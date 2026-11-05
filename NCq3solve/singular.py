"""
=====================================================================
  ENTIRE FILE ADDED BY CLAUDE (Singular back end)
  Script generation, safe execution (single or parallel), parsing of the
  @@ output protocol.
  Knows nothing about physics: equations, saturation and test ideals
  are supplied by pipeline.py.

  CHANGED (Claude, memory watchdog): rlimit/ulimit replaced by a polling
  watchdog (macOS does not enforce RLIMIT_AS). Three guards, all of which
  kill the whole Singular process tree:
    - per-process limit mem_gb,
    - shared budget for all Singular processes of this Python process,
    - reserve of free system memory.
  Output goes to a file (no pipe deadlock, partial output survives a kill).
=====================================================================
"""
import subprocess
import tempfile
import threading
import time
import os
import shutil
from concurrent.futures import ThreadPoolExecutor, as_completed

import psutil

SINGULAR = shutil.which('Singular')
GiB = 2**30
MIN_VICTIM = 0.25 * GiB      # shared guards never kill a process smaller than this


def poly_str(p):
    s = str(p)
    return s if s else '0'


def ideal_decl(name, polys):
    body = ',\n  '.join(poly_str(p) for p in polys) if polys else '0'
    return f'ideal {name} =\n  {body};\n'


def build_script(var_names, equations, extra=(), saturate=(), gb='modStd0', decompose='minAssGTZ',
                 tests=None, char=0, order='dp', cores=1):
    # modular algorithms only make sense over Q
    if char != 0 and gb in ('modStd', 'modStd0'):
        gb = 'slimgb'
    s = ['LIB "primdec.lib";\nLIB "modstd.lib";\nLIB "resources.lib";\n',
         f'setcores({int(cores)});\noption(redSB);\n',
         f'ring R = {char},({",".join(var_names)}),{order};\n']
    uniq, seen = [], set()
    for J in saturate:
        key = tuple(sorted(poly_str(p) for p in J))
        if key not in seen:
            seen.add(key)
            uniq.append(J)
    s.append(ideal_decl('I', list(equations) + list(extra)))
    s.append('int t0 = timer;\n')
    for k, J in enumerate(uniq):
        s.append(ideal_decl(f'J{k}', J))
        s.append(f'I = sat(I, J{k});\n')
    gbcall = {'std': 'std(I)', 'slimgb': 'slimgb(I)', 'modStd': 'modStd(I)', 'modStd0': 'modStd(I, 0)'}[gb]
    s.append(f'ideal G = {gbcall};\n')
    s.append('print("@@TGB " + string(timer - t0));\n')
    s.append('print("@@GBSIZE " + string(size(G)));\nprint("@@DIM " + string(dim(G)));\n')
    s.append('if (dim(G) == 0) { print("@@VDIM " + string(vdim(G))); }\n')
    if tests:
        for name, polys in tests.items():
            s.append(ideal_decl(f'T_{name}', polys))
        # whole chart ideal and its tests, printed before decomposing: the fallback on timeout
        s.append('if (dim(G) >= 0) {\n  print("@@GBGENS " + string(G));\n')
        for name in tests:
            s.append(f'  print("@@GBTEST {name} " + string(size(reduce(T_{name}, G)) == 0));\n')
        s.append('}\n')
    if decompose and decompose != 'none':
        # decomposing the unit ideal is a Singular error, hence the guard
        s.append('if (dim(G) >= 0) {\n')
        s.append(f'list P = {decompose}(G);\n')
        s.append('print("@@NCOMP " + string(size(P)));\n')
        s.append('int k; for (k = 1; k <= size(P); k++) {\n'
                 '  ideal Pk = std(P[k]);\n'
                 '  if (Pk[1] == 1) { print("@@EMPTY " + string(k)); continue; }\n'
                 '  print("@@COMP " + string(k) + " DIM " + string(dim(Pk)));\n'
                 '  print("@@GENS " + string(P[k]));\n')
        if tests:
            for name in tests:
                s.append(f'  print("@@TEST {name} " + string(size(reduce(T_{name}, Pk)) == 0));\n')
        s.append('}\n}\nelse { print("@@NCOMP 0"); }\n')
    s.append('print("@@TIME " + string(timer - t0));\nquit;\n')
    return ''.join(s)


# ---------------------------------------------------------------------
#  memory watchdog
# ---------------------------------------------------------------------
class MemoryGuard:
    """shared state of all running Singular processes; the largest one is killed when a shared limit is hit"""

    def __init__(self):
        self.lock = threading.Lock()
        self.live = {}            # pid -> current rss in bytes
        self.doomed = {}          # pid -> reason
        self.total_bytes = None
        self.reserve_bytes = None
        self.last_kill = 0.0
        self.configure()

    def configure(self, total_gb='auto', reserve_gb='auto'):
        vm = psutil.virtual_memory()
        if total_gb == 'auto':
            total_gb = 0.6 * vm.total / GiB
        if reserve_gb == 'auto':
            # macOS reports little 'available' memory even when idle (cache, compression):
            # the reserve is relative to what is available now, so it cannot fire at baseline
            reserve_gb = min(2.0, 0.5 * vm.available / GiB)
        with self.lock:
            self.total_bytes = None if total_gb is None else int(total_gb * GiB)
            self.reserve_bytes = None if reserve_gb is None else int(reserve_gb * GiB)

    def update(self, pid, rss):
        """records rss; returns the kill reason if this process has been chosen, else None"""
        with self.lock:
            self.live[pid] = rss
            if pid in self.doomed:
                return self.doomed.pop(pid)
            # cool-down: the memory of the last victim needs a moment to be returned to the system
            if time.time() - self.last_kill < 2.0:
                return None
            reason = None
            if self.total_bytes is not None and sum(self.live.values()) > self.total_bytes:
                reason = 'memory-budget'
            elif self.reserve_bytes is not None and psutil.virtual_memory().available < self.reserve_bytes:
                reason = 'memory-reserve'
            if reason is None:
                return None
            victim = max(self.live, key=self.live.get)
            if self.live[victim] < MIN_VICTIM:
                return None       # the pressure comes from elsewhere; killing a small run would not help
            self.live[victim] = 0
            self.last_kill = time.time()
            if victim == pid:
                return reason
            self.doomed[victim] = reason
            return None

    def release(self, pid):
        with self.lock:
            self.live.pop(pid, None)
            self.doomed.pop(pid, None)

    def status(self):
        with self.lock:
            return dict(n_running=len(self.live), rss_gb=sum(self.live.values()) / GiB,
                        total_gb=self.total_bytes and self.total_bytes / GiB,
                        reserve_gb=self.reserve_bytes and self.reserve_bytes / GiB,
                        available_gb=psutil.virtual_memory().available / GiB)


GUARD = MemoryGuard()


def configure_memory(total_gb='auto', reserve_gb='auto'):
    """total_gb: budget for all Singular processes together ('auto' = 60% of RAM, None = off);
    reserve_gb: kill the largest process when available system memory drops below this
    ('auto' = min(2, half of what is available now), None = off)"""
    GUARD.configure(total_gb, reserve_gb)


def _tree(proc):
    try:
        return [proc] + proc.children(recursive=True)
    except psutil.NoSuchProcess:
        return []


def _tree_rss(proc):
    total = 0
    for q in _tree(proc):
        try:
            total += q.memory_info().rss
        except psutil.NoSuchProcess:
            pass
    return total


def _kill_tree(proc):
    procs = _tree(proc)
    for q in procs:
        try:
            q.kill()
        except psutil.NoSuchProcess:
            pass
    psutil.wait_procs(procs, timeout=5)


def run(script, timeout=600, mem_gb=None, workdir=None, poll=0.25):
    if SINGULAR is None:
        raise RuntimeError('Singular not found in PATH')
    with tempfile.NamedTemporaryFile('w', suffix='.sing', delete=False, dir=workdir) as f:
        f.write(script)
        path = f.name
    t0, peak, killed = time.time(), 0, None
    limit = None if mem_gb is None else int(mem_gb * GiB)
    try:
        with tempfile.TemporaryFile('w+', dir=workdir) as out:
            # stdin closed: on a script error Singular would otherwise wait for terminal input forever
            p = subprocess.Popen([SINGULAR, '-q', '--no-warn', path], stdin=subprocess.DEVNULL,
                                 stdout=out, stderr=subprocess.STDOUT, start_new_session=True)
            proc = psutil.Process(p.pid)
            try:
                while True:
                    try:
                        p.wait(timeout=poll)
                        break
                    except subprocess.TimeoutExpired:
                        pass
                    rss = _tree_rss(proc)
                    peak = max(peak, rss)
                    if limit is not None and rss > limit:
                        killed = 'memory'
                    elif time.time() - t0 > timeout:
                        killed = 'timeout'
                    else:
                        killed = GUARD.update(p.pid, rss)
                    if killed:
                        _kill_tree(proc)
                        p.wait()
                        break
            finally:
                GUARD.release(p.pid)
                if p.poll() is None:          # interrupted from Python (e.g. Ctrl-C): no orphans
                    _kill_tree(proc)
            out.seek(0)
            text = out.read()
    finally:
        os.unlink(path)
    if killed and killed != 'timeout':
        text += f'\n[watchdog] killed: {killed}, peak {peak / GiB:.2f} GiB'
    status = 'timeout' if killed == 'timeout' else ('error' if killed else 'ok')
    res = parse(text, status, time.time() - t0)
    res['killed'] = killed
    res['peak_gb'] = round(peak / GiB, 3)
    return res


def run_jobs(jobs, workers=1, on_done=None):
    """runs zero-argument callables in a thread pool; results in input order, on_done(i, result) as each finishes"""
    workers = workers or max(1, (os.cpu_count() or 1) - 2)
    results = [None] * len(jobs)
    # threads suffice: each one only waits for its own Singular processes
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(job): i for i, job in enumerate(jobs)}
        for fut in as_completed(futures):
            i = futures[fut]
            try:
                results[i] = fut.result()
            except Exception as e:
                results[i] = dict(status='error', wall=0.0, components=[], raw_tail=repr(e),
                                  killed=None, peak_gb=0.0)
            if on_done is not None:
                on_done(i, results[i])
    return results


def run_many(scripts, workers=None, timeout=600, mem_gb=None, on_done=None):
    """runs independent scripts in parallel, one single-core Singular process per worker"""
    jobs = [lambda sc=sc: run(sc, timeout, mem_gb) for sc in scripts]
    return run_jobs(jobs, workers=workers, on_done=on_done)


def check_memory_guard(mem_gb=0.3, timeout=60):
    """run once on a new machine: a script that allocates without bound must come back killed='memory'.
    The shared guards are switched off during the test, so only the per-process limit is tested."""
    hog = 'list L; int i;\nfor (i = 1; i <= 10000000; i++) { L[i] = bigint(7)^200000; }\nquit;\n'
    vm = psutil.virtual_memory()
    print(f"system: total {vm.total / GiB:.1f} GiB, available {vm.available / GiB:.1f} GiB")
    saved = (GUARD.total_bytes, GUARD.reserve_bytes)
    GUARD.total_bytes = GUARD.reserve_bytes = None
    try:
        r = run(hog, timeout=timeout, mem_gb=mem_gb)
    finally:
        GUARD.total_bytes, GUARD.reserve_bytes = saved
    ok = r['killed'] == 'memory' and r['status'] == 'error' and r['peak_gb'] >= 0.9 * mem_gb
    print(f"memory guard {'OK' if ok else 'NOT WORKING'}: killed={r['killed']}, "
          f"peak={r['peak_gb']:.2f} GiB, wall={r['wall']:.1f}s")
    print(f"shared guards: {GUARD.status()}")
    return ok


def parse(text, status, wall):
    res = dict(status=status, wall=wall, components=[], raw_tail=text[-2000:])
    if status == 'ok' and ('? ' in text or '@@TIME' not in text):
        res['status'] = 'error'
    cur = None
    for line in text.splitlines():
        if line.startswith('@@GBSIZE'):
            res['gb_size'] = int(line.split()[1])
        elif line.startswith('@@DIM'):
            res['dim'] = int(line.split()[1])
        elif line.startswith('@@VDIM'):
            res['vdim'] = int(line.split()[1])
        elif line.startswith('@@NCOMP'):
            res['n_components_raw'] = int(line.split()[1])
        elif line.startswith('@@COMP'):
            parts = line.split()
            cur = dict(index=int(parts[1]), dim=int(parts[3]), tests={})
            res['components'].append(cur)
        elif line.startswith('@@GENS') and cur is not None:
            cur['generators'] = [g.strip() for g in line[len('@@GENS'):].split(',') if g.strip()]
        elif line.startswith('@@TEST') and cur is not None:
            _, name, val = line.split()
            cur['tests'][name] = (val == '1')
        elif line.startswith('@@GBGENS'):
            res['gb_generators'] = [g.strip() for g in line[len('@@GBGENS'):].split(',') if g.strip()]
        elif line.startswith('@@GBTEST'):
            _, name, val = line.split()
            res.setdefault('gb_tests', {})[name] = (val == '1')
        elif line.startswith('@@TGB'):
            res['cpu_s_gb'] = int(line.split()[1])
        elif line.startswith('@@TIME'):
            res['cpu_s_total'] = int(line.split()[1])
    return res
