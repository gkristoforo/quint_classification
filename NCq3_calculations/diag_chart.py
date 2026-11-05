"""Single-chart diagnosis: where does the time go, over Q and modulo p, with no parallel contention.

    python diag_chart.py --pattern 00,01,10,11,12 --chart 0 --exp p_full q_slimgb q_full --timeout 14400
"""
import argparse, json, os, re, shutil, signal, subprocess, sys, tempfile, time, resource
import re, os, time, json, sys; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

SINGULAR = shutil.which('Singular')
PRIME = 32003


def build_solver():
    # ADAPT: build exactly as in the KP run (same ansatz, correction, symmetries, extra equations),
    # otherwise --chart refers to a different chart than the one that timed out.
    from NCq3solve.families import family, symmetries_for
    from NCq3solve.pipeline import Solver
    ans, slice_eqs = family('KP_real')
    s = Solver(ans, symmetries=symmetries_for('KP_real'), extra_equations=slice_eqs, verbose=False)
    s.generate()
    return s


def poly_str(p):
    s = str(p)
    return s if s else '0'


def ideal_decl(name, polys):
    body = ',\n  '.join(poly_str(p) for p in polys) if polys else '0'
    return f'ideal {name} =\n  {body};\n'


def chart_tests(solver, chart):
    tests = {'twolocal': solver.two_local_forms or [solver.es.alg.zero],
             'vacuous': solver.vacuous_forms or [solver.es.alg.zero]}
    uniq = {tuple(sorted(map(str, J))): J for J in chart['saturate']}
    for k, J in enumerate(uniq.values()):
        tests[f'diag{k}'] = J
    return tests, list(uniq.values())


def build_script(var_names, gens, tests, sat_ideals, logpath, char, gb, decompose, refine, cores):
    s = ['LIB "primdec.lib";\nLIB "modstd.lib";\nLIB "resources.lib";\n',
         f'setcores({int(cores)});\nsystem("--ticks-per-sec", 1000);\n',
         f'string LOGF = ":a {logpath}";\n',
         'proc mark(string m) { write(LOGF, m + " cpu=" + string(timer) + " wall=" + string(rtimer)); }\n',
         f"ring R = {char},({','.join(var_names)}),dp;\noption(redSB);\n",
         ideal_decl('I', gens)]
    for name, polys in tests.items():
        s.append(ideal_decl(f'T_{name}', polys))
    s.append(f'mark("@@START ngens=" + string(size(I)));\n')
    # markers go to a file opened and closed per write, so they survive a SIGKILL
    for k, J in enumerate(sat_ideals):
        s.append(ideal_decl(f'J{k}', J) + f'I = sat(I, J{k});\nmark("@@SAT k={k} ngens=" + string(size(I)));\n')
    gbcall = {'slimgb': 'slimgb(I)', 'std': 'std(I)', 'modStd0': 'modStd(I, 0)', 'modStd': 'modStd(I)'}[gb]
    s.append(f'ideal G = {gbcall};\nmark("@@GB size=" + string(size(G)) + " dim=" + string(dim(G)));\n')
    if decompose != 'none':
        names = list(tests)
        flags = ''.join(f'{n}=" + string(t_{n}) + " ' for n in names)
        diag_any = ' || '.join(f't_{n}' for n in names if n.startswith('diag')) or '0'
        s.append(''.join(f'int t_{n};\n' for n in names))
        s.append('int genuine; int k; int m;\n')
        s.append('if (dim(G) == -1) { mark("@@EMPTY"); }\nelse {\n')
        s.append(f'list P = {decompose}(G);\nmark("@@DECOMP method={decompose} n=" + string(size(P)));\n')
        s.append('for (k = 1; k <= size(P); k++) {\n  ideal Pk = std(P[k]);\n'
                 '  if (Pk[1] == 1) { mark("@@PIECE k=" + string(k) + " empty=1"); continue; }\n')
        s.append(''.join(f'  t_{n} = (size(reduce(T_{n}, Pk)) == 0);\n' for n in names))
        s.append(f'  genuine = !({diag_any}) && !t_vacuous && !t_twolocal;\n')
        s.append(f'  mark("@@PIECE k=" + string(k) + " dim=" + string(dim(Pk)) + " {flags}genuine=" + string(genuine));\n')
        if refine and decompose == 'facstd':
            s.append('  if (genuine) {\n    list Q = minAssGTZ(Pk);\n'
                     '    mark("@@REFINE k=" + string(k) + " nprimes=" + string(size(Q)));\n'
                     '    for (m = 1; m <= size(Q); m++) {\n      ideal Qm = std(Q[m]);\n')
            s.append(''.join(f'      t_{n} = (size(reduce(T_{n}, Qm)) == 0);\n' for n in names))
            s.append(f'      genuine = !({diag_any}) && !t_vacuous && !t_twolocal;\n')
            s.append(f'      mark("@@PRIME k=" + string(k) + " m=" + string(m) + " dim=" + string(dim(Qm)) + " {flags}genuine=" + string(genuine));\n'
                     '      kill Qm;\n    }\n    kill Q;\n  }\n')
        s.append('  kill Pk;\n}\n}\n')
    s.append('mark("@@DONE");\nquit;\n')
    return ''.join(s)


EXPERIMENTS = {   # name: (characteristic, gb, decompose, refine)
    'p_gb':     (PRIME, 'slimgb', 'none', False),
    'p_full':   (PRIME, 'slimgb', 'minAssGTZ', False),
    'q_slimgb': (0, 'slimgb', 'none', False),
    'q_modstd': (0, 'modStd0', 'none', False),
    'q_full':   (0, 'modStd0', 'facstd', True),
}


def peak_child_gib():
    r = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    return r / 2**30 if sys.platform == 'darwin' else r / 2**20     # bytes on macOS, KiB on Linux


def run_singular(script, logpath, timeout):
    with tempfile.NamedTemporaryFile('w', suffix='.sing', delete=False) as f:
        f.write(script); path = f.name
    t0 = time.time()
    proc = subprocess.Popen([SINGULAR, '-q', '--no-warn', path], stdin=subprocess.DEVNULL,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, start_new_session=True)
    try:
        out, _ = proc.communicate(timeout=timeout); status = 'ok'
    except subprocess.TimeoutExpired:
        os.killpg(proc.pid, signal.SIGKILL)       # also kills modStd's forked workers
        out, _ = proc.communicate(); status = 'timeout'
    wall = time.time() - t0
    os.unlink(path)
    if status == 'ok' and ('? ' in out or proc.returncode != 0):
        status = 'error'
    marks = []
    if os.path.exists(logpath):
        for line in open(logpath):
            parts = line.split()
            kv = dict(p.split('=', 1) for p in parts[1:] if '=' in p)
            marks.append(dict(tag=parts[0], cpu_s=int(kv.pop('cpu')) / 1000, wall_s=int(kv.pop('wall')) / 1000, **kv))
    return dict(status=status, wall_s=round(wall, 1), peak_gib_so_far=round(peak_child_gib(), 2),
                last_stage=marks[-1]['tag'] if marks else None, marks=marks, tail=out[-1500:])


def summarise(name, r):
    print(f"\n== {name}: {r['status']}  wall={r['wall_s']}s  last stage={r['last_stage']}")
    prev = 0.0
    for m in r['marks']:
        extra = ' '.join(f'{k}={v}' for k, v in m.items() if k not in ('tag', 'cpu_s', 'wall_s'))
        print(f"   {m['tag']:<9} cpu={m['cpu_s']:9.1f}s (+{m['cpu_s'] - prev:8.1f})  wall={m['wall_s']:9.1f}s  {extra}")
        prev = m['cpu_s']
    if r['status'] == 'error':
        print('   Singular output tail:\n' + r['tail'])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--pattern', required=True, help='active blocks, e.g. 00,01,10,11,12')
    ap.add_argument('--chart', type=int, default=0)
    ap.add_argument('--exp', nargs='+', default=['p_full', 'q_slimgb', 'q_full'], choices=list(EXPERIMENTS))
    ap.add_argument('--timeout', type=float, default=4 * 3600, help='wall seconds per experiment')
    ap.add_argument('--saturate', action='store_true', help="saturate before the GB instead of 'test' mode")
    ap.add_argument('--cores', type=int, default=1)
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    if SINGULAR is None:
        sys.exit('Singular not found in PATH')

    S = [(int(b[0]), int(b[1])) for b in a.pattern.split(',')]
    solver = build_solver()
    charts = solver.stratum(S)
    chart = charts[a.chart]
    gens = list(solver.es.equations) + list(solver.extra) + list(chart['extra'])
    tests, sat_ideals = chart_tests(solver, chart)
    if a.saturate:
        tests = {k: v for k, v in tests.items() if not k.startswith('diag')}
    else:
        sat_ideals = []
    print(f"{chart['label']}  ({len(charts)} charts in this stratum)")
    print(f"  {len(solver.es.equations)} equations + {len(solver.extra)} family + {len(chart['extra'])} chart "
          f"= {len(gens)} generators in {len(solver.es.var_names)} variables; "
          f"{len(chart_tests(solver, chart)[1])} distinct diag ideals ({'saturated' if a.saturate else 'tests'}); "
          f"max degree {max(p.total_degree() for p in gens)}; {sum(len(p) for p in gens)} terms")

    report = dict(label=chart['label'], args=vars(a), experiments={})
    out = a.out or 'diagnostics/diag_' + re.sub(r'[^0-9a-z]+', '_', chart['label'].lower()).strip('_') + '.json'
    for name in a.exp:
        char, gb, dec, ref = EXPERIMENTS[name]
        logpath = os.path.abspath(f'diagnostics/.diag_{os.getpid()}_{name}.log')
        if os.path.exists(logpath): os.unlink(logpath)
        script = build_script(solver.es.var_names, gens, tests, sat_ideals, logpath, char, gb, dec, ref, a.cores)
        r = run_singular(script, logpath, a.timeout)
        os.unlink(logpath) if os.path.exists(logpath) else None
        summarise(name, r)
        report['experiments'][name] = r
        json.dump(report, open(out, 'w'), indent=1)       # rewritten after each experiment
    print(f'\nreport: {out}')


if __name__ == '__main__':
    main()
