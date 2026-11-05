import itertools as it
import json, os, time
import sympy as sp, numpy as np
import flint
from . import singular as sing
from .families import family
from .equations import EquationSystem
# from algebra import SparAlg

class Solver:
    def __init__(self, ansatz, htilde=None, symmetries=(), extra_equations=(), verbose=True):
        self.ans = ansatz
        self.d = ansatz.d
        self.es = EquationSystem(ansatz, htilde=htilde, verbose=verbose)
        self.symmetries = list(symmetries)
        self.results = {}
        self.verbose = verbose
        self.extra_sympy = list(extra_equations)
        self.nametable = self.es.nametable
        self.extra = []
        # self.log(f'symmetries of system {self.ans.name}: {self.symmetries}')

    def flint_to_sympy(self, eq, nametable):
        ctx = self.es.ctx
        eqdict = eq.to_dict()
        gens = ctx.gens()
        n = ctx.nvars()
        out = 0
        symbols = {sf: sp.Symbol(f'{nametable[sf]}', real=True) for sf in nametable.keys()}
        for k, c in eqdict.items():
            term = c
            for i in range(n):
                term *= symbols[str(gens[i])]**k[i]
            out += term
        return out

    def sympy_to_flint(self, eq):
        ctx = self.es.ctx
        syms = self.ans.syms
        poly = sp.Poly(sp.expand(eq), *syms)
        n = ctx.nvars()
        dct = dict()
        for m, c in poly.terms():
            if not sp.nsimplify(c).is_Rational:
                raise ValueError(f'c={c} is not rational')
            dct[tuple(list(m) + [0]*(n - len(m)))] = flint.fmpq(int(c.p), int(c.q))
        return ctx.from_dict(dct)

    def log(self, a):
        if self.verbose:
            print(f'[NCq3Solver] {a}', flush=True)

    def independent_linear(self, F):
        if not F: return []
        n = self.es.ctx.nvars()
        M = flint.fmpq_mat(len(F), n+1)     #+1 for the constant term
        for i, p in enumerate(F):
            for m, c in zip(p.monoms(), p.coeffs()):
                j = m.index(1) if sum(m) == 1 else n
                M[i, j] = c
        R, r = M.rref()
        out = []
        for i in range(r):
            dct = dict()
            for j in range(n):
                if R[i, j] != 0:
                    e = [0]*n
                    e[j] = 1
                    dct[tuple(e)] = R[i, j]
            if R[i, n] != 0:
                dct[(0,)*n] = R[i, n]
            out.append(self.es.ctx.from_dict(dct))
        return out

    def forms(self):
        es, A, d = self.es, self.es.alg, self.d
        h = es.density()
        self.block_forms = {}
        for (a, b) in it.product(range(d), repeat=2):
            F = []
            for (rw, cl), (re, im) in h.items():
                if rw[0] == a and rw[2] == b and rw[1] < cl[1]:
                    F += [p for p in (re, im) if not p.is_zero()]
            self.block_forms[(a, b)] = self.independent_linear(F)
        comps = A.canonical(h, 3)
        self.two_local_forms = [p for (rw, cl), v in comps[3].items() for p in v if not p.is_zero()]
        q = A.comm(A.pad(h, 0, 1), A.pad(h, 1, 0))
        comps_q = A.canonical(q, 4)
        # comps = {l: {(rw, cl): (re, im)}}, c is itself a dict
        self.vacuous_forms = [p for c in comps_q.values() for v in c.values() for p in v if not p.is_zero()]

    def generate(self, linear_reduce=True):
        eqs = self.es.generate()
        if linear_reduce:
            self.es.linear_reduce()
        self.extra = [self.sympy_to_flint(s) for s in self.extra_sympy]
        self.forms()
        return self.es.equations

    def active_candidates(self):
        return [k for k in it.product(range(self.d), repeat=2) if self.block_forms[k]]

    def consistent(self, S):
        def key(k):
            return tuple(sorted(map(str, self.block_forms[k])))
        act = {key(k) for k in S}
        return all(key(k) not in act for k in self.active_candidates() if k not in S)

    def orbit(self, S, do_log=False):
        orbit, frontier = {S}, [S]
        while frontier:
            X = frontier.pop()
            for g in self.symmetries:
                Y = frozenset(g['block'](*k) for k in X)
                if do_log: self.log(f'{set(X)} --> {set(Y)} by symmetry {g['name']}')
                if Y not in orbit:
                    if do_log: self.log(f'{set(Y)} was not in orbit')
                    orbit.add(Y)
                    frontier.append(Y)
                else:
                    if do_log: self.log(f'{set(Y)} already in orbit')
        return orbit

    def patterns(self, min_active=1, max_active=None):
        cands = self.active_candidates()
        max_active = max_active or len(cands)
        seen, reps = set(), []
        for r in range(min_active, max_active+1):
            for S in it.combinations(cands, r):
                S = frozenset(S)
                if not self.consistent(S): continue
                if S in seen: continue
                orbit = self.orbit(S)
                seen |= orbit
                reps.append((tuple(sorted(S)), len(orbit)))
        self.log(f'{len(reps)} unique activity patterns from {len(cands)} activatable blocks')
        return reps

    def stratum(self, S):
        S = sorted(S)
        inactive = [k for k in self.active_candidates() if k not in S]
        lin = [p for k in inactive for p in self.block_forms[k]]
        first = self.block_forms[S[0]]
        # we don't use k == S[0] so that the code instantly filters for redundancies in the active blocks (if block == first, then we can't use it for saturation later)
        others = [self.block_forms[k] for k in S if not self.block_forms[k] == first]
        charts = []
        for c in range(len(first)):
            eqs = lin + first[:c] + [first[c]-1]
            charts.append(dict(label=f'S={S} chart={c}', extra=eqs, saturate=others))
        return charts

    def system(self):
        return self.es.equations + self.extra

    # ==========================================================
    # everything below is just the wiring of the singular module
    # ==========================================================

    # ---- CHANGED (Claude): script builders split from execution, so charts can run in parallel ----
    def probe_script(self, chart, p, sat=False):
        return sing.build_script(self.es.var_names, self.system(), extra=chart['extra'],
                                 saturate=chart['saturate'] if sat else [], gb='slimgb',
                                 decompose='none', char=p)

    def solve_script(self, chart, gb='modStd0'):
        # facstd: cheap split into pieces that need not be prime; genuine pieces are refined afterwards
        return sing.build_script(self.es.var_names, self.system(), extra=chart['extra'], saturate=[],
                                 gb=gb, decompose='facstd', tests=self.tests_for(chart))

    def refine_script(self, comp, chart):
        # generator strings go back to Singular verbatim, no sympy round trip
        return sing.build_script(self.es.var_names, comp['generators'], gb='slimgb', decompose='minAssGTZ', tests=self.tests_for(chart))
    # ---- END CHANGED ----

    def probe_chart(self, chart, primes=(32003, 32009), timeout=120, sat_timeout=30, mem_gb=None):
        """single chart, sequential (kept for interactive use)"""
        def dims(sat, tmo):
            out = []
            for p in primes:
                r = sing.run(self.probe_script(chart, p, sat), timeout=tmo, mem_gb=mem_gb)   # CHANGED (Claude)
                out.append(r.get('dim') if r['status'] == 'ok' else r['status'])
            return out
        res = dict(label=chart['label'], dims_modp=dims(False, timeout))
        # 1 in I implies 1 in I:J^oo, so saturation is only tried if needed
        if all(x == -1 for x in res['dims_modp']):
            res['empty'] = True
            return res
        if chart['saturate']:
            res['dims_modp_sat'] = dims(True, sat_timeout)
            if all(x == -1 for x in res['dims_modp_sat']):
                res['empty'] = True
                return res
        res['empty'] = False
        return res

    def tests_for(self, chart):
        zero = [self.es.alg.zero]
        tests = {'twolocal': self.two_local_forms or zero, 'vacuous': self.vacuous_forms or zero}
        uniq = {tuple(sorted(map(str, J))): J for J in chart['saturate']}
        for k, J in enumerate(uniq.values()):
            tests[f'diag{k}'] = J
        return tests

    @staticmethod
    def classify(comp):
        t = comp.get('tests', {})
        if any(v for k, v in t.items() if k.startswith('diag')):
            return 'excluded'
        if t.get('vacuous'):
            return 'vacuous'
        if t.get('twolocal'):
            return 'two-local'
        return 'genuine'

    # ---- CHANGED (Claude): post-processing shared by solve_chart and the parallel solve_all ----
    def finish_solve(self, chart, res, gb):
        res['label'] = chart['label']
        res['gb_used'] = gb
        for comp in res['components']:
            comp['class'] = self.classify(comp)
        if res['status'] == 'ok':
            res.pop('raw_tail', None)
        return res

    def solve_chart(self, chart, gb='modStd0', timeout=600, mem_gb=None, refine=True):
        """facstd split, then minAssGTZ on every genuine piece, all within one time budget"""
        t0 = time.time()
        res = self.finish_solve(chart, sing.run(self.solve_script(chart, gb), timeout=timeout, mem_gb=mem_gb), gb)
        if res['status'] == 'timeout' and 'gb_generators' in res:
            # decomposition unfinished: keep the whole chart ideal as one unrefined piece
            piece = dict(index=0, dim=res['dim'], generators=res['gb_generators'], tests=res.get('gb_tests', {}),
                         refined=False, note='decomposition timed out; whole chart ideal, may be reducible')
            piece['class'] = self.classify(piece)
            res.update(status='partial', components=[piece])
            res.pop('raw_tail', None)
            return res
        if not refine or res['status'] != 'ok':
            return res
        comps, seen = [], set()
        for comp in res['components']:
            # trivial classes are already exact on a non-prime piece: T in the piece ideal means T vanishes on it
            if comp['class'] != 'genuine':
                comp['refined'] = False
                comps.append(comp)
                continue
            script = self.refine_script(comp, chart)
            left = timeout - (time.time() - t0)
            sub = sing.run(script, timeout=left, mem_gb=mem_gb) if left > 1 else None
            if sub is None or sub['status'] != 'ok':
                comp['refined'] = False
                comp['note'] = 'refinement ' + (sub['status'] if sub else 'skipped: no time left') + ', piece may not be prime'
                res['status'] = 'partial'
                comps.append(comp)
                continue
            for c in sub['components']:
                key = tuple(sorted(c.get('generators', [])))
                if key in seen:
                    continue
                seen.add(key)
                c['class'] = self.classify(c)
                c['refined'] = True
                comps.append(c)
        res['components'] = comps
        res['wall'] = time.time() - t0
        return res
    # ---- END CHANGED ----

    # ---- CHANGED (Claude): parallel solve_all, your out_json resumability kept ----
    def solve_all(self, patterns=None, probe=True, out_json=None, workers=None, primes=(32003, 32009),
                  probe_timeout=120, sat_timeout=30, gb='modStd0', timeout=600, mem_gb=None):
        """three parallel stages: modular probe, saturated modular probe (only where needed), exact solve"""
        patterns = patterns if patterns is not None else self.patterns()
        results = json.load(open(out_json)) if out_json and os.path.exists(out_json) else {}
        t_all = time.time()

        def save():
            if out_json:
                json.dump(results, open(out_json, 'w'), indent=1)

        charts = [(chart, orb) for S, orb in patterns for chart in self.stratum(S)]
        # resuming: finished charts are skipped, timeouts and errors are tried again
        todo = [k for k in range(len(charts))
                if results.get(charts[k][0]['label'], {}).get('status') not in ('ok', 'empty_modp')]
        if len(todo) < len(charts):
            self.log(f'resuming: {len(charts) - len(todo)} of {len(charts)} charts already in {out_json}')

        def probe_stage(idx, sat, tmo, key):
            jobs = [(k, p) for k in idx for p in primes]
            out = sing.run_many([self.probe_script(charts[k][0], p, sat) for k, p in jobs],
                                workers=workers, timeout=tmo, mem_gb=mem_gb)
            dims = {k: [] for k in idx}
            for (k, p), r in zip(jobs, out):
                dims[k].append(r.get('dim') if r['status'] == 'ok' else r['status'])
            remaining = []
            for k in idx:
                chart, orb = charts[k]
                entry = results.setdefault(chart['label'], {})
                entry.setdefault('probe', {})[key] = dims[k]
                if all(x == -1 for x in dims[k]):
                    entry.update(status='empty_modp', components=[], orbit_size=orb)
                    self.log(f"{chart['label']}: empty mod p")
                else:
                    remaining.append(k)
            save()
            return remaining

        if probe and todo:
            todo = probe_stage(todo, False, probe_timeout, 'dims_modp')
            sat_idx = [k for k in todo if charts[k][0]['saturate']]
            if sat_idx:
                emptied = set(sat_idx) - set(probe_stage(sat_idx, True, sat_timeout, 'dims_modp_sat'))
                todo = [k for k in todo if k not in emptied]
        self.log(f'{len(todo)} of {len(charts)} charts go to the exact solve '
                 f'({workers or "all"} workers; {time.time() - t_all:.1f}s so far)')

        # on_done runs in the main thread, so writing results and the JSON file here is safe
        def done(i, r):
            chart, orb = charts[todo[i]]
            r['orbit_size'] = orb
            r['probe'] = results.get(chart['label'], {}).get('probe')
            results[chart['label']] = r
            ng = sum(c['class'] == 'genuine' for c in r['components'])
            self.log(f"{chart['label']}: {r['status']} dim={r.get('dim')} comps={len(r['components'])} "
                     f"genuine={ng} ({r['wall']:.1f}s)")
            save()

        jobs = [lambda k=k: self.solve_chart(charts[k][0], gb, timeout, mem_gb) for k in todo]
        sing.run_jobs(jobs, workers=workers, on_done=done)
        self.log(f'solve_all finished in {time.time() - t_all:.1f}s')
        self.results = results
        return results
    # ---- END CHANGED ----

    def genuine_components(self, results=None):
        results = results if results is not None else self.results
        return [(label, c) for label, r in results.items() for c in r['components'] if c['class'] == 'genuine']

    # ---------- readable output ----------
    def pretty(self, generator_strings):
        loc = {v: sp.Symbol(v) for v in self.es.var_names}
        table = {sp.Symbol(v): sp.Symbol(nm) for v, nm in self.es.nametable.items()}
        return [sp.sympify(g.replace('^', '**'), locals=loc).xreplace(table) for g in generator_strings]

    def parse_generators(self, gens):
        loc = {v: sp.Symbol(v) for v in self.es.var_names}
        syms = [loc[v] for v in self.es.var_names]
        out = []
        for g in gens:
            poly = sp.Poly(sp.expand(sp.sympify(g.replace('^', '**'), locals=loc)), *syms)
            dct = {}
            for mon, c in poly.terms():
                c = sp.nsimplify(c)
                dct[tuple(mon)] = flint.fmpq(int(c.p), int(c.q))
            out.append(self.es.ctx.from_dict(dct))
        return out