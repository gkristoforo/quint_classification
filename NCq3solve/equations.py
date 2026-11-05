import itertools as it
import time
import flint
from sympy import re as sre, im as sim
from .algebra import SparAlg

class EquationSystem:
    def __init__(self, ansatz, htilde=None, q5_form = 'irf', order = 'degrevlex', verbose=True):
        self.ans = ansatz
        self.d = ansatz.d
        d = self.d
        names = [s.name for s in ansatz.syms]
        self.aux_names = []
        if htilde == 'diag-gauge':
            self.aux_names = [f'g{a}{c}' for a in range(d) for c in range(d)]
        elif htilde == 'irf':
            for a, b in it.product(range(d), repeat=2):
                for i in range(d):
                    self.aux_names.append(f't{a}{b}_u{i}')
                    for j in range(i+1, d):
                        self.aux_names.append(f't{a}{b}_u{i}{j}')
                        self.aux_names.append(f't{a}{b}_v{i}{j}')
        elif htilde != None:
            raise ValueError(f'not valid htilde format, htilde={htilde}')
        self.names = names + self.aux_names
        self.var_names = [f'x{k}' for k in range(len(self.names))]
        self.ctx = flint.fmpq_mpoly_ctx.get(tuple(self.var_names), order)
        self.alg = SparAlg(self.ctx, d)
        self.htilde = htilde
        self.q5_form = q5_form
        self.aux_dict = self.aux_density()[0]
        self.equations = None      
        self.verbose = verbose

    @property
    def nvars(self):
        return len(self.names)

    def const(self, q):
        return self.ctx.from_dict({(0,)*self.ctx.nvars(): q}) if q != 0 else self.alg.zero

    def var(self, k, q=1):
        a = [0 for i in range(self.ctx.nvars())]
        a[k] = 1
        return self.ctx.from_dict({tuple(a): q})

    def entry(self, c, lin):
        R = self.const(flint.fmpq(int(sre(c).p), int(sre(c).q)))
        I = self.const(flint.fmpq(int(sim(c).p), int(sim(c).q)))
        for k, v in lin.items():
            vre = sre(v)
            vim = sim(v)
            if vre != 0:
                rk = self.var(k, flint.fmpq(int(vre.p), int(vre.q)))
                R = R + rk
            if vim != 0:
                ik = self.var(k, flint.fmpq(int(vim.p), int(vim.q)))
                I = I + ik
        return (R, I)

    def density(self):
        D = dict()
        for (a, b, i, j), (const, lin) in self.ans.linear.items():
            D[((a, i, b), (a, j, b))] = self.entry(const, lin)
        return self.alg.clean(D)

    def aux_density(self):
        d = self.d
        n0 = len(self.ans.syms)
        D_aux = dict()
        if self.htilde == 'diag-gauge':
            for i, (a, c) in enumerate(it.product(range(d), repeat=2)):
                D_aux[(a, c), (a, c)] = (self.var(n0 + i), self.alg.zero)
            return D_aux, 2
        if self.htilde == 'irf':
            n = n0
            for a, b in it.product(range(d), repeat=2):
                for i in range(d):
                    D_aux[((a, i, b), (a, i, b))] = (self.alg.zero, self.var(n))
                    n += 1
                    for j in range(i+1, d):
                        u = self.var(n)
                        v = self.var(n + 1)
                        n += 2
                        D_aux[((a, i, b), (a, j, b))] = (-v, u)
                        D_aux[((a, j, b), (a, i, b))] = (v, u)
            return D_aux, 3
        return None, 0

    def generate(self):
        t0 = time.time()
        h = self.density()
        A = self.alg
        if self.q5_form == 'irf':
            q = A.comm(A.pad(h, 0, 1), A.pad(h, 1, 0))
            qn = 4
        else:
            q = A.comm(A.pad(h, 0, 2), A.add(A.pad(h, 1, 1), A.pad(h, 2, 0)))
            qn = 5
        # range(h) + range(q) - smallest overlap(1)
        N = 3 + qn - 1
        R = A.ti_comm(h, q, N)
        aux, an = self.aux_density()
        if aux is not None:
            if self.htilde == 'diag-gauge':
                HG = A.ti_comm(h, aux, 4)
                R = A.add(R, A.pad(A.ti_comm(h, HG, 6), 0, N - 6))
            else:
                R = A.add(R, A.pad(A.ti_comm(h, aux, 5), 0, N - 5))
        R = A.clean(R)
        t1 = time.time()
        Rcan = A.canonical(R, N)
        t2 = time.time()
        eqs = []
        seen = set()
        for l, c in Rcan.items():
            for (rw, cl), (re, im) in c.items():
                if (cl, rw) in c and (cl, rw) < (rw, cl):
                    continue
                for p in (re, im):
                    if p.is_zero():
                        continue
                    p = p / p.leading_coefficient()
                    key = str(p)
                    if key not in seen:
                        eqs.append(p)
                        seen.add(key)
        self.equations = eqs
        self.stats = dict(nvars=self.ctx.nvars(), n_eqs_raw=len(eqs), t_comms=t1-t0, 
                          t_canonical=t2-t1, n_terms_res=len(R), 
                          max_degree=max([p.total_degree() for p in eqs], default=0), n_monoms=len({m for p in eqs for m in p.monoms()}))
        if self.verbose:
            print(f'[NCq3 EqSystem] eq generation: {self.stats}')
        return eqs

    def linear_reduce(self):
        eqs = self.equations
        monos = sorted({m for p in eqs for m in p.monoms()}, reverse=True)
        col = {m: k for k, m in enumerate(monos)}
        M = flint.fmpq_mat(len(eqs), len(monos))
        for i, p in enumerate(eqs):
            for m, c in p.terms():
                M[i, col[m]] = c
        t0 = time.time()
        Rm, rank = M.rref()
        t1 = time.time()
        out = []
        for i in range(rank):
            out.append(self.ctx.from_dict({monos[j]: Rm[i, j] for j in range(len(monos)) if Rm[i, j] != 0}))
        self.equations = out
        if self.verbose:
            stats = dict(shape=(M.nrows(), M.ncols()), before=len(eqs),
                         after=rank, t_rref=t1-t0)              
            print(f'[NCq3 EqSystem] linear reduce: {stats}')
        return out

    @property
    def nametable(self):
        return dict(zip(self.var_names, self.names))