import sympy as sp
import numpy as np
from .operators import SPIN1, PAULI, LABELS, LABELS_DIP_REAL, LABELS_FULL, LABELS_REAL, KP_DICT
import itertools as it

class Ansatz:
    # hdict = {(a, b): some matrix}
    def __init__(self, hdict, d=3, check_hermitian=True, name=None, verbose=True, rebuild=False):
        self.hdict = {(a, b): sp.Matrix(M) for (a, b), M in hdict.items()}
        self.d = d
        self.name = name
        self.verbose = verbose
        syms = set()
        hkeys = self.hdict.keys()
        for a, b in it.product(range(d), repeat=2):
            if (a, b) not in hkeys:
                self.hdict[(a, b)] = sp.zeros(d, d)
            else: 
                m_ab = sp.Matrix(self.hdict[(a, b)])
                if m_ab != 0:
                    for i, j in it.product(range(d), repeat=2):
                        syms_ij = m_ab[i, j].free_symbols
                        syms = syms.union(syms_ij)
        for M in hdict.values():
            syms = syms.union(*[ij.free_symbols for ij in M])
        self.syms = sorted(syms, key=lambda s: s.name)
        for s in self.syms:
            if not s.is_real:
                raise ValueError(f'the symbol {s} is not real')
        if rebuild:
            return
        if check_hermitian:
            # this checks whether all blocks specified in hdict are Hermitian (the ansatz is forced to ensure hermiticity by declaring all blocks in hdict Hermitian)
            for (a, b), M in self.hdict.items():
                if sp.simplify(M - M.adjoint()) != sp.zeros(d, d):
                    # print('Hehe, it is not Hermitian, whatever\n', sp.latex(sp.simplify(M - M.adjoint())))
                    raise ValueError(f'the block {(a, b)} is not Hermitian, the matrix is: {M}')
        indices = {s: n for n, s in enumerate(self.syms)}
        self.linear = {}
        for (a, b), M in self.hdict.items():
            for i, j in it.product(range(self.d), repeat=2):
                eq = sp.expand(M[i, j])
                if eq == 0:
                    continue
                poly = sp.Poly(eq, *self.syms) if len(self.syms) != 0 else None
                if poly is not None and poly.total_degree() > 1:
                    raise ValueError(f'slot {(i, j)} in matrix {M} has degree {poly.total_degree()}')
                const = eq.subs({s: 0 for s in eq.free_symbols})
                lin_coeffs = {indices[s]: sp.nsimplify(sp.diff(eq, s)) for s in eq.free_symbols if sp.diff(eq, s) != 0}
                self.linear[(a, b, i, j)] = (sp.nsimplify(const), lin_coeffs)
        self.check_rational()
        if self.verbose:
            stats = dict(name=self.name, nsyms=self.nvars, d=self.d,
                        lin_forms=len(self.linear))
            print(f'[NCq3Ansatz] {stats}')

    @property
    def nvars(self):
        return len(self.syms)

    def check_rational(self):
        for (a, b, i, j), (c, lin_dict) in self.linear.items():
            for x in [c, *lin_dict.values()]:
                re, im = sp.re(x), sp.im(x)
                if not (re.is_Rational and im.is_Rational):
                    raise ValueError(f'a coefficient of M_{i}{j} in block {(a, b)} is not rational')
        return None
    
    def numeric_block(self, sym_vals):
        num_hdict = self.hdict.copy()
        for (a, b), M in num_hdict.items():
            for s in M.free_symbols:
                num_hdict[(a, b)] = num_hdict[(a, b)].subs(s, sym_vals[s]).evalf()
        for (a, b), M in num_hdict.items():
            num_hdict[(a, b)] = np.array(M, dtype=complex)
        return num_hdict

# --- builders needed for producing an Ansatz object ---

def hermitian_block(a, b, d=3):
    hblock = sp.zeros(d, d)
    for i, j in it.product(range(d), repeat=2):
        if i==j:
            u_ii = sp.symbols(f'u_{i}{i}^{a}{b}', real=True)
            hblock[i, i] = u_ii
        elif i < j:
            u_ij, v_ij = sp.symbols([f'u_{i}{j}^{a}{b}', f'v_{i}{j}^{a}{b}'], real=True)
            hblock[i, j] = (u_ij + sp.I *v_ij)
            hblock[j, i] = (u_ij - sp.I *v_ij)
    return hblock

def KP_block(a, b, real=False, d=3):
    hblock = sp.zeros(d, d)
    u_00, u_11, u_22 = [sp.Symbol(f'u_{i}{i}^{a}{b}', real=True) for i in range(d)]
    hblock += sp.diag(u_00, u_11, u_22)
    for (i, j) in KP_DICT[(a, b)]:
        u_ij, v_ij = sp.symbols([f'u_{i}{j}^{a}{b}', f'v_{i}{j}^{a}{b}'], real=True)
        P_ij = sp.zeros(d, d)
        if real:
            P_ij[i, j] = u_ij
            P_ij[j, i] = u_ij
        else:
            P_ij[i, j] = u_ij + sp.I * v_ij
            P_ij[j, i] = u_ij - sp.I * v_ij
        hblock += P_ij
    return hblock


def block_from_labels(a, b, labels, prefix='A', basis=SPIN1):
    d, _ = basis[labels[0]].shape
    syms = {l: sp.Symbol(f'{prefix}^{a}{b}_{l}', real=True) for l in labels}
    M = sp.zeros(d, d)
    for l in labels:
        M += syms[l] * basis[l]
    return M

def rebuild_hermitian_block(a, b, d=3, subsyms={}):
    hblock = sp.zeros(d, d)
    for i, j in it.product(range(d), repeat=2):
        if i==j:
            k = ('u', a, b, i, i)
            if k in subsyms.keys():
                hblock[i, i] = subsyms[k]
        elif i < j:
            k1, k2 = ('u', a, b, i, j), ('v', a, b, i, j)
            if k1 in subsyms.keys():
                hblock[i, j] += subsyms[k1]
                hblock[j, i] += subsyms[k1]
            if k2 in subsyms.keys():
                hblock[i, j] += sp.I * subsyms[k2]
                hblock[j, i] += -sp.I * subsyms[k2]
    return hblock

def rebuild_block_from_labels(a, b, labels, prefix='A', basis=SPIN1, subsyms=()):
    d, _ = basis[labels[0]].shape
    syms = {l: subsyms[(a, b, l)] for l in labels}
    M = sp.zeros(d, d)
    for l in labels:
        M += syms[l] * basis[l]
    return M

def hdict_builder(family='full', reflection=False, pattern=None, prefix='A', basis_name='SPIN1', rebuild=False, subsyms=()):
    fam_labels = {'full': LABELS_FULL, 'real': LABELS_REAL, 'dipolar': LABELS, 'dipolar_real': LABELS_DIP_REAL}
    # this should be a dynamic variable that is based on the dimension of the underlying Hilbert space 
    if basis_name == 'SPIN1':
        basis = SPIN1
        d = 3
    elif basis_name == 'PAULI':
        basis = PAULI
        d = 2
    hdict = dict()
    for a, b in it.product(range(d), repeat=2):
        (ka, kb) = (min(a, b), max(a, b)) if reflection else (a, b)
        if (ka, kb) in hdict.keys():
            hdict[(a, b)] = hdict[(ka, kb)]
            continue
        active = 1 if pattern is None else pattern(a, b)
        if family == 'hermitian':
            M = hermitian_block(ka, kb, d)
            if rebuild:
                M = rebuild_hermitian_block(ka, kb, d, subsyms)
            # if active == 0, the pattern does not let the block contain off-diagonal parts
            if not active:
                M = sp.diag(*[M[i, i] for i in range(d)])
        elif family in ['KP', 'KP_real']:
            if family == 'KP': M = KP_block(ka, kb, real=False, d=d)
            elif family == 'KP_real': M = KP_block(ka, kb, real=True, d=d)
        else:
            labels = fam_labels[family]
            # labels = fam_labels[family]
            if not active:
                diag_labels = [l for l in labels if basis[l].is_diagonal()]
                labels = diag_labels
            M = block_from_labels(ka, kb, labels, prefix=prefix, basis=basis)
            if rebuild:
                M = rebuild_block_from_labels(ka, kb, labels, prefix=prefix, basis=basis, subsyms=subsyms)
        hdict[(ka, kb)] = M
        hdict[(a, b)] = M
    return hdict

def restrict(ansatz, linear_restr_eqs, name=None, verbose=True):
    sol = sp.solve(list(linear_restr_eqs), ansatz.syms, dict=True)
    nvars0 = ansatz.nvars
    if not sol:
        raise ValueError('inconsistent set of linear restrictions')
    sub = sol[0]
    hdict = {ab: sp.expand(M.subs(sub)) for ab, M in ansatz.hdict.items()}
    ans_family = Ansatz(hdict, d=ansatz.d, name=name or ansatz.name + '|restricted')
    nvars1 = ans_family.nvars
    if verbose:
        print(f'[NCq3Ansatz] number of symbols: {nvars0} --> {nvars1}')
    return ans_family
