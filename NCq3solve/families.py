from .ansatz import Ansatz, hdict_builder, restrict
import itertools as it
import sympy as sp
import numpy as np

M_OF = [1, 0, -1]
BAR = [2, 1, 0]

"""
 The following families can be chosen for the classification:
 1. A, C, D: hardest families with no constraint
 2. B, B_R, B_RZ2, Bp: dipolar families
 3. E, E_B: spin-flip invariant models
 4. Z3, Z3_RR
"""

def family(name):
    extra = []
    extra_for_solver = []
    if name == 'A':         hd = hdict_builder('hermitian')
    elif name == 'B':       hd = hdict_builder('dipolar_real')
    elif name == 'Bp':      hd = hdict_builder('dipolar')
    elif name == 'B_R':     hd = hdict_builder('dipolar_real', reflection=True)
    elif name == 'B_RZ2':
        hd = hdict_builder('dipolar_real', reflection=True)
        def sym(a, b, l):
            return sp.Symbol(f'A^{min(a, b)}{max(a, b)}_{l}', real=True)
        can_keys = set()
        for a, b in it.product(range(3), repeat=2):
            a2, b2 = BAR[a], BAR[b]
            k1, k2 = (min(a, b), max(a, b)), (min(a2, b2), max(a2, b2))
            if (k1 in can_keys) or (k2 in can_keys):
                continue
            can_keys.add(k1)
            extra += [sym(a, b, '0') - sym(a2, b2, '0'), sym(a, b, 'x') - sym(a2, b2, 'x'), sym(a, b, 'z') + sym(a2, b2, 'z')]
    elif name == 'C':       hd = hdict_builder('real')
    elif name == 'D':       hd = hdict_builder('real', reflection=True)
    elif name in ['KP', 'KP_real']:
        if name == 'KP': hd = hdict_builder('KP')
        else:            hd = hdict_builder('KP_real')
        for (a, b, i) in [(0, 0, 0), (1, 0, 0), (2, 0, 0), (1, 0, 2), (2, 0, 1)]:
            s = sp.Matrix(hd[(a, b)])[i, i]
            assert isinstance(s, sp.Symbol), f'u^{a}{b}_{i}{i} must be a bare symbol'
            extra_for_solver.append(s)
    elif name == 'E':
        def pattern_E(a, b): return int(M_OF[a] + M_OF[b] == 0)
        hd = hdict_builder('real', pattern=pattern_E)
    elif name == 'E_B':
        def pattern_E(a, b): return int(M_OF[a] + M_OF[b] == 0)
        hd = hdict_builder('dipolar_real', pattern=pattern_E)
    elif name == 'Z3':
        # "a" can be fixed to 0, changing it to any other value would yield redundant equations
        hd = hdict_builder('hermitian')
        a = 0
        for b in range(3):
            M = hd[(a, b)]
            for k in (1, 2):
                ak, bk = (a+k)%3, (b+k)%3
                Mk = hd[ak, bk]
                for i, j in it.product(range(3), repeat=2):
                    e = sp.expand(Mk[(i+k)%3, (j+k)%3] - M[i, j])
                    extra += [sp.re(e), sp.im(e)]
    elif name == 'Z3_RR':
        spin = False
        if spin: hd = hdict_builder('real', reflection=True)
        else:
            hd = hdict_builder('hermitian', reflection=True)
            extra += [sp.im(e) for M in hd.values() for e in M]
        a = 0
        for b in range(3):
            M = hd[(a, b)]
            for k in (1, 2):
                Mk = hd[(a+k)%3, (b+k)%3]
                for i, j in it.product(range(3), repeat=2):
                    e = sp.expand(Mk[(i+k)%3, (j+k)%3] - M[i, j])
                    extra += [sp.re(e), sp.im(e)]
    else:
        raise NameError(f'symmetry family {name} is not defined')
    extra = [e for e in extra if e != 0]
    ans = Ansatz(hd, name=name)
    return (restrict(ans, extra, name if not extra else name + '|restricted') if extra else ans), extra_for_solver

def reflection(d):
    def refl(a, b): return (b, a)
    return dict(name='reflection', block=refl, mid=None)

def relabel(perm):
    def permute(a, b): return (perm[a], perm[b])
    d = len(perm)
    Pm = np.zeros((d, d))
    for c in range(d):
        Pm[perm[c], c] = 1
    return dict(name=f'relabel{tuple(perm)}', block=permute, mid=Pm)

def symmetries_for(name):
    if name in ['A', 'C', 'D']:
        symmetries = [reflection(3), *[relabel(list(p)) for p in it.permutations(range(3)) if list(p) != [0, 1, 2]]]
    elif name in ['Z3', 'Z3_RR']:
        symmetries = [reflection(3), relabel([0, 2, 1])]
    elif name in ['B', 'Bp', 'B_R', 'B_RZ2', 'E', 'E_B']:
        symmetries = [reflection(3), relabel([2, 1, 0])]
    elif name in ['KP', 'KP_real']:
        symmetries = [reflection(3), relabel([0, 2, 1])]
    else:
        raise NameError(f'symmetry family {name} is not defined')
    return symmetries


# --- REBUILD SOLVED FAMILIES ---
def rebuild_family(name, subsyms=()):
    extra = []
    if name == 'A':         hd = hdict_builder('hermitian', rebuild=True, subsyms=subsyms)
    elif name == 'B':       hd = hdict_builder('dipolar_real', rebuild=True, subsyms=subsyms)
    elif name == 'Bp':      hd = hdict_builder('dipolar', rebuild=True, subsyms=subsyms)
    elif name == 'B_R':     hd = hdict_builder('dipolar_real', reflection=True, rebuild=True, subsyms=subsyms)
    elif name == 'B_RZ2':
        hd = hdict_builder('dipolar_real', reflection=True, rebuild=True, subsyms=subsyms)
        def sym(a, b, l):
            return sp.Symbol(f'A^{min(a, b)}{max(a, b)}_{l}', real=True)
        can_keys = set()
        for a, b in it.product(range(3), repeat=2):
            a2, b2 = BAR[a], BAR[b]
            k1, k2 = (min(a, b), max(a, b)), (min(a2, b2), max(a2, b2))
            if (k1 in can_keys) or (k2 in can_keys):
                continue
            can_keys.add(k1)
            extra += [sym(a, b, '0') - sym(a2, b2, '0'), sym(a, b, 'x') - sym(a2, b2, 'x'), sym(a, b, 'z') + sym(a2, b2, 'z')]
    elif name == 'C':       hd = hdict_builder('real', rebuild=True, subsyms=subsyms)
    elif name == 'D':       hd = hdict_builder('real', reflection=True, rebuild=True, subsyms=subsyms)
    elif name == 'E':
        def pattern_E(a, b): return int(M_OF[a] + M_OF[b] == 0)
        hd = hdict_builder('real', pattern=pattern_E, rebuild=True, subsyms=subsyms)
    elif name == 'E_B':
        def pattern_E(a, b): return int(M_OF[a] + M_OF[b] == 0)
        hd = hdict_builder('dipolar_real', pattern=pattern_E, rebuild=True, subsyms=subsyms)
    elif name == 'Z3':
        # "a" can be fixed to 0, changing it to any other value would yield redundant equations
        hd = hdict_builder('hermitian', rebuild=True, subsyms=subsyms)
        a = 0
        for b in range(3):
            M = hd[(a, b)]
            for k in (1, 2):
                ak, bk = (a+k)%3, (b+k)%3
                Mk = hd[ak, bk]
                for i, j in it.product(range(3), repeat=2):
                    e = sp.expand(Mk[(i+k)%3, (j+k)%3] - M[i, j])
                    extra += [sp.re(e), sp.im(e)]
    elif name == 'Z3_RR':
        spin = False
        if spin: hd = hdict_builder('real', reflection=True, rebuild=True, subsyms=subsyms)
        else:
            hd = hdict_builder('hermitian', reflection=True, rebuild=True, subsyms=subsyms)
            extra += [sp.im(e) for M in hd.values() for e in M]
        a = 0
        for b in range(3):
            M = hd[(a, b)]
            for k in (1, 2):
                Mk = hd[(a+k)%3, (b+k)%3]
                for i, j in it.product(range(3), repeat=2):
                    e = sp.expand(Mk[(i+k)%3, (j+k)%3] - M[i, j])
                    extra += [sp.re(e), sp.im(e)]
    else:
        raise NameError(f'symmetry family {name} is not defined')
    extra = [e for e in extra if e != 0]
    ans = Ansatz(hd, name=name, rebuild=True)
    return ans
