import os, sys; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pathlib import Path
from NCq3solve.families import family, rebuild_family, symmetries_for
from NCq3solve.ansatz import hdict_builder
from NCq3solve.pipeline import Solver
from NCq3solve.operators import PROJECTORS
from flint import fmpq_mpoly
import sympy as sp
import itertools as it

name = sys.argv[1] if len(sys.argv) > 1 else 'Z3_RR'
htilde = sys.argv[2] if len(sys.argv) > 2 else None
if name in ['Z3', 'Z3_RR', 'hermitian']:
    ansatz_type = 'coordinate'
else:
    ansatz_type = 'spin'

s = Solver(family(name), htilde=htilde, symmetries=symmetries_for(name), verbose=True)
d = s.d
ctx = s.es.ctx
nametable = s.es.nametable
# s.generate()

# ---- check whether the symmetries actually apply ----
A = family('Z3_RR')
tau = sp.Matrix(3, 3, lambda i, j: 1 if i == (j+1) % 3 else 0)
# 1. ansatz level: the sum must be shift-invariant identically
M = sum((A.hdict[(a, b)] for a, b in it.product(range(3), repeat=2)), sp.zeros(3, 3))
assert sp.simplify(M - tau*M*tau.T) == sp.zeros(3, 3)
# 2. block level: h^{a+1,b+1} = tau h^{ab} tau^dag for every (a, b)
for a, b in it.product(range(3), repeat=2):
    assert sp.simplify(A.hdict[((a+1) % 3, (b+1) % 3)] - tau*A.hdict[(a, b)]*tau.T) == sp.zeros(3, 3)
# ---- end of check ----

# --- reading everything from the JSON file ---
import re
import sys
from pathlib import Path

MARKER = re.compile(r'"refined"\s*:\s*true')
def _find_matching_open(text: str, close_idx: int) -> int:
    """Walk backwards from a `]` and return the index of its balanced `[`."""
    depth = 0
    for i in range(close_idx, -1, -1):
        ch = text[i]
        if ch == "]":
            depth += 1
        elif ch == "[":
            depth -= 1
            if depth == 0:
                return i
    return -1

def extract_refined_blocks(path, include_brackets: bool = False, nested: bool = False, encoding: str = "utf-8"):
    text = Path(path).read_text(encoding=encoding)
    blocks = []
    seen_spans = set()
    for match in MARKER.finditer(text):
        close_idx = text.rfind("]", 0, match.start())
        if close_idx == -1:
            continue
        if nested:
            open_idx = _find_matching_open(text, close_idx)
        else:
            open_idx = text.rfind("[", 0, close_idx)
        if open_idx == -1:
            continue
        span = (open_idx, close_idx)
        if span in seen_spans:
            continue
        seen_spans.add(span)
        if include_brackets:
            blocks.append(text[open_idx : close_idx + 1])
        else:
            blocks.append(text[open_idx + 1 : close_idx])
    return blocks
# --- end of read JSON ---

def txt_to_flint(eq, ctx=ctx):
    s = str(eq)
    s = s.replace('"', '')
    s = s.replace(',', '')
    a = s.count('\t')
    while a > 0:
        s = s.replace('\t', '')
    poly = fmpq_mpoly(s, ctx)
    return poly

def flint_to_sympy(eqflint, nametable=nametable):
    eq, ctx = eqflint, eqflint.context()
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

# --- LINEARIZING EXPRESSIONS ---
def linearize_expr(expr, prefix='m', syms_dict=None):
    if syms_dict is None:
        syms_dict = dict()
    reverse_syms_dict = {v: k for k, v in syms_dict.items()}
    expanded = sp.expand(expr)
    if isinstance(expanded, sp.Add):
        terms = expanded.args
    else:
        terms = [expanded]
    out = 0
    for term in terms:
        coeff, subexpr = term.as_coeff_Mul()
        subexpr = sp.cancel(subexpr)
        if subexpr == 1:
            out += coeff
            continue
        if subexpr not in reverse_syms_dict:
            s = sp.symbols(f'{prefix}_{len(syms_dict) + 1}')
            syms_dict[s] = subexpr
            reverse_syms_dict[subexpr] = s
            out += coeff * s
        else:
            out += coeff * reverse_syms_dict[subexpr]
    return out, syms_dict

def linearize_matrix(A, syms_dict=None):
    r, c = A.shape
    A_lin = sp.zeros(r, c)
    if syms_dict is None:
        syms_dict = dict()
    for i in range(r):
        for j in range(c):
            lin_term, syms_dict = linearize_expr(A[i, j], prefix='m', syms_dict=syms_dict)
            A_lin[i, j] = lin_term
    return A_lin, syms_dict

def kron(X, Y):
    rx, cx = X.shape
    ry, cy = Y.shape
    M = sp.zeros(rx * ry, cx * cy)
    for i in range(0, rx):
        for j in range(0, cx):
            M[i*ry:(i+1)*ry, j*cy:(j+1)*cy] = X[i, j] * Y
    return M

def create_chain(mats):
    M = mats[0]
    for i in range(1, len(mats)):
        M = kron(M, mats[i])
    return M

# ---- END OF DEFINITIONS ----

ctx = s.es.ctx
nametable = s.es.nametable
symbols = {sf: sp.Symbol(f'{nametable[sf]}', real=True) for sf in nametable.keys()}
txt = []

gen_strs = extract_refined_blocks(f'/Users/kristof/main/research/integrable_systems/NCq3_calculations/results/{name}_{htilde}.json',
                                  include_brackets=False, nested=False)

for gen in gen_strs:
    txt_i = []
    with open('temp.txt', 'w') as file:
        file.write(gen)
    with open('temp.txt', 'r') as file:
        for line in file:
            str_to_app = line.strip()
            if str_to_app not in ['', '\n']:
                txt_i += [str_to_app]
    os.remove('temp.txt')
    txt.append(txt_i)

# with open(f'/Users/kristof/main/research/integrable_systems/NCq3_calculations/generators/{name}_comps.txt', 'r') as file:
    # for line in file:
    #     txt += [line.strip()]

sympy_eqs = []
for comp in txt:
    eq_i = []
    for eqstr in comp:
        eqsp = flint_to_sympy(txt_to_flint(eqstr))
        eq_i += [eqsp]
    sympy_eqs.append(eq_i)
print(f'all sympy equations are obtained')

for i, comp in enumerate(sympy_eqs):
    solution = sp.solve(comp, *[s for s in symbols.values()], dict=True)
    if type(solution) == list:
        solution = solution[0]
    print(f'component {i} is solved')

    free_syms = {}
    for s in symbols.values():
        if s not in solution.keys():
            solution[s] = s
            free_syms[s] = s

    hdict = A.hdict
    hdict_subs = {}
    print(f'solution: {solution}')
    index_solsyms = {}

    # if ansatz_type == 'coordinate':
    #     index_solsyms = {(str(k)[0], int(str(k)[5]), int(str(k)[6]), int(str(k)[2]), int(str(k)[3])): s for k, s in solution.items()}
    # elif ansatz_type == 'spin':
    #     index_solsyms = {(int(str(k)[2]), int(str(k)[3]), str(k)[5:]): s for k, s in solution.items()}

    for (a, b), M in hdict.items():
        M_ab = sp.zeros(d, d)
        for i, j in it.product(range(d), repeat=2):
            M_ab[i, j] = M[i, j].subs([(k, s) for k, s in solution.items()])
        hdict_subs[(a, b)] = M_ab
    print(f'hdict_{i}[0] before: {hdict[[k for k in hdict.keys()][0]]}')
    print(f'hdict_{i}[0] after: {hdict_subs[[k for k in hdict.keys()][0]]}')
    hblock = sp.zeros(d**3, d**3)
    with open(f'/Users/kristof/main/research/integrable_systems/NCq3_calculations/output_viewer/{name}_solutions.txt', 'a') as file:
        file.write(f'Genuine Solution no. {i} for family {name}: \n')
        file.write('\\begin{gather} \n')
        for (a, b), M in hdict_subs.items():
            string = 'h_' + str(i) + '^{' + str(a) + str(b) + '}=' + str(sp.latex(M)) + '\\\\ \n'
            file.write(string)
        file.write('\\end{gather} \n')

# hdict = {}
# if ansatz_type == 'coordinate':
#     index_solsyms = {(str(k)[0], int(str(k)[5]), int(str(k)[6]), int(str(k)[2]), int(str(k)[3])): s for k, s in solution.items()}
#     hdict = rebuild_family(name, index_solsyms).hdict
# elif ansatz_type == 'spin':
#     index_solsyms = {(int(str(k)[2]), int(str(k)[3]), str(k)[5:]): s for k, s in solution.items()}
#     hdict = rebuild_family(name, index_solsyms).hdict
# else:
#     raise ValueError(f'not proper ansatz type (coordinate/spin)')

"""
hblock = sp.zeros(d**3, d**3)
with open(f'/Users/kristof/main/research/integrable_systems/NCq3_calculations/output_viewer/{name}_9blocks.txt', 'w') as file:
    for (a, b), M in hdict_subs.items():
        string = 'h^{' + str(a) + str(b) + '}=' + str(sp.latex(M)) + '\\\\ \n'
        hblock += create_chain([PROJECTORS[a], M, PROJECTORS[b]])
        file.write(string)
        # print(string, sp.latex(M), '\\\\')

# --- Full h_j determination ---

hblock, newsyms = linearize_matrix(hblock, free_syms)

print(f'nametable: {nametable}')
# print('H_{\mathrm{gen}=10, \mathrm{lin}}=', sp.latex(hblock), '\\\\')
"""