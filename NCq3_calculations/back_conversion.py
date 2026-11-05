import re, os, time, json, sys; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pathlib import Path
from NCq3solve.families import family, symmetries_for
from NCq3solve.pipeline import Solver
from flint import fmpq_mpoly
import sympy as sp, numpy as np
import itertools as it

name = sys.argv[1] if len(sys.argv) > 1 else 'KP_real'
htilde = sys.argv[2] if len(sys.argv) > 2 else None
if name in ['Z3', 'Z3_RR', 'hermitian', 'KP', 'KP_real']:
    ansatz_type = 'coordinate'
else:
    ansatz_type = 'spin'

A, slice_eqs = family('KP_real')

S = Solver(A, htilde=htilde, symmetries=symmetries_for(name), extra_equations=slice_eqs, verbose=True)
d = S.d
ctx = S.es.ctx
nametable = S.es.nametable
# syms_ans = A.syms
symbols = {sf: sp.Symbol(f'{nametable[sf]}', real=True) for sf in nametable.keys()}

sym_vals = {s: np.random.uniform(-10, 10) for s in symbols.values()}

# ---- DEFINITIONS ----

# --- txt to sympy flow with an interstitial flint step ---
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

def timer(t0: float):
    """Print elapsed time since t0 (a time.time() timestamp) until Ctrl+C."""
    try:
        while True:
            h, rem = divmod(int(time.time() - t0), 3600)
            m, s = divmod(rem, 60)
            print(f"\rElapsed: {h:02d}:{m:02d}:{s:02d}", end="", flush=True)
            time.sleep(0.2)
    except KeyboardInterrupt:
        print()

# ---- END OF DEFINITIONS ----

with open(f'/Users/kristof/main/research/integrable_systems/NCq3_calculations/results/{name}_{htilde}.json', 'r') as file:
    data = json.load(file)

eqs_extracted = {}
for pattern, v in data.items():
    eqs_p = []
    for i, comp in enumerate(v['components']):
        if comp['refined']:
            eqs_extracted[(pattern, i)] = comp['generators']
sympy_eqs = {}
for key, eqs in eqs_extracted.items():
    eq_i = []
    for eqstr in eqs:
        eqsp = flint_to_sympy(txt_to_flint(eqstr))
        eq_i += [eqsp]
    sympy_eqs[key] = eq_i
print(f'all sympy equations are obtained')
print(f'symbols: {symbols}')

# ind_to_skip = [2, 3, 4] # Z3_RR, htilde = None
ind_to_skip = [2]

with open(f'/Users/kristof/main/research/integrable_systems/NCq3_calculations/output_viewer/{name}_{htilde}_solutions.txt', 'w+') as file:
    for n, (key, comp) in enumerate(sympy_eqs.items()):
        if n in ind_to_skip:
            continue
        # timer(time.time()) 
        solution = sp.solve(comp, *[s for s in symbols.values()], dict=True)
        if len(solution) == 0:
            print(f'component {n}.: {key} has no solutions, comp and equations:\n {comp}')
            continue
        if type(solution) == list:
            solution = solution[0]
            print
        print(f'component {n} is solved')

        free_syms = {}
        for s in symbols.values():
            if s not in solution.keys():
                solution[s] = s
                free_syms[s] = s

        hdict = A.hdict
        hdict_subs = {}
        # print(f'solution: {solution}')
        print(f'free symbols: {free_syms}')

        for (a, b), M in hdict.items():
            M_ab = sp.zeros(d, d)
            for i, j in it.product(range(d), repeat=2):
                M_ab[i, j] = M[i, j].subs([(k, s) for k, s in solution.items()])
            hdict_subs[(a, b)] = M_ab

        hblock = sp.zeros(d**3, d**3)
        file.write(f'Genuine Solution no. {n} for family {name}_{htilde} with pattern {key[0]}, comp {key[1]}: \n')
        file.write('\\begin{gather*} \n')
        for (a, b), M in hdict_subs.items():
            string = '\t h_{' + str(n) + '}^{' + str(a) + str(b) + '}=' + str(sp.latex(M, mat_str="pmatrix", mat_delim="")) + '\\\\ \n'
            file.write(string)
        file.write('\\end{gather*} \n')
        
        if htilde == 'diag-gauge':
            aux_dict = S.es.aux_dict
            aux_subs = {}
            for k, (vre, vim) in aux_dict.items():
                vre_sp = flint_to_sympy(vre)
                vsub = vre_sp.subs([(key, s) for key, s in solution.items()])
                aux_subs[k] = vsub
            file.write(f'htilde entries for Solution no. {n} in family {name} {htilde} with pattern {key[0]}, comp {key[1]}: \n')
            file.write('\\begin{gather*} \n')
            for (ac, ac), v in aux_subs.items():
                string = '\t \\gamma^{' + str(ac) + '}_{' + str(n) + '}=' + str(sp.latex(v, mat_str="pmatrix", mat_delim="")) + '\\\\ \n'
                file.write(string)
            file.write('\\end{gather*} \n')