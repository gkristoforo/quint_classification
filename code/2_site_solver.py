import numpy as np
import sympy as sp
from math import prod

S0 = sp.eye(2)
Sp = sp.Matrix([[0, 1], [0, 0]])
Sm = sp.Matrix([[0, 0], [1, 0]])
Sz = sp.Matrix([[1, 0], [0, -1]])
Sx = Sm + Sp

PAULI = {'0': S0, '+': Sp, '-': Sm, 'z': Sz, 'x': Sx}
LABELS = ['0', '+', '-', 'z']

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

def id_Llist(L):
    return [S0 for _ in range(L)]

def id_Lchain(L):
    return create_chain(id_Llist(L))

def embed_pair(density, i, j, L):
    dim = 2**L
    out_M = sp.zeros(dim, dim)
    for (a, b), coeff in density.items():
        if coeff == 0:
            continue
        slots = [S0 for _ in range(L)]
        slots[i] = PAULI[a]
        slots[j] = PAULI[b]
        M = coeff * create_chain(slots)
        # display(M)
        out_M += M
    return out_M

# 2 site Hamilton-sűrűségre ansatz
def ansatz_2site(gauge_type='I', ext_fields=True):
    syms = list()
    density = dict()
    zero_terms = [('0', '0')]
    if ext_fields == False:
        zero_terms = [('0', '0'), ('0', '+'), ('0', '-'), ('0', 'z'), ('+', '0'), ('-', '0'), ('z', '0')]
    for a in LABELS:
        for b in LABELS:
            if (a, b) in zero_terms:
                density[(a, b)] = 0
                continue
            elif gauge_type == 'I' and (a, b) in [('+', '+'), ('-', '-')]:
                density[(a, b)] = 0
            elif gauge_type == 'II' and (a, b) == ('+', '+'):
                density[(a, b)] = 1
            elif gauge_type == 'II' and (a, b) == ('-', '-'):
                density[(a, b)] = 0
            else:
                s = sp.symbols(f'A_{a}{b}')
                syms.append(s)
                density[(a, b)] = s
    if gauge_type == 'II':
        syms = [s for s in syms if s not in (density[('z', '+')], density[('z', '-')], density[('z', 'z')])]
        density[('z', '+')] = -density[('+', 'z')]
        density[('z', '-')] = -density[('-', 'z')]
        density[('z', 'z')] = density[('+', '-')] + density[('-', '+')]
    return density, syms

# density = {('pauli_a', 'pauli_b'): coeff}
def Q2_periodic(density, L):
    dim = 2 ** L
    if L < 2:
        raise Exception('chain is too short, at least 2 sites are needed')
    H = sp.zeros(dim, dim)
    # print(density)
    for i in range(L):
        hi = embed_pair(density, i, (i+1)%L, L)
        # display(hi)
        H += hi
    return H

def Q3_periodic(density, L):
    dim = 2**L
    Q3 = sp.zeros(dim, dim)
    for i in range(L):
        hi = embed_pair(density, i, (i+1)%L, L)
        hip1 = embed_pair(density, (i+1)%L, (i+2)%L, L)
        q3i = hi * hip1 - hip1 * hi
        Q3 += q3i
    return Q3

def eq_generator(density, L):
    dim = 2**L
    Q2 = Q2_periodic(density, L)
    Q3 = Q3_periodic(density, L)
    comm = Q2 * Q3 - Q3 * Q2
    # display(comm)
    eqs = set()
    for i in range(dim):
        for j in range(dim):
            e = sp.expand(comm[i, j])
            if e != 0:
                eqs.add(e)
    return list(eqs)

def eq_linear_eliminator(eqs, syms):
    eq_list = sorted(list(eqs), key=lambda x: (len(str(x)), str(x)), reverse=False)
    polys = [sp.Poly(e, *syms) for e in eq_list if sp.Poly(e, *syms) != 0]
    # print(len(polys))
    monoms_set = set()
    for p in polys:
        for monom in p.monoms():
            if sum([i for i in monom]) > 0:
                # expr = sp.Mul(*[g**e for g, e in zip(p.gens, term[0])])
                monoms_set.add(monom)
    monoms = sorted(list(monoms_set), key=lambda x: (len(str(x))), reverse=False)
    M = sp.zeros(len(polys), len(monoms))
    # print(M.shape)
    for i in range(len(polys)):
        p = polys[i]
        for j in range(len(monoms)):
            m = monoms[j]
            M[i, j] = p.coeff_monomial(m)
    _, pivots = M.T.rref()
    return [eq_list[i] for i in pivots if eq_list[i] != 0]


# --- systematically choosing independent subset of Hamiltonians ---
def linearize_expr(expr, prefix='m', syms_dict=None):
    if syms_dict is None:
        syms_dict = dict()
    reverse_syms_dict = {v: k for k, v in syms_dict.items()}
    expanded = sp.expand(expr)
    if isinstance(expanded, sp.Add):
        terms = expanded.args
    else:
        terms = [expanded]
    # linear_subexprs = []
    out = 0
    for term in terms:
        coeff, subexpr = term.as_coeff_Mul()
        subexpr = sp.cancel(subexpr)
        # print(subexpr, coeff)
        if subexpr == 1:
            out += coeff
            continue
        if subexpr not in reverse_syms_dict:
            # linear_subexprs.append(subexpr)
            s = sp.symbols(f'{prefix}_{len(syms_dict) + 1}')
            # m, i = prefix, str(counter)
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
            # print(lin_term, syms_dict)
            # print(f'ok, i={i}, j={j}')
    return A_lin, syms_dict

def construct_basis_matrix(A, syms=[]):
    r, c = A.shape
    out = []
    syms_dict = {s: s for s in syms}
    A_lin, syms_dict = linearize_matrix(A, syms_dict=syms_dict)
    syms = [a for a, b in syms_dict.items()]
    row = len(syms)
    for i in range(row):
        s = syms[i]
        pds = sp.diff(A_lin, s)
        if pds != sp.zeros(r, c):
            out.append(list(pds))
    return sp.Matrix(out) if out else sp.Matrix([])

def check_independence(A, B):
    Wa = construct_basis_matrix(A)
    Wb = construct_basis_matrix(B)
    rankA = Wa.rank()
    rankB = Wb.rank()
    # indep típusai (az erősebb marad a listában, vagy mindkettő, ha van független alterük):
    # 0 - A mátrix általánosabb, vagy ugyanolyan általánosak
    # 1 - B mátrix általánosabb
    # 2 - mindkét mátrix feszít ki független megoldást
    indep = 0
    if rankA == 0 or rankB == 0:
        # print(f'Az egyik mátrix a nullmátrix')
        pass
    else:
        Wab = sp.Matrix.vstack(Wa, Wb)
        rankAB = Wab.rank()
        if rankA == rankAB and rankB == rankAB:
            # print(f'A két mátrix ugyanazt a megoldásteret feszíti ki')
            pass
        elif rankA == rankAB and rankB < rankAB:
            # display(B)
            # print(f'függ H_1-től.')
            pass
        elif rankB == rankAB and rankA < rankAB:
            indep = 1
            # print(f'H_1 függ az alábbi mátrixtól:')
            # display(B)
        elif rankAB == rankA + rankB:
            indep = 2
            # print(f'A két mátrix teljesen független')
        else:
            indep = 2
            # print(f'A két mátrix független, de van közös alterük')
    return indep

"""
--- Gauge fixing relations ---
    The ansatz_2site function relies on the code below for using 
    Gauge types I. and II., but the relations are obtained by
    looking at the most general case, so it is not going to cause issues 
    that the equations are obtained later than the ansatz is being made.
"""
alpha, beta, gamma, delta = sp.symbols(['alpha', 'beta', 'gamma', 'delta'])

# the function takes the fact granted that the transformation matrix V has unit determinant
def embed_transformed_pair(density, i, j, L, V=None):
    dim = 2**L
    out_M = sp.zeros(dim, dim)
    for (a, b), coeff in density.items():
        if coeff == 0:
            continue
        slots = [S0 for _ in range(L)]
        if V != None:
            alpha, beta, gamma, delta = V[0, 0], V[0, 1], V[1, 0], V[1, 1]
            Vinv = sp.Matrix([[delta, -beta], [-gamma, alpha]])
            slots[i] = V * PAULI[a] * Vinv
            slots[j] = V * PAULI[b] * Vinv
        else:
            slots[i] = PAULI[a]
            slots[j] = PAULI[b]
        M = coeff * create_chain(slots)
        # display(M)
        out_M += M
    return out_M

h_dict, syms = ansatz_2site('not type I')
V = sp.Matrix([[alpha, beta], [gamma, delta]])
h = embed_transformed_pair(h_dict, 0, 1, 2, V=None)
h_prime = embed_transformed_pair(h_dict, 0, 1, 2, V=V)
# display(h)
# display(h_prime)

"""
This sadly works only in a .ipynb file:

print('Gauge fixing relations to be set to zero:')
print(f'The term corresponding to A\'_++:')
display(sp.expand(h_prime[0, 3]))
print(f'The term corresponding to A\'_--:')
display(sp.expand(h_prime[3, 0]))
"""

"""
---Solving the case of No External Coupling (NEC) case---
    Every interaction with the identity as one of the two site Pauli operators is
    excluded from the Hamiltonian (this would correspond to single site magnetic fields
    coupled to certain spin directions).
"""
L=5
h, syms = ansatz_2site(ext_fields=False)
print(syms)
eqs = eq_generator(h, L)
indep = sorted(eq_linear_eliminator(eqs, syms), key=lambda x: (len(str(x)), str(x)), reverse=False)

print(f"raw nonzero equations: {len(eqs)}  ->  independent: {len(indep)}")
with open('eqs_out2.txt', 'w') as file:
    for i in indep:
        string = '\\begin{dmath}' + str(sp.latex(i)) + '=0' + '\\end{dmath}' + '\n'
        file.write(string)
print("finished writing the file")
solution = sp.solve(indep, syms, dict=True)
print("finished solving the NEC case\n number of solutions: ", len(solution))

# ---analyzing the solution of No External Coupling (NEC) class---
n = len(solution)
print(f'number of solutions: {n}')
H_sols = list()
count = 1
for sol in solution:
    H = sp.zeros(4, 4)
    sol_syms = list(sol.keys())
    for s in sol_syms:
        label1 = str(s)[-2]
        label2 = str(s)[-1]
        site1 = PAULI[label1]
        site2 = PAULI[label2]
        hs = sol[s] * kron(site1, site2)
        H += hs
    # print(f'Hamiltonian of solution family number {count}:')
    # display(H)
    # print(sp.latex(H))
    H_sols.append(H)
    count += 1

n = len(H_sols)
indices = [i for i in range(0, n)]
i = 0
l = len(indices)
while i < l:
    index = indices[i]
    A = H_sols[index]
    j = i+1
    while j < l:
    # for j in range(i+1, l):
        indep = check_independence(A, H_sols[indices[j]])
        print(indices, i, j, indices[j], indep)
        if indep == 0:
            print(f'{index} kiütötte a(z) {indices[j]} mátrixot')
            indices.pop(j)
        elif indep == 1:
            indices.pop(i)
            break
        else:
            j += 1
        l = len(indices)
        print(f'len: {l}')
    i += 1
print(f'A végső indexek: {indices}')

with open('NEC_Hamiltonians.txt', 'w') as file:
    for i in range(len(indices)):
        H = H_sols[indices[i]]
        string = 'Solution family no. ' + str(i) + ': \n' + '\\begin{equation} H_{' + str(i) +'}=' + str(sp.latex(H)) + '\\end{equation}' + '\n'
        file.write(string)

