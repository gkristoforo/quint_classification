# custom függvények az utils csomagból
from utils import greedy_solver

import sympy as sp
from sympy import eye, zeros, Matrix, symbols, Rational, init_printing
from IPython.display import display

# Pauli-algebra a szokásos bázison
s0 = eye(2)
sP = Matrix([[0, 1], [0, 0]]) 
sM = Matrix([[0, 0], [1, 0]])   
sZ = Matrix([[1, 0], [0, -1]])  
PAULI = {'0': s0, '+': sP, '-': sM, 'z': sZ}
LABELS = ['0', '+', '-', 'z']

# row, column count a bemenet, triviális eseteket skippeljük, 
# majd A celláiba betesszük A_ij * B-t, slice-okkal
def kron(A, B):
    ra, ca = A.shape
    rb, cb = B.shape
    M = zeros(ra * rb, ca * cb)
    for i in range(ra):
        for j in range(ca):
            if A[i, j] != 0:
                M[i * rb:(i + 1) * rb, j * cb:(j + 1) * cb] = A[i, j] * B
    return M

# mats listában tárolt mátrixok összetenzorszorzása
def kron_chain(mats):
    M = mats[0]
    for m in mats[1:]:
        M = kron(M, m)
    return M

# sűrűségek legyártása
# sigma^a az i-edik, sigma^b a j-edik helyre az L-hosszú láncon
# Hdict tárolja a Hamilton kifejtésében szereplő együtthatókat
def embed_pair(Hdict, i, j, L):
    dim = 2 ** L
    out = zeros(dim, dim)
    for (a, b), coeff in Hdict.items():
        if coeff == 0:
            continue
        # slots = [s0] * L
        slots = [s0.copy() for _ in range(L)]
        slots[i] = PAULI[a]
        slots[j] = PAULI[b]
        out += coeff * kron_chain(slots)
    return out

# teljes lánc felépítése a sűrűségekből
def Q2_periodic(Hdict, L):
    dim = 2 ** L
    Q = zeros(dim, dim)
    for n in range(L):
        Q += embed_pair(Hdict, n, (n + 1) % L, L)
    return Q

# használjuk a boostos formulát: Q_3 = [B(H), Q_2] -> q_3jj+1j+2 = [q_2jj+1,q_2j+1j+2]
def Q3_periodic(Hdict, L):
    dim = 2 ** L
    Q = zeros(dim, dim)
    for n in range(L):
        H1 = embed_pair(Hdict, n, (n + 1) % L, L)
        H2 = embed_pair(Hdict, (n + 1) % L, (n + 2) % L, L)
        Q += H1 * H2 - H2 * H1
    return Q

# ALGEBRAI EGYSZERŰSÍTÉS
# gauge-type I.: kinullázzuk a 00, ++, -- elemeket, a maradékot paraméterként beletesszük az ansatzba
# gauge-type I.: ha ++, akkor 1, ha --, akkor 0, valamint további feltételek az offdiagonális elemekre (amitől gII-es lesz egyáltalán)
def make_ansatz(gauge_type="I"):
    A, syms = {}, []
    # zero_terms = [('0', '0'), ('0', '+'), ('0', '-'), ('0', 'z'), ('+', '0'), ('-', '0'), ('z', '0')]
    zero_terms = [('0', '0')]
    for a in LABELS:
        for b in LABELS:
            if (a, b) in zero_terms:
                A[(a, b)] = 0
            elif gauge_type == "I" and (a, b) in [('+', '+'), ('-', '-')]:
                A[(a, b)] = 0
            elif gauge_type == "II" and (a, b) == ('+', '+'):
                A[(a, b)] = 1
            elif gauge_type == "II" and (a, b) == ('-', '-'):
                A[(a, b)] = 0
            else:
                s = symbols(f"A_{a}{b}")
                A[(a, b)] = s
                syms.append(s)
    if gauge_type == "II":
        syms = [s for s in syms if s not in (A[('z', '+')], A[('z', '-')], A[('z', 'z')])]

        A[('z', '+')] = -A[('+', 'z')]
        A[('z', '-')] = -A[('-', 'z')]
        A[('z', 'z')] = A[('+', '-')] + A[('-', '+')]
    return A, syms

# A cikkben leírt manuális feltételgyártás helyett legyártunk egy akkora láncot, amely már tartója egy site-on ható nemtriviális sűrűségnek 
# Hdict helyett az ansatzot írjuk be
# Mindenhol, ahol nem teljesül by default a [,]_ij = 0, teljesítendő feltételként hozzá adjuk a komponenst az egyenletrendszerhez
def generate_equations(A, syms, L=5):
    Q2, Q3 = Q2_periodic(A, L), Q3_periodic(A, L)
    comm = Q2 * Q3 - Q3 * Q2
    dim = 2 ** L
    eqs = set()
    for i in range(dim):
        for j in range(dim):
            e = sp.expand(comm[i, j])
            if e != 0:
                eqs.add(e)
    return list(eqs)

# a kapott egyenletrendszert osztályozzuk a monomiáljai szerint (egy monomiál mintha egy elsőrendű új változó lenne), majd a lineárisan függő kifejezéseket elhagyjuk
def independent_subset(eqs, syms):
    sorted_eqs = sorted(list(eqs), key=lambda x: (len(str(x)), str(x)))
    polys = [sp.Poly(e, *syms) for e in sorted_eqs]
    monoms = set()
    for pol in polys:
        for monom in pol.monoms():
            if sum([abs(i) for i in monom]) > 0:
                monoms.add(monom)
    mon_list = sorted(list(monoms), key=lambda x: (len(str(x)), str(x)))
    M = sp.zeros(len(polys), len(mon_list))
    for n, p in enumerate(polys):
        for mon_tup, coeff in p.as_dict().items():
            if mon_tup in mon_list:
                m = mon_list.index(mon_tup)
                M[n, m] = coeff
    _, pivots = M.T.rref()
    return [sorted_eqs[i] for i in pivots]

# valamilyen random numerikus A mátrixszal L hosszú próbaláncon verifikáció
def verify_solution(A_numeric, L_check=6):
    Q2, Q3 = Q2_periodic(A_numeric, L_check), Q3_periodic(A_numeric, L_check)
    comm = Q2 * Q3 - Q3 * Q2
    dim = 2 ** L_check
    return all(comm[i, j] == 0 for i in range(dim) for j in range(dim))

# actual kód
# Gauge type (megadni).-es ansatz
gauge_type = "I"
A, syms = make_ansatz(gauge_type)
print(f"Type {gauge_type} ansatz: {len(syms)} free parameters")
print(syms)
eqs = generate_equations(A, syms, L=5)
indep = independent_subset(eqs, syms)

# megoldás, elég fapados
print(f"raw nonzero equations: {len(eqs)}  ->  independent: {len(indep)}")
with open('2_site/eqs_out.txt', 'w') as file:
    for i in indep:
        string = '\\begin{dmath}' + str(sp.latex(i)) + '=0' + '\\end{dmath}' + '\n'
        file.write(string)
print("started solving the system of equations")
# solutions = sp.solve(eqs, syms, dict=True)


# NO EXTERNAL HAMILTONIAN CONTRIBUTION
# Mivel AxI - IxA jön a kommutátorból, ez teleszkopikusan kiesik: ezeket a tagokat nem kell nézni
for k in [('0', '+'), ('0', '-'), ('0', 'z'), ('+', '0'), ('-', '0'), ('z', '0')]:
    A[k] = 0
branch_syms = [A[k] for k in [('+', '-'), ('-', '+'), ('+', 'z'), ('-', 'z'), ('z', '+'), ('z', '-'), ('z', 'z')]]
branch_eqs = generate_equations(A, branch_syms, L=5)
solutions = sp.solve(branch_eqs, branch_syms, dict=True)
print(f"\ncase-split branch solved: {len(solutions)} solution families found")

# Verifikáció
sol = solutions[0]
free = [s for s in branch_syms if s in sol.values() or s not in sol]
numeric_map = {s: Rational([3, -2, 5, 1][i % 4]) for i, s in enumerate(branch_syms) if s not in sol}
A_num = {k: (v.subs(sol).subs(numeric_map) if hasattr(v, "subs") else v) for k, v in A.items()}
A_num = {k: (v.subs(numeric_map) if hasattr(v, "subs") else v) for k, v in A_num.items()}
ok = verify_solution(A_num, L_check=6)
print(f"solution family #1: {sol}")
print(f"independently re-verified at L=6: {ok}")
