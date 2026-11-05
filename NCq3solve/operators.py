import sympy as sp
import itertools as it

"""
 ---- OPERATOR STRUCTURES ----
 The following operators need to be built into the program, divided into separate categories.
 Spin-1:
 1. Dipolar subspace: {S0, Sx, Sy, Sz}
 2. Real subspace: {S0, Sx, Sz, (Sx)^2, (Sz)^2, {Sx,Sz} }
 3. Full Hermitian subspace: {S0, Sx, Sy, Sz, (Sx)^2, (Sz)^2, {Sx,Sz}, {Sx,Sy}, {Sy,Sz}}
"""
im = sp.I
S0 = sp.eye(3)
Sx = sp.Matrix([[0, 1, 0,], [1, 0, 1], [0, 1, 0]])
Sy = im * sp.Matrix([[0, -1, 0,], [1, 0, -1], [0, 1, 0]])
Sz = sp.Matrix([[1, 0, 0,], [0, 0, 0], [0, 0, -1]])
Sp = Sx + im * Sy
Sm = Sx - im * Sy

# useful spin-1 operators and their combinations (double letters mean anticommutators)
SPIN1 = {'0': S0, 'y': Sy, 'x': Sx, 'z': Sz, '+': Sp, '-': Sm,
         'xx': Sx*Sx, 'zz': Sz*Sz, 'xz': Sx*Sz + Sz*Sx,
         'xy': Sx*Sy + Sy*Sx, 'yz': Sy*Sz + Sz*Sy}

# projectors onto the different basis states, augmented with offdiagonals
# PROJECTORS = dict()
# for i, j in it.product(range(3), repeat=2):
#     mat = sp.zeros(3, 3)
#     mat[i, j] = 1
#     PROJECTORS[(i, j)] = mat

# projectors onto the different basis states
PROJECTORS = dict()
for i in range(3):
    diags = [0, 0, 0]
    diags[i] = 1
    PROJECTORS[i] = sp.diag(*diags)

# diads for the Klobas-Prosen family
KP_DICT = {}
for a, b in it.product(range(3), repeat=2):
    if (a, b) in [(0, 0)]:
        KP_DICT[(a, b)] = [(0, 1), (0, 2), (1, 2)]
    elif (a, b) in [(0, 1), (1, 0), (1, 1)]:
        KP_DICT[(a, b)] = [(0, 1)]
    elif (a, b) in [(0, 2), (2, 0), (2, 2)]:
        KP_DICT[(a, b)] = [(0, 2)]
    else:
        KP_DICT[(a, b)] = [(1, 2)]


# these are the lists over which iterations will be carried out
LABELS = ['0', 'x', 'y', 'z']
LABELS_DIP_REAL = ['0', 'x', 'z']
LABELS_REAL = ['0', 'x', 'z', 'xx', 'zz', 'xz']
LABELS_FULL = ['0', 'x', 'y', 'z', 'xx', 'zz', 'xz', 'xy', 'yz']

# Pauli matrices for spin-1/2 validations
# this uses the same labels as the default spin-1 operators
PAULI = {'0': sp.eye(2), 'x': sp.Matrix([[0, 1], [1, 0]]), 'y': im * sp.Matrix([[0, -1], [1, 0]]), 'z': sp.Matrix([[1, 0], [0, -1]])}
