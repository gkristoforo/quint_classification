KÉRDÉS:
Az arxiv:2108.02053 cikkben, a III.53. egyenletben látható ansatz komplexitásában nagyon hasonlít a 2-site modellben láthatóra. Éppen azon dolgozom, hogy egy 16 paraméteres halmazt belátható időn belül, Gröbner-bázisos, polinomgyűrűs módszerrel meg tudjak oldani. Érdemes lehet ezt jobban megvizsgálni, hátha van egy bővebb paramétertér is, amit még klasszifikálni tudnék? Szerintem itt vannak érdekes dolgok.

---- 2-site classification ----

Current state of progress:
- ansatz via extended Pauli basis $(0, +, -, z)$
- set of equations generated via the Reshetikhin condition (proved by A. Hokkyo), $q_{3,j}=[h_j, h_{j+1}]$ from the $R$-matrix expansion
- boundary conditions: the support of the $[H, Q_3]$ commutator is 4, after several trials, it was found that an $L=5$ spin chain properly generates every constraint accounting for the boundary conditions between overlapping sites. Taking the matrix elements of this commutator, the system of equations is obatined

Skeleton of the algorithm:
1. Using the gauge freedom, eliminate the $00, ++, --$ components (there will be two solution families, in the latter, only 10 free symbols will be available initially).
2. Gaussian elimination via monomial linearization.
3. If $\alpha m_i + \beta m_j=0$ obtained somewhere, perform the substitution $m_j=-\frac{\alpha}{\beta}m_i$ for all $m_j$ ($m_{i,j}$ are monomials, the rhs. is 0 as the system of equations is, by construction, homogeneous). We look only for these types of substitutions as this introduces no radicals, exponential blowups in term count, or the untraceable branching of the equations.
4. Look for factorizable equations, if one is found, branch the solution (this avoids any sort of solution pruning).
5. Repeat steps 1-3. until no substitution, new branch, or linear dependence can be found.
6. Compute grevlex Gröbner basis (fastest), if needed, assign this for a better optimized algebra system.
7. Find dimension of the ideal and teh subset of independent variables, define a new base ring with the free variables as constants $\rightarrow$ the ideal with this ring will be zero dimensional, FGLM can obtain the lex basis, from which the solution can be easily (probably?) found.
8. Need to find a generating set for the set of solutions with respect to a change of basis.
