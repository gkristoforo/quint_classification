# The Klobas–Prosen (KP) family — search plan

*NCq3 framework · supervisor brief · 2 November 2026*

This plan follows the per-milestone loop of the development plan. It fixes the family, the facts already computed, the controls and the staged runs with acceptance tests. The design of the builder, of the gauge slice and of the new tests stays with the student; every number marked *(computed)* should be reproduced before it is relied on.

---

## 1. The family, derived from its two charges

### 1.1 Origin

Klobas and Prosen (J. Phys. A 55, 094003 (2022)) introduced two reversible three-state IRF automata, the species-preserving rule $\chi_p$ and the species-flipping rule $\chi_f$. State 0 is a vacancy, states 1 and 2 are two particle species. The rules differ only on the outer pair $00$. Read from their rule table, the middle-site map of $\chi_p$ is a transposition in every block except $(0,0)$:

| block $(a,b)$ | 01 | 10 | 11 | 02 | 20 | 22 | 12 | 21 | 00 |
|---|---|---|---|---|---|---|---|---|---|
| $\chi_p$ swaps | 0↔1 | 0↔1 | 0↔1 | 0↔2 | 0↔2 | 0↔2 | 1↔2 | 1↔2 | none ($\chi_f$: 1↔2) |

Among their elementary conserved quantities are the bond indicators of the ordered pairs $(1,2)$ and $(2,1)$. As one-site translation-invariant bond charges:

$$
N_{12}=\sum_j w_1(s_j,s_{j+1}),\qquad N_{21}=\sum_j w_2(s_j,s_{j+1}),\qquad
w_1(x,y)=\delta_{x1}\delta_{y2},\quad w_2(x,y)=\delta_{x2}\delta_{y1}.
$$

### 1.2 Selection rule

For a bond charge $Q_w=\sum_j w(s_j,s_{j+1})$ and an NCq3 density $h=\sum_{ab}P^a\otimes h^{ab}\otimes P^b$,

$$
[H,Q_w]=0\iff \Delta_{ab}(i,k)[w]\;h^{ab}_{ik}=0\quad\forall\,a,b,\ i<k,\qquad
\Delta_{ab}(i,k)[w]=w(a,i)+w(i,b)-w(a,k)-w(k,b).
$$

*Critical step.* The matrix element of $[Q_w,h_j]$ between $|a\,i\,b\rangle$ and $|a\,k\,b\rangle$ (sites $j,j+1,j+2$) is $\Delta_{ab}(i,k)\,h^{ab}_{ik}$, because only the two bonds touching site $j+1$ differ. A matrix element of $\sum_j[Q_w,h_j]$ between two configurations that differ only at site $m$ receives a contribution only from the term centred at $m$. Hence the sum vanishes iff every coefficient vanishes.

For the two charges:

$$
\Delta^{(1)}_{ab}(i,k)=\delta_{a1}\delta_{i2}+\delta_{i1}\delta_{b2}-\delta_{a1}\delta_{k2}-\delta_{k1}\delta_{b2},\qquad
\Delta^{(2)}_{ab}(i,k)=\delta_{a2}\delta_{i1}+\delta_{i2}\delta_{b1}-\delta_{a2}\delta_{k1}-\delta_{k2}\delta_{b1}.
$$

Edge by edge:

- **Edge $(0,1)$.** Terms with $i=0$ vanish. With $k=1$: $\delta_{a1}\delta_{k2}=0$, $\delta_{k1}\delta_{b2}=\delta_{b2}$, $\delta_{a2}\delta_{k1}=\delta_{a2}$, $\delta_{k2}\delta_{b1}=0$. So $\Delta^{(1)}=-\delta_{b2}$, $\Delta^{(2)}=-\delta_{a2}$. Allowed iff $a,b\neq 2$: blocks 00, 01, 10, 11.
- **Edge $(0,2)$.** Similarly $\Delta^{(1)}=-\delta_{a1}$, $\Delta^{(2)}=-\delta_{b1}$. Allowed iff $a,b\neq1$: blocks 00, 02, 20, 22.
- **Edge $(1,2)$.** $\Delta^{(1)}=\delta_{b2}-\delta_{a1}$, $\Delta^{(2)}=\delta_{a2}-\delta_{b1}$. Checking all nine pairs, both vanish only for $(0,0)$, $(1,2)$, $(2,1)$.

Total: $4+4+3=11$ allowed edges. All $\chi_p$ and $\chi_f$ transpositions are among them; the only addition is the full $(0,0)$ block. The rank of the 11 functionals is 4 *(computed)*, so the conserved bond charges are exactly the 3 trivial ones ($\alpha(a)-\alpha(i)+c$) plus $\mathrm{span}(w_1,w_2)$.

### 1.3 Definition

> **KP family:** the most general NCq3 density conserving $N_{12}$ and $N_{21}$,
> $$h^{ab}=D^{ab}+\sum_{(i,k)\in\mathcal E_{ab}}\big(z^{ab}_{ik}\,|i\rangle\langle k|+\bar z^{ab}_{ik}\,|k\rangle\langle i|\big),$$
> with $D^{ab}$ a full real diagonal and allowed edge sets
>
> | $\mathcal E_{ab}$ | blocks |
> |---|---|
> | $\{01,02,12\}$ | 00 |
> | $\{01\}$ | 01, 10, 11 |
> | $\{02\}$ | 02, 20, 22 |
> | $\{12\}$ | 12, 21 |

### 1.4 Facts *(computed)*

| quantity | value |
|---|---|
| unknowns, real version ($z\in\mathbb R$) | $27+11=38$ |
| unknowns, complex version | $27+22=49$ |
| phase-gauge rank $r=9-\dim W$ | 4 |
| gauge-invariant off-diagonal data (complex) | $22-4=18$ = 11 moduli + 7 fluxes |
| stabiliser in $S_3\times$ reflection | $\{1,\rho,\pi_{12},\rho\pi_{12}\}$, with $\pi_{12}$ = `relabel([0,2,1])` |
| activity-pattern orbits | 167 (orbit sizes sum to 511) |
| single-U(1) families containing KP | 3: the winding family (15 edges) and two symmetry-related 17-edge families conserving only $N_{12}$ or only $N_{21}$ |

So KP is the intersection of the winding family with the conservation of $N_{12}$ and $N_{21}$.

---

## 2. Exact reductions

### 2.1 Diagonal slice (always)

Five diagonal directions change neither $H$ nor $Q_5$:

$$
\delta h^{ab}_{ii}=c,\qquad
\delta h^{ab}_{ii}=\alpha(a)+\alpha(b)-2\alpha(i),\qquad
\delta h^{ab}_{ii}=\beta_1\big[w_1(a,i)-w_1(i,b)\big]+\beta_2\big[w_2(a,i)-w_2(i,b)\big].
$$

*Why they are exact.* The first shifts $H$ by a constant. The other two are telescoping shifts $h\mapsto h+g_{j,j+1}-g_{j+1,j+2}$, under which

$$
Q_5(h+\delta g)=Q_5(h)+[H,G],\qquad G=\sum_j g_{j,j+1}.
$$

For $g(a,c)=\alpha(a)-\alpha(c)$ the operator $G$ is identically zero (telescoping). For $G=\beta_1N_{12}+\beta_2N_{21}$, every member of the family commutes with $G$ by construction.

*Proof of the identity.* With $h'_j=h_j+g_j-g_{j+1}$:

$$
Q_5(h')=Q_5(h)+\sum_j\Big([h_j,g_{j+1}]-[h_j,g_{j+2}]+[g_j,h_{j+1}]-[g_{j+1},h_{j+1}]\Big)+\sum_j[g_j-g_{j+1},g_{j+1}-g_{j+2}].
$$

The last sum vanishes (all diagonal). $[h_j,g_{j+2}]=0$ and $[g_j,h_{j+1}]=0$, because in each case the only non-diagonal site of the density is not touched by $g$. Shifting the index, $-\sum_j[g_{j+1},h_{j+1}]=\sum_j[h_j,g_j]$. So $Q_5(h')=Q_5(h)+\sum_j([h_j,g_j]+[h_j,g_{j+1}])$, which equals $[H,G]$ because $[h_j,g_k]\neq0$ only for $k\in\{j,j+1\}$.

*Implementation rule.* The five generators are independent (rank 5, *computed*). Impose five homogeneous linear functionals whose $5\times5$ pairing with the generators is invertible. The diagonal drops from 27 to 22 unknowns, and every component loses a spurious 5-dimensional cylinder direction.

*Notes.*

- No set of five zeroed single entries is invariant under the stabiliser *(checked)*. An invariant slice is not needed: orbit reduction only claims $V_S\cong V_{gS}$ for the unfixed family, so run `check_symmetries` **before** slicing.
- If a symmetric slice is wanted anyway, one exists using character-weighted orbit sums; the isotypic multiplicities of the gauge directions are 2, 2, 1.
- Never fix the other six telescoping directions under `'none'`: they change $Q_5$ by a non-zero $[H,G]$.

### 2.2 Phase fixing (complex version, per stratum)

The controlled-phase gauge $U=e^{i\sum_j\varphi(s_j,s_{j+1})}$ acts as $h^{ab}_{ik}\mapsto e^{i\Delta_{ab}(i,k)[\varphi]}h^{ab}_{ik}$ and maps $[H,Q_5]$ to $U[H,Q_5]U^\dagger$, so it is an exact symmetry. In stratum $S$, choose $r_S\le4$ active edges with independent $\Delta_e$ and impose $\operatorname{Im}h_e=0$. No solution is lost: non-zero entries can be rotated real, zero entries satisfy the constraint already.

Place these equations in the chart's `extra` list (stratum-dependent; Singular eliminates linear equations cheaply). In the real version only sign gauges remain; expect sign-gauge copies of solutions and deduplicate them in triage.

### 2.3 Chart normaliser

Proposition 4.17 allows any active block to carry the normalisation $F_c=1$. Use an active block with the fewest forms, so that $(0,0)$ (3 forms) is the normaliser only when it is the sole active block. This gives **169** real charts instead of **335** with the sorted-first convention *(computed)*.

---

## 3. Controls *(computed)*

Values are $\max|[H,Q_5]_{ij}|$ on periodic chains (`verify_periodic`).

| point | expected | $L=7$ | $L=8$ |
|---|---|---|---|
| bond–site image of su(3) XX: hopping 1 on the six off-diagonal blocks, diagonal 0 | passes | 0 | 0 |
| same, unequal species amplitudes $t_2=0.37\,t_1$ | fails | 0.863 | 0.863 |
| $\chi_p$ generator (unit transposition generators) | — | 2.0 | 2.0 |
| $\chi_f$ generator | — | 2.0 | 2.0 |
| bXYZ on $\{0,1\}$ with state 2 inert, $(J_x,J_y,J_z)=(0.7,-0.4,1.1)$ | — | 6.556 | 6.556 |
| random complex KP member | fails | 55.4 | 45.0 |

For the random member, $[H,N_{12}]$ and $[H,N_{21}]$ vanish exactly at both lengths, confirming the selection rule.

**Positive control.** The $\tau$-covariant part of KP ($h^{a+1,b+1}=\tau h^{ab}\tau^\dagger$) maps under $d_j=s_{j+1}-s_j \bmod 3$ to two species hopping past vacancies without exchanging. It contains the image of the su(N) XX model of Maassarani and Mathieu (Nucl. Phys. B 517 (1998) 395). The image is genuine for the exact tests: $Q_5\neq0$ (not vacuous), and its density has a non-zero range-3 canonical part (not two-local). It has 4 range-2 local charges, so the reference's numeric verdict ladder would call it "suspect"; that is expected for an XX-type model whose spectrum is free-like. It tests the solver, not novelty. It lies in the stratum $S_\tau$ of the six off-diagonal blocks (real $\kappa=6$) and is rational, so the acceptance test is exact.

**Dictionary fact.** A qubit solution with an inert third state is not automatically an NCq3 solution (bXYZ row). Frozen-state embeddings must come out of the solver.

---

## 4. Cost

The predictor is the active kinetic count per stratum after phase fixing,

$$
\kappa_{\rm eff}(S)=\sum_{k\in S}m_k-r_S,
$$

with $m_k$ the number of real off-diagonal forms of block $k$ (no $r_S$ term in the real version). Working hypothesis from earlier anchors (to be calibrated): $\kappa_{\rm eff}\le9$ → minutes per chart; 10–12 → calibration zone; $\ge13$ → expect time-outs.

| version | $\kappa\le9$ | $\kappa=10$–$11$ | $\kappa\ge12$ | max |
|---|---|---|---|---|
| real | 163 patterns | 4 | 0 | 11 |
| complex | 109 | 26 | 32 | 18 |

Refining block $(0,0)$ by edge activity helps only partly in the complex version (29 sub-strata still at $\kappa\ge14$), because the hardest strata are dominated by single-edge blocks. Refine adaptively, on time-out only.

Budget rule with $W$ workers:

$$
T\approx\frac{N_{\rm chart}\,(2t_0+2\bar t_p)+\sum_{\text{non-empty}}t_{\rm solve}}{W}\le 8\ \text{h}.
$$

**Pre-flight before every run:** probe all charts, solve three non-empty ones, extrapolate, then commit.

---

## 5. Stages

| stage | what runs | size | budget | acceptance / gate |
|---|---|---|---|---|
| KP-0 | builder, `check_symmetries` on the unfixed family, U(1) check, exact substitution of controls | no Singular | 1 evening | XX point satisfies the generated equations exactly; $\chi_p$ point does not |
| KP-1 | real family, all 167 patterns | 169 charts, 22 diagonal unknowns | ≤ 2 h | XX lies on a component classed *genuine* in $S_\tau$; $\chi_p$ point lies on none |
| KP-2 | triage of KP-1 | — | 1 evening | every genuine component labelled |
| KP-3a | complex family, $\kappa_{\rm eff}\le11$ | 135 patterns | 1 night | pre-flight projects ≤ 8 h |
| KP-3b | complex family, $\kappa_{\rm eff}\ge12$ | 32 patterns | 1 night | refine on time-out; unfinished strata listed explicitly, never reported as empty |
| KP-4 | winding family, then the $N_{12}$-only family | separate plan | later | only after KP-3a is triaged |

**KP-1 acceptance in detail.** Decide whether the XX point lies on a component by evaluating that component's generators at the point; do not rely on the $\tau$-test. If the component containing XX is larger and not $\tau$-covariant, that is the most interesting outcome of the stage: a non-covariant deformation of a bond–site image. If XX is not found at all, there is a bug.

**Triage (KP-2), in order:**

1. Existing exact classes: excluded, vacuous, two-local.
2. New exact test — $\tau$-covariance: the real and imaginary parts of $h^{a+1,b+1}-\tau h^{ab}\tau^\dagger$. Pass ⇒ bond–site image; compare with two-species nearest-neighbour literature (XX, XXC-type deformations).
3. New exact test — frozen state $c$: all off-diagonal entries in row/column $c$ of every block.
4. Real points.
5. Periodic check at $L=7$ **and** $L=8$.
6. Range-2 local-charge count.

Deduplicating sign-gauge, relabelling and reflection copies requires the student's own criterion for "same model up to gauge" (ground rule 3).

**Trap checklist.**

- Homogeneity: every slice and phase equation is `= 0`; only the chart normalisation is `= 1` (C1).
- Symmetries verified on the unfixed family (C2).
- All 167 patterns solved, no `max_active` cut, so excluded components are counted at smaller patterns (C4).
- Re-probe empty charts with two further primes before final claims (C5).

---

## 6. Expectations

The real pass is cheap (comparable to family B, at most an order of magnitude more) and, unlike B, contains a known non-trivial solution, so it cannot come back silently empty. Most components will probably be bond–site images (XX type and diagonal deformations), vacuous models or fragmented models.

A genuinely new KP model is a genuine component that fails the $\tau$-test, has no frozen state, has real points, passes at $L=7$ and $L=8$, and has few local charges.

On the link to the automata: Prosen (arXiv:2608.13080) leaves a Yang–Baxter or Lax structure for these rules open. A classification of the KP family tests whether any integrable continuous-time Hamiltonian shares the two bond charges. It is a natural first test of that question, not an answer to it, because a Lax operator need not conserve these charges away from its special point.

---

## 7. Next step

Send the design in one message:

- the edge-mask builder and its output contract;
- where the diagonal slice and the per-stratum phase equations enter;
- the normaliser rule;
- the $\tau$-covariance and frozen-state tests.

**Hard case to state with the design:** the complex stratum with all nine blocks active, $\kappa_{\rm eff}=18$. Stress tests follow, then certification.

---

## References

- T. Gombor, B. Pozsgay, *Integrable spin chains and cellular automata with medium range interaction*, arXiv:2108.02053.
- K. Klobas, T. Prosen, J. Phys. A 55, 094003 (2022).
- T. Prosen, *Local and quasilocal conservation laws of three-state IRF cellular automata and their quantum deformations*, arXiv:2608.13080 (2026).
- Z. Maassarani, P. Mathieu, *The su(N) XX model*, Nucl. Phys. B 517 (1998) 395, arXiv:cond-mat/9709163.
- irfsolve implementation reference (16 Sept 2026) and pipeline deep dive (10 Oct 2026), project files.


# KP family — supervision session notes

*NCq3solve project · real Klobas–Prosen (KP) ansatz · 4 November 2026*

This file records the conversation that followed the KP family search plan (`KP_family_plan.md`, 2 Nov 2026). It keeps the derivations, code reviews, gates and open items in the order they were discussed. Numbers marked *(computed)* were checked in the container and should still be reproduced before they are relied on.

---

## Contents

1. Where to begin: derivation roadmap
2. The commutator $[Q_w,h_j]$ and the selection rule
3. Exact reductions for the real ansatz
4. Review of `check_symmetries`
5. Review of the KP builder
6. What to do now: Steps 1–7 to KP-1
7. The goal of Step 3 (diagonal slice)
8. The five diagonal translation symmetries, spelled out
9. Open items

---

## 1. Where to begin: derivation roadmap

**Assessment.** The Klobas–Prosen paper is not needed to understand the computation. The KP family is defined entirely by the two bond charges $N_{12}$ and $N_{21}$. The automata $\chi_p$, $\chi_f$ serve only as motivation and as diagnostic points. Start from the charges, not from the automata.

Steps 1–3 suffice for KP-0 (the builder), Steps 4–5 for the slice and phase equations, Step 6 for the acceptance test.

| step | content | critical point | check |
|---|---|---|---|
| 0 | rule table of $\chi_p$ | why the charges are bond charges | — |
| 1 | selection rule | isolating a single $h_j$ in the sum | 11 edges, rank 4, 5 conserved bond charges |
| 2 | family and unknowns | all $\chi_p,\chi_f$ transpositions lie in the masks | 38 real / 49 complex |
| 3 | stabiliser, orbits, charts | block map of $\rho\pi_{12}$ is $(a,b)\mapsto(\pi b,\pi a)$ | 167 orbits; 169 vs 335 charts |
| 4 | phase gauge (complex) | no solution lost when fixing $\operatorname{Im}h_e=0$ | rank 4; $18=11+7$ |
| 5 | diagonal slice | $[h_j,h_{j+2}]=0$ is load-bearing | 5 generators |
| 6 | bond–site map, XX control | $\tau$-covariance $\Leftrightarrow$ $H$ depends only on $d_j$ | control table values |

Burnside checks for Step 3:

$$
\frac{512+64+32+64}{4}-1=167,\qquad
\text{orbits with }(0,0)\text{ active}=\frac{256+32+16+32}{4}=84,
$$

$$
\text{sorted-first: }84\cdot3+83\cdot1=335,\qquad \text{fewest-forms: }166\cdot1+1\cdot3=169 .
$$

**Discrepancy in the plan.** Sec. 2.1 says "never fix the other **six** telescoping directions under `'none'`". The diagonal telescoping map $g\in\mathbb R^9\mapsto\delta D$ has a 1-dimensional kernel (constant $g$), so its image is 8-dimensional. Of these, 4 are harmless ($\alpha$: 2, $\beta$: 2), leaving **four** forbidden directions *(computed)*. This does not change any computation; the slice fixes 5 directions either way.

Not to be derived: the cost model (Sec. 4, empirical) and the symmetric-slice remark (isotypic multiplicities 2, 2, 1).

---

## 2. The commutator $[Q_w,h_j]$ and the selection rule

### 2.0 Conventions

$P^c=|c\rangle\langle c|$, $s_j\in\{0,1,2\}$, periodic chain with $L\ge4$. The projector-product image of a bond function $w$ is

$$
Q_w=\sum_jq_j,\qquad q_j=\sum_{x,y}w(x,y)\,P^x_jP^y_{j+1},
$$

so for $w_1(x,y)=\delta_{x1}\delta_{y2}$: $N_{12}=\sum_jP^1_jP^2_{j+1}$. The density is

$$
h_j=\sum_{a,b}P^a_j\,h^{ab}_{j+1}\,P^b_{j+2},\qquad h^{ab}_{ik}=\langle i|h^{ab}|k\rangle,\qquad H=\sum_jh_j .
$$

### 2.1 Which terms fail to commute

Every $q_m$ is diagonal. The only non-diagonal site of $h_j$ is $j+1$. So $[q_m,h_j]\neq0$ requires $q_m$ to act on $j+1$, i.e. $m\in\{j,j+1\}$:

$$
[Q_w,h_j]=[\,q_j+q_{j+1},\,h_j\,].
$$

### 2.2 The matrix element

On $|a\,i\,b\rangle$ (sites $j,j+1,j+2$): $q_j|aib\rangle=w(a,i)|aib\rangle$, $q_{j+1}|aib\rangle=w(i,b)|aib\rangle$. Define $E_{ab}(i)=w(a,i)+w(i,b)$.

$$
h_j|akb\rangle=\sum_{a'b'}\delta_{a'a}\delta_{b'b}\,|a\rangle\otimes h^{a'b'}|k\rangle\otimes|b\rangle=\sum_ih^{ab}_{ik}|aib\rangle .
$$

$$
\langle aib|(q_j+q_{j+1})h_j|akb\rangle=E_{ab}(i)h^{ab}_{ik},\qquad
\langle aib|h_j(q_j+q_{j+1})|akb\rangle=h^{ab}_{ik}E_{ab}(k).
$$

$$
\boxed{\;\langle aib|[Q_w,h_j]|akb\rangle=\Delta_{ab}(i,k)[w]\,h^{ab}_{ik},\qquad
\Delta_{ab}(i,k)[w]=w(a,i)+w(i,b)-w(a,k)-w(k,b)\;}
$$

Operator form (Hadamard product):

$$
[Q_w,h_j]=\sum_{a,b}P^a_j\Big(\sum_{i,k}\Delta_{ab}(i,k)\,h^{ab}_{ik}\,|i\rangle\langle k|_{j+1}\Big)P^b_{j+2}.
$$

Properties: $\Delta_{ab}(i,i)=0$ (the diagonal never enters), and $\Delta_{ab}(k,i)=-\Delta_{ab}(i,k)$, which pairs with $h^{ab}_{ki}=\overline{h^{ab}_{ik}}$, so $i<k$ suffices.

### 2.3 Specialisation to $N_{12}$

$$
\Delta^{(1)}_{ab}(i,k)=\delta_{a1}(\delta_{i2}-\delta_{k2})+\delta_{b2}(\delta_{i1}-\delta_{k1}).
$$

Operator cross-check: $[N_{12},h_j]=P^1_j[P^2_{j+1},h_j]+P^2_{j+2}[P^1_{j+1},h_j]$, using $\langle i|[P^c,h]|k\rangle=(\delta_{ci}-\delta_{ck})h_{ik}$.

- Block $(0,2)$, edge $(0,1)$: $\Delta^{(1)}=0+1\cdot(0-1)=-1$. Forbidden: $|0\,0\,2\rangle\to|0\,1\,2\rangle$ creates a $(1,2)$ bond.
- Block $(1,2)$, edge $(1,2)$: $\Delta^{(1)}=-1+1=0$ and $\Delta^{(2)}=0$. Allowed: $|1\,1\,2\rangle\leftrightarrow|1\,2\,2\rangle$ moves the $(1,2)$ bond.

### 2.4 From one term to the sum (critical step)

Let $C_j=[Q_w,h_j]$, so $[H,Q_w]=-\sum_jC_j$.

1. $C_j$ has zero diagonal, so it changes site $j+1$ only.
2. For configurations $s,s'$ differing only at site $m$ ($s_m=i\neq k=s'_m$), only $C_{m-1}$ contributes:
$$
\Big\langle s\Big|\sum_jC_j\Big|s'\Big\rangle=\Delta_{ab}(i,k)\,h^{ab}_{ik},\qquad a=s_{m-1},\ b=s_{m+1}.
$$
3. For $L\ge3$ every $(a,b,i,k)$ is realised, hence
$$
[H,Q_w]=0\iff\Delta_{ab}(i,k)[w]\,h^{ab}_{ik}=0\quad\forall\,a,b,\ i<k .
$$

Both steps rely on the zero diagonal of $C_j$.

*Exercise left open:* the edge table for $N_{12}$ alone (17 edges expected from the plan; to be checked against the hand count).

---

## 3. Exact reductions for the real ansatz

### 3.0 Setting

$x=(D,z)\in\mathbb R^{27}\times\mathbb R^{11}$; under `'none'` each equation $f_r$ is a homogeneous cubic; $V=\{f=0\}$.

A reduction is **exact** if it rests on $T$ with $f(Tx)=M_Tf(x)$ for all $x$, $M_T$ invertible. Two checks every time:

- **(E1)** the invariance is an identity on the whole family, not only on $V$;
- **(E2)** the representative rule meets every orbit and keeps the equations homogeneous (Prop. 4.17).

### 3.1 Scale: projective charts

$f(\lambda x)=\lambda^3f(x)$, so $V$ is a cone. Strata are an exact partition (inactive forms $=0$, active blocks saturated by $J_k$); the diagonal plays no role in $S$.

Chart $U_r$ of normaliser block $c$: $F_1=\dots=F_{r-1}=0$, $F_r=1$. Each punctured line meets $\bigcup_rU_r$ exactly once.

Charts per pattern $=m_c$. In real KP, $m_{00}=3$, all other $m_k=1$. Fewest-forms rule: 169 charts; sorted-first: 335. In `stratum`, saturate every active block except the normaliser and drop blocks whose form list equals the normaliser's.

### 3.2 Diagonal translations: the 5-dimensional slice

| generator | $\delta D^{ab}_{ii}$ | count | reason |
|---|---|---|---|
| $v_c$ | $1$ | 1 | $H\to H+tL\mathbb 1$ |
| $v_\alpha$ | $\alpha(a)+\alpha(b)-2\alpha(i)$ | 2 | $G\equiv0$ |
| $v_\beta$ | $\beta_1[w_1(a,i)-w_1(i,b)]+\beta_2[w_2(a,i)-w_2(i,b)]$ | 2 | $[H,G]\equiv0$ by edge masks |

Rank 5 *(computed)*. Telescoping image is 8-dimensional; 4 harmless, 4 forbidden *(computed)*.

**Slice.** Impose $\ell_r(D)=0$, $r=1..5$, with $M_{rs}=\ell_r(v_s)$. If $\det M\neq0$: $\ell(D+\sum t_sv_s)=\ell(D)+Mt=0\Rightarrow t=-M^{-1}\ell(D)$, unique. So $V\cong(V\cap\text{slice})\times\mathbb R^5$.

- Purely diagonal $\ell_r$ commute with strata, charts and saturations.
- Homogeneous ($=0$), so Prop. 4.17 survives (C1).
- Assert $\det M\not\equiv0\pmod p$ for probe primes.
- Gain: $27\to22$ diagonal unknowns; every component loses 5 spurious dimensions.

### 3.3 Discrete symmetries: 511 → 167 patterns

Stabiliser $\{1,\rho,\pi_{12},\rho\pi_{12}\}$.

- Relabelling: conjugation by $\Pi^{\otimes L}$; blocks $(a,b)\to(\pi a,\pi b)$, middle $\Pi h^{ab}\Pi^T$.
- Reflection: $\sum_j[h_j,h_{j+1}]\to\sum_k[\tilde h_{k+1},\tilde h_k]=-\tilde Q_5$; the sign is irrelevant.

Conditions: each element maps the real family into itself (`check_symmetries` on the **unfixed** family) and maps strata to strata. Slice and normaliser need not be invariant: the orbit claim is about unfixed $V_S$. Patterns with non-trivial stabiliser (e.g. $\{00\}$) give $g$-image copies inside one stratum — triage work.

### 3.4 Sign gauge: a duplication, not a reduction

$z_e\to(-1)^{\Delta_e[n]}z_e$, $n\in\{0,1\}^9$. $\mathbb F_2$-rank 4 *(computed)*, so up to 16 sign copies. Not quotientable algebraically ($z_e>0$ is an inequality); deduplicate in triage. A flip of $z_c$ combined with $\lambda=-1$ maps a chart point to $-UHU^\dagger$ in the same chart.

The real ansatz is a **restriction**: it sets the 7 gauge-invariant fluxes to $0$ or $\pi$.

### 3.5 Order of operations

| step | operation | requires |
|---|---|---|
| 1 | build the 38-unknown family | — |
| 2 | `check_symmetries` on the unfixed family | before any slice |
| 3 | 511 patterns → 167 representatives | step 2 passed |
| 4 | per pattern: inactive forms $=0$, normaliser, charts, saturations | homogeneous `extra` |
| 5 | diagonal slice | $\det M\neq0$, also mod $p$ |

Free unknowns per chart: $22+\sum_{k\in S}m_k-r$; maximum $32$.

---

## 4. Review of `check_symmetries`

The logic matches the reference (span the image with $n+2$ random members, least-squares membership with relative tolerance).

| severity | issue | fix |
|---|---|---|
| must fix | `g['mid'] == None` raises `ValueError` for every relabelling (elementwise comparison of an array, checked with NumPy 2.4) | `g['mid'] is None` |
| should fix | `lstsq` default `rcond` is version-dependent; basis has two zero singular values | `rcond=None` |
| should fix | unseeded global RNG, failures not reproducible | `rng = np.random.default_rng(seed)` |
| should fix | values drawn over `nametable`, basis sized by `nvars` | `assert len(self.nametable) == n` |
| optional | non-bijective block map silently overwrites keys | `assert len(new) == d*d` |

Notes: checking the generators $\rho,\pi_{12}$ suffices; "into" implies "onto" for injective linear maps; run on the unfixed family only; the method tests family invariance, not invariance of $[H,Q_5]=0$ (which is analytic).

**Acceptance test (with negative controls)** on the unfixed real KP family:

| symmetry | expected | reason |
|---|---|---|
| `reflection(3)` | `True` | swaps $N_{12}\leftrightarrow N_{21}$ |
| `relabel([0,2,1])` | `True` | swaps $N_{12}\leftrightarrow N_{21}$ |
| `relabel([1,0,2])` | `False` | maps $N_{12}\to N_{02}$ |
| `relabel([1,2,0])` | `False` | maps $N_{12}\to N_{20}$ |

---

## 5. Review of the KP builder

Code under review: `KP_DICT` in `operators.py`, `KP_block(a, b, d)` in `ansatz.py`, `'KP'` branches in the hdict builder and `families.py`.

**Verified** *(computed)*: 49 real symbols ($27+2\cdot11$); all blocks Hermitian; `KP_DICT` equals the masks computed from $\Delta^{(1)}=\Delta^{(2)}=0$.

| severity | issue |
|---|---|
| must fix | no real version: every edge carries $u$ and $v$; KP-1 needs 38 unknowns → add a `real` flag and a `'KP_R'` family |
| must fix | check that the `KP` branch honours `pattern` (or raises when it is passed); `reflection` works only because the masks are $(a,b)\leftrightarrow(b,a)$ symmetric |
| must fix | `d` parameter inconsistent with hard-coded `range(3)` → `assert d == 3` |
| design (deferred) | derive masks from the charges (`edge_mask(charges)`), keep `KP_DICT` as a unit test; needed for KP-4 |
| design (deferred) | reuse `hermitian_block` naming (`h{a}{b}_x{i}`, `_u{i}{j}`, `_v{i}{j}`) |
| design (deferred) | move `KP_DICT` out of `operators.py` |
| must add | `symmetries_for('KP')` returning `[reflection(3), relabel([0,2,1])]` |

Builder acceptance tests:

| test | expected |
|---|---|
| `nvars` of `'KP_R'` / `'KP'` | 38 / 49 |
| `edge_mask([w1, w2]) == KP_DICT` | `True` |
| `check_symmetries` with the four symmetries above | `True, True, False, False` |
| $\max|[H,N_{12}]|$, $\max|[H,N_{21}]|$ at $L=7,8$ | machine zero |

**Status reported by the student:** NCq3solve has `_orbit`; symmetries are tested numerically at `EquationSystem` initialisation; the edge-mask generalisation is deferred in favour of results.

---

## 6. What to do now: Steps 1–7 to KP-1

### Step 1 — Minimal builder fixes (≈30 min)

`real` flag and `'KP_R'`; `assert d == 3`; `symmetries_for` for both names. Optional insurance: a ten-line test recomputing the masks from the selection rule.

**Gate:** `nvars` 38/49; symmetry table `True, True, False, False`; 167 orbits from `_orbit`.

### Step 2 — Symmetry check vs. slice

Because the check runs at `EquationSystem` initialisation, it checks whatever family the system is built from.

- Slice via `restrict` → the check sees the sliced family and fails for $\pi_{12}$.
- Slice via `extra` → the check sees the unfixed family and passes. **Use this.**

A failed check must **raise**, not warn; otherwise `_orbit` merges inequivalent patterns silently.

### Step 3 — Diagonal slice (design first)

Deliver: the five generators as explicit 27-vectors (with the chosen basis of $\alpha$ mod constants), five purely diagonal homogeneous functionals $\ell_r$, the matrix $M_{rs}=\ell_r(v_s)$ entry by entry, $\det M$.

**Gate:** $\det M\neq0$ over $\mathbb Q$ and mod every probe prime; runtime assert recomputing $M$ from the code's own vectors.

Hint: single-entry functionals make $M$ the generators read off at five entries.

### Step 4 — Normaliser rule in `stratum`

Fewest forms, ties broken by sorted order; saturate all active blocks except the normaliser.

| setting | expected charts |
|---|---|
| 167 representatives, fewest-forms | 169 |
| 511 patterns, fewest-forms | $510+3=513$ |
| 511 patterns, sorted-first | $256\cdot3+255=1023$ |

### Step 5 — Exact control substitution (KP-0 gate)

Evaluate the **raw** FLINT equations of `'KP_R'` under `'none'` (no slice, no chart) at rational points.

| point | definition | expected |
|---|---|---|
| su(3) XX image | $z^{ab}=1$ on the single edge of each block with $a\neq b$; $z^{00}=0$; all $D=0$ | all $=0$ |
| unequal amplitudes | amplitude $37/100$ on the three blocks with edge $(0,2)$ | some $\neq0$ |
| $\chi_p$ generator | transposition permutation matrices, block $00=\mathbb 1$ | some $\neq0$ |

Reproduce the periodic control values ($0$, $0.863$, $2.0$ at $L=7,8$) and the U(1) check. Verify by hand first that every entry set to 1 lies on an allowed edge.

### Step 6 — Pre-flight

Probe all 169 charts; solve three non-empty ones including $\kappa=11$; extrapolate with

$$
T\approx\frac{N_{\rm chart}(2t_0+2\bar t_p)+\sum t_{\rm solve}}{W};
$$

commit only if $T\le2$ h. Trap list: C1 (homogeneity), C4 (no `max_active`), C5 (re-probe empty charts with two more primes).

### Step 7 — KP-1 run and acceptance

All 167 representatives. XX must lie on a *genuine* component of $S_\tau$, decided by evaluating the component's generators at the point after mapping it into the chart and slice ($t=-M^{-1}\ell(D)$, divide by the normaliser entry). $\chi_p$ lies on no component. If XX is not found, there is a bug — no triage before that.

**Order:** Steps 1, 2, 4 (mechanical) → Step 3 → Step 5 → Steps 6, 7.

---

## 7. The goal of Step 3 (diagonal slice)

Remove the five diagonal directions along which the KP equations are exactly invariant, so that Singular solves a smaller system and every returned component is free of five spurious dimensions, without losing any solution.

Constraints: homogeneous; purely diagonal; $\det M\neq0$ (also mod $p$); passed as `extra`. Gain: $27\to22$ diagonal unknowns per chart.

---

## 8. The five diagonal translation symmetries, spelled out

### 8.1 Statement

$T_t(D,z)=(D+\sum_{s=1}^5t_sv_s,\,z)$. Claim, as a polynomial identity:

$$
f(T_tx)=f(x)\qquad\forall\,x,\ \forall\,t\in\mathbb R^5 .
$$

Since canonical form is unique, it suffices to show $[H(T_tx),Q_5(T_tx)]=[H(x),Q_5(x)]$ as operators.

### 8.2 Lemma A — $Q_5=\sum_j[h_j,h_{j+1}]$ is exact

$h_j$ is non-diagonal only at $j+1$, $h_{j+2}$ only at $j+3$; they overlap only at $j+2$, where both are diagonal. So $[h_j,h_{j+2}]=0$. Diagonal shifts keep the NCq3 form, so this holds for every translated density.

### 8.3 Lemma B — telescoping diagonal shifts

If $\delta D^{ab}_{ii}=g(a,i)-g(i,b)$ and $g_j=\sum_{x,y}g(x,y)P^x_jP^y_{j+1}$, then on $|a\,i\,b\rangle$:

$$
\delta h_j|aib\rangle=\delta D^{ab}_{ii}|aib\rangle,\quad g_j|aib\rangle=g(a,i)|aib\rangle,\quad g_{j+1}|aib\rangle=g(i,b)|aib\rangle,
$$

so $\delta h_j=g_j-g_{j+1}$ and $H'=H+\sum_jg_j-\sum_jg_{j+1}=H$.

### 8.4 Lemma C — $Q_5(h')=Q_5(h)+[H,G]$

$$
\begin{aligned}
Q_5(h')=\;&\textstyle\sum_j[h_j,h_{j+1}] &&(1)\\
&+\textstyle\sum_j[h_j,g_{j+1}]-\sum_j[h_j,g_{j+2}] &&(2),(3)\\
&+\textstyle\sum_j[g_j,h_{j+1}]-\sum_j[g_{j+1},h_{j+1}] &&(4),(5)\\
&+\textstyle\sum_j[g_j-g_{j+1},\,g_{j+1}-g_{j+2}] &&(6)
\end{aligned}
$$

- (1) $=Q_5(h)$; (6) $=0$ (diagonal operators commute).
- (3) $=0$: overlap at $j+2$, both diagonal. (4) $=0$: overlap at $j+1$, both diagonal.
- (5) $=\sum_j[h_j,g_j]$ by antisymmetry and an index shift.

So $Q_5(h')=Q_5(h)+\sum_j([h_j,g_j]+[h_j,g_{j+1}])$. In $[H,G]=\sum_{j,k}[h_j,g_k]$ only $k\in\{j,j+1\}$ survive, giving the same sum. Consequence:

$$
[H',Q_5']=[H,Q_5]+[H,[H,G]] .\tag{$\ast$}
$$

### 8.5 The generators

- **$v_c$**: $\delta h_j=t\,\mathbb 1$; $H'=H+tL\mathbb 1$, $Q_5'=Q_5$, so $[H',Q_5']=[H,Q_5]$.
- **$v_\alpha$**: $g(x,y)=\alpha(x)-\alpha(y)$ gives $\delta D^{ab}_{ii}=\alpha(a)+\alpha(b)-2\alpha(i)$; $G=\sum_jA_j-\sum_jA_{j+1}=0$. Constant $\alpha$ gives $\delta D=0$, so $\alpha$ mod constants: 2 directions.
- **$v_\beta$**: $G=\beta_1N_{12}+\beta_2N_{21}$. **Critical step:** every $h^{ab}_{ik}$ with $\Delta\neq0$ is identically zero by the edge masks, so $[H,G]\equiv0$ as a polynomial identity and $(\ast)$ gives invariance. Imposing $[H,N_{12}]=0$ by equations instead would make the $\beta$-directions symmetries only on that zero set.

### 8.6 Combining

The five translations act only on $D$ and commute; each identity holds for every amplitude and every $x$. Hence $f(D+\sum t_sv_s,z)=f(D,z)$ for all $x,t$. The remaining 4 telescoping directions leave $[H,[H,G]]\neq0$ in $(\ast)$ and are not symmetries under `'none'`.

---

## 9. Open items

1. Builder fixes (Step 1) and their gates.
2. Confirm that a failed symmetry check raises.
3. Slice design: generators as 27-vectors, functionals $\ell_r$, explicit $M$ and $\det M$ — **send for review before coding**.
4. Normaliser rule in `stratum` and the chart-count gate.
5. Exact control substitution (XX, unequal amplitudes, $\chi_p$) and the U(1) check.
6. Hand count of the $N_{12}$-only edge table.
7. Deferred: edge-mask generalisation (needed for KP-4).