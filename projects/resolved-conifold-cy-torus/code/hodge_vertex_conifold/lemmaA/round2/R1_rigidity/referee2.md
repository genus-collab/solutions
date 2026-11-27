# Lemma A round 2: referee2 R1 (saved by the main session from the workflow journal, 2026-10-31)

Subagents cannot write report files, so this file holds the agent's structured return value: the fields first, then the full report text (report_text) verbatim.

## confirmed

- T1 (closed form of the polar-order relations). Proof s.2-s.7 re-derived by hand. An independent engine (Theorem D for all R, closed-form genus-0 kernel E0 = sigma Pi_n Pi^Q Xi_0(Q+{1,n})) reproduces T1 with 0 mismatches: to order 32 at tau = 2/7 (p = 10^9+7) and tau = -3/11 (p = 998244353); coefficient-wise in the one-leaf data to order 28 at tau = -1/5 and order 24 at tau = 2/7. Three deliberate perturbations are all detected.
- KP propagation is consistent: three different pair choices in Theorem D give identical multi-leaf values for all multisets at K = 8, 16 and 26.
- Flow theorem (Thm 5.1) and free energy on the slice (Thm 6.1): re-derived by hand; 7 random paths, orders 6-12, framings 2/7, -3/5, 5/3, -2/7, -1/5, 3; 816 coefficients, 0 mismatches.
- Branch-point identity (Lemma 4.1): re-derived; j <= 40 at tau = 2/7 and -5/13.
- Residue identity and Pfaff step (s.7(c)): re-derived; n <= 8 to t^12 at tau = 2/7 and -3/11.
- Curve identity FC (Thm 7.1) with C = W e^(tau^2 t)/H_1: re-derived; holds to t^13 at tau = 2/7, -3/5, 3, -1/5, -2/9. T1 needs only the leaf ratios b_j; the zeta of Theorem 11.1 cancels exactly.
- Leaf curve: b_j = A_j alpha^-j/c^2 = (-1)^j (Pi_j^2/j) t^(j-1) H_j/H_1^j follows from Theorem 11.1 by exact algebra. An independent rebuild of A_j with its own zeta matches q^(1..3) of s.11.4, the sum rule and W = alpha (275 comparisons, 5 values of x).
- Special-framing lemma (Lemma 8.1): re-derived; row spaces computed from Theorem D equal those of the p-conditions at all 311 framings with q <= 2G-2 and 2k0 < q, for G = 3..12.
- Forcing theorem (Thm 9.1): re-derived. A direct test using neither T1 nor Lemma 8.1 forces every one of the 119 tuples in its hypothesis range for G = 5..12. Both ingredients are needed: at G = 8, 0/12 are forced without the polar rows and 3/12 without the two-leaf rows.
- Odd polar orders vanish identically (kappa_k = 0) in every KP propagation run, consistent with T1 and the parity README.
- Zero count (Prop 10.1): Z1 confirmed by brute force over all T for D <= 16. Z1 + Z2 >= G-1 fails only at (G,D) = (4..7, 3) for G <= 1000; margins at G = 8, 9, 10 are 2, 1, 1. (B1), (B2) and the case split re-derived; (B2) checked to D = 6000 and (B1) at G = 300..2000 samples, 0 violations.
- (OL) for G >= 8 and Lemma A for all h (hence Corollary R): the induction (Theorem A; even D by T1/(f); odd D by (OL); Corollary E) uses each hypothesis only where it is available, and there is no circularity. Base G <= 7 from the data, Theorem F and Route E; the referee's pointwise scan closes (OL) at G = 5, 6, 7 independently (mod p).
- Side claim reproduced: TLR nullities in window G-1 with t >= 3 are 2, 3, 0, 0, 0 at G = 3..7; window G with t >= 3 gives 4, 7, 6, 1, 0; Route A's convention gives 11, 11, 7, 0 at G = 4..7.
- TR audit clean: no topological recursion of any kind in the author's chain. The KP hierarchy of the MV generating function (allowed) is the only integrable-systems input; the uniformising coordinate is used only as a change of variables.

## issues

### [1]

#### description

The dependency list of the all-genus statement is incomplete for (I1). Route D's single-vertex polar relation also rests on: the frame-B Givental-Teleman premise (Hodge-vertex graph sum = F_{h,0} of the resolved conifold x C^2); Theorem LC (F1, F3, L1-L3); Theorems A-C of conifold_regularity_note.md (used in Lemma 5.3); and the regularity at t = 0 of the colour-2 single-vertex term (Route D referee, issue 2).

#### location

summary; s.1 INPUT (I1); s.11; obstructions

#### severity

minor

#### suggested_fix

List these inputs explicitly next to (H) = Brini-Schuler Prop. 3.1, Lemma 5.3 and (V).

### [2]

#### description

Lagrange inversion is applied to f = 1/(C u_t), a t-series whose coefficients are Laurent series with positive degree growing with the t-order (f_k has degree k+1). The justification is one line. The claim is correct and was numerically confirmed by the referee (FC to t^13 at 5 framings), but the argument should be written out.

#### location

s.7(b), Lagrange inversion step of the curve identity

#### severity

minor

#### suggested_fix

State the order-by-order formal residue change of variables, or do the step analytically for small t as in s.7(c).

### [3]

#### description

Prop. 10.1 is an exact computer count for G <= 1000 plus analytic bounds for G > 300. (B1) actually holds for all G and D, not only for G > 300.

#### location

claims list ('Zero count', status proved); s.10

#### severity

minor

#### suggested_fix

Label it 'proved (computer-assisted)' and state (B1) in general.

### [4]

#### description

Pointwise forcing is decided by mod-p ranks. 'Nullity 0 mod p implies nullity 0 over Q' is valid, but 'y_D forced mod p' does not by itself imply forcing over Q: a rank drop of the reduced matrix without the y_D column gives a false positive. Side claims only; the main line uses the analytic forcing theorem.

#### location

s.12, pointwise scan 'closes (OL) for 5 <= G <= 12' and the single-level claims

#### severity

minor

#### suggested_fix

Label these as computer-verified evidence, or certify the ranks over Q (exact BigRat) at the framings used.

### [5]

#### description

Not an error, a strengthening. The minimal margin Z1 + Z2 - (G-1) is >= 1 for every 8 <= G <= 1000, and the analytic argument also works with G in place of G-1. So for G >= 8 the Hurwitz-component input 2(b) (Routes B/C) is not needed: the lambda_g input 2(a) alone, which gives G unknowns, suffices. 2(a) itself is needed (margin 1 at G = 9, 10).

#### location

s.11 dependencies (brief 2(b))

#### severity

minor

#### suggested_fix

Add the remark and drop Routes B/C from the dependencies of the G >= 8 step.

### [6]

#### description

The scope of the closed formula for F_G on the slice is stated imprecisely.

#### location

Theorem 6.1 / claims list, 'for every P (in range)'

#### severity

minor

#### suggested_fix

State it as an identity of formal power series, graded by level, for the KP-propagated F_G (any one-leaf data), valid for the true F_G on Taylor coefficients of level <= 2G-3.

### [7]

#### description

'Kflow_t - K_t = log C' should read log(W/C) = log H_1(1-e^-t) - tau^2 t. It is a constant in mu either way, so the conclusion is unaffected.

#### location

Corollary 7.2

#### severity

cosmetic

#### suggested_fix

Correct the formula.

### [8]

#### description

'(s+n)!/(s+n+1-m)' is missing a factorial.

#### location

s.1, proof of (I4)

#### severity

cosmetic

#### suggested_fix

Write (s+n)!/(s+n+1-m)!.

### [9]

#### description

In delta, a(1+tau) is the value of a at lambda = 1+tau, which is easily misread as a product. Lemma 8.1 uses m for two different indices. Theorem 11.1's 'tau' is T = 1 - e^-t, not the framing.

#### location

notation in s.5(iii), s.8, s.1

#### severity

cosmetic

#### suggested_fix

Disambiguate the notation.


## notes

I found no failure in the all-genus argument; verdict sound, with minor and cosmetic issues only. All checks used new independent Perl code (exact modular arithmetic, two primes) in G:\My Drive\Projects\Old\refGW\hodge_vertex_conifold\lemmaA\round2\R1_rigidity\referee2\, with outputs in referee2\out\; no wolframscript and no TR. The author's Mathematica checks (C1)-(C6) in R1_checks.wls are all covered by my independent Perl checks, so running them is optional. The main remaining risk is in the inputs, which are refereed elsewhere: Route A Theorem D for all R; Route D's unreduced single-vertex polar relation and the frame-B and Theorem LC framework behind it; Theorem 11.1. Strengthening: for G >= 8 the Hurwitz-component input brief 2(b) is not needed.

## refuted


## route

R1: rigidity of the KP system (TLR and variants), round 2 -- second independent referee of the all-genus claim (Lemma A for all h)

## verdict

sound

## report_text

# R1 (rigidity of the KP system), round 2: second referee report on the all-genus claim

- **Author's folder:** G:\My Drive\Projects\Old\refGW\hodge_vertex_conifold\lemmaA\round2\R1_rigidity\
- **Referee scripts:** G:\My Drive\Projects\Old\refGW\hodge_vertex_conifold\lemmaA\round2\R1_rigidity\referee2\
- **Referee outputs:** G:\My Drive\Projects\Old\refGW\hodge_vertex_conifold\lemmaA\round2\R1_rigidity\referee2\out\

**How I worked**
- **No TR.** I used or read no topological recursion of any kind: no spectral-curve recursion, loop equations, BKMP, dictionary P2 or TR scripts.
- **Inputs.** The only inputs were:
  - the KP hierarchy of the MV generating function, through Route A's Theorem D (allowed by the user ruling of 2026-10-31, and again in the relayed request);
  - exact algebra and Lagrange inversion.
- **No wolframscript** was run.
- **Code.** All code is new Perl with exact modular arithmetic (p = 10^9+7, second prime 998244353).
  - I wrote it from the *statements* of Route A Theorem D and leaf_structure_proof.md Theorem 11.1.
  - It shares no code with the author's kpmod.pm, kp_point.pl, r1lib.pm, ode_*.pl or t_*.pl.
  - The genus-0 kernel is used in closed form: E0(n;Q) = sigma Pi_n Pi^Q Xi_0(Q+{1,n}), from the genus-0 dispersionless Hirota identity. This differs from the exp-recursion used by the author and by Route A. I re-ran Route A's referee self-test (ref_kp_null.pl selftest, G = 5): 1772 relations, 0 mismatches.

## Verdict: SOUND

I tried to make the argument fail and could not.

- **Main line checked.** Every step was re-derived by hand: T1, flow theorem, curve identity, special-framing lemma, forcing theorem, zero count, then (OL) for G >= 8, then Lemma A for all h.
- **Numerical tests.** Each step was also tested with independent code at G = 7, 8, 9, 10 and beyond. Every test agreed with the author, and sensitivity tests show that the tests do detect errors.
- **Conclusion.** The all-genus statement follows from the project inputs the author lists:
  - Route A Theorems D, F, Corollary E, Theorem A;
  - Route D's single-vertex polar relation;
  - leaf_structure_proof Theorem 11.1;
  - brief 2(a),(b);
  - the exact data for h <= 4.
- **Missing from the list.** Route D's relation itself rests on some framework inputs that should be named (issue 1).
- **Other issues** are minor or cosmetic.

## 1. Steps asserted for all G, and my tests

| # | step (author's section) | my check | range | result |
|---|---|---|---|---|
| 1 | (I4) genus-0 data on the slice, Y_P(z) = W Y_0(z/W) | re-derived by hand: plane-forest count, and [u^n](1/y) = -sigma Pi_n/(n+1); also enters row 4 implicitly (W from the genus-0 equation against an engine whose kernel comes from Xi_0 directly) | -- | correct |
| 2 | (I3) Theorem D at every point of the slice; Lemma 2.1 | own engine rk.pm (Theorem D for all R); compared three different choices of the pair (al,be) per multiset | K = 8, 16 (914 multisets), 26 (11731 multisets) | identical, so the propagation is consistent |
| 3 | Lemma 4.1, branch point d[u^-j]_+(1+tau) = j Pi_j | ref_flow.pl branch | j <= 40, tau = 2/7, -5/13 | 0 mismatches |
| 4 | Theorems 5.1 and 6.1 (flow theorem; F_G on the slice) | engine F_G(P(t)), symbolic in the one-leaf data, against Res Psi [log W + log u(mu) - log u(psi_t(mu))] | 7 random paths (dense and sparse), orders t^6..t^12, tau = 2/7, -3/5, 5/3, -2/7, -1/5, 3 | 816 coefficients, 0 mismatches |
| 5 | residue identity, s.7(c) | re-derived by hand (below); ref_flow.pl res | n <= 8 to t^12, tau = 2/7, -3/11 | 0 mismatches |
| 6 | Theorem 7.1 (FC) | tau log((phi-1)/(l-1)) - (1+tau) log(phi/l) - K_t(phi) is independent of l and equals -log(W e^{tau^2 t}/H_1); phi(0) = 0 and phi(1) = 1 | to t^13, tau = 2/7, -3/5, 3, -1/5, -2/9 | 0 mismatches |
| 7 | T1 | relation_k from the full KP propagation on the leaf curve against the closed form | K = 32 at tau = 2/7 (p = 10^9+7) and tau = -3/11 (p = 998244353), random one-leaf data; symbolic in the one-leaf data at K = 28 (tau = -1/5), K = 24 (2/7), K = 16, 10, 8; K = 24 at tau = -2/9 | 0 mismatches. Covers T1 directly for every G <= 18. Three perturbations (no H_1 normalisation; T = t instead of 1-e^-t; one kernel entry doubled) give 66, 45, 28 mismatches |
| 8 | leaf curve b_j (Theorem 11.1 + string) | ref_leaf.pl rebuilds A_j with its own zeta (s*, U_1) and checks: q^(1), q^(2), q^(3) of s.11.4; the sum rule; b_j = A_j alpha^-j/c^2 = (-1)^j (Pi_j^2/j) t^(j-1) H_j/H_1^j; W = alpha | x = 2/7, -3/5, 5/3, 3, -2/9 | 275 comparisons, 0 mismatches |
| 9 | Lemma 8.1 | regularity rows computed from Theorem D (the engine), not from the y1 y2 divided-difference formula, compared with the p-conditions | every framing -k0/q with q <= 2G-2, 2k0 < q, for G = 3..12 (2, 5, 10, 15, 22, 31, 39, 50, 63, 74 framings; 311 in all) | identical row spaces everywhere |
| 10 | Theorem 9.1 (forcing) | direct test using neither T1 nor Lemma 8.1: polar rows from the full KP propagation on the leaf curve, regularity rows from the engine, y_L = 0 for L < D. Question asked: is y_D forced? | every tuple in the hypothesis range for G = 5..12 (1, 1, 4, 12, 12, 17, 36, 36 tuples; 119 in all) | all forced. Odd polar orders vanish identically in every run (kappa_k = 0). Sensitivity at G = 8: without the polar rows 0/12 forced, without the two-leaf rows 3/12 |
| 11 | Z1(D) (Route A, Cor. E) | brute force: maximum over all T at level D of the forced order (orders of the linear factors of the Pi's), halved at a = 2 | D <= 16 | equal to the closed formula (0,1,2,5,6,11,13,19,22,29,31,44 for D = 1..12) |
| 12 | Prop. 10.1 | direct enumeration for G <= 60 (ref_zeros.pl); inclusion-exclusion count for G <= 1000 (ref_margin.pl); (B1) and (B2) re-derived by hand; (B2) checked for D <= 6000; (B1) plus the full count at G = 300, 301, 500, 777, 1001, 1500, 2000 (all odd D) | as stated | failures only at (G,D) = (4,3), (5,3), (6,3), (7,3); margins at G = 8, 9, 10, 11, 12 are 2, 1, 1, 4, 3, always attained at D = 3; 0 violations of (B1) and (B2) |

**At G = 7, 8, 9, 10 specifically**

| G | T1 orders needed (k <= 2G-4) | Lemma 8.1 framings | tuples forced | margin at D = 3 |
|---|---|---|---|---|
| 7 | 10 | 22 | 4 of 4 | -1 (Z1 = 2, Z2 = 3, need 6) |
| 8 | 12 | 31 | 12 of 12 | 2 |
| 9 | 14 | 39 | 12 of 12 | 1 |
| 10 | 16 | 50 | 17 of 17 | 1 |

- T1 is checked to k = 32, which covers all four.
- G = 7 is not covered by the count, as the author says. It rests on Theorem F or on Route E's M_irr(7).
- My own pointwise scan closes G = 5, 6, 7 independently, using all special framings with q >= D+2, also outside the forcing range (ref_forcing_perD.pl, mod p, evidence-level):
  - G = 5, D = 3: 2 + 6 >= 4;
  - G = 6, D = 3: 2 + 6 >= 5;
  - G = 7, D = 3: 2 + 8 >= 6.

**Why I ran no real-data test.** In the range every Xi_h(R) with h <= 7 vanishes, so every relation is trivially true for the true values. The content of T1, the flow theorem, Lemma 8.1 and the forcing theorem is about KP-propagated values for arbitrary one-leaf data. That is what rows 2-10 test, with formal or random one-leaf data.

**Caveat on my own checks.** Mod-p agreement is evidence, not proof: the false-positive probability is about 1/p per coefficient. The proofs themselves were checked by hand (section 2).

## 2. Re-derivation (where I tried to break it)

### s.2 (Lemma 2.1)
- The coefficient [t^r z^-a w^-b] of (H) at P(t) is a finite combination of Theorem D's coefficients at (R, a, b) with D(R) <= r, since P_j = O(t^(j-1)). It is therefore valid when r + a + b - 2 <= 2G - 3.
- Coefficients with a = 1 or b = 1 are just the identities G_1n = V_n.
- (S) and (E) are exact.
- At a fixed level, the dependencies strictly lower the t-order. So the formal solution exists and is unique, and the true values agree with it up to level 2G-3. Correct.

### s.3
- Lemma 3.1 re-derived: A = (f o z) z'; Res A/(l_w - l) = -A_+(l_w); [A_- hatV]_+ = 0.
- (S'), (E') and (P') re-derived, including:
  - d_t Pcal|_l = d_t T|_z - (W'/W) h d_l Pcal;
  - b = N - (W'/W) a;
  - theta U = h d hatU.

### s.4
- [g]_+(1+tau) = j Res u^-j/(l(l-1)) = j binom(j tau + j - 1, j - 1).
- Hence N(1+tau) = -W' + W' = 0. Correct.

### s.5
- (i) [u^-j]_+(1) - [u^-j]_+(0) = binom(j tau + j, j) - binom(j tau + j - 1, j) = Pi_j.
- (ii) h = -(l+tau) - sigma/(l-1-tau) and N(1+tau) = 0 make hN a polynomial. Then d_t phi = -hN; the constant is fixed at l = 0 because h(0) = 0.
- (iii) In delta = W sigma c_R/(a(1+tau) (l-1-tau)), the symbol a(1+tau) is the *value* of a at l = 1+tau, so delta = O(1/l). With this, (S') and (E') check term by term, and uniqueness follows from the invertible Toeplitz operator [a.]_-. Correct.

### s.6
- d_t psi_t = -v_t(psi_t).
- Residues are invariant under mu = phi_t(l), and Res(R b) = 0. Correct.

### s.7
- **(a)** With zeta~ := tW/H_1, Phat_j = (-1)^j (Pi_j/j) t^-1 zeta~^j H_j for all j >= 1; the j = 1 term is -W.
  - FC, and hence T1, uses only the leaf ratios b_j = A_j alpha^-j/c^2.
  - The function zeta(t) of Theorem 11.1, which involves sqrt(1 - 4xT/c^2), cancels exactly, because C_j/c^2 = Pi_j^2/j.
  - So T1 does not depend on the sum rule.
- **(b)** Lagrange inversion: see issue 2.
- **(c)** The residue identity re-derived:
  - substitute Lambda = e^(t mu), then Lambda = 1 + (e^t - 1)s;
  - expand (1 + (e^t-1)s)^(-n tau - 1) binomially; the coefficients are (-1)^(n-1) binom(n tau + n - 1 + i, n - 1 + i) binom(n tau, i);
  - this gives (-1)^(n-1) t^(n-1) e^(-n sigma t) Pi_n 2F1(n + n tau, -n tau; n; 1 - e^t);
  - Pfaff (DLMF 15.8.1), with (1-z)^(-b) = e^(n tau t) and z/(z-1) = 1 - e^-t, gives e^(-n tau^2 t) Pi_n H_n(T).
  - The analytic-continuation step (real tau in (0,1), small t, contour |mu| = r with 1 < r < 2 pi/|t|) is legitimate.
- **(d)** The polynomiality argument for eps = 0 is correct.

### s-form
- Res Psi s/(1 - s l)^2 = -(1/s) Psi'(1/s).
- h(1/s) = -(1-s)/(s(1-(1+tau)s)).
- h(1 + 1/s) = -(1+s)/(s(1 - tau s)).
- Correct.

### s.8 (Lemma 8.1)
- [u^(qm)] g = (1/(qm)) [y^(qm-1)] g'(y)(1-y)^(m k0).
- Sum_r (-1)^r binom(m k0, r)(qm - r) y^r = mq (1-y)^(m k0 - 1)(1 - (1+tau0) y).
- The y1-transform is triangular with unit diagonal, and (m-1) k0 + q - 1 < qm. Correct.
- Xi_G(m,n) is finite at sigma0 != 0 for *every* R, because sigma^G Xi_G is a polynomial. So the conditions are legitimate.

### s.9 (Theorem 9.1)
- q R(1/l) = l X(l) with X = A r(l) + B r(l+1); checked term by term.
- The even part of X is O(l^-2G).
- The denominator of E divides (l(l^2-1))^A because B < A, and 3A < 2G gives E = 0.
- Matching the poles at +1 and -1 gives ord_0 pi >= A - B.
- deg s1 <= 2B - D - 2 < B, so s1 = 0.
- A != B, so pi_1 = 0. Correct.
- The proof actually gives more than stated: F = 0 mod y^(2G-1) at tau0, i.e. every y_L(tau0) with L <= 2G-3 vanishes.

### s.10
- **Multiplicities.** For a >= 3, sigma is a local coordinate at -c/a (1 + 2 tau0 != 0), and each pair {c, a-c} is counted once. At a = 2 the sigma-order is half the tau-order D - 1.
- **Distinctness.** All points are distinct: x -> x(1-x) is injective on (0, 1/2], and the reduced denominators differ (<= D+1 for Z1, >= D+2 for Z2).
- **(B1)** holds for all G and D. The counts for k0 = 1..4 are g-D, >= (g-D)/2, >= (2(g-D)-2)/3 and >= (g-D-2)/2, including the edge cases D = 3, 5 where the condition q > 2k0 matters.
- **(B2)** follows from:
  - #coprime pairs >= n^2 - Sum_{d>=2} floor(n/d)^2 >= (2 - pi^2/6) n^2;
  - m_a >= #{j >= 2 : a - 1 <= D/j}.
- **Case split.** The split at D + 1 = (7G-6)/24 is right, and so is the bound D >= 80.

## 3. Logic: hypotheses, base cases, circularity, brief s.3

**Chain at genus G >= 8, by inner induction on the level D**
- **D = 1:** Theorem A.
- **Even D:** with y_L = 0 for L < D, T1 at k = D gives relation_D = kappa_D y_D (beta_{D,D} = kappa_D). Then (I1) gives y_D = 0, and Pi_{D+1} is not identically 0.
- **Odd D >= 3:**
  - The forcing theorem gives Z2(G,D) simple zeros. It needs only Lemma A in genera < G (for T1, (I1) and Theorem D) and y_L = 0 for L < D.
  - Corollary E gives Z1(D) zeros. It needs Lemma A at genus G for levels < D, which is the inner induction hypothesis.
  - Brief 2(a),(b) give sigma^G Xi_G(D+1) = sigma Q with deg Q <= G-2.
  - Z1 + Z2 >= G - 1 at nonzero, distinct points, so Q = 0.
- **Then** Corollary E gives every T at level D.

**Hypotheses.** Each hypothesis is used only where it is available. Nothing at genus G and level >= D is assumed. Lemma 5.3 needs only genera < G, because a genus-G stable graph without legs that has an edge has no genus-G vertex.

**Base**
- h <= 4: exact data.
- h = 5, 6, 7: Route A's Theorem F (computer-verified; the mod-p rank-lifting argument is rigorous).
- Independently, Route E's direct engine gives M_irr(5), (6) and (7) complete.
- My pointwise scan is a third, evidence-level confirmation at G = 5, 6, 7.

**Circularity: none**
- (H) is Brini-Schuler Prop. 3.1 (degeneration plus Bryan-Pandharipande). It does not depend on Lemma A or on Corollary R.
- Corollary R at genus g needs M_irr(h) only for h < g (for the edge graphs) and the known polar part (for the single-vertex graph).
- So deriving Lemma A from (H) and then Corollary R from Lemma A is not circular.

**(OL) implies Lemma A (brief s.3)**
- The author's s.11 is exactly the brief s.3 induction: even D by (f) (equivalently T1 at k = D), odd D by (OL), then Corollary E.
- The final statement gives (OL) for every G >= 8 and every odd D in [3, 2G-3].
- With G <= 7 from Theorem F, the data and Route E, (OL) holds for all G >= 3.
- Hence Lemma A holds for all h, and Corollary R holds in all genera.

**Robustness remark (not needed by the author, but it removes a dependency)**
- The minimal margin Z1 + Z2 - (G-1) is >= 1 for every 8 <= G <= 1000 (ref_margin.pl, own inclusion-exclusion count).
- The analytic argument also works with G in place of G-1: (B1) covers D + 1 <= (7G-15)/24; otherwise D >= 87, and Z1(D) >= (24D + 39)/7 for D >= 81.
- So for G >= 8 the Hurwitz-component input, brief 2(b) (Routes B and C), is dispensable. The lambda_g input 2(a) alone puts sigma^G Xi_G(D+1) in span(sigma..sigma^G), with G unknowns and at least G zeros.
- 2(a) cannot be dropped: the margin is only 1 at G = 9 and G = 10.

## 4. Inputs and literature, with their exact hypotheses

**Route A, Theorem D (A_mv\report.md; referee verdict sound)**
- Hypotheses: G >= 2; Lemma A in all genera < G; R with parts >= 2; al, be >= 2; excess >= 1.
- The author's "Theorem D at every point P of the slice" is exactly the generating function over R of these relations: the genus-0 factor is exp(D(z)D(w)F_0)(P) = (Y_P(z)-Y_P(w))/(z-w). Lemma 2.1 uses only coefficients of level <= 2G-3.
- I re-checked that at order hbar^(2G) the shift operators of the Fay identity contribute F_G only through D(z)D(w)F_G and d_1 D(z) F_G. So the F_G-linear part is exactly the linearised dispersionless equation.
- All other terms (lower-genus products and pure genus-0 dispersive terms) vanish at excess >= 1 by excess homogeneity, given Lemma A in genera < G.
- Literature behind it:
  - Liu-Liu-Zhou, J. Diff. Geom. 65 (2003), and Okounkov-Pandharipande, Geom. Topol. 8 (2004), for the MV formula;
  - Takasaki 1984 and Orlov-Scherbin, TMP 128 (2001), for the KP tau-function property;
  - Takasaki-Takebe, Rev. Math. Phys. 7 (1995), for the dispersionless Hirota equations.
- Allowed by the user ruling.

**Route A, Corollary E, Theorem A, Theorem F**
- Corollary E: needs Lemma A in genera < G and at genus G for levels < D; used only for Z1, under (OL)'s hypotheses.
- Theorem A: unconditional.
- Theorem F: 5 <= G <= 20, computer-verified.

**(I1): Route D's single-vertex polar relation, unreduced form (brief 2(e))**
- The form used is [t^k] F_G(P(t)) = 0 for 1 <= k <= 2G-3, on the curve P_j(t) = (-1)^j (Pi_j/j) t^(j-1) H_j(1-e^-t)/H_1^j.
- This matches T_k(G) at leading order (weights prod (-1)^j eps_j/j).
- It also matches Theorem 11.1 after string elimination of the lam_1 leaves: b_j = A_j alpha^-j/c^2 with alpha = -A_1/c^2 (row 8). In Theorem 11.1 the symbol "tau" denotes T = 1 - e^-t, x is the conifold weight and the framing is -x/c.
- It rests on:
  - (H) = Brini-Schuler Prop. 3.1 = Thm A1. Its hypothesis alpha_i = eps_j holds on the project torus (Route D referee).
  - Lemma 5.3: needs Lemma A in genera < G, Theorems A-C of conifold_regularity_note.md, and Theorem LC.
  - (V): the colour-1 term, with the colour-2 single-vertex term regular at t = 0 (F1, F3).
  - The frame-B Givental-Teleman premise.

**Other inputs**
- Brief 2(a): Faber-Pandharipande lambda_g formula. Brief 2(b): Route C Theorem H and Route B Theorem B; dispensable for G >= 8 (section 3).
- **Pfaff transformation.** 2F1(a,b;c;z) = (1-z)^(-b) 2F1(c-a,b;c;z/(z-1)) (DLMF 15.8.1), here with a = n + n tau, b = -n tau, c = n >= 1 and z = 1 - e^t. It holds as an identity of power series in t, since z and z/(z-1) are both O(t).
- **Standard formal tools.** Lagrange inversion (one formal variable); invariance of formal residues under l -> l + (a t-small polynomial).
- **Elementary counting.** #coprime pairs in [1,n]^2 = 2 Phi(n) - 1 = Sum_d mu(d) floor(n/d)^2.

## 5. Issues

| # | severity | location | issue | suggested fix |
|---|---|---|---|---|
| 1 | minor | summary; s.1 (I1); s.11; obstructions | The dependency list of the all-genus statement is incomplete for (I1). Route D's relation also rests on the frame-B Givental-Teleman premise (Hodge-vertex graph sum = F_{h,0}), Theorem LC (F1, F3, L1-L3), Theorems A-C of conifold_regularity_note.md (used in Lemma 5.3), and the regularity of the colour-2 single-vertex term (Route D referee, issue 2). | List these explicitly. |
| 2 | minor | s.7(b) | Lagrange inversion is applied to f = 1/(C u_t), a t-series whose coefficients are Laurent series with positive degree growing with the t-order (f_k has degree k+1). It is justified in one line. | Write out the order-by-order formal residue change of variables, or do it analytically for small t as in (c). FC itself is confirmed to t^13 at 5 framings. |
| 3 | minor | claims list ("Zero count", status "proved"); s.10 | The proof is an exact computer count for G <= 1000 plus analytic bounds for G > 300. Also, (B1) holds for all G and D, not only for G > 300. | Label it "proved (computer-assisted)" and state (B1) in general. |
| 4 | minor | s.12, pointwise scan "closes (OL) for 5 <= G <= 12"; single-level claims | "Nullity 0 mod p implies nullity 0 over Q" is valid. But "y_D forced mod p" does not by itself imply forcing over Q: a rank drop of the reduced matrix without column D gives a false positive. This affects side claims only. | Label as computer-verified evidence, or certify the ranks over Q at the framings used. |
| 5 | minor (suggestion) | s.11 dependencies | Brief 2(b) is not needed for G >= 8 (section 3). | Add the remark, and drop Routes B/C from the dependencies of the G >= 8 step. |
| 6 | minor | Theorem 6.1, "for every P (in range)" | The scope is imprecise. | State it as an identity of formal power series, graded by level, for the KP-propagated F_G (any one-leaf data), valid for the true F_G on coefficients of level <= 2G-3. |
| 7 | cosmetic | Cor. 7.2 | "Kflow_t - K_t = log C" should read log(W/C) = log H_1(1-e^-t) - tau^2 t. It is a constant either way, so the conclusion is unaffected. | Correct the formula. |
| 8 | cosmetic | (I4) | "(s+n)!/(s+n+1-m)" is missing a factorial. | Write (s+n)!/(s+n+1-m)!. |
| 9 | cosmetic | notation | In delta (s.5(iii)), a(1+tau) is the value of a at lambda = 1+tau, which is easily misread as a product. Lemma 8.1 uses m for two different indices. Theorem 11.1's "tau" is T = 1 - e^-t. | Disambiguate. |

## 6. Side claims (s.12)

- **TLR nullities.** Reproduced exactly with my own engine (ref_tlr.pl):
  - window G-1, t >= 3: 2, 3, 0, 0, 0 at G = 3..7;
  - window G, t >= 3: 4, 7, 6, 1, 0 at G = 3..7;
  - Route A's convention (window G, t >= 2): 11, 11, 7, 0 at G = 4..7.
  - So in the author's window TLR already holds at G = 5 (checked for G <= 7).
  - Route A's TLR conjecture is stated with c = 3..2G-2 but was computed with c = 2..2G-2, which explains the different numbers.
- **Not re-run in full:** the single-level nullities and the G <= 12 scan. My per-level scan agrees that (OL) closes at G = 5, 6, 7.
- **Author's refuted ansatze (s.13).** Consistent with my row 6: with alpha in place of C the identity fails at t^2, because log H_1 - tau^2 t = O(t^2) is not zero.

## 7. Files (in referee2\)

**Scripts**
- rk.pm: modular library; Theorem D engine; leaf curve.
- ref_t1.pl: T1 test and consistency (`cons`) mode.
- ref_t1_pert.pl with rk_pert.pm: sensitivity runs.
- ref_flow.pl: modes branch, res, curve, path.
- ref_leaf.pl: leaf-curve checks.
- ref_forcing.pl: Lemma 8.1 and the forcing test.
- ref_forcing_sens.pl, ref_forcing_perD.pl: variants of the forcing test.
- ref_zeros.pl, ref_margin.pl, ref_bounds.pl: zero count and bounds.
- ref_tlr.pl: TLR side claims.

**Outputs (in out\)**
- T1: t1_K32_2_7.txt, t1_K32_m3_11.txt, t1_K28_m1_5_sym.txt, cons_K26_5_3.txt.
- Forcing: forcing_G8..G12.txt, forcing_summary.txt, forcing_perD_G5to7.txt.
- Zero count: zeros_G3to60.txt, margin_G400.txt, margin_G1000.txt, bounds_D6000_G2000.txt.
- Everything else: misc_checks.txt (branch, res, curve, path, leaf, sensitivity, TLR).

