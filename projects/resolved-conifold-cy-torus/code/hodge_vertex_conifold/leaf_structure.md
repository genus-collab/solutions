# Leaf structure of the t-expansion of the singular column (2026-10-30)

Task: second ingredient of the all-genus argument. For k = 0..3 decompose the order-t^k coefficients of
 (i) rho q_t - q_0,  (ii) rho f_{j,b} (legs, j = 1, 2, external power b),  (iii) rho^2 E11(psi, psi')  (loops / singular-singular edges)
into leaves lam_j(psi) = 1/(c - j psi), j = 0..k+1 (products lam_i(psi) lam_l(psi') for loops), fitting from the low psi-orders and
checking every remaining available psi-order (genus-4 factors hvc_factors_g4_Bm.mx, z^9 truncation).
Conventions: rc / hvc_lib.m, mode B, eps1 = x, eps2 = 1, c = 1 + x, singular root U1, t = Log Q,
rho(t) = (Delta_1 t^2/(-x c^3))^(-1/2), q_t(psi) = Sum_{d>=1} Dil_1(d; t) psi^d = -psi r_{0,1}(-psi; t)  (Dil_1(1) = -1 is the dilaton shift),
q_0 = c - c^2 lam_1 = c^2 (lam_0 - lam_1),  f_{j,b}(psi) = Sum_a E_{1j}(a, b) psi^a.
Script: hvc_leaf_structure.wls.

## Log

Step 0 (hvc_leaf_probe0.wls): factor file loads in 0.1 s; ZMAX = 9, so Dil_1(d) for d = 1..10 and E_{1j}(a, b) for a + b <= 8.
Delta_1 = -x c^3/t^2 [1 + (1-x)^2/c^2 t + (5 - 4x - 6x^2 - 4x^3 + 5x^4)/(12 c^4) t^2 + (1-x)^2/(12 c^2) t^3 + O(t^4)]  [exact].
Hence rho = 1 - (1-x)^2/(2c^2) t + O(t^2). The known O(t) row correction r^(1)(z) = (1-x)^2 z/(c (c+z)(c+2z)) gives
q^(1) = (1-x)^2/2 (lam_0 - 2 lam_1 + lam_2), and rho_1 q_0 = -(1-x)^2/2 (lam_0 - lam_1): the lam_0 parts cancel at order t
(by hand; checked by the script below).

## Step 1: leaf decompositions (hvc_leaf_structure.wls -> hvc_leaf_structure_out.txt, data in hvc_leaf_structure_data.m; 9 s)
Fit from psi-orders 0..k+1, then EVERY remaining available psi-order checked exactly (rational functions of x).
rho = 1 - (1-x)^2/(2c^2) t + (1 - 8x + 15x^2 - 8x^3 + x^4)/(6c^4) t^2 - (1-x)^2 (1 - 20x + 60x^2 - 20x^3 + x^4)/(24 c^6) t^3 + O(t^4).

(i) rho q_t - q_0 = Sum_k t^k delta_k  [checks: 9, 8, 7, 6 extra psi-orders at k = 0..3, all exact]
  k = 0: 0
  k = 1: ((1-x)^2/2) (lam_2 - lam_1)
  k = 2: [-(2 - 4x + 3x^2 - 4x^3 + 2x^4) lam_1 + 6 (1-x)^4 lam_2 - (2-x)^2 (1-2x)^2 lam_3] / (12 c^2)
  k = 3: (1-x)^2 [-(3 + 5x^2 + 3x^4)/72 lam_1 + (7 - 20x + 22x^2 - 20x^3 + 7x^4)/24 lam_2 - (2-x)^2 (1-2x)^2/8 lam_3 + (3-x)^2 (1-3x)^2/36 lam_4] / c^4
  NO lam_0 component at any order (the rho normalisation removes it exactly; sum of the coefficients = 0 because q_t(0) = 0).

(ii) legs rho f_{j,b}: for j = 1 and j = 2, every b with at least one check (7 - b - k extra psi-orders: b <= 6 - k, i.e. b <= 6, 5, 4, 3 at k = 0..3, same for j = 1, 2), the order-t^k part is Sum_{j'=1}^{k+1} (coef) lam_j', all checks exact; lam_0 coefficient = 0 in EVERY case (also in the
  0-check fits). t^0: E11-leg = -(x^3/2) c^(-b-1) lam_1, E12-leg = chi_b lam_1 (rank one, as in F3).
  Examples (b = 0): E11-leg t^1: (1-x)^2 x^3/c^3 [-(1/12) lam_1 + (2/3) lam_2];  E12-leg t^1: (1-x)^2 x^3/c^3 [-(1/2) lam_1 + 2 lam_2].
  E12 legs: the lam_1 coefficient at order t^k is chi_b times a b-independent factor ((1-x)^2/(2c^2) at k = 1, (2 - 4x - 3x^2 - 4x^3 + 2x^4)/(12 c^4)
  at k = 2), i.e. rank one for k = 1, 2; this FAILS at k = 3 (b = 0 vs b = 1 ratios differ: 89x^3 vs -55x^3 in the sextic), so it is not structural.

(iii) loops rho^2 E11(psi, psi') = Sum_k t^k Sum_{i,l} e^(k)_{il} lam_i(psi) lam_l(psi')  [checks: 41, 36, 29, 20 coefficient pairs (a, a'), a + a' <= 8, all exact]
  k = 0: e_11 = -x^3/2.
  k = 1: e_11 = -(1-x)^2 x^3/(2c^2), e_12 = e_21 = 2(1-x)^2 x^3/(3c^2), e_22 = 0.
  k = 2: e_11 = -x^3 (7 - 20x + 21x^2 - 20x^3 + 7x^4)/(24 c^4), e_12 = e_21 = (1-x)^4 x^3/c^4, e_13 = e_31 = -3(2-x)^2(1-2x)^2 x^3/(16 c^4), e_22 = -(1-x)^4 x^3/c^4; e_23 = e_33 = 0.
  k = 3: nonzero exactly for (i, l) in {11, 12, 13, 14, 22, 23} (and transposes); e_24 = e_33 = e_34 = e_44 = 0.
  So the support is exactly {i, l >= 1, (i-1) + (l-1) <= k}: total degree D <= k for the PAIR, which is what a singular-singular edge needs. No lam_0 row/column.
Top coefficients carry the Gauss-hypergeometric factors prod_{j=1}^{k} (k+1-j - j x)^2: (1-x)^2, (2-x)^2(1-2x)^2, (3-x)^2(1-x)^2(1-3x)^2, i.e.
the residues of (b+1)_k^2/(k!(a+b+2)_k) at z = -c/(k+1) (section 6 of conifold_regularity_note.md).

## Step 2: bookkeeping (stage 2 of hvc_leaf_structure.wls, run with argument "book"; 13 s in total)
Leaf+Xi prediction: write rho q_t = q_0 + Sum_{k>=1} t^k delta_k, rho f = Sum_k t^k f^(k), rho^2 E11 = Sum_k t^k E^(k) with the fitted leaf coefficients;
the coefficient of t^(2-2h+k) of V_h = (-x c^3/t^2)^(h-1) d_{rho f_1}..d_{rho f_n} [loops rho^2 E11] F_h(q_0 + delta) is then
(-x c^3)^(h-1) Sum (leaf coefficients) Xi_h(multiset), the multisets read off as monomials of
prod_legs (Sum t^k f^(k)) (Sum t^k E^(k))^loops exp(Sum t^k delta_k)  (lam_j -> commuting symbol LL[j]).
Xi_h evaluated exactly (lam_1 leaves stripped by the proved string identity Xi(1^p, R) = prod_{i<p} (J_R + i)/c^2 Xi(R); sanity-checked against
direct summation for (1,2), (1,1,2), (1,3), (1,2,2) at genus 2 and 3).
Direct: the vertex functional summed from the factor t-series themselves (Dil_1(d; t), E_1j(a, b; t), Delta_1(t)^(h-1); NO leaf decomposition).

Results (exact in x; every order agrees between leaf+Xi and direct):
 genus 2, coefficients of t^-2, t^-1, t^0, t^1:
  bare vertex: {gap_2, 0, -(7 + 40x + 190x^2 - 280x^3 - 184x^4 - 280x^5 + 190x^6 + 40x^7 + 7x^8)/(69120 x c^4), (1-x)^2 x (5 - 61x + 5x^2)/(3456 c^4)}
  E11-leg b = 0; E12-leg b = 0; E12-leg b = 1; one loop; E11-leg + E12-leg: {0, 0, nonzero, nonzero} in each case.
  => the leaf/Xi bookkeeping reproduces the REGULAR orders t^0, t^1 too (these involve Xi_2 at D = 2, 3, outside Conjecture M), a non-trivial
     check of the rho rescaling, the homogeneity and the decompositions.
 genus 3, coefficients of t^-4 .. t^-1 (first genus-3 test of Conjecture V):
  bare vertex: {gap_3, 0, 0, 0}, gap_3 = -(1 + x^2)(31 + 18x^2 + 31x^4)/(161280 x);
  E11-leg b = 0, 1; E12-leg b = 0, 1; one loop; two E12-legs (b = 0); E11-leg + E12-leg (b = 0): {0, 0, 0, 0} in every case.
 The maximal leaf degree occurring at order t^k is exactly D = k (k = 0..3): the bound D <= k is sharp.
 Truncation: these genus-3 configurations need psi-orders <= 8 on the legs/loops and <= 9 on the dilaton leaves, all inside the z^9 factors,
 so the genus-3 direct results are complete, not truncated.
 Explicit Conjecture M instances used (after stripping lam_1 by the string identity): Xi_2(2) = 0; Xi_3 at {2}, {3}, {4}, {2,2}, {2,3}, {2,2,2} = 0
 (printed by the "book" run); the first D = 2h-2 ones are nonzero (Xi_2(3), Xi_2(2,2); Xi_3(5), Xi_3(2,4), Xi_3(3,3), Xi_3(2,2,3)).

## Step 3: orders t^4, t^5 (genus-4 range; run "kt5", output hvc_leaf_structure_out_kt5.txt, data hvc_leaf_structure_data_kt5.m; 16 s)
Loops now with the RESTRICTED symmetric ansatz Sum_{1<=i<=l, i+l-2<=k} e_il (lam_i lam_l' + lam_l lam_i'), solved from all 45 ordered pairs a + a' <= 8
(overdetermined): consistent at every k = 0..5 (44, 43, 41, 39, 36, 33 extra equations satisfied), so E11 is symmetric with support D <= k and
no lam_0. (The full (k+2)^2 ansatz needs (a, a') up to (k+1, k+1), impossible beyond k = 3 at z^9.)
 dilaton k = 4, 5: lam_0 coefficient 0, leaves lam_1..lam_{k+1}, 5 resp. 4 extra psi-orders exact.
 legs k = 4: b = 0, 1, 2 (3, 2, 1 checks), k = 5: b = 0, 1 (2, 1 checks), j = 1 and 2: lam_0 coefficient 0, all checks exact.
 rho to O(t^5): coefficient of t^4 = (1 - 52x + 493x^2 - 1624x^3 + 2365x^4 - 1624x^5 + 493x^6 - 52x^7 + x^8)/(120 c^8).
Top coefficients (k = 1..5): dilaton a_{k+1}^(k) = (-1)^(k+1) prod_{j=1}^k (k+1 - j c)^2 / (k! (k+1)! c^(2k-2))
 [P_k = prod (k+1-jc)^2 = (1-x)^2, (2-x)^2(1-2x)^2, 4(3-x)^2(1-x)^2(1-3x)^2, ...], and the E11 corner e_{1,k+1}^(k) = ((k+1)^2/(k+2)) (x^3/c^2) a_{k+1}^(k)
 (4/3, 9/4, 16/5, 25/6, 36/7): the hypergeometric residues (b+1)_k^2/(k!(a+b+2)_k) at z = -c/(k+1), as expected from section 6 of the note.
Also checked: E21(b, a) = E12(a, b) and E11(a, b) = E11(b, a) for a + b <= 8 (so a singular vertex at either end of an edge sees the same leg).

## Step 4: conclusions

Lemma L (leaf structure; VERIFIED exactly in x at every available psi-order, k <= 5, not proved). For every k >= 0
 (L1) rho q_t - q_0 = Sum_{k>=1} t^k delta_k with delta_k in span(lam_1, ..., lam_{k+1});
 (L2) for every leg type (j, b): [t^k] rho f_{j,b} in span(lam_1, ..., lam_{k+1});
 (L3) [t^k] rho^2 E11(psi, psi') in span{lam_i(psi) lam_l(psi') : i, l >= 1, (i-1) + (l-1) <= k}.
 In words: no lam_0 component anywhere, no polynomial remainder, no higher poles, and total leaf degree D <= k (sharp).
 Equivalent large-z form of the lam_0-freeness: lim_{z->oo} z r_{0,1}(z; t) = c/rho(t) order by order in t, and E_{1j}(z, w) -> 0 as z -> oo
 at every order in t (fixed power of w). The rho normalisation is exactly what kills the lam_0 part of the dilaton input (at O(t):
 q^(1) = ((1-x)^2/2)(lam_0 - 2 lam_1 + lam_2) and rho_1 q_0 = -((1-x)^2/2)(lam_0 - lam_1)).

Proposition (bookkeeping; a proof, given its inputs). Let h >= 2. Assume L1-L3 through order t^(2h-3) and Xi_h(R) = 0 for every multiset R of
integers >= 2 with D(R) <= 2h-3 (the part of Conjecture M that survives the string reduction). Then Conjecture V holds at genus h.
 Proof. By homogeneity V_h = (-x c^3/t^2)^(h-1) d_{rho f_1}...d_{rho f_n} [rho^2 E11]^loops F_h(q_0 + delta). Expanding, [t^(2-2h+k)] V_h is a finite
 combination of Xi_h(j_1..j_r) with all j_i >= 1 (L1-L3) and D <= k (degrees add). Strip the lam_1's with the proved string identity
 Xi(1^p, R) = prod_{i<p} (J_R + i)/c^2 Xi(R): if R is empty the factor with i = 0 is J_R = 0 (Theorem C), otherwise Xi_h(R) = 0 by assumption as
 D(R) = D <= k <= 2h-3. So every polar coefficient vanishes, except the single term with r = 0, which occurs only for n = 0, no loops, k = 0:
 (-x c^3)^(h-1) Xi_h() = gap_h (F4). QED
 Graph form: in a coloured stable graph with an edge every singular vertex v of genus >= 2 has r_v >= 1 leaves, all with j >= 1, and an E11 edge
 between two singular vertices splits its D <= k_e between its ends (L3); so a non-zero term needs D_v >= 2h_v - 2 at each such v, the total
 t-order is >= Sum_v (2 - 2h_v) + Sum_v D_v >= 0, and the graph is regular (this is the Reduction Lemma of the note in leaf language;
 E21(b,a) = E12(a,b) checked).

Status at genus 2 and 3: the Conjecture M input needed is Xi_2(2) = 0 and Xi_3 = 0 at {2}, {3}, {4}, {2,2}, {2,3}, {2,2,2}: all exact computations.
 Hence Conjecture V at genus 2 and 3 is IMPLIED by Lemma L through t^1 resp. t^3, which is verified exactly at all available psi-orders.
 NO OBSTRUCTION: no lam_0 components, no remainders, D <= k everywhere. Direct confirmation (Step 2): genus-3 vertex polar parts {gap_3, 0, 0, 0}
 (bare) and 0 (7 configurations with legs/loops), and the leaf+Xi expansion reproduces the genus-2 vertex through t^1, including its regular
 t^0, t^1 coefficients.
What is NOT proved:
 (G1) Lemma L itself: verified to z^9 and k <= 5; a proof needs the closed form of the singular column (Euler integral / Gauss 2F1,
      note section 6 step 1). Range: dilaton to psi^10, E_{1j}(a, b) for a + b <= 8. Since every edge factor in a stable graph of genus g has
      a + b <= 3g - 4 (Sum_v dim_v = 3g - 3 - E), these are all the coefficients that occur in graphs of genus <= 4 (where a fit has fewer than
      k + 2 orders, a lam_0-free leaf representation of the available coefficients exists trivially). For the ABSTRACT vertex functional
      (arbitrary external power b, many half-edges) leg orders up to psi^(3h - 3 + n_half) are needed, beyond the data when b + 3h - 3 + n_half > 8.
 (G2) Conjecture M in general (only the instances above are needed at genus 2, 3).
 (G3) Genus 4: Lemma L is verified through t^5 (the full range k <= 2h-3 = 5), but the needed Conjecture M instances with r >= 3
      ({2,2,2}, {2,2,3}, {2,2,4}, {2,3,3}, {2,2,2,2}, {2,2,2,3}, {2,2,2,2,2}) are not yet computed (only r <= 2 at genus 4).

## Step 5 (2026-10-31): Lemma L proved to all orders
leaf_proof/leaf_structure_proof.md proves Lemma L, with F1 and F3, for every k >= 0 and at every psi-order, without topological recursion
(Theorem LC; summary in conifold_regularity_note.md section 7b). So (G1) is closed. The tables of Steps 1-3 are reproduced by the closed forms
    rho q_t = c + Sum_j A_j(t) lam_j,   A_j = (-1)^j C_j t^-1 zeta^j 2F1(j x/c, j x/c; j; tau)   (leaf_structure_proof.md Theorem 11.1),
    e_il(t) = (i l/(c^2 t^2)) Int_0^t Theta (t A_i)' (t A_l)'   (Theorem 11.4),   and the legs (9.7),
checked against hvc_leaf_structure_data*.m by hvc_leaf_closed_forms.wls. The top-coefficient formula of Step 3 is Corollary 11.3, now proved for
all k. Two qualifications of the observations above: the rank-one pattern of the lam_1 coefficients of the E_12 legs at k = 1, 2 is not
structural and fails at k = 3 (leaf_structure_proof.md s.11.3); "denominator exactly prod (c + j z)" holds for generic x only (s.5.5). (G3) is
settled by hvc_multiprobe_g4b.wls (proof_status_mixed_scale.md, Addendum (1)).
