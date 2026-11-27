# Referee report 1 on proof_M_residue.md (logical structure and cited results)

Date: 2026-10-31. Referee 1 (skeptical). Notes appended as the review proceeds.

## 0. Files read
proof_M_residue.md (all sections 0-9a), mixed_scale_findings.md, hvc_multiprobe.wls, hvc_multiprobe_out.txt,
hvc_multiprobe_g4.wls, hvc_multiprobe_g4_out.txt, hvc_which_class.wls, hvc_which_class_out.txt,
check_closed_formula_r1.pl, proof_M_hurwitz.md, proof_M_framing.md, conifold_regularity_note.md sections 1-7.

## 1. Lemma 1 (dilaton leaves -> K = exp(-Sum kappa_k/k)) and Corollary 1
- Pull-back of psi_i: correct (one forgotten point at a time; psi_{n+1} D_{i,n+1} = 0 already needs only psi^{>=1}
  on the forgotten point, psi^{>=2} is more than enough).
- Exponential formula: Sum_m 1/m! Sum_{sigma in S_m} prod_cycles W = exp(Sum_l W_l/l), W_l = Sum_k kappa_k [z^k] t(z)^l,
  t(z) = T(z)/z = -z/(1-z): exponent = Sum_k kappa_k [z^k](-log(1-t)) = Sum_k kappa_k [z^k] log(1-z) = -Sum kappa_k/k. CORRECT.
- Corollary 1: kappa_k = pi^*kappa_k + psi^k on Mbar_{h,1} (AC), so K = pi^*K_h (1-psi);
  (1-psi)/(1-j psi) = 1 + (j-1) Sum_b j^b psi^{b+1}, pi_* psi^{b+1} = kappa_b. CORRECT. (Gives Theorem C for r = 1 at once.)
- Homogeneity degree -(3h-3+2r) re-derived: CORRECT, so c = 1 is legitimate.

## 2. Lemma 2 (Lambda^vee(1) K = c(-R pi_* omega_log^{-1}))
- Chiodo's formula at r = 1 (= Mumford GRR with markings), AC kappa convention: the boundary term depends only on the
  node multiplicity q = 0, the psi term only on a_i, so s: 0 -> -1 changes only the kappa part by
  [B_{m+1}(-1)-B_{m+1}(0)]/(m+1)! = -(-1)^m/m!, and c(-R pi_*S) = exp(Sum (-1)^m (m-1)! ch_m) turns this into exp(-Sum kappa_m/m). CORRECT.
- s = a = 0: -R pi_* O = E^vee - O, c = Lambda^vee(1). CORRECT.
- R^0 = 0 (omega_log^{-1} negative on every component), so R^1 is locally free of rank 3g-3+n. CORRECT.
- The identification of R^1 pi_* omega_log^{-1} with T(-log) as a BUNDLE is stated, not proved (it is standard:
  dual of pi_*(omega (x) omega(D)) = Omega(log)). Only the Chern-class identity is used in the proofs; the bundle
  statement enters only the Gauss-Bonnet sanity check. Spot checks: Mbar_{0,4}: Int c_1 = -1 = chi(M_{0,4});
  Mbar_{1,1}: -1/24-1/24 = -1/12 = chi(M_{1,1}); genus 2, 3 data: 1/120, -1/252 (Harer-Zagier). OK.
- (2.1) and the sigma-decomposition: u2 + u3 = 1/sigma, u2 u3 = -1/sigma, so c_{u2}c_{u3} = prod(1 - alpha(1+alpha)/sigma). CORRECT.
  Genus-2 sigma form re-derived from hvc_which_class_out.txt (w2 = tau): matches 7 - (12j-28)sigma + (5j^2-26j+24)sigma^2.

## 3. Section 3: line-bundle mechanism, Lemma 3, Proposition 3
- 0 -> S(-x_i) -> S -> S|_{x_i} -> 0, x_i^*omega_log = O, x_i^*O(x_i) = L_i^{-1}: c(-R pi_* S_D) = c(T(-log))/prod prod(1-k psi_i). CORRECT
  (consistent with the psi-term of Chiodo's formula: B_{m+1}(a+1)-B_{m+1}(a) = (m+1)a^m).
- rank = 3h-3+r-J; deg < 0 on smooth curves iff J <= 2h-3+r. CORRECT. (Notation: N = 3h-3+r here, N = 2h-3+r in the Norlund formula.)
- Lemma 3: psi-polynomial of degree J - r_+, needs c_k with k >= N-J+r_+ > rank. CORRECT.
- Prop 3: on C_0^~ the bundle is omega^{-1}((j-1)x_1 - a - q - q'); vanishing at the e0 attachment points gives degree
  j+1-2(g0+l0+e0). Self-node preimages are NOT vanishing conditions (only gluing), but the bound is an injection into
  sections on C_0^~, so it is fine. Minor: "stability forces g0+l0+e0 >= 2" is false in general (stability gives
  2g0+e0+2l0 >= 2; g0 = 1, e0 = l0 = 0 is stable); it holds here because e0 = 0 means C = C_0 of genus g0+l0 = h >= 2. CORRECT.

## 4. Section 4: (4.1) and Proposition 4 (Hurwitz, j = 3, h >= 3)
- (4.1): (1-psi)(1-2psi)/(1-j psi) = 1 + (j-3)psi + (j-1)(j-2)psi^2/(1-j psi); the first two terms pair with c_N, c_{N-1}
  of V_2 (rank N-2) and vanish; c_{N-2}(V_2 (x) L^j) = Sum_b c_{N-2-b}(V_2)(j psi)^b. CORRECT, polynomial in j.
- pi_* S_3 = 0 (torsion free, generically 0 for h >= 3), so L^{-3} -> V_2 is a sheaf injection; fibre map = connecting map
  H^0(S_3|x_1) -> H^1(S_2), zero iff H^0(C,S_3) != 0. Z(s)_red = Z_3. Fulton's localized top Chern class lies in
  A_{dim - rank}(Z_3) = A_2(Z_3). CORRECT.
- Dimension bound: V irreducible surface in Z_3, V' a component of xi_Gamma^{-1}(V) over V; xi_Gamma^* psi_1 = p_{v0}^* psi;
  Int_V psi_1^2 = 0 if dim p_{v0}(V') <= 1. The section lives on C_0 only; eta in the stratum
  H(j-1, -2^{e0}, -1^{2l0}, -1^{d0}) (at least one pole since e0 >= 1), dim (projectivised) 2g0-3+n; coincidences of E
  with x_1, a, q or with each other only lower n. dim Y_0 <= j-1-e0 <= 1 for j = 3. Spot cases: (g0,e0,l0) = (0,1,1):
  Y_0 = M_{0,4} (dim 1); (1,1,0): 2-torsion locus (dim 1); (0,2,0): point. CORRECT.
  The citation "meromorphic strata have dimension 2g-2+n" should be made precise (Boissy 2015; Bainbridge-Chen-Gendron-
  Grushevsky-Moeller), but it is used correctly.

## 5. Section 7: components, Lemma 5, Proposition 6, Theorem 7
- e_1 = lambda_1 (ch_2 = 0), e_2 = lambda_2 - 6ch_3 (p_2 = p_4 = 0), e_h = lambda_h c(E), c(T(-log)) e_h = lambda_h K. CORRECT.
- Lemma 5 needs only c(V_2)(1-2psi) = pi^* c(T_h(-log)) as CHERN CLASSES: c(V_2) = Lambda^vee(1) K_{h,1}/((1-psi)(1-2psi))
  = pi^*(Lambda^vee(1)K_h)/(1-2psi) by Corollary 1 and Lemma 2 on Mbar_h. The bundle statement "V_2 = pi^*T_h(-log)/L^{-2}"
  is not proved but not needed. Degree 3h-3 part with c_{3h-3}(V_2) = 0 gives psi e(V_2) = -(1/2)pi^*c_{3h-3}. CORRECT.
- [Theta_2]_{3h-4} = (1/2) pi_* pi^*(top class) = 0 as a class (projection formula). Hence ANY degree <= 1 class pulled back
  from Mbar_h pairs to zero at j = 2 (not only lambda_1). CORRECT.
- Lambda_g component (b): re-derived. With K re-expanded as dilaton leaves, lambda_g formula -> genus-0 multinomial:
  Sum_m (N+m)!/m! [s^{N+m}] e^{Js}(1+s-e^s)^m = N! Res e^{Js}(s - (1+s-e^s))^{-N-1} = N! Res e^{Js}/(e^s-1)^{N+1}
  = (J-1)...(J-N). CORRECT. Genus 2 directly: Xi^{(2)}_2(j) = (j-1)(j-2) Int_{Mbar_2} lambda_2 kappa_1 = b_2 (j-1)(j-2). CORRECT.
- String/dilaton (a) re-derived at c = 1: lam_1 = 1 + psi - T, string derivative of T is T - psi; the m-sums telescope
  to Xi(1,j) = J Xi(j), Xi(0,j) = (J - r - 2h + 2) Xi(j). CORRECT.
- Theorem 7: irreducible points at h = 2 are only j = (2); reduction by string/dilaton is a correct induction
  (J' = 0 or D' = 2h-2 give the vanishing factor at the boundary). Xi_2(2) = 0 - 0/sigma + 0/sigma^2. CORRECT.

## 6. Does the argument "prove" something false? (WHICH CLASSES test)
Every proved step uses a property that fails for the counterexample classes:
- Lemma 2 needs dilaton scale = the weight carried by Lambda^vee: fails for KW (no Lambda^vee to absorb K) and for the
  CY triple with leaf scale 2 (K at scale 2, Lambda^vee at 1). Data: KW genus 2 at j = 2 gives -2/5760 != 0.
- The k = h cancellation uses e_h = lambda_h c(E), i.e. u2 + u3 = -u2 u3, which is exactly w1+w2+w3 = 0 (with u1 = 1).
- Quantitative test of the component machinery on NON-CY classes (scratchpad ref1_classes.pl, exact rationals):
  * genus 2, triple (1, w2, w3) generic, j = 2: components k = 0, 1 vanish (P3, P6), the surviving degree-2/3 parts give
    Xi(2) = -2 (u2+u3)(u2+u3+u2 u3)/5760 = -2(w2+w3)(w2+w3+1)/(5760 w2^2 w3^2); the stored numerator agrees at 9/9 points.
    It vanishes iff w2 + w3 = 0 (Mumford-trivial pair, back to the Hurwitz class) or w2 + w3 = -1 (CY): exactly the
    observed classes.
  * genus 2, Lambda(1)Lambda(w2), j = 2: prediction u2^2 (j-1) Int lambda_2 [Theta_2]_1 = -2 u2^2/5760; data -2/(5760 w2^2). Agrees.
  * genus 3, Lambda(1)Lambda(w2): u2^0 and u2^1 coefficients vanish at j = 2 (P3, P6) and at j = 3 (P4 and, empirically,
    the open case (3, k=1)); the u2^2, u2^3 coefficients are 288, 212 (j = 2) and 128, 64 (j = 3), non-zero.
So the proofs pass the test: they do not prove vanishing for any class where it fails.

## 7. Closed formula (6.1): evidence check (it is NOT claimed proved)
- Transcription of the data into sigma-form (check_closed_formula_r1.pl, psihat): re-derived by hand from the stored output.
  genus 2: from hvc_which_class_out.txt (w2 = tau), all w2-coefficients match. genus 3: from the CY control line, all
  seven w2-coefficients match. genus 4: from hvc_multiprobe_g4_out.txt with Xi_norm = c^{11} Xi, sigma = -x/(1+x)^2,
  num = Sum_k a_k (-x)^k (1+x)^{8-2k}: all 21 coefficients match. genus 1: Xi_1(j) = (j-2-1/sigma)/24 re-derived on Mbar_{1,1}.
- Sample size: the bracket has degree <= h in j and degree <= h in sigma (E^m is divisible by s^m, E/s has
  l^{2k}-coefficient of degree k-1 in s; formula is symmetric under e1 <-> e2), and the 5 tau values give 5 distinct
  sigma; so 5 x 9 points per genus prove the polynomial identity for h <= 4. Re-ran the script: 180 points, 0 mismatches.
- Re-checked: Hurwitz limit S^{j+2} e^{-jxS'/S} = S^2 H^j, S^2 H^j = H^j - (x/j)(H^j)', 1 + x^2(log S)'' = S^2,
  genus-3 Hurwitz bracket (80 - 126 j + 35 j^2)/2903040 reproduced from [x^6]H^j; 9a items 1-4 (lambda_g, Harer-Zagier,
  <tau_{3h-2}>, and the j^{3h-3} coefficient via <lambda_1 tau_{3g-3}> = g(g+4)/(5 24^g g!), re-derived from
  12 lambda_1 = kappa_1 - psi + delta) all correct. The residue form 6a(i) reduces to (6.1) by u = v - beta.
- "Read off from genus 1-3": E_2, E_3 were fitted; E_4 and all of genus 4 (about 14 independent numbers) are predictions.

## 8. Logical structure
Dependencies: Lemma 1 -> Cor 1 -> Lemma 2 -> (2.1) -> {Lemma 3 -> P3 -> (4.1) -> P4}; {7(a) -> Lemma 5 -> P6};
(b) lambda_g; (a) string/dilaton -> reduction -> T7. The conjectural (6.1) is never used in a "proved" item. No circularity.
Status section 9 matches the body.

## 9. Issues found (none fatal)
1. GAP 5 IS NOT A GAP. Nonvanishing outside the range follows for ALL h, r from the proved item (b): sigma^h Xi is a
   polynomial in sigma whose constant term is (-1)^h Xi^{(h)} = (-1)^h b_h (J-1)...(J-N), N = 2h-3+r, which is non-zero for
   every integer J outside [1, N]. So the "=>" direction of Conjecture M (as an identity in tau) is proved; (6.1) is not needed
   for it (the remark in section 6 about the tau -> 0 value of the bracket is this same fact).
   Also: with the reduction lemma of T7, the genus-3 irreducible points are r <= 3, J <= 3+r, and all of them are in the exact
   lattice tests of hvc_multiprobe_out.txt; so Conjecture M at genus 3 holds for all r (computer-assisted). Genus 4 would
   additionally need the irreducible points with r = 3, 4, 5.
2. Section 6 "Consequence": for the Theorem-C configuration with at least one phi leaf (r >= 2),
   Xi(0, 1^{r-1}) = 0 by the string factor j(j+1)...(j+r-2); "vanishes exactly when j in {1..2h-2}" is wrong at j = 0
   (the set is {0..2h-2}, in agreement with Conjecture M). Harmless for scales c/j, j >= 1.
3. Prop 3: "stability forces g0+l0+e0 >= 2" needs h >= 2 in the case e0 = 0; stability alone gives 2g0+e0+2l0 >= 2.
4. Bundle identifications (R^1 pi_* omega_log^{-1} = T(-log); V_2 = pi^*T_h(-log)/L^{-2}) are asserted, not proved; only
   Chern-class identities are used and those do follow. Rephrase or cite.
5. P4: cite the dimension of strata of meromorphic differentials precisely (Boissy; BCGGM): 2g-2+n with at least one pole.
6. Wording: "Hurwitz component at j = 2, 3 in every genus" holds for j = 3 only when h >= 3; "every tau" means tau not in
   {0, -1} (sigma = 0, where the class has poles).
7. Notation: N = 3h-3+r in sections 3-4 vs N = 2h-3+r in the Norlund formula; section 5's tag "[data h=2,3; identity]" for the
   gap formula, which is (6.1) at j = 0 and so conjectural for h >= 5.

## 10. Verdict
NOT REFUTED. The status "partial" is accurate. Proved items L1, L2 (as a Chern-class identity), L3, P3, P4 (h >= 3),
Lemma 5, P6, (a), (b) and T7 (genus 2, all r, all tau) are correct. Formula (6.1) is a conjecture backed by exact data for
h <= 4: the transcription of that data was re-checked and the sample is large enough to prove the identity for h <= 4.
The admitted gaps 1-4 are real. Gap 5 can be removed (issue 1).
