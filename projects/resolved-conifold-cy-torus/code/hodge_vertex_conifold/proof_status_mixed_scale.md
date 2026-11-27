# Proof status: mixed-scale vanishing (Conjecture M) and the all-genus regularity theorem

> **Addendum, 2026-10-31 (main session, after this synthesis).**
> (1) The seven missing genus-4 points {2,2,2}, {2,2,3}, {2,2,4}, {2,3,3}, {2,2,2,2}, {2,2,2,3}, {2,2,2,2,2} all vanish identically in x
>     (hvc_multiprobe_g4b.wls; the six control points at D = 6 are non-zero). So M_irr(4) holds for every tau, and Conjecture M is
>     established at genus 2, 3, 4; with Lemma C (verified through t^5) the leaf route gives the theorem through genus 4.
> (2) Formula (6.1) agrees with the exact Hodge functional identically in j and x at h = 2, 3, 4 (hvc_check61.wls).
> (3) P2 checked exactly at tau = 2 and tau = 3 (no global involution) through genus 3, r <= 3, with a general-tau exact TR code
>     (hvc_tr_tau.wls); Lemma D holds there in both directions.
> (4) The shear/escaping-branch-point heuristic for Lemma D predicts a tau-independent first moment Res y^(2g-1) omega_{g,1}; the exact
>     values depend on tau (hvc_moment_test.wls), so that route needs a correction at the logarithmic point v = oo.
> Current summary: conifold_regularity_note.md section 7a.
>
> **Addendum 2, 2026-10-31 (later): Lemma C proved to all orders, without topological recursion.** leaf_proof/leaf_structure_proof.md
> proves F1, F3 and Lemma C (L1, L2 for j = 1, 2 and every b, L3) for every k >= 0 and at every psi-order, together with E_21(b, a) = E_12(a, b),
> from the exact hypergeometric form of the singular row, rho r_{0,1} = (c/(c+z)) (tau/t) 2F1(1/z+1, 1/z+1; c/z+2; tau) Exp[-dphi/z], and an exact
> t-flow for E_11 (summary: conifold_regularity_note.md section 7b). With Theorem R (section 3) the all-genus regularity theorem is therefore
> conditional only on Lemma A (M_irr(h) for h >= 5); it holds through genus 4, genus 2 without computer input. Lemmas B and D are not needed.
> Two corrections to this file. (a) Theorem R is not the Reduction Lemma of conifold_regularity_note.md s.5a "in leaf language": it replaces
> it. The Reduction Lemma factorises an E_11 block between two singular vertices of genus >= 2 into columns, which is impossible with
> coefficients regular at t = 0 (the block has rank 1 at t = 0 and rank 2 for t != 0; leaf_structure_proof.md Remark 13.1), while Theorem R
> expands rho^2 E_11 in t-independent leaves. (b) Theorem R uses Lemma C also at k = 0 and the clause E_21(b, a) = E_12(a, b) (orientation of
> mixed edges in graphValueC); both are part of what is now proved.
>
> **Addendum 3, 2026-11-02 (round 6; supersedes the status statements below where they differ).** Lemma A = M_irr(h), for every h >= 2 and non-empty R, is now
> lemmaA/LEMMA_A_PROOF.md, Theorem L (s.7.4 there for the conifold convention), P modulo the Marino-Vafa formula (MV) and textbook facts, with no TR and no computer step; two round-6 referees, "sound with minor
> fixes", no major item. So Corollary R (leaf_proof/leaf_structure_proof.md s.13.2), a statement about the decoration sum, is proved modulo (MV) and textbook facts. Its identification with F_g is Proposition 13.2 there [repaired in round 6d; passes on the repaired paragraph: hodge_vertex/lp2_explore/ROUND6_STATUS_2026-11-02.md], P modulo the literature inputs listed in its heading and in (G3); genus 3 and 4 are no longer
> computer-assisted (their exact data are independent confirmations). Graphs with an edge need Lemma A only for h <= g - 1; the single-vertex graph also needs h = g. The value of the single pole is Theorem VI of the local-P2 draft (draft s.7.8;
> further inputs, beyond (MV) at mu = 1^a and textbook facts: lemmaA/LEMMA_A_PROOF.md Prop. 2.1, Lemma 3.2 and eq. (4.2) at k = 2, with the Frobenius formula; LLZ eq. (5) and the vertex-operator form of psi, psi^* are no longer inputs; two referee passes (round 5; round 6, D_VI); specialisation in leaf_structure_proof.md s.13.3 (5)), and [BS] Prop. 3.1 is now an independent check. Lemmas B and D, and "mod P2" below, are no longer needed for anything.


Synthesis of the proof routes of 2026-10-30/31 (proof_M_hurwitz.md, proof_M_residue.md, leaf_structure.md) after the referee
reports (referee_hurwitz_1.md, referee_residue_1.md). Written 2026-10-31. No wolframscript was run for this file; the small
checks quoted as [this file] were done in Perl/Math::BigRat (scripts in review_checks/status_*.pl).

User's question: does a single geometric leaf of a different scale vanish in the Theorem-C configuration?
The answer is in section 0.

Conventions as in mixed_scale_findings.md. Normalise c = 1 (Xi^(c) = c^-(3h-3+2r) Xi^(c=1), referee_hurwitz_1 s.1), so
w = (1, tau, -1-tau), sigma := tau(1+tau) = -w2 w3, lam_j = 1/(1 - j psi), T = -psi^2/(1 - psi), r = number of leaves,
J = Sum j_i, D = J - r, N := 2h-3+r. "Vanishes" means identically in tau (equivalently in x), unless a value of tau is named.
P2 = the spectral-curve dictionary (section 1.6); "mod P2" means: proved, assuming P2.

## 0. Summary, and the answer to the user's question

Theorem-C configuration with one leaf of another scale: Y_h(j; n) := Sum_m 1/m! < Omega_w lam_j phi^(n-1) T^m >_h, phi = lam_1.
The string identity (1.2) gives Y_h(j; n) = j(j+1)...(j+n-2) c^(-2(n-1)) Xi_h(j). So for j >= 1 the question does not depend
on n and is the one-leaf question Xi_h(j) = 0; for j = 0 and n >= 2 it vanishes trivially (factor j = 0).

Answer, for every genus h >= 2:
 * j = 1: vanishes (Theorem C).
 * j = 0, n = 1: Xi_h(0) = ((2-2h)/c^2) Z[h,0] with (-x c^3)^(h-1) Z[h,0] = gap_h != 0 (F4). j = 0, n >= 2: vanishes.
 * 2 <= j <= 2h-2: vanishes. PROVED: h = 2 (all tau); h = 3, j = 2, 3 (all tau, mod P2; new, section 1.7) and j = 2, 3, 4
   (computer-assisted, exact); every h at tau in {1, -2} (mod P2) and at tau = -1/2 for odd j (mod P2). Exact data: h = 4.
   The conjectural closed formula (6.1) (proof_M_residue.md) predicts it in every genus; it is now confirmed against exact
   topological recursion at h = 5 (tau = 1, tau = -1/2) and h = 6 (tau = -1/2) as polynomial identities in j (section 1.8).
 * j >= 2h-1: does NOT vanish (proved, every h); for each such j, Xi_h(j; tau) = 0 for at most 2h values of tau.
 * non-integer j (a generic scale c/j): does NOT vanish (proved, every h).
So a single geometric leaf of scale c/j vanishes in the Theorem-C configuration exactly for j in {1, ..., 2h-2}
(and trivially for j = 0 when n >= 2). These are exactly the scales of the poles z = -c/j, j <= k+1 <= 2h-2, that occur in the
polar orders t^k, k <= 2h-3, of a genus-h vertex (Fact G); this is why Conjecture M is the right input for regularity.

Status of Conjecture M (vanishing direction; non-vanishing is proved in general, 1.3):
   genus 2, all r, all tau ............ PROVED (T7, proof_M_residue.md s.7(d))
   genus 3, all r, all tau ............ PROVED, computer-assisted (exact lattice data cover every irreducible point)
   genus 4, r <= 2 .................... exact data; r = 3, 4, 5 irreducible points not computed
   all genera, all r, tau = 1 or -2 ... PROVED mod P2 (Theorem 3.1, proof_M_hurwitz.md s.3)
   all genera, odd parts, tau = -1/2 .. PROVED mod P2
   all genera, general tau ............ OPEN; reduced (1.7) to the degree bound for C_tau at h integer values of tau. [Superseded by Addendum 3: proved modulo (MV), all h, all r, all tau with sigma != 0.]
All-genus regularity (section 3): F1 + F3 + leaf structure (Lemma C = Lemma L of leaf_structure.md) + Conjecture M (irreducible points) => all-genus
theorem, by a proof with no further gap; no obstruction was found. F1, F3 and Lemma C are verified exactly (so proved, as
finite computations) only for the coefficients that occur up to genus 4; M is proved only as above. [Superseded by Addendum 2:
F1, F3 and Lemma C are proved to all orders.]

## 1. What is proved (after the referees' objections)

### 1.1 Reformulations [proved; referees re-derived all of them]
(a) Kappa form (Lemma 1 of proof_M_residue.md). Xi_h(j_1..j_r) = Int_{Mbar_{h,r}} Omega_w K prod_i 1/(1 - j_i psi_i),
K = exp(-Sum_k kappa_k/k). Proof: T has psi^(>=2), so on Mbar_{h,r+m} the psi_i of the probes may be replaced by pull-backs
(the correction divisors put a forgotten point on a genus-0 bubble where its psi^(>=2) is 0); Omega_w is a pull-back; and
Sum_m 1/m! pi_*(prod T(psi_{r+k})) = exp(Sum_k kappa_k [z^k](-log(1 - T(z)/z))) = exp(Sum_k kappa_k [z^k] log(1-z)).
(b) Corollary (r = 1). kappa_k = pi^*kappa_k + psi^k on Mbar_{h,1}, so K = pi^*K_h (1 - psi) and
   Xi_h(j) = (j-1) Sum_{b>=0} j^b Int_{Mbar_h} Omega_w K_h kappa_b       (kappa_0 = 2h-2).
(c) Log tangent bundle (Lemma 2, as a Chern-class identity). Lambda^vee(1) K = c(-R pi_* omega_log^{-1}) =: c(T(-log)).
Proof: Mumford/Chiodo GRR for omega_log^s(-Sum a_i x_i); passing from s = 0 to s = -1 changes only the kappa terms, by
Sum_m (-1)^m (m-1)! [B_{m+1}(-1) - B_{m+1}(0)]/(m+1)! kappa_m = -Sum_m kappa_m/m. (Only the Chern-class identity is used below;
the bundle identification is standard but not needed.) Here the leaf/dilaton scale c = w1 is what lets Lambda^vee(w1) absorb K.
(d) Twisting. For S_D = omega_log^{-1}(Sum j_i x_i): c(T(-log)) prod 1/(1 - j_i psi_i) = c(-R pi_* S_D) prod_i p_{j_i}(psi_i),
p_j(psi) = prod_{k=1}^{j-1}(1 - k psi) (p_0 = 1); rank(-R pi_* S_D) = 3h-3+r-J and deg S_D = J - (2h-2+r) on smooth curves,
so the range of Conjecture M is exactly "deg S_D < 0" (Lemma 3 of proof_M_residue.md).
(e) Components (s-decomposition). With alpha_i the Chern roots of E, u2 + u3 = 1/sigma, u2 u3 = -1/sigma:
   c_{u2}(E^vee) c_{u3}(E^vee) = prod_i (1 - alpha_i u2)(1 - alpha_i u3) = prod_i (1 - (alpha_i + alpha_i^2)/sigma), hence
   Xi_h(R; tau) = Sum_{k=0}^{h} (-1/sigma)^k Xi^(k)_h(R),  Xi^(k)_h(R) := Int c(T(-log)) e_k(alpha + alpha^2) prod lam_{j_i},
with tau-INDEPENDENT Xi^(k). k = 0 is the single class Lambda^vee(1) (tau -> oo, "Hurwitz component"); e_1 = lambda_1
(ch_2(E) = 0); e_h = lambda_h c(E), and c(T(-log)) e_h = lambda_h c(E) c(E^vee) K = lambda_h K (Mumford).

### 1.2 String and dilaton relations; reduction to irreducible points [proved]
For h >= 2 and any multiset R:  Xi(1, R) = J_R Xi(R),   Xi(0, R) = (D_R - 2h + 2) Xi(R)   (c = 1).
Proof. lam_1 = 1 + psi lam_1. Insert the constant 1 by the string equation (valid: Omega_w is pulled back, h >= 2): the probes
give S lam_j = (lam_j - lam_j(0))/psi = j lam_j, the dilaton leaves give T/psi = -psi lam_1; summed over m the latter is
-<psi lam_1 ...>, which cancels the second term psi lam_1 of lam_1: Xi(1, R) = J_R Xi(R). For lam_0 = 1 the string term is
J_R Xi(R) - B with B = Sum_m 1/m! <lam_R (psi lam_1) T^m>; write psi lam_1 = psi - T, apply the dilaton equation to psi and
reindex the T-sum: B = (2h-2+r) Xi(R), so Xi(0, R) = (J_R - r - 2h + 2) Xi(R).  QED
Iterating: Xi(1^p, R) = prod_{i=0}^{p-1} (J_R + i) Xi(R). In particular Xi(1^p) = 0 for p >= 1 (Theorem C) and
Y_h(j; n) = j(j+1)...(j+n-2) Xi_h(j) (section 0).
Reduction. If a point of the range contains j_i = 1, removing it keeps the rest in the range or produces the factor J' = 0;
if it contains j_i = 0, removing it keeps the rest in the range or produces the factor D' - 2h + 2 = 0. Hence the vanishing
direction of Conjecture M at genus h is equivalent to
   M_irr(h):  Xi_h(R) = 0 for every non-empty multiset R of integers >= 2 with D(R) <= 2h-3   (so 1 <= r <= 2h-3, 2r <= J <= 2h-3+r).
The empty multiset is excluded: Xi_h(empty) = Z[h,0] is the gap.
This is exactly the input of the leaf-structure proposition (section 3).

### 1.3 The lambda_g component; the non-vanishing direction [proved, all h, all r]
Xi^(h)_h(R) = b_h (J-1)(J-2)...(J-N),  N = 2h-3+r,  b_h = [z^(2h)] (z/2)/sin(z/2)   (any real j_i).
Proof. By 1.1(e), Xi^(h) = Int lambda_h K prod lam_{j_i}; re-expand K as dilaton leaves (1.1(a) backwards; lambda_h is a
pull-back) and use <lambda_h prod tau_{a_i}>_h = b_h <tau_0^(2h) prod tau_{a_i}>_0 with <prod_{i=1}^n tau_{a_i}>_0 = (n-3)!/prod a_i!.
Summing the probes (j^a) and m dilaton leaves (-Sum_{d>=2} psi^d) gives Sum_m (N+m)!/m! [s^(N+m)] e^(Js)(1 + s - e^s)^m
= N! [s^N] e^(Js) (1 - (1+s-e^s)/s)^(-N-1) = N! Res_s e^(Js) ds/(e^s - 1)^(N+1); with w = e^s - 1 this is
N! Res_w (1+w)^(J-1) dw/w^(N+1) = N! C(J-1, N).  QED  (Checked: x -> 0 data give b_2, b_3, b_4.)
Consequences. (i) sigma^h Xi_h(R) = Sum_k (-1)^k sigma^(h-k) Xi^(k)(R) is a polynomial in sigma of degree <= h with constant
term (-1)^h b_h (J-1)...(J-N). Outside the range this is non-zero, so Xi_h(R) is not identically zero: the NON-VANISHING
direction of Conjecture M holds for all h, r (referee_residue_1, issue 1), also for non-integer j_i (a generic leaf scale never
vanishes), and outside the range Xi_h(R; tau) = 0 for at most h values of sigma (2h values of tau).
(ii) Inside the range the constant term vanishes: sigma^h Xi_h(R) = sigma P_R(sigma) with deg P_R <= h-1.

### 1.4 Genus 2, all r, all tau [proved; T7 of proof_M_residue.md]
Self-contained proof of the two components it needs, valid in every genus h >= 2:
Lemma (V_2). On every stable curve of genus h >= 2, S_2 = omega_log^{-1}(2 x_1) has no non-zero section. Hence
V_2 := -R pi_* S_2 is a vector bundle of rank 3h-4 on Mbar_{h,1}, with c(V_2) = pi^* c(T_h(-log))/(1 - 2psi).
Proof. Sections vanish on the components without x_1 (negative degree), hence at the e0 nodes attaching the component C_0 of
x_1; on its normalisation (genus g0, l0 self-nodes) the twisted bundle has degree 3 - 2(g0 + l0 + e0) < 0, because
g0 + l0 + e0 >= 2 (if e0 = 0 then g0 + l0 = h; otherwise stability gives 2g0 + e0 + 2l0 >= 2). The Chern class follows from
1.1(d) (D = 2x_1) and 1.1(b),(c): c(T_{h,1}(-log)) = pi^* c(T_h(-log)) (1 - psi).  QED
Then Theta_2 := c(T_h(-log)) Sum_b 2^b kappa_b = pi_*(c(V_2)(1 - psi)) satisfies [Theta_2]_{3h-3} = pi_*(c_{3h-2}(V_2) - psi c_{3h-3}(V_2)) = 0
and [Theta_2]_{3h-4} = -pi_*(psi e(V_2)) = (1/2) pi_* pi^* c_{3h-3}(T_h(-log)) = 0 (degree 3h-3 part of c(V_2)(1 - 2psi) = pi^*c(T_h),
with c_{3h-3}(V_2) = 0). Since Xi^(k)_h(2) = Int_{Mbar_h} e_k(alpha+alpha^2) Theta_2 (1.1(b),(c)) and e_1 = lambda_1:
   Xi^(0)_h(2) = 0 (P3)  and  Xi^(1)_h(2) = 0 (P6)  for every h >= 2.
Genus 2: Xi_2(2) = Xi^(0)_2(2) - Xi^(1)_2(2)/sigma + Xi^(2)_2(2)/sigma^2 = 0 - 0 + b_2 (2-1)(2-2)/sigma^2 = 0, and R = {2} is the only
irreducible point at h = 2 (1.2: r <= 2h-3 = 1).  QED

### 1.5 Hurwitz component at j = 3 [proved for h >= 3; P4 of proof_M_residue.md, s.4]
Xi^(0)_h(3) = 0 for h >= 3: L^{-3} -> V_2 vanishes exactly on the jump locus Z_3 of S_3, Xi^(0)(3) = 2 Int_{[Z_3]^loc} psi_1^2
(Fulton's localized top Chern class), and psi_1^2 kills A_2(Z_3) because the strata of meromorphic differentials met there have
dimension <= j - 1 - e0 <= 1 (dimension 2g-2+n of meromorphic strata with a pole: Boissy; Bainbridge-Chen-Gendron-Grushevsky-
Moeller -- the precise citation still has to be added). Re-derived by the referee. Fails for h = 2 (Weierstrass divisor).

### 1.6 tau in {1, -2}: all genera, all r; tau = -1/2: odd parts [Theorem 3.1 of proof_M_hurwitz.md]
(a) TR theorem (unconditional). On C_1: x = v/(1+v)^2, y = -v, B = dv dv'/(v-v')^2, TR with respect to X = log x (v = -1 is a
double pole of x and a simple pole of dX, not a branch point), the correlators omega_{g,n} (g >= 2, n >= 1) have total x-degree
>= 2g-2+n at the zeros v = 0, oo of x. Proof: sigma(v) = 1/v is a global involution with X o sigma = X; the linear loop
equation makes omega_{g,n} odd in each variable (this needs only TR, not the Hodge dictionary); in the expansion of the
kernel at v_0 = oo, [t_0^(m+1)] omega_{g,n+1} = (1/2) Res_{v=1} h_m(v, 1/v) Brk/dX; the residue theorem moves this residue to
v = 0, oo, v_j (no pole at v = -1 for g >= 2), where the induction hypothesis gives degree >= 2g-2+n-m. The tau = -1/2 curve
(sigma(v) = -1-v) works the same way but only sees odd parts. Referee: every step correct; exact recursion confirms the
bound through genus 5 (and this file: genus 6, n = 1, both curves).
(b) Dictionary P2 (needed to translate (a) into Xi). For every tau: Xi_h(mu; tau) = sigma^(h-1) [x^mu] F^{C_tau}_{h,r} / prod_i
C(mu_i(tau+1) - 1, mu_i - 1), with C_tau: x = w/(1+w)^(tau+1), y = -1/w, expanded at w = 0 (normalisation of [x^mu]F as in
proof_M_hurwitz.md s.2 and review_checks/referee1_trC1.pl). Status: a theorem modulo BKMP for C^3
(L. Chen; J. Zhou) and the Eynard-Orantin variational formula, via the route of referee_hurwitz_1.md s.4 (the regular combined
deformation Y_B + e Y_0, rationality in e, homogeneity, which also explains the n-independent factor sigma^(1-g)); this route
is sketched, not yet written out. Exact agreement with the Hodge formulas: tau = 1 at (g,r) = (2,1),(2,2),(3,1),(3,2),(4,1);
tau = -1/2 (odd parts) at (2,2),(3,2),(4,1).
(c) Consequence (mod P2): Conjecture M, vanishing direction, at tau = 1 and tau = -2 for every h and r; at tau = -1/2 for every
h and every R with all parts odd. The argument does not over-reach: it gives no zero at j = 2h (absent for tau = 1), and it
uses CY + leaf scale = weight + a global involution, which all fail for the counterexample classes.

### 1.7 New consequences of 1.3-1.6 [this file]
(a) Lemma (h values of sigma suffice). Let R be in the range at genus h. If Xi_h(R; tau) = 0 at h values of tau with distinct
non-zero sigma -- or, more generally, if h linearly independent conditions on (Xi^(0)(R), ..., Xi^(h-1)(R)) are known -- then
Xi_h(R) = 0 identically. Proof: 1.3(ii); the h unknowns Xi^(0..h-1) are tau-independent.
(b) Proposition (genus 3, one leaf). Xi_3(2) = Xi_3(3) = 0 for all tau, mod P2 (no computer data needed).
Proof. j = 2: Xi^(0)(2) = Xi^(1)(2) = 0 (1.4), Xi^(3)(2) = b_3 (1)(0)(-1)(-2) = 0 (1.3); Theorem 3.1 at tau = 1 (sigma = 2) gives
Xi^(2)(2)/4 = 0. j = 3: Xi^(0)(3) = 0 (1.5), Xi^(3)(3) = 0 (1.3); tau = 1 gives -Xi^(1)(3)/2 + Xi^(2)(3)/4 = 0 and tau = -1/2
(sigma = -1/4, j odd) gives 4 Xi^(1)(3) + 16 Xi^(2)(3) = 0; the determinant is -9, so both vanish.  QED
In genus h the same accounting leaves h-2 unknown components at j = 2 against one TR relation and h-1 at j = 3 against two, so
it closes exactly at h = 3 and goes no further. (Genus 3 is in any case covered by the exact lattice data, 1.9.)
(c) Corollary (what is needed for all tau). Mod P2, the vanishing direction of Conjecture M for all h, r and tau follows from
  Lemma D: for every g >= 2, n >= 1 and every integer tau >= 1 (at genus g, any g of them suffice), the TR correlators of
  C_tau: x = v^tau/(1+v)^(tau+1), y = -v, have total x-degree >= 2g-2+n at the simple zero v = oo of x.
For integer tau >= 1 the sigma = tau(tau+1) are distinct and non-zero and the probe normalisations C(mu(tau+1)-1, mu-1) are
positive, so Lemma D at tau is equivalent (mod P2) to Conjecture M at that tau; then apply (a). Lemma D holds for tau = 1
(Theorem 3.1) and is consistent with the float recursion at tau = 2, genus 2 (proof_M_hurwitz.md s.2-3: zeros j = 1, 2 only);
at tau = 2, 3 the correlators have degree 1 at the OTHER zero v = 0, which is why the involution argument cannot be copied.

### 1.8 Exact TR test of the closed one-leaf formula (6.1) at genus 5 and 6 [this file]
(6.1): sigma^h Xi_h(j) = prod_{l=1}^{2h-2} (j-l) [lambda^(2h)] S(e1 lambda)^(1+j e2) S(e2 lambda)^(1+j e1), e1 = -tau, e2 = 1+tau.
It was read off from genus <= 3 Hodge data and proved as a polynomial identity for h <= 4 (referee_residue_1 s.7). Test against
the referee's exact recursions (copies review_checks/status_trC1.pl, status_trhalf.pl with an extended j-range; comparison
review_checks/status_cmp61.pl, TR = sigma^(1-h) Xi_h):
   tau = 1  (sigma = 2):    h = 2..6, j = 1..3h-1:   0 mismatches (genus 6: zeros j = 1..10, first non-zero -297675/33554432 at j = 11)
   tau = -1/2 (sigma = -1/4): h = 2..6, odd j = 1..33: 0 mismatches (17 values = deg_j + 1 at h = 6)
(run times: genus 6 at tau = 1 2.5 min, at tau = -1/2 12 s; genus 5 20 s / 4 s.) The number of j-values is at least
deg_j + 1 = 3h-1 for both sides, so these are identities in j. Hence, mod P2, (6.1) holds for h <= 6 at sigma in {0, 2, -1/4}
(sigma = 0: the lambda_g limit, proved for all h). A proof for a given h needs h+1 values of sigma, so (6.1) remains a conjecture
for h >= 5, now with parameter-free exact confirmations in genus 5 and 6 at two values of tau.

### 1.9 Other proved or computer-proved items (pointers)
- Genus 3, all r, all tau (computer-assisted): the irreducible points at h = 3 are {2}, {3}, {4}, {2,2}, {2,3}, {2,2,2}, all in
  the exact lattice tests of hvc_multiprobe_out.txt (Mathematica, exact in x); with 1.2 this is Conjecture M at genus 3
  (referee_residue_1, issue 1). The non-vanishing direction is 1.3.
- Genus 4: exact data for r <= 2 (hvc_multiprobe_g4_out.txt). Missing irreducible points: {2,2,2}, {2,2,3}, {2,2,4}, {2,3,3},
  {2,2,2,2}, {2,2,2,3}, {2,2,2,2,2}.
- Single class: Xi^(0)_g(2) = 0 for all g >= 1 by cut-and-join on p = X delta_1 (proof_M_hurwitz.md s.5; formal analytic
  continuation, harmless per referee); for g = 1 this is the extra zero j = 2h. The extra zero j = 2h for h >= 2 is explained only
  through (6.1) (identity S^2 H^j = H^j - (x/j)(H^j)', proof_M_residue.md 6a(iii)).
- Consistency identities of (6.1) in all genera (proof_M_residue.md 9a): lambda_g limit, Harer-Zagier at j = 0, the top two
  j-coefficients (DVV; Mumford's lambda_1 formula).

## 2. What remains open, as precise lemmas

[Status 2026-11-02: Lemma A is proved modulo (MV), lemmaA/LEMMA_A_PROOF.md Theorem L; Addendum 3. The route lists below are history.]
Lemma A (= M_irr(h), the vanishing direction of Conjecture M after the string/dilaton reduction). For h >= 4 and every tau:
Xi_h(R) = 0 for every non-empty multiset R of integers >= 2 with D(R) = Sum (j_i - 1) <= 2h-3.
 Equivalent (1.1(e), 1.3): Xi^(k)_h(R) = 0 for 0 <= k <= h-1 (k = h is proved). Equivalent mod P2 (1.7(c)): Lemma D at h integer
 values of tau. Known pieces in every genus: (j = 2; k = 0, 1), (j = 3; k = 0) for h >= 3, all of tau = 1, -2 and odd R at tau = -1/2 (these two mod P2).
 Sub-lemmas, each a self-contained target:
 A1 (one leaf, closed form). Formula (6.1) for all h. It gives Lemma A for r = 1 in both directions, but r = 1 is NOT enough
    for regularity: genus h needs r up to 2h-3 (e.g. {2,2,2} at genus 3, {2,2,2,2,2} at genus 4). No closed formula for r >= 2
    is known (no formula in J alone can satisfy both string and dilaton; the "flow the deviation" guess fails at second order).
 A2 (Hurwitz component). Xi^(0)_h(j) = 0 for 4 <= j <= 2h-2, and Xi^(0)_h(R) = 0 for irreducible R with r >= 2 (jump loci of
    omega_log^{-1}(Sum j_i x_i) on the boundary; for r >= 2 already a divisor of rational tails). Excess-intersection problem.
 A3 (intermediate components). Xi^(k)_h for 1 <= k <= h-1, except (k, j) = (1, 2).
Lemma B (= P2, the dictionary). For all tau, Xi_h(mu; tau) = sigma^(h-1) [x^mu] F^{C_tau}_{h,r} / prod C(mu_i(tau+1)-1, mu_i-1).
 Proof route (referee_hurwitz_1.md s.4): BKMP for the framed one-leg vertex (Chen; Zhou) at the background q_B; the EO
 variational formula along the REGULAR family Y_e = Y_B + e Y_0 (the separate lam_0, lam_1 deformations are irregular and
 must not be used); rationality of both sides in e; homogeneity in y plus the dilaton equation as e -> oo. To be written out,
 including the sign/normalisation conventions of [x^mu] and the n-independent factor sigma^(1-g).
Lemma C (= leaf structure L1-L3 to all orders, with F1 and F3). [PROVED 2026-10-31, leaf_proof/leaf_structure_proof.md Theorem LC, in the
 stronger form that Theorem R needs: also for k = 0 ([t^0](rho q_t) = q_0, [t^0](rho f_{j,b}) in span(lam_1), [t^0](rho^2 E_11) in
 span(lam_1 lam_1')), and with E_21(b, a) = E_12(a, b). See Addendum 2.] With rho(t) = (Delta_1 t^2/(-x c^3))^(-1/2), for every k >= 1:
 (L1) [t^k] rho q_t in span(lam_1, ..., lam_{k+1});   (L2) [t^k] rho f_{j,b} in span(lam_1, ..., lam_{k+1}) for j = 1, 2 and all b;
 (L3) [t^k] rho^2 E_11(psi, psi') in span{lam_i(psi) lam_l(psi') : i, l >= 1, (i-1) + (l-1) <= k};
 plus F1 (Delta_1 = -x c^3 t^-2 (1 + O(t)), Delta_2 regular and non-zero) and F3 (all frame-B edge and dilaton factors regular
 at t = 0). Large-z form: Fact G (only simple poles at z = -c/j, j <= k+1) plus vanishing at z = oo after the rho-normalisation
 (z r_{0,1}(z; t) -> c/rho(t); E_{1j}(z, w) -> 0 at fixed power of w). Verified exactly for k <= 5 at every psi-order available
 at z^9, which is every coefficient occurring in graphs of genus <= 4 (leaf_structure.md, G1). Natural proof: the thimble of the
 degenerate critical point is an Euler integral, tau^(a+b+1) B(a+1, b+1) 2F1(b+1, b+1; a+b+2; tau), whose k-th tau-coefficient
 has exactly the poles of Fact G (conifold_regularity_note.md s.6, step 1); what is missing is the passage from this exact
 integral to the formal z-asymptotics that define R', order by order in t, including the rho-normalisation that removes lam_0.
 The lam_0-freeness is essential: a lam_0 component at order t^k would produce Xi(0, R) = (D_R - 2h + 2) Xi(R), i.e. for R empty
 a non-zero multiple of the gap at order t^(2-2h+k).
Lemma D (TR form of Lemma A for integer tau). For every integer tau >= 1, g >= 2, n >= 1: the EO correlators of
 C_tau: x = v^tau/(1+v)^(tau+1), y = -v (TR with respect to X = log x; one simple branch point v = tau) have total x-degree
 >= 2g-2+n at v = oo. Known: tau = 1 (Theorem 3.1). For tau >= 2 there is no global involution, but x is a rational function of
 degree tau+1 and the GLOBAL linear loop equation holds: the trace of omega_{g,n} over the fibres of x is a holomorphic
 differential on P^1 (the only poles, at v = tau, cancel by the local linear loop equation; at v = 0, -1, where x ramifies but
 X is not a TR branch point, omega_{g,n} is holomorphic), hence zero:
     omega_{g,n}(v_oo(x), J) = - Sum_{l=1}^{tau} omega_{g,n}(v_l(x), J)     (v_l(x) -> 0 the tau sheets near v = 0).
 This replaces anti-invariance. What is missing is the analogue of the residue step: control of the quadratic loop equation
 summed over pairs of sheets, including the mixed pairs (v near tau, v' on a far sheet), with the kernel expanded at v_0 = oo.

Also open but not needed for regularity: the extra zero j = 2h of the single class for h >= 2; (6.1) itself beyond h <= 4; a
closed formula for r >= 2; the first non-vanishing layer at fixed tau (proof_M_hurwitz.md s.6).

## 3. Does Conjecture M plus the leaf structure give the all-genus regularity theorem?

Yes; the implication has no further gap. Precisely:

Theorem R (conditional). Fix g >= 2. Assume F1, F3, Lemma C through order t^(2h-3) for every 2 <= h <= g (only at the
psi-orders that occur in genus-g graphs). (a) If M_irr(h) holds for every 2 <= h <= g - 1, every coloured stable graph of genus g with
at least one edge is regular at t = 0 after summing over decorations (no vertex of such a graph has genus g: leaf_proof/leaf_structure_proof.md Lemma 13.0). (b) If M_irr(g) holds, the single-vertex graph has the single pole
(-x c^3)^(g-1) F_g(q_0) t^(2-2g). By Theorem VI of the local-P2 draft (draft s.7.8; specialisation h_w, w = (c, -x, -1) written out in leaf_structure_proof.md s.13.3 (5)), (-x c^3)^(g-1) F_g(q_0) = beta_{2g-2}(x, 1) = gap_g, the t^(2-2g) coefficient of P_g Li_{3-2g}(Q)/(eps1 eps2);
so F4 at n = 0 holds in genus g, with no input from F_g itself. Independent check: the conifold F_g is P_g Li_{3-2g}(Q)/(eps1 eps2) up to Q-independent terms ([BS] Prop. 3.1) and the only singular term of Li_{3-2g}(e^t) is (2g-3)! t^(2-2g), so by (a) and Proposition 13.2 the pole must be gap_g t^(2-2g). Assuming the hypotheses
for all g gives the all-genus theorem.
Proof. (i) Homogeneity. At a singular vertex v of genus h >= 2 with n_v half-edges, the vertex factor is
Delta_1^(h-1) d^(n_v) F_h(q_t)[f_1, ..., f_{n_v}], F_h the genus-h potential of the CohFT Omega_w
(the vertex class is the t-independent Omega_w in frame B, all t-dependence sits in Delta, Dil, E). The dilaton equation, in the form (H0) on the domain q_0 + t C[[t]][[psi]] of leaf_structure_proof.md s.13.0 (needed because rho q_t has psi-coefficient -rho), gives
F_h(lambda q) = lambda^(2-2h) F_h(q), so this equals (-x c^3/t^2)^(h-1) d^(n_v) F_h(rho q_t)[rho f_1, ..., rho f_{n_v}], because
rho^2 Delta_1 = -x c^3/t^2. (ii) Leaf expansion. By L1-L3, rho q_t = q_0 + Sum_{k>=1} t^k delta_k, and every leg, loop and
singular-singular edge is a t-series whose order-t^k coefficient is a combination of leaves lam_j (j >= 1) of total degree
D <= k (an E_11 edge between two singular vertices splits its D <= k_e between its ends; at an end attached to a singular vertex
of genus 0 or 1, which is not rescaled, divide by rho = 1 + O(t), which does not raise D). Taylor-expanding F_h(q_0 + delta),
the coefficient of t^(2-2h+k) of the vertex factor is a finite combination of Xi_h(R) with all parts j >= 1 and D(R) <= k.
(iii) Vanishing. Strip the parts equal to 1 with 1.2: Xi_h(1^p, R') = prod_{i<p}(J_{R'} + i) Xi_h(R'). If R' is empty and p >= 1
the factor J_{R'} = 0 kills it (Theorem C); if R' is non-empty, all its parts are >= 2 and D(R') = D(R), so it vanishes by
M_irr(h) whenever D(R) <= 2h-3. So a non-zero term needs D_v >= 2h_v - 2 at every singular vertex of genus >= 2 with r_v >= 1.
(iv) Counting. In a connected graph with an edge every vertex has a half-edge, so r_v >= 1. The t-order of a non-zero term is
>= Sum_{v singular, h_v >= 2} (2 - 2h_v + D_v) + (singular genus-0 vertices: +2 each; genus-1 and regular vertices: >= 0, by F1,
F3) >= 0. This proves (a). (v) The single-vertex graph has no edges and no loops: its only term with R empty is the order k = 0, which is
(-x c^3)^(g-1) F_g(q_0); all its other polar terms have R non-empty with D <= k <= 2g-3 and vanish by M_irr(g). The value is Theorem VI, as in (b).  QED
This replaces the Reduction Lemma of conifold_regularity_note.md s.5a (whose factorisation of an E_11 block between two singular
vertices of genus >= 2 into columns f_{1,b} is impossible with coefficients regular at t = 0, see the erratum there and
leaf_proof/leaf_structure_proof.md Remark 13.1; Theorem R needs no factorisation, because it expands rho^2 E_11 in t-independent leaves);
Conjecture V for one vertex with legs f_{j,b} and loops is the conclusion of steps (i)-(iii), leaf_structure.md step 4. The only properties of the singular idempotent it uses are
F1, L1-L3 and the Xi_h-zeros, and nothing about the regular idempotent beyond F1, F3.

Checks on the logic. (a) The leaves that occur at order t^k have integer scales j <= k+1, so only lattice points enter -- this
is why the "lattice, not hypersurface" nature of Conjecture M is harmless. (b) Only the vanishing direction of M is used.
(c) The ranges match exactly: polar orders k <= 2h-3 need D <= 2h-3, which is the range of M; the bound D <= k of Lemma C is
sharp (leaf_structure.md), and the first non-zero M-values (D = 2h-2) sit exactly at the order t^0. So there is no slack on
either side, but also no missing case. (d) Where it would break: a lam_0 component or a polynomial remainder in Lemma C (2.C),
or a failure of M_irr. M is sufficient, not shown to be necessary (polar coefficients are fixed combinations of Xi's).

Obstruction: none found. [Status 2026-11-02: see Addendum 3; the three items below are superseded.] What is unconditional:
 * genus 2 and 3: every input is a finite exact computation (F1, F3, L1-L3 from the stored z^9 factors; M_irr(2) proved, M_irr(3)
   exact), so regularity at genus <= 3 follows by this route (computer-assisted), independently of the direct graph sums;
 * genus 4: L1-L3 are verified through t^5 at all needed psi-orders; the route needs the 7 missing M_irr(4) points
   (1.9). The direct genus-4 graph sums (Stanage) already show regularity, so these points are a test of M, not of the theorem;
 * all genera: conditional on Lemma A (or Lemma D + Lemma B) and Lemma C. [Lemma C is proved (Addendum 2): conditional on Lemma A only.]

## 4. The most promising next step

1. Lemma D by global loop equations, prepared by an exact integer-tau TR code. Write an exact (Perl/BigRat) local TR for C_tau
   at integer tau >= 2: local coordinate zeta at v = tau with X = X(tau) - zeta^2/(2 tau(tau+1)), so that v(zeta) = tau + zeta +
   O(zeta^2) has rational coefficients; basis e_k = dxi_k with poles only at v = tau; the regular part of B in zeta enters the
   recursion as a B_{kl} matrix. Cost: the referee's tau = 1 code reached genus 6, n = 1 in 2.5 min. Uses: (a) by 1.7(a), exact
   zeros at h integer tau settle any given lattice point for ALL tau (mod P2): tau = 1, 2, ..., h settle every genus-h lattice point (genus 4:
   the 7 missing points; genus 5: all of M_irr(5)) without Mathematica; (b) it produces the data needed to find the global
   residue argument for Lemma D (which pair-of-sheets combinations have the degree bound at v = oo; the global linear loop
   equation of 2.D is the starting point).
2. In parallel, Lemma C from the Euler integral: derive the singular column r_{0,1}, r_{1,1} and E_{1j} order by order in t from
   tau^(a+b+1) B(a+1,b+1) 2F1(b+1, b+1; a+b+2; tau), showing simple poles at z = -c/j (j <= k+1) and decay at z = oo after the
   rho-normalisation. This is needed for every genus and is independent of M.
3. Cheap cross-check on Stanage (Mathematica/13.2.1): the 7 genus-4 irreducible points at 4 numeric values of tau (exact rational
   x, much cheaper than symbolic x); by 1.7(a) this proves them for all tau without using P2.
4. Write out Lemma B (referee route) so that every "mod P2" statement above becomes unconditional.
Not recommended as a first move: proving (6.1) (it only covers r = 1), or the framing cut-and-join flow from tau = 1 (by the variational formula the
tau-derivative is an integral of omega_{g,n+1} along the cut [-1, 0] of log(v/(1+v)), away from v = oo, where the induction
hypothesis gives no degree control; the compensating change of the expansion variable x would have to be tracked).

Scripts of this file (review_checks/, Perl + Math::BigRat, no Mathematica): status_check61.pl (evaluates (6.1)),
status_trC1.pl / status_trhalf.pl (copies of referee1_trC1.pl / referee1_trhalf.pl with env JEXTRA / JMAX for the j-range),
status_cmp61.pl (TR versus (6.1)). Example: JEXTRA=6 perl status_trC1.pl 11 6 7 > out.txt; perl status_cmp61.pl 1 out.txt.
