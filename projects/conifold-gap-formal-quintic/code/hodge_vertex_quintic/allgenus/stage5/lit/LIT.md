# Stage 5, task `lit`: literature pinning for (QRR), formal stationary phase and the chart mechanism (2026-11-14)

Worker note of the Stage-5 workflow (key `lit`). Continued in `LIT_2.md` (s.2 formal stationary phase, s.3 prior art for the chart
mechanism) and `LIT_3.md` (s.4 consequences for Stage 5, s.5 what was not read, s.6 record). AI-written; no human has read it. Labels:
PROVED, VERIFIED, CONJECTURAL, SKETCH, REFUTED. Notation of `../../MECHANISM_STAGE4.md`: colour 0, `u = 1 - 5^5 q`, `L = u^(-1/5)`,
`b_k = B_2k/(2k(2k-1))`, `h_w(z) = exp(-Sum_k b_k z^(2k-1) Sum_(a in w) a^(1-2k))`, `w_i = {zeta^i - zeta^j}_(j != i) U {-5 zeta^i}`.

## 0. Bottom line

* **(QRR) is in the literature, for the formal quintic itself, but only as a cited input.**
  * Guo-Janda-Ruan arXiv:1812.11908v1, eq. (17), inside the proof of Lemma 5.2, defines diagonal constants `C_alpha(z)` for the
    formal quintic (their "lambda-twisted theory", s.3.1 and eq. (9): twist parameter `t = 0`, `lambda_i = zeta^i lambda`). Their
    Birkhoff factorisation `S C_alpha^-1 = e^(u_alpha/z) R` gives `R(q=0) = C_alpha^-1 = h_(w_alpha)`: exactly the project's
    `diag(h_(w_i))`, to all orders in `z`, with the twist weight `-5 lambda_alpha`. GJR justify (17) only by citing Givental
    arXiv:math/0008067 and Coates-Givental arXiv:math/0110142.
  * Guo-Janda-Ruan arXiv:1709.07392v1, eqs. (19)-(20), in the proof of Lemma 4.6, state the same constant term for the theory with
    twist parameter `t` and base weights `lambda_i = xi^i lambda`. The twist weight there is `t - 5 lambda_i`; at `t = 0` it is the
    project's formula.
  * Lho-Pandharipande (arXiv:1702.06096 eq. (5); arXiv:1803.01409 eq. (5)) DEFINE the formal quintic. The second paper also gives the
    localisation vertex class, eq. (26). Neither paper writes down `R(q=0)`: their R_k are the I-function's own asymptotic series (for the quintic, the
    Zagier-Zinger series), which vanish at `q = 0` for `k >= 1`.
  * Givental proves the untwisted analogue (torus-equivariant theories with isolated fixed points): arXiv:math/0008067 Theorem 1 with
    its Remark; arXiv:math/0108100 Theorem 9.1. The proofs are by localisation and are partly sketched (s.1.3).
* **Sign of the twist weight: settled, three independent ways, all agreeing with the project** (s.1.4). Coates-Givental Cor. 1, the
  LP vertex class (26) with Mumford's relation, and GJR eq. (17) all put the O(5) weight into `h` as `-5 lambda_i`, i.e. as the
  factor `exp(+Sum_k b_k (z/(5 lambda_i))^(2k-1))`. My exact check: from LP (26) alone, `N_(2,0) = -19/9600` and
  `N_(3,0) = 1/224000`, the project's localisation values. The class identity holds modulo Mumford's relation in all degrees
  (g = 2, 3), and the weight `+5 lambda_i` gives `-97/28800` and `1/81000` instead (VERIFIED, exact over Q(zeta_5); s.1.5).
* **Formal stationary phase (LIT_2 s.2).** Citable, read statements exist for:
  * one- and several-variable Laplace / stationary phase (DLMF s.2.3; Etingof arXiv:2409.03117v2 Thms 2.16-2.17; Kirwin
    arXiv:0810.1700 Thm 1.1);
  * Wick and Feynman (Etingof Thms 3.1, 3.5; Polyak arXiv:math/0406251 Thms 2.7-2.8, Prop. 2.17);
  * Fubini for formal Gaussian integrals (Bar-Natan-Garoufalidis-Rozansky-Thurston arXiv:math/9801049 Prop. 2.13; diagrammatic);
  * Stirling's series (DLMF (5.11.1), (5.11.8); Olver as cited there).

  What is missing:
  * No source I read states a NUMBERED nonlinear change-of-variables lemma, a parameter lemma, or a delta lemma for FORMAL integrals
    with complex phases.
  * Hormander s.7.7, Olver, de Bruijn and Wong are books I could not open (LIT_2 s.5).
* **Prior art for the chart mechanism (LIT_2 s.3).** Nothing I read PROVES (Gamma) for the O(5)-twisted formal quintic.
  * The closest text is GJR arXiv:1812.11908v1 s.7. Its Remark 7.2 states the identity relating the formal quintic's LG oscillatory
    integral, the QRR constants `C_alpha` and the Zagier-Zinger series. That is the project's (Gamma3) modulo formal Fubini. It is
    stated without proof, and it is used for their Cor. 7.3-7.4 and Lemma 5.5.
  * The same s.7 PROVES (formal stationary phase evaluated at the orbifold point `L = 0`, Lemma 7.7, eq. (37)) a closed formula.
    For `m = 5`, `k = 1`, `lambda = 1` it is literally the project's `Ginf(z) = exp(Sum_l (-1)^(l+1) B_(5l+1)(1/5) z^(5l)/(l(5l+1)))`,
    which Stage 4 lists as CONJECTURAL.
  * Givental arXiv:math/0108100 s.10 gives the chart mechanism at `Q = 0` for untwisted equivariant CP^n as an outline: the chart is
    `x_0 = Q/(x_1...x_n)`, the integral factorises into Gamma integrals, and Lemma 10.4 is applied; continuity in Q is asserted.

## 1. (QRR): what each source says

### 1.1 Definition of the formal quintic (read in full: 1803.01409v2 PDF, 1702.06096 via its ar5iv rendering)

* **arXiv:1803.01409v2 (Lho-Pandharipande, "Holomorphic anomaly equations for the formal quintic").**
  * s.0.2, eq. (5): `Ntilde^GW_(g,d) = Int_[Mbar_g(P^4,d)]^vir e(R pi_*(S^-5))`. Here S is the universal subbundle `f^* O(-1)`, and
    footnote 3 says `S^-5 = (S^*)^5`.
  * The Euler class is defined after localisation, as `c_top(R^0)/c_top(R^1)` on each fixed locus, where both are bundles with
    moving weights.
  * The torus `T = (C^*)^5` has weights `lambda_0..lambda_4` on `H^0(P^4, O(1))`. The integral is homogeneous of degree 0 and
    becomes a rational number after `lambda_i = zeta^i`. Footnote 4 says `lambda_i = zeta^i lambda_0` may equally be used.
  * s.1.1 states the same convention as weights `-lambda_0..-lambda_m` on `C^(m+1)`, so the tangent weights at `p_j` are
    `lambda_j - lambda_i` (i != j).
* **arXiv:1702.06096 (Lho-Pandharipande, "Stable quotients and the holomorphic anomaly equation"), s.0.2, eq. (5).**
  * The same definition. The general twisted theories on P^m are in s.0.3, eq. (6), with the same weight convention and the
    specialisation `lambda_i = zeta^i lambda_0`.
  * The formal-quintic theorems are stated in s.0.6 (Thms 3, 4) and proved in 1803.01409.
  * So the Stage-4 citation "1702.06096 eq. (5)" is correct; 1803.01409 eq. (5) is the same formula. This supports Stage-4 K1.
* **Localisation vertex: 1803.01409v2 s.3.2, eq. (26).** The vertex class at `p_i` is
  `H_g^(p_i) = [e(E^* (x) T_(p_i)P^4) / e(T_(p_i)P^4)] * [(5 lambda_i) / e(E^* (x) (5 lambda_i))]`, with E the Hodge bundle. The
  first factor is the virtual normal bundle; the second is the twist, with weight `5 lambda_i` for O(5) at `p_i`.
* **What LP do NOT state.** In 1803.01409 s.2.4.4, eq. (18) gives the asymptotic form
  `Ibar|_(H=lambda_i) = e^(mu lambda_i/z) (R_0 + R_1 z/lambda_i + ...)`, with `R_k in C[[q]]`, derived in their reference [27] (= Zagier-Zinger).
  * Its examples `R_0 = L` and `R_1`, `R_2` multiples of `(L - L^5)` and `(1 - L^4)^2` vanish at `q = 0` for `k >= 1`.
  * Prop. 5 (attributed to Zagier-Zinger) gives `R_k in C[L^(+-1)]`. eq. (21) carries the same normalisation into `S_i(H^k)`.
  * This is the I-function's own normalisation (`Phi_s(0) = delta_(0,s)`, ZZ Thm 4(i); LIT_2 s.3.5), NOT the QRR-normalised
    R-matrix. LP never need `R(q=0)`: their vertices are Hodge integrals against (26), reduced to genus-0 data via
    Prop. 16 and Prop. 21-22.
  * The gamma worker's reading in Stage 4 (O11: "s.2.4.4 eq. (18) and Prop. 5") is therefore confirmed in its numbers. But eq. (18)
    is the ZZ normalisation; it does not contain the QRR constants.

### 1.2 The q = 0 R-matrix stated for the formal quintic (GJR; PDFs read)

* **arXiv:1812.11908v1** (Guo-Janda-Ruan, "Structure of higher genus Gromov-Witten invariants of quintic 3-folds"; only version, no
  journal reference on the abstract page).
  * s.3.1 calls the formal quintic of LP the "lambda-twisted theory": `t = 0` (allowed because `c_1(O(5))` is invertible
    equivariantly) and `lambda_i = zeta^i lambda`. Eq. (9) is the specialisation `(lambda_0..lambda_4, t) = (lambda, zeta lambda, ..., 0)`.
  * s.4.1 calls the theory a semisimple CohFT, reconstructed by Givental-Teleman.
  * Lemma 5.2 gives the quantum differential equations for S and R. Its proof uses the Birkhoff factorisation
    `S_(i alpha)(z) C_alpha(z)^-1 = e^(u_alpha/z) R_(i alpha)(z)`, citing [16, 12] (= math/0008067, math/0110142) for it,
    `C_alpha(z) = exp(Sum_(k>0, beta != alpha) b_k z^(2k-1)/(zeta^alpha lambda - zeta^beta lambda)^(2k-1))`
    `            * exp(Sum_(k>0) b_k z^(2k-1)/(-5 h_alpha)^(2k-1))`  (17),
    where `h_alpha = zeta^alpha lambda` is H at the fixed point.
  * Hence `R(q=0)_(alpha alpha) = C_alpha^-1 = h_(w_alpha)(z)`, with `w_alpha = {lambda_alpha - lambda_beta} U {-5 lambda_alpha}`: the
    project's (QRR), all orders in z, including the `z^(5k)` part.
  * The convention is Givental's (`S ~ R e^(u/z) C`, `R = 1 + O(z)`), the same as the project's (legs `R^-1(psi)`, edges
    `(eta^-1 - R^-1(psi) eta^-1 R^-1(psi')^T)/(psi + psi')`; `allgenus/model/ms_rows.py` docstring).
* **arXiv:1709.07392v1** (Guo-Janda-Ruan, genus-two quintic mirror theorem).
  * The main theory has twist parameter `t != 0` and a non-equivariant base. Its I-function, eq. (7), has numerator
    `prod_(j<=5d)(5H + jz - t)`.
  * In the proof of Lemma 4.6 GJR note that the standard q = 0 normalisation needs semisimplicity at q = 0. So they pass to the theory
    with `(C^*)^5` acting on the base, `lambda_i = xi^i lambda`, which is semisimple at q = 0.
  * eq. (19) relates the fundamental solution S to `Psi^-1 R e^(U/z)` through constant diagonal matrices, justified "by uniqueness of
    fundamental solutions of the quantum differential equation (see also [21, 13])" (= math/0008067, math/0110142).
  * eq. (20): `C(lambda, z) = diag(exp(Sum_(k>0, j != i) b_k z^(2k-1)/(lambda_i - lambda_j)^(2k-1)))` and
    `Gamma(t, z) = exp(Sum_(k>0) b_k z^(2k-1)/(t - 5H)^(2k-1))`. GJR say these two give the constant term of R.
  * Lemma 4.6 then writes the q = 0 diagonal of R through `z^3` with exactly these weights (twist weight `t - 5 xi^alpha lambda`, and
    `Sum_(j=1..4) 1/(1 - xi^j) = 2`, `Sum 1/(1 - xi^j)^3 = -1`). At `t = 0` this is `h_(w_alpha)` through `z^3`.
* **Status.** In both GJR papers (QRR) is an input justified by citation, not proved there. I read the statements; I did not find a
  separate proof in either paper.

### 1.3 The cited sources (read: math/0108100v2 and math/0008067v4 in PDF, math/0110142v2 in PDF and ar5iv)

* **Givental arXiv:math/0008067v4 ("Semisimple Frobenius structures at higher genus").**
  * s.2.2: the Hodge factors of the localisation vertices are absorbed by modifying the fundamental solution S. The new S is S times
    `exp(-Sum_k z^(2k-1) N_(2k-1)(1/chi^i) B_2k/(2k(2k-1)))`, where `N_l(1/chi^i)` is the Newton sum of the inverse TANGENT
    characters at the i-th fixed point.
  * Theorem 1 (s.2.2): with this S, formula (1) gives the higher-genus potentials of torus-equivariant GW theory with isolated fixed
    points.
  * The Remark after it: the new R(z) becomes the diagonal matrix of these compensating constants in the classical limit.
  * Example 6 checks genus 1. s.2.3 proves the Hodge-factor step by a Lemma, using Mumford's GRR ([24] there): both sides satisfy the
    same PDE system and agree at s = 0. Some identifications are referred to Givental's earlier papers [15, 16].
  * In the project's notation: `R(q=0)_i = exp(-Sum_k b_k z^(2k-1) Sum_(tangent a) a^(1-2k))`, i.e. `h` on the tangent weights.
* **Givental arXiv:math/0108100v2 ("Gromov-Witten invariants and quantization of quadratic hamiltonians").**
  * Theorem 9.1: in the torus-equivariant theory with isolated fixed points, the ancestor potential equals the semisimple formula,
    PROVIDED R is normalised at `Q = 0` by `R|_(Q=0) = exp diag(b_1..b_N)`. Here `b_i = Sum_k N^(i)_(2k-1) B_2k z^(2k-1)/(2k(2k-1))`,
    and `N^(i)_l` is the Newton sum of the inverse characters on the COTANGENT space at the fixed point.
  * Since cotangent characters are minus tangent weights and `1 - 2k` is odd, this equals the formula of math/0008067.
  * Remark 9.2(a): Thm 9.1 is equivalent (via 5.1) to Theorem 2 of math/0008067. Both are proved by fixed-point localisation, which
    relies on Graber-Pandharipande. The Bernoulli series offset the Hodge integrals via Thm 4.1 for X = point.
  * Thm 9.1 itself is not proved in math/0108100.
  * Theorem 4.1 (s.4): the Hodge potential with insertion `exp(Sum_k s_k ch_(2k-1)(E))` equals
    `exp(Sum_k B_2k/(2k)! s_k (z^(2k-1))^)` applied to the descendant potential. Remark 4.2 calls it a reformulation of
    Faber-Pandharipande and Mumford's GRR.
* **Coates-Givental arXiv:math/0110142v2 ("Quantum Riemann-Roch, Lefschetz and Serre").**
  * Theorem 1 (s.3, eq. (7)): for `c(.) = exp(Sum s_k ch_k(.))` it gives the twisted potential, up to explicit scalar factors, as
    `exp(Sum_(m>0, l>=0) s_(2m-1+l) B_2m/(2m)! (ch_l(E) z^(2m-1))^)` times a `1/z` term, applied to `D_0`.
  * Remark (2): for E trivial this is Givental's Thm 4.1, i.e. Mumford plus Faber-Pandharipande. The proof is in Appendix 1, by
    Mumford's GRR argument; its verifications are called tedious but straightforward.
  * s.4, eq. (9): the S^1-equivariant Euler class corresponds to `s_0 = ln lambda`, `s_k = (-1)^(k-1)(k-1)!/lambda^k`.
  * Corollary 1 (s.4): the Euler-class twist acts by `prod_i exp(Sum_(m>0) b_m ((z/(lambda + rho_i))^(2m-1))^)` (sign +), times the
    `1/z` factors, where `rho_i` are the Chern roots.
  * Cor. 3 (s.5) is the inverse Euler class (Serre). s.7, eq. (11) is the genus-0 Lagrangian-cone version (LIT_2 s.3.2).
* **Faber-Pandharipande arXiv:math/9810173 ("Hodge integrals and Gromov-Witten theory"; ar5iv rendering).**
  * s.1.1 states Mumford's GRR for ch(E) as an unnumbered "Theorem (Mumford)". Proposition 1 is its virtual-class version.
  * Proposition 2: `D_(2l-1) Z^X_E = 0`, the operator form of Givental's Thm 4.1.
  * Theorem 4: `Int_(Mbar_g) lambda_(g-1)^3 = |B_2g|/(2g) |B_(2g-2)|/(2g-2) 1/(2g-2)!`.
  * Mumford's relation `c(E) c(E^*) = 1` is unnumbered in s.2.1, with its consequence eq. (22).
* **Pandharipande arXiv:1712.02528 ("Cohomological field theory calculations"; ar5iv).** Definition 5 defines the R-matrix action.
  The unnumbered example after Theorem 7 (s.1.3) gives `R(z) = exp(-Sum b_k z^(2k-1))` for the CohFT `c(E)`. Footnote 11: that
  paper's R is Givental-Teleman's `R^-1`. In Givental's convention `c(E)` therefore has `R = exp(+Sum b_k z^(2k-1))`, which is the
  project's `h` at the single weight `a = -1`. This is consistent with the next subsection.

### 1.4 The sign of the twist weight (PROVED, given the statements of s.1.1-1.3)

**Claim.** In Givental's convention, the q = 0 vertex R-matrix at `p_i` of the formal quintic is `h_(w_i)`, with
`w_i = {lambda_i - lambda_j} U {-5 lambda_i}`. Equivalently, the O(5) factor is `exp(+Sum_k b_k (z/(5 lambda_i))^(2k-1))`, the
tangent factor is `exp(-Sum_k b_k Sum_(j != i) (z/(lambda_i - lambda_j))^(2k-1))`, and the twist has the opposite sign structure
to a tangent weight.

* **(a) LP (26), Mumford, Thm 4.1.**
  * Put `Lv(w) := e(E^* (x) L_w) = Sum_k (-1)^k lambda_k w^(g-k)`. Then (26) is
    `prod_(j != i) Lv(a_j)/a_j * (5 lambda_i)/Lv(5 lambda_i)`, with `a_j = lambda_i - lambda_j`.
  * Mumford's `c(E) c(E^*) = 1`, homogenised, reads `Lv(w) Lv(-w) = (-1)^g w^(2g)`. Hence
    `(5 lambda_i)/Lv(5 lambda_i) = (5 lambda_i)^(1-g) * Lv(-5 lambda_i)/(-5 lambda_i)^g`.
  * So `H_g^(p_i) = W_i^(g-1) prod_(a in w_i) Lv(a)/a^g`, with `W_i = prod_(j != i)(lambda_i - lambda_j)/(5 lambda_i)`. This is
    exactly the project's constant-map vertex (`FORMAL_QUINTIC_genus2.md` s.0: `(W_i^0)^(g-1) Hodge_(wt_i)`).
  * Next, `Lv(a)/a^g = prod_r (1 - x_r/a) = exp(-Sum_k (2k-2)! ch_(2k-1)(E) a^(1-2k))`, because `ch_2k(E) = 0` for `k >= 1`
    (Mumford; math/0108100 s.4).
  * Thm 4.1 of math/0108100 with `s_k = -(2k-2)! Sum_(a in w_i) a^(1-2k)` turns the insertion into the operator
    `exp(-Sum_k b_k Sum_(a in w_i) a^(1-2k) (z^(2k-1))^)`, i.e. into `h_(w_i)`.
  * The scalar `W_i^(g-1)` is the genus weight of the vertex (the metric / `Psi` part), not part of R.
* **(b) Coates-Givental Cor. 1.**
  * An Euler-class twist with Chern root `lambda + rho` contributes `exp(+Sum_m b_m (z/(lambda + rho))^(2m-1))`.
  * At `p_i` the T-equivariant root of O(5) is `5 lambda_i`. This is the localised substitution that GJR make when they cite [12].
  * The factor is `exp(+Sum b_m (z/(5 lambda_i))^(2m-1)) = exp(-Sum b_m z^(2m-1) (-5 lambda_i)^(1-2m))`.
  * The tangent part, `1/e(N^vir)`, carries the opposite sign: Givental Thm 9.1 and math/0008067 Thm 1.
* **(c) GJR arXiv:1812.11908v1 eq. (17)** gives the same thing explicitly: `(-5 h_alpha)` in `C_alpha` and `R(q=0) = C_alpha^-1`.
* **Against the project.**
  * `natural_form/referee/ref_rows.py`, `hcoef` (read-only), builds `w = [zeta^i - zeta^j] + [-5 zeta^i]` and
    `log h = -Sum_k b_k z^(2k-1) Sum_(a in w) a^(1-2k)`. Up to the overall factor, that is `C_alpha^-1` of GJR (17) at `lambda = 1`.
  * The Stage-4 statement "the sign of the weight `-5 zeta^i` is pinned by (H) at n = 1" is therefore confirmed independently, from
    the literature.
* **What (a)-(c) do not give.** That the R-matrix of the FULL formal-quintic CohFT restricts at q = 0 to these vertex R-matrices is
  quoted, not proved here.
  * The sources are GJR (17), and in the untwisted case Givental's localisation theorems.
  * The missing one-line argument: the q = 0 (degree-0) part of the theory is the direct sum of the five vertex theories (constant
    maps, `Mbar_(g,n)(P^4,0) = Mbar_(g,n) x P^4`). The R-matrix of a semisimple CohFT is unique given the CohFT. The second step
    is Teleman's classification (arXiv:0712.0160), which I did NOT open in this task.

### 1.5 Exact check (VERIFIED): `lit_vertex_check.py`, log `logs/lit_vertex_check.log`

New code: stdlib only, exact over Q(zeta_5), no project import, < 1 s.
* Inputs: LP (26) with `lambda_i = zeta^i`.
* Top-degree Hodge integrals are reduced by Mumford's relation to `lambda_1^(3g-3)`:
  * `g = 2`: `lambda_2 = lambda_1^2/2`;
  * `g = 3`: also `lambda_1 lambda_3 = lambda_1^4/8` and `lambda_3^2 = 0`.
* They are evaluated with FP Theorem 4: `Int_(Mbar_2) lambda_1^3 = 1/2880`, `Int_(Mbar_3) lambda_2^3 = 1/725760`.

| check | g = 2 | g = 3 |
|---|---|---|
| (A0) LP class vs project class as polynomials in FREE lambda's | differ (expected) | differ (expected) |
| (A) the same, modulo Mumford's relation, all degrees, i = 0..4 | PASS | PASS |
| (B) `N_(g,0) = Sum_i Int H_g^(p_i)` vs project's localisation value | `-19/9600` = `-19/9600` PASS | `1/224000` = `1/224000` PASS |
| (C) control: twist weight `+5 lambda_i` | `-97/28800` (differs) | `1/81000` (differs) |

* **What this shows.**
  * The project's numbers `F_2(0)`, `F_3(0)` (`FORMAL_QUINTIC_genus2.md` s.2) follow from LP's published definition.
  * The vertex identity of s.1.4(a) holds.
  * Already at genus 2, the opposite sign of the twist weight would change `F_2(0)`.
* **What it does not show.** It does not test the project's R-matrix code. The step from the Hodge class to `h` is Thm 4.1, a
  quoted theorem. `N_(3,0)` involves `ch_5(E)`, i.e. the `z^5` level of (26). But the comparison value is itself a localisation
  number, so this is NOT the critic's suggested test of (QRR) at `z^5` (the q = 0 Givental-Teleman graph sum on `Mbar_3` with
  `R = diag(h)`); that test was not done.
