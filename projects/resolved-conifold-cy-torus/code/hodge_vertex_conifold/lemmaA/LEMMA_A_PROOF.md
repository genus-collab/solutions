# Lemma A: all-genus proof of record

File: hodge_vertex_conifold/lemmaA/LEMMA_A_PROOF.md. Written 2026-11-02 (round 5, plan item P1 of lemmaA/round2/round2_status.md s.6). Edited in round 6 with the referee fixes of s.11.4 (round-5 text: LEMMA_A_PROOF_round5.md).

Purpose: the single written proof of Lemma A = M_irr(h) for non-empty R, for the weights (1, tau, -1-tau), usable by both the resolved conifold x C^2 proof and the local P2 x C^2 proof.

Tags: **[P]** proved in this document. **[C]** cited; s.9 says what was and was not confirmed. **[E]** evidence, used by no proof step. **[O]** open.

Provenance: the proof combines the round-5 reports MVmatch (s.2, s.8) and ThmL (s.3-s.6), each with two referee reports (verdict "sound with minor fixes", no major item). All referee fixes are applied; s.11 lists them. The writer re-derived s.2-s.6 line by line during assembly. The assembled text was refereed in round 6 under two lenses, logic and literature (L_logic) and independent computation (L_comp): both "sound with minor fixes", no major item; the fixes are applied (s.11.4, and O1 in s.10.2). (MV) itself stays cited and not re-proved (O3).

Rulings respected: no TR for any target theory, and no TR anywhere in the proof. No KP theorem is cited: the one integrable identity needed is derived in s.4. No finite-order computation fixes a formula or a constant: s.12 is evidence only.

---

## 0. Result

**Definitions** (full notation in s.1). tau is an indeterminate, sigma = tau(1+tau), Lambda^vee(w) = Sum_i (-1)^i lambda_i w^{h-i}, and

    Omega_h = Lambda^vee(1) Lambda^vee(tau) Lambda^vee(-1-tau) / (-sigma)^h,
    lam_j(psi) = 1/(1 - j psi),     T(psi) = -psi^2/(1 - psi).

For a finite multiset U of numbers, with r = #U:

    Xi_h(U) := Sum_{m>=0} (1/m!) Int_{Mbar_{h,r+m}} Omega_h prod_{j in U} lam_j(psi) prod_{k=1}^m T(psi_{r+k}),

a finite sum because T has psi-degree >= 2. D(U) := Sum_{j in U} (j - 1).

**Theorem L (Lemma A = M_irr(h)).** [P modulo the Marino-Vafa formula] For every h >= 2 and every **non-empty** multiset R of integers >= 2 with D(R) <= 2h - 3,

    Xi_h(R) = 0   in Q(tau).

**Corollary L1 (specialisation).** [P modulo (MV)] Write Lambda(t) = Sum_i t^i lambda_i, and let alpha_i be the Chern roots of the Hodge bundle. Then Xi_h(R; tau) = Sum_{k=0}^h (-1/sigma)^k Xi^(k)_h(R), with Xi^(k)_h(R) = Sum_m (1/m!) Int Lambda(-1) e_k(alpha + alpha^2) lam_R T^m independent of tau. For R as in Theorem L, every Xi^(k)_h(R) vanishes. Hence Xi_h(R; tau_0) = 0 at every tau_0 in a field of characteristic 0 (in particular at every complex tau_0) with sigma(tau_0) != 0, including the diagonal tau_0 = -1/2 and all other zeros of the factors Pi_k used in the proof.

**Corollary L2 (parts >= 1).** [P modulo (MV)] The same vanishing, in Q(tau) and at every such tau_0, holds for every non-empty multiset U of positive integers with D(U) <= 2h - 3. No part 0 is allowed.

**Corollary L3 (general scale).** [scale law P; vanishing P modulo (MV)] Let E, w_2, w_3 lie in a field of characteristic 0, with E + w_2 + w_3 = 0 and E w_2 w_3 != 0. Put Omega_w = prod_a Lambda^vee(w_a)/w_a^h for w = (E, w_2, w_3), lam^E_j = 1/(E - j psi), T^E = -psi^2/(E - psi). The corresponding integrals satisfy Xi^E_h(U) = E^{-(3h-3+2r)} Xi_h(U; tau = w_2/E), hence vanish in the same range. This covers the conifold and local P2 conventions (s.7.4, s.7.5).

**Excluded.** R empty: Xi_h(empty) is the gap, and part 0 reduces to it, Xi_h(U + 0) = (D(U) - 2h + 2) Xi_h(U) (s.7.2). See s.8.5 for what is known about Xi_h(empty). Nothing in the proof uses it.

**Inputs.** Exactly one non-textbook input:

- **(MV)** the Marino-Vafa formula, Liu-Liu-Zhou, arXiv:math/0306434v3, Introduction, eqs. (1)-(3) (called "the Marino-Vafa formula" there, and not one of their numbered theorems; their Theorems 1, 2 are the cut-and-join equations).

Textbook inputs: string and dilaton equations with classes pulled back along forgetful maps; the genus-0 multinomial formula for psi integrals; Frobenius formula and Murnaghan-Nakayama rule; the hook-length form of the principal specialisation of Schur functions; Lagrange inversion.

Not used: TR, Eynard, ABDKS, BKMP, any cited KP or 2D Toda theorem, the MV cut-and-join, Route A Theorems B, D, F and Cor. E, any conifold or local-P2 input, any computer output.

**Mechanism.** The pair Z_0 = Z_MV(p), Z_1 = Z_MV(e^{tau x} o p) satisfies one explicit bilinear equation P_k(H) Z_0 . Z_1 = 0 for each k >= 2. The shift enters through e^{beta Eu} with Eu = Sum k p_k d/dp_k and beta = tau x. At the dilaton point (y = u e^{-u}, w = 1/(1-u), w -> 0 is s = -1), Eu raises the order in w by one with coefficient delta = D - 2g + 2, so powers of Eu produce falling factorials (delta)_k. The lowest coefficient of the equation in w is an exact all-genus relation (Theorem M_k). Inside the Lemma A range a counting identity kills every term except a non-zero multiple of the target value.

---

## 1. Notation

**Field.** All identities hold over Q(tau) unless stated otherwise. Specialisation happens only in s.7.1.

- sigma = tau(1+tau).
- Pi_m = prod_{a=1}^{m-1} (m tau + a)/(m-1)!, so Pi_1 = 1, and Pi_m != 0 in Q[tau]. Pi_mu = prod_i Pi_{mu_i}.

**Hodge classes.** Lambda(t) = Sum_i t^i lambda_i, and Lambda^vee(w) = Sum_i (-1)^i lambda_i w^{g-i} = w^g Lambda(-1/w). So

    Omega_g = Lambda(-1) Lambda(-1/tau) Lambda(1/(1+tau)),   Omega_0 = 1.

Omega_g is pulled back under forgetful maps.

**Leaves.** lam_j = 1/(1 - j psi), lam_0 = 1, T = -psi^2/(1-psi), Tp := -T = psi^2/(1-psi) = lam_1 - 1 - psi.

**Multisets.** U, V are finite multisets, R is a multiset of integers >= 2. r_U = #U, J_U = Sum u, D(U) = J_U - r_U, delta_g(U) = D(U) - 2g + 2. U + V is the union, 1^a is a parts equal to 1. |Aut U| is the product of the factorials of multiplicities. y_U = prod_{u in U} y_u.

**Factorials.** (a)_k = a(a-1)...(a-k+1) falling, (a)^(k) = a(a+1)...(a+k-1) rising, both 1 at k = 0.

**Integrals.**

- For a multiset mu of positive integers, l = l(mu) >= 1: M_g(mu) = Int_{Mbar_{g,l}} Omega_g prod_i lam_{mu_i}(psi_i) if 2g-2+l > 0; and M_0(mu) := |mu|^{l-3} for all l >= 1. For l >= 3 this is the integral, by Int_{Mbar_{0,l}} prod psi_i^{d_i} = (l-3)!/prod d_i!. For l = 1, 2 it is the convention of LLZ: M_0((m)) = 1/m^2, M_0((m,n)) = 1/(m+n).
- For 2g-2+r > 0, r = r_U, N = 3g-3+r:

      Q_g(s;U) := Sum_{m>=0} (s^m/m!) Int_{Mbar_{g,r+m}} Omega_g lam_U Tp^{tensor m},

  a polynomial in s of degree <= N, since 2m <= 3g-3+r+m.
- X_g(U) := Xi_g(U) = Q_g(-1;U) for g >= 1 (U non-empty when g = 1). This is the definition, as T = -Tp.
- X_0(U) := Gamma(J_U)/Gamma(J_U - r_U + 3) for U non-empty with positive parts. Explicitly: 1/(J(J+1)) for r = 1, 1/J for r = 2, (J-1)(J-2)...(J-r+3) for r >= 3. X_0 is defined by this formula and is not a Xi. Lemma 3.3 shows why it plays the same role.

**Generating function.** q = e^x, p_k(q^rho) = 1/(q^{k/2} - q^{-k/2}) = 1/(2 sinh(kx/2)), f_2(nu) = sum of contents of nu, kappa_nu = Sum_i nu_i(nu_i - 2i + 1) = 2 f_2(nu), s_nu = Sum_mu chi_nu(mu) p_mu/z_mu, and

    Z(p) = Sum_nu e^{tau x f_2(nu)} s_nu(q^rho) s_nu(p),     F = log Z,   in Q[tau]((x))[[p_1, p_2, ...]].

Z has constant term 1 in p, so F is defined and has no constant term in p.

**Local symbols of s.4-s.5.** t_k = p_k/k are the times. H_k is the Hirota derivative in t_k. Sch_j is the elementary Schur polynomial, exp(Sum t_k z^k) = Sum_j Sch_j(t) z^j. Eu = Sum_k k p_k d/dp_k is the Euler operator, and beta = tau x.

The notations of the earlier reports differ. Dictionary: this document's y_k, theta, Eu, Z_0, Z_1, H_k, Sch_j, epsilon are R3's and ThmL's q_k, D, E, tau_0, tau_1, D_k, p_j, u (Hirota shift). MVmatch's C is F here.

---

## 2. The input: Marino-Vafa in real form

### 2.1 The cited statement [C]

LLZ, arXiv:math/0306434v3, Introduction. Let I_g(mu) = Int_{Mbar_{g,l}} Lambda^vee(1) Lambda^vee(-tau-1) Lambda^vee(tau) / prod_i (1 - mu_i psi_i), extended to l < 3 by I_0(mu) = |mu|^{l-3}, and

    C_{g,mu}(tau) = -( sqrt(-1)^{|mu|+l} / |Aut mu| ) [tau(tau+1)]^{l-1} Pi_mu I_g(mu).        (LLZ-def)

Then, with s_nu as above,

    Sum_{mu != empty} Sum_{g>=0} lambda^{2g-2+l} C_{g,mu}(tau) p_mu
        = log Sum_nu e^{sqrt(-1)(tau+1/2) kappa_nu lambda/2} V_nu(lambda) s_nu(p),            (LLZ-1)

    V_nu = prod_{a<b} sin[(nu_a-nu_b+b-a) lambda/2] / sin[(b-a) lambda/2]
           * prod_i prod_{v=1}^{nu_i} 1 / (2 sin[(v-i+l(nu)) lambda/2]).                      (LLZ-3)

Here the sum on nu has |nu| = |mu| at the coefficient of p_mu, and every power of lambda on the left has exponent 2g-2+l with g >= 0.

### 2.2 Real form

**Proposition 2.1.** [P modulo (LLZ-def), (LLZ-1), (LLZ-3)] Z(p) = Sum_nu e^{(tau+1/2) kappa_nu x/2} prod_{box in nu} (2 sinh(h x/2))^{-1} s_nu(p), and for every mu != empty with l = l(mu):

    [x^{2g-2+l}] [p_mu] F = sigma^{g+l-1} Pi_mu M_g(mu) / |Aut mu|     (g >= 0),         (2.1)

and all other powers of x vanish in [p_mu] F.

*Proof.*

(a) Hook form of V_nu. Let f(n) = 2 sin(n lambda/2), [n]! = f(1)...f(n), l = l(nu), beta_i = nu_i + l - i. Then nu_a - nu_b + b - a = beta_a - beta_b, prod_{v=1}^{nu_i} f(v+l-i) = [beta_i]!/[l-i]!, and prod_i [l-i]! = prod_{a<b} f(b-a). So (LLZ-3) reads V_nu = prod_{a<b} f(beta_a - beta_b) / prod_i [beta_i]!. Row-hook lemma: for fixed i the numbers {h_{i,j}}_j together with {beta_i - beta_k}_{k>i} are exactly {1, ..., beta_i}. (In the Maya diagram with occupied positions beta_k, the positions below beta_i are either occupied, by some beta_k with k > i, or empty; sliding the particle at beta_i down to an empty position removes a rim hook, and the lengths of these rim hooks are the hook lengths of row i.) Hence V_nu = prod_box 1/f(h). This is a multiset identity, valid for any non-vanishing f.

(b) Principal specialisation [textbook]. s_nu(1, t, t^2, ...) = t^{n(nu)} prod_box (1 - t^h)^{-1}. At t = q^{-1}, with the prefactor q^{-|nu|/2} from q^rho = (q^{-1/2}, q^{-3/2}, ...), Sum_box h = n(nu) + n(nu') + |nu| and kappa_nu = 2(n(nu') - n(nu)) give s_nu(q^rho) = q^{kappa_nu/4} prod_box 1/(2 sinh(hx/2)). Since e^{tau x f_2} q^{kappa/4} = e^{(tau+1/2) kappa x/2}, this is the first claim. It is an identity of rational functions of q^{1/2}, with p_k(q^rho) as defined in s.1.

(c) Removal of sqrt(-1) =: i. Work over Q(i), with x = i lambda, so lambda = -i x. Then e^{i(tau+1/2) kappa lambda/2} = e^{(tau+1/2) kappa x/2}, and 2 sinh(hx/2) = 2i sin(h lambda/2), so 1/(2 sin(h lambda/2)) = i/(2 sinh(hx/2)) and V_nu = i^{|nu|} prod_box 1/(2 sinh(hx/2)). By homogeneity i^{|nu|} s_nu(p) = s_nu((i^k p_k)_k). So the right side of (LLZ-1) is F(x; (i^k p_k)_k), whose [p_mu] coefficient is i^{|mu|} [p_mu] F.

(d) Signs. On the left of (LLZ-1), lambda^{2g-2+l} C_{g,mu} carries -i^{|mu|+l} (-i)^{2g-2+l} = -i^{|mu|} (i * (-i))^l (-i)^{2g-2} = -i^{|mu|} (-1)^{g-1} = i^{|mu|} (-1)^g. Cancelling i^{|mu|}:

    [x^{2g-2+l}] [p_mu] F = (-1)^g sigma^{l-1} Pi_mu I_g(mu) / |Aut mu|.

Finally Lambda^vee(1) Lambda^vee(tau) Lambda^vee(-1-tau) = (-sigma)^g Omega_g, so I_g(mu) = (-sigma)^g M_g(mu) and (-1)^g (-sigma)^g = sigma^g. In genus 0, Omega_0 = 1 and the extension of (LLZ-def) is M_0. QED

*Second derivation (not of record).* Okounkov-Pandharipande, Geom. Topol. 8 (2004), Theorem 1 (0.8) with a = tau, u = x, t = (-1, -1/a, 1/(a+1)) (the class is exactly Omega_g), binom((a+1)m, m) = (a+1) Pi_m(a), and the rescaling sqrt(a(a+1)) u, gives (2.1) again. Caution [P, hand check at mu = (1,1)]: the factor 1/|Aut mu| printed in OP (0.5), (0.8) is inconsistent with OP's own disconnected series at repeated parts: at mu = (1,1), (0.5) has left side u^{-4}(1 + ...) and right side (1/2) u^{-4} cosh u. OP's operator forms (1.24), (1.28) have no such factor, and the derivation uses that reading. The LLZ route is unaffected and is the route of record.

*Consistency in all genera, not used in the proof.* (i) tau = 0: Z = exp(Sum_d p_d/(d * 2 sinh(dx/2))), so [p_d] F = 1/(2d sinh(dx/2)) and no l >= 2 terms occur; on the Hodge side Lambda^vee(1) Lambda^vee(-1) = (-1)^g and Lambda^vee(0) = (-1)^g lambda_g, which matches (the (-1)^g above) in every genus. (ii) nu -> nu' gives Z(x; -1-tau; p) = Z(x; tau; ((-1)^{k-1} p_k)), matching Pi_m(-1-tau) = (-1)^{m-1} Pi_m(tau) and the invariance of Omega_g and sigma.

### 2.3 Scaled variables

Over Q(tau) put

    hbar = sigma x^2,   y_k = x sigma Pi_k p_k (k >= 1),   y := y_1 = x sigma p_1.

This uses Pi_k != 0 in Q(tau). The substitution rescales each monomial by a unit of Q(tau)((x)). Since x^{2g-2+l} sigma^{g+l-1} Pi_mu p_mu = hbar^{g-1} y_mu, and |Aut(R + 1^a)| = |Aut R| a! for R with parts >= 2, (2.1) is equivalent to

    F = Sum_{g>=0} hbar^{g-1} F_g,
    F_g = Sum_{R, parts >= 2, R != empty} y_R G^(g)_R(y) / |Aut R| + Sum_{a>=1} M_g(1^a) y^a / a!,    (2.2)

where for a multiset U of positive integers

    G^(g)_U(y) := Sum_{a>=0} M_g(U + 1^a) y^a / a!     (U non-empty, g >= 0; or U empty, g >= 2).   (2.3)

The empty-R part of F_g has no constant term, because F has none. Its a = 0 term would be Int_{Mbar_g} Omega_g, which is not in the series. G^(g)_empty with its a = 0 term is an auxiliary function used only in Cor. 3.5, where d/dy kills the constant. For all g >= 0,

    d/dy [ Sum_{a>=1} M_g(1^a) y^a/a! ] = G^(g)_(1),     d/dy G^(g)_U = G^(g)_{U+1}.          (2.4)

---

## 3. Dilaton-point calculus [P]

These are Route A Prop. 2.1, Cor. 2.2 and Lemma 3.1 (lemmaA/A_mv/report.md s.2, s.3.1), re-proved here so that this document stands alone, and extended uniformly to parts equal to 1.

Let u be a formal variable, y = u e^{-u} (an automorphism of the power series ring), w = 1/(1-u) in Q(tau)[[u]], s = w - 1, so u = s/(1+s). Let theta = y d/dy. Then theta = (w-1) d/du and d/du w^n = n w^{n+1}. Also 1/y = e^u w/(w-1).

w is transcendental over Q(tau): a polynomial relation in w would be a rational relation in u. So e^{uJ} Q(tau)[w, 1/w] is a free rank-one module over the Laurent polynomial ring in w, and coefficients in w are well defined. This purely formal fact replaces every analytic continuation to the second sheet in the earlier reports.

**Lemma 3.1 (operators).**

    (J + theta)(e^{uJ} w^n) = e^{uJ} ( (J - n) w^{n+1} + n w^{n+2} ),
    d/dy (e^{uJ} w^n) = e^{u(J+1)} ( J w^{n+1} + n w^{n+2} ),
    (d/dy) theta = (theta + 1) (d/dy).

*Proof.* Direct computation from the formulas above. QED

**Lemma 3.2 (dilaton leaves as unit legs).** Let U be a multiset of numbers (parts 0, 1, non-integers allowed), r = r_U, J = J_U, with 2g-2+r > 0, and put M_a = Int_{Mbar_{g,r+a}} Omega_g lam_U lam_1^a, G_U(y) = Sum_{a>=0} M_a y^a/a!. In Q(tau)[[s]], with u = s/(1+s):

    Q_g(s;U) = (1+s)^{-(2g-2+r)} e^{-uJ} G_U(u e^{-u}).                                       (3.1)

Equivalently G_U(u e^{-u}) = e^{uJ} w^{2g-2+r} Q_g(w-1;U), where Q_g is a polynomial of degree <= N = 3g-3+r, and X_g(U) = Q_g(-1;U) is its constant term in w.

*Proof.* Write Tp = lam_1 - 1 - psi and expand the m-th term of Q_g:

    Q_g = Sum_{a,b,c>=0} (s^a/a!) ((-s)^b/b!) ((-s)^c/c!) < Omega_g lam_U lam_1^a 1^b psi^c >.

Omega_g is pulled back along pi: Mbar_{g,n+1} -> Mbar_{g,n}, and psi_i = pi^* psi_i + D_{i,n+1}, psi_{n+1} D_{i,n+1} = 0, pi_* psi_{n+1} = 2g-2+n. So:

- The dilaton equation removes the c factors psi with total factor (e)^(c), where e = 2g-2+r+a+b.
- The string equation acts on legs by psi^n -> psi^{n-1}, killing psi^0. It sends lam_j to j lam_j. So removing one leg 1 multiplies by the sum of the labels, and removing the b legs 1 gives (J+a)^b.

All intermediate spaces are stable since 2g-2+r > 0. Resum: Sum_c ((-s)^c/c!) (e)^(c) = (1+s)^{-e}, and Sum_b ((-s)^b/b!) (J+a)^b (1+s)^{-b} = e^{-u(J+a)}. Therefore

    Q_g = Sum_a (s^a/a!) (1+s)^{-(2g-2+r+a)} e^{-u(J+a)} M_a = (1+s)^{-(2g-2+r)} e^{-uJ} Sum_a (u e^{-u})^a M_a/a!.

Each s-coefficient is a finite sum, so this holds in Q(tau)[[s]]. The degree bound comes from the definition of Q_g. No analytic continuation of G_U in y is used: X_g(U) is the value at s = -1 of a polynomial whose Taylor series at s = 0 is the right side. QED

*Laguerre form* (used only for s.8.2 and in the evidence). With L_n^(al)(x) = Sum_i (-1)^i binom(n+al, n-i) x^i/i!,

    [s^m] Q_g = Sum_{a<=m} (M_a/a!) (-1)^{m-a} L_{m-a}^{(2g-3+r+a)}(-J-a),
    X_g(U) = Sum_{a=0}^{N} ((-1)^a/a!) L_{N-a}^{(2g-2+r+a)}(-J-a) M_a.

The right side of the first line vanishes for m > N.

**Lemma 3.3 (genus 0).** Let U be non-empty with positive parts, r = r_U, J = J_U. Then e^{-uJ} G^(0)_U = w^{r-2} (polynomial in w) with constant term X_0(U).

*Proof.* G^(0)_U = Sum_a (J+a)^{a+r-3} y^a/a!, and J + theta multiplies the coefficient of y^a by J + a. So J + theta is injective on power series and raises r by one in this formula. By Lagrange inversion for u = y e^u, Sum_a (J+a)^{a-1} y^a/a! = e^{uJ}/J, which is r = 2. Lemma 3.1 then gives r >= 3, with leading coefficient prod_{n=1}^{r-3} (J - n). For r = 1 the solution of (J + theta) G = e^{uJ}/J is

    G^(0)_(J) = e^{uJ} ( w^{-1}/(J(J+1)) + 1/(J^2(J+1)) ),

by Lemma 3.1 and injectivity. QED

**Lemma 3.4 (uniform leading coefficient).** For all g >= 0, k >= 0 and every **non-empty** multiset U of positive integers:

    e^{-u J_U} (J_U + theta)^k G^(g)_U = w^{2g-2+r_U+k} ( (delta_g(U))_k X_g(U) + w * (polynomial in w) ).    (3.2)

*Proof.* k = 0 is Lemma 3.2 (g >= 1) or Lemma 3.3 (g = 0). The step k -> k+1 is Lemma 3.1 at n = 2g-2+r_U+l, where J - n = delta_g(U) - l, applied to the lowest term; higher terms stay higher. QED

**Corollary 3.5 (string).** Let g >= 1 and U non-empty with positive parts. Then X_g(U + 1) = J_U X_g(U). For g >= 2, X_g(1^v) = 0 for all v >= 1. Consequently, for U = U' + 1^v with U' having parts >= 2, and g >= 2:

    X_g(U) = (J_{U'})^(v) X_g(U') if U' is non-empty,     X_g(1^v) = 0.                         (S)

*Proof.* G_{U+1} = d/dy G_U, by (2.3). Apply the second line of Lemma 3.1 to (3.1): the coefficient of w^{2g-2+r+1} gives X_g(U+1) = J_U X_g(U). For g >= 2 and U empty, apply the same with the auxiliary G_empty (a = 0 term included; (3.1) holds for it, and d/dy kills the constant): J = 0 gives X_g(1) = 0. Iterate. QED

**Remark 3.6 (the empty-profile table).** R3 Lemma 2.3 is now a corollary and is not used separately. Put U = 1^v, v >= 1, so delta_g = 2-2g. In (3.2): X_0(1^v) = (v-1)!/2, with leading coefficients (v-1)!/2, (v-1)!, (v-1)!, 0, 0, ... at k = 0, 1, 2, 3, ...; at g = 1 the coefficient is (v-1)! X_1(1) at k = 0 and 0 for k >= 1, with X_1(1) = -(1+1/sigma)/24; at g >= 2 it is 0.

---

## 4. The bilinear equations P_k, for all k [P]

Self-contained. No integrable-systems theorem is cited; s.9 lists the literature that corroborates the result.

**Conventions.** v_i (i in Z) is the wedge basis, psi_i = v_i wedge, psi*_i = d/dv_i, |n> = v_{n-1} wedge v_{n-2} wedge ..., and F_n is the charge-n Fock space. |nu, n> has occupied modes nu_j + n - j, and is an orthonormal basis of F_n. Put alpha_k = Sum_j :psi_j psi*_{j+k}:. Then

    [alpha_k, psi_i] = psi_{i-k},   [alpha_k, psi*_i] = -psi*_{i+k},   [alpha_k, alpha_l] = k delta_{k+l,0},
    alpha_k |n> = 0 = <n| alpha_{-k}   (k > 0).

Put Gamma_+(t) = exp(Sum_{k>=1} t_k alpha_k), Gamma_-(t') = exp(Sum_{k>=1} t'_k alpha_{-k}), with p_k = k t_k and p'_k = k t'_k. xi(t,z) = Sum t_k z^k, and [1/z] is the time vector with k-th entry z^{-k}/k. A is a formal variable; t' are formal variables. All sums below are finite at fixed degree in t, t~, t'.

**(F1) Schur functions.** alpha_{-k} moves one particle up by k with sign (-1)^{height}. This is the Murnaghan-Nakayama rule, so alpha_{-mu}|n> = Sum_nu chi_nu(mu) |nu,n>. By the Frobenius formula,

    <n| Gamma_+(t) |nu,n> = s_nu(t),     Gamma_-(t') |n> = Sum_nu s_nu(t') |nu,n>.

Since the character table is invertible, F_n is spanned by the alpha_{-mu}|n>.

**(F2) Content weight.** Let W be diagonal with W|S> = A^{Lambda(S)}|S>, where Lambda(S) = Sum_{j in S, j>=0} f(j) - Sum_{j not in S, j<0} f(j) and f(j) = j(j+1)/2. Adding j to S changes Lambda by f(j) in both cases, so W psi_i W^{-1} = A^{f(i)} psi_i and W psi*_i W^{-1} = A^{-f(i)} psi*_i. Since f(nu_j + m) - f(m) at m = n-j equals nu_j(nu_j - 2j + 1)/2 + n nu_j, and Lambda(|0>) = Lambda(|1>) = 0, one has for n = 0, 1:

    Lambda(|nu,n>) = f_2(nu) + n |nu|.

Put Gop := W Gamma_-(t'). Then for n = 0, 1:

    Z_n(t) := <n| Gamma_+(t) Gop |n> = Sum_nu A^{f_2(nu) + n|nu|} s_nu(t') s_nu(t).

At A = e^beta, beta = tau x, t' = q^rho (that is p'_k = 1/(2 sinh(kx/2))): Z_0 = Z and Z_1 = Z(e^beta o p), where (e^beta o p)_k = e^{k beta} p_k. No normalising constant occurs, and the orientation e^{+tau x} is derived here, not fitted.

**(F3) Commutation.** S := Sum_i psi_i tensor psi*_i is a finite sum on any tensor of wedges. It commutes with Gop tensor Gop:

- for W termwise, by (F2);
- for Gamma_-, because Sum_i ([alpha_{-k}, psi_i] tensor psi*_i + psi_i tensor [alpha_{-k}, psi*_i]) = Sum_i (psi_{i+k} tensor psi*_i - psi_i tensor psi*_{i-k}) = 0.

On |0> tensor |1> only i = 0 survives, and S(|0> tensor |1>) = |1> tensor |0>. Hence S(Gop|0> tensor Gop|1>) = Gop|1> tensor Gop|0>.

**(F4) Vertex operators as shifts.** Put psi(z) = Sum psi_i z^i and psi*(z) = Sum psi*_i z^{-i}. Then [alpha_k, psi(z)] = z^k psi(z) and [alpha_k, psi*(z)] = -z^k psi*(z), so

    Gamma_+(t) psi(z) = e^{xi(t,z)} psi(z) Gamma_+(t),     Gamma_+(t) psi*(z) = e^{-xi(t,z)} psi*(z) Gamma_+(t).

The covector <1|psi(z) on F_0 is an eigenvector of right multiplication by alpha_{-k} with eigenvalue -z^{-k}, and it takes the value 1 on |0>. The covector <0|Gamma_+(-[1/z]) has the same two properties. By (F1),

    <1| psi(z) = <0| Gamma_+(-[1/z])   on F_0.

Likewise, with eigenvalue +z^{-k} and value 1 on |1>,

    <0| psi*(z) = <1| Gamma_+([1/z])   on F_1.

**(F5) Bilinear identity.** Note S = [z^0] psi(z) tensor psi*(z). Pair (F3) with <1|Gamma_+(t) tensor <0|Gamma_+(t~), move Gamma_+ past psi(z), psi*(z) by the first line of (F4), and apply the two covector identities to Gamma_+(t) Gop|0> and Gamma_+(t~) Gop|1>:

    [z^0] e^{xi(t - t~, z)} Z_0(t - [1/z]) Z_1(t~ + [1/z]) = Z_1(t) Z_0(t~).                      (4.1)

**(F6) The linear members.** Define Hirota derivatives by Z_0(b + eps) Z_1(b - eps) = exp(Sum eps_k H_k) Z_0 . Z_1 at base point b. Put t = b + eps, t~ = b - eps in (4.1):

    Sum_{j>=0} Sch_j(2 eps) Sch_j(-H~) exp(Sum eps_k H_k) Z_0 . Z_1 = exp(-Sum eps_k H_k) Z_0 . Z_1,   H~ = (H_k/k)_k.

Since Sch_j(2 eps) = 2 eps_j + (degree >= 2 in eps), the part linear in eps_k at eps = 0 gives H_k + 2 Sch_k(-H~) = -H_k on Z_0 . Z_1. That is:

    **P_k(H) Z_0 . Z_1 = 0,   P_k := H_k + Sch_k(-H~),   for every k >= 1.**                     (4.2)

Explicitly P_k = Sum_{alpha |- k} c_alpha H^alpha, where for alpha = (i^{m_i}):

    c_alpha = prod_i (-1/i)^{m_i}/m_i! + [alpha = (k)],    so c_(k) = (k-1)/k, which is != 0 for k >= 2.     (4.3)

P_1 = 0 identically. Examples: 2 P_2 = H_1^2 + H_2; 6 P_3 = 4 H_3 + 3 H_1 H_2 - H_1^3; 24 P_4 = H_1^4 + 8 H_1 H_3 - 6 H_1^2 H_2 + 3 H_2^2 + 18 H_4.

**Specialisation.** At a fixed t-monomial of weight d, the left side of (4.2) is a polynomial in t' of weight d + k with coefficients in Q[A, 1/A]. So (4.2) holds under the ring map A -> e^{tau x}, p'_k -> 1/(2 sinh(kx/2)), coefficientwise in Q[tau]((x))[[p]]. In particular it holds for Z_0 = Z, Z_1 = Z(e^beta o p).

---

## 5. The Toda relations at the dilaton point: Theorem M_k [P modulo (MV)]

No inductive hypothesis enters this section.

**Setting.** F_1 := log Z_1 = e^{beta Eu} F = Sum_{k>=0} beta^k Eu^k F / k!. Let Phi(eps) = F(b + eps) + F_1(b - eps) - F(b) - F_1(b). Since Z has constant term 1, (4.2) gives

    P_{k_0}(d/d eps) e^{Phi} |_{eps=0} = 0   for k_0 >= 2.

**Blocks.** Let alpha |- k_0 be a partition, viewed as a list of positions. For a monomial H^alpha, Faa di Bruno over set partitions pi of the positions (recall Phi(0) = 0) gives a sum of prod_{B in pi} d^B Phi|_0, with d^B = prod_{i in B} d/dt_i = prod_{i in B} (i x sigma Pi_i d/dy_i). And

    d^B Phi |_0 = Sum_{k>=0} c(|B|, k) (beta^k/k!) d^B Eu^k F,

with c(odd, 0) = 0, c(odd, k >= 1) = -1, c(even, 0) = 2, c(even, k >= 1) = 1.

**Coefficient of a block.** In the y-variables, Eu = Sum k y_k d/dy_k, which acts on y_{R'} G(y) as y_{R'} (J_{R'} + theta) G(y). Using (2.2), (2.4) and the third line of Lemma 3.1, the coefficient of y_{R_B}/|Aut R_B| (R_B with parts >= 2) in the (g,k) term of d^B beta^k/k! Eu^k hbar^{g-1} F_g is

    c(|B|,k) (prod_{i in B} i Pi_i) sigma^{|B|+g-1} (tau^k/k!) x^{n_B} (J_U + theta)^k G^(g)_U,
    U = R_B + B,     n_B := 2g - 2 + |B| + k.                                                     (5.0)

In detail: the derivatives d/dy_i, i >= 2, in d^B raise the profile from R_B to R' = R_B + {i in B : i >= 2}. With v = #{i in B : i = 1}, one has (d/dy)^v (J_{R'} + theta)^k G_{R'} = (J_{R'} + v + theta)^k G_{R'+1^v}. For R' empty, v = |B| >= 1, and (2.4) applies to the empty-profile part of F_g, with no constant term needed. The powers of x and sigma are x^{|B|} x^k x^{2g-2} and sigma^{|B|} sigma^{g-1}.

U is **non-empty**, because |B| >= 1. So Lemma 3.4 applies to every block. Neither the undifferentiated empty-profile function nor a singular v = 0 function ever occurs, and unit parts are treated like all other parts.

**Budget.** n_B >= 2g >= 0, because |B| + k >= 2 in every allowed case (odd |B| forces k >= 1, even |B| has |B| >= 2). The number of blocks is at most k_0. So at fixed x^n the sum over (g_B, k_B) is finite.

**Coefficient extraction.** The coefficient of y_R/|Aut R| in a product of blocks is the sum over labelled assignments of the parts of R to the blocks. So the coefficient of (x^n, y_R/|Aut R|) in the equation is a finite Q(tau)-combination of products of the series in (5.0), and it vanishes in Q(tau)[[y]].

**Single sector.** Substitute y = u e^{-u}. In every term, Sum_B J_{U_B} = J_R + k_0 and, by (3.2), Sum_B (w-order of B) = Sum_B (n_B + r_{R_B}) = n + r_R. So the coefficient is

    e^{u(J_R + k_0)} w^{n + r_R} * (polynomial in w),

which vanishes in Q(tau)[[u]]. Since e^{u(J_R+k_0)} is a unit and w is transcendental over Q(tau), every coefficient of that polynomial vanishes. In particular so does its constant term. This single-sector property holds because F contains all unstable genus-0 terms (l = 1, 2).

**Theorem M_k.** For all k_0 >= 2, n >= 0 and every multiset R of integers >= 2 (R empty allowed):

    Sum_{alpha |- k_0} c_alpha Sum_{pi} Sum_{R -> blocks} Sum_{(g_B, k_B) : Sum_B n_B = n} prod_{B in pi} phi_B = 0,     (5.1)

    phi_B = c(|B|, k_B) (prod_{i in B} i Pi_i) sigma^{|B|+g_B-1} (tau^{k_B}/k_B!) (delta_{g_B}(U_B))_{k_B} X_{g_B}(U_B),    U_B = R_B + B.

For k_0 = 2 this is, after multiplying by 2, R3's Theorem M: alpha = (2) gives the -2 Pi_2 line, the single block {1,1} gives the 2 sigma^{G+1} B_G term (k = 0) and the C-terms (k >= 1), using X_g(R+1+1) = J(J+1) X_g(R), and the two blocks {1},{1} give the Phi Phi products. Remark 3.6 reproduces R3's values Phi^0_1(empty) = tau, Phi^0_2(empty) = tau^2/2.

---

## 6. Proof of Theorem L [P modulo (MV)]

Induct on h >= 2 and, for fixed h, on d = 1, ..., 2h-3. Fix G = h and d, and assume:

- (H1) Theorem L in genera 2, ..., G-1;
- (H2) X_G(R') = 0 for all non-empty R' with parts >= 2 and D(R') < d.

No ordering inside a level is needed. By (S) these extend to multisets with parts >= 1:

- (H1') X_g(U) = 0 for 2 <= g < G, U non-empty, D(U) <= 2g-3;
- (H2') X_G(U) = 0 for U non-empty, D(U) < d.

Let V be a non-empty multiset with parts >= 2 and D(V) = d. Pick **any** part k_0 of V and put R = V - k_0, which may be empty; then D(R) = d - k_0 + 1. Take (5.1) with n = 2G. Only this even value of n is used.

**Count.** Since Sum_B D(U_B) = D(R) + k_0 - l(alpha), where l(alpha) = number of parts of alpha = Sum_B |B|, every term of (5.1) satisfies

    Sum_B (delta_{g_B}(U_B) - k_B) = D(R) + k_0 - Sum_B n_B = d + 1 - 2G <= -2 < 0.

So in every term some block has delta_B < k_B. Fix such a block, with data (g, k, U), and consider its factor phi_B.

- If 0 <= delta_B: then (delta_B)_k = 0, so phi_B = 0.
- If delta_B < 0: then D(U) <= 2g-3, and g >= 2 because delta_0(U) >= 2 and delta_1(U) >= 0 for non-empty U. Also g <= G, since n_B >= 2g and all n_B >= 0 sum to 2G.
  - If g < G: X_g(U) = 0 by (H1').
  - If g = G: n_B = 2G forces |B| + k = 2, and D(U) = D(R_B) + Sum_{i in B}(i-1) <= D(R) + k_0 - l(alpha) = d + 1 - l(alpha).
    - If l(alpha) >= 2: D(U) <= d - 1, so X_G(U) = 0 by (H2').
    - If l(alpha) = 1: alpha = (k_0), the only block is B = {k_0}, k = 1, R_B = R, so U = V. This is the target term.

So (5.1) reduces to its target term. By (4.3), c(1,1) = -1 and the block factor k_0 Pi_{k_0} sigma^G tau, it reads

    -(k_0 - 1) Pi_{k_0} sigma^G tau delta_G(V) X_G(V) = 0,      delta_G(V) = d + 2 - 2G <= -1.

The coefficient is a non-zero element of Q(tau), as k_0 >= 2. So X_G(V) = 0. This closes the induction. QED

**The base case by hand.** h = 2, d = 1, V = (2), R = empty, k_0 = 2, n = 4. (H1), (H2) are vacuous, and (H2') is Cor. 3.5.
- alpha = (2), c = 1/2, one block {2}, (g,k) in {(2,1), (1,3), (0,5)}: the last two die by (1)_3 = (3)_5 = 0, and the first gives Pi_2 sigma^2 tau X_2(2).
- alpha = (1,1), c = 1/2, one block {1,1}, (g,k) in {(2,0), (1,2), (0,4)}: X_2(1,1) = 0, (0)_2 = 0, (2)_4 = 0.
- alpha = (1,1), two blocks {1},{1}: each block has n_B = 2g - 1 + k. Genus 0 needs k <= 2 by (2)_k, so n_B <= 1; genus 1 dies by (0)_k, k >= 1; genus >= 2 dies by X_g(1) = 0. The total cannot reach 4.

So Pi_2 sigma^2 tau X_2(2) = 0.

**Uniformity.** Theorem M_k and Lemmas 3.1-3.4 are identities for every genus, inside and outside the range. Every hypothesis is used as an identity in Q(tau), and division by Pi_{k_0} (which vanishes at, for instance, tau = -1/2 when k_0 is even) occurs only in Q(tau). Specialisation comes only afterwards, in s.7.1.

---

## 7. Corollaries and interfaces

### 7.1 Corollary L1 [P modulo (MV)]

With Chern roots alpha_i of the Hodge bundle, (1 - alpha/tau)(1 + alpha/(1+tau)) = 1 - (alpha + alpha^2)/sigma, so

    Omega_h = Lambda(-1) Sum_{k=0}^h (-1/sigma)^k e_k(alpha + alpha^2).

e_k(alpha + alpha^2) is symmetric in the roots, hence a polynomial in the lambda classes, free of tau. So Xi_h(R;tau) = Sum_{k=0}^h (-1/sigma)^k Xi^(k)_h(R) with Xi^(k)_h(R) in Q for integer parts. As 1/sigma is transcendental over Q, Theorem L gives Xi^(k)_h(R) = 0 for each k. So the vanishing holds at every tau_0 in any field of characteristic 0 with sigma(tau_0) != 0, also where the proof's Pi_{k_0} vanishes. QED

### 7.2 Corollary L2 [P modulo (MV)] and the parts 1 and 0 [P]

Let U be non-empty with positive parts, and D(U) <= 2h-3. If U has only parts 1, X_h(U) = 0 by (S). Otherwise U = U' + 1^v with U' non-empty, parts >= 2, D(U') = D(U), and (S) with Theorem L gives the claim. The argument of s.7.1 applies to U: each Xi^(k)_h(U) = 0, so L2 also holds at every such tau_0.

For any multiset U of numbers with 2g-2+r > 0, and n = 2g-2+r:

- **Part 1.** G_{U+1} = d/dy G_U and Lemma 3.1 give Q_g(s;U+1) = (1+s)^2 d_s Q_g + [J + (1+s) n] Q_g, hence Xi_g(U+1) = J_U Xi_g(U).
- **Part 0.** Int Omega lam_U lam_1^a * 1 = (J+a) M_a by string, so G_{U+0} = (J + theta) G_U. Lemma 3.1 gives Q_g(s;U+0) = (J + s n) Q_g + s(1+s) d_s Q_g, hence Xi_g(U+0) = (J - n) Xi_g(U) = (D(U) - 2g + 2) Xi_g(U). (Note y d/dy = s(1+s)^2 d/ds.)

So part 0 is not covered by Theorem L: Xi_h(0) = (2-2h) Xi_h(empty). It is a lowering operator on D and reaches the gap.

### 7.3 Corollary L3: general scale [scale law P; vanishing P modulo (MV)]

For w = (E, w_2, w_3) with sum 0, tau := w_2/E gives w = E(1, tau, -1-tau), sigma = -w_2 w_3/E^2 != 0. Replace psi by psi/E. The degree-i part of Omega_w (normalised by w_a^h) is homogeneous of degree -i in w, so it scales as E^{-i}. The degree-b part of a leaf 1/(E - j psi) scales as E^{-1-b}. The degree-b part of -psi^2/(E - psi) scales as E^{1-b}. On Mbar_{h,r+m} the total degree is 3h-3+r+m, so the total exponent is -(3h-3+r+m) - r + m = -(3h-3+2r), independent of m. So Xi^E_h(U) = E^{-(3h-3+2r)} Xi_h(U; tau). QED

The normalisation by w_a^h matters: with the unnormalised product the exponent would be 3 - 2r. The vanishing is not affected. In this normalisation:

    Xi^E(U + 1) = (J_U/E^2) Xi^E(U),     Xi^E(U + 0) = ((D(U) - 2h + 2)/E^2) Xi^E(U),   lam^E_0 = 1/E.

These are (S_P2) and (L0_P2) of the local P2 draft.

The hypothesis that the **leaf scale E is one of the three weights** is essential ("mode B"). If w -> -w (mode A), the leaf scale is not a weight, and neither Theorem L nor its proof applies; the critics' controls at h = 2 confirm Xi(2) != 0 there [E].

### 7.4 Resolved conifold x C^2

The conifold convention (lemmaA/brief.md s.1) is w = (c, -x, -1), c = 1 + x, leaves 1/(c - j psi), T_c = -psi^2/(c - psi). So E = c, tau = -x/c, -1-tau = -1/c, sigma = -x/c^2. Cor. L3 gives Xi^(c)_h(R) = c^{-(3h-3+2r)} Xi_h(R; tau = -x/c) = 0 identically in x, and at every x with x c != 0. This is the Lemma A used by Theorem R and Corollary R of the conifold programme.

### 7.5 Local P2 x C^2

The local P2 convention (LP2 draft s.1.5, s.7) is w = (E, -eps1, -eps2), E = eps1 + eps2. So tau = -eps1/E, sigma = -eps1 eps2/E^2. Cor. L3 and Cor. L1 give Xi_h(R)(eps1, eps2) = 0 in Q(eps1, eps2) and at every complex (eps1, eps2) with eps1 eps2 E != 0, the diagonal eps1 = eps2 included. This proves Lemma A_P2(b) of the draft, including its extension to parts >= 1 (Cor. L2). The object is purely Hodge-side: the local P2 mirror curve and periods do not enter, so the ruling against target TR is not touched.

**What this document does not provide for either theory.**

- The identification of the Taylor coefficients of the graph sum with Xi_h(R). For local P2 this is written in LP2 draft s.1.5 (Fpot_h on the neighbourhood of q_0, dilaton-shifted form), s.7.1 (H0) and s.7.5(i)-(ii), tagged [P] there; the s.1.5 extension and (H0) had a pass in round 6 (D_rows, D_hodge). For the conifold it is leaf_structure_proof s.13.1(i), with the domain and (H0) in leaf_structure_proof s.13.0.
- The vertex identity for R empty (the gap): Theorem VI of the LP2 draft (s.7.6, proof s.7.8), equal to F4 at n = 0 for the conifold. It is a separate result with its own inputs (O5), and it uses nothing from Theorem L. (Referee D_VI's simplification of its Step 2 would use (4.2) at k = 2 of s.4 here, which does not depend on Theorem L, so no circularity arises.)
- The frame-B formula identifying the graph sum with the GW potential: Proposition 13.2_P2 (FRAMEB_PROOF.md; LP2 draft s.8.3, P modulo the literature and the draft statements listed there, among them eq:J=I (CCIT15 Thm 31) (T_max use confirmed by the user 2026-11-02, item O1 of FRAMEB_PROOF.md), assembled text refereed only as recorded there) and conifold Proposition 13.2 (leaf_structure_proof s.13.3).

**Range.** Graphs with an edge use Theorem L only for 2 <= h <= g - 1. The single-vertex graph also uses h = g.

---

## 8. Complements (proved, not used in s.3-s.6)

These state the whole matching of the MV series with the Hodge side. They live over Q(tau): the series below are literally undefined at zeros of Pi_k, and must not be evaluated there. Specialised tau is reached only through s.7.1.

### 8.1 All-profile, all-genus matching [P modulo (MV)]

In Q(tau)((x))[[s, P_2, P_3, ...]] put u = s/(1+s), x' = x/(1+s), p0_1 = s/(sigma x), p0_k = P_k/(sigma x Pi_k) (k >= 2), and

    FF(x,s;P) := log Sum_nu e^{-u|nu|} e^{(tau+1/2) kappa_nu x'/2} prod_box (2 sinh(h x'/2))^{-1} s_nu(p0).

**Proposition 8.1.** FF = Sum_{g>=0} Sum_R hbar^{g-1} QQ_g(s;R) P_R/|Aut R|, over multisets R of integers >= 2, where

    QQ_g(s;R) := (1+s)^{-(2g-2+r)} e^{-uJ} Sum'_{a>=0} M_g(R + 1^a) (u e^{-u})^a/a!,

and the prime omits (R,a) = (empty, 0). In particular only the powers x^{2g-2} occur, and sigma^{1-g} [x^{2g-2}] d_{P_{j_1}}...d_{P_{j_r}} FF |_{P=0} = QQ_g(s;R), with no Aut factor.

*Proof.* The map x' -> x/(1+s), p_1 -> e^{-u} s/(sigma x), p_k -> e^{-uk} P_k/(sigma x Pi_k) is a well-defined continuous ring map: p_1 lands in (s), p_k in (P), and each s-coefficient of x'^n = x^n (1+s)^{-n} is a monomial in x. Homogeneity gives s_nu(e^{-uk} p0_k) = e^{-u|nu|} s_nu(p0), so this map sends F to FF. Apply it to (2.1) at mu = R + 1^a, using 1/(1+s) = 1-u and s(1-u) = u. All x and sigma powers cancel except hbar^{g-1}. QED

### 8.2 Finite formula [P modulo (MV)]

For h >= 1, 2h-2+r > 0, R a non-empty multiset of integers >= 2, N = 3h-3+r, in Q(tau):

    Xi_h(R) = ( |Aut R| / (sigma^{h+r-1} Pi_R) ) Sum_{a=0}^{N} (-1)^a sigma^{-a} L_{N-a}^{(2h-2+r+a)}(-J-a) [x^{2h-2+r+a}] [p_R p_1^a] F.

This is (2.1) inserted into the Laguerre form of Lemma 3.2. For R empty and h >= 2, N = 3h-3: Xi_h(empty) = L_N^{(2h-2)}(0) <Omega_h> + sigma^{1-h} Sum_{a=1}^{N} (-1)^a sigma^{-a} L_{N-a}^{(2h-2+a)}(-a) [x^{2h-2+a}] [p_1^a] F, with <Omega_h> := Int_{Mbar_h} Omega_h (s.8.5); the a = 0 term does NOT carry the prefactor sigma^{1-h}.

### 8.3 Reformulation of Lemma A [P]

Theorem L is equivalent to: for each non-empty multiset R of integers >= 2, the series Sum_h hbar^{h-1} sigma^{1-h} [x^{2h-2}] d_{P_R} FF |_{P=0, s=-1}, taken coefficientwise as polynomials in s, is a Laurent polynomial in hbar of degree <= floor(D(R)/2). Indeed Xi_h(R) = 0 for all h with D(R) <= 2h-3 means that only h - 1 <= floor(D(R)/2) survive.

### 8.4 Table of QQ_g [P modulo (MV)]

Using Lagrange inversion for u = y e^u:

| (g, R) | QQ_g(s;R) | value at s = -1 |
|---|---|---|
| r >= 1, 2g-2+r > 0 | Q_g(s;R), polynomial of degree <= N | Xi_g(R) |
| g = 1, R = (j) | (j - 1 - 1/sigma + s)/24 | (j - 2 - 1/sigma)/24 |
| g = 0, r >= 3 | polynomial depending only on (J, r) | (J-1)...(J-r+3) |
| (0, (j,k)) | 1/(j+k) | 1/(j+k) |
| (0, (j)) | (1+s)/j^2 - s/(j(j+1)) | 1/(j(j+1)) |
| g >= 2, R empty | Q_g(s;empty) - (1+s)^{-(2g-2)} <Omega_g> | none: pole of order 2g-2 |
| (1, empty) | (1/24) [ log(1+s) + (1+1/sigma)(1/(1+s) - 1) ] | none: log and simple pole |
| (0, empty) | s + s^2/4 + s^3/(6(1+s)) | none: simple pole |

So QQ_g has a value at s = -1 for every non-empty R, and for R empty it has none in any genus. For R empty and g >= 2, the value Xi_g(empty) is recovered only after restoring the a = 0 term, as in s.8.2.

Derivations. (0,(j,k)): M_0 = (J+a)^{a-1}, so G = e^{Ju}/J. (0,(j)): Lemma 3.3. (0,empty): M_0 = a^{a-3}, (y d/dy)^2 G = u, G(0) = 0, so G = u - 3u^2/4 + u^3/6. (1,(j)): Omega_1 = 1 - (1+1/sigma) lambda_1 on Mbar_{1,n}, as -1/tau + 1/(1+tau) = -1/sigma and lambda_1^2 = 0, with <tau_1>_1 = <lambda_1>_1 = <tau_0 tau_2>_1 = 1/24; k >= 2 is excluded by degree. (1,empty): y d/dy G = y G_(1) = s Q_1(s;(1)) = s (s - 1/sigma)/24 by Lemma 3.2, and y d/dy = s(1+s)^2 d/ds, G(s=0) = 0. (0, r >= 3): A_m(X) = Sum_a (X+a)^{a+m} y^a/a! satisfies (X + theta) A_{m-1} = A_m, which gives Q^[m] = X Q^[m-1] + m s Q^[m-1] + s(1+s) d_s Q^[m-1] with Q^[-1] = 1/X, so Q^[m](-1) = (X-m) Q^[m-1](-1).

### 8.5 The empty profile (not part of Lemma A)

- **<Omega_g> from the MV series alone.** [P modulo (MV)] Lemma 3.2 at r = 0, J = 0 and the degree bound at m = N+1 give Sum_{a<=m} (M_a/a!) (-1)^{m-a} L_{m-a}^{(2h-3+a)}(-a) = 0. The a = 0 coefficient is (-1)^{N+1} L_{N+1}^{(2h-3)}(0) = (-1)^{N+1} binom(N+2h-2, N+1) != 0. So <Omega_h>, the a = 0 value of Lemma 3.2 at U = empty, is a finite combination of MV coefficients M_a, a >= 1.
- **Closed form.** [P modulo C: Mumford, Faber-Pandharipande] <Omega_g> = -A_g sigma^{1-g}, A_g = Int lambda_{g-2} lambda_{g-1} lambda_g. On w_1+w_2+w_3 = 0 the degree 3g-3 part of prod Lambda^vee(w_a) integrates to (-1)^{g-1} [ e_3(w) Int lambda_{g-1}^3 + (Sum_{a != b} w_a^2 w_b) A_g ]: the (3,0,0) terms die by lambda_g^2 = 0. Mumford gives lambda_{g-1}^2 = 2 lambda_{g-2} lambda_g, so Int lambda_{g-1}^3 = 2 A_g, and Sum_{a != b} w_a^2 w_b = -3 e_3(w). Hence Int prod Lambda^vee(w_a) = (-1)^g A_g e_3(w); at w = (1, tau, -1-tau), e_3 = -sigma, this gives the stated value. The closed form of A_g (as quoted in OP s.0.6) is not needed here.
- **Xi_h(empty).** Xi_h(empty) = Q_h(-1;empty) is the gap. Closed form: [C: LP2 draft Theorem VI_tau, s.7.8, TR-free] sigma^{h-1} Xi_h(empty) = beta_{2h-2}(-tau, 1+tau), beta_n(a,b) = (n-1)! [t^n] 1/(4 sinh(at/2) sinh(bt/2)); this is R2 Cor. 7 (Hodge-side TR, cross-check only). Non-vanishing: with S(y) = (y/2)/sinh(y/2) = Sum c_i y^{2i}, sign c_i = (-1)^i, c_i != 0, one has beta_{2h-2}(a,b) = ((2h-3)!/(ab)) Sum_{i+j=h} c_i c_j a^{2i} b^{2j}, whose terms all have sign (-1)^h at real a, b != 0. So Xi_h(empty) != 0 at every real tau with sigma != 0, hence in Q(tau). This non-vanishing, **unlike Theorem L, does not specialise**: at h = 2 one has Xi_2(empty) = -(24 + 28 s' + 7 s'^2)/5760 with s' = 1/sigma, which has zeros at sigma = -0.804..., -0.363...; both satisfy sigma < -1/4, so tau is complex there, in agreement with the real non-vanishing. Local P2 statements allowing complex weights must take this into account.
- **Part 0.** Xi_h(U + 0) = (D(U) - 2h + 2) Xi_h(U) is s.7.2. This identity needs no citation.

### 8.6 Variants, not used

- Proof (i) of R3: P_2 alone (2-insertion, Theorem M), then Route A Theorem B (MV cut-and-join) or Route A Cor. E (KP). More inputs and no gain.
- R2 (lemmaA/round2/R2_onepoint): ABDKS x-y swap plus Eynard 2014 on the Hodge curve C_tau, non-inductive. Allowed as Hodge-side TR. It stays a backup. Its literature and its normalisation "fixed against TR-free data for g <= 9" are unpinned, so it does not meet the all-orders ruling as written (critic 1 S3(b) of round 4).
- R1: not admissible under the all-orders ruling, since its base h <= 7 and zero count are computer-assisted.

---

## 9. Literature: what was and was not confirmed

| item | used for | status | confirmed | not confirmed |
|---|---|---|---|---|
| Liu-Liu-Zhou, J. Differential Geom. 65 (2003) 289-340; arXiv:math/0306434v3 | **the input** | [C] | Introduction read directly in arXiv v3 (author and referees of round 5) and in the JDG print (round 6, EMIS mirror, https://www.maths.tcd.ie/EMIS/journals/NYJM/jdg/p/2003/65-2-4.pdf), same numbering: C_{g,mu} with -sqrt(-1)^{|mu|+l}/|Aut mu| and the l < 3 extension (JDG pp. 290-291); eqs. (1), (2), (3) (p. 291); "we will call (1) the Marino-Vafa formula", Theorem 1 with (4), (5) (p. 292); Theorem 2 with (6) (p. 293); (8), (9), Prop. 2.1 (p. 294); kappa_nu. Theorems 1, 2 are the cut-and-join for R and C, so "LLZ Thm 1" for the formula is wrong; "eq. (1)" is confirmed in JDG and arXiv v3. Prop. 2.2 (tau = 0) agrees. Journal data confirmed on Project Euclid. | LLZ's proof not re-checked. |
| Okounkov-Pandharipande, Geom. Topol. 8 (2004) 675-699; arXiv:math/0307209v3 | second derivation of (2.1); independent published proof of the formula | [C], not load-bearing | Theorem 1 (0.8), (0.1), (0.5), s.0.6, (1.24), (1.25), (1.28) read directly | The printed 1/|Aut mu| in (0.5), (0.8) is inconsistent at repeated parts (s.2.2) |
| String and dilaton with pulled-back classes | Lemma 3.2 | textbook | Mechanism stated in Lemma 3.2 | No theorem number pinned |
| Genus-0 multinomial formula | M_0, s.1 | textbook | Stated in s.1 | No number pinned |
| Row-hook lemma | Prop 2.1(a) | textbook, sketch in place | LLZ (8)-(9) point to Macdonald pp. 10-11 | Macdonald I.1 Ex. 1 from memory |
| Principal specialisation of Schur functions | Prop 2.1(b) | textbook | Exponent computation done in place; LLZ eq. (12) (JDG p. 297) states s_nu(1,q,q^2,...) = q^{n(nu)}/prod(1-q^{h}) | Macdonald I.3 Ex. 2 and Stanley EC2 Cor. 7.21.3, both from memory |
| Frobenius; Murnaghan-Nakayama | s.4 (F1) | textbook | | No number pinned |
| Lagrange inversion | Lemma 3.3, s.8.4 | textbook | | |
| Jimbo-Miwa, Publ. RIMS 19 (1983) 943-1001 | corroboration of (4.2), not an input | not an input | Read in scan by the ThmL author: Thm 1.1-1.2; Appendix 1 p. 996, table "1st Modified KP hierarchy". Exact linear algebra: P_2..P_5 lie in the span of the transcribed entries | Equation number of the bilinear identity not looked up; no referee re-read the table |
| Harnad-Orlov, CMP 338 (2015); arXiv:1407.7800v9 | corroboration only | not an input | (2.2), (2.8)-(2.10), (5.1)-(5.4) read. (2.2) as printed has exponent e^{-xi(dt,z)}, which appears to be a sign slip: it fails exactly at level >= 2 on a generic content-product tau, while e^{+xi} holds and coincides with (4.1) | Journal version not compared |
| Okounkov, Math. Res. Lett. 7 (2000) 447-453 | no longer needed | not an input | Abstract read | "s.2.7" of R3 unconfirmed |
| Orlov-Scherbin TMP 128 (2001) 906-926; Ueno-Takasaki ASPM 4 (1984) 1-95; Takasaki ASPM 4 (1984) 139-163; Takebe LMP 21 (1991) 77-84 | no longer needed | not an input | Bibliographic data from Harnad-Orlov's list only | Texts and theorem numbers not read |
| Mumford relation; Faber-Pandharipande A_g | s.8.5 only | [C], outside Lemma A | A_g closed form as quoted in OP s.0.6. FP, arXiv:math/9810173v1, read: Theorem 4 and s.4.3 give Int lambda_{g-1}^3 = (|B_{2g}|/2g)(|B_{2g-2}|/(2g-2))/(2g-2)! for g >= 2, and that this equals 2 Int lambda_g lambda_{g-1} lambda_{g-2} from c(E)c(E^vee) = 1. Agrees with OP s.0.6 and s.8.5 | Journal numbering of FP not compared |
| Arbarello-Cornalba pushforward | only to rewrite the m-sum as K = exp(-Sum kappa_k/k); the m-sum is the definition of record here | textbook, not used | Rederived: 1 - s w/(1-w) = (1 - (1+s) w)/(1-w) gives K_s = exp(Sum ((1+s)^k - 1) kappa_k/k) | |

Web sources read in round 5: arXiv:math/0306434, projecteuclid.org/euclid.jdg/1090511689, arXiv:math/0307209, arXiv:1407.7800. In round 6: the JDG print of LLZ (EMIS mirror above), arXiv:math/9810173v1.

---

## 10. Status table and open points

### 10.1 Status

| result | status | rests on |
|---|---|---|
| Prop 2.1, (2.1) real form | P mod C | LLZ (1)-(3), textbook specialisation; all signs derived |
| (2.2) | P mod (MV) | restates (2.1); Pi_k != 0 in Q(tau) |
| (2.3), (2.4) | P | definition and formal |
| Lemma 3.1, Lemma 3.4, Cor 3.5, Remark 3.6 | P | Lemmas 3.2, 3.3 |
| Lemma 3.2 | P mod textbook | string, dilaton |
| Lemma 3.3 | P mod textbook | Lagrange inversion |
| (F1)-(F6), P_k for all k, (4.2), (4.3) | P mod textbook | Murnaghan-Nakayama, Frobenius |
| Theorem M_k (5.1), all k_0, n, R | P mod (MV) | s.2-s.4 |
| Theorem L | P mod (MV) | Theorem M_k, (S); no hypothesis beyond (H1), (H2) |
| Cor L1, L2, L3 | P mod (MV) | Theorem L; the scale law of L3 and the identities for parts 1, 0 of s.7.2 are P |
| s.7.4, s.7.5 transfers of the vanishing | P mod (MV) | Cor L1-L3 |
| Prop 8.1, s.8.2, s.8.3, s.8.4 | P mod (MV) | not used in Theorem L |
| s.8.5, closed form of <Omega_g> | P mod C | Mumford; outside Lemma A |
| s.8.5, closed form and non-vanishing of Xi_h(empty) | C | LP2 draft Theorem VI_tau (s.7.8); outside Lemma A |

**Precise claim.** Theorem L, for non-empty R, is proved to all orders, modulo the single cited theorem (MV) and the textbook facts of s.0. No step of s.2-s.6 is open.

### 10.2 Open points

- **O1 (refereeing): closed in round 6.** The assembled text (frozen copy LEMMA_A_PROOF_round5.md) was refereed under two lenses. L_logic: s.2-s.7 line by line by hand, including every sign, (F1)-(F6) and the induction of s.6, a search for hidden finite-order inputs and circularity, and the literature of s.9 read directly. L_comp: new code written from the document's formulas, plus hand checks of the same steps: (4.2), (5.1), the logic of s.6 term by term, Theorem L and its polynomial form in 1/sigma; these runs are evidence only. Both verdicts: sound with minor fixes, no major item. The confirmed fixes are in s.11.4. Their stated limits: the proof of (MV) was not re-checked (O3); s.8.6 and the runs of s.12 were not repeated; the interfaces listed in s.7.5 were not assessed. The referees had not seen this edited text, and the edits change no mathematics of s.2-s.6.
- **O2 (citation details; none is mathematical).** Macdonald and Stanley exercise numbers, from memory. No pinned theorem number for string/dilaton with pulled-back classes, for Frobenius or for Murnaghan-Nakayama. Harnad-Orlov's sign in the journal version.
- **O3 (input not re-proved).** The proof of (MV) in LLZ was not re-checked. (MV) has a second published proof (OP), with the caveat on the printed 1/|Aut mu|.
- **O4 (editing not done, outside this task).** The other files still carry the items of s.11.3.
- **O5 (outside Lemma A, unchanged by this document).**
  - Empty profile: the gap formula and (VI) are Theorem VI of the LP2 draft (s.7.8), TR-free, P modulo literature. In the round-5 draft that literature is LLZ eqs. (1), (5) and the vertex-operator form of psi, psi*; referee D_VI of round 6 (sound with minor fixes, no major item; second pass on Theorem VI) gives replacements that reduce it to LLZ eq. (1) at mu = 1^a plus textbook facts and (4.2) at k = 2 here. The draft's s.7.8 is the record: it uses LLZ eq. (1) only, at mu = 1^a, and neither eq. (5) nor the vertex-operator form. This closes plan P5 for j = 0, and critic 1 S3 of round 4 is answered there. R2 Cor. 7 is a cross-check only (draft s.7.6 (B-a)).
  - The identification of Theorem R's Taylor coefficients with Xi_h(R): see s.7.5 (local P2 written in the draft, conifold: leaf_structure_proof s.13.0, s.13.1).
  - lem:modrel (critic 1 S4) is closed by Proposition 2.12 of the LP2 draft. Prop. 13.2_P2 and the local P2 gap corollary are now in FRAMEB_PROOF.md (P modulo the literature and the draft statements listed there, among them eq:J=I (CCIT15 Thm 31) (T_max use confirmed by the user, item O1 of that file, which is not the O1 of s.10.2 here); assembled text refereed only as recorded there); the conifold Prop. 13.2 is leaf_structure_proof s.13.3.

---

## 11. Disposition of critic and referee items

### 11.1 Items now discharged

| item | how |
|---|---|
| Round 4 critic 1 S1 (MV-Xi matching only "normalisation exact" with data comparisons) | Prop 2.1 and Lemma 3.2, all genera, all profiles; no constant from data |
| Round 4 critic 1 S2 (no proof of record; P_k fixed by checks for k <= 4; input list incomplete; P1 missing) | Record = R3 proof (ii) rewritten here; s.4 derives P_k for all k; input list in s.0; this file is P1 |
| Round 4 critic 2 minor on s.7.4, s.8.2 (input list, proof of record) | s.0 |
| Round 4 lemmaA referee: MV matching belongs among open steps | Closed by s.2 |
| Round 4 lemmaA referee: hypotheses never used at specialised tau; direct source of the Laurent form; no part 0 | s.6 Uniformity; s.7.1; Cor L2 |
| round2_status s.5.2 item 1 (non-empty R) | In the statement |
| round2_status s.1.1 (a) (additive constant of G_empty) | s.2.3: the a = 0 term is explicit; R3 Lemma 2.1 needs the correction, which is not a new finding |
| round2_status s.1.1 (b) (v = 0 functions) | Moot: every block has U non-empty |
| round2_status s.1.1 (c) (single sector) | s.5, by a purely formal argument |
| round2_status s.1.1 (d) (odd n sampling) | (5.1) is proved for all n; only n = 2G is used; no sampling claim |
| round2_status s.5.2 item 2 (P_k, I_2 wording, dependencies) | s.4; the I_2 wording does not arise; proof (i) dependencies are not used |
| R3 Lemma 2.3 (empty-profile table, "referee should re-derive") | Remark 3.6, corollary of Lemma 3.4, not needed separately |

Item 2 of the round-4 draft's s.8.2 (backup LP2_REGULARITY_PROOF_draft_round4.md): its first bullet (2D Toda input) is replaced by s.4, its second bullet (LLZ/OP numbers) is settled against the JDG print and arXiv v3 (s.9), and its third bullet (matching) is closed by s.2.

### 11.2 Round-5 referee fixes applied here

From the MVmatch referees:
- False sentence on poles of QQ_g corrected: s.8.4 says that every empty-profile QQ_g lacks a value at s = -1, with pole order 2g-2 at g >= 2.
- Quantifier of the reformulation fixed: s.8.3.
- Overclaim of "discharged" removed: s.11.1 and O5 state exactly what is and is not closed.
- LLZ Prop. 2.1 not cited for the hook form: Prop 2.1(a) uses (LLZ-3) and its own derivation.
- <Omega_h> shown to be determined by MV data: s.8.5.
- E3 given as a cross-reference to round2_status s.1.1(a).
- Citation limits recorded: s.9.
- Scope note: Prop 8.1 and s.8.2 hold over Q(tau) only, stated at the head of s.8.

From the ThmL referees:
- Non-vanishing of Xi_h(empty) moved out of the theorem to s.8.5, with its source and with the warning that it does not specialise (exact h = 2 formula, complex zeros).
- The (2-2h) identity proved: s.7.2.
- The vacuous sub-case "T_B != R" of the induction deleted: s.6.
- Definition of G^(g)_U and the empty-profile part of F_g made explicit in all genera: s.2.3.
- Principal specialisation added to the textbook list, and the hook identity written in, not only cited.
- Notation cleaned: T no longer means both a multiset and a leaf class, D no longer means level, y d/dy and Hirota derivative, u no longer means both tree variable and Hirota shift, tau no longer means both framing and tau function.
- X_0 flagged as not a Xi: s.1.
- (S) stated for g >= 2, the only use.
- Harnad-Orlov worded as "appears to be" a slip: s.9.
- The label "P mod (MV)" restricted to the vanishing, and not extended to the Taylor identification: s.7.5, O5.

### 11.3 Errata and edits to make in other files (not done here)

- **E1.** lemmaA/A_mv/report.md line 177 misquotes LLZ: prefactor -(sqrt(-1))^l instead of -(sqrt(-1))^{|mu|+l}, and V_nu in [m] = q^{m/2} - q^{-m/2} instead of the printed 2 sin. The two slips cancel only together, because prod 1/[h] = (-i)^{|nu|} V_nu is exactly the rescaling p_k -> (-i)^k p_k. Its (1.1) is correct.
- **E1'.** "LLZ Thm 1" should read "LLZ eq. (1)": A_mv/report.md line 8, R3 report line 250, round2_status lines 26, 42 and 460.
- **E2.** Route A Cor. 2.3: "+ (unstable constant)" is not constant in s at fixed x. It is -Sum_{g>=2} A_g x'^{2g-2} plus the genus-0, 1 empty terms of s.8.4. "Xi is the genus-by-genus value at s = -1" fails for R empty in every genus: at g = 0, 1 there is no value, and at g >= 2 it needs the a = 0 term restored. "P_1 := s" means that parts equal to 1 are not separate variables.
- **E3.** R3 Lemma 2.1: G^(g)_0 lies in w^{2g-2}Q[w] only after adding the constant <Omega_g> (round2_status s.1.1(a)). Harmless there.
- **E4.** OP (0.5), (0.8): printed 1/|Aut mu|, s.2.2.
- **E5.** LLZ Prop. 2.1 as printed has 2^l (JDG p. 294) where (LLZ-3) gives one factor 2 per box.
- **E6.** R3 report lines 11, 34, 246: statement lacks "non-empty".
- **LP2 draft.** s.0.3, s.7.4(b), s.8.2 were edited in round 5, and "Nothing else" is gone. Remaining: s.8.1 label of Lemma A_P2(b) to "P mod (MV), LEMMA_A_PROOF.md", and the statements in s.0.3, s.7.4 and s.8.2 that this file has had no referee, now that O1 is closed; the last sentence of the s.7.4 bullet "Not provided by that file" (see O5); "JDG numbering not confirmed" in s.7.4, s.7.8, s.8.2 (see s.9); "LLZ eq. (13) states it" in s.7.8 should read eq. (12).
- **Conifold notes.** leaf_structure_proof s.0.4, s.13, G1; proof_status; conifold_regularity_note; proof_M_residue (round2_status s.5.5); RESUME.md and memory.

### 11.4 Round-6 referee fixes applied here

From L_logic and L_comp (scripts: scratchpad lp2_round6/L_logic/, lp2_round6/L_comp/; no mathematics of s.2-s.6 changed):
- Stale passages synced with the LP2 draft: s.7.5, s.8.5, s.10.1, O5, s.11.1-s.11.3.
- Xi_h(empty): closed form quoted from Theorem VI_tau, and a written proof of non-vanishing in Q(tau) from the signs of the coefficients of S: s.8.5.
- Letters in s.6: multiset V, level d, l(alpha), so that T, D and |.| keep one meaning each (s.11.2 claimed this before it was true); "m-th term" in Lemma 3.2.
- Status labels: L1, L2, the vanishing in L3 and (2.2) are P mod (MV), not P: s.0, s.7.1-7.3, s.10.1.
- Scope of specialisation: any field of characteristic 0, and parts >= 1: s.0, s.7.1, s.7.2.
- Sign (-1)^{N+1} of the a = 0 coefficient, and <Omega_h> in place of M_0: s.8.5. Placement of the prefactor sigma^{1-h} for R empty: s.8.2.
- Citations: JDG print of LLZ read, LLZ eq. (12), FP Theorem 4 and s.4.3, E5: s.9, O2, s.11.3.

---

## 12. Evidence [E] (no proof step uses any of this)

All exact unless stated. Paths are under <user home>/AppData/Local/Temp/claude/G--My-Drive-Projects-Old-refGW/37f24460-3601-4c85-b123-5156e19eb475/scratchpad/lp2_round5/ (session scratchpad; copy them into the project if they are to be kept).

| check | range | result | path |
|---|---|---|---|
| Fresh check of this document's formulas: Prop 2.1 Schur side, Laguerre form, (5.1), target coefficient | tau = 2/7, -5/3; G = 2, 3; k_0 <= 4; n = 2G, 2G+1; 48 relations each | 0 failures; target response exact 24/24; flipping c(odd,k) or beta -> -beta breaks 13/24 (the rest are terms that vanish anyway). X_2(3) = 1367/186624, X_3(5) = 205027/47029248 at tau = 2/7, equal to the reports' values | LEMMA_A/doc_check.py, doc_check_t27.out, doc_check_tm53.out |
| (2.1) against an independent Hodge engine (DVV, Mumford), all new code | author: g <= 4, 150 pairs (g, mu), 1387 comparisons; referee 2: g <= 5 at l <= 2, g <= 4 at l = 3, 595 comparisons, several tau including zeros of Pi | 0 mismatches; 9 kinds of control all fail | MVmatch/, MVmatch_ref2/ |
| LLZ (1)-(3) literally in lambda over Q(i) | 90 cases | match; the mixed misquote fails 66/90 | MVmatch/check3.py |
| Prop 8.1 directly from the Schur sum | g <= 3 (author), g <= 6 empty profile (referee) | 0 mismatches; x' = x, no e^{-u|nu|}, s -> -s, no Pi all fail | MVmatch/, MVmatch_ref2/testC*.py |
| Lemma 3.2 and operators | genus 0-4, non-integer leaves | 0 failures | MVmatch/, MVmatch_ref2/testB*.py |
| s.8.2 against the definition, h = 2, 3, 4; Schur-side zeros at genus 4-7 | Xi_4(7) = 193375/214990848 at tau = 2/7 non-zero, as expected out of range | 0 differences | MVmatch/check4.py, MVmatch_ref2/testD*.py |
| P_2..P_7 on generic content-product tau; P_k <= 10 at low weight | degree 8, 9 | 0 non-zero; swapped order, tau_{-1}, wrong c_(k) non-zero | ThmL/check1_toda.py, ThmL_ref2/t6_pk.py |
| P_k on Z_MV as Laurent series in x | tau = 2/7, -5/3; k <= 6, weight <= 8, order 30 | 114 coefficients, 0 failures | ThmL/check6_mvbilinear.py |
| (5.1), MV-side values | G <= 4 (author, 1638 relations), G <= 6 mod 2^61-1 (referee 1, 218 + 196 relations) | 0 failures | ThmL/check3_relation.py, ThmL_ref1/ref_m41.py |
| (5.1), Hodge-side values, MV not assumed | G <= 5, including tau = -1/3, a zero of Pi_3 | 186 relations, 0 failures; 5 kinds of control break it | ThmL_ref2/t2_*.py, t3_run.py |
| Target coefficient for every part k_0, including repeated | | 202/202, 124/124, 32/32 | ThmL/, ThmL_ref2/t7_resp.py |
| Mechanical run of the s.6 logic over all terms | all V for G <= 5 | 53.4 million terms, 0 survivors, one target per case | ThmL/check4_logic.py |
| Theorem L directly on the Hodge side | h = 3 (r <= 4), h = 4 (r <= 5, all in-range V), h = 5 (r <= 3), h = 6 (r <= 2), at tau = 2/7, -5/3, -1/2, 13/9 | 383 tests, all zero; just outside the range all 228 controls non-zero | ThmL_ref2/t1_thmL.py |
| Earlier rounds | genus-7 engines (194 zeros), R3 referees' mechanical logic G <= 12, Mathematica q-tables | see round2_status s.1 | lemmaA/round2/ |

Limits: the g = 6 direct Hodge comparison of (2.1) did not finish. The G = 6-8 mechanical logic run with r <= 2 did not finish. Neither is claimed.
