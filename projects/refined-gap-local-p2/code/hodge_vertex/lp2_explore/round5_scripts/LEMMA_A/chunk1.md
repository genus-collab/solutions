# Lemma A: all-genus proof of record

File: hodge_vertex_conifold/lemmaA/LEMMA_A_PROOF.md. Written 2026-11-02 (round 5, plan item P1 of lemmaA/round2/round2_status.md s.6).

Purpose: the single written proof of Lemma A = M_irr(h) for non-empty R, for the weights (1, tau, -1-tau), usable by both the resolved conifold x C^2 proof and the local P2 x C^2 proof.

Tags: **[P]** proved in this document. **[C]** cited; s.9 says what was and was not confirmed. **[E]** evidence, used by no proof step. **[O]** open.

Provenance: the proof combines the round-5 reports MVmatch (s.2, s.8) and ThmL (s.3-s.6), each with two referee reports (verdict "sound with minor fixes", no major item). All referee fixes are applied; s.11 lists them. The writer re-derived s.2-s.6 line by line during assembly. The assembled text has not yet had its own referee (open point O1, s.10).

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

**Corollary L1 (specialisation).** [P] Write Lambda(t) = Sum_i t^i lambda_i, and let alpha_i be the Chern roots of the Hodge bundle. Then Xi_h(R; tau) = Sum_{k=0}^h (-1/sigma)^k Xi^(k)_h(R), with Xi^(k)_h(R) = Sum_m (1/m!) Int Lambda(-1) e_k(alpha + alpha^2) lam_R T^m independent of tau. For R as in Theorem L, every Xi^(k)_h(R) vanishes. Hence Xi_h(R; tau_0) = 0 at every complex tau_0 with sigma(tau_0) != 0, including the diagonal tau_0 = -1/2 and all other zeros of the factors Pi_k used in the proof.

**Corollary L2 (parts >= 1).** [P] The same vanishing holds for every non-empty multiset U of positive integers with D(U) <= 2h - 3. No part 0 is allowed.

**Corollary L3 (general scale).** [P] Let E, w_2, w_3 lie in a field of characteristic 0, with E + w_2 + w_3 = 0 and E w_2 w_3 != 0. Put Omega_w = prod_a Lambda^vee(w_a)/w_a^h for w = (E, w_2, w_3), lam^E_j = 1/(E - j psi), T^E = -psi^2/(E - psi). The corresponding integrals satisfy Xi^E_h(U) = E^{-(3h-3+2r)} Xi_h(U; tau = w_2/E), hence vanish in the same range. This covers the conifold and local P2 conventions (s.7.4, s.7.5).

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

