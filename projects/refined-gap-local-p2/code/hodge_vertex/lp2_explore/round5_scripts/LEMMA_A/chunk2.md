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

*Proof.* Write Tp = lam_1 - 1 - psi and expand the k-th term of Q_g:

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

Induct on h >= 2 and, for fixed h, on D = 1, ..., 2h-3. Fix G = h and D, and assume:

- (H1) Theorem L in genera 2, ..., G-1;
- (H2) X_G(R') = 0 for all non-empty R' with parts >= 2 and D(R') < D.

No ordering inside a level is needed. By (S) these extend to multisets with parts >= 1:

- (H1') X_g(U) = 0 for 2 <= g < G, U non-empty, D(U) <= 2g-3;
- (H2') X_G(U) = 0 for U non-empty, D(U) < D.

Let T be a non-empty multiset with parts >= 2 and D(T) = D. Pick **any** part k_0 of T and put R = T - k_0, which may be empty; then D(R) = D - k_0 + 1. Take (5.1) with n = 2G. Only this even value of n is used.

**Count.** Since Sum_B D(U_B) = D(R) + k_0 - |alpha|, where |alpha| = number of parts of alpha = Sum_B |B|, every term of (5.1) satisfies

    Sum_B (delta_{g_B}(U_B) - k_B) = D(R) + k_0 - Sum_B n_B = D + 1 - 2G <= -2 < 0.

So in every term some block has delta_B < k_B. Fix such a block, with data (g, k, U), and consider its factor phi_B.

- If 0 <= delta_B: then (delta_B)_k = 0, so phi_B = 0.
- If delta_B < 0: then D(U) <= 2g-3, and g >= 2 because delta_0(U) >= 2 and delta_1(U) >= 0 for non-empty U. Also g <= G, since n_B >= 2g and all n_B >= 0 sum to 2G.
  - If g < G: X_g(U) = 0 by (H1').
  - If g = G: n_B = 2G forces |B| + k = 2, and D(U) = D(R_B) + Sum_{i in B}(i-1) <= D(R) + k_0 - |alpha| = D + 1 - |alpha|.
    - If |alpha| >= 2: D(U) <= D - 1, so X_G(U) = 0 by (H2').
    - If |alpha| = 1: alpha = (k_0), the only block is B = {k_0}, k = 1, R_B = R, so U = T. This is the target term.

So (5.1) reduces to its target term. By (4.3), c(1,1) = -1 and the block factor k_0 Pi_{k_0} sigma^G tau, it reads

    -(k_0 - 1) Pi_{k_0} sigma^G tau delta_G(T) X_G(T) = 0,      delta_G(T) = D + 2 - 2G <= -1.

The coefficient is a non-zero element of Q(tau), as k_0 >= 2. So X_G(T) = 0. This closes the induction. QED

**The base case by hand.** h = 2, D = 1, T = (2), R = empty, k_0 = 2, n = 4. (H1), (H2) are vacuous, and (H2') is Cor. 3.5.
- alpha = (2), c = 1/2, one block {2}, (g,k) in {(2,1), (1,3), (0,5)}: the last two die by (1)_3 = (3)_5 = 0, and the first gives Pi_2 sigma^2 tau X_2(2).
- alpha = (1,1), c = 1/2, one block {1,1}, (g,k) in {(2,0), (1,2), (0,4)}: X_2(1,1) = 0, (0)_2 = 0, (2)_4 = 0.
- alpha = (1,1), two blocks {1},{1}: each block has n_B = 2g - 1 + k. Genus 0 needs k <= 2 by (2)_k, so n_B <= 1; genus 1 dies by (0)_k, k >= 1; genus >= 2 dies by X_g(1) = 0. The total cannot reach 4.

So Pi_2 sigma^2 tau X_2(2) = 0.

**Uniformity.** Theorem M_k and Lemmas 3.1-3.4 are identities for every genus, inside and outside the range. Every hypothesis is used as an identity in Q(tau), and division by Pi_{k_0} (which vanishes at, for instance, tau = -1/2 when k_0 is even) occurs only in Q(tau). Specialisation comes only afterwards, in s.7.1.

---

## 7. Corollaries and interfaces

### 7.1 Corollary L1 [P]

With Chern roots alpha_i of the Hodge bundle, (1 - alpha/tau)(1 + alpha/(1+tau)) = 1 - (alpha + alpha^2)/sigma, so

    Omega_h = Lambda(-1) Sum_{k=0}^h (-1/sigma)^k e_k(alpha + alpha^2).

e_k(alpha + alpha^2) is symmetric in the roots, hence a polynomial in the lambda classes, free of tau. So Xi_h(R;tau) = Sum_{k=0}^h (-1/sigma)^k Xi^(k)_h(R) with Xi^(k)_h(R) in Q for integer parts. As 1/sigma is transcendental over Q, Theorem L gives Xi^(k)_h(R) = 0 for each k. So the vanishing holds at every tau_0 with sigma(tau_0) != 0, also where the proof's Pi_{k_0} vanishes. QED

### 7.2 Corollary L2 and the parts 1 and 0 [P]

Let U be non-empty with positive parts, and D(U) <= 2h-3. If U has only parts 1, X_h(U) = 0 by (S). Otherwise U = U' + 1^v with U' non-empty, parts >= 2, D(U') = D(U), and (S) with Theorem L gives the claim.

For any multiset U of numbers with 2g-2+r > 0, and n = 2g-2+r:

- **Part 1.** G_{U+1} = d/dy G_U and Lemma 3.1 give Q_g(s;U+1) = (1+s)^2 d_s Q_g + [J + (1+s) n] Q_g, hence Xi_g(U+1) = J_U Xi_g(U).
- **Part 0.** Int Omega lam_U lam_1^a * 1 = (J+a) M_a by string, so G_{U+0} = (J + theta) G_U. Lemma 3.1 gives Q_g(s;U+0) = (J + s n) Q_g + s(1+s) d_s Q_g, hence Xi_g(U+0) = (J - n) Xi_g(U) = (D(U) - 2g + 2) Xi_g(U). (Note y d/dy = s(1+s)^2 d/ds.)

So part 0 is not covered by Theorem L: Xi_h(0) = (2-2h) Xi_h(empty). It is a lowering operator on D and reaches the gap.

### 7.3 Corollary L3: general scale [P]

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

- The identification of the Taylor coefficients of the graph sum with Xi_h(R). For local P2 this is draft s.7.5(ii) with interface item (h), and it has rank-3 features with no conifold counterpart.
- The vertex identity for R empty (the gap): (VI) of the LP2 draft, F4 at n = 0 for the conifold.
- The frame-B formula identifying the graph sum with the GW potential.

**Range.** Graphs with an edge use Theorem L only for 2 <= h <= g - 1. The single-vertex graph also uses h = g.

