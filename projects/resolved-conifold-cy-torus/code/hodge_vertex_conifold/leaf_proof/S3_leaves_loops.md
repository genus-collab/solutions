# S3. Dilaton leaves (L1), the H-row, and the E_11 flow (L3; L2 for j = 1)

Section S3 of the all-order proof of Lemma C (leaf_proof/outline.md): Steps 4, 7 and 8 of the outline. Written 2026-10-31 and revised
the same day after a referee report. The revision tags each input by its source ([S1.d], [S2.4c], [S2.7], [S2.9c]) and states the
dependency chain explicitly. It cites facts proved in S2 instead of proving them again, states the ground field, and aligns the
notation with S2, S4 and leaf_structure.md (0.5).
No topological recursion and no spectral curves are used. Everything below is a complete argument given the four inputs of 0.4 and the
algebraic facts of S2 recalled where they are used. All of these are proved in S1 and S2 without any use of S3. Checks: by hand at
the orders t^0, t^1, t^2 against leaf_structure.md (4.1); by the exact rational script leaf_proof/S3_checks.pl (Perl + Math::BigRat;
output leaf_proof/S3_checks_out.txt: 1540 exact checks at x = 2/7, 3/5, 7/3, -5/9, no failure; 4.2). No wolframscript was run.

Contents. 0: setting, inputs, notation. 1: Step 4, the dilaton leaves (Theorem L1). 2: Step 7, the H-row in sigma form
(Proposition 4: leaf structure of sigma; Proposition 5: its projection on the regular eigen-covector). 3: Step 8, the E_11 flow
(Theorem 7), its integration, F3 for E_11, Theorem L3, and L2 for j = 1 (Corollary 8). 4: checks. 5: logical summary, and where
F1-F3 are proved.


## 0. Setting

### 0.1 Ground field, rings and conventions

- Ground field. eps1 = x, eps2 = 1, c = 1 + x. S1 and S2 treat x as a complex parameter with x != 0, -1, and they work with
  convergent germs in t. This section is purely formal. It is written over K := Q(x), with x an indeterminate. The arguments use
  only the field operations of K, the facts x != 0 and c != 0, the formal derivation d/dt, and partial fractions over an algebraic
  closure of K. They therefore apply verbatim with K replaced by C and x a fixed complex number, x != 0, -1.
  Passage between the two settings. Every t- and z-coefficient of the objects taken from S1 and S2 is a rational function of x with
  rational coefficients (S1 s.0 and Corollary 7.1; S2 (0.1), (0.4)). These objects are tau, s*, delta := (c^2 - 4 x tau)^(1/2),
  dphi, alg, rho, L_i, Delta_i, the rows r_{a,i}[m], the F_n, and hence G and the edge and dilaton factors. S1 and S2 prove their
  identities between such rational functions for every admissible complex x, and a non-zero rational function of x has only
  finitely many zeros. So these identities hold in Q(x), and the inputs of 0.4 hold over K. (S1 extends its one analytic input,
  Proposition 6.1, from real x > 0 to all x in the same way.)
  Conversely, run the same argument over C. Every upper bound proved below then holds for every complex x with x != 0, -1: pole
  orders, membership in Lam_m, the support bound of L3, and F3 for E_11; in short, L1-L3. Only the sharpness statements need x
  generic. The top coefficient (5) and the corner (17) are non-zero iff P_k(x) != 0 (and x != 0 for (17)), i.e. for x outside a
  finite set depending on k (1.6). Lemma 1(a) also notes that certain poles are exactly simple, because c is not rational in K;
  only "at most simple" is used.
- Every function of t below (tau, s*, delta, L_i, l, Delta_i, rho, Theta, ...) is an element of the field K((t)) of formal Laurent
  series. It is the expansion at t = 0 (Q = e^t = 1) of the corresponding germ of S2 (0.2), with U_1 the singular root (U_1(0) = 1)
  and U_2(0) = -1/x. tau := 1 - e^-t = t - t^2/2 + ... is in Q[[t]]. d/dt is the derivation of K((t)), extended coefficientwise to
  series in the auxiliary variables z, w, psi, psi'.
- [t^k] denotes the coefficient of t^k, taken coefficientwise in the auxiliary variables. It is applied to three kinds of object.
  (i) Elements of K((t)) (L_1, l, Theta, rho, rho^-1, ...).
  (ii) Series in K((t))[[psi]] or K((t))[[psi, psi']]. All of these lie in t^-2 K[[t]][[psi, psi']]: sigma is in
  t^-1 K[[t]][[psi]]^2 (Proposition 4), rho^2 E_11 is in t^-2 K[[t]][[psi, psi']] by (12), and q_t, g, Dg, Y are in K[[t]][[psi]]
  (= K[[psi]][[t]]). So products are expanded order by order with finite sums: [t^k](a f) = Sum_m [t^m] a [t^(k-m)] f.
  (iii) G and z G, which lie in K(z)[[t]] (Lemma 2). Here [t^k] gives an element of K(z).
  An element of K(z)[[t]] is compared with a series in z or psi only through (KF) of 0.4: G_k is regular at z = 0 and its Taylor
  series is [t^k](rho r_{0,1}).
  K((t))[[psi]] and K((t))[[psi, psi']] are integral domains (power series rings over a field). Non-zero factors such as psi,
  psi psi' and psi + psi' are cancelled in them without comment.
- Rational functions regular at 0 are identified with their Taylor series. This is an injective K-linear map, compatible with
  products and with the substitution z -> -psi. Thus lam_j(psi) = 1/(c - j psi) = Sum_{n>=0} j^n psi^n/c^(n+1) and lam_0 = 1/c.
  Lam_m := span_K(lam_1, ..., lam_m) has dimension m, because the lam_j, j >= 1, have the distinct poles c/j. In the R-matrix
  variable z = -psi one has lam_j(psi) = 1/(c + j z).
- The pairing on covectors is <a, b> := Sum_{i,j} eta^{ij} a_i b_j = x^2 (a_0 b_1 + a_1 b_0) - c x a_1 b_1. This is etainvSym of
  rc_lib.m, with eta^-1 = [[0, x^2], [x^2, -c x]] (S2 Lemma S2.1(a)). It is symmetric with constant coefficients, so
  d/dt <a, b> = <da/dt, b> + <a, db/dt>.
- v_1(z) = (r_{0,1}(z), r_{1,1}(z)) are the frame-B rows of the singular idempotent (rmatrixAt0, then twistRows with
  hCoeffs[w, -1, .]). The twist multiplies both rows by the same t-independent series h_w(z)^-1, and h_w(z) h_w(-z) = 1. So the
  H-row relation, the QDE and the symplectic identity have the same form in the point frame and in frame B (S2 (0.4), (0.5),
  (5.4)). Edges and dilaton leaves are those of evalQ: E_11(z, w) = Sum E_11(m, n) z^m w^n and Dil_1(d) = (-1)^d r_{0,1}[d-1]. In
  leaf language, E_11(psi, psi') means E_11(z = psi, w = psi') (no sign change), as in hvc_leaf_structure.wls and V8.

Elementary leaf identities (j, m >= 1), used throughout:

    (LAM1) [psi^n] lam_j = j^n/c^(n+1).
    (LAM2) -psi lam_j = (1 - c lam_j)/j.
    (LAM3) (c - m psi) lam_j = m/j + c (1 - m/j) lam_j;  in particular (c - j psi) lam_j = 1.
    (LAM4) For f = Sum_{j=1}^M f_j lam_j:  f = -(1/psi) Sum_j f_j/j + O(psi^-2) at psi = oo  (lam_j = -1/(j psi) - c/(j^2 psi^2) - ...),
           so f = O(psi^-2) iff Sum_j f_j/j = 0.
    (LAM5) lam_i lam_j = (i lam_i - j lam_j)/(c (i - j)) for i != j  (e.g. lam_1 lam_2 = (2 lam_2 - lam_1)/c).

Proof of (LAM2), (LAM3): (c - j psi) lam_j = 1 gives psi lam_j = (c lam_j - 1)/j, which is (LAM2); and
(c - m psi) lam_j = 1 + (j - m) psi lam_j = 1 + ((j - m)/j)(c lam_j - 1), which is (LAM3).

### 0.2 Two lemmas

(PF) Partial fractions. Let S be a finite subset of K and f in K(z) a rational function whose poles all lie in S and are of
order <= 1, and which is bounded at z = oo. Then f = f(oo) + Sum_{p in S} Res_{z=p}(f)/(z - p).
Proof. Call the difference h. It has no pole in an algebraic closure of K, and h(oo) = 0. A rational function without poles is a
polynomial, and a polynomial with limit 0 at oo is 0. QED
For S = {-c/j : 1 <= j <= M} (pairwise distinct and non-zero, as c != 0) this reads
    f = f(oo) + Sum_{j=1}^M beta_j/(c + j z),   beta_j = j Res_{z=-c/j} f = lim_{z -> -c/j} (c + j z) f(z).

(T) tau-orders versus t-orders. Let V be a K-vector space and f = Sum_n f_n tau^n in V[[tau]]. Since tau is in t + t^2 Q[[t]],
V[[tau]] = V[[t]] and
    [t^k] f = Sum_{n=0}^{k} u_{k,n} f_n,   u_{k,n} := [t^k] tau^n in Q,   u_{k,k} = 1,  u_{k,n} = 0 for n > k.
Conversely t = -log(1 - tau) is in tau + tau^2 Q[[tau]], which gives the inverse unitriangular relation. Hence, for every increasing
chain of subspaces V_0 c V_1 c V_2 c ... of V,
    ([tau^n] f in V_n for all n <= k)  <=>  ([t^n] f in V_n for all n <= k)   ("tau = 1 - e^-t preserves the filtration").
Likewise, multiplication by a series a in K[[t]] that does not involve the auxiliary variables preserves such a filtration:
[t^k](a f) = Sum_m a_m [t^(k-m)] f.

### 0.3 The closed form (notation of S1)

    F_n(z) := prod_{j=1}^n (1 + j z)^2 / ( n! z^n prod_{j=2}^{n+1} (c + j z) )   (F_0 = 1),     F(z, tau) := Sum_{n>=0} F_n(z) tau^n in K(z)[[tau]].

Coefficientwise F = 2F1(1/z + 1, 1/z + 1; c/z + 2; tau), because (1/z + 1)_n = z^-n prod_{j=1}^n (1 + j z) and
(c/z + 2)_n = z^-n prod_{j=2}^{n+1} (c + j z).
Put phi_tau(s) := log s + x log(1 - s) - log(1 - tau s). Its critical points satisfy x tau s^2 - c s + 1 = 0: multiply
phi_tau'(s) = 1/s - x/(1-s) + tau/(1 - tau s) by s (1-s)(1 - tau s). Let s*(tau) in K[[tau]] be the root with s*(0) = 1/c, i.e.
s* = (1 + x tau s*^2)/c = 1/c + x tau/c^3 + O(tau^2). Put dphi(tau) := phi_tau(s*) - phi_0(1/c). Since s* is critical,
d(dphi)/dtau = (d/dtau phi_tau)(s*) = s*/(1 - tau s*) (S1 s.7.4; S2 Lemma S2.2(e)). This lies in K[[tau]], and dphi(0) = 0.
So dphi is in tau K[[tau]], with dphi = tau/c + (1 + 2x) tau^2/(2 c^3) + O(tau^3). Next, alg := (1 - tau s*) (phi''_tau(s*)/
phi''_0(1/c))^(1/2) = (delta/c)^(1/2) (c + delta)/(2c) is in 1 + tau K[[tau]] (S1 s.0, s.7.5). Exp[-dphi/z] :=
Sum_l (-1)^l dphi^l/(l! z^l) converges tau-adically in K[1/z][[tau]]. We put

    X_m := [tau^m] Exp[-dphi/z]  in K[1/z],      H_n := [tau^n](F Exp[-dphi/z]) = Sum_{n1+n2=n} F_{n1} X_{n2}  in K(z),
    Psi_m := [tau^m](F Exp[-dphi/z] alg)  in K(z)     (the functions of S1 Corollary 7.1),
    G(z; t) := (c/(c+z)) (tau/t) F(z, tau) Exp[-dphi(tau)/z]  in K(z)[[t]]     (K(z)[[tau]] = K(z)[[t]] by (T)),    G_k := [t^k] G  in K(z).

The series g, g_k and Dg in the leaf variable are defined after (KF) in 0.4.

### 0.4 Inputs from S1 and S2; the key formula; order of proof

S3 uses four results of S1 and S2, cited by the tags in brackets. It also uses some algebraic facts of S2, recalled where they are
used (2.1, 2.2, 2.4, 2.5, 3.1, 3.2). None of these depends on S3.

[S1.d] S1 Theorem S1(d), in the precise form of S1 Theorem 7.2. It rests on S1 Theorem S1(b), (c) and Proposition 6.1. For every
     k >= 0,
         [t^k] r_{0,1}(z; t) = Sum_{j<=k} ([t^k] tau^j) (c/(c + z)) Psi_j(z).
     Every Psi_j has no pole at z = 0 (S1 Corollary 7.1; second, purely formal proof in S1 s.7.4, (7.3)). The Taylor coefficients
     in z of this rational function are the t^k-coefficients of the z-coefficients of r_{0,1}. In particular r_{0,1}(z; 0) =
     c/(c + z) (F2; S1 Corollary 7.3).
[S2.4c] S2 Proposition S2.4(c). rho = tau/(t alg) is a unit of K[[t]] with rho(0) = 1, and t^2 rho^2 Delta_1 = -x c^3 exactly (F1
     for Delta_1).
[S2.7] S2 Proposition S2.7. Let D_1 f := df/dt - (dDelta_1/dt)/(2 Delta_1) f. The H-row is defined (rmatrixCore; S2 (0.5)) by
         (H)  r_{1,1} = (L_1 + z D_1) r_{0,1},
     and the second row of the QDE holds in K((t))[[z]]:
         z D_1 r_{1,1} = L_2 (r_{1,1} - L_1 r_{0,1}).
     Equivalently (D_1 r_{1,1})[m-1] = L_2 (r_{1,1}[m] - L_1 r_{0,1}[m]) for m >= 1, the form tested in V5. The z^0 term of the right
     side vanishes because r_{1,1}[0] = L_1 r_{0,1}[0]. In matrix form z D_1 v_1 = (A - L_1) v_1, with A as in 2.4; this is the
     outline's z ds/dt = A s.
[S2.9c] S2 Proposition S2.9(c). K_1 = 1, i.e. <v_1(w), v_1(-w)> = Delta_1 in K((t))[[w]] (V6 for i = j = 1).

(KF) The key formula with regularity (S2 Proposition S2.4(d), made precise). G is in K(z)[[t]] (0.3), and rho r_{0,1}(z; t) is
in K[[t]][[z]]. For every k >= 0 the rational function G_k = [t^k] G has no pole at z = 0, and its Taylor series at z = 0 is
[t^k](rho r_{0,1}). In particular G_0 = c/(c + z).
Proof. By [S2.4c], rho alg = tau/t, so G = rho (c/(c + z)) F Exp[-dphi/z] alg = Sum_j rho (c/(c + z)) Psi_j tau^j in K(z)[[t]], and
    G_k = Sum_{m=0}^{k} ([t^m] rho) Sum_{j<=k-m} ([t^(k-m)] tau^j) (c/(c + z)) Psi_j.
This is a finite K-linear combination of rational functions regular at z = 0. By [S1.d] its Taylor series is
Sum_m ([t^m] rho) [t^(k-m)] r_{0,1} = [t^k](rho r_{0,1}). For k = 0, [t^0] rho = 1, [t^0] tau^j = delta_{j0} and
Psi_0 = F_0 X_0 alg(0) = 1. QED
Put g(psi; t) := (rho r_{0,1})(-psi; t) in K[[t]][[psi]]. By (KF), g_k := [t^k] g is the Taylor series of G_k(-psi); this is the
meaning of "g = G(-psi)" coefficientwise in t. Put Dg(psi; t) := d/dt (t g(psi; t)) = Sum_{k>=0} (k+1) g_k(psi) t^k.

Order of proof (no circularity). The dependencies are
    S1 Theorem S1(b)-(d), Corollary 7.1 (or s.7.4)  ->  S2 Proposition S2.4(c),(d)  ->  S2 Proposition S2.7  ->  S2 Proposition S2.9  ->  S3.
In detail:
- (KF) is [S1.d] plus [S2.4c].
- [S2.7] is proved in S2 from S1's closed form (S2 (S1.a), (S1.b)), Gauss's equation, and the operator identity S2 Lemma S2.6.
- [S2.9c] is proved in S2 from S2 Proposition S2.9(a) and (b). Part (a) says that K_1 is t-independent, by [S2.7] and
  self-adjointness. Part (b) derives F2 and [t^-1] sigma = (0, c) -- the statement of Proposition 4(a) below -- directly from (H)
  and S1's closed form.
So S2 uses nothing from S3.
Within S3, Proposition 4(a) uses (KF) at k = 0 (F2), L_1 = c/t + l (S2 Lemma S2.2(d)), and formula (6) of 2.1. Formula (6) needs
[S2.4c], and in fact only F1: with rho defined by the outline's formula, t^2 rho^2 Delta_1 = -x c^3 holds by definition. [S2.9c] is
used only in section 3.

### 0.5 Notation shared with other sections

For assembly with S1, S2, S4 and leaf_structure.md:
- a_j^(k) is the lam_j-coefficient of [t^k](rho q_t), equivalently the coefficient of 1/(c + j z) in [t^k](z G). With
  A_j(t) := Sum_k a_j^(k) t^k, these are the a_j^(k) and A_j(t) of S4 ((I2), Theorem CF.1), and a_{k+1}^(k) is the top coefficient
  of leaf_structure.md, Step 3. B_j := t A_j, as in S4 CF.2(b). S4 writes A_j(t) = (-1)^j C_j t^-1 kappa(t)^j 2F1(j x/c, j x/c; j; tau)
  with kappa := tau e^(dphi/c); kappa is not used here.
- g_j^(k) := -(j/c) a_j^(k) is the lam_j-coefficient of g_k. S4 does not name it; it writes g_k = -(1/c) Sum_j j a_j^(k) lam_j.
- e^(k)_il are the loop coefficients of leaf_structure.md (iii) (written e_il^(k) in S4 Theorem CF.4). delta_k := [t^k](rho q_t)
  (k >= 1) is as in leaf_structure.md (i). It is not to be confused with delta = (c^2 - 4 x tau)^(1/2) of S1, S2, or with the
  Kronecker delta_{k0}.
- Lam_m = span(lam_1, ..., lam_m) is V_m of S4 (N2).
- sigma = (sigma_0, sigma_1) denotes the two components, and sigma^(k) := [t^k] sigma. S4 uses the opposite convention: it writes
  sigma_k for [t^k] sigma and, in the proof of its Lemma 9.1, sigma^(1) for the second component.
- chi_i := (1, L_i), pi, A and the rank-one formula (10) are those of S2 Proposition S2.8. Dg := d/dt (t g) (written d/dt (t g) in
  S2.11 and S4). The letter omega (S2 Lemma S2.6: omega = (c + z)/tau - (1 + z)) is not used here.
- P_k := prod_{i=1}^{k} (k + 1 - i c)^2 (P_0 = 1), as in leaf_structure.md; C_j := P_{j-1}/(j! (j-1)! c^(2j-4)), as in S4 CF.1.
  S2's projector P_2 is not used here.
- K := Q(x) here; S2 uses K for C((z)). Y := t sigma, Ebar := t^2 rho^2 E_11 and Nbar := (psi + psi') Ebar are local to section 3.
  S2.10 uses Y for another series, and S1 s.7.4 uses Phi and N for other objects. In S2's proof of Proposition S2.9(c), (d),
  Y = kappa_0 refl(w).


## 1. Step 4: the dilaton leaves (L1)

### 1.1 The tau-coefficients

Lemma 1. Let n, m >= 0.
 (a) F_n has poles only at z = 0 (order exactly n) and at the points z = -c/j, 2 <= j <= n+1 (order <= 1; exactly 1 over Q(x),
     where c is not rational), and F_n(oo) = 1/(n+1).
 (b) X_m = [tau^m] Exp[-dphi/z] is a polynomial in 1/z of degree <= m, with X_m(oo) = delta_{m0}.
 (c) H_n = Sum_{n1+n2=n} F_{n1} X_{n2} has poles only at z = 0 (order <= n) and at z = -c/j, 2 <= j <= n+1 (order <= 1), and
     H_n(oo) = 1/(n+1).
Proof. (a) The numerator prod_{j=1}^n (1 + jz)^2 has degree 2n, value 1 at z = 0, and leading coefficient (n!)^2. The denominator
n! z^n prod_{j=2}^{n+1}(c + jz) has degree 2n and leading coefficient n! (n+1)!. It has a zero of order n at 0, and simple zeros at
the n points -c/j (2 <= j <= n+1), which are distinct and non-zero. A zero of the numerator can only lower a pole order: the zero
-1/i of 1 + i z coincides with -c/j iff c = j/i. So the poles are as stated, and F_n(oo) = (n!)^2/(n! (n+1)!) = 1/(n+1).
(b) dphi in tau K[[tau]] gives dphi^l in tau^l K[[tau]], hence X_m = Sum_{l=0}^{m} ((-1)^l/l!) z^-l [tau^m](dphi^l). This is a
polynomial in 1/z of degree <= m, whose constant term is [tau^m] dphi^0 = delta_{m0}.
(c) Product of (a) and (b). X_{n2} is regular away from 0 and has a pole of order <= n2 at 0, so pole orders add at 0 and do not
change elsewhere. At oo, F_{n1} X_{n2} -> delta_{n2,0}/(n1 + 1). QED

### 1.2 The t-coefficients of z G

Lemma 2. For k >= 0 put R_k(z) := [t^k](z G) = z G_k(z). Then

    R_k = (c z/(c+z)) Sum_{n=0}^{k} d_{k,n} H_n,     d_{k,n} := [t^k]((tau/t) tau^n) in Q,   d_{k,k} = 1.         (1)

So R_k is a rational function with poles only at z = 0 and at z = -c/j, 1 <= j <= k+1, of order <= 1 at each -c/j. It is bounded
at oo, and
    R_k(oo) = c delta_{k0}.                                                                                         (2)
Proof. z G = (c z/(c+z)) (tau/t) Sum_n H_n tau^n. As tau/t is in 1 + t Q[[t]] and tau^n in t^n (1 + t Q[[t]]), d_{k,n} = 0
for n > k and d_{k,k} = 1. This is (T) for the series (tau/t) F Exp[-dphi/z] in V[[tau]], V = K(z), and it gives the finite sum
(1).
Poles. By Lemma 1(c), each H_n with n <= k has poles only at 0 and at -c/j with 2 <= j <= n+1 <= k+1, simple at the latter. The
factor c z/(c+z) is regular except for a simple pole at -c (this is j = 1), and -c is not one of the points -c/j, j >= 2.
Multiplying by a function regular at a point does not raise the pole order there.
Limit. c z/(c+z) -> c and H_n -> 1/(n+1), so by (1)
    R_k(oo) = c Sum_{n=0}^{k} d_{k,n}/(n+1) = c [t^k] ( (tau/t) Sum_{n>=0} tau^n/(n+1) ) = c [t^k] ( -log(1 - tau)/t ) = c [t^k] 1 = c delta_{k0}.
The second equality holds because the terms n > k do not contribute to [t^k], and -log(1 - tau) = -log(e^-t) = t in Q[[t]]. QED

Remark (the coefficientwise limit). F(oo, tau) := Sum_n F_n(oo) tau^n = Sum_n tau^n/(n+1) = -log(1 - tau)/tau is the limit of
F(z, tau) coefficient by coefficient in tau, and only this is used. For each k, R_k is the finite combination (1) of the rational
functions H_0, ..., H_k, and the limit z -> oo is taken in each of them. No analytic statement about 2F1(1/z + 1, 1/z + 1; c/z + 2;
tau) as z -> oo at fixed tau is needed or claimed. For k >= 1 the limit vanishes because (tau/t) F(oo, tau) = (tau/t)(t/tau) = 1 is
t-independent. This is the precise sense in which the normalisation by rho removes lam_0. Equivalently, lim_{z -> oo} z r_{0,1}(z; t)
= c/rho(t), coefficientwise in t.

### 1.3 No pole at z = 0

By (KF), G_k is regular at z = 0, so R_k = z G_k is regular at 0 and R_k(0) = 0. In the closed form, the separate terms of (1) do
have poles at 0 (z H_n can have a pole of order n - 1). That they cancel is S1's regularity theorem, and nothing in this section
depends on how it is proved. Explicitly, every H_n is regular at z = 0, for either of two reasons.
 (a) S1 Corollary 7.1: Psi_m = [tau^m](F Exp[-dphi/z] alg) has no pole at z = 0, for every m. Since alg is z-independent with
     1/alg in 1 + tau K[[tau]],
         H_n = [tau^n]( (F Exp[-dphi/z] alg) (1/alg) ) = Sum_{m=0}^{n} ([tau^(n-m)] (1/alg)) Psi_m
     is a finite K-linear combination of Psi_0, ..., Psi_n, hence regular at 0. No induction is needed.
 (b) S1 s.7.4, (7.3): F Exp[-dphi/z] equals S1's transport solution, a series in z whose coefficients are analytic in tau. So the
     Laurent expansion at z = 0 of its tau^n-coefficient H_n has no negative powers. This proof is purely formal: it uses neither
     the Laplace identification nor the value of S1's constant A(z).
Then (1) shows directly that R_k is regular at 0 with R_k(0) = 0, since each term has the factor c z/(c + z). Conversely, the
regularity of all H_n gives that of all G_k by (T), since c/(c + z) is a unit at z = 0 and tau/t is a z-independent unit of
K[[t]]. So "every G_k is regular at 0" and "every H_n is regular at 0" are equivalent. (Checked by hand for n = 1, 2 in 4.1:
H_1 = (2x + c z)/(c (c + 2z)), and the z^-2 and z^-1 parts of H_2 = F_2 + F_1 X_1 + X_2 cancel.)

### 1.4 Partial fractions

Proposition 3 (L1 in the variable z). R_0 = c z/(c+z) = c - c^2/(c+z). For k >= 1 there are a_j^(k) in K, 1 <= j <= k+1, with

    R_k(z) = [t^k](z G) = Sum_{j=1}^{k+1} a_j^(k)/(c + j z),    a_j^(k) = j Res_{z=-c/j} R_k,    Sum_{j=1}^{k+1} a_j^(k) = 0.     (3)

Proof. For k = 0, (tau/t) F Exp[-dphi/z] = H_0 = 1 at t = 0, so G_0 = c/(c+z). For k >= 1, Lemma 2 and 1.3 show that the poles
of R_k lie in {-c/j : j <= k+1} and are simple, and R_k(oo) = 0. So (PF) gives (3) without constant term. Evaluating at z = 0
gives 0 = R_k(0) = (1/c) Sum_j a_j^(k). QED
We also write a_1^(0) := -c^2 (the coefficient of 1/(c+z) in R_0) and a_j^(k) := 0 for j > k+1. The sum rule then reads
Sum_j a_j^(k) = -c^2 delta_{k0}, as in S4 (I2).

### 1.5 Translation to psi and to q_t

Theorem L1. For every k >= 0,

    [t^0](rho q_t) = q_0 = c - c^2 lam_1,      [t^k](rho q_t) = delta_k := Sum_{j=1}^{k+1} a_j^(k) lam_j  in Lam_{k+1}   (k >= 1),

with the a_j^(k) of (3). Moreover g_0 = c lam_1, and for k >= 1

    g_k = Sum_{j=1}^{k+1} g_j^(k) lam_j = -(1/c) Sum_{j=1}^{k+1} j a_j^(k) lam_j  in Lam_{k+1},     g_k = O(psi^-2) at psi = oo.     (4)

Proof. By definition q_t(psi) = Sum_{d>=1} Dil_1(d; t) psi^d with Dil_1(d) = (-1)^d r_{0,1}[d-1] (DilN in evalQ), so
    q_t(psi) = Sum_{m>=0} (-1)^(m+1) r_{0,1}[m] psi^(m+1) = -psi r_{0,1}(-psi; t).
Multiply by rho, which does not involve psi, and use (KF): rho q_t(psi) = -psi (rho r_{0,1})(-psi; t) = -psi g(psi; t). Hence
[t^k](rho q_t)(psi) = (z G_k(z))|_{z=-psi} = R_k(-psi) as formal series (Taylor expansion commutes with z -> -psi). As
1/(c + j z) = lam_j(psi) at z = -psi, Proposition 3 gives the first two statements.
For (4): G_k = R_k/z and 1/(z (c + jz)) = (1/c)(1/z - j/(c + jz)), so by Sum_j a_j^(k) = 0,
    G_k(z) = (1/(c z)) Sum_j a_j^(k) - (1/c) Sum_j j a_j^(k)/(c + j z) = -(1/c) Sum_j j a_j^(k)/(c + j z).
Substitute z = -psi. Since R_k = O(1/z) at oo (no constant term), G_k = R_k/z = O(z^-2), i.e. g_k = O(psi^-2). Equivalently, by
(LAM4), Sum_j g_j^(k)/j = -(1/c) Sum_j a_j^(k) = 0. QED

Remarks. (i) g_j^(k) = -(j/c) a_j^(k) for all k >= 0; at k = 0, g_1^(0) = c = -(1/c) a_1^(0). (ii) Conversely, by (LAM2), L1 at
order k >= 1 is equivalent to "g_k in Lam_{k+1} and g_k = O(psi^-2)". Indeed -psi g_k = Sum_j g_j^(k)/j + Sum_j a_j^(k) lam_j, and
the constant Sum_j g_j^(k)/j must vanish. This is the form used in section 2. (iii) The three ingredients are separate: the pole
structure (Lemmas 1, 2; Fact G of conifold_regularity_note.md s.6), S1's regularity theorem (no pole at z = 0; 1.3), and the
coefficientwise limit at z = oo (no lam_0; Remark in 1.2).

### 1.6 The leaf coefficients as functions of t; the top coefficient

Put A_j(t) := Sum_{k>=0} a_j^(k) t^k in K[[t]] and z_j := -c/j. Then, coefficientwise in t (only j <= k+1 contribute to [t^k]),
    rho q_t(psi) = c + Sum_{j>=1} A_j(t) lam_j(psi),      g(psi; t) = -(1/c) Sum_{j>=1} j A_j(t) lam_j(psi).
The constant c occurs only at t^0: q_0 = c + a_1^(0) lam_1.

Proposition (residue and hypergeometric form). In K[[t]],

    A_1(t) = -c^2 (tau/t) F(-c, tau) e^(dphi(tau)/c),       F(-c, tau) = Sum_n ((x/c)_n)^2 tau^n/(n!)^2 = 2F1(x/c, x/c; 1; tau),
    A_j(t) = j (c z_j/(c + z_j)) (tau/t) e^(j dphi(tau)/c) Sum_{n >= j-1} Res_{z=z_j}(F_n) tau^n          (j >= 2),
           = (-1)^j C_j t^(j-1) ((tau/t) e^(dphi/c))^j 2F1(j x/c, j x/c; j; tau),     C_j := P_{j-1}/(j! (j-1)! c^(2j-4)).

The last line holds for j = 1 too (C_1 = c^2). Here P_k := prod_{i=1}^{k} (k + 1 - i c)^2 (P_0 = 1).
Proof. a_j^(k) = j Res_{z_j} R_k, and R_k is the finite sum (1).
Residue form, j >= 2. Only the F_{n1} with n1 >= j-1 have a pole at z_j, and it is simple. All other factors (c z/(c+z), the
d_{k,n}, the X_{n2}) are regular at z_j, so Res_{z_j}(u F_{n1}) = u(z_j) Res_{z_j} F_{n1} for each such factor u. Re-sum over k
(finitely many terms at each order) and use Sum_m X_m(z_j) tau^m = Exp[-dphi/z_j] = e^(j dphi/c) in K[[tau]].
j = 1. The pole at -c comes only from c z/(c+z), since F_n and X_m are regular at -c, and Res_{-c}(c z/(c+z)) = -c^2. So
A_1(t) = -c^2 (tau/t) Sum_n H_n(-c) tau^n, and Sum_n H_n(-c) tau^n = F(-c, tau) e^(dphi/c) because Sum_m X_m(-c) tau^m = e^(dphi/c).
Moreover F_n(-c) = prod_{i=1}^n (i - 1/c)^2/(n!)^2 and (x/c)_n = prod_{i=1}^n (i - 1/c).
Hypergeometric form, j >= 2, with n = j - 1 + m. Using 1 + i z_j = (j - i c)/j and c + i z_j = c (j - i)/j,
    Res_{z_j} F_n = (-1)^(j+1) prod_{i=1}^n (j - i c)^2 / ( j^2 n! c^(2n-1) (j-2)! (n+1-j)! ).
Next, prod_{i=j}^{j-1+m} (j - i c) = (-c)^m (j x/c)_m and n! = (j-1)! (j)_m, so
    Sum_{n>=j-1} Res_{z_j}(F_n) tau^n = (-1)^(j+1) P_{j-1}/(j^2 (j-1)! (j-2)! c^(2j-3)) tau^(j-1) 2F1(j x/c, j x/c; j; tau).
With c z_j/(c + z_j) = -c/(j-1) this gives the hypergeometric form. QED
Remarks. The Proposition is S4 Theorem CF.1, proved there by the same residue computation. S4 writes it as
A_j(t) = (-1)^j C_j t^-1 kappa(t)^j 2F1(j x/c, j x/c; j; tau) with kappa := tau e^(dphi/c), which is the last line above. Since
A_j(t) = O(t^(j-1)), it re-proves a_j^(k) = 0 for j > k+1. Nothing in L1-L3 depends on it. This section uses only its corollary (5),
for the sharpness statements, and the direct computation below proves (5) without it.

Corollary (top coefficient). For k >= 0,

    a_{k+1}^(k) = (-1)^(k+1) P_k/(k! (k+1)! c^(2k-2)),      g_{k+1}^(k) = (-1)^k P_k/((k!)^2 c^(2k-1)).                       (5)

This is the top-coefficient formula of leaf_structure.md, Step 3 (where the same symbol a_{k+1}^(k) is used, verified there for
k <= 5), and of S4 CF.3.
P_k is a non-zero polynomial in x, since each factor (k + 1 - i) - i x is non-zero. For every x with P_k(x) != 0, i.e. with
c not in {(k+1)/i : 1 <= i <= k}, the top leaf lam_{k+1} occurs at order t^k. So the bound in L1 is sharp for generic x.
Direct computation (k >= 1). Put z_0 := -c/(k+1). Only the term n = k of (1), and in H_k only F_k, has a pole at z_0. Indeed the
F_n with n < k and all X_m are regular there, and d_{k,k} = 1. So a_{k+1}^(k) = (k+1) Res_{z_0} R_k = (k+1) (c z_0/(c + z_0))
Res_{z_0} F_k. At z_0 one has c z_0/(c + z_0) = -c/k, 1 + i z_0 = (k + 1 - i c)/(k+1), z_0^k = (-c)^k/(k+1)^k,
prod_{i=2}^{k} (c + i z_0) = c^(k-1) (k-1)!/(k+1)^(k-1), and Res_{z_0} 1/(c + (k+1) z) = 1/(k+1). Hence
Res_{z_0} F_k = (-1)^k P_k/((k+1)^2 k! (k-1)! c^(2k-1)), and a_{k+1}^(k) is as in (5). For k = 0, a_1^(0) = -c^2 = -P_0 c^2. The
formula for g_{k+1}^(k) follows from g_j^(k) = -(j/c) a_j^(k).


## 2. Step 7: the H-row in sigma form

### 2.1 The H-row

By [S2.4c], Delta_1 = -x c^3 (t rho)^-2. Hence (dDelta_1/dt)/(2 Delta_1) = -(d(t rho)/dt)/(t rho), and

    D_1 f = (t rho)^-1 d/dt (t rho f).                                                                              (6)

This is S2 (5.2.1) with nu_1 = t rho/kappa_0. No square root of Delta_1 is needed, and only F1 is used. With rho defined by the
outline's formula, t^2 rho^2 Delta_1 = -x c^3 holds by definition, so the logarithmic derivative of Delta_1 is that of (t rho)^-2.
By (H), (6) and (KF):
    rho r_{1,1}(z) = rho (L_1 + z D_1) r_{0,1} = L_1 G + (z/t) d/dt (t G) = (L_1 + z/t) G + z dG/dt
(S2 Proposition S2.9(b); S2 Remark S2.7.2). Define, as in the outline and S2 Proposition S2.9(b),

    sigma(psi) := rho v_1(-psi) = (sigma_0(psi), sigma_1(psi)),   sigma_0 = g,   sigma_1 = L_1 g - (psi/t) Dg = (L_1 - psi/t) g - psi dg/dt.     (7)

The substitution z = -psi turns z/t into -psi/t and z dG/dt into -psi dg/dt. Both components of t sigma lie in K[[t]][[psi]],
because t L_1 is in K[[t]] (2.2).

### 2.2 The pole of L_1

From LofU, L(U) = -Q (x U + 1)/(1 - Q) = (x U + 1)/tau (Q = 1/(1 - tau)), and L_1 = L(U_1) = c/tau - x s* (S2 Lemma S2.2(b), (d)).
Since 1/tau = 1/t + 1/2 + t/12 - t^3/720 + ... (1/tau - 1/t is in Q[[t]]),

    L_1 = c/t + l(t),   l := c (1/tau - 1/t) - x s*  in K[[t]],   l_0 = c/2 - x/c = (1 + x^2)/(2c),   l_1 = c/12 - x^2/c^3.

### 2.3 The leaf structure of sigma

Proposition 4. (a) [t^-1] sigma = (0, c).
(b) For k >= 0, [t^k] sigma_0 = g_k and

    [t^k] sigma_1 = (c - (k+2) psi) g_{k+1} + Sum_{m=0}^{k} l_m g_{k-m}   in Lam_{k+1}.                                    (8)

Hence sigma = t^-1 (0, c) + Sum_{k>=0} t^k sigma^(k), with both components of sigma^(k) in Lam_{k+1}. Explicitly

    sigma^(0) = ( c lam_1, x lam_1 ),
    sigma^(1) = ( ((1-x)^2/(2c)) (lam_1 - 2 lam_2),  (x/(2c^2)) [ (1 - x + x^2) lam_1 - 2 (1-x)^2 lam_2 ] ).

Proof. From (7) and 2.2, [t^k](L_1 g) = c g_{k+1} + Sum_{m=0}^{k} l_m g_{k-m} and [t^k]((psi/t) Dg) = psi [t^(k+1)] Dg =
(k+2) psi g_{k+1}. This gives (8) for every k >= -1 (empty sum for k = -1). For k = -1, [t^-1] sigma_1 = (c - psi) g_0 =
(c - psi) c lam_1 = c by (LAM3), and [t^-1] sigma_0 = 0. This is (a).
Leaf structure, psi-form. Write g_{k+1} = Sum_{j=1}^{k+2} g_j^(k+1) lam_j. By Theorem L1 (k+1 >= 1) and (LAM4),
Sum_j g_j^(k+1)/j = 0. By (LAM3) and g_j^(k+1) = -(j/c) a_j^(k+1),
    (c - (k+2) psi) g_{k+1} = (k+2) Sum_j g_j^(k+1)/j + c Sum_j (1 - (k+2)/j) g_j^(k+1) lam_j = Sum_{j=1}^{k+1} (k+2-j) a_j^(k+1) lam_j.
The factor c - (k+2) psi turns the top leaf lam_{k+2} into the constant 1, and all the constants add up to
(k+2) Sum_j g_j^(k+1)/j = 0. The terms l_m g_{k-m} lie in Lam_{k-m+1}, which is contained in Lam_{k+1}. So, explicitly,
    [t^k] sigma_1 = Sum_{j=1}^{k+1} ( (k+2-j) a_j^(k+1) + Sum_{m=0}^{k} l_m g_j^(k-m) ) lam_j       (g_j^(n) := 0 for j > n+1).
Together with [t^k] sigma_0 = g_k = Sum_j g_j^(k) lam_j, this is the coefficient vector
S_j^(k) = (g_j^(k), (k+2-j) a_j^(k+1) + Sum_m l_m g_j^(k-m)) of S4 Lemma 9.1, (9.1).
Leaf structure, z-form (the outline's formulation). At psi = -z, (8) is (c + (k+2) z) G_{k+1}(z) + Sum_m l_m G_{k-m}(z).
G_{k+1} has at most a simple pole at z = -c/(k+2), and the factor c + (k+2) z kills it. Its remaining poles are simple, at -c/j
with j <= k+1. It is regular at 0 by (KF), and (c + (k+2) z) G_{k+1} = O(1/z) because G_{k+1} = O(z^-2) (Theorem L1). The G_{k-m}
are O(1/z), regular at 0, with simple poles at -c/j, j <= k-m+1. So (8) is a rational function with simple poles in
{-c/j : j <= k+1}, regular at 0 and vanishing at oo, and (PF) puts it in span{1/(c + j z) : 1 <= j <= k+1}.
Explicit values, k = 0. g_1 = ((1-x)^2/(2c))(lam_1 - 2 lam_2), by Theorem L1 with (a_1^(1), a_2^(1)) = ((1-x)^2/2)(-1, 1). Then
(c - 2 psi) g_1 = a_1^(1) lam_1 = -((1-x)^2/2) lam_1 and l_0 g_0 = (c^2/2 - x) lam_1, so
[t^0] sigma_1 = ((c^2 - (1-x)^2)/2 - x) lam_1 = x lam_1.
Explicit values, k = 1. (c - 3 psi) g_2 = 2 a_1^(2) lam_1 + a_2^(2) lam_2, with a_1^(2) = -(2 - 4x + 3x^2 - 4x^3 + 2x^4)/(12 c^2)
and a_2^(2) = (1-x)^4/(2c^2). Also l_0 g_1 = ((1+x^2)(1-x)^2/(4c^2))(lam_1 - 2 lam_2) and l_1 g_0 = ((c^4 - 12 x^2)/(12 c^2)) lam_1.
The lam_1-coefficient is
[-2(2 - 4x + 3x^2 - 4x^3 + 2x^4) + 3(1+x^2)(1-x)^2 + (1+x)^4 - 12x^2]/(12c^2) = 6x(1 - x + x^2)/(12c^2).
The lam_2-coefficient is (1-x)^4/(2c^2) - (1+x^2)(1-x)^2/(2c^2) = -x (1-x)^2/c^2. QED

Remarks. (i) sigma^(0)_1 = x lam_1 is the z-dependent part x/(c+z) of the t^0 term of the H-row in F2 of
conifold_regularity_note.md; the constant part is absorbed in (c/t)(1 + O(t)). In z: rho r_{1,1}(z) = c/t + x/(c + z) + O(t), as
in S4 Lemma 9.1.
(ii) sigma_1 is O(1/psi) but not O(psi^-2), since sigma^(0)_1 = x lam_1. The H-row components carry a lam_1 tail, which does not
matter for the leaf structure.
(iii) Inputs of (a): (KF) at k = 0 (F2), L_1 = c/t + l, and (6), i.e. [S2.4c] or just F1. Statement (a) is S2 Proposition
S2.9(b), which S2 proves independently from (H) and S1's closed form and then uses for K_1 = 1. S2 does not use S3 (0.4).
(iv) Proposition 4 is the input of Step 9 (S4 Lemma 9.1, which restates (8) with the coefficient vectors S_j^(k) above).

### 2.4 Spectral data of the QDE matrix (recalled from S2)

Let A := [[0, 1], [-x/tau, c/tau]] act on covectors (a_0, a_1); it is the transpose of the matrix of H * in the basis (1, H) (S2
Proposition S2.8(a)). Let chi_i := (1, L_i). The following facts are proved in S2.
- (S2 Lemma S2.2(d).) L_1 = c/tau - x s*, L_2 = x s*, L_1 + L_2 = c/tau, L_1 L_2 = x/tau, L_1 - L_2 = delta/tau. So
  L_2 - L_1 = -c/t + O(1) is invertible in K((t)).
- (S2 Proposition S2.8(b).) A chi_i = L_i chi_i.
- (S2 Proposition S2.8(c).) A is self-adjoint for < , >: eta^-1 A is symmetric, so <A a, b> = <a, A b>.
- (S2 Proposition S2.8(d), with S2 Proposition S2.3(b), (c).) <chi_i, chi_j> = delta_ij Delta_i, where Delta_i = 1/eta(e_i, e_i) is
  the code's DeltaCN. For {i, j} = {1, 2}, Delta_i = x^3 (L_i - L_j)/(x - c L_j).
- (S2 Proposition S2.3(d), (e).) x^2 - c x L_2 = x (x - c L_2) = -x^3 tau s*^2 = -x^3 t/c^2 + O(t^2), and Delta_2 = x^3/c + O(t)
  (F1 for Delta_2).
One further form of the metric is used:

    Delta_1 = (L_1 - L_2)(x^2 - c x L_1),       Delta_2 = (L_2 - L_1)(x^2 - c x L_2).                                  (9)

Proof. (x - c L_1)(x - c L_2) = x^2 - c x (L_1 + L_2) + c^2 L_1 L_2 = x^2 - c^2 x/tau + c^2 x/tau = x^2. So x^3/(x - c L_j) =
x (x - c L_i) = x^2 - c x L_i; insert this in Delta_i = x^3 (L_i - L_j)/(x - c L_j). QED
(With L_1 = c/t + O(1) and x^2 - c x L_1 = -c^2 x/t + O(1), (9) gives Delta_1 = -x c^3/t^2 + O(1/t) again. With
x^2 - c x L_2 = -x^3 t/c^2 + O(t^2), it gives Delta_2 = x^3/c + O(t).)
Finally, the rank-one formula (S2 Proposition S2.8(e)). For all covectors a, b, with pi(a) := <a, chi_2> = x^2 L_2 a_0 +
(x^2 - c x L_2) a_1,

    <(A - L_1) a, b> = ((L_2 - L_1)/Delta_2) pi(a) pi(b).                                                          (10)

S2 also writes pi(a) = x^3 s* (a_0 - tau s* a_1), and the right side of (10) as -(x^3/tau)(a_0 - tau s* a_1)(b_0 - tau s* b_1).
The forms agree because x^2 L_2 = x^3 s*, x^2 - c x L_2 = -x^3 tau s*^2, and (L_2 - L_1)/Delta_2 = (-delta/tau)/(x^3 delta s*^2)
= -1/(x^3 tau s*^2) (S2 Proposition S2.3(d)).
Self-adjointness is a K((t))-linear identity and (10) a K((t))-bilinear one. Both extend coefficientwise to covectors with entries in
K((t))[[psi]] and K((t))[[psi']], with values in K((t))[[psi, psi']]; this is how they are applied in 3.3.

### 2.5 The projection of sigma on the regular eigen-covector

Proposition 5 (S2.11, item (I2)). pi(sigma(psi)) = -(x^2 - c x L_2) (psi/t) Dg(psi), and x^2 - c x L_2 = -x^3 tau s*^2. Thus

    pi(sigma(psi)) = x^3 (tau/t) s*^2 psi Dg(psi)   in K[[t]][[psi]].

Proof. This is S2.11 (I2). In the present notation: sigma = g chi_1 + (sigma_1 - L_1 g)(0, 1), pi(chi_1) = <chi_1, chi_2> = 0,
pi((0, 1)) = x^2 - c x L_2, and sigma_1 - L_1 g = -(psi/t) Dg by (7). QED


## 3. Step 8: the E_11 flow; F3 for E_11; L3; L2 for j = 1

### 3.1 E_11 in terms of sigma

In evalQ,
    EdN[1,1][m, n] = Sum_{k=0}^{n} (-1)^k P_{m+k+1, n-k},     P_pq := -(-1)^(p+q) <v_1[p], v_1[q]> = [z^p w^q] P(z, w),
where v_1[p] := [z^p] v_1 = (r_{0,1}[p], r_{1,1}[p]) and P(z, w) := -<v_1(-z), v_1(-w)>. Since (z^p - (-w)^p)/(z + w) =
Sum_{i+j=p-1} z^i (-w)^j, EdN[1,1][m, n] is the coefficient of z^m w^n of the difference quotient
    E_11(z, w) = [P(z, w) - P(-w, w)]/(z + w).
(S4 (1.1)-(1.2) prove the same for every E_ij.) By [S2.9c], P(-w, w) = -<v_1(w), v_1(-w)> = -Delta_1. In the leaf variables
(psi, psi') = (z, w), rho^2 P(psi, psi') = -<sigma(psi), sigma(psi')>. Hence (S2 Proposition S2.9(d))

    rho^2 E_11(psi, psi') = [ rho^2 Delta_1 - <sigma(psi), sigma(psi')> ] / (psi + psi').                               (11)

Put Y(psi) := t sigma(psi), which is in K[[t]][[psi]]^2 by Proposition 4, and use [S2.4c]:

    Ebar(psi, psi') := t^2 rho^2 E_11(psi, psi') = Nbar(psi, psi')/(psi + psi'),      Nbar := -x c^3 - <Y(psi), Y(psi')>.          (12)

Nbar is in K[[t]][[psi, psi']], and by [S2.9c] and [S2.4c]
    Nbar(psi, -psi) = -x c^3 - t^2 rho^2 <v_1(-psi), v_1(psi)> = -x c^3 - t^2 rho^2 Delta_1 = 0.
Hence Nbar is divisible by psi + psi' in K[[t]][[psi, psi']]: in the variables u = psi, v = psi + psi' it vanishes at v = 0. So
Ebar is in K[[t]][[psi, psi']]. It is t^2 rho^2 times the code's difference quotient, since quotients by the non-zero-divisor
psi + psi' are unique. In particular rho^2 E_11 is in t^-2 K[[t]][[psi, psi']] (0.1).
Since < , > is symmetric, Ebar(psi, psi') = Ebar(psi', psi), i.e. E_11(a, b) = E_11(b, a); this is the symmetry observed in
leaf_structure.md, Step 3. Without [S2.9c] the subtracted term would be P(-w, w) = -K_1(w) Delta_1, and the symmetry would not be
manifest.

### 3.2 The QDE for Y

Lemma 6 (S2 Proposition S2.9(d)). -psi dY(psi)/dt = (A - L_1) Y(psi) in K((t))[[psi]]^2. Equivalently,
d sigma/dt = -sigma/t - (1/psi)(A - L_1) sigma.
Proof. The first identity is -z (t sigma(z))' = t (A - L_1) sigma(z) in the proof of S2 Proposition S2.9(d), with z renamed psi.
It is [S2.7] in the form z w' = (A - L_1) w, w := nu_1 v_1 = (t/kappa_0) rho v_1 (S2 Proposition S2.7), transported by z -> -z.
Here Y(psi) = kappa_0 w(-psi), and A, L_1 do not involve z.
For the second identity: since Y = t sigma, the first one reads t (A - L_1) sigma = -psi sigma - psi t dsigma/dt; divide by
-psi t. The division by psi is legitimate because (A - L_1) sigma vanishes at psi = 0: sigma(0) = rho v_1(0) = rho r_{0,1}[0] chi_1,
since r_{1,1}[0] = L_1 r_{0,1}[0] by (H), and (A - L_1) chi_1 = 0. QED
Sign: sigma is built from v_1(-psi), so z = -psi and the 1/z of the QDE becomes -1/psi. Only the first form is used below.

### 3.3 The flow

Theorem 7 (E_11 flow; V8). In K[[t]][[psi, psi']],

    d/dt Ebar(psi, psi') = Theta(t) Dg(psi) Dg(psi'),     Dg = d/dt (t g),                                              (13)
    Theta := (L_2 - L_1)(x^2 - c x L_2)^2/Delta_2 = x^2 - c x L_2 = Delta_2/(L_2 - L_1) = -x^3 tau s*^2 = x c d(dphi)/dt - x.     (14)

Proof. (i) The role of rho^2 Delta_1 = -x c^3/t^2. Since t^2 rho^2 Delta_1 is a constant ([S2.4c]), only the pairing term of Nbar
depends on t, and dNbar/dt = -<dY(psi)/dt, Y(psi')> - <Y(psi), dY(psi')/dt>. Multiply by psi psi', and use Lemma 6 and
self-adjointness (2.4):
    psi psi' dNbar/dt = psi' <(A - L_1) Y(psi), Y(psi')> + psi <Y(psi), (A - L_1) Y(psi')> = (psi + psi') <(A - L_1) Y(psi), Y(psi')>.
(ii) psi + psi' does not depend on t, so (psi + psi') dEbar/dt = dNbar/dt and
    (psi + psi') psi psi' dEbar/dt = (psi + psi') <(A - L_1) Y(psi), Y(psi')>.
Cancel psi + psi'. This is why the factor (psi + psi') disappears: the flows of the two arguments contribute 1/psi and 1/psi',
self-adjointness makes the two terms equal, and 1/psi + 1/psi' = (psi + psi')/(psi psi').
(iii) By (10) and Proposition 5, pi(Y(psi)) = t pi(sigma(psi)) = -(x^2 - c x L_2) psi Dg(psi), so
    psi psi' dEbar/dt = ((L_2 - L_1)/Delta_2) (x^2 - c x L_2)^2 psi psi' Dg(psi) Dg(psi').
Cancel psi psi'. This is (13) with Theta = (L_2 - L_1)(x^2 - c x L_2)^2/Delta_2, the definition in V8.
(iv) Simplification. By (9), Delta_2 = (L_2 - L_1)(x^2 - c x L_2), hence Theta = x^2 - c x L_2 = Delta_2/(L_2 - L_1). The last two
forms are S2.11 (I3). First, x^2 - c x L_2 = -x^3 tau s*^2 (2.4). Second, d(dphi)/dt = 1 - x s* = 1 - L_2 (S2 Lemma S2.2(e);
S1 s.7.4), so x c d(dphi)/dt - x = x c - x - c x L_2 = x^2 - c x L_2. QED
Expansion (tau = t - t^2/2 + ..., s* = 1/c + x tau/c^3 + ...):
    Theta = -x^3 t/c^2 + x^3 (1-x)^2 t^2/(2 c^4) - x^3 (1 - 8x + 12x^2 - 8x^3 + x^4) t^3/(6 c^6) + O(t^4),    Theta_0 = 0,
and Theta_k = (k+1) x c [t^(k+1)] dphi for k >= 1. This reproduces the printed Theta(t) and dphi(t) of
hvc_leaf_proof_check_out.txt through t^7 (4.2).
Remarks. (a) Theta = Delta_2/(L_2 - L_1): the loop flow is governed by the regular idempotent. Delta_2 is regular and non-zero,
while L_2 - L_1 has a simple pole, so Theta = O(t).
(b) Rank one. Before expanding in t, the right side Theta(t) Dg(psi; t) Dg(psi'; t) of (13) is a product of a function of
(t, psi) and a function of (t, psi'). It has rank one as a function of t because A - L_1 has rank one (S2 Proposition S2.8(e)).
Its t^k-coefficients are not of rank one in general. They are sums over k0 + k1 + k2 = k of rank-one terms
Theta_k0 (k1+1)(k2+1) g_k1(psi) g_k2(psi'), as in (16). For example, the t^1-coefficient of rho^2 E_11 has e^(1)_11 and
e^(1)_12 = e^(1)_21 non-zero and e^(1)_22 = 0 (4.1), a 2 x 2 matrix of rank 2. What rank one gives is that every loop coefficient is
built from products of the same leaf series g_k.

### 3.4 Integration from t = 0; F3 for E_11

By Proposition 4(a), [t^0] Y = (0, c), and <(0, c), (0, c)> = -c x c^2 = -x c^3. Hence [t^0] Nbar = -x c^3 + x c^3 = 0, and, as
psi + psi' is t-independent, [t^0] Ebar = 0. This is where K_1 = 1 enters. For a t-independent K_1(w) one would have
P(-w, w) = -K_1(w) Delta_1, and [t^0] Ebar = -x c^3 (K_1(psi') - 1)/(psi + psi') (the outline's integration constant).
(This constant also vanishes for a second reason. The numerator of the difference quotient vanishes at psi = -psi' by construction,
while its t^0-coefficient after multiplication by t^2 rho^2, namely -x c^3 (K_1(psi') - 1), does not depend on psi. This is S2's
proof of K_1 = 1, S2 Proposition S2.9(c), seen through (11).)
The kernel of d/dt on K[[t]] is K, and [t^(n+1)] Ebar = [t^n](dEbar/dt)/(n+1). So (13) and [t^0] Ebar = 0 determine Ebar:

    Ebar = Int_0^t Theta(u) Dg(psi; u) Dg(psi'; u) du  (termwise),   i.e.   [t^(n+1)] Ebar = [t^n](Theta Dg(psi) Dg(psi'))/(n+1).      (15)

As Theta_0 = 0, also [t^1] Ebar = 0. So Ebar is in t^2 K[[t]][[psi, psi']], and E_11 = (t rho)^-2 Ebar is in K[[t]][[psi, psi']],
since rho^-2 is in K[[t]] ([S2.4c]). So E_11 is regular at t = 0: this is F3 for E_11. (In the point frame and in frame A the
same edge factor has a double pole, conifold_regularity_note.md F3. In frame B the double pole of Delta_1 is exactly compensated.)

### 3.5 The leaf structure of the loops (L3)

Theorem L3. For every k >= 0,

    [t^k](rho^2 E_11)(psi, psi') = (1/(k+2)) Sum_{k0 >= 1, k1, k2 >= 0, k0 + k1 + k2 = k+1} Theta_{k0} (k1+1)(k2+1) g_{k1}(psi) g_{k2}(psi')    (16)
                                = Sum_{i, l >= 1, (i-1)+(l-1) <= k} e^(k)_{il} lam_i(psi) lam_l(psi'),
    e^(k)_{il} = (1/(k+2)) Sum_{k0 >= 1, k1 + k2 = k+1-k0} Theta_{k0} (k1+1)(k2+1) g_i^(k1) g_l^(k2),       e^(k)_{il} = e^(k)_{li},

where g_i^(m) := 0 for i > m+1. There is no lam_0 in either variable. The support bound is sharp for generic x. For
(i-1) + (l-1) = k, i.e. i = k1 + 1 and l = k2 + 1 with k1 + k2 = k,

    e^(k)_{k1+1, k2+1} = (-1)^(k+1) x^3 (k1+1)(k2+1) P_{k1} P_{k2} / ((k+2) (k1!)^2 (k2!)^2 c^(2k)),             (17)

which is non-zero whenever x P_{k1}(x) P_{k2}(x) != 0. At t = 0 the loop has rank one:
[t^0](rho^2 E_11) = (1/2) Theta_1 g_0(psi) g_0(psi') = -(x^3/2) lam_1(psi) lam_1(psi'). As rho(0) = 1, this gives
E_11(z, w)|_{t=0} = -(x^3/2) phi(z) phi(w), phi(z) = 1/(c - z), which is the E_11 part of F3 of conifold_regularity_note.md.
Proof. [t^k](rho^2 E_11) = [t^(k+2)] Ebar = [t^(k+1)](Theta Dg(psi) Dg(psi'))/(k+2) by (15). Inserting
Theta = Sum_{k0>=1} Theta_{k0} t^k0 and Dg = Sum_k (k+1) g_k t^k gives (16).
Support. Every term has k1 + k2 = k + 1 - k0 <= k. By Theorem L1, g_{k1} is in Lam_{k1+1} and g_{k2} in Lam_{k2+1}, so the term
lies in span{lam_i lam_l' : i <= k1+1, l <= k2+1}, whose total degree is (i-1) + (l-1) <= k1 + k2 <= k.
Symmetry. (16) is invariant under (psi, k1) <-> (psi', k2).
Corner. A product lam_{k1+1}(psi) lam_{k2+1}(psi') with k1 + k2 = k comes only from terms with k1' >= k1, k2' >= k2 and
k1' + k2' <= k, i.e. k1' = k1, k2' = k2 and k0 = 1. Its coefficient is (1/(k+2)) Theta_1 (k1+1)(k2+1) g_{k1+1}^(k1) g_{k2+1}^(k2).
Insert Theta_1 = -x^3/c^2 and (5). QED
Generating form. The lam_i-coefficient of Dg is d/dt(-(i/c) t A_i) = -(i/c) B_i', with B_i := t A_i (1.6). So
e_il(t) := Sum_k e^(k)_{il} t^k = (i l/(c^2 t^2)) Int_0^t Theta(u) B_i'(u) B_l'(u) du; this is S4 Theorem CF.4. The case k1 = 0 of
(17) is the corner e^(k)_{1,k+1} = ((k+1)^2/(k+2)) (x^3/c^2) a_{k+1}^(k) of leaf_structure.md, Step 3 (same symbol a_{k+1}^(k)
there; S4 Corollary CF.4').

### 3.6 rho versus rho^2; L2 for j = 1

In Theorem R (proof_status_mixed_scale.md s.3), every half-edge at a singular vertex of genus h >= 2 carries one factor rho. By the
dilaton equation F_h(lambda q) = lambda^(2-2h) F_h(q), one has
d^n F_h(q_t)[f_1, ..., f_n] = rho^(2h-2) d^n F_h(rho q_t)[rho f_1, ..., rho f_n], and rho^(2h-2) Delta_1^(h-1) = (-x c^3/t^2)^(h-1) by
[S2.4c] (Theorem R(i)). So:
- a loop, or an edge between two such vertices, carries rho^2 E_11;
- a leg, i.e. an E_11-edge to a vertex that is not rescaled (singular of genus 0 or 1), with external power b at the far end,
  carries rho f_{1,b};
- an E_11-edge between two non-rescaled vertices carries E_11.
All of them have the leaf structure:

Corollary 8. For every integer a and every k >= 0, [t^k](rho^a E_11) lies in span{lam_i(psi) lam_l(psi') : i, l >= 1,
(i-1) + (l-1) <= k}. In particular (L2 for j = 1), for every b >= 0,

    [t^k](rho f_{1,b}) = Sum_{i=1}^{k+1} ( Sum_{m=0}^{k} [t^m](rho^-1) Sum_{l} e^(k-m)_{il} l^b/c^(b+1) ) lam_i(psi)   in Lam_{k+1}.

Proof. rho^(a-2) is in K[[t]] by [S2.4c], with psi-independent coefficients. So [t^k](rho^a E_11) = Sum_m [t^m](rho^(a-2))
[t^(k-m)](rho^2 E_11), and the supports at the orders k - m <= k are contained in the support at order k. For the leg,
f_{1,b}(psi) = Sum_a E_11(a, b) psi^a = [psi'^b] E_11(psi, psi') and rho f_{1,b} = [psi'^b](rho E_11) (a = 1). By (LAM1),
[psi'^b] lam_l(psi') = l^b/c^(b+1), and the lam_i(psi) that occur have i - 1 <= k - (l - 1) <= k. QED
Explicitly, with rho^-1 = 1 + ((1-x)^2/(2c^2)) t + ((1 + 4x - 12x^2 + 4x^3 + x^4)/(12 c^4)) t^2 + O(t^3):

    [t^0](rho f_{1,b}) = -(x^3/2) c^(-b-1) lam_1,
    [t^1](rho f_{1,b}) = ((1-x)^2 x^3/(12 c^(b+3))) [ (2^(b+3) - 9) lam_1 + 8 lam_2 ].

At t^0 the b-dependence is only through c^-b. This is the rank-one observation of conifold_regularity_note.md s.6, step 3(i).


## 4. Checks

### 4.1 By hand, against the tables of leaf_structure.md

In this subsection a prime on a leaf series means its value at psi' (lam_1' = lam_1(psi'), g_1' = g_1(psi')), not a derivative.

Order t^0.
- L1: R_0 = c - c^2/(c+z), i.e. [t^0](rho q_t) = c - c^2 lam_1 = q_0.
- H-row: [t^-1] sigma = (0, c); sigma^(0) = (c lam_1, x lam_1) (F2).
- Theta_1 = -x^3/c^2 (printed Theta(t)).
- Loops: (1/2) Theta_1 c^2 lam_1 lam_1' = -(x^3/2) lam_1 lam_1'; table (iii), k = 0: e_11 = -x^3/2.
  Independently of the flow, from (12):
      [t^2] <Y(psi), Y(psi')> = x^2 [g_0 sigma^(0)_1' + sigma^(0)_1 g_0' + c (g_1 + g_1')] - c x [c (sigma^(1)_1 + sigma^(1)_1')
                                 + sigma^(0)_1 sigma^(0)_1']
                              = c x^3 lam_1 lam_1' + c x [x (g_1 + g_1') - c (sigma^(1)_1 + sigma^(1)_1')],
  and x g_1 - c sigma^(1)_1 = (x/(2c)) [(1-x)^2 - (1 - x + x^2)] lam_1 = -(x^2/(2c)) lam_1. Using
  lam_1 + lam_1' = (2c - psi - psi') lam_1 lam_1', this gives
      [t^2] <Y, Y'> = x^3 lam_1 lam_1' (c - (2c - psi - psi')/2) = (x^3/2)(psi + psi') lam_1 lam_1',
  hence [t^2] Ebar = -(x^3/2) lam_1 lam_1'.
- Legs: -(x^3/2) c^(-b-1) lam_1; table (ii), j = 1, t^0: -(1/2) x^3/(1+x)^(b+1), b = 0..7.

Order t^1.
- L1 from the closed form. tau/t = 1 - t/2 + O(t^2) and H_1 = F_1 + X_1 = (1+z)^2/(z (c+2z)) - 1/(c z) = (2x + c z)/(c (c+2z)), so
      G_1 = (c/(c+z)) (H_1 - 1/2) = (c/(c+z)) (4x - c^2)/(2c (c+2z)) = -(1-x)^2/(2 (c+z)(c+2z))      (c^2 - 4x = (1-x)^2),
  and R_1 = ((1-x)^2/2) (1/(c+2z) - 1/(c+z)). So delta_1 = ((1-x)^2/2)(lam_2 - lam_1), as in table (i). (No pole at z = 0, limit 0
  at oo, a_1^(1) + a_2^(1) = 0.) Also G_1 - rho_1 G_0 = (1-x)^2 z/(c (c+z)(c+2z)) with rho_1 = -(1-x)^2/(2c^2). This is the O(t)
  row correction [t^1] r_{0,1} of leaf_structure.md, Step 0, and of S1 s.7.5.
- g_1 = ((1-x)^2/(2c)) (lam_1 - 2 lam_2) = -((1-x)^2/2) lam_1 lam_2, which is O(psi^-2).
- Theta_2 = x^3 (1-x)^2/(2c^4) (printed Theta(t)). Also Theta_2 = 3 x c [t^3] dphi = 3 x c (1-x)^2 x^2/(6 c^5).
- H-row and symplecticity, with psi' = -psi:
      [t^1] <Y(psi), Y(-psi)> = x^2 c (g_0(psi) + g_0(-psi)) - c^2 x (sigma^(0)_1(psi) + sigma^(0)_1(-psi)) = 0,
      [t^2] <Y(psi), Y(-psi)> = c x^3 lam_1(psi) lam_1(-psi) + c x [x (g_1(psi) + g_1(-psi)) - c (sigma^(1)_1(psi) + sigma^(1)_1(-psi))] = 0.
  For the second, use lam_j(psi) + lam_j(-psi) = 2c/(c^2 - j^2 psi^2). The coefficient of 1/(c^2 - psi^2) is
  x^2 + x (1-x)^2 - x (1 - x + x^2) = 0, and that of 1/(c^2 - 4 psi^2) is -2x (1-x)^2 + 2x (1-x)^2 = 0. This tests (8) at k = 1,
  including l_0 and l_1, against K_1 = 1 (<Y(psi), Y(-psi)> = t^2 rho^2 Delta_1 = -x c^3).
- Loops, (16) at k = 1:
      (1/3) [Theta_2 g_0 g_0' + 2 Theta_1 (g_1 g_0' + g_0 g_1')] = ((1-x)^2 x^3/c^2) [-(1/2) lam_1 lam_1' + (2/3)(lam_1 lam_2' + lam_2 lam_1')].
  So e_11 = -(1-x)^2 x^3/(2c^2), e_12 = e_21 = 2 (1-x)^2 x^3/(3c^2) and e_22 = 0, as in table (iii), k = 1. (e_22 = 0 is forced by
  the support bound.)
- Legs. [t^1](rho f_{1,b}) = [t^1](rho^2 f_{1,b}) + rho^-1_1 [t^0](rho^2 f_{1,b}), where rho^-1_1 := [t^1](rho^-1) = (1-x)^2/(2c^2),
  [t^0](rho^2 f_{1,b}) = -(x^3/2) c^(-b-1) lam_1, and [t^1](rho^2 f_{1,b}) = c^(-b-1) [(e_11 + 2^b e_12) lam_1 + e_21 lam_2] with the
  k = 1 entries. The lam_1-coefficient is
      ((1-x)^2 x^3/c^(b+3)) (-1/2 + (2/3) 2^b - 1/4) = ((1-x)^2 x^3/(12 c^(b+3))) (2^(b+3) - 9),
  i.e. -1/12, 7/12, 23/12, 55/12, 119/12, 247/12, 503/12 times (1-x)^2 x^3/c^(b+3) for b = 0..6. The lam_2-coefficient is
  (2/3) (1-x)^2 x^3/c^(b+3). This is table (ii), j = 1, t^1, for every b with data.

Order t^2.
- L1 from the closed form. [t^2]((tau/t) H) = H_2 - H_1 + 1/6, where H_2 = F_2 + F_1 X_1 + X_2, X_1 = -1/(c z) and
  X_2 = -(1 + 2x)/(2 c^3 z) + 1/(2 c^2 z^2) (from dphi = tau/c + (1+2x) tau^2/(2c^3) + O(tau^3)).
  Regularity at 0. The z^-2 parts of F_2, F_1 X_1, X_2 are 1/(2c^2), -1/c^2, 1/(2c^2), and the z^-1 parts are (6 - 5/c)/(2c^2),
  -(2 - 2/c)/c^2, -(1+2x)/(2c^3). Both sums vanish, so H_2 (and H_1) are regular at z = 0, as S1 asserts.
  Residues. a_3^(2) = (5) at k = 2 = -(2-x)^2 (1-2x)^2/(12c^2). Next,
      a_2^(2) = [(c z/(c+z)) (c + 2z)(H_2 - H_1)]_{z=-c/2} = (-c) [-(1-x)^2 x^2/c^3 - (1-x)^2/c^3 + (1-x)^2/(2c)] = (1-x)^4/(2c^2).
  Finally a_1^(2) = [c z (H_2 - H_1 + 1/6)]_{z=-c}, with H_1(-c) = (1 + x^2)/c^2, F_2(-c) = x^2 (1+2x)^2/(4c^4),
  (F_1 X_1)(-c) = x^2/c^4 and X_2(-c) = 1/c^3. So 12 c^4 (H_2 - H_1 + 1/6)(-c) = 2 - 4x + 3x^2 - 4x^3 + 2x^4, and
  a_1^(2) = -(2 - 4x + 3x^2 - 4x^3 + 2x^4)/(12c^2).
  These are the three entries of table (i), k = 2, and a_1^(2) + a_2^(2) + a_3^(2) = 0.
- g_2 = (1/(12 c^3)) [(2 - 4x + 3x^2 - 4x^3 + 2x^4) lam_1 - 12 (1-x)^4 lam_2 + 3 (2-x)^2 (1-2x)^2 lam_3].
- Theta_3 = -x^3 (1 - 8x + 12x^2 - 8x^3 + x^4)/(6c^6) = 4 x c [t^4] dphi.
- Loops, (16) at k = 2:
      (1/4) [Theta_3 g_0 g_0' + 2 Theta_2 (g_1 g_0' + g_0 g_1') + 3 Theta_1 (g_2 g_0' + g_0 g_2') + 4 Theta_1 g_1 g_1'].
  The coefficients are:
      lam_1 lam_1': (x^3/(4c^4)) [-(1 - 8x + 12x^2 - 8x^3 + x^4)/6 + (1-x)^4 - (2 - 4x + 3x^2 - 4x^3 + 2x^4)/2 - (1-x)^4]
                    = -x^3 (7 - 20x + 21x^2 - 20x^3 + 7x^4)/(24 c^4);
      lam_1 lam_2': (1/4) (1-x)^4 x^3 (-1 + 3 + 2)/c^4 = (1-x)^4 x^3/c^4;
      lam_1 lam_3': -3 x^3 (2-x)^2 (1-2x)^2/(16 c^4);
      lam_2 lam_2': (1/4) 4 Theta_1 (1-x)^4/c^2 = -(1-x)^4 x^3/c^4;
      lam_2 lam_3', lam_3 lam_3': 0.
  These are all the entries of table (iii), k = 2.
- Legs, k = 2, b = 0: lam_1: x^4 (4 - 3x + 4x^2)/(48 c^5); lam_2: (1-x)^4 x^3/(3 c^5); lam_3: -3 (2-x)^2 (1-2x)^2 x^3/(16 c^5).
  These come from (1/c)(e_11 + e_12 + e_13, e_21 + e_22, e_31) + rho^-1_1 [t^1] + rho^-1_2 [t^0], and they match table (ii), j = 1,
  b = 0, t^2.

### 4.2 Exact rational checks (S3_checks.pl)

leaf_proof/S3_checks.pl (Perl + Math::BigRat; no Mathematica) parses the stored Mathematica tables hvc_leaf_structure_out_kt5.txt,
hvc_leaf_structure_out.txt and hvc_leaf_proof_check_out.txt, and evaluates them at x = 2/7, 3/5, 7/3, -5/9. It runs in about 10 s
and checks the following in exact rational arithmetic, with t-series truncated at t^10 and compared only where complete. Result:
1540 checks, 0 failures (output S3_checks_out.txt, with the count per type).
- Sections 0.3, 2.4, 3.3: the critical equation; L_1 + L_2 = c/tau; L_1 L_2 = x/tau; (x - c L_1)(x - c L_2) = x^2;
  x^2 - c x L_2 = -x^3 tau s*^2; (1 - tau) s*/(1 - tau s*) = 1 - x s*; Theta = x c d(dphi)/dt - x; Theta_0, Theta_1, Theta_2; the
  printed Theta(t) and dphi(t) (t^1..t^7) at nine values of t; rho = tau/(t alg) against the stored rho through t^5; l_0.
- Section 1: the closed form G, as a t-series at the rational points z = 3/11, -2/19, 5, -13/4, against the dilaton tables
  t^0..t^5 (z G_k = Sum_j a_j^(k)/(c + jz)); G_1; Sum_j a_j^(k) = 0 and g_k = O(psi^-2) for k = 1..5; no lam_0; the top
  coefficients (5) for k = 1..5; the residue form and the hypergeometric form of A_j(t), j = 1..6, t^0..t^5.
- Section 2: [t^0] Y = [t^-1] sigma = (0, c); sigma^(0)_1 = x lam_1; sigma^(1) (both components); the symplectic identity
  <Y(psi), Y(-psi)> = -x c^3 through t^5, with sigma built from (8) out of the tabulated g_k and l(t), at psi = 1/13, -2/17, 3/19;
  Proposition 5 through t^5.
- Section 3: (16) against the restricted loop tables t^0..t^5 (every entry, the support, symmetry, all corners (17)) and against
  the full loop matrices t^0..t^3 (including all their zeros); Ebar = O(t^2); [t^k](rho^2 E_11) computed from (12), with sigma
  built from the data, equal to (16) for k <= 3 at three pairs (psi, psi'); Corollary 8 against every tabulated j = 1 leg
  (b = 0..7, t^0..t^5 where the table has an entry: 33 entries per value of x) and the closed t^0, t^1 leg formulas.
Remark on the leg tables. Some entries are fitted from exactly k + 2 psi-orders ("checks 0"). They are still a genuine test of
Corollary 8 at the psi-orders 0..k+1, since the fit is a bijection onto the coefficients, but not of the leaf form beyond them.

### 4.3 Mathematica checks requested (not run; licence)

 M1. From the closed form G alone, with no stored factors: build v_1 = rho^-1 (G, (L_1 + z/t) G + z dG/dt), compute E_11 by the
     evalQ difference quotient for m + n <= 14 and t^0..t^8, and compare with (16). This extends V8 beyond the stored z^9 range.
 M2. From the same closed-form rows: [t^k](rho f_{1,b}) for b <= 12 and k <= 4, against Corollary 8. The stored leg tables stop at
     b <= 7 and have few psi-orders at large b.
 M3. The checks (iii)-(iv) of S1 s.8, on which 1.3(b) rests. Run the transport recursion of S1 s.7.4 for m <= 5 against
     [z^m](Exp[-dphi/z] F) to tau^10, symbolically in x, including S1's Phi_0 = 1/alg.


## 5. Logical summary

Inputs (0.4): [S1.d] and [S2.4c] (hence (KF)), [S2.7], [S2.9c], the definition (H), and the algebraic facts of S2 recalled in 2.1,
2.2, 2.4, 2.5, 3.1 and 3.2. Proved in this section:
- (PF), (T), Lemmas 1 and 2: unconditional. They use only dphi in tau K[[tau]] and the definition of F_n.
- Theorem L1 for all k (needs (KF)): g_k in Lam_{k+1} and g_k = O(psi^-2) for k >= 1. The residue / hypergeometric form of A_j(t)
  (= S4 CF.1) and the top coefficient (5) need only the closed form G and Proposition 3.
- Proposition 4 (needs (KF), (H), L_1 = c/t + l, and (6), i.e. [S2.4c] or just F1): sigma = t^-1 (0, c) + Sum_k t^k sigma^(k), with
  sigma^(k) in Lam_{k+1}^2. Part (a) is S2 Proposition S2.9(b), and it does not use [S2.9c]. Proposition 4 is the input of Step 9
  (S4 Lemma 9.1).
- (9), from S2 Proposition S2.3(c). Proposition 5 is S2.11 (I2), and (10) is S2 Proposition S2.8(e).
- Theorem 7: the E_11 flow (13), with Theta as in (14). It needs [S2.4c]; Lemma 6 (= S2 Proposition S2.9(d), i.e. [S2.7]);
  self-adjointness (S2 Proposition S2.8(c)); (10) and Proposition 5; and [S2.9c] through the form (12) of E_11. The
  t-independence of K_1 alone would suffice for the flow.
- 3.4 (needs [S2.9c] and Proposition 4(a)): Ebar(t = 0) = 0 and Ebar = O(t^2); F3 for E_11.
- Theorem L3 for all k: the explicit formula (16), the support bound (i-1) + (l-1) <= k, no lam_0, symmetry, sharpness (17) for
  generic x, and rank one at t^0 (F3 of the note).
- Corollary 8: L2 for j = 1, for all b and all k; [t^k](rho^a E_11) for every integer a.
- F3 for Dil_1. Dil_1(d) = (-1)^d [z^(d-1)] r_{0,1} = (-1)^d [z^(d-1)] (rho^-1 G) is in K[[t]], by (KF) and rho^-1 in K[[t]]
  ([S2.4c]). Directly: r^pt_{0,1} is in C{tau}[[z]] by S1 Theorem S1(b), hence so is r_{0,1} = h_w^-1 r^pt_{0,1}, since h_w is
  t-independent. The value at t = 0 is Dil_1(d; 0) = -c^(1-d), i.e.
  q_t|_{t=0} = q_0 = -c psi/(c - psi).
Where the inputs are used: (KF) in 1.3-1.5 and 2.1 (and through L1 everywhere); [S2.4c] in 0.4, 2.1, 3.1, 3.3(i), 3.4 and 3.6; (H)
in 2.1; [S2.7] in 3.2; [S2.9c] in 3.1 and 3.4.

Where F1-F3 are proved (the non-leaf input of Theorem R, proof_status_mixed_scale.md s.3):

    F1   Delta_1 = -x c^3 t^-2 (1 + O(t))          S2 Proposition S2.4(c) (= [S2.4c]); also (9) in 2.4
         Delta_2 regular, Delta_2(0) = x^3/c        S2 Proposition S2.3(e); S4 10.3
    F2   r_{0,1}(z; 0) = c/(c + z)                  S1 Corollary 7.3; S2 Proposition S2.9(b)
    F3   Dil_1                                      S1 Theorem S1(b), or (KF) + [S2.4c] (this section, s.5)
         E_11                                       this section, 3.4
         E_12                                       S4 Proposition 9.3(i)
         E_21                                       S4 10.6(b) and Proposition 9.5 (uses <v_1(w), v_2(-w)> = 0: S2 Proposition S2.10 with S4 Lemma 10.5)
         Dil_2, E_22                                S4 10.4 (Step 10)

What S4 takes from S3: S4 (I2) is Proposition 3 / Theorem L1; S4 (I4) is (7) with Proposition 4 (S4 Lemma 9.1 restates (8) with
explicit coefficient vectors); S4 (I6) is Theorem 7 with 3.4.
