# S2. Metric, eigenvalues, quantum differential equation, symplectic constant (Steps 3, 5, 6 of outline.md)

2026-10-31, revised the same day after the referee report. Section S2 of the all-order proof of Lemma C (leaf_proof/outline.md). No topological
recursion and no spectral curves. Conventions: outline.md s.0, i.e. rc_lib.m and hvc_lib.m at alpha0 = 0, alpha1 = eps1 = x, alpha2 = eps2 = 1,
mode B (s = -1). "rc N" and "hvc N" are line numbers in rc_lib.m and hvc_lib.m.
Checks. Every identity used downstream (in S3, S4 and Theorem R) has an exact rational spot check (Perl + Math::BigRat, no Mathematica):
 leaf_proof/S2_checks.pl for S2.1-S2.9 and S2.11 (output S2_checks_out.txt: 943 checks, 0 failures, series to order 22), and
 leaf_proof/S2_checks_reg.pl for the regular idempotent, S2.10 (output S2_checks_reg_out.txt: 420 checks, 0 failures; z^6, jets to order 4, t^8).
 The statements without a check of their own are immediate from the definitions: Lemma S2.1(d) (total symmetry of third derivatives), the module
 formalism (5.1), (5.2.1) (a logarithmic derivative) and (10.1), and the recursion in the proof of Prop S2.10 (only its conclusion is checked).

Changes in this revision. (1) The series e^{-dphi/z} is called Ephi; the letter g now always means g(psi) = G(-psi). (2) S2.10 now proves the second
QDE row at the regular idempotent (Lemmas S2.10.1-S2.10.3, Prop S2.10.4), with Gauss's equation formulated in a module (d/dtau on E_lambda, (5.1));
hence the orthogonality (Prop S2.10) is unconditional, and K_2 = 1 holds to all orders given S4 Prop. 10.6(a) (Prop S2.10.5). (3) Smaller corrections
in Prop S2.4(b), (c), Remark S2.7.4 and the header; the logical order S1 -> S2 -> S3 is stated explicitly.

Proved here, to all orders in z:
 - Step 3: Delta_1 = x^2 W2(U_1) = x^2 (1 - tau s*)^2 phi''_tau(s*)/tau^2 = -x c^3 alg^2/tau^2, so rho = tau/(t alg), and with the closed form of S1
   the KEY FORMULA (V3); the eigenvalue and metric identities (V1), (V2); F1 for Delta_1 and Delta_2.
 - Step 5: du_1/dt = L_1; the identity (V7), by hand, as an identity of differential operators; the second row of the QDE at the singular
   idempotent, from (V7) and the closed form of S1; the operator A (eigen-covectors, self-adjointness, normalisation, spectral decomposition).
 - Step 6: <v_1(w), v_1(-w)> = Delta_1 (K_1 = 1), from t-independence and the t -> 0 limit.
 - The regular idempotent (S2.10): the second row of the QDE at U_2 ((V5) for i = 2), from S1's formal Laplace calculus applied at the regular
   critical point s~ for fixed t != 0, without expansion at tau = 0; hence <v_1(w), v_2(-w)> = 0 (Prop S2.10), v_2(z; 0) = r_{0,2}(z; 0) (1, x/c), and,
   with S4 Prop. 10.6(a), <v_2(w), v_2(-w)> = Delta_2 (K_2 = 1). So (V6) holds in full, to all orders. The orthogonality and Row 2 at U_2 are not
   needed for L1-L3 or Step 8; they are needed for the transposition E_21(b, a) = E_12(a, b) and the regularity of E_21 (S4, Prop. 9.5), hence for
   Theorem R (ii) (graphValueC attaches EdN[c(v1), c(v2)] in the order of the edge list, so E_21 occurs whenever the regular vertex is listed first).
Input from other sections. (i) The closed form of the singular identity row (S1, Theorem S1(b)-(d)), in the form (S1.a), (S1.b) of S2.0. (ii) For S2.10
 only: S1's formal Laplace calculus, i.e. the definitions 2.2 and (D1)-(D3), Props. 2.2-2.4 ((P3), (P2), (P1)), which S1 states for an arbitrary
 nondegenerate critical point depending analytically on a parameter, the change of variables (1.1) and the amplitude identity (4.1). (iii) For K_2 = 1
 only (Prop S2.10.5(c), not used downstream): the value r_{0,2}(z; 0) = h_w(z)^-2 of S4, Prop. 10.6(a).
Logical order: S1 -> S2 -> S3, S4. S2 uses nothing from S3 or S4 except the input (iii), whose proof (S4 10.6(a): Laplace's method for a Beta integral
 at t = 0 and Stirling's series) uses nothing from S2. In particular the value [t^-1] sigma = (0, c) that enters the proof of K_1 = 1 is derived here,
 in Prop S2.9(b), from (H) and (S1.b). (S3, s.0.4 attributes this value to its Prop. 4(a), an independent derivation of the same fact; there is no
 circularity in either direction.)
Notation. The vector s(z) = (s_0, s_1) of outline.md Steps 5-6 (idempotent 1) is S(z) = (S_0, S_1) here (S2.5, S2.7), and S^(2) at idempotent 2 (S2.10);
 s* always denotes the critical point of phi_tau with s*(0) = 1/c (Euler coordinate s = (1 - U)/tau), s~ the other one, never a component.
 Ephi := e^{-dphi/z} (S2.0 (0.6)); the letter g is reserved for g(psi) := G(-psi) (Prop S2.9(b), S2.11 and outline (V8)). The amplitude (1 - tau s)^-1,
 which S1 calls g, is q_0 here (S2.10). u_i := W(U_i) (Lemma S2.2(e), S2.10.1(c)); f2r := f2/tau^2 at the regular critical point (Lemma S2.10.1(b)).
 kappa_1(z) (S2.5) is the t-independent constant of the outline's "s_0 = kappa(z) Q^{x/z} tau^{c/z+1} F". The phase constant -e^{x Log[-1]/z} of S1 (1.2),
 which S1 also calls kappa(z), is written kappa_ph here; it plays no role in S2.


## S2.0 Standing assumptions, rings, input

(0.1) x is a parameter with x != 0, -1, and c = 1 + x. All statements are identities of rational or analytic functions of x, so they hold for symbolic x.
 Q = e^t, tau = 1 - 1/Q = 1 - e^-t, and ' := d/dt = (1 - tau) d/dtau (dtau/dt = 1 - tau). The quantum product is written a * b.

(0.2) Germs. tau(t) = t - t^2/2 + ... is biholomorphic near t = 0, so O := C{t} = C{tau} (convergent power series) and its fraction field
 M := C({t}) = C({tau}) (germs of meromorphic functions) do not depend on the choice of variable; tau/t and t/tau are units of O.
   delta := c (1 - 4 x tau/c^2)^(1/2) in O (principal branch, delta(0) = c), so delta^2 = c^2 - 4 x tau;
   s* := 2/(c + delta) in O, s*(0) = 1/c;   s~ := 1/(x tau s*) in tau^-1 O.
 Since (c - delta)(c + delta) = 4 x tau, s* = (c - delta)/(2 x tau), i.e. 2 x tau s* = c - delta; s* and s~ are the two roots of
 x tau s^2 - c s + 1 (sum c/(x tau), product 1/(x tau)).

(0.3) Formal rings. K := C((z)); R := K((tau)) with the K-linear derivation d/dt = (1 - tau) d/dtau; R_0 := C[[z]]((tau)) = C[[z, tau]][1/tau].
 R_0 is a subring of R and (expanding in z) of C((tau))[[z]], compatibly with products and with d/dt. An element of M[[z]] whose z-coefficients
 have poles of bounded order at tau = 0 is an element of R_0.
 TRANSFER. An identity in R between elements of R_0 is an identity of formal double series: reading off z^m gives identities of formal Laurent
 series in tau, and when these are expansions of meromorphic germs (as for all rows below) identities of meromorphic functions of t near 0.
 "To all orders in z" always means this.
 The z-adic ring M((z)) (formal Laurent series in z, bounded below, with coefficients in M; d/dt acts coefficientwise) is used at the regular idempotent
 (S2.10), where no expansion at tau = 0 is needed; there "to all orders in z" means: coefficientwise in z, as identities in M.
 refl := the automorphism z -> -z of K, extended to R coefficientwise in tau (and to M((z)) coefficientwise in z); it commutes with d/dt.
 For f in R or M((z)), f(-z) := refl(f).
 (S1 uses "iota" for its Taylor-expansion map C{tau}[[z]] -> C[[z, tau]]; here that map is left implicit, as the inclusion O[[z]] in C[[z, tau]].)

(0.4) Rows and pairing. r_{a,i}(z) = Sum_{m>=0} r_{a,i}[m] z^m (a = 0 identity row, a = 1 H-row; i = 1, 2 idempotent), in frame B. The point-frame
 rows are r^pt_{a,i} = h_w(z) r_{a,i}: buildFactors (hvc 151-154) runs rmatrixAt0 and then twistRows (hvc 146-149), which multiplies every row
 by the t-independent series hCoeffs[w, -1] = h_w(z)^-1, w = (c, -x, -1) (hvc 52-55). The coefficients r_{a,i}[m] are elements of the quadratic algebra
 C(Q)[U]/(LeqU) evaluated at U = U_i (rc 121-132), i.e. rational functions of (tau, delta); they lie in M.
 v_i(z) := (r_{0,i}(z), r_{1,i}(z)) is a covector (its components pair with the flat basis 1, H), and
   <a, b> := Sum eta^{a'b'} a_{a'} b_{b'} = x^2 (a_0 b_1 + a_1 b_0) - c x a_1 b_1,
 the pairing of evalQ (rc 147) and of check (5) of hvc_leaf_proof_check.wls. Its matrix eta^-1 is constant in t (Lemma S2.1).

(0.5) H-row, the definition (rc 130-132; cDerN is the total Q-derivative, including dU_i/dQ through dUC, and Q d/dQ = d/dt):
 r_{1,i}[0] = L_i r_{0,i}[0], and r_{1,i}[m] = L_i r_{0,i}[m] + (d/dt - Delta_i'/(2 Delta_i)) r_{0,i}[m-1] for m >= 1, i.e.
   (H)   r_{1,i}(z) = (L_i + z D_i) r_{0,i}(z),   D_i := d/dt - Delta_i'/(2 Delta_i).
 (H) is K-linear in r_{0,i}, so it commutes with the twist by h_w(z)^-1: it holds in frame B and in the point frame.

(0.6) Input from S1 (Steps 1-2 of the outline); S2.10 uses in addition the general formal Laplace calculus of S1 (header, input (ii)).
 (S1.a) r_{0,1} is in O[[z]] (every z-coefficient is regular at t = 0); hence r_{0,1} is in C[[z, tau]], a subring of R_0. [Theorem S1(b).]
 (S1.b) In R (in fact in K[[tau]]), in frame B [Theorem S1(c), (d)],
          r_{0,1}(z; t) = (c/(c + z)) alg(tau) F(z, tau) Ephi(z, tau),     and r^pt_{0,1} = h_w(z) times the right-hand side, where
   F := 2F1(1/z + 1, 1/z + 1; c/z + 2; tau) = Sum_n F_n(z) tau^n,  F_n = prod_{j=1}^n (1 + j z)^2 / (n! z^n prod_{j=2}^{n+1} (c + j z)),  an element of
        C(z)[[tau]], a subring of K[[tau]] (F_0 = 1; F_n has a pole of order n at z = 0);
   dphi := phi_tau(s*) - phi_0(1/c), in tau O, with phi_tau(s) := log s + x log(1 - s) - log(1 - tau s);
   Ephi := e^{-dphi/z} := Sum_k (-dphi)^k/(k! z^k), in C[1/z][[tau]] (tau-adically convergent; [tau^n] Ephi is a polynomial of degree <= n in 1/z);
   alg := (1 - tau s*) (phi''_tau(s*)/phi''_0(1/c))^(1/2), in O, branch alg(0) = 1 (the ratio lies in 1 + tau O).
 F and Ephi are units of K[[tau]] (constant term 1, Ephi^-1 = e^{dphi/z} = refl(Ephi)), and z Ephi' = -dphi' Ephi (termwise differentiation of a
 tau-adically convergent series). Gauss's equation Gauss[F] = 0 (notation of Lemma S2.6) holds in K[[tau]]: its tau^n coefficient is (n + 1)(n + gamma) F_{n+1}
 - (n + alpha)(n + beta) F_n, and F_{n+1}/F_n = (1 + (n+1) z)^2/((n+1) z (c + (n+2) z)) = (n + alpha)(n + beta)/((n+1)(n + gamma)) for
 alpha = beta = 1/z + 1, gamma = c/z + 2 (gamma + n = (c + (n+2) z)/z is invertible in K). Setting tau = 0 in (S1.b) gives F2: r_{0,1}(z; 0) = c/(c + z).


## S2.1 Frobenius data at the torus point

Lemma S2.1. At alpha0 = 0, alpha1 = eps1 = x, alpha2 = eps2 = 1, in the flat basis (1, H):
 (a) eta = {{c/x^3, 1/x^2}, {1/x^2, 0}},   eta^-1 = {{0, x^2}, {x^2, -c x}}   (independent of t);
 (b) 1 * a = a for all a, and  H * H = -(x/tau) 1 + (c/tau) H;
 (c) Leq = L^2 - (c/tau) L + x/tau, Ssum = c/tau, Pprod = x/tau (rc 23-25);
 (d) the product is commutative, associative and Frobenius: eta(a * b, d) = eta(a, b * d).
Proof. (a) rc 16: Replp2 = -Expand[p (p + a0)] + p^2 = -a0 p, so p^2 is reduced by the relation p (p + a0) = 0, whose roots are the fixed points
 p = 0 and p = -a0 of P^1. The reduced cube exp3 (rc 17) agrees with (t0 + t1 p)^3 at both roots, so rc 18 is the localisation formula
   PrepCl = (1/6) [ t0^3/(a0 a1 a2 e1 e2) - (t0 - a0 t1)^3/(a0 (a0 + a1)(a0 + a2) e1 e2) ]
 (a cubic polynomial, so the Series truncation in rc 18 is exact). Its third derivatives are
   PrepCl_000 = (a0 + a1 + a2)/(a1 a2 (a0 + a1)(a0 + a2) e1 e2)   [since 1/(a1 a2) - 1/((a0+a1)(a0+a2)) = a0 (a0 + a1 + a2)/(a1 a2 (a0+a1)(a0+a2))],
   PrepCl_001 = 1/((a0 + a1)(a0 + a2) e1 e2),   PrepCl_011 = -a0/((a0 + a1)(a0 + a2) e1 e2),   PrepCl_111 = a0^2/((a0 + a1)(a0 + a2) e1 e2),
 i.e. the entries of etaSym (rc 82) and PrepCl_111. None has a0 in a denominator, so a0 = 0 may be set. eta_ij = Prep_0ij (rc 19). At the torus
 point a1 a2 = x, (a0 + a1)(a0 + a2) = x, e1 e2 = x and a0 + a1 + a2 = c, which gives eta; det eta = -1/x^4, and eta^-1 follows.
 (b) Prep = PrepCl + Li_3(e^{t1})/(e1 e2) (rc 21), and (d/dt1)^3 Li_3(e^{t1}) = e^{t1}/(1 - e^{t1}). With t1 = log Q (rc 22) and
 Q/(1 - Q) = 1/(1/Q - 1) = -1/tau:  Prep_111 = PrepCl_111 + Q/((1 - Q) x) = -1/(x tau); all other third derivatives are classical.
 cch[i, j, k] = Sum_m eta^{im} Prep_mjk (rc 22) are the structure constants, e_j * e_k = Sum_i cch[i, j, k] e_i (e_0 = 1, e_1 = H).
 Prep_0jk = eta_jk gives cch[i, 0, k] = delta_ik (1 is the unit), and
   cch[0, 1, 1] = eta^00 Prep_011 + eta^01 Prep_111 = x^2 (-1/(x tau)) = -x/tau,
   cch[1, 1, 1] = eta^10 Prep_011 + eta^11 Prep_111 = (-c x)(-1/(x tau)) = c/tau.
 (c) rc 23-25 with (b). (d) The algebra is generated by H, and eta(e_j * e_k, e_l) = Sum_i cch[i, j, k] eta_il = Prep_ljk is totally symmetric. QED


## S2.2 Critical points, eigenvalues, canonical coordinate

W(U) := x log(1 - Q U) - log U + log(1 - U) (spotU, rc 80, at a0 = 0, a1 = x, a2 = 1). Its derivatives are single valued:
   U W_U = -x Q U/(1 - Q U) - 1/(1 - U),   W2(U) := (U d/dU)^2 W = -x Q U/(1 - Q U)^2 - U/(1 - U)^2,   (partial_t W)(U) = -x Q U/(1 - Q U),
because -1 - U/(1 - U) = -1/(1 - U), U d/dU [Q U/(1 - Q U)] = Q U/(1 - Q U)^2 and U d/dU [U/(1 - U)] = U/(1 - U)^2.

Lemma S2.2. (a) U W_U = 0 iff x Q U^2 + Q (1 - x) U - 1 = 0. This is LeqU (rc 115) at the torus point; U_1 + U_2 = (x - 1)/x = SU and
 U_1 U_2 = -1/(Q x) = PU (rc 113).
 (b) At a critical point: (partial_t W)(U) = 1/(1 - U) = (x U + 1)/tau = L(U), where L(U) = -Q (x U + 1)/(1 - Q) is LofU (rc 81).
 (c) Put U = 1 - tau s (tau != 0). U is critical iff x tau s^2 - c s + 1 = 0. The critical points are
   U_1 = 1 - tau s* = (x - 1 + delta)/(2x) -> 1,     U_2 = 1 - tau s~ = 1 - 1/(x s*) -> -1/x     (t -> 0).
 U_1 is the root (S0 + sq)/2 of hvc_leaf_proof_check.wls (line 22; sq = delta/x there), and s* = (1 - U_1)/tau is the s* of (V2).
 (d) L_1 := L(U_1) = 1/(tau s*) = c/tau - x s*,   L_2 := L(U_2) = x s* = 1/(tau s~),
   L_1 + L_2 = c/tau,   L_1 L_2 = x/tau,   L_1 - L_2 = delta/tau.
 Hence L_1, L_2 are the roots of Leq, the eigenvalues of H *. L_1 = c/t + l(t) with l in O, l(0) = c/2 - x/c; L_2 is in O, L_2(0) = x/c.
 (e) u_1 := W(U_1) (any branch of the logarithms) satisfies du_1/dt = L_1. Explicitly,
   W(1 - tau s) = x t + c log tau + phi_tau(s) + x log(-1),   so   u_1 = x t + c log tau + dphi + const,   and
   dphi' = (1 - tau) s*/(1 - tau s*) = 1 - x s*,    L_1 = x + c (1 - tau)/tau + dphi',    L_1 - dphi' = c/tau - 1.
Proof. (a) U W_U = -[x Q U (1 - U) + (1 - Q U)]/((1 - Q U)(1 - U)), and x Q U (1 - U) + 1 - Q U = -(x Q U^2 + Q (1 - x) U - 1).
 (b) U W_U = 0 gives -x Q U/(1 - Q U) = 1/(1 - U). By (a), x U^2 = 1/Q - (1 - x) U, so (x U + 1)(1 - U) = 1 + (x - 1) U - x U^2 = 1 - 1/Q = tau.
 Finally -Q/(1 - Q) = 1/tau.
 (c) 1 - U = tau s and 1 - Q U = Q (1/Q - 1 + tau s) = -Q tau (1 - s). For tau != 0 the relation (x U + 1)(1 - U) = tau of (b) is equivalent to (a),
 and it becomes (c - x tau s) tau s = tau, i.e. x tau s^2 - c s + 1 = 0. By (0.2), 1 - U_1 = tau s* = (c - delta)/(2x) and tau s~ = 1/(x s*).
 (d) By (b), L(U) = 1/(1 - U) = 1/(tau s) at s = s*, s~. So L_2 = 1/(tau s~) = x s*, and L_1 = 1/(tau s*) = c/tau - x s* is the critical equation
 divided by tau s*. L_1 + L_2 = (s* + s~)/(tau s* s~) = (c/(x tau)) x = c/tau, L_1 L_2 = 1/(tau^2 s* s~) = x/tau, and L_1 - L_2 = (c - 2 x tau s*)/tau
 = delta/tau. By Lemma S2.1(b) the minimal polynomial of H * is Leq, so L_1, L_2 are its eigenvalues. 1/tau = 1/t + 1/2 + O(t) gives l.
 (e) u_1' = (partial_t W)(U_1) + W_U(U_1) dU_1/dt = L_1 by (b) (U_1 is analytic in t). For the explicit form: log(1 - Q U) = log(-1) + t + log tau
 + log(1 - s), log U = log(1 - tau s), log(1 - U) = log tau + log s, with branches fixed on a simply connected domain of small t != 0 (other branches
 change u_1 by a constant). Directly,
   phi'_tau(s) = 1/s - x/(1 - s) + tau/(1 - tau s) = (x tau s^2 - c s + 1)/(s (1 - s)(1 - tau s)),
 so s* is a critical point of phi_tau and d/dtau [phi_tau(s*(tau))] = (partial_tau phi_tau)(s*) = s*/(1 - tau s*); hence dphi' = (1 - tau) s*/(1 - tau s*).
 The identity (1 - tau) s* = (1 - x s*)(1 - tau s*) is once more x tau s*^2 - c s* + 1 = 0. Finally x + c (1 - tau)/tau + 1 - x s* = c/tau - x s* + (x + 1 - c)
 = L_1. QED


## S2.3 Idempotents and the metric: (V1), (V2), F1

Prop S2.3. (a) The columns of Psi^-1, Psi = {{1, L1}, {1, L2}} (rc 29), are (components with respect to (1, H))
   e_1 = (H - L_2)/(L_1 - L_2),    e_2 = (H - L_1)/(L_2 - L_1),
 and they are the idempotents: e_i * e_j = delta_ij e_i, H * e_i = L_i e_i, e_1 + e_2 = 1.
 (b) etacan = Psi^-T eta Psi^-1 = diag(eta(e_1, e_1), eta(e_2, e_2)), and the code's Delta (DeltaL rc 31, DeltaCN rc 128-129) is
 Delta_i = 1/eta(e_i, e_i): Delta_1 at U_1 (DeltaCN[1]) and Delta_2 = its conjugate at U_2 (DeltaCN[2]).
 (c) For {i, j} = {1, 2}:
   Delta_i = x^3 (L_i - L_j)/(x - c L_j) = x L_i (2x - c L_i) = x^2 W2(U_i).
 In particular Delta_1 = x^2 W2(U_1) is an identity in the quadratic algebra C(Q)[U]/(LeqU) (V1), and Delta_2 = x^2 W2(U_2).
 (d) Delta_1 = -x delta/(tau^2 s*^2),  Delta_2 = x^3 delta s*^2,  x - c L_2 = -x^2 tau s*^2,  x - c L_1 = -1/(tau s*^2),
   Delta_1 Delta_2 = -x^4 delta^2/tau^2,   1/Delta_1 + 1/Delta_2 = c/x^3 = eta(1, 1).
 (e) F1: Delta_2 is in O with Delta_2 = x^3/c - x^5 tau^2/c^5 + O(tau^3), so Delta_2(0) = x^3/c != 0 (and Delta_2'(0) = 0, as in S4 10.3);
 Delta_1 = -x c^3 t^-2 (1 + O(t)) (Prop S2.4(c)).
Proof. (a) Psi^-1 = (L_2 - L_1)^-1 {{L_2, -L_1}, {-1, 1}}: its first column (L_2, -1)/(L_2 - L_1) is (L_2 1 - H)/(L_2 - L_1) = e_1, the second is e_2.
 By Lemma S2.1(b) and S2.2(d), H * H = (L_1 + L_2) H - L_1 L_2 1, i.e. (H - L_1) * (H - L_2) = 0; hence H * (H - L_2) = L_1 (H - L_2), i.e. H * e_1 = L_1 e_1,
 and e_1 * e_2 = 0. As e_1 + e_2 = 1, e_1 * e_1 = e_1 * (e_1 + e_2) = e_1; likewise for e_2.
 (b) The (i, j) entry of Psi^-T eta Psi^-1 is eta(col_i, col_j) = eta(e_i, e_j) = eta(e_i * e_j, 1) (Lemma S2.1(d)), which vanishes for i != j.
 (c) First expression, Frobenius route: eta(e_i, e_i) = eta(e_i * e_i, 1) = eta(e_i, 1) = (eta(H, 1) - L_j eta(1, 1))/(L_i - L_j)
 = (1/x^2 - c L_j/x^3)/(L_i - L_j) = (x - c L_j)/(x^3 (L_i - L_j)).
 The code evaluates the quadratic form literally (rc 30-31: L1, L2 independent, then L2 -> Ssum - L1, then reduction modulo Leq):
 e_1^T eta e_1 = (c L_2^2/x^3 - 2 L_2/x^2)/(L_1 - L_2)^2 = L_2 (c L_2 - 2x)/(x^3 (L_1 - L_2)^2). The two agree on the critical locus, since
 L_2 (c L_2 - 2x) - (x - c L_2)(L_1 - L_2) = c L_1 L_2 - x (L_1 + L_2) = c x/tau - x c/tau = 0 (this is -1/x times <(1, L_1), (1, L_2)>, Prop S2.8(d)).
 Second expression: L_i (2x - c L_i)(x - c L_j) = L_i [2 x^2 - c x L_j - c x (L_i + L_j) + c^2 L_i L_j] = L_i (2 x^2 - c x L_j), because
 c x (L_i + L_j) = c^2 x/tau = c^2 L_i L_j; and L_i (2 x^2 - c x L_j) = 2 x^2 L_i - x^2 (L_i + L_j) = x^2 (L_i - L_j).
 Third expression: at a critical point -x Q U/(1 - Q U) = L (S2.2(b)), so Q U/(1 - Q U) = -L/x, 1/(1 - Q U) = 1 + Q U/(1 - Q U) = (x - L)/x and
 -x Q U/(1 - Q U)^2 = L (x - L)/x; and U = 1 - 1/L gives U/(1 - U)^2 = L (L - 1). Therefore
   W2(U) = L (x - L)/x - L (L - 1) = L (2x - c L)/x    at every critical point.
 The quadratic algebra C(Q)[U]/(LeqU) has the two characters "evaluate at U_1", "evaluate at U_2" (LeqU has distinct roots: delta != 0 near tau = 0),
 and an element vanishes iff both its values do. DeltaCN[1] = toC[DeltaL(LofU)] (rc 128-129) and toC[x^2 W2] have equal values at U_1 and at U_2,
 hence are equal: this is check (1) of hvc_leaf_proof_check.wls.
 (d) With L_1 = 1/(tau s*) and 2 x tau s* - c = -delta: Delta_1 = x L_1 (2x - c L_1) = x (2 x tau s* - c)/(tau^2 s*^2) = -x delta/(tau^2 s*^2).
 With L_2 = x s* and 2 - c s* = 2 - 2c/(c + delta) = delta s*: Delta_2 = x^2 s* (2x - c x s*) = x^3 s* (2 - c s*) = x^3 delta s*^2.
 x - c L_2 = x (1 - c s*) = -x^2 tau s*^2 by the critical equation; x - c L_1 = (x tau s* - c)/(tau s*) with x tau s* - c = (c - delta)/2 - c = -1/s*.
 The product is immediate, and 1/Delta_1 + 1/Delta_2 = eta(e_1 + e_2, 1) = eta(1, 1) = c/x^3 by (c).
 (e) delta and s* are in O with delta(0) = c, s*(0) = 1/c. With u := 4 x tau/c^2: delta/c = 1 - u/2 - u^2/8 + O(u^3) and
 (c + delta)/c = 2 - u/2 - u^2/8 + O(u^3), so delta s*^2 = 4 (delta/c)/(c ((c + delta)/c)^2) = (1/c)(1 - u^2/16 + O(u^3)) = 1/c - x^2 tau^2/c^5 + O(tau^3). QED


## S2.4 Step 3: the Hessian in the Euler coordinate; rho = tau/(t alg)

Prop S2.4. (a) For every small t != 0 and both critical points, with (s_1, s_2) := (s*, s~):
   W2(U_i) = (1 - tau s_i)^2 phi''_tau(s_i)/tau^2,     phi''_tau(s) = -1/s^2 - x/(1 - s)^2 + tau^2/(1 - tau s)^2.
 (b) phi''_0(1/c) = -c^3/x and (1 - tau s*)^2 phi''_tau(s*) = -delta/(x s*^2). Hence alg^2 = delta/(c^3 s*^2) and alg = (delta/c)^(1/2)/(c s*),
 a unit of O with alg(0) = 1.
 (c) Delta_1 = x^2 W2(U_1) = x^2 (1 - tau s*)^2 phi''_tau(s*)/tau^2 = -x c^3 alg^2/tau^2. Consequently
   rho := (Delta_1 t^2/(-x c^3))^(-1/2) = tau/(t alg),   a unit of O,  rho(0) = 1,   and   rho^2 Delta_1 = -x c^3/t^2   exactly,
 where the square root is the branch with rho(0) = 1, as in leaf_structure.md and hvc_leaf_proof_check.wls (line 36, check (2c)). (rho occurs in neither
 rc_lib.m nor hvc_lib.m: the graph factors never use it; it is defined by this formula, outline s.0.)
 Also Delta_1 t^2/(-x c^3) = (t/tau)^2 alg^2 is a unit of O with value 1 at t = 0: F1 for Delta_1.
 (d) With (S1.b): G := rho r_{0,1} = (c/(c + z)) (tau/t) F(z, tau) Ephi(z, tau), which is the KEY FORMULA (V3); G is in C[[z, t]] by (S1.a), and
 G(z; 0) = c/(c + z).
 (e) rho = 1 - (1-x)^2 t/(2 c^2) + O(t^2) and Delta_1 = -(x c^3/t^2)(1 + (1-x)^2 t/c^2 + O(t^2)), as in leaf_structure.md.
Proof. (a) Fix t != 0. By S2.2(e), W(1 - tau s) = phi_tau(s) + (a function of t), and s -> U = 1 - tau s is affine with dU/ds = -tau, so
 U d/dU = a(s) d/ds with a(s) = -(1 - tau s)/tau. Then (U d/dU)^2 W = a (a W_s)_s = a^2 W_ss + a a_s W_s. At s_i, W_s = phi'_tau(s_i) = 0 and
 W_ss = phi''_tau(s_i), while a^2 = (1 - tau s)^2/tau^2. (This is the usual tensor rule for a Hessian at a critical point.)
 (b) phi''_0(1/c) = -c^2 - x/(1 - 1/c)^2 = -c^2 - c^2/x = -c^3/x. By (a) and Prop S2.3(c), (d), (1 - tau s*)^2 phi''_tau(s*) = tau^2 W2(U_1)
 = tau^2 Delta_1/x^2 = -delta/(x s*^2). Hence, by the definition (0.6),
   alg^2 = (1 - tau s*)^2 phi''_tau(s*)/phi''_0(1/c) = delta/(c^3 s*^2),
 a unit of O with value c/(c^3/c^2) = 1 at tau = 0. alg is (1 - tau s*), a unit with value 1, times the square root with value 1 of the ratio
 phi''_tau(s*)/phi''_0(1/c) (an element of 1 + tau O); so alg is the square root of alg^2 with value 1 at 0. (delta/c)^(1/2)/(c s*) (principal root of
 delta/c in 1 + tau O) is also a square root of alg^2 with value 1 at 0. Two square roots of the same unit of O agree up to a sign, and the sign is fixed by
 the value at 0; hence alg = (delta/c)^(1/2)/(c s*).
 (c) The first equality is Prop S2.3(c), the second is (a), the third is (b) (tau^2 Delta_1/x^2 = -delta/(x s*^2) = -c^3 alg^2/x). Then
 rho^-2 = Delta_1 t^2/(-x c^3) = (t alg/tau)^2, and its square root with value 1 at t = 0 is t alg/tau.
 (d) rho r_{0,1} = (tau/(t alg)) (c/(c + z)) alg F Ephi. At tau = 0: F = 1, Ephi = 1, tau/t = 1.
 (e) tau = t - t^2/2 + O(t^3), delta/c = 1 - 2 x tau/c^2 + O(tau^2), c s* = 2/(1 + delta/c) = 1 + x tau/c^2 + O(tau^2), so alg = (delta/c)^(1/2)/(c s*)
 = 1 - 2 x tau/c^2 + O(tau^2) and rho = (tau/t)/alg = 1 + (2x/c^2 - 1/2) t + O(t^2) = 1 - (c^2 - 4x) t/(2 c^2) + O(t^2), with c^2 - 4x = (1 - x)^2.
 Delta_1 t^2/(-x c^3) = rho^-2. QED
Checks: S2_checks.pl part E compares the t-series of tau/(t alg) (to t^4) and of (t/tau)^2 alg^2 (to t^3) with the series of rho and Delta_1 in
leaf_structure.md (obtained there from the stored Delta_1), and dphi = Int (1 - x s*) dt, Theta = -x^3 tau s*^2 (to t^7) with the series printed in
hvc_leaf_proof_check_out.txt, at x = 2, 3/5, -3: all identical. Mathematica: (V2), (V3) and checks (2c), (3), (3b) of hvc_leaf_proof_check.wls.


## S2.5 Formal exponentials, square roots, t-independent constants

The outline writes s_a := e^{u_1/z} Delta_1^(-1/2) r_{a,1} and s_0 = kappa(z) Q^{x/z} tau^{c/z+1} F. Their precise meaning is as follows.

(5.1) Exponential symbols. For lambda in C((tau)) let E_lambda := R eps_lambda be a free R-module of rank one with the operator
   nabla (f eps_lambda) := (z f' + lambda f) eps_lambda.
 We write f e^{u/z} := f eps_{u'} for any (possibly multivalued) function u with u' = lambda. So e^{u/z} depends only on u', and u, u + const
 define the same symbol. The product (f eps_lambda)(h eps_mu) := f h eps_{lambda+mu} satisfies the Leibniz rule for nabla, since
 z (f h)' + (lambda + mu) f h = (z f' + lambda f) h + f (z h' + mu h); E_0 = R with nabla = z d/dt. If gamma is a unit of R with
 z gamma' = (lambda - mu) gamma, then f eps_lambda -> gamma f eps_mu is an isomorphism (E_lambda, nabla) -> (E_mu, nabla); this is the meaning of an
 equation e^{u/z} = gamma e^{v/z}. The substitution z -> -z acts by refl(f eps_lambda) := refl(f) eps_{-lambda}, and refl nabla = -nabla refl; the pairing
 E_lambda x E_{-lambda} -> R, (f eps_lambda, h eps_{-lambda}) -> f h, satisfies z (f h)' = (z f' + lambda f) h + f (z h' - lambda h).
 Q^{x/z} = e^{x t/z} and e^{+-dphi/z} are honest units of R (x t and dphi lie in tau O, so the exponential series converge tau-adically). The only
 genuinely formal factor in s_0 is tau^{c/z} = e^{c log(tau)/z}; together with the honest factor tau it carries the local exponent c/z + 1 at the
 regular singular point tau = 0 of the Picard-Fuchs operator (Remark S2.7.4).
 Gauss's operator on E_lambda. Put
   d/dtau (f eps_lambda) := (z (1 - tau))^-1 nabla (f eps_lambda) = (f_tau + lambda f/(z (1 - tau))) eps_lambda,    f_tau := df/dtau = (1 - tau)^-1 f'.
 It satisfies the Leibniz rule d/dtau (a Y) = a_tau Y + a d/dtau Y (a in R, Y in E_lambda), and on E_0 = R it is d/dtau. For Y in E_lambda put
   Gauss[Y] := tau (1 - tau) (d/dtau)^2 Y + (gamma - (alpha + beta + 1) tau) d/dtau Y - alpha beta Y     (alpha, beta, gamma as in Lemma S2.6).
 Every identity between differential operators in D = z d/dt (equivalently in d/dtau = (z (1 - tau))^-1 D) with coefficients in
 B0 := C[z, 1/z][tau, 1/tau, 1/(1 - tau)] that is proved from the Leibniz rule and D(tau) = z (1 - tau) alone holds verbatim on every E_lambda, with D := nabla;
 this is how Lemma S2.6 is applied on E_lambda.
 z-adic version (used only in S2.10). For lambda in M let E^M_lambda := M((z)) eps_lambda, with the same nabla, product, refl, pairing and d/dtau. Everything
 in (5.1) holds verbatim (only the Leibniz rule and the d/dt-linearity of the formulas were used), and B0 is a subring of M((z)) as well as of R.
(5.2) Square roots. Fix kappa_0 in C with kappa_0^2 = -x c^3, and put
   nu_1 := tau/(kappa_0 alg) = t rho/kappa_0,   an element of tau O whose quotient by tau is a unit.
 Then nu_1^2 = 1/Delta_1 (Prop S2.4(c)) and nu_1'/nu_1 = -Delta_1'/(2 Delta_1), so for every f in R, D_1 f = nu_1^-1 (nu_1 f)', i.e.
   L_1 + z D_1 = nu_1^-1 (z d/dt + L_1) nu_1    as operators on R.                                           (5.2.1)
 The sign of kappa_0 never matters: D_1 depends on Delta_1 only through its logarithmic derivative, and everything below involves nu_1 only through
 (5.2.1) or through nu_1^2 = 1/Delta_1.
(5.3) Dictionary. Put mu := x + c (1 - tau)/tau = c/tau - 1, so that "Q^{x/z} tau^{c/z}" = eps_mu. By S2.2(e), L_1 = mu + dphi', and gamma := e^{dphi/z}
 satisfies z gamma' = dphi' gamma, so E_{L_1} ~= E_mu. The outline's s_0 is S_0 := nu_1 r_{0,1} eps_{L_1} ("e^{u_1/z} Delta_1^(-1/2) r_{0,1}"), and by (S1.b)
 and the definition of nu_1 its image in E_mu is
   e^{dphi/z} nu_1 r_{0,1} eps_mu = e^{dphi/z} (tau/(kappa_0 alg)) (c/(c + z)) alg F Ephi eps_mu = kappa_1(z) tau F eps_mu,   kappa_1(z) := c/(kappa_0 (c + z)).
 This is the precise form of "s_0 = kappa(z) Q^{x/z} tau^{c/z+1} F" (outline, Step 5; this kappa is not S1's phase constant kappa_ph):
 kappa(z) = kappa_1(z) in frame B and h_w(z) kappa_1(z) in the point frame, a unit of K
 (t-independent). The additive constants of u_1 (the phase x log(-1), phi_0(1/c), the branches of the logarithms), which would give factors e^{const/z}
 that are not in K, do not occur at all.
(5.4) Why kappa does not matter. Rows 1 and 2 of the QDE (S2.7) are K-linear in r_{0,1}, and (H) commutes with multiplication by a unit kappa of K.
 Replacing r_{0,1} by kappa r_{0,1} replaces r_{1,1} by kappa r_{1,1}, preserves the second QDE row (Prop S2.7), and replaces K_1(z) (Prop S2.9) by
 kappa(z) kappa(-z) K_1(z). For the frame change kappa = h_w one has h_w(z) h_w(-z) = 1 (the exponent of h_w is odd in z), so K_1 and the symplectic
 relations are the same in frame B and in the point frame. The same holds at idempotent 2 (K_2, S2.10).


## S2.6 The identity (V7)

Lemma S2.6. Let D := z d/dt = z (1 - tau) d/dtau, omega := (c + z)/tau - (1 + z), and
   Gauss[f] := tau (1 - tau) f_tautau + (gamma - (alpha + beta + 1) tau) f_tau - alpha beta f,    alpha = beta = 1/z + 1,  gamma = c/z + 2.
 Then, as an identity of differential operators with coefficients in B0 = C[z, 1/z][tau, 1/tau, 1/(1 - tau)] (multiplication rule D a = a D + D(a)),
 hence on R, and, with D replaced by nabla and d/dtau as in (5.1), on every E_lambda and every E^M_lambda:
   (D + omega)^2 - (c/tau) (D + omega) + x/tau  =  ((1 - tau) z^2/tau) Gauss.                                  (S2.6.1)
 For P := tau eps_mu ("(1 - tau)^(-x/z) tau^(c/z+1)", 1 - tau = 1/Q) one has nabla(P f) = P (D + omega) f, because omega = mu + D(tau)/tau
 = x + (c + z)(1 - tau)/tau. So (S2.6.1) says [nabla^2 - (c/tau) nabla + x/tau](P f) = -((tau - 1) z^2/tau) P Gauss[f]: this is (V7).
 Module form: for every lambda and every Y in E_lambda (resp. E^M_lambda, with mu, lambda in M), in E_{mu+lambda} (resp. E^M_{mu+lambda}),
   nabla(P Y) = P (nabla + omega) Y,    [nabla^2 - (c/tau) nabla + x/tau](P Y) = ((1 - tau) z^2/tau) P Gauss[Y].        (S2.6.2)
Proof. (D + omega)^2 f = D^2 f + 2 omega D f + (D omega) f + omega^2 f, so the left side of (S2.6.1) applied to f is
   D^2 f + (2 omega - c/tau) D f + (D omega + omega^2 - (c/tau) omega + x/tau) f,
 with D f = z (1 - tau) f_tau and D^2 f = z^2 (1 - tau)^2 f_tautau - z^2 (1 - tau) f_tau. The coefficients are:
 * of f: D omega = -z (1 - tau)(c + z)/tau^2 and omega (omega - c/tau) = [(c + z)/tau - (1 + z)][z/tau - (1 + z)] = z (c + z)/tau^2 - (1 + z)(c + 2z)/tau + (1 + z)^2.
   The tau^-2 terms add up to z (c + z)(1 - (1 - tau))/tau^2 = z (c + z)/tau, and z (c + z) - (1 + z)(c + 2z) + x = -z^2 - 2z - (c - x) = -(1 + z)^2.
   So the coefficient is -(1 + z)^2/tau + (1 + z)^2 = -(1 - tau)(1 + z)^2/tau.
 * of f_tau: -z^2 (1 - tau) + z (1 - tau)[(c + 2z)/tau - 2 (1 + z)] = z (1 - tau)[(c + 2z)/tau - (2 + 3z)].
 * of f_tautau: z^2 (1 - tau)^2.
 On the right, ((1 - tau) z^2/tau) Gauss[f] has f_tautau-coefficient z^2 (1 - tau)^2, f_tau-coefficient (1 - tau) z^2 [(c/z + 2)/tau - (2/z + 3)]
 = z (1 - tau)[(c + 2z)/tau - (2 + 3z)], and f-coefficient -((1 - tau) z^2/tau)(1 + 1/z)^2 = -(1 - tau)(1 + z)^2/tau. All three agree. Only the Leibniz
 rule and D(tau) = z (1 - tau) were used (D f and D^2 f in terms of f_tau, f_tautau are the same computation), so this is an identity in B0<D>; on
 E_lambda and E^M_lambda the formulas nabla Y = z (1 - tau) d/dtau Y and nabla^2 Y = z^2 (1 - tau)^2 (d/dtau)^2 Y - z^2 (1 - tau) d/dtau Y hold by the
 definition of d/dtau and its Leibniz rule (5.1), so the same computation applies. For P: in E_mu, nabla(tau eps_mu) = (D(tau) + mu tau) eps_mu
 = P (z (1 - tau)/tau + mu) and z (1 - tau)/tau + mu = x + (c + z)(1 - tau)/tau = (c + z)/tau - (1 + z) = omega (x - c = -1); by Leibniz
 nabla(P f) = nabla(P) f + P D f = P (D + omega) f. (S2.6.2): by the Leibniz rule for products of module elements (5.1), nabla(P Y) = nabla(P) Y + P nabla(Y)
 = P (nabla + omega) Y; iterating, Pi(nabla)(P Y) = P Pi(nabla + omega)(Y) for every polynomial Pi in nabla with coefficients in B0 (multiplication by
 elements of B0 commutes with multiplication by P), and (S2.6.1) on E_lambda (resp. E^M_lambda) gives the second identity. QED
 (hvc_pf_gauss.wls checks the same identity with a symbolic F; S2_checks.pl part B checks (S2.6.1) exactly on tau^n, -3 <= n <= 7, and on a mixed
 Laurent polynomial, for x in {2, 3/5, -3, -7/4} and z in {2/7, -5/3, 4/9, 11}. (S2.6.2) on E^M_{lambda_2} is what S2_checks_reg.pl part P tests through
 Prop S2.10.4.)


## S2.7 Step 5: the QDE at the singular idempotent

Prop S2.7. Assume (S1.a), (S1.b), and let A := {{0, 1}, {-x/tau, c/tau}}. The following hold in R_0 (r_{0,1} in C[[z, tau]], r_{1,1} in tau^-1 C[[z, tau]]),
 hence, by TRANSFER, to all orders in z as identities of meromorphic functions of t near 0:
   Row 1:   r_{1,1} = (L_1 + z D_1) r_{0,1}                                     (the definition (H));
   Row 2:   (L_1 + z D_1) r_{1,1} = -(x/tau) r_{0,1} + (c/tau) r_{1,1},   equivalently   z D_1 r_{1,1} = L_2 (r_{1,1} - L_1 r_{0,1}),
            i.e. (D_1 r_{1,1})[m-1] = L_2 (r_{1,1}[m] - L_1 r_{0,1}[m]) for every m >= 1: exactly the form tested in (V5) for i = 1 (m <= 9 there).
 In matrix form: z D_1 v_1(z) = (A - L_1) v_1(z); equivalently z w' = (A - L_1) w for w(z) := nu_1 v_1(z), which lies in C[[z, tau]]^2; equivalently
   nabla S = A S   for   S := nu_1 v_1 eps_{L_1} = "e^{u_1/z} Delta_1^(-1/2) v_1"  in (E_{L_1})^2,
 which is the outline's "z ds/dt = A s"; in components nabla S_0 = S_1 (Row 1) and nabla S_1 = -(x/tau) S_0 + (c/tau) S_1 (Row 2).
Proof. (i) Conjugation by tau Ephi. For f in R, using z Ephi' = -dphi' Ephi, tau' = 1 - tau and S2.2(e),
   (z d/dt + L_1)(tau Ephi f) = tau Ephi [z f' + (L_1 - dphi' + z (1 - tau)/tau) f] = tau Ephi (D + omega) f,
 because L_1 - dphi' + z (1 - tau)/tau = c/tau - 1 + z/tau - z = omega.
 (ii) By (S1.b), nu_1 r_{0,1} = (tau/(kappa_0 alg)) (c/(c + z)) alg F Ephi = kappa_1 tau Ephi F. By (5.2.1) and (i),
   nu_1 r_{1,1} = (z d/dt + L_1)(nu_1 r_{0,1}) = kappa_1 tau Ephi (D + omega) F,
   nu_1 (L_1 + z D_1) r_{1,1} = (z d/dt + L_1)(nu_1 r_{1,1}) = kappa_1 tau Ephi (D + omega)^2 F.
 (iii) Hence, by (S2.6.1) and Gauss[F] = 0 (0.6),
   nu_1 [(L_1 + z D_1) r_{1,1} + (x/tau) r_{0,1} - (c/tau) r_{1,1}] = kappa_1 tau Ephi [(D + omega)^2 F - (c/tau)(D + omega) F + (x/tau) F]
                                                                  = kappa_1 Ephi (1 - tau) z^2 Gauss[F] = 0,
 and nu_1 is a unit of R. The second form of Row 2 follows from c/tau - L_1 = L_2 and x/tau = L_1 L_2. The first component of (A - L_1) v_1 is
 r_{1,1} - L_1 r_{0,1} = z D_1 r_{0,1} (Row 1), the second is -(x/tau) r_{0,1} + (c/tau - L_1) r_{1,1} = z D_1 r_{1,1} (Row 2), and z (nu_1 f)' = nu_1 z D_1 f.
 In E_{L_1}, nabla(nu_1 f eps_{L_1}) = nu_1 (L_1 + z D_1) f eps_{L_1}, which gives the last form.
 Membership: L_1 and Delta_1'/Delta_1 lie in tau^-1 O (Delta_1 is tau^-2 times a unit of O), so L_1 + z D_1 maps C[[z, tau]] into tau^-1 C[[z, tau]] and
 tau^-1 C[[z, tau]] into tau^-2 C[[z, tau]]; with (S1.a) every term lies in R_0, and TRANSFER (0.3) applies. w = nu_1 v_1 is in C[[z, tau]]^2 since nu_1 is
 in tau O. QED

Remark S2.7.1 (what is used). Given Row 1 and (S1.b), Row 2 is EQUIVALENT to Gauss[F] = 0, since kappa_1 Ephi (1 - tau) z^2 is a unit of R. The closed form
 enters only through nu_1 r_{0,1} = kappa_1 tau Ephi F with dphi' = 1 - x s*. No property of the Wick series beyond (S1.a), (S1.b) is used, and nothing
 about the regular idempotent.
Remark S2.7.2 (closed form of the singular H-row). From (ii), since F is a unit of K[[tau]],
   r_{1,1} = r_{0,1} [(D + omega) F]/F = r_{0,1} [ (c + z)/tau - (1 + z) + z (1 - tau) F_tau/F ];
 equivalently rho r_{1,1} = (L_1 + z/t) G + z G' (Prop S2.9(b)). The two agree: z G'/G = z (1 - tau)/tau - z/t + z (1 - tau) F_tau/F - dphi' by Prop S2.4(d), and
 L_1 - dphi' = c/tau - 1.
Remark S2.7.3 (F-form of Row 2). With Psi_1 := tau (D + omega) F = z tau (1 - tau) F_tau + (c + z - (1 + z) tau) F, a direct expansion gives
 z (1 - tau) dPsi_1/dtau - Psi_1 + x F = (1 - tau) z^2 Gauss[F]; so Row 2 is equivalent to z (1 - tau) dPsi_1/dtau = Psi_1 - x F (S2_checks.pl part D,
 to tau^22, for x in {2, 3/5, -3, -7/4} and several z, both signs).
Remark S2.7.4 (the two thimbles are the two Frobenius solutions at tau = 0; a remark, used nowhere). Let lambda = a/tau + (element of O) with a in C,
 and let f in R have tau-adic valuation r and leading coefficient f_r in K (non-zero). Then, by the definition of nabla, the operator
 Pi := nabla^2 - (c/tau) nabla + x/tau maps f eps_lambda to (z rho (z (rho - 1) - c) f_r tau^(r-2) + O(tau^(r-1))) eps_lambda, with rho := r + a/z:
 nabla(f eps_lambda) = (z rho f_r tau^(r-1) + O(tau^r)) eps_lambda, applied twice, and the term (x/tau) f is of order tau^(r-1). So a solution has rho = 0 or
 rho = 1 + c/z (on E_0 = R this is the indicial equation z^2 r (r - 1) - c z r = 0 of the Picard-Fuchs operator, exponents 0 and c/z + 1), and an
 isomorphism of (5.1) by a unit of K[[tau]] with constant term 1 changes neither r nor a. Now:
 - S_0 corresponds in E_mu, mu = c/tau - 1 (a = c), to kappa_1 tau F eps_mu (5.3), with F a unit of K[[tau]] (r = 1): rho = 1 + c/z.
 - S^(2)_0 := nu_2 r^pt_{0,2} eps_{L_2} (S2.10) is annihilated by Pi (proof of Prop S2.10.4, an identity in M((z)) between elements whose z-coefficients
   have bounded pole orders, hence also one in R by (0.3)), lies in E_{L_2} with L_2 in O (a = 0), and
   nu_2 r^pt_{0,2} is in O[[z]] with tau^0 coefficient nu_2(0) r^pt_{0,2}(z; 0), whose z^0 term nu_2(0) is non-zero (Lemma S2.10.1(d), (e); r = 0): rho = 0.
 So the singular thimble carries the exponent c/z + 1 and the regular one the exponent 0, as anticipated in conifold_regularity_note.md s.6.
 The same computation for Gauss's operator on E_lambda (5.1), where d/dtau (f eps_lambda) = (rho f_r tau^(r-1) + O(tau^r)) eps_lambda, gives the indicial
 polynomial rho (rho - 1 + gamma), exponents 0 and 1 - gamma: S1's Gauss-normalised singular expansion M = A F (an element of E_0 = R with r = 0) has
 rho = 0, and the regular expansion LE_2 of Lemma S2.10.3 has rho = -1 - c/z = 1 - gamma (Remark S2.10.3'(ii)).


## S2.8 The operator A

Prop S2.8. (a) A is the transpose of the matrix M_H = {{0, -x/tau}, {1, c/tau}} of H * in the basis (1, H) (columns: coordinates of H * 1, H * H).
 For a covector s with components s_b = s(phi_b), phi_0 = 1, phi_1 = H, (A s)_b = s(H * phi_b): A acts on covectors by s -> s(H * .).
 (b) chi_i := (1, L_i), i.e. the covector a_0 1 + a_1 H -> a_0 + a_1 L_i (the character of the algebra at idempotent i; chi_i(e_j) = delta_ij),
 satisfies A chi_i = L_i chi_i.
 (c) A is self-adjoint for < , >: eta^-1 A = {{-x^3/tau, c x^2/tau}, {c x^2/tau, x^2 - c^2 x/tau}} is symmetric, so <A a, b> = <a, A b>.
 (d) <chi_i, chi_j> = delta_ij Delta_i.
 (e) A - L_1 = (L_2 - L_1) P_2, where P_2 a := (pi(a)/Delta_2) chi_2 and
   pi(a) := <a, chi_2> = x^2 L_2 a_0 + (x^2 - c x L_2) a_1 = x^3 s* (a_0 - tau s* a_1).
 P_2 is the < , >-orthogonal projector onto C chi_2 along C chi_1, <P_2 a, b> = pi(a) pi(b)/Delta_2, and
   <(A - L_1) a, b> = (L_2 - L_1) pi(a) pi(b)/Delta_2 = -(x^3/tau) (a_0 - tau s* a_1)(b_0 - tau s* b_1).
Proof. (a) H * 1 = H and H * H = -(x/tau) 1 + (c/tau) H (Lemma S2.1(b)) give M_H, and (A s)_0 = s_1 = s(H * 1), (A s)_1 = -(x/tau) s_0 + (c/tau) s_1
 = s(H * H).
 (b) A (1, L) = (L, -x/tau + (c/tau) L), which is L (1, L) iff L^2 - (c/tau) L + x/tau = 0 (S2.2(d)). [Conceptually: chi_i is multiplicative, so
 (A chi_i)(a) = chi_i(H * a) = L_i chi_i(a).] chi_i(e_j) = delta_ij by Prop S2.3(a).
 (c) Multiply out with eta^-1 = {{0, x^2}, {x^2, -c x}}. [Conceptually: H * is eta-self-adjoint (Lemma S2.1(d)), M_H^T eta = eta M_H, hence
 eta^-1 M_H^T = M_H eta^-1, i.e. eta^-1 A is symmetric.]
 (d) Directly: <chi_i, chi_j> = x^2 (L_i + L_j) - c x L_i L_j. For i != j this is x^2 c/tau - c x x/tau = 0; for i = j it is 2 x^2 L_i - c x L_i^2
 = x L_i (2x - c L_i) = Delta_i (Prop S2.3(c)). Conceptually: a * e_j = chi_j(a) e_j (true for a = 1 and a = H, and linear in a), so
 a = Sum_j a * e_j = Sum_j chi_j(a) e_j, hence eta(a, e_i) = chi_i(a) eta(e_i, e_i), i.e. chi_i = eta(Delta_i e_i, .);
 the pairing of covectors given by eta^-1 is eta on the corresponding vectors, so <chi_i, chi_j> = eta(Delta_i e_i, Delta_j e_j) = delta_ij Delta_i.
 (e) Delta_1, Delta_2 are non-zero (Prop S2.3(d)) and chi_1, chi_2 is a < , >-orthogonal basis (L_1 != L_2). Hence a = Sum_i (<a, chi_i>/Delta_i) chi_i,
 i.e. 1 = P_1 + P_2 with P_i a := (<a, chi_i>/Delta_i) chi_i, and A = L_1 P_1 + L_2 P_2 by (b), so A - L_1 = (L_2 - L_1) P_2; <P_2 a, b> = (pi(a)/Delta_2)
 <chi_2, b> = pi(a) pi(b)/Delta_2. For pi: x^2 L_2 = x^3 s* and x^2 - c x L_2 = -x^3 tau s*^2 (Prop S2.3(d)). Last formula: (L_2 - L_1)/Delta_2
 = (-delta/tau)/(x^3 delta s*^2) = -1/(x^3 tau s*^2). QED
 Note that pi(chi_1) = x^3 s* (1 - tau s* L_1) = 0 is literally the formula L_1 = 1/(tau s*).


## S2.9 Step 6: the symplectic constant, K_1 = 1

Prop S2.9. Assume (S1.a), (S1.b).
 (a) K_1(z) := <v_1(z), v_1(-z)>/Delta_1 = <w(z), w(-z)> lies in C[[z, tau]] and is independent of t: K_1 is in C[[z]].
 (b) sigma(psi) := rho v_1(-psi) = (g(psi), (L_1 - psi/t) g(psi) - psi g(psi)'), where g(psi) := G(-psi) and ' = d/dt. Moreover t sigma(psi) is in
 C[[psi, t]]^2 and
   t sigma(psi)|_{t=0} = (0, c),  i.e.  [t^-1] sigma(psi) = (0, c);   and  g(psi)|_{t=0} = c/(c - psi)   (F2).
 (c) K_1(z) = 1, i.e. <v_1(w), v_1(-w)> = Delta_1 to all orders in w: (V6) for i = j = 1.
 (d) Consequences: <sigma(psi), sigma(-psi)> = rho^2 Delta_1 = -x c^3/t^2. In the edge formula (rc 147-148) E_11(z, w) = [P(z, w) - P(-w, w)]/(z + w),
 P(z, w) = -<v_1(-z), v_1(-w)>, the subtracted term is P(-w, w) = -Delta_1, so
   rho^2 E_11(psi, psi') = [rho^2 Delta_1 - <sigma(psi), sigma(psi')>]/(psi + psi'),   rho^2 Delta_1 = -x c^3/t^2;
 and, in C((psi))((t))^2,
   d sigma/dt = -sigma/t - (1/psi) (A - L_1) sigma.
Proof. (a) w = nu_1 v_1 is in C[[z, tau]]^2 (Prop S2.7) and nu_1^2 = 1/Delta_1, so K_1 = <w(z), w(-z)> is in C[[z, tau]]. Applying refl to z w' = (A - L_1) w
 gives -z w(-z)' = (A - L_1) w(-z) (A and L_1 do not involve z). eta^-1 is constant (Lemma S2.1(a)), so d/dt is a derivation of the pairing, and
   K_1' = <w(z)', w(-z)> + <w(z), w(-z)'> = (1/z) [ <(A - L_1) w(z), w(-z)> - <w(z), (A - L_1) w(-z)> ] = 0
 by Prop S2.8(c). Thus (1 - tau) dK_1/dtau = 0; 1 - tau is a unit, and the kernel of d/dtau on C[[z, tau]] is C[[z]].
 [Module form: refl S = nu_1 v_1(-z) eps_{-L_1} satisfies nabla(refl S) = -refl(nabla S) = -A refl S, so z <S, refl S>' = <A S, refl S> - <S, A refl S> = 0:
 "<s(z), s(-z)> is t-independent"; and <S, refl S> = nu_1^2 <v_1(z), v_1(-z)> = K_1.]
 (b) G = rho r_{0,1} is in C[[z, t]] with G(z; 0) = c/(c + z) (Prop S2.4(d)). By (H) and (5.2.1) with nu_1 = t rho/kappa_0,
   rho r_{1,1} = rho nu_1^-1 (z d/dt + L_1)(nu_1 r_{0,1}) = t^-1 (z d/dt + L_1)(t G) = (L_1 + z/t) G + z G'.
 By S2.2(d), t L_1 = c t/tau - x t s* is in O with value c at t = 0. So t rho r_{1,1} = (t L_1 + z) G + z t G' is in C[[z, t]], and at t = 0 it equals
 (c + z) c/(c + z) = c: the z-dependence cancels exactly, because of the geometric form c/(c + z) of F2. Substitute z = -psi.
 (c) w(z) = nu_1 v_1(z) = (t/kappa_0) rho v_1(z) = (t/kappa_0) sigma(-z), and refl(w)(z) = (t/kappa_0) sigma(z). Hence
   K_1(z) = <w(z), w(-z)> = (t^2/kappa_0^2) <sigma(-z), sigma(z)> = -(1/(x c^3)) <t sigma(-z), t sigma(z)>.
 By (b) the right side lies in C[[z, t]], with value at t = 0 equal to -(1/(x c^3)) <(0, c), (0, c)> = -(1/(x c^3)) c^2 eta^11 = -(1/(x c^3)) c^2 (-c x) = 1.
 By (a), K_1 is constant in t, so K_1 = K_1|_{t=0} = 1.
 (d) <sigma(psi), sigma(-psi)> = rho^2 <v_1(-psi), v_1(psi)> = rho^2 Delta_1 K_1(-psi) = rho^2 Delta_1 = -x c^3/t^2 (Prop S2.4(c)), and
 P(-w, w) = -<v_1(w), v_1(-w)> = -Delta_1 K_1(w) = -Delta_1; with z = psi, w = psi' this gives the formula for rho^2 E_11. Finally apply refl to
 z w' = (A - L_1) w and use refl(w)(z) = (t/kappa_0) sigma(z): -z (t sigma(z))' = t (A - L_1) sigma(z); rename z as psi and divide by -psi t. QED
 Only the t -> 0 data F2 and [t^-1] sigma = (0, c) enter the limit in (c); they are derived in (b) from (H) and (S1.b) alone, so S2 does not depend on
 Step 7 (S3). The t-independence (a) uses Row 2, i.e. the full closed form (S1.b).

Remark S2.9.1 (K_1 = 1 as a Wronskian identity for Gauss). By Prop S2.7(ii), v_1 = nu_1^-1 kappa_1 tau Ephi Phi with Phi := (F, (D + omega) F); moreover
 Ephi(z) Ephi(-z) = e^{-dphi/z} e^{dphi/z} = 1 and kappa_1(z) kappa_1(-z) = c^2/(kappa_0^2 (c^2 - z^2)) = -1/(x c (c^2 - z^2)). Hence
   K_1(z) = -tau^2 <Phi(z), Phi(-z)>/(x c (c^2 - z^2)),    so    K_1 = 1  <=>  tau^2 <Phi(z), Phi(-z)> = -x c (c^2 - z^2).
 The latter follows from Gauss[F] = 0 alone, without (S1.b): (D + omega) Phi = A Phi (Row 2 in F-form), refl gives (-D + omega(-z)) Phi(-z) = A Phi(-z), so by
 Prop S2.8(c) d/dt <Phi(z), Phi(-z)> = ((omega(-z) - omega(z))/z) <Phi(z), Phi(-z)> = -2 ((1 - tau)/tau) <Phi(z), Phi(-z)> and d/dt (tau^2 <Phi(z), Phi(-z)>) = 0;
 at tau = 0, tau Phi(z) -> (0, c + z), which gives -c x (c + z)(c - z). Given (S1.b) this is a second proof of (c) (S2_checks.pl part D: exact to tau^22
 at x in {2, 3/5, -3, -7/4} and several z).


## S2.10 The regular idempotent: second QDE row, orthogonality, K_2 = 1

Why a separate argument is needed. Row 1 at i = 2 is the definition (H). Row 2 at i = 2 cannot be obtained as in S2.7: the closed form (S1.b) concerns the
singular row only. Nor is the route of S1, Prop. 4.1 available, which multiplies the formal expansion by Exp[dphi/z] (f2(tau)/f2(0))^(-1/2) and expands at
tau = 0 in C((z))[[tau]]. At the regular critical point s~ -> oo, and phi_tau(s~) = -c log tau + (a function analytic at tau = 0), up to branch constants:
by Lemma S2.10.1(a) below, s~ = tau^-1 b with b := (c + delta)/(2x) in O, b(0) = c/x, so log s~ = -log tau + log b, log(1 - s~) = -log tau + log(tau - b), and
1 - tau s~ = U_2 is a unit. So e^{phi_tau(s~)/z} is tau^(-c/z) times a unit and is not an element of C((z))((tau)); d/dtau phi_tau(s~) = s~/(1 - tau s~) has a
simple pole; and f2 = phi''_tau(s~) has a double zero (Lemma S2.10.1). Instead we work coefficientwise in z at fixed small t != 0, in the z-adic ring M((z)) (0.3), and
formulate Gauss's equation in the module E^M_{lambda_2} of (5.1); its symbol eps_{lambda_2} records only the t-derivative lambda_2 of the critical value.
The proof of Row 2 at U_2 (Lemmas S2.10.2, S2.10.3, Prop S2.10.4) never expands anything at tau = 0; analyticity at tau = 0 is used only through units
of O and values at t = 0. Put r0 := min(1, |c|^2/(4|x|)) and D'' := {tau : 0 < |tau| < r0}. On |tau| < |c|^2/(4|x|),
delta = c (1 - 4 x tau/c^2)^(1/2) is analytic and delta/c has positive real part (S1, 1.2), so c + delta has no zero there.

Lemma S2.10.1 (the regular critical point). For tau in D'' (i.e. for small t != 0):
 (a) s~ = (c + delta)/(2 x tau) = 1/(x tau s*) is a critical point of phi_tau, and s~ is not in {0, 1, 1/tau}; tau s~ = (c + delta)/(2x) is in O, with value
     c/x at 0; 2 x tau s~ - c = delta; x s~ = L_1; U_2 = 1 - tau s~ = (x - 1 - delta)/(2x) = 1 - 1/(x s*) is a unit of O with U_2(0) = -1/x.
 (b) f2 := phi''_tau(s~) = delta/(x s~^2 U_2^2) = tau^2 Delta_2/(x^2 U_2^2). Hence f2 = tau^2 f2r with f2r := Delta_2/(x^2 U_2^2) a unit of O,
     f2r(0) = x^3/c; in particular s~ is a nondegenerate critical point for every tau in D''.
 (c) f0(tau) := phi_tau(s~(tau)) (any branch of the logarithms) has f0_tau := df0/dtau = s~/(1 - tau s~), an element of tau^-1 O with leading term -c/tau, and
       lambda_2 := df0/dt = (1 - tau) s~/(1 - tau s~) = 1 - x s~ = 1 - L_1,     mu + lambda_2 = L_2     (mu = c/tau - 1 = x + c (1 - tau)/tau, (5.3)).
     Correspondingly, u_2 := W(U_2) satisfies du_2/dt = L_2 and u_2 = x t + c log tau + f0 + const.
 (d) Fix a square root f2r^(1/2) in O (it exists since f2r(0) != 0; which of the two roots is irrelevant, as for kappa_0 in (5.2)), and put
     f2^(1/2) := tau f2r^(1/2) and
       nu_2 := tau/(x U_2 f2^(1/2)) = 1/(x U_2 f2r^(1/2)).
     Then nu_2 is a unit of O, nu_2^2 = 1/Delta_2, nu_2'/nu_2 = -Delta_2'/(2 Delta_2), and for every f in M((z)), in E^M_{L_2},
       nabla(nu_2 f eps_{L_2}) = nu_2 (L_2 + z D_2) f eps_{L_2}.                                                      (10.1)
 (e) (S4, 10.1-10.2.) v_2 is in O[[z]]^2: every r^pt_{0,2}[m], r_{0,2}[m] and r_{1,2}[m] is analytic at t = 0.
Proof. (a) By (0.2), s* s~ = 1/(x tau) and s* = 2/(c + delta), so s~ = (c + delta)/(2 x tau) and 2 x tau s~ - c = delta; s~ is the second root of the critical
 equation, i.e. a critical point of phi_tau (Lemma S2.2(c), (e)). s~ != 0 since s* s~ = 1/(x tau). If s~ = 1 or s~ = 1/tau, the critical equation
 x tau s^2 - c s + 1 = 0 gives x (tau - 1) = 0, resp. (tau - 1)/tau = 0, impossible for |tau| < 1. U_2 = 1 - (c + delta)/(2x) = (x - 1 - delta)/(2x); U_2 = 0
 would force delta^2 = (x - 1)^2, i.e. 4 x (1 - tau) = 0 (as c^2 - (x - 1)^2 = 4x), impossible; U_2(0) = (x - 1 - c)/(2x) = -1/x, and U_2 is analytic on
 |tau| < |c|^2/(4|x|), hence a unit of O. tau s~ = 1/(x s*) and x s~ = 1/(tau s*) = L_1 by Lemma S2.2(d).
 (b) The critical equation 1 - c s + x tau s^2 = 0 can be written 1 - s = x s (1 - tau s). So x/(1 - s~)^2 = 1/(x s~^2 U_2^2) and, as in S1 (1.3),
   U_2^2 f2 = [-(1 - tau s~)^2 - 1/x + tau^2 s~^2]/s~^2 = (2 tau s~ - c/x)/s~^2 = delta/(x s~^2)
 (-(1 - tau s)^2 + tau^2 s^2 = 2 tau s - 1 and 1 + 1/x = c/x), so f2 != 0 on D''. By Prop S2.4(a) at U_2 and Prop S2.3(c), tau^2 Delta_2/x^2 = tau^2 W2(U_2)
 = U_2^2 f2. Delta_2 and U_2 are units of O (Prop S2.3(e) and (a)), so f2r is one, with f2r(0) = (x^3/c)/(x^2 x^-2) = x^3/c.
 (c) phi'_tau(s~) = 0, so df0/dtau = (partial_tau phi_tau)(s~) = s~/(1 - tau s~) (as in Lemma S2.2(e); a change of branch changes f0 by a constant). tau s~ is
 in O with value c/x and U_2 is a unit with value -1/x, so f0_tau = tau^-1 (tau s~)/U_2 has leading term -c/tau. lambda_2 = (1 - tau) f0_tau, and
 (1 - tau) s = (1 - x s)(1 - tau s) holds at both roots of the critical equation (expand the right side and use x tau s^2 = c s - 1). So lambda_2 = 1 - x s~ = 1 - L_1,
 and mu + lambda_2 = c/tau - 1 + 1 - L_1 = L_2 (Lemma S2.2(d)). For u_2: u_2' = (partial_t W)(U_2) + W_U(U_2) dU_2/dt = L(U_2) = L_2 by Lemma S2.2(b), and the
 formula W(1 - tau s) = x t + c log tau + phi_tau(s) + x log(-1) of Lemma S2.2(e), at s = s~, gives u_2 = x t + c log tau + f0 + const, i.e. again
 L_2 = x + c (1 - tau)/tau + lambda_2.
 (d) nu_2^2 = tau^2/(x^2 U_2^2 f2) = 1/Delta_2 by (b), and nu_2 is a unit since U_2 and f2r are. The logarithmic derivative of nu_2^2 Delta_2 = 1 gives
 nu_2'/nu_2 = -Delta_2'/(2 Delta_2), and then nabla(nu_2 f eps_{L_2}) = (z nu_2 f' + z nu_2' f + L_2 nu_2 f) eps_{L_2} = nu_2 (L_2 + z D_2) f eps_{L_2}.
 (e) Let theta := U d/dU. theta log U = 1, theta log(1 - U) = 1 - (1 - U)^-1, theta log(1 - Q U) = 1 - (1 - Q U)^-1, and
 theta (1 - a U)^-n = n ((1 - a U)^-(n+1) - (1 - a U)^-n) for a in {1, Q}. So for k >= 1, W_k := theta^k W is a polynomial in (1 - U)^-1 and (1 - Q U)^-1 with
 coefficients in Z[x] (W_1 = x - x (1 - Q U)^-1 - (1 - U)^-1). At U = U_2(t) both are in O, since 1 - U_2 and 1 - Q U_2 both have the value c/x != 0 at t = 0;
 and W_2(U_2) = Delta_2/x^2 is a unit of O (Prop S2.3(c), (e)). By rc 121-123 and rmatrixAt0 (rc 139: no normalisation constants at alpha0 = 0),
 ratCN[0, 2][m] = cconj[wickCN[m]], whose value at U_1 is the value of wickCN[m] at U_2; so r^pt_{0,2}[m] is wickZ[[m+1]] with Wk[k] -> W_k(U_2), a polynomial in
 W_2(U_2)^-1, W_3(U_2), ..., W_{2m+2}(U_2), hence in O. The twist by h_w(z)^-1 is t-independent, so r_{0,2}[m] is in O, and by (H)
 r_{1,2}[m] = L_2 r_{0,2}[m] + D_2 r_{0,2}[m-1] is in O (L_2 in O, Delta_2 a unit of O). QED

Lemma S2.10.2 ((A2): the regular row in the Euler coordinate). Let Lr := Lnorm[phi_tau, (1 - tau s)^-1; s~] (S1, 2.2). For every tau in D'', coefficientwise in z,
   r^pt_{0,2}(z; t) = Lhat[phi_tau, (1 - tau s)^-1; s~] = U_2 Lr.                                                         (10.2)
 Hence Lr = U_2^-1 r^pt_{0,2} is in O[[z]] (Lemma S2.10.1(a), (e)).
Proof. This is S1, Prop. 3.1, steps (1)-(2), at the second root. Those steps use only that the root U is a nondegenerate critical point with U != 0 and that
 the change of coordinate has non-zero derivative at the critical point; at U_2 all of this holds for every tau in D''.
 (1) Fix tau in D'' and a branch y2 of log U_2 (U_2 != 0). Wy(y) := W(e^y) is analytic near y2, with Wy^(k)(y2) = W_k(U_2), Wy'(y2) = W_1(U_2) = 0 and
 Wy''(y2) = W_2(U_2) = Delta_2/x^2 != 0; the integrand e^{W/z} dU/U is e^{Wy(y)/z} dy, with amplitude 1. By Lemma S2.10.1(e) and S1 (D1),
 r^pt_{0,2} = Lnorm[Wy, 1; y2] = Lhat[Wy, 1; y2], to every order in z.
 (2) h(s) := log(1 - tau s), with the branch for which h(s~) = y2, is analytic near s~, and h'(s~) = -tau/U_2 != 0 (tau != 0). By S1 (1.1), valid for every s off
 {0, 1, 1/tau} with compatible branches, Wy o h = W(1 - tau s) = phi_tau + const near s~, and the amplitude 1 becomes (1 o h) h' = -tau (1 - tau s)^-1.
 S1 (P1)(iii) [Prop. 2.4], (D3) (a constant added to the phase changes nothing) and the linearity of Lnorm in the amplitude over functions of tau give
   Lhat[Wy, 1; y2] = Lhat[phi_tau, -tau (1 - tau s)^-1; s~] = Lhat[phi_tau, (1 - tau s)^-1; s~] = Lnorm[phi_tau, (1 - tau s)^-1; s~]/U_2^-1 = U_2 Lr,
 since Lhat[f, q; u0] = Lnorm[f, q; u0]/q(u0) (S1, 2.2, where the amplitude is called g) and (1 - tau s~)^-1 = U_2^-1. QED
 Both sides of (10.2) are, coefficientwise in z, holomorphic on D'': the left side is a rational function of U_2(Q) and Q, the right side a polynomial in
 1/f2, phi_tau^(k)(s~) = (k - 1)! [(-1)^(k-1) s~^-k - x (1 - s~)^-k + tau^k U_2^-k] and tau^a U_2^(-a-1) (S1 (D2)).

Lemma S2.10.3 ((B2): Gauss's equation at the regular critical point). Let
   Dop := d/dtau + f0_tau/z - f2_tau/(2 f2)       (f0_tau, f2 as in Lemma S2.10.1; f2_tau := df2/dtau),
 acting coefficientwise in z on M((z)), and alpha = beta = 1/z + 1, gamma = c/z + 2. Then
   [tau (1 - tau) Dop^2 + (gamma - (alpha + beta + 1) tau) Dop - alpha beta] Lr = 0      in M((z)).                     (10.3)
 Equivalently, with lambda_2 = 1 - x s~ and d/dtau, Gauss as in (5.1) on E^M_{lambda_2}: the element
   LE_2 := f2^(-1/2) Lr eps_{lambda_2}   satisfies   d/dtau LE_2 = f2^(-1/2) (Dop Lr) eps_{lambda_2}   and   Gauss[LE_2] = 0.       (10.4)
Proof. Fix tau0 in D'' and put s0 := s~(tau0). Then s0 is not in {0, 1, 1/tau0} (Lemma S2.10.1(a)), so f(tau, s) := phi_tau(s) (with branches of the logarithms
 chosen near (tau0, s0)) and the amplitude q_0(tau, s) := (1 - tau s)^-1 (called g in S1, 4.1) are analytic near (tau0, s0); s~(tau) is analytic with s~(tau0) = s0,
 d_s f(tau, s~(tau)) = 0, and f2(tau) = d_s^2 f(tau, s~(tau)) != 0 (Lemma S2.10.1(b)). So S1, Prop. 2.3 (P2) applies near tau0. Its f0' is
 (d_tau f)(tau, s~(tau)) = s~/(1 - tau s~) = f0_tau, and with delta_tau q := d_tau q + q (d_tau f)/z (S1, 4.1; d_tau f = s/(1 - tau s)) it says
   Dop Lnorm[phi_tau, q; s~] = Lnorm[phi_tau, delta_tau q; s~]
 near tau0, for every amplitude q analytic near (tau0, s0) and polynomial in z, 1/z. Apply it to q = q_0 and to q = delta_tau q_0 = (1 + 1/z) s (1 - tau s)^-2
 (analytic near (tau0, s0) because U_2(tau0) != 0); then Dop^2 Lr = Lnorm[phi_tau, delta_tau^2 q_0; s~] with delta_tau^2 q_0 = (1 + 1/z)(2 + 1/z) s^2 (1 - tau s)^-3.
 By the linearity of Lnorm over functions of tau and over C[z, 1/z] (S1 (D3)),
   [tau (1 - tau) Dop^2 + (gamma - (alpha + beta + 1) tau) Dop - alpha beta] Lr
        = Lnorm[phi_tau, tau (1 - tau) delta_tau^2 q_0 + (gamma - (alpha + beta + 1) tau) delta_tau q_0 - alpha beta q_0; s~].
 By the amplitude identity S1 (4.1) -- an identity of Laurent polynomials in z with coefficients rational in (s, tau), valid identically in (s, tau) and not
 only at a critical point: after multiplication by (1 - tau s)^3/(1 + 1/z) both sides equal (1/z)(-x tau s^2 + c s - 1) + (2 s - tau s - 1) -- the amplitude
 on the right is -(1 + 1/z)(k' + phi'_tau k/z) with k := s (1 - s)(1 - tau s)^-2 (' = d/ds), and k is analytic near s0. Its Lnorm vanishes by S1, Prop. 2.2 (P3).
 So every z-coefficient of the left side of (10.3) vanishes near tau0. As tau0 in D'' was arbitrary, it vanishes on D''. The z-coefficients of Lr, Dop Lr and
 Dop^2 Lr are finite combinations of the functions listed after Lemma S2.10.2 and of f0_tau, f2_tau/f2 and their tau-derivatives, all algebraic in tau and
 delta, holomorphic on D'' and meromorphic at 0, i.e. elements of M; so (10.3) holds in M((z)).
 (10.4): f2^(1/2) = tau f2r^(1/2) is in M and (f2^(-1/2))_tau = -(f2_tau/(2 f2)) f2^(-1/2); and lambda_2/(z (1 - tau)) = f0_tau/z by Lemma S2.10.1(c). By the
 definition of d/dtau on E^M_{lambda_2} (5.1), for every Y in M((z)),
   d/dtau (f2^(-1/2) Y eps_{lambda_2}) = [ (f2^(-1/2) Y)_tau + (f0_tau/z) f2^(-1/2) Y ] eps_{lambda_2} = f2^(-1/2) (Dop Y) eps_{lambda_2}.
 With Y = Lr this is the first claim; with Y = Lr and Y = Dop Lr, (d/dtau)^2 LE_2 = f2^(-1/2) (Dop^2 Lr) eps_{lambda_2}, so
 Gauss[LE_2] = f2^(-1/2) ([tau (1 - tau) Dop^2 + (gamma - (alpha + beta + 1) tau) Dop - alpha beta] Lr) eps_{lambda_2} = 0 by (10.3). QED
Remark S2.10.3'. (i) This is the computation of S1, Prop. 4.1 at the other critical point; what changes is only the last step. There the conjugating factor
 E = Exp[dphi/z] (f2(tau)/f2(0))^(-1/2) is an element of C((z))[[tau]] and Gauss's equation is then read in that ring; here the corresponding factor
 e^{f0/z} f2^(-1/2) is not (e^{f0/z} ~ tau^(-c/z), f2^(-1/2) ~ tau^-1), so it is kept as the symbol f2^(-1/2) eps_{lambda_2}, with f2^(-1/2) in M.
 (ii) The z^-2 coefficient of (10.3) is Lr[0] times tau (1 - tau) f0_tau^2 + (c - 2 tau) f0_tau - 1, the eikonal polynomial of S1, 7.4 at S' = f0_tau. It vanishes
 because (1 - tau) times it equals x (x tau s^2 - c s + 1) at s = s~ (S1, 7.4: the computation there uses only the critical equation). So f0_tau = s~/U_2 is
 the second root of the eikonal equation, the root -c/tau + O(1) that S1, 7.4 associates with the exponent 1 - gamma. In the language of Remark S2.7.4
 (LE_2 also lies in E_{lambda_2} over R, since f2^(-1/2) Lr is in tau^-1 O[[z]] by Lemma S2.10.2): lambda_2 = 1 - L_1 = -c/tau + O(1), and f2^(-1/2) Lr has
 tau-order exactly -1 (its tau^-1 coefficient f2r(0)^(-1/2) U_2(0)^-1 r^pt_{0,2}(z; 0) has z^0 term -x f2r(0)^(-1/2) != 0), so LE_2 has rho = -1 - c/z = 1 - gamma,
 while S1's M = A F (Gauss-normalised singular expansion) has rho = 0. These remarks are not used.

Prop S2.10.4 (Row 2 of the QDE at the regular idempotent). In M[[z]], indeed in O[[z]] after multiplication by tau,
   (L_2 + z D_2) r_{1,2} = -(x/tau) r_{0,2} + (c/tau) r_{1,2},    equivalently    z D_2 r_{1,2} = L_1 (r_{1,2} - L_2 r_{0,2}),
 i.e. (D_2 r_{1,2})[m-1] = L_1 (r_{1,2}[m] - L_2 r_{0,2}[m]) for every m >= 1 ((V5) for i = 2; m <= 9 there). With Row 1 (= (H)): z D_2 v_2 = (A - L_2) v_2, and
 w_2 := nu_2 v_2, an element of O[[z]]^2, satisfies z w_2' = (A - L_2) w_2. Both forms hold in frame B and in the point frame.
Proof. We work in the point frame and in E^M_{L_2}. Put S^(2)_a := nu_2 r^pt_{a,2} eps_{L_2} (a = 0, 1), the precise form of "e^{u_2/z} Delta_2^(-1/2) r^pt_{a,2}"
 (du_2/dt = L_2, Lemma S2.10.1(c); the module records only L_2).
 (i) Row 1. By (10.1) and (H), nabla S^(2)_0 = nu_2 (L_2 + z D_2) r^pt_{0,2} eps_{L_2} = S^(2)_1.
 (ii) By (10.2) and the definition of nu_2, nu_2 r^pt_{0,2} = nu_2 U_2 Lr = (tau/(x f2^(1/2))) Lr. With the product of (5.1) and L_2 = mu + lambda_2
 (Lemma S2.10.1(c)),
   S^(2)_0 = (1/x) (tau eps_mu)(f2^(-1/2) Lr eps_{lambda_2}) = (1/x) P LE_2,      P := tau eps_mu in E^M_mu (mu = c/tau - 1 in M).
 (iii) By (S2.6.2) with Y = LE_2 in E^M_{lambda_2}, and (10.4),
   [nabla^2 - (c/tau) nabla + x/tau] S^(2)_0 = (1/x) ((1 - tau) z^2/tau) P Gauss[LE_2] = 0.
 With (i): nabla S^(2)_1 = nabla^2 S^(2)_0 = (c/tau) S^(2)_1 - (x/tau) S^(2)_0. By (10.1), nabla S^(2)_1 = nu_2 (L_2 + z D_2) r^pt_{1,2} eps_{L_2}. The map
 f -> nu_2 f eps_{L_2} is injective (nu_2 is a unit), so (L_2 + z D_2) r^pt_{1,2} = (c/tau) r^pt_{1,2} - (x/tau) r^pt_{0,2} in M((z)). Both sides are in M[[z]],
 and multiplication by the t-independent h_w(z)^-1 (which commutes with D_2 and with (H)) gives the frame-B form. c/tau - L_2 = L_1 and x/tau = L_1 L_2
 (Lemma S2.2(d)) give the second form; its z^0 coefficient is 0 = L_1 (r_{1,2}[0] - L_2 r_{0,2}[0]), i.e. (H) at m = 0. The matrix forms follow as in Prop S2.7:
 the components of (A - L_2) v_2 are r_{1,2} - L_2 r_{0,2} = z D_2 r_{0,2} and -(x/tau) r_{0,2} + (c/tau - L_2) r_{1,2} = z D_2 r_{1,2}, and z (nu_2 f)' = nu_2 z D_2 f.
 Membership: v_2 and w_2 are in O[[z]]^2 (Lemma S2.10.1(d), (e)), and tau L_1 is in O. QED

Prop S2.10 (orthogonality). Assume Prop S2.7 and Prop S2.10.4 (so: given S1, header input (i), (ii)). Then <v_1(w), v_2(-w)> = 0 to all orders in w.
Proof. Y(z) := <v_1(z), v_2(-z)> is in M[[z]] (v_1, v_2 in M[[z]]^2; Prop S2.7 holds in M[[z]] by TRANSFER). Using v_i' = D_i v_i + (Delta_i'/(2 Delta_i)) v_i,
 the two QDEs (the second with z -> -z: -z D_2 v_2(-z) = (A - L_2) v_2(-z)) and Prop S2.8(c),
   z Y' = <(A - L_1) v_1(z), v_2(-z)> - <v_1(z), (A - L_2) v_2(-z)> + (z/2) (log(Delta_1 Delta_2))' Y = [(L_2 - L_1) + (z/2) (log(Delta_1 Delta_2))'] Y.
 Write Y = Sum_m Y_m z^m with Y_m in M. The z^0 coefficient gives (L_2 - L_1) Y_0 = 0, and for m >= 1 the z^m coefficient gives
 (L_2 - L_1) Y_m = Y'_{m-1} - (1/2) (log(Delta_1 Delta_2))' Y_{m-1}. Since L_2 - L_1 = -delta/tau is invertible in M, Y_0 = 0 and inductively Y_m = 0. QED
 This is the "e^{(u_1 - u_2)/z}" argument: (Delta_1 Delta_2)^(-1/2) Y eps_{L_1 - L_2} = <S^(1)(z), S^(2)(-z)> is nabla-flat, and a nabla-flat element f eps_lambda
 with f a formal power series in z and lambda = L_1 - L_2 invertible vanishes, since nabla(f eps_lambda) = 0 reads lambda f = -z f' (compare lowest z-orders).
 It is essential that f is a power series in z, bounded below: the symbol e^{-(u_1 - u_2)/z}, an infinite series in 1/z (and not even an element of R,
 since u_1 - u_2 contains c log tau), would formally give a non-zero flat section. (Y_0 = <chi_1, chi_2> = 0 is also direct.)
 Row 2 at i = 2 is genuinely needed here: Row 1 at i = 2 (the definition) only determines r_{1,2} from r_{0,2}; without Row 2 at i = 2 the derivative of
 v_2(-z) in Y' is not expressible through A and the recursion does not close. In fact the z^m coefficient of Y contains r_{0,2}[m] only through
 <chi_1, chi_2> r_{0,2}[m] = 0, so Y = 0 is a genuine constraint on the regular row, which the closed form of S1 (singular row only) cannot supply.
 Confirmations: (V5) for i = 2, m = 1..9, t^0..t^6 and (V6) for i != j, exact in Q, w^n for n <= 9 (Mathematica); S2_checks_reg.pl part P.

Prop S2.10.5 (the regular rows at t = 0; K_2). Assume Prop S2.10.4.
 (a) v_2(z; 0) = r_{0,2}(z; 0) (1, x/c), i.e. x r_{0,2}(z; 0) = c r_{1,2}(z; 0) (S4 Prop. 10.6(b), first proof).
 (b) K_2(z) := <v_2(z), v_2(-z)>/Delta_2 = <w_2(z), w_2(-z)> lies in C[[z]] (it is t-independent), and K_2(z) = r_{0,2}(z; 0) r_{0,2}(-z; 0).
 (c) Given S4 Prop. 10.6(a), r_{0,2}(z; 0) = h_w(z)^-2, one has K_2 = 1, i.e. <v_2(w), v_2(-w)> = Delta_2 to all orders in w. Together with Props S2.9(c) and
     S2.10: <v_i(w), v_j(-w)> = delta_ij Delta_i for all i, j, to all orders in w, i.e. (V6) in full.
Proof. (a) In the coefficient form of Prop S2.10.4, (D_2 r_{1,2})[m-1] = L_1 (r_{1,2}[m] - L_2 r_{0,2}[m]) (m >= 1), the left side is in O (Lemma S2.10.1(e);
 D_2 preserves O because Delta_2 is a unit of O). On the right, L_1 = c/t + l(t) with l in O (Lemma S2.2(d)) and r_{1,2}[m] - L_2 r_{0,2}[m] is in O; if the latter had
 a non-zero value at t = 0 the right side would have a simple pole. Hence r_{1,2}[m](0) = L_2(0) r_{0,2}[m](0) = (x/c) r_{0,2}[m](0) for m >= 1, and for m = 0
 by (H).
 (b) As in Prop S2.9(a): w_2 = nu_2 v_2 is in O[[z]]^2 and nu_2^2 = 1/Delta_2, so K_2 = <w_2(z), w_2(-z)> is in O[[z]]; refl applied to z w_2' = (A - L_2) w_2 gives
 -z w_2(-z)' = (A - L_2) w_2(-z), and K_2' = (1/z) [<(A - L_2) w_2(z), w_2(-z)> - <w_2(z), (A - L_2) w_2(-z)>] = 0 by Prop S2.8(c). An element of O with zero
 derivative is a constant, so K_2 = K_2|_{t=0} = <v_2(z; 0), v_2(-z; 0)>/Delta_2(0). By (a) this is r_{0,2}(z; 0) r_{0,2}(-z; 0) <(1, x/c), (1, x/c)>/Delta_2(0), and
 <(1, x/c), (1, x/c)> = x^2 (2x/c) - c x (x/c)^2 = x^3/c = Delta_2(0) (Prop S2.3(e)).
 (c) h_w(z)^-2 h_w(-z)^-2 = 1, because the exponent of h_w is odd in z. QED
 Remarks. (b) uses only Row 2 at U_2, and K_2 is the same in frame B and in the point frame (5.4). K_2 = 1 + O(z^2) is automatic (K_2 is even in z and
 K_2(0) = <chi_2, chi_2>/Delta_2 = 1); the vanishing of the higher coefficients is the content of S4 Prop. 10.6(a), an analytic input (Laplace's method for a
 Beta integral at t = 0 and Stirling's series), parallel to S1, Prop. 6.1 at the singular idempotent. K_2 = 1 is not used for Lemma C or Theorem R:
 E_22 and Dil_2 are regular whatever K_2 is (S4, 10.4). Confirmations: (V6) for i = j = 2, w^0..w^9 (Mathematica); S2_checks_reg.pl parts P and T
 (including (a) and r_{0,2}(z; 0) = h_w(z)^-2 to z^6).

Where Row 2 at U_2 and the orthogonality are used. Not in L1-L3 as stated (legs f_{j,b} = [w^b] E_1j(psi, w)): in Step 9 the term P_12(-w, w)
 = -<v_1(w), v_2(-w)> of E_12 enters only as a function of w alone; and not in Step 8, which uses only K_1 = 1 and the singular QDE (Prop S2.9(d)).
 They ARE used for E_21: (i) its regularity at t = 0 (part of F3) needs the order-t^-1 part, x r_{0,2}(z; 0) = c r_{1,2}(z; 0), which is Prop S2.10.5(a)
 (S4: [t^-1] rho E_21(z, w) = [f(z) - f(-w)]/(z + w) with f(z) = -c x (x r_{0,2}(-z; 0) - c r_{1,2}(-z; 0))); (ii) the transposition E_21(z, w) = E_12(w, z),
 i.e. E_21(b, a) = E_12(a, b), needs the full orthogonality (S4, Prop. 9.5: without it, E_21(z, w) - E_12(w, z) = [phi(z) - phi(-w)]/(z + w) with
 phi(u) = -<v_1(u), v_2(-u)>), and with it the leg seen by a singular vertex at the second end of an edge. Theorem R (ii) needs both, because graphValueC
 attaches EdN[c(v1), c(v2)][a_{v1}, a_{v2}] in the order of the edge list, so E_21 occurs whenever the regular vertex of a mixed edge is listed first.
 S4, Lemma 10.5 states Row 2 at U_2 and sketches the route of this section; Lemmas S2.10.2-S2.10.3 and Prop S2.10.4 supply the details, in particular
 (B2) = Lemma S2.10.3, which S4 asserts "together with Step 2(a),(b)".


## S2.11 Interface with Steps 7-10 (consequences used in S3, S4)

(I1) rho r_{1,1}(z) = (L_1 + z/t) G + z G' with L_1 = c/t + l(t), l in O, l(0) = c/2 - x/c (Prop S2.9(b), Lemma S2.2(d));
     sigma(psi) = (g, (L_1 - psi/t) g - psi g'), g(psi) = G(-psi), [t^-1] sigma = (0, c).
(I2) pi(sigma(psi)) = -(x^2 - c x L_2) (psi/t) d/dt (t g(psi)), with x^2 - c x L_2 = -x^3 tau s*^2 = O(t).
     Proof: pi(a) = x^3 s* (a_0 - tau s* a_1) (Prop S2.8(e)) and tau s* L_1 = 1 give pi(sigma) = x^3 s* [g - tau s* (L_1 - psi/t) g + tau s* psi g']
     = x^3 tau s*^2 psi (g/t + g') = x^3 tau s*^2 (psi/t) (t g)'.
(I3) Theta := (L_2 - L_1)(c x L_2 - x^2)^2/Delta_2 = -x^3 tau s*^2 = x c dphi' - x = -x^3 t/c^2 + O(t^2), and by Prop S2.8(e) and (I2)
       <(A - L_1) sigma(psi), sigma(psi')> = (L_2 - L_1) pi(sigma(psi)) pi(sigma(psi'))/Delta_2 = Theta (psi psi'/t^2) d/dt(t g(psi)) d/dt(t g(psi')).
     Proof: L_2 - L_1 = -delta/tau, c x L_2 - x^2 = x^3 tau s*^2, Delta_2 = x^3 delta s*^2; and x c (1 - x s*) - x = x^2 (1 - c s*) = -x^3 tau s*^2.
(I4) rho^2 E_11(psi, psi') = [-x c^3/t^2 - <sigma(psi), sigma(psi')>]/(psi + psi') and d sigma/dt = -sigma/t - (1/psi) (A - L_1) sigma (Prop S2.9(d)).
     With (I3) these give the E_11 flow (V8) of Step 8, written out there.
(I5) [= S4's (I5)] <v_1(w), v_2(-w)> = 0 to all orders in w (Prop S2.10, given Props S2.7 and S2.10.4).
(I6) [S4's Lemma 10.5 and Prop. 10.6(b)] z D_2 v_2 = (A - L_2) v_2 in O[[z]]^2 after multiplication by tau (Prop S2.10.4); v_2(z; 0) = r_{0,2}(z; 0) (1, x/c)
     (Prop S2.10.5(a)); K_2 = 1 given S4 Prop. 10.6(a) (Prop S2.10.5(c)).
S2_checks.pl part C checks (I1) (F2 and [t^-1] sigma_1 = c), the sigma form of Row 2, (I2) and K_1 = 1 (t^2 <sigma(psi), sigma(-psi)> = -x c^3) as exact
t-series to t^22 from the closed form (S1.b), at x in {2, 3/5, -3} and several rational psi; part A checks (I3) and Prop S2.8(e) at rational points.
S2_checks_reg.pl checks (I5) and (I6) without using any closed form: part P builds the code's rows at both critical points as exact jets at four points
tau0 != 0 and tests (A2), (B2), Row 2 at both idempotents and (V6) for all (i, j); part T builds the regular rows as t-series at t = 0 and tests Row 2 at U_2,
v_2(z; 0) = h_w(z)^-2 (1, x/c) and K_2 = 1 (details in S2.12).


## S2.12 Status

Each item: statement -- status -- independent confirmation.
 - Lemma S2.1 (eta, eta^-1, H * H = -(x/tau) 1 + (c/tau) H, Leq) -- proved -- S2_checks A.
 - Lemma S2.2 (critical points; L(U) = 1/(1 - U) = partial_t W; L_1 = 1/(tau s*) = c/tau - x s*, L_2 = x s*; L_1 + L_2 = c/tau, L_1 L_2 = x/tau;
   du_1/dt = L_1; dphi' = 1 - x s*, L_1 - dphi' = c/tau - 1) -- proved -- (V2), (V3) (d(dphi)/dt = 1 - x s*); S2_checks A, C, E.
 - Prop S2.3 (idempotents; Delta_i = x^3 (L_i - L_j)/(x - c L_j) = x L_i (2x - c L_i) = x^2 W2(U_i); Delta_1 = -x delta/(tau^2 s*^2),
   Delta_2 = x^3 delta s*^2; F1) -- proved -- (V1), (V2); S2_checks A; S2_checks_reg P, T.
 - Prop S2.4 (tensor rule W2 = (1 - tau s)^2 phi''_tau/tau^2 at critical points; alg^2 = delta/(c^3 s*^2); Delta_1 = -x c^3 alg^2/tau^2; rho = tau/(t alg))
   -- proved -- (V3), check (2c); S2_checks A, E (against leaf_structure.md and hvc_leaf_proof_check_out.txt).
 - KEY FORMULA G = rho r_{0,1} = (c/(c + z)) (tau/t) F Ephi -- proved given S1 (S1.b) -- (V3).
 - Lemma S2.6 = (V7), as an operator identity on R, on every E_lambda and every E^M_lambda, and its module form (S2.6.2) -- proved -- (V7) (hvc_pf_gauss.wls);
   S2_checks B; S2_checks_reg P (through Prop S2.10.4).
 - Prop S2.7 (Row 2 of the QDE at idempotent 1, all orders in z; z D_1 v_1 = (A - L_1) v_1) -- proved given S1 (S1.a, S1.b) -- (V5) for i = 1,
   m = 1..9, t^0..t^6; S2_checks C (sigma form, t^22), D (F-form, tau^22); S2_checks_reg P (code rows, jets).
 - Prop S2.8 (A = M_H^T; A chi_i = L_i chi_i; self-adjointness; <chi_i, chi_j> = delta_ij Delta_i; A - L_1 = (L_2 - L_1) P_2, <P_2 a, b> = pi(a) pi(b)/Delta_2;
   rank-one form of <(A - L_1) a, b>) -- proved -- S2_checks A.
 - Prop S2.9 (K_1 t-independent; F2 and [t^-1] sigma = (0, c); K_1 = 1; <sigma(psi), sigma(-psi)> = -x c^3/t^2; P(-w, w) = -Delta_1 in E_11; sigma flow)
   -- proved given S1 (S1.a, S1.b) -- (V6) for i = j = 1, w^0..w^9; S2_checks C (t^22), D (Wronskian form, tau^22); S2_checks_reg P.
 - Lemma S2.10.1 (s~, U_2 unit, f2 = tau^2 Delta_2/(x^2 U_2^2) = (x^3/c) tau^2 (1 + O(tau)), lambda_2 = 1 - x s~, mu + lambda_2 = L_2, nu_2 a unit of O, v_2 in O[[z]]^2)
   -- proved -- S2_checks_reg P (jets), T (t-series).
 - Lemma S2.10.2 ((A2): r^pt_{0,2} = U_2 Lnorm[phi_tau, (1 - tau s)^-1; s~]) -- proved given S1 (P1), (D1), (D3), (1.1) -- S2_checks_reg P (z^0..z^6, jets to
   order 4 at four points).
 - Lemma S2.10.3 ((B2): Gauss's equation for the regular expansion, in M((z)) and in E^M_{lambda_2}) -- proved given S1 (P2), (P3), (D3), (4.1) --
   S2_checks_reg P (z^-2..z^4, with a negative control: flipping the sign of f2_tau/(2 f2) makes it fail).
 - Prop S2.10.4 (Row 2 of the QDE at U_2, all orders in z) -- proved given Lemmas S2.10.2, S2.10.3 -- (V5) for i = 2, m = 1..9, t^0..t^6;
   S2_checks_reg P (code rows, jets at four points), T (t-series at t = 0 to t^8, with a negative control).
 - Prop S2.10 (<v_1(w), v_2(-w)> = 0 to all orders) -- proved given S1 -- (V6) for i != j, w^0..w^9; S2_checks_reg P.
 - Prop S2.10.5 ((a) v_2(z; 0) = r_{0,2}(z; 0) (1, x/c); (b) K_2 t-independent, = r_{0,2}(z; 0) r_{0,2}(-z; 0); (c) K_2 = 1) -- (a), (b) proved given S1;
   (c) proved given S4 Prop. 10.6(a) -- (V6) for i = j = 2, w^0..w^9; S2_checks_reg P, T (also r_{0,2}(z; 0) = h_w^-2 to z^6).
 - Interface (I1)-(I6), Theta = -x^3 tau s*^2 = x c dphi' - x -- proved ((I1), (I2), (I4) given S1; (I5), (I6) given S1, and K_2 = 1 also given S4 10.6(a)) --
   Theta as printed in hvc_leaf_proof_check_out.txt (S2_checks E); S2_checks A, C; S2_checks_reg P, T.
S2_checks_reg.pl in detail (output S2_checks_reg_out.txt, M = 6, J = 4, N = 8: 420 checks, 0 failures, 2.4 min). Part P (268 checks) at
 (x, tau0, s~(tau0)) = (2, 5/8, 2), (3/5, 16/27, 15/4), (-3, 5/12, 2), (1/3, 7/12, 6), exact jets in e = tau - tau0 to e^4 (e^2 after two derivatives):
 the data of Lemma S2.10.1, (V1) at both roots, (A2) for z^0..z^6, (B2) for z^-2..z^4 and its negative control, (V5) at both idempotents for z^1..z^6, and (V6)
 for all (i, j), w^0..w^6. Part T (152 checks) at x = 2, 3/5, -3, 1/3, t-series to t^8: L_2 = x s*, Delta_2 (the code's literal 1/(e_2^T eta e_2)) = x^2 W2(U_2)
 = x^3 delta s*^2 with Delta_2(0) = x^3/c and Delta_2'(0) = 0, f2 = Delta_2 tau^2/(x^2 U_2^2) = (x^3/c) t^2 + O(t^3), (V5) at U_2 for z^1..z^6 and its negative
 control (sign of Delta_2'/(2 Delta_2) flipped), r_{0,2}[m](0) = [z^m] h_w^-2 and r_{1,2}[m](0) = (x/c) [z^m] h_w^-2 for m <= 6, and K_2 = 1 for w^0..w^6.

Open points and cross-references for the referee.
 (O1) [resolved in this revision] Row 2 at the regular idempotent is proved in S2.10 (Lemmas S2.10.2, S2.10.3, Prop S2.10.4), coefficientwise in z at fixed
      small t != 0, with no expansion at tau = 0, from S1's formal Laplace calculus at the nondegenerate critical point s~ and the amplitude identity S1 (4.1).
      Consequently the orthogonality (Prop S2.10) holds unconditionally (given S1), and K_2 = r_{0,2}(z; 0) r_{0,2}(-z; 0) = 1 to all orders given S4 Prop.
      10.6(a) (Prop S2.10.5); so (V6) holds in full. The orthogonality and Row 2 at U_2 are used for E_21 (S4, Prop. 9.5 and 10.6(b)), hence in Theorem R (ii),
      not for L1-L3 or Step 8. K_2 = 1 is not used downstream.
 (O2) Consistency with S1 (checked against S1_thimble.md): the input (S1.b) is needed as an identity in K[[tau]] with both sides in C[[z, tau]], for the
      TRANSFER argument of S2.0. S1 states it in exactly this form: Theorem S1(c) in C((z))[[tau]] = K[[tau]], Theorem S1(b), (d) for the regularity in
      tau and in z. S1's formula W_2(U1) = U1^2 phi''_tau(s*)/tau^2 is Prop S2.4(a), since U_1 = 1 - tau s*. S2.10 uses S1's Props. 2.2-2.4 exactly as stated
      there (for an arbitrary nondegenerate critical point depending analytically on a parameter), at tau0 != 0 and s0 = s~(tau0).
 (O3) Suggested edits outside this file (not made here): S3, s.0.4 ("Logical order") should cite Prop S2.9(b) for [t^-1] sigma = (0, c); S1 should rename
      its phase constant kappa(z) of (1.2) (e.g. kappa_ph(z)), the name used here; S4, Lemma 10.5 can cite Prop S2.10.4 for (B2).
Suggested Mathematica checks (not run; licence):
 (M1) At x = 2 and x = 3/5: buildFactors[12, "B", -1, x] and then checks (4) and (5) of hvc_leaf_proof_check.wls ((V5) at both idempotents and (V6)) to z^12.
 (M2) Symbolic x and z: the Wronskian identity tau^2 <Phi(z), Phi(-z)> = -x c (c^2 - z^2), Phi = (F, (D + omega) F) with F the product formula, as a tau-series
      to tau^12 (Remark S2.9.1, the F-only proof of K_1 = 1).
 (M3) From hvc_factors_g4_Bm.mx: v_2(z; 0) = h_w(z)^-2 (1, x/c) to z^9 (Prop S2.10.5(a) with S4 10.6(a); with Prop S2.10.5(b) this is K_2 = 1 to z^9).
 (M4) Symbolic x: (10.3), i.e. [tau (1 - tau) Dop^2 + (gamma - (alpha + beta + 1) tau) Dop - alpha beta] Lr = 0 with Dop = d/dtau + s~/((1 - tau s~) z) - f2'/(2 f2),
      as Laurent series in tau at tau = 0 for z^-2..z^4, with Lr computed from the Wick formula (S1 (D1)) in the coordinate s at s~ = (c + delta)/(2 x tau)
      (Lr is in O[[z]] by Lemma S2.10.2; f0_tau = s~/(1 - tau s~) and f2'/f2 have simple poles at tau = 0).
