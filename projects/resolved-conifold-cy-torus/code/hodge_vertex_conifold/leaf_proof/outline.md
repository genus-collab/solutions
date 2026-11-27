# All-order proof of the leaf structure (Lemma C) -- outline and verified facts (2026-10-31)

No topological recursion is used anywhere. The route: the singular thimble is an Euler integral, so the singular identity row of the
frame-B R-matrix has an exact closed form in terms of Gauss 2F1; the H-row follows by the first row of the quantum differential
equation (QDE); the singular-singular edge E_11 obeys an exact t-flow whose right-hand side is a product of two copies of the
same leaf series; the E_12 legs follow by a rationality argument.

## 0. Conventions (from rc_lib.m, hvc_lib.m; read them)

- eps1 = x, eps2 = 1, c = 1 + x; torus alpha0 = 0, alpha1 = x, alpha2 = 1. Q = e^t, tau := 1 - e^-t = 1 - 1/Q (so dtau/dt = 1 - tau).
- LG model: W(U) = x Log[1 - Q U] - Log[U] + Log[1 - U] (spotU at alpha0 = 0), integrals Int e^{W/z} dU/U.
  Critical points: Q x U^2 + Q (1-x) U - 1 = 0. Singular root U1 -> 1 as t -> 0 (it sits between the log points U = 1 and U = 1/Q),
  regular root U2 -> -1/x.
- Flat basis {1, H}; eta = {{c/x^3, 1/x^2}, {1/x^2, 0}}, eta^-1 = {{0, x^2}, {x^2, -c x}}; on covectors <a, b> := Sum eta^{ab} a_a b_b
  = x^2 (a_0 b_1 + a_1 b_0) - c x a_1 b_1.
- H*-eigenvalue at a critical point: L(U) = -Q (x U + 1)/(1 - Q) (LofU); L_i = L(U_i).
- Point-frame identity row r^pt_{0,i}(z) = Sum_m w_m z^m: the Wick series (wickDirect in rc_lib.m) = the normalised formal
  stationary-phase (Laplace) expansion of Int e^{W/z} dU/U at U_i in the coordinate y = log U - log U_i, leading coefficient 1.
  At alpha0 = 0 there are no QRR normalisation constants (rmatrixAt0).
- Frame B (mode B, s = -1): every row is multiplied by h_w(z)^-1, h_w(z) = Exp[-Sum_k B_2k/(2k(2k-1)) z^(2k-1) Sum_a w_a^(1-2k)],
  w = (c, -x, -1) (twistRows[hCoeffs[w, -1, ...]] in hvc_lib.m). h_w is t-independent. One checks h_w(z) = Stirling tail of
  B(x/z, 1/z) = exp(Sum_k B_2k/(2k(2k-1)) z^(2k-1) (x^(1-2k) + 1 - c^(1-2k))).
- H-row (rmatrixCore): r_{1,i}(z) = L_i r_{0,i}(z) + z Delta_i^(1/2) d/dt (Delta_i^(-1/2) r_{0,i}(z)), d/dt = Q d/dQ (total).
- Delta_i = 1/eta(e_i, e_i) (DeltaCN), vertex weight Delta_i^(g-1).
- Dilaton leaves Dil_i(n) = (-1)^n r_{0,i}[n-1], so q_t(psi) := Sum_{d>=1} Dil_1(d; t) psi^d = -psi r_{0,1}(-psi; t).
- Edges (evalQ / EdN): E_ij(z, w) = Sum_{m,n} E_ij(m, n) z^m w^n = [P(z, w) - P(-w, w)]/(z + w),
  P(z, w) = -<v_i(-z), v_j(-w)>, v_i(z) := (r_{0,i}(z), r_{1,i}(z)).
- rho(t) := (Delta_1 t^2/(-x c^3))^(-1/2), lam_j(psi) := 1/(c - j psi), j >= 0 (lam_0 = 1/c is a constant).
- Legs f_{j,b}(psi) := Sum_a E_{1j}(a, b) psi^a = [w^b] E_1j(psi, w).

## 1. Statement to prove (Lemma C = Lemma L of leaf_structure.md)

For every k >= 0:
 (L1) [t^k] (rho q_t) = q_0 delta_{k0} + (element of span(lam_1, ..., lam_{k+1}) for k >= 1), q_0 = c - c^2 lam_1;
 (L2) [t^k] (rho f_{j,b}) in span(lam_1, ..., lam_{k+1}) for j = 1, 2 and every b >= 0;
 (L3) [t^k] (rho^2 E_11(psi, psi')) in span{lam_i(psi) lam_l(psi') : i, l >= 1, (i-1) + (l-1) <= k};
together with F1 (Delta_1 = -x c^3 t^-2 (1 + O(t)); Delta_2 regular, = x^3/c at t = 0), F2 (r_{0,1}(z; 0) = c/(c+z)) and F3
(all frame-B edge and dilaton factors regular at t = 0). This is exactly the input of Theorem R (proof_status_mixed_scale.md s.3)
besides Conjecture M.

## 2. Facts verified exactly (Mathematica; hvc_leaf_proof_check.wls -> hvc_leaf_proof_check_out.txt, hvc_pf_gauss.wls)

 (V1) Delta_1 = x^2 W2(U1), W2 := (U d/dU)^2 W -- an identity in the quadratic algebra (so also Delta_2 = x^2 W2(U2)).
 (V2) s* := (1 - U1)/tau satisfies x tau s*^2 - c s* + 1 = 0, s*(0) = 1/c; L_2 = x s*, L_1 = c/tau - x s*; L_1 + L_2 = c/tau,
      L_1 L_2 = x/tau; Delta_1 = x^3 (L_1 - L_2)/(x - c L_2), Delta_2 = x^3 (L_2 - L_1)/(x - c L_1). With delta := Sqrt[c^2 - 4 x tau]:
      s* = 2/(c + delta), x - c L_2 = -x^2 tau s*^2, Delta_1 = -x delta/(tau^2 s*^2).
 (V3) KEY FORMULA, checked against the stored factors for all t^k, k <= 6, and z^m, m <= 9:
        G(z; t) := rho(t) r_{0,1}(z; t) = (c/(c+z)) (tau/t) F(z, tau) Exp[-dphi(tau)/z],
        F(z, tau) := 2F1(1/z + 1, 1/z + 1; c/z + 2; tau) = Sum_n (tau/z)^n prod_{j=1}^n (1 + j z)^2 / (n! prod_{j=2}^{n+1} (c + j z)),
        dphi(tau) := phi_tau(s*) - phi_0(1/c),  phi_tau(s) := Log[s] + x Log[1 - s] - Log[1 - tau s].
      Equivalently r_{0,1} = (c/(c+z)) F Exp[-dphi/z] alg(tau) with alg := (1 - tau s*) (phi''_tau(s*)/phi''_0(1/c))^(1/2),
      phi''_0(1/c) = -c^3/x, and rho = tau/(t alg).
      dphi(t) = t/c - x^2 t^2/(2 c^3) + (1-x)^2 x^2 t^3/(6 c^5) - ...; d(dphi)/dt = 1 - x s*.
 (V4) [t^k] (z G) for k <= 12, from the closed form with x symbolic: regular at z = 0, limit c delta_{k0} at z = oo, denominator
      exactly prod_{j=1}^{k+1} (c + j z) (simple poles).
 (V5) Second QDE row, both idempotents, z^1..z^9, t^0..t^6: (D_i r_{1,i})[m-1] = L_j (r_{1,i}[m] - L_i r_{0,i}[m]), D_i = d/dt - Delta_i'/(2 Delta_i).
 (V6) Symplecticity, exact in Q, w^n for n <= 9: <v_i(w), v_j(-w)> = delta_ij Delta_i.
 (V7) PF <=> Gauss (exact, symbolic): with D := z d/dt = z (1 - tau) d/dtau,
        [D^2 - (c/tau) D + x/tau] ((1-tau)^(-x/z) tau^(c/z+1) F) = -((tau - 1) z^2/tau) (1-tau)^(-x/z) tau^(c/z+1) Gauss[F],
      Gauss[F] = tau(1-tau) F'' + (gamma - (alpha+beta+1) tau) F' - alpha beta F, alpha = beta = 1/z + 1, gamma = c/z + 2.
 (V8) E_11 FLOW, checked for all (m, n) with m + n <= 8 and t^0..t^6:
        d/dt [t^2 rho^2 E_11(z, w)] = Theta(t) d/dt(t g(z)) d/dt(t g(w)),   g(psi) := G(-psi),
        Theta := (L_2 - L_1)(c x L_2 - x^2)^2/Delta_2 = -x^3 tau s*^2 = x c (d(dphi)/dt) - x = -x^3 t/c^2 + O(t^2).

## 3. Proof outline (each step to be written rigorously)

Step 1 (rescaled thimble). v = 1 - U = tau s: e^{W/z} dU/U = (const) Q^{x/z} tau^{c/z+1} e^{phi_tau(s)/z} (1 - tau s)^-1 ds, where the
 constant is independent of s and t (a phase). The critical point s*(tau) of phi_tau is nondegenerate and analytic at tau = 0.
 Euler: J(z, tau) := Int_0^1 s^{1/z} (1-s)^{x/z} (1-tau s)^{-1/z-1} ds = B(1/z+1, x/z+1) F(z, tau).
Step 2 (formal Laplace lemma). (a) Coordinate invariance: for each fixed small t != 0, y = log U - log U1 is an analytic coordinate
 near s*, so r^pt_{0,1}(z; t) = N(z, tau) := the normalised formal Laplace expansion of J at s*(tau); N in C{tau}[[z]].
 (b) The formal expansion commutes with d/dtau and kills total s-derivatives, hence satisfies Gauss's equation (Euler's
 derivation of Gauss's equation is an integration by parts). (c) In C((z))[[tau]] the power-series solution of Gauss with value 1
 at tau = 0 is unique ((k+1)(k+gamma) F_{k+1} = (k+alpha)(k+beta) F_k with k + gamma invertible), so
 N = A(z) F(z, tau) Exp[-dphi/z] alg(tau), A(z) := normalised expansion of B(1/z+1, x/z+1) = h_w(z) c/(c+z) (Stirling; the
 Beta-function identity B(a+1,b+1) = B(a,b) ab/((a+b)(a+b+1)) gives the factor x/(c(c+z)) and the tail of B(x/z,1/z) is h_w).
 (d) N in C[[z]][[tau]] and A in 1 + z C[[z]], so every tau^k-coefficient of F Exp[-dphi/z] alg is a power series in z; being a
 rational function of z it has no pole at z = 0. Frame B: r_{0,1} = h_w^-1 N = (c/(c+z)) F Exp[-dphi/z] alg. (F2: tau = 0.)
Step 3 (metric). Delta_1 = x^2 W2(U1) = x^2 (1 - tau s*)^2 phi''_tau(s*)/tau^2 = -x c^3 alg^2/tau^2, so rho = tau/(t alg) and the KEY
 FORMULA (V3) holds. F1 follows.
Step 4 (L1). [tau^n] F has simple poles at z = -c/j (j = 2..n+1), a pole of order n at z = 0, and tends to 1/(n+1) at z = oo;
 [tau^m] Exp[-dphi/z] is a polynomial in 1/z. Hence [t^k](z G) is rational in z with at most simple poles at z = -c/j, j = 1..k+1
 (the j = 1 pole from c/(c+z)) and at z = 0; none at 0 by Step 2(d); bounded at oo with limit [t^k] c (tau/t) F(oo, tau) =
 [t^k] c (tau/t)(-log(1-tau)/tau) = c delta_{k0}. Partial fractions: [t^k](zG) in span{1/(c + j z): 1 <= j <= k+1}, no constant,
 for k >= 1. With psi = -z: rho q_t(psi) = -psi g(psi) and L1 follows. Equivalently: g_0 = c lam_1, and for k >= 1
 g_k := [t^k] g in span(lam_1..lam_{k+1}) with g_k = O(psi^-2) at oo.
Step 5 (QDE). u_1 := W(U1) has du_1/dt = L_1 (dW/dt = -x Q U/(1 - Q U) equals L(U) on the critical locus). Put
 s_a := e^{u_1/z} Delta_1^(-1/2) r_{a,1}. The H-row definition is z ds_0/dt = s_1. By Steps 1-3, s_0 = kappa(z) Q^{x/z} tau^{c/z+1}
 F(z, tau) (formal; kappa t-independent), so by (V7) s_0 satisfies D^2 s_0 - (c/tau) D s_0 + (x/tau) s_0 = 0, i.e. the second QDE
 row z ds_1/dt = -(x/tau) s_0 + (c/tau) s_1. So z ds/dt = A s with A = [[0, 1], [-x/tau, c/tau]] (covector form of H*), whose
 eigenvectors are (1, L_i) with eigenvalues L_i; A is self-adjoint for < , > (<(1,L_1), (1,L_2)> = x^2 (L_1+L_2) - c x L_1 L_2 = 0).
Step 6 (symplectic constant). d/dt <s(z), s(-z)> = 0 by Step 5, so K_1(z) := <v_1(z), v_1(-z)>/Delta_1 is t-independent. Its t -> 0
 limit is computed from sigma below: <sigma(psi), sigma(-psi)> = c^2 eta^{11}/t^2 + O(1/t) and rho^2 Delta_1 = -x c^3/t^2 give
 K_1 = c^2 eta^{11}/(-x c^3) = 1. (Also <v_1(w), v_2(-w)> = 0: a formal series P(z) with e^{(u_1-u_2)/z} P t-independent vanishes.)
Step 7 (H-row, sigma). rho r_{1,1}(z) = (L_1 + z/t) G + z dG/dt, L_1 = c/t + l(t) with l regular. In psi: sigma(psi) := rho v_1(-psi)
 = (g(psi), (L_1 - psi/t) g - psi dg/dt). [t^-1] sigma = (0, c). Orthogonality <(1,L_1),(1,L_2)> = 0 gives
 pi(sigma(psi)) := <sigma(psi), (1, L_2)> = -(x^2 - c x L_2) (psi/t) d/dt(t g(psi)), and x^2 - c x L_2 = -x^3 tau s*^2 = O(t).
Step 8 (E_11 flow and L3). By Step 6, rho^2 E_11(psi, psi') = [rho^2 Delta_1 - <sigma(psi), sigma(psi')>]/(psi + psi'), rho^2 Delta_1 = -x c^3/t^2.
 Step 5 gives d sigma/dt = -sigma/t - (1/psi)(A - L_1) sigma, and A - L_1 = (L_2 - L_1) P_2 with <P_2 a, b> = pi(a) pi(b)/Delta_2. Hence
   d/dt [t^2 rho^2 E_11] = t^2 <(A - L_1) sigma(psi), sigma(psi')>/(psi psi') = Theta(t) d/dt(t g(psi)) d/dt(t g(psi'))   (V8).
 Integration constant: [t^2 rho^2 E_11]_{t=0} = -x c^3 (K_1 - 1)/(psi + psi') = 0. Theta = O(t). Therefore
   [t^k] rho^2 E_11 = (1/(k+2)) Sum_{k0 >= 1, k0 + k1 + k2 = k+1} Theta_{k0} (k1+1) g_{k1}(psi) (k2+1) g_{k2}(psi'),
 which gives L3 with support (i-1) + (l-1) <= k1 + k2 <= k and no lam_0, F3 for E_11, and rank one at t^0:
 (Theta_1/2) c^2 lam_1 lam_1 = -(x^3/2) lam_1(psi) lam_1(psi'). L2 for j = 1: [w^b] lam_l(w) = l^b/c^(b+1), and rho^-1 = 1 + O(t).
Step 9 (E_12 legs). rho E_12(psi, w) = [-<sigma(psi), v_2(-w)> + (function of w alone)]/(psi + w). [w^b] of this equals
 Sum_{n=0}^b (-1)^n psi^(-n-1) [w^(b-n)](numerator). At order t^k the numerator is a combination of the two components of [t^k']
 sigma (each rational in psi with at most simple poles at psi = c/j, j <= k'+1 <= k+1, bounded at oo; the t^-1 part is the constant
 (0, c)) times functions of w. Hence [t^k][w^b] rho E_12 is rational in psi, with poles only at psi = c/j (simple, j <= k+1) and at
 psi = 0, and tends to 0 at oo; it is a power series in psi, so there is no pole at 0, i.e. it lies in span(lam_1..lam_{k+1}).
 The same argument at order t^-1 shows that E_12 is regular (F3).
 E_21 (added 2026-10-31, erratum E3 of leaf_structure_proof.md): E_21(b, a) = E_12(a, b) by the orthogonality of Step 6, which needs the QDE
 at the regular idempotent (leaf_structure_proof.md Theorem 8.4, Corollary 8.5, Proposition 9.4); already the regularity of E_21 needs the
 t^-1 part of it (Proposition 8.6(b)). Step 6's argument for the orthogonality is made rigorous by Corollary 8.5.
Step 10 (regular idempotent). U2 is a nondegenerate critical point away from the log points at t = 0, so r_{a,2}, Delta_2, L_2 are
 analytic at t = 0; Delta_2(0) = x^3/c != 0; Dil_2, E_22 regular (F1, F3).

## 4. Bonus closed forms (to be derived and checked)

- Leaf coefficients as exact functions of t: the coefficient of 1/(c + j z) in z G is j Res_{z=-c/j}(z G), and at z = -c/j the factor
  Exp[-dphi/z] is Exp[j dphi(tau)/c], an honest function of t; only the 2F1 terms with n >= j-1 have a pole there.
- Top coefficients (j = k+1 at order t^k) and the E_11 corner (leaf_structure.md, Step 3 tables) as checks.
- Theta = -x^3 tau s*^2 = -4 x^3 tau/(c + delta)^2.

## 5. What is NOT allowed / practical rules

- No topological recursion, no spectral curves.
- Do not run wolframscript (the Mathematica licence allows only ~2 kernels and the main session uses them). Perl with Math::BigRat is
  available for numeric spot checks at rational x (python is NOT available). List any Mathematica check you want in your output;
  it will be run afterwards.
- Never kill processes. Write files with the Write tool (heredocs with some characters fail in this Git Bash).
