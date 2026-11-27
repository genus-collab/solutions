# Referee 1 report on proof_M_hurwitz.md (logical structure and cited results)

Started 2026-10-31. Notes appended as the review proceeds.

## 0. Scope

Read: proof_M_hurwitz.md (complete), mixed_scale_findings.md, hvc_multiprobe_out.txt, hvc_multiprobe_g4_out.txt,
hvc_which_class_out.txt, hvc_multiprobe.wls, hvc_which_class.wls. No wolframscript used.

Claims under review (from proof_M_hurwitz.md section 7):
 (P1) kappa / Lambert / Hurwitz reformulations, s-decomposition Xi = Sum_k (-s)^(-k) Xi^(k).
 (P2) dictionary: base point q_0 = lam_0 - lam_1 <=> TR on C_tau: x = w/(1+w)^(tau+1), y = -1/w (derived, not proved).
 (P3) Theorem 3.1: tau = 1, all g >= 2, all n: TR correlators of C_1 have x-degree >= 2g-2+n at the zeros of x (conditional on P2).
 (P4) tau = -1/2 variant for all-odd parts.
 (P5) single class: Xi_g(2) = 0 for all g >= 1 by cut-and-join.

## 1. Theorem 3.1 (tau = 1): line-by-line check of the TR argument

Coordinates: zeta = (v-1)/(v+1) is odd under sigma(v) = 1/v; then x = (1 - zeta^2)/4, y = -(1+zeta)/(1-zeta),
B = dzeta dzeta'/(zeta-zeta')^2, dX = -2 zeta dzeta/(1-zeta^2), y(zeta) - y(-zeta) = -4 zeta/(1-zeta^2).
The zeros of x are zeta = +1 (v = oo, the pole of y) and zeta = -1 (v = 0); v = -1 is zeta = oo (a double pole of x,
a simple pole of dX, not a branch point).

Checked, and correct:
 - K(v0,v) = -dv0/(2 (v0-v)(v0-sigma v) dX(v)) and the expansion in 1/v0 with coefficients h_m(v, 1/v); (*) has the right sign.
 - Anti-invariance. It does NOT need the Hodge dictionary: omega_{g,n} (2g-2+n > 0) has poles only at zeta = 0, and the
   linear loop equation makes omega(zeta) + omega(-zeta) holomorphic there. This sum is then a global holomorphic
   differential on P^1, hence 0. So omega_{g,n} = Sum C[k] prod dzeta_i/zeta_i^(2k_i+2).
 - Poles of f_m: v = 1, 0, oo, v_j, sigma(v_j). There is none at v = -1 for (g,n) != (1,0): Brk is regular there and 1/dX
   vanishes. sigma^* f_m = f_m (h_m and Brk are symmetric and X o sigma = X), so Res_0 = Res_oo and Res_{sigma v_j} = Res_{v_j}.
 - Degrees: B(v,v_j) expanded for |t| < |t_j| has degree 0; B(sigma v, v_j) has degree >= 2; the stable factors have degree
   >= 2g'-2+n' by induction and anti-invariance; deg h_m >= -m; 1/dX = ((1+t)/(1-t)) t/dt has degree 0.
   Res_{v_j} B phi = d_{v_j} phi(v_j), and this preserves degree. Negative powers of t_j occur in the separate residues.
   They cancel in the sum, and the total-degree bookkeeping is unaffected.
 - Base cases g = 0, 1: holomorphy at x = 0 already gives degree >= 1 per variable, so no xi_k structure is needed.
 - n = 0 (omega_{g,1}) works the same way: [t0^(m+1)] = 0 for m < 2g-2.
 - The j_i = 0 entries are handled by the string/dilaton relations (a), which I re-derived. Xi(0,j') = ((D'-2h+2)/c^2) Xi(j').
   At J = 2h-3+r the prefactor is D'-2h+2 = 0, so the full vanishing range follows. Xi(1,j) = (J/c^2) Xi(j) was also re-derived.
 - The rescaling to c = 1 is legitimate: Xi^(c) = c^-(3h-3+2r) Xi^(c=1), uniformly in all terms (checked by degree count).
Verdict on section 3: the TR statement "omega^{C_1}_{g,n} has x-degree >= 2g-2+n at the zeros of x" is PROVED, for all g >= 2
and all n, with no Hodge input. The Hodge consequence depends only on the dictionary (P2).

## 2. Independent exact checks (my own Perl/BigRat code, not the authors' scripts)

Scripts (scratchpad): trC1.pl (C_1 in zeta = (v-1)/(v+1)) and trhalf.pl (C_{-1/2} in zeta = v + 1/2).
They implement the EO recursion directly as a coefficient recursion in the basis e_k = dzeta/zeta^(2k+2). This basis is
legitimate by the anti-invariance noted in section 1. Kernel: K = P(zeta^2) dz0/(zeta (z0^2 - zeta^2) dzeta), with
P = (1-zeta^2)^2/8 for C_1 and P = (zeta^2 - 1/4)/2 for C_{-1/2}. All omega_{g,n} come out symmetric.
My normalisation reproduces the authors' TR values exactly: [x^1]F_{1,1} = -1/16, Xi_2(3) = -1/512, Xi_3(5) = -3/4096.

(a) Degree bound (Theorem 3.1, pure TR statement), tau = 1: exact zeros for all mu with |mu| <= 2g-3+n, and non-zero at
    |mu| = 2g-2+n, for (g,n) = (2,1..4), (3,1..2), (4,1..2). No violation.
(b) Dictionary at tau = 1 (x = -1/2 in the (c,-x,-1) normalisation, c = 1/2). TR/Mathematica is constant on each (g,n):
    (2,1): 1/2; (2,2): 1/256; (3,1): 1/4; (3,2): 1/4096; (4,1): 1/16384.
    Each equals 2^(1-g) * c^(3h-3+2r) = (tau(tau+1))^(1-g) times the c-rescaling, independent of n.
    Genus 4: TR Xi(7..10) = -225/262144, -1823/262144, 669/327680, 93551/327680; the Mathematica values are 16384 times these.
(c) Dictionary at tau = -1/2 (x = 1, c = 2), odd parts: the ratio is constant, (2,2): -512, (3,2): 16384, (4,1): -131072,
    which is (-1/4)^(1-g) * 2^(3h-3+2r). Zeros exactly in the range, e.g. g = 3: (1,1), (1,3) vanish; (1,5), (3,3) do not.
Conclusion of section 2: the dictionary (P2) holds EXACTLY at tau = 1 and tau = -1/2 in every case I could compare (g <= 4, r <= 2).
This is far stronger evidence than the authors' float checks (their genus-3 tau = 1 check has abs. error 6e-7, rel. error 1e-3).
(d) Genus 5, tau = 1, n = 1: Xi_5(j) = 0 exactly for j = 1..8 and Xi_5(9) = -2205/1048576 != 0, as Theorem 3.1 predicts.

## 3. The other "proved" items

 - 1a (kappa form): correct. psi_i = pi^*psi_i + D_{i,n+1} and psi_{n+1} D = 0; the AC pushforward formula; p_k = -c^(-k)/k.
 - 1b (Lambert form): correct; [x^mu] xi_k = mu^(mu+k)/mu!. h = 1 by hand: Xi_1(j) = (j-1)/24 - 1/24 = (j-2)/24 (re-derived).
 - 1c (Hurwitz limit): correct. The analytic continuation the authors worry about is harmless. Along (c_0, c_1) = (s, -s)
   the point -psi + s(lam_0 - lam_1) = -(1+s)psi - s(psi^2 + ...) has no t_0 shift, so F is rational in s.
   X = -s e^s stays on the principal Lambert branch T = -s (X is monotone in T for T < 1). The string-flow signs check out:
   c_0 -> c_0 + a, c_mu -> c_mu e^(-mu a), and the probe lam_j picks up e^(jL).
   One genus-1 caveat: <Omega tau_0>_{1,1} = -1/24 breaks the string equation for F_1 itself. Only derivatives (W_j) and
   0*F_1 enter, so nothing breaks.
 - Section 4 (s-decomposition): correct; (tau-x)(-1-tau-x) = -(s - x - x^2), e_h(x+x^2) Lambda^vee(1) = lambda_h c(E)c(E^vee) = lambda_h.
 - 5, cut-and-join Xi_g(2) = 0: correct. On p = X delta_1 only the (i,j) = (1,1) join survives. The Euler operator is
   2g-2+2X d/dX, and dH/dp_2 = 2 W_2. Orders: (2g-2)H_g and 2X W_1 are O(L^(2-2g)), while X^2 W_2 = L^(3-2g) Xi(2) + O(L^(2-2g)).
   The argument is specific to the single class (ELSV + CJ) and agrees with the data (factor (j-2) at h = 1, 2, 3).

## 4. Is the dictionary (P2) justified? A rigorous route exists

The file's derivation is a sketch. It treats lam_0 and lam_1 separately, and each of those deformations gives a y with a
pole at the branch point: an irregular curve, for which the variational formula is not safe (cf. y = z + eps/z on the Airy
curve, where TR at fixed eps is not the eps-expansion). It also invokes "linearity of the background" (DOSS) without
spelling it out. A clean proof uses only the REGULAR combined deformation Omega = Y_0 dX:
 (1) BKMP (Chen; Zhou): TR(X, Y_B) = A_g B^n (CohFT_{g,n} at q_B = -kappa psi) in the xi_k basis.
 (2) EO variational formula along Y_e = Y_B + e Y_0, regular near e = 0, iterated. The Taylor coefficients in e at e = 0
     equal A_g B^n (d_{q_0})^k CohFT(q_B), because Int Lambda d xi_k = [psi^k](lam_0 - lam_1). These pairings are exact
     computations: Int_0^oo d xi_k = -delta_{k0}/tau and Res_0 x^-1 d xi_k = 1.
 (3) Both sides are rational in e. TR depends polynomially on Taylor coefficients of Y_e at the branch point, which are
     affine in e, divided by powers of Y_e'(a); the CohFT side has a pole only where the dilaton coefficient vanishes.
     So TR(Y_e) = A_g B^n CohFT(q_B + e B q_0) identically.
 (4) Homogeneity in Y and the dilaton equation, with e -> oo, give TR(X, Y_0) = A_g B^(2-2g) CohFT_{g,n}(q_0).
     The factor is independent of n, exactly as observed in section 2(b).
So P2 is a theorem modulo BKMP for C^3 and the EO variational formula, both published. The file should state it this way.

## 5. Minor issues (none breaks the claim)

 1. Section 2: the lam_0 and lam_1 deformations are separately irregular, as the authors themselves note. Use the
    combined, regular deformation and the rationality argument of section 4 instead of "linearity (DOSS)".
 2. Anti-invariance is derived from the xi_k form, which presupposes the dictionary. It follows from TR alone (linear loop
    equation plus a global involution of P^1), so Theorem 3.1 is an unconditional TR theorem. Only its translation into
    Xi uses P2. Saying so would make the logical dependence cleaner.
 3. tau = -1/2: x = w(1+w)^(-1/2) is not single-valued on P^1. The argument only uses dX and a local branch at v = oo;
    the sign convention x o sigma = -x should be stated with that branch. Also "h_m(v,-1-v) ~ v^m" is false for odd m:
    the leading coefficient is Sum_{i+k=m} (-1)^k = 0. Harmless, since only deg <= m is used. Here B(sigma v, v_j) has
    degree 0, not >= 2 as at tau = 1. Also harmless.
 4. The "iff two weights coincide" remark about global involutions is unproved. It is not used.
 5. The file's numerical checks are weaker than stated: genus 3 at tau = 1 has abs. error 6e-7 against values ~7e-4.
    My exact checks (section 2) supersede them.
 6. v = -1 is a double pole of x, hence a ramification point of x as a map. TR is taken with respect to X = log x, whose
    branch points are the zeros of dX only. This convention, standard for BKMP, should be stated explicitly.

## 6. Counterexample test (does the argument prove too much?)

The argument uses two properties that fail for the WHICH-CLASSES counterexamples:
 (i) The spectral-curve dictionary. It requires the CY condition and a leaf scale equal to a weight; only then are the
     probes x^j-coefficients of the framed curve.
 (ii) A global involution exchanging two "good" zeros of x, which needs tau in {1, -2} (or -1/2 for odd parts).
The trivial (KW) class, Lambda(1)Lambda(w2), the non-CY triple and the CY triple with leaf scale 2 have no such curve,
so the argument says nothing about them. It also does not overreach where it does apply: at tau = 1 it gives only
j <= 2h-2 and not the single-class extra zero j = 2h. That is consistent with the data: Xi_2(4; tau=1) = -11/1280 != 0.

## 7. Verdict

NOT REFUTED. What is proved:
 - Theorem 3.1 as a TR statement: on C_1 (x = v/(1+v)^2, y = -v), omega_{g,n} has x-degree >= 2g-2+n at the zeros of x,
   for all g >= 2 and n >= 1. Each step checks out, and it agrees with exact recursion through g = 5.
 - Translation to Conjecture M (vanishing direction, all h, r) at tau = 1 and -2, and at tau = -1/2 for all-odd parts.
   This is conditional on P2. P2 is not proved in the file, but it follows from BKMP plus the EO variational formula by
   the route in section 4, and it holds exactly in every case I compared (g <= 4, r <= 2).
 - The kappa/Lambert/Hurwitz forms, the s-decomposition and Xi_g(2) = 0 (single class) are correct.
Open, as the authors state: the single class beyond j = 1, 2; general tau; the extra zero j = 2h for h >= 2; the
non-vanishing direction.

Scripts (copies): review_checks/referee1_trC1.pl, referee1_trhalf.pl, referee1_cmpM.pl (Perl + Math::BigRat; no Mathematica).
Usage: "perl referee1_trC1.pl GN GMAX SMAX" (e.g. 7 4 5 for genus 4, n = 1), "perl referee1_trhalf.pl 7 4 5", "perl referee1_cmpM.pl".
