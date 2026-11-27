# S1. The singular row from the Euler integral (Steps 1-2 of outline.md, and F2)

Written 2026-10-31 for leaf_proof/outline.md; revised the same day after a referee report (verdict "sound"; the changes are
listed at the end of section 8). No topological recursion and no spectral curves are used. Every statement is marked [formal]
(an identity of formal power series, valid for every x with x != 0, c != 0) or [analytic] (a statement about convergent
integrals, with explicit hypotheses). The only analytic input on which the conclusions depend is the value of one constant
series A(z) (section 6); besides it, only the analyticity in tau of finitely many explicit algebraic functions is used (1.2,
1.3, Lemma 1.1, (D2)). Exact spot checks in Perl/Math::BigRat: leaf_proof/S1_checks.pl (checks C1-C6, 62 checks, all pass,
50 s; output leaf_proof/S1_checks_out.txt). The Mathematica facts V1-V8 of outline.md are quoted only as confirmations, never
used in a proof (the labelling remark in section 0 and the notational remark 7.6, which no proof in S1 uses, quote V1, which
S2 proves).

## 0. Conventions and the result

Conventions (rc_lib.m s.1, s.3; hvc_lib.m s.5; alpha0 = 0, alpha1 = eps1 = x, alpha2 = eps2 = 1):

- x is a parameter with x != 0 and c := 1 + x != 0 (when x is symbolic, identities are identities in Q(x)).
- Q = e^t, tau := 1 - 1/Q = 1 - e^-t; so Q = 1/(1 - tau), Q tau = Q - 1, dtau/dt = 1 - tau, tau = t - t^2/2 + O(t^3).
- W(U) := spotU = x Log[1 - Q U] - Log[U] + Log[1 - U]; theta := U d/dU; W_k := theta^k W (rational in U, Q for k >= 1).
- Critical points. theta W = -x Q U/(1 - Q U) - 1 - U/(1 - U); multiplying by (1 - Q U)(1 - U) gives
  -x Q U (1 - U) - (1 - Q U)(1 - U) - U (1 - Q U) = LeqU = Q x U^2 + Q (1 - x) U - 1. Its roots satisfy
  U1 + U2 = SsumQ = (x - 1)/x, U1 U2 = PprodQ = -1/(Q x), and the discriminant is ((x-1)^2 + 4x/Q)/x^2 = (c^2 - 4 x tau)/x^2.
  With delta := Sqrt[c^2 - 4 x tau], branch delta(0) = c,
      U1 = (x - 1 + delta)/(2x),   U2 = (x - 1 - delta)/(2x),   U1(0) = 1,   U2(0) = -1/x.
  Idempotent 1 is U1. The scripts take sq = ((1 + x)/x) Sqrt[1 + 4 x (1/Q - 1)/(1 + x)^2] = (c/x) (1 - 4 x tau/c^2)^(1/2),
  expanded in t, so sq = delta/x with delta(0) = c whatever the sign of c, and U1 = (S0 + sq)/2 (hvc_leaf_proof_check.wls,
  hvc_leaf_structure.wls). Only the consistency of this labelling matters, and it holds: the code computes every graph-sum
  factor (rows, Delta, edge and dilaton factors) in the quadratic algebra Q(x, Q)[U]/(LeqU) (elements C[a, b] = a + b U; the
  data of idempotent 2 are stored in conjugate encoding, via cconj), and Lemma C/L and all checks (hvc_leaf_proof_check.wls,
  leaf_structure.md) evaluate every factor at the same root U1, which is the root with the singular Delta (F1: Delta_1 =
  x^2 W_2(U1) by V1, = -x delta/(tau^2 s*^2) by (1.5)). Evaluation U -> U1 is a ring homomorphism to algebraic functions of Q,
  so cEval/cinv followed by U -> U1 is substitution of the values W_k(U1) (cconj followed by U -> U1 gives the values at U2).
- Point-frame identity row of idempotent 1. rmatrixAt0 calls rmatrixCore with ffc1 = {1, 0, 0, ...}, so
  ratCN[0,1][m] = wickCN[m] = wickZ[[m+1]] with Wk[k] -> W_k reduced in the quadratic algebra; evaluated at U = U1:
      r^pt_{0,1}(z; t) := Sum_m z^m wickZ[[m+1]] /. {Wk[k] -> W_k(U1(t))},
      wickZ[[m+1]] = Sum_{n=m}^{3m} (2n)! (-1)^n/(n! 2^n W_2^n) Sum_{p} prod_k W_k^(m_k)/((k!)^(m_k) m_k!),
  p running over the partitions of 2n into exactly n - m parts >= 3 (m_k = multiplicity of k in p). This is the
  coefficient of z^m of
      Sum_{n>=0} ((2n)!/(n! 2^n)) (-z/W_2)^n [v^(2n)] Exp[Sum_{k>=3} W_k v^k/(k! z)]
  (a term with j factors W_k and v-degree 2n has z-order n - j). The cut-off "parts <= 2 zmax + 2" in wickDirect never
  binds: at order z^m with j = n - m >= 1 factors the largest part is 2n - 3(j - 1) = 3m + 3 - n <= 2m + 2. Every monomial
  W_2^-n prod_i W_{k_i} of wickZ[[m+1]] has Sum_i k_i = 2n, i.e. weight 0 if W_k is given weight k (Lemma 2.6).
- Frame B. buildFactors[zmax, "B", -1, x] applies twistRows[hCoeffs[{c, -x, -1}, -1, zmax]], i.e. multiplies every row by
  Exp[-ex], ex := -Sum_{k>=1} B_2k/(2k(2k-1)) z^(2k-1) Sum_a w_a^(1-2k), w = (c, -x, -1). Since Sum_a w_a^(1-2k) =
  c^(1-2k) - x^(1-2k) - 1,
      h_w(z) := Exp[ex] = Exp[Sum_{k>=1} B_2k/(2k(2k-1)) z^(2k-1) (x^(1-2k) + 1 - c^(1-2k))]  in 1 + z Q(x)[[z]],
  and the frame-B row is r_{0,1}(z; t) = h_w(z)^-1 r^pt_{0,1}(z; t); h_w does not depend on t.
- F(z, tau) := Sum_{n>=0} F_n(z) tau^n, F_n := (alpha)_n (beta)_n/((gamma)_n n!), alpha = beta = 1/z + 1, gamma = c/z + 2.
  Since (1/z + 1)_n = z^-n prod_{j=1}^n (1 + j z) and (c/z + 2)_n = z^-n prod_{j=2}^{n+1} (c + j z),
      F_n(z) = z^-n prod_{j=1}^n (1 + j z)^2 / (n! prod_{j=2}^{n+1} (c + j z))  in Q(x)(z).
- phi_tau(s) := Log[s] + x Log[1 - s] - Log[1 - tau s];  s*(tau) := 2/(c + delta);  dphi(tau) := phi_tau(s*(tau)) - phi_0(1/c);
  alg(tau) := (1 - tau s*) (phi''_tau(s*)/phi''_0(1/c))^(1/2), the branch with alg(0) = 1. (' = d/ds on phi.)

Rings. C{tau}: convergent power series; O(D) (O(D')): functions analytic on the disc D (D') of 1.2. C((z)): formal Laurent
series. C(z) is embedded in C((z)) by Laurent expansion at z = 0 (an injective field homomorphism, used tacitly for F_n and
c/(c+z)). The map iota: C{tau}[[z]] -> C[[z]][[tau]] = C[[z, tau]], "Taylor-expand every z-coefficient", is an injective
ring homomorphism commuting with d/dtau and with multiplication by elements of C{tau}. Exp[+-dphi/z] := Sum_j (+-dphi/z)^j/j!
in C[1/z][[tau]] (tau-adically convergent since dphi = O(tau)); its tau^k coefficient is a polynomial in 1/z of degree <= k.
C((z))[[tau]] contains all of these.

**Theorem S1.** Let x != 0, c = 1 + x != 0.
 (a) [Step 1] Under U = 1 - tau s, e^{W/z} dU/U = kappa_ph(z) Q^(x/z) tau^(c/z+1) e^{phi_tau(s)/z} (1 - tau s)^-1 ds with a
     phase constant kappa_ph independent of s and t (it is not the constant kappa(z) of outline Step 5; see 7.6). s* is the
     critical point of phi_tau analytic at tau = 0 (s*(0) = 1/c, s* = (1 - U1)/tau); it is nondegenerate,
     phi''_tau(s*) = -delta/(x s*^2 (1 - tau s*)^2), phi''_0(1/c) = -c^3/x, and W_2(U1) = U1^2 phi''_tau(s*)/tau^2 =
     -delta/(x tau^2 s*^2). For every k >= 2, w_k := tau^k W_k(U1) is analytic at tau = 0 with w_k(0) = (-1)^k phi_0^(k)(1/c)
     (Lemma 1.1).
 (b) [Step 2(a)] For small t != 0, r^pt_{0,1}(z; t) = N(z, tau), the normalised formal Laplace expansion (section 2) of
     J(z, tau) = Int s^(1/z) (1-s)^(x/z) (1 - tau s)^(-1/z-1) ds at s*(tau). The z-coefficients of r^pt_{0,1} are analytic at
     t = 0 (by the homogeneity of the Wick formula they are polynomials in 1/w_2, w_3, w_4, ...), the identity extends to
     t = 0, where both sides equal A(z) := Lnorm[phi_0, 1; 1/c], and N lies in C{tau}[[z]]; so the t-Taylor series of the
     z-coefficients of r^pt_{0,1} are those of N(z, tau(t)).
 (c) [Step 2(b),(c)] In C((z))[[tau]]:  N(z, tau) = A(z) F(z, tau) Exp[-dphi(tau)/z] alg(tau),  with  A(z) = h_w(z) c/(c + z).
 (d) [Step 2(d)] For every k >= 0, [tau^k] (F Exp[-dphi/z] alg) and [tau^k] (F Exp[-dphi/z]) are rational functions of z
     without pole at z = 0, and in C[[z]][[t]] (and coefficientwise in z, as analytic functions of t near 0)
         r_{0,1}(z; t) = (c/(c + z)) F(z, tau) Exp[-dphi(tau)/z] alg(tau),   tau = 1 - e^-t.
     In particular r_{0,1}(z; 0) = c/(c + z) (F2) and r^pt_{0,1}(z; 0) = h_w(z) c/(c + z).
Status: (a) elementary; (b), (c) apart from the value of A, and (d) are [formal] ((b) uses in addition the analyticity in tau
of finitely many explicit algebraic functions); A = h_w c/(c+z) is [analytic] for real x > 0 (Laplace's method and Stirling's
series as z -> 0+) and extends to all x because both sides have coefficients in Q(x). Section 7.4 gives a second, purely
formal proof of the regularity in (d) by eikonal and transport equations.

## 1. The rescaled thimble (Step 1)

### 1.1 Change of variables
For tau != 0 put U = 1 - tau s. Then, using Q tau = Q - 1,
    1 - U = tau s,    1 - Q U = 1 - Q + Q tau s = (1 - Q)(1 - s) = -Q tau (1 - s),    dU/U = -tau ds/(1 - tau s).
Hence, on any simply connected domain of the s-plane avoiding s = 0, 1, 1/tau, with compatible branches of the logarithms,
    W(1 - tau s) = x Log[-Q tau (1 - s)] + Log[tau s] - Log[1 - tau s] = phi_tau(s) + c Log[tau] + x Log[Q] + x Log[-1],   (1.1)
    e^{W/z} dU/U = kappa_ph(z) Q^(x/z) tau^(c/z + 1) e^{phi_tau(s)/z} (1 - tau s)^-1 ds,    kappa_ph(z) := -e^{x Log[-1]/z},   (1.2)
and kappa_ph depends on z and on the branches only, not on s or t (for real x, z it is a phase up to sign). It is a constant
of the integrand only: it is not the t-independent constant kappa(z) of outline Step 5 (s_0 = kappa Q^{x/z} tau^{c/z+1} F),
which also absorbs e^{phi_0(1/c)/z}, the metric normalisation Delta_1^(-1/2) and the constant c/(c + z) of the frame-B row
(h_w c/(c + z) in the point frame); the relation is (7.4) in 7.6. (S2 s.5.3 refers to kappa_ph as "the phase constant kappa
of Theorem S1(a)".) In Pochhammer form
    e^{phi_tau(s)/z} (1 - tau s)^-1 = s^(1/z) (1 - s)^(x/z) (1 - tau s)^(-1/z - 1).
The logarithmic points U = 1, U = 1/Q, U = 0 of W are s = 0, s = 1, s = 1/tau. At t = 0 the first two collide (1 = 1/Q)
and U1 runs into them; in the coordinate s they stay at 0 and 1, the critical point stays between them (s*(0) = 1/c), and
U = 0 escapes to s = oo.

### 1.2 The critical point s*
    phi'_tau(s) = 1/s - x/(1 - s) + tau/(1 - tau s) = (x tau s^2 - c s + 1)/(s (1 - s)(1 - tau s)),
since (1 - s)(1 - tau s) - x s (1 - tau s) + tau s (1 - s) = 1 - c s + x tau s^2. So for tau != 0 the critical points of
phi_tau in C - {0, 1, 1/tau} are the roots of  x tau s^2 - c s + 1 = 0:
    s* = (c - delta)/(2 x tau) = 2/(c + delta),     s** = (c + delta)/(2 x tau) = 1/(x tau s*)
((c - delta)(c + delta) = 4 x tau). Under U = 1 - tau s: 1 - tau s* = 1 - (c - delta)/(2x) = (x - 1 + delta)/(2x) = U1 and
1 - tau s** = U2. So s* = (1 - U1)/tau (V2), and s** (pole at tau = 0) is the regular idempotent.
Analyticity. delta = c (1 - 4 x tau/c^2)^(1/2) (principal root) is analytic on D := {|tau| < |c|^2/(4|x|)} and delta/c has
positive real part there, so c + delta != 0 and s* is analytic on D, s*(0) = 1/c; equivalently, the implicit function theorem
applies since d/ds (x tau s^2 - c s + 1) = -c != 0 at (tau, s) = (0, 1/c). Series: s* = (1/c) Sum_n Cat_n (x tau/c^2)^n.
If s* = 1 or s* = 1/tau, the critical equation forces tau = 1; so on the disc D' := D intersected with {|tau| < 1} we have
s* not in {0, 1, 1/tau}. For real x > 0, (1 + x)^2/(4x) >= 1 and D' = {|tau| < 1}.
Identities at s = s* (consequences of the critical equation; used throughout):
 (a) c - 2 x tau s* = delta              [2 x tau s* = 4 x tau/(c + delta) = c - delta];
 (b) 1 - s* = x s* (1 - tau s*)          [this is 1 - c s* + x tau s*^2 = 0];
 (c) s*/(1 - tau s*) = (1 - x s*)/(1 - tau)   [s*(1 - tau) - (1 - x s*)(1 - tau s*) = -(x tau s*^2 - c s* + 1) = 0];
 (d) ds*/dtau = x s*^2/delta              [differentiate: x s*^2 + (2 x tau s* - c) s*' = 0, then (a)].

### 1.3 Nondegeneracy; W_2(U1) versus phi''_tau(s*)
phi''_tau(s) = -1/s^2 - x/(1 - s)^2 + tau^2/(1 - tau s)^2. By (b), x/(1 - s*)^2 = 1/(x s*^2 (1 - tau s*)^2), hence
    (1 - tau s*)^2 phi''_tau(s*) = [-(1 - tau s*)^2 - 1/x + tau^2 s*^2]/s*^2 = (2 tau s* - c/x)/s*^2 = -delta/(x s*^2)      (1.3)
by (a). So phi''_tau(s*) = -delta/(x s*^2 (1 - tau s*)^2) != 0 on D' (delta^2 = c^2 - 4 x tau != 0 on D), and at tau = 0
phi''_0(1/c) = -c^3/x. Consequently
    alg^2 = (1 - tau s*)^2 phi''_tau(s*)/phi''_0(1/c) = delta/(c^3 s*^2),
    alg = (delta/c)^(1/2)/(c s*) = (delta/c)^(1/2) (c + delta)/(2c);                                              (1.4)
the right side is analytic on D with value 1 at 0, so it is the branch of the definition.
W_2 in the coordinate y = Log[U] (theta = d/dy). Let h(s) := Log[1 - tau s], so h'(s) = -tau/(1 - tau s). For a function G
of U and the composite G(1 - tau s) = (G o exp o h)(s):
    d/ds G(1 - tau s) = -(tau/(1 - tau s)) (theta G)(1 - tau s),
    d^2/ds^2 G(1 - tau s) = (tau/(1 - tau s))^2 [(theta^2 G) - (theta G)](1 - tau s)
(differentiate the first line once more: d/ds(-tau/(1 - tau s)) = -tau^2/(1 - tau s)^2). Take G = W: by (1.1),
phi_tau(s) = W(1 - tau s) + const, so critical points correspond, and at s = s*, where W_1(U1) = 0,
    W_2(U1) = (1 - tau s*)^2 phi''_tau(s*)/tau^2 = U1^2 phi''_tau(s*)/tau^2 = -delta/(x tau^2 s*^2).      (1.5)
W_2(U1) has a double pole at tau = 0: the coordinate y degenerates there (h'(s*) = O(tau)), the coordinate s does not. (More
generally W_k(U1) has a pole of order k, compensated inside every single Wick monomial: Lemma 1.1, Remark 3.2.)
[C1 checks (1.5) exactly at (x, tau) = (2, 5/8) and (1/3, 7/12).] With V1 (Delta_1 = x^2 W_2(U1), Step 3) this is
Delta_1 = -x delta/(tau^2 s*^2) = -x c^3 alg^2/tau^2 (V2); in S1 this is used only in the remarks of section 0 (labelling) and
7.6 (notation), not in any proof.

### 1.4 Euler's integral and Gauss's equation [analytic; motivation only, not used in any proof below]
Let J(z, tau) := Int_0^1 s^(1/z) (1 - s)^(x/z) (1 - tau s)^(-1/z-1) ds (principal powers), under the hypotheses
    (H_E)  Re(1/z) > -1,  Re(x/z) > -1,  |tau| < 1        (e.g. z > 0, x > 0, |tau| < 1).
Euler's formula: J(z, tau) = B(1/z + 1, x/z + 1) F(z, tau).  (1.6)
Proof: put a := 1/z + 1, b := 1/z + 1, g := c/z + 2, so g - b - 1 = x/z. For |tau| < 1, (1 - tau s)^-a = Sum_n (a)_n (tau s)^n/n!
uniformly on [0, 1]; integrating termwise, Int_0^1 s^(b-1+n) (1-s)^(g-b-1) ds = B(b + n, g - b) = B(b, g - b) (b)_n/(g)_n, so
J = B(b, g - b) Sum_n (a)_n (b)_n/((g)_n n!) tau^n, and B(b, g - b) = B(1/z + 1, x/z + 1).
Gauss's equation by integration by parts. For general a, b, g let K := s^(b-1) (1 - s)^(g-b-1) (1 - tau s)^-a. Then
d_tau K = a s (1 - tau s)^-1 K, d_tau^2 K = a (a+1) s^2 (1 - tau s)^-2 K, and the exact identity
    [tau (1 - tau) d_tau^2 + (g - (a + b + 1) tau) d_tau - a b] K = -a d_s [s^b (1 - s)^(g-b) (1 - tau s)^(-a-1)]      (1.7)
holds: after division by a K (1 - tau s)^-2 both sides equal -b + g s + (b - a - 1) tau s + (a + 1 - g) tau s^2 (direct
expansion; the tau^2 s^2 terms cancel on the left). For a = b = 1/z + 1, g = c/z + 2 this is
    Gauss_z [s^(1/z) (1-s)^(x/z) (1 - tau s)^(-1/z-1)] = -(1/z + 1) d_s [s^(1/z+1) (1-s)^(x/z+1) (1 - tau s)^(-1/z-2)],   (1.8)
    Gauss_z := tau (1 - tau) d_tau^2 + (c/z + 2 - (2/z + 3) tau) d_tau - (1/z + 1)^2.
Under (H_E) and |tau| <= r < 1, |d_tau^j K| <= C_r s^(Re(1/z)) (1 - s)^(Re(x/z)), which is integrable, so J is analytic in tau
and can be differentiated under the integral sign; the bracket on the right of (1.8) is O(s^(Re(1/z)+1)) at s = 0 and
O((1 - s)^(Re(x/z)+1)) at s = 1, so the boundary terms vanish and Gauss_z J = 0. This also re-proves (1.6) by the recursion of
Lemma 5.1, now at a numerical z: the Taylor coefficients J_k of J at tau = 0 satisfy (k + 1)(k + gamma) J_{k+1} =
(k + alpha)(k + beta) J_k, and k + gamma != 0 for every k >= 0 because Re gamma = Re(1/z) + Re(x/z) + 2 > 0 under (H_E); so
J = J(z, 0) F with J(z, 0) = B(1/z + 1, x/z + 1).
Below, the formal expansion of J obeys Gauss's equation because of the pointwise identity (1.8) alone (section 4):
the formal functional kills total derivatives exactly, so no boundary terms and no hypothesis (H_E) are needed.

### 1.5 The rescaled derivatives tau^k W_k(U1)
Derivatives of the phase. For k >= 1,
    phi_tau^(k)(s) = (k - 1)! [(-1)^(k-1) s^-k - x (1 - s)^-k + tau^k (1 - tau s)^-k],                                (1.9)
since d^k/ds^k Log[s] = (-1)^(k-1) (k-1)! s^-k, d^k/ds^k Log[1 - s] = -(k-1)! (1 - s)^-k and d^k/ds^k Log[1 - tau s] =
-(k-1)! tau^k (1 - tau s)^-k. In particular phi_0^(k)(1/c) = (k - 1)! c^k ((-1)^(k-1) - x^(1-k)), as (1 - 1/c)^-k = c^k x^-k.
Derivatives of W. For a in {1, Q} put p_a := (1 - a U)^-1. Then theta Log[1 - a U] = -a U/(1 - a U) = 1 - p_a and
theta p_a = a U p_a^2 = p_a^2 - p_a (as a U = 1 - 1/p_a). Define polynomials P_0(p) := p, P_n(p) := (p^2 - p) P_{n-1}'(p)
(P_1 = p^2 - p, P_2 = 2p^3 - 3p^2 + p). Then theta^k Log[1 - a U] = -P_{k-1}(p_a) for k >= 2, by induction on k
(theta (1 - p_a) = -(p_a^2 - p_a) = -P_1(p_a), and theta(-P_{n-1}(p_a)) = -P_{n-1}'(p_a) (p_a^2 - p_a) = -P_n(p_a)), while
theta^k Log[U] = 0 for k >= 2 (theta Log[U] = 1). Hence
    W_k = -x P_{k-1}(p_Q) - P_{k-1}(p_1)     (k >= 2).                                                              (1.10)
By induction P_n has degree n + 1 with leading coefficient n!, and P_n(0) = 0 (P_0 = p, and P_n has the factor p^2 - p for
n >= 1); write P_{k-1}(p) = Sum_{i=1}^{k} a_{k-1,i} p^i, so a_{k-1,k} = (k - 1)!.

**Lemma 1.1 (rescaled derivatives).** For k >= 2 and tau in D' - {0} put w_k(tau) := tau^k W_k(U1). Then
    w_k = -Sum_{i=1}^{k} a_{k-1,i} tau^(k-i) [ x (-(1 - tau)/(1 - s*))^i + s*^(-i) ].                                 (1.11)
The right side is analytic on D', tau = 0 included, and
    w_k(0) = -(k - 1)! c^k (1 + (-1)^k x^(1-k)) = (-1)^k phi_0^(k)(1/c),      w_2 = -delta/(x s*^2)  (no zero on D').
Proof. For tau in D' - {0}, U1 = 1 - tau s* is not in {1, 1/Q} (tau s* != 0, and Q U1 = 1 would force s* = 1, which 1.2
excludes on D'), so W_k(U1) is defined (by (1.10), W_k has poles only at U = 1 and U = 1/Q). At U1, 1 - U1 = tau s*, and
1 - Q U1 = -Q tau (1 - s*) by 1.1 (with s = s*); with 1/Q = 1 - tau,
    tau p_1(U1) = 1/s*,     tau p_Q(U1) = -(1 - tau)/(1 - s*).
Multiplying (1.10) by tau^k and using tau^k P_{k-1}(p) = Sum_{i=1}^k a_{k-1,i} tau^(k-i) (tau p)^i gives (1.11). Only the powers
tau^(k-i) with 0 <= k - i <= k - 1 occur, and s*, 1/s*, 1/(1 - s*) are analytic on D' (1.2), so the right side of (1.11) is
analytic on D'. At tau = 0 only i = k survives; with s*(0) = 1/c, 1 - s*(0) = x/c and a_{k-1,k} = (k - 1)!,
    w_k(0) = -(k - 1)! [x (-c/x)^k + c^k] = -(k - 1)! c^k (1 + (-1)^k x^(1-k)),
and by (1.9), (-1)^k phi_0^(k)(1/c) = (-1)^k (k - 1)! c^k ((-1)^(k-1) - x^(1-k)) = -(k - 1)! c^k (1 + (-1)^k x^(1-k)). Finally
w_2 = tau^2 W_2(U1) = -delta/(x s*^2) by (1.5), and delta^2 = c^2 - 4 x tau has no zero on D. QED
So W_k(U1) = tau^-k w_k = O(tau^-k) and W_2(U1)^-1 = O(tau^2): the coordinate y degenerates at t = 0 exactly like -tau times a
good coordinate ((3.4) in Remark 3.2). [C6: (1.11) against tau^k W_k(U1) computed from (1.10) at (x, tau) = (2, 5/8) for
k <= 12 (w_2 = -25/4, W_2(U1) = -16 there); w_k(0) for k <= 12 and w_2 = -delta/(x s*^2) as tau-series at x = 2, 1/3, -5/2.]

## 2. Formal Laplace expansions (Step 2: foundations) [formal, except Lemma 2.5]

### 2.1 The formal Gaussian functional
Let A be the set of formal sums K = Sum_{p >= 0, q in Z} K_pq v^p z^q (coefficients complex numbers, or analytic functions
of parameters) such that, with the weight wt(p, q) := p/2 + q,
 (A1) wt is bounded below on the support of K, and (A2) for every N only finitely many (p, q) with wt <= N occur.
A is a ring (weights add under multiplication; by (A1), (A2) every coefficient of a product is a finite sum, and the product
again satisfies (A1), (A2)). It contains C[[v]][z, 1/z] (a monomial v^p z^q with q >= -J has weight >= -J, and weight <= N
forces p <= 2(N + J), q <= N), and Exp[R(v)/z] := Sum_j R^j/(j! z^j) for every R in v^3 C[[v]]: (R/z)^j has weight >= j/2,
so only j <= 2N contribute to weights <= N. d/dv is a derivation of A (weights drop by 1/2), multiplication by v raises
weights by 1/2, and d/dv Exp[R/z] = (R'/z) Exp[R/z].
For a != 0 define the Gaussian functional (the formal moments of a centred Gaussian of variance -z/a)
    <<K>>_a := Sum_{p even, q} K_pq (p - 1)!! (-z/a)^(p/2) z^q   in C((z))      ((-1)!! := 1).
Each power z^w receives contributions only from monomials of weight w, finitely many by (A2), and w is bounded below by (A1);
<<.>>_a is linear over C((z)) and over functions of parameters, and maps elements of weight >= N into z^N C[[z]].
**Lemma 2.1.** For K in A:
 (G0) scaling: for lambda != 0, <<K(lambda v)>>_{lambda^2 a} = <<K(v)>>_a;
 (G1) Gaussian integration by parts: <<dK/dv>>_a = -(a/z) <<v K>>_a;
 (G2) variance derivative: if a = a(tau) and K depends analytically on tau (coefficientwise, with support independent of tau),
      d/dtau <<K>>_{a(tau)} = <<dK/dtau>>_a + (a'/(2z)) <<v^2 K>>_a + (a'/(2a)) <<K>>_a.
Proof. Each side is <<L K>> for an operator L on A that shifts weights by a bounded amount, so every z-coefficient of either
side involves finitely many monomials of K and it suffices to take K = v^p (and, for (G2), to use the product rule for the
tau-dependence of the coefficients). Put m_p := <<v^p>>_a, so m_{2n} = (2n-1)!! (-z/a)^n, m_odd = 0.
(G0): lambda^(2n) (2n-1)!! (-z/(lambda^2 a))^n = m_{2n}. (G1): for p = 2n + 1, <<p v^(2n)>> = (2n+1)(2n-1)!! (-z/a)^n =
(2n+1)!! (-z/a)^n = -(a/z) (2n+1)!! (-z/a)^(n+1) = -(a/z) m_{2n+2}; for p even both sides vanish. (G2): d/dtau (-z/a)^n =
-n (a'/a) (-z/a)^n, and (a'/(2z)) m_{2n+2} + (a'/(2a)) m_{2n} = m_{2n} (a'/(2a)) [-(2n+1) + 1] = -n (a'/a) m_{2n}; odd p: 0 = 0. QED

### 2.2 Definition
Data: a point u0, a phase f analytic near u0 (only its Taylor series matters) with f'(u0) = 0 and f2 := f''(u0) != 0
(nondegenerate critical point), and an amplitude g in O_{u0}[z, 1/z] (a Laurent polynomial in z with coefficients analytic
near u0). Write f(u0 + v) = f(u0) + f2 v^2/2 + R_f(v), R_f(v) := Sum_{k>=3} f^(k)(u0) v^k/k!. Define
    Lnorm[f, g; u0] := << g(u0 + v) Exp[R_f(v)/z] >>_{f2}   in C((z))     (Gaussian-normalised formal Laplace expansion),
    Lhat[f, g; u0]  := Lnorm[f, g; u0]/g(u0)    for g in O_{u0} with g(u0) != 0   (normalised: leading coefficient 1),
    Lfull[f, g; u0] := e^{f(u0)/z} (2 pi z)^(1/2) (-f2)^(-1/2) Lnorm[f, g; u0]     (unnormalised; a formal symbol, with a
                       chosen root (-f2)^(-1/2)).
Lfull is the formal counterpart of Int e^{f/z} g du over a path through u0 (Lemma 2.5). Immediate properties:
 (D1) Wick form: Lnorm[f, g; u0] = Sum_{n>=0} ((2n)!/(n! 2^n)) (-z/f2)^n [v^(2n)] (g(u0 + v) Exp[R_f(v)/z]). For g = 1 the
      coefficient of z^m is wickZ[[m+1]] with Wk[k] -> f^(k)(u0) (compare section 0).
 (D2) For g in O_{u0}, Lnorm[f, g; u0] lies in C[[z]] with constant term g(u0), and its z^m coefficient is a polynomial with
      rational coefficients in 1/f2, f^(3)(u0), ..., f^(2m+2)(u0), g(u0), ..., g^(2m)(u0). [Both factors have weight >= 0;
      the weight-m part involves g-monomials v^a with a/2 <= m and R-monomials f^(k) v^k/(k! z) of weight k/2 - 1 <= m.]
 (D3) Lnorm depends only on the Taylor coefficients of f and g at u0 (so a translation u = u0 + v changes nothing), is
      unchanged when a constant is added to f, and is linear in g over C[z, 1/z] and over functions independent of u.

### 2.3 Total derivatives
**Proposition 2.2 (P3).** For every k in O_{u0}[z, 1/z]:  Lnorm[f, k' + f' k/z; u0] = 0.
(The amplitude k' + f' k/z is the one of the total derivative: e^{f/z} (k' + f' k/z) = d/du (e^{f/z} k).)
Proof. Put K(v) := k(u0 + v) Exp[R_f(v)/z] in A. Since f'(u0 + v) = f2 v + R_f'(v),
    (k' + f' k/z)(u0 + v) Exp[R_f(v)/z] = dK/dv + (f2/z) v K,
and <<dK/dv>>_{f2} + (f2/z) <<v K>>_{f2} = 0 by (G1). QED

### 2.4 Derivative in a parameter, with a moving expansion point
**Proposition 2.3 (P2).** Let f(tau, s) and g(tau, s) be analytic near (tau0, s0) (g may be a Laurent polynomial in z with
such coefficients), let s*(tau) be analytic with s*(tau0) = s0, d_s f(tau, s*(tau)) = 0 and f2(tau) := d_s^2 f(tau, s*(tau)) != 0.
Put f0(tau) := f(tau, s*(tau)) and N(tau) := Lnorm[f(tau, .), g(tau, .); s*(tau)]. Then every z-coefficient of N is
analytic in tau, and
    d_tau N + (f0'/z - f2'/(2 f2)) N = Lnorm[f(tau, .), d_tau g + g d_tau f/z; s*(tau)],      (2.1)
where f0' = (d_tau f)(tau, s*(tau)) (by the critical-point condition) and f2' is the total derivative of f2. Equivalently,
choosing (-f2(tau))^(-1/2) continuously and using d_tau e^{f0/z} = (f0'/z) e^{f0/z},
    d_tau Lfull[f, g; s*(tau)] = Lfull[f, d_tau g + g d_tau f/z; s*(tau)]:
the expansion of the integrand e^{f/z} g commutes with d/dtau although the expansion point moves.
Proof. Analyticity follows from (D2), as tau -> d_s^k f(tau, s*(tau)), d_s^k g(tau, s*(tau)) and 1/f2(tau) are analytic.
Step 1 (shifted variable). Let ft(tau, v) := f(tau, s*(tau) + v), gt(tau, v) := g(tau, s*(tau) + v). By (D3)
N(tau) = Lnorm[ft(tau, .), gt(tau, .); 0], and d_v ft(tau, 0) = 0 for all tau. Write ft = f0(tau) + f2(tau) v^2/2 + Rt(tau, v)
and Kt := gt Exp[Rt/z] in A. By (G2), with d_tau Kt = (d_tau gt + gt d_tau Rt/z) Exp[Rt/z],
    d_tau N = <<(d_tau gt + gt d_tau Rt/z) Exp[Rt/z]>> + (f2'/(2z)) <<v^2 Kt>> + (f2'/(2 f2)) N.
Since d_tau ft = f0' + f2' v^2/2 + d_tau Rt,
    Lnorm[ft, d_tau gt + gt d_tau ft/z; 0] = <<(d_tau gt + gt d_tau Rt/z) Exp[Rt/z]>> + (f0'/z) N + (f2'/(2z)) <<v^2 Kt>>
                                          = d_tau N + (f0'/z - f2'/(2 f2)) N.
Step 2 (back to s). d_tau gt = (d_tau g + s*' d_s g)(tau, s* + v) and d_tau ft = (d_tau f + s*' d_s f)(tau, s* + v), so
    d_tau gt + gt d_tau ft/z = [(d_tau g + g d_tau f/z) + s*'(tau) (d_s g + (d_s f) g/z)](tau, s* + v).
The last term is s*'(tau) times an amplitude of the form k' + f' k/z (with k = g(tau, .)), whose Lnorm vanishes by (P3).
Translating back by (D3) gives (2.1). The unnormalised form follows from d_tau [e^{f0/z} (-f2)^(-1/2) N] =
e^{f0/z} (-f2)^(-1/2) [d_tau N + (f0'/z - f2'/(2 f2)) N]. QED
In (2.1), f0'/z comes from the moving critical value e^{f0/z}, -f2'/(2 f2) from the Gaussian normalisation (-f2)^(-1/2),
and the motion of the expansion point is a total derivative.

### 2.5 Coordinate invariance
**Proposition 2.4 (P1).** Let h be analytic near v0 with h(v0) = u0 and h1 := h'(v0) != 0, and f, g as in 2.2 with g in O_{u0}.
Then f o h has a nondegenerate critical point at v0, (f o h)(v0) = f(u0), (f o h)''(v0) = f2 h1^2, and
    Lnorm[f o h, (g o h) h'; v0] = h1 Lnorm[f, g; u0].      (2.2)
Exactly what is invariant:
 (i)  the critical value e^{f(u0)/z};
 (ii) the unnormalised expansion: Lfull[f o h, (g o h) h'; v0] = Lfull[f, g; u0] when the roots are related by
      (-(f o h)''(v0))^(1/2) = h1 (-f2)^(1/2): the Jacobian h' in the transformed amplitude is compensated exactly by the
      Jacobian h1 in the Gaussian normalisation (for real data and h1 > 0 both roots are the positive ones; h1 < 0 is the
      orientation reversal of the path);
 (iii) the normalised expansion, with no convention at all: Lhat[f o h, (g o h) h'; v0] = Lhat[f, g; u0], because the
      leading coefficient transforms in the same way, ((g o h) h')(v0) = h1 g(u0).
Proof. For eps in C let h_eps(v) := u0 + h1 (v - v0) + eps (h(v) - u0 - h1 (v - v0)); then h_eps(v0) = u0, h_eps'(v0) = h1 for
every eps, h_0 is affine and h_1 = h. On a neighbourhood V of v0 (uniform for |eps| <= 2), h_eps' has no zero and h_eps(V)
lies in the domain of f, g. Put f_eps := f o h_eps, g_eps := (g o h_eps) h_eps', n(eps) := Lnorm[f_eps, g_eps; v0]. Every f_eps
has a critical point at v0 with critical value f(u0) and f_eps''(v0) = f2 h1^2, both independent of eps; the Taylor
coefficients of f_eps, g_eps at v0 are polynomials in eps, so by (D2) each z-coefficient of n(eps) is a polynomial in eps.
(2.1) with tau = eps, s*(eps) = v0, f0' = 0, f2' = 0 gives
    n'(eps) = Lnorm[f_eps, d_eps g_eps + g_eps d_eps f_eps/z; v0].
Let eta := d_eps h_eps = h(v) - u0 - h1 (v - v0) and xi := eta/h_eps' (analytic on V). Then d_eps f_eps = (f' o h_eps) eta =
xi f_eps', and d_eps g_eps = (g' o h_eps) eta h_eps' + (g o h_eps) eta' = (xi g_eps)' (because xi g_eps = eta (g o h_eps)). So the
amplitude is (xi g_eps)' + f_eps' (xi g_eps)/z, a total derivative, and n'(eps) = 0 by (P3); hence n(1) = n(0). Finally, with
w = v - v0, f_0(v0 + w) = f(u0 + h1 w) has R_{f_0}(w) = R_f(h1 w) and f_0''(v0) = h1^2 f2, so by (G0) with lambda = h1
    n(0) = h1 <<g(u0 + h1 w) Exp[R_f(h1 w)/z]>>_{h1^2 f2} = h1 <<g(u0 + v) Exp[R_f(v)/z]>>_{f2} = h1 Lnorm[f, g; u0].
(ii) and (iii) are restatements of (2.2) together with (f o h)(v0) = f(u0). QED

### 2.6 Analytic meaning (Laplace's method) [analytic]
**Lemma 2.5.** Let I be an open real interval, f and g real-analytic on I, u0 in I with f'(u0) = 0, f''(u0) < 0, and assume
sup_{I - (u0 - e, u0 + e)} f < f(u0) for every e > 0 and Int_I e^{f/z0} |g| du < oo for some z0 > 0. Then for every M, as z -> 0+,
    Int_I e^{f(u)/z} g(u) du = e^{f(u0)/z} (2 pi z/(-f''(u0)))^(1/2) [ Sum_{m<M} ([z^m] Lnorm[f, g; u0]) z^m + O(z^M) ].
Proof. Localisation: fix e > 0 small and eta > 0 with f <= f(u0) - eta outside (u0 - e, u0 + e). There, for 0 < z < z0,
e^{f/z} = e^{f/z0} e^{f (1/z - 1/z0)} <= e^{f/z0} e^{(f(u0) - eta)(1/z - 1/z0)}, so the outer part is O(e^{(f(u0) - eta)/z}) =
e^{f(u0)/z} O(z^M) for every M. Morse coordinate: f(u) - f(u0) = (u - u0)^2 q(u) with q real-analytic and q(u0) = f''(u0)/2 < 0;
shrinking e, q < 0 on (u0 - e, u0 + e), and w(u) := (u - u0) (2 q(u)/f''(u0))^(1/2) is real-analytic, strictly increasing,
w(u0) = 0, w'(u0) = 1, with f(u) = f(u0) + f''(u0) w^2/2. Let u(w) be the inverse map, [w_-, w_+] (w_- < 0 < w_+) the image
interval and G(w) := g(u(w)) u'(w). Watson: with alpha := -f''(u0) > 0 and G(w) = Sum_{p<2M} G_p w^p + O(|w|^(2M)) on
[w_-, w_+], each Int_{w_-}^{w_+} e^{-alpha w^2/(2z)} w^p dw equals the integral over R (namely (2 pi z/alpha)^(1/2) (p-1)!!
(z/alpha)^(p/2) for p even, 0 for p odd) up to O(e^{-b/z}), b > 0, and Int_R e^{-alpha w^2/(2z)} |w|^(2M) dw = O(z^(M + 1/2)).
So the inner integral is e^{f(u0)/z} (2 pi z/alpha)^(1/2) [Sum_{n<M} G_{2n} (2n-1)!! (z/alpha)^n + O(z^M)], and
Sum_n G_{2n} (2n-1)!! (z/alpha)^n = <<G>>_{-alpha} = Lnorm[f''(u0) w^2/2, G; 0] (purely quadratic phase, R = 0). Finally (2.2)
with h = u(w), h1 = 1, and (D3) for the constant f(u0) give Lnorm[f o u, (g o u) u'; 0] = Lnorm[f, g; u0]. QED
(Lemma 2.5 is proved in full above; what follows is a cross-reference only. The same coefficient formula is Hormander's L_j
in the stationary phase theorem, ALPDO I, Thm 7.7.5, which is stated for integrals of u e^{i omega F} with Im F >= 0; the
Laplace case is its special case omega = 1/z, F := -i (f - f(u0)), for which e^{f/z} = e^{f(u0)/z} e^{i omega F} and
Im F = f(u0) - f >= 0 near u0. Lemma 2.5 is used only for the Beta integral in section 6; everything else in sections 2-5 and
7 is formal.)

### 2.7 Homogeneity of the Wick coefficients
For m >= 0 write wickZ_m(W_2, W_3, ...) := wickZ[[m+1]] (section 0); it is a polynomial in W_2^-1, W_3, ..., W_{2m+2} with
rational coefficients (at order z^m no part exceeds 2m + 2, section 0).
**Lemma 2.6 (homogeneity).** Every monomial of wickZ_m is a rational multiple of W_2^-n prod_{i=1}^{n-m} W_{k_i} with all
k_i >= 3 and Sum_i k_i = 2n. Consequently, for every lambda != 0, as an identity of rational functions of the W_k,
    wickZ_m(lambda^2 W_2, lambda^3 W_3, ..., lambda^(2m+2) W_{2m+2}) = wickZ_m(W_2, W_3, ..., W_{2m+2}).                 (2.3)
In particular, for f as in 2.2 and lambda != 0: Lnorm[f(u0 + lambda .), 1; 0] = Lnorm[f, 1; u0].
Proof. By section 0, the z^m coefficient is Sum_n (2n)! (-1)^n/(n! 2^n W_2^n) Sum_p prod_k W_k^(m_k)/((k!)^(m_k) m_k!), p running
over the partitions of 2n into exactly n - m parts >= 3; so its monomials are W_2^-n prod_{k in p} W_k with Sum_{k in p} k = 2n.
Under W_k -> lambda^k W_k such a monomial is multiplied by lambda^(-2n) lambda^(2n) = 1. For the second statement put
f_lambda(v) := f(u0 + lambda v); then f_lambda'(0) = 0 and f_lambda^(k)(0) = lambda^k f^(k)(u0), so by (D1) the z^m coefficients
of Lnorm[f_lambda, 1; 0] and of Lnorm[f, 1; u0] are wickZ_m at the arguments lambda^k f^(k)(u0) and f^(k)(u0), which agree by
(2.3). (Equivalently: (G0) with K(v) = Exp[R_f(v)/z], since R_{f_lambda}(v) = R_f(lambda v) and f_lambda''(0) = lambda^2 f2.) QED
In words: if W_k is given the weight k, every Wick monomial has weight -2n + Sum_i k_i = 0, so the Wick formula is blind to
linear rescalings of the coordinate. This is what makes the degenerate coordinate y = Log U harmless (Proposition 3.1(ii),
Remark 3.2).

## 3. The Wick row is the normalised expansion of J (Step 2(a)) [formal]

**Proposition 3.1.** (i) For every tau in D' - {0} (i.e. every small t != 0), coefficientwise in z,
    r^pt_{0,1}(z; t) = N(z, tau) := Lhat[phi_tau, (1 - tau s)^-1; s*(tau)] = (1 - tau s*) Lnorm[phi_tau, (1 - tau s)^-1; s*(tau)].   (3.1)
(ii) For every m >= 0 and every tau in D' - {0},
    [z^m] r^pt_{0,1}(z; t) = wickZ_m(W_2(U1), W_3(U1), ...) = wickZ_m(w_2(tau), w_3(tau), ..., w_{2m+2}(tau)),             (3.2)
a polynomial in 1/w_2, w_3, ..., w_{2m+2} with rational coefficients. Hence the z^m coefficient of r^pt_{0,1}, an algebraic
function of Q, extends analytically to all of D', tau = 0 included, and its value at tau = 0 is A_m := [z^m] A(z), where
A(z) := Lnorm[phi_0, 1; 1/c].
(iii) The z^m coefficient N_m(tau) of N is analytic on D', and N(z, 0) = A(z). So, with the extension of (ii), (3.1) holds on
all of D' (tau = 0 included), N lies in O(D')[[z]], a subring of C{tau}[[z]], and the Taylor series at t = 0 of the z^m
coefficient of r^pt_{0,1} (which is what Series[... /. Q -> Exp[t], {t, 0, M}] in evalQseries and in the check scripts
computes) is that of N_m(tau(t)). So (3.1) holds in C[[z]][[t]].
Proof. (1) The code's row is a normalised expansion in y = Log[U]. Fix tau in D' - {0}. Then U1 = 1 - tau s* != 0, so
y1 := Log[U1] (any branch) is defined and Wy(y) := W(e^y) is analytic near y1, with Wy^(k)(y1) = (theta^k W)(U1) = W_k(U1).
Wy'(y1) = W_1(U1) = 0 and Wy''(y1) = W_2(U1) = -delta/(x tau^2 s*^2) != 0 by (1.5). In this coordinate the integrand is
e^{W/z} dU/U = e^{Wy(y)/z} dy, amplitude 1. By (D1) and the description of wickZ in section 0,
    r^pt_{0,1}(z; t) = Lnorm[Wy, 1; y1] = Lhat[Wy, 1; y1]
(to every order the code keeps; its truncation at z^ZMAX is the truncation of this series).
(2) Change of coordinate; proof of (i). h(s) := Log[1 - tau s], with the branch for which h(s*) = y1, is analytic near s*, and
h'(s*) = -tau/(1 - tau s*) != 0 because tau != 0. By (1.1), Wy o h = phi_tau + const near s*, and the amplitude 1 becomes
(1 o h) h' = -tau (1 - tau s)^-1. By (P1)(iii), (D3) (constant in the phase) and linearity (the factor -tau),
    Lhat[Wy, 1; y1] = Lhat[phi_tau, -tau (1 - tau s)^-1; s*] = Lhat[phi_tau, (1 - tau s)^-1; s*] = N(z, tau).
(3) Proof of (ii). The first equality in (3.2) is the definition of the code's row (section 0). By (2.3) with lambda = tau
(tau != 0), wickZ_m(W_k(U1)) = wickZ_m(tau^k W_k(U1)) = wickZ_m(w_k), the second equality. By Lemma 1.1 every w_k is analytic
on D' and w_2 has no zero there, so the right side of (3.2) is analytic on D'. At tau = 0, w_k(0) = (-1)^k phi_0^(k)(1/c)
(Lemma 1.1), and (2.3) with lambda = -1, followed by (D1) (phi_0'(1/c) = 0, phi_0''(1/c) = -c^3/x != 0, 1.3), gives
    wickZ_m(w_k(0)) = wickZ_m((-1)^k phi_0^(k)(1/c)) = wickZ_m(phi_0^(k)(1/c)) = [z^m] Lnorm[phi_0, 1; 1/c] = A_m.
No coordinate change, i.e. no use of (P1), occurs in (3).
(4) Proof of (iii). By (D2), N_m(tau) is a polynomial in 1/phi''_tau(s*), in phi_tau^(k)(s*) (k <= 2m + 2), in
d_s^a (1 - tau s)^-1 at s* = a! tau^a (1 - tau s*)^(-a-1) (a <= 2m), and in the factor 1 - tau s*. By (1.9) all these are
analytic on D', because s* is analytic there, s* is not in {0, 1, 1/tau}, and phi''_tau(s*) != 0 (1.2, 1.3). At tau = 0,
s*(0) = 1/c and (1 - 0 s)^-1 = 1, so N(z, 0) = Lnorm[phi_0, 1; 1/c] = A(z). By (i) the analytic functions of (ii) and N_m agree
on D' - {0}, and by (ii) and N(z, 0) = A(z) they also agree at tau = 0 (directly; continuity would give the same). Finally
t -> tau = 1 - e^-t is analytic with tau(0) = 0, so the z^m coefficient of r^pt_{0,1} is an analytic function of t near 0,
equal to N_m(tau(t)), and the two Taylor series coincide. QED

**Remark 3.2 (why the degenerate coordinate is harmless; a uniform version).** In the coordinate y the inputs of the Wick
formula blow up at t = 0: W_k(U1) = tau^-k w_k = O(tau^-k) and W_2(U1)^-1 = O(tau^2) (Lemma 1.1); also h'(s*) = O(tau) in
step (2). Nevertheless every single Wick monomial is analytic at t = 0: a monomial W_2^-n prod_i W_{k_i} of wickZ_m has
Sum_i k_i = 2n (Lemma 2.6), so it equals w_2^-n prod_i w_{k_i}; the factor tau^(2n) coming from W_2^-n and the factor
tau^(-2n) coming from the numerator cancel inside the monomial. No cancellation between different monomials, and no
reparametrisation invariance, is needed for the analyticity; (P1) is needed only to identify the row with the expansion N of
the Euler integral J.
Geometrically, y is -tau times a coordinate that is good uniformly near tau = 0. Let
    ytilde := -Log[1 - tau s]/tau = Sum_{n>=1} tau^(n-1) s^n/n,
    s = sigma_tau(ytilde) := (1 - e^{-tau ytilde})/tau = Sum_{n>=1} (-tau)^(n-1) ytilde^n/n!
(so ytilde = s at tau = 0). sigma is entire in (tau, ytilde), sigma_0 = id, d sigma_tau/d ytilde = e^{-tau ytilde} =
1 - tau sigma_tau(ytilde), and U = 1 - tau s = e^{-tau ytilde}, i.e. y = -tau ytilde and d/d ytilde = -tau theta. By (1.1),
    W(e^{-tau ytilde}) = phi_tau(sigma_tau(ytilde)) + c Log[tau] + x Log[Q] + x Log[-1],                                   (3.3)
and the first term on the right is analytic in (tau, ytilde) near (0, 1/c) (there s = sigma_tau(ytilde) is close to 1/c, so
s, 1 - s and 1 - tau s do not vanish): it has no singularity at tau = 0. Put ytilde1(tau) := -Log[U1]/tau (principal branch;
analytic near tau = 0, with ytilde1(0) = s*(0) = 1/c and sigma_tau(ytilde1) = s*). Differentiating (3.3) k >= 2 times at
ytilde1, with d/d ytilde = -tau theta on the left,
    (-tau)^k W_k(U1) = (-1)^k w_k(tau) = d^k/d ytilde^k [phi_tau o sigma_tau](ytilde1(tau))     (k >= 2),                  (3.4)
which re-proves Lemma 1.1 near tau = 0 (at tau = 0 the right side is phi_0^(k)(1/c), as sigma_0 = id). For k = 1 both sides
of the differentiated (3.3) vanish: W_1(U1) = 0, and ytilde1 is a critical point of phi_tau o sigma_tau since phi'_tau(s*) = 0.
Uniform version of Proposition 3.1. Let tau be small (tau = 0 allowed) and choose y1 := -tau ytilde1 = Log[U1]. For tau != 0,
Lemma 2.6 in the form Lnorm[f(u0 + lambda .), 1; 0] = Lnorm[f, 1; u0] with f = Wy, u0 = y1, lambda = -tau (note
Wy(y1 - tau v) = W(e^{-tau (ytilde1 + v)})), then (3.3) and (D3), give
    r^pt_{0,1}(z; t) = Lnorm[Wy, 1; y1] = Lnorm[W(e^{-tau .}), 1; ytilde1] = Lnorm[phi_tau o sigma_tau, 1; ytilde1].
The right side is defined for every small tau, tau = 0 included: phi_tau o sigma_tau has a critical point at ytilde1 with second
derivative phi''_tau(s*) U1^2 = w_2 != 0 (-c^3/x at tau = 0), and its z-coefficients are analytic in tau by (D2), because its
Taylor coefficients at ytilde1(tau) are derivatives of a function analytic in (tau, ytilde), evaluated at the analytic point
ytilde1(tau), and 1/w_2 is analytic. Now (P1) with
h = sigma_tau, v0 = ytilde1, u0 = s*, h1 = e^{-tau ytilde1} = U1 and g = (1 - tau s)^-1, for which (g o h) h' =
e^{-tau ytilde}/(1 - tau sigma_tau(ytilde)) = 1, gives at every small tau, tau = 0 included (where sigma_0 = id and both sides
are A(z)),
    Lnorm[phi_tau o sigma_tau, 1; ytilde1] = U1 Lnorm[phi_tau, (1 - tau s)^-1; s*] = N(z, tau).
So for tau != 0 the code's row is the restriction of the family Lnorm[phi_tau o sigma_tau, 1; ytilde1], which is analytic at
tau = 0 and equal to N there as well: this is Proposition 3.1 near tau = 0 without a continuation step. [C6 checks (3.4) for
k <= 12, and wickZ_m(w_k) = Lnorm[phi_tau o sigma_tau, 1; ytilde1] = N mod (tau^7, z^5) as tau-series, with value
A(z) = h_w c/(c+z) at tau = 0, at x = 2, 1/3, -5/2.]
[C1: (3.1) through z^5 in exact rationals at (x, tau) = (2, 5/8) (Q = 8/3, U1 = 3/4, s* = 2/5, delta = 2) and (1/3, 7/12)
(Q = 12/5, U1 = 1/2, s* = 6/7, delta = 1), the y-side computed exactly as the code does (W_k = theta^k W at U1, amplitude 1);
at (2, 5/8) both sides are 1, -97/512, 20809/524288, 1164725/268435456, ...]

## 4. Gauss's equation for the formal expansion (Step 2(b)) [formal]

### 4.1 The amplitude identity
Put f := phi_tau, g := (1 - tau s)^-1, and for an amplitude q(tau, s) (Laurent polynomial in z with coefficients analytic
in (tau, s)) let  delta_tau q := d_tau q + q (d_tau f)/z,  so that e^{f/z} delta_tau q = d_tau (e^{f/z} q). Since
d_tau f = s/(1 - tau s),
    delta_tau g = (1 + 1/z) s (1 - tau s)^-2,     delta_tau^2 g = (1 + 1/z)(2 + 1/z) s^2 (1 - tau s)^-3.
Dividing (1.8) by e^{f/z} = s^(1/z) (1 - s)^(x/z) (1 - tau s)^(-1/z) gives, with k := s (1 - s)(1 - tau s)^-2,
    tau (1 - tau) delta_tau^2 g + (gamma - (alpha + beta + 1) tau) delta_tau g - alpha beta g = -(1 + 1/z) (k' + f' k/z).   (4.1)
This is an identity between Laurent polynomials in z with coefficients rational in (s, tau), and it can be checked directly
without reference to (1.8): since f' k = (x tau s^2 - c s + 1)(1 - tau s)^-3, after multiplication by (1 - tau s)^3/(1 + 1/z)
both sides equal (1/z)(-x tau s^2 + c s - 1) + (2 s - tau s - 1). [C4 checks (4.1) exactly at 6 random rational points
(x, s, tau, z).]

### 4.2 Transfer to the formal expansion
**Proposition 4.1.** Let L(tau) := Lnorm[phi_tau, (1 - tau s)^-1; s*(tau)] (an element of C{tau}[[z]]) and
    M(z, tau) := Exp[dphi(tau)/z] (phi''_tau(s*)/phi''_0(1/c))^(-1/2) iota(L(tau))   in C((z))[[tau]]
(the root with value 1 at tau = 0). Then, with ' = d/dtau,
    Gauss[M] := tau (1 - tau) M'' + (gamma - (alpha + beta + 1) tau) M' - alpha beta M = 0   in C((z))[[tau]],
and M(z, 0) = Lnorm[phi_0, 1; 1/c] = A(z) (notation of Proposition 3.1), an element of 1 + z Q(x)[[z]]. In terms of N:
M = Exp[dphi/z] N/alg.
Proof. Here f0(tau) = phi_tau(s*(tau)), so f0' = dphi', and f2(tau) = phi''_tau(s*(tau)). Let Dop := d/dtau + dphi'/z -
f2'/(2 f2), acting on C{tau}((z))-valued expressions. (P2), applied to each z-power of an amplitude q, says
Dop Lnorm[q] = Lnorm[delta_tau q]; applying it twice and using (D3) (linearity over C[tau][z, 1/z]),
    [tau (1 - tau) Dop^2 + (gamma - (alpha + beta + 1) tau) Dop - alpha beta] L = Lnorm[LHS of (4.1)]
                                                                                = -(1 + 1/z) Lnorm[k' + f' k/z] = 0
by (P3). Let E(tau) := Exp[dphi/z] (f2(tau)/f2(0))^(-1/2) in C[1/z][[tau]]; in C((z))[[tau]], d/dtau (E Y) = E Dop Y because
E' = E (dphi'/z - f2'/(2 f2)) (iota commutes with d/dtau and with multiplication by convergent functions of tau, and
Exp[dphi/z]' = (dphi'/z) Exp[dphi/z]). Hence Gauss[E iota(L)] = E iota([tau(1 - tau) Dop^2 + ...] L) = 0.
At tau = 0: s*(0) = 1/c, (1 - 0 s)^-1 = 1 and dphi(0) = 0, so M(z, 0) = Lnorm[phi_0, 1; 1/c], which lies in 1 + z Q(x)[[z]]
by (D2) since phi_0^(k)(1/c) = (k - 1)! c^k ((-1)^(k-1) - x^(1-k)) is in Q(x) (1.9). Finally N = (1 - tau s*) L and
alg = (1 - tau s*)(f2(tau)/f2(0))^(1/2), so M = Exp[dphi/z] N/alg. QED
No boundary term occurs at any stage: the formal functional kills total derivatives exactly (P3), the formal shadow of
the vanishing boundary terms of 1.4, valid for all x and z without (H_E).

## 5. Uniqueness and the closed form (Step 2(c)) [formal]

**Lemma 5.1 (uniqueness).** If Y = Sum_k Y_k tau^k in C((z))[[tau]] satisfies Gauss[Y] = 0 (alpha = beta = 1/z + 1,
gamma = c/z + 2, c != 0), then Y = Y_0 F.
Proof. With Y' = Sum (k+1) Y_{k+1} tau^k, Y'' = Sum (k+2)(k+1) Y_{k+2} tau^k, the coefficient of tau^k of Gauss[Y] is
    (k+1) k Y_{k+1} - k (k-1) Y_k + gamma (k+1) Y_{k+1} - (alpha + beta + 1) k Y_k - alpha beta Y_k
        = (k+1)(k + gamma) Y_{k+1} - (k + alpha)(k + beta) Y_k.
Now k + gamma = (c + (k+2) z)/z is a unit of C((z)) (inverse (z/c)(1 + (k+2) z/c)^-1 in C[[z]]) and k + 1 != 0, so
Y_{k+1} = (k + alpha)(k + beta) Y_k/((k+1)(k + gamma)), whence Y_k = Y_0 (alpha)_k (beta)_k/((gamma)_k k!) = Y_0 F_k. QED
So the solutions of Gauss's equation in C((z))[[tau]] form the line C((z)) F; the second local solution, of exponent
1 - gamma = -c/z - 1 (the one that would carry tau^(-c/z)), is not in this ring.

**Theorem 5.2.** In C((z))[[tau]]:
    N(z, tau) = A(z) F(z, tau) Exp[-dphi(tau)/z] alg(tau),     A(z) = Lnorm[phi_0, 1; 1/c].      (5.1)
Proof. By Proposition 4.1 and Lemma 5.1, M = M(z, 0) F = A F. Multiply by alg Exp[-dphi/z] (alg is in C{tau}, and
Exp[-dphi/z] is the inverse of Exp[dphi/z] in C[1/z][[tau]]) and use M = Exp[dphi/z] N/alg. QED
With Proposition 3.1: r^pt_{0,1}(z; t) = A(z) F(z, tau) Exp[-dphi(tau)/z] alg(tau) in C((z))[[tau]], tau = 1 - e^-t; the left
side lies in C[[z]][[tau]]. [C3: (5.1) checked mod (tau^5, z^5) at x = 2 and x = 1/3, N computed from the Wick formula in
the coordinate s with tau-series coefficients, A taken from Stirling (Proposition 6.1); also the identities (a), (c) of 1.2,
(1.3) and (1.4) as tau-series. C6: the left side computed as the code does, in y = Log U, via (3.2), equals N mod
(tau^7, z^5) at x = 2, 1/3, -5/2.]

## 6. The constant A(z) = h_w(z) c/(c + z)

**Proposition 6.1.** A(z) := Lnorm[phi_0, 1; 1/c] = h_w(z) c/(c + z) in Q(x)[[z]], with h_w the series of hvc_lib.m for
w = (c, -x, -1) (section 0). In words: the normalised Laplace expansion of the Beta integral B(1/z + 1, x/z + 1) =
Int_0^1 e^{phi_0(s)/z} ds is the Stirling quotient, i.e. the Stirling tail h_w of B(x/z, 1/z) times the elementary factor
c/(c + z) coming from B(a + 1, b + 1) = B(a, b) a b/((a + b)(a + b + 1)).
Proof. [analytic for x > 0, then algebraic in x]
(1) Reduction to x > 0. Both sides have z-coefficients in Q(x): the left by (D2), because phi_0^(k)(1/c) =
(k - 1)! c^k ((-1)^(k-1) - x^(1-k)) and 1/phi_0''(1/c) = -x/c^3; the right by definition (c/(c + z) = Sum_n (-z/c)^n).
Two rational functions of x that agree on (0, oo) are equal. Fix a real x > 0 from now on.
(2) Laplace. For real z > 0, B(1/z + 1, x/z + 1) = Int_0^1 s^(1/z) (1 - s)^(x/z) ds = Int_0^1 e^{phi_0(s)/z} ds. On I = (0, 1),
phi_0 = Log[s] + x Log[1 - s] is real-analytic, phi_0'' = -1/s^2 - x/(1 - s)^2 < 0, phi_0' = (1 - c s)/(s (1 - s)) vanishes
only at s = 1/c, and phi_0 -> -oo at both ends; so phi_0 is strictly concave with its unique maximum at 1/c, the
supremum of phi_0 off (1/c - e, 1/c + e) is max(phi_0(1/c - e), phi_0(1/c + e)) < phi_0(1/c), and Int_0^1 e^{phi_0} ds =
B(2, x + 1) < oo. Lemma 2.5 (g = 1, z0 = 1), with -phi_0''(1/c) = c^3/x and phi_0(1/c) = -Log c + x Log(x/c) = x Log x - c Log c:
    B(1/z + 1, x/z + 1) = e^{(x Log x - c Log c)/z} (2 pi z x/c^3)^(1/2) [Sum_{m<M} A_m z^m + O(z^M)]   (z -> 0+), A = Sum A_m z^m.
(3) Stirling. Gamma(u + 1) = u Gamma(u) gives B(a + 1, b + 1) = B(a, b) a b/((a + b)(a + b + 1)); for a = 1/z, b = x/z,
a b/((a + b)(a + b + 1)) = (x/z^2)/((c/z)((c + z)/z)) = x/(c (c + z)). Stirling's series with remainder (the series is DLMF
5.11.1; for the remainder see DLMF s.5.11(ii): for real positive argument the remainder after the term k = K is bounded in
absolute value by the first neglected term, B_{2K+2}/((2K+2)(2K+1) u^(2K+1)), and has the same sign; also Whittaker-Watson
s.12.33): for u -> +oo and every K,
    Log Gamma(u) = (u - 1/2) Log u - u + (1/2) Log(2 pi) + Sum_{k=1}^K B_2k/(2k (2k-1) u^(2k-1)) + O(u^(-2K-1)).
Apply it to Log B(1/z, x/z) = Log Gamma(1/z) + Log Gamma(x/z) - Log Gamma(c/z) (all three arguments -> +oo). The terms -u
cancel (1/z + x/z = c/z); the coefficient of Log z in the (u - 1/2) Log u terms is -(1/z - 1/2) - (x/z - 1/2) + (c/z - 1/2)
= 1/2, and the remaining ones give (x Log x - c Log c)/z - (1/2) Log x + (1/2) Log c. With u^(1-2k) = z^(2k-1) w^(1-2k) for
u = w/z,
    B(1/z, x/z) = e^{(x Log x - c Log c)/z} (2 pi c z/x)^(1/2) Exp[Sum_{k=1}^K B_2k/(2k (2k-1)) z^(2k-1) (1 + x^(1-2k) - c^(1-2k))
                  + O(z^(2K+1))].
The exponent is the truncation of Log h_w (section 0) up to O(z^(2K+1)), and the exponential of (polynomial + O(z^(2K+1)))
is the truncation of the formal exponential + O(z^(2K+1)). Multiplying by x/(c (c + z)) (analytic at z = 0) and using
(2 pi c z/x)^(1/2) x/c = (2 pi z x/c^3)^(1/2) c:
    B(1/z + 1, x/z + 1) = e^{(x Log x - c Log c)/z} (2 pi z x/c^3)^(1/2) [ (h_w(z) c/(c + z))_{<M} + O(z^M) ].
(4) The prefactors in (2) and (3) are the same positive function of z, so dividing by it, the two power series are
asymptotic expansions of the same function as z -> 0+; asymptotic power series expansions are unique, so they coincide. QED
Remarks. (a) z^1 by hand: (D1) gives A_1 = f4/(8 f2^2) - 5 f3^2/(24 f2^3) with f2 = -c^3/x, f3 = 2 c^3 (x^2 - 1)/x^2,
f4 = -6 c^4 (x^3 + 1)/x^3, i.e. A_1 = (x^2 - 11 x + 1)/(12 c x); and [z^1] h_w c/(c + z) = (1/12)(1/x + 1 - 1/c) - 1/c =
(x^2 - 11 x + 1)/(12 c x). (b) C2 checks Proposition 6.1 exactly through z^7 at x = 2, 1/3, 5 and at x = -3/7 (outside
x > 0, as the rationality argument predicts). (c) This completes the F2 sketch of conifold_regularity_note.md s.1 to all
orders in z: point frame h_w c/(c + z), frame B c/(c + z). (d) A purely formal proof would have to reproduce Stirling's
series inside the formal calculus; the asymptotic argument is the natural one and is the only analytic step of S1. (e) That
the code's row at t = 0 equals A(z) = Lnorm[phi_0, 1; 1/c] is formal (Proposition 3.1(ii)); only the evaluation of A uses
analysis.

## 7. Regularity, the frame-B row and F2 (Step 2(d)) [formal, given Proposition 6.1]

**Corollary 7.1 (regularity).** For every k >= 0, (a) Psi_k(z) := [tau^k] (F Exp[-dphi/z] alg) and (b) H_k(z) :=
[tau^k] (F Exp[-dphi/z]) are rational functions of z, with coefficients in Q(x), without pole at z = 0.
Proof. (a) Psi_k = Sum_{i+j+l=k} F_i(z) ([tau^j] Exp[-dphi/z]) ([tau^l] alg) is a finite sum of products of elements of Q(x)(z),
of Q(x)[1/z] and of Q(x) (alg and dphi have Taylor coefficients in Q(x), by (1.4), 1.2(c) and the binomial series), so it
lies in Q(x)(z). By (5.1), Psi_k = [tau^k] (N/A) in C((z)): A is a unit of C[[z]] (constant term 1, (D2)) and N lies in
C[[z]][[tau]] (Proposition 3.1), so the Laurent expansion of Psi_k at z = 0 has no negative powers, i.e. Psi_k is regular at 0.
(Only A in 1 + z C[[z]] is used here, not its value.) (b) alg lies in 1 + tau Q(x)[[tau]] and does not depend on z, so 1/alg
does too, and H_k = Sum_{j=0}^{k} Psi_j [tau^(k-j)] (1/alg) is a Q(x)-linear combination of Psi_0, ..., Psi_k. QED
Part (b) is the form used in Step 4 (S3, (I1)): there G = rho r_{0,1} = (c/(c + z)) (tau/t) F Exp[-dphi/z] contains
F Exp[-dphi/z] without alg (rho = tau/(t alg), Step 3). Section 7.4 gives (b) directly: Exp[-dphi/z] F = Phi in O(D)[[z]] (7.3).

**Theorem 7.2 (frame-B row).** In C[[z]][[t]], and coefficientwise in z as analytic functions of t near 0 (made precise at the
end of the proof),
    r_{0,1}(z; t) = (c/(c + z)) F(z, tau) Exp[-dphi(tau)/z] alg(tau),     tau = 1 - e^-t.
Precisely: for every k, [t^k] r_{0,1}(z; t) = Sum_{j<=k} ([t^k] tau(t)^j) (c/(c + z)) Psi_j(z), a rational function of z regular
at z = 0, whose Taylor coefficients in z are the t^k Taylor coefficients of the z-coefficients of r_{0,1}.
Proof. r_{0,1} = h_w^-1 r^pt_{0,1} (frame B; h_w^-1 is t-independent, so it commutes with the t-expansion) = h_w^-1 N
(Proposition 3.1) = h_w^-1 A F Exp[-dphi/z] alg (5.1) = (c/(c + z)) F Exp[-dphi/z] alg (Proposition 6.1). The passage from
tau to t is the triangular substitution tau = t - t^2/2 + ..., with rational coefficients, and preserves regularity.
Coefficientwise meaning. The right side is a priori an element of C((z))[[tau]], so its z-coefficients are not yet functions
of t. By Corollary 7.1 it lies in C[[z]][[tau]], and the chain above says that it equals iota(h_w^-1 N), where h_w^-1 N lies in
O(D')[[z]] (Proposition 3.1(iii); h_w^-1 is in C[[z]], so each z-coefficient of h_w^-1 N is a finite combination of the N_m).
So each z-coefficient of the right side is the Taylor series at 0 of a function analytic on the disc D', hence converges on D',
and its sum is the corresponding z-coefficient of h_w^-1 N = r_{0,1}. By (7.3) below the right side also equals
(c/(c + z)) Phi alg with Phi in O(D)[[z]], so these sums continue analytically to D. Hence for every m and every t with tau(t)
in D' (in particular every small t), [z^m] r_{0,1}(z; t) = [z^m] ((c/(c + z)) Phi(z, tau(t)) alg(tau(t))). QED
Consistency: multiplied by rho = tau/(t alg) (Step 3) this is the KEY FORMULA V3, verified against the stored genus-4 factors
for t^k, k <= 6, z^m, m <= 9. Independently of the factor files, C3 checks [tau^0] h_w^-1 N = c/(c + z) and
[tau^1] h_w^-1 N = (1-x)^2 z/(c (c + z)(c + 2z)) through z^4 (the latter is the O(t) correction of leaf_structure.md, Step 0).

**Corollary 7.3 (F2).** r_{0,1}(z; 0) = c/(c + z), and r^pt_{0,1}(z; 0) = A(z) = h_w(z) c/(c + z).
Proof. t = 0 means tau = 0, and F(z, 0) = 1, dphi(0) = 0, alg(0) = 1. (The equality r^pt_{0,1}(z; 0) = A(z) is also
Proposition 3.1(ii), obtained there directly from the homogeneity of the Wick formula, without (P1).) QED

### 7.4 Purely formal alternative for Corollary 7.1: eikonal and transport equations
This argument uses neither the Laplace identification (sections 2.5-3, 6) nor the value of A; it shows directly that the
z-expansion of Exp[-dphi/z] F exists and has coefficients analytic in tau.
Write F = Exp[S/z] Phi with S in tau C{tau}. From F' = e^{S/z}(Phi' + S' Phi/z), F'' = e^{S/z}(Phi'' + 2 S' Phi'/z + S'' Phi/z +
S'^2 Phi/z^2) and gamma = c/z + 2, alpha + beta + 1 = 2/z + 3, alpha beta = 1/z^2 + 2/z + 1 one finds, as an identity of
operators on C((z))[[tau]],
    Exp[-S/z] Gauss[Exp[S/z] Phi] = z^-2 Eik(S') Phi + z^-1 (P Phi' + R Phi) + T Phi,      (7.1)
    Eik(S') := tau (1 - tau) S'^2 + (c - 2 tau) S' - 1,      P := 2 tau (1 - tau) S' + c - 2 tau,
    R := tau (1 - tau) S'' + (2 - 3 tau) S' - 2,             T := tau (1 - tau) d^2/dtau^2 + (2 - 3 tau) d/dtau - 1.
Eikonal equation. Eik(S') = 0 is quadratic in S'; its root regular at tau = 0 has S'(0) = 1/c, the other root is
-c/tau + O(1) (the exponent 1 - gamma). Claim: S := dphi solves it, with
    d(dphi)/dtau = s*/(1 - tau s*) = (1 - x s*)/(1 - tau),     d(dphi)/dt = (1 - tau) d(dphi)/dtau = 1 - x s*.
Indeed d(dphi)/dtau = (d_tau phi_tau)(s*) + phi'_tau(s*) s*' = s*/(1 - tau s*) (envelope; d_tau phi_tau(s) = s/(1 - tau s)),
which equals (1 - x s*)/(1 - tau) by 1.2(c); and dtau/dt = 1 - tau. Explicitly, with S' = (1 - x s)/(1 - tau), s = s*,
    (1 - tau) Eik = tau (1 - x s)^2 + (c - 2 tau)(1 - x s) - (1 - tau) = x^2 tau s^2 - c x s + (c - 1) = x (x tau s^2 - c s + 1) = 0
(expand: tau - 2 x tau s + x^2 tau s^2 + c - c x s - 2 tau + 2 x tau s - 1 + tau; and c - 1 = x). S'(0) = 1/c: the regular root.
The coefficients. P = d Eik/d S' at the root, so P^2 = (c - 2 tau)^2 + 4 tau (1 - tau) = c^2 - 4 x tau, and P(0) = c: P = delta.
Directly: P - c + 2 tau = 2 tau (1 - tau) S' = 2 tau (1 - x s*), and delta - c + 2 tau = 2 (1 - c s* + tau s*)/s* = 2 tau (1 - x s*)
(delta = 2/s* - c, 1 - c s* = -x tau s*^2). For R: differentiating (1 - tau) S' = 1 - x s* gives tau (1 - tau) S'' = tau S' - x tau s*',
so R = 2 (1 - tau) S' - x tau s*' - 2 = -x s* (2 + x tau s*/delta) (by 1.2(d)) = -x (c + 3 delta)/(delta (c + delta))
(x tau s* = (c - delta)/2, s* = 2/(c + delta)); R(0) = -2x/c. [C5 checks all of this at (x, tau) = (2, 5/8): P = 2, R = -9/5.]
With Eik = 0, multiplying (7.1) by z:
    z Exp[-S/z] Gauss[Exp[S/z] Phi] = P Phi' + R Phi + z T Phi.      (7.2)
Transport equations. For Phi = Sum_{m>=0} z^m Phi_m(tau), (7.2) vanishes iff P Phi_m' + R Phi_m + T Phi_{m-1} = 0 for all
m >= 0 (Phi_{-1} := 0). On D, P = delta has no zero and R is analytic (c + delta != 0), so there are unique analytic solutions
with Phi_0(0) = 1 and Phi_m(0) = 0 (m >= 1):
    Phi_0 = Exp[-Int_0^tau R/P] = 1/alg = 2 c (c/delta)^(1/2)/(c + delta),
    Phi_m(tau) = -Phi_0(tau) Int_0^tau T[Phi_{m-1}](u)/(P(u) Phi_0(u)) du      (m >= 1).
[Phi_0 = 1/alg: by (1.4), (Log alg)' = delta'/(2 delta) + delta'/(c + delta) = delta' (c + 3 delta)/(2 delta (c + delta)), and with
delta' = -2x/delta this is -x (c + 3 delta)/(delta^2 (c + delta)) = R/P.] So Phi lies in O(D)[[z]], a subring of C{tau}[[z]], and
Phi(z, 0) = 1. By (7.1)-(7.2), Yt := Exp[dphi/z] iota(Phi) in C((z))[[tau]] satisfies Gauss[Yt] = 0 and Yt(z, 0) = 1, so
Lemma 5.1 gives Yt = F, i.e.
    Exp[-dphi/z] F = Phi   in C((z))[[tau]],   with [z^0] Phi = 1/alg.      (7.3)
Hence Psi_k = [tau^k](Phi alg) lies in C[[z]], which is Corollary 7.1 again. (7.3) says more: the z-coefficients of
Exp[-dphi/z] F, a priori formal in tau, are analytic on the whole disc D = {|tau| < |c|^2/(4|x|)}; the boundary point
tau = c^2/(4x) is where s* and s** coalesce (delta = 0). Comparison with (5.1): N = A (Phi alg). The Laplace route fixes the
constant A, the transport route the tau-dependence; consistently [z^0] N = A_0 ([z^0] Phi) alg = 1.

### 7.5 Closed forms and the first orders
    delta = (c^2 - 4 x tau)^(1/2),   s* = 2/(c + delta) = (c - delta)/(2 x tau) = (1/c) Sum_n Cat_n (x tau/c^2)^n,   U1 = 1 - tau s*,
    phi''_tau(s*) = -delta/(x s*^2 (1 - tau s*)^2),   W_2(U1) = -delta/(x tau^2 s*^2),   alg = (delta/c)^(1/2) (c + delta)/(2c),
    d(dphi)/dtau = s*/(1 - tau s*) = (1 - x s*)/(1 - tau),   d(dphi)/dt = 1 - x s*.
Series: 1 - x s* = 1/c - x^2 tau/c^3 - 2 x^3 tau^2/c^5 + O(tau^3); with tau = t - t^2/2 + t^3/6 + O(t^4) and c^2 - 4x = (1 - x)^2,
    dphi(t) = t/c - x^2 t^2/(2 c^3) + (1 - x)^2 x^2 t^3/(6 c^5) + O(t^4),
in agreement with V3 and the series printed in hvc_leaf_proof_check_out.txt; alg = 1 - 2 x tau/c^2 + O(tau^2).
Order tau^1, by hand: F_1 = (1 + z)^2/(z (c + 2z)), [tau^1] Exp[-dphi/z] = -1/(c z), [tau^1] alg = -2x/c^2, so
    Psi_1 = (1 + z)^2/(z (c + 2z)) - 1/(c z) - 2x/c^2 = (2x + c z)/(c (c + 2z)) - 2x/c^2 = (1 - x)^2 z/(c^2 (c + 2z)):
the poles at z = 0 of F_1 and of Exp[-dphi/z] cancel, and [t^1] r_{0,1} = [tau^1] r_{0,1} = (1 - x)^2 z/(c (c + z)(c + 2z)).
(The pole structure of the Psi_k at z = -c/j and the behaviour at z = oo are the subject of Step 4, not of this section.)

### 7.6 The constant kappa(z) of outline Step 5 (notation only; not used in S1)
Outline Step 5 writes s_0 := e^{u_1/z} Delta_1^(-1/2) r_{0,1} = kappa(z) Q^{x/z} tau^{c/z+1} F with a t-independent kappa(z).
This kappa is not the phase constant kappa_ph(z) = -e^{x Log[-1]/z} of (1.2). Their relation, with the exponentials read as
formal symbols (S2 s.5.1 makes them precise): u_1 := W(U1), and (1.1) at s = s* gives
    u_1 = phi_tau(s*) + c Log[tau] + x Log[Q] + x Log[-1] = phi_0(1/c) + dphi(tau) + c Log[tau] + x Log[Q] + x Log[-1];
by V1 (Delta_1 = x^2 W_2(U1), proved in S2, Step 3), (1.5) and (1.4), Delta_1 = -x delta/(tau^2 s*^2) = -x c^3 alg^2/tau^2, so
Delta_1^(-1/2) = tau/(kappa_0 alg) with kappa_0^2 = -x c^3 (either root, the same one in Delta_1^(-1/2) and in kappa); and
r_{0,1} is given by Theorem 7.2. In the product the factors e^{dphi/z} and alg cancel:
    s_0 = e^{(phi_0(1/c) + x Log[-1])/z} (c/(kappa_0 (c + z))) Q^{x/z} tau^{c/z+1} F     (frame B).                         (7.4)
So kappa(z) = e^{(phi_0(1/c) + x Log[-1])/z} c/(kappa_0 (c + z)) = -kappa_ph(z) e^{phi_0(1/c)/z} c/(kappa_0 (c + z)) in frame B,
and h_w(z) times this in the point frame (r^pt_{0,1} = h_w r_{0,1}). The sign -1 in front of kappa_ph is the orientation
factor of dU/U = -tau ds/(1 - tau s), which belongs to the integrand (1.2) but not to s_0. S2 (s.5.1-5.3) uses exponential
symbols that depend only on the t-derivative of the exponent, so the additive constant phi_0(1/c) + x Log[-1] of u_1 drops out
there and its constant is kappa_1(z) := c/(kappa_0 (c + z)) (this is why the phase constant of (1.2) is not called kappa_1
here). In every normalisation kappa is t-independent, which is all that Step 5 uses.

## 8. Status, dependencies, checks

Elementary (explicit algebra, and analyticity of explicit functions of tau on D or D'): 1.1-1.3, 1.5 (Lemma 1.1).
Formal (no analytic input; valid for every x with x != 0, c != 0): Lemma 2.1, the definitions 2.2, Propositions 2.2-2.4
(P3, P2, P1), Lemma 2.6, Proposition 3.1 and Remark 3.2 (their analytic content is only the analyticity of finitely many
explicit algebraic functions of tau on D': Lemma 1.1 for the code's side, (D2) and 1.2-1.3 for N), the identity (4.1),
Proposition 4.1, Lemma 5.1, Theorem 5.2, Corollary 7.1, Theorem 7.2 and Corollary 7.3 (given Proposition 6.1), and 7.4
entirely.
Analytic: 1.4 (Euler's integral and Gauss's equation for J under (H_E); motivation only, used nowhere), Lemma 2.5 (Laplace's
method for real data, z -> 0+), Proposition 6.1 (proved for real x > 0 from Lemma 2.5 and Stirling's series with remainder,
then extended to all x because both sides have coefficients in Q(x)).
Logical chain for Theorem S1: (a) = 1.1-1.3 + Lemma 1.1; (b) = Prop. 3.1 [(i): P1; (ii): Lemma 1.1 + Lemma 2.6 + (D1);
(iii): (D2)]; (c) = Prop. 4.1 [(4.1), P2, P3] + Lemma 5.1 + Prop. 6.1 [Lemma 2.5 <- P1]; (d) = Cor. 7.1 [(c), Prop. 3.1, (D2)]
or 7.4 [(7.1), Lemma 5.1] + Thm 7.2 + Cor. 7.3.
Handed on to Steps 3-10: the closed form r_{0,1} = (c/(c + z)) F Exp[-dphi/z] alg (Theorem 7.2, in C[[z]][[t]] and
coefficientwise as analytic functions of t near 0), the regularity at z = 0 of every tau^k coefficient of F Exp[-dphi/z] alg
and of F Exp[-dphi/z] (Corollary 7.1(a), (b); (b) is the form used in Step 4), F2 (Corollary 7.3), and the closed forms of
1.2, 1.3, 1.5 and 7.5 (s*, delta, phi''_tau(s*), W_2(U1), w_k, alg, d(dphi)/dt = 1 - x s*). Not proved here: V1
(Delta_1 = x^2 W_2(U1)), hence F1 and rho = tau/(t alg), which are Step 3 (S2).
Exact checks run for this section (leaf_proof/S1_checks.pl, perl + Math::BigRat, 50 s, 62 checks, all pass; output
leaf_proof/S1_checks_out.txt):
 C1 (3.1) at two points with t != 0, through z^5, the y-side computed as the code does, and (1.5) there;
 C2 Proposition 6.1 through z^7 at x = 2, 1/3, 5, -3/7;
 C3 (5.1) mod (tau^5, z^5) at x = 2, 1/3, no negative z-powers on the right (Corollary 7.1), the frame-B row at tau^0 and tau^1,
    and the identities 1.2(c), (1.3), (1.4), x tau s*^2 - c s* + 1 = 0, phi'_tau(s*) = 0 as tau-series;
 C4 the amplitude identity (4.1) at 6 random rational points (x, s, tau, z); draws at which (4.1) is undefined, or trivially
    0 = 0 (z = -1), are discarded and redrawn (10 draws for 6 points);
 C5 eikonal equation, P = delta, the closed form of R, s*' = x s*^2/delta and (Log alg)' = R/P at (x, tau) = (2, 5/8);
 C6 at x = 2, 1/3, -5/2, as tau-series mod tau^7: Lemma 1.1 (w_k(0) = (-1)^k phi_0^(k)(1/c) for k = 2..12, w_2 = -delta/(x s*^2),
    and (1.11) against tau^k W_k(U1) from (1.10) at (x, tau) = (2, 5/8)); (3.4) for k <= 12; the code's row in y rewritten by
    homogeneity, wickZ_m(w_k), equal to N mod (tau^7, z^5) with value A(z) = h_w c/(c+z) at tau = 0 (Proposition 3.1(ii),
    (iii)); the uniform version Lnorm[phi_tau o sigma_tau, 1; ytilde1] = N (Remark 3.2).
Suggested Mathematica checks (symbolic in x; not run here, wolframscript being reserved for the main session):
 (i) C3 symbolically, tau^0..tau^8, z^0..z^9, with N from the Wick formula in the coordinate s (not from the stored factors);
 (ii) Proposition 6.1 through z^15;
 (iii) the transport recursion of 7.4 for m <= 5 against [z^m] (Exp[-dphi/z] F) as tau-series to tau^10, including Phi_0 = 1/alg;
 (iv) (4.1), (7.1), P = delta and the closed form of R by Together/PolynomialRemainder modulo x tau s^2 - c s + 1;
 (v) Lemma 1.1 for k <= 20: Together[tau^k W_k(U1)], with W_k = (U d/dU)^k spotU at a0 = 0, a1 = x, a2 = 1 and
     U1 = (x - 1 + delta)/(2x), is analytic at tau = 0 with value (-1)^k phi_0^(k)(1/c) = (-1)^k (k-1)! c^k ((-1)^(k-1) - x^(1-k));
 (vi) V3 from freshly built factors (not the stored factor file): setupWick[12]; rmatrixAt0[{a1 -> x, a2 -> 1, e1 -> x, e2 -> 1}];
     twistRows[hCoeffs[{1 + x, -x, -1}, -1, 12]]; ratCN[0,1][m] evaluated at U1 as a t-series to t^10 against
     (c/(c + z)) F(z, tau) Exp[-dphi/z] alg, symbolically in x, m <= 12.
Changes in this revision (referee report of 2026-10-31, verdict "sound"):
 - The remark after Proposition 3.1 said that the individual Wick monomials in the coordinate y are not obviously bounded at
   t = 0. That was wrong: every monomial has weight 0 under W_k -> lambda^k W_k (Lemma 2.6, new) and therefore equals a monomial
   in the w_k = tau^k W_k(U1), which are analytic at tau = 0 with w_k(0) = (-1)^k phi_0^(k)(1/c) (Lemma 1.1, new; with (1.9),
   (1.10)). The remark is replaced by Remark 3.2. Proposition 3.1 now proves the analyticity directly (part (ii)) and obtains
   r^pt_{0,1}(z; 0) = A(z) without (P1); Remark 3.2 adds the uniform version in ytilde = -Log[1 - tau s]/tau, with no
   continuation step. The proof of Proposition 3.1(i) (fixed t != 0) is unchanged.
 - Theorem 7.2: the clause "coefficientwise in z for every small t" is now proved (end of the proof, via Proposition 3.1(iii),
   and via (7.3) for the larger disc D). Corollary 7.1(b) is added (the tau^k coefficients of F Exp[-dphi/z], without alg, are
   regular at z = 0; this is the form Step 4 uses).
 - (1.2): the phase constant is renamed kappa_ph (the referee's suggested name kappa_1 is already used in S2 s.5.3 for
   c/(kappa_0 (c + z))), and its relation to the constant kappa(z) of outline Step 5 is derived in 7.6, (7.4).
 - Section 0: the labelling of idempotent 1 as U1 is justified (all factors live in one quadratic algebra and are evaluated at
   the same root U1, delta(0) = c, the root with the singular Delta).
 - Citations: DLMF s.5.11(ii) for the remainder in Stirling's series; Hormander's Thm 7.7.5 applies to the Laplace case via
   F = -i (f - f(u0)), omega = 1/z. Lemma 2.5 is proved in full; the citations are cross-references only.
 - 1.4: the re-proof of (1.6) by the recursion of Lemma 5.1 at numerical z now states k + gamma != 0 (Re gamma > 0 under (H_E)).
 - S1_checks.pl: C4 draws until 6 valid points have been tested (before, 2 of the 6 draws were skipped); new check C6.
