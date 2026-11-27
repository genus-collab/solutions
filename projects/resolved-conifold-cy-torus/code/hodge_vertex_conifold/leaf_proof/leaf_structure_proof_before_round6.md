# Leaf structure of the singular column to all orders: proof of Lemma C with F1-F3

Resolved conifold x C^2 at the conifold point Q = 1, fibrewise torus alpha0 = 0, alpha1 = eps1 = x, alpha2 = eps2 = 1, Hodge-vertex frame B.
2026-10-31. This is the assembled reference text of the all-order proof outlined in leaf_proof/outline.md. It is built from the four
section files leaf_proof/S1_thimble.md (Steps 1-2), S2_metric_qde.md (Steps 3, 5, 6 and the second QDE row at the regular idempotent),
S3_leaves_loops.md (Steps 4, 7, 8) and S4_legs_closed.md (Steps 9, 10 and the closed forms). Each section was refereed (verdict
"sound") and revised after its report. Here every result is stated once, in one notation, with complete proofs, in logical order. The
section files remain the working files; Appendix B maps their labels to the labels used here, and Appendix A lists the notation.

No topological recursion and no spectral curves are used anywhere. No step of a proof depends on a computer calculation. The exact
computations listed in section 12 (Mathematica, run before and after the assembly; Perl with Math::BigRat, run for the sections) are
confirmations only. Revised 2026-10-31 after the completeness critic and the Mathematica checks of s.12.2(e) (all passed).

Contents. 0 Summary. 1 Conventions (the code's definitions). 2 Plan of the proof. 3 The singular identity row from the Euler integral.
4 Metric, eigenvalues, QDE and the symplectic constant at the singular idempotent. 5 Dilaton leaves (L1). 6 The H-row in leaf form.
7 The E_11 flow (L3, L2 for j = 1). 8 The regular idempotent. 9 The E_12 and E_21 legs (L2 for j = 2). 10 Proof of the main theorems;
F1-F3. 11 Closed forms. 12 What is proved and what is computer-verified. 13 Consequence: the all-genus regularity theorem is conditional
only on M_irr. 14 Remaining gaps, errata, open points. Appendix A notation; Appendix B concordance with the section files.


## 0. Summary

### 0.1 The objects
In frame B a coloured stable graph of the conifold contributes (conifold_regularity_note.md, preamble)
    (1/|Aut|) Sum_{colourings} Sum_{decorations} prod_v Delta_{c_v}^(g_v - 1) Int Hodge_w psi^a prod Dil_{c_v}(d) prod_e E_{c c'}(a_e, b_e),
with vertex weights Delta_i, dilaton factors Dil_i(d) and edge factors E_ij(m, n), all built from the frame-B R-matrix rows r_{0,i}
(identity row) and r_{1,i} (H-row) of the two idempotents i = 1, 2. Section 1.2 states their definitions exactly as rc_lib.m and
hvc_lib.m compute them. Idempotent 1 is the critical point U_1 of the superpotential that runs into the logarithmic singularity at the
conifold point t = log Q = 0 (the singular idempotent); idempotent 2 stays away from it (the regular idempotent). Put c := 1 + x,
    rho(t) := (Delta_1 t^2/(-x c^3))^(-1/2) (rho(0) = 1),   lam_j(psi) := 1/(c - j psi),   Lam_m := span(lam_1, ..., lam_m),
q_t(psi) := Sum_{d>=1} Dil_1(d; t) psi^d (the dilaton series), f_{j,b}(psi) := [w^b] E_1j(psi, w) (the legs) and
q_0 := c - c^2 lam_1 = q_t|_{t=0}.

### 0.2 The main theorem
**Theorem LC (Lemma C with F1-F3, to all orders in t and in z).** Let x be a complex number with x != 0, -1, or an indeterminate,
and c = 1 + x. All factors are the frame-B factors of section 1.2, at their true idempotent (s.1.2 (C3)). Then:
 (F1) Delta_1 = -x c^3/(t^2 rho^2), where rho = tau/(t alg) (tau = 1 - e^-t, alg as in s.1.5) is analytic at t = 0 with
      rho(0) = 1 (equivalently, t^2 Delta_1/(-x c^3) is a unit of C{t} with value 1); so Delta_1 = -(x c^3/t^2)(1 + (1-x)^2 t/c^2 + O(t^2)).
      Delta_2 = x^3 delta s*^2 is analytic at t = 0, with Delta_2(0) = x^3/c != 0 and Delta_2'(0) = 0.
 (F2) r_{0,1}(z; 0) = c/(c + z), and rho r_{1,1}(z; t) = c/t + x/(c + z) + O(t).
 (F3) Every frame-B dilaton factor Dil_i(d) (i = 1, 2, d >= 1) and every edge factor E_ij(m, n) (i, j in {1, 2}, m, n >= 0) is
      analytic at t = 0.
 (L1) [t^0](rho q_t) = q_0 = c - c^2 lam_1, and [t^k](rho q_t) lies in Lam_{k+1} for every k >= 1.
 (L2) For j in {1, 2}, every b >= 0 and every k >= 0, [t^k](rho f_{j,b}) lies in Lam_{k+1}. Moreover E_21(b, a) = E_12(a, b), so the
      leg seen by a singular vertex is f_{2,b} for either orientation of a mixed edge; and [t^k] f_{j,b} lies in Lam_{k+1} as well.
 (L3) For every k >= 0, [t^k](rho^2 E_11(psi, psi')) lies in span{lam_i(psi) lam_l(psi') : i, l >= 1, (i-1) + (l-1) <= k}, and
      E_11(a, b) = E_11(b, a). The same support bound holds for [t^k](rho^e E_11) for every integer e.
Every membership is an identity of formal power series in psi (and psi'), valid at every psi-order (that is, to all orders in z) and
every t-order; the coefficients are rational functions of x. The bounds are sharp for generic x: if c is not in {(k+1)/i : 1 <= i <= k},
the top leaf lam_{k+1} occurs in [t^k](rho q_t); if c is not in {(m+1)/i : 1 <= i <= m <= k}, every entry of the anti-diagonal
(i-1) + (l-1) = k of L3 is non-zero (Corollaries 11.3, 11.5: by (7.6) the entry (k1+1, k2+1) is a non-zero multiple of x P_{k1} P_{k2}, and
P_m vanishes at c = (m+1)/i for every m <= k, e.g. e^(2)_22 = 0 at x = 1). Equivalent large-z form: coefficientwise in t, z r_{0,1}(z; t) -> c/rho(t) and E_1j(z, w) -> 0 as z -> oo at fixed order
in w, in the sense of the rational functions of z of Lemma 5.2 and Corollary 9.5.

### 0.3 Further results proved on the way
**Theorem Q (quantum differential equation and symplecticity).** (i) For i = 1, 2 the second row of the quantum differential equation
holds to all orders in z: z D_i v_i = (A - L_i) v_i, with v_i = (r_{0,i}, r_{1,i}), D_i = d/dt - Delta_i'/(2 Delta_i) and
A = [[0, 1], [-x/tau, c/tau]] (Proposition 4.6, Theorem 8.4). (ii) <v_i(w), v_j(-w)> = delta_{ij} Delta_i for all i, j, to all orders in w
(Proposition 4.8, Corollaries 8.5, 8.7). (iii) v_2(z; 0) = h_w(z)^-2 (1, x/c) (Proposition 8.6).

**Theorem CF (closed forms; section 11).** In frame B,
    rho r_{0,1}(z; t) = (c/(c + z)) (tau/t) 2F1(1/z + 1, 1/z + 1; c/z + 2; tau) Exp[-dphi(tau)/z]                     (KEY FORMULA),
    rho q_t(psi) = c + Sum_{j>=1} A_j(t) lam_j(psi),   A_j(t) = (-1)^j C_j t^-1 zeta(t)^j 2F1(j x/c, j x/c; j; tau),
    d/dt [t^2 rho^2 E_11(psi, psi')] = Theta(t) d/dt(t g(psi)) d/dt(t g(psi')),   Theta = -x^3 tau s*^2,   g(psi) = -(1/c) Sum_j j A_j lam_j,
with C_j, zeta, s*, dphi explicit (s.1.5); hence closed forms of every leaf and loop coefficient, of the top and subleading coefficients
for all k, and of the t = 0 values of E_12 and E_22.

### 0.4 Consequence
**Corollary R (section 13).** Fix g >= 2 and assume M_irr(h) for 2 <= h <= g (proof_status_mixed_scale.md s.1.2). Then every coloured
stable graph of genus g with at least one edge is regular at t = 0 after summing over decorations, and the single-vertex graph has the
single pole gap_g t^(2-2g). M_irr(2) is proved, and M_irr(3) and M_irr(4) are exact computations (proof_status_mixed_scale.md, addendum
(1) and s.1.9). So the theorem holds through genus 4 (computer-assisted, independently of the direct graph sums), and in every genus it
is conditional only on M_irr(h), h >= 5, i.e. on Lemma A of proof_status_mixed_scale.md s.2. Lemmas B (P2) and D, the
topological-recursion route to Lemma A, are not needed.

### 0.5 Inputs beyond algebra
One analytic input enters Theorem LC: the value Abeta(z) = h_w(z) c/(c + z) of the normalised Laplace expansion of the Beta integral at
t = 0 (Proposition 3.15). It is proved from Laplace's method (Lemma 3.8, proved here in full) and Stirling's series with remainder for
real positive argument (DLMF 5.11.1 and s.5.11(ii), cited), first for real x > 0 and then for all x by rationality. A second analytic
input of the same kind, r^pt_{0,2}(z; 0) = h_w(z)^-1 (Proposition 8.6(a)), is used only for K_2 = 1 and for the t = 0 values at the
regular idempotent, not for Theorem LC or Corollary R. Everything else is formal: identities of formal power series, plus the
analyticity in tau of finitely many explicit algebraic functions. The regularity at z = 0 of the tau-coefficients of the singular row
has, in addition, a purely formal proof (section 3.10) that uses neither the Laplace identification nor the value of Abeta.

## 1. Conventions

### 1.1 Parameters and variables
- x is a complex number with x != 0, -1, or an indeterminate; c := 1 + x (so x != 0 and c != 0 throughout).
- Q = e^t and tau := 1 - 1/Q = 1 - e^-t. Hence Q = 1/(1 - tau), Q tau = Q - 1, dtau/dt = 1 - tau, tau = t - t^2/2 + t^3/6 + O(t^4) and
  t = -log(1 - tau).
- Derivatives. On functions of t (at fixed auxiliary variables z, w, psi, psi') f' := df/dt, and f_tau := df/dtau, so f' = (1 - tau) f_tau.
  On functions written with the Euler variable s (phases phi_tau, amplitudes) a prime or a superscript (k) means d/ds, and d_tau is the
  partial derivative in tau at fixed s. Dz := z d/dt = z (1 - tau) d/dtau (the operator D of outline V7).
- [t^k], [z^m], [psi^a], [w^b] are coefficients, taken coefficientwise in the other variables; X^(k) := [t^k] X (k may be negative).

### 1.2 The factors, as rc_lib.m and hvc_lib.m compute them
(rc_lib.m s.1, s.3 and hvc_lib.m s.5, at alpha0 = 0, alpha1 = e1 = x, alpha2 = e2 = 1, mode B, s = -1.)
 (C1) Superpotential (spotU): W(U) := x log(1 - Q U) - log U + log(1 - U); theta := U d/dU; W_k := theta^k W for k >= 1. The W_k are
      single valued, rational in U and Q.
 (C2) Critical points: U W_U = 0 iff LeqU := Q x U^2 + Q (1 - x) U - 1 = 0 (Lemma 4.2). Its roots satisfy U_1 + U_2 = (x - 1)/x (SsumQ),
      U_1 U_2 = -1/(Q x) (PprodQ), and its discriminant is ((x - 1)^2 + 4x/Q)/x^2 = (c^2 - 4 x tau)/x^2. U_1 -> 1 and U_2 -> -1/x as t -> 0
      (explicit forms in s.1.3).
 (C3) Quadratic algebra and conjugate encoding. Every factor is computed in QA := Q(x, Q)[U]/(LeqU), with elements C[a, b] = a + b U.
      ev_i : C[a, b] -> a + b U_i are ring homomorphisms into algebraic functions of Q, jointly injective (LeqU has distinct roots), and
      cconj[C[a, b]] := C[a + b Ssum, -b] is U -> Ssum - U, so ev_1 o cconj = ev_2 and ev_2 o cconj = ev_1. rc_lib.m stores every quantity
      of idempotent 2 as cconj of the corresponding quantity of idempotent 1 (ratCN[0,2] = cconj[wickCN], LCN[2] = cconj[LCN[1]],
      DeltaCN[2] = cconj[DeltaCN[1]]), and all factors, of either idempotent, are evaluated at U_1 (the checks take
      sq = (c/x)(1 - 4 x tau/c^2)^(1/2) as a t-series and U1 = (S0 + sq)/2, hvc_leaf_proof_check.wls line 22). In this text the value of
      a factor of idempotent 2 always means its true value, ev_1 of the stored element. (Remark 8.8 shows what evaluation at U_2 gives.)
 (C4) Flat basis {1, H}. Frobenius data (Lemma 4.1): eta = {{c/x^3, 1/x^2}, {1/x^2, 0}}, eta^-1 = {{0, x^2}, {x^2, -c x}}, both constant in
      t. Covectors a = (a_0, a_1) (components paired with 1 and H) are paired by
          <a, b> := Sum eta^{a'b'} a_{a'} b_{b'} = x^2 (a_0 b_1 + a_1 b_0) - c x a_1 b_1        (etainvSym; the pairing of evalQ).
      It is symmetric with constant coefficients, so d/dt <a, b> = <a', b> + <a, b'>.
 (C5) Eigenvalue of H*: L(U) := -Q (x U + 1)/(1 - Q) = (x U + 1)/tau (LofU); L_i := L(U_i).
 (C6) Delta_i := 1/eta(e_i, e_i) (DeltaCN[i]; e_i the idempotents, Proposition 4.3). The vertex weight is Delta_i^(g_v - 1).
 (C7) Point-frame identity row (rmatrixAt0 calls rmatrixCore with ffc1 = ffc2 = (1, 0, 0, ...): no QRR constants at alpha0 = 0):
          r^pt_{0,i}(z; t) := Sum_{m>=0} z^m wickZ_m(W_2(U_i), W_3(U_i), ...),
          wickZ_m := Sum_{n=m}^{3m} ((2n)! (-1)^n/(n! 2^n W_2^n)) Sum_p prod_k W_k^(m_k)/((k!)^(m_k) m_k!)        (wickDirect),
      p running over the partitions of 2n into exactly n - m parts >= 3, m_k the multiplicity of k in p. Equivalently wickZ_m is the
      coefficient of z^m of Sum_{n>=0} ((2n)!/(n! 2^n)) (-z/W_2)^n [v^(2n)] Exp[Sum_{k>=3} W_k v^k/(k! z)]: a term with j factors W_k
      and v-degree 2n has z-order n - j. The cut-off "parts <= 2 ZMAX + 2" of wickDirect never binds: at order z^m with j = n - m >= 1
      factors the largest part is 2n - 3(j - 1) = 3m + 3 - n <= 2m + 2. Every monomial W_2^-n prod_i W_{k_i} of wickZ_m has
      Sum_i k_i = 2n (Lemma 3.9).
 (C8) Frame B (buildFactors[zmax, "B", -1, x], i.e. twistRows[hCoeffs[{c, -x, -1}, -1, zmax]]): every row is multiplied by h_w(z)^-1,
          r_{a,i}(z; t) := h_w(z)^-1 r^pt_{a,i}(z; t),    h_w(z) := Exp[Sum_{k>=1} bk_k St_k z^(2k-1)],
          bk_k := B_2k/(2k (2k - 1)),   St_k := 1 + x^(1-2k) - c^(1-2k).
      Indeed hCoeffs builds ex := -Sum_k bk_k z^(2k-1) Sum_a w_a^(1-2k) with w = (c, -x, -1), Sum_a w_a^(1-2k) = c^(1-2k) - x^(1-2k) - 1 = -St_k,
      and multiplies by Exp[s ex] = Exp[-ex] = h_w^-1 (s = -1). h_w is t-independent, lies in 1 + z Q(x)[[z]], and h_w(z) h_w(-z) = 1
      (its exponent is odd in z). By Proposition 3.15, h_w is the Stirling tail of the Beta function B(x/z, 1/z).
 (C9) H-row (rmatrixCore; cDerN is the total Q-derivative, including dU_i/dQ through dUC, and Q d/dQ = d/dt):
          (H)   r_{1,i}(z) = (L_i + z D_i) r_{0,i}(z),     D_i := d/dt - Delta_i'/(2 Delta_i),
      i.e. r_{1,i}[0] = L_i r_{0,i}[0] and r_{1,i}[m] = L_i r_{0,i}[m] + D_i r_{0,i}[m-1] for m >= 1. (H) is linear over C[[z]] with
      t-independent coefficients, so it holds in frame B and in the point frame. We write v_i(z) := (r_{0,i}(z), r_{1,i}(z)) (a covector)
      and v_i[p] := [z^p] v_i.
 (C10) Dilaton factors (DilN): Dil_i(n) := (-1)^n r_{0,i}[n-1] (n >= 1). Hence
          q_t(psi) := Sum_{d>=1} Dil_1(d; t) psi^d = Sum_{m>=0} (-1)^(m+1) r_{0,1}[m] psi^(m+1) = -psi r_{0,1}(-psi; t).
      The code's decorations carry dilaton leaves of degree d >= 2 only (graphs.m line 54, hvc_lib.m line 124: Min[#] >= 2). The term d = 1 is
      the dilaton shift: Dil_1(1) psi = -r_{0,1}[0] psi = -psi at every t, since r_{0,1}[0] = 1. So graphValueC computes
      Sum_m (1/m!) < ... T_t^m > with T_t := q_t + psi, which is Fpot_h(q_t) in dilaton-shifted coordinates (s.13.1 (i)); this identification is
      what makes the homogeneity step of Theorem R (rescaling the shift -psi to -rho psi) match the code. The same holds at idempotent 2.
 (C11) Edge factors (evalQ): pre[[p+1, q+1]] = Sum_{a,b} -(-1)^(p+q) eta^{ab} rN[a,i][p] rN[b,j][q] and
          EdN[i,j][m, n] = Sum_{k=0}^{n} (-1)^k P^ij_{m+1+k, n-k},     P^ij_pq := -(-1)^(p+q) <v_i[p], v_j[q]>.
 (C12) Legs and leaf variables. f_{j,b}(psi) := Sum_a E_1j(a, b) psi^a = [w^b] E_1j(psi, w). In leaf language E_11(psi, psi') means
      E_11(z = psi, w = psi'), without change of sign, as in hvc_leaf_structure.wls and outline V8.
 (C13) rho := (Delta_1 t^2/(-x c^3))^(-1/2), the branch with rho(0) = 1 (as in leaf_structure.md and hvc_leaf_proof_check.wls, check (2c);
      rho occurs in neither library). Leaves lam_j(psi) := 1/(c - j psi) (j >= 0; lam_0 = 1/c); in the R-matrix variable z = -psi,
      lam_j = 1/(c + j z).
 (C14) Scope. The objects of the theorems are the formal series in z (and w) defined by (C7)-(C11) for all m, n. For every ZMAX the code
      computes exactly their truncation.

**Lemma 1.1 (the edge factors as generating functions).** For any two rows (no symplecticity is assumed), with
P_ij(z, w) := Sum_{p,q} P^ij_pq z^p w^q = -<v_i(-z), v_j(-w)>,
    E_ij(z, w) := Sum_{m,n} E_ij(m, n) z^m w^n = [P_ij(z, w) - P_ij(-w, w)]/(z + w).                                            (1.1)
The numerator vanishes on z = -w, so the quotient lies in C[[z, w]], coefficientwise in t.
Proof. For p >= 1, (z^p - (-w)^p)/(z + w) = Sum_{l=0}^{p-1} z^l (-w)^(p-1-l); for p = 0 the quotient is 0. So the right side of (1.1) is
Sum_{p,q} P^ij_pq w^q Sum_{l=0}^{p-1} z^l (-w)^(p-1-l). Its coefficient of z^m w^n collects l = m, p = m + 1 + k and q = n - k with
k = p - 1 - m in [0, n], with sign (-1)^k: this is EdN. Only entries with p + q = m + n + 1 <= ZMAX occur, so (1.1) is exact for every
entry the code computes. QED

### 1.3 Algebraic functions of tau
- delta := c (1 - 4 x tau/c^2)^(1/2), with the principal root of the bracket. It is analytic on the disc D := {|tau| < |c|^2/(4|x|)},
  delta/c has positive real part there, delta(0) = c, and delta^2 = c^2 - 4 x tau. (For c < 0 this is not the principal root of
  c^2 - 4 x tau.)
- s* := 2/(c + delta) = (c - delta)/(2 x tau) = (1/c) Sum_n Cat_n (x tau/c^2)^n, analytic on D, s*(0) = 1/c ((c - delta)(c + delta) = 4 x tau).
- s~ := (c + delta)/(2 x tau) = 1/(x tau s*); s* and s~ are the two roots of x tau s^2 - c s + 1 (sum c/(x tau), product 1/(x tau)).
  tau s~ = (c + delta)/(2x) is analytic with value c/x at 0.
- U_1 = 1 - tau s* = (x - 1 + delta)/(2x) -> 1 and U_2 = 1 - tau s~ = 1 - 1/(x s*) = (x - 1 - delta)/(2x) -> -1/x (Lemma 4.2(c)). U_1 is the
  root (S0 + sq)/2 of the scripts, whatever the sign of c, and it is the root with the singular Delta (Proposition 4.3(d)).
- Discs: D as above; D' := D intersected with {|tau| < 1}; D'' := D' - {0} = {0 < |tau| < r0}, r0 := min(1, |c|^2/(4|x|)). For real x > 0,
  (1 + x)^2/(4x) >= 1 and D' = {|tau| < 1}.

### 1.4 Rings and general principles
- O := C{t} = C{tau} (convergent power series; tau(t) is biholomorphic near 0 and tau/t, t/tau are units of O); M := C({t}) = C({tau}),
  the field of meromorphic germs; O(D), O(D') the functions analytic on D, D'.
- Kz := C((z)) (formal Laurent series); C(z) is embedded in Kz by Laurent expansion at z = 0. Rz := Kz((tau)), with the Kz-linear
  derivation d/dt = (1 - tau) d/dtau; Rz0 := C[[z]]((tau)) = C[[z, tau]][1/tau], a subring of Rz and (expanding in z) of C((tau))[[z]],
  compatibly with products and d/dt. iota : O[[z]] -> C[[z, tau]], "Taylor-expand every z-coefficient", is an injective ring homomorphism
  commuting with d/dtau and with multiplication by elements of O. refl : z -> -z (coefficientwise) commutes with d/dt; f(-z) := refl(f).
  M((z)) denotes formal Laurent series in z, bounded below, with coefficients in M (d/dt coefficientwise).
- TRANSFER. An identity in Rz between elements of Rz0 is an identity of formal double series: reading off z^m gives identities of formal
  Laurent series in tau, and, when these are expansions of meromorphic germs (as for all rows), identities of meromorphic functions of t
  near 0. "To all orders in z" means this; in M((z)) it means "coefficientwise in z, as identities in M".
- (P) Passage between complex x and x indeterminate. Kx := Q(x). The objects of this text (tau, s*, delta, dphi, alg, rho, L_i, Delta_i,
  Theta, the rows r_{a,i}[m], h_w, and the edge and dilaton factors) are built by field operations, d/dt and the expansions of
  (1 - 4 x tau/c^2)^(1/2) and of logarithms of units from x and tau. They are analytic in (x, tau) near every (x0, 0) with x0 != 0, -1
  (the only denominators at tau = 0 are powers of x and c: U_1(0) = 1, U_2(0) = -1/x, 1 - U_2(0) = c/x, delta(0) = c), and their
  t-, z- and psi-coefficients are rational functions of x with rational coefficients, regular at every x0 != 0, -1. The F_n, G_k and
  H_n lie in Kx(z). Hence an identity between such coefficients proved for all complex x in a set with an accumulation point holds in
  Q(x), and conversely an identity or a bound proved over Q(x) specialises to every x != 0, -1. Sections 5-7 are written over Kx and
  apply verbatim over C for every fixed x != 0, -1. Only the sharpness statements need x generic (prod_{m<=k} P_m(x) != 0).
- Leaf spaces. Lam_m := span(lam_1, ..., lam_m) has dimension m (the poles c/j are distinct); Lam_m^+ := (constants) + Lam_m. A rational
  function f of psi lies in Lam_m^+ iff its poles are simple, lie in {c/j : 1 <= j <= m}, and f is bounded at oo; it lies in Lam_m iff
  moreover f(psi) -> 0 as psi -> oo. A rational function regular at 0 is identified with its Taylor series at 0 (an injective linear map,
  compatible with products and with the substitution z -> -psi); thus lam_j(psi) = Sum_{n>=0} j^n psi^n/c^(n+1).
- Elementary leaf identities (j, m >= 1):
    (LAM1) [psi^n] lam_j = j^n/c^(n+1);
    (LAM2) psi lam_j = (c lam_j - 1)/j, i.e. -psi lam_j = (1 - c lam_j)/j;
    (LAM3) (c - m psi) lam_j = m/j + c (1 - m/j) lam_j; in particular (c - j psi) lam_j = 1;
    (LAM4) for f = Sum_{j=1}^{M} f_j lam_j one has f = -(1/psi) Sum_j f_j/j + O(psi^-2) at psi = oo, so f = O(psi^-2) iff Sum_j f_j/j = 0;
    (LAM5) lam_i lam_j = (i lam_i - j lam_j)/(c (i - j)) for i != j;
    (LAM6) lam_j(psi) - lam_j(-w) = j (psi + w) lam_j(psi) lam_j(-w).
  Proofs: (c - j psi) lam_j = 1 gives (LAM2); (c - m psi) lam_j = 1 + (j - m) psi lam_j = 1 + ((j - m)/j)(c lam_j - 1) gives (LAM3);
  lam_j = -1/(j psi) - c/(j^2 psi^2) + ... gives (LAM4); (LAM5) and (LAM6) are identities of partial fractions.
- (PF) Partial fractions. Let S be a finite set of points and f a rational function whose poles all lie in S, of order <= 1, with f
  bounded at oo. Then f = f(oo) + Sum_{p in S} Res_{z=p}(f)/(z - p). [The difference h has no pole and h(oo) = 0; a rational function
  without poles is a polynomial, and a polynomial with limit 0 is 0.] For S = {-c/j : 1 <= j <= M} (distinct and non-zero since c != 0):
      f = f(oo) + Sum_{j=1}^{M} beta_j/(c + j z),     beta_j = j Res_{z=-c/j} f = lim_{z -> -c/j} (c + j z) f(z).
- (T) tau-orders versus t-orders. For a vector space V and f = Sum_n f_n tau^n in V[[tau]], V[[tau]] = V[[t]] and
      [t^k] f = Sum_{n=0}^{k} u_{k,n} f_n,    u_{k,n} := [t^k] tau^n in Q,   u_{k,k} = 1,   u_{k,n} = 0 for n > k,
  with the inverse unitriangular relation from t = -log(1 - tau). So for every increasing chain V_0 c V_1 c ..., ([tau^n] f in V_n for
  n <= k) iff ([t^n] f in V_n for n <= k). Multiplication by a series in t that does not involve the auxiliary variables also preserves
  such a filtration: [t^k](a f) = Sum_m a_m [t^(k-m)] f.

### 1.5 The closed-form ingredients
- Euler phase: phi_tau(s) := log s + x log(1 - s) - log(1 - tau s). phi_0(1/c) = x log x - c log c and phi_0''(1/c) = -c^3/x (Lemma 3.1).
- dphi(tau) := phi_tau(s*(tau)) - phi_0(1/c), principal branches near (tau, s) = (0, 1/c); dphi lies in tau O, with
  dphi_tau = s*/(1 - tau s*) = (1 - x s*)/(1 - tau), dphi' = 1 - x s* (Lemma 3.1(b)), and
      dphi = tau/c + (1 + 2x) tau^2/(2 c^3) + O(tau^3) = t/c - x^2 t^2/(2 c^3) + (1 - x)^2 x^2 t^3/(6 c^5) + O(t^4).
- alg(tau) := (1 - tau s*) (phi_tau''(s*)/phi_0''(1/c))^(1/2), the branch with alg(0) = 1; alg = (delta/c)^(1/2)/(c s*)
  = (delta/c)^(1/2)(c + delta)/(2c) (Lemma 3.1(c)).
- Gauss's hypergeometric series. alpha = beta := 1/z + 1, gamma := c/z + 2, and
      F(z, tau) := 2F1(alpha, beta; gamma; tau) = Sum_{n>=0} F_n(z) tau^n,   F_n = (alpha)_n (beta)_n/((gamma)_n n!)
               = prod_{j=1}^{n} (1 + j z)^2 / ( n! z^n prod_{j=2}^{n+1} (c + j z) )   in Kx(z),
  using (1/z + 1)_n = z^-n prod_{j=1}^{n} (1 + j z) and (c/z + 2)_n = z^-n prod_{j=2}^{n+1} (c + j z). Gauss's operator:
      Gauss[Y] := tau (1 - tau) Y_tautau + (gamma - (alpha + beta + 1) tau) Y_tau - alpha beta Y.
- Ephi := Exp[-dphi/z] := Sum_l (-dphi)^l/(l! z^l), an element of C[1/z][[tau]] (tau-adically convergent since dphi = O(tau)); its tau^n
  coefficient X_n is a polynomial of degree <= n in 1/z. Ephi^-1 = refl(Ephi) = Exp[dphi/z] and z Ephi' = -dphi' Ephi.
- H_n := [tau^n](F Ephi) = Sum_{n1+n2=n} F_{n1} X_{n2} and Psi_n := [tau^n](F Ephi alg), elements of Kx(z).
- G(z; t) := (c/(c + z)) (tau/t) F(z, tau) Ephi(z, tau) in Kx(z)[[t]] (Kx(z)[[tau]] = Kx(z)[[t]] by (T)); G_k := [t^k] G. By the KEY FORMULA
  (Proposition 4.4) G = rho r_{0,1}. g(psi; t) := (rho r_{0,1})(-psi; t), g_k := [t^k] g, Dg := d/dt (t g) = Sum_k (k + 1) g_k t^k.
- QDE data: A := [[0, 1], [-x/tau, c/tau]] (acting on covectors (a_0, a_1) -> (a_1, -(x/tau) a_0 + (c/tau) a_1)); chi_i := (1, L_i);
  pi(a) := <a, chi_2>; mu := c/tau - 1 = x + c (1 - tau)/tau; omega := (c + z)/tau - (1 + z).
- Theta := (L_2 - L_1)(c x L_2 - x^2)^2/Delta_2 (the definition in V8; Theta = -x^3 tau s*^2 by Theorem 7.2).
- sigma(psi) := rho v_1(-psi) = (sigma_0(psi), sigma_1(psi)) (two components) and sigma^(k) := [t^k] sigma, k >= -1.
- Leaf coefficients: q^(k) := [t^k](rho q_t) (written delta_k in leaf_structure.md and proof_status_mixed_scale.md); a_j^(k) := the
  lam_j-coefficient of q^(k) (a_1^(0) := -c^2; a_j^(k) := 0 for j > k + 1); g_j^(k) := -(j/c) a_j^(k), the lam_j-coefficient of g_k;
  A_j(t) := Sum_k a_j^(k) t^k and Ahat_j := t A_j; e^(k)_il := the lam_i(psi) lam_l(psi')-coefficient of [t^k](rho^2 E_11) and
  e_il(t) := Sum_k e^(k)_il t^k.
- P_k := prod_{i=1}^{k} (k + 1 - i c)^2 (P_0 = 1); C_j := P_{j-1}/(j! (j-1)! c^(2j-4)) (C_1 = c^2); zeta(t) := tau e^(dphi(tau)/c).

## 2. Plan of the proof and logical order

The outline's ten steps become the following chain of statements; an arrow means "is used in the proof of". Forward references occur
only for explicit low-order values and for sharpness (Proposition 4.4(d), the k = 1 values in Proposition 6.1, the anti-diagonal (7.6)),
which use section 11; section 11 uses only statements that do not depend on these, so there is no circularity.

 Step 1  Euler coordinate, critical points, Hessians, rescaled derivatives ........................... s.3.1-3.3 (Lemmas 3.1, 3.2)
 Step 2  formal Laplace calculus (P1, P2, P3; homogeneity; Laplace's method) ........................ s.3.4-3.5 (Lemma 3.3 - Lemma 3.9)
         -> the code's Wick row equals the normalised expansion N of the Euler integrand (Prop. 3.10)
         -> Gauss's equation for the formal expansion (Prop. 3.12); uniqueness (Lemma 3.13)
         -> N = Abeta F Ephi alg (Thm 3.14); Abeta = h_w c/(c + z) [the analytic input] (Prop. 3.15)
         -> regularity at z = 0 of every tau^k-coefficient (Cor. 3.16; purely formal second proof: Prop. 3.19)
         -> frame-B row r_{0,1} = (c/(c + z)) F Ephi alg (Thm 3.17); F2 (Cor. 3.18)
 Step 3  Frobenius data, eigenvalues, metric; Delta_1 = x^2 W_2(U_1) (V1) -> rho = tau/(t alg), F1,
         KEY FORMULA G = rho r_{0,1} (V3) .............................................................. s.4.1-4.4 (Lemmas 4.1, 4.2, Props 4.3, 4.4)
 Step 5  operator identity V7 (Lemma 4.5) + Gauss + KEY FORMULA -> second QDE row at U_1 (Prop. 4.6); the operator A (Prop. 4.7)
 Step 6  K_1 is t-independent; its t -> 0 limit from F2 and the H-row -> K_1 = 1 (Prop. 4.8)
 Step 4  partial fractions of [t^k](z G) -> L1 (Thm 5.4) ........................................................ section 5
 Step 7  the H-row in leaf form: sigma = t^-1 (0, c) + Sum_k t^k sigma^(k), sigma^(k) in Lam_{k+1}^2 (Prop. 6.1); pi(sigma) (Prop. 6.2)
 Step 8  the E_11 flow (Thm 7.2) from Steps 5-7; integration with K_1 = 1 -> F3 for E_11, L3 (Thm 7.3), L2 for j = 1 (Cor. 7.4)
 Step 10 regular idempotent: analyticity of all its rows, Delta_2, Dil_2, E_22 (Lemmas 8.1, 8.2); the second QDE row at U_2 by the
         formal Laplace calculus at s~ (Thm 8.4) -> orthogonality (Cor. 8.5), v_2(z; 0) proportional to (1, x/c) (Prop. 8.6(b)),
         K_2 = 1 (Cor. 8.7, with the analytic input Prop. 8.6(a))
 Step 9  rationality argument -> L2 for j = 2 and regularity of E_12 (Prop. 9.2, using only the analyticity of the regular rows);
         orthogonality -> E_21 = E_12^T and E_21 regular (Prop. 9.4)
 Assembly: Theorem LC and Theorem Q (section 10); closed forms (section 11); Theorem R (section 13).

Dependencies on the code. Only the definitions (C1)-(C14) are used. In particular nothing about the code's numerical output is used.

Two remarks on the logic. (1) The only place where the frame matters is Proposition 3.15 combined with (C8): the twist h_w^-1 cancels
exactly the Stirling tail h_w of the singular row and leaves the geometric series c/(c + z) (F2); every later statement that involves
the value c/(c + z) (L1 at t^0, [t^-1] sigma = (0, c), K_1 = 1, F3 for E_11, the rank-one t^0 loop) rests on it. (2) The second QDE row
at the regular idempotent is not needed for L1, L3, for L2 with j = 1, 2 as legs f_{j,b} = [w^b] E_1j(psi, w), or for the regularity of
E_12; it is needed for the regularity of E_21 and for the transposition E_21(b, a) = E_12(a, b), which the graph sum uses because
graphValueC attaches EdN[c(v1), c(v2)][a_{v1}, a_{v2}] in the order of the edge list.

## 3. The singular identity row from the Euler integral (Steps 1-2; F2)

Statements are marked [formal] (identities of formal series, valid for every x with x != 0, c != 0) or [analytic] (statements about
convergent integrals, with explicit hypotheses).

### 3.1 The Euler coordinate
For tau != 0 put U = 1 - tau s. Using Q tau = Q - 1,
    1 - U = tau s,    1 - Q U = 1 - Q + Q tau s = (1 - Q)(1 - s) = -Q tau (1 - s),    dU/U = -tau ds/(1 - tau s).
Hence, on every simply connected domain of the s-plane avoiding s = 0, 1, 1/tau, with compatible branches of the logarithms,
    W(1 - tau s) = x log(-Q tau (1 - s)) + log(tau s) - log(1 - tau s) = phi_tau(s) + c log tau + x log Q + x log(-1),           (3.1)
    e^{W/z} dU/U = kappa_ph(z) Q^(x/z) tau^(c/z + 1) e^{phi_tau(s)/z} (1 - tau s)^-1 ds,      kappa_ph(z) := -e^{x log(-1)/z}.       (3.2)
kappa_ph depends on z and on the branches only, not on s or t. It is a constant of the integrand; the t-independent constant of the
outline's Step 5 is a different one (Remark 4.10). In Pochhammer form e^{phi_tau(s)/z} (1 - tau s)^-1 = s^(1/z) (1 - s)^(x/z)
(1 - tau s)^(-1/z - 1). The logarithmic points U = 1, U = 1/Q, U = 0 of W are s = 0, s = 1, s = 1/tau. At t = 0 the first two collide
(1 = 1/Q) and U_1 runs into them; in the coordinate s they stay at 0 and 1, the critical point stays between them (s*(0) = 1/c), and
U = 0 escapes to s = oo.

### 3.2 Critical points, Hessians, the tensor rule
**Lemma 3.1.** Let tau != 0.
 (a) phi_tau'(s) = 1/s - x/(1 - s) + tau/(1 - tau s) = (x tau s^2 - c s + 1)/(s (1 - s)(1 - tau s)). So the critical points of phi_tau in
     C - {0, 1, 1/tau} are s* and s~, and 1 - tau s* = U_1, 1 - tau s~ = U_2; in particular s* = (1 - U_1)/tau (V2). On D', s* is not in
     {0, 1, 1/tau}; on D'', s~ is not in {0, 1, 1/tau}, and U_2 is not in {0, 1, 1/Q}.
 (b) At s = s* (and, for (ii) and (iii), at s = s~):
     (i) c - 2 x tau s* = delta;   (ii) 1 - s = x s (1 - tau s);   (iii) s/(1 - tau s) = (1 - x s)/(1 - tau);
     (iv) ds*/dtau = x s*^2/delta;   (v) 1 - c s* = -x tau s*^2.
     Consequently dphi_tau = s*/(1 - tau s*) = (1 - x s*)/(1 - tau) and dphi' = 1 - x s* (envelope formula).
 (c) Hessians. phi_tau''(s) = -1/s^2 - x/(1 - s)^2 + tau^2/(1 - tau s)^2, and at both critical points
        (1 - tau s)^2 phi_tau''(s) = (2 tau s - c/x)/s^2,                                                                    (3.3)
     so (1 - tau s*)^2 phi_tau''(s*) = -delta/(x s*^2) and (1 - tau s~)^2 phi_tau''(s~) = delta/(x s~^2); phi_0''(1/c) = -c^3/x. Hence both
     critical points are nondegenerate on D'', and
        alg^2 = (1 - tau s*)^2 phi_tau''(s*)/phi_0''(1/c) = delta/(c^3 s*^2),     alg = (delta/c)^(1/2)/(c s*) = (delta/c)^(1/2)(c + delta)/(2c),   (3.4)
     a unit of O with alg(0) = 1.
 (d) Tensor rule. For a function G of U and the composite G(1 - tau s),
        d/ds G(1 - tau s) = -(tau/(1 - tau s)) (theta G)(1 - tau s),   d^2/ds^2 G(1 - tau s) = (tau/(1 - tau s))^2 [theta^2 G - theta G](1 - tau s).
     Since phi_tau = W(1 - tau s) + const (3.1) and W_1(U_i) = 0,
        W_2(U_i) = U_i^2 phi_tau''(s_i)/tau^2     ((s_1, s_2) := (s*, s~)),
        W_2(U_1) = -delta/(x tau^2 s*^2),    W_2(U_2) = delta/(x tau^2 s~^2) = x delta s*^2.                                     (3.5)
     W_2(U_1) has a double pole at tau = 0: the coordinate y = log U degenerates there (d/ds log(1 - tau s) = O(tau)), the coordinate s
     does not.
Proof. (a) (1 - s)(1 - tau s) - x s (1 - tau s) + tau s (1 - s) = 1 - c s + x tau s^2. Under U = 1 - tau s, by s.1.3, 1 - tau s* = U_1 and
1 - tau s~ = U_2. If s* = 1 or s* = 1/tau, the critical equation forces tau = 1, which is excluded on D'; s* = 0 contradicts the constant
term 1. The same holds for s~ on D''. U_2 = 0 would force delta^2 = (x - 1)^2, i.e. 4x (1 - tau) = 0 (since c^2 - (x - 1)^2 = 4x); and
U_2 in {1, 1/Q} means s~ in {0, 1}, by 1 - Q U_2 = -Q tau (1 - s~).
(b) (i): 2 x tau s* = 4 x tau/(c + delta) = c - delta. (ii) is the critical equation 1 - c s + x tau s^2 = 0 rewritten. (iii):
s (1 - tau) - (1 - x s)(1 - tau s) = -(x tau s^2 - c s + 1). (iv): differentiate the critical equation, x s*^2 + (2 x tau s* - c) s*_tau = 0,
and use (i). (v) is the critical equation. Envelope: d/dtau [phi_tau(s*(tau))] = (d_tau phi_tau)(s*) + phi_tau'(s*) s*_tau =
s*/(1 - tau s*), since d_tau phi_tau(s) = s/(1 - tau s); then (iii), and d/dt = (1 - tau) d/dtau.
(c) By (ii), x/(1 - s)^2 = 1/(x s^2 (1 - tau s)^2), so (1 - tau s)^2 phi_tau''(s) = [-(1 - tau s)^2 - 1/x + tau^2 s^2]/s^2, and
-(1 - tau s)^2 + tau^2 s^2 = 2 tau s - 1, 1 + 1/x = c/x. With (i), 2 tau s* - c/x = -delta/x; with 2 x tau s~ = c + delta,
2 tau s~ - c/x = delta/x. At tau = 0: phi_0''(1/c) = -c^2 - x/(1 - 1/c)^2 = -c^2 - c^2/x = -c^3/x. delta^2 = c^2 - 4 x tau has no zero on D,
so the Hessians do not vanish on D''. For (3.4): alg is (1 - tau s*), a unit with value 1, times the root with value 1 of an element of
1 + tau O, so alg is the square root of alg^2 = delta/(c^3 s*^2) with value 1 at 0; (delta/c)^(1/2)/(c s*) (principal root of
delta/c in 1 + tau O) is also such a root, and two square roots of a unit of O with the same value at 0 coincide. Finally 1/(c s*) =
(c + delta)/(2c).
(d) dU/ds = -tau and theta = U d/dU give d/ds = -(tau/U) theta. Differentiating once more, d/ds(-tau/(1 - tau s)) = -tau^2/(1 - tau s)^2,
which gives the second formula. At a critical point theta W = 0, so d^2/ds^2 [W(1 - tau s)] = (tau/U)^2 W_2, which is (3.5) by (c). For
W_2(U_2) use tau s~ = 1/(x s*). QED
[Mathematica: (V1)-(V2), checks (1), (3), (3b) of hvc_leaf_proof_check.wls; Perl: S1_checks.pl C1 (W_2(U_1) = -16 at
(x, tau) = (2, 5/8)), C3; S2_checks.pl part A.]

### 3.3 The rescaled derivatives tau^k W_k(U_1)
For k >= 1,
    phi_tau^(k)(s) = (k - 1)! [(-1)^(k-1) s^-k - x (1 - s)^-k + tau^k (1 - tau s)^-k],   phi_0^(k)(1/c) = (k - 1)! c^k ((-1)^(k-1) - x^(1-k)),   (3.6)
since (1 - 1/c)^-k = c^k x^-k. For a in {1, Q} put p_a := (1 - a U)^-1. Then theta log(1 - a U) = 1 - p_a and theta p_a = p_a^2 - p_a. Define
polynomials Pol_0(p) := p and Pol_n(p) := (p^2 - p) Pol_{n-1}'(p) (Pol_1 = p^2 - p, Pol_2 = 2p^3 - 3p^2 + p). By induction
theta^k log(1 - a U) = -Pol_{k-1}(p_a) for k >= 2, while theta^k log U = 0 for k >= 2. Hence
    W_k = -x Pol_{k-1}(p_Q) - Pol_{k-1}(p_1)      (k >= 2).                                                                   (3.7)
Pol_n has degree n + 1, leading coefficient n!, and Pol_n(0) = 0; write Pol_{k-1}(p) = Sum_{i=1}^{k} pc_{k-1,i} p^i (pc_{k-1,k} = (k - 1)!).

**Lemma 3.2 (rescaled derivatives).** For k >= 2 and tau in D'' put wr_k(tau) := tau^k W_k(U_1). Then
    wr_k = -Sum_{i=1}^{k} pc_{k-1,i} tau^(k-i) [ x (-(1 - tau)/(1 - s*))^i + s*^(-i) ].                                          (3.8)
The right side is analytic on D', tau = 0 included, and
    wr_k(0) = -(k - 1)! c^k (1 + (-1)^k x^(1-k)) = (-1)^k phi_0^(k)(1/c),     wr_2 = -delta/(x s*^2) (no zero on D').
Proof. For tau in D'', U_1 = 1 - tau s* is not in {1, 1/Q} (tau s* != 0, and Q U_1 = 1 would force s* = 1), so W_k(U_1) is defined (by (3.7)
W_k has poles only at U = 1 and U = 1/Q). At U_1, 1 - U_1 = tau s* and 1 - Q U_1 = -Q tau (1 - s*); with 1/Q = 1 - tau,
tau p_1(U_1) = 1/s* and tau p_Q(U_1) = -(1 - tau)/(1 - s*). Multiplying (3.7) by tau^k and writing tau^k Pol_{k-1}(p) =
Sum_i pc_{k-1,i} tau^(k-i) (tau p)^i gives (3.8). Only tau^(k-i) with 0 <= k - i <= k - 1 occur, and s*, 1/s*, 1/(1 - s*) are analytic on
D' (Lemma 3.1(a)), so the right side is analytic on D'. At tau = 0 only i = k survives; with s*(0) = 1/c, 1 - s*(0) = x/c and pc_{k-1,k} = (k - 1)!,
wr_k(0) = -(k - 1)! [x (-c/x)^k + c^k] = -(k - 1)! c^k (1 + (-1)^k x^(1-k)), which equals (-1)^k phi_0^(k)(1/c) by (3.6). Finally
wr_2 = tau^2 W_2(U_1) = -delta/(x s*^2) by (3.5). QED
So W_k(U_1) = tau^-k wr_k = O(tau^-k) and W_2(U_1)^-1 = O(tau^2): the coordinate y = log U behaves at t = 0 like -tau times a good
coordinate (Remark 3.11). [S1_checks.pl C6: (3.8) against tau^k W_k(U_1) from (3.7) at (x, tau) = (2, 5/8) for k <= 12; wr_k(0) for
k <= 12 and wr_2 = -delta/(x s*^2) as tau-series at x = 2, 1/3, -5/2.]

### 3.4 Formal Laplace expansions [formal, except Lemma 3.8]
Let Acal be the set of formal sums K = Sum_{p >= 0, q in Z} K_pq v^p z^q (coefficients complex numbers, or analytic functions of
parameters) such that, with the weight wt(p, q) := p/2 + q, (A1) wt is bounded below on the support of K, and (A2) for every N only
finitely many (p, q) with wt <= N occur. Acal is a ring (weights add; every coefficient of a product is a finite sum). It contains
C[[v]][z, 1/z] and Exp[Rr(v)/z] := Sum_j Rr^j/(j! z^j) for every Rr in v^3 C[[v]] ((Rr/z)^j has weight >= j/2). d/dv is a derivation of
Acal lowering weights by 1/2, multiplication by v raises them by 1/2, and d/dv Exp[Rr/z] = (Rr'/z) Exp[Rr/z]. For a != 0 the Gaussian
functional (formal moments of a centred Gaussian of variance -z/a) is
    <<K>>_a := Sum_{p even, q} K_pq (p - 1)!! (-z/a)^(p/2) z^q   in Kz       ((-1)!! := 1).
Each power z^w receives contributions only from monomials of weight w, finitely many by (A2), and w is bounded below by (A1); <<.>>_a is
linear over Kz and over functions of parameters, and maps elements of weight >= N into z^N C[[z]].

**Lemma 3.3 (Gaussian calculus).** For K in Acal:
 (G0) scaling: <<K(lambda v)>>_{lambda^2 a} = <<K(v)>>_a for lambda != 0;
 (G1) integration by parts: <<dK/dv>>_a = -(a/z) <<v K>>_a;
 (G2) variance derivative: if a = a(tau) and K depends analytically on tau (coefficientwise, with tau-independent support),
      d/dtau <<K>>_{a(tau)} = <<dK/dtau>>_a + (a_tau/(2z)) <<v^2 K>>_a + (a_tau/(2a)) <<K>>_a.
Proof. Each side is <<Lop K>> for an operator Lop on Acal that shifts weights by a bounded amount, so every z-coefficient involves finitely
many monomials of K and it suffices to take K = v^p (for (G2) with the product rule for the coefficients). Put m_p := <<v^p>>_a, so
m_{2n} = (2n - 1)!! (-z/a)^n and m_odd = 0. (G0): lambda^(2n) (2n - 1)!! (-z/(lambda^2 a))^n = m_{2n}. (G1): for p = 2n + 1,
<<p v^(2n)>> = (2n + 1)!! (-z/a)^n = -(a/z) (2n + 1)!! (-z/a)^(n+1) = -(a/z) m_{2n+2}; for p even both sides vanish. (G2):
d/dtau (-z/a)^n = -n (a_tau/a)(-z/a)^n, and (a_tau/(2z)) m_{2n+2} + (a_tau/(2a)) m_{2n} = m_{2n} (a_tau/(2a)) [-(2n + 1) + 1]. QED

**Definition 3.4.** Data: a point u0, a phase f analytic near u0 (only its Taylor series matters) with f'(u0) = 0 and f2 := f''(u0) != 0,
and an amplitude am in O_{u0}[z, 1/z] (a Laurent polynomial in z with coefficients analytic near u0). Write
f(u0 + v) = f(u0) + f2 v^2/2 + Rr_f(v), Rr_f(v) := Sum_{k>=3} f^(k)(u0) v^k/k!, and define
    Lnorm[f, am; u0] := << am(u0 + v) Exp[Rr_f(v)/z] >>_{f2}   in Kz                (Gaussian-normalised formal Laplace expansion),
    Lhat[f, am; u0]  := Lnorm[f, am; u0]/am(u0)   for am in O_{u0}, am(u0) != 0       (normalised: leading coefficient 1),
    Lfull[f, am; u0] := e^{f(u0)/z} (2 pi z)^(1/2) (-f2)^(-1/2) Lnorm[f, am; u0]   (unnormalised; a formal symbol, with a chosen root).
Lfull is the formal counterpart of Int e^{f/z} am du over a path through u0 (Lemma 3.8). Immediate properties:
 (D1) Wick form: Lnorm[f, am; u0] = Sum_{n>=0} ((2n)!/(n! 2^n)) (-z/f2)^n [v^(2n)] (am(u0 + v) Exp[Rr_f(v)/z]). For am = 1 the coefficient of z^m
      is wickZ_m with W_k -> f^(k)(u0) (compare (C7)).
 (D2) For am in O_{u0}, Lnorm[f, am; u0] lies in C[[z]] with constant term am(u0), and its z^m coefficient is a polynomial with rational
      coefficients in 1/f2, f^(3)(u0), ..., f^(2m+2)(u0), am(u0), ..., am^(2m)(u0). [Both factors have weight >= 0; the weight-m part involves
      monomials v^a of am with a/2 <= m and terms f^(k) v^k/(k! z) of weight k/2 - 1 <= m.]
 (D3) Lnorm depends only on the Taylor coefficients of f and am at u0 (so a translation u = u0 + v changes nothing), is unchanged when a
      constant is added to f, and is linear in am over C[z, 1/z] and over functions independent of u.

**Proposition 3.5 (P3: total derivatives).** For every km in O_{u0}[z, 1/z]: Lnorm[f, km' + f' km/z; u0] = 0. (The amplitude
km' + f' km/z is that of a total derivative: e^{f/z}(km' + f' km/z) = d/du (e^{f/z} km).)
Proof. Put K(v) := km(u0 + v) Exp[Rr_f(v)/z] in Acal. Since f'(u0 + v) = f2 v + Rr_f'(v), (km' + f' km/z)(u0 + v) Exp[Rr_f(v)/z] =
dK/dv + (f2/z) v K, and <<dK/dv>>_{f2} + (f2/z) <<v K>>_{f2} = 0 by (G1). QED

**Proposition 3.6 (P2: derivative in a parameter with a moving expansion point).** Let f(tau, s) and am(tau, s) be analytic near
(tau0, s0) (am may be a Laurent polynomial in z with such coefficients), let sc(tau) be analytic with sc(tau0) = s0,
d_s f(tau, sc(tau)) = 0 and f2(tau) := d_s^2 f(tau, sc(tau)) != 0. Put f0(tau) := f(tau, sc(tau)) and Nn(tau) := Lnorm[f(tau, .),
am(tau, .); sc(tau)]. Then every z-coefficient of Nn is analytic in tau, and
    d_tau Nn + (f0_tau/z - f2_tau/(2 f2)) Nn = Lnorm[f(tau, .), Dam am; sc(tau)],     Dam am := d_tau am + am (d_tau f)/z,        (3.9)
where f0_tau = (d_tau f)(tau, sc(tau)) (by the critical-point condition) and f2_tau is the total derivative of f2. Equivalently, with
(-f2(tau))^(-1/2) chosen continuously, d_tau Lfull[f, am; sc(tau)] = Lfull[f, Dam am; sc(tau)]: the expansion of the integrand e^{f/z} am
commutes with d/dtau although the expansion point moves.
Proof. Analyticity follows from (D2), as d_s^k f(tau, sc(tau)), d_s^k am(tau, sc(tau)) and 1/f2(tau) are analytic.
Step 1 (shifted variable). Let ft(tau, v) := f(tau, sc(tau) + v), amt(tau, v) := am(tau, sc(tau) + v). By (D3) Nn = Lnorm[ft, amt; 0], and
d_v ft(tau, 0) = 0 for all tau. Write ft = f0(tau) + f2(tau) v^2/2 + Rt(tau, v) and Kt := amt Exp[Rt/z] in Acal. By (G2), with
d_tau Kt = (d_tau amt + amt d_tau Rt/z) Exp[Rt/z],
    d_tau Nn = <<(d_tau amt + amt d_tau Rt/z) Exp[Rt/z]>> + (f2_tau/(2z)) <<v^2 Kt>> + (f2_tau/(2 f2)) Nn.
Since d_tau ft = f0_tau + f2_tau v^2/2 + d_tau Rt,
    Lnorm[ft, d_tau amt + amt d_tau ft/z; 0] = <<(d_tau amt + amt d_tau Rt/z) Exp[Rt/z]>> + (f0_tau/z) Nn + (f2_tau/(2z)) <<v^2 Kt>>
                                            = d_tau Nn + (f0_tau/z - f2_tau/(2 f2)) Nn.
Step 2 (back to s). d_tau amt = (d_tau am + sc_tau d_s am)(tau, sc + v) and d_tau ft = (d_tau f + sc_tau d_s f)(tau, sc + v), so
d_tau amt + amt d_tau ft/z = [Dam am + sc_tau (d_s am + (d_s f) am/z)](tau, sc + v). The last term is sc_tau times an amplitude of the form
km' + f' km/z (km = am(tau, .)), whose Lnorm vanishes by (P3). Translating back by (D3) gives (3.9). The unnormalised form follows from
d_tau [e^{f0/z} (-f2)^(-1/2) Nn] = e^{f0/z} (-f2)^(-1/2) [d_tau Nn + (f0_tau/z - f2_tau/(2 f2)) Nn]. QED
In (3.9), f0_tau/z comes from the moving critical value e^{f0/z}, -f2_tau/(2 f2) from the Gaussian normalisation (-f2)^(-1/2), and the
motion of the expansion point is a total derivative.

**Proposition 3.7 (P1: coordinate invariance).** Let hc be analytic near v0 with hc(v0) = u0 and hc1 := hc'(v0) != 0, and f, am as in
Definition 3.4 with am in O_{u0}. Then f o hc has a nondegenerate critical point at v0, (f o hc)(v0) = f(u0), (f o hc)''(v0) = f2 hc1^2, and
    Lnorm[f o hc, (am o hc) hc'; v0] = hc1 Lnorm[f, am; u0].                                                                   (3.10)
Exactly what is invariant: (i) the critical value e^{f(u0)/z}; (ii) the unnormalised expansion, Lfull[f o hc, (am o hc) hc'; v0] =
Lfull[f, am; u0], when the roots are related by (-(f o hc)''(v0))^(1/2) = hc1 (-f2)^(1/2) (the Jacobian hc' in the amplitude is
compensated by hc1 in the Gaussian normalisation; for real data and hc1 > 0 both roots are positive, hc1 < 0 is the orientation reversal
of the path); (iii) the normalised expansion, with no convention at all: Lhat[f o hc, (am o hc) hc'; v0] = Lhat[f, am; u0], because
((am o hc) hc')(v0) = hc1 am(u0).
Proof. For eps in C let hc_eps(v) := u0 + hc1 (v - v0) + eps (hc(v) - u0 - hc1 (v - v0)); then hc_eps(v0) = u0 and hc_eps'(v0) = hc1 for every eps,
hc_0 is affine and hc_1 = hc. On a neighbourhood V of v0 (uniform for |eps| <= 2) hc_eps' has no zero and hc_eps(V) lies in the domain of f, am.
Put f_eps := f o hc_eps, am_eps := (am o hc_eps) hc_eps' and n(eps) := Lnorm[f_eps, am_eps; v0]. Each f_eps has a critical point at v0 with
critical value f(u0) and f_eps''(v0) = f2 hc1^2, independent of eps; the Taylor coefficients of f_eps, am_eps at v0 are polynomials in eps,
so by (D2) each z-coefficient of n(eps) is a polynomial in eps. (3.9) with tau = eps, sc(eps) = v0, f0_eps = 0, f2_eps = 0 gives
n'(eps) = Lnorm[f_eps, d_eps am_eps + am_eps d_eps f_eps/z; v0]. Let eta := d_eps hc_eps = hc(v) - u0 - hc1 (v - v0) and xi := eta/hc_eps'.
Then d_eps f_eps = (f' o hc_eps) eta = xi f_eps', and d_eps am_eps = (am' o hc_eps) eta hc_eps' + (am o hc_eps) eta' = (xi am_eps)' (because
xi am_eps = eta (am o hc_eps)). So the amplitude is (xi am_eps)' + f_eps' (xi am_eps)/z, a total derivative, and n'(eps) = 0 by (P3); hence
n(1) = n(0). Finally, with w = v - v0, f_0(v0 + w) = f(u0 + hc1 w) has Rr_{f_0}(w) = Rr_f(hc1 w) and f_0''(v0) = hc1^2 f2, so by (G0) with
lambda = hc1, n(0) = hc1 <<am(u0 + hc1 w) Exp[Rr_f(hc1 w)/z]>>_{hc1^2 f2} = hc1 Lnorm[f, am; u0]. (ii) and (iii) restate (3.10). QED

**Lemma 3.8 (Laplace's method) [analytic].** Let I be an open real interval, f and am real-analytic on I, u0 in I with f'(u0) = 0,
f''(u0) < 0, and assume sup_{I - (u0 - e, u0 + e)} f < f(u0) for every e > 0 and Int_I e^{f/z0} |am| du < oo for some z0 > 0. Then for
every M, as z -> 0+,
    Int_I e^{f(u)/z} am(u) du = e^{f(u0)/z} (2 pi z/(-f''(u0)))^(1/2) [ Sum_{m<M} ([z^m] Lnorm[f, am; u0]) z^m + O(z^M) ].
Proof. Localisation: fix e > 0 small and eta > 0 with f <= f(u0) - eta outside (u0 - e, u0 + e). There, for 0 < z < z0,
e^{f/z} = e^{f/z0} e^{f (1/z - 1/z0)} <= e^{f/z0} e^{(f(u0) - eta)(1/z - 1/z0)}, so the outer part is O(e^{(f(u0) - eta)/z}) = e^{f(u0)/z} O(z^M) for
every M. Morse coordinate: f(u) - f(u0) = (u - u0)^2 qq(u) with qq real-analytic and qq(u0) = f''(u0)/2 < 0; shrinking e, qq < 0 on
(u0 - e, u0 + e), and wm(u) := (u - u0)(2 qq(u)/f''(u0))^(1/2) is real-analytic and strictly increasing, wm(u0) = 0, wm'(u0) = 1, with
f(u) = f(u0) + f''(u0) wm^2/2. Let u(wm) be the inverse map, [w_-, w_+] (w_- < 0 < w_+) the image interval and Gw(wm) := am(u(wm)) u'(wm).
Watson: with al := -f''(u0) > 0 and Gw(wm) = Sum_{p<2M} Gw_p wm^p + O(|wm|^(2M)) on [w_-, w_+], each Int_{w_-}^{w_+} e^{-al wm^2/(2z)} wm^p dwm
equals the integral over R ((2 pi z/al)^(1/2) (p - 1)!! (z/al)^(p/2) for p even, 0 for p odd) up to O(e^{-b/z}) with b > 0, and
Int_R e^{-al wm^2/(2z)} |wm|^(2M) dwm = O(z^(M + 1/2)). So the inner integral is e^{f(u0)/z} (2 pi z/al)^(1/2) [Sum_{n<M} Gw_{2n} (2n - 1)!!
(z/al)^n + O(z^M)], and Sum_n Gw_{2n} (2n - 1)!! (z/al)^n = <<Gw>>_{-al} = Lnorm[f''(u0) wm^2/2, Gw; 0]. Finally (3.10) with hc = u(wm),
hc1 = 1, and (D3) for the constant f(u0) give Lnorm[f o u, (am o u) u'; 0] = Lnorm[f, am; u0]. QED
(The lemma is proved in full above. Cross-reference only: the same coefficient formula is the L_j of Hormander's stationary-phase theorem,
ALPDO I, Thm 7.7.5, which applies to the Laplace case with omega = 1/z and F := -i (f - f(u0)), Im F = f(u0) - f >= 0 near u0. Lemma 3.8 is
used only for Beta integrals, Propositions 3.15 and 8.6(a).)

**Lemma 3.9 (homogeneity of the Wick coefficients).** Every monomial of wickZ_m is a rational multiple of W_2^-n prod_{i=1}^{n-m} W_{k_i}
with all k_i >= 3 and Sum_i k_i = 2n. Consequently, for every lambda != 0, as an identity of rational functions of the W_k,
    wickZ_m(lambda^2 W_2, lambda^3 W_3, ..., lambda^(2m+2) W_{2m+2}) = wickZ_m(W_2, W_3, ..., W_{2m+2}),                         (3.11)
and Lnorm[f(u0 + lambda .), 1; 0] = Lnorm[f, 1; u0] for f as in Definition 3.4.
Proof. By (C7) the monomials of the z^m coefficient are W_2^-n prod_{k in p} W_k with Sum_{k in p} k = 2n; under W_k -> lambda^k W_k such a
monomial is multiplied by lambda^(-2n) lambda^(2n) = 1. For the second statement put f_lambda(v) := f(u0 + lambda v); then f_lambda'(0) = 0 and
f_lambda^(k)(0) = lambda^k f^(k)(u0), so by (D1) the two z^m coefficients are wickZ_m at the arguments lambda^k f^(k)(u0) and f^(k)(u0), which
agree by (3.11). (Equivalently: (G0) for K = Exp[Rr_f/z].) QED
In words: giving W_k the weight k, every Wick monomial has weight -2n + Sum k_i = 0, so the Wick formula is blind to linear rescalings of
the coordinate. This is why the degenerate coordinate y = log U is harmless.

### 3.5 The code's row is the normalised expansion of the Euler integrand (Step 2(a)) [formal]
**Proposition 3.10.** (i) For every tau in D'' (every small t != 0), coefficientwise in z,
    r^pt_{0,1}(z; t) = N(z, tau) := Lhat[phi_tau, (1 - tau s)^-1; s*(tau)] = U_1 Lnorm[phi_tau, (1 - tau s)^-1; s*(tau)].              (3.12)
(ii) For every m >= 0 and tau in D'',
    [z^m] r^pt_{0,1}(z; t) = wickZ_m(W_2(U_1), W_3(U_1), ...) = wickZ_m(wr_2(tau), wr_3(tau), ..., wr_{2m+2}(tau)),                   (3.13)
a polynomial in 1/wr_2, wr_3, ..., wr_{2m+2} with rational coefficients. Hence the z^m coefficient of r^pt_{0,1}, an algebraic function
of Q, extends analytically to all of D', tau = 0 included, and its value at tau = 0 is Abeta_m := [z^m] Abeta(z), where
    Abeta(z) := Lnorm[phi_0, 1; 1/c].
(iii) The z^m coefficient N_m of N is analytic on D' and N(z, 0) = Abeta(z). So (3.12) holds on all of D', N lies in O(D')[[z]], and the
Taylor series at t = 0 of the z^m coefficient of r^pt_{0,1} (what Series[... /. Q -> Exp[t], {t, 0, M}] in evalQseries and in the check
scripts computes) is that of N_m(tau(t)): (3.12) holds in C[[z]][[t]].
Proof. (1) The code's row is a normalised expansion in y = log U. Fix tau in D''. Then U_1 != 0, so y1 := log U_1 (any branch) is defined
and Wy(y) := W(e^y) is analytic near y1, with Wy^(k)(y1) = W_k(U_1), Wy'(y1) = 0 and Wy''(y1) = W_2(U_1) = -delta/(x tau^2 s*^2) != 0 by
(3.5). In this coordinate the integrand is e^{W/z} dU/U = e^{Wy(y)/z} dy, amplitude 1. By (D1) and (C7), r^pt_{0,1}(z; t) =
Lnorm[Wy, 1; y1] = Lhat[Wy, 1; y1] (to every order; the code's truncation at z^ZMAX is the truncation of this series).
(2) Change of coordinate; proof of (i). hc(s) := log(1 - tau s), with the branch for which hc(s*) = y1, is analytic near s*, and
hc'(s*) = -tau/U_1 != 0. By (3.1), Wy o hc = phi_tau + const near s*, and the amplitude 1 becomes (1 o hc) hc' = -tau (1 - tau s)^-1. By
(P1)(iii), (D3) and linearity (the factor -tau), Lhat[Wy, 1; y1] = Lhat[phi_tau, -tau (1 - tau s)^-1; s*] = Lhat[phi_tau, (1 - tau s)^-1; s*]
= N(z, tau); and Lhat = Lnorm/(1 - tau s*)^-1 = U_1 Lnorm.
(3) Proof of (ii). The first equality in (3.13) is the definition (C7). By (3.11) with lambda = tau, wickZ_m(W_k(U_1)) = wickZ_m(wr_k). By
Lemma 3.2 every wr_k is analytic on D' and wr_2 has no zero there, so the right side of (3.13) is analytic on D'. At tau = 0,
wr_k(0) = (-1)^k phi_0^(k)(1/c), and (3.11) with lambda = -1 followed by (D1) (phi_0'(1/c) = 0, phi_0''(1/c) = -c^3/x != 0) gives
wickZ_m(wr_k(0)) = wickZ_m(phi_0^(k)(1/c)) = [z^m] Lnorm[phi_0, 1; 1/c] = Abeta_m. No coordinate change, i.e. no use of (P1), occurs here.
(4) Proof of (iii). By (D2), N_m is a polynomial in 1/phi_tau''(s*), in phi_tau^(k)(s*) (k <= 2m + 2), in d_s^a (1 - tau s)^-1 at s* =
a! tau^a U_1^(-a-1) (a <= 2m), and in U_1. By (3.6) all of these are analytic on D', since s* is analytic there, s* is not in
{0, 1, 1/tau}, and phi_tau''(s*) != 0 (Lemma 3.1). At tau = 0, s*(0) = 1/c and (1 - 0 s)^-1 = 1, so N(z, 0) = Abeta(z). By (i) the analytic
functions of (ii) and N_m agree on D'', and by (ii) and N(z, 0) = Abeta(z) also at tau = 0. Finally t -> tau = 1 - e^-t is analytic with
tau(0) = 0. QED

**Remark 3.11 (the degenerate coordinate; a uniform version without continuation).** In y the inputs of the Wick formula blow up at
t = 0 (W_k(U_1) = tau^-k wr_k, W_2(U_1)^-1 = O(tau^2), hc'(s*) = O(tau)), yet every single Wick monomial is analytic there: a monomial
W_2^-n prod W_{k_i} has Sum k_i = 2n, so it equals wr_2^-n prod wr_{k_i}; no cancellation between monomials is needed. Geometrically y is
-tau times a coordinate that is good uniformly near tau = 0. Let
    ytilde := -log(1 - tau s)/tau = Sum_{n>=1} tau^(n-1) s^n/n,    s = sigma_tau(ytilde) := (1 - e^{-tau ytilde})/tau = Sum_{n>=1} (-tau)^(n-1) ytilde^n/n!
(ytilde = s at tau = 0). sigma is entire in (tau, ytilde), sigma_0 = id, d sigma_tau/d ytilde = e^{-tau ytilde} = 1 - tau sigma_tau(ytilde),
U = e^{-tau ytilde}, i.e. y = -tau ytilde and d/d ytilde = -tau theta. By (3.1),
    W(e^{-tau ytilde}) = phi_tau(sigma_tau(ytilde)) + c log tau + x log Q + x log(-1),                                       (3.14)
and the first term on the right is analytic in (tau, ytilde) near (0, 1/c). With ytilde1(tau) := -log(U_1)/tau (principal branch;
ytilde1(0) = 1/c, sigma_tau(ytilde1) = s*), differentiating (3.14) k >= 2 times gives
    (-tau)^k W_k(U_1) = (-1)^k wr_k(tau) = d^k/d ytilde^k [phi_tau o sigma_tau](ytilde1(tau)),                                   (3.15)
which re-proves Lemma 3.2 near tau = 0. For tau != 0, Lemma 3.9 (lambda = -tau), (3.14) and (D3) give r^pt_{0,1}(z; t) = Lnorm[Wy, 1; y1]
= Lnorm[phi_tau o sigma_tau, 1; ytilde1]. The right side is defined for every small tau, tau = 0 included (second derivative wr_2 != 0,
-c^3/x at tau = 0; coefficients analytic in tau by (D2)). (P1) with hc = sigma_tau, v0 = ytilde1, u0 = s*, hc1 = U_1 and am = (1 - tau s)^-1,
for which (am o hc) hc' = e^{-tau ytilde}/(1 - tau sigma_tau(ytilde)) = 1, gives Lnorm[phi_tau o sigma_tau, 1; ytilde1] =
U_1 Lnorm[phi_tau, (1 - tau s)^-1; s*] = N(z, tau) at every small tau, tau = 0 included. [S1_checks.pl C6: (3.15) for k <= 12, and
wickZ_m(wr_k) = Lnorm[phi_tau o sigma_tau, 1; ytilde1] = N mod (tau^7, z^5), with value Abeta = h_w c/(c + z) at tau = 0; C1: (3.12)
through z^5 at (x, tau) = (2, 5/8) and (1/3, 7/12), the y-side computed exactly as the code does, e.g. 1, -97/512, 20809/524288,
1164725/268435456, ... at (2, 5/8).]

### 3.6 Gauss's equation for the formal expansion (Step 2(b)) [formal]
The amplitude identity. Put f := phi_tau, am0 := (1 - tau s)^-1, and Dam q := d_tau q + q (d_tau f)/z as in (3.9), so that
e^{f/z} Dam q = d_tau (e^{f/z} q). Since d_tau f = s/(1 - tau s),
    Dam am0 = (1 + 1/z) s (1 - tau s)^-2,      Dam^2 am0 = (1 + 1/z)(2 + 1/z) s^2 (1 - tau s)^-3,
and, with kE := s (1 - s)(1 - tau s)^-2,
    tau (1 - tau) Dam^2 am0 + (gamma - (alpha + beta + 1) tau) Dam am0 - alpha beta am0 = -(1 + 1/z) (kE' + f' kE/z).                 (3.16)
This is an identity of Laurent polynomials in z with coefficients rational in (s, tau), valid identically in (s, tau) (not only at a
critical point): since f' kE = (x tau s^2 - c s + 1)(1 - tau s)^-3, after multiplication by (1 - tau s)^3/(1 + 1/z) both sides equal
(1/z)(-x tau s^2 + c s - 1) + (2 s - tau s - 1). [S1_checks.pl C4: 6 random rational points (x, s, tau, z).] (3.16) is Euler's
integration by parts (Remark 3.20) written for amplitudes.

**Proposition 3.12.** Let L(tau) := Lnorm[phi_tau, (1 - tau s)^-1; s*(tau)] (in O(D')[[z]], by Proposition 3.10(iii)) and
    NG(z, tau) := Exp[dphi(tau)/z] (phi_tau''(s*)/phi_0''(1/c))^(-1/2) iota(L)   in Kz[[tau]]
(root with value 1 at tau = 0). Then Gauss[NG] = 0 in Kz[[tau]] (d/dtau-derivatives), NG(z, 0) = Abeta(z), an element of
1 + z Q(x)[[z]], and NG = Exp[dphi/z] N/alg.
Proof. Here f0(tau) = phi_tau(s*(tau)), so f0_tau = dphi_tau, and f2(tau) = phi_tau''(s*(tau)). Let Dop := d/dtau + dphi_tau/z - f2_tau/(2 f2),
acting on O((z))-valued expressions. (P2), applied to each z-power of an amplitude q, says Dop Lnorm[q] = Lnorm[Dam q]; applying it twice
and using (D3) (linearity over C[tau][z, 1/z]),
    [tau (1 - tau) Dop^2 + (gamma - (alpha + beta + 1) tau) Dop - alpha beta] L = Lnorm[LHS of (3.16)] = -(1 + 1/z) Lnorm[kE' + f' kE/z] = 0
by (P3). Let E := Exp[dphi/z] (f2(tau)/f2(0))^(-1/2) in C[1/z][[tau]]. In Kz[[tau]], d/dtau (E Y) = E Dop Y, because
E_tau = E (dphi_tau/z - f2_tau/(2 f2)) (iota commutes with d/dtau and with multiplication by convergent functions of tau, and
(Exp[dphi/z])_tau = (dphi_tau/z) Exp[dphi/z]). Hence Gauss[E iota(L)] = E iota([tau (1 - tau) Dop^2 + ...] L) = 0. At tau = 0:
s*(0) = 1/c, (1 - 0 s)^-1 = 1 and dphi(0) = 0, so NG(z, 0) = Lnorm[phi_0, 1; 1/c] = Abeta(z), which lies in 1 + z Q(x)[[z]] by (D2) and
(3.6). Finally N = U_1 L and alg = U_1 (f2(tau)/f2(0))^(1/2), so NG = Exp[dphi/z] N/alg. QED
No boundary term occurs at any stage: the formal functional kills total derivatives exactly (P3), the formal shadow of the vanishing
boundary terms of Euler's integration by parts, valid for all x and z.

### 3.7 Uniqueness and the closed form (Step 2(c)) [formal]
**Lemma 3.13 (uniqueness).** If Y = Sum_k Y_k tau^k in Kz[[tau]] satisfies Gauss[Y] = 0, then Y = Y_0 F.
Proof. The coefficient of tau^k of Gauss[Y] is (k + 1) k Y_{k+1} - k (k - 1) Y_k + gamma (k + 1) Y_{k+1} - (alpha + beta + 1) k Y_k
- alpha beta Y_k = (k + 1)(k + gamma) Y_{k+1} - (k + alpha)(k + beta) Y_k. Now k + gamma = (c + (k + 2) z)/z is a unit of Kz (c != 0) and
k + 1 != 0, so Y_{k+1} = (k + alpha)(k + beta) Y_k/((k + 1)(k + gamma)), whence Y_k = Y_0 F_k. QED
So the solutions of Gauss's equation in Kz[[tau]] form the line Kz F; the second local solution, of exponent 1 - gamma = -c/z - 1, is not
in this ring. (The same recursion shows that F itself satisfies Gauss[F] = 0: F_{n+1}/F_n = (n + alpha)(n + beta)/((n + 1)(n + gamma)).)

**Theorem 3.14.** In Kz[[tau]]:
    N(z, tau) = Abeta(z) F(z, tau) Ephi(z, tau) alg(tau).                                                                      (3.17)
Proof. By Proposition 3.12 and Lemma 3.13, NG = NG(z, 0) F = Abeta F. Multiply by alg Ephi (alg is in O, Ephi is the inverse of
Exp[dphi/z] in C[1/z][[tau]]) and use NG = Exp[dphi/z] N/alg. QED
[S1_checks.pl C3: (3.17) mod (tau^5, z^5) at x = 2, 1/3, with N from the Wick formula in the coordinate s; C6: the left side computed as
the code does, in y = log U, via (3.13), mod (tau^7, z^5) at x = 2, 1/3, -5/2.]

### 3.8 The constant Abeta(z)
**Proposition 3.15.** Abeta(z) := Lnorm[phi_0, 1; 1/c] = h_w(z) c/(c + z) in Q(x)[[z]]. In words: the normalised Laplace expansion of the
Beta integral B(1/z + 1, x/z + 1) = Int_0^1 e^{phi_0(s)/z} ds is the Stirling tail h_w of B(x/z, 1/z) times the elementary factor c/(c + z)
coming from B(a + 1, b + 1) = B(a, b) a b/((a + b)(a + b + 1)).
Proof. [analytic for x > 0, then algebraic in x.]
(1) Reduction to x > 0. Both sides have z-coefficients in Q(x): the left by (D2) and (3.6), with 1/phi_0''(1/c) = -x/c^3; the right by
definition (c/(c + z) = Sum_n (-z/c)^n). Two rational functions of x that agree on (0, oo) are equal. Fix a real x > 0.
(2) Laplace. For real z > 0, B(1/z + 1, x/z + 1) = Int_0^1 s^(1/z) (1 - s)^(x/z) ds = Int_0^1 e^{phi_0(s)/z} ds. On I = (0, 1), phi_0 =
log s + x log(1 - s) is real-analytic, phi_0'' = -1/s^2 - x/(1 - s)^2 < 0, phi_0' = (1 - c s)/(s (1 - s)) vanishes only at 1/c, and
phi_0 -> -oo at both ends; so phi_0 is strictly concave with its unique maximum at 1/c, its supremum off (1/c - e, 1/c + e) is
max(phi_0(1/c - e), phi_0(1/c + e)) < phi_0(1/c), and Int_0^1 e^{phi_0} ds = B(2, x + 1) < oo. Lemma 3.8 (am = 1, z0 = 1), with
-phi_0''(1/c) = c^3/x and phi_0(1/c) = x log x - c log c, gives
    B(1/z + 1, x/z + 1) = e^{(x log x - c log c)/z} (2 pi z x/c^3)^(1/2) [Sum_{m<M} Abeta_m z^m + O(z^M)]     (z -> 0+).
(3) Stirling. Gamma(u + 1) = u Gamma(u) gives B(a + 1, b + 1) = B(a, b) a b/((a + b)(a + b + 1)); for a = 1/z, b = x/z this factor is
(x/z^2)/((c/z)((c + z)/z)) = x/(c (c + z)). Stirling's series with remainder (DLMF 5.11.1; for real positive argument the remainder after
the term k = K is bounded in absolute value by the first neglected term B_{2K+2}/((2K + 2)(2K + 1) u^(2K+1)) and has its sign, DLMF
s.5.11(ii); Whittaker-Watson s.12.33): for u -> +oo and every K,
    log Gamma(u) = (u - 1/2) log u - u + (1/2) log(2 pi) + Sum_{k=1}^{K} bk_k u^(1-2k) + O(u^(-2K-1)).
Apply it to log B(1/z, x/z) = log Gamma(1/z) + log Gamma(x/z) - log Gamma(c/z) (all three arguments -> +oo). The terms -u cancel; the
coefficient of log z in the (u - 1/2) log u terms is -(1/z - 1/2) - (x/z - 1/2) + (c/z - 1/2) = 1/2, and the remaining ones give
(x log x - c log c)/z - (1/2) log x + (1/2) log c. With u^(1-2k) = z^(2k-1) w^(1-2k) for u = w/z,
    B(1/z, x/z) = e^{(x log x - c log c)/z} (2 pi c z/x)^(1/2) Exp[Sum_{k=1}^{K} bk_k St_k z^(2k-1) + O(z^(2K+1))].
The exponent is the truncation of log h_w (C8) up to O(z^(2K+1)), and the exponential of (polynomial + O(z^(2K+1))) is the truncation
of the formal exponential plus O(z^(2K+1)). Multiplying by x/(c (c + z)) (analytic at z = 0) and using (2 pi c z/x)^(1/2) x/c =
(2 pi z x/c^3)^(1/2) c,
    B(1/z + 1, x/z + 1) = e^{(x log x - c log c)/z} (2 pi z x/c^3)^(1/2) [ (h_w(z) c/(c + z))_{<M} + O(z^M) ].
(4) The prefactors in (2) and (3) are the same positive function of z, so the two power series are asymptotic expansions of the same
function as z -> 0+; asymptotic power series are unique, so they coincide. QED
Check at z^1, by hand: (D1) gives Abeta_1 = f4/(8 f2^2) - 5 f3^2/(24 f2^3) with f2 = -c^3/x, f3 = 2 c^3 (x^2 - 1)/x^2, f4 = -6 c^4 (x^3 + 1)/x^3,
i.e. Abeta_1 = (x^2 - 11x + 1)/(12 c x); and [z^1](h_w c/(c + z)) = (1/12)(1/x + 1 - 1/c) - 1/c = (x^2 - 11x + 1)/(12 c x). [S1_checks.pl
C2: through z^7 at x = 2, 1/3, 5 and at x = -3/7, outside x > 0, as the rationality argument predicts.] A purely formal proof would have
to reproduce Stirling's series inside the formal calculus; none is given. That the code's row at t = 0 equals Abeta (Proposition 3.10(ii))
is formal; only the evaluation of Abeta uses analysis.

### 3.9 Regularity, the frame-B row, F2 (Step 2(d)) [formal, given Proposition 3.15]
**Corollary 3.16 (regularity).** For every k >= 0, Psi_k = [tau^k](F Ephi alg) and H_k = [tau^k](F Ephi) are rational functions of z with
coefficients in Q(x), without pole at z = 0.
Proof. Psi_k = Sum_{i+j+l=k} F_i X_j [tau^l] alg is a finite sum of products of elements of Q(x)(z), Q(x)[1/z] and Q(x) (alg and dphi have
Taylor coefficients in Q(x)), so it lies in Q(x)(z). By (3.17), Psi_k = [tau^k](N/Abeta) in Kz: Abeta is a unit of C[[z]] (constant term 1)
and N lies in C[[z]][[tau]] (Proposition 3.10), so the Laurent expansion of Psi_k at z = 0 has no negative powers. Only Abeta in
1 + z C[[z]] is used here, not its value. alg is in 1 + tau Q(x)[[tau]] and does not depend on z, so H_k = Sum_{j<=k} Psi_j [tau^(k-j)](1/alg)
is a Q(x)-linear combination of Psi_0, ..., Psi_k. QED

**Theorem 3.17 (the frame-B singular row).** In C[[z]][[t]], and coefficientwise in z as analytic functions of t near 0,
    r_{0,1}(z; t) = (c/(c + z)) F(z, tau) Ephi(z, tau) alg(tau),     tau = 1 - e^-t.                                             (3.18)
Precisely: for every k, [t^k] r_{0,1}(z; t) = Sum_{j<=k} ([t^k] tau^j) (c/(c + z)) Psi_j(z), a rational function of z regular at z = 0, whose
Taylor coefficients in z are the t^k Taylor coefficients of the z-coefficients of r_{0,1}.
Proof. r_{0,1} = h_w^-1 r^pt_{0,1} (C8; h_w^-1 is t-independent) = h_w^-1 N (Proposition 3.10) = h_w^-1 Abeta F Ephi alg (3.17)
= (c/(c + z)) F Ephi alg (Proposition 3.15). The passage from tau to t is the triangular substitution (T) and preserves regularity.
Coefficientwise meaning: the right side is a priori an element of Kz[[tau]]; by Corollary 3.16 it lies in C[[z]][[tau]], and the chain
above says that it equals iota(h_w^-1 N) with h_w^-1 N in O(D')[[z]] (Proposition 3.10(iii)). So each z-coefficient of the right side is
the Taylor series of a function analytic on D', and its sum is the corresponding z-coefficient of r_{0,1}. By (3.20) below the right side
also equals (c/(c + z)) Phi alg with Phi in O(D)[[z]]. QED
[Mathematica (V3), check (2b) of hvc_leaf_proof_check.wls: rho times (3.18), against the stored genus-4 factors, all t^k with k <= 6 and
z^m with m <= 9. Independently of the factor files, S1_checks.pl C3 checks [tau^0] h_w^-1 N = c/(c + z) and
[tau^1] h_w^-1 N = (1 - x)^2 z/(c (c + z)(c + 2z)) through z^4.]

**Corollary 3.18 (F2).** r_{0,1}(z; 0) = c/(c + z), and r^pt_{0,1}(z; 0) = Abeta(z) = h_w(z) c/(c + z).
Proof. t = 0 means tau = 0, and F(z, 0) = 1, dphi(0) = 0, alg(0) = 1. (r^pt_{0,1}(z; 0) = Abeta is also Proposition 3.10(ii).) QED

### 3.10 A purely formal proof of the regularity: eikonal and transport equations
This uses neither the Laplace identification (s.3.4-3.5, s.3.8) nor the value of Abeta; it shows directly that the z-expansion of Ephi F
exists and has coefficients analytic in tau. Write F = Exp[Sph/z] Phi with Sph in tau O. From F_tau = e^{Sph/z}(Phi_tau + Sph_tau Phi/z),
F_tautau = e^{Sph/z}(Phi_tautau + 2 Sph_tau Phi_tau/z + Sph_tautau Phi/z + Sph_tau^2 Phi/z^2) and gamma = c/z + 2, alpha + beta + 1 = 2/z + 3,
alpha beta = 1/z^2 + 2/z + 1, as operators on Kz[[tau]],
    Exp[-Sph/z] Gauss[Exp[Sph/z] Phi] = z^-2 Eik(Sph_tau) Phi + z^-1 (P_tr Phi_tau + R_tr Phi) + T_tr Phi,                         (3.19)
    Eik(p) := tau (1 - tau) p^2 + (c - 2 tau) p - 1,      P_tr := 2 tau (1 - tau) Sph_tau + c - 2 tau,
    R_tr := tau (1 - tau) Sph_tautau + (2 - 3 tau) Sph_tau - 2,      T_tr := tau (1 - tau) d^2/dtau^2 + (2 - 3 tau) d/dtau - 1.
Eikonal equation. Eik(p) = 0 is quadratic in p; its root regular at tau = 0 has p(0) = 1/c, the other root is -c/tau + O(1) (exponent
1 - gamma). Sph := dphi solves it, with Sph_tau = (1 - x s*)/(1 - tau) (Lemma 3.1(b)): with s = s*,
    (1 - tau) Eik = tau (1 - x s)^2 + (c - 2 tau)(1 - x s) - (1 - tau) = x^2 tau s^2 - c x s + (c - 1) = x (x tau s^2 - c s + 1) = 0,
and Sph_tau(0) = 1/c. The coefficients: P_tr = (d Eik/dp)(Sph_tau), so P_tr^2 = (c - 2 tau)^2 + 4 tau (1 - tau) = c^2 - 4 x tau and
P_tr(0) = c, i.e. P_tr = delta (directly: P_tr - c + 2 tau = 2 tau (1 - x s*) and delta - c + 2 tau = 2 (1 - c s* + tau s*)/s* = 2 tau (1 - x s*)).
Differentiating (1 - tau) Sph_tau = 1 - x s* gives tau (1 - tau) Sph_tautau = tau Sph_tau - x tau s*_tau, so
R_tr = 2 (1 - tau) Sph_tau - x tau s*_tau - 2 = -x s* (2 + x tau s*/delta) = -x (c + 3 delta)/(delta (c + delta)), R_tr(0) = -2x/c.
[S1_checks.pl C5 at (x, tau) = (2, 5/8): P_tr = 2, R_tr = -9/5.] With Eik = 0, z Exp[-Sph/z] Gauss[Exp[Sph/z] Phi] = P_tr Phi_tau +
R_tr Phi + z T_tr Phi.

**Proposition 3.19.** For Phi = Sum_m z^m Phi_m, the transport equations delta (Phi_m)_tau + R_tr Phi_m + T_tr Phi_{m-1} = 0 (Phi_{-1} := 0) have
unique solutions analytic on D with Phi_0(0) = 1, Phi_m(0) = 0 (m >= 1):
    Phi_0 = Exp[-Int_0^tau R_tr/delta] = 1/alg,    Phi_m(tau) = -Phi_0(tau) Int_0^tau T_tr[Phi_{m-1}](u)/(delta(u) Phi_0(u)) du,
and in Kz[[tau]]
    Ephi F = Phi,   Phi in O(D)[[z]],   [z^0] Phi = 1/alg.                                                                       (3.20)
Proof. On D, delta has no zero and R_tr is analytic (c + delta != 0), which gives existence and uniqueness. Phi_0 = 1/alg: by (3.4),
(log alg)_tau = delta_tau (c + 3 delta)/(2 delta (c + delta)) with delta_tau = -2x/delta, i.e. -x (c + 3 delta)/(delta^2 (c + delta)) = R_tr/delta.
By (3.19), Yt := Exp[dphi/z] iota(Phi) satisfies Gauss[Yt] = 0 and Yt(z, 0) = 1, so Yt = F by Lemma 3.13. QED
Hence Psi_k = [tau^k](Phi alg) lies in C[[z]], which is Corollary 3.16 again, and (3.20) says more: the z-coefficients of Ephi F are
analytic on the whole disc D = {|tau| < |c|^2/(4|x|)}, whose boundary point tau = c^2/(4x) is where s* and s~ coalesce (delta = 0).
Comparison with (3.17): N = Abeta (Phi alg). The Laplace route fixes the constant Abeta, the transport route the tau-dependence.
[S1 referee script referee_scripts/s1ref/ref2.pl R5: the recursion against [z^m](Ephi F) for m <= 3 mod tau^7 at x = 2, 1/3, -5/2.]

### 3.11 First orders; the Euler integral
Closed forms: delta = (c^2 - 4 x tau)^(1/2), s* = 2/(c + delta), U_1 = 1 - tau s*, phi_tau''(s*) = -delta/(x s*^2 U_1^2), W_2(U_1) =
-delta/(x tau^2 s*^2), alg = (delta/c)^(1/2)(c + delta)/(2c), dphi_tau = s*/(1 - tau s*) = (1 - x s*)/(1 - tau), dphi' = 1 - x s*. Series:
1 - x s* = 1/c - x^2 tau/c^3 - 2 x^3 tau^2/c^5 + O(tau^3); with tau = t - t^2/2 + t^3/6 and c^2 - 4x = (1 - x)^2,
dphi(t) = t/c - x^2 t^2/(2c^3) + (1 - x)^2 x^2 t^3/(6c^5) + O(t^4) (as printed in hvc_leaf_proof_check_out.txt); alg = 1 - 2x tau/c^2 + O(tau^2).
Order tau^1, by hand: F_1 = (1 + z)^2/(z (c + 2z)), X_1 = -1/(c z), [tau^1] alg = -2x/c^2, so
    Psi_1 = (1 + z)^2/(z (c + 2z)) - 1/(c z) - 2x/c^2 = (2x + c z)/(c (c + 2z)) - 2x/c^2 = (1 - x)^2 z/(c^2 (c + 2z)):
the poles at z = 0 of F_1 and of Ephi cancel, and [t^1] r_{0,1} = [tau^1] r_{0,1} = (1 - x)^2 z/(c (c + z)(c + 2z)), the O(t) row correction
of leaf_structure.md, Step 0.

**Remark 3.20 (Euler's integral; motivation, used nowhere).** Under (H_E): Re(1/z) > -1, Re(x/z) > -1, |tau| < 1, the integral
J(z, tau) := Int_0^1 s^(1/z) (1 - s)^(x/z) (1 - tau s)^(-1/z-1) ds equals B(1/z + 1, x/z + 1) F(z, tau) (expand (1 - tau s)^-a, a = 1/z + 1,
and use Int_0^1 s^(b-1+n)(1 - s)^(g-b-1) ds = B(b, g - b) (b)_n/(g)_n with b = 1/z + 1, g = c/z + 2). For general a, b, g and
Kk := s^(b-1)(1 - s)^(g-b-1)(1 - tau s)^-a the exact identity
    [tau (1 - tau) d_tau^2 + (g - (a + b + 1) tau) d_tau - a b] Kk = -a d_s [s^b (1 - s)^(g-b) (1 - tau s)^(-a-1)]
holds (after division by a Kk (1 - tau s)^-2 both sides equal -b + g s + (b - a - 1) tau s + (a + 1 - g) tau s^2). Under (H_E) the
boundary terms vanish and Gauss J = 0; the recursion of Lemma 3.13 at numerical z (k + gamma != 0 since Re gamma > 0 under (H_E))
re-proves J = J(z, 0) F. The formal expansion obeys Gauss's equation because of the pointwise identity (3.16) alone, with no boundary
terms and no hypothesis (H_E). This is the precise form of conifold_regularity_note.md s.6, step 1 ("the thimble is the exponent-0
solution of Gauss's equation"), and it completes the F2 sketch of that note's s.1 to all orders in z.

## 4. Metric, eigenvalues, QDE and the symplectic constant at the singular idempotent (Steps 3, 5, 6)

"rc N" refers to line N of rc_lib.m.

### 4.1 Frobenius data at the torus point
**Lemma 4.1.** At alpha0 = 0, alpha1 = e1 = x, alpha2 = e2 = 1, in the flat basis (1, H):
 (a) eta = {{c/x^3, 1/x^2}, {1/x^2, 0}} and eta^-1 = {{0, x^2}, {x^2, -c x}}, independent of t;
 (b) 1 * a = a for all a, and H * H = -(x/tau) 1 + (c/tau) H;
 (c) Leq = L^2 - (c/tau) L + x/tau, Ssum = c/tau, Pprod = x/tau (rc 23-25);
 (d) the product is commutative, associative and Frobenius: eta(a * b, d) = eta(a, b * d).
Proof. (a) rc 16: Replp2 = -Expand[p (p + a0)] + p^2 = -a0 p, so p^2 is reduced by p (p + a0) = 0, whose roots are the fixed points
p = 0, -a0 of P^1. The reduced cube exp3 (rc 17) agrees with (t0 + t1 p)^3 at both roots, so rc 18 is the localisation formula
PrepCl = (1/6)[t0^3/(a0 a1 a2 e1 e2) - (t0 - a0 t1)^3/(a0 (a0 + a1)(a0 + a2) e1 e2)] (cubic, so the Series truncation is exact). Its third
derivatives are PrepCl_000 = (a0 + a1 + a2)/(a1 a2 (a0 + a1)(a0 + a2) e1 e2), PrepCl_001 = 1/((a0 + a1)(a0 + a2) e1 e2),
PrepCl_011 = -a0/((a0 + a1)(a0 + a2) e1 e2), PrepCl_111 = a0^2/((a0 + a1)(a0 + a2) e1 e2) (the entries of etaSym, rc 82); none has a0 in a
denominator, so a0 = 0 may be set. eta_ij = Prep_0ij (rc 19). At the torus point a1 a2 = (a0 + a1)(a0 + a2) = e1 e2 = x and
a0 + a1 + a2 = c, which gives eta; det eta = -1/x^4. (b) Prep = PrepCl + Li_3(e^{t1})/(e1 e2) (rc 21), (d/dt1)^3 Li_3(e^{t1}) =
e^{t1}/(1 - e^{t1}), and with t1 = log Q, Q/(1 - Q) = -1/tau: Prep_111 = Q/((1 - Q) x) = -1/(x tau); the other third derivatives are classical.
cch[i, j, k] = Sum_m eta^{im} Prep_mjk (rc 22) are the structure constants, e_j * e_k = Sum_i cch[i, j, k] e_i (e_0 = 1, e_1 = H). Prep_0jk = eta_jk
gives cch[i, 0, k] = delta_{ik}, and cch[0, 1, 1] = x^2 (-1/(x tau)) = -x/tau, cch[1, 1, 1] = (-c x)(-1/(x tau)) = c/tau. (c) rc 23-25 with (b).
(d) The algebra is generated by H, and eta(e_j * e_k, e_l) = Prep_ljk is totally symmetric. QED

### 4.2 Critical points, eigenvalues, canonical coordinates
U W_U = -x Q U/(1 - Q U) - 1/(1 - U), W_2 = -x Q U/(1 - Q U)^2 - U/(1 - U)^2, and (partial_t W)(U) = -x Q U/(1 - Q U).
**Lemma 4.2.** (a) U W_U = 0 iff x Q U^2 + Q (1 - x) U - 1 = 0 (LeqU, rc 115).
 (b) At a critical point (partial_t W)(U) = 1/(1 - U) = (x U + 1)/tau = L(U).
 (c) With U = 1 - tau s (tau != 0), U is critical iff x tau s^2 - c s + 1 = 0; the critical points are U_1 = 1 - tau s* and U_2 = 1 - tau s~.
 (d) L_1 = 1/(tau s*) = c/tau - x s*, L_2 = x s* = 1/(tau s~), L_1 + L_2 = c/tau, L_1 L_2 = x/tau, L_1 - L_2 = delta/tau. L_1, L_2 are the
     roots of Leq, the eigenvalues of H*. L_1 = c/t + l(t) with l := c (1/tau - 1/t) - x s* in O, l(0) = c/2 - x/c = (1 + x^2)/(2c),
     l_1 = c/12 - x^2/c^3; L_2 is in O with L_2(0) = x/c.
 (e) u_i := W(U_i) (any branch) satisfies du_i/dt = L_i. Explicitly u_1 = x t + c log tau + phi_tau(s*) + x log(-1), so with
     mu = c/tau - 1 = x + c (1 - tau)/tau: L_1 = mu + dphi', L_1 - dphi' = c/tau - 1; and u_2 = x t + c log tau + phi_tau(s~) + x log(-1),
     L_2 = mu + (1 - x s~) (d/dt phi_tau(s~) = 1 - x s~ = 1 - L_1, Lemma 3.1(b)(iii) at s~).
Proof. (a) U W_U = -[x Q U (1 - U) + (1 - Q U)]/((1 - Q U)(1 - U)), and x Q U (1 - U) + 1 - Q U = -(x Q U^2 + Q (1 - x) U - 1).
(b) U W_U = 0 gives -x Q U/(1 - Q U) = 1/(1 - U). By (a), x U^2 = 1/Q - (1 - x) U, so (x U + 1)(1 - U) = 1 + (x - 1) U - x U^2 = 1 - 1/Q = tau;
and -Q/(1 - Q) = 1/tau. (c) 1 - U = tau s and 1 - Q U = -Q tau (1 - s); (x U + 1)(1 - U) = tau becomes (c - x tau s) tau s = tau.
(d) By (b), L(U) = 1/(1 - U) = 1/(tau s) at s = s*, s~. So L_2 = 1/(tau s~) = x s*, and L_1 = 1/(tau s*) = c/tau - x s* is the critical equation
divided by tau s*. L_1 + L_2 = (s* + s~)/(tau s* s~) = c/tau, L_1 L_2 = 1/(tau^2 s* s~) = x/tau, L_1 - L_2 = (c - 2 x tau s*)/tau = delta/tau. By
Lemma 4.1(b) the minimal polynomial of H* is Leq. 1/tau = 1/t + 1/2 + t/12 + O(t^3) and x s* = x/c + x^2 t/c^3 + O(t^2) give l(0), l_1.
(e) u_i' = (partial_t W)(U_i) + W_U(U_i) U_i' = L_i by (b). The explicit form is (3.1) with log Q = t; then d/dt phi_tau(s*) = dphi' by the
envelope formula, and x + c (1 - tau)/tau + dphi' = c/tau - 1 + 1 - x s* = L_1. At s~ the same envelope computation gives
d/dt phi_tau(s~) = (1 - tau) s~/(1 - tau s~) = 1 - x s~, and x s~ = 1/(tau s*) = L_1, so mu + 1 - x s~ = c/tau - L_1 = L_2. QED

### 4.3 Idempotents and the metric: (V1), (V2), F1
**Proposition 4.3.** (a) The columns of Psi^-1, Psi = {{1, L1}, {1, L2}} (rc 29), are the idempotents e_1 = (H - L_2)/(L_1 - L_2),
e_2 = (H - L_1)/(L_2 - L_1): e_i * e_j = delta_{ij} e_i, H * e_i = L_i e_i, e_1 + e_2 = 1.
 (b) etacan = Psi^-T eta Psi^-1 = diag(eta(e_1, e_1), eta(e_2, e_2)), and the code's Delta (DeltaL rc 31, DeltaCN rc 128-129) is
     Delta_i = 1/eta(e_i, e_i) at U_i (DeltaCN[2] = cconj[DeltaCN[1]]).
 (c) For {i, j} = {1, 2}:
        Delta_i = x^3 (L_i - L_j)/(x - c L_j) = x L_i (2x - c L_i) = (L_i - L_j)(x^2 - c x L_i) = x^2 W_2(U_i).                    (4.1)
     In particular Delta_1 = x^2 W_2(U_1) is an identity in the quadratic algebra QA (V1), and so is Delta_2 = x^2 W_2(U_2).
 (d) Delta_1 = -x delta/(tau^2 s*^2), Delta_2 = x^3 delta s*^2, x - c L_2 = -x^2 tau s*^2, x - c L_1 = -1/(tau s*^2), (x - c L_1)(x - c L_2) = x^2,
     Delta_1 Delta_2 = -x^4 delta^2/tau^2, 1/Delta_1 + 1/Delta_2 = c/x^3 = eta(1, 1).
 (e) F1: Delta_2 is in O with Delta_2 = x^3/c - x^5 tau^2/c^5 + O(tau^3), so Delta_2(0) = x^3/c != 0 and Delta_2'(0) = 0; Delta_1 is tau^-2 times a
     unit of O, Delta_1 = -x c^3 t^-2 (1 + O(t)) (Proposition 4.4).
Proof. (a) Psi^-1 = (L_2 - L_1)^-1 {{L_2, -L_1}, {-1, 1}}; its first column (L_2, -1)/(L_2 - L_1) is (L_2 1 - H)/(L_2 - L_1) = e_1, the second e_2.
By Lemma 4.1(b) and 4.2(d), (H - L_1) * (H - L_2) = 0, hence H * e_1 = L_1 e_1 and e_1 * e_2 = 0; with e_1 + e_2 = 1, e_i * e_i = e_i.
(b) The (i, j) entry of Psi^-T eta Psi^-1 is eta(e_i, e_j) = eta(e_i * e_j, 1) = 0 for i != j.
(c) First expression (Frobenius route): eta(e_i, e_i) = eta(e_i * e_i, 1) = eta(e_i, 1) = (eta(H, 1) - L_j eta(1, 1))/(L_i - L_j) =
(x - c L_j)/(x^3 (L_i - L_j)). The code evaluates the quadratic form literally (rc 30-31: L1, L2 independent, then L2 -> Ssum - L1, then
reduction modulo Leq): e_1^T eta e_1 = L_2 (c L_2 - 2x)/(x^3 (L_1 - L_2)^2). The two agree on the critical locus, since
L_2 (c L_2 - 2x) - (x - c L_2)(L_1 - L_2) = c L_1 L_2 - x (L_1 + L_2) = c x/tau - x c/tau = 0. Second expression: L_i (2x - c L_i)(x - c L_j) =
L_i [2x^2 - c x L_j - c x (L_i + L_j) + c^2 L_i L_j] = L_i (2x^2 - c x L_j), because c x (L_i + L_j) = c^2 x/tau = c^2 L_i L_j, and
L_i (2x^2 - c x L_j) = 2x^2 L_i - x^2 (L_i + L_j) = x^2 (L_i - L_j). Third expression: (x - c L_1)(x - c L_2) = x^2 - c x (L_1 + L_2) + c^2 L_1 L_2 = x^2,
so x^3/(x - c L_j) = x (x - c L_i). Fourth expression: at a critical point -x Q U/(1 - Q U) = L (Lemma 4.2(b)), so 1/(1 - Q U) = (x - L)/x and
-x Q U/(1 - Q U)^2 = L (x - L)/x; and U = 1 - 1/L gives U/(1 - U)^2 = L (L - 1). Therefore W_2(U) = L (x - L)/x - L (L - 1) = L (2x - c L)/x.
QA has the two characters ev_1, ev_2, jointly injective; DeltaCN[1] = toC[DeltaL(LofU)] (rc 128-129) and toC[x^2 W_2] have equal values at
U_1 and at U_2, hence are equal (check (1) of hvc_leaf_proof_check.wls).
(d) With L_1 = 1/(tau s*) and 2 x tau s* - c = -delta: Delta_1 = x (2 x tau s* - c)/(tau^2 s*^2) = -x delta/(tau^2 s*^2). With L_2 = x s* and
2 - c s* = 2 delta/(c + delta) = delta s*: Delta_2 = x^3 s* (2 - c s*) = x^3 delta s*^2. x - c L_2 = x (1 - c s*) = -x^2 tau s*^2 (Lemma 3.1(b)(v));
x - c L_1 = (x tau s* - c)/(tau s*) with x tau s* - c = (c - delta)/2 - c = -1/s*. The rest follows, and 1/Delta_1 + 1/Delta_2 =
eta(e_1 + e_2, 1) = c/x^3.
(e) With u := 4 x tau/c^2, delta/c = 1 - u/2 - u^2/8 + O(u^3) and c s* = 2/(1 + delta/c) = 1 + u/4 + u^2/8 + O(u^3), so
delta s*^2 = (1/c)(1 - u^2/16 + O(u^3)) = 1/c - x^2 tau^2/c^5 + O(tau^3). QED

### 4.4 Step 3: the Hessian in the Euler coordinate; rho; the KEY FORMULA
**Proposition 4.4.** (a) Delta_1 = x^2 W_2(U_1) = x^2 U_1^2 phi_tau''(s*)/tau^2 = -x c^3 alg^2/tau^2 (by (4.1), (3.5), (3.4)).
 (b) rho = tau/(t alg), a unit of O with rho(0) = 1 (the branch of (C13)), and
        t^2 rho^2 Delta_1 = -x c^3   exactly;                                                                                  (4.2)
     Delta_1 t^2/(-x c^3) = (t/tau)^2 alg^2 is a unit of O with value 1 at t = 0 (F1 for Delta_1).
 (c) KEY FORMULA (V3). With (3.18),
        G(z; t) := rho r_{0,1}(z; t) = (c/(c + z)) (tau/t) F(z, tau) Ephi(z, tau);                                                (4.3)
     G is in C[[z, t]] (Theorem 3.17), G(z; 0) = c/(c + z), and every G_k = [t^k] G is a rational function of z without pole at z = 0, whose
     Taylor series is [t^k](rho r_{0,1}).
 (d) rho = 1 - (1 - x)^2 t/(2c^2) + (1 - 8x + 15x^2 - 8x^3 + x^4) t^2/(6c^4) + O(t^3) and Delta_1 = -(x c^3/t^2)(1 + (1 - x)^2 t/c^2 + O(t^2)),
     as in leaf_structure.md (Step 1 and Step 0).
Proof. (a), (b): rho^-2 = Delta_1 t^2/(-x c^3) = (t alg/tau)^2, and its square root with value 1 at t = 0 is t alg/tau. (c) rho r_{0,1} =
(tau/(t alg)) (c/(c + z)) F Ephi alg. The rational form of G_k: by (b), G = Sum_j rho (c/(c + z)) Psi_j tau^j, so G_k = Sum_{m<=k} ([t^m] rho)
Sum_{j<=k-m} ([t^(k-m)] tau^j) (c/(c + z)) Psi_j, a finite combination of functions regular at z = 0 (Corollary 3.16) whose Taylor series is
Sum_m ([t^m] rho) [t^(k-m)] r_{0,1} = [t^k](rho r_{0,1}) by Theorem 3.17. G_0 = c/(c + z) since Psi_0 = 1. (d) tau = t - t^2/2 + O(t^3),
alg = 1 - 2x tau/c^2 + O(tau^2), so rho = (tau/t)/alg = 1 + (2x/c^2 - 1/2) t + O(t^2) = 1 - (c^2 - 4x) t/(2c^2) + O(t^2), c^2 - 4x = (1 - x)^2;
the t^2 term follows from the branch-free form of Corollary 11.2(d). QED
[Mathematica (V2), (V3), checks (2b), (2c), (3), (3b); Perl: S2_checks.pl part E (rho and Delta_1 against leaf_structure.md),
S4_checks.pl ("rho: three expressions"), and the spot check referee_scripts/assemble/spot.pl.]

### 4.5 The operator identity (V7)
**Lemma 4.5.** Let Bop := C[z, 1/z][tau, 1/tau, 1/(1 - tau)], and consider the ring of differential operators generated by Bop and one symbol
Dd subject to Dd a = a Dd + Dd(a), with Dd(a) := z a' (so Dd(tau) = z (1 - tau)), and put d/dtau := (z (1 - tau))^-1 Dd. Then, with omega as in
s.1.5 and Gauss written with this d/dtau,
    (Dd + omega)^2 - (c/tau)(Dd + omega) + x/tau = ((1 - tau) z^2/tau) Gauss.                                                     (4.4)
In particular (Dd = Dz) for every f in Rz: [(Dz + omega)^2 - (c/tau)(Dz + omega) + x/tau] f = ((1 - tau) z^2/tau) Gauss[f]. For any other operator
Dd' that satisfies the same Leibniz rule, Dd'(a Y) = z a' Y + a Dd' Y, (4.4) holds with Dd' in place of Dd and (z (1 - tau))^-1 Dd' in place of
d/dtau. (V7) is the case Dd = Dz conjugated by Pf := (1 - tau)^(-x/z) tau^(c/z+1): formally Dz(Pf f) = Pf (Dz + omega) f, since
omega = mu + Dz(tau)/tau with mu = x + c (1 - tau)/tau the logarithmic t-derivative of (1 - tau)^(-x/z) tau^(c/z) times z.
Proof. (Dd + omega)^2 f = Dd^2 f + 2 omega Dd f + Dd(omega) f + omega^2 f, so the left side applied to f is
Dd^2 f + (2 omega - c/tau) Dd f + (Dd(omega) + omega^2 - (c/tau) omega + x/tau) f, with Dd f = z (1 - tau) f_tau and
Dd^2 f = z^2 (1 - tau)^2 f_tautau - z^2 (1 - tau) f_tau (f_tau := d/dtau f). Coefficients:
 * of f: Dd(omega) = -z (1 - tau)(c + z)/tau^2 and omega (omega - c/tau) = z (c + z)/tau^2 - (1 + z)(c + 2z)/tau + (1 + z)^2. The tau^-2 terms add
   up to z (c + z)/tau, and z (c + z) - (1 + z)(c + 2z) + x = -(1 + z)^2, so the coefficient is -(1 + z)^2/tau + (1 + z)^2 = -(1 - tau)(1 + z)^2/tau;
 * of f_tau: -z^2 (1 - tau) + z (1 - tau)[(c + 2z)/tau - 2 (1 + z)] = z (1 - tau)[(c + 2z)/tau - (2 + 3z)];
 * of f_tautau: z^2 (1 - tau)^2.
On the right, ((1 - tau) z^2/tau) Gauss[f] has f_tautau-coefficient z^2 (1 - tau)^2, f_tau-coefficient (1 - tau) z^2 [(c/z + 2)/tau - (2/z + 3)]
= z (1 - tau)[(c + 2z)/tau - (2 + 3z)], and f-coefficient -((1 - tau) z^2/tau)(1 + 1/z)^2 = -(1 - tau)(1 + z)^2/tau. Only the Leibniz rule and
Dd(tau) = z (1 - tau) were used, which gives the statement for Dd'. For Pf: z (log Pf)' = x + c (1 - tau)/tau + z (1 - tau)/tau = omega. QED
[hvc_pf_gauss.wls (Mathematica, symbolic F); S2_checks.pl part B (exactly on tau^n, -3 <= n <= 7, and on a mixed Laurent polynomial, for
four x and four z).]

### 4.6 Step 5: the second row of the QDE at the singular idempotent
**Proposition 4.6.** In Rz0 (r_{0,1} in C[[z, tau]], r_{1,1} in tau^-1 C[[z, tau]]), hence by TRANSFER to all orders in z as identities of
meromorphic functions of t near 0:
    Row 1:  r_{1,1} = (L_1 + z D_1) r_{0,1}    (the definition (H));
    Row 2:  (L_1 + z D_1) r_{1,1} = -(x/tau) r_{0,1} + (c/tau) r_{1,1},  equivalently  z D_1 r_{1,1} = L_2 (r_{1,1} - L_1 r_{0,1}),
i.e. (D_1 r_{1,1})[m-1] = L_2 (r_{1,1}[m] - L_1 r_{0,1}[m]) for every m >= 1 (the form tested in V5). In matrix form z D_1 v_1 = (A - L_1) v_1, and
for vn := t rho v_1, an element of C[[z, t]]^2,
    z vn' = (A - L_1) vn.                                                                                                        (4.5)
Proof. By (4.2), Delta_1'/(2 Delta_1) = -(t rho)'/(t rho), so D_1 f = (t rho)^-1 (t rho f)', i.e.
    L_1 + z D_1 = (t rho)^-1 (z d/dt + L_1)(t rho)    as operators on Rz.                                                            (4.6)
(i) Conjugation. For f in Rz, using z Ephi' = -dphi' Ephi, tau' = 1 - tau and Lemma 4.2(e),
    (z d/dt + L_1)(tau Ephi f) = tau Ephi [z f' + (L_1 - dphi' + z (1 - tau)/tau) f] = tau Ephi (Dz + omega) f,
because L_1 - dphi' + z (1 - tau)/tau = c/tau - 1 + z/tau - z = omega.
(ii) By (4.3), t rho r_{0,1} = t G = (c/(c + z)) tau Ephi F. By (4.6) and (i),
    t rho r_{1,1} = (z d/dt + L_1)(t G) = (c/(c + z)) tau Ephi (Dz + omega) F,
    t rho (L_1 + z D_1) r_{1,1} = (z d/dt + L_1)(t rho r_{1,1}) = (c/(c + z)) tau Ephi (Dz + omega)^2 F.
(iii) Hence, by (4.4) and Gauss[F] = 0 (Lemma 3.13),
    t rho [(L_1 + z D_1) r_{1,1} + (x/tau) r_{0,1} - (c/tau) r_{1,1}] = (c/(c + z)) tau Ephi [(Dz + omega)^2 - (c/tau)(Dz + omega) + x/tau] F
                                                                 = (c/(c + z)) Ephi (1 - tau) z^2 Gauss[F] = 0,
and t rho is a unit of Rz. The second form uses c/tau - L_1 = L_2 and x/tau = L_1 L_2. The components of (A - L_1) v_1 are r_{1,1} - L_1 r_{0,1}
= z D_1 r_{0,1} (Row 1) and -(x/tau) r_{0,1} + (c/tau - L_1) r_{1,1} = z D_1 r_{1,1} (Row 2), and z (t rho f)' = t rho z D_1 f, which gives (4.5).
Membership: L_1 and Delta_1'/Delta_1 lie in tau^-1 O, so L_1 + z D_1 maps C[[z, tau]] into tau^-1 C[[z, tau]] and tau^-1 C[[z, tau]] into
tau^-2 C[[z, tau]]; with Theorem 3.17 every term lies in Rz0. vn = (t G, (t L_1 + z) G + z t G') lies in C[[z, t]]^2, since t L_1 is in O. QED
Given Row 1 and (4.3), Row 2 is equivalent to Gauss[F] = 0, since (c/(c + z)) Ephi (1 - tau) z^2 is a unit of Rz. No property of the Wick
series beyond Theorem 3.17 is used, and nothing about the regular idempotent. From (ii), since F is a unit of Kz[[tau]],
    r_{1,1} = r_{0,1} [(c + z)/tau - (1 + z) + z (1 - tau) F_tau/F],    equivalently    rho r_{1,1} = (L_1 + z/t) G + z G'           (4.7)
(the two agree: z G'/G = z (1 - tau)/tau - z/t + z (1 - tau) F_tau/F - dphi', and L_1 - dphi' = c/tau - 1). In F-form, with
PsiF1 := tau (Dz + omega) F = z tau (1 - tau) F_tau + (c + z - (1 + z) tau) F, Row 2 reads z (1 - tau) (PsiF1)_tau = PsiF1 - x F, because
z (1 - tau) (PsiF1)_tau - PsiF1 + x F = (1 - tau) z^2 Gauss[F]. [Mathematica (V5) for i = 1, m = 1..9, t^0..t^6; S2_checks.pl part C (sigma
form, t^22) and D (F-form, tau^22); the referee script referee_scripts/ref_S2_jets.pl from the code's rows at generic t0 != 0.]

### 4.7 The operator A
**Proposition 4.7.** (a) A is the transpose of the matrix M_H = {{0, -x/tau}, {1, c/tau}} of H* in the basis (1, H) (columns: coordinates of
H * 1, H * H); A acts on covectors by s -> s(H * .).
 (b) chi_i = (1, L_i) (the character a_0 1 + a_1 H -> a_0 + a_1 L_i; chi_i(e_j) = delta_{ij}) satisfies A chi_i = L_i chi_i.
 (c) A is self-adjoint for < , >: eta^-1 A = {{-x^3/tau, c x^2/tau}, {c x^2/tau, x^2 - c^2 x/tau}} is symmetric, so <A a, b> = <a, A b>.
 (d) <chi_i, chi_j> = delta_{ij} Delta_i.
 (e) A - L_1 = (L_2 - L_1) Pr_2 with Pr_2 a := (pi(a)/Delta_2) chi_2 the < , >-orthogonal projector onto C chi_2 along C chi_1, where
        pi(a) := <a, chi_2> = x^2 L_2 a_0 + (x^2 - c x L_2) a_1 = x^3 s* (a_0 - tau s* a_1);
     and, for all covectors a, b,
        <(A - L_1) a, b> = ((L_2 - L_1)/Delta_2) pi(a) pi(b) = -(x^3/tau)(a_0 - tau s* a_1)(b_0 - tau s* b_1).                     (4.8)
Proof. (a) H * 1 = H and H * H = -(x/tau) 1 + (c/tau) H; (A s)_0 = s_1 = s(H * 1), (A s)_1 = -(x/tau) s_0 + (c/tau) s_1 = s(H * H).
(b) A (1, L) = (L, -x/tau + (c/tau) L), which is L (1, L) iff Leq(L) = 0. (c) Multiply out; conceptually H* is eta-self-adjoint (Lemma 4.1(d)).
(d) <chi_i, chi_j> = x^2 (L_i + L_j) - c x L_i L_j: for i != j this is x^2 c/tau - c x x/tau = 0; for i = j it is 2x^2 L_i - c x L_i^2 = Delta_i
by (4.1). (e) chi_1, chi_2 is an orthogonal basis with non-zero norms, so a = Sum_i (<a, chi_i>/Delta_i) chi_i, 1 = Pr_1 + Pr_2 and
A = L_1 Pr_1 + L_2 Pr_2; <Pr_2 a, b> = pi(a) pi(b)/Delta_2. x^2 L_2 = x^3 s* and x^2 - c x L_2 = -x^3 tau s*^2 (Proposition 4.3(d)) give pi;
(L_2 - L_1)/Delta_2 = (-delta/tau)/(x^3 delta s*^2) = -1/(x^3 tau s*^2) gives the last form. QED
Note pi(chi_1) = x^3 s* (1 - tau s* L_1) = 0 is literally L_1 = 1/(tau s*). Self-adjointness is K((t))-linear and (4.8) bilinear; both extend
coefficientwise to covectors with entries in series rings in psi, psi'. [S2_checks.pl part A.]

### 4.8 Step 6: the symplectic constant K_1 = 1
**Proposition 4.8.** (a) K_1(z) := <v_1(z), v_1(-z)>/Delta_1 = -(1/(x c^3)) <vn(z), vn(-z)> lies in C[[z, t]] and is independent of t: K_1 is in C[[z]].
 (b) sigma(psi) := rho v_1(-psi) = (g, (L_1 - psi/t) g - psi g'), with g(psi) = G(-psi) and ' = d/dt at fixed psi. t sigma(psi) is in C[[psi, t]]^2,
        sigma^(-1) = [t^-1] sigma(psi) = (0, c),     g(psi)|_{t=0} = c/(c - psi) = c lam_1(psi)   (F2).
 (c) K_1 = 1, i.e. <v_1(w), v_1(-w)> = Delta_1 to all orders in w ((V6) for i = j = 1).
 (d) Consequences: <sigma(psi), sigma(-psi)> = rho^2 Delta_1 = -x c^3/t^2; in (1.1) the subtracted term is P_11(-w, w) = -Delta_1, so
        rho^2 E_11(psi, psi') = [ -x c^3/t^2 - <sigma(psi), sigma(psi')> ]/(psi + psi');                                           (4.9)
     and -psi d/dt (t sigma(psi)) = (A - L_1)(t sigma(psi)), i.e. d sigma/dt = -sigma/t - (1/psi)(A - L_1) sigma.
Proof. (a) (t rho)^2 = -x c^3/Delta_1 gives the second expression, which lies in C[[z, t]] by Proposition 4.6. Applying refl to (4.5) gives
-z vn(-z)' = (A - L_1) vn(-z) (A and L_1 do not involve z). eta^-1 is constant, so
    z <vn(z), vn(-z)>' = <(A - L_1) vn(z), vn(-z)> - <vn(z), (A - L_1) vn(-z)> = 0
by Proposition 4.7(c). Thus (1 - tau) dK_1/dtau = 0, and the kernel of d/dtau on C[[z, tau]] is C[[z]].
(b) By (C9), (4.6) and (4.3), rho r_{1,1} = rho (t rho)^-1 (z d/dt + L_1)(t G) = (L_1 + z/t) G + z G'. By Lemma 4.2(d), t L_1 is in O with value c,
so t rho r_{1,1} = (t L_1 + z) G + z t G' is in C[[z, t]], and at t = 0 it equals (c + z) c/(c + z) = c: the z-dependence cancels exactly because
of the geometric form c/(c + z) of F2. Substitute z = -psi.
(c) vn(z) = t rho v_1(z) = t sigma(-z), so K_1(z) = -(1/(x c^3)) <t sigma(-z), t sigma(z)>. By (b) this lies in C[[z, t]] with value at t = 0
equal to -(1/(x c^3)) <(0, c), (0, c)> = -(1/(x c^3)) c^2 (-c x) = 1. By (a), K_1 = 1.
(d) <sigma(psi), sigma(-psi)> = rho^2 Delta_1 K_1(-psi), and P_11(-w, w) = -<v_1(w), v_1(-w)> = -Delta_1 K_1(w); with z = psi, w = psi' and
rho^2 P_11(psi, psi') = -<sigma(psi), sigma(psi')> this gives (4.9). Finally (4.5) at z = -psi, with vn(-psi) = t sigma(psi), is
-psi (t sigma(psi))' = (A - L_1)(t sigma(psi)); divide by -psi t. The division by psi is legitimate: v_1(0) = r_{0,1}[0] chi_1 by (H), and
(A - L_1) chi_1 = 0, so (A - L_1) sigma vanishes at psi = 0. QED
Only the t -> 0 data F2 and sigma^(-1) = (0, c) enter the limit in (c), and they come from (H) and (4.3) alone; the t-independence (a) uses
Row 2, i.e. the full closed form. [Mathematica (V6) for i = j = 1, w^0..w^9; S2_checks.pl part C (t^2 <sigma(psi), sigma(-psi)> = -x c^3
to t^22); referee scripts referee_scripts/ref_S2_key.pl, ref_S2_leaves.pl.]

**Remark 4.9 (K_1 = 1 as a Wronskian identity for Gauss).** By Proposition 4.6(ii), vn = (c/(c + z)) tau Ephi PhiF with PhiF := (F, (Dz + omega) F).
Since Ephi(z) Ephi(-z) = 1 and (c/(c + z))(c/(c - z)) = c^2/(c^2 - z^2),
    K_1(z) = -tau^2 <PhiF(z), PhiF(-z)>/(x c (c^2 - z^2)),    so    K_1 = 1  iff  tau^2 <PhiF(z), PhiF(-z)> = -x c (c^2 - z^2).
The latter follows from Gauss[F] = 0 alone: (Dz + omega) PhiF = A PhiF (Row 2 in F-form), refl gives (-Dz + omega(-z)) PhiF(-z) = A PhiF(-z), so
d/dt <PhiF(z), PhiF(-z)> = ((omega(-z) - omega(z))/z) <PhiF(z), PhiF(-z)> = -2 ((1 - tau)/tau) <PhiF(z), PhiF(-z)>, and d/dt (tau^2 <...>) = 0; at
tau = 0, tau PhiF(z) -> (0, c + z), which gives -c x (c + z)(c - z). Given (4.3) this is a second proof of Proposition 4.8(c). [S2_checks.pl
part D, exact to tau^22.]

**Remark 4.10 (the flat sections of the outline, and the constants kappa).** For lambda in C((tau)) let Mod_lambda := Rz eps_lambda be free
of rank one with nabla(f eps_lambda) := (z f' + lambda f) eps_lambda, and write f e^{u/z} := f eps_{u'} for any (possibly multivalued) u with
u' = lambda; so e^{u/z} depends only on u', and u, u + const define the same symbol. Products eps_lambda eps_mu := eps_{lambda + mu} satisfy the
Leibniz rule for nabla; if gamma0 is a unit of Rz with z gamma0' = (lambda - mu) gamma0, then f eps_lambda -> gamma0 f eps_mu is an isomorphism
(this is the meaning of e^{u/z} = gamma0 e^{v/z}); refl(f eps_lambda) := refl(f) eps_{-lambda}. Q^{x/z} = e^{x t/z} and e^{+-dphi/z} are honest
units of Rz; the only genuinely formal factor in the outline's s_0 is tau^{c/z} = e^{c log(tau)/z}. eps_mu (mu = c/tau - 1) is the symbol
"Q^{x/z} tau^{c/z}". Fix kappa_0 with kappa_0^2 = -x c^3 and put nu_1 := t rho/kappa_0, so nu_1^2 = 1/Delta_1. The outline's
s = e^{u_1/z} Delta_1^(-1/2) v_1 is Sm := nu_1 v_1 eps_{L_1}; Proposition 4.6 says nabla Sm = A Sm ("z ds/dt = A s"). Since L_1 = mu + dphi',
Mod_{L_1} ~= Mod_mu via f eps_{L_1} -> e^{dphi/z} f eps_mu, and Sm_0 -> kappa_1(z) tau F eps_mu with kappa_1(z) := c/(kappa_0 (c + z)): this is
"s_0 = kappa(z) Q^{x/z} tau^{c/z+1} F" with kappa = kappa_1 in frame B and h_w kappa_1 in the point frame, t-independent. Read with honest
exponentials (u_1 = phi_0(1/c) + dphi + c log tau + x log Q + x log(-1), from (3.1)), s_0 = e^{(phi_0(1/c) + x log(-1))/z} kappa_1(z)
Q^{x/z} tau^{c/z+1} F = -kappa_ph(z) e^{phi_0(1/c)/z} kappa_1(z) Q^{x/z} tau^{c/z+1} F in frame B; the sign -1 is the orientation factor of
dU/U = -tau ds/(1 - tau s), which belongs to the integrand (3.2) and not to s_0. The additive constants of u_1 never enter any proof.
Replacing r_{0,1} by kappa r_{0,1} with kappa a unit of Kz replaces r_{1,1} by kappa r_{1,1}, preserves Row 2 and replaces K_1(z) by
kappa(z) kappa(-z) K_1(z); for the frame change kappa = h_w, h_w(z) h_w(-z) = 1, so K_1 and the symplectic relations are the same in both
frames. Exponents (a remark, used nowhere): for lambda = a/tau + O(1) and f of tau-valuation r with leading coefficient f_r, the operator
nabla^2 - (c/tau) nabla + x/tau maps f eps_lambda to (z ex (z (ex - 1) - c) f_r tau^(r-2) + O(tau^(r-1))) eps_lambda, ex := r + a/z; so solutions of this second-order
equation have ex = 0 or ex = 1 + c/z (on Rz the indicial equation z^2 r (r - 1) - c z r = 0 of the Picard-Fuchs operator). Sm_0 has
ex = 1 + c/z (a = c, r = 1); the regular flat section Sm^(2) of Remark 8.10 has ex = 0. The singular thimble carries the exponent c/z + 1,
the regular one the exponent 0, as anticipated in conifold_regularity_note.md s.6.

## 5. Step 4: the dilaton leaves (L1)

Sections 5-7 are written over Kx = Q(x) (see s.1.4 (P)); every function of t is an element of Kx((t)), the expansion at t = 0 of the
corresponding germ. Products of series in Kx((t))[[psi]] or Kx((t))[[psi, psi']] are expanded order by order with finite sums (all of them
lie in t^-2 Kx[[t]][[psi, psi']]); elements of Kx(z)[[t]] (such as G and z G) are compared with series in z or psi only through
Proposition 4.4(c): G_k is regular at z = 0 with Taylor series [t^k](rho r_{0,1}), and g_k is the Taylor series of G_k(-psi).

### 5.1 The tau-coefficients
**Lemma 5.1.** Let n, m >= 0.
 (a) F_n has poles only at z = 0 (order exactly n) and at the points z = -c/j, 2 <= j <= n + 1 (order <= 1; exactly 1 over Q(x), where c is not
     rational), and F_n(oo) = 1/(n + 1).
 (b) X_m = [tau^m] Ephi is a polynomial in 1/z of degree <= m with X_m(oo) = delta_{m0}.
 (c) H_n = Sum_{n1+n2=n} F_{n1} X_{n2} has poles only at z = 0 (order <= n) and at z = -c/j, 2 <= j <= n + 1 (order <= 1), and H_n(oo) = 1/(n + 1).
Proof. (a) The numerator prod_{j=1}^{n} (1 + j z)^2 has degree 2n, value 1 at 0 and leading coefficient (n!)^2. The denominator
n! z^n prod_{j=2}^{n+1} (c + j z) has degree 2n, leading coefficient n! (n + 1)!, a zero of order n at 0 and simple zeros at the n distinct
non-zero points -c/j. A zero of the numerator can only lower a pole order (the zero -1/i of 1 + i z meets -c/j iff c = j/i). Hence the poles
are as stated and F_n(oo) = (n!)^2/(n! (n + 1)!) = 1/(n + 1). (b) dphi in tau Kx[[tau]] gives dphi^l in tau^l Kx[[tau]], so
X_m = Sum_{l=0}^{m} ((-1)^l/l!) z^-l [tau^m](dphi^l), a polynomial in 1/z of degree <= m with constant term [tau^m] dphi^0 = delta_{m0}.
(c) X_{n2} is regular away from 0 with a pole of order <= n2 at 0; at oo, F_{n1} X_{n2} -> delta_{n2,0}/(n1 + 1). QED

### 5.2 The t-coefficients of z G
**Lemma 5.2.** For k >= 0 put R_k(z) := [t^k](z G) = z G_k(z). Then
    R_k = (c z/(c + z)) Sum_{n=0}^{k} d_{k,n} H_n,      d_{k,n} := [t^k]((tau/t) tau^n) in Q,   d_{k,k} = 1.                           (5.1)
So R_k is a rational function with poles only at z = 0 and at z = -c/j, 1 <= j <= k + 1, of order <= 1 at each -c/j; it is bounded at oo and
    R_k(oo) = c delta_{k0}.                                                                                                     (5.2)
Proof. z G = (c z/(c + z)) (tau/t) Sum_n H_n tau^n, with tau/t in 1 + t Q[[t]] and tau^n in t^n (1 + t Q[[t]]); this is (T) and gives the finite
sum (5.1). Poles: by Lemma 5.1(c) each H_n with n <= k has poles only at 0 and at -c/j with 2 <= j <= n + 1 <= k + 1, simple at the latter; the
factor c z/(c + z) adds a simple pole at -c (this is j = 1), which is not one of the -c/j with j >= 2. Limit: c z/(c + z) -> c and H_n -> 1/(n + 1),
so R_k(oo) = c Sum_{n<=k} d_{k,n}/(n + 1) = c [t^k]((tau/t) Sum_{n>=0} tau^n/(n + 1)) = c [t^k](-log(1 - tau)/t) = c [t^k] 1 = c delta_{k0},
since the terms n > k do not contribute to [t^k] and -log(1 - tau) = t. QED
The limit is coefficientwise: F(oo, tau) := Sum_n F_n(oo) tau^n = -log(1 - tau)/tau is the limit of F coefficient by coefficient in tau, and
only this is used; no statement about 2F1(1/z + 1, 1/z + 1; c/z + 2; tau) as z -> oo at fixed tau is needed or claimed. For k >= 1 the limit
vanishes because (tau/t) F(oo, tau) = 1 is t-independent: this is the precise sense in which the normalisation by rho removes lam_0.
Equivalently lim_{z -> oo} z r_{0,1}(z; t) = c/rho(t), coefficientwise in t.

### 5.3 No pole at z = 0
By Proposition 4.4(c), G_k is regular at z = 0, so R_k = z G_k is regular at 0 with R_k(0) = 0. In (5.1) the separate terms do have poles at 0
(z H_n can have a pole of order n - 1); that they cancel is the regularity theorem of section 3. Explicitly every H_n is regular at z = 0, for
either of two reasons: (a) Corollary 3.16: H_n = Sum_{m<=n} ([tau^(n-m)](1/alg)) Psi_m with every Psi_m regular at 0 (no induction is
needed); (b) Proposition 3.19: F Ephi = Phi has z-coefficients analytic in tau, so the Laurent expansion at z = 0 of its tau^n-coefficient H_n
has no negative powers (purely formal: neither the Laplace identification nor the value of Abeta is used). Conversely the regularity of
all H_n gives that of all G_k by (T), since c/(c + z) is a unit at z = 0 and tau/t is a z-independent unit of Kx[[t]]: "every G_k is regular at
0" and "every H_n is regular at 0" are equivalent. By hand: H_1 = F_1 + X_1 = (1 + z)^2/(z (c + 2z)) - 1/(c z) = (2x + c z)/(c (c + 2z)); for
H_2 = F_2 + F_1 X_1 + X_2 (X_2 = -(1 + 2x)/(2 c^3 z) + 1/(2 c^2 z^2)), the z^-2 parts of F_2, F_1 X_1, X_2 are 1/(2c^2), -1/c^2, 1/(2c^2) and the
z^-1 parts are (6 - 5/c)/(2c^2), -(2 - 2/c)/c^2, -(1 + 2x)/(2c^3); both sums vanish.

### 5.4 Partial fractions
**Proposition 5.3 (L1 in the variable z).** R_0 = c z/(c + z) = c - c^2/(c + z). For k >= 1 there are a_j^(k) in Kx, 1 <= j <= k + 1, with
    R_k(z) = [t^k](z G) = Sum_{j=1}^{k+1} a_j^(k)/(c + j z),     a_j^(k) = j Res_{z=-c/j} R_k,     Sum_{j=1}^{k+1} a_j^(k) = 0.          (5.3)
Proof. For k = 0, (tau/t) F Ephi = H_0 = 1 at t = 0, so G_0 = c/(c + z). For k >= 1, Lemma 5.2 and s.5.3 show that the poles of R_k lie in
{-c/j : j <= k + 1} and are simple, and R_k(oo) = 0, so (PF) gives (5.3) without constant term; evaluating at z = 0 gives
0 = R_k(0) = (1/c) Sum_j a_j^(k). QED
With a_1^(0) := -c^2 and a_j^(k) := 0 for j > k + 1, the sum rule reads Sum_j a_j^(k) = -c^2 delta_{k0}.

### 5.5 Translation to psi: Theorem L1
**Theorem 5.4 (L1, all k).** For every k >= 0,
    q^(0) = [t^0](rho q_t) = q_0 = c - c^2 lam_1,      q^(k) = [t^k](rho q_t) = Sum_{j=1}^{k+1} a_j^(k) lam_j  in Lam_{k+1}   (k >= 1),
with the a_j^(k) of (5.3). Moreover g_0 = c lam_1, and for k >= 1
    g_k = Sum_{j=1}^{k+1} g_j^(k) lam_j = -(1/c) Sum_{j=1}^{k+1} j a_j^(k) lam_j  in Lam_{k+1},     g_k = O(psi^-2) at psi = oo.       (5.4)
Proof. By (C10) q_t(psi) = -psi r_{0,1}(-psi; t). Multiply by rho, which does not involve psi: rho q_t(psi) = -psi g(psi; t), and
[t^k](rho q_t)(psi) = R_k(-psi) as formal series (Taylor expansion commutes with z -> -psi). As 1/(c + j z) = lam_j(psi) at z = -psi,
Proposition 5.3 gives the first two statements. For (5.4): G_k = R_k/z and 1/(z (c + j z)) = (1/c)(1/z - j/(c + j z)), so by the sum rule
G_k(z) = (1/(c z)) Sum_j a_j^(k) - (1/c) Sum_j j a_j^(k)/(c + j z) = -(1/c) Sum_j j a_j^(k)/(c + j z); substitute z = -psi. Since R_k = O(1/z) at oo,
G_k = O(z^-2), i.e. g_k = O(psi^-2); equivalently, by (LAM4), Sum_j g_j^(k)/j = -(1/c) Sum_j a_j^(k) = 0. QED
Remarks. (i) g_j^(k) = -(j/c) a_j^(k) for all k >= 0 (at k = 0, g_1^(0) = c). (ii) Conversely, by (LAM2), L1 at order k >= 1 is equivalent to
"g_k in Lam_{k+1} and g_k = O(psi^-2)": -psi g_k = Sum_j g_j^(k)/j + Sum_j a_j^(k) lam_j, and the constant must vanish. This is the form used in
section 6. (iii) The three ingredients are separate: the pole structure (Lemmas 5.1, 5.2), the regularity theorem (s.5.3), and the
coefficientwise limit at z = oo (Lemma 5.2). The pole structure proved here is the upper-bound half of Fact G of conifold_regularity_note.md
s.6: poles only at z = -c/j with j <= k + 1, all simple. That every such pole occurs is proved only for the top pole and generic x
(Corollary 11.3), is checked for k <= 12 (V4), and fails at special x (at x = 1, [t^1](z G) = ((1 - x)^2/2)(1/(c + 2z) - 1/(c + z)) = 0).
Theorem R needs only the upper bound. (iv) F3 for Dil_1: Dil_1(d) = (-1)^d [z^(d-1)] r_{0,1} =
(-1)^d [z^(d-1)](rho^-1 G) is in Kx[[t]] (Proposition 4.4); its value at t = 0 is -c^(1-d), i.e. q_t|_{t=0} = q_0 = -c psi/(c - psi).
[Mathematica check (7) of hvc_leaf_proof_check.wls: [t^k](z G) for k <= 12, exact in x: regular at z = 0, limit c delta_{k0}, denominator
exactly prod_{j=1}^{k+1} (c + j z) (V4); S3_checks.pl and S4_checks.pl: the closed form against the dilaton tables t^0..t^5 of
leaf_structure.md and partial fractions at rational z for k <= 11.]


## 6. Step 7: the H-row in leaf form

### 6.1 The H-row
By (4.6) and (4.7),
    rho r_{1,1}(z) = (L_1 + z/t) G + z dG/dt,    sigma(psi) = rho v_1(-psi) = (sigma_0, sigma_1),   sigma_0 = g,
    sigma_1 = L_1 g - (psi/t) Dg = (L_1 - psi/t) g - psi dg/dt                                                                    (6.1)
(z/t becomes -psi/t and z dG/dt becomes -psi dg/dt under z = -psi). Both components of t sigma lie in Kx[[t]][[psi]], because t L_1 is in
Kx[[t]]. Only F1 enters (6.1): with rho defined by (C13), t^2 rho^2 Delta_1 = -x c^3 holds by definition, so the logarithmic derivative of
Delta_1 is that of (t rho)^-2; no square root of Delta_1 is needed. By Lemma 4.2(d),
    L_1 = c/t + l(t),   l = c (1/tau - 1/t) - x s* in Kx[[t]],   l_0 = c/2 - x/c = (1 + x^2)/(2c),   l_1 = c/12 - x^2/c^3.          (6.2)

### 6.2 The leaf structure of sigma
**Proposition 6.1.** (a) sigma^(-1) = [t^-1] sigma = (0, c).
(b) For k >= 0, [t^k] sigma_0 = g_k and
    [t^k] sigma_1 = (c - (k + 2) psi) g_{k+1} + Sum_{m=0}^{k} l_m g_{k-m}   in Lam_{k+1}.                                                (6.3)
Hence sigma = t^-1 (0, c) + Sum_{k>=0} t^k sigma^(k) with both components of sigma^(k) in Lam_{k+1}; explicitly
    sigma^(k) = Sum_{j=1}^{k+1} S_j^(k) lam_j,     S_j^(k) := ( g_j^(k),  (k + 2 - j) a_j^(k+1) + Sum_{m=0}^{k} l_m g_j^(k-m) )           (6.4)
(g_j^(n) := 0 for j > n + 1). In particular
    sigma^(0) = (c lam_1, x lam_1),
    sigma^(1) = ( ((1 - x)^2/(2c)) (lam_1 - 2 lam_2),  (x/(2c^2)) [ (1 - x + x^2) lam_1 - 2 (1 - x)^2 lam_2 ] ).
Proof. By (6.1) and (6.2), [t^k](L_1 g) = c g_{k+1} + Sum_{m=0}^{k} l_m g_{k-m} and [t^k]((psi/t) Dg) = psi [t^(k+1)] Dg = (k + 2) psi g_{k+1}; this
gives (6.3) for every k >= -1 (empty sum for k = -1). For k = -1, [t^-1] sigma_1 = (c - psi) g_0 = (c - psi) c lam_1 = c by (LAM3), and
[t^-1] sigma_0 = 0: this is (a). Leaf structure: write g_{k+1} = Sum_{j=1}^{k+2} g_j^(k+1) lam_j. By Theorem 5.4 (k + 1 >= 1) and (LAM4),
Sum_j g_j^(k+1)/j = 0. By (LAM3) and g_j^(k+1) = -(j/c) a_j^(k+1),
    (c - (k + 2) psi) g_{k+1} = (k + 2) Sum_j g_j^(k+1)/j + c Sum_j (1 - (k + 2)/j) g_j^(k+1) lam_j = Sum_{j=1}^{k+1} (k + 2 - j) a_j^(k+1) lam_j:
the factor c - (k + 2) psi turns the top leaf lam_{k+2} into a constant, and all the constants add up to (k + 2) Sum_j g_j^(k+1)/j = 0. The terms
l_m g_{k-m} lie in Lam_{k-m+1}, contained in Lam_{k+1}. Together with [t^k] sigma_0 = g_k this is (6.4). (In z: (6.3) at psi = -z is
(c + (k + 2) z) G_{k+1} + Sum_m l_m G_{k-m}; the factor c + (k + 2) z kills the simple pole of G_{k+1} at -c/(k + 2); the result is regular at 0,
O(1/z) at oo, with simple poles in {-c/j : j <= k + 1}, so (PF) puts it in span{1/(c + j z) : j <= k + 1}.)
k = 0: g_1 = ((1 - x)^2/(2c))(lam_1 - 2 lam_2) by Theorem 5.4 with (a_1^(1), a_2^(1)) = ((1 - x)^2/2)(-1, 1); then (c - 2 psi) g_1 = a_1^(1) lam_1 =
-((1 - x)^2/2) lam_1 and l_0 g_0 = (c^2/2 - x) lam_1, so [t^0] sigma_1 = ((c^2 - (1 - x)^2)/2 - x) lam_1 = x lam_1.
k = 1: (c - 3 psi) g_2 = 2 a_1^(2) lam_1 + a_2^(2) lam_2 with a_1^(2) = -(2 - 4x + 3x^2 - 4x^3 + 2x^4)/(12 c^2) and a_2^(2) = (1 - x)^4/(2c^2)
(Corollary 11.3, s.11.4); l_0 g_1 = ((1 + x^2)(1 - x)^2/(4c^2))(lam_1 - 2 lam_2) and l_1 g_0 = ((c^4 - 12 x^2)/(12 c^2)) lam_1. The lam_1-coefficient
is [-2(2 - 4x + 3x^2 - 4x^3 + 2x^4) + 3(1 + x^2)(1 - x)^2 + (1 + x)^4 - 12x^2]/(12 c^2) = 6x (1 - x + x^2)/(12 c^2), and the lam_2-coefficient is
(1 - x)^4/(2c^2) - (1 + x^2)(1 - x)^2/(2c^2) = -x (1 - x)^2/c^2. QED
Remarks. (i) sigma^(0)_1 = x lam_1 is the z-dependent part x/(c + z) of the H-row at t^0, so rho r_{1,1}(z) = c/t + x/(c + z) + O(t), the second half
of F2 (conifold_regularity_note.md s.1). (ii) sigma_1 is O(1/psi) but not O(psi^-2) (sigma^(0)_1 = x lam_1); this does not matter for the
leaf structure. (iii) Part (a) is the same statement as Proposition 4.8(b), proved there from (H) and (4.3); here it comes out of the leaf
expansion, using only F2 (Corollary 3.18), (6.2) and (6.1). [S3_checks.pl: sigma^(-1), sigma^(0), sigma^(1) and the symplectic identity
<t sigma(psi), t sigma(-psi)> = -x c^3 through t^5, sigma built from (6.3); S4_checks.pl ("sigma"); ref_s3_indep.pl.]

### 6.3 The projection of sigma on the regular eigen-covector
**Proposition 6.2.** pi(sigma(psi)) = -(x^2 - c x L_2) (psi/t) Dg(psi), and x^2 - c x L_2 = -x^3 tau s*^2. Thus
    pi(sigma(psi)) = x^3 (tau/t) s*^2 psi Dg(psi)   in Kx[[t]][[psi]].
Proof. sigma = g chi_1 + (sigma_1 - L_1 g)(0, 1); pi(chi_1) = <chi_1, chi_2> = 0, pi((0, 1)) = x^2 - c x L_2, and sigma_1 - L_1 g = -(psi/t) Dg by (6.1).
x^2 - c x L_2 = x (x - c L_2) = -x^3 tau s*^2 by Proposition 4.3(d). QED
(Equivalently, with pi(a) = x^3 s* (a_0 - tau s* a_1) and tau s* L_1 = 1: pi(sigma) = x^3 s* [g - tau s* (L_1 - psi/t) g + tau s* psi g'] =
x^3 tau s*^2 (psi/t)(t g)'.)

## 7. Step 8: the E_11 flow; F3 for E_11; L3; L2 for j = 1

### 7.1 E_11 in terms of sigma
By Lemma 1.1 and Proposition 4.8(d), in the leaf variables (psi, psi') = (z, w),
    rho^2 E_11(psi, psi') = [ rho^2 Delta_1 - <sigma(psi), sigma(psi')> ]/(psi + psi'),      rho^2 Delta_1 = -x c^3/t^2.
Put Y(psi) := t sigma(psi), an element of Kx[[t]][[psi]]^2 by Proposition 6.1, and
    Ebar(psi, psi') := t^2 rho^2 E_11(psi, psi') = Nbar(psi, psi')/(psi + psi'),     Nbar := -x c^3 - <Y(psi), Y(psi')>.                (7.1)
Nbar is in Kx[[t]][[psi, psi']], and Nbar(psi, -psi) = -x c^3 - t^2 rho^2 <v_1(-psi), v_1(psi)> = -x c^3 - t^2 rho^2 Delta_1 = 0 (Proposition 4.8(c),
(4.2)). So Nbar is divisible by psi + psi' (in the variables u = psi, v = psi + psi' it vanishes at v = 0), Ebar is in Kx[[t]][[psi, psi']], and it
is t^2 rho^2 times the code's difference quotient (quotients by the non-zero-divisor psi + psi' are unique). In particular rho^2 E_11 lies in
t^-2 Kx[[t]][[psi, psi']]. Since < , > is symmetric, Ebar(psi, psi') = Ebar(psi', psi), i.e. E_11(a, b) = E_11(b, a) (observed in leaf_structure.md,
Step 3). Without K_1 = 1 the subtracted term would be -K_1(w) Delta_1 and the symmetry would not be manifest.

### 7.2 The QDE for Y
**Lemma 7.1.** -psi dY(psi)/dt = (A - L_1) Y(psi) in Kx((t))[[psi]]^2; equivalently d sigma/dt = -sigma/t - (1/psi)(A - L_1) sigma.
Proof. This is Proposition 4.8(d) (the QDE (4.5) at z = -psi, with vn(-psi) = Y(psi); A and L_1 do not involve z). Sign: sigma is built from
v_1(-psi), so z = -psi and the 1/z of the QDE becomes -1/psi. QED

### 7.3 The flow
**Theorem 7.2 (E_11 flow; V8).** In Kx[[t]][[psi, psi']],
    d/dt Ebar(psi, psi') = Theta(t) Dg(psi) Dg(psi'),     Dg = d/dt (t g),                                                          (7.2)
    Theta = (L_2 - L_1)(x^2 - c x L_2)^2/Delta_2 = x^2 - c x L_2 = Delta_2/(L_2 - L_1) = -x^3 tau s*^2 = x c dphi' - x = -4 x^3 tau/(c + delta)^2.   (7.3)
Proof. (i) Since t^2 rho^2 Delta_1 is constant, only the pairing term of Nbar depends on t: dNbar/dt = -<dY(psi)/dt, Y(psi')> - <Y(psi), dY(psi')/dt>.
Multiply by psi psi' and use Lemma 7.1 and self-adjointness (Proposition 4.7(c)):
    psi psi' dNbar/dt = psi' <(A - L_1) Y(psi), Y(psi')> + psi <Y(psi), (A - L_1) Y(psi')> = (psi + psi') <(A - L_1) Y(psi), Y(psi')>.
(ii) psi + psi' is t-independent, so (psi + psi') dEbar/dt = dNbar/dt, and after cancelling psi + psi' (a non-zero-divisor),
psi psi' dEbar/dt = <(A - L_1) Y(psi), Y(psi')>. This is why the factor psi + psi' disappears: the flows of the two arguments contribute 1/psi
and 1/psi', self-adjointness makes the two terms equal, and 1/psi + 1/psi' = (psi + psi')/(psi psi').
(iii) By (4.8) and Proposition 6.2, pi(Y(psi)) = t pi(sigma(psi)) = -(x^2 - c x L_2) psi Dg(psi), so
psi psi' dEbar/dt = ((L_2 - L_1)/Delta_2)(x^2 - c x L_2)^2 psi psi' Dg(psi) Dg(psi'). Cancel psi psi': this is (7.2) with the V8 definition of Theta.
(iv) By (4.1), Delta_2 = (L_2 - L_1)(x^2 - c x L_2), hence Theta = x^2 - c x L_2 = Delta_2/(L_2 - L_1); x^2 - c x L_2 = -x^3 tau s*^2 (Proposition 4.3(d));
dphi' = 1 - x s* = 1 - L_2, so x c dphi' - x = x c - x - c x L_2 = x^2 - c x L_2; and s*^2 = 4/(c + delta)^2. QED
Expansion (tau = t - t^2/2 + ..., s* = 1/c + x tau/c^3 + ...):
    Theta = -x^3 t/c^2 + x^3 (1 - x)^2 t^2/(2 c^4) - x^3 (1 - 8x + 12x^2 - 8x^3 + x^4) t^3/(6 c^6) + O(t^4),   Theta_0 = 0,
and Theta_k = (k + 1) x c [t^(k+1)] dphi for k >= 1, in agreement with the printed Theta(t) and dphi(t) of hvc_leaf_proof_check_out.txt
through t^7. Remarks. (a) Theta = Delta_2/(L_2 - L_1): the loop flow is governed by the regular idempotent; Delta_2 is regular and non-zero
while L_2 - L_1 has a simple pole, so Theta = O(t). (b) Rank one: before expanding in t, the right side of (7.2) is a product of a function of
(t, psi) and a function of (t, psi'), because A - L_1 has rank one (Proposition 4.7(e)). Its t^k-coefficients are not of rank one in general;
they are sums of rank-one terms Theta_{k0}(k1 + 1)(k2 + 1) g_{k1}(psi) g_{k2}(psi') (Theorem 7.3); e.g. the t^1-coefficient of rho^2 E_11 has
e^(1)_11, e^(1)_12 = e^(1)_21 non-zero and e^(1)_22 = 0, a matrix of rank 2. What rank one gives is that every loop coefficient is built from
products of the same leaf series g_k. [Mathematica (V8), check (6) of hvc_leaf_proof_check.wls: every (m, n) with m + n <= 8, t^0..t^6;
referee_scripts/ref_s3_indep.pl: the code's edge formula against the integrated flow for m + n <= 8, t^0..t^6 at four rational x.]

### 7.4 Integration from t = 0; F3 for E_11
By Proposition 6.1(a), [t^0] Y = (0, c) and <(0, c), (0, c)> = -c x c^2 = -x c^3. Hence [t^0] Nbar = -x c^3 + x c^3 = 0 and [t^0] Ebar = 0. This
is where K_1 = 1 enters: for a t-independent K_1 one would have P_11(-w, w) = -K_1(w) Delta_1 and [t^0] Ebar = -x c^3 (K_1(psi') - 1)/(psi + psi').
(The constant also vanishes for a second reason: the numerator of the difference quotient vanishes at psi = -psi' by construction, while its
t^0-coefficient after multiplication by t^2 rho^2, -x c^3 (K_1(psi') - 1), does not depend on psi; this is the proof of Proposition 4.8(c) seen
through (7.1).) The kernel of d/dt on Kx[[t]] is Kx and [t^(n+1)] Ebar = [t^n](dEbar/dt)/(n + 1), so (7.2) and [t^0] Ebar = 0 give
    Ebar = Int_0^t Theta(u) Dg(psi; u) Dg(psi'; u) du  (termwise),   i.e.   [t^(n+1)] Ebar = [t^n](Theta Dg(psi) Dg(psi'))/(n + 1).        (7.4)
As Theta_0 = 0, also [t^1] Ebar = 0. So Ebar is in t^2 Kx[[t]][[psi, psi']], and E_11 = (t rho)^-2 Ebar is in Kx[[t]][[psi, psi']] since rho^-2 is in
Kx[[t]] (Proposition 4.4): E_11 is regular at t = 0. This is F3 for E_11. (In the point frame and in frame A the same factor has a double pole,
conifold_regularity_note.md F3; in frame B the double pole of Delta_1 is exactly compensated.)

### 7.5 The leaf structure of the loops (L3)
**Theorem 7.3 (L3, all k).** For every k >= 0,
    [t^k](rho^2 E_11)(psi, psi') = (1/(k + 2)) Sum_{k0 >= 1, k1, k2 >= 0, k0 + k1 + k2 = k + 1} Theta_{k0} (k1 + 1)(k2 + 1) g_{k1}(psi) g_{k2}(psi')     (7.5)
                               = Sum_{i, l >= 1, (i-1) + (l-1) <= k} e^(k)_il lam_i(psi) lam_l(psi'),
    e^(k)_il = (1/(k + 2)) Sum_{k0 >= 1, k1 + k2 = k + 1 - k0} Theta_{k0} (k1 + 1)(k2 + 1) g_i^(k1) g_l^(k2),     e^(k)_il = e^(k)_li,
with g_i^(m) := 0 for i > m + 1. There is no lam_0 in either variable. On the anti-diagonal (i - 1) + (l - 1) = k, i.e. i = k1 + 1, l = k2 + 1 with
k1 + k2 = k,
    e^(k)_{k1+1, k2+1} = (-1)^(k+1) x^3 (k1 + 1)(k2 + 1) P_{k1} P_{k2}/((k + 2) (k1!)^2 (k2!)^2 c^(2k)),                                 (7.6)
which is non-zero whenever x P_{k1}(x) P_{k2}(x) != 0. At t = 0 the loop has rank one: [t^0](rho^2 E_11) = (1/2) Theta_1 g_0(psi) g_0(psi') =
-(x^3/2) lam_1(psi) lam_1(psi'), hence E_11(z, w)|_{t=0} = -(x^3/2) lam_1(z) lam_1(w) (rho(0) = 1), the E_11 part of F3 of
conifold_regularity_note.md.
Proof. [t^k](rho^2 E_11) = [t^(k+2)] Ebar = [t^(k+1)](Theta Dg(psi) Dg(psi'))/(k + 2) by (7.4); insert Theta = Sum_{k0>=1} Theta_{k0} t^k0 and
Dg = Sum_k (k + 1) g_k t^k. Support: every term has k1 + k2 = k + 1 - k0 <= k; by Theorem 5.4, g_{k1} is in Lam_{k1+1} and g_{k2} in Lam_{k2+1},
so the term lies in span{lam_i lam_l' : i <= k1 + 1, l <= k2 + 1}, of total degree (i - 1) + (l - 1) <= k1 + k2 <= k. Symmetry: (7.5) is
invariant under (psi, k1) <-> (psi', k2). Anti-diagonal: a product lam_{k1+1}(psi) lam_{k2+1}(psi') with k1 + k2 = k comes only from terms with
k1' >= k1, k2' >= k2 and k1' + k2' <= k, i.e. k1' = k1, k2' = k2 and k0 = 1; its coefficient is (1/(k + 2)) Theta_1 (k1 + 1)(k2 + 1)
g^(k1)_{k1+1} g^(k2)_{k2+1}. Insert Theta_1 = -x^3/c^2 and g^(k)_{k+1} = (-1)^k P_k/((k!)^2 c^(2k-1)) (Corollary 11.3). QED
Generating form: the lam_i-coefficient of Dg is -(i/c) Ahat_i' (Theorem 11.1), so e_il(t) = (i l/(c^2 t^2)) Int_0^t Theta(u) Ahat_i'(u) Ahat_l'(u) du
(Theorem 11.4). The case k1 = 0 of (7.6) is the corner e^(k)_{1,k+1} = ((k + 1)^2/(k + 2))(x^3/c^2) a_{k+1}^(k) of leaf_structure.md, Step 3.
[S3_checks.pl: (7.5) against the restricted loop tables t^0..t^5 (every entry, support, symmetry, all corners) and the full loop matrices
t^0..t^3; S4_checks.pl: every loop entry k <= 5, corners k <= 6, anti-diagonal i, l <= 5.]

### 7.6 Powers of rho; L2 for j = 1
In Theorem R every half-edge at a singular vertex of genus h >= 2 carries one factor rho (section 13). So a loop, or an edge between two such
vertices, carries rho^2 E_11; a leg, i.e. an E_11-edge to a vertex that is not rescaled (singular of genus 0 or 1), with external power b at the
far end, carries rho f_{1,b}; an E_11-edge between two non-rescaled vertices carries E_11.
**Corollary 7.4.** For every integer e and every k >= 0, [t^k](rho^e E_11) lies in span{lam_i(psi) lam_l(psi') : i, l >= 1, (i-1) + (l-1) <= k}. In
particular (L2 for j = 1), for every b >= 0,
    [t^k](rho f_{1,b}) = Sum_{i=1}^{k+1} ( Sum_{m=0}^{k} [t^m](rho^-1) Sum_{l} e^(k-m)_il l^b/c^(b+1) ) lam_i(psi)   in Lam_{k+1},       (7.7)
and [t^k] f_{1,b} lies in Lam_{k+1} as well.
Proof. rho^(e-2) is in Kx[[t]] (Proposition 4.4) with psi-independent coefficients, so [t^k](rho^e E_11) = Sum_m [t^m](rho^(e-2)) [t^(k-m)](rho^2 E_11),
and the supports at the orders k - m <= k are contained in the support at order k. For the leg, f_{1,b}(psi) = [psi'^b] E_11(psi, psi') and
rho f_{1,b} = [psi'^b](rho E_11) (e = 1); by (LAM1), [psi'^b] lam_l(psi') = l^b/c^(b+1), and the lam_i(psi) that occur have i - 1 <= k. QED
Explicitly, with rho^-1 = 1 + ((1 - x)^2/(2c^2)) t + ((1 + 4x - 12x^2 + 4x^3 + x^4)/(12 c^4)) t^2 + O(t^3),
    [t^0](rho f_{1,b}) = -(x^3/2) c^(-b-1) lam_1,      [t^1](rho f_{1,b}) = ((1 - x)^2 x^3/(12 c^(b+3))) [ (2^(b+3) - 9) lam_1 + 8 lam_2 ].
[t^1](rho f_{1,b}) = [t^1](rho^2 f_{1,b}) + [t^1](rho^-1) [t^0](rho^2 f_{1,b}) with [t^1](rho^2 f_{1,b}) = c^(-b-1)[(e^(1)_11 + 2^b e^(1)_12) lam_1 + e^(1)_21 lam_2]
and [t^0](rho^2 f_{1,b}) = -(x^3/2) c^(-b-1) lam_1: the lam_1-coefficient is ((1 - x)^2 x^3/c^(b+3))(-1/2 + (2/3) 2^b - 1/4), the lam_2-coefficient
(2/3)(1 - x)^2 x^3/c^(b+3). At t^0 the b-dependence is only through c^-b: the rank-one observation of conifold_regularity_note.md s.6,
step 3(i). [S3_checks.pl: every tabulated j = 1 leg, b = 0..7, t^0..t^5 (33 entries per x); S4_checks.pl ("E11 legs table", 524 comparisons).]

## 8. The regular idempotent (Step 10, and the second QDE row at U_2)

### 8.1 The conjugate encoding
**Lemma 8.1.** (i) cconj is a ring automorphism of QA (the Galois conjugation U -> Ssum - U); it commutes with cadd, cmul, cinv, with cscale by
elements of Q(x, Q), with cDerN and with twistRows. (ii) Consequently, at alpha0 = 0: ratCN[a,2][m] = cconj[ratCN[a,1][m]] (a = 0, 1, all m,
before and after twistRows), LCN[2] = cconj[LCN[1]], DeltaCN[2] = cconj[DeltaCN[1]], dlogDeltaN[2] = cconj[dlogDeltaN[1]], and evalQ gives exactly
    EdN[2,2][m, n] = cconj[EdN[1,1][m, n]],   EdN[2,1][m, n] = cconj[EdN[1,2][m, n]],   DilN[2][n] = cconj[DilN[1][n]].
(iii) The true values at idempotent 2 are ev_1 of the stored elements: r^pt_{0,2}(z; t) = Sum_m z^m wickZ_m(W_2(U_2), W_3(U_2), ...), and
r_{1,2} = (L_2 + z D_2) r_{0,2} with d/dt the total derivative along U_2(t) and Delta_2 = ev_2(DeltaCN[1]).
Proof. (i) ev_1, ev_2 are ring homomorphisms, jointly injective, and ev_1 o cconj = ev_2, ev_2 o cconj = ev_1 (U_1 + U_2 = Ssum); this gives the
compatibility with the ring operations. cDerN[C[a, b]] = C[d_Q a, d_Q b] + b dUC with dUC = toC[-D[LeqU, Q]/D[LeqU, U]], and
ev_i(dUC) = dU_i/dQ (implicit differentiation), so ev_i o cDerN = (d/dQ) o ev_i. For {i, i'} = {1, 2}: ev_i(cDerN cconj C) = d/dQ ev_{i'}(C) =
ev_{i'}(cDerN C) = ev_i(cconj cDerN C). twistRows multiplies the rows of both idempotents by the same scalar series. (ii) rmatrixAt0 uses
ffc1 = ffc2 = (1, 0, 0, ...), so ratCN[0,2][m] = cconj[wickCN[m]] = cconj[ratCN[0,1][m]]; applying (i) to the H-row formula of rmatrixCore,
whose only scalar is Q, gives ratCN[1,2][m] = cconj[ratCN[1,1][m]]. etainvCur contains no U, so the edge and dilaton formulas of evalQ commute
with cconj. (iii) ev_1 o cconj = ev_2, and ev_2 commutes with total d/dQ. QED

### 8.2 Analyticity at U_2
**Lemma 8.2.** (i) For k >= 1, W_k is a polynomial in (1 - U)^-1 and (1 - Q U)^-1 with coefficients in Z[x] (W_1 = x - x (1 - Q U)^-1 -
(1 - U)^-1). U_2(t) is analytic at t = 0 with U_2(0) = -1/x, and 1 - U_2(0) = 1 - Q(0) U_2(0) = c/x != 0; so every W_k(U_2(t)) is analytic at t = 0,
and W_2(U_2) = x delta s*^2 is a unit of O with W_2(U_2(0)) = x/c.
 (ii) (R) Every r^pt_{0,2}[m], r_{0,2}[m] and r_{1,2}[m] is analytic at t = 0; r_{0,2}[0] = 1 and r_{1,2}[0] = L_2.
 (iii) Delta_2 = x^3 delta s*^2 is analytic at t = 0 with Delta_2(0) = x^3/c and Delta_2'(0) = 0 (F1), Delta_2 = x^2 W_2(U_2) = x^2 U_2^2 phi_tau''(s~)/tau^2,
     and Delta_1 Delta_2 = -x^4 delta^2/tau^2.
 (iv) Dil_2(n) = (-1)^n r_{0,2}[n-1] and every E_22(m, n) are analytic at t = 0 (F3 for Dil_2, E_22).
Proof. (i) theta log U = 1, theta log(1 - U) = 1 - (1 - U)^-1, theta log(1 - Q U) = 1 - (1 - Q U)^-1, and theta (1 - a U)^-n = n ((1 - a U)^-(n+1) -
(1 - a U)^-n) for a in {1, Q}. U_2 is a simple root of LeqU at Q = 1 ((x U + 1)(U - 1) = 0 has the simple roots 1 and -1/x, distinct since c != 0),
hence analytic in t. W_2(U_2) = x delta s*^2 by (3.5); at t = 0, W_2 = -c U/(1 - U)^2 at U = -1/x equals x/c. (ii) By Lemma 8.1(iii) r^pt_{0,2}[m] is a
polynomial in W_2(U_2)^-1, W_3(U_2), ..., W_{2m+2}(U_2), hence analytic; the twist h_w^-1 is t-independent; and r_{1,2}[m] = L_2 r_{0,2}[m] +
D_2 r_{0,2}[m-1] with L_2 = x s* in O and Delta_2 a unit of O. (iii) Proposition 4.3(c)-(e) and (3.5). (iv) By Lemma 1.1, E_22(m, n) is a polynomial
with constant coefficients in the entries of v_2[p], p <= m + n + 1. QED
[S4_checks_r1.pl B1 (Delta_2 four ways, Delta_2'(0) = 0) at four x; S2_checks_reg.pl part T.]

### 8.3 The data at the regular critical point
Why a separate argument is needed. Row 1 at i = 2 is the definition (H). Row 2 at i = 2 cannot be obtained as in Proposition 4.6, because the
closed form concerns the singular row only; nor does the route of Proposition 3.12 apply, which multiplies the formal expansion by
Exp[dphi/z] (f2/f2(0))^(-1/2) and expands at tau = 0 in Kz[[tau]]: at s~ -> oo, e^{phi_tau(s~)/z} is tau^(-c/z) times a unit, d/dtau phi_tau(s~)
has a simple pole and phi_tau''(s~) a double zero. Instead we work coefficientwise in z at fixed small t != 0, with no expansion at tau = 0;
analyticity at tau = 0 enters only through units of O and values at t = 0.
**Lemma 8.3.** For tau in D'' (small t != 0):
 (a) s~ = (c + delta)/(2 x tau) = 1/(x tau s*) is a nondegenerate critical point of phi_tau, s~ is not in {0, 1, 1/tau}; tau s~ is in O with value
     c/x; 2 x tau s~ - c = delta; x s~ = L_1; U_2 = 1 - tau s~ is a unit of O, U_2(0) = -1/x.
 (b) f2 := phi_tau''(s~) = delta/(x s~^2 U_2^2) = tau^2 Delta_2/(x^2 U_2^2) = tau^2 f2r, with f2r := Delta_2/(x^2 U_2^2) a unit of O, f2r(0) = x^3/c.
 (c) f0(tau) := phi_tau(s~(tau)) (any branch) has f0_tau = s~/(1 - tau s~) in tau^-1 O, with leading term -c/tau, and
        f0' = (1 - tau) s~/(1 - tau s~) = 1 - x s~ = 1 - L_1 = L_2 - mu.                                                          (8.1)
     u_2 := W(U_2) satisfies du_2/dt = L_2 and u_2 = x t + c log tau + f0 + const.
Proof. (a) Lemma 3.1(a), s.1.3, Lemma 4.2(d); U_2(0) = (x - 1 - c)/(2x) = -1/x, and U_2 has no zero on D'' (Lemma 3.1(a)). (b) (3.3) at s~ and (3.5)
with Delta_2 = x^2 W_2(U_2) (Proposition 4.3(c)); f2r is a unit since Delta_2 and U_2 are, with f2r(0) = (x^3/c)/(x^2 x^-2). (c) Envelope formula and
Lemma 3.1(b)(iii) at s~: (1 - x s~)(1 - tau s~) = 1 - (x + tau) s~ + x tau s~^2 = (c - x - tau) s~ = (1 - tau) s~. tau s~ is in O with value c/x and
U_2 is a unit with value -1/x, so f0_tau = tau^-1 (tau s~)/U_2 has leading term -c/tau. The last equality in (8.1) is Lemma 4.2(e). QED

### 8.4 The second QDE row at the regular idempotent
Fix tau0 in D'' and an open disc Db around tau0 whose closure lies in D''. Fix on Db branches of log s, log(1 - s), log(1 - tau s) near
s = s~(tau) (possible since s~, 1 - s~ and 1 - tau s~ do not vanish on the closure of Db); then phi_tau(s) is analytic near every point
(tau, s~(tau)), tau in Db. O_Db denotes the analytic functions on Db and O_Db((z)) the Laurent series in z, bounded below, with coefficients
in O_Db (d/dt = (1 - tau) d/dtau coefficientwise). Put, with f0 and f2 of Lemma 8.3,
    Dt := d/dt + f0'/z - f2'/(2 f2),      Dop := Dt/(1 - tau) = d/dtau + f0_tau/z - f2_tau/(2 f2).

**Theorem 8.4 (Row 2 of the QDE at U_2).** In M[[z]] (in O[[z]] after multiplication by tau),
    (L_2 + z D_2) r_{1,2} = -(x/tau) r_{0,2} + (c/tau) r_{1,2},   equivalently   z D_2 r_{1,2} = L_1 (r_{1,2} - L_2 r_{0,2}),                (8.2)
i.e. (D_2 r_{1,2})[m-1] = L_1 (r_{1,2}[m] - L_2 r_{0,2}[m]) for m >= 1 (the form tested in V5 for i = 2). With Row 1 (the definition (H)) this is
z D_2 v_2 = (A - L_2) v_2, and the same holds for the point-frame rows h_w v_2.
Proof. Step A (the regular row as a normalised Laplace expansion in the Euler coordinate). For every tau in Db,
    r^pt_{0,2}(z; t) = U_2 L2(tau),     L2(tau) := Lnorm[phi_tau, (1 - tau s)^-1; s~(tau)]   in O_Db[[z]].                                (8.3)
This is Proposition 3.10, steps (1)-(2), with (U_1, s*) replaced by (U_2, s~). Let y2 be a branch of log U_2 and Wy(y) := W(e^y); Wy is analytic
near y2 with Wy^(k)(y2) = W_k(U_2), Wy'(y2) = 0 and Wy''(y2) = W_2(U_2) = x delta s*^2 != 0; the integrand is e^{Wy/z} dy with amplitude 1, so by
(D1) and Lemma 8.1(iii), r^pt_{0,2} = Lnorm[Wy, 1; y2] = Lhat[Wy, 1; y2]. hc(s) := log(1 - tau s), with the branch for which hc(s~) = y2, is
analytic near s~ with hc'(s~) = -tau/U_2 != 0. By (3.1), Wy(hc(s)) = phi_tau(s) + C(t) near s~ with C(t) = c log tau + x t + x log(-1) independent of
s, and the amplitude 1 becomes -tau (1 - tau s)^-1. By (P1)(iii), (D3) and linearity, Lhat[Wy, 1; y2] = Lhat[phi_tau, (1 - tau s)^-1; s~] = U_2 L2,
since the amplitude takes the value 1/U_2 at s~. By (D2) the z-coefficients of L2 are polynomials in 1/f2, in phi_tau^(k)(s~) (3.6) and in
tau^a U_2^(-a-1), so they lie in O_Db.
Step B (conjugation). As operators on O_Db((z)),
    (L_2 + z D_2) o U_2 = U_2 o (z Dt + omega).                                                                                  (8.4)
Indeed (L_2 + z D_2)(U_2 Y) = U_2 [z Y' + L_2 Y + z (U_2'/U_2 - Delta_2'/(2 Delta_2)) Y]. By Lemma 8.3(b), Delta_2 = x^2 U_2^2 f2/tau^2, so
Delta_2'/(2 Delta_2) = U_2'/U_2 + f2'/(2 f2) - tau'/tau with tau'/tau = (1 - tau)/tau; and L_2 = f0' + mu by (8.1). So the bracket equals
z Y' + f0' Y - z (f2'/(2 f2)) Y + (mu + z (1 - tau)/tau) Y = (z Dt + omega) Y, since mu + z (1 - tau)/tau = omega. Only the logarithmic derivative of
Delta_2 enters, so no square root is involved. (8.4) is where V1 at U_2 (Delta_2 = x^2 W_2(U_2), Proposition 4.3(c)) is used: if Delta_2 were
x^2 W_2(U_2) times a non-constant p(t), the right side of (8.4) would contain the extra term -(z/2)(p'/p) U_2 Y, and Step D would not give Row 2.
Step C (Gauss's equation for L2). With alpha = beta = 1/z + 1, gamma = c/z + 2,
    Gauss_Dop[L2] := tau (1 - tau) Dop^2 L2 + (gamma - (alpha + beta + 1) tau) Dop L2 - alpha beta L2 = 0   in O_Db((z)).                   (8.5)
(P2) in the normalised form (3.9) applies to f(tau, s) = phi_tau(s), to the analytic nondegenerate critical point s~(tau) and to every amplitude
q(tau, s) that is a Laurent polynomial in z with coefficients analytic near (tau1, s~(tau1)), tau1 in Db: Dop Lnorm[phi_tau, q; s~] =
Lnorm[phi_tau, Dam q; s~], with the f0 and f2 of s~. Apply it to q = am0 = (1 - tau s)^-1 and to q = Dam am0, and use (D3):
    Gauss_Dop[L2] = Lnorm[ tau (1 - tau) Dam^2 am0 + (gamma - (alpha + beta + 1) tau) Dam am0 - alpha beta am0 ; s~ ]
                  = -(1 + 1/z) Lnorm[ kE' + phi_tau' kE/z ; s~ ] = 0,
by the amplitude identity (3.16), valid identically in (s, tau), and (P3). The amplitudes am0, Dam am0, Dam^2 am0 and kE are analytic near s~
because s~ is not in {0, 1, 1/tau}; that s~ -> oo as tau -> 0 plays no role, since Db does not contain 0.
Step D (transfer by V7). z Dt = Dz + a0 with the function a0 := f0' - z f2'/(2 f2) satisfies the Leibniz rule z Dt (a Y) = z a' Y + a z Dt Y, and
(z (1 - tau))^-1 z Dt = Dop. So Lemma 4.5 (with Dd' = z Dt) gives
    (z Dt + omega)^2 - (c/tau)(z Dt + omega) + x/tau = ((1 - tau) z^2/tau) Gauss_Dop   on O_Db((z)).                                      (8.6)
By (8.3), (8.4) applied twice, (8.6) and (8.5),
    [(L_2 + z D_2)^2 - (c/tau)(L_2 + z D_2) + x/tau] r^pt_{0,2} = U_2 [(z Dt + omega)^2 - (c/tau)(z Dt + omega) + x/tau] L2
                                                                = U_2 ((1 - tau) z^2/tau) Gauss_Dop[L2] = 0.
h_w does not depend on t, so the same identity holds for r_{0,2}; with r_{1,2} = (L_2 + z D_2) r_{0,2} it is (8.2) on Db.
Conclusion. Every z-coefficient of (8.2) is an identity of analytic functions on Db; tau0 was arbitrary, so it holds on D''. All its terms are
meromorphic at t = 0 (Lemma 8.2; x/tau and c/tau have simple poles), so it holds in M. The second form uses c/tau - L_2 = L_1 and x/tau = L_1 L_2.
The components of (A - L_2) v_2 are r_{1,2} - L_2 r_{0,2} = z D_2 r_{0,2} (Row 1) and -(x/tau) r_{0,2} + (c/tau - L_2) r_{1,2} = z D_2 r_{1,2} (Row 2). QED
Why the normalised forms suffice: the code's row is itself normalised, so (P1)(iii) identifies it with U_2 L2 exactly, with no convention; the
t-dependent factors that the Lfull form of (P1) carries (e^{C/z}, the Jacobian, and (-f2)^(-1/2) in place of Delta_2^(-1/2)) are accounted for in
(8.4), through mu, U_2, tau'/tau and Lemma 8.3(b). [Mathematica (V5) for i = 2, m = 1..9, t^0..t^6; S2_checks_reg.pl part P (Step A for z^0..z^6
and Step C for z^-2..z^4 as exact jets at four points tau0 != 0, with a negative control) and part T (Row 2 as t-series at t = 0, with a negative
control); S4_checks_r1.pl A2 (Step A, z^0..z^5, at three points t != 0) and B3 (Row 2 at U_2, z^1..z^6, t^0..t^6, four x).]

### 8.5 Orthogonality
**Corollary 8.5 (orthogonality).** <v_1(w), v_2(-w)> = 0 in M[[w]], to all orders in w.
Proof. Y12(z) := <v_1(z), v_2(-z)> lies in M[[z]]. Use z D_1 v_1 = (A - L_1) v_1 (Proposition 4.6), (8.2) with z replaced by -z
(-z D_2 v_2(-z) = (A - L_2) v_2(-z); A, L_2, D_2 do not involve z), v_i' = D_i v_i + (Delta_i'/(2 Delta_i)) v_i, the t-independence of eta^-1 and
self-adjointness (Proposition 4.7(c)):
    z Y12' = <(A - L_1) v_1(z), v_2(-z)> - <v_1(z), (A - L_2) v_2(-z)> + (z/2)(Delta_1'/Delta_1 + Delta_2'/Delta_2) Y12
           = [ (L_2 - L_1) + (z/2)(Delta_1'/Delta_1 + Delta_2'/Delta_2) ] Y12.
Write Y12 = Sum_m Y_m z^m with Y_m in M. The z^0 coefficient gives (L_2 - L_1) Y_0 = 0, and for m >= 1 the z^m coefficient gives
(L_2 - L_1) Y_m = Y_{m-1}' - (1/2)(Delta_1'/Delta_1 + Delta_2'/Delta_2) Y_{m-1}. Since L_2 - L_1 = -delta/tau is a unit of M, Y_0 = 0 and inductively
Y_m = 0. QED
It is essential that Y12 is a power series in z, bounded below: the formal symbol e^{-(u_1 - u_2)/z}, an infinite series in 1/z (and not even an
element of Rz, since u_1 - u_2 contains c log tau), would formally give a non-zero flat section (this is the outline's "e^{(u_1 - u_2)/z} P
t-independent" argument made rigorous). Y_0 = <chi_1, chi_2> = 0 is also direct. Row 2 at U_2 is genuinely needed: Row 1 only determines
r_{1,2} from r_{0,2}, and the z^m coefficient of Y12 contains r_{0,2}[m] only through <chi_1, chi_2> r_{0,2}[m] = 0, so Y12 = 0 is a genuine
constraint on the regular row, which the closed form of the singular row cannot supply. [Mathematica (V6) for i != j, w^0..w^9;
S2_checks_reg.pl part P; referee_scripts/ref_S2_jets.pl and s4ref/ref_indep.pl.]

### 8.6 The regular rows at t = 0
**Proposition 8.6.** (a) [analytic input] r^pt_{0,2}(z; 0) = h_w(z)^-1, equivalently r_{0,2}(z; 0) = h_w(z)^-2.
 (b) (Uses Theorem 8.4.) x r_{0,2}(z; 0) = c r_{1,2}(z; 0), i.e. v_2(z; 0) = r_{0,2}(z; 0) (1, x/c) = h_w(z)^-2 (1, x/c). Equivalently [t^1] r_{0,2}(z; t) = 0,
     i.e. d/dt [Delta_2^(-1/2) r_{0,2}] = 0 at t = 0.
 (c) (Uses (a), (b).) E_22(z, w)|_{t=0} = (x^3/c) [1 - h_w(z)^2 h_w(w)^2]/(z + w), and Dil_2(n)|_{t=0} = (-1)^n [z^(n-1)] h_w(z)^-2 (by (a) alone).
Proof. (a) By Lemma 8.2 each r^pt_{0,2}[m] is analytic at t = 0, and its value there is wickZ_m evaluated on the t = 0 data (theta^k W_0)(-1/x),
W_0(U) := c log(1 - U) - log U (spotU at Q = 1). By (D1), r^pt_{0,2}(z; 0) = Lnorm[Wy_0, 1; y_c] with Wy_0(y) := W_0(e^y) and y_c a branch of
log(-1/x). We identify this series with Lemma 3.8, first for real x > 0.
 (i) The phase. Write y_c = y0 + i pi with y0 := -log x, and put Fb(y) := y - c log(1 + e^y) for real y. Since e^{y + i pi} = -e^y,
     Wy_0(y + i pi) = -Fb(y) - i pi, so the Taylor coefficients of order >= 1 of Wy_0 at y_c are those of -Fb at y0. In (D1) the replacement
     (f, z) -> (-f, -z) leaves (-z/f2)^n and Rr_f/z unchanged, so Lnorm[-f, am; u0](z) = Lnorm[f, am; u0](-z). With (D3),
        r^pt_{0,2}(z; 0) = Lnorm[-Fb, 1; y0](z) = Lnorm[Fb, 1; y0](-z).                                                           (8.7)
     Equivalently: at z = -1/lambda and U = -e^y, W_0/z = lambda Fb(y) + const, so the Wick series in z is the Laplace series in 1/lambda = -z.
 (ii) Laplace. Fb is real-analytic on I = R, Fb'(y) = (1 - x e^y)/(1 + e^y) vanishes only at y0, Fb' > 0 on (-oo, y0) and < 0 on (y0, oo);
     Fb''(y0) = -c e^{y0}/(1 + e^{y0})^2 = -x/c < 0 and Fb(y0) = x log x - c log c. So sup_{|y - y0| >= e} Fb = max(Fb(y0 - e), Fb(y0 + e)) < Fb(y0)
     and Int_R e^{Fb(y)} dy = Int_0^oo (1 + u)^(-c) du = 1/x < oo (u = e^y). Lemma 3.8 (am = 1, z0 = 1) gives, for every M, as z' -> 0+,
        Int_R e^{Fb(y)/z'} dy = e^{Fb(y0)/z'} (2 pi z' c/x)^(1/2) [ Sum_{m<M} ([z'^m] Lnorm[Fb, 1; y0]) z'^m + O(z'^M) ],
     and with u = e^y, lambda = 1/z' the left side is Int_0^oo u^(lambda - 1)(1 + u)^(-c lambda) du = B(lambda, x lambda) (DLMF 5.12.3).
 (iii) Stirling (as in Proposition 3.15(3)) for log B(lambda, x lambda) = log Gamma(lambda) + log Gamma(x lambda) - log Gamma(c lambda): the terms
     -u cancel, the terms in log lambda add up to -(1/2) log lambda, and for every K
        B(lambda, x lambda) = (2 pi/lambda)^(1/2) x^(x lambda - 1/2) c^(1/2 - c lambda) exp( Sum_{k=1}^{K} bk_k St_k z'^(2k-1) + O(z'^(2K+1)) ),
     whose prefactor equals e^{Fb(y0)/z'} (2 pi z' c/x)^(1/2) = x^(x lambda) c^(-c lambda) (2 pi c/(lambda x))^(1/2). By uniqueness of asymptotic
     expansions, Lnorm[Fb, 1; y0](z') = exp(Sum_k bk_k St_k z'^(2k-1)) = h_w(z'), and by (8.7) r^pt_{0,2}(z; 0) = h_w(-z) = h_w(z)^-1.
 (iv) All x. Both sides have z-coefficients in Q(x): on the left by (D2), since Fb' = 1 - c sg with sg(y) := e^y/(1 + e^y), sg' = sg (1 - sg), so for
     k >= 2 Fb^(k) = -c sg^(k-1) is a polynomial in sg with coefficients in Z[c], with sg(y0) = 1/c and 1/Fb''(y0) = -c/x; on the right by definition.
     Frame B multiplies by a further h_w^-1.
 Check at z^1, by hand: at U = -1/x, Q = 1, W_2 = x/c, W_3 = -x (1 - x)/c^2, W_4 = x (1 - 4x + x^2)/c^3, and (D1) gives W_4/(8 W_2^2) - 5 W_3^2/(24 W_2^3)
 = -(1 + x + x^2)/(12 x c) = -(1/12) St_1 = [z^1] h_w^-1. At idempotent 1 the same method gives the opposite sign of the tail:
 Abeta = h_w c/(c + z) (Proposition 3.15).
(b) First proof, from Theorem 8.4 alone. In (D_2 r_{1,2})[m-1] = L_1 (r_{1,2}[m] - L_2 r_{0,2}[m]) the left side is analytic at t = 0 (Lemma 8.2),
while L_1 = c/t + O(1) has a simple pole; so r_{1,2}[m](0) = L_2(0) r_{0,2}[m](0) = (x/c) r_{0,2}[m](0) for m >= 1, and for m = 0 by (H). Second
proof, from Corollary 8.5: put w = -z and multiply by rho, <sigma(z), v_2(z)> = 0; its t^-1 coefficient is <(0, c), v_2(z; 0)> =
c x (x r_{0,2}(z; 0) - c r_{1,2}(z; 0)) (by <(0, c), b> = c x (x b_0 - c b_1)). Equivalent forms: by (H), r_{1,2}[m](0) - (x/c) r_{0,2}[m](0) =
(D_2 r_{0,2}[m-1])(0) = (r_{0,2}[m-1])'(0), because L_2(0) = x/c and Delta_2'(0) = 0; and D_2 f = Delta_2^(1/2) (Delta_2^(-1/2) f)'.
(c) By (1.1), P_22(z, w)|_{t=0} = -h_w(-z)^-2 h_w(-w)^-2 <(1, x/c), (1, x/c)> = -(x^3/c) h_w(z)^2 h_w(w)^2, since <(1, x/c), (1, x/c)> =
x^2 (2x/c) - c x (x/c)^2 = x^3/c and h_w(u) h_w(-u) = 1; and P_22(-w, w)|_{t=0} = -x^3/c. QED
[S2_checks_reg.pl part T and S4_checks_r1.pl B4 ((a), (b) and [t^1] r_{0,2} = 0 for z^0..z^6 at four x), B5 ((c) for m + n <= 5).]

**Corollary 8.7 (K_2 = 1; V6 in full).** K_2(z) := <v_2(z), v_2(-z)>/Delta_2 lies in O[[z]], is t-independent, and K_2 = r_{0,2}(z; 0) r_{0,2}(-z; 0) = 1.
With Proposition 4.8(c) and Corollary 8.5: <v_i(w), v_j(-w)> = delta_{ij} Delta_i for all i, j, to all orders in w. (Not used for Theorem LC.)
Proof. vn_2 := nu_2 v_2 with nu_2 := 1/(x U_2 f2r^(1/2)) (a unit of O, nu_2^2 = 1/Delta_2, nu_2'/nu_2 = -Delta_2'/(2 Delta_2)) lies in O[[z]]^2 and
satisfies z vn_2' = (A - L_2) vn_2 (Theorem 8.4). As in Proposition 4.8(a), K_2 = <vn_2(z), vn_2(-z)> has zero t-derivative, so
K_2 = <v_2(z; 0), v_2(-z; 0)>/Delta_2(0) = r_{0,2}(z; 0) r_{0,2}(-z; 0) (x^3/c)/(x^3/c) by Proposition 8.6(b), and h_w(z)^-2 h_w(-z)^-2 = 1. QED
K_2 = 1 + O(z^2) is automatic (K_2 is even in z, K_2(0) = <chi_2, chi_2>/Delta_2 = 1); the vanishing of the higher coefficients is the content of
the analytic input (a). E_22 and Dil_2 are regular whatever K_2 is. [Mathematica (V6) for i = j = 2; S4_checks_r1.pl B7 (z^0..z^6, t^0..t^6).]

**Remark 8.8 (E_22 at t = 0 is not of rank one; an erratum to conifold_regularity_note.md F3).** (1) Write h_w(z)^2 = exp(Sum_k s_k z^(2k-1)),
s_k := 2 bk_k St_k (s_1 = St_1/6, s_2 = -St_2/180, St_1 = (1 + x + x^2)/(x c)). Expanding Proposition 8.6(c),
(1 - h_w(z)^2 h_w(w)^2)/(z + w) = -s_1 - (s_1^2/2)(z + w) - s_2 (z^2 - z w + w^2) - (s_1^3/6)(z + w)^2 + (degree >= 3), so at t = 0
E_22(0, 0) = -(x^3/c) s_1 = -x^2 (1 + x + x^2)/(6 c^2), E_22(1, 0) = E_22(0, 1) = -(x^3/(2c)) s_1^2, E_22(1, 1) = (x^3/c)(s_2 - s_1^3/3), and
    [E_22(0,0) E_22(1,1) - E_22(0,1) E_22(1,0)]_{t=0} = (x^3/c)^2 (s_1^4/12 - s_1 s_2) = (x^6/c^2) St_1 (72 St_2 + 5 St_1^3)/77760,
a non-zero rational function of x (3/4096 at x = 1, where E_22(0,0) = -1/8 happens to equal E_11(0,0)). (2) conifold_regularity_note.md F3 states
"E_11(z, w) = E_22(z, w) = -(x^3/2) phi(z) phi(w)" at t = 0. Its E_22 part is the E_11 table read through the conjugate encoding:
hvc_vertex_identity.wls line 14 prints val0[EdN[2,2][m, n], U2] = ev_2(cconj EdN[1,1][m, n]) = ev_1(EdN[1,1][m, n]) (Lemma 8.1(ii)), i.e.
-(x^3/2)/c^(m+n+2). Evaluated at U1 it gives (c). The graph sums are not affected: in graphValueC, cconj maps the term of a colouring to the term
of the complementary colouring, so the sum over all colourings is cconj-invariant (cSymQ in orbitTerm) and has the same value at U_1 and U_2.
Theorem A of the note uses only F1 and the regularity of the factors; Theorem B uses the t = 0 rank-one forms of E_11 and E_12; the Reduction
Lemma uses Conjecture V (every t-order of the singular legs and loops). None of them uses E_22 beyond its regularity, so the E_22 erratum affects
none of them. (The Reduction Lemma has a gap of its own, for E_11 edges between two singular vertices of genus >= 2: Remark 13.1. Corollary R
does not use it.) (3) The same artifact sits in the note's claim (F3) that E_22 has a double pole in the point frame and in frame A: line 24 of
hvc_asymptotics.wls and hvc_asymptotics_g4.wls evaluated EdN[2,2] at U2 and so printed the orders of E_11. Read at U1, E_22 is O(t^0) in the point
frame and O(t^2) in frame A (there v_2(z; 0) = (1, x/c) and [t^1] v_2 = 0, so E_22|_{t=0} = 0); only E_11 has the double pole. By Lemma 8.2 and
the t-independence of the frame changes, E_22 is analytic at t = 0 in every frame.
[referee_scripts/s4ref/ref_indep.pl: agreement with (c) at 21 of 21 entries for each of six x, with the note's value at 0 of 21 (1 of 21 at
x = 1); S4_checks_r1.pl B5 at x = 3/7, 5/3, -5/2, 2/9. Mathematica, hvc_leaf_proof_check2.wls (M14), symbolic x: EdN[2,2] = cconj EdN[1,1] exactly,
(c) for every m + n <= 8, E_22(0,0), the minor and its value 3/4096 at x = 1; read at U2 the table is -(x^3/2)/c^(m+n+2), the E_11 table;
(Q2): the orders in the point frame and in frame A.]

**Remark 8.9 (the regular column is determined by the singular data; not needed below).** Corollary 8.5 at order t^k (k >= -1), the H-row in the
form x r_{0,2} - c r_{1,2} = (x - c L_2) r_{0,2} - c z D_2 r_{0,2}, <(0, c), b> = c x (x b_0 - c b_1), and x - c L_2 = -x^2 tau s*^2 = O(t) give
    c (k + 2) z [t^(k+2)] r_{0,2} = [t^(k+1)]((x - c L_2) r_{0,2}) + c z [t^(k+1)](r_{0,2} Delta_2'/(2 Delta_2)) + (1/(c x)) Sum_{k1=0}^{k} <sigma^(k1), v_{2,k-k1}>,
v_{2,k2} := [t^(k2)] v_2(z), whose right side involves r_{0,2} only up to order t^(k+1). So, given r_{0,2}(z; 0) = h_w^-2, the whole regular column, and
with it every E_12 leg and E_22 to all orders, is determined recursively by the closed-form singular data and by L_2, Delta_2; e.g. [t^1] r_{0,2} = 0.

**Remark 8.10 (Theorem 8.4 in the module language; the formal form).** In the notation of Remark 4.10, with M((z))-modules Mod^M_lambda :=
M((z)) eps_lambda and d/dtau := (z (1 - tau))^-1 nabla on them (Leibniz rule d/dtau (a Y) = a_tau Y + a d/dtau Y; on Mod^M_0 it is d/dtau): put
mu_2 := f0' = 1 - x s~, f2^(1/2) := tau f2r^(1/2), LE_2 := f2^(-1/2) L2 eps_{mu_2} and Pf := tau eps_mu. Then mu + mu_2 = L_2 (8.1), and Step A says
Sm^(2)_0 := nu_2 r^pt_{0,2} eps_{L_2} = (1/x) Pf LE_2; Step C says d/dtau LE_2 = f2^(-1/2)(Dop L2) eps_{mu_2} and Gauss[LE_2] = 0 (this is S2's (B2));
Step D is Lemma 4.5 in module form, nabla(Pf Y) = Pf (nabla + omega) Y and [nabla^2 - (c/tau) nabla + x/tau](Pf Y) = ((1 - tau) z^2/tau) Pf Gauss[Y];
so nabla Sm^(2) = A Sm^(2). Read with honest exponentials, e^{u_2/z} Delta_2^(-1/2) r^pt_{0,2} = kappa_2(z) Q^{x/z} tau^{c/z+1} NG_2, with
NG_2 := e^{f0/z} (-f2)^(-1/2) L2 the full (unnormalised) expansion of the Euler integrand at s~ and the t-independent constant
kappa_2(z) := (-x^2)^(-1/2) e^{x log(-1)/z} = -(-x^2)^(-1/2) kappa_ph(z); in this form Step A is (P1)(ii) and Step C is (P2) in Lfull form. The z^-2
coefficient of (8.5) is L2[0] times Eik(f0_tau) (3.19), which vanishes because (1 - tau) Eik equals x (x tau s^2 - c s + 1) at s = s~: f0_tau = s~/U_2 is
the second root of the eikonal equation, -c/tau + O(1), associated with the exponent 1 - gamma. Accordingly f2^(-1/2) L2 has tau-order exactly -1
and LE_2 has exponent -1 - c/z = 1 - gamma, while Sm^(2) has exponent 0 (Remark 4.10).

## 9. Step 9: the E_12 and E_21 legs (L2 for j = 2; F3 for E_12, E_21)

### 9.1 Extraction of a leg
**Lemma 9.1.** Let E be in C[[psi, w]] and Num := (psi + w) E. For every b >= 0, in C((psi)),
    [w^b] E = Sum_{n=0}^{b} (-1)^n psi^(-n-1) [w^(b-n)] Num.                                                                       (9.1)
Proof. Put E_m := [w^m] E in C[[psi]] and E_{-1} := 0. Then [w^m] Num = psi E_m + E_{m-1}, so Sum_{n=0}^{b} (-1)^n psi^(-n-1)(psi E_{b-n} + E_{b-n-1})
= Sum_{n=0}^{b} (-1)^n psi^(-n) E_{b-n} - Sum_{n=1}^{b+1} (-1)^n psi^(-n) E_{b-n} = E_b. QED
Equivalently (9.1) is the expansion of Num/(psi + w) in C((psi))[[w]] (where psi + w is invertible, 1/(psi + w) = Sum_n (-1)^n w^n psi^(-n-1)); C[[psi, w]]
is contained in C((psi))[[w]]. Being linear, (9.1) holds coefficientwise in t.

### 9.2 L2 for j = 2 and the regularity of E_12
**Proposition 9.2.** Assume (R) of Lemma 8.2(ii) (every v_2[m] analytic at t = 0). Then:
 (i) rho E_12 and E_12 are O(t^0): every E_12(a, b) is regular at t = 0 (F3 for E_12);
 (ii) for all k >= 0 and b >= 0, [t^k](rho f_{2,b}) lies in Lam_{k+1}.
Proof. Put N12(psi, w) := rho P_12(psi, w) = -<sigma(psi), v_2(-w)>. Since rho depends on t only, (1.1) gives
    (psi + w) rho E_12(psi, w) = N12(psi, w) - N12(-w, w).                                                                       (9.2)
By Proposition 6.1 and (R), N12 is a Laurent series in t starting at t^-1, and each coefficient is a finite sum:
    [t^k] N12(psi, w) = -Sum_{k1=-1}^{k} <sigma^(k1)(psi), v_{2,k-k1}(-w)>,     v_{2,k2}(-w) := [t^(k2)] v_2(-w), a pair in C[[w]].           (9.3)
(i) [t^-1] N12(psi, w) = -<(0, c), v_{2,0}(-w)> does not depend on psi, so the t^-1 coefficient of the right side of (9.2) is zero; C[[psi, w]] has no
zero divisors, so [t^-1](rho E_12) = 0, and lower orders do not occur. So rho E_12 = O(t^0), and E_12 = rho^-1 (rho E_12) = O(t^0).
(ii) Fix k, b >= 0. Apply (9.1) to E = [t^k](rho E_12) and Num = [t^k](N12(psi, w) - N12(-w, w)):
    [t^k](rho f_{2,b})(psi) = Sum_{n=0}^{b} (-1)^n psi^(-n-1) Ym_{b-n}(psi),                                                            (9.4)
    Ym_m(psi) := [w^m][t^k](N12(psi, w) - N12(-w, w)) = -(-1)^m Sum_{k1=-1}^{k} <sigma^(k1)(psi), v_{2,k-k1}[m]> - ym_m,
where v_{2,k2}[m] := [t^(k2)] v_2[m] is in C^2 ([w^m] v_2(-w) = (-1)^m v_2[m]) and ym_m := [w^m][t^k] N12(-w, w) is a constant. By Proposition 6.1
both components of every sigma^(k1) (-1 <= k1 <= k) lie in Lam_{k1+1}^+, contained in Lam_{k+1}^+; hence Ym_m is the Taylor series of a rational
function Ymh_m in Lam_{k+1}^+. Consider Rat(psi) := Sum_{n=0}^{b} (-1)^n psi^(-n-1) Ymh_{b-n}(psi).
 (a) Away from psi = 0, Rat has at most simple poles, located at the points c/j with 1 <= j <= k + 1.
 (b) Rat(psi) -> 0 as psi -> oo, because each Ymh is bounded at oo and psi^(-n-1) -> 0.
 (c) The Laurent expansion of Rat at psi = 0 is the right side of (9.4), i.e. [t^k](rho f_{2,b}), which lies in C[[psi]] because E_12(psi, w) is a
     power series in psi (Lemma 1.1). Hence the principal part of Rat at 0 vanishes: the apparent pole of order <= b + 1 cancels.
So Rat - Sum_{j=1}^{k+1} Rat_j lam_j, Rat_j := lim_{psi -> c/j}(c - j psi) Rat(psi), has no poles in C and tends to 0 at oo, hence is 0:
[t^k](rho f_{2,b}) = Rat lies in Lam_{k+1}. QED
Explicitly, since psi^(-n-1) lam_j = (j/c)^(n+1) lam_j + (a polynomial in 1/psi without constant term), Rat_j = Sum_{n=0}^{b} (-1)^n (j/c)^(n+1) [lam_j] Ymh_{b-n}.
What is used about the regular rows: only (R). The rows v_2 enter (9.4) through the constants v_{2,k2}[m] and ym_m; their z-dependence, their values,
symplecticity and the orthogonality play no role. Analyticity is exactly what the degree bound needs: if some v_2[m] had a pole of order p at
t = 0, (9.3) would bring sigma^(k+p) into the order t^k, i.e. lam_j with j up to k + p + 1.

### 9.3 The leaf form of the E_12 edge
**Proposition 9.3.** Put S_j(t) := Sum_{k>=j-1} t^k S_j^(k) with the coefficient vectors (6.4); S_j = O(t^(j-1)) (closed form in Corollary 11.2(f)). Then
    sigma(psi) = (0, c/t) + Sum_{j>=1} S_j(t) lam_j(psi),                                                                            (9.5)
    rho E_12(psi, w) = -Sum_{j>=1} j <S_j(t), v_2(-w; t)> lam_j(psi) lam_j(-w),     lam_j(-w) = 1/(c + j w),                                (9.6)
identities in C[[psi, w]]((t)) in which only j <= k + 1 contribute at order t^k (by (R)). Consequently, for all k and b,
    [t^k](rho f_{2,b})(psi) = -Sum_{j=1}^{k+1} j lam_j(psi) Sum_{k1=j-1}^{k} < S_j^(k1), [w^b] ( v_{2,k-k1}(-w)/(c + j w) ) >.               (9.7)
Proof. (9.5) is Proposition 6.1. Insert it in N12. The term (0, c/t) contributes the same function of w to N12(psi, w) and N12(-w, w), so by (LAM6)
N12(psi, w) - N12(-w, w) = -Sum_j <S_j, v_2(-w)>(lam_j(psi) - lam_j(-w)) = -(psi + w) Sum_j j <S_j, v_2(-w)> lam_j(psi) lam_j(-w). Divide by psi + w
((9.2)); (9.7) is the coefficient of t^k w^b. QED
In (9.6) the cancellation at psi = 0 of step (c) in the proof of Proposition 9.2 is explicit, and the subtracted term N12(-w, w), which the orthogonality would set to
zero, drops out without being evaluated. Proposition 9.2(ii) is the statement "j <= k + 1" in (9.7).

### 9.4 The transposition E_21 = E_12^T; regularity of E_21
**Proposition 9.4.** Assume Corollary 8.5 (orthogonality). Then E_21(z, w) = E_12(w, z), i.e. E_21(b, a) = E_12(a, b) for all a, b. Hence E_21 is
regular at t = 0, and an edge joining a singular vertex to a regular one gives the singular vertex the leg f_{2,b} (b the power at the regular end)
whichever end comes first in the edge list (graphValueC attaches EdN[c(v1), c(v2)][a_{v1}, a_{v2}] to the edge (v1, v2)).
Proof. ph12(u) := P_12(-u, u) = -<v_1(u), v_2(-u)> = 0, so E_12(z, w) = P_12(z, w)/(z + w). By the symmetry of < , >, P_21(z, w) = -<v_2(-z), v_1(-w)> =
P_12(w, z) and P_21(-w, w) = P_12(w, -w) = ph12(-w) = 0. Hence E_21(z, w) = P_12(w, z)/(z + w) = E_12(w, z). QED
What the orthogonality is needed for. Without it, E_21(z, w) - E_12(w, z) = [ph12(z) - ph12(-w)]/(z + w), and by (9.1) the singular-side leg of an
E_21 edge would be [z^b] E_21(z, psi) = rho^-1 [Sum_{n=0}^{b} (-1)^n psi^(-n-1) [z^(b-n)] N12(psi, z) - (-1)^b psi^(-b-1) N12(psi, -psi)], which contains
N12(psi, -psi) = -<sigma(psi), v_2(psi)>, a pairing of sigma with the full power series v_2(psi) (at t = 0, h_w(psi)^-2 (1, x/c), a divergent Stirling
series); the rationality argument of Proposition 9.2 does not apply to it. With the orthogonality the leg is f_{2,b} and is rational. Even the
regularity of E_21 needs the t^-1 part of the orthogonality: since rho P_21(z, w) = -<v_2(-z), sigma(w)> and sigma^(-1) = (0, c),
    [t^-1](rho E_21)(z, w) = [fe(z) - fe(-w)]/(z + w),     fe(z) := -c x (x r_{0,2}(-z; 0) - c r_{1,2}(-z; 0)),
and no lower powers of t occur (Lemma 8.2). This vanishes iff fe is constant; the constant is fe(0) = -c x (x - c L_2(0)) = 0 (r_{0,2}[0] = 1,
r_{1,2}[0] = L_2, L_2(0) = x/c). So E_21 is regular iff x r_{0,2}(z; 0) = c r_{1,2}(z; 0), which is Proposition 8.6(b) (from Theorem 8.4 alone).
[Mathematica: E21(b, a) = E12(a, b) for a + b <= 8 (leaf_structure.md, Step 3); referee_scripts/s4ref/ref_indep.pl.]

### 9.5 All legs; unrescaled ends
**Corollary 9.5.** For all k, b >= 0 and j in {1, 2}: [t^k] f_{j,b} and [t^k](rho f_{j,b}) lie in Lam_{k+1}; [t^k](rho^e E_11(psi, psi')) lies in
span{lam_i(psi) lam_l(psi') : i, l >= 1, (i-1) + (l-1) <= k} for every integer e; all these series are O(t^0).
Proof. j = 2: Proposition 9.2 and multiplication by rho^-1 in Kx[[t]] ((T): [t^k](u X) = Sum_i u_i [t^(k-i)] X, and the spaces grow with k). j = 1 and
the loops: Corollary 7.4. QED

### 9.6 Use in Theorem R
At a singular vertex v of genus h_v >= 2, homogeneity rescales every half-edge at v exactly once (one factor rho per half-edge at v; section 13).
The factors attached to v depend on the other end of the edge:
 - the other end is a regular vertex: the edge factor E_12 (or E_21) gives v the leg rho f_{2,b}, b the power at the regular end (Proposition 9.2;
   either orientation by Proposition 9.4);
 - the other end is a singular vertex v' of genus 0 or 1 (not rescaled): the edge factor is rho E_11 = rho^-1 (rho^2 E_11), and for each power b' at v'
   it gives v the leg rho f_{1,b'} (Corollary 7.4), of leaf degree D <= k at order t^k; either orientation, since E_11 is symmetric;
 - the other end is a singular vertex v' of genus >= 2 (also rescaled): the edge factor is rho^2 E_11, which splits its D <= k between the two ends
   (Theorem 7.3);
 - a self-loop at v: the factor is rho^2 E_11 (Theorem 7.3).
At a vertex v' of genus 0 or 1 itself only regularity is needed: Delta_1^(h-1) is O(t^2) for h = 0 and equals 1 for h = 1, and every factor attached
to v' is O(t^0) (F3, section 10). This is precisely the parenthetical remark of Theorem R (ii).

## 10. Proof of Theorems LC and Q; F1-F3 in one table

Proof of Theorem LC.
 F1: t^2 rho^2 Delta_1 = -x c^3 and rho = tau/(t alg) with rho(0) = 1: Proposition 4.4(b), (d). Delta_2 = x^3 delta s*^2 analytic with Delta_2(0) = x^3/c,
     Delta_2'(0) = 0: Proposition 4.3(d), (e), Lemma 8.2(iii).
 F2: r_{0,1}(z; 0) = c/(c + z): Corollary 3.18. rho r_{1,1} = c/t + x/(c + z) + O(t): Proposition 6.1(a), (b) (sigma^(-1) = (0, c), sigma^(0)_1 = x lam_1)
     at psi = -z.
 F3: see the table below.
 L1: Theorem 5.4.
 L2: j = 1 by Corollary 7.4, j = 2 by Proposition 9.2(ii); the unrescaled legs [t^k] f_{j,b} by Corollary 9.5; E_21(b, a) = E_12(a, b) by Proposition 9.4.
 L3: Theorem 7.3; symmetry s.7.1; all powers rho^e by Corollary 7.4.
 Sharpness: Corollaries 11.3 and 11.5. Large-z form: Lemma 5.2 (lim z r_{0,1} = c/rho) and Lam_{k+1} vanishing at oo (Corollary 9.5). QED
Proof of Theorem Q. (i) Proposition 4.6 and Theorem 8.4. (ii) Proposition 4.8(c), Corollary 8.5, Corollary 8.7. (iii) Proposition 8.6. QED

F3, factor by factor. For F_{g,0} graphValueC uses only deltaPow (vertex weights), DilN and EdN; LfN occurs only with marked points.
    factor  proof of regularity at t = 0                                     value at t = 0
    Dil_1   Theorem 5.4, remark (iv) (q_t = rho^-1 (rho q_t))             q_0(psi) = -c psi/(c - psi)   (Dil_1(d; 0) = -c^(1-d))
    Dil_2   Lemma 8.2(iv)                                                 (-1)^n [z^(n-1)] h_w(z)^-2   (Proposition 8.6(a))
    E_11    s.7.4 (t^2 rho^2 E_11 = O(t^2), rho^-2 in Kx[[t]])               -(x^3/2) lam_1(psi) lam_1(psi')   (Theorem 7.3)
    E_12    Proposition 9.2(i)                                            lam_1(psi) chi12(w), chi12(w) = -x^3 h_w(w)^2/(c + w)   (s.11.3(b))
    E_21    Proposition 8.6(b) (Theorem 8.4) for the t^-1 part;          E_21(b, a) = E_12(a, b)
            Proposition 9.4 (Corollary 8.5) for E_21 = E_12^T
    E_22    Lemma 8.2(iv)                                                 (x^3/c) [1 - h_w(z)^2 h_w(w)^2]/(z + w)   (Proposition 8.6(c); not rank one)
Remark. By (9.5) the t^-1 part of sigma is the psi-independent covector (0, c). Hence the leg factors LfN[1][n] = (-1)^n r_{1,1}[n] are regular for
n >= 1, while LfN[1][0] = L_1 has a simple pole; they do not enter F_{g,0}.

Where each external input is used. Proposition 3.15 (the analytic input) enters through Theorem 3.17 into the KEY FORMULA (4.3) and from there into
everything that uses c/(c + z): F2, L1 at t^0 (q_0), sigma^(-1) = (0, c), K_1 = 1, [t^0] Ebar = 0, F3 for E_11, the t^0 rank-one loop. Lemma 3.8 is
used only for Beta integrals. V1 at U_1 (Proposition 4.3(c)) enters Proposition 4.4 (rho = tau/(t alg)); V1 at U_2 enters (8.4). The formal
Laplace calculus (P1-P3, homogeneity) enters Proposition 3.10, Proposition 3.12 and Theorem 8.4.

## 11. Closed forms (section 4 of the outline)

Nothing in Theorem LC depends on this section except the sharpness statements and the explicit values quoted as examples.

### 11.1 The dilaton leaf coefficients
**Theorem 11.1.** For j >= 1 put
    C_j := P_{j-1}/(j! (j-1)! c^(2j-4)) = c^2 [(1 - j/c)_{j-1}]^2/(j! (j-1)!)    (C_1 = c^2),
    zeta(t) := tau e^(dphi(tau)/c) = c (1 - U_1) U_1^(-(1-x)/c) = t + (1/c^2 - 1/2) t^2 + O(t^3),
    A_j(t) := (-1)^j C_j t^-1 zeta(t)^j 2F1(j x/c, j x/c; j; tau),   so A_j = O(t^(j-1)).
Then a_j^(k) = [t^k] A_j(t) for all j >= 1, k >= 0. Equivalently, as formal series in t with finitely many j at each order,
    z G(z; t) = c + Sum_{j>=1} A_j(t)/(c + j z),    rho q_t(psi) = c + Sum_{j>=1} A_j(t) lam_j(psi),    g(psi) = -(1/c) Sum_{j>=1} j A_j(t) lam_j(psi).
Proof. Fix j >= 1 and k >= 0. By Proposition 5.3, [t^k](z G) has at -c/j at most a simple pole, contributing a_j^(k)/(c + j z), so
a_j^(k) = j Res_{z=-c/j} [t^k](z G). Write z G = (tau/t) Sum_n tau^n (c z/(c + z)) Ephi F_n; [t^k](z G) is a finite sum over n <= k of terms
(scalar)(c z/(c + z)) z^-p F_n(z) (each [t^m] Ephi is a polynomial in 1/z, and evaluation at z = -c/j commutes with taking t-coefficients).
j = 1: F_n is regular at -c (its poles are at 0 and -c/l, l >= 2), so Res_{z=-c}[t^k](z G) = [t^k](-c^2 (tau/t) e^(dphi/c) F(-c, tau)), and at z = -c,
alpha = 1 - 1/c = x/c and gamma = 1: F(-c, tau) = 2F1(x/c, x/c; 1; tau) ((x/c)_n = prod_{i=1}^{n} (i - 1/c), Sum_m X_m(-c) tau^m = e^(dphi/c)).
j >= 2: c z/(c + z) and z^-p are regular at -c/j; F_n is regular there for n < j - 1 and has a simple pole for n >= j - 1, from the single factor
c + j z. At z = -c/j, 1 + l z = (j - l c)/j, c + l z = c (j - l)/j, z^n = (-c/j)^n and prod_{l=2, l != j}^{n+1} (j - l) = (j - 2)! (-1)^(n+1-j) (n + 1 - j)!, so
    Res_{z=-c/j} F_n = (-1)^(j-1) prod_{l=1}^{n} (j - l c)^2/(j^2 c^(2n-1) n! (j-2)! (n+1-j)!) = (-1)^(j-1) c [(1 - j/c)_n]^2/(j^2 n! (j-2)! (n+1-j)!),
using prod_{l=1}^{n} (j - l c) = (-c)^n (1 - j/c)_n. Reindex n = j - 1 + m and use (1 - j/c)_{j-1+m} = (1 - j/c)_{j-1} (j x/c)_m (since
1 - j/c + j - 1 = j x/c) and (j - 1 + m)! = (j - 1)! (j)_m:
    Sum_{n>=j-1} tau^n Res_{z=-c/j} F_n = (-1)^(j-1) c [(1 - j/c)_{j-1}]^2 tau^(j-1)/(j^2 (j-2)! (j-1)!) 2F1(j x/c, j x/c; j; tau)
(the limit formula DLMF 15.2.3_5 for the regularised 2F1 at gamma = 2 - j). The remaining factors at z = -c/j are c z/(c + z) = -c/(j - 1) and
Ephi = e^(j dphi/c). Multiplying by j and using j (j - 1)(j - 2)! = j!, a_j^(k) = [t^k] A_j. The three generating functions follow from Proposition 5.3,
Theorem 5.4 and g_j^(k) = -(j/c) a_j^(k). QED
(For fixed small t, z -> z G(z; t) is meromorphic near -c/j and A_j(t) = j Res_{z=-c/j} z G(z; t) as functions; the 2F1 series converges uniformly on
a small circle around -c/j, and its poles accumulate only at z = 0. Not needed.)

**Corollary 11.2.** (a) Sum rule: Sum_{j>=1} A_j(t) = -c^2, i.e. Sum_{j<=k+1} a_j^(k) = -c^2 delta_{k0} (no pole of G at z = 0); equivalently
    Sum_j (-1)^j C_j zeta^j 2F1(j x/c, j x/c; j; tau) = -c^2 t = c^2 log(1 - tau).
    Check at x -> 0: c = 1, C_j = 1/j, 2F1 = 1, s* = 1, zeta = tau/(1 - tau), and Sum_j (-1)^j zeta^j/j = -log(1 + zeta) = log(1 - tau).
 (b) With Ahat_j := t A_j = (-1)^j C_j zeta^j 2F1(j x/c, j x/c; j; tau) (a function of tau), t z G(z) = -c log(1 - tau) + Sum_j Ahat_j(tau)/(c + j z).
 (c) e^(dphi) = (c s*)^c (1 - tau s*)^(x-1) = (c s*)^c U_1^(x-1): by the critical equation 1 - s* = x s* (1 - tau s*), dphi = log s* + x log(1 - s*)
     - log(1 - tau s*) + log c - x log(x/c) = c log(c s*) + (x - 1) log U_1 (principal branches, quantities near 1). Hence zeta = c (1 - U_1) U_1^(-(1-x)/c).
 (d) With u := 4 x tau/c^2 and principal roots of 1 - u (near 1),
        rho = (tau/t) 2/((1 - u)^(1/4) (1 + (1 - u)^(1/2))),     alg = (1 - u)^(1/4) (1 + (1 - u)^(1/2))/2,
     a form without branch ambiguity, also for c < 0 (it follows from rho^2 = (tau/t)^2 s*^2 c^3/delta, delta = c (1 - u)^(1/2), s* = 2/(c (1 + (1 - u)^(1/2)))
     and rho(0) = 1). Written as (tau/t) s* (c^3/delta)^(1/2) the root must be c (1 - u)^(-1/4); for c < 0 the principal root of c^3/delta has value
     |c| = -c at t = 0 and would give rho(0) = -1. Expansion: 2/((1 - u)^(1/4)(1 + (1 - u)^(1/2))) = 1 + u/2 + (11/32) u^2 + O(u^3), so
        rho = 1 - (1 - x)^2 t/(2c^2) + (1 - 8x + 15x^2 - 8x^3 + x^4) t^2/(6c^4) - (1 - x)^2 (1 - 20x + 60x^2 - 20x^3 + x^4) t^3/(24 c^6)
              + (1 - 52x + 493x^2 - 1624x^3 + 2365x^4 - 1624x^5 + 493x^6 - 52x^7 + x^8) t^4/(120 c^8) + O(t^5),
        Delta_1 = -(x c^3/t^2) [1 + (1 - x)^2 t/c^2 + (5 - 4x - 6x^2 - 4x^3 + 5x^4) t^2/(12 c^4) + (1 - x)^2 t^3/(12 c^2) + O(t^4)].
 (e) Theta = -x^3 tau s*^2 = -4 x^3 tau/(c + delta)^2 = x c (1 - x s*) - x (c s* - 1 = x tau s*^2) = -(x^3/c^2) t + (1 - x)^2 x^3 t^2/(2c^4) + O(t^3).
 (f) The H-component. By g = -(1/c) Sum_j j A_j lam_j, (LAM2) and Sum_j (t A_j)' = -c^2 (from (a)),
        sigma(psi) = (0, c/t) + Sum_{j>=1} S_j(t) lam_j(psi),
        S_j(t) = ( -(j/c) A_j,  A_j' - ((j - 1)/t) A_j - (j/c) l A_j ) = t^(j-1) ( -(j/c) Atil_j,  Atil_j' - (j/c) l Atil_j ),
     with Atil_j := t^(1-j) A_j = (-1)^j C_j (zeta/t)^j 2F1(j x/c, j x/c; j; tau) analytic, Atil_j(0) = (-1)^j C_j; the coefficient of t^k is (6.4);
     S_1(0) = (c, x) = c (1, x/c). [Derivation: -(psi/t) g = (1/t) Sum_j A_j lam_j - (1/(c t)) Sum_j A_j, -psi g' = Sum_j A_j' lam_j - (1/c) Sum_j A_j', and the
     constants add up to -(1/(c t)) Sum_j (t A_j)' = c/t.]
 (g) Dilaton consistency (a relation between Theorem 11.1 and (d), not used in any proof). Dil_1(1) = -r_{0,1}[0] = -1 at every t, so
     [psi^1](rho q_t) = -rho; with [psi^1] lam_j = j/c^2 this reads
        Sum_{j>=1} j A_j(t) = -c^2 rho(t),   i.e.   Sum_{j<=k+1} j a_j^(k) = -c^2 [t^k] rho   for every k.
     [Mathematica, hvc_leaf_proof_check3.wls (Q3): t^0..t^12, symbolic x.]

**Corollary 11.3 (top and subleading coefficients, all k).** For every k >= 0, resp. k >= 1,
    a_{k+1}^(k) = (-1)^(k+1) C_{k+1} = (-1)^(k+1) P_k/((k + 1)! k! c^(2k-2)),     g_{k+1}^(k) = (-1)^k P_k/((k!)^2 c^(2k-1)),
    a_k^(k) = (-1)^k k C_k (1 - x)^2/(2 c^2).
P_k = prod_{m=1}^{k} (k + 1 - m - m x)^2: (1 - x)^2, (2 - x)^2 (1 - 2x)^2, 4 (3 - x)^2 (1 - x)^2 (1 - 3x)^2, ...; this is the top-coefficient formula of
leaf_structure.md, Step 3 (verified there for k <= 5), now proved for all k.
Proof. A_j = (-1)^j C_j t^(j-1) (zeta/t)^j 2F1(j x/c, j x/c; j; tau) with (zeta/t)(0) = 1 and 2F1 = 1 + (j x^2/c^2) tau + O(tau^2): the first claim is
j = k + 1 at t^k; g_{k+1}^(k) = -((k + 1)/c) a_{k+1}^(k). For the second, [t^1] of (zeta/t)^k 2F1(k x/c, k x/c; k; tau) is
k (1/c^2 - 1/2) + k x^2/c^2 = k (2 + 2x^2 - c^2)/(2 c^2) = k (1 - x)^2/(2 c^2). QED
Direct computation of the top coefficient (k >= 1), without Theorem 11.1: put z0 := -c/(k + 1). Only the term n = k of (5.1), and in H_k only F_k,
has a pole at z0 (d_{k,k} = 1), so a_{k+1}^(k) = (k + 1)(c z0/(c + z0)) Res_{z0} F_k with c z0/(c + z0) = -c/k, 1 + i z0 = (k + 1 - i c)/(k + 1),
z0^k = (-c)^k/(k + 1)^k, prod_{i=2}^{k} (c + i z0) = c^(k-1) (k - 1)!/(k + 1)^(k-1) and Res_{z0} 1/(c + (k + 1) z) = 1/(k + 1); hence
Res_{z0} F_k = (-1)^k P_k/((k + 1)^2 k! (k - 1)! c^(2k-1)). For k = 0, a_1^(0) = -c^2 = -P_0 c^2.
Sharpness: C_{k+1} != 0 unless c is in {(k + 1)/m : 1 <= m <= k}, a finite set; for generic x the top leaf lam_{k+1} occurs at order t^k, and the
bound j <= k + 1 of L1 is sharp.

### 11.2 The loops
**Theorem 11.4.** With Ahat_j = t A_j,
    rho^2 E_11(psi, psi') = Sum_{i, l >= 1} e_il(t) lam_i(psi) lam_l(psi'),
    e_il(t) = (i l/(c^2 t^2)) Int_0^t Theta(u) Ahat_i'(u) Ahat_l'(u) du = -(i l x^3/(c^2 t^2)) Int_0^tau v (1 - v) s*(v)^2 (dAhat_i/dv)(dAhat_l/dv) dv,
and e_il = O(t^((i-1) + (l-1))).
Proof. By Theorem 11.1, Dg(psi) = d/dt (t g(psi)) = -(1/c) Sum_j j Ahat_j' lam_j(psi). Insert this in (7.2) and integrate from 0 ((7.4)). For the second
form substitute v = tau: dt = dv/(1 - v), Ahat' = (1 - v) dAhat/dv and Theta = -x^3 v s*^2. The order follows from Theta = O(t) and Ahat_j' = O(t^(j-1)). QED

**Corollary 11.5 (anti-diagonal and corner, all k).** Since Theta = -(x^3/c^2) t + O(t^2) and Ahat_j' = j a_j^(j-1) t^(j-1) + O(t^j), for i + l = k + 2:
    e_il^(k) = -x^3 i^2 l^2 a_i^(i-1) a_l^(l-1)/(c^4 (i + l)) = -(-1)^(i+l) x^3 i^2 l^2 C_i C_l/(c^4 (i + l))     (= (7.6)).
With a_1^(0) = -c^2 this gives the corner e_{1,k+1}^(k) = ((k + 1)^2/(k + 2))(x^3/c^2) a_{k+1}^(k) of leaf_structure.md, Step 3, now proved for all k. The
anti-diagonal entries are non-zero for generic x, so the support bound (i-1) + (l-1) <= k of L3 is sharp entry by entry. Examples (as in
hvc_leaf_structure_out_kt5.txt): e_22^(2) = -(1 - x)^4 x^3/c^4, e_23^(3) = 3 (2 - x)^2 (1 - x)^2 (1 - 2x)^2 x^3/(10 c^6), e_33^(4) = -3 (2 - x)^4 (1 - 2x)^4 x^3/(32 c^8),
e_24^(4) = -4 (3 - x)^2 (1 - x)^4 (1 - 3x)^2 x^3/(27 c^8).

### 11.3 The legs
 (a) E_11 legs. [t^k](rho f_{1,b}) = Sum_{k1+k2=k} [t^(k1)](rho^-1) Sum_{i,l} e_il^(k2) (l^b/c^(b+1)) lam_i(psi): no lam_0, and its top coefficient (lam_{k+1} at
     t^k) is e_{k+1,1}^(k)/c^(b+1) = ((k + 1)^2/(k + 2)) x^3 a_{k+1}^(k)/c^(b+3), for all k and b. At t^0 the leg is -(x^3/2) c^(-b-1) lam_1.
 (b) E_12 at t = 0 (rank one): E_12(psi, w)|_{t=0} = lam_1(psi) chi12(w) with chi12(w) := -x^3 r_{0,2}(-w; 0)/(c + w) = -x^3 h_w(w)^2/(c + w). Proof: at t^0,
     (9.6) has only j = 1 with S_1^(0) = (c, x) = c (1, x/c), and <(1, x/c), b> = (x^3/c) b_0 for every b, so [t^0](rho E_12) = -x^3 r_{0,2}(-w; 0) lam_1(psi)
     lam_1(-w) and [t^0] E_12 = [t^0](rho E_12); then Proposition 8.6(a) (the orthogonality is not needed). Low orders: chi12_0 = -x^3/c,
     chi12_1 = -x^2 (1 - 5x + x^2)/(6c^2), chi12_2 = -x (1 - 10x + 63x^2 - 10x^3 + x^4)/(72 c^3), as in leaf_structure.md (ii).
 (c) E_12 top coefficients (all k, b), given Proposition 8.6(a), (b):
        [lam_{k+1}][t^k](rho f_{2,b}) = (k + 1)^2 (x^3/c^2) a_{k+1}^(k) [w^b] ( h_w(w)^2/(c + (k + 1) w) ).
     Proof: in (9.7) only j = k + 1 and k1 = k contribute, giving -(k + 1) <S_{k+1}^(k), [w^b](v_{2,0}(-w)/(c + (k + 1) w))>; v_{2,0}(-w) = h_w(w)^2 (1, x/c),
     <S, (1, x/c)> = (x^3/c) S_0, and S_0 = -((k + 1)/c) a_{k+1}^(k) by (6.4). Explicitly:
        b = 0: (k + 1)^2 x^3 a_{k+1}^(k)/c^3;   b = 1: (k + 1)^2 x^2 a_{k+1}^(k) (1 - (6k + 5) x + x^2)/(6 c^4);
        b = 2: (k + 1)^2 x a_{k+1}^(k) [(1 + x + x^2)^2 - 12 (k + 1) x (1 + x + x^2) + 72 (k + 1)^2 x^2]/(72 c^5);
     for k = 0, 1, 2 the last bracket is 1 - 10x + 63x^2 - 10x^3 + x^4, 1 - 22x + 267x^2 - 22x^3 + x^4, 1 - 34x + 615x^2 - 34x^3 + x^4, as in the tables.
     The lower coefficients [lam_j][t^k](rho f_{2,b}), j <= k, involve v_{2,k2} with k2 >= 1 ((9.7)); they are determined by Remark 8.9, but no closed
     form in b is claimed. In particular nothing forces the rank-one pattern of the lam_1 coefficients seen at k = 1, 2 (leaf_structure.md (ii)), and
     it fails at k = 3.

### 11.4 Low orders, by hand (all in agreement with leaf_structure.md)
 - rho, Delta_1, dphi, Theta: Corollary 11.2(d), (e) and s.1.5; dphi = t/c - x^2 t^2/(2c^3) + (1 - x)^2 x^2 t^3/(6c^5) - x^2 (1 - 8x + 12x^2 - 8x^3 + x^4) t^4/(24 c^7) + O(t^5).
 - q^(1) = ((1 - x)^2/2)(lam_2 - lam_1): G_1 = (c/(c + z))(H_1 - 1/2) = -(1 - x)^2/(2 (c + z)(c + 2z)), R_1 = ((1 - x)^2/2)(1/(c + 2z) - 1/(c + z)); also
   G_1 - rho_1 G_0 = (1 - x)^2 z/(c (c + z)(c + 2z)) with rho_1 = -(1 - x)^2/(2c^2), the O(t) row correction of leaf_structure.md, Step 0.
 - q^(2) = [-(2 - 4x + 3x^2 - 4x^3 + 2x^4) lam_1 + 6 (1 - x)^4 lam_2 - (2 - x)^2 (1 - 2x)^2 lam_3]/(12 c^2): a_3^(2) = -C_3 by Corollary 11.3;
   a_2^(2) = [(c z/(c + z))(c + 2z)(H_2 - H_1)]_{z=-c/2} = (1 - x)^4/(2c^2); a_1^(2) = [c z (H_2 - H_1 + 1/6)]_{z=-c} with H_1(-c) = (1 + x^2)/c^2,
   F_2(-c) = x^2 (1 + 2x)^2/(4c^4), (F_1 X_1)(-c) = x^2/c^4, X_2(-c) = 1/c^3, so 12 c^4 (H_2 - H_1 + 1/6)(-c) = 2 - 4x + 3x^2 - 4x^3 + 2x^4; the three sum to 0.
   Equivalently from Theorem 11.1: zeta/t = 1 + (1/c^2 - 1/2) t + (1/6 - 1/(2c^2) + (1 - x^2)/(2c^4)) t^2 + ..., 2F1(x/c, x/c; 1; tau) =
   1 + (x^2/c^2) t + (x^2 (1 + 2x)^2/(4c^4) - x^2/(2c^2)) t^2 + ...
 - q^(3) = (1 - x)^2 [-(3 + 5x^2 + 3x^4)/72 lam_1 + (7 - 20x + 22x^2 - 20x^3 + 7x^4)/24 lam_2 - (2 - x)^2 (1 - 2x)^2/8 lam_3 + (3 - x)^2 (1 - 3x)^2/36 lam_4]/c^4
   (top a_4^(3) = C_4 and subleading a_3^(3) = -3 C_3 (1 - x)^2/(2c^2) by Corollary 11.3; a_2^(3), a_1^(3) checked by S4_checks.pl).
 - g_1 = ((1 - x)^2/(2c))(lam_1 - 2 lam_2) = -((1 - x)^2/2) lam_1 lam_2;
   g_2 = (1/(12 c^3)) [(2 - 4x + 3x^2 - 4x^3 + 2x^4) lam_1 - 12 (1 - x)^4 lam_2 + 3 (2 - x)^2 (1 - 2x)^2 lam_3].
 - Theta_1 = -x^3/c^2, Theta_2 = x^3 (1 - x)^2/(2c^4) = 3 x c [t^3] dphi, Theta_3 = -x^3 (1 - 8x + 12x^2 - 8x^3 + x^4)/(6c^6) = 4 x c [t^4] dphi.
 - Loops: k = 0: e_11 = -x^3/2 (also from (7.1) directly: [t^2]<Y(psi), Y(psi')> = (x^3/2)(psi + psi') lam_1 lam_1', using x g_1 - c sigma^(1)_1 =
   -(x^2/(2c)) lam_1 and lam_1 + lam_1' = (2c - psi - psi') lam_1 lam_1'); k = 1: e_11 = -(1 - x)^2 x^3/(2c^2), e_12 = e_21 = 2 (1 - x)^2 x^3/(3c^2), e_22 = 0;
   k = 2: e_11 = -x^3 (7 - 20x + 21x^2 - 20x^3 + 7x^4)/(24 c^4), e_12 = e_21 = (1 - x)^4 x^3/c^4, e_13 = e_31 = -3 (2 - x)^2 (1 - 2x)^2 x^3/(16 c^4),
   e_22 = -(1 - x)^4 x^3/c^4; k = 3: non-zero exactly for (i, l) in {11, 12, 13, 14, 22, 23} and transposes.
 - Symplecticity at low orders (a test of (6.3) with l_0, l_1 against K_1 = 1): [t^1]<Y(psi), Y(-psi)> = 0 and [t^2]<Y(psi), Y(-psi)> = 0, using
   lam_j(psi) + lam_j(-psi) = 2c/(c^2 - j^2 psi^2): the coefficient of 1/(c^2 - psi^2) is x^2 + x (1 - x)^2 - x (1 - x + x^2) = 0, that of 1/(c^2 - 4 psi^2) is
   -2x (1 - x)^2 + 2x (1 - x)^2 = 0.
 - E_11 legs: t^0: -(x^3/2) c^(-b-1) lam_1; t^1: ((1 - x)^2 x^3/(12 c^(b+3)))[(2^(b+3) - 9) lam_1 + 8 lam_2] (i.e. -1/12, 7/12, 23/12, 55/12, 119/12, 247/12,
   503/12 times (1 - x)^2 x^3/c^(b+3) for b = 0..6 on lam_1); t^2, b = 0: lam_1: x^4 (4 - 3x + 4x^2)/(48 c^5), lam_2: (1 - x)^4 x^3/(3c^5),
   lam_3: -3 (2 - x)^2 (1 - 2x)^2 x^3/(16 c^5).
 - E_12 legs, b = 0, t^1: (1 - x)^2 x^3/c^3 [-(1/2) lam_1 + 2 lam_2].

## 12. What is proved and what is computer-verified

### 12.1 Status of every statement
Proved (complete arguments in this text; [f] = formal, [e] = elementary algebra plus analyticity of explicit algebraic functions of tau,
[a] = analytic):
 - Lemma 1.1 (edge formula) [e]; Lemmas 3.1, 3.2 [e]; Lemma 3.3, Definition 3.4, Propositions 3.5-3.7 (P3, P2, P1), Lemma 3.9 (homogeneity) [f];
   Lemma 3.8 (Laplace's method) [a, proved in full]; Proposition 3.10 and Remark 3.11 [f + e]; Proposition 3.12, Lemma 3.13, Theorem 3.14 [f];
   Proposition 3.15 [a: Lemma 3.8 + Stirling's series with remainder (cited: DLMF 5.11.1, s.5.11(ii); Whittaker-Watson s.12.33) for real x > 0,
   then rationality in x]; Corollary 3.16, Theorem 3.17, Corollary 3.18 [f, given Proposition 3.15]; Proposition 3.19 [f; existence of solutions
   of linear ODEs with analytic coefficients on a disc].
 - Lemmas 4.1, 4.2, Propositions 4.3, 4.4 (V1, V2, F1, rho, KEY FORMULA V3) [e]; Lemma 4.5 (V7) [e]; Propositions 4.6-4.8 (QDE at U_1, the
   operator A, K_1 = 1) [f, given section 3]; Remarks 4.9, 4.10 [f].
 - Sections 5, 6, 7 (L1, the sigma leaf structure, the E_11 flow V8, L3, L2 for j = 1, F3 for E_11) [f, given sections 3-4].
 - Lemmas 8.1-8.3 [e]; Theorem 8.4 (QDE at U_2) [f: the formal Laplace calculus at s~ and (3.16), coefficientwise at fixed t != 0];
   Corollary 8.5 (orthogonality) [f]; Proposition 8.6(a) [a: same method as Proposition 3.15], (b), (c) [f]; Corollary 8.7 (K_2 = 1) [f, given 8.6(a)];
   Remarks 8.8-8.10.
 - Section 9 (L2 for j = 2, F3 for E_12 and E_21, transposition) [f, given sections 6 and 8]; section 10; section 11 [f].
 - Proposition 13.2 (the frame-B formula, input (1) of Theorem R) [cited: Teleman, Givental, Pandharipande-Pixton-Zvonkine, Graber-Pandharipande,
   Coates-Givental, Mumford, Faber-Pandharipande; the identities (13.2), (13.3) [e]]; Lemma 13.3 (normalisation constants) [a: Lemma 3.8 + Stirling's
   series, as Proposition 3.15].
Not proved here and not needed: see section 14.

### 12.2 Computer checks (confirmations only; none is used in a proof)
(a) Mathematica, run before this proof was written, from the stored frame-B genus-4 factors hvc_factors_g4_Bm.mx (ZMAX = 9):
 hvc_leaf_proof_check.wls -> hvc_leaf_proof_check_out.txt (112 s): (1) V1 as an identity in QA; (2a) no pole at z = 0 in [t^k] G, k <= 6; (2b) the KEY
 FORMULA V3 for all t^k, k <= 6, and z^m, m <= 9; (2c) rho = tau/(t alg); (3), (3b) V2; (4) V5 at both idempotents, m = 1..9, t^0..t^6; (5) V6, exact in Q,
 w^n for n <= 9; (6) V8 at every (m, n) with m + n <= 8, t^0..t^6; (7) V4: [t^k](z G) for k <= 12, exact in x (regular at 0, limit c delta_{k0},
 denominator exactly prod_{j=1}^{k+1} (c + j z)); the printed dphi(t) and Theta(t) to t^7.
 hvc_pf_gauss.wls: V7 with a symbolic F, and Delta_1 = x^2 W_2 modulo the critical equation.
 hvc_leaf_structure.wls -> hvc_leaf_structure_out.txt, hvc_leaf_structure_out_kt5.txt (data hvc_leaf_structure_data.m, _kt5.m): the leaf tables of
 leaf_structure.md (dilaton k <= 5, legs, loops), E21(b, a) = E12(a, b) and E11(a, b) = E11(b, a) for a + b <= 8.
(b) Perl + Math::BigRat, exact rational arithmetic, in leaf_proof/ (scripts and outputs kept there):
    S1_checks.pl    C1-C6, 62 checks, 0 failures, 50 s (S1_checks_out.txt): Proposition 3.10 at t != 0 with the y-side computed as the code does;
                    Proposition 3.15 through z^7 at four x (one with x < 0); (3.17) mod (tau^5, z^5); (3.16) at 6 random points; the eikonal/transport
                    data; Lemma 3.2, (3.15) and the uniform version of Remark 3.11.
    S2_checks.pl    parts A-E, 943 checks, 0 failures, series to order 22 (S2_checks_out.txt): the identities of section 4 at rational points;
                    (4.4) on Laurent monomials; closed form -> H-row -> Row 2 -> K_1 = 1 as t-series; Gauss, the F-form of Row 2 and the Wronskian
                    identity as tau-series; rho, Delta_1, dphi, Theta against the Mathematica outputs.
    S2_checks_reg.pl parts P, T, 420 checks, 0 failures, 2.4 min (S2_checks_reg_out.txt): Lemma 8.3, Steps A and C of Theorem 8.4 as exact jets
                    at four points tau0 != 0 (with a negative control), V5 at both idempotents and V6 for all (i, j); Row 2 at U_2 as t-series (with a
                    negative control), v_2(z; 0) = h_w^-2 (1, x/c), K_2 = 1.
    S3_checks.pl    1540 checks at x = 2/7, 3/5, 7/3, -5/9, 0 failures (S3_checks_out.txt): the closed form against all dilaton tables t^0..t^5; the
                    residue and hypergeometric forms of A_j; sigma; symplecticity; (7.5) against all loop tables; Corollary 7.4 against all j = 1 legs.
    S4_checks.pl    2592 comparisons at x = 3/7, 2/5, 5/3, 7/2 (and 1314 at x = -1/3, 11/5), 0 mismatches (S4_checks_out.txt): section 11 against
                    the tables of leaf_structure.md (rho to t^5, q^(k) for k <= 5, all loops k <= 5, all E_11 legs, t^0 and top E_12 legs), and
                    partial fractions of the KEY FORMULA at rational z for k <= 11.
    S4_checks_r1.pl 901 comparisons, 0 mismatches, 163 s (S4_checks_r1_out.txt): from the code's definitions only (no closed form of the rows):
                    Step A at s~ at three points with t != 0; V1 at U_2 and (8.1); Row 2 at U_2 (z^1..z^6, t^0..t^6); Proposition 8.6(a), (b), (c);
                    the minor of Remark 8.8; Corollary 11.2(d) for c < 0; K_2 = 1.
(c) Independent referee scripts (Perl + Math::BigRat), all passing; written in the session scratchpad and archived, with their outputs, in
 leaf_proof/referee_scripts/ (paths below are relative to it; ref_indep.pl reads the project folder through its $ROOT):
    s1ref\ref1.pl, s1ref\ref2.pl      the code's point-frame row rebuilt from wickDirect as exact tau-series against h_w c/(c + z) F Ephi alg
                                    mod (tau^7, z^6) at x = 2, 1/3, 3/5, -5/2; rho, Delta_1 and the top leaf coefficients; the transport recursion; Lemma 3.2 for k <= 12.
    ref_S2_jets.pl (200 checks), ref_S2_key.pl (438), ref_S2_leaves.pl (1081): V1, V5 at both idempotents and V6 for all (i, j) at generic t0 != 0;
                                    the code's frame-B singular row as exact t-series against the KEY FORMULA (negative control: twisting by h_w
                                    instead of h_w^-1 fails); dilaton, K_1 = 1, loops, flow, legs.
    ref_s3_indep.pl (2868 checks)    at x = 3/7, -2/5, 5/2, 4/3: regularity of [t^k](F Ephi) and [t^k] G for k <= 8 with the negative z-powers kept;
                                    the code's evalQ formula for E_11 against the integrated flow (7.5), m + n <= 8, t^0..t^6; all stored tables.
    s3rev_spot.pl (162 checks)       residues, top coefficients, corners for k <= 6, the t^1 leg formula for b <= 7.
    s4ref\ref_indep.pl (7752 comparisons at x = 3/7, 5/3, -1/3, 11/5, -5/2, 1): all frame-B factors rebuilt from the code's definitions; roots,
                                    V1 at both roots, Proposition 8.6, Row 2 at U_2, orthogonality, K_1 = K_2 = 1, E_21(b, a) = E_12(a, b), E_22 at t = 0
                                    (agreeing with Proposition 8.6(c) at 21 of 21 entries per x, with the note's F3 at 0 of 21, 1 of 21 at x = 1), (6.4) and
                                    (9.7) against EdN.
(d) Spot checks of this assembly (referee_scripts/assemble/spot.pl, spot2.pl, spot3.pl; 135 checks, 0 failures): the unified forms of alg (including the
 branch-free one), rho to t^4, rho^-1 to t^2, the bracket of Delta_1 to t^3, Delta_2 = (L_2 - L_1)(x^2 - c x L_2) and its tau^2 term, the five forms of
 Theta, the three forms of dphi and its t^4 term, the two forms of zeta, l_0, l_1, and A_1..A_4 at low orders, at x = 2, 3/5, -5/2, 1/3, 7/3.
 The critic's spot check referee_scripts/crit_rank.pl (the E_11 block determinant of Remark 13.1 from the fitted leaf data) is archived there too.
(e) Mathematica, run after the assembly (2026-10-31); every check passed. Outputs in the project folder.
    hvc_leaf_proof_check2.wls 6 5 -> hvc_leaf_proof_check2_out.txt (63 s; stored factors, symbolic x): Q2 (point frame and frame A, genus-3
        factors: E_22 at U1 is O(t^0), resp. O(t^2), while EdN[2,2] read at U2 shows the double pole of E_11); M14 (conjugate encoding exact for
        rows, Dil, E_22, E_21; Proposition 8.6(c) for m + n <= 8; E_22(0,0); the minor of Remark 8.8); M15 (E_21 = E_12^T exactly, m + n <= 8; all
        Dil_i, E_ij of order >= 0); M8 (v_2(z; 0) = h_w^-2 (1, x/c) to z^9, [t^1] r_{0,2} = 0, Delta_2 = x^3 delta s*^2 to t^10); Q1 (Remark 13.1:
        det of the E_11 block over psi-orders {0, 1} is 0, 0, (1-x)^4 x^6/(18 c^10) at t^0, t^1, t^2); the branch-free rho of Corollary 11.2(d);
        M10b (Theorem 11.4 against EdN[1,1], m + n <= 8, t^0..t^6); M13 ((9.6) against EdN[1,2], m + n <= 8, t^0..t^5).
    hvc_leaf_proof_check3.wls -> hvc_leaf_proof_check3_out.txt (62 s; no stored factors): M5 ((3.16), (3.19), eikonal, P_tr, R_tr); M6 (Lemma 3.2,
        k <= 20); Q3 (Corollary 11.2(g), t^12); M7 (Remark 4.9 to tau^14, symbolic x and z); M4 (Proposition 3.19, m <= 5, tau^10); M3
        (Proposition 3.15 and Proposition 8.6(a) through z^15, from the code's wickDirect); M9 ((8.5) at z^-2..z^5 as Laurent series in tau,
        L2 computed in the coordinate sigma = tau s); M2 ((3.17) with N from the Wick formula with amplitude in the coordinate s, tau^8 z^9
        at x = 2, 1/3, -5/2, 3/7 and tau^5 z^6 symbolically); M10, M11 (E_11 built from the KEY FORMULA alone by the code's difference quotient,
        m + n <= 14, t^-2..t^8: no polar terms, (7.5) both lines, (7.6), Theorem 11.4; at x = 2, 3/7, -5/2).
    hvc_leaf_proof_check4.wls x0 12 8 -> hvc_leaf_proof_check4_out_2.txt, _3_5.txt (40 s each): M1 and M16 at rational x. The frame-B factors
        rebuilt from scratch by the code (setupWick[12], rmatrixAt0, twistRows) at x = 2 and x = 3/5: the KEY FORMULA for m <= 12, t^0..t^8;
        Row 2 of the QDE at both idempotents, m = 1..12, t^0..t^8; V6 exactly in Q for all i, j, w^n with n <= 12.
    M12 is covered by checks (2b), (7) of hvc_leaf_proof_check.wls and (a), (b) of hvc_leaf_closed_forms.wls.

### 12.3 Mathematica checks requested by the sections (all run, 2026-10-31: s.12.2(e))
None of these is needed for a proof; each extends a confirmation. The list is kept as the record of what was asked.
 M1  Fresh factors, not the stored file: setupWick[12]; rmatrixAt0[{a1 -> x, a2 -> 1, e1 -> x, e2 -> 1}]; twistRows[hCoeffs[{1 + x, -x, -1}, -1, 12]];
     ratCN[0,1][m] at U1 as a t-series to t^10 against (c/(c + z)) F Ephi alg, symbolic x, m <= 12; then checks (4), (5) of hvc_leaf_proof_check.wls
     (V5 at both idempotents, V6) to z^12 (symbolic x, or at x = 2 and x = 3/5).
 M2  (3.17) symbolically, tau^0..tau^8, z^0..z^9, with N from the Wick formula in the coordinate s.
 M3  Proposition 3.15 through z^15, symbolic x.
 M4  The transport recursion of Proposition 3.19 for m <= 5 against [z^m](Ephi F) to tau^10, including Phi_0 = 1/alg.
 M5  (3.16), (3.19), P_tr = delta and the closed form of R_tr by Together/PolynomialRemainder modulo x tau s^2 - c s + 1.
 M6  Lemma 3.2 for k <= 20: Together[tau^k W_k(U1)] analytic at tau = 0 with value (-1)^k phi_0^(k)(1/c).
 M7  Remark 4.9: tau^2 <PhiF(z), PhiF(-z)> = -x c (c^2 - z^2), symbolic x and z, to tau^12.
 M8  From hvc_factors_g4_Bm.mx: v_2(z; 0) = h_w(z)^-2 (1, x/c) to z^9 and [t^1] ratCN[0,2][m] = 0 (evaluated at U1); Delta_2 = x^3 delta s*^2 exactly.
 M9  (8.5) at symbolic x as Laurent series in tau for z^-2..z^4, with L2 from the Wick formula at s~ = (c + delta)/(2 x tau).
 M10 E_11 from the closed form alone via the evalQ difference quotient for m + n <= 14, t^0..t^8, against (7.5); Theorem 11.4 against EdN[1,1]
     (m + n <= 8, t^0..t^6) and loopResR in hvc_leaf_structure_data_kt5.m.
 M11 Corollary 7.4 for b <= 12, k <= 4 (the stored leg tables stop at b <= 7).
 M12 Theorem 11.1 symbolically against the lam_j coefficients of [t^k](rho z r_{0,1}) (k <= 6), and the sum rule to t^12.
 M13 (9.6) against EdN[1,2][a, b] (stored rows, evaluated at U1), a + b <= 8, t^0..t^5 (the referee's Perl check covers a + b <= 5).
 M14 Confirmation of Remark 8.8: zeroC[cadd[EdN[2,2][m, n], cscale[-1, cconj[EdN[1,1][m, n]]]]] and the same for EdN[2,1] against cconj[EdN[1,2]];
     val0[EdN[2,2][m, n], U1] against Proposition 8.6(c) for m + n <= 8 (E_22(0, 0)|_{t=0} = -x^2 (1 + x + x^2)/(6c^2)); line 14 of hvc_vertex_identity.wls
     re-run with U1 in place of U2.
 M15 E_21(m, n) = E_12(n, m) exactly (m + n <= 8), and no negative t-powers in DilN[2], EdN[2,2], EdN[1,2], EdN[2,1], all evaluated at U1.
 M16 (Optional.) V5 for i = 2 and V6 for i != j to z^9, t^0..t^8, symbolic x.

## 13. Consequence: the all-genus regularity theorem is conditional only on M_irr

Notation of proof_status_mixed_scale.md: Fpot_h denotes the genus-h potential of the Hodge CohFT Omega_w in dilaton-shifted coordinates (called F_h
there), Xi_h(j_1, ..., j_r) := d_{lam_{j_1}} ... d_{lam_{j_r}} Fpot_h(q_0) (a finite sum of Hodge integrals at each order, since the lam_j are power
series in psi), J_R := Sum j_i, D(R) := Sum (j_i - 1), and
    M_irr(h):  Xi_h(R) = 0 for every multiset R of integers >= 2 with D(R) <= 2h - 3.
In this section g denotes the genus.

### 13.1 Theorem R (proof_status_mixed_scale.md s.3)
**Theorem R.** Fix g >= 2. Assume F1, F3, Lemma C through order t^(2h-3) for every 2 <= h <= g (at the psi-orders that occur in genus-g graphs), and
M_irr(h) for every 2 <= h <= g. Then every coloured stable graph of genus g with at least one edge is regular at t = 0 after summing over decorations,
and the single-vertex graph has the single pole (-x c^3)^(g-1) Fpot_g(q_0) t^(2-2g). Since the conifold F_g is P_g Li_{3-2g}(Q)/(eps1 eps2) up to
Q-independent terms, whose only singular term is (2g-3)! t^(2-2g), this pole is gap_g t^(2-2g), and F4 holds in genus g as a consequence.
Proof (proof_status_mixed_scale.md s.3; recalled with the references of this text). (i) Homogeneity. At a singular vertex v of genus h >= 2 with
n_v half-edges the vertex factor is Delta_1^(h-1) d^(n_v) Fpot_h(q_t)[f_1, ..., f_{n_v}]: the vertex class is the t-independent Omega_w in frame B,
and all t-dependence sits in Delta, Dil and E. The dilaton equation gives Fpot_h(lambda q) = lambda^(2-2h) Fpot_h(q), so
d^n Fpot_h(q)[f_1, ..., f_n] = lambda^(2h-2) d^n Fpot_h(lambda q)[lambda f_1, ..., lambda f_n]; with lambda = rho and (4.2) the vertex factor is
(-x c^3/t^2)^(h-1) d^(n_v) Fpot_h(rho q_t)[rho f_1, ..., rho f_{n_v}]. (ii) Leaf expansion. By L1, rho q_t = q_0 + Sum_{k>=1} t^k q^(k) with q^(k) in
Lam_{k+1}; by L2, L3, Proposition 9.4 and Corollary 9.5 (section 9.6), every leg, loop and singular-singular edge attached to v is a t-series whose
order-t^k coefficient is a combination of leaves lam_j (j >= 1) of total degree D <= k (an E_11 edge between two rescaled vertices splits its
D <= k between its ends; at an end attached to a singular vertex of genus 0 or 1, which is not rescaled, the factor rho^-1 = 1 + O(t) does not raise
D). Taylor-expanding Fpot_h(q_0 + Sum_k t^k q^(k)), the coefficient of t^(2-2h+k) of the vertex factor is a finite combination of Xi_h(R) with all parts
j >= 1 and D(R) <= k. (iii) Vanishing. Strip the parts equal to 1 by the string identity Xi_h(1, R) = (J_R/c^2) Xi_h(R) (proof_status_mixed_scale.md
s.1.2 in the normalisation c = 1): Xi_h(1^p, R') = prod_{i<p} (J_{R'} + i) c^(-2p) Xi_h(R'). If R' is empty and p >= 1, the factor J_{R'} = 0 kills it
(Theorem C); if R' is non-empty, all its parts are >= 2 and D(R') = D(R) <= k, so Xi_h(R') = 0 by M_irr(h) whenever k <= 2h - 3. So a non-zero polar
term needs D_v >= 2h_v - 2 at every singular vertex of genus >= 2 with r_v >= 1 leaves. (iv) Counting. In a connected graph with an edge every vertex
has a half-edge, so r_v >= 1. The t-order of a non-zero term is >= Sum_{v singular, h_v >= 2} (2 - 2h_v + D_v) + 2 (number of singular genus-0
vertices) + 0 (genus-1 singular and all regular vertices, by F1 and F3) >= 0. (v) The single-vertex graph has no edges: its only term with R empty is
the order k = 0, (-x c^3)^(g-1) Fpot_g(q_0) t^(2-2g); all its other polar terms have R non-empty with D <= k <= 2g - 3 and vanish. Comparing the sum
over all graphs with the known polar part of F_g identifies this term with gap_g t^(2-2g). QED

### 13.2 What this text adds
**Corollary R.** Fix g >= 2 and assume M_irr(h) for every 2 <= h <= g. Then every coloured stable graph of genus g with at least one edge is regular at
t = 0 after summing over decorations, the single-vertex graph has the single pole gap_g t^(2-2g), and F4 holds in genus g. In particular:
 - genus 2: unconditional, and without any computer input (M_irr(2) is proved: proof_status_mixed_scale.md s.1.4);
 - genus 3 and 4: unconditional, computer-assisted (M_irr(3): the exact lattice data hvc_multiprobe_out.txt; M_irr(4): hvc_multiprobe_g4_out.txt for
   r <= 2 and hvc_multiprobe_g4b.wls for the seven points with r >= 3; proof_status_mixed_scale.md s.1.9 and addendum (1)); this is independent of
   the direct genus-4 graph sums, which already show regularity;
 - genus g >= 5: conditional only on M_irr(h) for 5 <= h <= g, i.e. on Lemma A of proof_status_mixed_scale.md s.2 (its case h = 4 is now settled).
Proof. Theorem LC gives F1, F3 and Lemma C (L1, L2 for j = 1, 2 and every b, L3) for every k and at every psi-order, so the hypotheses of Theorem R
other than M_irr hold in every genus. The orientation issue of a mixed edge (graphValueC attaches EdN[c(v1), c(v2)] in the order of the edge list,
so E_21 occurs whenever the regular vertex is listed first) is settled by Proposition 9.4: E_21(b, a) = E_12(a, b), and E_21 is regular. The
unrescaled ends are covered by Corollary 9.5. QED
Three remarks on the logic. (a) The leaves that occur at order t^k have integer scales j <= k + 1 (Theorem LC), so only lattice points of
Conjecture M enter; only its vanishing direction is used. (b) The ranges match exactly: the polar orders k <= 2h - 3 need D <= 2h - 3, the range of
M_irr(h); the bound D <= k of Lemma C is sharp (Corollaries 11.3, 11.5) and the first non-zero values of M (D = 2h - 2) sit exactly at the order t^0,
so there is no slack and no missing case. (c) The route uses neither topological recursion nor the dictionary P2: Lemmas B and D of
proof_status_mixed_scale.md s.2 are one possible route to Lemma A and are not needed for Corollary R.
(d) Theorem LC is strictly stronger than Lemma C as stated in proof_status_mixed_scale.md s.2 in two ways that Theorem R uses: it includes k = 0
([t^0](rho f_{j,b}) in Lam_1, [t^0](rho^2 E_11) in span(lam_1 lam_1'), [t^0](rho q_t) = q_0), which was previously supplied by the rank-one limits of
the note's F3 (whose E_22 part is wrong, Remark 8.8), and it includes E_21(b, a) = E_12(a, b) (Proposition 9.4), which the orientation of mixed
edges in graphValueC needs.

**Remark 13.1 (Corollary R rests on Theorem R, not on the Reduction Lemma).** The Reduction Lemma of conifold_regularity_note.md s.5a writes the
finite E_11 block of an edge between two distinct singular vertices as Sum_alpha f_alpha f_alpha^T with f_alpha in the column space of the legs
f_{1,b}, and then applies Conjecture V, whose legs are such columns; for this the coefficients c_alpha in f_alpha = E c_alpha must be regular at
t = 0. No such factorisation exists once the block has two rows: E = E C E with C regular makes E C idempotent, so rank E(t) = trace(E(t) C(t))
would be constant near t = 0, but the block has rank 1 at t = 0 (Theorem 7.3) and rank 2 for small t != 0:
    det( E_11(a, a') )_{a, a' in {0, 1}} = (1 - x)^4 x^6 t^2/(18 c^10) + O(t^3)
(from the loop coefficients of s.11: the t^2 term is -(e^(1)_12)^2 + e^(0)_11 e^(2)_22 in leaf form; it vanishes at x = 1). Such blocks first occur
at genus 4 (two genus-2 vertices joined by an edge, both of colour 1, psi-powers 0..4 at each end). Theorem R does not factorise the block: it
expands rho^2 E_11 in the t-independent leaves lam_i(psi) lam_l(psi') with coefficients e^(k)_il regular in t (Theorem 7.3), and splits the leaf
degree D <= k between the two ends (s.9.6, s.13.1 (ii)). So the gap is in the note's Reduction Lemma, not in Theorem R or Corollary R; the
sentence of proof_status_mixed_scale.md s.3 calling Theorem R "the Reduction Lemma ... written in leaf language" should be read in this sense.
[Mathematica: hvc_leaf_proof_check2.wls (Q1) from the stored factors, symbolic x, also the t^3 term; referee_scripts/crit_rank.pl from the
fitted leaf data at x = 2, 1/3, 3/7, -5/2 and x = 1.]

### 13.3 The inputs of Theorem R that are not part of Lemma C
Inputs (2)-(5) are proved in the cited notes, not here; input (1) is Proposition 13.2 below: (1) the frame-B Hodge-vertex graph-sum formula
for the conifold free energies (conifold_regularity_note.md,
preamble; hvc_lib.m), the basis of the whole project, with the vertex class Omega_w t-independent; (2) the dilaton equation for Fpot_h (homogeneity);
(3) Theorem C (conifold_regularity_note.md s.3a, by the string equation); (4) the string/dilaton reduction of Conjecture M to M_irr
(proof_status_mixed_scale.md s.1.2); (5) the known polar part of the conifold F_g (the definition of gap_g in F4).

**Proposition 13.2 (the frame-B formula).** Let T = (C^*)^5 act on Z = Tot(O(-1)^2 + O^2 -> P^1) with tangent weights wt_+ = {alpha0, alpha1,
alpha2, -eps1, -eps2} and wt_- = {-alpha0, alpha0 + alpha1, alpha0 + alpha2, -eps1, -eps2} at the fixed points p_+, p_- of the zero section
(rc_lib.m line 84; [BS] App. A), and let T_F be the subtorus alpha0 = 0, alpha_i = eps_i, which fixes the zero section P^1 pointwise, with normal
weights (eps1, eps2, -eps1, -eps2). For g >= 2 let F_g(Q) := Sum_{d>=0} Q^d GW^{Z,T_F}_{g,d[P^1]}, defined by virtual localisation on the proper
fixed loci M_{g,0}(Z, d)^{T_F} = M_{g,0}(P^1, d) ([BS] Def. 2.2); F_g is homogeneous of degree 2g - 2 in (eps1, eps2), so nothing is lost by setting
eps1 = x, eps2 = 1. Then F_g lies in Q(x)(Q), and in Q(x)(Q)
    F_g = Sum_Gamma |Aut Gamma|^-1 Sum_c Sum_dec prod_v Delta_{c_v}^(g_v - 1) Int Omega_w psi^a prod Dil_{c_v}(d) prod_e E_{c c'}(a_e, b_e)   (13.1)
over the connected stable graphs Gamma of genus g without legs and the colourings c: V -> {1, 2}, with the decorations of s.0.1 and (C10) (exponents
a >= 0 on the half-edges; at each vertex an ordered m_v-tuple of dilaton degrees d >= 2 with weight 1/m_v!, nd in vertexOptionsH4, hvc_lib.m lines
115-126; integrals over M_{g_v, n_v + m_v}); Omega_w := prod_{a=1}^3 Lambda^vee_{g_v}(w_a)/w_a^{g_v} with w = (c, -x, -1); and Delta_i, Dil_i, E_ij
the frame-B factors (C6), (C8)-(C11), computed in QA (C3). At the true idempotents (C3) every coloured term is analytic on D'' and meromorphic at
tau = 0 (Lemmas 3.1(a), 8.2, 8.3, Proposition 4.3, F1, F3), so (13.1) also holds on D'' and for the Laurent series at t = 0; Omega_w is
t-independent. Notation: for a multiset v, h_v(z) := Exp[-Sum_k bk_k z^(2k-1) Sum_a v_a^(1-2k)], so that h_{u+v} = h_u h_v (union of multisets),
h_{-v} = h_v^-1 and h_v(z) h_v(-z) = 1; h_w is that of (C8).

**Lemma 13.3 (normalisation constants).** Let the weights be indeterminate and let r^pt_{0,+-}(z; Q) := Sum_m z^m wickZ_m(W_2(U_+-), W_3(U_+-), ...)
be the normalised expansion (C7) of e^{W/z} dU/U (W = spotU, rc_lib.m line 80) at the two roots of LeqU (line 115), labelled by U_+ ->
-alpha2/alpha0 and Q U_- -> alpha0/(alpha0 + alpha1) as Q -> 0. Its z-coefficients are analytic at Q = 0, with
    r^pt_{0,+}(z; 0) = h_{(alpha2, -alpha0-alpha2, alpha0)},      r^pt_{0,-}(z; 0) = h_{(-alpha0, -alpha1, alpha0+alpha1)}.
Hence h_{wt_+}/r^pt_{0,+}(z; 0) = h_{wt_-}/r^pt_{0,-}(z; 0) = h_nw with nw := (alpha1, alpha0 + alpha2, -eps1, -eps2), i.e. ffcLimit = ffcClosed
(rc_lib.m lines 91-109) to all orders in z (previously verified only to z^6, comment at line 105).
Proof. At U_+, theta^k log(1 - Q U) = O(Q), so the W_k(U_+) tend to those of f = a log U + b log(1 - U), (a, b) = (-alpha2, alpha0 + alpha2), at its
critical point a/(a + b). At U_-, V := Q U_- is analytic at Q = 0, theta_U = theta_V, theta log(1 - V/Q) = V/(V - Q) -> 1 and theta^k log(1 - V/Q)
-> 0 for k >= 2, so the W_k(U_-) tend to those of a log V + b log(1 - V), (a, b) = (alpha0, alpha1). In both cases W_2 -> -a (a + b)/b != 0, so by
(C7) and (D1) each z-coefficient tends to that of Lnorm[f, 1] in y = log U at the critical point. This limit is h_{(-a, -b, a+b)}: both sides have
z-coefficients rational in (a, b) ((D2); h by definition), so let a, b > 0. Then f(y) = a y + b log(1 - e^y) is strictly concave on (-oo, 0) (f'' =
-b e^y/(1 - e^y)^2), tends to -oo at both ends, and Int e^{f/z0} dy = B(a/z0, b/z0 + 1) < oo (DLMF 5.12.1), so Lemma 3.8 (am = 1) gives Int_0^1
U^(a/z) (1 - U)^(b/z) dU/U = B(a/z, b/z + 1) = e^{f*/z} (2 pi z b/(a (a + b)))^(1/2) [Sum_{m<M} ([z^m] Lnorm) z^m + O(z^M)] with f* := a log a + b
log b - (a + b) log(a + b). Stirling's series (DLMF 5.11.1, s.5.11(ii)) for log B(a/z, b/z), as in Proposition 3.15(3), and B(a/z, b/z + 1) = B(a/z,
b/z) b/(a + b), whose factor b/(a + b) does not depend on z, give the same prefactor times h_{(-a, -b, a+b)} + O(z^M), and asymptotic expansions are
unique. Substituting the two pairs (a, b) gives the limits, and h_{wt_+} h_{(-alpha2, alpha0+alpha2, -alpha0)} = h_{wt_-} h_{(alpha0, alpha1,
-alpha0-alpha1)} = h_nw because pairs v, -v cancel. QED (Check: (a, b) = (-1, c), the T_F phase at Q = 1 and U = -1/x, gives h_{(1, -c, x)} =
h_w^-1, Proposition 8.6(a).)

Proof of Proposition 13.2. Kw := Q(alpha0, alpha1, alpha2, eps1, eps2).
(a) The CohFT and its reconstruction. A T-fixed stable map to Z has image in the compact T-invariant curves, i.e. in the zero section, so the
T-fixed loci of M_{g,n}(Z, d) are those of the proper M_{g,n}(P^1, d), and R pi_* f^* Nb, Nb := O(-1)_{alpha1} + O(-1)_{alpha2} + O_{-eps1} +
O_{-eps2}, is moving there. Hence, by [GP99] on M_{g,n}(P^1, d), the T-equivariant classes Omega^T_{g,n} of Z (defined by localisation, [BS]
Notation 2.1, Def. 2.2) are those of the e_T(R pi_* f^* Nb)^-1-twisted T-equivariant theory of P^1 in the sense of [CG07] s.1, and they form a CohFT
with flat unit over Kw[[Q]] for eta(a, b) = Int_{P^1} a b/e_T(Nb) = Int_Z a b: the GW classes of P^1 do, and the twisting class is multiplicative
along the forgetful and gluing maps ([CG07] App. 1, (iv)-(vi): R pi_* f^* Nb is pulled back along the forgetful maps, and restricts to a boundary
divisor as the sum of the two sides minus ev^* Nb), which replaces the Poincare pairing by eta. The genus-0 potential is PrepCl + Li_3(e^{t1})/(eps1
eps2) (rc_lib.m lines 18-21; [BS] App. A): the classical part is localisation (proof of Lemma 4.1(a)); for d > 0, H^0(C, f^*O(-1)) = 0 gives
M_{0,0}(Z, d) = M_{0,0}(P^1, d) x C^2 with obstruction bundle R^1 pi_* f^*(O(-1) + O(-1)), whose rank is the virtual dimension, so, localising on
C^2, GW_{0,d} = C(0, d)/(eps1 eps2) with C(0, d) = 1/d^3 ([FP00] p. 3, the Aspinwall-Morrison formula), and <H, H, H>_{0,d} = d^3 GW_{0,d} by the
divisor equation. At Q = 0 the quantum product is the product of H_T(Z) (x) Kw = Kw e_+ + Kw e_- (e_p|_q = delta_pq), so the CohFT is semisimple. By
Teleman's classification ([T12] s.1.9, Theorem (2); s.8.1, Proposition 8.1(ii); s.8.2; applied over the field Kw((Q))), Omega^T = R.omega for the
unit-preserving action of [PPZ] Defs. 2.2, 2.13 on omega(e_i^{(x)n}) := Delta_i^(g-1), where Psi R e^{U/z} is a unitary fundamental solution of the
QDE ([G01] s.1.3, (3)-(5)); two such R differ by right multiplication by Exp[diag Sum_k a_k z^(2k-1)] with constant a_k ([G01] s.1.3, Proposition
(c); its proof is formal and holds over Kw((Q)) with d/dt = Q d/dQ). Integrated over M_g, R.omega is Givental's graph sum ([G01] (1)): vertices
Delta_i^(g_v - 1) Int psi^a, dilaton leaves T^i_d (d >= 2; T^i_0 = T^i_1 = 0, [G01] (7)), edges Delta_i^(1/2) Delta_j^(1/2) V^ij_kl ([G01] (1),
(6)). With the rows v_i := Delta_i^(1/2) e^{-u_i/z} S^i = e_i^* o R^-1(-z) (S^i the i-th solution, in the flat basis {1, H}) these are the factors
of graphValueC: the vertex weights (C6), Dil_i(d) = T^i_d (C10) and E_ij(k, l) = Delta_i^(1/2) Delta_j^(1/2) V^ij_kl (C11); and by symplecticity
Lemma 1.1 reads E_ij(z, w) = [delta_ij Delta_i - <v_i(-z), v_j(-w)>]/(z + w), the edge of [PPZ] Def. 2.2. So F^T_g := Int_{M_g} Omega^T_{g,0} is the
point-frame sum for the rows of R.
(b) The rows. Define R_gen by its rows v_i = e_i^* o R_gen^-1(-z) := (r_{0,i}, r_{1,i}), i = +-, with r_{0,i} := h_nw r^pt_{0,i} (Lemma 13.3) and
r_{1,i} := (L_i + z D_i) r_{0,i} ((H), L_i = LofU(U_i), rc_lib.m line 81): the construction of rmatrixAtClosed (line 138) with the weights kept
indeterminate. (i) Frobenius data. From PrepCl + Li_3(e^{t1})/(eps1 eps2), H * H = (L_+ + L_-) H - L_+ L_- 1 with L_+ + L_- = -(alpha0 + (alpha1 +
alpha2) Q)/(1 - Q) and L_+ L_- = -alpha1 alpha2 Q/(1 - Q). At a root U_i of LeqU, L_i = -alpha1 Q U_i/(1 - Q U_i) = du_i/dt for u_i := W(U_i), and
U_i/(1 - U_i) = (L_i - alpha2)/(alpha0 + alpha2); with e_i = (H - L_j)/(L_i - L_j) and sig := alpha0 + alpha1 + alpha2 this gives the generic form
of (V1),
    Delta_i = 1/eta(e_i, 1) = alpha1 alpha2 (alpha0 + alpha1)(alpha0 + alpha2) eps1 eps2 (L_i - L_j)/(alpha1 alpha2 - sig L_j)
            = eps1 eps2 (alpha0 alpha1 alpha2 + 2 alpha1 alpha2 L_i - sig L_i^2) = alpha1 (alpha0 + alpha2) eps1 eps2 W_2(U_i),          (13.2)
which is (4.1) on T_F. ([BS] App. A states s^2 d_s^2 W = alpha1 (alpha0 + alpha2) eps1 eps2 Delta_+- at the critical points, with the constant on
the wrong side; its own formula for Delta_+- agrees with (13.2).) (ii) Row 2. By (13.2) and du_i/dt = L_i, L_i + z D_i = z (d/dt + u_i'/z -
W_2(U_i)'/(2 W_2(U_i))), which is z times the operator of Proposition 3.6 (P2) for the phase Wy(y) := W(e^y) at y_i = log U_i, in which the code's
row is the normalised expansion with amplitude 1 (Proposition 3.10, step (1)). With d/dt = Q d/dQ at fixed U and PF := z^2 (Q - 1) (d/dt)^2 - z
(alpha0 + (alpha1 + alpha2) Q) d/dt + alpha1 alpha2 Q = (Q - 1)[(z d/dt)^2 - (L_+ + L_-) z d/dt + L_+ L_-], which is (Q - 1) times the operator of
[BS] (A.2), Euler's integration by parts is the identity in U and Q (multiply out)
    PF e^{W/z} = theta(km e^{W/z}),      km := -z alpha1 Q (1 - U)/(1 - Q U).   (13.3)
So (P2) twice, (D3) and (P3) give [(L_i + z D_i)^2 - (L_+ + L_-)(L_i + z D_i) + L_+ L_-] r_{0,i} = (Q - 1)^-1 h_nw Lnorm[Wy, theta km + km W_1/z;
y_i] = 0 at every Q with distinct roots not in {0, 1, 1/Q}: this is Row 2 of the QDE, z D_i v_i = (A - L_i) v_i with A^T the matrix of H * (Theorem
Q(i); the argument of Theorem 8.4, Steps C-D, with Gauss's operator replaced by PF). (iii) Q = 0. U_+ and V := Q U_- are analytic at Q = 0, every
W_k(U_i) is a polynomial in (1 - U_i)^-1 and (1 - Q U_i)^-1 (as in Lemma 8.2(i)), analytic at Q = 0, W_2(U_i) has the non-zero limit -a (a + b)/b of
Lemma 13.3, L_+ -> 0, L_- -> -alpha0, and Delta_+- -> e(wt_+-) by (13.2). So v_i lies in Kw[[Q]][[z]]^2, and by Lemma 13.3 and Q d/dQ = 0 at Q = 0,
v_i(z)|_{Q=0} = h_{wt_i}(z) (1, L_i(0)): R_gen|_{Q=0} = diag(h_{wt_+}, h_{wt_-}). (iv) Symplecticity. H * is eta-self-adjoint, eta^-1 is constant
and L_+ - L_- is a unit of Kw[[Q]], so the proofs of Corollary 8.5 and Proposition 4.8(a) apply: Y := <v_+(z), v_-(-z)> satisfies z Y' = [(L_- -
L_+) + (z/2)(Delta_+'/Delta_+ + Delta_-'/Delta_-)] Y, hence Y = 0, and <v_i(z), v_i(-z)>/Delta_i is Q-independent, hence equal to its value
h_{wt_i}(z) h_{wt_i}(-z) = 1 at Q = 0. So Psi R_gen e^{U/z} is a unitary fundamental solution of the QDE, and R = R_gen Exp[diag Sum_k a_k z^(2k-1)]
by [G01] Proposition (c).
(c) The constants. For a multiset v put Hodge_v := prod_a Lambda^vee(v_a)/v_a^g. Since ch_2k(E) = 0 by Mumford's relation c(E) c(E^vee) = 1 ([M83];
[FP00] s.1.2, s.2.1), Hodge_v = Exp[-Sum_k (2k-2)! ch_{2k-1}(E) Sum_a v_a^(1-2k)]; by Mumford's formula ([M83]; [FP00] s.1.1) a scalar Exp[Sum_k b_k
z^(2k-1)] acts on 1 by Exp[Sum_k b_k ((2k)!/B_2k) ch_{2k-1}(E)] ([T12] s.8.2; [G01] s.2.3, Lemma), so Hodge_v = h_v.1. At Q = 0 only constant maps
contribute: by [GP99] on M_{g,n} x {p} (N^vir = T_p Z - E^vee (x) T_p Z), Omega^T_{g,n}(e_p, ..., e_p)|_{Q=0} = e(wt_p)^(g-1) Hodge_{wt_p}, mixed
insertions vanish, and Delta_p|_{Q=0} = e(wt_p). The CohFT R|_{Q=0}.omega has components Delta_p^(g-1) (h_{wt_p} s_p).1 with s_p := Exp[Sum_k a_k^p
z^(2k-1)], and comparing with Hodge_{wt_p} = h_{wt_p}.1 gives s_p.1 = 1 ([PPZ] Prop. 2.14). Now ch_{2k-1}(E) != 0 on M_k: for k >= 2, (2k-1)!
ch_{2k-1}(E) = (-1)^(k-1) lambda_{k-1} lambda_k by c(E) c(E^vee) = 1, and lambda_{k-1}^3 = 2 lambda_k lambda_{k-1} lambda_{k-2} has non-zero
integral ([FP00] Thm. 4); for k = 1, Int_{M_{1,1}} lambda_1 = 1/24. Induction on k gives a^p = 0, so R = R_gen, and F^T_g in Kw[[Q]] is the
point-frame sum for R_gen. (For compact Hamiltonian targets with isolated fixed points and 1-dimensional orbits this normalisation of R is [G01]
s.2.2, Theorem 1 and the Remark after it, proved by localisation; Z has such orbits but is not compact. [Z16] Thms. 5.4, 5.7, 6.2, a preprint,
states it for GKM orbifolds, via [T12] and orbifold quantum Riemann-Roch.)
(d) Specialisation to T_F. The factors of (b) are elements of QA over Kw(Q) whose coordinates are regular along Lie(T_F) at generic Q: the two roots
of LeqU are distinct (discriminant ((x - 1)^2 + 4x/Q)/x^2) and avoid {0, 1, 1/Q}, W_2 and Delta are finite and non-zero at both, eta^-1 is
polynomial ([BS] (A.1)), and h_nw -> h_{(eps1, eps2, -eps1, -eps2)} = 1. So the colour-summed total, cconj-invariant hence in Kw(Q), restricts to
the point-frame sum for R_Z := the rows of rmatrixAt0 (ffc = 1; rc_lib.m line 139). On the GW side, by (a), GW^T_{g,d} = Int_{[M_{g,0}(P^1, d)]^vir}
e_T(R pi_* f^* Nb)^-1 over a proper space on which T acts through alpha0, and the fibre characters of Nb are alpha1, alpha2, -eps1, -eps2 up to
multiples of alpha0; expanding e_T^-1 in alpha0 puts GW^T_{g,d} in Q[alpha1, alpha2, eps1, eps2][1/(alpha1 alpha2 eps1 eps2)][[alpha0]], with image
GW^{T_F}_{g,d} under alpha0 -> 0, alpha_i -> eps_i. Writing the total as num/den with num, den polynomials in the weights and Q, den|_{T_F} != 0,
and applying this ring map to den F^T_g = num gives F_g = num|/den| in Q(x)(Q): the point-frame sum for R_Z.
(e) Frame B. By (c), Omega_w = Hodge_w = h_w.1 in the conventions of [PPZ] Def. 2.13, so Omega_w (x) omega = h_w.omega (h_w is scalar) and, by the
group property [PPZ] Prop. 2.14, R'.(Omega_w (x) omega) = (R' h_w).omega. Integrated over M_g, the left side is the Hodge-vertex sum (13.1) for R'
(vertex classes Delta_i^(g_v - 1) Omega_w; translation leaves Dil_i and edges E_ij from the rows of R'), and the right side is the point-frame sum
for R' h_w, which is F_g for R' h_w = R_Z by (d). So R' = R_Z h_w^-1. Multiplying the code's rows by h_w(z)^s is R -> R h_w^s (v_i = e_i^* o
R^-1(-z) and h_w(z) h_w(-z) = 1), so this is (C8) with s = -1, and (13.1) follows. QED

Remarks on Proposition 13.2. (i) s = -1 is derived in (e), not fitted. It agrees with the empirical choice (hvc_lib.m lines 5-6; hvc.log line 31: in
genus 2, s = -1 reproduces F_2 in modes A and B, s = +1 in neither). Omega_w = h_w.1 is checked by ppzScalar (hvc.log lines 6-20: genus 2 at
symbolic w, genus 3 at two numerical w); ppzScalar takes R^-1(psi) = h_w(psi)^-1 (hvc_lib.m lines 101-105), i.e. code rows h_w(z), and "rows h_w^-1"
in hvc.log line 6 refers to R^-1(psi). (13.1) holds for every w, with w-dependent frame-B factors; w = (c, -x, -1) is the choice for which F2 holds
(s.2, remark (1)). (ii) h_w is scalar (the same at both idempotents), t-independent and h_w(z) h_w(-z) = 1, so R' is again a unitary solution of the
QDE (the ambiguity of [G01] Proposition (c); on T_F, Theorem Q proves both properties of the frame-B rows directly), h_w.omega = Omega_w (x) omega,
and all t-dependence of (13.1) sits in Delta_i, Dil_i, E_ij. (iii) h_nw is c_+ of [BS] App. A (where Theta(z/v) = h_{(v)}^-1); at alpha0 = 0, nw
lists the fibre weights of O(-1)^2 + O^2, and on T_F h_nw = 1 because nw = (v, -v) (equivalently Lambda^vee(e) Lambda^vee(-e) = (-e^2)^g; README
lines 18-22), so rmatrixAt0 needs no constants. (iv) The passage to T_F is the algebraic step (d), not an analytic continuation in Q from large
radius: on T_F the point Q = 0 is not semisimple (both roots of LeqU -> oo, both L_i -> 0), and the stable-graph terms have poles there that cancel
only in the total (README lines 44-47). For x > 0 the discriminant (c^2 - 4 x tau)/x^2 has no zero for Q in (0, 1], but the endpoint Q = 0, where
the constants are fixed at generic weights, is degenerate on T_F. (13.1) is an identity in Q(x)(Q), valid for every x != 0, -1, so x > 0 is not
needed; the discriminant vanishes only at tau = c^2/(4x), outside D'. (v) The translation of [PPZ] Def. 2.13 for R', z(1 - R'^-1(z) 1), has no z^0
and z^1 terms; its colour-i part is Sum_{d>=2} Dil_i(d) psi^d (C10), and the omitted -psi is the dilaton shift.

Evidence (exact, Q and x symbolic). Genus 2 and 3: the frame-B sums equal F_{g,0} of ref/resolvedconifold_results.m (hvc.log lines 107-109 and
1917-1919; hvc_results_g2_Bm.m, hvc_results_g3_Bm.m: totalOK, gapOK). That reference is the generic-torus point-frame Givental sum of rc_lib.m with
the ffcLimit constants (equal to ffcClosed by Lemma 13.3), reconstructed from 97 integer weight samples (resolvedconifold_genus3_checks.log line 1)
and specialised, so this comparison tests (d), (e) and the constants, not (a) or (b). Steps (a)-(c) are tested independently against [BS]: the
degree-1 localisation formula of [BS] s.3.1 in genus 1, 2 (checks log lines 6-7), and [BS] Prop. 3.1 (F_g - F_g(0) = P_g Li_{3-2g}(Q)/(eps1 eps2) on
T_F, all genera, proved by degeneration) in genus 2, 3 (lines 12-13). Genus 4 (27172 decorated graphs; hvc.log lines 1961-1965 and 1970;
stanage_genus4/hvc_results_g4_Bm.m: totalOK, gapOK): total - P_4 Li_{-5}(Q)/(eps1 eps2) is Q-independent and equals -(x^6 + 1)/87091200. This is
F_4(0): for g >= 2, F_g(0) = -(eps1^(2g-2) + eps2^(2g-2)) |B_2g B_{2g-2}|/(4g (2g-2) (2g-2)!) (F0closed, hvc_lib.m line 178, labelled there a
prediction). Indeed [GP99] on M_g x P^1 (fixed part e(E^vee (x) T_P1) = (-1)^g (lambda_g - 2 lambda_{g-1} h), h the point class; moving part Nb -
E^vee (x) Nb) and Lambda^vee(u) Lambda^vee(-u) = (-1)^g u^(2g) leave, in degree 3g - 3, F_g(0) = (-1)^(g+1) (eps1^(2g-2) + eps2^(2g-2)) Int_{M_g}
lambda_g (2g-3)! ch_{2g-3}(E); c(E) c(E^vee) = 1 gives lambda_g (2g-3)! ch_{2g-3}(E) = (-1)^g lambda_g lambda_{g-1} lambda_{g-2} and lambda_{g-1}^2
= 2 lambda_g lambda_{g-2}, and [FP00] Thm. 4 gives the value. So the genus-4 comparison is with a proven F_4.

References. [BS] A. Brini, Y. Schuler, Refined Gromov-Witten invariants, arXiv:2410.00118 (v1, 2024). [CG07] T. Coates, A. Givental, Quantum
Riemann-Roch, Lefschetz and Serre, Ann. of Math. (2) 165 (2007) 15-53 (math/0110142). [FP00] C. Faber, R. Pandharipande, Hodge integrals and
Gromov-Witten theory, Invent. Math. 139 (2000) 173-199 (math/9810173). [G01] A. Givental, Semisimple Frobenius structures at higher genus, Internat.
Math. Res. Notices 2001, no. 23, 1265-1286 (math/0008067). [GP99] T. Graber, R. Pandharipande, Localization of virtual classes, Invent. Math. 135
(1999) 487-518. [M83] D. Mumford, Towards an enumerative geometry of the moduli space of curves, in Arithmetic and Geometry II, Progr. Math. 36,
Birkhauser 1983, 271-328. [PPZ] R. Pandharipande, A. Pixton, D. Zvonkine, Relations on M_{g,n} via 3-spin structures, J. Amer. Math. Soc. 28 (2015)
279-309 (1303.1043). [T12] C. Teleman, The structure of 2D semi-simple field theories, Invent. Math. 188 (2012) 525-588 (0712.0160). [Z16] Z. Zong,
Equivariant Gromov-Witten theory of GKM orbifolds, arXiv:1604.07270 (preprint; not used in the proof). DLMF: NIST Digital Library of Mathematical
Functions, 5.11.1, s.5.11(ii), 5.12.1.

## 14. Remaining gaps, errata, open points

Gaps (stated as precisely as possible; none of them is inside the proof of Theorem LC, except where (G2) says so).
 (G1) The all-genus regularity theorem remains conditional on M_irr(h) for h >= 5 (Lemma A). Nothing in this text bears on it.
 (G2) Analytic inputs. Proposition 3.15 (Abeta = h_w c/(c + z)) is the one non-formal step in the proof of Theorem LC: it fixes the frame-B constant
      c/(c + z), on which F2, L1 at t^0, sigma^(-1) = (0, c), K_1 = 1 and hence L2, L3 and F3 for E_11 depend. It is proved from Laplace's method
      (Lemma 3.8, proved in full) and Stirling's series with an explicit remainder bound for real positive argument, which is cited (DLMF 5.11.1 and
      s.5.11(ii); Whittaker-Watson s.12.33), not re-proved; the passage from x > 0 to all x is by rationality. No purely formal proof of Proposition
      3.15 is given. (The regularity at z = 0 of every tau^k-coefficient does have a purely formal proof, Proposition 3.19.) Proposition 8.6(a)
      (r^pt_{0,2}(z; 0) = h_w^-1) rests on the same two inputs; it is used only for K_2 = 1 and for the t = 0 values at the regular idempotent
      (E_22, Dil_2, the E_12 legs at t = 0 and their top coefficients), not for Theorem LC or Corollary R. Lemma 13.3 (s.13.3) rests on
      the same two inputs; it serves input (1) of Theorem R (Proposition 13.2), not Theorem LC.
 (G3) Closed as a write-up item. The premise that the frame-B Hodge-vertex graph sum of rc_lib.m/hvc_lib.m computes the conifold free energies
      (input (1) of s.13.3) is now Proposition 13.2, proved from: the CohFT of Z at indeterminate weights (virtual localisation, as a twisted theory
      of P^1) and Teleman's classification; the identification of its R-matrix with the code's Laplace rows times ffcClosed (generic V1 (13.2) and
      Euler's identity (13.3), proved there; Lemma 13.3 proves ffcClosed = ffcLimit to all orders, previously verified only to z^6) and the fixing
      of its constants at Q = 0 by the constant maps; the restriction to T_F as an identity in Q(x)(Q), where the constants become h_{(eps1, eps2,
      -eps1, -eps2)} = 1 (no analytic continuation from large radius, which is not a semisimple point on T_F); and the frame change by Mumford's
      formula and the group property of the Givental action, which forces s = -1. Cited, not re-proved: Teleman (Invent. Math. 188, 2012: s.1.9
      Theorem (2), s.8.1-8.2); Givental (IMRN 2001: s.1.3 Proposition (c), formulas (1), (6), (7), s.2.3 Lemma); Pandharipande-Pixton-Zvonkine (JAMS
      28, 2015: Defs. 2.2, 2.13, Prop. 2.14); Graber-Pandharipande (virtual localisation); Coates-Givental (Ann. of Math. 165, 2007: twisted
      theories, App. 1 (iv)-(vi)); Mumford's formula and relation c(E) c(E^vee) = 1 and Faber-Pandharipande (Invent. Math. 139, 2000: s.1.1-s.2.1,
      Thm. 4, C(0, d) = 1/d^3); Brini-Schuler Def. 2.2. Lemma 13.3 uses the analytic input of (G2) (Lemma 3.8, Stirling's series). Givental's
      localisation theorem (compact targets) and Zong's preprint arXiv:1604.07270 are quoted only for comparison. Inputs (2)-(5) of s.13.3 are
      proved in the cited notes and are not re-examined here ((5): Brini-Schuler Prop. 3.1, proved by degeneration).
 (G4) Theorem LC is a statement about the formal series defined by the code's formulas (s.1.2 (C7)-(C11)); the code computes exactly their truncation
      at z^ZMAX, and every Mathematica confirmation is limited to the stored z^9 factors (t-orders as listed in s.12.2). The proof itself has no such
      limitation.
 (G5) Not proved and not needed: a closed form in b of the lower E_12 leg coefficients (lam_j with j <= k at order t^k; they are determined
      recursively by Remark 8.9); the rank-one pattern of the lam_1 coefficients of the E_12 legs at k = 1, 2 (leaf_structure.md (ii)) is not
      structural and fails at k = 3. The statement R~'_11 - 1 = O(t^2) in F2 of conifold_regularity_note.md is not re-derived. The sharpness
      statements hold only for x outside a finite set depending on k.
 (G6) The Mathematica checks M1-M16 of s.12.3 have all been run and passed (s.12.2(e); M1 and M16 at x = 2, 3/5 with factors rebuilt to z^12).
 (G7) The independent referee scripts of s.12.2(c) and the assembly spot checks of s.12.2(d) are archived in leaf_proof/referee_scripts/.

Errata to other files (found while proving; applied 2026-10-31 after the Mathematica confirmation of s.12.2(e)):
 (E1) conifold_regularity_note.md, F3: "E_11(z, w) = E_22(z, w) = -(x^3/2) phi(z) phi(w)" replaced by "E_11(z, w) = -(x^3/2) phi(z) phi(w), and
      E_22(z, w) = (x^3/c) [1 - h_w(z)^2 h_w(w)^2]/(z + w), which is not of rank one" (Proposition 8.6(c), Remark 8.8); and "in the point frame and in
      frame A the same-idempotent edge factors E_11, E_22 have double poles" replaced by the statement that only E_11 has one (Remark 8.8(3)).
      Theorems A and B of the note are unaffected; for the Reduction Lemma see (E5).
 (E2) hvc_vertex_identity.wls line 14, and hvc_asymptotics.wls, hvc_asymptotics_g4.wls line 24: EdN[2,2] is now evaluated at U1, not at U2
      (conjugate encoding, Lemma 8.1; Remark 8.8(2), (3)).
 (E3) outline.md, Step 9: added "E_21: E_21(b, a) = E_12(a, b) by the orthogonality of Step 6, which needs the QDE at the regular idempotent (Theorem 8.4,
      Corollary 8.5, Proposition 9.4); already the regularity of E_21 needs the t^-1 part of it (Proposition 8.6(b))." Outline Step 6's argument for
      the orthogonality ("a formal series P(z) with e^{(u_1-u_2)/z} P t-independent vanishes") is made rigorous by Corollary 8.5.
 (E5) conifold_regularity_note.md s.5a (Reduction Lemma): the factorisation of an E_11 block between two singular vertices of genus >= 2 is not
      possible with coefficients regular at t = 0 (Remark 13.1); a note added there says so and points to Theorem R, which is the valid
      reduction. The note's s.7 ("by the Reduction Lemma") now cites Theorem R.
 (E6) proof_status_mixed_scale.md: s.2 Lemma C marked proved, with k >= 0 and the E_21 clause (s.13.2 (d)); the s.3 sentence on the Reduction
      Lemma qualified as in Remark 13.1; addendum (2026-10-31) with the status after this text.
 (E4) Notation in the section files (resolved here; see Appendix A): kappa (S1's phase constant kappa_ph vs the outline's kappa(z) vs S4's first
      zeta), omega (S2's omega(z, tau) vs zeta/t in S3 s.1.6 and S4_checks.pl), g (S1's amplitude vs g(psi) = G(-psi)), sigma_k (S3: component; S4:
      t-order), K (S2: C((z)); S3: Q(x)), Y (S3: t sigma; S2/S4: the orthogonality pairing), A (S1's constant A(z) vs the matrix A), delta_k (leaf
      decomposition vs delta), P (projector, prefactor, transport coefficient, polynomials, P_k), q_0 (S2.10's amplitude vs q_0 = c - c^2 lam_1).

Open points that are not gaps of the proof: the closed forms of section 11 are complete for the singular data (all leaf and loop coefficients
as explicit functions of t); for the regular data only the t = 0 values and the recursive description of Remark 8.9 are available.


## Appendix A. Notation (with the names used in the section files and other notes)

The column "where" gives the subsection (s.N.M) or the statement where a symbol is defined.

    symbol           meaning                                                                  where     section files / other notes
    x, c             eps1 = x, eps2 = 1, c = 1 + x                                              1.1
    Q, t, tau        Q = e^t, tau = 1 - e^-t                                                    1.1
    f', f_tau, Dz    d/dt, d/dtau, Dz = z d/dt                                                  1.1       Dz = D (outline V7, S2.6)
    W, theta, W_k    superpotential spotU, U d/dU, theta^k W                                    1.2 (C1)  W2 = W_2
    U_1, U_2         singular and regular critical points                                       1.3       U1, U2
    delta            c (1 - 4 x tau/c^2)^(1/2), delta(0) = c                                    1.3       (not the leaf decomposition)
    s*, s~           critical points of phi_tau, s*(0) = 1/c, s~ = 1/(x tau s*)                 1.3       s~ = s** (S1), s_2 (S4)
    D, D', D''       discs |tau| < |c|^2/(4|x|); intersected with |tau| < 1; minus 0              1.3
    L_i, l           eigenvalues of H*; L_1 = c/t + l                                           1.2 (C5), 4.2
    Delta_i, rho     metric 1/eta(e_i, e_i); rho = (Delta_1 t^2/(-x c^3))^(-1/2)                1.2 (C6), (C13)
    h_w, bk_k, St_k  frame-B twist, B_2k/(2k(2k-1)), 1 + x^(1-2k) - c^(1-2k)                    1.2 (C8)  beta_k, S_k (S4 10.6)
    r_{a,i}, v_i     frame-B rows, (r_{0,i}, r_{1,i}); point frame r^pt                         1.2       ratCN[a,i]
    D_i              d/dt - Delta_i'/(2 Delta_i)                                                1.2 (C9)
    <a, b>           x^2 (a_0 b_1 + a_1 b_0) - c x a_1 b_1                                      1.2 (C4)
    q_t, q_0, q^(k)  dilaton series, its t = 0 value, [t^k](rho q_t)                            1.2, 1.5  q^(k) = delta_k (leaf_structure.md)
    f_{j,b}          legs [w^b] E_1j(psi, w)                                                    1.2 (C12)
    lam_j, Lam_m     1/(c - j psi), span(lam_1, ..., lam_m)                                     1.2, 1.4  Lam_m = V_m (S4)
    phi_tau, dphi    Euler phase, phi_tau(s*) - phi_0(1/c)                                      1.5
    alg              (1 - tau s*) (phi_tau''(s*)/phi_0''(1/c))^(1/2)                             1.5
    F, F_n           2F1(1/z + 1, 1/z + 1; c/z + 2; tau) and its tau-coefficients                1.5
    alpha, beta, gamma Gauss parameters 1/z + 1, 1/z + 1, c/z + 2                               1.5
    Ephi, X_m        Exp[-dphi/z], its tau-coefficients                                         1.5       E_m (S3 first version)
    H_n, Psi_n       [tau^n](F Ephi), [tau^n](F Ephi alg)                                       1.5
    G, g, Dg         rho r_{0,1}, G(-psi), d/dt(t g)                                            1.5
    a_j^(k), A_j     lam_j-coefficient of q^(k), Sum_k a_j^(k) t^k                              1.5       b_j^(k), b_j(t) (S3 first version)
    g_j^(k)          lam_j-coefficient of g_k, -(j/c) a_j^(k)                                   1.5
    Ahat_j, Atil_j   t A_j, t^(1-j) A_j                                                         1.5, 11.2 B_j (S3, S4), A~_j (S4)
    e^(k)_il, e_il   loop coefficients                                                          1.5       e_il^(k) (S4)
    P_k, C_j, zeta   prod (k+1-ic)^2, P_{j-1}/(j!(j-1)! c^(2j-4)), tau e^(dphi/c)               1.5       zeta = kappa(t) (S4 first version); zeta/t = omega (S3 1.6, S4_checks.pl)
    A, chi_i, pi     QDE matrix, (1, L_i), <., chi_2>                                           1.5, 4.7
    mu, omega        c/tau - 1, (c + z)/tau - (1 + z)                                           1.5
    Theta            (L_2 - L_1)(c x L_2 - x^2)^2/Delta_2 = -x^3 tau s*^2                        1.5, 7.2
    sigma, sigma^(k) rho v_1(-psi) (components sigma_0, sigma_1), [t^k] sigma                   1.5, 6    S4 writes sigma_k for [t^k] sigma
    S_j^(k), S_j(t)  coefficient vectors of sigma                                               6.2, 9.3
    Abeta(z)         Lnorm[phi_0, 1; 1/c] = h_w c/(c + z)                                        3.5, 3.15 A(z) (S1)
    N, L, NG         Lhat[phi_tau, (1 - tau s)^-1; s*], Lnorm[...], Exp[dphi/z] N/alg            3.5, 3.6  M (S1 Prop. 4.1)
    wr_k             tau^k W_k(U_1)                                                             3.3       w_k (S1)
    Pol_n, pc        polynomials of (3.7) and their coefficients                                3.3       P_n, a_{n,i} (S1)
    am, am0, kE, Dam amplitude, (1 - tau s)^-1, s(1-s)(1-tau s)^-2, d_tau + (d_tau f)/z            3.4, 3.6  g, k, delta_tau (S1); q_0 (S2.10)
    hc               coordinate change in (P1)                                                  3.4       h (S1)
    P_tr, R_tr, T_tr transport coefficients (P_tr = delta)                                      3.10      P, R, T (S1 7.4)
    Sph, Phi         eikonal phase (= dphi), transport solution                                 3.10      S (S1 7.4)
    kappa_ph         -e^{x log(-1)/z}, phase constant of (3.2)                                   3.1       kappa (S1 first version)
    kappa_0, kappa_1 kappa_0^2 = -x c^3, c/(kappa_0 (c + z))                                     4.10      kappa(z) of outline Step 5 = kappa_1
    nu_1, nu_2, vn   t rho/kappa_0; 1/(x U_2 f2r^(1/2)); t rho v_1                               4.6-4.10, 8.7  w (S2)
    PhiF             (F, (Dz + omega) F)                                                        4.9       Phi (S2 Remark S2.9.1)
    Pf, Mod_lambda   tau eps_mu; Rz eps_lambda                                                   4.5, 4.10 P, E_lambda (S2)
    K_1, K_2         symplectic constants                                                       4.8, 8.7
    Y, Ebar, Nbar    t sigma, t^2 rho^2 E_11, the numerator of (7.1)                            7.1
    Y12              <v_1(z), v_2(-z)>                                                          8.5       Y (S2, S4)
    f0, f2, f2r      phi_tau(s~), phi_tau''(s~), f2/tau^2                                       8.3
    L2, Dt, Dop      Lnorm at s~, d/dt + f0'/z - f2'/(2 f2), Dt/(1 - tau)                        8.4       Lr (S2), L (S4)
    N12              rho P_12(psi, w)                                                           9.2       N (S4)
    ph12, fe         P_12(-u, u), the t^-1 function of E_21                                     9.4       phi, f (S4)
    chi12            t = 0 leg of E_12                                                          11.3      chi (S4, note F3)
    Fb, sg           y - c log(1 + e^y), e^y/(1 + e^y)                                           8.6       F, ell (S4)
    Acal, <<.>>_a    weighted formal series in (v, z); Gaussian functional                      3.4
    Lnorm, Lhat, Lfull Gaussian-normalised, leading-coefficient-normalised and full formal expansions  Def. 3.4
    Sm, Sm^(2), mu_2 flat sections nu_i v_i eps_{L_i}; mu_2 = f0' = 1 - x s~                    4.10, 8.10 S, S^(2), lambda_2 (S2)
    Db, O_Db         disc around tau0 with closure in D''; analytic functions on it                 8.4       B, O_B (S4)
    Kx, Kz, Rz, Rz0  Q(x), C((z)), C((z))((tau)), C[[z]]((tau))                                  1.4       K (S2 = C((z)); S3 = Q(x)), R, R_0 (S2)
    Fpot_h, Xi_h     genus-h potential of Omega_w, its leaf derivatives                         13        F_h (proof_status)


## Appendix B. Concordance with the section files and the outline

Here "s.N.M" is a subsection and an unprefixed label is a statement or an equation; the left side of each entry uses the labels of the
section file, the right side those of this text.

Outline steps: Step 1 = s.3.1-3.3; Step 2(a) = Proposition 3.10; 2(b) = Proposition 3.12; 2(c) = Lemma 3.13, Theorem 3.14, Proposition 3.15;
2(d) = Corollary 3.16, Theorem 3.17, Corollary 3.18 (and Proposition 3.19); Step 3 = Propositions 4.3, 4.4; Step 4 = section 5; Step 5 = Lemma 4.5,
Propositions 4.6, 4.7; Step 6 = Proposition 4.8 and Corollary 8.5; Step 7 = section 6; Step 8 = section 7; Step 9 = section 9; Step 10 = section 8;
outline s.4 = section 11. Verified facts: V1 = Proposition 4.3(c); V2 = Lemma 4.2, Proposition 4.3(d); V3 = (4.3); V4 = section 5; V5 = Proposition
4.6, Theorem 8.4; V6 = Proposition 4.8(c), Corollaries 8.5, 8.7; V7 = Lemma 4.5; V8 = Theorem 7.2.

S1_thimble.md: s.1.1 = s.3.1; s.1.2, s.1.3 = Lemma 3.1; s.1.4 = Remark 3.20; s.1.5 and Lemma 1.1 = s.3.3 and Lemma 3.2; Lemma 2.1 = Lemma 3.3;
s.2.2 = Definition 3.4; Prop. 2.2 (P3) = Prop. 3.5; Prop. 2.3 (P2) = Prop. 3.6; Prop. 2.4 (P1) = Prop. 3.7; Lemma 2.5 = Lemma 3.8; Lemma 2.6 = Lemma 3.9;
Prop. 3.1 = Prop. 3.10; Remark 3.2 = Remark 3.11; (4.1) = (3.16); Prop. 4.1 = Prop. 3.12; Lemma 5.1 = Lemma 3.13; Thm 5.2 = Thm 3.14; Prop. 6.1 =
Prop. 3.15; Cor. 7.1 = Cor. 3.16; Thm 7.2 = Thm 3.17; Cor. 7.3 = Cor. 3.18; s.7.4 = s.3.10 (Prop. 3.19); s.7.5 = s.3.11; s.7.6 = Remark 4.10;
Theorem S1 = Lemmas 3.1, 3.2, Props 3.10, 3.12, Thm 3.14, Prop. 3.15, Thm 3.17.
S2_metric_qde.md: Lemma S2.1 = Lemma 4.1; Lemma S2.2 = Lemma 4.2; Prop. S2.3 = Prop. 4.3; Prop. S2.4 = Prop. 4.4; s.S2.5 = Remark 4.10; Lemma S2.6 =
Lemma 4.5; Prop. S2.7 = Prop. 4.6 (Remarks S2.7.2, S2.7.3 = (4.7) and the F-form after it; Remark S2.7.4 = the exponent part of Remark 4.10);
Prop. S2.8 = Prop. 4.7; Prop. S2.9 = Prop. 4.8; Remark S2.9.1 = Remark 4.9; Lemma S2.10.1 = Lemmas 8.2, 8.3; Lemmas S2.10.2, S2.10.3 = Steps A and C
of Theorem 8.4 (module form: Remark 8.10); Prop. S2.10.4 = Theorem 8.4; Prop. S2.10 = Cor. 8.5; Prop. S2.10.5 = Prop. 8.6(b) and Cor. 8.7;
s.S2.11 (I1)-(I6) = (4.7), Prop. 6.2, (7.3) with (4.8), (4.9) with Lemma 7.1, Cor. 8.5, Theorem 8.4 with Prop. 8.6(b) and Cor. 8.7.
S3_leaves_loops.md: (LAM1)-(LAM5), (PF), (T) = s.1.4; Lemma 1 = Lemma 5.1; Lemma 2 = Lemma 5.2; s.1.3 = s.5.3; Prop. 3 = Prop. 5.3; Theorem L1 =
Theorem 5.4; the Proposition of s.1.6 = Theorem 11.1; the Corollary (5) of s.1.6 = Corollary 11.3; (6) = (4.6); (7) = (6.1); Prop. 4 = Prop. 6.1;
(9) = (4.1); (10) = (4.8); Prop. 5 = Prop. 6.2; (11) = (4.9); (12) = (7.1); Lemma 6 = Lemma 7.1; Theorem 7 = Theorem 7.2; (15) = (7.4); Theorem L3 =
Theorem 7.3; (16), (17) = (7.5), (7.6); Cor. 8 = Cor. 7.4; s.4.1 = s.11.4.
S4_legs_closed.md: (1.1)-(1.2) = Lemma 1.1; Lemma 9.1 = Prop. 6.1 with (6.4); Lemma 9.2 = Lemma 9.1; Prop. 9.3 = Prop. 9.2; Prop. 9.4 = Prop. 9.3;
Prop. 9.5 = Prop. 9.4; Cor. 9.6 = Cor. 9.5; "Use in Theorem R" = s.9.6; s.10.1, s.10.2, s.10.4 = Lemma 8.2 (and Lemma 8.1); s.10.3 = Lemma 8.2(iii)
with Prop. 4.3(d), (e); Lemma 10.5 = Theorem 8.4 ((10.11)-(10.18) = (8.2)-(8.6)); Remark 10.5.1 = Remark 8.10 and the paragraph after Theorem 8.4;
Cor. 10.5' = Cor. 8.5; Prop. 10.6 = Prop. 8.6 ((10.19) = (8.7)); Remark 10.6.1 = Remark 8.8 and Lemma 8.1; Remark 10.6' = Remark 8.9; Remark 10.6'' =
Cor. 8.7; s.10.7 = section 10; CF.1 = Thm 11.1; CF.2 = Cor. 11.2; CF.3 = Cor. 11.3; CF.4 = Thm 11.4; CF.4' = Cor. 11.5; CF.5 = s.11.4; CF.6 = s.11.3.
