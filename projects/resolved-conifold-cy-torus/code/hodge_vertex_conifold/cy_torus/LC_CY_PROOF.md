# Theorem LC_CY: leaf structure of the singular column at generic Calabi-Yau weights (resolved conifold x C^2)

2026-11-03. Resolved conifold x C^2 at the conifold point Q = 1, torus weights (a0, a1, a2, e1, e2) with fixed-point tangent weights
wt_+ = (a0, a1, a2, -e1, -e2) and wt_- = (-a0, a0 + a1, a0 + a2, -e1, -e2), restricted to the CY 4-torus T_CY : E := a0 + a1 + a2 = e1 + e2,
Hodge-vertex frame B with class w = (E, -e1, -e2). This text extends Theorem LC of leaf_proof/leaf_structure_proof.md (2026-11-02; below
"the proof of record", cited as [PR] with its section and label numbers) from the fibrewise torus T_F (a0 = 0, a1 = e1 = x, a2 = e2 = 1,
c = 1 + x) to generic CY weights. It follows [PR] s.0-s.11 section by section. For every statement it says whether the statement of [PR]
is used verbatim (then it is cited, with the substituted constants quoted, and not re-proved), is a substitution (then the substituted
statement is written out and re-derived in full), or needs a new argument (then written in full). Finite-order computations are
evidence only and never proof steps; they are listed in s.13.3. No topological recursion is used. Four slips in intermediate displays (s.4.6 F-form of Row 2, Thm 11.4_CY factor (1-v), s.11.4 sign of g_1, Cor 12.1 leg exponent) were found while transcribing to LaTeX (hodge_vertex_conifold/note/) and corrected in place on 2026-11-03; none is used by a later step. Refereed twice (AI; no major item; the record is
s.13.1, whose fixes M1-M6 of 2026-11-03 are applied); no human has read this text.

Notation. A := a0 + a2 = E - a1, B := a0 + a1 = E - a2, u := (E - 2 a1)(E - 2 a2), tau := 1 - 1/Q = 1 - e^-t, t := log Q.
Non-degeneracy (ND): a1 a2 A B e1 e2 E != 0. The point of T_F in this notation: a0 = 0, A = a2 = 1, B = a1 = x, E = c, u = -(1 - x)^2.
The dictionary from [PR] to this text is: (x, 1, 1) -> (a1, a2, A) in the Euler phase, c -> E, x/tau -> a1 a2/tau, c/tau -> E/tau - a0,
-x c^3 -> -e1 e2 E^3, x^3 (the regular Delta_2(0) c) -> e1 e2 a1 a2 A B/E, class (c, -x, -1) -> (E, -e1, -e2), rows h_w^-1 r^pt -> h_w^-1 h_nw r^pt.

Contents. 0 Statement. 1 Conventions and factors at generic weights. 2 Plan. 3 The singular identity row from the Euler integral.
4 Metric, eigenvalues, QDE, symplectic constant. 5 Dilaton leaves (L1). 6 The H-row in leaf form. 7 The E_11 flow (F3, L3, L2 for j = 1).
8 The regular idempotent. 9 The legs (L2 for j = 2). 10 Assembly. 11 Closed forms (Theorem CF_CY). 12 The three-torus a0 = 0 by
rescaling. 13 Status, what was not re-derived, finite-order checks.

## 0. Statement

### 0.1 The objects
Exactly as [PR] s.0.1, with the factors of section 1 below: a coloured stable graph contributes
    (1/|Aut|) Sum_colourings Sum_decorations prod_v Delta_{c_v}^(g_v - 1) Int Hodge_w psi^a prod Dil_{c_v}(d) prod_e E_{c c'}(a_e, b_e),
built from the frame-B rows r_{0,i} (identity row) and r_{1,i} (H-row) at the two idempotents. Idempotent 1 is the critical point U_1 of the
superpotential that runs into the logarithmic singularity at t = 0 (singular idempotent), idempotent 2 stays away from it (regular). Put
    rho(t) := (Delta_1 t^2/(-e1 e2 E^3))^(-1/2) (branch with rho(0) = 1),   lam_j(psi) := 1/(E - j psi),   Lam_m := span(lam_1, ..., lam_m),
    q_t(psi) := Sum_{d>=1} Dil_1(d; t) psi^d,   f_{j,b}(psi) := [w^b] E_1j(psi, w),   q_0 := E - E^2 lam_1.

### 0.2 The main theorem
**Theorem LC_CY.** Let (a0, a1, a2, e1, e2) be complex numbers on T_CY with (ND), or indeterminates with E = e1 + e2 imposed. All factors are
the frame-B factors of section 1 (rows R' = R_Z h_w^-1, w = (E, -e1, -e2), R_Z the rows with the Q = 0 constant h_nw of [PR] Lemma 13.3),
at their true idempotents. Then:
 (F1) Delta_1 = -e1 e2 E^3/(t^2 rho^2), where rho = tau/(t alg) (alg as in s.1.5) is analytic at t = 0 with rho(0) = 1; so
      Delta_1 = -(e1 e2 E^3/t^2)(1 - u t/E^2 + O(t^2)). Delta_2 = e1 e2 B delta s* L_2 is analytic at t = 0 with Delta_2(0) = e1 e2 a1 a2 A B/E != 0
      and Delta_2'(0) = 0.
 (F2) r_{0,1}(z; 0) = E/(E + z), and rho r_{1,1}(z; t) = E/t + a1 a2/(E + z) + O(t).
 (F3) Every frame-B dilaton factor Dil_i(d) (i = 1, 2, d >= 1) and every edge factor E_ij(m, n) (i, j in {1, 2}, m, n >= 0) is analytic at t = 0.
 (L1) [t^0](rho q_t) = q_0 = E - E^2 lam_1, and [t^k](rho q_t) lies in Lam_{k+1} for every k >= 1.
 (L2) For j in {1, 2}, every b >= 0 and every k >= 0, [t^k](rho f_{j,b}) lies in Lam_{k+1}. Moreover E_21(b, a) = E_12(a, b), and [t^k] f_{j,b} lies
      in Lam_{k+1} as well.
 (L3) For every k >= 0, [t^k](rho^2 E_11(psi, psi')) lies in span{lam_i(psi) lam_l(psi') : i, l >= 1, (i - 1) + (l - 1) <= k}, and E_11(a, b) = E_11(b, a).
      The same support bound holds for [t^k](rho^e E_11) for every integer e.
Every membership is an identity of formal power series in psi (and psi'), valid at every psi-order and every t-order; the coefficients are
rational functions of the weights, regular on (ND). Sharpness: with P_k := prod_{i=1}^{k} ((k + 1) a2 - i E)((k + 1) A - i E), the top leaf
lam_{k+1} occurs in [t^k](rho q_t) iff P_k != 0, i.e. iff E is not in {(k + 1) a2/i, (k + 1) A/i : 1 <= i <= k}; every entry of the anti-diagonal
(i - 1) + (l - 1) = k of L3 is non-zero iff P_{i-1} P_{l-1} != 0 (Corollaries 11.3, 11.5). Large-z form: coefficientwise in t, z r_{0,1}(z; t) -> E/rho(t)
and E_1j(z, w) -> 0 as z -> oo at fixed order in w, in the sense of Lemma 5.2 and Corollary 9.5.

Remark 0.1 (where the CY condition enters). No step of the proof below uses E = e1 + e2: the frame-B rows are h_{(a1, A, -E)} r^pt and do not
involve e1, e2 at all; Delta_i, E_ij and the pairing carry e1 e2 as an overall factor. Theorem LC_CY therefore holds for every weight point
with (ND), with w := (E, -e1, -e2), E := a0 + a1 + a2. The CY condition is what makes w a sum-zero triple whose first weight is the leaf scale E,
which is what Theorem L and Theorem VI of the local-P2 programme need in Corollary R ([PR] s.13); so it is needed for the regularity theorem
built on Theorem LC_CY, not for Theorem LC_CY itself. (Evidence, not a step: off T_CY the reports in cy_torus_explore/ find the same row
structure and real polar terms in the graph sums.)

### 0.3 Further results proved on the way
**Theorem Q_CY.** (i) For i = 1, 2 the second row of the QDE holds to all orders in z: z D_i v_i = (A - L_i) v_i, with v_i = (r_{0,i}, r_{1,i}),
D_i = d/dt - Delta_i'/(2 Delta_i) and A = [[0, 1], [-a1 a2/tau, E/tau - a0]] (Proposition 4.6; Theorem 8.4 or [PR] Prop. 13.2(b)(ii)).
(ii) <v_i(w), v_j(-w)> = delta_ij Delta_i for all i, j (Proposition 4.8; Corollaries 8.5, 8.7). (iii) v_2(z; 0) = h_{(a1, a2, A, B, -E, -E)}(z) (1, a1 a2/E)
(Proposition 8.6; the scalar factor is the analytic input of s.8.6(a), the direction (1, a1 a2/E) is formal).

**Theorem CF_CY (closed forms; section 11).** In frame B,
    rho r_{0,1}(z; t) = (E/(E + z)) (tau/t) 2F1(1 + a2/z, 1 + A/z; 2 + E/z; tau) Exp[-dphi(tau)/z]                              (KEY FORMULA),
    rho q_t(psi) = E + Sum_{j>=1} A_j(t) lam_j(psi),   A_j(t) = (-1)^j C_j t^-1 zeta(t)^j 2F1(j B/E, j a1/E; j; tau),
    C_j := P_{j-1}/(j! (j - 1)! E^(2j-4)),   zeta := tau e^(dphi/E),
    d/dt [t^2 rho^2 E_11(psi, psi')] = Theta(t) Dg(psi) Dg(psi'),   Theta = e1 e2 (a1 a2 - E L_2) = -e1 e2 B tau s* L_2,   g = -(1/E) Sum_j j A_j lam_j,
with closed forms of every leaf and loop coefficient and of the t = 0 values of E_12, E_22 (s.11).

### 0.4 Consequences not proved here
Corollary R_CY (every coloured stable graph with an edge regular at t = 0 after summing decorations, the single-vertex graph with the single
pole (-e1 e2 E^3)^(g-1) Fpot_g(q_0) t^(2-2g) = beta_{2g-2}(e1, e2) t^(2-2g)) follows from Theorem LC_CY exactly as Corollary R follows from Theorem
LC in [PR] s.13 (Theorem R with c -> E and -x c^3 -> -e1 e2 E^3; M_irr = Theorem L of lemmaA/LEMMA_A_PROOF.md at the sum-zero triple w; Theorem VI by
homogeneity in (e1, e2)); the identification of the graph sum with F_g on T_CY is [PR] Proposition 13.2 restricted to the hyperplane (parts
(a)-(c), (e) are stated there at indeterminate weights; part (d) is not needed when a0 != 0). Neither is written out here.

### 0.5 Inputs beyond algebra
One analytic input enters Theorem LC_CY: [PR] Proposition 3.15 (Laplace's method, [PR] Lemma 3.8, and Stirling's series with remainder), used
here only through its rescaling (Proposition 3.15_CY). A second input of the same kind, the Q = 0 constant h_nw of [PR] Lemma 13.3 (same method),
enters the definition of the rows and the constant line of s.3.9. A third, r^pt_{0,2}(z; 0) = h_{(a2, -E, B)} (Proposition 8.6(a), again by the
method of [PR] Lemma 13.3), is used only for the t = 0 values at the regular idempotent and for K_2 = 1, not for Theorem LC_CY. Everything
else is formal, plus the analyticity in tau of finitely many explicit algebraic functions.

## 1. Conventions

### 1.1 Parameters and variables
Weights as above with (ND); Kw := Q(a0, a1, a2, e1, e2), and K := Kw restricted to T_CY, i.e. Q(a1, a2, e1, e2) with a0 := e1 + e2 - a1 - a2
(the frame-B rows lie in Q(a0, a1, a2)-coefficient series and never see e1, e2 separately, so nothing depends on the choice). Q = e^t,
tau = 1 - 1/Q, derivatives, [t^k], [z^m] and X^(k) as in [PR] s.1.1 (verbatim).

### 1.2 The factors at generic weights
This replaces [PR] s.1.2 (C1)-(C14). The generic data are taken from [PR] Prop. 13.2(b)(i) and Lemma 13.3 (the code rc_lib.m, rmatrixAtClosed, was read
afterwards in the referee pass and re-read for this revision; it agrees with (C1'), (C2'), (C4'), (C5'), (C7') line by line, s.13.1): L(U) = -a1 Q U/(1 - Q U) = du/dt, U/(1 - U) = (L - a2)/A, and the Q -> 0 phases -a2 log U + A log(1 - U)
at U_+ and a0 log V + a1 log(1 - V) at V = Q U_-. These determine the superpotential up to an additive constant:
 (C1') W(U) := a1 log(1 - Q U) - a2 log U + A log(1 - U); theta := U d/dU; W_k := theta^k W (k >= 1), single valued, rational in U and Q.
       (On T_F: x log(1 - Q U) - log U + log(1 - U), [PR] (C1).) Indeed partial_t W = -a1 Q U/(1 - Q U) = L(U), and the critical equation
       below gives U W_U = 0 iff L = a2 + A U/(1 - U), i.e. U/(1 - U) = (L - a2)/A.
 (C2') Critical points: U W_U = -a1 Q U/(1 - Q U) - a2 - A U/(1 - U) = 0 iff
           LeqU := B Q U^2 - (a0 + Q (a1 - a2)) U - a2 = 0
       (multiply by (1 - Q U)(1 - U): -a1 Q U (1 - U) - a2 (1 - Q U)(1 - U) - A U (1 - Q U) = -[Q U^2 (a1 - a2 + A) - U (a1 Q - a2 Q - a2 + A) - a2], and
       a1 - a2 + A = B, a1 Q - a2 Q - a2 + A = a0 + Q (a1 - a2)). Roots: U_1 + U_2 = (a0 + Q (a1 - a2))/(B Q) (Ssum), U_1 U_2 = -a2/(B Q) (Pprod),
       discriminant (a0 + Q (a1 - a2))^2 + 4 a2 B Q = Q^2 delta^2 with delta^2 of s.1.3. U_1 -> 1 and U_2 -> -a2/B as t -> 0 (s.1.3, Lemma 8.2).
 (C3') Quadratic algebra QA := K(Q)[U]/(LeqU) and conjugate encoding U -> Ssum - U, exactly as [PR] (C3); the stored quantities of idempotent 2
       are cconj of those of idempotent 1, and the scalar series h_nw, h_w^-1 below are the same at both idempotents (Lemma 13.3 of [PR]), so
       cconj commutes with the twists. The true value of a factor of idempotent 2 is ev_1 of its stored element.
 (C4') Flat basis {1, H}. Frobenius data (Lemma 4.1): eta = (1/(A B e1 e2)) [[E/(a1 a2), 1], [1, -a0]], eta^-1 = e1 e2 [[a0 a1 a2, a1 a2], [a1 a2, -E]],
       constant in t. Covectors a = (a_0, a_1) are paired by
           <a, b> := e1 e2 [ a0 a1 a2 a_0 b_0 + a1 a2 (a_0 b_1 + a_1 b_0) - E a_1 b_1 ],
       symmetric with constant coefficients, so d/dt <a, b> = <a', b> + <a, b'>. (On T_F: x^2 (a_0 b_1 + a_1 b_0) - c x a_1 b_1.)
 (C5') Eigenvalue of H*: L(U) := -a1 Q U/(1 - Q U); L_i := L(U_i).
       The code's LofU (rc_lib.m line 81) is -(Q B U + Q a2)/(1 - Q) = (B U + a2)/tau (as 1 - Q = -Q tau), linear in U, not -a1 Q U/(1 - Q U); the difference is
       (B U + a2)/tau + a1 Q U/(1 - Q U) = -LeqU/(tau (1 - Q U)) (numerator: (B U + a2)(1 - Q U) + a1 Q U tau = -B Q U^2 + (B - a1) U + Q (a1 - a2) U + a2 = -LeqU,
       using a1 Q U tau = a1 U (Q - 1) and B - a1 = a0), and 1 - Q U is a unit of QA (LeqU(1/Q) = (B - a0)/Q - a1 = -a1 tau != 0 for tau != 0). So the two agree at
       both roots and are equal in QA: the L_i, and Delta_i = DeltaL(L_i) (line 128), are the same.
 (C6') Delta_i := 1/eta(e_i, e_i) (Proposition 4.3). Vertex weight Delta_i^(g_v - 1).
 (C7') Point-frame identity row: r^pt_{0,i}(z; t) := Sum_{m>=0} z^m wickZ_m(W_2(U_i), W_3(U_i), ...), wickZ_m exactly as [PR] (C7) (verbatim; it is a
       universal polynomial in W_2^-1, W_3, ..., W_{2m+2} with rational coefficients, every monomial W_2^-n prod W_{k_i} having Sum k_i = 2n).
       The row with the Q = 0 constant (rmatrixAtClosed; [PR] Prop. 13.2(b), Lemma 13.3):
           r_{0,i}(z; t) := h_nw(z) r^pt_{0,i}(z; t),   nw := (a1, A, -e1, -e2),   h_v(z) := Exp[-Sum_{k>=1} bk_k z^(2k-1) Sum_a v_a^(1-2k)],
       bk_k := B_2k/(2k (2k - 1)); h_{u+v} = h_u h_v, h_{-v} = h_v^-1, h_v(z) h_v(-z) = 1. h_nw is the same scalar at both idempotents ([PR] Lemma 13.3,
       stated and proved there at indeterminate weights; it needs a0 != 0 for the Q = 0 limit, see s.12 for a0 = 0). On T_F, nw = (x, 1, -x, -1) and h_nw = 1.
 (C8') Frame B: every row is multiplied by h_w(z)^-1 with w := (E, -e1, -e2):
           r'_{a,i}(z; t) := h_w(z)^-1 r_{a,i}(z; t) = h_w^-1 h_nw r^pt_{a,i} = h_{(a1, A, -E)}(z) r^pt_{a,i}(z; t),
       because (a1, A, -e1, -e2) + (-E, e1, e2) = (a1, A, -E) as multisets. So the frame-B rows depend on (a0, a1, a2) only. Below the prime is
       dropped: r_{a,i} denotes the frame-B row unless "point frame" or "R_Z" is said. On T_F, h_{(x, 1, -c)} = h_{(c, -x, -1)}^-1 = h_w^-1, [PR] (C8).
 (C9') H-row: r_{1,i}(z) = (L_i + z D_i) r_{0,i}(z), D_i := d/dt - Delta_i'/(2 Delta_i) (verbatim [PR] (C9); linear over C[[z]] with t-independent
       coefficients, so it holds in every frame). v_i(z) := (r_{0,i}(z), r_{1,i}(z)).
 (C10')-(C12') Dilaton factors Dil_i(n) := (-1)^n r_{0,i}[n-1], q_t = -psi r_{0,1}(-psi; t); edge factors EdN and Lemma 1.1 (E_ij(z, w) =
       [P_ij(z, w) - P_ij(-w, w)]/(z + w), P_ij(z, w) = -<v_i(-z), v_j(-w)>); legs f_{j,b}: verbatim [PR] (C10)-(C12) and Lemma 1.1, whose proof uses no
       weights. Dil_1(1) = -1 at every t since r_{0,1}[0] = 1 (wickZ_0 = 1 and h_{(a1,A,-E)}[0] = 1).
 (C13') rho := (Delta_1 t^2/(-e1 e2 E^3))^(-1/2) with rho(0) = 1; lam_j(psi) := 1/(E - j psi) (lam_0 = 1/E); in z = -psi, lam_j = 1/(E + j z).
 (C14') Scope: the formal series in z, w defined by (C7')-(C11') for all m, n.

### 1.3 Algebraic functions of tau
- delta^2 := (E + a0 tau)^2 - 4 A B tau = E^2 - (E^2 + u) tau + a0^2 tau^2 (the two forms agree since 4 A B = 4 a1 a2 + 4 a0 E and
  2 a0 E + 4 a1 a2 = E^2 + u; indeed u = E^2 - 2 E (a1 + a2) + 4 a1 a2 = E^2 - 2 E (E - a0) + 4 a1 a2). delta := the root analytic near tau = 0 with
  delta(0) = E; it is analytic on the disc D := {|tau| < r_D}, r_D := the smallest modulus of a root of delta^2 (r_D is finite under (ND): for a0 != 0 delta^2 is
  quadratic in tau, and at a0 = 0, where E^2 + u = 4 a1 a2 != 0, delta^2 = E^2 - 4 a1 a2 tau).
  (On T_F: delta^2 = c^2 - 4 x tau, [PR] s.1.3.)
- s* := (E + a0 tau - delta)/(2 B tau) = 2 A/(E + a0 tau + delta), analytic on D with s*(0) = A/E: the denominator E + a0 tau + delta vanishes only
  where (E + a0 tau)^2 = delta^2, i.e. 4 A B tau = 0, i.e. tau = 0, where it equals 2 E != 0. s~ := (E + a0 tau + delta)/(2 B tau) = A/(B tau s*);
  s*, s~ are the two roots of N(s) := B tau s^2 - (E + a0 tau) s + A (sum (E + a0 tau)/(B tau), product A/(B tau)); tau s~ is analytic with value E/B at 0.
  N'(s*) = 2 B tau s* - (E + a0 tau) = -delta, N'(s~) = +delta. Also delta = E + a0 tau - 2 B tau s*.
- U_1 = 1 - tau s* -> 1 and U_2 = 1 - tau s~ -> 1 - E/B = -a2/B (Lemma 4.2(c)). U_1 is the root with the singular Delta (Proposition 4.3(d)).
- Discs: D' := D intersected with {|tau| < 1}; D'' := D' - {0}.

### 1.4 Rings and general principles
Verbatim [PR] s.1.4 with Kx replaced by Kw (or K): O, M, Kz, Rz, Rz0, iota, refl, TRANSFER, the leaf spaces Lam_m (dimension m because the
poles E/j are distinct, E != 0), (LAM1)-(LAM6) with c -> E (they are identities of partial fractions in which c is a free constant), (PF) with
S = {-E/j : 1 <= j <= M} (distinct and non-zero since E != 0), and (T). The passage (P) between complex and indeterminate weights: every object
below (tau, s*, delta, dphi, alg, rho, L_i, Delta_i, Theta, the rows, h_{(a1,A,-E)}, the factors) is built by field operations, d/dt and the
expansions of the square roots delta, alg (units at tau = 0) and of logarithms of units from the weights and tau; at tau = 0 the only
denominators are a1, a2, A, B, E, e1, e2 (s*(0) = A/E, 1 - s*(0) = a1/E, U_2(0) = -a2/B, 1 - U_2(0) = E/B, phi_0''(A/E) = -E^3/(a1 A),
Delta_2(0) = e1 e2 a1 a2 A B/E, the Gauss denominators E + j z). So the t-, z-, psi-coefficients are rational functions of the weights, regular
on (ND), and an identity proved over Kw holds at every complex point of (ND) and conversely. Sections 5-7 are written over Kw.

### 1.5 The closed-form ingredients
- Euler phase: phi_tau(s) := A log s + a1 log(1 - s) - a2 log(1 - tau s). phi_0(A/E) = A log(A/E) + a1 log(a1/E), phi_0''(A/E) = -E^3/(a1 A) (Lemma 3.1).
- dphi(tau) := phi_tau(s*(tau)) - phi_0(A/E), principal branches near (tau, s) = (0, A/E); dphi lies in tau O, with
      dphi_tau = a2 s*/(1 - tau s*) = (A - B s*)/(1 - tau),   dphi' = A - B s* = a2 - L_2   (Lemma 3.1(b)),
      dphi = (a2 A/E) tau + O(tau^2) = (a2 A/E) t - (a1 a2 A B/(2 E^3)) t^2 + O(t^3)   (the t^2 term from Theta_1, s.7.3).
- alg(tau) := (1 - tau s*) (phi_tau''(s*)/phi_0''(A/E))^(1/2), the branch with alg(0) = 1; alg^2 = a1 A delta U_1/(E^3 s* (1 - s*)) = A delta (E - B tau s*)/(E^3 s*) (Lemma 3.1(c), (3.4), (3.4')).
- Gauss's hypergeometric series: alpha := 1 + a2/z, beta := 1 + A/z, gamma := 2 + E/z, F(z, tau) := 2F1(alpha, beta; gamma; tau) = Sum_n F_n(z) tau^n,
      F_n = (alpha)_n (beta)_n/((gamma)_n n!) = prod_{j=1}^{n} (a2 + j z)(A + j z) / ( n! z^n prod_{j=2}^{n+1} (E + j z) )   in Kw(z),
  using (1 + v/z)_n = z^-n prod_{j=1}^{n} (v + j z) and (2 + E/z)_n = z^-n prod_{j=2}^{n+1} (E + j z). Gauss[Y] := tau (1 - tau) Y_tautau + (gamma - (alpha + beta + 1) tau) Y_tau - alpha beta Y,
  with alpha + beta + 1 = 3 + (a2 + A)/z and alpha beta = 1 + (a2 + A)/z + a2 A/z^2.
- Ephi := Exp[-dphi/z] in C[1/z][[tau]]; X_n := [tau^n] Ephi, a polynomial of degree <= n in 1/z. H_n := [tau^n](F Ephi), Psi_n := [tau^n](F Ephi alg).
- G(z; t) := (E/(E + z)) (tau/t) F Ephi in Kw(z)[[t]]; G_k := [t^k] G; g(psi; t) := G(-psi; t); Dg := d/dt (t g).
- QDE data: A := [[0, 1], [-a1 a2/tau, E/tau - a0]] on covectors, (a_0, a_1) -> (a_1, -(a1 a2/tau) a_0 + (E/tau - a0) a_1); chi_i := (1, L_i);
  pi(a) := <a, chi_2>; mu := E/tau - A; omega := (E + z)/tau - (A + z) = mu + z (1 - tau)/tau.
- Theta := Delta_2/(L_2 - L_1) = e1 e2 (a1 a2 - E L_2) = -e1 e2 B tau s* L_2 = e1 e2 (E dphi' - a2 A) (Theorem 7.2).
- sigma(psi) := rho v_1(-psi) = (sigma_0, sigma_1); sigma^(k) := [t^k] sigma. Leaf coefficients q^(k) := [t^k](rho q_t), a_j^(k), g_j^(k) := -(j/E) a_j^(k),
  A_j(t) := Sum_k a_j^(k) t^k, Ahat_j := t A_j, e^(k)_il and e_il(t): as [PR] s.1.5 with c -> E.
- P_k := prod_{i=1}^{k} ((k + 1) a2 - i E)((k + 1) A - i E) (P_0 = 1); C_j := P_{j-1}/(j! (j - 1)! E^(2j-4)) (C_1 = E^2); zeta(t) := tau e^(dphi(tau)/E).
  On T_F: P_k = prod (k + 1 - i c)^2, [PR] s.1.5.

## 2. Plan of the proof and logical order
The chain of [PR] s.2 is kept step by step and label by label (3.1 -> 3.1_CY etc.; the suffix is dropped when no confusion is possible).
Verbatim statements are cited; the substituted ones are re-derived in s.3-s.9; new lines are: the rescaling of Proposition 3.15 and the
constant line of s.3.9 (where h_nw enters), the identity (delta + tau L_2)(1 - s*) = a1 U_1 replacing [PR]'s 1 - s = x s (1 - tau s)
(Lemma 3.1(b)), the form R_tr = delta (log alg)_tau in s.3.10, and the generic Theta and pi of s.4.7, s.6.3, s.7.3. Row 2 of the QDE at the
regular idempotent and the symplectic relations may alternatively be taken from [PR] Prop. 13.2(b)(ii), (iv) (s.8.4, s.8.5: what those cover is
stated there). Nothing about the code's numerical output is used.

## 3. The singular identity row from the Euler integral (Steps 1-2; F2)

Statements are [formal] (identities of formal series, valid for every weight point with (ND)) or [analytic].

### 3.1 The Euler coordinate [substitution]
For tau != 0 put U = 1 - tau s. As in [PR] s.3.1, 1 - U = tau s, 1 - Q U = (1 - Q)(1 - s) = -Q tau (1 - s), dU/U = -tau ds/(1 - tau s). Hence, with
compatible branches,
    W(1 - tau s) = a1 log(-Q tau (1 - s)) - a2 log(1 - tau s) + A log(tau s) = phi_tau(s) + E log tau + a1 log Q + a1 log(-1),           (3.1)
    e^{W/z} dU/U = kappa_ph(z) Q^(a1/z) tau^(E/z + 1) e^{phi_tau(s)/z} (1 - tau s)^-1 ds,   kappa_ph(z) := -e^{a1 log(-1)/z},          (3.2)
since a1 + A = E. In Pochhammer form e^{phi_tau/z} (1 - tau s)^-1 = s^(A/z) (1 - s)^(a1/z) (1 - tau s)^(-a2/z - 1). The logarithmic points U = 1, 1/Q, 0
are s = 0, 1, 1/tau; at t = 0 the first two collide and U_1 runs into them, while in s the critical point stays at s*(0) = A/E strictly
between 0 and 1 when a1, A > 0. (On T_F this is [PR] (3.1), (3.2).)

### 3.2 Critical points, Hessians, the tensor rule [substitution]
**Lemma 3.1_CY.** Let tau != 0 and N(s) := B tau s^2 - (E + a0 tau) s + A.
 (a) phi_tau'(s) = A/s - a1/(1 - s) + a2 tau/(1 - tau s) = N(s)/(s (1 - s)(1 - tau s)). The critical points of phi_tau in C - {0, 1, 1/tau} are s* and s~,
     and 1 - tau s* = U_1, 1 - tau s~ = U_2. On D', s* is not in {0, 1, 1/tau}; on D'', s~ is not in {0, 1, 1/tau} and U_2 is not in {0, 1, 1/Q}.
 (b) At s = s* (and, for (ii)-(iv), also at s = s~, with delta -> -delta in (i), (iii)):
     (i) 2 B tau s* - (E + a0 tau) = -delta;   (ii) (1 - tau) a2 s = (A - B s)(1 - tau s), equivalently (E - B tau s)(1 - s) = a1 (1 - tau s);
     (iii) ds*/dtau = s* L_2/delta with L_2 := B s* - a0;   (iv) A - E s = -tau s (B s - a0).
     Consequently dphi_tau = a2 s*/(1 - tau s*) = (A - B s*)/(1 - tau) and dphi' = A - B s* (envelope formula).
 (c) Hessians. At a root s_c of N, phi_tau''(s_c) = N'(s_c)/(s_c (1 - s_c)(1 - tau s_c)); so
        phi_tau''(s*) = -delta/(s* (1 - s*) U_1),   phi_tau''(s~) = delta/(s~ (1 - s~) U_2),   phi_0''(A/E) = -E^3/(a1 A).                       (3.3)
     Both critical points are nondegenerate on D'', and
        alg^2 = U_1^2 phi_tau''(s*)/phi_0''(A/E) = a1 A delta U_1/(E^3 s* (1 - s*)),                                                    (3.4)
     a unit of O with alg(0) = 1. Equivalently, by (b)(ii) in its second form (U_1/(1 - s*) = (E - B tau s*)/a1),
        alg^2 = A delta (E - B tau s*)/(E^3 s*),   E - B tau s* = delta + tau L_2,                                                              (3.4')
     and E - B tau s* has no zero on D: if B tau s* = E then N(s*) = s* (B tau s* - E - a0 tau) + A = -a0 E/B + A = a1 a2/B != 0. Since delta and s* are units
     of O(D) as well (s.1.3), alg^2 is a unit of O(D) and 1/alg is analytic on D.
 (d) Tensor rule (verbatim [PR] Lemma 3.1(d), no weights): d/ds G(1 - tau s) = -(tau/(1 - tau s)) (theta G)(1 - tau s), and at a critical point
        W_2(U_i) = U_i^2 phi_tau''(s_i)/tau^2,   W_2(U_1) = -delta U_1/(tau^2 s* (1 - s*)),   W_2(U_2) = delta U_2/(tau^2 s~ (1 - s~)).         (3.5)
Proof. (a) A (1 - s)(1 - tau s) - a1 s (1 - tau s) + a2 tau s (1 - s) = A - (A + a1) s - A tau s + a2 tau s + (A + a1 - a2) tau s^2 = N(s), using
A + a1 = E, a2 - A = -a0, A + a1 - a2 = B. The roots of N are s*, s~ (s.1.3), and 1 - tau s* = U_1, 1 - tau s~ = U_2 by definition. s = 0 is not a
root (N(0) = A != 0); N(1) = B tau - E - a0 tau + A = (B - a0) tau - a1 = a1 (tau - 1), so a root equals 1 only if tau = 1 (a1 != 0); N(1/tau) =
B/tau - E/tau - a0 + A = (B - E)/tau + a2 = a2 (1 - 1/tau), so a root equals 1/tau only if tau = 1 (a2 != 0). On D'' both roots are finite (B tau != 0).
U_2 = 0 means s~ = 1/tau, U_2 in {1, 1/Q} means s~ in {0, 1} (1 - Q U_2 = -Q tau (1 - s~)).
(b) (i) is N'(s*) = -delta (s.1.3). (ii) First form: (1 - tau) a2 s - (A - B s)(1 - tau s) = a2 s - a2 tau s - A + A tau s + B s - B tau s^2 =
-[B tau s^2 - (B + a2) s + (a2 - A) tau s + A] = -N(s), since B + a2 = E and a2 - A = -a0. Second form: by (i), delta = E + a0 tau - 2 B tau s*, and
tau L_2 = B tau s* - a0 tau, so E - B tau s* = delta + tau L_2, and (E - B tau s)(1 - s) - a1 (1 - tau s) = E - E s - B tau s + B tau s^2 - a1 + a1 tau s
= A - E s + B tau s^2 - (B - a1) tau s = N(s). (iii) Differentiate N(s*(tau), tau) = 0: B s*^2 + 2 B tau s* s*_tau - a0 s* - (E + a0 tau) s*_tau = 0, i.e.
N'(s*) s*_tau = -(B s*^2 - a0 s*) = -s* L_2; with (i), s*_tau = s* L_2/delta. (iv) is N(s) = 0 rearranged: A - E s = a0 tau s - B tau s^2.
Envelope: d/dtau [phi_tau(s*(tau))] = (d_tau phi_tau)(s*) + phi_tau'(s*) s*_tau = (d_tau phi_tau)(s*) = a2 s*/(1 - tau s*), since
d_tau phi_tau(s) = a2 s/(1 - tau s); then (ii), and d/dt = (1 - tau) d/dtau.
(c) Since phi' = N/(s (1 - s)(1 - tau s)) and N(s_c) = 0, phi''(s_c) = N'(s_c)/(s_c (1 - s_c)(1 - tau s_c)) (the derivative of the denominator
multiplies N(s_c) = 0). N'(s*) = -delta, N'(s~) = delta give (3.3); at tau = 0, phi_0''(A/E) = -A (E/A)^2 - a1 (E/a1)^2 = -E^2 (1/A + 1/a1) = -E^3/(a1 A).
delta has no zero on D, s*, 1 - s*, U_1 are units of O (values A/E, a1/E, 1), so phi''(s*) is a unit of O, and (3.4) follows: alg^2 =
U_1^2 (-delta/(s* (1 - s*) U_1)) (-a1 A/E^3) = a1 A delta U_1/(E^3 s* (1 - s*)), with value E (a1 A/E^3) (E/A)(E/a1) = 1 at tau = 0; alg is its square root
with value 1 (two square roots of a unit of O with the same value at 0 coincide). (3.4') is (3.4) with U_1/(1 - s*) = (E - B tau s*)/a1, the second form of (b)(ii);
the value of N at s = E/(B tau) is direct.
(d) Verbatim [PR]: d/ds = -(tau/U) theta, theta W = 0 at a critical point, so d^2/ds^2 [W(1 - tau s)] = (tau/U)^2 W_2; with (3.1) and (3.3) this is (3.5). QED
[Evidence: check1_identities.py (s.13.3): every identity of (a)-(d) at 25 random rational points with rational delta, both roots.]

Remark 3.1. On T_F, (b)(ii) is 1 - s = x s (1 - tau s) ([PR] Lemma 3.1(b)(ii)) and (3.4) is delta/(c^3 s*^2). Generically there is no factorisation
of N(s) of the form (1 - s) - x s (1 - tau s); the two forms of (ii) are what replaces it in s.3.10 and s.4.

### 3.3 The rescaled derivatives tau^k W_k(U_1) [substitution]
For k >= 1, with (3.6) of [PR] generalised,
    phi_tau^(k)(s) = (k - 1)! [ (-1)^(k-1) A s^-k - a1 (1 - s)^-k + a2 tau^k (1 - tau s)^-k ],   phi_0^(k)(A/E) = (k - 1)! E^k [ (-1)^(k-1) A^(1-k) - a1^(1-k) ].   (3.6)
With p_a := (1 - a U)^-1 (a in {1, Q}), theta log(1 - a U) = 1 - p_a, theta p_a = p_a^2 - p_a, Pol_0(p) := p, Pol_n(p) := (p^2 - p) Pol_{n-1}'(p), and
theta^k log(1 - a U) = -Pol_{k-1}(p_a) for k >= 2, theta^k log U = 0 for k >= 2 (verbatim [PR] s.3.3). Hence
    W_k = -a1 Pol_{k-1}(p_Q) - A Pol_{k-1}(p_1)   (k >= 2),                                                                            (3.7)
Pol_n of degree n + 1, leading coefficient n!, Pol_n(0) = 0; Pol_{k-1}(p) = Sum_{i=1}^{k} pc_{k-1,i} p^i, pc_{k-1,k} = (k - 1)!.
**Lemma 3.2_CY.** For k >= 2 and tau in D'' put wr_k(tau) := tau^k W_k(U_1). Then
    wr_k = -Sum_{i=1}^{k} pc_{k-1,i} tau^(k-i) [ a1 (-(1 - tau)/(1 - s*))^i + A s*^(-i) ],                                                (3.8)
analytic on D' (tau = 0 included), with
    wr_k(0) = -(k - 1)! E^k [ (-1)^k a1^(1-k) + A^(1-k) ] = (-1)^k phi_0^(k)(A/E),   wr_2 = -delta U_1/(s* (1 - s*)) (no zero on D').
Proof. As [PR] Lemma 3.2: at U_1, 1 - U_1 = tau s* and 1 - Q U_1 = -Q tau (1 - s*), so tau p_1(U_1) = 1/s* and tau p_Q(U_1) = -(1 - tau)/(1 - s*);
multiply (3.7) by tau^k. s*, 1/s*, 1/(1 - s*) are analytic on D' (Lemma 3.1(a)). At tau = 0 only i = k survives, with s*(0) = A/E and 1 - s*(0) = a1/E:
wr_k(0) = -(k - 1)! [a1 (-E/a1)^k + A (E/A)^k] = -(k - 1)! E^k [(-1)^k a1^(1-k) + A^(1-k)], which is (-1)^k phi_0^(k)(A/E) by (3.6). wr_2 = tau^2 W_2(U_1) is
(3.5). QED

### 3.4 Formal Laplace expansions [verbatim]
[PR] s.3.4 (the ring Acal, the Gaussian functional, Lemma 3.3 (G0)-(G2), Definition 3.4 with (D1)-(D3), Proposition 3.5 (P3), Proposition 3.6
(P2), Proposition 3.7 (P1), Lemma 3.8 (Laplace's method) and Lemma 3.9 (homogeneity)) contains no weights and is used verbatim. Two
consequences of Definition 3.4 that are used below and are immediate from (D1):
 (D4) rescaling of the phase: for a constant lambda != 0 and am independent of z, Lnorm[lambda f, am; u0](z) = Lnorm[f, am; u0](z/lambda), since
      Rr_{lambda f}/z = Rr_f/(z/lambda) and (-z/(lambda f2))^(p/2) = (-(z/lambda)/f2)^(p/2);
 (D5) the parameters: by (D2) and (3.6), the z^m coefficient of Lnorm[phi_0, 1; A/E] is a polynomial with rational coefficients in
      1/phi_0''(A/E) = -a1 A/E^3 and in the phi_0^(k)(A/E) of (3.6), hence an element of Q(a1, A) = Q(a1, A, E), homogeneous of degree -m in the weights.

### 3.5 The code's row is the normalised expansion of the Euler integrand [verbatim after (C1')]
**Proposition 3.10_CY.** (i) For every tau in D'' , coefficientwise in z,
    r^pt_{0,1}(z; t) = N(z, tau) := Lhat[phi_tau, (1 - tau s)^-1; s*(tau)] = U_1 Lnorm[phi_tau, (1 - tau s)^-1; s*(tau)].                      (3.12)
(ii) [z^m] r^pt_{0,1}(z; t) = wickZ_m(W_2(U_1), ...) = wickZ_m(wr_2, ..., wr_{2m+2}), a polynomial in 1/wr_2, wr_3, ... with rational coefficients, hence
analytic on all of D', with value at tau = 0 equal to Abeta_m := [z^m] Abeta(z), Abeta(z) := Lnorm[phi_0, 1; A/E].
(iii) N lies in O(D')[[z]], N(z, 0) = Abeta(z), and (3.12) holds in C[[z]][[t]].
Proof. Word for word the proof of [PR] Proposition 3.10 with (3.1), (3.5), Lemma 3.2_CY and (3.6) in place of their T_F versions: the code's row
is Lnorm[Wy, 1; y1] in y = log U with Wy''(y1) = W_2(U_1) != 0 (by (3.5), delta U_1 != 0 on D''); the coordinate change hc(s) = log(1 - tau s) has
hc'(s*) = -tau/U_1 != 0, Wy o hc = phi_tau + const by (3.1), and (P1)(iii), (D3) give (i); (ii) is Lemma 3.9 with lambda = tau, then Lemma 3.2_CY, then
Lemma 3.9 with lambda = -1 and (D1) at tau = 0 (phi_0'(A/E) = 0, phi_0''(A/E) = -E^3/(a1 A) != 0); (iii) is (D2) with the analyticity of s* on D',
s* not in {0, 1, 1/tau}, phi''(s*) != 0, and s*(0) = A/E, (1 - 0 s)^-1 = 1. QED
Remark 3.11 of [PR] (the uniform coordinate ytilde = -log(1 - tau s)/tau) carries over verbatim with (3.1).

### 3.6 Gauss's equation for the formal expansion [substitution]
The amplitude identity. Put f := phi_tau, am0 := (1 - tau s)^-1, Dam q := d_tau q + q (d_tau f)/z, so e^{f/z} Dam q = d_tau (e^{f/z} q). Since
d_tau f = a2 s/(1 - tau s),
    Dam am0 = alpha s (1 - tau s)^-2,   Dam^2 am0 = alpha (alpha + 1) s^2 (1 - tau s)^-3   (alpha = 1 + a2/z),
and with kE := s (1 - s)(1 - tau s)^-2,
    tau (1 - tau) Dam^2 am0 + (gamma - (alpha + beta + 1) tau) Dam am0 - alpha beta am0 = -alpha (kE' + f' kE/z).                        (3.16)
Proof of (3.16). It is an identity of Laurent polynomials in z with coefficients rational in (s, tau), valid identically (not only at a critical
point). Multiply by (1 - tau s)^3/alpha; since f' kE = N(s)/(1 - tau s)^3 (Lemma 3.1(a)) and kE' = [(1 - 2s)(1 - tau s) + 2 tau s (1 - s)]/(1 - tau s)^3,
the claim is
    tau (1 - tau)(alpha + 1) s^2 + (gamma - (alpha + beta + 1) tau) s (1 - tau s) - beta (1 - tau s)^2 = -(1 - 2s)(1 - tau s) - 2 tau s (1 - s) - N(s)/z.
Coefficient of 1/z on the left: tau (1 - tau) a2 s^2 + (E - (a2 + A) tau) s (1 - tau s) - A (1 - tau s)^2. Expanding,
  tau a2 s^2 - tau^2 a2 s^2 + E s - E tau s^2 - (a2 + A) tau s + (a2 + A) tau^2 s^2 - A + 2 A tau s - A tau^2 s^2
  = -A + [E + (A - a2) tau] s + (a2 - E) tau s^2 + (-a2 + a2 + A - A) tau^2 s^2 = -A + (E + a0 tau) s - B tau s^2 = -N(s).
Coefficient of z^0 on the left: 2 tau (1 - tau) s^2 + (2 - 3 tau) s (1 - tau s) - (1 - tau s)^2 = 2 tau s^2 - 2 tau^2 s^2 + 2 s - 3 tau s - 2 tau s^2 + 3 tau^2 s^2
  - 1 + 2 tau s - tau^2 s^2 = 2 s - tau s - 1; on the right -(1 - 2s)(1 - tau s) - 2 tau s (1 - s) = -1 + tau s + 2 s - 2 tau s^2 - 2 tau s + 2 tau s^2 = 2 s - tau s - 1. QED
[Evidence: check1_identities.py, 30 random rational points (weights, s, tau, z).] (3.16) is Euler's integration by parts for the integrand
s^(beta-1) (1 - s)^(gamma-beta-1) (1 - tau s)^-alpha written for amplitudes ([PR] Remark 3.20 with a = alpha, b = beta, g = gamma); on T_F it is [PR] (3.16).

**Proposition 3.12_CY.** Let L(tau) := Lnorm[phi_tau, (1 - tau s)^-1; s*(tau)] in O(D')[[z]] and NG := Exp[dphi/z] (phi_tau''(s*)/phi_0''(A/E))^(-1/2) iota(L).
Then Gauss[NG] = 0 in Kz[[tau]], NG(z, 0) = Abeta(z) in 1 + z Kw[[z]], and NG = Exp[dphi/z] N/alg.
Proof. Verbatim [PR] Proposition 3.12: (P2) twice with f0 = phi_tau(s*), f0_tau = dphi_tau, f2 = phi''(s*), (D3), the amplitude identity (3.16) and
(P3) give [tau (1 - tau) Dop^2 + (gamma - (alpha + beta + 1) tau) Dop - alpha beta] L = -alpha Lnorm[kE' + f' kE/z] = 0; conjugation by
E := Exp[dphi/z] (f2/f2(0))^(-1/2) turns Dop into d/dtau. At tau = 0, NG(z, 0) = Lnorm[phi_0, 1; A/E] = Abeta, in 1 + z Kw[[z]] by (D5). N = U_1 L and
alg = U_1 (f2/f2(0))^(1/2). QED

### 3.7 Uniqueness and the closed form [substitution in one constant]
**Lemma 3.13_CY.** If Y in Kz[[tau]] satisfies Gauss[Y] = 0 then Y = Y_0 F.
Proof. [PR] Lemma 3.13 with the recursion (k + 1)(k + gamma) Y_{k+1} = (k + alpha)(k + beta) Y_k, where k + gamma = (E + (k + 2) z)/z is a unit of
Kz because E != 0 (this is the only place where the parameters enter). QED
So the solutions form the line Kz F; the second local exponent 1 - gamma = -1 - E/z is not in this ring. (The exponents of Gauss's equation
at tau = 0 are 0 and 1 - gamma; those of the Picard-Fuchs operator z^2 (Q - 1) theta^2 - z (a0 + (a1 + a2) Q) theta + a1 a2 Q at Q = 1 are 0 and
1 + E/z: only E enters at Q = 1. Not used.)

**Theorem 3.14_CY.** In Kz[[tau]]: N(z, tau) = Abeta(z) F(z, tau) Ephi(z, tau) alg(tau).                                                  (3.17)
Proof. Verbatim [PR] Theorem 3.14 from Proposition 3.12_CY and Lemma 3.13_CY. QED

### 3.8 The constant Abeta(z) [substitution: Proposition 3.15 rescaled]
**Proposition 3.15_CY.** Abeta(z) := Lnorm[phi_0, 1; A/E] = h_{(E, -a1, -A)}(z) E/(E + z) in Q(a1, A)[[z]]. In words: the normalised Laplace
expansion of the Beta integral B(A/z + 1, a1/z + 1) = Int_0^1 s^(A/z) (1 - s)^(a1/z) ds = Int_0^1 e^{phi_0(s)/z} ds is the Stirling tail
h_{(E, -a1, -A)} of B(A/z, a1/z) times the elementary factor E/(E + z) from B(a + 1, b + 1) = B(a, b) a b/((a + b)(a + b + 1)).
Proof. Put x := a1/A (x != 0, -1 by (ND): a1 != 0 and a1 + A = E != 0) and ph(s) := log s + x log(1 - s), the Euler phase of [PR] at tau = 0 with
parameter x; its critical point is 1/(1 + x) = A/E. Then phi_0 = A ph, so by (D4) Abeta(z) = Lnorm[ph, 1; 1/(1 + x)](z/A). [PR] Proposition 3.15,
an identity in Q(x)[[z']] valid for every x != 0, -1, gives Lnorm[ph, 1; 1/(1 + x)](z') = h_{(1 + x, -x, -1)}(z') (1 + x)/(1 + x + z'). At z' = z/A the
exponent of h is -Sum_k bk_k (z/A)^(2k-1) Sum_a w_a^(1-2k) = -Sum_k bk_k z^(2k-1) Sum_a (A w_a)^(1-2k) with A w = (E, -a1, -A), i.e. h_{(E, -a1, -A)}(z); and
(1 + x)/(1 + x + z/A) = E/(E + z). Both sides lie in Q(a1, A)[[z]] by (D5). QED
Check at z^1, by hand (as in [PR]): by (D1), Abeta_1 = f4/(8 f2^2) - 5 f3^2/(24 f2^3) with f_k := phi_0^(k)(A/E) from (3.6): f2 = -E^3/(a1 A),
f3 = 2 E^3 (a1^2 - A^2)/(a1 A)^2, f4 = -6 E^4 (a1^3 + A^3)/(a1 A)^3. So f4/(8 f2^2) = -3 (a1^3 + A^3)/(4 a1 A E^2), 5 f3^2/(24 f2^3) = -5 (a1^2 - A^2)^2/(6 a1 A E^3), and
with a1^3 + A^3 = E (a1^2 - a1 A + A^2), (a1^2 - A^2)^2 = E^2 (a1 - A)^2: Abeta_1 = [-9 (a1^2 - a1 A + A^2) + 10 (a1 - A)^2]/(12 a1 A E) = (a1^2 - 11 a1 A + A^2)/(12 a1 A E).
On the right, [z^1](h_{(E,-a1,-A)} E/(E + z)) = (1/12)(1/a1 + 1/A - 1/E) - 1/E = (E^2 - 13 a1 A)/(12 a1 A E) = (a1^2 - 11 a1 A + A^2)/(12 a1 A E). On T_F both are
(x^2 - 11 x + 1)/(12 c x), [PR] s.3.8. [Evidence: check2_rows.py: r^pt_{0,1}(z; 0) = h_{(E,-a1,-A)} E/(E + z) through z^11 at two generic CY points, from the
Wick row.]
Remark 3.15. Rationality replaces the passage x > 0 -> all x of [PR]: here nothing new is proved analytically; the only analytic input is [PR]
Proposition 3.15 as it stands.

### 3.9 Regularity, the frame-B row, F2 [substitution; one new line: the constant]
**Corollary 3.16_CY (regularity).** For every k >= 0, Psi_k = [tau^k](F Ephi alg) and H_k = [tau^k](F Ephi) are rational functions of z with
coefficients in Kw, without pole at z = 0.
Proof. Verbatim [PR] Corollary 3.16: Psi_k is in Kw(z) (F_i in Kw(z), X_j in Kw[1/z], alg with Taylor coefficients in Kw); by (3.17)
Psi_k = [tau^k](N/Abeta), Abeta a unit of C[[z]] and N in C[[z]][[tau]], so no negative powers of z; H_k is a Kw-linear combination of Psi_0, ..., Psi_k. QED

**The constant line (new).** Multisets: nw + (E, -a1, -A) = (a1, A, -e1, -e2, E, -a1, -A) = (E, -e1, -e2) + (a1, -a1) + (A, -A), so
    h_nw(z) h_{(E, -a1, -A)}(z) = h_{(E, -e1, -e2)}(z) = h_w(z),   equivalently   h_w^-1 h_nw h_{(E, -a1, -A)} = 1.                             (3.9a)
This is where the class w = (E, -e1, -e2) is forced: it is the unique multiset class (up to pairs v, -v) for which the frame-B constant is the
geometric series E/(E + z) (h_v depends on v only through the odd power sums Sum_a v_a^(1-2k), k >= 1, of the reciprocals, and h_v = 1 iff all of them vanish
iff the multiset {1/v_a} equals its negative, by Newton's identities (the power sums of a finite multiset determine it, and the even ones of {1/v_a} and {-1/v_a}
agree); so h_{w'} = h_w iff h_{w' - w} = 1 iff w' and w agree up to pairs v, -v). The same scalar h_nw multiplies both idempotents ([PR] Lemma 13.3), so the frame-B rows of the regular idempotent carry the
non-trivial constant h_{(a1, A, -E)} r^pt_{0,2}(z; 0) (s.8.6). On T_F, h_nw = 1 and (3.9a) is h_w = h_{(c, -x, -1)}, [PR] Proposition 3.15.

**Theorem 3.17_CY (the frame-B singular row).** In C[[z]][[t]], and coefficientwise in z as analytic functions of t near 0,
    r_{0,1}(z; t) = (E/(E + z)) F(z, tau) Ephi(z, tau) alg(tau),   tau = 1 - e^-t;                                                          (3.18)
precisely, [t^k] r_{0,1}(z; t) = Sum_{j<=k} ([t^k] tau^j) (E/(E + z)) Psi_j(z), a rational function of z regular at z = 0.
Proof. r_{0,1} = h_w^-1 h_nw r^pt_{0,1} (C7'), (C8') = h_w^-1 h_nw N (Proposition 3.10_CY) = h_w^-1 h_nw Abeta F Ephi alg (3.17) = (E/(E + z)) F Ephi alg by
Proposition 3.15_CY and (3.9a). The coefficientwise meaning is as in [PR] Theorem 3.17 (Corollary 3.16_CY, Proposition 3.10_CY(iii), (T)). QED
[Evidence: check2_rows.py, KEY FORMULA block: h_{(a1,A,-E)} times the Wick row in y = log U equals (E/(E + z)) F Ephi alg through tau^6 z^9 and through
tau^8 z^11 at two generic CY points; the negative z-powers of the right side cancel.]

**Corollary 3.18_CY (F2).** r_{0,1}(z; 0) = E/(E + z). Also r^pt_{0,1}(z; 0) = Abeta = h_{(E, -a1, -A)} E/(E + z) and, for the R_Z row with constant,
h_nw r^pt_{0,1}(z; 0) = h_{(E, -e1, -e2)} E/(E + z).
Proof. tau = 0: F = 1, dphi = 0, alg = 1. QED

### 3.10 A purely formal proof of the regularity: eikonal and transport equations [substitution; re-derived]
Write F = Exp[Sph/z] Phi with Sph in tau O. With gamma = E/z + 2, alpha + beta + 1 = (a2 + A)/z + 3, alpha beta = a2 A/z^2 + (a2 + A)/z + 1, as operators on Kz[[tau]],
    Exp[-Sph/z] Gauss[Exp[Sph/z] Phi] = z^-2 Eik(Sph_tau) Phi + z^-1 (P_tr Phi_tau + R_tr Phi) + T_tr Phi,                                   (3.19)
    Eik(p) := tau (1 - tau) p^2 + (E - (a2 + A) tau) p - a2 A,   P_tr := 2 tau (1 - tau) Sph_tau + E - (a2 + A) tau,
    R_tr := tau (1 - tau) Sph_tautau + (2 - 3 tau) Sph_tau - (a2 + A),   T_tr := tau (1 - tau) d^2/dtau^2 + (2 - 3 tau) d/dtau - 1.
Derivation: F_tau = e^{Sph/z} (Phi_tau + Sph_tau Phi/z), F_tautau = e^{Sph/z} (Phi_tautau + 2 Sph_tau Phi_tau/z + Sph_tautau Phi/z + Sph_tau^2 Phi/z^2); collect the
powers of z: z^-2: tau (1 - tau) Sph_tau^2 + (E - (a2 + A) tau) Sph_tau - a2 A; z^-1: 2 tau (1 - tau) Sph_tau Phi_tau + tau (1 - tau) Sph_tautau Phi + (E - (a2 + A) tau) Phi_tau
+ (2 - 3 tau) Sph_tau Phi - (a2 + A) Phi; z^0: tau (1 - tau) Phi_tautau + (2 - 3 tau) Phi_tau - Phi.
Eikonal equation. Eik(p) = 0 is quadratic in p with discriminant (E - (a2 + A) tau)^2 + 4 a2 A tau (1 - tau) = E^2 - [2 E (a2 + A) - 4 a2 A] tau + (A - a2)^2 tau^2
= delta^2, because A - a2 = a0 and 2 E (a2 + A) - 4 a2 A = 2 E (E - a1 + a2) - 4 a2 (E - a1) = 2 E (E - a1 - a2) + 4 a1 a2 = E^2 + u. Its root regular at
tau = 0 has p(0) = a2 A/E, the other root is -E/tau + O(1) (exponent 1 - gamma). Sph := dphi solves it: with s = s* and p := dphi_tau = (A - B s)/(1 - tau)
(Lemma 3.1(b)),
    (1 - tau) Eik(p) = tau (A - B s)^2 + (E - (a2 + A) tau)(A - B s) - a2 A (1 - tau)
                     = B^2 tau s^2 + B tau s (a2 + A - 2 A) + A tau (A - a2 - A + a2) + A (E - a2) - E B s = B [B tau s^2 - a0 tau s - E s + A] = B N(s) = 0,
and p(0) = (A - B A/E) = A a2/E. Coefficients: P_tr = (d Eik/dp)(Sph_tau), so P_tr^2 = discriminant = delta^2; directly, P_tr = 2 tau (A - B s*) + E - (a2 + A) tau
= E + (A - a2) tau - 2 B tau s* = E + a0 tau - 2 B tau s* = delta (Lemma 3.1(b)(i)). For R_tr: differentiating (1 - tau) Sph_tau = A - B s* gives
tau (1 - tau) Sph_tautau = tau Sph_tau - B tau s*_tau, so
    R_tr = 2 (1 - tau) Sph_tau - B tau s*_tau - (a2 + A) = 2 (A - B s*) - a2 - A - B tau s*_tau = a0 - 2 B s* - B tau s* L_2/delta,   R_tr(0) = a0 - 2 B A/E - 0 = -(2 a1 a2 + a0 E)/E,
using s*_tau = s* L_2/delta (Lemma 3.1(b)(iii)) and a0 - 2 B A/E = (a0 E - 2 A B)/E = -(2 a1 a2 + a0 E)/E (A B = a1 a2 + a0 E). On T_F: R_tr = -x s* (2 + x tau s*/delta),
R_tr(0) = -2x/c, [PR] s.3.10. With Eik = 0, z Exp[-Sph/z] Gauss[Exp[Sph/z] Phi] = delta Phi_tau + R_tr Phi + z T_tr Phi.

**Proposition 3.19_CY.** For Phi = Sum_m z^m Phi_m, the transport equations delta (Phi_m)_tau + R_tr Phi_m + T_tr Phi_{m-1} = 0 (Phi_{-1} := 0) have unique
solutions analytic on D with Phi_0(0) = 1, Phi_m(0) = 0 (m >= 1): Phi_0 = Exp[-Int_0^tau R_tr/delta] = 1/alg, Phi_m = -Phi_0 Int_0^tau T_tr[Phi_{m-1}]/(delta Phi_0),
and in Kz[[tau]]
    Ephi F = Phi,   Phi in O(D)[[z]],   [z^0] Phi = 1/alg.                                                                               (3.20)
Proof. On D, delta has no zero, s* is analytic (s.1.3) and so is R_tr; this gives existence, uniqueness and analyticity. By (3.19) with Eik = 0,
Yt := Exp[dphi/z] iota(Phi) satisfies Gauss[Yt] = 0 and Yt(z, 0) = 1, so Yt = F by Lemma 3.13_CY, i.e. Ephi F = Phi. It remains to identify Phi_0 with 1/alg,
i.e. to show R_tr = delta (log alg)_tau; since both Phi_0 and 1/alg are analytic on D with value 1 at 0 (1/alg by (3.4')) and Phi_0 is the unique solution of delta Phi_0' + R_tr Phi_0 = 0, it
suffices to verify this identity. By (3.4), 2 (log alg)_tau = delta_tau/delta + U_1'/U_1 - s*_tau/s* + s*_tau/(1 - s*) (prime = d/dtau). From
delta = E + a0 tau - 2 B tau s*, delta_tau = a0 - 2 B s* - 2 B tau s*_tau = R_tr - B tau s*_tau, so the claim 2 R_tr/delta = 2 (log alg)_tau becomes, after multiplying
by delta and inserting s*_tau = s* L_2/delta, U_1' = -s* - tau s*_tau = -s* (delta + tau L_2)/delta,
    -s* (delta + tau L_2)/U_1 - L_2 + L_2 s*/(1 - s*) - delta_tau = 2 B tau s* L_2/delta
(the right side is 2 R_tr - 2 delta_tau = 2 B tau s* L_2/delta), i.e. (inserting delta_tau = a0 - 2 B s* - 2 B tau s* L_2/delta, whose last term cancels the right side)
    -s* (delta + tau L_2)/U_1 - L_2 + L_2 s*/(1 - s*) - a0 + 2 B s* = 0.
Since -L_2 - a0 + 2 B s* = B s*, dividing by s* this reads (delta + tau L_2)/U_1 = L_2/(1 - s*) + B = (L_2 + B (1 - s*))/(1 - s*) = (B - a0)/(1 - s*) = a1/(1 - s*),
i.e. (delta + tau L_2)(1 - s*) = a1 U_1, which is Lemma 3.1(b)(ii) (second form). QED
[Evidence: check1_identities.py: (log alg^2)_tau = 2 R_tr/delta at 25 rational points, both roots.]
Hence Psi_k = [tau^k](Phi alg) lies in C[[z]], Corollary 3.16_CY again, and the z-coefficients of Ephi F are analytic on the whole disc D, up to the
nearest zero of delta^2 (where s* and s~ coalesce). Comparison with (3.17): N = Abeta (Phi alg); the Laplace route fixes Abeta, the transport route the
tau-dependence. (The identification Phi_0 = 1/alg also follows formally from Proposition 3.12_CY, whose z^0 coefficient is [z^0](N/Abeta)/alg = 1/alg.)

### 3.11 First orders
F_1 = (a2 + z)(A + z)/(z (E + 2z)), X_1 = -dphi_1/z with dphi_1 := [tau^1] dphi = dphi_tau(0) = a2 A/E, so
    H_1 = F_1 + X_1 = [E (a2 + z)(A + z) - a2 A (E + 2z)]/(E z (E + 2z)) = [E z + E (a2 + A) - 2 a2 A]/(E (E + 2z)):
the poles at z = 0 cancel. [tau^1] alg = -(E^2 + u)/(2 E^2) (from (3.4): log alg^2 = log(delta/E) + log U_1 - log(E s*/A) - log(E (1 - s*)/a1) with first-order
terms -(E^2 + u)/(2 E^2), -A/E, -a1 a2/E^2, +A a2/E^2 in tau, using s*_tau(0) = s*(0) L_2(0)/E = A a1 a2/E^3; their sum is -(E^2 + u)/(2 E^2) - (A B + a1 a2)/E^2
= -(E^2 + u)/E^2, since A B + a1 a2 = 2 a1 a2 + a0 E = (E^2 + u)/2). So [tau^1] r_{0,1} = (E/(E + z)) (H_1 - (E^2 + u)/(2 E^2)) and [t^1] r_{0,1} = [tau^1] r_{0,1}.
On T_F: H_1 = (2x + c z)/(c (c + 2z)), [tau^1] alg = -2x/c^2, [PR] s.3.11.

Remark 3.20 (Euler's integral; motivation, used nowhere). Under Re(A/z) > -1, Re(a1/z) > -1, |tau| < 1, J(z, tau) := Int_0^1 s^(A/z) (1 - s)^(a1/z) (1 - tau s)^(-a2/z - 1) ds
= B(A/z + 1, a1/z + 1) F(z, tau) (Euler's integral with a = alpha, b = beta, c = gamma), and Gauss J = 0 by the exact identity of [PR] Remark 3.20 with
(a, b, g) = (alpha, beta, gamma). The formal expansion obeys Gauss's equation because of (3.16) alone.

## 4. Metric, eigenvalues, QDE and the symplectic constant at the singular idempotent (Steps 3, 5, 6)

### 4.1 Frobenius data [substitution; the generic data are those of [PR] Prop. 13.2(b)(i)]
**Lemma 4.1_CY.** In the flat basis (1, H):
 (a) eta = (1/(A B e1 e2)) [[E/(a1 a2), 1], [1, -a0]] and eta^-1 = e1 e2 [[a0 a1 a2, a1 a2], [a1 a2, -E]], independent of t;
 (b) 1 * a = a, and H * H = -(a1 a2/tau) 1 + (E/tau - a0) H;
 (c) Leq(L) = L^2 - (E/tau - a0) L + a1 a2/tau, Ssum = E/tau - a0, Pprod = a1 a2/tau;
 (d) the product is commutative, associative and Frobenius.
Proof. (a) The genus-0 classical potential is PrepCl = (1/6)[t0^3/(a0 a1 a2 e1 e2) - (t0 - a0 t1)^3/(a0 A B e1 e2)] (localisation at the two fixed points,
[PR] proof of Lemma 4.1(a), stated there at generic weights), with third derivatives PrepCl_000 = E/(a1 a2 A B e1 e2), PrepCl_001 = 1/(A B e1 e2),
PrepCl_011 = -a0/(A B e1 e2), PrepCl_111 = a0^2/(A B e1 e2); eta_ij = PrepCl_0ij gives eta. det eta = (-a0 E - a1 a2)/(a1 a2 (A B e1 e2)^2) = -1/(a1 a2 A B e1^2 e2^2),
since a0 E + a1 a2 = A B; inverting, eta^-1 = -a1 a2 A B e1^2 e2^2 [[-a0, -1], [-1, E/(a1 a2)]]/(A B e1 e2) = e1 e2 [[a0 a1 a2, a1 a2], [a1 a2, -E]].
(b) Prep = PrepCl + Li_3(e^{t1})/(e1 e2), Prep_111 = PrepCl_111 - 1/(e1 e2 tau) (as in [PR]: (d/dt1)^3 Li_3(e^{t1}) = Q/(1 - Q) = -1/tau). Structure constants
cch[i, 1, 1] = Sum_m eta^{im} Prep_m11: cch[0, 1, 1] = e1 e2 [a0 a1 a2 (-a0/(A B e1 e2)) + a1 a2 (a0^2/(A B e1 e2) - 1/(e1 e2 tau))] = -a1 a2/tau;
cch[1, 1, 1] = e1 e2 [a1 a2 (-a0/(A B e1 e2)) - E (a0^2/(A B e1 e2) - 1/(e1 e2 tau))] = -a0 (a1 a2 + a0 E)/(A B) + E/tau = E/tau - a0. (This agrees with
[PR] Prop. 13.2(b)(i): L_+ + L_- = -(a0 + (a1 + a2) Q)/(1 - Q) = E/tau - a0 and L_+ L_- = -a1 a2 Q/(1 - Q) = a1 a2/tau, with Q/(1 - Q) = -1/tau,
1/(1 - Q) = -(1 - tau)/tau.) (c), (d) as in [PR]. QED
On T_F: eta^-1 = {{0, x^2}, {x^2, -c x}}, H * H = -(x/tau) 1 + (c/tau) H. [Evidence: check1_identities.py (4).]

### 4.2 Critical points, eigenvalues, canonical coordinates [substitution]
U W_U = -a1 Q U/(1 - Q U) - a2 - A U/(1 - U), W_2 = -a1 Q U/(1 - Q U)^2 - A U/(1 - U)^2, (partial_t W)(U) = -a1 Q U/(1 - Q U).
**Lemma 4.2_CY.** (a) U W_U = 0 iff LeqU = B Q U^2 - (a0 + Q (a1 - a2)) U - a2 = 0 ((C2')).
 (b) At a critical point (partial_t W)(U) = a2 + A U/(1 - U) = L(U).
 (c) With U = 1 - tau s (tau != 0), U is critical iff N(s) = 0; the critical points are U_1 = 1 - tau s*, U_2 = 1 - tau s~.
 (d) L_1 = A/(tau s*) - a0, L_2 = B s* - a0 = A/(tau s~) - a0; L_1 + L_2 = E/tau - a0, L_1 L_2 = a1 a2/tau, L_1 - L_2 = delta/tau. L_1, L_2 are the roots of
     Leq, the eigenvalues of H*. L_1 = E/t + l(t) with l := A/(tau s*) - a0 - E/t in O, l_0 = E/2 - a0 - a1 a2/E; L_2 is in O with L_2(0) = a1 a2/E.
 (e) u_i := W(U_i) satisfies du_i/dt = L_i. Explicitly u_1 = a1 t + E log tau + phi_tau(s*) + a1 log(-1), so with mu := E/tau - A:
     L_1 = mu + dphi', L_1 - dphi' = E/tau - A; and L_2 = mu + (a2 - L_1), since d/dt phi_tau(s~) = A - B s~ = a2 - L_1.
Proof. (a) is (C2'). (b) U W_U = 0 gives -a1 Q U/(1 - Q U) = a2 + A U/(1 - U). (c) By (3.1), W(1 - tau s) = phi_tau(s) + const for tau != 0, and d/ds = -(tau/U) theta
(Lemma 3.1(d)) with tau/U != 0 near the critical points, so U = 1 - tau s is critical iff s is a critical point of phi_tau, i.e. N(s) = 0 (Lemma 3.1(a));
equivalently LeqU(1 - tau s) = Q tau N(s) as polynomials in s (compare the s^2 coefficients B Q tau^2 and B tau, and the roots). (d) By (b) and 1 - U_i = tau s_i: L_i = a2 + A (1 - tau s_i)/(tau s_i) =
A/(tau s_i) - a0. With s* s~ = A/(B tau), A/(tau s~) = B s*. Sum: A (s* + s~)/(tau s* s~) - 2 a0 = (E + a0 tau)/tau - 2 a0 = E/tau - a0. Product:
(A/(tau s*) - a0)(B s* - a0) = A B/tau - a0 (L_1 + a0) - a0 (L_2 + a0) + a0^2 = (A B - a0 E)/tau = a1 a2/tau. Difference: A/(tau s*) - B s* =
(A - B tau s*^2)/(tau s*) = (E + a0 tau - 2 B tau s*)/tau = delta/tau, by N(s*) = 0 and Lemma 3.1(b)(i). By Lemma 4.1(b) the minimal polynomial of H* is
Leq. Expansion: 1/tau = 1/t + 1/2 + O(t), s* = A/E + s*_1 tau + O(tau^2) with s*_1 = s*_tau(0) = s*(0) L_2(0)/delta(0) = (A/E)(a1 a2/E)/E = A a1 a2/E^3
(Lemma 3.1(b)(iii), L_2(0) = B A/E - a0 = (A B - a0 E)/E = a1 a2/E); so A/(tau s*) = (E/tau)(1 - (E/A) s*_1 tau + O(tau^2)) = E/t + E/2 - a1 a2/E + O(t), and
l_0 = E/2 - a0 - a1 a2/E. (e) u_i' = (partial_t W)(U_i) = L_i by (b). The explicit form is (3.1) with log Q = t, so u_1' = a1 + E (1 - tau)/tau + d/dt phi_tau(s*), and
d/dt phi_tau(s*) = dphi' = A - B s* (Lemma 3.1(b)); since a1 + E (1 - tau)/tau = E/tau - A + (a1 + A - E) = mu, this gives L_1 = mu + dphi'. Directly:
mu + dphi' = E/tau - B s* = A/(tau s*) - a0 = L_1, because N(s*) = 0 divided by tau s* reads B s* - (E + a0 tau)/tau + A/(tau s*) = 0. At s~, Lemma 3.1(b)(ii) at the second root gives d/dt phi_tau(s~) = (1 - tau) a2 s~/(1 - tau s~) = A - B s~, and B s~ = A/(tau s*) = L_1 + a0, so
A - B s~ = a2 - L_1 and mu + a2 - L_1 = E/tau - a0 - L_1 = L_2. QED
On T_F: L_1 = 1/(tau s*) = c/tau - x s*, L_2 = x s*, l_0 = c/2 - x/c, mu = c/tau - 1, [PR] Lemma 4.2.

### 4.3 Idempotents and the metric: (V1), (V2), F1 [substitution]
**Proposition 4.3_CY.** (a) e_1 = (H - L_2)/(L_1 - L_2), e_2 = (H - L_1)/(L_2 - L_1) are the idempotents (verbatim [PR] 4.3(a)).
 (b) Delta_i = 1/eta(e_i, e_i) (verbatim).
 (c) For {i, j} = {1, 2}:
        Delta_i = a1 a2 A B e1 e2 (L_i - L_j)/(a1 a2 - E L_j) = e1 e2 (L_i - L_j)(a1 a2 - E L_i) = e1 e2 (a0 a1 a2 + 2 a1 a2 L_i - E L_i^2) = a1 A e1 e2 W_2(U_i),   (4.1)
     identities in QA (V1 at both roots). This is [PR] (13.2).
 (d) Delta_1 = -a1 A e1 e2 delta U_1/(tau^2 s* (1 - s*)) = -e1 e2 E^3 alg^2/tau^2;  Delta_2 = e1 e2 B delta s* L_2;  a1 a2 - E L_2 = -B tau s* L_2;
     (a1 a2 - E L_1)(a1 a2 - E L_2) = a1 a2 A B;  Delta_1 Delta_2 = -e1^2 e2^2 a1 a2 A B delta^2/tau^2 (from the second form of (4.1) and L_1 - L_2 = delta/tau);  1/Delta_1 + 1/Delta_2 = E/(a1 a2 A B e1 e2) = eta(1, 1).
 (e) F1: Delta_2 is in O with Delta_2(0) = e1 e2 a1 a2 A B/E != 0 and Delta_2'(0) = 0; Delta_1 is tau^-2 times a unit of O, Delta_1 = -e1 e2 E^3 t^-2 (1 + O(t)).
Proof. (c) First form: eta(e_i, e_i) = eta(e_i, 1) = (eta(H, 1) - L_j eta(1, 1))/(L_i - L_j) = [1/(A B e1 e2) - L_j E/(a1 a2 A B e1 e2)]/(L_i - L_j) =
(a1 a2 - E L_j)/(a1 a2 A B e1 e2 (L_i - L_j)). Second form: (a1 a2 - E L_1)(a1 a2 - E L_2) = a1^2 a2^2 - E a1 a2 (E/tau - a0) + E^2 a1 a2/tau = a1 a2 (a1 a2 + a0 E) = a1 a2 A B,
so a1 a2 A B/(a1 a2 - E L_j) = a1 a2 - E L_i. Third form: (L_i - L_j)(a1 a2 - E L_i) = a1 a2 L_i - E L_i^2 - a1 a2 L_j + E L_i L_j, and -a1 a2 L_j = a1 a2 L_i - a1 a2 (E/tau - a0),
E L_i L_j = E a1 a2/tau; the tau^-1 terms cancel and a0 a1 a2 + 2 a1 a2 L_i - E L_i^2 remains. Fourth form: at a critical point L = -a1 Q U/(1 - Q U), so
1/(1 - Q U) = 1 + Q U/(1 - Q U) = (a1 - L)/a1 and -a1 Q U/(1 - Q U)^2 = L (a1 - L)/a1; and U/(1 - U) = (L - a2)/A (Lemma 4.2(b)), 1/(1 - U) = (L + a0)/A, so
A U/(1 - U)^2 = (L - a2)(L + a0)/A. Hence a1 A W_2(U) = A L (a1 - L) - a1 (L - a2)(L + a0) = -(A + a1) L^2 + a1 (A - a0 + a2) L + a0 a1 a2 = -E L^2 + 2 a1 a2 L + a0 a1 a2.
The code's literal evaluation of the quadratic form and the identification in QA by the two characters ev_1, ev_2 are as in [PR] (c).
(d) The first form of Delta_1 is (4.1) with (3.5); the second is (3.4). Delta_2 = e1 e2 (L_2 - L_1)(a1 a2 - E L_2) with L_2 - L_1 = -delta/tau and
a1 a2 - E L_2 = a1 a2 - E (B s* - a0) = A B - E B s* = B (A - E s*) = -B tau s* L_2 (Lemma 3.1(b)(iv)) gives Delta_2 = e1 e2 B delta s* L_2. The rest is (c).
(e) Delta_2 = e1 e2 B delta s* L_2 is a product of elements of O, with value e1 e2 B E (A/E)(a1 a2/E) = e1 e2 a1 a2 A B/E at 0. Derivative at 0: with
delta_tau(0) = -(E^2 + u)/(2E), s*_tau(0) = A a1 a2/E^3 =: s_1, (L_2)_tau(0) = B s_1,
    (delta s* L_2)_tau(0) = -(E^2 + u)/(2E) (A/E)(a1 a2/E) + E s_1 (a1 a2/E) + E (A/E) B s_1 = (A a1 a2/E^3) [ -(E^2 + u)/2 + a1 a2 + A B ] = 0,
because A B + a1 a2 = 2 a1 a2 + a0 E = (E^2 + u)/2. So Delta_2'(0) = 0 (this was an observation in the feasibility reports; it is proved here). Delta_1 is
Proposition 4.4. QED
On T_F: Delta_1 = -x delta/(tau^2 s*^2), Delta_2 = x^3 delta s*^2, [PR] Prop. 4.3(d). [Evidence: check1_identities.py (3), (6).]

### 4.4 Step 3: the Hessian in the Euler coordinate; rho; the KEY FORMULA [substitution]
**Proposition 4.4_CY.** (a) Delta_1 = a1 A e1 e2 W_2(U_1) = a1 A e1 e2 U_1^2 phi_tau''(s*)/tau^2 = -e1 e2 E^3 alg^2/tau^2 (by (4.1), (3.5), (3.4)).
 (b) rho = tau/(t alg), a unit of O with rho(0) = 1, and t^2 rho^2 Delta_1 = -e1 e2 E^3 exactly (4.2); Delta_1 t^2/(-e1 e2 E^3) = (t/tau)^2 alg^2 is a unit with
     value 1 (F1 for Delta_1).
 (c) KEY FORMULA. G(z; t) := rho r_{0,1}(z; t) = (E/(E + z)) (tau/t) F(z, tau) Ephi(z, tau) (4.3); G is in C[[z, t]], G(z; 0) = E/(E + z), and every G_k is a
     rational function of z without pole at z = 0 whose Taylor series is [t^k](rho r_{0,1}).
 (d) rho = 1 + u t/(2 E^2) + O(t^2), Delta_1 = -(e1 e2 E^3/t^2)(1 - u t/E^2 + O(t^2)).
Proof. (a), (b), (c) verbatim [PR] Proposition 4.4 with (3.18) and Corollary 3.16_CY. (d) alg = 1 - ((E^2 + u)/(2 E^2)) tau + O(tau^2) (s.3.11), tau/t =
1 - t/2 + O(t^2), so rho = (tau/t)/alg = 1 + [(E^2 + u)/(2 E^2) - 1/2] t + O(t^2) = 1 + u t/(2 E^2) + O(t^2); rho^-2 = 1 - u t/E^2 + O(t^2). QED
On T_F: rho = 1 - (1 - x)^2 t/(2c^2) + O(t^2), [PR] 4.4(d) (u = -(1 - x)^2). [Evidence: check2_rows.py: rho_1 = u/(2E^2); KEY FORMULA as in s.3.9.]

### 4.5 The operator identity (V7) [substitution; re-derived]
**Lemma 4.5_CY.** In the ring of differential operators generated by Bop := C[z, 1/z][tau, 1/tau, 1/(1 - tau)] and one symbol Dd with Dd a = a Dd + z a' (so
Dd(tau) = z (1 - tau)), with d/dtau := (z (1 - tau))^-1 Dd, omega = (E + z)/tau - (A + z) and Gauss with the parameters of s.1.5,
    (Dd + omega)^2 - (E/tau - a0)(Dd + omega) + a1 a2/tau = ((1 - tau) z^2/tau) Gauss.                                                     (4.4)
The same holds for any Dd' with the same Leibniz rule. (V7) is the case Dd = Dz conjugated by Pf := (1 - tau)^(-a1/z) tau^(E/z + 1): Dz(Pf f) = Pf (Dz + omega) f,
since z (log Pf)' = a1 + E (1 - tau)/tau + z (1 - tau)/tau = mu + z (1 - tau)/tau = omega.
Proof. Write S := E/tau - a0. For f in Rz, (Dd + omega)^2 f - S (Dd + omega) f + (a1 a2/tau) f = Dd^2 f + (2 omega - S) Dd f + [Dd(omega) + omega (omega - S) + a1 a2/tau] f,
with Dd f = z (1 - tau) f_tau, Dd^2 f = z^2 (1 - tau)^2 f_tautau - z^2 (1 - tau) f_tau.
 * Coefficient of f. Dd(omega) = -z (1 - tau)(E + z)/tau^2. omega - S = z/tau - (a2 + z) (since (E + z) - E = z and -(A + z) + a0 = -(a2 + z)), so
   omega (omega - S) = z (E + z)/tau^2 - [(E + z)(a2 + z) + z (A + z)]/tau + (A + z)(a2 + z). The tau^-2 terms add up to z (E + z)/tau, and
   z (E + z) - (E + z)(a2 + z) - z (A + z) + a1 a2 = -[a2 E + (a2 + A) z + z^2 - a1 a2] = -[a2 A + (a2 + A) z + z^2] = -(A + z)(a2 + z) (a2 E - a1 a2 = a2 A).
   So the coefficient is -(A + z)(a2 + z)/tau + (A + z)(a2 + z) = -(1 - tau)(A + z)(a2 + z)/tau.
 * Coefficient of f_tau: -z^2 (1 - tau) + z (1 - tau)(2 omega - S) = z (1 - tau) [2 omega - S - z] = z (1 - tau) [(E + 2z)/tau - (2A - a0 + 3z)] = z (1 - tau) [(E + 2z)/tau - (a2 + A + 3z)].
 * Coefficient of f_tautau: z^2 (1 - tau)^2.
On the right, ((1 - tau) z^2/tau) Gauss[f] has f_tautau-coefficient z^2 (1 - tau)^2, f_tau-coefficient (1 - tau) z^2 [(2 + E/z)/tau - (3 + (a2 + A)/z)] =
z (1 - tau) [(E + 2z)/tau - (a2 + A + 3z)], and f-coefficient -((1 - tau) z^2/tau)(1 + a2/z)(1 + A/z) = -(1 - tau)(a2 + z)(A + z)/tau. QED
On T_F: a2 = A = 1, [PR] Lemma 4.5. [Evidence: check1_identities.py (2): (4.4) on tau^n, -3 <= n <= 6, at 12 random points (weights, z).]

### 4.6 Step 5: the second row of the QDE at the singular idempotent [substitution]
**Proposition 4.6_CY.** In Rz0 (r_{0,1} in C[[z, tau]], r_{1,1} in tau^-1 C[[z, tau]]), hence to all orders in z as identities of meromorphic functions of t:
    Row 1: r_{1,1} = (L_1 + z D_1) r_{0,1};   Row 2: (L_1 + z D_1) r_{1,1} = -(a1 a2/tau) r_{0,1} + (E/tau - a0) r_{1,1},  equivalently  z D_1 r_{1,1} = L_2 (r_{1,1} - L_1 r_{0,1}).
In matrix form z D_1 v_1 = (A - L_1) v_1, and for vn := t rho v_1 in C[[z, t]]^2, z vn' = (A - L_1) vn (4.5).
Proof. As [PR] Proposition 4.6. By (4.2), L_1 + z D_1 = (t rho)^-1 (z d/dt + L_1)(t rho) (4.6). (i) Conjugation: (z d/dt + L_1)(tau Ephi f) = tau Ephi [z f' +
(L_1 - dphi' + z (1 - tau)/tau) f] = tau Ephi (Dz + omega) f, because L_1 - dphi' = E/tau - A = mu (Lemma 4.2(e)) and mu + z (1 - tau)/tau = omega. (ii) By (4.3),
t rho r_{0,1} = (E/(E + z)) tau Ephi F, so t rho r_{1,1} = (E/(E + z)) tau Ephi (Dz + omega) F and t rho (L_1 + z D_1) r_{1,1} = (E/(E + z)) tau Ephi (Dz + omega)^2 F.
(iii) By (4.4) and Gauss[F] = 0, t rho [(L_1 + z D_1) r_{1,1} + (a1 a2/tau) r_{0,1} - (E/tau - a0) r_{1,1}] = (E/(E + z)) Ephi (1 - tau) z^2 Gauss[F] = 0. The second
form uses E/tau - a0 - L_1 = L_2 and a1 a2/tau = L_1 L_2. Membership as in [PR]: t L_1 is in O. QED
From (ii): r_{1,1} = r_{0,1} [(E + z)/tau - (A + z) + z (1 - tau) F_tau/F], equivalently rho r_{1,1} = (L_1 + z/t) G + z G' (4.7). In F-form, with
PsiF1 := tau (Dz + omega) F = z tau (1 - tau) F_tau + (E + z - (A + z) tau) F, Row 2 reads z (1 - tau) (PsiF1)_tau = a2 PsiF1 - a1 a2 F (the coefficient a2 is -[omega - z(1 - tau)/tau - (E/tau - a0)] = A - a0 = L_2 + dphi'; corrected 2026-11-03, the earlier text had coefficient 1).
[Evidence: check2_rows.py: Row 2 at U_1 from the Wick-built rows, z^m <= 9, t^-3..t^2, at two CY points.]

### 4.7 The operator A [substitution]
**Proposition 4.7_CY.** (a) A is the transpose of the matrix M_H = [[0, -a1 a2/tau], [1, E/tau - a0]] of H* in the basis (1, H); A acts on covectors by s -> s(H * .).
 (b) chi_i = (1, L_i) satisfies A chi_i = L_i chi_i.
 (c) A is self-adjoint for < , >: the bilinear form <A ., .> has the symmetric matrix e1 e2 [[-(a1 a2)^2/tau, a1 a2 E/tau], [a1 a2 E/tau, a1 a2 + a0 E - E^2/tau]].
 (d) <chi_i, chi_j> = delta_ij Delta_i.
 (e) A - L_1 = (L_2 - L_1) Pr_2, Pr_2 a := (pi(a)/Delta_2) chi_2, with
        pi(a) := <a, chi_2> = e1 e2 [ a1 a2 (a0 + L_2) a_0 + (a1 a2 - E L_2) a_1 ] = e1 e2 B s* ( a1 a2 a_0 - tau L_2 a_1 ),
     and, for all covectors a, b,
        <(A - L_1) a, b> = ((L_2 - L_1)/Delta_2) pi(a) pi(b) = pi(a) pi(b)/Theta,   Theta := Delta_2/(L_2 - L_1) = e1 e2 (a1 a2 - E L_2) = -e1 e2 B tau s* L_2.      (4.8)
Proof. (a) H * 1 = H, H * H = -(a1 a2/tau) 1 + (E/tau - a0) H (Lemma 4.1(b)). (b) A (1, L) = (L, -a1 a2/tau + (E/tau - a0) L) = L (1, L) iff Leq(L) = 0.
(c) With M := eta^-1 = e1 e2 [[p, q], [q, r]], p = a0 a1 a2, q = a1 a2, r = -E, and (A a)_{a'} = Sum_c A_{a' c} a_c, <A a, b> = Sum_{a', b'} M^{a' b'} (A a)_{a'} b_{b'} = a^T (A^T M) b,
so the form (a, b) -> <A a, b> has matrix A^T M = e1 e2 [[-a1 a2 q/tau, -a1 a2 r/tau], [p + (E/tau - a0) q, q + (E/tau - a0) r]]; it is symmetric iff -a1 a2 r/tau = p + (E/tau - a0) q, i.e. a1 a2 E/tau = a0 a1 a2 + a1 a2 E/tau - a0 a1 a2, true. (Conceptually
H* is eta-self-adjoint.) (d) <chi_i, chi_j> = e1 e2 [a0 a1 a2 + a1 a2 (L_i + L_j) - E L_i L_j]; for i != j this is e1 e2 [a0 a1 a2 + a1 a2 (E/tau - a0) - E a1 a2/tau] = 0; for
i = j it is e1 e2 [a0 a1 a2 + 2 a1 a2 L_i - E L_i^2] = Delta_i by (4.1). (e) chi_1, chi_2 is an orthogonal basis with non-zero norms, A = L_1 Pr_1 + L_2 Pr_2. For pi:
a0 a1 a2 + a1 a2 L_2 = a1 a2 (a0 + L_2) = a1 a2 B s*, and a1 a2 - E L_2 = -B tau s* L_2 (Proposition 4.3(d)). (L_2 - L_1)/Delta_2 = 1/Theta by definition of Theta, and
Delta_2 = e1 e2 (L_2 - L_1)(a1 a2 - E L_2) gives the two other forms of Theta. QED
pi(chi_1) = e1 e2 B s* (a1 a2 - tau L_1 L_2) = 0 is literally L_1 L_2 = a1 a2/tau. On T_F: pi(a) = x^3 s* (a_0 - tau s* a_1), Theta = -x^3 tau s*^2, [PR] 4.7.
[Evidence: check1_identities.py (3): self-adjointness and (4.8) on random covectors.]

### 4.8 Step 6: the symplectic constant K_1 = 1 [substitution]
**Proposition 4.8_CY.** (a) K_1(z) := <v_1(z), v_1(-z)>/Delta_1 = -(1/(e1 e2 E^3)) <vn(z), vn(-z)> lies in C[[z, t]] and is t-independent.
 (b) sigma(psi) := rho v_1(-psi) = (g, (L_1 - psi/t) g - psi g'), t sigma(psi) in C[[psi, t]]^2, sigma^(-1) = (0, E), g(psi)|_{t=0} = E/(E - psi) = E lam_1(psi) (F2).
 (c) K_1 = 1, i.e. <v_1(w), v_1(-w)> = Delta_1 to all orders in w.
 (d) <sigma(psi), sigma(-psi)> = rho^2 Delta_1 = -e1 e2 E^3/t^2; rho^2 E_11(psi, psi') = [-e1 e2 E^3/t^2 - <sigma(psi), sigma(psi')>]/(psi + psi') (4.9); and
     d sigma/dt = -sigma/t - (1/psi)(A - L_1) sigma.
Proof. Verbatim [PR] Proposition 4.8 with the substituted constants: (a) from (4.5), refl and self-adjointness (Proposition 4.7(c)); (b) rho r_{1,1} =
(L_1 + z/t) G + z G', t L_1 in O with value E (Lemma 4.2(d)), so t rho r_{1,1} = (t L_1 + z) G + z t G' has value (E + z) E/(E + z) = E at t = 0: the z-dependence
cancels exactly because of the geometric form E/(E + z) of F2. (c) K_1(z) = -(1/(e1 e2 E^3)) <t sigma(-z), t sigma(z)> has value -(1/(e1 e2 E^3)) <(0, E), (0, E)> =
-(1/(e1 e2 E^3)) e1 e2 (-E) E^2 = 1 at t = 0, and is t-independent by (a). (d) as in [PR]; the division by psi is legitimate since (A - L_1) chi_1 = 0. QED
[Evidence: check2_rows.py: K_1 = 1 for w^m <= 9 (11), t^0..t^2 (t^4), from the Wick-built rows.]

**Remark 4.9 (K_1 = 1 as a Wronskian identity).** vn = (E/(E + z)) tau Ephi PhiF, PhiF := (F, (Dz + omega) F), so K_1 = -tau^2 <PhiF(z), PhiF(-z)>/(e1 e2 E (E^2 - z^2)),
and K_1 = 1 iff tau^2 <PhiF(z), PhiF(-z)> = -e1 e2 E (E^2 - z^2). This follows from Gauss[F] = 0 alone: (omega(-z) - omega(z))/z = -2 (1 - tau)/tau, so
d/dt (tau^2 <PhiF(z), PhiF(-z)>) = 0, and at tau = 0, tau PhiF -> (0, E + z), giving e1 e2 (-E)(E + z)(E - z). (Verbatim [PR] Remark 4.9 with c -> E, -x c -> -e1 e2 E.)

**Remark 4.10 (flat sections; constants kappa).** Verbatim [PR] Remark 4.10 with Q^{x/z} -> Q^{a1/z}, tau^{c/z} -> tau^{E/z}, mu = E/tau - A, kappa_0^2 = -e1 e2 E^3,
kappa_1 := E/(kappa_0 (E + z)); the indicial exponents of nabla^2 - (E/tau - a0) nabla + a1 a2/tau on Rz at tau = 0 are 0 and 1 + E/z (the a0 and a1 a2 terms are
of lower tau-order). Not used.

## 5. Step 4: the dilaton leaves (L1) [verbatim after substitution; every step re-checked for hidden uses of a2 = A]

Sections 5-7 are written over Kw (s.1.4 (P)); G and z G are compared with series in z or psi only through Proposition 4.4(c).

### 5.1 The tau-coefficients
**Lemma 5.1_CY.** Let n, m >= 0.
 (a) F_n has poles only at z = 0 (order exactly n) and at z = -E/j, 2 <= j <= n + 1 (order <= 1), and F_n(oo) = 1/(n + 1).
 (b) X_m = [tau^m] Ephi is a polynomial in 1/z of degree <= m with X_m(oo) = delta_m0.
 (c) H_n = Sum_{n1+n2=n} F_{n1} X_{n2} has poles only at z = 0 (order <= n) and at z = -E/j, 2 <= j <= n + 1 (order <= 1), and H_n(oo) = 1/(n + 1).
Proof. (a) The numerator prod_{j=1}^{n} (a2 + j z)(A + j z) has degree 2n, value (a2 A)^n != 0 at z = 0 (so the pole at 0 has order exactly n) and leading
coefficient (n!)^2. The denominator n! z^n prod_{j=2}^{n+1} (E + j z) has degree 2n, leading coefficient n! (n + 1)!, a zero of order n at 0 and simple zeros at
the n distinct non-zero points -E/j (E != 0). A zero of the numerator (at -a2/i or -A/i) can only lower a pole order; the poles at -E/j are therefore of
order <= 1 (exactly 1 over Kw, and at every weight point with E not in {j a2/i, j A/i : 1 <= i <= n, 2 <= j <= n + 1}). F_n(oo) = (n!)^2/(n! (n + 1)!) = 1/(n + 1),
independent of the weights. (b), (c) verbatim [PR] Lemma 5.1: dphi in tau Kw[[tau]], X_m = Sum_{l<=m} ((-1)^l/l!) z^-l [tau^m](dphi^l); F_{n1} X_{n2} -> delta_{n2,0}/(n1 + 1) at oo. QED
The difference from T_F: there the numerator was prod (1 + j z)^2, with double zeros at -1/j; here it has simple zeros at -a2/j and -A/j. Nothing in the
argument uses the multiplicity or the location of these zeros, only the degree and the leading coefficient.

### 5.2 The t-coefficients of z G
**Lemma 5.2_CY.** For k >= 0 put R_k(z) := [t^k](z G) = z G_k(z). Then
    R_k = (E z/(E + z)) Sum_{n=0}^{k} d_{k,n} H_n,   d_{k,n} := [t^k]((tau/t) tau^n) in Q,   d_{k,k} = 1.                                           (5.1)
So R_k is a rational function with poles only at z = 0 and at z = -E/j, 1 <= j <= k + 1, of order <= 1 at each -E/j; it is bounded at oo and R_k(oo) = E delta_k0 (5.2).
Proof. Verbatim [PR] Lemma 5.2 with c -> E: (T) gives (5.1); the factor E z/(E + z) adds the simple pole j = 1; R_k(oo) = E Sum_{n<=k} d_{k,n}/(n + 1) =
E [t^k]((tau/t) Sum_n tau^n/(n + 1)) = E [t^k](-log(1 - tau)/t) = E delta_k0. QED
The limit is coefficientwise in tau (F(oo, tau) := Sum_n tau^n/(n + 1) = -log(1 - tau)/tau for every weight point, because F_n(oo) = 1/(n + 1) does not depend on
the weights); for k >= 1 it vanishes because (tau/t) F(oo, tau) = 1: the normalisation by rho removes lam_0 at generic weights for the same reason as on T_F.
Equivalently lim_{z->oo} z r_{0,1}(z; t) = E/rho(t), coefficientwise in t.

### 5.3 No pole at z = 0
By Proposition 4.4(c), G_k is regular at z = 0, so R_k = z G_k is regular at 0 with R_k(0) = 0. In (5.1) the separate terms have poles at 0 (z H_n up to order
n - 1); their cancellation is the regularity theorem of section 3, for either of two reasons: (a) Corollary 3.16_CY, which rests on (3.17) and the fact
that Abeta is a unit of C[[z]] (its value is not used); (b) Proposition 3.19_CY, purely formal. By hand: H_1 = [E z + E (a2 + A) - 2 a2 A]/(E (E + 2z)) (s.3.11).
Nothing here uses a2 = A.

### 5.4 Partial fractions
**Proposition 5.3_CY (L1 in the variable z).** R_0 = E z/(E + z) = E - E^2/(E + z). For k >= 1 there are a_j^(k) in Kw, 1 <= j <= k + 1, with
    R_k(z) = Sum_{j=1}^{k+1} a_j^(k)/(E + j z),   a_j^(k) = j Res_{z=-E/j} R_k,   Sum_{j=1}^{k+1} a_j^(k) = 0.                                          (5.3)
Proof. Verbatim [PR] Proposition 5.3: for k >= 1 the poles of R_k are simple and lie in {-E/j : j <= k + 1} (Lemma 5.2_CY, s.5.3), R_k(oo) = 0 (5.2), so (PF)
gives (5.3) without constant term (no polynomial part: a rational function with simple poles that tends to 0 at oo is the sum of its principal parts);
evaluating at z = 0, 0 = R_k(0) = (1/E) Sum_j a_j^(k) (E != 0). QED
With a_1^(0) := -E^2 and a_j^(k) := 0 for j > k + 1, the sum rule reads Sum_j a_j^(k) = -E^2 delta_k0.

### 5.5 Translation to psi: Theorem L1
**Theorem 5.4_CY (L1, all k).** q^(0) = E - E^2 lam_1, and for k >= 1, q^(k) = Sum_{j=1}^{k+1} a_j^(k) lam_j in Lam_{k+1}. Moreover g_0 = E lam_1 and for k >= 1
    g_k = Sum_{j=1}^{k+1} g_j^(k) lam_j = -(1/E) Sum_j j a_j^(k) lam_j in Lam_{k+1},   g_k = O(psi^-2) at psi = oo.                                        (5.4)
Proof. Verbatim [PR] Theorem 5.4: rho q_t(psi) = -psi g(psi; t), [t^k](rho q_t) = R_k(-psi), 1/(E + j z) = lam_j at z = -psi; (5.4) from G_k = R_k/z,
1/(z (E + j z)) = (1/E)(1/z - j/(E + j z)) and the sum rule; g_k = O(psi^-2) since R_k = O(1/z). QED
Remarks (as [PR]). (i) g_j^(k) = -(j/E) a_j^(k). (ii) L1 at order k >= 1 is equivalent to "g_k in Lam_{k+1} and g_k = O(psi^-2)" (by (LAM2) with E). (iii) The
three ingredients are the pole structure (Lemmas 5.1_CY, 5.2_CY), the regularity theorem (s.5.3) and the coefficientwise limit at z = oo. That every pole
-E/j, j <= k + 1, actually occurs is proved for the top pole when P_k != 0 (Corollary 11.3). (iv) F3 for Dil_1: Dil_1(d) = (-1)^d [z^(d-1)](rho^-1 G) is in Kw[[t]],
with value -E^(1-d) at t = 0, i.e. q_t|_{t=0} = q_0 = -E psi/(E - psi).
[Evidence: check2_rows.py: [t^k](z G) prod_{j<=k+1}(E + j z) has z-degree exactly k for 1 <= k <= 6 (through z^9) and k <= 8 (through z^11), from the
Wick-built row; the critic's leaf_high.py reports the same for k <= 10 at another point. The A_j closed form reproduces [t^k](z G) for k <= 5 (7).]

## 6. Step 7: the H-row in leaf form [verbatim after substitution]

### 6.1 The H-row
By (4.6), (4.7): rho r_{1,1}(z) = (L_1 + z/t) G + z dG/dt, sigma(psi) = rho v_1(-psi) = (g, (L_1 - psi/t) g - psi dg/dt) (6.1). Both components of t sigma lie in
Kw[[t]][[psi]] since t L_1 is in Kw[[t]]. Only F1 enters (6.1): t^2 rho^2 Delta_1 = -e1 e2 E^3 by definition of rho. By Lemma 4.2(d),
    L_1 = E/t + l(t),   l = A/(tau s*) - a0 - E/t in Kw[[t]],   l_0 = E/2 - a0 - a1 a2/E.                                                        (6.2)
The structural fact used in s.6.2 is that the coefficient of the pole of L_1, t L_1(0) = E, equals the leaf scale E of lam_j = 1/(E - j psi) and the
first weight of w. It holds because L_1 = A/(tau s*) - a0 with s*(0) = A/E.

### 6.2 The leaf structure of sigma
**Proposition 6.1_CY.** (a) sigma^(-1) = (0, E). (b) For k >= 0, [t^k] sigma_0 = g_k and
    [t^k] sigma_1 = (E - (k + 2) psi) g_{k+1} + Sum_{m=0}^{k} l_m g_{k-m}   in Lam_{k+1}.                                                        (6.3)
Hence sigma = t^-1 (0, E) + Sum_{k>=0} t^k sigma^(k), both components of sigma^(k) in Lam_{k+1}; explicitly
    sigma^(k) = Sum_{j=1}^{k+1} S_j^(k) lam_j,   S_j^(k) := ( g_j^(k), (k + 2 - j) a_j^(k+1) + Sum_{m=0}^{k} l_m g_j^(k-m) ).                                 (6.4)
In particular sigma^(0) = (E lam_1, a1 a2 lam_1) = E lam_1 (1, L_2(0)), and g_1 = (u/(2E))(2 lam_2 - lam_1).
Proof. Verbatim [PR] Proposition 6.1 with c -> E: [t^k](L_1 g) = E g_{k+1} + Sum_m l_m g_{k-m}, [t^k]((psi/t) Dg) = (k + 2) psi g_{k+1}; for k = -1,
[t^-1] sigma_1 = (E - psi) g_0 = (E - psi) E lam_1 = E by (LAM3). Leaf structure: by Theorem 5.4_CY and (LAM4), Sum_j g_j^(k+1)/j = 0; by (LAM3) with E,
(E - (k + 2) psi) g_{k+1} = (k + 2) Sum_j g_j^(k+1)/j + E Sum_j (1 - (k + 2)/j) g_j^(k+1) lam_j = Sum_{j<=k+1} (k + 2 - j) a_j^(k+1) lam_j: the factor E - (k + 2) psi turns the
top leaf lam_{k+2} into a constant, and the constants add up to 0. The terms l_m g_{k-m} lie in Lam_{k-m+1}.
k = 0, by hand. From s.3.11 and (4.3) (tau/t = 1 - t/2 + O(t^2)), G_1 = (E/(E + z))(H_1 - 1/2) (= [t^1] r_{0,1} + (u/(2E^2)) G_0, since G = rho r_{0,1} and
rho_1 = u/(2E^2)) = [2 E z + 2 E (a2 + A) - 4 a2 A - E^2 - 2 E z]/(2 (E + z)(E + 2z)) = u/(2 (E + z)(E + 2z)),
because 2 E (a2 + A) - 4 a2 A - E^2 = -(E - 2 a2)(E - 2 A) = -(E - 2 a2)(2 a1 - E) = u. So R_1 = (u/2)[1/(E + z) - 1/(E + 2z)], (a_1^(1), a_2^(1)) = (u/2)(1, -1),
q^(1) = (u/2)(lam_1 - lam_2), g_1 = -(1/E)(a_1^(1) lam_1 + 2 a_2^(1) lam_2) = (u/(2E))(2 lam_2 - lam_1). Then (E - 2 psi) g_1 = a_1^(1) lam_1 = (u/2) lam_1 (by (LAM3):
(E - 2 psi) lam_2 = 1, (E - 2 psi) lam_1 = 2 - E lam_1), l_0 g_0 = (E^2/2 - a0 E - a1 a2) lam_1, and [t^0] sigma_1 = (u/2 + E^2/2 - a0 E - a1 a2) lam_1 = a1 a2 lam_1,
since (u + E^2)/2 = a0 E + 2 a1 a2. QED
Remarks. (i) sigma^(0)_1 = a1 a2 lam_1 is the z-dependent part a1 a2/(E + z) of the H-row at t^0: rho r_{1,1} = E/t + a1 a2/(E + z) + O(t), the second half of
F2; on T_F, x/(c + z). (ii) The two cancellations in the k = 0 computation, 2 E (a2 + A) - 4 a2 A - E^2 = u and u/2 + l_0 E = a1 a2, are the generic form of
[PR]'s (1 - x)^2 and x; they use only E = a1 + A = a2 + B and A B = a1 a2 + a0 E. [Evidence: check2_rows.py: F2 for all z^m <= 9 (11).]

### 6.3 The projection of sigma on the regular eigen-covector
**Proposition 6.2_CY.** pi(sigma(psi)) = -Theta (psi/t) Dg(psi) = e1 e2 B tau s* L_2 (psi/t) Dg(psi) in Kw[[t]][[psi]], Theta = e1 e2 (a1 a2 - E L_2).
Proof. sigma = g chi_1 + (sigma_1 - L_1 g)(0, 1); pi(chi_1) = 0, pi((0, 1)) = e1 e2 (a1 a2 - E L_2) = Theta (Proposition 4.7(e)), and sigma_1 - L_1 g = -(psi/t) Dg by (6.1). QED
On T_F: x^3 (tau/t) s*^2 psi Dg, [PR] Proposition 6.2.

## 7. Step 8: the E_11 flow; F3 for E_11; L3; L2 for j = 1 [verbatim after substitution; bookkeeping re-done]

### 7.1 E_11 in terms of sigma
By Lemma 1.1 and Proposition 4.8(d), rho^2 E_11(psi, psi') = [rho^2 Delta_1 - <sigma(psi), sigma(psi')>]/(psi + psi'), rho^2 Delta_1 = -e1 e2 E^3/t^2. Put Y(psi) := t sigma(psi)
in Kw[[t]][[psi]]^2 and
    Ebar(psi, psi') := t^2 rho^2 E_11(psi, psi') = Nbar/(psi + psi'),   Nbar := -e1 e2 E^3 - <Y(psi), Y(psi')>.                                            (7.1)
Nbar(psi, -psi) = -e1 e2 E^3 - t^2 rho^2 Delta_1 = 0 (Proposition 4.8(c), (4.2)), so Nbar is divisible by psi + psi', Ebar is in Kw[[t]][[psi, psi']], rho^2 E_11 is in
t^-2 Kw[[t]][[psi, psi']], and Ebar(psi, psi') = Ebar(psi', psi), i.e. E_11(a, b) = E_11(b, a). (Verbatim [PR] s.7.1.)

### 7.2 The QDE for Y
**Lemma 7.1_CY.** -psi dY/dt = (A - L_1) Y in Kw((t))[[psi]]^2. (Proposition 4.8(d); verbatim.)

### 7.3 The flow
**Theorem 7.2_CY (E_11 flow).** In Kw[[t]][[psi, psi']],
    d/dt Ebar(psi, psi') = Theta(t) Dg(psi) Dg(psi'),   Dg = d/dt (t g),                                                                        (7.2)
    Theta = Delta_2/(L_2 - L_1) = e1 e2 (a1 a2 - E L_2) = -e1 e2 B tau s* L_2 = e1 e2 (E dphi' - a2 A).                                                (7.3)
Proof. Steps (i), (ii) verbatim [PR] Theorem 7.2: psi psi' dEbar/dt = <(A - L_1) Y(psi), Y(psi')> by Lemma 7.1 and self-adjointness. (iii) By (4.8) and
Proposition 6.2, <(A - L_1) Y(psi), Y(psi')> = pi(Y(psi)) pi(Y(psi'))/Theta with pi(Y(psi)) = t pi(sigma(psi)) = -Theta psi Dg(psi); so psi psi' dEbar/dt =
Theta psi psi' Dg(psi) Dg(psi'); cancel psi psi'. (iv) The forms of Theta are Proposition 4.7(e), 4.3(d), and dphi' = a2 - L_2 gives a1 a2 - E L_2 = a1 a2 - E a2 + E dphi'
= E dphi' - a2 A (a2 E - a1 a2 = a2 A). QED
Expansion. Theta_0 = e1 e2 (E dphi_1 - a2 A) = 0 since dphi_1 = a2 A/E; Theta = e1 e2 E dphi' - e1 e2 a2 A gives Theta_k = (k + 1) e1 e2 E [t^(k+1)] dphi for k >= 1, and from
the explicit form Theta_1 = -e1 e2 B s*(0) L_2(0) = -e1 e2 B (A/E)(a1 a2/E) = -e1 e2 a1 a2 A B/E^2 (tau = t + O(t^2)); hence [t^2] dphi = -a1 a2 A B/(2 E^3).
On T_F: Theta = -x^3 tau s*^2 = x c dphi' - x, Theta_1 = -x^3/c^2, [t^2] dphi = -x^2/(2 c^3), [PR] (7.3). Remarks (a), (b) of [PR] hold verbatim: Theta = O(t) because
Delta_2 is regular and L_2 - L_1 has a simple pole; the right side of (7.2) is of rank one before expanding in t. [Evidence: check2_rows.py: (7.2) for
m + n <= 3, t^0..t^1 (t^3), from the Wick-built rows; Theta_1.]

### 7.4 Integration from t = 0; F3 for E_11
[t^0] Y = (0, E) (Proposition 6.1(a)) and <(0, E), (0, E)> = e1 e2 (-E) E^2 = -e1 e2 E^3, so [t^0] Nbar = 0 and [t^0] Ebar = 0: this is where K_1 = 1 enters (for a
t-independent K_1 one would have [t^0] Ebar = -e1 e2 E^3 (K_1(psi') - 1)/(psi + psi')). Hence, as in [PR] s.7.4,
    Ebar = Int_0^t Theta(u) Dg(psi; u) Dg(psi'; u) du  (termwise),   [t^(n+1)] Ebar = [t^n](Theta Dg(psi) Dg(psi'))/(n + 1).                                  (7.4)
As Theta_0 = 0, [t^1] Ebar = 0; Ebar is in t^2 Kw[[t]][[psi, psi']] and E_11 = (t rho)^-2 Ebar is in Kw[[t]][[psi, psi']]: E_11 is regular at t = 0 (F3 for E_11).

### 7.5 The leaf structure of the loops (L3)
**Theorem 7.3_CY (L3, all k).** For every k >= 0,
    [t^k](rho^2 E_11)(psi, psi') = (1/(k + 2)) Sum_{k0>=1, k1+k2 = k+1-k0} Theta_{k0} (k1 + 1)(k2 + 1) g_{k1}(psi) g_{k2}(psi') = Sum_{(i-1)+(l-1) <= k} e^(k)_il lam_i(psi) lam_l(psi'),   (7.5)
    e^(k)_il = (1/(k + 2)) Sum_{k0>=1, k1+k2 = k+1-k0} Theta_{k0} (k1 + 1)(k2 + 1) g_i^(k1) g_l^(k2),   e^(k)_il = e^(k)_li,
no lam_0 in either variable; on the anti-diagonal i = k1 + 1, l = k2 + 1, k1 + k2 = k,
    e^(k)_{k1+1, k2+1} = (-1)^(k+1) e1 e2 a1 a2 A B (k1 + 1)(k2 + 1) P_{k1} P_{k2}/((k + 2) (k1!)^2 (k2!)^2 E^(2k)),                                      (7.6)
non-zero iff P_{k1} P_{k2} != 0. At t = 0 the loop has rank one: [t^0](rho^2 E_11) = (1/2) Theta_1 g_0(psi) g_0(psi') = -(e1 e2 a1 a2 A B/2) lam_1(psi) lam_1(psi'), hence
E_11(z, w)|_{t=0} = -(e1 e2 a1 a2 A B/2) lam_1(z) lam_1(w).
Proof. Verbatim [PR] Theorem 7.3 from (7.4), Theorem 5.4_CY and (7.3): support (i - 1) + (l - 1) <= k1 + k2 <= k; symmetry; the anti-diagonal comes only from
k0 = 1, with coefficient (1/(k + 2)) Theta_1 (k1 + 1)(k2 + 1) g^(k1)_{k1+1} g^(k2)_{k2+1}, and g^(k)_{k+1} = -((k + 1)/E) a_{k+1}^(k) = (-1)^k P_k/((k!)^2 E^(2k-1)) (Corollary 11.3),
Theta_1 = -e1 e2 a1 a2 A B/E^2; collecting: (-1)^(k1+k2) (-1) e1 e2 a1 a2 A B (k1 + 1)(k2 + 1) P_{k1} P_{k2}/((k + 2)(k1!)^2 (k2!)^2 E^(2 + 2k1 - 1 + 2k2 - 1)) and 2 k1 + 2 k2 = 2k. QED
On T_F: e1 e2 a1 a2 A B -> x^3, [PR] (7.6). [Evidence: check2_rows.py: L3 bidegree bound for k <= 2 on psi, psi' orders <= 4; E_11(0, 0) and E_11(1, 0) at t = 0; the critic's
leaf_high.py (read, not re-run): bidegree bound for k <= 6, m + n <= 14.]

### 7.6 Powers of rho; L2 for j = 1
**Corollary 7.4_CY.** For every integer e and k >= 0, [t^k](rho^e E_11) lies in span{lam_i(psi) lam_l(psi') : (i - 1) + (l - 1) <= k}; in particular (L2 for j = 1)
    [t^k](rho f_{1,b}) = Sum_{i=1}^{k+1} ( Sum_{m=0}^{k} [t^m](rho^-1) Sum_l e^(k-m)_il l^b/E^(b+1) ) lam_i(psi)   in Lam_{k+1},                              (7.7)
and [t^k] f_{1,b} lies in Lam_{k+1}. (Verbatim [PR] Corollary 7.4: rho^(e-2) in Kw[[t]], [psi'^b] lam_l(psi') = l^b/E^(b+1).)
Explicitly [t^0](rho f_{1,b}) = -(e1 e2 a1 a2 A B/2) E^(-b-1) lam_1 for every b (rank one at t = 0).

## 8. The regular idempotent (Step 10, and the second QDE row at U_2)

### 8.1 The conjugate encoding [verbatim]
**Lemma 8.1_CY.** [PR] Lemma 8.1 holds with the generic LeqU and Ssum, Pprod of (C2'): cconj is a ring automorphism of QA commuting with the ring
operations, with cDerN and with multiplication by scalar series; since h_nw and h_w^-1 are the same scalar series at both idempotents (Lemma 13.3 of [PR]:
h_{wt_+}/r^pt_{0,+}(z; 0) = h_{wt_-}/r^pt_{0,-}(z; 0) = h_nw), the stored quantities of idempotent 2 are cconj of those of idempotent 1 before and after the
twists, and the true values at idempotent 2 are ev_1 of the stored elements: r^pt_{0,2}(z; t) = Sum_m z^m wickZ_m(W_2(U_2), ...), r_{1,2} = (L_2 + z D_2) r_{0,2}.
(This is the only place where "the same scalar at both colours" is needed; with a colour-dependent constant the conjugate encoding would fail.)

### 8.2 Analyticity at U_2 [substitution]
**Lemma 8.2_CY.** (i) For k >= 1, W_k is a polynomial in p_1 = (1 - U)^-1 and p_Q = (1 - Q U)^-1 with coefficients in Z[a1, a2, A] (W_1 = a1 (1 - p_Q) - a2 + A (1 - p_1),
and (3.7) for k >= 2). U_2(t) is analytic at t = 0 with U_2(0) = -a2/B, and 1 - U_2(0) = 1 - Q(0) U_2(0) = E/B != 0; so every
W_k(U_2(t)) is analytic at t = 0, and W_2(U_2) is a unit of O with W_2(U_2(0)) = a2 B/E.
 (ii) (R) Every r^pt_{0,2}[m], r_{0,2}[m], r_{1,2}[m] is analytic at t = 0; r_{0,2}[0] = 1, r_{1,2}[0] = L_2.
 (iii) Delta_2 = e1 e2 B delta s* L_2 = a1 A e1 e2 W_2(U_2) is analytic with Delta_2(0) = e1 e2 a1 a2 A B/E and Delta_2'(0) = 0 (F1).
 (iv) Dil_2(n) and every E_22(m, n) are analytic at t = 0 (F3 for Dil_2, E_22).
Proof. (i) theta log(1 - a U) = 1 - p_a, theta p_a^n = n (p_a^(n+1) - p_a^n). At Q = 1, LeqU = B U^2 - (B - a2) U - a2 = (B U + a2)(U - 1), with the simple roots 1 and
-a2/B (distinct since -a2/B = 1 iff E = 0); so U_2 is a simple root at Q = 1 and analytic in t. At Q = 1 the superpotential is W_0(U) = E log(1 - U) - a2 log U,
W_2 = theta^2 W_0 = -E U/(1 - U)^2, and at U = -a2/B: -E (-a2/B)/(E/B)^2 = a2 B/E. (ii) r^pt_{0,2}[m] is a polynomial in W_2(U_2)^-1, W_3(U_2), ..., hence analytic;
h_nw, h_w^-1 are t-independent; L_2 = B s* - a0 in O, Delta_2 a unit. (iii) Proposition 4.3(d), (e); the value a1 A e1 e2 a2 B/E agrees. (iv) Lemma 1.1. QED
On T_F: U_2(0) = -1/x, W_2(U_2(0)) = x/c, Delta_2(0) = x^3/c, [PR] Lemma 8.2.

### 8.3 The data at the regular critical point [substitution]
**Lemma 8.3_CY.** For tau in D'':
 (a) s~ = (E + a0 tau + delta)/(2 B tau) = A/(B tau s*) is a nondegenerate critical point of phi_tau, not in {0, 1, 1/tau}; tau s~ is in O with value E/B;
     2 B tau s~ - (E + a0 tau) = delta; B s~ = L_1 + a0; U_2 = 1 - tau s~ is a unit of O with U_2(0) = -a2/B.
 (b) f2 := phi_tau''(s~) = delta/(s~ (1 - s~) U_2) = tau^2 Delta_2/(a1 A e1 e2 U_2^2) = tau^2 f2r, f2r := Delta_2/(a1 A e1 e2 U_2^2) a unit of O with
     f2r(0) = Delta_2(0)/(a1 A e1 e2 U_2(0)^2) = (a1 a2 A B/E) B^2/(a1 A a2^2) = B^3/(a2 E) (on T_F: x^3/c).
 (c) f0(tau) := phi_tau(s~(tau)) has f0_tau = a2 s~/(1 - tau s~) in tau^-1 O with leading term -E/tau (tau s~ -> E/B, U_2 -> -a2/B, so f0_tau ~ a2 (E/B)/(tau (-a2/B)) = -E/tau), and
        f0' = (1 - tau) a2 s~/(1 - tau s~) = A - B s~ = a2 - L_1 = L_2 - mu.                                                                      (8.1)
     u_2 := W(U_2) satisfies du_2/dt = L_2 and u_2 = a1 t + E log tau + f0 + const.
Proof. (a) s.1.3, Lemma 3.1(a), Lemma 4.2(d). (b) (3.3) at s~ and (3.5) with Delta_2 = a1 A e1 e2 W_2(U_2) (Proposition 4.3(c)). (c) Envelope formula and Lemma
3.1(b)(ii) at the root s~ (the identity (1 - tau) a2 s = (A - B s)(1 - tau s) holds at both roots, being equivalent to N(s) = 0); B s~ = A/(tau s*) = L_1 + a0 so
A - B s~ = a2 - L_1, and L_2 - mu = L_2 - E/tau + A = a2 - L_1 by L_1 + L_2 = E/tau - a0. QED
On T_F: f0' = 1 - x s~ = 1 - L_1, [PR] (8.1).

### 8.4 The second QDE row at the regular idempotent [substitution; alternatively cited]
**Theorem 8.4_CY (Row 2 at U_2).** In M[[z]],
    (L_2 + z D_2) r_{1,2} = -(a1 a2/tau) r_{0,2} + (E/tau - a0) r_{1,2},   equivalently   z D_2 r_{1,2} = L_1 (r_{1,2} - L_2 r_{0,2}),                        (8.2)
i.e. z D_2 v_2 = (A - L_2) v_2, for the frame-B rows, the R_Z rows and the point-frame rows alike.
Proof. Steps A-D of [PR] Theorem 8.4 with the generic data. Step A: r^pt_{0,2}(z; t) = U_2 L2(tau), L2 := Lnorm[phi_tau, (1 - tau s)^-1; s~] in O_Db[[z]] (Proposition
3.10_CY, steps (1)-(2), at s~, using (3.1) and W_2(U_2) != 0 on D''). Step B: with Dt := d/dt + f0'/z - f2'/(2 f2), (L_2 + z D_2) o U_2 = U_2 o (z Dt + omega) (8.4), because
Delta_2 = a1 A e1 e2 U_2^2 f2/tau^2 gives Delta_2'/(2 Delta_2) = U_2'/U_2 + f2'/(2 f2) - (1 - tau)/tau, and L_2 = f0' + mu (8.1), mu + z (1 - tau)/tau = omega. Step C: Gauss_Dop[L2] = 0
on O_Db((z)) by (P2), (D3), the generic amplitude identity (3.16) (valid identically in (s, tau), hence also near s~) and (P3). Step D: Lemma 4.5_CY with Dd' = z Dt gives
(z Dt + omega)^2 - (E/tau - a0)(z Dt + omega) + a1 a2/tau = ((1 - tau) z^2/tau) Gauss_Dop, so [(L_2 + z D_2)^2 - (E/tau - a0)(L_2 + z D_2) + a1 a2/tau] r^pt_{0,2} = 0; the
t-independent scalars h_nw, h_w^-1 commute with everything, and (8.2) follows on Db, hence on D'' and in M by Lemma 8.2_CY. The second form uses E/tau - a0 - L_2 = L_1
and a1 a2/tau = L_1 L_2. QED
Citation. [PR] Proposition 13.2(b)(ii) proves Row 2 at both roots at indeterminate weights by the same argument with Gauss's operator replaced by the
Picard-Fuchs operator and (3.16) by Euler's identity (13.3) PF e^{W/z} = theta(km e^{W/z}), "at every Q with distinct roots not in {0, 1, 1/Q}". That statement
is pointwise in Q and in the weights and does not use Q = 0, so it covers every weight point with (ND), a0 = 0 included, and both colours; it can be cited in place
of Proposition 4.6_CY and Theorem 8.4_CY. (The two routes are the same calculation in two coordinates.)
[Evidence: for colour 1, check2_rows.py (Row 2); for colour 2, the reports' T2 tests (cy_torus_explore), not re-run here.]

### 8.5 Orthogonality [verbatim]
**Corollary 8.5_CY.** <v_1(w), v_2(-w)> = 0 in M[[w]], to all orders in w.
Proof. Verbatim [PR] Corollary 8.5: Y12 := <v_1(z), v_2(-z)> satisfies z Y12' = [(L_2 - L_1) + (z/2)(Delta_1'/Delta_1 + Delta_2'/Delta_2)] Y12 by Proposition 4.6_CY, (8.2),
the t-independence of eta^-1 and self-adjointness (Proposition 4.7(c)); L_2 - L_1 = -delta/tau is a unit of M, so Y_0 = 0 and inductively Y_m = 0. QED
Citation. [PR] Proposition 13.2(b)(iv) proves <v_+(z), v_-(-z)> = 0 and <v_i(z), v_i(-z)>/Delta_i = h_{wt_i}(z) h_{wt_i}(-z) = 1 at indeterminate weights, the second by
evaluation at Q = 0. The first is the ODE argument above (no Q = 0 needed; covers a0 = 0). The second uses the semisimple point Q = 0, so it needs a0 != 0: at a0 = 0
both roots run to infinity at Q = 0 ([PR] Remark (iv) after Proposition 13.2). For K_1 = 1 this text uses Proposition 4.8_CY (the t -> 0 limit at Q = 1), which covers
a0 = 0; for K_2 = 1, Corollary 8.7_CY below (analytic input) or, at a0 != 0, [PR] 13.2(b)(iv); by rationality in the weights (s.1.4 (P)) an identity of the
t-, z-coefficients valid for a0 != 0 extends to a0 = 0 on (ND). K_2 is not needed for Theorem LC_CY.

### 8.6 The regular rows at t = 0 [substitution; (a) is an analytic input]
**Proposition 8.6_CY.** (a) [analytic input] r^pt_{0,2}(z; 0) = h_{(a2, -E, B)}(z); equivalently the frame-B value r_{0,2}(z; 0) = h_{(a1, A, -E)} h_{(a2, -E, B)} =
     h_{(a1, a2, A, B, -E, -E)}(z) =: hreg(z), with hreg(z) hreg(-z) = 1.
 (b) (Uses Theorem 8.4_CY.) a1 a2 r_{0,2}(z; 0) = E r_{1,2}(z; 0), i.e. v_2(z; 0) = r_{0,2}(z; 0) (1, a1 a2/E) = hreg(z) (1, a1 a2/E). Equivalently [t^1] r_{0,2}(z; t) = 0.
 (c) (Uses (a), (b).) E_22(z, w)|_{t=0} = (e1 e2 a1 a2 A B/E) [1 - hreg(z)^-1 hreg(w)^-1]/(z + w), and Dil_2(n)|_{t=0} = (-1)^n [z^(n-1)] hreg(z).
Proof. (a) By Lemma 8.2_CY, r^pt_{0,2}(z; 0) = Lnorm[Wy_0, 1; y_c], Wy_0(y) := W_0(e^y) with W_0 = E log(1 - U) - a2 log U = a log U + b log(1 - U), (a, b) := (-a2, E),
at the critical point U = a/(a + b) = -a2/B. This is exactly the limit computed in [PR] Lemma 13.3 ("This limit is h_{(-a, -b, a+b)}": both sides rational in (a, b),
Laplace's method for a, b > 0 on the integral Int_0^1 U^(a/z) (1 - U)^(b/z) dU/U = B(a/z, b/z + 1), Stirling's series, uniqueness of asymptotic expansions), which
needs only a b (a + b) != 0, here a2 E B != 0. So r^pt_{0,2}(z; 0) = h_{(a2, -E, B)}. Frame B multiplies by h_{(a1, A, -E)} (C8'), and the union (a1, A, -E, a2, -E, B) is
(a1, a2, A, B, -E, -E). On T_F this is h_{(1, -c, x)} h_{(x, 1, -c)} = h_w^-2, [PR] Proposition 8.6(a).
(b) First proof: in (D_2 r_{1,2})[m-1] = L_1 (r_{1,2}[m] - L_2 r_{0,2}[m]) the left side is analytic at t = 0 (Lemma 8.2_CY) and L_1 = E/t + O(1) has a simple pole, so
r_{1,2}[m](0) = L_2(0) r_{0,2}[m](0) = (a1 a2/E) r_{0,2}[m](0) for m >= 1, and for m = 0 by (H). Second proof: from Corollary 8.5_CY, <sigma(z), v_2(z)> = 0, whose t^-1
coefficient is <(0, E), v_2(z; 0)> = e1 e2 E (a1 a2 r_{0,2}(z; 0) - E r_{1,2}(z; 0)) (by <(0, E), b> = e1 e2 [a1 a2 E b_0 - E^2 b_1]). Equivalent form: by (H),
r_{1,2}[m](0) - (a1 a2/E) r_{0,2}[m](0) = (D_2 r_{0,2}[m-1])(0) = (r_{0,2}[m-1])'(0), because L_2(0) = a1 a2/E and Delta_2'(0) = 0 (Proposition 4.3(e)).
(c) By Lemma 1.1, P_22(z, w)|_{t=0} = -hreg(-z) hreg(-w) <(1, a1 a2/E), (1, a1 a2/E)>, and <(1, a1 a2/E), (1, a1 a2/E)> = e1 e2 [a0 a1 a2 + 2 a1 a2 (a1 a2/E) - E (a1 a2/E)^2]
= e1 e2 a1 a2 (a0 E + a1 a2)/E = e1 e2 a1 a2 A B/E = Delta_2(0); P_22(-w, w)|_{t=0} = -Delta_2(0) since hreg(w) hreg(-w) = 1. QED
Remark 8.8 of [PR] (E_22 at t = 0 is not of rank one) carries over: (1 - hreg(z)^-1 hreg(w)^-1)/(z + w) is not a product; only the regularity of E_22 is used.

**Corollary 8.7_CY (K_2 = 1).** K_2(z) := <v_2(z), v_2(-z)>/Delta_2 is in O[[z]], t-independent and equal to r_{0,2}(z; 0) r_{0,2}(-z; 0) = hreg(z) hreg(-z) = 1.
(Verbatim [PR] Corollary 8.7 with Theorem 8.4_CY and Proposition 8.6_CY(a), (b). Not used for Theorem LC_CY.)

Remark 8.9 of [PR] (the regular column is determined recursively by the singular data) and Remark 8.10 (module language) carry over verbatim with the
substitutions of s.1.5; a1 a2 r_{0,2} - E r_{1,2} = (a1 a2 - E L_2) r_{0,2} - E z D_2 r_{0,2} with a1 a2 - E L_2 = -B tau s* L_2 = O(t). Not used.

## 9. Step 9: the E_12 and E_21 legs (L2 for j = 2; F3 for E_12, E_21) [verbatim]

### 9.1 Extraction of a leg
Lemma 9.1 of [PR] ([w^b] E = Sum_{n<=b} (-1)^n psi^(-n-1) [w^(b-n)] Num for Num = (psi + w) E) has no weights; verbatim.

### 9.2 L2 for j = 2 and the regularity of E_12
**Proposition 9.2_CY.** Assume (R) (Lemma 8.2_CY(ii)). Then (i) E_12(a, b) is regular at t = 0 (F3 for E_12); (ii) [t^k](rho f_{2,b}) lies in Lam_{k+1} for all k, b >= 0.
Proof. Verbatim [PR] Proposition 9.2 with (0, c) -> (0, E), c/j -> E/j: N12 := rho P_12 = -<sigma(psi), v_2(-w)> is a Laurent series in t starting at t^-1 whose t^-1
coefficient -<(0, E), v_{2,0}(-w)> does not depend on psi (so [t^-1](rho E_12) = 0); for (ii), by Proposition 6.1_CY both components of every sigma^(k1), -1 <= k1 <= k,
lie in Lam_{k1+1}^+ and the rows v_2 enter only through the constants v_{2,k2}[m] and ym_m; the rational function Rat(psi) has simple poles only at E/j, j <= k + 1,
tends to 0 at oo, and its principal part at 0 vanishes because E_12 is a power series in psi; so Rat is in Lam_{k+1}. QED
What is used about the regular rows: only (R). The t = 0 value is s.11.3(b).

### 9.3 The leaf form of the E_12 edge
**Proposition 9.3_CY.** With S_j(t) := Sum_{k>=j-1} t^k S_j^(k) (6.4):
    sigma(psi) = (0, E/t) + Sum_{j>=1} S_j(t) lam_j(psi),   rho E_12(psi, w) = -Sum_{j>=1} j <S_j(t), v_2(-w; t)> lam_j(psi) lam_j(-w),   lam_j(-w) = 1/(E + j w),      (9.6)
and [t^k](rho f_{2,b})(psi) = -Sum_{j=1}^{k+1} j lam_j(psi) Sum_{k1=j-1}^{k} <S_j^(k1), [w^b](v_{2,k-k1}(-w)/(E + j w))> (9.7).
(Verbatim [PR] Proposition 9.3: the term (0, E/t) cancels between N12(psi, w) and N12(-w, w), and (LAM6) with E.)

### 9.4 The transposition E_21 = E_12^T; regularity of E_21
**Proposition 9.4_CY.** Assume Corollary 8.5_CY. Then E_21(z, w) = E_12(w, z), E_21 is regular at t = 0, and an edge between a singular and a regular vertex gives the
singular vertex the leg f_{2,b} for either orientation. (Verbatim [PR] Proposition 9.4.)
Without orthogonality: [t^-1](rho E_21)(z, w) = [fe(z) - fe(-w)]/(z + w) with fe(z) := -<v_2(-z; 0), (0, E)> = -e1 e2 E [a1 a2 r_{0,2}(-z; 0) - E r_{1,2}(-z; 0)], which
vanishes iff fe is constant, i.e. iff Proposition 8.6_CY(b) (the constant is fe(0) = -e1 e2 E (a1 a2 - E L_2(0)) = 0). So already the regularity of E_21 needs the
t^-1 part of the orthogonality, i.e. Row 2 at U_2.

### 9.5 All legs; unrescaled ends
**Corollary 9.5_CY.** For all k, b >= 0 and j in {1, 2}: [t^k] f_{j,b} and [t^k](rho f_{j,b}) lie in Lam_{k+1}; [t^k](rho^e E_11) lies in the span of (7.5) for every e;
all these series are O(t^0). (Verbatim [PR] Corollary 9.5.)

### 9.6 Use in Theorem R
Verbatim [PR] s.9.6: at a singular vertex of genus >= 2 every half-edge carries one factor rho; E_12/E_21 give the leg rho f_{2,b}; an edge to a non-rescaled singular
vertex gives rho f_{1,b'} of leaf degree <= k at order t^k; an edge between two rescaled singular vertices, or a self-loop, carries rho^2 E_11, which splits its D <= k
between its ends. At genus-0 and genus-1 singular vertices only regularity is needed: Delta_1^(h-1) is O(t^2) for h = 0 and 1 for h = 1.

## 10. Proof of Theorems LC_CY and Q_CY; F1-F3 in one table

Proof of Theorem LC_CY.
 F1: t^2 rho^2 Delta_1 = -e1 e2 E^3, rho = tau/(t alg), rho(0) = 1: Proposition 4.4_CY(b), (d). Delta_2 = e1 e2 B delta s* L_2 analytic with Delta_2(0) = e1 e2 a1 a2 A B/E,
     Delta_2'(0) = 0: Proposition 4.3_CY(d), (e), Lemma 8.2_CY(iii).
 F2: r_{0,1}(z; 0) = E/(E + z): Corollary 3.18_CY. rho r_{1,1} = E/t + a1 a2/(E + z) + O(t): Proposition 6.1_CY(a), (b) at psi = -z.
 F3: table below.
 L1: Theorem 5.4_CY.  L2: j = 1 by Corollary 7.4_CY, j = 2 by Proposition 9.2_CY(ii); unrescaled legs by Corollary 9.5_CY; E_21(b, a) = E_12(a, b) by Proposition 9.4_CY.
 L3: Theorem 7.3_CY; symmetry s.7.1; all powers rho^e by Corollary 7.4_CY.
 Sharpness: Corollaries 11.3_CY, 11.5_CY. Large-z form: Lemma 5.2_CY and Corollary 9.5_CY. QED
Proof of Theorem Q_CY. (i) Proposition 4.6_CY and Theorem 8.4_CY (or [PR] 13.2(b)(ii)). (ii) Proposition 4.8_CY, Corollary 8.5_CY, Corollary 8.7_CY. (iii) Proposition 8.6_CY. QED

F3, factor by factor (frame B; values at t = 0; hreg := h_{(a1, a2, A, B, -E, -E)}, kappa_2 := e1 e2 a1 a2 A B):
    factor  regularity                                             value at t = 0
    Dil_1   Theorem 5.4_CY (iv)                                    q_0 = -E psi/(E - psi)   (Dil_1(d; 0) = -E^(1-d))
    Dil_2   Lemma 8.2_CY(iv)                                       (-1)^n [z^(n-1)] hreg(z)   (Proposition 8.6_CY(a))
    E_11    s.7.4                                                  -(kappa_2/2) lam_1(psi) lam_1(psi')   (Theorem 7.3_CY)
    E_12    Proposition 9.2_CY(i)                                  lam_1(psi) chi12(w), chi12(w) = -kappa_2 hreg(w)^-1/(E + w)   (s.11.3(b))
    E_21    Proposition 8.6_CY(b) for the t^-1 part; 9.4_CY        E_21(b, a) = E_12(a, b)
    E_22    Lemma 8.2_CY(iv)                                       (kappa_2/E) [1 - hreg(z)^-1 hreg(w)^-1]/(z + w)   (Proposition 8.6_CY(c))
On T_F: kappa_2 = x^3, hreg = h_w^-2, [PR] s.10.

Where each external input is used. [PR] Proposition 3.15 (through Proposition 3.15_CY and the constant line (3.9a), which also uses [PR] Lemma 13.3 for the
definition of the rows) enters through Theorem 3.17_CY into the KEY FORMULA and everything that uses E/(E + z): F2, L1 at t^0, sigma^(-1) = (0, E), K_1 = 1,
[t^0] Ebar = 0, F3 for E_11, the rank-one t^0 loop. V1 at U_1 (Proposition 4.3_CY(c)) enters Proposition 4.4_CY; V1 at U_2 enters (8.4). The formal Laplace calculus
enters Propositions 3.10_CY, 3.12_CY and Theorem 8.4_CY. The condition (ND): E != 0 (leaf scale, Gauss recursion, distinct roots at Q = 1), a1 A != 0 (nondegenerate
Euler critical point, Proposition 3.15 at x = a1/A), a2 != 0 (s* != 1/tau; U_2(0) != 0), B != 0 (U_2(0) finite, Delta_2(0) != 0), e1 e2 != 0 (eta invertible, Delta != 0).

## 11. Closed forms: Theorem CF_CY

Nothing in Theorem LC_CY depends on this section except the sharpness statements and the explicit values quoted as examples.

### 11.1 The dilaton leaf coefficients [substitution; residues re-derived]
**Theorem 11.1_CY.** For j >= 1 put C_j := P_{j-1}/(j! (j - 1)! E^(2j-4)) = E^2 (1 - j a2/E)_{j-1} (1 - j A/E)_{j-1}/(j! (j - 1)!) (C_1 = E^2), zeta(t) := tau e^(dphi(tau)/E) = t + O(t^2),
    A_j(t) := (-1)^j C_j t^-1 zeta(t)^j 2F1(j B/E, j a1/E; j; tau),   A_j = O(t^(j-1)).
Then a_j^(k) = [t^k] A_j(t) for all j >= 1, k >= 0; equivalently z G = E + Sum_j A_j/(E + j z), rho q_t = E + Sum_j A_j lam_j, g = -(1/E) Sum_j j A_j lam_j.
Proof. a_j^(k) = j Res_{z=-E/j} [t^k](z G) (Proposition 5.3_CY), z G = (tau/t) Sum_n tau^n (E z/(E + z)) Ephi F_n, and [t^k](z G) is a finite sum over n <= k of terms
(scalar)(E z/(E + z)) z^-p F_n; evaluation at z = -E/j commutes with taking t-coefficients.
j = 1: F_n is regular at z = -E (its poles are at 0 and -E/l, l >= 2), so Res_{z=-E} [t^k](z G) = [t^k](-E^2 (tau/t) e^(dphi/E) F(-E, tau)), and at z = -E: alpha = 1 - a2/E = B/E,
beta = 1 - A/E = a1/E, gamma = 2 - 1 = 1, so F(-E, tau) = 2F1(B/E, a1/E; 1; tau); Sum_m X_m(-E) tau^m = e^(dphi/E). Hence A_1 = -E^2 t^-1 zeta 2F1(B/E, a1/E; 1; tau) = (-1)^1 C_1 t^-1 zeta 2F1.
j >= 2: E z/(E + z) and z^-p are regular at -E/j; F_n is regular there for n < j - 1 and has a simple pole for n >= j - 1 from the single factor E + j z. At z = -E/j:
a2 + l z = (j a2 - l E)/j, A + l z = (j A - l E)/j, z^n = (-E/j)^n, prod_{l=2, l != j}^{n+1} (E + l z) = (E/j)^(n-1) (j - 2)! (-1)^(n+1-j) (n + 1 - j)!, Res 1/(E + j z) = 1/j; so
    Res_{z=-E/j} F_n = (-1)^(j-1) E (1 - j a2/E)_n (1 - j A/E)_n/(j^2 n! (j - 2)! (n + 1 - j)!),
using prod_{l=1}^{n} (j a2 - l E) = (-E)^n (1 - j a2/E)_n and the same with A. Reindex n = j - 1 + m: (1 - j a2/E)_{j-1+m} = (1 - j a2/E)_{j-1} (j B/E)_m, since 1 - j a2/E + j - 1 = j (E - a2)/E = j B/E;
(1 - j A/E)_{j-1+m} = (1 - j A/E)_{j-1} (j a1/E)_m; (j - 1 + m)! = (j - 1)! (j)_m. So
    Sum_{n>=j-1} tau^n Res_{z=-E/j} F_n = (-1)^(j-1) E (1 - j a2/E)_{j-1} (1 - j A/E)_{j-1} tau^(j-1)/(j^2 (j - 2)! (j - 1)!) 2F1(j B/E, j a1/E; j; tau).
The remaining factors at z = -E/j are E z/(E + z) = -E/(j - 1) and Ephi = e^(j dphi/E). Multiplying by j and using j (j - 1)(j - 2)! = j!,
a_j^(k) = [t^k] { (-1)^j E^2 (1 - j a2/E)_{j-1} (1 - j A/E)_{j-1}/(j! (j - 1)!) t^-1 zeta^j 2F1 }, and E^2 (1 - j a2/E)_{j-1} (1 - j A/E)_{j-1} = P_{j-1}/E^(2j-4), because
prod_{i=1}^{j-1} (j a2 - i E) = (-E)^(j-1) (1 - j a2/E)_{j-1} and likewise for A. QED
On T_F: A_j = (-1)^j C_j t^-1 zeta^j 2F1(j x/c, j x/c; j; tau), [PR] Theorem 11.1 (B/E = a1/E = x/c there). [Evidence: check2_rows.py: A_j against [t^k](z G), k <= 5 (7).]

**Corollary 11.2_CY.** (a) Sum rule: Sum_j A_j(t) = -E^2, i.e. Sum_j (-1)^j C_j zeta^j 2F1(j B/E, j a1/E; j; tau) = E^2 log(1 - tau) (no pole of G at z = 0).
 (b) t z G(z) = -E log(1 - tau) + Sum_j Ahat_j(tau)/(E + j z), Ahat_j := t A_j a function of tau.
 (c) Theta = e1 e2 (E dphi' - a2 A) (7.3), Theta_1 = -e1 e2 a1 a2 A B/E^2, Theta_k = (k + 1) e1 e2 E [t^(k+1)] dphi (k >= 1).
 (d) The H-component: sigma = (0, E/t) + Sum_j S_j(t) lam_j with S_j = ( -(j/E) A_j, A_j' - ((j - 1)/t) A_j - (j/E) l A_j ), S_1(0) = (E, a1 a2) = E (1, a1 a2/E)
     (verbatim [PR] Corollary 11.2(f) with c -> E, using Sum_j (t A_j)' = -E^2).
 (e) Dilaton consistency: Sum_j j A_j(t) = -E^2 rho(t) (from Dil_1(1) = -1; not used).
Not re-derived here: the branch-free forms of rho and alg ([PR] 11.2(d)) and the closed form of e^(dphi) ([PR] 11.2(c)); generically e^(dphi) = (E s*/A)^A (E (1 - s*)/a1)^(a1) U_1^(-a2)
with 1 - s* = a1 U_1/(E - B tau s*) (Lemma 3.1(b)(ii)), which is all that is needed to define zeta.

**Corollary 11.3_CY (top and subleading coefficients).** For k >= 0, resp. k >= 1,
    a_{k+1}^(k) = (-1)^(k+1) C_{k+1} = (-1)^(k+1) P_k/((k + 1)! k! E^(2k-2)),   g_{k+1}^(k) = (-1)^k P_k/((k!)^2 E^(2k-1)),   a_k^(k) = (-1)^(k+1) k C_k u/(2 E^2).
Proof. A_j = (-1)^j C_j t^(j-1) (zeta/t)^j 2F1 with (zeta/t)(0) = 1; the first claim is j = k + 1 at t^k; g_{k+1}^(k) = -((k + 1)/E) a_{k+1}^(k). For the second,
[t^1]((zeta/t)^k 2F1(k B/E, k a1/E; k; tau)) = k (dphi_1/E - 1/2) + (k B/E)(k a1/E)/k = k [a2 A/E^2 - 1/2 + a1 B/E^2] = k (2 a2 A + 2 a1 B - E^2)/(2 E^2) = -k u/(2 E^2),
since 2 a2 A + 2 a1 B - E^2 = 2 E (a1 + a2) - 4 a1 a2 - E^2 = -u. QED
Direct computation of the top coefficient without Theorem 11.1 (as [PR]): z0 := -E/(k + 1); only n = k in (5.1) and only F_k in H_k has a pole at z0; E z0/(E + z0) = -E/k,
a2 + i z0 = ((k + 1) a2 - i E)/(k + 1), A + i z0 = ((k + 1) A - i E)/(k + 1), z0^k = (-E)^k/(k + 1)^k, prod_{i=2}^{k} (E + i z0) = E^(k-1) (k - 1)!/(k + 1)^(k-1), so
Res_{z0} F_k = (-1)^k P_k/((k + 1)^2 k! (k - 1)! E^(2k-1)) and a_{k+1}^(k) = (k + 1)(-E/k) Res_{z0} F_k = (-1)^(k+1) P_k/((k + 1)! k! E^(2k-2)).
Sharpness: C_{k+1} != 0 iff P_k != 0 iff E is not in {(k + 1) a2/i, (k + 1) A/i : 1 <= i <= k}. On T_F: a_k^(k) = (-1)^k k C_k (1 - x)^2/(2c^2), [PR] 11.3.

### 11.2 The loops
**Theorem 11.4_CY.** rho^2 E_11(psi, psi') = Sum_{i,l>=1} e_il(t) lam_i(psi) lam_l(psi'), e_il(t) = (i l/(E^2 t^2)) Int_0^t Theta(u) Ahat_i'(u) Ahat_l'(u) du
= -(i l e1 e2 B/(E^2 t^2)) Int_0^tau v (1 - v) s*(v) L_2(v) (dAhat_i/dv)(dAhat_l/dv) dv (the factor (1 - v) comes from du = dv/(1 - v) and Ahat' = (1 - v) dAhat/dv; corrected 2026-11-03), and e_il = O(t^((i-1)+(l-1))). (Verbatim [PR] Theorem 11.4 with Dg = -(1/E) Sum_j j Ahat_j' lam_j and
Theta = -e1 e2 B tau s* L_2.)
**Corollary 11.5_CY (anti-diagonal and corner).** For i + l = k + 2:
    e_il^(k) = Theta_1 i^2 l^2 a_i^(i-1) a_l^(l-1)/(E^2 (i + l)) = -e1 e2 a1 a2 A B i^2 l^2 a_i^(i-1) a_l^(l-1)/(E^4 (i + l)) = -(-1)^(i+l) e1 e2 a1 a2 A B i^2 l^2 C_i C_l/(E^4 (i + l))   (= (7.6));
with a_1^(0) = -E^2 the corner is e_{1,k+1}^(k) = ((k + 1)^2/(k + 2)) (e1 e2 a1 a2 A B/E^2) a_{k+1}^(k). Example k = 1: e_12^(1) = e_21^(1) = -2 e1 e2 a1 a2 A B u/(3 E^2), e_22^(1) = 0.
Non-zero iff P_{i-1} P_{l-1} != 0, so the support bound of L3 is sharp entry by entry at generic weights.
(Proof: Theta = Theta_1 t + O(t^2), Ahat_j' = j a_j^(j-1) t^(j-1) + O(t^j), integrate. On T_F: e_12^(1) = 2 (1 - x)^2 x^3/(3 c^2), [PR] s.11.4.)

### 11.3 The legs
 (a) E_11 legs: [t^k](rho f_{1,b}) = Sum_{k1+k2=k} [t^(k1)](rho^-1) Sum_{i,l} e_il^(k2) (l^b/E^(b+1)) lam_i(psi); top coefficient (lam_{k+1} at t^k)
     e_{k+1,1}^(k)/E^(b+1) = ((k + 1)^2/(k + 2)) e1 e2 a1 a2 A B a_{k+1}^(k)/E^(b+3). At t^0: -(e1 e2 a1 a2 A B/2) E^(-b-1) lam_1.
 (b) E_12 at t = 0 (rank one): E_12(psi, w)|_{t=0} = lam_1(psi) chi12(w), chi12(w) := -e1 e2 a1 a2 A B r_{0,2}(-w; 0)/(E + w) = -e1 e2 a1 a2 A B hreg(w)^-1/(E + w).
     Proof: at t^0, (9.6) has only j = 1 with S_1^(0) = (E, a1 a2) = E (1, a1 a2/E), and <(1, a1 a2/E), b> = e1 e2 [(a0 a1 a2 + a1 a2 a1 a2/E) b_0 + (a1 a2 - E a1 a2/E) b_1] =
     (e1 e2 a1 a2 A B/E) b_0 for every b (a0 E + a1 a2 = A B); so [t^0](rho E_12) = -E (e1 e2 a1 a2 A B/E) r_{0,2}(-w; 0) lam_1(psi) lam_1(-w), [t^0] E_12 = [t^0](rho E_12), then
     Proposition 8.6_CY(a) with hreg(-w) = hreg(w)^-1. On T_F: -x^3 h_w(w)^2/(c + w).
 (c) E_12 top coefficients, given Proposition 8.6_CY(a), (b): [lam_{k+1}][t^k](rho f_{2,b}) = (k + 1)^2 (e1 e2 a1 a2 A B/E^2) a_{k+1}^(k) [w^b](hreg(w)^-1/(E + (k + 1) w)).
     Proof: in (9.7) only j = k + 1, k1 = k contribute: -(k + 1) <S_{k+1}^(k), [w^b](v_{2,0}(-w)/(E + (k + 1) w))>, v_{2,0}(-w) = hreg(w)^-1 (1, a1 a2/E), <S, (1, a1 a2/E)> =
     (e1 e2 a1 a2 A B/E) S_0, S_0 = -((k + 1)/E) a_{k+1}^(k). The lower coefficients (j <= k) involve v_{2,k2} with k2 >= 1 and are determined recursively (Remark 8.9); no
     closed form in b is claimed.

### 11.4 Low orders, by hand (all re-derived above)
 - dphi = (a2 A/E) t - (a1 a2 A B/(2 E^3)) t^2 + O(t^3); alg = 1 - ((E^2 + u)/(2 E^2)) tau + O(tau^2); rho = 1 + (u/(2 E^2)) t + O(t^2); Delta_1 = -(e1 e2 E^3/t^2)(1 - u t/E^2 + O(t^2));
   Delta_2(0) = e1 e2 a1 a2 A B/E, Delta_2'(0) = 0; l_0 = E/2 - a0 - a1 a2/E; Theta_1 = -e1 e2 a1 a2 A B/E^2.
 - q^(1) = (u/2)(lam_1 - lam_2); G_1 = u/(2 (E + z)(E + 2z)); g_1 = (u/(2E))(2 lam_2 - lam_1) = +(u/2) lam_1 lam_2 (by (LAM5); sign corrected 2026-11-03); G_1 - rho_1 G_0 = [t^1] r_{0,1} = (E/(E + z))(H_1 - (E^2 + u)/(2E^2)).
 - q^(k) top and subleading: a_{k+1}^(k) = (-1)^(k+1) C_{k+1}, a_k^(k) = (-1)^(k+1) k C_k u/(2 E^2) (Corollary 11.3_CY); e.g. a_3^(2) = -C_3 = -P_2/(12 E^2) = -(3 a2 - E)(3 a2 - 2E)(3 A - E)(3 A - 2E)/(12 E^2),
   a_2^(2) = -2 C_2 u/(2 E^2) = -C_2 u/E^2 = u^2/(2 E^2), since C_2 = P_1/2 = (2 a2 - E)(2 A - E)/2 = (E - 2 a2)(2 a1 - E)/2 = -u/2; and a_1^(2) = -a_2^(2) - a_3^(2) by the sum rule. On T_F: a_2^(2) = (1 - x)^4/(2 c^2), a_3^(2) = -(2 - x)^2 (1 - 2x)^2/(12 c^2), [PR] s.11.4.
 - sigma^(0) = (E lam_1, a1 a2 lam_1); [t^0] Y = (0, E).
 - Loops: e_11^(0) = -e1 e2 a1 a2 A B/2; k = 1: e_12 = e_21 = -2 e1 e2 a1 a2 A B u/(3 E^2), e_22 = 0 (e_11^(1) needs Theta_2 = 3 e1 e2 E [t^3] dphi, not computed here).
 - E_11 legs at t^0: -(e1 e2 a1 a2 A B/2) E^(-b-1) lam_1. E_12 at t^0: lam_1(psi) chi12(w).
 Not re-derived: [PR] 11.4's t^2, t^3 tables, sigma^(1), Theta_2, Theta_3, the t^1 E_11 leg formula in b, and the explicit E_12 leg polynomials in b.

## 12. The three-torus a0 = 0 with a1, a2, e1 free, by rescaling

**Corollary 12.1 (a0 = 0).** Let a0 = 0 and a1 a2 E e1 e2 != 0 (then A = a2, B = a1, E = a1 + a2), and put x := a1/a2 (so x != 0, -1) and c := 1 + x = E/a2. Let the superscript
TF denote the frame-B objects of [PR] at parameter x. Then, for i = 1, 2 and all t:
    r_{0,i}(z; t) = r^TF_{0,i}(z/a2; t),   r_{1,i}(z; t) = a2 r^TF_{1,i}(z/a2; t),   L_i = a2 L_i^TF,   Delta_i = kappa Delta_i^TF,   kappa := e1 e2 a2^4/a1,
    <a, b> = kappa <a^TF, b^TF>_TF for a = (a_0, a2 a_1^TF), b likewise,   Dil_i(d) = a2^(1-d) Dil^TF_i(d),   E_ij(m, n) = kappa a2^(-1-m-n) E^TF_ij(m, n),
    rho = rho^TF,   lam_j(psi) = a2^-1 lam^TF_j(psi/a2),   q_t(psi) = a2 q^TF_t(psi/a2),   f_{j,b}(psi) = kappa a2^(-1-b) f^TF_{j,b}(psi/a2) (exponent corrected 2026-11-03, from E_1j(a,b) = kappa a2^(-1-a-b) E^TF(a,b)).
Consequently Theorem LC_CY at a0 = 0 is Theorem LC of [PR] at x = a1/a2 read through these rescalings; no new row analysis is needed.
Proof. (1) The superpotential: at a0 = 0, W = a1 log(1 - Q U) - a2 log U + a2 log(1 - U) = a2 W^TF with W^TF := x log(1 - Q U) - log U + log(1 - U) ([PR] (C1)); so W_k = a2 W_k^TF, the
critical points coincide (LeqU = a2 LeqU^TF), and by (D1) and (D4) (Lnorm[lambda f, 1](z) = Lnorm[f, 1](z/lambda), applied in the coordinate y = log U) r^pt_{0,i}(z) = r^{pt,TF}_{0,i}(z/a2).
(2) The twist: (C8') gives r_{0,i} = h_{(a1, a2, -E)}(z) r^pt_{0,i}, and h_{(a1, a2, -E)}(z) = h_{(x, 1, -c)}(z/a2) = h^TF_w(z/a2)^-1 (rescaling the multiset by a2 rescales z;
h_{(x, 1, -c)} = h_{(c, -x, -1)}^-1 by h_{-v} = h_v^-1), which is the T_F twist of [PR] (C8) at z/a2. So r_{0,i}(z) = r^TF_{0,i}(z/a2).
(3) L = -a1 Q U/(1 - Q U) = a2 L^TF. D_i is scale-invariant. r_{1,i}(z) = (a2 L_i^TF + z D_i) r^TF_{0,i}(z/a2) = a2 (L_i^TF + (z/a2) D_i) r^TF_{0,i}(z/a2) = a2 r^TF_{1,i}(z/a2).
(4) Pairing: at a0 = 0, <a, b> = e1 e2 [a1 a2 (a_0 b_1 + a_1 b_0) - E a_1 b_1]; with a_1 = a2 a_1^TF, b_1 = a2 b_1^TF and a1 a2 = x a2^2, E = c a2:
<a, b> = e1 e2 a2^3 [x (a_0 b_1^TF + a_1^TF b_0) - c a_1^TF b_1^TF] = (e1 e2 a2^3/x) <a^TF, b^TF>_TF, since <a, b>_TF = x^2 (a_0 b_1 + a_1 b_0) - c x a_1 b_1 ([PR] (C4)); e1 e2 a2^3/x = kappa.
Delta_i = e1 e2 (L_i - L_j)(a1 a2 - E L_i) = e1 e2 a2^3 (L_i^TF - L_j^TF)(x - c L_i^TF) = (e1 e2 a2^3/x) Delta_i^TF by [PR] (4.1), Delta^TF_i = x (L_i - L_j)(x - c L_i).
(5) Dil_i(d) = (-1)^d r_{0,i}[d-1] = a2^(1-d) Dil^TF_i(d). Lemma 1.1: P_ij(z, w) = -<v_i(-z), v_j(-w)> = kappa P^TF_ij(z/a2, w/a2), so E_ij(z, w) = [P_ij(z, w) - P_ij(-w, w)]/(z + w) =
(kappa/a2) E^TF_ij(z/a2, w/a2), i.e. E_ij(m, n) = kappa a2^(-1-m-n) E^TF_ij(m, n). (6) rho: Delta_1 t^2/(-e1 e2 E^3) = kappa Delta^TF_1 t^2/(-e1 e2 a2^3 c^3) = Delta_1^TF t^2/(-x c^3), so rho = rho^TF.
lam_j(psi) = 1/(a2 c - j psi) = a2^-1/(c - j psi/a2). q_t(psi) = -psi r_{0,1}(-psi) = a2 (-psi/a2) r^TF_{0,1}(-psi/a2) = a2 q^TF_t(psi/a2). f_{j,b}(psi) = [w^b] E_1j(psi, w) =
(kappa/a2) a2^-b [w'^b] E^TF_1j(psi/a2, w'). (7) Each statement of Theorem LC_CY is now the corresponding statement of [PR] Theorem LC at x: membership in Lam_{k+1} is
preserved by psi -> psi/a2 and the overall constants; F1, F2 (E/(E + z) = c/(c + z/a2), E/t + a1 a2/(E + z) = a2 [c/t + x/(c + z/a2)]), F3 are preserved. QED
Remarks. (i) [PR] Theorem LC needs x != 0, -1, i.e. a1 != 0 and E != 0; a2 != 0 for the rescaling; e1 e2 != 0 for kappa != 0: exactly (ND) at a0 = 0 (A = a2, B = a1).
(ii) The QRR constant at a0 = 0 is h_nw = h_{(a1, a2, -e1, -e2)}, not 1 unless a_i = e_i; it is the frame change h_w^-1 h_nw = h_{(a1, a2, -E)} that is e-independent. On T_F both
are trivial. (iii) The vertex class Omega_w, w = (E, -e1, -e2) = a2 (c, -e1/a2, -e2/a2), differs from the T_F class (c, -x, -1) unless e1 = a1; this concerns Theorem R
(Theorem L at tau_w = -e1/E), not Theorem LC_CY. (iv) The direct proof of sections 3-10 covers a0 = 0 as well (the point Q = 0 is never used there; h_nw enters only
through the multiset identity (3.9a)); the only statements of this text that need a0 != 0 are the citations of [PR] Lemma 13.3 for the definition of the constant
and of [PR] 13.2(b)(iv) for K_2 = 1, which extend to a0 = 0 by rationality (s.1.4 (P)) once the factors are defined there.
[Evidence: check2_rows.py: frame-B row at (0, a1, a2, e1, e2) equals the fibrewise frame-B row at x = a1/a2 with z -> z/a2 through tau^6 z^9 (tau^8 z^11), and
h_{(a1, a2, -E)}(z) = h_{(c, -x, -1)}(z/a2)^-1.]

## 13. Status

### 13.1 Status of every statement
P = proved here to all orders (complete argument in this text, or a verbatim citation of a weight-free statement of [PR]); P mod = proved modulo the named inputs;
E = evidence only.
 - s.1 (C1')-(C14'): definitions. (C1') is taken from [PR] Prop. 13.2(b)(i) and Lemma 13.3; that it is the code's generic superpotential is now confirmed by reading
   rc_lib.m (hodge_vertex_conifold/rc_lib.m, read by the hand referee and re-read for this revision): line 80 spotU = a1 Log[1 - Q U] - a2 Log[U] + (a0 + a2) Log[1 - U]
   is (C1') verbatim; line 115 LeqU = Q (a0 + a1) U^2 - (a0 + Q (a1 - a2)) U - a2 is (C2'); line 82 etaSym is (C4'); lines 106-109 ffcClosed =
   Exp[-Sum_k (B_2k/(2k (2k - 1))) z^(2k-1) (a1^(1-2k) + (a0 + a2)^(1-2k) - e1^(1-2k) - e2^(1-2k))] is h_nw, nw = (a1, A, -e1, -e2), in the sign convention of (C7');
   line 138 rmatrixAtClosed passes the same list ff to both roots; line 81 LofU = -(Q (a0 + a1) U + Q a2)/(1 - Q) equals L(U) of (C5') in QA (remark under (C5')),
   and line 128 evaluates Delta at L = LofU [P].
   (C7') uses [PR] Lemma 13.3 (h_nw, same at both roots) [P mod: Lemma 13.3 = Laplace + Stirling, at a0 != 0; a0 = 0 by rationality].
 - Lemmas 3.1_CY, 3.2_CY [P, elementary]. s.3.4 [P, verbatim]. Proposition 3.10_CY [P]. (3.16) [P]. Proposition 3.12_CY, Lemma 3.13_CY, Theorem 3.14_CY [P].
   Proposition 3.15_CY [P mod [PR] Proposition 3.15 (Laplace's method and Stirling's series, cited there)]. Corollary 3.16_CY, Theorem 3.17_CY, Corollary 3.18_CY (F2)
   [P mod Prop. 3.15_CY and Lemma 13.3]. Proposition 3.19_CY [P, formal; the identity R_tr = delta (log alg)_tau proved by hand].
 - Lemmas 4.1_CY, 4.2_CY, Propositions 4.3_CY (incl. Delta_2'(0) = 0), 4.4_CY, Lemma 4.5_CY [P, elementary]. Propositions 4.6_CY, 4.7_CY, 4.8_CY [P mod section 3].
 - Sections 5, 6, 7 [P mod sections 3-4]: L1, the sigma leaf structure, the flow, F3 for E_11, L3, L2 for j = 1.
 - Lemmas 8.1_CY-8.3_CY [P]. Theorem 8.4_CY [P; alternatively [PR] 13.2(b)(ii)]. Corollary 8.5_CY [P]. Proposition 8.6_CY(a) [P mod the method of [PR] Lemma 13.3 =
   Laplace + Stirling]; (b), (c) [P mod (a) for the values]. Corollary 8.7_CY [P mod 8.6(a)].
 - Section 9 [P mod sections 6, 8]. Section 10 [P mod the above]. Theorem LC_CY: P modulo [PR] Proposition 3.15 and [PR] Lemma 13.3 (both Laplace + Stirling), i.e.
   the same analytic inputs as [PR] Theorem LC (the generic superpotential (C1') is the code's, see the first item).
 - Theorem CF_CY: Theorem 11.1_CY, Corollaries 11.2_CY(a)-(e), 11.3_CY, Theorem 11.4_CY, Corollary 11.5_CY, s.11.3(a) [P mod Theorem LC_CY]; s.11.3(b), (c) and the E_22,
   Dil_2 values [P mod Proposition 8.6_CY(a)].
 - Corollary 12.1 (a0 = 0) [P mod [PR] Theorem LC].
 - Theorem Q_CY [P mod section 3 and, for (iii), 8.6(a)].
 - Remark 0.1 (no use of E = e1 + e2) [P: inspection of the proof; E for the failure of regularity off T_CY].
 - Corollary R_CY, Proposition 13.2 on T_CY, the rigidity corollary of the feasibility reports: not written here (s.0.4).
No step failed. No step uses a2 = A or a_i = e_i; the places where [PR] used them (1 - s = x s (1 - tau s); the double zeros (1 + j z)^2 of F_n; x^2 - c x L_2 = -x^3 tau s*^2;
h_nw = 1) are replaced by Lemma 3.1(b)(ii) in its two forms, the degree/leading-coefficient argument of Lemma 5.1_CY, a1 a2 - E L_2 = -B tau s* L_2, and the constant line (3.9a).

Referee record (2026-11-03). Two referee passes, both AI, both with no major item. (1) A hand referee re-derived s.1, s.3-s.10, Corollary 12.1 and Remark 0.1
line by line against [PR]; verdict: sound with minor fixes. Its items M1-M6 are applied in this revision: M1 the first expression for G_1 in the k = 0 computation
of Proposition 6.1_CY (s.6.2; the printed intermediate was [t^1] r_{0,1} = -u z/(E (E + z)(E + 2z)), the numerator and the value were right); M2 the sign of
2 B tau s* L_2/delta in the display of the proof of Proposition 3.19_CY (s.3.10; the final line was right); M3 the parenthesis on r_D in s.1.3 (the case a0 = 0,
E^2 + u = 0 is excluded by (ND)); M4 the form (3.4') of alg^2, which makes the analyticity of 1/alg on D in Proposition 3.19_CY manifest; M5 the code record above
(rc_lib.m read) and the remark on LofU under (C5'); M6 the one-line reason for the uniqueness of the class w in s.3.9. Every changed line was re-derived by hand
before the change. (2) An independent numerical referee with its own F_p code (p = 2^127 - 1 and 2^89 - 1; truncated Laurent series in t to t^14, z^23) tested the
KEY FORMULA, F1-F3, L1 (k <= 13), L2 (k <= 9), L3 (k <= 9), the flow (7.2), Theorem Q_CY, the regular-colour values (s.8), Corollary 12.1 at a0 = 0, and Remark 0.1
at non-CY points, at nine weight points including a1 = a2, u = 0 and e1 = e2; all passed, and four controls failed as required: the class (E + 1, -e1, -e2) in place
of w; no h_nw; A -> a2 in the rows; lam_j = 1/(E + j psi). It also observed that at u = 0 the top leaf lam_{k+1} is absent from [t^k](rho q_t) for odd k, and that at
E = 9 a2/7 the lam_9 coefficient of [t^8](rho q_t) vanishes, as the sharpness statement predicts (for odd k and E = 2 a2, resp. E = 2 a1, P_k contains the factor
(k + 1) a2 - ((k + 1)/2) E, resp. (k + 1) A - ((k + 1)/2) E; P_8 contains 9 a2 - 7 E). Both passes are evidence only, not proof steps. No human has read this text.

### 13.2 What was not re-derived
 - [PR] s.3.4 (formal Laplace calculus), Lemma 3.8, Lemma 3.9: cited as weight-free.
 - [PR] Proposition 3.15 and Lemma 13.3 (analytic inputs): cited; only their rescaling/multiset consequences are derived here.
 - [PR] Lemma 8.1 (conjugate encoding in the code): cited. That rc_lib.m's rmatrixAtClosed uses the same scalar at both roots is now read in the code (line 138,
   s.13.1); that this scalar is the Q = 0 constant h_nw is [PR] Prop. 13.2(b) and Lemma 13.3, not re-proved.
 - The two referee passes of s.13.1: their derivations and F_p code were not re-run here; only the fixes M1-M6 were re-derived.
 - [PR] Prop. 13.2(b)(ii), (iv): read, not re-proved; their scope (generic weights; (iv) needs a0 != 0 at Q = 0) is recorded in s.8.4, s.8.5.
 - [PR] s.11.2(c), (d) closed forms of e^(dphi), rho, alg without branch ambiguity; [PR] s.11.4's higher tables; sigma^(1), Theta_2, Theta_3, l_1.
 - Section 13 of [PR] (Theorem R, Corollary R, Proposition 13.2 at generic weights) and the rigidity corollary: not re-derived or restated beyond s.0.4.
 - The genus-0 potential PrepCl + Li_3 at generic weights is quoted from [PR] Lemma 4.1 / Prop. 13.2(a).

### 13.3 Finite-order checks run for this text (evidence only; fresh python with exact Fractions, no sympy; scripts and outputs in the scratch folder
cy_torus/lc_writer/: check1_identities.py, check1_out.txt, check2_rows.py, check2_out.txt, check2_out_T8Z11.txt)
 check1_identities.py (random rational weights on T_CY, a0 != 0, and random rational s, tau, z):
  - (3.16)_CY amplitude identity times (1 - tau s)^3, 30 points;
  - (4.4)_CY operator identity on tau^n, -3 <= n <= 6, 12 points;
  - at 25 rational points where delta^2 is a rational square, both roots: N(s) = 0; L = A/(tau s) - a0 from the superpotential; L1 + L2, L1 L2, tau (L1 - L2) = delta;
    a1 a2 - E L2 = -B tau s L2; (1 - tau) a2 s/(1 - tau s) = A - B s; the Hessian f''(s) = -delta/(s (1 - s) U); W_2(U) = U^2 f''/tau^2; Delta = e1e2 (L1 - L2)(a1 a2 - E L1)
    = a1 A e1 e2 W_2 = -e1 e2 E^3 alg^2/tau^2; the eikonal identity (1 - tau) Eik = B N(s); P_tr = -N'(s); <chi_1, chi_2> = 0, <chi_1, chi_1> = Delta_1; the three forms of Theta;
    Delta_2 = e1 e2 B delta s* L_2; self-adjointness of A and (4.8) on random covectors; L1 - dphi' = E/tau - A; tau p_1 = 1/s, tau p_Q = -(1 - tau)/(1 - s); s_tau = s L2/delta;
    (log alg^2)_tau = 2 R_tr/delta; the u identities; eta eta^-1 = 1; the multiset identity (3.9a) as power sums; the first-order values (u, l_0, a1 a2, s*'(0), Delta_2'(0) = 0, alg_1).
    All pass. (One line of the script, a mis-typed simplified form of R_tr with s^2 in place of s, failed and was corrected in the text; the identity R_tr = delta (log alg)_tau,
    computed from the definition of R_tr, passed.)
 check2_rows.py T Z (T = tau/t order, Z = z order), at the CY points (a0, a1, a2, e1) = (5/3, -4/7, 13/6, 9/11) and (-2/5, 7/4, 3, 1/6) with e2 = E - e1, run at (6, 9) and (8, 11):
  - the Wick row in y = log U, computed from the rescaled wr_k(tau) through the universal Wick polynomials (partitions of 2n into n - m parts >= 3), times h_{(a1, A, -E)}(z),
    equals (E/(E + z)) F Ephi alg (KEY FORMULA / Theorem 3.17_CY) coefficientwise, and the negative z-powers of the right side cancel;
  - r^pt_{0,1}(z; 0) = h_{(E, -a1, -A)} E/(E + z);
  - L1: [t^k](z G) prod_{j<=k+1} (E + j z) has z-degree exactly k for 1 <= k <= 6 (Z = 9), resp. k <= 8 (Z = 11), with G = rho r_{0,1} built from the Wick row;
  - the A_j closed form reproduces [t^k](z G) - E delta_k0 for k <= 5 (7); [t^1](z G) = (u/2)(1/(E + z) - 1/(E + 2z)); rho_1 = u/(2 E^2);
  - F2: rho r_{1,1} = E/t + a1 a2/(E + z) + O(t) for all z^m <= Z, with r_{1,1} built from (H) and Delta_1 = a1 A e1 e2 wr_2/tau^2;
  - Row 2 of the QDE at U_1, z^m <= Z, t^-3..t^(T-4);
  - K_1 = <v_1(z), v_1(-z)>/Delta_1 = 1, w^m <= Z, t^0..t^(T-4);
  - E_11(0, 0), E_11(1, 0) have no negative t-powers and the t^0 values -e1 e2 a1 a2 A B/(2 E^2), -e1 e2 a1 a2 A B/(2 E^3);
  - Theta(0) = 0, Theta_1 = -e1 e2 a1 a2 A B/E^2; the flow d/dt (t^2 rho^2 E_11) = Theta Dg Dg for m + n <= 3, t^0..t^(T-5);
  - L3 bidegree bound for k <= 2 on psi, psi' orders <= 4;
  - a0 = 0: the frame-B row at (0, 5/2, -3/7, 4/9, e2) equals the fibrewise frame-B row at x = a1/a2 with z -> z/a2, and h_{(a1, a2, -E)}(z) = h_{(c, -x, -1)}(z/a2)^-1.
    All pass at (6, 9) and at (8, 11) (check2_out.txt, check2_out_T8Z11.txt).
 Read but not re-run: cy_torus_explore/exact_test/t1_basics.py (QRR constants at Q = 0), critic/leaf_high.py (L1 to k = 10, L3 to k = 6 at another point), critic/misc.py.
