# S4. The E_12 legs, the regular idempotent, and closed forms (Steps 9, 10 and section 4 of outline.md)

2026-10-31, revision 1. This revision answers the referee report of 2026-10-31 (verdict "sound"); the changes are listed in S4.9. No
topological recursion and no spectral curves are used. Conventions: outline.md s.0 (rc_lib.m, hvc_lib.m; mode B, s = -1, eps1 = x,
eps2 = 1, c = 1 + x, alpha0 = 0; c != 0 and x != 0 throughout). No wolframscript was run for this section. The exact rational spot checks
(Perl + Math::BigRat) are leaf_proof/S4_checks.pl (output S4_checks_out.txt) and, for the parts added in revision 1,
leaf_proof/S4_checks_r1.pl (output S4_checks_r1_out.txt); see S4.5.

Summary.
(1) Step 9. L2 holds for j = 2 to all orders, and E_12 is regular at t = 0 (Prop. 9.3). The only property of the regular idempotent this
    uses is that its rows are analytic at t = 0. The same computation gives the leaf form of the whole edge (Prop. 9.4):
        rho E_12(psi, w) = -Sum_{j>=1} j <S_j(t), v_2(-w; t)> lam_j(psi) lam_j(-w),    S_j(t) = O(t^(j-1)), given explicitly.
    The graph sum also needs the transposition E_21(b, a) = E_12(a, b), because a singular vertex can be the second end of an edge. It
    follows from the orthogonality <v_1(w), v_2(-w)> = 0 (Prop. 9.5). Without the orthogonality the rationality argument of Prop. 9.3
    does not apply to the singular-side leg of an E_21 edge; with it, that leg is f_{2,b} and is rational. The orthogonality needs the QDE
    at the regular idempotent, which Lemma 10.5 proves.
(2) Step 10. The following are analytic at t = 0: all frame-B rows at U_2, L_2 = x s*, Delta_2 = x^3 delta s*^2 (with Delta_2(0) = x^3/c and
    Delta_2'(0) = 0), Dil_2 and E_22. Lemma 10.5 proves Row 2 of the QDE at U_2 to all orders in z: it is the argument of S1 Prop. 4.1 at
    the second critical point s_2 of the Euler phase, together with V1 at U_2 and S2's Lemma S2.6. With S2's Prop. S2.10 this gives the
    orthogonality (Corollary 10.5'). At t = 0, v_2(z; 0) = h_w(z)^-2 (1, x/c) (the Beta integral through S1's Lemma 2.5 and Stirling, and
    Row 2 at U_2). E_22 at t = 0 is (x^3/c)[1 - h_w(z)^2 h_w(w)^2]/(z + w), which is not of rank one. The rank-one value in
    conifold_regularity_note.md F3 is the E_11 table, read through the conjugate encoding of idempotent 2 in rc_lib.m (Remark 10.6.1).
    This completes F1 and F3 (see 10.7).
(3) Closed forms. Every dilaton leaf coefficient is an explicit function of t (Theorem CF.1):
        rho q_t(psi) = c + Sum_{j>=1} A_j(t) lam_j(psi),     A_j(t) = (-1)^j C_j t^-1 zeta(t)^j 2F1(j x/c, j x/c; j; tau),
        C_j = prod_{m=1}^{j-1} (j - m c)^2 / (j! (j-1)! c^(2j-4)),    zeta = tau e^(dphi(tau)/c) = c (1 - U_1) U_1^(-(1-x)/c).
    The following consequences hold for all k: the top-coefficient formula of leaf_structure.md (now proved) and the subleading coefficient;
    the E_11 corner formula (proved) and the whole anti-diagonal of L3, so the bounds in L1 and L3 are sharp; closed forms of rho, of Theta
    and of every loop coefficient e_il(t); the t = 0 form of E_12; and the top coefficients of all E_11 and E_12 legs.
(4) S4 predicts the following tables of leaf_structure.md, and all of them are reproduced exactly at x = 3/7, 2/5, 5/3, 7/2 (2592 exact
    comparisons; also at x = -1/3 and 11/5): rho to t^5, delta_k for k <= 5, all loop coefficients for k <= 5, all E_11-leg coefficients,
    and the t^0 and top E_12-leg coefficients. The parts rewritten or added in revision 1 (Lemma 10.5, Prop. 10.6, Remarks 10.6.1 and
    10.6'', CF.2(d)) are checked from the code's definitions by S4_checks_r1.pl (901 exact comparisons at x = 3/7, 5/3, -5/2, 2/9 and at three
    points with t != 0).


## S4.0 Notation and inputs

(N1) t-expansions. U_1(t) and U_2(t) are the two roots of Q x U^2 + Q (1-x) U - 1 = 0. At Q = 1 the equation is (x U + 1)(U - 1) = 0, with the
     simple roots 1 and -1/x (distinct because c != 0). By the implicit function theorem both roots are analytic functions of t near 0
     (Q = e^t). For k >= 1, W_k := (U d/dU)^k W is a polynomial in (1-U)^-1 and (1-QU)^-1 with coefficients in Z[x] (see 10.1). The rows r_{a,i}[m],
     L_i, Delta_i, rho, Dil_i and E_ij are built from these, from W_2^-1, and from logarithmic t-derivatives of Delta_i. So every coefficient
     (at fixed powers of z, w, psi) is meromorphic at t = 0. At U_1 the factors (1 - U_1)^-1 = (tau s*)^-1 and
     (1 - Q U_1)^-1 = ((e^t - 1)(s* - 1))^-1 have simple poles; this is where the singular behaviour of idempotent 1 comes from. [t^k] X
     denotes the coefficient of t^k in the Laurent expansion, taken coefficientwise in the formal variables. For a generating function we
     write X = Sum_k t^k X^(k) with X^(k) in C[[psi, w]].
(N2) Leaf spaces. As spaces of rational functions of psi, put V_m := span(lam_1, ..., lam_m) (V_0 = 0) and V_m^+ := C + V_m. By partial
     fractions, a rational function f lies in V_m^+ iff its poles are simple and lie in {c/j : 1 <= j <= m} and f is bounded at infinity.
     It lies in V_m iff in addition f(psi) -> 0 as psi -> oo. A rational function regular at psi = 0 is identified with its Taylor series
     at 0 (this map is injective). In particular lam_j(psi) = Sum_{a>=0} j^a psi^a / c^(a+1).
(N3) Three elementary identities, used repeatedly:
        psi lam_j = (c lam_j - 1)/j,    lam_j(psi) - lam_j(-w) = j (psi + w) lam_j(psi) lam_j(-w),    (c - m psi) lam_j = m/j + c (1 - m/j) lam_j.
(N4) <a, b> = x^2 (a_0 b_1 + a_1 b_0) - c x a_1 b_1 (outline s.0). Two special covectors: for every covector b,
        <(0, 1), b> = x (x b_0 - c b_1),    <(1, x/c), b> = (x^3/c) b_0.
     A covector proportional to (1, x/c) therefore only sees the identity component of the other argument.

Inputs from S1-S3, exactly as stated in outline.md:
(I1) [Steps 2-3, V3] G := rho r_{0,1} = (c/(c+z)) (tau/t) F(z, tau) Exp[-dphi(tau)/z] as a formal series in t. Each t-coefficient is a
     rational function of z, regular at z = 0 (Step 2(d)), identified with its Taylor series in z.
(I2) [Step 4] For k >= 0, [t^k](z G) = c delta_{k0} + Sum_{j=1}^{k+1} a_j^(k)/(c + j z) for some numbers a_j^(k), with a_1^(0) = -c^2. Since
     z G vanishes at z = 0, Sum_{j} a_j^(k) = -c^2 delta_{k0} (the "sum rule"). In psi = -z this reads
     rho q_t(psi) = -psi g(psi) = c + Sum_{k,j} t^k a_j^(k) lam_j(psi). We set a_j^(m) := 0 for m < j - 1.
(I3) [Step 3, F1] rho is analytic at t = 0 with rho(0) = 1, so rho^-1 and rho^-2 are power series in t with constant term 1. Also
     rho^2 Delta_1 = -x c^3/t^2.
(I4) [Step 7] sigma(psi) := rho v_1(-psi) = (g, (L_1 - psi/t) g - psi dg/dt), with d/dt at fixed psi, and L_1 = c/tau - x s* = c/t + l(t)
     with l analytic. Since c/tau = c/t + c/2 + O(t), l(0) = c/2 - x/c = (1 + x^2)/(2c).
(I5) [Step 6; S2 Prop. S2.10 with Lemma 10.5] Orthogonality: <v_1(w), v_2(-w)> = 0. It is proved in 10.5 (Corollary 10.5'). It is used ONLY
     in Prop. 9.5 (transposition), in Remark 10.6' and in the last sentence of Remark 10.6'' (V6 for i != j). 10.6(b),(c) and CF.6(c) need only
     its t^-1 part, which also follows from Lemma 10.5 directly (10.6(b), first proof).
(I6) [Step 8] d/dt [t^2 rho^2 E_11(psi, psi')] = Theta(t) d/dt(t g(psi)) d/dt(t g(psi')) with Theta = -x^3 tau s*^2, and t^2 rho^2 E_11 = O(t^2).
     This gives L3, L2 for j = 1 and the regularity of E_11 at t = 0; E_11(a, b) = E_11(b, a) follows from K_1 = 1 (Step 6).
(I7) [V2; S1 s.0, S2 (0.2), Lemma S2.2, Prop. S2.3] u := 4 x tau/c^2 and delta := c (1 - u)^(1/2), with the principal root of 1 - u (which
     is near 1). So delta^2 = c^2 - 4 x tau and delta(0) = c; for c < 0 this is not the principal root of c^2 - 4 x tau.
     x tau s*^2 - c s* + 1 = 0, s* = 2/(c + delta), hence x tau s* = (c - delta)/2. Also L_2 = x s*, L_1 = c/tau - x s*,
     Delta_i = x^3 (L_i - L_j)/(x - c L_j) for i != j, and Delta_1 = -x delta/(tau^2 s*^2).
Inputs used only in 10.5 and 10.6(a):
(I8) [S1 s.1, s.2, s.4, s.6] The Euler coordinate: for tau != 0, W(1 - tau s) = phi_tau(s) + c log tau + x t + x log(-1) (S1 (1.1)), with
     phi_tau(s) = log s + x log(1 - s) - log(1 - tau s). The formal Laplace calculus of S1 s.2: the Gaussian-normalised expansion
     Lnorm[f, g; u0] at a nondegenerate critical point u0 and Lhat := Lnorm/g(u0), with the Wick form (D1), the regularity statement (D2) and
     the invariances (D3); P1 = S1 Prop. 2.4 (coordinate invariance; used here in the normalised form (iii), which involves no choice of
     square root); P2 = S1 Prop. 2.3 (derivative in a parameter with a moving expansion point; used in the normalised form (2.1));
     P3 = S1 Prop. 2.2 (Lnorm kills total derivatives); and the amplitude identity S1 (4.1) (Euler's integration by parts in amplitude
     form, an identity of rational functions of (s, tau)). S1 states P1-P3 for an arbitrary nondegenerate critical point depending
     analytically on a parameter. Analytic input: S1 Lemma 2.5 (Laplace's method for real data) and S1 Prop. 6.1(1), (3), (4) (rationality
     in x, Stirling's series with remainder, uniqueness of asymptotic power series).
(I9) [S2] Lemma S2.2 (on the critical locus L(U) = (x U + 1)/tau = 1/(1 - U); L_1 = 1/(tau s*), L_2 = x s* = 1/(tau s~) with s~ = 1/(x tau s*));
     Prop. S2.3(c), (d) (Delta_i = x^2 W_2(U_i) for i = 1 AND i = 2, proved by hand: this is V1 at both roots; Delta_2 = x^3 delta s*^2);
     Prop. S2.4(a) (W_2(U_i) = (1 - tau s_i)^2 phi''_tau(s_i)/tau^2 at both critical points); Lemma S2.6 (= V7, proved by hand as an identity
     of differential operators); Prop. S2.7 (the QDE at U_1, z D_1 v_1 = (A - L_1) v_1, from S1's closed form); Prop. S2.8(c) (the matrix
     A = [[0, 1], [-x/tau, c/tau]] is self-adjoint for < , >); Prop. S2.10 (the QDE at both idempotents implies the orthogonality).
     M := C({t}) is the field of meromorphic germs at t = 0 (S2 (0.2)); the rows lie in M[[z]]^2.


## S4.1 The edge factors as the code computes them

In evalQ (rc_lib.m s.3):
    pre[[p+1, q+1]] = Sum_{a,b=0,1} -(-1)^(p+q) etainvCur[[a+1, b+1]] rN[a,i][p] rN[b,j][q],
    EdN[i,j][m, n]  = Sum_{k=0}^{n} (-1)^k pre[[m+k+2, n-k+1]]          (m + n <= ZMAX - 1).
Here etainvCur = eta^-1 is t-independent. With v_i[p] := (r_{0,i}[p], r_{1,i}[p]) this reads
    E_ij(m, n) = Sum_{k=0}^{n} (-1)^k P^ij_{m+1+k, n-k},    P^ij_{pq} := -(-1)^(p+q) <v_i[p], v_j[q]>.                           (1.1)
Only entries with p + q = m + n + 1 <= ZMAX occur, so (1.1) is exact for every (m, n) the code computes. In generating functions,
    P_ij(z, w) := Sum_{p,q} P^ij_{pq} z^p w^q = -<v_i(-z), v_j(-w)>,
    E_ij(z, w) := Sum_{m,n} E_ij(m, n) z^m w^n = [P_ij(z, w) - P_ij(-w, w)] / (z + w).                                         (1.2)
Proof of (1.2). For p >= 1, (z^p - (-w)^p)/(z + w) = Sum_{l=0}^{p-1} z^l (-w)^(p-1-l); for p = 0 the quotient is 0. So the right side of (1.2)
is Sum_{p,q} P_pq w^q Sum_{l=0}^{p-1} z^l (-w)^(p-1-l). Its coefficient of z^m w^n collects l = m, p = m+1+k and q = n-k with k = p-1-m in
[0, n], with sign (-1)^k. This is (1.1). QED
The numerator of (1.2) vanishes on z = -w by construction, so the quotient lies in C[[z, w]] (coefficientwise in t) for ANY pair of rows.
Thus (1.2) uses no symplecticity. The legs are f_{j,b}(psi) = [w^b] E_1j(psi, w) and the dilaton leaves are Dil_i(n) = (-1)^n r_{0,i}[n-1] (DilN).
(The stored rows of idempotent 2 are in conjugate encoding, and all factors are evaluated at U_1; see 10.2 and Remark 10.6.1.)

What symplecticity (V6: <v_i(w), v_j(-w)> = delta_ij Delta_i; proved for all i, j by S2 Prop. S2.9(c), Corollary 10.5' and Remark 10.6'')
adds:
 (S-a) The subtracted term becomes P_ij(-w, w) = -delta_ij Delta_i. This gives the usual form E_ij = [delta_ij Delta_i - <v_i(-z), v_j(-w)>]/(z + w).
 (S-b) Transposition E_ji(z, w) = E_ij(w, z) (Prop. 9.5). For i != j only the orthogonality (I5) is needed.
 (S-c) Nothing for L2 with j = 2 or for the regularity of E_12 (Prop. 9.3). There the subtracted term P_12(-w, w) enters only as a function
       of w alone, and its value is irrelevant.


## S4.2 Step 9: the E_12 legs

Lemma 9.1 (sigma order by order). Write sigma = Sum_{k>=-1} t^k sigma_k. Then sigma_{-1} = (0, c). For k >= 0 both components of sigma_k
lie in V_{k+1}; explicitly
    g_k = -(1/c) Sum_{j=1}^{k+1} j a_j^(k) lam_j,
    sigma_k = Sum_{j=1}^{k+1} S_j^(k) lam_j,    S_j^(k) := ( -(j/c) a_j^(k),  (k+2-j) a_j^(k+1) - (j/c) Sum_{i=0}^{k} l_i a_j^(k-i) ),         (9.1)
where l(t) = Sum_i l_i t^i.
Proof. By (I2), the sum rule and (N3),
c delta_{k0} + Sum_j a_j^(k) lam_j = Sum_j a_j^(k) (lam_j - 1/c) = (psi/c) Sum_j j a_j^(k) lam_j.
Dividing -psi g_k = c delta_{k0} + Sum_j a_j^(k) lam_j by -psi gives g_k, with g_0 = c lam_1. Write sigma^(1) for the second component.
By (I4) with L_1 = c/t + l,
    [t^k] sigma^(1) = c g_{k+1} - psi g_{k+1} - (k+1) psi g_{k+1} + Sum_{i=0}^{k} l_i g_{k-i} = (c - (k+2) psi) g_{k+1} + Sum_{i=0}^{k} l_i g_{k-i}.
For k = -1 this is (c - psi) c lam_1 = c, while the first component has no t^-1 part. For k >= 0, (N3) gives
    (c - (k+2) psi) g_{k+1} = -(1/c) Sum_{j=1}^{k+2} j a_j^(k+1) [ (k+2)/j + c (1 - (k+2)/j) lam_j ]
                            = -((k+2)/c) Sum_j a_j^(k+1) + Sum_{j=1}^{k+2} (k+2-j) a_j^(k+1) lam_j.
The constant vanishes by the sum rule (k+1 >= 1), and the term j = k+2 carries the factor 0. Collecting terms gives (9.1). QED
So for k >= 0, sigma_k has no constant term and has poles only at c/j with j <= k+1. The factor c - (k+2) psi removes the pole at c/(k+2)
that g_{k+1} carries. (This is Proposition 4 of S3_leaves_loops.md; the explicit coefficient vectors S_j^(k) are added here. They agree
with S3's sigma^(1).) Check: with a_1^(1) = -(1-x)^2/2 and l(0) = (1+x^2)/(2c), sigma_0 = (c, x) lam_1. In z this is F2:
rho r_{1,1}(z) = c/t + x/(c+z) + O(t). The t^-1 part of sigma is the psi-independent covector (0, c).

Lemma 9.2 ([w^b] extraction). Let E be in C[[psi, w]] and Num := (psi + w) E. For every b >= 0, in C((psi)),
    [w^b] E = Sum_{n=0}^{b} (-1)^n psi^(-n-1) [w^(b-n)] Num.                                                                  (9.2)
Proof. Put E_m := [w^m] E in C[[psi]] and E_{-1} := 0. Then [w^m] Num = psi E_m + E_{m-1}, so
Sum_{n=0}^{b} (-1)^n psi^(-n-1) (psi E_{b-n} + E_{b-n-1}) = Sum_{n=0}^{b} (-1)^n psi^(-n) E_{b-n} - Sum_{n=1}^{b+1} (-1)^n psi^(-n) E_{b-n} = E_b. QED
Equivalently, (9.2) is the expansion of Num/(psi + w) in C((psi))[[w]], i.e. with 1/(psi + w) = Sum_n (-1)^n w^n psi^(-n-1). This is legitimate
because C[[psi, w]] is contained in C((psi))[[w]], where psi + w is invertible. Being linear, (9.2) holds coefficientwise in t.

Proposition 9.3 (L2 for j = 2; regularity of E_12). Assume (I2)-(I4) (through Lemma 9.1) and
    (R) for every m, the regular row v_2[m] = (r_{0,2}[m], r_{1,2}[m]) is analytic at t = 0 (proved in 10.2).
Then:
 (i) rho E_12 and E_12 are O(t^0), i.e. every E_12(a, b) is regular at t = 0 (F3 for E_12);
 (ii) for all k >= 0 and b >= 0, [t^k] (rho f_{2,b}) lies in V_{k+1}.
Proof. Put N(psi, w) := rho P_12(psi, w) = -<sigma(psi), v_2(-w)>. Since rho depends on t only, (1.2) gives
    (psi + w) rho E_12(psi, w) = N(psi, w) - N(-w, w).                                                                       (9.3)
By Lemma 9.1 and (R), N is a Laurent series in t that starts at t^-1, and each coefficient is a finite sum:
    [t^k] N(psi, w) = -Sum_{k1=-1}^{k} <sigma_{k1}(psi), v_{2,k-k1}(-w)>,    v_{2,k2}(-w) := [t^(k2)] v_2(-w), a pair in C[[w]].     (9.4)
(i) [t^-1] N(psi, w) = -<(0, c), v_{2,0}(-w)> does not depend on psi, so the t^-1 coefficient of the right side of (9.3) is zero. Since C[[psi, w]]
has no zero divisors, [t^-1](rho E_12) = 0, and lower orders do not occur. So rho E_12 = O(t^0), and E_12 = rho^-1 (rho E_12) = O(t^0) by (I3).
(ii) Fix k >= 0 and b >= 0. Apply (9.2) to E = [t^k](rho E_12) and Num = [t^k](N(psi, w) - N(-w, w)):
    [t^k] rho f_{2,b}(psi) = Sum_{n=0}^{b} (-1)^n psi^(-n-1) Y_{b-n}(psi),                                                    (9.5)
    Y_m(psi) := [w^m] [t^k] (N(psi, w) - N(-w, w)) = -(-1)^m Sum_{k1=-1}^{k} <sigma_{k1}(psi), v_{2,k-k1}[m]> - y_m.
Here v_{2,k2}[m] := [t^(k2)] v_2[m] is in C^2 (note [w^m] v_2(-w) = (-1)^m v_2[m]), and y_m := [w^m] [t^k] N(-w, w) is a constant. By Lemma 9.1
both components of every sigma_{k1} (-1 <= k1 <= k) lie in V_{k1+1}^+, which is contained in V_{k+1}^+. Hence Y_m is the Taylor series of a
rational function Yhat_m in V_{k+1}^+. Consider the rational function
    R(psi) := Sum_{n=0}^{b} (-1)^n psi^(-n-1) Yhat_{b-n}(psi).
 (a) Away from psi = 0, R has at most simple poles, located at the points c/j with 1 <= j <= k+1 (psi^(-n-1) is regular there).
 (b) R(psi) -> 0 as psi -> oo, because each Yhat is bounded at oo and psi^(-n-1) -> 0.
 (c) The Laurent expansion of R at psi = 0 is the right side of (9.5) (expand each Yhat_{b-n} by its Taylor series). So it equals
     [t^k] rho f_{2,b}, which lies in C[[psi]] because E_12(psi, w) is a power series in psi (1.2). Hence the principal part of R at 0
     vanishes: the apparent pole of order <= b+1 at psi = 0 cancels, and R is regular at 0.
By (a)-(c), consider R - Sum_{j=1}^{k+1} R_j lam_j with R_j := lim_{psi -> c/j} (c - j psi) R(psi). It has no poles in C, so it is a
polynomial; it tends to 0 at infinity, so it is 0. Thus [t^k] rho f_{2,b} = R lies in V_{k+1}. QED
The lam_j coefficients of R are explicit. Since psi^(-n-1) lam_j = (j/c)^(n+1) lam_j + (a polynomial in 1/psi without constant term),
    R_j = Sum_{n=0}^{b} (-1)^n (j/c)^(n+1) [lam_j] Yhat_{b-n}.

What is used about the regular rows. Only (R). The rows v_2 enter (9.5) through the constants v_{2,k2}[m] and y_m. Their z-dependence, their
values, symplecticity and the orthogonality (I5) play no role. Analyticity is exactly what the degree bound needs: if some v_2[m] had a pole of
order p at t = 0, (9.4) would bring sigma_{k+p} into the order t^k, i.e. lam_j with j up to k+p+1.

Proposition 9.4 (leaf form of the E_12 edge; second proof). Put S_j(t) := Sum_{k>=j-1} t^k S_j^(k). This is a formal series with
S_j = O(t^(j-1)); a closed form is given in CF.2(f). Then
    sigma(psi) = (0, c/t) + Sum_{j>=1} S_j(t) lam_j(psi),                                                                      (9.6)
    rho E_12(psi, w) = -Sum_{j>=1} j <S_j(t), v_2(-w; t)> lam_j(psi) lam_j(-w),     lam_j(-w) = 1/(c + j w).                    (9.7)
These are identities in C[[psi, w]]((t)), and at order t^k only j <= k+1 contribute (by (R)). Consequently, for all k and b,
    [t^k] rho f_{2,b}(psi) = -Sum_{j=1}^{k+1} j lam_j(psi) Sum_{k1=j-1}^{k} < S_j^(k1), [w^b] ( v_{2,k-k1}(-w) / (c + j w) ) >.        (9.8)
Proof. (9.6) is (9.1). Insert it in N. The term (0, c/t) contributes the same function of w to N(psi, w) and to N(-w, w), so by (N3)
    N(psi, w) - N(-w, w) = -Sum_j <S_j, v_2(-w)> (lam_j(psi) - lam_j(-w)) = -(psi + w) Sum_j j <S_j, v_2(-w)> lam_j(psi) lam_j(-w).
Dividing by psi + w (9.3) gives (9.7), and (9.8) is its coefficient of t^k w^b. QED
In (9.7) the cancellation at psi = 0 of Prop. 9.3(c) is explicit, since psi + w divides each lam_j(psi) - lam_j(-w). The subtracted term
N(-w, w), which symplecticity would set to zero, drops out without being evaluated. Prop. 9.3(ii) is the statement "j <= k+1" in (9.8). The
explicit form is used for the closed forms in CF.6.

Proposition 9.5 (transposition). Assume (I5). Then E_21(z, w) = E_12(w, z), i.e. E_21(b, a) = E_12(a, b) for all a, b. Hence E_21 is regular at
t = 0. An edge joining a singular vertex to a regular one therefore gives the singular vertex the leg f_{2,b} (b the power at the regular end)
whichever end comes first in the code's edge list. (graphValueC attaches EdN[c(v1), c(v2)][a_{v1}, a_{v2}] to the edge (v1, v2).)
Proof. phi(u) := P_12(-u, u) = -<v_1(u), v_2(-u)> = 0 by (I5), so E_12(z, w) = P_12(z, w)/(z + w). By the symmetry of < , >,
P_21(z, w) = -<v_2(-z), v_1(-w)> = P_12(w, z) and P_21(-w, w) = P_12(w, -w) = phi(-w) = 0. Hence E_21(z, w) = P_12(w, z)/(z + w) = E_12(w, z). QED
What (I5) is needed for. Without (I5) one only has E_21(z, w) - E_12(w, z) = [phi(z) - phi(-w)]/(z + w). By (9.2) the singular-side leg of an
E_21 edge is then
    [z^b] E_21(z, psi) = rho^-1 [ Sum_{n=0}^{b} (-1)^n psi^(-n-1) [z^(b-n)] N(psi, z) - (-1)^b psi^(-b-1) N(psi, -psi) ].
This contains N(psi, -psi) = -<sigma(psi), v_2(psi)>, which pairs sigma with the full power series v_2(psi). At t = 0 that series is
h_w(psi)^-2 (1, x/c) (Prop. 10.6), a divergent Stirling series, so the rationality argument of Prop. 9.3 does not apply to this expression. With
(I5) the leg is f_{2,b} by the transposition, and it is rational by Prop. 9.3.
Even the regularity of E_21 needs the order-t^-1 part of (I5). Since rho P_21(z, w) = -<v_2(-z), sigma(w)> and sigma_{-1} = (0, c), (N4) gives
    [t^-1] rho E_21(z, w) = [f(z) - f(-w)]/(z + w),     f(z) := -c x (x r_{0,2}(-z; 0) - c r_{1,2}(-z; 0)),
and no lower powers of t occur (the rows v_2 are analytic, 10.2). This vanishes iff f(z) = f(-w) in C[[z, w]], i.e. iff f is constant. The
constant is f(0) = -c x (x r_{0,2}[0](0) - c r_{1,2}[0](0)) = -c x (x - c L_2(0)) = 0, because r_{0,2}[0] = 1, r_{1,2}[0] = L_2 r_{0,2}[0] = L_2
and L_2(0) = x/c. So E_21 is regular iff f = 0, i.e. iff x r_{0,2}(z; 0) = c r_{1,2}(z; 0). That is 10.6(b), the t^-1 part of (I5); it also
follows from Row 2 at U_2 alone (10.6(b), first proof). So for E_21 the second QDE row at the regular idempotent is genuinely needed: for the
regularity of E_21 directly, and for the transposition through (I5). Both are proved: Lemma 10.5 and Corollary 10.5' (S2 Prop. S2.10).

Corollary 9.6 (unrescaled ends; the legs in Theorem R (ii)). For all k, b >= 0 and j in {1, 2}:
    [t^k] f_{j,b} and [t^k] (rho f_{j,b}) lie in V_{k+1};
    [t^k] (rho^e E_11(psi, psi')) lies in span{ lam_i(psi) lam_l(psi') : i, l >= 1, (i-1) + (l-1) <= k } for e = 0, 1, 2.
All these series are O(t^0).
Proof. For j = 2 this is Prop. 9.3. For j = 1, rho f_{1,b} = rho^-1 [psi'^b](rho^2 E_11) and [psi'^b] lam_l(psi') = l^b/c^(b+1). So by L3 (I6),
[t^k] rho f_{1,b} is a combination of lam_i(psi) with i - 1 <= k. The remaining cases follow on multiplying by rho^-1 or rho^-2 (I3): for u in C[[t]],
[t^k](u X) = Sum_{i=0}^{k} u_i [t^(k-i)] X, and the spans grow with k. QED

Use in Theorem R (proof_status_mixed_scale.md s.3, step (ii)). At a singular vertex v of genus h_v >= 2, homogeneity rescales every half-edge
at v exactly once (one factor rho per half-edge at v). The factors attached to v depend on the other end of the edge:
 - the other end is a regular vertex: the edge factor E_12 (or E_21) gives v the leg rho f_{2,b}, b the power at the regular end (Prop. 9.3;
   either orientation by Prop. 9.5);
 - the other end is a singular vertex v' of genus 0 or 1, which is not rescaled: the edge factor is rho E_11 = rho^-1 (rho^2 E_11), and for
   each power b' at v' it gives v the leg rho f_{1,b'} (L2 for j = 1, Cor. 9.6), with D <= k at the v-end; either orientation, since E_11
   is symmetric;
 - the other end is a singular vertex v' of genus >= 2, which is also rescaled: the edge factor is rho^2 E_11, which splits its D <= k
   between the two ends (L3);
 - a self-loop at v: the factor is rho^2 E_11 (L3).
At a vertex v' of genus 0 or 1 itself only regularity is needed: Delta_1^(h-1) is O(t^2) for h = 0 and equals 1 for h = 1, and every factor
attached to v' is O(t^0) (F3, 10.7). This is precisely the parenthetical remark in Theorem R (ii); nothing else about the legs is used there.


## S4.3 Step 10: the regular idempotent

10.1 The critical point U_2 and the W_k. Put theta := U d/dU. Then theta log U = 1, theta log(1-U) = 1 - (1-U)^-1, theta log(1-QU) = 1 - (1-QU)^-1,
and theta (1-U)^-n = n ((1-U)^-(n+1) - (1-U)^-n), and likewise with QU in place of U. Hence W_1 = x - x (1-QU)^-1 - (1-U)^-1, and for k >= 1, W_k
is a polynomial in (1-U)^-1 and (1-QU)^-1 with coefficients in Z[x]. In particular W_k is regular wherever U is neither 1 nor 1/Q. By (N1),
U_2(t) is analytic at t = 0 with U_2(0) = -1/x; explicitly U_2 = 1 - tau s_2 = 1 - (c + delta)/(2x), where s_2 = (c + delta)/(2 x tau) is the
second root of x tau s^2 - c s + 1. Since -1/x != 1 = 1/Q(0), every W_k(U_2(t)) is analytic at t = 0. Moreover W_2 = -x Q U/(1-QU)^2 - U/(1-U)^2,
so at t = 0 its value is -c U/(1-U)^2 at U = -1/x:
    W_2(U_2(0)) = x/c != 0.
For all small t, W_2(U_2) = x delta s*^2 (S2 Prop. S2.3(c), (d); (10.13) below), which does not vanish.

10.2 The rows; proof of (R). rmatrixAt0 (alpha0 = 0, hence no QRR constants) sets ratCN[0,2][m] = cconj[wickCN[m]]. This is the algebra element
whose value at U_1 (evC) equals the value at U_2 of wickZ[[m+1]] with Wk[k] -> W_k. By wickDirect, wickZ[[m+1]] is a polynomial in
W_2^-1, W_3, ..., W_{2m+2}, so by 10.1 each point-frame row r^pt_{0,2}[m] is analytic at t = 0. Frame B (twistRows[hCoeffs[w, -1, ZMAX]],
w = (c, -x, -1)) multiplies each row by the t-independent series h_w(z)^-1, which preserves analyticity. The H-row of rmatrixCore is the same
formula for i = 1 and i = 2:
    r_{1,2}[m] = L_2 r_{0,2}[m] + d/dt r_{0,2}[m-1] - r_{0,2}[m-1] Delta_2'/(2 Delta_2),      d/dt = Q d/dQ.
Here cDerN differentiates the algebra element including the motion of the root (C[a, b] -> C[a', b'] + b dU/dQ). Applied to a cconj'd element,
whose value at U_1 is a function evaluated at U_2, it therefore returns the total derivative along U_2(t). The formula is linear over C[[z]]
with t-independent coefficients, so it commutes with the twist and holds in frame B. With L_2 = x s* (analytic, L_2(0) = x/c) and 10.3,
every r_{1,2}[m] is analytic at t = 0. This proves (R). Note r_{0,2}[0] = 1 and r_{1,2}[0] = L_2.
Encoding. rc_lib.m computes in the quadratic algebra QA := Q(x, Q)[U]/(LeqU) and stores every idempotent-2 quantity in conjugate encoding:
the stored element is the one whose value at U_1 is the idempotent-2 value. All factors, of either idempotent, are therefore to be evaluated
at U_1. Remark 10.6.1 shows what evaluation at U_2 produces instead.

10.3 Delta_2. DeltaCN[2] = cconj[DeltaCN[1]] is DeltaL at L = L_2, and DeltaL = 1/eta(e_1, e_1). Indeed, in rc_lib.m s.1,
etacan = Psi^-T eta Psi^-1 with Psi = {{1, L1}, {1, L2}}, and the first column of Psi^-1 is the idempotent e_1 = (H - L_2)/(L_1 - L_2) in the flat
basis. By the Frobenius property,
    eta(e_i, e_i) = eta(e_i, 1) = [eta(H, 1) - L_j eta(1, 1)]/(L_i - L_j) = (x - c L_j)/(x^3 (L_i - L_j)),
which is the Delta formula in (I7). The remaining identities of (I7) are elementary: L(U) = -Q (xU + 1)/(1 - Q) = (xU + 1)/tau because
tau = 1 - 1/Q; with U_1 = 1 - tau s*, U_2 = 1 - tau s_2 and s* + s_2 = c/(x tau) this gives L_1 = c/tau - x s* and L_2 = x s*. Using (I7),
2 x tau s* - c = -delta and x tau - c^2 + c x tau s* = -(c + delta)^2/4, so
    Delta_2 = x^3 (L_2 - L_1)/(x - c L_1) = x^3 (2 x tau s* - c)/(x tau - c^2 + c x tau s*) = 4 x^3 delta/(c + delta)^2 = x^3 delta s*^2.     (10.1)
delta is analytic at tau = 0 with delta(0) = c != 0. Hence Delta_2 is analytic at t = 0, and
    Delta_2(0) = x^3/c,    Delta_2'(0) = 0    (because delta s*^2 = 1/c - x^2 tau^2/c^5 + O(tau^3)).
For comparison, Delta_1 = -x delta/(tau^2 s*^2), so Delta_1 Delta_2 = -x^4 delta^2/tau^2.

10.4 Dil_2 and E_22. Dil_2(n) = (-1)^n r_{0,2}[n-1] (DilN). By (1.1), E_22(m, n) is a polynomial with constant coefficients in the entries of
v_2[p], p <= m+n+1. Both are therefore analytic at t = 0 by 10.2.

10.5 The QDE at U_2 and the orthogonality (I5).
Step 5 of the outline is the QDE at idempotent 1 (S2 Prop. S2.7, from S1's closed form). Step 6 needs it at both idempotents: S2 Prop. S2.10
derives <v_1(w), v_2(-w)> = 0 from the QDE at U_1 and at U_2. There is no closed form at U_2, and none is needed. The proof of S1 Prop. 4.1
(Gauss's equation for the formal expansion of the Euler integrand) works verbatim at the second critical point s_2 of phi_tau, for tau near any
fixed tau0 != 0. The passage to tau-series at tau0 = 0 (the map iota) and the uniqueness Lemma 5.1, which S1 needs at s*, are not used.

Standing assumptions for 10.5. Put r0 := min(1, |c|^2/(4|x|)) and let 0 < |tau| < r0 (S1's domain D' minus 0). There delta is analytic and
non-zero with c + delta != 0 (S1 1.2), and
    s_2 := (c + delta)/(2 x tau) = 1/(x tau s*),    U_2 = 1 - tau s_2 = 1 - (c + delta)/(2x),    f2 := phi''_tau(s_2) = delta/(x s_2^2 U_2^2).
s_2 is a root of x tau s^2 - c s + 1 = 0, and s_2 is not in {0, 1, 1/tau}: s_2 = 0 contradicts the constant term 1, and s_2 = 1 or s_2 = 1/tau
forces tau = 1. Hence U_2 is not in {0, 1, 1/Q} (1 - Q U_2 = -Q tau (1 - s_2), S1 1.1), and f2 != 0. The formula for f2 is the computation of
S1 (1.3) at the other root: at any root s of the critical equation one has 1 - s = x s (1 - tau s), hence (1 - tau s)^2 phi''_tau(s) =
(2 tau s - c/x)/s^2, and 2 tau s_2 - c/x = delta/x. The rows r_{a,2}[m] are analytic on the whole disc |tau| < r0: by 10.1-10.3 they are
polynomials in the W_k(U_2), in 1/W_2(U_2) = 1/(x delta s*^2), in L_2 = x s* and in t-derivatives of these and of log Delta_2 = log(x^3 delta s*^2),
and U_2 is not in {1, 1/Q} on that disc (at tau = 0, U_2 = -1/x).
Let B be an open disc with closure in the punctured disc, and fix on B branches of log s, log(1 - s) and log(1 - tau s) near s = s_2(tau)
(possible because s_2, 1 - s_2 and 1 - tau s_2 do not vanish on the closure of B). Then phi_tau(s) is analytic near every point (tau, s_2(tau)),
tau in B. O_B denotes the analytic functions on B and O_B((z)) the Laurent series in z with finitely many negative powers and coefficients in
O_B; d/dt = (1 - tau) d/dtau acts coefficientwise, and identities in O_B((z)) are identities coefficientwise in z. Put
    f0(tau) := phi_tau(s_2(tau)),   mu := x + c (1 - tau)/tau,   omega := (c + z)/tau - (1 + z),   D_2 := d/dt - Delta_2'/(2 Delta_2).
Here ' = d/dt; omega is the omega of S2 Lemma S2.6, and mu is the mu of S2 (5.3). f0 depends on the branches only through an additive
constant, so f0' does not depend on them.

Lemma 10.5 (Row 2 of the QDE at the regular idempotent). In M[[z]],
    (L_2 + z D_2) r_{1,2} = -(x/tau) r_{0,2} + (c/tau) r_{1,2}.                                                      (10.11)
Equivalently z D_2 r_{1,2} = L_1 (r_{1,2} - L_2 r_{0,2}), i.e. (D_2 r_{1,2})[m-1] = L_1 (r_{1,2}[m] - L_2 r_{0,2}[m]) for m >= 1 (the form
tested in V5 for i = 2). Together with the H-row r_{1,2} = (L_2 + z D_2) r_{0,2} (Row 1, the definition, 10.2) this is
    z D_2 v_2 = (A - L_2) v_2,     A = [[0, 1], [-x/tau, c/tau]].                                                    (10.12)
The same holds for the point-frame rows h_w v_2. In the module language of S2.5, (10.12) says nabla S^(2) = A S^(2) for S^(2) := nu_2 v_2 eps_{L_2}
with nu_2^2 = 1/Delta_2, because nabla(nu_2 f eps_{L_2}) = nu_2 (L_2 + z D_2) f eps_{L_2} (as in S2 (5.2.1)). This is the outline's
"z ds/dt = A s" for s = e^{u_2/z} Delta_2^(-1/2) v_2.

Inputs. Each is proved by hand where it is cited.
 (a) [10.2] r^pt_{0,2}(z; t) = Sum_m z^m wickZ[[m+1]]|_{W_k -> W_k(U_2(t))}, r_{a,2} = h_w^-1 r^pt_{a,2}, and the H-row r_{1,2} = (L_2 + z D_2) r_{0,2},
     where d/dt is the total derivative along U_2(t).
 (b) [V1 at U_2: S2 Prop. S2.3(c); the tensor rule: S2 Prop. S2.4(a) at U_2] Delta_2 = x^2 W_2(U_2) and W_2(U_2) = U_2^2 f2/tau^2. Hence
        Delta_2 = x^2 U_2^2 f2/tau^2,                                                                                    (10.13)
     consistent with (10.1): with tau s_2 = 1/(x s*) one gets W_2(U_2) = delta/(x tau^2 s_2^2) = x delta s*^2.
 (c) [S2 Lemma S2.2(b), (d); the envelope formula] L_2 = (x U_2 + 1)/tau = 1/(tau s_2) = c/tau - x s_2 (the last equality is the critical
     equation divided by tau s_2). Since s_2 is critical, d f0/dtau = (d_tau phi_tau)(s_2) = s_2/(1 - tau s_2). Hence
        f0' = (1 - tau) s_2/(1 - tau s_2) = 1 - x s_2 = L_2 - mu.                                                        (10.14)
     The middle equality holds because (1 - x s_2)(1 - tau s_2) = 1 - (x + tau) s_2 + x tau s_2^2 = (c - x - tau) s_2 = (1 - tau) s_2 (the
     critical equation). The last one holds because c/tau - x s_2 - (1 - x s_2) = c/tau - 1 = mu.
 (d) [(I8)] S1 (1.1); P1(iii), P2 (2.1), P3, (D1)-(D3) and the amplitude identity S1 (4.1), at the nondegenerate critical point s_2(tau).
 (e) [(I9)] S2 Lemma S2.6 (V7).

Proof. Fix tau0 with 0 < |tau0| < r0 and a disc B around tau0 as above.
Step A (the regular row is a normalised Laplace expansion in the Euler coordinate). For every tau in B,
    r^pt_{0,2}(z; t) = U_2 L(tau),     L(tau) := Lnorm[phi_tau, (1 - tau s)^-1; s_2(tau)]   in O_B[[z]].                     (10.15)
This is S1 Prop. 3.1, steps (1)-(2), with (U_1, s*) replaced by (U_2, s_2).
(1) Let y_2 be a branch of log U_2 and Wy(y) := W(e^y). Then Wy is analytic near y_2, Wy^(k)(y_2) = W_k(U_2), Wy'(y_2) = W_1(U_2) = 0 and
Wy''(y_2) = W_2(U_2) = x delta s*^2 != 0. In the coordinate y the integrand e^{W/z} dU/U is e^{Wy/z} dy, with amplitude 1. By (D1) and the
description of wickZ (S1 s.0), r^pt_{0,2} = Lnorm[Wy, 1; y_2] = Lhat[Wy, 1; y_2].
(2) h(s) := log(1 - tau s), with the branch for which h(s_2) = y_2, is analytic near s_2, and h'(s_2) = -tau/U_2 != 0 because tau != 0. By S1
(1.1), Wy(h(s)) = phi_tau(s) + C(t) near s_2 with C(t) = c log tau + x t + x log(-1) independent of s, and the amplitude 1 becomes
(1 o h) h' = -tau (1 - tau s)^-1. By P1(iii), by (D3) (a constant added to the phase changes nothing) and by linearity (the constant factor -tau),
    Lhat[Wy, 1; y_2] = Lhat[phi_tau, -tau (1 - tau s)^-1; s_2] = Lhat[phi_tau, (1 - tau s)^-1; s_2] = U_2 L,
since the amplitude (1 - tau s)^-1 takes the value 1/U_2 at s_2. By (D2) the z-coefficients of L are polynomials in 1/f2, in the derivatives
phi_tau^(k)(s_2) and in the Taylor coefficients tau^a U_2^(-a-1) of the amplitude, so they are analytic on B (as in S1 Prop. 3.1(3)).
Step B (conjugation). Put
    Dt := d/dt + f0'/z - f2'/(2 f2),     Dop := Dt/(1 - tau) = d/dtau + (d f0/dtau)/z - (d f2/dtau)/(2 f2).
Dop is the operator of the proof of S1 Prop. 4.1, at s_2. Then, as operators on O_B((z)),
    (L_2 + z D_2) o U_2 = U_2 o (z Dt + omega),                                                                    (10.16)
where U_2 stands for multiplication by U_2. Indeed, for Y in O_B((z)),
    (L_2 + z D_2)(U_2 Y) = U_2 [ z Y' + L_2 Y + z (U_2'/U_2 - Delta_2'/(2 Delta_2)) Y ].
By (10.13), Delta_2'/(2 Delta_2) = U_2'/U_2 + f2'/(2 f2) - tau'/tau, with tau'/tau = (1 - tau)/tau; and L_2 = f0' + mu by (10.14). So the bracket
equals z Y' + f0' Y - z (f2'/(2 f2)) Y + (mu + z (1 - tau)/tau) Y = (z Dt + omega) Y, since mu + z (1 - tau)/tau = x + (c + z)(1 - tau)/tau = omega.
Only the logarithmic derivative of Delta_2 enters, so no square root and no branch is involved. (10.16) is where V1 at U_2 is used: if
Delta_2 were x^2 W_2(U_2) times a non-constant function p(t), the right side of (10.16) would contain the extra term -(z/2)(p'/p) U_2 Y, and
Step D would not give Row 2.
Step C (Gauss's equation for L). With alpha = beta = 1/z + 1 and gamma = c/z + 2,
    Gauss_Dop[L] := tau (1 - tau) Dop^2 L + (gamma - (alpha + beta + 1) tau) Dop L - alpha beta L = 0   in O_B((z)).    (10.17)
This is the first half of the proof of S1 Prop. 4.1, at s_2 instead of s*. P2 in the normalised form (2.1) applies to f(tau, s) = phi_tau(s),
to the analytic critical point s_2(tau) (f2 != 0), and to every amplitude q(tau, s) that is a Laurent polynomial in z with coefficients analytic
near (tau1, s_2(tau1)), tau1 in B. With delta_tau q := d_tau q + q (d_tau phi_tau)/z it says Dop Lnorm[phi_tau, q; s_2] = Lnorm[phi_tau, delta_tau q; s_2];
the critical value f0 and the Hessian f2 in Dop are those of s_2. Apply it to q = g := (1 - tau s)^-1 and then to q = delta_tau g, and use
(D3) (linearity over functions of tau and over C[z, 1/z]):
    Gauss_Dop[L] = Lnorm[ tau (1 - tau) delta_tau^2 g + (gamma - (alpha + beta + 1) tau) delta_tau g - alpha beta g ; s_2 ]
                 = -(1 + 1/z) Lnorm[ d_s k + (d_s phi_tau) k/z ; s_2 ] = 0,      k := s (1 - s)(1 - tau s)^-2.
The second equality is the amplitude identity S1 (4.1). It is an identity of rational functions of (s, tau), valid wherever both sides are
defined, in particular near s_2. The last equality is P3 (S1 Prop. 2.2). The amplitudes g, delta_tau g, delta_tau^2 g and k are analytic near
s_2 because s_2 is not in {0, 1, 1/tau}. That s_2 -> oo as tau -> 0 plays no role, since B does not contain 0.
Step D (transfer by V7). S2 Lemma S2.6 is an identity in the ring of differential operators generated by functions of (tau, z) and one operator
D, subject to the Leibniz rule D a = a D + D(a) with D(tau) = z (1 - tau); in it d/dtau stands for (z (1 - tau))^-1 D. Its proof uses nothing
else (S2 notes that it therefore holds with D replaced by nabla on every E_lambda). The operator z Dt = D + a0, with the function
a0 := f0' - z f2'/(2 f2), satisfies the same Leibniz rule, z Dt (a Y) = D(a) Y + a z Dt Y, and (z (1 - tau))^-1 z Dt = Dop. Hence S2 Lemma S2.6 gives
    (z Dt + omega)^2 - (c/tau)(z Dt + omega) + x/tau = ((1 - tau) z^2/tau) Gauss_Dop   on O_B((z)).                      (10.18)
(Directly: (z Dt + omega)^2 Y = (z Dt)^2 Y + 2 omega z Dt Y + D(omega) Y + omega^2 Y and (z Dt)^2 = z^2 (1 - tau)^2 Dop^2 - z^2 (1 - tau) Dop; the
comparison of coefficients in S2's proof is then literally the same, with f_tau and f_tautau replaced by Dop Y and Dop^2 Y.) By (10.15), (10.16)
applied twice, (10.18) and (10.17),
    [(L_2 + z D_2)^2 - (c/tau)(L_2 + z D_2) + x/tau] r^pt_{0,2} = U_2 [(z Dt + omega)^2 - (c/tau)(z Dt + omega) + x/tau] L
                                                                = U_2 ((1 - tau) z^2/tau) Gauss_Dop[L] = 0.
h_w does not depend on t, so it commutes with L_2 + z D_2, and the same identity holds for r_{0,2} = h_w^-1 r^pt_{0,2}. With
r_{1,2} = (L_2 + z D_2) r_{0,2} it is (10.11) on B.
Conclusion. Every z-coefficient of (10.11) is an identity of analytic functions on B. Since tau0 was arbitrary, it holds on the punctured disc
0 < |tau| < r0. All its terms are meromorphic at t = 0 (the r_{a,2}[m], L_2 and Delta_2'/Delta_2 are analytic there, and x/tau, c/tau have simple
poles), so it holds in M, i.e. for the Laurent expansions at t = 0. The second form uses c/tau - L_2 = L_1 and x/tau = L_1 L_2 (S2 Lemma
S2.2(d)). The first component of (A - L_2) v_2 is r_{1,2} - L_2 r_{0,2} = z D_2 r_{0,2} (Row 1), and the second is
-(x/tau) r_{0,2} + (c/tau - L_2) r_{1,2} = z D_2 r_{1,2} (Row 2); this is (10.12). QED

Remark 10.5.1 (the formal form; S1 Prop. 2.4(ii) and S2.10 (A2), (B2)). Formally, with u_2 := W(U_2) = f0 + C(t), Step A and (10.13) say
    e^{u_2/z} Delta_2^(-1/2) r^pt_{0,2} = K(z) Q^{x/z} tau^{c/z + 1} M_2,      M_2 := e^{f0/z} (-f2)^(-1/2) L,
with the t-independent constant K(z) = (-x^2)^(-1/2) e^{x log(-1)/z}. It contains the phase kappa(z) = -e^{x log(-1)/z} of S1 (1.2) and a root
of -1/x^2. M_2 is the full (unnormalised) formal expansion of the Euler integrand at s_2. In this form the statement is S1 Prop. 2.4(ii) (the
Lfull form, with the square roots chosen continuously on B), with the t-dependent prefactor e^{C/z} = -kappa(z) Q^{x/z} tau^{c/z} of S1 (1.1)-(1.2)
and the Jacobian -tau (1 - tau s)^-1, combined with V1 at U_2, which replaces (-W_2(U_2))^(-1/2) by a constant times Delta_2^(-1/2). Step C says
Gauss[M_2] = 0, and Step D is V7 applied to P M_2 with P = (1 - tau)^(-x/z) tau^(c/z+1) (outline, Step 5). The exponential symbols are not
elements of a ring, and (10.16) is their precise content: their logarithmic derivatives are the terms f0'/z and mu/z contained in z Dt + omega,
while the Jacobian and the Hessian normalisation contribute U_2 and z (1 - tau)/tau - z f2'/(2 f2). S2.10 (A2), (B2) say the same in the module
language E_lambda of S2.5: (A2) is (10.15), and (B2) is (10.17) transported to E_{lambda_2} with lambda_2 = f0' = 1 - x s_2.
Why the normalised forms suffice. The code's row is itself normalised (leading coefficient 1, no Gaussian prefactor), so P1(iii) identifies it
with U_2 L exactly, with no convention. The t-dependent factors that the Lfull form of P1 carries (e^{C/z}, the Jacobian, and (-f2)^(-1/2) in
place of Delta_2^(-1/2)) are then accounted for in (10.16): through mu, through U_2 and tau'/tau, and through (10.13). Likewise P2 is used in
the normalised form (2.1), which S1 Prop. 2.3 states as equivalent to its unnormalised form; S1 Prop. 4.1 uses it in the same way. The two
routes are therefore equivalent, and the normalised one avoids choosing square roots.

Corollary 10.5' (orthogonality (I5); this is S2 Prop. S2.10). <v_1(w), v_2(-w)> = 0 in M[[w]].
Proof. Y(z) := <v_1(z), v_2(-z)> lies in M[[z]], since v_1 and v_2 lie in M[[z]]^2. The QDE at U_1 (S2 Prop. S2.7: z D_1 v_1 = (A - L_1) v_1) and
(10.12) with z replaced by -z (A, L_2 and D_2 do not involve z: -z D_2 v_2(-z) = (A - L_2) v_2(-z)) are used together with
v_i' = D_i v_i + (Delta_i'/(2 Delta_i)) v_i, the t-independence of eta^-1 and the self-adjointness <A a, b> = <a, A b> (S2 Prop. S2.8(c)):
    z Y' = <(A - L_1) v_1(z), v_2(-z)> - <v_1(z), (A - L_2) v_2(-z)> + (z/2)(Delta_1'/Delta_1 + Delta_2'/Delta_2) Y
         = [ (L_2 - L_1) + (z/2)(Delta_1'/Delta_1 + Delta_2'/Delta_2) ] Y.
Write Y = Sum_m Y_m z^m with Y_m in M. The z^0 coefficient gives (L_2 - L_1) Y_0 = 0, and for m >= 1 the z^m coefficient gives
(L_2 - L_1) Y_m = Y_{m-1}' - (1/2)(Delta_1'/Delta_1 + Delta_2'/Delta_2) Y_{m-1}. Since L_2 - L_1 = -delta/tau is a unit of M, Y_0 = 0, and
inductively Y_m = 0 for all m. QED
That Y is a power series in z, bounded below, is essential; the formal symbol e^{-(u_1 - u_2)/z}, an infinite series in 1/z, would be a non-zero
solution (S2.10). Y_0 = <chi_1, chi_2> = 0 is also direct (S2 Prop. S2.8(d)). Confirmations, not used: V6 (exact in Q, w^n for n <= 9), V5 for
i = 2 (m <= 9, t^0..t^6), and S4_checks_r1.pl B3 (Row 2 at U_2 from the code's definitions, z^1..z^6, t^0..t^6, four values of x).

10.6 The regular rows at t = 0.
Proposition 10.6.
 (a) r_{0,2}(z; 0) = h_w(z)^-2; equivalently, in the point frame, r^pt_{0,2}(z; 0) = h_w(z)^-1.
 (b) (Uses Lemma 10.5.) r_{1,2}(z; 0) = (x/c) h_w(z)^-2, so v_2(z; 0) = h_w(z)^-2 (1, x/c). Equivalently, [t^1] r_{0,2}(z; t) = 0, i.e.
     d/dt [Delta_2^(-1/2) r_{0,2}(z; t)] = 0 at t = 0.
 (c) (Uses (b).) E_22(z, w)|_{t=0} = (x^3/c) [1 - h_w(z)^2 h_w(w)^2]/(z + w), and (by (a) alone) Dil_2(n)|_{t=0} = (-1)^n [z^(n-1)] h_w(z)^-2.
Proof. (a) By 10.2, each r^pt_{0,2}[m] is analytic at t = 0, and its value there is wickZ[[m+1]] evaluated on the t = 0 data W_k(U_2(0)) at Q = 1,
i.e. on (theta^k W_0)(-1/x) with W_0(U) := c log(1 - U) - log U (spotU at Q = 1). By (D1), r^pt_{0,2}(z; 0) = Lnorm[Wy_0, 1; y_c] with
Wy_0(y) := W_0(e^y) and y_c a branch of log(-1/x). We identify this series with S1's analytic lemma, first for real x > 0.
 (i) The phase. Write y_c = y0 + i pi with y0 := -log x, and put F(y) := y - c log(1 + e^y) for real y. Since e^{y + i pi} = -e^y,
     Wy_0(y + i pi) = c log(1 + e^y) - y - i pi = -F(y) - i pi, so the Taylor coefficients of Wy_0 at y_c of order k >= 1 are those of -F at y0.
     In (D1) the replacement (f, z) -> (-f, -z) leaves (-z/f2)^n and R_f/z unchanged, so Lnorm[-f, g; u0](z) = Lnorm[f, g; u0](-z).
     With (D3),
        r^pt_{0,2}(z; 0) = Lnorm[-F, 1; y0](z) = Lnorm[F, 1; y0](-z).                                                (10.19)
     Equivalently: at z = -1/lambda and U = -e^y one has W_0/z = lambda F(y) + const, so the Wick series in z is the Laplace series of
     Int e^{lambda F(y)} dy in the variable 1/lambda = -z.
 (ii) Laplace. F is real-analytic on I = R. F'(y) = (1 - x e^y)/(1 + e^y) vanishes only at y0, with F' > 0 on (-oo, y0) and F' < 0 on (y0, oo);
     F''(y0) = -c e^{y0}/(1 + e^{y0})^2 = -x/c < 0 and F(y0) = x log x - c log c. So F is strictly increasing on (-oo, y0] and strictly
     decreasing on [y0, oo), and sup_{|y - y0| >= e} F = max(F(y0 - e), F(y0 + e)) < F(y0) for every e > 0. Moreover
     Int_R e^{F(y)} dy = Int_0^oo (1 + u)^(-c) du = 1/x < oo (u = e^y). S1 Lemma 2.5 (g = 1, z0 = 1) gives, for every M, as z' -> 0+,
        Int_R e^{F(y)/z'} dy = e^{F(y0)/z'} (2 pi z' c/x)^(1/2) [ Sum_{m<M} ([z'^m] Lnorm[F, 1; y0]) z'^m + O(z'^M) ].
     With u = e^y and lambda = 1/z' the left side is Int_0^oo u^(lambda - 1) (1 + u)^(-c lambda) du = B(lambda, x lambda) (DLMF 5.12.3;
     c lambda - lambda = x lambda > 0).
 (iii) Stirling. Apply Stirling's series with remainder (S1 Prop. 6.1(3); DLMF 5.11.1 with the error bound 5.11.10) to
     log B(lambda, x lambda) = log Gamma(lambda) + log Gamma(x lambda) - log Gamma(c lambda); all three arguments tend to +oo. The terms -u
     cancel, the terms in log lambda add up to -(1/2) log lambda, and for every K
        B(lambda, x lambda) = (2 pi/lambda)^(1/2) x^(x lambda - 1/2) c^(1/2 - c lambda) exp( Sum_{k=1}^{K} beta_k S_k z'^(2k-1) + O(z'^(2K+1)) ),
        beta_k := B_2k/(2k(2k-1)),   S_k := 1 + x^(1-2k) - c^(1-2k),   z' = 1/lambda.
     The prefactors agree: e^{F(y0)/z'} (2 pi z' c/x)^(1/2) = x^(x lambda) c^(-c lambda) (2 pi c/(lambda x))^(1/2). By the uniqueness of
     asymptotic power series (S1 Prop. 6.1(4)), Lnorm[F, 1; y0](z') = exp(Sum_k beta_k S_k z'^(2k-1)). By (10.19),
        r^pt_{0,2}(z; 0) = exp(-Sum_k beta_k S_k z^(2k-1)) = h_w(z)^-1,
     because log h_w(z) = Sum_k beta_k S_k z^(2k-1) (outline s.0; Sum_a w_a^(1-2k) = c^(1-2k) - x^(1-2k) - 1 for w = (c, -x, -1)).
 (iv) All x. Both sides have z-coefficients in Q(x), as in S1 Prop. 6.1(1). On the left this follows from (D2), which involves only F^(k)(y0)
     for k >= 2: F' = 1 - c ell with ell(y) := e^y/(1 + e^y), ell' = ell (1 - ell), so for k >= 2 F^(k) = -c ell^(k-1) is a polynomial in ell
     with coefficients in Z[c], and ell(y0) = 1/c; also 1/F''(y0) = -c/x. On the right it holds by definition. Two rational functions of x
     that agree on (0, oo) are equal.
 Frame B multiplies by a further h_w^-1. Check at z^1, by hand: at U = -1/x, Q = 1 one has W_2 = x/c, W_3 = -x (1-x)/c^2, W_4 = x (1 - 4x + x^2)/c^3,
 and (D1) gives W_4/(8 W_2^2) - 5 W_3^2/(24 W_2^3) = -(1 + x + x^2)/(12 x c) = -(1/12) S_1 = [z^1] h_w^-1. At idempotent 1 the same method
 (S1 Prop. 6.1) gives the opposite sign of the tail: A(z) = h_w(z) c/(c+z).
(b) First proof, from Row 2 at U_2 alone. The second form of (10.11) reads (D_2 r_{1,2})[m-1] = L_1 (r_{1,2}[m] - L_2 r_{0,2}[m]). The left side is
analytic at t = 0 by 10.2, while L_1 = c/t + O(1) has a simple pole. Hence r_{1,2}[m](0) = L_2(0) r_{0,2}[m](0) = (x/c) r_{0,2}[m](0) for every m.
Second proof, from (I5) (Corollary 10.5'). Put w = -z in (I5) and multiply by rho: <sigma(z), v_2(z)> = 0. Its t^-1 coefficient is
<(0, c), v_{2,0}(z)> = c x (x r_{0,2}(z; 0) - c r_{1,2}(z; 0)) = 0 (N4).
The equivalent forms. By the H-row, r_{1,2}[m](0) - (x/c) r_{0,2}[m](0) = (D_2 r_{0,2}[m-1])(0) = (d/dt r_{0,2}[m-1])(0), because L_2(0) = x/c and
Delta_2'(0) = 0 (see 10.3); and D_2 f = Delta_2^(1/2) d/dt (Delta_2^(-1/2) f).
(c) Insert v_{2,0} into (1.2). By (N4), P_22(z, w)|_{t=0} = -(x^3/c) h_w(-z)^-2 h_w(-w)^-2 = -(x^3/c) h_w(z)^2 h_w(w)^2, using h_w(z) h_w(-z) = 1
(the exponent of h_w is odd in z). Also P_22(-w, w)|_{t=0} = -x^3/c. The formula for Dil_2 follows from (a). QED

Remark 10.6.1 (E_22 at t = 0, and conifold_regularity_note.md F3).
 (1) E_22 at t = 0 is not of rank one, unlike E_11 and E_12 (Step 8, CF.6(b)). Write h_w(z)^2 = exp(Sum_k s_k z^(2k-1)) with s_k := 2 beta_k S_k
     (beta_k, S_k as in (a)(iii)); so s_1 = S_1/6, s_2 = -S_2/180, and S_1 = (1 + x + x^2)/(x c). Expanding (c) to second order,
        (1 - h_w(z)^2 h_w(w)^2)/(z + w) = -s_1 - (s_1^2/2)(z + w) - s_2 (z^2 - z w + w^2) - (s_1^3/6)(z + w)^2 + (terms of degree >= 3),
     so at t = 0: E_22(0, 0) = -(x^3/c) s_1 = -x^2 (1 + x + x^2)/(6 c^2), E_22(1, 0) = E_22(0, 1) = -(x^3/(2c)) s_1^2 and
     E_22(1, 1) = (x^3/c)(s_2 - s_1^3/3). Hence
        [ E_22(0,0) E_22(1,1) - E_22(0,1) E_22(1,0) ]_{t=0} = (x^3/c)^2 (s_1^4/12 - s_1 s_2) = (x^6/c^2) S_1 (72 S_2 + 5 S_1^3)/77760,
     a non-zero rational function of x (it equals 3/4096 at x = 1, where E_22(0,0) = -1/8 happens to coincide with E_11(0,0)).
 (2) F3 of conifold_regularity_note.md states "E_11(z, w) = E_22(z, w) = -(x^3/2) phi(z) phi(w)" at t = 0. Its E_22 part is the E_11 table,
     read through the conjugate encoding of 10.2. The cause is in hvc_vertex_identity.wls, not in the mathematics:
     (i) cconj[C[a, b]] = C[a + b Ssum, -b] is U -> Ssum - U on QA = Q(x, Q)[U]/(LeqU). The evaluations ev_i: C[a, b] -> a + b U_i(Q) are ring
         homomorphisms, jointly injective (LeqU has distinct roots, delta != 0), and ev_1 o cconj = ev_2, ev_2 o cconj = ev_1 (U_1 + U_2 = Ssum).
         So cconj is a ring automorphism (the Galois conjugation), commuting with cadd, cmul, cinv and with cscale by elements of Q(x, Q).
     (ii) cconj commutes with cDerN. Indeed cDerN[C[a, b]] = C[d_Q a, d_Q b] + b dUC with ev_i(dUC) = dU_i/dQ (dUC = toC[-D[LeqU, Q]/D[LeqU, U]],
         implicit differentiation), so ev_i o cDerN = d/dQ o ev_i (the total derivative along the root). For {i, i'} = {1, 2} this gives
         ev_i(cDerN cconj C) = d/dQ ev_{i'}(C) = ev_{i'}(cDerN C) = ev_i(cconj cDerN C). twistRows multiplies the rows of both idempotents by the
         same scalar series, so it also commutes with cconj.
     (iii) At alpha0 = 0, rmatrixAt0 uses ffc1 = ffc2 = (1, 0, 0, ...), so ratCN[0,2][m] = cconj[ratCN[0,1][m]]; LCN[2] = cconj[LCN[1]] and
         DeltaCN[2] = cconj[DeltaCN[1]], hence dlogDeltaN[2] = cconj[dlogDeltaN[1]] by (i), (ii). Applying (i), (ii) to the H-row formula of
         rmatrixCore (whose only scalar is Q) gives ratCN[1,2][m] = cconj[ratCN[1,1][m]] for every m, and this persists after twistRows.
         etainvCur contains no U, so evalQ gives exactly
         EdN[2,2][m, n] = cconj[EdN[1,1][m, n]], EdN[2,1][m, n] = cconj[EdN[1,2][m, n]] and DilN[2][n] = cconj[DilN[1][n]].
     (iv) The true value of every factor is its value at U_1 (see 10.2); hvc_vertex_identity.wls evaluates E11 and E12 in this way (val0[., U1]).
         Its line 14 prints val0[EdN[2,2][m, n], U2] = ev_2(cconj EdN[1,1][m, n]) = ev_1(EdN[1,1][m, n]), i.e. the E_11 table -(x^3/2)/c^(m+n+2).
         Evaluated at U_1 it gives (c). The graph sums are not affected: in graphValueC, cconj maps the term of a colouring c to the term of the
         colouring 3 - c, so the sum over all colourings is cconj-invariant, i.e. of the form C[a, 0] (orbitTerm records this as cSymQ), and its
         value a is the same at U_1 and at U_2.
     So F3 of the note should read E_22|_{t=0} = (x^3/c) [1 - h_w(z)^2 h_w(w)^2]/(z + w) (proposed correction in S4.8). The referee's independent
     recomputation from the code's definitions agrees with (c) at 21 of 21 entries for each of six values of x, and with the note's value at 0 of
     21 entries (1 of 21 at x = 1); S4_checks_r1.pl B5 confirms this at x = 3/7, 5/3, -5/2, 2/9 (and at x = 1: 1 of 21, minor 3/4096). Requested
     check R5 (S4.7) is therefore a confirmation, not an open question. Nothing in Lemma C or in Theorem R depends on E_22 beyond its regularity,
     and Theorems A and B and the Reduction Lemma of the note use only E_11 and E_12 at t = 0 and the regularity of the other factors.
Remark 10.6' (the regular rows are determined by the singular ones; not needed below). Combine three facts: (I5) at order t^k, k >= -1
(the sum is empty for k = -1),
    <(0, c), v_{2,k+1}(z)> = -Sum_{k1=0}^{k} <sigma_{k1}(z), v_{2,k-k1}(z)>;
the H-row in the form x r_{0,2} - c r_{1,2} = (x - c L_2) r_{0,2} - c z (d/dt - Delta_2'/(2 Delta_2)) r_{0,2}; and x - c L_2 = -x^2 tau s*^2 = O(t).
Together they give
    c (k+2) z [t^(k+2)] r_{0,2} = [t^(k+1)]((x - c L_2) r_{0,2}) + c z [t^(k+1)](r_{0,2} Delta_2'/(2 Delta_2)) + (1/(c x)) Sum_{k1=0}^{k} <sigma_{k1}, v_{2,k-k1}>,
whose right side involves r_{0,2} only up to order t^(k+1). Hence, given r_{0,2}(z; 0) = h_w^-2, the whole regular column (and with it every
E_12 leg and E_22 to all orders) is determined recursively by the closed-form singular data (S4.4) and by L_2 and Delta_2. For example
[t^1] r_{0,2} = 0.
Remark 10.6'' (K_2 = 1; not needed below). K_2(z) := <v_2(z), v_2(-z)>/Delta_2 lies in M[[z]], and its coefficients are analytic at t = 0. It is
t-independent: by (10.12) and its z -> -z form, z <D_2 v_2(z), v_2(-z)> + z <v_2(z), D_2 v_2(-z)> = <(A - L_2) v_2(z), v_2(-z)> -
<v_2(z), (A - L_2) v_2(-z)> = 0 (S2 Prop. S2.8(c)), so with v_2' = D_2 v_2 + (Delta_2'/(2 Delta_2)) v_2 one gets
d/dt <v_2(z), v_2(-z)> = (Delta_2'/Delta_2) <v_2(z), v_2(-z)>, i.e. K_2' = 0 (as in the proof of S2 Prop. S2.9(a)). Hence K_2 = K_2|_{t=0}, and by
10.6(b), (N4) and Delta_2(0) = x^3/c,
    K_2(z) = <h_w(z)^-2 (1, x/c), h_w(-z)^-2 (1, x/c)>/(x^3/c) = h_w(z)^-2 h_w(-z)^-2 = 1.
With K_1 = 1 (S2 Prop. S2.9(c)) and Corollary 10.5' this proves V6 for all i, j: <v_i(w), v_j(-w)> = delta_ij Delta_i, so (S-a) of S4.1 holds
for every pair of idempotents. (S4_checks_r1.pl B7 confirms K_2 = 1 from the code's rows, z^0..z^6, t^0..t^6.)

10.7 F1 and F3, complete.
F1. Delta_1 = -x delta/(tau^2 s*^2) = -x c^3/(t^2 rho^2) = -(x c^3/t^2) (1 + (1-x)^2 t/c^2 + O(t^2)), by (I3), (I7) and CF.2(d). Delta_2 = x^3 delta s*^2
is analytic with Delta_2(0) = x^3/c != 0 (see 10.3).
F3. For F_{g,0}, graphValueC uses the factors deltaPow (vertex weights), DilN and EdN; LfN occurs only with marked points. All frame-B
dilaton and edge factors are analytic at t = 0:
    factor   proof of regularity                                        value at t = 0
    Dil_1    q_t = rho^-1 (rho q_t), rho q_t = q_0 + O(t) (L1, (I2))    q_0(psi) = -c psi/(c - psi) = -psi + T(psi), T = -psi^2/(c - psi)
    Dil_2    10.4                                                       (-1)^n [z^(n-1)] h_w(z)^-2 (Prop. 10.6)
    E_11     (I6): t^2 rho^2 E_11 = O(t^2); rho^-2 in C[[t]]            -(x^3/2) lam_1(psi) lam_1(psi')
    E_12     Prop. 9.3(i)                                               lam_1(psi) chi(w), chi(w) = -x^3 h_w(w)^2/(c + w) (CF.6(b))
    E_21     10.6(b) (Row 2 at U_2, Lemma 10.5); Prop. 9.5 ((I5),       E_21(b, a) = E_12(a, b)
             Corollary 10.5')
    E_22     10.4                                                       (x^3/c) [1 - h_w(z)^2 h_w(w)^2]/(z + w) (10.6(c); not rank one, Remark 10.6.1)
Remark. By (9.6), the t^-1 part of sigma is the psi-independent covector (0, c). Hence the leg factors LfN[1][n] = (-1)^n r_{1,1}[n] are regular
for n >= 1, while LfN[1][0] = L_1 has a simple pole. They do not enter F_{g,0}.


## S4.4 Closed forms (section 4 of the outline)

Notation. The function zeta(t) := tau e^(dphi(tau)/c) of CF.1 was called kappa in the first version of this section. It is renamed because
kappa(z) is the phase constant of S1 (1.2) and of the outline, Step 5 (and kappa_0, kappa_1 are constants of S2.5). The ratio zeta/t is called
omega in S3 s.1.6 and in S4_checks.pl; that omega is unrelated to the omega(z, tau) of S2 Lemma S2.6 used in 10.5.

Theorem CF.1 (the dilaton leaf coefficients). For j >= 1 put
    C_j := prod_{m=1}^{j-1} (j - m c)^2 / ( j! (j-1)! c^(2j-4) ) = c^2 [(1 - j/c)_{j-1}]^2 / ( j! (j-1)! )        (C_1 = c^2),
    zeta(t) := tau e^(dphi(tau)/c) = c (1 - U_1) U_1^(-(1-x)/c) = t + (1/c^2 - 1/2) t^2 + O(t^3)        (second form: CF.2(c)),
    A_j(t) := (-1)^j C_j t^-1 zeta(t)^j 2F1(j x/c, j x/c; j; tau),       so A_j = O(t^(j-1)).
Then a_j^(k) = [t^k] A_j(t) for all j >= 1 and k >= 0. Equivalently, as formal series in t with finitely many j at each order,
    z G(z; t) = c + Sum_{j>=1} A_j(t)/(c + j z),     rho q_t(psi) = c + Sum_{j>=1} A_j(t) lam_j(psi),     g(psi) = -(1/c) Sum_{j>=1} j A_j(t) lam_j(psi).
Proof. Fix j >= 1 and k >= 0. By (I2), [t^k](z G) is a rational function of z whose only possible singularity at z = -c/j is a simple pole,
contributing a_j^(k)/(c + j z). Hence a_j^(k) = j Res_{z=-c/j} [t^k](z G). Write
    z G = (tau/t) Sum_{n>=0} tau^n (c z/(c+z)) e^(-dphi/z) F_n(z),      F_n(z) := prod_{l=1}^{n} (1 + l z)^2 / ( n! z^n prod_{l=2}^{n+1} (c + l z) ).
In Pochhammer form F_n = (alpha)_n^2/((gamma)_n n!) with alpha = 1/z + 1 and gamma = c/z + 2. Since tau = O(t) and dphi = O(t), [t^k](z G) is a
finite sum over n <= k of terms (scalar) (c z/(c+z)) z^(-p) F_n(z): each [t^m] e^(-dphi/z) is a polynomial in 1/z, and evaluating it at z = -c/j
commutes with taking t-coefficients. We take residues term by term.
 j = 1. F_n is regular at z = -c, since its poles are at 0 and at -c/l with l >= 2. Hence
    Res_{z=-c} [t^k](z G) = [t^k] ( (c z)|_{z=-c} (tau/t) e^(dphi/c) F(-c, tau) ) = [t^k] ( -c^2 (tau/t) e^(dphi/c) 2F1(x/c, x/c; 1; tau) ),
 because at z = -c, alpha = 1 - 1/c = x/c and gamma = 1.
 j >= 2. Here c z/(c+z) and z^(-p) are regular at -c/j. For n < j - 1, F_n is regular there. For n >= j - 1 it has a simple pole, coming from the
 single factor c + j z = j (z + c/j). At z = -c/j one has 1 + l z = (j - l c)/j, c + l z = c (j - l)/j and z^n = (-c/j)^n, and
 prod_{l=2, l != j}^{n+1} (j - l) = (j-2)! (-1)^(n+1-j) (n+1-j)!. Hence
    Res_{z=-c/j} F_n = (-1)^(j-1) prod_{l=1}^{n} (j - l c)^2 / ( j^2 c^(2n-1) n! (j-2)! (n+1-j)! ) = (-1)^(j-1) c [(1 - j/c)_n]^2 / ( j^2 n! (j-2)! (n+1-j)! ),
 using prod_{l=1}^{n} (j - l c) = (-c)^n (1 - j/c)_n. Reindex n = j - 1 + m, and use (1 - j/c)_{j-1+m} = (1 - j/c)_{j-1} (j x/c)_m (because
 1 - j/c + j - 1 = j x/c) and (j-1+m)! = (j-1)! (j)_m:
    Sum_{n>=j-1} tau^n Res_{z=-c/j} F_n = (-1)^(j-1) c [(1 - j/c)_{j-1}]^2 tau^(j-1) / ( j^2 (j-2)! (j-1)! ) 2F1(j x/c, j x/c; j; tau).
 Only the terms with n >= j - 1 contribute, and their sum is this single 2F1; it is the limit formula DLMF 15.2.3_5 for the regularised
 function 2F1(a, b; gamma; tau)/Gamma(gamma) at gamma = 2 - j. The remaining factors at z = -c/j are c z/(c+z) = -c/(j-1) and e^(-dphi/z) = e^(j dphi/c).
 Multiplying by j and using j (j-1) (j-2)! = j!:
    a_j^(k) = [t^k] ( (-1)^j c^2 [(1 - j/c)_{j-1}]^2 / (j! (j-1)!) (tau^j/t) e^(j dphi/c) 2F1(j x/c, j x/c; j; tau) ) = [t^k] A_j(t).
 The case j = 1 is the same formula (C_1 = c^2, empty Pochhammer symbol).
The three generating-function identities follow from (I2) and Lemma 9.1. QED
Analytic remark (not needed). For fixed small t, z -> z G(z; t) is meromorphic near -c/j, and A_j(t) = j Res_{z=-c/j} z G(z; t) as functions. The
2F1 series converges uniformly on a small circle around -c/j, and its poles (gamma(z) in {0, -1, -2, ...}) accumulate only at z = 0.

CF.2 Corollaries.
 (a) Sum rule: Sum_{j>=1} A_j(t) = -c^2 as formal series in t, i.e. Sum_{j<=k+1} a_j^(k) = -c^2 delta_{k0}. This is the absence of a pole of G at
     z = 0 (Step 2(d)). Equivalently Sum_j (-1)^j C_j zeta^j 2F1(j x/c, j x/c; j; tau) = -c^2 t = c^2 log(1 - tau). Check in the limit x -> 0: there
     c = 1, C_j = 1/j, 2F1 = 1, s* = 1 and zeta = tau/(1 - tau), and the identity becomes Sum_j (-1)^j zeta^j/j = -log(1 + zeta) = log(1 - tau).
 (b) Everything is a function of tau. With B_j := t A_j = (-1)^j C_j zeta^j 2F1(j x/c, j x/c; j; tau),
        t z G(z) = -c log(1 - tau) + Sum_j B_j(tau)/(c + j z),    since t = -log(1 - tau).
 (c) e^(dphi) = (c s*)^c (1 - tau s*)^(x-1) = (c s*)^c U_1^(x-1). Proof: the critical equation gives 1 - s* = x s* (1 - tau s*), so
     dphi = log s* + x log(1-s*) - log(1 - tau s*) + log c - x log(x/c) = c log(c s*) + (x-1) log(1 - tau s*). All logarithms are of quantities
     tending to 1 and are taken on the principal branch, consistent with dphi(0) = 0. Hence zeta = c tau s* U_1^(-(1-x)/c) = c (1 - U_1) U_1^(-(1-x)/c).
     Also d(dphi)/dt = 1 - x s* (V3), equivalently d(dphi)/dtau = s*/(1 - tau s*).
 (d) With u = 4 x tau/c^2 and principal roots of 1 - u (which is near 1),
        rho = (tau/t) 2/( (1 - u)^(1/4) (1 + (1 - u)^(1/2)) ),     alg = tau/(t rho) = (1 - u)^(1/4) (1 + (1 - u)^(1/2))/2.
     This form has no branch ambiguity, also for c < 0. It follows from rho^2 = -x c^3/(t^2 Delta_1) and Delta_1 = -x delta/(tau^2 s*^2), which
     give rho^2 = (tau/t)^2 s*^2 c^3/delta, together with delta = c (1 - u)^(1/2), s* = 2/(c (1 + (1 - u)^(1/2))) and rho(0) = 1 (S2 Prop. S2.4(c)).
     (Written as (tau/t) s* (c^3/delta)^(1/2), the root must be c (1 - u)^(-1/4), whose value at t = 0 is c. For c < 0 the principal root of
     c^3/delta has the value |c| = -c there and would give rho(0) = -1; a form such as c^(3/2)/delta^(1/2) depends on how c^(3/2) is read.)
     Expansion: 2/((1 - u)^(1/4)(1 + (1 - u)^(1/2))) = 1 + u/2 + (11/32) u^2 + O(u^3), and
     rho = 1 - (1-x)^2 t/(2c^2) + (1 - 8x + 15x^2 - 8x^3 + x^4) t^2/(6c^4) + O(t^3) (CF.5).
 (e) Theta = -x^3 tau s*^2 = -4 x^3 tau/(c + delta)^2 = x c (1 - x s*) - x = -(x^3/c^2) t + (1-x)^2 x^3 t^2/(2 c^4) + O(t^3). The middle identity uses
     c s* - 1 = x tau s*^2.
 (f) The H-component. By (I4), g = -(1/c) Sum_j j A_j lam_j, psi lam_j = (c lam_j - 1)/j (N3) and Sum_j (t A_j)' = -c^2 (from (a)):
        sigma(psi) = (0, c/t) + Sum_{j>=1} S_j(t) lam_j(psi),
        S_j(t) = ( -(j/c) A_j,  A_j' - ((j-1)/t) A_j - (j/c) l A_j ) = t^(j-1) ( -(j/c) A~_j,  A~_j' - (j/c) l A~_j ),
     with A~_j := t^(1-j) A_j = (-1)^j C_j (zeta/t)^j 2F1(j x/c, j x/c; j; tau) analytic and A~_j(0) = (-1)^j C_j. The coefficient of t^k is (9.1).
     In particular S_1(0) = (c, x) = c (1, x/c).

CF.3 Top and subleading coefficients (all k). For every k >= 0 and every k >= 1 respectively,
    a_{k+1}^(k) = (-1)^(k+1) C_{k+1} = (-1)^(k+1) prod_{m=1}^{k} (k+1 - m c)^2 / ( (k+1)! k! c^(2k-2) ),
    a_k^(k)     = (-1)^k k C_k (1-x)^2/(2 c^2).
The first line is the formula of leaf_structure.md Step 3, previously verified for k <= 5 and now proved for all k. Note
prod_m (k+1 - m c)^2 = prod_m (k+1 - m - m x)^2, which gives (1-x)^2, (2-x)^2(1-2x)^2, 4(3-x)^2(1-x)^2(1-3x)^2, ...
Proof. A_j = (-1)^j C_j t^(j-1) (zeta/t)^j 2F1(j x/c, j x/c; j; tau), with (zeta/t)(0) = 1 and 2F1 = 1 + (j x^2/c^2) tau + O(tau^2). The first
claim is the case j = k+1 at t^k. For the second, [t^1] of (zeta/t)^k 2F1(k x/c, k x/c; k; tau) equals
k (1/c^2 - 1/2) + k x^2/c^2 = k (2 + 2x^2 - c^2)/(2c^2) = k (1-x)^2/(2c^2). QED
Sharpness. C_{k+1} != 0 unless x lies in the finite set {(k+1)/m - 1 : 1 <= m <= k}. So for generic x, lam_{k+1} really occurs at order t^k, and
the bound j <= k+1 in L1 is sharp.

Theorem CF.4 (closed form of the loops). With B_j = t A_j,
    rho^2 E_11(psi, psi') = Sum_{i,l>=1} e_il(t) lam_i(psi) lam_l(psi'),
    e_il(t) = (i l/(c^2 t^2)) Int_0^t Theta(u) B_i'(u) B_l'(u) du = -(i l x^3/(c^2 t^2)) Int_0^tau v (1-v) s*(v)^2 (dB_i/dv) (dB_l/dv) dv,
and e_il = O(t^((i-1)+(l-1))).
Proof. By CF.1, d/dt (t g(psi)) = -(1/c) Sum_j j B_j' lam_j(psi). Insert this in (I6) and integrate from 0; the integration constant is 0 because
t^2 rho^2 E_11 = O(t^2). For the second form substitute v = tau, so dt = dv/(1-v), B' = (1-v) dB/dv and Theta = -x^3 v s*^2. The order follows
from Theta = O(t) and B_j' = O(t^(j-1)). QED
Corollary CF.4' (corner, anti-diagonal; all k). Since Theta = -(x^3/c^2) t + O(t^2) and B_j' = j a_j^(j-1) t^(j-1) + O(t^j), for i + l = k + 2:
    e_il^(k) = -x^3 i^2 l^2 a_i^(i-1) a_l^(l-1) / ( c^4 (i + l) ) = -(-1)^(i+l) x^3 i^2 l^2 C_i C_l / ( c^4 (i + l) ).
With a_1^(0) = -c^2 this gives the corner
    e_{1,k+1}^(k) = ((k+1)^2/(k+2)) (x^3/c^2) a_{k+1}^(k),
which is the formula of leaf_structure.md Step 3, now proved for all k. The same follows directly from the Step 8 sum: only
(k0, k1, k2) = (1, 0, k) reaches lam_1(psi) lam_{k+1}(psi'). The anti-diagonal entries are non-zero for generic x, so the support bound
(i-1) + (l-1) <= k of L3 is sharp entry by entry. Examples: e_22^(2) = -(1-x)^4 x^3/c^4, e_23^(3) = 3 (2-x)^2 (1-x)^2 (1-2x)^2 x^3/(10 c^6),
e_33^(4) = -3 (2-x)^4 (1-2x)^4 x^3/(32 c^8), e_24^(4) = -4 (3-x)^2 (1-x)^4 (1-3x)^2 x^3/(27 c^8), all as in hvc_leaf_structure_out_kt5.txt.

CF.5 Explicit low orders, by hand, against leaf_structure.md Step 1.
 - rho. From CF.2(d), with u = 4 x tau/c^2: 2/((1 - u)^(1/4)(1 + (1 - u)^(1/2))) = 1 + u/2 + (11/32) u^2 + O(u^3). Multiplying by
   tau/t = 1 - t/2 + t^2/6 gives [t^1] rho = -1/2 + 2x/c^2 = -(1-x)^2/(2c^2) and
   [t^2] rho = 1/6 - 2x/c^2 + 11 x^2/(2c^4) = (1 - 8x + 15x^2 - 8x^3 + x^4)/(6c^4). Both agree with the table.
 - delta_1. a_1^(1) = -(1-x)^2/2 (CF.3, k = 1) and a_2^(1) = C_2 = (1-x)^2/2, so delta_1 = ((1-x)^2/2)(lam_2 - lam_1).
 - delta_2. a_3^(2) = -C_3 = -(2-x)^2 (1-2x)^2/(12 c^2) and a_2^(2) = 2 C_2 (1-x)^2/(2c^2) = (1-x)^4/(2c^2). For a_1^(2) directly:
   zeta/t = 1 + (1/c^2 - 1/2) t + (1/6 - 1/(2c^2) + (1-x^2)/(2c^4)) t^2 and 2F1(x/c, x/c; 1; tau) = 1 + (x^2/c^2) t + (x^2 (1+2x)^2/(4c^4) - x^2/(2c^2)) t^2.
   This gives [t^2] (zeta/t) 2F1 = (2 - 4x + 3x^2 - 4x^3 + 2x^4)/(12 c^4), so a_1^(2) = -(2 - 4x + 3x^2 - 4x^3 + 2x^4)/(12 c^2). All three agree with
   the table, and their sum vanishes (sum rule).
 - delta_3. Top a_4^(3) = C_4 = (3-x)^2 (1-x)^2 (1-3x)^2/(36 c^4); subleading a_3^(3) = -3 C_3 (1-x)^2/(2c^2) = -(2-x)^2 (1-x)^2 (1-2x)^2/(8 c^4).
   Both agree with the table; a_2^(3) and a_1^(3) are checked by S4_checks.pl.
 - Loops, k = 1. From CF.4 (or the Step 8 sum with g_1 = ((1-x)^2/(2c)) (lam_1 - 2 lam_2) and Theta_2 = (1-x)^2 x^3/(2c^4)):
   e_11^(1) = -(1-x)^2 x^3/(2c^2), e_12^(1) = 2 (1-x)^2 x^3/(3c^2), e_22^(1) = 0. For k = 2: e_13^(2) = -3 (2-x)^2 (1-2x)^2 x^3/(16 c^4) (corner) and
   e_22^(2) = -(1-x)^4 x^3/c^4 (anti-diagonal). All agree with the table.

CF.6 The legs.
 (a) E_11 legs. [t^k] rho f_{1,b} = Sum_{k1+k2=k} (rho^-1)_{k1} Sum_{i,l} e_il^(k2) (l^b/c^(b+1)) lam_i(psi). This has no lam_0, and its top
     coefficient (lam_{k+1} at t^k) is e_{k+1,1}^(k)/c^(b+1) = ((k+1)^2/(k+2)) x^3 a_{k+1}^(k)/c^(b+3), for all k and b. At t^0 the leg is
     -(x^3/2) c^(-b-1) lam_1, as in leaf_structure.md (ii).
 (b) E_12 at t = 0 (rank one). E_12(psi, w)|_{t=0} = lam_1(psi) chi(w) with
        chi(w) = -x^3 r_{0,2}(-w; 0)/(c + w) = -x^3 h_w(w)^2/(c + w).
     Proof: at t^0, (9.7) has only j = 1 with S_1^(0) = (c, x) = c (1, x/c), and by (N4) <(c, x), b> = x^3 b_0 for every b. So
     [t^0] rho E_12 = -x^3 r_{0,2}(-w; 0) lam_1(psi) lam_1(-w), and [t^0] E_12 = [t^0] (rho E_12). Then use 10.6(a); (I5) is not needed. Low
     orders: chi_0 = -x^3/c, chi_1 = -x^2 (1 - 5x + x^2)/(6c^2), chi_2 = -x (1 - 10x + 63x^2 - 10x^3 + x^4)/(72 c^3), as in the tables.
 (c) E_12 top coefficients (all k, b). Assume 10.6(b) (Lemma 10.5 or (I5)). Then
        [lam_{k+1}] [t^k] rho f_{2,b} = (k+1)^2 (x^3/c^2) a_{k+1}^(k) [w^b] ( h_w(w)^2 / (c + (k+1) w) ).
     Proof: in (9.8) only j = k+1 and k1 = k contribute, giving -(k+1) <S_{k+1}^(k), [w^b](v_{2,0}(-w)/(c + (k+1) w))>. By 10.6(b),
     v_{2,0}(-w) = h_w(w)^2 (1, x/c). By (N4), <S, (1, x/c)> = (x^3/c) S_0 for every S, and S_0 = -((k+1)/c) a_{k+1}^(k) by (9.1). For k = 0 this
     is (b). Explicitly:
        b = 0: (k+1)^2 x^3 a_{k+1}^(k)/c^3;
        b = 1: (k+1)^2 x^2 a_{k+1}^(k) (1 - (6k+5) x + x^2)/(6 c^4);
        b = 2: (k+1)^2 x a_{k+1}^(k) [ (1 + x + x^2)^2 - 12 (k+1) x (1 + x + x^2) + 72 (k+1)^2 x^2 ]/(72 c^5).
     For k = 0, 1, 2 the bracket is 1 - 10x + 63x^2 - 10x^3 + x^4, 1 - 22x + 267x^2 - 22x^3 + x^4 and 1 - 34x + 615x^2 - 34x^3 + x^4, as in the
     tables. The lower coefficients [lam_j][t^k] rho f_{2,b} with j <= k involve v_{2,k2} for k2 >= 1 (9.8). They are determined by Remark 10.6',
     but no closed form in b is claimed. In particular nothing forces the rank-one pattern of the lam_1 coefficients seen at k = 1, 2
     (leaf_structure.md (ii)), and it fails at k = 3.


## S4.5 Verification (Perl + Math::BigRat, no Mathematica)

leaf_proof/S4_checks.pl (unchanged; it calls zeta/t "omega"). At a rational x the script computes the following exactly, as truncated power
series in t (to t^11): tau; s* (from the critical equation); U_1; dphi in three ways (Int (1 - x s*) dt, c log(c s*) + (x-1) log U_1, and the
definition); zeta; A_j for j <= 12 (CF.1); rho in three ways (CF.2(d), the V2 formula x^3 (L_1 - L_2)/(x - c L_2) for Delta_1, and tau/(t alg));
Theta in four ways (CF.2(e) and the V8 definition); Delta_2 in two ways (10.1); l(t) and S_j (CF.2(f)); e_il for i, l <= 7 (CF.4); and h_w(w)^2.
It parses the Mathematica tables hvc_leaf_structure_out.txt, hvc_leaf_structure_out_kt5.txt and hvc_leaf_proof_check_out.txt and evaluates them
at the same x. Run at x = 3/7, 2/5, 5/3, 7/2 (output S4_checks_out.txt): 2592 exact comparisons, no mismatch. A further run at x = -1/3 and 11/5
gave 1314 comparisons, no mismatch. The comparisons, as counts over the four main values of x:
    critical equation of s*                                                      48
    dphi: three expressions agree, and agree with the printed expansion to t^7     128
    sum rule CF.2(a), t^0..t^11                                                   48
    partial fractions: [t^k](z G) computed directly from the KEY FORMULA at 13-14
      rational z equals c delta_k0 + Sum_j A_j[k]/(c + j z), k <= 11              636
    delta_k table (k <= 5, including lam_0 = 0)                                  108
    top coefficient CF.3 (k <= 11)                                                48
    subleading coefficient CF.3 (k <= 11)                                         44
    rho: three expressions to t^10, and the table to t^5                         112
    Theta: four expressions to t^10, and the printed expansion to t^7            164
    Delta_2: (10.1) both ways to t^11, and Delta_2(0) = x^3/c                     52
    sigma: l(0) = (1+x^2)/(2c), S_j = O(t^(j-1)), S_1(0) = (c, x)                 72
    loops: every entry of the restricted-ansatz table, k <= 5 (CF.4)             136
    loops: e_il^(k) = 0 outside (i-1)+(l-1) <= k, k <= 5                          88
    loops: corner, k <= 6, and anti-diagonal, i, l <= 5 (CF.4')                  28 + 60
    E_11 legs: every coefficient of every fitted j = 1 leg, k <= 5 (CF.6(a))     524
    E_12 legs: lam_0 = 0; t^0 values chi_b, b <= 7; top coefficients CF.6(c)
      for every fitted (b, k), k <= 5                                            132 + 32 + 132
These are checks of the closed forms and of their consistency with the Mathematica fits from the stored z^9 factors. They are not part of the
proofs. The partial-fraction test checks the residue computation of CF.1 independently of the tables.

leaf_proof/S4_checks_r1.pl (revision 1; output S4_checks_r1_out.txt, 163 s). It checks the parts rewritten in revision 1 without using any closed
form of the rows: the rows at U_2 are built from the code's definitions (the Wick row wickZ in the W_k(U_2(t)), computed as the Gaussian
moments (D1); the twist h_w^-1; the H-row of rmatrixCore; Delta_i as the literal quadratic form 1/eta(e_i, e_i); the edge formula (1.1)). Counts:
    Part A, at (x, tau) = (2, 5/8), (1/3, 7/12), (-3, -1/4) (t != 0; delta = 2, 1, -1):
      roots of LeqU at U_1, U_2; s_2 critical and not in {0, 1, 1/tau}; W_1(U_2) = 0                                12
      Step A of Lemma 10.5: the code's Wick row at U_2 equals U_2 Lnorm[phi_tau, (1 - tau s)^-1; s_2], z^0..z^5         18
      W_2(U_2) = U_2^2 f2/tau^2 = x delta s*^2, f2 = delta/(x s_2^2 U_2^2), Delta_2 = x^2 W_2(U_2) = x^3 delta s*^2 (10.13)    12
      (1 - tau) s_2/(1 - tau s_2) = 1 - x s_2 = L_2 - mu, L_2 = x s* (10.14)                                              3
    Part B, t-series at x = 3/7, 5/3, -5/2, 2/9:
      roots; W_1(U_2) = 0; L_2 = x s*; Delta_2 = x^3 delta s*^2 = x^2 W_2(U_2) = x^2 U_2^2 f2/tau^2; Delta_2'(0) = 0   280
      f0' = L_2 - mu, in the form d/dt [phi_tau(s_2) + c log tau] = L_2 - x                                          32
      Row 2 at U_2 (Lemma 10.5), z^1..z^6, t^0..t^6                                                                 168
      10.6(a): r^pt_{0,2}(z; 0) = h_w^-1, r_{0,2}(z; 0) = h_w^-2, z^0..z^6                                           28
      10.6(b): r_{1,2}(z; 0) = (x/c) r_{0,2}(z; 0), and [t^1] r_{0,2} = 0, z^0..z^6                                    28
      E_22 at t = 0 from (1.1) against 10.6(c), m + n <= 5                                                          84
      E_22(0,0) = -x^2 (1 + x + x^2)/(6c^2); the 2 x 2 minor of Remark 10.6.1, non-zero                                8
      CF.2(d): rho from the quadratic form of Delta_1 against the branch-free form, t^0..t^6 (c < 0 at x = -5/2)      32
      K_2 = 1 (Remark 10.6''): <v_2(z), v_2(-z)> = Delta_2, z^0..z^6, t^0..t^6                                         196
901 exact comparisons, no mismatch. The value -(x^3/2)/c^(m+n+2) of the note's F3 agrees with E_22 at 0 of 21 entries for each of the four x;
a separate run at x = 1 gives 1 of 21 and the minor 3/4096. The referee's independent recomputation (report of 2026-10-31; 7752 comparisons
at six values of x, including the orthogonality, K_1 = K_2 = 1, (9.1) and (9.8) against EdN) agrees with every claim of S4 that it tests.


## S4.6 Status of the claims (inputs (I1)-(I9) as stated in S4.0)

 - No input from S1-S3 beyond the elementary algebra of V2 (re-derived in 10.3): the edge formula (1.2); 10.1-10.4 (analyticity of the regular
   rows, which is (R); Delta_2 = x^3 delta s*^2 with Delta_2(0) = x^3/c and Delta_2'(0) = 0; Dil_2; E_22); the analysis of the conjugate encoding
   in Remark 10.6.1(2), which uses only rc_lib.m.
 - Given 10.2, (D1)-(D3), S1 Lemma 2.5 and S1 Prop. 6.1(1), (3), (4) (Laplace's method for real data, Stirling's series with remainder,
   rationality in x, uniqueness of asymptotic power series): Prop. 10.6(a).
 - Given (I1)-(I4): Lemma 9.1; Prop. 9.3 (L2 for j = 2 to all orders, and the regularity of E_12), which uses about v_2 only (R); Prop. 9.4 (the
   leaf form (9.7)); Theorem CF.1; CF.2; CF.3 (top and subleading coefficients for all k); CF.6(b) (E_12 at t = 0).
 - Given in addition (I6): Cor. 9.6 and its use in Theorem R (ii); Theorem CF.4 and Cor. CF.4' (corner and anti-diagonal for all k); CF.6(a).
 - Lemma 10.5 (Row 2 at the regular idempotent, to all orders in z; this is S2's open point O1). Its inputs: 10.2 (the code's rows and H-row at
   U_2); V1 at U_2 (S2 Prop. S2.3(c), by hand) and the tensor rule (S2 Prop. S2.4(a) at U_2, by hand), i.e. (10.13); L_2 = 1/(tau s_2) (S2 Lemma
   S2.2) and the envelope formula, i.e. (10.14); S1 (1.1) (the Euler coordinate); S1 Prop. 2.4 (P1, normalised form (iii)), Prop. 2.3 (P2,
   normalised form (2.1), applied twice), Prop. 2.2 (P3), the amplitude identity S1 (4.1) and (D1)-(D3), all at the nondegenerate critical point
   s_2(tau) for tau near some tau0 != 0 (S1 states them for an arbitrary analytic family of nondegenerate critical points); S2 Lemma S2.6
   (= V7, proved by hand).
 - Corollary 10.5' (the orthogonality (I5)), given Lemma 10.5, S2 Prop. S2.7 (the QDE at U_1, from S1's closed form) and S2 Prop. S2.8(c)
   (self-adjointness of A); this is S2 Prop. S2.10.
 - Given Lemma 10.5 and 10.6(a): 10.6(b) (both proofs), hence the regularity of E_21; 10.6(c) and Remark 10.6.1(1); CF.6(c); Remark 10.6''
   (K_2 = 1, which with K_1 = 1 and Corollary 10.5' gives V6 for all i, j).
 - Given (I5): Prop. 9.5 (E_21 = E_12^T); Remark 10.6'.
 - F1 and F3 are complete (see 10.7), given all of the above.
 - Verified only (S4.5): the agreement of the closed forms with the tables of leaf_structure.md (k <= 5), at six rational values of x; Row 2 at
   U_2, 10.6(a)-(c), the minor of Remark 10.6.1, K_2 = 1 and CF.2(d) from the code's definitions, at four rational values of x and at three
   points with t != 0.

## S4.7 Mathematica checks requested (not run; licence)

 R1. CF.1 symbolically: [t^k] A_j against the lam_j coefficients of [t^k](rho z r_{0,1}) from hvc_factors_g4_Bm.mx, k <= 6, exact in x; and the
     sum rule CF.2(a) to t^12.
 R2. CF.4 symbolically against EdN[1,1] (m + n <= 8, t^0..t^6) and against loopResR in hvc_leaf_structure_data_kt5.m.
 R3. (9.7): build -Sum_j j <S_j(t), v_2(-w)> lam_j(psi) lam_j(-w) from the stored ratCN[a,2] (conjugate encoding, i.e. evaluated at U1) and the
     S_j of CF.2(f), and compare with EdN[1,2][a, b] (evaluated at U1), a + b <= 8, t^0..t^5. (The referee's Perl check covers a + b <= 5.)
 R4. 10.6(a),(b) and Remark 10.6': the t-series of ratCN[0,2][m] and ratCN[1,2][m], evaluated at U1, at t = 0 against h_w^-2 and (x/c) h_w^-2
     (m <= 9), and [t^1] ratCN[0,2][m] = 0. Also (10.1), Delta_2 = x^3 delta s*^2, exactly.
 R5. (A confirmation; the cause of the discrepancy with the note is identified in Remark 10.6.1(2).) (a) With hvc_factors_g3_Bm.mx (or g4) loaded:
     zeroC[cadd[EdN[2,2][m, n], cscale[-1, cconj[EdN[1,1][m, n]]]]] for m + n <= ZMAX - 1 (zeroC as in hvc_vertex_identity2.wls), and the same for
     EdN[2,1] against cconj[EdN[1,2]]. (b) Table[val0[EdN[2,2][m, n], U1], ...] against (x^3/c)[z^m w^n](1 - h_w(z)^2 h_w(w)^2)/(z + w), exact in x,
     m + n <= 8; in particular E_22(0, 0)|_{t=0} = -x^2 (1+x+x^2)/(6c^2). (c) Line 14 of hvc_vertex_identity.wls re-run with U1 in place of U2.
 R6. E_21(m, n) = E_12(n, m) exactly in Q (m + n <= 8), and no negative t-powers in DilN[2], EdN[2,2], EdN[1,2], EdN[2,1] (the F3 table), all
     evaluated at U1.
 R7. (Optional.) V5 for i = 2 (Row 2 at U2) and V6 for i != j extended to z^9, t^0..t^8, symbolic in x: a confirmation of Lemma 10.5 and
     Corollary 10.5' beyond the Perl range (z^6, t^6).

## S4.8 Open points and cross-references

 - S2_metric_qde.md (revised), s.S2.10 ("Where the orthogonality is used"), S2.12 and (O1), now say what S4 needs: the orthogonality is used for
   E_21(b, a) = E_12(a, b) and for the regularity of E_21, and Row 2 at the regular idempotent is needed for the graph sum, through S4. Lemma 10.5
   is the argument sketched in S2.10, written out: (A2) is Step A (10.15), in the normalised form; (B2) is Steps B-D ((10.16)-(10.18); Remark 10.5.1);
   Corollary 10.5' is S2 Prop. S2.10. S2.12 and (O1) already refer to Lemma 10.5 for the proof; S2's remark that K_2 = 1 beyond
   K_2 = 1 + O(z^2) needs a further argument is answered by Remark 10.6''.
 - Proposed one-line addition to outline.md, Step 9 (not applied; outside this file): "E_21: E_21(b, a) = E_12(a, b) by the orthogonality of
   Step 6, which needs the QDE at the regular idempotent (S4, Prop. 9.5, Lemma 10.5, Corollary 10.5'); already the regularity of E_21 needs the
   t^-1 part of it (S4, 10.6(b))."
 - Proposed correction of conifold_regularity_note.md F3 (not applied; outside this file): replace "E_11(z, w) = E_22(z, w) = -(x^3/2) phi(z) phi(w)"
   by "E_11(z, w) = -(x^3/2) phi(z) phi(w) and E_22(z, w) = (x^3/c)[1 - h_w(z)^2 h_w(w)^2]/(z + w), which is not of rank one", and in
   hvc_vertex_identity.wls, line 14, evaluate EdN[2,2] at U1 (Remark 10.6.1). Theorems A and B and the Reduction Lemma of the note are unaffected.
 - Notation for the merged text: zeta(t) of CF.1 is t times the omega of S3 s.1.6 and of S4_checks.pl; kappa stays reserved for the phase
   constant of S1 (1.2) (kappa_0, kappa_1 for S2.5); omega(z, tau) is S2's (Lemma S2.6, used in 10.5). S3's b_j(t) is A_j(t).
 - K_2 := <v_2(z), v_2(-z)>/Delta_2 = 1 to all orders (Remark 10.6''), so V6 holds for all i, j. It is not needed: E_22 and Dil_2 are regular
   whatever K_2 is.

## S4.9 Changes in revision 1 (referee report of 2026-10-31)

 1. Lemma 10.5 (minor issue on its proof): rewritten as a parallel of S1 Prop. 4.1 at s_2(tau), tau near tau0 != 0, in four steps. Step A is
    S1 Prop. 3.1 at U_2 (normalised P1). Step B is the conjugation identity (10.16), which uses V1 at U_2 and the tensor rule (10.13) and the
    envelope formula (10.14). Step C is S1 (4.1), P2 twice and P3 at s_2. Step D is S2 Lemma S2.6 in conjugated form. Remark 10.5.1 relates the
    proof to the Lfull form of S1 Prop. 2.4(ii) (prefactor Q^(x/z) tau^(c/z), Jacobian -tau (1 - tau s)^-1, phase kappa(z)) and to S2.10 (A2),
    (B2). Inputs listed in S4.0 (I8), (I9) and S4.6.
 2. Orthogonality (minor): the argument with the formal symbol e^((u_1 - u_2)/z) is replaced by the computation of S2 Prop. S2.10 (Corollary 10.5').
 3. 10.6(a) (minor): the analytic input is now S1 Lemma 2.5 on I = R in y = log u, with F(y) = y - c log(1 + e^y), and S1 Prop. 6.1(1), (3), (4).
    The sign identity Lnorm[-f](z) = Lnorm[f](-z), (10.19), shows that the Wick series in z is the Laplace series in 1/lambda = -z. The z^1
    coefficient is checked by hand.
 4. E_22 and the note (minor): Remark 10.6.1 identifies the cause of the note's F3 claim (conjugate encoding; EdN[2,2] = cconj[EdN[1,1]] exactly,
    evaluated at U2 in hvc_vertex_identity.wls) and proves that E_22 at t = 0 is not of rank one (a 2 x 2 minor). R5 becomes a confirmation.
    The corrections of the note and of the script are proposed in S4.8 (they are outside this file).
 5. Branches (cosmetic): (I7) fixes delta(0) = c; CF.2(d) gives rho and alg in a form without branch ambiguity.
 6. S4.8 (cosmetic): the cross-reference to the revised S2 replaces the out-of-date bullet; a line for outline Step 9 is proposed (not applied).
 7. S4.6 (cosmetic): the inputs of Lemma 10.5 are complete (V1 at U_2, S1 (4.1), P3; V7 is cited as S2 Lemma S2.6).
 8. "Use in Theorem R" (cosmetic): the bullets are now organised by the other end of the edge (regular; singular of genus 0 or 1; singular of
    genus >= 2; self-loop) and no longer overlap.
 9. Notation (cosmetic): kappa(t) is renamed zeta(t); alpha_j of CF.2(f) is renamed A~_j (alpha is a Gauss parameter in 10.5); beta_j of
    Prop. 9.3 is renamed R_j.
10. Wording (cosmetic): Summary (1) and the paragraph after Prop. 9.5 (without (I5) the rationality argument does not apply; with (I5) the leg
    is f_{2,b}; "f constant implies f = 0" uses r_{1,2}[0] = L_2 and L_2(0) = x/c).
11. New check script S4_checks_r1.pl (901 exact comparisons, S4.5).
12. New (not requested, not needed): Remark 10.6'' proves K_2 = 1 to all orders from Lemma 10.5 and 10.6(b), so V6 holds for all i, j. This
    answers the remark in S2.10 / (O1) that K_2 = 1 beyond O(z^2) would need a further argument.
