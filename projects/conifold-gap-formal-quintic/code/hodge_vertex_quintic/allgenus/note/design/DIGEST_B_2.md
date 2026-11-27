# DIGEST_B_2 (key read_B), s.2 part 1: statements S1-S19 (PF4/PF5 and (CYC); Theorem U; the model rows)

Continuation of `DIGEST_B.md` (abbreviations there). Format: label proposal (prefix by topic: `cyc` the cyclic systems and (CYC),
`gau` Theorem U and the gauge, `mod` the model, `van` vanishing solutions and (G), `sym` symplecticity, `res` residue and degree
structure, `hat` hatted transport), type, statement in final form with hypotheses and scope, source, status, corrections, dependencies.

## s.2 Statements

### (A) PF4, PF5 and (CYC)

**S1** `lem:cyc-pf5` (lem). As operators with coefficients in `Q[u]`: `D^5 - q prod_(k=1..5)(5D + k) = D o PF4`; hence, `z` being
constant, `PF5 = z^5 D PF4 - 1`. Uses `D o q = q o (D + 1)`, i.e. `q (5D + 5) = 5 D q`.
Source: CYC Cor. proof; PB s.2.2; CF2 s.8.2. Status: PROVED (elementary). Deps: D2.

**S2** `lem:cyc-frob` (lem; "Frobenius structure at the conifold"). (i) PF4 maps `Q[[u]]` into itself; `PF4 u^n = Sum_(s=-3..1) c_s(n)
u^(n+s)` with `c_(-3)(n) = n(n-1)^2(n-2)` and `c_s(n) = 0` whenever `n + s < 0`, `n <= 2`; so the coefficients of `u^0, u^1, u^2` are
free and determine a power-series solution (exponents 0, 1, 1, 2). (ii) For a power-series solution `g`, `r(g) := PF4(g l) - l PF4(g)`
lies in `u^(-3) Q[[u]]`, and `g l + h` (`h` in `Q[[u]]`) is a solution iff `[u^m] r(g) = 0` for `m = -3, -2, -1`. (iii) On
`span(f_0, f_1, f_2)` the obstruction matrix (rows `u^-3, u^-2, u^-1`) is `[[-2,0,0],[17/5,0,0],[-7/5,-7/5,2]]`, kernel `Q (0, 1, 7/10)`:
the log solution is unique up to scale and up to adding power-series solutions, and `g = f_1 + (7/10) f_2` is canonical up to scale.
Equivalently, with `PF4 u^x = Sum_j P_j(x) u^(x+j)`, `P_(-3) = x(x-1)^2(x-2)`, `P_(-2) = -4x^4 + 10x^3 - (47/5)x^2 + (17/5)x`:
`O_(-3) = -2 b_0`, `O_(-2)` automatic (`P_(-3)'(1) = 0`), `gamma = -P_(-2)'(1)/P_(-3)'(2) = (7/5)/2 = 7/10`; and `gamma` is the unique value
with `B(f_0, f_1 + gamma f_2) = -7/5 + 2 gamma = 0`. Source: CYC s.1; CF s.2.2-s.2.3 (hand); CA s.3 (closed form). Status: PROVED.
Correction: CYC s.1 "holds exactly for `g = f_1 + (7/10) f_2` (computed)" claimed more than its code searched (`span(f_1, f_2)` only); replaced
by the 3x3 matrix (RS2 s.5 item 1 = CF-E1 + CA-K3). The uniqueness is needed only for the words "the vanishing period", not by S9.
Deps: D2, D3.

**S3** `lem:cyc-const` (lem; CYC "Lemma 0"). `{f in R : Df = 0} = Q`. Hypothesis used: `D l = 1 - 1/u` is not in `D(K)` (control: with
`D l' = 1 - u = D(-u)`, `l' + u` is a non-constant constant; CF s.2.1). Source: CYC s.1; CF s.2.1. Status: PROVED. Deps: D3.

**S4** `lem:cyc-adj` (lem; C1). `PF4* = PF4`. Criterion (CF s.3, CA s.3): an order-4 operator is formally self-adjoint iff
`A_3 = 2 D A_4` and `A_1 = D A_2 - D^3 A_4`; PF4 satisfies both. Source: CYC s.2; CF s.3; CA s.3. Status: PROVED. Deps: D2.

**S5** `lem:cyc-form` (lem; C2). (i) `g L f - f L* g = D P(f,g)` for all `f, g` (identity of differential polynomials, any order-4 `L`).
(ii) For `L = PF4` and `f, g` in S: `B(f,g) := P(f,g)` lies in Q; B is bilinear and alternating. (iii) `P(f,g) = Sum_(a,b) D^a f P_ab D^b g`
with `P_ab = 0` for `a + b >= 4` and `P_(3-j,j) = (-1)^j u` (explicit matrix D17; `B(f,g) = (2/5)(1-u) W_01 + (1-u) W_02 + u W_12 - u W_03`,
`W_ab = D^a f D^b g - D^b f D^a g`); `det P = u^4`. (iv) `Bmat = Wr^T P Wr`; in the basis `(f_0, f_1, f_2, y_3)`
`Bmat = [[0,-7/5,2,-9/10],[7/5,0,0,-1],[-2,0,0,0],[9/10,1,0,0]]`; `Pf(Bmat) = -B(f_0,f_2) B(f_1,y_3) = -(2)(-1) = 2` (as
`B(f_1,f_2) = B(f_2,y_3) = 0`), `det Bmat = 4`; hence B is non-degenerate and `det Wr != 0`; in fact `det Wr = -2/u^2` (Abel's identity
`D det Wr = -(A_3/A_4) det Wr` and S3). Source: CYC s.3; CF s.4.1-s.4.4 (hand values); CA s.3. Status: PROVED. Correction: "computed
exactly from the solutions" replaced by the hand values (RS2 s.5 item 2 = CF-E2). Deps: S3, S4, D17.

**S6** `lem:cyc-transv` (lem; C3, "Griffiths transversality"). For every Q-basis `(y_i)` of S: `Sum_(i,j) (Bmat^(-1))_ij D^a y_i D^b y_j =
(P^(-1))_ab`; it vanishes for `a + b <= 2`, and `(P^(-1))_(0,3) = 1/u`. In particular `T_0 = T_1 = T_2 = 0`, `T_3 = 1/u`.
Explicit `P^(-1)` (CF s.5): `[[0,0,0,1/u],[0,0,-1/u,-(1-u)/u^2],[0,1/u,0,(2/5)(1-u)/u^2],[-1/u,(1-u)/u^2,-(2/5)(1-u)/u^2,0]]`.
Identity in `Frac(R)`. Source: CYC s.4; CF s.5. Status: PROVED. Deps: S5.

**S7** `lem:cyc-iso` (lem; C4, "formal Picard-Lefschetz"). `B(f, g) = 0` for every power-series solution `f`; in particular
`B(y_0, g) = 0` in every conifold frame; `B(y_3, g) = 1`. So `span(y_0, g)` is isotropic for every power-series `y_0`: the frame
condition is vacuous. Source: CYC s.5; CF2 s.6; CA s.3. Status: PROVED. Remarks: proves FQ2 s.1 "Cf preserves the Wronskian relation" and
"`B(g, analytic) = 0` automatically" (RS2 s.5 item 44; CA-K5). Deps: S5, D3.

**S8** `lem:cyc-sg` (lem; C5, "special geometry"). Let `y_0` be a power-series solution not proportional to `g`, `y_1 := g`, `t := y_1/y_0`,
`th := (Dt)^(-1) D`. There are `y_2, y_3'` in S with `B(y_0, y_3') = B(y_1, y_2) = 1`, `B(y_0, y_2) = B(y_1, y_3') = B(y_2, y_3') = 0` such that,
with `phi_i := y_i/y_0` (`phi_0 = 1`, `phi_1 = t`, `phi_3 := y_3'/y_0`) and `Y := th^2 phi_2`: (i) `th^2 phi_3 = -t Y`, hence
`th^2(y/y_0)` lies in `Y span_Q(1, t)` for every `y` in S; (ii) `Y = 1/(u y_0^2 (Dt)^3)`. Source: CYC s.6; CF2 s.7. Status: PROVED.
Corrections: (a) add (CF-E7, RS2 s.5 item 7): `tau_1 = tau_2 = 0`, `tau_3 = 1/(u y_0^2 Dt^3)` (CYC's `tau_k := Sum (Bmat^-1)_ij phi_i
th^k phi_j`) hold for every unit `y_0` and every `t` with `Dt` a unit, since S6 is basis-free; isotropy enters only through the adapted
basis, which exists iff `B(y_0, y_1) = 0`. (b) CA-K6 (RS2 s.5 item 9): `Y` is normalised by the Darboux form of B; the Yukawa coupling of
FQ2 s.1 is `5Y` (`eta = 5 delta`); the row coefficient `c_1 = Y Dt` carries no 5 (a 5 there breaks S10); the `u` is `A_4`. (c) Rename
CYC's `f_i := y_i/y_0` (clash with the Frobenius basis; CF3 E6). Deps: S6, S7.

**S9** `thm:cyc-main` (thm; "Theorem CYC"). Let `y_0` be a non-zero power-series solution of PF4 not proportional to `g` (so `t = g/y_0`
is non-constant), `g` taken at any non-zero scale. Then, as operators with coefficients in `Q((u))`,
`PF4 = M_cyc(y_0, g) = (1/y_0) D (1/Dt) D (u y_0^2 (Dt)^2) D (1/Dt) D (1/y_0)`. More generally the identity holds for every pair
`(y_0, y_1)` in `S cap Q((u))` spanning a B-Lagrangian plane (these are exactly the planes containing g). Scope: `y_0(0) != 0` is NOT
used (every non-zero element of K is a unit); it is the condition for the Theorem-U unit shape `val(c) = (0,-1,0,0,0)` (e.g.
`(0,-3,0,1,1)` for `y_0 = f_1`, `(-2,-1,-2,2,2)` for `y_0 = f_2`). The large-radius frame `(I_0, I_1)` and the ring N of fq_ring are not
covered and not needed. `M_cyc` is invariant under `y_0 -> mu y_0`, `g -> lambda g`, `y_1 -> y_1 + a y_0`. In the programme's frames
(D4) this is `PF4 = (1/varpi_0) D (1/Dt_c) D (u varpi_0^2 (Dt_c)^2) D (1/Dt_c) D (1/varpi_0)`.
Source: CYC s.0, s.7; scope RS1 A1 (CA-K2 adopted over CF-E5); CF2 s.8.1 (scale invariance); CA s.2.5. Status: PROVED; AI-refereed twice
(CF, CA), human check pending (MS5 s.5 item 1.2). Corrections: RS2 s.5 item 5 (scope), item 6 (notation: `M` -> `M_cyc`, `tau` -> `sigma`,
`P` jet matrix -> `Pi_ab`, `y_3` vs `y_3'`). Deps: S2-S8.

**S10** `cor:cyc-system` (cor; "the cyclic system is PF5"). Let `c = (Dt, L^5/(y_0^2 (Dt)^2), Dt, y_0, y_0)`, `L^5 = 1/u`, with `(y_0, t)` as in
S9 and `c_a` invertible. In any differential extension of `(Q((u)), D)` in which `z` is constant and the `c_a` are invertible, the system
`zD Y_a = c_a Y_(a+1)` (`a` in Z/5) holds iff `Psi := y_0 Y_0` solves PF5 and `Y_(a+1) = z D Y_a/c_a` for `a = 0..3`. The rings used by the
programme: `Q(zeta)((s))((z))` for the WKB rows (`zD + zeta^i L` is the image of `zD` under the automorphism `D -> D + zeta^i L/z`);
`L((u))`, `L = Q((z))`, for `P_can Q^m` and `Y^C_Q`; the hatted series via `z = s^4 zhat`. The identities lie in `Q[[u]][D]` for frames
with `y_0(0) != 0`, so they act on every differential `Q[[u]]`-module with z constant (CF version). Key step: `(1/c_4) D (1/c_3) D (1/c_2)
D (1/c_1) D (1/c_0) D o (1/y_0) = (1/y_0) D M_cyc`, so the system reads `z^5 D M_cyc Psi = Psi`; with S1 and S9, `PF5 Psi = 0`.
Both directions are used: "cyclic system => PF5" in P(b) step 5 (PB s.2.8), the glue (PB s.3), POLY Lemma 2 for `Phi^can`, EV Lemma N
(PS-C9); "PF5 => closing relation `4 -> 0`" for the actual rows (UR s.1, s.4; hence R0, S20) and for PB's `Q^SP`. The existence of
log-free WKB rows at all orders is a separate fact, derived inside the glue (PS s.4.2; RS1 A4), not part of this corollary.
Source: CYC s.0 Cor.; CF2 s.8.2; CA s.2.8; rings RS2 s.5 item 8 (CF-E8 + CA-K4). Status: PROVED. Deps: S1, S9.

**S11** `rmk:cyc-robust` (rmk; for the writers' choice). (a) In the programme frame `(alpha, beta)`, the project's own row coefficients
pushed through the conifold map equal `1/c_a` (`Cf(kinv[a]) = 1/c_a`, `Cf(C_1) = D t_c`; CA s.2.2), i.e. S10 is the statement the
programme uses (CA s.0). (b) Theorem U's unit shape does not detect a wrong `gamma` or a non-isotropic pair: the identification of the
Theorem-U system with the quintic's rests on `gamma = 7/10` exactly (S2) and S7 (CA-K3; RS1 A10). (c) The order of the `c_a` matters
(`c_0 = c_2`, `c_3 = c_4` make `M_cyc` palindromic; swapping `c_1, c_2` breaks the closing relation; CA s.2.3). (d) Deformations break
exactly the lemma they should (CA2 s.5 table). (e) SKETCH (CA2 s.5 Remark): the same proof works for every self-adjoint order-4 operator
with `A_k` in `Q[[u]]`, `A_4 = u` times a unit, `det Bmat != 0` (middle factor `A_4 y_0^2 Dt^2`); the Wronskian step is not proved in general
(RS2 s.6.3). Not used. (f) Classical special geometry (CYC s.9 names Strominger, Candelas-de la Ossa, arXiv:hep-th/9204035, NOT read).

### (B) Theorem U

**S12** `lem:gau-hyp` (lem; "the two systems satisfy the hypotheses of Theorem U"). In every frame `(alpha, beta)` of the programme family:
`b_0 = b_2 = t_c'` and `b_3 = b_4 = -varpi_0/(1-u)` are units (`t_c'(0) = varpi_0(0) = 1`), `u b_1` is a unit; the pulled-back model
`bt = (v', -v'/v, v', -v', -v')` has the same shape; `prod b_a = -1/(u(1-u)^5)` and `prod bt_a = -v'^5/v`, equal iff
`v'^5/v = 1/(u(1-u)^5)` iff (with `v = u + O(u^2)`) `v^(4/5) = u^(4/5)/tau`, i.e. `v = u tau^(-5/4)` (D7). Both systems are self-dual and
have residue `r = rt = -1`. Necessity (UR s.3 Rmk (c)): `P_0` must conjugate `A` into `At`, so equal products are necessary; for a change
of variable `v(u)` this is the equality of the canonical one-forms, which fixes `v` given `v = u + O(u^2)`.
Source: PB s.2.8; PS s.3; UR s.3 Rmk (c), s.4; GI s.1; CA2 s.4 (unit shape iff `y_0(0) != 0`). Status: PROVED.
Corrections: "v unique given v(0) = 0" -> unique given `v = u + O(u^2)` (NFS s.4; PS s.3); the normalisations `t_c'(0) = varpi_0(0) = 1`
hold in the programme family only (CA-K2; RS2 s.5 items 27-28). CA-K7 (RS2 s.5 item 13): Theorem U itself needs no (CYC); (CYC) supplies
(i) that the actual rows solve the system to which U is applied (UR s.4) and (ii) that every solution of it solves PF5. Deps: D4-D8, D10.

**S13** `thm:gau-U` (thm; "Theorem U"). Let `A`, `At` be cyclic systems with one simple pole (D10) with `prod_a b_a = prod_a bt_a`. Then
there is `P = Sum_(k>=0) z^k P_k(u)`, `P_k` in `gl_5(O)`, `P_0` diagonal with unit entries, such that `z P' = A P - P At` (so `Y = P Yt`
maps solutions to solutions). Every such `P` with coefficients in `gl_5(O)` and `P_0` invertible is obtained from one of them by
`P -> P c(z)`, `c` in `C[[z]]^x` (`Q[[z]]^x` over Q). If one restricts to the Z/5-graded form `P_k = Sum_b p_(k,b) E_(b+k,b)`, the freedom is
`c(z^5)`. Source: UR s.3 (statement); GI s.1 (quoted over `Q[[u]]`). Status: PROVED in UR s.3; AI-audited only (NFS s.1 item 2), not in any
referee round (RS1 s.3 row (U)); human check pending (MS5 s.5 item 1.6). It is input (U) of the main theorem (MS5 s.2).
Remarks in UR: (a) the proof uses only the shape (the model is "versal" for such systems); (b) geometric reading: the lattice `+ O Y_a`
is the structure sheaf of the spectral curve `y^5 = 1/B(u)`, exponents `j/5` = ramification weights (a remark, not used); (d) computation
(evidence, s.5 E2). Deps: S14-S16, S17 (graded clause).

**S14** `lem:gau-comm` (lem; UR s.3 (ii), "the commutant"). `y := At^(-1) = Sum_a bt_a^(-1) E_(a+1,a)` is holomorphic, `y^5 = 1/B`,
`B := prod bt_a = unit/u`. `At` is regular semisimple over K, its commutant in `gl_5(K)` is `K[y]`, and `gl_5(K) = K[y] (+) im(ad At)`.
`K[y] cap gl_5(O) = O[y]`. Write `pr` for the projection onto `K[y]`. Grading: `X = Sum_a xi_a E_(a,a+m)` has degree `m` (`1 <= m <= 5`,
degree 5 = degree 0); `W_a^(m) := bt_a .. bt_(a+m-1)`, `W_a^(5) = B`; then `pr(X) = (1/5)[Sum_a xi_a/W_a^(m)] At^m =
(1/5)[Sum_a xi_a B/W_a^(m)] y^(5-m)`. For holomorphic X the coefficient of `y^j` in `pr(X)` has at most a simple pole (`j = 1..4`) and is
holomorphic for `j = 0`. Source: UR s.3 (ii). Status: PROVED (UR); several steps stated without proof (s.3 below). Deps: D10.

**S15** `lem:gau-lift` (lem; UR's "Lemma 1"; rename: clashes with GI Lemma 1). If `X` in `gl_5(O)` and `pr(X) = 0`, there is `Pi` in
`gl_5(O)` with `[At, Pi] = X`, unique up to `O[y]`. Source: UR s.3 (iii). Status: PROVED. Deps: S14.

**S16** `lem:gau-scal` (lem; UR's "Lemma 2"; rename: clashes with POLY Lemma 2). Let `D_0 = diag(d_a)`, `d_(a+1)/d_a = bt_a/b_a`, and
`Dl := D_0^(-1) D_0'`. For `Z` in `gl_5(O)` the equation `pr(c' + Dl c) = -pr(Z)` has a solution `c = Sum_(j=0..4) c_j y^j` in `O[y]`, unique
up to `c_0 -> c_0 + const exp(-Int tr Dl/5)`. Facts used: `pr((y^j)') = -(j/5)(log B)' y^j = (j/(5u) + hol) y^j`; `pr(Dl y^j) =
(tr Dl/5) y^j`; the `y^j` component is `u c_j' + (j/5 + u h_j) c_j = -u phi_j` with `h_j`, `u phi_j` in O (regular singular, exponent
`-j/5`, pivots `n + j/5 != 0` for `j = 1..4`; regular for `j = 0`). Source: UR s.3 (iv). Status: PROVED. Deps: S14.

**S17** `prop:gau-graded` (prop; "graded gauges and covariance"). Both systems are covariant under `z -> zeta z` combined with
`diag(zeta^(-a))`; hence the regular gauge can be taken graded, and for a graded gauge the constant diagonal series of S20 satisfies
`D_i(z) = D_0(zeta^(-i) z)`; the graded regular gauges are `P c(z^5)`. For a graded P, `P_11(0,z) = Sum_q p_(5q,1)(0) z^(5q)` is a unit of
`Q[[z^5]]`, so `P_can := P P_11(0,z)^(-1)` is graded (D12). Source: UR s.4 (one sentence), s.3 (statement of U); GI s.1, s.3 Cor. (a);
PS s.4.5 ("confirmed", with `h_(w_i)(z) = h_(w_0)(zeta^(-i) z)` from `w_i = zeta^i w_0`). Status: PROVED in the sources by one sentence
each (see s.3). Deps: S13, D5, D8.

**S18** `prop:mod-rows` (prop; UR "Proposition S"; rename: clashes with GI Prop. S). (a) The model system (D8) holds for
`Y = (I_0, -I_(-1), -I_3, I_2, I_1)` by (Mod1), (Mod2) ("three lines"). (b) With `Y_b = e^(U/z) sigma^(c_b-1) F_b(t)` (D18) the model system is
equivalent to `F_1 = T_(-1) F_0`, `F_2 = -T_(-2) F_1`, `F_3 = T_2 F_2`, `F_4 = -T_1 F_3`, `F_0 = -T_0 F_4`; for each colour it has a unique
solution with `F_0(0) = 1` (the recursion for `[t^n] F_0` has pivot `-4n`); `F_b(0) = (1,-1,-1,1,1)` (colour 0); the general WKB column with
the same exponential is the homogeneous one times a constant series, which is homogeneous only if constant ("uniqueness half", used by
(S2) of PB s.2.6, independent of (H): PS s.1, RS1 s.3 row P(b)). `F_0 = 1 - (3/20)t + (9/800)t^2 + (11571/80000)t^3 - (2415609/6400000)t^4
+ (47735919/128000000)t^5 + (52989974037/25600000000)t^6 - (7749624763881/512000000000)t^7 + ..`. (c) Given (H) (now a consequence of
`D = 1`, ST2 s.2 "(F_g) implies (H)"), the conifold symbol (D19) is the model: `R^_b^(i)(zhat) = zeta^(i(b+1)) F_b(zeta^(-i) zhat)`.
Source: UR s.1 "Model", s.2; LG s.0 item 4 (closed recursion `g^(a+1)_n - g^a_n = ((4n - 3 - c_a)/5) g^a_(n-1)`). Status: (a), (b) PROVED;
(c) PROVED given (H). Deps: D8, D18, D19.

**S19** `rmk:mod-limit` (rmk). (a) Scaling limit (UR s.1, PROVED): with `z = u^(4/5) zhat`, in the u-form of PF5 the terms `u th^5`,
`3 z th^4`, `1` are of order `u^0` and all others smaller by integer powers of u; the limit system is `b -> (1, -1/u, 1, -1, -1)` = the model
system; local exponents `0, 1, 1, 2, 3`. (b) Mellin-Barnes form (UR s.2; "PROVED, elementary; not machine-checked beyond the ODE"):
`Int_0^inf exp(-xi^4/4 - k/xi) dxi = (1/2 pi i) Int Gamma(s) Gamma((s+1)/4) 4^((s-3)/4) k^(-s) ds`, `k = v/z^(5/4)`; five Gamma factors with the
indicial polynomial `th(th-1)^2(th-2)(th-3)`. (c) The model is a point of the Hurwitz space of `x^4 + a x^2 + b x + c + u/(x - d)` on
`a = b = d = 0` (UR s.1; LG s.0 item 5, s.5, VERIFIED locally; arXiv:2211.12259 QUOTED). Not needed by the chain; candidates for remarks.
