# DIGEST_B_4 (key read_B), s.3: proofs, step by step

Continuation of `DIGEST_B_3.md`. For each statement of s.2: the proof as numbered steps, with where each step is written.
Marks: **[W]** written out in the source named; **[W*]** written, but in one or two terse sentences (transcribable as is); **[S]** only
stated or sketched in the sources (justification missing: a `\draftnote` unless the authors supply it); **[V]** only VERIFIED;
**[M]** missing from every source read.

## s.3 Proofs

**S1** (PF5 = z^5 D PF4 - 1). 1. `D o q = q o (D + 1)`, so `q(5D + 5) = 5 D q` [W] CYC Cor.; PB s.2.2. 2. `D o PF4 = D^5 - 5 D q prod_(k<=4)
(5D+k) = D^5 - q prod_(k<=5)(5D+k)` [W] CF2 s.8.2. 3. `z` constant: `(zD)^5 = z^5 D^5`, `5zD + kz = z(5D + k)` [W] CF2 s.8.2.

**S2** (Frobenius structure, gamma). 1. `PF4 u^n = Sum_(s=-3..1) c_s(n) u^(n+s)`, `c_(-3)(n) = n(n-1)^2(n-2)` from `A_4 D^4 + A_3 D^3` only
[W] CF s.2.2 (all five `c_s(n)` symbolic in n); closed form `P_j(x)` [W] CA s.3. 2. `c_s(n) = 0` for `n + s < 0`, `n <= 2`: PF4 maps
`Q[[u]]` to itself, triangular recursion, `u^0, u^1, u^2` free [W] CF s.2.2. 3. `r(f) = Sum_k A_k Sum_(m>=1) C(k,m) D^(k-m) f D^m l`; `D l =
-u^-1 + 1`, `D^2 l = -u^-2 + u^-1`, `D^3 l = -2u^-3 + 3u^-2 - u^-1`, `D^4 l = -6u^-4 + 12u^-3 - 7u^-2 + u^-1`; so `r(f)` is in `u^-3 Q[[u]]`
[W] CF s.2.3. 4. `g l + h` is a solution iff `[u^m] r(g) = 0`, `m = -3,-2,-1` (PF4 of a power series cannot produce these terms) [W] CYC s.1;
CF s.2.3. 5. The obstruction matrix by hand: `[u^-3] r(f) = -6f(0) + 4f(0) = -2f(0)`; `[u^-1] r(f_1) = 0 - 7/5 + 2 - 2 = -7/5`;
`[u^-1] r(f_2) = 0 + 0 + 4 - 2 = 2`; the `u^-2` entries of the `f_1`, `f_2` columns are 0, using `Df_1 = -1 + u + (19/20)u^2`, `D^2 f_1 = -1 -
(9/10)u`, `D^3 f_1 = 9/10`, `Df_2 = -2u - (17/10)u^2`, `D^2 f_2 = 2 + (7/5)u`, `D^3 f_2 = -7/5` [W] CF s.2.3; equivalently
`O_m = Sum_j P_j'(m-j) b_(m-j)` [W] CA s.3. 6. Kernel `Q(0,1,7/10)`, hence `gamma = 7/10`; equality with the isotropy condition [W] CA s.3.

**S3** (Lemma 0). 1. In K: `Df = -(1-u) f' = 0` gives `f' = 0`, `f` in Q [W] CYC s.1. 2. For `f = Sum_(k<=n) a_k l^k`, `n >= 1`, the `l^n`
coefficient of `Df` is `D a_n`, so `a_n` in `Q^x` [W]. 3. The `l^(n-1)` coefficient `D a_(n-1) + n a_n (1 - 1/u) = 0`; dividing by the unit
`1-u`: `a_(n-1)' = -n a_n/u`, impossible since a derivative in `Q((u))` has no `u^-1` term [W] CF s.2.1 (makes the division explicit).

**S4** (C1). 1. `PF4 = D^4 - 5^5 q o P(D)`, `P(D) = prod_(k=1..4)(D + k/5)` [W] CYC s.2. 2. `(q o P(D))* = P(-D) o q = q o P(-D-1)` [W].
3. `P(-D-1) = prod(-(D + (5-k)/5)) = P(D)` [W]. 4. `(D^4)* = D^4` [W]. Alternative: the two-condition criterion [W] CF s.3; CA s.3.

**S5** (C2). (i) 1. With `H = A_k g`: `D[(-1)^j D^j H D^(k-1-j) f] = T_j - T_(j+1)`, `T_j = (-1)^j D^j H D^(k-j) f`; the sum telescopes to
`H D^k f - (-D)^k(H) f` [W] CYC s.3; CF s.4.1. (ii) 2. On S: `D P(f,g) = g PF4 f - f PF4* g = 0`, so `P(f,g)` in Q by S3 [W]. 3. Alternation:
`Q := P(f,g) + P(g,f)` has `DQ = 0` as a differential polynomial; at top total order `s` the coefficient of `f_(a+1) g_b` in `DQ` is
`q_ab + q_(a+1,b-1)` and that of `f_0 g_(s+1)` is `q_(0,s)`; successively all `q_ab` with `a + b = s` vanish, a contradiction, so `Q = 0` [W] CYC
s.3; CF s.4.2. (iii) 4. Order count (`k - 1 <= 3`; order-3 part `(-1)^j u D^j g D^(3-j) f` from `k = 4`) [W] CYC s.3; the full jet matrix
from `P_ab = Sum_(k=a+b+1..4) (-1)^(k-1-a) C(k-1-a,b) D^(k-1-a-b) A_k` [W] CF s.4.3; `det P = u^4` (anti-diagonal `u,-u,u,-u`, the longest
permutation of four letters is even) [W] CF s.4.3, CA s.3. (iv) 5. `Bmat = Wr^T P Wr` ((iii) on the basis) [W]. 6. Hand values: at `u = 0`,
for power series `B(f,g) = (2/5)(f Dg - g Df) + (f D^2 g - g D^2 f)`, `(f, Df, D^2 f)(0) = (1,0,0), (0,-1,-1), (0,0,2)` for `f_0, f_1, f_2`:
`b_01 = -7/5`, `b_02 = 2`, `b_12 = 0`; with `Y_k` the `l^0` part of `D^k y_3` (`Y_0 = h`, `Y_1 = -1 + (3/10)u + ..`, `Y_2 = u^-1 - 9/10 + ..`,
`Y_3 = u^-2 - (7/5)u^-1 + 87/100 + ..`): `b_13 = 0 + 1 - 1 - 1 = -1`, `b_23 = 0` [W] CF s.4.4. 7. `Pf = -b_02 b_13 = 2`, `det Bmat = 4 =
(det Wr)^2 u^4`, so `det Wr != 0` [W] CF s.4.4; CYC s.3. The entry `b_03 = -9/10` is [V] only (CF S3; CYC s.3) and is not needed (it
multiplies `b_12 = 0`). Optional: Abel's identity gives `D(u^2 det Wr) = 0`, so `det Wr = C u^-2` by S3 [W] CF s.4.4; `C = -2` [V] (`C^2 = 4`
follows from 7).

**S6** (C3). 1. `P^(-1) = Wr Bmat^(-1) Wr^T` in `Frac(R)` (from S5 (iv)) [W] CYC s.4. 2. With the exchange matrix J: `PJ` is upper
triangular with diagonal `P_(c,3-c)`, so is `(PJ)^(-1) = J P^(-1)` with diagonal `1/P_(c,3-c)`; hence `(P^(-1))_(ba) = 0` for `a + b < 3`
and `(P^(-1))_(0,3) = 1/u` [W] CYC s.4; CF s.5 (explicit inverse, `P P^(-1) = 1`).

**S7** (C4). 1. `l -> l + 1` defines a K-algebra automorphism `sigma` of the polynomial ring R, inverse `l -> l - 1` [W] CF2 s.6.
2. `D sigma = sigma D`: both are sigma-twisted derivations, agree on K and on `l` (`D(l+1) = D l = sigma(D l)`) [W] CF2 s.6 (CYC s.5 gives the
one-line reason "D l in K"). 3. `sigma f_i = f_i`, `sigma y_3 = y_3 + g`, so `sigma(S) = S` [W]. 4. `sigma P(f,f') = P(sigma f, sigma f')`
(coefficients in K); B is in Q, fixed by sigma; so `B(f, y_3) = B(f, y_3 + g)`, i.e. `B(f,g) = 0` [W] CYC s.5.

**S8** (C5). 1. `span(y_0, y_1)` isotropic (S7 and alternation), hence Lagrangian in the 4-dimensional symplectic space `(S, B)` [W] CYC
s.6. 2. Adapted basis: `x -> B(., x)|Lam` from S onto `Lam*` (kernel `Lam^perp = Lam`); choose `y_3'`, `y_2`; replace `y_2` by
`y_2 - B(y_2, y_3') y_0`; then `Bmat` is the signed permutation J and `Bmat^(-1) = -Bmat` [W] CF2 s.7 (CYC: "linear algebra gives").
3. `s_r := Sum (Bmat^-1)_ij phi_i D^r phi_j`, `tau_k := Sum (Bmat^-1)_ij phi_i th^k phi_j`, `s_0 = tau_0 = 0` [W] CYC s.6. 4. Leibniz on
`y_j = y_0 phi_j`: `T_1 = y_0^2 s_1`, `T_2 = 2 y_0 (Dy_0) s_1 + y_0^2 s_2`, `T_3 = 3 y_0 (D^2 y_0) s_1 + 3 y_0 (Dy_0) s_2 + y_0^2 s_3`; with
`D = Dt th`: `s_1 = Dt tau_1`, `s_2 = D(Dt) tau_1 + (Dt)^2 tau_2`, `s_3 = D^2(Dt) tau_1 + 3 Dt D(Dt) tau_2 + (Dt)^3 tau_3` [W] CYC s.6; CF2 s.7.
5. S6 (`T_1 = T_2 = 0`, `T_3 = 1/u`), `y_0`, `Dt` units of K: `tau_1 = tau_2 = 0`, `tau_3 = 1/(u y_0^2 (Dt)^3)` [W]. 6. In the adapted basis
`tau_k = -th^k phi_3 + phi_3 th^k phi_0 - t th^k phi_2 + phi_2 th^k phi_1`; for `k >= 2` `th^k 1 = th^k t = 0`; `tau_2 = -th^2 phi_3 - tY = 0`
gives (i); `tau_3 = -th^3 phi_3 - t th Y = (Y + t th Y) - t th Y = Y` gives (ii) [W] CYC s.6. 7. Where isotropy enters (only in 6) [W] CF2 s.7
(sentence to add, CF-E7).

**S9** (Theorem CYC). 1. `D (1/Dt) D = Dt th^2`, so `M_cyc = (Dt/y_0) th^2 o (u y_0^2 (Dt)^3) o th^2 o (1/y_0)` [W] CYC s.7; CF2 s.8.1.
2. For `y = Sum c_i y_i` (adapted basis) S8 gives `u y_0^2 (Dt)^3 th^2(y/y_0) = c_2 - c_3 t`, which `th^2` kills: `M_cyc y = 0` on S [W].
3. The `D^4` coefficient of `M_cyc` is `(1/y_0)(1/Dt)(u y_0^2 (Dt)^2)(1/Dt)(1/y_0) = u = A_4` [W]. 4. `E := M_cyc - PF4 = Sum_(k<=3) e_k D^k`
kills `(f_0, f_1, f_2, y_3)`, so `(e_0..e_3) Wr = 0`; `Wr` invertible (S5) and R a domain give `E = 0` [W] CYC s.7; CF2 s.8.1. 5. Scope
(RS1 A1): every step takes place in K or `K[l]`, where non-zero elements are units; S8 needs `y_0`, `Dt` invertible in K only; step 4 needs
only `det Wr != 0` [W] RS1 A1; for a general Lagrangian pair `(y_0, y_1)` of `S cap K` the proof is CYC's verbatim with `y_1` for `g`
(B has rank 2 on `S cap K` with radical `Q g`, so these planes contain g) [W*] CA2 s.4; RS2 s.6.3. 6. Scale invariance: the five factors
scale by `(1/mu)(mu/nu)(nu^2)(mu/nu)(1/mu) = 1` [W] CF2 s.8.1; CA s.2.5.

**S10** (Corollary). 1. The four definitions and the closing relation give `Y_0 = z^5 (1/c_4) D (1/c_3) D (1/c_2) D (1/c_1) D (1/c_0) D Y_0`
[W] CYC Cor. 2. `(1/c_4) D (1/c_3) D (1/c_2) D (1/c_1) D (1/c_0) D o (1/y_0) = (1/y_0) D M_cyc` (insert `c`) [W] CF2 s.8.2. 3. Hence the system
reads `z^5 D M_cyc Psi = Psi`; by S1, `PF5 Psi = 0` reads `z^5 D PF4 Psi = Psi`; apply S9 [W] CYC Cor. 4. "iff": the same identity read
in both directions [W] CF2 s.8.2. 5. The rings: the identities lie in `K[D]` (in `Q[[u]][D]` when `y_0(0) != 0`), so they hold in every
differential extension with z constant and the `c_a` invertible [W*] RS2 s.5 item 8 (one sentence to add). PB s.2.3 adds that (CYC) is
equivalent to `PF4 = M_cyc` because `c_4 = varpi_0` and `D o O = 0` forces `O = 0` [W*].

**S12** (hypotheses of U). 1. Units and the simple pole of `b_1`, `bt_1` [W] PB s.2.8. 2. `prod bt = -v'^5/v`; equality with `-1/(u(1-u)^5)`
iff `v'^5/v = 1/(u(1-u)^5)`; with `v = u + O(u^2)` this picks the fifth root `v' = v^(1/5) u^(-1/5)/(1-u)` and integrates to `v^(4/5) =
u^(4/5)/tau` [W] PS s.3; PB s.2.8; `Int_0^u t^(-1/5)(1-t)^(-1) dt = (5/4) u^(4/5) 2F1(1,4/5;9/5;u)` termwise [W*] RN s.3(a).
3. Unit shape iff `y_0(0) != 0` (examples `f_1`, `f_2`) [W] CA2 s.4. 4. Necessity [W*] UR s.3 Rmk (c).

**S13-S16** (Theorem U). 1. *Reduction* (UR s.3 (i)): `D_0 = diag(d_a)`, `d_(a+1)/d_a = bt_a/b_a`, consistent around the cycle by equal
products, unit entries, `A D_0 = D_0 At`; `P = D_0 Ph` gives `z Ph' = [At, Ph] - z Dl Ph`; order `z^0`: `[At, Ph_0] = 0`; order `z^(k+1)`:
`[At, Ph_(k+1)] = Ph_k' + Dl Ph_k` [W]. 2. `y = At^(-1)` holomorphic (`1/bt_1 = O(u)`), `y^5 = 1/B` [W] (ii). 3. "`At` is regular semisimple
over K, so its commutant in `gl_5(K)` is `K[y]` and `gl_5(K) = K[y] + im(ad At)` (direct)" [S] (no justification given). 4. `K[y] cap
gl_5(O) = O[y]` ("for each `j <= 4` some entry of `y^j` is a unit"; the disjoint supports of the `y^j` are implicit) [W*]. 5. The
projection formula `pr(X) = (1/5)[Sum_a xi_a/W_a^(m)] At^m` [S] (stated; the trace-form derivation is not written). 6. Pole bound: `B/W_a^(m)`
is a product of `5 - m` weights, at most one of them `bt_1` [W]. 7. S15 (UR Lemma 1): `Pi = Sum_a pi_a E_(a,a+d)`, `pi_a = rho_a W_a^(d)`;
`[At, Pi] = Sum_a (rho_(a+1) - rho_a) W_a^(d+1) E_(a,a+d+1)`; `rho_(a+1) - rho_a = sg_a := xi_a/W_a^(d+1)`, holomorphic, `sg_a(0) = 0` when the
window `{a..a+d}` contains 1; `Sum sg_a = 0` is `pr(X) = 0`; holomorphy of `pi_a` iff `rho_a` holomorphic and `rho_a(0) = 0` for
`a = 2-d..1`, which reduce to the single condition `rho_1(0) = 0`, fixed by the common additive function; two solutions differ by `O[y]`
[W] (iii). 8. S16 (UR Lemma 2): `pr((y^j)') = -(j/5)(log B)' y^j` ("every weight occurs in j of the five windows") [W*]; `pr(Dl y^j) =
(tr Dl/5) y^j` [S] (stated); the `y^j` component `u c_j' + (j/5 + u h_j) c_j = -u phi_j`, `h_j`, `u phi_j` in O by step 6; pivots `n + j/5`
non-zero for `j = 1..4`, regular equation for `j = 0` [W] (iv). 9. *Induction* (v): `Ph_0`: `Pi_0 = 0`, S16 with `Z = 0` gives `Ph_0 =
rho_0(u) 1`, `rho_0 = const exp(-Int tr Dl/5)`; given `Ph_k = Pi_k + c_k`, choose `c_k` by S16 with `Z = Pi_k' + Dl Pi_k`, then `X_(k+1) =
Ph_k' + Dl Ph_k` is holomorphic with `pr(X_(k+1)) = 0` and S15 gives `Pi_(k+1)` [W]. 10. *Uniqueness clause*: "the only freedom is one
constant per order, which is `P -> P c(z)`" [S] (that every regular gauge with `P_0` invertible arises from the construction, and that the
constants are exactly realised by `P c(z)`, is not written). 11. *Graded clause* ("the freedom is `c(z^5)`") and existence of a graded
gauge [S]: UR s.4 one sentence (covariance); GI s.1 "it can be taken graded". 12. *Field*: UR states U over `C[[u]]` (or convergent
series), GI quotes it over `Q[[u]]`; the rationality of the construction (including the constant in step 9) is not discussed [S].
Overall: PROVED in UR s.3; audited (NFS item 2) without recorded gaps; steps 3, 5, 8 (second fact), 10, 11, 12 need `\draftnote`s or
authors' justification (s.7). Evidence for existence and uniqueness: s.5 E2.

**S17** (graded gauges; `D_i(z) = D_0(zeta^-i z)`). 1. `S A S^(-1) = zeta A` for `S = diag(zeta^(-a))`, and the same for `At`, so
`S P(u, zeta z) S^(-1)` is again a regular gauge; averaging (or the uniqueness clause) gives a graded one [M] (UR s.4 says only "both
systems are covariant under `z -> zeta z` combined with `diag(zeta^-a)`, so P can be taken graded"). 2. Covariance of both column sets:
quintic `Q^(i)_(a,n) = zeta^(i(a+1-n)) Q^(0)_(a,n)` [W*] CS s.1 (from `P_a^(i)(z) = zeta^(ia) P_a^(0)(zeta^-i z)` and `Delta_i^(1/2)`); model:
S18 (c) form `zeta^(i(b+1)) F_b(zeta^-i t)` [W] UR s.2; EV s.2.2. 3. With P graded, `(P Q^m)^(i)_a(z) = zeta^(i(a+1)) (P Q^m)^(0)_a(zeta^-i z)`,
hence `D_i(z) = D_0(zeta^(-i) z)` [S] (UR s.4 one sentence; PS s.4.5 "confirmed"; the three-line computation is not written).

**S18** (model rows). (a) 1. (Mod1), (Mod2) from differentiation under the integral and the exact form `d(x^m e^(lambda/z))` [W] UR s.1;
the five rows (the `b = 1` row uses (Mod2) at `m = 0`) [W] PS s.1 (ii), GR s.2. (b) 2. Chain rule: `z d/dv[e^(U/z) sigma^(c-1) F(t)] =
e^(U/z) sigma^(c-2) T_(c-1) F` (`U' = -1/sigma`, `sigma' = 1/(5 sigma^4)`), giving the T-system [W] PS s.1 (iv) (UR s.2 only states it).
3. Uniqueness: `F_0` enters the closing relation at `t^(n+1)` with coefficient `-(1/5) Sum_i (k_i - 4n) = 4n` (`Sum k_i = 0`) [W] PS s.1 (v)
(UR: "pivot `-4n`"; the sign depends on the side; immaterial). 4. Homogeneous WKB columns differ by `d(z) = Sum d_k sigma^(4k) t^k`,
homogeneous only if `d_k = 0`, `k >= 1` [W] UR s.2. (c) 5. Given (H), the `u^0` terms of the hatted rows satisfy the limit system (S19) and
are homogeneous, hence equal the model rows [W] UR s.2 proof.

**S19** (a) The order count in the u-form of PF5 at `z = u^(4/5) zhat` [W] UR s.1. (b) [W*] UR s.2 (one paragraph; "not machine-checked
beyond the ODE").
