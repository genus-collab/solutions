# DIGEST_B_3 (key read_B), s.2 part 2: statements S20-S38 (R0; the vanishing solutions and (G); symplecticity; residues)

Continuation of `DIGEST_B_2.md`. Same format and abbreviations (`DIGEST_B.md`).

### (C) The gauge and the WKB columns

**S20** `thm:gau-R0` (thm; "R0"). Assume (U). Let `Q` be the quintic's WKB columns (QRR-normalised at `q = 0`; rows built by S10, so they
satisfy the closing relation in the conifold frame) and `Q^m(v(u))` the model's homogeneous WKB columns (S18) at `v = v(u)`. Then
`P_can Q^m(v(u))` is a formal WKB solution of the quintic system whose exponentials are those of the quintic (the choice of `v(u)`:
`-(5/4) zeta^i v^(4/5) = u_i - u_i^c`, D6-D7). By formal WKB uniqueness (the roots `zeta^i L` are distinct; the transport operator is first
order with kernel `Q L`; PS-C10, RS2 s.5 item 26): `Q D(z) = P_can Q^m(v(u))` with `D = diag(D_i(z))` constant, `D_i` in `Q[[z]]`, the SAME
exponential `e^(u_i/z)` stripped on both sides (PS s.4.1; UR's `e^((u_i - u_i^c)/z)` differs by `e^(u_i^c/z)`, which is not a power series in
z). Covariance (S17): `D_i(z) = D_0(zeta^(-i) z)`. Hatted: `Qhat(u,zhat) D = Ghat(u,zhat) R^(tau zhat)`, the right side a power series in u
whose coefficients are power series in zhat. Colour 0, row 0, hatted: `Phi^QRR(w,u) D_0(s^4 w) = Phi^can(w,u)` with
`Phi^can := C_0 x (colour-0 entry of row 0 of the hatted P_can Q^m(v(u)))` (the factor `1/C_0` is common; PB s.3). Ring of comparison:
`Q[L,1/L][[z]]` inside `Q((s))[[w]]` (PS s.4.2, C8). Source: UR s.4; EV s.2.3 steps 1-5; PB s.3; PS s.4.1-s.4.2, s.4.5.
Status: PROVED given (U) (with (CYC) for the actual rows, S10). The value `D = 1` is the connection theorem (ST5 S1-S6; PROVED modulo (U)
and (QRR); other digests). Superseded labels: UR s.4 "(H) Purity" (CONJECTURAL) and "Canonical normalisation (CONJECTURAL)" are now
consequences of `D = 1`; the note must not present them as conjectures. Deps: S10, S13, S17, S18.

### (D) Residue relation, vanishing solutions, (G)

**S21** `lem:res-R1` (lem; "(R1)"). For every regular gauge P between cyclic systems with one simple pole, `A = u^(-1) r E_12 + ..`,
`At = u^(-1) rt E_12 + ..` (`r, rt != 0`): `r E_12 P(0) = rt P(0) E_12`, i.e. row 2 of `P(0,z)` is `P_22(0,z) e_2^T`, column 1 is
`P_11(0,z) e_1`, and `r P_22(0,z) = rt P_11(0,z)`; here `r = rt = -1`, so `P_22(0,z) = P_11(0,z)`. Reading: `P(0,z)` commutes with the residue;
at the conifold point the `H^2`-row of the quintic is proportional to the `H^2`-row of the model, and the `H`-row of the model enters only
the `H`-row of the quintic. Source: UR s.5; GI s.3. Status: PROVED (order `u^(-1)` of the gauge equation). Deps: D11.

**S22** `lem:van-model` (lem; GI "Lemma V"). `I^C_m` (D13) is the formal residue `Res_(x=0) x^m e^(-(x^4/4 + v/x)/z) dx`; (i) `z d/dv I^C_m =
-I^C_(m-1)`; (ii) `m z I^C_(m-1) - I^C_(m+3) + v I^C_(m-2) = 0` (residue of the exact form `d(x^m e^(lambda/z))`); hence `Y^C` solves the model
system (the `b = 1` step uses (ii) at `m = 0`: `I_3 = v I_(-2)`), and `Y^C(0,z) = -e_1` (only `I^C_(-1)` has a `v^0` term).
Source: GI s.2; GR s.2. Status: PROVED. Deps: D8, D13.

**S23** `lem:van-shape` (lem). `Y^C_b` is a sum of monomials `v^j z^(-n)` with `n = 5i + c_b + 1`, `5j = 4n + c_b + 1`; at order `u^N` of
`Y^C(v(u),z)` only `n <= (5N - c_b - 1)/4` occur, so `Y^C(v(u),z)` lies in `Q[z^(-1)][[u]]^5`; leading terms `Y^C = (-v/z, -1, -v^4/(24z^4),
-v^3/(6z^3), v^2/(2z^2)) + ..`, `[z^(-1)] Y^C_0 = -v`; homogeneity `Y^C_b(v,z) = sigma^(c_b+1) y_b(z/sigma^4)`, `y_b(t) := Y^C_b(1,t) = eps_b
Sum_i (-1/(4t))^i/i! (-1/t)^(c_b+4i+1)/(c_b+4i+1)!`. Source: GI s.2; GH s.1-s.2; TH s.1.3. Status: PROVED. Deps: S22.

**S24** `lem:van-quintic` (lem; GI "Lemma Q"). For the quintic system (and verbatim for any cyclic one-pole system, with `T_A :=
Int_0^u b_0` in place of `t_c`): (a) `y^(n)_a = 0` unless `n = 1 - a mod 5`, and then `ord_u y^(n)_a = (4n + c_a + 1)/5` exactly, with the same
leading coefficient as `[z^(-n)] Y^C_a(v(u))`; (b) every integrand is a power series (the pole of `b_1` multiplies `y^(n)_2`, non-zero only
for `n = 4 mod 5`, of order `(4n+4)/5 >= 4`; GI writes `>= 1`), so no logarithm arises; (c) `Y^C_Q` lies in `Q[z^(-1)][[u]]^5` (order
`>= 4n/5` at `z^(-n)`), solves the quintic system, `Y^C_Q(0) = -e_1`; (d) `[z^(-1)] Y^C_Q = -t_c e_0`; also `[u^0](Y^C_Q)_0 = 0` (stated by
PS-C7 as used in P(b) step 7). The normalisation `t_c(0) = 0`, `t_c'(0) = 1` refers to the programme family (CA-K2; RS2 s.5 item 28).
Source: GI s.2; GR s.2; PS s.3; TH s.1.3. Status: PROVED. Deps: D5, D9, D13, S23, S25 (for the leading coefficients).

**S25** `lem:van-uniq` (lem; GI "Lemma 1"; rename: clashes with UR Lemma 1). Let `B = u^(-1) B_(-1) + Sum_(j>=0) u^j B_j` with constant
(z-free) matrices `B_j`. If `Y` in `L[[u]]^5` satisfies `z dY/du = B Y` and `Y(0) = 0`, then `Y = 0`. For `B_(-1) = r E_12` (nilpotent)
`det((n+1)z - B_(-1)) = ((n+1)z)^5`, a unit of `L = Q((z))`. Source: GI s.2; GR s.2. Status: PROVED. Deps: none.

**S26** `rmk:van-monodromy` (rmk). Over L the regular solutions are parametrised by `Y(0)` with `Y(0)_2 = 0`; the fifth solution is
`log(u) Y^C_Q + W`, W regular, `W_2(0) = z`; so `Y^C_Q` spans the image of the logarithm of the monodromy (the vanishing cycle), and the
same for `Y^C`. Not used by any proof; it explains why the initial value `-e_1` (and not "pure 1/z") characterises the vanishing solution.
Source: GI s.2 Remark; GR s.2 (W_2(0) = z re-derived; the characterisation). Status: remark (re-derived by GR).

**S27** `lem:van-diag` (lem; "the diagonal of the gauge"). For a graded regular gauge: order `z^0` gives `p_(0,a+1)/p_(0,a) = bt_a/b_a`;
order `z^1`, diagonal: `p_(0,a)' = b_a p_(1,a) - bt_(a-1) p_(1,a-1)`; with the weights `1/p_(0,a)` the right sides cancel, so `prod_a p_(0,a)` is
constant, `= 1` for `P_can`; with equal products (`v'^5 = v/(u(1-u)^5)`) this gives `p_(0,0)^5 = prod_(j<4)(b_j/bt_j)^(4-j) = v/(u varpi_0^5)`,
hence `p_(0,0)(u) = tau^(-1/4)/varpi_0` and `Ghat(u,0)_00 = 1/varpi_0` (`=: 1/xhat`). Also `g_0 g_3 = 1` (GR s.2; used in GH s.2).
Source: GI s.5 "The diagonal"; GR s.2. Status: PROVED. Used by P(b) step 5 (`ph^can_0 = 1`, PB s.2.8) and by GH (T). Deps: S13, S17, D7.

**S28** `thm:van-G` (thm; the gauge identity "(G)"; do NOT call it "Theorem G": clash with the Gap Theorem, S01 C9). For every regular
gauge P from the pulled-back model `At` to the quintic `A`: `P(u,z) Y^C(v(u),z) = P_11(0,z) Y^C_Q(u,z)` in `L[[u]]^5`.
Corollaries: (a) for graded P, `P_11(0,z)` is a unit of `Q[[z^5]]`, `P_can = P P_11(0,z)^(-1)` is a graded regular gauge, and it is the
unique regular gauge with `P Y^C = Y^C_Q` (TH s.1.4 property (A): `P Y^C(v(u))` has no positive power of z); (b) for any regular gauge,
`P Y^C` has no positive power of z iff `P_11(0,z)` is constant (the `u^0` term of `P Y^C` is `-P_11(0,z) e_1`); (c) universality: for any
pair of cyclic systems with one simple pole and equal products, the regular gauge maps the vanishing solution of the second (iterated
integrals, value `-e_1`) to `P_11(0,z)` times that of the first, whose `[z^(-1)]` term is `-Int_0^u b_0 e_0`.
Source: GI s.3; GR s.2 ("the whole content is (R1) + Lemma 1"); PB s.2.8 (the four ingredients); S01 E1. Status: PROVED given Theorem U
(only for the existence of a regular gauge); AI-refereed (GR); human check pending (MS5 s.5 item 1.6). Deps: S13, S21, S22-S25.

**S29** `prop:sym-main` (prop; "symplecticity of the regular gauge"; GI calls it "Proposition S": rename, clash with UR Prop. S).
For every regular gauge between two self-dual cyclic one-pole systems with equal products and `r = rt`:
`P(u,-z)^T eta P(u,z) = c(z) c(-z) eta`, `c(z) := P_11(0,z)`. Consequences: `P_can(-z)^T eta P_can(z) = eta`; the symplectic regular gauges are
`P_can c(z)` with `c(z) c(-z) = 1`, for graded ones `c = +-exp(a_1 z^5 + a_3 z^15 + a_5 z^25 + ..)`; symplecticity fixes only the even part
of `log c` and cannot see the odd multiples of `z^5` (e.g. `c_5`): only (G) singles out `P_can`. Source: GI s.4; GR s.2; S01 E2; ST2 s.2
step 3; TH s.1.4 (Sp). Status: PROVED given Theorem U (its uniqueness clause for the pair `(At, At)`). Deps: S13, S21, D16.

**S30** `cor:sym-Gstar` (cor; "(G*)"). For every regular gauge `Y^C_Q(-z)^T eta P(z) = c(z) Y^C(-z)^T eta`; for `P_can`,
`Y^C_Q(-z)^T eta P_can(z) = Y^C(-z)^T eta`. Source: GI s.4; GR s.2 (negative controls supplied, GR s.4 item 3). Status: PROVED. Deps: S28, S29.

**S31** `prop:res-allorder` (prop; "residue relations to all orders"). Write `P_can = Sum_n u^n P^(n)(z)`, `Y^C(v(u),z) = Sum_n u^n
Ychi_n(z)` (GI writes `eta_n`: rename, clash with the pairing), `Y^C_Q = Sum_n u^n theta_n(z)`, `Ychi_0 = theta_0 = -e_1`. Then
`(G)_n: P^(n) e_1 = Sum_(j<n) P^(j) Ychi_(n-j) - theta_n` and `(G*)_n: 5 e_2^T P^(n) = Sum_(j<n) theta_(n-j)(-z)^T eta P^(j) - Ychi_n(-z)^T eta`:
the `H`-column and the `H^2`-row of every `[u^n] P` are fixed by lower orders. `n = 0`: (R1) canonically (`P(0,z) e_1 = e_1`,
`e_2^T P(0,z) = e_2^T`). `n = 1`: `P^(1)_(a1) = (delta_(a0) - P^(0)_(a0))/z`, `e_2^T P^(1) = (e_3^T P^(0) - e_3^T)/z`; so `P^(1)_21 = 0`
((R2)), `P^(1)_31 = -pi/z` ((R3b)), `P^(1)_20 = pi/z` ((R3a)), `pi(z) := P_30(0,z)`. `n = 2`: `P^(2)_21 = -pi/z^2` ((R3c)). `n = 3, 4`: explicit
row-2 relations (GI s.5), the `n = 4` one containing `(1 - P_22(0,z))/(24 z^4)`, the first place where the normalisation enters. Data in
frame (0,0): `Ychi_1 = theta_1 = -e_0/z`, `Ychi_2 = -(5/9) e_0/z + e_4/(2z^2)`, `theta_2 = -(7/10) e_0/z + e_4/(2z^2)`; slot 2 starts at `u^4`
with `-(1/24) z^(-4)`. Source: GI s.5; GR s.2. Status: PROVED given Theorem U. Not needed by the chain (structure). Deps: S28, S30.

**S32** `lem:res-R23` (lem; UR "(R2), (R3)"). For every regular gauge (here `b_2(0) = bt_0(0) = 1`): (R2) `P_21 = O(u^2)` at every order in z
(order `u^0`, entry (2,1): `b_2(0) P_31(0) - bt_0(0) P_20(0) = 0` by (R1)); (R3) with `pi(z) := P_30(0,z)`: `P_20 = u pi/z + O(u^2)`,
`P_31 = -u pi/z + O(u^2)`, `P_21 = -u^2 pi/z^2 + O(u^3)`: the block rows (2,3) x columns (0,1) is `pi(z)(u/z, 1)^T (1, -u/z)` to leading order,
the bivector `pi(z)(1 + (u/z)H)(x)(1 - (u/z)H) = -pi(z) z^(-2) e(z)(x)e(-z)`, `e(z) = z 1 + uH` (D20). Consequence: no `H(x)H` in `G_0`, `G_1`
(forced). Hatted: the `H`-part of `G_2` is `alpha zhat (H - zhat 1)(x)(H + zhat 1)`, `alpha = -127/15750` (frame (0,0)); the sign differs
from the flat form because the Fourier vector of H is minus the flat class (`N_1 = -1`), NOT by an ordering convention (NFS s.2 item 8;
dictionary in VT `v3_fourier`). Canonical values: `pi(z) = (127/3150) z^3 + (1416472237/725359556250000) z^8 - ..`; frame `(1/3,-2/7)`:
`pi = -(239/450) z^3 + ..`; UR's `z^8` value `-77141089/362679778125000` belongs to the stored normalisation (difference `c_5 (127/3150)`;
GI s.5, GR rf_misc). Frame dependence: the `H(x)H` coefficient of `[u^2 zhat] G` is `-127/15750 - (2/5) beta` (NFS item 13).
Source: UR s.5; GI s.5; NFS s.1 items 11, 13, s.2 item 8. Status: PROVED given U (NFS: PROVED-mod(Theorem U)); values VERIFIED.
The reading "two Euler-dilaton vectors" is a remark. Deps: S21, S31.

**S33** `thm:res-deg` (thm; UR "Theorem D"; rename: clashes with `D`, `D(z)`, CP's Theorem D). Given U and `D = 1` (UR: and (H)):
`[u^n] G` is a polynomial in zhat of degree at most `floor((5n+2)/4) = n + floor((n+2)/4)` (`n = 4q + j`: `5q + (0,1,3,4)_j`); the a priori
grading bound is `floor((5n+4)/4)`; the improvement uses (R1), (R2) (`eps >= 2` for `(a,b) = (2,1)`, `eps >= 1` if `a = 2` or `b = 1`,
`5 !| k`). Equivalent "rule W" form: `G - 1` is a formal power series in z with coefficients regular at `u = 0` in the frame
`(1, uH, H^2, H^3, H^4)` (RN s.2; NFS item 10: as a bound equivalent to U + (R1) + (R2) + the normalisation). Equality: VERIFIED only, and
attained only in generic frames (`beta = -127/6300` gives `deg G_2 = 2`; NFS s.2 item 3). Source: UR s.5; RN s.2; NFS s.1 items 10, 13,
s.4. Status: bound PROVED; equality VERIFIED (UR's "Theorem D equality PROVED" is an overclaim, NFS s.4). Use: finiteness of graph terms
needs only `mu(b,k) >= (4k-4)/5` (GH s.1); Theorem D refines it. Deps: S20, S21, S32.

**S34** `prop:res-orderu` (prop; "the order-u propagator"). With `bt_a/b_a = 1 + eps_a u + ..`: `eps = (e, -2e-2f, e, f, f)`, `e = 2v_2 - 2
gamma_t`, `f = 2v_2 - 1 - w_1` (`v = u + v_2 u^2`, `v_2 = 5/9`); `p_(1,b)(0) = (-e-f, 0, -e-f, e+2f, e+2f)` = `(8/45, 0, 8/45, -1/15, -1/15)` in
frame (0,0) [`(-7/45, 0, -7/45, -1/15, -1/15)` in `(1/3,-2/7)`, VT v2_run]; `G_1 = zhat X_1`, `Xhat(0,1) = Xhat(1,0) = (e+2f)/5 = -(2 gamma_g -
4/3)/5 = -1/75` (frame independent), `Xhat(2,4) = Xhat(4,2) = (e+f)/5 = -8/225 + w_1/5` (string direction, frame dependent); no `H(x)H`.
The leaf coefficient `-p_(1,4)(0) = 2 gamma_t + 2 w_1 - 4/3 = 1/15` in EVERY frame (GR s.4 item 1; S01 C6). Source: UR s.5, s.6.2; NFS
item 12; `ur_chain_check.log` (X1). Status: UR "PROVED formula, VERIFIED"; the derivation is not written in UR (NFS item 12: "PROVED by
hand from the G-equation mod the hatted R-equation", not located in the sources read). Not used by the chain.

**S35** `cor:gau-F` (cor; UR "Corollary F", with RN's companion facts). (a) Given `D = 1` (UR: given (H)): `R~(zhat;u) = G(u,zhat)
M_0(tau(u) zhat)`, `G = Psi(u)^(-1) Ghat Psi(0)`, `Ghat` the hatted `P_can`; `[u^n] G` polynomial in zhat; `M_0` on the RIGHT; the same
`M_0(tau zhat)` in every conifold frame, frame-dependent G. (b) Uniqueness of the rescaling (RN s.3(a)): if `R~ M_0(c zhat)^(-1)` and
`R~ M_0(c' zhat)^(-1)` were both polynomial at every u-order then `c = c'`; and `tau` is the canonical coordinate,
`u_i/z = -(5/4) zeta^i/(tau(u) zhat)` (RN: PROVED; NFS item 3: PROVED-mod(hatted R-equation to all orders, i.e. PF + (H))). (c) `M_0` on the
left fails at every order (RN s.3(c), VERIFIED). (d) `G(zhat)^T G(-zhat) = 1` since `R~` and `M_0(tau zhat)` are symplectic (RN s.3(d), PROVED).
Source: UR s.4; RN s.3; NFS s.1 items 1, 3. Place after the connection theorem. Deps: S20, S28, S29, D = 1.

**S36** `prop:gau-actual` (prop; "the actual gauge"). (a) Given (H): the actual R-matrix defines `P^act = P_can c^act(z^5)`,
`c^act(z) = P^act_11(0,z)`; `G^act`, hence `P^act`, is symplectic, so `c^act(z) c^act(-z) = 1`, `c^act = exp(a_1 z^5 + a_3 z^15 + ..)`: the
`z^10`, `z^20` coefficients are forced; "the actual R-matrix is in the canonical normalisation" means `a_1 = a_3 = .. = 0`. (b) Reformulation
(ST2b s.6.1, from S28 (b), given (H)): (N) holds iff `P^act Y^C(v(u)) = Y^C_Q`. (c) Superseded by `D = 1` (MS5 s.2): `P^act = P_can`.
Source: GI s.6; GR s.2, s.4 item 2; ST2b s.6.1. Status: (a), (b) PROVED given (H); windows VERIFIED (s.5 E7). Deps: S28, S29, S20.

**S37** `lem:hat-transport` (lem; GH). (i) hatted (G): `Ghat(u,zhat) Yhat^mod = Yhat^Q` (the same terms regrouped), `Yhat^Q = u x (series in u
with coefficients in Q[[zhat^-1]])`: no positive power of zhat; (ii) closed form of `phi_(i,m)` (D21), one Fourier mode `r = 2 - m`,
`phi(tau_m(Hhat^c)) = (-1)^m Ghat(u,0)_cc [zhat^(-m)] Yhat^mod_c`; (iii) Lemma L (uses only symplecticity of the rows of R'):
`Sum_(i,m) phi_(i,m) E'(m,n)_ij = [psi'^(n+1)] Sum_ab h_a(psi') eta^ab Q'^(j)_b(-psi')`, `h(x) := [Q'(x) f(x)]_(>0)`, for any covector;
(iv) (E): for `P_can`, `Sum_(i,m) phi_(i,m) E'(m,n)_ij = 0` for all `n, j`, all u-orders, every graph of every genus; (v) (T):
`Pi(T') := Sum_(i,d>=2) phi_(i,d) T'^(i)_d = vhat - t_c`, `vhat = u tau^(-3/2)/xhat`, leading `(p_(1,4)(0)/2) u^2 = -u^2/30` in every frame.
Corrections: "frame dependence" of the leaf bookkeeping withdrawn (GR s.4 item 1; S01 C6); dimension caveat concerns the vertex genus h,
not g (GR item 7); "needs the non-linear genus-0 extension" superseded (S01 s.2). Source: GH s.0-s.5; GR s.2, s.4. Status: (i)-(v)
PROVED given U for `P_can`. The transport lemmas of record for the chain are TH s.3.1 (Lemmas E, T, F; other digest).
Deps: S23, S27, S28, S29.

**S38** `rmk:ur-other` (remarks; for DIGEST_A coordination and for remarks). (a) Conjecture N as recorded in UR s.6.5 (hatted graph sum,
`M = M_0(tau zhat)`, `R' = Psi(u) G`; every stable graph with an edge `O(u^(2g-2))`; vertex term `c_g t_c^(2-2g)` + regular,
`c_g = (-1)^(g-1) B_2g/(2g(2g-2)) 5^(1-g)`; flat form `Omega^quintic = P.(Omega^LG at v(u))`); UR's CONJECTURAL is superseded (MS5 s.2).
(b) Counting (UR s.6.1): a genus-g graph with E edges is `O(u^(E+2-2g))`; only `G_1..G_(2g-3)` matter. (c) UR Prop. E (Euler-dilaton
first-order response) is PROVED only for genus 2, 3 (NFS s.4); not needed. (d) The Euler-dilaton element is not an exact symmetry
(genus-1 anomaly `-1/8`; VERIFIED p89). (e) UR s.8 exact statements 1-2 (PROVED, elementary; D6 and the function `W(x;u)`); the rest of
UR s.8 is SKETCH, and the Stage-4 phase is a different phase with the same critical data (ST4b K6). (f) UR s.7: the analogy with the
Hodge trick. (g) The single-leaf bookkeeping `-(g-1) c_g/15`, observed in UR s.6.2, is derived later (S01 E7; other digest).
