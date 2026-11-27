# Outline, part 4: s03 the two cyclic systems, (CYC), Theorem U, the gauge identity, symplecticity, R0 and the colours

## s03_gauge — target 10 pp

Purpose: prove that the quintic's rows form a cyclic one-pole system (CYC), relate it to the pulled-back model by the regular gauge of
Theorem U, single out the canonical gauge by the gauge identity (G), prove its symplecticity, define the normalisation series `D` by
comparing WKB columns in colour 0 (R0), put the rows into the differential ring, and transfer R0 to every colour, so that (F_g) is
defined before s07.

Covers: DB s.1 D1-D21, s.2 S3-S17, S20-S32, s.3, s.7; DD1 lem:pol-HPF, lem:pol-poly, lem:glu-wkb, thm:glu-main (a)-(d), cor:glu-colours;
CYC s.0-s.9; CF, CF2, CF3; CA, CA2; RS1 A1, A2, A4; RS2 s.5 items 1-16, 24-28, 44; UR s.1-s.5; GI s.1-s.5; GR s.2, s.4; PS s.3, s.4.1-s.4.5;
PB s.2.3, s.2.8, s.3; POLY s.1-s.2; PR s.2; EVID s.2.1, s.2.3; ER s.2, s.4; G2 (all pieces) with G2R, G2R_2, G2R_3; G1_3 s.2.10 (c);
GS s.2.5, s.2.9. Placement: the model rows are in s02 (`prop:set-modrows`); HPF and the polynomiality of power-series solutions are here,
because R0 needs them. Local notation of `sec:gau-U` (declared at the top of s03; `NOTATION_2.md` s.9): `\mathcal B^{(m)}_a`, `\mathcal Y`,
`dlog Theta_U`, `\mathfrak d`, `pr`, `chi^U_j`, `P^\natural`, `P^{[delta]}`, `Rot_i`, `zeta^(gr)`, `\Gau(A, A~)`, the field `\Bbbk`.
Subsections: `sec:gau-cyc`, `sec:gau-systems`, `sec:gau-U`, `sec:gau-van`, `sec:gau-symp`, `sec:gau-R0`.

### sec:gau-cyc — the cyclic form of PF4, (CYC)
* `def:gau-conc` (defn; DB D3, D17): `KK = Q((u))` with `D`; `KK[ell]`, `D ell = 1 - 1/u`; `y_3 = varpi_1 ell + rho_3`; solution space
  `SS = span(f_0, f_1, f_2, y_3)` (no claim that it is all of `KK[ell]`'s solutions, CF2 s.6); formal adjoint; concomitant `CC(f,g)`;
  jet matrix `JJ_ab`; Gram matrix `CCmat`, Wronskian `\Wr`; the automorphism `sh: ell -> ell + 1`; `\partial_t := (Dt)^-1 D`; sums `TT_k`;
  `M_cyc(y_0, y_1) := (1/y_0) D (1/Dt) D (u y_0^2 (Dt)^2) D (1/Dt) D (1/y_0)`, `t = y_1/y_0`. Renamings (RS2 item 6; CF3 E6): CYC's `tau`
  -> `sh`, `M` -> `M_cyc`, `P` -> `JJ`, `f_i = y_i/y_0` -> `phi_i`, `tau_k` -> `TT~_k`, `B` -> `CC`.
* `lem:gau-const` (lem; DB S3): `{f in KK[ell] : Df = 0} = Q`. Proof DB s.3 S3 steps 1-3 [W] (CYC s.1; CF s.2.1 makes the division by the
  unit `1 - u` explicit).
* `lem:gau-adj` (lem; DB S4): `PF4* = PF4`. Proof S4 steps 1-4 [W] (CYC s.2); the two-condition criterion as an alternative (CF s.3).
* `lem:gau-conc` (lem; DB S5): (i) Lagrange identity `g L f - f L* g = D CC(f,g)`; (ii) on `SS`, `CC` is `Q`-valued, bilinear, alternating;
  (iii) the jet matrix, `JJ_ab = 0` for `a + b >= 4`, `JJ_(3-j,j) = (-1)^j u`, `det JJ = u^4`; (iv) `CCmat = Wr^T JJ Wr`; `Pf(CCmat) = 2`,
  `det CCmat = 4`; `CC` non-degenerate; `det Wr != 0`. Proof S5 steps 1-7 [W] (CYC s.3; CF s.4.1-s.4.4, hand values `CC(f_0,f_1) = -7/5`,
  `CC(f_0,f_2) = 2`, `CC(f_1,f_2) = 0`, `CC(f_1,y_3) = -1`, `CC(f_2,y_3) = 0`). Corrections: RS2 item 2 (hand values replace "computed"). The
  entry `CC(f_0, y_3)` is computed only and multiplies `0` in the Pfaffian (write `*`, say it is not needed). Optional: `D(u^2 det Wr) = 0`,
  so `det Wr = C u^-2`, `C^2 = 4` (Abel; `lem:gau-const`); the value `C = -2` only in appendix C.
* `lem:gau-transv` (lem; DB S6, "transversality"): `Sum (CCmat^-1)_ij D^a y_i D^b y_j = (JJ^-1)_ab`, zero for `a + b <= 2`, `(JJ^-1)_(0,3) =
  1/u`; `TT_0 = TT_1 = TT_2 = 0`, `TT_3 = 1/u`, in `Frac(KK[ell])`. Proof S6 steps 1-2 [W] (CYC s.4; CF s.5 explicit inverse).
* `lem:gau-iso` (lem; DB S7, formal Picard-Lefschetz): `CC(f, varpi_1) = 0` for every power-series solution `f`; `CC(y_3, varpi_1) = 1`; every
  `span(y_0, varpi_1)` with `y_0` a power series is isotropic. Proof S7 steps 1-4 [W] (CF2 s.6; CYC s.5).
* `lem:gau-confexp` (lem): the conifold expansion `\Cexp_(alpha,beta)` of `def:set-confexp` is a well-defined ring map commuting with `D`
  (it respects PF4 and the concomitant relation `CC(C_0, I_1) = 0`, whose image `CC(varpi_0, varpi_1) = 0` holds by `lem:gau-iso`). Source
  FQ2 s.1 (statement); proof RS2 s.5 item 44, CA-K5 [W*], G5_2 s.2.2 [G]; expand in two sentences.
* `lem:gau-sg` (lem; DB S8, special geometry): for `y_0` a power series not proportional to `varpi_1`, `t = varpi_1/y_0`: adapted basis
  `(y_0, varpi_1, y_2, y_3')`, `phi_i = y_i/y_0`, `Yuk := \partial_t^2 phi_2`; (i) `\partial_t^2 phi_3 = -t Yuk`, so `\partial_t^2(y/y_0)` lies in
  `Yuk span(1, t)`; (ii) `Yuk = 1/(u y_0^2 (Dt)^3)`. Proof S8 steps 1-7 [W] (CYC s.6; CF2 s.7). Corrections RS2 item 7 (CF-E7: `TT~_1 = TT~_2 =
  0`, `TT~_3 = 1/(u y_0^2 (Dt)^3)` for every unit `y_0`; isotropy enters only through the adapted basis) and item 9 (CA-K6: `Yuk` normalised by
  the Darboux form; the Yukawa coupling of `def:set-genus0` is `5 Yuk`; the row coefficient `kappa_1 = Yuk Dt` has no 5).
* `thm:gau-cyc` (thm; DB S9, "(CYC)"): for `y_0` a non-zero power-series solution not proportional to `varpi_1` (any scale of `varpi_1`),
  `PF4 = M_cyc(y_0, varpi_1)` as operators over `KK`; more generally for every pair spanning a `CC`-Lagrangian plane of `SS cap KK` (these
  planes contain `varpi_1`). Scope: `y_0(0) != 0` is NOT needed (RS1 A1); invariance under `y_0 -> mu y_0`, `varpi_1 -> lambda varpi_1`, `y_1 ->
  y_1 + a y_0`. In the conifold frames: `PF4 = M_cyc(varpi_0, varpi_1)`. Proof S9 steps 1-6 [W] (CYC s.7; CF2 s.8.1); step 5 scope [W] RS1
  A1, general Lagrangian pairs [W*] CA2 s.4.
* `cor:gau-cycsys` (cor; DB S10): with `kappa =
  (Dt, L^5/(y_0^2 (Dt)^2), Dt, y_0, y_0)`, `L^5 = 1/u`, in any differential extension of `(KK, D)` with `z` constant and the `kappa_a`
  invertible, `zD Y_a = kappa_a Y_(a+1)` (`a` in Z/5) iff `y_0 Y_0` solves PF5 and `Y_(a+1) = zD Y_a/kappa_a` (`a = 0..3`); rings used later:
  `Q(zeta)((s))((z))`, `Q((z))((u))`, hatted series. Proof S10 steps 1-5 [W] (CYC Cor.; CF2 s.8.2), rings [W*] RS2 item 8 (one sentence
  to add). Both directions are used (DB S10). Log-free WKB rows are not part of it (RS1 A4): they are `thm:gau-R0` (iii).
* `rmk:gau-cyc`: dropped (D8).

### sec:gau-systems — the quintic system and the pulled-back model
* `def:gau-onepole` (defn; DB D10, D11; G2 s.2): `z d_u Y = A Y`, `A = Sum_a b_a E_(a,a+1)`, `b_a` units of `\Bbbk[[u]]` for `a != 1`, `u b_1`
  a unit (the one-pole shape; G2's "(H)", never printed as (H)); residue `r = (u b_1)(0)`; self-dual `b_0 = b_2`, `b_3 = b_4`; (EP) equal
  products `Bprod := prod_a b_a = prod_a b~_a`; the gauge equation `z P' = A P - P A~`; its solution module `\Gau(A, A~)` in `gl_5(\Bbbk[[u]])[[z]]`
  (no condition on `P_0`; a `\Bbbk[[z]]`-module); a regular gauge is an element with `P_0` invertible; graded `P_k = Sum_b p_(k,b) E_(b+k,b)`.
* `prop:gau-quintic` (prop; DB D5, S12; REV R1 item 4): in every conifold frame the rows satisfy `zD Y_a = kappa_a Y_(a+1)` with `kappa =
  (D t_c, L^5/(varpi_0^2 (D t_c)^2), D t_c, varpi_0, varpi_0)`; `b_a = -kappa_a/(1-u) = (t_c', -1/(u varpi_0^2 (1-u)^3 t_c'^2), t_c',
  -varpi_0/(1-u), -varpi_0/(1-u))` is a self-dual one-pole system with `r = -1`, `prod b_a = -1/(u(1-u)^5)`; the unit shape holds iff
  `y_0(0) != 0` (CA2 s.4). Proof GS s.2.9 [G]: `\Cexp` is a differential homomorphism (`lem:gau-confexp`) with `\Cexp(C_0) = varpi_0`,
  `\Cexp(L) = 1/s`, `\Cexp(C_1) = \Cexp(\mathsf W)/varpi_0^2 = D t_c` (`\Cexp(\mathsf W) = varpi_0 D varpi_1 - varpi_1 D varpi_0`, G5_2 s.2.2);
  applied to the forward relations of `prop:set-cyclicLR` it gives the four relations with `kappa = \Cexp(C_1 kappa^LR)`; the closing
  relation is `cor:gau-cycsys` (or `\Cexp` of `cor:set-cyclicPF5` (b)). Also UR s.1; PB s.2.8; CA s.2.1-s.2.2; RS2 items 27-28 (`t_c'(0) =
  varpi_0(0) = 1` in the conifold frames only).
* `prop:gau-equalprod` (prop; S12): the model pulled back by `v(u)`: `b~ = (v', -v'/v, v', -v', -v')`, `prod b~ = -v'^5/v`; equal products iff
  `v'^5/v = 1/(u(1-u)^5)` iff, for `v = u + O(u^2)`, `v = u tau(u)^(-5/4)`; necessity of (EP): a solution of the gauge equation with `P_0 != 0`
  forces equal products (order `z^0`: `A P_0 = P_0 A~`, hence `A^5 P_0 = P_0 A~^5`; G2_3 s.5.1 (e)). Proof DB s.3 S12 steps 1-4 [W]/[W*] (PS
  s.3; PB s.2.8; RN s.3(a)). Correction: unique given `v = u + O(u^2)`, not `v(0) = 0` (NFS s.4; PS s.3).

### sec:gau-U — Theorem U
* `thm:gau-U` (thm; S13; DEC D7, D15): statement G2_3 s.5.1 (a)-(d) [G], over any field `\Bbbk` of characteristic 0: (a) there is
  `P^\natural` in `\Gau(A, A~)` with `P^\natural_0 = \mathfrak d Theta_U`, diagonal with unit entries; (b) for any such `P^\natural` (G2R_2 s.3 R2)
  the map `c -> P^\natural c(z)` is a bijection `\Bbbk[[z]] -> \Gau(A, A~)`: every solution with coefficients in `gl_5(\Bbbk[[u]])` is
  `P^\natural c(z)`, its `z^0`-term is `c(0) \mathfrak d Theta_U`, it is regular iff `c(0) != 0`, and the regular gauges are `P c(z)`, `c` in
  `\Bbbk[[z]]^x`; (c) graded regular gauges exist, with freedom `c(z^5)`; (d) if all `b_a`, `b~_a` lie in `Q((u))`, `P^\natural` can be chosen
  over `Q`, and `P^\natural c` has rational coefficients iff `c` lies in `Q[[z]]`. Proof: reduction G2_3 s.5.2, existence s.5.3 (order of the
  choices), leading term and uniqueness s.5.4, field s.5.5 [G]. Renames (G2R_2 s.3 C5): the free constant (G2's `gamma`) -> the coefficient
  sans `c_n`; `beta = u Bprod` -> a local letter; the field -> `\Bbbk`; `S` -> `zeta^(-gr)`; `V(A, A~)` -> `\Gau(A, A~)`; `e_c`, `h_j`, `phi_j`,
  `rho_a`, `rho^U_a` local; `kappa^U_j` -> `chi^U_j`. Optional: index shifts `m = 1..5` with `gl^(5) = gl^(0)` (G2R_2 s.3 R1). Every earlier
  `\draftnote` (G1-G5, G7, `C[[u]]`) is deleted. Evidence (existence, uniqueness): appendix C only.
* `lem:gau-comm` (lem; S14): G2_2 s.3 (a)-(g) [G]: the commutant `\Bbbk((u))[\mathcal Y]`, `\mathcal Y = A~^-1`; `gl_5 = \Bbbk((u))[\mathcal Y] (+) im ad A~`;
  the projection `pr`; integrality and the pole bound (f); (g) `\Bbbk((u))[\mathcal Y]` is a field (Eisenstein); regular semisimplicity is
  true but not needed.
* `lem:gau-lift` (lem; S15 = UR "Lemma 1"): G2_2 s.4.1 [G]. `lem:gau-scal` (lem; S16 = UR "Lemma 2"): G2_2 s.4.2 [G] with the two
  projection facts `pr((\mathcal Y^j)') = -(j/5)(log Bprod)' \mathcal Y^j`, `pr(dlog Theta_U \mathcal Y^j) = (tr dlog Theta_U/5) \mathcal Y^j`.
* `lem:gau-lead` (lem; new): the leading term of a solution of the gauge equation, G2_3 s.5.4 [G]; placed before the uniqueness clause.
* `lem:gau-degree` (lem; new): the Z/5-degree of the gauge equation, G2_4 s.6.1 [G]; placed before `prop:gau-graded`.
* `prop:gau-graded` (prop; S17): graded regular gauges exist; they are `P c(z^5)`; for graded `P`, `P_11(0,z)` is a unit of `Q[[z^5]]`. Proof
  G2_4 s.6.2 [G] (step 1 by `lem:gau-degree`: no `zeta`, no averaging); one sentence: gradedness is equivalent to `zeta^(-gr)`-covariance
  and to invariance under the colour average (G2_4 s.6.4 (i)); the covariance statements of step 2 [W*] CS s.1, [W] UR s.2.
* `rmk:gau-Uremarks` (rmk; kept under D8, reduced): (a) what the proof uses (G2_3 s.5.6); the role of the simple pole: for two systems
  without a pole every index resonates and the solution module has rank 5 (the one-line reason of G2R_2 s.3 R3; the example is computed,
  appendix C; not used); the spectral-curve reading is not used.

### sec:gau-van — (R1), the diagonal, the vanishing solutions, (G)
* `lem:gau-R1` (lem; S21): for every regular gauge, `r E_12 P(0) = r~ P(0) E_12`; column 1 of `P(0,z)` is `P_11(0,z) e_1`, row 2 is
  `P_22(0,z) e_2^T`, `r P_22(0,z) = r~ P_11(0,z)`; here `P_22(0,z) = P_11(0,z)`. Proof S21 steps 1-2 [W] (UR s.5; GI s.3).
* `lem:gau-diag` (lem; S27): for a graded regular gauge `prod_a p_(0,a)` is constant, with the general value `p_(0,1)(0)^5 prod_a rho^U_a`
  (G2_4 s.7); under (N0) (`b_a(0) = b~_a(0)` for `a != 1`, `r = r~`; true in the conifold frames, G2 s.2) it is `1` for `P_can` and
  `p_(0,a)(0) = 1` for all `a`; `p_(0,0)(u) = tau^(-1/4)/varpi_0`, `Ghat(u,0)_00 = 1/varpi_0`. Proof S27 steps 1-6 [W]/[W*] (GI s.5; GR s.2),
  step 4 = G2_4 s.7 [G]; the second proof from `lem:gau-lead` (G2_4 s.7, "Second proof") may replace steps 1-5. Drop `g_0 g_3 = 1`
  (undefined `g_a`; REV_2 M13 (f)). Delete `\draftnote{G8}`.
* `def:gau-van` (defn; DB D13): `I^C_m(v,z) = Res_(p=0) p^m e^(-W_v(p)/z) dp`, `Y^C := (I^C_0, -I^C_(-1), -I^C_3, I^C_2, I^C_1)`; the quintic
  `Y^C_Q := Sum_n z^-n y^(n)`, `y^(0) = -e_1`, `y^(n+1)_a = Int_0^u b_a y^(n)_(a+1) du`; both characterised by the initial value `-e_1` (GR s.2).
* `lem:gau-vanmod` (S22 = GI Lemma V), `lem:gau-vanshape` (S23: the closed form `[z^-n] Y^C_b = eps_b (-1)^n v^(n-q)/(4^q q! (n-q)!)`, `q =
  floor(n/5)`, for `n = nu_b + 1 mod 5`, else 0; used by `lem:tra-conv` (c), G3R K2), `lem:gau-vanq` (S24 = GI Lemma Q; also for any one-pole
  system with `T_A := Int_0^u b_0` in place of `t_c`; (d) `[z^-1] Y^C_Q = -t_c e_0`, `[u^0](Y^C_Q)_0 = 0`), `lem:gau-uniq` (S25 = GI Lemma 1).
  Proofs DB s.3 S22-S25 [W]; S24 step 4 [W*], step 7 [W*] (PS-C7). Correction S24 (b): the pole of `b_1` meets `y^(n)_2` of order
  `(4n+4)/5 >= 4` (GI writes `>= 1`). `rmk:gau-monodromy`: dropped (D8).
* `thm:gau-G` (thm; S28, the gauge identity "(G)"; never "Theorem G"): for every regular gauge `P` from `A~` (pulled-back model) to `A`
  (quintic), `P(u,z) Y^C(v(u), z) = P_11(0,z) Y^C_Q(u,z)` in `Q((z))[[u]]^5`. Cor. (a): `P_can := P P_11(0,z)^-1` is a graded regular gauge, the
  unique regular gauge with `P Y^C = Y^C_Q` (proof G2_4 s.8 [G]: `thm:gau-U` (b), (G) for `P_can`, `lem:gau-vanq` (c)); (b) `P Y^C` has no
  positive power of `z` iff `P_11(0,z)` is constant; (c) the same for any pair of one-pole systems with equal products. Proof S28 steps 1-6
  [W] (GI s.3; GR s.2). Delete `\draftnote{G9}`.
* `def:gau-can` (defn; DB D12): the canonical gauge `P_can`, `P_11(0,z) = 1`; graded and colour covariant (G2_4 s.6.3, s.6.4 (ii)).

### sec:gau-symp — symplecticity
* `prop:gau-symp` (prop; S29; GI's "Proposition S", renamed): for every regular gauge between self-dual one-pole systems with equal products
  and `r = r~`, `P(u,-z)^T eta P(u,z) = c(z) c(-z) eta`, `c(z) := P_11(0,z)`; `P_can` is symplectic; the symplectic regular gauges are `P_can
  c(z)` with `c(z)c(-z) = 1`; symplecticity sees only the even part of `log c`. Proof S29 steps 1-5 [W] (GI s.4; GR s.2); step 3 cites
  `thm:gau-U` (b) for the pair `(A~, A~)`: every solution of `z T' = A~ T - T A~` in `gl_5(Q[[u]])[[z]]` is `t(z) 1` (G2_3 s.5.4, consequence 2;
  G2_6 s.11 item 6; invertibility of `T_0` is not needed).
* `cor:gau-Gstar` (cor; S30): `Y^C_Q(-z)^T eta P_can(z) = Y^C(-z)^T eta`. Proof [W] GI s.4. `rmk:gau-res`: dropped (D8; values in appendix C).

### sec:gau-R0 — HPF, formal WKB uniqueness, colour covariance, R0, the rows in the ring, all colours
* `lem:gau-HPF` (lem; DD1 lem:pol-HPF): with the operators `\cO_e := 1 + zhat(1-u)((4/5) theta_zhat - e/5 - theta_u)` and the iterates
  `\Phi^{[j]} := \cO_(-j) \Phi^{[j-1]}` (REV R8 (c): not `A_e`, `Phi_j`), `L e^(u_0/z) Phi(zhat, u)` solves PF5 iff (HPF) holds; `[zhat^(n+1)] HPF =
  (1-u)[(4n - 5 theta_u) ph_n + r_n]`, `r_n` from `ph_0..ph_(n-1)`; for `n = 5m` the coefficient `frak f_m := [u^(4m)] ph_(5m)` is free and
  `[u^(4m)] r_n = 0` is required. Proof DD1 s.3.1 steps 1-5 [W] (POLY s.1; PR s.2; PB s.2.2).
* `lem:gau-poly` (lem; DD1 lem:pol-poly): a solution of HPF with `ph_0 = 1` and every `[zhat^n] Phi` in `Q[[u]]` has every `ph_n` a polynomial in
  `u`; such solutions are `{Phi^(0) c(z) : c in 1 + z^5 Q[[z^5]]}` once one exists. Proof DD1 s.3.1 steps 1-6 [W]; step 6 [T] (POLY s.2).
* `lem:gau-wkb` (lem; formal WKB uniqueness; DD1 lem:glu-wkb; G2R_2 s.3 C2): for every colour `i`, over any field containing `Q(zeta)`: two
  colour-`i` WKB columns of PF5 with the same exponential `e^(u_i/z)` and leading term `L` differ by a constant series `\mathsf k(z)` in `1 +
  z Q(zeta)[[z]]` (`1 + z Q[[z]]` for `i = 0`), in `Q[[q]]` and in formal series in `z` over `Q((s))` (REV_2 M13 (c): not `d(z)`). Proof:
  colour 0 [T] PS-C10, PS s.4.2 (b), ER s.4 (`T_1(x^m) = zeta^-i (1-m)(x^(m-1) - x^(m+4))`: the exponent is a simple root of the symbol, the
  transport operator is first order with kernel `Q L`; induction on the first differing coefficient); colour `i`: the substitution `z ->
  zeta^-i z` maps `PF5^(u_0)` to `PF5^(u_i)` (every `z`-monomial of PF5 has degree 0 or 5; GS s.2.5), or the colour-`i` transport operator
  `T(L g) = 5 zeta^(4i) D g` (G2R_2 s.2.10).
* `lem:gau-colmod` (lem; new; G2_5 s.9.3 [G]): the colour covariance of the model columns: the colour-`i` T-system (`\mathcal T^(i)_k`,
  `Fs_0(0) = 1`) has the solution `zeta^(ib) Fs_b(zeta^-i t)`, unique (pivot `4n zeta^(4i)`), so the rows `Rhat^mod,(i)` of `def:set-M0` are a
  WKB solution with exponential `e^(U_i/z)` and `Q^(m,i)_b(u, z) = zeta^(i(b+1)) Q^(m,0)_b(u, zeta^-i z)` at `v = v(u)` (G2R_2 s.2.9).
* `lem:gau-colrows` (lem; new; G2_5 s.9.4 Steps 1-5 [G] with G2R_2 s.3 C3 and `psi_i = L Phi_i`): inputs `prop:set-cyclicLR` (the cyclic
  shape only), `D u_i = zeta^i L`, the WKB form of row 0 (`cor:set-cyclicPF5` (b)) with `Phi_i(z, 1) = h_(w_i)(z)` (`thm:set-qrr`, G1_4 s.3.1),
  `lem:gau-wkb`, `lem:set-hrot`; statement: `Phi_i(zhat, u) = Phi_0(zeta^-i zhat, u)` and `Q^(i)_(a,n) = zeta^(i(a+1-n)) Q^(0)_(a,n)` as
  identities of formal series in `q` over `Q(zeta)`. Step 2 reads "the inverse of the substitution of Step 1" (C3). Step 3 is where (QRR)
  enters s03. These are the covariance clauses of `def:set-rows` (GS s.1 item (a)).
* `thm:gau-R0` (thm; DB S20 with DD1 thm:glu-main (a)-(d)), colour 0: (i) `P_can Q^m(v(u))` is a formal WKB solution matrix of the quintic
  system with exponentials `e^(u_i/z)` (`lem:set-canon`); hatted `Ghat(u,zhat) Rhat^mod(tau zhat)`, a power series in `u` whose coefficients are
  power series in `zhat` (`mu(b,k) >= 0`, `lem:mai-Ghat`); (ii) `Phi^can := varpi_0 x` (colour 0, row 0 of the hatted `P_can Q^m`) solves HPF
  (`cor:gau-cycsys`), `ph^can_0 = 1` (`lem:gau-diag`), so every `ph^can_n` is a polynomial (`lem:gau-poly`) and `A^can_n` lies in `Q[L, 1/L]`;
  (iii) the actual colour-0 column (WKB form: `cor:set-cyclicPF5` (b)) satisfies `A^can = A^act D_0(z)` in `Q[[q]]` (`lem:gau-wkb`), so
  `A^act_n` lies in `Q[L, 1/L]`; (iv) in `Q[L, 1/L][[z]]`, embedded in `Q((s))[[zhat]]`, `Q^(0) D_0 = P_can Q^(m,0)(v(u))` with `D_0` in `Q[[z]]`
  constant in `u`, the SAME exponential `e^(u_0/z)` stripped on both sides; row 0: `Phi^act(zhat, u) D_0(s^4 zhat) = Phi^can(zhat, u)`. Proof:
  UR s.4; EVID s.2.3 steps 1-5; PS s.4.1-s.4.2 (a)-(d) [W]; PB s.3 [W]. Corrections PS-C8 (ring of comparison; same exponential), RS1 A4
  (log-free actual rows are derived), RS2 items 13, 24. No normalisation of the actual rows is used; colour 0 uses no (QRR).
* `cor:gau-rowsNN` (cor; new; REV R1 item 3; G1_3 s.2.10 (c) [G] with G2_5 s.9.4 Step 6): the recursion coefficients `(C_1 kappa^LR_a)^-1 =
  (1/C_1, C_0^2 C_1^2 L^-5, 1/C_1, 1/C_0)` (`a = 0..3`), `D L = (L^6 - L)/5`, `thm:gau-R0` (iii) and `lem:gau-colrows` put every row coefficient
  `Q^(i)_(a,n)` into `\Nring (x) Q(zeta)`, free of `\logq`; hence edges, leaves and graph terms lie there (G5_2 s.2.2), and the covariance
  clauses survive the conifold expansion and hatting. Sources RS1 A4-A5, PS s.4.2 (c).
* `cor:gau-colours` (cor; G2_5 s.9.1-s.9.2 [G] with G2R_2 s.2.9, s.3 C4): for every `i`, `D_i(z) := D_0(zeta^-i z)` satisfies `Q^(i) D_i =
  P_can Q^(m,i)(v(u))` in `Q(zeta)[L, 1/L][[z]]` inside the domain `Q(zeta)((s))[[z]]`, and is the unique such series; hence `Q D = P_can Q^m`,
  `D = diag(D_i)` (`eq:gau-R0`); `D_0` in `Q[[z]]`, `D_i` in `Q(zeta)[[z]]` (in `Q[[z]]` for `i != 0` iff `D_0` in `Q[[z^5]]`); `D_0 = 1` implies `D = 1`.
  Proof: the three-line computation with `lem:gau-colmod`, `lem:gau-colrows`, `cor:gau-rowsNN` and the gradedness of `P_can`, read without
  its first equality (G2R_2 s.2.9). Delete `\draftnote{G6}`.
* `def:gau-Fg` (defn): (F_g) is `D_0 = 1 mod z^(3g-2)`; at genus `g` only `P` through `z^(3g-3)` enters (ST2 s.2; FR (I2)).

Not stated (superseded; rule 5): UR's "(H) purity", "canonical normalisation (CONJECTURAL)", GI s.6's conjecture; they are consequences
of `thm:gam-D` (`cor:mai-natfact`, `cor:mai-actual`). Length 10 pp.
