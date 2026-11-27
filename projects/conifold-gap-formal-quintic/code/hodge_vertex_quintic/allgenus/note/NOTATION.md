# One notation for the note (part 1 of 2)

Working file for the section writers; not part of the note. Revised after the gap phase (record: `design/DESIGN_CHANGES.md`). Part 2
(`NOTATION_2.md`): tables 6-9, the renaming dictionary source -> note, the fixed conventions, and the macros of `NOTE/macros.tex`. Every
recurring symbol has a macro; a symbol used in one proof only may be declared at the top of the section file (with `\providecommand`
and the section prefix if it could clash) and must be said to be local.

**Rules of precedence.** (1) The letters of the task statement are kept: `c_g`, `tau(u)`, `M = M_0(tau(u) zhat)`, `t_c`, `K`, `K_0`,
`c(z; L)`, `D_0` (as sans-serif `D`), (QRR), (F_g), `Theorem U`, `Theorem P`, `Theorem LMI`, `Theorem Gamma`. (2) Otherwise the letter of the
source of record is kept when it does not clash. (3) A clash is resolved by the tables below; a renamed symbol is never used under its old
name; the dictionary of part 2 lists every renaming. (4) Dummy indices: `a, b` in Z/5 for flat components (rows), `i, j` in Z/5 for colours,
`k, m, n` for orders, `g` for the genus of a graph, `g_v` for the genus of a vertex (REV_2 M13 (a)), `h` for the genus index of the model
(`ML_K(h)`, `F^m_h`, `c_h`, `omega_(h,n)`); the series `h(z)`, `h_n`, `h_w` keep their meaning and never meet the genus index in one formula.

## 1. Variables, points, operators, rings

| note | macro | meaning | defined | sources |
|---|---|---|---|---|
| `q` | `q` | parameter of the formal quintic (degree variable in the I-function) | `def:set-genus0` | q |
| `D` | `\Dq` | `q d/dq = -(1-u) d/du` | `def:set-genus0` | D (never the normalisation series) |
| `\partial_T` | `\partial_T` | `C_1^-1 D`, the derivative along the mirror coordinate `T` (s02) | `prop:set-cyclicLR` | d_T (G1) |
| `u` | `u` | `1 - 5^5 q`; conifold point `u = 0`, large radius `u = 1` | `def:set-points` | u (TRM's `u` = residue: rename `v`) |
| `s` | `s` | `u^(1/5)` | `def:set-points` | s (other `s`: renamed, part 2) |
| `L` | `L` | `u^(-1/5) = (1 - 5^5 q)^(-1/5) = 1/s`; in s08-s09 the independent variable | `def:set-points`, `def:gam-c` | L, x (ER, G5) |
| `zeta` | `\zeta` | `e^(2 pi i/5)` | s02 | zeta, xi (GJR17) |
| `z` | `z` | descendant variable, constant for `D` | s02 | z (XY's curve coordinate: `p`) |
| `zhat` | `\zh` | `z/s^4 = L^4 z`, the hatted variable | `def:set-points` | zhat, w (PB, POLY, PS) |
| `theta_u`, `theta_zhat` | `\theta_u`, `\theta_{\zh}` | Euler operators `u d/du`, `zhat d/dzhat` | `lem:gau-HPF` | th_u, th_w; UR's `th = z d/du` is written `z d_u` |
| `\partial_t` | `\partial_t` | `(Dt)^-1 D` in s03 | `def:gau-conc` | CYC's `th` |
| `\Cstar` | `\Cstar` | `C \ {0, 1, zeta, .., zeta^4}` | `def:gam-c` (s08) | CP's `Lam` |
| `q(L)` | `q(L)` | `(1 - L^-5)/5^5` | `def:gam-c` (s08) | q(L) |
| `\AL` | `\AL` | `Q[L, L^-1, (L^5 - 1)^-1]` with `D L = (L^6 - L)/5` | `prop:gam-rational`, `thm:con-Pa` | A_L (G5) |
| `\jmath` | `\jmath` | the embedding `\AL -> Q((s))`, `L -> s^-1` | `thm:con-Pa` | iota (G5; clashes with CGYZ26's `iota`) |
| `KK` | `\bbK` | `Q((u))` with `D` (s03) | `def:gau-conc` | CYC's `K`, UR's `K = O[1/u]` |
| `\Bbbk` | `\Bbbk` | a field of characteristic 0 (Theorem U) | `thm:gau-U` | k (G2) |
| `ell` | `\ell` | indeterminate with `D ell = 1 - 1/u` (`ell = log u`) | `def:gau-conc` | CYC's `l` |

## 2. The formal quintic, genus-0 data, rows, (QRR)

| note | macro | meaning | defined | sources |
|---|---|---|---|---|
| `lambda_i` | `\lambda_i` | torus weights; set to `zeta^i` (`lambda = 1`) | `def:set-fq` | lambda_i |
| `p_i` | `p_i` | fixed points of `P^4` (always with an index; the model coordinate is `p`) | `def:set-fq` | p_i |
| `N_(g,d)`, `F_g`, `\Qn` | `N_{g,d}`, `F_g`, `\Qn` | invariants; `F_g = Sum_(d>=0) N_(g,d) \Qn^d`, `\Qn` the degree variable (`\Qn = e^(T(q))` along the mirror map) | `def:set-fq` | N_{g,d}, F_g, Q (REV R8 (g)) |
| `I_i(q,z)` | `I_i` | I-function at `p_i` | `def:set-genus0` | I_i |
| `C_0`, `I_1`, `T`, `C_1` | `C_0, I_1, T, C_1` | large-radius periods, mirror map `T = I_1/C_0`, `C_1 = DT` | `def:set-genus0` | C0, C1 (in a conifold frame write `varpi_0`, `D t_c`) |
| `\logq` | `\logq` | the formal `log q` (`D \logq = 1`); its conifold image `log(1-u) = -\Xi` | `def:set-confexp` | Lam (G5), ell (G1) |
| PF5, PF4 | `\PFfive`, `\PFfour` | the fifth- and fourth-order operators | `def:set-genus0` | PF5, PF4 |
| `kappa^LR_a` | `\kappa^{\mathrm{LR}}_a` | `(1, L^5/(C_0^2 C_1^3), 1, C_0/C_1, C_0/C_1)`, `H * H^a = kappa^LR_a H^(a+1)` | `prop:set-cyclicLR` | kappa (G1), c_a |
| `mu_i` | `\mu_i` | eigenvalues `zeta^i L/C_1` of `H*` (local to `cor:set-cyclicPF5`) | `cor:set-cyclicPF5` | mu_i (G1) |
| `\Wr` | `\Wr` | Wronskian (of PF4 solutions; at large radius `= L^10`) | `lem:set-wronskian`, `def:gau-conc` | W (G1), Wr |
| `Delta_i` | `\Delta_i` | `1/eta(e_i, e_i) = C_0^2 zeta^(3i)/L^2` | `cor:set-cyclicPF5` | Delta_i |
| `e_i`, `e~_i` | `e_i`, `\tilde e_i` | idempotents, normalised canonical basis | `cor:set-cyclicPF5`, `def:set-frob` | e_i (standard basis vectors are bold `\mathbf e_a`) |
| `\bbH` | `\bbH` | the state space, flat basis `H^a` (REV_2 M13 (e)) | `def:set-frob` | H = Q^5 |
| `W_i^0` | `W^0_i` | `zeta^(3i)` | `def:set-frob` | W_i^0 |
| `eta` | `\eta` | `eta(H^a, H^b) = 5 delta_(a+b=3 mod 5)`, `eta^-1 = eta/25` | `def:set-frob` | eta (GI's `eta_n`: not used) |
| `\TFT` | `\TFT` | the TFT of the formal quintic, `\TFT_(g,n)(e~_i, ..) = Delta_i^((2g-2+n)/2)` (REV R6) | `def:set-frob` | omega |
| `u_i`, `u_i^c` | `u_i`, `u_i^{\mathrm c}` | canonical coordinates, `D u_i = zeta^i L`; value at `u = 0` | `cor:set-cyclicPF5`, `lem:set-canon` | u_i |
| `R(z)` | `R` | the R-matrix of the formal quintic, Givental's convention | `def:set-rows` | R |
| `Q_a^(i)` | `\Qr^{(i)}_a` | rows `eta(H^a, R(z) e~_i)` (calligraphic `Q`) | `def:set-rows` | Q, P/Delta^(1/2) (never blackboard) |
| `Phi_i`, `ph_n` | `\Phi_i`, `\mathrm{ph}_n` | `Phi_i := C_0 Q^(i)_0/(zeta^i L)` (a definition), its coefficients | `def:set-rows` | Phi, ph_n; `psi_i = L Phi_i` |
| `A_n` | `A_n` | flat amplitude `[z^n](C_0 Q^(0)_0) = L^(4n+1) ph_n(L^-5)` | `def:set-rows` | A_n (psi_n in ER) |
| `Phi^act`, `A^act` | `\Phiact` | the actual (QRR-normalised) colour-0 column | `thm:gau-R0` | Phi^QRR, A^QRR |
| `bk_k` | `\bk_k` | `B_2k/(2k(2k-1))` (sans-serif b) | `thm:set-qrr` | b_k |
| `w_i` (bold) | `\bw_i` | `{zeta^i - zeta^j}_(j != i) U {-5 zeta^i}` | `thm:set-qrr` | w_i, wt_i |
| `h_w`, `h`, `h_n` | `\hw`, `h`, `h_n` | `exp(-Sum bk_k z^(2k-1) Sum_(a in w) a^(1-2k))`; `h = h_(w_0)`; `h_n = [z^n] h` | `thm:set-qrr` | h_w, Theta (FQ2) |
| `T_k` | `\mathsf T_k` | `5^(1-2k) + Sum_(j>=1) (zeta^j - 1)^(1-2k)` | `thm:gam-main` | T_k (CP3) |

## 3. Frames, conifold expansion, the natural variable

| note | macro | meaning | defined | sources |
|---|---|---|---|---|
| `f_0, f_1, f_2` | `f_0, f_1, f_2` | Frobenius basis at `u = 0`, `f_k = u^k + O(u^3)` | `lem:set-frobbasis` | f_k (CYC's `f_i = y_i/y_0` is `phi_i`) |
| `varpi_1` | `\varpi_1` | vanishing period `f_1 + (7/10) f_2`, scale `u + O(u^2)` | `lem:set-frobbasis` | g (never: clash with the genus) |
| `gamma` | `\gamma` | the number `7/10` | `lem:set-frobbasis` | gamma, gamma_g (G2's free constant: sans `c_n`) |
| `rho_3` | `\rho_3` | power-series part of the log solution | `lem:set-frobbasis` | h (CYC) |
| `(alpha, beta)` | `\alpha, \beta` | frame parameters, `varpi_0 = f_0 + alpha varpi_1 + beta f_2` | `def:set-frames` | (alpha, beta); `(a, b)` labels only in appendix C |
| `varpi_0` | `\varpi_0` | conifold-frame period | `def:set-frames` | varpi_0, C0 (frame), y_0 (CYC) |
| `t_c` | `\tc` | `varpi_1/varpi_0` | `def:set-frames` | t_c, t |
| `gamma_t`, `w_1` | `\gamma_t`, `w_1` | `[u^2] t_c`, `[u^1] varpi_0 = alpha` | `def:set-frames` | same |
| `\Nring` | `\Nring` | the presented differential ring (generators `L^(+-1)`, `\logq`, `c_0..c_3`, `i_0..i_2`; `c_0`, `\mathsf W` inverted) | `def:set-confexp` | NN (G5; not DP's sans `N`) |
| `\mathsf W` | `\mathsf W` | `c_0^2 + c_0 i_1 - c_1 i_0`, `ev(\mathsf W) = C_0^2 C_1` | `def:set-confexp` | W (G5) |
| `ev` | `\ev` | evaluation `\Nring -> ` functions of `q` | `def:set-confexp` | ev |
| `\Cexp_(alpha,beta)` | `\Cexp_{\alpha,\beta}` | the conifold expansion `\Nring -> Q((s))` (never the code name) | `def:set-confexp` | Cf (code name) |
| `tau(u)` | `\tau(u)` | `1/2F1(1, 4/5; 9/5; u)`. **`tau` means only this.** | `def:set-tau` | tau (CYC's automorphism: `sh`; `tau_k(H^b)`: `psi^k H^b`) |
| `v(u)` | `v(u)` | `u tau^(-5/4)`; also the model's parameter `v` | `def:set-tau` | v |
| `varsigma` | `\varsigma` | `v^(1/5) = s tau^(-1/4)` | `def:set-tau` | sigma (UR): the Mellin variable keeps `sigma` |
| `V` | `V` | `tau^(-5/4)` (s08) | `lem:con-cfu` | V (PB's residue `V` is `VV`) |

## 4. Cyclic systems, gauges, normalisation

| note | macro | meaning | defined | sources |
|---|---|---|---|---|
| `kappa_a` | `\kappa_a` | coefficients of `zD Y_a = kappa_a Y_(a+1)` in a frame, `kappa = \Cexp(C_1 kappa^LR)` | `cor:gau-cycsys`, `prop:gau-quintic` | c_a (UR, CYC, D1.4) |
| `b_a`, `b~_a` | `b_a`, `\tilde b_a` | coefficients of `z d_u Y = A Y`, `b_a = -kappa_a/(1-u)`; pulled-back model | `def:gau-onepole` | b_a, bt_a |
| `A`, `A~` | `A`, `\tilde A` | `Sum_a b_a E_(a,a+1)` | `def:gau-onepole` | A, At |
| `Bprod` | `\Bprod` | `prod_a b_a` (calligraphic B); (EP): `Bprod = prod_a b~_a` | `def:gau-onepole` | B (UR; never the Bergman kernel or the concomitant) |
| `r`, `r~` | `r` | residues | `def:gau-onepole` | r, rt |
| (N0) | | `b_a(0) = b~_a(0)` for `a != 1` and `r = r~` (true in the conifold frames) | `lem:gau-diag` | (N0) (G2) |
| `\Gau(A, A~)` | `\Gau(A,\tilde A)` | all solutions of the gauge equation in `gl_5(\Bbbk[[u]])[[z]]` | `def:gau-onepole` | V(A, A~) (G2) |
| `P`, `P_can` | `P`, `\Pc` | regular gauge, canonical gauge | `def:gau-onepole`, `def:gau-can` | P, P_can |
| `c(z)` (sans) | `\cn(z)` | the freedom `P -> P c(z)`; `c(z) := P_11(0,z)` in (G); its coefficients `\cn_n` | `thm:gau-U`, `thm:gau-G` | c(z); G2's `gamma` |
| local to `sec:gau-U` | | `\mathcal B^{(m)}_a` (partial products), `\mathcal Y = A~^-1`, `\Theta_U = diag(d^U_a)`, `dlog Theta_U`, `\mathfrak d`, `pr`, `chi^U_j`, `P^\natural`, `P^{[delta]}`, `Rot_i`, `zeta^(gr)` | proof of `thm:gau-U` | G2 s.2 (`kappa^U_j` -> `chi^U_j`; `S` -> `zeta^(-gr)`) |
| `Ghat` | `\Ghat` | hatted gauge `Lambda P(u, s^4 zhat) Lambda'^-1` (local diagonal matrices) | `lem:mai-Ghat` | Ghat |
| `G(u, zhat)` | `G` | `Psi(u)^-1 Ghat Psi_0` (only in `cor:mai-natfact`) | `cor:mai-natfact` | G |
| `Y^C`, `Y^C_Q` | `\YC`, `\YCQ` | vanishing solutions (model, quintic) | `def:gau-van` | Y^C, Y^C_Q |
| `I^C_m` | `\IC_m` | `Res_(p=0) p^m e^(-W_v(p)/z) dp` | `def:gau-van` | I^C_m |
| `D(z)` (sans) | `\Dn(z) = \operatorname{diag}(\Dn_i)` | normalisation series of R0; `D_0` in `Q[[z]]`, `D_i` in `Q(zeta)[[z]]` | `thm:gau-R0`, `cor:gau-colours` | D(z), D_i, D_0 |
| `\mathsf k(z)` | `\mathsf k(z)` | the constant series of formal WKB uniqueness (local) | `lem:gau-wkb` | d(z) (REV_2 M13 (c)) |
| (F_g) | `\Fg{g}` | `D_0 = 1 mod z^(3g-2)` | `def:gau-Fg` | (F_g); (F) is never used alone |
| `Phi^can`, `Phi^SP` | `\Phican`, `\PhiSP` | canonical and full-Stirling colour-0 WKB functions | `thm:gau-R0`, `def:con-PhiSP` | same |
| `\cO_e`, `\Phi^{[j]}` | `\cO_e`, `\Phi^{[j]}` | HPF operators `1 + zhat(1-u)((4/5) theta_zhat - e/5 - theta_u)` and iterates `\cO_(-j) \Phi^{[j-1]}` | `lem:gau-HPF` | A_e, Phi_j (REV R8 (c)) |
| `frak f_m` | `\fc_m` | free coefficient `[u^(4m)] ph_(5m)` | `lem:gau-HPF` | f_m, f_(5m) (EVID), x_m (GAM) |
| `CC(f,g)`, `JJ`, `CCmat` | `\Conc`, `\Jet` | concomitant, jet matrix, Gram matrix | `def:gau-conc` | B(f,g), P_ab, Bmat |
| `sh` | `\sh` | automorphism `ell -> ell + 1` | `def:gau-conc` | tau (CYC) |
| `TT_k`, `TT~_k` | `\Upsilon_k`, `\widetilde\Upsilon_k` | transversality sums (with `D`, with `\partial_t`) | `def:gau-conc` | T_k, tau_k (CYC) |
| `M_cyc` | `\Mcyc` | the factorised operator | `def:gau-conc` | M (CYC) |
| `Yuk` | `\mathrm{Yuk}` | `\partial_t^2 phi_2` (special geometry) | `lem:gau-sg` | Y (CYC) |
| `m_5` | `m_5` | the scalar of the other normalisation (value in appendix C only) | `sec:sgn-norm` | m5, c_5 |

## 5. The model (Landau-Ginzburg side)

| note | macro | meaning | defined | sources |
|---|---|---|---|---|
| `W_v(p)` | `\Wv(p)` | `p^4/4 + v/p` | `def:set-model` | W, lambda = -(x^4/4 + v/x), x^4 + v/x |
| `p` | `p` | model coordinate (LG variable and curve coordinate) | `def:set-model` | x (GI, UR), p (ST2), z (XY) |
| `I_m(v,z)` | `I_m` | `Int p^m e^(-W_v(p)/z) dp` | `def:set-model` | I_m |
| `nu_a` | `\nu_a` | charges `(0,-1,3,2,1)`, `nu_a + nu_(3-a) = 2` | `def:set-model` | c_a (hatted exponents) |
| `eps_b` | `\varepsilon_b` | signs `(1,-1,-1,1,1)`, `Y_b = eps_b I_(nu_b)` | `def:set-model` | eps_b, N_a |
| `bmod_b` | `b^{\mathrm m}_b` | `(1, -1/v, 1, -1, -1)` | `prop:set-modrows` | bmod |
| `Fs_b(t)` | `\Fs_b(t)` | model row series, `Fs_0 = 1 - (3/20)t + ..`, `Fs_b(0) = eps_b` | `prop:set-modrows` | F_b(t) |
| `T_k` (model) | `\mathcal T_k`, `\mathcal T^{(i)}_k` | `-zeta^i + (t/5)(k - 4t d/dt)` | `prop:set-modrows` | T_k (UR), calT^(i)_k (G2) |
| `Rhat^mod`, `Psi_0` | `\Rhat^{\mathrm m}`, `\Psi_0` | `zeta^(i(b+1)) Fs_b(zeta^-i t)`; its value at `t = 0` | `def:set-M0` | R^, Psi0 |
| `M_0` | `\Mo` | `Psi_0^-1 Rhat^mod(t)` (Teleman form); `M_0(-t)^T M_0(t) = 1` | `def:set-M0`, `lem:set-M0sympl` | M_0 |
| `Q^m` | `\Qr^{\mathrm m}` | model WKB columns at `v(u)` | `thm:gau-R0` | Q^m |
| `Omega^mod_v`, `\TFT^m_v` | `\Omega^{\mathrm{mod}}_v`, `\TFT^{\mathrm m}_v` | model CohFT; model TFT, defined by `Deltahat_i^(1/2) := zeta^(4i) varsigma` | `def:set-modcohft` | Omega^mod, omega^m |
| `\varkappa` | `\varkappa` | the dilaton and TFT rescaling scalar (`\varkappa` in `1 + m`) | `lem:set-dilaton`, `def:flh-outer` | kappa (G1, S15; REV R8 (i)) |
| `Fbar_h`, `\dd^n Fbar_h` | `\overline F_h` | dilaton-shifted ancestor potential; its `n`-th derivative (roman `\dd`) | `lem:set-dilaton` | Fb_h, d^n Fbar (the dilaton factor `d` stays italic) |
| `F^m_h`, `c_h` | `\Fm_h`, `c_h` | model free energy `c_h v^(2-2h)`; constant | `thm:mod-gap`, `cor:mod-trF` | F_h, Phi_h |
| `wt(k,a)` | `\wt(k,a)` | `(4k + nu_a - 4)/5` | `lem:set-modhom` | w(k,a), beta (TRM) |
