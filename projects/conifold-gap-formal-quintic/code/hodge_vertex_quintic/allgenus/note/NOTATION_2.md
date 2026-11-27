# One notation for the note (part 2 of 2)

Continuation of `NOTATION.md` (rules of precedence and tables 1-5 there).

## 6. Tower, slice, transport, graph sums

| note | macro | meaning | defined | sources |
|---|---|---|---|---|
| `\cH`, `\cH_+` | `\cH`, `\cH_+` | graded setting: `\bbH((z^-1))`, `\bbH[z]` | `conv:tra-forms` | Hc, Hc_+ |
| `\cH^u`, `\cH^u_+`, `\cH^u_-` | `\cH^u`, .. | u-adic setting: `\bbH (x) Q((z))[[u]]`, `\bbH (x) Q[[u, z]]`, `\bbH (x) z^-1 Q[z^-1][[u]]` | `conv:tra-forms` | HH^u (G3) |
| "summable", "of finite type" | | families and maps of the summation lemma | `conv:tra-forms`, `lem:tra-sum` | G3 s.2.2 |
| `\Theta = 1 + N` | `\Theta` | the operator of the shift rule on `(Pi_i(q))_i` (local letters of `sec:tra-tower`; `N` is not DP's sans `N`) | `lem:tra-theta` | Theta, N (G3) |
| `\Omega_s(f,g)` | `\Omega_{\mathrm s}` | `Res_z eta(f(-z), g(z)) dz` | `conv:tra-forms` | Omega (never the CohFT letter) |
| `psi^k H^b` | `\psi^k H^b` | an insertion; as a vector `z^k H^b` | `conv:tra-forms` | tau_k(H^b) |
| `Pi_i` | `\Pi_i` | vanishing tower; `Pi_0` the vanishing covector | `def:tra-tower` | Pi_i, PiC, V^(k) = (-1)^(k+1) Pi_(k-1) |
| `K`, `K_0` | `\sK`, `\sKo` | the slice and its kernel part (calligraphic K) | `def:tra-slice` | K, K_vc, K_0, T_v, Lambda_v |
| `mu` | `\mu` | the local scalar of the affine slice, `Pi_0(-z 1 + n + mu H) = v + mu` | `def:tra-slice` | s (REV R8 (a)) |
| `frame(k,a)` | `\fr(k,a)` | the frames of `K_0` (ten trivial ones, `k < nu_a + 1`) | `lem:tra-frames` | frame(k,a) |
| `e` (sans) | `\se` | Euler-dilaton vector `z 1 + v H = frame(1,0)` | `lem:tra-frames` | e, A1 (UR) |
| `P`, (Sp), (A) | | dressing, symplecticity, annihilation | `def:tra-dressing` | same |
| `T_P`, `d`, `R` | `T_P`, `d`, `R` | transported vanishing coordinate, dilaton factor `(P_0)_33` (italic), `P_0^-1 P` | `def:tra-dressing` | same |
| `E_R`, `q_T`, `qt` | `E_R`, `q_T`, `\tilde q_T` | edge, leaf point, `d q_T` | `def:tra-edge` | same |
| `\ttr_R(z)` | `\ttr_R` | PPZ's translation `z 1 - z R^-1(z) 1` ([PPZ15] Def. 2.13) | `def:set-graphsum` | T(z) (REV_2 M13 (b); `T` is the mirror map) |
| `delta_P` | `\delta_P` | `d T_P - v`, in `u^2 Q[[u]]` | `lem:tra-T` | s (RD2) |
| `vhat` | `\hat v` | `v/d` | `rmk:red-leaf` | vhat |
| `Cont_Gamma`, `F^P_g`, `F_Gamma` | `\mathrm{Cont}_\Gamma` | contraction (including `|Aut Gamma|^-1`), flat sum, flat term for `P_can` | `def:set-graphsum`, `def:red-flatsum` | same |
| `g_v`, `n_v` | `g_v`, `n_v` | genus and valence of a vertex | `def:set-graphsum` | h_v (REV_2 M13 (a)) |
| ML_K(h) | `\ML{h}` | the Model Lemma on `K` at genus `h` (one argument) | `def:red-ML` | ML_K, ML_Lambda, (ML) |
| `T_Gamma`, `G_Gamma` | `T_\Gamma`, `G_\Gamma` | outer term, hatted term | `def:set-outer` | same |

## 7. Topological recursion on the model curve

| note | macro | meaning | defined | sources |
|---|---|---|---|---|
| `x`, `y` | `x`, `y` | `x = W_v(p)`, `y = p` | `def:mod-curve` | x = W, z^4 + v/z, kappa W |
| `omega_(h,n)` | `\omega_{h,n}` | TR correlators, engine convention; `omega_(0,2)` the Bergman kernel | `conv:mod-kernel` | omega, W^(g)_k (EO), B |
| `varsigma_b` | `\varsigma_b` | branch points `zeta^-b varsigma` | `def:mod-curve` | c_b, c_j |
| `iota_b` | `\iota_b` | the local involution of `x` at `varsigma_b` | `conv:mod-kernel` | sigma, sigma_b (REV_2 M13 (d)) |
| `pbar` | `\bar p` | `v/p` | `def:mod-curve` | q (DP) |
| `xi_a`, `xi_(a,d)`, `D_x` | `\xi_a`, `\xi_{a,d}`, `D_x` | dual descendant functions, `D_x = d/dx` | `def:mod-xi` | xif_a, D |
| `q_phi` | `q_\phi` | ancestor vector of a leg | `def:mod-ancvec` | q_Phi |
| `e^TR_j`, `\bar e_j` | `e^{\mathrm{TR}}_j`, `\bar e_j` | the basis of DNOPS (2-30); the frozen basis `e^TR_j/sqrt(-5)` | `thm:mod-tr` | e_j, et_j (G4) |
| `vartheta` | `\vartheta` | the Laplace variable of (2-31), (2-32) (local to `sec:mod-tr`) | `prop:mod-Rtr` | zeta (DNOPS), lam (G4) |
| `L_c`, `G_m` | `\mathsf L_c`, `\mathsf G_m` | the formal Laplace functional, Gaussian moments | `lem:mod-laplace` | L_c, G_m (G4) |
| `-5 kappa` | | a number (DNOPS's TFT constant at `x = kappa W_v`); never the leg map | `prop:mod-Rtr` | Lambda_kappa (G4), Lambda (TRM) |
| `phi^b_k`, `w_b` | `\phi^b_k`, `w_b` | legs and their weights `w_b = nu_b + 1` | `def:mod-legs` | Phi^b_k |
| `varrho_b` | `\varrho_b` | leg-dictionary constants `(-1, 1, -1, -3, -2)` | `prop:mod-LD` | c_b (LD) |
| `V_w`, `FF`, `Lambda` | `V_w`, `\mathbb F`, `\Lambda` | weight spaces, their sum, the leg map | `def:mod-weights` | S_w, F, Lambda |
| `LL`, `MM`, `NN` | `\cL`, `\mathsf M`, `\mathsf N` | leg class `Q[v][p] + Q(v)[x]`, DP's module and its quotient (not `\Nring`) | `def:mod-weights` | L, M, N (DP) |
| `d(w)`, `e(w)` | `\mathrm d(w)`, `\mathrm e(w)` | dimension counts | `def:mod-weights` | d(w), e(w) |
| `XX` | `\cX` | variation field `p/(p^5 - v) d/dp` | `lem:mod-fixedx` | xi (LL) |
| `Phi_ydx` | `\Phi_{y\mathrm dx}` | primitive of `y dx` | `conv:mod-kernel` | Phi |
| `frak r` | `\mathfrak r` | a ray from 0 to oo avoiding the branch points | `def:mod-pair` | gamma (LL) |
| `\Sigma_h` | `\Sigma_h` | `Sum_r (r-1)! c_(h,r) = -4^h B_2h/(2h)` (local to `prop:mod-period`) | `prop:mod-period` | S_h (G4) |
| `F^TR_h` | `\FTR_h` | TR free energies (engine convention) | `conv:mod-kernel` | F_h^TR |

## 8. Hatted theory, connection problem, Gamma identity, formal stationary phase

| note | macro | meaning | defined | sources |
|---|---|---|---|---|
| `Hhat^a`, `Qhat`, `xhat`, `1hat`, `\hat\TFT` | `\hat H^a`, `\hat{\Qr}`, `\hat x`, `\hat 1`, `\hat\TFT` | hatted classes, rows, TFT datum, unit, hatted TFT | `def:set-hat` | same; omega^ |
| `M`, `R'` | `M`, `R'` | `M_0(tau(u) zhat)`, `Q M^-1` | `def:set-outer` | M_nat, R' |
| `F_g^\Nring`, `F_g^cf` | `F_g^{\Nring}`, `F_g^{\mathrm{cf}}` | the genus-`g` sum in the ring; its conifold expansion | `def:mai-Fcf` | F_g^NN, F_g^cf (G5) |
| `Psi(u)`, `R~` | `\Psi(u)`, `\tilde R` | `Qhat(u, zhat = 0)`, `Psi^T eta^-1 Qhat` (only s10) | `cor:mai-natfact` | Psi, R~ |
| `Xi` | `\Xi` | `-log(1-u)` | `lem:con-form` | Lam (PB, POLY), -Lam' (G5) |
| `WW`, `GG`, `frak a` | `\cW`, `\cG`, `\fa` | flat phase `Xi/x + GG(x)`, its regular part, amplitude | `lem:con-form` | W, G, at |
| `WWhat`, `frak ahat` | `\hat\cW`, `\hat\fa` | hatted phase and amplitude | `lem:con-form` (c) | What, B(X) |
| `VV`, `\hat\cV` | `\cV`, `\hat\cV` | the residue `Res_(x=0) e^(-WW/z) frak a dx`; its hatted relabelling | `def:con-PhiSP`, `lem:con-relab` | V (PB), HH (REV R8 (e)) |
| `Y`, `Y_c`, `pi_m` | `Y`, `Y_c`, `\pi_m` | uniformised variable, critical point `tau^(-1/4)`, reduction coefficients | `lem:con-cfu`, `lem:con-red` | same |
| `PF5^L` | `\PFfive^{L}` | PF5 with `zD` replaced by `zD + L` | `thm:con-Pa` | PF5^L (G5) |
| `c(z;L)`, `c_n(L)` | `c(z;L)` | full-Stirling saddle series; coefficients ALWAYS with their argument | `def:gam-c` (s08) | same |
| `SS_0`, `SS_1`, `SS_2k`, `A_St` | `\mathsf S_0`, .. | phase, half-density, Stirling parts, amplitude | `def:gam-c` (s08) | S_0, S_1, S_2k |
| `lambda` | `\lambda` | the logarithm of `-q` making `L` critical (local to s08-s09) | `def:gam-c` | lambda |
| `p_k(L)` | `p_k(L)` | local power sums `Sum_j (zeta^j - L)^-k` | `def:gam-c` | P_k(L) (REV R8 (f)) |
| `[5]_L` | `[5]_L` | `1 + L + L^2 + L^3 + L^4` | `rmk:gam-sevenvar` | Lam1 |
| `s_*, a_*, M_*, phi, FF` | `s_*, a_*, M_*, \varphi, \mathcal F` | LMI data | `def:gam-lmidata` | s*, alpha0, M0, g, F |
| `Q` (LMI) | `\mathsf Q` | chart parameter `e^(-lambda/a_*)` | `def:gam-lmidata` | Q |
| `Ginf` | `\Ginf` | the orbifold value | `thm:gam-ginf` | Ginf |
| `Gamma_St`, `St_m`, `Poch_M`, `beta_n(theta)` | `\Gst` | formal Stirling form, normalised Stirling series, Pochhammer, Bernoulli-polynomial coefficients | `def:fsp-stirling` | same; DLMF shift `h` -> `theta` |
| `SPn`, `SPh` | `\SPn`, `\SPh` | normalised stationary phase; divided by its constant term | `def:fsp-wick` | N (CP), SP^norm (FSL), Nhat |
| `R_n(A)`, `E_(zC)` | `\cR_n(A)`, `E_{zC}` | weighted ring, Gaussian functional | `def:fsp-ring` | same |
| `frak b` | `\mathfrak b` | Fubini factor | `lem:fsp-fubini` | beta |
| `RR_res` | `\cR_{\mathrm{res}}` | residue ring | `lem:fsp-res` | L (FSL) |

## 9. Dictionary source -> note (every renaming; the writers search this table before transcribing)

Inputs: "(S1)" (mirror theorem) -> (M); "(S2)" (Teleman) -> (T); "(S3)" (implicit function theorem) -> not an input (DEC Q2: `lem:fsp-ift`);
the input names printed are exactly (Q), (QRR), (M), (G0), (T). CP's "(S1)" -> `def:fsp-ring`/`lem:fsp-univ`; the sources' "(S2)" (model
identity) -> `lem:con-S2`; both never printed (REV R2).
Letters: `c_a` (hatted exponents) -> `nu_a`; `c_a` (cyclic, D-form) -> `kappa_a`; G1's `kappa` -> `kappa^LR`; `c_b` (leg dictionary) -> `varrho_b`;
`c_b`, `c_j` (branch points) -> `varsigma_b`; `c(z)` (Theorem-U freedom, `P_11(0,z)`) -> sans `c(z)`; G2's `gamma` -> sans `c_n`; `D(z)`, `D_i`, `D_0`
-> sans `D`; UR's `D_0 = diag(d_a)` -> `Theta_U`; `g` (vanishing period) -> `varpi_1`; CYC's `h` -> `rho_3`; `f_m` -> `frak f_m`; CYC's `f_i` -> `phi_i`;
CYC's `tau`, `tau_k`, `T_k` -> `sh`, `Upsilon~_k`, `Upsilon_k`; `tau_k(H^b)` -> `psi^k H^b`; `B = prod b_a` -> `Bprod`; `B(f,g)`, `Bmat`, `P_ab` (CYC) ->
`CC(f,g)`, `CCmat`, `JJ_ab`; CYC's `M` -> `M_cyc`; Bergman `B` -> `omega_(0,2)`; `b_k` (Stirling) -> sans `b_k`; `K`, `K_vc`, `K_0` (slice) ->
calligraphic `K`, `K_0`; `K = Q((u))` -> blackboard `K`; `K_q` -> sans `K_q`; `L` (leg class) -> calligraphic `L`; GI's `L = Q((z))` -> written out;
FSL's ring `L` -> `RR_res`; CP's `Lam` -> `\Cstar`; PB's `Lam` -> `Xi`; G5's `Lam` -> `\logq`, `Lam'` -> `-\Xi`; G1's `ell` -> `\logq`; G1's `d_T` ->
`\partial_T`; G1's `W` (Wronskian) -> `\Wr`; G5's `W` -> `\mathsf W`; G5's `NN` -> `\Nring`; `Cf` -> `\Cexp_(alpha,beta)`; G5's `iota` -> `\jmath`;
G5's `A_L` -> `\AL`; TR's `Lambda = -5 kappa` -> the number `-5 kappa`; RD's `Lambda_v` -> `-z + K`; CP's `N` -> `SPn`; `N_a` -> `eps_a`; DP's `M`, `N`
-> sans `M`, `N`; PB's `V` (residue) -> `VV`; PB's `HH` -> `\hat\cV`; theory's `V^(k)` -> `(-1)^(k+1) Pi_(k-1)`; PB's `W`, `G(x)`, `at`, `What`,
`B(X)` -> `WW`, `GG`, `frak a`, `WWhat`, `frak ahat`; `w` (hatted variable) -> `zhat`; `w_i` -> bold `w_i`; `w(k,a)` -> `wt(k,a)`; `F_b(t)` -> sans
`F_b(t)`; `F_h`, `Phi_h` (model) -> `F^m_h`; TR primitive `Phi` -> `Phi_ydx`; `Phi^b_k` -> `phi^b_k`; UR's `sigma = v^(1/5)` -> `varsigma`; the TR
involution `sigma` -> `iota_b`; RD's `s = dT_P - v` -> `delta_P`; the scalar of the affine slice -> `mu`; PB's `s(z)` -> local sans `s(z)`; LMI's
`s*`, `alpha0`, `M0`, `g`, `F`, `Q` -> `s_*`, `a_*`, `M_*`, `phi`, calligraphic `F`, sans `Q`; FSL's `beta` -> `frak b`; LL's `xi`, `gamma` -> `XX`,
`frak r`; `xif_a` -> `xi_a`; GI's `x` (LG variable) -> `p`; S15's `x = s varpi_0` -> `s varpi_0`; S15's and G1's `kappa` (TFT scalar) -> `\varkappa`;
`e` (Euler-dilaton) -> sans `e`; GI's `eta_n`, `theta_n` -> not used; CYC's `Y` (Yukawa) -> `Yuk`; Hessians `H` -> `Hess`; the state space `H = Q^5` ->
`\bbH`; `omega` (TFT) -> `\TFT`; `omega^` -> `\hat\TFT`; `omega^m_v` -> `\TFT^m_v`; PPZ's `T(z)` -> `\ttr_R(z)`; `h_v` (vertex genus) -> `g_v`; the
`n`-th derivative `d^n Fbar` -> `\dd^n Fbar`; `A_e`, `Phi_j` (HPF) -> `\cO_e`, `\Phi^{[j]}`; DD1's `d(z)` (WKB constant) -> `\mathsf k(z)`; `c_10 =
c_5^2/2` -> `[z^10] D_0 = ([z^5] D_0)^2/2`; CP's `Lam1` -> `[5]_L`; CP's 7-variable phase `Phi` -> local `SS_7`; UR's lift `Pi` -> local
calligraphic `Z`; DLMF shift `h` -> `theta`; CP's `P_k(L)` -> local `p_k(L)`; S15's "(U)" -> "(unit)"; FSL's "(U)" -> `lem:fsp-univ`.
Gap notes: G2's `k` (field) -> `\Bbbk`; `V(A, A~)` -> `\Gau(A, A~)`; `S = diag(zeta^-a)` -> `zeta^(-gr)`; `kappa^U_j` -> `chi^U_j`; `beta = u Bprod` ->
a local letter; `e_c`, `h_j`, `phi_j`, `rho_a`, `rho^U_a` -> local; `psi_i` -> `L Phi_i`; G3's `s(z)`, `sigma(q)`, `pr`, `iota = (x + psi')^-1`, `AA` ->
local letters; G4's `lam` -> `vartheta`, `et_j` -> `\bar e_j`, `e_j` -> `e^TR_j`, `S_h` -> `\Sigma_h`, `C_a` -> local; G5's `c(eps)`, `Delta`, `g`,
`psi_eps`, `e(eps)`, `w_a`, `m_a`, `lambda`, `mu`, `sigma`, `M`, `P`, `Pi`, `E`, `rho_n`, `Lambda`, `Lambda'` -> local letters (declared).
Names: "Theorem G" -> (G) `thm:gau-G` or the Gap Theorem `thm:mod-gap`; "Theorem S15", "s.1.5" -> the comparison theorem `thm:flh-main`;
"D(ii)", "D(iii)", "LD" -> `prop:mod-dii`, `thm:mod-lattice`, `prop:mod-LD`; "M1".."M10" -> the statements of s06 by name; "Lemma 0" ->
`lem:tra-frames` (theory) or `lem:gau-const` (CYC); "Lemma E, T, F" -> note names (`lem:tra-E`, `lem:tra-T`, `lem:tra-F`); G4's "Lemma T"
(trace formula) -> a local name; PB's and RL's "Lemma T(b)" -> `lem:fsp-transl` (b); "Lemma K" -> dropped; "Lemma C" -> `lem:mai-coord`;
"Lemma V, Q, 1" (GI) -> `lem:gau-vanmod`, `lem:gau-vanq`, `lem:gau-uniq`; "Lemma 1, 2" (UR) -> `lem:gau-lift`, `lem:gau-scal`; "Lemma 1, 2"
(POLY) -> `lem:con-u0`, `lem:gau-poly`; G1's "Lemma A".."Lemma F" -> `lem:set-grading`, `lem:set-unitfrob`, `lem:set-birkhoff`,
`lem:set-periods`, `lem:set-wronskian`, `lem:set-order3`; G2's "5.1" -> `thm:gau-U`; G3's "FL", "CV" -> `lem:tra-frames`, `lem:tra-conv`; G5's
"Lemma 3.A", "3.B" -> `lem:mai-Psi`, `lem:mai-Ghat`; G5's "FL4 at the generic point" -> `lem:fsp-generic`; "Prop. S" (UR) -> `prop:set-modrows`;
"Prop. S" (GI) -> `prop:gau-symp`; "Corollary F" -> `cor:mai-natfact`; "Theorem D" (UR) -> dropped; "Theorem D" (CP) -> `thm:gam-D`; GAM's
"Theorem A(i)" -> never printed; "Lemma CFU" -> `lem:con-cfu` (the uniformising map); "(A2)", "(H)" (UR's purity and G2's one-pole shape),
"(N)", "E4", "R0..R7", "S1..S9", "H1..H9", "(L0)..(L4)", G4's "(L1)..(L4)", "(1.1)" (G5) (tags) -> never printed; "ML1/ML2" -> never printed;
(Q1)-(Q3) of LL -> [Eynard:2007kz] Thm 5.1, `thm:mod-Lreg`, `thm:mod-gap`. Bibliography: `LPfq` (and G1's "LP18" for arXiv:1803.01409) -> `LP20`.

## 10. Conventions fixed for the whole note

1. **R-matrix:** Givental's `S ~ R e^(u/z)`, `R = 1 + O(z)`, symplectic `R(-z)^T eta R(z) = eta` (flat form); [PPZ15] Def. 2.2 action with legs
   `R^-1(psi)`; translation `\ttr_R(z) = z 1 - z R^-1(z) 1` ([PPZ15] Def. 2.13), also for the dilaton leaves of TR (`prop:mod-Rtr`);
   `Cont_Gamma` includes `|Aut Gamma|^-1`.
2. **`D = q d/dq`**; `D = -(1-u) d/du`; `zD = -(1-u) z d_u`; `D = C_1 \partial_T`.
3. **Colour 0** is the critical point `sigma = L` of the Mellin-Barnes summand and the fixed point `p_0`; colour `i` = `zeta^i L` = `p_i`;
   covariance `Phi_i(zhat, u) = Phi_0(zeta^-i zhat, u)` (`lem:gau-colrows`), `D_i(z) = D_0(zeta^-i z)`, `h_(w_i)(z) = h(zeta^-i z)`.
4. **Hatted versus flat:** flat `(u, z)`, classes `H^a`, rows `Q`; hatted `(u, zhat)`, `zhat = z/s^4`, `Hhat^a = s^(1-nu_a) H^a`, `xhat = varpi_0`.
5. **Model normalisation:** `W_v(p) = p^4/4 + v/p`; LG integrals of `e^(-W_v/z)`; TR curve `x = W_v`, `y = p`; `H = d/dv`; `e = z 1 + v H`; the
   model TFT is DEFINED by `Deltahat_i^(1/2) := zeta^(4i) varsigma`; the TR TFT is `Delta_j = -5 varsigma_j^2` with the uniform orientation
   `dw_j/dp(varsigma_j) = sqrt(-5) varsigma_j`.
6. **Lowered components** for solutions and gauges; vector form only in `lem:tra-E`, `lem:tra-T` (`P^vec = J P J`).
7. **Square roots** `zeta^(1/2) := zeta^3`, `Delta_i^(1/2) = zeta^(4i) C_0/L`, `(W_i^0)^(-1/2) = zeta^i`.
8. **QRR sign:** twist weight `-5 zeta^i`; `h_1 = -3/20`; (QRR) used in the form `Phi_i(z, 1) = h_(w_i)(z)`.
9. **TR kernel:** engine convention (`conv:mod-kernel`); EO's printed kernel = ABDKS's = minus the engine's.
10. **Normalisations:** canonical `P_11(0,z) = 1`; `varpi_1 = u + O(u^2)`; `varpi_0(0) = 1` (conifold frames); `v = u + O(u^2)`; `Fs_0(0) = 1`;
    `ph_0 = 1`; `c_0(L) = 1`; normalised stationary-phase series have constant term 1; `d(0) = 1` for `P_can`.
11. **"Regular":** a model correlator is regular if its `v`-exponent is `>= 0`; a graph term is regular if it lies in `Q[[u]]`; never "zero".
12. **Logarithms and branches:** `log(a/z) := log a - log z`; `Gamma_St(1 + w) := w Gamma_St(w)`; branches are part of the symbols only.
13. **Pairings:** legs in `LL` are paired by `Sum_b Res_(varsigma_b)`; the circle pairing only for legs in `Q[v][p]`; log legs on a ray.
14. **Vertex tensor:** the homogeneous form `d^(2g_v-2+n) \dd^n Fbar_(g_v)(d q_T)`; it is the Taylor expansion at the base point when `d(0) = 1`.

## 11. Macros of `NOTE/macros.tex` (loaded by `main.tex`; in addition to the generic letters `\bbA .. \cZ`, `\dd`, `\ri`, `\re`, `\de`, `\fps`)

`\Dq` (D), `\zh` (\hat z), `\tc` (t_{\mathrm c}), `\Pc` (P_{\mathrm{can}}), `\Dn` (\mathsf D), `\cn` (\mathsf c), `\bk` (\mathsf b), `\bw`
(\mathbf w), `\hw` (h_{\bw}), `\Bprod` (\mathcal B), `\YC` (Y^{\mathrm C}), `\YCQ` (Y^{\mathrm C}_{\mathrm Q}), `\IC` (I^{\mathrm C}), `\Qr`
(\mathcal Q), `\Ghat` (\hat G), `\Rhat` (\hat R), `\Mo` (M_0), `\Fs` (\mathsf F), `\Fm` (F^{\mathrm m}), `\FTR` (F^{\mathrm{TR}}), `\Wv` (W_v),
`\sK` (\mathcal K), `\sKo` (\mathcal K_0), `\fr` (\mathrm{frame}), `\se` (\mathsf e), `\ML{h}` (\mathrm{ML}_{\mathcal K}(h), one argument),
`\cH`, `\cL`, `\cX`, `\cW`, `\cG`, `\cV`, `\cO`, `\cR` (calligraphic), `\fa` (\mathfrak a), `\fc` (\mathfrak f), `\Conc` (\mathcal C), `\Jet` (\mathcal J),
`\Mcyc` (M_{\mathrm{cyc}}), `\sh` (\mathrm{sh}), `\Phican`, `\PhiSP`, `\Phiact`, `\SPn` (\mathrm{SP}^{\circ}), `\SPh` (\widehat{\mathrm{SP}}),
`\Gst` (\Gamma_{\mathrm{St}}), `\Ginf` (G_\infty), `\Cstar` (\mathbb C_\ast), `\wt`, `\Res` (operators), `\Hess`, `\PFfour` (\mathrm{PF}_4),
`\PFfive` (\mathrm{PF}_5), `\Fg{g}` ((\mathrm F_{g}), one argument), and the new ones of the revision: `\Qn` (\mathscr Q), `\logq`
(\mathfrak l), `\Cexp` (\mathscr C), `\Nring` (\mathscr N), `\TFT` (\boldsymbol\omega), `\ttr` (\mathfrak t), `\bbH` (\mathbb H), `\AL`
(\mathbb A_L), `\Gau` (\mathrm{Gau}), `\Wr` (\mathrm{Wr}). Theorem environments as in the template (`thm`, `lem`, `prop`, `cor`, `defn`, `rmk`,
`notation`, `example`, `conv`), numbered within sections; `\draftnote` exactly as in the template. Bibliography `\bibliographystyle{amsalpha}`,
`\bibliography{refs}` (`NOTE/refs.bib`).
