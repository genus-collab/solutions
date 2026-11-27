# Outline, part 8: s09 the Gamma identity; s10 the main theorem assembled, universality, remarks

## s09_gamma — target 5 pp

Purpose: prove Theorem Gamma, the value `c(z; 1) = h(z)` of the full-Stirling saddle series of `def:gam-c` (defined in s08), by local
Mellin inversion; transfer it to the hatted series; conclude `D = 1` and (F_g) for all `g`; by-products at `q = oo` (`Ginf`, the free
coefficients) and remarks.

Covers: DD2 s.1 D1-D9, s.2 GAM-3..GAM-18, s.3A, s.4, s.6, s.7 items 1-12; DD1 prop:pol-Ginf, rmk:glu-K8, thm:glu-D; CP s.0-s.1; CP3
s.3.1-s.3.11; CP4 s.4; RF, RF2 (RF-C1..C6); RA s.1, s.3; LF (LF-C1..C6); LA (LA-C1..C7); GEN s.3.2; LIT2 s.3; ST5 s.1-s.3; ST5b K1, K2, K6-K10,
O1-O6, O14, O18; RS1 A3, A7-A10; RS2 items 29-42; G5 s.1.2, s.1.6. Subsections `sec:gam-series`, `sec:gam-lmi`, `sec:gam-thm`, `sec:gam-rmk`. The
proof of record of Theorem Gamma is LMI; the 7-variable chart is a remark.

### sec:gam-series
* Recall `def:gam-c` and `prop:gam-rational` (s08, `sec:con-form`; DEC Q6) and the second Mellin-Barnes form `lem:con-form` (d); nothing new
  is defined here. `rmk:gam-motivation`: dropped (D8).

### sec:gam-lmi — local Mellin inversion
* `lem:gam-deltasrc` (lem; GAM-3, "the delta lemma with a source"): phase `P(s) + a_* s m`, `P = Sum_(k>=2) p_k s^k` (exactly this form), amplitude
  `e^(c m) A(s, z)`: `0` is non-degenerate, `det = -a_*^2`, and exactly `N = exp(P(delta)/z) A(delta, z)`, `delta = -c z/a_*`. Proof DD2 s.3A steps
  1-5 [W] (RA s.3.2; LF s.2; LA s.1.4); the moment formula of step 3 with "`= 0` for `r < a`" (LF s.2).
* `def:gam-lmidata` (defn; DD2 D7): `s_*` in `C`, `a_*` in `C^x` (source `alpha0`), `M_*` in `Z` (source `M0`), `phi` holomorphic near `s_*` (source
  `g`), `phi_k := phi^(k)(s_*)/k!`, `FF` in `O_(s_*)[[z]]` (source `F`) with `FF(s_*, 0) != 0`; `a(sigma) := a_*(sigma - s_*)` EXACTLY linear; a
  punctured disc `U_*` with `a_* + (L - s_*) phi''(L) != 0` (LF-C3 = LA-C3); `lambda(L) := -phi'(L) - a_* log a(L)`; `S := lambda sigma + phi + a log
  a - a`; `A := a^(M_* - 1/2) Poch_(M_*)(z/a) exp(Sum_k bk_k (z/a)^(2k-1)) FF`; `c^FF(z; L)` the normalised series at `L`; `delta := -M_* z/a_*`, `phit(x)
  := phi(s_* + x) - phi(s_*) - phi'(s_*) x`, chart parameter `\mathsf Q := e^(-lambda/a_*) = a(L) e^(phi'(L)/a_*)`.
* `thm:gam-lmi` (thm; GAM-4, Theorem LMI): (i) each `c^FF_n(L)`, `L` in `U_*`, depends only on `L`; (ii) it extends holomorphically to `s_*`;
  (iii) `c^FF(z; s_*) = exp(phit(delta)/z) FF(s_* + delta, z)/FF(s_*, 0)`, re-expanded in `z`. Scope (LA-C6): `phi`, `FF` holomorphic at `s_*` and
  free of `lambda` (holomorphic dependence on `\mathsf Q` allowed, the value is that of the `\mathsf Q = 0` data); `a` exactly linear; `FF(s_*, 0) !=
  0` is a genuine hypothesis. Proof DD2 s.3A steps 0-3 [W] (RA s.3.3 with LF s.1, s.3.2 and LA s.1.1-s.1.4): step 1 the two-variable formal
  integrand (LA-C1), `lem:fsp-cov` (translation, `chi = exp`), `lem:fsp-stirling` (`a = a(sigma)`, `M = M_*`, measure `dx`), `lem:fsp-fubini`;
  step 2 the chart `nu = mu - lambda/a_*` (rename `nu` locally: it clashes with `nu_a`), `\mathsf Q` single-valued (LF s.3.2), the family of
  critical points by the formal implicit function theorem `lem:fsp-ift` (DEC Q2; LA s.2.1), the identification of critical points (LF-C1:
  written out), `\mathsf Q^(M_*)` drops out of the normalisation (LF-C2), `lem:fsp-param` (i); step 3 `lem:gam-deltasrc` with `P = phit`, `c = M_*`.
  Full input list LA-C2 (including the `A[[z]]`-linearity of `SPn`, `def:fsp-wick`).
* `rmk:gam-lmi-words` (rmk): the saddle series continued to the colour point is the normalised value of the regular part at the first pole
  `sigma = s_* + delta` of `Gamma(M_* + a/z)` (RA s.3.1).

### sec:gam-thm — Theorem Gamma and the value of `D`
* `thm:gam-main` (thm; GAM-5, Theorem Gamma): (i) `c_n` in `\AL`; (ii) no `c_n` has a pole at `L = 1` (`q = 0`) and `c(z; 1) = h(z)`; (iii)
  `c_n(zeta L) = zeta^-n c_n(L)`, so no pole at the `zeta^i`, `c_n` in `Q[L, L^-1]`, exponents `= 4n mod 5`, `c(z; zeta^i) = h_(w_i)(z)`; (iv) no pole
  at `L = 0` and `c(z; 0) = Ginf(z)`, the series of [GJR18, Lemma 7.7, (37)] (`thm:gam-ginf`; REV_2 M4), so `c_n` in `Q[L]` (`cor:gam-poly`).
  Proofs: (ii) LMI at `s_* = 1`, `a_* = -1`, `M_* = 0` [W] LF s.4.1, LA s.1.5, with `FF = (5 sigma)^(1/2) prod_(j>=1) (zeta^j - sigma)^(-1/2)
  exp(Sum_k z^(2k-1) bk_k [(5 sigma)^(1-2k) + Sum_(j>=1) (zeta^j - sigma)^(1-2k)])`, the regular part of `A_St` (REV_2 M10; LF-C4 = LA-C4); value
  `FF(1, z)/FF(1, 0) = exp(Sum_k bk_k T_k z^(2k-1)) = h(z)`, `T_k = 5^(1-2k) + Sum_(j>=1)(zeta^j - 1)^(1-2k)`; (iii) [W] CP3 s.3.9, RF2 s.4 (or
  LMI at `zeta^i`). Only (ii) is used by the chain. Exactly one Gamma factor degenerates at each of `L = 1, zeta^i, 0` (LA s.1.5). Numbering
  differs from CP's (DD2 s.7 item 1). `cor:gam-colours`: dropped (D8; `cor:gau-colours` gives the colours on the gauge side).
* `prop:gam-hatflat` (prop; GAM-7): `c_n(L) = L^(4n) ph^SP_n(L^-5)` as rational functions of `L` (this is `thm:con-Pa` (iv); G5 s.1.2, s.1.6);
  hence `ph^SP_n(1) = c_n(1)` and, with `thm:gam-main` (ii), `Phi^SP(zhat; 1) = h(zhat)`. Proof: cite `thm:con-Pa` (iv); the separate route of
  CP3 s.3.10 (with RF-C2) is not printed. No circularity: `thm:con-Pa` does not use this proposition.
* `cor:gam-equiv` (cor): `Phi^SP(z; 1) = h(z)` iff `c(z; 1) = h(z)`; with `cor:con-A2`: `D = 1` iff (Gamma).
* `thm:gam-D` (thm; Theorem C of s01; DD1 thm:glu-D): `D_0 = 1` in `Q[[z]]`, `D = 1`, and (F_g) holds for every `g >= 2`. Proof: `thm:con-glue`,
  `thm:con-Pb`, `prop:gam-hatflat`, `thm:gam-main` (ii), `cor:gau-colours` [W] ST5 s.2 S1-S7; RS2 s.6.2. Inputs (REV_2 M3 with DEC Q2): (QRR)
  through `thm:con-glue`, appendix A (with `lem:fsp-ift` through `thm:gam-lmi`), and s03, s08. Correction ST5b O8 (CP's longer route through
  Lemma N is not used).
* `rmk:gam-K8` (rmk; DD1 rmk:glu-K8): given `thm:con-Pa`, regularity at `L = 1` is automatic; the content of (Gamma) is the value; LMI proves
  the regularity independently (PB P8; ST5b K8).
* `cor:gam-poly` (cor; GAM-8): `c_n(L) = Sum_k a_(n,k) L^(4n-5k)`, `4n - 5k >= 0`; with `prop:gam-hatflat`, `[u^(4m)] ph^SP_(5m) = c_(5m)(0)`. Proof [W] LF
  s.4.3. Correction RF-C3/RS2 item 38 (polynomiality of `ph^SP` is `thm:con-Pa` (iv), not the chart).
* `thm:gam-ginf` (thm; GAM-9): `c(z; 0) = Ginf(z) = exp(Sum_(l>=1) (-1)^(l+1) B_(5l+1)(1/5) z^(5l)/(l(5l+1)))`; LMI at `s_* = 0`, `a_* = 5`, `M_* = 1`,
  `delta = -z/5`, `\mathsf Q = (-q)^(-1/5)`; normalisation of LF-C5. Proof DD2 s.3A steps (1)-(5) [W] LF s.4.2; LA s.1.5, s.2.3, using
  `lem:fsp-dlmf5118` (formal; LA-C5) and `B_k(theta+1) - B_k(theta) = k theta^(k-1)` [DLMF 24.4.1, from 24.2.3].
* `cor:gam-free` (cor; GAM-10): `frak f^SP_m = c_(5m)(0) = [z^(5m)] Ginf`; with `thm:con-Pb`, `Sum_m frak f^can_m z^(5m) = Ginf(z)` (no (QRR) needed);
  `frak f^act = frak f^can` by `thm:gam-D`; `[z^5] Ginf = B_6(1/5)/6 = 4537/3937500`. Proof [W] RA s.3.4; LF s.4.3. Prior art: [GJR18] Lemma 7.7,
  (37) (`m = 5`, `k = 1`) proves the formula for their Landau-Ginzburg series; the note does NOT claim that series equals `c(z; L)` (LIT2 s.3.6).

### sec:gam-rmk — remarks (at most about one page; D8)
* `rmk:gam-sevenvar` (rmk; GAM-11; kept under D8, started as statement plus sketch, REV_2 M19): `Omega_7`, critical point `(L, 5L, zeta^j -
  L)`, `det H_7 = 1/(L^2(1 - L^5)^2)`, two orders of integration, the 5-variable chart with `1/det H_ch = -L^2 [5]_L^2` (RF-C1), the product of
  Euler integrals at `L = 1`; warning (FSL3 s.8): the 7-variable stationary phase is not defined at `L = 1`, the value is read off the chart.
  Not used in any proof: `\draftnote{sketch only; the authors decide whether to include the full argument (then with lem:fsp-delta)}`.
* `rmk:gam-blind` (rmk; GAM-14, trimmed): (Gamma) does not see the measure of the twist factor (`Gamma(5s)` gives the same value at `L = 1`);
  the measure `dy` enters through `thm:con-Pb` and `thm:con-Pa` (iv); at `L = 0` the twist factor is singular and `M_* = 1` matters. The
  unwritten evaluation at `L = 0` with `M_* = 0` is not printed (its planned `\draftnote` goes).
* `rmk:gam-necessity` (rmk; GAM-15): (a) with proof (GEN s.3.2 Lemma 4 [W], an analytic statement about actual critical points); (b) the
  one-Gamma identity `c^FF(z; L) = C(z/a(L))`, `C = 1` (LF s.1 [W]); (c) the falsifier cases only in appendix C (rule 3).
* `rmk:gam-Mcomplex` (rmk; GAM-13; D8: one sentence): LMI for `M_*` in `C` with `exp(Sum_n beta_n(M_*) (z/a)^n)` in place of the Pochhammer factor;
  not used; `\draftnote{not proved here}`. `lem:fsp-L3C` is dropped.
* `rmk:gam-general` (rmk; GAM-16, short): LMI applies to any Mellin-Barnes summand with exactly one degenerate factor with a linear argument
  and `lambda`-free regular data; the examples only in appendix C; the (QRR) of those theories is not quoted: no claim about them.
* `rmk:gam-prior` (rmk; GAM-17; REV R4): [MR1901075] s.10 (untwisted `CP^n`, the chart `x_0 = Q/(x_1..x_n)`, 10.4 Lemma; continuity in `Q`
  asserted); [GJR18] Rem. 7.2 states the analogous identity for their Landau-Ginzburg integral without proof (display L39), Lemma 7.7 (37);
  [MR2553377] Lemma 4.20 (exact expansion, at `z = -1`); [MR2276766] s.7 (11). `\draftnote{describe [MR1901075, s.10] neutrally (outline or
  proof)}`. "We are not aware of a proof of (Gamma) in the literature" (DL s.2.3 item 6). `rmk:gam-Iexp`: dropped (D8).

Length 5 pp.

## s10_main — target 4 pp

Purpose: assemble Theorem A and Corollary B, record the consequences of `D = 1` for the natural form, prove universality in flat form,
and collect the remarks on graph-locality and on the hypotheses.

Covers: MS5 s.1-s.2; ST2 s.2 steps 5-8 and the Corollary; DA A9-A12, A14-A17, A20-A21, A25-A27; DB S33-S36; DC1 UNI-1..UNI-9; DC2 s.2.5-
s.2.6; UR s.4-s.6; RN s.3; NF s.1-s.2, s.6; BR s.1, s.3; BV s.1; RD2 s.3.3-s.3.4; UN1, UN2, URF; DL s.2.1, s.7 items 4, 9; G5_2 s.2, G5_3
s.3 with G5R s.2, G5R_2 s.3. Subsections `sec:mai-thm`, `sec:mai-nat`, `sec:mai-univ`, `sec:mai-rmk`.

### sec:mai-thm
* `thm:mai-main` (thm; Theorem A): OUTLINE s.0.2 (REV R6 form). Proof: `thm:gam-D` gives (F_g); `thm:flh-main` gives `T_Gamma = F_Gamma`;
  `cor:red-can` with `thm:mod-MLK` gives `F_Gamma` in `Q[[u]]` for every edge graph and `F_vertex = c_g t_c^(2-2g) + Q[[u]]` [W] ST2 s.2 step 8;
  MS5 s.2. Close with the conditionality sentence of OUTLINE s.0.2 and the remark that the statement is per graph (`rmk:mai-graphlocal`).
* `def:mai-Fcf` (defn; new; G5_2 s.2.2 [G] with G5R F2.2): `F_g^\Nring := Sum_Gamma Cont_Gamma[E_(R~), q_(R~); \TFT]` in `\Nring (x) Q(zeta)` (no
  further `|Aut Gamma|^-1`), `R~ := Psi^-1 Q = Psi^T eta^-1 Q`, `Psi := Q(z = 0)`, `Psi^T eta^-1 Psi = 1` in `\Nring`; the edge written without
  division, `E_X(psi, psi') := (N_X(psi, psi') - N_X(-psi', psi'))/(psi + psi')`, `N_X := 1 - X^-1(psi) X^-1(psi')^T` (PPZ's edge for symplectic `X`,
  [PPZ15] Rem. 2.1); the Galois group of `Q(zeta)` permutes the colours (print the covariance formula, G5R s.2.3), so `F_g^\Nring` lies in
  `\Nring`; `ev(F_g^\Nring) = F_g` by (T) with `thm:set-qrr`; `F_g^cf := \Cexp_(alpha,beta)(F_g^\Nring)` in `Q((s))`, the conifold expansion of `F_g`.
* `cor:mai-gap` (cor; Corollary B; DEC D3): in every conifold frame `F_g^cf = c_g t_c^(2-2g) + Q[[u]]` (since `Q[[u]] = Q[[t_c]]`: no fractional
  power, no `t_c^k`, `3 - 2g <= k <= -1`). Proof G5_2 s.2.4 [G] with G5R s.2.3: Step 1 `\Cexp(R~)` is symplectic (from `\Cexp(Q) = P_can Q^m`,
  `thm:gau-R0` (iv) with `cor:gau-colours` and `thm:gam-D`; `prop:gau-symp`; `Q^m(-z)^T eta^-1 Q^m(z) = 1`, which needs the line `Lambda'^-1 eta^-1
  Lambda'^-1 = eta^-1` (`nu_a + nu_(3-a) = 2`) and `lem:set-M0sympl`; then `Psi Psi^T = eta`); Step 2 `lem:flh-frame` (no extra `|Aut Gamma|^-1`,
  F2.2); Steps 3-4 Theorem A and the regrouping `R'.(M.\TFT) = (R' M).\TFT` by [PPZ15] Def. 2.13 and Prop. 2.14 (Prop. 2.4 for the plain
  action; REV_2 M2); finite sum. Local letters for G5's `Lambda`, `Lambda'` (not the leg map). The complex-frames remark of G5_2 s.2.5 is
  dropped (D8); the scaling count (with F2.3: the count of `lem:flh-comm`'s proof) goes into `rmk:mai-constant`. Delete the `\draftnote`.
* `rmk:mai-constant` (rmk; D10): `c_g` is the constant of [CGYZ26, Conj. 1.5] in the same normalisation of the vanishing period, and the two
  conifold expansions agree on the period ring: G5_2 s.2.6 [G] (dictionary, differences 1-6), pointer of G5R F2.4 ([CGYZ26] Ex. 1.3; s.4.2,
  Prop.-Def. 4.5 (1)); the scaling count of G5_2 s.2.5 explains "the same normalisation"; `\draftnote{numbers from arXiv:2609.13921v2, read in
  the HTML rendering; to be re-read in the PDF}`. Do not import the sign of [CGYZ26] Def. 1.8. Values `c_2..c_5`: `rmk:set-cg`.

### sec:mai-nat — consequences of `D = 1` for the natural form
* `cor:mai-natfact` (cor; DB S35 (a), (d), S18 (c); DA A10, A11): statement G5_3 s.3.3 [G]: assume `D = 1` (`thm:gam-D`); in every conifold
  frame (a) `Qhat(u, zhat) = Ghat(u, zhat) Rhat^mod(tau(u) zhat)`, `Ghat` a power series in `u` whose coefficients are polynomials in `zhat` (degree
  `<= floor((5n+4)/4)` at `u^n`); (b) `Qhat(0, zhat) = Rhat^mod(zhat)`: the conifold symbol is the model's, so `M_0 = R~(zhat; 0)` (NF's definition)
  and `Psi(0) = Psi_0`; (c) `R~(zhat; u) = G(u, zhat) M_0(tau(u) zhat)` with `M_0` on the RIGHT, `G = Psi(u)^-1 Ghat Psi_0`, `G(zhat)^T G(-zhat) = 1`, the
  same `M_0(tau zhat)` in every frame. Definitions (G5_3 s.3.2): `Psi(u) := Qhat(u, zhat = 0)`, `R~ := Psi^-1 Qhat = Psi^T eta^-1 Qhat`. Lemmas:
  `lem:mai-Psi` (G5_3 s.3.2 Lemma 3.A: `Psi(u)_(aj) = zeta^(j(a+1)) w_a(u)` (rename `w_a`, `m_a`), `Psi(0) = Psi_0`, `Psi^T eta^-1 Psi = 1`; cite the
  `Delta^(1/2)` convention for `Q = zeta^j X`, G5R_2 s.3.2) and `lem:mai-Ghat` (Lemma 3.B: `Ghat := Lambda P_can(u, s^4 zhat) Lambda'^-1` with
  local diagonal matrices, entries `tau^((1-nu_b)/4) Sum_k zhat^k u^(mu(b,k)) p_(k,b)(u)`, `mu(b,k) = (4k + nu_b - nu_(b+k))/5 >= max(0, (4k-4)/5)`).
  Proof G5_3 s.3.3 [G] with the `Lambda'` line and the `u = 0` clause (G5R_2 s.3.2); inputs `thm:gau-R0`, `cor:gau-colours`, `thm:gam-D`,
  `lem:gau-R1`, `lem:gau-diag`, `lem:set-canon`, `prop:gau-symp`, `lem:set-M0sympl`, `prop:set-cyclicLR`. (a), (b) are not tests of `D = 1`
  (G5R_2 s.3.3). "`M_0` on the left fails" is computed only (appendix C). Delete the `Psi(u)` `\draftnote`.
* `cor:mai-actual` (cor; DB S36): the actual gauge is the canonical one, `P^act = P_can` (`c^act = 1`; `P^act Y^C(v(u)) = Y^C_Q`). Proof [W] GI s.6
  and ST2b s.6.1 with `thm:gam-D`.
* `rmk:mai-apriori` (rmk; DA A9; UR s.6.1): every outer edge and leaf is `O(u)`; a genus-`g` outer graph with `E` edges is `O(u^(E+2-2g))`; only
  `G_1, .., G_(2g-3)` enter. [W] UR s.6.1; NF s.6. Not a step of the proof.
* `rmk:mai-tau` (rmk; kept under D8, labelled unproved): uniqueness of the rescaling `tau` (DA A12; DB S35 (b)): "not used; not proved here",
  `\draftnote{not proved here; the authors decide whether to keep this remark}`; the genus-by-genus selection of `tau` is evidence (appendix C).
* `rmk:mai-vertexterm` (rmk; DA A25): the leafless vertex is `c_g (u tau^(-3/2)/xhat)^(2-2g)`; the dilaton leaves convert it to `t_c` (`rmk:red-leaf`);
  `t_c` enters only as `Pi_0` of the leaf point. Dropped (D8): `rmk:mai-degree`, `rmk:mai-G1`; the Euler-dilaton remark only in appendix C.

### sec:mai-univ — universality; sec:mai-rmk — remarks
* `lem:mai-pole` (lem; DC2 lem:uni-pole): a regular gauge to the pulled-back model needs `prod b = -v'^5/v`; a pole forces `ord_0 Bprod = -1` and
  `v'(0)^4 = -res_0 Bprod`. Proof [W] BR s.1; BV s.1 (its converse paragraph is not used).
* `lem:mai-coord` (lem; "Lemma C" of BR s.3; DC1 UNI-2): a change of `Bprod` is a change of coordinate `u = phi(w)`, `phi = v_X^-1 o v_q`; (a)-(e)
  as DC2 `lem:uni-coord`. Proof [W] BR s.3. The hatted-encoding count (`lam^(10(g-1))`) is one sentence.
* `cor:mai-univ` (cor; DC2 `cor:uni-flat`): for a self-dual one-pole system `A` with the model's product and its canonical gauge, `r = -1/lambda^2`,
  `lambda = f(0)|g(0)|` (S01 C10): every edge graph of the flat sum is regular and the vertex graph is `c_g (T_A/lambda)^(2-2g) +` regular, `T_A =
  Int_0^u b_0`; B-perturbed systems by `lem:mai-coord`. Proof [W] RD2 s.3.3, TRF s.2, with `thm:red-main`, `thm:mod-MLK`. Flat form only (the
  hatted-engine form is NOT claimed). Correction BVRF 3 / ST2b K6 (`|c0| = m^10/|f(0) g(0)|`).
* `rmk:mai-annih` (rmk; DC1 UNI-3): `thm:red-main` applies to annihilating non-gauge dressings (flat form); no hatted statement.
* `rmk:mai-stored` (rmk; DC1 UNI-4): the gauge of the other normalisation, `P_can(1 + m z^5 + ..)`, is regular, symplectic only through `z^9`,
  violating (A); flatness is not sufficient by construction. Proof of the algebra [W] UN1 s.2.3; URF s.2.1; numbers in appendix C.
* `rmk:mai-graphlocal` (rmk; DC1 UNI-7; REV_2 M6): the computations of appendix C show that in a non-canonical normalisation six genus-4 edge
  graphs become singular while the total is unchanged; so no argument about total potentials can prove Theorem A; the proof is graph-local
  and uses the odd part of `log c` fixed by (G).
* `rmk:mai-necessity` (rmk; DC1 UNI-8): none of the hypotheses (A), flat regularity can be dropped in general (examples in appendix C); for (Sp)
  no example exists (without it the edge is not a power series); flatness in `u` is not a hypothesis. `lem:mai-delta6`: dropped (D8).
* `rmk:mai-analogy` (rmk; DA A20; short): the factorisation has the shape of the Hodge factorisation for the resolved conifold; the counterpart
  of the Hodge factor is `M_0(tau zhat)` together with the regular gauge (NF s.2; NFS s.4: "natural" is not "forced by the graph sum", DA A19).

Length 4 pp.
