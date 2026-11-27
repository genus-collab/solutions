# Outline, part 9: appendices A (formal stationary phase), B (signs and normalisations), C (evidence); page budget; check-list

## a_formalsp — appendix A, formal stationary phase — target 7 pp

Purpose: state and prove, over any commutative `Q`-algebra and purely formally, every stationary-phase lemma used in s08-s09, including
Stirling's series as the stationary phase of Euler's integral, the formal DLMF 5.11.8 and the formal implicit function theorem.

Covers: DD2 s.1 D2, D5, s.2B (SP-*), s.3B, s.4 (normalisations, phases, Stirling, branches); FSL s.0-s.10 (base text, as ST5b s.8 item 2
(vii) prescribes), CP2 s.2.1-s.2.8 (shorter proofs where named), RL s.1-s.3 with RL-C1..C8, LF s.4.2, LA s.2.3, RS1 A3, A7; G5 s.1.3.
Rule: descriptive names only (no "FL2", no "(U)"; PB's FL1-FL7 and FSL's FL1-FL4 never appear). Every lemma is proved. The non-formal
inputs of the appendix are named in its first paragraph (REV_2 M11 (a)): the exact Gaussian integral inside the alternative proof of the
perturbation lemma (the proof of record, CP2 s.2.2, is algebraic), and the branch part of `lem:fsp-refl` (used only in `lem:con-form` (d)).
Subsections `sec:fsp-def`, `sec:fsp-basic`, `sec:fsp-change`, `sec:fsp-stirling`, `sec:fsp-res`.

* `def:fsp-ring` (defn; SP-DEF): the weighted completed ring `R_n(A)` (`wt(t^alpha z^j) = |alpha|/2 + j`, finitely many terms of weight `<= W`),
  the Gaussian expectation `<t^alpha>_C` and `E_(zC)` (FSL s.1.2-s.1.3); `A` any commutative `Q`-algebra with `H` in `GL_n(A)`, `C = -H^-1`.
  `lem:fsp-ring` (FSL Lemma 1.1) [W FSL s.1.2].
* `def:fsp-wick` (defn): the normalised stationary-phase functional `SPn_c[F; a] := E_(zC)[exp(F_3/z) a(c + t)]` (FSL Def. 1.2 = CP's `N`,
  RL-C7), `A[[z]]`-linear in the amplitude (LA-C2 = RS2 item 35, an input of LMI; REV_2 M11 (b)); the normalised series `SPh` (divided by its
  constant term `a(c)`); the full symbol `(2 pi z)^(n/2) det(-H)^(-1/2) e^(F(c)/z) SPn` (FSL Def. 1.3; branches as part of the datum). Phase
  conventions `e^(F/z)` versus `e^(-W/zhat)` (FSL s.1.4).
* `lem:fsp-univ` (lem; SP-U; FSL's "Lemma U", renamed): `[z^k] SPn` is a polynomial over `Q` in finitely many Taylor coefficients (phase
  orders `3..2k+2`, `A_m` orders `<= 2(k-m)`, the sharp CP (F1) bound) and the entries of `C`; ring maps and parameters pass through [W CP s.1.2
  (F1); FSL s.1.5]. `lem:fsp-const` (SP-CONST, trivial).
* `lem:fsp-homog` (lem; SP-H, RL's Lemma H, missing in FSL: RL-C1): `SPn[mu F, a(., z)](z) = SPn[F, a(., mu z')](z')`, `z' = z/mu`, with the symbol
  rule `det(-mu H)^(1/2) := mu^(n/2) det(-H)^(1/2)` [W RL s.3.2; CP s.1.2 (F3)].
* `lem:fsp-ibp` (lem; SP-IBP): exact forms have zero stationary phase, including convergent infinite sums of exact forms (used by
  `thm:con-Pa`) [W CP2 s.2.1; FSL2 s.3.1].
* `lem:fsp-pg` (lem; SP-PG): the Gaussian perturbation identity and its infinitesimal form; the case `K_1 = K(y) - H_xx` depending on outer
  variables (RL-C5) [W CP2 s.2.2, algebraic proof of record; FSL s.2 (G5), (G5') alternative]; RL-C4 (wording of the ideal).
* `lem:fsp-prod` (lem; SP-PROD) [W CP2 s.2.3; FSL (G4)].
* `lem:fsp-ift` (lem; SP-IFT, formal implicit function theorem) [W FSL2 s.3.2]; RL-C4. It replaces the holomorphic implicit function theorem
  in LMI step 2 (DEC Q2) and is used in `thm:con-Pa` (G5 s.1.3 Step 1).
* `lem:fsp-transl` (lem; SP-T, translation (a) `p` in `(qA)^n`, (b) `p` in `(zA[[z]])^n`) [W FSL2 s.3.2]. PB's and RL's "Lemma T(b)" is (b) (REV_2
  M11 (c)); never "Lemma T" (that name is `lem:tra-T`).
* `cor:fsp-movebase` (cor; SP-T', moving the base point) [W FSL2 s.3.2; CP2 s.2.6 (iii)].
* `lem:fsp-param` (lem; SP-PAR, parameters (i)-(iv); (iv) the substitution `eps -> z`, RL-C4) [W FSL2 s.4; CP2 s.2.6].
* `lem:fsp-dq` (lem; SP-DQ: `d/dq` of the stationary phase; the iterated clause `P(z D_q) SP = SP[P(sigma) ..]` when `D_q F = sigma`, needing
  `D_q a = 0`, RL-C2) [W FSL2 s.4.5; CP2 s.2.6 (iii')].
* `lem:fsp-generic` (lem; new; G5 s.1.3 Steps 1-2 [G]): the stationary phase at a generic base point: for the phase `S + eps sigma` with
  `S''(L)` a unit of the differential ring `\AL`, the critical point given by `lem:fsp-ift` is `exp(eps D)(L)`, and `D` acts on normalised saddle
  series through the Taylor data at the moving critical point; used by `thm:con-Pa` (i).
* `lem:fsp-cov` (lem; SP-COV, change of variables, analytic, formal or parameter-dependent) with `lem:fsp-piola` (Piola identity) [W FSL2 s.5;
  CP2 s.2.4]; RL-C3 (`P_s := G o psi_s - G(0) - y^T H y/2`).
* `lem:fsp-cfu` (lem; SP-CFU: the change of variables over `Q[[u]]` for `Y = X psi(X^5; u)` at `X = 1`, also at `X = zeta^-i`) [W FSL2 s.5.4].
* `lem:fsp-fubini` (lem; SP-FUB: iterated stationary phase with the factor `frak b(w) = (det H_vv(w)/det H_vv(w_0))^(-1/2)`, the Schur complement,
  symbols) [W CP2 s.2.5; FSL3 s.6 with the "Rings" paragraph]; RL-C3, RL-C5.
* `lem:fsp-delta` (lem; SP-L2) [W CP2 s.2.7; FSL3 s.7]: written only if `rmk:gam-sevenvar` carries its full proof (it is its only user).
* `def:fsp-stirling` (defn; DD2 D2): `bk_k`, `Gamma_St(Y) := exp((Y - 1/2) log Y - Y + (1/2) log 2 pi + Sum_k bk_k Y^(1-2k))`, `log(a/z) := log a -
  log z`, `Gamma_St(1 + w) := w Gamma_St(w)`, `Poch_M`, `St_m(w)`, `beta_n(theta) := (-1)^(n+1) B_(n+1)(theta)/(n(n+1))` (DLMF's shift `h` -> `theta`).
* `lem:fsp-stirfe` (lem; SP-FE: `Gamma_St(Y + 1) = Y Gamma_St(Y)` formally; `St_(m+1) = (1 + mw) St_m`) [W FSL4 s.9.1 (iv), s.9.2].
* `lem:fsp-stirling` (lem; SP-L3 (i)-(v): Stirling's series as the formal stationary phase of Euler's integral with `x^m`, both measures) [W FSL4
  s.9.3, Proof A, purely formal, proof of record (ST5b K2)]; `rmk:fsp-classical` (rmk; one sentence, not used): Proof B through DLMF 5.11.1,
  2.3.8 and Laplace's method (CP2 s.2.8; FSL4 s.9.4); RF-C6 (the Pochhammer part is formal).
* `lem:fsp-dlmf5118` (lem; SP-5118, the formal DLMF 5.11.8): `log Gamma_St(w + theta) - [(w + theta - 1/2) log w - w + (1/2) log 2 pi] = Sum_n
  beta_n(theta) w^-n` in `Q[theta][[1/w]]` [W LF s.4.2 = LA s.2.3], from `lem:fsp-stirfe` and [DLMF] 24.4.1 (from 24.2.3). Only `thm:gam-ginf`.
* `lem:fsp-refl` (lem; SP-REFL, formal reflection) [W RL s.3.3]; the branch part is not formal (used in `lem:con-form` (d)).
* `lem:fsp-res` (lem; SP-RES: the residue ring `RR_res` (FSL's `L`), `Res := [X^-1]`; closure (FSL Lemma 10.1), exact forms have zero residue,
  invariance under `X -> X + eps`) [W FSL4 s.10.1-s.10.3]; PS-C2. The relabelling is `lem:con-relab` (s08).
* `rmk:fsp-analytic` (rmk; SP-ANAL; one sentence): for real data with negative definite Hessian the full symbol is the Laplace expansion; used
  in no proof; NO citation of Hoermander; the optional context citations are dropped (D8). `lem:fsp-L3C`: dropped (D8).

Length 7 pp.

## b_signs — appendix B, sign conventions and normalisations — target 3 pp

Purpose: one place for every sign and normalisation dictionary, so that the body can refer to it. Proves nothing; every entry is a
convention or a dictionary with its derivation in the body. Covers: DA s.4 C1-C20; DB s.4.1-s.4.2; DC1 s.4 CV1-CV17; DC2 s.4 C1-C12; DD1
s.4.1-s.4.12; DD2 s.4; DL s.4; OUTLINE_2 s.0.4; GS_3 s.5 ("OUTLINE_9"), G4 s.1, G4_2 s.7.
* `sec:sgn-R`: Givental's convention `S ~ R e^(u/z)`, `R = 1 + O(z)`; [PPZ15] Def. 2.2 uses `R^-1` on legs and states that Givental and Teleman
  use the same convention (p.20); DNOPS Remark 2.1 (Givental's convention is the transpose of theirs); the translation `\ttr_R` (Def. 2.13 /
  Prop. 2.12) and the dilaton-leaf sign `+z(1 - R^-1(z) 1)` (DNOPS s.2.3 prints the opposite sign; DOSS s.2.3 agrees with PPZ; `prop:mod-Rtr`);
  `Cont_Gamma` includes `|Aut Gamma|^-1`; the colour-basis edge `E_ij` and the leaves `Dil_i(d)` (REV_2 M8; DA D9, C6-C7), obtained from the
  naturality of `def:set-graphsum`; lowered versus vector form (`conv:tra-forms`).
* `sec:sgn-var`: `D = q d/dq = -(1-u) d/du`; `z d_u`; `\partial_T = C_1^-1 D`; Euler operators `theta_u`, `theta_zhat`; `zhat = L^4 z`; colour `i` = critical
  point `zeta^i L` = fixed point `p_i`; the three points; hatted versus flat (`def:set-hat`, `lem:flh-hat`); `xhat = varpi_0 = 1 + alpha u + ..`.
* `sec:sgn-qrr`: the twist weight `-5 zeta^i` (minus the fibre weight for an Euler-class twist), `h_w`, `h_1 = -3/20`; tangent versus cotangent
  characters ([Giv01] Thm 1 versus [MR1901075] 9.1); Bernoulli conventions (only `B_2k` and `B_k(theta)` enter); [ZZ08] and [LP20] normalise
  asymptotic series by `Phi_s(0) = delta_(0,s)`, which is not (QRR).
* `sec:sgn-model` (`tab:sgn-model`): the model normalisations: `W_v` here; DNOPS (8-1) at `a = (0,0,0,0,4v)`, `f = 4 W_v` (DEC D13); the curve
  `x = z^4 + V/z` of [ABDKS25] (7.2) and XY, with the factor `4^(2h-2+n)` (G4 F6); the LG form `x^4 + u_LG/x`, `u_LG = 4u`, `z_LG = -4z`,
  `eps = -4 zeta^-i zhat`; the TR TFT `Delta_j = -5 varsigma_j^2` with the uniform orientation `dw_j/dp(varsigma_j) = sqrt(-5) varsigma_j` and the
  frozen basis; the model TFT `Deltahat_i^(1/2) = zeta^(4i) varsigma` (a definition); the factor `(-5)^(1-h)`; Hurwitz-space names (DNOPS (8-1)
  `n = 4`; [CvdLPS23] `H_(0;1,4)`); `H = d/dv`, `e = z 1 + v H`, `<H>_1 = -1/12`, `<e>_1 = 1/8`.
* `sec:sgn-tr` (`tab:sgn-tr`): the kernel dictionary of `conv:mod-kernel` (EO's printed kernel = ABDKS's = minus the engine's); `F_h =
  -F^print_h`; this `omega_(0,3)` versus EO (4-20); EO's dilaton equation (4-30) in the engine form `(2 - 2g - n)`; `Lambda^print = -1`; the
  formal Laplace reading of DNOPS (2-31), (2-32) (`lem:mod-laplace`).
* `sec:sgn-legs`: LD (`varrho_b`) as the dictionary of record; the literature-style dictionary (`tau_0(1) <-> -z`, step `+Int dx`) is not used;
  string sign `+d/dx`; the residue pairing for `L`-legs (`def:mod-pair`).
* `sec:sgn-tower`: `Pi_i` only; three equivalent forms of the weight exponent (RD1 s.1.2; MSG1 s.1; TRM s.3); "regular" never means "zero".
* `sec:sgn-norm`: the canonical normalisation `P_11(0,z) = 1`; the other normalisation used in appendix C, defined by `p_(5q,0)(0) = 0` (`q >= 1`),
  `P = P_can(1 + m_5 z^5 + ..)`; the value of `m_5` only in appendix C (REV_2 M7).
* `sec:sgn-sp`: phase conventions; measures `dy` and `dx/x`; normalised versus full stationary phase; the branch rules (DD2 s.4).
* `sec:sgn-lit` (`tab:sgn-lit`): conventions of the literature (DL s.4): LG variables of [GJR18] s.7 (`x_i = L - lambda_i`, measure `dx/x`) versus
  the charts here (`x_j = zeta^j - L`); [MR2553377] Lemma 4.20 at `z = -1`; gap constants ([HKQ09] (1.2); [CGYZ26] Conj. 1.5 with `(-1)^(g-1)`,
  while [CGYZ26] Def. 1.8 has `(-1)^g`: do not import it; [Br25] (1.5)); [MR1408320] Cor. 11.7 prints `l_alpha` for `lambda_alpha` (Thm 11.6 is used).

## c_evidence — appendix C, evidence — target 5 pp

Purpose: report the finite computations that test the statements of the note (ranges, fields, controls, caveats); it proves nothing and
nothing in the body depends on it. Covers: DA s.5 V1-V13, s.7 items 23, 25; DB s.5 E1-E12; DC1 s.5; DC2 s.5 E1-E20, W6; DD1 s.5 E1-E16; DD2
s.5 E1a-E14; DL s.5; the evidence sections of the gap notes and referee reports (G1_3 s.2.12, G1R s.1.2, s.2.2, G1R_2 s.6, s.9; G2_6 s.10,
G2R_3 s.4; G3_4 s.4, G3R_2 s.12; G4_2 s.8, G4_3 s.4, G4_4 s.5, G4_5 s.5, G4R s.1.4, G4R_2 s.3; G5_4 s.5, G5R s.1, G5R_2 s.3.3, s.5; GS_3 s.8).
Rules (rule 3, rule 5): content, ranges, fields (`Q`, `Q(zeta_5)`, `F_p` with the primes written as numbers) and controls only; no script,
log, file or folder names; the other normalisation named by its condition; caveats stated plainly; nothing here is cited by a proof.
* `sec:evi-calib`: calibration of the graph sum against localisation (DA V1; FQ2 s.2: `N_(2,0) = -19/9600`, `N_(1,1) = -305/6`, `N_(1,2) =
  -535775/12`, `N_(2,1) = 2129/5760`, `N_(2,2) = 203401/144`, `F_3(0) = 1/224000`; [LP20] (26) and the weight-`+5` control; the true quintic's
  `F_2(0) = -5/144`); the genus-0 data: (M) at order `q` at generic weights and at the zeta point, the literal reading of Cor. 11.7 failing by
  `+600`; `kappa^LR` to `q^30` with four failing controls; the degree-1 invariants `2875 = 5 [q] kappa_1`, `-3250 = 5 [q] kappa_3`; the Yukawa check.
* `sec:evi-graphs`: genus 2-5 directly (DA V2-V7): graph counts, both frames `(0,0)`, `(1/3,-2/7)`, primes `2^89-1`, `2^61-1` (`2^521-1` once),
  totals and polar parts, content accounting, the controls, the independent second implementation; caveats of V7; `c_g` and the gap of the
  total do not depend on `M`, so they are not evidence for the factorisation (DA A26, V13).
* `sec:evi-hodge`: the linear Hodge frames fail at genus 2 (DA A18, V4), with the genus-2 tables V2-V3.
* `sec:evi-normal`: `D = 1` and (Gamma) in finite windows (DD1 E1-E16; DD2 E1a-E7; DA V10), `Phi^SP = Phi^can`, the residue identity,
  `lem:con-S2` (REV R2), the uniformising map, the free coefficients (`m <= 25` exact, `m <= 100` two primes), the degree equality (`n != 2`),
  the construction of Theorem U, the finite route; `thm:con-Pa` (i)-(iv) at five points, also at `u = 32` and `u = -1/32`, with the control
  `g = q`, and at twelve random points mod `2^61 - 1`; the degree bound and its two obstructing controls (`n <= 30`, `n <= 40`); the checks of
  the natural factorisation (`M_0` on the left fails; only the symplecticity of `G` detects `D = 1`); failing controls; independence caveats.
* `sec:evi-gauge`: (CYC) in many frames (DB E1), (R1)-(R3), (G), symplecticity, Lemmas V, Q, the diagonal, the windows actual = canonical,
  model rows, degrees and the factorisation windows (DB E2-E12; DA V8-V9); Theorem U built along its proof (frames `(0,0)`, `(1/3,-2/7)`,
  `z^0..z^15`), the values `pi(z)`, `p_(1,b)(0)`; the uniqueness clause mod `2^61 - 1` including non-graded gauges, and the rank-5 control
  without the pole; `M_0` symplectic to `t^24` with a non-self-adjoint control; the colour covariance of the model (to `t^10`) and of the actual
  row 0 (`z^0..z^6`, `q^0..q^8`, over `Q(zeta)`) with the controls `h_(w_0)` and PF5 `+ q z^2`; the colour-`i` transport operator on `L^m`.
* `sec:evi-model`: TR = frozen theory, `F_2..F_5`, (1.3) of [Hock25] for `g <= 12`, `v -> 0` limits, closed formulas, legs (DC1 EV-TR1..EV-XY2;
  DC2 E5-E11); `R_TR = M_0` exactly over `Q(zeta_5)` through order 10 (`kappa = 1, 4, -1`; `varsigma = 1, 2`); the dilaton-leaf sign decided by
  TR at `(1,1)` (four points, `kappa = 1, 4`; the printed opposite sign fails); `-1/60` per branch point; the period identity (`h <= 5`) and
  `\Sigma_h` (`h <= 15`); the lattice theorem to weight 400 (exact) and 1000 (two primes), "16200 raw-leg kernel vectors"; the Model Lemma lists
  (326 genuine cancellations); the off-`K` counterexamples of `rmk:mod-offK` and the response `Delta_g(5)`.
* `sec:evi-transport`: Lemma E on random symplectic `R`, Lemma T and Cor. E' on the canonical gauge, `delta_P = u^2/30 + (127/3150) u^3 + ..` in
  both frames, Theorem R at genus 2, the comparison at genus 3 (108 + 4 non-zero coefficients) and its controls; the ten trivial frames
  (`k < 40`; control `k <= nu_a + 1`).
* `sec:evi-univ`: universality (perturbed systems, `Bprod`-variation, `lem:mai-coord` numerically), the non-canonical normalisation (six
  genus-4 graphs singular, total unchanged, responses per unit `m_5`; the value of `m_5`), the necessity examples (DC1 EV-U1..EV-U6; DC2
  E16-E18; DA A17).
* `sec:evi-lmi`: LMI on general data, the falsifiers of `rmk:gam-necessity` (c), the examples of `rmk:gam-general`, the formal lemmas on
  random data (DD2 E8-E13).

## Page budget (target 60-85; this design: 84)
s01 5 | s02 10 | s03 10 | s04 6 | s05 3 | s06 12 | s07 4 | s08 10 | s09 5 | s10 4 | A 7 | B 3 | C 5. Remarks at most about two pages in all
(D8). If the total exceeds 85 pp, drop in this order: `sec:con-finite`, `rmk:gam-sevenvar`, `rmk:gam-general`, `rmk:mai-analogy`,
`rmk:gau-Uremarks`.

## Writer check-list (every section)
1. Read OUTLINE s.0 (rules and the file plan), the section's part of this outline, `NOTATION.md` (+ `_2`), `LABELS.md` (+ `_2`), the digest
   items named, then the sources and gap notes named at each step (the gap note is the text of record for the steps it closes; apply the
   corrections of its referee report named next to it, and the adjudications of GS s.2).
2. Apply the corrections listed at each statement; record each applied correction and each `\draftnote` in `DRAFT_NOTES.md` (section,
   label, source tag): the tags never enter the `.tex` file.
3. Use only the labels of `LABELS.md` (+ `_2`); new auxiliary labels with the section prefix; refer to other sections only through them.
4. No finite computation in a proof; numbers only as definitions, illustrations or (appendix C) evidence.
5. State hypotheses with scope (rule 7); state the conditionality where an input enters (rule 8).
6. Use the literature exactly as `BIBPLAN.md` pins it; do not cite unread statements; use root keys.
7. No step of the chain carries a `\draftnote`; the planned `\draftnote`s are those named in this outline (all outside the chain, OUTLINE
   s.0.2); each says in mathematical words what is missing and names no project file (rule 5).
