# Gap phase synthesis (key gaps/critic), part 3 of 3: OUTLINE s04-App. C, NOTATION, LABELS, BIBPLAN, errata, questions, re-runs

Continuation of `GAPS_SYNTHESIS_2.md` (s.3 conditionality, s.4 open items, s.5 first half). Same conventions.

## 5 (continued). Changes to the OUTLINE files, s04 to the appendices

**OUTLINE_5.md (s04, s05).**
* `conv:tra-forms`: the two settings of G3 s.2.1 (graded `Q[v]` with Givental's `H((z^-1))`, `HH_+ = H[z]`; u-adic `HH^u = H (x) Q((z))[[u]]`,
  `HH^u_+ = H (x) Q[[u, z]]`, `HH^u_- = H (x) z^-1 Q[z^-1][[u]]`, splitting, `Omega_s`), and G3 s.2.2 (summable families, maps of finite type,
  the summation lemma). The `K_0` condition reads `z^-1 Q[z^-1][[u]]` in the u-adic setting (G3 s.3 item 6).
* `sec:tra-tower`: a lemma for `Theta = 1 + N` and its inverse `Sum_r (-N)^r` in both settings (G3 part 2 Step 1; G3R K8 (b)).
* `lem:tra-frames`: statement G3 part 2 s.1 (a)-(d) with K1; proof Steps 0-5 (Step 2 replaces the source's "top power"); trivial frames
  exactly `k < nu_a + 1` (ten); delete the `[S]` `\draftnote`. `cor:tra-slice` as G3 Cor. (ii).
* `lem:tra-conv`: parts (a)-(d) (G3 part 3 s.1-s.3) and the added (e) vertex resummation, (f) dilaton re-expansion with G3R K3, (g) the
  two-slot frame expansion (G3 part 2 s.3; or in `cor:tra-edges` step 5); K2 (cite `lem:gau-vanshape`), K7, K10 (print the dependencies:
  `def:set-modcohft` homogeneity = `lem:set-modhom`, `lem:set-dilaton`, `lem:gau-diag`, ML_K(g), `lem:tra-T` (d), `v'(0) != 0`); delete the
  "argued" `\draftnote`.
* `lem:tra-E`: statement and proof G3 part 4 s.2 in `AA = Q((x))[[u, psi']]` with `f` in `HH^u_-`; step 0 (divisibility over `Q[[u]]`) also
  serves `def:tra-edge`; K5 for the side remark; delete the O-8 `\draftnote`. `cor:tra-edges` step 3: the line `P^vec J^C = eta^-1 (P Y^C)` (K6).
* `lem:tra-F` step 1 and `def:red-flatsum`: define the vertex tensor by the homogeneous form `d^(2h-2+n) d^n Fbar_h(d q_T)`; it equals the
  Taylor expansion at the base point iff `d(0) = 1` (`lem:tra-conv` (f)), which holds for `P_can` (s.2.3).
* `thm:red-main` steps 8-9: cite `lem:tra-conv` (e) or its finite route. `rmk:tra-K`: the decay hypothesis is `lem:tra-conv` (c).
* The `\draftnote{D1 ..}` on `Pi_i` through `S(z)` is moot under DECISIONS D1 (define from `Y^C`); at most a remark.

**OUTLINE_6.md (s06).**
* R3 everywhere: (8-1) with `n = 4` at `a = (0,0,0,0,4v)`, `f = 4 W_v`, critical points `p^5 = v`, values `5 p_c^4`; R8 (d) renames.
* `thm:mod-tr`: the quoted statement is s.3.1's DOSS entry (G4R R-3); step 3 is the TFT dictionary `Delta_j = -5 varsigma_j^2` from (2-30)
  with the uniform orientation `dw_j/dp(varsigma_j) = sqrt(-5) varsigma_j` (G4 F10), the frozen basis and the factor `(-5)^(h-1)` (s.2.10);
  the transfer `kappa = 4 -> 1` by the new `lem:mod-scale` and the bookkeeping `[R(./c).omega]_(h,n) = Sum_D c^-D [R.omega]^(D)_(h,n)` (G4R R-2).
* `prop:mod-Rtr`: statement G4 part 2 s.1 with R-1 (dilaton leaf `psi (1hat - M_0^-1(psi/(kappa varsigma^4)) 1hat)`), `R_TR(w) =
  M_0(w/(kappa varsigma^4))`; new `lem:mod-laplace` (F1; its "Meaning" sentence as a remark, R-7); proof s.3-s.5 with `lem:set-M0sympl`;
  delete the `\draftnote{D6/G2 ..}` and the "uniqueness of homogeneous R-matrices"/"sign of `kappa`" text (F7).
* `def:mod-pair`: legs in `LL = Q[v][p] + Q(v)[x]` are paired by `Sum_b Res_(varsigma_b)`; the circle `|p| = R` equals it only for legs in
  `Q[v][p]` (F2; `Res_(p=0) x^5 omega_(1,1) = -15/4`); log legs on a ray. Same wording in `thm:mod-Lreg`.
* `lem:mod-xfun`: proof G4 part 3 s.2 (Lemma T, the trace formula for residues; discs disjoint from `J`, R-8); the EO pins (F9) optional;
  delete the `\draftnote{G3 ..}`.
* `thm:mod-Lreg`, `lem:mod-abdks`, `prop:mod-period`, `lem:mod-vanper`: G4 part 4 with R-4 (the printed edge weight (7.7); drop the optional
  bound `K = 6h - 4 + 2n`); Thm 7.1 applied to `x = W_v` (W5 dissolves; XY's `x = z^4 + V/z` only in appendix B with `4^(2h-2+n)`, F6);
  route 2 by the exact identity (F5); delete the `\draftnote{W5 ..}`.
* `thm:mod-V` base `h = 1`: new `lem:mod-w11` and the evaluation (G4 part 5 s.2-s.3) with R-6; cross-checks optional; delete the
  `\draftnote{D5/G1 ..}`.

**OUTLINE_7.md (s07, s08).**
* `lem:flh-dil`: applies `lem:set-dilaton` with `\varkappa = d`, `d(0) = 1` (s.2.3).
* `thm:con-Pa`: statement G5 s.1.2 (i)-(iv), proofs s.1.3-s.1.5, with F1.1 ((1.1) for amplitudes in `Q[sigma, z]`), F1.3 (over `A_L(zeta)`, then
  Galois invariance), F1.4 (local letters); F1.2: move `def:gam-c` and `prop:gam-rational` into `sec:con-form`, or state (i) for the s08
  Mellin-Barnes form and identify it with `c(z; L)` by `lem:con-form` (d) (s.7 Q6); a new appendix-A lemma "FL4 at the generic point" (G5
  s.1.3 Steps 1-2); delete the domain `\draftnote`.
* `prop:con-deg`: (L0)-(L3) (G5 s.4.2) and (L4) (G5R s.4.2) in step 3; step 5 as G5 s.4.3 with "`deg rho_n <= 4n` by (L4)"; step 4 no longer
  needs `lem:gau-poly`; F4.2 (attribution: completes ER s.4, s.7 item 3); F4.3 (write `L`; rename `M`, `P`, `Pi`, `E`, `rho_n`); delete the
  resonance `\draftnote`.

**OUTLINE_8.md (s09, s10).**
* `prop:gam-hatflat`: cites `thm:con-Pa` (iv); CP3 s.3.10 at most a remark.
* `cor:mai-gap`: definitions G5 s.2.2 (`NN`, `Cf`, `F_g^NN` with the division-free edge, `F_g^cf`) with F2.1, F2.2; proof G5 s.2.4 with the
  colour-covariance formula and the `Lambda'`-isometry line (G5R s.2.3); cite PPZ15 Def. 2.13, Prop. 2.14 (Prop. 2.4 for the plain action);
  the two frame remarks optional (scaling: the count of `lem:flh-comm`'s proof, F2.3); delete the `\draftnote`.
* `rmk:mai-constant`: dictionary and differences 1-6 (G5 s.2.6), pointer F2.4; keep `\draftnote{CGYZ26 v2; re-read in the PDF}`.
* `cor:mai-natfact`: `Psi(u) := Qhat(u, 0)`, `R~ := Psi^T eta^-1 Qhat`; Lemmas 3.A (citing the `Delta^(1/2)` convention), 3.B; proof G5 s.3.3
  with the `Lambda'` line and the `u = 0` clause (G5R s.3.2); rename `w_a`, `m_a`; delete the `Psi(u)` `\draftnote`. `rmk:con-why`: only the
  symplecticity of `G` detects `D = 1` (computed, G5R s.3.3; appendix C).
* `thm:mai-main`: the conditionality paragraph of s.3.10.

**OUTLINE_9.md (appendices).** A: "FL4 at the generic point" after `lem:fsp-dq` (and `lem:mod-laplace` if the editor prefers it here). B:
conventions in one place: Givental's convention for `R` and PPZ's for the action; the dilaton-leaf sign (PPZ15 Def. 2.13); kernels (EO =
ABDKS = minus the engine); the orientation `sqrt(-5) varsigma_j`; the TFT normalisations `Delta_j = -5 varsigma_j^2`, `Deltahat_i^(1/2) =
zeta^(4i) varsigma`; XY's normalisation with `4^(2h-2+n)`. C: the phase's VERIFIED computations with ranges (s.8), numbers only there.

## 6. NOTATION, LABELS, BIBPLAN, SOURCE_ERRATA

**NOTATION** (register; local unless stated):
* Input names (M), (G0), (T), (IFT) (never (S1)-(S3)); `\varkappa` for the dilaton/TFT rescaling scalar (R8 (i)); `NN`, `Cf`, the formal
  `log q` (`Lam`; its image `log(1-u)` written `-Xi`, table 8) (R8 (h)).
* s02: in `prop:set-cyclicLR` the derivative in `T` (G1's `d_T`; not s03's `d_t`), the formal `log q` (G1's `ell`; s03's `ell` is `log u`),
  `Mop`, the periods `c_k` and `phi_k`, `Fhat`, `ctilde_m`, `X_a`, `s_a`, `g`, `W`, `w` (each clashes, G1R s.2.3); in `lem:set-M0sympl`,
  `lem:set-modhom`: `P_ij`, `B`, `S_ij`, `Lam`, `Pi`.
* s03: `\mathcal B_a^{(m)}`, `\mathcal Y = A~^-1`, `dlog Theta_U`, `\mathfrak d`, `pr`, `lambda_m`, `kappa^U_j` (rename: `kappa_a` taken), `P^\natural`,
  `P^{[delta]}`, `Rot_i`, `zeta^(gr)` (instead of `T_i`, `S`); renames of G2R C5: the free constant `gamma` (NOTATION's `gamma = 7/10`), `beta =
  u Bprod` (frame parameter), the field `k` (`\Bbbk`), `V(A, A~)` (e.g. `\mathrm{Gau}(A, A~)`), `e_c`, `h_j`, `phi_j`, `rho_a`, `rho^U_a`; `psi_i`
  -> `L Phi_i`.
* s04: `\cH^u`, `\cH^u_+`, `\cH^u_-`, `AA`, `iota = (x + psi')^-1`, `Theta`, `N`, "summable", "of finite type"; the scalar of `def:tra-slice` (R8 (a)).
* s06: the Laplace variable (DNOPS's `zeta`, G4's `lam`), `L_c`, `G_m`, the number `-5 kappa` (G4's `Lambda_kappa`; not the leg map), the frozen
  basis `et_j`, `C_a`, `beta_00`, `alpha`, `y_k` of `lem:mod-w11`; `varsigma_b` for branch points (R8 (d)).
* s08: `A_L = Q[L, L^-1, (L^5 - 1)^-1]`, the embedding `L -> 1/s` (G5's `iota`; clashes with CGYZ26's `iota` in s10), `PF5^L`; G5's `c(eps)`,
  `Delta`, `g`, `psi_eps`, `e(eps)` (F1.4); in `prop:con-deg`: `L`, and new letters for `M`, `P`, `Pi`, `E`, `rho_n` (F4.3).
* s10: `F_g^NN`, `F_g^cf`, `Psi(u)`, `R~`, `Ghat`, `G` (table 8 reserves some); G5's `w_a`, `m_a`, `lambda`, `mu`, `sigma` (rename).
* NOTATION_2 s.10: add that `Deltahat_i^(1/2) := zeta^(4i) varsigma` is a definition; the residue pairing for `L`-legs; PPZ's `T` for dilaton leaves.

**LABELS** (register; prefix as shown): s02 `prop:set-cyclicLR`, its corollary (G1 proposes `cor:set-cyclicPF5`), `lem:set-M0sympl` (one label
for G1 s.5.1 = G4 s.5; drop G1's `lem:set-M0symp`), `lem:set-modhom`; s03 `lem:gau-lead`, `lem:gau-degree`, `lem:gau-colmod`, `lem:gau-colrows`,
`cor:gau-rowsNN`; s04 a label for the `Theta` lemma (e.g. `lem:tra-theta`) and, if separate, the two-slot lemma (e.g. `lem:tra-twoslot`);
s06 `lem:mod-laplace`, `lem:mod-scale`, `lem:mod-w11`; appendix A one label for "FL4 at the generic point" (e.g. `lem:fsp-generic`). LABELS
l.12 and l.32: (S1) -> (M) (R2). No existing label changes meaning.

**BIBPLAN.**
* MR1408320 (Givental, alg-geom/9603021v2; root): from "optional, s01" to cited in s02: s.9 set-up (p.23-27), Prop. 11.4, Thm 11.6; Cor. 11.7
  only for its misprint; eprint and version note; `\draftnote{IMRN numbering not compared}`.
* DuninBarkowski:2012bw (DOSS, arXiv:1211.4021; root): from optional to cited (s06, (Q)): s.2.3 (graph weights, vertex weight, dilaton-leaf
  weight), Thm 4.1. Dunin-Barkowski:2015caa: s.2.3, s.2.6 (2-30)-(2-35), Prop. 4.3 (1), the sentence before Thm 8.1, (8-1) at `a = (0,0,0,0,4v)`;
  Thm 8.1 only as context; Prop. 5.2 dropped.
* New key for Arbarello-Cornalba, alg-geom/9406008v1 (e.g. `AC96`): formula (1.7), `kappa_0 = 2g - 2 + n` (p.5); journal data to confirm.
* PPZ15: add Def. 0.3 (iii) and the convention sentence after Def. 2.2 (p.20); Prop. 2.14 also in s02 and s10.
* MR2917177: add 8.2 Proposition, (8.6), 8.7 Proposition (i), the Riemann-Roch identity of p.40, 8.10 Proposition (context only); the use is
  G1 s.3.2 (not "8.4 and 8.7 (ii)" alone).
* GJR18: add s.4.1 (the twisted CohFT) and s.5.2 eq. (14) (QDE) for (G0); s.5.1 (13) only as corroboration (misprinted).
* LPfq (1803.01409v2): add (10), and s.2.4.3 (17), (19)-(20) (corroboration); s.3.2 (26) also in the (QRR) derivation. G1's "LP18" is LPfq.
* MR2276766: Thm 2 only as the modern formulation; Cors. 5, 6 non-equivariant on `X`; no Cor. 5.1 exists.
* MR1728879, MR717614: unchanged roles; to be read for O1 before (QRR) is promoted.
* ABDKS25: add (7.6), (7.7) (with `1/hbar^2`), (7.10). EO: add Thm 4.5 (4-25), App. A (A-1), (A-5), (A-25) (optional pins). CGYZ26: add s.4.2
  Prop.-Def. 4.5 (1) (to confirm in the PDF). Hock25: optional.
* Optional, not read: Guo-Ross 1703.06955v1 (Compositio Math. 155 (2019); Props. 6.1-6.2, corroboration; new key); LLY97 alg-geom/9712011;
  Witten 1991 (for `1/24`, only if the genus-1 log term is printed); R1's candidate MR1653024 (not needed with the pin on MR1408320).

**SOURCE_ERRATA** (entries for the main session; the file is the editor's). External: (1) Givental alg-geom/9603021v2, Cor. 11.7 (2):
`l_alpha` for `lambda_alpha` (Thm 11.6 is right); (2) DNOPS 1509.06954v2 s.2.3: dilaton-leaf sign (PPZ15 Def. 2.13 and DOSS s.2.3 are right;
journal version not checked); (3) GJR18 (13): `tau-dot = H + q d/dq(I_1/I_0)` (intended `C_1 H *`). Project sources: (4) XYR l.37-38: the ABDKS
kernel equals EO's printed kernel (both minus the engine's; conclusion unchanged); (5) RD1 s.1.3, DC1 D9: ten trivial frames, not eight (DC1
D10 is right); (6) TRF s.2: `Q((x^-1))[[psi1]]` is right only in the hatted setting; the flat Lemma E needs `Q((x))[[u, psi']]` and `f` in
`HH^u_-`; (7) PS s.4.1, DB D14, DD1 D1.8: `D_i` in `Q(zeta)[[z]]`, in `Q[[z]]` only for `i = 0`; (8) OU3 s.3: `Psi^T eta^-1 Psi = 1`, not `Psi eta^-1
Psi^T = eta^-1`; OU3 s.6d: "orthogonal" -> isometry; (9) PB s.2.2 "on `0 < u < 1`" -> the comparison in `Q((s))`; CP3 s.3.10's "as rational
functions" is circular inside P(a); (10) ER s.4 step 5: "W" undefined; the step needs (L1) and (L4); (11) PLAN_2026-11-03: CG has no Cor.
5.1; (12) MS5 s.2: "Teleman Prop. 8.5(ii) (s.8.2)" is the ar5iv label of v3 8.7 Proposition (ii); (13) the "O-8" of `lem:tra-E` is DIGEST_C1_6
s.7 O-8, not STAGE2_SYNTHESIS_2 s.5 O8. Gap-note corrections to apply while transcribing (not errata): G1 s.1.2 and s.3.2 Step 5; G2 s.9.4
Step 2; G3 part 2 s.1 (c), part 4 s.3 Rmk 2; G4 s.1, part 2 s.1 (b), s.4, s.6 (iii), s.7, part 4 s.1; G5 (1.1), s.2.2, s.2.4, s.2.5.

## 7. Questions for the main session

* **Q1 (QRR).** Keep it quoted (recommended; the conditionality is then exactly s.3, and G1's Steps 1-5 appear as a remark with one
  `\draftnote` on O1), or promote it to a proposition: then read Mumford's relation and a non-vanishing statement (O1, s.2.2) and cite the
  convention identity (columns of the CohFT R-matrix solve the QDE: GJR18 (14) with Lemma 5.2; sign fixed by `h_1 = ph_1 = -3/20`; G1R s.4.6).
* **Q2 (IFT).** Drop it by using `lem:fsp-ift` in LMI step 2 (recommended).
* **Q3.** One label for the symplecticity of `M_0` (`lem:set-M0sympl` recommended).
* **Q4 (D4).** Route 2 of the Gap Theorem is complete with ABDKS alone; making it the proof of record removes [Hock25] from (Q).
* **Q5.** Adopt the homogeneous definition of the vertex tensor in `def:red-flatsum` (G3-2f).
* **Q6.** Placement for `thm:con-Pa` (i) (G5R F1.2).
* **Q7.** A second source for (M) (LLY97 or an exposition; to be read first).

## 8. Re-runs and own checks (evidence only; nothing above rests on them except where a sign or a misprint was decided)

* **Re-runs** (copies in the scratch folder, `py -I -B`; `critic/rerun_record.log`): `g1_cyclic_check.py 30 12`, `ref_g1_loc.py`,
  `ref_g1_model.py`, `ref_wkb_run.py 6 10`, `r3_frames.py 60`, `g4_w11.py`, `g4_period.py 5`, `rf4_check_b.py`, `ref5_item4.py 40`,
  `ref5_item1.py 12 12`: every log is identical to the stored one (two stored G5 referee logs also contain timing lines). What they decide:
  `kappa_a` to `q^30` with four failing controls; (M) at order `q` and the `+600` of the literal Cor. 11.7; `<H,H>_(0,2,1) = 2875 = 5 [q] kappa_1`,
  `<H^3,H^4>_(0,2,1) = -3250 = 5 [q] kappa_3`; `M_0` symplectic to `t^24`; the colour covariance of row 0 with the `h_(w_0)` control; ten trivial
  frames; `-1/60` per branch point; `S_h`; the dilaton-leaf sign by TR at `(1,1)`; (L0)-(L4) with both obstructing controls; `thm:con-Pa`
  (a)-(c) with the control `g = q`.
* **Own checks** (`critic/critic_checks.py`, log `critic/critic_checks.log`; exact; under one second): A, the rewriting of Givental's Thm 11.6
  (controls: non-CY degrees, literal Cor. 11.7); B, the trivial frames from the support rule of `Pi_0` alone (control `k <= nu_a + 1`); C,
  `[z^j] PF5^(u_i) = zeta^(-ij) [z^j] PF5^(u_0)` and the transport operator `5 zeta^(4i) D` on `L^m` over `Q[c]/(c^5 - 1)` (control PF5 `+ q z^2`);
  D, (L0)-(L4) by applying `P` to monomials (controls `P + z^2 x^3`, `P + z x^6`). All pass; every control fails.
* Not re-run: `g2_unique.py`, `ref_unique_run.py`, `rf4_check_a.py` (minutes each; re-run by their referees), the G3 Lemma E scripts and the
  G5 provers' scripts (re-run identically by G3R and G5R).

## 9. Files and reading record

Written (all new): `NOTE/gaps/GAPS_SYNTHESIS.md`, `GAPS_SYNTHESIS_2.md`, `GAPS_SYNTHESIS_3.md`; `NOTE/gaps/critic/critic_checks.py`,
`critic_checks.log`, `rerun_record.log`. The scratch folder holds copies only. Read: every file of `NOTE/gaps/G1..G5` named in s.0 of the
prover notes and referee reports; `design/OUTLINE_REVIEW.md`, `design/DECISIONS.md`; `OUTLINE.md`; `OUTLINE_3.md` items 1-11; `OUTLINE_4.md`
`sec:gau-systems`-`sec:gau-R0`; the planned `\draftnote`s of `OUTLINE_3`-`_9` (by search); `BIBPLAN.md` part 1; the label list of `LABELS(_2)`;
`NOTATION_2.md` s.9 (search); `stage2/flat/S15_PROOF.md` l.37. No literature was opened by the critic: every citation status above rests on the
readings recorded by the provers and referees, which I checked for consistency against each other and against the design files.
