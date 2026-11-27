# Outline, part 2: s.0.3 sources and gap notes, s.0.4 conflicts settled by the design, s.0.5 decisions

## 0.3 Source abbreviations (paths relative to HVQ = `formalquintic/hodge_vertex_quintic/`; AG = HVQ/`allgenus/`; NOTE = AG/`note/`)

| abbr. | file | abbr. | file |
|---|---|---|---|
| NF | HVQ/`NATURAL_FORM.md` | FQ2 | HVQ/`FORMAL_QUINTIC_genus2.md` |
| OU2, OU3 | HVQ/`OU_EXTENSION_genus2.md`, `_genus3.md` | CS | HVQ/`CONIFOLD_SYMBOL_genus2.md` |
| NFS | HVQ/`natural_form/SYNTHESIS.md` | UR | HVQ/`natural_form/theory/UNIFORM_REDUCTION.md` |
| RN | HVQ/`natural_form/referee/REFEREE_NAT.md` | LG | HVQ/`natural_form/lg/LG_MODEL.md` |
| MEC, S01 | AG/`MECHANISM.md`, AG/`SYNTHESIS_STAGE01.md` | MS2, MS4, MS5 | AG/`MECHANISM_STAGE2.md`, `_STAGE4.md`, `_STAGE5.md` |
| GI, GH | AG/`gauge/GAUGE_IDENTITY.md`, `GAUGE_IDENTITY_HAT.md` | GR | AG/`gauge/referee/REFEREE.md` (gauge referee) |
| RD1, RD2, RD3 | AG/`theory/REDUCTION_DRAFT.md`, `_2`, `_3` | TRF | AG/`theory/referee/REFEREE.md` |
| MSG1, MSG2 | AG/`model/MODEL_SINGULAR.md`, `_2` | MRF | AG/`model/referee/REFEREE.md` |
| UN1, UN2 | AG/`univ/UNIVERSALITY.md`, `UNIVERSALITY_E3_E5.md` | URF | AG/`univ/referee/REFEREE.md` |
| TRM, TRRF | AG/`tr/TR_MODEL.md`, `tr/referee/REFEREE.md` | XY, XYR | AG/`xyswap/XYSWAP_MODEL.md`, `xyswap/referee/REFEREE.md` |
| ST2, ST2b | AG/`stage2/STAGE2_SYNTHESIS.md`, `_2` | DP, DR | AG/`stage2/dict/DICT_PROOF.md`, `dict/referee/REFEREE.md` |
| LL, LL2, FB | AG/`stage2/logleg/LOGLEG_PROOF.md`, `_2`, `FALSIFIER_B.md` | LR | AG/`stage2/logleg/referee/REFEREE.md` |
| S15, FR | AG/`stage2/flat/S15_PROOF.md`, `flat/referee/REFEREE.md` | BV, BVT, BR, BR4 | AG/`stage2/bvar/BVAR.md`, `BV_TABLES.md`, `referee/REFEREE.md`, `REFEREE_S4.md` |
| ST4, ST4b | AG/`stage4/STAGE4_SYNTHESIS.md`, `_2` | EVID, ER | AG/`stage4/evid/EVID.md`, `evid/referee/REFEREE.md` |
| POLY, PR | AG/`stage4/poly/POLY.md`, `poly/referee/REFEREE.md` | GAM, GP, GMR | AG/`stage4/gamma/GAMMA.md`, `GAMMA_PROOFS.md`, `gamma/referee/REFEREE.md` |
| ST5, ST5b | AG/`stage5/STAGE5_SYNTHESIS.md`, `_2` | CYC | AG/`stage5/critic/CYC_PROOF.md` |
| CF, CF2, CF3 | AG/`stage5/referees/cyc_formal/REFEREE.md`, `_2`, `_3` | CA, CA2 | AG/`stage5/referees/cyc_adversary/REFEREE.md`, `_2` |
| RS1, RS2, RS3 | AG/`stage5/referees/REFEREES_SYNTHESIS.md`, `_2`, `_3` | PB, PB2 | AG/`stage5/pb/PB_AUDIT.md`, `_2` |
| PS, PS2 | AG/`stage5/referees/pb_steps/REFEREE.md`, `_2` | CP, CP2, CP3, CP4 | AG/`stage5/chart/CHART_PROOF.md`, `_2`, `_3`, `_4` |
| RF, RF2 | AG/`stage5/chart/referee_formal/REFEREE.md`, `_2` | RA, RA2 | AG/`stage5/chart/referee_adversary/REFEREE.md`, `_2` |
| LF, LF2 | AG/`stage5/referees/lmi_formal/REFEREE.md`, `_2` | LA, LA2 | AG/`stage5/referees/lmi_adversary/REFEREE.md`, `_2` |
| FSL, FSL2..FSL5 | AG/`stage5/lemmas/FORMAL_SP_LEMMAS.md`, `_2`..`_5` | RL, RL2 | AG/`stage5/lemmas/referee/REFEREE.md`, `_2` |
| GEN, GEN2 | AG/`stage5/general/GENERAL.md`, `_2` | LIT, LIT2, LIT3 | AG/`stage5/lit/LIT.md`, `_2`, `_3` |

Two clashes in the digests' abbreviations: DIGEST_D1 writes GR for the gamma referee (here GMR) while DIGEST_B writes GR for the gauge
referee (here GR); DIGEST_D1 and DIGEST_D2 write RS for `REFEREES_SYNTHESIS.md` (here RS1), so "RS A1-A10" means RS1 s.2.

**Gap notes** (final form of every step they close; pointer = piece + section, e.g. `G2_5 s.9.4`). The pieces and their sections:

| note | pieces (NOTE/`gaps/..`) and sections | referee report (sections) |
|---|---|---|
| G1 (s02) | `G1/GAPS_G1.md` s.0-s.1 ((M)); `_2` s.2.0-s.2.7 (`prop:set-cyclicLR`: statement, inputs, Lemmas A-F); `_3` s.2.8-s.2.12 (proof, corollaries, literature, evidence); `_4` s.3 (`thm:set-qrr`), s.4 (`lem:set-dilaton`); `_5` s.5 (`M_0`, model CohFT), s.6-s.7 | `G1/referee/REFEREE.md` s.0-s.2 (G1R); `_2` s.3-s.10 (G1R_2) |
| G2 (s03) | `G2/GAPS_G2.md` s.0-s.2; `_2` s.3-s.4 (commutant, lifting, scalar lemma); `_3` s.5 (Theorem U); `_4` s.6-s.8 (graded gauges, `lem:gau-diag`, `thm:gau-G` Cor. (a)); `_5` s.9 (colours); `_6` s.10-s.13 | `G2/referee/REFEREE.md` s.0-s.2.8 (G2R); `_2` s.2.9-s.3 (G2R_2); `_3` s.4-s.7 (G2R_3) |
| G3 (s04) | `G3/GAPS_G3.md` s.0-s.3 (conventions); `_2` (frame lemma; two-slot lemma); `_3` (`lem:tra-conv`); `_4` (`lem:tra-E`) | `G3/referee/REFEREE.md` s.0-s.5 (G3R); `_2` s.6-s.14 (G3R_2; s.13 corrections K1-K10) |
| G4 (s06) | `G4/GAPS_G4.md` s.0-s.2; `_2` (`prop:mod-Rtr`, `lem:mod-laplace`, `lem:set-M0sympl`, `lem:mod-scale`); `_3` (`lem:mod-xfun`); `_4` (`thm:mod-Lreg`, second proof of the Gap Theorem); `_5` (`thm:mod-V` base `h = 1`, `lem:mod-w11`; findings F1-F10) | `G4/referee/REFEREE.md` s.0-s.2 (G4R); `_2` s.1-s.6 (G4R_2; s.5 corrections R-1..R-8) |
| G5 (s08, s10) | `G5/GAPS_G5.md` s.0-s.1 (`thm:con-Pa`); `_2` s.2 (`cor:mai-gap`, the constant); `_3` s.3 (`cor:mai-natfact`), s.4 (`prop:con-deg`); `_4` s.5-s.8 | `G5/referee/REFEREE.md` s.0-s.2 (G5R); `_2` s.3-s.7 (G5R_2; s.6 consolidated) |
| GS | `GAPS_SYNTHESIS.md` (s.0 verdict, s.1 status, s.2 adjudications); `_2` (s.3 conditionality, s.4 open, s.5 OUTLINE changes); `_3` (s.5-s.6 changes, s.7 questions, s.8 checks) | - |

**Supersession order** (final form = the last item): MEC s.1, s.3, s.4 -> MS2 s.1 -> MS4 s.2 -> MS5 s.2 -> the gap notes; inside each
stage the synthesis (ST2b, ST4b, ST5b, RS1-RS3) overrides the worker notes; inside the gap phase the referee report overrides the gap
note where it corrects it, and GS s.2 decides every dispute (always for the referee). NF's ARGUED/CONJECTURAL labels and NFS's
PROVED-mod((H), ..) are superseded by MS5 s.2 (the dependence on (H) and on the canonical normalisation disappears through `D = 1`).

## 0.4 Conflicts settled by the design (record; apply silently in the note; listed in `DRAFT_NOTES.md`)

1. **Eynard-Orantin numbers** (v4 PDF of math-ph/0702045): Thm 5.1 closes s.5.1 with (5-14) (`g + k > 1`) and (5-15) (`g >= 1`); the
   temperatures are in s.5.3, (5-27)-(5-30); (5-37)-(5-38) are the pole-moduli variations of s.5.3; Thm 4.4 is (4-23)-(4-24); Thm 4.5 is
   (4-25); Cor. 4.1 is (4-26)-(4-29), proved in App. A with (A-25); Thm 4.7 is (4-30); kernel (4-3), recursion (4-9), `F^(g)` (4-10).
2. **PPZ numbering** (v3 PDF of 1303.1043): Prop. 2.12 (`T_a = z[R(1) - 1]`, `T_b = z[1 - R^-1(1)]`), Def. 2.13 (`R.Omega = R T Omega`,
   `T(z) = z 1 - z R^-1(z) 1`; in the note `\ttr_R(z)`), Prop. 2.14 (group action). The leaf translation of the note is Def. 2.13's.
   PPZ state that they use the inverse of the R-matrix in Def. 2.2 and that Givental and Teleman use the same convention (p.20).
3. **Teleman numbering**: v3 PDF (Theorem (2) in (1.15); 8.2 Proposition; 8.4 Proposition; (8.6); 8.7 Proposition (i)-(ii); 8.10
   Proposition), never the ar5iv labels ("Prop. 8.5(ii) (s.8.2)" is v3's 8.7 (ii)).
4. **(QRR) source**: GJR18 (17) (ST5b K9, O11); the sign `-5 h_alpha` is read in the authors' LaTeX (DL L27). math/0108100 Thm 9.1 is the
   untwisted statement and 1702.06096 (5) the definition of the theory.
5. **Vanishing tower**: the note uses only `Pi_i`, defined from `Y^C` (DEC D1); theory's `V^(k) = (-1)^(k+1) Pi_(k-1)`.
6. **TR kernel orientation**: engine convention = DNOPS (2-28) = Hock (2.1); the printed EO kernel (4-3), (4-9) EQUALS the ABDKS kernel
   (1.1) with `omega^(0)_(1,0) = -y dx`, and both are minus the engine's: `omega_(h,n) = (-1)^n W^EO_(h,n) = (-1)^n omega^ABDKS_(h,n)`,
   `F_h = -F^print_h`, `Lambda^print = -1` (ST2b K1; LR s.2; G4 s.1; G4_5 s.6 F4). Appendix B states it once.
7. **Leg dictionary of record**: LD, not XY's (ST2b K10); string sign `+d/dx` (S01 C8); Euler-dilaton edge leg `-(16/45) z^9 + ..` (S01 C7).
8. **(CYC) scope**: any power-series `y_0` not proportional to `varpi_1` (RS1 A1, CA-K2 over CF-E5); the frame labels relabelled (RS1 A2).
   **Placement**: (CYC) opens s03 (a fair reading of "a set-up lemma", MS5 s.5 item 2; REV_2 M18); do not move it into s02.
9. **The other normalisation of the computations**: `P = P_can (1 + m_5 z^5 + ..)`, defined by `p_(5q,0)(0) = 0`; its value of `m_5` only
   in appendix C (REV_2 M7).
10. **`H = d/dv`, Euler-dilaton vector `e = z 1 + v H`** (S01 C3); `<H>_1 = -1/12`, `<e>_1 = +1/8` (values in appendix B only).
11. **`lambda` in the universality corollary**: `r = -1/lambda^2`, `lambda = f(0)|g(0)|` (S01 C10); `|c0| = m^10/|f(0) g(0)|` (ST2b K6).
12. **Exponentials in R0**: the same `e^(u_i/z)` is stripped on both sides (PS-C8); UR's `e^((u_i - u_i^c)/z)` is not used (on amplitudes the
    two describe the same identity, G2R_2 s.2.9). PB s.2.1's `e^(-What(1)/zhat) = e^(u_0/z)` holds only up to `e^(u_0^c/z)`; s08 expresses
    Theorem P's normalisation through `Phi^SP` (normalised series), never through the exponential.
13. **Rows**: lowered components `Q_a^(i)(z) = eta(H^a, R(z) e~_i)`; `Phi_i := C_0 Q^(i)_0/(zeta^i L)` is a DEFINITION (G1_4 s.3.1); the
    WKB form of colour `i` is `cor:set-cyclicPF5` (b), the covariance clauses are `lem:gau-colrows` (not CS s.1).
14. **Ordering of the connection problem**: s08 proves `D_0 = Phi^SP(z,1)/h(z)`; s09 proves Theorem Gamma and closes with `thm:gam-D`.
15. **GAM-2** (two Mellin-Barnes summands, same normalised series) is `lem:con-form` (d); `def:gam-c` and `prop:gam-rational` are placed in
    `sec:con-form` (s08) before it (DEC Q6), so that `thm:con-Pa` (i) has no forward reference; s09 uses them.
16. **Name of S15**: "the comparison theorem" (`thm:flh-main`).
17. **Lemma name clashes** (S01 C9; DC1 O-1): theory "Lemma 0" -> `lem:tra-frames`; CYC "Lemma 0" -> `lem:gau-const`; UR "Lemma 1"/"2" ->
    `lem:gau-lift`/`lem:gau-scal`; GI "Lemma 1" -> `lem:gau-uniq`; two "Prop. S" -> `prop:set-modrows` (UR) and `prop:gau-symp` (GI);
    "Theorem G" -> (G) `thm:gau-G` and the Gap Theorem `thm:mod-gap`; "Lemma C" -> `lem:mai-coord`; G1's Lemmas A-F -> `lem:set-grading`
    .. `lem:set-order3`; G4's "Lemma T" (trace formula for residues) -> a local name, never "Lemma T" (`lem:tra-T`); G5's Lemmas 3.A,
    3.B -> `lem:mai-Psi`, `lem:mai-Ghat`; "Theorem D" (UR) -> dropped; CP's "Theorem D" -> `thm:gam-D`; "(A2)", "(H)" never printed.
18. **Definition of `M_0`**: the model's Teleman-form R-matrix (`def:set-M0`); `cor:mai-natfact` proves that it is the conifold symbol of
    the quintic (UR s.2; G5_3 s.3.3 (b)), which recovers NF's definition.
19. **Flat and hatted outer terms**: `T_Gamma = F_Gamma` plus `G_Gamma = u^(2g-2) T_Gamma`; Theorem A is stated for `T_Gamma` (REV R6).
20. **R0 placement**: colour 0 of R0 (`thm:gau-R0`) is in s03; the colours `i != 0` follow from colour 0 by `lem:gau-colmod`,
    `lem:gau-colrows` and gradedness (G2R_2 s.2.9), in `cor:gau-colours`; the glue keeps only the evaluation at `u = 1`.
21. **Dilaton-leaf sign**: PPZ15 Def. 2.13's translation, `+z(1 - R^-1(z) 1)`; DNOPS s.2.3 prints the opposite sign (external misprint;
    DOSS s.2.3 and TR at `(1,1)` agree with PPZ) (G4R s.1.3 R-1; GS s.2.6).
22. **Pin of (M)**: Givental's Theorem 11.6 with the two-line rewriting for `s_i`; Cor. 11.7 prints `l_alpha` for `lambda_alpha` (G1R s.1.1;
    GS s.2.1).
23. **Trivial frames**: `frame(k,a) = z^k H^a` iff `k < nu_a + 1`, ten of them (G3 s.3 item 7; G3R_2 s.13).
24. **The quoted TR statement** is DOSS Thm 4.1 via DNOPS s.2.6 at `f = 4 W_v`, not "DNOPS Thm 8.1" (GS s.2.8; G4R s.1.3 R-3).
25. **`D_i` lies in `Q(zeta)[[z]]`**; in `Q[[z]]` only for `i = 0` (G2_6 s.11 item 3).
26. **`d(0) = 1` is proved** (`d = (P_0)_33 = 1/p_(0,0)`, `p_(0,a)(0) = 1` under (N0), GS s.2.3). The formula `d = tau^(1/4) varpi_0` is not
    needed by the chain; it agrees with `p_(0,0) = tau^(-1/4)/varpi_0` of `lem:gau-diag` (DB S27 steps 5-6), never cite S15's remark for it.
27. **Vertex tensor**: the homogeneous form `d^(2g_v-2+n) d^n Fbar(d q_T)` (DEC Q5; G3_3 s.5); it equals the Taylor expansion iff `d(0) = 1`.
28. **Bibliography key**: arXiv:1803.01409 is `LP20` (the digests and gap notes call it `LPfq`; G1 calls it "LP18"; GS s.2.11 (i)).

## 0.5 Decisions (settled; `design/DECISIONS.md`; each applied at its place)

* **D1**: `Pi_i`, `K`, `K_0` defined from `Y^C` (no `\draftnote`; the descendant reading is dropped under D8). **D2**: closed by G1:
  (M) pinned (G1 s.1 with GS s.2.1, DEC Q7), `prop:set-cyclicLR` proved (G1_2, G1_3). **D3**: Corollary B kept, proved (G5_2 s.2).
* **D4**: replaced by **Q4**: the second proof (ABDKS Thm 7.1 with the identity `S_h = -4^h B_2h/(2h)`) is the proof of record of the Gap
  Theorem; Hock's Cor. 1.2 becomes `rmk:mod-hock`; [Hock25] leaves (Q). **D5**: closed (G4_5). **D6**: closed (G4_2). **D7, D15**:
  closed (G2): Theorem U and `P_can` over any field of characteristic 0, hence over `Q`.
* **D8**: remarks: s.0 rule 10. **D9**: Theorem A, Corollary B, Theorem C, Theorem Gamma; "the comparison theorem". **D10**:
  `rmk:mai-constant` kept with G5_2 s.2.6. **D11**: pending; a `\draftnote` in s01 (companion notes; [Br25, s.1.1] sentence).
  **D12**: title and abstract of `design/TITLE_ABSTRACT.md` are drafts for the authors. **D13**: one model normalisation `W_v`; the
  others only in appendix B (with DNOPS's point `a = (0,0,0,0,4v)`). **D14**: frames `(alpha, beta)` in `Q^2`.
* **Q1**: (QRR) stays quoted; G1_4 s.3.2 Steps 1-5 become `rmk:set-qrrderiv` with one `\draftnote` (O1). **Q2**: (IFT) dropped; LMI
  step 2 uses `lem:fsp-ift`. **Q3**: one label `lem:set-M0sympl`. **Q5**: homogeneous vertex tensor in `def:red-flatsum`; `d(0) = 1` stated
  where the dilaton lemma is applied. **Q6**: `def:gam-c`, `prop:gam-rational` in `sec:con-form`. **Q7**: (M) pinned to Thm 11.6; LLY97
  only if the writer of s02 reads it and pins a statement (otherwise no second source and no `\draftnote`).
