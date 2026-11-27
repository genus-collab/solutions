# Gap phase synthesis (key gaps/critic), part 1 of 3: verdict, status of R1 and R7, adjudications

Author-facing; not part of the note; nothing here enters a `.tex` file. Object, read in full: the prover notes `NOTE/gaps/G1..G5/GAPS_G*.md`,
the referee reports `NOTE/gaps/G*/referee/REFEREE*.md`, their scripts and logs; `design/OUTLINE_REVIEW.md` R1 and R7; checked against
`OUTLINE.md` s.0-s.0.2, `OUTLINE_3.md` .. `OUTLINE_8.md` (the entries named), `NOTATION(_2)`, `LABELS(_2)`, `BIBPLAN(_2)`, `design/DECISIONS.md`.
Pieces: this file (s.0 verdict, s.1 status of every item, s.2 adjudications); `GAPS_SYNTHESIS_2.md` (s.3 conditionality of the main
theorem, s.4 what is still open, s.5 changes to the OUTLINE files); `GAPS_SYNTHESIS_3.md` (s.6 changes to NOTATION, LABELS, BIBPLAN and
SOURCE_ERRATA, s.7 questions for the main session, s.8 re-runs and own checks, s.9 files). Labels: PROVED, VERIFIED, CONJECTURAL, SKETCH,
OPEN. Own code and logs: `NOTE/gaps/critic/`. No existing file was modified. AI-written; no human has read it.
Abbreviations: G1..G5 = the five clusters; "G2 s.9.4" = the prover's note of G2, its s.9.4; "G2R s.2.10" = the G2 referee report, its s.2.10.

## 0. Verdict

* **Every chain item of R7 and of R1 is closed, with one qualification.** Each closure is a written argument (PROVED) or a citation whose
  statement a prover or referee read in the full text (CITED); every one was refereed. Of the 37 item verdicts of the five referees, 21 are
  *accepted* and 16 *accepted-with-corrections*; none is *rejected*, none was reclassified OPEN. I agree with every verdict (s.2 adjudicates
  the disputes). The qualification: item 2, `thm:set-qrr` (b), is a derivation of the quoted input (QRR) and still lacks two pins (O1); as
  long as (QRR) is quoted (from GJR18 Lemma 5.2 with (17), read by G1), it is not a step of the proof of Theorem A.
* **The table "steps of the chain left as `\draftnote`" required by R7 is empty after this phase.** The main theorem (`thm:mai-main`,
  Theorem A) is conditional exactly on the quoted theorems (Q), on (QRR), on (M), (G0), (T) and, unless the formal IFT is used, (IFT), each
  with a pinned statement (s.3). Compared with the design, (S1) is now the pinned (M) and (Q) changes in one entry (DOSS Thm 4.1 replaces
  "DNOPS Thm 8.1"; s.2.8).
* **What remains open** (s.4), none of it in the conditionality: (O1) `thm:set-qrr` (b), the derivation of (QRR) from the CohFT, uses two
  classical inputs nobody read in the full text (Mumford's relation; non-vanishing of `ch_(2k-1)(E)`); it matters only if the main session
  promotes (QRR) from a quoted input to a proposition. (O2) "the remarks" of R7's last line were not assigned to any cluster. (O3)
  bibliographic confirmations (journal numberings, one PDF re-reading), which are not mathematical gaps.
* **Repairs made inside this phase** (all written, all adopted here): G5's identity (1.1) holds only for amplitudes in `Q[sigma, z]` (F1.1);
  `prop:con-deg` needs the degree bound (L4) (F4.1); G3's (f2) needed the `m`-grouped argument (K3); G4's dilaton leaf had the misprinted sign
  of DNOPS s.2.3 (R-1); G1's pin of (M) was a misprinted corollary (Cor. 11.7 -> Thm 11.6). Each was found by the referee and is decided in
  the referee's favour.
* **Re-runs** (s.8): ten scripts of the provers and referees, copied to a scratch folder, reproduce their stored logs; my own four groups of
  exact checks pass, and their controls fail.

## 1. Status of every item (R1 with its fix items, then R7 in its order)

Resolution = the prover's; verdict = the referee's (A accepted, AC accepted-with-corrections); status = after this synthesis.

| # | item (outline label) | resolution (where) | verdict | status; what the editor must carry |
|---|---|---|---|---|
| R1.1 | `prop:set-cyclicLR`: `H * H^a = kappa_a H^(a+1)`, `kappa^LR = (1, L^5/(C_0^2 C_1^3), 1, C_0/C_1, C_0/C_1)` | PROVED, all orders in `q`, from (M) and (G0) (G1 s.2.0-s.2.8, Lemmas A-F) | A | **closed**; corollaries `Delta_i`, `D u_i`, "`C_0 v_0` solves PF5 for every QDE solution" (G1 s.2.10 (a)-(b)) |
| R1.2 | (S1) rewritten as what is used: the J-function | CITED: [MR1408320] Thm 11.6 (CY case `n=4, r=1, l_1=5, lambda'=0`) with the two-line rewriting for `s_i` | AC (pin Thm 11.6, not Cor. 11.7) | **closed** (s.2.1) |
| R1.3 | `cor:gau-rowsNN` (rows in the differential ring) | PROVED given `thm:gau-R0` (iii) (G1 s.2.10 (c)) | A | **closed**; place after `thm:gau-R0` |
| R1.4 | `prop:gau-quintic`: conifold expansion sends `C_1 kappa^LR_a` to `kappa_a`, `C_1` to `D t_c` | not assigned to a cluster | - | **closed here** (PROVED, s.2.9) |
| R1.5 | row C1' in the chain table | design | - | editor (s.5) |
| 1 | `thm:set-qrr` (a): `R(q=0) = diag(h_(w_i))` -> `A^QRR_n(q=0) = h_n` | PROVED (G1 s.3.1), with `Phi_i := C_0 Q^(i)_0/(zeta^i L)` as the definition | A | **closed** (all colours) |
| 2 | `thm:set-qrr` (b): uniqueness, `R(q=0)` is the vertex R-matrix | CITED + Step 1 PROVED (G1 s.3.2, Steps 1-5) | AC (C3b-1..C3b-4) | **closed except two pins** (O1, s.2.2); outside the conditionality while (QRR) is quoted |
| 3 | `lem:set-dilaton` | PROVED: PPZ15 Def. 0.3 (iii), Prop. 2.12, Def. 2.13; Arbarello-Cornalba (1.7) (G1 s.4) | AC (C4-1..C4-4) | **closed**; `d(0) = 1` is proved (s.2.3); `Int psi_1 = 1/24` only for the genus-1 log term, not used by the chain |
| 4 | `def:set-M0` symplecticity | PROVED twice, same argument (G1 s.5.1; G4 s.5 `lem:set-M0sympl`) | A (G1R); "correct" inside G4 item 1 (G4R s.1.2) | **closed**; one lemma (s.2.4) |
| 5 | `def:set-modcohft` homogeneity | PROVED (G1 s.5.2: degree count; Z/5 selection) | A (naturality sentence) | **closed**; define `omega^m_v` by its canonical data (s.2.10) |
| 6 | `thm:gau-U` steps 3, 5, 8, 10, 11, 12 (and 4, 6) | PROVED (G2 s.3-s.6; new `lem:gau-lead`, `lem:gau-degree`), over any field of characteristic 0 | A (all eight) | **closed** (D7, D15) |
| 7 | `prop:gau-graded` step 1 [M] | PROVED (G2 s.6.1-s.6.2; no `zeta`, no averaging) | A | **closed** |
| 8 | `lem:gau-diag` step 4 (G8) | PROVED (G2 s.7; the same line in G5 s.3.3 Rmk (i)) | A | **closed**, with the scope (N0) (programme family) |
| 9 | `thm:gau-G` Cor. (a) (G9) | PROVED (G2 s.8) | A | **closed** |
| 10 | `cor:gau-colours` (G6), with `lem:gau-colmod`, `lem:gau-colrows` | PROVED given G1 items 1, 2, 3(a) and `lem:gau-wkb` over `Q(zeta)` (G2 s.9) | AC (C1-C5) | **closed** (s.2.5): the G1 inputs are refereed; `lem:gau-wkb` over `Q(zeta)` follows by substitution |
| 11 | `lem:tra-frames` step 4 [S] | PROVED: the operator `Theta = 1 + N` is invertible (G3 part 2) | AC (K1, K8) | **closed**; ten trivial frames, not eight |
| 12 | `lem:tra-conv` (every "argued" item) | PROVED (a)-(d), added (e)-(g); the decay in (c) FOUND in GI s.2 (`lem:gau-vanshape`) | 2a, 2d, 2e, 2g A; 2b, 2c, 2f AC | **closed** (K2, K3, K4, K7, K10) |
| 13 | `lem:tra-E` step 2 (O-8) | PROVED: ring `Q((x))[[u, psi']]`, hypothesis `f` in `HH^u_-` (G3 part 4) | AC (K5, K6) | **closed** |
| 14 | `prop:mod-Rtr` (D6) | PROVED mod the quoted DOSS/DNOPS data: `X(lam) = M_0(t)^T`, `R_TR(w) = M_0(w/(kappa varsigma^4))`; transfer `kappa = 4 -> 1` (G4 part 2) | AC (R-1 required; R-2, R-3) | **closed** after R-1 (s.2.6); the outline's sign was right |
| 15 | `lem:mod-xfun` (G3) | PROVED, self-contained (trace over the sheets); EO v4 Thm 4.4, (A-25) as optional pins | A | **closed**; residue pairing for `LL` (F2) |
| 16 | `thm:mod-Lreg`, `lem:mod-abdks`, `prop:mod-period`, `lem:mod-vanper` (W5; route 2) | PROVED mod ABDKS Thm 7.1 and EO Thm 4.2; exact identity `S_h = -4^h B_2h/(2h)` replaces the Borel step | AC (R-4) | **closed**; W5 dissolves; route 2 is complete, so Hock is optional (D4) |
| 17 | `thm:mod-V` base `h = 1` (D5) | PROVED by a finite derivation by hand (`lem:mod-w11`; `-1/60` per branch point) | A (R-6) | **closed** |
| 18 | `thm:con-Pa` step 5 (the domain) | PROVED: flat identity in `A_L`, comparison in `Q((s))` via `L -> 1/s` (G5 s.1) | AC (F1.1 required, F1.2) | **closed** after F1.1 (s.2.7) |
| 19 | `cor:mai-gap` (Corollary B; D3) | PROVED given Theorem A; CITED PPZ15 Def. 2.13, Prop. 2.14; regrouping FOUND in `natural_form/genus4_referee/review/REVIEW.md` s.2.3 | AC (F2.1, F2.2) | **closed** |
| 20 | `cor:mai-natfact` (`Psi(u)`) | PROVED: `Psi(u) := Qhat(u, 0)`, `R~ := Psi^T eta^-1 Qhat` (G5 s.3) | AC (minor) | **closed** |
| 21 | `prop:con-deg` step 5 (resonance) | PROVED: (L1) left factor `1 - x^5` with ER's valuation lemma (G5 s.4) | AC (F4.1 required) | **closed** after (L4) (s.2.7) |
| 22 | "the remarks" (R7, outside the chain) | not assigned | - | **OPEN** (O2): rule 10, D8 |

Items raised during the phase that are not R7 items:
* (a) The clause `Phi_i(zhat, u) = Phi_0(zeta^-i zhat, u)` of `def:set-rows` needs a proof once `Phi_i` is defined by row 0 (G1R s.3):
  it is Steps 1-3 of `lem:gau-colrows` (G2 s.9.4; G2R s.2.10). **Closed.**
* (b) `def:set-graphsum`'s `\draftnote` (DA C6: the colour-basis form of PPZ Def. 2.2 not written in the sources): it is the naturality
  `A_*(R.omega) = (A R A^-1).(A_* omega)` for an isometry `A`, immediate from PPZ15 Defs. 2.2, 2.13 (G1R s.7). **Closed by one sentence.**
* (c) The vertex tensor for dressings with `d(0) != 1` is not a Taylor series (G3 s.3 item 2, part 3 s.5): a change of definition in
  `def:red-flatsum`, not a gap of the chain (for `P_can`, `d(0) = 1`).

## 2. Adjudications (prover against referee, and between clusters)

**2.1 The pin of (M) (G1 item 1).** G1 pins Givental's Cor. 11.7 and transcribes its factor (2) with `Sum_alpha lambda_alpha`; the referee
shows (font-level extraction, `ref_g1_pdfglyph.log`) that the PDF prints the Latin `l_alpha`, under which the factor at the formal quintic
point is `C_0 exp(5 g_1/(z C_0))` and (M) fails (+600 at order `q`; control A1 of `ref_g1_loc.log`, re-run identically). **The referee is
right.** Theorem 11.6 prints `lambda` correctly, and the rewriting is two lines: with `s^*_i(t) = e^(lambda_i t/hbar) Z^*_i(e^t)`, `s_i(T) =
e^(lambda_i T/hbar) Z_i(e^T)`, Thm 11.6 gives `s^*_i = f E_i e^(-lambda_i (T - t)/hbar) s_i(T)`, and the exponent is `(f hbar)^-1 [g_1 Sum_alpha
lambda_alpha - Sum_a lambda'_a g_(l_a) + lambda_i g_1 (Sum_a l_a - (n+1))]`, whose last term vanishes in the CY case. I re-derived it and
checked it on random rational data with two failing controls (non-CY degrees; the literal reading) (`critic_checks.log`, group A).
Decision: cite "[MR1408320] Theorem 11.6, Calabi-Yau case `n = 4`, `r = 1`, `l_1 = 5`, `lambda' = 0`, rewritten for `s_i`"; one clause says
that Cor. 11.7 prints `l_alpha` for `lambda_alpha` (external misprint, SOURCE_ERRATA). The translation (G1 s.1.3) and the hypotheses (s.1.4)
stand; the referee's degree-1 localisation confirms (M) at order `q` at generic weights and at the zeta point.

**2.2 The uniqueness step of (QRR) (G1 item 3(b)).** Steps 1-5 are correct mathematics (both agree). The referee's four citation points:
* C3b-1 (the GJR18 R-matrix solves the QDE): closed by Teleman v3 8.7 Proposition (ii) (`R = E exp h`, `E` solving the system of 8.2 (ii)),
  read by the referee.
* C3b-3 (Teleman's `T_z^-1 exp h T_z (I_A)` is `exp(h).omega` in PPZ's convention): closed by PPZ15 p.20 ("the same convention is used by
  Givental and Teleman") plus the referee's genus-1 computation `r_1 = a_1/12` (a sign check that would detect `z -> -z`).
* C3b-2 (Mumford's relation `ch_(2k)(E) = 0`) and C3b-4 (`ch_(2k-1)(E) != 0` in some `H^*(Mbar_(g,n))` for every `k >= 1`): **open**. The
  referee is right that G1's attribution of the injectivity to Teleman 8.10 is inaccurate and that Madsen-Weiss is the wrong name.
* Most promising route (not executed; nothing read): Newton's identities with Mumford's relation give `ch_(2g-1)(E) = (-1)^(g-1)
  lambda_(g-1) lambda_g/(2g-1)!` on `Mbar_g` (I re-checked `g = 2`: `lambda_1^2 = 2 lambda_2` gives `ch_3 = -lambda_1 lambda_2/6`); then a
  pinned non-zero evaluation of `Int_(Mbar_g) lambda_(g-2) lambda_(g-1) lambda_g` (for `k = g >= 2`) and `Int_(Mbar_(1,1)) lambda_1 = 1/24`
  (`k = 1`) give exactly what Step 5 needs (one `(g, n)` per `k`). Both inputs are classical Hodge-integral statements; MR1728879 (root key,
  already in BIBPLAN for Mumford's relation, s.2.1) is the first paper to read; no statement number is given here (rule 4).
* Consequence: the conditionality is unchanged while (QRR) is quoted (it is quoted in the design, with GJR18 (17), GJR17 (19)-(20), Giv01,
  MR2276766 as sources). Steps 1-5 then document how (QRR) follows from the CohFT; the decision whether to promote (QRR) is s.7 Q1.

**2.3 `d(0) = 1` (G1 C4-3; G3 K10).** The G1 referee calls `d(0) = 1` "an input of the chain" until it is derived from S15 (I4); the G3
referee lets it rest on `lem:gau-diag` step 4 (G8). **It is proved:** `d = (P_0)_33 = p_(0,3)` (S15 s.1.1-1.2; with `P_0^T eta P_0 = eta` at
`z^0`, `p_(0,3) = 1/p_(0,0)`), and G2 s.7 (accepted) proves `p_(0,a)(0) = 1` for every `a` for `P_can` under (N0), which holds in the programme
family. So `d(0) = 1` there; the S15 (I4) route and the formula `d = tau^(1/4) varpi_0` (VERIFIED in the sources) are not needed.

**2.4 Two proofs of the symplecticity of `M_0` (G1 s.5.1; G4 s.5).** The same argument (constancy of `Y(-z)^T eta^-1 Y(z)` for the self-adjoint
cyclic system, `bm_(a-1) = bm_(3-a)`, then homogeneity), refereed by both referees; G1R VERIFIED it to `t^24` with a gauge-equivalent but
non-self-adjoint control, re-run identically. One lemma in s02; one label (s.6). G4's F3 ("closes the R7 item") is right for symplecticity;
the companion item (homogeneity of the model CohFT) is closed by G1 s.5.2, not by G4 (G4R s.4, F3).

**2.5 `cor:gau-colours` (G2 item 5).** G2's one-word "PROVED" and G2R's "conditional on G1 (pending its referee)" are both correct for the
moment each was written; G1R has since accepted G1's items 1 (with corrections), 2 and 3(a), so the item is **closed**. G2R is right that the
parenthetical of Step 2 ("Step 1 with `i` replaced by `-i`") is wrong (C3). The cleanest form of Steps 1-2: the substitution `z -> zeta^-i z`
turns `PF5^(u_0)` into `PF5^(u_i)` exactly (every `z`-monomial of PF5 has degree 0 or 5), so it is a bijection between colour-0 and colour-`i`
WKB amplitudes, and `lem:gau-wkb` holds in every colour over any field containing `Q(zeta)`; equivalently, the colour-`i` transport operator
is `T(L g) = 5 zeta^(4i) D g` (G2R s.2.10). I re-derived the `z^1` part by hand and checked both statements exactly on `L^m`, with the control
PF5 `+ q z^2` (`critic_checks.log`, group C). G2's field correction (`D_i` in `Q(zeta)[[z]]`; in `Q[[z]]` iff `D_0` in `Q[[z^5]]`) is right.

**2.6 The dilaton-leaf sign (G4 item 1, R-1).** G4 calls the difference between DNOPS s.2.3 (`-psi (Id - R^-1(psi)) e_1`) and the outline
(`z(1hat - M_0^-1(z) 1hat)`) a convention; G4R calls it a misprint of DNOPS. **G4R is right**, for two independent reasons. (i) TR decides it:
the `(1,1)` graph sum with DNOPS's own `R^-1_1` matches the TR `omega_(1,1)` at all four points for `kappa = 1, 4` only with the PPZ sign
(`rf4_check_b.log`, re-run identically). (ii) With legs `R^-1(psi)` the translation is fixed by the unit, not by a convention: PPZ15 Prop. 2.12
needs `T(z) = z 1 - z R^-1(z) 1`; in rank one, G1R's genus-1 computation gives `r_1 (Int kappa_1 - Int psi_1) + r_1/2 = a_1/24`, i.e. `r_1 =
a_1/12 = a_1 B_2/2!` as Teleman 8.7 (i) requires, whereas `-T` would give `5 r_1/12 = a_1/24`. DOSS's own dilaton-leaf weight has the PPZ sign
(read by G4R). Nothing in the chain changes: the outline already uses PPZ's `T`.

**2.7 The two mathematical repairs of G5.** (F1.1) The counterexample is two lines: for `g = q`, `(zD + L)(L q c) = q (zD + L)(L c) + z (Dq) L c`
and `Dq = q`, so (1.1) fails by `z q L c`; Step 3 applies (1.1) only to `sigma^m` and `prod_k (5 sigma + kz)`, in `Q[sigma, z]`, so the proof stands
with the restricted statement. (F4.1) The bound `deg RHS_n <= 4n + 5` is needed and is not implied by (L0)-(L3); the referee's weight proof
(`x`: 1, `z`: -4, `theta`: 0; `M^5` and `Pi` share the weight-5 part `M_1^5`; all weights are `0 mod 5`) is correct. I re-derived it and checked
(L0)-(L4) in a third representation (P applied to monomials), with the controls `P + z^2 x^3` (breaks L1) and `P + z x^6` (breaks L4)
(`critic_checks.log`, group D). G5's other corrections (F2.1: the primitive of 1; F2.2: `|Aut Gamma|` counted twice, confirmed in
`stage2/flat/S15_PROOF.md` l.37, where `Cont_Gamma` already contains `|Aut Gamma|^-1`) are right.

**2.8 The operative quoted TR theorem (G4 R-3).** What the chain uses is: for this compact curve with simple branch points, TR equals the
graph sum (2-35) with `R` from (2-31) and the TFT `Delta_j^((2h-2+n)/2)`. In DNOPS this is the s.2.6 recollection of [DOSS] Thm 4.1 with the
compatibility of their Prop. 4.3 (1) and the sentence before Thm 8.1; the second half of Thm 8.1 (Frobenius manifold) is not used (F7). G4R read
DOSS (arXiv:1211.4021) s.2.3 and Thm 4.1 in the source. **Adopted**: (Q) names [DuninBarkowski:2012bw] Thm 4.1 and DNOPS s.2.6, Prop. 4.3 (1)
at `f = 4 W_v`, `a = (0,0,0,0,4v)` (R3), instead of "DNOPS Thm 8.1".

**2.9 R1 fix item 4 (not assigned; closed here, PROVED).** `Cf` is a differential homomorphism (`lem:gau-confexp`) with `Cf(C_0) = varpi_0`,
`Cf(L) = 1/s`, and `Cf(C_1) = Cf(W)/varpi_0^2 = D t_c` since `Cf(W) = varpi_0 D varpi_1 - varpi_1 D varpi_0 = varpi_0^2 D t_c` (G5 s.2.2,
recomputed by G5R). Applied to the forward relations of `prop:set-cyclicLR` it gives `z D Y_a = kappa_a Y_(a+1)` with `kappa = Cf(C_1 kappa^LR) =
(D t_c, 1/(u varpi_0^2 (D t_c)^2), D t_c, varpi_0, varpi_0)` (the outline's `L^5` is `1/u`). With `D = -(1-u) d/du`, `b_a = -kappa_a/(1-u) = (t_c',
-1/(u varpi_0^2 (1-u)^3 t_c'^2), t_c', -varpi_0/(1-u), -varpi_0/(1-u))`, which is G2 s.2's quintic system, and `prod_a b_a = -1/(u(1-u)^5)` as in
`prop:gau-quintic`. The closing relation in the frame is `cor:gau-cycsys` (or `Cf` of G1 s.2.10 (b)).

**2.10 The model TFT (G1 s.5.2 scope; G4 F7, F10).** `def:set-modcohft` both names `omega^m_v` "the Hurwitz Frobenius manifold of DNOPS (8-1)"
and gives `Deltahat_i^(1/2) = zeta^(4i) varsigma`. G1's homogeneity proof uses only the second; G4 shows the Frobenius-manifold
identification is not needed. Recommendation: **define** `omega^m_v` by `Deltahat_i^(1/2) := zeta^(4i) varsigma` (unit `Sum_i Deltahat_i^(-1/2)
e~_i`). The TR TFT of (2-30) is `Delta_j = -5 varsigma_j^2` (G4 s.4); in the frozen basis `e_j/sqrt(-5)` with the pairing rescaled, the TR
theory is the model theory with the genus-`h` part multiplied by `(-5)^(h-1)`, which is the factor between `F^TR_h = B_2h v^(2-2h)/(2h(2h-2))`
and `c_h v^(2-2h)`, `c_h = (-5)^(1-h) B_2h/(2h(2h-2))`. This is `thm:mod-tr` step 3, pinned by DOSS's vertex weight (R-3).

**2.11 Smaller points, each decided for the referee.** G3 K3 (I re-derived the supplement: `[u^o] T_m = c_0^m` times a non-zero polynomial in
`m`, the degrees `N_(m_2) - 1` being distinct); G3 K2 (NOTATION_2 s.9 maps GI's Lemma V to `lem:gau-vanmod`; the closed form is
`lem:gau-vanshape`); G2 C5 and every notation clash listed by the referees; G4 R-2, R-4..R-8; G5 F1.2-F1.4, F2.3-F2.5, F4.2-F4.3.
Two points of my own: (i) G1 writes "LP18" for arXiv:1803.01409v2, which BIBPLAN keys `LPfq` (BIBPLAN's `LP18` is 1702.06096v4): every
"LP18 (10), (17), (19)-(20), (26)" of G1 is `LPfq`. (ii) (QRR) now also enters s03: Step 3 of `lem:gau-colrows` uses the `q = 0` values of
the colours `i != 0`; OUTLINE s.0.1 ("(QRR) enters only in s02 and in the glue") and `thm:gau-R0`'s last sentence must say so (colour 0 of
`thm:gau-R0` still uses no (QRR)).
