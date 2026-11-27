# Record of the design revision (key editor), part 1 of 2: scope, R1-R8, M1-M21, D1-D15, Q1-Q7

Author-facing; not part of the note; nothing here enters a `.tex` file. The design files were revised in place so that the writers can
work from them alone. Inputs applied: `design/OUTLINE_REVIEW.md` (R1-R8), `OUTLINE_REVIEW_2.md` s.2 (M1-M21, where no later decision
conflicts), `design/DECISIONS.md` (D1-D15, Q1-Q7), `gaps/GAPS_SYNTHESIS_2.md` s.3 and s.5, `GAPS_SYNTHESIS_3.md` s.5-s.6, the gap notes
G1-G5 and their referee reports. Part 2 (`DESIGN_CHANGES_2.md`): the gap-synthesis changes file by file, the final inputs, the editor's own
changes, open points for the main session, files.

**Files revised:** `OUTLINE.md` .. `OUTLINE_9.md` (rewritten), `NOTATION.md`, `NOTATION_2.md`, `LABELS.md`, `LABELS_2.md`, `BIBPLAN.md`,
`BIBPLAN_2.md` (rewritten), `design/TITLE_ABSTRACT.md` (rewritten), `design/SOURCE_ERRATA.md` (s.5 appended, items 51-77). **Files
created:** `OUTLINE_10.md` (s02 items 10-17: `OUTLINE_3.md` would otherwise exceed the piece limit), this record and its part 2, and
`main.tex`, `preamble.tex`, `macros.tex`, `refs.bib`, `sections/*.tex` (placeholders). Not touched: the digests, the gap notes, the reviews,
`design/DECISIONS.md`, every file outside `NOTE`.

**Pointer convention introduced** (OUTLINE header): `G1_3 s.2.8` = `gaps/G1/GAPS_G1_3.md` s.2.8; `G1R_2 s.4.2` = the referee report piece;
`GS`, `GS_2`, `GS_3` = the synthesis; `REV`, `REV_2`, `DEC`; mark [G] = written out in a gap note in final form. Every proof plan of a step
that a gap note closes now points to that file and section, with the referee corrections named next to it.

## 1. Required changes R1-R8 (`OUTLINE_REVIEW.md` s.1): all applied

| item | what was done | where |
|---|---|---|
| R1.1 | `prop:set-cyclicLR` added with statement, local lemmas A-F (`lem:set-grading`..`lem:set-order3`) and proof from G1_2-G1_3; no `\draftnote` (proved in G1, accepted by G1R) | OUTLINE_3 item 4; LABELS s02 |
| R1.2 | (S1) rewritten as (M) = the J-function statement, pinned to MR1408320 Thm 11.6 (DEC Q7); (G0) stated next to it | OUTLINE s.0.2; OUTLINE_3 item 3; TITLE_ABSTRACT |
| R1.3 | `def:set-confexp` points to `cor:gau-rowsNN` (new, after `thm:gau-R0`); FQ2 s.4 moved to appendix C | OUTLINE_10 item 11; OUTLINE_4 `sec:gau-R0` |
| R1.4 | the step "`\Cexp` sends `C_1 kappa^LR` to `kappa`, `C_1` to `D t_c`" added to `prop:gau-quintic` (GS s.2.9) | OUTLINE_4 `prop:gau-quintic` |
| R1.5 | row C1' in the chain table | OUTLINE s.0.2 |
| R2 | (S1)-(S3) renamed (M), (T); (S3) dropped (DEC Q2); (G0) added; dictionary entries for CP's "(S1)" and the sources' "(S2)"; `sec:evi-normal` writes `lem:con-S2` | OUTLINE s.0.2; NOTATION_2 s.9; LABELS scheme; OUTLINE_9 C |
| R3 | DNOPS at `f = 4 W_v`, `a = (0,0,0,0,4v)`, critical points `p^5 = v`, values `5 p_c^4`; the transfer `kappa = 4 -> 1` now by `lem:mod-scale` (G4_2 s.6); the quoted statement is DOSS Thm 4.1 via DNOPS s.2.6 (GS s.2.8) | OUTLINE s.0.2 (Q); OUTLINE_6 `thm:mod-tr`; OUTLINE_10 item 15 |
| R4 | "the analogous identity for their Landau-Ginzburg integral" in s01 items 2, 5 and `rmk:gam-prior`; SOURCE_ERRATA #51 | OUTLINE_3; OUTLINE_8 |
| R5 | CGYZ26 Thm 2 for the quintic only; the other targets dropped (option "drop"); SOURCE_ERRATA #52 | OUTLINE_3 item 2 |
| R6 | Theorem A stated for `T_Gamma := Cont_Gamma[E_(R'), q_(R'); M.\TFT]`, `R' := Q M^-1`; `\TFT`, `\hat\TFT` added to NOTATION | OUTLINE s.0.2; OUTLINE_10 item 16; NOTATION tables 2, 8 |
| R7 | the table "steps of the chain left as `\draftnote`" added; after the gap phase it reads "none"; O1-O4 listed outside the chain | OUTLINE s.0.2 (last paragraph) |
| R8 (a) | `s` -> `delta_P` (`lem:tra-T`, `prop:tra-can`); the slice scalar -> `mu` | OUTLINE_5; NOTATION_2 table 6 |
| R8 (b) | `c_10 = c_5^2/2` -> `[z^10] D_0 = ([z^5] D_0)^2/2` | OUTLINE_7 `thm:flh-main` |
| R8 (c) | HPF letters `A_e`, `Phi_j` -> `\cO_e`, `\Phi^{[j]}` | OUTLINE_4 `lem:gau-HPF`; NOTATION table 4 |
| R8 (d) | branch points `varsigma_b`; `Lambda = -5` -> the number `-5`; `Lambda(s)` -> `Lambda(X)` | OUTLINE_6 |
| R8 (e) | `C_0` -> `varpi_0` in `prop:con-res`, `thm:con-Pb` step 7; `psi_n`, `x` -> `A_n`, `L` in `prop:con-deg`; `HH` -> `\hat\cV`; `beta(Y)`, `Gamma` local | OUTLINE_7 |
| R8 (f) | `P_k(L)` -> local `p_k(L)` | OUTLINE_7 `def:gam-c` (now in s08) |
| R8 (g) | the degree variable `\Qn` | OUTLINE_3 item 1; NOTATION table 2 |
| R8 (h) | `\Nring` (the ring), `\Cexp` (the expansion map), `\logq` (formal `log q`) | OUTLINE_10 item 11; NOTATION tables 2-3 |
| R8 (i) | `\varkappa` for the dilaton/TFT scalar | OUTLINE_7 `def:flh-outer`; OUTLINE_10 item 15; NOTATION table 5 |

## 2. Recommended changes M1-M21 (`OUTLINE_REVIEW_2.md` s.2)

| item | applied? | how, or why not |
|---|---|---|
| M1 | superseded by DEC Q2 | (S3) dropped; LMI step 2 uses `lem:fsp-ift`. M1's observation that no source writes the formal route in full is recorded as open point P1 (part 2 s.7) |
| M2 | yes | `cor:mai-gap` cites PPZ15 Def. 2.13, Prop. 2.14 (Prop. 2.4 for the plain action), with (T) and the proved statements (GS_2 s.3.9) |
| M3 | yes, adapted to Q2 | `thm:gam-D` inputs: (QRR) through `thm:con-glue`, appendix A (with `lem:fsp-ift` through `thm:gam-lmi`), s03, s08 |
| M4 | yes | `thm:int-gamma` and `thm:gam-main` (iv) attribute `Ginf` to [GJR18, Lemma 7.7, (37)] |
| M5 | yes | (a) leg-free; (b) `p^4/4 + v/p`; (c) "two-parameter rational family"; (d) the mirror theorem stays (TITLE_ABSTRACT) |
| M6 | yes | "the computations at genus 2 of appendix C show ..." (s01 item 1); "the computations of appendix C show ..." (s01 item 3, `rmk:mai-graphlocal`) |
| M7 | yes | `m_5` value only in appendix C; `sec:sgn-norm` keeps the defining condition (OUTLINE_2 s.0.4 item 9; OUTLINE_9 B) |
| M8 | yes | colour-basis edge and leaves in appendix B; `def:set-graphsum` gets the naturality sentence (GS s.1 item (b)) instead of the `\draftnote` |
| M9 | yes | "the conifold frames of `def:set-frames`"; "programme" added to the words never printed (rule 5) |
| M10 | yes | `FF` spelled out in `thm:gam-main` (ii) (OUTLINE_8) |
| M11 | yes | (a) both non-formal inputs named; (b) `A[[z]]`-linearity in `def:fsp-wick`; (c) "Lemma T(b)" -> `lem:fsp-transl` (b) in NOTATION_2 s.9 |
| M12 | superseded by DEC Q4 | the second proof is the proof of record; W5 dissolved (G4_4 s.1; F6) |
| M13 | yes | (a) vertex genera `g_v`; (b) the translation `\ttr_R`; (c) `\mathsf k(z)`; (d) `iota_b`; (e) the state space `\bbH`; (f) `g_0 g_3 = 1` dropped; (g) dictionary entries |
| M14 | yes | s02 items assigned to `sec:set-fq` .. `sec:set-N`; new labels registered |
| M15 | yes | `thm:mod-var` in the two ranges of EO (5-14), (5-15) |
| M16 | superseded by GS s.2.8 | the (Q) entry is DOSS Thm 4.1 via DNOPS s.2.6 (2-30)-(2-35), which contains M16's list |
| M17 | yes | s.0.1 "Keep visible": (Q) also in s10 |
| M18 | yes | OUTLINE_2 s.0.4 item 8 (placement of (CYC)) |
| M19 | yes | `rmk:gam-sevenvar` as statement plus sketch; page budget revised (84 pp) |
| M20 | yes | key `LPfq` -> `LP20` (v2 of 18 April 2020, checked on the arXiv page); no (S3) entry (Q2) |
| M21 | yes | SOURCE_ERRATA #51-#54 |

## 3. Decisions D1-D15 and Q1-Q7 (`design/DECISIONS.md`)

| item | where applied |
|---|---|
| D1 | `def:tra-tower` from `Y^C`, no `\draftnote`; `rmk:tra-descendant` dropped under D8 (OUTLINE_5; LABELS) |
| D2 | closed: (M) with its pin, `prop:set-cyclicLR` proved (OUTLINE_3 items 3-4; OUTLINE_2 s.0.5) |
| D3 | `cor:mai-gap` kept, with `def:mai-Fcf` and the proof of G5_2 s.2.4 (OUTLINE_8) |
| D4 | replaced by Q4 (below) |
| D5, D6 | closed by G4_5, G4_2; their `\draftnote`s deleted (OUTLINE_6) |
| D7, D15 | closed by G2; Theorem U over any field of characteristic 0 (OUTLINE_4) |
| D8 | rule 10 rewritten; every optional remark not kept is marked "dropped (D8)" in LABELS; kept: the 7-variable chart (sketch), the necessity examples, LMI for complex `M_*` (one sentence), the weight lemma (one sentence), the finite route (`lem:con-N`, `thm:con-conn` only), the uniqueness of `tau` (unproved, labelled) and the non-optional remarks explaining mechanism or necessity |
| D9 | names as decided (OUTLINE_2 s.0.5; LABELS s01) |
| D10 | `rmk:mai-constant` (and `rmk:int-constant`) with G5_2 s.2.6; the comparison `\draftnote` replaced by the version `\draftnote` |
| D11 | default `\draftnote` in s01 (OUTLINE_3) |
| D12 | TITLE_ABSTRACT marked as drafts; `main.tex` carries the working title with a `\draftnote` |
| D13 | one model normalisation; the others only in appendix B, with `a = (0,0,0,0,4v)` |
| D14 | frames `(alpha, beta)` in `Q^2` in Theorem A and Corollary B; G5_2 s.2.5's complex-frames remark dropped (D8) |
| Q1 | (QRR) quoted; `rmk:set-qrrderiv` (new) with Steps 1-5 of G1_4 s.3.2, C3b-1, C3b-3, C3b-4 and one `\draftnote` (O1); `rmk:set-qrrvertex` dropped |
| Q2 | (IFT) dropped from the inputs; LMI step 2 by `lem:fsp-ift` (OUTLINE_8; OUTLINE_9 A); see P1 in part 2 |
| Q3 | one label `lem:set-M0sympl` (OUTLINE_10 item 14); G1's `lem:set-M0symp` not used |
| Q4 | `prop:mod-period` is the proof of record of `thm:mod-gap`; new `rmk:mod-hock`; [Hock25] removed from (Q) and from the conditionality |
| Q5 | `def:red-flatsum` with the homogeneous vertex tensor; `d(0) = 1` stated at `lem:tra-F`, `lem:tra-conv` (f), `lem:flh-dil` |
| Q6 | `def:gam-c`, `prop:gam-rational` placed in `sec:con-form` (labels kept; LABELS exception); `lem:con-form` (d) refers to them |
| Q7 | (M) pinned to Thm 11.6; LLY97 conditional (BIBPLAN key `LLY97`, cited only if read) |
