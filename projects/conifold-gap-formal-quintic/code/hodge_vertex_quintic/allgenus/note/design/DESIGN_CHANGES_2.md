# Record of the design revision (key editor), part 2 of 2: gap-synthesis changes, inputs, editor's changes, open points, files

Continuation of `DESIGN_CHANGES.md`.

## 4. The design changes of `GAPS_SYNTHESIS_2.md` s.5 and `GAPS_SYNTHESIS_3.md` s.5-s.6: all applied

| synthesis bullet | applied in |
|---|---|
| OUTLINE.md: rule 8; s.0.1; chain table (C1 inputs, C1', C6, C9); final inputs; R7 table "none"; main theorem with `(alpha, beta)` in `Q^2` | OUTLINE.md s.0, s.0.1, s.0.2 |
| OUTLINE_2: D2, D3, D5, D6, D7, D15 closed with pointers; conflicts "DNOPS dilaton sign", "Cor. 11.7", "ten trivial frames" | OUTLINE_2 s.0.4 items 21-28, s.0.5 |
| s01 `sec:int-inputs`: the conditionality paragraph; no list of open steps | OUTLINE_3 s01 item 4 |
| `def:set-genus0`: `Delta_i`, `D u_i` as consequences | OUTLINE_3 item 2; `cor:set-cyclicPF5` |
| `thm:set-mirror`: (M) with pin and rewriting; (G0) next to it; Yukawa sentence to appendix C | OUTLINE_3 item 3 |
| `prop:set-cyclicLR` with Lemmas A-F, the G1R sentence, the corollary, local letters | OUTLINE_3 item 4; LABELS (`cor:set-cyclicPF5`, `lem:set-grading`..`lem:set-order3`) |
| `def:set-rows`: `Phi_i`, `A_n` as definitions; covariance clauses cite `lem:gau-colrows` | OUTLINE_3 item 7 |
| `thm:set-qrr`: `\draftnote` (a) deleted; (b) as remark (Q1); "8.7 Proposition (ii)" | OUTLINE_3 item 8 |
| `rmk:set-qrrvertex` subsumed (keep one) | dropped; `rmk:set-qrrderiv` kept |
| `def:set-graphsum`: naturality sentence | OUTLINE_3 item 9 |
| `def:set-M0` + `lem:set-M0sympl`; `def:set-modcohft` defined by `Deltahat`, `lem:set-modhom`; `lem:set-dilaton` with C4-1, C4-2, K4 | OUTLINE_10 items 14-15 |
| `def:set-confexp`: `NN` as in G5 s.2.2, `Cf` with the automorphism remark, forward reference | OUTLINE_10 item 11 |
| OUTLINE_4: `prop:gau-quintic` step; `thm:gau-U` statement, `lem:gau-lead`, `lem:gau-degree`, renames, optional shifts; the three lemmas; `prop:gau-graded`; `lem:gau-diag` with (N0); `thm:gau-G` Cor. (a); `prop:gau-symp` step 3; `lem:gau-wkb` every colour; `thm:gau-R0` (iii)-(iv); `cor:gau-rowsNN`, `lem:gau-colmod`, `lem:gau-colrows`; `cor:gau-colours` | OUTLINE_4 (all) |
| OUTLINE_5: `conv:tra-forms` two settings, summable families; the `Theta` lemma; `lem:tra-frames`; `lem:tra-conv` (a)-(f) and the two-slot lemma; `lem:tra-E` ring; `cor:tra-edges` step 3; vertex tensor; `thm:red-main` steps 8-9; `rmk:tra-K`; D1 | OUTLINE_5; LABELS (`lem:tra-sum`, `lem:tra-theta`, `lem:tra-twoslot`); `rmk:tra-K` dropped (D8) |
| OUTLINE_6: R3; `thm:mod-tr` (DOSS entry, TFT dictionary, `lem:mod-scale`, bookkeeping); `prop:mod-Rtr` with R-1, `lem:mod-laplace`; `def:mod-pair`; `lem:mod-xfun`; `thm:mod-Lreg` and route 2; `thm:mod-V` base `h = 1`, `lem:mod-w11` | OUTLINE_6; LABELS_2 |
| OUTLINE_7: `lem:flh-dil`; `thm:con-Pa` (G5 s.1.2, F1.1, F1.3, F1.4, Q6, FL4 lemma); `prop:con-deg` with (L4) | OUTLINE_7 |
| OUTLINE_8: `prop:gam-hatflat` cites `thm:con-Pa` (iv); `cor:mai-gap` definitions and proof; `rmk:mai-constant`; `cor:mai-natfact`; `rmk:con-why`; `thm:mai-main` conditionality | OUTLINE_8; OUTLINE_7 (`rmk:con-why`) |
| OUTLINE_9: appendix A "FL4 at the generic point" (`lem:fsp-generic`); B conventions in one place; C the phase's computations | OUTLINE_9 |
| NOTATION (GS_3 s.6) | NOTATION tables 1-5; NOTATION_2 tables 6-8, s.9, s.10, s.11 |
| LABELS (GS_3 s.6) | LABELS (+ `_2`): every proposed label registered; "(S1)" -> "(M)" in the scheme |
| BIBPLAN (GS_3 s.6) | BIBPLAN (+ `_2`): MR1408320 and DuninBarkowski:2012bw cited; AC96 new; PPZ15, MR2917177, GJR18, LP20, MR2276766, ABDKS25, EO, CGYZ26, Hock25 updated; optional and conditional keys |
| SOURCE_ERRATA (GS_3 s.6, items (1)-(13) and the gap-note corrections) | SOURCE_ERRATA #55-#74 |

## 5. Final inputs and conditionality (GS_2 s.3)

Printed in OUTLINE s.0.2: (Q) as GS_2 s.3.1 with Hock removed (DEC Q4); (QRR) as s.3.2 (every colour); (M) as s.3.3; (G0) as s.3.4; (T)
as s.3.5; (IFT) not an input (DEC Q2); the standard formula [AC96, (1.7)] of s.3.7 is NAMED in the conditionality sentence (the editor's
completion: `lem:set-dilaton` (i) is used by the chain through `lem:tra-F` and `lem:flh-dil`, and rule 8 asks for the exact list); the
"not used" list of s.3.8; the corollaries of s.3.9; the sentence of s.3.10 with these two adjustments. The R7 table reads "none".

## 6. Changes made by the editor (consistency of the design; each traced to its source)

1. **`Cont_Gamma` includes `|Aut Gamma|^-1`** everywhere (S15 s.1.1, confirmed by G5R F2.2 and GS s.2.7); `def:red-flatsum` and `def:mai-Fcf`
   carry no further factor (the earlier outline wrote `F^P_g = Sum |Aut Gamma|^-1 Cont_Gamma`).
2. **Eigenvalues of `H*`**: `zeta^i L/C_1` (G1_3 s.2.10 (a)); the earlier outline had `zeta^i L` (SOURCE_ERRATA #75).
3. **`d = tau^(1/4) varpi_0`**: removed from (I4) and from rule 3's list of computed-only numbers; it is not needed and agrees with
   `lem:gau-diag` (`p_(0,0) = tau^(-1/4)/varpi_0`, DB S27 steps 5-6), while GS s.2.3 calls the S15 remark VERIFIED only (SOURCE_ERRATA #76).
4. **R0 for the colours**: `thm:gau-R0` is stated in colour 0; the matrix identity `eq:gau-R0` and the colours `i != 0` are in
   `cor:gau-colours` (G2R_2 s.2.9), and the order of `sec:gau-R0` is `lem:gau-wkb`, `lem:gau-colmod`, `lem:gau-colrows`, `thm:gau-R0`,
   `cor:gau-rowsNN`, `cor:gau-colours`, so that Step 6 of G2_5 s.9.4 (membership in the ring) follows `cor:gau-rowsNN` without circularity.
5. **New labels for statements the gap notes name without a label**: `thm:set-G0` ((G0)), `def:mai-Fcf` (G5_2 s.2.2), `rmk:mod-hock`
   (route 1, DEC Q4), `lem:tra-sum` (G3 s.2.2), and the local lemmas of `prop:set-cyclicLR`.
6. **`lem:set-dilaton`**: the genus-1 logarithmic term is not printed (GS_2 s.5 offered "pinned or dropped"; dropped, it is not used).
7. **`lem:con-cfu`**: the attribution `\draftnote` is moved to `DRAFT_NOTES.md` (no `\draftnote` in the chain); "the uniformising map".
8. **`thm:mod-tr`**: the remark on the metric of DNOPS (2-7) is not carried (the route of G4 identifies the TFT through (2-30) and the frozen
   basis); DNOPS (2-7) leaves BIBPLAN.
9. **`def:set-fq`**: `F_g` as a function of `q` through `\Qn = e^(T(q))`; the convention `\draftnote` is replaced by the definition (G5_2 s.2.2).
10. **`rmk:gam-blind`** trimmed (the unwritten evaluation at `L = 0` is not printed); **`rmk:mai-tau`**, **`rmk:gam-prior`**,
    **`rmk:gam-sevenvar`**, **`rmk:gam-Mcomplex`**: their `\draftnote`s reworded without project references (rule 5 applies to draftnotes).
11. **Hypothesis names**: G2's "(H)" is not printed ("(H)" is never printed); (N0) and (EP) may be printed; `(F_g)` typeset with a roman F.
12. **Section titles and the file plan** fixed in OUTLINE s.0 rule 12 (used by the placeholders).
13. **Notation macros** collected in `NOTE/macros.tex` (NOTATION_2 s.11); `\ML{h}` and `\Fg{g}` take one argument.
14. **BIBPLAN_2 in final form**: `note` carries the version statement only; provenance in `annote` (the earlier BIBPLAN_2 asked for this
    move "before typesetting"); bookkeeping words removed from the annotes; the optional appendix-A context entries are dropped (D8).
15. **Page budget**: 84 pp (s02 +2, s03 +1, s04 +1, s06 +1, s08 +1; appendix A -1; remarks trimmed by D8).
16. **The model TFT**: defined by `Deltahat_i^(1/2) := zeta^(4i) varsigma` (GS s.2.10); its identification with the Hurwitz Frobenius
    manifold of DNOPS (8-1) is not printed as a remark in s02 (GS s.2.10 allowed one; D8; G4_2 s.7: not used); the name stays in appendix B.
17. **Digest pointers** written "(DB S3)" etc. where a bare "(S3)" could be read as the old input name (OUTLINE_4 `sec:gau-cyc`).

## 7. Open points for the main session (not applied; none changes the conditionality as decided)

* **P1 (LMI step 2 and DEC Q2).** REV_2 M1 observed that no source writes the formal route of LMI step 2 in full: LA s.2.1 says only
  that "a formal IFT over `C[[Q]]` would do", and DD2 GAM-4 step 2 uses the holomorphic family of critical points `p(Q)`, `|Q| < r`, its
  values at `Q != 0` (the identification LF-C1) and the holomorphy of `N_n(Q)` (for (ii)). With `lem:fsp-ift` alone one obtains `p(Q)` as a
  formal power series; the identification at `Q != 0` and the holomorphic extension need convergence. The design follows Q2; the writer
  of s09 must check that step 2 can be written with `lem:fsp-ift` (for instance through the rationality of `c_n`, `prop:gam-rational`).
  If it cannot, rule 1 forces a `\draftnote` in the chain, and the main session must either restore the holomorphic implicit function
  theorem as a cited standard result (with a pinned textbook statement) or supply the formal argument.
* **P2.** AC96: journal data from a reader's memory only (BIBPLAN "points to confirm"); the arXiv v1 is the cited version.
* **P3.** The DOSS version: arXiv:1211.4021 has one version (checked on the arXiv page), so "v1" is pinned; G4R read it without naming
  the version.
* **P4.** CGYZ26: the comparison (G5_2 s.2.6) was read in the HTML rendering; a PDF re-reading is pending (O3).
* **P5.** `rmk:set-qrrderiv` (O1): the route of GS s.2.2 (Newton's identities and a pinned Hodge integral) goes to `DRAFT_NOTES.md`.

## 8. Files written by the editor (all inside `NOTE`)

Design: `OUTLINE.md` .. `OUTLINE_10.md`, `NOTATION.md`, `NOTATION_2.md`, `LABELS.md`, `LABELS_2.md`, `BIBPLAN.md`, `BIBPLAN_2.md`,
`design/TITLE_ABSTRACT.md`, `design/SOURCE_ERRATA.md` (appended), `design/DESIGN_CHANGES.md`, `design/DESIGN_CHANGES_2.md`. LaTeX:
`main.tex`, `preamble.tex`, `macros.tex`, `refs.bib`, and the placeholders `sections/s01_intro.tex` .. `s10_main.tex`, `sections/a_formalsp.tex`,
`sections/b_signs.tex`, `sections/c_evidence.tex` (each: its `\section` line and its label; the writers overwrite them).
