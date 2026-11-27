# Trim pass 2, group P1 (s01, s02, s10 + Complements, appendix A1, appendix E)

Files edited: `sections/s01_intro*.tex`, `s02_setup*.tex`, `a1_cyclicLR.tex`, `s10_main*.tex` (including the
Complements files `s10_main_3.tex`, `s10_main_4.tex`), `b_signs*.tex`. No other file was touched.

## Pages

Measured by compiling only these sections against the full `main.aux` (references resolved), in two ways:
(a) normal page layout (start of section to start of the next, in pages); (b) "content length": the same text on one
very tall page, divided by the normal text height (no page-break slack). Before = the 127-page version
(`allgenus_trim_trim_check` copy of the note); after = this pass.

| section | before (a) | after (a) | before (b) | after (b) |
|---|---|---|---|---|
| s01 Introduction | 5.42 | 5.20 | 5.46 | 5.22 |
| s02 Setup | 15.98 | 14.46 | 15.67 | 14.31 |
| s10 Main theorem | 4.14 | 3.84 | 4.01 | 3.82 |
| s11 Complements | 3.66 | 2.69 | 3.62 | 2.67 |
| A Large-radius cyclic structure | 2.62 | 2.52 | 2.54 | 2.54 |
| E Signs and normalisations | 5.31 | 4.73 | 5.31 | 4.69 |
| total | 37.14 | 33.44 | 36.61 | 33.25 |

Saving: 3.7 pp (a), 3.4 pp (b). The full build in the compile folder, which has the other groups' current files
too, gives 118 pages with no undefined references. No overfull hbox is left in these files. The b_signs box of 25.6pt
(colour-basis edge) is fixed by a display, and the full build reports no others from them. The two that remain,
12.9pt and 12.0pt, are in `a2_cyclicPF4.tex`, which is not in this group.

## What was done

* s02. I tightened the connecting prose and the definitions (fq, genus-zero data, Frobenius data, rows, graph sum,
  frames, conifold expansion, tau, model, modcohft, hatted theory, outer graph sum, N(g), rmk:set-cg) and turned short
  displays into inline formulas. I removed these derived displays and sentences, which no proof uses:
  * the u-form of PF5 (s08 states its own copy);
  * `eta^{-1}=eta/25`;
  * the frame coefficients `w_1`, `gamma_t`;
  * the sentence "uniqueness half";
  * the Phiact pointer.
  The paragraph "The symplecticity of R'" in def:set-outer now cites lem:flh-inputs (iii) instead of re-deriving that
  identity. The pin of (M) is kept step for step, only compressed.
* s01. I tightened the setting, the gap, the strategy and the plan. The Inputs, Conditionality and Where-the-inputs-
  enter paragraphs are unchanged in content.
* s10.
  * I tightened the proof of lem:mai-Psi, the justification after def:mai-Fcf and the proof of cor:mai-gap.
  * The sentence after cor:mai-gap no longer repeats the list of proved statements. It keeps "proved statements
    only" and the full list of PPZ15 items, which s01 refers to.
  * The input lists and the conditionality paragraph are unchanged.
* Complements.
  * **I cut `lem:mai-coord` (the coordinate-change lemma, "Lemma C"), its proof and the paragraph after it,
    about 0.75 pp.** No proof uses it: cor:mai-univ does not cite it, and nothing else does. The authors may want
    it back; it is in `archive_full_162pp/` and in the 127-page copy.
  * The paragraph after cor:mai-natfact is condensed, and rmk:mai-stored and rmk:mai-annih are tightened.
* A1. The closing remark on "every step is an identity in Q[[q]]" moved into the intro sentence, shortened.
* E.
  * sec:sgn-qrr is merged into sec:sgn-var. Its twist paragraph now points to rmk:set-qrrsign, and its asymptotic-
    series paragraph was dropped, because the same reference is in tab:sgn-lit.
  * The "Points" paragraph points to def:set-points.
  * sec:sgn-legs points to prop:mod-LD and def:mod-pair; the residue and circle pairing is in s06.
  * Weights and regularity are merged, the list of normalisations is shortened, and the R-matrix-action and LG-chart
    rows of tab:sgn-lit are shortened.

All theorem, proposition, lemma and corollary statements are unchanged, and every proof in the closure of
thm:mai-main and cor:mai-gap is kept complete. All `\draftnote`s are kept: s01 has 5 and s10 has 3, as before.

## Removed labels

* `lem:mai-coord` (s10_main_4)
* `sec:sgn-qrr` (b_signs)

Neither is referenced anywhere in `sections/*.tex` or `main.tex`, so no reference outside these files needs fixing.
The `.md` planning files mention `lem:mai-coord`: OUTLINE_2/8/9, LABELS_2, review/rf1/REVIEW*.md. These are records
and were not edited.
