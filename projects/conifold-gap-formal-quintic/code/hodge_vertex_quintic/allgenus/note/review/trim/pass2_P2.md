# Trim pass 2, group P2 (s08, s09, appendix D, appendix F)

The build is in the scratch directory (pdflatex, bibtex, pdflatex x3). It used the other groups' files as they were at
the start of this pass. Result: 127 -> 124 pages. All three pages come from the files of this group. There are no
undefined references and no overfull boxes left in the P2 files.

## Pages (start page of each section, whole document)

| section | before | after | length before -> after |
|---|---|---|---|
| s08 The connection problem | 71 | 71 | 12 -> 10 |
| s09 The Gamma identity | 83 | 81 | 7 -> 7 (about 0.7 p shorter) |
| App. D Formal stationary phase | 107 | 105 | 11 -> 10 |
| App. F Evidence | 123 | 121 | about 1.5 -> about 2 (float placement only, no content cut) |
| total | 127 | 124 | -3 |

## Removed items (no proof of the closure of thm:mai-main / cor:mai-gap cites them; grep of all sections)

* `lem:con-u0` (the conifold value) with its proof: it was cited nowhere.
* `prop:con-deg` (the degree bound) with its proof and `eq:con-PL`: it was cited only by "uses neither (QRR) nor
  prop:con-deg" before thm:con-Pb, which now reads "It does not use (QRR)." The file `sections/s08_connection_2.tex` held
  only this proposition and is deleted; `s08_connection.tex` now inputs `s08_connection_3`. The polynomiality that
  thm:con-Pa (iv) needs comes from `lem:gau-poly`, as before.
* `lem:con-weight` with its proof and the sentence before it: it was cited only by the parenthesis "(also by ...)" in the
  proof of prop:con-res, which is removed. The proof already obtains the ring from V in z^{-1}Q[z^{-1}][[u]].
* `rmk:con-710` (said "not used"). `cor:con-frame` (frame independence): it was cited only in c_evidence_2, whose sentence
  now cites thm:con-Pb and def:con-PhiSP directly.
* `rmk:gam-necessity` (said "None of (a)-(c) is used"): the clause citing it after thm:gam-lmi is removed. This departs
  from TRIM_PLAN_2 s.7, which kept it in sec:gam-rmk. That subsection is renamed "Prior art" and keeps `rmk:gam-prior` with
  its draftnote.
* `rmk:fsp-classical` (said "None of this is used"): the parenthetical citing it in the proof of thm:gam-ginf now reads
  "(no asymptotic statement is used)".
* Two sentences: the closing sentence of appendix D (about lem:con-relab), and the counterexample sentence at the end of
  the proof of lem:fsp-generic (iii).

No statement of a theorem, proposition, lemma or corollary changed. No \draftnote was removed, because none of the
removed items held one.

## Tightened prose (outside statements)

* The s08 introduction. The subsection sec:con-problem is renamed "Conventions"; its label is kept, and the clause
  defining the flat amplitude A_n, used nowhere in s08, is dropped. The genus remark rmk:con-genus is shortened but keeps
  the same content and label.
* The s09 sec:gam-series paragraph; one redundant clause in the proof of thm:gam-main (ii); the "bound deg c_n <= 4n,
  not claimed" sentence after cor:gam-poly.
* The appendix D introduction, the proof of lem:fsp-ring (i), and the example parenthesis in the proof of
  lem:fsp-fubini.

## Overfull boxes fixed

* s08_connection.tex: three boxes, in the conventions paragraph (the PF5 formula is now its own sentence) and in
  rmk:con-genus (reworded).
* s08_connection_4.tex: one box, in Step 1 of the proof of thm:con-Pa (i) (reworded).

## References to removed labels outside the P2 files

* .tex files: none. A grep of every section, main.tex, preamble.tex and macros.tex finds no reference to con-u0,
  con-deg, con-weight, con-710, con-frame, gam-necessity, fsp-classical or eq:con-PL.
* Markdown planning files (not edited): LABELS_2.md, OUTLINE_7.md, OUTLINE_8.md, OUTLINE_9.md, NOTATION_2.md,
  DRAFT_NOTES_2.md (lines 44-45 and 53-54 describe these cuts as "not done" or "kept"), BIBPLAN.md, and the reviews and
  trim records (ASSEMBLY.md s.5, T6.md, TRIM_PLAN*.md, CHECK.md, REVIEW_SYNTHESIS*.md).
