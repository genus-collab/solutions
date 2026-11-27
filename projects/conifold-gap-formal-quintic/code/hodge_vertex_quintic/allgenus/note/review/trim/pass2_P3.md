# Trim pass 2, group P3 (s03, a2, s04, a3, s05, s06, s07)

Build: scratch folder `allgenus_trim2_trim2_P3`, pdflatex, bibtex, pdflatex x3. The full note now has 119 pages (other groups
were editing at the same time). No error, no undefined reference or citation, no multiply defined label. BibTeX gives one
warning, as before ("missing pages in ABDKS25"). No overfull hbox is left in the P3 files.

## Pages

Section lengths are fractional. Each was measured from `\pdfsavepos` marks before every `\input` of `main.tex`, in a copy of
`main.tex` that exists only in the scratch folder. The "before" column is the first build of this pass, before any edit.
The "after" column is the final build. Expect about +-0.15 pp of noise from page breaking, and from other groups' edits
moving the section starts. s04, s05, s07, a2 and a3 were not tightened, so their small changes are that noise.

| part | before | after | change |
|---|---|---|---|
| 3 Cyclic systems, Theorem U and the canonical gauge | 14.68 (pp. 23-38) | 13.91 (pp. 21-35) | -0.77 |
| 4 The vanishing tower and transport | 9.05 | 9.00 | (noise) |
| 5 Reduction to the model | 2.05 | 2.04 | (noise) |
| 6 The model: topological recursion | 16.28 (pp. 49-65) | 15.02 (pp. 46-61) | -1.26 |
| 7 Flat and hatted graph sums | 5.39 | 5.53 | (noise) |
| App. B cyclic form of PF4 | 3.19 | 3.26 | (noise) |
| App. C u-adic convergence | 3.02 | 3.07 | (noise) |
| **s03 + s06** | **30.96** | **28.93** | **-2.03** |

## What was removed (no proof uses any of it)

* Whole statements outside the proof closure of thm:mai-main / cor:mai-gap, none referenced by any proof:
  `cor:gau-Gstar` (the dual identity), `cor:mod-basis`, `cor:mod-primary`.
* Subsection heading `sec:gau-symp`. `prop:gau-symp` now sits at the end of `sec:gau-van`, whose title now reads
  "The residue relation, the vanishing solutions, (G) and symplecticity".
* Expository paragraphs and sentences with no role in any proof, in s03:
  * the use of cor:gau-cycsys;
  * the characterisation of the vanishing solutions after lem:gau-uniq;
  * the remark that colour 0 uses no (QRR);
  * the remark that gradedness is essential;
  * a duplicate pointer to the appendix;
  * the motivation sentence in def:gau-Fg.
* The same in s06:
  * the "for orientation only" examples in conv:mod-kernel, and two pointer sentences there;
  * the circle-versus-residue caveat in def:mod-pair;
  * the exponent remark after thm:mod-Lreg;
  * the remark after thm:mod-var;
  * the termwise-application pointer after thm:mod-MLK;
  * the reading of the diagonal of (2-31);
  * the sentence on how the three sources place the dilaton leaves;
  * the parentheticals "(It is the lattice theorem ...)" and "(The decomposition is also lem:tra-frames (a))".
* Proof steps were shortened, with the logic unchanged:
  * prop:mod-Rtr Step 3 now cites lem:gau-colmod (i), (ii) instead of repeating its computation;
  * thm:gau-G (c) now cites the proof of lem:gau-vanq (a), and that proof now names lem:gau-vanmod explicitly;
  * the HPF proof writes PF5 as `uM^5-1-(1-u) sum sigma_k z^k M^{5-k}`, and r_1 inline;
  * the proof of thm:gau-U (c) pointer moved into the proof of (b)-(d);
  * the remark after lem:mod-w11 (the coordinate-free form of Loc) moved into its proof;
  * the citation sentence after lem:mod-xfun moved into its proof;
  * the parentheticals in thm:mod-var and in the proof of thm:mod-V (i) were compressed.
* Overfull boxes fixed:
  * s03_gauge_6 (cor:gau-Gstar proof; now removed, 22.6 pt);
  * s03_gauge_6 (HPF proof, 3.2 pt);
  * s03_gauge_7 (lem:gau-wkb (i), 23.8 pt, now a display);
  * the colour-rotation item in s03_gauge_2, which overflowed after a rewording and is now a display.

No statement of a theorem, proposition, lemma or corollary was changed; the only differences that
`check_scripts/closure.py` reports for P3 files are in definitions and conventions, and only in their explanatory sentences.
The draftnote count in the P3 files is unchanged (1, in rmk:mod-hock, kept; the remark was shortened). No conditionality
text was touched.

## References to removed labels

* Inside the P3 files: `s04_transport_2.tex` had one sentence citing `cor:mod-basis`, which is removed. `s03_gauge.tex`
  cited `sec:gau-symp` in the section introduction; it now cites only `sec:gau-van`.
* In section files outside P3: **none**. Grep over `main.tex` and `sections/*.tex` for `cor:gau-Gstar`, `cor:mod-basis`,
  `cor:mod-primary` and `sec:gau-symp` is empty.
* The planning files `OUTLINE_4.md`, `OUTLINE_5.md`, `OUTLINE_6.md`, `LABELS.md` and `LABELS_2.md` still list these labels.
  They are not part of the build and were not edited.
* `rmk:mod-hock` and `rmk:mod-offK` are referenced from `s01_intro_3.tex` and `c_evidence_3.tex`, so they were kept
  (shortened). `cor:mod-onepoint` is referenced from `c_evidence_3.tex`, so it was kept.
