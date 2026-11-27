# Consistency fixes I-2 .. I-9, I-12, I-15 (key fix/consistency)

Author-facing record; not part of the note. Source of the issues: `review/REVIEW_SYNTHESIS.md` and `_2` (s.4, s.6, s.7, s.10),
read together with the outputs of the lemma writer and its referee (`review/fix_lemma/LEMMA.md`, `REFEREE.md`), which I checked in the
files: `lem:set-pushins` stands in `s02_setup_4.tex` with its proof; s01 and s10 already carried its PPZ pins (§0.2 (1), (8), Def. 0.3,
the conventions after Defs. 2.2 and 2.5). Files edited: `main.tex`, `macros.tex`, `sections/s01_intro.tex`, `s01_intro_3.tex`,
`s04_transport_3.tex`, `s05_reduction.tex`, `s06_model.tex`, `s06_model_7.tex`, `s06_model_8.tex`, `s08_connection_4.tex`,
`s08_connection_6.tex`, `s09_gamma_3.tex`; one line each, for consistency (see s.3), in `s02_setup_4.tex`, `s10_main.tex`,
`s10_main_2.tex`. `refs.bib` unchanged (I-15). No design, digest, gap-note or review file of others was touched. Line numbers below are
those after the edits.

## 1. I-2: the role of (T) for Theorems A and C

**Choice: add (T), as the identification of R, to Theorem C's inputs** (option 1). Reasons, from the body:
* `def:set-rows` (s02_setup_4 l.9) defines R as "the R-matrix of its semisimple cohomological field theory [GJR18, Section 4.1]". That
  citation is part of the pin of (T) ("[MR2917177, Theorem (2) in (1.15)] with [GJR18, §4.1, Theorem 4.1]"), and the paragraph after
  `thm:set-GT` already says that (T) "identifies the R-matrix whose rows enter the outer graph sum". The object whose rows define the
  series D of Theorem C (`thm:gau-R0` (iii)-(iv), `cor:gau-colours`) is the same R, so the same identification is needed there.
* Option 2 is not supported by the body. No statement of the note says that R is characterised by (b) of `def:set-rows`, R(0) = 1
  and (QRR). The body proves only the uniqueness of a WKB amplitude up to a constant series (`lem:gau-wkb` (ii)), used colour by
  colour in `lem:gau-colrows` and `thm:gau-R0`. The existence of a series with (b) and the value (QRR) is not established anywhere
  independently of the citation of [GJR18, §4.1]. Writing the characterisation would be a new argument, and it would turn (T) into an
  interpretation-only input for Theorem A, against the conditionality sentence fixed by OUTLINE s.0.2 and rule 8.
* What the proofs use of R was checked in the files. For Theorem C (`thm:gam-D` <- `thm:con-glue`, `thm:gau-R0`, `lem:gau-colrows`,
  `cor:gau-rowsNN`, `cor:gau-colours`) it is (b) of `def:set-rows` (through `cor:set-cyclicPF5` (b) and the forward relations of
  `prop:set-cyclicLR`), R(0) = 1 (the leading terms `eta(1, e~_i) = Delta_i^(-1/2)`) and (QRR) (Step 3 of `lem:gau-colrows`, the glue).
  The symplecticity (a) is not used (`rmk:con-genus`, `s07_comparison_2` l.63, `s02_setup_8` l.38 say so). The note now states this for
  Theorem C, which makes the role of (T) there exactly "identifies R".

Edits:
* `s01_intro_3.tex` l.13, item (T): "In `thm:int-main,thm:int-D` it identifies that R-matrix, whose rows define T_Gamma and D; in
  `cor:int-gap` it identifies the total with F_g" (was: Theorem A only).
* `s01_intro_3.tex` l.21, Theorem C: "`thm:int-D`, like `thm:int-main`, concerns the rows of the R-matrix that (T) identifies
  (`thm:set-GT`); its proof uses this matrix only through `def:set-rows` (b), R(0) = 1 and (QRR), the last through the glue ... and,
  for the colours i != 0, through `lem:gau-colrows` ..., and it uses the genus-zero data ... (M) and (G0) ...". The rest unchanged.
* `s09_gamma_3.tex` l.92, the inputs paragraph after `thm:gam-D`: new sentence "They are the rows of the R-matrix of the cohomological
  field theory of the formal quintic, which the reconstruction (T) identifies (`thm:set-GT`), as in Theorem A; the proof uses this
  matrix only through `def:set-rows` (b), R(0) = 1 and (QRR)."
* `s02_setup_4.tex` l.151 (one-line consistency change required by I-2): "In Theorems A and C it identifies the R-matrix whose rows
  enter the outer graph sum (`sec:set-hat`) and define the normalisation series D (`thm:gau-R0`); in Corollary B ...".
* Not changed: the conditionality sentences of s01 l.19 and s10 l.35 (Theorem A), and s10 "Where the inputs enter" for (T).

## 2. The other issues

* **I-3 (abstract, `main.tex` l.18).** "in every genus" -> "in every genus g >= 2". The last two sentences now read: "The results hold
  modulo the following inputs from the literature, and every other step is proved: results on topological recursion, applied only to
  the model curve; results of Pandharipande, Pixton and Zvonkine on the action of R-matrices, including the formula for the pull-back
  of psi classes along the map forgetting a point; the value of the R-matrix at q = 0 given by quantum Riemann-Roch; the genus-0 mirror
  theorem and the genus-0 axioms of the twisted theory; Teleman's reconstruction; and a push-forward formula of Arbarello and Cornalba
  for products of psi classes. We do not treat the quintic threefold itself." So all inputs stand on one footing ("hold modulo"), as
  in s01/s10; "reconstruction" replaces "classification"; [AC96, (1.7)] (the push-forward formula, per REFEREE.md) and the lemma input
  [PPZ15, (8)] (the pull-back formula, in `lem:set-pushins`) are named in words. The `\draftnote` is kept.
* **I-4 (s01 lists).**
  - `s01_intro_3.tex` l.9, (Q), PPZ list: "(5)" added after "Def. 0.3 (iii)" (the group of symplectic series; PPZ's own name, REFEREE.md).
  - `s01_intro_3.tex` l.15, new sentence after "The definitions use ...": "The inputs are read with the following definitions and
    conventions of their sources: the normalisation of the J-function [MR2276766, Section 6, (10)] in (M) (`thm:set-mirror`); the free
    energies [EO, Def. 4.3, (4-10)] (`conv:mod-kernel`, `cor:mod-trF`) and the set-up [EO, (5-1), (5-4), (5-13)] of [EO, Theorem 5.1]
    (`thm:mod-var`); and the prescription [DOSS, §3.1, §3.5], by which the integrals [DNOPS, (2-31), (2-32)] are read as formal
    Laplace series (`lem:mod-laplace`, `prop:mod-Rtr`)." Each pin was checked at its place of use: s02_setup l.97, s06_model l.20-26,
    s06_model_8 l.39, s06_model_2 l.23. They are put here, not into the (Q) list, so that the list of (Q) printed in the s06 introduction
    (`s06_model.tex` l.7, not in my edit set) stays identical to s01's; the "Apart from the inputs above ..." sentence (l.23) is now exact.
  - `s01_intro_3.tex` l.17, "Where the inputs enter", PPZ places: `thm:set-G0` ("as part of (G0)", PPZ Def. 0.3) and `lem:set-modhom`
    (PPZ Defs. 2.2, 2.13, s02_setup_7 l.37) added.
  - `s10_main.tex` l.33 (one-line consistency change required by I-4): the same two labels added to s10's PPZ places, so that the two
    "Where the inputs enter" paragraphs agree again.
  - Cor. B lists: "(5)" added to the PPZ pins in `s01_intro_3.tex` l.21 and `s10_main_2.tex` l.40 (Step 2 uses "X symplectic" in the
    sense of `def:set-graphsum`); see I-5 for the other addition. The two lists stay identical (s01 writes "Remark", s10 "Rem.").
* **I-5 (pin changes ratified).** EO (4-30) and (5-15): already outside (Q) and in the "cited but not used" list of s01 l.23; s06 l.7
  says "not among the inputs (Q)"; s10 refers to s01's numbers; appendix B lists them in the convention table only. Nothing to change.
  PPZ Prop. 2.8: in (Q) with "Props. 2.4, 2.8 and 2.9 enter only through the proof of Prop. 2.14" (s01 l.9) and in the s10_main_2 l.31
  parenthesis; but the Cor. B lists named "Prop. 2.14" alone. Both lists (`s01_intro_3.tex` l.21; `s10_main_2.tex` l.40, a one-line
  consistency change required by I-5) now read "Prop. 2.14, with Props. 2.4, 2.8 and 2.9 through its proof". Checked against the
  public ar5iv text of arXiv:1303.1043: the proof of Prop. 2.14 invokes Prop. 2.9 and Prop. 2.8 explicitly, and combines the two
  R-actions, which is Prop. 2.4 (left group action of the action of Def. 2.2).
* **I-6 (`s01_intro.tex` l.16).** "(def:set-rows, prop:con-deg)" -> "(def:set-rows; by `thm:gau-R0` (ii)-(iii) it is a combination of
  the canonical amplitudes A^can_m, m <= n, to which `prop:con-deg` applies)". Checked: (ii) gives Phi^can solving HPF with
  coefficients in Q[[u]] and ph_0 = 1, so `prop:con-deg` gives A^can_m in L Q[L] of degree <= 4m+1; (iii) and its proof give
  A^act_n = sum_(m<=n) A^can_m [z^(n-m)] D_0^(-1).
* **I-7 (`macros.tex` l.1).** The parenthesis naming the design files is deleted.
* **I-8 (`s06_model_7.tex` l.99).** Draftnote now: "The closed forms printed in [Hock25, §3.2] are not used; only Cor. 1.2 there is
  cited."
* **I-9 (`s01_intro_3.tex` l.34).** "For thm:int-gamma we are not aware of a proof in the literature; the analogous identity ... is
  stated without proof in [GJR18, Remark 7.2], and ..." -> "For thm:int-gamma the note cites [GJR18, Remark 7.2], where the analogous
  identity for the Landau-Ginzburg integral of Guo, Janda and Ruan is stated without proof, and we do not claim ...". This keeps the
  REV R4 wording of rule 11 and drops the near-priority clause.
* **I-12 ("in" for the vertex law).** `s05_reduction.tex` l.59 (`thm:red-main` (b)): "Cont_vertex in c_g T_P^(2-2g) + Q[[u]]"; l.94
  (`cor:red-can`): "the vertex graph satisfies F_vertex in c_g t_c^(2-2g) + Q[[u]]"; for consistency also the display of Steps 8-9,
  l.85, and the same coset statement in `lem:tra-conv` (e3), `s04_transport_3.tex` l.24 (a file of my edit set through I-15). Not
  changed: `s10_main_3.tex` l.74 ("the full vertex term is ...", a remark; s10 is outside the allowed edits) and `s10_main_4.tex`
  l.29, 44, 56, which already say "lies in" or name a vertex law.
* **I-15 (typography).** All five overfull boxes removed, without changing any formula:
  - `s04_transport_3.tex` l.59 (8.6 pt): "... is (1/m_2!)<X y^(m_2)>_h times the binomial series sum_(m_1) ...", the phrase "the
    binomial series" being the one (f1) uses next.
  - `s06_model.tex` l.58 (the alignment, 1.5 pt): `\qquad` -> `\quad` before the list of triples in the display of `def:mod-xi`.
  - `s06_model_8.tex` l.20 (1.9 pt): "Then (2 pi i)^(-1) oint ... w equals -<(X dphi)_+, w>."
  - `s08_connection_4.tex` l.16-19 (2.2 pt): the formula for hat W'(X) set as a display, unchanged; "which is the formula stated and
    vanishes at X = 1" kept.
  - `s08_connection_6.tex` l.27 (3.9 pt): "(d) Let Delta_b in Q[[u]][[zhat]], b in Z/5, be series that contain only powers zhat^k with
    k = -b mod 5. If sum_b ... = 0, then every Delta_b vanishes."
  - ABDKS25: the arXiv page gives only "Journal of the European Mathematical Society, 10.4171/JEMS/1615, 2025", and the journal page
    (ems.press, via the DOI) says "published online first", with no volume, issue or pages. So `refs.bib` is left as it is, and the
    BibTeX warning "missing pages in ABDKS25" remains.

## 3. One-line changes in the protected files (as the task allows)

| file, line | issue | change |
|---|---|---|
| `s02_setup_4.tex` l.151 | I-2 | "In Theorems A and C it identifies the R-matrix whose rows enter the outer graph sum and define D (`thm:gau-R0`)" |
| `s10_main.tex` l.33 | I-4 | `thm:set-G0` (as part of (G0)) and `lem:set-modhom` added to the places where [PPZ15] enters |
| `s10_main_2.tex` l.40 | I-4, I-5 | Cor. B's PPZ list: "(5)" and "Prop. 2.14, with Props. 2.4, 2.8 and 2.9 through its proof" |

`s06_model_2.tex`, `s06_model_3.tex`, the other s02 parts and the other s10 parts are unchanged.

## 4. Build

Scratch folder `allgenus_fix_fix_consistency/` (`base/`: the note before the edits; `build/`: after). pdflatex, bibtex, pdflatex x3,
then two more pdflatex runs after the last small edits; `build/` was checked identical to NOTE (all `.tex` and `refs.bib`).
* Before: 162 pages; 5 overfull boxes (8.6, 1.5, 1.9, 2.2, 3.9 pt); underfull vboxes on two pages; BibTeX "missing pages in ABDKS25".
* After: 162 pages; no error; no undefined or multiply defined reference or citation; no LaTeX or package warning; **no overfull box**;
  BibTeX warning unchanged. One underfull vbox remains (badness 10000, page 7). On that page the inputs list of s01 now runs across the
  page break, and the page has almost no stretchable glue, because `\parskip` is a fixed 4pt (`preamble.tex`). The two underfull vboxes
  of the old build (pages 2 and 154) are gone. Shortening the new s01 text did not move the break; I left the preamble alone.
  The two underfull hboxes (badness 10000 in appendix C, 1215 in the bibliography) were in the old build too.
* NOTE/main.pdf was not rebuilt (compiling only in the scratch folder).

## 5. Open, for others (files outside this edit set)

* Design files are now behind on these points (I-16): OUTLINE s.0.2, item (T) ("In Theorem A it identifies that R-matrix (C1); in
  Corollary B, the total") should also name Theorem C; the PPZ pin list there lacks (5), Prop. 2.8 and the §0.2 (1)/(8) pins;
  `design/TITLE_ABSTRACT.md` still has the old abstract and its rule-8 check says [AC96, (1.7)] is not named in the abstract, which is
  no longer so; BIBPLAN should record the reading pins MR2276766 §6 (10), EO Def. 4.3 (4-10) and (5-1), (5-4), (5-13), DOSS §3.1, §3.5,
  and PPZ (5).
* DRAFT_NOTES(_2): the text of the `s06_model_7.tex` draftnote changed (I-8); the number of draftnotes is unchanged by this task.
* `s06_model.tex` l.7 (the list of (Q) in the s06 introduction) could name the reading pins of s01 l.15 as well. Not needed: they are
  not in (Q).
* Not done here: I-1 (done by the lemma writer), I-10, I-11, I-13 (author decisions and length), I-14 (stated in the note), and
  REVIEW_SYNTHESIS_2's I-12 on letters with several meanings (the task's I-12 is the "in" form only).
* Human check, in addition to REVIEW_SYNTHESIS_2 s.8: the statement in s01 l.21 and s09 l.92 that the proof of Theorem C uses R only
  through (b), R(0) = 1 and (QRR) (I traced `thm:gam-D` <- `thm:con-glue`, `thm:gau-R0`, `lem:gau-colrows`, `cor:gau-rowsNN`,
  `cor:gau-colours`, `cor:set-cyclicPF5` (b), `prop:set-cyclicLR` (iii)).
