# Final assembly after the fix round (key fix/assembler)

This file is for the authors and is not part of the note.

Inputs, all checked against the files:
* the lemma writer and its referee (`review/fix_lemma/LEMMA.md`, `REFEREE.md`);
* the consistency fixes (`review/fix_consistency/FIXES.md`);
* `review/REVIEW_SYNTHESIS.md` and `_2`;
* the previous assembly (`review/assembler/ASSEMBLY.md`).

Files written:
* the comment lines of 43 `.tex` files (s.2);
* `NOTE/DRAFT_NOTES.md` and `_2`;
* `NOTE/main.pdf`;
* this file.

## 1. Build

**Where.** The scratch folder `allgenus_fix_fix_assembler/final/`, a clean copy of `main.tex`, `preamble.tex`, `macros.tex`,
`refs.bib` and `sections/`. After the edits it is byte-identical to NOTE. The run was pdflatex, bibtex, then pdflatex three
times. The first two passes after bibtex asked for reruns (undefined references while the citations were being read, then
"Label(s) may have changed"), and the third pass is clean.

**Result.**
* 162 pages (159 before the fix round; `lem:set-pushins` added about 1.5 pages).
* No error, no undefined reference or citation, no multiply defined label, no LaTeX or package warning, no overfull box.
* A parser of the `\input` chain finds 398 labels, 253 labelled statements and 172 proofs, with no duplicate label and no
  reference to an undefined label.

**What remains, all cosmetic:**
* an underfull vbox on p. 7, where the s01 inputs list runs across a page break and `\parskip` is a fixed 4pt;
* an underfull hbox on p. 153, in a table of appendix C (`c_evidence.tex` l.65-66);
* an underfull hbox on p. 162, the Teleman entry of the bibliography;
* the BibTeX warning "missing pages in ABDKS25" (JEMS, published online first, without volume or pages).

`NOTE/main.pdf` is this build: 162 pages, 2195033 bytes.

## 2. Scan of the `.tex` files, comments included

**Searched for:**
* bookkeeping words: stage, workflow, agent, referee, verdict, AI, programme, prover, cluster, script, notebook, scratch, file,
  design, outline, digest, synthesis, errata, decision, review, assembler, checked, verified, numerically, tested, confirmed,
  status tags;
* item tags of the reviews and the design;
* months and dates;
* file names with an extension, paths, e-mail addresses and URLs;
* design-file names and project names;
* personal data.

**Found and fixed: file names in the structural comments.** 47 comment lines in 43 files named `.tex` files of the note. For
example:
* `% s02_setup_4. Section s02, part 4 of 8:` now reads `% Section 2, part 4 of 8:`;
* `Continued in sections/s01_intro_2.tex.` now reads `Continued in part 2.`;
* `main.tex`, `preamble.tex` and `macros.tex` no longer name each other in their header comments;
* in `macros.tex`, `(s08)` and `(s06)` now read `(Section 8)` and `(Section 6)`.

Only comments changed, so the PDF has the same pages and the same size. The `\input` commands necessarily keep the part
names.

**Checked and left as they are:**
* "checked" and "tested" occur only in appendix C and in one pointer to it (`rmk:con-genus`). This is allowed by rule 3.
* "the sign is confirmed in genus one" (`rmk:set-qrrderiv`, Step 3) is a computation written out in the text.
* "bookkeeping" occurs only in its mathematical sense: `rmk:con-genus`, the leaf remark in s05, `s07_comparison_2.tex` l.96
  and `a_formalsp_5.tex` l.55.
* "The finite route" is the title of `sec:con-finite`.
* arXiv versions appear only in the two pin notes. There is no date, path, e-mail address, design-file name, script name or
  personal data. The author line is a placeholder.

**Outside the `.tex` files** (not in this task's edit set, not edited): the printed note of the Teleman entry in `refs.bib`
adds "(the ar5iv rendering renumbers sections 1 and 8)", which is a remark about how the paper was read. The other printed
notes are version pins ("numbers cited from this version").

## 3. Draftnotes

There are **16** (17 before; the note of `def:set-graphsum` went with `lem:set-pushins`): 7 citations to confirm, 8 author
decisions and 1 step to check (`DRAFT_NOTES.md` s.1-s.3).

**None lies on the chain of Theorem A** (`DRAFT_NOTES.md` s.0):
* The proof closure of `thm:mai-main` has 182 statements and contains no draftnote host.
* With the citations inside statements added, it has 197 statements, and the only host reached is `rmk:mai-constant`,
  through the closing pointer of `def:set-confexp`.
* No proof of the 172 cites appendix C or a draftnoted remark.
* The same holds from `cor:mai-gap`, `thm:gam-D` and `thm:gam-main`.

`lem:set-pushins` is now on the chain (C9, the reading of the dilaton sign), so [PPZ15, (8)] is a pin of (Q).

## 4. Remaining author decisions

1. **Title** (`main.tex:12`) and the author line.
2. **Abstract** (`main.tex:19`). Since I-3 it has the inputs of s01/s10; `design/TITLE_ABSTRACT.md` is behind.
3. **D11** (`s01_intro.tex:38`): whether to cite the companion notes, and how to relate the note to [Br25, §1.1].
4. **Version re-checks:**
   * the latest versions of [CGYZ26] (arXiv:2609.13921) and [Br25] (arXiv:2509.19298) (`s01_intro.tex:40`);
   * the [CGYZ26] numbers, read in the v2 HTML, to be re-read in the PDF (`s01_intro_2.tex:74`, `s10_main_2.tex:55`);
   * the journal numbering of Teleman (T) and Givental (M) (`s01_intro.tex:12`, `s01_intro_3.tex:11`).
5. **Remarks.**
   * `rmk:gam-sevenvar` (a sketch), `rmk:gam-Mcomplex` and `rmk:mai-tau` (unproved), `rmk:gam-prior` (neutral wording) and
     `rmk:mai-analogy` (an unrecorded computation): keep, prove, record or drop each.
   * `rmk:set-qrrderiv` (Mumford's relation is unpinned) and `rmk:tra-consistency` (the string equation is unpinned): pin or
     drop.
   * `rmk:mod-hock`: delete its note when the remark is accepted.
   * In any case, cut the remarks to the budget of D8.
6. **Pins to ratify:**
   * [PPZ15, (8)], or a textbook reference for $\psi_i=\pi^*\psi_i+D_i$, $\psi_iD_i=0$;
   * the deviations of I-5: EO (4-30) and (5-15) out of (Q), PPZ Prop. 2.8 into it;
   * the choice of I-2: (T) identifies $R$ in Theorem C as well.
7. **Length and a possible split** (REVIEW_SYNTHESIS_2 s.9). The note has 162 pages against the design's 84 (OUTLINE_9), and its
   remarks fill about 10 pages against about 2 (D8). The recommendation is to split it into two papers:
   * (I) the universal flat theory and the model: Theorem U, the gauge identity, transport, Theorem R, the model and the Model
     Lemma by topological recursion on the model curve, and universality in flat form (`cor:mai-univ`). Its inputs are (Q) and
     [AC96, (1.7)] only. About 60 pages.
   * (II)+(III) Theorem Gamma (about 25 pages, with no quoted input) together with the formal quintic (about 70 pages, with
     inputs (M), (G0), (T), (QRR) and [PPZ15]).

   A third paper, on local Mellin inversion, is possible. If the note stays one text, the recommendation is:
   * move the large-radius cyclic structure, the cyclic form of `PF4`, the u-adic convergence, the CGYZ26 dictionary and
     s10.2-10.3 to appendices or a complements section;
   * drop `sec:con-finite`, the remarks beyond D8 and most of appendix C, keeping its tables.

   That would give about 110-120 pages.
8. **Human checks of REVIEW_SYNTHESIS_2 s.8**, items 1-11, before further editing. Item 10 now reduces to [PPZ15, (8)] and to
   reading PPZ Def. 2.2 as a push-forward along §0.2 (1).
9. **Design files** (I-16): see `DRAFT_NOTES_2.md` s.6.
10. **Draftnotes:** switch `\showdraftnotesfalse` only when all 16 are gone.
