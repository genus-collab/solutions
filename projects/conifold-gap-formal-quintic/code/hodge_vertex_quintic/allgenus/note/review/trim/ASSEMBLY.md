# Trim assembly (key trim/assembler)

Inputs: `TRIM_PLAN.md`, `TRIM_PLAN_2.md`, reports `T1.md` ... `T6.md`. Output: `note/main.pdf`, **127 pages**.

## 1. What the assembler changed

* `main.tex`: after `\input{sections/s10_main}` the line `\input{sections/s10_main_3}` (section 11, Complements); after
  `\appendix`, in the planned order, `\input{sections/a1_cyclicLR}`, `\input{sections/a2_cyclicPF4}`,
  `\input{sections/a3_uconv}`, then the old `a_formalsp`, `b_signs`, `c_evidence`. The header comment now says that
  Complements is input separately. Nothing else in `main.tex` (title, abstract) changed.
* `sections/s10_main_2.tex`, `rmk:mai-necessity`: one `\draftnote` added at its end (its first sentence cites computations
  that the trimmed evidence appendix no longer records; same situation as the note T1 added in `rmk:mai-graphlocal`). No
  other section file was edited: the groups had already rewritten every reference to a removed label.
* `DRAFT_NOTES.md`, `DRAFT_NOTES_2.md` rewritten for the 12 remaining draftnotes (16 before; 6 removed with cut remarks,
  2 added).

## 2. Dangling references

Grep over `main.tex` and `sections/*.tex` for every removed label (the 22 cut remarks of TRIM_PLAN s.4, `sec:con-finite`,
`lem:con-N`, `thm:con-conn`, `sec:mai-rmk`, the six removed `sec:evi-*` subsections, `lem:mod-vanper` removed by T5): **no
occurrence**, neither as `\label` nor in `\cref`/`\Cref`/`\ref`/`\eqref` (multi-label references included). So no reference
had to be fixed by the assembler. Labels moved to appendices are unchanged, and all their references resolve.

## 3. Build (scratch folder, pdflatex, bibtex, pdflatex x3, clean aux)

No error; no undefined reference or citation; no multiply defined label; no rerun request. BibTeX: one warning, "missing
pages in ABDKS25" (as before). 11 overfull hboxes (3-51 pt; listed in `DRAFT_NOTES_2.md` s.5), none in the earlier
162-pp build: a copy-editing point. All section files are below 250 lines (largest: `s04_transport_4.tex`, 156). No project
bookkeeping in any `.tex` file (grep for dates, review/design/gap identifiers, group keys: empty); no text names an
appendix by letter.

## 4. Pages per section (from `main.toc`; a section starting mid-page is counted from that page)

| part | pages (first-last) | pages | full version | plan target |
|---|---|---|---|---|
| title, abstract, contents | 1 | 1 | 1 | 1 |
| 1 Introduction | 2-7 | 5 | 8.5 | 5 |
| 2 The formal quintic, the model and the outer graph sum | 7-23 | 16 | 23 | 13 |
| 3 Cyclic systems, Theorem U and the canonical gauge | 23-38 | 15 | 23 | 14 |
| 4 The vanishing tower and transport | 38-47 | 9 | 13 | 8 |
| 5 Reduction to the model | 47-49 | 2 | 3.5 | 2.5 |
| 6 The model: topological recursion on the model curve | 49-65 | 16 | 20.5 | 15 |
| 7 Flat and hatted graph sums | 65-71 | 6 | 6.5 | 5 |
| 8 The connection problem | 71-83 | 12 | 15.5 | 10.5 |
| 9 The Gamma identity and the normalisation | 83-90 | 7 | 11 | 8 |
| 10 The main theorem | 90-94 | 4 | 5.7 | 4 |
| 11 Complements (old 10.2-10.3) | 94-98 | 4 | 4.3 | 3 |
| App. A The large-radius cyclic structure (new) | 98-100 | 2 | - | 3 |
| App. B The cyclic form of the fourth-order operator (new) | 100-104 | 4 | - | 4 |
| App. C $u$-adic convergence (new) | 104-107 | 3 | - | 2.5 |
| App. D Formal stationary phase (old A) | 107-118 | 11 | 11 | 9 |
| App. E Signs and normalisations (old B) | 118-123 | 5 | 5 | 5 |
| App. F Evidence: finite computations (old C) | 123-124 | 1 | 10 | 3 |
| References | 124-127 | 4 | 2 | 2 |
| **total** | | **127** | **162** | **117** |

Integer page differences overstate short parts (each is rounded up to the next start page); the groups' fractional figures
are in `T1.md` ... `T6.md`.

## 5. Over budget

127 pp against the target 117 and the ceiling 120 (the request: about 110-120). The excess sits where the groups reported
it: s02 (+3), s03 (+1), s06 (+1), s08 (+1.5), appendix D (+2), s01/Complements (+1). Each group stopped because what remains
is statement text, which the binding rules keep verbatim. Options that need the authors' approval, none a step of a proof:
* T6: cut `prop:con-deg` with `lem:con-u0` (about 1.6 pp; used by no proof since `sec:con-finite` is gone) and
  `lem:con-weight` (0.25 pp).
* T5: move the proof of `prop:mod-Rtr` (about 2 pp) to an appendix (saves nothing by itself), or shorten the explanatory
  sentences inside `conv:mod-kernel` and `def:mod-pair`.
* T4: shorten expository sentences inside `conv:tra-forms` and `def:flh-outer` (ii), or compress routine parts of the proof of
  `lem:tra-conv` ((b1), (b2), (b4), (f2)).
* T2: move the pin of (M) (about 1 p) into appendix A, or shorten the explanatory parts of the long definitions.
* Typographic: the overfull boxes and the reference list (4 pp) are untouched; a denser bibliography style would save about
  one page.
With the T6 cut and about one page from each of s02, s06 and appendix D, the note reaches about 122; reaching 120 needs edits
inside definitions or conventions.

## 6. Open points for the authors (also in `DRAFT_NOTES_2.md`)

* `rmk:mai-graphlocal` and `rmk:mai-necessity` (and the sentence of `s01_intro_2` pointing to the former) cite computations
  that are only in `archive_full_162pp/` (draftnotes in both remarks).
* D8 vs the removal list (TRIM_PLAN_2 s.7 item 1): `rmk:gam-sevenvar`, `rmk:gam-Mcomplex`, `sec:con-finite`, `rmk:mai-tau`
  are cut as the task requires.
* MR1728879 is no longer cited (BibTeX omits it); the entry stays in `refs.bib`.
