# Draft notes (part 1 of 2)

Every `\draftnote` in the trimmed note (one paper, 119 pp, after the second trim pass and its check, 2026-11-15), by category, with file,
line and nearest label. Files are in `sections/` unless `main.tex` is named. The notes are printed in red in the PDF
(`\draftnote` of `preamble.tex`). There are **12** of them, counted by a raw search for `\draftnote{` over `main.tex` and
`sections/*.tex` and by a parser that follows the `\input` chain from `main.tex` (68 files; the definition in `preamble.tex`
is not counted). Remove each note when its point is settled. Part 2 (`DRAFT_NOTES_2.md`): the notes removed or added by the
trim, the points recorded without a draftnote, and the corrections that the design files need. The full version (16 notes)
is documented in `archive_full_162pp/`.

## 0. The chain of Theorem A

**None of the 12 notes lies on the chain of Theorem A or Corollary B.**

* 2 are in `main.tex` (title, abstract), 2 in the literature paragraph `sec:int-gap` of the introduction, and 2 are
  bibliographic notes at the first citations of (T) (`s01_intro.tex:12`) and (M) (`s01_intro_3.tex:11`). The inputs are the
  pinned statements of the named arXiv versions (`sec:int-inputs`), so these two notes concern only journal numbering.
* 6 are in remarks: `rmk:int-constant`, `rmk:mod-hock`, `rmk:gam-prior`, `rmk:mai-constant`, `rmk:mai-graphlocal`,
  `rmk:mai-necessity`. None of these remarks is cited in a proof:
  - `rmk:int-constant` is cited in the introduction and in the signs appendix (`b_signs_2`);
  - `rmk:mod-hock` in the introduction and in the evidence appendix;
  - `rmk:gam-prior` and `rmk:mai-necessity` nowhere;
  - `rmk:mai-graphlocal` in the introduction (`s01_intro_2`) and in the remark `rmk:flh-fail` (`s07_comparison_2`);
  - `rmk:mai-constant` in the introduction, the signs appendix and by the pointer in the last sentence of
    `def:set-confexp` (`s02_setup_6`), which is not used to define the expansion.
* The trim moved proofs to appendices (A: `sec:set-lrproof`, B: `sec:gau-cycproof`, C: `sec:tra-convproof`) with labels
  unchanged; none of the moved blocks contains a draftnote. No draftnote stands in a statement or a proof.

The proof-closure check (`review/trim/CHECK.md`, rerun after the second pass in `review/trim/CHECK2.md`) confirms that every
statement of the old closure of Theorem A and Corollary B is still present with its proof; no draftnote lies in a statement
or proof of that closure. The second pass removed and added no draftnote; only line numbers changed (updated above).

## 1. Citations to confirm or to pin (6)

| File:line | Nearest label | Note | Remark |
|---|---|---|---|
| s01_intro.tex:12 | sec:int-setting | Teleman [MR2917177], first citation: "numbers from arXiv v3; journal numbering not compared" | pin of (T); the only copy |
| s01_intro_3.tex:11 | sec:int-inputs, item (M) | Givental 1996 [MR1408320, Theorem 11.6]: "numbers from arXiv v2; IMRN numbering not compared" | pin of (M); the only copy |
| s01_intro.tex:32 | sec:int-gap | "Re-check the latest versions of arXiv:2609.13921 and arXiv:2509.19298 before submission." | [CGYZ26] and [Br25] |
| s01_intro_2.tex:72 | rmk:int-constant | "Numbers from arXiv:2609.13921v2, read in the HTML rendering; to be re-read in the PDF." | [CGYZ26] numbers |
| s10_main_2.tex:41 | rmk:mai-constant | the same note for the same numbers of [CGYZ26] | the remark after Corollary B (the dictionary itself is now in `sec:sgn-lit`) |
| s06_model_7.tex:87 | rmk:mod-hock | "The closed forms printed in [Hock25, §3.2] are not used; only Cor. 1.2 there is cited." | [Hock25] is not an input |

## 2. Author decisions (6)

| File:line | Nearest label | Decision |
|---|---|---|
| main.tex:12 | (title) | Working title, to be settled. The author line "Author(s) to be added" is a placeholder too. |
| main.tex:19 | (abstract) | Draft abstract, to be settled (left unchanged by the trim). |
| s01_intro.tex:30 | sec:int-gap | Whether to cite the companion notes, and how to relate the note to the sentence of [Br25, §1.1] on $R$-matrix methods. |
| s09_gamma_4.tex:72 | rmk:gam-prior | Describe [MR1901075, §10] neutrally: outline or proof. |
| s10_main_2.tex:46 | rmk:mai-graphlocal | **New (trim).** The genus-4 computation (six edge graphs singular in the normalisation of `rmk:mai-stored`) is no longer among the tables of the evidence appendix: record it there, in an ancillary file, or drop the first sentence. |
| s10_main_2.tex:51 | rmk:mai-necessity | **New (trim assembly).** The computations behind the first sentence (dressings violating (A), dressings not regular in $u$) are not among the tables of the evidence appendix: same choice. |

## 3. Steps to check (0)

The only step-to-check note of the full version (`rmk:tra-consistency`, the string equation for the model) went with the
remark, which the trim cut. No step to check remains.

Continued in `DRAFT_NOTES_2.md`.
