# Draft notes (part 1 of 2)

Every `\draftnote` in the note, by category, with file, line and nearest label. Line numbers are those of the version after the
fix round (lemma `lem:set-pushins`, consistency fixes I-2..I-9, I-12, I-15, final assembly); files are in `sections/` unless
`main.tex` is named. The notes are printed in red in the PDF (`\draftnote` of `preamble.tex`). There are **16** of them, counted
both by a raw search for `\draftnote{` and by a parser that follows the `\input` chain from `main.tex` (the definition in
`preamble.tex` is not counted). Remove each note when its point is settled. Part 2 (`DRAFT_NOTES_2.md`): the notes removed or
changed, the points recorded without a draftnote, and the corrections that the design files need.

## 0. The chain of Theorem A

**None of the 16 notes lies on the chain of Theorem A.** This was checked on the files as they now stand, as follows.

* Proof closure. Starting from the proof of `thm:mai-main` (`s10_main.tex` l.25-29, with its dependency paragraph), every
  label cited in a proof was followed to its statement and that statement's proof, recursively. A proof is attributed to the
  statement it follows, or to the statements named in `\begin{proof}[Proof of ...]`. The closure has 182 labelled statements.
  Among them are `lem:set-pushins` (through `thm:mod-tr` and `prop:mod-Rtr` (ii)), `def:set-graphsum` (through `lem:tra-F`)
  and the quoted `thm:set-mirror`, `thm:set-G0`. The quoted (T), `thm:set-GT`, is reached through the definitions
  `def:set-frob` and `def:set-rows` (the $R$-matrix of the cohomological field theory), not through a proof; it is an input of
  the conditionality. **No draftnote stands inside any statement or proof of the closure.** The only remarks and conventions in it are
  `conv:mod-kernel`, `conv:tra-forms` and `rmk:fsp-classical` (a pointer from `thm:gam-ginf`), and none of them carries a
  draftnote.
* Statement edges. Adding the citations made inside statements enlarges the closure to 197 statements, because definitions
  point ahead. The only draftnote host reached is `rmk:mai-constant`, through the last sentence of `def:set-confexp`
  (`s02_setup_5.tex` l.64: "The comparison with the conifold expansion of [CGYZ26] is `rmk:mai-constant`"). That is a pointer,
  not used to define the expansion.
* Sections cited by the proofs of the closure: `sec:tra` ("the transport lemmas of"), `sec:gau-U` (the colour rotation),
  `sec:fsp` (formal integrands) and `sec:con-finite` ("Neither ... is used"). They contain no draftnote in a lemma, definition or
  proof. The one note in `sec:tra` stands in `rmk:tra-consistency`, a remark that no statement cites. The sections cited only
  from statements are `sec:con`, `sec:gam`, `sec:gau-van`, `sec:set-rows`, `sec:sgn`, `sec:sgn-tr`, `sec:sgn-lit`,
  `sec:evi-calib`, `sec:evi-model` and `tab:sgn-tr`. They are cited for letters, conventions, dictionaries or "not used"
  examples. Their draftnoted remarks (`rmk:gam-sevenvar`, `rmk:gam-Mcomplex`, `rmk:gam-prior`, `rmk:set-qrrderiv`) are cited
  by no proof.
* No proof cites appendix C, and no draftnote stands in appendix A, B or C.

Of the 16 notes:
* 2 are in `main.tex` (the title and the abstract) and 2 in the literature paragraph `sec:int-gap` of the introduction.
* 2 are bibliographic notes, in the introduction, at the first citations of two quoted inputs: Teleman's reconstruction (T)
  (`s01_intro.tex:12`) and Givental's mirror theorem (M) (`s01_intro_3.tex:11`). They concern only the numbering of the
  journal versions. The inputs are the pinned statements of the named arXiv versions (`sec:int-inputs`: "Statement numbers
  refer to the arXiv versions named in the bibliography"), so neither note is a step of the chain.
* 10 are in remarks, and none of these remarks is cited in a proof:
  - `rmk:tra-consistency` is cited nowhere;
  - `rmk:int-constant` and `rmk:mai-analogy` are cited only in the introduction;
  - `rmk:set-qrrderiv` is cited in the introduction, in the opening paragraph of s02 ("which no proof uses") and in the
    remark `rmk:set-qrrsign`;
  - `rmk:mod-hock` is cited in the introduction and in appendix C;
  - `rmk:gam-sevenvar` is cited in appendices B and C;
  - `rmk:gam-Mcomplex` and `rmk:mai-tau` are cited in appendix C;
  - `rmk:gam-prior` is cited in `rmk:gam-sevenvar`;
  - `rmk:mai-constant` is cited in the introduction, in appendix B and by the pointer in `def:set-confexp`.

The same check, run from `cor:mai-gap` (Corollary B), `thm:gam-D` (Theorem C) and `thm:gam-main` (Theorem Gamma), gives the same
result. The proof closures have 185, 109 and 52 statements, and none contains a draftnote host. With statement edges, again
only `rmk:mai-constant` is reached, through the same pointer. The draftnote in `def:set-graphsum` on which Corollary B depended
(REVIEW_SYNTHESIS_2, I-1) is gone. Its argument is now `lem:set-pushins` (`s02_setup_4.tex` l.109-140), which has a proof.

Two statements made in the previous version of this file are superseded:
* "the comparison of the two forms of the translation is not used by the chain";
* "Corollary B rests on a draftnoted sentence".

`lem:set-pushins` is on the chain of Theorem A, in the reading of the dilaton sign of the input of C9 (`prop:mod-Rtr` (ii),
`thm:mod-tr` Step 1). So [PPZ15, (8)] is now one of the pins of (Q). The quoted works are unchanged.

## 1. Citations to confirm or to pin (7)

| File:line | Nearest label | Note | Remark |
|---|---|---|---|
| s01_intro.tex:12 | sec:int-setting | Teleman [MR2917177, Theorem (2) in (1.15)], first citation: "numbers from arXiv v3; journal numbering not compared" | pin of (T); OUTLINE rule 4 (the only copy) |
| s01_intro_3.tex:11 | sec:int-inputs, item (M) | Givental 1996 [MR1408320, Theorem 11.6], first citation: "numbers from arXiv v2; IMRN numbering not compared" | pin of (M); the only copy |
| s01_intro.tex:40 | sec:int-gap | "Re-check the latest versions of arXiv:2609.13921 and arXiv:2509.19298 before submission." | [CGYZ26] and [Br25]; literature paragraph |
| s01_intro_2.tex:74 | rmk:int-constant | "Numbers from arXiv:2609.13921v2, read in the HTML rendering; to be re-read in the PDF." | [CGYZ26] Defs. 1.1-1.2, Ex. 1.3, Rem. 1.4, Prop.-Def. 4.5 (1), Conj. 1.5 |
| s10_main_2.tex:55 | rmk:mai-constant | the same note for the same numbers of [CGYZ26] | the remark after Corollary B (the gap for rescaled periods) |
| s02_setup_4.tex:77 | rmk:set-qrrderiv | Mumford's relation $\mathrm{ch}_{2k}(\mathbb E)=0$ and the non-vanishing of $\mathrm{ch}_{2k-1}(\mathbb E)$ "are used here without a pinned statement" | the genus-one values $\int\psi_1=\int\lambda_1=1/24$ of Step 3 are also unpinned; the remark is used in no proof |
| s06_model_7.tex:99 | rmk:mod-hock | "The closed forms printed in [Hock25, §3.2] are not used; only Cor. 1.2 there is cited." (reworded in the fix round, I-8) | [Hock25] is not an input (DEC Q4); the remark's degree count at $\infty$ and its value at $p=0$ are sketched and used in no proof |

## 2. Author decisions (8)

| File:line | Nearest label | Decision |
|---|---|---|
| main.tex:12 | (title) | Working title, to be settled (`design/TITLE_ABSTRACT.md`). The author line "Author(s) to be added" is a placeholder too. |
| main.tex:19 | (abstract) | Draft abstract, to be settled. Since I-3 it states the inputs of s01/s10 on one footing, with $g\geq2$. |
| s01_intro.tex:38 | sec:int-gap | Whether to cite the companion notes, and how to relate the note to the sentence of [Br25, §1.1] on $R$-matrix methods (DECISIONS D11). |
| s09_gamma_4.tex:88 | rmk:gam-sevenvar | Sketch only. Either include the full argument, which needs the formal delta lemma (the exact restriction to the constraint, with its Jacobian) that appendix A does not contain, or keep the sketch. |
| s09_gamma_4.tex:107 | rmk:gam-Mcomplex | Not proved: the proof needs the formal stationary phase of Euler's integral with a non-integer exponent, which appendix A does not contain. Keep it as an unproved remark, prove it, or drop it. |
| s09_gamma_4.tex:117 | rmk:gam-prior | Describe [MR1901075, §10] neutrally: decide whether it is presented there as an outline or as a proof (OUTLINE rule 11). |
| s10_main_3.tex:68 | rmk:mai-tau | Not proved; keep the remark or drop it. Appendix C records the computation at genus 5 modulo $u^8$ only. |
| s10_main_4.tex:79 | rmk:mai-analogy | The computation behind the last sentence (factorisations other than the natural form that make every edge graph of genus $\leq3$ regular) is not recorded in appendix C. Record it there with its range, or drop the sentence. |

## 3. Steps to check (1)

These are places where an argument is missing from the text or is given only inside the note.

| File:line | Nearest label | Note | Where it is used |
|---|---|---|---|
| s04_transport_3.tex:70 | rmk:tra-consistency | The remark uses the string equation for the ancestor correlators of the model, a standard consequence of its flat unit (`lem:set-dilaton` (ii)). The note says this "is not stated, proved or pinned to a reference in the note". | Nowhere: the remark is cited by no statement. The sources use the string equation without proof (rf4 O1). To close: state and pin the string equation, or drop the remark. |

The step-to-check note of `def:set-graphsum` (the push-forward translation against insertions) was removed. It is closed by
`lem:set-pushins` (part 2, s.4).

Continued in `DRAFT_NOTES_2.md`.
