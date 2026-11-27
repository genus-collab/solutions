# Draft notes (part 2 of 2)

Continuation of `DRAFT_NOTES.md`: s.0 is the chain of Theorem A, and s.1-s.3 list the 16 notes by category.

## 4. Notes removed, added or changed

| Place | Label | Note | Action | Reason |
|---|---|---|---|---|
| s02_setup_4.tex (was l.104) | def:set-graphsum | the push-forward translation against insertions | removed (fix round, I-1) | the argument is now `lem:set-pushins` (`s02_setup_4.tex` l.109-140), with a full proof from [PPZ15, (8)]; see below |
| s06_model_7.tex:99 | rmk:mod-hock | "were not re-checked" changed to "are not used; only Cor. 1.2 there is cited" | reworded (fix round, I-8) | process wording |
| s05_reduction.tex | rmk:red-leaf | "No derivation of the $u^2$-coefficient of $\delta_{P_{\mathrm{can}}}$ is written ..." | removed (rf4 F5) | $t_c-\hat v=u^2/30+O(u^3)$ is derived from `lem:tra-T` (d), `lem:gau-diag`, `prop:tra-can`, `def:set-frames`, `def:set-tau`, `lem:set-frobbasis` |
| s02_setup_4.tex, Step 1 | rmk:set-qrrderiv | second copy of the Teleman v3 note | removed (rf2 S10) | OUTLINE rule 4: the note stands at the first citation, `s01_intro.tex:12` |
| s02_setup.tex | thm:set-mirror | second copy of the Givental v2 note | removed (review assembly) | rule 4: the first citation of [MR1408320] is item (M) of `sec:int-inputs` (`s01_intro_3.tex:11`) |
| s10_main_4.tex:79 | rmk:mai-analogy | the genus-3 computation is not recorded in appendix C | added (rf1 S10-15) | a remark pointed to a computation that appendix C does not contain |

No note was added in the fix round, so the count went from 17 to 16.

The new `lem:set-pushins` (the translation: push-forward and insertion) is cited at the places that used the removed sentence:
* `thm:set-GT` "In particular", and the paragraph after it (`s02_setup_4.tex` l.144, l.151);
* the third property after `def:mai-Fcf` (`s10_main.tex` l.83);
* `cor:mai-gap`, Steps 1 and 2 (`s10_main_2.tex` l.25, l.31) and its list of inputs (l.40);
* the sign paragraph of s06 (`s06_model_2.tex` l.55), `prop:mod-Rtr` (ii) (`s06_model_3.tex` l.31) and `thm:mod-tr` Step 1
  (l.35);
* the lists of inputs in s01 and s02 (`s01_intro_3.tex` l.9, l.17, l.21; `s02_setup.tex` l.16; `s10_main.tex` l.33).

## 5. Points recorded without a draftnote (for the authors)

None of these is a gap in the chain of Theorem A.

* **Pins to ratify.** The conditionality is unchanged at the level of the quoted works.
  - [PPZ15, (8)] is the pull-back formula for $\psi$ classes, used for $d=1,2$ in the proof of `lem:set-pushins` (i). It is a
    formula recalled inside PPZ's proof of Prop. 2.12, and it now stands in (Q) (`s01_intro_3.tex` l.9). Ratify it, or
    replace it by a textbook reference for $\psi_i=\pi^*\psi_i+D_i$ and $\psi_iD_i=0$; that would need a new bibliography
    entry. [AC96, (1.7), (1.9)] cannot serve, since both are push-forward statements.
  - Added to the PPZ lists in the fix round: §0.2 (1), (5), Def. 0.3, and the conventions stated after Defs. 2.2 and 2.5.
  - [PPZ15, Prop. 2.8] is in (Q) through the proof of Prop. 2.14. Both lists of Corollary B read "Prop. 2.14, with Props. 2.4,
    2.8 and 2.9 through its proof" (I-5).
  - [Eynard:2007kz, Thm 4.7, (4-30)] and the variation (5-15) of [Eynard:2007kz, Thm 5.1] stand under "cited but not used".
  - The reading pins are named in the sentence "The inputs are read with ..." of `s01_intro_3.tex` l.15, not in (Q) (I-4):
    [MR2276766, §6, (10)]; [Eynard:2007kz, Def. 4.3, (4-10)] and (5-1), (5-4), (5-13); [DuninBarkowski:2012bw, §3.1, §3.5].
* **The role of (T)** (I-2, option 1). (T) identifies the $R$-matrix in Theorems A and C alike (`s01_intro_3.tex` l.13, l.21;
  `s09_gamma_3.tex` l.92; `s02_setup_4.tex` l.151). Human check: the sentence that the proof of Theorem C uses this matrix
  only through `def:set-rows` (b), $R(0)=1$ and (QRR).
* **DNOPS §2.3** prints its action without the factor $1/|\mathrm{Aut}\,\Gamma|$. The note reads it as PPZ's action ("edges as
  in [PPZ15, Def. 2.2]"), which supplies the factor, and DOSS has it. Optionally, say so next to the sign remark of (Q).
* **Human check 10 of REVIEW_SYNTHESIS_2 s.8** now reduces to two points: [PPZ15, (8)], and the reading of the contribution
  of PPZ's Def. 2.2 as the push-forward along the gluing map of PPZ §0.2 (1).
* **Optional pointers to `lem:set-pushins`**, not needed for Theorem A: `s07_comparison.tex` l.38 (`def:flh-outer`),
  `s02_setup_8.tex` l.28 (`def:set-outer`) and `s04_transport_4.tex` l.40.
* **Earlier referee points, unchanged.**
  - rf2:
    - (a) The translation (d) of the pin of (M) uses the string equation for stable maps without a pin; it belongs to (M).
    - (b) `def:set-fq`: the weights of the virtual normal bundle stay invertible at $\lambda_i=\zeta^i$. This is part of
      the definition quoted from [LP18, LP20], and the note does not write it.
    - (c) In `def:set-rows`, the properties after "This means" do not determine $R$; (T) and (QRR) fix it.
  - rf3:
    - $\Phi^\circ$ of `lem:gau-poly` (iii) is used nowhere.
    - `def:gau-van` (ii) cites `lem:gau-vanq` (b) ahead of it, which is not circular.
    - A few local letters inside single proofs reuse global ones.
  - rf4:
    - Three forward references in statements, by design: `conv:tra-forms`, `lem:tra-conv` (e3) and `def:tra-dressing` (ii).
    - Example (i) of `def:tra-dressing` leaves $v_A(u)$ implicit; `cor:mai-univ` makes it precise.
  - rf6: whether $\Dn_0(z)\Dn_0(-z)=1$ can be proved without Theorem Gamma is open, and it is not needed (`rmk:con-genus`).
  - rf7: appendix A writes the Jacobian matrices with the letter $D$; only `lem:fsp-generic` uses $\Dq$.
  - rf8: the response $\Delta_g(5)$ listed in OUTLINE_9 is not in appendix C.
* **Letters with several meanings**, each declared where it is used: $q$ (the parameter; a point of the loop space in
  `def:set-graphsum` and `lem:set-pushins`; a point of the curve in `conv:mod-kernel`; the vectors $q_\phi$), $J$, $x$,
  $\lambda$, $\mu$, $\xi$. A global renaming of the loop-space point would remove the most frequent clash.
* **Vertex law.** `s10_main_3.tex` l.74 (a remark) still writes "the full vertex term is $c_g T^{2-2g}+\mathbb Q[[u]]$",
  while the theorems now write "$\in$". This is harmless.
* **Wording.** "bookkeeping" occurs in the printed text only in its mathematical sense: the remark titles "Genus bookkeeping"
  (`rmk:con-genus`) and "the leaf bookkeeping" (`s05_reduction.tex` l.108), and `s07_comparison_2.tex` l.96,
  `a_formalsp_5.tex` l.55. "The finite route" (`sec:con-finite`) is a subsection title. Neither was changed.
* **Typography** (final build): no overfull box. Three boxes remain:
  - one underfull vbox on page 7, where the s01 inputs list runs across a page break and `\parskip` is a fixed 4pt;
  - two underfull hboxes, one in a table of appendix C (`c_evidence.tex` l.65-66) and one in the bibliography (Teleman).

  BibTeX warns "missing pages in ABDKS25": JEMS lists the paper as published online first, without volume or pages.

## 6. Corrections for the design files

These files were not edited: the design, digests and gap notes are outside the assembly.

* **NOTATION.**
  - The hatting matrices are $\Lambda_s:=\diag(s^{1-\nu_a})$ and $\Lambda_\varsigma:=\diag(\varsigma^{1-\nu_b})$, and bare
    $\Lambda$ is the leg map.
  - The HPF operators are $\mathscr O_e$.
  - The free constants of `rmk:gau-Uremarks` are $\mathsf a_j$.
  - The Mellin variable is $\sigma$ in every Mellin--Barnes summand.
  - The dictionary entry $x=s\varpi_0$ is now written $s\varpi_0$.
  - The moduli spaces are written $\overline{\mathcal M}$.
* **LABELS.** Add `lem:set-pushins` (`s02_setup_4.tex`). LABELS_2: add the display label `eq:tra-leaf` (`s04_transport_4.tex`).
* **BIBPLAN.**
  - PPZ15: the pins (8), §0.2 (1), (5), Def. 0.3, the conventions after Defs. 2.2 and 2.5, and "Prop. 2.14, with Props. 2.4,
    2.8 and 2.9 through its proof".
  - The reading pins of s.5 (MR2276766, EO, DOSS).
  - DOSS: §3.1 and §3.5, after `lem:mod-laplace`.
  - LP20: §2.4.3 (17) and §2.4.4 (18)-(20).
  - Hock25: also cited in `tab:sgn-tr`.
  - EO: Thm 4.7 (4-30) and (5-15) are recalled, not inputs.
* **OUTLINE s.0.2.**
  - Item (T) should name Theorem C as well (I-2).
  - The PPZ pin list lacks (5), Prop. 2.8, §0.2 (1) and (8).
  - Step C9 now uses `lem:set-pushins`, so (8) is a pin of (Q) for Theorem A.
  - Row C10: the comparison theorem uses PPZ Def. 0.3, Rem. 2.1 and (6), Defs. 2.2, 2.5, 2.13, Prop. 2.12 and the dilaton
    lemma, not Props. 2.4 and 2.9 (rf1 X7).
* **design/TITLE_ABSTRACT.md** still holds the old abstract. Its rule-8 check says that [AC96, (1.7)] is not named in the
  abstract, which is no longer so (I-3).
* **The derived value $u^2/30$** of `rmk:red-leaf` is still treated as computed in OUTLINE s.0 rule 3, s.0.2 list O2,
  OUTLINE_5 and GAPS_SYNTHESIS_2 s.4 (O-9). The note now derives it (rf4 F5).
* **design/SOURCE_ERRATA.md** (rf6 O1, rf3 F1).
  - Some sources carry the symplecticity of the rows in $q$ over to the conifold expansion, which is not proved; the note
    proves $R'(-z)^T\eta^{-1}R'(z)=M(z)(\Dn(z)\Dn(-z))^{-1}M(z)^{-1}$ instead.
  - Some sources take the actual colour-0 amplitude to be rational before (QRR) is used; the note proves this in
    `thm:con-glue`.
* **DECISIONS D8** (remarks about two pages) is not met: there are 38 remarks, about 10 pages (REVIEW_SYNTHESIS_2, I-10).
