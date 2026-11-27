# Assembly after the review (key review/assembler)

This file is for the authors and is not part of the note. Inputs: the eight reports `review/rf1..rf8/REVIEW.md`, `REVIEW_2.md`,
read in full, together with NOTATION (+ `_2`), OUTLINE s.0, GAPS_SYNTHESIS_2 s.3 and the section files.

What was done:
* (1) The cross-section changes of the reports were applied by local edits to 36 section files (s.1). Items that other
  referees had already applied were checked (s.2), and the items not applied are given with reasons (s.3).
* (2) The note compiles cleanly (s.5).
* (3) `NOTE/DRAFT_NOTES.md` (+ `_2`) lists every `\draftnote`.
* (4) The final PDF is `NOTE/main.pdf`.

The edits change no mathematics: no statement, no proof step and no label. The only change of wording in the conditionality
is item J below.

## 1. Edits, by file (report item in brackets)

* `s01_intro.tex`: $\overline M\to\overline{\mathcal M}$ (rule 6).
* `s01_intro_2.tex`:
  - after the display $\Qr\Dn=\Pc\Qr^{\mathrm m}$: "that $\Dn_0$ has rational coefficients is part of `thm:con-glue`" [rf3 X3,
    rf6 X2; `thm:gau-R0` (iii) now gives $\Dn_0$ over $\bbQ(\zeta)$ only];
  - the summand of Theorem Gamma is written in $\sigma$ [rf8 4].
* `s01_intro_3.tex`:
  - [J] (Q): [EO, Thm 4.7 (4-30)] and the variation (5-15) moved to the "cited but not used" paragraph [rf5 s.3 item 1: no proof
    uses them];
  - "PPZ15 enters ... in `sec:tra-dress`" instead of `def:tra-edge` [rf4: Prop. 2.12 and Def. 2.13 now sit after the
    definition];
  - LP20 pin "s.2.4.3, (17); s.2.4.4, (19)-(20)" [rf8 16].
* `s02_setup.tex`:
  - `lem:set-pf5` stated on any differential $\bbQ$-algebra containing $q$ with $\Dq q=q$ [rf3 X4; the proof is unchanged, and
    `cor:gau-cycsys` uses it in rings containing $\bbQ((u))$];
  - the duplicate Givental v2 draftnote in `thm:set-mirror` removed (rule 4; the first citation is s01).
* `s02_setup_3.tex`: LP20 pin, as in s01 [rf8 16].
* `s02_setup_4.tex`: `def:set-graphsum`: the loop-space point $q$ is declared [rf8 11].
* `s02_setup_8.tex`: `def:set-hat`: no letter $x$ for $s\varpi_0$ [rf8 7; NOTATION dictionary "S15's $x=s\varpi_0\to s\varpi_0$"].
* `s03_gauge.tex`: `def:gau-conc` (iii): $\det\Wr_{\mathbf x}=\Wr(x_0,\dots,x_3)$ of `lem:set-wronskian` [rf8 1].
* `s03_gauge_2.tex`: `thm:gau-cyc` and its proof: the scalar $\lambda\to\mu'$ [rf8 6].
* `s03_gauge_4.tex`: `rmk:gau-Uremarks`: $\cn_j\to\mathsf a_j$ [rf8 8; $\cn_n$ are the coefficients of $\cn(z)$].
* `s03_gauge_6.tex`, `s08_connection.tex`: HPF operators `\cO_e` $\to$ $\mathscr O_e$ [rf8 9; `\cO` is also $\cO(5)$].
* `s03_gauge_7.tex` (`thm:gau-R0` (i), (ii)), `s08_connection_8.tex` (`lem:con-N`): $\Lambda,\Lambda'\to\Lambda_s,\Lambda_\varsigma$,
  as in s07 and s10 [rf3 X6, rf6 O2; bare $\Lambda$ is the leg map].
* `s06_model.tex`:
  - "Two further statements of [EO], not among the inputs (Q), ..." (consistent with J);
  - `conv:mod-kernel`: $q$ is a point of the curve, and the parameter $q$ does not occur in `sec:mod` [rf8 11].
* `s06_model_3.tex`: $J^{(j)}_m$ declared local [rf8 14].
* `s06_model_10.tex`: `\texorpdfstring` in the title of `sec:mod-MLK` (removes the two hyperref warnings).
* `s07_comparison.tex`: $\overline{\mathcal M}$.
* `s07_comparison_2.tex`:
  - no letter $x$ in the hatted data [rf8 7];
  - number of leaf insertions $L\to n_{\mathrm l}$ in the proof of `lem:flh-hat` and in `thm:flh-main` Step 1 [rf1 X6;
    $L=u^{-1/5}$];
  - $\overline{\mathcal M}$.
* `s08_connection_3.tex`, `s08_connection_4.tex`, `s09_gamma.tex`, `s09_gamma_4.tex`, `b_signs_2.tex`, `c_evidence_2.tex`,
  `c_evidence_3.tex`: the Mellin variable is $\sigma$. Every $s$ with $\sigma=zs$ is replaced by $\sigma/z$ [rf8 4; NOTATION
  table 3 "the Mellin variable keeps $\sigma$"; $s=u^{1/5}$]. The substitution is exact: $\Gamma(\zeta^j/z-s)=\Gamma((\zeta^j-\sigma)/z)$,
  $\Gamma(1+s-\zeta^j/z)=\Gamma(1+(\sigma-\zeta^j)/z)$, $(-q)^s=(-q)^{\sigma/z}$, $y^{5s}=y^{5\sigma/z}$.
* `s08_connection_6.tex`: `lem:con-S2` and part (a) of the next proof: the variable $\eta\to\tilde p$ [rf8 3; $\eta$ is the pairing].
* `s09_gamma_3.tex`: after `thm:gam-D`: (QRR) enters in two places, through `thm:con-glue` (colour 0) and through
  `cor:gau-colours` (`cor:gau-rowsNN` (iii), `lem:gau-colrows`, colours $i\neq0$) [rf1 X1].
* `s10_main.tex`: `sec:tra-dress` [rf4].
* `s10_main_2.tex`:
  - `rmk:mai-constant`: $A_n\in\bbQ(\zeta)[L,L^{-1}]$ by `thm:gau-R0` (iii), and rational by `thm:con-glue` [rf3 X2, rf6 X2];
  - $\overline{\mathcal M}$.
* `s10_main_3.tex`: `rmk:mai-apriori`: the gauge factor is written $G(u,\zh)$ next to $G_\Gamma$ [rf8 12].
* `s10_main_4.tex`: `cor:mai-univ`: $\lambda$ declared local [rf8 6].
* `a_formalsp.tex`:
  - `def:fsp-wick`: $\eta\to\tilde t$ [rf8 3];
  - the opening paragraph declares $\xi$ as the variable of generating functions [rf8 13].
* `a_formalsp_5.tex`: `def:fsp-stirling`: $\mu$ local to the subsection [rf8 15]. The first placement gave a 55 pt overfull
  line, so the declaration was moved to the end of the definition.
* `b_signs.tex`: $J_k$ declared [rf8 14].
* `c_evidence.tex`:
  - the splitting $E_R=E_M+M^{-1}E_{R'}M^{-T}$ is an identity because $M$ and, by `thm:gam-D`, $R'$ are symplectic [rf2 cross-section];
  - $\overline{\mathcal M}$.
* `c_evidence_2.tex`:
  - natural factorisation: the symplecticity of $G$ is the only conclusion that fails for the perturbation tested, which breaks
    $\Dn_0(z)\Dn_0(-z)=1$; it is blind to $\exp(\text{odd})$ [rf1 X4, rf6 X3];
  - "Independence": the odd coefficients would be independent only if $\Dn_0(z)\Dn_0(-z)=1$ were known a priori
    (`rmk:con-genus`) [rf6 X3].
* `c_evidence_3.tex`: `sec:evi-univ`: "the polar part of the total does not move", as the sources say [rf1 X3].

Labels: none added, removed or renamed. Statements and definitions whose wording changed:
* `lem:set-pf5` (wider scope; its uses in s02_setup_3 and s03_gauge_2 are covered);
* `thm:gau-cyc`, `thm:gau-R0` (i), `lem:con-N`, `lem:con-S2`, `def:fsp-wick`, `def:set-hat`, `cor:mai-univ`, `def:fsp-stirling`
  (letters only);
* `def:gau-conc`, `def:set-graphsum`, `conv:mod-kernel` (a declaration added);
* `thm:int-gamma` (the summand written in $\sigma$).

Every use in other sections was checked, and none depends on the changed wording. All files stay below 250 lines; the longest
has 139.

## 2. Items found already applied (checked, no edit)

* rf1 X2, the pointer to $m_5$: applied in appendix B by rf8, and the values are now in `sec:evi-univ`.
* rf1 X4, `rmk:con-why`: applied by rf6.
* rf3:
  - X1, the glue over $\bbQ(\zeta)$: applied by rf6, and the statement now proves $\Dn_0\in1+z\bbQ[[z]]$;
  - X5: no edit was needed;
  - X7 ($\PFfive^L=\PFfive^{(0)}$) and X8 (index $r$): applied by rf6.
* rf8:
  - item 2 ($\Lambda_\kappa$): rf5;
  - item 5 (index $q$): rf6;
  - item 10 (discs $\bbD$): rf7.
* rf2's table and rf6 X1: `def:set-outer`, `lem:flh-inputs` (iii)-(iv), `thm:flh-main` Steps 1-2, `def:mai-Fcf`, `cor:mai-gap` and
  `rmk:con-genus` agree on the symmetrised edge and on the defect formula $R'(-z)^T\eta^{-1}R'(z)=M(z)(\Dn(z)\Dn(-z))^{-1}M(-z)^T$.
* rf4, rf5, rf7: the uses of the reworded statements in other sections are unaffected. These are `lem:tra-T` (c),
  `lem:tra-frames` (c), `lem:tra-conv` (c), `def:tra-edge`, the s06 items, and `thm:gam-D` with its frame scope.

## 3. Not applied, with reasons

* rf1 X5 (record the genus-3 classification in appendix C): the computation would have to be transcribed from the sources with
  its range. This is left to the authors, under the draftnote in `rmk:mai-analogy`.
* rf5 item 3 (an optional precision in the dilaton row of appendix B): rf8 found the row right as written.
* rf7's optional sentence in s01 (LMI step 2 needs no implicit function theorem): the s01 sentence is true as written.
* rf8 items 11 and 14 (the letters $q$ and $J$) and the letter $x$ of `lem:con-form`:
  - NOTATION sanctions these uses ($q_T$, $q_\phi$, $J^{\mathrm C}$, the flat phase $\Xi/x+\cG(x)$);
  - a global renaming would touch many statements;
  - the local meanings are now declared where they are used.
* The design files (NOTATION, LABELS, BIBPLAN, OUTLINE, SOURCE_ERRATA) were not edited. The corrections they need are in
  `DRAFT_NOTES_2.md` s.6.

## 4. Settled during the assembly: the open question of rf2 S11 for step C9

rf2 asked whether `thm:mod-tr` (C9) uses the comparison of the push-forward translation with insertions (the new draftnote in
`def:set-graphsum`). The quoted theorem is [Dunin-Barkowski:2015caa, s.2.6, (2-35)], an integral over
$\overline{\mathcal M}_{g,k}$ of the class $\alpha=R.\TFT^{\mathrm{TR}}$. That class is formed with the action of its s.2.3,
which decorates each dilaton leaf by $-\psi(1-R^{-1}(\psi))e_1$ and forgets the dilaton points with the forgetful map. This is
PPZ's push-forward form; it was read in the public arXiv text of 1509.06954 (ar5iv rendering). So C9 compares classes with
classes, and the comparison is not used there. The draftnote is not on the chain of Theorem A; it is used by Corollary B
(`DRAFT_NOTES.md` s.0).

## 5. Compilation

Scratch folder `allgenus_review_review_assembler/` (copies of `main.tex`, `preamble.tex`, `macros.tex`, `refs.bib`, `sections/`).
Each build ran pdflatex, bibtex and pdflatex twice, with a fourth pass when a rerun was requested.

| build | pages | errors | undefined refs/cites | multiply defined | warnings | overfull boxes |
|---|---|---|---|---|---|---|
| before the assembly (state after rf1-rf8) | 159 | 0 | 0 | 0 | 2 hyperref (title of `sec:mod-MLK`), 1 rerun after three passes | 6 |
| final | 159 | 0 | 0 | 0 | none (BibTeX: "missing pages in ABDKS25") | 5, all present before |

The final PDF was copied to `NOTE/main.pdf`. No other build product is in NOTE.

## 6. Files written

* `NOTE/DRAFT_NOTES.md`, `NOTE/DRAFT_NOTES_2.md`, this file, and `NOTE/main.pdf`.
* The 36 section files of s.1:
  - `a_formalsp`, `a_formalsp_5`;
  - `b_signs`, `b_signs_2`;
  - `c_evidence`, `_2`, `_3`;
  - `s01_intro`, `_2`, `_3`;
  - `s02_setup`, `_3`, `_4`, `_8`;
  - `s03_gauge`, `_2`, `_4`, `_6`, `_7`;
  - `s06_model`, `_3`, `_10`;
  - `s07_comparison`, `_2`;
  - `s08_connection`, `_3`, `_4`, `_6`, `_8`;
  - `s09_gamma`, `_3`, `_4`;
  - `s10_main`, `_2`, `_3`, `_4`.
* `main.tex`, `preamble.tex`, `macros.tex` and `refs.bib` are unchanged.
