# APPLIED_REF: referee fixes REF_1 .. REF_5 applied to the note (2026-11-18)

Order of application: REF_5, REF_2, REF_1, REF_4, REF_3. Edits by exact string replacement only, in `sections/*.tex`
(`main.tex` not touched: no finding concerns it). One paragraph = one line kept; no new label, macro or bib key.
Originals of all sources were snapshotted in the session scratchpad (`reffix/orig/`) before the first edit.

Totals: 30 findings (REF_5: 6, REF_2: 4, REF_1: 5, REF_4: 6, REF_3: 9), all applied; 41 string replacements, of which
40 verbatim from the reports and 1 ADAPTED; 0 SKIPPED. 14 files modified.

## REF_5 (`thm:mod-var` (ii), `cor:mai-univ`)
- F1.1 APPLIED `s06_model_8.tex` (proof (i): alpha, alpha' "letters of that paper").
- F1.2 APPLIED `s06_model_8.tex` ("the assumption before (5-4)").
- F1.2 (same pin, intro) ADAPTED `s01_intro_3.tex` l.15: the report gives no old/new fragment for this place; the pin
  `\cite[(5-1), (5-4), (5-13)]` was changed to `\cite[(5-1), the assumption before (5-4), (5-13)]` (only (5-4) pin of the paragraph).
- F1.3 APPLIED `s06_model_8.tex` (item (ii) now opens "(Not used.)"; theorem title unchanged, as instructed).
- F2.1 APPLIED `s10_main_4.tex` l.21 (lambda := -b_0(0)b_3(0)/v_A'(0)^2) and l.28 (quintic sentence, sign-convention sentence).
- F2.2 APPLIED `s10_main_4.tex` l.21 (replaces "(compare \cref{lem:mai-pole})").
- F2.3 APPLIED `s10_main_4.tex` l.31 (`lem:gau-vanshape` added).
- F2.4 no change required by the report; nothing done.

## REF_2 (`prop:mod-Rtr` (iii), `eq:mod-edge`)
- F1 APPLIED `s06_model_2.tex` l.57; APPLIED `s06_model_3.tex` l.29; APPLIED `s01_intro_3.tex` l.9.
- F2 APPLIED `s06_model_2.tex` l.57 (coordinate pin; merged with F1 on the same line, both present);
  APPLIED `s01_intro_3.tex` l.15 (end of the paragraph; merged with REF_5 F1.2 on the same line).
- F3 APPLIED `s06_model_3.tex` l.21 ("symmetric" removed).
- F4 APPLIED `s06_model_3.tex` l.25 (double Laurent series clause).
- Section 3, last bullet (optional addition to the table row of `b_signs_2.tex` l.25): not a finding, no text given; not done.

## REF_1 (`lem:set-Runique`)
- W1 APPLIED `s02_setup_4.tex` (first sentences of the proof).
- W2 APPLIED `s02_setup_4.tex` (existence and uniqueness; "unit of Q(zeta)[[q]]").
- W3 APPLIED `s02_setup_4.tex` (diagonal of N).
- W4 APPLIED `s02_setup_4.tex` (`def:set-rows`, last sentence); APPLIED `s02_setup_5.tex` (`thm:set-GT` parenthesis).
- W5 APPLIED `s01_intro_3.tex` l.13 and l.21 (`R|_{z=0}=1`), only these two places.

## REF_4 (`eq:con-bz`)
- W1 APPLIED `s08_connection_4.tex` l.83 (T_z on k[[t,z]]).
- W2 APPLIED `s08_connection_4.tex` l.83 (tensor product over Q, injectivity) and l.61 (`\Bbbk` in the list of local letters).
- W3 APPLIED `s08_connection_4.tex` l.83 (rho a unit).
- W4 APPLIED `s08_connection_4.tex` l.93 (w^{-1} log as a power series).
- W5 APPLIED `s08_connection_4.tex` l.95 (parenthesis on lambda replaced by the coefficient comparison).
- W6 APPLIED `s08_connection_5.tex` l.4 (a_0 defined and listed as local).

## REF_3 (membership, formal change of frame)
- F1 APPLIED `b_signs_2.tex` l.108.
- F2 APPLIED `b_signs_2.tex` l.106.
- F3 APPLIED `b_signs_2.tex` l.100, both replacements (wording "as read in version 2 of the preprint" kept).
- F4 APPLIED `s01_intro_3.tex` l.30.
- F5 APPLIED `s01_intro.tex` l.16.
- F6 APPLIED `b_signs_2.tex` l.108: step (2) sentence, step (3) definition of B_5, and `$\mathsf B_5$` added to the list of
  local letters (inside the F1 replacement text, after `$\mathsf B$`).
- F7 APPLIED `s02_setup_6.tex` l.20.
- F8 APPLIED `s01_intro_3.tex` l.32.
- F9 APPLIED `s10_main_2.tex` l.40; APPLIED `b_signs.tex` l.6.

## Checks
- Per edited file, before/after: differences of `{`/`}` equal, `$` difference even, `(`/`)` equal, `\begin`/`\end` equal,
  line count unchanged, CRLF kept, no `\label` added; also per changed line. All 14 files pass
  (script `reffix/counts.py` in the session scratchpad). Not compiled.
- Macros and keys used by the new text all occur already in the note: `\coef`, `\fps`, `\Dq`, `\bbH`, `\bw`, `\Bbbk`, `\AL`,
  `\Pc`, `\tc`, `\Bprod`, `\ev`, `\Cexp`, `\Nring`; bib keys `PPZ15`, `CGYZ26`, `Eynard:2007kz`, `DuninBarkowski:2012bw`,
  `Dunin-Barkowski:2015caa`, `MR2917177`; labels `lem:gau-vanshape`, `def:red-flatsum`, `sec:sgn-R`, `sec:sgn-lit`,
  `eq:mod-edge`, `prop:mod-Rtr`, `def:fsp-stirling`, `thm:gam-D`, `cor:set-cyclicPF5`, `def:set-frob`, `lem:set-Runique` exist.

## Lambda consistency check (`cor:mai-univ` and proof, after F2.1): CONSISTENT
With lambda := -b_0(0)b_3(0)/v_A'(0)^2:
- r = -1/lambda^2: even in lambda, unchanged (r = -v_A'(0)^4/(b_0(0)b_3(0))^2).
- P_A := lambda^{-1} P^A_can, so (P_A)_11(0,z) = lambda^{-1}; tilde r / r = (-1)/(-1/lambda^2) = lambda^2; the symplectic
  factor (tilde r/r) P_11(0,z) P_11(0,-z) = 1. Even in lambda.
- P_A Y^C(v_A) = lambda^{-1} Y^C_A, hence T_{P_A} = T_A/lambda (`lem:gau-vanq` (d)); item (ii) c_g (T_A/lambda)^{2-2g} is even.
- Quintic: `prop:gau-quintic` (ii) gives b_0(0) = 1, b_3(0) = -1; v'(0) = 1 (`s05_reduction.tex` l.84); so lambda = 1,
  P_A = P_can, T_{P_A} = T_A = t_c, and the corollary is literally `cor:red-can` with `thm:mod-MLK`, as the new sentence says.
- No other file uses this lambda (grep on `b_0(0)b_3(0)`, `T_A/\lambda`, `cor:mai-univ`: `s10_main_4.tex` only).
No inconsistency found; no mathematics changed beyond the report's text.

## Doubts and remarks (nothing fixed here)
1. `s10_main_4.tex` l.28, new sentence of F2.1: "since $R$, $E_P$ and $\tilde q_T$ are unchanged, $d$ changes sign": the letters
   $R$ and $d$ are not introduced in `sec:mai-univ`; they are those of `def:tra-dressing` / `def:tra-edge` ($R^{\mathrm{vec}}$, $d$).
   A reference to `def:tra-edge` could be added. Text applied exactly as given.
2. `s09_gamma_3.tex` l.82 still has "$R(0)=1$" in the sense z = 0 (same clash as REF_1 W5); not among the two places named,
   left unchanged.
3. REF_5 F1.2 in the intro is an adaptation (see above); the pin now reads "the assumption before (5-4)" inside a list of
   equation numbers.
4. REF_3 F3: the pins Def. 1.8 (2) and Ex. 1.3 of the preprint come from an automatic reading of the HTML (report, section 2);
   to be confirmed on the PDF. REF_3 F1 states that admissibility is asserted in the preprint and not proved here.
5. REF_3 F6: `$\mathsf B_5$` is declared local twice (list of local letters, and "(local)" in step (2)), as the report asks.
6. REF_1 W4 makes `def:set-rows` say that `lem:set-Runique` "uses only (M) and (G0)"; this is the lemma's hypothesis, a forward
   reference as before.
7. `APPLIED_1.md` doubt 9 (a_0 local to the proof of (i)) is settled by REF_4 W6.
