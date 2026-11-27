# APPLIED_4: log of the application of AMEND_4 (2026-11-18)

Applied to `sections/*.tex` and `refs.bib` with exact-string edits (Edit tool only), after AMEND_1, AMEND_2, AMEND_3
(`APPLIED_1.md`, `APPLIED_2.md`, `APPLIED_3.md`). Nothing compiled. `main.tex`, `macros.tex`, `preamble.tex` unchanged.
No new label, macro or bibliography key. All AI work; no source paper was opened, so every corrected pin rests on the
readings reported in AMEND_4.

Counts, items A4-01..A4-32: APPLIED 18, ADAPTED 9 (A4-01, 04, 05, 07, 16, 17, 19, 21, 27), ALREADY DONE 2 (A4-02,
A4-06), NOT APPLIED as duplicates 2 (A4-03, A4-09), SKIPPED 1 (A4-24, no text). Rejected R1-R8: not applied.
Counts, pin table rows 1-18: APPLIED 12, ADAPTED 2 (rows 6, 7), ALREADY DONE 4 (rows 9, 15, 17, 18; for 17 and 18 a
remaining part was applied). NO LONGER IN THE TEXT: none of rows 1-18.

## P1

| item | status | file |
|---|---|---|
| A4-01 | ADAPTED. (a): the sentence had been changed by A2-22 to "Then at ..., with $\Qn=1$ and $t=T(q)H$ (read through the divisor equation, with $\re^{dT}=(\re^T)^d$ as in def:set-genus0; see (d) below),". The whole clause was replaced by the A4-01 text (which now gives this reading as a definition); the A2-22 parenthesis is thereby superseded, and its pointer is kept as "Then at $\lambda_i=\zeta^i$ (see also (d) of the translation below),". (b), (c), (d): applied as given. | s02_setup.tex |
| A4-02 | ALREADY DONE by the first pass (C3.1): the proof of `thm:mod-var` has "(ii) This part is not used. ..." from `prop:mod-period`; `b_signs_2.tex` has "(5-14)]". Verified, not touched. | s06_model_8.tex, b_signs_2.tex |
| A4-03 | NOT APPLIED: duplicate of AMEND_1 C1 | - |
| A4-04 | (a) APPLIED in `conv:mod-kernel`; in the table of `b_signs_2.tex` ADAPTED: the cell read "$\omega_{h,n}=(-1)^nW^{\mathrm{EO}}_{h,n}=(-1)^n\omega^{\mathrm{ABDKS}}_{h,n}$" (a chain), now "$\omega_{h,n}=(-1)^n\omega^{\mathrm{ABDKS}}_{h,n}$; $\omega_{h,n}=(-1)^nW^{\mathrm{EO}}_{h,n}$ for $2h-2+n\geq0$". (b), (c), (d) APPLIED. | s06_model.tex, s06_model_7.tex, b_signs_2.tex |
| A4-05 | (a) APPLIED. (b) ADAPTED: AMEND_1 A1.5 had already replaced the false clause by "\Cref{sec:gau} uses no quoted input beyond (M), (G0) and (QRR), which enter through the rows of sec:set (flatness) and in lem:gau-colrows (the value);". That sentence is kept (flatness is part of (QRR)), and the pointers of A4-05 were added inside its two parentheses: "(flatness, def:set-rows (b), in thm:gau-R0 (iii)--(iv))" and "(the value; hence in cor:gau-rowsNN (iii) and cor:gau-colours)". The second change ("and its appendix \cref{sec:set-lrproof}") APPLIED as given. | s02_setup.tex, s01_intro_3.tex |
| A4-06 | ALREADY DONE in substance by A2-23 (same finding): the parenthesis already lists the other denominators and proves that none vanishes. The A4-06 version of the parenthesis was NOT substituted, so as not to undo A2-23; see doubt 1. | s02_setup.tex |
| A4-07 | ADAPTED: in `def:set-graphsum` and `lem:set-pushins` A2-25 had already put "with a perfect symmetric pairing $\eta$"; only " (its Gram matrix is invertible)" was added. `def:flh-outer`: applied as given. | s02_setup_4.tex, s02_setup_5.tex, s07_comparison.tex |
| A4-08 | APPLIED (inserted at the beginning of Step 2, with the parenthesis on the derivatives of $\Fm_h$ after it) | s06_model_10.tex |
| A4-09 | NOT APPLIED: duplicate of AMEND_1 C4-C6 | - |
| A4-10 | APPLIED | b_signs.tex |

## P2

| item | status | file |
|---|---|---|
| A4-11 | APPLIED (both places, inside the texts of A4-12) | s01_intro.tex, s02_setup.tex |
| A4-12 | APPLIED (three places) | s01_intro.tex, s02_setup.tex, s01_intro_3.tex |
| A4-13 | APPLIED (both places); the optional sentence "(d) follows from (b), (c) and the genus-zero topological recursion relation." added at the end of (G0)(d), after the sentence of A4-01 (c) | s02_setup.tex, s01_intro_3.tex |
| A4-14 | APPLIED | s01_intro.tex |
| A4-15 | APPLIED (three places) | s01_intro_2.tex, s09_gamma_3.tex, s09_gamma_4.tex |
| A4-16 | ADAPTED: A3-27 had inserted a clause in the sentence on the free energy ("For $h\geq 2$, on a curve where $\omega_{h,1}$ has zero residues ... (for the model curve, lem:mod-xfun), the free energy is defined by ..."). That clause is kept. Two edits: the sentence "This is the recursion of [DNOPS, (2-26)--(2-29)], with the kernel (2-28), to which the quoted theorem of thm:mod-tr refers." inserted after "Here $q$ is a point of the curve."; and "... $=y\,\dd x$. This is the convention of [DNOPS (2-26)--(2-29)] and of [Hock25 (2.1)]." replaced by "... $=y\,\dd x$; only this definition is used (it is that of [Hock25 (2.1)]; [DNOPS] defines no free energy)." Intro list: applied as given. `b_signs_2.tex` l.19 unchanged, as the item says. | s06_model.tex, s01_intro_3.tex |
| A4-17 | ADAPTED (intro): the old clause was followed by " (\cref{lem:mod-laplace,prop:mod-Rtr})."; the new text contains this reference itself, so the old trailing parenthesis was removed with the clause (otherwise it would occur twice). `s06_model_2.tex`: appended as given; the symbol $R_{\mathrm{TR}}$ is the one used in `prop:mod-Rtr`. | s01_intro_3.tex, s06_model_2.tex |
| A4-18 | APPLIED (both places) | s06_model_4.tex, s01_intro_3.tex |
| A4-19 | ADAPTED. "Def.~0.3 (iii);" in (Q): after AMEND_1 D8 the text is "Def.~0.3, in particular (iii);", which already pins the whole definition; not touched (AMEND_1 has priority). "(p.~20)" added after "Def.~2.2, with the convention stated after it". `s10_main_2.tex`: "Rem.~2.1 and (6);" was already there (AMEND_1 B4.5/D8). | s01_intro_3.tex |
| A4-20 | APPLIED (both changes); grep confirms that "(4-24)" and Theorem 2 of MR2276766 are cited nowhere else | s01_intro_3.tex |
| A4-21 | ADAPTED. Sentence of the introduction: ALREADY DONE by A3-53 ("uses no quoted input; the identity [DLMF 24.4.1] ... is derived from their generating function in the proof of lem:fsp-dlmf5118"), not replaced. Optional citation applied: "(the formal counterpart of \cite[5.11.8]{DLMF})" added at the end of the statement of `lem:fsp-dlmf5118`, after the reformulation of \eqref{eq:fsp-5118} (the item says "after \eqref{eq:fsp-5118}"; the equation is cited nowhere by \eqref, so the statement of the lemma is the place). `refs.bib`, DLMF annote: "Cited: Eqs. 24.2.3 (definition), 24.4.1 (derived in the text), 5.11.8 (comparison)" (the two annotes of the item merged). | a_formalsp_5.tex, refs.bib |
| A4-22 | APPLIED (bib note; `s06_model_7.tex` l.5 and l.20; (Q) in the introduction) | refs.bib, s06_model_7.tex, s01_intro_3.tex |
| A4-23 | APPLIED | s06_model_7.tex |
| A4-24 | SKIPPED: no text; a task for the author (journal numbering of Givental and Teleman). The draft notes stay. | - |

## P3

| item | status | file |
|---|---|---|
| A4-25 | APPLIED (both places) | s03_gauge_4.tex, s03_gauge_8.tex |
| A4-26 | APPLIED, minimal change (the parenthesis after "is summable, and" in (e2)); the fuller change $H^{\otimes j}$ not made | a3_uconv.tex |
| A4-27 | ADAPTED. First replacement APPLIED. Second: the old sentence no longer exists; A2-19 had replaced it by "with $X\coloneqq R'(0)^{-1}R'$ ... so the outer graph sum is the leg-free part of the unit-preserving action $X.(M.\TFT)$ on the colour space (lem:set-pushins (iii))", which removes the incorrect "$R'.(M.\TFT)$". Not touched. | s02_setup_8.tex |
| A4-28 | APPLIED (the two places named; the same words in `a_formalsp.tex` l.46 are not named by the item and are unchanged) | s01_intro_2.tex, b_signs_2.tex |
| A4-29 | APPLIED (after the sentence defining "regular at the conifold point") | s02_setup_6.tex |
| A4-30 | APPLIED (before "An element of $\Nring\otimes\bbQ(\zeta)$") | s02_setup_6.tex |
| A4-31 | APPLIED | s06_model_9.tex |
| A4-32 | APPLIED (both insertions, in the "Sources for the value" of (QRR) as rewritten by AMEND_1 A1.1) | s01_intro_3.tex |

## Pin table, rows 1-18

| row | status | file |
|---|---|---|
| 1 GJR18 Sec. 3.1, (9) | APPLIED: "\S3.1; \S3.2, (9)" (A4-11) | s01_intro.tex, s02_setup.tex |
| 2 ZZ08 as source of the specialisation | APPLIED (A4-12) | s01_intro.tex, s02_setup.tex, s01_intro_3.tex |
| 3 GJR18 Sec. 5.2, (14) | APPLIED (A4-13) | s02_setup.tex, s01_intro_3.tex |
| 4 Giv01 / MR2276766 "gives" (QRR) | APPLIED (A4-14) | s01_intro.tex |
| 5 GJR18 Lemma 7.7, (37) | APPLIED (A4-15), three places | s01_intro_2.tex, s09_gamma_3.tex, s09_gamma_4.tex |
| 6 DNOPS (2-26)-(2-29) as "the convention" | ADAPTED (A4-16) | s06_model.tex, s01_intro_3.tex |
| 7 DOSS Sec. 3.1, 3.5 "the prescription" | ADAPTED (A4-17): now (3-57), cf. (3-64) | s01_intro_3.tex, s06_model_2.tex |
| 8 EO App. A, (A-1), (A-5) | APPLIED (A4-18) | s06_model_4.tex, s01_intro_3.tex |
| 9 EO (5-15) | ALREADY DONE (first pass, A4-02). (5-15) remains only in the list "enter no proof" and in "(... \cite[(5-15)], which is not used)". | s06_model_8.tex, b_signs_2.tex |
| 10 EO (4-23)-(4-24); MR2276766 Theorem 2 | APPLIED (A4-20) | s01_intro_3.tex |
| 11 EO dictionary, range | APPLIED (A4-04 (a); table cell adapted) | s06_model.tex, b_signs_2.tex |
| 12 EO Thm 4.7, (4-30) | APPLIED (A4-04 (d)) | b_signs_2.tex |
| 13 ABDKS numbering in print | APPLIED (bib note, A4-22) | refs.bib |
| 14 ABDKS Thm 7.1 hypotheses | APPLIED (A4-04 (c), A4-22) | s06_model_7.tex |
| 15 PPZ Def. 0.3 (iii) | ALREADY DONE in (Q) by AMEND_1 D8 ("Def.~0.3, in particular (iii)"); not touched. See doubt 2 for `s10_main_2.tex`. | s01_intro_3.tex |
| 16 PPZ "p. 20" | APPLIED: "(p.~20)" named in (Q) (A4-19); the two citations in the body are unchanged, as the table says | s01_intro_3.tex |
| 17 PPZ Rem. 2.1 without (6) | ALREADY DONE in `s10_main_2.tex` (AMEND_1). One further occurrence found by grep and corrected: `s10_main.tex` l.71, `\cite[Rem.~2.1]{PPZ15}` -> `\cite[Rem.~2.1, (6)]{PPZ15}` (the form used in `s02_setup_4.tex`, `s02_setup_8.tex`). | s10_main.tex |
| 18 DLMF 24.4.1 as input; annote | sentence ALREADY DONE (A3-53); annote APPLIED (A4-21) | refs.bib |

Rows 19-26: no change asked (author's check); row 21 is with AMEND_1 A2.

## Remaining occurrences of old pins (grep of `main.tex` and `sections/` after the edits)

None in `main.tex`. In `sections/`:
* `s10_main_2.tex` l.36: "Def.~0.3~(iii);" in the list after Corollary B. Put there by A3-24 (opposite direction to row 15). Not undone; doubt 2.
* `s01_intro.tex` l.32: "the value at the orbifold point is \cite[Lemma~7.7, (37)]{GJR18}". Not one of the three places of A4-15 and it does not say "proved"; unchanged. The reading "right side of (37) as a series in $z$" is not repeated there.
* `s02_setup_8.tex` l.20: "\cite[Def.~0.3~(iii)]{PPZ15}" is a use of axiom (iii) itself; correct, unchanged.
* `s02_setup_4.tex` l.92, `s02_setup_5.tex` l.37: "\cite[p.~20]{PPZ15}"; row 16 says OK.
* `b_signs_2.tex` l.19: "equal to the kernel of \cite[(2-26)--(2-29)]{Dunin-Barkowski:2015caa} and of \cite[(2.1)]{Hock25}"; A4-16 says correct (it speaks of the kernel).
* `s01_intro_3.tex` l.23: "(5-15)" in the list of citations that enter no proof (now true).
* `s02_setup_4.tex` l.52, `b_signs_2.tex` l.91: \cite[Corollary~1]{MR2276766} in comparisons of the sign; not the "gives (QRR)" of row 4.
No occurrence remains of: "\S3.1, (9)", ZZ08 as the specialisation, "(14)]{GJR18}", "\S3.1, \S3.5", "(4-24)",
"Theorem~2]{MR2276766}", "(as in \cite[App.~A", "simple critical points of $x$", "5.11.1", "2.3.8".

## Checks made by script after the edits

No duplicate label; every `\cref`, `\Cref`, `\eqref` target exists; every citation key is in `refs.bib`; in every file
of `sections/` and in `main.tex` the counts of braces, parentheses, brackets, `$`, `\[`/`\]`, `\bigl`/`\bigr`,
`\Bigl`/`\Bigr`, `\begin`/`\end` are balanced exactly as in `archive_pre_amend_2026-11-18/`; each line has balanced
braces and an even number of `$`; braces of `refs.bib` balanced. One paragraph per line kept; `s02_setup.tex` gained
4 lines (the new display of A4-01 (a) with the sentence after it). Commands in the new text: all already in use
(`\Gst`, `\Ginf`, `\SPh`, `\ML`, `\Fm`, `\lser`, `\coloneqq`, `\mathsf`).

## Files modified by this pass (25)

`refs.bib`; sections: a3_uconv, a_formalsp_5, b_signs, b_signs_2, s01_intro, s01_intro_2, s01_intro_3, s02_setup,
s02_setup_4, s02_setup_5, s02_setup_6, s02_setup_8, s03_gauge_4, s03_gauge_8, s06_model, s06_model_2, s06_model_4,
s06_model_7, s06_model_9, s06_model_10, s07_comparison, s09_gamma_3, s09_gamma_4, s10_main (.tex).

## Doubts and points for the author (nothing was changed on their account)

1. A4-06 against A2-23. The parenthesis of `def:set-fq` is the A2-23 text (weights $(a\lambda_i+b\lambda_j)/d-\lambda_k$,
   node-smoothing weights $(\lambda_i-\lambda_j)/d_1+(\lambda_i-\lambda_k)/d_2$, classes $(\lambda_i-\lambda_j)/d-\psi$;
   "a proper convex combination of two distinct fifth roots of unity has modulus $<1$"). The A4-06 text differs: it
   opens with "No denominator of the localisation formula vanishes there.", treats $m\in\{i,j\}$ explicitly, and has
   no node-smoothing weights. The author may prefer one of the two; I did not merge them.
2. Def. 0.3 of [PPZ15]. Row 15 / A4-19 want "Def. 0.3" (all of (i)-(iii) are used); A3-24 put "Def.~0.3~(iii)" into
   the list after Corollary B (`s10_main_2.tex` l.36), and (Q) reads "Def.~0.3, in particular (iii)". The two lists
   therefore still differ in wording, and the one after Corollary B is the narrower.
3. A4-01. (i) The statement of `thm:set-mirror` now cites `thm:set-G0` (b), a quoted theorem stated after it (forward
   reference between two quoted inputs). (ii) In (G0)(d) the degree-0 term is written $\eta(xy,w)$ with the classical
   product $xy$; typed as given. (iii) The optional sentence "(d) follows from (b), (c) and the genus-zero topological
   recursion relation" stands inside a theorem marked "quoted". (iv) I rechecked the formula for $J(\tau H,z)$: the
   degree-0 part of the sum is $z\re^{\tau H/z}$ and the descendant divisor equation gives the factor $\re^{\tau H/z}$
   on each degree; no error found.
4. A4-16, introduction: "the recursion [DNOPS, (2-26)--(2-29)], in which [DOSS, Theorem 4.1] is recalled" is typed as
   given; (Q) in the same subsection says that the theorem is recalled in [DNOPS, \S2.6, (2-30)--(2-32), (2-34),
   (2-35)]. The two descriptions should be reconciled. The kernel number (2-28) is a single-referee pin.
5. A4-17. In `s06_model_2.tex` the new sentence speaks of $\dd\xi^i$ and $R_{\mathrm{TR}}$ at l.23, before they are
   introduced (in `prop:mod-Rtr`, l.43). In the introduction the symbols $\mathsf G_{-2}$, $\vartheta$, $w$ are used
   without explanation (they are those of `lem:mod-laplace`). Rechecked: $\int w^{-2}\re^{-\vartheta w^2/2}\dd w=-\sqrt{2\pi\vartheta}$
   by one integration by parts, and $\mathsf G_{-2}=-\vartheta$ in \eqref{eq:mod-laplace}. Pins (3-57), (3-64): single referee.
6. A4-23. The sentence now reads "\cite[(7.6), (7.7)]{ABDKS25}, whose closed form $W_{i,j}=\dots$, printed in
   \cite[(7.7)]{ABDKS25}, expands as ..."; "whose" has no clear antecedent, and the display above it defines $W_{i,j}$
   by an exponential. Rechecked to order $\hbar^2$: the exponential formula, the closed form and the series agree.
7. A4-22. "all at finite points" and "the swapped data $(\bbP^1,p,\Wv)$, whose recursion is trivial because $\dd p$
   has no zeros": typed as given; Theorem 1.14 and the preamble of Sec. 7.1 of [ABDKS25] are single-referee pins.
8. A4-21. The new citation [DLMF, 5.11.8] rests on one referee's term-by-term comparison to order 9; it is a comparison,
   not an input. The annote of the bib entry was reworded by me (merge of the two versions of the item).
9. A4-08. The claim that the derivatives of $\Fm_h$ which occur are in $\bbQ[v,v^{-1}]$, and the reference
   `lem:set-modhom` for the shape of the correlators, are typed as given; not rechecked.
10. A4-10. "$B_1=-\frac12$ occurs only inside the proof of `lem:fsp-stirfe`": typed as given, not rechecked by grep
    of every use of $B_1$.
11. A4-12. `s02_setup.tex` now attributes the specialisation to [LP18, Sec. 0.2] and [LP20, Sec. 0.2]; that ZZ08 has
    no torus weights rests on one referee.
12. A4-30. The target $\bbQ(\zeta)\lser s$ of the extended $\Cexp_{\alpha,\beta}$ is typed as given.
13. A4-05 (b). The merged sentence says that flatness is used "in `thm:gau-R0` (iii)--(iv)"; this is the pointer of
    A4-05 moved under (QRR). The sentence earlier in the same paragraph lists flatness "in `prop:set-cyclicLR` (iii)
    and `cor:set-cyclicPF5` (b)"; both are true if `thm:gau-R0` uses the rows through these, which I did not recheck.
14. A4-26. The parenthesis stands between "is summable, and" and a displayed formula; it reads a little heavily.
