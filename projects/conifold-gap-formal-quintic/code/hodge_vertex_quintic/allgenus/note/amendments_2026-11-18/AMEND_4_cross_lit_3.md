# AMEND_4: cross-cutting and literature reports (part 3 of 3): P2 continued, P3, rejections, ledger

Conventions as in `AMEND_4_cross_lit.md`. "Old" texts were read in the current `sections/*.tex` and `refs.bib`.

## 4. P2 items (continued)

### A4-19 (P2) The PPZ items are listed differently in three places [X5-F4, third bullet; X5-F10]
Checked: `s01_intro_3.tex` l.9 "Def.~0.3 (iii);" while `s02_setup_8.tex` l.6 uses "\cite[Def.~0.3~(i)--(iii)]{PPZ15}"
and `s10_main_2.tex` l.36 lists "Def.~0.3; (5); Rem.~2.1; Defs.~2.2 and~2.5"; "p.~20" is cited at `s02_setup_4.tex`
l.79 and `s02_setup_5.tex` l.37 and is the convention after Def. 2.2 (L1: pp. 20-21). Agree.
`s01_intro_3.tex` l.9: `Def.~0.3 (iii);` -> `Def.~0.3;` and
`Def.~2.2, with the convention stated after it;` -> `Def.~2.2, with the convention stated after it (p.~20);`.
`s10_main_2.tex` l.36: `Rem.~2.1;` -> `Rem.~2.1 and (6);`.
(E6 rewrites the end of l.9, "Props. 2.4, 2.8 and 2.9 enter only through ...", and adds Props 2.3, 2.7, Lemma 2.11.)

### A4-20 (P2) The list of citations that "enter no proof" names two items the body no longer cites [X5-F10]
Checked: `s01_intro_3.tex` l.23 "\cite[Theorem~4.4, (4-23)--(4-24); ...]{Eynard:2007kz}, \cite[(26)]{LP20},
\cite[Theorem~2]{MR2276766},". grep: the body cites "(4-23)" only (`s06_model_4.tex` l.16) and never Theorem 2 of
MR2276766. Agree. `(4-23)--(4-24);` -> `(4-23);` and delete `\cite[Theorem~2]{MR2276766}, `.

### A4-21 (P2) Theorem Gamma (b) uses no quoted input; the DLMF annotation is stale [X5-F7 = L6-m6]
Checked: `s01_intro_3.tex` l.21 "\Cref{thm:int-gamma} uses no quoted input, except the identity \cite[24.4.1]{DLMF} of
Bernoulli polynomials in part (b) (\cref{thm:gam-ginf})."; the proof of lem:fsp-dlmf5118 (`a_formalsp_5.tex` l.69)
derives the identity from the generating function ("the generating functions of \cite[24.2.3]{DLMF} at $\theta+1$ and
$\theta$ differ by $\xi\re^{\theta\xi}$, so ..."). `refs.bib`, entry DLMF: annote "Cited: Eqs. 24.2.3, 24.4.1; in a
remark 5.11.1, 5.11.8, 2.3.8"; grep finds none of the last three in `sections/`. Agree. Old: the sentence quoted. New:
```latex
\Cref{thm:int-gamma} uses no quoted input: in part (b) (\cref{thm:gam-ginf}) the Bernoulli polynomials are defined by the generating function \cite[24.2.3]{DLMF} (\cref{def:fsp-stirling}), and the difference identity \cite[24.4.1]{DLMF} is derived from it in the proof of \cref{lem:fsp-dlmf5118}.
```
`refs.bib`, DLMF: `annote = {Cited: Eqs. 24.2.3 (definition), 24.4.1 (derived in the text)}`. Optionally, after
\eqref{eq:fsp-5118}, add "(the formal counterpart of \cite[5.11.8]{DLMF})"; L6 checked the agreement term by term to
order 9. Then the annote reads "Cited: Eqs. 24.2.3, 24.4.1, 5.11.8".

### A4-22 (P2) ABDKS: the printed numbering differs; the hypotheses of Theorem 7.1 are paraphrased incompletely [L4-F1, L4-F2]
Checked: `refs.bib` ABDKS25: `note = {arXiv:2212.00320v3; numbers cited from this version}` with the JEMS DOI;
`s06_model_7.tex` l.5 "for a rational function $x(z)$ on $\bbP^1$ such that $\dd x$ has only simple zeros";
l.20 "The standing assumptions of \cite[\S1.1]{ABDKS25} (simple critical points of $x$, $\dd y\neq 0$ there) hold".
L4 read v3 and the JEMS text. Agree (single referee for the JEMS numbers).
`refs.bib`, ABDKS25, field note:
`note = {arXiv:2212.00320v3; numbers cited from this version. In the journal version (7.6) is (7.4), (7.9) is (7.5), (7.7) and (7.10) are unnumbered displays, and (1.13) is (1.7); Theorems 1.14 and 7.1 and (1.1) keep their numbers}`.
`s06_model_7.tex` l.5: old `such that $\dd x$ has only simple zeros,` new
`such that $\dd x$ has only simple zeros, all at finite points (the theorem is printed for an arbitrary rational $x$; the condition is that of the preamble of \cite[\S7.1]{ABDKS25} and of \cite[Theorem~1.14]{ABDKS25}, from which the theorem is derived),`.
`s06_model_7.tex` l.20: old `The standing assumptions of \cite[\S1.1]{ABDKS25} (simple critical points of $x$, $\dd y\neq 0$ there) hold for` new
`The standing assumptions of \cite[\S1.1]{ABDKS25} (a compact curve with its Bergman kernel; the zeros of $\dd x$ simple; $y$ regular with $\dd y\neq0$ there) and the hypotheses of \cite[Theorem~1.14]{ABDKS25} for the swapped data $(\bbP^1,p,\Wv)$, whose recursion is trivial because $\dd p$ has no zeros, hold for`.
`s01_intro_3.tex` l.9: `with the kernel (1.1)]{ABDKS25}` -> `with the kernel (1.1) and the assumptions of \S1.1]{ABDKS25}`.

### A4-23 (P2) The expansion of $W_{i,j}$ is stated without source; its closed form is part of (7.7) [X5-F6, c9]
Checked: `s06_model_7.tex` l.15 "\cite[(7.6), (7.7)]{ABDKS25}, and $W_{i,j}=\sum_{k\geq 0}(\hbar^2/4)^k(w_i+w_j)^{2k}w_iw_j(z_i-z_j)^{-2k-2}$".
L4 (2.3): (7.7) itself prints the closed form; X5 and L4 both verified the expansion. Agree. Old: `and $W_{i,j}=\sum_{k\geq 0}` New:
`whose closed form $W_{i,j}=w_iw_j\big/\bigl((z_i-z_j)^2-\hbar^2(w_i+w_j)^2/4\bigr)$, printed in \cite[(7.7)]{ABDKS25}, expands as $W_{i,j}=\sum_{k\geq 0}`.

### A4-24 (P2) Statement numbers of (T) and (M) were not compared with the journal versions [X5-F11]
Checked: draft notes at `s01_intro_3.tex` l.11 ("numbers from arXiv v2; IMRN numbering not compared") and in
`refs.bib` MR2917177 ("journal numbering not compared"). The sentence "Statement numbers refer to the arXiv versions
named in the bibliography" (l.7) makes the citations consistent, so no change is forced. Before submission: compare
Theorem 11.6, Prop. 11.4, Cor. 11.7 and Section 9 of Givental with IMRN 1996, and Theorem (2) in (1.15) and (after E4)
Props 7.7, 8.2, 8.4, 8.7 of Teleman with Invent. Math. 188; then replace each draft note by the journal numbers in the
bib note. The audit could not do this (journal versions not available to L5, L6).

## 5. P3 items (optional)

### A4-25 (P3) Two mutual citations are resolved only by the order of the parts [X1-F3]
Checked: `s03_gauge_4.tex` l.22 "Part (c) is \cref{prop:gau-graded}."; `s03_gauge_8.tex` l.30 "this is
\cref{cor:gau-rowsNN}~(i) for colour $0$". Both well founded (X1 (a2), (b)-17).
l.22: -> `Part (c) is \cref{prop:gau-graded}, proved below from (a), (b) and (d).`
l.30: -> `this is \cref{cor:gau-rowsNN}~(i) for colour $0$ (whose proof uses (iii), not (iv))`.

### A4-26 (P3) $H^j$ for "j insertions of H" clashes with the flat class $H^j$ [X2-m2]
Checked: `a3_uconv.tex` l.26-28 ((e2), (e3): `\langle(y^\circ)^lH^j\rangle_g`, `\langle H^j\rangle_g`,
`\langle(y^\circ)^0H^0\rangle_g`) and l.53-61. Agree; this is where $c_g$ is produced. Minimal change, l.26, after
"is summable, and": insert `(here and below $H^j$ under a bracket denotes $j$ insertions of $H=H^1$, not the flat class $H^j$)`.
Fuller change: write `H^{\otimes j}` in the brackets of l.26-28 and l.53-61.

### A4-27 (P3) def:set-outer: "sum over all colourings" and "$R'.(M.\TFT)$" [X2-m3; L1-M3]
Checked: `s02_setup_8.tex` l.50: "the vertex tensors are the full tensors of $M.\TFT$, so $T_\Gamma$ contains the sum over
all colourings of $\Gamma$" ($M.\TFT$ is not colour-diagonal); and "the outer graph sum is the leg-free part of
$R'.(M.\TFT)$ in the sense of \cref{def:flh-outer}" ($R'(0)=\Psi_{\mathrm Q}\neq1$ maps the colour space to $\bbH$, so
this is not the action of PPZ Def. 2.13, and def:flh-outer defines contractions only). Agree with both.
Old: `so $T_\Gamma$ contains the sum over all colourings of $\Gamma$.` New:
`with all their components in the colour basis ($M.\TFT$ is not diagonal in the colours, so there is no separate sum over colourings).`
Old: `and the outer graph sum is the leg-free part of $R'.(M.\TFT)$ in the sense of \cref{def:flh-outer}.` New:
`and the outer graph sum is the sum of the contractions of \cref{def:flh-outer}~(iii); it equals the leg-free part of $(\Psi_{\mathrm{cf}})_*\bigl(X.(M.\TFT)\bigr)$ with $X\coloneqq\Psi_{\mathrm{cf}}^{-1}R'$ (\cref{cor:mai-gap}, Step~2).`
(E8 edits other sentences of the same definition.)

### A4-28 (P3) "normalised stationary phase" names $\SPn$ and $\SPh$ [X3-m4]
Checked: `a_formalsp.tex` l.38-42 calls $\SPn$ "the normalised stationary phase"; `s01_intro_2.tex` l.53 "the normalised
stationary-phase series at $\sigma=L$" (with $c_0=1$, so $\SPh$); `b_signs_2.tex` l.65 "Only normalised series are
canonical". Agree; each place is unambiguous in context.
`s01_intro_2.tex` l.53: `the normalised stationary-phase series at` -> `the normalised stationary-phase series, divided by its constant term ($\SPh$ of \cref{def:fsp-wick}), at`.
`b_signs_2.tex` l.65: `Only normalised series are canonical` -> `Only the series $\SPh$ are canonical`.

### A4-29 (P3) "regular" has three meanings [X4-m6]
Checked: `s02_setup_6.tex` l.20 ("regular at the conifold point": expansion in $\bbQ\fps u$); `s02_setup_7.tex` l.56
("A correlator is \emph{regular} if $e\geq0$"); `s03_gauge_2.tex` l.9 ("A \emph{regular gauge} is ... with $P_0$
invertible"). Agree; no statement mixes them. `s02_setup_6.tex` l.20, after the sentence defining "regular at the
conifold point", add: `(This includes the rationality of the coefficients. A \emph{regular gauge} (\cref{def:gau-onepole}) and a \emph{regular} model correlator (\cref{lem:set-modhom}~(c)) are different notions.)`

### A4-30 (P3) The $\bbQ(\zeta)$-linear extension of $\Cexp_{\alpha,\beta}$ is named only in the proof of cor:mai-gap [X4-m7, second point]
Checked: `s10_main_2.tex` l.11 "Write $\Cexp$ for $\Cexp_{\alpha,\beta}$, extended $\bbQ(\zeta)$-linearly"; def:set-confexp
(`s02_setup_6.tex` l.20) already speaks of the conifold expansion of elements of $\Nring\otimes\bbQ(\zeta)$. Agree.
`s02_setup_6.tex` l.20, before "An element of $\Nring\otimes\bbQ(\zeta)$": insert
`The $\bbQ(\zeta)$-linear extension $\Nring\otimes\bbQ(\zeta)\to\bbQ(\zeta)\lser s$ is denoted by the same letter.`

### A4-31 (P3) $\mathrm{reg}_{00}=\frac16\{t;\chi\}(0)$ is stated without proof [X5-F6, c11]
Checked: `s06_model_9.tex` l.2. I recomputed it: for $t=\chi+a\chi^2+b\chi^3+\cdots$,
$t'(\chi_1)t'(\chi_2)-\bigl((t(\chi_1)-t(\chi_2))/(\chi_1-\chi_2)\bigr)^2=(b-a^2)(\chi_1-\chi_2)^2+O(3)$. Agree. After
"is the Schwarzian derivative of $t$ as a function of $\chi$" add:
`(both sides are $t_3/t_1-(t_2/t_1)^2$ for $t=t_1\chi+t_2\chi^2+t_3\chi^3+\cdots$, by expanding $t'(\chi_1)t'(\chi_2)/(t(\chi_1)-t(\chi_2))^2$)`.

### A4-32 (P3) (QRR), "Sources": say what each source gives [L5-F8, second location]
Checked: `s01_intro_3.tex` l.10 "Sources: \cite[(17), ...]{GJR18}; \cite[(19)--(20), ...]{GJR17}; \cite[\S2.2, Theorem~1 and
the Remark after it]{Giv01}; \cite[Theorem~1, Corollary~1]{MR2276766}." Accurate as a list. Optional: after
`at $t=0$]{GJR17}` insert `, which state it for this theory`, and after the last item add
`, which give the tangent factor and the Euler-class twist separately`.

## 6. Rejected (R)

- R1 [L1-M4] "as in Givental and Teleman \cite[p.~20]{PPZ15}" (`s02_setup_4.tex` l.79, `b_signs.tex` l.14). The
  quotation of PPZ is accurate (L1), and the note does not rest on it: the identification of the action is (T), see E4.
- R2 [L5-F7] GJR18 Theorem 4.1 mixes $\Omega^\lambda$ and $R_\tau$. L5 itself asks for no change; thm:set-GT states the
  coherent reading, and E4 restates (T).
- R3 [L6, remark in 1.1] Cite a complete account of the mirror theorem (LLY97 is in `refs.bib`, uncited). The entry's
  own annote says it is cited only if a statement is read and pinned; no referee read it. Author's choice.
- R4 [X4-m3, side remark] Corollary B for all complex $(\alpha,\beta)$ by Zariski density. A new statement, checked
  by one referee only and outside the scope of the note (frames are rational by def:set-frames).
- R5 [X4-m7, first point] "$\Nring\lser z\otimes\bbQ(\zeta)$" in prop:gau-quintic (i) (`s03_gauge_2.tex` l.20). Not an
  error: the statement about the forward relations holds for Laurent series in $z$, and is more general than for rows.
- R6 [X4-m7, third point] def:set-rows, "$A_n=L^{4n+1}\mathrm{ph}_n(L^{-5})$" (`s02_setup_4.tex` l.19). It is a
  definition in $\bbQ(\zeta)\fps q$ with $u=L^{-5}$ and is consistent; the Laurent-polynomial statement is correctly
  placed in thm:gau-R0 (iii).
- R7 [X5-F6, c10] Lagrange-Buermann in prop:mod-period Step 4 (`s06_model_7.tex` l.71). The text states the formula
  with its hypotheses and its proof in a parenthesis ("change of variables and one integration by parts"); X5 checked
  the result for $h\leq6$. Nothing to add.
- R8 [L6-m6, part] Renaming the label lem:fsp-dlmf5118. Labels are not printed; renaming risks broken references.
  The optional citation of DLMF 5.11.8 in A4-21 meets the point.

Disagreements with the audit, in short: R4-R8 above; in A4-02 I choose one of the three repairs (prove (ii) from
prop:mod-period); in A4-16 the repairs of X5 and L3 contradict each other and both are needed in part; X5's closure of
F1 (embedding the rings in $\bbC$) is not adopted (E6 already says so).

## 7. Ledger: the inputs listed in `s01_intro_3.tex` (l.7-23) and `s10_main.tex` (l.31), `s10_main_2.tex` (l.36), against what the proofs use (X5)

| used in a proof, or in a statement a proof uses | where | in the lists? | amendment |
|---|---|---|---|
| def:set-rows (b): the columns of $R$ solve the quantum differential equation | s02_setup_4 l.24; s03_gauge_7; s03_gauge_8; s08_connection_8 | no (called a convention) | see E4 |
| all-genus divisor equation; all-genus axioms with flat unit | s02_setup_5 l.42 | no ((G0) is genus 0) | see E4 |
| descendant string and divisor equations in genus 0; meaning of $t=T(q)H$, $\Qn=1$ | s02_setup l.80, l.87 (d); a1_cyclicLR l.13 | partly ((G0)(b) is primary) | A4-01 |
| (QRR) inside Section 3 | s03_gauge_7 l.62 | listed, but l.17 denies it | A4-05 (b) |
| (M), (G0) in the appendix of Section 2 | a1_cyclicLR | "only in sec:set" | A4-05 (b) |
| PPZ over rings other than Q; Props 2.3, 2.7, Lemma 2.11 | s02_setup_4 l.79; s10_main_2 l.27 | no | see E6 |
| PPZ Def. 0.3 (i)-(ii); "p. 20"; "(6)" in s10_main_2 | s02_setup_8 l.6; s02_setup_4 l.79; s10_main_2 l.36 | inconsistent | A4-19 |
| DNOPS (2-26)-(2-29) (the recursion of the quoted theorem) | s06_model l.20 | no | A4-16 |
| DNOPS (2-33), Prop. 4.3 (1), (8-1), Sec. 2.3 | s06_model_3 l.23, l.28 | listed, wrong | see E5 |
| DOSS Sec. 3.1, 3.5 (reading of the Laplace integrals) | s06_model_2 l.23 | listed as "prescription" | A4-17 |
| EO (5-15) | s06_model_8 l.39 | listed as "enters no proof" | A4-02 |
| EO App. A, (A-1), (A-5) | s06_model_4 l.12 | no | A4-18 |
| ABDKS Sec. 1.1, Theorem 1.14 (hypotheses of Theorem 7.1) | s06_model_7 l.5, l.20 | no | A4-22 |
| DLMF 24.2.3 (definition); 24.4.1 | a_formalsp_5 l.9, l.69 | 24.4.1 listed as an input, not needed | A4-21 |
| listed but not cited: EO (4-24); MR2276766 Theorem 2 | s01_intro_3 l.23 | residue | A4-20 |
| computations not in sec:evi behind three remarks | s01_intro_2 l.25; s10_main_2 l.46, l.51 | "no finite computation enters a proof" stays true | A4-09 |

Everything else in X5's ledger (Section 2, about 40 rows) matches the lists: every other citation inside a proof of
Theorems A, C, Gamma (a) or Corollary B is a declared item. X5's unproved assertions c1-c20: c1 is E6; c3, c4 are E4 and
A4-01; c9, c11 are A4-23, A4-31; the others were verified by X5 and need no change. X1: no circularity; nothing to amend
beyond A4-05 (a) and A4-25. X2: all interfaces and the constant $c_g$ match; only A4-08, A4-26, A4-27.
X3: no sign or normalisation error; only A4-10, A4-28 (and E5, E8).

## 8. What I could not verify

- I did not open any external source: every "corrected pin" rests on what the literature referees report having read
  in the arXiv versions. Single-referee pins: GJR18 Cor. 7.6 (A4-15); DNOPS kernel number (2-28) (A4-16); DOSS (3-57),
  (3-64) (A4-17); ABDKS Theorem 1.14, preamble of Sec. 7.1 and the JEMS numbers (A4-22); EO (4-7) (A4-04).
  Two referees agree on GJR18 (9) (A4-11).
- Rows 19-26 of the table of pins: journal numbering of Teleman and Givental; Hock25; LP20 (17)-(20), (26); GR19;
  MR1901075; MR2553377; HKQ09; Br25; CvdLPS23; CGYZ26. None was read by the audit.
- In A4-17 the symbol of the TR matrix in the new sentence must be the one used at `s06_model_2.tex` l.43; in A4-27 the
  letters $\Psi_{\mathrm{cf}}$, $X$ are those of lem:mai-Psi and cor:mai-gap. Nothing was compiled.
- A4-06: the list of denominators follows L6 and the standard localisation formula; I checked the non-vanishing, not the
  completeness of the list against a written-out formula in the sources.
