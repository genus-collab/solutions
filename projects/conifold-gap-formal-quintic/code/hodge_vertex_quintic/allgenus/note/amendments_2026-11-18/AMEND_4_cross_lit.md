# AMEND_4: amendments from the cross-cutting and literature reports (part 1 of 3)

Sources: `audit_2026-11-18/reports/X1.md` ... `X5.md` (circularity, interfaces, conventions, rings, ledger) and
`L1.md` ... `L6.md` (PPZ, Eynard-Orantin, DOSS/DNOPS, ABDKS, QRR/Teleman/G0, mirror theorem and others).
Paths are relative to the note folder. Line numbers are those of the current `sections/*.tex` (one paragraph per line);
every "old" text below was read in the current TeX. Nothing in the note was edited.

Parts: this file = scope, counts, table of pins, P1 items; `_2.md` = P2 items; `_3.md` = P3 items, rejections,
ledger of inputs (X5) against `s01_intro_3.tex` / `s10_main.tex`, and what could not be verified.

Classes. P1 = a statement or proof step is wrong or incomplete as written. P2 = wrong or imprecise citation, pin or
cross-reference. P3 = wording or notation, optional. R = rejected, with the reason.

## 0. Findings left to the other agent (not treated here)

| finding | where ruled |
|---|---|
| X1-F4, X5-F1, L1-M1 (PPZ over rings other than Q; Props 2.3, 2.7, Lemma 2.11 inside Prop. 2.12) | see E6 |
| X2-m1, X4-m5 (the TFT in Theorem A is the expanded one); X3-m1 (hat-Delta); X3-m3 (eta^{-1} in the outer edge) | see E8 |
| X2-m4, X4-m2, X4-m3, X4-m4, X5-F8 (what Corollary B states; injectivity of ev; CGYZ comparison) | see E7 |
| X3-m5, L3-F1, L3-F3, L3-F4, L1-M5 ((2-31) index reading, (8-1) tuple, 1/abs(Aut), Prop. 4.3, dilaton sign) | see E5 |
| L5-F1, L5-F2, X5-F5 first two bullets (flatness (b); Teleman's hypotheses; all-genus divisor equation; log-free ring) | see E4 |
| X4-m1 (Galois descent cites lem:gau-colrows instead of thm:con-glue) | README, routine list |

Three items of the README list "false as written, used by no proof" come from my reports and have no E-ruling; they are
treated here (A4-02, A4-03, A4-09). If the other agent also covers them, keep one version.

## 1. Counts

| class | number | items |
|---|---|---|
| P1 | 10 | A4-01 ... A4-10 |
| P2 | 14 | A4-11 ... A4-24 |
| P3 | 8 | A4-25 ... A4-32 |
| R | 8 | R1 ... R8 |

Merged duplicates: X1-F2 = X5-F2; X5-F3 = L2-m1; X3-m2 = L6-m5; L5-F4 = L6-m2; X5-F7 = L6-m6 = part of X5-F10;
X5-F4 (first bullet) and L3-F5 (opposite repairs, reconciled in A4-16); X5-F4 (last bullet) = L3-F2.

## 2. Table of pins (statement numbers, equation numbers, versions)

Status: OK = located by a literature referee in the version named in `refs.bib`; WRONG = the number or attribution is
wrong; IMPRECISE = right place, wrong description or missing range; UNVERIFIED = no referee compared it.
Only rows with a status other than plain OK are listed in full; the last rows record what was confirmed in bulk.

| # | citation in the note | where | status | corrected pin |
|---|---|---|---|---|
| 1 | GJR18, Sec. 3.1, (9) | s01_intro l.10; s02_setup l.24 | WRONG | (9) is the first display of Sec. 3.2; the cycle at lambda_i = zeta^i lambda is the last (unnumbered) display of Sec. 3.1. Pin: "Sec. 3.1; Sec. 3.2, (9)" (A4-11) |
| 2 | ZZ08 as the source of the specialisation lambda_i = zeta^i | s01_intro l.10; s02_setup l.24; s01_intro_3 l.15 | WRONG | arXiv:0710.0889v1 has no torus weights; it studies the series only. Pin the specialisation to LP18 Sec. 0.2, LP20 Sec. 0.2, GJR18 Sec. 3.1; cite ZZ08 for the series of lem:set-periods (A4-12) |
| 3 | GJR18, Sec. 5.2, (14) for (G0)(d) | s02_setup l.102; s01_intro_3 l.12 | IMPRECISE | (14) is the q-direction equation for S; the quantum differential equation is the unnumbered display before (14) (A4-13) |
| 4 | Giv01 Sec. 2.2 Thm 1; MR2276766 Cor. 1 "gives" (QRR) | s01_intro l.12 | IMPRECISE | neither states it for this theory; the specialised statement is GJR18 (17), GJR17 (19)-(20) (A4-14) |
| 5 | GJR18, Lemma 7.7, (37), "the value proved there" | s01_intro_2 l.68; s09_gamma_3 l.28; s09_gamma_4 l.72 | IMPRECISE | the formula is the right side of (37) read as a series in z; the coefficients r are introduced in GJR18 Cor. 7.6 through (-z)^l (A4-15) |
| 6 | DNOPS (2-26)-(2-29) as "the convention" of the free energy | s06_model l.20 | WRONG | (2-26)-(2-29) define the recursion only (kernel (2-28)); no free energy there. Also not in the list of inputs although the kernel identification is used (A4-16) |
| 7 | DOSS Sec. 3.1, 3.5 "the prescription" | s01_intro_3 l.15 | IMPRECISE | DOSS states no rule for polar integrands; consistent formulas are DOSS (3-57), (3-64) (A4-17) |
| 8 | EO App. A, (A-1), (A-5) | s06_model_4 l.12 (inside a proof) | IMPRECISE | not in the lists of s01_intro_3; make it "cf." and add to the non-proof list (A4-18) |
| 9 | EO (5-15) | s06_model_8 l.39 (proof); b_signs_2 l.24 | WRONG (ledger) | declared "enters no proof" at s01_intro_3 l.23 but used for thm:mod-var (ii) (A4-02) |
| 10 | EO (4-23)-(4-24); MR2276766 Theorem 2 | s01_intro_3 l.23 | WRONG (residue) | the body cites (4-23) only and never Theorem 2 of MR2276766 (A4-20) |
| 11 | EO dictionary omega_{h,n} = (-1)^n W_n^(h) | s06_model l.21; b_signs_2 l.20 | IMPRECISE | holds for 2h-2+n >= 0; W_1^(0) = 0 is EO (4-7) (A4-04) |
| 12 | EO Thm 4.7, (4-30) | b_signs_2 l.23 | IMPRECISE | range 2h-2+n > 0 missing; (4-30) has residues at the p_i and a delta term (A4-04) |
| 13 | ABDKS (7.6), (7.7), (7.9), (7.10) | s01_intro_3 l.9; s06_model_7 l.10-15 | OK in v3; differ in print | JEMS: (7.6) = (7.4), (7.9) = (7.5), (7.7) and (7.10) unnumbered, (1.13) = (1.7); Thm 7.1, Thm 1.14, (1.1) unchanged (A4-22) |
| 14 | ABDKS Thm 7.1 "for x with only simple zeros of dx" | s06_model_7 l.5, l.20 | IMPRECISE | printed for "arbitrary rational x"; the condition is in the preamble of Sec. 7.1 and in Thm 1.14; (7.9) fails at (g,m) = (0,1) (A4-04, A4-22) |
| 15 | PPZ Def. 0.3 (iii) | s01_intro_3 l.9 | IMPRECISE | (i)-(iii) are used (s02_setup_8 l.6); pin "Def. 0.3" (A4-19) |
| 16 | PPZ "p. 20" | s02_setup_4 l.79; s02_setup_5 l.37 | OK, not itemised | it is the convention after Def. 2.2 (pp. 20-21), already in (Q); say so (A4-19) |
| 17 | PPZ Rem. 2.1 without (6) | s10_main_2 l.36 | IMPRECISE | "Rem. 2.1 and (6)", as in s01_intro_3 l.9 (A4-19) |
| 18 | DLMF 24.4.1 as a quoted input; annote "5.11.1, 5.11.8, 2.3.8" | s01_intro_3 l.21; refs.bib DLMF | WRONG (ledger) | 24.4.1 is derived in the proof of lem:fsp-dlmf5118; 24.2.3 is the definition; the three numbers of the annote are not cited (A4-21) |
| 19 | Teleman, Theorem (2) in (1.15) | s01_intro_3 l.13; s02_setup_5 l.42 | OK in arXiv v3; content: see E4 | exists in v3; journal numbering not compared (UNVERIFIED); the later rendering renumbers (v3 Props 7.7, 8.2, 8.4, 8.7 = later 7.5, 8.1, 8.3, 8.5) |
| 20 | Givental, Thm 11.6, Sec. 9 pp. 23-27, Prop. 11.4, Cor. 11.7 | s01_intro_3 l.11; s02_setup l.85 | OK in arXiv v2 | IMRN numbering not compared (UNVERIFIED); keep the draft note until compared |
| 21 | DNOPS Prop. 4.3 (1); (8-1) with a = (0,0,0,0,4v); (2-33); Sec. 2.3 | s01_intro_3 l.9; s06_model_3 l.28; b_signs l.14, l.52, l.58 | WRONG | see E5: Prop. 4.3 (2); a_1 = ... = a_4 = 0, a_5 = 4v; (2-33) misprinted (use DOSS (4-1)); Sec. 2.3 omits 1/abs(Aut) |
| 22 | Hock25 (2.1), (1.1), Cor. 1.2, Sec. 3.2 | s06_model l.20; s06_model_7 l.87 | UNVERIFIED | not read by any referee; enters no proof |
| 23 | LP20 (26); Sec. 2.4.3 (17); Sec. 2.4.4 (18), (19)-(20) | s02_setup_4 l.52; a1_cyclicLR l.107; b_signs_2 l.89, l.92; s01_intro_3 l.23 | UNVERIFIED | not read (L5 says so); note that b_signs_2 l.92 cites "(18)" while s01_intro_3 l.23 lists "(19)-(20)": check which is meant |
| 24 | GR19 Props 6.1-6.2; MR1901075 Sec. 10, 10.4 Lemma, 9.1 Theorem; MR2553377 Lemma 4.20; MR2276766 Sec. 7 (11); HKQ09 (1.2); Br25 (1.5); CvdLPS23 | a1_cyclicLR l.107; s09_gamma_4 l.72; b_signs_2 l.88-94; b_signs l.58 | UNVERIFIED | comparisons only; enter no proof |
| 25 | CGYZ26 Conj. 1.5, Thm 2, Defs 1.1, 1.2, 1.8, Rem. 1.9 | s10_main_2 l.40; s01_intro_3 l.30-32; b_signs_2 | UNVERIFIED | the draft note at s10_main_2 l.41 says they were read in the HTML rendering only |
| 26 | AC96 (1.7), p. 5 | s01_intro_3 l.15; s02_setup_8 l.20 | OK in arXiv v1 | journal data of the bib entry to be confirmed (its annote says so) |

Confirmed in bulk, no change: PPZ Sec. 0.2 (1), (5), Rem. 2.1 and (6), Defs 2.2, 2.5, 2.13, Props 2.4, 2.8, 2.9, 2.12
with (8), 2.14 (v3; (5) numbers the series, the symplectic condition is the display after it); EO Def. 4.2 (4-9),
kernel (4-3), Thms 4.1 (4-20), 4.2, 4.4, 4.6, Cor. 4.1, Def. 4.3 (4-10), (5-1), (5-4), (5-13) with footnote 8, Thm 5.1
(5-14), Sec. 5.3 (5-27)-(5-30), (A-25) (v4); DOSS Thm 4.1, Sec. 2.3 (2-22)-(2-26) (v1); DNOPS Sec. 2.6 (2-30)-(2-32),
(2-34), (2-35), the sentence before Thm 8.1, Remark 2.1 (v2); ABDKS Thm 7.1, (1.1), Sec. 1.1 (v3); GJR18 Sec. 4.1,
Thm 4.1, Lemma 5.2 (16)-(17), Remark 7.2 (v1); GJR17 Lemma 4.6 (19)-(20) (v1); Giv01 Sec. 2.2 Thm 1 and Remark (v4);
MR2276766 Thm 1, Cor. 1, Sec. 6 (10) (v2); LP18 Sec. 0.2 (5) (v4); LP20 Sec. 0.2 (5), (10) (v2); DLMF 24.2.3, 24.4.1
(version 1.2.8); Givental Cor. 11.7 does print l_alpha for lambda_alpha.
Number of pins to correct in the text: 18 (rows 1-18, of which rows 9, 10, 18 are ledger corrections); rows 19-26
need a check by the author but no change now; row 21 is with E5.

## 3. P1 items

### A4-01 (P1) J and the quantum product at t = T(q)H, Q = 1 are not defined by substitution [L6-m3; X5-F5, third bullet]
Checked: `s02_setup.tex` l.76-82 defines J(t,z) as a formal series in t and then puts "t = T(q)H" with
T = \logq + \tilde I_1/C_0; \logq is not topologically nilpotent, so the substitution is not literally defined. The
pin (l.87, step (d)) and `a1_cyclicLR.tex` l.13 use the divisor equation to give it a meaning, and (G0)(b) (l.100)
states the string and divisor equations for primary insertions only. Agree with both referees. The part of X5-F5
about thm:set-GT is with E4.

(a) `s02_setup.tex` l.80 (thm:set-mirror). Old:
```latex
where $H^\vee_a$ is the $\eta$-dual basis and the sum runs over the $(n,d)$ for which $\overline{\mathcal M}_{0,n+1}(\bbP^4,d)$ is defined. Then at $\lambda_i=\zeta^i$, with $\Qn=1$ and $t=T(q)H$,
```
New:
```latex
where $H^\vee_a$ is the $\eta$-dual basis and the sum runs over the $(n,d)$ for which $\overline{\mathcal M}_{0,n+1}(\bbP^4,d)$ is defined. For $t=\tau H$ the divisor equation (\cref{thm:set-G0}~(b)) gives $J(\tau H,z)=\re^{\tau H/z}\bigl(z+\sum_{d\geq1}(\Qn\re^{\tau})^d\sum_a\langle H^a/(z-\psi)\rangle_{0,1,d}H^\vee_a\bigr)$. Since $T=\logq+\tilde I_1/C_0$ contains the formal logarithm, $J(T(q)H,z)$ at $\Qn=1$ is \emph{defined} by this formula, with $\re^{dT}\coloneqq(\re^T)^d$, $\re^T=q\exp(\tilde I_1/C_0)$ (\cref{def:set-genus0}), and $\re^{TH/z}\coloneqq\sum_{k\geq0}(TH)^k/(k!\,z^k)$:
\[
J(TH,z)\coloneqq\re^{TH/z}\Bigl(z+\sum_{d\geq1}\re^{dT}\sum_{a\in\bbZ/5}\Bigl\langle\frac{H^a}{z-\psi}\Bigr\rangle_{0,1,d}H^\vee_a\Bigr).
\]
Then at $\lambda_i=\zeta^i$,
```
(b) `s02_setup.tex` l.100 ((G0)(b)). Old:
```latex
\item[(b)] The string and divisor equations hold: $E_{0,n+1,d}$ is the pull-back of $E_{0,n,d}$ along the map forgetting a point, and $\int_dH=d$.
```
New:
```latex
\item[(b)] The string and divisor equations hold, also with descendant insertions at the other points (with the usual correction terms, from $\psi_j=\pi^*\psi_j+D_j$): $E_{0,n+1,d}$ is the pull-back of $E_{0,n,d}$ along the map $\pi$ forgetting a point, and $\int_dH=d$.
```
(c) `s02_setup.tex` l.106 ((G0)(d)), append after "and $J_0=J/z$ (string equation).":
```latex
At $t=T(q)H$, $\Qn=1$ these series are understood as in \cref{thm:set-mirror}: by (b), $\langle\!\langle x,y,w\rangle\!\rangle(TH)\coloneqq\eta(xy,w)+\sum_{d\geq1}\re^{dT}\langle x,y,w\rangle_{0,3,d}$, and $J_a(TH,z)$ is $\partial J/\partial t^a$ at $t=\tau H$ with $(\Qn\re^{\tau})^d$ replaced by $\re^{dT}$ and $\re^{\tau H/z}$ by $\re^{TH/z}$.
```
(d) In the pin, l.87 step (d), "so by the divisor equation $J(TH,z)=\dots$" becomes "so the definition of $J(TH,z)$ reads".

### A4-02 (P1) thm:mod-var (ii) is proved from EO (5-15), declared "enters no proof" [X5-F3 = L2-m1; L2-m3]
Checked: `s06_model_8.tex` l.34 (statement (ii)), l.39 (proof: "(5-14) for $g+k>1$ and (5-15) for $g\geq 1$");
`s01_intro_3.tex` l.23 lists (5-15) among citations that enter no proof. (ii) is cited nowhere (grep: thm:mod-var is
used only in prop:mod-V1, part (i)). Agree. Of the three repairs offered I recommend keeping (ii) and proving it from
prop:mod-period, which precedes it and does not depend on it; then l.23 becomes true as it stands.
`s06_model_8.tex` l.39, replace the proof by:
```latex
\begin{proof}
(i) This is \cite[Thm.~5.1]{Eynard:2007kz}, (5-14) for $g+k>1$, for the temperature variation of \cite[\S5.3, (5-27)--(5-30)]{Eynard:2007kz}, translated by \cref{conv:mod-kernel}. There a variation at fixed $x$ is $\delta_\Omega y|_x\,\dd x=-\Omega$ \cite[(5-1)]{Eynard:2007kz} with $\Omega(p)=\int_{\partial\Omega}B(p,q)\Lambda(q)$ \cite[(5-13)]{Eynard:2007kz}, and (5-14) reads $\delta_\Omega W^{(h)}_n(J)|_{x}=\int_{\partial\Omega}W^{(h)}_{n+1}(J,q)\Lambda(q)$; on $\bbP^1$ there is no period matrix and no parameter $\kappa$, so the covariant variation $D_\Omega$ of \cite[Thm.~5.1]{Eynard:2007kz} is $\delta_\Omega$. By \cref{lem:mod-fixedx}, $\Omega=-\partial_v(y\,\dd x)|_x=-\int_{\mathfrak r}\omega_{0,2}(\cdot,q)$: the path $\partial\Omega=\mathfrak r$ joins the poles $0$, $\infty$ of $y\,\dd x$, with $\Lambda=-1$. With $\omega_{h,n}=(-1)^nW^{(h)}_n$ this gives $\partial_v|_x\omega_{h,n}=(-1)^{n}(-1)^{n+1}(-1)\int_{\mathfrak r}\omega_{h,n+1}=\int_{\mathfrak r}\omega_{h,n+1}$. The hypotheses hold: [keep the last sentence of the present proof unchanged].

(ii) This part is not used. By \cref{prop:mod-period}~(iii), $\partial_v\FTR_h=-\frac{B_{2h}}{2h}v^{1-2h}$, and by \cref{prop:mod-period}~(ii), $\int_{\mathfrak r}\omega_{h,1}=-\int_\infty^0\omega_{h,1}=-\frac{B_{2h}}{2h}v^{1-2h}$. (It is also the translation of \cite[(5-15)]{Eynard:2007kz}, which is not used.)
\end{proof}
```
`b_signs_2.tex` l.24: `\cite[Theorem~5.1, (5-14)--(5-15)]{Eynard:2007kz}` -> `\cite[Theorem~5.1, (5-14)]{Eynard:2007kz}`.
Alternative (shorter): delete item (ii) and its proof; nothing cites it.

### A4-03 (P1) rmk:set-qrrsign, last sentence: "forced, independently of (QRR)" is false for the actual rows [L5-F3]
Checked: `s02_setup_4.tex` l.52, last sentence: "That $\mathrm{ph}_1=-\frac3{20}$ is forced, independently of (QRR), is
part of \cref{lem:gau-HPF}." lem:gau-HPF (`s03_gauge_6.tex` l.65) gives $(4-5\theta_u)\mathrm{ph}_1=-\frac35$, whose
solutions in $\bbQ(\zeta)\fps q$ are $-\frac3{20}+\mathsf c\,u^{4/5}$ ($u^{4/5}=L^{-4}\in\bbQ\fps q$); by thm:gau-R0 (iv)
$\mathsf c=-[z^1]\Dn_0$, which vanishes by (QRR) (thm:con-glue), not by lem:gau-HPF. Agree. Remark only; no proof uses it.
Old: the sentence quoted above. New:
```latex
For every solution with coefficients in $\bbQ\fps u$, in particular for $\Phican$, the coefficient of $\zh$ is forced to be $-\frac3{20}$, independently of (QRR) (\cref{lem:gau-HPF}); so the twist sign of (QRR) is the one compatible with $\Dn_0=1$ at order $z$. For the actual rows this is not an independent derivation: there $\mathrm{ph}_1=-\frac3{20}-\mathsf d_1u^{4/5}$ with $\mathsf d_1=\coef{z}\Dn_0$ (\cref{thm:gau-R0}~(iv)), and $\mathsf d_1=0$ uses (QRR) (\cref{thm:con-glue}).
```

### A4-04 (P1) Missing ranges: three statements are false or undefined at the excluded index [L2-m2, L2-m4, L3-F6, L4-F3]
All four checked in the current text; agree with each.
(a) `s06_model.tex` l.21 (conv:mod-kernel). Old:
```latex
$\omega_{h,n}=(-1)^nW^{(h)}_n=(-1)^n\omega^{(h)}_{n,0}$ and $\FTR_h=-F^{(h)}$,
```
New:
```latex
$\omega_{h,n}=(-1)^n\omega^{(h)}_{n,0}$ for all $(h,n)$, $\omega_{h,n}=(-1)^nW^{(h)}_n$ for $2h-2+n\geq0$ (at $(h,n)=(0,1)$, $W^{(0)}_1=0$ \cite[(4-7)]{Eynard:2007kz}), and $\FTR_h=-F^{(h)}$,
```
and in `b_signs_2.tex` l.20 add "for $2h-2+n\geq0$" after `$\omega_{h,n}=(-1)^nW^{\mathrm{EO}}_{h,n}$`.
(b) `s06_model.tex` l.66 (thm:mod-tr). Old: `Let $v\neq 0$ and $2h-2+n>0$.` New: `Let $v\neq 0$, $n\geq1$ and $2h-2+n>0$.`
(for $n=0$ there is no $\omega_{h,0}$; the theorem is used with $n\geq1$ only).
(c) `s06_model_7.tex` l.5. Old: `with the Bergman kernel and $\omega^{(0)}_{1,0}=-y\,\dd x$ are`
New: `with the Bergman kernel and $\omega^{(0)}_{1,0}=-y\,\dd x$ are, for $(g,m)\neq(0,1)$,`
(at (0,1) the right side of (7.9) is 0; the note uses (7.9) for $2h-2+n>0$ and (7.10) for $h\geq1$).
(d) `b_signs_2.tex` l.23 (table row "dilaton equation"): add "for $2h-2+n>0$" after the displayed equation, and in the
source column replace "with the factor $2h+n-2$ for $W^{\mathrm{EO}}$" by "with the factor $2h+n-2$ for
$W^{\mathrm{EO}}$; the residues at the $p_i$ and the term $\delta_{g,0}\delta_{k,1}y\,\dd x$ printed there drop out in
this range".
