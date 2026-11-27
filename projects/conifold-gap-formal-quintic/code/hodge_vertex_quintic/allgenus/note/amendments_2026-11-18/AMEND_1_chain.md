# Suggested amendments 1: the proof chain (judgements E1-E8 and the README lists)

Part 1 of 4. Source: `audit_2026-11-18/README.md`, `audit_2026-11-18/verify/E1.md` ... `E8.md`, and the referee
reports cited there. Nothing in `main.tex` or `sections/` was changed; this file only proposes changes.

How to read an entry. `file lN` is the line of the current file (each paragraph of the note is one line). "old" quotes
the current text (a fragment long enough to be found; `...` marks an omission inside one line), "new" is the
replacement, in the note's macros. Every "old" fragment was checked against the current file. "Checked" says what was
re-derived for this list, and how. Labels proposed here that do not yet exist: `lem:set-Runique` (optional),
`eq:mod-edge` (in a new part (iii) of `prop:mod-Rtr`), `eq:con-factor`. No new macro and no new bibliography entry is needed,
except in the optional item D9.

Groups: A must-fix (A1-A3), B routine missing steps (B1-B5), C false as written and used by no proof (C1-C6),
D other wording (D1-D9), then the resulting statement of the inputs. Parts: this file A1; `_2` A2, A3; `_3` B; `_4` C, D
and the inputs.

Where this list departs from the judges (details at the entries):
* A1: E4 proposes to *redefine* R by flatness and its value at q = 0. Recommended here: keep the definition, declare
  flatness as part of (QRR), and (optionally) add the uniqueness statement as a lemma. Same logical content, smaller edit.
* B2: E1 repairs the proof of `eq:con-bz` in a ring with seven formal logarithm symbols. A shorter repair exists that
  uses no symbols, no branches and not `lem:con-form` (d): only `lem:fsp-stirfe` in the form St_m = Poch_m St_0.
* A2: on the disagreement E5 against U09/L3, see A2 (the two sides are talking about different conditions).

---------------------------------------------------------------------------------------------------------------------

## A. Must fix

### A1. Flatness of the R-matrix is an input, not a convention; all-genus axioms in (T)  [E4; L5-3, U01-4, U02-1, X5-6]

Confirmed in the text.
* `s02_setup_4.tex` l.9-15 (`def:set-rows`): "be the $R$-matrix of the semisimple cohomological field theory ... in
  Givental's convention: (a) ... (b) for every colour $i$ the column $R(z)\tilde e_i\re^{u_i/z}$ solves the quantum
  differential equation ... (c) ..." and l.15 names only the value at q = 0 and the graph sum as quoted inputs.
* `s01_intro_3.tex` l.10: (QRR) is "The value of the $R$-matrix at the large-radius point ... in Givental's convention
  (the columns ... solve the quantum differential equation, and the legs ... carry $R^{-1}(\psi)$)".
* (b) is used in the proof at `s02_setup_4.tex` l.24 ("\cref{def:set-rows}~(b) and the argument of ..."), which gives
  `prop:set-cyclicLR` (iii) and the last sentence of `cor:set-cyclicPF5` (b); these feed `lem:gau-colrows` Steps 1 and 4,
  `thm:gau-R0` (iii), `cor:gau-rowsNN`, `thm:con-glue` and `lem:mai-Psi`. Both statements open with "Assume (M) and (G0)."
* `s02_setup_5.tex` l.42 (`thm:set-GT`): "(which, by the divisor equation, is the theory at the point $t=T(q)H$ with
  $\Qn=1$)"; (G0) (b) gives the divisor equation in genus 0 only.

Checked. (1) The pins: in the stored text of arXiv:1812.11908v1 (`audit_2026-11-18/work/L5/html/`), Lemma 5.2 has the
equation (16) `(z q d/dq + L_alpha) R_{i alpha-bar} = sum_j A_i^j R_{j alpha-bar}` and (17) follows in its proof; in the
stored text of arXiv:0712.0160v3 the statements are numbered 7.7, 8.2, 8.4, 8.7 as E4 says. The translation of (16) into
`prop:set-cyclicLR` (iii) (E4 3.1) was re-derived: R_{i alpha-bar} = eta(H^i, R e-bar_alpha) is the row Q^{(alpha)}_i, and
A H^i = C_1 kappa^LR_i H^{i+1}. (2) E4's claim that (b) and the value at q = 0 determine R, and that (a) follows, was
re-derived by hand (proof written out in A1.16). I did not re-read Givental, Coates-Givental or GJR17.

Disagreement with E4, and recommendation. E4 (5)(1) replaces the definition of R by "the unique series with (b) and
the value at q = 0". That is correct, but it changes the object named in Theorems A and C and every sentence that says
"the R-matrix of the formal quintic". The same gain is obtained by (i) declaring (b) as part of (QRR) and (ii) adding the
uniqueness as a lemma (A1.16, optional). Then (QRR) = "the R-matrix of (T) is flat and has the value diag(h) at q = 0",
and symplecticity need not be quoted. The edits A1.1-A1.15 are needed in either case.

**A1.1** `s01_intro_3.tex` l.10, the whole item (QRR). Old: `\item[(QRR)] \emph{The value of the $R$-matrix at the
large-radius point}: $R(q=0)=\diag(h_{\bw_0},\dots,h_{\bw_4})$ in the normalised canonical basis at $q=0$, in Givental's
convention (the columns ... solve the quantum differential equation, and the legs of the graph sum carry $R^{-1}(\psi)$)
(\cref{thm:set-qrr}). It is used for every colour, ... Sources: ... \cite[Theorem~1, Corollary~1]{MR2276766}.` New:

```latex
\item[(QRR)] \emph{The $R$-matrix along the line: flatness and the value at the large-radius point} (\cref{thm:set-qrr}). The $R$-matrix $R$ of (T) has coefficients in $\End(\bbH)\otimes\bbQ(\zeta)\fps q$ and satisfies: (\emph{flatness}) for every colour $i$ the column $R\,\tilde e_i\,\re^{u_i/z}$ solves the quantum differential equation $z\partial_Tv=H*v$ along the line (\cref{def:set-rows}~(b)); (\emph{value}) $R(q=0)=\diag(h_{\bw_0},\dots,h_{\bw_4})$ in the normalised canonical basis at $q=0$. Flatness is a theorem about the $R$-matrix of the cohomological field theory, not a convention; the convention (Givental's) is only that the legs of the graph sum carry $R^{-1}(\psi)$. Flatness is used in \cref{prop:set-cyclicLR}~(iii) and in the last sentence of \cref{cor:set-cyclicPF5}~(b), on which \cref{lem:gau-colrows}, \cref{thm:gau-R0}~(iii), \cref{cor:gau-rowsNN}, \cref{thm:con-glue} and \cref{lem:mai-Psi} rest. The value is used for every colour, in the form $\Phi_i(z,1)=h_{\bw_i}(z)$ (\cref{eq:set-qrr}): in \cref{lem:gau-colrows} and, for colour $0$, in the glue \cref{thm:con-glue}. Sources for flatness: \cite[Lemma~5.2, (16)]{GJR18}, which is \cref{prop:set-cyclicLR}~(iii) at $\lambda=1$; \cite[\S4.2, (12)--(13)]{GJR17}; in general \cite[Section~1, Proposition (a)--(c)]{Giv01} and \cite[Props.~7.7, 8.2~(ii), 8.4, 8.7~(ii)]{MR2917177}. Sources for the value: \cite[(17), in the proof of Lemma~5.2]{GJR18}; \cite[(19)--(20), in the proof of Lemma~4.6, at $t=0$]{GJR17}; \cite[\S2.2, Theorem~1 and the Remark after it]{Giv01}; \cite[Theorem~1, Corollary~1]{MR2276766}.
```

(The clause "in \cref{def:set-rows}" of the old list of uses is dropped: `def:set-rows` uses R(0) = 1, not the value.)

**A1.2** `s01_intro_3.tex` l.13, the whole item (T). Old: `\item[(T)] \emph{Teleman's reconstruction}
\cite[Theorem~(2) in (1.15)]{MR2917177} with \cite[\S4.1, Theorem~4.1]{GJR18}: the cohomological field theory of the
formal quintic is the graph sum of \cref{def:set-graphsum} with the $R$-matrix of \cref{thm:set-qrr} (\cref{thm:set-GT}).
In \cref{thm:int-main,thm:int-D} it identifies ...; in \cref{cor:int-gap} it identifies the total with $F_g$.` New:

```latex
\item[(T)] \emph{Teleman's reconstruction, with the axioms of the twisted theory in all genera} (\cref{thm:set-GT}). (i) At $\lambda_i=\zeta^i$ the classes $\Omega_{g,n}$ of \cref{def:set-fq}, in all genera, form a cohomological field theory with flat unit $1$ for $\eta$ \cite[\S4.1]{GJR18}, \cite[Def.~0.3]{PPZ15}, and satisfy the divisor equation: for $2g-2+n>0$ the part of degree $d$ of $\pi_*\Omega_{g,n+1}(\gamma_1,\dots,\gamma_n,H)$ is $d$ times that of $\Omega_{g,n}(\gamma_1,\dots,\gamma_n)$, where $\pi$ forgets the last point; so the theory at $\Qn=\re^{T(q)}$ is the theory at the point $t=T(q)H$ with $\Qn=1$. (ii) This theory, which is semisimple (\cref{cor:set-cyclicPF5}~(a)), is the graph sum of \cref{def:set-graphsum}, $\Omega=R.\TFT$, for an $R$-matrix $R$ with coefficients in $\End(\bbH)\otimes\bbQ(\zeta)\fps q$, symplectic and with $R(0)=1$ \cite[Theorem~(2) in (1.15)]{MR2917177}, \cite[\S4.1, Theorem~4.1]{GJR18}; this is the $R$-matrix of \cref{def:set-rows}, to which (QRR) refers. In \cref{thm:int-main,thm:int-D} (T) identifies that $R$-matrix, whose rows define $T_\Gamma$ and $\Dn$; in \cref{cor:int-gap} it identifies the total with $F_g$.
```

**A1.3** `s01_intro_3.tex` l.17 ("Where the inputs enter"). Two replacements.
Old: `(QRR) at the three places named above;` New: `(QRR) at the places named above: flatness in \cref{prop:set-cyclicLR}~(iii) and \cref{cor:set-cyclicPF5}~(b), the value in \cref{lem:gau-colrows} and \cref{thm:con-glue};`
Old: `\Cref{sec:gau} uses no quoted input beyond the genus-zero data of \cref{sec:set};` New: `\Cref{sec:gau} uses no quoted input beyond (M), (G0) and (QRR), which enter through the rows of \cref{sec:set} (flatness) and in \cref{lem:gau-colrows} (the value);`
(The old sentence was already inexact: `lem:gau-colrows` Step 3, in `sec:gau`, says that (QRR) enters there.)

**A1.4** `s01_intro_3.tex` l.19 and, identically, `s10_main.tex` l.33 ("Conditionality"). Old: `the value (QRR) of the
$R$-matrix at the large-radius point, the genus-$0$ mirror theorem (M) and the genus-$0$ axioms (G0) of the twisted
theory, Teleman's reconstruction (T),` New: `the flatness of the $R$-matrix and its value at the large-radius point
(QRR), the genus-$0$ mirror theorem (M) and the genus-$0$ axioms (G0) of the twisted theory, Teleman's reconstruction
with the axioms of the twisted theory in all genera (T),`

**A1.5** `s01_intro_3.tex` l.21. Old: `\Cref{thm:int-D} uses the $R$-matrix that (T) identifies only through
\cref{def:set-rows}~(b), $R(0)=1$ and (QRR), and the genus-zero data` New: `\Cref{thm:int-D} uses the $R$-matrix that
(T) identifies only through $R(0)=1$, its coefficient ring $\bbQ(\zeta)\fps q$ and (QRR), that is, flatness
(\cref{def:set-rows}~(b)) and the value at $q=0$, and the genus-zero data`

**A1.6** `s02_setup_4.tex` l.9 and l.15 (`def:set-rows`).
Old (l.9): `be the $R$-matrix of the semisimple cohomological field theory of the formal quintic along the line $t=T(q)H$ \cite[Section~4.1]{GJR18}, in Givental's convention:`
New: `be the $R$-matrix of the semisimple cohomological field theory of the formal quintic along the line $t=T(q)H$ \cite[Section~4.1]{GJR18}, that is, the $R$-matrix of the quoted input (T), \cref{thm:set-GT}. It has the following properties: (a) is part of the quoted input (T), (b) is part of the quoted input (QRR), \cref{thm:set-qrr}, and (c) is Givental's convention:`
Old (l.15): `Its value at $q=0$ is the quoted input (QRR), \cref{thm:set-qrr}, and its graph sum is the cohomological field theory of the formal quintic, the quoted input (T), \cref{thm:set-GT}.`
New: `The flatness (b) and the value of $R$ at $q=0$ form the quoted input (QRR), \cref{thm:set-qrr}; that the graph sum of $R$ is the cohomological field theory of the formal quintic is the quoted input (T), \cref{thm:set-GT}, which includes (a). Property (b) is not a convention: it says that the $R$-matrix of the cohomological field theory is the asymptotic solution of the quantum differential equation, and it is used in \cref{prop:set-cyclicLR}~(iii).`
(With A1.16 add: `By \cref{lem:set-Runique}, (b) and the value at $q=0$ determine $R$, and (a) follows from them.`)

**A1.7** `s02_setup_4.tex` l.29-39 (`thm:set-qrr` and the sentence after it).
Old (l.29): `\begin{thm}[The $R$-matrix at $q=0$, (QRR); quoted]` New: `\begin{thm}[Flatness of the $R$-matrix and its value at $q=0$, (QRR); quoted]`
Old (l.36): `In Givental's convention (\cref{def:set-rows}), $R(q=0)=\diag(h_{\bw_0},\dots,h_{\bw_4})$ in the basis $\tilde e_i(q=0)$:`
New: `The $R$-matrix of \cref{def:set-rows} has coefficients in $\End(\bbH)\otimes\bbQ(\zeta)\fps q$ and satisfies \cref{def:set-rows}~(b); and $R(q=0)=\diag(h_{\bw_0},\dots,h_{\bw_4})$ in the basis $\tilde e_i(q=0)$:`
Old (l.39): `This is \cite[proof of Lemma~5.2, (17)]{GJR18} and \cite[proof of Lemma~4.6, (19)--(20)]{GJR17} at $t=0$; both refer to ... \cite[Theorem~1, Corollary~1]{MR2276766}. It is used only in the following form, for every colour $i$:`
New: `Flatness is \cite[Lemma~5.2, (16)]{GJR18}: there $R_{i\bar\alpha}=\eta(H^i,R\,\bar e_\alpha)$ is our row $\Qr^{(\alpha)}_i$, $L_\alpha=\zeta^\alpha L$ at $\lambda=1$, and the matrix on the right is that of $C_1H*$, so that (16) is literally \cref{prop:set-cyclicLR}~(iii); see also \cite[\S4.2, (12)--(13)]{GJR17}, and in general \cite[Section~1, Proposition (a)--(c)]{Giv01} and \cite[Props.~7.7, 8.2~(ii), 8.4, 8.7~(ii)]{MR2917177}. The value at $q=0$ is \cite[proof of Lemma~5.2, (17)]{GJR18} and \cite[proof of Lemma~4.6, (19)--(20)]{GJR17} at $t=0$; both refer to \cite[Section~2.2, Theorem~1 and the Remark after it]{Giv01} and \cite[Theorem~1, Corollary~1]{MR2276766}. The value is used only in the following form, for every colour $i$:`

**A1.8** `s02_setup_2.tex` l.7 (`prop:set-cyclicLR`). Old: `Assume (M) and (G0). Along the line` New: `Assume (M) and (G0), and for (iii) the flatness part of (QRR), \cref{def:set-rows}~(b). Along the line`

**A1.9** `s02_setup_3.tex` l.5 (`cor:set-cyclicPF5`). Old: `Assume (M) and (G0).` New: `Assume (M) and (G0), and for the last sentence of (b) the flatness part of (QRR), \cref{def:set-rows}~(b).`

**A1.10** `s02_setup_5.tex` l.42 (`thm:set-GT`). Old: `The cohomological field theory of the formal quintic along the
line, i.e.\ the classes $\Omega_{g,n}$ of \cref{def:set-fq} with $\Qn=\re^{T(q)}$ (which, by the divisor equation, is the
theory at the point $t=T(q)H$ with $\Qn=1$), equals $R.\TFT$, with $R$ of \cref{def:set-rows} and $\TFT$ of
\cref{def:set-frob} \cite[Theorem~(2) in (1.15)]{MR2917177}, \cite[Section~4.1, Theorem~4.1]{GJR18}.` New:

```latex
At $\lambda_i=\zeta^i$ the classes $\Omega_{g,n}$ of \cref{def:set-fq}, in all genera, form a cohomological field theory with flat unit $1$ for $\eta$ over $\bbQ(\zeta)\fps{\Qn}$ \cite[Section~4.1]{GJR18}, \cite[Def.~0.3]{PPZ15}, and satisfy the divisor equation in all genera: for $2g-2+n>0$ the part of degree $d$ of $\pi_*\Omega_{g,n+1}(\gamma_1,\dots,\gamma_n,H)$ is $d$ times that of $\Omega_{g,n}(\gamma_1,\dots,\gamma_n)$, $\pi$ forgetting the last point (as in genus $0$, \cref{thm:set-G0}~(b): $E_{g,n+1,d}$ is the pull-back of $E_{g,n,d}$ and $\int_dH=d$). Hence the classes $\Omega_{g,n}$ with $\Qn=\re^{T(q)}$ are the theory at the point $t=T(q)H$ with $\Qn=1$, the cohomological field theory of the formal quintic along the line. It equals $R.\TFT$, with $\TFT$ of \cref{def:set-frob} and an $R$-matrix $R=1+O(z)$ with coefficients in $\End(\bbH)\otimes\bbQ(\zeta)\fps q$, symplectic as in \cref{def:set-rows}~(a) \cite[Theorem~(2) in (1.15)]{MR2917177}, \cite[Section~4.1, Theorem~4.1]{GJR18}; this $R$ is the $R$-matrix of \cref{def:set-rows}.
```

(The rest of the theorem, "In particular, for $g\geq2$, ...", is unchanged.)

**A1.11** `s08_connection_8.tex` l.19, last sentence of the proof of `thm:con-glue`. Old: `The inputs used here are
(QRR) and the actual rows with their WKB form, which rest on (M) and (G0) (\cref{thm:set-mirror,thm:set-G0}) through
\cref{cor:set-cyclicPF5} and \cref{thm:gau-R0};` New: `The inputs used here are the value at $q=0$ of (QRR) and the
actual rows with their WKB form, which rest on (M), (G0) (\cref{thm:set-mirror,thm:set-G0}) and the flatness part of
(QRR) (\cref{def:set-rows}~(b)) through \cref{cor:set-cyclicPF5} and \cref{thm:gau-R0};`

**A1.12** `s09_gamma_3.tex` l.82. Old: `\cref{thm:con-glue} also uses the rows and their WKB form, which rest on (M)
and (G0).` New: `\cref{thm:con-glue} also uses the rows and their WKB form, which rest on (M), (G0) and the flatness part
of (QRR) (\cref{def:set-rows}~(b)).`

**A1.13** `s10_main.tex` l.31 ("Where the inputs enter"). Old: `(QRR) enters through the rows (\cref{thm:set-qrr}),
\cref{lem:gau-colrows} and the glue \cref{thm:con-glue}; (M) and (G0) through
\cref{thm:set-mirror,prop:set-cyclicLR,cor:set-cyclicPF5};` New: `(QRR) enters through the rows: its flatness part
(\cref{def:set-rows}~(b)) in \cref{prop:set-cyclicLR}~(iii) and \cref{cor:set-cyclicPF5}~(b), its value
(\cref{thm:set-qrr}) in \cref{lem:gau-colrows} and the glue \cref{thm:con-glue}; (M) and (G0) through
\cref{thm:set-mirror,prop:set-cyclicLR,cor:set-cyclicPF5};`
(`s10_main.tex` l.33 is A1.4.)

**A1.14** `main.tex` l.18 (abstract, list of inputs). Old: `the value of the $R$-matrix at $q=0$ given by quantum
Riemann--Roch; the genus-$0$ mirror theorem and the genus-$0$ axioms of the twisted theory; Teleman's reconstruction;`
New: `the fact that the $R$-matrix is the asymptotic solution of the quantum differential equation, and its value at
$q=0$ given by quantum Riemann--Roch; the genus-$0$ mirror theorem and the genus-$0$ axioms of the twisted theory;
Teleman's reconstruction, with the axioms of the twisted theory in all genera;`

**A1.15** `b_signs.tex` l.12 (appendix B, "Givental's convention"). Old: `$R(z)=1+O(z)$ is taken in Givental's
convention, $S\sim R\,\re^{u/z}$: the columns $R(z)\tilde e_i\,\re^{u_i/z}$ solve $z\de_Tv=H*v$ along $t=T(q)H$, and
$R(-z)^T\eta R(z)=\eta$ in the flat basis (\cref{def:set-rows});` New: `$R(z)=1+O(z)$ is taken in Givental's convention,
$S\sim R\,\re^{u/z}$, with legs $R^{-1}(\psi)$. That the columns $R(z)\tilde e_i\,\re^{u_i/z}$ solve $z\de_Tv=H*v$ along
$t=T(q)H$ is not part of the convention but the quoted property (QRR) (\cref{def:set-rows}~(b), \cref{thm:set-qrr}); and
$R(-z)^T\eta R(z)=\eta$ in the flat basis (\cref{def:set-rows}~(a));`

**A1.16 (optional, recommended)** New lemma after the proof at `s02_setup_4.tex` l.23-25 (before l.27). It makes the
input (QRR) minimal (symplecticity is then proved) and shows that "the R-matrix with (b) and the value diag(h)" is one
definite object, which is what E4 (5)(1) wants.

```latex
\begin{lem}[Flatness and the value at $q=0$ determine $R$]
\label{lem:set-Runique}
Assume (M) and (G0). There is exactly one series $R(z)=1+\sum_{k\geq1}R_kz^k$, $R_k\in\End(\bbH)\otimes\bbQ(\zeta)\fps q$, which satisfies \cref{def:set-rows}~(b) and $R(0,z)\tilde e_i(0)=h_{\bw_i}(z)\tilde e_i(0)$ for every $i$; and it satisfies \cref{def:set-rows}~(a).
\end{lem}

\begin{proof}
Write $R(z)\tilde e_i=\sum_jY_{ji}(z)\tilde e_j$, $Y=1+\sum_{k\geq1}Y_kz^k$ with entries in $\bbQ(\zeta)\fps q$, and $\Dq\tilde e_i=\sum_j\Theta_{ji}\tilde e_j$ ($\Theta$ a local letter). As $\eta(\tilde e_i,\tilde e_j)=\delta_{ij}$ is constant, $\Theta$ is antisymmetric, and its entries lie in $q\bbQ(\zeta)\fps q$, because $\Dq$ maps $\bbQ(\zeta)\fps q$ into $q\bbQ(\zeta)\fps q$. Since $C_1H*\tilde e_j=\zeta^jL\,\tilde e_j$ (\cref{cor:set-cyclicPF5}~(a)), property (b), multiplied by $C_1$, reads $z\Dq Y+z\Theta Y=[\Dq U,Y]$ with $\Dq U\coloneqq\diag(\Dq u_j)=\diag(\zeta^jL)$, that is,
\[
(\Dq Y_k+\Theta Y_k)_{ji}=(\zeta^j-\zeta^i)L\,(Y_{k+1})_{ji}\qquad(k\geq0,\ Y_0=1).
\]
For $j\neq i$ this determines the off-diagonal entries of $Y_{k+1}$ from $Y_k$, as $(\zeta^j-\zeta^i)L$ is a unit; they lie in $q\bbQ(\zeta)\fps q$, because the left side does. For $j=i$, with $k+1$ in place of $k$, it reads $\Dq(Y_{k+1})_{ii}=-\sum_{l\neq i}\Theta_{il}(Y_{k+1})_{li}\in q\bbQ(\zeta)\fps q$, which determines $(Y_{k+1})_{ii}$ up to its constant term; the equation with $j=i$ and $k=0$ holds, as $\Theta_{ii}=0$. So the solutions $Y$ correspond to the families of constant terms of the diagonal entries, $Y(q=0)$ is diagonal, and the condition at $q=0$ fixes it as $\diag(h_{\bw_i})$. For (a) put $N(z)\coloneqq Y(-z)^TY(z)$ (a local letter). From the equation and its transpose at $-z$, using $\Theta^T=-\Theta$, $z\Dq N=[\Dq U,N]$. For $j\neq i$ the lowest coefficient in $z$ of $N_{ji}$ is annihilated by the unit $(\zeta^j-\zeta^i)L$, so $N_{ji}=0$; and $\Dq N_{ii}=0$, so $N_{ii}$ is its value $h_{\bw_i}(-z)h_{\bw_i}(z)=1$ at $q=0$, the exponent of $h_{\bw_i}$ being odd in $z$ \eqref{eq:set-hw}. Hence $Y(-z)^TY(z)=1$, which is (a) in the orthonormal basis $\tilde e_i$.
\end{proof}
```

Checked: the recursion and the identity z D N = [D U, N] were derived by hand; `\eqref{eq:set-hw}` is stated after
this place in the file (in `thm:set-qrr`), so either move the lemma after `lem:set-hrot` (l.60) or keep the forward
reference. If the lemma is adopted, the words "symplectic and" may be dropped from the quoted part of A1.2 and A1.10.

Not changed, and why: `s01_intro.tex` l.16 and `s01_intro_2.tex` l.36, l.46 name (QRR) for the value only, which stays
correct; (G0) is unchanged (the all-genus statements go into (T), as E4 proposes; putting them into (G0) would make
the name "genus-zero axioms" wrong).
