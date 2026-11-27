# AMEND_4: cross-cutting and literature reports (part 2 of 3): P1 continued, P2

Conventions as in `AMEND_4_cross_lit.md`. "Old" texts were read in the current `sections/*.tex`.

## 3. P1 items (continued)

### A4-05 (P1) Two sentences on the logical order and on the inputs are false as written [X1-F1; X1-F2 = X5-F2]
(a) `s02_setup.tex` l.14. Checked: "Four facts used in the definitions are proved later, independently of this
section: \cref{lem:gau-confexp}, ..." Each of the four uses Section 2 (lem:gau-confexp uses def:set-confexp;
cor:gau-rowsNN uses thm:gau-R0, hence cor:set-cyclicPF5; lem:gau-colrows uses def:set-rows, thm:set-qrr; thm:gam-D
uses all of them). What X1 verified is that none uses the definitions that cite it. Agree.
Old: `are proved later, independently of this section:` New: `are proved later, without using the definitions that cite them:`
(b) `s01_intro_3.tex` l.17. Checked: "\Cref{sec:gau} uses no quoted input beyond the genus-zero data of \cref{sec:set};"
while `s03_gauge_7.tex` l.62 says "This step is the only place in this section where (QRR) enters", and the same
paragraph of the introduction lists lem:gau-colrows among the places of (QRR). Agree. Old: the clause quoted. New:
```latex
\Cref{sec:gau} uses, beyond the genus-zero data of \cref{sec:set}, only (QRR), in \cref{lem:gau-colrows} (hence in \cref{cor:gau-rowsNN}~(iii) and \cref{cor:gau-colours}), and the rows of \cref{def:set-rows}~(b), in \cref{thm:gau-R0}~(iii)--(iv);
```
E4 changes the status of def:set-rows (b) from "convention" to a quoted property and also touches this line; apply the
two changes together. In the same line, old `(M) and (G0) enter only in \cref{sec:set} (` new
`(M) and (G0) enter only in \cref{sec:set} and its appendix \cref{sec:set-lrproof} (` (X5, Section 2: the proofs of
prop:set-cyclicLR (i)-(ii) are in `a1_cyclicLR.tex`).

### A4-06 (P1) def:set-fq: the reason why the localisation sum can be specialised at lambda_i = zeta^i is partial [L6-m4]
Checked: `s02_setup.tex` l.24: "(The second term, $\bbE^\vee\otimes\cO(5)|_{p_i}$ on the contracted components, has an
invertible top Chern class there, as $5\zeta^i\neq0$.)" Only one family of denominators is treated. The others are
non-zero for the reason L6 gives (a fifth root of unity is not on the segment joining two others; I checked the list
of weights: for $m=i$ the weight is $b(\zeta^j-\zeta^i)/e$, zero only for the removed zero weight). Agree. Old: the
parenthesis quoted. New:
```latex
(No denominator of the localisation formula vanishes there. The second term, $\bbE^\vee\otimes\cO(5)|_{p_i}$ on the contracted components, has an invertible top Chern class, as $5\zeta^i\neq0$. The other denominators are products of non-zero rational multiples of $\zeta^i-\zeta^j$ ($j\neq i$), of classes $(\zeta^i-\zeta^j)/e-\psi$, and of the non-zero weights $(a\zeta^i+b\zeta^j)/e-\zeta^m$, $a+b=e$, $a,b\geq0$, of $H^0(f^*T_{\bbP^4})$ on an edge of degree $e$ over the line $p_ip_j$; these do not vanish at $\lambda=\zeta$, because for $m\notin\{i,j\}$ the point $\zeta^m$ is not on the segment $[\zeta^i,\zeta^j]$, and for $m\in\{i,j\}$ the weight is a non-zero multiple of $\zeta^i-\zeta^j$.)
```

### A4-07 (P1) "pairing" must be a perfect pairing over a ring [L1-M2]
Checked: `s02_setup_4.tex` l.70 "on a free module $V$ of finite rank with pairing $\eta$ over a commutative
$\bbQ$-algebra", the same words in `s02_setup_5.tex` l.5 (lem:set-pushins) and `s07_comparison.tex` l.19
(def:flh-outer: "of finite rank with pairing $\eta$"); all three use $\eta^{-1}$. Agree. In the three places, old
`with pairing $\eta$` new `with a perfect symmetric pairing $\eta$ (its Gram matrix is invertible)`.
The repair text of E6 for the last sentence of def:set-graphsum already says "eta perfect"; this item covers the
openings of the three statements. In all uses $\eta$ is $5$ times a permutation matrix, the identity, or $-\frac15$.

### A4-08 (P1) thm:mod-MLK: the passage from complex v to the formal variable v is not said [X2-m5]
Checked: `s06_model_10.tex` l.12-24. ML(h) (def:red-ML) is an identity over $\bbQ[v]$ with $v$ formal; thm:mod-tr,
thm:mod-Lreg and thm:mod-V are stated for complex $v\neq0$; the proof does not connect them. Agree; one sentence closes
it. Insert at the beginning of Step 2 (l.17), before "By \cref{thm:mod-tr}, multilinearity ...":
```latex
\Cref{thm:mod-tr,thm:mod-Lreg,thm:mod-V} are stated for a complex number $v\neq0$, and $\ML{h}$ for the variable $v$. This makes no difference: the correlator $\langle X_1\cdots X_n\rangle_h$ is a finite $\bbQ[v]$-combination of monomials $Cv^e$, $C\in\bbQ$, $e\in\bbZ$ (\cref{lem:set-modhom}), every term on the right below is the value at $v$ of an element of $\bbQ[v,v^{-1}]$, and an identity between elements of $\bbQ[v,v^{-1}]$ that holds for all complex $v\neq0$ holds identically.
```
(The derivatives $(\Fm_h)^{(n)}$ that occur, of order $\geq1$ for $h=1$ and $\geq3$ for $h=0$, are in $\bbQ[v,v^{-1}]$.)

### A4-09 (P1) Three assertions rest on computations that are not in the evidence appendix [X5-F9; README list 4]
Checked: `s01_intro_2.tex` l.25 "it is not implied by a statement about the singular part of the sum
(\cref{rmk:mai-graphlocal})"; `s10_main_2.tex` l.46 (rmk:mai-graphlocal) "computations (\cref{sec:evi}) show that six
edge graphs of genus $4$ become singular ..." with a draft note saying the computation is not in the tables; l.51
(rmk:mai-necessity) likewise. None is a proof step. Agree. Either record the computations in `c_evidence*.tex`
(then no text change), or:
(a) `s01_intro_2.tex` l.25. Old: `The statement is about each graph separately; it is not implied by a statement about the singular part of the sum (\cref{rmk:mai-graphlocal}).`
New: `The statement is about each graph separately; a statement about the singular part of the sum does not address the individual terms (see \cref{rmk:mai-graphlocal}).`
(b) `s10_main_2.tex` l.46. Old: `computations (\cref{sec:evi}) show that six edge graphs of genus $4$ become singular while the singular part of the total is unchanged. So the singular part of the total does not determine those of the graph terms, and no argument that uses the total only through its singular part can prove Theorem~A. The proof is graph-local,`
New: `a computation that is not recorded in \cref{sec:evi} indicates that six edge graphs of genus $4$ become singular while the singular part of the total is unchanged. If so, the singular part of the total does not determine those of the graph terms. In any case the proof of Theorem~A is graph-local,`
(c) `s10_main_2.tex` l.51. Old: `Computations (\cref{sec:evi}) give dressings` New: `Computations that are not recorded in \cref{sec:evi} indicate that there are dressings`;
old `: neither hypothesis of \cref{thm:red-main} can be dropped in general.` new `; on this evidence neither hypothesis of \cref{thm:red-main} can be dropped in general.`
The two draft notes then go.

### A4-10 (P1) sec:sgn: "only the $B_{2k}$ enter" is false for Theorem Gamma (b) [X3-m2 = L6-m5]
Checked: `b_signs.tex` l.34, last sentence: "only the $B_{2k}$, $k\geq1$, enter, so the conventions $\xi/(\re^\xi-1)$ and
$\xi/(1-\re^{-\xi})$ give the same formulas." The polynomials $B_{5l+1}(\frac15)$ of $\Ginf$ and $B_n(\theta)$ of
lem:fsp-dlmf5118 depend on the convention. The mathematics is right (def:fsp-stirling fixes the convention). Agree.
Old: the sentence quoted (from "only the"). New:
```latex
in $h_{\bw}$, $\Gst$ and Stirling's series only the numbers $B_{2k}$, $k\geq1$, enter, so there the conventions $\xi/(\re^\xi-1)$ and $\xi/(1-\re^{-\xi})$ give the same formulas ($B_1=-\frac12$ occurs only inside the proof of \cref{lem:fsp-stirfe}). The polynomials $B_n(\theta)$ of \cref{lem:fsp-dlmf5118} and the values $B_{5l+1}(\frac15)$ in $\Ginf$ (\cref{thm:gam-ginf}) depend on the convention, which is that of \cite[24.2.3]{DLMF}; the other one would replace $B_n(\theta)$ by $B_n(\theta+1)=B_n(\theta)+n\theta^{n-1}$.
```

## 4. P2 items

### A4-11 (P2) GJR18 "Section 3.1, (9)" [L5-F4 = L6-m2]
Checked: `s01_intro.tex` l.10 "(see \cite[\S3.1, (9)]{GJR18})"; `s02_setup.tex` l.24 "see \cite[Section~3.1, (9)]{GJR18}".
Two referees located (9) in Sec. 3.2 of arXiv v1 (p. 14). Agree.
`s01_intro.tex` l.10: `\cite[\S3.1, (9)]{GJR18}` -> `\cite[\S3.1; \S3.2, (9)]{GJR18}`.
`s02_setup.tex` l.24: `\cite[Section~3.1, (9)]{GJR18}` -> `\cite[Section~3.1; Section~3.2, (9)]{GJR18}`.

### A4-12 (P2) The specialisation is not in ZZ08 [L6-m1]
Checked: `s01_intro.tex` l.10 "at the specialisation $\lambda_i=\zeta^i$, $\zeta=\re^{2\pi\ri/5}$, of \cite{ZZ08}";
`s02_setup.tex` l.24 "This is the specialisation of \cite{ZZ08};" `s01_intro_3.tex` l.15 "The definitions use ...,
\cite{ZZ08} and \cite{MR1666787}." L6 searched the full text of arXiv:0710.0889v1: no torus weights. Agree (single
referee; the attribution comes from GJR18 p. 3).
`s01_intro.tex` l.10: old `of \cite{ZZ08} (see \cite[\S3.1, (9)]{GJR18})` new
`(see \cite[\S3.1; \S3.2, (9)]{GJR18}; the resulting hypergeometric series is studied in \cite{ZZ08})`.
`s02_setup.tex` l.24: old `This is the specialisation of \cite{ZZ08}; see \cite[Section~3.1, (9)]{GJR18}.` new
`This is the specialisation of \cite[Section~0.2]{LP18} and \cite[Section~0.2]{LP20}; see \cite[Section~3.1; Section~3.2, (9)]{GJR18}. The resulting hypergeometric series (\cref{lem:set-periods}) is studied in \cite{ZZ08}.`
`s01_intro_3.tex` l.15: old `\cite[\S0.2, (5), (10)]{LP20}, \cite{ZZ08} and \cite{MR1666787}.` new
`\cite[\S0.2, (5), (10)]{LP20} and \cite{MR1666787}; \cite{ZZ08} is cited for the series of \cref{lem:set-periods} only.`

### A4-13 (P2) (G0)(d): GJR18 (14) is not the quantum differential equation quoted [L5-F5]
Checked: `s02_setup.tex` l.102 "The quantum differential equation \cite[Section~5.2, (14)]{GJR18}"; `s01_intro_3.tex`
l.12 "the quantum differential equation \cite[\S5.2, (14)]{GJR18}". (14) is $z\Dq S+S\cdot H=A\cdot S$; the general
equation is the display before it. Agree (single referee). In both places:
`(14)]{GJR18}` -> `the display before (14), of which (14) is the form in the $q$-direction]{GJR18}`.
Optionally add in l.106: "(d) follows from (b), (c) and the genus-zero topological recursion relation."

### A4-14 (P2) Introduction: Giv01 and Coates-Givental do not state (QRR) for this theory [L5-F8]
Checked: `s01_intro.tex` l.12 "at $q=0$ quantum Riemann--Roch \cite[\S2.2, Theorem~1]{Giv01}, \cite[Corollary~1]{MR2276766},
\cite[(17)]{GJR18} gives". thm:set-qrr's own sentence (`s02_setup_4.tex` l.39) is accurate. Agree. Old: the clause
quoted. New:
```latex
at $q=0$ quantum Riemann--Roch, in the form \cite[(17)]{GJR18}, \cite[(19)--(20)]{GJR17} (which rest on \cite[\S2.2, Theorem~1]{Giv01} for the tangent weights and on \cite[Corollary~1]{MR2276766} for the twist), gives
```

### A4-15 (P2) GJR18 Lemma 7.7, (37): say "the series displayed on the right of (37)" [L5-F6]
Checked: `s01_intro_2.tex` l.68 "takes at $L=0$ the value proved in \cite[Lemma~7.7, (37)]{GJR18} for the
Landau--Ginzburg series there"; `s09_gamma_3.tex` l.28 "this is the formula (37) of \cite[Lemma~7.7]{GJR18} ...,
proved there for their Landau--Ginzburg integral"; `s09_gamma_4.tex` l.72 "\cite[Lemma~7.7, (37)]{GJR18} proves its
value at the orbifold point". L5: the right side of (37) at $m=5$, $k=1$, $\lambda=1$ is $\Ginf(z)$, but the
coefficients $r_{\alpha,k;l}$ of the left side are introduced in GJR18 Cor. 7.6 through $(-z)^l$, so "the value proved
there" is $\Ginf(z)$ or $\Ginf(-z)$ according to the reading. No proof depends on it. Agree, with the caveat that only
one referee read Cor. 7.6.
`s01_intro_2.tex` l.68: old `takes at $L=0$ the value proved in \cite[Lemma~7.7, (37)]{GJR18} for the Landau--Ginzburg series there` new
`takes at $L=0$ the series displayed on the right of \cite[Lemma~7.7, (37)]{GJR18}, read as a series in $z$, which is obtained there for the Landau--Ginzburg series of that paper`.
`s09_gamma_3.tex` l.28: old `this is the formula (37) of \cite[Lemma~7.7]{GJR18} at $m=5$, $k=1$ and equivariant parameter $1$ in the notation there, proved there for their Landau--Ginzburg integral` new
`this is the series displayed on the right of (37) of \cite[Lemma~7.7]{GJR18} at $m=5$, $k=1$ and equivariant parameter $1$ in the notation there, read as a series in $z$ (the coefficients on the left of (37) are introduced in \cite[Corollary~7.6]{GJR18} through powers of $-z$); it is obtained there for their Landau--Ginzburg integral`.
`s09_gamma_4.tex` l.72: `proves its value at the orbifold point` -> `gives its value at the orbifold point`.

### A4-16 (P2) conv:mod-kernel: DNOPS (2-26)-(2-29) is the recursion, not the free energy; and it is not declared [X5-F4, first bullet; L3-F5]
Checked: `s06_model.tex` l.20: "For $h\geq 2$ the free energy is defined by ... . This is the convention of
\cite[(2-26)--(2-29)]{Dunin-Barkowski:2015caa} and of \cite[(2.1)]{Hock25}." The two referees pull in opposite
directions: X5 wants (2-26)-(2-29) added to the declared conventions (the identification of the note's recursion with
the one of the quoted theorem is used in thm:mod-tr); L3 points out that (2-26)-(2-29) contain no free energy. Both
are right; the sentence is attached to the wrong object.
Old: `Here $q$ is a point of the curve. For $h\geq 2$ the free energy is defined by $(2-2h)\FTR_h:=\sum_b\Res_{\varsigma_b}\Phi_{y\dd x}\,\omega_{h,1}$, where $\dd\Phi_{y\dd x}=y\,\dd x$. This is the convention of \cite[(2-26)--(2-29)]{Dunin-Barkowski:2015caa} and of \cite[(2.1)]{Hock25}.`
New:
```latex
Here $q$ is a point of the curve. This is the recursion of \cite[(2-26)--(2-29)]{Dunin-Barkowski:2015caa}, with the kernel (2-28), to which the quoted theorem of \cref{thm:mod-tr} refers. For $h\geq 2$ the free energy is defined by $(2-2h)\FTR_h:=\sum_b\Res_{\varsigma_b}\Phi_{y\dd x}\,\omega_{h,1}$, where $\dd\Phi_{y\dd x}=y\,\dd x$; only this definition is used (it is that of \cite[(2.1)]{Hock25}; \cite{Dunin-Barkowski:2015caa} defines no free energy).
```
`s01_intro_3.tex` l.15, in the list of conventions, after `(\cref{thm:mod-var});` insert
`the recursion \cite[(2-26)--(2-29)]{Dunin-Barkowski:2015caa}, in which \cite[Theorem~4.1]{DuninBarkowski:2012bw} is recalled (\cref{conv:mod-kernel});`.
`b_signs_2.tex` l.19 is correct as it stands (it speaks of the kernel). Hock25 (2.1) was not read by the audit.

### A4-17 (P2) "the prescription [DOSS, 3.1, 3.5]" [L3-F2; X5-F4, last bullet]
Checked: `s01_intro_3.tex` l.15 "and the prescription \cite[\S3.1, \S3.5]{DuninBarkowski:2012bw}, by which the integrals
\cite[(2-31), (2-32)]{Dunin-Barkowski:2015caa} are read as formal Laplace series"; `s06_model_2.tex` l.23 gives the
reading without a reason. L3: DOSS states no rule for polar integrands; (3-57) is the termwise Gaussian rule and
(3-64) the value of the polar part. Agree. `s01_intro_3.tex` l.15, old: the clause quoted. New:
```latex
and the reading of the integrals \cite[(2-31), (2-32)]{Dunin-Barkowski:2015caa} as formal Laplace series (\cref{lem:mod-laplace,prop:mod-Rtr}): termwise Gaussian integration of the holomorphic part, as in \cite[(3-57)]{DuninBarkowski:2012bw}, and the value $\mathsf G_{-2}=-\vartheta$ for the polar term $w^{-2}\dd w$ (cf.\ \cite[(3-64)]{DuninBarkowski:2012bw})
```
`s06_model_2.tex` l.23, append to "... are read as $\sqrt{2\pi/\vartheta}\,\mathsf L_c[\omega]$.":
```latex
For the holomorphic part this is the termwise Gaussian integral. The only polar term that occurs is $w^{-2}\dd w$, in $\dd\xi^i$ at its own branch point; as $w^{-2}\re^{-\vartheta w^2/2}\dd w$ has no residue, its integral along any path avoiding $0$ is $-\sqrt{2\pi\vartheta}=\sqrt{2\pi/\vartheta}\,\mathsf G_{-2}$, the only value for which $R_{\mathrm{TR}}^{-1}=1+O(\vartheta^{-1})$.
```
(Use the note's symbol for the TR matrix at that place; E5 also edits this neighbourhood.)

### A4-18 (P2) An undeclared pointer to EO App. A inside the proof of lem:mod-xfun [X5-F4, second bullet]
Checked: `s06_model_4.tex` l.12, end: "hence $\mathrm{Tr}\,\omega=0$ by meromorphy (as in \cite[App.~A, (A-1), (A-5)]{Eynard:2007kz})."
The argument is complete without it. Agree. Old: `by meromorphy (as in \cite[App.~A, (A-1), (A-5)]{Eynard:2007kz}).`
New: `because a meromorphic form on $\bbP^1$ that vanishes on a non-empty open set vanishes (cf.\ \cite[App.~A, (A-1), (A-5)]{Eynard:2007kz}, not used).`
`s01_intro_3.tex` l.23: `(A-25);` -> `(A-1), (A-5), (A-25);`.
