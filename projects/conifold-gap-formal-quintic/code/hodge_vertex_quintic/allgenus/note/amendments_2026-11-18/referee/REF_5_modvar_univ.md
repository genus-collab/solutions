# REF_5: `thm:mod-var` (ii) and `cor:mai-univ` after the amendments of 2026-11-18

Adversarial referee report, read-only on the sources. Files read: `sections/s06_model.tex` (l.10-46),
`s06_model_2.tex` (lem:mod-scale), `s06_model_7.tex`, `s06_model_8.tex`, `s06_model_9.tex` (l.30-50),
`s06_model_10.tex`, `s03_gauge_2..6.tex`, `s04_transport_4.tex` (l.1-116), `s05_reduction.tex`, `s10_main_4.tex`,
the archived `s10_main_4.tex` and `s06_model_8.tex`, and the amendment items A4-02, A3-14, A3-25, A3-62.

| piece | verdict |
|---|---|
| 1. `thm:mod-var` (ii) and its proof | **SOUND** (two wording items, both in the recollection of EO in part (i)) |
| 2. `cor:mai-univ`, `lem:mai-pole` | **SOUND AFTER FIXES** (one sign slip introduced by A3-62 in the paragraph after the corollary; wording) |

---------------------------------------------------------------------------------------------------------------

## Piece 1. `thm:mod-var` (ii), `s06_model_8.tex` l.34 and l.41

### 1.1 What was checked

Statement: for $h\geq2$, $\partial_v\FTR_h=\int_{\mathfrak r}\omega_{h,1}$.

Proof, re-derived.
- `prop:mod-period` (iii): $\FTR_h=\frac{B_{2h}}{2h(2h-2)}v^{2-2h}$, so
  $\partial_v\FTR_h=(2-2h)\frac{B_{2h}}{2h(2h-2)}v^{1-2h}=-\frac{B_{2h}}{2h}v^{1-2h}$. Correct.
- `prop:mod-period` (ii): $\int_\infty^0\omega_{h,1}=\frac{B_{2h}}{2h}v^{1-2h}$ on any path off the branch points.
  $\mathfrak r$ runs from $0$ to $\infty$ (`def:mod-pair`, and the bracket $[\,\cdot\,]_{q=0}^{q=\infty}$ in
  `lem:mod-fixedx`), so $\int_{\mathfrak r}=-\int_\infty^0$. Correct sign; the two sides agree.
- Genus range. $h\geq2$ is what both inputs need ((iii) is stated for $h\geq2$; $\FTR_h$ is defined in
  `conv:mod-kernel` only for $h\geq2$). $h=0,1$ are rightly excluded from (ii): there is no $\FTR_0,\FTR_1$ in the
  convention. (With the ad hoc $\FTR_1=-\frac1{12}\log v$ of `thm:mod-V` (ii) the identity also holds at $h=1$:
  `prop:mod-period` (ii) gives $\int_{\mathfrak r}\omega_{1,1}=-B_2/(2v)=-1/(12v)$, the value found independently
  in the proof of `thm:mod-V`. Nothing to change.)

No circularity. `prop:mod-period` (l.50-74 of `s06_model_7.tex`) is proved from `eq:mod-710`, `def:mod-xi`,
`thm:mod-Lreg`, `lem:mod-scale` (ii), `def:mod-pair` and the definition of $\FTR_h$; all precede `thm:mod-var`,
and none cites it. Step 5 was re-derived: $\Phi_{y\dd x}=p^5/5-v\log p$ ($y\,\dd x=(p^4-v/p)\dd p$),
$\langle p^5\rangle_h$ has exponent $2-2h-1+1<0$ and lies in $\bbQ[v]$, hence vanishes,
and $\FTR_h=\frac{-v}{2-2h}\cdot\frac{B_{2h}}{2h}v^{1-2h}$. Correct.

"This part is not used": true. All citations of `thm:mod-var` in `main.tex`, `sections/*.tex`:
`s06_model_8.tex` l.87 (proof of `prop:mod-V1`, explicitly part (i)); `s01_intro_3.tex` l.15 (the set-up
(5-1), (5-4), (5-13) of EO Thm 5.1, i.e. part (i)); `b_signs_2.tex` l.24 (table row, (5-14), part (i)).
`thm:mod-V` (ii), $h\geq2$, $k=1$, which is the same identity, is proved separately from the definition of
$\FTR_h$ and homogeneity (`s06_model_9.tex` l.43), not from `thm:mod-var` (ii). `s01_intro_3.tex` l.23 lists
(5-15) among the citations that enter no proof: now true.

Part (i) is unaffected: the edit only split off the reference to (5-15). Its sign chain was re-derived:
$\Omega=-\partial_v(y\,\dd x)|_x=\dd p/p=\int_{\mathfrak r}B(\cdot,q)\Lambda(q)$ with $\Lambda=-1$;
$\partial_v|_x\omega_{h,n}=(-1)^n\cdot(-1)\cdot(-1)^{n+1}\int_{\mathfrak r}\omega_{h,n+1}=+\int_{\mathfrak r}\omega_{h,n+1}$.
The closing parenthesis of (ii), "it is also the translation of (5-15)", is correct:
(5-15) reads $D_\Omega F^{(g)}=-\int_{\partial\Omega}W^{(g)}_1\Lambda$, and with $\FTR_h=-F^{(h)}$,
$\omega_{h,1}=-W^{(h)}_1$, $\Lambda=-1$ this is $\partial_v\FTR_h=+\int_{\mathfrak r}\omega_{h,1}$.

### 1.2 The pins of EO

Read: the stored TeX source `audit_2026-11-18/work/L2/src/x/invariance3.tex` (Section 5 of the paper, complete)
and the text extract `audit_2026-11-18/work/L2/dl/eo_v4.txt` (for the printed equation numbers).
Not read: the page images `L2/img/*.png`; the appendix proof of Thm 5.1 (`appvariation.tex`) beyond its first lines.

Confirmed: (5-1) $\delta_\Omega y|_x\,\dd x=-\Omega$; (5-13) $\Omega(p)=\int_{\partial\Omega}B(p,q)\Lambda(q)$ with
the footnote on hard edges; Theorem 5.1 with (5-14) for $g+k>1$ and (5-15) for $g\geq1$, both for the covariant
variation $D_\Omega=\delta_\Omega+\tr(\kappa\,\delta_\Omega\tau\,\kappa\,\partial/\partial\kappa)$, which is
$\delta_\Omega$ in genus $0$; (5-27) $\Omega=-\dd S_{\alpha,\alpha'}=\int_\alpha^{\alpha'}B(p,q)$,
$\partial\Omega=[\alpha,\alpha']$, $\Lambda=1$. With $\alpha=\infty$, $\alpha'=0$:
$\int_\infty^0B(p,q)=\dd p/p=\Omega$. So the sentence added by A3-14 is mathematically right.

### 1.3 Findings

**F1.1 (wording).** `s06_model_8.tex` l.39. The paper's letters in (5-27) are $\alpha$ and $\alpha'$, not
$\alpha$ and $\beta$ ($\beta$ occurs in (5-28) as a running index).
- old: `printed from $\alpha$ to $\beta$ with $\Lambda=1$, is taken with $\alpha=\infty$, $\beta=0$, that is,`
- new: `printed from $\alpha$ to $\alpha'$ with $\Lambda=1$ (letters of that paper), is taken with $\alpha=\infty$, $\alpha'=0$, that is,`

**F1.2 (wording).** `s06_model_8.tex` l.39, and the same pin in `s01_intro_3.tex` l.15. The hypothesis
"$\Omega/\dd y$ has no pole at the branch points" is the sentence printed between (5-3) and (5-4); (5-4) is the
Rauch formula derived under it.
- old: `is regular at the branch points \cite[(5-4)]{Eynard:2007kz};`
- new: `is regular at the branch points \cite[the assumption before (5-4)]{Eynard:2007kz};`

**F1.3 (wording, optional).** l.29, title "The variational formula, in the form used", while (ii) opens with
"This part is not used". Either leave it, or
- old: `\item For $h\geq 2$: $\partial_v\FTR_h=\int_{\mathfrak r}\omega_{h,1}$.`
- new: `\item (Not used.) For $h\geq 2$: $\partial_v\FTR_h=\int_{\mathfrak r}\omega_{h,1}$.`

**Remark (no change required).** (ii) does not need the Bernoulli numbers: by `lem:mod-scale` (ii),
$\FTR_h$ is homogeneous of degree $2-2h$ in $v$, and $(2-2h)\FTR_h=v\langle-\log p\rangle^{\mathrm{TR}}_h$ as in
Step 5 of `prop:mod-period`; hence $\partial_v\FTR_h=(2-2h)\FTR_h/v=\int_{\mathfrak r}\omega_{h,1}$. This is the
argument of `thm:mod-V` (ii); (ii) duplicates that statement at $k=1$. The present proof is correct as it stands.

### 1.4 Exact computation

Script: `<user home>\AppData\Local\Temp\claude\d7c1b608d03d194a\37f24460-3601-4c85-b123-5156e19eb475\scratchpad\ref5\period.py`
(Python, `fractions`, Laurent series at $p=0$). It evaluates $\Phi^{\mathrm A}_h=\sum_{r\geq1}D_x^{r-1}(a_r/x')$
of Step 1 of `prop:mod-period` directly from `eq:mod-710`, without Steps 3-4, for $h=1,\dots,5$ and $v=1,3$.

| $h$ | $\Phi^{\mathrm A}_h(0)$ at $v=1$ | $-B_{2h}/(2h)$ | at $v=3$ | no negative powers of $p$ |
|---|---|---|---|---|
| 1 | $-1/12$ | $-1/12$ | $-1/36$ | yes |
| 2 | $1/120$ | $1/120$ | $1/3240$ | yes |
| 3 | $-1/252$ | $-1/252$ | $-1/61236$ | yes |
| 4 | $1/240$ | $1/240$ | $1/524880$ | yes |
| 5 | $-1/132$ | $-1/132$ | $-1/2598156$ | yes |

So $\int_{\mathfrak r}\omega_{h,1}=\Phi^{\mathrm A}_h(0)=-\frac{B_{2h}}{2h}v^{1-2h}=\partial_v\FTR_h$ for
$h=2,\dots,5$, with the $v$-scaling $v^{1-2h}$. The constants $c_h=(-5)^{1-h}B_{2h}/(2h(2h-2))$ come out as
$1/1200$, $1/25200$, $1/180000$, $1/660000$ for $h=2,\dots,5$.
Scope of the test: it checks `prop:mod-period` (ii) and the consistency of (ii), (iii) with `thm:mod-var` (ii);
it takes `eq:mod-710` as recalled in the note and is not an independent run of the recursion `eq:mod-kernel`
at genus $\geq2$ (not attempted).

---------------------------------------------------------------------------------------------------------------

## Piece 2. `cor:mai-univ` and `lem:mai-pole`, `s10_main_4.tex` l.6-32

### 2.1 The residue formula

Self-duality (`def:gau-onepole`): $b_0=b_2$, $b_3=b_4$, so $\Bprod=b_0^2b_1b_3^2$ and $b_1=\Bprod/(b_0^2b_3^2)$.
$r:=(ub_1)(0)=(u\Bprod)(0)/(b_0(0)b_3(0))^2$, and $(u\Bprod)(0)=-(uv_A'^5/v_A)(0)=-v_A'(0)^4$
(this is `lem:mai-pole` (ii) with $y=v_A$, $n=1$). Hence $r=-v_A'(0)^4/(b_0(0)b_3(0))^2=-1/\lambda^2$ with
$\lambda=b_0(0)b_3(0)/v_A'(0)^2$. Correct. What it needs is all in the hypotheses: $b_0,b_3$ are units of
$\bbQ\fps u$ (one-pole shape), so $b_0(0)b_3(0)\in\bbQ^\times$; $v_A'(0)\in\bbQ^\times$; the residue is that of
`def:gau-onepole`. $\lambda\in\bbQ^\times$ holds. (The hypothesis "residue $r$" is now only a name; fine.)

### 2.2 Do the cited results apply to a general such $A$?

- $\tilde A_{v_A}$ has the one-pole shape for every $v_A$ with $v_A'(0)\neq0$, residue $\tilde r=-1$, is
  self-dual, with product $-v_A'^5/v_A$ (the computation of `prop:gau-equalprod` (i), which does not use
  $v'(0)=1$ except for the values $\tilde b_a(0)$). So (EP) holds by hypothesis.
- `thm:gau-U` and `prop:gau-graded` assume exactly: one-pole shape for both systems and (EP), over a field of
  characteristic $0$; rationality in (d) and in `prop:gau-graded` (i) needs all $b_a,\tilde b_a\in\bbQ\lser u$.
  I read the proofs (`lem:gau-comm`, `lem:gau-lift`, `lem:gau-scal`, `lem:gau-lead`, `lem:gau-degree`): none
  uses the normalisation (N0), $r=\tilde r$, self-duality, or any datum of the quintic. All hypotheses are
  implied by those of the corollary. Nothing is missing.
- The construction in the statement: a rational graded regular gauge $P$ exists (`prop:gau-graded` (i));
  $P_{11}(0,z)\in\bbQ\fps{z^5}^\times$ ((iii)); $PP_{11}(0,z)^{-1}$ is graded regular ((ii)), rational;
  uniqueness since regular gauges are $P\cn(z)$ and $(P\cn)_{11}(0,z)=P_{11}(0,z)\cn(z)$ (`thm:gau-U` (b)). Correct.
- `lem:gau-R1`, `lem:gau-vanq` (a first part, b, c, d), `lem:gau-uniq`, `thm:gau-G` (c) are stated and proved
  for general one-pole systems over $\bbQ\fps u$. (A3-25 had left `thm:gau-G` (c) unexamined; it is fine.)
- `prop:gau-symp` with $r\neq\tilde r$: re-derived. $S:=P(-z)^T\eta P(z)$ gives $\eta^{-1}S=t(z)1$ by
  `thm:gau-U` (b) for $(\tilde A,\tilde A)$; $S_{21}(0,z)=5P_{22}(0,-z)P_{11}(0,z)$ and
  $rP_{22}(0,z)=\tilde rP_{11}(0,z)$, so $t=(\tilde r/r)P_{11}(0,z)P_{11}(0,-z)$. For $P_A$:
  $\tilde r/r=\lambda^2$, $P_{11}(0,z)=\lambda^{-1}$, $t=1$. Correct.
- (A) and $T_{P_A}=T_A/\lambda$: from `thm:gau-G` (c) and `lem:gau-vanq` (c), (d). Correct.
- `thm:red-main` is stated for a dressing at any $v(u)$ with $v(0)=0$, $v'(0)\neq0$, with (Sp), (A), and does
  not use the gauge equation; `def:red-flatsum` is in the homogeneous form valid for $d(0)\neq1$. With
  `thm:mod-MLK` ($\ML h$ for all $h\geq0$) it gives (i) and (ii) for every $g\geq2$. Correct.

Consistency checks done by hand for $P_A$ (leading coefficients from `lem:gau-diag`, $\rho^U_a$):
$p_{0,0}p_{0,3}=1$ at $u=0$ (needed by (Sp) at $z^0$); $d\,T_{P_A}\equiv v_A$ and
$T_{P_A}\equiv p_{0,0}v_A$ modulo $u^2$ (`lem:tra-T` (c), (d)). All three hold identically in
$b_0(0),b_3(0),v_A'(0)$.

`lem:mai-pole`: statement and proof correct ($A^5=\Bprod\,1$, $\tilde A_y^5=-(y'^5/y)1$,
$\ord_0(y'^5/y)=4n-5$, $\Res_{u=0}\Bprod=-y'(0)^4$ for $n=1$). The corollary uses only the residue computation of
(ii); "compare" undersells it but is not wrong.

### 2.3 Comparison with the archived text

Changes: (a) "over $\bbQ\fps u$" added; (b) "assume $r=-1/\lambda^2$ with $\lambda\in\bbQ^\times$" replaced by the
derived formula with $\lambda$ now *defined*; (c) the existence argument spelled out; (d) in the paragraph after
the corollary, the sentence on the system $(b_0,\Bprod/(b_0^2b_3^2),b_0,b_3,b_3)$ with
"$\lambda=\pm b_0(0)b_3(0)/v_A'(0)^2$" replaced by "Replacing $\lambda$ by $-\lambda$ changes nothing: the vertex
law is even in $\lambda$". The proof is unchanged.
Lost: the freedom of the sign of $\lambda$, which the old text used silently for the quintic ($\lambda=1$). See F2.1.

### 2.4 Findings

**F2.1 (error, minor; introduced by A3-62).** `s10_main_4.tex` l.21 and l.28. With the new definition, the
formal quintic has $b_0(0)=1$, $b_3(0)=-1$, $v'(0)=1$ (`prop:gau-quintic` (ii)), so $\lambda=-1$ and
$P_A=-\Pc$, not $\Pc$. The sentence "the corollary is \cref{cor:red-can} with \cref{thm:mod-MLK}" is then not
literally true, and the justification offered for the sign, "the vertex law is even in $\lambda$", covers only
the formula in (ii), not the identity of the flat sums of $P$ and $-P$ (which is true: $R$, $E_P$,
$\tilde q_T$ are unchanged, $d\mapsto-d$, each vertex tensor changes by $(-1)^{n_v}$ and $\sum_vn_v$ is even).
The pulled-back model itself also has $b_0(0)b_3(0)/v_A'(0)^2=-1$. Fix: put the sign into the definition.
- l.21 old: `with $\lambda\coloneqq b_0(0)b_3(0)/v_A'(0)^2\in\bbQ^\times$`
- l.21 new: `with $\lambda\coloneqq-b_0(0)b_3(0)/v_A'(0)^2\in\bbQ^\times$`
- l.28 old: `$r=-1$ and $T_A=\tc$ (\cref{lem:gau-vanq}~(d)), and the corollary is \cref{cor:red-can} with \cref{thm:mod-MLK}. Replacing $\lambda$ by $-\lambda$ changes nothing: the vertex law is even in $\lambda$.`
- l.28 new: `$r=-1$, $\lambda=1$ ($b_0(0)=1$, $b_3(0)=-1$, $v'(0)=1$), $P_A=\Pc$ and $T_A=\tc$ (\cref{lem:gau-vanq}~(d)), and the corollary is \cref{cor:red-can} with \cref{thm:mod-MLK}. The sign of $\lambda$ is a convention: the proof uses only $\lambda^2=-1/r$, and the flat sums of $P$ and of $-P$ agree term by term, since $R$, $E_P$ and $\tilde q_T$ are unchanged, $d$ changes sign, and the vertex tensors of \cref{def:red-flatsum} change by $(-1)^{n_v}$ with $\sum_vn_v$ even.`

  (The proof at l.31 needs no change: it uses $P_{11}(0,z)=\lambda^{-1}$ and $\tilde r/r=\lambda^2$ only.)

**F2.2 (wording).** l.21: the reason given for $r=-1/\lambda^2$ names only one of its two ingredients, and
`lem:mai-pole` (ii) is the other, not a comparison.
- old: `because $b_1=\Bprod/(b_0^2b_3^2)$ (compare \cref{lem:mai-pole})`
- new: `because $b_1=\Bprod/(b_0^2b_3^2)$ by self-duality and $(u\Bprod)(0)=-v_A'(0)^4$ (\cref{lem:mai-pole}~(ii) with $y=v_A$)`

**F2.3 (wording).** l.31: `lem:gau-uniq` needs $\YC(v_A(u),z)\in\bbQ\lser z\fps u^5$, which is
`lem:gau-vanshape` (valid for every $v\in u\bbQ\fps u$, though worded for $v(u)$).
- old: `(\cref{lem:gau-vanmod}, the chain rule and \cref{lem:gau-uniq}, as in the proof of \cref{thm:gau-G})`
- new: `(\cref{lem:gau-vanmod}, the chain rule, \cref{lem:gau-vanshape}, whose last claim holds for every $v\in u\bbQ\fps u$, and \cref{lem:gau-uniq}, as in the proof of \cref{thm:gau-G})`

**F2.4 (wording).** l.21: `thm:gau-U` (d) is not needed once everything is over $\bbQ$ (a regular gauge is by
definition over the ground field); harmless. Also "self-dual one-pole system" is the note's shorthand for
"self-dual cyclic system with one simple pole" of `def:gau-onepole`; used already at l.4. No change required.

### 2.5 Conclusion for Piece 2

The new statement is correct as mathematics for either sign of $\lambda$: residue formula, existence and
uniqueness of $P^A_{\mathrm{can}}$, symplecticity of $P_A$, (A), $T_{P_A}=T_A/\lambda$, and the conclusion for
every $g\geq2$ all follow from the hypotheses as cited. The only defect is the sign convention of $\lambda$
against the quintic specialisation in the following paragraph (F2.1), with its incomplete justification.
