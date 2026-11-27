# Suggested amendments, part 3: minor findings of units U09-U18 (file 2 of 3: P2, and P3 for Sections 6-7)

Conventions as in file 1 (`AMEND_3_units_09_18.md`). Paths are relative to `sections/`.
"Not checked against the source" means: the referee read the cited paper, I did not (no source PDF was kept).

---------------------------------------------------------------------------------------------------------------------

## P2

### A3-13 [P2] `thm:mod-var` (ii) rests on [EO, (5-15)], which the note declares unused — `s06_model_8.tex` l.31-35, l.39; `b_signs_2.tex` l.24
- Reports: U11_A F1, U11_B m1 (README "false as written, used by no proof").
- Current (l.34): "\item For $h\geq 2$: $\partial_v\FTR_h=\int_{\mathfrak r}\omega_{h,1}$." Current (l.39): "(5-14) for $g+k>1$ and
  (5-15) for $g\geq 1$". `s01_intro_3.tex` l.9 lists only "Theorem~5.1, (5-14)" in (Q), and l.23 lists "(5-15)" among the
  citations that "enter no proof".
- Judgement: agree. (ii) is cited nowhere (grep of `thm:mod-var`: only `prop:mod-V1`, part (i)); the identity it states
  is obtained independently in `thm:mod-V` (ii) (`s06_model_9.tex` l.43) from `prop:mod-period` and `thm:mod-gap`.
- Change (recommended, delete (ii)):
  - l.31-35: `Let $J=(p_1,\dots,p_n)$ be points off the branch points. \begin{enumerate} \item For $n\geq 1$ and $h+n>1$: ... \item For $h\geq 2$: ... \end{enumerate}` ->
    `Let $J=(p_1,\dots,p_n)$ be points off the branch points, $n\geq1$ and $h+n>1$. Then $\partial_v|_x\,\omega_{h,n}(J)=\int_{q\in\mathfrak r}\omega_{h,n+1}(J,q)$.`
  - l.39: `(5-14) for $g+k>1$ and (5-15) for $g\geq 1$, for the` -> `(5-14) for $g+k>1$, for the`;
    `and (5-14), (5-15) read $\delta_\Omega W^{(h)}_n(J)|_{x}=\int_{\partial\Omega}W^{(h)}_{n+1}(J,q)\Lambda(q)$, $\delta_\Omega F^{(h)}=-\int_{\partial\Omega}W^{(h)}_1\Lambda$.` ->
    `and (5-14) reads $\delta_\Omega W^{(h)}_n(J)|_{x}=\int_{\partial\Omega}W^{(h)}_{n+1}(J,q)\Lambda(q)$.`;
    `With $\omega_{h,n}=(-1)^nW^{(h)}_n$ and $\FTR_h=-F^{(h)}$ this gives $\partial_v|_x\omega_{h,n}=\ldots=\int_{\mathfrak r}\omega_{h,n+1}$ and $\partial_v\FTR_h=\int_{\mathfrak r}W^{(h)}_1\Lambda=\int_{\mathfrak r}\omega_{h,1}$.` ->
    `With $\omega_{h,n}=(-1)^nW^{(h)}_n$ this gives $\partial_v|_x\omega_{h,n}=(-1)^{n}(-1)^{n+1}(-1)\int_{\mathfrak r}\omega_{h,n+1}=\int_{\mathfrak r}\omega_{h,n+1}$.`
  - `s06_model_8.tex` l.85: `\cref{thm:mod-var}~(i)` -> `\cref{thm:mod-var}`. `b_signs_2.tex` l.24: `(5-14)--(5-15)` -> `(5-14)`.
- Alternative (keep (ii)): add "(5-15)" to (Q) in `s01_intro_3.tex` l.9 and remove it from the list at l.23.

### A3-14 [P2] `thm:mod-var`, proof: [EO, Thm. 5.1] is about the covariant variation — `s06_model_8.tex` l.39
- Reports: U11_B m2. Not checked against the source.
- Current: "This is \cite[Thm.~5.1]{Eynard:2007kz}, (5-14) \ldots, for the temperature variation of
  \cite[\S5.3, (5-27)--(5-30)]{Eynard:2007kz}, translated by \cref{conv:mod-kernel}."
- Judgement: agree that the pin should say how the quoted statement specialises; the mathematics is not in doubt
  (U11_B checked the orientation of the path).
- Change: after that sentence insert
  `(The theorem is stated there for the covariant variation $D_\Omega$ of \cite[(5-9)]{Eynard:2007kz}; on a curve of genus $0$ there are no cycle terms and $D_\Omega=\delta_\Omega$. The path of \cite[(5-27)]{Eynard:2007kz}, printed from $\alpha$ to $\beta$ with $\Lambda=1$, is taken with $\alpha=\infty$, $\beta=0$, that is, along $\mathfrak r$ from $0$ to $\infty$ with $\Lambda=-1$.)`

### A3-15 [P2] `conv:mod-kernel`: the free energy is not part of [DNOPS, (2-26)-(2-29)] — `s06_model.tex` l.20
- Reports: U09_B M3. Not checked against the source.
- Current: "For $h\geq 2$ the free energy is defined by $(2-2h)\FTR_h:=\ldots$, where $\dd\Phi_{y\dd x}=y\,\dd x$. This is the
  convention of \cite[(2-26)--(2-29)]{Dunin-Barkowski:2015caa} and of \cite[(2.1)]{Hock25}."
- Judgement: agree (the sentence follows the free energy, so "This" reads as covering it; l.21 already says
  $\FTR_h=-F^{(h)}$ of [EO, (4-10)]).
- Change: `This is the convention of` -> `The recursion \eqref{eq:mod-kernel} is in the convention of`, and append
  `; the free energy is not defined in \cite{Dunin-Barkowski:2015caa}, and ours is compared with \cite{Eynard:2007kz} below`.

### A3-16 [P2] `def:flh-outer` (ii): convergence attributed to `lem:tra-conv` (f) — `s07_comparison.tex` l.24
- Reports: U12_B m2, U12_A F7 (second half).
- Current: "so these expansions converge $u$-adically (\cref{lem:set-dilaton}, \cref{lem:tra-conv}~(f))."
- Judgement: agree. `lem:tra-conv` (f) (`a3_uconv.tex` l.29) is stated for the model correlators with
  $X_i\in\cH^u_+$; for a general $\Omega$ over $\cA$ the statement is `lem:set-dilaton` (iii) (`s02_setup_8.tex` l.11).
- Change: `(\cref{lem:set-dilaton}, \cref{lem:tra-conv}~(f))` ->
  `(\cref{lem:set-dilaton}~(iii); for the flat model correlators this is \cref{lem:tra-conv}~(f))`

### A3-17 [P2] After `def:flh-outer`: identification with the term of $O.\Omega$ needs $g\geq2$ and `lem:set-pushins` (iii) — `s07_comparison.tex` l.38
- Reports: U12_A F7 (first half).
- Current: "and $\mathrm{Cont}_\Gamma[E_O,q_O;\Omega]$ is the term of $\Gamma$ in the leg-free part of the unit-preserving action
  $O.\Omega$ (\cref{def:set-graphsum})."
- Judgement: agree. The statement that inserted translation points equal pushed-forward ones is `lem:set-pushins`
  (iii) (`s02_setup_5.tex` l.19), for $g\geq2$ and $R$ symplectic with $R(0)=1$. Not used in a proof of the unit.
- Change: `is the term of $\Gamma$ in the leg-free part of the unit-preserving action $O.\Omega$ (\cref{def:set-graphsum}).` ->
  `is, for $\Gamma$ of genus $g\geq2$, the term of $\Gamma$ in $\int_{\overline{\mathcal M}_g}(O.\Omega)_{g,0}$ for the unit-preserving action $O.\Omega$ (\cref{def:set-graphsum}, \cref{lem:set-pushins}~(iii)).`

### A3-18 [P2] `thm:flh-main`, Step 2: $E_{R''}$ is not literally the edge of [PPZ, Def. 2.2] — `s07_comparison_2.tex` l.80
- Reports: U12_B m5.
- Current: "so the two correction terms of its numerator vanish, and $E_{R''}$ is the edge of \cite[Def.~2.2]{PPZ15}, a
  symmetric power series \cite[Rem.~2.1 and (6)]{PPZ15}."
- Judgement: agree. $R''=\Pc\Psi^{\mathrm m}$ maps the colour space to $\bbH$ and $R''(0)=\Psi_{\mathrm Q}\neq1$; the edge
  meant is that of `def:flh-outer` (iii).
- Change: `and $E_{R''}$ is the edge of \cite[Def.~2.2]{PPZ15}, a symmetric power series \cite[Rem.~2.1 and (6)]{PPZ15}.` ->
  `and $E_{R''}$ is the edge $E_O$ of \cref{def:flh-outer}~(iii) for $O=R''$, a power series by the argument of \cite[Rem.~2.1 and (6)]{PPZ15} recalled there.`

### A3-19 [P2] `lem:flh-inputs` (iv): $d$ is defined in `def:tra-dressing`, not `conv:tra-forms`; vector form of $P_0$ — `s07_comparison_2.tex` l.56, l.60
- Reports: U12_A F5.
- Current (l.56): "maps $H^0$ to $dH^0$ (\cref{conv:tra-forms}: $d=(P_0)_{33}$ is the $(0,0)$ entry of the vector form)".
  Current (l.60): "$1_{\TFT}=\Psi_{\mathrm Q}^{-1}H^0=(\Psi^{\mathrm m})^{-1}P_0^{-1}H^0=d^{-1}(\Psi^{\mathrm m})^{-1}H^0$".
- Judgement: agree. $d:=(P_0)_{33}=(P_0^{\mathrm{vec}})_{00}$ is set in `s04_transport_4.tex` l.14. At l.60, $P_0$ acts
  on lowered components; the result is the same as for $P_0^{\mathrm{vec}}$ because the lowered components of $H^0$ are
  $5e_3$ and $(P_0)_{33}=d$. Wording only.
- Change (l.56): `(\cref{conv:tra-forms}: $d=(P_0)_{33}$ is the $(0,0)$ entry of the vector form)` ->
  `(\cref{def:tra-dressing}: $d=(P_0)_{33}=(P_0^{\mathrm{vec}})_{00}$, the vector form being that of \cref{conv:tra-forms})`
- Change (l.60, optional): after `$=d^{-1}(\Psi^{\mathrm m})^{-1}H^0$` add `(in lowered components $H^0$ is $5e_3$ and $(P_0)_{33}=d$)`.

### A3-20 [P2] `lem:con-key`: $v=uV$ uses `lem:con-cfu` — `s08_connection_7.tex` l.9
- Reports: U14_A m5.
- Current: "and $v^e=u^eV^e$ (\cref{def:gau-van,def:set-tau})."
- Judgement: agree. `def:set-tau` gives $v(u)=u\,\tau^{-5/4}$; that $V=\tau^{-5/4}$ is `lem:con-cfu`.
- Change: `(\cref{def:gau-van,def:set-tau})` -> `(\cref{def:gau-van}; $v(u)=u\,\tau^{-5/4}$ by \cref{def:set-tau}, and $V=\tau^{-5/4}$ by \cref{lem:con-cfu})`

### A3-21 [P2] `thm:gam-main`, proof: (iv) also depends on (ii), through (iii) — `s09_gamma_3.tex` l.35
- Reports: U15_B M5 (U15_A M1 confirms that (iii) is proved by rotation and (ii)).
- Current: "for the statements of (iv) these proofs use only (i), (iii), \cref{thm:gam-lmi} and \cref{lem:fsp-dlmf5118}".
- Judgement: agree; no circle, (ii) and (iii) are proved before (iv) is used.
- Change: `use only (i), (iii), \cref{thm:gam-lmi}` -> `use only (i), (iii) (and through it (ii)), \cref{thm:gam-lmi}`

### A3-22 [P2] `thm:gam-main` (iv): the parenthesis attributes the [GJR18] formula to `thm:gam-ginf` — `s09_gamma_3.tex` l.28
- Reports: U15_A M2, U15_B M6. Neither referee nor I compared (37) of [GJR18, Lemma 7.7] with the source.
- Current: "this is the formula (37) of \cite[Lemma~7.7]{GJR18} at $m=5$, $k=1$ and equivariant parameter $1$ in the notation
  there, proved there for their Landau--Ginzburg integral (\cref{thm:gam-ginf})."
- Judgement: agree on the wording. The comparison enters no proof (`s09_gamma_4.tex` l.12 says so).
- Change: `proved there for their Landau--Ginzburg integral (\cref{thm:gam-ginf}).` ->
  `proved there for their Landau--Ginzburg integral; the value $c(z;0)=\Ginf(z)$ is proved here independently (\cref{thm:gam-ginf}), and the comparison with \cite{GJR18} enters no proof.`
  The pin itself (m, k, parameter) remains to be read against the PDF by the author.

### A3-23 [P2] Introduction: $\mathrm{Cont}_\Gamma$ "as in [PPZ15, Defs. 2.2, 2.13]" for a general $X$ — `s01_intro_2.tex` l.7
- Reports: U18_A m5.
- Current: "with the leaf point $q_X$ of $X$, as in \cite[Defs.~2.2, 2.13]{PPZ15}, including the factor $|\Aut\Gamma|^{-1}$
  (\cref{def:set-graphsum,def:flh-outer})."
- Judgement: agree. Theorem A uses $X=R'$, with $R'(0)\neq1$ and symplecticity not known beforehand; its edge is that of
  `def:set-outer`, which is the PPZ edge only for symplectic $X$ with $X(0)=1$.
- Change: `as in \cite[Defs.~2.2, 2.13]{PPZ15}, including the factor $|\Aut\Gamma|^{-1}$ (\cref{def:set-graphsum,def:flh-outer})` ->
  `as in \cite[Defs.~2.2, 2.13]{PPZ15} when $X$ is symplectic with $X(0)=1$, and with the edge of \cref{def:set-outer} in general, including the factor $|\Aut\Gamma|^{-1}$ (\cref{def:set-graphsum,def:set-outer,def:flh-outer})`

### A3-24 [P2] The two lists of PPZ items disagree — `s10_main_2.tex` l.36 against `s01_intro_3.tex` l.9
- Reports: U18_A m9. (The clause on Props. 2.4, 2.8, 2.9 in both lists is rewritten by E6; do this together.)
- Current (l.36): "\cite[\S0.2, (1); Def.~0.3; (5); Rem.~2.1; Defs.~2.2 and~2.5, \ldots]{PPZ15}"; (Q) has
  "Def.~0.3 (iii)" and "Remark~2.1 and (6)".
- Judgement: agree; harmless.
- Change (l.36): `Def.~0.3; (5); Rem.~2.1;` -> `Def.~0.3~(iii); (5); Rem.~2.1 and (6);`

### A3-25 [P2] `cor:mai-univ`: existence of $P^A_{\mathrm{can}}$ referred to `def:gau-can`, which is stated for the quintic pair — `s10_main_4.tex` l.21
- Reports: U18_A m10 (second point), U18_B F5 (first point).
- Current: "which exists, is unique, is graded and has rational coefficients, as in \cref{def:gau-can} (\cref{thm:gau-U},
  \cref{prop:gau-graded})".
- Judgement: agree, mildly: the two general statements are already named, but "as in \cref{def:gau-can}" routes through
  `thm:gau-G`, stated for the quintic pair. The item numbers below were read in `s03_gauge_4.tex` l.38-42 and
  `s03_gauge_2.tex` l.70-75. (The use of `thm:gau-G` (c) in the proof, l.31, for a general pair was not examined.)
- Change: `as in \cref{def:gau-can} (\cref{thm:gau-U}, \cref{prop:gau-graded})` ->
  `by the argument of \cref{def:gau-can}: a rational graded regular gauge $P$ exists (\cref{prop:gau-graded}~(i)), $P_{11}(0,z)$ is a unit of $\bbQ\fps{z^5}$ (\cref{prop:gau-graded}~(iii)), so $P\,P_{11}(0,z)^{-1}$ is a graded regular gauge (\cref{prop:gau-graded}~(ii)), and it is the only regular gauge with this normalisation and has rational coefficients (\cref{thm:gau-U}~(b), (d))`

---------------------------------------------------------------------------------------------------------------------

## P3 (Sections 6 and 7)

### A3-26 [P3] $\kappa$ "a rational number" — `s06_model_2.tex` l.2
- Reports: U09_A F2, U09_B M2. Agree: `lem:mod-scale` (i) is stated for "constants $\kappa,\mu\neq0$" and applied in (ii) with
  $\kappa=k^4$, $k\in\bbC^\times$.
- Change: `$\kappa\neq 0$ is a rational number and` -> `$\kappa\neq 0$ is a complex constant and`

### A3-27 [P3] "applied here, and only here, to the curve $x=\Wv(p)$, $y=p$" — `s06_model.tex` l.5
- Reports: U09_B M5. Agree: the quoted identification is applied to $(\bbP^1,4\Wv,p)$ (`s06_model_3.tex` l.28).
  ((Q) in `s01_intro_3.tex` l.9 is rewritten by E5.)
- Change: `to the curve $x=\Wv(p)$, $y=p$ on $\bbP^1$, $v\neq 0$;` ->
  `to the curve $x=\Wv(p)$, $y=p$ on $\bbP^1$, $v\neq 0$, and to its rescalings $(\kappa\Wv,\mu p)$ (\cref{lem:mod-scale});`

### A3-28 [P3] Free energy of a general curve; $f\in\bbQ(v)[X]$ at a fixed $v$ — `s06_model.tex` l.20, `s06_model_4.tex` l.4
- Reports: U09_A F4. Agree; only the model curve is used.
- Change (l.20): `For $h\geq 2$ the free energy is defined by` ->
  `For $h\geq 2$, on a curve where $\omega_{h,1}$ has zero residues at the branch points, so that the choice of $\Phi_{y\dd x}$ does not matter (for the model curve, \cref{lem:mod-xfun}), the free energy is defined by`
- Change (`s06_model_4.tex` l.4): `(for instance $f\in\bbQ(v)[X]$, $v\neq 0$ fixed)` ->
  `(for instance $f\in\bbQ(v)[X]$ with coefficients regular at the fixed $v\neq0$)`

### A3-29 [P3] Complex $v$ against the formal root of `def:set-modcohft` — `s06_model_3.tex` l.32
- Reports: U09_A F5. Agree.
- Change: `($\hat\Delta_j^{1/2}=\varsigma_j$, \cref{def:set-modcohft})` ->
  `($\hat\Delta_j^{1/2}=\varsigma_j$, \cref{def:set-modcohft}, the formal root $\varsigma$ there being specialised to the fifth root of $v$ fixed in \cref{def:mod-curve})`

### A3-30 [P3] "By (i), $1/\varsigma_k=\zeta^k/\varsigma$" — `s06_model_3.tex` l.19
- Reports: U09_B M8. Agree: this is `def:mod-curve` ($\varsigma_k=\zeta^{-k}\varsigma$); (i) is what gives the left side.
- Change: `By (i), $1/\varsigma_k=\zeta^k/\varsigma$ and $\sum_k\zeta^{-kb}=5\delta_{b,0}$, the left side` ->
  `By (i), and as $1/\varsigma_k=\zeta^k/\varsigma$ (\cref{def:mod-curve}) and $\sum_k\zeta^{-kb}=5\delta_{b,0}$, the left side`

### A3-31 [P3] Codomain "$\{\text{legs}\}$" of the leg map — `s06_model_5.tex` l.4
- Reports: U10_A m-1. Agree.
- Change: `$\Lambda\colon\mathbb F\to\{\text{legs}\}$,` -> `$\Lambda\colon\mathbb F\to\mathsf M$ (the module defined below; \cref{lem:mod-closure}~(iii)),`

### A3-32 [P3] `thm:mod-Lreg`: pairing of legs with $\bbQ(v)$-coefficients — `s06_model_7.tex` l.40
- Reports: U10_A m-2. Agree; downstream only $\bbQ[v][p]+\bbQ[v][x]$ occurs.
- Change: `the residue pairing lies in $\bbQ[v]$; the components` ->
  `the residue pairing, defined for $v$ off the poles of the coefficients, is the restriction of an element of $\bbQ[v]$; the components`

### A3-33 [P3] "$(\cL/\bbQ(v)[x])_w$": the quotient is graded only through `lem:mod-classL` (iv) — `s06_model_6.tex` l.43
- Reports: U10_A m-4, U10_B m2. Agree.
- Change: `is bijective, and both sides have dimension $\mathrm d(w)$;` ->
  `is bijective, and both sides have dimension $\mathrm d(w)$; here $\cL/\bbQ(v)[x]$ is graded through its identification with $\bbQ[p,\bar p]/\bbQ[v,x]$ (\cref{lem:mod-classL}~(iv));`

### A3-34 [P3] "$\mathcal R^b_n$ \ldots a monomial in $v$" can be $0$ — `s06_model_5.tex` l.26, l.30
- Reports: U10_A m-7. Agree.
- Change (both places): `a monomial in $v$` -> `a rational multiple of a power of $v$, possibly $0$`

### A3-35 [P3] `prop:mod-period`, Step 1: index set of the finite sum — `s06_model_7.tex` l.61
- Reports: U10_B m1. Agree: $r+1=2h+j$, so $j$ is fixed by $r$.
- Change: `the finite sum over $1\leq j\leq h$ and $k_1+\dots+k_j=h$ with $r=2h-1+j$` ->
  `the finite sum over the ordered $j$-tuples of positive integers with $k_1+\dots+k_j=h$, where $j:=r-2h+1$ (the sum is empty unless $1\leq j\leq h$)`

### A3-36 [P3] `cor:mod-onepoint`: "its $v$-exponent" when the correlator vanishes — `s06_model_10.tex` l.28
- Reports: U11_A F5. Agree.
- Change: `If $\langle X\rangle_h$ is singular, that is, if its $v$-exponent is negative, then` ->
  `If the $v$-exponent $e=1-2h+w/5$ of $\langle X\rangle_h$ (\cref{lem:set-modhom}) is negative, then`

### A3-37 [P3] `lem:mod-fixedx`: "every family of meromorphic forms" — `s06_model_8.tex` l.9
- Reports: U11_B m4. Agree.
- Change: `For every family of meromorphic forms,` -> `For every family of meromorphic forms depending holomorphically on $v$,`

### A3-38 [P3] `rmk:mod-offK`: the clause on two-point correlators is not a statement — `s06_model_10.tex` l.37
- Reports: U11_B m5. Agree ("$\Fm_h(V(q))$" has no defined meaning here). Unused.
- Change: `is in general not $\Pi_0(X)(\Fm_h)'$, nor are those of two-point correlators given by $\Fm_h(V(q))$ for a genus-independent $V$; both facts are known only by computation` ->
  `is in general not $\Pi_0(X)(\Fm_h)'$; this is known only by computation`

### A3-39 [P3] What the proof of `thm:flh-main` uses; (unit) for a series not known to be symplectic — `s07_comparison.tex` l.38
- Reports: U12_A F1, F3. Agree with both: the proof uses `lem:flh-inputs` (ii) and the exact symplecticity of
  $\Pc\Psi^{\mathrm m}$, and (unit) involves only the invertible constant term.
- Change: `satisfies (unit) with $\varkappa=1$, and` ->
  `satisfies (unit) with $\varkappa=1$ (the condition involves only the constant term, the isometry $\Psi_{\mathrm Q}$, and makes sense without symplecticity), and`
- Change: `under $\Fg{g}$ it is symplectic modulo $z^{3g-2}$ (\cref{lem:flh-inputs}~(iii)), which is all that the proof of \cref{thm:flh-main} uses.` ->
  `under $\Fg{g}$ it agrees modulo $z^{3g-2}$ with the symplectic series $\Pc\Psi^{\mathrm m}$ (\cref{lem:flh-inputs}~(ii), (iii)), which is all that the proof of \cref{thm:flh-main} uses.`

### A3-40 [P3] $\Psi_{\mathrm Q}$ denotes both the matrix over $\Nring\otimes\bbQ(\zeta)$ and its expansion — `s07_comparison_2.tex` l.56
- Reports: U12_A F4 (also U18_A m6, U18_B F3; `s10_main.tex` l.35 acknowledges it). Agree; optional.
- Change: `So, as an identity in $\Nring\otimes\bbQ(\zeta)$, $(\Psi_{\mathrm Q})_{aj}=` ->
  `So the representative in $\Nring\otimes\bbQ(\zeta)$ of $(\Psi_{\mathrm Q})_{aj}$ (denoted $\Psi_{\mathrm Q}$ in \cref{lem:mai-Psi}, its expansion being $\Psi_{\mathrm{cf}}$ there) is`

### A3-41 [P3] `rmk:flh-fail` (i), (ii) — `s07_comparison_2.tex` l.100
- Reports: U12_B m6, U12_A F8. Agree; remarks only. In (i) the scalar case forces $\cn=\Dn_0\in\bbQ\fps{z^5}$.
- Change (i): `If $\Dn=\cn(z)1$ is scalar with $\cn(z)\cn(-z)=1$, then` ->
  `If $\Dn=\cn(z)1$ is scalar (then $\cn=\Dn_0$ lies in $\bbQ\fps{z^5}$, by $\Dn_i(z)=\Dn_0(\zeta^{-i}z)$ and \cref{thm:con-glue}) with $\cn(z)\cn(-z)=1$, then`
- Change (ii): `but different graph terms.` -> `but in general different graph terms (no example is recorded here).`
