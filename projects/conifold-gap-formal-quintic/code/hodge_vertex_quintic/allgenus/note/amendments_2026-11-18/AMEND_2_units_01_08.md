# Suggested amendments, part 2: minor findings of units U01-U08 (file 1 of 2: P1 and P2)

Source: `audit_2026-11-18/reports/` U01_C, U02_C, U03_A, U03B, U04_A, U04_B, U05_A, U05B, U06_A, U06_B, U07_A, U07B,
U08_A, U08_B. Findings ruled on in `verify/E1..E8.md` or in the README lists "Must fix" / "Routine missing steps" are
skipped (list at the end of file 2). Two items of the README list "False as written, used by no proof" that sit in
these units (`rmk:set-qrrsign`, `cor:gau-cycsys` (ii)) are included as A2-01, A2-02; drop them if part 1 has them.

Method: duplicates between lenses merged; every place was opened in the current TeX (paths relative to `sections/`,
line numbers are those of the current files, one paragraph per line); "Current" quotes the text now on the page;
"Change" is old -> new in the note's macros. No `.tex` file was edited. Nothing was compiled.

Classes: P1 = statement or proof step wrong or incomplete as written; P2 = citation, pin or cross-reference;
P3 = wording or notation, optional; R = rejected. This file: P1 (A2-01..11), P2 (A2-12..20). File 2: P3, R, skipped.

---------------------------------------------------------------------------------------------------------------------

## P1

### A2-01 [P1] `rmk:set-qrrsign`, last sentence false for the rows — `s02_setup_4.tex` l.52
- Reports: U01_C F1, U05B F5 (also L5; README "false as written").
- Current: "That $\mathrm{ph}_1=-\frac3{20}$ is forced, independently of (QRR), is part of \cref{lem:gau-HPF}."
- Judgement: agree. `lem:gau-HPF` (iii) assumes $\mathrm{ph}_n\in\bbQ\fps u$. The $\mathrm{ph}_1$ of `def:set-rows` lies in
  $\bbQ\fps q$, where $(4-5\theta_u)\mathrm{ph}_1=-\frac35$ has the solutions $-\frac3{20}+c\,u^{4/5}$, and
  $u^{4/5}=(1-5^5q)^{4/5}\in\bbQ\fps q$ (checked: $\theta_uu^{4/5}=\frac45u^{4/5}$). No proof uses the remark.
- Change: replace the sentence by
  `For a solution of \eqref{eq:gau-HPF} with $\mathrm{ph}_0=1$ and coefficients in $\bbQ\fps u$, such as $\Phican$, the value $\mathrm{ph}_1=-\frac3{20}$ is forced (\cref{lem:gau-HPF}~(iii)). For the rows, where $\mathrm{ph}_1\in\bbQ\fps q$ is fixed by the equation only up to a multiple of $u^{4/5}=(1-5^5q)^{4/5}$, the value at $q=0$ comes from (QRR); so $-5\zeta^i$ is the sign compatible with $\Dn_0=1$ (\cref{thm:gam-D}).`

### A2-02 [P1] `cor:gau-cycsys` (ii), literal reading false — `s03_gauge.tex` l.82
- Reports: U03B M4 (README "false as written").
- Current: "the same holds with $z\Dq$ replaced by $z\Dq+w$ in \eqref{eq:gau-system} and in $\PFfive$."
- Judgement: agree. Read literally the recursion $Y_{a+1}=z\Dq Y_a/\kappa_a$ of (i) stays unshifted, and that
  equivalence is false. The proof (`a2_cyclicPF4.tex` l.122, automorphism $\Dq\mapsto\Dq+w/z$ applied to the whole
  chain) proves the version with the recursion shifted, which is the one used (`prop:gau-quintic` (i), `thm:gau-R0` (ii)).
- Change: `... replaced by $z\Dq+w$ in \eqref{eq:gau-system} and in $\PFfive$.` ->
  `... replaced by $z\Dq+w$ in \eqref{eq:gau-system}, in $\PFfive$ and in the recursion $Y_{a+1}=z\Dq Y_a/\kappa_a$.`

### A2-03 [P1] `prop:gau-quintic` (i): the closing-relation clause asserts what is not yet proved — `s03_gauge_2.tex` l.20, l.31
- Reports: U04_A M1, U04_B Finding 1, U03B M5 (three reports, one finding).
- Current (l.20): "... to the relations $(z\Dq+\zeta^iL)\Cexp(X_a)=\kappa_a\Cexp(X_{a+1})$, $L=s^{-1}$; the closing relation is then given by \cref{cor:gau-cycsys}."
  Current (l.31): "... it carries each forward relation to the relation with $\kappa_a$; the closing relation is \cref{cor:gau-cycsys}~(ii) with $w=\zeta^is^{-1}$."
- Judgement: agree. `cor:gau-cycsys` (ii) is an equivalence; for arbitrary $X_a$ with the forward relations the closing
  relation fails in general, and for the rows the needed $\PFfive$-equation of the expansion is available only later.
  All later uses quote (i) for $a\leq3$ (`s03_gauge_8.tex` l.30 "for $a\leq3$", l.43), so nothing downstream changes.
- Change (l.20): `; the closing relation is then given by \cref{cor:gau-cycsys}.` ->
  `; by \cref{cor:gau-cycsys}~(ii) with $w=\zeta^is^{-1}$, the closing relation $(z\Dq+\zeta^iL)\Cexp(X_4)=\kappa_4\Cexp(X_0)$ then holds if and only if $\varpi_0\Cexp(X_0)$ is annihilated by $\PFfive$ with $z\Dq$ replaced by $z\Dq+\zeta^is^{-1}$. For the rows this is proved in \cref{thm:gau-R0}~(iv) (colour $0$) and \cref{cor:gau-colours} (all colours).`
- Change (l.31): `; the closing relation is \cref{cor:gau-cycsys}~(ii) with $w=\zeta^is^{-1}$.` ->
  `; the statement on the closing relation is \cref{cor:gau-cycsys}~(ii) with $w=\zeta^is^{-1}$, $\cR=\bbQ(\zeta)\lser s\lser z$ and $Y_a=\Cexp(X_a)$ ($\kappa_4=\varpi_0=y_0$).`

### A2-04 [P1] Stripping the formal exponential: from "$\re^{u_i/z}zL\Phi_i$ solves $\PFfive$" to $\PFfive^{(i)}(L\Phi_i)=0$ — `s02_setup_4.tex` l.24; used at `s03_gauge_7.tex` l.62, `s03_gauge_8.tex` l.28
- Reports: U01_C F7, U06_A F1, U06_B F4.
- Current (l.24): "For the WKB form, $v=R(z)\tilde e_i\re^{u_i/z}$, with the formal exponential $\Dq\re^{u_i/z}=z^{-1}(\Dq u_i)\re^{u_i/z}$, solves $z\partial_Tv=H*v$, so $C_0v_0=\zeta^i\re^{u_i/z}L\Phi_i$ solves $\PFfive$ by \cref{cor:set-cyclicPF5}~(b); multiply by the constant $z\zeta^{-i}$."
  Current (`s03_gauge_7.tex` l.62): "By \cref{cor:set-cyclicPF5}~(b), $L\Phi_i$ is a colour-$i$ WKB amplitude in the sense of \cref{lem:gau-wkb}".
- Judgement: agree, incomplete as written. `lem:gau-wkb` defines an amplitude by $\PFfive^{(i)}\cA=0$; the cited
  statement is about a product with an exponential in an unnamed algebra (a power series in $z$ times a series in
  $z^{-1}$ with $\logq$). The step closes in two lines once the algebra is named (I re-checked the computation below).
- Change (l.24): replace the quoted sentence by
  `For the WKB form let $\mathsf B\coloneqq\bbQ(\zeta)\fps q\lser z[\mathsf E,\mathsf E^{-1}]$ with $\Dq\mathsf E\coloneqq z^{-1}\zeta^iL\,\mathsf E$ (a local letter: $\mathsf E$ stands for $\re^{u_i/z}$, $\Dq u_i=\zeta^iL$), a differential $\bbQ\fps q$-algebra in which $z$ is a constant. By \cref{def:set-rows}~(b), $v\coloneqq R(z)\tilde e_i\,\mathsf E$ solves $z\partial_Tv=H*v$ over $\mathsf B$, as $z\partial_T(w\mathsf E)=(z\partial_Tw+\mu_iw)\mathsf E$; so $C_0v_0=\zeta^iL\Phi_i\,\mathsf E$ solves $\PFfive$ by \cref{cor:set-cyclicPF5}~(b); multiply by the constant $z\zeta^{-i}$. Moreover $\Dq(f\mathsf E)=\bigl((\Dq+z^{-1}\zeta^iL)f\bigr)\mathsf E$ gives $\PFfive(f\mathsf E)=(\PFfive^{(i)}f)\,\mathsf E$ for $f\in\bbQ(\zeta)\fps q\lser z$, where $\PFfive^{(i)}$ is $\PFfive$ with $z\Dq$ replaced by $z\Dq+\zeta^iL$; as $\mathsf E$ is a unit of $\mathsf B$, $\PFfive^{(i)}(L\Phi_i)=0$.`
- Change (`s03_gauge_7.tex` l.62): `By \cref{cor:set-cyclicPF5}~(b), $L\Phi_i$ is` -> `By \cref{cor:set-cyclicPF5}~(b), in the form $\PFfive^{(i)}(L\Phi_i)=0$ obtained in the proof after \cref{def:set-rows}, $L\Phi_i$ is`.
- Change (`s03_gauge_8.tex` l.28): `$A^{\mathrm{act}}$ is a colour-$0$ WKB amplitude by \cref{cor:set-cyclicPF5}~(b),` -> `$A^{\mathrm{act}}=L\Phi_0$ is a colour-$0$ WKB amplitude by \cref{cor:set-cyclicPF5}~(b) in the form $\PFfive^{(0)}(L\Phi_0)=0$ (proof after \cref{def:set-rows}),`.

### A2-05 [P1] `lem:gau-conc` (iv) uses $\Conc(f,\varpi_1)=0$ before it is stated — `a2_cyclicPF4.tex` l.50
- Reports: U03_A M1, U03B M1.
- Current: "then the constant $\Conc(f,y_3)=\Conc(f,\varpi_1)\ell+\sum_{a,b}\Dq^af\,\Jet_{ab}Y_b$ is the coefficient of $u^0$ of the second term."
- Judgement: agree. The vanishing is stated only in `lem:gau-iso` (ii). No circularity: it follows on the spot.
- Change: replace by
  `then $\Conc(f,y_3)=\Conc(f,\varpi_1)\ell+\sum_{a,b}\Dq^af\,\Jet_{ab}Y_b$. The left side lies in $\bbQ$ by (ii), and the right side has degree at most $1$ in the indeterminate $\ell$ over $\bbK$; hence $\Conc(f,\varpi_1)=0$ for $f\in\mathcal S\cap\bbK$, and the second term is a constant of $\bbK$, equal to its coefficient of $u^0$.`
  (`\mathcal S` is the solution space of `def:gau-conc`, as written in this appendix.)

### A2-06 [P1] `lem:gau-HPF` (ii): the "if and only if" omits the coefficient of $\zh^0$ — `s03_gauge_6.tex` l.65
- Reports: U05_A F2, U05B F2.
- Current: "depends only on $\mathrm{ph}_0,\dots,\mathrm{ph}_{n-1}$. So \eqref{eq:gau-HPF} holds if and only if $(4n-5\theta_u)\mathrm{ph}_n=-r_n$ for all $n\geq0$"
- Judgement: agree. (ii) computes only the coefficients of $\zh^{n+1}$, $n\geq0$. The $\zh^0$ coefficient is
  $\Phi^{[5]}_0-\mathrm{ph}_0=0$ by the recursion of step (4) at $m=0$; cancelling $1-u$ is allowed in $\bbQ\lser s$.
- Change: `depends only on $\mathrm{ph}_0,\dots,\mathrm{ph}_{n-1}$. So \eqref{eq:gau-HPF} holds` ->
  `depends only on $\mathrm{ph}_0,\dots,\mathrm{ph}_{n-1}$; the coefficient of $\zh^0$ of the left side is $\Phi^{[5]}_0-\mathrm{ph}_0=0$, as every $\mathscr O_e$ preserves the coefficient of $\zh^0$. So, $1-u$ being non-zero in the field $\bbQ\lser s$, \eqref{eq:gau-HPF} holds`
- Optional (P3, same lemma, l.71-73): the proof steps are numbered (1)-(5) while the items are (i)-(iii) ("using (2)",
  "(5) (iii) follows"). Rename the steps `\emph{Step 1.}` ... `\emph{Step 5.}`.

### A2-07 [P1] `lem:tra-twoslot`, (i) => (ii): summability of the regrouped families is not what the first paragraph shows — `s04_transport_2.tex` l.77; `s04_transport.tex` l.30
- Reports: U07_A M2, U07B N6.
- Current (l.77): "regrouping the double family over $(\alpha,(l,b))$, summable by the first paragraph, gives $E=\sum_\alpha\fr_\alpha(\psi)\otimes G_\alpha(\psi')$".
  Current (`s04_transport.tex` l.30): "Every regrouping of a summable family (a partition into subfamilies) consists of summable families, and the sum of their sums is the total sum, because coefficientwise this is a statement about finite sums".
- Judgement: agree. The first paragraph treats $(E^{\alpha\beta}\fr_\alpha\otimes\fr_\beta)$, not the family with a
  flat second slot; and `lem:tra-sum` is then applied to $(\fr_\alpha\otimes G_\alpha)_\alpha$, whose summability
  (a family of sub-sums) is not stated. Both close at once.
- Change (`s04_transport.tex` l.30): `consists of summable families, and the sum of their sums is the total sum, because coefficientwise this is a statement about finite sums` ->
  `consists of summable families; the family of their sums is summable (for a fixed coefficient index and $N$ the finitely many members outside $u^N\bbQ\fps{u}$ lie in finitely many subfamilies, and every other sub-sum lies in $u^N\bbQ\fps{u}$); and the sum of their sums is the total sum, because coefficientwise and modulo $u^N$ this is a statement about finite sums`
- Change (l.77): `regrouping the double family over $(\alpha,(l,b))$, summable by the first paragraph, gives` ->
  `regrouping the double family $\bigl(E^{ab}_{kl}\,\fr_\alpha(\psi)\otimes\psi'^lH^b\bigr)$ over $(\alpha,(l,b))$, summable by the computation of the first paragraph in the first slot (the second slot is a flat monomial), gives the summable family $(\fr_\alpha\otimes G_\alpha)_\alpha$ (\cref{conv:tra-forms}) and`

### A2-08 [P1] `thm:red-main`, Step 3: the unit $d^{2g_w-2+n_w}$ of the flat vertex tensor is not covered by (b4) — `s05_reduction.tex` l.73
- Reports: U08_B m5 (digest U08-14).
- Current: "By \cref{lem:tra-conv}~(b4), applied to the edges $E_P$, ... and to the leaves $dT_{\geq 2}\in z^2\cH^u_+\cap\sK$, the contraction along $\Gamma$ is the sum of the summable family ..."
- Judgement: agree. (b4) has "a vertex factor as in (b3)", i.e. $\sum_m\frac1{m!}\langle X_1\cdots X_ny^m\rangle_h$, while
  `def:red-flatsum` has $d^{2g_w-2+n_w}$ times it. Harmless: the product of the units is $d^{2g-2}$
  ($\sum_w(2g_w-2+n_w)=2g-2$, checked).
- Change: insert after `... and to the leaves $dT_{\geq 2}\in z^2\cH^u_+\cap\sK$,`:
  `the vertex tensors of \cref{def:red-flatsum} being the vertex factors of \cref{lem:tra-conv}~(b3) with $y=dT_{\geq2}$ times the units $d^{\,2g_w-2+n_w}$ of $\bbQ\fps{u}$, of product $d^{\,2g-2}$, which multiplies the contraction and every member of the family,`

### A2-09 [P1] `lem:tra-conv` (b4), (e2) invoke (b2)/(b3) outside the standing $n\geq1$ of (b) — `a3_uconv.tex` l.21, l.24
- Reports: U08_A m2, U08_B m4.
- Current (l.21): "(b4) Let $\Gamma$ be a leg-free stable graph. ... at each vertex $w$ a vertex factor as in (b3)".
  Current (l.24): "The family $(\frac{\delta^j}{l!\,j!}\langle(y^\circ)^lH^j\rangle_g)_{l,j\geq 0}$, every member expanded in frames by (b2), is summable".
- Judgement: agree (scope of the statement only). (b) assumes $n\geq1$; the vertex graph has $n=0$ and so has the member
  $(l,j)=(0,0)$. (b4) is used once, for edge graphs (`s05_reduction.tex` l.73; grep confirms no other use).
- Change (l.21): `(b4) Let $\Gamma$ be a leg-free stable graph.` -> `(b4) Let $\Gamma$ be a leg-free stable graph with at least one edge, so that every vertex carries $n_w\geq1$ half-edges (\cref{lem:red-genus}).`
- Change (l.24): `every member expanded in frames by (b2)` -> `every member with $(l,j)\neq(0,0)$ expanded in frames by (b2)`

### A2-10 [P1] `def:red-ML`: the "in short" form is false for $h=0,1$ in low order — `s05_reduction.tex` l.39
- Reports: U08_A m3, U08_B m2.
- Current: "; in short, $\overline F_h|_{-z+\sK}=\Fm_h\circ\Pi_0+\mathrm{Reg}_h$ with a regular remainder $\mathrm{Reg}_h$."
- Judgement: agree. For $h=0$ the Taylor coefficients of order $n\leq2$ of $\overline F_0$ vanish while those of
  $\Fm_0\circ\Pi_0$ contain $\log v$; likewise order $0$ for $h=1$. The precise sentence before it is right and is
  the one used.
- Change: append to the sentence, before the full stop:
  `, an identity of the Taylor coefficients at the base point of order $n>2-2h$ (all orders for $h\geq2$, $n\geq1$ for $h=1$, $n\geq3$ for $h=0$)`

### A2-11 [P1] Section 2, opening: "proved later, independently of this section" is false — `s02_setup.tex` l.14
- Reports: U01_C F4, U06_A F6, U06_B F2.
- Current: "Four facts used in the definitions are proved later, independently of this section: \cref{lem:gau-confexp}, \cref{cor:gau-rowsNN}, \cref{lem:gau-colrows}, and the symplecticity of $R'$ (from \cref{thm:gam-D}, whose proof does not use $R'$)."
- Judgement: agree (a sentence about the logical order, not a proof step). The four facts use Section 2
  (`def:set-rows`, `prop:set-cyclicLR`, `cor:set-cyclicPF5`, `thm:set-qrr`, `lem:set-hrot`, ...). What holds is
  non-circularity.
- Change: `are proved later, independently of this section:` -> `are proved later, by arguments that do not use the definitions in which they are invoked, so that no circle arises:`

---------------------------------------------------------------------------------------------------------------------

## P2

### A2-12 [P2] Pairing attributed to `def:set-frob` — `s04_transport.tex` l.15
- Reports: U01_C F9, U07_A M1, U07B N5.
- Current: "and the pairing $\eta(H^a,H^b) = 5\,\delta_{a+b\equiv 3}$ of \cref{def:set-frob}". Agree: the formula is in `def:set-fq`.
- Change: `of \cref{def:set-frob}` -> `of \cref{def:set-fq}`

### A2-13 [P2] References to items that are not in `def:tra-tower` — `s04_transport_2.tex` l.26, l.28
- Reports: U07_A M1, U07B N5.
- Current (l.26): "$\Pi_i(q)\in\bbQ v^{(w-4i)/5}$ (\cref{def:tra-tower})"; (l.28): "at $i = k$ by (ii) of \cref{def:tra-tower}".
- Judgement: agree. Items (i)-(iii) and the homogeneity statement are in the unnumbered paragraphs after the definition
  (`s04_transport.tex` l.54-72).
- Change (l.26): `(\cref{def:tra-tower})` -> `(the paragraph on homogeneity after \cref{def:tra-tower})`
- Change (l.28): `by (ii) of \cref{def:tra-tower}` -> `by the value (ii) after \cref{def:tra-tower}`

### A2-14 [P2] Order bound attributed to `lem:tra-frames` (c) — `s04_transport_2.tex` l.75
- Reports: U07_A M1, U07B N5.
- Current: "for $a' = 1$ the $u$-order is at least $(4(k-k')+1)/5$ (\cref{lem:tra-frames}~(c))". Agree: it is the exponent of (b)
  with $\nu_a\geq0$ and $\ord_uv\geq1$.
- Change: `(\cref{lem:tra-frames}~(c))` -> `(the exponent of \cref{lem:tra-frames}~(b), with $\nu_a\geq0$ and $\ord_uv\geq1$)`

### A2-15 [P2] Definition of $\overline F_g$ cited to `lem:set-dilaton` — `a3_uconv.tex` l.57
- Reports: U08_B m3.
- Current: "which is $\overline F_g(\tilde q)$ by the definition of $\overline F_g$ (\cref{lem:set-dilaton})". Agree.
- Change: `(\cref{lem:set-dilaton})` -> `(\cref{def:set-graphsum,def:set-modcohft})`

### A2-16 [P2] Dependency sentences say "only" but are incomplete — `a3_uconv.tex` l.35; `s05_reduction.tex` l.7
- Reports: U08_A m1, U08_B m1.
- Current (a3 l.35): "The parts depend on the following statements only: (b) and (e) on the homogeneity of \cref{lem:set-modhom}; ..."
  Current (s05 l.7): "The proof uses only the results of \cref{sec:tra}, no gauge equation and no datum of the quintic beyond $T_P$."
- Judgement: agree. (a), (b2)-(b4), (e) use `lem:tra-sum`, `lem:tra-frames`, `cor:tra-slice`, `lem:tra-twoslot`; (c) uses
  `lem:gau-vanshape`. I checked by grep that `s04_transport.tex` and `s04_transport_2.tex` do not use `lem:tra-conv`
  in a proof (l.42 only lists it), so there is no circle.
- Change (a3 l.35): `The parts depend on the following statements only:` ->
  `Beyond \cref{conv:tra-forms}, \cref{lem:tra-sum,lem:tra-frames,cor:tra-slice,lem:tra-twoslot} (which precede \cref{sec:tra-dress} and do not use this lemma), \cref{lem:gau-vanshape} in (c) and the definitions of \cref{def:set-graphsum,def:set-modcohft}, the parts depend on the following statements only:`
- Change (s05 l.7): `uses only the results of \cref{sec:tra},` -> `uses only the results of \cref{sec:tra} and \cref{lem:tra-conv},`

### A2-17 [P2] `def:gau-can`: existence cites a statement that presupposes a graded regular gauge — `s03_gauge_6.tex` l.25
- Reports: U05B F4.
- Current: "It exists and is unique: $\Pc=P\,P_{11}(0,z)^{-1}$ for any graded regular gauge $P$ (\cref{thm:gau-G}~(a)),". Agree.
- Change: `for any graded regular gauge $P$ (\cref{thm:gau-G}~(a)),` -> `for any graded regular gauge $P$ (one exists by \cref{prop:gau-graded}~(i); \cref{thm:gau-G}~(a)),`

### A2-18 [P2] `lem:set-modhom` (c): step (1) of `lem:set-M0sympl` proves a neighbouring statement — `s02_setup_7.tex` l.65
- Reports: U02_C F5.
- Current: "Replacing $\varsigma$ by $\zeta^{-1}\varsigma$ replaces, by step~(1) of the proof of \cref{lem:set-M0sympl}, the colour-$i$ column of \cref{def:set-M0} by the colour-$(i+1)$ column".
- Judgement: agree, wording. Step (1) (l.24) applies $\varsigma\mapsto\zeta^{-i}\varsigma$ to colour $0$; the claim is the
  same computation applied to colour $i$.
- Change: `replaces, by step~(1) of the proof of \cref{lem:set-M0sympl},` -> `replaces, by the computation of step~(1) of the proof of \cref{lem:set-M0sympl} ($U_i\mapsto U_{i+1}$, $t\mapsto\zeta^{-1}t$, $\varsigma^{\nu_b-1}\mapsto\zeta^{b+1}\varsigma^{\nu_b-1}$),`

### A2-19 [P2] `def:set-outer`: "leg-free part of $R'.(M.\TFT)$ in the sense of `def:flh-outer`" — `s02_setup_8.tex` l.50
- Reports: U02_C F4.
- Current: "$E_{R'}$ is the edge of \cite[Def.~2.2]{PPZ15} \cite[Rem.~2.1, (6)]{PPZ15}, and the outer graph sum is the leg-free part of $R'.(M.\TFT)$ in the sense of \cref{def:flh-outer}."
- Judgement: agree in part. The report says `def:flh-outer` defines nothing for a map between two modules; it does
  (item (ii): edge $E_O$ and leaf point $q_O$ for $O\colon\mathbb V\to\mathbb V'$ with (unit)). What it does not define
  is an action "$R'.$" or a "leg-free part" when $\mathbb V'\neq\mathbb V$ (the paragraph after it treats only
  $\mathbb V'=\mathbb V$, $O(0)=1$). With $X\coloneqq R'(0)^{-1}R'$ one has $E_{R'}=E_X$ and $q_{R'}=q_X$ (checked:
  $R'(0)$ is an isometry mapping $1_{\TFT}$ to $1$). No proof uses the sentence.
- Change: `and the outer graph sum is the leg-free part of $R'.(M.\TFT)$ in the sense of \cref{def:flh-outer}.` ->
  `and $E_{R'}$, $q_{R'}$ are the edge and the leaf point of \cref{def:flh-outer}~(ii) for $O=R'$; with $X\coloneqq R'(0)^{-1}R'$ (symplectic, $X(0)=1$) they are $E_X$ and $q_X$, so the outer graph sum is the leg-free part of the unit-preserving action $X.(M.\TFT)$ on the colour space (\cref{lem:set-pushins}~(iii)).`

### A2-20 [P2] Pointer "the equation of $\Phi_0$ is `lem:gau-HPF`" — `s02_setup_4.tex` l.27
- Reports: U05B F5 (second half).
- Current: "the equation of $\Phi_0$ is \cref{lem:gau-HPF}."
- Judgement: agree, minor. `lem:gau-HPF` is stated for $\mathrm{ph}_n\in\bbQ\lser s$; for $\Phi_0$ (coefficients in
  $\bbQ(\zeta)\fps q$) items (i), (ii) carry over, item (iii) does not (see A2-01). Used in no proof.
- Change: `the equation of $\Phi_0$ is \cref{lem:gau-HPF}.` -> `the equation of $\Phi_0$ is that of \cref{lem:gau-HPF}~(i), (ii), which hold in the same way for coefficients in $\bbQ(\zeta)\fps q$; item (iii) there concerns coefficients in $\bbQ\fps u$ only.`

(Continued in `AMEND_2_units_01_08_2.md`: P3, R, skipped findings, what was not verified.)
