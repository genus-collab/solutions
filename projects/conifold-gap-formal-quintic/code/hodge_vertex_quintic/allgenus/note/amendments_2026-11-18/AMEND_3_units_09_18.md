# Suggested amendments, part 3: minor findings of units U09-U18 (file 1 of 3: P1)

Source: `audit_2026-11-18/reports/` U09_A, U09_B, U10_A, U10_B, U11_A, U11_B, U12_A, U12_B, U13_A, U13_B, U14_A,
U14_B, U15_A, U15_B, U16_C, U17_C, U1617B, U18_A, U18_B. Findings ruled on in `verify/E1..E8.md` or in the README
lists "Must fix" / "Routine missing steps" are skipped (table at the end of file 3). Items of the README list "False
as written, used by no proof" that sit in these units (`thm:mod-var` (ii); `rmk:mai-graphlocal`, `rmk:mai-necessity`,
`s01_intro_2.tex` l.25) are included as A3-13 and A3-11; drop them if another part has them.

Method: duplicates between the lenses merged; every place was opened in the current TeX (paths relative to
`sections/`; line numbers are those of the current files, one paragraph per line); "Current" quotes the text now on
the page; "Change" is old -> new in the note's macros. No `.tex` file was edited. Nothing was compiled or computed.

Classes: P1 = statement or proof step wrong or incomplete as written (even if the repair is one phrase);
P2 = citation, pin or cross-reference; P3 = wording or notation, optional; R = rejected.
This file: P1 (A3-01..12). File 2: P2 (A3-13..25) and P3 (A3-26..41). File 3: P3 (A3-42..63), R, skipped.

---------------------------------------------------------------------------------------------------------------------

## P1

### A3-01 [P1] `prop:mod-LD`: $q_\phi$ and "TR leg" used for legs outside the class of `def:mod-ancvec` — `s06_model_4.tex` l.34
- Reports: U09_A F1, U09_B M4.
- Current (`s06_model.tex` l.59): "For a leg $\phi\in\cL$, or $\phi=-\log p$, or a $\bbQ[v]$-linear combination of these,
  put $q_\phi(z):=\ldots$". Current (l.34): "Hence $q_{\phi^b_k}=-\varrho_b\,z^kH^b$, and the TR leg of $v^mz^kH^b$ (the
  insertion $v^m\psi^kH^b$, \cref{conv:tra-forms}) is $-\varrho_b^{-1}v^m\phi^b_k$."
- Judgement: agree. For instance $\phi^1_1=-x\log p+p^4/16-v/p$ is not in $\cL=\bbQ[v][p]+\bbQ(v)[x]$, so $q_{\phi^1_1}$
  is not defined. Only the one-slot pairing with $\dd\xi_{a,d}$ is needed, and \eqref{eq:mod-LD} defines it. The uses are
  inside the class: in `thm:mod-MLK` Step 1 the legs are $-\Lambda(Y_i)\in\bbQ[p,\bar p]\subset\cL$ and $-\log p$.
- Change (l.34): `Hence $q_{\phi^b_k}=-\varrho_b\,z^kH^b$, and the TR leg of` ->
  `For the legs $\phi^b_k$ and their $\bbQ[v]$-linear combinations, which need not lie in the class of \cref{def:mod-ancvec}, we define $q_\phi(z):=\sum_{a,d}\langle\phi,\dd\xi_{a,d}\rangle H^az^d$ by the pairing \eqref{eq:mod-LD}; where both definitions apply they agree (for $-\log p$ by the last display of \cref{def:mod-pair}). Hence $q_{\phi^b_k}=-\varrho_b\,z^kH^b$, and, in this sense, the TR leg of`
  and add at the end of the statement:
  `\Cref{thm:mod-tr} is applied only to legs in $\cL$ and to $-\log p$ (proof of \cref{thm:mod-MLK}).`

### A3-02 [P1] `prop:mod-Rtr`, Step 3: `prop:set-modrows` (b) is stated over $\bbQ(\zeta)$ and applied over $\bbC$ — `s06_model_3.tex` l.13
- Reports: U09_A F3.
- Current: "By the homogeneity above, $F_b=\varepsilon_bJ_{\nu_b}(1;\varsigma^4\vartheta)$ is a power series in $t$, and the
  computation of \cref{prop:set-modrows}~(b) shows that $(F_b)_b$ solves the $T$-system of colour $j$."
- Judgement: agree. `prop:set-modrows` (b) (`s02_setup_6.tex` l.66) takes $F_b\in\bbQ(\zeta)\fps t$; here the
  coefficients contain $(-5\kappa)^{-1/2}$ and a complex $\varsigma$. The computation and the pivot argument are formal.
- Change: after that sentence insert
  `(Here $F_b\in\bbC\fps t$. The computation and the uniqueness statement of \cref{prop:set-modrows}~(b) are formal and hold verbatim over $\bbC$ in place of $\bbQ(\zeta)$: the recursion has the pivot $4n\zeta^{4j}\neq0$.)`

### A3-03 [P1] `def:mod-pair`, `prop:mod-V1`: iterated ray integrals need a joint bound — `s06_model.tex` l.41, `s06_model_8.tex` l.85
- Reports: U09_B M7, U11_A F2 and F3, U11_B m3.
- Current (l.41): "These integrals converge and do not depend on the rays: in each variable $\omega_{h,n}$ is holomorphic
  off the branch points (the paragraph after \cref{conv:mod-kernel}), and its residues there vanish".
  Current (l.85): "and the integrand is holomorphic in $v$, bounded on the rays by a constant times
  $|\dd q|/(1+|q|)^2$ uniformly near $v_0$."
- Judgement: agree. Holomorphy "in each variable" does not give convergence of an iterated integral, the exchange of
  the order, or differentiation under the integral sign; the bound at l.85 is stated for one ray variable and without
  a reason. Both close with the form of `lem:mod-abdks` (proved in `s06_model_7.tex` without pairings).
- Change (l.41): `These integrals converge and do not depend on the rays: in each variable $\omega_{h,n}$ is holomorphic off the branch points (the paragraph after \cref{conv:mod-kernel}), and its residues there vanish` ->
  `These integrals converge, may be taken in any order, and do not depend on the rays. Indeed $\omega_{h,n}$ is a meromorphic multidifferential on $(\bbP^1)^n$ whose polar divisor lies in $\bigcup_{i,b}\{p_i=\varsigma_b\}$: another component would give, for generic values of the other variables, a pole in one variable off the branch points, against the paragraph after \cref{conv:mod-kernel}. So $\omega_{h,n}$ is jointly holomorphic on products of rays and near their end points, hence bounded there by a constant times $\prod_i|\dd p_i|/(1+|p_i|)^2$; in each variable its residues at the branch points vanish`
- Change (l.85): `bounded on the rays by a constant times $|\dd q|/(1+|q|)^2$ uniformly near $v_0$.` ->
  `bounded on the rays by a constant times $\prod_j|\dd q_j|/(1+|q_j|)^2$ uniformly near $v_0$: by \cref{lem:mod-abdks}, $\omega_{h,n}=N\prod_i\dd p_i/\prod_i(p_i^5-v)^K$ with $N\in\bbQ[v][p_1,\dots,p_n]$, and holomorphy at $p_i=\infty$ for every $v\neq0$ gives $\deg_{p_i}N\leq5K-2$. The same bound justifies the exchanges of integrals (Fubini) below.`

### A3-04 [P1] `lem:mod-modv`: two steps of "the map" are not on the page — `s06_model_6.tex` l.34
- Reports: U10_A m-3, U10_B m3.
- Current: "as $\bar p\equiv -p^4/4$ modulo $\bbQ[v,x]$, $\phi^2_k\equiv(-1)^{k+1}p^{4k+4}/\big(4\prod_{i=1}^k(4+4i)\big)$; and
  by \cref{lem:mod-closure}~(ii) the $\log p$-component of $\phi^1_k$ is $(-1)^kx^k/k!$ modulo $v$."
- Judgement: agree. The congruence $\bar p\equiv-p^4/4$ holds for $\phi^2_0$ only; its propagation through $I$ is not
  argued, and for $\phi^1_k$ closure (ii) treats only the part $\log p\cdot g(x)$. I checked the replacement below by
  hand from `lem:mod-closure` (i), (ii): $I(\bar p)=-vp^3/3-\bar p^2/2$, $\phi^2_1\equiv p^8/32$, in agreement with the
  printed formula ($4\prod_{i=1}^k(4+4i)=4^{k+1}(k+1)!$), and $\phi^1_1=-x\log p+p^4/16-\bar p$.
- Change: replace the quoted text by
  `again by \cref{lem:mod-closure}~(i), $I(\bar p^c)\equiv-\bar p^{c+1}/(c+1)$ mod $v\mathsf M$ for $c\geq1$ (the other term is $-vp^3\bar p^{c-1}/(4-c)$, or $-v^4\log p$ if $c=4$), so $\phi^2_k\equiv(-1)^k\bar p^{k+1}/(k+1)!$ mod $v\mathsf M$, which is log-free; and $x^n\equiv p^{4n}/4^n+\bar p^n$ modulo $p\bar p\,\bbQ[p,\bar p]$ gives $\bar p^n\equiv-p^{4n}/4^n$ in $\mathsf N/v\mathsf N$, hence $\phi^2_k\equiv(-1)^{k+1}p^{4k+4}/\big(4\prod_{i=1}^k(4+4i)\big)$. The $\log p$-component of $\phi^1_k$ is $(-1)^kx^k/k!$ modulo $v$, by induction on $k$: for the part $\log p\cdot g(x)$ this is \cref{lem:mod-closure}~(ii) with $G$ without constant term (the terms $v^c\log p$ of $\int G\,\dd p/p$ then have $c\geq1$), and by \cref{lem:mod-closure}~(i) $I$ produces $\log p$ from $\bbQ[p,\bar p]$ only with a factor $v$.`

### A3-05 [P1] `lem:mod-ibp` (i): the radius of the circle is not quantified — `s06_model_8.tex` l.20
- Reports: U11_A F4 (U11_B m4 for the unused hypothesis).
- Current: "Then $(2\pi\ri)^{-1}\oint_{|p|=R}\phi\,\mathrm{Lie}_{\cX}w$ equals $-\langle(\cX\dd\phi)_+,w\rangle$."
- Judgement: agree. The proof ("has all its poles inside the circle") needs $R>|\varsigma|$; for a smaller circle the
  statement is false. The hypothesis "zero residues" is not used in (i) or (ii); it is harmless and may stay.
- Change: `Then $(2\pi\ri)^{-1}\oint_{|p|=R}` -> `Then, for $R>|\varsigma|$, $(2\pi\ri)^{-1}\oint_{|p|=R}`

### A3-06 [P1] `def:flh-outer` (i): the symmetry required of $E$ is not the one the edges have — `s07_comparison.tex` l.23, l.33
- Reports: U12_B m1.
- Current (l.23): "a bivector series $E\in(\mathbb V\otimes\mathbb V)\fps{\psi,\psi'}$ that is symmetric under the
  exchange of its two slots". Current (l.33): "so $E_O$ is a symmetric power series in $\psi,\psi'$".
- Judgement: agree. $E_O$, $E_{R'}$, $E_{\Pc}$ satisfy $E(\psi,\psi')^T=E(\psi',\psi)$ (`s04_transport_4.tex` l.33;
  `def:set-outer`, `s02_setup_8.tex` l.50, says "together with $\psi$ and $\psi'$"), not symmetry of the tensor alone.
- Change (l.23): `that is symmetric under the exchange of its two slots` ->
  `that is symmetric under the exchange of its two slots together with $\psi$ and $\psi'$, $E(\psi,\psi')^T=E(\psi',\psi)$`
- Change (l.33): `so $E_O$ is a symmetric power series in $\psi,\psi'$` ->
  `so $E_O$ is a power series in $\psi,\psi'$, symmetric in the sense of (i)`

### A3-07 [P1] `def:flh-outer` (iii): "thus $O^{-1}(z)=O(-z)^*$" needs equal ranks — `s07_comparison.tex` l.25
- Reports: U12_A F2.
- Current: "Let $(\mathbb V',\eta',1')$ be a second such module and $O(z)=\sum_{k\geq0}O_kz^k$, \ldots a symplectic series:
  $O(-z)^*O(z)=1$, \ldots; thus $O^{-1}(z)=O(-z)^*$."
- Judgement: agree. $O(-z)^*O(z)=1$ gives a left inverse only; it is two-sided when the ranks agree. All uses have rank 5.
- Change: `be a second such module and` -> `be a second such module, of the same rank as $\mathbb V$, and`

### A3-08 [P1] `prop:gam-rational`, step (5): the order-1 coefficient of $\mathsf S_0$ is not rational — `s08_connection_3.tex` l.50
- Reports: U13_A m3, U13_B m2.
- Current: "(5) The Taylor coefficients of order at least $1$ of $\mathsf S_0$, $\mathsf S_1$, $\mathsf S_{2k}$ are rational
  in $\sigma$ and carry no branch".
- Judgement: agree. $\mathsf S_0'=\lambda+5\log(5\sigma)-\sum_j\log(\zeta^j-\sigma)$ carries the branches; it vanishes at
  $L$ only by the choice of $\lambda$, as the next sentence of the proof itself says. Step (2) uses orders $\geq3$ only.
- Change: `of order at least $1$ of $\mathsf S_0$, $\mathsf S_1$, $\mathsf S_{2k}$ are rational` ->
  `of order at least $2$ of $\mathsf S_0$, and of order at least $1$ of $\mathsf S_1$ and $\mathsf S_{2k}$, are rational`

### A3-09 [P1] `lem:con-id` (a), (b): the square-root rules of the full symbols are not reconciled — `s08_connection_6.tex` l.35
- Reports: U14_A m2, U14_B m1.
- Current (end of the proof of (b)): "The exponents $m=-1,\dots,3$ are the $\nu_b$; insert (a) and divide by
  $N_{\zh}=(2\pi\zh/5)^{1/2}\re^{-\hat\cW(1)/\zh}$ (\cref{def:con-PhiSP})."
- Judgement: agree. By `def:fsp-wick` an identity of full symbols includes a rule for $\det(-\mathsf H)^{1/2}$. Three
  rules are chained ( `lem:fsp-cfu` for $X\to Y$; `lem:fsp-cov` and `lem:fsp-homog` in (a); `def:con-PhiSP` at $X=1$ ),
  and they agree only if $Y'(1)Y_c=(1-u)^{-1/2}$, which is not stated. It is true (two lines, below; I re-derived it:
  $W_V''(Y_c)=3Y_c^2+2V/Y_c^3=5Y_c^2$, and $\hat\cW''(1)=5/(1-u)$ is the value used in `def:con-PhiSP`).
- Change: append to the proof of (b)
  `The square roots agree: \cref{lem:fsp-cfu} takes $\hat\cW''(1)^{1/2}=Y'(1)\,W_V''(Y_c)^{1/2}$, (a) takes $W_V''(Y_c)^{1/2}=\sqrt5\,Y_c$, and \cref{def:con-PhiSP} takes $\hat\cW''(1)^{1/2}=\sqrt5\,(1-u)^{-1/2}$. Differentiating $\hat\cW=W_V\circ Y$ twice at $X=1$, where $W_V'(Y_c)=0$, gives $5/(1-u)=5Y_c^2\,Y'(1)^2$; as $Y_c$, $Y'(1)$ and $(1-u)^{-1/2}$ lie in $1+u\bbQ\fps u$, where a square root with constant term $1$ is unique, $Y'(1)\,Y_c=(1-u)^{-1/2}$.`
  (Alternative of U14_B: state (a), (b) for $\SPn$ only; then $\hat\fa(1;0)=(1-u)^{-1/2}$ and the same identity are used.)

### A3-10 [P1] `lem:fsp-dq`: the holomorphic clause is not a meaningful statement as worded, and is unused — `a_formalsp_3.tex` l.55, l.59
- Reports: U16_C M4, U1617B N2.
- Current (l.55): "For holomorphic data with a holomorphic family $c(q)$ the same identities hold for the coefficient
  functions." Current (l.59): "The holomorphic version follows from \cref{lem:fsp-param}~(ii) at every point."
- Judgement: agree. $\mathcal S_q[F;a]$ contains $\exp((F(c)-F_0(c_0))/z)$, which has no expansion in $z$ with holomorphic
  coefficients unless $F(c(q))$ is constant; the proof is one sentence. Both uses of the lemma (`lem:con-S2` over
  $\bbQ\fps{v-1}$, `lem:fsp-generic` over $\AL\fps\varepsilon$; grep of `lem:fsp-dq`) are formal.
- Change: delete both sentences. (I do not follow U1617B N2 further: the clauses for $\mathrm D$ and $P(z\mathrm D)$ are
  proved on the page and harmless; keep them.)

### A3-11 [P1] Assertions resting on computations that are not recorded — `s10_main_2.tex` l.46, l.51; `s01_intro_2.tex` l.25
- Reports: U18_A m7, U18_B F4 (README "false as written": the references to \cref{sec:evi}).
- Current (l.46): "In the normalisation of \cref{rmk:mai-stored}, computations (\cref{sec:evi}) show that six edge graphs
  of genus $4$ become singular while the singular part of the total is unchanged. So the singular part of the total
  does not determine those of the graph terms, and no argument that uses the total only through its singular part can
  prove Theorem~A." Current (l.51): "Computations (\cref{sec:evi}) give dressings that violate (A), \ldots, whose flat
  sums have singular edge graphs: neither hypothesis of \cref{thm:red-main} can be dropped in general."
  Current (intro l.25): "it is not implied by a statement about the singular part of the sum (\cref{rmk:mai-graphlocal})."
- Judgement: agree. The draft notes in both remarks say that the computations are not in the evidence appendix, so
  "(\cref{sec:evi})" points to nothing, and three sentences state as fact what the note does not record. No proof uses
  them. Either record the data (then no change), or:
- Change (l.46): the two sentences -> `In the normalisation of \cref{rmk:mai-stored}, a genus-$4$ computation, which is not recorded in \cref{sec:evi} and enters no proof, indicates that six edge graphs become singular while the singular part of the total is unchanged. If so, the singular part of the total does not determine those of the graph terms, and an argument that uses the total only through its singular part cannot prove Theorem~A.`
- Change (l.51): first sentence -> `Computations that are not recorded in \cref{sec:evi} and enter no proof indicate that dressings violating (A), among them $P^\sharp$ of \cref{rmk:mai-stored}, and dressings that are not regular in $u$ in the flat frame have flat sums with singular edge graphs; so we expect that neither hypothesis of \cref{thm:red-main} can be dropped.`
- Change (intro l.25): `it is not implied by a statement about the singular part of the sum (\cref{rmk:mai-graphlocal}).` ->
  `we know no way to deduce it from a statement about the singular part of the sum (\cref{rmk:mai-graphlocal}).`

### A3-12 [P1] Commentary after `lem:flh-inputs`: not all identities at $z^0$ are proved in $\Nring\otimes\bbQ(\zeta)$ — `s07_comparison_2.tex` l.63
- Reports: U12_B m4.
- Current: "This is why (iii) carries the factor $\Dn$ and why the identities at $z^0$ in (iv) are proved in
  $\Nring\otimes\bbQ(\zeta)$ itself."
- Judgement: agree. In the proof of (iv) (l.56-60) only the unit identity is proved in $\Nring\otimes\bbQ(\zeta)$;
  $\Psi_{\mathrm Q}=P_0\Psi^{\mathrm m}$ and the isometry are read off from (ii) at $z^0$, in $\bbQ(\zeta)\lser s$.
- Change: `and why the identities at $z^0$ in (iv) are proved in $\Nring\otimes\bbQ(\zeta)$ itself.` ->
  `and why the unit identity $\sum_j\Delta_j^{-1/2}(\Psi_{\mathrm Q})_{aj}=5\,\delta_{a3}$ in (iv) is proved in $\Nring\otimes\bbQ(\zeta)$ itself; the other identities at $z^0$ in (iv), $\Psi_{\mathrm Q}=P_0\Psi^{\mathrm m}$ and the isometry, are identities in $\bbQ(\zeta)\lser s$ obtained from (ii).`
