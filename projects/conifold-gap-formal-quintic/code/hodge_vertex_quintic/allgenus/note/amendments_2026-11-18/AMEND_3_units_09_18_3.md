# Suggested amendments, part 3: minor findings of units U09-U18 (file 3 of 3: P3 for Sections 8-10 and appendix D, R, skipped)

Conventions as in file 1 (`AMEND_3_units_09_18.md`). Paths are relative to `sections/`.

---------------------------------------------------------------------------------------------------------------------

## P3 (Sections 8, 9, 10, appendix D, introduction)

### A3-42 [P3] `lem:con-form` (a): "the critical points are $x^5=u$" is an identity of symbols — `s08_connection_3.tex` l.61
- Reports: U13_A m4, U13_B m4. Agree; (a) is used only through (c), and Theorem P(a) (ii) is proved by derivatives.
- Change: `the critical points are $x^5=u$.` ->
  `formally, the critical points are $x^5=u$ (an identity of symbols, with $\log(1-x^5)$ and $\log(1-u)$ as series; it is used below only through (c)).`

### A3-43 [P3] `def:gam-c`: the colour-$i$ series — `s08_connection_3.tex` l.23
- Reports: U13_A m6. Agree; not used in the unit.
- Change: `with the logarithm of $-q$ chosen so that $\zeta^iL$ is critical.` ->
  `with a disc around $\zeta^iL$, and with the branches and the logarithm of $-q$ chosen there as for $i=0$, so that $\zeta^iL$ is critical.`

### A3-44 [P3] `def:con-PhiSP`: $\mathrm{ph}^{\mathrm{SP}}_n\in\bbQ\fps u$ uses two more facts — `s08_connection_4.tex` l.35
- Reports: U13_B m6. Agree (both facts are on the same line of the text).
- Change: `By \cref{lem:con-form}~(c) and \cref{lem:fsp-univ}, every` ->
  `By \cref{lem:con-form}~(c) and \cref{lem:fsp-univ}, and since the covariance $(1-u)/5$ and the unit $\hat\fa(1;0)=(1-u)^{-1/2}$ lie in $\bbQ\fps u$, every`

### A3-45 [P3] `thm:con-Pa` (ii): semicolon in $\PhiSP(s^{-4}z;s^5)$ — `s08_connection_4.tex` l.52
- Reports: U13_A m5 (iv). Agree: `def:con-PhiSP` writes $\PhiSP(\zh,u)$.
- Change: `$\jmath(c)(z)=\PhiSP(s^{-4}z;s^5)$` -> `$\jmath(c)(z)=\PhiSP(s^{-4}z,s^5)$`

### A3-46 [P3] `lem:con-relab`: the letter $k$ is both the power of $\Xi/(xz)$ and the Stirling index — `s08_connection_6.tex` l.51
- Reports: U14_A m1. Agree (l.51 uses $k$ for both inside one sentence; l.52-53 use $k$ for the power).
- Change: `$\bk_k\,c_{k,P}\,(xz)^{2k-1}x^{5P}$, with $c_{k,P}:=5^{1-2k}\delta_{P,0}-5\binom{2k-2+5P}{5P}$` ->
  `$\bk_\kappa\,c_{\kappa,P}\,(xz)^{2\kappa-1}x^{5P}$ ($\kappa\geq1$ a local index), with $c_{\kappa,P}:=5^{1-2\kappa}\delta_{P,0}-5\binom{2\kappa-2+5P}{5P}$`

### A3-47 [P3] `lem:con-cfu`: "its inverse is of the same kind" — `s08_connection_5.tex` l.41
- Reports: U14_B m4. Agree; the form of the inverse needs the equivariance, two lines. (Checked: a monomial
  $u^ay^{5d}$ of the inverse becomes $u^{a+d}Y^{5d}$ after $y=sY$, which gives $\deg X_k\leq k$.)
- Change: `a change of variables in $\bbQ\fps{u,x}$ tangent to the identity; its inverse is of the same kind, and returning to $Y$ gives the stated form` ->
  `a change of variables in $\bbQ\fps{u,x}$ of the form $x\,(1+g(u,x^5))$ with $g(0,0)=0$. Its compositional inverse exists in $\bbQ\fps{u,y}$ and commutes with $y\mapsto\zeta y$, so it has the form $y\,(1+h(u,y^5))$ with $h(0,0)=0$; a monomial $u^ay^{5d}$ of $h$ is $u^{a+d}Y^{5d}$ at $y=sY$, and returning to $X$, $Y$ gives the stated form`

### A3-48 [P3] `lem:con-id` (c), `prop:con-res`: coefficient field of the gauge — `s08_connection_6.tex` l.22, `s08_connection_7.tex` l.18
- Reports: U14_B m5. Agree: `def:gau-onepole` allows any field $\Bbbk$, and (c) asserts $G_b\in\bbQ\fps u\fps\zh$.
- Change (both places): `a graded regular gauge` -> `a graded regular gauge with coefficients in $\bbQ\fps u$`

### A3-49 [P3] "\cref{thm:gam-lmi} applies at each of them" — `s09_gamma.tex` l.12
- Reports: U15_A M1, U15_B M1. Agree: true, but at $\zeta^i$ it is not what the proof does.
- Change: `and \cref{thm:gam-lmi} applies at each of them.` ->
  `and \cref{thm:gam-lmi} applies at each of them; it is used at $L=1$ and $L=0$, and the points $\zeta^i$ are reached by rotation (\cref{thm:gam-main}~(iii)).`

### A3-50 [P3] `thm:gam-lmi`, Step 1: "without its factor $z^{M_*+a(\sigma)/z}$" — `s09_gamma_2.tex` l.27
- Reports: U15_B M2. Agree: with $x=\re^\chi$ the integrand is that of
  $\int_0^\infty x^{M_*+a/z-1}\re^{-x/z}\dd x=z^{M_*+a/z}\,\Gamma(M_*+a/z)$, so the factor is present, not removed.
- Change: `motivated by the Euler integral of $\Gamma(M_*+a(\sigma)/z)$ without its factor $z^{M_*+a(\sigma)/z}$.` ->
  `motivated by Euler's integral $\int_0^\infty x^{M_*+a/z-1}\re^{-x/z}\,\dd x=z^{M_*+a/z}\,\Gamma(M_*+a/z)$, $a=a(\sigma)$, in the variable $x=\re^\chi$.`

### A3-51 [P3] `thm:gam-lmi`, Step 0: "the Taylor coefficients \ldots of the $\cF_m/\cF_0$" — `s09_gamma_2.tex` l.21
- Reports: U15_A M4. Agree: what enters is $\cF_m^{(j)}(L)/(j!\,\cF_0(L))$; `def:gam-lmidata` (`s09_gamma.tex` l.39) has
  "$\cF_0$ has no zero on $\bbD$".
- Change: `the Taylor coefficients of $\varphi$ and of the $\cF_m/\cF_0$ are holomorphic on all of $\bbD$.` ->
  `the Taylor coefficients of $\varphi$, and those of the $\cF_m$ divided by $\cF_0(L)$, are holomorphic on all of $\bbD$, where $\cF_0$ has no zero.`

### A3-52 [P3] `cor:gam-free`: what $\Phiact$ means in the conifold expansion — `s09_gamma_4.tex` l.56
- Reports: U15_A M5, U15_B M7. Agree; terse, correct, not on the main chain.
- Change: `(by \cref{thm:gam-D}, see the proof)` ->
  `(by \cref{thm:gam-D}, see the proof; here $\Phiact$ is expanded through the representatives $A^{\mathrm{act}}_n\in\bbQ(\zeta)[L,L^{-1}]$ of \cref{thm:gau-R0}~(iii))`

### A3-53 [P3] Conditionality paragraph: [DLMF, 24.4.1] is derived, not quoted — `s01_intro_3.tex` l.21
- Reports: U15_A M7. Agree: `a_formalsp_5.tex` l.69 derives it from the generating function [DLMF, 24.2.3].
- Change: `uses no quoted input, except the identity \cite[24.4.1]{DLMF} of Bernoulli polynomials in part (b) (\cref{thm:gam-ginf}).` ->
  `uses no quoted input; the identity \cite[24.4.1]{DLMF} of Bernoulli polynomials used in part (b) is derived from their generating function in the proof of \cref{lem:fsp-dlmf5118}.`

### A3-54 [P3] Draft notes left in the text of these units (editorial; no text proposed)
- Reports: U10_A m-5, U10_B m5, U15_A M6, U15_B M4, U18_A m8.
- Places: `s06_model_7.tex` l.87 (`rmk:mod-hock`); `s09_gamma_4.tex` l.72 (`rmk:gam-prior`: outline or proof in
  [Givental, Sec. 10]); `s10_main_2.tex` l.41 (`rmk:mai-constant`: numbers of [CGYZ26] read in an HTML rendering),
  l.46, l.51 (see A3-11). Each needs a decision by the author; the comparison with [CGYZ26] could not be checked here.

### A3-55 [P3] `lem:fsp-ring` (iii): "perfect matchings of the multiset" — `a_formalsp.tex` l.25
- Reports: U16_C M5. Agree: read literally $\langle t_1^4\rangle_C$ would be $C_{11}^2$, not $3C_{11}^2$; the recursion fixes it.
- Change: `the sum over the perfect matchings of the multiset $\{1^{\alpha_1},\dots,n^{\alpha_n}\}$ of` ->
  `the sum over the perfect matchings of $|\alpha|$ labelled points, $\alpha_i$ of them carrying the index $i$, of`

### A3-56 [P3] `cor:fsp-movebase`, proof: the Hessian hypothesis of `lem:fsp-ift` — `a_formalsp_3.tex` l.23
- Reports: U16_C M7. Agree: `lem:fsp-ift` (`a_formalsp_2.tex` l.96) asks for $\Hess F(0)\in GL_n(B)$.
- Change: `gives $c$: $\nabla F(c_0)=\nabla F_1(c_0)\in(q)^n$.` ->
  `gives $c$: $\nabla F(c_0)=\nabla F_1(c_0)\in(q)^n$, and $\Hess F(c_0)\in\mathsf H_0+q\,\mathfrak{gl}_n(A)$ is invertible.`

### A3-57 [P3] `lem:fsp-fubini` (ii): $w_0$ is non-degenerate for $F^\flat$ — `a_formalsp_4.tex` l.56
- Reports: U17_C F4 (second half; the first half is in E2). Agree: (iii) takes $\SPn_{w_0}[F^\flat;\cdot]$.
- Change: `and $\det\mathsf H=\det\mathsf H_{vv}\det\mathsf S$.` ->
  `and $\det\mathsf H=\det\mathsf H_{vv}\det\mathsf S$; hence $\det\mathsf S\in A^\times$, and $w_0$ is a non-degenerate critical point of $F^\flat$.`

### A3-58 [P3] `lem:fsp-cov`: typography of the formal map — `a_formalsp_4.tex` l.19
- Reports: U17_C F5 (last point). Agree.
- Change: `$\phi\in(c+(y-c')A\fps{y-c'})^n$` -> `$\phi\in c+\bigl((y-c')A\fps{y-c'}\bigr)^n$`

### A3-59 [P3] `cor:mai-gap`, Step 2: the reason for $E_{R'}=E_X$ — `s10_main_2.tex` l.27
- Reports: U18_A m3. Agree: `lem:flh-frame` is about a push-forward; the direct reason is one line and needs no
  symplecticity. (The parenthesis on [PPZ15] in the same paragraph is rewritten by E6.)
- Change: `so $E_{R'}=E_X$ and $q_{R'}=q_X$, by the computation of the proof of \cref{lem:flh-frame} for the isometry $A=\Psi_{\mathrm{cf}}$.` ->
  `so $E_{R'}=E_X$ and $q_{R'}=q_X$: $R'^{-1}=X^{-1}\Psi_{\mathrm{cf}}^{-1}$, and $\Psi_{\mathrm{cf}}^{-1}\otimes\Psi_{\mathrm{cf}}^{-1}$ maps the inverse bivector of $(\bbH,\eta)$ to that of the colour space, $\Psi_{\mathrm{cf}}$ being an isometry.`

### A3-60 [P3] Introduction, Theorem C: the frame is not quantified in the first sentence — `s01_intro_2.tex` l.46
- Reports: U18_A m4. Agree: $\Dn_0$ is defined through $\Pc$; the body (`thm:gam-D`) says "In every conifold frame".
- Change: `$\Dn_0=1$ in $\bbQ\fps{z}$. Hence` -> `In every conifold frame, $\Dn_0=1$ in $\bbQ\fps{z}$. Hence`
  and, in the next sentence, `Equivalently, in every conifold frame $\Qr=` -> `Equivalently, $\Qr=`

### A3-61 [P3] `cor:mai-natfact` (c): "$G$ depends on the frame" is asserted without proof — `s10_main_3.tex` l.35
- Reports: U18_A m6 (last sentence). Agree; unused.
- Change: `$G$ depends on the frame.` -> `$G$ is built from the canonical gauge of the frame and may depend on it.`

### A3-62 [P3] `cor:mai-univ`: the hypothesis on $r$ is automatic; the field is not stated — `s10_main_4.tex` l.21, l.28
- Reports: U18_A m10 (first point), U18_B F5 (second point). Agree, for a system over $\bbQ\fps u$: self-duality
  (`def:gau-onepole`: $b_0=b_2$, $b_3=b_4$) gives $b_1=\Bprod/(b_0^2b_3^2)$, and $\Res_{u=0}(-v_A'^5/v_A)=-v_A'(0)^4$, so
  $r=-v_A'(0)^4/(b_0(0)b_3(0))^2$; the paragraph at l.28 says the same. The statement does not say over which field
  $A$ is defined (`def:gau-onepole` allows any $\Bbbk$), and $\lambda\in\bbQ^\times$ needs $\bbQ$.
- Change (l.21): `let $A$ be a self-dual one-pole system (\cref{def:gau-onepole}) with product $-v_A'^5/v_A$ and residue $r$; assume $r=-1/\lambda^2$ with $\lambda\in\bbQ^\times$ (a letter local to this corollary and its proof).` ->
  `let $A$ be a self-dual one-pole system over $\bbQ\fps u$ (\cref{def:gau-onepole}) with product $-v_A'^5/v_A$ and residue $r$. Then $r=-1/\lambda^2$ with $\lambda\coloneqq b_0(0)b_3(0)/v_A'(0)^2\in\bbQ^\times$ (a letter local to this corollary and its proof), because $b_1=\Bprod/(b_0^2b_3^2)$.`
- Change (l.28): the whole sentence `For a system $b=(b_0,\Bprod/(b_0^2b_3^2),b_0,b_3,b_3)$ with $\Bprod=-v_A'^5/v_A$ one has \ldots; the vertex law is even in $\lambda$.` ->
  `Replacing $\lambda$ by $-\lambda$ changes nothing: the vertex law is even in $\lambda$.`

### A3-63 [P3] `rmk:mai-stored`: "(G) fixes the odd part of $\log\cn$" — `s10_main_4.tex` l.41
- Reports: U18_B F6. Agree, mildly: by the sentence before, (A) forces $\cn=1$, not only the odd part.
- Change: `(G) fixes the odd part of $\log\cn$.` ->
  `(G) fixes $\cn$ altogether, in particular the odd part of $\log\cn$, which symplecticity leaves free.`
  (The same phrase in `rmk:mai-graphlocal`, `s10_main_2.tex` l.46, "the odd part of $\log\cn$, which is fixed by (G)
  \ldots and not by symplecticity", is accurate and can stay.)

---------------------------------------------------------------------------------------------------------------------

## R (no change proposed)

| Finding | Place | Reason |
|---|---|---|
| U10_A m-6, U10_B m4 | Section 6: $\mathcal B_n$ / `\Bprod`, $\mathcal R^b_n$ / `\cR`, $I(F)$ / $I_m$, $B$, $w$ | Local letters, each declared where used; renaming would touch many lines for no gain. |
| U11_A F6 | `s06_model_9.tex` l.45-47: $k=\varsigma$, $\rho_b$ / $\varrho_b$, $r$, $t$, $\lambda$ | Same; the one real clash of a letter inside a sentence in these units is A3-46. |
| U11_B m4 (first half) | `lem:mod-ibp`: hypothesis "zero residues" unused | True but harmless, and it holds for the forms the lemma is applied to. |
| U12_A F6, U12_B m7 | Section 7: $E$ (edge / number of edges), $N$, $\cA$ | Local; the text says "a local letter" at the definitions. |
| U12_A F9 | "the flat sum of $\Pc.\Omega^{\mathrm{mod}}$" | The definition (`def:red-flatsum`) is cited at each use. |
| U13_A m5 (i)-(iii), U13_B m5 | Section 8: $u\,\dd/\dd u$ at fixed $z$, $c_0(L)$ / $c_0\in\Nring$, $\lambda$ | Declared local in the text; (iv) is A3-45. |
| U14_B m6 | `thm:con-glue` (a), `s08_connection_8.tex` l.19: "so the constant series is $\Dn_0$" | Redundant with the definition of $\Dn_0$ in `thm:gau-R0`, but correct, and it is where $\Dn_0\in\bbQ(\zeta)\fps z$ is tied to $A^{\mathrm{act}}_n\in\bbQ(\zeta)[L,L^{-1}]$. |
| U14_B notation notes | $\mathsf s$ next to $s$; two meanings of $G_b$ | Both declared local. |
| U15_A M3, U15_B M3 | Section 9: $\delta$, $\lambda$, $\mathsf Q$, $\phi$, $c$ | Declared local. |
| U17_C F5 (letters) | Appendix D: $\beta$, $\rho$, $B$, $\psi$, $\chi$, $J$ | Local to single lemmas; the typographical point is A3-58. |
| U18_A m6 (letters) | Section 10: $\tilde R$, $G$, $T$, $\Psi_{\mathrm Q}$ | Acknowledged in the text (`s10_main.tex` l.35); see A3-40 for the one optional clarification. |
| U1617B N2 | `lem:fsp-dq` ($\mathrm D$, $P(z\mathrm D)$ clauses), `lem:fsp-param` (iii) unused | Proved on the page; deleting them is a matter of taste. Only the holomorphic clause is removed (A3-10). |
| U18_A m8 | `rmk:mai-constant`: numbers of [CGYZ26] | Already flagged by the draft note in the text; nothing to amend until the PDF is read (A3-54). |

---------------------------------------------------------------------------------------------------------------------

## Skipped: ruled on elsewhere

| Finding | Subject | See |
|---|---|---|
| U09_B M1 | [DNOPS, (2-33)] misprinted; Step 1 of `thm:mod-tr` | E5 (README must-fix 2) |
| U09_B M6 | [DNOPS, Prop. 4.3] (1) against (2); the sentence before Thm. 8.1 | E5 |
| U09_B M9 | dilaton sign "not a convention" (`s06_model_2.tex` l.55) | E5 |
| U12_B m3 | hatted unit $\hat1=s1$ in `lem:flh-hat` | E8 |
| U13_A m1, U13_B m1 | `eq:con-bz` proved "factor by factor" | E1 (README routine) |
| U13_A m2, U13_B m3 | $\mathsf S_0(L)$, $A_{\mathrm{St}}(L,0)\notin\AL$ in `def:gam-c` and the proof of `thm:con-Pa` (i) | E1 |
| U14_A m3, U14_B m2 | ring of the second homomorphism, `thm:con-Pa` (ii) Steps 3-4 | E1 |
| U14_B m3 | `lem:fsp-cov` applied to $\mathsf S_0$ over $\bbQ\lser s$ | E1 |
| U14_A m4, U17_C F3, U1617B N3 | `lem:fsp-res` (iii), exponential under $X\mapsto X+\chi$ | E3 (README routine) |
| U16_C M1, M2, M3, M6 | $E_{zC}$ with $q$-dependent $C$; `lem:fsp-param` (i), (iv); linearity of $E^t_{zC}$ | E2 (README routine) |
| U17_C F1, F2, F4 (first half) | `lem:fsp-cov`/`cfu` "every amplitude"; $P_\pm$ over $\Bbbk$; $E^v$ and $\re^{F^\flat_{\geq3}/z}$ | E2 |
| U1617B N1 | `lem:fsp-generic` asks for a phase over $\AL$; $\mathsf S_0(L)\notin\AL$ | E1 (closed by replacing $\mathsf S_0$ by $\mathsf S_0-\mathsf S_0(L)$; no change to `a_formalsp_3.tex` l.64 is then needed) |
| U18_A m1, U18_B F1 | Galois descent in `def:mai-Fcf` cites `lem:gau-colrows` | README routine (X4, U06) |
| U18_A m2, U18_B F2 | what Corollary B is a statement about; injectivity of $\ev$ | E7 (README must-fix 3) |
| U18_B F3 | "$\TFT$ the TFT of `def:set-frob`" without "expanded in the frame" | E8 |
| U18_A m9 (part) | which propositions of [PPZ15] enter, over which ring | E6; the residual mismatch of the two lists is A3-24 |

## Counts

P1: 12 (A3-01..12). P2: 13 (A3-13..25). P3: 38 (A3-26..63). R: 13 rows. Skipped (see E# / README): 16 rows.
