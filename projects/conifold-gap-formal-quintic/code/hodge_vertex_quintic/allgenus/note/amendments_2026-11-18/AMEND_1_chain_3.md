# Suggested amendments 1: the proof chain. Part 3 of 4 (B1-B5)

Conventions as in `AMEND_1_chain.md`.

## B. Routine missing steps

### B1. `lem:fsp-res`: the exponential under X -> X + chi; coefficients of P_+-  [E3; E2 (U17-4); U14 m4, U17 F3]

Confirmed in the text. `a_formalsp_5.tex` l.94, end of (iii): "The map $\Omega\mapsto\Omega(X+\chi)$ is a ring
homomorphism, equal on Laurent polynomials to the substitution ...; so for $\Omega=\re^{-f(Y)/\zh}\beta(Y)$ it is
$\re^{-f(X+\chi)/\zh}\beta(X+\chi)$, expanded termwise in $u$." Proof, l.98, last sentence: "The map is multiplicative by
the Leibniz rule for higher derivatives, and it is the binomial expansion on $X^{\pm1}$." Nothing determines the image of
$\re^{-f_0/\zh}$, which is not a Laurent polynomial and not u-adically small. l.90 (i): "$P_\pm$ polynomials without
constant term" (the bound in the proof needs coefficients in the field, not in its Laurent series).

Checked. The statement is true. E3 closes it with complete Bell polynomials (Section 3) or with an ODE in a parameter
(alternative (a)); I re-derived the second, which is shorter, and use it below: d/dt T_t(Omega) = chi T_t(Omega'),
uniqueness of the solution of dZ/dt = chi T_t(g') Z, and the substitution t = 1. The only use is `prop:con-res`
(f_0 = 1/Y + Y^4/4 over Q), which needs exactly this statement (E3 Section 4).

**B1.1** `a_formalsp_5.tex` l.90. Old: `$P_\pm$ polynomials without constant term;` New: `$P_\pm\in X\Bbbk[X]$ (polynomials over $\Bbbk$ without constant term);`

**B1.2** `a_formalsp_5.tex` l.94. Old: `so for $\Omega=\re^{-f(Y)/\zh}\beta(Y)$ it is $\re^{-f(X+\chi)/\zh}\beta(X+\chi)$, expanded termwise in $u$.`
New: `and for $\Omega=\re^{-f/\zh}\beta$ as in (ii), $\beta\in\Bbbk[\zh,\zh^{-1}][X,X^{-1}]\fps u$, it is $\re^{-f(X+\chi)/\zh}\beta(X+\chi)$, where $f(X+\chi)$ and $\beta(X+\chi)$ are the substitutions and the exponential is read with the convention of (ii): $f(X+\chi)=f_0(X)+\bigl(f(X+\chi)-f_0(X)\bigr)$, the second summand lying in $u\Bbbk[X,X^{-1}]\fps u$.`

**B1.3** `a_formalsp_5.tex` l.98 (proof). In (ii), old: `and $\frac{\dd}{\dd X}(\re^{-f/\zh}Q)=\re^{-f/\zh}(Q'-\zh^{-1}f'Q)$.`
New: `and $\frac{\dd}{\dd X}(\re^{-f/\zh}Q)=\re^{-f/\zh}(Q'-\zh^{-1}f'Q)$, because $E:=\re^{-f_0/\zh}$ satisfies $E'=-\zh^{-1}f_0'E$: for each of the two factors $\sum_kp^k/k!$ of (i), $p=-P_\pm/\zh$, differentiate termwise, the coefficients being finite sums.`
In (iii), old: `The map is multiplicative by the Leibniz rule for higher derivatives, and it is the binomial expansion on $X^{\pm1}$.` New:

```latex
The map $\mathsf T\colon\Omega\mapsto\Omega(X+\chi)$ is multiplicative by the Leibniz rule for higher derivatives, linear over $\Bbbk\lser{\zh^{-1}}\fps u$ and $u$-adically continuous, and it is the binomial expansion on $X^{\pm1}$; so it is the substitution $X\mapsto X+\chi$ on $\Bbbk\lser{\zh^{-1}}[X,X^{-1}]\fps u$, and $\mathsf T(\re^{-f_1/\zh})=\re^{-f_1(X+\chi)/\zh}$. It remains to find $\mathsf T(E)$ for $E=\re^{-f_0/\zh}$, which is neither a Laurent polynomial nor $u$-adically small. Put $g:=-f_0/\zh$, so that $E'=g'E$, and for a variable $\mathsf t$ let $\mathsf T_{\mathsf t}(\Omega):=\sum_j\Omega^{(j)}\mathsf t^j\chi^j/j!\in\cR_{\mathrm{res}}\fps u\fps{\mathsf t}$ (the letters $\mathsf T$, $\mathsf t$, $Z$ are local). It is multiplicative, and $\de_{\mathsf t}\mathsf T_{\mathsf t}(\Omega)=\chi\,\mathsf T_{\mathsf t}(\Omega')$. Hence $Z=\mathsf T_{\mathsf t}(E)$ satisfies $\de_{\mathsf t}Z=\chi\,\mathsf T_{\mathsf t}(g')\,Z$ and $Z|_{\mathsf t=0}=E$. So does $Z=E\exp\bigl(\mathsf T_{\mathsf t}(g)-g\bigr)$, where $\mathsf T_{\mathsf t}(g)-g=\sum_{j\geq1}g^{(j)}\mathsf t^j\chi^j/j!$ has no constant term in $\mathsf t$, so that the exponential is defined. A solution $Z=\sum_nZ_n\mathsf t^n$ is determined by $Z_0$, since $(n+1)Z_{n+1}=\sum_{i\leq n}\coef{\mathsf t^i}\bigl(\chi\,\mathsf T_{\mathsf t}(g')\bigr)Z_{n-i}$. Therefore $\mathsf T_{\mathsf t}(E)=E\exp(\mathsf T_{\mathsf t}(g)-g)$. On both sides the coefficient of $\mathsf t^n$ lies in $u^n\cR_{\mathrm{res}}\fps u$, so $\mathsf t=1$ may be substituted: $\mathsf T(E)=\re^{-f_0(X)/\zh}\,\re^{-(f_0(X+\chi)-f_0(X))/\zh}$. By multiplicativity, $\Omega(X+\chi)=\re^{-f_0(X)/\zh}\bigl(\re^{-(f(X+\chi)-f_0(X))/\zh}\beta(X+\chi)\bigr)$ for $\Omega=\re^{-f/\zh}\beta$, which is $\re^{-f(X+\chi)/\zh}\beta(X+\chi)$ in the convention of (ii).
```

---------------------------------------------------------------------------------------------------------------------

### B2. `eq:con-bz` (Theorem P(a) (i), Step 3): the factors are not in the ring named  [E1; U13-1/10]

Confirmed in the text. `s08_connection_4.tex` l.79: "which embeds in $\AL(\zeta)\fps{t}\fps{z}$ ...; we prove
\eqref{eq:con-bz} there, factor by factor ... Write the Stirling-expanded second summand as the product of the six
factors $\re^{\sigma\log q/z}$, $\Gst(1+5\sigma/z)$ and $\Gst(1+(\sigma-\zeta^j)/z)^{-1}$". None of these factors, nor
the multipliers $\prod_k(5\sigma/z+k)$, lies in that ring (logarithms, negative powers of z). Also l.61: "$\mathsf S_0$ is
the formal phase at $\sigma=L$ given by its Taylor data, which lie in $\AL$" (the constant $\mathsf S_0(L)$ and
$A_{\mathrm{St}}(L,0)$ do not).

Departure from E1. E1's repair (Section 3, R1-R5) works in a polynomial ring over k[[t]]((z)) in seven logarithm
symbols and proves an intertwining of two ring homomorphisms. A shorter proof stays inside power series: per factor
rho, the three honest series U_rho (Taylor part of rho log rho - rho), l_rho = log(rho/rho_0) and the Stirling tail
satisfy one identity (dagger below), which is the statement St_m = Poch_m St_0 of `lem:fsp-stirfe` at w = z/rho. It works
with the first Mellin-Barnes form (the definition of c in `def:gam-c`), so the passage to the second summand, and with
it `lem:con-form` (d) and the branch argument, is no longer used in this proof. Both repairs prove the same identity.

Checked by hand: (dagger); the decomposition of S_0, S_1, S_2k of `def:gam-c` into the six factors (the linear terms of
S_0 cancel by the definition of lambda); the product of the multipliers, with prod_j (zeta^j - x) = 1 - x^5 and
prod_j (zeta^j - L)/(5L)^5 = -q. E1's own symbolic checks (work/E1) confirm the end result to t^4 z^6.

**B2.1** `s08_connection_4.tex` l.61. Old: `In this proof $\mathsf S_0$ is the formal phase at $\sigma=L$ given by its
Taylor data, which lie in $\AL$ (\cref{def:gam-c}), and $\sigma=L+t$. For $\chi\in\bbQ[\sigma,z]$ put
$c\langle\chi\rangle(z;L):=\SPn_L[\mathsf S_0;\chi A_{\mathrm{St}}]/A_{\mathrm{St}}(L,0)$, so that $c=c\langle1\rangle$;`
New: `In this proof $\sigma=L+t$, $\mathsf S_0$ stands for the formal phase $\mathsf S_0(L+t)-\mathsf S_0(L)\in t^2\AL\fps t$
at $\sigma=L$ (\cref{def:gam-c}; the constant does not matter, \cref{lem:fsp-const}), and
$a_0:=A_{\mathrm{St}}/A_{\mathrm{St}}(L,0)=\exp\bigl(\mathsf S_1(L+t)-\mathsf S_1(L)+\sum_kz^{2k-1}\mathsf S_{2k}(L+t)\bigr)\in\AL\fps t\fps z$
for the normalised amplitude, both given by the Taylor data of \cref{def:gam-c}. For $\chi\in\bbQ[\sigma,z]$ put
$c\langle\chi\rangle(z;L):=\SPn_L[\mathsf S_0;\chi a_0]$, so that $c=c\langle1\rangle$;`
and at the end of the same line, old: `The letters $\chi$, $a$, $b$, $b_z$ are local to this proof.` new: `The letters
$\chi$, $a$, $a_0$, $b$, $b_z$, $\mathsf E$, $\mathsf T_z$, $\rho$, $\mathsf l_\rho$, $\mathsf U_\rho$, $\mathsf m_\rho$, $w$ are local to this proof.`

**B2.2** `s08_connection_4.tex` l.74-79: replace Step 3 (from `\emph{Step 3 (the shift).}` to the end of the proof,
keeping `\end{proof}`) by:

```latex
\emph{Step 3 (the shift).} Put $a:=(\sigma^5-1)a_0$, and let $\mathsf T_z$ be the continuous ring endomorphism $f(t,z)\mapsto f(t+z,z)$ of $\AL\fps{t,z}$. The translate of $a$ by $z$ (\cref{lem:fsp-transl}~(b)) is $b_z=\re^{\mathsf E}\,\mathsf T_z(a)$ with
\[
\mathsf E:=\bigl(\mathsf S_0(L+t+z)-\mathsf S_0(L+t)\bigr)/z=\sum_{n\geq1}\mathsf S_0^{(n)}(L+t)\,z^{n-1}/n!\in\AL\fps{t,z},
\]
which has no constant term, because $\mathsf S_0'(L)=0$. We claim
\begin{equation}
\label{eq:con-bz}
b_z=q\prod_{k=1}^5(5\sigma+kz)\,a_0\qquad\text{in }\AL\fps{t,z},\quad\sigma=L+t .
\end{equation}
The ring $\AL\fps{t,z}$ embeds in $\Bbbk\fps{t,z}$, $\Bbbk:=\AL\otimes\bbQ(\zeta)$, in which $5L$ and the $\zeta^j-L$ are units; we prove \eqref{eq:con-bz} there. Let $\rho(\sigma)=\alpha\sigma+\beta$ with $\alpha\in\bbZ$, $\beta\in\bbQ(\zeta)$ and $\rho_0:=\rho(L)\in\Bbbk^\times$, so that $\rho(L+t)=\rho_0+\alpha t$, and put
\[
\mathsf l_\rho:=\log(1+\alpha t/\rho_0)\in t\Bbbk\fps t,\qquad\mathsf U_\rho:=\rho\,\mathsf l_\rho-\alpha t\in t^2\Bbbk\fps t,\qquad\mathsf m_\rho:=\sum_{k\geq1}\bk_k\,z^{2k-1}\rho^{1-2k},\qquad w:=z/\rho\in z\Bbbk\fps{t,z};
\]
these are the Taylor parts at $\sigma=L$ of $\log\rho$ (orders $\geq1$) and of $\rho\log\rho-\rho$ (orders $\geq2$), and the Stirling tail $\mu(\rho/z)$ of \cref{def:fsp-stirling}. Then
\begin{equation}
\label{eq:con-factor}
\frac{\mathsf T_z\mathsf U_\rho-\mathsf U_\rho}z+(\mathsf T_z-1)\bigl(-\tfrac12\mathsf l_\rho+\mathsf m_\rho\bigr)=\alpha\,\mathsf l_\rho+\log\mathrm{Poch}_\alpha(z/\rho)\qquad\text{in }\Bbbk\fps{t,z}.
\end{equation}
Indeed $\mathsf T_z\rho=\rho\,(1+\alpha w)$, so $\mathsf T_z\mathsf l_\rho=\mathsf l_\rho+\log(1+\alpha w)$ and $(\mathsf T_z\mathsf U_\rho-\mathsf U_\rho)/z=\alpha\,\mathsf l_\rho+(w^{-1}+\alpha)\log(1+\alpha w)-\alpha$, while $\mathsf m_\rho=\sum_k\bk_kw^{2k-1}$ and $\mathsf T_z\mathsf m_\rho=\sum_k\bk_kw^{2k-1}(1+\alpha w)^{1-2k}$. So the left side of \eqref{eq:con-factor} is $\alpha\,\mathsf l_\rho+\log\mathrm{St}_\alpha(w)-\log\mathrm{St}_0(w)$, and $\mathrm{St}_\alpha=\mathrm{Poch}_\alpha\,\mathrm{St}_0$ in $\bbQ\fps w$ by \cref{lem:fsp-stirfe}; the substitution $w=z/\rho$ is a continuous ring homomorphism $\bbQ\fps w\to\Bbbk\fps{t,z}$.

Apply \eqref{eq:con-factor} to $\rho_{\mathrm N}:=5\sigma$ ($\alpha=5$) and to $\rho_j:=\zeta^j-\sigma$ ($\alpha=-1$), $j=0,\dots,4$. By \cref{def:gam-c}, $\mathsf S_0(L+t)-\mathsf S_0(L)=\mathsf U_{\rho_{\mathrm N}}+\sum_j\mathsf U_{\rho_j}$ (the terms linear in $t$ cancel by the definition of $\lambda$, that is, $\mathsf S_0'(L)=0$), $\mathsf S_1(L+t)-\mathsf S_1(L)=\frac12\mathsf l_{\rho_{\mathrm N}}-\frac12\sum_j\mathsf l_{\rho_j}$ and $\sum_kz^{2k-1}\mathsf S_{2k}(L+t)=\mathsf m_{\rho_{\mathrm N}}+\sum_j\mathsf m_{\rho_j}$. Summing \eqref{eq:con-factor} over the six factors and adding $(\mathsf T_z-1)\mathsf l_{\rho_{\mathrm N}}=\log(1+5z/\rho_{\mathrm N})$, with $\mathrm{Poch}_5(w)(1+5w)=\prod_{k=1}^5(1+kw)$ and $\mathrm{Poch}_{-1}(w)=(1-w)^{-1}$,
\[
\mathsf E+(\mathsf T_z-1)\log a_0=5\,\mathsf l_{\rho_{\mathrm N}}-\sum_j\mathsf l_{\rho_j}+\sum_{k=1}^5\log\Bigl(1+\frac{kz}{\rho_{\mathrm N}}\Bigr)-\sum_j\log\Bigl(1-\frac z{\rho_j}\Bigr).
\]
Exponentiate, using $\re^{\mathsf l_\rho}=\rho/\rho_0$, $\prod_j(\zeta^j-x)=1-x^5$ and $\prod_j(\zeta^j-L)/(5L)^5=-q$ (\cref{def:gam-c}):
\[
\re^{\mathsf E}\,\frac{\mathsf T_z(a_0)}{a_0}=\frac{\prod_{k=1}^5(5\sigma+kz)}{(5L)^5}\cdot\frac{\prod_j(\zeta^j-L)}{\prod_j(\zeta^j-\sigma-z)}=\frac{q\prod_{k=1}^5(5\sigma+kz)}{(\sigma+z)^5-1}.
\]
Multiplying by $\mathsf T_z(\sigma^5-1)\,a_0=((\sigma+z)^5-1)\,a_0$ gives \eqref{eq:con-bz}. By \cref{lem:fsp-transl}~(b), with $p=z$ and $a\in\AL\fps t\fps z$, $\SPn_L[\mathsf S_0;b_z]=\SPn_L[\mathsf S_0;a]$, that is, $c\langle q\prod_k(5\sigma+kz)\rangle=c\langle\sigma^5-1\rangle$; by Step~2, $\PFfive^{L}(L\,c)=0$.
```

The label `eq:con-bz` is kept (its display changes from `A_{\mathrm{St}}/A_{\mathrm{St}}(L,0)` to `a_0`, the same
element); `eq:con-factor` is new. In Step 1 (l.63) the two occurrences of `\chi A_{\mathrm{St}}/A_{\mathrm{St}}(L,0)` may
be written `\chi a_0`; nothing else changes there. `lem:con-form` (d) stays as a statement (it is still what identifies
the two Mellin-Barnes forms) but is no longer cited in this proof.

---------------------------------------------------------------------------------------------------------------------

### B3. Appendix D: scope of the Gaussian functional and of four lemmas  [E2; U16-2/3/5, U17-1/2]

Confirmed in the text: `a_formalsp.tex` l.16 ("coefficientwise in $q$", with no definition of E_{zC} for a
q-dependent C on $\cR_n(A_0)\fps q$); `a_formalsp_3.tex` l.36 ("lie in a subring $A'\subset A$"), l.42 ((iv) names neither
the critical point nor the class of a); `a_formalsp_4.tex` l.19 ("for every amplitude $a$"), l.39 ("for every amplitude
$\beta$ at $Y_0$ over $A$"), l.77 ("the linearity of $E^v$ over functions of $w$ and $z$").

Checked. The counterexample to `lem:fsp-param` (i) over Z (coefficient 15/2 for F = -t^2/2 + t^3) was recomputed.
The definition below (E2's "U16 form") agrees with the old one on $\cR_n(A)$. E2's list of downstream uses (all inside
the restricted scope) was read, not re-derived line by line.

**B3.1** `a_formalsp.tex` l.16 (`def:fsp-ring`, "Parameters"). Old: `Over $A=A_0\fps{q}$ we work in
$\cR_n(A_0)\fps{q}\supset\cR_n(A)$, coefficientwise in $q$.` New:

```latex
Over $A=A_0\fps{q}$ we work in $\cR_n(A_0)\fps{q}\supset\cR_n(A)$, coefficientwise in $q$, in the following sense. For a symmetric matrix $C$ over $A$, which may depend on $q$, and $P=\sum_dq^dP_d$ with $P_d\in\cR_n(A_0)$, put $E_{zC}[P]:=\sum_dq^dE_{zC}[P_d]\in A_0\lser z\fps q$, where $E_{zC}[P_d]\in A\lser z$; on $\cR_n(A)$ this is the functional defined above. Modulo $q^N$ it is the Gaussian functional over $B_N:=A_0[q]/(q^N)$ on $\cR_n(B_N)=\cR_n(A_0)\fps q/(q^N)$. Identities in $\cR_n(A_0)\fps q$, in particular between exponentials, are meant modulo every $q^N$, in $\cR_n(B_N)$, where an element divisible by $q$ is nilpotent and $\re^{Y+X}=\re^Y\re^X$ for $Y$ of weight $\geq1/2$ and $X$ nilpotent. Likewise, for an amplitude $a\in A_0\fps{x-c_0}\lser z\fps q$ and a critical point $c\in c_0+(qA)^n$, the series $\SPn_c[F;a]\in A_0\lser z\fps q$ of \cref{def:fsp-wick} is defined modulo every $q^N$, where $a(c+t)$ is an amplitude over $B_N$.
```

**B3.2** `a_formalsp_2.tex` l.69 (`lem:fsp-pg` (iii)) and l.78 (its proof). l.69, old: `coefficientwise in $q$.` new:
`coefficientwise in $q$ (in the sense of \cref{def:fsp-ring}).` l.78, old: `Over $A_0\fps{q}$ the same holds at each power
of $q$, to which finitely many $k$ contribute.` new: `Over $A_0\fps q$ argue modulo $q^N$: over $B_N=A_0[q]/(q^N)$ the
ideal $(q)$ is nilpotent, (i) applies with $R=B_N$ and $I=(q)$, and the monomial argument above gives the identity on
$\cR_n(B_N)$.`

**B3.3** `a_formalsp_3.tex` l.23 (proof of `cor:fsp-movebase`). Old: `in $\cR_n(A_0)\fps{q}$` (before the display)
new: `in $\cR_n(A_0)\fps{q}$, that is, modulo every $q^N$ (\cref{def:fsp-ring}; both sides are the exponential of
$(F(c+t)-F_0(c_0)-\frac12t^{\mathrm T}\mathsf H_0t)/z$),`

**B3.4** `a_formalsp_3.tex` l.36 (`lem:fsp-param` (i)). Old: `lie in a subring $A'\subset A$` New: `lie in a $\bbQ$-subalgebra $A'\subset A$`

**B3.5** `a_formalsp_3.tex` l.42 (`lem:fsp-param` (iv)). Old: `Let $g\in A\fps{x-c}$, $g(c)=0$, $\varepsilon$ a
variable` New: `Let $c$ be a non-degenerate critical point of $F$ over $A$, $a\in A\fps{x-c}\lser z$ an amplitude (free of
$\varepsilon$), $g\in A\fps{x-c}$, $g(c)=0$, $\varepsilon$ a variable`

**B3.6** `a_formalsp_4.tex` l.19 (`lem:fsp-cov`). Old: `and for every amplitude $a$` New: `and for every amplitude $a\in A\fps{x-c}\lser z$`
`a_formalsp_4.tex` l.39 (`lem:fsp-cfu`). Old: `for every amplitude $\beta$ at $Y_0$ over $A$` New: `for every amplitude $\beta\in A\fps{Y-Y_0}\lser{\zh}$`
`s08_connection_5.tex` l.8 (the use in Theorem P(a) (ii), Step 2). Old: `for every amplitude $a$.` New: `for every amplitude $a\in A\fps{\sigma-L}\lser z$.`
(All uses have amplitudes of this kind; E2 lists them.)

**B3.7** `a_formalsp.tex` l.27 (`lem:fsp-ring` (iv)) and l.33 (proof). l.27, old: `and $E_{zC}[fP]=fE_{zC}[P]$ for $f\in A\lser{z}$.`
new: `$E_{zC}[fP]=fE_{zC}[P]$ for $f\in A\lser{z}$, and $E^t_{zC}[gP]=g\,E^t_{zC}[P]$ for $g\in\cR_m(A)$, a series in the outer variables and $z$, and $P\in\cR_{n+m}(A)$.`
l.33, add at the end of the proof of (iv): `For $g\in\cR_m(A)$ let $g_\gamma$ be its monomials; the family $(g_\gamma P)_\gamma$ is summable with sum $gP$, and $E^t_{zC}[g_\gamma P]=g_\gamma E^t_{zC}[P]$ by the definition on monomials.`
`a_formalsp_4.tex` l.77 (proof of `lem:fsp-fubini`). Old: `the linearity of $E^v$ over functions of $w$ and $z$,` New: `the linearity of $E^v$ over $\cR_{n_2}(A)$, the series in $w$ and $z$ (\cref{lem:fsp-ring}~(iv)),`
(`n_2` is the number of the variables $w$ in `lem:fsp-fubini`.)

---------------------------------------------------------------------------------------------------------------------

### B4. Propositions of [PPZ15] over rings other than Q  [E6; X1-F4, X5-F1, L1-M1]

Confirmed in the text. `s02_setup_4.tex` l.79, last sentence of `def:set-graphsum`: "The definitions and propositions
of \cite{PPZ15} used in the note are identities between finite sums, polynomial in the coefficients; stated there over
$\bbQ$, they hold verbatim over every commutative $\bbQ$-algebra." `s10_main_2.tex` l.27: "the proofs of their
Props.~2.4, 2.12 and 2.14, and of Props.~2.8 and~2.9, ... hold verbatim over $\bbQ(\zeta)\lser{s}$". As E6 says,
polynomiality alone is not a reason (the data are constrained), and the list of proofs is incomplete.

Checked. E6's argument (B) was re-derived: over a Q-algebra the series R = 1 + O(z) with R(-z)^T R(z) = 1 are exactly
exp(r) with r_k^T = (-1)^{k+1} r_k, so the data (y, r) are free; rational points with y in (Q^x)^r are Zariski dense; at
them [PPZ15] applies as stated. With it only the *statements* of Props. 2.12 and 2.14 over Q are used, so Props. 2.4,
2.8, 2.9 leave the list of inputs (simpler than E6's alternative of adding 2.3, 2.7, Lemma 2.11 and Section 0.3).
I did not re-read [PPZ15]; that Props. 2.12 and 2.14 are the only propositions used is E6's count of the 70 citations.

**B4.1** `s02_setup_4.tex` l.79. Old: `The definitions and propositions of \cite{PPZ15} used in the note are identities
between finite sums, polynomial in the coefficients; stated there over $\bbQ$, they hold verbatim over every commutative
$\bbQ$-algebra.` New:

```latex
The definitions of \cite{PPZ15} used in the note are finite sums; they make sense over every commutative $\bbQ$-algebra $A$, for $V$ free of finite rank with a perfect pairing and classes in $H^*(\overline{\mathcal M}_{g,n},\bbQ)\otimes_\bbQ A$. Two propositions of \cite{PPZ15}, stated there over $\bbQ$, are used over other rings: \cite[Prop.~2.12]{PPZ15} ($R.\Omega$ is a cohomological field theory with the unit of $\Omega$) and \cite[Prop.~2.14]{PPZ15} (the unit-preserving action is a left group action). They are applied only to a topological field theory $\TFT_y$ on $A^r$ with an orthonormal basis $e_i$, $(\TFT_y)_{g,n}(e_i,\dots,e_i)=y_i^{-(2g-2+n)}$, zero on mixed colours, with the unit $\sum_iy_ie_i$, $y_i\in A^\times$, and to series $R=1+O(z)$ with $R(-z)^TR(z)=1$, that is, $R=\exp(r)$ with $r=\sum_{k\geq1}r_kz^k$ and $r_k^T=(-1)^{k+1}r_k$. In a $\bbQ$-basis of $H^*(\overline{\mathcal M}_{g,n},\bbQ)$, for fixed $(g,n)$, each conclusion used is a finite family of identities between Laurent polynomials over $\bbQ$ in the $y_i$ and in finitely many free coordinates of the $r_k$. \cite{PPZ15} gives them at all points with $y\in(\bbQ^\times)^r$ and rational $r_k$; as $\bbQ$ is infinite they are identities of Laurent polynomials, and so they hold over $A$.
```

**B4.2** `s02_setup_5.tex` l.20 (`lem:set-pushins` (iii), last sentence). Old: `By \cite[Prop.~2.12]{PPZ15}, $R.\Omega$
is again a cohomological field theory with the unit $1$, to which (iii) applies.` New: `By \cite[Prop.~2.12]{PPZ15},
$R.\Omega$ is again a cohomological field theory with the unit $1$, to which (iii) applies; over a ring other than
$\bbQ$ this is used only for the topological field theories of \cref{def:set-graphsum}, last paragraph.`

**B4.3** `s10_main_2.tex` l.27, the parenthesis. Old: `(The results of \cite{PPZ15} are stated over $\bbQ$; the proofs
of their Props.~2.4, 2.12 and 2.14, and of Props.~2.8 and~2.9, which the proof of Prop.~2.14 uses, are identities between
finite sums, polynomial in the coefficients of the $R$-matrices and of the CohFT, and hold verbatim over
$\bbQ(\zeta)\lser{s}$, \cref{def:set-graphsum}.)` New: `(The results of \cite{PPZ15} are stated over $\bbQ$;
\cite[Props.~2.12 and~2.14]{PPZ15} hold over $\bbQ(\zeta)\lser s$ for the topological field theory $\Cexp(\TFT)$, whose
unit has the components $y_i=\Cexp(\Delta_i^{-1/2})$, units of $\bbQ(\zeta)\lser s$, in the orthonormal basis
$\tilde e_i$, and for the symplectic series $M$ and $X$, by the last paragraph of \cref{def:set-graphsum}.)`

**B4.4** `s01_intro_3.tex` l.9, the list of [PPZ15] in (Q). Delete `Prop.~2.4; ` and `Props.~2.8 and~2.9; ` from the
list. Old (last sentence of the item): `Props.~2.4, 2.8 and~2.9 enter only through the proof of Prop.~2.14, and (8) only
through \cref{lem:set-pushins}.` New: `Props.~2.12 and~2.14 are used as statements over $\bbQ$ and transferred to the
other coefficient rings in \cref{def:set-graphsum}; (8) enters only through \cref{lem:set-pushins}.`

**B4.5** `s10_main_2.tex` l.36. Old: `Prop.~2.14, with Props.~2.4, 2.8 and~2.9 through its proof]{PPZ15}` New: `Prop.~2.14]{PPZ15}`

---------------------------------------------------------------------------------------------------------------------

### B5. Galois descent of F_g^N: rationality of the colour-0 rows is cited from the wrong place  [X4 m1; U06]

Confirmed in the text. `s10_main.tex` l.81: "By \cref{lem:gau-colrows}, $\Qr^{(i)}_{a,n}=\zeta^{i(a+1-n)}\Qr^{(0)}_{a,n}$
with $\Qr^{(0)}_{a,n}\in\Nring$". `lem:gau-colrows` gives the covariance only; `thm:gau-R0` (iii) and `cor:gau-rowsNN`
(i) give $A^{\mathrm{act}}_n\in\bbQ(\zeta)[L,L^{-1}]$, that is $\Qr^{(0)}_{a,n}\in\Nring\otimes\bbQ(\zeta)$. Rationality
needs $\Dn_0\in1+z\bbQ\fps z$, which is `thm:con-glue`.

Checked: `thm:gau-R0` (iii) gives A^act_n = sum_{m<=n} A^can_m [z^{n-m}] D_0^{-1} with A^can_m in Q[L, L^{-1}] ((ii));
`thm:con-glue` states D_0 in 1 + zQ[[z]] and precedes Section 10; no circularity (its proof does not use F_g^N).

**B5.1** `s10_main.tex` l.81. Old: `By \cref{lem:gau-colrows}, $\Qr^{(i)}_{a,n}=\zeta^{i(a+1-n)}\Qr^{(0)}_{a,n}$ with
$\Qr^{(0)}_{a,n}\in\Nring$, and $\Delta_i^{1/2}=\zeta^{4i}C_0/L$;` New: `By \cref{cor:gau-rowsNN}~(iii), which rests on
\cref{lem:gau-colrows}, $\Qr^{(i)}_{a,n}=\zeta^{i(a+1-n)}\Qr^{(0)}_{a,n}$; and the colour-$0$ representatives lie in
$\Nring$: by \cref{thm:con-glue}, $\Dn_0\in1+z\bbQ\fps z$, so
$A^{\mathrm{act}}_n=\sum_{m\leq n}A^{\mathrm{can}}_m\coef{z^{n-m}}\Dn_0^{-1}$ lies in $\bbQ[L,L^{-1}]$
(\cref{thm:gau-R0}~(ii), (iii)), and the representatives of \cref{cor:gau-rowsNN}~(i) are obtained from
$A^{\mathrm{act}}_n(L)\,c_0^{-1}$ by a recursion with coefficients in $\Nring$. Moreover
$\Delta_i^{1/2}=\zeta^{4i}C_0/L$;`

**B5.2 (optional)** `s03_gauge_8.tex` l.36 (`cor:gau-rowsNN` (i)), add at the end: `(After \cref{thm:con-glue},
$A^{\mathrm{act}}_n\in\bbQ[L,L^{-1}]$, and the colour-$0$ representatives lie in $\Nring$.)`
