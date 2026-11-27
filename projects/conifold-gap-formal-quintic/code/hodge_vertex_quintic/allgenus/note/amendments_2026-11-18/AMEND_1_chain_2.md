# Suggested amendments 1: the proof chain. Part 2 of 4 (A2, A3)

Conventions as in `AMEND_1_chain.md`. Bibliography keys: DOSS = `DuninBarkowski:2012bw`, DNOPS = `Dunin-Barkowski:2015caa`.

### A2. The edge condition of the quoted TR theorem: DNOPS (2-33) is misprinted  [E5; U09-6, L3-2/3/4, X3-5]

Confirmed in the text.
* `s06_model_3.tex` l.23 (proof of `thm:mod-tr`, Step 1): "gives for a compact curve with simple branch points, $\dd y\neq 0$
  there, and the Bergman kernel: if the matrix $R$ of (2-31) is symplectic in the orthonormal basis and (2-32) holds, then"
  (the edge condition is not among the hypotheses);
* l.28: "the compatibility condition \cite[Prop.~4.3 (1)]{...} holds, since the curve is rational and $f$ has a pole at
  $\infty$ (...); and by the sentence preceding \cite[Thm.~8.1]{...} this curve satisfies (2-30)--(2-33)."
* `s01_intro_3.tex` l.9 (Q): "the compatibility condition \cite[Prop.~4.3 (1)]{...} and the sentence preceding
  \cite[Theorem~8.1]{...}, applied to $f=4\Wv$, that is, to \cite[(8-1)]{...} with $n=4$ and $a=(0,0,0,0,4v)$".

The disagreement, and the recommendation. Two different conditions are involved.
* The *compatibility condition* DNOPS Prop. 4.3 is equivalent to (2-32). The note proves (2-32) directly
  (`prop:mod-Rtr` (ii)), so Prop. 4.3 is redundant. U09 (M6), L3 (F4) and E5 all say this; there is no disagreement.
* The *edge condition* (2-33) relates the bidifferential B to R. It is independent of (2-32) and is a hypothesis of the
  quoted theorem (in DOSS Theorem 4.1 the local curve is *defined* by (4-1)). The proof of (2-32) does not make it
  redundant. L3 did not notice the misprint and accepted "always satisfied if the curve is compact"; U09 (M1) saw the
  misprint and relied on the correct statement being known (Eynard), which is not a declared input; E5 proves it.
Recommendation: follow E5. Drop Prop. 4.3 and the sentence before Theorem 8.1 from the inputs, and prove the edge
condition for the model curve in the note (a new part (iii) of `prop:mod-Rtr`, 15 lines). The alternative, to add
[Eynard, arXiv:1110.2949, App. B] to (Q), adds an input and a bibliography entry for something with a short proof.

Checked. E5's three-step proof was re-derived by hand: the coefficient of D(B), its regularity on the diagonal (value
f'''/6), the cancellation of the simple poles at the branch points, the principal part and (dp/dw_k)^2 = -1/x''(c_k);
the diagonal case D(dw dw'/(w-w')^2) = dw dw'/(w^2 w'^2); and the signs of the two Laplace transforms against the
definition of R_TR in `prop:mod-Rtr`. The misprint itself and the DOSS numbers (2-22)-(2-26), (4-1)-(4-3) were only
spot-checked for existence in the stored text extracts (`audit_2026-11-18/work/L3/`); E5's numerics were not rerun.

**A2.1** `s06_model_2.tex` l.51, in `prop:mod-Rtr`: add a third item after item (ii) (before `\end{enumerate}`).

```latex
\item (the edge condition) For a double series $\Omega=\sum_{m,n\geq-M}a_{mn}\,w_i^m\,w_j'^{\,n}\,\dd w_i\,\dd w_j'$ in the coordinates $w_i$ at $\varsigma_i$ in the first point and $w_j$ at $\varsigma_j$ in the second, with $a_{-1,n}=a_{m,-1}=0$, put $(\mathsf L_i\otimes\mathsf L_j)[\Omega](\vartheta_1,\vartheta_2):=\sum_{m,n}a_{mn}\mathsf G_m(\vartheta_1)\mathsf G_n(\vartheta_2)$, and let $B^{\mathrm{reg}}_{ij}:=\omega_{0,2}-\delta_{ij}\,\dd w_i\,\dd w_i'/(w_i-w_i')^2$, which is holomorphic at $(\varsigma_i,\varsigma_j)$. Then, with $\psi_a:=\vartheta_a^{-1}$,
\begin{equation}
\label{eq:mod-edge}
(\psi_1+\psi_2)\,(\mathsf L_i\otimes\mathsf L_j)\bigl[B^{\mathrm{reg}}_{ij}\bigr]=\delta_{ij}-\sum_kR_{\mathrm{TR}}^{-1}(\psi_1)^i{}_k\,R_{\mathrm{TR}}^{-1}(\psi_2)^j{}_k\qquad\text{for all }i,j .
\end{equation}
This is the edge condition \cite[(4-1), with (3-2) and (3-57)]{DuninBarkowski:2012bw}: the expansion of $\omega_{0,2}$ at the branch points is the edge of $R_{\mathrm{TR}}$. The form printed in \cite[(2-33)]{Dunin-Barkowski:2015caa} omits the term $\delta_{ij}$ and is not used.
```

**A2.2** `s06_model_3.tex` l.19, proof of `prop:mod-Rtr`: add after the paragraph "(ii) ..." (before `\end{proof}`).

```latex
(iii) Put $f:=1/x_\kappa'(p)=p^2/(\kappa(p^5-v))$ and, for a symmetric bidifferential $\Omega$, $\mathsf D\Omega:=\dd_1(\Omega/\dd x_{\kappa,1})+\dd_2(\Omega/\dd x_{\kappa,2})$, where $\Omega/\dd x_{\kappa,a}$ is a function of the $a$-th point (local letters). \emph{Step A:} $\mathsf D\omega_{0,2}=\sum_k\dd\xi^k(p_1)\,\dd\xi^k(p_2)$. Indeed
\[
\mathsf D\omega_{0,2}=\Bigl[\frac{f'(p_1)+f'(p_2)}{(p_1-p_2)^2}+\frac{2\bigl(f(p_2)-f(p_1)\bigr)}{(p_1-p_2)^3}\Bigr]\dd p_1\,\dd p_2 .
\]
Fix $p_2\neq0,\infty$ off the branch points. As a $1$-form in $p_1$ this is regular at $p_1=p_2$ (the bracket tends to $f'''(p_2)/6$) and holomorphic at $p_1=0$ and $p_1=\infty$ ($f$ has a double zero at $0$ and is $O(p^{-3})$ at $\infty$). At $p_1=\varsigma_k$, where $f=r_k/(p-\varsigma_k)+O(1)$ with $r_k:=1/x_\kappa''(\varsigma_k)$, the simple poles of the two terms cancel, and the principal part is $-r_k\,\dd p_1\,\dd p_2/\bigl((p_1-\varsigma_k)^2(\varsigma_k-p_2)^2\bigr)$. By Step~1, $\dd\xi^k(p_1)\,\dd\xi^k(p_2)=\bigl(\tfrac{\dd p}{\dd w_k}(\varsigma_k)\bigr)^2\dd p_1\,\dd p_2/\bigl((\varsigma_k-p_1)^2(\varsigma_k-p_2)^2\bigr)$, which is holomorphic in $p_1$ off $\varsigma_k$ and has the same principal part there, because $x_\kappa=x_\kappa(\varsigma_k)-w_k^2/2$ gives $x_\kappa''(\varsigma_k)\bigl(\tfrac{\dd p}{\dd w_k}(\varsigma_k)\bigr)^2=-1$. The difference of the two sides is a holomorphic $1$-form in $p_1$ on $\bbP^1$, hence zero. \emph{Step B.} If $\Omega$ is holomorphic at $(\varsigma_i,\varsigma_j)$, then $(\mathsf L_i\otimes\mathsf L_j)[\mathsf D\Omega]=-(\vartheta_1+\vartheta_2)(\mathsf L_i\otimes\mathsf L_j)[\Omega]$: apply \cref{lem:mod-laplace}~(i) in the first point to $F=\Omega/\dd x_{\kappa,1}$, which has at most a simple pole at $\varsigma_i$ while $F\,\dd x_{\kappa,1}=\Omega$ has no residue, coefficientwise in the second point, and likewise in the second point. \emph{Step C.} By the definition of $R_{\mathrm{TR}}$, $\mathsf L_j[\dd\xi^k]=-\vartheta\,R_{\mathrm{TR}}^{-1}(\vartheta^{-1})^j{}_k$. For $i\neq j$, $B^{\mathrm{reg}}_{ij}=\omega_{0,2}$, and Steps~A and~B give $-(\vartheta_1+\vartheta_2)(\mathsf L_i\otimes\mathsf L_j)[\omega_{0,2}]=\vartheta_1\vartheta_2\sum_kR_{\mathrm{TR}}^{-1}(\psi_1)^i{}_kR_{\mathrm{TR}}^{-1}(\psi_2)^j{}_k$; divide by $-\vartheta_1\vartheta_2$. For $i=j$ use the coordinate $w=w_i$, in which $1/(\dd x_\kappa/\dd w)=-1/w$: the formula of Step~A with $-1/w$ in place of $f$ gives $\mathsf D\bigl(\dd w\,\dd w'/(w-w')^2\bigr)=\dd w\,\dd w'/(w^2w'^2)$, as $(w^{-2}+w'^{-2})/(w-w')^2-2/(ww'(w-w')^2)=1/(w^2w'^2)$. So $\mathsf DB^{\mathrm{reg}}_{ii}=\sum_k\dd\xi^k\,\dd\xi^k-\dd w\,\dd w'/(w^2w'^2)$, and by Step~B and $\mathsf G_{-2}=-\vartheta$, $-(\vartheta_1+\vartheta_2)(\mathsf L_i\otimes\mathsf L_i)[B^{\mathrm{reg}}_{ii}]=\vartheta_1\vartheta_2\bigl(\sum_kR_{\mathrm{TR}}^{-1}(\psi_1)^i{}_kR_{\mathrm{TR}}^{-1}(\psi_2)^i{}_k-1\bigr)$; divide by $-\vartheta_1\vartheta_2$. (At $\psi_2=-\psi_1$, \eqref{eq:mod-edge} gives the symplecticity of (i) again.)
```

**A2.3** `s06_model_3.tex` l.23 (proof of `thm:mod-tr`, Step 1). Old: `if the matrix $R$ of (2-31) is symplectic in
the orthonormal basis and (2-32) holds, then for $2h-2+n>0$` New: `if the matrix $R$ of (2-31) is symplectic in the
orthonormal basis, (2-32) holds, and the expansion of the Bergman kernel at the branch points is the edge of $R$, in the
form \cite[(4-1)]{DuninBarkowski:2012bw}, that is \eqref{eq:mod-edge} (the printed \cite[(2-33)]{Dunin-Barkowski:2015caa}
omits the term $\delta_{ij}$), then for $2h-2+n>0$`

**A2.4** `s06_model_3.tex` l.28. Old: `We apply it to $f:=4\Wv$, the superpotential \cite[(8-1)]{Dunin-Barkowski:2015caa}
with $n=4$, $a=(0,0,0,0,4v)$: the critical points $p^5=v$ are simple with distinct critical values, $\dd y\neq 0$, and the
shift $p_{\mathrm{Dub}}=p$ is trivial ($a_1=0$); the compatibility condition \cite[Prop.~4.3 (1)]{Dunin-Barkowski:2015caa}
holds, since the curve is rational and $f$ has a pole at $\infty$ ($\dd y/\dd f=p^2/(4(p^5-v))$ has poles only at the zeros
of $\dd f$); and by the sentence preceding \cite[Thm.~8.1]{Dunin-Barkowski:2015caa} this curve satisfies (2-30)--(2-33).
The hypotheses on $R$ are \cref{prop:mod-Rtr} at $\kappa=4$, which gives` New:

```latex
We apply it to $(\bbP^1,4\Wv,p,\omega_{0,2})$ (for orientation, $4\Wv$ is the superpotential \cite[(8-1)]{Dunin-Barkowski:2015caa} with $n=4$, $a_1=\dots=a_4=0$ and $a_5=4v$): the critical points $p^5=v$ are simple with pairwise distinct critical values, and $\dd y\neq0$ there (\cref{def:mod-curve}). The three hypotheses are \cref{prop:mod-Rtr} at $\kappa=4$: symplecticity is (i), (2-32) is (ii), and the edge condition is (iii); none of them is taken from the literature. (The compatibility condition \cite[Prop.~4.3]{Dunin-Barkowski:2015caa} is equivalent to (2-32); it holds by item (2) there, as $\dd y/\dd(4\Wv)=p^2/(4(p^5-v))$ has poles only at the zeros of $\dd\Wv$, and is not used.) \Cref{prop:mod-Rtr} gives
```

Optional simplification (not needed): with (iii) the quoted theorem can be applied directly at kappa = 1; Step 2
("transfer to kappa = 1") then becomes a consistency check and may be cut, with `\alpha^1`, `\Delta^1_j`, `\xi^{(1,j)}`
defined from `prop:mod-Rtr` at kappa = 1. Recommended: keep Step 2, it is correct and costs nothing.

**A2.5** `s01_intro_3.tex` l.9, item (Q). Old: `Topological recursion, applied only to the model curve $x=\Wv(p)$,
$y=p$ (\cref{sec:mod}):` New: `Topological recursion, applied only to the model curve $x=\Wv(p)$, $y=p$ and its
rescalings $x=\kappa\Wv$ (\cref{sec:mod}):`
Old: `the theorem \cite[Theorem~4.1]{DuninBarkowski:2012bw}, with its graph weights \cite[\S2.3]{DuninBarkowski:2012bw},
as recalled in \cite[\S2.6, (2-30)--(2-35)]{Dunin-Barkowski:2015caa}, with the action of
\cite[\S2.3]{Dunin-Barkowski:2015caa} (its dilaton leaves taken with the sign of \cite[Def.~2.13]{PPZ15}), the
compatibility condition \cite[Prop.~4.3 (1)]{Dunin-Barkowski:2015caa} and the sentence preceding
\cite[Theorem~8.1]{Dunin-Barkowski:2015caa}, applied to $f=4\Wv$, that is, to \cite[(8-1)]{Dunin-Barkowski:2015caa} with
$n=4$ and $a=(0,0,0,0,4v)$, whose critical points $p^5=v$ have the pairwise distinct critical values $5p_c^4$;` New:

```latex
the theorem \cite[Theorem~4.1, with (4-1)--(4-3)]{DuninBarkowski:2012bw}, with its graph weights \cite[\S2.3, (2-22)--(2-26)]{DuninBarkowski:2012bw} (the dilaton leaves (2-23) and the factor $1/|\Aut|$ of (2-26)), as recalled in \cite[\S2.6, (2-30)--(2-32), (2-34), (2-35)]{Dunin-Barkowski:2015caa}, with matrices in the left-action convention of \cite[Remark~2.1]{Dunin-Barkowski:2015caa} (upper index the row), and with the action of \cite[\S2.3]{Dunin-Barkowski:2015caa} read with the factor $1/|\Aut\Gamma|$ and with the dilaton leaves of \cite[Def.~2.13]{PPZ15}, that is, as the action of \cite[Defs.~2.2, 2.13]{PPZ15}; it is applied to $x=4\Wv$, $y=p$, whose critical points $p^5=v$ have the pairwise distinct critical values $5p_c^4$. Its hypotheses, namely symplecticity, \cite[(2-32)]{Dunin-Barkowski:2015caa} and the edge condition \cite[(4-1)]{DuninBarkowski:2012bw}, are proved for this curve in \cref{prop:mod-Rtr}; the form \cite[(2-33)]{Dunin-Barkowski:2015caa} of the edge condition is misprinted (it omits the term $\delta_{ij}$) and is not used, nor are \cite[Prop.~4.3]{Dunin-Barkowski:2015caa} and the sentence preceding \cite[Theorem~8.1]{Dunin-Barkowski:2015caa};
```

**A2.6** `s06_model_2.tex` l.43 (`prop:mod-Rtr`). Old: `$R_{\mathrm{TR}}^{-1}(\vartheta^{-1})^j{}_i:=-\vartheta^{-1}\mathsf L_{\varsigma_j}[\dd\xi^i]$ (upper index the row).`
New: `$R_{\mathrm{TR}}^{-1}(\vartheta^{-1})^j{}_i:=-\vartheta^{-1}\mathsf L_{\varsigma_j}[\dd\xi^i]$ (upper index the row: the left-action convention of \cite[Remark~2.1]{Dunin-Barkowski:2015caa}, in which (2-32) is $R^{-1}$ applied to the unit \cite[Section~3]{Dunin-Barkowski:2015caa}).`

**A2.7** `s06_model_2.tex` l.55. Old: `With legs $R^{-1}(\psi)$ the sign is not a convention: it is the one for which
the action preserves the unit (\cref{sec:sgn-R}).` New: `With legs $R^{-1}(\psi)$ the sign is fixed by the requirement
that the action preserve the unit: it is that of \cite[(2-23)]{DuninBarkowski:2012bw} and of
\cite[Prop.~2.12]{PPZ15} (\cref{sec:sgn-R}).`

**A2.8** `b_signs.tex` l.14 (appendix B). Old: `The action of \cite[\S2.3]{Dunin-Barkowski:2015caa} is the same (legs
$R^{-1}(\psi)$ acting on column vectors, the row index being the outgoing colour);` New: `The action of
\cite[\S2.3]{Dunin-Barkowski:2015caa} is the same up to the factor $1/|\Aut\Gamma|$, omitted there and present in
\cite[(2-26)]{DuninBarkowski:2012bw} and \cite[Def.~2.2]{PPZ15}, and up to the sign of the dilaton leaves (legs
$R^{-1}(\psi)$ acting on column vectors, the row index being the outgoing colour);`

**A2.9** `b_signs.tex` l.52 and l.58 (table of appendix B). l.52, old: `at $a=(0,0,0,0,4v)$, $f=4\Wv$;` new:
`at $(a_1,\dots,a_5)=(0,0,0,0,4v)$, $f=4\Wv$;` and third column, old: `the quoted theorem of \cref{thm:mod-tr},
transferred to $x=\Wv$ by \cref{lem:mod-scale}` new: `the curve to which the quoted theorem is applied in
\cref{thm:mod-tr}, transferred to $x=\Wv$ by \cref{lem:mod-scale}; \cite[\S8]{Dunin-Barkowski:2015caa} itself is not used`.
l.58, old: `on the line $a=(0,0,0,0,4v)$` new: `on the line $(a_1,\dots,a_5)=(0,0,0,0,4v)$`.

No change needed: `s01_intro_3.tex` l.15 (the Laplace prescription), l.19 and `s10_main.tex` l.31, l.33 ("as recalled
in DNOPS") stay correct; `b_signs_2.tex` l.25 stays correct.

---------------------------------------------------------------------------------------------------------------------

### A3. What is proved is a formal change of frame, not an analytic continuation  [E7; X4 m2-m4, U18]

Confirmed in the text.
* `main.tex` l.18: "Near the conifold point $u=1-5^5q=0$ we factorise the $R$-matrix as $R=R'\,\Mo(\tau(u)\zh)$".
* `s02_setup_3.tex` l.28: "near the conifold point they are Laurent series in $s$, ...; the passage is the conifold
  expansion of \cref{def:set-confexp}."
* `b_signs_2.tex` l.106, end: "Consequently, for every pair of the family of \cite[Rem.~1.4]{CGYZ26} ..., the formal
  quintic has the conifold gap with the constant of \cite[Conj.~1.5]{CGYZ26}."; l.108: "whether $F_g^{\Nring}$ lies in
  the algebra of \cite[Def.~1.1]{CGYZ26} is not needed."

Checked. E7's membership argument (2.3 (iii)) was re-derived: N_0 = B[x] with W = c_0^2 + c_0 pi_1; pi_2 =
(DW - 2 c_0 c_1)/c_0; L -> zeta L commutes with D and multiplies Q^{(0)}_{a,n} by zeta^{a+1-n} once A^act_n = L^{4n+1}
ph_n(L^{-5}) (Theorem C); the action on Delta_i^{1/2} and Psi_Q. E7's proof that ev is injective (Beukers-Heckman,
Schlesinger) was read but not checked against the sources; it is not needed and is kept out of the text (see D9).
One caution on wording: "it is not an analytic continuation" is itself a statement (E7 2.2 proves it from the
irreducibility of the monodromy of PF4, which the note does not prove). The texts below say only that no analytic
continuation is involved or claimed.

**A3.1** `main.tex` l.18 (abstract), three replacements.
Old: `Near the conifold point $u=1-5^5q=0$ we factorise the $R$-matrix as $R=R'\,\Mo(\tau(u)\zh)$, where`
New: `We expand the theory at the conifold point $u=1-5^5q=0$ in a conifold frame: the graph sum is written in a differential ring generated by the periods, and a regular period and the vanishing period at $u=0$ are substituted for the two periods at $q=0$. This is a formal change of frame (the holomorphic limit in the conifold frame); no analytic continuation is involved. We factorise the expanded $R$-matrix as $R=R'\,\Mo(\tau(u)\zh)$, where`
Old: `the term of every leg-free stable graph of genus $g$ with at least one edge is regular at the conifold point,`
New: `in the graph sum of $R'$ on the theory obtained from the topological part by the action of $\Mo(\tau(u)\zh)$, the term of every leg-free stable graph of genus $g$ with at least one edge is regular at the conifold point,`
Old: `Summing over graphs gives the conifold gap for the formal quintic in every genus $g\geq2$.`
New: `Summing over graphs gives the conifold gap for this expansion of the genus-$g$ potential of the formal quintic, for every $g\geq2$.`
(The second replacement answers E7's side remark: the terms of the Givental--Teleman sum of R itself are not regular.)

**A3.2** `s02_setup_3.tex` l.28 (`def:set-points`). Old: `near the conifold point they are Laurent series in $s$, with
$\Dq=-\frac{1-s^5}5s^{-4}\frac{\dd}{\dd s}$; the passage is the conifold expansion of \cref{def:set-confexp}.` New:
`near the conifold point the note uses their conifold expansions (\cref{def:set-confexp}), which are Laurent series in
$s$, with $\Dq=-\frac{1-s^5}5s^{-4}\frac{\dd}{\dd s}$. The conifold expansion substitutes the periods of a conifold frame
for $C_0$ and $I_1$ in a representative in the ring $\Nring$; it is a formal change of frame, and no analytic
continuation is involved.`

**A3.3** `s01_intro_3.tex` l.31, new item after item 2 of "What is not claimed".

```latex
\item Nothing is claimed about the analytic continuation of the series $F_g(q)$ to a neighbourhood of $u=0$. The conifold expansion $\Cexp_{\alpha,\beta}$ is a formal change of isotropic frame on the ring $\Nring$ (\cref{def:set-confexp}), applied to the graph-sum representative $F_g^{\Nring}$ of $F_g$ (\cref{def:mai-Fcf}); this is the sense of ``conifold expansion'' and of ``conifold gap'' in \cref{cor:int-gap}.
```

**A3.4** `b_signs_2.tex` l.106 and l.108 (`sec:sgn-lit`). Two versions; (a) is recommended, (b) is the fallback if the
definitions of [CGYZ26, Def. 1.1-1.2], which the note recalls from an HTML rendering (draft note at `rmk:mai-constant`),
turn out to differ from the recollection at l.100.

(a) l.106 unchanged. l.108, old: `\emph{Differences.} The rows involve $L\notin\bbQ(q)$, yet $F_g^{\mathrm{cf}}$ is a
Laurent series in $u$ by Corollary~B; whether $F_g^{\Nring}$ lies in the algebra of \cite[Def.~1.1]{CGYZ26} is not
needed.` New:

```latex
\emph{Membership.} The comparison needs that $F_g^{\Nring}$ lies in the differential subring $\mathscr S$ of $\Nring$ generated by $c_0^{\pm1}$ and $(\mathsf W/c_0^2)^{\pm1}$ over $\bbQ[L^5,L^{-5}]$, on which the two expansions agree; then $F_g=\ev(F_g^{\Nring})$ lies in the algebra of \cite[Def.~1.1]{CGYZ26}, and $F_g^{\mathrm{cf}}$ is its conifold expansion in the sense of \cite[Def.~1.2]{CGYZ26}. This holds (the letters $\mathscr S$, $\mathsf B$, $x$, $\pi_j$, $\varrho$ are local). (1) $F_g^{\Nring}$ is free of $\logq$ (\cref{cor:gau-rowsNN}~(ii)) and fixed by the automorphisms $i_j\mapsto i_j-kc_j$, $\logq\mapsto\logq+k$, $k\in\bbQ$ (\cref{def:set-confexp}). With $x:=i_0/c_0$ and $\pi_j:=i_j-xc_j$ ($j=1,2$), the $\logq$-free subring of $\Nring$ is the polynomial ring $\mathsf B[x]$ over $\mathsf B:=\bbQ[L^{\pm1}][c_0^{\pm1},c_1,c_2,c_3,\pi_1,\pi_2][\mathsf W^{-1}]$, because $\mathsf W=c_0^2+c_0\pi_1$; the automorphisms fix $\mathsf B$ and send $x$ to $x-k$, so $F_g^{\Nring}\in\mathsf B$. (2) Let $\varrho$ be the automorphism of $\Nring\otimes\bbQ(\zeta)$ with $L\mapsto\zeta L$ which fixes $\zeta$ and the other generators; it commutes with $\Dq$, since $\Dq L=L(L^5-1)/5$ and $\Dq c_3$, $\Dq i_2$ involve $L$ only through $L^5$. By \cref{thm:gam-D}, $A^{\mathrm{act}}_n=A^{\mathrm{can}}_n=L^{4n+1}\mathrm{ph}^{\mathrm{can}}_n(L^{-5})$ (\cref{thm:gau-R0}~(ii), (iii)), so $\varrho$ multiplies $\Qr^{(0)}_{0,n}=A^{\mathrm{act}}_nc_0^{-1}$ by $\zeta^{1-n}$ and, by the recursion of \cref{cor:gau-rowsNN}~(i), whose coefficients it fixes, $\Qr^{(0)}_{a,n}$ by $\zeta^{a+1-n}$. Hence $\varrho(\Qr^{(i)}_{a,n})=\Qr^{(i+1)}_{a,n}$, $\varrho(\Delta_i^{1/2})=\Delta_{i+1}^{1/2}$ and $\varrho\bigl((\Psi_{\mathrm Q})_{aj}\bigr)=(\Psi_{\mathrm Q})_{a,j+1}$ (\cref{lem:mai-Psi}~(i)): $\varrho$ permutes the colours and fixes $F_g^{\Nring}$, which therefore contains only powers of $L^5$. (3) $\bbQ[L^{\pm5}][c_0^{\pm1},c_1,c_2,c_3,\pi_1,\pi_2][\mathsf W^{-1}]\subset\mathscr S$, as $c_j=\Dq^jc_0$, $\pi_1=(\mathsf W-c_0^2)/c_0$ and $\pi_2=(\Dq\mathsf W-2c_0c_1)/c_0$.

\emph{Differences.} The rows involve $L\notin\bbQ(q)$, yet $F_g^{\mathrm{cf}}$ is a Laurent series in $u$ by Corollary~B.
```

(b) l.106, old: `the formal quintic has the conifold gap with the constant of \cite[Conj.~1.5]{CGYZ26}.` new: `the
formal quintic has the conifold gap in the sense of Corollary~B, with the constant of \cite[Conj.~1.5]{CGYZ26}.` l.108,
old: `whether $F_g^{\Nring}$ lies in the algebra of \cite[Def.~1.1]{CGYZ26} is not needed.` new: `whether
$F_g^{\Nring}$ lies in the differential subring of $\Nring$ on which the two expansions agree is not used here; it is
what would identify $F_g^{\mathrm{cf}}$ with the conifold expansion of the function $F_g$ in the sense of
\cite[Def.~1.2]{CGYZ26}.`

In version (a) the membership uses Theorem C (through `thm:gam-D`), so it must stay in the appendix, after Section 9;
it is used by no proof.
