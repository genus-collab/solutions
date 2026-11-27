# Suggested amendments, part 2: minor findings of units U01-U08 (file 2 of 2: P3, R, skipped)

Continues `AMEND_2_units_01_08.md` (conventions and classes are explained there). Paths are relative to `sections/`.

---------------------------------------------------------------------------------------------------------------------

## P3 (wording and notation; optional)

### A2-21 [P3] `def:set-confexp`: vanishing of the constant $\Conc(C_0,I_1)$ is justified tersely — `s02_setup_6.tex` l.14
- Reports: U02_C F2.
- Current: "and this constant vanishes, as its part of order $q^0$ involves only the leading coefficient of $\PFfour$ and the jets $(1,0,0,0)$, $(\logq,1,0,0)$ of $C_0$, $I_1$."
- Judgement: the claim and the reason are correct (re-derived: at $q^0$ the jet matrix is the antidiagonal
  $(-1,1,-1,1)$, and $f_0\Jet_{03}g_3=0$), but why only the leading coefficient survives is not said. Not P1: a reason is given.
- Change: replace the quoted clause by
  `and this constant (the constants of $\Dq$ on $\bbQ\fps q[\logq]$ are $\bbQ$, proof of \cref{lem:set-periods}) vanishes: the coefficients of $\PFfour$ below the leading one are multiples of $1-u=5^5q$, so by \cref{lem:gau-conc}~(iii) the part of order $q^0$ of $\Conc(C_0,I_1)$ is $\sum_a(-1)^{a+1}f_ag_{3-a}$ for the jets $f=(1,0,0,0)$, $g=(\logq,1,0,0)$ of $C_0$, $I_1$ at order $q^0$, which is $-f_0g_3=0$.`

### A2-22 [P3] Meaning of functions "at $\Qn=1$, $t=T(q)H$" — `s02_setup.tex` l.78; `s02_setup_2.tex` l.3; `a1_cyclicLR.tex` l.5
- Reports: U01_C F5, U03_A M3, U03B M3.
- Current (`s02_setup.tex` l.78): "Then at $\lambda_i=\zeta^i$, with $\Qn=1$ and $t=T(q)H$,"; (`s02_setup_2.tex` l.3): "Along the line $t=T(q)H$, with $\Qn=1$, functions of $t$ are functions of $q$"; (`a1_cyclicLR.tex` l.5): "Every step is an identity in $\bbQ\fps q$ or $\bbQ\fps q[\logq]$".
- Judgement: agree that the sums over $n$ are not literally sums in $\bbQ\fps q[\logq]$ ($\sum_n(dT)^n/n!$ is not a
  polynomial in $\logq$). The intended reading is in the text, but late: Translation (d) and `lem:set-grading` both say
  "by the divisor equation", and `def:set-genus0` defines $\re^T$. Hence P3, one sentence.
- Change (`s02_setup_2.tex` l.3): after `functions of $t$ are functions of $q$` insert
  `: by the divisor equation (G0)(b) the insertions of $TH$ are resummed, so that the functions lie in $\bbQ[T]\fps{\re^T}$, and $T\mapsto T(q)$, $\re^T\mapsto q\exp(\tilde I_1/C_0)$ (\cref{def:set-genus0}) maps them to $\bbQ\fps q[\logq]$, intertwining $\partial_T$ with $C_1^{-1}\Dq$`
- Change (`s02_setup.tex` l.78): `with $\Qn=1$ and $t=T(q)H$,` -> `with $\Qn=1$ and $t=T(q)H$ (read through the divisor equation, with $\re^{dT}=(\re^T)^d$ as in \cref{def:set-genus0}; see (d) below),`
- Change (`a1_cyclicLR.tex` l.5): `Every step is an identity in $\bbQ\fps q$ or $\bbQ\fps q[\logq]$;` -> `Every step is an identity over $\bbQ\fps q$ or $\bbQ\fps q[\logq]$ (in \cref{lem:set-periods} with the further variables $\mathfrak w$, $z^{-1}$, in \cref{lem:set-order3} over $\bbQ\lser q$);`

### A2-23 [P3] `def:set-fq`: regularity of the localisation formula at $\lambda_i=\zeta^i$ — `s02_setup.tex` l.24
- Reports: U01_C F6.
- Current: "(The second term, $\bbE^\vee\otimes\cO(5)|_{p_i}$ on the contracted components, has an invertible top Chern class there, as $5\zeta^i\neq0$.)"
- Judgement: agree that this is not the whole justification; the other denominators are also non-zero (I re-checked the
  convexity argument of the report). The definition is the one of the cited sources, so optional.
- Change: append inside the parenthesis: ` The other denominators of the localisation formula are products of the weights $(a\lambda_i+b\lambda_j)/d-\lambda_k$ ($a+b=d$, non-zero before the specialisation), $(\lambda_i-\lambda_j)/d_1+(\lambda_i-\lambda_k)/d_2$ and of classes $(\lambda_i-\lambda_j)/d-\psi$ ($j,k\neq i$); none vanishes at $\lambda_i=\zeta^i$, as a proper convex combination of two distinct fifth roots of unity has modulus $<1$.`

### A2-24 [P3] "$R(0)$" for the value at $z=0$ next to $R(q=0)$, $R(0,z)$ — `s02_setup_3.tex` l.33; `s02_setup_4.tex` l.17 (against l.36)
- Reports: U01_C F3.
- Current: "By \cref{thm:set-GT}, since $R(0)=1$," and "$\mathrm{ph}_0=1$, because $R(0)=1$"; `thm:set-qrr` (l.36) writes "$R(q=0)$" and "$R(0,z)$". Agree, harmless.
- Change (both places): `$R(0)=1$` -> `$R|_{z=0}=1$`

### A2-25 [P3] `def:set-graphsum`: three tacit hypotheses — `s02_setup_4.tex` l.70, l.75-79
- Reports: U01_C F8 (a)-(c). (F8 (d), the "verbatim over every $\bbQ$-algebra" sentence, is E6.)
- Current: "on a free module $V$ of finite rank with pairing $\eta$"; "$\overline F^\Omega_h(q)\coloneqq\sum_n\frac1{n!}\langle(q+z1)^n\rangle_h$"; "$\mathrm{Cont}_\Gamma[E,q;\Omega]$" for a general $E$.
- Judgement: agree, all three hold in every use.
- Change (l.70, also `s02_setup_5.tex` l.5): `with pairing $\eta$` -> `with a perfect symmetric pairing $\eta$`
- Change (l.70): `$\overline F^\Omega_h(q)\coloneqq\sum_n\frac1{n!}\langle(q+z1)^n\rangle_h$` -> `$\overline F^\Omega_h(q)\coloneqq\sum_{n:\,2h-2+n>0}\frac1{n!}\langle(q+z1)^n\rangle_h$`
- Change (end of l.75, before the display): `the contraction` -> `and a bivector series $E$ symmetric under the exchange of its two slots (together with $\psi$, $\psi'$), the contraction`

### A2-26 [P3] `def:gau-onepole`: "with $P_0$ invertible"; "If $\Bbbk$ contains $\zeta$" — `s03_gauge_2.tex` l.9, l.61
- Reports: U04_A M2, M3.
- Judgement: agree; by `thm:gau-U` (b) invertibility over $\Bbbk\fps u$ and over $\Bbbk\lser u$ coincide.
- Change (l.9): `with $P_0$ invertible` -> `with $P_0\in\mathrm{GL}_5(\Bbbk\fps u)$`
- Change (l.61): `If $\Bbbk$ contains $\zeta$,` -> `If $\Bbbk$ contains a primitive fifth root of unity $\zeta$,`

### A2-27 [P3] Proof of `thm:gau-U`: two one-line omissions — `s03_gauge_4.tex` l.9, l.24
- Reports: U04_A M4.
- Current (l.9): "so \cref{lem:gau-lift}, applied to each shift component, gives"; (l.24): "by field operations, derivatives and $\int_0^u$".
- Judgement: agree; both facts are on the page (`lem:gau-comm` (e) "preserves the shift"; `lem:gau-scal` (d)).
- Change (l.9): `applied to each shift component,` -> `applied to each shift component ($\operatorname{pr}$ preserves the shift, \cref{lem:gau-comm}~(e), so $\operatorname{pr}$ of each component vanishes),`
- Change (l.24): `field operations, derivatives and $\int_0^u$;` -> `field operations, derivatives, $\int_0^u$ and the binomial series $\mathfrak d$ (\cref{lem:gau-scal}~(d));`

### A2-28 [P3] The letter $\cn$ has two meanings — `s03_gauge_5.tex` l.20 against `s03_gauge_6.tex` l.25, l.30
- Reports: U05_A F1, U05B F1.
- Current (`s03_gauge_5.tex` l.20): "where $\cn(0)\in\Bbbk^\times$ is as in \cref{thm:gau-U}~(b)"; (`s03_gauge_6.tex` l.30): "$\cn(z):=P_{11}(0,z)$".
- Judgement: agree; $P_{11}(0,0)=\cn_U(0)\,d^U_1(0)$, equal to $\cn_U(0)$ only under (N0). No statement is affected.
- Change (`s03_gauge_5.tex` l.20): `where $\cn(0)\in\Bbbk^\times$ is as in \cref{thm:gau-U}~(b)` -> `where $\cn(0)=p_{0,0}(0)\in\Bbbk^\times$ is the constant of \cref{thm:gau-U}~(b) (\cref{lem:gau-lead}; not the series $P_{11}(0,z)$ of \cref{prop:gau-symp})`

### A2-29 [P3] `lem:gau-vanq` (a): "the pulled-back model" for a general $A$ — `s03_gauge_5.tex` l.81
- Reports: U05_A F3, U05B F3.
- Current: "where $\tilde A$ is the pulled-back model,". By the standing convention (l.7) $\tilde A$ is the model pulled back
  along $v=v(u)$, so the statement is defined; only $A$ is general here. Clarify.
- Change: `where $\tilde A$ is the pulled-back model,` -> `where $\tilde A$ is the model pulled back along $v=v(u)$ (as in the standing convention; any $v\in u\bbQ\fps u$ with $v'(0)\neq0$ would do),`

### A2-30 [P3] `lem:gau-wkb`: the derivation on $\Bbbk\lser s$ — `s03_gauge_7.tex` l.20
- Reports: U06_A F3.
- Current: "or $(\Bbbk\lser s,\Dq,s^{-1})$; in both, $\Dq L=(L^6-L)/5$ and $5^5q=1-L^{-5}$".
- Judgement: agree in part. $q$ is fixed by the quoted clause; $\Dq$ on $\Bbbk\lser s$ is not named. The other half of
  the finding ("(iii) needs $\zeta\in\Bbbk$") is rejected: the first line already assumes $\Bbbk\supseteq\bbQ(\zeta)$ for $i\neq0$.
- Change: `$(\Bbbk\lser s,\Dq,s^{-1})$;` -> `$(\Bbbk\lser s,\Dq,s^{-1})$ with $\Dq=-\frac{1-s^5}5s^{-4}\frac{\dd}{\dd s}$ (\cref{def:set-points});`

### A2-31 [P3] `thm:gau-R0`, `def:gau-Fg`: referent, image, frame — `s03_gauge_8.tex` l.11, l.13, l.70
- Reports: U06_A F2, F5; U06_B F3.
- Judgement: agree, all harmless (the proof at l.28 already says "the image of a unique Laurent polynomial";
  $\Dn_0=1$ in every frame by `thm:gam-D`).
- Change (l.11): `; this is a power series in $u$` -> `; the product $\Ghat(u,\zh)\Rhat^{\mathrm m}(\tau\zh)$ is a power series in $u$`
- Change (l.13): `and $A^{\mathrm{act}}_n\in\bbQ(\zeta)[L,L^{-1}]$ for every $n$.` -> `and every $A^{\mathrm{act}}_n$ is the image of a unique element of $\bbQ(\zeta)[L,L^{-1}]$ under $L\mapsto(1-5^5q)^{-1/5}$, again written $A^{\mathrm{act}}_n$.`
- Change (l.70): `$\Fg{g}$ is the statement $\Dn_0\equiv1\pmod{z^{3g-2}}$;` -> `$\Fg{g}$ is the statement $\Dn_0\equiv1\pmod{z^{3g-2}}$ for the series $\Dn_0$ of \cref{thm:gau-R0} in the conifold frame considered;`

### A2-32 [P3] `lem:tra-theta`: name the domain of $\Theta$ — `s04_transport.tex` l.87
- Reports: U07_A M3, U07B N6. The domain is given implicitly in l.83 ("a sequence ... in $A$, finitely supported in the graded setting").
- Change: `and $\Theta$ is bijective,` -> `and $\Theta$ is bijective on the set of these sequences,`

### A2-33 [P3] `lem:tra-coeff`: "every coefficient" — `s04_transport_2.tex` l.58
- Reports: U07_A M4, U07B N6.
- Change: `and every coefficient $f_j$ of $\fr(k,a)$ has` -> `and every non-zero coefficient $f_j$ of $\fr(k,a)$ has`

### A2-34 [P3] `def:tra-dressing`, example (i): citations — `s04_transport_4.tex` l.14
- Reports: U07_A M6, U07B N7. The example is true and enters no proof; two facts are uncited.
- Change: `(\cref{thm:gau-G}~(c), \cref{prop:gau-symp})` -> `(\cref{prop:gau-graded}~(iii), \cref{thm:gau-G}~(c) with \cref{lem:gau-vanq}~(c), \cref{prop:gau-symp}; the pulled-back model has the one-pole shape for every $v$ with $v'(0)\neq0$, by the computation of \cref{prop:gau-equalprod}~(i))`

### A2-35 [P3] "neither (Sp) nor (A) is used here"; description of `lem:tra-T` — `s04_transport_4.tex` l.40; `b_signs.tex` l.22
- Reports: U07B N1, N2.
- Current (l.40): "(the $T_b$ of \cite[Prop.~2.12]{PPZ15}); neither (Sp) nor (A) is used here. Thus $E_P=E_{R^{\mathrm{vec}}}$ ...";
  (`b_signs.tex` l.22): "The condition (A) is stated in lowered form, \cref{lem:tra-E,lem:tra-T} in vector form."
- Change (l.40): `neither (Sp) nor (A) is used here.` -> `this computation uses neither (Sp) nor (A).`
- Change (`b_signs.tex` l.22): `\cref{lem:tra-E,lem:tra-T} in vector form.` -> `\cref{lem:tra-E} and the leaf point of \cref{lem:tra-T} in vector form; parts (b), (d) of \cref{lem:tra-T} are computed in lowered form.`

### A2-36 [P3] `lem:tra-conv` (c): $\YC$ is not an element of $\cH^u$; "this section" — `a3_uconv.tex` l.22
- Reports: U08_A m4, U08_B m7.
- Current: "Hence $\YC(v(u))$ and $J^{\mathrm C}$ lie in $\cH^u$ and ..., and every pairing and product used in this section".
- Change: `Hence $\YC(v(u))$ and $J^{\mathrm C}$ lie in $\cH^u$` -> `Hence $\YC(v(u))\in\bbQ[z^{-1}]\fps{u}^5$, $J^{\mathrm C}=\eta^{-1}\YC$ lies in $\cH^u$`; `used in this section` -> `used in \cref{sec:tra}`

### A2-37 [P3] `thm:red-main`: "If moreover"; repeated hypothesis — `s05_reduction.tex` l.56, l.59
- Reports: U08_A m7, U08_B m8.
- Judgement: (b) is proved from $\ML{g}$ alone (Step 4). The statement is true as it stands; stating it sharply costs nothing.
  The repetition "$v(0)=0$, $v'(0)\neq0$" (l.56) is harmless: leave.
- Change (l.59): `\item[(b)] If moreover $\ML{g}$ holds,` -> `\item[(b)] If $\ML{g}$ holds,` (in `cor:red-can` keep "If moreover", which is correct there).

### A2-38 [P3] `def:flh-outer` (ii): citation of `lem:tra-conv` (f) outside its setting — `s07_comparison.tex` l.24 (reported by U08_B m6)
- Current: "so these expansions converge $u$-adically (\cref{lem:set-dilaton}, \cref{lem:tra-conv}~(f))".
- Judgement: agree in part: `lem:set-dilaton` (iii), which covers the ring, is already cited; (f) is stated for model correlators only.
- Change: `(\cref{lem:set-dilaton}, \cref{lem:tra-conv}~(f))` -> `(\cref{lem:set-dilaton}~(iii); for the model, \cref{lem:tra-conv}~(f))`

### A2-39 [P3] Notation: the few clashes worth a change
- Reports: U02_C F7, U03_A M2, U03B M2, U04_A M3, U04_B Finding 2, U05_A F2 (letters), U07_A M5, U08_A m6, U08_B m7.
- Judgement: no clash creates an ambiguity in a proof (all reports agree). Worth changing, because the letter is a global macro:
  - $D\coloneqq3h-3+n$ against `\Dq` $=D$: `a3_uconv.tex` l.10 (and its uses, with $D_w$ at l.49), `s02_setup_7.tex` l.63.
    Change: `D := 3h-3+n` -> `\mathsf D_{h,n} := 3h-3+n` (uses accordingly); l.63 `The degree is $D=` -> `The degree is $\mathsf D=`.
  - `macros.tex` l.61, comment of `\cR`: "weighted rings (appendix A)"; the weighted rings are in the stationary-phase
    appendix, and `\cR` is also the generic differential ring of `def:gau-conc`, `cor:gau-cycsys`, `lem:gau-wkb`. Comment only.
- The remaining clashes are listed under R-4.

---------------------------------------------------------------------------------------------------------------------

## R (rejected: finding wrong, already met by the text, or not worth a change)

- **R-1** U01_C F8 (e): `def:set-rows` (c) "is vacuous for the leg-free sums". True, but (c) fixes the convention used by
  `lem:set-pushins` and Section 4; leave.
- **R-2** U08_A m5 (digest U08-5): "$\wt(k,b)\geq3/5$" in `lem:tra-conv` (d) (`a3_uconv.tex` l.23, l.51) is unused. Correct and
  harmless; deleting it is not worth a renumbering of the proof.
- **R-3** U07B N3: `lem:tra-F` (`s04_transport_4.tex` l.138) lists "the $(P^{-1})_k$ and $d^{\pm1}$" although only $P_0^{-1}$
  and non-negative powers of $d$ occur. An over-inclusive list is not false. U07B N4 ("annihilates the vanishing period",
  "$q_T\in\sK$" and "$\tilde q_T\in-z+\sK$" carry the same content): terminology, justified by `cor:tra-edges`, `cor:tra-slice`.
- **R-4** Notation clashes left alone (local, declared or separated by font): $s$ as index in `lem:set-pushins`; $T$, $t$
  overloaded; `\cP`, $Y_b$, $J$, $q_{ab}$ in Appendices A, B; primes on $\mathcal Z'$, $\cn'$; $Y$, $Z$ in `lem:gau-comm`;
  $r$ as residue and as residue class; $R^{(n)}$ in the proof of `thm:gau-U` (b); `\Bprod` for $A$ and $\tilde A$ (equal under (EP));
  $r_n$, $e$ in `lem:gau-HPF`; $L$, $A$, $N$, $\alpha,\beta$, $R$ in Section 4; $N$, $c$ in `lem:tra-conv` (f);
  `\mathrm{Cont}_\Gamma` of `def:red-flatsum` against `\mathrm{Cont}_\Gamma[E,q;\Omega]` (the factor $d^{2g_v-2+n_v}$ is displayed in the definition).
- **R-5** U06_A F3, second half: "`lem:gau-wkb` (iii) needs $\zeta\in\Bbbk$". The statement already begins
  "Let $\Bbbk\supseteq\bbQ(\zeta)$ be a field (for $i=0$, any field $\Bbbk\supseteq\bbQ$)".
- **R-6** U02_C F4, as worded ("`def:flh-outer` defines no action of a map from the colour space to $\bbH$"): it does define
  the edge and the leaf point of such a map (item (ii)). The residual point is kept as A2-19.
- **R-7** Remarks the reports themselves mark as "not findings": the open entry $\Conc(f_0,y_3)=-\frac9{10}$ and the sign
  $u^2\det\Wr_{\mathbf f}=-2$ (U03); "$z$ invertible" in `cor:gau-cycsys` (ii) used only to write $w/z$ (U03); the
  elimination of the input \cite[(1.7)]{AC96} through `lem:set-pushins` (i) (U02, a possible simplification, not a
  correction); `lem:tra-T` (a) stated for all $f$ (U07); `s04_transport_4.tex` l.33, "`\]` In the variables" on one source line.

---------------------------------------------------------------------------------------------------------------------

## Skipped (ruled on elsewhere)

| Finding | Where | See |
|---|---|---|
| U01_C F2 (digest U01-4): `def:set-rows` (a), (b), coefficient ring are undeclared inputs | `s02_setup_4.tex` l.7-15 | E4, README must-fix 1 |
| U02_C F1 (U02-1): all-genus divisor equation in (T) | `s02_setup_5.tex` l.42 | E4 |
| U01_C F8 (d): PPZ "hold verbatim over every commutative $\bbQ$-algebra" | `s02_setup_4.tex` l.79 | E6 |
| U02_C F6 (U02-4): the conifold expansion of $F_g$ depends on the representative | `s02_setup_6.tex` l.20 | E7 |
| U02_C F3 (U02-6): meaning of $\dd^n\overline F_h(q/\varkappa)$ for formal $y$ | `s02_setup_8.tex` l.11-28 | E8 |
| U06_A F4, U06_B F1 (U06-8, U06-10): `def:mai-Fcf` cites `lem:gau-colrows` for rationality | `s10_main.tex` l.81 | README routine list |
| U08_B m8, side remark: `rmk:mai-necessity`, "dressings that are not regular in $u$" | `s10_main_2.tex` l.51 | README "false as written" |

## Not verified by me

- The cited external sources (GJR18, LP20, PPZ15, AC96, Givental) were not opened; no amendment above changes a pin.
- A2-03 relies on the reports' reading that `cor:gau-colours` gives the closing relation for all colours; I read its
  statement ("the expanded rows of every colour form a WKB solution") and proof, not every lemma behind it.
- A2-23: the list of localisation weights follows the report and the standard graph-localisation formula; I checked the
  two non-vanishing arguments, not the list against a source.
- The reports' machine checks (series to finite order) were not rerun; no wolframscript or compilation was used.
- Line numbers were read from the current files on 2026-11-18; the new LaTeX was not compiled.
