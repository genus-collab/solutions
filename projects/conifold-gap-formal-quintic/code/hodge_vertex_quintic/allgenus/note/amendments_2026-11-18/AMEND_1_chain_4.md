# Suggested amendments 1: the proof chain. Part 4 of 4 (C, D, and the inputs after the amendments)

Conventions as in `AMEND_1_chain.md`.

## C. False as written, used by no proof

(Two more statements that are false as worded, `lem:fsp-param` (i) and `lem:fsp-res` (i), are B3.4 and B1.1.)

### C1. `rmk:set-qrrsign`: ph_1 is not forced independently of (QRR)  [U01_C F1, L5 F3, U05B F5]

Confirmed: `s02_setup_4.tex` l.52, last sentence: "That $\mathrm{ph}_1=-\frac3{20}$ is forced, independently of (QRR), is
part of \cref{lem:gau-HPF}." `lem:gau-HPF` (iii) forces it only "If every $\mathrm{ph}_n$ lies in $\bbQ\fps u$"; the
actual rows have coefficients in power series in q. Checked: by `lem:gau-HPF` (ii) with r_1 = 3/5 the equation is
(4 - 5 theta_u) ph_1 = -3/5, with the solutions -3/20 + k u^{4/5}, and u^{4/5} = (1 - 5^5 q)^{4/5} is a power series in q.

Old: `That $\mathrm{ph}_1=-\frac3{20}$ is forced, independently of (QRR), is part of \cref{lem:gau-HPF}.`
New: `The value $h_1=-\frac3{20}$ comes from (QRR). It agrees with \cref{lem:gau-HPF}~(iii), which forces $\mathrm{ph}_1=-\frac3{20}$ for a solution with coefficients in $\bbQ\fps u$; but this is not an independent check of the sign: over $\bbQ\fps q$ the equation $(4-5\theta_u)\mathrm{ph}_1=-\frac35$ of \cref{lem:gau-HPF}~(ii) also has the solutions $-\frac3{20}+k\,u^{4/5}$, $u^{4/5}=(1-5^5q)^{4/5}$ (the freedom of \cref{lem:gau-wkb}~(ii)), and that the actual $\mathrm{ph}_1$ is the constant one is part of Theorem~C, which uses (QRR).`

### C2. `cor:gau-cycsys` (ii), literal reading; the closing relation in `prop:gau-quintic` (i)  [U03B M4, M5; U04 A/B]

Confirmed: `s03_gauge.tex` l.82: "the same holds with $z\Dq$ replaced by $z\Dq+w$ in \eqref{eq:gau-system} and in
$\PFfive$." Read literally the relations $Y_{a+1}=z\Dq Y_a/\kappa_a$ of (i) stay unshifted, and that statement is false;
the proof (`a2_cyclicPF4.tex` l.122, the automorphism applied to "the identities of (i)") gives the shifted one, which
is how (ii) is used. `s03_gauge_2.tex` l.20: "the closing relation is then given by \cref{cor:gau-cycsys}." (it is an
equivalence, and the PF5 equation for the expansions is known only later).

**C2.1** `s03_gauge.tex` l.82. Old: `with $z\Dq$ replaced by $z\Dq+w$ in \eqref{eq:gau-system} and in $\PFfive$.`
New: `with $z\Dq$ replaced by $z\Dq+w$ in \eqref{eq:gau-system}, in the relations $Y_{a+1}=z\Dq Y_a/\kappa_a$ and in $\PFfive$.`

**C2.2** `s03_gauge_2.tex` l.20. Old: `the closing relation is then given by \cref{cor:gau-cycsys}.` New: `by
\cref{cor:gau-cycsys}~(ii) with $w=\zeta^is^{-1}$ the closing relation
$(z\Dq+\zeta^iL)\Cexp(X_4)=\kappa_4\Cexp(X_0)$ is then equivalent to the equation $\PFfive$, with $z\Dq$ replaced by
$z\Dq+\zeta^iL$, for $\varpi_0\Cexp(X_0)$; for the rows it is obtained in \cref{thm:gau-R0}~(iv) and
\cref{cor:gau-colours}.`
`s03_gauge_2.tex` l.31 (proof). Old: `the closing relation is \cref{cor:gau-cycsys}~(ii) with $w=\zeta^is^{-1}$.` New:
`the statement on the closing relation is \cref{cor:gau-cycsys}~(ii) with $w=\zeta^is^{-1}$ and $(y_0,y_1)=(\varpi_0,\varpi_1)$.`
(No proof uses the closing relation of the expanded actual rows before `thm:gau-R0` (iv), which takes a <= 3.)

### C3. `thm:mod-var` (ii) is proved from [EO (5-15)], declared to enter no proof  [L2 m1]

Confirmed: `s06_model_8.tex` l.39: "This is \cite[Thm.~5.1]{Eynard:2007kz}, (5-14) for $g+k>1$ and (5-15) for $g\geq 1$";
`s01_intro_3.tex` l.23 lists "(5-15)" among the citations that "enter no proof", and (Q) names only (5-14). Part (ii) is
cited nowhere (`thm:mod-var` is cited only in the proof at `s06_model_8.tex` l.85, for (i)). Checked: (ii) follows from
`prop:mod-period` (ii) and (iii), which do not use `thm:mod-var`: d/dv of B_{2h} v^{2-2h}/(2h(2h-2)) is
-(B_{2h}/2h) v^{1-2h}, and the ray of `def:mod-pair` runs from 0 to infinity.

**C3.1** `s06_model_8.tex` l.39, four replacements in the proof.
Old: `This is \cite[Thm.~5.1]{Eynard:2007kz}, (5-14) for $g+k>1$ and (5-15) for $g\geq 1$, for the temperature variation`
New: `Part (i) is \cite[Thm.~5.1]{Eynard:2007kz}, (5-14) for $g+k>1$, for the temperature variation`
Old: `and (5-14), (5-15) read $\delta_\Omega W^{(h)}_n(J)|_{x}=\int_{\partial\Omega}W^{(h)}_{n+1}(J,q)\Lambda(q)$, $\delta_\Omega F^{(h)}=-\int_{\partial\Omega}W^{(h)}_1\Lambda$.`
New: `and (5-14) reads $\delta_\Omega W^{(h)}_n(J)|_{x}=\int_{\partial\Omega}W^{(h)}_{n+1}(J,q)\Lambda(q)$.`
Old: `With $\omega_{h,n}=(-1)^nW^{(h)}_n$ and $\FTR_h=-F^{(h)}$ this gives $\partial_v|_x\omega_{h,n}=(-1)^{n}(-1)^{n+1}(-1)\int_{\mathfrak r}\omega_{h,n+1}=\int_{\mathfrak r}\omega_{h,n+1}$ and $\partial_v\FTR_h=\int_{\mathfrak r}W^{(h)}_1\Lambda=\int_{\mathfrak r}\omega_{h,1}$.`
New: `With $\omega_{h,n}=(-1)^nW^{(h)}_n$ this gives $\partial_v|_x\omega_{h,n}=(-1)^{n}(-1)^{n+1}(-1)\int_{\mathfrak r}\omega_{h,n+1}=\int_{\mathfrak r}\omega_{h,n+1}$.`
Add at the end of the proof: `Part (ii) is not used below. It follows from \cref{prop:mod-period}: by (iii) there, $\partial_v\FTR_h=-\frac{B_{2h}}{2h}v^{1-2h}$, and by (ii) there, $\int_{\mathfrak r}\omega_{h,1}=-\int_\infty^0\omega_{h,1}=-\frac{B_{2h}}{2h}v^{1-2h}$, the ray running from $0$ to $\infty$. (It is also \cite[(5-15)]{Eynard:2007kz}, which is not used.)`

**C3.2** `b_signs_2.tex` l.24. Old: `\cite[Theorem~5.1, (5-14)--(5-15)]{Eynard:2007kz}` New: `\cite[Theorem~5.1, (5-14)]{Eynard:2007kz}`
(`s01_intro_3.tex` l.23 then stays correct. Simpler alternative: delete part (ii) of the theorem.)

### C4-C6. Remarks that rest on computations not in the evidence appendix  [U18_A m7, U08_B m8; the draft notes]

Confirmed: `s10_main_2.tex` l.46 (`rmk:mai-graphlocal`): "computations (\cref{sec:evi}) show that six edge graphs of
genus $4$ become singular while the singular part of the total is unchanged" with a draft note saying the computation
is not in the tables; l.51 (`rmk:mai-necessity`): "Computations (\cref{sec:evi}) give dressings that violate (A) ... and
dressings that are not regular in $u$ in the flat frame" (also a contradiction in terms: a dressing is regular by
`def:tra-dressing`), with the same draft note; `s01_intro_2.tex` l.25: "it is not implied by a statement about the
singular part of the sum (\cref{rmk:mai-graphlocal})." I could not check the computations (none was rerun; they are
not in `sections/c_evidence*.tex`). The texts below are correct whether or not the computations are recorded later;
if they are recorded in an ancillary file, cite it in place of "not recorded".

**C4** `s10_main_2.tex` l.46. Old: `In the normalisation of \cref{rmk:mai-stored}, computations (\cref{sec:evi}) show
that six edge graphs of genus $4$ become singular while the singular part of the total is unchanged. So the singular
part of the total does not determine those of the graph terms, and no argument that uses the total only through its
singular part can prove Theorem~A.` New: `In the normalisation of \cref{rmk:mai-stored}, a computation that accompanied
the work, not recorded in \cref{sec:evi}, found six edge graphs of genus $4$ singular while the singular part of the
total was unchanged. If so, the singular part of the total does not determine those of the graph terms, and no argument
that uses the total only through its singular part can prove Theorem~A; this is not used.` (Then delete the draft note.)

**C5** `s10_main_2.tex` l.51. Old: `Computations (\cref{sec:evi}) give dressings that violate (A), among them
$P^\sharp$ of \cref{rmk:mai-stored}, and dressings that are not regular in $u$ in the flat frame, whose flat sums have
singular edge graphs: neither hypothesis of \cref{thm:red-main} can be dropped in general.` New: `Computations that
accompanied the work, not recorded in \cref{sec:evi}, found singular edge graphs in the flat sums of gauges that violate
(A), among them $P^\sharp$ of \cref{rmk:mai-stored}, and of gauges that are not regular in $u$ in the flat frame (and so
are not dressings in the sense of \cref{def:tra-dressing}). This indicates, but does not prove here, that neither
hypothesis of \cref{thm:red-main} can be dropped in general.` (Then delete the draft note.)

**C6** `s01_intro_2.tex` l.25. Old: `The statement is about each graph separately; it is not implied by a statement
about the singular part of the sum (\cref{rmk:mai-graphlocal}).` New: `The statement is about each graph separately; we
do not know how to obtain it from a statement about the singular part of the sum (\cref{rmk:mai-graphlocal}).`

---------------------------------------------------------------------------------------------------------------------

## D. Other wording changes (E8, and the companion points of E1, E4, E7)

All confirmed in the text; E8's readings (the bivector, reading A of the hatted unit, the expanded TFT) were re-derived
for D1 (matrix of the inverse bivector in lowered components is eta; N = 1 - R'(-psi)^T eta^{-1} R'(-psi') for
symplectic R') and read for D2-D5. None changes a proof.

**D1** (E8-1) `s02_setup_8.tex` l.46 (`def:set-outer`). Old: `where $1=\sum_i\tilde e_i\otimes\tilde e_i$ is the inverse
bivector of the colour space,` New: `where $1=\sum_i\tilde e_i\otimes\tilde e_i$ is the inverse bivector of the colour
space and $\eta^{-1}\in\bbH\otimes\bbH$ is the inverse bivector of $(\bbH,\eta)$ (in lowered components its matrix is
$\eta$, whereas the matrix $\eta^{-1}=\eta/25$ is the pairing on lowered components; so for symplectic $R'$, in lowered
components, $N(\psi,\psi')=1-R'(-\psi)^T\eta^{-1}R'(-\psi')$),`
`s07_comparison_2.tex` l.78 (`thm:flh-main`, Step 1). Old: `With $N(\psi,\psi')\coloneqq1-\bigl(R'^{-1}(\psi)\otimes R'^{-1}(\psi')\bigr)\eta^{-1}$, the edge`
New: `With $N(\psi,\psi')\coloneqq1-\bigl(R'^{-1}(\psi)\otimes R'^{-1}(\psi')\bigr)\eta^{-1}$, $\eta^{-1}$ being the inverse bivector of $(\bbH,\eta)$ as in \cref{def:set-outer}, the edge`

**D2** (E8-2) `b_signs.tex` l.32. Old: `$\hat\Delta_i^{1/2}=\zeta^{4i}\hat x$ with $\hat x=\varpi_0$,` New: `the hatted
TFT $\hat\TFT_{g,n}(\tilde e_i,\dots,\tilde e_i)=(\zeta^{4i}\hat x)^{2g-2+n}$ with $\hat x=\varpi_0$,` (the letter
$\hat\Delta_i^{1/2}$ is defined as $\zeta^{4i}\varsigma$ for the model TFT in `def:set-modcohft`).

**D3** (E8-3) `s02_setup_8.tex` l.11 (`lem:set-dilaton` (iii)). After `$y$ being a formal variable.` add: `The left
sides below are defined coefficientwise in $y$: writing $q/\varkappa+z1=y/\varkappa+\epsilon\,z1$ with
$\epsilon\coloneqq1-\varkappa^{-1}$, every monomial of degree $a$ in $y$ receives a series in $\epsilon$, which converges
in $\mathsf R$; for $y\in z^2V\fps z$ only finitely many $a$ contribute, and the result is the Taylor expansion of
\cref{def:flh-outer}~(ii).`

**D4** (E8-4) `s07_comparison_2.tex` l.9 (before the proof of `lem:flh-hat`). Old: `$\hat R'(\zh)\coloneqq\hat\Qr(\zh)\hat M(\zh)^{-1}=\Lambda_sR'(s^4\zh)$, the hatted term is`
New: `$\hat R'(\zh)\coloneqq\hat\Qr(\zh)\hat M(\zh)^{-1}=\Lambda_sR'(s^4\zh)$, which is the map $R'(s^4\zh)$ written in the lowered components with respect to the $\hat H^a$ (in these, $\hat1=s1$ has the components $5\,\delta_{a3}$; $\Lambda_s$ is not applied to vectors as an isometry of $\bbH$), the hatted term is`

**D5** (E8-5) The TFT in Theorem A is the expanded one. Three places.
`s10_main.tex` l.11. Old: `$\TFT$ the TFT of \cref{def:set-frob} (\cref{def:set-outer}).` New: `$\TFT$ the TFT of \cref{def:set-frob} expanded in the frame, $\TFT_{g,n}(\tilde e_i,\dots,\tilde e_i)=(\zeta^{4i}s\varpi_0)^{2g-2+n}$ (\cref{def:set-outer}).`
`s01_intro_2.tex` l.15. Old: `where $\TFT$ is the topological field theory of the formal quintic (\cref{def:set-frob}).` New: `where $\TFT$ is the topological field theory of the formal quintic (\cref{def:set-frob}), expanded in the frame.`
`s02_setup_8.tex` l.46. Old: `the classes of the cohomological field theory $M.\TFT$ on the colour space \cite[Def.~2.13]{PPZ15},` New: `the classes of the cohomological field theory $M.\TFT$ on the colour space \cite[Def.~2.13]{PPZ15}, where $\TFT$ is the TFT of \cref{def:set-frob} expanded in the frame, $\TFT_{g,n}(\tilde e_i,\dots,\tilde e_i)=(\zeta^{4i}s\varpi_0)^{2g-2+n}$,`

**D6** (E1, U14-5) `s08_connection_5.tex` l.4 (proof of Theorem P(a) (ii)). After `computed over $A$.` add: `As in the
proof of (i), $\mathsf S_0$ stands for the formal phase $\mathsf S_0(L+t)-\mathsf S_0(L)$ and the amplitude is $a_0$, both
with Taylor data in $\AL$, mapped to $A$ by $\jmath$.` (So `lem:fsp-cov` in Step 2 is applied to a phase over A.)

**D7** (E1, U14-2/3) `s08_connection_5.tex` l.10 (Step 3). Old: `and expansion of $(1-s^5X^5)^{-1}=\sum_Pu^PX^{5P}$ in
$\bbQ[X]\fps{u}$ followed by expansion at $X=1$ coefficientwise in $s$.` New: `and expansion of
$(1-s^5X^5)^{-1}=\sum_Ps^{5P}X^{5P}$ in $\bbQ[X,X^{-1}]\lser s$ followed by expansion at $X=1$ coefficientwise in $s$, a
ring homomorphism $\bbQ[X,X^{-1}]\lser s\to A\fps{X-1}$ which commutes with $\dd/\dd X$.` (The first ring does not
contain $s^{-1}$ and $X^{-1}$, which are generators of the source.)
`s08_connection_5.tex` l.16 (Step 4). Old: `these identities hold by the two-homomorphism argument of Step 3.` New:
`these identities hold by the two-homomorphism argument of Step~3, applied to the sums over $j$, which lie in $\bbQ[s,X,(1-s^5X^5)^{-1}]$.`

**D8** (U18_A m9) The two lists of [PPZ15] should agree. `s10_main_2.tex` l.36, old: `Def.~0.3; (5); Rem.~2.1;` new:
`Def.~0.3; (5); Rem.~2.1 and (6);` and `s01_intro_3.tex` l.9, old: `Def.~0.3 (iii);` new: `Def.~0.3, in particular (iii);`
(with B4.4 and B4.5 the lists then differ only by the items not needed for Corollary B).

**D9** (E7 W5; optional, not recommended) `s10_main.tex` l.81: "Whether $\ev$ is injective on $\Nring$ is neither
needed nor proved." E7 and X4 argue that ev is injective (the differential Galois group of PF4 is Sp_4, by
Beukers-Heckman and Schlesinger). The sentence is accurate as it stands. Adding "(it holds, by ...)" would bring in two
external theorems and a new bibliography entry for a fact that no proof uses; if the author wants it, it belongs in a
remark marked as not used, and `s02_setup_8.tex` l.52 ("the injectivity of $\ev$ is not claimed") stays as it is.

Not taken over from the judges: E1's remark on "the terms in $\log z$ cancel before and after the shift" (the sentence
disappears with B2.2); E5's optional pin for (8-1) (the citation becomes orientation only, A2.4).

---------------------------------------------------------------------------------------------------------------------

## Consequence: the declared inputs after the amendments

The list of inputs has the same five letters. Two grow in wording to say what was already used ((QRR), (T)); (Q)
loses five items and gains two pins. No new source is added. Exact new statements:

**(Q)** Quoted theorems.
* [Eynard-Orantin]: unchanged (Def. 4.2, (4-9) with the kernel (4-3); Theorems 4.2, 4.6; Theorem 5.1, (5-14), with
  Section 5.3, (5-27)-(5-30)). After C3, (5-15) is used nowhere, as l.23 of `s01_intro_3.tex` already says.
* [DOSS] Theorem 4.1 with (4-1)-(4-3), its graph weights (2-22)-(2-26) (dilaton leaves (2-23), factor 1/|Aut| in (2-26)),
  as recalled in [DNOPS] Section 2.6, (2-30)-(2-32), (2-34), (2-35), matrices in the left-action convention of [DNOPS]
  Remark 2.1, the action being that of [PPZ15] Defs. 2.2, 2.13; applied to x = 4 W_v, y = p. Its three hypotheses
  (symplecticity, (2-32), the edge condition (4-1)) are proved in the note (`prop:mod-Rtr` (i)-(iii)).
  Removed: [DNOPS] (2-33) (misprinted), Prop. 4.3, the sentence before Theorem 8.1, and (8-1) as an input.
* [ABDKS25] Theorem 7.1: unchanged.
* [PPZ15]: Section 0.2 (1); Def. 0.3; (5); Remark 2.1 and (6); Def. 2.2 and Def. 2.5 with the conventions stated after
  them; Prop. 2.12, with the formula (8) in its proof; Def. 2.13; Prop. 2.14. Props. 2.12 and 2.14 are used as statements
  over Q and transferred to the rings Q(zeta)[vs, vs^{-1}] and Q(zeta)((s)) by a density argument given in
  `def:set-graphsum`. Removed: Props. 2.4, 2.8, 2.9 (no proof of [PPZ15] is re-used).

**(QRR)** The R-matrix along the line: flatness and the value at the large-radius point. The R-matrix R of (T) has
coefficients in End(H) (x) Q(zeta)[[q]] and (flatness) its columns R e~_i e^{u_i/z} solve the quantum differential
equation z d_T v = H * v [GJR18, Lemma 5.2, (16); GJR17, Section 4.2, (12)-(13); Giv01, Section 1, Proposition (a)-(c);
Teleman, Props. 7.7, 8.2 (ii), 8.4, 8.7 (ii)]; (value) R(q = 0) = diag(h_{w_0}, ..., h_{w_4}) in the normalised canonical
basis at q = 0 [GJR18, (17); GJR17, (19)-(20); Giv01, Section 2.2, Theorem 1 and Remark; Coates-Givental, Theorem 1,
Corollary 1]. New relative to the current text: flatness, which was called a convention. With the optional lemma A1.16
these two properties determine R and imply that it is symplectic.

**(M)** unchanged. **(G0)** unchanged (genus zero only).

**(T)** Teleman's reconstruction, with the axioms of the twisted theory in all genera. (i) At lambda_i = zeta^i the
classes Omega_{g,n} of `def:set-fq`, in all genera, form a cohomological field theory with flat unit for eta
[GJR18, Section 4.1] and satisfy the divisor equation in all genera, so that the theory at Q = e^{T(q)} is the theory
at t = T(q)H with Q = 1. (ii) This (semisimple) theory equals R.omega for an R-matrix R = 1 + O(z) over Q(zeta)[[q]],
symplectic [Teleman, Theorem (2) in (1.15); GJR18, Section 4.1, Theorem 4.1]; this R is the R-matrix of `def:set-rows`.
New relative to the current text: the all-genus axioms and the all-genus divisor equation in (i).

**[AC96, (1.7)]** unchanged.

Conditionality, as it should now read (A1.4): "Theorem A holds modulo the quoted theorems (Q) of [EO], [DOSS] (as
recalled in [DNOPS]), [ABDKS25] and [PPZ15] with the statement numbers above, the flatness of the R-matrix and its
value at the large-radius point (QRR), the genus-0 mirror theorem (M) and the genus-0 axioms (G0) of the twisted
theory, Teleman's reconstruction with the axioms of the twisted theory in all genera (T), and the formula
[AC96, (1.7)]; every other step is proved in the note."

Which input is used where, after the amendments.
* Theorem C (`thm:int-D`) and the proof of Theorem A use of the R-matrix only: R(0) = 1, the ring Q(zeta)[[q]],
  flatness and the value at q = 0, that is (QRR), together with (M), (G0); Theorem A also uses (Q) and [AC96]. In both,
  (T) only names the object ("the R-matrix of the formal quintic"). With A1.16 that object is determined by (QRR).
* Corollary B uses in addition (T) itself (ev(F_g^N) = F_g), including its part (i).
* Theorem Gamma: unchanged (no quoted input, except [DLMF 24.4.1] in part (b)).
* What is proved, in words (A3): the gap for the conifold expansion of the graph-sum representative F_g^N, a formal
  change of frame; nothing about the analytic continuation of F_g. With A3.4 (a) this expansion is also the conifold
  expansion of the function F_g in the sense recalled from [CGYZ26, Def. 1.2].

Counts: A 3 items (A1: 15 edits + 1 optional lemma; A2: 9 edits; A3: 4 edits, one with two versions); B 5 items
(B1: 3 edits; B2: 2; B3: 7; B4: 5; B5: 1 + 1 optional); C 6 items (8 edits); D 9 items (8 with edits, D9 optional and
not recommended).

Not verified for this list: the statements of the sources beyond the spot checks named at A1 and A2 (in particular
[GJR17], [Giv01], [Coates-Givental], [PPZ15], [CGYZ26] were not re-read); the unrecorded computations behind C4-C6;
E7's proof that ev is injective. Nothing was compiled; the new LaTeX uses only macros of `macros.tex` and commands
already used in the sections quoted (`\coloneqq`, `\fps`, `\lser`, `\coef`, `\de`, `\dd`, `\re`, `\cref`).
