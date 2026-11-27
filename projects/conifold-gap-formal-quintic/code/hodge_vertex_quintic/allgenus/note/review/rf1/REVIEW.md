# Review rf1: `s01_intro*.tex`, `s10_main*.tex` (part 1 of 2)

Part 2 (`REVIEW_2.md`): the chain of Theorem A step by step, Corollary B and its constant, the conditionality, what remains
open, and every change that affects other sections.

## 0. Verdict

| section | verdict | one line |
|---|---|---|
| `s01_intro` (3 parts) | **sound-after-fixes** | statements agree with the body; inputs list corrected (PPZ pins, where the inputs enter, Theorem C inputs); the "not used" list corrected (ZZ08) and made exact as a statement; two overclaims removed |
| `s10_main` (4 parts) | **sound-after-fixes** | assembly of Theorem A and Corollary B correct and complete; two definitional/ordering gaps closed (forward use of `lem:mai-Psi` in a proof; non-symmetric edge in `def:mai-Fcf`); one wrong remark corrected (what tests `D = 1`); several scope and notation slips fixed |

No `\draftnote` was added on a step of the chain of Theorem A (none was needed). One `\draftnote` was added outside the chain,
in `rmk:mai-analogy` (a computation cited there is not recorded in appendix C; s.3, S10-15). Compile (scratch folder, current
files of all sections, pdflatex, bibtex, pdflatex x2): no error, no undefined or multiply defined reference, no overfull box in
`s01`/`s10`; 159 pages (the page count also reflects concurrent edits of other sections).

## 1. What was checked against what

* Every statement of s01 against its body version: `thm:int-main` / `thm:mai-main`, `def:set-N`, `eq:set-cg`, `rmk:set-cg`;
  `cor:int-gap` / `cor:mai-gap`, `def:mai-Fcf`; `thm:int-D` / `thm:gam-D`, `cor:mai-actual`; `thm:int-gamma` / `thm:gam-main`
  (i), (ii), (iv), `thm:gam-ginf`, `def:gam-c`; `rmk:int-constant` / `rmk:mai-constant`.
* Every `\cref` of s01 and s10 (138 labels, script over the `.aux`): existence, and for proofs that the target is earlier, an
  appendix lemma or a citation. Every citation key and pin against BIBPLAN, DIGEST_L(L1-L48), DIGEST_L_4 s.2 and the body
  (script over all `\cite[...]` of the note).
* The literature paragraphs of `sec:int-gap` against DIGEST_L_4 s.2.1-s.2.3 and L1-L19 (CGYZ26 Thm 2 only for the quintic
  threefold; CGYZ26 Rem. 1.9 for the formal quintic; GJR18 Rem. 7.2 "analogous", stated without proof; FCCZ26 not used; no
  withdrawn preprint; no priority or novelty claim beyond L43-L48).
* The inputs and the conditionality against OUTLINE s.0.2 ("Exact final list"), GAPS_SYNTHESIS_2 s.3, DECISIONS (D/Q items),
  and against what the body actually uses (pins scanned in every section file).
* s10 against the sources of record: MS5/ST2 step 8 (assembly), G5_2 s.2 with G5R s.2 (Corollary B, `rmk:mai-constant`), G5_3
  s.3 with G5R_2 s.3 (`lem:mai-Psi`, `lem:mai-Ghat`, `cor:mai-natfact`), DC1 UNI-1..UNI-8, DC2 s.2.6 (universality), DA A9,
  A12, A19-A20, A25 (remarks), OUTLINE_8 (s10 plan).
* Arithmetic redone by hand: `c_2..c_5`; `(-1)^(g-1) B_2g = |B_2g|`; `2F1(1,4/5;9/5;u) = Sum 4u^n/(5n+4)`; `z/varsigma^4 =
  tau zhat`; `Psi_Q`, `Psi(u)` entries and orthonormality; `mu(b,k)` and the zero set `{k = 0} U {(1,1)}`; the scaling
  endomorphism of `rmk:mai-constant` on `c_j`, `i_j`, `W`, `C_1` and on the rows; the Galois descent; `lem:mai-pole`,
  `lem:mai-coord` (a)-(e); the residue/`lambda` relation of `cor:mai-univ`; the critical values `5 p_c^4` of `4 W_v`.
* Web (public, generic): arXiv:1303.1043v3 (PPZ15) abstract page and HTML: Prop. 2.14 = "the unit-preserving R-matrix action is
  a left group action"; its proof uses Def. 2.13, Prop. 2.9 and Prop. 2.8 (used for S01-4, S10-3).

## 2. Issues in s01 (location; problem; fix)

* **S01-1** (`s01_intro_2.tex`, paragraph after Theorem A). "it is stronger than any statement about their sum" is false
  (Theorem A does not determine the regular part of `F_g^cf`), and "no statement about the total potential alone can imply"
  Theorem A goes beyond the computation cited (sec:evi-univ; DC1 UNI-7, MECHANISM s.2: the *polar part* of the total does
  not move). Fixed: "With (T) it implies the statement about their sum, Corollary B; it is not implied by a statement about
  the singular part of the sum: ..." Also "the hatted outer sum of `def:set-hat`" -> formed with the data of `def:set-hat`
  (`def:set-outer`), where `G_Gamma` is defined.
* **S01-2** (`s01_intro_2.tex`, after Theorem Gamma). "formed by parts (ii) and (iv)": the intro statement also contains part
  (i) (`c_n` in `A_L`) and only the first assertion of (iv). Fixed. The non-claim "we do not claim that the series of GJR18
  coincides with `c(z;L)`" cited `cor:gam-free`, which does not contain it (it is the sentence after `thm:gam-ginf`). Fixed:
  `\cref{thm:gam-ginf}`; `cor:gam-free` now cited for what it proves (the free coefficients), "not used in the proof of
  Theorem A".
* **S01-3** (`s01_intro_3.tex`, "Where the inputs enter"). "the comparison theorem of sec:flh uses [PPZ15, Defs. 2.2, 2.5,
  Props. 2.4, 2.9]": the s07 text uses Def. 0.3, Rem. 2.1 and (6), Defs. 2.2, 2.5, 2.13, Prop. 2.12 and the dilaton lemma, and
  neither Prop. 2.4 nor Prop. 2.9 (2.4 appears only in `rmk:flh-fail` (iii), as *not* needed; 2.9 nowhere). Fixed. The list
  of places where PPZ15 enters missed `sec:mod` (`thm:mod-tr`, PPZ Defs. 2.2, 2.5, 2.13, Prop. 2.12) and `def:set-outer`.
  Fixed.
* **S01-4** (`s01_intro_3.tex`, (Q), PPZ pins). Corollary B transfers PPZ's *proofs* to `Q(zeta)((s))`; PPZ's proof of
  Prop. 2.14 uses Props. 2.8 and 2.9 (checked in the v3 HTML), so Prop. 2.8 belongs to the pins as much as 2.9. Added
  "Props. 2.8 and 2.9" and "Props. 2.4, 2.8 and 2.9 enter only through the proof of Prop. 2.14". This deviates from the
  printed pin list of OUTLINE s.0.2 at the level of pins only; the conditionality ("PPZ15 quoted") is unchanged (to ratify).
* **S01-5** (`s01_intro_3.tex`, "In particular ..."). (a) Corollary B's PPZ inputs listed only Def. 2.13, Prop. 2.14; the
  proof also uses Prop. 2.12 (unit of `M.Cexp(TFT)`) and Rem. 2.1, Def. 2.2 (the edge). Fixed, consistent with s10. (b)
  "Theorem C uses (QRR), through the glue": `D = 1` for all colours needs `D_i(z) = D_0(zeta^-i z)` (`cor:gau-colours`), whose
  proof uses the colour covariance of the actual rows (`cor:gau-rowsNN` (iii), `lem:gau-colrows` Step 3, i.e. (QRR) for
  `i != 0`); and the rows rest on (M), (G0). Fixed; `cor:mai-actual` added.
* **S01-6** (`s01_intro_3.tex`, "not used" list). `[ZZ08]` was said to corroborate `prop:set-cyclicLR`; in the body ZZ08 is
  cited only for the specialisation (definitional, s02) and a convention (s.B), and the corroboration sentence
  (`s02_setup_3`) cites GR19 and LP20 s.2.4.3 only. Removed; LP20 pinned as in the body. Added the citations used only in
  remarks (Mumford's relation, LP20 (26), the Teleman s.8 statements of `rmk:set-qrrderiv`/`rmk:set-qrrsign`) and Givental
  Cor. 11.7 (misprinted form). Added one sentence classifying every other citation (further sources of inputs, comparisons,
  conventions, context; none enters a proof), so that the paragraph is exact as a statement; it now speaks of Theorem A
  and Corollary B. (Scan: EO Thm 4.1 (4-20), App. A (A-1),(A-5) and DNOPS (2-26)-(2-29), Prop. 5.2 are comparisons or
  conventions in s06; `lem:mod-genus0` (ii) is proved from the recursion, not from EO Thm 4.1.)
* **S01-7** (`s01_intro_3.tex`, What is not claimed (2)). "Nothing is claimed for frames other than those of
  `def:set-frames`" contradicted `rmk:mai-constant`, which proves the gap for the rescaled pairs `(k_0 varpi_0, k_1 varpi_1)`,
  `k_i` in `C^x` (G5_2 s.2.5). Fixed ("and the rescalings of their two periods treated in `rmk:mai-constant`").
* **S01-8** (`s01_intro.tex`, Hodge frames). "Frames of Hodge type do not suffice here:" was unconditional, the evidence being
  genus 2 only (REV_2 M6). Fixed: "at least at genus 2".
* **S01-9** (`s01_intro.tex`, TR paragraph). DOSS Thm 4.1 is applied to `4 W_v` and transferred by scaling (rule 8, s06
  opening), not to `W_v`; added "applied to `4W_v`, which is the same curve up to scaling".
* **S01-10** (`s01_intro.tex`, Givental s.10). The clause "after which the point `q = 0` is a regular point of the phase" is
  not in DIGEST_L_4 s.2.2(4); aligned with it and with `rmk:gam-prior`: "a change of chart in which the critical point
  survives at `q = 0` and the integrand there factorises into Euler integrals".
* **S01-11** (`s01_intro_2.tex`, `s01_intro_3.tex`). "unique up to a constant series in `z`" -> "unique up to a scalar factor
  `c(z)` independent of `u`" (`thm:gau-U` (b): `P -> P c(z)`, `c` a unit of `k[[z]]`).
* **S01-12** (`s01_intro.tex`). "in every genus" -> "in every genus `g >= 2`".

Checks of s01 that passed without change: the statements of Theorems A, C, Gamma and Corollary B agree with the body (the body
now also uses "`T_Gamma` in `c_g t_c^(2-2g) + Q[[u]]`", S10-16); the values `c_2..c_5`; `rmk:int-constant` against L8, G5_2
s.2.6 and G5R F2.4 (Ex. 1.3; s.4.2, Prop.-Def. 4.5 (1); Rem. 1.4 with rational `c_i`; Conj. 1.5 constant; Def. 1.8 sign not
imported); the (QRR), (M), (G0), (T), [AC96] items and the definitional citations against OUTLINE s.0.2 and GS_2 s.3; the
conditionality sentence (verbatim OUTLINE s.0.2, with [AC96, (1.7)]); the rule-4 `\draftnote`s at the first Teleman and
Givental 1996 citations; the D11 `\draftnote`; the CGYZ26/Br25 version `\draftnote`s; no project bookkeeping in the text.

## 3. Issues in s10 (location; problem; fix)

* **S10-1** (order; `def:mai-Fcf`, proof of `cor:mai-gap`). Both used `Psi_Q^T eta^-1 Psi_Q = 1` and the explicit `Psi_Q` of
  `lem:mai-Psi`, stated and proved later (sec:mai-nat): a forward dependence of a proof. Fixed: `lem:mai-Psi` moved before
  `def:mai-Fcf` (same label; number 10.5 -> 10.2); part (ii) now for `Psi_cf := Cexp(Psi_Q)` and `Psi(u) := Lambda_s Psi_cf`
  (both orthonormality relations proved); sec:mai-nat identifies `Qhat(u, zhat = 0) = Psi(u)` and defines `R~(zhat; u)`.
* **S10-2** (`def:mai-Fcf`, the edge). `E_X = (N_X(psi,psi') - N_X(-psi',psi'))/(psi+psi')` is symmetric only if `X` is
  symplectic, and in `NN (x) Q(zeta)` the symplecticity of `R~` is not known (G5R F2.3; `rmk:con-genus`; the paragraph after
  `lem:flh-inputs`), while `Cont_Gamma` is defined for symmetric bivectors only (`def:set-graphsum`, `def:flh-outer` (i)):
  `F_g^NN` depended on edge orientations. Fixed with the symmetric exact-division form `(N - N(-psi',psi')/2 -
  N(psi,-psi)/2)/(psi+psi')`, equal to PPZ's edge for symplectic `X` (the same form as the current `def:set-outer`); the
  `ev` property now says why `ev(E_R~)` is PPZ's edge (R symplectic over `Q[[q]]`).
* **S10-3** (proof of `cor:mai-gap`). Local `Psi := Cexp(Psi_Q)` clashed with `Psi(u)` (NOTATION_2) -> `Psi_cf`. "which is the
  case `A = Psi` of `lem:flh-frame`" misattributed (that lemma pushes a CohFT forward; here only the edge/leaf algebra) ->
  "by the computation of the proof of `lem:flh-frame`". Added why the edge commutes with ring maps (exact division). The
  coefficient-ring remark now names Props. 2.8, 2.9 (used by PPZ's proof of 2.14), as in G5_2 s.2.4. Long inline formulas
  turned into displays (removes two overfull lines, one pre-existing).
* **S10-4** (paragraph after `cor:mai-natfact`). "only the symplecticity of `G` tests the full identity `D = 1`" is wrong:
  (a), (b), polynomiality are invariant under `D_0 -> D_0 c(z)`, and the factor `R~ M_0(tau zhat)^-1 = G c(z)^-1` is
  symplectic iff `c(z)c(-z) = 1`: symplecticity sees only the even part of `log D_0`, a factor `exp(a)`, `a` odd, is
  invisible (`rmk:con-why`, `rmk:mai-stored`, `lem:flh-inputs` (iii)). Rewritten; the computation (shifted free coefficient,
  rows not symplectic) is described as such. The final text does not use the withdrawn claim "R' symplectic without `F_g`".
* **S10-5** (`rmk:mai-stored`). "The values of `m_5` and `m_10` ... are in sec:evi-univ": `m_10` is recorded nowhere, `m_5`
  in the opening of sec:evi. Fixed.
* **S10-6** (`rmk:mai-graphlocal`; sentence after Theorem A). "no argument about total potentials, such as the gap of the
  total or an integrable hierarchy for the total, can prove Theorem A" exceeds the computation (only the singular part of
  the total is unchanged). Narrowed to arguments that use the total only through its singular part.
* **S10-7** (sec:mai-univ opening). "The proof of Theorem A uses the formal quintic only through its cyclic system" is false
  for the whole proof (comparison theorem, `D = 1`); restricted to its flat half (`cor:red-can` with `thm:mod-MLK`).
* **S10-8** (sec:mai-univ opening). `A~_y` was defined only for `y'(0) != 0`, but `lem:mai-pole` uses `y` of any order. Fixed.
* **S10-9** (`lem:mai-coord` (d)). "a dressing `P` of `X`" is undefined -> "of the model at `v_X(u)`", with `P(phi(w),z)` a
  dressing at `v_q(w)`.
* **S10-10** (after `lem:mai-coord`). "systems with other products are coordinate changes ..." holds only for products
  `-v_X'^5/v_X`, `v_X` in `uQ[[u]]` (residue minus a rational fourth power). Fixed.
* **S10-11** (`cor:mai-univ`). "the canonical gauge from `A~_(v_A)` to `A`" was not defined for `A` other than the quintic
  system (`def:gau-can` is stated for it only). Defined in the statement (existence, uniqueness, grading, rationality by
  `thm:gau-U`, `prop:gau-graded`); in the proof, `Y^C(v_A(u))` identified with the vanishing solution of `A~_(v_A)`.
* **S10-12** (`rmk:mai-constant`). `2g-2+L` used `L` (global `u^(-1/5)`) for the number of leaves -> `Sum_v l_v`; the
  endomorphism `\mathcal S` clashed with the solution space `\mathcal S` of s03 and appendix A -> local `\mathfrak k`.
* **S10-13** (`rmk:mai-apriori`). `L` as a count -> `m`; `G_1..G_(2g-3)` (clash with `G_Gamma`) -> `[u^k]G`; added that the
  vertex classes modulo `u^(2g-2)` also enter.
* **S10-14** (`rmk:mai-tau`). "`vartheta = tau`, the canonical coordinate": `tau` is the rescaling given by the canonical
  coordinates. "known only from the computations of sec:evi" for every `g`: appendix C records it at genus 5, modulo `u^8`
  (sec:evi-graphs). Fixed.
* **S10-15** (`rmk:mai-analogy`). "z divided by the canonical coordinate" holds up to the factor `-(5/4) zeta^i` (fixed); "the
  computations of sec:evi show that through genus 3 the graph sums alone do not single it out" points to a computation that
  appendix C does not record (DIGEST_A A19, genus-3 classification). `\draftnote` added (record it or drop the sentence).
* **S10-16** (statements). Theorem A (ii), Corollary B, `cor:mai-univ` (ii), `rmk:mai-annih`: "`= c t^(2-2g) + Q[[u]]`" ->
  "`in`", as in s01.
* **S10-17** (proof of Theorem A). A paragraph naming the main dependencies of the cited statements was added, so that the
  chain C1-C13 of OUTLINE s.0.2 is visible in the assembly (s.4 of part 2).
* **S10-18** ("Where the inputs enter", after Theorem A). Made consistent with s01: (M), (G0) also through `thm:set-mirror`,
  `cor:set-cyclicPF5`; PPZ also through `def:set-modcohft`, `def:set-outer`, `thm:mod-tr`; DOSS "as recalled in" DNOPS.
* **S10-19** (`cor:mai-natfact` (c), proof). `Psi^-T Psi^-1` -> `Psi(u)^-T Psi(u)^-1`; the symplecticity of the expanded rows
  used there is noted as the case `D = 1` of `lem:flh-inputs` (iii).
