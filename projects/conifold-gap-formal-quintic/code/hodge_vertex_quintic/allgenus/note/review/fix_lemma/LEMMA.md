# Fix of review item I-1: the lemma `lem:set-pushins` (key fix/lemma)

Author-facing record; not part of the note. Object: review item I-1 (`review/REVIEW_SYNTHESIS.md` s.0, s.3;
`REVIEW_SYNTHESIS_2.md` s.7 I-1, s.8 item 10, s.10 recommendation 1). Files edited (only these):
`sections/s02_setup_4.tex`, `sections/s02_setup.tex`, `sections/s06_model_2.tex`, `sections/s06_model_3.tex`,
`sections/s10_main.tex`, `sections/s10_main_2.tex`, `sections/s01_intro_3.tex` (input lists only). `refs.bib` unchanged.

## 1. What was done

* **New lemma** `lem:set-pushins` ("The translation: push-forward and insertion"), Lemma 2.22 of the build, in s02.3
  (`s02_setup_4.tex`) directly after `def:set-graphsum` and before `thm:set-GT`. Statement and complete proof (about 1.5 pp):
  - (i) *class identity.* On `Mbar_(h,n+m)` (2h-2+n > 0), `p_m` forgetting the last m points, `k_j >= 0`, `c_s >= 1`:
    `prod_j psi_j^(k_j) prod_s psi_(n+s)^(c_s) = prod_j (p_m^* psi_j)^(k_j) prod_s psi_(n+s)^(c_s)` in `H^*(Mbar_(h,n+m))`.
  - (ii) *vertex.* For a CohFT `Omega` [PPZ15, Def. 0.3] and `t` in `z^2 V[[z]]`, the translation `t Omega` of [PPZ15, Def. 2.5]
    times `prod psi_j^(k_j)` equals `sum_m (1/m!) p_(m*)( Omega_(h,n+m)(a, t(psi), .., t(psi)) prod psi_j^(k_j) )` with the classes
    `psi_j` of `Mbar_(h,n+m)` (class level); integrating, `int (t Omega)_(h,n)(a) prod psi_j^(k_j) = sum_m (1/m!) <psi^k a ... t^m>_h`,
    i.e. `d^n Fbar_h(q)` at `q = -z1 + t` when `1` is a unit. Finite sums.
  - (iii) *graph.* For `Omega` with unit and `R = 1 + O(z)` symplectic, `g >= 2`, every leg-free stable graph: `Cont_Gamma[E_R, q_R;
    Omega]` is the term of `Gamma` in `int_(Mbar_g) (R.Omega)_(g,0)`, `R.Omega` of [PPZ15, Def. 2.13]; hence the integral is the graph
    sum of genus g; by [PPZ15, Prop. 2.12] `R.Omega` again has the unit `1`, so (iii) applies to it.
* **Proof** (all steps in sentences, as in the draftnote's argument and the review's suggested route, "induction on the number
  of forgotten points"):
  - (i), m = 1: [PPZ15, (8)] (the formula for the map forgetting one point, `pi^*(psi_j^d) = psi_j^d - delta_j pi^*(psi_j^(d-1))`,
    `delta_j` the divisor with a genus-0 component carrying only the markings j and N+1) gives `pi^* psi_j = psi_j - delta_j`
    (d = 1); comparing d = 2 with the square of d = 1 gives `psi_j delta_j = 0`; the relabelling automorphism exchanging j and
    N+1 gives `psi_(N+1) delta_j = 0`; so `psi_j psi_(N+1) = pi^*(psi_j) psi_(N+1)`, and (*) for monomials. Induction on m:
    `p_m = p_(m-1) o pi`, (*) twice and (i) for m-1.
  - (ii): Def. 2.5 with the convention stated after it, projection formula, (i) (every translation point carries `psi^c`,
    `c >= 2`), `int p_* W = int W`; finiteness by degree (`2m` against `3h-3+n+m`).
  - (iii): Def. 2.13 (R.Omega = action of Def. 2.2 on the translation by `t_R`); the contribution of Def. 2.2 is the push-forward
    along the gluing map [PPZ15, s.0.2, (1)]; the edge of Def. 2.2 (iii) is `E_R` (inverse R, [PPZ15, p. 20]); expansion of the
    edges into monomials, Fubini on `prod_v Mbar_(g_v,n_v)`, (ii) at every vertex (`q_R + z1 = t_R`).
  - The hypothesis actually needed in (i) is `c_s >= 1`; the translation has `c_s >= 2` ([PPZ15, Def. 2.5]).
* **Draftnote deleted** in `def:set-graphsum`; the sentence is now: the PPZ translation pushes the points forward (classes `psi`
  of `Mbar_(g_v,n_v)`), the vertex tensors insert them (classes of `Mbar_(g_v,n_v+m)`); "The two forms agree: for g >= 2,
  `Cont_Gamma[E_R,q_R;Omega]` is the term of `Gamma` in `int (R.Omega)_(g,0)` (`lem:set-pushins` (iii) below)". The definition of
  `R.Omega` (legs, edges, translation; PPZ p. 20 convention) is kept. Draftnotes in the `.tex` files: 17 -> 16.
* **Citations of the lemma added**:
  - `thm:set-GT` "In particular": "since `F_g = int Omega_(g,0)` (`def:set-fq`), `lem:set-pushins` (iii) gives ..."; the paragraph
    after the theorem says the second assertion follows from the first by (iii) (TFT with unit 1, `R(0) = 1`, `def:set-rows` (a)).
  - third property after `def:mai-Fcf` (`s10_main.tex` l.83): the graph sum with this R-matrix is the integral of the genus-g
    class of its unit-preserving action on `TFT` (`lem:set-pushins` (iii), unit `1_TFT`), and this integral is the potential by (T).
  - `cor:mai-gap` Step 1 (`s10_main_2.tex` l.25): `Cexp(TFT)` is a TFT over `Q(zeta)((s))` with unit `Cexp(1_TFT)` (the TFT axioms
    are polynomial identities in the `Delta_i^(+-1/2)`; this sentence is new, it was presupposed), then (iii) gives
    `int (Cexp(R~).Cexp(TFT))_(g,0)`.
  - `cor:mai-gap` Step 2 (l.31): the right-hand side by (iii) applied to `M.Cexp(TFT)` (unit by [PPZ15, Prop. 2.12]) and to the
    symplectic `X`, `X(0) = 1`; closing sentence (l.40): PPZ list now `s.0.2, (1); Rem. 2.1; Defs. 2.2 and 2.5; (8); Prop. 2.12;
    Def. 2.13; Prop. 2.14`, and `lem:set-pushins` added to the proved statements. "and proved statements only" is now true.
  - s06 sign paragraph (`s06_model_2.tex` l.55): the three forms of the dilaton leaves are described (DOSS s.2.3: marked points, a
    vertex is an integral of `psi` classes over the space with these points; DNOPS s.2.3: `psi` classes of legs and edges on that
    space, then the dilaton points are forgotten by push-forward; PPZ Defs. 2.2, 2.5: forget first, `psi` classes on the smaller
    space); by (ii) at every vertex the last two give the same classes, whose integrals are the sums over the dilaton leaves of the
    vertex integrals of the first form, when the translation is the series of the dilaton weights; so DOSS's weights fix the sign.
  - proof of `prop:mod-Rtr` (ii) (`s06_model_3.tex` l.31): a DOSS dilaton leaf is a marked point of its vertex space carrying
    `psi^k`, `k >= 2`; by (ii) at every vertex these give the same vertex integrals as the translation of Def. 2.5 by the series.
* **Input lists**: s01 (Q) list for PPZ15 now `s.0.2, (1); Def. 0.3 (iii); Remark 2.1 and (6); Def. 2.2, with the convention stated
  after it; Prop. 2.4; Def. 2.5, with the convention stated after it; Props. 2.8 and 2.9; Prop. 2.12, with the formula (8) in its
  proof; Def. 2.13; Prop. 2.14`, "(8) only through `lem:set-pushins`"; "Where the inputs enter" (s01 l.17) adds `lem:set-pushins`
  and `prop:mod-Rtr` (through the lemma); s01 l.21 (inputs of Corollary B) as in s10 l.40, naming the lemma; s10 l.33 adds
  `lem:set-pushins` "(in the reading of the topological-recursion input, `prop:mod-Rtr` (ii))"; s02 l.16 (inputs of s02) names
  [PPZ15, (8)]. The two conditionality sentences (s01 l.19, s10 l.35) are unchanged: (8) is a pin inside the quoted work PPZ15 of
  (Q), and s10 refers to "the statement numbers listed in `sec:int-inputs`".

## 2. Pins, read in the full texts

* **PPZ15** (arXiv:1303.1043v3; the PDF in the earlier scratch folders, 468647 bytes, identical in the three copies; text and
  HTML renderings): (8) p. 27, in the proof of Prop. 2.12: `p^*(psi_h^d) = psi_h^d - Delta_(h,n+1) p^*(psi_h^(d-1))`, "the
  well-known formula", `Delta_(h,n+1)` "the divisor in the moduli space at the vertex carrying h corresponding to curves with a
  genus 0 component carrying only the markings h and n+1"; Def. 2.5 p. 23 with the sentence after it (`Omega(..T(psi_i)..) =
  sum_(k>=2) psi_i^k Omega(..T_k..)`); Rem. 2.6 (finiteness; not cited, finiteness proved directly); Prop. 2.12 p. 27; Def. 2.13
  p. 30; the convention sentence p. 20 ("We use the inverse of the R-matrix in all of our formulas in Definition 2.2 ... the same
  convention is used by Givental and Teleman"); s.0.2 (1) p. 3, the gluing map `xi_Gamma` (degree `|Aut Gamma|`, footnote 1).
  Also noted: footnote 2 of PPZ (`kappa_i = pi_*(psi_(n+1)^(i+1))`), the AC convention; and the proof of Prop. 2.9 (p. 25) uses
  the same degree count for genus-0 vertices with `kappa`-legs.
* **AC96** (alg-geom/9406008v1; read in the TeX source of the arXiv e-print, which carries the equation numbers):
  - (1.7): `pi_n*(psi_1^(a_1) .. psi_(n-1)^(a_(n-1)) psi_n^(a_n+1)) = psi_1^(a_1) .. psi_(n-1)^(a_(n-1)) kappa_(a_n)`, the push-forward
    formula (as already used in `lem:set-dilaton`; "kappa_0 = 2g - 2 + n" stands just after the definition of the kappa classes);
  - (1.9): the string equation `pi_n*(psi_1^(a_1) .. psi_(n-1)^(a_(n-1))) = sum_(a_j>0) ... psi_j^(a_j - 1) ...`;
  - (1.10): `kappa_a = pi_n^*(kappa_a) + psi_n^a`, the kappa pull-back formula; (1.12): the iterated push-forward
    `(pi_(k+1)..pi_n)_*(psi_1^(a_1)..psi_k^(a_k) psi_(k+1)^(a_(k+1)+1)..psi_n^(a_n+1)) = psi_1^(a_1)..psi_k^(a_k) R_(a_(k+1)..a_n)`.
  - AC96 does not state the divisor formula `psi_i = pi^* psi_i + D_i` (the section divisors `D_i` are only defined in s.1, before
    the definition of `psi_i`, `K` and `kappa_i`).
  - **Why no AC96 formula enters the lemma.** (1.7), (1.9), (1.10), (1.12) are statements about push-forwards of `psi`/`kappa`
    monomials. With the projection formula they give the comparison only when the CohFT class is pulled back from `Mbar_(h,n)`, i.e.
    for a TFT (degree-0 classes). Corollary B Step 2 needs the comparison for the CohFT `M.Cexp(TFT)`, whose classes are not
    pulled back, and the s06 reading needs it as a class identity (DNOPS's form); by Poincare duality this requires the class
    identity (i), i.e. the divisor formula with the vanishing of `psi` on the genus-0 component. The review's alternative "pin
    [AC96, (1.9)]" rested on a misremembered number: (1.9) is the string equation. The pin used is [PPZ15, (8)], inside a work
    that is already a quoted input; no new bibliography entry.
* **DOSS** (1211.4021v1, TeX source of the e-print, s.2.3 "Expressions in terms of graphs"): dilaton leaves `L^bullet` with labels
  `k > 1`, weight `[z^(k-1)](-(R(-z))^i_1)` (eq. unmarkedLeaf); a vertex with n ordinary leaves and half-edges and m dilaton leaves
  is `hbar^(g-1) int_(Mbar_(g,n+m)) psi_1^(k_1) .. psi_(n+m)^(k_(n+m))` (eq. vertices). So the dilaton leaves are kept as marked
  points (insertion form), as the review said. Sign: `sum_(k>=2) z^k [z^(k-1)](-R^(-1)(z)1) = z1 - zR^(-1)(z)1 = t_R` (DOSS's
  `R(-z)` = DNOPS/PPZ `R^(-1)(z)`): PPZ's sign.
* **DNOPS** (1509.06954v2, source in the scratch folder of an earlier review): s.2.3 decorates each dilaton leaf by
  `-psi(Id - R^(-1)(psi))e_1`, takes the `psi` classes of leaves and edges on `Mbar_(g(v),k(v))` (the valency `k(v)` includes the
  dilaton leaves), then `f_Gamma` first forgets the dilaton points and then sews; `alpha'_(g,k) = sum_Gamma (f_Gamma)_* alpha_Gamma`.
  This is the "second form" of the s06 paragraph; it equals PPZ's form by (ii) (class level).

## 3. Exact hypotheses (as stated)

* (i): only `2h-2+n > 0` (all `Mbar_(h,n+m')` stable), `c_s >= 1`; statement in `H^*(Mbar_(h,n+m); Q)`, extended linearly.
* (ii): `Omega` a CohFT [PPZ15, Def. 0.3] (unit not needed; it only names `d^n Fbar_h`) on a free module of finite rank over a
  commutative Q-algebra (the extension of the PPZ definitions to such rings is the sentence of `def:set-graphsum`); `t` in
  `z^2 V[[z]]` as in [PPZ15, Def. 2.5].
* (iii): `Omega` with unit, `R = 1 + O(z)` symplectic (`def:set-graphsum`), `g >= 2`, leg-free graphs.
* Uses: `thm:set-GT` (`TFT` on `H` over `Q(zeta)[[q]]`, unit `1 = H^0`); `def:mai-Fcf` (colour basis, unit `1_TFT`); `cor:mai-gap`
  Step 1 (`Cexp(TFT)` over `Q(zeta)((s))`), Step 2 (`M.Cexp(TFT)`); `prop:mod-Rtr` (ii) and the s06 paragraph (vertex form (ii),
  `TFT^TR` with its unit).

## 4. Effect on the chain and the conditionality

* Corollary B: the step of I-1 is now proved (lemma), from [PPZ15] definitions already in (Q) plus [PPZ15, (8)].
* Theorem A: the lemma now enters C9 through the proof of `prop:mod-Rtr` (ii) (the reading of the TR input's dilaton sign);
  REVIEW_SYNTHESIS s.0 ("One nuance ...") is resolved by it. The only quoted item added to Theorem A is the pin [PPZ15, (8)],
  listed in (Q); the quoted works and the conditionality sentences are unchanged. AC96 enters only through `lem:set-dilaton`,
  as before.

## 5. Build (scratch folder `allgenus_fix_fix_lemma/`, copy of main.tex, preamble.tex, macros.tex, refs.bib, sections/)

pdflatex, bibtex, pdflatex x3 (a clean build needs the third pass for "Label(s) may have changed", as in the earlier builds):
162 pages (was 159; the lemma takes about 1.5 pp); no error; no undefined or multiply defined reference or citation; no
LaTeX/package warning; the same 5 overfull boxes as before (8.6, 1.5, 1.9, 2.2, 3.9 pt); BibTeX "missing pages in ABDKS25".
All edited files stay below 250 lines (s02_setup_4.tex: 153).

## 6. Left for others (files outside this task's edit set)

* `DRAFT_NOTES.md` s.0 item (def:set-graphsum) and the table row `s02_setup_4.tex:104`, `DRAFT_NOTES_2.md` row: the draftnote is
  gone (16 remain); the statement "not used by the chain" there is superseded (s.4 above).
* `LABELS.md`: add `lem:set-pushins` (s02.3, after `def:set-graphsum`). `BIBPLAN.md`: PPZ15 pins (8) and s.0.2 (1).
* Pointers that rely on the identification through `def:set-graphsum` and could cite the lemma directly: `s07_comparison.tex` l.38
  (`def:flh-outer`, "is the term of Gamma ..."), `s02_setup_8.tex` l.28 (`def:set-outer`, "the outer graph sum is the leg-free part
  of R'.(M.TFT)", for `R'(0)` an isometry this also needs `lem:flh-frame`), `s04_transport_4.tex` l.40. Optional: in `thm:mod-tr`
  Step 1 (`s06_model_3.tex` l.35) a pointer to the paragraph after `prop:mod-Rtr` for DNOPS's placement of the `psi` classes.
* REVIEW_SYNTHESIS_2 s.8 item 10 (human check of the comparison): the proof is now in the note (Lemma 2.22); the human check
  reduces to [PPZ15, (8)] and the reading of Def. 2.2 as a push-forward along (1).
