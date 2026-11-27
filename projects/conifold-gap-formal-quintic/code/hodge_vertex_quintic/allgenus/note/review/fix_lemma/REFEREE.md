# Referee report: the lemma `lem:set-pushins` and its uses (key fix/lemma_referee)

Author-facing record; not part of the note. Object: the lemma written for review item I-1 (record `review/fix_lemma/LEMMA.md`),
i.e. `lem:set-pushins` (Lemma 2.22 of the build, `sections/s02_setup_4.tex`) and its uses: `def:set-graphsum`, `thm:set-GT` and
the paragraph after it, the third property after `def:mai-Fcf` (`s10_main.tex` l.83), `cor:mai-gap` Steps 1-2 and the closing
input sentence (`s10_main_2.tex`), the s06 sign paragraph (`s06_model_2.tex` l.55), the proof of `prop:mod-Rtr` (ii)
(`s06_model_3.tex` l.31), the input lists (`s01_intro_3.tex` l.9, 17, 21; `s02_setup.tex` l.16; `s10_main.tex` l.33).

## Verdict: SOUND-AFTER-FIXES

No mathematical error. The statement is correct, every step of the proof is correct and complete, every pin was checked in
the full text, the hypotheses hold at every use, and there is no circularity. I made three small corrections (s.7). One was
in the proof: (*) was proved for `N := n` but used with `N = n+m-1`. Two were about exactness: the input lists of Corollary B,
and a missing pointer in `thm:mod-tr` Step 1. The previous agent's account in LEMMA.md is accurate in every point I checked.

## 1. Sources (downloaded afresh into my scratch folder; arXiv e-print sources and PDFs)

* **PPZ15** = arXiv:1303.1043v3. I read the TeX source (`3spin2015.tex`; `\setcounter{section}{-1}`, so the introduction is
  s.0) and the PDF pages 3, 20, 27 and 30. I recounted the numbering from the source. The `proposition` counter is shared by
  remark, definition and lemma, numbered within the section. The `eqnarray` (3) has one numbered line.
  - s.0.2 "Stable graphs", **(1)** = `xi_Gamma: Mbar_Gamma -> Mbar_(g,n)`, p.3; footnote 1: "The degree of xi_Gamma is |Aut(Gamma)|".
  - **Def. 0.3** = CohFT (axioms (i), (ii)); with (iii), a CohFT with unit.
  - **(5)** = `R(z) = 1 + R_1 z + ...`, with the symplectic condition `R(z)R^*(-z) = 1` in the display right after it. PPZ
    themselves say "R an element of the group (5)" on p.20, so the note's pin `[PPZ15, (5)]` for the symplectic group is
    PPZ's own usage.
  - **Rem. 2.1** contains **(6)**, the power-series property for `R^{-1}`.
  - **Def. 2.2**, p.20: (i) Omega at the vertices, (ii) `R^{-1}(psi_l)` on the legs, (iii) the edge
    `(eta^{-1} - R^{-1}(psi')eta^{-1}R^{-1}(psi'')^t)/(psi'+psi'')`, and `(R Omega)_(g,n) = sum_Gamma |Aut Gamma|^{-1} Cont_Gamma`.
    The convention sentence on p.20 reads "We use the inverse of the R-matrix in all of our formulas in Definition 2.2 ... the
    same convention is used by Givental and Teleman".
  - Prop. 2.3 is the CohFT property. **Prop. 2.4** is the left group action.
  - **Def. 2.5**: `(T Omega)_(g,n)(v) = sum_m (1/m!) (p_m)_* Omega_(g,n+m)(v, T(psi_(n+1)), ..., T(psi_(n+m)))`, with
    `T in z^2 V[[z]]`. The sentence after it reads `Omega(..T(psi_i)..) = sum_(k>=2) psi_i^k Omega(..T_k..)`.
  - Rem. 2.6 is finiteness. Prop. 2.7: `T Omega` is a CohFT. **Prop. 2.8**: translations form an abelian action.
    **Prop. 2.9**: `T_a R Omega = R T_b Omega`. Rem. 2.10. Lemma 2.11 is the fibre-square lemma.
  - **Prop. 2.12**, p.27: `T_a R Omega = R T_b Omega` with `T_b(z) = z[1 - R^{-1}(1)](z)` is a CohFT with unit 1. Its proof
    contains **(8)**, p.27: `p^*(psi_h^d) = psi_h^d - Delta_(h,n+1) p^*(psi_h^(d-1))`, "the well-known formula", where
    `Delta_(h,n+1)` is the divisor of curves with a genus-0 component carrying only the markings h and n+1.
  - **Def. 2.13**, p.30: `T(z) = z.1 - zR^{-1}(z)(1) in z^2V[[z]]`, `R.Omega = R T Omega`. **Prop. 2.14**, p.30: the
    unit-preserving action is a left group action.
  - All statement numbers, pages and conventions cited by the lemma and its uses are correct.
* **AC96** = alg-geom/9406008v1. I read the TeX source (equation tags with `\leqno`) and PDF p.5.
  - s.1 defines the section divisors `D_i` and `psi_i`, `K`, `kappa_i = pi_*(K^(i+1))`. The text "Notice that
    `kappa_0 = 2g-2+n`" is on p.5.
  - (1.5) gives `kappa` against `kappa~`.
  - **(1.7)** is the push-forward formula
    `pi_n*(psi_1^a1..psi_(n-1)^a(n-1) psi_n^(an+1)) = psi_1^a1..psi_(n-1)^a(n-1) kappa_an`, the classes `psi_j` on the left
    being those of the larger space.
  - (1.8) is the restriction of kappa to strata. (1.9) is the string equation. (1.10) is `kappa_a = pi_n^*(kappa_a) + psi_n^a`.
    (1.11) is the dualizing-sheaf comparison. (1.12) is the iterated push-forward.
  - The formula `psi_i = pi^* psi_i + D_i` is **not** stated.
  - The use of (1.7) in `lem:set-dilaton` (a = 0, kappa_0 on p.5) is literal and correct. AC96 is not used by the new lemma.
* **DOSS** = arXiv:1211.4021v1, TeX source, s.2.3 "Expressions in terms of graphs".
  - Dilaton leaves (labels `k > 1`) carry `[z^(k-1)](-(R(-z))^i_1)` (eq. `unmarkedLeaf`).
  - Ordinary leaves carry `R(-z)` applied to the input (eq. `markedLeaf`).
  - The edges are `(delta - R(-z)R(-w)^T)/(z+w)`.
  - A vertex with m dilaton leaves is `hbar^(g-1) int_(Mbar_(g,n+m)) psi_1^k1 ... psi_(n+m)^k(n+m)` (eq. `vertices`).
  - The sum carries `1/|Aut(gamma)|`, with the leaves inside gamma (eq. `graphs-Givental`).
* **DNOPS** = arXiv:1509.06954v2, TeX source, s.2.3.
  - Each dilaton leaf is decorated by `-psi(Id - R^{-1}(psi))e_1`.
  - The classes psi of legs and edges live on `Mbar_(g(v),k(v))`, where `k(v)` includes the dilaton leaves. `f_Gamma` first
    forgets the dilaton points and then sews: `alpha'_(g,k) = sum_Gamma (f_Gamma)_* alpha_Gamma`.
  - s.2.4-2.6 are the superpotential, TR, and TR via CohFTs.

## 2. The statement

* (i) A class identity on `Mbar_(h,n+m)`, `2h-2+n>0`. The hypothesis needed is `c_s >= 1`; it is correct as stated. It is the
  standard fact `psi_j = pi^*psi_j + D_(j,N+1)` together with `psi_(N+1) D_(j,N+1) = 0`, iterated.
* (ii) A class-level and an integrated identity for a CohFT without a unit hypothesis and for `t in z^2V[[z]]`. The two sides
  match PPZ Def. 2.5 and the convention after it literally. The last sentence is the definition of `d^n Fbar_h` in
  `def:set-graphsum` (`q + z1 = t`). Correct.
* (iii) For `Omega` with unit and `R = 1 + O(z)` symplectic, `g >= 2`, leg-free graphs:
  `Cont_Gamma[E_R, q_R; Omega] = |Aut Gamma|^{-1} int` of the Def. 2.2 contribution of `Gamma` to `R(t_R Omega) = R.Omega`.
  The closing sentence ("again a CohFT with unit, to which (iii) applies") is exactly Prop. 2.12 with `T_b = t_R`. Correct.

## 3. The proof, step by step

* (i), m = 1.
  - (8) with d = 1 gives `pi^*psi_j = psi_j - delta_j`.
  - (8) with d = 2 gives `pi^*(psi_j^2) = psi_j^2 - delta_j psi_j + delta_j^2`. Comparing this with `(psi_j - delta_j)^2`
    gives `psi_j delta_j = 0`. I checked this algebra: `-delta psi = -2 psi delta`.
  - The transposition of j and N+1 fixes `delta_j` and exchanges `psi_j` and `psi_(N+1)`, so `psi_(N+1) delta_j = 0`, hence
    `psi_j psi_(N+1) = pi^*(psi_j) psi_(N+1)`.
  - Replacing one factor at a time needs a single `psi_(N+1)` (`c >= 1`). This gives (*). Correct.
  - Using (8) for d = 1, 2 is within its statement: PPZ state it for the pull-back of every power `psi_h^d`, d >= 1. It also
    applies on a single vertex space `Mbar_(h,N+1) -> Mbar_(h,N)`.
* (i), induction step.
  - `p_m = p_(m-1) o pi`. Apply (*) with `N = n+m-1`, `c = c_m`, then (i) for m-1, then (*) again with `k_j = 0`.
  - I re-derived every equality. Correct.
  - The text set `N := n` and later used (*) with `N = n+m-1`. (*) holds for every N with `2h-2+N > 0`, because its proof
    never used `N = n`, but as written it had been proved only for `N = n`. Fixed (s.7, item 1).
* (ii).
  - Def. 2.5 with its convention, then the projection formula `p_*(Z) b = p_*(Z p^*b)`.
  - (i) applies to each term of `Z_m`, since every `c_s >= 2 >= 1`. Then `int p_* W = int W`.
  - Degree bound: `2m <= 3h-3+n+m` iff `m <= 3h-3+n`, and `c_s <= dim`. Correct.
* (iii).
  - Def. 2.13 gives `R.Omega = R(t_R Omega)`.
  - The Def. 2.2 contribution is the push-forward along `xi_Gamma` of (1). This is a reading: PPZ place the classes on the
    vertex spaces and call the result a class on `Mbar_(g,n)`. PPZ's proof of Prop. 2.3 confirms it: it pulls graph
    contributions back along `xi_Phi` "according to the rule given in Section 0.3".
  - The edge of Def. 2.2 (iii) equals `E_R` of `def:set-graphsum`. In components,
    `(R^{-1}(psi) (x) R^{-1}(psi'))eta^{-1} = R^{-1}(psi) eta^{-1} R^{-1}(psi')^t`, and `E_R` is symmetric under the swap,
    so orientation is irrelevant.
  - The push-forward preserves integrals. Expand into monomials (finitely many survive), apply Fubini on `prod_v Mbar_(gv,nv)`
    (self-loops included), then (ii) at every vertex with `q_R + z1 = t_R`. This gives `Cont_Gamma`, factor
    `|Aut Gamma|^{-1}` included, as in Def. 2.2. Correct.
* PPZ's own proof of Prop. 2.9 ("we add to every vertex ... kappa-legs") uses the same comparison tacitly. The lemma does not
  depend on that proof.

## 4. The sign of the dilaton translation

* PPZ Def. 2.13: `T(z) = z1 - zR^{-1}(z)1`; Prop. 2.12: `T_b = z[1 - R^{-1}(1)]`. The note's `t_R = z1 - zR^{-1}(z)1` and
  `q_R = R^{-1}(z)(-z1) = -z1 + t_R` agree.
* DOSS s.2.3: ordinary leaves carry `R(-z)`, the operator that PPZ and DNOPS write `R^{-1}(z)` on legs. The edges agree with
  this reading.
  - So `sum_(k>=2) z^k [z^(k-1)](-R(-z)1) = z(1 - R(-z)1) = z1 - zR^{-1}(z)1`, which is PPZ's sign.
  - DNOPS prints `-psi(Id - R^{-1}(psi))e_1`, the opposite sign.
* The note states exactly this:
  - the s06 paragraph;
  - `prop:mod-Rtr` (ii): `psi 1hat - psi R_TR^{-1}(psi) 1hat`;
  - the s01 (Q) list: "its dilaton leaves taken with the sign of [PPZ15, Def. 2.13]".
* The lemma gives the two comparisons the reading needs:
  - (ii), class level: DNOPS's form, with the classes psi on the spaces with the dilaton points followed by push-forward,
    equals PPZ's push-forward form.
  - (ii), integrated: DOSS's insertion form gives the same integrals.
  - The `1/m!` of Def. 2.5 is the orbit count of DOSS's `1/|Aut gamma|`, with the leaves included. "The sums over the
    dilaton leaves" in the s06 paragraph covers this.

## 5. Hypotheses at each use

* `thm:set-GT`: `TFT` on `H` over `Q(zeta)[[q]]` with unit `1`; `R` symplectic (`def:set-rows` (a)), `R(0) = 1`;
  `F_g = int Omega_(g,0)` (`def:set-fq`, checked); g >= 2. OK.
* Third property after `def:mai-Fcf`:
  - `ev(R~)` is the matrix of `R` in the normalised canonical basis. It is symplectic for the identity pairing and
    `ev(R~)(0) = 1`.
  - `TFT` in the colour basis has unit `1_TFT`. So `ev(R~).TFT` is `R.TFT` written in that basis, and (iii) and then (T)
    apply.
  - `ev` commutes with `E_X`, `q_X` and the vertex tensors. OK.
* `cor:mai-gap` Step 1:
  - `Cexp(R~)` is symplectic (shown) with `Cexp(R~)(0) = 1`.
  - `Cexp(TFT)` is a TFT with unit. I checked that the axioms are identities of Laurent monomials in `Delta_i^(1/2)`:
    gluing, `Delta_i^(-1/2) Delta_i^((2g-1+n)/2) = Delta_i^((2g-2+n)/2)` and `Omega_(0,3)(e~_i, e~_j, 1) = delta_ij`.
  - For symplectic X, the edge of `def:mai-Fcf` equals `E_X` of `def:set-graphsum`. OK.
* `cor:mai-gap` Step 2: `M` is symplectic with `M(0) = 1`, so `M.Cexp(TFT)` is a CohFT with unit `Cexp(1_TFT)` (Prop. 2.12).
  `X` is symplectic with `X(0) = 1`. OK.
* `prop:mod-Rtr` (ii) and the s06 paragraph:
  - `Omega = TFT^TR` (a TFT; no unit is needed for (ii)), over a Q-algebra (C).
  - `t = psi(1 - R_TR^{-1}(psi))1hat` lies in `psi^2 V[[psi]]`, since `R_TR(0) = 1`.
  - Legs and half-edges are main points with arbitrary monomials, as (ii) allows. OK.
* Extension of scalars: (i) holds in `H^*(Mbar; Q)`, and (ii) and (iii) are linear. The lemma is proved over every commutative
  Q-algebra, as its statement says.

## 6. Circularity

* The lemma uses only:
  - PPZ: Def. 0.3, (1), Defs. 2.2, 2.5, 2.13 with their conventions, (8), Prop. 2.12;
  - the definitions in `def:set-graphsum`.
* `def:set-graphsum` refers forward to the lemma only to state the agreement; its definitions do not depend on it.
* The users of the lemma do not occur in its proof: `thm:set-GT`, `def:mai-Fcf`, `cor:mai-gap`, `prop:mod-Rtr` (ii) and the s06
  paragraph. No cycle.

## 7. Changes made (only in files of the previous agent's edit set)

1. `sections/s02_setup_4.tex`, proof of (i), first sentence. It now reads "Let first m=1. For N>=0 with 2h-2+N>0 let
   pi: Mbar_(h,N+1) -> Mbar_(h,N) forget the last point, ... (N, pi and delta_j are local letters)". At the end of the m=1 step
   it reads "For N=n, pi=p_1 and (*) is (i) for m=1". So (*) is proved for every N, as the induction step uses it with
   `N = n+m-1`. Before, the text read "write pi := p_1 and N := n". The file still has 153 lines.
2. Corollary B input lists, `sections/s01_intro_3.tex` l.21 and `sections/s10_main_2.tex` l.40 (kept identical):
   - added `Def. 0.3`: Step 1 now invokes "a TFT with unit", Prop. 2.12 and the lemma need it, and Theorem A's (Q) list has it;
   - added "with the conventions stated after them" to "Defs. 2.2 and 2.5". The lemma uses the reading of `T(psi)` stated
     after Def. 2.5 and the `R^{-1}` convention of p.20; (T) is matched to PPZ's action through the p.20 sentence.
   - Theorem A's (Q) list already had both items. The conditionality, at the level of quoted works, is unchanged.
3. `sections/s06_model_3.tex`, `thm:mod-tr` Step 1. The parenthetical on the action of DNOPS s.2.3 now ends "...; by
   `lem:set-pushins` (ii) its classes are those of the push-forward form of [PPZ15, Def. 2.5], see the paragraph after
   `prop:mod-Rtr`".
   - This is the place where DNOPS's action is equated with `R.omega^TR`; before, the link went only through "the hypotheses on R
     are `prop:mod-Rtr`".
   - This makes the reading of the TR input in C9 self-contained, as REVIEW_SYNTHESIS_2 I-1 asked.

No other file was touched; `refs.bib` is unchanged. The lemma statement, the proofs of (ii) and (iii) and all the other uses are
unchanged: they are correct as written.

## 8. Build (scratch folder `allgenus_fix_fix_lemma_referee/build/`)

* Copied `main.tex`, `preamble.tex`, `macros.tex`, `refs.bib` and `sections/` after the edits, and removed the old auxiliary
  files.
* Ran pdflatex, bibtex and pdflatex x3. The third pass was needed for "Label(s) may have changed", as in earlier builds.
* Result: 162 pages; no error; no undefined or multiply defined reference or citation; no LaTeX or package warning.
  - The same 5 overfull boxes as before.
  - BibTeX "missing pages in ABDKS25" (pre-existing).
* Rendered pages 21-22 (the lemma) and 77 (`thm:mod-tr` Step 1) were inspected. `NOTE/main.pdf` was not rebuilt.

## 9. Open points for the authors and for others (no change made)

* **The pin [PPZ15, (8)]** is a "well-known formula" recalled inside the proof of Prop. 2.12, now in (Q).
  - It is correct and is used only for d = 1, 2.
  - The authors should ratify the pin. Alternatively they may replace it by a textbook reference for
    `psi_i = pi^*psi_i + D_i` and `psi_i D_i = 0`, which would need a new bibliography entry.
  - AC96 cannot serve: (1.7) and (1.9) are push-forward statements. With the projection formula they give the comparison only
    for classes pulled back from the smaller space (TFTs). Corollary B Step 2 (`M.Cexp(TFT)`) and the DNOPS reading need the
    class identity.
* **Human check (REVIEW_SYNTHESIS_2 s.8 item 10)** now reduces to two points:
  - (8);
  - the reading of PPZ Def. 2.2's `Cont_Gamma` as `xi_Gamma*` of the placed classes. This is implicit in PPZ, and supported by
    footnote 1 and by the proof of Prop. 2.3.
* **DNOPS s.2.3** prints `alpha'_(g,k) = sum_Gamma (f_Gamma)_* alpha_Gamma` without the factor `1/|Aut Gamma|`.
  - The note reads that action as PPZ's ("edges as in [PPZ15, Def. 2.2]"), which supplies the factor; DOSS has it.
  - Optional: the authors may say so in the (Q) item, next to "its dilaton leaves taken with the sign of [PPZ15, Def. 2.13]".
* **Files outside this edit set** (as in LEMMA.md s.6):
  - DRAFT_NOTES(_2): the `def:set-graphsum` draftnote is gone;
  - LABELS: add `lem:set-pushins`;
  - BIBPLAN: PPZ15 pins (8) and s.0.2 (1);
  - optional pointers at `s07_comparison.tex` l.38, `s02_setup_8.tex` l.28, `s04_transport_4.tex` l.40. These rely on
    `def:set-graphsum`, which now points to the lemma; they are not needed for Theorem A, since `T_Gamma` and `F_Gamma` are
    contractions (remark after `thm:flh-main`, item 3).
