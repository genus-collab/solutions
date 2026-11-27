# Final review synthesis (key review/critic), part 2 of 2

Continuation of `REVIEW_SYNTHESIS.md` (verdict s.0, method s.1, chain s.2, dependency graph s.3, conditionality s.4, binding
rules s.5). File:line refer to `NOTE/sections/` unless `main.tex` or `macros.tex` is named.

## 6. The four statements across abstract, s01 and body

| statement | abstract (main.tex l.18) | s01 | body | verdict |
|---|---|---|---|---|
| Theorem A | every `g >= 2`, every frame of a two-parameter rational family; edge graphs regular; vertex graph `c_g t_c^(2-2g)` + regular; `c_g = (-1)^(g-1) B_2g/(2g(2g-2) 5^(g-1))` | `thm:int-main` (s01_intro_2 l.9-23), same data and constant, also written `abs(B_2g)/(...)` | `thm:mai-main` (s10_main l.9-21), same, plus "the hatted terms are `G_Gamma = u^(2g-2) T_Gamma`" (`lem:flh-hat`); `def:set-N` is the same statement for one `g`; `eq:set-cg` | consistent |
| Corollary B | "the conifold gap ... in every genus" | `cor:int-gap` (l.27-34), `(alpha, beta)` in `Q^2`, `g >= 2` | `cor:mai-gap` (s10_main_2 l.1-8), "every conifold frame" | consistent; abstract should say `g >= 2` |
| Theorem C | "equality of two normalisations ... canonical at the conifold point, QRR at `q = 0`" | `thm:int-D` (l.44-47) = `thm:gam-D` + `cor:mai-actual` ("Equivalently, `Q = Pcan Q^m`") | `thm:gam-D` (s09_gamma_3 l.83-90), `cor:mai-actual` (s10_main_3 l.49-56) | consistent |
| Theorem Gamma | "identity of formal stationary phase for a Mellin-Barnes integral, giving the value at `q = 0` of its saddle series" | `thm:int-gamma` (l.53-68): `c_n` in `A_L`; (a) `c(z;1) = h(z)`; (b) `c(z;0) = G_inf(z)`; summand in `sigma` | `thm:gam-main` (i), (ii), first assertion of (iv); `thm:gam-ginf`; `def:gam-c` | consistent ("lambda = 1 in the notation there" = "equivariant parameter 1") |

The constant `c_g` appears identically in five places (abstract, `eq:set-cg`, `thm:int-main`, `thm:mai-main`, `thm:mod-gap`);
`c_2..c_5` recomputed. Theorem A (ii) and `def:set-N` write "`T_Gamma` in `c_g t_c^(2-2g) + Q[[u]]`", while `thm:red-main` (b) and
`cor:red-can` still write "`= c_g T^(2-2g) + Q[[u]]`" (s05_reduction l.59, l.94): harmless, but one form should be used.

## 7. Remaining issues

Categories: [M] mathematics or rigour, [K] conditionality, [C] citation, [R] binding rule, [P] presentation.
"Before submission" marks what should be fixed before the note leaves the authors; none of the items is a gap in Theorem A.

* **I-1 [M, before submission] Corollary B uses a step written only inside a `\draftnote`.** `def:set-graphsum`
  (s02_setup_4 l.94-108; the note at l.104) identifies `Cont_Gamma[E_R, q_R; Omega]` with the term of `Gamma` in
  `Int_(Mbar_g) (R.Omega)_(g,0)`; this compares the push-forward form of the translation of [PPZ15, Def. 2.5] with insertions.
  Used by: `thm:set-GT` "In particular" (s02_setup_4 l.111-116); the third property after `def:mai-Fcf` (s10_main l.83,
  `ev(F_g^NN) = F_g`); `cor:mai-gap` Step 1, last sentence, and Step 2, "integrate over `Mbar_g`" (s10_main_2 l.25, l.31). Not used
  by Theorem A as written: `T_Gamma` and `F_Gamma` are contractions, and `thm:mod-tr` quotes DNOPS (2-35), an integral of the class
  `R.omega^TR` formed with the push-forward action of DNOPS s.2.3 (confirmed in the public ar5iv text of 1509.06954 during this
  review). One nuance: the sign of the dilaton leaves used in C9 is taken from DOSS s.2.3 (s06_model_2 l.55, s06_model_3 l.31),
  whose graph weights `[z^(k-1)](-R(-z))` act on leaves kept as marked points, i.e. in the insertion form (ar5iv text of
  1211.4021), while DNOPS s.2.3 prints the opposite sign. The note also justifies the sign by unit preservation (s06_model_2
  l.55, `sec:sgn-R`). Reading DOSS's insertion-form theorem as DNOPS's push-forward class with PPZ's sign is exactly the comparison
  of this item; the lemma recommended below would therefore also make the reading of the TR input in C9 self-contained. The
  argument in the draftnote is correct: `psi_j - p^* psi_j` is a sum of divisors `D_(j,S)` with `j` and a
  set `S` of translation points on a rational component `Mbar_(0,|S|+2)`; the translation `t_R` lies in `z^2 V[[z]]`, so the classes
  `psi` of `S` have degree at least `2|S| > |S| - 1`, the dimension of that component, and every correction vanishes. Until it is a
  lemma, s01_intro_3 l.21 and s10_main_2 l.40 ("and statements proved in the note") overstate. Fix: a two-line lemma in s02.3
  (induction on the number of forgotten points), cited at the three places and in the sign paragraph of s06; or pin
  [AC96, (1.9)] and add it to the inputs (OUTLINE rule 8).
* **I-2 [K] The role of (T) differs between Theorem A and Theorem C.** s01_intro_3 l.13 says that in Theorem A (T) "identifies
  that R-matrix"; `def:set-rows` (s02_setup_4 l.9-14) defines R as the R-matrix of the cohomological field theory. Theorem C is a
  statement about the rows of the same R, but its list of inputs (s01_intro_3 l.21; s09_gamma_3 l.92) does not name (T). Both
  readings are coherent: R is already determined by (b) of `def:set-rows`, `R(0) = 1` and (QRR) (the QDE fixes R up to a constant
  diagonal series in `z`, which the full value at `q = 0` fixes), and the proofs of C6, C11, C12 use only these. Choose one: add
  "(T), which identifies R" to Theorem C, or say once in `def:set-rows` that R is characterised by (b) and (QRR) and that (T) is
  used only to identify it with the R-matrix of the CohFT (then for Theorem A it gives the meaning of `T_Gamma`, and for
  Corollary B the total).
* **I-3 [K] The abstract does not state the same inputs as s01 and s10.** main.tex l.18: no [AC96, (1.7)]; "Teleman's
  classification" (s01, s10: "Teleman's reconstruction"); (QRR) singled out as the condition and the other inputs called "cited
  theorems", while s01/s10 put all on one footing ("holds modulo"); "in every genus" without `g >= 2`. The abstract is a draft
  (`\draftnote`); align it when it is settled.
* **I-4 [C] Pinned items used on the chain but not listed in s01** (none is a further theorem; add them to (Q) or to the sentence
  "The definitions use ..." of s01_intro_3 l.15, for exactness):
  - [DuninBarkowski:2012bw, s.3.1, s.3.5], the formal-Laplace prescription for the integrals of DNOPS (2-31), (2-32), which
    `prop:mod-Rtr` uses (s06_model_2 l.23; rf5 s.3 item 2 asked for it in BIBPLAN);
  - [MR2276766, Section 6, (10)], the normalisation of the J-function in `thm:set-mirror` (s02_setup l.97);
  - [PPZ15, (5)], the symplectic condition in `def:set-graphsum` (s02_setup_4 l.95);
  - [Eynard:2007kz, Def. 4.3 (4-10)] (the free energies, `conv:mod-kernel`, `cor:mod-trF`) and (5-1), (5-4), (5-13) (the set-up of
    Thm 5.1 read in `thm:mod-var`).
  Also: "Where the inputs enter" (s01_intro_3 l.17) omits `lem:set-modhom` (PPZ Defs. 2.2, 2.13, s02_setup_7 l.37) and
  `thm:set-G0` (PPZ Def. 0.3, as part of (G0)) from the places where PPZ15 enters.
* **I-5 [C] Two pin deviations from GAPS_SYNTHESIS_2 s.3.1 await ratification** (EO (4-30), (5-15) out of (Q); PPZ Prop. 2.8 into
  (Q)). Both are correct (part 1, s.4.3). The conditionality at the level of the quoted works is unchanged.
* **I-6 [P] s01_intro l.16**: "the coefficient `A_n = [z^n](C_0 Q^(0)_0)` is a polynomial in `L` of degree at most `4n+1`
  (`def:set-rows`, `prop:con-deg`)". `prop:con-deg` concerns solutions of HPF with coefficients in `Q[[u]]`; for the actual
  amplitude the statement follows from `thm:gau-R0` (ii)-(iii) (`A^act_n` is a combination of the `A^can_m`, `m <= n`) with
  `prop:con-deg` applied to `Phi^can`. Cite `thm:gau-R0` instead of (or with) `def:set-rows`.
* **I-7 [R] macros.tex l.1**: the comment names the design files "NOTATION.md, NOTATION_2.md" (rule 5: no file names in a `.tex`
  file). Delete the parenthesis. The structural comments at the top of the section parts name only the note's own files; they may
  be stripped before submission.
* **I-8 [R] s06_model_7 l.99**: the draftnote "The closed forms printed in [Hock25, s.3.2] were not re-checked" describes a
  process; say "are not used; only Cor. 1.2 is cited".
* **I-9 [R, authors] s01_intro_3 l.34**: "For `thm:int-gamma` we are not aware of a proof in the literature" is close to a priority
  statement (rule 11). Keep only if the authors want it.
* **I-10 [P] DEC D8 is not met**: 38 remarks, about 10.3 pages (target: about 2). Largest: `rmk:mai-constant` 1.2 pp,
  `rmk:set-qrrderiv` 1.2 pp (with a draftnote, used by no proof), `rmk:gam-necessity`, `rmk:flh-fail`, `rmk:gam-sevenvar` (about
  0.5 pp each). The OUTLINE_9 drop list (`sec:con-finite`, `rmk:gam-sevenvar`, `rmk:gam-general`, `rmk:mai-analogy`,
  `rmk:gau-Uremarks`) is untouched.
* **I-11 [P] Length**: 159 pages against the design's 84 (OUTLINE_9); s02 21 pp (budget 10), s03 24 (10), s04 13 (6), s06 20 (12),
  s08 15 (10), s09 11 (5), s10 10 (4), appendix A 11 (7), C 10 (5). See s.9.
* **I-12 [P] Residual letters with several meanings**, each declared where used: `q` (parameter; loop-space point; point of the
  curve; `q_phi`), `J`, `x`, `lambda`, `mu`, `xi` (DRAFT_NOTES_2 s.5). Acceptable for a draft; a global renaming of the loop-space
  point (e.g. `\mathbf q`) would remove the most frequent clash.
* **I-13 [P] Draftnotes**: 17, all off the chain (part 1, s.3). Four must be resolved before submission whatever the authors
  decide on the remarks: the title and the abstract (main.tex l.12, l.19), D11 (s01_intro l.38), and `def:set-graphsum` (I-1);
  and the version re-checks of CGYZ26 and Br25 (s01_intro l.40, s01_intro_2 l.74, s10_main_2 l.55).
* **I-14 [M, for the reader] What Corollary B is a statement about.** The conifold expansion of `F_g` is defined as that of its
  graph-sum representative `F_g^NN` in the presented ring (`def:mai-Fcf`); the injectivity of `ev` is not claimed, and the
  agreement with the expansion of [CGYZ26] holds on the subring generated by `C_0^(+-1)`, `C_1^(+-1)` over `Q[L^5, L^-5]`, while
  membership of `F_g^NN` there is "neither needed nor proved" (`rmk:mai-constant` (2)-(3)). This is stated in the note; a reader
  comparing with [CGYZ26, Conj. 1.5] should be aware that Corollary B is a statement about this representative.
* **I-15 [P] Typography**: five overfull boxes of at most 8.6 pt (s04_transport_3 l.59-60, s06_model_8 l.20-21, s08_connection_4
  l.16-17, s08_connection_6 l.27-28, one alignment); BibTeX "missing pages in ABDKS25" (an arXiv entry).
* **I-16 [P, outside the note]** The design files are behind the note (DRAFT_NOTES_2 s.6: NOTATION tables 3-4, LABELS_2
  `eq:tra-leaf`, BIBPLAN pins of DOSS, LP20, Hock25, EO; OUTLINE s.0.2 row C10 and list (Q); the derived `u^2/30`;
  SOURCE_ERRATA). Not a defect of the note.

## 8. What a human reader must check first

All statement numbers of the quoted inputs were read, in this project, by AI agents only. In order of risk:

1. **The quoted TR inputs of s06 in the original papers**, with the conventions of `conv:mod-kernel`: DOSS Thm 4.1 as recalled in
   DNOPS s.2.6 (2-30)-(2-35), the action of DNOPS s.2.3 (dilaton leaves forgotten by push-forward), and above all the **sign of the
   dilaton leaves**: DNOPS s.2.3 prints `-psi(I - R^-1(psi))e_1`, DOSS s.2.3 weights `[z^(k-1)](-R(-z))` on kept leaves, PPZ Def. 2.13
   translates by `+psi(1 - R^-1(psi)1)`; the note uses PPZ's sign (s06_model_2 l.55) and the pin of (Q) says so. DNOPS
   Prop. 4.3 (1) and the sentence before its Thm 8.1, at `f = 4W_v` = DNOPS (8-1) with `n = 4`, `a = (0,0,0,0,4v)`; DOSS s.3.1,
   s.3.5; ABDKS Thm 7.1 (7.6), (7.7) with `1/hbar^2`, (7.9), (7.10) (which drives the Gap Theorem through `eq:mod-710`); EO Def. 4.2
   (4-9), kernel (4-3), Thms 4.2, 4.6, 5.1 (5-14) with s.5.3. Then the signs `omega = (-1)^n W^EO`, `F^TR = -F^(h)` and the factor
   `(-5)^(1-h)` of `cor:mod-trF`, which fix the sign and the power of 5 in `c_g`.
2. **(QRR) and its sign**: `R(q = 0) = diag(h_(w_i))` with the weight `-5 zeta^i` (`rmk:set-qrrsign`: `h_1 = -3/20`, not `-11/60`), in
   Givental's convention (columns `R e~_i e^(u_i/z)` solve the QDE; legs `R^-1(psi)`), against [GJR18, (17)], [GJR17, (19)-(20)],
   [Giv01, s.2.2], [MR2276766, Thm 1, Cor. 1]. The identity `D = 1` would fail for the other sign.
3. **(M) and (T) as pinned**: [MR1408320, Thm 11.6] in the case `n = 4`, `r = 1`, `l_1 = 5`, `lambda' = 0` with s.9 and Prop. 11.4;
   [MR2917177, Theorem (2) in (1.15)] with [GJR18, Thm 4.1]; the journal numbering (two draftnotes).
4. **The transfer between the two points** (`thm:gau-R0` (iii)-(iv), `thm:con-glue`): `A^can = A^act D_0` proved in
   `Q(zeta)[[q]]` and read in `Q(zeta)((s))` through `Q(zeta)[L, 1/L]`, then evaluated at `L = 1`. This is where the normalisation
   at `q = 0` reaches `u = 0`.
5. **Theorem P** (`thm:con-Pa`, `thm:con-Pb`, s08.2-8.3) and `prop:gam-hatflat`: the canonical colour-0 function equals the
   full-Stirling saddle series, and `c_n(L) = L^(4n) ph^SP_n(L^-5)`.
6. **Theorem LMI and Theorem Gamma (ii)** (s09.2-9.3): the continuation of the saddle series to the degenerate point along the
   explicit family `p(L)`; only elementary holomorphy and appendix A.
7. **The comparison theorem** (s07): Step 1's degree count and the window `z^(3g-2)`; `lem:flh-comm`, `lem:flh-dil`,
   `lem:flh-frame`, `lem:flh-inputs` (iii)-(iv).
8. **Theorem R and the transport lemmas** (s05, s04): the u-adic summability of `lem:tra-conv` and the leaf point `lem:tra-T`.
9. **The Model Lemma's assembly** (s06.3-6.6: lattice theorem, L-regularity, Theorem V, `thm:mod-MLK`).
10. **The comparison of the two forms of the translation** (I-1), for Corollary B.
11. **The literature paragraphs** (`sec:int-gap`, `rmk:int-constant`, `rmk:mai-constant`): CGYZ26 numbers (read in the v2 HTML),
    CGL26, Br25, FCCZ26, GJR18 Rem. 1.6, 1.9 and 7.2.

## 9. Length: split, and what can move to an appendix

The note is a monograph-length text with three logically separable parts:
* (I) **the universal flat theory and the model**: Theorem U, the gauge identity and symplecticity for one-pole cyclic systems
  (s03.2-3.5), the tower, slice and transport (s04), Theorem R (s05), the model and its CohFT (s02.5, `lem:set-dilaton`), the Model
  Lemma by TR on the model curve (s06), universality in flat form (`cor:mai-univ`, s10.3). Main theorem: `cor:mai-univ`. Inputs:
  (Q) and [AC96, (1.7)] only; nothing about the quintic. About 60 pages.
* (II) **Theorem Gamma**: `def:gam-c`, `prop:gam-rational` (s08.2), Theorem LMI, `thm:gam-main`, `thm:gam-ginf`, the first
  assertion of `cor:gam-poly` (s09) and appendix A; `prop:gam-hatflat`, the rest of `cor:gam-poly` (which uses `thm:con-Pa` (iv)),
  `thm:gam-D` and `cor:gam-free` stay with (III). Self-contained, no quoted input. About 25 pages.
* (III) **the formal quintic**: s02.1-2.4, 2.6-2.7, s03.1 (CYC), s03.6, s07, s08 (Theorem P, glue), s10.1-10.2, appendix C, citing
  (I) and (II). Inputs (M), (G0), (T), (QRR), [PPZ15]. About 70 pages.

**Recommendation: split into two papers**, (I) and (II)+(III), or three if a short paper on local Mellin inversion is wanted. The cut
is clean: s06 uses from s02-s04 only the model and its CohFT, the vanishing solutions, the tower and the slice (OUTLINE
s.0.1); s05 uses only s04 and the dilaton lemma; the s03 results that `cor:mai-univ` needs (Theorem U, `prop:gau-graded`,
`lem:gau-R1`, the proof of `prop:gau-symp`, `thm:gau-G`, `lem:gau-vanq`) are stated for, or (as its proof shows) apply to,
general one-pole systems (`def:gau-can` would have to be restated for a general pair, as `cor:mai-univ` already does); and
`prop:tra-can` becomes the case `A` = quintic system of `cor:mai-univ`. If the note stays one text, move:
* to appendices: the proof of the large-radius cyclic structure (s02.1, `lem:set-grading` .. `lem:set-order3`, about 6 pp); the
  cyclic form of `PF4` with the concomitant (s03.1, about 6 pp); the u-adic convergence machinery (`lem:tra-conv`, about 4 pp);
  the CGYZ26 dictionary of `rmk:mai-constant` into `sec:sgn-lit`; s10.2-10.3 (not used by Theorem A) into a final "Complements"
  section;
* out of the note: `sec:con-finite` (unused), the remarks beyond D8 (in particular `rmk:set-qrrderiv`, `rmk:gam-sevenvar`,
  `rmk:gam-Mcomplex`, `rmk:gam-general`, `rmk:mai-tau`, `rmk:mai-analogy`, `rmk:tra-consistency`), and most of appendix C (keep
  the tables; the rest can be an ancillary file).
This would bring a single text to about 110-120 pages.

## 10. Recommendations (in order)

1. Write the comparison of I-1 as a lemma in s02.3 and cite it in `thm:set-GT`, after `def:mai-Fcf` and in `cor:mai-gap`, and in
   the sign paragraph of s06 (s06_model_2 l.55, s06_model_3 l.31), where it turns DOSS's insertion-form weights into DNOPS's
   push-forward class with PPZ's sign; delete the draftnote. (Or pin [AC96, (1.9)] and list it as an input.)
2. Make the role of (T) the same for Theorems A and C (I-2).
3. Settle title and abstract; give the abstract the inputs of s01 (I-3).
4. Ratify I-5; complete the lists of s01 with the items of I-4; update BIBPLAN and the other design files (I-16).
5. Before further editing, have a human check items 1-3 of s.8 against the original papers.
6. Decide the split (s.9); in any case cut the remarks to about two pages and drop `sec:con-finite`.
7. Small fixes: I-6, I-7, I-8, I-9; one form ("in") for the vertex law in `thm:red-main` (b) and `cor:red-can`.
8. Resolve the author-decision draftnotes (DRAFT_NOTES.md s.2) and the version re-checks; switch `\showdraftnotesfalse` only when
   all 17 are gone.

Files: this file, `REVIEW_SYNTHESIS.md`. Scratch (not part of the project; may be cleared):
`<user home>/AppData/Local/Temp/claude/d7c1b608d03d194a/37f24460-3601-4c85-b123-5156e19eb475/scratchpad/allgenus_review_review_critic/`
(build and `tools/` scripts of part 1, s.1).
