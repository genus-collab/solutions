# Review of the note design (key design_critic), part 2: recommended changes, spot checks, the eight checks

Continuation of `OUTLINE_REVIEW.md` (verdict s.0, required changes R1-R8 in s.1). Same abbreviations.

## 2. Recommended changes (the writers or the assembler may apply them; none blocks the writing)

* **M1. (S3), the implicit function theorem.** BIBPLAN has no entry for it, so it is the only cited input without a pinned statement
  (rule 4). LA s.2.1 says only that "a formal IFT over C[[Q]] would do"; no source writes that route. Keep (S3) (renamed, R2) as the
  input of the proof of record of `thm:gam-lmi`, pin a textbook statement (authors' choice, `\draftnote`), and replace "the formal IFT
  `lem:fsp-ift` also suffices, LA s.2.1" in OUTLINE s.0.2 by "a formal route is suggested in LA s.2.1, not written".
* **M2. Composition of R-actions in `cor:mai-gap`.** The design cites [PPZ15] Prop. 2.4 alone. FR s.0 uses Prop. 2.4 with Prop. 2.9
  (translations), and the unit-preserving action is a group action by Prop. 2.14 (in the v3 text; BIBPLAN key table lists it). Cite
  Prop. 2.14 (or 2.4, 2.8, 2.9) and list it, with (S2) and `lem:gau-confexp`, as the inputs of Corollary B in `sec:int-inputs`.
* **M3. `thm:gam-D`, input sentence.** "Inputs: (QRR) (in `thm:con-glue`) and the appendix-A lemmas; nothing else beyond s03" omits
  s08 (`thm:con-Pa`, `thm:con-Pb`) and (S3) through `thm:gam-lmi`. Write "(QRR) through `thm:con-glue`, (S3) through `thm:gam-lmi`,
  appendix A, and s03, s08".
* **M4. Theorem Gamma (iv) in s01.** `thm:int-gamma` states `c(z; 0) = Ginf(z)` as a headline result. Attribute the formula in the
  statement to [GJR18, Lemma 7.7, (37)] (DL L47: "Not new as a formula"; ST5b K7), so that the theorem does not read as a new formula
  (rule 11). The note's content there is the identification of its own series' value at `L = 0`.
* **M5. Title and abstract** (`design/TITLE_ABSTRACT.md`). (a) "every stable graph of genus g with at least one edge" -> "every
  leg-free stable graph ..." (as in s.0.2). (b) "the Landau-Ginzburg model x^4 + v/x" -> the note's normalisation `p^4/4 + v/p`
  (D13), or "x^4/4 + v/x". (c) "in every conifold frame" -> "in every conifold frame of a two-parameter rational family", since D14
  fixes that scope and CGYZ26 Rem. 1.4 admits a wider one (DL L8). (d) After R1, "the genus-0 mirror theorem" stays in the list.
* **M6. Evidence wording (rule 3).** s01 item 1 ("why the linear Hodge frames ... fail here"): DA A18 is VERIFIED at genus 2 only;
  write it as "computations at genus 2 (appendix C) show ...". s01 item 3 and `rmk:mai-graphlocal` ("a statement about the total cannot
  imply it", "no argument about total potentials can prove Theorem A"): the six singular genus-4 graphs are computed (DC1 UNI-7; MEC
  s.2); say "the computations of appendix C show that ...".
* **M7. `m5` in appendix B.** `sec:sgn-norm` prints `m5 = -43397/807975000`, a computed number (S01 C11; rule 3: only appendix C or a
  remark saying it is computed and unused). Keep the defining condition `p_(5q,0)(0) = 0` in B, move the value to C.
* **M8. Colour-basis graph sum.** The colour-basis edge `E_ij` and the leaves `Dil_i(d)` of `def:set-graphsum` (DA D9, C6-C7) are code
  conventions used only by the computations; the chain uses the PPZ form. Put the colour-basis form in appendix B or C; the
  `\draftnote` "derivation ... not written: DA C6" then leaves s02.
* **M9. "Programme family".** The outline says "programme family" throughout (s.0.2, rule 7, D14, `prop:gau-quintic`, `thm:con-Pb`).
  Fix the printed name ("the conifold frames of `def:set-frames`") and add "programme" to the words never printed (rule 5).
* **M10. `thm:gam-main` (ii).** "`FF` is `e^(SS_1)` without its `j = 0` factor" is RS2 item 32's shorthand; read literally it drops the
  Stirling series that the value formula needs. Spell out `FF = (5 sigma)^(1/2) prod_(j>=1) (zeta^j - sigma)^(-1/2) exp(Sum_k
  z^(2k-1) bk_k [(5 sigma)^(1-2k) + Sum_(j>=1) (zeta^j - sigma)^(1-2k)])`, the regular part of `A_St` (LF s.4.1, LA s.1.5).
* **M11. Appendix A.** (a) The purpose line ("the only non-formal input in the whole appendix is the exact Gaussian integral ...")
  contradicts `lem:fsp-refl` ("the branch part is not formal"); name both. (b) Add the `A[[z]]`-linearity of `SPn` in the amplitude to
  `def:fsp-wick` (LA-C2 = RS2 item 35, an input of LMI). (c) Map PB's and RL's "Lemma T(b)" to `lem:fsp-transl` (b) in NOTATION_2 s.9:
  "Lemma T" is also the note's name for `lem:tra-T`.
* **M12. D4 and W5.** "Route 2 as a second proof, so that Hock is optional" holds only once W5 is closed: `prop:mod-period` uses the
  closed formula on `x = z^4 + V/z` (XY s.2.3) exactly as `thm:mod-Lreg` does. Say so in D4 and at `thm:mod-gap`.
* **M13. Second-tier notation (NOTATION).** (a) `h` is both the genus of a vertex (`h_v`, rule (4)) and the series `h(z)`, `h_n`
  (table 2): use `g_v` for vertex genera, as the conifold note does. (b) `T` is the mirror map (table 2) and the translation `T(z) = z 1 -
  z R^-1(z) 1` (s.10 item 1): rename one. (c) `d(z)` in `lem:gau-wkb` against the dilaton factor `d` of s04: write `D_0` or a local
  letter. (d) The TR involution in `conv:mod-kernel` ("`Int_(sigma p)^p`") against the Mellin variable `sigma`. (e) `H = Q^5` (the
  state space, `conv:tra-forms`) against the class `H`. (f) `g_0 g_3 = 1` in `lem:gau-diag` uses an undefined `g_a` (GI; used only by
  GH s.2): define or drop. (g) Add to NOTATION_2 s.9: POLY's Lemma 1, 2 -> `lem:con-u0`, `lem:gau-poly`; GAM's "Theorem A(i)" never
  printed (it would clash with the note's Theorem A).
* **M14. Labels.** Assign the thirteen items of s02 to the registered subsections `sec:set-fq` .. `sec:set-N` (OUTLINE_3 names none);
  register the labels created by R1 (`prop:set-cyclicLR`, `cor:gau-rowsNN`).
* **M15. Scope in `thm:mod-var`.** "`2h - 2 + n >= 0`, `(h,n) != (0,0)`" includes `(1,0)`, which is EO (5-15) (`g >= 1`), not (5-14)
  (`g + k > 1`); state the two ranges as EO does.
* **M16. (Q) entry of DNOPS.** List (2-26)-(2-35) with Thm 8.1 in OUTLINE s.0.2 (as BIBPLAN does): the proofs use (2-32) in
  `prop:mod-Rtr` and (2-35) in `thm:mod-tr`.
* **M17. OUTLINE s.0.1.** "the quoted (Q) enter in s02, s04, s06, s07": add s10 (`cor:mai-gap`, PPZ and (S2)).
* **M18. (CYC) placement.** MS5 s.5 item 2 and ST5b s.8 item 2 (vi) say "(CYC) as a set-up lemma"; the design opens s03 with it. That
  is a fair reading; record it in OUTLINE s.0.4 so that the assembler does not move it into s02.
* **M19. Page budget.** 80 pp is inside 60-85, but s06 (DP, DR, LL, LL2, XY, TRM in 11 pp) and appendix A (FSL's five pieces, about 600
  dense lines, in 8 pp) are likely to run over by 3-5 pp; start `rmk:gam-sevenvar` as statement plus sketch ("not used"), which is
  already second in the drop order.
* **M20. BIBPLAN details.** The key `LPfq` breaks the initials-plus-year style of rule 1 (e.g. `LP20`, with the journal point kept
  open); add an entry for the pinned (S3) once chosen (M1).
* **M21. `design/SOURCE_ERRATA.md` additions:** DL L46 and DL s.2.2 paragraph 4 ("equivalent identity"; R4); DL s.2.2 paragraph 1
  (Fermat hypersurfaces under CGYZ26 Thm 2; R5); DA D2 ("the all-orders cyclic form ... is (CYC)"; R1); DC2 `thm:trm-ident` and TRM s.1
  (hypotheses in DNOPS's parameter `4v`; R3).

## 3. Spot checks (statement in the design against its source; arithmetic by hand or by a stdlib script, s.5)

| # | label | what was checked | source | result |
|---|---|---|---|---|
| 1 | `rmk:set-cg` | `c_2..c_5 = 1/1200, 1/25200, 1/180000, 1/660000` from `c_g` | DA A16; NF s.1 item 7 | correct |
| 2 | `rmk:set-qrrsign` | `h_1 = -3/20` (`Sum_(a in w_0) 1/a = 2 - 1/5`), weight `+5` gives `-11/60` | DA A5; MS5 s.2 | correct |
| 3 | `cor:gam-free` | `[z^5] Ginf = B_6(1/5)/6 = 4537/3937500` | DD2 GAM-10; LF s.4.3 | correct |
| 4 | `def:set-tau`, `lem:set-canon` | `tau zhat = -(5/4) zeta^i z/(u_i - u_i^c)`; `v^(4/5) = u^(4/5)/tau` from `v'^5/v = 1/(u(1-u)^5)` | DB D6-D7; UR s.1 | correct |
| 5 | `prop:gau-equalprod` | `prod b~ = -v'^5/v`; uniqueness given `v = u + O(u^2)`, not `v(0) = 0` | DB S12; NFS s.4; PS s.3 | correct, correction applied |
| 6 | `prop:gau-quintic` | `prod kappa = L^5`, `prod b = -1/(u(1-u)^5)`, `r = -1` (uses `varpi_0(0) = 1`, `D t_c(0) = -1`) | DB D5; RS2 items 27-28 | correct; see R1 for the input |
| 7 | `lem:gau-conc` | `Pf(CCmat) = -CC(f_0,f_2) CC(f_1,y_3) = 2`; `CC(f_0,y_3)` multiplies `0` | RS2 item 2; CF s.4.4 | correct |
| 8 | `lem:set-frobbasis` | obstruction matrix `[[-2,0,0],[17/5,0,0],[-7/5,-7/5,2]]`, kernel `Q(0,1,7/10)` | RS2 item 1 | identical |
| 9 | `thm:gau-cyc` | scope: any power-series `y_0` not proportional to `varpi_1`; `y_0(0) != 0` only for Theorem U | RS1 A1 (CA-K2) | correct |
| 10 | `cor:gau-cycsys` | rings named; log-free rows not part of it | RS2 item 8; RS1 A4 | correct |
| 11 | `lem:gau-vanq` (b) | pole of `b_1` meets `y^(n)_2` of order `(4n+4)/5 >= 4` | DB S24; GI s.2 | correct |
| 12 | `lem:gau-diag` | `p_(0,0) = tau^(-1/4)/varpi_0`; consistent with `d = tau^(1/4) varpi_0` and `(P_0)_00 (P_0)_33 = 1` | DB S27; S15 s.3 (I4) | consistent |
| 13 | `lem:gau-R1`, `thm:gau-G` | (R1) and (G) as statements; uniqueness clause flagged | DB S21, S28; GR s.2 | correct |
| 14 | `lem:tra-frames` | the five frames, the exponent `(4(k-j) + nu_a + 1)/5`; `frame(6,0)` from the values of `Pi_0` | DC1 D7, D9 | correct |
| 15 | `rmk:tra-consistency` | `<e>_1 = 5/24 - 1/12 = 1/8`, with `lem:set-dilaton`'s `5/24` | DC1 A5; S01 C3 | consistent |
| 16 | `prop:set-modrows` | `bmod = (1, -1/v, 1, -1, -1)` from (Mod1)-(Mod2) and `I_3 = v I_(-2)` | DB S18; UR s.1 | correct |
| 17 | `thm:mod-tr` | hypotheses of DNOPS Thm 8.1 in the note's `v` | TRM s.1-s.2.4; DC2 `thm:trm-ident` | wrong parameter: R3 |
| 18 | `thm:mod-gap`, `cor:mod-trF` | `(-5)^(1-h) B_2h/(2h(2h-2)) = c_h`; same formula for `x = W` in its own `v` | TRM s.0 item 3, l.153 | consistent |
| 19 | `thm:mod-V` (ii), `lem:mod-genus0` | `(F^TR_0)''' = 1/v = <(-log p)^3>_0`; `(F^TR_1)' = -1/(12v)` | LL s.2.4; LL2 s.3 | consistent |
| 20 | `thm:mod-lattice`, `prop:mod-LD` | statements, scope `Q[v]`, constants `varrho = (-1,1,-1,-3,-2)` | DC2 s.2.1-s.2.2; DP | identical |
| 21 | `def:gau-Fg`, `thm:flh-main` | (F_g) is `D_0 = 1 mod z^(3g-2)`; window `z^(3g-3)`; genus 3 sees `z^5` only | ST2 s.2 l.81-83; MS4 s.1 | correct |
| 22 | `thm:flh-main` | restated as `T_Gamma = F_Gamma` plus `G_Gamma = u^(2g-2) T_Gamma` | S15 s.4 display | faithful |
| 23 | `lem:con-form` (c) | `WWhat(1) = (5/4)/tau`, `WWhat''(1) = 5/(1-u)`, critical point `X = 1` | POLY s.5; PB s.2.1 | correct |
| 24 | `lem:con-cfu` | `V = tau^(-5/4)`, `Y_c = tau^(-1/4)`: `W_V(Y_c) = (5/4) V^(4/5) = WWhat(1)` | PB s.2.4 | consistent |
| 25 | OUTLINE s.0.4 item 12 | `e^(-WWhat(1)/zhat) = e^((u_0 - u_0^c)/z)`, so `e^(u_0/z)` only up to `e^(u_0^c/z)` | DD1 s.4.6; PS-C8 | correct |
| 26 | `thm:con-Pb` step 7 | `[u^1]` of `varpi_0 s(z) (Y^C_Q)_0` is `-s(z)/z` (orders of `lem:gau-vanq` (a)), so `s = 1` | PS-C7; PB s.2.9 | correct |
| 27 | `thm:con-glue` | `thm:gau-R0` (iv) at `u = 1` with `A^act_n(q=0) = h_n` | PS s.4.2 (e); PB s.3 | correct |
| 28 | `prop:con-deg`, `cor:gam-poly` | `deg ph_n <= floor(4n/5)` iff `A_n` in `L Q[L]`; `[u^(4m)] ph_(5m) = c_(5m)(0)` | ER s.4; LF s.4.3 | consistent |
| 29 | `thm:gam-main` (ii) | LMI at `s_* = 1`, `a_* = -1`, `M_* = 0`; `exp(Sum bk_k T_k z^(2k-1)) = h_(w_0)` | LF s.4.1; LA s.1.5 | correct (see M10) |
| 30 | `thm:gam-main` (iii) | exponents `= 4n mod 5`, `c_n(zeta L) = zeta^-n c_n(L)` | CP3 s.3.9 | consistent |
| 31 | `thm:gam-ginf` | `s_* = 0`, `a_* = 5`, `M_* = 1`, `delta = -z/5`, `Q = (-q)^(-1/5)` | LF s.4.2; LF-C5 | consistent |
| 32 | `thm:gam-lmi` | hypotheses LF-C1..C4, LA-C1..C6 (punctured disc, `FF(s_*,0) != 0`, linear `a`) | RS2 items 29-37 | applied |
| 33 | `thm:set-qrr` | `R(q=0) = C_alpha^-1 = h_(w_alpha)` from GJR18 (17), weight `-5 h_alpha` | DL L27 | correct |
| 34 | `lem:set-pf5`, `def:set-genus0` | `5q prod_(k=1..4)(5D+k)` gives DA D2's coefficients; `PF5 = z^5 D PF4 - 1` | DB D2; CF2 s.8.2 | correct |
| 35 | `def:set-rows` | `Q^(i)_(a,n) = zeta^(i(a+1-n)) Q^(0)_(a,n)` from `P^(i)_a(z) = zeta^(ia) P^(0)_a(zeta^-i z)` | DA D6; CS s.1 | consistent |
| 36 | `prop:gau-quintic` input | the forward relations of the actual rows | UR s.1; CA s.2.1; FQ2 s.5 | missing: R1 |
| 37 | s01 positioning | GJR18 Rem. 7.2 wording; CGYZ26 Thm 2 scope | ST5b s.8; DL L8 | R4, R5 |
| 38 | `cor:mai-gap` | composition of R-actions | FR s.0; PPZ v3 Props 2.4, 2.9, 2.14 | M2 |

## 4. The eight checks of the task, item by item

1. **The chain closes, without circularity.** Order: s02 (definitions; quoted (S1), (S2), (QRR)) -> s03 ((CYC), Theorem U, (R1), the
   vanishing solutions, (G), symplecticity, HPF, polynomiality, WKB uniqueness, R0) -> s04 (transport, from s03) -> s05 (Theorem R, ML
   as hypothesis) -> s06 (Model Lemma, nothing about the quintic) -> s07 (comparison, (F_g) as hypothesis) -> s08 (P(a) from
   `lem:gau-poly` and appendix A, no Theorem U; P(b) from s03 and appendix A; the glue from `thm:gau-R0` (iv) and (QRR)) -> s09 (LMI;
   `prop:gam-hatflat` uses `thm:con-Pa`, which does not use it, DD1 s.7 item 7; `thm:gam-D`) -> s10. The only forward use is
   `def:set-confexp` -> `thm:gau-R0` (iii), removed by R1 item 3. Inputs: complete after R1, with R7's list; nothing superfluous (Lemma
   N only in `sec:con-finite`; DLMF only for `Ginf`; Mumford, MR1728879, MR1666787 only in remarks and appendix C; Hock optional, M12;
   (S3) used by LMI, M1). The items named by the task are present: (CYC) C2; the formal lemmas (inputs of `thm:con-Pa`, `thm:con-Pb`,
   `prop:gam-hatflat`, `thm:gam-lmi`); Lemma H (`lem:fsp-homog`, in `thm:con-Pa`, `lem:con-S2`, `prop:gam-hatflat`); formal WKB
   uniqueness (`lem:gau-wkb`, in `thm:gau-R0` (iii)-(iv)); the ring-of-comparison sentence (`thm:gau-R0`, PS-C8).
2. **Final corrected forms.** 38 spot checks (s.3): 33 correct, consistent or faithful; one wrong parameter (R3), one missing input
   (R1), two positioning sentences (R4, R5), one citation (M2). The adjudications are applied wherever I looked: RS2 items 1, 2, 8, 13,
   24-37; RS1 A1, A4, A6; ST5b K3, K8, K9, K12, O8, O10; FR C1-C6 via ST2b O14-O20; NFS s.4 (`v` uniqueness, Theorem D bound only).
3. **TR only on the model curve.** Every TR statement is in s06 and on `x = W_v`, `y = p` (DNOPS at `f = 4 W_v`, R3; ABDKS on `x = z^4 +
   V/z`; Hock and EO Thms 4.4, 4.7, 5.1 on `x = W_v`). Family statements (TRM s.4.1's `P(z) + v/(z - a)`, XY's four-parameter family,
   TR-5's Hurwitz space) are remarks that no proof uses; `rmk:mod-residue` is recommended omitted. `rmk:mod-H2` cites CvdLPS23 Thm 2,
   proved by TR on `x = z^4 + 1/z`, the model curve. FCCZ26 is context only. No TR on a mirror curve of the (formal) quintic.
4. **Positioning and novelty.** No novelty is claimed for the formal lemmas, `Ginf` or the chart mechanism, and nothing about the true
   quintic (DL s.2.3 items 1-10 are carried into `sec:int-notclaimed`). To fix: R4, R5; attribution in `thm:int-gamma` (M4); the
   abstract's scope (M5).
5. **Notation.** NOTATION resolves the clashes it lists, but the outline's chain formulas still use renamed source letters (R8 (a)-(f),
   (i)), and NOTATION lacks the generating variable, the differential ring, the conifold-expansion map, `log q`, `omega`, `omega^`
   (R6, R8 (g)-(h)). Second tier: M13.
6. **Labels.** A stdlib script (s.5) finds every `type:prefix-name` label of OUTLINE*.md in LABELS*.md, except three digest labels used
   as source pointers (`prop:pol-Ginf`, `rmk:mlk-dilaton`, `rmk:thp-Ai`), which are not note labels. Complete for the outline's
   cross-references; add the labels of R1 and the s02 subsection map (M14).
7. **Binding rules.** No stage, workflow, agent, verdict, status tag or file name is planned for the `.tex` files; evidence only in
   appendix C; errata only in `design/SOURCE_ERRATA.md`. To fix: the computed `m5` in appendix B (M7); evidence-based sentences in s01
   and `rmk:mai-graphlocal` (M6); FQ2 cited for an all-orders fact (R1 item 3); the word "programme" (M9).
8. **Page budget and pointers.** 80 pp, balanced; s06 and appendix A tight (M19). Every statement carries a digest item, a source file
   and section and a mark [W]/[W-c]/[S]/[M]/[V]. Two places send the writer to sources that no digest read (`def:set-modcohft`: LG s.3,
   UR s.1, MSG1 s.1; `cor:mai-natfact`: UR s.4, OU3 s.3); both carry a `\draftnote`.

Decisions D1-D15: I agree with the recommended options of D1, D3, D4 (with M12), D9, D13 (with R3's scaling sentence), D14 and D15;
D2 grows with R1; D5, D6, D7 stay with the authors and belong to R7's list.

## 5. My own checks (nothing in the project was modified)
* Eynard-Orantin and PPZ: read in the designer's text extractions of the arXiv PDFs (header "arXiv:math-ph/0702045v4 7 May 2007").
  Theorem 4.4 = (4-23)-(4-24); Corollary 4.1 = (4-26)-(4-29); Theorem 4.7 = (4-30); Theorem 5.1 = (5-14) ("For g + k > 1") and (5-15)
  ("for g >= 1"), just before s.5.2; s.5.3 "Variations with respect to the moduli" holds the temperatures (5-27)-(5-30) and later
  (5-37)-(5-38); (5-13) with footnote 8 (hard edges excluded). PPZ v3: Remark 2.1, Definition 2.2 ("We use the inverse of the R-matrix
  in all of our formulas in Definition 2.2"), Propositions 2.3, 2.4, Definition 2.5, Propositions 2.7-2.9, 2.12, Definition 2.13,
  Proposition 2.14. OUTLINE_2 s.0.4 items 1-2 and SOURCE_ERRATA item 1 are right; DL L28 is wrong where they say.
* Arithmetic (stdlib, exact fractions): `c_2..c_5`; `B_6(1/5)/6 = 4537/3937500`; `h_1 = -3/20`, `-11/60`; `WWhat(1) = (5/4)/tau` to `u^7`.
* Labels: every `thm|prop|lem|cor|def|rmk|conv|eq|sec|tab:prefix-name` of OUTLINE*.md compared with LABELS*.md (s.4 item 6).
* Root bibliography: the 16 root keys of BIBPLAN each occur once in `refGW/refs.bib`; none of the new-key papers (PPZ, Givental 2001
  IMRN, GJR, LP, ZZ, CGYZ, Brini 2025, HKQ, HK, Hock, ABDKS, Coates-Iritani 2021, ...) is there; `MR1653024` is (R1 item 2).
* Sources read beyond the digests: MS5 (all); UR s.1-s.4; FQ2 s.1, s.5; CA s.0-s.3; CF2 s.8.2; ST5b s.5-s.8; RS1 s.2-s.3; RS2 s.5-s.6;
  S15 s.4-s.5; TRM (headings, s.0, s.2.1-s.2.4 lines on the normalisation); LIT s.1.5, LIT_2 s.3.6, LIT_3 s.4; LA s.0, s.2.1; FR s.0;
  ST2 s.2 and ST2b s.6.3, ST4b s.6.3 (the write-up plans); the conifold note's OUTLINE s.0, LABELS header and `\draftnote` macro.
