# Review of the note design (key design_critic), part 1: verdict and required changes

Author-facing; not part of the note; nothing here enters a `.tex` file. Object, read in full: `NOTE/OUTLINE.md` .. `OUTLINE_9.md`,
`NOTATION.md` (+`_2`), `LABELS.md` (+`_2`), `BIBPLAN.md` (+`_2`), `design/TITLE_ABSTRACT.md`, `design/SOURCE_ERRATA.md`; checked against
the digests `design/DIGEST_*.md` and, where named, the sources (abbreviations of `OUTLINE_2.md` s.0.3; digests DA, DB, DC1, DC2, DD1,
DD2, DL as in `OUTLINE.md`). No design file was edited. Part 2 (`OUTLINE_REVIEW_2.md`): recommended changes (s.2), the spot-check
record (s.3), the eight checks of the task item by item (s.4), my own checks (s.5).

## 0. Verdict: **ready-after-required-changes**

* The chain C1-C13 of OUTLINE s.0.2 is the chain of MS5 s.2, ST5b s.2 and RS2 s.6.2 in final corrected form. Every input that an
  earlier stage had missed is at its place: (CYC) `thm:gau-cyc` with the scope of RS1 A1; the formal lemmas of appendix A; Lemma H
  `lem:fsp-homog`; formal WKB uniqueness `lem:gau-wkb`; the ring-of-comparison sentence and "the same exponential" in `thm:gau-R0`
  (PS-C8); the model identity `lem:con-S2`, the relabelling and the key lemma; L3(iv) `lem:fsp-stirfe`; FSL's Lemma T(b)
  `lem:fsp-transl` (b). Lemma N's `K_q` formula appears only in the optional `sec:con-finite` (ST5b K3). No circularity (s.4 item 1).
* Topological recursion is confined to s06 and to the model curve; family statements are remarks that no proof uses (s.4 item 3).
* BIBPLAN reuses every root key correctly: all 16 root keys it names exist in `refGW/refs.bib`, and none of the new-key papers is
  there. The designer's corrections of the Eynard-Orantin and PPZ numbers are right (s.5).
* Eight required changes (s.1). R1 is a genuine missing link of the chain: the flat-basis quantum product at large radius, which
  no source derives and which (S1) as worded does not supply. R2-R8 are transcription, scope, wording and notation corrections;
  none changes a statement of the sources. After R1-R8 the design is ready for the writers.

## 1. Required changes

### R1. The large-radius cyclic structure of the rows is missing from the chain
**Where.** OUTLINE s.0.2 (row C1; bullet (S1)); OUTLINE_3 s02 item 2 (`thm:set-mirror`), item 5 (`def:set-rows`), item 8
(`def:set-confexp`: "rows, edges, leaves and graph terms lie in `NN` (FQ2 s.1 "Rows"; FQ2 s.4)"); OUTLINE_4 `prop:gau-quintic`
("the rows satisfy `zD Y_a = kappa_a Y_(a+1)` ... (by `cor:gau-cycsys`)").

**Why.**
* `cor:gau-cycsys` (DB S10) gives the cyclic system iff `y_0 Y_0` solves PF5 AND the four forward relations `Y_(a+1) = zD Y_a/kappa_a`,
  `a = 0..3`, hold. (CYC) supplies only the closing relation `4 -> 0` for rows built by those four relations (RS2 s.5 item 13 =
  CA-K7; DB S10: "PF5 => closing relation ... for the actual rows (UR s.1, s.4)").
* For the actual R-matrix the forward relations are the flat-basis quantum product `H * H^a = kappa_a H^(a+1)`, `kappa^LR = (1,
  L^5/(C_0^2 C_1^3), 1, C_0/C_1, C_0/C_1)` (DA D2), i.e. `X_0 = psi_i/C_0`, `X_(a+1) = (1/(C_1 kappa_a))(zD + zeta^i L) X_a`.
* No source derives it. UR s.1 and CA s.2.1 read it off the construction of the project's rows ("fq_rows.Rows builds `X_0 =
  psi_i/c_0` and `X_(a+1) = kinv[a] (zD + zeta^i L) X_a`"); DB D5 records it as the definition of the rows; FQ2 s.5 lists
  "kappa_3 = kappa_4 = C_0/C_1 beyond q^12 evidence" among the points not checked (VERIFIED to `q^12` only; FQ2 s.1 "Tier-1
  agreement with the I-function q-series to q^12"). DA D2's "the all-orders cyclic form used by the programme is (CYC)" conflates
  the closing relation (PROVED) with the forward coefficients (not proved). `kappa_3` is a genuinely equivariant 3-point
  invariant, not the quintic's.
* The design's (S1) ("the quantum differential equation of the formal quintic is PF5 and its fundamental solution is given by the
  I-function after the mirror map") concerns the scalar equation and does not give `kappa_1`, `kappa_3`.
* The all-orders membership of rows, edges, leaves and graph terms in the differential ring (needed to expand `R` in a frame in
  `def:set-outer`, and for `cor:mai-gap`) is attributed to FQ2 s.1, s.4, a genus-2 computation (rule 3). In the sources it follows
  from `A^act_n` in `Q[L, 1/L]` (PS s.4.2 (c); RS1 A4-A5) together with the forward relations. The design proves the first in
  `thm:gau-R0` (iii), in s03: so `def:set-confexp` also has a forward dependency.
* The conifold note has the analogous fact as a proved proposition (its `LABELS.md`, s02 row "Frobenius data in the flat basis: eta,
  the product ... (proof from the genus-0 potential)").

**Fix.**
1. Add to s02, after `thm:set-mirror`, a statement `prop:set-cyclicLR` (register it in `LABELS.md`): in the flat basis
   `H * H^a = kappa_a H^(a+1)` with `kappa^LR` as above; equivalently the rows of `def:set-rows` satisfy `X_0 = psi_i/C_0`,
   `X_(a+1) = (1/(C_1 kappa_a))(zD + zeta^i L) X_a`, `a = 0..3`. Statement from DA D2, D6; DB D5; UR s.1. Proof:
   `\draftnote{the sources build the rows by this recursion (UR s.1; CA s.2.1) and check kappa_3 = kappa_4 = C_0/C_1 only to q^12
   (FQ2 s.1, s.5); a derivation from the genus-0 mirror theorem, or a pinned citation stating the flat-basis product, is needed
   (decision D2)}`. A possible route, for `DRAFT_NOTES.md` only (it is in no source; rule 1): the Z/5 grading of the specialised theory
   (`<H^a, H^b, H^c>_(0,d)` has lambda-degree `a + b + c - 3`, and at `lambda_i = zeta^i lambda` only `e_5` of the symmetric functions
   survives), the unit (`kappa_0 = 1`), the Frobenius property (`kappa_0 = kappa_2`, `kappa_3 = kappa_4`), the Yukawa coupling
   `5 L^5/(C_0^2 C_1^3)` from (S1), and `prod_a C_1 kappa_a = L^5` (the symbol; CA s.2.1).
2. Rewrite (S1) where it is printed (OUTLINE s.0.2, `sec:int-inputs`, `thm:set-mirror`, TITLE_ABSTRACT) as what is used: the
   J-function of the formal quintic is `I/C_0` after the mirror map `T = I_1/C_0`, and, if the pin chosen in D2 states it, the
   flat-basis product of `prop:set-cyclicLR`. Extend D2 to cover both; add the root key MR1653024 (DA s.7 item 2) to its candidates.
3. In `def:set-confexp` replace "rows, edges, leaves and graph terms lie in `NN` (FQ2 s.1, s.4)" by a forward reference to a new
   corollary after `thm:gau-R0` (`cor:gau-rowsNN`, register it): by `thm:gau-R0` (iii) and `prop:set-cyclicLR` the rows lie in the
   differential ring, hence so do edges, leaves and graph terms. Sources RS1 A4-A5, PS s.4.2 (c). FQ2 s.4 goes to appendix C.
4. Add to the proof plan of `prop:gau-quintic` the step that the conifold expansion sends `C_1 kappa_a^LR` to `kappa_a` (`a = 0..3`)
   and `C_1` to `D t_c` (CA s.2.2, with its hand computation of the image of `W = c_0^2 C_1`); the frame rows then satisfy the four
   forward relations, and the closing relation is `cor:gau-cycsys`.
5. Add a row C1' (`prop:set-cyclicLR`; input (S1); `\draftnote`) to the table of OUTLINE s.0.2 and to the list of R7.

### R2. The printed input names (S1)-(S3) collide with source labels
**Where.** OUTLINE s.0.2 ("(S) standard results cited: (S1) .., (S2) .., (S3) .."), s.0.1, OUTLINE_3 s01 item 4, `LABELS.md` l.12,
`NOTATION_2.md` s.9, TITLE_ABSTRACT "Checks".

**Why.** The writers transcribe sources in which "(S2)" is the model identity (PB s.2.6; PS s.1; DD1 `lem:thp-S2`; the design's own
`lem:con-S2`, "source "(S2)"", and "the residue identity, (S2), the uniformising map" in `sec:evi-normal`), "(S1)" is CP's ring
lemma (CP s.0 table "(S1), (PG), (L0), (L1a) .."; RS2 s.6.1), and S1-S9 is ST5's chain. NOTATION_2 s.9 maps the source "(S2)" to
`lem:con-S2` and declares "S1..S9" never printed, while s.0.2 prints (S1)-(S3) with other meanings. Rule 6.

**Fix.** Rename the three standard inputs with letters no source uses, e.g. (M) the genus-0 mirror theorem, (T) Teleman's
classification, (IFT) the holomorphic implicit function theorem, in every file above. Add to NOTATION_2 s.9: CP's "(S1)" ->
`def:fsp-ring`/`lem:fsp-univ`, the sources' "(S2)" -> `lem:con-S2`, both never printed. In `sec:evi-normal` write `lem:con-S2`.

### R3. The quoted DNOPS theorem is stated at the wrong point of the Hurwitz space
**Where.** OUTLINE s.0.2 (Q) bullet ("[Dunin-Barkowski:2015caa] Thm 8.1, superpotential (8-1), `n = 4`, `a = (0,0,0,0,v)`");
OUTLINE_3 `def:set-modcohft` ("(8-1), `n = 4`, at `a = (0,0,0,0,v)`"); OUTLINE_6 `thm:mod-tr` ("`n = 4`, `a = (0,0,0,0,v)`; critical
points `p^5 = v/4`, critical values `5 p_c^4`").

**Why.** With the note's normalisation `x = W_v(p) = p^4/4 + v/p` (D13; `def:mod-curve`: branch points `p^5 = v`), DNOPS (8-1) with
`n = 4` is `f = p^4 + a_2 p^2 + a_3 p + a_4 + a_5/(p - a_1)`, so `x = f/4` with `a = (0,0,0,0,4v)`. TRM s.2.1: "kappa = 4: x = lambda =
p^4 + v/p with v = 4u"; TRM s.2.3: "R_TR(z) = M_0(z/kappa)"; TRM s.2.4: factor `(-5 kappa)^(g-1)` (TRM l.110-114: `(-20)^(g-1)` and
`4^-D` at kappa = 4). The hypothesis check copied from TRM s.1 and DC2 `thm:trm-ident` (`p^5 = v/4`) is in DNOPS's parameter
`v_D = 4v`. Read with the note's `v` it is false by a factor 4 (rule 7).

**Fix.** Quote "(8-1) with `n = 4` at `a = (0,0,0,0,4v)`, i.e. `f = 4 W_v`; critical points `p^5 = v`, critical values `f(p_c) =
5 p_c^4`, pairwise distinct" in all three places, and add to the proof plan of `thm:mod-tr` the transfer from `x = 4 W_v` to
`x = W_v` (TRM s.2.1-s.2.4: `R_TR(z) = M_0(z/kappa)`, `Lambda = -5 kappa`, kappa = 4 -> 1). It is the DNOPS case of DC2 s.7 W5,
written on the CohFT side in TRM s.2.4, so no new `\draftnote` is needed. Keep "the same curve up to scaling" (rule 8), naming it.

### R4. GJR18 Remark 7.2 must not be called "an equivalent identity"
**Where.** OUTLINE_3 s01 item 5 ("[GJR18, Rem. 7.2 states an equivalent identity without proof]"); s01 item 2, which transcribes DL
s.2.2 paragraph 4 ("an equivalent identity for the Landau-Ginzburg integral of the formal quintic is stated without proof in [GJR18,
Remark 7.2]"); OUTLINE_8 `rmk:gam-prior` (to be worded the same way).

**Why.** The identification of GJR's Landau-Ginzburg saddle series with `c(z; L)` is "written neither in GJR nor in the project"
(LIT_2 s.3.6). ST5b s.8 item 2 and MS5 s.3.6: "the analogue of (Gamma3) for their LG integral"; DD2 s.7 item 7: "Do not claim that
GJR's LG series equals c(z; L)"; the design's own `cor:gam-free` says the note does not claim it. "Equivalent" asserts an equivalence
that the note does not prove (rules 1, 11). The word comes from DL L46 and LIT_3 s.4 item 3.

**Fix.** Write "the analogous identity for their Landau-Ginzburg integral is stated without proof in [GJR18, Remark 7.2]" in s01
items 2 and 5 and in `rmk:gam-prior`. Add DL L46 and DL s.2.2 paragraph 4 to `design/SOURCE_ERRATA.md`.

### R5. CGYZ26 Theorem 2 concerns the quintic only
**Where.** OUTLINE_3 s01 item 2 ("Positioning paragraphs 1-2 of DL s.2.2, transcribed"). DL s.2.2 paragraph 1 ends "... and for the
quintic and the other one-parameter Fermat hypersurfaces in [CGYZ26, Theorem 2], by master-space localisation".

**Why.** DL L8: Theorem 2 is the quintic ("the conifold gap property holds for the quintic threefold"); the Fermat hypersurfaces of
degree 6, 8, 10 are extensions in CGYZ26 s.5, which DL did not read (DL s.0.2: "Not read: s.2 beyond Thm 3, s.5"). That cites an
unread statement under the number of another (rule 4).

**Fix.** Transcribe paragraph 1 with "for the quintic in [CGYZ26, Theorem 2]"; mention the other targets only as "see also [CGYZ26,
s.5]" with `\draftnote{CGYZ26 s.5 not read}`, or drop them. Record in `design/SOURCE_ERRATA.md`.

### R6. The statement of Theorem A must name the unhatted outer term
**Where.** OUTLINE s.0.2 (main theorem: "In the outer graph sum `R'.(M.omega^)` of genus `g`: (i) the term `T_Gamma` ..."); OUTLINE_3
`def:set-outer`; OUTLINE_8 `thm:mai-main`; s01 `thm:int-main`.

**Why.** In the sources `omega^` is the hatted TFT (DA D13: "`M.omega^` (`omega^` the hatted TFT)"), and the term of the hatted outer
sum is `G_Gamma = u^(2g-2) T_Gamma`; `T_Gamma` is the genuine term (DA D18; S15 s.4: `T_Gamma = Cont[E_R', q_R'; M_nat.omega^Q]` with the
quintic's TFT `omega^Q`, then `G_Gamma = u^(2g-2) T_Gamma` by its Lemma 3). "The term of `R'.(M.omega^)` lies in `Q[[u]]`" is a weaker
statement than "`T_Gamma` lies in `Q[[u]]`". NOTATION defines neither `omega` nor `omega^`. Nor does `def:set-outer` say which matrix
`R` is: the comparison needs the flat-basis rows (`R' = Q M^-1`, `R'_0 = P_0 Psi^m`, a non-unital constant part, `def:flh-outer`;
S15 s.4), not the canonical-frame `R = 1 + O(z)` of `def:set-rows`.

**Fix.** In all four places: `T_Gamma := Cont_Gamma[E_(R'), q_(R'); M.omega]`, `omega` the TFT of `def:set-frob`, `R' := Q M^-1` with
`Q` the matrix of the lowered rows of `def:set-rows` expanded in the frame; `G_Gamma` the same contraction with the hatted data of
`def:set-hat`; `G_Gamma = u^(2g-2) T_Gamma` (`lem:flh-hat`). Add `omega` and `omega^` to NOTATION tables 2 and 8.

### R7. The conditionality sentence needs the list of the open steps of the chain
**Where.** OUTLINE s.0.2 "Exact final list of inputs"; s.0 rule 8; s01 `sec:int-inputs`; s10 `thm:mai-main` ("Close with the
conditionality paragraph").

**Why.** Rule 8 requires the conditionality to be stated exactly; rule 1 makes every missing step a `\draftnote`. The design leaves
about twenty chain steps as `\draftnote`s, scattered over nine files. While they are open, "conditional only on (Q), (QRR) and
(S1)-(S3)" is not what the note proves, and the assembler has no list to check it against.

**Fix.** Add to OUTLINE s.0.2 a table "steps of the chain left as `\draftnote`" and require `sec:int-inputs` and the closing paragraph
of `thm:mai-main` to carry one `\draftnote` naming them by label (deleted once they are closed). From the design: `thm:set-qrr` (a)
the passage to `A^QRR_n(q = 0) = h_n`, (b) the uniqueness step; `prop:set-cyclicLR` (R1); `lem:set-dilaton` (no pinned citation);
`def:set-M0` symplecticity and `def:set-modcohft` homogeneity (if not written in LG s.3, UR s.1, MSG1 s.1); `thm:gau-U` steps 3, 5, 8,
10, 11, 12 (G1-G5, G7; D7, D15); `prop:gau-graded` step 1 [M]; `lem:gau-diag` step 4 (G8); `thm:gau-G` Cor. (a) (G9);
`cor:gau-colours` (G6); `lem:tra-frames` step 4 [S]; `lem:tra-conv` (every item "argued" in the sources); `lem:tra-E` step 2 (O-8);
`prop:mod-Rtr` (D6, G2); `lem:mod-xfun` (G3); `thm:mod-Lreg` (W5; also route 2 of the Gap Theorem); `thm:mod-V` base `h = 1` (D5, G1);
`thm:con-Pa` step 5 (the domain). Outside the chain of Theorem A, listed apart: `cor:mai-gap` (D3), `cor:mai-natfact` (`Psi(u)`),
`prop:con-deg` step 5, the remarks.

### R8. Notation clashes inside chain statements and proofs (rule 6)
The outline's own formulas use source letters that NOTATION renames, or letters that clash with note symbols. Fix each:
* (a) s04 `lem:tra-T` (c)-(d), `prop:tra-can`, appendix C `sec:evi-transport`: `s` for `d T_P - v` -> `delta_P` (NOTATION_2 table 6);
  the scalar in `def:tra-slice` ("`Pi_0(-z 1 + n + s H) = v + s`") -> a declared local letter (clash with `s = u^(1/5)`).
* (b) s07 `thm:flh-main`: "`c_10 = c_5^2/2`" (`c_5 = 1/660000` is the genus-5 constant) -> "`[z^10] D_0 = ([z^5] D_0)^2/2`".
* (c) s03 `lem:gau-HPF`: `A_e` and `Phi_j := A_(-j) Phi_(j-1)` clash with the flat amplitude `A_n` and the colour-`i` WKB function
  `Phi_i` of `def:set-rows` -> new letters, registered in NOTATION table 2.
* (d) s06 `def:mod-curve`, `def:mod-pair`, `conv:mod-kernel`, `thm:mod-tr` ("`Delta^TR_j = Lambda c_j^2`, `Lambda = -5`"),
  `lem:mod-genus0`: branch points `c_b`, `c_j` -> `varsigma_b` (NOTATION_2 table 7), and `Lambda` -> the number `-5` (NOTATION_2 s.9;
  `Lambda` is the leg map of `def:mod-weights` in the same section); `prop:mod-dii` "`Lambda(s)` ... `s` in `K_0`" -> another letter.
* (e) s08 `prop:con-res` ("`VV = C_0 (P^SP Y^C)_0` in every frame") and `thm:con-Pb` step 7 ("`C_0(0) = 1`"): `C_0` -> `varpi_0`
  (NOTATION table 2: `C_0` only at large radius); `prop:con-deg` "`psi_n` ... `x = L`" -> `A_n`, `L`; `lem:con-relab` `HH` (the loop
  space `\cH` of s04) -> e.g. `\hat{\cV}`; `lem:con-red` `beta(Y)` (frame parameter) and the exact amplitude `Gamma` -> local letters.
* (f) s09 `def:gam-c`: `P_k(L)` -> the local `p_k(L)` of NOTATION_2 s.9.
* (g) s02 `def:set-fq`: `F_g = Sum N_(g,d) Q^d`; NOTATION has no letter for the generating variable, and `Q` is already `\bbQ`, `\Qr`
  and the sans `Q` of LMI -> add a symbol (with the `\draftnote` of `def:set-fq`).
* (h) `NN` is both the differential ring of `def:set-confexp` and DP's module `\mathsf N` of `def:mod-weights`; NOTATION has no symbol
  for the ring, for the conifold-expansion map (its source name is a code name, rule 5) or for the formal `log q` ("`Lam`") -> add.
* (i) s07 `def:flh-outer`: the TFT-rescaling scalar `kappa` clashes with `kappa_a` of s03 -> e.g. `\varkappa`, in NOTATION table 8.
