# Review rf1: `s01_intro*.tex`, `s10_main*.tex` (part 2 of 2)

Continuation of `REVIEW.md` (verdict s.0, checks s.1, issues S01-1..S01-12 and S10-1..S10-19).

## 4. The assembly of Theorem A against the chain of OUTLINE s.0.2

The proof of `thm:mai-main` cites C12, C10, C9, C8, C7 directly; the new dependency paragraph (S10-17) names C11, C6, C5, C4,
C3, C2, C1'. Each row: the statement cited, its hypotheses as printed, and whether they match at the place of use.

| step | statement used | hypotheses as printed | match |
|---|---|---|---|
| C1 | `def:set-rows` (rows of the CohFT R-matrix, Givental's convention), `thm:set-qrr` (QRR), `thm:set-GT` (T) | quoted | yes: `T_Gamma` is built from the rows of `def:set-rows`; (T) identifies that R |
| C1' | `prop:set-cyclicLR`, `cor:set-cyclicPF5` | (M), (G0) | yes: used through `prop:gau-quintic`, `cor:gau-rowsNN`, `thm:gau-R0`, `lem:mai-Psi` |
| C2 | `thm:gau-cyc`, `cor:gau-cycsys` | every frame | yes (closing relation in `prop:gau-quintic`) |
| C3 | `prop:gau-quintic`, `prop:gau-equalprod`, `thm:gau-U` | one-pole shape, (EP), char 0 | yes: (ii) of both props for `v = u tau^(-5/4)` and `(alpha, beta)` in `Q^2` |
| C4 | `lem:gau-R1`, `lem:gau-diag` (N0), `thm:gau-G`, `def:gau-can` | conifold frame | yes; (N0) holds in the frames of `def:set-frames` |
| C5 | `prop:gau-symp` | self-dual, (EP), `r = r~` | yes (`r = r~ = -1`) |
| C6 | `thm:gau-R0` (iv), `cor:gau-colours`, `cor:gau-rowsNN` | conifold frame | yes; colours `i != 0` use (QRR) through `lem:gau-colrows` |
| C7 | `prop:tra-can` | `(alpha, beta)` in `Q^2`, `v(0) = 0`, `v'(0) = 1` | yes: dressing with (Sp), (A), `T_Pcan = t_c`, `d(0) = 1` |
| C8 | `cor:red-can` (from `thm:red-main`) | `ML(h)`, `h <= g-1` (edges); `ML(g)` (vertex) | yes, supplied by C9 |
| C9 | `thm:mod-MLK` | (Q) only | yes: `ML(h)` for every `h >= 0` |
| C10 | `thm:flh-main` | `g >= 2`, (F_g), frame fixed, `M` exact in `u` | yes: (F_g) from C12; same `T_Gamma` (`def:set-outer`) |
| C11 | `thm:con-Pa`, `thm:con-Pb` (every conifold frame), `thm:con-glue` | (QRR) in the glue | yes |
| C12 | `thm:gam-main` (ii), `prop:gam-hatflat`, `thm:gam-D` | appendix A | yes: `D_0 = 1`, `D = 1`, (F_g) for all `g` |
| C13 | `thm:mai-main` | none new | yes |

No circularity: `thm:gam-D` uses neither `R'` nor `T_Gamma` (stated in `def:set-outer`); `thm:con-Pa` (iv) does not use
`prop:gam-hatflat` or `thm:gam-main` (stated in s09); `lem:mai-Psi` uses only `prop:set-cyclicLR`, `cor:gau-rowsNN`,
`lem:gau-confexp`, `def:set-frob`, `def:set-M0`. The frame: `D_0` is defined per frame (`thm:gau-R0`), the glue and Theorem P
hold in every frame, so `thm:gam-D` holds in every frame (and `cor:con-frame` makes `D_0` frame independent); Theorem A's proof
fixes the frame first. Regularity: `F_Gamma` lies in `Q[[u]]` (rational coefficients), hence so does `T_Gamma = F_Gamma`, which
a priori lies in `Q(zeta)((s))`. No proof in s10 cites appendix C; the remarks that mention it say that the computation is
not used.

## 5. Corollary B, its constant, `rmk:mai-constant`

* `def:mai-Fcf`: ring statements checked: `F_g^NN` is a polynomial over `Q` in finitely many row coefficients, `Delta_i^(+-1/2)`
  and `zeta` (now with the symmetric edge, S10-2); Galois descent (`sigma_k(Q^(i)_(a,n)) = Q^(ki)_(a,n)`,
  `sigma_k(Delta_i^(1/2)) = Delta_(ki)^(1/2)`, `sigma_k(R~_ij) = R~_(ki,kj)`, colourings permuted) gives `F_g^NN` in `NN`;
  `ev(F_g^NN) = F_g` by (T), R symplectic over `Q[[q]]`; independence of the primitive of 1 (G5R F2.1); no double `|Aut|`
  (G5R F2.2).
* Proof of `cor:mai-gap`: Step 1 (ring map; `Cexp(Q) = Pcan Q^m` from `D = 1`; `Cexp(R~)` symplectic from `prop:gau-symp`,
  `lem:set-M0sympl`, the `Lambda_varsigma` line and `lem:mai-Psi`; the edge is PPZ's), Step 2 (`X = Cexp(R~) M^-1`, `X(0) = 1`,
  symplectic; PPZ Prop. 2.14 for `Cexp(R~) = X M`; unit by Prop. 2.12; `E_(R') = E_X`, `q_(R') = q_X` because `Psi_cf` is a
  constant isometry mapping `1_TFT` to 1), Steps 3-4 (Theorem A, finitely many graphs): correct, and as in G5_2 s.2.4 with
  G5R s.2.3. Corollary B therefore holds exactly in the frames of Theorem A, with the conditionality of Theorem A.
* The constant: `c_g = (-1)^(g-1) B_2g/(2g(2g-2) 5^(g-1)) = |B_2g|/(2g(2g-2) 5^(g-1))`; `c_2 = 1/1200`, `c_3 = 1/25200`,
  `c_4 = 1/180000`, `c_5 = 1/660000` recomputed ( `B_4 = -1/30`, `B_6 = 1/42`, `B_8 = -1/30`, `B_10 = 5/66`). It is the
  constant of [CGYZ26, Conj. 1.5] (DIGEST_L L8; G5R F2.4) in the normalisation `varpi_1 = u + O(u^2)`.
* Scaling argument of `rmk:mai-constant` (G5_2 s.2.5 with G5R F2.3, F2.5) rechecked: the local endomorphism `\mathfrak k`
  sends `c_j -> k_0 c_j`, `i_j -> k_1 i_j + (k_1 - k_0)(l c_j + j c_(j-1))`; `Cexp' = Cexp o \mathfrak k` on `c_j` and `i_j`;
  `\mathfrak k(W) = k_0 k_1 W`, `\mathfrak k(C_1) = (k_1/k_0) C_1`; rows scale by `diag(k_0^-1, k_1^-1, k_1, k_0, 1)` (an
  `eta^-1`-isometry); every term scales by `k_0^(2g-2)`; result `c_g k_1^(2g-2) t_c'^(2-2g) + k_0^(2g-2) Q[[u]]`. Correct.
* Dictionary with CGYZ26 (Defs. 1.1-1.2, Ex. 1.3, Rem. 1.4, Prop.-Def. 4.5 (1), Conj. 1.5; Def. 1.8 not imported) and the
  differences (1)-(5): as in G5_2 s.2.6 and G5R s.2.4; the `\draftnote` on the HTML reading is kept (O3).

## 6. The conditionality paragraphs

`sec:int-inputs` and the closing paragraph after `thm:mai-main` state exactly: (Q) [Eynard:2007kz], [DuninBarkowski:2012bw] as
recalled in [Dunin-Barkowski:2015caa] (at `f = 4 W_v`), [ABDKS25], [PPZ15] with pins; (QRR); (M) Givental alg-geom/9603021v2
Thm 11.6; (G0); (T) Teleman Theorem (2) in (1.15) with GJR18 Thm 4.1; [AC96, (1.7)]; "every other step is proved in the note".
This is GAPS_SYNTHESIS_2 s.3 with DECISIONS Q2 (no IFT), Q4 (Hock out of (Q)), Q7 (Thm 11.6), plus [AC96] as OUTLINE s.0.2
requires. No other input is used by a proof (s.1 scan; S01-6). Topological recursion only on `x = W_v(p)`, `y = p` (s06
applies DOSS to `4 W_v`); no finite computation enters a proof. The only change is at pin level (PPZ Prop. 2.8, S01-4).

## 7. What remains open (none is a step of the chain)

* **O-a.** `rmk:mai-analogy`: the new `\draftnote` (the genus-3 computation behind "the graph sums alone do not single it
  out" is not recorded in appendix C). Either appendix C records it (DIGEST_A A19: kinematic plane, invisible orders 4-6,
  isolated rays, genus 3), or the sentence is dropped.
* **O-b.** `rmk:mai-tau`: unproved (existing `\draftnote`, remark only).
* **O-c.** Bibliographic (existing `\draftnote`s): CGYZ26 numbers read in the v2 HTML (`rmk:int-constant`, `rmk:mai-constant`);
  re-check arXiv:2609.13921 and 2509.19298; Teleman v3 and Givental v2 numbering; title and abstract; D11.
* **O-d.** The PPZ pin list (Prop. 2.8 added) to be ratified against OUTLINE s.0.2 by the editor.

## 8. Changes that affect other sections, and issues found there (not edited; other agents' files)

Changes made here that other sections see:
* `lem:mai-Psi` moved inside s10 before `def:mai-Fcf` (label unchanged, number 10.5 -> 10.2). Used outside s10 only by
  `def:set-outer` (s02, "`lem:mai-Psi` (i)"), which remains valid.
* `def:mai-Fcf`: the edge is now the symmetric exact-division form, the same as the current `def:set-outer` (s02) and
  `thm:flh-main` Step 1 (s07); `cor:gau-rowsNN` (ii) and `def:set-confexp` (s02/s03), which refer to `def:mai-Fcf`, are
  unaffected.
* s01 now states: PPZ enters `thm:mod-tr` and `def:set-outer`; the comparison theorem uses PPZ Def. 0.3, Rem. 2.1 and (6),
  Defs. 2.2, 2.5, 2.13, Prop. 2.12; Props. 2.4, 2.8, 2.9 enter only through Prop. 2.14 (Corollary B).
* All other labels referenced from other sections (`thm:mai-main`, `cor:mai-gap`, `def:mai-Fcf`, `cor:mai-natfact`,
  `cor:mai-actual`, `lem:mai-Ghat`, `lem:mai-coord`, `cor:mai-univ`, `rmk:mai-apriori`, `rmk:mai-tau`, `rmk:mai-stored`,
  `rmk:mai-graphlocal`, `rmk:mai-annih`, `rmk:mai-constant`) keep their content; checked at each place of use.

Issues found in other sections (for their owners):
* **X1** `s09_gamma_3.tex`, paragraph after `thm:gam-D`: "(QRR) enters the proof of `thm:gam-D` only through
  `thm:con-glue`". It also enters through `cor:gau-colours` (`D_i = D_0(zeta^-i z)`), whose proof uses `cor:gau-rowsNN` (iii),
  i.e. `lem:gau-colrows` Step 3 ((QRR) for `i != 0`).
* **X2** `b_signs_2.tex` (sec:sgn-norm): "`m_5` whose value is recorded in sec:evi-univ": it is in the opening paragraph of
  sec:evi (`c_evidence.tex`), not in sec:evi-univ.
* **X3** `c_evidence_3.tex` (sec:evi-univ): "while the other 372 terms and the total are unchanged": the sources (MECHANISM
  s.2; DIGEST_C1 UNI-7) say that the *polar part* of the total does not move; s01 and s10 now say "the singular part".
* **X4** `s08_connection_8.tex` (`rmk:con-why`, last sentence) and `c_evidence_2.tex` (the natural-factorisation paragraph):
  "only the symplecticity of `G` detects `D = 1`": symplecticity sees only `D_0(z)D_0(-z)` (the even part of `log D_0`); the
  perturbation tested breaks that product. Suggested wording: "is the only conclusion that fails for the perturbation tested,
  which breaks `D_0(z)D_0(-z) = 1`".
* **X5** `c_evidence*.tex`: no record of the genus-3 classification cited by `rmk:mai-analogy` (O-a).
* **X6** `s07_comparison.tex`/`_2.tex` (`lem:flh-hat` proof, `thm:flh-main` Step 1): `L` is used for the number of leaf
  insertions, clashing with `L = u^(-1/5)` (NOTATION table 1).
* **X7** OUTLINE s.0.2, row C10, gives the inputs of the comparison as "[PPZ15] Defs. 2.2, 2.5, Props. 2.4, 2.9"; the s07 text
  uses Def. 0.3, Rem. 2.1 and (6), Defs. 2.2, 2.5, 2.13, Prop. 2.12 (s01 now follows the text).
* **X8** `def:gau-can` (s03) is stated for the quintic system only; `cor:mai-univ` now defines the canonical gauge of a general
  pair itself, so nothing is needed in s03.
* **X9** During this review the owners of s02 and s07 changed `def:set-outer` and `lem:flh-inputs` (iii) (the expanded rows
  are symplectic iff `D_0(z)D_0(-z) = 1`; the large-radius symplecticity is not transferred). s10 is consistent with the new
  versions: no statement of s10 uses the symplecticity of the expanded rows without `D = 1`.
