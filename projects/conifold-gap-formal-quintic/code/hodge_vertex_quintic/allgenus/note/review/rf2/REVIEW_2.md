# Review rf2: `sections/s02_setup*.tex` (part 2 of 2)

Continuation of `REVIEW.md` (verdict s.0, issues S1-S15 s.1).

## 2. What remains open (all marked in the text; none is in the conditionality)

1. **S11, the new `\draftnote` in `def:set-graphsum`** (push-forward translation of [PPZ15, Def. 2.5] versus insertions in
   `Cont_Gamma`). Mathematically complete (argument inside the draftnote); open only as a written step. Possible closures, for the
   authors: (a) write it as a short lemma (two lines, induction on the number of forgotten points); (b) pin the one-point comparison
   `psi_j = pi^* psi_j + D_(j,n+1)`, which G1 s.1.3 item 4 reports as [AC96, (1.9)] (read in the gap phase) - but [AC96, (1.9)] is not
   among the inputs of OUTLINE s.0.2 (only (1.7)), so this would add an input (rule 8). Chain status: Theorem A as written in s05/s07
   compares contractions with contractions and does not use it; Corollary B uses it (`thm:set-GT`, `def:mai-Fcf`, `cor:mai-gap` Steps
   1-2); and **C9 may use it**: `thm:mod-tr` (s06) reads the graph formula of [DOSS, Thm 4.1] (as recalled in DNOPS s.2.6), whose dilaton
   leaves carry `psi^k` at the vertex, as the action of [PPZ15, Defs. 2.2, 2.13] (`def:set-graphsum`); if those leaves are insertions
   (not push-forwards), the same comparison is used there. To be settled by the reviewer of s06.
2. The existing `\draftnote`s: Givental arXiv v2 numbering (first citation of [MR1408320], `thm:set-mirror`); O1 in `rmk:set-qrrderiv`
   (Mumford's relation and the non-vanishing of `ch_(2k-1)(E)`; the genus-one values `Int psi_1 = Int lambda_1 = 1/24` of Step 3 are also
   unpinned, but the remark is not used).
3. Not changed, noted: (a) the translation (d) of the pin of (M) uses `pi_* psi_0^k = psi_0^(k-1)` for stable maps (the string
   equation), unpinned, as accepted by G1R s.1.2; it belongs to the quoted input (M). (b) `def:set-fq` asserts that the localisation
   value can be specialised at `lambda_i = zeta^i`; the note justifies only the invertibility of the `R^1`-term (from G1 s.1.4); the
   virtual normal bundle weights are also non-zero there (a convex combination of two distinct fifth roots of unity is not a third one),
   but this is part of the definition quoted from [LP18, LP20] and is not written; no change. (c) `def:set-rows`: "This means" lists
   properties (a)-(c) that do not determine `R` uniquely (right multiplication by `exp` of an odd constant diagonal series); `R` is the
   R-matrix of the CohFT by (T) and its value at `q = 0` is (QRR), as the definition says next; left as is.

## 3. Changes that affect other sections

| change in s02 | where it is used | status / action |
|---|---|---|
| `def:set-outer`: `E_(R')` is the symmetrised division-free edge; symplecticity of `R'` derived from `cor:gau-colours`, `prop:gau-symp`, `lem:set-M0sympl` and `thm:gam-D` | s07 `def:flh-outer`, `lem:flh-hat`, `lem:flh-inputs` (iii), `thm:flh-main` Steps 1-2, `rmk:flh-fail` (v); s10 `thm:mai-main`, `cor:mai-gap`; s01 `thm:int-main`; App. C (`E_R = E_M + M^-1 E_(R') M^-T`, valid when `R'` is symplectic) | s07 already reads (current text): "`T_Gamma = Cont_Gamma[E_(R'), q_(R'); M.TFT]`, with the edge of `def:set-outer`" and `lem:flh-inputs` (iii) "`R'` symplectic iff `D_0(z) D_0(-z) = 1`; modulo `z^(3g-2)` under (F_g)" - consistent (same identity, `M(-z)^T = M(z)^-1`). `lem:flh-hat`'s rule `E_(R^') = s^4 E_(R')(s^4 psi, s^4 psi')` holds for the new edge (each `N` term scales the same way). `rmk:flh-fail` (v) (current text) says the symplecticity of `R'` is equivalent to `D_0(z) D_0(-z) = 1` and holds modulo `z^(3g-2)` under (F_g) - consistent. |
| same: the constant term of `R'` justified by `lem:mai-Psi` (i) in `Nring (x) Q(zeta)` | s07 `lem:flh-inputs` (iv), s10 | s07 (iv) now has the same ring computation (current text) - consistent |
| `lem:set-dilaton` (iii): hypothesis extended to `Q(zeta)((s))`, `varkappa` in `1 + uQ[[u]]` | s07 `lem:flh-dil` (applies it over `A = Q(zeta)((s))`), s04 `lem:tra-conv` (f), `lem:tra-F` | s07 may now cite (iii) directly; its own sentence "It applies over `A`, because ..." stays correct |
| `def:set-fq` now contains `eta(H^a, H^b) = 5 delta`, `eta^(cc') = delta/5` | `def:set-frob` keeps the display; other sections cite `def:set-frob` for `eta` | no action |
| `prop:set-cyclicLR` (iii), last sentence of `cor:set-cyclicPF5` (b): proofs moved after `def:set-rows` | s03 (`lem:gau-colrows`, `cor:gau-rowsNN`, `thm:gau-R0`), s08 `thm:con-glue`, s10 `lem:mai-Psi` cite the statements | statements and labels unchanged; no action |
| `cor:set-cyclicPF5` (a): sentence on additive constants of `u_i` | `lem:set-canon`, s03 `thm:gau-R0` (i) | no action |
| new `\draftnote` S11 in `def:set-graphsum` | s06 `thm:mod-tr`; s04 `def:tra-edge`; s10 `def:mai-Fcf`, `cor:mai-gap`; s01 `sec:int-inputs` (list of draftnotes, if any) | s06 reviewer: see s.2 item 1; s01: not an input, not a chain step of Theorem A as written |
| Teleman draftnote removed from `rmk:set-qrrderiv` | s01 carries it at the first citation | no action |

Related places in other sections found during the review (not edited):
* s08 `rmk:con-genus` (current text) already says that `D_0(z) D_0(-z) = 1` "would" halve the conditions and that the transfer of the
  `q`-series symplecticity is not available - consistent with S1. `rmk:con-why` (current text: "symplecticity (which could at most give
  `D_0(z) D_0(-z) = 1` ...)") - consistent.
* s10 `def:mai-Fcf` (current text) uses the same symmetrised edge `E_X` and states why; `cor:mai-gap` Step 1 derives the symplecticity
  of `Cexp(R~)` from `D = 1` - consistent with S1.
* Labels: none added, none removed, none renamed. Every `\cref` in an s02 proof now points backwards, to a pinned citation, or (in the
  two definitions `def:set-confexp`, `def:set-outer`, not in proofs) forward to statements whose proofs do not use s02's later objects.

## 4. Checks on record (finite, exact, evidence for this review only; nothing enters the note)

Scratch scripts in the review's scratch folder (`checks/illus.py`, `checks/model.py`, Python stdlib):
* `f_0 = 1 + (2/625)u^3`, `f_1 = u - (19/60)u^3`, `f_2 = u^2 + (37/30)u^3`, `varpi_1 = u + (7/10)u^2 + (41/75)u^3`; in the frame `(0,0)`
  `t_c = u + (7/10)u^2 + (41/75)u^3 + (9/20)u^4`; `tau = 1 - (4/9)u`, `v = u + (5/9)u^2`: all as printed.
* T-system of colour 0: `Fs_0 = 1 - (3/20)t + (9/800)t^2 + (11571/80000)t^3`, `Fs_b(0) = eps_b`, closing relation to `t^6`;
  `Rhat^m(-t)^T eta^-1 Rhat^m(t) = 1` numerically through `t^5` (all colours).
* By hand: `PF4` coefficients (`e_k` of `1..4` and `1..5`), `P_(-3)`, `P_(-2)`, `P_(-1)'(0) = -7/5`, the obstruction matrix and kernel;
  `D^k ell`; `eta`; `Delta_i`; square roots; `h_1 = -3/20` versus `-11/60`; Givental's exponent `(f hbar)^-1[g_1 Sum lambda - Sum lambda' g +
  lambda_i g_1 (Sum l - (n+1))]`; the genus-one sign check `r_1 = a_1/12`; Mumford-relation rewriting `H_g = W_i^(g-1) Prod h(a)/a^g`;
  the degree count of `lem:set-modhom`; `M_0(zeta^-1 t) = Pi^-1 M_0(t) Pi`; the binomial step of `lem:set-dilaton` (iii); the reversion
  at `g = 2`; `c_2..c_5`; the automorphism `i_j -> i_j - k c_j` commutes with `D` (incl. `i_3 -> i_3 - k c_3` from the concomitant relation).
* Compilation: `pdflatex`, `bibtex`, `pdflatex` x2 in the scratch copy with the current versions of all sections: no error, no
  undefined reference, no overfull box in s02.

## 5. Reading record

Read in full: all eight s02 parts; OUTLINE (s.0-s.0.2), OUTLINE_3, OUTLINE_10, NOTATION (+`_2`), LABELS (s02 rows), BIBPLAN, DECISIONS,
GAPS_SYNTHESIS s.0-s.2, GAPS_SYNTHESIS_2 s.3-s.5, G1 parts 1-5, G1R parts 1-2; in other sections the statements cited from s02 (s03
`def:gau-conc`, `lem:gau-const`, `lem:gau-adj`, `lem:gau-conc`, `lem:gau-colrows`, `thm:gau-R0`, `cor:gau-rowsNN`, `cor:gau-colours`,
`def:gau-Fg`; s07 parts 1-2; s08 `rmk:con-problem`, `rmk:con-genus`, `thm:con-glue`, `rmk:con-why`; s09 `thm:gam-main`, `thm:gam-D`; s10
`thm:mai-main`, `lem:mai-Psi`, `def:mai-Fcf`, `cor:mai-gap`). External: PPZ15 (1303.1043v3) Def. 2.5, Remark 2.6, Prop. 2.7 (text
extracted from the arXiv PDF, p.23-24) - no statement relating the translation by push-forward to insertions was found in Section 2.
