# Decisions on the design's open points D1-D15 (main session, 2026-11-15)

Author-facing; not part of the note. Answers to `OUTLINE*.md` s.0.5 (D1-D15), taking `design/OUTLINE_REVIEW.md` (R1-R8) and
`OUTLINE_REVIEW_2.md` (M1-M21) into account. Items marked **gap** are handed to the gap-closing phase (`note/gaps/`); the design is
revised only after that phase, by one editor who applies R1-R8, the recommended M-items it agrees with, these decisions and the gap
results.

| item | decision |
|---|---|
| D1 | Define `Pi_i`, `K`, `K_0` directly from `Y^C` (the recommended option). |
| D2 | **gap** (cluster G1): pin the genus-0 mirror theorem for the formal quintic with a statement read in full text, and derive or pin the large-radius flat-basis product `prop:set-cyclicLR` (R1). Name the input (M), not (S1) (R2). |
| D3 | Keep Corollary B (the gap of the total), with its derivation written out (**gap**, cluster G5). |
| D4 | Route 1 (Hock, Cor. 1.2) as the proof of record of the Gap Theorem, route 2 (ABDKS + the Borel identity) as a second proof, so that Hock is optional; with M12. |
| D5 | **gap** (cluster G4): the base `h = 1` of Theorem V written out (the exact computation by hand, or the vanishing-period route with the transfer to `x = W_v`). |
| D6 | **gap** (cluster G4): agreement of the TR R-matrix and the model R-matrix to all orders, by an argument written out (e.g. homogeneity plus the uniqueness of the homogeneous model column, UR Prop. S, or a pinned uniqueness theorem). |
| D7 | **gap** (cluster G2): every step of Theorem U stated without proof in the sources is proved or pinned. |
| D8 | Keep only remarks that a proof uses, or that explain the mechanism or the necessity of a hypothesis: the 7-variable chart, the necessity examples (failing controls), LMI for complex `M0` (one sentence), the weight lemma (one sentence), the finite route (short), the uniqueness of `tau` (remark, labelled as such if unproved). Drop the other optional remarks unless a proof needs them. At most about two pages of remarks in total. |
| D9 | Names as proposed: Theorem A (main theorem), Corollary B (gap of the total), Theorem C, Theorem Gamma; S15 becomes "the comparison theorem". |
| D10 | Keep the remark comparing `c_g` with the constant of [CGYZ26, Conj. 1.5], with a `\draftnote` to compare the definitions of the conifold expansion. |
| D11 | Author decision, left open: a `\draftnote` in s01 asking whether to cite the authors' companion notes and how to relate the note to the sentence of [Br25, s.1.1] on R-matrix methods. |
| D12 | The working title and abstract of `design/TITLE_ABSTRACT.md` stand as drafts, marked for the authors. |
| D13 | One model normalisation, `W_v(p) = p^4/4 + v/p`; the others only in appendix B (with R3's DNOPS point `a = (0,0,0,0,4v)`). |
| D14 | Frame scope of the main theorem: the programme family `(alpha, beta)` in `Q^2`. |
| D15 | Theorem U and `P_can` over Q, proved in cluster G2 with D7 (**gap**). |

Literature (DIGEST_L L8-L9): the introduction cites [CGYZ26, Theorem 2] for the gap of the quintic threefold itself and [CGYZ26,
Remark 1.9] for the formal quintic, with the wording of DIGEST_L_4 s.2 and the corrections R4-R5. No claim beyond DIGEST_L_4 s.2.1.

## Answers to the gap critic's questions (`gaps/GAPS_SYNTHESIS_3.md` s.7), main session, 2026-11-15

| item | decision |
|---|---|
| Q1 | Keep (QRR) a quoted input (GJR18 Lemma 5.2 with (17); GJR17 (19)-(20); Giv01; MR2276766). G1's Steps 1-5 appear as a remark in s02 with one `\draftnote` naming O1 (Mumford's relation and the non-vanishing of `ch_(2k-1)(E)`, to be read before promoting it). |
| Q2 | Drop (IFT): LMI step 2 uses the formal `lem:fsp-ift` of appendix A. |
| Q3 | One label for the symplecticity of `M_0`: `lem:set-M0sympl`. |
| Q4 | Route 2 (ABDKS Thm 7.1 with the Borel identity, `S_h = -4^h B_2h/(2h)`) is the proof of record of the Gap Theorem; route 1 ([Hock25] Cor. 1.2) becomes a remark. [Hock25] leaves (Q). This replaces D4. |
| Q5 | Adopt the homogeneous definition of the vertex tensor in `def:red-flatsum` (G3-2f); state `d(0) = 1` where the dilaton lemma is applied. |
| Q6 | Move `def:gam-c` and `prop:gam-rational` into `sec:con-form` (s08), so that `thm:con-Pa` (i) has no forward reference; s09 then uses them. |
| Q7 | (M) is pinned to [MR1408320] Thm 11.6 (with the rewriting for `s_i`; Cor. 11.7 only as the misprinted form). Add LLY97 as a second source only if the writer of s02 reads it and pins a statement; otherwise no second source and no `\draftnote`. |
| D11 | Pending the authors' answer: a `\draftnote` in s01 (default). |

## Ratifications after the review and fix phases (main session, 2026-11-15)

| item | decision |
|---|---|
| R-PPZ8 | Ratified: [PPZ15, (8)] (the formula in the proof of Prop. 2.12, used for d = 1, 2) is the pin of `lem:set-pushins` (i); AC96 (1.7), (1.9) do not give the class identity. PPZ15 was already in (Q). |
| R-I5 | Ratified: EO (4-30) and (5-15) are not inputs (no proof uses them); PPZ Prop. 2.8 is an input through the proof of Prop. 2.14. |
| R-I2 | Ratified: (T) identifies the R-matrix of `def:set-rows` for Theorems A and C alike; Theorem C's proof uses R only through `def:set-rows` (b), `R(0) = 1` and (QRR). |
