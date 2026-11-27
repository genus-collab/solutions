# DIGEST_A, part 7: s.7 Open issues for the writers

Abbreviations: `DIGEST_A.md` header. Each item: the issue; sources; what to do; suggested placement.
Addendum to s.0: besides the reading listed there, ST2 s.1 table rows M4 and M8, LIT.md l.75 and l.97, LIT_2.md l.183 and GAMMA.md
l.148 were seen as one-line search hits only.

## s.7.1 Conditionality and quoted inputs (binding rules 4 and 8)

1. **Definition of `F_g` and its citation.** FQ2 s.1 cites arXiv:1702.06096 s.0.2 for the formal quintic; MS4 s.2 and the Stage-4 evid
   referee cite eq. (5) of the same paper; a Stage-4 note calls arXiv:1803.01409 "Lho-Pandharipande, formal quintic" (GAMMA.md l.148).
   The generating variable of `F_g` (`Q = e^T` with `T = I_1/I_0`) and the inclusion of `d = 0` are implicit (FQ2 s.2: `F_2(0) = N_{2,0}`).
   Do: define `N_{g,d}` and `F_g` in the set-up with one pinned citation; `\draftnote{statement number to confirm: 1702.06096 s.0.2 or
   eq. (5)}`. Placement: §2 (set-up), first definition.
2. **Genus-0 mirror theorem.** Every row is built from the I-function and PF (FQ2 s.1 "Genus-0 data (all verified exactly)"; D6). The
   identification of the formal quintic's quantum differential equation with PF (equivalently `J = I` after the mirror map) is not in
   the list (Q), and (QRR) is stated only as the `q = 0` normalisation (MS5 s.2). The (QRR) sources rest on arXiv:math/0008067 and
   arXiv:math/0110142 Thm 1, Cor. 1 (MS5 s.2); LIT_2 l.183 mentions Givental's equivariant mirror theorem (search hit only). Do: settle
   with the literature digest whether this input is covered by the (QRR) quotes; if not, name it in the hypothesis list of the main
   theorem (rule 8). Candidate root keys: MR1653024 (Givental, mirror theorem for toric complete intersections), MR2276766
   (Coates-Givental). `\draftnote` in §1 and at the definition of the rows until settled.
3. **Givental-Teleman reconstruction for the formal quintic** (A6). The identity "`F_g` = graph sum with the QRR-normalised R" is part of
   (QRR)'s package through Teleman's uniqueness, "read only through a summariser" (MS5 s.2); MS5 s.5 item 3 asks for "Teleman s.8.2
   verbatim". Do: state it as a quoted theorem with arXiv:0712.0160v3 Theorem (2) (s.1.9), Prop. 8.5(ii) (s.8.2), plus the R-action
   conventions of arXiv:1303.1043; `\draftnote{Teleman s.8.2 to be read verbatim}`. Placement: §2, after the rows.
4. **`Lambda^vee` and the Hodge class** (D8). The convention is in FRAMEB A3.2 (not read). Needed only for the background remark on
   Hodge frames and for the vertex-level derivation of (QRR). Do: take the conifold note's definition if identical; otherwise define it.
5. **Project documents behind the code conventions.** The docstrings cite "FRAMEB s.0.1-0.2, A2-A3", "FRAMEB 0.1", "FRAMEB A3.2",
   "FRAMEB C Step 1", "draft (1.1)": project notes, not literature. Do: state D6, D7, D9 as definitions in the note; never cite them.

## s.7.2 Conventions to fix before writing (NOTATION.md)

6. **Row normalisation and index placement** (C4). Choose one of `R_{ai}`, `P`, `Q`, `Qhat` as the named rows; check OU2 s.1's
   `[z^n] eta(H^a, R(z) e~_i)` against doc:fq_rows and doc:ref_rows. Placement: §2 and the sign appendix.
7. **`Psi(u)`** (D12). My reading `Psi_{aj}(u) = Qhat_{a,0}^(j)(u)` comes from doc:ref_rows's Teleman form; confirm with UR s.4
   (`G = Psi(u)^-1 Ghat Psi(0)`) and OU3 s.3 (`Psi eta^-1 Psi^T = eta^-1`) before writing.
8. **`xhat`.** OU3 s.1's `xhat = 1 + O(u^3)` is the frame `(0,0)` value; in general `xhat = varpi_0^cf(u) = 1 + alpha u + ..`.
9. **Frame scope of the theorem.** Evidence: two frames at genus 2-5 (second frame one prime at genus 3 and for genus-5 regularity),
   seven frames at low order (NF s.1, NFS item 13). The chain claims nothing frame-specific as far as the entry notes show ((CYC) in every
   frame, MS5 s.3.4; the leaf bookkeeping frame independent, S01 C6; `gamma_t + w_1 = 7/10`), and MS5 s.3.4 widens (CYC) to any
   power-series `varpi_0` not proportional to `g`, with `varpi_0(0) != 0` needed for Theorem U. Do: confirm with the Theorem U, S15 and
   Theorem R digests that no step fixes `(alpha, beta)`; then state "for every conifold frame". Placement: §1 statement.
10. **Model normalisation, `H`, `e`** (C9, C10). One superpotential normalisation for the whole note (proposal: `W = p^4/4 + v/p`,
    `H = d/dv`, `e = tau_1(1) + v tau_0(H)`), and one sentence giving the dictionary to `J_k` and `eps = -4 zeta^-i zhat`.
11. **Letter clashes** (C21): `s`, `tau`, `alpha`, `c`, `K`, `G`, `R~`, `M`, `L`, `D`, `H`, `P`, `Q`, `T`, `N`, `W`, theorem names.
    The renaming advice of S01 C9 ("(G)", "Gap Theorem", `K_vc`) is already adjudicated.
12. **Hatted variables as a definition.** Theorem S15 uses "Lemma 3's description of the hatted variables", which matches the code
    docstrings (ST2 s.2 step 6). Do: write D12 (with A7) as the note's definition and check it against the S15 digest.

## s.7.3 Gaps, sketches, unstated corollaries

13. **Gap of the total** (A15): VERIFIED `g = 2..5`; for all `g` it follows from the main theorem, A6 and the `M`-independence of the
    total, but no source read states it as a corollary. Do: state it with exactly these dependencies, or omit it;
    `\draftnote{corollary not in the sources as such}`. Placement: after the main theorem.
14. **Input (b) of the symbol identification** (A11): "the conifold-frame map respects D on the row subring, a standard stationary-phase
    lemma" has no citation in the sources read. Check UR s.2 and the formal stationary-phase digest; `\draftnote` if not covered.
15. **Uniqueness of `tau`** (A12): "last step to be written out" (NFS s.3 item 5). Not used by the main theorem. Do: remark with the
    evidence of V5-V7 only, or a `\draftnote` naming the missing step; do not prove it.
16. **Degree bound of `G_n`** (A10): through the rule-W equivalence, "index bookkeeping ... unrefereed" (NFS s.5, s.1 item 10). It is used
    only in A9's count of `zhat`-orders. Do: cite UR's Theorem D (bound only, NFS s.4) as proved in the note, or drop the count.
17. **Each hatted factor's power of `s`** (A7 step 1) is stated in OU2 s.1 with a one-line reason; write it out from D12.
18. **Elementary steps stated without derivation**: the fifth-order equation and PF4 from the I-function (A2 step 3); the indicial
    polynomial and the Frobenius recursion (A3 steps 1-2). Write them out (they are elementary) or cite.

## s.7.4 Presentation, scope and evidence wording

19. **"Hodge-vertex"**: for the formal quintic the linear Hodge frames fail (A18); the vertices of the theorem are the model CohFT
    `V_tau` (D15). Proposal: call the vertex factor "the model factor" and explain the analogy with the Hodge factor in one remark (A20).
20. **Coloured versus uncoloured**: the theorem is per uncoloured outer graph (the colour sum is inside the vertex tensors, D13, C19).
21. **Literature without pinned statements**: "HKQ" (FQ2 s.1, s.4), the "usual conifold-gap value", Hurwitz Frobenius manifolds, the
    Chester-Friedman-Ursell map. Rule 4: give a pinned statement or drop; `\draftnote` for the authors. Root refs.bib has no HKQ entry;
    for Hurwitz spaces Dubrovin:1994hc exists but no source cites it.
22. **Superseded statuses**: NF's ARGUED and CONJECTURAL items (factorisation, Theorem U, model, (H), normalisation, model gap,
    bookkeeping, Conjecture N for `g >= 6`) are superseded by MS5 and S01 s.2 "Record notes now outdated"; NFS's PROVED-mod((H), ..)
    items lose the (H) dependence through `D_0 = 1`. The note uses none of these labels (rule 5).
23. **Evidence caveats** that the appendix must state plainly (V6-V7, V13): genus-5 arithmetic probabilistic; the second frame at one
    prime (genus 3; genus-5 regularity); `tau_6` and `zhat^10..zhat^12` of `M_0` uncontrolled; the regrouping identity checked against
    genuine vertex tensors only at genus 2-3; the second implementation shares the split-sum method and the genus-5 graph list; `c_g` and
    the total's gap are `M`-independent; `F_3(0)` from the project's R-matrix graph sum not computed (MS5 s.5 item 3).
24. **No conflict found between the sources on the statement itself**: NF s.1 item 7, ST2 s.2 step 8 and MS5 s.2 state the same
    Conjecture N(g); the differences are in the inputs (s.1.5 identification dropped, MS2; (QRR) made explicit, MS4; (F_g) discharged
    and (CYC) added, MS5), and in the (QRR) source (MS4 s.2 superseded by MS5 s.2, correction O11).
25. **Upgrades recorded in the sources that the note must reflect**: the degree bound of `ph_n` (PROVED, MS4 s.3) replaces FQ2 s.3's
    "observed `= 4n` for `n <= 7`" as the input of the all-genus pole bound (A7); the single-leaf bookkeeping is derived (S01 E7); the
    model gap holds for all `g` modulo quoted theorems (S01 E8).

26. **Novelty and prior art** (MS5 s.3.6): "No paper read proves (Gamma). The write-up should cite these and claim no more novelty"
    (Givental arXiv:math/0108100 s.10, Lemma 10.4; GJR arXiv:1812.11908 Remark 7.2, Lemma 7.7, eq. (37); Iritani arXiv:0903.1463v3
    Lemma 4.20). For the formal stationary-phase appendix: "No source read by lit states these formal lemmas for complex data with
    statement numbers. The project's proofs are self-contained" (MS5 s.3.2). Placement: §1 and the (Gamma) section.
27. **The polynomial structure of `F_g`**: FQ2 s.4 places `F_2` in "the Lho-Pandharipande shape
    `C[L^(+-1)][A_2, A_4, A_6, C_0^(+-1), C_1^(-1), K_2]` with `C_0^(2g-2)` homogeneity" without a statement number or paper. MS5 s.5 item 2
    lists "the polynomial structure" first in the write-up. Do: pin the LP statement (arXiv:1702.06096 or arXiv:1803.01409) with the
    digest that owns the polynomial structure; `\draftnote` until then.

## s.7.5 Suggested placement (write-up plan ST2b s.6.3 as updated by MS5 s.5 item 2)

* §1 Statement: A14 with the hypotheses (Q), (QRR) and the input questions of items 2-3; A16; A17; the motivation A18-A20; scope A24;
  A15 if kept.
* §2 Set-up (before the cyclic systems and Theorem U): D1-D11, A1-A6 (A5, A6 as quoted statements), D9 with the conventions C2-C7,
  D12-D13, A7-A9; D16-D17 (the model), A11; the objects D19-D20 as forward references.
* Natural-form material (with the gauge sections): D14-D15, A10, A12 (remark), A21-A22 (optional), A25.
* §8 Main theorem: the chain of A14 (proofs from the other digests), with (QRR) entering only in the glue.
* Appendix, evidence: V1-V13, A18 data, A26, A27.
* Appendix, sign conventions: C1-C20; NOTATION.md: C21.
