# DIGEST_B_6 (key read_B): s.6 literature, s.7 open issues for the writers

Continuation of `DIGEST_B_5.md`. Same abbreviations (`DIGEST_B.md`).

## s.6 Literature used by the sources of this key

The objects of this key (Theorem U, (G), symplecticity, (CYC), the model rows) are proved inside the project from the definitions;
no step of their proofs cites a paper. The items below are what the sources name, with what they say about it.

| reference | statement number | what is used, and where | read in full text by the sources? |
|---|---|---|---|
| arXiv:2211.12259 (Hurwitz Frobenius manifold with one simple and one higher-order pole; rational reductions of KP) | UR, LG: none; MEC s.1, TH s.0 item 4: "Thm 2" (the `H^2` tower, via hypermaps) | UR s.6.5: a suggested route to the model's gap; LG s.0 item 5, NFS item 6: identification of the model with `H_(0;1,4)` and the tau-function statement, QUOTED. Not used by any statement of this key | no: UR s.6.5 "not read in detail", UR s.10 "was not read"; the model worker's reading of Thm 2 was not checked by me |
| Chester, Friedman, Ursell, Proc. Camb. Phil. Soc. 53 (1957) 599-611 | none | UR s.8: the uniform-reduction mechanism, SKETCH, analogy only ("Theorem U is its formal, ODE-level version and does not depend on this sketch") | no |
| CWI report ir.cwi.nl/pub/1473/1473D.pdf; the Laguerre example (doi 10.1142/9789814656054_0022) | none | UR "Sources (web)": context for uniform expansions | no |
| arXiv:math/0108100 (Givental, quantization) | none in GH | GH s.3 "Interpretation (QUOTED formalism, not used in the proofs)": covariance of quantized linear Hamiltonians, the edge-level reading of Lemma L. Other digests use its Thm 9.1, s.10, Thm 4.1 | GH: no |
| arXiv:hep-th/9204035 (Ceresole, D'Auria, Ferrara, Lerche, Louis, "Picard-Fuchs equations and special geometry"); Strominger; Candelas-de la Ossa (no arXiv numbers given) | none | CYC s.9: (CYC) is classical special geometry for a self-adjoint fourth-order operator in an isotropic frame; ST5b s.8 item 5: optional citation, "not needed since (CYC) is proved" | no: CYC s.9 "I have NOT read these sources and cite no statement from them"; CF3 s.12, CA2 s.8: not read |
| arXiv:1702.06096 s.0.2, eq. (5) and arXiv:1803.01409 (Lho-Pandharipande) | s.0.2; eq. (5) | definition of the formal quintic, hence of PF5 and the rows (FQ2 s.1; ST4b K1: cite both). DIGEST_A's subject | the critic of ST4b read the abstract pages only |
| arXiv:1303.1043 | Def. 2.2, Prop. 2.4, Def. 2.5, Prop. 2.9 | the R-matrix action on which the graph sums and `P.Omega^mod` rest (TH s.1.5; ST2 s.2 table: located in the PDF by the flat referee). Other digests | yes (flat referee, per ST2) |

No DLMF statement is needed here: `Int_0^u t^(-1/5)(1-t)^(-1) dt = (5/4) u^(4/5) 2F1(1,4/5;9/5;u)` is a termwise identity (RN s.3(a)), and
formal WKB uniqueness is proved in one line (PS s.4.2 (b), C10). The standard linear-algebra facts of Theorem U step 3 (s.3) carry no
citation in the sources.

## s.7 Open issues for the writers (conflicts, gaps, `\draftnote`s, placement)

**Placement.** ST2b s.6.3 §2 "The two cyclic systems, Theorem U, the canonical gauge, (G), symplecticity"; ST5b s.8 item 2 (vi) and MS5
s.5 item 2: "(CYC) as a set-up lemma". Proposed order inside §2, compatible with the dependencies of s.2:
1. PF4, PF5 (S1), the Frobenius basis and `gamma` (S2), conifold frames (D4);
2. (CYC): S3-S9, the Corollary S10 (remarks S11 optional);
3. the model, its rows and the uniqueness half of UR Prop. S (S18 (a), (b)); the scaling limit as a remark (S19 (a));
4. cyclic systems with one simple pole; S12 (hypotheses for the quintic and the pulled-back model; `v(u)`, `tau`);
5. Theorem U (S13 with S14-S16) and graded gauges (S17);
6. (R1) (S21) and the diagonal (S27);
7. vanishing solutions (S22-S25; S26 as a remark), the gauge identity (G) (S28), the canonical gauge (D12);
8. symplecticity (S29), (G*) (S30); the residue relations (S31, S32) only if the transport section needs them;
9. R0 (S20): at the end of §2 or at the start of the connection section; it must be formulated as the glue's writer uses it (PB s.3).
Remarks after the connection theorem (`D = 1`): S18 (c) (the symbol is the model), S35 (Corollary F), S36; S33 and S34 optional.

**Gaps in written proofs (each becomes a `\draftnote` unless the authors supply the step).**
* G1 (Theorem U, step 3; UR s.3 (ii)). "At is regular semisimple over K, so its commutant is K[y] and `gl_5(K) = K[y] (+) im(ad At)`":
  no justification in UR. Needed: why `At` is regular semisimple (its characteristic polynomial is `x^5 - B`, `B` of valuation -1) and the
  decomposition for semisimple elements. Not in any source read; do not invent: `\draftnote`.
* G2 (Theorem U, step 5). The projection formula `pr(X) = (1/5)[Sum_a xi_a/W_a^(m)] At^m` is stated; its derivation (trace form) is not
  written. `\draftnote`.
* G3 (Theorem U, step 8). `pr(Dl y^j) = (tr Dl/5) y^j` stated without derivation. `\draftnote`.
* G4 (Theorem U, step 10). The uniqueness clause ("the only freedom is one constant per order, which is `P -> P c(z)`") is one sentence;
  it carries S29 (symplecticity uses the uniqueness for `(At, At)`) and S28 Cor. (a). The induction (every regular gauge arises from the
  construction; the constants are realised by `P c(z)`) must be written or flagged. `\draftnote`.
* G5 (Theorem U, graded clause; S17 step 1). Existence of a graded gauge from the covariance `S A S^-1 = zeta A`, `S = diag(zeta^-a)`,
  and "the freedom is `c(z^5)`": one sentence each (UR s.3, s.4; GI s.1). `\draftnote`.
* G6 (S17 step 3). `D_i(z) = D_0(zeta^-i z)`: one sentence in UR s.4, "confirmed" in PS s.4.5; the three-line computation from the
  covariance of both column sets (CS s.1; UR s.2) and the gradedness of `P_can` is not written. `\draftnote`.
* G7 (Theorem U, field). UR works over `C[[u]]` (or convergent series), GI quotes U over `Q[[u]]`; that the construction stays over Q is not
  discussed. State the theorem over the ring the note uses and add a `\draftnote` for the authors.
* G8 (S27 step 4). "`prod p_(0,a) = 1` canonically" needs `p_(0,a)(0) = 1` for all a; one line not written.
* G9 (S28 Cor. (a)). Uniqueness of the canonical gauge among regular gauges with `P Y^C = Y^C_Q` is implicit (uses S13's uniqueness).
* G10 (S34). The order-u propagator formula has no written derivation (UR s.5 gives the result; NFS item 12's "PROVED by hand" not
  located). Not needed: omit, or values in the evidence appendix.
* G11 (S35 (b)). The uniqueness of `tau` rests on RN's (R-eq), whose derivation is only VERIFIED. Not needed: remark with `\draftnote`.
* Everything in (CYC) is written (RS2 s.6.2: no finite identity remains); its one computed-only number `B(f_0, y_3) = -9/10` is not
  needed (s.3 S5 step 7): print Bmat with that entry marked, or omit it.

**Superseded statements and labels (do not state as conjectures in the note).**
* UR s.0 item 1 "Cyclic form (VERIFIED)": now S9-S10, PROVED (RS2 s.5 item 45).
* UR s.4 "(H) Purity" and "Canonical normalisation (CONJECTURAL)", GI s.6 "(H) and `c^act = 1` to all orders: CONJECTURAL", UR s.6.5
  "Conjecture N ... CONJECTURAL for `g >= 4`", MEC s.3 open inputs: superseded by `D = 1` (MS5 s.2; RS2 s.6.2). UR Prop. S (c) and
  Corollary F become corollaries of `D = 1`.
* UR "Theorem D equality PROVED", "Proposition E PROVED", "v unique given v(0) = 0": overclaims (NFS s.4): bound only; genus 2, 3 only;
  unique given `v = u + O(u^2)`.
* GH "frame dependence" of the leaf bookkeeping: withdrawn; `-p_(1,4)(0) = 1/15` and `p_(1,4)(0)/2 = -1/30` in every frame (GR s.4 item 1;
  S01 C6). GH "needs the non-linear genus-0 descendant extension": superseded (S01 s.2).
* GI s.0 item 5 / GH s.1 "the production R' of genus 2-5 IS the canonical gauge at the orders used": an inference from the recorded
  agreement of production and referee rows, not a test (GR s.4 item 2; S01 s.7).

**Conflicts between sources and their settlements (for the record; all settled).**
* Scope of (CYC): CF-E5 / ST5b s.9 ("`y_0(0) != 0`") versus CA-K2 (not used): CA adopted (RS1 A1).
* Frame labels in the critic's checks (RS1 A2): relabel; harmless.
* Log-free rows as an input of S10: both CF and CA right about different statements; supplied inside the glue (RS1 A4).
* `Cf` and C4: no conflict (RS1 A5).
* Dictionary `H^a <-> x^a dx` (RN side remark) versus `x^(c_a) dx` (LG, UR): LG/UR right; RN's item withdrawn (NFS s.2 item 1).
* Existence of a polynomial G to all orders: RN "CONJECTURAL" (written before U) versus UR "PROVED given (H)": PROVED-mod((H), U) (NFS s.2
  item 2); now modulo U and `D = 1`.
* Frame independence of degrees and rule W (RN) versus special frames (verifier): the degree formula is a bound attained in generic
  frames (NFS s.2 item 3, item 6).
* Euler-dilaton sign: `N_1 = -1`, not an ordering convention (NFS s.2 item 8).
* Sign of `c_5`: unified (S01 C11; s.4.2).
* GI's (G*) "failure" in a first attempt: a floating-point bug in the test, not a wrong statement (GI s.8; GR s.4 item 4, reconstruction
  caveat). Irrelevant for the note.

**Statements to define with care.**
* The vanishing solutions by their initial value `-e_1` (D13; GR s.2), not by "only negative powers of z".
* R0 (S20) with the same exponential on both sides, `D` in `Q[[z]]`, and the ring of comparison (PS s.4.1-s.4.2).
* Conifold frames: the programme family `varpi_0 = f_0 + alpha g + beta f_2` for the normalisations `varpi_0(0) = t_c'(0) = 1`
  (CA-K2; RS2 s.5 items 27-28); (CYC) needs less (S9); Theorem U needs `y_0(0) != 0` (S12).
* "Conjecture N" as UR s.6.5 words it ("every stable graph with an edge is `O(u^(2g-2))`") refers to the hatted engine normalisation
  `G_Gamma = u^(2g-2) F_Gamma` (ST2 s.2 step 6): regular `F_Gamma` iff `O(u^(2g-2))`. Coordinate with DIGEST_A.
* Theorem U is an input "(U)" in MS5 s.2, but binding rule (2) requires the note to PROVE it: the main theorem is then conditional only on
  (Q) and (QRR) (rule 8), with the gaps G1-G7 as `\draftnote`s until the authors close them.

**Rule (5) items (keep out of the note; DRAFT_NOTES.md only).** Review history (AI-audited, refereed once/twice), human checks pending
(MS5 s.5 items 1.2 (CYC) and 1.6 (Theorem U, (G))), code names (`fq_conifold`, `fq_rows`, `ur_gauge`, `gi_*`, `cyc_*`), the "stored"
normalisation, primes, logs, frame lists of the codes, the label `Cf`.

**For NOTATION.md.** Resolve the clashes of s.4.3, above all `c` (hatted exponents vs the cyclic coefficients: both occur in §2 on the
same page), `D` (`q d/dq`, `D(z)`, UR's `D_0`), `tau` (CYC's automorphism), `M` and `P` (CYC), the two "Prop. S", "Theorem G", the two
"Lemma 1" and "Lemma 2", `(R1)`, `eta_n`. Fix `H = d/dv`, `e = tau_1(1) + v tau_0(H)`, `m5` (s.4.2).

**For LABELS.md.** The labels proposed in s.2 use the topic prefixes `cyc`, `gau`, `mod`, `van`, `sym`, `res`, `hat`; the section prefix
chosen by the design replaces them.

**Bookkeeping of this digest.** 38 statements (S1-S38). The end of s.3 (S20-S38) is in `DIGEST_B_5.md`, not in `DIGEST_B_4.md` as the
piece table of `DIGEST_B.md` announces. In `DIGEST_B_2.md` S18 (b), "PS s.1" refers to PS s.1 (iv)-(v) (the re-derivation of the
T-system and of the pivot `4n`), which states that this half of Prop. S does not use (H).
