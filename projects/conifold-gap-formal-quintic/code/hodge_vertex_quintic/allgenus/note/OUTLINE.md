# Outline of the note (design; no LaTeX here)

Working file for the section writers; not part of the note. Revised by the design editor after the gap-closing phase; every change is
recorded in `design/DESIGN_CHANGES.md` (+ `_2`). Pieces: this file (s.0 rules and the section-to-file plan, s.0.1 dependency graph,
s.0.2 logical chain, inputs, conditionality); `OUTLINE_2.md` (s.0.3 sources and gap notes, s.0.4 conflicts settled, s.0.5 decisions);
`OUTLINE_3.md` (s01; s02 items 1-9); `OUTLINE_10.md` (s02 items 10-17); `OUTLINE_4.md` (s03); `OUTLINE_5.md` (s04, s05); `OUTLINE_6.md`
(s06); `OUTLINE_7.md` (s07, s08); `OUTLINE_8.md` (s09, s10); `OUTLINE_9.md` (appendices A, B, C; page budget; writer check-list).
Notation `NOTATION.md` (+ `_2`); labels `LABELS.md` (+ `_2`); bibliography `BIBPLAN.md` (+ `_2`); errata (never in the note)
`design/SOURCE_ERRATA.md`; decisions `design/DECISIONS.md`; digests `design/DIGEST_*.md`.

Item references. Digests: `[DA A14]` = DIGEST_A statement A14; `[DB S13]`; `[DC1 B3]`; `[DC2 thm:mlk-main]` (digest label, not a note label);
`[DD1 thm:thp-b]`; `[DD2 GAM-4]`; `[DL L24]`. Gap notes (OUTLINE_2 s.0.3): `G1_3 s.2.8` = `NOTE/gaps/G1/GAPS_G1_3.md`, its section 2.8
(the sections of a gap note run across its pieces; every pointer names the piece); `G1R_2 s.4.2` = `NOTE/gaps/G1/referee/REFEREE_2.md`
s.4.2; `GS s.2.1`, `GS_2 s.3`, `GS_3 s.5` = `NOTE/gaps/GAPS_SYNTHESIS.md`, `_2`, `_3`; `REV R1`, `REV_2 M5` = `design/OUTLINE_REVIEW.md`, `_2`;
`DEC D8`, `DEC Q2` = `design/DECISIONS.md`. Marks: [W] written out; [W-c]/[W*]/[T] written tersely (expand from the same place); [G] written
out in a gap note (transcribe it with the referee corrections named next to it); [S], [M], [V] as in the digests. No step of the chain
of the main theorem carries [S] or [M] any more (s.0.2); [V] never enters a proof.

## 0. Binding rules

1. **Transcribe mathematics, invent none.** Every statement and proof step comes from the sources of OUTLINE_2 s.0.3 in their final
   corrected form: ST2b s.4-s.5, ST4b s.4-s.5, ST5b s.5-s.6, RS1 s.2, RS2 s.5, S01 s.3 and s.7, NFS s.2 and s.4, the referee corrections
   named at each statement, and the gap notes `G1`-`G5` with the corrections of their referee reports and the adjudications of GS s.2.
   For every step that a gap note closes, the gap note is the text of record; its pointer names piece and section. A step missing
   from all these sources is written as the surrounding text plus `\draftnote{...}` saying exactly what is missing; it is never filled.
   A suggestion that is not in the sources ([DW] in the digests, "suggestion" here) appears only in a `\draftnote` or in `DRAFT_NOTES.md`.
2. **Every all-orders proof step is written out in sentences.** The note is self-contained apart from the cited literature. Proved in
   the note: (CYC); the large-radius cyclic structure of the rows (`prop:set-cyclicLR`); Theorem U over any field of characteristic 0;
   (R1); the vanishing solutions; (G); the symplecticity of `P_can` and of `M_0`; the homogeneity of the model CohFT; the dilaton lemma;
   the colour lemmas and R0; the transport lemmas (frames, convergence, Lemmas E, T, F); Theorem R; the model results (the R-matrix of the
   curve, functions of `x`, the leg dictionary, the residue identity, the lattice theorem, L-regularity, the Gap Theorem, Theorem V, the
   Model Lemma); the comparison theorem; the polynomial structure; Theorem P; the glue; Theorem LMI; Theorem Gamma; Corollary B; the
   natural factorisation; the formal stationary-phase lemmas. A step marked [W-c]/[W*]/[T] is expanded from the same source place.
3. **Finite computations appear only in appendix C** (evidence), with ranges, fields and primes. No proof cites them; no proof says
   "checked" or "verified". This includes every computation marked VERIFIED in the gap notes and their referee reports. A number known
   only by computation (`t_c - vhat = u^2/30 + ..`, the order-`u` propagator, the concomitant entry `CC(f_0, y_3)`, `m_5`, the values
   `c_n(2)`) appears only in appendix C or in a remark that says it is computed and is not used.
4. **Citations with pinned statement numbers**, in the arXiv version named in `BIBPLAN.md`. Teleman: v3 PDF numbering, with
   `\draftnote{numbers from arXiv v3; journal numbering not compared}` at the first citation; Givental 1996 [MR1408320]: arXiv v2 numbering,
   with `\draftnote{numbers from arXiv v2; IMRN numbering not compared}` at the first citation. Never a number read from memory
   (Hoermander), never the withdrawn arXiv:2610.38473, never unpublished lecture notes without the authors' decision, never Givental's
   Cor. 11.7 as the pin of (M), never "DNOPS Thm 8.1" as the topological-recursion input (s.0.2).
5. **No project bookkeeping:** no stages, workflows, agents, referees, verdicts, AI, status tags (PROVED, VERIFIED, SKETCH, ...),
   correction or item tags (K.., O.., A.., C.., F.., R.., M.., D.., Q.., RF-C.., PS-C.., C3b-.., R-..; R1-R8, M1-M21, D1-D15, Q1-Q7,
   O1-O4), names of clusters or notes (G1..G5, "gap note", "cluster", "prover", "route 1/2": write "second proof"), the word "programme"
   (REV_2 M9), file, script, log, folder or code names (`fq_rows`, `ref_rows`, `ur_gauge`, `Cf` as a code name, ...), dates of writing.
   The same holds inside `\draftnote{...}`: a draftnote states what is missing in mathematical words and names no project file. The other
   normalisation used in the computations appears only in appendix C, named by its defining condition. Errata go only to
   `design/SOURCE_ERRATA.md`; the note has no errata appendix. Hypothesis names (Sp), (A), (F_g), (G), (unit), (N0), (EP) may be printed.
6. **One notation** (`NOTATION.md`), plain prose, macros from `NOTE/main.tex` (`NOTE/macros.tex`); new macros only at the top of a section
   file, with `\providecommand` and the section prefix. A renamed symbol is never used under its source name. Dummy letters inside one
   proof may be local (declare them).
7. **Hypotheses with their exact scope:** coefficient rings (`Q[v]`, `Q[[u]]`, `Q((z))[[u]]`, `Q[L, 1/L][[z]]`, `Q(zeta)[L, 1/L][[z]]`,
   `A_L`, ...); frames: the conifold frames of `def:set-frames`, `varpi_0 = f_0 + alpha varpi_1 + beta f_2`, `(alpha, beta)` in `Q^2`
   (DEC D14), wherever `varpi_0(0) = t_c'(0) = 1` or the normalisation (N0) of `lem:gau-diag` is used; any power-series `y_0` not
   proportional to `varpi_1` for (CYC); `y_0(0) != 0` for the unit shape of Theorem U; Theorem U over any field of characteristic 0;
   `D_0` in `Q[[z]]`, `D_i` in `Q(zeta)[[z]]`; genus ranges (`g_v <= g - 1` versus `g_v = g`); per outer graph, not only the total;
   `M = M_0(tau(u) zhat)` exact in `u`; formal versus analytic ("as formal series", "as rational functions of `L`"); the window (F_g);
   `d(0) = 1` stated wherever the dilaton lemma is applied (DEC Q5); legs of `L`-type with poles at `p = 0` paired by residues.
8. **Conditionality stated exactly** (s.0.2): the main theorem is conditional only on the inputs (Q), (QRR), (M), (G0), (T) and the
   formula [AC96, (1.7)], each with a pinned statement number. Say so in s01 (`sec:int-inputs`), in s10 (closing paragraph of
   `thm:mai-main`) and at each place where an input enters. Topological recursion is used ONLY on the model curve `x = W_v(p) =
   p^4/4 + v/p`, `y = p`: the quoted theorem of [DuninBarkowski:2012bw] (as recalled in [Dunin-Barkowski:2015caa] s.2.6) is applied to
   `f = 4 W_v`, "the same curve up to scaling", and transferred to `x = W_v` by `lem:mod-scale`; [ABDKS25] Thm 7.1 is applied to
   `x = W_v` itself. Results the sources state for families of curves appear at most in a remark that no proof uses. TR is NEVER used
   on a mirror curve of the quintic or of the formal quintic. No analytic theorem enters the chain beyond those cited: Stirling's
   series by the formal Proof A (appendix A), DLMF 5.11.8 by its formal version, the implicit function theorem of LMI step 2 by the
   formal `lem:fsp-ift` (DEC Q2).
9. **Read-only material:** the draft with Yannik (`refGW/BS_v2_2026-11-02/`), the root bibliography `refGW/refs.bib`, the conifold
   note folder, every project source, the digests, the gap notes and the reviews. The note has its own `NOTE/refs.bib` (built from
   `BIBPLAN.md`, with the ROOT key for every paper in the root file).
10. **Remarks** either carry their proof, written from the sources, or say "not used; not proved here" with a `\draftnote`. By DEC D8
    the note keeps only remarks that a proof uses or that explain the mechanism or the necessity of a hypothesis, at most about two
    pages in all; the remarks dropped are listed in `LABELS.md` (+ `_2`) as "dropped (D8)" and are not written. SKETCH material of the
    sources does not enter as an argument: TRM s.6; XY s.4-s.5; UR s.8 beyond its exact statements 1-2; BV s.1's converse paragraph;
    Prop. I beyond the model-outwards direction; Theorem R(c); the hatted-engine form of universality; TRM's sign-of-`kappa` argument;
    `rmk:mlk-dilaton`'s Euler-dilaton action.
11. **Literature positioning** only as in DL s.2.1-s.2.3 with REV R4 ([GJR18, Remark 7.2] states "the analogous identity for their
    Landau-Ginzburg integral" without proof; never "an equivalent identity") and REV R5 ([CGYZ26, Theorem 2] concerns the quintic; the
    other targets of CGYZ26 s.5 are not cited). The note does not claim the gap of the true quintic, priority, novelty of the formal
    lemmas, of `Ginf` ([GJR18, Lemma 7.7, (37)]) or of the chart mechanism ([MR1901075, s.10]).
12. **Files, titles, labels.** One file per section as below; a long section becomes `sNN_name.tex`, `sNN_name_2.tex`, ..., each part
    ending with `\input{sections/sNN_name_K}` of the next, and `main.tex` inputs only the first part. Labels only from `LABELS.md`
    (+ `_2`); a writer adds auxiliary labels only with its own prefix and refers to other sections only through registered labels.

| file (`NOTE/sections/`) | `\section` title (provisional) | label |
|---|---|---|
| `s01_intro.tex` | Introduction | `sec:int` |
| `s02_setup.tex` | The formal quintic, the model and the outer graph sum | `sec:set` |
| `s03_gauge.tex` | Cyclic systems, Theorem U and the canonical gauge | `sec:gau` |
| `s04_transport.tex` | The vanishing tower and transport | `sec:tra` |
| `s05_reduction.tex` | Reduction to the model | `sec:red` |
| `s06_model.tex` | The model: topological recursion on the model curve | `sec:mod` |
| `s07_comparison.tex` | Flat and hatted graph sums | `sec:flh` |
| `s08_connection.tex` | The connection problem | `sec:con` |
| `s09_gamma.tex` | The Gamma identity and the normalisation | `sec:gam` |
| `s10_main.tex` | The main theorem | `sec:mai` |
| `a_formalsp.tex` (appendix A) | Formal stationary phase | `sec:fsp` |
| `b_signs.tex` (appendix B) | Signs and normalisations | `sec:sgn` |
| `c_evidence.tex` (appendix C) | Evidence: finite computations | `sec:evi` |

### 0.1 Global dependency graph of the sections ("uses")

    s02 set-up: literature only: (M), (G0), (T), (QRR), [PPZ15], [AC96]; proves prop:set-cyclicLR, cor:set-cyclicPF5,
        lem:set-M0sympl, lem:set-modhom, lem:set-dilaton; forward references (no circularity): lem:gau-colrows (covariance
        clauses of def:set-rows), lem:gau-confexp (the conifold expansion is a ring map), cor:gau-rowsNN (rows in the ring)
    s03 cyclic systems, Theorem U, (G), symplecticity, R0, colours: s02 (lem:set-pf5, prop:set-cyclicLR, cor:set-cyclicPF5,
        def:set-rows, thm:set-qrr for the colours i != 0, def:set-frames, lem:set-frobbasis, def:set-confexp, def:set-tau,
        lem:set-canon, def:set-model, prop:set-modrows, def:set-M0, lem:set-hrot)
    s04 tower, slice, frames, transport: s02 (def:set-model, def:set-modcohft, lem:set-modhom, lem:set-dilaton, def:set-graphsum),
        s03 (lem:gau-vanmod, lem:gau-vanshape, lem:gau-vanq, lem:gau-uniq, lem:gau-diag, thm:gau-G, prop:gau-symp, def:gau-can)
    s05 Theorem R: s04 (all), s02 (lem:set-dilaton); the Model Lemma enters only as a hypothesis
    s06 the model: s02 (def:set-model, prop:set-modrows, def:set-M0, lem:set-M0sympl, def:set-modcohft, lem:set-modhom,
        lem:set-dilaton), s03 (lem:gau-vanmod, lem:gau-vanshape), s04 (def:tra-tower, def:tra-slice, lem:tra-frames, lem:tra-coeff);
        nothing about the quintic
    s07 flat versus hatted: s02 (def:set-hat, def:set-outer, prop:set-modrows, def:set-modcohft, lem:set-dilaton, lem:set-M0sympl),
        s03 (prop:gau-symp, thm:gau-R0, cor:gau-colours, lem:gau-diag), s05 (def:red-flatsum); (F_g) only as a hypothesis
    s08 connection problem: s02, s03 (thm:gau-cyc, cor:gau-cycsys, thm:gau-U, thm:gau-G, lem:gau-diag, lem:gau-HPF, lem:gau-poly,
        lem:gau-wkb, thm:gau-R0), appendix A
    s09 Gamma identity: appendix A, s08 (def:gam-c, prop:gam-rational, lem:con-form, def:con-PhiSP, thm:con-Pa, thm:con-Pb,
        thm:con-glue), s03 (cor:gau-colours)
    s10 main theorem: s05 (thm:red-main, cor:red-can), s06 (thm:mod-MLK), s07 (thm:flh-main), s09 (thm:gam-D); Corollary B also
        s02-s03 (def:set-confexp, lem:gau-confexp, cor:gau-rowsNN), (T), [PPZ15]; universality s03-s06
    appendix A: self-contained; used by s08, s09
    appendix B: conventions; referenced by all; proves nothing
    appendix C: evidence; points to all; proves nothing; nothing points to it from a proof

Keep visible in the text: topological recursion enters only in s06; (QRR) enters in s02 (`thm:set-qrr`, `def:set-rows`), s03
(`lem:gau-colrows`, Step 3) and s08 (`thm:con-glue`, colour 0); colour 0 of `thm:gau-R0` uses no (QRR); the quoted (Q) enter in s02
([PPZ15] in `def:set-graphsum`, `lem:set-dilaton`), s04 ([PPZ15] in `def:tra-edge`), s06 (EO, DOSS through DNOPS, ABDKS), s07 ([PPZ15])
and s10 ([PPZ15] and (T) in `cor:mai-gap`); (M) and (G0) enter only in s02; (T) in s02 and s10. Nothing in s02-s05, s07-s09 uses TR.

### 0.2 The logical chain of the main theorem

**Main theorem** (`thm:mai-main`, restated as `thm:int-main`; Theorem A). Let `g >= 2`, let `(alpha, beta)` in `Q^2` give a conifold
frame (`def:set-frames`), `varpi_1 = u + O(u^2)` the vanishing period and `t_c = varpi_1/varpi_0`. Let `M := M_0(tau(u) zhat)`, exactly
in `u`, and `R' := Q M^-1`, where `Q` is the matrix of the lowered rows of `def:set-rows` expanded in the frame (`def:set-confexp`). For a
leg-free stable graph `Gamma` of genus `g` let `T_Gamma := Cont_Gamma[E_(R'), q_(R'); M.TFT]` be its outer term (`TFT` the TFT of
`def:set-frob`; `def:set-outer`). Then (i) `T_Gamma` lies in `Q[[u]]` for every `Gamma` with at least one edge; (ii) for the vertex graph
`T_Gamma = c_g t_c^(2-2g) + Q[[u]]`, `c_g = (-1)^(g-1) B_2g 5^(1-g)/(2g(2g-2))`. The hatted terms are `G_Gamma = u^(2g-2) T_Gamma`
(`lem:flh-hat`). [DA A14; MS5 s.2; NF s.1 item 7; ST2 s.2 step 8; REV R6]

| step | statement | where (section, label) | inputs | proof of record |
|---|---|---|---|---|
| C1 | (M), (G0); the rows of the formal quintic, their WKB form, (QRR) for every colour; the potential is the graph sum | s02 `thm:set-mirror`, `thm:set-G0`, `def:set-rows`, `thm:set-qrr`, `def:set-graphsum`, `thm:set-GT` | (M), (G0), (QRR), (T), [PPZ15] | G1 s.1.1-s.1.4 with GS s.2.1; G1_4 s.3.1 [G] |
| C1' | flat-basis product `H * H^a = kappa_a H^(a+1)`; eigenvalues, `Delta_i`, `D u_i`; `C_0 v_0` solves PF5 | s02 `prop:set-cyclicLR`, `cor:set-cyclicPF5` | (M), (G0) | G1_2 s.2.0-s.2.7, G1_3 s.2.8-s.2.10 (a)-(b) [G] |
| C2 | (CYC): `PF4 = M_cyc(varpi_0, varpi_1)` in every frame; the cyclic system closes on PF5 | s03 `thm:gau-cyc`, `cor:gau-cycsys` | none | DB S3-S10 [W] (CYC; CF, CF2; CA) |
| C3 | the quintic system (from C1') and the pulled-back model at `v(u) = u tau^(-5/4)` satisfy the hypotheses; Theorem U; graded gauges | s03 `prop:gau-quintic`, `prop:gau-equalprod`, `thm:gau-U`, `lem:gau-comm`, `lem:gau-lift`, `lem:gau-scal`, `lem:gau-lead`, `lem:gau-degree`, `prop:gau-graded` | C1' | GS s.2.9; G2_2 s.3-s.4, G2_3 s.5, G2_4 s.6 [G] |
| C4 | (R1); the diagonal; vanishing solutions; (G); the canonical gauge | s03 `lem:gau-R1`, `lem:gau-diag`, `def:gau-van`, `lem:gau-vanmod`, `lem:gau-vanshape`, `lem:gau-vanq`, `lem:gau-uniq`, `thm:gau-G`, `def:gau-can` | none | DB S21-S28 [W]; G2_4 s.6.3, s.7, s.8 [G] |
| C5 | symplecticity `P_can(-z)^T eta P_can(z) = eta` | s03 `prop:gau-symp` | none | DB S29 [W]; G2_3 s.5.4 consequence 2 [G] |
| C6 | R0 in colour 0 (`D_0` constant in `Q[[z]]`); rows in the differential ring; colour covariance; all colours, `D_i(z) = D_0(zeta^-i z)` | s03 `lem:gau-HPF`, `lem:gau-poly`, `lem:gau-wkb`, `lem:gau-colmod`, `lem:gau-colrows`, `thm:gau-R0`, `cor:gau-rowsNN`, `cor:gau-colours` | C1', (QRR) at `q = 0` (colours `i != 0`) | DD1, PS s.4 [W]; G1_3 s.2.10 (c); G2_5 s.9.1-s.9.4 with G2R_2 s.2.9-s.2.10, s.3 C2-C4 [G] |
| C7 | transport: frames, convergence, edges of `P_can` in `K_0 (x) K_0`, leaf point, finiteness | s04 `conv:tra-forms`, `lem:tra-sum`, `lem:tra-theta`, `lem:tra-frames`, `lem:tra-conv`, `lem:tra-twoslot`, `lem:tra-E`, `cor:tra-edges`, `lem:tra-T`, `lem:tra-F`, `prop:tra-can` | [PPZ15] Defs. 2.2, 2.13 | DC1 [W]; G3 s.2, G3_2 s.1-s.3, G3_3 s.1-s.5, G3_4 s.1-s.2 with G3R_2 s.13 [G] |
| C8 | Theorem R for the flat sum `P_can.Omega^mod_(v(u))`, conditional on the Model Lemma for `h <= g` | s05 `def:red-flatsum`, `thm:red-main`, `cor:red-can` | none | DC1 C1 [W]; G3_3 s.4-s.5 [G] |
| C9 | the Model Lemma on `K` for every `h >= 0` | s06 `thm:mod-MLK` via `thm:mod-tr`, `prop:mod-Rtr`, `lem:mod-laplace`, `lem:mod-scale`, `lem:mod-xfun`, `prop:mod-LD`, `prop:mod-dii`, `thm:mod-lattice`, `lem:mod-abdks`, `thm:mod-Lreg`, `prop:mod-period`, `thm:mod-gap`, `lem:mod-w11`, `thm:mod-V` | (Q): DOSS Thm 4.1 via DNOPS s.2.6; EO; ABDKS Thm 7.1 | DC2 [W]; G4_2 s.1-s.6, G4_3 s.1-s.2, G4_4 s.1-s.4, G4_5 s.1-s.3 with G4R s.1.3, G4R_2 s.5 [G] |
| C10 | comparison: under (F_g), `T_Gamma = F_Gamma` for every leg-free genus-`g` graph | s07 `thm:flh-main` | [PPZ15] Defs. 2.2, 2.5, Props. 2.4, 2.9 | S15 s.2-s.4 with FR [W] |
| C11 | Theorem P: `Phi^SP = Phi^can`; the glue `h(z) D_0(z) = Phi^can(z, 1)` | s08 `thm:con-Pa`, `thm:con-Pb`, `thm:con-glue` | (QRR) (glue only) | G5 s.1.2-s.1.5 with G5R s.1 [G]; DD1 s.3.3 [W]; PS s.4.2 (e) with G1_4 s.3.1 |
| C12 | Theorem Gamma `c(z; 1) = h(z)` by LMI; `c_n(L) = L^(4n) ph^SP_n(L^-5)`; hence `D_0 = 1`, `D = 1`, (F_g) for all `g` | s09 `thm:gam-lmi`, `thm:gam-main` (ii), `prop:gam-hatflat`, `thm:gam-D` | appendix A (`lem:fsp-ift` in LMI step 2) | DD2 GAM-4, GAM-5, GAM-7 [W]; `thm:con-Pa` (iv) |
| C13 | assembly: C8 + C9 give the flat statement; C10 + C12 transfer it graph by graph | s10 `thm:mai-main` | none new | ST2 s.2 step 8; MS5 s.2 [W] |

**Exact final list of inputs** (printed once in s01, `sec:int-inputs`; keys of `BIBPLAN.md`; GS_2 s.3 with DEC Q2, Q4).
* **(Q) quoted theorems** (TR only on the model curve):
  [Eynard:2007kz] (math-ph/0702045v4): Def. 4.2 (4-9) with the kernel (4-3); Thm 4.2; Thm 4.6; Thm 4.7 (4-30); Thm 5.1 (5-14)-(5-15)
  with s.5.3 (5-27)-(5-30). Optional pins, not needed: Thm 4.4 (4-23)-(4-24), Cor. 4.1 (4-26)-(4-29), (A-25) (`lem:mod-xfun` is proved
  without them). [DuninBarkowski:2012bw] (1211.4021v1) Thm 4.1 with its graph weights (s.2.3), as recalled in [Dunin-Barkowski:2015caa]
  (1509.06954v2) s.2.6 (2-30)-(2-35), with the action of its s.2.3 (dilaton leaves with the sign of [PPZ15] Def. 2.13), the compatibility
  of its Prop. 4.3 (1) and the sentence before its Thm 8.1; at `f = 4 W_v`, i.e. (8-1) with `n = 4`, `a = (0,0,0,0,4v)` (critical points
  `p^5 = v`, critical values `5 p_c^4`, pairwise distinct). [ABDKS25] (2212.00320v3) Thm 7.1: (7.6), (7.7) (with its factor `1/hbar^2`),
  (7.9), (7.10); the kernel (1.1). [PPZ15] (1303.1043v3): Def. 0.3 (iii); Rem. 2.1 and (6); Def. 2.2 with the convention sentence after
  it (p.20); Prop. 2.4; Def. 2.5; Prop. 2.9; Prop. 2.12; Def. 2.13; Prop. 2.14.
* **(QRR)**: `R(q = 0) = diag(h_(w_0), .., h_(w_4))` in Givental's convention (the columns `R e~_i e^(u_i/z)` solve the QDE; legs
  `R^-1(psi)`), used for every colour in the form `Phi_i(z, 1) = h_(w_i)(z)` (G1_4 s.3.1): in `def:set-rows`, `lem:gau-colrows` Step 3
  and `thm:con-glue` (colour 0). Sources: [GJR18] (1812.11908v1) eq. (17) in the proof of Lemma 5.2; [GJR17] (1709.07392v1) eqs. (19)-(20)
  at `t = 0`; [Giv01] s.2.2 Thm 1 and the Remark after it; [MR2276766] Thm 1, Cor. 1. Its derivation from the CohFT is `rmk:set-qrrderiv`
  (not used; one `\draftnote`, DEC Q1).
* **(M)**: [MR1408320] (alg-geom/9603021v2) Theorem 11.6, Calabi-Yau case `n = 4`, `r = 1`, `l_1 = 5`, `lambda' = 0`, rewritten for `s_i`
  (GS s.2.1), with the set-up of s.9 (p.23-27) and Prop. 11.4; at `lambda_i = zeta^i`: `J(T(q) H, z) = I(q, z)/C_0(q)`, `T = I_1/C_0`.
  Hypotheses: genus-0 convexity of `O(5)`, polynomiality in `lambda` (G1 s.1.4). Cor. 11.7 only as the misprinted form (DEC Q7).
* **(G0)**: the genus-0 twisted classes form a CohFT with unit for `eta` ([GJR18] s.4.1; [PPZ15] Def. 0.3), polynomial and
  `S_5`-symmetric in `lambda` (proper moduli spaces; `E_(0,n,d)` a vector bundle of rank `5d + 1`), with the string, divisor and splitting
  axioms, the genus-0 TRR, and the quantum differential equation ([GJR18] s.5.2 eq. (14)) (G1_2 s.2.1).
* **(T)**: [MR2917177] Theorem (2) in paragraph (1.15) with [GJR18] s.4.1, Thm 4.1: the CohFT of the formal quintic is the graph sum of
  `def:set-graphsum` with the R-matrix of `thm:set-qrr`. In Theorem A it identifies that R-matrix (C1); in Corollary B, the total.
* **Standard formula** (read in the full text): [AC96] (alg-geom/9406008v1) formula (1.7) and `kappa_0 = 2g - 2 + n` (p.5), in
  `lem:set-dilaton` (i). Definitional citations: [LP18] s.0.2 (5); [LP20] s.0.2 (5), (10); [ZZ08] (the specialisation); [MR1666787].
* **Not used by the chain** (may be cited elsewhere): [Hock25] Cor. 1.2 (`rmk:mod-hock`, first proof of the Gap Theorem; DEC Q4); the
  holomorphic implicit function theorem (DEC Q2); `Int_(Mbar_(1,1)) psi_1 = 1/24` (the genus-1 logarithmic term is not printed);
  Mumford's relation, [MR1728879], [MR717614] (only `rmk:set-qrrderiv`); [MR2276766] Thm 2 (modern form of (M)); [GJR18] (13), [GR19]
  Props. 6.1-6.2, [LP20] (17), (19)-(20), [ZZ08] (corroboration of `prop:set-cyclicLR`); [CGYZ26] (`rmk:mai-constant`); DLMF 24.4.1
  (only `thm:gam-ginf`); DNOPS Prop. 5.2, the second half of DNOPS Thm 8.1, uniqueness of homogeneous R-matrices (dropped, G4_2 s.7).
* **Corollaries.** Corollary B (`cor:mai-gap`): Theorem A, (T), [PPZ15] Def. 2.13 and Prop. 2.14, and proved statements only.
  `cor:mai-natfact`: nothing beyond the inputs of Theorem A. `thm:con-Pa`, `prop:con-deg`: no quoted input. (GS_2 s.3.9)

**Conditionality sentence** (for `sec:int-inputs` and the closing paragraph of `thm:mai-main`; GS_2 s.3.10 with DEC Q2, Q4 and the
standard formula): "Theorem A holds modulo the quoted theorems (Q) of [Eynard:2007kz], [DuninBarkowski:2012bw] (as recalled in
[Dunin-Barkowski:2015caa]), [ABDKS25] and [PPZ15] with the statement numbers above, the value (QRR) of the R-matrix at the large-radius
point, the genus-0 mirror theorem (M) and the genus-0 axioms (G0) of the twisted theory, Teleman's reconstruction (T), and the formula
[AC96, (1.7)]; every other step is proved in the note."

**Steps of the chain left as `\draftnote`** (REV R7): **none**. No `\draftnote` listing open chain steps is printed. Outside the chain
(none enters the conditionality; GS_2 s.4): O1 `rmk:set-qrrderiv` (Mumford's relation; non-vanishing of `ch_(2k-1)(E)`; one
`\draftnote`); O2 the remarks kept under D8 that carry a `\draftnote` (`rmk:red-leaf`, `rmk:mai-tau`, `rmk:mod-hock`, `rmk:gam-Mcomplex`,
`rmk:gam-sevenvar`, `rmk:gam-prior`); O3 bibliographic and editorial points (`BIBPLAN.md`, "points to confirm"; the version `\draftnote`s of
rule 4, of `rmk:int-constant`/`rmk:mai-constant` and of `sec:int-gap`; the title; D11); O4 not claimed and not needed (injectivity of `ev` on the ring of
`def:set-confexp`; `deg ph_n = floor(4n/5)` as an equality; polynomiality of `P_can(0, z)`; regular semisimplicity of `A~`).
