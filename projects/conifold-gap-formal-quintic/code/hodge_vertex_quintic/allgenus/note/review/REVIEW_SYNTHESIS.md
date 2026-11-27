# Final review synthesis (key review/critic), part 1 of 2

Author-facing; not part of the note. Object: the assembled note (`main.tex`, 61 section files, `refs.bib`) after the eight
reviews `review/rf1..rf8` and the assembly (`review/assembler/ASSEMBLY.md`, `DRAFT_NOTES.md` + `_2`). Part 1: verdict, method,
the chain of Theorem A, the dependency graph, the conditionality, the binding rules. Part 2 (`REVIEW_SYNTHESIS_2.md`): the four
statements across abstract, s01 and body; every remaining issue with its location; the points a human reader must check first;
length and split; recommendations. No `.tex` file was edited.

## 0. Verdict

* **Theorem A: sound as assembled, with exactly the stated conditionality.** The chain C1-C13 written in the proof of
  `thm:mai-main` (s10_main.tex l.25-29) is complete: every cited statement exists, is proved earlier in the body, in appendix A,
  or is one of the quoted inputs (M), (G0), (T), (QRR), (Q); the hypotheses match at every place of use (s.2). The dependency
  graph of the labels cited in proofs has no genuine cycle (s.3). No step of the chain rests on a `\draftnote` (the only one
  inside a statement of the chain, in `def:set-graphsum`, concerns a sentence that the chain does not use), no proof cites
  appendix C, and no proof rests on a finite computation. One nuance concerns the reading of the TR input in C9: its form is
  that of DNOPS (push-forward), its sign that of DOSS (insertion form) and PPZ; bridging the two is the same comparison as in
  `def:set-graphsum` (part 2, I-1, and s.8 item 1).
* **Corollary B: sound, but one step is written only inside a `\draftnote`** (the identification of the contraction
  `Cont_Gamma[E_R, q_R; Omega]` with the term of `Gamma` in `Int (R.Omega)_(g,0)`, `def:set-graphsum`). The sentences "Corollary B uses
  ... and statements proved in the note" (s01_intro_3.tex l.21, s10_main_2.tex l.40) are true only if that argument is counted as
  proved. The argument in the note is correct (part 2, item I-1); it must become a lemma before submission.
* **Theorem C and Theorem Gamma: sound.** Theorem C's list of inputs omits (T), which Theorem A lists for the same R-matrix
  (part 2, I-2). Theorem Gamma uses no quoted input (its DLMF 24.4.1 is derived in the note from the generating function).
* **Conditionality: s01 and s10 agree word for word and match GAPS_SYNTHESIS_2 s.3 with DECISIONS Q2, Q4, Q7**, up to two pin-level
  changes that I checked and recommend ratifying (s.4.3). The abstract is looser (no [AC96, (1.7)], "classification" for
  "reconstruction"; part 2, I-3). The body uses no further theorem; a few pinned definitions/conventions used on the chain are
  missing from the lists of s01 (part 2, I-4).
* **Binding rules: met in the printed text**, with two small slips (a `.tex` comment naming design files; process wording in one
  `\draftnote`) and one borderline priority sentence (s.5).
* **Presentation: too long.** 159 pages against the design's 84 (OUTLINE_9); remarks fill about 10 pages against about 2 (DEC D8).
  A split into two papers is recommended (part 2, s.9).

Overall: ready for a careful human reading of the mathematics. Before submission: I-1 (lemma), I-2, I-3, the author decisions of
`DRAFT_NOTES.md` s.2, and the human checks of part 2 s.8.

## 1. What was read and done

* Read in full: `main.tex`, `macros.tex`, the s01 and s10 parts; `review/rf1..rf8/REVIEW*.md` (verdicts, issues, open items,
  cross-section changes); `ASSEMBLY.md`; `DRAFT_NOTES.md` (+ `_2`); OUTLINE s.0-s.0.2; GAPS_SYNTHESIS_2 s.3-s.4; DECISIONS;
  OUTLINE_9 (page budget). Read in the body: the statement of every label that the proof of `thm:mai-main` and its dependency
  paragraph cite, and the proofs at the interfaces listed at the end of s.2; the other proofs are covered by rf1-rf8 and by the
  scripts below (existence, order, cycles, citations, wording).
* Build in the scratch folder `allgenus_review_review_critic/` (pdflatex, bibtex, pdflatex x3): 159 pages; no error, no undefined
  or multiply defined label or citation, no LaTeX/package warning; 5 overfull boxes (8.6, 1.5, 1.9, 2.2, 3.9 pt); BibTeX "missing
  pages in ABDKS25". This agrees with ASSEMBLY.md s.5.
* Scripts (stdlib Python, scratch folder `tools/`, for this review only): `deps.py` parses the 61 files in `\input` order (397
  labels, 252 statement environments), attributes every proof to its statement (including `\begin{proof}[Proof of \cref{X}]`),
  and records the labels cited in statements and in proofs; `report.py`, `paths.py`, `path1.py`: unknown labels, forward
  references from proofs, strongly connected components (Tarjan), closures of `thm:mai-main` and `cor:mai-gap`, paths to
  remarks, draftnotes and appendix C; `cites.py`: every `\cite` with its pin and enclosing statement; `scan.py`, `proofscan.py`:
  rule-5 and rule-3 words; `sizes.py`: page estimates; `conifold.py`, `conifold2.py`: independent recomputation of normalisations
  (s.5, last item).

## 2. The chain of Theorem A as written in s10

The proof of `thm:mai-main` cites C12 (`thm:gam-D`), C10 (`thm:flh-main`), C9 (`thm:mod-MLK`), C8 (`cor:red-can`), C7
(`prop:tra-can`), `lem:flh-hat`, and names the rest in its dependency paragraph. Checked row by row:

| step | statements used | where proved | hypotheses at the place of use | status |
|---|---|---|---|---|
| C1 | `def:set-rows`, `thm:set-qrr` (eq:set-qrr), `def:set-graphsum`, `thm:set-GT` | s02.2-2.3; (QRR), (T) quoted | `T_Gamma` built from the rows of R; only `Cont_Gamma` of `def:set-graphsum` is used, not its draftnoted sentence | ok |
| C1' | `prop:set-cyclicLR` (i)-(iii), `cor:set-cyclicPF5` | s02.1-2.2 from (M), (G0), `lem:set-pf5` | (iii) and the end of cor. (b) are proved after `def:set-rows` (a joint proof); part order (i),(ii) -> cor.(a),(b) -> (iii) | ok |
| C2 | `thm:gau-cyc`, `cor:gau-cycsys` | s03.1 | every frame | ok |
| C3 | `prop:gau-quintic`, `prop:gau-equalprod`, `thm:gau-U`, `prop:gau-graded` | s03.2-3.3 | one-pole shape, (EP), char 0; `v = u tau^(-5/4)`, `(alpha, beta)` in `Q^2` | ok |
| C4 | `lem:gau-R1`, `lem:gau-diag`, `thm:gau-G`, `def:gau-can` | s03.4 | (N0) holds in the frames (`varpi_0(0) = t_c'(0) = 1`) | ok |
| C5 | `prop:gau-symp` | s03.5 | `r = r~ = -1` | ok |
| C6 | `thm:gau-R0` (iii)-(iv), `cor:gau-colours`, `cor:gau-rowsNN`, `lem:gau-colrows` | s03.6 | `D_0` in `1 + zQ(zeta)[[z]]` only; colour 0 without (QRR); colours `i != 0` use (QRR) (Step 3 of `lem:gau-colrows`) | ok |
| C7 | `prop:tra-can` (via `cor:tra-edges`, `lem:tra-T`, `lem:tra-F`, `lem:tra-twoslot`) | s04 | dressing at `v(u)`, `v(0) = 0`, `v'(0) = 1`, (Sp), (A), `T_Pcan = t_c`, `d(0) = 1`: the hypotheses of `thm:red-main` | ok |
| C8 | `thm:red-main` (a), (b); `cor:red-can` | s05 | ML(h), `h <= g-1` (edge graphs, all vertex genera incl. 0, 1); ML(g) (vertex) | ok |
| C9 | `thm:mod-MLK` (via `thm:mod-tr`, `prop:mod-Rtr`, `cor:mod-trF`, `prop:mod-LD`, `thm:mod-lattice`, `thm:mod-Lreg`, `thm:mod-gap` <- `prop:mod-period`, `thm:mod-V` <- `thm:mod-var` (i), `cor:mod-Ksplit`) | s06 | ML(h) for every `h >= 0`; inputs (Q) only; TR on `x = W_v` (DOSS at `f = 4W_v`, scaling `lem:mod-scale`) | ok |
| C10 | `thm:flh-main` (with `lem:flh-inputs` (i)-(iv), `lem:flh-comm`, `lem:flh-dil`, `lem:flh-frame`, `lem:flh-hat`) | s07 | `g >= 2`, (F_g), M exact in u, frame fixed; Step 1 degree count `3g-3-E-n_l` re-derived; window `z^(3g-2)` | ok |
| C11 | `thm:con-Pb`, `thm:con-glue`, `thm:con-Pa` (iv) (inside `prop:gam-hatflat`) | s08 | glue: (QRR) at colour 0, now also proves `D_0` in `1 + zQ[[z]]` (rf3 O1 closed) | ok |
| C12 | `thm:gam-lmi`, `thm:gam-main` (ii), `prop:gam-hatflat`, `thm:gam-D` | s09, appendix A | `D_0 = Phi^SP(z,1)/h(z)` exactly; LMI step 2 uses the explicit family `p(L)`, no implicit function theorem; `lem:fsp-ift` enters through `lem:fsp-fubini` (LMI Step 1) and `lem:fsp-generic` (Theorem P(a)); its use for the uniqueness sentence of LMI Step 2 is not needed | ok |
| C13 | `thm:mai-main` | s10 | fixes `g` and the frame first; `F_Gamma` rational, so `T_Gamma = F_Gamma` lies in `Q[[u]]` | ok |

Interfaces re-checked by hand: the identity `Q D = Pcan Q^m` (`thm:gau-R0` (iv), `cor:gau-colours`) and the defect formula
`R'(-z)^T eta^-1 R'(z) = M(z) (D(z)D(-z))^-1 M(z)^-1` (`lem:flh-inputs` (iii); ASSEMBLY.md s.2 writes `M(-z)^T` for `M(z)^-1`, the
same matrix since `M(-z)^T M(z) = 1`); Steps 1-4 of `thm:flh-main`; Steps 1-10 of `thm:red-main`; the glue (e) (evaluation at
`L = 1` of `A^can = A^act D_0` in `Q(zeta)[L, 1/L][[z]]`); `prop:gam-hatflat` -> `thm:gam-D`; `prop:mod-period` Steps 4-5 (the
Borel identity `Sigma_h = -4^h B_2h/(2h)` by Lagrange-Buermann, and `F^TR_h = B_2h v^(2-2h)/(2h(2h-2))`); `cor:mod-trF`
(`F^m_h = (-5)^(1-h) F^TR_h`), giving `c_h`.

## 3. The dependency graph

* No label is cited that does not exist (also: the build has no undefined reference).
* **Cycles, proof edges only:** four strongly connected components, none a circular argument.
  - `prop:set-cyclicLR` <-> `cor:set-cyclicPF5`: the joint proof "of `prop:set-cyclicLR` (iii) and of the end of
    `cor:set-cyclicPF5` (b)" (s02_setup_4.tex l.27) after `def:set-rows`; at the level of parts the order is acyclic (row C1').
  - `lem:con-red` <-> `lem:con-id`: `lem:con-red` only points ahead ("Uniqueness of the `pi_m` is not needed: `lem:con-id` below
    determines them").
  - `lem:fsp-generic` <-> `thm:con-Pa`: the appendix lemma only points to its application ("The lemma is applied in
    `thm:con-Pa` (i)"); appendix A stays self-contained.
  - `thm:gam-main` <-> `thm:gam-ginf`, `cor:gam-poly` (with `prop:gam-hatflat`): part (iv) of `thm:gam-main` is proved later by
    `thm:gam-ginf` and `cor:gam-poly`, which use only (i), (iii), `thm:gam-lmi`, `lem:fsp-dlmf5118`, `prop:gam-hatflat`, and
    `prop:gam-hatflat` uses only (ii) and `thm:con-Pa` (iv). No part uses (iv); the chain uses (ii) only.
* **Cycles with statement edges added:** larger components appear, all created by forward pointers inside statements (e.g.
  `def:set-N` -> `thm:mai-main`; `def:set-outer` -> `lem:mai-Psi`, `thm:gam-D`; `def:gau-Fg` -> `thm:flh-main`, `thm:gam-D`;
  `def:red-flatsum` -> `thm:flh-main`; `def:set-confexp` -> `rmk:mai-constant`; `def:tra-dressing` -> `rmk:mai-annih`;
  `cor:gam-equiv` <-> `cor:con-A2`). I read each one in a definition of the chain: none is used to define the object. In
  particular `def:set-outer` cites `lem:mai-Psi` (i) only for the remark that `R'(0)` is an isometry mapping `1_TFT` to 1, which
  `lem:flh-inputs` (iv) proves independently, and cites `thm:gam-D` for the symplecticity of `R'`, which no proof uses before it.
* **Forward references from proofs into the body:** 16 occurrences in 7 proofs, all pointers or deferred parts:
  `prop:set-cyclicLR` and `cor:set-cyclicPF5` -> `def:set-rows` ("(iii) is proved after ..."); `thm:gau-U` -> `lem:gau-lead` (an
  artefact: the proof of (b), (d) stands after the lemma); `thm:gau-R0` -> `cor:gau-rowsNN` ("this is `cor:gau-rowsNN` (i) for
  colour 0"); `lem:con-red` -> `lem:con-id`; `thm:con-glue` -> `sec:con-finite` ("is not used"); `thm:gam-main` -> `thm:gam-ginf`,
  `cor:gam-poly` (part (iv)). Forward references from proofs into appendix A: 53, allowed.
* **Closure of `thm:mai-main` (proof edges):** 194 nodes. Statements without a proof environment: `thm:set-mirror`,
  `thm:set-G0`, `thm:set-GT`, i.e. the quoted (M), (G0), (T) (in the whole note only these and the four s01 restatements). Remarks
  reached: `rmk:fsp-classical` only, as a pointer from `thm:gam-ginf` (part (iv), off the chain). Draftnote hosts reached:
  `def:set-graphsum` only (its draftnoted sentence is not used, s.2 row C1). Appendix C: not reached by proof edges; with the
  statement edges of definitions, only the pointer "(an example is computed in `sec:evi-model`; not used)" in `def:mod-pair`.
  No proof in the note cites a remark carrying a `\draftnote` (the only proof-to-remark edge is `thm:gam-ginf` ->
  `rmk:fsp-classical`, without a draftnote). `cor:mai-gap` adds `def:set-graphsum` directly (I-1).

## 4. The conditionality

### 4.1 s01, s10 and the abstract

* s01 `sec:int-inputs` (s01_intro_3.tex l.7-23) and s10 (s10_main.tex l.33-35) state the same inputs, and the conditionality
  sentence is verbatim the one of OUTLINE s.0.2, including [AC96, (1.7)]: (Q) [Eynard:2007kz], [DuninBarkowski:2012bw] as recalled
  in [Dunin-Barkowski:2015caa] (at `f = 4W_v`), [ABDKS25], [PPZ15] with pins; (QRR); (M) Givental Thm 11.6; (G0); (T) Teleman
  Theorem (2) in (1.15) with [GJR18, Thm 4.1]; [AC96, (1.7)]. A word diff of the two sentences shows only "above" against "listed
  in `sec:int-inputs`", and s10 adds one sentence (TR only on the model curve; no finite computation in the proof). The "Where
  the inputs enter" paragraphs agree (s10 omits
  `cor:mai-gap`, which is not Theorem A, and writes `sec:flh` where s01 writes `thm:flh-main`).
* The abstract (main.tex l.18) names (QRR) as the condition and the others as "cited theorems": the genus-0 mirror theorem, the
  genus-0 axioms, "Teleman's classification", the PPZ action, results on TR (model curve only). It omits [AC96, (1.7)] and frames
  the inputs differently from s01/s10. Not contradictory; to be aligned when the draft abstract is settled (part 2, I-3).

### 4.2 Against GAPS_SYNTHESIS_2 s.3 and DECISIONS

Same quoted works. Hock25 has left (Q) (DEC Q4) and is cited only in `rmk:mod-hock`, `conv:mod-kernel` (a convention) and
`tab:sgn-tr`. (IFT) is dropped (DEC Q2): LMI step 2 needs no implicit function theorem (rf7 s.1, confirmed by reading
s09_gamma_2.tex l.48-60). (M) is Thm 11.6; Cor. 11.7 only as the misprinted form (DEC Q7). (QRR) stays quoted (DEC Q1);
`rmk:set-qrrderiv` is used by no proof.

### 4.3 The two pin-level deviations (for ratification)

* [Eynard:2007kz, Thm 4.7 (4-30)] and the variation (5-15) of Thm 5.1 moved from (Q) to "cited but not used": correct.
  `cor:mod-trF` uses only the definition `(2-2h)F^TR_h := Sum_b Res Phi_ydx omega_(h,1)` (the form of EO Def. 4.3 (4-10)) and
  mentions (4-30) in a parenthesis; `thm:mod-var` (ii) "is used in no proof" (s06_model_8.tex l.42); `thm:mod-V` uses
  `thm:mod-var` (i) = (5-14) only (s06_model_9.tex l.45). Appendix B still lists both, as conventions.
* [PPZ15, Prop. 2.8] added to (Q): correct; it enters only through PPZ's proof of Prop. 2.14 (Corollary B, s10_main_2.tex l.31;
  and Step 5 of `rmk:set-qrrderiv`, which no proof uses).

### 4.4 The body uses no other input

All citations inside statements and proofs of the closure of `thm:mai-main` and `cor:mai-gap` (212 nodes with definitional
edges) were listed by key and pin. Every one is a pinned input, a definitional or conventional citation, or a comparison marked
as such ("This is the computation of [EO, App. A]"; "We use neither ... nor [DNOPS, Prop. 5.2]"; "[MR2276766, Theorem 2] ... is not
used"). No further theorem enters. Pinned items used on the chain but absent from s01's lists are given in part 2, I-4.

## 5. Binding rules (OUTLINE s.0)

* Rule 5: no hit in the `.tex` files for stage, workflow, agent, referee, verdict, programme, prover, gap note, cluster, route,
  AI, status tags, dates, item tags, code names; no withdrawn preprint. "verified" and "numerically": none. "checked": five times
  in appendix C (allowed) and once in a draftnote (s06_model_7.tex l.99, "were not re-checked": process wording, I-8). File
  names: only in comments; macros.tex l.1 names the design files "NOTATION.md, NOTATION_2.md" (I-7); the part comments of the
  section files name the note's own parts (structural; may be stripped before submission).
* Rule 3: no proof cites appendix C; the 26 pointers to `sec:evi*` in s01-s10 stand in the introduction, remarks, explanatory
  paragraphs, `sec:con-finite`, or a "not used" parenthesis. A scan of all proof environments for computation wording found only
  all-order symbolic steps ("computed by `lem:fsp-cov`", "computed in (iii)").
* Rule 8: TR only on `x = W_v(p)`, `y = p` and its scalings in s06; DOSS applied at `f = 4W_v`; no mirror curve; [FCCZ26] not used.
* Rule 4: the Teleman v3 and Givental v2 draftnotes stand at the first citations (s01_intro.tex l.12, s01_intro_3.tex l.11).
* Rule 11: s01_intro_3.tex l.34, "For `thm:int-gamma` we are not aware of a proof in the literature", is close to a priority
  claim (I-9). GJR18 Rem. 7.2 is described as "analogous ... stated without proof"; CGYZ26 Thm 2 only for the quintic.
* Rule 10 / DEC D8: 38 remarks, about 10.3 pages (I-10).
* Rule 2: the seven theorem-type environments without a proof are the four s01 restatements and the quoted (M), (G0), (T);
  every other one has a written proof, and the scan of the proofs found no step delegated to a computation (the line-by-line
  check of each proof is the referees' rf1-rf8; I re-read the proofs of the chain's interfaces listed in s.2).
* `\draftnote`s: 17, as DRAFT_NOTES.md says (2 main.tex, 3 s01_intro, 1 each s01_intro_2, s01_intro_3, s04_transport_3,
  s06_model_7, s10_main_2, s10_main_3, s10_main_4, 2 s02_setup_4, 3 s09_gamma_4); its classification is right.
* Independent spot checks (exact rational arithmetic; evidence for this review only): Frobenius basis of `PF4` at `u = 0`,
  `f_0 = 1 + (2/625)u^3`, `f_1 = u - (19/60)u^3`, `f_2 = u^2 + (37/30)u^3`; the log solution forces `varpi_1 = f_1 + (7/10)f_2 =
  u + (7/10)u^2 + (41/75)u^3` (conditions at `u^-3, u^-2, u^-1`: `e = 0`, `2c - 7/5 = 0`), as in `lem:set-frobbasis`;
  `tau = 1 - (4/9)u - (50/567)u^2 + ...`; `delta_Pcan = tau^(1/4) varpi_1 - u tau^(-5/4) = u^2/30 + O(u^3)` (`rmk:red-leaf`);
  `c_2..c_5 = 1/1200, 1/25200, 1/180000, 1/660000`; `B_6(1/5)/6 = 4537/3937500` (`cor:gam-free`); `Sigma_1 = -1/3`.

Continued in `REVIEW_SYNTHESIS_2.md`.
