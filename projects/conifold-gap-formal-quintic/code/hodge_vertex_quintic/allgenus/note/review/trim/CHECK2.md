# Trim check, second pass (key trim2/check)

Compared: `note/` after the three tighteners (`pass2_P1.md`, `pass2_P2.md`, `pass2_P3.md`) against
`note/archive_full_162pp/`. I copied the scripts of `check_scripts/` to the scratch folder and ran them there
(`closure.py`, `sdiff.py`, `cond.py`; Python stdlib, `py -I -B`). Method as in `CHECK.md` s.1.

## Verdict

**The chain of Theorem A and Corollary B is intact. The note is at 119 pages and builds cleanly.**

* Every statement of the old proof closure is present with the same label set, environment and proof, except the three
  items in 1a. Those three were cited in the old closure only by "not used" or "also by" remarks.
* No theorem, lemma, proposition or corollary statement differs from the archive, once one restoration is made (1b).
* The conditionality and the input lists are unchanged.
* No bookkeeping appears in any `.tex` file.
* Seven minimal restorations were made, listed in 1b-1d. Together they cost about one page (118 -> 119).

## 1. Proof-dependency closure

| | statements | closure (proof refs) | closure (proof+statement refs) |
|---|---|---|---|
| archive | 253 | 185 | 198 |
| after pass 2 | 217 | 180 | 194 |

**a. Missing from the old closure:**
* `prop:con-deg`: it was cited only by "uses neither (QRR) nor prop:con-deg" before `thm:con-Pb`, and by the "not used"
  sentence in `thm:con-glue`.
* `lem:con-weight`: it was cited only by "(also by ...)" in the proof of `prop:con-res`. That proof derives the ring from
  $V$ directly.
* `rmk:fsp-classical`: a "None of this is used" remark, cited in a parenthesis in the proof of `thm:gam-ginf`. The fact
  that the proof uses, DLMF 24.4.1, is still in the proof of `lem:fsp-dlmf5118`.

The other old-only nodes are not lost:
* `def:set-fq` and `def:mod-ancvec` are present but no longer reached through proof references.
* `cor:mai-univ` is present; it was reached only through the statement of `rmk:mai-annih`.

No new node is in the closure, no proof was lost, and no label set or environment changed.

**b. Statement text.** One lemma statement had changed.
* **Restored:** `lem:set-canon` (`s02_setup_6.tex`). The clause "which the note uses only through the equation $\Dq
  u_i=\zeta^iL$ and which are determined by it up to additive constants" is back verbatim.

After this, a check over all 253 archive statements finds no thm/lem/prop/cor difference.

I read the word diffs of all definitions and conventions that changed (24) and of the three remarks. They are rewording,
inline displays, and removed explanation, orientation examples or pointers. The removed derived facts are each available
elsewhere:
* $\eta^{-1}=\eta/25$ is in `conv:tra-forms`;
* the symplecticity argument of `def:set-outer` is `lem:flh-inputs` (iii), whose proof is complete;
* the u-form of $\PFfive$ is stated in `s08_connection.tex`;
* $|\Aut\Gamma|^{-1}$ is in the formula for $\mathrm{Cont}_\Gamma$.

Restored minimally:
* **`def:set-genus0`:** "are consequences of (M) and (G0)", for the eigenvalues, the canonical coordinates and the
  $\Delta_i$ (conditionality).
* **`conv:mod-kernel`:** "our $\omega_{0,3}$ is the negative of [EO, Thm. 4.1 (4-20)]" (a stated sign fact).
* **`def:mod-pair`:** "moving a ray across $\varsigma_b$ changes the integral by $2\pi\ri$ times such a residue". This is
  the step behind independence of the rays.
* **`rmk:set-cg`:** "$c_g>0$", with its reason $(-1)^{g-1}B_{2g}=|B_{2g}|$.

**c. Proofs.** I compared by hand the proofs whose citations dropped in pass 2: `cor:mai-gap`, `cor:mod-split`,
`thm:gau-G`, `prop:con-res`, `thm:gam-ginf`, `thm:gam-lmi` and `lem:gau-vanq`. Three needed restoration:
* **`lem:gau-vanq`** (`s03_gauge_5.tex`): restored to the archive proof. The compressed proof had lost:
  * the expansion $b_a=b^\circ_au^{-\delta_{a1}}(1+O(u))$;
  * the base case $n=0$;
  * $c^{(n+1)}_a\neq0$, which is needed for "exactly";
  * the finiteness in each $u$-order;
  * $(Y^{\mathrm C}_A)_0(0)=0$, which is needed for (d).

  Two pass-2 changes are kept: the explicit `lem:gau-vanmod`, and the moved sentence "In a conifold frame the hypothesis of
  (a) holds (prop:gau-equalprod)". The proof of `thm:gau-G` refers to the first.
* **`thm:gau-G`** (`s03_gauge_6.tex`): restored to the archive proof. The compressed proof had dropped the computation
  $z\de_u(P\tilde Y)=A(P\tilde Y)$, the constancy of $P_{11}(0,z)$ in $u$, and the converse of (b).
* **`cor:mai-gap`** (`s10_main_2.tex`), three sentences restored:
  * that the classes of $\Cexp(\TFT)$ are those of $\TFT$ with $\Cexp(\Delta_i^{\pm1/2})$;
  * the reason why the PPZ15 proofs hold over $\bbQ(\zeta)\lser s$ (identities between finite sums, polynomial in the
    coefficients);
  * Step 3 in full: edge graphs in $\bbQ\fps u$, the vertex graph, finiteness.

  The quotation of the symplecticity from `lem:flh-inputs` (iii) dates from pass 1 and was checked in `CHECK.md`.

The remaining dropped citations are only removed alternatives or pointers:
* `cor:mod-split`: "is also `lem:tra-frames` (a)";
* `prop:con-res`: "also by";
* `thm:gam-ginf`: `rmk:fsp-classical`;
* `thm:gam-lmi`: the uniqueness paragraph was already cut in pass 1, and the archive marked it "not needed below".

**d. No proof environment:** `thm:set-G0`, `thm:set-mirror` and `thm:set-GT`. These are the quoted inputs, which have no
proof in the archive either.

## 2. Conditionality (`cond.py`)

* **Abstract:** identical. `main.tex` differs from the archive only in the `\input` lines and one comment.
* **s10:** the paragraphs "Where the inputs enter" and "Conditionality" are identical. In the sentence after `cor:mai-gap`,
  the list of the PPZ15 items, "proved statements only" and "its conditionality is that of Theorem A" are unchanged. Only
  the parenthetical list of proved statements was dropped (P1); it is not an input list.
* **s01:** the "Conditionality." paragraph and the input list are identical. The "Where the inputs enter" paragraph still
  carries the three clauses restored in `CHECK.md` (sec:gau; transport lemmas with PPZ15 Defs. 2.2, 2.13 and the dilaton
  lemma, and Theorem R adds no input; Theorem P(a) and `thm:int-gamma` (a)). The `thm:int-D` clause is kept.

## 3. Bookkeeping and draftnotes

* I grepped `main.tex`, `preamble.tex`, `macros.tex` and `sections/*.tex` for dates, review, design, gaps, `.md`, todo,
  cluster, archive, trim, decision, Stanage, OUTLINE, LABELS and session. The only hits are mathematical ("Genus
  bookkeeping", Step 4, $D_w$).
* **Draftnotes:** 12 (2 in `main.tex` and 10 in sections), the same as after pass 1. Pass 2 removed none and added none,
  and their texts are unchanged.

## 4. Build, references, pages

The build ran in the scratch folder: pdflatex, bibtex, then pdflatex three times. **119 pages.** It has no error, no
undefined reference or citation, no multiply defined label and no rerun request. BibTeX gives one warning, "missing pages
in ABDKS25", as before.

**Overfull hboxes:** none. The three that remained in `a2_cyclicPF4.tex` (50.7, 12.9 and 12.0 pt) are fixed:
* the Jet sum is now a three-line `aligned` display;
* the sentence in (b) of `lem:gau-sg` was reworded;
* the closing relation in `cor:gau-cycsys` (i) is now a display.

**Dangling references:** none in `.tex`. Removed labels: `lem:mai-coord`, `sec:sgn-qrr`, `lem:con-u0`, `prop:con-deg`,
`eq:con-PL`, `lem:con-weight`, `rmk:con-710`, `cor:con-frame`, `rmk:gam-necessity`, `rmk:fsp-classical`,
`cor:gau-Gstar`, `cor:mod-basis`, `cor:mod-primary` and `sec:gau-symp`. Some `.md` planning files still name them (see
the pass reports); those files were not edited. This check removed no label.

All section files are under 250 lines; the largest is 156. The PDF was copied to `note/main.pdf`. `DRAFT_NOTES.md` and
`DRAFT_NOTES_2.md` are updated with the page count, line numbers, the cuts and typography.

| part | pages | count |
|---|---|---|
| title, abstract, contents | 1 | 1 |
| 1 Introduction | 2-6 | 5 |
| 2 Formal quintic, model, outer graph sum | 7-20 | 14 |
| 3 Cyclic systems, Theorem U, canonical gauge | 21-34 | 14 |
| 4 Vanishing tower and transport | 35-43 | 9 |
| 5 Reduction to the model | 44-45 | 2 |
| 6 The model: topological recursion | 46-60 | 15 |
| 7 Flat and hatted graph sums | 61-66 | 6 |
| 8 The connection problem | 67-76 | 10 |
| 9 Gamma identity and normalisation | 77-83 | 7 |
| 10 The main theorem | 84-87 | 4 |
| 11 Complements | 88-90 | 3 |
| A Large-radius cyclic structure | 91-93 | 3 |
| B Cyclic form of PF4 | 94-96 | 3 |
| C u-adic convergence | 97-99 | 3 |
| D Formal stationary phase | 100-110 | 11 |
| E Signs and normalisations | 111-114 | 4 |
| F Evidence | 115-117 | 3 |
| References | 118-119 | 2 |

Counts are by starting page, so each section's count includes its shared boundary page. The target of at most 120 pages
is met. The restorations of s.1 cost about one page against the 118-page tightened build, and the about-117 ideal is not
reached.
