# Trim check (key trim/check)

Compared: the trimmed note `note/` against `note/archive_full_162pp/`. The scripts are in `review/trim/check_scripts/`
(Python stdlib). Run them with `py -I -B closure.py <NOTE>`, `sdiff.py <NOTE> stmt|proof <labels>`, `cond.py <NOTE>`.

## Verdict

**The chain of Theorem A and Corollary B is intact.** Every statement in the old dependency closure is still in the
trimmed note. It has the same label set and environment, and every one that had a proof still has a proof. No theorem,
lemma, proposition or corollary statement changed. Definitions were tightened in wording only: every defining clause and
stated fact is kept. One pointer was restored (see 1c). One sentence of s01 on where the inputs enter had become wrong and
was restored, and a second one was completed (see 2). After the restorations: **127 pages**, a clean build.

## 1. Proof-dependency closure (`closure.py`)

Method: follow the `\input` chain from `main.tex` (comments stripped). Parse the environments thm/lem/prop/cor/defn/rmk/
notation/example/conv/proof. A label inside an environment stands for that environment. A proof `[Proof of \cref{X}]`
belongs to X; any other proof belongs to the statement before it. Then take the closure from `thm:mai-main` and
`cor:mai-gap` over `\cref/\Cref/\ref/\eqref` (multi-label references included). This was run in two modes: references in
proofs only, and references in proofs and statements.

| | statements | closure (proof refs) | closure (proof+statement refs) |
|---|---|---|---|
| archive | 253 | 185 | 198 |
| trimmed | 228 | 182 | 196 |

**a. Missing, proof lost, or label or environment changed: none.** All 25 removed statements lie outside the old closure.
The new closure contains no statement that is not in the old one.

**b. Statement text changed (whitespace aside).**
* thm/lem/prop/cor: **none**.
* Definitions in the closure: `def:mai-Fcf`, `def:set-confexp`, `def:set-fq`, `def:set-frob`, `def:set-graphsum`,
  `def:set-modcohft`, `def:set-model`, `def:set-outer`. I read each word diff. All changes are rewording or removed
  explanation, and every defining formula and stated fact is still there. Specific points:
  * `def:set-frob`: now cites `def:set-fq` for $\eta(H^a,H^b)=5\delta$. That derivation is still in `def:set-fq`, with
    $\eta^{-1}=\eta/25$.
  * `def:set-modcohft`: the sentence "identification with a Hurwitz Frobenius manifold is not used" is removed.
  * `def:mai-Fcf`: the parenthetical on why the symmetric edge form is needed is shortened.
* Remarks in the closure (proof+statement mode only): `rmk:mai-annih`, `rmk:mai-constant` and `rmk:set-cg` were shortened.
  All the values and the comparison are kept. `rmk:mai-constant` now points to `sec:sgn-lit`.
* Labels in running text that the closure cites (equations): no text changed. `sec:*` labels naturally differ.
  `sec:con-finite` (removed) was cited in the old closure only by a sentence of the proof of `thm:con-glue` saying that it
  is *not* used. The same holds for `prop:con-deg` and for `sec:evi-calib`.

**c. Proofs citing fewer labels than before.** I checked each case by hand against the archive. All are compressions of
the same argument, removed parenthetical pointers, or "not used" remarks:
* `cor:mai-gap`: the symplecticity argument in Step 1 is now quoted from `lem:flh-inputs`(iii) with $\Dn=1$
  (`thm:gam-D`). That lemma states the identity, precedes Section 10, and does not depend on `cor:mai-gap`. The deduction
  $\Cexp(\Qr)$ symplectic follows in one line from `lem:set-M0sympl`. Steps 3-4 were merged.
* `thm:con-Pa`: Step 1 now applies `lem:fsp-generic` (statement unchanged; its proof uses `lem:fsp-ift` and
  `lem:fsp-dq`). It checks the hypotheses: $S'(L)=0$, $S''=-1/\delta$, $A_{\rm St}(\sigma,0)^2=-\mathsf r$, and
  $b\in\bbA_\sigma\fps z$.
* `thm:con-glue`: steps (b)-(e) were merged; `prop:con-deg` and `sec:con-finite` were cited only as "not used".
* `thm:gam-lmi`: the uniqueness paragraph is removed. The old text said that it "is not needed below".
* `thm:mod-V` (`thm:mod-var`), `thm:red-main` (`conv:tra-forms`), `lem:fsp-stirling` (`lem:fsp-cov`, `lem:fsp-homog`),
  `lem:fsp-generic` (`def:gam-c`, `thm:con-Pa`): the cited items were "not used", "also follows from" or "is applied in"
  remarks.
* `thm:mod-MLK` (`cor:mod-trF`, `def:mod-ancvec`), `thm:mod-tr` (`lem:set-pushins`), `lem:gau-confexp`, `lem:mod-genus0`,
  `lem:set-grading`, `lem:set-unitfrob`, `cor:set-cyclicPF5`, `lem:set-pushins`, `prop:mod-LD`, `prop:set-cyclicLR`,
  `cor:gau-rowsNN`: only a pointer to a definition was removed, or routine steps were compressed. The word diffs show
  every step still stated.
* **Restored** (`s02_setup_8.tex`, `def:set-outer`): "and $\Psi^{\mathrm m}_{v(u)}$ the isometry of
  \cref{def:set-modcohft}". This returns that definition to the archive wording for the symbol it uses.

**d. Proofs with no proof environment:** `thm:set-G0`, `thm:set-mirror` and `thm:set-GT`. These are the quoted inputs
(G0), (M) and (T), and they have no proof in the archive either.

## 2. Conditionality (`cond.py`)

* **Abstract:** identical. `main.tex` differs from the archive only in the `\input` lines and one header comment.
* **s10 (`sec:mai-thm`):** the "Where the inputs enter" and "Conditionality" paragraphs are identical. The sentence after
  `cor:mai-gap` that lists the [PPZ15] items and says "its conditionality is that of Theorem A" is identical.
* **s01:** the "Conditionality." paragraph and the input list (Q), (QRR), (M), (G0), (T), [AC96 (1.7)] are identical.
  Two sentences were removed from it: the pointer to the cut `rmk:set-qrrderiv`, and "each input is stated again where
  it enters".
* **Fixed in s01** (`s01_intro_3.tex` l.17, "Where the inputs enter"): the trimmed sentence read "sec:gau, Theorem R,
  Theorem P(a) and part (a) of thm:int-gamma use no quoted input beyond the genus-zero data". For Theorem R this is
  false: through the transport lemmas it uses [PPZ15, Defs. 2.2, 2.13] and the dilaton lemma, which rests on [AC96]. The
  archive content is restored in one sentence:
  * sec:gau uses no quoted input beyond the genus-zero data;
  * the transport lemmas use [PPZ15, Defs. 2.2, 2.13] and the dilaton lemma, and Theorem R adds no input;
  * Theorem P(a) and part (a) of thm:int-gamma use no quoted input.
* **Completed in s01** (l.21): the trimmed sentence on `thm:int-D` no longer said that it uses the matrix identified by
  (T) and the genus-zero data resting on (M) and (G0). Both are restored in one clause.
* The s01 "Where" paragraph still omits the detailed [PPZ15]/[AC96] places. It now points to `sec:mai-thm`, whose
  paragraph is complete and unchanged, so no claim is lost.

## 3. Bookkeeping in `.tex` files

I grepped `main.tex`, `preamble.tex`, `macros.tex` and `sections/*.tex` (comments included) for dates, review, design and
gap identifiers, group keys, decision numbers, archive, todo, cluster and `.md`. Nothing was found. The only hits are
mathematical: `fontenc` T1, Step 4, and the remark title "Genus bookkeeping", which is also in the archive.

**Draftnotes:** 16 in the archive and 12 now. The six removed draftnotes were inside the removed remarks `rmk:set-qrrderiv`,
`rmk:tra-consistency`, `rmk:gam-sevenvar`, `rmk:gam-Mcomplex`, `rmk:mai-tau` and `rmk:mai-analogy`. Two were added, in
`rmk:mai-graphlocal` and `rmk:mai-necessity`, pointing to computations that are no longer tabulated. Every other
draftnote is kept verbatim.

## 4. Build and pages

Built in the scratch folder: pdflatex, bibtex, then pdflatex three times.
* **127 pages.** No error, no undefined reference or citation, no multiply defined label, no rerun request.
* BibTeX: one warning ("missing pages in ABDKS25"), as before.
* The PDF was copied to `note/main.pdf`.
* All section files are under 250 lines (largest 156).
* Labels removed by this check: none.

The note is still over the 110-120 target (127 pages). The remaining options are in `ASSEMBLY.md` s.5. None of them
concerns the chain.
