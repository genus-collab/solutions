# Fidelity check of the abridged version against the long version

Date: 2026-11-20. Scope: `main_summary.tex` + `summary/*.tex` (abridged) against
`main.tex` + `sections/*.tex` and the built `main.aux` (long, 110 pp, read-only).

## Result

- No statement of the abridged version was found wrong or stronger than the long version.
- One proof sketch contained a wrong formula (Section 3) and one pointer was misattributed (Section 6).
- Nine caveats were missing, the main ones being the steps that the long version itself flags as
  still to be checked (flop identity, restriction to subtori, three-torus, Stirling step).
- Twelve wording corrections.
- All were fixed in the abridged files. Rebuild: 33 pages, latexmk exit 0, no LaTeX warning
  (no undefined reference or citation, no multiply defined label), largest overfull box 4.57pt.

## What was compared

| Abridged file | Long file | Environments | Sketches | Read by |
|---|---|---|---|---|
| abstract, `s1_intro` | `main.tex`, `s01_intro` | 5 (Theorems A-D, corollary) | - | main checker, line by line |
| `s2_setup_a`, `s2_setup_b` | `s02_setup` | 20 | - | main checker, in full |
| `s3_rows` | `s03_rows` | 12 | 6 | delegated reader, in full |
| `s4_leaves_a`, `s4_leaves_b` | `s04_leaves`, `s04b_closedforms` | 15 | 14 | delegated reader, in full |
| `s5_identification` | `s05_identification` | 11 | 8 | delegated reader, in full |
| `s6_regularity` | `s06_regularity` | 7 | 3 | delegated reader, in full |
| `s7_rigidity` | `s07_rigidity` | 10 | 10 | main checker, in full |
| `s8_evidence` | `a_evidence` | 6 check items | - | main checker, in full |

Totals: 80 statement environments, 41 sketches, about 90 running-text claims, 153 pointers (151 numbered, 2 to Appendix A).
For Sections 3-6 the comparison was made by four separate readers, one per section, each reading
the whole long section (statements and proofs); their edits were reviewed afterwards as diffs, and
the cited theorems and the scale law of Section 6 were re-compared against the long text.

Checks made:
1. Each statement against the pointed long statement (hypotheses, conclusion, signs, constants,
   ranges, rings and fields; merged or shortened statements still true and not stronger).
2. Conditionality: inputs, "which result uses which", "what is not claimed", and every draft note
   of the long version (all 40 were read).
3. Each sketch against the long proof.
4. Each `[L, ...]` number against `main.aux`, by script and by content.
5. Front matter.

## Discrepancies and changes

Severity: W = wrong, S = stronger than proved, C = missing caveat, T = wording.

### Introduction (`s1_intro.tex`)

1. **C** First paragraph said "the complete proofs are in [L]". The long version carries draft
   notes marking several steps as written out beyond its sources and to be checked.
   Changed to "the proofs are in [L]", with a sentence pointing to the list in the final section.
2. **C** Input (ii). Long `s01_intro` l.109 notes that Teleman's paper is cited with the arXiv
   numbering, not compared with the journal. Sentence added.
3. **C** Input (iii). Long `s01_intro` l.110 (also `s04_leaves` l.47, `s06_regularity` l.8):
   the numbering of the companion paper is provisional. Sentence added.
4. **T** "What is not claimed" (ii): "it depends on alpha [L, Remark 7.16]". In the long version
   this rests on computed values at genus 2, 3, 4 (Remark 7.16, Appendix A.7), not on a proof.
   Changed to "it depends on alpha, as its computed values in genus 2, 3, 4 show".

Theorems A-D, the corollary, the four inputs, "which result uses which" and the four
"not claimed" items are otherwise word for word those of the long introduction.

### Section 2

No discrepancy. All 20 definitions and lemmas agree with Definitions 2.1-2.26 and Lemmas
2.5-2.24; the merged lemma (2.11 + 2.12) and the condensed Lemma 2.23 are true as written.

### Section 3 (`s3_rows.tex`)

5. **W** (in a sketch) Sketch of Lemma 3.19 / Proposition 3.20 said that Stirling's series
   "gives the same prefactor" as Laplace's method. In the long proofs Stirling gives
   `e^{f/z} (2 pi z (a+b)/(ab))^{1/2}`; it becomes the Laplace prefactor only after the shift
   factor `b/(a+b)`, respectively `A alpha_1/(E(E+z))`. Sketch corrected.
6. **T** Section opening listed Laplace's method and Stirling's series as the analytic inputs;
   the Beta integral (DLMF 5.12.1), used in the long proof and listed in input (i), was missing.
   Added.

### Section 4 (`s4_leaves_b.tex`)

7. **T** Remark on the Calabi-Yau condition said Theorem A holds "at every weight point with
   (ND)"; by Definition 2.1 such a point lies on the Calabi-Yau torus, so the off-torus
   statement was empty as written. Changed to "of the five-torus", as in the long introduction.
8. **C** After Corollary 4.45: "the direct proof of Theorem A covers alpha_0 = 0 as well".
   The long remark gives the reason and a restriction: the point Q = 0 is never used, and the
   description of `h_nw` as the constant at Q = 0 (Lemma 5.12) needs alpha_0 != 0. Added.

### Section 5 (`s5_identification.tex`)

9. **C** The section did not say that the argument fails on the three-torus (alpha_0 = 0),
   as the long version does after Proposition 5.1 and in Remark 5.18. Sentence added, with a
   new pointer [L, Remark 5.18].
10. **C** Lemma 5.2 was said to rest on two citations; the long proof also uses the
    fundamental-class and splitting axioms of the virtual class "in standard form", with no
    reference (draft note, long l.78). Added.
11. **T** Idea of the proof of Lemma 5.3: localisation on the factor C^2 and the divisor
    equation (cited) were omitted. Added.
12. **T** "is a multiple of lambda_k lambda_{k-1}" -> "non-zero multiple" (Lemma 5.6).
13. **T** Final sketch: "the listed denominators vanish only over ..." -> "vanish, at finite Q,
    only over ..." (Lemma 5.16(ii) treats Q = infinity separately).

### Section 6 (`s6_regularity.tex`)

14. **C** After the cited vanishing theorem: "its proof is not reproduced" under a header
    pointing to [L, Theorem 6.8] could be read as "the proof is in [L]". Changed to "neither
    here nor in [L]".
15. **T** "Part (iv) is the only place where the Calabi-Yau condition is used". The long text
    says the same, but the cited vertex theorem (6.9(ii)) also has the hypothesis
    E = eps_1 + eps_2. Reworded to name both; this departs from the long sentence on purpose.
16. **T** Sketch, step (ii): the span for (L3) lacked `i,l >= 1`, and the self-loop case was
    not distinguished. Corrected.
17. **W** (pointer) "first values of beta and instances of M_irr are in [L, Remarks 6.11, 6.14
    and 6.15]": Remark 6.14 is the specialisation to weight points. Pointer moved to the
    sentence it belongs to.
18. **T** "use Theorem A and input (iii) only" -> "Theorem A and the two cited theorems of the
    companion paper (input (iii)) only".

### Section 7 (`s7_rigidity.tex`)

19. **C** Flop identity, part (ii) (long Proposition 7.15(iii)). The long proof has two draft
    notes (l.356, l.362): the weights of the two small resolutions are stated without
    reference, and the identification of the factors at 1/Q = 0 with those of the flop at its
    own Q' = 0 "is not spelled out in the source" and "should be checked". The abridged sketch
    presented the step without reservation. A reservation was added to the sketch, saying
    that part (i) is not affected.

All ten statements (Lemmas 7.2, 7.4, 7.6-7.9, 7.11, 7.17, Proposition 7.10, Corollary 7.18)
agree with the long version; the sketches follow the long proofs.

### Section 8 (`s8_evidence.tex`)

20. **C** The section "What remains conditional" did not mention the flagged steps. New
    subsection "Steps of the long version still to be checked", six items:
    Stirling step (long `s03` l.566); the moving-weights list and the projection-formula step
    (`s05` l.74, l.453); restriction to a general subtorus (`s07` l.70); sign conventions of
    the d = 1 localisation formula (`s07` l.240); the flop identification (`s07` l.356, l.362);
    the ring map alpha_0 -> 0 and the restriction for d >= 1 on the three-torus (`s07` l.380,
    l.383). It states which of these Theorems A-D depend on. The paragraph on conditionality
    now also records the provisional numbering, the arXiv locators, that the results of the
    paper containing the conjecture are referred to by content, and the virtual-class axioms.
21. **T** Check (ii): "Other classes or the point frame give singular edge graphs". Appendix
    A.3: this holds for the class at colour 1; the class at colour 2 is not constrained.
    Corrected.
22. **T** Check (iii): "fails for every g >= 1" -> "for every g >= 1 tested".
23. **T** Added, from Appendix A.9 and A.6: the alpha_0-adic step and the regularity at
    Q = infinity were tested only through the totals and the rational reconstruction.

Counts: W 2 (one sketch formula, one pointer), S 0, C 9, T 12.

## Pointers

All 151 numbered `[L, ...]` pointers resolve in `main.aux` to a statement of the right type, and each
points to the statement it is used for (after fix 17). Two were added: [L, Remark 5.18] in
Section 5 and a second [L, Appendix A] in Section 8. Merged pointers checked by content:
Lemmas 2.11 and 2.12; Definitions 3.5 and 3.7; Propositions 3.8-3.10; Lemmas 5.5-5.7 (5.6 is
the Chern-character lemma, part (ii) of the merged statement); Lemmas 6.3 and 6.4; Lemmas 7.6
and 7.7; Proposition 7.15 (parts (i) and (iii)).

## Front matter

- Author: "The GENUS collaboration" (author line changed on 2026-11-23, both PDFs rebuilt and the
  title block rechecked). No person's name as author, no email, no `\draftnote`.
- The introduction ends with `\subsection*{Acknowledgements}` and the required sentence.
- Abstract: same claims as the long abstract, with "state, and sketch the proofs of".

## Build

`latexmk -pdf -interaction=nonstopmode -halt-on-error main_summary.tex`: exit 0, 33 pages.
No LaTeX warning in `main_summary.log`. Overfull boxes: 4.57pt (introduction, statements
paragraph, present before the check) and 1.22pt. All summary files are under 250 lines
(largest 196, the macro file; largest section file 166).

## What could not be checked

- The cited papers themselves: the locators and the content of the theorems of the companion
  paper, of Teleman's paper and of the paper containing the conjecture were compared with the
  long version only.
- The mathematics of the long version, except where a check was cheap: exact-rational tests of
  the closed form of the leaf coefficients through t^6 and of the loop flow through t^5 at two
  weight points (Section 4), a first-order check of the Beta-phase lemma (Section 3), and hand
  checks of the scale law and of the constants of Lemma 5.12.
- The numbers of the finite-order checks, beyond their agreement with Appendix A.
- Sections 3-6 were not re-read in full by the main checker (see the table above).

## Observations on the long version (not edited)

- `s06_regularity.tex` l.186: "it is the only place" where the Calabi-Yau condition is used,
  against the hypothesis of Theorem 6.9(ii) (l.196) and Remark 4.2.
- `s04_leaves.tex` l.47: "every weight point with (ND)" for a statement off the Calabi-Yau
  torus, against Definition 2.1.
- `s03_rows.tex` l.14 and l.512: the analytic inputs are listed without the Beta integral,
  which the proof of Lemma 3.19 cites.
- `s05_identification.tex` l.410: the draft note asks whether rule (P2) is stated with a
  general parameter; it is (Proposition 3.9), so that note can be closed.
- Remark 4.14: the parenthesis on the pole of `z H_n` reads as contradicting the next sentence.
