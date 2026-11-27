# Fidelity check of the abridged version against the long version

Date: 2026-11-20. Scope: `main_summary.tex` + `summary/*.tex` (abridged, "[A]") against
`main.tex` + `sections/*.tex`, `main.aux`, `main.toc` (long, "[L]", read-only, not modified).

## 1. What was compared

Every file of [A] was read in full next to the corresponding file of [L] (all of
`s01`-`s11` and `a_evidence`, about 4300 lines).

- 54 numbered statements of [A] (theorems, propositions, lemmas, corollaries, definitions,
  remarks), each against the statement of [L] it points to: hypotheses, conclusion, signs,
  constants, index ranges, rings and fields.
- 91 displayed formulas of [A], each against the display of [L].
- 20 proof sketches of [A], each against the proof in [L] (does the sketch describe the
  argument actually used; is it marked as a sketch or idea).
- The running text of every section: about 45 paragraphs carrying claims (dependencies,
  "uses only", "does not use", values at the conifold point, summaries of lemmas quoted
  with a pointer only).
- All pointers `[L, ...]` (114 before the fixes, 122 after), each against `main.aux`.
- Front matter, the list of inputs, "what is not claimed", the methods not used.

A few constants were recomputed by hand as a control on the transcription:
`r_1, r_2, r_3` from the recursion; `beta_2` and the genus-2 leading term; `Delta_j(0)`;
the value of the edge at colour 0 at `t_c = 0`; the constant `-pi^2/2` in the Barnes
argument; the exponent count in the scale law.

## 2. Result in one paragraph

No statement of [A] on local P^2 was found to be wrong, and none is stronger than the
statement of [L] it points to. All pointers match `main.aux`. The defects found are in
the treatment of conditionality: the resolved conifold was presented as proved, the
outstanding checks flagged in [L] were not mentioned, and the summary of the evidence
was slightly stronger than the appendix. They are fixed. Counts: wrong 0; stronger than
proved 2; missing caveat 3; wording 3.

## 3. Discrepancies and fixes

### D1. Resolved conifold corollary presented as proved (stronger than proved)
- Where: `s11_final.tex`, first paragraph and Corollary 11.1 of [A] (= [L, Corollary 11.4]).
- [A] before: "[L, Theorem 11.1] ... is stated in [L] with a description of its proof";
  "With these, the proof of ... gives the following", then the corollary with no
  hypothesis.
- [L]: Section 11 opens with a note that its text is to be decided and that the rank-2
  leaf structure "is only stated". Theorem 11.1 has no proof in [L]; the remark after it
  lists ingredients of a proof "not reproduced here", with a note that the pointer is to
  be replaced. The proof of Corollary 11.4 starts from Theorem 11.1.
- Fix: the paragraph now says that Section 11 of [L] is not in final form and that
  Theorem 11.1 is stated without proof; the corollary is titled "conditional" and begins
  "Assume the rank-2 leaf structure [L, Theorem 11.1]"; a sentence after it says that it
  is not proved as long as that theorem is not. The last bullet of the section says the
  same.

### D2. "Complete proofs are in the long version" without exception (missing caveat)
- Where: abstract and `\thanks` in `main_summary.tex`; "About this version" in
  `s01a_intro.tex`; last paragraph of `s01b_intro.tex`.
- [L]: complete for local P^2; not for Section 11 (see D1).
- Fix: abstract: the conifold counterpart "is conditional on a leaf-structure statement
  whose proof is given in neither version", and "complete proofs of the results on local
  P^2 are in the long version"; `\thanks` likewise; "About this version" names the
  exception and points to the final section; the plan says that the counterpart is
  conditional and that the final section lists the outstanding checks.

### D3. Status and inputs of [L, Proposition 11.7] and [L, Theorem 11.1] (missing caveat)
- Where: `s11_final.tex`.
- [A] before: "its proof follows Sections 4, 6 in rank 2", no inputs named.
- [L]: the proof is written "in the terms of" the rank-3 proofs; two notes say that
  step (b) relies on a Laplace calculus "in a separate notation", that steps (a), (c)
  go through the rank-3 proofs, that it is undecided whether to write them out, and that
  the phase at generic weights is to be compared with [BS]. Inputs listed in [L]:
  Tel12, GP99, CG07, PPZ15, Mumford's formula and FP00 Theorem 4, the Aspinwall-Morrison
  formula (FP00, p. 3), Stirling's series (DLMF 5.11.1, 5.12.1); for Theorem 11.1,
  Stirling's series with remainder (DLMF 5.11.1, Section 5.11(ii)).
- These inputs are not in the list of inputs of the introduction (in [L] either, since
  they are not used for the main statements).
- Fix: both the status and the inputs are now stated where the two results are mentioned,
  and the "Quoted theorems" bullet says that these inputs are used only there.

### D4. Outstanding checks flagged in [L] not mentioned (missing caveat)
- Where: `s11_final.tex`, "What the results are conditional on".
- [L] has draft notes saying: numbers of CCIT15 and Tel12 are from the arXiv versions,
  journal versions not compared; the formulas of the PPZ15 action are written from the
  way they are used and are to be compared with the print, including the sign in the
  symplectic condition; Barnes's second lemma was confirmed only through secondary
  statements; pointers to Mac95 and Sta99 are from memory; no theorem numbers pinned for
  the Frobenius formula, the Murnaghan-Nakayama rule, Frobenius' transposition identity,
  string and dilaton equations with pulled-back classes; a reference for the axioms of
  the virtual class is to be added; HK10, CI18, ZZ08 were not consulted, LP18 and
  BFGW21 were read in part only.
- Fix: new bullet "Checks which [L] marks as outstanding", listing these in eight
  sentences, without reproducing the notes. The sentence "No hypothesis of the main
  statements is left unproved in [L]" is kept, now restricted explicitly to Theorems
  1.1, 1.2 and Corollaries 1.3, 1.4.

### D5. Summary of the evidence stronger than the appendix (stronger than proved)
- Where: `s11_final.tex`, "Finite-order evidence".
- [A] before: "each check ... is accompanied by controls"; "The checks cover each link of
  the argument".
- [L, Appendix A]: defines a control but several items have none; Section A.7 lists
  controls that do not discriminate (rank-one sign controls at odd genus; first order in
  z cannot separate R from R*), and says that the class-level identity was paired only
  with psi classes and powers of lambda_1.
- Fix: "many are accompanied by controls"; "The checks concern the following links";
  the two non-discriminating controls and the limited class-level pairing are stated;
  "at two weight pairs" added to the per-colouring item.

### D6. "Quoted theorems" bullet imprecise (wording)
- Where: `s11_final.tex`.
- [A] before: Theorem 1.1 "rests in addition on" the mirror theorem etc. and "the
  special-function inputs of Lemma 6.4" (read literally: also on Mariño-Vafa and Barnes).
- [L, Sections 1.3, 1.5]: the frame formula uses neither the Mariño-Vafa formula nor
  Barnes's second lemma; Barnes enters only through `t_cf = t_c` in the corollaries.
- Fix: the bullet now gives the three dependency sets separately, and names
  BFGW21 Lemma 4.2 as used only for statements involving A.

### D7. Section 3 opening omits [L, Lemma 3.14] (wording)
- Where: `s03_rows.tex`, first paragraph: "These results are all that Section 7 uses".
- [L]: the list includes the closed forms of rows 1 and 2 (Lemma 3.14), which
  Section 7 uses; [A] Section 7 says so itself.
- Fix: Lemma 3.14 added to the sentence.

### D8. String identity quoted without its hypothesis (wording)
- Where: `s08_vanishing.tex`, sentence quoting [L, Corollary 8.13].
- [L]: "Let g >= 1 and U non-empty with positive parts."
- [A] before: "for g >= 1" only (X_1 of the empty multiset is not defined).
- Fix: hypothesis added.

## 4. Checked and found in agreement (no change)

- Introduction: Theorems 1.1, 1.2, Corollaries 1.3, 1.4 are word for word those of [L].
  The list of inputs has the same seven items with the same pinpoint references. "What
  is not claimed" is identical. The paragraph on methods not used (no topological
  recursion on the mirror curve, no remodelling, the only relation so called is the
  genus-zero psi-class relation, none of the conifold results of [BS], no finite-order
  computation as a proof step) is complete relative to [L].
- Conditional hypotheses kept: Theorem 10.8 and Corollary 10.9 of [L] keep the
  hypothesis M_irr(h) in [A], with the sentence that Lemma 10.7(b) discharges it.
- Statements given "in part", "condensed", "combined": each restriction is true as
  written. Examples: Lemma 8.10 + 8.12 combined is restricted to g >= 1 in its first
  half, which avoids the separate genus-0 definition; Theorem 5.1(a) omits the explicit
  pairing formula but keeps the F-derivative of the edge; Theorem 4.1 omits the
  structural form of (A3); Corollary 3.16 omits (Form).
- [BFGW21, Lemma 4.2] is named where it is used (after Definition 6.5 of [A]), with the
  statement that the frame formula does not use it.
- The constant `GW_{g,0}` is attributed in [A] to Theorem 5.1(b) where [L] points to
  Proposition 5.8; part (b) reproduces the formula of that proposition, so this is exact.
- Sketches: all 20 are headed "Sketch", "Idea of the proof" or "Proof, given the results
  above" (the last one, for the frame formula, is the complete proof of [L]). Each
  follows the route of [L].
- Pointers: all 122 match `main.aux` (theorem-like numbers, equation numbers (5.1) and
  (9.4), section numbers 1.5, 2.1, 3.2-3.4, 11.1, A.7).
- Front matter: author "The GENUS collaboration" (author line changed on 2026-11-23, both PDFs
  rebuilt and the title block rechecked); the introduction ends with
  `\subsection*{Acknowledgements}` and the prescribed sentence; no draft-note macro, no
  name of a person as author, no e-mail address in `main_summary.tex` or `summary/*.tex`.

## 5. Not checked, or not fixable here

- The cited literature was not consulted. The comparison is [A] against [L] only; the
  correctness of [L] itself was not refereed.
- The numerical values in the evidence appendix were not recomputed; only the summary in
  [A] was compared with the text of the appendix.
- The comparisons with [BS] ("the display of [BS] at genus 2", "the row recursion of
  [BS]") are carried over from [L], where they are tracked by labels in draft notes;
  [BS] was not opened.
- The bibliography of [A] is produced from `refs.bib`, which must not be edited. Four
  entries print working notes: "Publication data to be added" ([BS]), "Release to be
  added" ([DLMF]), "Source of Barnes's second lemma" ([Bar10]), and "numbers cited from
  this version" (several). These appear in the reference list of `main_summary.pdf`.
  Removing them needs a change of `refs.bib` or a separate bibliography file.
- One overfull box of 2.09pt remains (`s02a_setup.tex`, lines 126-128); it was there
  before and is under the 5pt limit.

## 6. Files changed

- `main_summary.tex`: abstract (two sentences), `\thanks`.
- `summary/s01a_intro.tex`: "About this version" (one clause).
- `summary/s01b_intro.tex`: last sentence of the plan.
- `summary/s03_rows.tex`: one sentence of the opening paragraph.
- `summary/s08_vanishing.tex`: hypothesis of the string identity.
- `summary/s11_final.tex`: rewritten (25 lines).

No other existing file was modified. The long version and `refs.bib` are untouched;
`main.pdf` and `main.aux` keep their time stamps.

## 7. Build

`latexmk -pdf -interaction=nonstopmode -halt-on-error main_summary.tex`: exit 0, 32 pages.
No undefined reference or citation, no multiply defined label, no BibTeX warning, one
overfull box (2.09pt). All files of `summary/` are under 250 lines (largest 225).
