# Trim plan, part 2 (key trim/plan)

Continues `TRIM_PLAN.md` (s.0-5 there). Line numbers are those of the current files in `note/sections/` (1-based).

## 6. Groups, files owned, budgets, instructions

Ownership is disjoint; a group edits only its own files, plus the new appendix file it creates. A reference into another group's
file that must change is reported to the assembler, not edited.

| group | files owned | creates | now (pp) | target (pp) |
|---|---|---|---|---|
| T1 | s01_intro*, s10_main*, b_signs* | none (Complements stays in s10_main_3, s10_main_4) | 8.5 + 10 + 5 = 23.5 | 5 + 7 + 5 = 17 |
| T2 | s02_setup* | a1_cyclicLR.tex | 23 | 13 + 3 = 16 |
| T3 | s03_gauge* | a2_cyclicPF4.tex | 23 | 14 + 4 = 18 |
| T4 | s04_transport*, s05_reduction*, s07_comparison* | a3_uconv.tex | 13 + 3.5 + 6.5 = 23 | 8 + 2.5 + 2.5 + 5 = 18 |
| T5 | s06_model* | none | 20.5 | 15 |
| T6 | s08_connection*, s09_gamma*, a_formalsp*, c_evidence* | none | 15.5 + 11 + 11 + 10 = 47.5 | 10.5 + 8 + 9 + 3 = 30.5 |
| (front matter, references) | main.tex, refs.bib (assembler) | | 3 | 3 |
| **total** | | | 162 | 117.5 |

Instructions common to all groups (from the binding rules):
* Every statement on the chain of Theorem A and Corollary B keeps its statement, hypotheses, label and a complete proof in
  sentences; a proof may move to an appendix, never be deleted. Compress routine steps (one sentence citing the lemma used),
  never skip a step. The conditionality ((Q) with its pins, (QRR), (M), (G0), (T), [AC96 (1.7)]) is unchanged.
* Lists of inputs are stated in s01 (1.4) and s10 (10.1) only. Elsewhere a proof names the inputs it uses in one clause, without
  restating the list; drop sentences of the form "this is the only place where an input enters ... neither ... nor ... is used"
  unless they are the only place that records which input a step uses.
* No repeated statements: a result stated in s01 and in the body keeps its full form in the body; the body does not restate what
  the introduction or an earlier section says, and section intros are at most a short paragraph.
* Keep every `\draftnote` except the seven inside cut remarks (TRIM_PLAN s.4). No project bookkeeping in `.tex` files (no
  review, gap, design or decision identifiers, no dates, comments included); header comments naming appendices by letter are
  rewritten without the letter.
* Each file stays below 250 lines; a new file starts with a one-line header comment and inputs nothing unless it is split.
* Report every reference to a removed label, and every reference you rewrote in another group's file's favour, to the assembler.

Per group:
* **T1.** s01: shorten 1.1-1.2 and 1.6 to about half, keep 1.3 (statements) and 1.4 (inputs, complete); rewrite 1.7 for the new
  layout (Complements, appendices A1-A6); drop the pointers to rmk:set-qrrderiv,
  rmk:mai-analogy, sec:evi-hodge, sec:evi-univ, sec:evi-calib. s10: M4, M5 of TRIM_PLAN s.3; cut rmk:mai-apriori,
  rmk:mai-tau, rmk:mai-vertexterm, rmk:mai-analogy; shorten the proof of cor:mai-gap where it repeats s02.3; the list of inputs of
  10.1 stays. b_signs: receive the dictionary in sec:sgn-lit; drop the rmk:gam-sevenvar pointer in tab:sgn-lit (b_signs_2:98)
  and the sec:evi-univ pointer (b_signs_2:60); tighten B.1-B.9 without dropping any convention used by a proof.
* **T2.** M1; cut rmk:set-qrrderiv (and its pointer in rmk:set-qrrsign, s02_setup_4:58, and in the section intro, s02_setup:16),
  rmk:set-reversion; tighten 2.1 (the genus-zero data), 2.3 (the graph sum: keep thm:set-GT and lem:set-pushins complete), 2.4, 2.5.
* **T3.** M2; cut rmk:gau-Uremarks; replace the pointer to sec:evi-gauge (s03_gauge:15) by sec:evi or drop it; tighten 3.3-3.6
  (3.6 is 6 pp: compress the routine colour-covariance and WKB steps).
* **T4.** M3; cut rmk:tra-consistency, rmk:red-hyp, rmk:red-leaf, rmk:flh-dilfactor; shorten rmk:flh-fail (drop its pointer to
  rmk:flh-dilfactor at s07_comparison_2:105); drop the pointer to rmk:mai-apriori after lem:flh-hat (s07_comparison_2:28), keeping
  the pointer to cor:mai-natfact; tighten 4.1 and 4.3.
* **T5.** prose only: 6.1 (7 pp: conventions and the TR dictionary) and 6.5 (Theorem V, 4 pp) are the main sources; shorten
  rmk:mod-hock to two sentences with its draftnote (s.7 item 2); keep rmk:mod-offK.
* **T6.** cut sec:con-finite, rmk:con-problem, rmk:con-model, rmk:con-points, rmk:con-why, rmk:gam-lmi-words, rmk:gam-K8,
  rmk:gam-sevenvar, rmk:gam-blind, rmk:gam-Mcomplex, rmk:gam-general, rmk:fsp-analytic; shorten rmk:con-genus, rmk:con-710 (one
  sentence), rmk:gam-necessity, rmk:gam-prior (keep its draftnote); sec:gam-rmk keeps rmk:gam-necessity and rmk:gam-prior; cut
  appendix C to its tables (TRIM_PLAN s.4); tighten 8.3 (Theorem P, 7 pp), 9.2 (LMI, 5 pp) and A.2-A.4.

## 7. Conflicts and points for the assembler

1. **D8 vs the removal list.** D8 keeps "the 7-variable chart", "LMI for complex M0 (one sentence)", "the finite route (short)" and
   "the uniqueness of tau (remark)". The task and REVIEW_SYNTHESIS_2 s.9 remove rmk:gam-sevenvar, rmk:gam-Mcomplex,
   sec:con-finite and rmk:mai-tau. This plan follows the task. If the authors want D8 kept, the cheapest fallback is one sentence
   each, without proof and marked as unproved where the source remark was (rmk:gam-Mcomplex and rmk:mai-tau are unproved,
   their draftnotes say so): about 0.3 pp in all.
2. **D4.** D4 makes Hock's formula the proof of record of the Gap Theorem; the current proof of thm:mod-gap uses prop:mod-period
   (iii) and cor:mod-trF, and rmk:mod-hock is the Hock route. The review does not list rmk:mod-hock; this plan keeps it, shortened.
3. **Statements containing references to cut items.** def:set-confexp cites rmk:mai-constant and def:tra-dressing cites
   rmk:mai-annih: both remarks are kept so that no definition changes. No kept statement cites any cut label (checked for all
   cut labels; the only statement citing one is thm:con-conn, cut itself).
4. **Optional further cut (not in the review's list).** With sec:con-finite gone, prop:con-deg and lem:con-u0 (s08.1, about 0.7 pp
   with prose) are used only in text and remarks (s01_intro:16, s08_connection:7, s08_connection_4:58, s08_connection_7:52,
   s08_connection_8:19, :35, c_evidence_2:34) and, by the explicit sentences of thm:con-glue and after thm:con-Pb, by no proof.
   Cutting them needs the assembler's approval; T6 does not cut them by default.
5. **Abstract.** Settled by the authors (D12); left unchanged and not counted as a saving.
6. Remarks kept: about 2.3 pp against D8's "about two pages"; T6 and T1 shorten first if the total is over.

## 8. References to cut labels, complete (file:line, enclosing block, owner of the file)

From a grep of `\cref`, `\Cref`, `\ref`, `\eqref` over `sections/*.tex` and `main.tex`, multi-label references included.
* rmk:set-qrrderiv: s01_intro_3:10, :23 (two), :32 (T1); s02_setup:16 (intro), s02_setup_4:58 (rmk:set-qrrsign) (T2).
* rmk:mod-hock (kept, listed for completeness): s01_intro_3:23 (T1), c_evidence_3:25 (T6, cut text).
* rmk:mai-analogy: s01_intro:26 (T1).
* rmk:red-leaf: s10_main_3:74 (rmk:mai-vertexterm, cut) (T1); c_evidence_3:32 (T6).
* rmk:flh-dilfactor: s07_comparison_2:105 (rmk:flh-fail) (T4).
* rmk:con-model: s08_connection:25 (rmk:con-problem, cut) (T6).
* rmk:con-why: s10_main_3:47 (text after cor:mai-natfact) (T1); c_evidence_2:38 (T6).
* rmk:gam-sevenvar: b_signs_2:98 (tab:sgn-lit) (T1); c_evidence_3:56, :58 (T6).
* rmk:gam-blind: c_evidence_3:52. rmk:gam-Mcomplex: c_evidence_3:48 (two), :58. rmk:gam-general: c_evidence_3:54 (T6).
* rmk:gam-prior (kept): s09_gamma_4:88 is inside rmk:gam-sevenvar (cut).
* rmk:mai-apriori: s07_comparison_2:28 (T4); c_evidence:72 (T6). rmk:mai-tau: c_evidence:74 (T6).
* rmk:fsp-analytic: a_formalsp:6 (T6).
* rmk:set-reversion, rmk:gau-Uremarks, rmk:tra-consistency, rmk:red-hyp, rmk:con-problem, rmk:con-points, rmk:gam-lmi-words,
  rmk:gam-K8, rmk:mai-vertexterm: no references.
* sec:con-finite: s08_connection_7:52, s08_connection_8:19 (proof of thm:con-glue), s09_gamma_3:92, c_evidence_2:17 (T6).
  lem:con-N: s08_connection_8:63, :69 (inside the cut subsection). thm:con-conn: none.
* sec:mai-rmk (dissolved): s10_main:4 (T1). sec:mai-nat, sec:mai-univ: s10_main:4 (the section intro: now point to sec:cmp).
* Removed appendix C subsections: see TRIM_PLAN s.4 (all outside references listed there).
* Pointers to moved items need no change (labels kept), except tab:sgn-lit's pointer to rmk:mai-constant (b_signs_2:99), which
  should point to the dictionary's new place in sec:sgn-lit.

## 9. Per-file pages now (start page, pages), from the scratch build

s01_intro 2 2 | s01_intro_2 4 2 | s01_intro_3 6 3 | s02_setup 9 4 | s02_setup_2 13 2 | s02_setup_3 15 3 | s02_setup_4 18 5 |
s02_setup_5 23 3 | s02_setup_6 26 2 | s02_setup_7 28 2 | s02_setup_8 30 2 | s03_gauge 32 3 | s03_gauge_2 35 3 | s03_gauge_3 38 4 |
s03_gauge_4 42 3 | s03_gauge_5 45 4 | s03_gauge_6 49 3 | s03_gauge_7 52 3 | s04_transport 55 3 | s04_transport_2 58 3 |
s04_transport_3 61 3 | s04_transport_4 64 4 | s05_reduction 68 3 | s06_model 71 3 | s06_model_2 74 1 | s06_model_3 75 3 |
s06_model_4 78 2 | s06_model_5 80 2 | s06_model_6 82 2 | s06_model_7 84 3 | s06_model_8 87 2 | s06_model_9 89 2 |
s06_model_10 91 1 | s07_comparison 92 2 | s07_comparison_2 94 4 | s08_connection 98 2 | s08_connection_2 100 2 |
s08_connection_3 102 1 | s08_connection_4 103 3 | s08_connection_5 106 2 | s08_connection_6 108 2 | s08_connection_7 110 2 |
s08_connection_8 112 2 | s09_gamma 114 2 | s09_gamma_2 116 3 | s09_gamma_3 119 2 | s09_gamma_4 121 4 | s10_main 125 2 |
s10_main_2 127 3 | s10_main_3 130 2 | s10_main_4 132 3 | a_formalsp 135 1 | a_formalsp_2 136 3 | a_formalsp_3 139 2 |
a_formalsp_4 141 2 | a_formalsp_5 143 3 | b_signs 146 2 | b_signs_2 148 3 | c_evidence 151 4 | c_evidence_2 155 2 |
c_evidence_3 157 5 (page differences of file starts; a file starting mid-page counts from that page).
