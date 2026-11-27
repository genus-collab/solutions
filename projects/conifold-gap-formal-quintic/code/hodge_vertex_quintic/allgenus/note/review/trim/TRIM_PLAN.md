# Trim plan: one paper of 110-120 pages (key trim/plan)

Inputs: `review/REVIEW_SYNTHESIS_2.md` s.9 and s.10 item 6, `design/DECISIONS.md` D8 (and D4, D10). Continued in
`TRIM_PLAN_2.md` (per-file page table, the full list of references to rewrite, checks for the assembler).

## 0. Measurement

The current sources (162 pp) were compiled in the scratch folder of this task, with a `\phantomsection\label{zzf:<file>}` put at
the top of every section file of the scratch copy only (the note itself is unchanged); pages were read from `main.aux` and
`main.toc`. A block's size in pages is its number of characters divided by 5243 (the mean over the section files). Pages per section:

| part | pages now | first-last page | subsections (start page) |
|---|---|---|---|
| title, abstract, contents | 1 | 1 | |
| s01 intro | 8.5 | 1-9 | 1.1 (2), 1.2 (3), 1.3 (4), 1.4 (6), 1.5-1.6 (8), 1.7 (9) |
| s02 setup | 23 | 9-32 | 2.1 (10, 8 pp), 2.2 (18), 2.3 (20), 2.4 (23), 2.5 (26), 2.6 (30), 2.7 (32) |
| s03 gauge | 23 | 32-55 | 3.1 (33, 5.4 pp), 3.2 (38), 3.3 (40), 3.4 (45), 3.5 (48), 3.6 (49, 6 pp) |
| s04 transport | 13 | 55-68 | 4.1 (56, 8 pp, of which lem:tra-conv 2.8), 4.2 (64), 4.3 (65) |
| s05 reduction | 3.5 | 68-71 | 5.1 (68), 5.2 (70) |
| s06 model | 20.5 | 71-91 | 6.1 (72, 7 pp), 6.2 (79), 6.3 (82), 6.4 (84), 6.5 (87, 4 pp), 6.6 (91) |
| s07 comparison | 6.5 | 92-98 | 7.1 (92), 7.2 (96) |
| s08 connection | 15.5 | 98-113 | 8.1 (99), 8.2 (102), 8.3 (105, 7 pp), 8.4 (112), 8.5 (113, 1 p) |
| s09 gamma | 11 | 114-124 | 9.1-9.2 (114, 5 pp), 9.3 (119), 9.4 (123, 1.5 pp) |
| s10 main | 10 | 125-134 | 10.1 (125, 5 pp), 10.2 (130), 10.3 (132), 10.4 (134) |
| A formal stationary phase | 11 | 135-146 | |
| B signs | 5 | 146-151 | |
| C evidence | 10 | 151-160 | tables about 3 pp |
| references | 2 | 161-162 | |

Remarks: 38 labelled remarks, about 9.7 pp in total (s.4 below).

## 1. Page budget (target 117; ceiling 120)

| part | now | target | of which |
|---|---|---|---|
| front matter | 1 | 1 | abstract left as is (D12; title and abstract are the authors') |
| s01 intro | 8.5 | 5 | inputs listed once here (1.4); plan 1.7 rewritten for the new layout |
| s02 setup | 23 | 13 | 3.5 pp moved to A1; rmk:set-qrrderiv (1.1 pp) and rmk:set-reversion cut |
| s03 gauge | 23 | 14 | about 4 pp of 3.1 moved to A2; rmk:gau-Uremarks cut |
| s04 transport | 13 | 8 | proof of lem:tra-conv (about 2.5 pp) moved to A3; rmk:tra-consistency cut |
| s05 reduction | 3.5 | 2.5 | rmk:red-hyp, rmk:red-leaf cut |
| s06 model | 20.5 | 15 | prose only (rmk:mod-hock shortened, see s.4) |
| s07 comparison | 6.5 | 5 | rmk:flh-dilfactor cut, rmk:flh-fail shortened |
| s08 connection | 15.5 | 10.5 | sec:con-finite (1 p) and four remarks cut |
| s09 gamma | 11 | 8 | six remarks cut (1.1 pp), two shortened |
| s10 main (10.1 + remarks) | 5.7 | 4 | CGYZ26 dictionary moved to B; rmk:mai-analogy cut |
| s11 Complements (old 10.2-10.3) | 4.3 | 3 | rmk:mai-tau, rmk:mai-vertexterm, rmk:mai-apriori cut |
| A1 large-radius cyclic structure (new) | - | 3 | from s02.1 |
| A2 the cyclic form of PF4 (new) | - | 4 | from s03.1 |
| A3 u-adic convergence (new) | - | 2.5 | from s04.1 |
| A4 formal stationary phase (old A) | 11 | 9 | prose only; rmk:fsp-analytic cut |
| A5 signs (old B) | 5 | 5 | gains the CGYZ26 dictionary (+0.5), tightened elsewhere |
| A6 evidence (old C) | 10 | 3 | the six tables with one short paragraph each |
| references | 2 | 2 | |
| **total** | **162** | **117** | |

Structural cuts and moves save about 15 pp (remarks 5.5, sec:con-finite 1, appendix C 7, shortened remarks 1, the dictionary
moved and shortened 0.5); the other 30 pp come from tightening prose (s.6). Moves to appendices save nothing by themselves.

## 2. New layout, files and the lines of main.tex

Order of the appendices: by first use in the body, then the old appendices.

```
\input{sections/s01_intro}
...                                   (s02-s09 unchanged)
\input{sections/s10_main}             % s10_main_2 no longer inputs s10_main_3
\input{sections/s10_main_3}           % new final section "Complements" (old 10.2-10.3); s10_main_3 inputs s10_main_4
\appendix
\input{sections/a1_cyclicLR}          % new (T2): proof of prop:set-cyclicLR
\input{sections/a2_cyclicPF4}         % new (T3): the cyclic form of PF4
\input{sections/a3_uconv}             % new (T4): proof of lem:tra-conv
\input{sections/a_formalsp}
\input{sections/b_signs}
\input{sections/c_evidence}
\bibliographystyle{amsalpha}
\bibliography{refs}
```

The Complements section stays in the files `s10_main_3`, `s10_main_4` (group T1), so no new file is needed for it: T1 puts
`\section{Complements}\label{sec:cmp}` at the top of `s10_main_3`, keeps `sec:mai-nat`, `sec:mai-univ` as its subsections and
removes `\input{sections/s10_main_3}` from the end of `s10_main_2`. The main-file comment "this file inputs only the first part of
every section" must then be adjusted by the assembler. New section labels (prefix of the owning group): `sec:set-lrproof` (A1),
`sec:gau-cycproof` (A2), `sec:tra-convproof` (A3), `sec:cmp` (Complements). Appendix letters change (old A, B, C become D, E,
F): every `.tex` header comment that names an appendix by letter ("Appendix A, part 2 of 5", ...) is rewritten without the letter.
No reference in the text uses a letter (grep for `Appendix~[A-C]` is empty); all go through `\cref`.

## 3. Moves (the group owning the source file creates the target file)

M1 (T2) `s02_setup_2` l.23-83 and `s02_setup_3` l.3-49 -> `a1_cyclicLR.tex`, in this order: lem:set-grading, lem:set-unitfrob,
lem:set-periods, lem:set-birkhoff, lem:set-wronskian, lem:set-order3, the proof of prop:set-cyclicLR (i)-(ii) with its closing
paragraph (l.49). Stay in s02.1: the paragraph before prop:set-cyclicLR, the statement of prop:set-cyclicLR, the sentence on the
superscript LR (l.21), one sentence "(i) and (ii) are proved in \cref{sec:set-lrproof}; (iii) after \cref{def:set-rows}", and
cor:set-cyclicPF5 with its proof. External reference into the moved block: lem:set-periods from s02_setup_5 l.58 (label kept, fine).

M2 (T3) s03.1 (`s03_gauge` l.18-112, `s03_gauge_2` l.1-126) -> `a2_cyclicPF4.tex`. Stay in s03.1 (about 1.5 pp): the subsection
intro (shortened), def:gau-conc with eq:gau-Mcyc (needed for the statements and cited from s02_setup_5, b_signs, c_evidence),
lem:gau-confexp with its proof (cited 8 times from s02, s03.2, s03.6, s10), and the statements of thm:gau-cyc (eq:gau-cyc) and
cor:gau-cycsys (eq:gau-system), each followed by "Proof in \cref{sec:gau-cycproof}". Move, in this order: lem:gau-const,
lem:gau-adj, lem:gau-conc, lem:gau-transv, lem:gau-iso, lem:gau-sg (statements and proofs), then the proof of thm:gau-cyc and the
proof of cor:gau-cycsys. Labels kept; external references to lem:gau-conc (s02_setup l.16, s02_setup_5 l.58, s03_gauge_7 l.70)
and lem:gau-adj (s02_setup l.16, s02_setup_5 l.58) need no change.

M3 (T4) lem:tra-conv, statement (`s04_transport_3` l.4-32, about 1.2 pp) and proof (l.34-64), -> `a3_uconv.tex`, unchanged
and with its label; s04.1 keeps two or three sentences saying what the lemma gives and pointing to \cref{lem:tra-conv} in
\cref{sec:tra-convproof} (lem:tra-conv is cited 20 times from s02, s04-s07; all references stay valid). rmk:tra-consistency
(l.68-71) is cut, not moved. If T4 prefers to keep the statement in s04.1, the budget of s04 rises by 1 p and must be found in s04.

M4 (T1) the CGYZ26 dictionary of rmk:mai-constant (`s10_main_2` l.43 ff.) -> `sec:sgn-lit` (`b_signs_2`, row or paragraph next
to tab:sgn-lit, which already refers to rmk:mai-constant at l.99: change that pointer to the new paragraph). rmk:mai-constant
stays in s10.1 as a short remark (D10): the comparison of c_g with [CGYZ26, Conj. 1.5], the values pointer to rmk:set-cg, one
sentence "the dictionary is in \cref{sec:sgn-lit}", and its `\draftnote` (s10_main_2 l.55). Its references from def:set-confexp
(s02_setup_5 l.64, inside a definition: keep the label so the definition is untouched), s01_intro_2 l.74, s01_intro_3 l.23, l.31
stay valid. The reference to sec:evi-calib in it (l.54) becomes sec:evi (s.4, C-cut).

M5 (T1) old 10.2-10.3 become the Complements section (s.2). rmk:mai-graphlocal and rmk:mai-necessity (old 10.4, about Theorem
A) move to the end of 10.1, after cor:mai-gap; subsection sec:mai-rmk is dissolved (its only reference: s10_main l.4, the section
intro, rewritten). Forward references from the body into Complements stay valid (cor:mai-natfact from def:set-hat, s07_comparison_2
l.28; lem:mai-Ghat from thm:gau-R0 is a pointer only, Ghat is defined in thm:gau-R0 itself).

## 4. Cuts

Remarks cut (with their `\draftnote`s where noted). Every reference to a cut label is listed with its file and line in
`TRIM_PLAN_2.md` s.2; the owner of the referring file rewrites it (drop the parenthetical, or point to a kept item).

| label | file (owner) | pp | draftnote inside (goes) | references to rewrite (file:line) |
|---|---|---|---|---|
| rmk:set-qrrderiv | s02_setup_4 (T2) | 1.13 | l.77 | s01_intro_3:10, :23 (x2), :32; s02_setup:16; s02_setup_4:58 (in rmk:set-qrrsign) |
| rmk:set-reversion | s02_setup_8 (T2) | 0.19 | | none |
| rmk:gau-Uremarks | s03_gauge_4 (T3) | 0.29 | | none (it cites sec:evi-gauge) |
| rmk:tra-consistency | s04_transport_3 (T4) | 0.30 | l.70 | none |
| rmk:red-hyp | s05_reduction (T4) | 0.21 | | none |
| rmk:red-leaf | s05_reduction (T4) | 0.27 | | s10_main_3:74 (in rmk:mai-vertexterm, cut); c_evidence_3:32 (cut text) |
| rmk:flh-dilfactor | s07_comparison_2 (T4) | 0.18 | | s07_comparison_2:105 (in rmk:flh-fail) |
| rmk:con-problem | s08_connection (T6) | 0.14 | | none |
| rmk:con-model | s08_connection (T6) | 0.11 | | s08_connection:25 (in rmk:con-problem, cut) |
| rmk:con-points | s08_connection_2 (T6) | 0.16 | | none |
| rmk:con-why | s08_connection_8 (T6) | 0.30 | | s10_main_3:47 (T1); c_evidence_2:38 (T6) |
| rmk:gam-lmi-words | s09_gamma_2 (T6) | 0.14 | | none |
| rmk:gam-K8 | s09_gamma_3 (T6) | 0.08 | | none |
| rmk:gam-sevenvar | s09_gamma_4 (T6) | 0.46 | l.88 | b_signs_2:98 (tab:sgn-lit, T1); c_evidence_3:56, :58 (cut text) |
| rmk:gam-blind | s09_gamma_4 (T6) | 0.20 | | c_evidence_3:52 (cut text) |
| rmk:gam-Mcomplex | s09_gamma_4 (T6) | 0.08 | l.107 | c_evidence_3:48 (x2), :58 (cut text) |
| rmk:gam-general | s09_gamma_4 (T6) | 0.13 | | c_evidence_3:54 (cut text) |
| rmk:mai-apriori | s10_main_3 (T1) | 0.24 | | s07_comparison_2:28 (T4); c_evidence:72 (T6) |
| rmk:mai-tau | s10_main_3 (T1) | 0.16 | l.68 | c_evidence:74 (T6) |
| rmk:mai-vertexterm | s10_main_3 (T1) | 0.16 | | none |
| rmk:mai-analogy | s10_main_4 (T1) | 0.26 | l.79 | s01_intro:26 (T1) |
| rmk:fsp-analytic | a_formalsp (T6) | 0.06 | | a_formalsp:6 (T6) |

Subsection cut: sec:con-finite (`s08_connection_8` l.39-72, T6, 1.05 pp), with lem:con-N and thm:con-conn (their only references
are inside it). References to sec:con-finite: s08_connection_7:52 (after thm:con-Pb), s08_connection_8:19 (proof of thm:con-glue,
the sentence "Neither the finite route of ... nor": keep "the degree bound prop:con-deg is not used"), s09_gamma_3:92 (after
thm:gam-D), c_evidence_2:17 (table tab:evi-windows, T6: drop the pointer). None of them is a step of a proof.

Appendix C cut (T6): keep sec:evi with a two-sentence intro (the present first sentence of l.5 "Nothing in this appendix is a step
of a proof ..." stays), and the subsections that carry tables: sec:evi-graphs (tab:evi-g2..g5), sec:evi-normal (tab:evi-windows),
sec:evi-model (tab:evi-model), each table with at most one short paragraph. Remove the subsections sec:evi-calib, sec:evi-hodge,
sec:evi-gauge, sec:evi-transport, sec:evi-univ, sec:evi-lmi. References to them outside appendix C, all rewritten to
`\cref{sec:evi}` or dropped: sec:evi-hodge s01_intro:26; sec:evi-univ s01_intro_2:25, s10_main_4:61, :69, :74, b_signs_2:60;
sec:evi-calib s01_intro_3:30, s10_main_2:54; sec:evi-gauge s03_gauge:15 (and s03_gauge_4:97, inside rmk:gau-Uremarks, cut);
sec:evi-transport s05_reduction:118 (inside rmk:red-leaf, cut); sec:evi-lmi s09_gamma_4:102 (in rmk:gam-necessity, kept: rewrite),
:93 and :112 (inside cut remarks). The full appendix stays in `archive_full_162pp/` (an ancillary file is an author decision).

## 5. Remarks kept (D8), target at most about 2.2 pp in total

rmk:int-constant (s01, shorten to 3 lines + its draftnote l.74), rmk:set-qrrsign (used in the proof at s02_setup_4:73),
rmk:set-cg (values of c_g), rmk:mod-hock (shorten to two sentences + its draftnote l.99; kept because of D4, see s.7),
rmk:mod-offK (necessity of K), rmk:flh-fail (necessity; shorten to about 0.3 pp), rmk:con-genus (cited from s07_comparison_2:96;
shorten), rmk:con-710 (the weight lemma: one sentence), rmk:gam-necessity (cited in the proof of thm:gam-lmi, s09_gamma_2:18;
shorten to about 0.3 pp), rmk:gam-prior (attribution; shorten, keep its draftnote l.117), rmk:mai-constant (short, s.3 M4),
rmk:mai-annih (cited inside def:tra-dressing, s04_transport_4:14, so kept to leave the definition untouched; two sentences),
rmk:mai-stored (necessity: a non-canonical normalisation), rmk:mai-graphlocal (cited after thm:mai-main and from s01),
rmk:mai-necessity, rmk:fsp-classical (cited in a proof, s09_gamma_4:34). Sum after shortening: about 2.3 pp.
