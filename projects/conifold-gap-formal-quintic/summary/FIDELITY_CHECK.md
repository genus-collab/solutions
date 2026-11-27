# Fidelity check of the abridged version against the long version

Date: 2026-11-20. Reference: `main.tex` + `sections/*.tex` (127 pp, `main.aux` of 17:12, newer than
every long source). Checked: `main_summary.tex` + `summary/*.tex` (32 pp).
The long version, `refs.bib` and all other existing files were not modified.

## 1. What was compared

| Abridged | Long | How |
|---|---|---|
| Abstract, front matter, Section 1 (`s01_intro*`) | abstract, `s01_intro`, `_2`, `_3`; `b_signs_2` (E.9); `s10_main_2` (Remark 10.5) | read in full, both sides |
| Section 5 (`s05_model*`) | Section 6, all ten parts, with every proof | read in full, both sides |
| Section 7 (`s07_connection*`) | Section 8, all parts, with every proof; statements of Appendix D (`a_formalsp*`) | read in full; Appendix D statements only |
| Section 8 (`s08_gamma`) | Section 9, all four parts, with every proof | read in full, both sides |
| Section 2 (`s02_setup*`) | Section 2 (8 parts, proofs), Appendices A, B | delegated reader; every reported item re-verified here |
| Section 3 (`s03_gauge*`) | Section 3 (8 parts, proofs), Appendix B | delegated reader; every reported item re-verified here |
| Sections 4, 6, 9, 10 | Sections 4, 5, 7, 10, 11, Appendices C, F | delegated reader; every reported item re-verified here |
| Pointers `\Lr`, `\LL` | `main.aux` | script: all 359 entries of `longrefs.tex` equal `main.aux`; all 178 labels used exist |

About 800 statements, formulas and claims were compared (intro about 60, Section 5 about 55,
Section 7 about 70, Section 8 about 40, Sections 2, 3, and 4/6/9/10 about 180, 165 and 225).

## 2. Result in one paragraph

No abridged theorem, proposition, lemma, corollary, definition or displayed formula is stronger than,
or different from, its long counterpart. Two items were wrong (one item pointer, one index slip in a
proof); seven places lacked a caveat that the long version states; eleven were pointer or
undefined-notation matters. All twenty were fixed. Sections 5, 7 and 8, whose sketches the writer had
drafted from truncated proofs, agree with the full long proofs step by step; nothing had to be
corrected there except one conditionality sentence and one pointer.

## 3. Discrepancies and fixes

Severity: W = wrong, S = stronger than proved, C = missing caveat, T = wording / pointer / notation.

| # | Sev | Abridged location | Long version | Change made |
|---|---|---|---|---|
| 1 | W | `s01_intro_2`, Inputs: "[AC96, (1.7)] is used in the dilaton lemma (Lemma (i))" | In the abridged merged lemma the dilaton equation is item (ii) (item (i) is homogeneity); the long numbering (i) was copied | now "the dilaton equation of Lemma (ii)" |
| 2 | W | `s03_gauge_2`, proof of the symplecticity identity: "its entry (2,1) at u=0 gives t(z)=cn(z)cn(-z)", "its" = T = t(z)1 | Long proof takes the entry of S = eta T: S_21(0,z) = 5 cn(-z) cn(z) = 5 t(z); the entry (2,1) of the scalar matrix T is 0 | "by (R1), r = r~ and eta_{a1} = 5 delta_{a2}, the entry (2,1) of eta T at u=0 gives 5t(z) = 5cn(-z)cn(z)" |
| 3 | C | `s01_intro`, Remark on the constant: comparison with [CGYZ26] without the admissibility caveat | E.9 and Rem. 10.5: F_g^cf is the conifold expansion of [CGYZ26, Def. 1.2] only for an admissible pair, "asserted in [CGYZ26, Ex. 1.3, Rem. 1.4] and not proved here", under the inputs of Theorem C | sentence added, with pointer [L, Appendix E.9] |
| 4 | C | `s02_setup_2`, outer graph sum, hatted term: hatted unit omitted, `s^{-1}.TFT` undefined at that point | Def. of the hatted theory: unit 1^ = H^0^ = s1, TFT^_{g,n} = s^{-(2g-2+n)} TFT_{g,n}; without 1^ the relation G = u^{2g-2} T fails | both written out |
| 5 | C | `s03_gauge`, quintic cyclic system (i): "the forward relations of Prop. (iii)" (stated there for all a in Z/5, for the rows), X_a undefined | Forward relations only for a = 0..3, for elements X_a of N((z)) x Q(zeta); the closing relation a = 4 is separate, and for the rows is obtained later | range, the ring of X_a and "for the rows it is obtained in Theorem (normalisation series)" added |
| 6 | C | `s06_comparison`, opening: R' presented as an instance of the symplectic O | "(unit) involves only the constant term ... The definitions do not decide whether R' is symplectic" | clause added: condition on the constant term only; R' not assumed symplectic; T_Gamma uses the edge of the outer-sum definition |
| 7 | C | `s09_main`: "the graph sum is the potential by (T)" | "by (T), with the normalisation (QRR)" | "(T), with the normalisation (QRR)" |
| 8 | C | `s09_main`, after Theorem A: "the inputs enter as described in Section 1.4"; Corollary B without conditionality | "Theorem A holds modulo (Q), (QRR), (M), (G0), (T) and [AC96,(1.7)]"; "its conditionality is that of Theorem A" | explicit list for Theorem A and Corollary B |
| 9 | C | `s08_gamma`, after Theorem C: "(QRR) enters through the glue and the colour covariance" | also: the rows and their WKB form rest on (M), (G0) and flatness; the R-matrix of (T) enters only through R(0)=1 and (QRR) | sentence completed ("the value of (QRR) ...; the rows ... rest on (M), (G0) and the flatness part of (QRR) ...") |
| 10 | T | `s02_setup`: `e^{T(q)}` used, never defined (T contains the formal log q) | Def. 2.2: e^T := q exp(I~_1/C_0) | definition added |
| 11 | T | `s02_setup`: I(q,z) pointed to Theorem 2.4 only | I(q,z) is defined in Def. 2.2 | pointer [L, Def. 2.2, Thm 2.4] |
| 12 | T | `s02_setup`: sketch of the large-radius structure covers (iii) but points to Appendix A only | (iii) is proved after Def. 2.10 | pointer extended |
| 13 | T | `s02_setup`: concomitant identity with P on both sides pointed to Lemma B.3 | B.3 has the adjoint P*; the stated form needs Lemma B.2 (P* = P) | Lemma B.2 added to the pointer |
| 14 | T | `s03_gauge_2`, canonical gauge: "colour covariant" undefined (Rot_i dropped) | Rot_i Pc = Pc, equivalent to gradedness (Prop. 3.14 (iv)) | explicit formula zeta^{i(a-b)} Pc_ab(u, zeta^{-i} z) = Pc_ab(u,z) |
| 15 | T | `s03_gauge_2`, sketch of the normalisation theorem: "the hatted gauge" undefined | G^ := Lambda_s Pc(u,s^4 z^) Lambda^{-1}, exponents mu(b,k) >= 0 integers | "in the hatted variables the entries of Pc carry non-negative integer powers of u only" |
| 16 | T | `s03_gauge`, sketch of Theorem U: the degree argument for (c) is in none of the pointed statements | Lemma 3.13 | added to the pointer |
| 17 | T | `s06_comparison`, last remark: "(A) only if cn = 1" pointed to Remark 7.8 only | that sentence is in Remark 11.7 (stored gauge) | pointer extended |
| 18 | T | `s10_evidence`: "the one-point law" undefined; "was found to fail" | Cor. 6.38; Appendix F.3 has one example (psi^5 H, h = 3, v = 1) | law written out; "fails in an example"; pointer to Cor. 6.38 |
| 19 | T | `s04_transport_2`: "edge graph" undefined; `t_alpha` undefined in the sketch of Theorem R | Def. of the flat sum; Lemma C.1 (b2) | both defined in place |
| 20 | T | `s05_model`, conventions: q_phi defined, Def. 6.5 not pointed | Def. 6.5 (ancestor vector of a leg) | added to the pointer |

Severity count: W 2, S 0, C 7, T 11.

## 4. Checks that found no discrepancy

- **Front matter.** Author "The GENUS collaboration" (author line changed on 2026-11-23, both PDFs rebuilt
  and the title block rechecked); the introduction ends with
  `\subsection*{Acknowledgements}` and the prescribed sentence; `\showdraftnotesfalse`, no draft note,
  no personal name as author, no email address in any abridged file.
- **Inputs.** (Q), (QRR), (M), (G0), (T) and [AC96, (1.7)] are complete relative to the long
  Section 1.4: (QRR) = flatness and the value at q = 0, with both lists of sources; (T) = (i) the CohFT
  axioms in all genera with the divisor equation and (ii) Teleman's reconstruction; (Q) lists every
  quoted statement of the long (Q) (the same items of the four sources), says that symplecticity,
  (2-32) and the edge condition (4-1) are proved for the model curve in [L, Prop. 6.9], and that
  (2-33) is not used. Omitted relative to the long (Q) are only reading conventions, for which the
  abridged text points to [L, Section 1.4].
- **What is not claimed.** All six items are present verbatim, including the formal change of frame
  (no analytic continuation) and the quintic threefold; they are repeated at the end of Section 10.
- **Conditionality in the body.** Theorem A, Corollary B and Theorem C are stated as conditional;
  Theorem Gamma, Theorem U and Theorem P(a) as using no quoted input, as in the long version;
  Theorem P(b) "does not use (QRR)", as in the long version. Finite computations appear only in
  Section 10, introduced by "Nothing in this paragraph and the next two is a step of a proof";
  all ranges, graph counts, primes and frames equal the long tables.
- **Section 5 against the long Section 6.** Conventions, the TR theorem with the factor (-5)^{h-1},
  the leg dictionary (rho, the legs, q of Phi, of -log p, of the Euler-dilaton vector), the log
  coefficients, the lattice theorem (abridged (a), (b), (c) = long (b), (c), Cor. 6.24), L-regularity,
  the Gap Theorem (signs of the period and of F^TR_h recomputed), the variational formula, Theorem V
  with F^TR_0 = (1/2) v^2 log v and F^TR_1 = -(1/12) log v, and the proof of the Model Lemma: all
  agree. Each sketch follows the long proof (Laplace series and T-system uniqueness, kappa = 4 then
  scaling; [ABDKS25, (7.9), (7.10)]; induction on the number of log legs; local formula for
  omega_{1,1}).
- **Section 7 against the long Section 8.** The saddle series, the 1-form, W^, a^, Phi^SP, the
  residue V, Theorem P(a) (i)-(iv), the uniformising map, the reduction, the model identity, the
  two functionals, Theorem P(b) and its proof (order u^1), the glue and its corollary: all agree.
  The sketch of P(a)(i) is the long argument (shift identity, translation by z, functional
  equation). The summary of Appendix D in Section 7.1 agrees with Definitions D.1, D.3, D.20 and
  Lemmas D.4, D.5, D.7, D.11, D.14, D.17, D.19, D.21, D.22 as stated.
- **Section 8 against the long Section 9.** The delta lemma, LMI data, Theorem LMI (i)-(iii) and its
  three-step sketch, Theorem Gamma (i)-(iv) and the proof of (ii), Theorem C and its proof: all agree.

## 5. Observations not acted on

- In both versions the corollary on the fifth-order equation, part (b), reads
  "e^{u_i/z} z L Phi_i solves PF5"; the factor z is harmless (a constant) but is probably not
  intended. It is in the long source (`sections/s02_setup_3.tex`), which is read-only here; the
  abridged text copies it.
- One connecting remark of the abridged version has no counterpart sentence in the long version
  ("the product of the coefficients, which does not depend on the frame, fixes tau: this is the
  origin of the natural form", Section 3); it is a consequence of two long propositions and was kept.
- BibTeX reports one warning, "missing pages in ABDKS25"; it comes from `refs.bib`, which was not
  to be touched, and was present before this check.

## 6. What could not be checked

- The cited literature: statement numbers and content of the external sources were not opened.
- The mathematics of the long version itself: it was read for fidelity, not re-derived (a few
  constants and signs were recomputed in passing).
- Sections 2, 3, 4, 6, 9 and 10 were compared by three delegated readers working in parallel.
  Each discrepancy they reported was re-verified against the long source before it was fixed, but
  the absence of further discrepancies in those sections rests on their reading, not on a second
  independent one.
- Appendix D: statements only (not its proofs); Appendix E: only E.9 and the passages pointed to.
- The rendered PDF was checked through the build log and a text extraction (no "??"), not page by
  page by eye.

## 7. Build

`latexmk -pdf -interaction=nonstopmode -halt-on-error main_summary.tex` in the note folder:
exit 0, 32 pages, no undefined reference or citation, no multiply defined label, no overfull box
(0 reported), `main.pdf` and `main.aux` untouched (time stamps unchanged). Every summary file is
under 250 lines (largest: `longrefs.tex` 182, `s02_setup.tex` 127).
