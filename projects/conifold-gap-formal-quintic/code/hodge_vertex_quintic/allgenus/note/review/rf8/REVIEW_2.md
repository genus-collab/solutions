# Review rf8, part 2: appendix C fixes, open points, effects on other sections, cross-section list

Continuation of `REVIEW.md` (verdict s.0, method s.1, appendix B fixes s.2.1).

### 2.2 Appendix C (line numbers before the fix)

1. **Normalisation slip** (`c_evidence_3.tex` l.12, `tab:evi-model` row 1): "topological recursion on `x = W_v` is `(-5)^(1-h)` times the
   model CohFT" is inverted: `thm:mod-tr` gives `omega_(h,n) = (-5)^(h-1) <..>` (TRM/TRRF: "TR x `(-5)^(1-g)` = model"). Corrected, both
   directions stated.
2. **R-matrix of the curve** (row 2): "is `M_0`" -> "of `x = kappa W_v` is `M_0(z/(kappa varsigma^4))`" (`prop:mod-Rtr`, eq:mod-Rtr; the
   check ran at `kappa = 1, 4, -1`, `varsigma = 1, 2`, G4_2 s.8, G4R s.1.4).
3. **Misattribution** (`c_evidence.tex` l.76): "[the second implementation] detected each of 26 deliberate errors". NF s.5: the 26
   mutants were caught by the genus-2/3 regression (23) and the genus-4 checks (3: `tau_4`, `tau_5`, orders `zhat >= 7`). Rewritten.
4. **Counts attributed to the wrong objects** (l.72): "1138 vanish by genuine cancellation, between 31 and 4525 contributions each"
   - NF s.6 census: of 3142 coefficients, 1655 lie in rows that receive pieces, 517 of these receive none at that order, 1138 cancel;
   "min 31, median 347, max 4525" is **per cancelling graph**. Rewritten.
5. **Genus-5 controls** (l.74): "every graph computed is singular" -> the 14 *edge* graphs with at most two edges (the vertex graph is
   "spoilt", not "singular"); the `tau + u^7` control was also restricted to `E <= 2` and spoils the vertex graph at `t_c^-1` only;
   "149 genus-4 graphs" -> "149 edge graphs". The uniqueness of `tau` mod `u^8` is now derived as the source does it (valuation
   `k + E - 1` by multilinearity; the lowest-order response of a one-edge graph is `c` times a `k`-independent non-zero number,
   measured at `k = 1, 7`), with the ansatz letter of `rmk:mai-tau`.
6. **"Exact in u"** (l.17) held only at genus 4-5 (NF s.5: no truncation of `tau`); at genus 2-3 `M` was taken mod `(zhat^7, u^4)`
   (V5, UR s.6.2), which suffices. Stated precisely; the regrouping sentence (which read as if the regrouping *produced* the
   Witten-Kontsevich vertices) rewritten with the splitting formula `E_R = E_M + M^-1(psi) E_R' M^-1(psi')^T` (re-derived; NF s.5).
7. **One check, not two** (NF s.5 review, NF s.6, A26): at genus 4-5 the tables give totals, so "the vertex graph carries the polar
   part" is the same check as "all edge graphs regular". Added. Letter: the hatted term is `G_Gamma` (was `G`, which in s.C.4 is the gauge
   factor of `cor:mai-natfact`).
8. **Genus-5 independent implementation** (l.76): "1474 of the 4555 regrouped graphs ... 19591 modulo `2^61-1`": the 19591 contributions
   are on 460 of these graphs (NF s.6); "regrouped graphs" -> graphs of the plain sum; the additivity identity spelled out; the
   common-mode caveat of NF s.6 (the second implementation shares the regrouping; the two primes test one implementation) added.
9. **Values pointed to from s05 and s10 were missing.** `rmk:red-leaf` (s05_reduction l.113) cites `sec:evi-transport` for
   `t_c - vhat = u^2/30 + O(u^3)`; C gave only `delta_P`. Added `t_c - vhat = delta_P/d`, `vhat = v/d`, with the frame-independence
   reason (`T_Pcan = t_c`, `prop:tra-can`; `d = tau^(1/4) varpi_0`, `lem:gau-diag`; DIGEST_C1_5 CV14). `rmk:mai-stored` (s10_main_4 l.61)
   says "the values of `m_5` and `m_10` ... are in sec:evi-univ"; C gave `m_5` elsewhere and `m_10` nowhere. Added to `sec:evi-univ`:
   `m_5 = -43397/807975000`, `m_10 = -765368951400317/1305647201250000000` (GAUGE_IDENTITY s.3; univ referee s.2.1: exact, both frames),
   `2 m_10 - m_5^2 = -18223115111903/15543419062500000` (GAUGE_IDENTITY s.4; re-checked to agree). The value of `m_5` was removed from
   the opening paragraph (pointer instead).
10. **Ranges stated more widely than the sources** (`c_evidence_2.tex` l.47, l.49): the `y_0(0) = 0` frames of the 17-frame check ran
    to `u^68` only (DB E1); the uniqueness window `#{k <= 15}` holds for the quintic pair, the random one-pole systems used
    `#{k <= 10}` (G2_6 s.10 vs G2R_3 s.4.2). Corrected; the two further cyclic implementations now carry their ranges (`u^70`, ten frames;
    `u^60`, four). `det Wr` -> `det Wr_f` of `def:gau-conc` (s03 uses `\Wr` for the matrix).
11. **Letter clashes in C**: `p_(5q,0)` and `z^(5q)`, `q <= 12` (parameter `q` as index) -> `n`; the LMI falsifier "a perturbation `eta`"
    (the pairing) -> words; "the delta lemma" named a lemma that the note does not contain (`lem:fsp-delta` was not written) -> "the
    formal delta lemma invoked in `rmk:gam-sevenvar` (not stated in `sec:fsp`)".
12. **Hodge paragraph** (`c_evidence_2.tex` l.6): "(tab:evi-g2)" was given for the coefficients of `s^-10..s^-7` per coloured graph;
    the table holds only the colour sums of `s^-10` (`[t_c^-2]`). Pointer made exact; "change of frame" -> "change of Hodge frame".
13. Smaller corrections: scope of the calibration (the `N_(1,d)` come from the genus-1 formula of the R-action, not from the graph
    sum of `def:set-graphsum`, which is stated for `g >= 2`; now said); LP20 vertex check "pins" -> "supports" (DL s.5 wording); degree-one calibration sentence made precise (FQ2 s.2:
    "the d = 1 genus split of the localisation matches the graph sum"); `tab:evi-g2` caption: frame `(0,0)` for both columns (FQ2 s.3:
    Cf at `alpha = beta = 0`); genus-4 bad-normalisation graphs named (bridges `(1,3)`, `(2,2)` at `u^4`; loop on `(3)` and three two-edge
    graphs at `u^5`, UNIVERSALITY l.102); LMI random cases ("twelve ... and a case with `phi'' = 0`", non-integer `M_*` flagged as outside
    `def:gam-lmidata`, four more mod-`p` cases, two with `a_* = sqrt(-1)`; LA2 s.4.2 P); twelve theories: `P^1`, `P^2` have two weight choices
    each (GEN s.4.1); the 17416 tuples are of legs in `L` and `-log p`; positive controls defined in the opening paragraph; "the two frames
    used throughout" -> "used most"; table columns ragged-right; "lemma Lemma", "theorem Theorem" duplications removed.

Checked and right as written: every other number of C (calibration values, predictions, `N_(3,0)`, `-5/144` (Faber-Pandharipande
constant, recomputed), instanton numbers; (M) at order `q`, `+600`; `kappa^LR` to `q^30` and the four controls; `2875`, `-3250`; all graph
counts and totals; the second-frame totals; content counts; caveats; `tab:evi-windows` (DD1 E1-E16, DD2 E1a-E7, V10); Theorem P,
degree, Theorem U, natural factorisation, normalisation controls, independence; the gauge, model, transport, universality and LMI
paragraphs not listed above).

## 3. Open points (none on the chain)

* O-rf8-1 (bibliographic, for the assembler and BIBPLAN): LP20's key table should read "s.2.4.3 Prop. 4, (17); s.2.4.4 (18)-(20)"; the
  appendix now cites (18). s02_setup_3 l.49 cites "[LP20, Section 2.4.3, (17), (19)-(20)]": (19)-(20) are in s.2.4.4 (see s.5 item 8).
* O-rf8-2 (outline item not carried over, minor): OUTLINE_9 lists "the response `Delta_g(5)`" for `sec:evi-model`; the sources give it
  only by pointer (DC1 EV-M3); not added.
* Hock25 is cited in `tab:sgn-tr` ((2.1)); BIBPLAN restricts it to `rmk:mod-hock` and `conv:mod-kernel`; the table repeats
  `conv:mod-kernel`'s dictionary, so this is harmless; the assembler may add "App. B" to BIBPLAN's "where".

## 4. Changes that affect other sections

* No label added, removed or renamed; no statement or definition used elsewhere changed.
* Pointers into C now satisfied: `rmk:red-leaf` (s05) and `rmk:mai-stored` (s10) need no change.
* `sec:sgn-norm` now writes the other normalisation as `P_can c(z)`, `c = 1 + m_5 z^5 + m_10 z^10 + ..` (as `rmk:mai-stored`), condition
  `p_(5n,0)(0) = 0`, `n >= 1`.

## 5. Cross-section list (read-only; for the assembler)

1. `\Wr` with two meanings: determinant in s02_setup_3 l.5 (`Wr(p_0..p_3) := det(..) = L^10`, as NOTATION) and **matrix** in
   s03_gauge l.37, l.80, l.84, s03_gauge_2 l.99 (`det Wr_x`).
2. `Lambda`: s06_model_3 l.3 "`Lambda_kappa := -5 kappa` (a number)" against NOTATION_2 s.7 ("`-5 kappa` ... never the leg map") and
   the leg map `Lambda` of s06_model_5 l.6.
3. `eta` used as a variable: s08_connection_6 l.5, l.9, l.33 (`lem:con-S2`: `W_1(eta)`, `SPn_(eta=1)[..; eta^nu_b]`) and a_formalsp l.42
   (`t = z^(1/2) eta`), against the pairing `eta`.
4. `s` as Mellin variable (`sigma = zs`): s08_connection_3 l.30, s09_gamma l.23, s09_gamma_4 l.93, against `s = u^(1/5)` (`def:set-points`).
5. `q` as an integer index: s08_connection_8 l.47-72 (`\mathsf K_q = [z^(5q)](..)`), in the same subsection as the parameter `q` (l.19).
6. `lambda`: torus weights `lambda_i` (s02); `lambda := log(-q)` (s08-s09, declared local); `lambda in Q^x` with `r = -1/lambda^2` in
   `cor:mai-univ` (s10_main_4 l.41, l.51; OUTLINE_2 s.0.4 item 11 allows it but NOTATION has no entry); scalars in s03_gauge_2 l.88,
   l.99; weights `lambda_0..lambda_(N-1)` in `rmk:gam-necessity` (s09_gamma_4 l.98).
7. `x`: `def:set-hat` (s02_setup_8 l.15) and s07_comparison_2 l.9 define `x := zeta^(-4i) Cexp(Delta_i^(1/2)) = s varpi_0`, while
   `x = W_v(p)` (s06) and `x` is the flat variable of `lem:con-form` (s08).
8. `\cn`: s03_gauge_4 l.97 writes free constants as `\cn_j` (`\cn_j exp(-Int chi^U_j)`), while NOTATION reserves `\cn_n` for the
   coefficients of `c(z)`.
9. `\cO`: the line bundle `O(5)` and the HPF operators `\cO_e` (s03_gauge_6 l.11) share the macro and glyph.
10. `D` (italic) for discs in `def:gam-lmidata` (s09_gamma l.65) and s09_gamma_3 l.37, while `D = q d/dq` (`\Dq`) everywhere else.
11. `q` with four meanings: the parameter; a point of the loop space (`def:set-graphsum`, `\overline F(q)`; `def:tra-tower`); the
    integration point of the TR kernel (`conv:mod-kernel`); `q_phi` (`def:mod-ancvec`). Typography does not separate them.
12. `G`: hatted term `G_Gamma` and gauge factor `G(u, zhat)` in the same remark (`rmk:mai-apriori`, s10_main_3 l.85-88).
13. `xi`: dummy of the Bernoulli generating functions in `def:fsp-stirling` (a_formalsp_5 l.10) vs `xi_a`, `xi_(a,d)` (s06).
14. `J`: J-function (s02), vector `J`, `J^C` (`conv:tra-forms`), `J^(j)_m` (s06_model_3 Step 2), `J_k` (LG dictionary, App. B).
15. `mu`: `mu_i` (s02, local), the slice scalar `mu` (s04), `mu(Y)` (`def:fsp-stirling`), `mu` in `lem:fsp-homog`, `mu` local in LMI step 2
    (s09_gamma_2 l.42), scalars in s03_gauge_2 l.88: all local, but `def:fsp-stirling`'s `mu(Y)` is not declared local.
16. LP20 pin in s02_setup_3 l.49: "Section 2.4.3, (17), (19)-(20)": (19)-(20) are in s.2.4.4 (arXiv v2 PDF, p.17).
17. Consistent (checked, no action): `d := (P_0)_33 = (P_0^vec)_00` (s04 vs App. B); `P^vec = eta^-1 P eta`; `Delta_i^(1/2) = zeta^(4i) C_0/L`;
    `<e>_1 = 1/8` pointer from s04_transport_3 l.70; the TR dictionary of `conv:mod-kernel` vs `tab:sgn-tr`; the dilaton equation of
    s06_model_4 l.62 vs `tab:sgn-tr`; `T_Pcan = t_c` (`prop:tra-can`) vs C; `delta_P` frame independence; `def:set-hat` vs `sec:sgn-var`.

## 6. Compilation

Scratch folder `allgenus_review_review_rf8` (copies only): pdflatex, bibtex, pdflatex x2: no error, no undefined reference or citation,
154 pages; the 15 overfull boxes are those of the baseline run (none in B or C); BibTeX: the known "missing pages in ABDKS25".
Files edited: `sections/b_signs.tex`, `b_signs_2.tex`, `c_evidence.tex`, `c_evidence_2.tex`, `c_evidence_3.tex` (85, 103, 80, 53, 58 lines).
