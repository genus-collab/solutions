# Review rf4, part 2 of 2: line-by-line record, effects on other sections, compile record

Continues `REVIEW.md` (verdict s.0, issues F1-F11 s.1, open items O1-O5 s.2).

## 1. What was checked, statement by statement

Read in full: OUTLINE.md (s.0, s.0.1, s.0.2), OUTLINE_5.md (s04, s05), DECISIONS.md, NOTATION.md/_2 (tables 6, s.9-s.10),
LABELS.md/_2 (s04, s05 rows, equation labels), BIBPLAN.md (PPZ15 row); gaps/G3/GAPS_G3.md-_4.md, G3/referee/REFEREE.md, _2.md
(K1-K10); GAPS_SYNTHESIS.md s.2.3, s.2.11, GAPS_SYNTHESIS_3.md s.1 (s04/s05 instructions); DIGEST_C1.md (D1-D19),
DIGEST_C1_2.md (A1-A6, B1-B7, C1-C6), DIGEST_C1_4.md (proof steps), DIGEST_C1_5.md (CV1-CV18), DIGEST_C1_6.md (O-2..O-9);
sources HVQ/allgenus/theory/REDUCTION_DRAFT.md (RD1 s.1-s.2), REDUCTION_DRAFT_2.md (RD2 s.3), theory/referee/REFEREE.md
(TRF s.0-s.2). The statements of other sections used by s04-s05 were read in the note itself: s02 (`def:set-frob`,
`lem:set-frobbasis`, `def:set-frames`, `def:set-tau`, `def:set-model`, `def:set-graphsum`, `def:set-modcohft`, `lem:set-modhom`,
`lem:set-dilaton`, `eq:set-cg`), s03 (`prop:gau-equalprod`, `thm:gau-U`, the definition of regular/graded gauges,
`lem:gau-R1`, `lem:gau-diag`, `def:gau-van`, `lem:gau-vanmod`, `lem:gau-vanshape`, `lem:gau-vanq`, `thm:gau-G`, `def:gau-can`,
`prop:gau-symp`), s06 (`def:mod-weights`, `cor:mod-basis`, `thm:mod-MLK` recall of ML), s07, s10 (every use of s04/s05 labels).

**1.1 s04 part 1.** `conv:tra-forms`: `eta = 5 J`, `eta^-1 = eta/25`; (Sp) for `P` iff for `P^vec` (re-derived:
`P^vec(-z)^T eta P^vec(z) = (1/25) eta P(-z)^T eta P(z) eta`, and `eta^3 = 25 eta`); lowered/vector forms; `4k + w_b =
5 wt(k,b) + 5`; antisymmetry and isotropy of `Omega_s`; the description of `Q((z))[[u]]` and of `HH^u_-`; summable families
and maps of finite type as G3 s.2.2. `lem:tra-sum`: proof as G3 s.2.2, correct. `def:tra-tower`, shift rule, values: the closed
form agrees with `lem:gau-vanshape`; `Pi_0(z^m H^b) = (-1)^(m+1) [z^-m] Y^C_b` re-derived; `nu_b + 1 = 1 - b mod 5`; (i)-(iii)
and the seven example values re-derived by hand. Homogeneity: `5(n - floor(n/5)) - 4n = nu_b + 1` because `0 <= nu_b + 1 <= 4`.
`eq:tra-pairing`: `Omega_s(z^-k J^C, q) = (-1)^k [z^(k-1)] <Y^C(-z), q(z)> = (-1)^(k+1) Pi_(k-1)(q)`, re-derived; it matches
the sources' `V^(k) = (-1)^(k+1) Pi_(k-1)` (DC1_5 CV1). `lem:tra-theta`: telescoping, nilpotency on finite supports, u-order
gain `>= 4`, weight bookkeeping; correct (G3_2 Step 1, G3R s.2.2). `def:tra-slice`: the equivalent description of `K_0` with
`z^-1 Q[z^-1][[u]]` (G3 s.3 item 6) is right (an element of `Q((z))[[u]]` without non-negative powers has polynomial
coefficients in `z^-1`).

**1.2 s04 part 2.** `lem:tra-frames` (a)-(d): proof = G3_2 Steps 1-5 with K1 (now stated precisely, F6). Support rule
`j = k - nu_a - 1 mod 5`, exponent `>= 1`. Trivial frames: exactly ten (`k < nu_a + 1`), proof re-checked (`j_0 + 5r >= k + 1`;
`f_(j_0) = eps_a v^(nu_a+1)/(nu_a+1)!`). The five examples re-derived from the recursion `f_i = -Pi_0(z^(k-i) H^a) - Sum_r
pi_r v^(4r) f_(i+5r)` (e.g. `frame(6,0)`: `f_5 = v`, `f_0 = v^5/480 - v^5/96 = -v^5/120`; `frame(2,0) = z^2 1 + v z H` in
`rmk:tra-consistency`). `cor:tra-slice` (`-z1` lies in `K`, so `-z + K = K`), `lem:tra-coeff` (`wt(k,a) - wt(j,1) =
(4(k-j) + nu_a + 1)/5`), `lem:tra-twoslot` (G3_2 s.3; the summability bound and both directions) correct.

**1.3 s04 part 3, `lem:tra-conv`.** (a) = `lem:tra-frames` (c) with the u-order of the frame coefficients. (b1): degree bound,
`wt >= wt(0,1) = -1`, `e_0 = min(0, (2-2h-n) ord_u v)` (K7: the printed bound uses homogeneity, and the dependency line says so).
(b2): the bound `e_0 + Sum max(0, (4(kappa - D) + 1)/5)`, iterated summation lemma, regrouping. (b3), (b4): as G3_3 s.2.
(c): as G3_3 s.2 with K2 (closed form `lem:gau-vanshape`, residue solution `lem:gau-vanmod`). (d): `wt(k,b) >= 3/5` for `k >= 2`,
`Sum k_i <= 3h - 3 + n - m`. (e): G3_3 s.4; (e2) needs `ord_u v = 1` for negative exponents (`v'(0) != 0`), the bound
`2 - 2g + j + 4l/5`; (e3) the difference `<H^j>_g - (F^m_g)^(j)` is homogeneous of negative exponent and in `Q[v]`, hence 0;
binomial resummation; the shorter argument's tail has u-order `>= g`. (f1): rising factorial `(N)_(m_1)`, `N = 2h-2+n+m_2 >= 1`,
every intermediate step stable; `Sum binom(N+m_1-1, m_1) c^(m_1) = (1-c)^-N = kappa^N`. (f2): the polynomial-in-`m` argument of
G3R_2 s.8 (K3), re-derived: `[u^o] T_m = c_0^m Sum_(m_2 in S) c_0^(-m_2) a_(m_2) binom(N_(m_2) + m - m_2 - 1, N_(m_2) - 1)`, a
non-zero polynomial since the degrees `N_(m_2) - 1` are distinct. The dependency line (K10) is printed. `rmk:tra-consistency`
(a)-(c): every identity re-derived (string shift `Pi_i([q/z]_+) = Pi_(i+1)(q)`; `<se>_h = 0` for `h >= 2`; `<se H^m>_1 = 0`);
it keeps its `\draftnote` (O1).

**1.4 s04 part 4.** `def:tra-dressing`: grading convention identical to `def:gau-can`; `(P_0^T eta P_0)_03 = 5 p_00 p_03`, so
`p_00 p_03 = 1`; `d = (P_0)_33 = (P_0^vec)_00`. Divisibility (G3_4 Step 0): `N(psi, -psi) = 0` from (Sp) inverted, the
substitution `psi' = xi - psi`, uniqueness, symmetry; `E_P` from `P^vec` because `P_0^vec` is an isometry. `eq:tra-leaf` (new,
F1): `P^vec = P_0^vec R^vec`, `(P_0^vec)^-1 1 = d^-1 1`. `lem:tra-E` (G3_4 s.2): the ring `Q((psi))[[u, psi']]`; `psi + psi'`
a unit with inverse `Sum (-psi')^j psi^(-j-1)`; Step 3 (no residue), Step 4 (`beta = R^-1(psi)^T alpha`), Step 5
(`R^-1(psi)^T = eta R(-psi) eta^-1`, `beta = eta (R f)(-psi)`, `beta_(>=0)(-psi') = eta (R f)_(>=0)(psi')`), the final
transpose `C = -R(psi')^-1 (R f)_(>=0)(psi')`: every sign and transpose re-derived; the side remark uses the K5 ring
`Q((psi'))((psi))[[u]]`. `cor:tra-edges`: `(R^vec z^-k J^C)_(>=0) = z^-k (P_0^vec)^-1 eta^-1 (Y')_(>=k)` with K6
(`P^vec J^C = eta^-1 (P Y^C)`); both directions of the equivalence. `lem:tra-T`: (a) from `eta (P^vec)^-1(z) = P^vec(-z)^T eta`;
(b) re-derived: `Pi_i(q_T) = (-1)^i Omega_s(z^(-i-1) eta^-1 Y', -z1) = (-1)^(i+1) [z^(i-1)] Y'_0`, using `(eta^-1 Y')^T eta 1 =
Y'_0`; it agrees with RD2's `V^(k)(q_T) = (-1)^(k+1) [z^(k-2)] Y'_0(-z)` under `V^(k) = (-1)^(k+1) Pi_(k-1)`; `i = 1` by the
grading (`Y'_0` has only `z^j`, `j = -1 mod 5`); (c) as F1, F2; (d) `T_P = -Sum_k p_(k,-k) [z^(-1-k)] Y^C_(-k)`, the `k = 0`
term `p_00 v`, the others of v-exponent `>= 2` (DC1 O-3 (c): the exclusion of the other `u^1` terms is written). `lem:tra-F`:
Step 3's inversion `[psi^k psi'^l] E = Sum_(j=0)^l (-1)^j [psi^(k+1+j) psi'^(l-j)] N` re-derived; Step 1 with `d(0) = 1` stated
(DEC Q5), `h = 1` needs `n >= 1` (DC1 O-3 (d)). `prop:tra-can`: every hypothesis traced to s02/s03 (`thm:gau-U` (d) rationality,
regular gauge = invertible `P_0`, grading `def:gau-can`, (Sp) `prop:gau-symp` with `cn = P_11(0,z) = 1`, (A) from
`thm:gau-G` (a) and `lem:gau-vanq` (c), `T = t_c` from `lem:gau-vanq` (d), `d(0) = 1` from `lem:gau-diag` under (N0)).

**PPZ15 numbering (F9).** Public arXiv HTML rendering of 1303.1043 (latest version): Def. 2.2 = the R-matrix action (vertices
`Omega`, legs `R^-1(psi)`, edges `(eta^-1 - R^-1(psi') eta^-1 R^-1(psi'')^t)/(psi' + psi'')`), followed by the sentence that the
inverse of R is used on the legs (Givental's and Teleman's convention); Def. 2.5 = translation action `T Omega`; Prop. 2.12 =
`T_a R Omega = R T_b Omega` with `T_b(z) = z(1 - R^-1(z)) 1`, a CohFT with unit; Def. 2.13 = `R.Omega := R T Omega` with
`T(z) = z1 - z R^-1(z) 1`; Prop. 2.14 = left group action. This matches BIBPLAN and every PPZ15 pin of s04-s05. (The HTML
rendering is of the latest arXiv version; BIBPLAN names v3, which I could not open as text; no discrepancy was visible.)

**1.5 s05.** `def:red-flatsum`: homogeneous vertex tensor (DEC Q5, G3 s.3 item 2), `|Aut|^-1` included (DESIGN_CHANGES_2
item 1); the vertex graph `eq:red-vertex`. `def:red-ML`: as DC1 D14 and as recalled in s06_model_10 l.5 (consistent); the
function-form equivalence holds (polarisation over `Q[v]`; for `h >= 2` the `n = 0` coefficient vanishes by the definition of
`F^m_h`; for `h = 0, 1` the excluded orders are those with `2h - 2 + n <= 0`). `lem:red-genus`: Betti-number count, a tree with an
edge has two leaves; the added sentence on unstable genus-0 vertices is correct (`2m > n + m - 3` for `n <= 2`). `thm:red-main`
Steps 1-10 re-checked against RD2 s.3.2, DC1_4 C1, G3_3 s.4: Step 5 needs `n >= 1` at every vertex (Step 2) so that one
insertion is a frame with `Pi_0 = 0`; Step 6 is `lem:tra-conv` (b4) with the hypotheses checked (edges with property (i) by
`cor:tra-edges`, leaves in `z^2 HH^u_+ cap K` by `lem:tra-T` (c)); Steps 8-9 are `lem:tra-conv` (e) with `qt_T` in
`-z1 + z^2 HH^u_+`, `Pi_i(qt_T) = 0`, `delta_P` in `u^2 Q[[u]]`; Step 10 `d^(2g-2) c_g (d T_P)^(2-2g) = c_g T_P^(2-2g)`. ML is used
only for `h <= g` (MSG2's "every `h`" corrected), the gauge equation is not used, and no claim that each hypothesis is needed
is made. `cor:red-can`, `rmk:red-hyp` (C2 (i)-(v)), `rmk:red-leaf` (F5) correct. Binomial expansion in `rmk:red-leaf`
re-derived (`-(g-1) c_g 2 vhat^(1-2g) (t_c - vhat)`, error `O(u^(4-2g))`); the genus-2 numbers of appendix C are consistent
with it (`-7/6000 + 1/900 = -1/18000 = -(g-1) c_2/15`).

**Binding rules.** No project bookkeeping in the five files (grep for stages, referees, statuses, file or code names, "route",
"tested", "checked", "verified": none left); every citation pinned (only [PPZ15] is cited, Def. 2.2, Prop. 2.12, Def. 2.13);
no proof cites a computation; every `\cref` inside a proof points to s02, s03 or an earlier item of s04-s05 (scanned
mechanically after the edits: no target in s06-s10 or the appendices); topological recursion is not used; the conditionality
(the PPZ15 definitions as part of (Q); the dilaton lemma through [AC96, (1.7)]) is stated in the s04 and s05 introductions.

## 2. Gap-note corrections: where they stand in the edited files

K1 `lem:tra-frames` (c) graded clause (F6) - done. K2 `lem:tra-conv` (c) cites `lem:gau-vanshape` and `lem:gau-vanmod` - done.
K3 (f2) polynomial-in-`m` argument - done. K4 homogeneous vertex tensor in `def:red-flatsum` - done. K5 side remark ring - done.
K6 `cor:tra-edges` line `P^vec J^C = eta^-1 (P Y^C)` - done. K7 `e_0` via homogeneity, dependency printed - done. K8 (b) the
`Theta` lemma `lem:tra-theta` - done. K10 dependency line of `lem:tra-conv` - done. DC1 O-3 (a)-(g): all written (frames for
infinite sums, duality step `lem:tra-twoslot`, exclusion of the other `u^1` terms in `lem:tra-T` (d), `h = 1` in `lem:tra-F`,
genus bound `lem:red-genus`, u-adic convergence `lem:tra-conv` (b), resummation `lem:tra-conv` (e)). O-6 ("monomials in
`v(u)`") - done. O-7 (vector form fixed once) - done. O-8 (ring of Lemma E) - done. O-9 (`rmk:red-leaf`) - now proved (F5).

## 3. Changes that affect other sections or the design files

* `eq:tra-leaf` (new auxiliary display label, s04_transport_4): defines `T_(>=2)` for every dressing with (Sp). Used in s04 and
  s05 only. No other section uses `T_(>=2)` (grep).
* `lem:tra-T` (c) reworded (F1, F2); same content except the corrected `y°` clause. External uses: s10_main_3 l.98 (cites
  `lem:tra-T` (d) and `lem:tra-T`), c_evidence_3 l.32 (cites `lem:tra-T`): unaffected.
* `lem:tra-frames` (c) reworded (F6). External uses: s06_model_6 l.59 (part (a)), l.77 (`cor:mod-basis`, graded clauses (b),
  (c)): consistent with the new wording.
* `lem:tra-conv` (c) reworded (F8). External uses are of part (f) only (s02_setup_7 l.72, s07_comparison l.24,
  s07_comparison_2 l.93): unaffected.
* `def:tra-edge`: the sentence defining `Omega^P` replaced (F3); `Omega^P` occurs nowhere else. `E_R`, `E_P`, `q_T`, `qt_T` and
  the divisibility statement are unchanged; `def:tra-edge` still cites [PPZ15, Def. 2.2], and [PPZ15, Prop. 2.12, Def. 2.13]
  now sit in the paragraph after it (same subsection `sec:tra-dress`). External uses: s01_intro_3 l.17 and s10_main l.33
  ("[PPZ15] enters in `def:tra-edge`": still true; the owners may write `sec:tra-dress`), s07_comparison l.38,
  s07_comparison_2 l.50, l.93 (edge and leaf point `q_T`: unchanged), s10_main_4 l.74 (the edge is a power series only under
  (Sp): unchanged). The phrase "the flat sum of `P.Omega^mod_(v(u))`" used in s01, s05, s07, s10 is now defined in `def:tra-edge`.
* `rmk:red-leaf` (F5): the `\draftnote` is gone and the `u^2/30` is derived. Design files that still list it as computed or as
  carrying a `\draftnote`: OUTLINE.md s.0 rule 3 (the example "`t_c - vhat = u^2/30 + ..`" of a number known only by
  computation), OUTLINE.md s.0.2 list O2 (`rmk:red-leaf` among the remarks with a `\draftnote`), OUTLINE_5.md (`lem:tra-T` step 7
  "[V]: appendix C only"; `rmk:red-leaf` `\draftnote`), GAPS_SYNTHESIS_2.md s.4 (planned `\draftnote` of `rmk:red-leaf`, O-9).
  For the editor; I did not edit them. s10_main_3 l.98 already derives `d = tau^(1/4) varpi_0`, `vhat = u tau^(-3/2)/varpi_0`
  and points to `rmk:red-leaf`: consistent. Appendix C (c_evidence_3 l.32) keeps the computed series as evidence: consistent.
* Unit notation (F10): s04-s05 now write `1` as s02, s07 and appendix B do. No macro was added or changed.

## 4. Compile record; files

Scratch folder `<user home>/AppData/Local/Temp/claude/d7c1b608d03d194a/37f24460-3601-4c85-b123-5156e19eb475/scratchpad/
allgenus_review_review_rf4/` (copies of main.tex, preamble.tex, macros.tex, refs.bib, sections/): baseline pdflatex, bibtex,
pdflatex, pdflatex: no error, no undefined reference; after the edits pdflatex twice (three times in all): no error, no
undefined reference, no LaTeX warning, 153 pages. Edited: `NOTE/sections/s04_transport.tex` (108 lines), `_2` (84), `_3` (73),
`_4` (152), `NOTE/sections/s05_reduction.tex` (119). Written: this file and `REVIEW.md`. Nothing else in NOTE was touched.
Evidence-only script (not used by any statement): `rf4_scripts/delta_check.py` in the scratchpad (exact series of
`tau^(1/4) varpi_1 - u tau^(-5/4)`: `u^2/30 + (127/3150) u^3`).
