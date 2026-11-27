# Review rf3: section s03 (`sections/s03_gauge.tex` .. `s03_gauge_7.tex`), part 1 of 2

Task key `review/rf3`. Author-facing; not note text. Part 2 (`REVIEW_2.md`): open items, every change that affects other
sections, the compile record, the files changed.

## 0. Verdict: SOUND AFTER FIXES

* Every definition, statement and proof of s03 was checked line by line against OUTLINE_4 (s03), OUTLINE s.0 (rules) and s.0.2
  (chain), NOTATION (+ `_2`), LABELS, DECISIONS, the gap notes G2 (`GAPS_G2.md` .. `_6`, with G2R, G2R_2, G2R_3), G1_3 s.2.10,
  G1_4 s.3.1, GAPS_SYNTHESIS_2 s.3, the digests DB (S20-S38) and the s02 statements that s03 cites (`def:set-genus0`, `lem:set-pf5`,
  `prop:set-cyclicLR`, `cor:set-cyclicPF5`, `def:set-points`, `def:set-frob`, `def:set-rows`, `thm:set-qrr`, `lem:set-hrot`,
  `lem:set-frobbasis`, `def:set-frames`, `def:set-confexp`, `def:set-tau`, `lem:set-canon`, `def:set-model`, `prop:set-modrows`,
  `def:set-M0`). All hand computations were redone on paper (list in s.3).
* One statement-level gap on a step of the chain (C6), closed by correcting the statement, with no loss for the chain: F1 below
  (rationality of the actual colour-0 column was asserted in `thm:gau-R0` (iii) without proof; it needs (QRR)). No `\draftnote`
  was added; none is needed.
* Two scope errors (F2, F3): operator identities stated over `KK = Q((u))` but used in rings that do not contain `KK`
  (`Q[[q]][logq]`, the ring `\Nring`), and the converse case for `lem:set-pf5`. Fixed in s03.
* The rest are notation clashes, one remark that carried only a sketch (rule 10), and wording (F4-F12).
* No circularity: every `\cref` in a proof of s03 points to an earlier statement of s03 or s02 (checked one by one); forward
  references occur only in statements, remarks or definitions, as names (`def:gau-can` in `lem:gau-diag`, `lem:mai-Ghat` in
  `thm:gau-R0` (i), `thm:con-glue` in the two new remarks, `thm:flh-main`, `thm:gam-D` in `def:gau-Fg`, `lem:gau-vanq` (b) in
  `def:gau-van` (ii) for well-definedness).
* No finite computation is used by a proof; `sec:evi-gauge` is cited only in the introduction and in `rmk:gau-Uremarks`. No
  bookkeeping words, file names or tags in the `.tex` files. No citation in s03 (s03 uses no quoted literature directly; the
  inputs enter through s02, as the section's convention paragraph says).
* Labels: none changed, none added. All seven files stay below 250 lines (114, 130, 126, 100, 165, 89, 101).

## 1. Method

Each part was read in full; every formula was re-derived (s.3). The gap note G2 is the text of record for `sec:gau-U`, the
graded gauges, `lem:gau-diag` step 4, `thm:gau-G` (a) and the colours; its referee corrections C1-C5, R1-R3 were checked against the
text (all applied, except R3, which is now applied by F6). Theorem U was compared step by step with G2_2 s.3-s.4, G2_3 s.5.1-s.5.6
and G2_4 s.6-s.8; the colour lemmas with G2_5 s.9.1-s.9.4 and G2R_2 s.2.9-s.2.10, s.3 C2-C4.

## 2. Findings and fixes

**F1 (significant; chain step C6; `thm:gau-R0` (iii)-(iv), `cor:gau-rowsNN` (i), `cor:gau-colours`; s03_gauge_7.tex l.38-39,
53-58, 63, 84, 93, 96).**
* *Issue.* (iii) asserted `A^act_n in Q[[q]]`, hence (by `lem:gau-wkb` (ii) in its rational form) `D_0 in 1 + z Q[[z]]` and
  `A^act_n in Q[L, 1/L]`. But the rows are defined over `Q(zeta)` (`def:set-rows`: `R_k in End(H) (x) Q(zeta)[[q]]`), and the properties
  of `R` used in s03 (symplecticity, the QDE, the WKB form `cor:set-cyclicPF5` (b)) fix the colour-0 amplitude only up to a constant
  series in `1 + z Q(zeta)[[z]]` (`lem:gau-wkb` (ii) over `Q(zeta)`). Rationality needs (QRR) at `q = 0` (`A^act(q=0) = h(z)`, rational):
  then `D_0 = Phi^can(z,1)/h(z)` (the glue), or a Galois argument. The sources do not prove it: PS s.4.2 (b) and DB S20 take the
  "QRR-normalised" column to be rational, G2_5 s.9.2 "Fields" asserts it, G1_3 s.2.10 (c) takes `A^act_n in Q[L, 1/L]` from
  `thm:gau-R0` (iii). With the note's stated scope "colour 0 uses no (QRR)" the claim was unjustified.
* *Fix (a correction of the statement, not a `\draftnote`).* (iii), (iv) and their proofs now read over `Q(zeta)`:
  `A^act_n in Q(zeta)[[q]]`, `D_0 in 1 + z Q(zeta)[[z]]` (`lem:gau-wkb` (ii) with `\Bbbk = Q(zeta)`), `A^act_n in Q(zeta)[L, 1/L]`, the
  comparison in `Q(zeta)[L,1/L][[z]]` inside `Q(zeta)((s))[[z]]` (the injectivity argument is unchanged); the representatives of
  the actual rows lie in `\Nring (x) Q(zeta)` (as `cor:gau-rowsNN` already said). A remark after the theorem says why `D_0` is only
  known over `Q(zeta)` there and that its rationality, equivalently that of every `A^act_n`, follows from (QRR) at `q = 0` in the
  glue `thm:con-glue`, and is not used before. `cor:gau-colours` now states `D_i in 1 + z Q(zeta)[[z]]` with
  `[z^n] D_i = zeta^(-in) [z^n] D_0`; the clause "for `i != 0`, `D_i in Q[[z]]` iff `D_0 in Q[[z^5]]`" moved to the remark after the
  corollary, conditional on the rationality from `thm:con-glue`.
* *Effect on the chain.* None: (F_g), `thm:flh-main`, `thm:con-glue` and `thm:gam-D` hold verbatim over `Q(zeta)`, and the glue
  itself yields `D_0 = Phi^can(z,1)/h(z) in Q[[z]]`. The matching edit of the proof of `thm:con-glue` (s08) is listed in REVIEW_2 s.5.

**F2 (moderate; scope; `def:gau-conc` (iii), `lem:gau-adj`, `lem:gau-conc` (i); s03_gauge.tex l.32-36, 59-63, 73).**
* *Issue.* The adjoint, the concomitant and `lem:gau-adj`, `lem:gau-conc` (i) were stated for `P in KK[D]` and functions in
  differential rings extending `KK = Q((u))`. They are used in `Q[[q]][logq]` (`def:set-confexp`, s02: `D CC(C_0, I_1) = 0`) and in
  `\Nring` (`lem:gau-confexp`: the relation defining `D i_2`; `cor:gau-rowsNN`: the `logq`-coefficient `Sum c_a JJ_ab c_b = 0`), rings
  that do not contain `Q((u))`.
* *Fix.* `def:gau-conc` (iii) defines `P^*` and `CC_P` for coefficients in any differential ring and records that `PF4`, `PF4^*`, `CC`
  have coefficients in `Q[u]` (`D u = -(1-u)`) and act on every differential ring containing it (`KK[ell]`, `Q[[q]][logq]` with
  `u = 1 - 5^5 q`, `\Nring` with `u = L^-5`, `Q((s))`). `lem:gau-adj` is stated in `Q[u][D]` (its proof uses only
  `D o q = q o (D + 1)`, `q = (1-u)/5^5`). `lem:gau-conc` (i) holds for coefficients in any differential ring `\cR_0`. The proof of
  `lem:gau-confexp` now says that `CC` has its coefficients in `Q[u]`, `u = L^-5`. Proofs unchanged otherwise.

**F3 (moderate; scope; proof of `cor:gau-cycsys` (i); s03_gauge_2.tex l.125).** `lem:set-pf5` is stated "on any differential
`Q[[q]]`-algebra", but it is applied in a differential ring `\cR` extending `KK`, which is not a `Q[[q]]`-algebra (a power series in
`q = (1-u)/5^5` does not converge `u`-adically). *Fix:* a parenthesis in the proof: `lem:set-pf5` is an identity of operators with
coefficients in `Q[q][z]` whose proof uses only `D o q = q o (D + 1)`, which holds in `\cR`. Recommendation for s02: REVIEW_2 s.5 X4.

**F4 (minor; one notation, rule 6; `lem:gau-confexp`, s03_gauge_2.tex l.35-49).** The lemma used `I := c_0 logq + i_0` (s02 calls this
element `I_1`, `def:set-confexp`) and the local `lambda_0 := log(1-u)`, which clashes with the torus weight `lambda_0` (s02 has
`Xi = -log(1-u)`). Replaced by `I_1` and `-Xi`; every formula re-derived (`\Cexp(i_j) = D^j(varpi_1 + Xi varpi_0)`, `D Xi = -1`,
`\Cexp(I_1) = varpi_1`, `\Cexp(D i_2) = D^3(varpi_1 + Xi varpi_0)`). Added the half-sentence that the reduction to generators of
the localisation needs (`delta(x^-1) = -\Cexp(x)^-2 delta(x)` on `L`, `c_0`, `\mathsf W`).

**F5 (minor; local letters clashing with global symbols).**
* `lem:gau-conc` proof (i) (s03_gauge.tex l.89): `H`, `T_j` (hyperplane class, mirror map) -> `w`, sans `t_j`, declared local.
* `lem:gau-vanq` proof (s03_gauge_5.tex l.97): `beta_a`, `beta~_a` (frame parameter `beta`; the clash G2R_2 C5 removed for
  `beta = u Bprod`) -> `b^circ_a`, `b~^circ_a`, declared local.
* `lem:gau-HPF` proof (3) (s03_gauge_6.tex l.30): `e_j` (clashes with the index `e` of `\cO_e` in the same proof and with the
  idempotents) -> `sigma_j`, declared local, `sigma_0 = 1` added.
* `prop:gau-graded` (iii) (s03_gauge_4.tex l.86): summation index `q` (the parameter) -> `m`.
* `lem:gau-colmod` (s03_gauge_6.tex l.77-86): `\Fs^(i)` -> `F^(i)`, the letter of `prop:set-modrows` (b); "and `\Fs_b := \Fs^(0)_b`"
  (a redefinition) -> "thus `\Fs_b = F^(0)_b`".
* `rmk:gau-Uremarks`: the new local `E_j` (matrix units `E_ab`) written `\mathcal E_j`.

**F6 (minor; rule 10; `rmk:gau-Uremarks`, s03_gauge_4.tex l.97).** "The same arguments give five free constants at every order:
the solution module has rank 5" was the one-line sketch of G2R_2 s.3 R3, which the referee asked to label SKETCH unless the
five-constant version of `lem:gau-lead` is written. The remark now carries that argument: without the pole every `chi^U_j` is
integral, so each equation of `lem:gau-scal` (b) has one free constant `c_j`, the proof of `lem:gau-lead` gives the leading terms
`Theta_U Sum_j c_j E_j Y^j`, the construction of `thm:gau-U` (a) started from `E_j Y^j` gives five solutions, and the proof of
`thm:gau-U` (b) shows they form a basis. Still marked "not used".

**F7 (minor; scope of a statement; `lem:gau-diag`, s03_gauge_5.tex l.24).** "For the quintic system ... this gives
`p_(0,0) = tau^(-1/4)/varpi_0`" did not say for which gauge; it needs `p_(0,1)(0) = 1` (for `P_can c(z)` the value is
`c(0) tau^(-1/4)/varpi_0`). Now "and a regular gauge with `p_(0,1)(0) = 1` (for instance `P_can`)". All uses elsewhere are for
`P_can` (s04_transport_4 l.151, s07_comparison_2 l.56, 60, s08_connection_7 l.52, s10_main_3 l.74): unaffected.

**F8-F12 (wording).**
* F8 `lem:gau-uniq` (s03_gauge_5.tex l.80): "monic polynomial" -> "polynomial of degree 5 with leading coefficient `(n+1)^5 != 0`".
* F9 `thm:gau-G` (a) (l.110, 119): the dummy gauge of the uniqueness clause clashed with the given graded `P`; now `P'`, and the
  proof says "Conversely, let `P'` be any regular gauge ...".
* F10 paragraph after `cor:gau-cycsys` (s03_gauge_2.tex l.128): the gauges applied to the model columns were listed under
  `Q((z))((u))`; they live in `Q(zeta)((s))[[z]]` (`thm:gau-R0` (i)). Corrected; the case `w = 0` in `Q((z))((u))` is the one of
  solutions without exponential (the vanishing solutions).
* F11 `lem:gau-transv` proof (s03_gauge_2.tex l.13): "which by the same identity does not depend on `x`" -> the reason
  (under `x -> x M` the Gram matrix becomes `M^T CC_x M`).
* F12 `def:gau-conc` (ii): "with `log u` replaced by `ell`" (`lem:set-frobbasis` already uses `ell`) -> "of `lem:set-frobbasis` (iii),
  whose indeterminate `ell` is the one above"; after `lem:gau-colrows`: "the two covariance clauses of `def:set-rows`" (there are no
  such clauses in `def:set-rows`) -> "the colour covariance of the rows announced after `def:set-rows`".

Continued in `REVIEW_2.md` (s.3 checks found correct, s.4 open items, s.5 cross-section changes, s.6 compile record, s.7 files).
