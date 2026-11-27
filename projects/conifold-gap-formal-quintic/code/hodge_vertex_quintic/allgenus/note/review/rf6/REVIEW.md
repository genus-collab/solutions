# Review rf6: s07 (`s07_comparison.tex`, `_2`) and s08 (`s08_connection.tex` .. `_8`), part 1 of 2

Task key `review/rf6`. Author-facing; not note text. Part 2 (`REVIEW_2.md`): what was checked and found correct, open items,
every change that affects other sections, the compile record, the files changed.

## 0. Verdict: SOUND AFTER FIXES (both sections)

* s07 (C10, the comparison theorem) and s08 (C11, Theorem P and the glue) were read line by line against OUTLINE s.0 (rules), s.0.2
  (chain), OUTLINE_7, NOTATION (+`_2`), LABELS (+`_2`), BIBPLAN, DECISIONS, the gap note G5 (`GAPS_G5.md` s.1, `GAPS_G5_3.md` s.4) with
  `G5/referee/REFEREE.md` (F1.1-F1.4, F2.x) and `REFEREE_2.md` (F4.1-F4.3), GAPS_SYNTHESIS_2 s.3, the digest DC2 (`DIGEST_C2_3.md`
  s.2.4, `_5`), the source referee `stage2/flat/referee/REFEREE.md` s.1, and every statement of s02-s05, s09, s10, appendix A that the
  two sections cite. Every formula was re-derived by hand (list in part 2, s.1).
* One gap on a step of the chain (C10), the known issue: the unconditional symplecticity of the outer matrix `R'` in
  `lem:flh-inputs` (iii). It is closed without a `\draftnote` and without circularity (F1). A second instance of the same transfer
  problem, in `lem:flh-inputs` (iv) (the unit and isometry of `Psi_Q`), is closed the same way (F2). In s08 the a priori identity
  `D_0(z) D_0(-z) = 1` of `rmk:con-genus` rested on the same transfer; the remark is corrected (F8). None of these enters the chain
  any more.
* One gap on chain step C11 created by a concurrent correction of s03 (rf3 F1: `thm:gau-R0` (iii) now gives `D_0` only over
  `Q(zeta)`): the proof of the glue `thm:con-glue` assumed rationality. Rewritten over `Q(zeta)`; the glue now proves that `D_0` and
  the actual amplitudes are rational (F7).
* The rest: one scope slip in Theorem P(a) (F9), small statement slips (F3-F6, F10, F11), notation and typography (F12).
* No `\draftnote` was added or removed; none is needed. No label was added, renamed or removed. No finite computation is used by a
  proof (appendix C is cited only in remarks, as "not used"). No bookkeeping words, file or code names in the `.tex` files.
* No circularity: every `\cref` in a proof of s07/s08 points to s02-s06, to an earlier statement of the same section, or to appendix A
  (checked one by one). Forward references occur only in statements, remarks and comments (`thm:gam-D`, `cor:mai-*`, `rmk:mai-*`,
  `thm:gam-ginf`, `cor:gam-equiv`, `sec:evi-*`). The proof of `thm:gam-D` (s09) uses `thm:con-glue`, `thm:con-Pb`,
  `prop:gam-hatflat`, `cor:gau-colours` and nothing of s07 (grep: no `flh` label and no `R'` in s09).

## 1. The known issue: symplecticity of `R'` (settled consistently with rf2)

**F1 (significant; chain C10; `lem:flh-inputs` (iii), `s07_comparison_2.tex` l.36-54, 63; uses in l.38 of part 1, l.80, 96, 106).**
* *Issue.* Old (iii): "Without (F_g), `R'` is symplectic", proved from `def:set-rows` (a) (`Q(-z)^T eta^-1 Q(z) = 1`) and
  `M(-z)^T M(z) = 1`. But `def:set-rows` (a) is an identity of power series in `q`; `Q` in s07 is the conifold expansion of the
  representatives of the rows in `Nring (x) Q(zeta)` (`cor:gau-rowsNN`). An identity between evaluations under `ev` is an identity in
  the ring only where `ev` is injective, which is not claimed (OUTLINE O4). The same unjustified transfer is in the sources: FR s.1
  ("holds for the ACTUAL `R'` with or without (F), since `R~` and `M_nat` are symplectic"), DC2 `rmk:flh-fail` (5), and OUTLINE_7
  (`thm:flh-main` bullet, "[z^10] D_0 = ([z^5] D_0)^2/2 is forced by `prop:gau-symp`"). G5R s.2.3 already obtains the symplecticity of
  the expanded rows only from `D = 1`.
* *Fix (a corrected statement with a proof, not a `\draftnote`).* New (iii): `P_can` and `P_can Psi^m` are symplectic, and
  `R'(-z)^T eta^-1 R'(z) = M(z) (D(z) D(-z))^-1 M(z)^-1`; hence `R'` is symplectic iff `D_0(z) D_0(-z) = 1`, and under (F_g) it is
  symplectic modulo `z^(3g-2)`. Proof: `R' = S K`, `S = P_can Psi^m` symplectic (by `prop:gau-symp`, `eta^-1 = eta/25`, (i)),
  `K = M D^-1 M^-1` (proof of (ii)), `M(-z)^T = M(z)^-1`, `D` diagonal, `D_i(z) = D_0(zeta^-i z)`. A paragraph after the lemma says why the
  `q`-series identity is not used, that the injectivity of `ev` is neither claimed nor needed, and that a posteriori (`thm:gam-D`)
  `R' = P_can Psi^m` is symplectic; `thm:flh-main` does not use this.
* *The proof of `thm:flh-main` now uses only the window.* Step 1 states the coefficient count for the edge of `def:set-outer` as now
  defined by rf2 (symmetrised exact-division numerator `N - N(-psi',psi')/2 - N(psi,-psi)/2`): the coefficient of `psi^a psi'^b`
  involves `N` only in total degree `a + b + 1`. Step 2: `R'' := P_can Psi^m` is symplectic, so its correction terms vanish and
  `E_(R'')` is the edge of [PPZ15, Def. 2.2]; `R' = R'' mod z^(3g-2)` gives equality of the numerators in total degree `<= 3g-3`,
  hence of the edge coefficients of total degree `<= 3g-4` (the only ones that enter) and of the leaves mod `z^(3g-1)`.
* *Consistency with rf2.* rf2 (S1 of `review/rf2/REVIEW.md`; `def:set-outer`, `s02_setup_8.tex`) settled the same point the same way:
  symmetric division-free edge as the definition (valid for every `R'` with `R'(0)` invertible), the same defect formula
  `M(z) diag((D_i(z) D_i(-z))^-1) M(-z)^T` (equal to mine, as `M(-z)^T = M(z)^-1`), symplecticity modulo `z^(3g-2)` under (F_g) and
  exactly by `thm:gam-D`, no circularity. s10 (rf1, S10-2, S10-19) uses the same edge and cites the case `D = 1` of (iii).
* *Other uses fixed in s07:* the paragraph after `def:flh-outer` (part 1, l.38: "Whether `R'` is symplectic is not decided by the
  definitions ..."); the intro (l.4: symplecticity of `R'` used only modulo `z^(3g-2)`); the paragraph after `thm:flh-main` (l.96,
  F4); `rmk:flh-fail` (v) (l.106: equivalent to `D_0(z)D_0(-z) = 1`, window under (F_g), exact by `thm:gam-D`; not inferred from large
  radius).

**F2 (moderate; chain C10; `lem:flh-inputs` (iv), `s07_comparison_2.tex` l.56-60).**
* *Issue.* "As `R(0) = 1`, `Q^(i)_a(0) = eta(H^a, e~_i)`, so `Psi_Q` maps `e~_i` to `e~_i` in `H`, hence `1_TFT` to `1`": a statement about
  `q`-series, transferred to the conifold expansion in the same way as in F1. The isometry of `Psi_Q` was not argued.
* *Fix.* The constant terms are computed from the representatives of `cor:gau-rowsNN` (i) at `n = 0`: `(Psi_Q)_(aj) = zeta^(j(a+1)) m_a`,
  `m = (L/C_0, L^2/(C_0 C_1), C_0 C_1/L^2, C_0/L, 1)`, starting from `A^act_0 c_0^-1 = L/C_0`; with `Delta_j^(-1/2) = zeta^j L/C_0`,
  `Sum_j Delta_j^(-1/2) (Psi_Q)_(aj) = 5 delta_(a3) = eta(H^a, H^0)` in `Nring (x) Q(zeta)`, hence after expansion: (unit) with `varkappa = 1`.
  The isometry of `Psi_Q` follows from `Psi_Q = P_0 Psi^m` (z^0 of (ii)) with `P_0^T eta P_0 = eta` (`prop:gau-symp` at `z^0`). This
  avoids the forward reference to `lem:mai-Psi` (s10, now placed in `s10_main.tex`), whose part (i) is the same computation.
  The agreement remark `d = tau^(1/4) varpi_0` now names `p_00 p_03 = 1`.

## 2. Further findings and fixes in s07

**F3 (minor; accuracy; intro l.4, `rmk:flh-dilfactor` l.115).** "The only property of `P_can` it uses is its symplecticity": the proof
also uses that `P_0` is diagonal (graded) and `d(0) = 1` (`lem:gau-diag`), through (iv) and `lem:flh-dil`. Both places now say so.

**F4 (minor; the paragraph after `thm:flh-main`, l.96).** "`[z^10] D_0 = ([z^5] D_0)^2/2` is forced by `D_0(z) D_0(-z) = 1`" presented the
identity as known. Now conditional ("if moreover `D_0(z)D_0(-z) = 1`, which by (iii) is the symplecticity of `R'`"); the genus at which
`z^10` first enters is named (genus 5, `10 <= 3g-3`); both conditions hold a posteriori by `thm:gam-D` and are not used.

**F5 (minor; `rmk:flh-fail`, l.102-107).** (i) "If `D` is a scalar series `c(z)1`, the proof identifies `T_Gamma` with the flat term of
`P_can c(z)^-1`": the flat sum is defined for dressings with (Sp), so `c(z)c(-z) = 1` is needed (`prop:gau-symp`); added, and "scalar"
explained ("`D_0` a series in `z^5`", `cor:gau-colours`). (ii) The regrouping in `cor:mai-gap` uses `D = 1`; said. (vi) The vertex
CohFT is that of the model up to `Psi^m` *and the rescaling by `d`* ((i), (iv)).

**F6 (minor; scope; `lem:flh-dil` proof, part 1 l.66).** `lem:set-dilaton` (iii) was applied over `Q(zeta)((s))`, not a complete local
ring. rf2 (S8) extended the hypothesis of (iii) to `Q(zeta)((s))` with `varkappa in 1 + u Q[[u]]`; the sentence now cites that case.

## 3. Findings and fixes in s08

**F7 (significant after rf3 F1; chain C11; `thm:con-glue`, `s08_connection_8.tex` l.8-20, and `cor:con-A2` proof l.28).**
* *Issue.* rf3 corrected `thm:gau-R0` (iii)-(iv) and `cor:gau-colours`: without (QRR), `D_0` is known only in `1 + z Q(zeta)[[z]]` and
  `A^act_n` only in `Q(zeta)[L, 1/L]` (the rows are defined over `Q(zeta)`; `lem:gau-wkb` (ii) over `Q(zeta)`). The glue's proof still
  read (b) "In `Q[[q]]`", (c) `A^act_n in Q[L,1/L]`, (d) `D_0 in Q[[z]]`, (e) evaluation `Q[L,1/L] -> Q` (rf3 REVIEW_2 X1, O1).
* *Fix.* Statement: "In `Q(zeta)[[z]]`, `h(z) D_0(z) = Phi^can(z,1)`; hence `D_0 in 1 + z Q[[z]]`, and, by `thm:con-Pb`,
  `D_0 = Phi^SP(z,1)/h(z)`." Proof (a)-(e) over `Q(zeta)` (WKB uniqueness with `\Bbbk = Q(zeta)`, embedding
  `Q(zeta)[L,1/L][[z]] -> Q(zeta)((s))[[zhat]]`, evaluation `Q(zeta)[L,1/L] -> Q(zeta)`); the value at `L = 1` of the Laurent polynomial
  `A^act_n` is its `q^0` coefficient because `(1 - 5^5 q)^(-1/5) = 1 + O(q)` (now said); since `h` and `Phi^can(z,1)` are rational and
  `h(0) = 1`, `D_0` is rational, and `A^act = A^can D_0^-1` has coefficients in `Q[L, 1/L]`. The displayed equation `eq:con-glue` now
  carries only `h D_0 = Phi^can(z,1)`; its two users (`cor:con-A2`, `thm:con-conn` (iii)) were checked, and `cor:con-A2` now cites
  `thm:con-Pb` as well. (QRR) still enters exactly once, here.

**F8 (moderate; `rmk:con-genus`, `s08_connection.tex` l.28-31; `rmk:con-why` l.35 of part 8).**
* *Issue.* "Symplecticity halves the conditions: `D_0(z) D_0(-z) = 1`. Indeed, the rows satisfy `Q(-z)^T eta^-1 Q(z) = 1`, because `R` is
  symplectic": the transfer of F1 again (source: PR C2, "from `prop:gau-symp` for the actual and the canonical gauge"; the actual gauge
  is symplectic only if the expanded rows are).
* *Fix.* The remark now proves the equivalence `Q(-z)^T eta^-1 Q(z) = (D(-z)D(z))^-1` (colour by colour, from `cor:gau-colours`,
  `prop:gau-symp`, the model columns), i.e. symplecticity of the expanded rows (of `R'`, `lem:flh-inputs` (iii)) iff
  `D_0(z)D_0(-z) = 1`; says that carrying the `q`-series identity to the expansion would need it between the representatives in
  `Nring (x) Q(zeta)`, which is not claimed; so the halving is not available a priori and is not used (`D_0 = 1` is proved directly,
  `thm:gam-D`). The genus bookkeeping (`z^n` enters (F_g) first at `g = ceil((n+3)/3)`) is kept. `rmk:con-why` now says that
  symplecticity "could at most give `D_0(z)D_0(-z) = 1`, and only through the symplecticity of the expanded rows".
* Its last sentence ("only the symplecticity of `G` detects `D = 1`, known by computation") contradicted the remark's own invariance
  statement (symplecticity is blind to `exp(a)`, `a` odd; also rf1 S10-4). Replaced by the proved statement: a scalar change of the rows by
  `c(z)^-1`, `c in 1 + z^5 Q[[z^5]]`, multiplies `G` by `c(z)^-1`, which keeps `G(u,0) = G(0,zhat) = 1` and polynomiality and violates only
  the symplecticity of `G`, exactly when `c(z)c(-z) != 1`; the computations only illustrate this.

Continued in `REVIEW_2.md` (F9-F12; checks; open items; cross-section changes; compile record; files).
