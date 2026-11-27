# Gap closing, cluster G5, part 4: evidence, suggestions for the editor, dependencies, sources read

Continuation of `GAPS_G5.md`, `_2`, `_3`. Same rules. Everything in this part is author-facing.

## 5. Evidence (VERIFIED; exact; own code, Python 3.14 stdlib, `py -I -B`, single processes; logs print no paths)
Scripts and logs in this folder (`logs/`). No project module was imported or run; each script re-derives what it needs.

| script | log | what is checked | result |
|---|---|---|---|
| `g5_degree_bound.py 30` | `logs/g5_degree_bound_N30.log` | item 4, over `Q`, colour 0: (A) `P_0 = 0` on `x^m`, `-25 <= m <= 40`; (B) `P_j(x^m)` divisible by `1 - x^5`, `j = 1..5`; (C) `P_j(x Q[x]) ⊂ x Q[x]` (`1 <= m <= 40`); (D) `Pi(x^m)` in `x^6 Q[x][z]`; (T) `P_1 x^m = (1-m)(x^(m-1) - x^(m+4))`; (R) for `n <= 30` and two choices of free coefficients (0, and `frak f_m = (-1)^m (2m+3)/7`): `RHS_n = (1 - x^5) rho_n`, `rho_n` in `x Q[x]` (trivial resonance), `A_n` in `x Q[x]` with grade `4n+1 mod 5`, the equation `[z^(n+1)] P psi = 0`; (H) `A_n = x^(4n+1) ph_n(x^-5)` with `ph_n` from an independent HPF recursion in `u` (POLY s.1), `deg ph_n <= floor(4n/5)` | all PASS; `ph_5 = 47735919/128000000 - (203511/250000) u + (548547/1000000) u^2 - (6789/62500) u^3 + frak f_1 u^4`, as POLY s.1 prints |
| (controls) | same log | (X1) a lower amplitude with a term of valuation `<= 0` (`A_0 = x + x^-4`, same grade): `rho_n` leaves `x Q[x]` from `n = 1` and the resonance is OBSTRUCTED at `n = 5, 10`; (X2) twist factors `(1,2,3,4,10)`: the valuation property fails | both fail as predicted |
| `g5_flat_hatted.py 9` | `logs/g5_flat_hatted_N9.log` | item 1: the FLAT full-Stirling saddle series `c_n(L_0)`, computed by one-variable Wick contraction from the flat data at `sigma = L_0`, against `L_0^(4n) ph^SP_n(L_0^-5)` from HPF, `L_0 = 2, 3, -2, 1/2, 5/3` (`u = 1/32, 1/243, -1/32, 32, 243/3125`), `n <= 9`; `c_5(L_0) - L_0^20 ph^(0)_5(L_0^-5) = 4537/3937500 = B_6(1/5)/6` at all five points (`cor:gam-free`) | all PASS (also at `u = 32` and `u < 0`, outside `0 < u < 1`) |
| (controls) | same log | (X1) no half-density factor; (X2) wrong sign of the twist term; (X3) scaling `L_0^(4n-1)` | all fail at `n = 1` |
| `g5_natfact_checks.py 12` | `logs/g5_natfact_checks_K12.log` | item 3, over `Q(zeta)`: `Psi_0^T eta^-1 Psi_0 = 1`; OU3's `Psi eta^-1 Psi^T = eta^-1` FALSE (erratum 3.4); row 0 of `Psi_0` = `(zeta^(-4j))` (the hatted unit at `u = 0`); `mu(b,k)` integral, `>= max(0,(4k-4)/5)`, zero only at `k = 0` and `(1,1)` (`k <= 40`); the model T-system: closing relation, `Fs_b(0) = eps_b`, `Fs_0 = 1 - 3t/20 + 9t^2/800 + 11571 t^3/80000 + ..`; `Rhat^mod(0) = Psi_0`; `M_0(-t)^T M_0(t) = 1` through `t^12`; control: a perturbed `Fs_0` breaks symplecticity | all PASS; control fails |

Run times: below 3 s each. These are evidence only; no proof of s.1-s.4 uses them.

## 6. Suggestions for the editor (no design file was changed)
1. **`thm:con-Pa`**: adopt the statement of s.1.2 (parts (i)-(iv)) and the proof of s.1.3-s.1.5; delete its `\draftnote` (domain).
   `prop:gam-hatflat`: cite `thm:con-Pa` (iv) for the identity of rational functions; CP3 s.3.10 becomes an alternative. Add to the
   inputs of `thm:con-Pa`: `lem:fsp-univ`, `lem:fsp-ift`, `lem:fsp-cov`, `lem:fsp-homog`, `prop:gam-rational`. A one-paragraph formal
   lemma ("FL4 at the generic point", s.1.3 Step 1-2) may go to appendix A after `lem:fsp-dq` (the editor chooses the label).
2. **`cor:mai-gap`**: print the definition of s.2.2 (in `def:set-confexp` or at the corollary) and the proof of s.2.4; replace the
   `\draftnote`; cite [PPZ15] Def. 2.13 and Prop. 2.14 (M2). The frame scope is that of Theorem A (D14); the two remarks of s.2.5 (complex
   frames, scalings) are optional, proved, and the scaling remark explains "the same normalisation of the vanishing period".
3. **`rmk:mai-constant`**: s.2.6 replaces `\draftnote{compare def:set-confexp with CGYZ26 Defs. 1.1-1.2, Rem. 1.4}`; keep the version
   `\draftnote` (arXiv v2) and add "checked against the HTML rendering; to be re-read in the PDF".
4. **`cor:mai-natfact`**: define `Psi(u) := Qhat(u, 0)` and `R~ := Psi^-1 Qhat = Psi^T eta^-1 Qhat` (NOTATION table 8 already reserves the
   letters); print Lemma 3.A (explicit `Psi`, orthonormality) and the proof of s.3.3; delete the `\draftnote` on `Psi(u)`. Add
   `lem:gau-R1` and `lem:set-canon` to its inputs.
5. **`lem:gau-diag` step 4 (G8)**: the line is written in s.3.3, Remark (i) (it belongs to another cluster; offered, not owned).
6. **`prop:con-deg`**: adopt s.4.2-s.4.3; delete the `\draftnote` of step 5; step 4 no longer needs `lem:gau-poly`. The proposition stays
   outside the main chain, but it now also gives an independent proof of the "resonances for every choice" clause of `lem:gau-poly`.
7. **`design/SOURCE_ERRATA.md`** (for the main session): OU3 s.3, "`Psi eta^-1 Psi^T = eta^-1`" -> "`Psi^T eta^-1 Psi = 1`" (s.3.4);
   PB s.2.2, "On `0 < u < 1`" -> the comparison in `Q((s))` (s.1); ER s.4, step 5's sentence is correct but its "W" is undefined (s.4).
8. **R7's list** (OUTLINE_REVIEW): the entries `thm:con-Pa` step 5, `cor:mai-gap`, `cor:mai-natfact` (`Psi(u)`), `prop:con-deg` step 5 can be
   removed once the editor applies 1, 2, 4, 6; `cor:mai-gap` remains conditional on Theorem A and on R1 (s.7).

## 7. Dependencies on other clusters (what these resolutions assume)
* **R1 / D2 (cluster G1): `prop:set-cyclicLR` and `cor:gau-rowsNN`.** Used in s.2.2 (the rows, hence `F_g^NN`, lie in `NN`), in Lemma 3.A
  (the `z^0` rows via the forward relations) and in the scaling remark of s.2.5. If R1 stays a `\draftnote`, Corollary B and the
  explicit form of `Psi(u)` carry the same `\draftnote` (Theorem A already does, through `def:set-outer`).
* **Symplecticity of `M_0`** (`def:set-M0`; R7's list): used in s.2.4 Step 1 and s.3.3 (c). VERIFIED through `t^12` here (s.5).
* **`thm:gau-R0` (iv), `thm:gam-D`, `prop:gau-symp`, `lem:gau-R1`, `prop:gau-graded`, `thm:gau-U`** (clusters on s03, s09; D7, D15) as stated in
  the outline; `lem:flh-frame`, `lem:flh-comm` (s07) as written in S15 s.2.
* **(T) `thm:set-GT` and (QRR) `thm:set-qrr`** (quoted; the `\draftnote`s (a), (b) of `thm:set-qrr` belong to another cluster).
* **Appendix A**: `lem:fsp-univ`, `lem:fsp-ift`, `lem:fsp-dq`, `lem:fsp-transl` (b), `lem:fsp-stirfe`, `lem:fsp-cov`, `lem:fsp-homog`,
  `lem:fsp-param` (only for the holomorphic variant).
Nothing here uses topological recursion, Theorem P(b) or the finite route.

## 8. Sources read
* Design: `OUTLINE.md`, `OUTLINE_2.md`, `OUTLINE_3.md` (s02), `OUTLINE_4.md` (s03 parts), `OUTLINE_7.md`, `OUTLINE_8.md`, `NOTATION.md` (+`_2`),
  `LABELS.md` (+`_2`, searched), `design/OUTLINE_REVIEW.md` (+`_2` s.2-s.4), `design/DECISIONS.md`; digests DD1 (`_2`, `_4`, `_5` parts),
  DA (`_2`, `_3`, `_4`, `_7` parts), DB (`_2`-`_5` parts), DL (`_2`, `_4`, `_7` parts).
* Project sources: `allgenus/stage4/poly/POLY.md` (all); `allgenus/stage5/pb/PB_AUDIT.md` s.0-s.2.3; `allgenus/stage5/lemmas/FORMAL_SP_LEMMAS.md`,
  `_2.md` (all); `.../lemmas/referee/REFEREE.md` s.3.2-s.3.4; `allgenus/stage5/chart/CHART_PROOF_3.md` s.3.10-s.3.11;
  `allgenus/stage5/referees/pb_steps/REFEREE.md` s.4; `.../cyc_adversary/REFEREE.md` s.2.7-s.2.8, `_2.md` K5; `allgenus/stage4/evid/referee/REFEREE.md` s.4;
  `allgenus/stage2/flat/S15_PROOF.md` s.5; `natural_form/theory/UNIFORM_REDUCTION.md` s.0-s.4; `natural_form/referee/REFEREE_NAT.md` s.3(c)-(d)
  and the docstring of `ref_rows.py`; `natural_form/genus4_referee/review/REVIEW.md` s.2; `NATURAL_FORM.md` s.6; `OU_EXTENSION_genus3.md` s.3,
  s.6d (passages on `Psi`); `FORMAL_QUINTIC_genus2.md` s.1-s.5.
* Literature (public, arXiv): [PPZ15] arXiv:1303.1043 (current version v3, 2015-03-18), s.2: Def. 2.2, Props. 2.3, 2.4, Def. 2.5,
  Props. 2.7-2.9, 2.12, Def. 2.13, Prop. 2.14, read in the ar5iv HTML rendering of the current version through a text-extraction tool
  (the designer checked the same numbers in the v3 PDF, `OUTLINE_2.md` s.0.4 item 2). [CGYZ26] arXiv:2609.13921v2 (2026-10-31), s.1 and the
  passages listed in s.2.6, in the arXiv HTML rendering through a text-extraction tool. No personal data was sent.
