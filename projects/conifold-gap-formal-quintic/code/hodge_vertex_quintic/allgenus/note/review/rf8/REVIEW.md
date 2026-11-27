# Review rf8: appendix B (`b_signs.tex`, `b_signs_2.tex`) and appendix C (`c_evidence.tex`, `_2`, `_3`)

Reviewer's report for the assembler (part 1 of 2; part 2: `REVIEW_2.md`, the read-only cross-section list). Not part of the note.

## 0. Verdict

| section | verdict | why |
|---|---|---|
| appendix B (`sec:sgn`) | **sound-after-fixes** | every dictionary entry re-derived or checked against the body (s02-s10, app. A) and the pinned literature; one wrong statement-number pin (LP20), two notation clashes, one false sentence on the constant term of `\SPn`, several imprecisions; all fixed |
| appendix C (`sec:evi`) | **sound-after-fixes** | every number and range compared with the sources (DIGEST_A/B/C1/C2/D1/D2/L s.5, the gap notes' evidence sections, NF s.5-s.6); one normalisation slip (the `(-5)` power of the TR/model comparison), several misattributions and overstatements (mutation tests, cancellation counts, frame windows, genus-5 controls), missing values that s05/s10 point to; all fixed |

Neither appendix carries a step of the chain of Theorem A. No `\draftnote` was added; no step of the chain is affected. The
appendices prove nothing, and no proof of s02-s10 cites appendix C: each of the 26 references to `sec:evi*` from s01-s10 lies in a
remark, in the introduction, in an explanatory paragraph after a proof (s10_main_3 l.47) or in s08's `sec:con-finite`, which the main
proof does not use (checked by grep with the enclosing environment).

## 1. How the check was done

* Appendix B: each convention compared with the place where the body defines or uses it (`def:set-rows`, `def:set-graphsum`,
  `def:set-frob`, `def:set-points`, `def:set-hat`, `lem:flh-hat`, `thm:set-qrr`, `rmk:set-qrrsign`, `lem:set-hrot`, `conv:tra-forms`,
  `def:tra-dressing`, `lem:tra-E`, `lem:tra-T`, `def:tra-tower`, `def:gau-van`, `def:gau-onepole`, `def:gau-can`, `lem:gau-diag`,
  `lem:set-modhom`, `conv:mod-kernel`, `def:mod-curve`, `thm:mod-tr`, `prop:mod-Rtr`, `cor:mod-trF`, `def:mod-legs`, `prop:mod-LD`,
  `def:mod-pair`, `thm:mod-var`, `prop:set-modrows`, `def:set-M0`, `lem:con-S2`, `def:fsp-wick`, `def:fsp-stirling`, `lem:fsp-stirling`,
  `lem:fsp-cov`, `lem:fsp-homog`, `def:gam-c`, `def:gam-lmidata`, `rmk:gam-sevenvar`, `rmk:mai-stored`, `cor:mai-univ`). Every
  formula was re-derived by hand: the colour-basis edge and leaves from `def:set-graphsum` (including `Delta_i^(1/2) Q^(i)_0(0) = 1`),
  `P^vec = eta^-1 P eta`, `h_1 = -3/20` and the `+5` value `-11/60`, the branch points and critical values of `W_v`, the
  `4^(2h-2+n)` scaling, the LG dictionary `J_(nu_b)(-4 zeta^-i t) = eps_b zeta^-(i(b+1)) Rhat^m` (via `lem:fsp-homog` and `lem:con-S2`),
  the leg-dictionary consequences (unit `p`, dilaton `-phi^0_1`, `H -> -log p`, `e -> -p^5/5`, the second dictionary's `-(4/5)p^5`,
  `-(16/45)p^9 + (V/5)p^4`), `Pi_0(H) = 1`, `Pi_0(z1) = -v`, the support rule `b + k = 1 mod 5`, `<e>_1 = 5/24 - 1/12 = 1/8`, the EO
  sign dictionary (kernel, `omega_(0,3)`, `F_h`, dilaton equation), the critical equation `u sigma^5 = 1`, `Q = -q` resp. `(-q)^(-1/5)`.
* Literature, public sources (no personal data sent): **arXiv math-ph/0702045v4 PDF** (EO): (4-3) is `dE_q(p) = (1/2) Int_q^(qbar) B`
  with footnote 3 (opposite to earlier papers), (4-9), (4-10), (4-20) and Thm 4.7 (4-30) `= (2g + k - 2) W_k + delta..`, "for `k = 0`,
  `g >= 2` it holds by definition if `W_0 = -F`": the table `tab:sgn-tr` is right as written. **arXiv 1803.01409v2 PDF** (LP20):
  (10) is "the T-weight on the fiber over `p_i` of `O(5)` is `5 lambda_i`" (pin right); **s.2.4.3 contains only Prop. 4 and (17); (18),
  (19), (20) are in s.2.4.4** (see 2.1 below).
* Appendix C: every number, range, field, prime and control compared with DIGEST_A_5-_6 (V1-V13), DIGEST_B_5 (E1-E12), DIGEST_C1_5
  (EV-*), DIGEST_C2_6 (E1-E20), DIGEST_D1_5 (E1-E16, controls, caveats), DIGEST_D2_6 (E1a-E14), DIGEST_L_7 s.5, the gap notes
  (G1_3 s.2.12, G1R s.1.2, s.2.2, G1R_2 s.6, s.9; G2_6 s.10, G2R_3 s.4; G3_4 s.4, G3R_2 s.12; G4_2 s.8, G4_3 s.4, G4_4 s.5, G4_5 s.5,
  G4R s.1.4, G4R_2 s.3; G5_4 s.5, G5R s.1.4, G5R_2 s.3.3, s.5; GS_3 s.8) and the project notes they cite where a digest was
  ambiguous (NATURAL_FORM s.5-s.6, FORMAL_QUINTIC_genus2 s.1-s.3, GAUGE_IDENTITY s.3-s.4, univ/UNIVERSALITY, univ/referee,
  stage5 lmi_adversary REFEREE_2 s.4, stage5 general GENERAL s.4). Arithmetic re-checked (Python, exact): the genus-2 column sums
  `1/1200`, `0`; `c_2..c_7`; `F^TR_2..F^TR_5`; the graph counts `42, 379, 4555`; `229, 2891, 1663`; the leading coefficients
  `-7/6000, -1/9000, -7/300000, -7/825000` and the second-frame ones `-11/189000, -11/900000, -1/225000`; the response sums `0, 0`;
  `2 m_10 - m_5^2`; `delta_P = u^2/30 + (127/3150) u^3` from `t_c = u + (7/10)u^2 + (41/75)u^3` and `tau = 1/2F1(1,4/5;9/5;u)`; the
  falsifier values `1/16, 9/512, 75/8192, 3675/524288`; the off-K values; primality of every modulus printed (all prime).
* Binding rules: no file, script, log, folder, code, cluster, stage, agent, referee or status word in either appendix (grep);
  no `\draftnote`; no new macro; labels unchanged.

## 2. Issues found and fixed (file, line before the fix; what was done)

### 2.1 Appendix B

1. **Wrong pin** (`b_signs.tex` l.54, `b_signs_2.tex` l.95): "[LP20, s.2.4.3]" for the normalisation of asymptotic series. In the v2
   PDF this is (18) in **s.2.4.4** (s.2.4.3 = Prop. 4 and (17)); (18) gives `R_0 = L`, `R_1 = (3/20)(L - L^5)`, `R_2 = (9L/800)(1-L^4)^2`,
   i.e. `delta_(0,k)` at `q = 0`. Changed to `\cite[\S2.4.4, (18)]{LP20}` in both places, with the three values. **BIBPLAN must add
   "s.2.4.4 (18)" to LP20** (see REVIEW_2 s.4).
2. **Notation clash, letter of the gauge** (`b_signs.tex` l.24-28): `P^{(i)}_a := Delta_i^(1/2) Q^(i)_a`, the source letter that NOTATION
   renamed (rule 6), next to the gauge `P` of the next paragraph. Removed: the edge and leaves are written with `Delta_i^(1/2) Q^(i)_a`;
   the derivation now states `e~_i (x) e~_j = Delta_i^(1/2) Delta_j^(1/2) e_i (x) e_j` and `Delta_i^(1/2) Q^(i)_0(0) = Delta_i^(1/2) eta(1, e~_i) = 1`.
3. **Index clash** (`b_signs.tex` l.28): `Dil_i(d)`, `d >= 2`, two lines before the dilaton factor `d`. Index renamed `k`.
4. **New symbol** (`b_signs.tex` l.30): `\mathsf J` (permutation matrix) clashed with the J-function, the vector `J`, `J^C` of
   `conv:tra-forms` and `J_k` of `tab:sgn-model`. Replaced by `P^vec = eta^-1 P eta` exactly as in `conv:tra-forms`.
5. **Frame-letter clash** (`b_signs_2.tex` l.50): `beta_a := 4 - nu_a` (TRM's letter, renamed `wt` in NOTATION) collides with the frame
   parameter `beta`. Replaced by the explicit vector `(4 - nu_a) = (4,5,1,2,3)`.
6. **False sentence** (`b_signs_2.tex` l.58): "normalised stationary-phase series have constant term 1 (def:fsp-wick)": by
   `def:fsp-wick` the normalised stationary phase `\SPn` has constant term `a_0(c)`; only `\SPh` has constant term 1. Corrected.
7. **Index clash and dangling pointer** (`b_signs_2.tex` l.60): `p_(5q,0)(0) = 0 for q >= 1` used the parameter `q` as an index, and
   "`m_5` whose value is recorded in sec:evi-univ" was false (the value stood in the opening paragraph of C). Index renamed `n`; the
   sentence now names `c(z) = 1 + m_5 z^5 + m_10 z^10 + ..` (`rmk:mai-stored`) and appendix C (`sec:evi-univ`) now gives both values.
8. **Mellin variable `s`** (`b_signs_2.tex` l.68, l.74): `Gamma(1 + 5s)` with `s` the Mellin variable, while `s = u^(1/5)` in `sec:sgn-var`.
   Added "the letter `s` is the Mellin variable, `sigma = zs` (sec:con-form), not `u^(1/5)`" and "(with `sigma = zs`)".
9. **Phases incomplete** (`b_signs_2.tex` l.66): the flat convention `e^(-W/z)` of `sec:con` (`lem:con-form`, `def:con-PhiSP`) was missing,
   and `S_0` was attributed to all of s08. Rewritten after `def:fsp-wick`'s last sentence.
10. **Hessian letter** (`b_signs_2.tex` l.70, l.72): `det(-H)` with the class letter `H`; app. A writes `\mathsf H`. Aligned.
11. **Overgeneral claim** (`b_signs_2.tex` l.72): "The normalised series do not depend on the branches (prop:gam-rational)": that
    proposition is about `c(z; L)` only. Restricted to `c(z; L)`.
12. **Pairing claim** (`b_signs_2.tex` l.37): "for legs with a pole at `p = 0` it does not [equal the circle pairing]" is stronger than
    `def:mod-pair` ("differs by the residue at 0, not zero in general"). Corrected.
13. **String sign** (`b_signs_2.tex` l.35): `+d/dx` is stated nowhere in the body (the string equation is not proved, cf. the
    `\draftnote` of `rmk:tra-consistency`, s04_transport_3 l.70). Marked "an orientation; the string equation is used in no proof".
14. **Seven-variable chart** (`b_signs_2.tex` l.96): "`dy` for `Gamma(1 + 5/z)`" is the value at `sigma = 1` only, and "(a sign flip)" was
    attached to `y = 5L` instead of `x_j`. Now "`dy` for the factor `Gamma(1 + 5 sigma/z)` and `dx_j/x_j` for `Gamma((zeta^j - sigma)/z)`;
    `x_j = zeta^j - L` (the opposite sign to `L_alpha - lambda_i`)", as in `rmk:gam-sevenvar`.
15. Precision: `K` -> `K_b` (as `conv:mod-kernel`); `Lambda` of the EO variation pinned to EO (5-13) ("a letter of that paper");
    the sign anchors now say what `conv:mod-kernel` says (`x = -p^2/2`; `F^TR_2 = -1/240`); `J_k` row: "at its critical point `x = 1`,
    `\SPn_(x=1)[lambda; x^k](z')`"; Bernoulli paragraph aligned with `def:fsp-stirling` (numbers vs DLMF polynomials); PPZ "Rem. 2.1,
    (6)" as in `def:set-graphsum`; vertex weight written `Delta_i^(g_v-1) = TFT_(g_v,n_v)(e_i,..,e_i)`; `tab:sgn-lit` floated above its
    run-in heading: the sentence now precedes the table (cosmetic).

Checked and right as written (no change): Givental's convention and `R^-1 = R(-z)^*`; PPZ Def. 2.2/2.13/Prop. 2.12, p.20; DNOPS
s.2.3 and Rem. 2.1 (G4_2 s.7, G4R s.1.3 R-1); `|Aut|^-1`; all of `sec:sgn-var`; the twist rule and [LP20, (10)]; the characters
paragraph; `tab:sgn-model` rows 1, 2, 4-9 (G4R s.1.1 R3; `prop:mod-Rtr`; Step 4 of `thm:mod-tr`); the EO dictionary; `sec:sgn-legs`;
`sec:sgn-tower` (with `def:tra-tower`, `def:gau-van`, `lem:set-modhom`, `conv:tra-forms`); `q and -q`; exponentials; the gap constants.
