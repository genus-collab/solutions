# Review rf3: section s03, part 2 of 2

Continuation of `REVIEW.md` (verdict s.0, method s.1, findings and fixes F1-F12 s.2).

## 3. Checked and found correct (no change)

Every item below was re-derived by hand; none relies on a computation.

**`sec:gau-cyc` ((CYC)).**
* `lem:gau-const`: the coefficient of `ell^(n-1)` gives `a_(n-1)' = -n a_n/u`, impossible (no `u^-1` in a derivative; characteristic 0).
* `lem:gau-adj`: `P_4(-D-1) = P_4(D)`; the two-condition criterion re-derived (coefficient of `D^3` in `P^*`: `-A_3 + 4 D A_4`; of `D^1`:
  `-A_1 + 2 D A_2 - 2 D^3 A_4`; the conditions at `D^2`, `D^0` follow); `2 D A_4 = A_3`, `D A_2 - D^3 A_4 = A_1`.
* `lem:gau-conc`: telescoping; alternation (top total degree `s` of `CC(f,g) + CC(g,f)`); the jet formula and all ten entries
  (`JJ_00 = 0`, `JJ_01 = (2/5)(1-u)`, `JJ_02 = 1-u`, `JJ_03 = -u`, `JJ_11 = 0`, `JJ_12 = u`, and antisymmetry); `det JJ = u^4` (reversal of
  four letters is even). Hand values: `(f, Df, D^2 f)(0)` for `f_0, f_1, f_2`; `CC(f_0,f_1) = -7/5`, `CC(f_0,f_2) = 2`, `CC(f_1,f_2) = 0`;
  `D^k ell`; the orders of `Y_0..Y_3`; the cancellation of `u^-1` in the first bracket; `CC(f, y_3) = (Df)(0)`, so `CC(f_1,y_3) = -1`,
  `CC(f_2,y_3) = 0`; `Pf = 0 - 2(-1) + * . 0 = 2`, `det = 4`; Abel's identity `D(u^2 det Wr) = 0`, `(u^2 det Wr)^2 = 4`. The entry
  `CC(f_0, y_3)` is written `*` and is not needed; the sign of `u^2 det Wr` is not claimed (rule 3).
* `lem:gau-transv` (exchange-matrix triangularity, `(JJ^-1)_03 = 1/u`), `lem:gau-iso` (twisted Leibniz rule for `sh`;
  `CC(y_3, varpi_1) = 1`; kernel `(0, 1, 7/10)`; `gamma' = 7/10 = gamma` is the only isotropic value), `lem:gau-sg` (the expansions of
  `Upsilon_k`; `D^2`, `D^3` in terms of `partial_t`; the Darboux basis, `G^-1 = -G`; `Upsilon~_3 = Yuk`), the remark on `Yuk`
  (`Yuk` is the expansion of `kappa^LR_1`; `eta(H*H, H) = 5 kappa^LR_1`; `kappa_1 = Yuk D t_c`).
* `thm:gau-cyc`: `D o (Dt)^-1 o D = (Dt) partial_t^2`; leading coefficient `u = A_4`; the Wronskian argument; the scaling factors
  `mu^-1, mu/lambda, lambda^2, mu/lambda, mu^-1`. Scope: any power-series `y_0` not proportional to `varpi_1`, `y_0(0) != 0` not used
  (RS1 A1), any scale of `varpi_1`, every basis of an isotropic plane; `PF4 = M_cyc(varpi_0, varpi_1)` in every conifold frame.
* `cor:gau-cycsys`: the operator identity with the five `kappa_a^-1`; both directions; the automorphism `D -> D + w/z` of `\cR[D]`.

**`sec:gau-systems`, `sec:gau-U` (Theorem U).**
* `def:gau-onepole`: the entries of `eta A` and `A^T eta`; self-duality iff `b_0 = b_2`, `b_3 = b_4`.
* `prop:gau-quintic`: the representatives of `C_1 kappa^LR`; the `b_a`; `r = -1`; product `-1/(u(1-u)^5)`; (iii) the unit shape iff
  `y_0(0) != 0`. `prop:gau-equalprod`: `b~`; (ii) by uniqueness of roots in `1 + uQ[[u]]` and `d/du (s^4/tau) = (4/5) s^-1 (1-u)^-1`;
  (iii) `A^5 = Bprod 1` on the leading non-zero `P_n` (a strengthening of the outline's `P_0 != 0`, as G2R s.2.4 allows).
* `lem:gau-comm` (a)-(g), `lem:gau-lift`, `lem:gau-scal`, the proof of `thm:gau-U` (a), `lem:gau-lead`, the proof of (b), (d),
  `lem:gau-degree`, `prop:gau-graded`: identical in substance with G2_2 s.3-s.4, G2_3 s.5.1-s.5.5, G2_4 s.6.1-s.6.4, with the referee
  corrections of G2R_2 s.3 applied (C2: every field `\Bbbk`; C5: the free constant is sans `c`, `u Bprod = Bprod^circ`, the field
  `\Bbbk`, `zeta^(gr)`, `\Gau`; R1: shifts `m = 1..5`; R2: "for every `P^nat` with `P^nat_0 = frak d Theta_U`"). Re-derived: the Gram
  determinant `5 (5 Bprod)^4`; the pivots `n + j/5`; the window indices `2-m, .., 1` in `lem:gau-lift`; `frak d' = -(1/5) tr(dlog Theta_U)
  frak d`; the order of choices in the construction; `Gau(A~, A~) = \Bbbk[[z]] 1`.

**`sec:gau-van`, `sec:gau-symp`.**
* `lem:gau-R1`; `lem:gau-diag` (general formula `c(0)^5 det Theta_U(0) = p_(0,1)(0)^5 prod rho^U_a`; the exponent bookkeeping `t_c'`: 0,
  `v'`: -10, `varpi_0`: -5, `(1-u)`: -10; `p_(0,0) = tau^(-1/4)/varpi_0`; (N0) stated as a hypothesis, true in the conifold frames).
* `def:gau-van`, `lem:gau-vanmod` (all five rows; (Mod2) at `m = 0`; `Y^C(0,z) = -e_1`), `lem:gau-vanshape` (exponent, sign,
  denominator, the five lowest terms, homogeneity, `Y^C(v(u),z) in Q[z^-1][[u]]^5`), `lem:gau-uniq`, `lem:gau-vanq` (`e(n,a)`;
  `nu_(a+1) = nu_a - 1 + 5 delta_(a1)`; the correction of DB S24 (b): the pole meets an order `>= 4`; `[z^-1] Y^C_Q = -t_c e_0`).
* `thm:gau-G` (c) via (R1) and `lem:gau-uniq`; (a) uniqueness via `thm:gau-U` (b), (d) (G2_4 s.8); (b); `def:gau-can`
  (existence, uniqueness, graded, rational, `Rot_i P_can = P_can`, `p_(0,a)(0) = 1`); `prop:gau-symp` (the equation for `S`,
  `T = eta^-1 S` solves the gauge equation of `(A~, A~)`, entry `(2,1)` with `eta_(a1) = 5 delta_(a2)`; graded freedom
  `+-exp(a_1 z^5 + a_3 z^15 + ..)`); `cor:gau-Gstar`.

**`sec:gau-R0`.**
* `lem:gau-HPF`: `(zD + L)(s^e F) = s^(e-1) \cO_e F`; the iterates; the elementary symmetric functions `(3, 17/5, 9/5, 274/625,
  24/625)`; the coefficient `4n - 5 theta_u` of `ph_n`; `Phi^flat_(j,1) = (1-u) j(j+1)/10`; `r_1 = 3/5`; `ph_1 = -3/20`.
* `lem:gau-poly` (i)-(iii), including the triangular system of the free coefficients; `lem:gau-wkb` (the `z^0` coefficient
  vanishes; `T^w_1 f = zeta^(4i)(5 L^-1 D f + (L^-1 - L^4) f)`, so `T^w_1(L g) = 5 zeta^(4i) D g`; the recursion; (iii) homogeneity of
  degree 5 of the `z`-monomials, `zeta^(-5i) = 1`).
* `lem:gau-colmod` (scaling identity, the five relations, pivot `4n zeta^(4i)`); `lem:gau-colrows` Steps 1-5 (G2_5 s.9.4 with G2R_2
  C3: Step 2 uses the inverse substitution through `lem:gau-wkb` (iii)); (QRR) only in Step 3.
* `thm:gau-R0` (i)-(iv) (after F1): the same exponential `e^(u_i/z)` is stripped on both sides (PS-C8), `U_i(v(u)) = u_i - u_i^c`;
  `mu(b,k)` an integer `>= max(0, (4k-4)/5)`; `ph^can_0 = 1` from `lem:gau-diag`; the injectivity of `L -> (1 - 5^5 q)^(-1/5)`.
  `cor:gau-rowsNN` (the `logq`-free subring is `D`-stable; antisymmetry of `JJ`); `cor:gau-colours` (the three-line computation read
  without its first equality, G2R_2 s.2.9; uniqueness by the non-zero-divisor `Q^(i)_0`); `def:gau-Fg`.

**Binding rules.** Rule 2 (all-orders steps in sentences): yes. Rule 3: no proof uses a computation. Rule 5: no bookkeeping in
the `.tex` files. Rule 6: after F4, F5. Rule 7: frames, `y_0`, `\Bbbk`, the fields of `D_0` (now `Q(zeta)` before the glue) and `D_i`.
Rule 8: (QRR) enters s03 only in `lem:gau-colrows` Step 3, (M) and (G0) only through s02, as the convention paragraph says. Rule 10:
after F6.

## 4. Open items and observations (no change made)

* O1 (priority): the proof of `thm:con-glue` (s08) still assumes the rationality that `thm:gau-R0` (iii) no longer asserts; the
  required edit is X1 below. Until it is made, the rationality of `D_0` is asserted, not proved, in s08 (it is not used anywhere
  else before `thm:gam-D`).
* O2: `Phi^circ` is defined in `lem:gau-poly` (iii) and used nowhere in the note.
* O3: s03 occupies 24 pages of the compiled note (pages 29-52); the design budget is 10 pp (OUTLINE_9).
* O4: `def:gau-van` (ii) cites `lem:gau-vanq` (b) for the well-definedness of the iterated integrals (forward, not circular).
* O5: some one-proof local letters reuse global ones (`s`, `q_ab` in `lem:gau-conc` (ii); the plane `Lambda` in `lem:gau-iso`,
  `lem:gau-sg`; `S`, `T` in `prop:gau-symp`; `mu` in `thm:gau-G` (b)); declared or plainly local, left as they are.

## 5. Changes that affect other sections (not edited here)

* X1 (needed) `s08_connection_8.tex` l.18, proof of `thm:con-glue`: (b) "In `Q[[q]]`" -> over `Q(zeta)[[q]]` (`A^act` is a WKB
  amplitude over `Q(zeta)[[q]]`, `lem:gau-wkb` (ii) with `\Bbbk = Q(zeta)`); (c) `A^act_n in Q(zeta)[L, 1/L]`; (d) the constant series
  `D_0 in 1 + z Q(zeta)[[z]]`, the embedding `Q(zeta)[L,1/L][[z]] -> Q(zeta)((s))[[zhat]]`; (e) evaluation at `L = 1` on
  `Q(zeta)[L,1/L]` gives `h D_0 = Phi^can(z,1)` in `Q(zeta)[[z]]`; as `h` and `Phi^can(z,1)` are rational and `h(0) = 1`, `D_0 =
  Phi^can(z,1)/h(z) in Q[[z]]`, hence `A^act = A^can D_0^-1` is rational. The statement of `thm:con-glue` stays as it is.
* X2 (wording) `s10_main_2.tex` l.44: "It fixes the flat amplitudes `A_n in Q[L, L^-1]` (`thm:gau-R0` (iii))": `thm:gau-R0` (iii) now
  gives `Q(zeta)[L, L^-1]`, which is all the argument needs (the endomorphism fixes `L`); write `Q(zeta)[L,L^-1]` or add `thm:con-glue`.
* X3 (wording) `s01_intro_2.tex` l.40-42: the display ends with `D_0 in Q[[z]]` and cites `thm:gau-R0`, `cor:gau-colours`; true, but
  the rationality is now proved in `thm:con-glue`: add that label or drop the clause. NOTATION table 4 and OUTLINE s.0 rule 7
  ("`D_0` in `Q[[z]]`") remain true as final statements.
* X4 (scope, recommended) `s02_setup.tex` l.76, `lem:set-pf5`: "on any differential `Q[[q]]`-algebra" is narrower than its use in
  `cor:gau-cycsys` (rings containing `Q((u))`); suggest "on any differential ring containing `q` with `D q = q` in which `z` is a
  constant". s03 now records this in the proof of `cor:gau-cycsys` (F3), so no gap remains.
* X5 (no change needed) `s02_setup_5.tex` l.58, `def:set-confexp`: its use of `lem:gau-conc` (i) and `lem:gau-adj` in
  `Q[[q]][logq]` is now covered by their statements (F2).
* X6 (notation) `\Lambda`, `\Lambda'` (s03 `thm:gau-R0` (i), s08 `lem:con-N`, NOTATION table 4) versus `\Lambda_s`, `\Lambda_varsigma`
  (s10 `lem:mai-Ghat`, `cor:mai-natfact`, s10_main_2 l.17): one notation for the hatting matrices should be chosen.
* X7 (notation) `PF5^(0)` (s03 `lem:gau-HPF`, `lem:gau-wkb`) and `PF5^L` (s08 `thm:con-Pa`; NOTATION_2 table 8) are the same operator;
  s08 could say `PF5^L = PF5^(0)` of `lem:gau-wkb`, or one symbol be used.
* X8 (notation) s08 `lem:con-N`, `thm:con-conn` use `q` as an index (`z^(5q)`, `K_q`), clashing with the parameter `q`; s03 renamed its
  analogous index to `m` (F5).
* Statements of s03 whose wording changed: `def:gau-conc` (iii), `lem:gau-adj`, `lem:gau-conc` (i) (wider scope); `lem:gau-confexp`
  (`I_1`, `Xi`); `lem:gau-diag` (scope of the last display); `thm:gau-G` (a) (dummy `P'`); `lem:gau-colmod` (`F^(i)`); `thm:gau-R0`
  (iii), (iv), `cor:gau-rowsNN` (i), `cor:gau-colours` (field `Q(zeta)`). Every use of these labels in other sections was read
  (grep over `sections/`); apart from X1-X3 none depends on the changed parts.

## 6. Compile record

Scratch copies of `main.tex`, `preamble.tex`, `macros.tex`, `refs.bib`, `sections/` in
`<user home>/AppData/Local/Temp/claude/d7c1b608d03d194a/37f24460-3601-4c85-b123-5156e19eb475/scratchpad/allgenus_review_review_rf3/`
(`base/` unmodified, `mod/` with the fixes); pdflatex, bibtex, pdflatex, pdflatex, each with a timeout.
* Baseline: 153 pages, no error, no undefined or multiply defined reference, 15 overfull boxes (one in s03, s03_gauge_7 l.63).
* After the fixes: 155 pages, no error, no undefined or multiply defined reference, 14 overfull boxes, none in s03.

## 7. Files changed

`NOTE/sections/s03_gauge.tex`, `s03_gauge_2.tex`, `s03_gauge_4.tex`, `s03_gauge_5.tex`, `s03_gauge_6.tex`, `s03_gauge_7.tex`
(`s03_gauge_3.tex` unchanged). New: `NOTE/review/rf3/REVIEW.md`, `NOTE/review/rf3/REVIEW_2.md`.
