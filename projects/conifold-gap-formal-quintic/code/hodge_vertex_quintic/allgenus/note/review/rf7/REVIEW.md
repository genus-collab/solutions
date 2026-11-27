# Review rf7: s09 (the Gamma identity and the normalisation) and appendix A (formal stationary phase)

Author-facing; not part of the note. Key `review/rf7`. Files refereed and edited: `NOTE/sections/s09_gamma.tex`, `_2`, `_3`, `_4`;
`NOTE/sections/a_formalsp.tex`, `_2`, `_3`, `_4`, `_5`. Continuation: `REVIEW_2.md` (appendix A lemma by lemma, the list of edits,
the changes that touch other sections, compilation). Nothing outside these files and this folder was modified.

**Read for this review.** OUTLINE s.0-s.0.2, OUTLINE_2 (s.0.3-s.0.5), OUTLINE_8 (s09), OUTLINE_9 (appendix A), DECISIONS (D1-D15,
Q1-Q7), DESIGN_CHANGES (+ _2, in particular open point P1), NOTATION (+ _2), LABELS_2 (s09, appendix A), BIBPLAN; the digest DD2
(DIGEST_D2 .. _5: D1-D10, GAM-1..GAM-18, SP-*, proofs 3A-3B, conventions s.4); the sources FSL, FSL2, FSL3, FSL4 (FORMAL_SP_LEMMAS
s.0-s.10), RL2 s.5 (the corrections RL-C1..C8), CP2 (CHART_PROOF_2, all of s.2), RA s.3 (Theorem LMI and (L2')), LF s.0-s.4
(LF-C1..C6), LA s.0-s.3 and LA2 s.4-s.8 (LA-C1..C7), G5 s.1 (`thm:con-Pa`, the generic point); the s08 statements that s09 uses
(`def:gam-c`, `prop:gam-rational`, `lem:con-form`, `def:con-PhiSP`, `thm:con-Pa`, `thm:con-Pb`, `thm:con-glue`, `cor:con-A2`) and
`thm:set-qrr`, `lem:set-hrot`, `lem:gau-HPF`, `thm:gau-R0`, `cor:gau-colours`, `def:gau-Fg`.

## 0. Verdict: **sound-after-fixes**

* Every statement of s09 and of appendix A was checked line by line against its source in final corrected form and for internal
  correctness. I found **no mathematical error and no gap**, in the chain or outside it. The fixes are precision of statements
  (frame scope of `thm:gam-D`, definedness of the actual free coefficients, the exact meaning of "Stirling form" in the LMI amplitude,
  the dependency sentence of `thm:gam-main` (iv)), one notation clash (the disc `D` against the operator `D = q d/dq`), one citation
  outside BIBPLAN (DLMF 5.11.8 in s09), one unproved display (now proved in one line), wording of three remarks and of one sentence of
  appendix A ("proof of record").
* **No `\draftnote` was added.** The chain step C12 (Theorem Gamma by LMI, then `D = 1`) carries no `\draftnote`. The three existing
  `\draftnote`s of s09 (`rmk:gam-sevenvar`, `rmk:gam-Mcomplex`, `rmk:gam-prior`) are outside the chain, planned by OUTLINE_8 and worded
  without project references; appendix A has none.
* **The known point (LMI step 2 with the formal implicit function theorem; DESIGN_CHANGES_2 open point P1) is resolved: the note's
  step 2 is complete and correct, and in fact needs no implicit function theorem at all** (s.1 below). Neither the holomorphic IFT nor
  a pinned textbook statement has to be restored; the conditionality of Theorem A stays exactly as in GAPS_SYNTHESIS_2 s.3 and DEC Q2.

## 1. Theorem LMI, step 2: the formal route (the point the task names)

The worry of P1 was: with `lem:fsp-ift` one obtains the family of chart critical points only as a formal series in `Q`; identifying
it with the critical points of step 1 at `Q != 0` and continuing the coefficients holomorphically to `L = s_*` seem to need convergence.
The note's version avoids the problem entirely, because it never parametrises by `Q`:

1. **The family is explicit and parametrised by `L`.** `p(L) := (L, mu_c(L))`, `mu_c(L) = -phi'(L)/a_*`, is defined and holomorphic on the
   whole disc, `L = s_*` included. That `p(L)` is a critical point of `Psi_(Q(L))` for every `L` in the disc is a two-line computation
   (`d_sigma Psi = phi' + a_* mu = 0`; `d_mu Psi = a(L) - Q(L) e^(mu_c(L)) = a(L) - a(L) = 0`, using `Q(L) = a(L) e^(phi'(L)/a_*)`). Its
   Hessian `[[phi''(L), a_*], [a_*, -a(L)]]` has determinant `-a_*(a_* + (L - s_*) phi''(L))`, non-zero on the disc by the LMI hypothesis.
   I re-derived both by hand. Source: LA s.1.3 ("for `Q = Q(L)` the point `(L, -g'(L)/alpha0)` is critical ... it is the translate of
   `(L, nu_c(L))`"), LA2 s.4.1 ("`Q e^(mu_c) = alpha0 e` exactly, `mu_c = -g'(L)/alpha0`").
2. **Identification with step 1 (LF-C1) for `L != s_*`, i.e. `Q != 0`.** The translation `chi = mu - lambda/a_*` sends `p(L)` to
   `(L, chi_L)`, `chi_L = log a(L)` the branch used in `lambda(L)`: `mu_c - lambda/a_* = -phi'/a_* + phi'/a_* + log a(L)`. This is a direct
   computation (the converse direction of LF s.3.2, written in LA s.1.3); no uniqueness statement and no IFT are used. With
   `lem:fsp-cov` (Jacobian 1), `lem:fsp-const` (the constant `lambda s_*`) and the cancellation of the constant `Q^(M_*)` in `SPh`
   (LF-C2) this gives (eq:gam-chart) at every `L` in the punctured disc.
3. **Holomorphy and the value at `L = s_*`.** In the coordinates `s = sigma - L`, `m = mu - mu_c(L)` the chart phase is
   `const + [phi(L+s) - phi(L) - phi'(L)s] + a_* s m - a(L)(e^m - 1 - m)` and the normalised amplitude is
   `e^(M_* m) FF(L + s, z)/FF_0(L)`: every Taylor datum is an explicit holomorphic function of `L` on the whole disc, and so are the entries
   of `-H(L)^-1`. By `lem:fsp-univ` each coefficient of the chart series is a polynomial over Q in finitely many of them, hence
   holomorphic on the disc, and equal to `c^FF_n` off `s_*` by (eq:gam-chart). At `L = s_*` the same universal polynomials evaluate the
   data of `Psi_0` at `(s_*, mu_*)`, which gives (eq:gam-atQ0); step 3 (`lem:gam-deltasrc`) then gives the value. This is LF s.3.2's
   "explicitly, by (F1) the chart coefficients are polynomials in the Taylor data ..., holomorphic at `L = s_*` by inspection".
4. **Where `lem:fsp-ift` enters.** Only (a) inside `lem:fsp-fubini` (i) in step 1 (the inner critical point), and (b) in the sentence
   that `p` is the unique formal family of chart critical points through `(s_*, mu_*)` over `C[[L - s_*]]`. I checked (b): the gradient
   at `(s_*, mu_*)` is `(0, -Q(L) e^(mu_*))`, in the maximal ideal; the Hessian there is invertible (constant term of the determinant
   `-a_*^2`); the Taylor series of `p` is a formal zero of the gradient because Taylor expansion commutes with composition. (b) is correct
   but is not used afterwards; the note now says so explicitly (edit E6 of REVIEW_2), and the list of inputs after the proof says where
   `lem:fsp-ift` is used (edit E8).

Conclusion: P1 needs no `\draftnote`; parametrising by `L` (not by `Q`) is the formal-enough route, and it is written in full in the note.
The only analytic facts used are elementary (derivatives, sums, products and quotients by zero-free functions of holomorphic functions
are holomorphic); in `thm:gam-main` (ii) and `thm:gam-ginf` one also uses that a rational function with a holomorphic extension at a
point has no pole there. This is the "elementary" continuation of FSL3 s.8 warning 2, now without any IFT.

## 2. s09, statement by statement

| item | checked against | result |
|---|---|---|
| `sec:gam-series` (recall of `def:gam-c`, `prop:gam-rational`, `lem:con-form` (d); measures; the degenerate factor at `L = 1, zeta^i, 0`) | DD2 D6, GAM-1, GAM-2, GAM-5; LA s.1.5; s08 | correct |
| `lem:gam-deltasrc` (L2'): Hessian, `C_ss = 0`, `C_sm = -1/a_*`, `C_mm = P''(0)/a_*^2`, the moment formula with "`= 0` for `r < a`", the resummation, finiteness per order, the substitution `s -> delta` | RA s.3.2; LF s.2; LA s.1.4; DD2 3A | correct; I re-derived the moment formula from the generating function (`<s^a m^(a+2j)> = r! C_sm^a C_mm^j/(2^j j!)`, equal to the printed form) |
| `def:gam-lmidata` (LF-C3 = LA-C3; `phi_k`; `FF_0` zero-free and `a_* + (sigma - s_*) phi'' != 0` on a disc, equivalent after shrinking to the outline's punctured-disc form; `a` exactly linear; `Q` single-valued) | DD2 D7; RA s.3.1; LF s.1 | correct; "Stirling form of `Gamma(M_* + a/z)`" made exact (E2) |
| `thm:gam-lmi` statement (i)-(iii), the value formula and its re-expansions | RA s.3.1; DD2 GAM-4 | correct |
| Step 0 (Taylor coefficients `(-1)^k a_*/(k(k-1)(L - s_*)^(k-1))` of `a log a - a`; branch-free normalised amplitude; `-1/S''(L)`) | LF s.1; LA s.1.1 | correct |
| Step 1 (formal integrand, LA-C1; inner critical point `chi(sigma)`; Fubini factor `(a(sigma)/a(L))^(-1/2)`; `lem:fsp-cov` with `chi -> e^chi`, Jacobian `a(sigma)`; `lem:fsp-stirling` over `C[[sigma - L]]` with the unit `a(sigma)`, `M = M_*`, measure `dx`; outer amplitude `a(L)^(1/2) A`) | RA s.3.3; LF s.3.1; LA s.1.2 | correct; the linearity used is over `C[[sigma - L]][[z]]`, not `C[[z]]` (E5) |
| Step 2 and holomorphy | s.1 above | correct, complete |
| Step 3 (`Psi_0 = phi(s_*) + phit(s) + a_* s m` exactly; source `e^(M_* m)`; constant term `e^(M_* mu_*) FF(s_*, 0)`) | RA s.3.3; LF s.3.3 | correct; the constant `phi(s_*)` now removed by name (E7) |
| `rmk:gam-lmi-words` | RA s.3.1 | correct (first pole of `Gamma(M_* + a/z)` at `sigma = s_* + delta`) |
| `thm:gam-main` (ii): data `s_* = 1`, `a_* = -1`, `M_* = 0`; `phi`; `FF` with `(5 sigma)^(1/2)` (LF-C4 = LA-C4 = REV_2 M10); `S = SS_0`, `A = A_St`, `lambda` = the `lambda` of `def:gam-c`; value `exp(Sum bk_k T_k z^(2k-1)) = h` | LF s.4.1; LA s.1.5 | correct; `T_k = -Sum_(a in w_0) a^(1-2k)` and `h_1 = -3/20` re-checked |
| `thm:gam-main` (iii): Taylor data at `zeta L`, the substitution `t = zeta t'`, `lem:fsp-homog` with the unit `zeta`, the values at `zeta^i`, exponents `= 4n mod 5` | CP3 s.3.9 via DD2 3A; RF2 s.4 | correct; the three Taylor formulas re-derived by hand |
| `thm:gam-main` (iv) and its proof sentence | DD2 GAM-5, GAM-8, GAM-9 | correct; the dependency sentence was too short (E10) |
| `prop:gam-hatflat`, `cor:gam-equiv`, `rmk:gam-K8` | G5 s.1.2 (`thm:con-Pa` (iv)); DD2 GAM-7; DD1 rmk:glu-K8 | correct, no circularity (`thm:con-Pa` uses neither) |
| `thm:gam-D` | ST5 s.2 via OUTLINE_8; `thm:con-glue`, `thm:con-Pb`, `cor:gau-colours`, `def:gau-Fg` | correct; frame scope added (E12) |
| `thm:gam-ginf`: `s_* = 0`, `a_* = 5`, `M_* = 1`, `delta = -z/5`, `Q^5 = (-q)^(-1)`; `phit_j(delta)/z = (w + theta) log(1 + theta/w) - theta`; `lem:fsp-dlmf5118`; `Sum_j zeta^(-jn)`; `5 beta_(5l)(1/5)` | LF s.4.2 (LF-C5); LA s.1.5, s.2.3 | correct; each line re-derived |
| the symbol display after `thm:gam-ginf` | LF s.4.2 (C5) | correct, but stated without proof: one-line proof added (E14) |
| `cor:gam-poly`, `cor:gam-free` | LF s.4.3; RA s.3.4; RF-C3; RS2 item 38 | correct; `B_6(1/5)/6 = 4537/3937500` re-computed exactly; the degree sentence and `frak f^act` made precise (E15, E16) |
| remarks `rmk:gam-sevenvar`, `-blind`, `-necessity`, `-Mcomplex`, `-general`, `-prior` | DD2 GAM-11, -13..-17; RF-C1; GEN s.3.2 Lemma 4 | correct; `det H_7 = 1/(L^2(1 - L^5)^2)` (Schur complement `SS_0''(L)`), `1/det Hess psi = -L^2 [5]_L^2`, `det Hess psi(1) = -1/25`, the constraint `x_0 = -q y^5/(x_1..x_4)` and the chart phase re-derived; (a) of `rmk:gam-necessity` re-proved; wording of `-blind`, `-necessity` (a), `-general` (E17-E19) |

Citations of s09 against BIBPLAN: [GJR18] Rem. 7.2, Lemma 7.7 (37); [MR1901075] s.10, 10.4 Lemma; [MR2553377] Lemma 4.20 (at `z = -1`,
the DL reading); [MR2276766] s.7 (11); [DLMF] 24.4.1. All pinned. The one exception was "[DLMF, 5.11.8]" in the proof of
`thm:gam-ginf`, which BIBPLAN allows only in `rmk:fsp-classical`; removed (E13). Section citations now use `\S` as elsewhere (E20).
No finite computation is cited by a proof; the remarks point to appendix C only with "not used" or "known only by computation".
