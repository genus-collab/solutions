# Review rf6: s07 and s08, part 2 of 2

Continuation of `REVIEW.md` (verdict s.0; the known issue F1-F2 s.1; s07 F3-F6 s.2; s08 F7-F8 s.3).

## 0. Further findings and fixes in s08

**F9 (minor; scope; Theorem P(a) (i), Step 4, `s08_connection_4.tex` l.80-85).** `a := (sigma^5 - 1) A_St` and `eq:con-bz`
(`b_z = q prod_k (5 sigma + kz) A_St`) were claimed "in `A_L[[t]][[z]]`", but `A_St` contains the half-density `e^(S_1)`, and
`A_St(L, 0) = (5L/(1-L^5))^(1/2)` is not in `A_L`. Fixed: `a := (sigma^5 - 1) A_St/A_St(L,0)`, `eq:con-bz` with `A_St/A_St(L,0)` on the
right (both sides in `A_L[[t]][[z]]`, the Taylor data being those listed after `def:gam-c`); and the exponent
`(S(L+t+z) - S(L+t))/z` is noted to lie in `A_L[[t,z]]` without constant term (`S'(L) = 0`). The conclusion
`c<q prod(5 sigma + kz)> = c<sigma^5 - 1>` is unchanged (G5 s.1.3 Step 4; G5R F1.3 was already applied).

**F10 (minor; `prop:con-deg` statement, `s08_connection_2.tex` l.5).** "The value of `L^-1 A_n` at `L = 0` is `f_(n/5)` if `5 | n` and 0
otherwise" is false at `n = 0` (`A_0 = L`). Now "for `n >= 1`", with `A_0 = L` stated.

**F11 (minor; `rmk:con-points`, l.39).** "for a colour-0 WKB function with coefficients in `Q[[u]]`" lacked "solving HPF" (in the
section's conventions a WKB function is any series with `ph_0 = 1`); added. The actual amplitude is not such a series a priori; the
remark now says that its coefficients are Laurent polynomials in `L` (`thm:gau-R0` (iii), `thm:con-glue`).

**F12 (notation, typography).**
* `lem:con-N`, `thm:con-conn` and the closing paragraph: the index `q` of `K_q`, `z^(5q)`, `u^(4q)` clashed with the parameter `q` (rf3
  X8); renamed `r` (declared local). `K_r - f_r` depends on `f_1, .., f_(r-1)`; proofs unchanged.
* Conventions of s08 (l.12): `PF5^L` is now identified with `PF5^(0)` of `lem:gau-HPF`, `lem:gau-wkb` (rf3 X7).
* `thm:con-Pb`: the domain is stated, `Phi^SP = Phi^can` in `Q[[u]][[zhat]]` (G5 s.1.2).
* `lem:con-u0` proof: `r_n(0) := [u^0] r_n` defined (was implicit), which also removed a 20 pt overfull line.
* `def:con-PhiSP`: the "equivalently" formula displayed (46 pt overfull line removed); `thm:con-glue` display split (1 pt).
* Section title of `sec:con-form`: `\texorpdfstring` for `$1$` and `$\PhiSP$` (removes the hyperref warnings).
* Three overfull boxes of 2-4 pt remain (s07 part 2 l.9; s08 part 4 l.16; s08 part 6 l.27); cosmetic.

## 1. Checked and found correct (no change)

Every item below was re-derived by hand; corrections already applied in the sources (G5R F1.1-F1.4, F4.1-F4.3; FR C1-C6) were
confirmed in the text.
* s07: `def:flh-outer` (i)-(iv) (finiteness for `q + z1 in z^2 V[[z]]`; the dilaton series `(1-eps)^-(2h-2+n)`; gluing of `varkappa.Omega`,
  exponent additive; unit `varkappa^-1 1`); `lem:flh-frame`; `lem:flh-dil` (identity `d^n Fbar^(k.Omega)(q) = k^(2h-2+n) (d^n Fbar)(k q)` and
  the dilaton step at `k q`); `lem:flh-comm` (`Sum_v (2g_v - 2 + n_v) = 2g - 2 + n`); `lem:flh-hat` (edges `s^(4(a+b+1))`, leaves `s^(4k-3)`,
  vertices `s^(4D_v - (2g_v-2+n_v+l_v))`, dimension count `3g - 3 - E + L`, total `N = 10(g-1)`); `lem:flh-inputs` (i), (ii), (iv)
  (`d = (P_0)_33 = (P_0^vec)_00`, TFT fixed by its unit, `TFT = d.TFT^m`, `d = tau^(1/4) varpi_0` consistent with `lem:gau-diag`);
  `thm:flh-main` Step 1 bounds (`a + b <= 3g - 4`, leaves `<= 3g - 2`), Step 2 (`E_(R'')`, `q_(R'')`, `q_(R'') + z d^-1 1 = O(z^2)`), the six-line
  chain of Step 3, Step 4 against `lem:tra-conv` (f1) and `def:red-flatsum` (`d q_T + z1 = d y`); `rmk:flh-dilfactor`.
* s08, `sec:con-problem`: `eq:con-chain` (`(zD + L)(s^e F) = s^(e-1) O_e F`); `PF5^L(s^-1 Phi) = s^-1 HPF(Phi)` with the `e_m`;
  `lem:con-u0` (the triple-product identity at `(10,-4,-3)` and `(0,1,2)`, `O^0_(-5) - 3 zhat = O^0_10`, the T-system closing relation,
  `[D, v] = z`); `prop:con-deg` Steps 1-7 (`P_0 = 0`; left factor; `(5(zD+L) + kz) L^m`; `P_1 L^m = (1-m)(1-L^5) L^(m-1)`; weights; the
  induction with `deg R_n <= 4n + 5`, `r~_n in L Q[L]`, trivial resonance; `r~_n = L^(4n) r_n(L^-5)`).
* `sec:con-form`: `e^lambda = -q(L)`; `S_0'' = 5/(sigma(1 - sigma^5))`; the Stirling identity with `(2 pi)^3 z^2` (log z coefficient 2);
  generating function of `p_k`, `p_k in (1-L^5)^-k Q[L]`; the Taylor data of `S_0`, `S_1`, `S_2k`; `prop:gam-rational` (1)-(5) (branch
  changes compensated by `lambda`; `S_0` shifts by `2 pi i Sum_j zeta^j n_j`); `lem:con-form` (a) (`S(1/x) = -WW(x)`, `GG`), (b)
  (`-sqrt5 (1-x^5)^(-1/2) dx`, the `5 Sum_P binom(2k-2+5P, 5P)` sums, degree bound), (c) (`WWhat`, `WWhat'(X)`, `WWhat''(1) = 5/(1-u)`,
  `WWhat(1) = (5/4)/tau`), (d) (reflection, parity of `e`, and the direct route); `def:con-PhiSP` (Gaussian factor times `ahat(1;0)` is
  `(2 pi zhat/5)^(1/2)`; `VV in z^-1 Q[z^-1][[u]]`, z-power `<= -1 - 5l`).
* `sec:con-P`: `thm:con-Pa` (i) Steps 1-4 (`sigma_eps' = sigma_eps(sigma_eps^5 - 1)/5`; the square-root bookkeeping `= sigma_eps/L`;
  `eq:con-shiftD`; `PF5^L(Lc) = L(c<sigma^5 - 1> - c<q prod(5 sigma + kz)>)`; the shift multipliers `q prod_k (5 sigma + kz)/((sigma+z)^5 - 1)`),
  consistent with `lem:fsp-generic` (iii); (ii) Steps 1-5 (equivariance, `phi = L/X`, phases via `5uX^4/(1 - uX^5)`, half-density
  `((1-u)/(1-uX^5))^(1/2)`, Stirling parts, `mu = s^4`); (iii), (iv). `lem:con-cfu` (`R_1 = 1/2 + T^2/18`, `V_1 = 5/9`, `Y_1 = (1+T)/18`,
  `deg Y_k <= k`, `V = tau^(-5/4)`, the inverse); `lem:con-red` (reduction relation, grading `a - b + 1 = 1 mod 5`); `lem:con-S2` (the model
  system for `(J_0, -J_-1, -J_3, J_2, J_1)`, the chain rule `T_k`, symbols `(2 pi z/5)^(1/2) varsigma^(m-1)`); `lem:con-id` (a)-(d) and the
  final clause; `lem:con-relab` (exponent map, flat and hatted fibres, determinant 5); `lem:con-key` (`beta = j - 5i' - m - 1`,
  `E = i + l - (4 beta - 1)/5`); `lem:con-weight`; `prop:con-res` (`p = sY`, model residues `s^-(m+1) I^C_m`, signs and powers of `s` cancel);
  `thm:con-Pb` Steps 5-7 (`[u^1] VV = -1/z`, `[u^1](Y^C_Q)_0 = -1/z`, `s = 1`); `rmk:con-710` (`-(1/2 + 1/5)/z = -(7/10)/z`).
* `sec:con-glue`, `sec:con-finite`: `cor:con-A2`; `lem:con-N` (`mu(b,5r) = 4r`, the triangular relation for `K_r`); `thm:con-conn` (i)-(iii).
* Citations: s07 cites only [PPZ15] (Def. 0.3, Rem. 2.1 and (6), Defs. 2.2, 2.5, 2.13, Prop. 2.12; Props. 2.4, 2.14 only in
  `rmk:flh-fail` (iii), as not needed), all pinned in BIBPLAN; s08 cites nothing (no quoted input in `thm:con-Pa`, `prop:con-deg`,
  `thm:con-Pb`; (QRR) through `thm:set-qrr` in the glue only). Conditionality statements of both sections agree with GAPS_SYNTHESIS_2 s.3.

## 2. Open items (none on the chain)

* O1 (design records, not the note): OUTLINE_7 (`thm:flh-main`: "`[z^10] D_0 = ([z^5] D_0)^2/2` is forced by `prop:gau-symp`";
  `rmk:con-genus`: "`D_0(z) D_0(-z) = 1`, from `prop:gau-symp` for the actual and the canonical gauge, PR C2"), DIGEST_C2_3
  (`rmk:flh-fail` (5): "true for the actual `R'` with or without (F)") and FR s.1 rest on the untransferred `q`-series identity. For
  `design/SOURCE_ERRATA.md` (not edited here).
* O2 (notation, editor): the hatting matrices are `Lambda`, `Lambda'` in s03 and `lem:con-N` (NOTATION table 4) and `Lambda_s`,
  `Lambda_varsigma` in s07 and s10 (local letters there); one notation should be chosen (rf3 X6).
* O3: whether `D_0(z) D_0(-z) = 1` can be proved without Theorem Gamma (equivalently, the symplecticity of the expanded rows from that
  of `R` at large radius) is open; it is not needed. Stated in `rmk:con-genus` as "not available a priori".
* O4: local letters reused within one proof: `S` (= `P_can Psi^m`) in the proof of `lem:flh-inputs` (iii) (also `S` the phase in s08);
  `r` (index of `K_r`) next to the HPF remainders `r_n`; declared local, left.

## 3. Changes that affect other sections (not edited there)

* X1 `lem:flh-inputs` (iii) changed meaning (no unconditional symplecticity). Uses read: `s02_setup_8.tex` (`def:set-outer`, rf2's new
  paragraph cites (i) and derives the same defect formula: consistent; the old sentence citing (iii) for unconditional symplecticity is
  already gone); `s10_main_2.tex` (cites (iv)), `s10_main_3.tex` (cites (i), (iv); rf1 S10-19 cites "the case `D = 1` of (iii)"):
  consistent. The proof of `thm:flh-main` (Step 1) now relies on the symmetrised edge of `def:set-outer`; if that definition is changed
  again, Step 1 must follow.
* X2 `thm:con-glue` statement now "In `Q(zeta)[[z]]` ... hence `D_0 in 1 + z Q[[z]]`". Users: `s09_gamma_3.tex` l.89 (`thm:gam-D`:
  "`D_0 = Phi^SP(z,1)/h(z)` exactly in `Q[[z]]`": consistent), the remark after `thm:gau-R0` and after `cor:gau-colours` (s03, point to the
  glue for rationality: consistent), `s10_main.tex` ("Where the inputs enter"), `s01_intro_3.tex` l.46: consistent. rf3 X2 (`s10_main_2`
  l.44, "`A_n in Q[L,1/L]` (`thm:gau-R0` (iii))") and X3 (`s01_intro_2` l.40-42) should cite `thm:con-glue` for the rationality.
* X3 `c_evidence_2.tex`, paragraph "Independence": "only the odd coefficients `a_j` of `log D_0` are independent, because
  `D_0(z) D_0(-z) = 1`" uses the identity as known; it holds only a posteriori (`rmk:con-genus`). Suggest "would be independent if
  `D_0(z)D_0(-z) = 1`, i.e. if the expanded rows were symplectic". Same file, "The natural factorisation": "only the symplecticity of `G`
  detects `D = 1` (`rmk:con-why`)" should read "detects `c(z)c(-z) != 1`" (`rmk:con-why` as corrected, F8).
* X4 `s01_intro_3.tex` l.17 (PPZ pins of `sec:flh`): already corrected by rf1 (S01-3) to the list of s.1 above; consistent.
* X5 Labels: none added, renamed or removed in s07/s08. Equation labels kept (`eq:con-glue` now displays one identity; its users
  checked).

## 4. Compile record

Scratch folder `<user home>/AppData/Local/Temp/claude/d7c1b608d03d194a/37f24460-3601-4c85-b123-5156e19eb475/scratchpad/allgenus_review_review_rf6/`
(copies of `main.tex`, `preamble.tex`, `macros.tex`, `refs.bib`, `sections/`, refreshed from NOTE after the concurrent edits of the other
reviews). pdflatex, bibtex, pdflatex, pdflatex (MiKTeX, `-interaction=nonstopmode -halt-on-error`, each with a timeout): no error, no
undefined reference or citation, no multiply defined label, no rerun request; 159 pages. Baseline before my edits (06:30 copy): 154 pages,
no undefined reference. In s07/s08: three overfull boxes of 2-4 pt; no hyperref warning.

## 5. Files changed (NOTE/sections only)

`s07_comparison.tex` (intro; paragraph after `def:flh-outer`; `lem:flh-dil` proof), `s07_comparison_2.tex` (`lem:flh-inputs` (iii),
(iv) and the paragraph after it; `thm:flh-main` Steps 1-2; paragraph after the theorem; `rmk:flh-fail` (i), (ii), (v), (vi);
`rmk:flh-dilfactor`), `s08_connection.tex` (conventions; `rmk:con-genus`; `lem:con-u0` proof), `s08_connection_2.tex`
(`prop:con-deg` statement; `rmk:con-points`), `s08_connection_3.tex` (subsection title), `s08_connection_4.tex` (`def:con-PhiSP` display;
`thm:con-Pa` (i) Step 4), `s08_connection_7.tex` (`thm:con-Pb` statement), `s08_connection_8.tex` (`thm:con-glue`; `cor:con-A2` proof;
`rmk:con-why`; `lem:con-N`, `thm:con-conn` index). Line counts: 81, 115; 59, 42, 80, 88, 61, 58, 78, 72 (all below 250). New files:
`NOTE/review/rf6/REVIEW.md`, `REVIEW_2.md`.
