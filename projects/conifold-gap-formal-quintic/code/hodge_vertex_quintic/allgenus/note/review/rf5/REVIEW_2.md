# Review rf5, part 2: statement-by-statement record, what remains open, changes affecting other sections, compilation

Continuation of `REVIEW.md`. Notation as in the note (ASCII renderings).

## 1. Statement-by-statement record (all re-derived by hand unless said otherwise)

* `conv:mod-kernel`. The kernel equals DNOPS (2-28) (`(1/2) Int_z^(sigma z) omega_(0,2) / (omega_(0,1)(sigma z) - omega_(0,1)(z))`:
  both the orientation and the denominator flip) and Hock's kernel after (2.1); EO's `dE_q(p)/omega(q)` ((4-1), (4-3), (4-9)) and
  ABDKS (1.1) with `omega^(0)_(1,0) = -y dx` are equal and opposite to it; `omega_(h,n) = (-1)^n W^(h)_n`, `F^TR_h = -F^(h)` (EO s.1
  also states `W_0^(g) = -F^(g)`). EO's (4-20) is `+Res BBB/(dx dy)`, so the note's `omega_(0,3)` is its negative. Anchors: on
  `x = -p^2/2` the convention gives `omega_(0,3) = dp^3/(p_0 p_1 p_2)^2` and `omega_(1,1) = (1/8) dp/p^4`; for `x = p + 1/p` the
  argument of `prop:mod-period` (pure part at `p = 0`, the polynomial part of `x` being irrelevant there) gives `F^TR_2 = -1/240`.
* `def:mod-curve`, `def:mod-pair`, `def:mod-xi`, `def:mod-ancvec`: facts checked (`x'' = 5 varsigma_b^2`, fibre polynomial, weights,
  `xi_(a,d)` orders, `<-log p, d xi_(a,d)> = xi_(a,d)(oo) - xi_(a,d)(0)`, the residue-theorem identity for log legs). `q_phi` is a
  finite sum for every leg of `def:mod-ancvec` (degree count at `oo`), and in any case the correlators vanish for `Sum d > 3h-3+n`.
* `lem:mod-laplace` (i)-(iv) on monomials (including `m = -1` in (i), the excluded `m = -2`); (iii): the sign of `w_j(v)(kp') =
  k^2 w_j(v')(p')` is fixed by the derivative at the branch point with one root for all `v`; the paragraph on the diagonal entries
  (`G_(-2) = -vartheta` gives the `1` of `R^-1 = 1 + O(vartheta^-1)`).
* `lem:mod-scale` (i) induction on `2h - 2 + n` (exponents of products, `omega_(0,2)` exponent 0), free energy through `Phi_ydx ->
  kappa mu Phi_ydx`; (ii) `m_k^* W_(k^5 v) = k^4 W_v`, naturality, rays to rays; the exponent `e`; `e` not an integer forces `C = 0`.
* `prop:mod-Rtr`. Step 1 (`xi^i`, IBP, holomorphy of `xi^i dx_kappa` at `varsigma_j`); Step 2 (`-kappa/(-5 kappa) = 1/5`); Step 3:
  (Mod2) `m z J_(m-1) - J_(m+3) + v J_(m-2) = 0`, (Mod1) via (iii) and (iv) with the elimination of `J_(m+4)`, giving `z d_v J_m =
  -J_(m-1) + J_m/varsigma_j`; the five relations `z d_v Y_b = bmod_b Y_(b+1)` one by one (`Y_2 = -v Ytilde_(-2)`); initial values
  `sqrt(-5 kappa) F_b(0) = eps_b zeta^(j(b+1))`; the colour rotation of the T-system (the closing `4 -> 0` by `zeta^(5j) = 1`);
  uniqueness by the pivot `4n zeta^(4j)` of `prop:set-modrows` (b) (valid over any field containing `Q(zeta)`); the table `m_b =
  (0,-1,-2,2,1)`. Step 4 (`varsigma_i^(1-m_b) varsigma^(m_b-1) = zeta^(-i(b+1))`, `(Psi_0^-1)_(ib) = eps_b zeta^(-i(b+1))/5`); Step 5
  (`R^-1 = M_0(-psi/(kappa varsigma^4))^T = M_0(psi/(kappa varsigma^4))^-1` by `lem:set-M0sympl`); (ii) both sides of (2-32) equal
  `Rhat^(m,(i))_0(t)/(sqrt(-5 kappa) varsigma)`; the unit `Sum bar e_k/varsigma_k`; the dilaton sign (REVIEW item 8).
* `thm:mod-tr`. Step 1: DNOPS s.2.6 recalls DOSS's theorem with (2-30)-(2-35) as cited; (2-33) holds for compact curves with
  meromorphic `dx` (DNOPS s.2.6) and, for the Hurwitz curves, by the sentence before Thm 8.1; (8-1) is `p^n + a_2 p^(n-2) + .. + a_n +
  a_(n+1)/(p - a_1)`, which at `n = 4`, `a = (0,0,0,0,4v)` is `4 W_v`; `f'' (p_c) = 20 p_c^2`, values `5 p_c^4`; Prop. 4.3 (REVIEW
  item 9). Step 3: `w^4 = 2 w^1`, `xi^(4) = xi^(1)/2`, `(Delta^4)^(1/2) = 2 (Delta^1)^(1/2)`, the TFT factor `2^(2h-2+n)` (vertices
  `2^(2h_v-2+n_v)`, dilaton leaves `2^-1`), the degree bookkeeping `c^-D` (legs `c^-k`, dilaton terms `c^-k`, edges `c^-(a+b+1)` with
  `E_(R(./c)) = c^-1 E_R(./c, ./c)`), total `4^(2-2h-n)`, and `lem:mod-scale` (i) at `(4,1)`. Step 4: `TFT^1(bar e_j^n) = (-5)^(h-1)
  varsigma_j^(2h-2+n)`, metric `-delta/5`, the scaling lemma for `lambda^(h-1) Omega` (splitting axioms, unit axiom, `R` with
  `lambda^(-1)` metric, `Sum_v (h_v - 1) + #E = h - 1`). Step 5: `H^(3-a)/5` dual to `H^a`, `nu_(3-a) = 2 - nu_a`, the partial
  fractions, `nu_a + a = (0,0,5,5,5)`, `eps_(3-a) = (1,-1,-1,1,1)`, and `xi'_a = xi_a` for all five `a`.
* `lem:mod-xfun`: the fibre polynomial, `Sum_i B(p^i(X), q') = dX dx(q')/(X - x(q'))^2`, the contour form of the recursion with discs
  disjoint from `J` (R-8), `Tr omega = 0` on an open set, the trace formula for residues at `X_b` (three simple points, distinct
  critical values), consequences (no diagonal poles; Fubini).
* `prop:mod-LD`: Steps 2-6; the Step 5 values `1, -1, 1, 3, 2` from `m = (4,0,1,2,3)`, `n = (0,-1,-2,2,1)`, `eps' = (-1,1,1,-1,-1)`,
  `alpha = (0,1,0,0,0)`; `I(p) = -p^5/5 + v log p`; `q_(Phi_ydx) = -z H^0`; the leg of `e` is `-p^5/5`.
* `cor:mod-trF`: `<Phi_ydx>^TR_h = (-5)^(h-1) <q_Phi>_h = -(-5)^(h-1) (2h - 2) F^m_h` with `lem:set-dilaton` (i) at `n = 0`.
* `lem:mod-closure` (i)-(iii) (exceptional terms, no `(log p)^2` because `G(x(p))'` has no `p^-1` term, `I(Q[v,x]) in Q[v,x]`);
  `lem:mod-genfn` (the telescoping `dB_n/dp = x^n phi_0'/n!`, `B` in `Q[v]`); `lem:mod-resid` (`[z^-n] I^C_m` against the explicit
  series of `def:gau-van`; `varrho = (-1,1,-1,-3,-2)`, the case `b = 2` through `Res W^n dW = 0`); `prop:mod-dii`; the remark on
  `e` (`Pi_0(z 1) + v Pi_0(H) = -v + v = 0`).
* `lem:mod-count`, `lem:mod-classL` (algebraic independence of `p, x`; the presentation `Q[v,x][T]/(T^5 - 4xT + 4v) = Q[x,T]`; the
  descent on the most negative exponent; (ii)-(iv), the count `d(w)` including `w = 0`), `lem:mod-modv` (bases of source and target
  mod `v`, the triangular matrix in weight `4n`), `thm:mod-lattice` (graded Nakayama for the cokernel; dimensions `d(w) + e(w)` on
  both sides; (b), (c), and the failure over `Q(v)`), `cor:mod-split`, `cor:mod-Ksplit` (`Pi_0` maps `V_w` into `Q v^(w/5)` by
  `def:tra-tower`), `cor:mod-basis` (`lem:tra-frames` (b), (c)).
* `lem:mod-abdks` against ABDKS Thm 7.1 as printed ((7.6), (7.7) with `1/hbar^2` and its series form, (7.9)); finiteness of the
  terms at fixed `hbar^(2h)`; the division argument by monic divisors. `thm:mod-Lreg`: residues at `oo`, the binomial expansion,
  iteration. `prop:mod-period`: Step 1 (`r = 2h - 1 + j`, `1 <= j <= h`, the `r = 0` term vanishes); Step 2 (degree `5 - 10h`); Step 3
  (pure part in `bar p`, corrections raise the order by `>= 3`); Step 4 (`C_(h,r) = ((-1)^j/j!) [x^h] a^j`, binomial resummation,
  `artanh`, Lagrange-Buermann for Laurent `H`, `u coth u`); check at `h = 1`: `Sigma_1 = -1/3`, and directly `Phi^A_1 = D_x(W''/(24 W'))`
  has the value `-1/(12 v)` at `p = 0`; Step 5. `lem:mod-vanper`, `thm:mod-gap`: `c_2 = 1/1200`, `c_3 = 1/25200`, `c_5 = 1/660000`.
* `lem:mod-fixedx`, `lem:mod-ibp` (the division remainder term is `O(dp/p^3)` at `oo`; `iota_X w` vanishes at `0` and `oo`).
  `thm:mod-var` against EO s.5: (5-1) `delta_Omega y|_x dx = -Omega`; the assumption "`Omega/dy` has no pole at branch points" stated
  with the Rauch formula (5-4); (5-13); Thm 5.1 = (5-14) for `g + k > 1`, (5-15) for `g >= 1`; the temperature variation (5-27)-(5-30)
  (`Omega = Int_alpha^alpha' B`, `Lambda = 1`; the note takes the opposite orientation with `Lambda = -1`); `D_Omega = delta_Omega`
  on a genus-0 curve; the sign bookkeeping `(-1)^n (-1)^(n+1) (-1)^(n+1) = 1`. The sign of (i) agrees with `lem:mod-genus0` (iii) at
  `(0,2)`, that of (ii) with `thm:mod-V` (ii) at `h >= 2` (derived without it).
* `lem:mod-genus0` (i)-(iv) (partial fractions `p/(p^5 - v) = Sum r_b/(p - varsigma_b)`, the Lie derivative identity, the residue at
  `t = 0`, `Int_r e_(b,1) = -1/varsigma_b`, the four pairings, `<(-log p)^3>_0 = 1/v`); `prop:mod-V1` (differentiation under the
  integrals, `d_v|_p = d_v|_x + Lie_X`, the sign of the `X`-term); `lem:mod-w11` (expansion of `K_c`, `omega_(0,2)(q, iota q)`, the
  coefficient of `chi^-1`, the Schwarzian value of `reg_00`; (ii) exactness, (iii) Fubini); `thm:mod-V`: (i) induction with the bases
  `(0,1)`, `(0,2)`; (ii) `h >= 2`, `h = 0` (`(F^TR_0)''' = 1/v`), `h = 1`: the expansion `x(1 + t) = 5/4 + (5/2)t^2 + (5/4)t^4 - t^5 +
  ..`, `chi = t + t^3/4 - t^4/5 + ..`, `y_3 = -1/4`, `{t; chi}(0) = -3/2`, `reg_00 = -1/4`, `eta_0 = -1`, `eta_2 = -1/4`, `Loc_1 =
  -(1/10)(12 - 1 - 3)/48 = -1/60`, rotation invariance of the coordinate-free form of `Loc`, homogeneity `v^-1`.
* `thm:mod-MLK` Steps 1-5 (`F^m_0 = (-5) F^TR_0 = -(5/2) v^2 log v` and `F^m_1 = F^TR_1` agree with `def:red-ML`), `cor:mod-onepoint`
  (exponent `1 - 2h + w/5` on both terms), `cor:mod-primary`, `rmk:mod-offK` (computed facts, labelled, unused: OUTLINE rule 3).

## 2. What remains open in s06

* Nothing in the chain of Theorem A. No `\draftnote` was added.
* Outside the chain (unchanged): the `\draftnote` of `rmk:mod-hock` (OUTLINE s.0.2, O2). The remark's evaluation at `p = 0` and the
  degree count at `oo` are sketched, as the remark now says; it is used in no proof.

## 3. Changes and findings that affect other sections (not edited by this review)

1. **s01 `sec:int-inputs` and s10 (closing paragraph of `thm:mai-main`), design OUTLINE s.0.2 / GS_2 s.3.1.** The (Q) list names
   EO Thm 4.7 (4-30) and EO Thm 5.1 (5-15), but no proof of the note uses either (REVIEW item 3 (b)); s06 now says so. The
   conditionality stated in s01 and s10 is therefore over-complete (true, but not "exact" in the sense of OUTLINE rule 8). Proposed:
   move these two items to the list "cited but not used" of `sec:int-inputs`, or keep them and accept the over-statement. Owners: s01,
   s10 and the design editor.
2. **BIBPLAN** (DuninBarkowski:2012bw): new pin `s.3.1, s.3.5` (the local Bergman kernel with its singular part and its Laplace
   transform), used in the paragraph after `lem:mod-laplace` (REVIEW item 1). BIBPLAN lists only s.2.3 and Thm 4.1 for DOSS.
3. **Appendix B** (`b_signs_2.tex`, row "dilaton equation"): consistent with s06. Optional precision for its owner: EO state Thm 4.7
   for `n >= 1`; the case `n = 0` is the definition (4-10). The row "`q_(phi^b_k) = -varrho_b psi^k H^b`" (l. 35) means the same as the
   new wording `-varrho_b z^k H^b` of `prop:mod-LD`.
4. **Labels and statements.** No label of s06 was added, removed or renamed, and no statement changed its content. Places in other
   sections that cite s06 (checked, none affected): s01 (`thm:mod-tr`, `prop:mod-Rtr`, `prop:mod-LD`, `thm:mod-lattice`,
   `thm:mod-Lreg`, `thm:mod-V`, `thm:mod-gap`, `thm:mod-MLK`, `lem:mod-xfun`, `rmk:mod-hock`), s02 (`thm:mod-gap`, `cor:mod-trF`), s04
   (`def:mod-weights`, `cor:mod-basis`), s05 and s10 (`thm:mod-MLK`, `thm:mod-tr`, `def:mod-weights`, `thm:mod-gap`), s07
   (`def:mod-weights`), appendix B (`prop:mod-Rtr`, `lem:mod-abdks`, `def:mod-curve`, `lem:mod-scale`, `thm:mod-tr`, `cor:mod-trF`,
   `prop:mod-LD`, `thm:mod-gap`, `conv:mod-kernel`, `lem:mod-laplace`, `def:mod-legs`, `def:mod-pair`, `thm:mod-var`), appendix C
   (`prop:mod-LD`, `thm:mod-MLK`, `prop:mod-Rtr`, `thm:mod-V`, `prop:mod-V1`, `prop:mod-dii`, `rmk:mod-hock`, `thm:mod-lattice`,
   `rmk:mod-offK`, `cor:mod-onepoint`, `def:mod-pair`, `prop:mod-period`, `thm:mod-tr`, `lem:mod-w11`). No other section refers to an
   item of an s06 statement, so the corrected item reference of `def:mod-curve` affects no one.
5. **NOTATION_2**: no change needed; s06 now uses the number `-5 kappa` (not `Lambda_kappa`) as NOTATION_2 prescribes.

## 4. Compilation

Scratch folder `allgenus_review_review_rf5/` (copies of `main.tex`, `preamble.tex`, `macros.tex`, `refs.bib`, `sections/`):
pdflatex, bibtex, pdflatex, pdflatex, each with exit code 0; no LaTeX error, no undefined reference or citation (as in the
baseline compile of the unmodified note, made first). The only remaining overfull boxes in s06 are two of less than 2 pt, present
before. (Other sections were being edited concurrently by their reviewers; the final compile used the state of all files at that
time.)

## 5. Files written

* `NOTE/review/rf5/REVIEW.md`, `REVIEW_2.md` (this report).
* Edited: `NOTE/sections/s06_model.tex`, `s06_model_2.tex`, `s06_model_3.tex`, `s06_model_4.tex`, `s06_model_5.tex`,
  `s06_model_6.tex`, `s06_model_7.tex`, `s06_model_8.tex`, `s06_model_9.tex` (`s06_model_10.tex` unchanged).
