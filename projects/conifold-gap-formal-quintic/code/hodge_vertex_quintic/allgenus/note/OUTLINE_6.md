# Outline, part 6: s06 the model (topological recursion on x = W_v, the R-matrix of the curve, the leg dictionary, the lattice,
# L-regularity, the Gap Theorem, Theorem V, the Model Lemma on K)

## s06_model — target 12 pp

Purpose: prove the Model Lemma on `K` in every genus, by identifying the model CohFT with topological recursion on the model curve
(the quoted theorem of DOSS, as recalled by DNOPS, at `f = 4 W_v`, transferred to `x = W_v`), translating insertions of `K` into legs
(leg dictionary, log coefficients, lattice theorem) and controlling their pairings (L-regularity, the Gap Theorem, Theorem V).

Covers: DC2 s.1 D1-D11, s.2.1-s.2.3, s.3 (parts 4-5), s.4 C1-C13, s.7 G1-G8, W1-W11, L1-L5; DC1 s.2 MOD, TR, s.3 MOD/TR, s.4 CV2-CV5,
CV8-CV12, CV17, s.7 O-10..O-13, O-17, O-19; ST2 s.1 (M1-M10), s.2 step 7; ST2b K1-K3, K10-K11, O1-O13, O26-O27; TRM s.0-s.4.3; TRRF;
XY s.0-s.3; XYR; DP; DR; LL; LL2; LR; MSG1 s.1, s.3-s.5; MSG2 s.8-s.9; MRF items 1-5, 9-10; G4 (all pieces) with G4R, G4R_2 (R-1..R-8),
G4_5 s.6 (F1-F10); GS s.2.6, s.2.8, s.2.10.
First sentence of the section: topological recursion is applied here, and only here, to the curve `x = W_v(p)`, `y = p` (the quoted
theorem is applied to `f = 4 W_v`, the same curve up to scaling, and transferred by `lem:mod-scale`; [ABDKS25] Thm 7.1 is applied to
`x = W_v` itself); nothing about the quintic enters (ST2 s.1). Statements that the sources prove for families of curves are written for
the model curve only (rule 8); no family version is printed. Subsections `sec:mod-tr`, `sec:mod-legs`, `sec:mod-lattice`, `sec:mod-Lreg`,
`sec:mod-V`, `sec:mod-MLK`.

### sec:mod-tr — conventions; the model as topological recursion
* `conv:mod-kernel` (convention; DC2 C1-C2; DC1 CV9; ST2b K1; LR s.2; G4 s.1): the engine kernel `(1/2) Int_(iota_b q)^q omega_(0,2)(p_0, .)/
  ((y(q) - y(iota_b q)) dx(q))`, `iota_b` the local involution of `x` at `varsigma_b` (REV_2 M13 (d): not `sigma`), `omega_(0,1) = y dx`,
  `omega_(0,2) = dp_1 dp_2/(p_1 - p_2)^2`, equal to [Dunin-Barkowski:2015caa] (2-26)-(2-29) and [Hock25] (2.1); the printed kernel of
  [Eynard:2007kz] (4-3), (4-9) EQUALS that of [ABDKS25] (1.1) with `omega^(0)_(1,0) = -y dx`, and both are the opposite (G4_5 s.6 F4):
  `omega_(h,n) = (-1)^n W^EO_(h,n) = (-1)^n omega^ABDKS_(h,n)`, `F_h = -F^print_h`, `Lambda^print = -1` for the temperature variation; this
  `omega_(0,3)` is minus that of [Eynard:2007kz] Thm 4.1 (4-20); `(2 - 2h) F_h = Sum_b Res_(varsigma_b) Phi_ydx omega_(h,1)`, `d Phi_ydx = y dx`.
  DNOPS's (2-31), (2-32) are read as formal Laplace series (`lem:mod-laplace`; F1). Stated once here and in appendix B. Sign anchors
  (classical, not used in proofs): Airy intersection numbers; `F_2 = -1/240` for `x = p + 1/p`.
* `def:mod-curve` (defn; DC2 D1): `x = W_v(p)`, `y = p`, `v != 0`; `dx = (p^5 - v) dp/p^2`; branch points `varsigma_b = zeta^-b varsigma`
  (REV R8 (d)), `x''(varsigma_b) = 5 varsigma_b^2`, critical values `x(varsigma_b) = (5/4) zeta^b varsigma^4`, pairwise distinct, `dy != 0`;
  weights `p:1, v:5, x:4, pbar := v/p : 4`; homogeneity `p -> kp, v -> k^5 v`: `omega_(h,n) -> k^(5(2-2h-n)) omega_(h,n)` (TRRF s.3.4; also by
  `lem:mod-scale`, G4_2 s.6). Other normalisations only in appendix B.
* `def:mod-pair` (defn; DC2 D3; G4_3 s.1 "Precision", F2): legs in the leg class `LL = Q[v][p] + Q(v)[x]` are paired by the residue pairing
  `<phi, w> := Sum_b Res_(varsigma_b) phi w`; for legs in `Q[v][p]` it equals the pairing on a circle `|p| = R > |varsigma|`, for legs with a pole
  at `p = 0` it does not (`Res_(p=0) x^5 omega_(1,1) = -15/4`); log legs `-log p` are paired on a ray `frak r` from 0 to oo avoiding the branch
  points (independent of the ray); mixed pairings.
* `def:mod-xi` (defn; DC2 D4): `xi_a = eps'_a v^(alpha_a) p^(m_a)/(p^5 - v)`, `(m_a, eps'_a, alpha_a) = (4,-1,0), (0,1,1), (1,1,0), (2,-1,0),
  (3,-1,0)`; `D_x = d/dx = p^2/(p^5 - v) d/dp`; `xi_(a,d) := D_x^d xi_a`. The signs `eps'_a` are not the row signs `eps_b` (rename locally).
* `def:mod-ancvec` (defn; DC2 D5): `q_phi(z) := Sum_(a,d) <phi, d xi_(a,d)> H^a z^d`.
* `thm:mod-tr` (thm; M1, DC2 thm:trm-ident): for `v != 0`, `2h - 2 + n > 0`: `omega_(h,n) = (-5)^(h-1) Sum <prod_i psi^(d_i) H^(a_i)>_h prod_i d
  xi_(a_i, d_i)`, the brackets being the ancestor correlators of `Omega^mod_v`; equivalently `<q_(phi_1) .. q_(phi_n)>_h = (-5)^(1-h) <phi_1 ..
  phi_n>^TR_h`. Proof: (1) the quoted theorem of OUTLINE s.0.2 ((Q): DOSS Thm 4.1 via DNOPS s.2.6 (2-30)-(2-35), Prop. 4.3 (1) and the
  sentence before Thm 8.1) at `f = 4 W_v`, i.e. DNOPS (8-1) with `n = 4`, `a = (0,0,0,0,4v)` (REV R3): critical points `p^5 = v`, `f''(p_c) =
  20 p_c^2 != 0`, critical values `5 p_c^4` pairwise distinct; the compatibility holds (`dy/dx` has poles only at the branch points, no residue
  at oo; G4R s.1.1); `p_Dub = p` since `a_1 = 0`; (2) `prop:mod-Rtr` at `kappa = 4`; (3) the TFT dictionary `Delta_j = -5 varsigma_j^2` from (2-30)
  with the uniform orientation `dw_j/dp(varsigma_j) = sqrt(-5) varsigma_j` (F10), the frozen basis `\bar e_j := e^TR_j/sqrt(-5)` with the pairing
  rescaled, and the factor `(-5)^(h-1)` (GS s.2.10; DOSS's vertex weight, G4R s.1.3 R-3); (4) the transfer `kappa = 4 -> 1` by `lem:mod-scale`
  and the bookkeeping `[R(./c).TFT]_(h,n) = Sum_D c^-D [R.TFT]^(D)_(h,n)` (G4_2 s.6 with R-2; never "push-forwards preserve degree"); (5) the
  identification with `Omega^mod_v` (`def:set-modcohft`; TRM s.2.4). Corrections TRRF s.6 ("any point where `x` is simply ramified"); author
  names -> citations (O-19). Not used (G4_2 s.7; F7): the second half of DNOPS Thm 8.1, DNOPS Prop. 5.2, the uniqueness of homogeneous
  R-matrices, LG s.3's canonical-frame recursion, TRM's sign-of-`kappa` argument.
* `prop:mod-Rtr` (prop; TR-2; DEC D6): statement G4_2 s.1 [G] with G4R s.1.3 R-1: for `kappa != 0`, `v != 0`, DNOPS's coordinates `x_kappa =
  x_kappa(varsigma_j) - w_j^2/2` with `dw_j/dp(varsigma_j) = sqrt(-5 kappa) varsigma_j` (one root for all `j`), and (2-31), (2-32) read by
  `lem:mod-laplace`: (a) `R^-1(vartheta^-1)^j_i = M_0(t)_ij`, `t = -1/(kappa varsigma^4 vartheta)` (`vartheta` the Laplace variable, DNOPS's
  `zeta`, G4's `lam`), so by `lem:set-M0sympl` `R_TR(w) = M_0(w/(kappa varsigma^4))` in the basis `e^TR_j` of (2-30); (b) (2-32) holds; in the
  frozen basis the unit is `1hat = Sum_k \bar e_k/varsigma_k`, and `omega_(0,1) = y dx` produces the dilaton leaves `psi(1hat -
  M_0^-1(psi/(kappa varsigma^4)) 1hat)`, the translation of [PPZ15] Def. 2.13 (DNOPS s.2.3 prints the opposite sign: R-1; GS s.2.6). Proof
  G4_2 s.3-s.4 [G]. Delete `\draftnote{D6/G2 ..}` and the "uniqueness of homogeneous R-matrices" and "sign of `kappa`" text.
* `lem:mod-laplace` (lem; new; G4_2 s.2 [G]): the formal Laplace functional `L_c[omega] := Sum_m a_m G_m(vartheta)` for 1-forms with zero
  residue at a simple zero `c` of `dx_kappa`, Gaussian moments `G_(2n) = (2n-1)!! vartheta^-n`, `G_(-2n) = (-1)^n vartheta^n/(2n-1)!!`, `G_m = 0`
  for odd `m != -1`; its properties (integration by parts, coordinate sign, homogeneity, derivative in `vartheta`), numbered (i)-(iv) (G4's
  (L1)-(L4) are not printed as tags); the finite part at `i = j` gives the `delta_ij` of `R^-1 = I + O(vartheta^-1)`. The "Meaning" sentence
  (Watson's lemma) is motivation, not part of the definition (R-7).
* `lem:mod-scale` (lem; new; G4_2 s.6 [G]): for constants `kappa, mu != 0` the TR data of `(P^1, kappa x, mu y, omega_(0,2))` are `(kappa
  mu)^(2-2h-n) omega_(h,n)` (`2h - 2 + n > 0`) and `(kappa mu)^(2-2h) F_h` (`h >= 2`); induction on `2h - 2 + n`.
* `lem:mod-xfun` (lem; M6): statement and proof G4_3 s.1-s.2 [G]: for `2h - 2 + n > 0`, `n >= 1`, `p_2..p_n` off the branch points and `f`
  holomorphic near the branch points (`f` in `Q(v)[x]` in the application), `Res_(varsigma_b) f(x) omega_(h,n)(., p_2, ..) = 0` for each `b`:
  the trace of `omega_(h,n)` over the five sheets of `x` vanishes (the trace of the kernel is zero), then the trace formula for residues
  (G4's "Lemma T"; rename locally) with the discs `U_b` disjoint from `{p_2, .., p_n}` (R-8) and the distinct critical values; so the leg
  class acts modulo `Q(v)[x]` for the residue pairing. Optional pins [Eynard:2007kz] Thm 4.4 (4-23), (A-25) (G4_3 s.3). Delete `\draftnote{G3}`.

### sec:mod-legs — the leg dictionary, the residue identity
* `def:mod-legs` (defn; DC2 D8): `phi^b_0 = p, log p, v/p, p^3, p^2` (`b = 0..4`); `phi^b_(k+1) = I(phi^b_k)`, `I(F) := -Int F dx` termwise in
  `p`, constant 0 (DR s.2); `D_x phi_(k+1) = -phi_k`; `wt(phi^b_k) = w_b + 4k`, `w_b = nu_b + 1 = (1,0,4,3,2)`. Source `Phi^b_k` -> `phi^b_k`.
* `prop:mod-LD` (prop; LD): `<phi^b_k, d xi_(a,d)> = -varrho_b delta_ab delta_dk`, `varrho = (-1, 1, -1, -3, -2)`; hence `q_(phi^b_k) = -varrho_b
  psi^k H^b`; the TR leg of `v^m psi^k H^b` is `-varrho_b^-1 v^m phi^b_k`; unit leg `y = p`; dilaton leg `Int y dx = -phi^0_1` with `q = -z 1`;
  log leg `-log p` with `q = H`; `e` has leg `-p^5/5`. Proof DC2 s.3 LD steps 1-7 [W] (ST2 s.2 step 7.2); step 6 [W-]: print the table
  `(kappa_b, n_b) = (1,0), (1,-1), (-v,-2), (3,2), (2,1)` (rename `kappa_b` locally) and `Sum_c c^N = 5 v^(N/5) [5 | N]`.
* `cor:mod-trF` (cor): `F^m_h = (-5)^(1-h) F^TR_h` (`h >= 2`), by the dilaton equation on both sides (`lem:set-dilaton`; [Eynard:2007kz] Thm
  4.7 (4-30) through `conv:mod-kernel`) and LD at `(b,k) = (0,1)` [W] ST2 M1. `lem:mod-dil`: dropped (D8).
* `lem:mod-genfn` (Lemma G), `lem:mod-resid` (the residue identity `R^b_n = varrho_b Pi_0(psi^n H^b)`, `R^b_n := Res_(p=0) W_v^n d phi^b_0/n!`),
  `prop:mod-dii` (the log coefficients of the legs, source "D(ii)"): the `log p`-coefficient of `phi^b_k` is `varrho_b Sum_i Pi_i(psi^k H^b)
  (-x)^i/i!`; `Lambda(X)` is log-free iff `X` lies in `K_0` (REV R8 (d): not `s`). Proof DC2 s.3 steps 1-8 [W] (DR s.3, proof of record; DP
  s.1); inputs `lem:gau-vanmod`, `def:tra-tower`, `lem:mod-closure` (state `lem:mod-closure` first). Corrections O6, O26, DR Correction 4(ii).

### sec:mod-lattice — the lattice theorem (source "D(iii)")
* `def:mod-weights` (defn; DC2 D9-D10): `V_w` (source `S_w`); `FF = (+)_w V_w`; the leg map `Lambda(v^m psi^k H^b) := varrho_b^-1 v^m phi^b_k`;
  `Q[p, pbar]`, the leg class `LL`; `MM := Q[p, pbar] (+) log p . Q[v,x]`, `NN := MM/Q[v,x]` (DP's sans `M`, `N`; not the ring `\Nring` of s02);
  `d(w) = floor(w/5) + 1 - [5 | w]`, `e(w) = #{(m,n) : 5m + 4n = w}`.
* `lem:mod-count` (DP Lemma 0) [W DP s.2]; `lem:mod-classL` (DP Lemma 1 (a)-(d)) [W DP s.3; DR Correction 6]; `lem:mod-closure` (DP Lemma 2
  (i)-(iii)) [W DP s.4; DR Correction 6]; `lem:mod-modv` (DP Lemma 3) [W DP s.5].
* `thm:mod-lattice` (thm; D(iii)): for every `w >= 0`: (a) `V_w` has `d(w)` triples with `b != 1` and `e(w)` with `b = 1`; (b) for `X` in `K_0 cap
  V_w`, `Lambda(X)` lies in `Q[p, pbar] = Q[v][p] + Q[v][x]`, exactly; (c) `Lambda : K_0 cap V_w -> (LL/Q(v)[x])_w` is bijective, both of
  dimension `d(w)`; (d) `Lambda : V_w -> NN_w` is bijective, `NN` free over `Q[v]`. Scope: `Q[v]` coefficients (over `Q(v)` it fails: `H^2/v` has
  leg `1/p`). Proof DC2 s.3 steps 1-4 [W] (DP s.1-s.5; DR s.2, s.6; graded Nakayama). Corrections O1/DR Corr. 1, O2/DR Corr. 2, K11.
* `cor:mod-split` (DP Cor. 4) [W DP s.5; proved without D(ii), DP s.6; G3R_2 s.11, K9]; `cor:mod-Ksplit` (DP Cor. 5) [W DP s.9]: on `K` the
  only legs beyond `Q[p, pbar]` are multiples of `log p`.
* `cor:mod-basis` (cor; DC2 prop:trm-coeff): the `v^m frame(k,a)` of weight `w` form a basis of `K_0 cap V_w` and `dim K_0 cap V_w = d(w)`: from
  `lem:tra-frames` in the graded setting and the count of `thm:mod-lattice` (a) alone (G3_2 s.2 Corollary (i); G3 s.3 item 1); ST2 M5 [W-c]
  as a second reference. `rmk:mod-h2leg`: dropped (D8).

### sec:mod-Lreg — L-regularity and the Gap Theorem
* `lem:mod-abdks` (lem; TR-7): G4_4 s.1-s.2 [G] with R-4: [ABDKS25] Thm 7.1 applies to `x = W_v`, `y = p` itself (`W_v` rational, `dx` with five
  simple zeros, `y` regular with `dy != 0`; G4 F6); with `omega^engine = (-1)^n omega^ABDKS`: `omega_(h,n) = N_(h,n) prod dp_i/prod (p_i^5 - v)^K`,
  `N_(h,n)` in `Q[v][p_1, .., p_n]` (divisibility in `Q[v][p]`, the `p_i`, `p_i - p_j` being monic); the edge weight (7.7) with its factor
  `1/hbar^2`; the optional bound `K = 6h - 4 + 2n` is dropped. XY's `x = z^4 + V/z` only in appendix B (factor `4^(2h-2+n)`). Correction XYR
  s.2 (Remark 7.3 concerns the discriminant; `V = 0` is a boundary point).
* `thm:mod-Lreg` (thm; M7, "L-regularity"): G4_4 s.3 [G]: for `2h - 2 + n > 0` every residue pairing of `omega_(h,n)` with legs in `Q[v][p]` (all
  `n >= 1` slots) lies in `Q[v]`; legs in `Q[v][x]` pair to zero (`lem:mod-xfun`). The "11 forced zeros" are not evidence (XYR s.3). Delete
  `\draftnote{W5}`.
* `prop:mod-period` (prop; TR-8; proof of record of the Gap Theorem, DEC Q4): G4_4 s.4 [G]: for `h >= 1`: (a) `omega_(h,1) = d Phi_h`, `Phi_h`
  rational, holomorphic at `p = 0` and `p = oo`; (b) `Int_oo^0 omega_(h,1) = (B_2h/(2h)) v^(1-2h)` (any path avoiding the branch points); (c) for
  `h >= 2`, `F^TR_h = B_2h v^(2-2h)/(2h(2h-2))`. Proof Steps 1-5 from [ABDKS25] (7.10): exactness, `p = oo`, the monomial count at `p = 0`, the
  exact identity `\Sigma_h = -4^h B_2h/(2h)` (binomial series and Lagrange-Buermann inversion; it replaces the "Borel identity", F5), the
  contour lemma; (c) uses `thm:mod-Lreg` for `<p^5>_h`. Local letters for G4's `S_h` (`\Sigma_h`), `c_(h,r)`, `epsl`, `t`, `phit`.
* `lem:mod-vanper` (lem): `<-log p>_h = Int_(frak r) omega_(h,1) = Sum_b Res_(varsigma_b) (-log p) omega_(h,1) = -(B_2h/(2h)) v^(1-2h)` (G4_4 s.4
  Step 5 [G]; orientation XYR s.2).
* `thm:mod-gap` (thm; "Gap Theorem", TR-4, M8; never "Theorem G"): on the model curve, for `h >= 2`, `F^TR_h = B_2h v^(2-2h)/(2h(2h-2))` (engine
  convention), hence `F^m_h = c_h v^(2-2h)`, `c_h = (-5)^(1-h) B_2h/(2h(2h-2))` (`cor:mod-trF`). Proof of record: `prop:mod-period` (c).
* `rmk:mod-hock` (rmk; new; DEC Q4): the first proof, by [Hock25] Cor. 1.2, (1.3) applied to `x = W_v` (b-points `0` and `oo`; the convention
  of Cor. 1.2 stated, LR F5, O13), TRM s.4.1 steps 2-5, TRRF s.3.2 (Lagrange inversion `A_g = -B_2g/(2g)`); not used.
  `\draftnote{the closed forms printed in [Hock25, s.3.2] were not re-checked; only Cor. 1.2 is cited}`. `rmk:mod-residue`: dropped.

### sec:mod-V — Theorem V (log legs)
* `lem:mod-fixedx` (lem) [W LL s.2.1; LR s.3]: `d_v p|_x = -p/(p^5 - v)`, `d_v|_x = d_v|_p - Lie_\cX`, `\cX := p/(p^5 - v) d/dp`, `d_v(y dx)|_x =
  -dp/p = Int_(frak r) omega_(0,2)(., q)`. `lem:mod-ibp` (lem) [W LL s.2.2; LR s.3]: (a) L-slot `oint phi Lie_\cX w = -<(\cX dphi)_+, w>`; (b) log
  slot `Int_(frak r) Lie_\cX w = 0`.
* `thm:mod-var` (quoted, in the form used; REV_2 M15): in the two ranges of [Eynard:2007kz] Thm 5.1: (5-14) `d_v|_x omega_(h,n)(J) =
  Int_(frak r) omega_(h,n+1)(J, q)` for `n >= 1`, `h + n > 1`; (5-15) `d_v F_h = Int_(frak r) omega_(h,1)` for `h >= 1` (used for `h >= 2`); with the
  temperature variation of s.5.3 (5-27)-(5-30), through `conv:mod-kernel` (`Lambda^print = -1`; LR s.2). Hypotheses as checked (LL s.2.3; LR
  s.2): simple branch points with `dy != 0`, variation at fixed `x`, `Omega/dy` regular at the branch points (EO (5-4)), the path joins the two
  poles of `y dx` avoiding the branch points (EO (5-13), footnote 8), convergence. The case `(0,1)` is the hypothesis, `(0,2)` is
  `lem:mod-genus0`. Corrections F1/O7 (s.5.3, not s.5.1); W9 (`Omega = -d_v(y dx)|_x = +dp/p` by EO (5-1)).
* `lem:mod-genus0` (lem) [W LL s.2.4; LR s.3]: `d_v|_x omega_(0,2) = Sum_b r_b e_(b,1)(p_1) e_(b,1)(p_2)`, `omega_(0,3) = -Sum_b e_(b,1)^3/x''(varsigma_b)`
  (step 3 [W-]), `Int omega_(0,3) = d_v|_x omega_(0,2)`; the four explicit (0,3) pairings; `<(-log p)^3>_0 = 1/v`. Correction F1/O13 (this
  `omega_(0,3)` is minus EO Thm 4.1 (4-20)). `prop:mod-V1` (prop; (V1)) [W LL2 s.3; LR s.3]: the recursion with `+` before the `\cX`-term.
* `lem:mod-w11` (lem; new; G4_5 s.2 [G] with R-6): the local formula for `omega_(1,1)` at a point off the branch points, as a sum of local
  operators applied to `omega_(0,2)(p_0, .)` (G4's `C_a`: rename locally), and the reason why `omega_(1,1)` has zero residues and is
  holomorphic at `p = 0` (R-6).
* `thm:mod-V` (thm; Theorem V): (i) pairings with `a >= 1` L-legs in `Q[v][p]` and `k` log legs lie in `Q[v]`; (ii) `<(-log p)^k>_h =
  (F^TR_h)^((k))(v)`, `F^TR_1 = -(1/12) log v`, `F^TR_0 = (1/2) v^2 log v` (`k >= 3`); (iii) (V1). Proof DC2 s.3 steps 1-7: (i) induction on `k`
  [W LL2 s.3], bases `(h,a) = (0,1), (0,2)` by the explicit (0,3) pairings (W7); (ii) base `h >= 2` by LR F2 without the variational formula
  [W]; base `h = 0` [W] `lem:mod-genus0`; base `h = 1` G4_5 s.2-s.3 [G]: `lem:mod-w11` and the finite evaluation (`-1/60` per branch point at
  `v = 1`, five equal contributions by rotation, total `-1/12`; homogeneity gives `<-log p>_1 = -1/(12 v)`); the two cross-checks of G4_5 s.4
  are not printed. Delete `\draftnote{D5/G1 ..}`. Inputs: [Eynard:2007kz] Thm 5.1, [ABDKS25] Thm 7.1, `thm:mod-gap` (O11).
  `rmk:mod-dilroute`: dropped (D8).

### sec:mod-MLK — the Model Lemma on K
* `thm:mod-MLK` (thm): ML_K(h) (`def:red-ML`) holds for every `h >= 0`. Proof DC2 s.3 steps 1-6 [W] (ST2 s.2 steps 7.1-7.6; LL2 s.3 Corollary
  with LR F3): `thm:mod-tr`, `prop:mod-LD`, `prop:mod-dii`, `cor:mod-Ksplit` with `thm:mod-lattice` (b) and `lem:mod-classL` (a), `lem:mod-xfun`,
  `thm:mod-Lreg`, `thm:mod-V`, `thm:mod-gap`; insertions of Theorem R are `Q[[u]]`-combinations of elements of `K cap V_w` (`lem:tra-coeff`),
  applied termwise (`lem:tra-conv` (b)). Corrections LR F3/O9 (closed by `thm:mod-lattice` (b)), O5, O27, DR Correction 5, K11.
* `cor:mod-onepoint` (cor; MOD-1 on `K`): for `X` in `K` homogeneous with `<X>_h` singular, `<X>_h = Pi_0(X) (F^m_h)'(v)`. [W] from the theorem.
* `cor:mod-primary` (cor; MOD-5): for `h >= 1` every primary correlator vanishes unless all insertions are `H`, and `<H^n>_h = (F^m_h)^((n))`.
  [W] MSG2 s.9(c); MRF item 9.
* `rmk:mod-offK` (rmk; kept: necessity of `K`): off `K` the one-point law fails and no genus-independent non-linear `V` exists (MOD-1, MOD-2);
  the counterexamples (`<psi^5 H>_3 = 53/1209600`, the Hessians) are computed: appendix C only.
* Dropped (D8): `rmk:mod-toda`, `rmk:mod-H2`; MOD-6 (first-order response `Delta_g(5)`) only in appendix C, no label.

Length 12 pp.
