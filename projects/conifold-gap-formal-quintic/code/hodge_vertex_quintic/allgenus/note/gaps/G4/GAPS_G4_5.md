# GAPS_G4, part 5: item 4, base `h = 1` of `thm:mod-V` (D5, DIGEST item G1) -- PROVED; findings for the editor

Conventions: GAPS_G4.md s.1. Sources searched: LL2 s.3 (proof of (ii): "`h = 1`: by homogeneity `<-log p>_1 = c/v` with one rational
number `c`, computed exactly over `Q(zeta5)`: `c = -1/12`"), LR s.3 (reproduced mod two primes), XY s.2.3 (vanishing period; `F_1` via the
variational formula "with its sign fixed by consistency"), XYR s.2 (`F_1` "weakened"), LG s.5.4 (genus-one formula on the CohFT side).
No source writes the number by hand. Below: a finite exact derivation from the definition of the recursion (no quoted theorem), and two
independent cross-checks.

## 1. Statement

**(base `h = 1` of `thm:mod-V` (ii)).** On `x = W_v`, `y = p`, `v != 0`: `<-log p>_1 := Int_(frak r) omega_(1,1) = -1/(12 v)`, i.e.
`<-log p>_1 = (F^TR_1)'(v)` with `F^TR_1 = -(1/12) log v` (`F^TR_1` is defined up to a constant by this; with (V1) at `a = 0` it gives
`<(-log p)^k>_1 = (F^TR_1)^((k))(v)` for all `k >= 1`).

## 2. A local formula (`lem:mod-w11`, PROVED from the definition; new lemma for `sec:mod-V`, or inline in the proof)

Let `a` be a branch point, `t := p - a`, and `zeta_a` a local coordinate with `x = x(a) + alpha zeta_a^2` (`alpha != 0` a constant; then
`sigma_a : zeta_a -> -zeta_a`). Write `y = Sum_k y_k zeta_a^k`, `B = [ (zeta_1 - zeta_2)^-2 + beta(zeta_1, zeta_2) ] dzeta_1 dzeta_2` near `(a, a)`
(`beta` holomorphic and symmetric), and `beta_00 := beta(0,0) = (1/6) {t; zeta_a}(0)`, `{t; zeta} = t'''/t' - (3/2)(t''/t')^2` (the regular
part of `t'(z_1) t'(z_2)/(t(z_1) - t(z_2))^2` on the diagonal). For a 1-form `eta = Sum_m eta_m zeta_a^m dzeta_a` holomorphic at `a` put
  `C_a[eta] := -(1/(4 alpha y_1)) [ beta_00 eta_0 + eta_2/12 - (y_3/(4 y_1)) eta_0 ]`.
**Lemma.** (i) `omega_(1,1)(p_0) = Sum_a C_a[B(p_0, .)]` for `p_0` off the branch points (`eta = B(p_0, .)` as a form in the second slot).
(ii) `<-log p>_1 = Int_(frak r) omega_(1,1) = Sum_a C_a[-dq/q]`.
*Proof.* (i) The recursion (`conv:mod-kernel`) gives `omega_(1,1)(p_0) = Sum_a Res_(q -> a) K(p_0, q) omega_(0,2)(q, sigma_a q)`. Expand
`B(p_0, xi) = dp_0 Sum_m B_m(p_0) zeta^m dzeta` at `xi` near `a`. Then `Int_(sigma q)^q B(p_0, .) = 2 dp_0 Sum_(m even) B_m zeta^(m+1)/(m+1)`,
`y(q) - y(sigma q) = 2 Sum_(k odd) y_k zeta^k`, `dx = 2 alpha zeta dzeta`, so
`K = (dp_0/dzeta) (4 alpha y_1)^-1 [B_0 zeta^-1 + (B_2/3 - B_0 y_3/y_1) zeta + O(zeta^3)]`, and `omega_(0,2)(q, sigma q) = B(zeta, -zeta) =
-[(4 zeta^2)^-1 + beta_00 + O(zeta^2)] dzeta^2` (`dzeta_2 = d(-zeta)`; `beta(zeta, -zeta) = beta_00 + O(zeta^2)` by symmetry). The `zeta^-1`
coefficient of the product is `-(4 alpha y_1)^-1 [B_0 beta_00 + (B_2/3 - B_0 y_3/y_1)/4]` times `dp_0`, i.e. `C_a[B(p_0, .)]`.
(ii) Each residue is a contour integral over a small circle around `a`; `frak r` avoids the circles and `omega_(1,1)` is `O(dp/p^2)` at `oo`
and holomorphic at `0`, so `Int_(frak r)` may be taken inside: `Int_(p_0 in frak r) B(p_0, q) = dq [-1/(p_0 - q)]_(p_0 = 0)^(oo) = -dq/q` for
`q` off the ray; `C_a` is linear and continuous in the coefficients. QED
*Check of the convention.* On the Airy curve `x = zeta^2`, `y = zeta` (`alpha = y_1 = 1`, `y_3 = beta = 0`, `B_2(p_0) = 3/p_0^4`) the lemma gives
`omega_(1,1) = -(1/16) dp/p^4`, and `lem:mod-scale` with `kappa = -1/2` gives `(1/8) dp/p^4` on `x = -p^2/2`, which is
`<tau_1>_1 (3!!) = 3/24` (Witten-Kontsevich), the sign anchor of `conv:mod-kernel` (TRM s.3 regression).

## 3. The evaluation (finite, exact)

Take `v = 1`, `a = varsigma_0 = 1`. Then `x(1 + t) = 5/4 + (5/2) t^2 + 0 t^3 + (5/4) t^4 - t^5 + O(t^6)`, so `alpha = 5/2` and
`zeta = t (1 + t^2/2 - (2/5) t^3 + O(t^4))^(1/2) = t + t^3/4 - t^4/5 + O(t^5)`, `t = zeta - zeta^3/4 + zeta^4/5 + O(zeta^5)`.
* `y = 1 + t`: `y_1 = 1`, `y_2 = 0`, `y_3 = -1/4`.
* `t'(0) = 1`, `t''(0) = 0`, `t'''(0) = -3/2`: `{t; zeta}(0) = -3/2`, `beta_00 = -1/4`.
* `-dq/q = -t'(zeta) dzeta/(1 + t(zeta)) = -(1 - (3/4) zeta^2 + ..)(1 - zeta + zeta^2 + ..) dzeta = -(1 - zeta + zeta^2/4 + O(zeta^3)) dzeta`:
  `eta_0 = -1`, `eta_2 = -1/4`.
* `C_1[-dq/q] = -(1/10) [ (-1/4)(-1) + (-1/4)/12 - (-1/4)(-1)/4 ] = -(1/10) (12/48 - 1/48 - 3/48) = -(1/10)(1/6) = -1/60`.
*Five equal contributions.* `rho_b(p) := zeta^-b p` maps `1` to `varsigma_b` and satisfies `x o rho_b = zeta^b x`, `y o rho_b = zeta^-b y`,
`rho_b^* B = B`, `rho_b^*(dq/q) = dq/q`; so `rho_b` carries the local data at `1` to those at `varsigma_b` (`sigma_b o rho_b = rho_b o sigma_0`;
the kernel `K` is invariant: `(y(q) - y(sigma q)) dx(q)` picks up `zeta^-b zeta^b = 1`), whence `C_(varsigma_b)[-dq/q] = C_1[-dq/q] = -1/60`.
Total at `v = 1`: **`<-log p>_1 = 5 (-1/60) = -1/12`**.
*All `v`.* With `p = varsigma p'` (`varsigma^5 = v`, any root): `x_v(varsigma p') = varsigma^4 x_1(p')`, `y = varsigma y'`, so the pull-back of
`omega^(v)_(1,1)` is `varsigma^(5(2-2-1)) omega^(1)_(1,1) = v^-1 omega^(1)_(1,1)` (`lem:mod-scale` with the `p`-scaling; TRRF s.3.4), and the
ray `frak r` pulls back to a ray; the pairing does not depend on the ray (`omega_(1,1)` has zero residues). Hence `<-log p>_1 = -1/(12 v)`.
QED. (Inputs: `conv:mod-kernel`, `def:mod-curve`, `def:mod-pair`, `lem:mod-scale`; nothing quoted.)

## 4. Two independent cross-checks (written; not needed for the proof)

* **Vanishing period** (`prop:mod-period` at `h = 1`, G4_4 s.4; mod ABDKS Thm 7.1): `<-log p>_1 = -(B_2/2) v^-1 = -1/(12 v)`; in XY's
  normalisation `<-log z>_1 = -1/(12 V)` and the transfer factor `4^(2h-2+n) = 4` with `V = 4v` give the same. This is the route the design
  (D5) named; it is a second proof.
* **Closed form.** `omega_(1,1) = -(5/8) p^4 (p^10 + 18 v p^5 + 6 v^2)/(p^5 - v)^4 dp` (ABDKS (7.10) via G4_4 s.3 (iii); it also equals
  `Sum_a C_a[B(p_0, .)]`, checked exactly at 8 points, `g4_w11.log` (B)). With `u = p^5` the ray integral is `-(1/8) Int_0^oo (u^2 + 18 v u +
  6 v^2)/(u - v)^4 du` (`p^4 dp = du/5`; the image path avoids `u = v`); in `w = u - v` the integrand is `w^-2 + 20 v w^-3 + 25 v^2 w^-4`
  (no residue, so the path does not matter), with primitive `G(w) = -1/w - 10 v/w^2 - (25/3) v^2/w^3`; `G(oo) = 0`, `G(-v) = (1 - 10 +
  25/3)/v = -(2/3)/v`, so the integral from `w = -v` to `oo` is `(2/3)/v` and `<-log p>_1 = -(1/8)(2/3)/v = -1/(12 v)`.
* (Not used: the CohFT side, Givental's genus-one formula, LG s.5.4, gives `dF_1 = -(1/12) du/u`; it needs the shifted-family reading.)

## 5. Evidence (VERIFIED; `g4_w11.py` -> `g4_w11.log`, exact)

* (A) `C_(varsigma_b)[-dq/q] = -1/60` for each `b` at `v = 1` (exact over `Q(zeta5)`, local data computed by series reversion), total `-1/12`;
  controls: opposite ray orientation `+1/12`; `beta_00` dropped `1/24`; `y_3` dropped `-11/96`.
* (B) `Sum_a C_a[B(p_0, .)]` = the closed form at `p_0 = 3, -2, 1/2, 7/3`, `v = 1` and `v = 32` (8 exact equalities).
* (D) the closed-form integral `= -1/(12 v)` at `v = 1, 32, -3`.

## 6. Findings for the editor (all clusters may read; nothing here changes a statement of the sources)

* **F1** (item 1): state once (in `conv:mod-kernel` or appendix A) that DNOPS's integrals (2-31), (2-32) are formal Laplace series
  (`lem:mod-laplace`), with the finite part at `i = j` (equivalently DOSS's "`delta_ij` + regular part"); otherwise `R^-1 = I + O(z)` is not
  even defined for the diagonal entries.
* **F2** (item 2): `def:mod-pair`, `thm:mod-Lreg`, `lem:mod-xfun`: the L-leg pairing for `LL = Q[v][p] + Q(v)[x]` is the residue pairing
  `Sum_b Res_(varsigma_b)`; the circle `|p| = R` equals it only for legs holomorphic at `p = 0`.
* **F3** (item 1, cross-cluster): `lem:set-M0sympl` (G4_2 s.5) proves the symplecticity of `M_0` from `prop:set-modrows` alone; it closes
  the R7 entry "`def:set-M0` symplecticity" and belongs in s02.
* **F4** (source erratum, for `design/SOURCE_ERRATA.md`): XYR s.2 says the ABDKS kernel is "minus the kernel of arXiv:math-ph/0702045";
  in fact ABDKS (1.1) with `omega_(1,0) = -y dx` EQUALS the printed EO kernel (4-3), (4-9) (both `(1/2) Int_z^(sigma z) B/((y(z) - y(sigma z))
  dx)`), and both are minus the engine kernel. XYR's conclusion `omega^code = (-1)^n omega^A` is right.
* **F5** (item 3): route 2 needs no Borel summation; print Step 4 of G4_4 s.4 (the identity `S_h = 4^h A_h = -4^h B_2h/(2h)` by the binomial
  series and Lagrange inversion); the Borel integral, if kept, is a remark.
* **F6** (item 3, W5): apply ABDKS Thm 7.1 to `x = W_v` itself; XY's `x = z^4 + V/z` only in appendix B with the factor `4^(2h-2+n)`.
* **F7** (item 1): the following drop out of the chain: the uniqueness of homogeneous R-matrices, DNOPS Prop. 5.2, TRM's sign-of-`kappa`
  argument, LG s.3's canonical-frame recursion, and the second half of DNOPS Thm 8.1. The `\draftnote{D6/G2 ..}` of OUTLINE_6 can go.
* **F8** (item 4): the `\draftnote{D5/G1 ..}` can go: base `h = 1` is the finite derivation of s.3 (inputs listed there).
* **F9** (item 2): pin EO by the PDF numbers: Thm 4.2; Thm 4.4 = (4-23), (4-24); Thm 4.5 = (4-25); Cor. 4.1 = (4-26)-(4-29), proof (A-25).
* **F10** (item 1): the uniform orientation `dw_j/dp = sqrt(Lambda) varsigma_j` is the one matching `Psi_0` and the TFT dictionary
  `Delta_j = Lambda varsigma_j^2` of `thm:mod-tr` step 3; state it where (2-30) is used.
* **Nothing in G4 is left OPEN.** Inputs that remain quoted: DNOPS Thm 8.1 (+ s.2.3, s.2.6, proof of Thm 5.1), ABDKS Thm 7.1, EO Thm 4.2
  (optionally Thm 4.4); used from s02 of the design: `prop:set-modrows` (b) (its derivation is written in the sources, DB S18 [W]).
