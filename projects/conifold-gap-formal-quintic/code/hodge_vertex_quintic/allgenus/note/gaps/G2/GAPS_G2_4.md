# GAPS_G2 part 4: graded gauges and `prop:gau-graded` (G5, item 2), `lem:gau-diag` step 4 (G8), `thm:gau-G` Cor. (a) (G9)

Continuation of `GAPS_G2_3.md`. Setting and local notation: `GAPS_G2.md` s.2.

## 6. Graded gauges

### 6.1 `lem:gau-degree` (new; proposed label): the Z/5-degree of the gauge equation
For a series `P = Sum_k z^k P_k` (`P_k` in `gl_5(k((u)))`) and `delta` in Z/5 put `P^[delta] := Sum_k z^k Sum_((a,b): k+b-a = delta mod 5)
(P_k)_(ab) E_(ab)`, so `P = Sum_delta P^[delta]` (direct), and `P` is **graded** (`P_k = Sum_b p_(k,b) E_(b+k,b)` for all `k`, as in
`def:gau-onepole`) iff `P = P^[0]`. Let `GEop(P) := z P' - A P + P A~`. Then `GEop(P^[delta]) = GEop(P)^[delta+1]`. In particular,
if `P` solves `GE(A, A~)`, so does every `P^[delta]`.
*Proof.* Give `z^k f(u) E_(ab)` the degree `k + b - a` mod 5. `d/du` does not change it and the factor `z` adds 1, so `z P'` raises the
degree of each homogeneous term by 1. `A (f z^k E_(cd)) = b_(c-1) f z^k E_(c-1,d)` has degree `k + d - c + 1`, and
`(f z^k E_(cd)) A~ = b~_d f z^k E_(c,d+1)` has degree `k + d + 1 - c`. So `GEop` maps homogeneous series of degree `delta` to
homogeneous series of degree `delta + 1`, and by additivity `GEop(P)^[delta+1] = GEop(P^[delta])`. QED

### 6.2 Existence of graded gauges and the graded freedom (items 1g, 2; G5)
* (i) *Existence.* Let `P^nat` be as in 5.1 (a). Its `z^0`-term `frak d Theta_U` is diagonal, hence of degree 0. By 6.1, `P^gr :=
  (P^nat)^[0]` solves `GE(A, A~)`, and `P^gr_0 = P^nat_0` is invertible: `P^gr` is a graded regular gauge, over `Q` if `P^nat` is.
  (Equivalently, run 5.3 in graded form: take `c_k` in `k[[u]] calY^j`, `j = k mod 5`, and the shift-`(-k)` lift; the free constant
  `gamma` then occurs only at `k = 0 mod 5`, the "resonances only at `k = 0 mod 5`" of UR s.3 Rmk (d).)
* (ii) *Freedom.* Let `P` be a graded regular gauge. By 5.1 (b) the regular gauges are `P c(z)`, `c` in `k[[z]]^x`. Write
  `c = Sum_(r in Z/5) c^(r)`, `c^(r) := Sum_(n = r mod 5) c_n z^n`. Since `P` has degree 0 and `z^n` adds `n` to the degree,
  `P c^(r)` is homogeneous of degree `r`, so `(P c)^[r] = P c^(r)`. Hence `P c` is graded iff `P c^(r) = 0` for `r != 0`, iff `c^(r) = 0`
  for `r != 0` (injectivity in 5.1 (b)), iff `c` lies in `k[[z^5]]`. So the graded regular gauges are exactly `P c(z^5)`, `c` in
  `k[[z^5]]^x`: this is the graded clause of Theorem U.
* (iii) *The normalising entry.* For graded `P` the entry `(P_k)_(11)` vanishes unless `1 = 1 + k mod 5`, so
  `P_11(0, z) = Sum_(q>=0) (P_(5q))_(11)(0) z^(5q)` lies in `k[[z^5]]`; its constant term `(P_0)_(11)(0) = c(0) frak d(0) d^U_1(0)` is
  non-zero (5.1 (b)), so `P_11(0, z)` is a unit of `k[[z^5]]` (rational if `P` is). This is `prop:gau-graded`'s last clause and GI
  s.3 Cor. (a)'s first sentence.

### 6.3 The canonical gauge (`def:gau-can`)
For a graded regular gauge `P` put `P_can := P P_11(0, z)^-1`. By 6.2 (iii) and (ii) it is a graded regular gauge, and
`(P_can)_11(0, z) = 1` because the factor is a scalar. It does not depend on `P` (two graded gauges differ by a unit `c(z^5)`, which the
normalisation removes). It is the unique regular gauge, graded or not, with `P_11(0, z) = 1`: if `P'` is a regular gauge with
`P'_11(0, z) = 1`, then `P' = P_can c(z)` by 5.1 (b), so `c(z) = P'_11(0, z) = 1`. Over `Q` its coefficients are rational (5.5).
*Evidence* (s.10, `g2_gauge.log`, `g2_G.log`): the construction of 5.3 in graded form, frames `(0,0)` and `(1/3,-2/7)`, reproduces
the canonical values recorded in GI s.5 and VT `v2_run.log`.

### 6.4 Colour covariance: equivalent forms of gradedness (item 2: "`S A S^-1 = zeta A`"; item 5: "colour covariance of `P_can`")
Let `k` contain `Q(zeta)`, `S := diag(zeta^-a)`, `T_i := diag(zeta^(i(a+1))) = zeta^i S^-i`, and `(Rot_i P)(u, z) := T_i P(u, zeta^-i z)
T_i^-1 = S^-i P(u, zeta^-i z) S^i` (matrices), `(Rot_i Y)(u, z) := T_i Y(u, zeta^-i z)` (vectors).
* (i) The entry `(a, b)` of `[z^k] Rot_i P` is `zeta^(i(a+1) - ik - i(b+1)) (P_k)_(ab) = zeta^(-i(k+b-a)) (P_k)_(ab)`, so
  `Rot_i P = Sum_delta zeta^(-i delta) P^[delta]`. Hence `P` is graded iff `Rot_1 P = P` iff `Rot_i P = P` for all `i`. UR's formulation
  ("covariant under `z -> zeta z` combined with `diag(zeta^-a)`") is `Rot_(-1) P(u, z) = S P(u, zeta z) S^-1 = P(u, z)`: the same
  condition. `(1/5) Sum_i Rot_i P = P^[0]` is the averaging of `prop:gau-graded` step 1; 6.1 makes the averaging (and `zeta`)
  unnecessary.
* (ii) `S A S^-1 = zeta A`, `S A~ S^-1 = zeta A~` (`(S A S^-1)_(a,a+1) = zeta^(-a) b_a zeta^(a+1)`); equivalently `A T_i = zeta^i T_i A`.
  Hence `Rot_i` maps solutions of `z d_u Y = A Y` to solutions: `z d_u (Rot_i Y) = zeta^i T_i [(zeta^-i z) d_u Y](u, zeta^-i z) =
  zeta^i T_i A Y(u, zeta^-i z) = A Rot_i Y`; the same for `A~`. Also `Rot_i (P Y) = (Rot_i P)(Rot_i Y)`.
* (iii) (**colour covariance of `P_can`**) `P_can` is graded (6.3), so `Rot_i P_can = P_can`, i.e.
  `T_i P_can(u, zeta^-i z) T_i^-1 = P_can(u, z)` and `Rot_i(P_can Y) = P_can Rot_i(Y)` for every vector series `Y` for which the products
  are defined. (The non-graded regular gauges `P_can c(z)` are covariant iff `c` is in `k[[z^5]]`.)

## 7. `lem:gau-diag` step 4 (item 3, G8)

**Statement.** Let `P` be a graded regular gauge between one-pole systems with (H), (EP), `P_0 = diag(p_(0,a))`. Then
`prod_a p_(0,a) = c(0)^5 det Theta_U(0) = p_(0,1)(0)^5 prod_a rho^U_a`, `rho^U_a := d^U_a(0)/d^U_1(0) = prod_(c=1..a-1) (b~_c/b_c)(0)`
(read `(b~_1/b_1)(0) = r~/r`), a constant. Under the programme normalisation (N0) of s.2 (`b_a(0) = b~_a(0)` for `a != 1`, `r = r~`):
every `(b~_a/b_a)(0) = 1`, so `p_(0,a)(0) = p_(0,1)(0)` for all `a`; for `P_can`, `p_(0,1)(0) = [z^0] (P_can)_11(0, z) = 1`, hence
`p_(0,a)(0) = 1` for every `a` and `prod_a p_(0,a) = 1`.

**Proof (GI's route, with the missing line).** GI s.5 shows that `prod_a p_(0,a)` is constant in `u` (diagonal of order `z^1`, [W]).
Its value is `prod_a p_(0,a)(0)`. Order `z^0` gives `b_a p_(0,a+1) = p_(0,a) b~_a`; at `u = 0` this is `p_(0,a+1)(0) = p_(0,a)(0)
(b~_a/b_a)(0)` for `a != 1`, and after multiplication by `u` it is `r p_(0,2)(0) = r~ p_(0,1)(0)` for `a = 1`. Going round the cycle
from `a = 1` gives `p_(0,a)(0) = p_(0,1)(0) rho^U_a`, which is the general formula; under (N0) every factor is 1. For `P_can` the
`z^0` coefficient of `(P_can)_11(0, z) = 1` is `p_(0,1)(0) = 1`. QED
**Second proof (from `lem:gau-lead`; also replaces GI s.5 steps 1-5).** By 5.1 (b), `P_0 = gamma frak d Theta_U` with `gamma = c(0)`.
Hence `prod_a p_(0,a) = det P_0 = gamma^5 frak d^5 det Theta_U = gamma^5 det Theta_U(0)`, constant, since `frak d^5 = det Theta_U(0)/
det Theta_U`. For `P_can`: `1 = p_(0,1)(0) = gamma d^U_1(0) = gamma (b~_0/b_0)(0)`, so `gamma = 1` under (N0), and then
`prod_a p_(0,a) = det Theta_U(0) = prod_a d^U_a(0) = 1`. Moreover `p_(0,0) = frak d = (det Theta_U)^(-1/5)` (under (N0)), and
`det Theta_U = prod_(a=0..4) prod_(c<a) (b~_c/b_c) = prod_(c=0..3) (b~_c/b_c)^(4-c)`: so `p_(0,0)^5 = prod_(c<4) (b_c/b~_c)^(4-c)`, GR s.2's
formula (S27 step 5), from which `p_(0,0) = tau^(-1/4)/varpi_0` follows as in GR s.2 (with `v'^5 = v/(u (1-u)^5)`). QED
*Scope.* The value 1 uses (N0), i.e. `t_c'(0) = varpi_0(0) = 1`, `v = u + O(u^2)` and the scale `varpi_1 = u + O(u^2)` (programme
family). Rescaling `varpi_0` by a constant keeps the product 1 (the ratios at 0 become `(w, 1, w, 1/w, 1/w)`, product of the
`rho^U_a` = 1); rescaling `varpi_1` changes `r` and the product. *Evidence:* s.10 (`g2_gauge.log`: `p_(0,a)(0) = 1`, product 1 to
`u^30`, both frames; control `b~_0 -> 2 b~_0`, `b~_1 -> b~_1/2` gives `1/16 = prod rho^U_a`, as the formula predicts).

## 8. `thm:gau-G` Cor. (a) (item 4, G9)

**Statement.** `P_can` is a graded regular gauge from `A~` (pulled-back model) to `A` (quintic), and it is the unique regular gauge
`P` with `P(u, z) Y^C(v(u), z) = Y^C_Q(u, z)` in `Q((z))[[u]]^5`.
**Proof.** Graded: 6.3. Existence: (G) (`thm:gau-G`, GI s.3 [W]) for `P_can` gives `P_can Y^C = (P_can)_11(0, z) Y^C_Q = Y^C_Q`.
Uniqueness: let `P` be a regular gauge with `P Y^C = Y^C_Q`. By 5.1 (b), `P = P_can c(z)` with `c` in `Q[[z]]^x` (both are regular
gauges over `Q`; 5.5). Multiplication by the scalar `c(z)` commutes with everything, and all products lie in `Q((z))[[u]]^5` (at each
order in `u` only finitely many terms meet: `P_k` has `u`-coefficients, `Y^C(v(u), z)` lies in `Q[z^-1][[u]]^5` by `lem:gau-vanshape`).
So `Y^C_Q = P Y^C = c(z) P_can Y^C = c(z) Y^C_Q`, i.e. `(c(z) - 1) Y^C_Q = 0`. The `u^0`-coefficient of `Y^C_Q` is `-e_1`
(`lem:gau-vanq` (c)), so `(c(z) - 1)(-e_1) = 0` in `Q((z))^5` and `c = 1`. QED
*Variant without (G) for `P_can`:* (G) for `P` itself gives `P Y^C = P_11(0, z) Y^C_Q`; comparing `u^0`-terms with `Y^C_Q` gives
`P_11(0, z) = 1`, and 6.3 gives `P = P_can`. Either way the uniqueness clause (5.1 (b)) is the only input beyond (G).
*Related (Cor. (b), [W*] in GI):* `P Y^C` has no positive power of `z` iff `P_11(0, z)` is a constant: if it is the constant `c_0`, then
`P Y^C = c_0 Y^C_Q` lies in `Q[z^-1][[u]]^5`; conversely the `u^0`-term of `P Y^C` is `-P_11(0, z) e_1` by (R1). Hence the regular gauges
with this property are exactly `Q^x P_can`.
*Evidence* (s.10, `g2_G.log`): (G) for the constructed `P_can` in all slots `z^-13..z^3`, mod `u^10`, both frames; for
`P_can (1 + z^5/7)`: `P Y^C = (1 + z^5/7) Y^C_Q`, and the `[z^5 u^0]` slot is `-1/7` (the failing control of the uniqueness).
