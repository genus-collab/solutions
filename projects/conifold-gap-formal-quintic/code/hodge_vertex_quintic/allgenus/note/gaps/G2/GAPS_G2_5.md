# GAPS_G2 part 5: `cor:gau-colours` (item 5, G6): `D_i(z) = D_0(zeta^-i z)`

Continuation of `GAPS_G2_4.md`. Setting and local notation: `GAPS_G2.md` s.2; `Rot_i`, `T_i` as in s.6.4.

## 9. The colours of `D`

### 9.1 Statement (proposed final form)
Let `D = diag(D_0, .., D_4)` be defined by `thm:gau-R0` (iv): `QQ D = P_can QQ^m` for the amplitude matrices (rows `a`, columns =
colours), the same exponential `e^((u_i - u_i^c)/z)` being stripped from `P_can` times the model's colour-`i` WKB solution at `v(u)` and
from `e^(-u_i^c/z)` times the quintic's colour-`i` WKB solution (PS-C8). Here `QQ^(i)` is the colour-`i` column of the actual rows
(`def:set-rows`, QRR-normalised, conifold expansion) and `QQ^(m,i)` that of the model rows `varsigma^(nu_b - 1) Rhat^mod,(i)_b(z/varsigma^4)`
at `v = v(u)` (`def:set-M0`). Then:
* `D_i(z) = D_0(zeta^-i z)` for every `i` in Z/5; hence `D_0 = 1` implies `D = 1`.
* `D_0` lies in `Q[[z]]`; `D_i` lies in `Q(zeta)[[z]]`, and for `i != 0` it lies in `Q[[z]]` iff `D_0` lies in `Q[[z^5]]`.
(`QQ` stands for the calligraphic `\Qr` of NOTATION table 2; `QQ^m` for `\Qr^{\mathrm m}`.)

### 9.2 The three-line computation, with its three inputs
Ring: `R_zeta := Q(zeta)((s))[[z]]` (`s = u^(1/5)`), an integral domain, which contains the amplitudes of every colour and
`P_can = Sum z^k P_k(u)` (`u = s^5`); for colour 0 it contains the comparison ring `Q[L, 1/L][[z]]` of `thm:gau-R0` (iv). On
`R_zeta^5`, `Rot_i Y(u, z) = T_i Y(u, zeta^-i z)` is additive and satisfies `Rot_i(Y d(z)) = Rot_i(Y) d(zeta^-i z)` for scalar `d` in
`Q(zeta)[[z]]`, and `Rot_i(P Y) = (Rot_i P)(Rot_i Y)` (6.4 (ii)).
* (C1) `QQ^(m,i) = Rot_i QQ^(m,0)` (model columns; `lem:gau-colmod`, 9.3, PROVED).
* (C2) `QQ^(i) = Rot_i QQ^(0)` (actual columns; `lem:gau-colrows`, 9.4, PROVED from the inputs listed there).
* (C3) `Rot_i P_can = P_can` (6.4 (iii), PROVED: `P_can` is graded).
Then, for each `i`:
    `QQ^(i) D_i(z) = P_can QQ^(m,i) = P_can Rot_i(QQ^(m,0)) = Rot_i(P_can QQ^(m,0)) = Rot_i(QQ^(0) D_0(z)) = Rot_i(QQ^(0)) D_0(zeta^-i z)
     = QQ^(i) D_0(zeta^-i z)`,
using R0 for colour `i`, (C1), (C3) with `Rot_i(P Y) = (Rot_i P)(Rot_i Y)`, R0 for colour 0, the scalar rule, and (C2). The row-0
entry of `QQ^(i)` has `z^0`-term `zeta^i L/varpi_0 != 0` (`ph_0 = 1`, `def:set-rows`, conifold frame), so it is not a zero divisor in
the domain `R_zeta`, and `D_i(z) = D_0(zeta^-i z)`. If `D_0 = 1` then every `D_i = 1`. QED
*Exponents.* The computation is on amplitudes only. That the stripped exponentials match colour by colour is what makes (C1), (C2)
statements about WKB solutions of the right colour: `Rot_i` sends `e^(U_0/z)` to `e^(U_0/(zeta^-i z)) = e^(U_i/z)` (`U_i = -(5/4)
zeta^i varsigma^4 = zeta^i U_0`) and `e^((u_0 - u_0^c)/z)` to `e^((u_i - u_i^c)/z)` (`u_i - u_i^c = -(5/4) zeta^i u^(4/5)/tau`,
`lem:set-canon`). *Fields.* Colour 0: the actual colour-0 amplitudes, `P_can` and `Fs_b` are rational, so `D_0` is in `Q[[z]]`
(already in `thm:gau-R0` (iv)). Then `[z^n] D_i = zeta^(-in) [z^n] D_0`, which gives the second bullet of 9.1.
*Why gradedness is essential.* For the non-graded regular gauge `P' = P_can c(z)` the series defined by `QQ D' = P' QQ^m` is
`D'_i = D_i c(z)`, and `D'_i(z) = D'_0(zeta^-i z)` holds for all `i` iff `c(z) = c(zeta^-i z)`, i.e. iff `c` is in `k[[z^5]]`: the
colour relation characterises the graded gauges among the regular ones.

### 9.3 `lem:gau-colmod` (new; proposed label): the colour covariance of the model columns (C1)
**Statement.** For each `i`, the T-system of `prop:set-modrows` (b) with `calT^(i)_k := -zeta^i + (t/5)(k - 4t d/dt)`,
`F_1 = calT_-1 F_0, F_2 = -calT_-2 F_1, F_3 = calT_2 F_2, F_4 = -calT_1 F_3, F_0 = -calT_0 F_4`, `F_0(0) = 1`, has the solution
`F^(i)_b(t) = zeta^(ib) Fs_b(zeta^-i t)` (`Fs_b` the colour-0 solution). Hence `Rhat^mod,(i)_b(t) = zeta^(i(b+1)) Fs_b(zeta^-i t) =
zeta^i F^(i)_b(t)` (`def:set-M0`) are the rows of a WKB solution of the model system with exponential `e^(U_i/z)`, and at `v = v(u)`:
`QQ^(m,i)_b(u, z) = zeta^(i(b+1)) QQ^(m,0)_b(u, zeta^-i z)`, i.e. (C1).
**Proof.** (1) `t d/dt` commutes with the substitution `t -> zeta^-i t`, and `t = zeta^i (zeta^-i t)`; hence for every series `G`:
`calT^(i)_k [G(zeta^-i t)] = zeta^i [calT^(0)_k G](zeta^-i t)`. (2) Put `G_b(t) := zeta^(ib) Fs_b(zeta^-i t)`. By (1):
`G_1 = zeta^i [calT^(0)_-1 Fs_0](zeta^-i t) = calT^(i)_-1 G_0`; `G_2 = -zeta^(2i) [calT^(0)_-2 Fs_1](zeta^-i t) = -calT^(i)_-2 G_1`; in the
same way `G_3 = calT^(i)_2 G_2` and `G_4 = -calT^(i)_1 G_3`; finally `G_0 = Fs_0(zeta^-i t) = -[calT^(0)_0 Fs_4](zeta^-i t) =
-zeta^(-i) calT^(i)_0 [Fs_4(zeta^-i t)] = -zeta^(-5i) calT^(i)_0 G_4 = -calT^(i)_0 G_4`, and `G_0(0) = 1`. By the uniqueness clause of
`prop:set-modrows` (b) (the pivot at `t^(n+1)` is `4n zeta^(4i) != 0`), `F^(i) = G`. (3) By `prop:set-modrows` (b),
`e^(U_i/z) varsigma^(nu_b - 1) F^(i)_b(z/varsigma^4)` solves the model system; so does its multiple by the constant `zeta^i`. As
`varsigma` depends on `u` only, `QQ^(m,i)_b = varsigma^(nu_b - 1) zeta^(i(b+1)) Fs_b(zeta^-i z/varsigma^4) = zeta^(i(b+1)) QQ^(m,0)_b(u,
zeta^-i z)`. QED
*Remark.* Equivalently: the model system is a weighted cyclic shift with weights depending on `v` only, so `Rot_i` maps its solutions
to solutions (6.4 (ii)), and it sends the colour-0 WKB solution to a colour-`i` one with the same normalisation at `t = 0` up to the
factor `zeta^i` that `def:set-M0` builds in. UR s.2 and EV s.2.2 state the formula `zeta^(i(b+1)) F_b(zeta^-i t)`; the covariance
computation was not written. *Evidence:* `g2_colours.log` (all five colours to `t^10` over `Q(zeta)`; control: the normalised rows
`zeta^(i(b+1)) Fs_b(zeta^-i t)` are NOT the T-system solution, they are `zeta^i` times it).

### 9.4 `lem:gau-colrows` (new; proposed label): the colour covariance of the actual columns (C2)
**Inputs (all from s02; none is proved in this cluster).**
* (i) *Cyclic shape* (part of `prop:set-cyclicLR`, R1, cluster G1): the actual flat solution of colour `i`, `Y^(i) = e^(u_i/z) X^(i)`
  (rows `X^(i)_a`), satisfies `z D Y_a = kappa'_a Y_(a+1)` (`a` in Z/5) with coefficients `kappa'_a` (at large radius `C_1 kappa^LR_a`)
  that are invertible and independent of `z` and of `i`. Only this shape is used, not the values of the `kappa'_a`.
* (ii) *Canonical coordinates* (`def:set-genus0`): `D u_i = zeta^i L`. Only `D u_i` enters: every statement below is about amplitudes,
  i.e. about solutions of the conjugated equations (`D` replaced by `D + (D u_i)/z`), so the integration constants of the `u_i` play no
  role (changing `u_i` by a constant multiplies the solution, not the amplitude, by `e^(c/z)`).
* (iii) *Row 0 and (QRR)* (`thm:set-mirror` (M), cluster G1, and `thm:set-qrr` in its colour-`i` form): `X^(i)_0 = psi_i/C_0`, where
  `e^(u_i/z) z psi_i` is a formal WKB solution of PF5, `psi_i = L + O(z)`, and `psi_i(q = 0, z) = h_(w_i)(z)` (at `q = 0`:
  `L = C_0 = 1`; `P^(i)_a(q = 0) = zeta^(ia) h_(w_i)(z)`, DA D6).
* (iv) `lem:gau-wkb` (formal WKB uniqueness, colour 0) and `lem:set-hrot` (`h_(w_i)(z) = h_(w_0)(zeta^-i z)`).
**Statement.** Under (i)-(iv): `psi_i(q, z) = psi_0(q, zeta^-i z)`, `X^(i)_a(q, z) = zeta^(ia) X^(0)_a(q, zeta^-i z)`, and for the
normalised rows `QQ^(i)_a = zeta^i X^(i)_a` (`def:set-rows`: `QQ = P/Delta^(1/2)`, `Delta_i^(1/2) = zeta^(4i) C_0/L`, `P^(i)_a = (C_0/L)
X^(i)_a`): `QQ^(i)_a(z) = zeta^(i(a+1)) QQ^(0)_a(zeta^-i z)`, i.e. `QQ^(i)_(a,n) = zeta^(i(a+1-n)) QQ^(0)_(a,n)`, an identity of
coefficients in the differential ring (tensored with `Q(zeta)`); applying the conifold expansion gives (C2).
**Proof.** *Step 1 (PF5 is invariant under `z -> zeta z`).* In `PF5 = (zD)^5 - 1 - q prod_(k=1..5)(5zD + kz)` every monomial that
contains `z` has degree exactly 5 in `z`; so the substitution `z -> zeta^-i z` multiplies it by `zeta^(-5i) = 1`. The conjugated
operator `PF5^(u) := e^(-u/z) PF5 e^(u/z)` is `PF5` with `zD` replaced by `zD + Du`; with `zeta^-i z D + L = zeta^-i (zD + zeta^i L)` and
(ii), the substitution `z -> zeta^-i z` turns `PF5^(u_0)` into `PF5^(u_i)`. Hence `z psi_0(q, zeta^-i z)` (= `zeta^i` times the
substituted amplitude `zeta^-i z psi_0(zeta^-i z)`) is the amplitude of a colour-`i` WKB solution of PF5, with leading term `z L`.
*Step 2 (uniqueness in colour `i`).* `lem:gau-wkb` is stated for colour 0; the substitution `z -> zeta^i z` maps colour-`i` WKB
amplitudes of PF5 to colour-0 ones (Step 1 with `i` replaced by `-i`) and constant series to constant series. So
`psi_0(q, zeta^-i z) = psi_i(q, z) d(z)` with `d` in `1 + z Q(zeta)[[z]]`.
*Step 3 ((QRR) fixes `d`).* At `q = 0`: `h_(w_0)(zeta^-i z) = h_(w_i)(z) d(z)`. By `lem:set-hrot` the left side is `h_(w_i)(z)`, a unit
of `Q(zeta)[[z]]`; so `d = 1`.
*Step 4 (rows).* By (i), `X^(i)_(a+1) = (1/kappa'_a)(zD + zeta^i L) X^(i)_a` (`z D (e^(u_i/z) X) = e^(u_i/z)(zD + Du_i) X`). Induction on
`a`: `X^(i)_0 = psi_i/C_0 = X^(0)_0(zeta^-i z)`; if `X^(i)_a(z) = zeta^(ia) X^(0)_a(zeta^-i z)` then `X^(i)_(a+1)(z) = zeta^(ia)
(1/kappa'_a)(zD + zeta^i L)[X^(0)_a(zeta^-i z)] = zeta^(i(a+1)) [(1/kappa'_a)((zeta^-i z) D + L) X^(0)_a](zeta^-i z) =
zeta^(i(a+1)) X^(0)_(a+1)(zeta^-i z)`, as `kappa'_a` does not depend on `z`.
*Step 5 (normalisation).* `QQ^(i)_a = zeta^(-4i) X^(i)_a = zeta^i X^(i)_a` and `QQ^(0) = X^(0)`, so `QQ^(i)_a(z) = zeta^(i(a+1))
QQ^(0)_a(zeta^-i z)`. *Step 6 (conifold expansion).* The coefficient identity `QQ^(i)_(a,n) = zeta^(i(a+1-n)) QQ^(0)_(a,n)` is preserved
by any ring map acting on coefficients (the conifold expansion of `def:set-confexp`, `lem:gau-confexp`) and by hatting. This step
needs the coefficients of the actual rows to lie in the differential ring on which the conifold expansion is defined: for row 0 this
is `thm:gau-R0` (iii) (colour 0; colour `i` then by Steps 1-3), for rows `a >= 1` it follows from row 0 by the forward relations of
(i) (this is R1's proposed `cor:gau-rowsNN`, R1 item 3). QED
*Status.* PROVED given (i)-(iv). (iv) is proved in the note (`lem:gau-wkb` [T] in another cluster; `lem:set-hrot` one line); (iii)
contains the quoted inputs (M) and (QRR) (the passage from `R(q=0) = diag(h_(w_i))` to `psi_i(q=0) = h_(w_i)` is `thm:set-qrr`'s
`\draftnote` (a), the same for every colour); (i) is the open item R1 of cluster G1 (only its shape part is needed); (ii) is in the
sources. The covariance itself was VERIFIED numerically in the sources (FQ2 s.1; CS s.1); it was nowhere derived.
*Evidence here:* `g2_colours.log` checks the multiset identity `w_i = zeta^i w_0` behind `lem:set-hrot` (control: `{zeta^i - zeta^j}
U {-5}` is not covariant). No new computation of the actual WKB columns was made (s.10).
