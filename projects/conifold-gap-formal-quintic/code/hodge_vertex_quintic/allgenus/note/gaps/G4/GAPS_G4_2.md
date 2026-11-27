# GAPS_G4, part 2: item 1, `prop:mod-Rtr` (D6, DIGEST item G2) -- PROVED

Conventions: GAPS_G4.md s.1. Here `kappa != 0` is a rational number (the cases used are `kappa = 4`, DNOPS's `f = 4 W_v` at `a =
(0,0,0,0,4v)` of R3, and `kappa = 1`, the note's curve), `x_kappa := kappa W_v`, `Lambda := -5 kappa` (a number), `sqrt(Lambda)` one fixed
square root, `lam` the Laplace variable of DNOPS (2-31) (their `zeta`). Sources searched first: TRM s.2.3 (z^8 check, uniqueness
"standard, not re-proved", sign of `kappa` "ARGUED"), TRRF claim 2 and correction 4 (z^10 check; "a direct route also exists: xit^i x' =
-(1/c_i) Sum_k c_i^k p^(2-k)"), LG s.3-s.4, UR Prop. S. No source writes an all-orders proof; the route below is TRRF's "direct route",
carried out. It does not use the uniqueness of homogeneous R-matrices, DNOPS Prop. 5.2, or TRM's sign-of-`kappa` argument.

## 1. Statement (proposed wording for the note)

**Proposition (`prop:mod-Rtr`).** Let `kappa != 0`, `v != 0`, and choose the local coordinates of DNOPS s.2.6 at the branch points
of `x_kappa` by `x_kappa = x_kappa(varsigma_j) - w_j^2/2`, `dw_j/dp (varsigma_j) = sqrt(Lambda) varsigma_j` (the same root for all `j`).
Read (2-31) and (2-32) as formal Laplace series (`lem:mod-laplace` below). Then, as formal series in `lam^-1`, for all `i, j`:
  (a) `R^-1(lam^-1)^j_i = M_0(t)_(ij)`, `t = -1/(kappa varsigma^4 lam)`; with `w := lam^-1` and `M_0(-t)^T = M_0(t)^-1` (`lem:set-M0sympl`):
      **`R_TR(w) = M_0(w/(kappa varsigma^4))`** in the basis `e_j` of (2-30); for `kappa = 1`, `v = 1`: `R_TR = M_0`;
  (b) (2-32) holds; both sides equal `Lambda^(-1/2) varsigma^-1 Rhat^mod_0^(i)(t)`. In the frozen basis `et_j := e_j/sqrt(Lambda)` the unit
      `Sum_k Delta_k^(-1/2) e_k` is `1hat = Sum_k et_k/varsigma_k`, and `omega_(0,1) = y dx` produces the dilaton leaves of DNOPS s.2.3,
      `-psi (I - M_0^-1(psi/(kappa varsigma^4))) 1hat`.
Consequently (DNOPS Thm 8.1 with s.2.6, at `kappa = 4`, plus the transfer of s.6 below) `omega_(h,n)` of `x = W_v` is given by (2-35) with
`alpha = M_0(./varsigma^4).omega^TR`, `omega^TR` the TFT with `Delta_j = Lambda_1 varsigma_j^2`, `Lambda_1 = -5`.

## 2. Formal Laplace series (`lem:mod-laplace`, PROVED, elementary; new lemma for `sec:mod-tr` or appendix A)

Let `c` be a simple zero of `dx_kappa`, `w` a coordinate at `c` with `x_kappa = x_kappa(c) - w^2/2`. For a 1-form `omega` meromorphic at
`c` with `Res_c omega = 0`, write `omega = Sum_(m >= -M) a_m w^m dw` and put `L_c[omega] := Sum_m a_m G_m(lam)`, with the Gaussian moments
  `G_(2n) := (2n-1)!! lam^-n` (`n >= 0`),   `G_(-2n) := (-1)^n lam^n/(2n-1)!!` (`n >= 1`),   `G_m := 0` (`m` odd, `m != -1`).
*Meaning.* `sqrt(2 pi/lam) L_c[omega]` is the asymptotic expansion, as `lam -> +oo`, of `Int omega e^(-lam w^2/2)` over the real `w`-line
deformed near `w = 0` into `Im w > 0` (Watson's lemma; the deformation does not matter because `Res_c omega = 0`). This is the reading of
the integrals in DNOPS (2-31), (2-32): for `i != j` and for (2-32) the integrand is holomorphic at `c`; for `i = j` its polar part is
`dw/w^2`, `G_(-2) = -lam`, which produces exactly the `delta_ij` of `R^-1 = I + O(lam^-1)` (DNOPS s.2.3 requires `R = I + O(z)`), i.e.
DOSS's prescription "delta_ij from the singular part `dw dw'/(w - w')^2` plus the transform of the regular part". (Recorded as finding
F1 in G4_5: the note should state this reading once, in `conv:mod-kernel` or appendix A.)
Properties (each checked on the monomials `w^m`, which span; no analysis is involved):
* (L1) *Integration by parts.* If `F` is meromorphic at `c` and `F dx_kappa` has zero residue at `c`, then `L_c[dF] = -lam L_c[F dx_kappa]`.
  (`F = w^m`, `m != -2`: `dF = m w^(m-1) dw`, `F dx_kappa = -w^(m+1) dw`, and `m G_(m-1) = lam G_(m+1)`: for `m = 2n+1 >= 1` both are
  `(2n+1)!! lam^-n`; for `m = -(2n+1) <= -1` both are `(-1)^n lam^(n+1)/(2n-1)!!`; for even `m` both vanish.)
* (L2) *Coordinate.* `w` is fixed by `x_kappa` up to `w -> -w`, which changes `L_c` by the factor `-1`.
* (L3) *Homogeneity.* Under `p = k p'`, `v = k^5 v'`: `x_kappa,v(k p') = k^4 x_kappa,v'(p')`, `varsigma_j(v) = k varsigma_j(v')`, the
  normalised coordinates satisfy `w_j(v)(k p') = k^2 w_j(v')(p')` (compare derivatives at the branch point), and therefore
  `L_(varsigma_j(v))[p^m dp](lam) = k^(m-1) L_(varsigma_j(v'))[p'^m dp'](k^4 lam)` (substitute in `[w^(2n)](p^m dp/dw)`).
* (L4) *Derivative in `lam`.* For `omega` holomorphic at `c`: `lam d/dlam L_c[omega] = lam L_c[(x_kappa - x_kappa(c)) omega] + L_c[omega]/2`
  (on `w^(2n) dw`: `-n G_(2n) = -(lam/2) G_(2n+2) + G_(2n)/2`).

## 3. Proof of (a)

*Step 1 (integration by parts).* `B(w_i, .)/dw_i |_(w_i = 0) = (dp/dw_i)(varsigma_i) dq/(varsigma_i - q)^2 = d xi^i`, `xi^i(q) :=
-1/(sqrt(Lambda) varsigma_i (q - varsigma_i))` (this is (2-34) with `xi^i(oo) = 0`). The function `xi^i` has at most a simple pole at
`varsigma_j` (only for `j = i`), so `xi^i dx_kappa` is holomorphic there and (L1) applies:
  `R^-1(lam^-1)^j_i = -(2 pi lam)^(-1/2) sqrt(2 pi/lam) L_j[d xi^i] = -lam^-1 L_j[d xi^i] = L_j[xi^i dx_kappa]`.
*Step 2 (partial fractions).* `dx_kappa = kappa (p^5 - v) dp/p^2` and `(p^5 - v)/(p - varsigma_i) = Sum_(k=0..4) varsigma_i^k p^(4-k)`
(`varsigma_i^5 = v`), so `xi^i dx_kappa = -(kappa/(sqrt(Lambda) varsigma_i)) Sum_(k=0..4) varsigma_i^k p^(2-k) dp`, and with
`J^(j)_m := L_j[p^m dp]`, `It^(j)_m := sqrt(Lambda) J^(j)_m` and `-kappa/Lambda = 1/5`:
  **`R^-1(lam^-1)^j_i = (1/5) Sum_(k=0..4) varsigma_i^(k-1) It^(j)_(2-k)`.**                                                   (1)
*Step 3 (the five integrals are the model rows).* Put `z := -1/(kappa lam)` (so `lam x_kappa = -W_v/z`) and `t := z/varsigma^4`.
* (Mod2) for the `J`: (L1) with `F = p^m` (holomorphic at `varsigma_j != 0`) gives `m J_(m-1) = -lam kappa (J_(m+3) - v J_(m-2))`, i.e.
  `m z J_(m-1) - J_(m+3) + v J_(m-2) = 0`; in particular `J_3 = v J_(-2)`.
* (Mod1) for `Yt_m := e^(-W_v(varsigma_j)/z) J_m` (local letter, indexed by the power `m`): by (L3), `5 v d_v J_m = (m-1) J_m + 4 lam
  d_lam J_m`; by (L4) with `x_kappa -
  x_kappa(varsigma_j) = kappa (p^4/4 + v/p) - (5/4) kappa varsigma_j^4`, `lam d_lam J_m = lam kappa (J_(m+4)/4 + v J_(m-1)) - (5/4) lam kappa
  varsigma_j^4 J_m + J_m/2`; eliminating `J_(m+4)` by (Mod2) at `m+1` and using `lam kappa = -1/z`, `varsigma_j^4 = v/varsigma_j`:
  `z d_v J_m = -J_(m-1) + J_m/varsigma_j`. Since `d_v [W_v(varsigma_j(v))] = 1/varsigma_j` (`W_v'(varsigma_j) = 0`): **`z d_v Yt_m = -Yt_(m-1)`**.
* Hence `Y_b := eps_b Yt_(nu_b)` (for `b = 2`: `Y_2 = -Yt_3 = -v Yt_(-2)`) satisfy `z d_v Y_b = bmod_b Y_(b+1)` (`b` in `Z/5`): this is the derivation of
  `prop:set-modrows` (a) from (Mod1), (Mod2), now valid for these formal series. Write `Y_b = e^(U_j/z) varsigma^(nu_b - 1) F_b`, `U_j =
  -W_v(varsigma_j) = -(5/4) zeta^j varsigma^4`, `F_b := eps_b varsigma^(1-nu_b) J^(j)_(nu_b)`. By (L3) (`k = varsigma`), `F_b` depends on
  `varsigma^4 lam = -1/(kappa t)` only, i.e. it is a power series in `t`. So `(F_b)` solves the colour-`j` T-system of
  `prop:set-modrows` (b) (the computation there: `z d_v [e^(U_j/z) varsigma^(nu_b - 1) F_b(t)] = e^(U_j/z) varsigma^(nu_b - 2) T_(nu_b - 1) F_b`).
* Initial values: `J^(j)_m = varsigma_j^m (dp/dw_j)(varsigma_j) + O(lam^-1) = varsigma_j^(m-1)/sqrt(Lambda) + ..`, and `(varsigma_j/varsigma)^(nu_b - 1)
  = zeta^(-j(nu_b - 1)) = zeta^(j(b+1))` (`nu_b - 1 = -(b+1)` mod 5), so `sqrt(Lambda) F_b(0) = eps_b zeta^(j(b+1)) = Rhat^mod_b^(j)(0)`.
* `Rhat^mod^(j)` solves the colour-`j` T-system: with `s := zeta^-j t`, `T^(j)_k [zeta^(j(b+1)) Fs_b(s)] = zeta^(j(b+2)) (T^(0)_k Fs_b)(s)`
  (`t d/dt = s d/ds`), so the colour-0 relations `Fs_(b+1) = +-T_k Fs_b` become the colour-`j` ones (`zeta^(5j) = 1` closes `4 -> 0`).
* Uniqueness (`prop:set-modrows` (b), pivot `4n`: the closing equation `Fs_0 = -T_0 Fs_4` reads `4 n f_n + (terms in f_(n-4..n-1)) = 0` at
  order `t^(n+1)`; for colour `j` the pivot is `4 n zeta^(4j)`) gives **`sqrt(Lambda) F_b = Rhat^mod_b^(j)(t)`**, i.e. `It^(j)_(m_b) = eps_b varsigma^(m_b - 1) Rhat^mod_b^(j)(t)` for
  `(m_0, .., m_4) = (0, -1, -2, 2, 1)` (`m_b = -b` mod 5; for `b = 2` use `It_(-2) = It_3/v`).
*Step 4 (assembly).* Insert into (1) with `k = 2 - m_b`: `varsigma_i^(1-m_b) varsigma^(m_b - 1) = zeta^(-i(1 - m_b)) = zeta^(-i(b+1))`, so
  `R^-1(lam^-1)^j_i = (1/5) Sum_b eps_b zeta^(-i(b+1)) Rhat^mod_b^(j)(t) = Sum_b (Psi_0^-1)_(ib) Rhat^mod_b^(j)(t) = M_0(t)_(ij)`. 
*Step 5 (inverse).* As matrices (row = upper index): `R^-1(w) = M_0(-w/(kappa varsigma^4))^T = M_0(w/(kappa varsigma^4))^-1` by
`lem:set-M0sympl`. Hence `R_TR(w) = M_0(w/(kappa varsigma^4))`. QED (a). [One could instead use that the DNOPS matrix is in the group
of their s.2.3 (`R(z) R^*(-z) = I`; their s.5 writes `R^-1(xi,u) = R(-xi,u)^T`; they state that (2-33) "is always satisfied if the curve
is compact and the differential dx is meromorphic"); then the symplecticity of `M_0` would follow. Not needed: the proof of record uses
`lem:set-M0sympl`, which is elementary.]

## 4. Proof of (b)

`Delta_k^(-1/2) = dy/dw_k(0) = 1/(sqrt(Lambda) varsigma_k)` (2-30). Left side of (2-32): by (a),
`Sum_k M_0(t)_(ki)/(sqrt(Lambda) varsigma_k) = (1/(sqrt(Lambda) varsigma)) Sum_b eps_b [(1/5) Sum_k zeta^(k - k(b+1))] Rhat^mod_b^(i)(t) =
Rhat^mod_0^(i)(t)/(sqrt(Lambda) varsigma)` (`Sum_k zeta^(-kb) = 5 delta_(b0)`). Right side: `sqrt(lam/2 pi) sqrt(2 pi/lam) L_i[dp] = J^(i)_0 =
It^(i)_0/sqrt(Lambda) = Rhat^mod_0^(i)(t)/(sqrt(Lambda) varsigma)` (Step 3, `b = 0`). Equal, to all orders. The unit statement is then (a)
applied to `Sum_k Delta_k^(-1/2) e_k = Sum_k et_k/varsigma_k`. QED (b).

## 5. `lem:set-M0sympl` (PROVED; new lemma, belongs in s02 next to `def:set-M0`; closes the R7 item "def:set-M0 symplecticity")

**Lemma.** `M_0(-t)^T M_0(t) = 1` (hence also `M_0(t) M_0(-t)^T = 1`), i.e. `M_0` lies in Givental's group (DNOPS s.2.3: `R = I + O(z)`,
`R(z) R^*(-z) = I` in the orthonormal basis).
*Proof.* Let `Y = (Y^(i)_b)` (row `b`, column `i`), `Y^(i)_b(v, z) := e^(U_i/z) varsigma^(nu_b - 1) Rhat^mod_b^(i)(z/varsigma^4)` on a simply
connected `v`-domain with a branch of `varsigma` (formal: `e^(U_i/z)` times a series in `z`). By `prop:set-modrows` (b) and the colour
rotation of Step 3, `z d_v Y = A Y`, `A := Sum_b bmod_b E_(b,b+1)`; substituting `z -> -z`, `-z d_v Y(v,-z) = A Y(v,-z)`. Put `P :=
Y(v,-z)^T eta^-1 Y(v,z)`. Then `d_v P = z^-1 Y(v,-z)^T (eta^-1 A - A^T eta^-1) Y(v,z) = 0`, because `(eta^-1 A)_(ac) = (1/5) bmod_(3-a)
delta_(c, 4-a)` is symmetric: `bmod_(3-a) = bmod_(a-1)` for all `a` (`bmod = (1, -1/v, 1, -1, -1)`). On `eta`-pairs `nu_a + nu_b = 2`, so
`P_(ij) = e^((U_j - U_i)/z) f_(ij)(t)`, `f_(ij)(t) := (Rhat^mod(-t)^T eta^-1 Rhat^mod(t))_(ij) = Sum_n f_(ij,n) t^n`, `f_(ij,n)` constants.
Write `f_(ij) = Sum_n f_(ij,n) varsigma^(-4n) z^n` and `Delta := U_j - U_i`. `d_v P = 0` means `Delta' f + z d_v f = 0` termwise in `z`: the
`z^n` coefficient is `Delta' f_(ij,n) varsigma^(-4n) + f_(ij,n-1) d_v(varsigma^(-4(n-1))) = 0`. For `i != j`,
`Delta' = -(zeta^j - zeta^i) varsigma^-1 != 0`, so `f_(ij,0) = 0` and inductively all `f_(ij,n) = 0`. For `i = j`, `Delta = 0` and
`d_v(varsigma^(-4n)) != 0` gives `f_(ii,n) = 0` for `n >= 1`; `f_(ii,0) = (Psi_0^T eta^-1 Psi_0)_(ii) = 1`
(`(Psi_0^T eta^-1 Psi_0)_(ij) = (1/5) Sum_b eps_b eps_(3-b) zeta^((i-j)(b+1)) = delta_ij`, as `eps_b eps_(3-b) = 1`). So
`Rhat^mod(-t)^T eta^-1 Rhat^mod(t) = 1`, and with `Psi_0 Psi_0^T = eta` (equivalent to `Psi_0^T eta^-1 Psi_0 = 1`):
`M_0(-t)^T M_0(t) = Rhat^mod(-t)^T (Psi_0 Psi_0^T)^-1 Rhat^mod(t) = 1`. QED
(The same argument is Dubrovin's/Givental's for the stationary-phase solutions; here everything is the explicit cyclic system. It also
re-proves LG s.3's VERIFIED symplecticity "through zhat^40".)

## 6. The transfer `kappa = 4 -> 1` (R3; PROVED)

**`lem:mod-scale` (PROVED).** For constants `kappa, mu != 0` let `omega^(kappa,mu)_(h,n)`, `F^(kappa,mu)_h` be the TR data of
`(P^1, kappa x, mu y, B)`. Then `omega^(kappa,mu)_(h,n) = (kappa mu)^(2-2h-n) omega_(h,n)` (`2h - 2 + n > 0`) and `F^(kappa,mu)_h =
(kappa mu)^(2-2h) F_h` (`h >= 2`). *Proof.* Same branch points and involutions; the kernel becomes `K/(kappa mu)`; `omega_(0,2) = B` is
unchanged (exponent `0 = 2 - 0 - 2`); each term of the recursion for `omega_(h,n+1)` is one kernel times data whose exponents add up to
`2 - 2h - n` (both for `omega_(h-1,n+2)` and for the products, the `omega_(0,2)` factors included), so induction on `2h - 2 + n` gives
`(kappa mu)^(1-2h-n)`; for `F_h`, `Phi_ydx -> kappa mu Phi_ydx`. QED (TR is natural under biholomorphisms of the curve; with `p = k p'`,
`v = k^5 v'` this gives the homogeneity `omega^(v) -> k^(5(2-2h-n)) omega^(v')` of `def:mod-curve`, TRRF s.3.4.)
**Transfer.** DNOPS Thm 8.1 with s.2.6 at `f = 4 W_v` (hypotheses: R3; TRM s.1; TRRF claim 1) gives (2-35) for `omega^4` with `alpha^4 =
R^4.omega^4_TFT`, `R^4(w) = M_0(w/(4 varsigma^4))` (s.3 at `kappa = 4`), `Delta^4_j = Lambda_4 varsigma_j^2`. Take `sqrt(Lambda_4) =
2 sqrt(Lambda_1)`; then `w^4_j = 2 w^1_j`, `xi^(4,i) = xi^(1,i)/2`, `d/dx_4 = (1/4) d/dx_1`, `Delta^4_j = 4 Delta^1_j`, and in the orthonormal
bases `e^4_j <-> e^1_j`: (i) `omega^4_TFT,(h,n) = 2^(2h-2+n) omega^1_TFT,(h,n)` (the only non-zero values are `omega_(h,n)(e_j,..,e_j) =
Delta_j^((2h-2+n)/2)`), with units `e^4_1 = e^1_1/2`; (ii) `R^4(w) = R^1(w/4)`, `R^1(w) := M_0(w/varsigma^4)`. (iii) Therefore `alpha^4_(h,n) =
2^(2h-2+n) Sum_D 4^(-D) [alpha^1_(h,n)]_D`, `alpha^1 := R^1.omega^1_TFT`, `[.]_D` the part of cohomological degree `D`: in the graph sum of
DNOPS s.2.3 the vertices give `2^(Sum_v (2h_v - 2 + n_v)) = 2^(2h - 2 + n + #dil)` and the dilaton leaves (decorated by the unit) give
`2^(-#dil)`; `R^1(./4)` multiplies a leaf or dilaton-leaf decoration `psi^k` by `4^-k` and an edge coefficient of `psi'^a psi''^b` by
`4^-(a+b+1)` (`(I (x) I - R^-1(psi'/4) (x) R^-1(psi''/4))/(psi' + psi'') = (1/4) E^1(psi'/4, psi''/4)`), i.e. every class of degree `d` by
`4^-d` (push-forwards preserve degree). (iv) In (2-35) only `D = 3h - 3 + n - Sum d_l` contributes and leg `l` gives `4^(-d_l)/2`, so
`omega^4_(h,n) = 2^(2h-2+n) 4^-(3h-3+n) 2^-n RHS^1 = 4^(2-2h-n) RHS^1`, `RHS^1 := Sum Int alpha^1(e^1_(i_1), ..) prod psi^(d_l)
d((d/dx_1)^(d_l) xi^(1,i_l))`. By `lem:mod-scale`, `omega^1_(h,n) = RHS^1`: (2-35) holds on `x = W_v` with `alpha^1 = M_0(./varsigma^4).
omega^1_TFT`, `Delta^1_j = -5 varsigma_j^2`. (This is TRM s.2.4's `Lambda^(g-1) kappa^-D` bookkeeping, done here on the TR side.)

## 7. Remarks for the editor

* What is used from DNOPS: the sentence of s.8 just before Thm 8.1 (the spectral curve `(CP^1, f(p,a), p, dp_1 dp_2/(p_1 - p_2)^2)` "satisfies
  equations (2-30)-(2-33)" and hence stores the correlators of a CohFT via (2-35)), with Thm 8.1 itself (whose proof runs through Thm 5.1:
  the compatibility of their Thm 4.1 holds, so the CohFT reconstructed from (2-30)-(2-32) exists and satisfies (2-35)); the definitions
  (2-30), (2-31), (2-32), (2-34); the action of s.2.3.
  The second half of Thm 8.1 (that this CohFT is the Frobenius manifold's) is NOT needed: the note identifies `alpha` with
  `Omega^mod_v` directly through (a). So "homogeneity of `R_TR`" (DNOPS Prop. 5.2), the uniqueness of homogeneous R-matrices (TRM s.2.3;
  the suggestion MR2917177 8.4 in OUTLINE_6), the canonical-frame recursion `(5/4)[C_(n+1), Z] = (V - n) C_n` (LG s.3, VERIFIED only) and
  TRM's sign-of-`kappa` argument all drop out of the chain. The suggested route "homogeneous model column unique (UR Prop. S, pivot `4n`)"
  is exactly Step 3: the rows of `R_TR` ARE the homogeneous model columns.
* Conventions that must be stated once (with `def:set-modcohft`): the R-action is DNOPS s.2.3 = [PPZ15] Def. 2.2 (leaves `R^-1(psi)` on
  column vectors, row index = output colour); DNOPS Remark 2.1: Givental's convention is the transpose. The orientation choice
  `dw_j/dp(varsigma_j) = sqrt(Lambda) varsigma_j` (uniform in `j`) is the one for which (a) holds as printed; another choice conjugates by a
  diagonal sign matrix together with `e_j -> -e_j` and leaves the CohFT unchanged.
* The dilaton leaves: DNOPS decorate them by `-psi (I - R^-1(psi)) e_1`; the outline's `z(1hat - M_0^-1(z) 1hat)` is the same translation
  with the opposite sign convention for the leaf: fix it once in `conv:mod-kernel`/appendix A.
* TRM s.2.3's VERIFIED statement "`X(z) = M_0(z/kappa)^-1 = M_0(-z/kappa)^T`, z^0..z^8" is (a) at `v = 1` with TRM's `z = lam^-1`.

## 8. Evidence (VERIFIED; `g4_rtr.py 10` -> `g4_rtr_10.log`, exact over `Q(zeta5)`, 10 s)

* `Fs_0 = 1 - (3/20)t + (9/800)t^2 + (11571/80000)t^3 - ..` from the T-system with pivot `4n`; closing relation checked to `t^10`.
* `M_0(-t)^T M_0(t) = 1` through `t^10` (0 failures of 275 entries).
* For `(kappa, varsigma) = (1,1), (4,1), (-1,1), (1,2), (4,2)`: (2-31) computed from the local expansion of `xi^i` in `w_j` with the finite
  part at `i = j` (independent of Steps 1-4): `X = M_0(t)^T`, `t = -1/(kappa varsigma^4 lam)`, 0 mismatches in 275 coefficients each; the
  integration-by-parts form (1) with `It^(j)_m` computed from its definition: 0 mismatches; (2-32): 0 mismatches.
* Controls (must fail, do fail): `X = M_0(t)` (no transpose): 180 mismatches; `X = M_0(-t)^T`: 125; `kappa` dropped from `t` (for
  `kappa = 4, -1`): 250, 125.
* Not evidence for anything beyond the stated orders; the proof is s.3-s.6.
