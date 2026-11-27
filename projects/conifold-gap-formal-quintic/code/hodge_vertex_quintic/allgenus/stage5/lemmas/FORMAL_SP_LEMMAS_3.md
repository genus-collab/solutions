# Stage 5 / lemmas, part 3: Fubini (L1b), the delta lemma (L2), the interface for the chart argument

Continuation of `FORMAL_SP_LEMMAS.md` and `_2.md`. Same notation and labels.

## 6. Lemma L1b: formal Fubini (PROVED)

**Lemma L1b.** Variables `(x, y)`, `n + m` of them. Let `F in A[[x - x_c, y - y_c]]` have the critical point `(x_c, y_c)` with Hessian
`H = [[H_xx, H_xy], [H_yx, H_yy]]`, where `H_xx in GL_n(A)` and `S := H_yy - H_yx H_xx^-1 H_xy in GL_m(A)` (given `H_xx` invertible, `S`
invertible <=> `H` invertible, `det H = det H_xx det S`). Let `a(x, y; z)` be an amplitude and `A_y := A[[y - y_c]]`.
1. *(Inner critical point.)* There is a unique `xi(y) in (x_c + (y - y_c) A_y)^n` with `d_x F(xi(y), y) = 0`; the inner Hessian
   `K(y) := d_x^2 F(xi(y), y)` lies in `GL_n(A_y)`, `K(y_c) = H_xx`. So `x -> F(x, y)` has a nondegenerate critical point over `A_y`.
2. *(Inner data.)* `G(y) := F(xi(y), y) in A_y` (inner critical value), `Delta_in(y) := det(-K(y)) in A_y^x`, and the inner series
   `b(y; z) := SP^norm_(xi(y))[F(., y); a(., y)] in A_y((z))` (Definition 1.2 over the coefficient ring `A_y`).
3. *(Outer phase.)* `grad G(y_c) = 0`, `Hess G(y_c) = S`, and `Delta := det(-H) = Delta_in(y_c) det(-S)`.
4. *(Fubini.)*  `SP^norm_((x_c, y_c))[F; a] = SP^norm_(y_c)[ G ; (Delta_in(y)/Delta_in(y_c))^(-1/2) b(y; z) ]`  (binomial series).

*Full form.* `SP[e^(F/z) a dx dy] = SP_(y_c)[ e^(G/z) (2 pi z)^(n/2) Delta_in(y)^(-1/2) b(y; z) dy ]`, with the branch rules
`Delta_in(y)^(1/2) := Delta_in(y_c)^(1/2) (Delta_in(y)/Delta_in(y_c))^(1/2)` and `Delta^(1/2) := Delta_in(y_c)^(1/2) det(-S)^(1/2)`.
So: the inner SP is taken at the inner critical point `xi(y)`; the inner critical value `G(y)` becomes the outer PHASE; the
y-dependent part of the inner Gaussian factor, `(Delta_in(y)/Delta_in(y_c))^(-1/2)`, and the inner series `b(y; z)` become the outer
AMPLITUDE; the constants `(2 pi z)^(n/2) Delta_in(y_c)^(-1/2)` multiply. Every splitting of the variables into an inner block with
invertible Hessian and the rest computes the same total SP; so two orders of integration give equal results (chart step (c)).
*Proof.* (1) IFT of s.3.2 over `B = A_y`, `I = (y - y_c)`, for `x' -> F(x_c + x', y)`: `d_x F(x_c, y) in I^n`, and the x-Hessian at
`(x_c, y)` is invertible over `A_y` (its value at `y_c` is `H_xx`).
(3) `grad G(y) = d_y F(xi(y), y)` (the `d_x F` term vanishes), so `grad G(y_c) = 0`. Differentiating `d_x F(xi(y), y) = 0` at `y_c` gives
`D xi(y_c) = -H_xx^-1 H_xy`, hence `Hess G(y_c) = H_yy + H_yx D xi(y_c) = S`. The determinant formula is Schur's.
(4) *Step 1 (straighten the inner critical point).* `Theta(x', y) := (xi(y) + x', y)` has `Theta(0, y_c) = (x_c, y_c)`,
`D Theta = [[I, D xi], [0, I]]`, `det D Theta = 1`. By L1a the left side equals `SP^norm_((0, y_c))[F o Theta; a o Theta]`. The right side
is the same for `(F, a)` and `(F o Theta, a o Theta)`: the inner SP of `F o Theta` at `x' = 0` has the same Taylor data in `x' = x - xi(y)`
as the inner SP of `F` at `xi(y)`, and `G`, `Delta_in` do not change. So we may assume `xi = 0`, and translate `y_c = 0`. Then
`d_x F(0, y) = 0` identically, hence `H_xy = 0`, `S = H_yy`, `G(y) = F(0, y)`.
*Step 2 (Gaussian factorisation).* Write `F(x, y) = G(y) + x^T K(y) x/2 + R(x, y)`, `R in (x)^3` over `A_y`. The total Hessian is
`diag(H_xx, H_yy)`; the total cubic part is `F_3 = G_3(y) + x^T K_1(y) x/2 + R(x, y)` with `G_3 := G - y^T H_yy y/2` and `K_1 := K - H_xx`
(entries in `(y)`, of weight `>= 1/2`). By (G4), then (G5') in the x-variables with `y, z` as coefficients:
`SP^norm_0[F; a] = E_(-zH_yy^-1)[ e^(G_3/z) E_(-zH_xx^-1)[ e^(x^T K_1 x/(2z)) e^(R/z) a ] ]`
`              = E_(-zH_yy^-1)[ e^(G_3/z) (det K(y)/det H_xx)^(-1/2) E_(-zK(y)^-1)[ e^(R/z) a ] ]`.
The inner functional is `b(y; z)` (Definition 1.2 for `x -> F(x, y) - G(y)`: Hessian `K(y)`, cubic part `R(., y)`), and
`det K/det H_xx = Delta_in(y)/Delta_in(0)`. The remaining expression is `SP^norm_0[G; (Delta_in/Delta_in(0))^(-1/2) b]`. QED
*Rings.* The inner functional `E_(-zK(y)^-1)` (covariance in `GL_n(A_y)`) maps `R_(n+m)(A)` to `R_m(A)` without lowering the weight, and
it produces the same element as the inner SP computed in `A_y((z))`: the coefficient of each `y^beta z^k` is the same finite sum.
*Remark (uniformity).* An analytic Fubini for Laplace integrals needs the inner expansion to be uniform in the outer variables. In
the formal setting this is replaced by the fact that the inner data `xi, K, R, a` are power series in `y - y_c`; nothing else is used.

## 7. Lemma L2: the formal delta lemma (PROVED)

**Lemma L2.** Variables `sigma`, `t` (one each) and `w` (`m` of them). `Phi(sigma, t, w) = sigma G(t, w) + f(t, w)` with
`G, f in A[[t - t_c, w - w_c]]`, and an amplitude `a(t, w; z)` NOT involving `sigma`. Let `(sigma_c, t_c, w_c)` be critical:
`G = 0`, `sigma_c G_t + f_t = 0`, `sigma_c G_w + f_w = 0` at `(t_c, w_c)`. Assume `g_c := G_t(t_c, w_c) in A^x`.
1. The `(sigma, t)`-block of `Hess Phi` is `[[0, g_c], [g_c, sigma_c G_tt + f_tt]]`, of determinant `-g_c^2`, invertible.
2. Over `A_w := A[[w - w_c]]` there is a unique `tau(w) in t_c + (w - w_c) A_w` with `G(tau(w), w) = 0`. The inner critical point is
   `(sigma*(w), tau(w))`, `sigma*(w) := -f_t/G_t` at `(tau(w), w)`; the inner critical value is `f_res(w) := f(tau(w), w)`, the restriction
   of `f` to `{G = 0}`; and `Delta_in(w) = -G_t(tau(w), w)^2`.
3. *(Exactness.)*  `SP^norm_((sigma*(w), tau(w)))[ Phi(., ., w) ; a(., w) ] = a(tau(w), w; z)`  exactly: no correction at any order.
4. *(With L1b.)* If the total Hessian is invertible (equivalently `Hess f_res(w_c)` is, by L1b(3)), then
   `SP^norm_((sigma_c, t_c, w_c))[Phi; a] = SP^norm_(w_c)[ f_res ; a(tau(w), w; z) g_c / G_t(tau(w), w) ]`,
   and in full form `SP[e^(Phi/z) a dsigma dt dw] = SP_(w_c)[ e^(f_res/z) (2 pi z) (-G_t(tau(w), w)^2)^(-1/2) a(tau(w), w; z) dw ]`,
   with `(-G_t(tau(w), w)^2)^(-1/2) := (-g_c^2)^(-1/2) g_c/G_t(tau(w), w)`: restriction to `G = 0` with the Jacobian `1/G_t`, to all orders.
*Proof.* (1) Direct (`Phi` is linear in `sigma`). (2) `tau` is the fixed point of `t -> t - G(t, w)/g_c` (the contraction argument of
s.3.2); then `d_t Phi = 0` fixes `sigma*(w)` and `d_sigma Phi = G = 0`. The value at the critical point is `sigma* . 0 + f(tau, w)`.
(3) Fix `w` (coefficients in `A_w`). Let `T(g) in g A_w[[g]]` be the compositional inverse of `t' -> G(tau(w) + t', w)`, so that
`G(tau + T(g), w) = g`. Change variables by `(s, g) -> (sigma, t) = (sigma*(w) + s, tau(w) + T(g))` (L1a over `A_w`): it maps `(0, 0)` to the
inner critical point; the Jacobian is `T'(g) = 1/G_t(tau + T(g), w)`, normalised `j(g) := G_t(tau, w)/G_t(tau + T(g), w)`, `j(0) = 1`. The
phase becomes `(sigma* + s) g + f(tau + T(g), w) = f_res(w) + s g + fh(g)`, `fh(g) := f(tau + T(g), w) - f(tau, w) + sigma* g`, with
`fh(0) = 0` and `fh'(0) = f_t/G_t + sigma* = 0`. The Hessian in `(s, g)` is `[[0, 1], [1, h]]`, `h := fh''(0)`, so
`C = -[[0, 1], [1, h]]^-1 = [[h, -1], [-1, 0]]`:  **`C_gg = 0`**. The cubic part `fh(g) - h g^2/2` and the amplitude
`ah(g) := a(tau + T(g), w; z) j(g)` depend on `g` alone (the term `s g` belongs to the quadratic form). Hence
`SP^norm = E_zC[ e^((fh(g) - h g^2/2)/z) ah(g) ]` is a combination of the moments `<g^k>_C`, and `<g^k>_C = 0` for `k >= 1`: every perfect
matching of `k` copies of the index `g` contains the factor `C_gg = 0`. Only `k = 0` survives: `SP^norm = ah(0) = a(tau(w), w; z)`. QED
(4) L1b with inner block `(sigma, t)` and outer `w`: the outer phase is `f_res`; the outer amplitude is
`(Delta_in(w)/Delta_in(w_c))^(-1/2) a(tau(w), w; z) = (rho^2)^(-1/2) a = rho^-1 a`, `rho := G_t(tau(w), w)/g_c in 1 + (w - w_c) A_w` (in a
Q-algebra, `u^2 = v^2` with `u, v in 1 + (w - w_c)A_w` implies `u = v`, since `u + v` is a unit). QED
*Meaning.* `int_(iR) e^(sigma G/z) dsigma = 2 pi i z delta(G)` and `int delta(G(t)) a dt = a(tau)/|G_t(tau)|`; the formal factor
`(2 pi z)(-g_c^2)^(-1/2) = 2 pi z/(+-i g_c)` is this, with the orientation of the sigma-contour fixing the sign.
*In chart step (b)* (STAGE4_SYNTHESIS s.1.4): `t = x_0`, `w = (y, x_1, .., x_4)`, `G = log(-q) + 5 log y - Sum_(j=0..4) log x_j`,
`G_(x_0) = -1/x_0`, a unit iff `x_(0,c) = 1 - L != 0`; `tau(w) = X_0 := -q y^5/prod_(j>=1) x_j`; the amplitude `z^-1 prod_(j=0..4) x_j^-1` gives
`a(tau(w), w) g_c/G_t(tau(w), w) = z^-1 x_(0,c)^-1 prod_(j>=1) x_j^-1`: the Jacobian factor `X_0` cancels the `1/X_0` of the measure.

## 8. Interface: what the chart argument needs (for `stage5/chart/`)

The chart argument (STAGE4_SYNTHESIS s.1.4) is a statement about formal SPs over a field `K` containing `zeta` and `L` (e.g.
`Q(zeta)(L)`, or `C` with numerical `L`), and over power-series rings in `sigma - L`, `w - w_c`, `q`. The lemmas apply as follows.

| step | lemma | coefficient ring | hypotheses to be checked by the chart proof |
|---|---|---|---|
| (a) inner `(x_0..x_4, y)` at fixed `sigma` | L1b; (G4) for the product of the six inner integrals; L3 with `a = zeta^j - sigma` (`dx/x`, `m = 0`) and `a = 5 sigma` (`dx`, `m = 0`) | `A_sigma = K[[sigma - L]]` | `zeta^j - L != 0` (`j = 0..4`) and `L != 0`; the outer Hessian in `sigma` (the Schur complement) is a unit; the outer amplitude `Delta_in(sigma)^(-1/2) prod(Stirling series)` is the full-Stirling amplitude of GAMMA s.2 / POLY s.5 up to constants |
| (b) inner `(sigma, x_0)` at fixed `w` | L2 (with item 4) | `A_w = K[[w - w_c]]` | `x_(0,c) = 1 - L != 0`, i.e. `L != 1`; the amplitude is free of `sigma` (the measure `dsigma dy prod dx_j/x_j`) |
| (c) equality of (a) and (b) | L1b, twice | `K` | the 7 x 7 Hessian is invertible at `(L, 5L, zeta^j - L)`; both sides are then the same full SP. Prefactors: `(2 pi z)^(7/2)` on both sides; the critical values agree |
| (d) `q = 0` | L1c (i)-(ii); (G4); L3 (`a = 5`, `dy`; `a = zeta^j - 1`, `dx_j/x_j`; `m = 0`) | chart variables over `K'[[q]]` or holomorphic in `q` | the chart phase, its critical point `(5L, zeta^j - L)` and its Hessian are analytic at `q = 0`, Hessian `diag(-1/(zeta^j - 1), -1/5)` at `q = 0` |
| (e) q-expansion | L1c (iii) (T-q) with `delta = q y^5/prod x_j` at the `q = 0` critical point; (G4); L3 (v) with `m = 5d` (`y`) and `m = -d` (`x_j`) | `K'[[q]]` | NEGATIVE `m` is needed for the `x_j`; L3 is proved for all `m in Z` |

*Two warnings.* (1) The constant prefactor of the inner `(sigma, x_0)` integral contains `(-g_c^2)^(-1/2)`, `g_c = -1/(1 - L)`: it is
proportional to `1 - L`. The 7-variable full SP is therefore not defined at `L = 1` (its critical point leaves the domain of
`log x_0`), and the comparison (c) must be made at `L != 1` and continued coefficientwise; the normalised chart series is unaffected.
(2) The continuation "a rational function of `L` that agrees on a punctured neighbourhood of `L = 1` with a holomorphic function is
regular at `L = 1`" is not a lemma of this note; it is elementary (identity theorem), but it needs the chart coefficients to be
holomorphic in `L` near `L = 1`, which is L1c(i) applied to the chart data with `q = (1 - L^-5)/5^5` as parameter.
