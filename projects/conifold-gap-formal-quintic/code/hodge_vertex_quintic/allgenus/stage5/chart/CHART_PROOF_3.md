# Stage 5 / chart, part 3: the chart argument, steps (a)-(e); covariance; proof of Theorem D (s.3)

Continuation of `CHART_PROOF.md` (s.0-s.1) and `CHART_PROOF_2.md` (s.2, lemmas). `L` is a fixed point of
`Lam = C \ {0, 1, zeta, .., zeta^4}`, `q = q(L) = (1 - L^-5)/5^5`, `lambda` the logarithm of `-q` of s.1.6, and
`Lam1(L) := prod_(j=1..4) (zeta^j - L) = 1 + L + L^2 + L^3 + L^4` (so `prod_(j=0..4) (zeta^j - L) = 1 - L^5 = (1 - L) Lam1(L)`).

## 3.1 Step (a), first half: the 7-variable formal integrand; the powers of z cancel identically in sigma

*Motivation (not used in the proofs).* For `Re w > 0`, `Gamma(w) = Int_0^oo t^w e^-t dt/t`; with `t = x/z` this is
`z^(-w) Int_0^oo x^w e^(-x/z) dx/x`. Apply it to the six factors of the summand, `sigma = z s`:
* `Gamma(1 + 5 sigma/z) = z^(-1 - 5 sigma/z) Int e^((5 sigma log y - y)/z) dy`;
* `Gamma((zeta^j - sigma)/z) = z^(-(zeta^j - sigma)/z) Int e^(((zeta^j - sigma) log x_j - x_j)/z) dx_j/x_j`, `j = 0..4`.

The product of the six powers of `z` is `z^E` with `E = -1 - (1/z)[5 sigma + Sum_j (zeta^j - sigma)] = -1 - (1/z) Sum_j zeta^j = -1`,
**identically in `sigma`** (an identity in `Q(zeta)[sigma, 1/z]`: the `sigma`-coefficient is `5 - 5 = 0`, the constant is
`-Sum_j zeta^j = 0`). This is what makes the 7-variable phase free of `log z`: a `sigma`-dependent power `z^(c sigma/z)` would add
`c sigma log z` to the phase and move the critical point by a `log z`-dependent amount. VERIFIED symbolically (`ch_zpow.py` (i));
controls: dropping the factor `j = 0`, or `Gamma(1 + 4s)` for `Gamma(1 + 5s)`, give `sigma`-dependent exponents.

**Definition.** `Omega_7 := z^-1 e^(Phi/z) (x_0 x_1 x_2 x_3 x_4)^-1 d sigma dy dx_0 .. dx_4`, with
`Phi(sigma, y, x) := lambda sigma - y + 5 sigma log y + Sum_(j=0..4) (-x_j + (zeta^j - sigma) log x_j)`, where `log y` is a branch near
`5L` and `log x_j` near `zeta^j - L` (s.4). It is a formal integrand in 7 variables (s.1.1); `(-q)^s = e^(lambda sigma/z)`.

## 3.2 The critical point and every Hessian used, in both orders

* **Critical point** `p(L) = (sigma, y, x_0, .., x_4) = (L, 5L, zeta^0 - L, .., zeta^4 - L)`:
  `d_sigma Phi = lambda + 5 log y - Sum_j log x_j` vanishes at `p(L)` by the definition of `lambda` (same branches);
  `d_y Phi = -1 + 5 sigma/y`, `d_(x_j) Phi = -1 + (zeta^j - sigma)/x_j` vanish there. Exponentiated, the first equation is
  `-q (5L)^5 = prod_j (zeta^j - L) = 1 - L^5`, i.e. `u L^5 = 1`.
* **Hessian** `H_7` at `p(L)` (rows and columns `sigma, y, x_0..x_4`): `H_(sigma sigma) = 0`, `H_(sigma y) = 1/L`,
  `H_(sigma x_j) = -1/(zeta^j - L)`, `H_(y y) = -1/(5L)`, `H_(x_j x_j) = -1/(zeta^j - L)`, all others 0. With the rows multiplied by
  `L(1 - L^5)`, `5L`, `zeta^j - L` it becomes a polynomial matrix of determinant 5 (Bareiss over `Q(zeta)[L]`, `ch_hess.py`), so
  `det H_7 = 1/(L^2 (1 - L^5)^2)`, non-zero on `Lam`.
* **Order (a)**, `v = (y, x_0..x_4)` first: `H_vv = diag(-1/(5L), -1/(zeta^j - L))`, `det H_vv = 1/(5L(1 - L^5)) != 0`; Schur
  complement in `sigma`: `S_sigma = 0 - (1/L)^2/(-1/(5L)) - Sum_j (zeta^j - L)^-2/(-(zeta^j - L)^-1) = 5/L + Sum_j 1/(zeta^j - L)
  = 5/(L(1 - L^5)) = S_0''(L)`, using `Sum_j 1/(zeta^j - L) = 5L^4/(1 - L^5)`. Product: `1/(L^2 (1 - L^5)^2) = det H_7`.
  At any `sigma` near `L` the inner Hessian `diag(-1/(5 sigma), -1/(zeta^j - sigma))` stays invertible.
* **Order (b)**, `v = (sigma, x_0)` first: `H_vv = [[0, -1/(1-L)], [-1/(1-L), -1/(1-L)]]`, `det H_vv = -1/(1 - L)^2 != 0` (needs
  `L != 1`). Only the `sigma`-column couples to `w = (y, x_1..x_4)`, so the Schur complement is the rank-one update
  `H_ww - (1 - L) hh^T`, `h = (1/L, -1/(zeta^1 - L), .., -1/(zeta^4 - L))`, which is EXACTLY the Hessian of the chart phase `psi`
  (s.3.4) at `p_ch(L)` (entrywise; `ch_hess.py`). In log coordinates, `D H_ch D = diag(-5L, -(zeta^j - L)) + (L - 1) v v^T` with
  `D = diag(5L, zeta^j - L)`, `v = (5, -1, -1, -1, -1)`; its determinant is the constant `-25` (Bareiss), so
  `det H_ch = -25/(25 L^2 Lam1(L)^2) = -(1 - L)^2/(L^2 (1 - L^5)^2)`. Product: `det H_vv det H_ch = 1/(L^2 (1 - L^5)^2) = det H_7`.
  **`det H_ch = -1/(L^2 Lam1(L)^2)` vanishes only at `L = 0, zeta, .., zeta^4`; it is `-1/25` at `L = 1`.**
* VERIFIED (`ch_hess.py`): all entries, both orders and the direct determinant agree exactly at eight rational `L`; the two
  constant determinants (5 and -25) are polynomial identities in `L`. Controls: `q -> -q` moves the critical point.

## 3.3 Step (a), second half: order (a) computes c(z; L)

Apply (Fub) (L1b) to `Omega_7` with `v = (y, x_0..x_4)` and `w = sigma`. Here `v_c(sigma) = (5 sigma, zeta^j - sigma)` and
`phi(sigma) = Phi(v_c(sigma), sigma) = S_0(sigma)` (s.1.6). At fixed `sigma` the inner phase is
`lambda sigma + [5 sigma log y - y] + Sum_j [(zeta^j - sigma) log x_j - x_j]` and the inner amplitude is `z^-1 . 1 . prod_j x_j^-1`: a
product (L0) of six one-variable integrands. By (L3), with `a = 5 sigma`, `M = 1` (measure `dy`) and `a = zeta^j - sigma`, `M = 0`
(measure `dx_j/x_j`):

    B(sigma, z) = z^-1 exp(Sum_k b_k (z/(5 sigma))^(2k-1)) prod_j (zeta^j - sigma)^-1 exp(Sum_k b_k (z/(zeta^j - sigma))^(2k-1))
                = z^-1 prod_j (zeta^j - sigma)^-1 exp(Sum_k z^(2k-1) S_2k(sigma)).

With `det H_vv(sigma) = 1/(5 sigma prod_j (zeta^j - sigma))`, `beta(sigma) = [5 sigma prod_j (zeta^j - sigma)/(5L prod_j (zeta^j - L))]^(1/2)`,
and `beta B = z^-1 (5L (1 - L^5))^(-1/2) (5 sigma)^(1/2) prod_j (zeta^j - sigma)^(-1/2) exp(Sum_k z^(2k-1) S_2k) = const . z^-1 A_St(sigma, z)`,
since `e^(S_1) = (5 sigma)^(1/2) prod_j (zeta^j - sigma)^(-1/2)`. Hence by (Fub) and (F2)

    N[Phi, z^-1 prod_j x_j^-1; p(L)] = const . z^-1 N[S_0, A_St; L],      so      Nhat_7(z; L) = c(z; L).

*Bookkeeping in symbol form.* The inner SP is `(2 pi z)^3 det(-H_vv)^(-1/2) e^(phi/z) B`; the six factors `z^(a/z + M)` of (L3)
multiply to `z^(+1)` identically in `sigma` (`ch_zpow.py` (ii)) and cancel the `z^-1` of `Omega_7`; the result is exactly the
Stirling-expanded summand `(2 pi)^3 z^2 exp(S_0/z + S_1 + Sum z^(2k-1) S_2k)` of s.1.6. So "order (a) = full-Stirling SP" holds as
symbols, not only for normalised series. Valid for every `L` in `Lam` (it needs `5L != 0` and `zeta^j - L != 0` for all `j`, in
particular `L != 1`). VERIFIED directly: the 7-variable Wick series of `Omega_7` (full 7x7 covariance) equals `c(z; L)` exactly over
`Q(zeta)` through `z^4` at `L = 2`, and mod two primes through `z^5` at `L = 3` (`ch_chart.py`).

## 3.4 Step (b): order (b), the (sigma, x_0)-integral by (L2), and the chart integrand

Write `Phi = sigma G + f` with `G := lambda + 5 log y - Sum_(j=0..4) log x_j` and `f := -y + Sum_j (-x_j + zeta^j log x_j)`; the amplitude
`z^-1 (x_0 prod_(j>=1) x_j)^-1` does not involve `sigma`. Apply (L2) with `w = (y, x_1, .., x_4)`: `d_(x_0) G = -1/x_0 != 0`,
* `X_0(w) = exp(lambda + 5 log y - Sum_(j>=1) log x_j) = -q y^5/prod_(j>=1) x_j`, `sigma_c(w) = 1 - X_0(w)`;
* inner series EXACTLY `N = z^-1 (X_0(w) prod_(j>=1) x_j)^-1`; `det H_vv(w) = -X_0(w)^-2`, so `beta(w) = X_0(w)/X_0(w_*)` with
  `w_* = p_ch(L)` and `X_0(w_*) = 1 - L`; hence `beta N = z^-1 (1 - L)^-1 prod_(j>=1) x_j^-1`;
* inner critical value `f(X_0(w), w) = lambda + psi(w)`, where (using `log X_0 = lambda + 5 log y - Sum_(j>=1) log x_j` and
  `-X_0 = q y^5/prod x_j`)

      psi(y, x_1..x_4) := 5 log y - y + Sum_(j=1..4) (-x_j + (zeta^j - 1) log x_j) + q y^5/(x_1 x_2 x_3 x_4).

**The chart integrand** is `Omega_ch := e^(psi/z) dy prod_(j=1..4) dx_j/x_j`. In symbol form the `(sigma, x_0)`-SP of `Omega_7` is
`2 pi z (-X_0^-2)^(-1/2) e^((lambda + psi)/z) z^-1 (X_0 prod x_j)^-1 = 2 pi kappa (-q)^(1/z) e^(psi/z) prod_(j>=1) x_j^-1`, `kappa = +-i`:
the factor `X_0` of the Gaussian factor (the "Jacobian") cancels the `1/x_0` of `dx_0/x_0`. By (Fub) and (F2):

    N[Phi, z^-1 prod_j x_j^-1; p(L)] = z^-1 (1 - L)^-1 N[psi, prod_(j>=1) x_j^-1; p_ch(L)],      so      Nhat_7(z; L) = Nch(z; L).

*The chart critical point.* `p_ch(L) = (5L, zeta^1 - L, .., zeta^4 - L)`, where `K := q y^5/prod x_j = L - 1`:
`d_y psi = 5/y - 1 + 5K/y` and `d_(x_j) psi = -1 + (zeta^j - 1 - K)/x_j` vanish there. The phase `psi` depends polynomially on `q`,
and `q(L)`, `p_ch(L)` are holomorphic at `L = 1` (where `q = 0`, `K = 0`, `p_ch(1) = (5, zeta^j - 1)`); the Hessian is invertible
for `L` outside `{0, zeta, .., zeta^4}` (s.3.2). Order (b) itself needs `L != 1` (`X_0 = 1 - L != 0`).
VERIFIED: the 5-variable Wick series of `Omega_ch` equals `c(z; L)` exactly over `Q(zeta)` through `z^3` at `L = 2, -2/3` and mod two
primes through `z^6` at `L = 3, 5/4, -1/2` (`ch_chart.py`; the Stage-4 critic: through `z^5` against `R_n(L)`). Controls: without the
Jacobian `X_0` it fails at `z^1`; with `dy/y` at `z^2`.

## 3.5 Step (c): Fubini, comparison of the normalised series, continuation to L = 1

* **Comparison.** Both orders compute the SAME object `N[Phi, z^-1 prod_j x_j^-1; p(L)]` (by (Fub), an identity of series `N`, so
  no prefactor has to be tracked by hand). Hence, for every `L` in `Lam`: **`c(z; L) = Nhat_7(z; L) = Nch(z; L)`** (3.3, 3.4).
  In the language of symbols: both orders produce `(2 pi z)^(7/2) det(-H_7)^(-1/2) e^(Phi(p)/z) z^-1 N`, with the same critical
  value (`S_0(L) = Phi(p(L)) = lambda + psi(p_ch(L))`), the same power `z^(5/2)` and non-zero constants (`det(-H_7)^(-1/2)` with the
  branch convention of (Fub)); every non-series prefactor is a constant times a power of `z`, and the normalised series agree.
* **Rationality of `Nch`** ((L1c)(i) with parameter `L`). The Taylor coefficients of `psi` at `p_ch(L)` are `5 (-1)^(k+1)/(k (5L)^k)`,
  `(zeta^j - 1)(-1)^(k+1)/(k (zeta^j - L)^k)`, and `(L - 1)` times the Taylor coefficients of `(1 + eta/(5L))^5 prod_j (1 + xi_j/(zeta^j - L))^-1`;
  the amplitude has those of `prod_j (1 + xi_j/(zeta^j - L))^-1`; the covariance is `-adj(H_ch)/det H_ch` with
  `1/det H_ch = -L^2 Lam1(L)^2` a polynomial. So `Nch_n in Q(zeta)[L, L^-1, (zeta^j - L)^-1 : j = 1..4]`: **regular at `L = 1`**
  (Theorem Gamma (ii)).
* **Continuation.** `c_n` (s.1.6) and `Nch_n` are rational functions of `L` that coincide on `Lam`; so they are equal as rational
  functions. `Nch_n` has no pole at `L = 1`, hence neither has `c_n`, and `c_n(1) = Nch_n(1)`. No a priori polynomiality of `c_n` and
  no uniformity in `q` is used: the continuation acts on coefficients (rational functions), not on asymptotic expansions.

## 3.6 Step (d): q = 0, the factorisation and the exact match with h_(w_0)

At `L = 1`: `q = 0`, `K = 0`, and the chart phase is `psi_0 = [5 log y - y] + Sum_(j=1..4) [(zeta^j - 1) log x_j - x_j]`, critical point
`p_ch(1) = (5, zeta^j - 1)`, Hessian `diag(-1/5, -1/(zeta^1 - 1), .., -1/(zeta^4 - 1))` (`ch_hess.py`), amplitude `prod_(j>=1) x_j^-1`
with measure `dy`. By (L1c)(i), `Nch(z; 1)` is the normalised series of this integrand. By (L0) it is the product of five one-variable
series, and by (L3) (`a = 5`, `M = 1`, measure `dy`; `a = zeta^j - 1`, `M = 0`, measure `dx_j/x_j`; `Poch_1 = Poch_0 = 1`):

    Nch(z; 1) = exp( Sum_k b_k z^(2k-1) T_k ),   T_k := 5^(1-2k) + Sum_(j=1..4) (zeta^j - 1)^(1-2k).

These are the Euler integrals of `Gamma(1 + 5/z)` and `Gamma((zeta^j - 1)/z)`. Since `1 - 2k` is odd,
`-Sum_(a in w_0) a^(1-2k) = -(-5)^(1-2k) - Sum_j (1 - zeta^j)^(1-2k) = T_k`. Therefore `Nch(z; 1) = exp(-Sum_k b_k z^(2k-1) Sum_(a in w_0) a^(1-2k))
= h_(w_0)(z)`, with `b_k = B_2k/(2k(2k-1))` as in QRR. With 3.5: **`c(z; 1) = h(z)`, Theorem Gamma (iii).** QED
* *The twist weight.* The factor `Gamma(1 + 5/z)` contributes `+5^(1-2k)`, i.e. the weight `-5` of `w_0`. The weight `+5` would need
  `Gamma(1 - 5/z)` and changes the `z^1` coefficient of `log h` from `b_1 (1/5 - 2) = -3/20 = h_1` to `b_1(-1/5 - 2) = -11/60`.
* VERIFIED: `T_k = -Sum_(a in w_0) a^(1-2k)` exactly in `Q(zeta)` for `k <= 30`, and `exp(Sum b_k T_k z^(2k-1)) = h` (Newton sums) through
  `z^59`, the product of the five normalised Stirling series in `Q(zeta)[[z]]` through `z^25`, `h_1 = -3/20`, `h_5 = -257071/8064000000`
  (`ch_stirling.py`); the chart SP at `L = 1` equals `h` exactly through `z^7` and mod two primes through `z^6` (`ch_chart.py`).
  Controls fail at `z^1`: twist `+5`; `1 - zeta^j` for `zeta^j - 1`; `Gamma(1 + 5/z)` dropped; `b_k -> -b_k`.

## 3.7 Step (e): consistency with the I-function, coefficientwise in q (PROVED given (L1c)(iii); not used above)

Near `L = 1` write the chart phase as `psi_q = psi_0 + q K_1`, `K_1 := y^5/(x_1 x_2 x_3 x_4)`, with the `q`-independent amplitude
`A = prod_(j>=1) x_j^-1`, and let `p_0 = p_ch(1)`. By (L1c)(iii), coefficientwise in `q`,
`e^(-psi_0(p_0)/z) SP_(p_ch(q))[e^(psi_q/z) A] = Sum_d (q^d/(d! z^d)) e^(-psi_0(p_0)/z) SP_(p_0)[e^(psi_0/z) K_1^d A]`.
The `d`-th term factorises (L0) into Euler integrals: `y^(5d)` with `dy` (`a = 5`, `M = 5d + 1`) and `x_j^(-d)` with `dx_j/x_j`
(`a = zeta^j - 1`, `M = -d`, the negative case of (L3)). With `Gamma_St(w + 1) = w Gamma_St(w)` the ratio to the `d = 0` term is

    (q/z)^d (1/d!) z^(5d) prod_(k=1..5d) (5/z + k) . prod_(j>=1) z^(-d) / prod_(i=1..d) ((zeta^j - 1)/z - i)
        = q^d (1/d!) prod_(k=1..5d) (5/z + k) / prod_(j>=1) prod_(i=1..d) ((1 - zeta^j)/z + i) =: q^d I_d(z)

(the powers of `z` cancel, `ch_zpow.py` (iv); the signs `(-1)^(4d) = 1`). `I_d` is the `d`-th Frobenius coefficient of the colour-0
I-function, `z q^(1/z) Sum_d q^d prod_(k<=5d)(5 + kz)/prod_j prod_(k<=d)(1 - zeta^j + kz)` (FORMAL_QUINTIC_genus2 s.1, `lambda_i = 1`).
The `d = 0` term has normalised series `h(z)` (3.6). The prefactors of the moving critical point combine to
`[det(-H_ch(q))/det(-H_ch(0))]^(-1/2) A(p_ch(q))/A(p_0) = (L Lam1(L)/5)(5/Lam1(L)) = L`. With 3.5 this gives, in `Q((z))[[q]]`,

    L e^(Dpsi(q)/z) c(z; L(q)) = h(z) Sum_(d>=0) q^d I_d(z),        Dpsi(q) := psi_q(p_ch(q)) - psi_0(p_0),

with `L(q) = (1 - 5^5 q)^(-1/5)` and `D_q Dpsi = q (d_q psi_q)(p_ch) = K = L - 1`, so `log q + Dpsi` is (up to a constant) the colour-0
canonical coordinate `u_0` (`D_q u_0 = L`). So the chart series is `h(z)` times the I-function coefficientwise in `q`; at `q^0` this is
step (d) again; it is the form (Gamma3) of SYN4 s.1.3 (identification with arXiv:0710.0889 eq. (5) not re-checked here).
VERIFIED (`ch_eexp.py 6 10`): all 98 non-zero coefficients `q^d z^n`, `d <= 6`, `-d <= n <= 10`, exactly over Q, with `c_n(L)`
from the sigma-form (fitted as polynomials in `L`, fits confirmed on further points), `h` from Newton sums and `I_d` from the product;
`D_q Dpsi = L - 1` through `q^6`. Controls fail: the factor `L` dropped, `e^(Dpsi/z)` dropped (first at `q^1`), twist `+5` in `h`
(`q^0 z^1`). The `(L1c)(iii)` mechanism alone is VERIFIED in a one-variable instance (`ch_lemmas.py`).

## 3.8 Relation to Givental's chart computation (arXiv:math/0108100 v2, s.9-10, read 2026-11-14)

* *What is there.* Thm 9.1: for a torus action with isolated fixed points, the R-matrix of the equivariant theory, normalised by
  `R(z)|_(Q=0) = exp diag(b_1, .., b_N)`, `b_i(z) = Sum_k N^(i)_(2k-1) (B_2k/2k) z^(2k-1)/(2k-1)`, `N^(i)_l` the Newton sums of the inverse
  characters of the COTANGENT space at the fixed point `w_i`, gives the total ancestor potential. Section 10 (pp. 13-15): Thm 10.1, the
  oscillating integral `I = Int_(Gamma in {x_0..x_n = Q}) e^((x_0+..+x_n)/z) x_0^(lambda_0/z)..x_n^(lambda_n/z) dx_0..dx_n/d(x_0..x_n)` solves
  the equivariant quantum differential equation of `CP^n`; Example 10.3, the formal stationary-phase expansion; then the chart
  `(x_1, .., x_n)` with `x_0 = Q/(x_1..x_n)`, which still contains the critical point `P^(0)` when `Q = 0`; at `Q = 0` the integral
  factorises into one-dimensional integrals; Lemma 10.4 (Stirling); the result `ln R^00 = Sum_k (B_2k/2k) z^(2k-1)/(2k-1) Sum_(a != 0) (lambda_a - lambda_0)^(1-2k)`,
  which "coincides with `b_0(z)` in 9.1". The text calls itself an outline (details elsewhere); it states no Fubini or delta lemma.
* *What is different here.* (1) The theory is twisted by `O(5)`: the extra variable `y` (the factor `Gamma(1 + 5s)`), whose Stirling
  factor supplies the twist weight `-5` of `w_0`. (2) The object of (Gamma) is the ONE-variable full-Stirling series `c(z; L)`; s.3.3
  connects it with the 7-variable integral by (Fub). (3) Givental integrates over the constraint surface; here the constraint
  `G = log(-q y^5/prod x_j) = 0` is produced by the Mellin-Barnes variable `sigma` through (L2): `sigma` is exactly the Lagrange
  multiplier `P` of Givental's p. 14 ("constrained extrema ... with the Lagrange multiplier P"). (4) The continuation to `q = 0` is
  done on coefficients (rational functions of `L`), not on asymptotic expansions.

## 3.9 The Z/5 covariance (Theorem Gamma (iv), s.0.2) (PROVED)

* From `P_k(zeta L) = Sum_j (zeta^j - zeta L)^-k = zeta^-k P_k(L)` and the Taylor data of s.1.6:
  `[t^(k+1)] (S_0(zeta L + t) - S_0(zeta L)) = zeta^-k [t^(k+1)](S_0(L + t) - S_0(L))`, `[t^k] (S_1(zeta L + t) - S_1(zeta L)) = zeta^-k [t^k](...)_L`,
  `[t^r] S_2k(zeta L + t) = zeta^(1-2k-r) [t^r] S_2k(L + t)`.
* Substitute `t = zeta t'` ((L1a), constant Jacobian `zeta`). In `t'` the phase is `zeta` times the phase at `L`, the normalised
  `e^(S_1)`-part is the one at `L`, and `z^(2k-1) S_2k` becomes `(zeta^-1 z)^(2k-1) S_2k(L + t')`. So the integrand at `zeta L` is
  `e^(Phi_L/(zeta^-1 z)) A_L(t', zeta^-1 z)` up to a constant, and by (F3): `c(z; zeta L) = c(zeta^-1 z; L)`, i.e.
  `c_n(zeta L) = zeta^-n c_n(L)`; by iteration `c(z; zeta^i L) = c(zeta^-i z; L)`.
* *Colour `i`.* For the parameter `q(L) = q(zeta^i L)` the critical point `sigma = zeta^i L` is the colour-0 point of the parameter
  `zeta^i L`. Hence `c^(i)(z; L) = c(z; zeta^i L) = c(zeta^-i z; L)`.
* *Consequences.* `c_n` is regular at `L = zeta^i` with `c(z; zeta^i) = h(zeta^-i z)`; and `h_(w_i)(z) = h_(w_0)(zeta^-i z)`
  because `Sum_(a in w_i) a^(1-2k) = zeta^(i(1-2k)) Sum_(a in w_0) a^(1-2k)`. By Theorem Gamma (i) the only possible poles of `c_n` are
  `L = 0` and `L^5 = 1`; none at `L^5 = 1`, so `c_n in Q[L, L^-1]`, with exponents `= 4n mod 5` (from `c_n(zeta L) = zeta^-n c_n(L)`).
* VERIFIED (`ch_sigma.py 10`): `c_n(zeta^i L) = zeta^(-in) c_n(L)` exactly in `Q(zeta)` for `i = 1..4`, `n <= 10`, `L = 2, -2/3`;
  the control `zeta^(+in)` fails.

## 3.10 The hatted form: `c_n(L) = L^(4n) ph^SP_n(L^-5)` (PROVED)

POLY s.5 uses `x = 1/sigma = X/L` and `z = L^-4 w`. The map `sigma = L/X` is analytic near `X = 1` ((L1a); its Jacobian is removed by
normalisation), and `S_0(L/X)/z = L^4 S_0(L/X)/w`, so by (F3) (applied to the phase `L^4 S_0(L/X)` with `mu = L^-4`) the normalised series in `z` of the sigma-form equals the
normalised series in `w` of the pulled-back form, evaluated at `w = L^4 z`. The pulled-back form is POLY's hatted integrand up to a
constant (s.1.7), and `Phi^SP(w; u) = Sum_n w^n ph^SP_n(u)` is its normalised series. By (L1c)(i) both sides are rational in the
parameter; `ph^SP_n` is a polynomial (Theorem P(a)). Hence `Sum_n c_n(L) z^n = Sum_n ph^SP_n(L^-5) L^(4n) z^n` as rational functions,
and at `L = 1`: `c_n(1) = ph^SP_n(1)`. (Numerically this is the critic's E3: `c_5(2)` from `u = 1/32`, reproduced in `ch_sigma.py`.)

## 3.11 Proof of Theorem D

By Theorem Gamma (iii) and s.3.10, `ph^SP_n(1) = h_n` for every `n`, i.e. Phi^SP(w; 1) = h(w) (poly's "(A2)"). By Theorem P(b) (R5),
`Phi^can = Phi^SP`, so `ph^can_n(1) = h_n` for every `n`. By the connection form (EVID s.2.4; R3, which rests on Lemma N, R2), `D = 1`.
In particular `D_0 = 1` (and `D_i(z) = D_0(zeta^-i z) = 1`). Hence (F_g) for every `g >= 2`. QED
(This is R7 of SYN4, `D = 1 <=> (Gamma)`, read from right to left, with (Gamma) now proved.)
