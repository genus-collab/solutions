# Uniform reduction of the formal quintic at the conifold point: why R~(zhat; u) = G(u, zhat) M_0(tau(u) zhat)

Date 2026-11-10. Folder `natural_form/theory/` (new files only; nothing else in the project was modified). Conventions of
OU_EXTENSION_genus3.md s.1 and CONIFOLD_SYMBOL_genus2.md: formal quintic of LP18, `zeta` a primitive 5th root of unity, `u = 1 - 5^5 q`,
`s = u^(1/5) = 1/L`, `zhat = z/s^4`, `c = (0, -1, 3, 2, 1)`, hatted rows `Qhat_a = s^(1-c_a) Q_a`, conifold frame `(alpha, beta) = (0, 0)`
unless stated, flat index `a` in Z/5, `E_ab` the matrix units. Status labels: **PROVED** (argument given here), **VERIFIED** (exact
computation, orders and primes stated), **CONJECTURAL**, **SKETCH**. Everything is AI-checked only; no human has read it.

## 0. Summary

1. **Cyclic form (VERIFIED).** In the flat frame the rows of the formal quintic satisfy the cyclic system
   `th Y_a = b_a(u) Y_{a+1}` (`a` in Z/5, `th = z d/du`) with
   `b = (t_c', -1/(u varpi_0^2 (1-u)^3 t_c'^2), t_c', -varpi_0/(1-u), -varpi_0/(1-u))`, `prod b_a = -1/(u (1-u)^5)`: four entries are units
   at `u = 0` and `b_1` (the Yukawa coupling `H * H = b_1 H^2`) has a simple pole. Checked on the s-series rows including the closing
   relation `4 -> 0`, all colours, `z^1..z^6`, both primes (s.1).
2. **Scaling limit = an LG model (PROVED).** For `u -> 0` at fixed `zhat` the system becomes `z d/dv Y_b = bmod_b Y_{b+1}`, `bmod = (1, -1/v, 1,
   -1, -1)`, which is exactly the Gauss-Manin system of `lambda(x) = -(x^4/4 + v/x)` with primitive form `dx` in the basis
   `(dx, -dx/x, -v dx/x^2, x^2 dx, x dx) <-> (1, H, H^2, H^3, H^4)`; scalar equation `(v th^5 + 3 z th^4 + 1) Y_0 = 0`. (The same model as
   `x^4 + u/x` of the hand check after constant rescalings, `v -> c v` with `c^4 = -1/4` in the scalar equation. In the conventions used
   here the dictionary is `H^a <-> x^(-a) dx` modulo `x^5 = v`, i.e. the five flat directions `d_c, d_u, d_d, d_a, d_b` of the Hurwitz
   space of `x^4 + a x^2 + b x + c + u/(x - d)`; in particular `H <-> dx/x`, the derivative of `lambda` in the residue. It is NOT
   `H^a <-> x^a dx`: with the colour `i` attached to the critical value `-(5/4) zeta^i sigma^4`, i.e. `x_c = zeta^-i sigma`, the form
   `x^a dx` has the phase `zeta^(i(1-a))` and the class `H^b` the phase `zeta^(i(b+1))`.)
3. **The symbol is the model R-matrix (PROVED given the purity of the hatted rows; VERIFIED all five rows, all colours, `zhat^0..zhat^6`,
   both primes).** `M_0` is the homogeneous WKB solution of the model system: `R^_b^(i)(zhat) = zeta^(i(b+1)) F_b(zeta^-i zhat)` with `F_b`
   computed over Q from the model alone; the identity row is the stationary phase of `Int exp(-(x^4/4 + v/x)/z) dx` (checked through
   `t^10`). By a Mellin-Barnes transform the model integral is `Int Gamma(s) Gamma((s+1)/4) (...)^s ds` (a Meijer function of type
   `G^{5,0}_{0,5}`, a `0F4` equation in `v^4/z^5`): the direct generalisation of the Gamma integral behind the Hodge factor (s.2, s.7).
4. **Theorem U (PROVED, s.3).** Two cyclic systems `th Y_a = b_a Y_{a+1}`, `th Y_a = bt_a Y_{a+1}` whose coefficients are units at `u = 0`
   except for a simple pole of `b_1`, `bt_1`, and with `prod b_a = prod bt_a`, are related by a formal gauge `P(u, z) = Sum_k z^k P_k(u)`
   whose coefficients are REGULAR at `u = 0` (power series), `P_0` diagonal and invertible; `P` is unique up to `P -> P c(z)`, `c` a
   constant scalar series. Applied to the quintic and to the model pulled back by `v(u)`: the condition `prod b = prod bt` (equality of
   the canonical 1-forms `du_i`) is `v^(-1/5) dv = u^(-1/5) du/(1-u)`, i.e. `v = u tau^(-5/4)` with `tau = 1/2F1(1, 4/5; 9/5; u)`. This is why
   the natural variable is `z` over the canonical-coordinate difference: `z/v^(4/5) = tau zhat`, and no other change of variable admits
   even a leading-order gauge.
5. **The factorisation (PROVED modulo the purity (H) of the hatted rows, which is VERIFIED mod `z^7`, `u^6`).** `R~(zhat; u) = G(u, zhat)
   M_0(tau zhat)` with `G = Psi(u)^-1 Ghat Psi(0)`, `Ghat_ab = tau^((1-c_b)/4) Sum_k zhat^k u^(mu(b,k)) p_{k,b}(u) delta_{a,b+k}` the hatted
   form of `P`. **Independent check (VERIFIED):** `P` computed over Q from the periods alone (no R-matrix input) reproduces the project's
   `G` in all 25 entries, `zhat^0..zhat^6`, `u^0..u^5`, at both primes, with ONE constant `c_5 = 43397/807975000` of the freedom
   `P -> P (1 + c_5 z^5)`; and that constant is exactly the one for which `P_11(0, z) = P_22(0, z) = 1 + O(z^10)` (canonical normalisation,
   s.4).
6. **Degrees (PROVED as upper bounds, VERIFIED as equalities for n <= 12).** `deg_zhat [u^n] G = floor((5n + 2)/4) = n + floor((n + 2)/4)`:
   `0, 1, 3, 4, 5, 6, 8, 9, 10, 11, 13, 14, 15, ...`. The a priori grading bound is `floor((5n + 4)/4)`; the improvement comes from the
   residue relations (R1)-(R2) of s.5.
7. **The `H (x) H` component (PROVED).** `P(0, z)` commutes with the residue `E_12` of the connection, so column 1 and row 2 of `P(0, z)`
   vanish off the diagonal; the entry `(2, 1)` (the flat-frame term `z^k H (x) H`, `k = 1 mod 5`) vanishes to SECOND order at `u = 0`. Hence no
   `H (x) H` in `G_0`, `G_1`: forced, not accidental. Moreover the block `rows (2, 3) x columns (0, 1)` is, to leading order, the rank-one
   matrix `pi(z) (u/z, 1)^T (1, -u/z)` with `pi(z) = P_30(0, z) = (127/3150) z^3 + ...`, i.e. in hatted variables the `H`-part of `G_2` is
   `alpha zhat (H - zhat 1) (x) (H + zhat 1)`, `alpha = -127/15750`: two insertions of the Euler-dilaton vector.
8. **Graph sum (VERIFIED, production runs of today: frame (0, 0) at both primes, frame (1/3, -2/7) at 2^89-1).** With `M = M_0(tau(u) zhat)`
   and nothing else, all 41 edge graphs of
   genus 3 are `O(u^4)` and all six edge graphs of genus 2 are `O(u^2)`: regular; the vertex term carries `(1/25200) t_c^-4` resp.
   `(1/1200) t_c^-2`. Nothing has to be added to the vertex factor through genus 3. The first-order correction is the point `D_1 = -X_1` of
   the kinematic plane Pi of OU_EXTENSION_genus3 s.6d, with `X_1` explicit in the periods: its `H^3 (x) H^4` coefficient `-(2 gamma_g - 4/3)/5 =
   -1/75` is frame independent and its string-direction coefficient `-8/225 + w_1/5` is not (formula confirmed in the second frame, s.5).
   Partial mechanism (PROVED given the primary-sector inputs, VERIFIED both primes):
   the Euler-dilaton part is invisible at first order on every one-edge graph with a vertex of genus >= 2. Negative result: it is NOT an
   exact symmetry (genus-1 anomaly `-1/8`).
9. **Conjecture N (s.6).** For all `g >= 2`, with `M = M_0(tau zhat)` (the LG model `-(x^4/4 + v/x)` in its natural variable) and
   `R' = Psi(u) G` (the regular gauge), every stable graph with edges is regular at `u = 0`. Supported at `g = 2, 3`; open points and the
   genus-4 test are listed in s.6.5.

## 1. The two systems

**Quintic.** fq_rows.py: `X_{a+1} = (1/(C_1 kappa_a)) (z D + zeta^i L) X_a` (`a = 0..3`), `kappa = (1, L^5/(C_0^2 C_1^3), 1, C_0/C_1, C_0/C_1)`,
`D = q d/dq = -(1-u) d/du`, PF operator `M^5 - 1 - q prod_{k=1..5}(5M + kz)`, `M = zD`. With `Y_a = e^{u_i/z} X_a` (`D u_i = zeta^i L`) and
`c = (C_1, L^5/(C_0^2 C_1^2), C_1, C_0, C_0)` (`prod c_a = L^5`):

    z D Y_a = c_a Y_{a+1}   (a in Z/5),      i.e.   th Y_a = b_a Y_{a+1},   b_a = -c_a/(1 - u),   th = z d/du.

The closing relation `a = 4 -> 0` is the fifth-order equation (`H^5 = 1` in the formal quintic); it is not used by fq_rows and was
checked here (ur_model_check.py (A): all colours, `z^1..z^6`, >= 32 s-orders, both primes). In the conifold frame `C_0 -> varpi_0^cf`, `C_1 ->
D t_c = -(1-u) t_c'`, which gives the `b_a` of s.0 item 1; `varpi_0 = 1 + O(u^3)`, `t_c = u + (7/10) u^2 + ...` in frame `(0, 0)`, and the
hypotheses "units, one simple pole" hold in every conifold frame `(alpha, beta)`. The normalised rows `Q_a^{(i)} = zeta^i X_a` satisfy the
same system. In matrix form `th Y = A Y`, `A = Sum_a b_a E_{a,a+1}`, a weighted cyclic shift; its eigenvalues are the `du_i/du = -zeta^i
u^(-1/5)/(1-u)`, so `u_i - u_i^c = -zeta^i (5/4) u^(4/5) 2F1(1, 4/5; 9/5; u) = -zeta^i (5/4) u^(4/5)/tau(u)`.

**Scaling limit (PROVED).** Put `z = u^(4/5) zhat`. Then `th ~ u^(-1/5)`, and in `P = u M^5 - (1-u)[3 z M^4 + (17/5) z^2 M^3 + (9/5) z^3
M^2 + (274/625) z^4 M + (24/625) z^5] - 1`, `M = -(1-u) th`, the terms `u th^5`, `3 z th^4`, `1` are of order `u^0` and every other term
(including all commutators `[th, u] = z`) is smaller by integer powers of `u`. So the limit equation is `(u th^5 + 3 z th^4 + 1) Y_0 = 0`,
and for the system: `b -> (1, -1/u, 1, -1, -1)` (`t_c'(0) = 1`, `varpi_0(0) = 1`). Local exponents at `u = 0`: `0, 1, 1, 2, 3`.

**Model.** For `I_m = Int x^m exp(lambda/z) dx`, `lambda = -(x^4/4 + v/x)`: `z d/dv I_m = -I_{m-1}` and, from `Int d(x^m e^{lambda/z}) = 0`,
`z m I_{m-1} - I_{m+3} + v I_{m-2} = 0`. Hence with `(Y_0, .., Y_4) = (I_0, -I_{-1}, -v I_{-2} = -I_3, I_2, I_1)`:

    z d/dv Y_b = bmod_b Y_{b+1},   bmod = (1, -1/v, 1, -1, -1),     (v th^5 + 3 z th^4 + 1) Y_0 = 0     (PROVED, three lines).

So the limit of the quintic system IS the Gauss-Manin system of this LG model (a point of the Hurwitz space `M_{0;3,0}` of functions
`x^4 + a x^2 + b x + c + u/(x - d)`, on the locus `a = b = d = 0`, with `c` irrelevant), in the frame of the five flat directions:
`1 <-> d_c lambda`, `H <-> d_u lambda = 1/x`, `H^2 <-> d_d lambda = u/x^2`, `H^3 <-> x^2`, `H^4 <-> x` (up to constants). The model is
homogeneous: `x -> k x`, `v -> k^5 v`, `z -> k^4 z`; critical points `x^5 = v`, critical values `-(5/4) x_c^4`, `Delta_i = lambda''(x_i)
~ x_i^2`, so the weights `Y_b ~ sigma^(c_b - 1) x (series in z/sigma^4)`, `sigma = v^(1/5)`, are exactly the hatted weights
`c = (0, -1, 3, 2, 1)` of the quintic (this is where `c` comes from).

## 2. The conifold symbol is the R-matrix of the model

**Proposition S.** Let `Y_b = e^{U/z} sigma^(c_b - 1) F_b(t)`, `t = z/sigma^4`, `U = -(5/4) zeta^i sigma^4`. The model system is equivalent to

    F_1 = T_{-1} F_0,  F_2 = -T_{-2} F_1,  F_3 = T_2 F_2,  F_4 = -T_1 F_3,  F_0 = -T_0 F_4,     T_k = -zeta^i + (t/5)(k - 4 t d/dt),

and has, for each colour, a unique solution with `F_0(0) = 1` (the recursion for `[t^n] F_0` has the non-zero pivot `-4n`). Its rows are the
conifold symbol: `R^_b^(i)(zhat) = zeta^(i(b+1)) F_b(zeta^-i zhat)` (`F_b` of colour 0), `F_b(0) = (1, -1, -1, 1, 1)`,

    F_0 = 1 - (3/20) t + (9/800) t^2 + (11571/80000) t^3 - (2415609/6400000) t^4 + (47735919/128000000) t^5 + (52989974037/25600000000) t^6
          - (7749624763881/512000000000) t^7 + ...

*Proof.* The hatted rows are power series in `u` at fixed `zhat` (hypothesis (H), s.4; asserted by `HatRowsU` mod `z^7`, `u^6`), their `u^0`
terms satisfy the limit system (s.1), and a WKB solution of the model system that is homogeneous (a function of `t` alone times
`sigma^(c_b-1)`) is unique up to a constant: the general WKB column is the homogeneous one times a constant series `d(z) = Sum d_k
sigma^(4k) t^k`, which is homogeneous only if `d_k = 0` for `k >= 1`. QED (PROVED given (H)).
*Verification* (ur_model_check.py (B), log ur_model_check.log): all 5 x 5 x 7 numbers `R^_{b,n}^(i)`, both primes: PASS; the stationary
phase of `Int exp(-(x^4/4 + v/x)/z) dx` at the critical point `x = 1`, `v = 1` equals `F_0` through `t^10` (exact, Q): PASS.

**Mellin-Barnes form (PROVED, elementary; not machine-checked beyond the ODE).** `Int_0^inf exp(-xi^4/4 - k/xi) d xi = (1/2 pi i) Int
Gamma(s) Gamma((s+1)/4) 4^((s-3)/4) k^(-s) ds`, `k = v/z^(5/4)`; by the multiplication formula the integrand has five Gamma factors
`Gamma(s/4) Gamma((s+1)/4)^2 Gamma((s+2)/4) Gamma((s+3)/4)`: the indicial polynomial `th(th-1)^2(th-2)(th-3)` of the model equation. `M_0` is the
large-argument (exponentially small/large branches, one per colour) expansion of this function and of its four companions.

## 3. Theorem U: the regular gauge

Let `O = C[[u]]` (or convergent series), `K = O[1/u]`. A *cyclic system with one simple pole* is `th Y = A Y`, `A = Sum_a b_a E_{a,a+1}`
(`a` in Z/5), with `b_a` in `O^x` for `a != 1` and `u b_1` in `O^x`.

**Theorem U.** Let `A`, `At` be two such systems with `prod_a b_a = prod_a bt_a`. Then there is `P = Sum_{k>=0} z^k P_k(u)`, `P_k` in
`gl_5(O)`, `P_0` diagonal with unit entries, such that `Y = P Yt` maps solutions to solutions, i.e. `z P' = A P - P At`. Every such `P`
with coefficients in `gl_5(O)` and `P_0` invertible is obtained from one of them by `P -> P c(z)`, `c` in `C[[z]]^x`. If one restricts to
the Z/5-graded form `P_k = Sum_b p_{k,b} E_{b+k,b}`, the freedom is `c(z^5)`.

*Proof.* (i) *Reduction.* `D_0 = diag(d_a)` with `d_{a+1}/d_a = bt_a/b_a` is consistent around the cycle (equal products), has unit
entries, and `A D_0 = D_0 At`. Writing `P = D_0 Ph` the equation becomes `z Ph' = [At, Ph] - z Dl Ph`, `Dl = D_0^-1 D_0'` diagonal in
`gl_5(O)`; order `z^(k+1)`: `[At, Ph_{k+1}] = Ph_k' + Dl Ph_k`, order `z^0`: `[At, Ph_0] = 0`.
(ii) *The commutant.* `y := At^-1 = Sum_a bt_a^-1 E_{a+1,a}` is holomorphic (`1/bt_1 = O(u)`), `y^5 = 1/B`, `B = prod bt_a = unit/u`. `At` is
regular semisimple over `K`, so its commutant in `gl_5(K)` is `K[y]` and `gl_5(K) = K[y] + im(ad At)` (direct). `K[y]` meets `gl_5(O)` in
`O[y]` (for each `j <= 4` some entry of `y^j` is a unit). Call `pr` the projection onto `K[y]`. Grading: `X = Sum_a xi_a E_{a,a+m}` has
degree `m` (`1 <= m <= 5`, degree 5 = degree 0), like `At^m = Sum_a W_a^(m) E_{a,a+m}`, `W_a^(m) = bt_a bt_{a+1} .. bt_{a+m-1}`, `W_a^(5) = B`.
Then `pr(X) = (1/5) [Sum_a xi_a / W_a^(m)] At^m = (1/5) [Sum_a xi_a B/W_a^(m)] y^(5-m)`. Since `B/W_a^(m)` is a product
of `5 - m` weights, at most one of them `bt_1`: **for holomorphic `X` the coefficient of `y^j` in `pr(X)` has at most a simple pole
(`j = 1..4`) and is holomorphic for `j = 0`.**
(iii) *Lemma 1.* If `X` in `gl_5(O)` and `pr(X) = 0`, there is `Pi` in `gl_5(O)` with `[At, Pi] = X`, unique up to `O[y]`. Indeed for
`Pi = Sum_a pi_a E_{a,a+d}` of degree `d` (0 <= d <= 4) with `pi_a = rho_a W_a^(d)` one finds `[At, Pi] = Sum_a (rho_{a+1} - rho_a)
W_a^(d+1) E_{a,a+d+1}`. So `rho_{a+1} - rho_a = sg_a := xi_a/W_a^(d+1)`, which is holomorphic and vanishes at `u = 0` whenever the window
`{a, .., a+d}` contains the pole index 1; `Sum sg_a = 0` is `pr(X) = 0`, and the `rho_a` are determined up to a common additive function.
`pi_a` is holomorphic iff `rho_a` is and `rho_a(0) = 0` for the `d` values of `a` whose window `{a, .., a+d-1}` contains the pole index,
`a = 2-d, .., 1`; for `a = 2-d, .., 0` the longer window contains it too, so `sg_a(0) = 0` and these `d` conditions are the single condition
`rho_1(0) = 0`, which the additive function fixes (no condition for `d = 0`). Two holomorphic solutions differ by an element of `O[y]`.
(iv) *Lemma 2.* For `Z` in `gl_5(O)` the equation `pr(c' + Dl c) = -pr(Z)` has a solution `c = Sum_j c_j y^j` in `O[y]`, unique up to
`c_0 -> c_0 + const x exp(-Int tr Dl/5)`. Indeed `pr((y^j)') = -(j/5) (log B)' y^j = (j/(5u) + hol) y^j` (every weight occurs in `j` of the
five windows) and `pr(Dl y^j) = (tr Dl/5) y^j`, so the equation is `u c_j' + (j/5 + u h_j) c_j = -u phi_j` with `h_j`, `u phi_j` in `O` by
(ii): a regular singular equation with exponent `-j/5`. For `j = 1..4` the pivots `n + j/5` never vanish: a unique power-series
solution. For `j = 0` the equation is regular: a one-parameter family.
(v) *Induction.* `Ph_0`: `Pi_0 = 0`, and Lemma 2 with `Z = 0` gives `Ph_0 = rho_0(u) 1`, `rho_0 = const x exp(-Int tr Dl/5)`. Given `Ph_k =
Pi_k + c_k` with `Pi_k` in `gl_5(O)` a solution of the order-`k` equation (Lemma 1), choose `c_k` in `O[y]` by Lemma 2 with `Z = Pi_k' +
Dl Pi_k`; then `X_{k+1} = Ph_k' + Dl Ph_k` is holomorphic with `pr(X_{k+1}) = 0` and Lemma 1 gives `Pi_{k+1}`. The only freedom is one
constant per order, which is `P -> P c(z)`. QED

**Remarks.** (a) The proof uses only the shape of the two systems; the model is *versal* for cyclic systems with one simple pole.
(b) Geometrically (ii) says that the lattice `E = + O Y_a` is the structure sheaf of the spectral curve `y^5 = 1/B(u)` (smooth, totally
ramified over `u = 0`, generator `Y_1`, `Y_a ~ y^(c_a + 1)`), and the exponents `j/5` are the ramification weights. (c) Necessity of
`prod b = prod bt`: `P_0` must conjugate `A` into `At`. For a change of the independent variable `v(u)` this is the equality of
the canonical one-forms and fixes `v = u tau^(-5/4) = u + (5/9) u^2 + (220/567) u^3 + ...` uniquely (given `v(0) = 0`). (d) *Computation*
(ur_gauge.py, exact over Q, inputs: `varpi_0^cf`, `t_c`, `tau` only): `p_{k,b}` for `k <= 16` to `O(u^18)` at least; every one is a power
series; the residues of the ODEs are `(0, 1/5, -3/5, -2/5, -1/5)` for `k = (0, 1, 2, 3, 4) mod 5`, congruent to `k/5` mod 1 (non-integral
except for `k = 0 mod 5`), as in Lemma 2. PASS (ur_gauge_run.log).

## 4. From the gauge to the factorisation; the natural variable; the normalisation

Apply Theorem U to `A` (quintic, any conifold frame) and `At = v'(u) A^mod(v(u))` (the model pulled back by `v = u tau^(-5/4)`). The model
rows at `v(u)` are `Q^m_b = sigma^(c_b-1) R^_b(z/sigma^4)` with `sigma = s tau^(-1/4)`, so `z/sigma^4 = tau zhat`: **the model R-matrix is a
function of `z` over its canonical coordinate, and the pull-back makes that the canonical coordinate difference of the quintic**
(`tau zhat = -(5/4) zeta^i z/(u_i - u_i^c)`). `P Q^m` is a formal solution of the quintic system with the exponentials `e^{(u_i - u_i^c)/z}`,
hence equals the quintic WKB columns up to a constant diagonal series: `Q D(z) = P Q^m(v(u))`. Both systems are covariant under
`z -> zeta z` combined with `diag(zeta^-a)`, so `P` can be taken graded and then `D_i(z) = D_0(zeta^-i z)`. In hatted variables

    Qhat(u, zhat) D = Ghat(u, zhat) R^(tau zhat),    Ghat_ab = tau^((1-c_b)/4) Sum_k zhat^k u^(mu(b,k)) p_{k,b}(u) delta_{a,b+k},
    mu(b, k) = (4k + c_b - c_{b+k})/5  (an integer >= 0),

and the right-hand side is a power series in `u` whose coefficients are power series in `zhat`.

**(H) Purity.** The hatted rows of the formal quintic (normalised by the QRR constants `h_{wt_i}(z)` at `q = 0`) are power series in `u`.
Equivalently `D` is a series in `z^5` (then scalar by covariance, and absorbed into `P`). VERIFIED mod `z^7`, `u^6` (the assertion in
`fq3_msum.HatRowsU`; also `Symbol.val_ok`). Not proved to all orders: it is a statement about the constants of integration at `q = 0`
(a connection formula between the large-radius and the conifold normalisations of the WKB series); the oscillatory-integral picture
of s.8 makes it plausible. Without (H) everything below holds with `M_0(tau zhat)` replaced by `M_0(tau zhat) D(z)^-1`.

**Corollary F (PROVED given (H)).** `R~(zhat; u) = G(u, zhat) M_0(tau(u) zhat)`, `G = Psi(u)^-1 Ghat(u, zhat) Psi(0)`, where `Ghat` is the
hatted form of a regular gauge; `[u^n] G` is a polynomial in `zhat` (s.5). The same holds in every conifold frame `(alpha, beta)` with the
same `M_0(tau zhat)` and a frame-dependent `G`.

**Verification (ur_gauge_run.py, log ur_gauge_run.log).** `P` was computed from `(varpi_0^cf, t_c, tau)` over Q (`k <= 16`) and compared with
`Ghat^proj := Qhat(u, zhat) R^(tau zhat)^-1` built from the project's rows and symbol:

| check | result (both primes) |
|---|---|
| Z/5 grading of `Ghat^proj` | PASS |
| `zhat^0..zhat^4`, all 25 entries, `u^0..u^5`: `Ghat^proj` = gauge, no free constant | PASS |
| `zhat^5`, `zhat^6` through `u^3` | PASS |
| `zhat^5`, `zhat^6` through `u^5` with `P -> P (1 + c_5 z^5)`, `c_5` fitted on the entry (0, 0) at `u^4 zhat^5` | PASS, `c_5 = 43397/807975000` |
| `c_5 = -p_{5,1}(0) = -p_{5,2}(0)` in the normalisation `p_{5,0}(0) = 0` | PASS |

The last line says that the gauge of the ACTUAL R-matrix satisfies `P_11(0, z) = P_22(0, z) = 1 + O(z^10)` (recall `P_11(0, z) = P_22(0, z)`
always, (R1) below), while its other diagonal entries at `z^5`, `u = 0` are `(43397/807975000, ., ., 283/2835000, -12659/201993750)` for
`b = 0, 3, 4`.

**Canonical normalisation (CONJECTURAL; verified at the only order available, `z^5`, one rational number, both primes).** The regular
gauge normalised by `P_11(0, z) = 1` is the gauge of the actual R-matrix with `D = 1`. Heuristic: `P_22(0, z)` compares the two `H^2`-rows
at the conifold point, i.e. the coefficients of `log u` (the pairing with the vanishing cycle) of the WKB columns, and in the model the
row `-x^3 dx` is exact at `v = 0`.

**Consequence (practical).** Near the conifold point the R-matrix of the formal quintic is determined, to any order in `z`, by the
periods (through `P`) and the universal series `F_b` of s.2. For genus 4 one needs `z^9`: `ur_gauge.gauge(9, ..)` and
`ur_model_check.model_rows(9)` give it without any new WKB computation (given (H) and the constant `c_5`).

## 5. Structure of the gauge: grading, residue relations, degrees

**Grading (PROVED).** `c_a = -a + 5 e_a`, `e = (0, 0, 1, 1, 1)`, gives `mu(b, 5q + r) = 4q + mu(b, r)` with

| r | a = b + r for b = 0..4 | mu(b, r), b = 0..4 | after (R1), (R2) | min |
|---|---|---|---|---|
| 0 | 0 1 2 3 4 | 0 0 0 0 0 | 0 0 0 0 0 | 0 |
| 1 | 1 2 3 4 0 | 1 0 1 1 1 | 1 2 1 1 1 | 1 |
| 2 | 2 3 4 0 1 | 1 1 2 2 2 | 2 2 2 2 2 | 2 |
| 3 | 3 4 0 1 2 | 2 2 3 3 2 | 2 3 3 3 3 | 2 |
| 4 | 4 0 1 2 3 | 3 3 4 3 3 | 3 4 4 4 3 | 3 |

`u^(mu(b,k))` is the least power of `u` at which the flat-frame term `z^k E_{b+k,b}` can appear in hatted variables. The only entry with
`mu = 0` at `k >= 1` is `(a, b) = (2, 1)`, `k = 1`: the flat-frame term `z H (x) H` (the entry `(a, b)` is the bivector `H^(3-a) (x) H^b`).

**Residue relations (PROVED).** Write `A = r E_12/u + A_0 + u A_1 + ..`, `At = rt E_12/u + At_0 + ..` (`r = rt = -1` here), `P = P(0) + u P^(1) + u^2
P^(2) + ..` (coefficients series in `z`). The gauge equation `z P' = A P - P At` gives
- order `u^-1`: `r E_12 P(0) = rt P(0) E_12`, i.e. **(R1)** `P_2b(0, z) = 0` for `b != 2`, `P_a1(0, z) = 0` for `a != 1`, `r P_22(0, z) = rt P_11(0, z)`:
  at the conifold point the `H^2`-row of the quintic is proportional to the `H^2`-row of the model, and the `H`-row of the model enters only
  the `H`-row of the quintic;
- order `u^0`, entry `(2, 1)`: the residue terms have no `(2, 1)` entry and `(A_0 P(0) - P(0) At_0)_21 = b_2(0) P_31(0) - bt_0(0) P_20(0) = 0` by (R1),
  so **(R2)** `P_21 = O(u^2)` at every order in `z`;
- the same orders on the entries `(2, 0)`, `(3, 1)` and order `u^1` on `(2, 1)`: with `pi(z) := P_30(0, z)` and `b_2(0) = bt_0(0) = 1`,
  **(R3)** `P_20 = u pi/z + O(u^2)`, `P_31 = -u pi/z + O(u^2)`, `P_21 = -u^2 pi/z^2 + O(u^3)`: the block `rows (2, 3) x columns (0, 1)` is
  `pi(z) (u/z, 1)^T (1, -u/z)` to leading order in each entry, i.e. the bivector `pi(z) (1 + (u/z) H) (x) (1 - (u/z) H)`.

VERIFIED exactly over Q for `k <= 14` (ur_chain_check.py): (R1), (R2), and (R3) with `pi = (127/3150) z^3 - (77141089/362679778125000) z^8 - ..`.
In hatted variables `1 +- (u/z) H = s^-1 (1hat +- Hhat/zhat)`, so the `z^3` term is `-(127/3150) u^2 zhat (Hhat + zhat)(x)(Hhat - zhat)` up to
the ordering convention of the two factors: this is the observed block of `[u^2] G` (`Xhat(3,3) = Xhat(3,4) = -Xhat(4,3) = -127/15750` at
`zhat`, `zhat^2`, `zhat^2`, and `Xhat(4,4) = +127/15750` at `zhat^3`; `127/15750 = (1/5)(127/3150)`), i.e. `alpha zhat (H - zhat 1)(x)(H + zhat 1)`
in the Fourier convention of genus3_natural.py. `z 1 +- u H` is the Euler-dilaton vector of the homogeneous model (s.6.4).

**Answer to (c).** The absence of `H (x) H` from `G_1` (and from `G_0`) is FORCED: it is (R2). The component first appears in `G_2`, tied by
(R3) to the `1 (x) 1` term of `zhat^3` and the `H (x) 1 - 1 (x) H` terms of `zhat^2`.

**Theorem D (degrees; PROVED as a bound given Theorem U and (H)).** `[u^n] G` is a polynomial in `zhat` of degree at most
`floor((5n + 2)/4) = n + floor((n + 2)/4)`: for `n = 4q + j`, degree `5q + (0, 1, 3, 4)_j`.
*Proof.* The entry `(b + k, b)` of `[zhat^k] Ghat` is `u^(mu(b,k) + eps)` times a power series, `eps = ord_u p_{k,b}`; by (R1), (R2) `eps >= 2` for
`(a, b) = (2, 1)` and `eps >= 1` if `a = 2` or `b = 1` (and `5` does not divide `k`). The fourth column of the table gives the first order `n(k) =
4q + (0, 1, 2, 2, 3)_r` at which `zhat^k` can occur (`k = 5q + r >= 1`), and `Psi(u)` is a power series in `u`. QED
*Sharpness* (VERIFIED by the recursion for `k <= 16`, i.e. `n <= 12`; for `n <= 4` also by the production data): no further vanishing occurs;
degrees `0, 1, 3, 4, 5, 6, 8, 9, 10, 11, 13, 14, 15` for `n = 0..12`. Top-degree components: `n = 1 mod 4`: the two Pi-type propagators `1 (x)
H^2 + H^2 (x) 1`, `H^3 (x) H^4 + H^4 (x) H^3` at `zhat^(5q+1)`; `n = 2 mod 4`: `1 (x) 1` at `zhat^(5q+3)`; `n = 3 mod 4`: `1 (x) H^4`, `H^4 (x) 1` at
`zhat^(5q+4)`; `n = 0 mod 4`: the diagonal (`H^(3-a) (x) H^a`) terms at `zhat^(5q)`.

**The order-u propagator (PROVED formula, VERIFIED).** With `bt_a/b_a = 1 + eps_a u + ..`: `eps = (e, -2e - 2f, e, f, f)`, `e = 2 v_2 - 2 gamma_t =
-13/45`, `f = 2 v_2 - 1 = 1/9` (`v = u + v_2 u^2`, `v_2 = 5/9`; `t_c = u + gamma_t u^2`, `gamma_t = 7/10`), and

    p_{1,b}(0) = (-e - f, 0, -e - f, e + 2f, e + 2f) = (8/45, 0, 8/45, -1/15, -1/15),      G_1 = zhat X_1,

`X_1` having the Fourier components `Xhat(2,4) = Xhat(4,2) = -(1/5)(8/45) = -8/225` and `Xhat(0,1) = Xhat(1,0) = (1/5)(-1/15) = -1/75`. So `X_1` is
explicit in the first-order data of the periods. In a general conifold frame (`varpi_0 = 1 + w_1 u + ..`, `t_c = u + gamma_t u^2 + ..`,
`gamma_t + w_1 = gamma_g = 7/10` the `u^2` coefficient of the vanishing period) `e = 2 v_2 - 2 gamma_t`, `f = 2 v_2 - 1 - w_1`, hence

    Xhat(0,1) = (e + 2f)/5 = -(2 gamma_g - 4/3)/5 = -1/75   (frame INDEPENDENT),      Xhat(2,4) = (e + f)/5 = -8/225 + w_1/5   (frame dependent).

**Check in a second frame (VERIFIED, 2^89-1):** for `(alpha, beta) = (1/3, -2/7)` (`w_1 = 1/3`) the formula gives `Xhat(2,4) = 7/225`, `Xhat(0,1) =
-1/75`, which are the values printed by the production run in that frame (logs/genus3_natural_89_nat_ab1o3_m2o7.log). So the choice of
conifold frame moves `X_1` along the string direction of the plane Pi only (the direction that is an exact gauge direction by the
string equation), and the other coefficient is the invariant `2 gamma_g - 4/3` that converts the model coordinate into `t_c` (s.6.2).

## 6. Consequences for the graph sum

### 6.1 Counting

Factorise the hatted R-matrix as `R = R' M` with `M = M_nat := M_0(tau(u) zhat)` and `R' = Psi(u) G(u, zhat)`. The vertex CohFT is `V_tau =
M_0(tau zhat).omega^`: the model theory in its natural variable (degree-`d` classes of the frozen theory `M_0.omega^` times `tau^d`). Edges
`(1 - G(-psi)^T G(-psi'))/(psi + psi')` and leaves `z(1 - G^-1 1)` are `O(u)`. A genus-g decoration is `u^(2-2g + Sum_e n_e + Sum_l n_l)`
times a power series, `n_e, n_l >= 1` the u-orders taken from the edge and leaf factors. Hence: a graph with `E` edges is `O(u^(E + 2 - 2g))`;
graphs with `E >= 2g - 2` are regular for free; for the others only the finitely many terms with `Sum n <= 2g - 3` matter, and they
involve `G_1, .., G_{2g-3}` only, i.e. by Theorem D zhat-orders `<= floor((10g - 13)/4)` of the gauge (4 at genus 3, 6 at genus 4).

### 6.2 What the genus-3 run says (VERIFIED: frame (0, 0) at both primes, frame (1/3, -2/7) at 2^89-1; logs/genus3_natural_*_nat*.log)

With `M = M_nat` exactly (no dressing; `M` taken mod `u^4`, `zhat^7`): **all 41 edge graphs of genus 3 are `O(u^4)`, i.e. regular**; the total
equals the point-frame total; the vertex term is `u^-4 (1/25200 - u/9000 + 407 u^2/3780000 - 151 u^3/3780000 + ..) = (1/25200) t_c^-4 + 0
t_c^-3 + 0 t_c^-2 + 0 t_c^-1 + ..`. **Genus 2 with the same `M`**: the six edge graphs are `O(u^2)` (`[u^2]`: 479/756000, 1/2250, -127/2016000,
0, 0, 0), the vertex term is `u^-2 (1/1200 - 7u/6000 + ..)`. The same holds in the frame `(1/3, -2/7)` (2^89-1: 41 edge graphs `O(u^4)`,
vertex term `u^-4 (1/25200 - 11u/189000 - ..) = (1/25200) t_c^-4 + 0 + 0 + 0`, genus 2 regular), with the same `M_0(tau zhat)` and the
frame-dependent `G`, as Corollary F predicts. So **nothing has to be added to the vertex factor** through genus 3: the
`H (x) H` part of `G` may stay in `R'`. Of the conditions of OU_EXTENSION_genus3 s.6d: the `N_1`-conditions (4 linear, 5 quadrics, 9 cubics)
hold because `Delta = N_1^nat - N_1^taut = -M_0^-1 (zhat X_1) M_0` lies in the kinematic plane Pi (`X_1` has no `H (x) H` by (R2)); here
`N_1^nat = tau_1 M_0^-1 zhat dM_0/dzhat`, `tau_1 = -4/9`, is the generator of the rescaling of `zhat`. The remaining nine numbers (`[u^2]`, `[u^3]` of
the two one-edge graphs, `[u^3]` of the five two-edge graphs), which in s.6d were solved for with free parameters `N_2`, `N_3`, vanish for
the canonical `N_2^nat`, `N_3^nat` with no parameter to adjust: this is the new evidence. So the natural form is ONE point of the
solution set of s.6d (first order: the point `D_1 = -X_1` of Pi; not the tautological point, not pure-high, not one of the 17 rays),
singled out by an independent principle.

*The vertex term is not the bare model potential.* By homogeneity the leafless vertex is `c_g (u tau^(-3/2)/xhat)^(2-2g)` (`c_2 = 1/1200`,
`c_3 = 1/25200`), whose `u^1` coefficient is `-(3g - 3)(4/9) c_g`, while the gap needs `-(2g - 2)(7/10) c_g`; the difference `-(g - 1) c_g/15`
(`-1/18000` at genus 2, `-1/189000` at genus 3, from the stored values) is supplied by the dilaton leaves of `R'`, and `1/15 = -p_{1,4}(0)` is
the coefficient of the unit row of the gauge at order `z` (OBSERVED from the run values; the single-leaf formula was not derived). In
a general conifold frame `-p_{1,4}(0) = 2 gamma_t + 2 w_1 - 4/3` (`w_1 = [u^1] varpi_0`), which is exactly what the gap requires of such a
leaf (`-(2g-2)(gamma_t + w_1) + (3g-3)(4/9) = -(g-1)(2 gamma_t + 2 w_1 - 4/3)`), so this bookkeeping is frame independent.

### 6.3 Mechanism: what is understood and what is not

(i) *First order* (kinematic, OU_EXTENSION_genus3 s.6d): `u^1` edges are psi-free propagators in `Sym^2 span(1, H^2, H^3, H^4)`, and the
frozen vertex theory has genus `>= 1` primary correlators only in the direction `H`.
(ii) *The Euler-dilaton part* (**Proposition E**; PROVED from the string and dilaton equations and the primary sector `<H^n>_h = d_t^n
(c_h t^(2-2h))` at `t = -1`, `h >= 2`, which is VERIFIED for `(h, n) = (2, 0..3), (3, 0)` only; VERIFIED directly, both primes,
ur_euler_linear.py). For the left dressing `exp(eps n)`, `n(zhat) = zhat (H - zhat)(x)(H + zhat)`, the edge is, up to an overall sign,
`eps [H(x)H + (psi - psi')(H(x)1 - 1(x)H) - (psi^2 - psi psi' + psi'^2) 1(x)1] + O(eps^2)` and the first-order response of a one-edge graph is
(up to that sign)
- loop on a vertex of genus `h`: `(1/2)(<H,H>_h - (4h - 4)<H>_h + (2h - 2)(2h - 3) F_h) = 0` for `h >= 2`;
- bridge `h_1 - h_2`: `<e>_{h_1} <e>_{h_2}` with `<e>_h = <H>_h - <tau_1(1)>_h = 0` for `h >= 2` and `<e>_1 = 1/12 - 5/24 = -1/8`.
Measured (eps^1 coefficients, frozen theory): genus 3: loop on `g = 2`: 0, bridge `g1 - g2`: 0; genus 2: loop on `g = 1`: `1/16 = -(1/2)(-1/8)`,
bridge `g1 - g1`: `-1/128 = -(1/2)(1/8)^2`, vertex graph `-7/128` (the three add up to zero). Controls: the opposite relative sign
`zhat (H + zhat)(x)(H - zhat)` gives `-1/150`, `-7/7200` at genus 3, the bare `zhat H(x)H` gives `-1/400 = -(1/2)<H,H>_2`, `-1/7200`. So at every
genus `g >= 3` the `alpha u^2` term of (R3) drops out of the `[u^2]` coefficient of the one-edge graphs.
(iii) **Negative result.** `n` is NOT an exact symmetry of the frozen vertex theory: for `M = exp(eps n) M_0` all genus-2 edge graphs and the
genus-3 graphs 2, 5 are non-zero at finite `eps` (ur_euler_check_89.log; e.g. `eps^2` coefficients 7/128, 7/512 at genus 3). The
obstruction is the genus-1 anomaly `<e>_1 = -1/8` (and, at second order, pairs of Euler-dilaton edges). The apparent invisibility of
`zhat^3 (H - zhat)(x)(H + zhat)` in the same log is trivial (not Z/5-covariant: the colour sums vanish).
(iv) *Not explained:* the vanishing of the remaining `[u^2]`, `[u^3]` contributions (the non-`H` parts of `G_2`, `G_3`, products `G_1 G_1`,
leaves, the `tau`-dependence of the vertices). There is no all-order mechanism yet.

### 6.4 The Euler-dilaton vector

In the flat frame `e(z) = z 1 + u H`. For the homogeneous model `F_h = c_h v^(2-2h)`, `H = d_v`: `<tau_1(1) X> + v <H X> = (2h - 2 + n + v d_v) <X>`,
which vanishes on the primary `H`-sector. (R3) says that the part of the regular gauge that is not suppressed by a power of `u` in the
`H`-directions is `-pi(z) z^-2 e(z) (x) e(-z)`: the residue structure of the connection at the conifold point produces exactly the
combination that the homogeneity of the model theory annihilates (up to the genus-1 anomaly). This is the structural reason to expect
regularity beyond the orders tested, and the place where it could fail.

### 6.5 Conjecture and programme

**Conjecture N (natural factorisation).** For every `g >= 2`, in the hatted graph sum of the formal quintic with `M = M_0(tau(u) zhat)`,
`tau = 1/2F1(1, 4/5; 9/5; u)`, `M_0` the R-matrix of the LG model `-(x^4/4 + v/x)`, `dx` (s.2), and `R' = Psi(u) G(u, zhat)`, `G` the hatted regular
gauge of s.3-4 (canonical normalisation), every stable graph with at least one edge is `O(u^(2g-2))`, i.e. regular at the conifold
point; equivalently the vertex term (with its leaves) is `c_g t_c^(2-2g)` + regular, `c_g = (-1)^(g-1) B_2g/(2g(2g-2)) 5^(1-g)`.
Flat-frame form: `Omega^quintic = P.(Omega^LG at v(u))` with `P` regular at `u = 0`, and all graphs with edges regular.
Status: VERIFIED at `g = 2` and `g = 3` (frame (0, 0): two primes; frame (1/3, -2/7): 2^89-1). CONJECTURAL for `g >= 4`.

What would have to be proved:
1. (H) to all orders, and the canonical normalisation (s.4). [verified mod `z^7`]
2. *The gap of the model*: the genus-g potential of the LG model on the symmetric locus is `c_g v^(2-2g)` with the Bernoulli constants
   (the power is homogeneity; the constant is verified at `g = 2, 3`), `F_1 = -(1/12) log v`. A route: the total descendant potential of
   the Hurwitz Frobenius manifold of functions with one simple pole is a tau function of a rational reduction of KP (arXiv:2211.12259,
   not read in detail); its restriction to the residue direction should be the Gaussian/c = 1 string free energy.
3. *Primary-sector lemma*: at genus `>= 1` the primary correlators of the model at the symmetric point vanish unless all insertions are
   `H` (i.e. `F_g` of the Hurwitz Frobenius manifold depends on the residue `u` only, to the relevant order). [verified for `(1, <= 4)`,
   `(2, <= 3)` in s.6d]
4. *Edge/leaf structure lemma* (the analogue of Lemmas A and C of the conifold programme): every decoration with `Sum n_e + Sum n_l <=
   2g - 3` cancels within its graph. Available inputs: Z/5 selection, string, Euler-dilaton (Proposition E), item 3, Theorem D (which
   components of `G_n` carry the top psi-degrees: Pi-type, `1 (x) 1`, `1 (x) H^4`, diagonal), (R1)-(R3). Missing: any control of the
   descendant correlators of the model beyond these Ward identities, and of the genus-1 anomaly in graphs with more edges.

*Genus 4 as a test.* It needs `G_1..G_5` (degrees 1, 3, 4, 5, 6; all computed) and `M_0` mod `zhat^10` (s.2), 379 graphs; new features: the
diagonal `u^4 zhat^5` terms of `G_4` (values in s.4; `0` on the `H`, `H^2` entries), the second Euler-dilaton term only at `u^6` (`pi_8`), and
graphs in which an Euler-dilaton edge sits on a genus-1 vertex together with other half-edges (the anomaly). If Conjecture N fails
there, the natural repair is a left dressing `M = A M_0(tau zhat)` by the offending low-order part of `G` (the `--variant` mechanism of
genus3_natural.py); the data through genus 3 give no indication that this is needed.

## 7. The analogy with the Hodge trick

| | resolved conifold x C^2, local P^2 (Hodge vertex) | formal quintic at the conifold point |
|---|---|---|
| model integral | `Int_0^inf exp(-(x - w log x)/z) dx/x = z^(w/z) Gamma(w/z)`, one for each normal weight `w` | `Int exp(-(x^4/4 + v/x)/z) x^m dx`, `m = 0, -1, 3, 2, 1` |
| LG model, primitive form | `lambda = x - w log x`, `dx/x`; one critical point `x = w` | `lambda = x^4/4 + v/x`, `dx`; five critical points `x^5 = v` coalescing on the pole |
| Frobenius manifold | rank 1 per colour (trivial); CohFT = Hodge classes `Lambda^vee(w)` of the normal bundle (GW of C^3) | rank 5, colour-mixing: Hurwitz space `M_{0;3,0}` on its Z/5-symmetric locus |
| Mellin-Barnes kernel | `Gamma(s)` | `Gamma(s) Gamma((s+1)/4)` = five Gamma factors with exponents `0, 1/4, 1/4, 1/2, 3/4` |
| model R-matrix | Stirling series `exp(-+ Sum_k B_2k/(2k(2k-1)) (z/w)^(2k-1))`: scalar, exp(odd), Bernoulli numbers | `M_0(z/sigma^4)`: the series `F_b` of s.2 (recursion with pivot `-4n`); symplectic, Z/5-covariant, not exp(odd), not diagonal |
| natural variable | `z/w` (homogeneity `x, w, z -> k x, k w, k z`); `w` is constant, so `M` is constant in `q` | `z/v^(4/5)` = `z` over the canonical coordinate difference; it blows up at the conifold, so `M` must depend on `u`: `M_0(tau(u) zhat)` |
| change of variable | none | `v = u tau^(-5/4)`, `tau = 1/2F1(1, 4/5; 9/5; u)`, forced by the canonical one-forms |
| remaining factor `R'` | `R_Z h_w^-1` | `Psi(u) G(u, zhat)`: hatted form of the gauge `P(u, z)`, regular at `u = 0` in the flat frame, computable from the periods; `[u^n] G` polynomial of degree `floor((5n+2)/4)` |
| vertex term | Hodge integrals; carries the whole gap | model potential `c_g v^(2-2g)` dressed by the dilaton leaves of `R'`; carries the whole gap (g = 2, 3) |
| edge graphs | regular (genus 2-4 for local P^2; conifold: programme of Lemmas A, C) | regular at genus 2, 3 (41 graphs at genus 3); Conjecture N |
| why a scalar factor fails / works | weights are constants; colours do not interact at the singular point | the five colours coalesce (`u_i - u_j -> 0`); only a colour-mixing model can absorb the singular part (CONIFOLD_SYMBOL s.1) |

In both cases the trick is the same uniform-asymptotics statement: *R = (regular gauge) x (R-matrix of the universal local model in its
natural variable)*, the model being the simplest oscillatory integral with the same singular structure: one critical point next to a
logarithmic singularity (Gamma), respectively five critical points coalescing onto a simple pole (the 2-saddle analogue `x + v/x` is
the Bessel-type uniform expansion known from Laguerre polynomials, and is the GUE/NLS Frobenius manifold `M_{0;0,0}`). What the formal
quintic lacks so far is the analogue of the closed Hodge-integral evaluation of the vertices: the all-genus theory of the model (s.6.5,
items 2-3).

## 8. Oscillatory-integral picture (SKETCH, with two exact statements)

*Exact statement 1 (PROVED, elementary).* Let `Phi(m) = Sum_j zeta^j log(m - zeta^j) = -5 Int_0^x t^3/(1 - t^5) dt`, `x = 1/m`. At the points
`u m^5 = 1` (`m = zeta^i L`): `Phi = -(5/4) x^4 2F1(1, 4/5; 9/5; u) = u_i - u_i^c`, and `d Phi/du = -m/(1-u) = du_i/du`. So `tau` is literally
`x_c^4` over (4/5 of) the critical value.
*Exact statement 2 (PROVED).* The one-variable function `W(x; u) = (5/x)[1 - ((1 - x^5)/(1 - u))^(1/5)] - 5 Int_0^x t^3/(1 - t^5) dt` has
`dW/dx = 0` exactly at `x^5 = u`, critical values `u_i - u_i^c`, a simple pole at `x = 0` with residue `5(1 - (1-u)^(-1/5)) = -u + ..`, and
`W = -(u/x + x^4/4) + ` higher order.
*Interpretation (SKETCH; the form of the mirror is recalled, not checked against a reference).* The I-function of the formal quintic is
the stationary phase of `Int exp(W_5/z)`, `W_5 = Sum_j (x_j + zeta^j log x_j) - y` on `prod x_j = q y^5` (equivariant mirror of P^4 with the
Laplace variable `y` of the twist); its critical points are `x_j = m - zeta^j`, `y = 5m`, `u m^5 = 1`, which run to infinity at the
conifold point along the curve on which `W_5` restricts to `W(1/m; u)`. Integrating out the four transverse (Morse) directions leaves a
one-dimensional oscillatory integral with phase equivalent to `W` and an amplitude regular at `x = 0`. `W` is a deformation of
`-(x^4/4 + v/x)` inside the functions with one simple pole; the versal family is the Hurwitz family `x^4 + a x^2 + b x + c + u/(x - d)`,
and Z/5-covariance keeps the quintic on `a = b = d = 0`: a change of variable `x -> X(x, u)` brings `W` to `-(X^4/4 + v(u)/X)` (the
critical values force `v`), and reducing the amplitude modulo exact forms in the Brieskorn lattice of the model gives the regular
gauge `P`. This is the Chester-Friedman-Ursell mechanism; Theorem U is its formal, ODE-level version and does not depend on this sketch.
It also suggests why (H) holds: the QRR constants are the Stirling series of the Gamma factors of the same integral.

## 9. Files and checks (all in natural_form/theory/, single process, exact)

| file | content | result |
|---|---|---|
| `ur_model_check.py`, `.log` | (A) cyclic system incl. the closing relation on the s-series rows; (B) model series `F_b` over Q vs the symbol, all rows/colours, `zhat^0..6`; stationary phase through `t^10` | ALL PASS, both primes |
| `ur_gauge.py`, `ur_gauge_run.py`, `.log` | Theorem U as a recursion over Q (`k <= 16`); regularity; first u-orders; degrees for `n <= 12`; comparison with the project's `G` (`zhat^0..6`, `u^0..5`); `c_5` | ALL PASS, both primes |
| `ur_chain_check.py`, `.log` | (R1), (R2), (R3) for `k <= 14`, the formula for `p_{1,b}(0)` | ALL PASS (exact, Q) |
| `ur_c5.py` | `c_5` at `2^61-1` against the rational value | True |
| `ur_euler_linear.py`, `_89.log`, `_61.log` | Proposition E: first-order response of the Euler-dilaton element, two controls | ALL PASS, both primes |
| `ur_euler_check.py`, `_89.log` | exact invisibility of `exp(eps n)`: REFUTED (two intended FAIL lines; s.6.3 (iii)) | negative result, `2^89-1` only |

Read-only inputs: the project modules (`fq_*`, `fq3_*`), `logs/genus3_natural_{89,61}_nat.log` and
`logs/genus3_natural_89_nat_ab1o3_m2o7.log` (the production runs, which finished while this note was written). No production file was
written or modified.

## 10. Caveats

- AI-checked only. Theorem U and the residue relations are proved here in full but have not been refereed; the computations confirm
  them for `k <= 16` (regularity) and against the project's independent R-matrix for `zhat^k`, `k <= 6`, `u^n`, `n <= 5`.
- (H) (purity of the hatted rows, `D = 1`) and the canonical normalisation `P_11(0, z) = 1` are verified only mod `z^7` (the latter is a
  single number). Proposition S and Corollary F are conditional on (H).
- The sharpness of Theorem D (no accidental vanishing) is a finite computation (`n <= 12`).
- Proposition E uses the primary sector of the frozen vertex theory, known only for the `(g, n)` listed; its direct verification covers
  the one-edge graphs of genus 2 and 3 only. The single-leaf bookkeeping of s.6.2 is an observation.
- Conjecture N rests on genus 2 and 3 (frame (0, 0) at two primes, a second frame at one prime). No mechanism is known for the cancellations beyond first order,
  and the genus-1 anomaly of the Euler-dilaton vector shows that the model is not simply invariant under the low-order part of the
  gauge. Genus 4 is the first test of the diagonal `u^4 zhat^5` terms.
- The identification of the model with a point of a Hurwitz Frobenius manifold is at the level of the Gauss-Manin system on the
  one-parameter locus; nothing was checked about the full five-dimensional Frobenius structure, about the gap constants of the model
  beyond `g = 2, 3` (from the frozen theory), or about the literature statement quoted in s.6.5 (arXiv:2211.12259 was not read).
- The mirror `W_5` of s.8 is recalled, not derived; only the two exact statements there are claimed.
- Sign and ordering conventions of bivectors (first/second tensor factor, overall sign of an edge) are those of genus3_natural.py for
  quoted data and are otherwise fixed only up to the stated ambiguities; the flat-frame matrix identities (R1)-(R3) are unambiguous.

## Sources (web)

- arXiv:2211.12259 (Hurwitz Dubrovin-Frobenius structure for one simple and one higher-order pole; rational reductions of KP):
  https://www.arxiv.org/abs/2211.12259
- Chester, Friedman, Ursell, "An extension of the method of steepest descents", Proc. Camb. Phil. Soc. 53 (1957) 599-611; a modern
  instance with two saddle points coalescing onto a simple pole (Laguerre polynomials): https://worldscientific.com/doi/10.1142/9789814656054_0022
- Uniform Airy-type expansions of integrals (CWI report): https://ir.cwi.nl/pub/1473/1473D.pdf
