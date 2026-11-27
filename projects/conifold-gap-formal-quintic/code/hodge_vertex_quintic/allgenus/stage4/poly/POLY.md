# Stage 4 / poly: D_0 = 1 from the polynomial structure of the hatted WKB functions (2026-11-14)

Folder `allgenus/stage4/poly/` (new files only; no project module imported or modified). AI-derived and AI-checked only; no human has
read it. Labels: PROVED (argument written here, or quoted with arXiv number and statement), VERIFIED (exact computation, range stated),
CONJECTURAL, SKETCH, REFUTED. Conventions: `natural_form/theory/UNIFORM_REDUCTION.md` (UR), `allgenus/gauge/GAUGE_IDENTITY.md` (GI);
colour 0 throughout (colour i: `Phi^(i)(w,u) = Phi^(0)(w/zeta^i, u)`, PROVED: HPF depends on `w/zeta^i` only and `h_(w_i)(z) = h_(w_0)(z/zeta^i)`).
`w := zhat`, `z = u^(4/5) w`, `c = (0,-1,3,2,1)`, `Lam := -log(1-u)`, `G*(y)` = Stirling remainder (`log G*(y) ~ Sum_k b_k y^(1-2k)`, `b_k = B_2k/(2k(2k-1))`).

## 0. Bottom line

1. **The naive reduction is vacuous (PROVED).** `ph_n(0)` equals the model coefficient for EVERY n (the `u^0` part of the recursion is the
   model recursion, pivot `4n != 0`). The model never constrains the free coefficient `f_m := [u^(4m)] ph_(5m)`; s.2.
2. **The correct reduced criterion (PROVED, s.3).** With the triangular first-row expansion `Phi/C0 = Sum_k w^k g_k(u) F_(-k)(tau w)`,
   `[z^5q] P_11(0,z) = [u^4q](g_(5q) + g_(5q-1))`. So (N) is one explicit linear identity per `n = 5q`, involving `ph_0..ph_5q`, `varpi_0`,
   `tau` and the model rows to order `u^4q`. It reproduces `m5 = -43397/807975000` (= `p'_(4,1)(0)` = `-p^can_(5,0)(0)`) and E4.
3. **New evidence (VERIFIED, exact over Q, frame (0,0)).** `D == 1 mod z^71`: (H) for all `n <= 70`, `5 not| n` (56 identities) and
   `c^act = 1` at `z^5..z^70` (14 identities; `a_5 = a_7 = .. = a_13 = 0`). Hence (F_g) holds for `g <= 24` (was `g <= 9`). s.4.
4. **One integrand carries both normalisations (s.5).** The Mellin-Barnes summand of the I-function, Stirling-expanded, is in `x = 1/sigma`
   `exp(-(Lam/x + G(x))/z) at(x,z) dx`, `G(x) = Sum_(m>=1) x^(5m-1)/(m(5m-1))` (PROVED), amplitude independent of `u`. Its formal
   stationary phase (SP) at the critical point is the QRR-normalised WKB function (VERIFIED: all coefficients, `w^0..w^25`, two primes,
   including `f_1..f_5`); its formal residue at the pole `x = 0` is the quintic vanishing solution `varpi_0 (Y^C_Q)_0` (VERIFIED, `u^1..u^14`,
   exact). In hatted variables the phase is `W^(X;u) = Sum_m (u^(m-1)/m)[1/X + X^(5m-1)/(5m-1)]`, critical points `X^5 = 1` fixed.
5. **Theorem P (PROVED modulo Theorem U and GI Theorem G, s.6).** `Phi^SP = Phi^can`: the SP normalisation is the canonical conifold
   normalisation. The proof is formal (Chester-Friedman-Ursell change of variables, reduction modulo exact forms, residue = vanishing).
6. **Verdict.** D_0 = 1 is REDUCED to the single large-radius identity **(A2): `Phi^SP(z, u=1) = h_(w_0)(z)`** ("QRR = stationary
   phase of the MB summand at q = 0"), which no longer involves the conifold. (A2) is VERIFIED directly through `w^25` (two primes) and,
   via Theorem P and item 3, through `w^70` (exact). (A2) is NOT proved here (s.7).

## 1. HPF and the recursion (PROVED; re-derived from the PF operator, not copied)

`D = q d/dq = -(1-u) d/du`, `3125 q = 1-u`, `s = u^(1/5)`, `z = s^4 w`. On `F = s^e Fh(w,u)`: `(zD + L)F = s^(e-1) A_e Fh`,
`A_e = 1 + w(1-u)((4/5)theta_w - e/5 - theta_u)` (chain rule; `L = 1/s`, colour 0). With `psi = s^-1 Phi`, `Phi_j = A_(-j) Phi_(j-1)`, the PF
operator `u M^5 - 1 - (1-u)(3zM^4 + (17/5)z^2M^3 + (9/5)z^3M^2 + (274/625)z^4M + (24/625)z^5)` (the `e_k(1/5,..,5/5)`) gives, after dividing by `s^-1`,

    (HPF)  Phi_5 - Phi_0 - (1-u)[3 w Phi_4 + (17/5) u w^2 Phi_3 + (9/5) u^2 w^3 Phi_2 + (274/625) u^3 w^4 Phi_1 + (24/625) u^4 w^5 Phi_0] = 0

(`z^m M^(5-m) psi = s^-1 u^(m-1) w^m Phi_(5-m)`). Coefficients: `[w^m]Phi_j = [w^m]Phi_(j-1) + (1-u)((4(m-1)+j)/5 - theta_u)[w^(m-1)]Phi_(j-1)`,
`[w^m]Phi_0 = ph_m`. Every raising carries a left factor `(1-u)`, so `[w^(n+1)] HPF = (1-u) [ (4n - 5 theta_u) ph_n + r_n ]` with

    r_n = Sum_(j=1..5) ((4n+j)/5 - theta_u) Phi^b_(j-1,n) - 3 Phi^b_(4,n) - (17/5) u Phi_(3,n-1) - (9/5) u^2 Phi_(2,n-2) - (274/625) u^3 Phi_(1,n-3) - (24/625) u^4 ph_(n-4),

`Phi^b_(j,n)` = `[w^n]Phi_j` computed with `ph_n := 0` (the `ph_n` part sums to `(4n + 3 - 5theta) ph_n - 3 ph_n`). So
`[u^k] ph_n = [u^k] r_n / (5k - 4n)` for `5k != 4n`; for `n = 5m` the coefficient `[u^4m] r_n` must vanish (resonance) and `f_m = [u^4m] ph_5m` is free.
VERIFIED (`pp_core.py`): `ph_n(0)` = UR's `F_0` (`1, -3/20, 9/800, 11571/80000, ..`), and for the QRR choice `ph_n(1) = h_n` holds automatically
for all `5 not| n`, `n <= 70` (the referee's `h_auto`, now exact over Q and 70 orders). Example: `ph_5 = 47735919/128000000 - (203511/250000) u +
(548547/1000000) u^2 - (6789/62500) u^3 + f_1 u^4`, `h_5 = -257071/8064000000`, so `f_1 = 4537/3937500`.

## 2. Polynomial structure

* **Lemma 1 (PROVED).** `ph_n(0)` is fixed by `ph_0(0), .., ph_(n-1)(0)` for every `n >= 1` (constant term of `(4n - 5 theta)ph_n` is `4n ph_n(0)`;
  `theta_u` kills constants and the `u^(m-1)` brackets vanish at `u = 0`): the model recursion of UR s.2. **So "ph_n(0) = model" carries no
  information about the normalisation**: the free term `f_m u^4m` vanishes at `u = 0`.
* **Lemma 2 (PROVED).** A solution of HPF with `[w^n]Phi in Q[[u]]` has polynomial coefficients; the solutions are `Phi^(0) c(z)`,
  `c in 1 + z^5 Q[[z^5]]`, `z^5m = u^4m w^5m`, where `Phi^(0)` has all `f_m = 0`. (The operator `4n - 5theta` is diagonal on `u^k`; `r_n` is a
  polynomial if the lower `ph` are; multiplication by a function of z commutes with HPF.) The resonance conditions hold for every choice
  of the `f_m` as soon as ONE power-series solution exists (Theorem P(a) or Theorem U provide one). VERIFIED for `n <= 70`.
* **Lemma 3 (VERIFIED, n <= 70; not proved).** `deg ph_n <= floor(4n/5)` (equality except `n = 2`). Hence `f_m` is the LEADING coefficient:
  in original variables the `z^n` amplitude is `u^(-(4n+1)/5) ph_n(u)`, so `F_inf(z) := Sum_m f_m z^5m` is the termwise limit of
  `u^(1/5) x amplitude` at `u = infinity` (the orbifold point q = infinity). Large radius is `u = 1` (`z = w` there), the conifold `u = 0`.
* **Lemma 4 (PROVED).** For the QRR solution, `Phi^QRR = Phi^(0) F_inf^QRR(z)` and `F_inf^QRR(z) = h_(w_0)(z) / Phi^(0)(z, 1)` (evaluate at u = 1,
  where `z = w`). VERIFIED: `F_inf = 1 + (4537/3937500) z^5 - (177640111151/31007812500000) z^10 + ..` (log coefficients: `z^5`: `(1/1260)(4537/3125)`,
  `z^15` has the denominator of `B_16/240`, `z^10, z^20` only powers of 2 and 5).
* **The leading operator.** On polynomials `(4n - 5theta)` has kernel `Q u^(4n/5)` (5 | n) or 0. The `u = 0` datum is the constant term, the
  `u = 1` datum the coefficient sum: `ph_n(1) - f = Sum_k [u^k] r_n /(5k - 4n)` (`f = 0` if `5 not| n`). This map alone does not relate the
  model to the QRR constants; the relation needs the conifold gauge (s.3) or an integral representation (s.5-6).

## 3. What the conifold normalisation fixes (PROVED)

Hatted first row (UR s.4, D = 1 or any normalisation): `(1/C0) Phi = Sum_k w^k g_k(u) F_(-k mod 5)(tau w)`, `C0 = varpi_0`,
`g_k = u^mu(-k,k) tau^((1-c_b)/4) p_(k,b)(u)`, `b = -k mod 5`, `mu(-k,k) = (4k + c_(-k))/5 = 4q + (0,1,2,3,3)_r` for `k = 5q+r` (a first draft
used `4q + r`, which is wrong at `r = 4`; the divisibility check caught nothing because `p_(4,1)(0) = 0` by (R1), but `c` came out wrong).
`F_b(0) = (1,-1,-1,1,1)`, so `g_N` is determined triangularly from `ph_0..ph_N`; divisibility of `g_N` by `u^mu` is the regularity of `P`.
Row 0 of `z P' = A P - P At`: `P_1. = (z P_0.' + P_0. At)/b_0`; entry (1,1) with `b_0(0) = t_c'(0) = 1`, `bt_0(0) = v'(0) = 1`:

    c(z) := P_11(0, z) = P_00(0, z) + z (d/du) P_01(0, z),     [z^5q] c = p_(5q,0)(0) + p'_(5q-1,1)(0) = [u^4q] (g_5q + g_(5q-1))

(`p_(5q-1,1)(0) = 0` by (R1)). **(N) <=> `[u^4q](g_5q + g_(5q-1)) = 0` for all q >= 1, for the QRR polynomials.** This is the precise
"polynomial connection identity per n = 5q"; it is linear in `f_q` with coefficient 1 (`[u^4q w^5q](Phi/C0) = f_q + ..`).
**m5 explained.** `p_(4,1)` does not depend on the normalisation (only `z^k`, `k < 5`); `p'_(4,1)(0) = -43397/807975000`. The canonical
condition is `p_(5,0)(0) = -p'_(4,1)(0) = 43397/807975000` (= UR's `c_5`, found here from the QRR `f_1`); the stored normalisation
`p_(5q,0)(0) = 0` gives `c = 1 + p'_(4,1)(0) z^5 + .. = 1 + m5 z^5 + ..`. VERIFIED exactly (control in every `pp_run_norm` log).
Equivalent form (Lemma 4): (N) <=> `F_inf^QRR = F_inf^can := 1/c^(0)(z)`, `c^(0)` the `c` of `Phi^(0)` (control: VERIFIED to `z^70`).

## 4. Evidence: `c^act` through `z^70` (VERIFIED)

`pp_run_norm.py Q`: QRR polynomials to `w^5Q`, first-row extraction mod `u^(4Q+3)`, exact over Q, frame (0,0). Q = 14 (121 s):
(H) 56/56; `deg ph_n <= floor(4n/5)`; all divisibility checks; `[z^5q] c^act = 0` for `q = 1..14`; controls: `Phi^(0)` gives `c = 1/F_inf^QRR`,
the stored normalisation gives `m5` and `c != 1`. `p^act_(5q,0)(0)` agrees with UR/GI (`q = 1`: 43397/807975000). So `D == 1 mod z^71`, i.e.
(F_g) for every `g <= 24`, conditional (like E4) on (U) and the UR conventions, which this run re-tests through the reproduced numbers.
Not done: a second conifold frame and a second prime-free re-implementation of the first-row formula (the formula is two lines, s.3).

## 5. The integral: one phase, one amplitude, two functionals

**Phase (PROVED).** The colour-0 I-function is `Sum_d q^s Gamma(1+5s)/prod_j Gamma(1+s-zeta^j/z)`, `s = (zeta^0 + dz)/z`. With
`Gamma(1+y) = sqrt(2pi) y^(y+1/2) e^-y G*(y)`, `sigma = s z`: the exponents give `S(sigma)/z`, `S = sigma log q + 5 sigma log 5 sigma -
Sum_j (sigma - zeta^j) log(sigma - zeta^j)` (the `log z` terms cancel: `5 sigma = Sum_j (sigma - zeta^j)`). For `x = 1/sigma`, using
`Sum_j (1/x - zeta^j) log(1 - zeta^j x) = Sum_m x^(5m-1)/(m(5m-1))` (cyclotomic sums), `S(1/x) = -W(x;u)`, `W = Lam/x + G(x)`. Critical points
`dW/dx = -(Lam + log(1 - x^5))/x^2 = 0 <=> x^5 = u`. Hatted (`x = sX`, `z = s^4 w`): `W = s^4 W^(X;u)` with W^ as in s.0; `W^'(1) = 0`,
`W^(1) = (5/4)/tau`, `W^''(1) = 5/(1-u)`, `[t^k] W^(1+t) = Sum_m u^(m-1)[(-1)^k + C(5m-1,k)/(5m-1)]/m` rational in u, poles only at u = 1.
**Amplitude (PROVED).** Summand x `dsigma/z` = `alpha(x,z) e^(-W/z) dx`, `alpha = -(2pi)^-2 sqrt5 z at`, `at = (1-x^5)^(-1/2) G*(5/(xz)) /
prod_j G*((1-zeta^j x)/(xz))`; `log at = -(1/2)log(1-x^5) + Sum_k b_k (xz)^(2k-1)[5^(1-2k) - 5 Sum_P C(2k-2+5P,5P) x^5P]`. It does not depend on
u; in hatted variables it is a power series in u (`x z = u X w`), and `at(0, z) = 1` exactly.
**SP (VERIFIED).** `Phi^SP := SP_(X=1)[e^(-W^/w) B^ dX] / sqrt(2pi w/5)` (`B^` = `at` in hatted variables; the prefactors
`(2pi)^(-3/2) z^(3/2) e^(u_0/z) L` match those of the `d = 0` term at q = 0 exactly). `pp_run_sp.py 25 22`: `[w^n]Phi^SP == ph_n^QRR` in all
coefficients `u^0..u^21`, `n <= 25`, at 2^61-1 and 2^89-1, including `f_1..f_5`; vanishing above `u^floor(4n/5)`.
**Residue (VERIFIED).** `V := Res_(x=0) e^(-W/z) at dx = Sum_(k>=1) (-Lam/z)^k/k! [x^(k-1)](e^(-G/z) at)`. `pp_res.py 14`: `V = varpi_0 (Y^C_Q)_0`
in every coefficient `u^1..u^14` (exact; `Y^C_Q` built independently from the periods by Lemma Q of GI). E.g. `[u^2]`: `-1/(2z)` from `Lam` and
`-1/(5z)` from the Stirling term `b_1 = 1/12` give `-7/(10z) = -gamma/z`: the Bernoulli number produces the conifold period.
**CFU (PROVED; VERIFIED k < 10).** There are unique `Y = X(1 + Sum_(k>=1) u^k Y_k(X^5))`, `Y_k in Q[X^5]`, and `V(u)` with
`W^(X;u) = V/Y + Y^4/4`; order `u^k` times X reads `(X^5 - 1) Y_k = R_k(X^5) - V_k` (solvable iff `V_k = R_k(1)`), and the critical value forces
`V = ((4/5) W^(1))^(5/4) = tau^(-5/4)`, i.e. `v = u tau^(-5/4)` (UR's pull-back, now from the integral). `Y_1 = (1 + X^5)/18`.

## 6. Theorem P: the stationary-phase normalisation is the canonical one

**Theorem P.** (a) `Phi^SP` solves HPF and `[w^n]Phi^SP` is a polynomial in u (PROVED). (b) `Phi^SP = Phi^can`, the hatted first row of
`P_can Q^m(v(u))` (PROVED modulo Theorem U (UR s.3) and GI Theorem G; the formal lemmas used are standard and written in outline).
*Proof of (a).* `PF[q^s phi(s)] = g(s) - g(s+1)`, `g(s) = q^s ((zs)^5 - 1) phi(s)`, by `((z(s+1))^5 - 1) phi(s+1) = z^5 prod_k (5s+k) phi(s)`.
With the Stirling-expanded `phi` (asymptotic expansion of an exact identity, so the recursion holds formally) and `sigma = zs`,
`g(sigma) - g(sigma + z) = -d/dsigma Sum_(k>=1) z^k g^(k-1)/k!`; the formal SP of a total derivative is 0 (Gaussian integration by parts).
So `SP_(sigma=L)` is a formal solution of PF with exponent `e^(u_0/z)`, and `Phi^SP` solves HPF. In hatted variables every ingredient is in
`Q[[u]]` (s.5: `[t^k]W^(1+t)`, `1/W^''(1) = (1-u)/5`, `B^`), so `[w^n]Phi^SP in Q[[u]]`, a polynomial by Lemma 2. QED
*Proof of (b).* 1. (CFU, s.5) `Y(X;u)` with `W^ = V/Y + Y^4/4`, `V = tau^(-5/4)`; it maps `X = zeta^-i` to the model critical points, `0` to `0`.
2. (Reduction) `B^ dX = beta(Y) dY`, `beta in Q[Y, w][[u]]` (polynomial in Y and w at each u-order). Modulo the exact forms
`d(Y^m e^(-f/w))`, `f = V/Y + Y^4/4`: `Y^(m+3) == m w Y^(m-1) + V Y^(m-2)`, so `beta dY == Sum_(m=-1..3) pi_m(u,w) Y^m dY`, `pi_m in Q[w][[u]]`.
3. (Functionals) Formal SP at a nondegenerate critical point and the formal residue at `Y = 0` (coefficient of `Y^-1`, at each u-order an
entire series in `1/w`) vanish on exact forms and are invariant under the formal change of variables of step 1 (for the residue: termwise
in u, `Res F(X + u h)(1 + u h') dX = Res F dX` since the correction is `Sum_j (u^j/j!) d(F^(j-1) h^j)`).
4. (Identification) The five SPs of `Y^m e^(-f/w) dY` are the model's hatted WKB rows (UR s.1-2: rows `(I_0, -I_-1, -I_3, I_2, I_1)`,
`x = sY`, critical point `Y_c = V^(1/5)`, `sqrt(5/f''(Y_c)) = tau^(1/4)`, `F_b(tau w)`, weights `tau^((1-c_b)/4)`); the five SPs of `B^ dX` at
`X = zeta^-i` are the colours of `Phi^SP` (`W^(zeta X) = W^(X)/zeta`, `B^(zeta X; w) = B^(X; zeta w)`). Since `R^(tau w)` is invertible, `(pi_m)` is,
up to the factor `C0` and the signs of step 4, the hatted first row of the gauge `P^SP := Q^SP (Q^m)^-1`.
5. (Regularity) By (a) and Lemma 2, `Phi^SP = Phi^can s(z)` with `s in 1 + z^5 Q[[z^5]]`; so `P^SP = P_can s(z)` is a regular gauge (Theorem U)
and GI Theorem G gives `P^SP Y^C(v(u)) = s(z) Y^C_Q`.
6. (Residue) Steps 2-4 applied to the residue: `V := Res_(x=0) e^(-W/z) at dx = kappa(u) . (P^SP Y^C)_0 = kappa(u) s(z) (Y^C_Q)_0` in
`Q((z))[[u]]`, `kappa` a function of u only (powers of s, tau, C0 from hatting).
7. (Order u^1) `V = -u/z + O(u^2)` exactly (only `k = 1`; `G = O(x^4)`, `at(0,z) = 1`) and `(Y^C_Q)_0 = -u/z + O(u^2)` (GI Lemma Q: at `u^1` only
`n = 1`). Hence `kappa(0) s(z) = 1` for all z: `s` is constant, `s = s(0) = 1`, and `Phi^SP = Phi^can`. QED
*Consistency.* Step 6 with `s = 1`, `kappa = 1` predicts `V = varpi_0 (Y^C_Q)_0` in all orders: VERIFIED through `u^14` (s.5). This tests the
conventions of steps 4-6 independently of the gauge computation of s.3.

**Corollary.** `D = 1  <=>  Phi^QRR = Phi^SP  <=>  (A2) Phi^SP(z, 1) = h_(w_0)(z)` (the QRR WKB solution is the unique one with value `h` at
u = 1, and `Phi^SP` is a polynomial solution). (A2) contains (H) (`5 not| n`) and (N) (`5 | n`). It is a statement at q = 0 alone.

## 7. What remains: (A2), "QRR = stationary phase at q = 0"

* **Status.** VERIFIED directly for `n <= 25` (two primes, `pp_run_sp`), and for `n <= 70` via Theorem P + s.4 (exact). Not proved.
* **Equivalent forms.** `F_inf^SP(z) = h_(w_0)(z) / Phi^(0)(z, 1)` (Lemma 4; a u = 1 versus u = infinity identity for the polynomials
  with zero free coefficients); for `5 not| n` it is exactly the identity `ph_n^poly(1) = h_n` (the unique polynomial solution).
* **Mechanism (SKETCH).** As `u -> 1` the saddle `sigma = L -> 1` meets the colour-0 factor `1/Gamma(1 + (sigma-1)/z) = 1/d!`; the SP
  coefficients stay finite (they are polynomials) and should tend to the `d = 0` term, whose Stirling expansion is exactly
  `(2pi)^(-3/2) z^(3/2) e^(c/z) h_(w_0)(z)` (checked: the prefactors agree, s.5). Toy version: the formal SP of `Sum_d lambda^d/d!` at
  `d = lambda` is `e^lambda` with no corrections, and its value at `lambda = 0` is the `d = 0` term.
* **Routes.** (i) Laplace + Euler-Maclaurin for the exact `I_0(q)`, `0 < q < 5^-5`, `z -> 0+`, plus uniformity as `q -> 0` (analytic, allowed
  tools; the Poisson terms `k != 0` are exponentially small away from q = 0). (ii) A formal two-variable SP keeping the colour-0 Gamma factor
  exact through its Hankel integral. (iii) A literature theorem identifying the q = 0 R-matrix of the twisted theory with the SP of the
  mirror integral (not located or checked here; no statement number is claimed).
* **Why this is progress.** The conifold, the gauge, the periods and the model have disappeared from the open statement: what is left
  compares two normalisations of the same formal WKB solution at the large-radius point.

## 8. Files, runs, caveats

| file | content | run / log | result |
|---|---|---|---|
| `pp_core.py` | HPF recursion over Q, QRR constants (Newton sums of `1/(1-zeta^j)`), free coefficients | imported | - |
| `pp_series.py` | u-series over Q, periods `f0, f2` of the PF operator, `tau`, model rows `F_b` (own solver) | imported | - |
| `pp_gauge.py` | first-row extraction, `c(z)` (s.3) | imported | - |
| `pp_run_norm.py` | s.4 with controls | `logs/pp_run_norm_Q{4,7,14}.log` | PASS to `z^70` (121 s) |
| `pp_sp.py`, `pp_run_sp.py` | formal SP of the MB summand mod p | `logs/pp_run_sp_N25_K22.log` | PASS `w^0..w^25`, two primes |
| `pp_res.py` | residue versus `varpi_0 (Y^C_Q)_0` (own Lemma-Q iteration) | `logs/pp_res_14.log` | PASS `u^1..u^14` |
| `pp_cfu.py` | CFU map, `V = tau^(-5/4)` | `logs/pp_cfu_10.log` | PASS `k < 10` |

All runs `py -I -B`, stdlib, working directory in the session scratch folder, single processes (at most two at a time), none left running.
**Caveats.** AI-only. Theorem P rests on Theorem U and Theorem G (AI-audited) and on standard formal lemmas (SP and residue vanish on
exact forms and are invariant under formal changes of variables) that are stated, not proved in detail. s.4 uses one frame and relies on the
UR conventions (validated by reproducing `c_5`, `m5` and E4). Lemma 3 is only VERIFIED. (A2), hence D_0 = 1, is not proved.
