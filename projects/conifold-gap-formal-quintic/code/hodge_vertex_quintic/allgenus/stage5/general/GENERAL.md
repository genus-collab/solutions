# Stage 5, task `general`: a falsifier of the chart ARGUMENT for (Gamma) (2026-11-14)

Folder `allgenus/stage5/general/` (new files only). AI-written; no human and no referee has read it. Labels: PROVED (argument
written here, or quoted with statement number actually read), VERIFIED (exact computation, range and field stated), CONJECTURAL,
SKETCH, REFUTED. Notation as in the brief and `../../stage4/STAGE4_SYNTHESIS.md` (SYN4). Continuation: `GENERAL_2.md`.
Code (stdlib, `py -I -B`, run from the session scratch folder, copied here with their logs): `g5_core.py` (engine), `g5_data.py`,
`g5_theory.py`, `g5_valid.py`, `g5_t1.py`, `g5_t2.py`, `g5_t3.py`, `g5_t6.py`, `g5_t6b.py`, `g5_t6c.py`; logs in `logs/`.

## 0. Bottom line

* **Nothing is REFUTED.** Every prediction of the chart argument was confirmed, in every theory and colour tested; and in every
  tested case where a hypothesis of the argument fails, the conclusion fails too (poles), or there is no colour at all. Every
  failure in the logs is a deliberate control, and each control fails in the predicted WAY (pole versus wrong value).
* **(T1)** New method: `c(z; 1 + eps)` expanded directly as a Laurent series in `eps` at every order `z^n`. VERIFIED exactly over
  Q: no negative power of `eps` and `[eps^0] c_n = h_n` for all `n <= 60` (s.4; validated against an independent code, s.1).
* **(T2)** O(5)-twisted P^4 with distinct rational weights (four sets, one with `Sum lambda_j = 0`): the chart argument applies
  verbatim (the Euler substitutions leave the sigma-independent factor `z^(-1 - Sum lambda_j/z)`), and for EVERY colour the
  singularity at `sigma = lambda_i` is removable with value `h_(w_i)`, `w_i = {lambda_i - lambda_j} U {-5 lambda_i}`. VERIFIED
  exactly, `n <= 24`, 20 colour cases.
* **(T3)** P^1, P^2 (untwisted), local P^2 (two lifts), P^1 with O(-1)+O(-1) (two lifts), P^4 with O(3), O(4), O(6) and O(4) with a
  shifted lift: every prediction made from the twisted I-function holds (s.4.1), with the QRR signs: an Euler-class twist
  contributes minus its fibre weight, an inverse-Euler-class twist plus its fibre weight (a tangent weight of the total space).
* **The argument is local and more general than written** (GENERAL_2 s.5). Only ONE ingredient matters: the factor that
  degenerates at the colour point must be an exact `Gamma_St`, and it must be the only one. The other factors may be replaced by an
  ARBITRARY regular WKB function (random data, even powers of `z` included): the singularity stays removable, and the value is
  the regular part at the first pole of the singular Gamma factor. VERIFIED exactly (`n <= 20`), with a written SKETCH of the
  generalised chart argument. This gives a sharper and simpler statement than (Gamma) (s.3.3).
* **New prediction confirmed: the q = oo chart.** Solving the delta restriction for the twist variable instead of `x_0` predicts the
  value at `L = 0`: `c(z; 0) = Ginf(z)`, the CONJECTURAL closed form of Stage 4 (the shift 1/5 comes from the measure `dy`). VERIFIED
  exactly for `n <= 20`, also for general weights (`exp(Sum_n beta_n(1/5) z^n Sum_j lambda_j^-n)`). So Ginf is the same mechanism.
* **Where it must fail, it fails.** Coincident weights (two or three singular factors): poles of order exactly `n` at `z^n`; for
  two, the leading pole coefficients are Hankel's coefficients of `K_0`. Two singular factors of different slope (`lambda_i = 0` at
  `q = oo`): poles. Elementary and Bernoulli parts from different Gammas: poles. A vanishing twist weight: no colour exists.
* **(T5)** The critic's 5-variable Wick test of the chart step, `cr4_chart.py 6`: ALL CHECKS PASSED at `2^61-1` and `2^89-1` (742 s).
* **Status of (Gamma)** is unchanged by this task: CONJECTURAL with a SKETCH proof. What changed: the argument has now been tested
  as an argument (its hypotheses one at a time) and its general form is VERIFIED in every case tried.

## 1. The method of T1: Laurent expansion in eps (PROVED)

*Setting.* A Mellin-Barnes summand `prod_f Gamma(e_f + (c_f sigma + d_f)/z) Q^(sigma/z)`, `sigma = z s`, every Gamma replaced by
`Gamma_St(e + Y) = exp((Y + e - 1/2) log Y - Y + (1/2) log 2pi + Sum_(n>=1) beta_n(e) Y^-n)`,
`beta_n(e) = (-1)^(n+1) B_(n+1)(e)/(n(n+1))` (DLMF 5.11.8, read 2026-11-14; for `e = 0, 1` this is 5.11.1 with `beta_(2k-1) = b_k`,
`beta_(2k) = 0`). `c(z; sigma_c)` = formal SP at a critical point `sigma_c`, divided by its leading term (Gaussian factor,
`exp(S_0/z)`, `exp(A_0)` and constants), so `c_0 = 1`; the Stirling terms at `sigma_c` are kept (as in Stage 4).

*Observation.* No derivative of order >= 2 of the phase involves `Q`. So `c(z; sigma_c)` is a function of `sigma_c` alone, and each
`c_n` is a rational function of `sigma_c`. The colour point is `s*` where one factor, the singular one, has `c_* s* + d_* = 0`.

**Lemma 1 (scaling; PROVED).** Put `sigma = s* + eps (1 + T)`, `sigma_c = s* + eps`, `nu = z/eps`, `U = 1 + T`. Then
* the phase times `nu`, minus its value and linear part at `T = 0`, is
  `F(T; eps) = M c_* [U log U - T] + Sum_(m>=2) Phi_m eps^(m-1) [U^m - 1 - m T]`, with `Phi(s*+x) = Sum Phi_m x^m` the regular phase
  (`M` = number of singular factors; `M = 1` except in the degenerate tests);
* the amplitude is `g = C(eps) U^(M(e_* - 1/2)) exp(M Sum_n beta_n(e_*) c_*^-n nu^n U^-n) B(nu eps, eps U)`, with
  `B(z, x) = exp(A_0(s*+x) - A_0(s*) + Sum_n z^n A_n(s*+x))` built from the regular factors and `C = 1/B(0, eps)`;
* every coefficient of `F` and `g` (in `T`, resp. `T` and `nu`) lies in `K[[eps]]`, and `a := -F''(0) = -M c_* - eps Phi''(s*+eps)`
  is a unit of `K[[eps]]`;
* the normalised SP in `T` with parameter `nu`, `ct(nu, eps)`, satisfies `ct(z/eps, eps) = c(z; s* + eps)` term by term.

*Proof.* `sigma - sigma_c = eps T` is linear. The phase terms `(c sigma + d) log(c sigma + d)` of the singular factor give
`(c_*/nu) U log(c_* eps U/z)`; the part `(c_*/nu) U log(c_* eps/z)` is linear in `T` and drops out. Every regular term is a power
series in `x = eps U`. A Wick term of order `m` is a product of Taylor coefficients times `<t^(2m)> = (2m-1)!! (z/A)^m`, where
`A = -S_0''(sigma_c) = a/eps`; and `eps^(2m) <T^(2m)>_nu = eps^(2m) (2m-1)!! (nu/a)^m = (2m-1)!! (z/A)^m`. The Taylor coefficients
scale the same way (`d^k/dT^k = eps^k d^k/dsigma^k`), so the two expansions agree term by term. □

**Corollary.** `c_n(s* + eps) = eps^-n ct_n(eps)` with `ct_n := [nu^n] ct in K[[eps]]`. So the pole order at `eps = 0` is at most `n`
(a priori), the singularity is removable iff `ct_n = O(eps^n)` for all `n`, and then `c_n(s*) = [eps^n] ct_n`.

*Computation* (`g5_core.engine`): Lagrange-Buermann. With `F = -(a/2) T^2 Theta(T)`,
`ct = Sum_m (2m-1)!! (nu/a)^m [T^(2m)] g Theta^(-(2m+1)/2)` (Morse substitution `W = T Theta^(1/2)`). Writing `g` in powers of `U`
reduces everything to `L_m(e) = [T^(2m)] U^e Theta^(-(2m+1)/2)` for `2N + 1` exponents `e`. Cost `O(N^5)` operations in `K`.
Exact over Q, or mod p. This is not one of the Stage-4 methods (interpolation in `L`, SP at fixed `L` or `u`, WKB recursions, the
5-variable Wick test): it works AT the singular point and shows the cancellation of the poles directly.

*The layers of the cancellation.* At `eps = 0` the regular part drops out of `F` and `g` (its terms carry `eps^(>=1)`), so the
leading pole coefficients `ct_n(0)` are those of the singular factor ALONE: `ct_n(0) = 0` for all `n >= 1` is the one-Gamma
identity (formal Mellin inversion of `Gamma`). VERIFIED (`g5_t6c`, L0): `ct_n(0) = 0` for `1 <= n <= 20`; without its Stirling series the
values are `1/12, 25/288, 9041/51840` (so at `z^1` the residue is `1/12 - b_1 = 0`). The lower poles `eps^(-n+1), .., eps^-1` of
`c_n` involve the regular part; their cancellation is the content of (Gamma) and of s.3.3.

*Validation* (`g5_valid.py 6 24`, log `logs/g5_valid_6_24.log`; independent fixed-point SP code, no `eps`):
* quintic, `eps`-truncation `5N = 30`: the engine's series vanishes beyond `eps^(5n)` at `z^n`, `n <= 6` (it must: `c_n in Q[L]` of
  degree `<= 4n`, PROVED in Stage 4); the resulting polynomial equals the fixed-point SP at `L = 2, -2/3, 3/2`; and
  `c_5(2) = 11495493197/31500`, the value printed in SYN4 s.3 (E3);
* twisted P^4, weights W1, colour 0 (non-polynomial `c_n`): the truncated series at `eps = 1/1000` agrees with the fixed-point SP
  to relative `10^-27` or better, `n <= 6`.

## 2. The summands from the twisted I-functions (PROVED, elementary; VERIFIED by residues)

*Input (standard, quoted; statement numbers not checked).* Equivariant `P^(N-1)`, weights `lambda_j`; at the fixed point `p_i`
(`H = lambda_i`) the degree-`d` coefficient of the I-function is `M_d / prod_j prod_(k=1..d) (lambda_i - lambda_j + k z)`, with
`M_d = prod_(k=1..ld) (l lambda_i + mu + k z)` for a convex `O(l)` twisted by the Euler class, `c_1 = lH + mu` (for O(5), `mu = 0`
this is the project's formula, `FORMAL_QUINTIC_genus2.md` s.1), and `M_d = prod_(k=0..ld-1) (-l lambda_i + mu - k z)` for a concave
`O(-l)` twisted by the inverse Euler class, `c_1 = -lH + mu` (local P^2: `l = 3`; resolved conifold: two copies of `l = 1`).

**Lemma 2 (PROVED).** The summand `a(s) = prod_j Gamma(lambda_j/z - s) prod_conv Gamma(1 + l s + mu/z) prod_conc Gamma(l s - mu/z) Q^s`,
`Q = (-1)^(N + l_conc) z^(Sum_f c_f) q` (`c_f = -1` for the base factors, `c_f = l` for each twist; `l_conc` = sum of the
concave `l`) has `Res_(s = lambda_i/z + d) a / Res_(s = lambda_i/z) a = q^d I_d(p_i)` for all `d >= 0`.
*Proof.* At `s = lambda_i/z + d`: `Gamma(lambda_i/z - s)` has residue ratio `(-1)^d/d! = (-1)^d z^d / prod_(k<=d) (k z)`;
for `j != i`, `Gamma(x - d)/Gamma(x) = (-1)^d z^d / prod_(k<=d)(lambda_i - lambda_j + k z)` with `x = (lambda_j - lambda_i)/z`; a convex
factor gives `prod_(k=1..ld)(l lambda_i + mu + k z) z^-ld`; a concave one `prod_(k=0..ld-1)(l lambda_i - mu + k z) z^-ld`
`= (-1)^(ld) z^-ld M_d`. The signs `(-1)^(Nd)` and `(-1)^(l_conc d)` and the powers `z^((N - Sum l) d)` are cancelled by `Q^d`. □
VERIFIED (`residue_check`, exact, `z = 1/7`, `d <= 3`, every colour) for every theory of T2 and T3: all ratios equal 1.

**Lemma 3 (z-powers; PROVED).** Write each factor as an Euler integral and rescale `x_f -> x_f/z`. The factor `Gamma(e + Y_f)`,
`Y_f = (c_f sigma + d_f)/z`, contributes `z^(-e - Y_f)`. The total exponent of `z` is `-Sum_f e_f - Sum_f d_f/z - (Sum_f c_f) sigma/z`;
`Q^s` contributes `+(Sum_f c_f) sigma/z`. So the only `z`-power left is the sigma-INDEPENDENT `z^(-Sum e_f - Sum d_f/z)`, common to both
orders of integration. For the twisted P^4 with `O(5)` it is `z^(-1 - Sum_j lambda_j/z)`: when `Sum lambda_j != 0` this is the only
change to the quintic bookkeeping (VERIFIED as a symbolic check per weight set, `g5_t2` lines (Z)). □
*Remark.* With the naive `(-q)^s` in a non-CY theory, `(Sum_f c_f) sigma log z / z` survives in the phase, and the critical point
in `sigma` depends on `log z`. Since `c` depends on `sigma_c` only, this changes only the relation between `q` and `sigma_c`.

## 3. The chart argument in general form, its hypotheses, its predictions

### 3.1 The argument (SYN4 s.1.4, for any Gamma-product summand; SKETCH, as there)
Factors `f` with data `(c_f, d_f, e_f)`; the singular factor `*` has `c_* s* + d_* = 0`; put `w_f := c_f s* + d_f` for `f != *`.
* Euler integrals (Lemma 3 for the `z`-powers): integrand `exp(Phi_tot/z) prod_f x_f^(e_f) dx_f/x_f dsigma`, with
  `Phi_tot = sigma G + Sum_f (d_f log x_f - x_f)`, `G = log Q + Sum_f c_f log x_f`. Linear in `sigma`; the amplitude is sigma-free.
* (a) `x` first, then `sigma`: Stirling (L3) for every factor, then the SP in `sigma`: this is `c(z; sigma_c)`.
* (b) `(sigma, x_*)` first: the delta lemma (L2) restricts to `G = 0`, i.e. `x_*^(c_*) = Q^-1 prod_(f != *) x_f^(-c_f)`. With the
  measure `dx_*/x_*` the Jacobian is the constant `1/|c_*|`; the insertion `x_*^(e_*)` becomes `prod_(f != *) x_f^(-c_f e_*/c_*)`
  times a constant. Since `-d_*/c_* = s*`, the chart integrand is
  `exp([Sum_(f != *) (w_f log x_f - x_f) - X_*(x)]/z) prod_(f != *) x_f^(h_f) dx_f/x_f`, `h_f := e_f - c_f e_*/c_*`,
  `X_* = kappa prod_f x_f^(-c_f/c_*)`, `kappa = Q^(-1/c_*)`. For the base singular factor `kappa = Q`; at `q = oo`, `kappa = Q^(-1/5)`.
* (c) The chart critical point `x_f = w_f + O(kappa)` is analytic in `kappa`, and non-degenerate at `kappa = 0` iff every `w_f != 0`
  (Hessian `diag(-1/w_f)`). At the critical point `x_* = c_* eps`, so `eps <-> kappa` is analytic and invertible near 0.
  Fubini (L1) identifies (a) and (b) for `eps != 0`; hence each `c_n(s* + eps)` is analytic at `eps = 0`.
* (d) At `kappa = 0` the chart is `prod_(f != *) Gamma(h_f + w_f/z)`. **Prediction (P):** removable singularity at `s*`, value
  `c(z; s*) = prod_(f != *) exp(Sum_n beta_n(h_f) (z/w_f)^n)` (the regular factors at the first pole `sigma = s* - e_* z/c_*`).
* Base singular factor (`c_* = -1`, `e_* = 0`, the colour-`i` point `s* = lambda_i`): `h_f = e_f`, and (P) is `h_(w_i)` with
  `w_i = {-w_f}`: tangent weights `lambda_i - lambda_j`; convex `O(l)`: `-(l lambda_i + mu)`; concave `O(-l)`: `mu - l lambda_i`.

### 3.2 Hypotheses, and what happens without them
* (C1) exactly one factor degenerates at `s*`; (C2) it is an exact `Gamma_St(e_* + Y)` (elementary part and Bernoulli part of the
  SAME Gamma); (C3) `w_f != 0` for every other factor; (C4) the formal lemmas (L1)-(L3).
* Not used: the CY condition (Lemma 3), `Sum lambda_j = 0`, `|c_*| = 1`, `e_* = 0`, the sign of `q`, and (s.3.3) the Gamma form of
  the regular factors.
* (C1) fails for coincident weights (`lambda_j = lambda_i`: the chart variable `x_j` has `w_j = 0`) and for two singular factors of
  different slope (`lambda_i = 0` at `q = oo`). (C2) fails for perturbed Bernoulli coefficients or a stray factor `(sigma - s*)^alpha`.
  (C3) fails at `q = 0` exactly when a twist weight vanishes; then (C1) fails as well, and no colour exists:

**Lemma 4 (PROVED).** Let `lambda_j != lambda_i` for `j != i`, and let the argument of the twist factor `Gamma(e + l s + d/z)`
vanish at `sigma = lambda_i` (`l lambda_i + d = 0`, `l >= 1`). Then no critical point tends to `lambda_i` as `q -> 0`.
*Proof.* Near `sigma = lambda_i` the critical equation reads `log(+-Q) + l log(l(sigma - lambda_i)) - log(lambda_i - sigma) + O(1) = 0`
(the other factors are regular and non-zero there), i.e. `(l - 1) log|sigma - lambda_i| = -log|Q| + O(1)`. For `l >= 2` the right side
tends to `+oo`, contradicting `sigma -> lambda_i`; for `l = 1` it reads `log|Q| = O(1)`, impossible as `Q -> 0`. □

### 3.3 The general form: local Mellin inversion (LMI)
**Statement** (CONJECTURAL; VERIFIED in GENERAL_2 s.5). Let `B(sigma; z) = exp(Phi(sigma)/z + Sum_(k>=0) z^k A_k(sigma))` with `Phi`,
`A_k` power series in `sigma - s*`, and let `c(z; sigma_c)` be the normalised SP of `Gamma_St(e_* + c_*(sigma - s*)/z) B Q^(sigma/z) dsigma`.
Then `c_n(s* + eps) in K[[eps]]` for every `n`, and `c(z; s*) = exp(Lam(z) - Lam(0))`, with
`Lam(z) = [Phi(s* + delta) - Phi(s*)]/z + Sum_(k>=0) z^k A_k(s* + delta)`, `delta = -e_* z/c_*`:
**the normalised value of the regular part at the first pole of the singular Gamma factor.**
*SKETCH of proof.* Represent `B` as a formal Laplace integral, `B(sigma) = SP Int exp(sigma t/z) Bhat(t; z) dt` (`Bhat` the formal
Legendre-Laplace dual; it exists if `Phi''(s*) != 0`), and the singular factor by its Euler integral. The phase is linear in `sigma`.
The delta lemma restricts to `c_* log x = -(log Q + t)`, i.e. `x(t) = kappa e^(-t/c_*)`. The remaining `t`-integral has the phase
`s* t - x(t) + Phihat(t)` and the insertion `x(t)^(e_*) = kappa^(e_*) e^(delta t/z)`, analytic in `kappa`. At `kappa = 0` it is the
Laplace integral of `Bhat` at `sigma = s* + delta`, and Laplace inversion gives `B(s* + delta)`. If `Phi''(s*) = 0`: every coefficient
`[nu^n eps^r] ct` is a polynomial in the Taylor coefficients of `Phi` and `A_k` (`a = -c_* - eps Phi''` is a unit for every `Phi`),
so the identity extends by Zariski density. For a Gamma-product `B` this is s.3.1: `Gamma(e + Y) = Int exp((e + Y) u - e^u) du` is a
Laplace integral in `u = log x`, and `sigma` enters only through `Sum_f c_f u_f`.
*Lemmas needed:* (L1), (L2), (L3) for the singular factor only, and (L4) formal Laplace inversion for WKB functions. None is
written out here. *What it gives:* (Gamma) for the formal quintic is the case where `B` is the five regular Gamma factors (value
`h` by (L3) for them). Regularity at `q = 0` is a property of the singular factor alone.

## 4. Results of T1-T4 (all exact over Q unless stated)

| test | content | range | result | log (`logs/`) |
|---|---|---|---|---|
| T1 | quintic colour 0: (R) no `eps^-k` in `c_n(1+eps)`; (V) `[eps^0] c_n = h_n` | `n <= 40`, trunc. `eps^40` | PASS (105 s) | `g5_t1_40_40_Q.log` |
| T1 | the same | `n <= 60`, trunc. `eps^60` | PASS (750 s) | `g5_t1_60_60_Q.log` |
| T2 | twisted P^4, W1 = (1,3,-2,5/2,-7/3), W2 = (2,-1,3,-5,1) (sum 0), W3 = (1,1/2,1/3,1/4,1/5), W4 = (7,-3,4/3,-11/2,2/5); all 20 colours; (Z), (I), (R), (V) | `n <= 24` | PASS (330 s) | `g5_t2_24.log` |
| T3 | 12 theories, 40 colours (s.4.1) | `n <= 20` | PASS (258 s) | `g5_t3_20.log` |

**4.1 T3: predictions (made from the I-function before the runs) and outcome.** Weights rational; `w_i` as listed in the log.
* P^1 (weights (1,-2) and (3,1/2)) and P^2 ((1,3,-2) and (2,-1,1/3)), untwisted, `Q = (+-) z^-N q`: `w_i = {lambda_i - lambda_j}`.
* local P^2, `Gamma(3s - mu/z)`, `Q = (+-) q`, lifts `c_1 = -3H` and `-3H + 1/2`: `w_i = {lambda_i - lambda_j} U {mu - 3 lambda_i}`.
* P^1 with O(-1)+O(-1), `Gamma(s - mu_1/z) Gamma(s - mu_2/z)`, lifts `(0,0)` and `(1/3,-2/5)`: `w_i = {lambda_i - lambda_j, mu_1 - lambda_i, mu_2 - lambda_i}`.
* P^4 with O(k), k = 3, 4, 6, `Gamma(1 + ks)`, `Q = -z^(k-5) q`, weights W1: `w_i = {lambda_i - lambda_j} U {-k lambda_i}`; and O(4) with
  `c_1 = 4H + 2/3`: twist weight `-(4 lambda_i + 2/3)`. For O(6) one critical point runs to infinity as `q -> 0`; it is no colour.
* Outcome: (I), (R), (V) PASS for all 40 colours; each control (twist weights with the opposite sign in h, resp. one tangent weight
  + 1/11 in h) fails at `z^1`. The QRR sign rule (Euler class: minus the fibre weight; inverse Euler class: plus the fibre weight)
  is the one the chart predicts; I did not re-derive it from the QRR papers.

**4.2 T4 controls** (quintic: `g5_t1`, `n <= 12`; T2: W1, every colour, `n <= 10`). Each fails in the predicted WAY.

| control | prediction | observed |
|---|---|---|
| twist weight `+5` instead of `-5` in h (quintic; T2 every colour) | (R) holds, value wrong | value fails at `z^1` |
| one weight `lambda_k + 1/11` in h only (T2, every colour) | (R) holds, value wrong | value fails at `z^1` |
| Stirling series of `Gamma(1+5s)` dropped | (R) holds (regular factor), value wrong | value fails at `z^1` |
| Stirling series of the singular factor dropped | pole | pole at `z^1` |
| `b_1 + 1/7` in every factor (Stage-4 control) | pole (singular factor perturbed) | pole at `z^1` |
| `b_3 + 1/7` in the singular factor only | pole | pole at `z^5` |
| `b_3 + 1/7` in the regular factors only | (R) holds; value = perturbed regular tail | as predicted, exactly, `n <= 12` |
| positive control: `Gamma(5s)` for `Gamma(1+5s)` (`dy/y` measure) | (R) and (V) still hold | both hold, `n <= 12` |
