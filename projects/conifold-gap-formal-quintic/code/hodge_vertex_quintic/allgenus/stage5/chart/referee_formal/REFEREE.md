# Referee report: Stage 5 / chart, "the chart argument for (Gamma), written out in full" (2026-11-14)

Folder `allgenus/stage5/chart/referee_formal/` (task key `chart/referee_formal`). New files only. AI-written: no human has read this
report or the note it referees. Labels as in the brief: PROVED, VERIFIED, CONJECTURAL, SKETCH, REFUTED; in addition CONFIRMED means
"I checked this argument line by line and found it correct". Part 2: `REFEREE_2.md` (colours, Theorem D, corrections, code, runs).
Notation of the note: `Lam = C \ {0, 1, zeta^j}`, `q(L) = (1 - L^-5)/5^5`, `Lam1(L) = 1 + L + .. + L^4`, `h = h_(w_0)`.

**Read.**
* In full: `CHART_PROOF.md`, `_2`, `_3`, `_4`.
* The worker's logs: all 16 scanned. None contains a FAIL line; every control is printed as "PASS control"; each ends "ALL CHECKS
  PASSED". Five were read in full, and eight were compared line by line with my re-runs (REFEREE_2 s.7).
* The worker's code: `ch_phases.py` (Taylor data), the head of `ch_wick.py`, and the set-up of `ch_chart.py`.
* Stage 4: `MECHANISM_STAGE4.md`; `STAGE4_SYNTHESIS.md` (all, incl. s.1.2 and s.1.4) and `_2`; `GAMMA.md`; `POLY.md` s.2-s.8;
  `EVID.md` s.2.1-2.4; poly referee `REFEREE.md` s.2-s.5; the docstring of `critic/cr4_sp.py`.
* Not read: `GAMMA_PROOFS.md`, the gamma and evid referee reports, UR, GI (only as quoted in SYN4).
* Literature, re-read by me: DLMF s.5.11 and s.2.3 (statement numbers); arXiv:math/0108100v2 pp. 12-15 (s.9.1, s.10).

## 0. Verdict: SOUND, WITH (MINOR) CORRECTIONS

* **Theorem Gamma (i)-(iv): CONFIRMED.** The chart argument is a complete and correct formal proof. Its only non-formal inputs are
  DLMF 5.11.1 (Stirling's series) and DLMF 2.3.8 (Watson's lemma), used once, in (L3). I re-read both: the numbers and the
  statements are as quoted. No gap. In detail:
  * every formula, prefactor, power of `z`, Hessian and determinant of s.1.6 and s.3.1-s.3.6 is right (by hand, and exactly in code);
  * the hypotheses of (L1b), (L2), (L3) hold at every `L` in `Lam` in both orders, and nothing in the argument needs uniformity in `L`
    near `L = 1` (s.3 below);
  * the continuation to `L = 1` is an identity of rational functions, and it is correct;
  * order (a) gives exactly the `c(z; L)` of GAMMA s.2 = POLY s.5 = `cr4_sp`;
  * the twist weight in step (d) is `-5`, with the right sign;
  * the Z/5 covariance and the reduction of the other colours to colour 0 are correct.
* **Theorem D: CONFIRMED as an assembly.** `D_0 = 1` follows from Theorem Gamma (iii), s.3.10, R5 (Theorem P(b)), R3 with R2, and R0.
  The list of hypotheses in s.0.3 is complete. One addition: EVID's "u = 1 convention of ref_rows", which the note folds into (QRR).
  Theorem D is PROVED only relative to these Stage-4 inputs and (QRR). None of them is re-proved here, and none is human-checked.
* **Corrections** C1-C6 (REFEREE_2 s.6). They concern wording and precision; none changes a conclusion. The two that matter for a
  reader: C1, a false sentence ("det H_ch vanishes only at ..."; it never vanishes, it has poles there), and C2 (s.3.10: the
  Jacobian of `sigma = L/X` is not constant; it is part of POLY's amplitude, and only its value at `X = 1` drops out).
* **Independent evidence** (new code, stdlib, no import of the worker's or any Stage-4 module, failing controls; REFEREE_2 s.7):
  * (Gamma) directly at `L = 1`, with no assumption on the shape of `c_n`: `c_n` is expanded at `L = 1 + eps` in truncated Laurent
    series with rigorous precision tracking. All poles `eps^-n .. eps^-1` cancel and `c_n(1) = h_n`, for every `n <= 40` (exact, Q).
  * POLY's hatted integrand with `u` formal: `deg ph_n <= floor(4n/5)` and `ph_n(1) = h_n` for `n <= 30`; `L^(4n) ph_n(L^-5) = c_n(L)`
    (s.3.10); `ph_5` equals the QRR polynomial printed in GAMMA s.1.
  * The chart and the 7-variable integral against `c(z; L)`, with two new Gaussian-expectation algorithms (not Isserlis):
    * exact in `Q(zeta)`: through `z^3` for the chart and `z^2` for the 7-variable integral;
    * mod two new primes: through `z^5` for the chart (`L = -3, 5/3`; and `= h` at `L = 1`) and `z^3` for the 7-variable integral.
      `z^5` is the first order that sees (N).
    * All Hessian and determinant identities hold exactly.
  * The lemmas (L1a), (L1b) with the factor `beta`, (L2), (L3) and (S1), on random exact data; and step (e) through `q^6 z^12`.
  * Eight of the worker's scripts re-run from a scratch copy: the logs are identical to the stored ones, timings stripped.

## 1. Section 1 (definitions; the object c(z; L))

* **s.1.2 (Wick functional).** The range `V <= 2(n-m)` and the degree `2(n-m+V)` are right: a term `Phi_3^V A_m` of total degree `D`
  carries `z^(D/2 - V + m)`. (F1)-(F3) CONFIRMED; (F3) by `eta = mu^(-1/2) eta'`. (F4) is a remark, used in no proof.
* **s.1.6.** `e^lambda = prod_j (zeta^j - L)/(5L)^5 = (1 - L^5)/(5^5 L^5) = -q(L)` CONFIRMED. `S_0'`, `S_0'' = 5/(sigma(1 - sigma^5))`,
  the Taylor data and the generating function `Sum_k P_k t^(k-1) = 5(L+t)^4/(1 - (L+t)^5)` CONFIRMED by hand.
  * The displayed identity `e^(lambda sigma/z) Gamma_St(1 + 5 sigma/z) prod_j Gamma_St((zeta^j - sigma)/z) = (2 pi)^3 z^2 exp(S_0/z + S_1 + ...)`
    is CONFIRMED. The `(log z)/z` coefficient is `-(5 sigma + Sum_j (zeta^j - sigma)) = 0`, the `log z` coefficient is `-1/2 + 5/2 = 2`,
    and the constant is `(2 pi)^(6/2)`.
  * Theorem Gamma (i), `c_n in Q[L, L^-1, (1 - L^5)^-1]`: CONFIRMED (`P_k in (1 - L^5)^-k Q[L]`, `-1/S_0''(L) = -L(1 - L^5)/5`).
* **s.1.7 (identification with the Stage-4 objects).**
  * GAMMA s.2 uses the same integrand, measure `d sigma` and normalisation `c_0 = 1`. Identical.
  * POLY s.5 uses the `1/Gamma(1 + s - zeta^j/z)` form, `q^s` and the hatted variables.
    * The reflection argument is right: `log Gamma_St(y) + log Gamma_St(1-y) = log 2 pi + (y - 1/2)(log y - log(-y))`, the `b_k` terms
      cancel, and `e` is even and forced to be 0.
    * A shorter route: the two flat phases have the same second derivative, `S'' = 5/sigma - Sum_j 1/(sigma - zeta^j) = S_0''`. The
      normalised amplitudes agree, because `(zeta^j - sigma)^(1-2k) = -(sigma - zeta^j)^(1-2k)`. And `L` is critical for both by the
      choice of the logarithm.
    * Re-derived from the summand by me: the hatted phase `d/dX [L^4 S(L/X)] = (1/(uX^2)) log((1 - uX^5)/(1 - u))`, and the amplitude
      `(1 - uX^5)^(-1/2) exp(Sum_k b_k (uXw)^(2k-1) [5^(1-2k) - 5 Sum_P C(2k-2+5P, 5P) u^P X^(5P)])`. The latter includes the Jacobian
      `d sigma/dX` (see C2). They agree with POLY s.5 and with `cr4_sp.py`.
    * VERIFIED (`rf_hat.py`, (H3)): `L^(4n) ph_n(L^-5) = c_n(L)` for `n <= 12` at `L = 2, -2/3, 3`.
  * The claims on `gamma_interp.py`, `cr4_chart.py` (prefactor `2 pi kappa (-q)^(1/z)`, "Jacobian |x_0|" = `kappa X_0`), on `q` versus
    `-q` and on `dx/x` versus `dx` are CONFIRMED.
  * A remark on `dy/y`. With `Gamma(5s)` in place of `Gamma(1+5s)`, `c(z; L)` changes from `z^2` on at `L != 1`. At `L = 1` the series
    is exactly unchanged: by (L3), `Poch_0 = Poch_1 = 1`, i.e. `Gamma(1+w) = w Gamma(w)` changes only the elementary factor. VERIFIED
    to `z^6` (`rf_sigma`). The measure test is therefore meaningful only at `q != 0`, as the critic had noted.

## 2. Section 2 (the formal lemmas): each proof checked

| lemma | verdict | remarks | test (`rf_lemmas.py 7 6`; controls fail) |
|---|---|---|---|
| (S1) | CONFIRMED | `<(H eta)_i P + d_i P>_C = 0` from `HC = -1`; the `z^(1/2)` bookkeeping is right | `N[Phi, (d_v Phi) h + z d_v h] = 0`, 6 seeds, `z^7`; control `(d_v Phi) h` alone |
| (PG) | CONFIRMED | `dC/ds = C delta C`; `<QP> = (1/2) tr(delta C) <P> + (1/2) Sum (C delta C)_ij <d_i d_j P>`; `d log d(s)/ds = (1/2) tr(C delta)`; `I`-adic convergence | through (L1b) |
| (L0) | CONFIRMED | block-diagonal covariance | - |
| (L1a) | CONFIRMED | linear part plus homotopy `chi_t` tangent to the identity; `d_t Phi_t = O(|x|^3)`; Lie-derivative step; (S1) applies at every `t` | nonlinear `chi`, 6 seeds; control: Jacobian dropped |
| (L1b) | CONFIRMED | reduction to `v_c = 0` by a Jacobian-1 shift; `H = diag(H_vv, S)`; three exponents of degree `>= 3`; (PG) in `C[eta_w][[z^(1/2)]]` gives `beta` and `-H_vv(w)^-1`; inner and outer expansions commute (finitely many terms per order) | inner series over `Q[[w]]`, inner point critical identically in `w`, 6 seeds, `z^7`; controls: `beta` dropped; inner point frozen |
| (L1c) | CONFIRMED | (i)-(ii) from (F1). (iii'): translation to a fixed point, the infinitesimal (PG) for the Hessian, the `p'` terms killed by (S1). (iii) by iteration; each `t^d`-coefficient lies in `z^(n/2) C((z))` | through step (e) (s.3.7) |
| (L2) | CONFIRMED | `C_(xi xi) = 0`; each monomial `s^a xi^b` of `Phi_3^V A_m` has `a <= b/2`; a non-zero pairing needs `b <= a`; so `a = b = 0`. `det = -g_1^2`, critical value `f(X_0)` | `N = A(X_0, z)` exactly, random `G`, `f`, `A`, 6 seeds; controls: `s xi` in `A`, `s^2 xi` in the phase |
| (L3) | CONFIRMED mod DLMF 5.11.1 and 2.3.8 | scaling `x = at` + (F3) reduces it to rational identities; Stirling with `(w)_M`; Laplace via the Morse substitution and Watson's lemma (`mu = 2`, `lam = 1`; odd terms cancel between `u > 0` and `u < 0`; `2 Gamma(k + 1/2) 2^(k - 1/2) = sqrt(2 pi) (2k-1)!!`); uniqueness of asymptotic expansions | `M = -3..5`: exact over Q to `z^18` (`a = 1`) and `z^12` (`a = 7/3, -2/5`); over `Q(zeta)` at `a = zeta - 2` to `z^8`; controls fail at `z^3`, `z^1`, `z^1` |

*(L3), a simplification (optional).* The Pochhammer part is purely formal. For `phi = a log x - x`, (S1) with `h = x^M` reads
`N[phi, (a/x - 1) x^M + z M x^(M-1)] = 0`, i.e. `N[phi, x^M] = (a + Mz) N[phi, x^(M-1)]` for every integer `M`. This is exactly
`a^M Poch_(M+1)(z/a) = (a + Mz) a^(M-1) Poch_M(z/a)`. So the analytic input (DLMF) is needed for one value of `M` only, the Stirling
series `N[a log x - x, 1; a] = exp(Sum_k b_k (z/a)^(2k-1))` itself. For that remaining case the worker's remark in s.2.8 stands: a
fully formal proof would need a shift `a -> a + z`.

## 3. Section 3, steps (a)-(e)

### 3.1 Powers of z and prefactors (s.3.1, s.3.3)

* `Gamma(1 + 5 sigma/z) = z^(-1 - 5 sigma/z) Int e^((5 sigma log y - y)/z) dy` (measure `dy`, since `Gamma(1+w) = Int t^w e^-t dt`) and
  `Gamma((zeta^j - sigma)/z) = z^(-(zeta^j - sigma)/z) Int e^(((zeta^j - sigma) log x - x)/z) dx/x`. CONFIRMED.
* The six exponents sum to `-1 - (1/z)[5 sigma + Sum_j (zeta^j - sigma)] = -1`, identically in `sigma`. The inner Stirling factors give
  `z^(5 sigma/z + 1) prod_j z^((zeta^j - sigma)/z) = z^(+1)`. Together with the `z^-1` of `Omega_7`, the symbol of the order-(a) inner SP
  is exactly the Stirling-expanded summand of s.1.6. CONFIRMED.
* (`ds = d sigma/z` would add one more `z^-1` to `Omega_7` if one insisted on the measure `ds`. It is irrelevant: only normalised
  series are compared, and (Fub) is an identity of `N`-series with no prefactor to track. The note says so in s.3.5.)

### 3.2 Critical point, Hessians, determinants (s.3.2): all CONFIRMED by hand and VERIFIED exactly (`rf_hess.py`)

* `p(L) = (L, 5L, zeta^j - L)`. The `sigma`-equation exponentiated is `-q (5L)^5 = prod_j (zeta^j - L) = 1 - L^5`, i.e. `u L^5 = 1`.
* `H_7`: `H_ss = 0`, `H_sy = 1/L`, `H_(s x_j) = H_(x_j x_j) = -1/(zeta^j - L)`, `H_yy = -1/(5L)`, all others 0. Read off the Taylor data
  that I built from the definition. The `y` and `x_j` linear terms are computed and vanish. The `s`-linear term
  `lambda + 5 log 5L - Sum_j log(zeta^j - L)` is the `sigma`-equation, which is checked in its exponentiated form above.
* `det H_7 = 1/(L^2 (1 - L^5)^2)`.
  * Order (a): block `diag(-1/(5L), -1/(zeta^j - L))`, `det = 1/(5L(1 - L^5))`. Schur complement
    `5/L + Sum_j 1/(zeta^j - L) = 5/(L(1 - L^5)) = S_0''(L)`, using `Sum_j 1/(zeta^j - L) = 5L^4/(1 - L^5)`.
  * Order (b): `H_(sigma,x0) = [[0, -1/(1-L)], [-1/(1-L), -1/(1-L)]]`, `det = -1/(1 - L)^2`, inverse `[[1-L, -(1-L)], [-(1-L), 0]]`. Hence
    the Schur complement is `H_ww - (1 - L) h h^T`.
  * I recomputed the four kinds of chart second derivatives at `p_ch(L)` (`yy`, `y x_j`, `x_j x_j`, `x_i x_j`): `(4L - 5)/(5L^2)`, `(1-L)/(L(zeta^j - L))`,
    `(2L - 1 - zeta^j)/(zeta^j - L)^2` and `(L-1)/((zeta^i - L)(zeta^j - L))`. They equal that Schur complement entry by entry.
* `det H_ch = -1/(L^2 Lam1(L)^2)`, by the matrix determinant lemma. In log coordinates `det(diag(-5L, -(zeta^j - L)) + (L-1) v v^T)
  = -5L Lam1 (1 + (L-1)(-5/L + Lam1'/Lam1)) = -25`, with `(L-1) Lam1 = L^5 - 1` and `(L-1) Lam1' = 5L^4 - Lam1`.
* VERIFIED exactly in `Q(zeta)` at eight `L`, and as polynomial identities by evaluation at 40 points (more than the degree bounds
  11 and 5): the scaled `det H_7` is the constant 5 and the log-coordinate determinant is -25. At `L = 1` the chart Hessian is
  `diag(-1/5, -1/(zeta^j - 1))`, with determinant `-1/25`. Controls fail: `+q`; a wrong closed form.
* Also used (my derivation, VERIFIED as `C H = -1` in `rf_chart.py`): `C_ch = -H_ch^-1 = diag(5L, zeta^j - L) + (L(L^5-1)/5) u u^T`
  with `u = (-5,1,1,1,1)`, and `C_7 = diag(0, 5L, zeta^j - L) + (L(L^5-1)/5) u u^T` with `u = (-1,-5,1,1,1,1,1)`.

### 3.3 Order (a) (s.3.3): CONFIRMED

* Inner critical point `(5 sigma, zeta^j - sigma)`, `phi = S_0` with consistent branches.
* (L0) + (L3) apply with `a = 5 sigma, M = 1` (amplitude `y^0`, measure `dy`) and `a = zeta^j - sigma, M = 0` (amplitude `x_j^-1`).
  This gives `B = z^-1 prod_j (zeta^j - sigma)^-1 exp(Sum_k z^(2k-1) S_2k)`.
* `det H_vv(sigma) = 1/(5 sigma prod_j (zeta^j - sigma))` and `beta B = (5L(1 - L^5))^(-1/2) z^-1 e^(S_1) exp(Sum_k z^(2k-1) S_2k)`, so
  `Nhat_7 = c(z; L)`. As a cross-check, the leading coefficients of the two orders agree: `z^-1/(1 - L^5)` from both.
* (L1b) needs `H_vv` invertible at `p` and at nearby `sigma`. It is, for every `L` in `Lam` (`5L != 0`, `zeta^j - L != 0`, `j = 0..4`).
* VERIFIED: the 7-variable SP equals `c(z; L)` exactly over `Q(zeta)` through `z^2` (`L = 1/2`) and mod two new primes through `z^3`
  (`L = -3`). The expectation algorithm is independent of the worker's (s.7).

### 3.4 Order (b) and the chart (s.3.4): CONFIRMED

* `Phi = sigma G + f` with `G = lambda + 5 log y - Sum_(j=0..4) log x_j` and `f = -y + Sum_j (-x_j + zeta^j log x_j)` (checked by
  expanding). The amplitude has no `sigma`, so (L2) applies: `G(p) = 0` and `d_(x0) G = -1/(1-L) != 0` for `L != 1`.
* `X_0(w) = -q y^5/prod_(j>=1) x_j` and `sigma_c = -(d_x0 f)/(d_x0 G) = 1 - X_0` (with `d_x0 f = -1 + 1/x_0`).
* `beta = (det H_vv(w)/det H_vv(w_*))^(-1/2) = X_0(w)/X_0(w_*)` on the branch with value 1 at `w_*`. So
  `beta N_inner = z^-1 (1 - L)^-1 prod_(j>=1) x_j^-1`, which is "the Jacobian `X_0` cancels `1/x_0`".
* `f(X_0(w), w) = lambda + psi(w)`, from `log X_0 = lambda + 5 log y - Sum_(j>=1) log x_j` and `-X_0 = q y^5/prod x_j`. The chart phase
  `psi` and its critical point `p_ch(L)`, with `K = L - 1`, are CONFIRMED by hand and in code (`K` computed from `q`).
* VERIFIED: chart SP `= c(z; L)` exactly over `Q(zeta)` through `z^3` (`L = 1/2`). Mod two new primes it holds through `z^5` at
  `L = -3, 5/3`, and the chart SP is `h` through `z^5` at `L = 1`. Controls fail as they should:
  * without the Jacobian, at `z^1`;
  * with `dy/y`, at `z^2`;
  * with `q -> -q`, `p_ch(L)` is no longer critical.

### 3.5 Fubini near L = 1, and the continuation (s.3.5): CONFIRMED, no gap

* **Fubini near `L = 1`.** (Fub) is applied at each FIXED `L` in `Lam`. It is an identity between formal series whose coefficients
  are polynomials in finitely many Taylor coefficients at `p(L)` and in the entries of inverse Hessians. No integral, estimate,
  radius of convergence or uniformity in `L` enters.
  * As `L -> 1`, the domain of the branch of `log x_0` shrinks (`x_0 = 1 - L`). So do the neighbourhood on which `X_0(w)` exists and
    the inner Hessian `-1/X_0^2`. This is invisible to an identity of Taylor coefficients at one point.
  * The intermediate quantities of order (b) that are singular at `L = 1`, namely `beta`, `det H_(sigma,x0)` and the leading
    coefficient `z^-1/(1 - L^5)` of `N_7`, enter only through `N_7`. They cancel in `Nhat_7`. `Nch` contains no factor `(1 - L)^-1`.
* **Continuation.**
  * `c_n` lies in `Q[L, L^-1, (1-L^5)^-1]` (s.1.6) and `Nch_n` lies in `Q(zeta)[L, L^-1, (zeta^j - L)^-1 : j >= 1]`. The latter needs the
    Taylor data of `psi`, the normalised amplitude, `adj H_ch` and `1/det H_ch = -L^2 Lam1^2` (a polynomial).
  * They agree on the infinite set `Lam`, so they are equal as rational functions. Hence `c_n` is regular at 1, and
    `c_n(1) = Nch_n(1)`, which is `N_n` of the `q = 0` data (the formula (F1) evaluated at `L = 1`, or (L1c)(i)).
  * VERIFIED independently of any fit: the Laurent expansion of `c_n` at `L = 1 + eps` has no negative powers and constant term `h_n`
    for all `n <= 40`. Control: `b_1 + 1/7` in the `j = 0` Stirling term leaves a pole `eps^-1` in `c_1`.

### 3.6 Step (d) and the sign of the twist weight (s.3.6): CONFIRMED

* At `q = 0`, `psi_0` is a sum of one-variable phases. (L0) + (L3), with `M = 1` for `y` (measure `dy`) and `M = 0` for `x_j`
  (`dx_j/x_j`), give `exp(Sum_k b_k z^(2k-1) T_k)`, `T_k = 5^(1-2k) + Sum_(j>=1) (zeta^j - 1)^(1-2k) = -Sum_(a in w_0) a^(1-2k)`
  (odd exponent).
* The sign: `Gamma(1 + 5/z)` contributes `+5^(1-2k) = -(-5)^(1-2k)`, i.e. the weight `-5` of `w_0` (the I-function numerator
  `prod (5 + kz)` = twist by `O(5)`).
* Check: `h_1 = b_1 (1/5 - 2) = -3/20 = ph_1` (forced by HPF, no free parameter). The weight `+5` would give `-11/60`.
* Givental's convention agrees (pp. 14-15: `ln R^00 = Sum_k (B_2k/2k) z^(2k-1)/(2k-1) Sum_(alpha != 0) (lambda_alpha - lambda_0)^(1-2k)`,
  with cotangent characters). With `lambda_alpha = zeta^alpha` this is the tangent part of `h_(w_0)` with the same sign.
* VERIFIED: `T_k` rational and `= -Sum_(w_0) a^(1-2k)` in `Q(zeta)`, and `h` from the five Stirling factors equals the Newton-sum `h`,
  through `z^40`. Controls fail at `z^1`: twist `+5`; the `Gamma(1+5s)` Stirling term dropped.

### 3.7 Step (e) (s.3.7; not used for Theorems Gamma, D): CONFIRMED

* The prefactor `[det(-H_ch(q))/det(-H_ch(0))]^(-1/2) A(p_ch)/A(p_0) = (L Lam1/5)(5/Lam1) = L` is right.
* The `q^d` ratio `(1/d!) prod_(k<=5d) (5/z + k)/prod_(j>=1) prod_(i<=d) ((1 - zeta^j)/z + i)` is right (`M = 5d+1` and `M = -d` in (L3)).
  The `j = 0` factor of the I-function gives the `d! z^d`.
* `D_q Dpsi = K = L - 1` by the envelope theorem.
* VERIFIED (`rf_eexp.py 6 12`): all 112 non-zero coefficients `q^d z^n`, `d <= 6`, `-d <= n <= 12`, exactly over Q. Here `c_n` comes
  from the HATTED polynomials of `rf_hat` (not from a fit of the sigma-form) and `I_d` from
  `prod_(j>=1) (1 - zeta^j + iz) = ((1+iz)^5 - 1)/(iz)`. Controls fail: `L` dropped, `e^(Dpsi/z)` dropped, twist `+5`.
