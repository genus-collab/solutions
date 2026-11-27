# Stage 5 / chart, adversarial referee, part 2: computations, comparison, limits, record (2026-11-14)

Continuation of `REFEREE.md` (s.0 verdict, s.1 hand check, s.2 the attacks, s.3 the two-variable chart and Theorem LMI).
Same labels and abbreviations. AI-written; unrefereed.

## 4. My computations (new code, new method; exact over Q or modulo two new primes)

### 4.1 Method (`ra_sp.py`; PROVED as an algorithm, it is CP s.1.2 after a change of variables)
* *Integrand.* A Mellin-Barnes-type summand at a colour point `s*`: one singular factor `Gamma_St(M0 + alpha0 (sigma - s*)/z)`
  (power `count`), groups of regular factors `Gamma_St(M + (alpha sigma + beta)/z)` (Stirling coefficients and the exponent `-1/2`
  can be overridden), an optional extra phase `g(sigma)` and an extra amplitude `F(sigma, z)`. All `log(a/z)` read as `log a`.
* *Scaling.* `L = s* + eps`, `sigma = L + eps tau`, `z = eps z'`. Every Taylor datum in `(tau, z')` is then a power series in `eps`.
  The singular factor becomes `eps`-free: its phase coefficients are `alpha0 (-1)^k/(k(k-1))`, its amplitude
  `(1 + tau)^(-1/2) exp(Sum_k b_k (z'/alpha0)^(2k-1) (1 + tau)^(1-2k))`. The quadratic coefficient is a unit of `F[[eps]]`.
* The Wick series in `z'` therefore has coefficients `cp_n(eps)` in `F[[eps]]`, and `c_n(s* + eps) = eps^(-n) cp_n(eps)`. So `c_n` is
  regular at `s*` iff `[eps^k] cp_n = 0` for all `k < n`, and then `c_n(s*) = [eps^n] cp_n`.
* *Computation.* The direct Wick sum of CP s.1.2 over `R = F[eps]/(eps^(N+1))`, `F = Q` or `F_p`. A ring `RNum` evaluates at one
  `eps` instead (then `c_n = eps^-n N_n`).
* *Prediction.* `predict()` evaluates the residue formula of Theorem LMI from the same factor data at `eps = 0`. Closed forms are
  computed separately where they exist.
* *Independence.* Written before I read `stage5/general`. It imports nothing from the worker, the general agent or Stage 4. The
  quintic's `zeta` enters by two routes that share nothing:
  * route `m`: the `e_j = zeta^j - 1` are the roots of `m(s) = ((1+s)^5 - 1)/s`, and
    `Sum_r t^(r-1) Sum_j (e_j - eps)^(-r) = -m'(eps + t)/m(eps + t)`, over any field;
  * route `zeta`: an explicit `zeta` in `F_p`.
* Primes (`ra_core.py`): `PR1 = 999999999999999631`, `PR2 = 999999999999999999661`. Both are `= 1 mod 10`, deterministic
  Miller-Rabin; they are neither the worker's nor Stage 4's.

### 4.2 Engine checks (log `logs/ra_run_18_24.log`)
* E1: `c_1..c_5(2) = -12/5, 72/25, 354278/625, -73040834/3125, 11495493197/31500`, the worker's `ch_sigma` and `ch_chart` values and
  the Stage-4 value `c_5(2)`. Also `c_1..c_3(-2/3) = -4/135, 8/18225, 4067314/110716875`, the worker's values.
* E2: the Laurent data at `L = 1` give `c_1 = -(3/20)(1+eps)^4` and `c_2 = (9/800)(1+eps)^8` in all available `eps`-orders.
* E3: routes `m` and `zeta` give identical Laurent data mod PR1 and PR2 through `z^24`.
* E4 (`ra_poly.py 8 6`, log `logs/ra_poly_8_6.log`). For each `n <= 8`, fit `c_n = Sum_k a_k L^(4n-5k)` to the Taylor data at
  `L = 1`. The fit is confirmed on 8 to 21 further coefficients, and its value at `L = 2` equals the independently evaluated
  `c_n(2)`. (`a_1 = 0` for `n = 2`: `deg ph_2 = 0`, as in SYN4.)
* C0: `h_1..h_5 = -3/20, 9/800, -269/80000, 2823/6400000, -257071/8064000000`.

### 4.3 Predictions: every one holds (`ra_run.py 18 24`, 594 s; exact `z^18` over Q unless stated)
"Regular" means no negative power of `eps` in any `c_n`, `n <= 18`. "Value" means `c_n(s*)` = the residue formula (and the closed
form, where one is given) for all `n <= 18`. "CP" means inside CP's 7-variable argument as written; "LMI" means it needs s.3.

| id | integrand (colour 0 unless stated) | predicted value at `s*` | covered by | result |
|---|---|---|---|---|
| C1, C2 | the quintic | `h` | CP | regular, value: exact `z^18`; mod PR1, PR2 `z^24` |
| C3 | the quintic at `sigma* = zeta^i`, `i = 1, 2` (explicit `zeta`) | `h(zeta^-i z)` | CP s.3.9 | mod PR1, PR2 `z^24` |
| P1 | `Gamma(M + 5s)`, `M = 0, 2, 3` | `Poch_M(z/5) h` | CP (measure `y^(M-1) dy`) | PASS (3) |
| P2 | `Gamma(M0 + (1 - sigma)/z)`, `M0 = 1, -1, 2` | residue at `sigma = 1 + M0 z`; closed forms `Poch_6(z/5) (5/m(z)) h` and `Poch_-4(z/5) h` for `M0 = +-1` | CP (chart amplitude `X_0^M0`) | PASS (3) |
| P3 | `Gamma(M + (zeta^j - sigma)/z)` for all `j >= 1`, `M = 1, 2, -1` | `h`, `(m(-z)/5) h`, `(5/m(z)) h` | CP | PASS (3) |
| P4 | `Gamma(1 + d s)`, `d = 0, 3, 7` (non-CY; `d = 0` is equivariant `P^4`) | `exp(Sum b_k z^(2k-1) [d^(1-2k) + Sum_j e_j^(1-2k)])` | CP | PASS (3) |
| P7 | rational weights `(1; 3, -1/2, 5/2, -2)`, `d = 5`; `(1; 2, -3)`, `d = 2` | product of the regular Stirling factors | CP | PASS (2) |
| P9 | `Gamma(2(1 - sigma)/z)` for the `j = 0` factor | `h` | CP (square-root chart) | PASS |
| P2b | singular exponent `+1/2`, i.e. `Gamma(1 + w) = w Gamma(w)` | as `M0 = 1` | CP | PASS |
| P10 | `1/Gamma((1 - sigma)/z)` | residue of `Gamma(1 + (sigma - 1)/z)` at `sigma = 1 - z` (reflection) | CP after reflection | PASS |
| P5 | amplitude times `e^(2 sigma)`, `1/(3 + sigma)`, `(2 + sigma)`, `e^(z sigma^2)` | `h`, `h`, `h`, `e^z h` | LMI | PASS (4) |
| P6 | phase plus `(3/7) x^2`, `(-2/5) x^3`, `x^2 - x^3 + 2x^4`, `x = sigma - 1` | `h` | LMI | PASS (3) |
| P8 | ONE singular factor, random `g`, `F`; `(alpha0, M0)` in `(-1,0), (-1,1), (-1,-1), (2,0), (2,1), (-1/3,0), (-1/3,2)` | LMI formula | LMI | PASS (7) |
| I1-I4 | a REGULAR factor perturbed: `y` with `b_1 + 1/7`; `x` with `b_2 + 3`; `x` Stirling dropped; `y` exponent `+1/2` | LMI formula with the perturbed data | LMI | PASS (4) |
| O1 | the quintic at `L = 0` (`q = oo`): singular factor `Gamma(1 + 5 sigma/z)`, `alpha0 = 5`, `M0 = 1` | `Ginf(z)` (Bernoulli polynomials; `f_1, f_2, f_3`) | LMI (s.3.4) | PASS |

The two-variable chart itself (`ra_2var.py 9`, 169 s, exact) was tested against the one-variable series. The 2-variable Wick
series (2x2 covariance, Isserlis moments) of the chart of s.3.3 equals the normalised 1-variable full-Stirling series. This holds
at `eps = 1/2` and `eps = -3/7`, and as an `eps`-series, through `z^9`. Cases: the quintic (`M0 = 0, 1`), the orbifold point and two
random core-lemma cases. Control: the 1-variable form without the singular Stirling series differs.

### 4.4 Falsifiers: the identity fails, and so does a hypothesis of the argument (exact, `z^18`)

| id | modification | failing step of CP's argument | observed |
|---|---|---|---|
| X1 | Stirling series of the singular factor dropped | s.3.3: order (a) produces the FULL Stirling series of `x_0` by (L3), so the modified summand is no reduction of any `Omega_7` | pole at `z^1`; pole of order `n` at every `z^n` |
| X2 | singular `b_1 + 1/7` | the same ((L3)) | pole at `z^1` |
| X3 | `b_1 + 1/7` in all six factors (Stage-4 control) | the same, through the singular factor | pole at `z^1` |
| X5 | singular exponent `0` for `-1/2` (no Gamma function) | the same | pole at `z^1` |
| X7 | singular `b_3 + 1/7` | the same | first pole at `z^5`, order 5 |
| X8 | singular `b_2 + 1/7` | the same | first pole at `z^3`, order 3 |
| X4 | `Gamma((1 - sigma)/z)^2` | `Omega_8` exists, and (a), (b) hold for `L != 1`. Step (c) fails: the second variable `x_0'` has coefficient `1 - 1 = 0` in front of `log x_0'`, and its critical point `x_0' = 1 - L` tends to 0. The chart data have poles at `L = 1` (two degenerate factors, against LMI) | pole at `z^1`; order `n` at `z^n` |

Pattern: a perturbation of the singular `b_k` first shows at `z^(2k-1)`, as a pole of the top order `2k - 1`. This is the leading
layer `[eps^0] cp_n`, which comes from the singular factor alone. No falsifier passed. No predicted case failed.

### 4.5 (L2') on random data (`ra_l2src.py 8 8`, exact, `z^8`, 8 seeds)
* The full 2-variable Wick series equals `exp(P(delta)/z) A(delta, z)` for random `P`, `A`, `alpha0`, and `c` in `{0, +-1, 2, -3, 1/5}`.
* Controls: the wrong sign of `delta` is rejected; a pure `m^3` term in the phase breaks the identity.
* Smaller run `ra_l2src.py 6 4`.

## 5. Comparison with `stage5/general` (read after my runs; GEN, GEN_2)
* *Independent and convergent.*
  * GEN found the same Laurent method in `eps` (T1). Its engine is Lagrange-Buermann, mine the direct Wick sum, and the codes are
    disjoint.
  * GEN found the same local statement (GEN s.3.3, "LMI").
  * GEN found the same `q = oo` prediction `c(z; 0) = Ginf` (GEN_2 s.5.1 A).
  * GEN found the same failure pattern: perturbing the singular factor gives poles; perturbing the regular factors changes only the
    value; coincident weights give poles of order exactly `n`.
  * There is no contradiction anywhere between GEN's results and mine.
* *What GEN has that I do not.*
  * Larger ranges for (Gamma) itself (exact `z^60`).
  * The theories T2-T3 (twisted `P^4` at all colours, local `P^2`, `P^1` with `O(-1)+O(-1)`, `O(k)` twists).
  * Non-integer `e_*`.
  * The Hankel `K_0` coefficients of the double pole.
  * A re-run of the critic's `cr4_chart.py 6`.
* *What this report adds.*
  * A PROOF of LMI (s.3.3). GEN labels LMI CONJECTURAL with a SKETCH that needs a formal Laplace dual of the regular part (it
    requires `Phi''(s*) != 0`) and a Zariski-density step; my chart needs neither.
  * The delta lemma with a source (L2'), PROVED and tested.
  * Ginf PROVED given Theorem D's inputs, and polynomiality `c_n in Q[L]` from the chart at `L = 0` (REFEREE s.3.4).
  * The identity-of-Gamma cases P2b, P10.
  * `sigma`-dependent amplitude and phase perturbations of the quintic itself (P5, P6).
  * Colours via an explicit `zeta` in `F_p` (C3).
  * The two-variable chart tested against the 1-variable series (s.4.3).
* *Alignment.* My s.3 covers integer `M0`. For non-integer `M0` (GEN's `e_*`), replace (L3) by its analogue with DLMF 5.11.8 (same
  proof); (L2') already allows any `c`.

## 6. What I did not check; limits
* The Stage-4 inputs of Theorem D: R0 (Theorem U), (R1) of UR, R2-R3 (Lemma N, the connection form), R5 (Theorem P with GI
  Theorem G, the CFU map, the Brieskorn reduction), the residue part of (FL) and Stirling's functional equation. Also (QRR): the
  chart computes `c(z; 1)` and never uses (QRR); that the actual `q = 0` R-matrix is `diag(h_(w_i))` is an input.
* CP s.3.8 (Givental): not read.
* The worker's multi-variable codes were not re-run. My tests are one- and two-variable. GEN re-ran the critic's 5-variable test.
* All my computations are finite. The ranges are:
  * exact `z^18` over Q (C, P, I, O, X);
  * `z^24` mod PR1, PR2 (C2, C3, E3);
  * the 2-variable chart through `z^9`;
  * (L2') through `z^8`.
* Theorem LMI and (L2') are my proofs: written out, not refereed.

## 7. Record

**Files** (all new, in this folder; each under 250 lines; stdlib only; import no project, worker or agent module):

| file | content |
|---|---|
| `REFEREE.md`, `REFEREE_2.md` | this report |
| `ra_core.py` | fields Q and F_p, primes, Bernoulli numbers, generalised binomials, the rings `RSym` (`F[eps]/eps^K`) and `RNum` |
| `ra_sp.py` | integrand specification (`Sing`, `Group`, extra `g`, `F`), Taylor data in `(tau, z')`, Wick sum, residue prediction |
| `ra_models.py` | the quintic and its modifications, colours `zeta^i`, rational weights, the core lemma, the orbifold point; power-sum routes |
| `ra_run.py` | sections E, C, P, I, X, O (s.4.2-s.4.4) |
| `ra_2var.py` | the two-variable chart of s.3.3 against the 1-variable series |
| `ra_l2src.py` | (L2') on random data |
| `ra_poly.py` | E4: the polynomial shape of `c_n` read off its Laurent data at `L = 1` |
| `logs/*.log` | all runs below (relative paths only; every line PASS; "identity FAILS as expected" lines are PASSed falsifier checks) |

**Runs** (`py -I -B`, working directory in the session scratch folder; at most three of my processes at a time):

| command | result | time |
|---|---|---|
| `ra_run.py 4 5`, `5 6`, `8 10`, `12 16` | ALL CHECKS PASSED | 1 s, 1 s, 9 s, 68 s |
| `ra_run.py 18 24` | ALL CHECKS PASSED (56 checks) | 594 s |
| `ra_2var.py 3`, `6`, `9` | ALL CHECKS PASSED | 1 s, 22 s, 169 s |
| `ra_l2src.py 6 4`, `8 8` | ALL CHECKS PASSED | under 2 s |
| `ra_poly.py 8 6` | ALL CHECKS PASSED | 1 s |

**My own errors, found and fixed (reported plainly).**
1. *Harness bug.* The first draft of `ra_run.py` compared the numeric `L = -2/3` values without dividing by `eps^n`, and printed
   FAIL for E1. The numbers themselves were right: `4/81 = (-5/3)(-4/135)`.
2. *Two misclassified falsifiers.* "Singular exponent `+1/2`" and "`1/Gamma`" were run as falsifiers and showed no pole. They are
   `Gamma(1 + w) = w Gamma(w)` and the reflection formula, i.e. covered cases. They are now positive tests P2b and P10, each
   compared with its own prediction.
3. *Blind check.* At `z^4` the `b_3` falsifier is blind (`b_3` enters at `z^5`); the check is now skipped below `z^5`.
4. *Bad test point.* `ra_2var.py` first used `eps = 1` at the orbifold point, which is `L = 1`, a singular point there. It now
   uses `eps = 1/2, -3/7`.
5. Cosmetics: an unreachable branch of the source term in `ra_2var.py` and the title line of `ra_l2src.py`.

The first-draft logs (`ra_run_4_5.log` with the two FAIL lines of items 1 and 3, `ra_2var_3.log` before item 4,
`ra_l2src_6_4.log` before item 5) were overwritten by re-runs of the final code. That log of the first run of `ra_l2src.py 6 4`
was deleted by me before its re-run. No other file was deleted or changed, and nothing outside this folder was written except
scratch output.

**Compliance.** Only new files, only in this folder. No existing project file was modified, moved or deleted. No `__pycache__`.
No process I did not start was signalled. Web: DLMF s.5.11 and s.2.3 only, with no personal data in any request. Logs and notes
contain no absolute path and no user name.
