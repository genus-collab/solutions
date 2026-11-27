# Stage 5, task `general`, part 2: the hypotheses one at a time, T5, the verdict (2026-11-14)

Continuation of `GENERAL.md` (s.0 bottom line, s.1 the eps-method, s.2 summands, s.3 the general chart argument, s.4 T1-T4).
Same labels; AI-written, unrefereed.

## 5. Testing the argument: its hypotheses one at a time (`g5_t6.py`, `g5_t6b.py`, `g5_t6c.py`; exact over Q, `n <= 20`)

Each line: what the chart argument (s.3.1-s.3.3) predicts, then what happened. Logs `logs/g5_t6_20.log`, `logs/g5_t6c_20.log`.

### 5.1 The argument applied at OTHER points (it must work there too)
* **A. Formal quintic at `L -> 0` (`q = oo`).** Solve `G = 0` for the twist variable `y` (`c_* = 5`, `e_* = 1`, `s* = 0`). Then
  `y_0 = (prod_j x_j/(-q))^(1/5)` and the measure `dy` inserts `y_0`, i.e. `x_j^(1/5)`: the chart at `q = oo` is
  `prod_j Gamma(zeta^j/z + 1/5)`. Prediction: `c(z; 0) = prod_j exp(Sum_n beta_n(1/5) (z/zeta^j)^n)`
  `= exp(Sum_l (-1)^(l+1) B_(5l+1)(1/5) z^(5l)/(l(5l+1))) = Ginf(z)`, the closed form of SYN4 s.1.1 (CONJECTURAL there).
  **Observed:** no pole at `L = 0` and `c_n(0) = [z^n] Ginf` for all `n <= 20` (exact). The UNshifted tail
  `prod_j exp(Sum b_k (z/zeta^j)^(2k-1))` fails at `z^5`: the shift 1/5 is needed and is the one the measure gives.
  With (Gamma) and `c = R` (SYN4 s.1.3) this is the statement `f_m = [z^(5m)] Ginf` of Stage 4, now with a mechanism.
* **B. Twisted P^4, general weights, at `q = oo`** (all five critical points tend to 0 when `prod lambda_j != 0`). Prediction
  `exp(Sum_n beta_n(1/5) z^n Sum_j lambda_j^-n)`. **Observed:** exact agreement, no pole, `n <= 20`, W1, W2, W3.
* **C. Local P^2 at `q = oo`** (three critical points tend to `mu/3`; singular factor `Gamma(3s - mu/z)`, `e_* = 0`). Prediction:
  NO shift, value `exp(Sum_k b_k z^(2k-1) Sum_j (lambda_j - mu/3)^(1-2k))`. **Observed:** agreement, `mu = 0` and `1/2`, `n <= 20`.
* **D. A shifted singular factor** `Gamma(1/3 + (lambda_i - sigma)/z)` in the twisted P^4 (not a GW summand; a test of the shift
  rule `h_f = e_f - c_f e_*/c_*`, here `e_f + c_f/3`). **Observed:** no pole, value with the shifts, colours 0 and 2; the unshifted
  value fails at `z^1`.
* **E. Orbifold-type singular factor**, weighted P(1,2): `Gamma(3/z - s) Gamma(-1/z - 2s)`; at `s* = -1/2` the singular factor has
  `c_* = -2` (the root `x_*^2 = ...` in the chart). Prediction: no pole, unshifted regular tail (`e_* = 0`). **Observed:** agreement at
  both colour points (`c_* = -1` and `-2`), `n <= 20`.

### 5.2 The regular factors need not be Gamma factors (the argument is local)
* **F.** Regular part = random WKB data: `Phi''`, `A_0'`, and `A_n` for EVERY `n` (even `n` included, absent for `e_f in {0, 1}`),
  random rationals of small height; singular factor `Gamma_St((s* - sigma)/z)`. Prediction (LMI, s.3.3): no pole, value
  `exp(Sum_n z^n A_n(s*))`. **Observed:** three trials, `n <= 20`, exact agreement.
* **F2.** The same with `(c_*, e_*) = (-1, 1/3), (-1, 1), (2, 1/3), (5, 1)`: value = normalised `B` at the first pole
  `delta = -e_* z/c_*`. **Observed:** agreement in all four, `n <= 20`.
* **F3.** No regular phase (`Phi = 0`; the Legendre step of the SKETCH is unavailable, the Zariski-density step covers it).
  **Observed:** no pole and value `exp(Sum z^n A_n(s*))`, two trials, `n <= 20`.
* **F'.** Controls with the same random regular data: `b_1 + 1/7` in the singular factor gives a pole at `z^1`; an extra term
  `(1/7) Y^-2` in the singular Stirling series gives a pole at `z^2`.

### 5.3 Where a hypothesis fails (the argument says: no conclusion; the data say: it really fails)
* **G (C1 fails: coincident weights).** P^1 with `lambda = (1, 1)`, P^2 with `(1, 1, -2)`, twisted P^4 with `(1, 1, 3, -2, 5/2)`:
  the singular factor is `Gamma_St((s* - sigma)/z)^2`; P^2 with `(1, 1, 1)`: cube. **Observed:** poles of order EXACTLY `n` at every
  `z^n`, `1 <= n <= 20`, in all four. For the square the `eps^-n` coefficient of `c_n` is `((2n-1)!!)^2/(n! 16^n)` (`n <= 20`,
  exact), independent of the regular factors: Hankel's coefficients of `K_0(x)` at `x = 2 y_c`, as `Int Gamma(y)^2 X^-y dy` is a `K_0`
  (standard Mellin pair; used only as an explanation). For the cube: `1/9, 2/81, 14/2187, 8/19683, -440/177147, ...`.
* **H1 (C1 fails: two singular factors of different slope).** Twisted P^4 with `lambda_0 = 0` at `q = oo`: `Gamma(-sigma/z)` and
  `Gamma(1 + 5 sigma/z)` both degenerate at `sigma = 0` (four critical points, `sigma^4 ~ 1/q`). **Observed:** poles from `z^1` on.
* **H2 (C2 fails).** Elementary part of `Gamma_St(1/2 + Y)` with the Bernoulli series of `Gamma_St(Y)` (equivalently a stray
  `(sigma - s*)^(1/2)` in `B`): pole at `z^1`. Positive control: elementary part of `Gamma_St(1 + Y)` with the series of `Gamma_St(Y)`,
  which IS `Gamma_St(1 + Y)` (`beta_n(1) = beta_n(0)`): no pole.
* **C3 fails (a vanishing twist weight).** By Lemma 4 (PROVED) the colour does not exist; nothing to test.
* Also from T1 (GENERAL s.4.2): perturbing the singular factor's Bernoulli series gives poles; perturbing only the regular ones
  keeps regularity and changes the value exactly as LMI predicts.

**Reading.** The removable singularity is a property of ONE factor: in every case tried it holds when (C1)-(C3) hold and fails
when (C1) or (C2) fails; the value is then fixed by the regular part at the first pole of that factor. This is exactly what the chart argument (in the general form
of s.3.3) says, including the shift rule and the `q = oo` chart.

## 6. T5: the critic's Wick test of the chart step at `z^6`

`stage4/critic/cr4_chart.py` and `cr4_core.py` copied unchanged to the scratch folder (they write next to themselves), run there
as `py -I -B cr4_chart.py 6`; the log is copied here as `logs/t5_cr4_chart_6.log`. ALL CHECKS PASSED at `2^61-1` and `2^89-1`
(742 s, real 12 min 22 s):
* `L = 1`: the 5-variable chart SP factorises into Stirling series and equals `h_(w_0)` through `z^6`;
* `L = 2, 3/2, -2/3`: the chart SP equals `R_n(L) = L^(4n) ph_n(L^-5)` for `n <= 6`; the stored normalisation differs at `z^5` by
  exactly `m5` (so `z^5` sees (N));
* controls: no Jacobian `|x_0|` fails at `z^1`; measure `dy/y` fails from `z^2` on.
This tests steps (b)-(c) of SYN4 s.1.4 one order beyond Stage 4 (VERIFIED mod two primes; not exact).

## 7. Verdict: which statement is supported in which generality

| # | statement | generality | status |
|---|---|---|---|
| S1 | (Gamma): `c_n(1 + eps)` has no pole and `[eps^0] c_n = h_n` | formal quintic, colour 0 | VERIFIED exactly `n <= 60` (new method); proof SKETCH (SYN4 s.1.4), unrefereed |
| S2 | the analogue at every colour, value `h_(w_i)`, `w_i = {lambda_i - lambda_j} U {-5 lambda_i}` | O(5)-twisted P^4, distinct non-zero rational weights, `Sum lambda_j` arbitrary | VERIFIED exactly `n <= 24`, 4 weight sets x 5 colours; chart argument applies verbatim (Lemma 3, PROVED); CONJECTURAL |
| S3 | the analogue with the chart-predicted QRR weights (tangent weights; Euler twist: minus, inverse Euler twist: plus the fibre weight) | P^1, P^2; local P^2; P^1 with O(-1)+O(-1); P^4 with O(3), O(4), O(6); lifts included | VERIFIED exactly `n <= 20`, 12 theories, 40 colours; summands PROVED (Lemma 2); CONJECTURAL in general |
| S4 | LMI: one exact singular Gamma factor + ANY regular WKB part: no pole, value = regular part at the first pole | any `(c_*, e_*)`, any regular data | VERIFIED exactly `n <= 20` (random data, 9 cases); SKETCH (s.3.3); CONJECTURAL |
| S5 | `q = oo` chart: `c(z; 0) = Ginf(z)`; general weights `exp(Sum beta_n(1/5) z^n p_-n)` | formal quintic; twisted P^4 | VERIFIED exactly `n <= 20`; SKETCH (same chart, `x_* = y`) |
| S6 | (C1), (C2) cannot be dropped: poles in every tested case where one fails; no colour when a twist weight vanishes | coincident weights (2, 3), two slopes, mismatched factor, perturbed Bernoulli | VERIFIED in these cases (`n <= 20`), not in general; Lemma 4 PROVED |

Not supported or not tested:
* nothing here PROVES (Gamma) or LMI: all ranges are finite, and (L1)-(L4) are not written out;
* the (QRR) input itself: the tests compare with the QRR formula `h_(w_i)`; they do not compute the `q = 0` R-matrix of the
  theories independently. The sign rule is the one the chart predicts and the project uses (`-5 zeta^i`, pinned by (H) at `n = 1`);
* weights in a number field other than Q (the quintic's `zeta` enters only through rational power sums), irrational or
  symbolic weights, and non-hypergeometric I-functions;
* `n > 20` for S3-S6, `n > 24` for S2; `n > 60` exactly for S1 (the `eps`-method at `z^80` would take an estimated 1 h or more).

**Consequence for the programme.** The chart argument survives a falsification test of the ARGUMENT, not only of the identity:
its conclusion holds in every case where its hypotheses hold, it fails exactly where (C1) or (C2) fails, and its two
side-predictions (the shift rule, the `q = oo` chart giving Ginf) hold. Its general form S4 suggests a shorter write-up for
Stage 5: prove LMI once, from (L1)-(L4), in one Laplace variable plus one Euler variable, and obtain (Gamma) (and Ginf) as the
special cases `B` = the regular Gamma factors at `q = 0` (resp. at `q = oo`). This is a suggestion, not done here.

## 8. Self-criticism: limits, failures, bugs

* Everything above is a finite exact computation, except Lemmas 1-4 (PROVED, elementary). No statement about all orders is PROVED.
* Independence: the engine and the builders are new and import no project or agent module. They share with Stage 4 only the
  definitions of `c` and `h` and the Bernoulli numbers. `g5_valid` compares with my own fixed-point SP code and with ONE number
  printed in SYN4 (`c_5(2)`).
* The random tests (F, F2, F3: 9 cases) use rationals of small height; they test the statement generically, not every datum.
* The values are compared with the QRR FORMULA; an error common to the formula and to the chart prediction (e.g. a global sign
  convention of the twist) would not be seen. The project pins the sign independently ((H) at `n = 1`, Stage 4).
* For S2-S6 the comparison is with predictions I derived myself (s.3); the predictions were written into the code before the
  runs and none was adjusted afterwards. The only change after a run was the threshold of check V2 (next item).
* Bugs found during the work, none affecting a reported result:
  * the positive control `Gamma(5s)` first crashed (a keyword clash in `g5_data.quintic_groups`); fixed, then run;
  * `g5_valid` V2 first used a wrong tolerance formula and printed FAIL although the agreement was `10^-27`; I replaced the
    formula by a fixed threshold `10^-20`; the compared numbers did not change;
  * my first draft of `residue_check` used `Q = (-1)^N z^(Sum c) q`; re-deriving it before any run gave `(-1)^(N + l_conc)` for
    concave twists. No run used the draft; the sign of `Q` does not enter any `c_n` anyway (s.1).
* After the first runs I let the engine take several singular factors (for H1). The change is backward compatible; all record
  logs were produced again with the final code (s.9).
* Not done: T1 beyond `z^60`; a mod-p extension; a re-derivation of the QRR sign from the QRR literature; (L1)-(L4) written out.

## 9. Record

**Files** (all new, each under 250 lines; code stdlib only, imports no project or agent module):

| file | content |
|---|---|
| `GENERAL.md`, `GENERAL_2.md` | this note |
| `g5_core.py` | the `eps`-engine of s.1 (Lagrange-Buermann in `T`, Laurent coefficients in `eps`), exact Q or F_p |
| `g5_data.py` | Bernoulli numbers and polynomials, `beta_n(h)`; groups of regular factors; Newton power sums for the quintic; predictions |
| `g5_theory.py` | Gamma-factor descriptions of the theories (s.2), the residue check against the I-function, one generic test run |
| `g5_valid.py` | validation of the engine by an independent fixed-point SP (s.1) |
| `g5_t1.py`, `g5_t2.py`, `g5_t3.py` | tests T1 (+ quintic controls), T2 (+ controls), T3 (+ controls) |
| `g5_t6.py`, `g5_t6b.py`, `g5_t6c.py` | the hypothesis tests of s.5 (A-E; F, F2, G; F3, H1, H2, L0) |
| `logs/*.log` | the logs below (relative paths only); `logs/t5_cr4_chart_6.log` is the critic's log from the T5 run |

**Runs** (final code, `py -I -B`, working directory in the session scratch folder; times include up to three runs in parallel):

| command | log | result | time |
|---|---|---|---|
| `g5_valid.py 6 24` | `g5_valid_6_24.log` | ALL CHECKS PASSED | 1 s |
| `g5_t1.py 40 40 Q 12` | `g5_t1_40_40_Q.log` | ALL CHECKS PASSED | 105 s |
| `g5_t1.py 60 60 Q 6` | `g5_t1_60_60_Q.log` | ALL CHECKS PASSED | 750 s |
| `g5_t2.py 24 10` | `g5_t2_24.log` | ALL CHECKS PASSED | 330 s |
| `g5_t3.py 20` | `g5_t3_20.log` | ALL CHECKS PASSED | 258 s |
| `g5_t6.py 20` | `g5_t6_20.log` | ALL CHECKS PASSED | 116 s |
| `g5_t6c.py 20` | `g5_t6c_20.log` | ALL CHECKS PASSED | 19 s |
| critic `cr4_chart.py 6` (unchanged copy) | `t5_cr4_chart_6.log` | ALL CHECKS PASSED | 742 s |

Smaller smoke runs (`z^4` to `z^12`) stayed in the scratch folder. Every FAIL-type outcome in these logs is a control, reported as
the PASS of a check that the control fails in the predicted way.

**Compliance.** No existing file was modified, moved or deleted; everything written into the project is in `stage5/general/`.
The critic's two scripts were copied to the scratch folder and run there unchanged. At most three of my processes ran at a time
(one short-lived waiting loop that would have been a fourth was stopped at once); I signalled no process I did not start. Web
access: one public page, DLMF s.5.11 (statement 5.11.8); no personal data in any request. Logs print relative paths only; `-B`
prevented `__pycache__`.
