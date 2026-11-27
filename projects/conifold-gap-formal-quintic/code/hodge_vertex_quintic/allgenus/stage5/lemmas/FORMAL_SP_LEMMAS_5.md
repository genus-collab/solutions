# Stage 5 / lemmas, part 5: checks, caveats, record

Continuation of `FORMAL_SP_LEMMAS.md`, `_2.md`, `_3.md`, `_4.md`. Same labels.

## 11. Exact checks (VERIFIED; own code, every check with failing controls)

All scripts are new, stdlib only, import no project module, and were run as `py -I -B <script> ARGS` with the working directory in
the session scratch folder; each writes `logs/<name>_<args>.log` next to itself (relative paths only). Exact arithmetic over Q
(`fractions.Fraction`) throughout; random test data with small integer coefficients and fixed seeds.

### 11.1 Files
| file | content |
|---|---|
| `sp_core.py` | graded series with weight truncation (doubled weight: 1 per integration variable, 2 per parameter, 2 per power of z); Wick moments by the recursion (W), also with series-valued covariance; `SP^norm` (Definition 1.2) |
| `sp_iter.py` | total SP; inner critical point by the formal IFT; inner Taylor data, inner covariance `-K(y)^-1`, inner series `b(y; z)`; iterated SP (right side of L1b) |
| `sp_rand.py`, `chk_util.py` | random data; logging, controls, first differing order |
| `chk_fubini.py` | L1b on random polynomial phases; controls |
| `chk_change.py` | L1a for random polynomial changes of variables; FL1 for random exact forms; controls |
| `chk_delta.py` | L2 with and without outer variables; controls |
| `chk_param.py` | L1c(iii) (T-q), L1c(iv) (T-z), Lemma T(b) (translation by `-z`), FL4, FL3 for the SP with moving critical points; controls |
| `chk_stirling.py` | L3 by Wick's theorem in one variable against `St_m` built from Bernoulli numbers; (iv); the a-shift identity of s.9.3; homogeneity; controls |
| `chk_res.py` | FL2 and FL3 for the formal residue of s.10; controls |

### 11.2 Runs and results (all logs end with ALL CHECKS PASSED)
| log | what is compared | range | controls (each fails, as it must) |
|---|---|---|---|
| `chk_fubini_8_3.log` (84 s) | total SP = iterated SP, shapes `(nx, ny) = (1,1), (1,2), (2,1)` | `z^0..z^8`, 3 random phases each (degree <= 5, resp. 4) | no factor `(det K(y)/det K(0))^(-1/2)`; inner covariance frozen at `y = 0`; inner expansion at `x = 0` with the linear term dropped: each first differs at `z^1` |
| same log | shape `(2, 2)` (4 variables) | `z^0..z^6`, 3 phases | same three controls, `z^1` |
| `chk_fubini_8_2`, `_6_1`, `_4_1` | the same at smaller sizes | | |
| `chk_change_8_2.log` (14 s) | L1a: `SP_x[f; a] = SP_y[f o phi; (a o phi) det Dphi/det J]`, random cubic maps | `n = 1, 2`: `z^8`; `n = 3`: `z^6`; 2 instances each | without `det Dphi`; phase pulled back by `J` only: `z^1` |
| same log | FL1: SP of `Sum_i (d_i h_i + z^-1 (d_i f) h_i)` vanishes | same ranges | without the `z^-1 (d_i f) h_i` term: non-zero from `z^2` (`n = 1`), `z^0` (`n = 2, 3`) |
| `chk_delta_8_2.log` (26 s) | L2(3): SP over `(sigma, t)` = `a(0; z)` exactly | `z^8`, 2 instances | amplitude `a + s(1 + t)` (depends on sigma): `z^1` |
| same log | L2(4): SP over `(sigma, t, w)` = `SP_w[f(tau(w), w); a(tau(w), w) g_t/G_t(tau(w), w)]` | 1 outer variable: `z^8`; 2: `z^6`; 2 instances | without the Jacobian factor: `z^1`; sigma-dependent amplitude: `z^1` |
| `chk_param_8.log` (3 s) | T-q (moving critical point vs perturbation), FL4 (times z), T-z (`eps -> z`), T(b) (translation by `-z`), `n = 1, 2`; FL3 for the SP (`Y = X(1 + u h(X))`) | all `q^d z^k` (`u^d z^k`) with `d + k <= 8` | no determinant ratio; no `(d_q F)/z` term; no `e^g`; no exponential factor in T(b); no `Y'(X)`: all differ |
| `chk_stirling_40_12.log` (2 s), `chk_stirling_80_12.log` (11 s) | L3: `SP^norm_a[F; (x/a)^e] = St_(e+1)(z/a)` for `e = -3..4`, i.e. `m = -2..5` (`dx/x`) and `m = -3..4` (`dx`); (iv) for `m = -2..4`; the a-shift identity | `a = 1`: `z^40`, resp. `z^80`; `a = 3/7, -2`: `z^12` (homogeneity) | `b_2 -> -b_2`: `z^3`; `dx/x` versus `dx` at `m = 1`: `z^1`; `St_0` used for `m = 2`: `z^1` |
| `chk_res_6_14.log` | FL2 (3 random `Q`, `V = 1 + O(u)` random); FL3 for the residue (random CFU-type map, random `beta`) | `u^0..u^6`, `w^0..w^-14` | residue of `e^(-f/w) dX` is non-zero (`-1/w` at `u^0`); pull-back without `Y'(X)` differs |

The task asked for: L1b in 2 and 3 variables to `z^8` (done, plus 4 variables to `z^6`); L1a for a random analytic change of variables
(random polynomial maps, `n = 1, 2, 3`); L2 with random `f` (done, with up to 2 outer variables); L3 through `z^40` for `m = 0..4` and
both measures (done for `m = -3..5`, through `z^80` at `a = 1`); FL1-FL4 on examples (done; FL4 with moving critical points).
`St_0 = 1 + w/12 + w^2/288 - 139 w^3/51840 - 571 w^4/2488320 + ..` (the classical Stirling coefficients) is printed in the L3 logs.

### 11.3 A failure of my own, found and fixed
* The first run of the FL4 check (`chk_param.py 5`, before the fix) reported FAIL for `n = 1, 2`; T-q, T-z and FL3 passed.
* Diagnosis (scratch scripts, not kept): the perturbative form of FL4 passed, and T-q with the amplitude `z^-1 (d_q F) a + d_q a` failed in
  the harness, while T-q with `a` passed. The amplitude has terms of NEGATIVE weight (`z^-1 x`), and `sp_core.sp_norm` had pre-truncated
  the phase data `F_3/z` at the target weight. A term of weight `W + 1/2` of `F_3/z` times an amplitude term of weight `-1/2` lands at
  weight `W`, so the harness dropped needed terms. The same applies to a series-valued covariance.
* Fix: `sp_norm` now keeps phase terms up to weight `W - (minimal amplitude weight)`, and refuses (assertion) a negative-weight amplitude
  with a series-valued covariance. The FL4 check is run in the equivalent form multiplied by `z` (amplitude `(d_q F) a + z d_q a`, weight
  `>= 0`). All logs on record were produced after the fix; the FL1 check (amplitude weight `-1/2`, constant covariance, polynomial phase
  with all phase terms kept) was not affected.
* This was a harness error, not a counterexample: FL4 is a corollary of T-q (s.4.5), which passed throughout.

### 11.4 What the checks do not test
* Coefficients: all random data are over Q. The proofs hold over any Q-algebra; complex coefficients (`zeta`, numerical `L`) are not
  exercised by these scripts (the critic's `cr4_chart.py` works over F_p with `zeta in F_p`).
* The checks test the identities on polynomial phases of bounded degree; the proofs cover arbitrary formal power series (by Lemma U each
  coefficient involves only finitely many Taylor coefficients, so polynomial phases are no restriction for a given order).
* The analytic meaning (Laplace expansions) is not tested, and not needed.
* s.10.4 (the residue in the original variables) is not tested here; the poly referee's `rf_res` tests it numerically.

## 12. Scope, caveats, self-audit

* **AI-written, unrefereed.** No human has read the proofs. Every PROVED item has its argument written out in s.1-s.10; the only quoted
  statements are DLMF s.2.1(iii), (2.3.13)-(2.3.16) and 5.11.1, used only in Proof B of L3 (Proof A does not use them).
* **The analytic route was not used.** The polynomial-identity principle appears once, in (G5), for an EXACT Gaussian integral; no
  asymptotic theorem and no uniformity estimate is needed for L1a-L1c, L2, FL1-FL4. I did not read Hormander; nothing here cites it.
* **What "formal SP" means.** The definition needs the critical point and the inverse Hessian in the coefficient ring (e.g. `L` in a
  field `K` with `L != 0, 1, zeta^j`), and uses only Taylor series at that point. "Analytic" data enter through their Taylor series;
  convergence is never used. A statement about the actual integral is a different statement (s.1.4, Meaning).
* **Square-root branches.** Only normalised series are canonical. Each identity of full SPs is stated with an explicit branch rule
  (L1a, L1b, L2(4)); the project normalises by `c_0 = 1`, so constant factors never matter, but z-dependent symbols do; they are explicit.
* **Hypotheses the users must check.** Units: the Hessian at the critical point (L1a, L1c), the inner block and the Schur complement
  (L1b), `G_t` (L2), `a` (L3). For the chart argument these are listed in s.8, with two warnings (the 7-variable SP is not defined at
  `L = 1`; step (e) needs negative `m`).
* **SKETCH items.** s.10.4 (the relabelling between the original and the hatted residue in Theorem P(b) step 6) is a SKETCH; it is
  bookkeeping outside FL2-FL3, supported by the poly referee's numerical check through `u^30`.
* **Possible referee questions, answered in the text.** (1) Is `SP^norm` independent of choices? It depends only on the Taylor data at
  `c` (Definition 1.2); changing coordinates is L1a, moving terms between phase and amplitude is L1c(iv). (2) Why may `eps` be set to `z`?
  Because of the weight bound (L1c(iv)(2)). (3) Why is the Fubini proof free of uniformity issues? s.6, Remark. (4) Why is L2 exact?
  `C_gg = 0` (s.7). (5) Is `Gamma_St(y+1) = y Gamma_St(y)` proved or assumed? Proved, formally (s.9.2).
* **Not done here.** The chart argument (`stage5/chart/`); Theorem P(b) itself; Theorem U; (QRR). These lemmas make the "standard formal
  lemmas" of R5 (POLY s.6, poly referee C6) and of the chart outline explicit; whether they are applied correctly is for those proofs.
* **Not covered.** gamma's (R3) (GAMMA_PROOFS s.1: the ORIGINAL-variable residue is invariant under `x -> x/(1 - zx)`, i.e.
  `sigma -> sigma - z`). That change of variables is not small in `u`, so FL3 (s.10.3) does not apply. Nor is it small in `z` alone:
  each derivative of `e^(-G(x)/z)` brings a factor `1/z` that the factor `z` of `eps = z x^2/(1 - zx)` only cancels, so the pull-back
  series converges at best x-adically at fixed `u`-order. A proof needs its own ring; it is not attempted here. (R3) is used by gamma's
  Theorem A (the Picard-Fuchs equation of the residue), not by Theorem P(b).

### 12.1 Compliance record
* New files only, all in `allgenus/stage5/lemmas/` (5 notes, 10 scripts, 15 logs); no existing file was modified, moved or deleted.
* Code: Python 3.14 stdlib, `py -I -B`, working directory in the session scratch folder; no `__pycache__`; no project module imported;
  at most 3 of my processes at a time, all finished. Scratch-only files (two debugging scripts, run transcripts) stayed in the scratch
  folder. No wolframscript, no Mathematica, no cluster jobs.
* Web: three DLMF pages (s.2.1, s.2.3, s.5.11), read to check statements and equation numbers; no personal data in any request.
