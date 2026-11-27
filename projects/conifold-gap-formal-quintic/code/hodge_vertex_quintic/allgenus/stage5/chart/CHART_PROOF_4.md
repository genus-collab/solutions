# Stage 5 / chart, part 4: branches, what is not proved, checks, files (s.4-s.7)

Continuation of `CHART_PROOF.md`, `CHART_PROOF_2.md`, `CHART_PROOF_3.md`.

## 4. Every choice of branch, and why the normalised series does not depend on it

* **(B1) `lambda = log(-q)`.** It is forced by the critical point: given the branches (B2), `lambda := -5 log(5L) + Sum_j log(zeta^j - L)`.
  Another logarithm `lambda + 2 pi i k`, `k != 0`, gives `d_sigma Phi = 2 pi i k + O(|x - p(L)|)` near `p(L)`: no critical point there
  at all. The same holds in the sigma-form (s.1.6) and, for colour `i`, at `sigma = zeta^i L`.
* **(B2) `log y`, `log x_j` near `p(L)`; `log(5 sigma)`, `log(zeta^j - sigma)` on `Delta`.** Changing them by `2 pi i n_y`, `2 pi i n_j`
  changes `d_sigma Phi` by `2 pi i (5 n_y - Sum_j n_j)`, which (B1) compensates in `lambda`; then `Phi` changes by the constant
  `2 pi i Sum_j zeta^j n_j`. Only the critical value moves; `N` does not see constants (F2).
* **(B3) `log x_0` at `X_0 = 1 - L`.** In s.3.4 the identity `log X_0(w) = lambda + 5 log y - Sum_(j>=1) log x_j` is the equation
  `G(X_0(w), w) = 0` itself, with the branches fixed in `Phi`; it holds identically near `w_* = p_ch(L)`.
* **(B4) Square roots in Gaussian factors.** `det(-H)^(-1/2)` appears only in SP symbols, never in `N`. In (Fub) only
  `beta(w) = (det H_vv(w)/det H_vv(w_0))^(-1/2)` enters, with the branch `beta(w_0) = 1` (the ratio is close to 1 near `w_0`). In order
  (a) `beta = [5 sigma prod_j (zeta^j - sigma)/(5L prod_j (zeta^j - L))]^(1/2)`; in order (b) `beta = X_0(w)/X_0(w_*)`, no root at all.
  In the order-(b) symbol `(-X_0^-2)^(-1/2) = kappa X_0` with `kappa = +-i`: a constant.
* **(B5) `e^(S_1)`** contains `(5 sigma)^(1/2) prod_j (zeta^j - sigma)^(-1/2)`; `c` is normalised by `A_St(L, 0) = e^(S_1(L))`, and the
  normalised amplitude `e^(S_1(sigma) - S_1(L))` is the branch equal to 1 at `L` (a binomial series): no choice survives.
* **(B6) (L3).** `a^(M - 1/2)` and `log a` occur only in the symbol; the series `N = a^(M-1) Poch_M exp(...)` has an integer power of `a`.
* **(B7) Stirling forms with argument `a/z`.** `log(a/z) := log a - log z`; the `log z` parts cancel identically in `sigma` (s.1.6,
  `ch_zpow.py`), and the `log a` are the branches of (B2).
* **(B8) GAMMA versus POLY form** (s.1.7): the odd integers `eps_j`, `eps` and the even `e` are fixed by the critical point (`e = 0`).
* **(B9) The colour.** `L = u^(-1/5)` would require a fifth root; here `L` is the independent variable and `q(L) = (1 - L^-5)/5^5` is
  single-valued, so no root is taken. Colour 0 near `q = 0` means `L` near 1; `sigma = zeta^i L` is colour `i` (s.3.9).
* **(B10) Step (e).** `L(q) = (1 - 5^5 q)^(-1/5)` is the formal power series with `L(0) = 1`, `log L` the one with value 0; `Dpsi` is a
  power series in `q` with rational coefficients (no constant logarithm survives in the difference).

So `c(z; L)`, `Nhat_7(z; L)` and `Nch(z; L)` depend on `L` only.

## 5. What I could not prove, or did not check (labelled honestly)

1. **Stage-4 inputs of Theorem D** (R0, R2, R3, R5) are quoted, not re-proved. They are AI-written and AI-refereed; Theorem U is
   AI-audited only; GI Theorem G is AI-refereed. No human has checked any of them. Theorem D is PROVED only given them.
2. **(QRR)** is an input. arXiv:1702.06096 eq. (5) (the definition of the formal quintic) was not read by me. arXiv:math/0108100
   Thm 9.1, which I read, is the UNTWISTED equivariant statement `R|_(Q=0) = exp diag(b_i)` (Newton sums over the cotangent characters);
   the twisted version with the weights `w_i` (including `-5 zeta^i`) is the project's convention (FORMAL_QUINTIC_genus2 s.0), not a
   statement I have located. Its `z^(5k)` part is exactly what (N) tests. Theorem Gamma itself does not use (QRR): it computes `c(z; 1)`.
3. **(FL), residue part.** "The formal residue kills exact forms and is invariant under formal changes of variables, termwise in `u`"
   and Stirling's functional equation `Gamma_St(y+1) = y Gamma_St(y)` in `Q[[1/y]]` are used by Theorem P but are not written out
   here (elementary; the gamma referee VERIFIED the functional equation through `v^60`). The SP part of (FL) is (S1), (L1a), (L1c)(iii').
4. **(L3)** is PROVED modulo two classical analytic theorems, DLMF 5.11.1 (Stirling's series) and DLMF 2.3.8 (Watson's lemma), whose
   statements I read on the DLMF site (2026-11-14). A purely formal proof is not given (s.2.8). (L3) is VERIFIED exactly to `z^40`.
5. **The lemma proofs** (S1), (PG), (L0)-(L2) are mine: written out, not refereed, and not compared with the parallel `stage5/lemmas`
   agent. Remark (F4) cites Hoermander's Theorem 7.7.5 from memory; it is used in no proof.
6. **Literature.** Read: DLMF 5.11.1, 5.11.3-5.11.4, 2.3(ii)-(iii); arXiv:math/0108100 v2 s.9-10 (pp. 11-15; s.3.8). Not read by me:
   arXiv:0710.0889 eq. (5) (for the reading of step (e) as (Gamma3)); arXiv:1702.06096 and arXiv:1803.01409 (item 2).
7. **Not claimed:** any statement about actual integrals or analytic asymptotics in `q`; polynomiality `c_n in Q[L]` (the chart gives
   `Q[L, L^-1]`; polynomiality is not needed; VERIFIED for `n <= 16` in `ch_eexp.py`, PROVED elsewhere per SYN4 s.1.1).
8. Every modular computation is evidence at the stated primes only.

## 6. Checks (own exact code, with failing controls)

All scripts are new, stdlib only, run as `py -I -B script.py ARGS` with the working directory in the session scratch folder; they
import no project module and no Stage-4 module. Shared modules: `ch_core.py` (exact `Q(zeta_5)`, `F_p`, `Q`, `Q[[w]]/(w^D)`; Bernoulli
numbers; `h` by Newton sums), `ch_wick.py` (graded polynomials; the functional `N` of s.1.2 by the Isserlis recursion; a one-variable
version), `ch_phases.py` (Taylor data of the sigma-form, of `Omega_7` at `p(L)` and of `Omega_ch` at `p_ch(L)`, with the control
variants), `ch_hess_sym.py` (polynomials over `Q(zeta)[L]`, Bareiss). Primes: `p1 = 4611686018427388081`,
`p2 = 1208925819614629174706261` (both `= 1 mod 10`; Miller-Rabin with the bases 2..41 is a proof below `3.3e24`); they are NOT the
Stage-4 critic's primes. Logs in `logs/`, relative paths only. Every line of every log is PASS; every FAIL is a deliberate control.

| script and arguments (log) | what is checked | field, range | time |
|---|---|---|---|
| `ch_zpow.py` | s.3.1, s.1.6: exponent of `z` of the six Euler substitutions `= -1`; inner Stirling factors of order (a) give `z^(+1)`; coefficient of `log z` = 2; `q = 0` chart and every `q^d` term of (e) carry `z^(+1)`. Controls: factor `j = 0` dropped, `Gamma(1 + 4s)` | symbolic, `Q(zeta)[sigma, 1/z]` | < 1 s |
| `ch_hess.py` (+ `ch_hess_sym.py`) | s.3.2: critical equation with `q = (1 - L^-5)/5^5`; gradients; `K = L - 1`; Hessian entries; `det H_7` directly and by both orders; Schur complement of order (b) = chart Hessian entrywise; closed forms; log-coordinate form; `q = 0` chart Hessian `diag(-1/5, -1/(zeta^j - 1))`, det `-1/25`. Symbolic: `det = 5`, `det = -25`, `prod (zeta^j - L) = 1 - L^5`, `Sum_j prod_(i != j) = 5L^4`. Controls: `q -> -q`, `+q`, wrong closed form | exact `Q(zeta)`, 8 rational `L`; `Q(zeta)[L]` | < 1 s |
| `ch_stirling.py 30` | s.3.6: `T_k = -Sum_(w_0) a^(1-2k)`, `k <= 30`; `exp(Sum b_k T_k z^(2k-1)) = h` through `z^59`; product of the five normalised Stirling series = `h` through `z^25`; `h_1`, `h_5`. Controls (fail at `z^1`): twist `+5`, `1 - zeta^j`, `Gamma(1 + 5/z)` dropped, `b_k -> -b_k` | exact | < 1 s |
| `ch_euler.py 40 40` (`ch_euler_40.log`) | (L3): `x^m`, `m = -3..3`, measures `dx/x` and `dx`, `a = 1, -2/7, 3/2`, through `z^40`; DLMF 5.11.4 coefficients. Controls: `B_2k/(2k)` (`z^3`), Pochhammer dropped (`z^1`), `1/Gamma` (`z^1`) | exact Q | 16 s |
| `ch_sigma.py 10` | s.1.6-1.7, s.3.9: two assemblies of the sigma-data agree; `c_n` rational; `c_1 = -(3/20)L^4`, `c_2 = (9/800)L^8`; stored `c_5(2) = 11495493197/31500`; `c_n(zeta^i L) = zeta^(-in) c_n(L)`, `i = 1..4`, `n <= 10`. Control: `zeta^(+in)` | exact `Q(zeta)` | 3 s |
| `ch_chart.py 3 4 6 5 7` | s.3.3-3.6: chart SP at `L = 1` = `h` (exact `z^7`; mod p1, p2 `z^6`); chart SP = `c(z; L)` (exact `z^3` at `L = 2, -2/3`; mod p1, p2 `z^6` at `L = 3, 5/4, -1/2`); 7-variable SP of `Omega_7` = `c(z; L)` (exact `z^4` at `L = 2`; mod p1, p2 `z^5` at `L = 3`). Controls: no Jacobian (`z^1`), `dy/y` (`z^2`), `Omega_7` with `dy/y` against `c` (fails) and against the `Gamma(5s)` sigma-form (passes) | exact `Q(zeta)`; `F_p1`, `F_p2` | 2137 s |
| `ch_lemmas.py 6 4` | s.2: (S1); (L1a); (L1b) with inner dimension 1 and 2 (inner series over `Q[[w]]/(w^15)`); (L2); four seeds, through `z^6`; (L1c)(iii) with the moving critical point, `q^0..q^4`, `z^-d..z^6`, and against (L3). Controls: no Jacobian, no `beta`, `sigma` in the amplitude, `sigma^2` in the phase, no critical-value exponential, no Hessian ratio | exact Q | 44 s |
| `ch_eexp.py 6 10` | s.3.7: `L e^(Dpsi/z) c(z; L(q)) = h(z) Sum_d q^d I_d(z)`, all 98 non-zero coefficients, `d <= 6`, `-d <= n <= 10`; `D_q Dpsi = L - 1`; `c_n` polynomial in `L` with exponents `4n - 5k` (`n <= 16`, fits confirmed on 3 further points). Controls: `L` dropped, `e^(Dpsi/z)` dropped, twist `+5` in `h` | exact Q | 21 s |

Smaller runs with the same verdicts (logs kept): `ch_euler.py 8 6`, `ch_chart.py 1 1 2 2 2` and `3 2 4 3 6`, `ch_lemmas.py 3 1` and `5 3`,
`ch_eexp.py 2 3` and `4 6`, `ch_sigma.py 8`.

**What the checks add to the proof.** Nothing logically: Theorem Gamma is proved in s.3. They test the bookkeeping that a written proof
can get wrong: the measures (`dy` versus `dy/y`), the Jacobian of (L2), the sign of `q`, the twist weight, the normalisations, and the
lemmas themselves on random data. New relative to Stage 4: the 7-variable integral itself (both orders against the total), exactness
over `Q(zeta)` of the chart step, the chart against `c(z; L)` directly (the critic compared with `R_n(L)`), order `z^6`, and the
coefficientwise agreement with the I-function (step (e)).

**Development record (honest).** Bugs found and fixed before any reported run: (1) in `ch_lemmas.py` the `-y` term was missing from the
linear Taylor coefficient of the (L1c) test; the critical-point assertion stopped the run. (2) `ch_eexp.py` first evaluated the
sigma-data over Q (which has no `zeta`) and crashed; a control then raised a KeyError (missing negative powers of `z`); both fixed.
(3) Two scripts had leftover lines in a draft (`ch_stirling.py`, `ch_hess.py`), removed before the first run. No check other than a
deliberate control ever failed. The `ch_euler_40.log` is from `40 40`, which overwrote the earlier `40 20` run (same verdicts).

## 7. Files and compliance

* Notes: `CHART_PROOF.md` (s.0-s.1), `CHART_PROOF_2.md` (s.2), `CHART_PROOF_3.md` (s.3), `CHART_PROOF_4.md` (s.4-s.7).
* Code: `ch_core.py`, `ch_wick.py`, `ch_phases.py`, `ch_hess_sym.py` (modules); `ch_zpow.py`, `ch_hess.py`, `ch_stirling.py`,
  `ch_euler.py`, `ch_sigma.py`, `ch_chart.py`, `ch_lemmas.py`, `ch_eexp.py` (scripts); logs in `logs/`.
* Every file is new and under 250 lines; no existing file anywhere was modified, moved or deleted; nothing was written outside
  `stage5/chart/` except throw-away output in the session scratch folder. At most three of my processes ran at a time; none is left
  running. No `__pycache__` (`-B`). Logs and notes contain no absolute path and no user name.
* Web: DLMF s.5.11 and s.2.3 (nist.gov), the arXiv abstract page of math/0108100 and its PDF (pages 11-15 read); no personal data in
  any request.
