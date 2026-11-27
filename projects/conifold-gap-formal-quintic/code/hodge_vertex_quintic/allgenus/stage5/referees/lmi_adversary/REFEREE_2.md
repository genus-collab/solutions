# Stage 5 / referees / lmi_adversary, part 2: computations, comparison, corrections, limits, record (2026-11-14)

Continuation of `REFEREE.md` (s.0 verdict, s.1 the proof re-derived, s.2 the attacks, s.3 the falsifier map). Same labels
and abbreviations. AI-written; unrefereed.

## 4. My computations (new code; exact over Q, or modulo two new primes)

### 4.1 Method (`lm_core.py`, `lm_wick.py`; an algorithm, PROVED: it is the definition of CP s.1.2)
* **No rescaling.** At `L = s* + e`, every Taylor datum of the one-variable integrand is an element of `LS`, the Laurent series in
  `e` with ABSOLUTE precision: all coefficients of `e^k` with `k < prec` are exact. The singular phase coefficients
  `alpha (-1)^k/(k(k-1)) e^(1-k)` and the singular amplitude (`(1 + t/e)^E`, `(z/(alpha (e + t)))^n`) are exact Laurent polynomials.
  Polynomial regular data are exact. Truncated Taylor data (the quintic) carry their true precision. `C = -1/(2 p_2)` is a
  series, since `D = alpha0 + e g''(L)` is a unit. So "no negative power of `e`" is an exact statement once the final precision
  is `>= 1`; this is checked in every run. No a priori pole bound and no scaling lemma is used.
* **Wick sum.** `E = exp(Sum_k p_k s^(k-2) tau^k) Sum_(m,j) a_(m,j) s^(2m+j) tau^j` and `N_n = Sum_r (2r-1)!! C^r [s^(2n) tau^(2r)] E`.
  The vertex part comes from the exponential recursion in the weight (its terms are `tau^(w+2V)`), the amplitude from a graded
  exponential in `(z, t)`.
* **Two-variable chart** (`lm_chart.py`). At the critical point belonging to `L` one has `Q e^(mu_c) = alpha0 e` exactly and
  `mu_c = -g'(L)/alpha0`, so all chart data are polynomials in `e`. The covariance is `C_ss = -e/D`, `C_sm = -1/D`,
  `C_mm = 2 G_2/(alpha0 D)`. Phase and amplitude are separable; the moments `<eta_s^a eta_m^b>` are counted by cross-pairings.
  At `e = 0` this is (L2') itself.
* **Fields.** Q, or `F_p` with `p1 = 4611686018427387761` (the largest prime below `2^62` that is `1 mod 20`) and
  `p2 = 1180591620717411303161` (the same below `2^70`). Both are found by a deterministic Miller-Rabin test, and both carry
  `sqrt(-1)` and a primitive fifth root of unity. They are not RA's primes (`999999999999999631`, `999999999999999999661`) and
  not Stage 4's (`2^61 - 1`, `2^89 - 1`).
* **Independence.** Written from the definitions before I opened RA's or GEN's code. It imports no project or agent module and
  shares only the definitions (`c`, `h`, `b_k`, `beta_n`). Three different methods now exist: RA (scaled Wick sum, truncated in
  `eps`), GEN (Lagrange-Buermann in `T`) and this one (unscaled Laurent-series Wick sum at `L`, with precision tracking).

### 4.2 Engine checks, LMI predictions, proof steps (`lm_run_eps.py`; logs `lm_run_eps_10_Q`, `_12_p1`, `_14_p2`)

| id | check | result |
|---|---|---|
| E1 | quintic at `L = 2`: `c_1..c_5 = -12/5, 72/25, 354278/625, -73040834/3125, 11495493197/31500` (the stored values) | PASS |
| E2 | quintic `c_1(1+e) = -(3/20)(1+e)^4`, `c_2 = (9/800)(1+e)^8` (28 and 27 `e`-coefficients); `h_1..h_5` | PASS |
| E3 | (L3), `z^max(N,10)`: `M = -3..4` against both forms; `M = 1/3, -7/5, 5/2, 2/7` against the Bernoulli form; generic `M` in `F_p` (2 each); controls: Bernoulli part of `M = 0` for `M = 1/3`, and `Poch` with `(1 - z)`, both fail | PASS |
| E4 | one-Gamma identity, `g = 0`, `F = 1`: `c_n = 0` identically, 5 cases (`M = 1/3` included) | PASS |
| P | LMI on random data (deg `g` = 6; `F = exp(f)`, `f` of degree 3 in `z`, 5 in `x`): 12 cases `(alpha0, M0)` = `(-1,0)`, `(-1,1)`, `(-1,-2)`, `(2,0)`, `(3,1)`, `(-1/3,2)`, `(5,1)`, `(1/2,-1)`, `(-7/4,3)`, `(-1,1/3)`, `(2/3,-2/5)`, `(-3,7/2)`; plus `g''(s*) = 0`. Mod p, 4 more: `alpha0 = sqrt(-1)` with `M0 = 0, 2`, and random `alpha0` with generic `M0` | PASS |
| S1 | steps 1+2: the 1-variable `c_n(s* + e)` equals the 2-variable chart `N_n(e)` coefficientwise (up to 20 `e`-coefficients per `n` at `z^14`); 4 cases, `M0 = 1/3` included | PASS |
| S2 | step 3: the chart has no negative power of `e`, and its value at `e = 0` is the (L2') formula | PASS |
| S3 | control: a pure `(2/7) m^3` term breaks (L2') | fails as predicted |
| S4 | control: singular `b_1 + 1/7` breaks step 1 (1-variable != chart) and gives poles of order `n` | fails as predicted |

Ranges: exact `z^10` over Q (45 checks, 152 s); `z^12` mod `p1` (51 checks, 35 s); `z^14` mod `p2` (51 checks, 70 s).

### 4.3 Falsifiers (`lm_run_x.py 12`; exact over Q, `n <= 12`; 37 checks, 492 s)
* X1, `Gamma(a0/z)^2`: pole of order exactly `n` at every `z^n`. Leading coefficients at `alpha0 = -1`: `1/16, 9/512, 75/8192,
  3675/524288, ..`, i.e. GEN's `((2n-1)!!)^2/(n! 16^n)`; at `alpha0 = 2` they carry the predicted factor `(-1/alpha0)^n` (my scaling).
* X2 (`Gamma(a0/z) Gamma(-a0/z)`): pole orders `0, 2, 4, .., 24`. X3 (slopes 1, 2) and X3b (`Gamma(a0/z) Gamma(1 + a0/z)`): orders `n`.
* X4 (`alpha1 = 1/3`, `ell` in `{0, 1, 5/2}`): regular at fixed `ell`. For `M0 = 0`, `[e^1] c_n` depends on `ell` for every
  `n = 1..12`, the `[e^0]` values do not, and they equal `F(s*, z)/F(s*, 0)`. For `M0 = 1`, `[e^0] c_n` depends on `ell` for every `n`.
* X5: `d c_1(s*)/d ell = 2/3` exactly. X5+: regular, value of `F_0`, for `M0 = 0, 1, -1`.
* X6, `alpha0 = -1, 2`: first poles `(n0, order, coefficient)` = `(1, 1, -1/7)`, `(2, 2, 1/7)`, `(3, 3, 2/3)`; `(1, 1, 1/14)`, `(2, 2, 1/28)`,
  `(3, 3, -1/12)`, all `eta/alpha0^n0`. Reason: the perturbation adds exactly `eta alpha0^(-n0) e^(-n0)` to `c_n0`, because its
  `t`-dependence enters only at higher orders of `z`. X6m: `-1/9`, `1/18` (`= 1/(9 alpha0)`).
* X7: `-1/8, 3/8, 1/24` (`alpha0 = -1`) and `1/16, -3/16, -1/48` (`alpha0 = 2`), all `(1/24 + beta_1(M) - E^2/2)/alpha0`, derived from
  the five order-`z` Wick terms. X7+: exponent `+1/2` with `M = 0` is regular, with the value of `M0 = 1`.
* X8: leading `1, 2, 6, 24, 120, 720` (`alpha0 = -1`) and `-1/2, 1/2, -3/4, 3/2, -15/4, 45/4` (`alpha0 = 2`), i.e. `n! (-1/alpha0)^n`.
* X9+: regular, with the value of `M0 + 1`. X9 (`(alpha0, M0) = (-1, 1), (3, 2), (2, -1)`): only simple poles, and
  `[e^-1] c_n = (-M0/alpha0) V_(n-1)` for all `n <= 12`.
* X10: with `M0 = 1`, the `delta = 0` value and the wrong-sign value both differ from the series.

### 4.4 The quintic (`lm_run_q.py`; `z^12` exact, 9 checks, 156 s; `z^16` mod `p1`, `p2`, 13 checks each, about 60 s)
* Q1: no pole at `L = 1`, and `c_n(1) = h_n`. The quintic's data enter as Taylor series at `sigma = 1` of one rational function
  `-(4s^3 + 3s^2 + 2s + 1)/(1 + s + .. + s^4) = Sum_(j>=1) 1/(zeta^j - s)`, truncated at `x^(3N+12)` with tracked precision.
* Q2: with `Gamma(5s)` the value is still `h`; `c_n` differs from Q1 exactly for `n >= 2`.
* Q3: `Gamma(1 + (1 - sigma)/z)` is regular, its value is the LMI value with `delta = z`, and it differs from `h` at `z^1`.
* Q4: at `L = 0`, no pole, and `c_n(0) = [z^n] Ginf` (the Bernoulli closed form) = the LMI value with `delta = -z/5`.
* Q5: at `L = 0` with `Gamma(5s)`, the value `exp(5 Sum_(5|2k-1) b_k z^(2k-1))` is not `Ginf`.
* Q6 (mod p, explicit `zeta`): colours 1 and 2 are regular, with value `h(zeta^-i z)`; the control (value `h(z)`) fails.

## 5. Comparison with GEN s.3.3 (+ GEN2 s.5) and with RA2 s.4.4

* **GEN's statement is RA's,** provided `Lam(0)` is read as the constant term of `Lam(z)` (its `z -> 0` limit
  `-Phi'(s*) e_*/c_* + A_0(s*)`). The linear part of `Phi` at `s*` contributes only that constant. Read literally, `Lam(0)` is undefined.
* **GEN's SKETCH** (a formal Legendre-Laplace dual, which needs `Phi''(s*) != 0`, plus Zariski density) is superseded by RA's chart,
  which needs neither; my P case with `g''(s*) = 0` holds. GEN's label "CONJECTURAL" (SYN5 O2) is obsolete.
* **The shift rule** `h_f = e_f - c_f e_*/c_*` (GEN s.3.1(d), GEN2 s.5.1 D) is a corollary of LMI and Lemma F5.11.8. A regular
  factor `Gamma_St(e_f + Y_f)` at `sigma = s* + delta` is `Gamma_St(e_f - c_f e_*/c_* + w_f/z)`, re-expanded by F5.11.8. GEN's T2, T3
  and GEN2 5.1 A-E are thus instances of LMI (one degenerate factor with a linear argument, `q`-free regular data); I did not
  re-check GEN's Lemma 2 (the summands).
* **GEN's G, H1, H2** are my X1, X3 and X6m; I reproduce GEN's `K_0` numbers exactly. GEN's F, F2, F3 (random regular data, all
  `(c_*, e_*)`) are my P cases, extended to complex `alpha0` and generic `M0`.
* **RA2 s.4.4 inside LMI's hypotheses:** X1 (singular Stirling dropped), X2, X7, X8 (singular `b_k` perturbed), X3 (all six `b_1`)
  and X5 (exponent 0) violate "A is the exact Stirling form", and fail at step 1. X4 (the square) violates "g holomorphic" and
  fails at step 2. RA's pattern "a perturbation of `b_k` first shows at `z^(2k-1)` with top order" is my X6 formula, which
  includes the coefficient `eta/alpha0^(2k-1)`; RA's X5 is my X7 with `E = 0`, coefficient `1/(8 alpha0)`. No conflict anywhere.

## 6. Corrections (presentation and bookkeeping; none affects validity)

* **C1** (RA s.3.3 step 1): introduce `Omega_2` as a formal integrand. The Euler-integral motivation carries `z^(M0 + a0/z)`, whose
  `sigma`-dependent part, if kept, would shift `lambda` by `alpha0 log z` (REFEREE s.1.2).
* **C2** (RA s.3.1, "PROVED from ..."): the list of inputs should also contain (F1)-(F2), the translation case of (L1a), the
  holomorphic IFT (CP s.1.4) and the `C[[z]]`-linearity of `N` (used for the `nu`-free factor `F`).
* **C3** (RA s.3.1): "for `L != s*` near `s*`": the punctured disc on which `g` and the `F_m` are holomorphic and
  `alpha0 + (L - s*) g''(L) != 0` (equivalently `S''(L) != 0`, `det H_2 != 0`).
* **C4** (RA s.3.4, (Gamma)): "the `j >= 1` part of `e^(S_1)`" should read "`e^(S_1)` without its `j = 0` term `a0^(-1/2)`, i.e.
  `(5 sigma)^(1/2) prod_(j>=1) (zeta^j - sigma)^(-1/2)`".
* **C5** (RA s.3.4, `Ginf`): DLMF 5.11.8 via uniqueness of asymptotic expansions is correct. Lemma F5.11.8 (REFEREE s.2.3) makes
  the step formal.
* **C6** (RA s.0 item 4 and s.3.4 "What this says"): "ANY holomorphic phase and amplitude" should say: holomorphic at `s*`,
  independent of `lambda` (holomorphic dependence on `Q` is allowed, s.2.5), with `a0` exactly linear and `F(s*, 0) != 0`. The
  last condition is a real hypothesis, not only a normalisation (X9). `M0` may be any complex number with the
  Bernoulli-polynomial form (REFEREE s.2.3).
* **C7** (SYN5 s.0, s.2 S2, K6): LMI and (L2') are "PROVED (u)" there. Now they are PROVED and refereed once (this report, AI).
  Lemma F5.11.8, (L3-C) and the `M0`-in-C corollary are new here and unrefereed.
* **For the write-up:** (Gamma) is blind to the measure of the twist factor (REFEREE s.2.4). Say where the correct `dy` is used
  (Theorem P(b), the hatted/flat identity), so that a reader does not take (Gamma) for a test of it.

## 7. What I did not check; limits

* I did not check CP3, CP4 (the 7-variable chart), the referee_formal reports, the inputs of Theorem D, (QRR), P(b) or (CYC).
* For the lemmas LMI quotes, I checked their statements and their use, not their proofs (refereed by RF and RL). FSL4's Proof A
  was read, not re-refereed; (L3-C) and F5.11.8 rest on its (i) and (iv).
* All computations are finite: exact `z^10` (random LMI data, steps), `z^12` (quintic, falsifiers); modulo two primes `z^12`, `z^14`
  (random data), `z^16` (quintic and colours). The random data are rationals of height `<= 7` (or uniform in `F_p`); the tests are
  generic, not exhaustive.
* Some falsifier facts were observed, not predicted: the pole orders of X2, X3, X3b; the X8 coefficient `n! (-1/alpha0)^n`, explained
  afterwards by the inverse Mellin transform `E_1` of `Gamma(w)/w` and DLMF 6.12.1 (an explanation, not a proof); and the X4 limit.

## 8. Record

**Files** (new, in this folder; code stdlib only; no project or agent module imported; each file under 250 lines):

| file | content |
|---|---|
| `REFEREE.md`, `REFEREE_2.md` | this report |
| `lm_core.py` | field Q or `F_p`; `LS` (Laurent series in `e` with absolute precision); `exp` of `LS` |
| `lm_wick.py` | Bernoulli numbers and polynomials, `beta_n`, graded exponential, vertex table, the one-variable Wick engine |
| `lm_data.py` | LMI integrands at `L = s* + e` (singular factors with overrides, `g`, `F`, extra powers and terms), the LMI value |
| `lm_chart.py` | the two-variable chart of LMI step 2 and its Wick sum |
| `lm_fals.py` | falsifier builders: non-linear singular argument, `lambda`-dependent and `Q`-dependent amplitudes |
| `lm_models.py` | the quintic in LMI form at `s* = 1`, `zeta^i`, `0`; the fixed point `L = 2`; `h`, `Ginf`; random data |
| `lm_harness.py` | logging, checks, comparisons, primes, roots of unity |
| `lm_run_eps.py`, `lm_run_x.py`, `lm_run_q.py` | sections E, P, S; X; Q |
| `logs/*.log` | the 7 final runs below; `logs/DRAFT1_lm_run_x_4.log` (first draft, 2 FAIL lines, item 2 below) |

**Final runs** (`py -I -B`, working directory in the session scratch folder, at most three of my processes at a time):
`lm_run_eps.py 10 Q` (45, 152 s), `12 p1` (51, 35 s), `14 p2` (51, 70 s); `lm_run_x.py 12` (37, 492 s); `lm_run_q.py 12 Q`
(9, 156 s), `16 p1` (13, 63 s), `16 p2` (13, 59 s). That is 219 checks, ALL CHECKS PASSED in each log; every failure-type
outcome is a labelled falsifier or control, reported as the PASS of the predicted failure. The logs print relative paths only.

**My own errors, found and fixed (reported plainly).**
1. *A wrong prediction (X8).* For `F = 1/(sigma - s*)` I first predicted the leading pole coefficient `(-1/alpha0)^n`, from a
   wrong expansion of `1/w` in inverse falling factorials. The smoke run showed `n!`. The recorded check uses the observed formula
   and is labelled "observed".
2. *A wrong sub-claim (X4).* My first X4 check also required the `ell`-dependence of `[e^1] c_n` to be affine in `ell`. It is not, so
   two FAIL lines appear in `DRAFT1_lm_run_x_4.log` (kept). The main claim, `ell`-dependence, held there. I removed the sub-claim
   and made the draft's observation, an `ell`-free limit for `M0 = 0` and `ell`-dependence at `e^0` for `M0 = 1`, into checks.
3. *A hang.* My comparison helper looped up to `2^40` on two exact series. Fixed; the stuck process, mine, was stopped.
4. *A reused prime.* My first `p2` (the largest prime below `10^21` that is `1 mod 20`) turned out to be RA's `PR2`. I replaced it
   by a prime below `2^70` before the final runs.
5. *File size.* `lm_core.py` was first written as one piece of about 306 lines, against the 250-line rule. I split it into
   `lm_core.py` and `lm_wick.py` (pure refactor) and re-ran every recorded log from the split code. The pre-split development logs
   (all PASS, apart from item 2) stayed in the scratch folder.
6. *A stray file.* One shell command began with a stray `cat > /tmp/...`, which blocked on stdin and created an empty file in
   `/tmp`. I stopped my process and deleted that empty file. Nothing else was written outside the scratch folder and this folder.

**Compliance.** Only new files, only in this folder. No existing project file was modified, moved or deleted. No `__pycache__`.
No process I did not start was signalled. No Mathematica, no cluster. Web: DLMF s.24.2, s.24.4 and s.6.12 (public pages), with
no personal data in any request. Logs and notes contain no absolute path and no user name.
