# Stage 5 / referees / lmi_formal, part 2: computations, limits, record (2026-11-14)

Continuation of `REFEREE.md` (s.0 verdict and the list C1-C6, s.1 statement, s.2 (L2'), s.3 the proof, s.4 the
consequences). Same labels and abbreviations. AI-written; unrefereed.

## 5. My computations (new code, exact over Q; one float check)

### 5.1 Method and independence
* *Independence.* All code is new, Python stdlib only, and imports only its own modules `lf_*.py`. I did not read or import
  RA's `ra_*.py`, the main session's `independent_check/lmi_check.py`, or any worker, Stage-4 or project module.
* *The functional `N`* (`lf_wick.py`) is CP s.1.2 implemented from its definition: with `xi = t eta`, `t = z^(1/2)`, the
  normalised integrand is `exp(E)`, `E = Sum_(k>=3) ph_k t^(k-2) eta^k + Sum logA_(d,m) t^(d+2m) eta^d`; `exp` by the recursion
  `r f_r = Sum_k k E_k f_(r-k)`; then Isserlis moments. Two-variable engines: `wick_2sep` (data whose cubic-and-higher phase and
  amplitude separate, `G(s) + K(m)`, `A_s A_m`, which is the case for every chart of LMI; FULL 2x2 covariance) and `wick_2gen`
  (arbitrary data). The moment table is checked against the Isserlis recursion (L4).
* *No `zeta` in the quintic data.* The Taylor data at rational `L` are rational: `P_r(x) = Sum_(j=0..4) (zeta^j - x)^-r` from
  `5 y^4/(1 - y^5)`, `p_r(x)` (`j = 1..4`) from `-Phi5'(y)/Phi5(y)`, `y = x + t`. `h = h_(w_0)` is computed from its DEFINITION in
  the field Q(zeta_5) (class `Cyc5`, inverse by the norm), a route that shares nothing with the power sums (Q0 links them).
* *Three ways to see regularity at `s*`:* (i) the chart = `c^F` at `L != s*`, with the chart coefficients holomorphic by
  inspection (s.3.2); (ii) Method A: `sigma = L + eps tau`, `z = eps w` ((L1a), (F3)), over `Q[eps]/(eps^(N+1))`, so
  `c^F_n(s* + eps) = c'_n(eps)/eps^n`, regular iff `[eps^k] c'_n = 0` for `k < n`, value `[eps^n] c'_n`; (iii) for the quintic, a
  GLOBAL fit: `c_n(L)` exactly at 32 rational points, fitted by a Laurent polynomial with exponents `4n - 5k` including two
  slots below 0 and two above `4n`, confirmed on three further points (no expansion at any special point).

### 5.2 Checks (every line of every log is PASS)

| id | what | range, field | script (log) |
|---|---|---|---|
| Q0 | `P_r = p_r + (1 - x)^-r`; `P_r(0) = 5 [5 \| r]`; `p_r(x)` = explicit sums in Q(zeta_5) at `x = 1, 2` | `r <= 40`, `r <= 25`; exact | `lf_quintic 14 60` |
| Q1 | `h_1..h_5 = -3/20, 9/800, -269/80000, 2823/6400000, -257071/8064000000` (= RA_2 C0); `h_1`, `h_2` as quoted in CP | exact | same |
| Q2 | chart at `s* = 1` = `c(z; L)` at `L = 2, -2/3, 1/2, 3, -1, 5/4`; `det H = -g''(L) a0(L) - 1 = -a0(L) S_0''(L)`; `c_1..c_5(2)` = RA_2 E1 | `z^14`, exact | same |
| Q3 | chart at `L = 1` (`Q = 0`) = `h` from its definition; `det H = -1` | `z^60` | same |
| Q4 | chart at `s* = 0` = `c(z; L)` at `L = 1/2, -1/3, 2, -2` (`z^14`); chart at `L = 0` = `Ginf`, `det H = -25` (`z^60`); `[z^5] Ginf = B_6(1/5)/6 = 4537/3937500` | exact | same |
| L1 | inner step: `N[a nu - e^nu, e^(M nu)] = a^M Poch_M(z/a) exp(Sum b_k (z/a)^(2k-1))`, `a = 1, -2/7, 3/2, 5`, `M = -3..3` | `z^12` | `lf_lmi 12 12` |
| L2 | `g = 0`, `F = 1`: `c^F(z; L) = 1` exactly; 5 values of `alpha0`, `M0 = -2..3`, 2 values of `L` | `z^12` | same |
| L3 | ten random cases (`alpha0` in `{-1, 2, -1/3, 5/2, 5, -7/4}`, `M0` in `{-1, 0, 1, 2}`, quartic `g`, seven-term `log F`): (a) chart = `c^F` at `eps0 = 1/3, -2/5`; (b) chart at `Q = 0` = LMI formula; (c) Method A: no pole, value = formula | `z^12` | same |
| L4 | moment table = Isserlis recursion (degrees `<= 24`, 4 covariances); the generating-function step (`r, a <= 30`); (L2') by the generic engine, `c = 0, 1, -1, 2, 1/5, -3/2` | `z^12` | same |
| L5 | floats: Newton from `(s*, mu_*)` at `Q(L) = a0(L) e^(g'(L)/alpha0)` lands on `(L, -g'(L)/alpha0)` (error `<= 1e-10`, six cases, three `L` each); `Q` is branch-independent | double precision | same |
| M1 | quintic as LMI data at `s* = 1`: `c_n(1 + eps)` regular, `c(z; 1) = h` | `z^14` | `lf_lmi2 14` |
| M2 | quintic as LMI data at `s* = 0`: `c_n(eps)` regular, `c(z; 0) = Ginf` | `z^14` | same |
| M3 | entire `g` (`g_k = r^k/k!`), `log F = log(1 + x/3) + Sum_(m=1,2) (m/2) z^m/(1 + x/2)`, `M0 = 0, 1, -1, 2` | `z^14` | same |
| M4 | beyond LMI: `Gamma(h + a0/z)` with the Stirling series of DLMF 5.11.8, `h = 1/5, 1/3, -1/2, 7/4`: regular, value with `delta = -h z/alpha0`; for integer `h = M` that amplitude is the Poch form (`M = -2..3`) | `z^14` | same |
| F | quintic: `c_n in Q[L]`, exponents `4n - 5k in [0, 4n]`; `c_n(1) = h_n`; `c_n(0) = [z^n] Ginf` (`m <= 6` for `n = 5m`); `c_1 = -(3/20) L^4`, `c_2 = (9/800) L^8` (= CP s.1.7) | `n <= 30`, exact | `lf_fit 30` |
| D1-D3 | 5.11.1 at `w + h` re-expanded = 5.11.8 as polynomials in `h` (`w^-60`); `B_k(h+1) - B_k(h) = k h^(k-1)` (`k <= 61`); LMI value at `s* = 0` = `Ginf` (`z^80`) | exact | `lf_dlmf 60 80` |

### 5.3 Failing controls (each must fail; each did, at the stated place)

| control | where it fails |
|---|---|
| chart with `Q e^(2 mu)` for `Q e^mu` (quintic `L = 2, -2/3`; random case 1) | chart differs from `c^F` at `z^1` |
| `c(z; L)` with the singular (`j = 0`) Stirling series dropped, against the chart | `z^1` |
| twist weight `+5` in `h` (chart at `L = 1`; fit) | `z^1`; `c_n(1) != h_n` |
| orbifold chart, Method A or value with `M0 = 0` (no shift, i.e. `Gamma(5 s)`) against `Ginf` | `z^5` (still regular, as LMI predicts) |
| Bernoulli polynomials at `0`, `-1/5`, `4/5` instead of `1/5` | `z^5`, `z^5`, `z^10` |
| L1: `b_k -> B_2k/(2k)`; `Poch_2` dropped | `z^3`; `z^1` |
| L2, Method A: singular `b_1 + 1/7`; exponent `M0` for `M0 - 1/2`; Poch dropped (`M0 = 2, -1`); singular Stirling dropped | value `!= 1` at `z^1`; pole `eps^-1` in `c_1` |
| L3: wrong-sign `delta` (all eight cases with `M0 != 0`) | `z^1` |
| L4: wrong-sign `delta`; a pure `m^3` term; an `s^2 m` term (`c = 1, 1/5`) | identity fails |
| D1: sign `(-1)^(k+1)`; `(h + 1/2)` for `(h - 1/2)` | identity fails |
| F: singular Stirling dropped (`c_1` admits no Laurent fit); M4: `B_k(0)` for `B_k(1/5)` | fit not confirmed; pole at `z^1` |

### 5.4 Agreement with the stored numbers
`h_1..h_5` = RA_2 C0; `c_1..c_5(2)` = RA_2 E1 (and `c_5(2)` = CP s.1.7, the Stage-4 value); `c_1 = -(3/20) L^4`, `c_2 = (9/800) L^8`
= CP s.1.7; `[z^5] Ginf = B_6(1/5)/6 = 4537/3937500` = SYN4 s.1.1. My ranges extend RA's: chart = one-variable series `z^14` (RA
`z^9`); (L2') `z^12` (RA `z^8`); `c(z; 0) = Ginf` through `z^30` by the fit (RA `z^18` by a Laurent expansion).

## 6. What I did not check; limits
* The lemmas themselves. I read CP s.1.4, (F1)-(F3), (L1a), (L1b), (L1c), (L3) and their proofs in CP2, and checked that each
  use in RA s.3.3 meets its hypotheses (REFEREE s.3). I did not referee those proofs (RF and RL did, per SYN5). (L3) as proved
  in CP2 s.2.8 rests on DLMF 5.11.1 and 2.3.8; I did not read FSL's formal proof (SYN5 K2).
* DLMF was read through the fetch tool's summariser, not verbatim; the statements of 5.11.1 and 5.11.8 agree with CP2 s.2.8
  and RA s.3.4. My formal proof of 5.11.8 (REFEREE s.4.2) uses Stirling's functional equation in `Q[[1/y]]`, which I did not
  re-prove (it is an (FL) item of CP s.0.3).
* The Stage-4 inputs used after `c(z; 0) = Ginf` (CP s.3.10, P(a), P(b), Theorem D) were not audited.
* All computations are finite. The global fit (F) is evidence: it presumes a Laurent shape in a window of exponents, confirmed
  exactly on three further points; its control shows that it detects a pole. L5 is a floating-point check.
* Non-integer `M0` (M4) is outside Theorem LMI; it is VERIFIED only (REFEREE C6).

## 7. Record

**Files** (all new, in this folder; each under 250 lines; Python 3.14 stdlib only):

| file | content |
|---|---|
| `REFEREE.md`, `REFEREE_2.md` | this report (196 and 115 lines) |
| `lf_core.py` | Bernoulli numbers and polynomials, `b_k`, binomials, series, `TS = Q[eps]/(eps^K)`, `Cyc5 = Q(zeta_5)`, power sums |
| `lf_wick.py` | the functional `N`: one variable, two variables (separated, generic), moments, covariance |
| `lf_models.py` | Taylor data: quintic one-variable series and the two charts (`s* = 1`, `s* = 0`); general LMI data (one variable, chart, Method A, prediction) |
| `lf_log.py` | logging (relative paths only), checks and failing controls |
| `lf_quintic.py`, `lf_lmi.py`, `lf_lmi2.py`, `lf_fit.py`, `lf_dlmf.py` | the checks of s.5 |
| `logs/*.log` | the runs below |

**Runs** (`py -I -B <script> ARGS`, working directory the session scratch folder; at most three of my processes at a time):

| command | checks | time |
|---|---|---|
| `lf_quintic.py 4 10`, `8 30`, `14 60` | 35, 36, 36: ALL CHECKS PASSED | 0.3 s, 3.0 s, 34.1 s |
| `lf_lmi.py 4 4`, `8 8`, `12 12` | 95 each: ALL CHECKS PASSED | 0.8 s, 8.0 s, 40.6 s |
| `lf_lmi2.py 6`, `14` | 16 each: ALL CHECKS PASSED | 1.1 s, 42.8 s |
| `lf_fit.py 6`, `20`, `30` | 6 each: ALL CHECKS PASSED | 0.1 s, 11.7 s, 87.5 s |
| `lf_dlmf.py 10 20`, `60 80` | 7 each: ALL CHECKS PASSED | 0.0 s, 0.4 s |

**Development record (honest).** No check and no positive test failed in any run; every FAIL-type line in the logs is a
labelled control that failed as it must. Before the recorded runs I ran one throw-away smoke test in the scratch folder
(engine, (L3) in the `nu`-form, `C = 1`, `c(z; 2)` through `z^5`), not kept. After the first runs I added one check to
`lf_quintic.py` (the stored `c_1..c_5(2)`) and re-ran its three runs, which overwrote my own three earlier logs (same verdicts,
one more check for `N >= 5`). Since `lf_fit`, `lf_lmi2` and `lf_dlmf` import `h_definition` and `ginf` from `lf_quintic`, I then
re-ran all their runs from the final code as well (overwriting my own logs; same verdicts). Every log is from the final code.
One background shell wrapper of mine reported exit code 1: a `tail` option error AFTER the two Python runs it wrapped, which
both finished with ALL CHECKS PASSED (their logs are complete).

**Compliance.** Only new files, only in this folder (plus throw-away output in the session scratch folder). No existing project
file was modified, moved or deleted; no project, worker, referee or critic code was read or imported. No `__pycache__` (`-B`). Every process I started has finished; none of
anyone else's was touched. No wolframscript, no Mathematica, no cluster job. Web: DLMF s.5.11 (two reads through the fetch
tool), no personal data in any request. Logs and notes contain no absolute path and no user name.
