# Referee report, part 2: citations, corrections, code, own checks, record

Continuation of `REFEREE.md` (s.0 verdict, s.1-s.2 the proofs, s.3 the statements against the chart proof and Theorem P(b)).

## 4. Citations (dlmf.nist.gov, pages s.2.1, s.2.3 and s.5.11 opened 2026-11-14)

| cited in s.9.4 | the note's quote | the source | verdict |
|---|---|---|---|
| DLMF s.2.3(iii), (2.3.13)-(2.3.16) | Laplace's method for `int_a^b e^(-x p(t)) q(t) dt`, minimum at `a`; hypotheses (a)-(c); conclusion (2.3.15); `b_s` from (2.3.16) | as quoted. Also: "Again lambda and mu are positive constants. Also p_0 > 0"; the minimum is put at the left endpoint "without loss of generality"; (c) absolute convergence for all sufficiently large `x` | matches |
| DLMF 5.11.1 | `Ln Gamma(y) ~ (y - 1/2) ln y - y + (1/2) ln 2 pi + Sum_k B_2k/(2k(2k-1) y^(2k-1))`, `abs(ph y) <= pi - delta` | as quoted (`Ln` the general logarithm, `delta > 0`); for real `y > 0` it is the real `ln Gamma` | matches |
| DLMF s.2.1(iii) | uniqueness; addition, multiplication; "substitution, logarithms, and powers" permitted | "Asymptotic expansions of the forms (2.1.14), (2.1.16) are unique"; "Substitution, logarithms, and powers are also permissible" | matches |

* **Hypotheses of Proof B, checked.**
  * On `[1, oo)`: `p' = 1 - 1/x > 0` on `(1, oo)`; `mu = 2`, `lambda = 1`, `p_0 = 1/2`, `q_0 = 1`. The integral converges absolutely for `y > 0`.
  * On `(0, 1]`, after reflection: `p` is decreasing on `(0, 1)`, and `x^(y-1) e^(-yx)` is integrable at 0.
  * I re-derived `b_s = (-1)^s g_s/2` there with the orientation (`dx/dv < 0` on that side). The odd `s` cancel and the even `s` double.
  * Exponentiating the expansion of 5.11.1 is a substitution of an expansion without constant term into the exponential series, so it
    is covered by s.2.1(iii).
  * Only the case `a = 1` is analytic; the general case follows formally from (H). So no uniformity in `a` is involved.
* The coefficients `1, 1/12, 1/288, -139/51840, -571/2488320, 163879/209018880` of DLMF 5.11.4 are reproduced by all three of my
  methods: Wick, Lagrange-Buermann, and `St_0` by its definition.
* The note cites nothing else and says so (no Hoermander). No arXiv statement is cited.

## 5. Corrections (the precise list)

* **C1 (missing statement; trivial).** Add the z-rescaling Lemma H (REFEREE s.3.2). The chart proof uses it as (F3). Theorem P(b)
  uses it in step 4 (`w -> tau w`) and in P(a) (flat to hatted). The note has only the special case (H) of s.9.3.
* **C2 (FL4, missing hypothesis).** The clause "`P(zD_q) SP[..] = SP[P(sigma) ..]` when `D_q F = sigma`" needs `D_q a = 0`; otherwise
  extra terms `z D_q a` appear at each step. The Mellin-Barnes amplitude (the Stirling factors) is `q`-independent, so the use in
  Theorem P(a) is fine (VERIFIED, `rl_fl4`: the control "P(sigma) = sigma + z" fails).
* **C3 (constants).** As written, `e^(P_s/z)` is not an element of `R_n(A)` when `G(0) != 0`. Harmless.
  * L1a, Step 3: `P_s := G o psi_s - G(0) - y^T H y/2`.
  * L1b, Step 2: `G_3 := G - G(0) - y^T H_yy y/2`.
* **C4 (wording).**
  * s.4(iv)(2): "every term of `e^(eps g/z)` has weight `>= 1/2`" should read "every term of `eps g/z`".
  * s.3.2, formal IFT: "the bracket being quadratic in `(c, c')` times `(c - c')`" should read "a matrix with entries in `I`, times
    `(c - c')`".
  * (G5'): "an ideal of topologically nilpotent elements (..., or of weight `>= 1/2`)". The last case is not an ideal of `A`; say
    "entries in the outer variables, of weight `>= 1/2` in the weighted ring".
* **C5 ((G5'), scope).** The statement is for `P` in `R_n(A)`. L1b applies it with `K_1 = K(y) - H_xx` depending on the outer
  variables and with `P` in `R_(n+m)(A)`. The proof is the same (substitute `k = K_1(y)`; both sides are continuous in the weight
  topology), but the statement should include this case. REFEREE s.2 explains why it is not the same as working in `R_n(A[[y]])`.
* **C6 (ledger of R5).** s.0 ("Who uses what") and s.12 list the poly referee's four lemmas (C6) as the formal input of Theorem P(b).
  The complete list is in REFEREE s.3.4. It adds:
  * L3(iv);
  * Lemma T(b);
  * Lemma H;
  * the relabelling of s.10.4, which is a SKETCH.
* **C7 (naming).** "`SP^norm`" is normalised by the Gaussian (`E_zC[1] = 1`; its constant term is `a(c)`). The project's "normalised
  series" has `c_0 = 1`. Say once, in s.1.4, that `SP^norm` is CHART_PROOF's `N`, not its `Nhat`.
* **C8 (what the checks test).**
  * The FL2 check of `chk_res` is close to a tautology: it takes the residue of an explicit X-derivative with the same routine. It
    tests the implementation of the exponentials, not the lemma, which is trivial.
  * All worker data are polynomial phases with small integer coefficients, so log-type (Mellin-Barnes) data are not exercised. My
    checks in s.7 supply them.
* *Not a correction:* the bound of Lemma U is not sharp; the remark (s.1.5) is correct as an upper bound.

## 6. Code review and reproduction

* **Read in full:** `sp_core.py`, `sp_iter.py`, `sp_rand.py`, `chk_util.py` and the six `chk_*.py`.
  * Weight truncation uses the doubled weight: 1 per integration variable, 2 per parameter, 2 per power of `z`.
  * The fix of s.11.3 is correct. `sp_norm` keeps phase data up to `W2 - wmin`, and it asserts a constant covariance when the amplitude
    has negative weight.
  * In `chk_change` the FL1 amplitude is truncated by `c.mul` before the shift by `z^-1`. In principle that could drop terms of doubled
    weight in `(2K, 2K + 2]`. It cannot here: `deg((d_i f) h) <= 7 < 2K` for `K = 6, 8`. No effect.
  * The controls are genuine: each perturbs exactly one ingredient of the lemma.
* **Reproduction.** I re-ran all 15 stored logs from a scratch copy with the stored arguments. Every log is identical to the stored
  one once timing figures are removed:
  * `chk_change` (8 2, 4 1), `chk_delta` (8 2, 4 1), `chk_param` (8, 5), `chk_res` (6 14, 4 10);
  * `chk_stirling` (80 12, 40 12, 10 6);
  * `chk_fubini` (8 3, 8 2, 6 1, 4 1).
  Every stored log ends with `ALL CHECKS PASSED`; there is no `FAIL` line, and no log prints an absolute path.
* **The worker's account of its own failure** (s.11.3) is consistent with the code. The current `sp_norm` contains the
  negative-weight guard, and FL4 is checked in the form multiplied by `z`.

## 7. Own checks (VERIFIED: exact over Q, new code `rl_core.py`, every check with failing controls)

The engine `rl_core.py` computes `SP^norm` by three methods, all different from the worker's `sp_core`:
* Wick, in the scaled variables `t = s eta`, `z = s^2`, with the exponential computed by its power-series recurrence in `s`;
* the same, with Gaussian moments from the generating function `exp(xi^T C xi/2)` instead of the recursion (W);
* in one variable, Lagrange-Buermann (Morse coordinate).
Coefficients are rationals, or truncated Laurent series in one parameter (with a precision guard).

### 7.1 A toy chart: the s.8 interface run end to end (`rl_toy.py`; `rl_toy_6_14_52.log`, 299 s; also `rl_toy_4_10_40.log`)

The integrand is `e^(Phi/z) dsigma dy dx_0 dx_1/(x_0 x_1)` with
`Phi = lam sigma - y + k sigma log y + Sum_(j=0,1) (-x_j + (c_j - sigma) log x_j)`, `k = 2`, `(c_0, c_1) = (1, -2)`. Its critical
point is `(L, 2L, 1 - L, -2 - L)`. This is the quintic's 7-variable integrand with two colours instead of five and twist 2 instead of
5. No lemma uses `Sum c_j = 0`. The toy (Gamma) predicted by the chart argument reads
`c(z; c_0) = h_toy(z) = exp(Sum_k b_k z^(2k-1) [(k c_0)^(1-2k) + (c_1 - c_0)^(1-2k)]) = 1 + z/72 + z^2/10368 - 2731 z^3/11197440 + ..`.

| check | lemmas exercised | range | result |
|---|---|---|---|
| (T) total 4-variable SP (full 4x4 covariance) = (A) the sigma-SP of `e^(S_0/z) A_St`, `A_St = e^(S_1) e^(Sum z^(2k-1) S_2k)` | L1b, (G4), L3 | `z^0..z^6`, `L = 3, -3/2` | PASS |
| (A) by Wick = (A) by Lagrange-Buermann | L1a (Morse) | same | PASS |
| (T) = (B), the chart SP in `(y, x_1)` with phase `f(X_0(w), w)` and amplitude `1/x_1` | L2(3)-(4), L1b | same | PASS |
| the measure `dy/y` throughout: (T) = (A) = (B) | same | same | PASS |
| the chart at `L = c_0` (`q = 0`) = `h_toy` | L1c(ii), (G4), L3 | `z^0..z^14` | PASS |
| (D) `c_n(L)` as Laurent series in `d = L - c_0`: no negative power, and `c_n(c_0) = [z^n] h_toy` | the continuation of s.8 (2) | `n <= 14`; cut-offs 54 and 60 agree on `d^(<=2)` | PASS |
| (D) for the `dy/y` version: regular, value `h_toy` (the chart predicts this, since `St_1 = St_0`) | same | `n <= 14` | PASS |

Every control fails, as it must:
* (A) without `(Delta_in(sigma)/Delta_in(L))^(-1/2)`: first differs at `z^1`;
* (B) without the Jacobian `X_0`: `z^1`;
* (B) with `dy/y` against (T) with `dy`: `z^1`;
* (D) with `b_k -> -b_k`: a pole at `n = 1`;
* (D) with the colour-0 Stirling terms dropped: a pole at `n = 1`.

So regularity at `L = c_0` is not automatic. It holds exactly for the full-Stirling amplitude, as the chart argument predicts.

### 7.2 L3 and Proof A (`rl_l3.py`; `rl_l3_40_16_8.log`, 5 s)

| check | range | result |
|---|---|---|
| L3(i) at `a = 1`: SP by Lagrange-Buermann against `St_m` by its definition, `m = -6..6` (both measures) | `z^0..z^40` | PASS |
| L3(i) by Wick at `a = 3/7, -5/2`, `m = -4, 0, 3` (this tests (H)) | `z^0..z^16` | PASS |
| (iv): `rho = 0`; the three displayed identities of s.9.2 | through `Y^-42` | PASS |
| `St_(m+1) = (1 + m w) St_m`, `m = -6..6` | `w^40` | PASS |
| Step A2 (T-z): Wick SP at the moving point `1 + eps` over `Q[[eps]]`, times `e^(..)(1+eps)^(1/2)`, then `eps -> z`, equals `SP^norm_1[F; 1]` | `z^0..z^8` | PASS |
| L1a with the map `x = e^u`, `m = -3, -1, 0, 2, 3` | `z^0..z^16` | PASS |

Every control fails, as it must:
* `St_1` claimed for the amplitude `x^1`: `z^1`;
* `b_3 + 1/1000` in `rho`: `Y^-6`;
* T-z without `(Delta_eps/Delta)^(-1/2)`: `z^1`;
* `eps -> 2z`: `z^1`;
* the Jacobian `e^u` omitted: `z^1`.

### 7.3 FL4 and the shift of Theorem P(a) (`rl_fl4.py`; `rl_fl4_8_8.log`, 4 s)

The data are the toy sigma-phase at `L = 3` with its full-Stirling amplitude. The q-dependence is `lam -> lam + eps`, so
`d/deps = D_q`. The moving critical point is computed by my own formal IFT.

| check | range | result |
|---|---|---|
| `z d/deps S_0 = S_1` (`sigma` inserted, at the moving critical point) | all `z^k eps^d` with `k + d <= 8`, `d <= 7` | PASS |
| `(z d/deps)^2 S_0 = S_2` (`sigma^2` inserted) | `k + d <= 8`, `d <= 6` | PASS |
| Lemma T(b) with `p = +z`: `SP[e^((S_0(L+t+z) - S_0(L+t))/z) A_St(L+t+z)] = SP[A_St]` | `z^0..z^8` | PASS |

Every control fails, as it must:
* without `e^((F_eps(c) - F_0(c_0))/z)`;
* without `(Delta(eps)/Delta(0))^(-1/2)`;
* the claim `P(sigma) = sigma + z`;
* T(b) without the exponential factor (`z^1`).

The comparison region `k + d <= K` is forced: `N_k` is computed only for `k <= K`, and `X^n = O(eps^n)` couples `z^-n` to `eps^n`.

### 7.4 A blind control, and my own failures

* **A blind control** (`rl_blind.py`; `rl_blind_9.log`). In the toy, replacing the chart amplitude `1/x_1` by `1/(X_0 x_1)`, which is
  proportional to `y^-2` (dropping the Jacobian), gives the SAME normalised chart series at `L = 2`, through `z^9`.
  * At `L = 3`, `7/4` and `-3/2` the two series differ at `z^1`.
  * I did not find the cause, presumably an identity special to that parameter value.
  * So a control run at a single parameter value can fail to fail. The final toy runs use `L = 3` and `-3/2`.
  * The worker's random-data controls are not affected.
  * The Stage-4 critic's chart controls (`L = 2, 3/2, -2/3`, quintic data) fail at `z^1` as reported, so they are not blind.
* **My own failures.**
  * A first version of `rl_core` raised a `KeyError` on vanishing moments; that run aborted and produced no result. Fixed.
  * A preliminary toy run at `L = 2` reported the blind control above as FAIL. Its log was discarded; s.7.4 records it.
* **Not tested by me:**
  * the quintic data themselves (`Q(zeta)`, 7 variables), which are covered by the chart agent's and the critic's runs;
  * Theorem P(b) as a whole;
  * the relabelling of s.10.4;
  * the chart's branch bookkeeping;
  * complex coefficients.

## 8. What survives

| item | status after refereeing |
|---|---|
| Lemma 1.1, Lemma U, (G1)-(G5), (G5') | PROVED (C3-C5: wording and scope) |
| FL1, IFT, Lemma T, Corollary T', L1c (i)-(iv) | PROVED; L1c(iv) VERIFIED in its L3 instance by a direct Wick computation |
| FL4 | PROVED; the iterated clause needs `D_q a = 0` (C2); VERIFIED on Mellin-Barnes-type data |
| L1a (with Piola), FL3 for the SP | PROVED; VERIFIED with a non-polynomial map |
| L1b, L2 | PROVED; VERIFIED in the chart configuration (toy) |
| L3 (i)-(v), Proof A | PROVED (formal); VERIFIED by a third method |
| L3, Proof B | PROVED modulo DLMF (2.3.15)-(2.3.16), 5.11.1, s.2.1(iii) (quotes checked) |
| FL2, FL3 for the residue | PROVED (trivial) |
| s.10.4, the relabelling | SKETCH (label correct); written out in PB_AUDIT s.2.7, not refereed here |
| Lemma H (z-rescaling), the reflection identity | PROVED here (REFEREE s.3.2-s.3.3) |
| s.8, the interface | correct; matches the chart proof; the toy chart is VERIFIED end to end |
| everything above | AI-written and AI-refereed; no human check |

## 9. Files, runs, compliance

* **New files only**, all in `stage5/lemmas/referee/`:
  * `REFEREE.md`, `REFEREE_2.md`;
  * `rl_core.py`, `rl_toy.py`, `rl_l3.py`, `rl_fl4.py`, `rl_blind.py`;
  * `logs/` with 5 logs: `rl_toy_4_10_40`, `rl_toy_6_14_52`, `rl_l3_40_16_8`, `rl_fl4_8_8`, `rl_blind_9`.
* No existing file was modified, moved or deleted. I re-ran the worker's scripts only from a scratch copy.
* **Runs.** `py -I -B` from the session scratch folder (the copies of these scripts), stdlib only, no project or worker module
  imported. At most two of my processes ran at a time, and none is left running. The logs print relative paths only.
* **Web.** Three DLMF pages, opened to check statements and equation numbers. No personal data was sent in any request.
