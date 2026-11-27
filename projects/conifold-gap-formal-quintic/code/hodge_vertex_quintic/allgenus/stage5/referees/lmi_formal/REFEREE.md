# Stage 5 / referees / lmi_formal: line-by-line referee of Theorem LMI and (L2') (2026-11-14)

Folder `allgenus/stage5/referees/lmi_formal/` (task key `referees/lmi_formal`); new files only. AI-written by the referee; no
human has read this report. Labels as in the brief: PROVED, VERIFIED, CONJECTURAL, SKETCH, REFUTED.
Abbreviations: RA = `../../chart/referee_adversary/REFEREE.md` (s.3 is the object of this report), RA_2 = its `REFEREE_2.md`;
CP, CP2, CP3, CP4 = `../../chart/CHART_PROOF.md`, `_2`, `_3`, `_4`; SYN5 = `../../STAGE5_SYNTHESIS.md` (+ `_2`); GEN =
`../../general/GENERAL.md`. Continuation: `REFEREE_2.md` (s.5 computations, s.6 limits, s.7 record).

**Read.** RA in full (s.3 line by line) and RA_2 in full; CP s.0-s.1; CP2 s.2.1-s.2.8 (every lemma RA uses, with proofs);
CP3 s.3.1-s.3.6, s.3.9-s.3.10; CP4 s.4-s.5; SYN5 parts 1-2 (K6, K7, O6, O14); SYN4 s.1.1 (the `Ginf` line) and SYN4b K2;
GEN s.3.1-s.3.3 (the conjectural form of LMI). Web (public page, no personal data): DLMF s.5.11, through the fetch tool's
summariser (not verbatim): 5.11.1 (Stirling's series, `z -> oo` in `|ph z| <= pi - delta`) and 5.11.8,
`Ln Gamma(z + h) ~ (z + h - 1/2) ln z - z + (1/2) ln(2 pi) + Sum_(k>=2) (-1)^k B_k(h)/(k(k-1) z^(k-1))`, same conditions,
`h` in C fixed (a generalisation note says the earlier restriction `h in [0, 1]` was unnecessary). **Not read:** the `ra_*.py`
scripts, the main session's `independent_check/lmi_check.py` (I read only that folder's README), FSL and RF. My code
(REFEREE_2 s.5) imports none of them.

## 0. Verdict: **sound-with-corrections**

Theorem LMI and (L2') are true as stated, and RA s.3.3 proves LMI from CP's (L1a), (L1b), (L1c)(i), (L3), CP s.1.4 and (F2)
with no gap. Every hypothesis of every lemma is met where it is used (s.3). The statement is right for every `M0` in Z,
including the Pochhammer part and the shift `delta = -M0 z/alpha0` (s.1, s.3.3). The consequences hold: (Gamma) with
`s* = 1, alpha0 = -1, M0 = 0`; `c(z; 0) = Ginf(z)`, with DLMF 5.11.8 used legitimately as a formal identity and the shift `+1/5`
right; `c_n in Q[L]` (s.4). All corrections are presentational; none changes a claim:

* **C1** (step 2, identification of the critical points). The bullet "For Q != 0 this is the critical point of step 1" is
  correct but compressed; its "S'(L) = g'(L) + alpha0 mu_c = 0" holds only with the branch `log a0(L) := mu_c - lambda/alpha0`.
  The full argument is in s.3.2 below (it needs no uniqueness statement).
* **C2** (step 2). "The constants drop out of Nhat (F2)": (F2) removes the phase constant `lambda s*`; the amplitude constant
  `Q^M0` drops out because `Nhat = N/N_0` (for `Q != 0`).
* **C3** (s.3.1). `g_k := g^(k)(s*)/k!` is used but not defined. "For `L != s*` near `s*`" should say why `c^F` is defined
  there: `S''(L) = g''(L) + alpha0/(L - s*)` is non-zero on a punctured neighbourhood.
* **C4** (s.3.4, (Gamma) bullet). "F their amplitudes (the j >= 1 part of e^(S_1) ...)": `F` must contain `(5 sigma)^(1/2)`;
  read "`e^(S_1)` without its `j = 0` factor `(1 - sigma)^(-1/2)`".
* **C5** (s.3.4, orbifold bullet). "The normalised value of `prod_j Gamma_St((zeta^j - sigma)/z)` at `sigma = -z/5`" needs its
  convention stated (division by the non-series part `(2 pi)^(5/2) z^(3/2) e^(C/z)` times a constant); with it the bullet is
  right (s.4.2).
* **C6** (RA_2 s.5, "Alignment"). "For non-integer `M0` ... same proof" is a SKETCH, not a proof: it needs the analogue of (L3)
  for non-integer `M` (DLMF 5.11.8 in place of 5.11.1 and the Pochhammer factor). VERIFIED by me for four non-integer
  values through `z^14` (REFEREE_2 s.5, M4).

Status of the pieces (my labels):

| item | label |
|---|---|
| (L2') with its proof (RA s.3.2) | PROVED (re-derived by hand, s.2); VERIFIED, generic two-variable engine, `z^12` |
| Theorem LMI (RA s.3.1, s.3.3) | PROVED given CP s.1.4, (F2), (L1a), (L1b), (L1c)(i), (L3); (L3) is CP's (PROVED modulo DLMF 5.11.1, 2.3.8) |
| LMI on test data: `g = 0`, `F = 1`, `M0 = -2..3`; ten random polynomial and four entire `(g, F)`, `M0 in {-1, 0, 1, 2}` | VERIFIED exactly: `z^12` (chart and one-variable series), `z^14` (Method A) |
| (Gamma) from LMI (`c(z; 1) = h`, `c(z; zeta^i) = h(zeta^-i z)`) | PROVED given LMI; VERIFIED `z^30` (global fit), `z^60` (chart at `L = 1`) |
| `c(z; 0) = Ginf(z)` | PROVED given LMI and DLMF 5.11.8 (formal); VERIFIED `z^30` (global fit), `z^14` (Laurent), `z^80` (value step) |
| `c_n in Q[L]`, exponents `4n - 5k >= 0` | PROVED given LMI at `s* = 0` and Theorem Gamma (i), (iv); VERIFIED `n <= 30` |
| non-integer `M0` (RA_2 s.5) | SKETCH; VERIFIED (four values of `h`, `z^14`) |

## 1. The statement (RA s.3.1)

Data: `s*`, `alpha0 != 0`, `M0` in Z, `g` holomorphic near `s*`, `F` in `O_(s*)[[z]]`, `F(s*, 0) != 0`, `a0 = alpha0 (sigma - s*)`.
* `S'(sigma) = lambda + g'(sigma) + alpha0 log a0(sigma)` (the derivative of `a0 log a0 - a0` is `alpha0 log a0`), so
  `lambda = -g'(L) - alpha0 log a0(L)` makes `L` critical. Correct.
* `S''(L) = g''(L) + alpha0^2/a0(L) = [g''(L) a0(L) + alpha0^2]/a0(L)`. Non-zero for `0 < |L - s*|` small (C3). The higher
  Taylor coefficients `(a0 log a0)^(k)/k! = (-1)^k alpha0/(k(k-1)(L - s*)^(k-1))` (`k >= 2`) carry no branch; so `c^F(z; L)`
  depends on `L` only, not on `lambda` or on the branch of `log a0`.
* `A` is CP2 s.2.8's symbol `z^(a/z+M) Gamma_St(a/z + M) = (2 pi z)^(1/2) a^(M-1/2) e^((a log a - a)/z) Poch_M(z/a) exp(Sum b_k (z/a)^(2k-1))`
  with `a = a0`, `M = M0`, stripped of `(2 pi z)^(1/2)` and of the phase, times `F`. Correct, both signs of `M0`
  (`Poch_M(e) = prod_(i<M) (1 + i e)`, resp. `prod_(i<=-M) (1 - i e)^-1`, from `Gamma(w + M) = (w)_M Gamma(w)`).
* `S` contains no `log z`; `c^F` is a formal object. Its identification with an actual Mellin-Barnes summand needs the
  cancellation of the `sigma`-dependent powers of `z` (CP s.1.6 for the quintic). RA s.3.4 relies on CP s.1.6; fine.
* The value: `gt(delta)/z = Sum_(k>=2) g_k (-M0/alpha0)^k z^(k-1)` and `F(s* + delta, z) = Sum_n z^n Sum_(m+j=n) F_m^(j)(s*) (-M0/alpha0)^j/j!`
  are power series in `z` with finite sums at each order. The constant term is 1, as it must be.
* *A structural check.* For `g = 0`, `F = 1`, (L1a) with `sigma - s* = (L - s*) t` and (F3) with `mu = a0(L)` give
  `c^F(z; L) = C(z/a0(L))` for one universal series `C` (depending on `M0` only). Regularity at `s*` then forces `C = 1`, and
  LMI asserts exactly this. By hand, `N_1 = 1/12 - 3/8 + 1/2 - 5/24 = 0` (`M0 = 0`) and `1/12 + 1/8 + 0 - 5/24 = 0` (`M0 = 1`).
  VERIFIED `C = 1` exactly through `z^12` for `M0 = -2..3` (REFEREE_2 s.5, L2). This is the formal shadow of Mellin
  inversion, `(1/2 pi i) Int Gamma(M0 + a) x^(-a) da = x^M0 e^(-x)`, whose stationary-phase series has no corrections.
* *A second structural check* (my derivation, not in RA): for `g = 0`, `F = 1 + (sigma - s*)`, `M0 = 1`, `alpha0 = 1`, the
  identity `C = 1` for `M0 = 1, 2` gives `c^F(z; L) = 1 - z/(1 + L - s*)`, whose value at `s*` is `1 - z = F(s* + delta)/F(s*)`
  with `delta = -z`. The sign of `delta` is right.

## 2. (L2') (RA s.3.2): re-derived by hand

* Critical point 0: `d_s = P'(s) + alpha0 m`, `d_m = alpha0 s` vanish at 0 since `P = O(s^2)`. `H = [[2 p_2, alpha0], [alpha0, 0]]`,
  `det H = -alpha0^2 != 0`. `C = -H^-1 = [[0, -1/alpha0], [-1/alpha0, 2 p_2/alpha0^2]]`: `C_ss = 0`, `C_sm = -1/alpha0`,
  `C_mm = 2 p_2/alpha0^2`. Correct.
* `Phi_3 = Sum_(k>=3) p_k s^k` (the term `alpha0 s m` is quadratic and lies in `H`), so after `s = z^(1/2) eta_s`,
  `m = z^(1/2) eta_m` the integrand is `e^(J eta_m) Pi(eta_s)`, `J = c z^(1/2)`. Correct.
* Moments: with `C_ss = 0` each `eta_s` is matched to a distinct `eta_m` (`r!/(r - a)!` injections), the remaining `r - a`
  copies of `eta_m` are matched among themselves (`(r - a - 1)!!` ways), so
  `<eta_m^r eta_s^a> = r!/(r-a)! (r-a-1)!! C_sm^a C_mm^((r-a)/2)` for `r - a` even and `>= 0`, else 0. With
  `(2j-1)!!/(2j)! = 1/(2^j j!)` the sum over `r` is `(J C_sm)^a e^(J^2 C_mm/2)`. Correct (RA's formula, with "`= 0` for `r < a`"
  implicit). VERIFIED exactly for `r, a <= 30`, and the whole moment table against the Isserlis recursion (REFEREE_2 s.5, L4).
* Termwise in `z^(1/2)`: the coefficient of `z^(R/2)` of `e^(J eta_m) Pi` is a finite sum; `Pi(J C_sm)` is evaluated at
  `J C_sm = O(z^(1/2))`. `z^(1/2) J C_sm = -c z/alpha0 = delta` and `J^2 C_mm/2 = p_2 c^2 z/alpha0^2 = p_2 delta^2/z`; with
  `Pi(J C_sm) = exp(Sum_(k>=3) p_k delta^k/z) A(delta, z)` this is `exp(P(delta)/z) A(delta, z)`. EXACT: correct.
* `c = 0`: `N = A(0, z)`, which is CP2 (L2) with `G = alpha0 s`, `f = P`, `X_0 = 0`. The remark is right.
* Controls (REFEREE_2 s.5, L4): the formula with the wrong sign of `delta` fails; a pure `m^3` term or an `s^2 m` term added
  to the phase breaks the identity, as it must (the proof uses that the phase is EXACTLY `P(s) + alpha0 s m`).

## 3. The proof (RA s.3.3), step by step

### 3.1 Step 1: two variables
* `Omega_2` is a formal integrand in the sense of CP s.1.1 (phase and amplitude holomorphic near `(L, log a0(L))`); the
  Euler integral only motivates it. Nothing analytic is used.
* Inner variable `nu` at fixed `sigma` near `L`: `d_nu = a0 - e^nu`, so `nu_c(sigma) = log a0(sigma)` (a branch continuous in
  `sigma`), Hessian `-a0(sigma) != 0` because `a0(L) != 0`. In `m = nu - nu_c` the phase is `const - a0 Sum_(k>=2) m^k/k!`.
* (L1a) with `chi(nu) = e^nu` (biholomorphic near `nu_c`, `J = e^nu`, `J(nu_c) = a0`):
  `N[a0 nu - e^nu, e^(M0 nu); nu_c] = a0 N[a0 log x - x, x^(M0-1); a0] = a0^M0 Poch_M0(z/a0) exp(Sum b_k (z/a0)^(2k-1))` by (L3)
  (`a = a0` in C^x, `M = M0` in Z, measure `dx`). `F` does not depend on `nu` and factors out. Correct. VERIFIED directly in
  the `nu`-form, without (L1a): `a = 1, -2/7, 3/2, 5`, `M = -3..3`, `z^12` (REFEREE_2 s.5, L1).
* (L1b)(i): `H_vv = -a0(L)` is invertible; `phi(sigma) = Phi(nu_c(sigma), sigma) = S(sigma)`. The total Hessian is
  `[[g''(L), alpha0], [alpha0, -a0(L)]]` (`a0'' = 0`), determinant `-g''(L) a0(L) - alpha0^2`, which tends to `-alpha0^2` as
  `L -> s*`; it equals `det H_vv . S''(L)` (Schur, (L1b)(i)). Correct; VERIFIED for the quintic at six `L` (REFEREE_2 s.5, Q2).
* (L1b)(iii): `beta = (a0(sigma)/a0(L))^(-1/2)` (branch `beta(L) = 1`), and `beta B = a0(L)^(1/2) A`. RA omits the constant
  `a0(L)^(1/2)`; it is harmless (`Nhat`). Hence `Nhat[Omega_2; (L, nu_c(L))] = c^F(z; L)` for `L != s*` near `s*`. Correct.

### 3.2 Step 2: the chart
* `nu = mu - lambda/alpha0` (at fixed `L`) is a translation, Jacobian 1. `a0 nu = alpha0 (sigma - s*) mu - lambda sigma + lambda s*`
  and `e^nu = Q e^mu`, so the phase is `lambda s* + Psi_Q` and the amplitude `Q^M0 e^(M0 mu) F`. Correct (C2 for "(F2)").
* *Branch independence of `Q`.* `log a0 -> log a0 + 2 pi i k` changes `lambda` by `-2 pi i k alpha0`, so `-lambda/alpha0` by `2 pi i k`.
  Explicitly `Q(L) = a0(L) e^(g'(L)/alpha0)`: single-valued, holomorphic at `s*`, `Q(s*) = 0`, `Q'(s*) = alpha0 e^(g'(s*)/alpha0) != 0`.
  So `Q` is a holomorphic coordinate at `L = s*`; RA states the equivalent inverse statement. Correct.
* *The family.* At `Q = 0` the critical point is `(s*, mu_*)`, `mu_* = -g'(s*)/alpha0`, with Hessian `[[g''(s*), alpha0], [alpha0, 0]]`,
  determinant `-alpha0^2`; CP s.1.4 gives the holomorphic non-degenerate family `p(Q) = (sigma_c(Q), mu_c(Q))`, `|Q| < r`.
  The critical equations are `g'(sigma) + alpha0 mu = 0` and `alpha0 (sigma - s*) = Q e^mu`, so
  `sigma_c(Q) = s* + Q e^(mu_*)/alpha0 + O(Q^2)` and `sigma_c'(0) != 0`. Correct.
* *The identification (C1), written out.* For `0 < |Q| < r` put `L := sigma_c(Q)`; then `L != s*`, since
  `alpha0 (L - s*) = Q e^(mu_c) != 0`. Take any branch of `log Q` and put `lambda := -alpha0 log Q`, `nu_c := mu_c + log Q`. Then
  `e^(nu_c) = a0(L)` (so `nu_c` is a branch of `log a0(L)`), and the `sigma`-equation gives `lambda = -g'(L) - alpha0 nu_c`: the
  `lambda` of s.3.1 for this branch. The translation maps `(L, nu_c)` to `p(Q)` and `Omega_2` to `(Psi_Q, Q^M0 e^(M0 mu) F)`.
  Step 1 applies at `(L, nu_c)` (`L != s*`; the total Hessian is the chart Hessian at `p(Q)`, non-degenerate). So
  `Nhat[Psi_Q, e^(M0 mu) F; p(Q)] = c^F(z; sigma_c(Q))`. Since `sigma_c` maps the disc biholomorphically onto a neighbourhood
  `U` of `s*`, with `0 -> s*`, we get `c^F_n = chart_n o sigma_c^-1` on `U \ {s*}`. No uniqueness statement is needed. The
  converse, that `(L, -g'(L)/alpha0)` is the chart critical point at `Q(L)`, is VERIFIED numerically (REFEREE_2 s.5, L5).
* *(L1c)(i).* `Psi_Q` and `e^(M0 mu) F` are holomorphic in `(sigma, mu, Q)` near `(s*, mu_*, 0)`, so every `N_n(Q)` is
  holomorphic. `N_0(Q) = e^(M0 mu_c(Q)) F(sigma_c(Q), 0) != 0` near 0, so `Nhat_n(Q)` is holomorphic, and at `Q = 0` it is
  `Nhat_n` of the `Q = 0` data. Composed with `sigma_c^-1`: `c^F_n` extends holomorphically to `s*`. Correct.
* *Explicitly*, by (F1) the chart coefficients are polynomials in the Taylor data of `g`, `F` at `L`, in `Q e^(mu_c) = a0(L)` and
  in the entries of `-H^-1`, with `1/det H = -1/(g''(L) a0(L) + alpha0^2)`: holomorphic at `L = s*` by inspection.
  VERIFIED: chart = `c^F` exactly at `L != s*` (ten random cases at two `L`; the quintic at ten `L`; REFEREE_2 s.5).

### 3.3 Step 3: the value
* With `s = sigma - s*`, `m = mu - mu_*`: `Psi_0 = g(s*) + gt(s) + alpha0 s m`. The linear terms `g'(s*) s + alpha0 mu_* s` cancel,
  there is no pure `m` term and no `s^k m` with `k >= 2`. So (L2') applies with `P = gt`, `c = M0`, `A = F(s* + s, z)`;
  `e^(M0 mu) = e^(M0 mu_*) e^(M0 m)`. `N = e^(M0 mu_*) exp(gt(delta)/z) F(s* + delta, z)`, `N_0 = e^(M0 mu_*) F(s*, 0)`. Correct.
* **The `M0 != 0` question.** The Pochhammer factor enters only in step 1, through (L3), as an exact identity; at `Q = 0` the
  same `M0` reappears as the source `e^(M0 m)`, and (L2') turns it into the shift `delta = -M0 z/alpha0`, the first pole of
  `Gamma(M0 + a0/z)`. The statement is right for every `M0` in Z. VERIFIED: `M0 = -2..3` with `g = 0`, `F = 1` (s.1);
  `M0 in {-1, 0, 1, 2}` with random polynomial data (ten cases) and with entire data (four cases); the wrong sign of `delta`
  fails at `z^1` in every case with `M0 != 0` (REFEREE_2 s.5, L2, L3, M3).

## 4. The consequences (RA s.3.4)

### 4.1 (Gamma)
* `s* = 1`, `alpha0 = -1`, `M0 = 0`, `delta = 0`. `g = [5 sigma log 5 sigma - 5 sigma] + Sum_(j>=1) [(zeta^j - sigma) log(zeta^j - sigma) - (zeta^j - sigma)]`,
  `F = (5 sigma)^(1/2) prod_(j>=1) (zeta^j - sigma)^(-1/2) exp(Sum_k b_k z^(2k-1) [(5 sigma)^(1-2k) + Sum_(j>=1) (zeta^j - sigma)^(1-2k)])` (C4).
  Then `S = S_0`, `A = A_St`, and RA's `lambda` is CP's: `c^F = c` near `L = 1`. `g` is holomorphic and `F(1, 0) != 0` at 1.
* Value `F(1, z)/F(1, 0) = exp(Sum_k b_k T_k z^(2k-1)) = h(z)`, with `T_k = -Sum_(a in w_0) a^(1-2k)` (odd exponent). Correct.
  Here `Q = e^lambda = -q`: the chart parameter is `-q` itself and `Q = 0` is `q = 0`.
* `s* = zeta^i`: `a0 = zeta^i - sigma`, `alpha0 = -1`; value `exp(Sum_k b_k z^(2k-1) [(5 zeta^i)^(1-2k) + Sum_(j != i) (zeta^j - zeta^i)^(1-2k)])`
  `= h_(w_i)(z) = h(zeta^-i z)`. Correct.
* VERIFIED: the chart at `L = 1` equals `h` (computed from its definition in Q(zeta_5)) through `z^60`; Method A at `s* = 1`
  through `z^14`; the global fit gives `c_n(1) = h_n` for `n <= 30` (REFEREE_2 s.5, Q3, M1, F).

### 4.2 The orbifold point `L = 0`; DLMF 5.11.8; the shift by 1/5
* At `sigma = 0` only `Gamma(1 + 5 sigma/z)` degenerates: `a0 = 5 sigma`, `alpha0 = 5`, `M0 = 1` (`(5 sigma)^(1/2) = a0^(M0 - 1/2)`,
  `Poch_1 = 1`); the five factors `j = 0..4` are regular at 0 and give `g`, `F`; `lambda` is CP's; `Q = (-q)^(-1/5)`; `delta = -z/5`.
* *The value as a product of shifted Gamma functions (C5).* With `Gamma_St(a/z) := (2 pi)^(1/2) exp((a/z - 1/2)(log a - log z) - a/z + Sum b_k (z/a)^(2k-1))`
  and `a_j = zeta^j - sigma`: `e^(g/z) F = (2 pi)^(-5/2) z^(-5/2) z^(-5 sigma/z) prod_j Gamma_St(a_j/z)` (`Sum_j zeta^j = 0`).
  At `sigma = delta`: `a_j(delta)/z = zeta^j/z + 1/5`, `z^(-5 delta/z) = z`, `e^(g(delta)/z) = e^(g(0)/z) e^(-g'(0)/5) e^(gt(delta)/z)`.
  So the LMI value is `prod_j Gamma_St(zeta^j/z + 1/5)/[(2 pi)^(5/2) z^(3/2) e^(g(0)/z) e^(-g'(0)/5) F(0, 0)]`. The shift is `h = +1/5`.
* *DLMF 5.11.8 at `w = zeta^j/z`, `h = 1/5`* (`log w = log zeta^j - log z`): the `log z` terms give `z^(-5(h - 1/2)) = z^(3/2)`; the
  `1/z` terms give `e^(Sum_j zeta^j log zeta^j/z) = e^(g(0)/z)`; the constants give `exp((h - 1/2) Sum_j log zeta^j) = e^(-g'(0)/5) F(0, 0)`
  (`g'(0) = -Sum_j log zeta^j`, `F(0, 0) = prod_j zeta^(-j/2)`, same branches). Everything cancels and
  `c(z; 0) = exp(Sum_(n>=1) (-1)^(n+1) B_(n+1)(1/5) z^n Sum_j zeta^(-jn)/(n(n+1)))`; `Sum_j zeta^(-jn) = 5 [5 | n]`, and `n = 5l`
  gives `(-1)^(l+1) B_(5l+1)(1/5) z^(5l)/(l(5l+1))`: **`c(z; 0) = Ginf(z)`**. Correct.
* *5.11.8 as a formal identity.* RA's argument is valid: for fixed `h`, the term-by-term re-expansion of 5.11.1 at `w + h` in
  powers of `1/w` (finitely many terms per order) is an asymptotic expansion, in the scale `w^-n` along `w -> +oo`, of
  `Ln Gamma(w + h) - (w + h - 1/2) ln w + w - (1/2) ln 2 pi`, and so is 5.11.8; uniqueness identifies the coefficients
  (polynomials in `h`). *A purely formal proof* (mine; PROVED given Stirling's functional equation `Gamma_St(y + 1) = y Gamma_St(y)`
  in `Q[[1/y]]`, an (FL) item of CP s.0.3): let `D_n(h)` be the difference of the two coefficients of `w^-n`. Under
  `h -> h + 1` both sides change by the series of `ln(1 + h/w)` (left: the functional equation at `y = w + h`; right:
  `B_k(h + 1) - B_k(h) = k h^(k-1)`), so `D_n(h + 1) = D_n(h)`. At `h = 0` both are 5.11.1 (`B_k(0) = B_k`, `B_odd = 0` for `k >= 3`),
  so `D_n(0) = 0`. A polynomial with these two properties vanishes on all integers `>= 0`, hence identically.
  VERIFIED: as polynomials in `h` through `w^-60`; the difference equation for `k <= 61` (REFEREE_2 s.5, D1, D2).
* VERIFIED, `c(z; 0) = Ginf`: the value step through `z^80` (D3); the chart at `Q = 0` through `z^60` (Q4); Method A at `s* = 0`
  through `z^14` (M2); the global fit, `c_(5m)(0) = [z^(5m)] Ginf` for `m <= 6` (F). Controls: `h = 0` and `h = -1/5` fail at
  `z^5`; `h = 4/5` fails only at `z^10` (`B_n(1 - h) = (-1)^n B_n(h)` agrees at `l = 1`); `M0 = 0` (no shift) fails at `z^5`.

### 4.3 `c_n in Q[L]`, the degree bound, `f^SP`
* **`c_n in Q[L]`: PROVED.** Theorem Gamma (i) (CP s.1.6: `c_n in Q[L, L^-1, (1 - L^5)^-1]`), regularity at `L^5 = 1` (Theorem
  Gamma (iii)-(iv), or LMI at `s* = zeta^i`) and at `L = 0` (LMI at `s* = 0`) leave no pole. Z/5 (CP s.3.9) puts the exponents in
  `4n + 5Z`. Neither P(a) nor SYN4b K2 is used; RA R4 is right. VERIFIED for `n <= 30` (fit with slots for poles; F).
* *Degree bound.* With CP s.3.10, `c_n(L) = L^(4n) ph^SP_n(L^-5)`, the polynomiality of `c_n` gives that `ph^SP_n` has top
  `u`-degree `<= floor(4n/5)` and `[u^(4m)] ph^SP_(5m) = c_(5m)(0)`, with no further input. That `ph^SP_n` has no negative power of `u`
  (`deg c_n <= 4n`) is Theorem P(a) (used by CP s.3.10) or the regularity of the hatted series at `u = 0`. RA cites CP s.3.10,
  so the bullet is correct; its "without Theorem P(a)" applies to `c_n in Q[L]` only. VERIFIED `deg c_n <= 4n`, `n <= 30`.
* `f^SP = f^can = f^QRR`: as SYN5 K7 says, `c(z; 0) = Ginf` itself needs no input of Theorem D; `f^can` needs P(b); `f^QRR` needs
  `D = 1`. RA's "given the inputs of Theorem D" is correct but more than needed (SYN5 O14).
