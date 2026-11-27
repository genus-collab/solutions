# Stage 5 / referees / lmi_adversary: an attempt to break Theorem LMI (2026-11-14)

Folder `allgenus/stage5/referees/lmi_adversary/` (task key `referees/lmi_adversary`); new files only. AI-written by an
adversarial referee; no human has read it. Labels: PROVED, VERIFIED, CONJECTURAL, SKETCH, REFUTED (as in the brief).
RA, RA2 = `../../chart/referee_adversary/REFEREE.md`, `_2.md`; CP, CP2 = `../../chart/CHART_PROOF.md`, `_2.md`;
FSL4 = `../../lemmas/FORMAL_SP_LEMMAS_4.md`; GEN, GEN2 = `../../general/GENERAL.md`, `_2.md`; SYN5 = `../../STAGE5_SYNTHESIS.md` (+ `_2`).
Continuation: `REFEREE_2.md` (s.4 computations, s.5 comparison with GEN s.3.3 and RA2 s.4.4, s.6 corrections, s.7 limits,
s.8 record). Notation of the brief; `x = sigma - s*`, `a0 = alpha0 x`, `e = L - s*`, `delta = -M0 z/alpha0`.

**Read in full:** SYN5 (both parts); RA, RA2; GEN, GEN2; CP s.0-s.1; CP2 (all lemmas); FSL4 s.9 (L3, Proof A). Context:
POLY s.7. Not read: CP3, CP4 (the 7-variable chart, not used by LMI), the referee_formal reports, the other FSL parts.
Web (public pages, no personal data): DLMF s.24.2, s.24.4, s.6.12 through the fetch tool (renderings, not verbatim
typography): 24.2.3 (generating function of `B_n(x)`), 24.2.5 (`B_n(x) = Sum binom(n,k) B_k x^(n-k)`), 24.4.1
(`B_n(x+1) - B_n(x) = n x^(n-1)`), 6.12.1 (`E_1(z) ~ e^-z/z (1 - 1!/z + 2!/z^2 - ..)`).

## 0. Verdict: **sound**

Theorem LMI (RA s.3.1) is true as stated, and RA's proof (s.3.2-s.3.3) is complete, given CP's (L1a), (L1b), (L1c)(i), (L3)
(or FSL's formal (L3)), (F1)-(F2) and the holomorphic implicit function theorem (CP s.1.4). I could not break it. The list:

1. **Every step re-derived** (s.1) against the exact statements of the lemmas it quotes. No error. The deduction of (Gamma)
   (RA s.3.4) is correct: at `L = 1` exactly one factor degenerates, `Gamma(0 + (1 - sigma)/z)`, and the value is `h` by algebra.
2. **Every attack of the brief fails** (s.2): no hidden analytic step; several degenerate factors are excluded by "g
   holomorphic"; non-integer `M0` is excluded by the statement and in fact also TRUE (proved here); the measures are
   encoded by `M0`, and only the singular factor's measure matters at `L = 1`; holomorphy in `Q` is legitimate exactly because
   `a0` is linear in `sigma` and the regular data do not see `lambda`; no branch survives; no series in `z` can intervene.
3. **Falsifiers** (s.3): 27 integrands outside the hypotheses, each with a named hypothesis and a named failing step, each
   failing as LMI's own mechanism predicts (pole order, leading coefficient, residue, log-dependence; the pole orders of X2,
   X3, X3b and the X8 coefficient were observed, not predicted, s.3). Also 9 dangerous-looking positive controls that satisfy the hypotheses, or reduce to them
   by a Gamma identity, and hold (X5+, X7+, X9+, Q2, Q5).
4. **Computations** (REFEREE_2 s.4; own engine, no rescaling, Laurent-series coefficients in `e` with tracked precision; a
   separate two-variable chart engine): exact over Q to `z^10` (random LMI data, steps 1+2 and 3 separately), `z^12` (quintic,
   falsifiers); modulo two new primes to `z^12`, `z^14` (complex `alpha0`, generic `M0`) and `z^16` (quintic, colours 1, 2
   with an explicit `zeta`). All final logs: ALL CHECKS PASSED.
5. **Added here (PROVED, s.2.3):** (L3-C), i.e. (L3) for every `M` in C with Bernoulli polynomials, from FSL's integer case by
   interpolation; hence LMI for every `M0` in C; and DLMF 5.11.8 as a formal identity, so RA's `Ginf` needs no asymptotics.
6. **Scope (PROVED, s.2.5):** LMI holds verbatim when `g` and `F` depend holomorphically on `Q`; it FAILS when they depend on
   `lambda` (`log Q`), or when the singular argument is non-linear in `sigma` (falsifiers X4, X5).
7. **Corrections** (REFEREE_2 s.6, C1-C7): presentation and bookkeeping only; none affects the validity of LMI or of (Gamma).

## 1. The statement and the proof, re-derived

**1.1 Statement.** `s*` in C, `alpha0 != 0`, `M0` in Z, `g` holomorphic near `s*`, `F` in `O_(s*)[[z]]`, `F(s*, 0) != 0`; for `L != s*`
near `s*`, `lambda(L) = -g'(L) - alpha0 log a0(L)`, `S = lambda sigma + g + a0 log a0 - a0`,
`A = a0^(M0 - 1/2) Poch_M0(z/a0) exp(Sum_k b_k (z/a0)^(2k-1)) F`, `c^F(z; L) = Nhat[S, A; L]`. Here `S''(L) = g''(L) + alpha0/e`, so
`-1/S''(L) = -e/(alpha0 + e g''(L))` is holomorphic at `e = 0`; by (F1) every `c^F_n` is holomorphic on a punctured disc
(where `g`, `F_m` are holomorphic and `alpha0 + e g''(L) != 0`) with an isolated singularity at `s*`, a priori a pole.

**1.2 Step 1 ((Fub) with inner variable `nu`).** Checked: `p = (L, nu_c(L))`, `nu_c = log a0(L)` (the branch used in `lambda`), is
critical for `Phi = lambda sigma + g + a0 nu - e^nu` (`d_nu`: `a0 - e^nu`; `d_sigma`: `lambda + g' + alpha0 nu`); `H_vv = -a0(L) != 0`;
`Phi_(sigma sigma) = g''` because `a0'' = 0`, so `det H = -a0(L) g''(L) - alpha0^2 != 0` near `s*`. Inner `N`: (L1a) with `chi = exp`,
`J = e^nu`, `J(p') = a`, gives `N[a nu - e^nu, e^(M0 nu)] = a N[a log x - x, x^(M0-1); a] = a^M0 Poch_M0(z/a) exp(..)` by (L3); `F` is a
`nu`-free factor and `N` is linear over `C[[z]]`. Outer: `phi = S`, `beta = (a0(sigma)/a0(L))^(-1/2)`, `beta B = a0(L)^(1/2) A`, so
`Nhat[Omega_2; p] = Nhat[S, A; L]`: no `z`-dependent factor appears. *Remark (correction C1):* `Omega_2` is a formal integrand.
The real Euler integral of `Gamma(M0 + a0/z)` (after `x -> x/z`) carries `z^(M0 + a0/z)`, whose `sigma`-dependent part
`exp(alpha0 x (log z)/z)` is NOT in `Omega_2` (it would shift `lambda` by `alpha0 log z`). Harmless: `c^F` is defined without it,
and in the quintic the total `(log z)/z` coefficient vanishes identically (CP s.1.6).

**1.3 Step 2 (the chart).** `nu = mu - lambda/alpha0` is an exact identity of functions:
`lambda sigma + alpha0 x (mu - lambda/alpha0) = lambda s* + alpha0 x mu` because `a0` is EXACTLY `alpha0 (sigma - s*)`;
`e^nu = Q e^mu`; `e^(M0 nu) = Q^M0 e^(M0 mu)`. (F2) and `Nhat` remove `lambda s*` and `Q^M0`. After this, `lambda` enters only through
`Q = e^(-lambda/alpha0) = a0(L) e^(g'(L)/alpha0)`: single-valued, holomorphic in `L`, simple zero at `s*`. At `Q = 0` the critical
point is `(s*, mu_*)`, `mu_* = -g'(s*)/alpha0`, with Hessian determinant `-alpha0^2`. For `Q = Q(L)` the point `(L, -g'(L)/alpha0)` is
critical (`alpha0 e - Q(L) e^(-g'(L)/alpha0) = 0`; `g'(L) + alpha0 mu = 0`), tends to `(s*, mu_*)`, hence is the implicit-function
point; it is the translate of `(L, nu_c(L))` (`nu_c + lambda/alpha0 = -g'(L)/alpha0`, same branch). (L1c)(i) gives `Nhat_n(Q)`
holomorphic at `Q = 0`; the normalisation needs `N_0(Q) = e^(M0 mu_c) F(sigma_c, 0) != 0`, true near `Q = 0` since `F(s*, 0) != 0`.
So `c^F_n(L) = Nhat_n(Q(L))` on a punctured disc: the singularity is removable.

**1.4 Step 3 ((L2')).** Re-derived: `C = [[0, -1/alpha0], [-1/alpha0, 2 p_2/alpha0^2]]`; `<eta_m^r eta_s^a> = r!/(r-a)! (r-a-1)!! C_sm^a C_mm^((r-a)/2)`
(each `eta_s` pairs with an `eta_m`); the resummation `(J C_sm)^a e^(J^2 C_mm/2)`; `z^(1/2) J C_sm = delta`, `J^2 C_mm/2 = p_2 delta^2/z`.
With `m = mu - mu_*`, `Psi_0 = g(s*) + gt(x) + alpha0 x m` exactly: no pure-`m` term and no `x^2 m` term, as (L2') requires. Hence
`c^F(z; s*) = exp(gt(delta)/z) F(s* + delta, z)/F(s*, 0)`.

**1.5 (Gamma) (RA s.3.4).** At `s* = 1` the data of CP s.1.6 are LMI data with `alpha0 = -1` (`a0 = 1 - sigma`) and `M0 = 0`: the
`j = 0` factor is `Gamma((1 - sigma)/z)`, coming from `1/Gamma(1 + s - 1/z)` by reflection (CP s.1.7, RL s.3.3), and its part
of `e^(S_1)` is `a0^(-1/2)`. `g` collects the other five phases (holomorphic at 1: `5 != 0`, `zeta^j - 1 != 0`), and
`F = (5 sigma)^(1/2) prod_(j>=1) (zeta^j - sigma)^(-1/2) exp(Sum_k z^(2k-1) b_k [(5 sigma)^(1-2k) + Sum_(j>=1) (zeta^j - sigma)^(1-2k)])`.
Then `c(z; 1) = F(1, z)/F(1, 0) = h`, because `1 - 2k` is odd. At each of `L = 1`, `zeta^i`, `0` exactly one factor degenerates
(`j = 0`, `j = i`, and the twist factor `Gamma(1 + 5s)`, respectively), so LMI applies at all three points.

## 2. The attacks of the brief

**2.1 Hidden analytic step.** None in LMI's own argument. The inputs: (L3), via DLMF 5.11.1 and 2.3.8 (CP: legitimate, it
identifies two rational sequences) or purely formally (FSL4 s.9.3, refereed); the holomorphic IFT (CP s.1.4; a formal IFT
over `C[[Q]]` would do); the removable singularity, which is elementary (`c^F_n` equals an explicit holomorphic function of `L`
on a punctured disc; for the quintic `c_n` is rational, so it has no pole). RA's `Ginf` used DLMF 5.11.8 through uniqueness of
asymptotic expansions (legitimate); s.2.3 replaces it by a formal proof.

**2.2 Several Gamma factors degenerating at `s*`.** Excluded by the hypotheses: a second degenerate factor contributes
`b log b - b` (`b = beta0 x`) to `g`, which is not holomorphic at `s*`; a stray `x^(1/2)` or `1/x` violates "F holomorphic".
With two Euler variables the chart phase at `Q = 0` is `g + alpha0 x mu_1 + beta0 x mu_2`, whose Hessian
`[[g'', alpha0, beta0], [alpha0, 0, 0], [beta0, 0, 0]]` has rank 2, so (L1c)(i) is unavailable. Observed (X1-X3b): poles, with
order exactly `n` and leading coefficient `((2n-1)!!)^2/(n! 16^n) (-1/alpha0)^n` for the square (GEN's `K_0` numbers at
`alpha0 = -1`; the `alpha0`-scaling is my prediction, VERIFIED at `alpha0 = 2`). For `Gamma(a/z) Gamma(-a/z)` the phases cancel,
the saddle is no longer pinned to `s*`, and the poles have order `2n`. In the quintic only one factor degenerates at each
colour point and at `L = 0`.

**2.3 Non-integer `M0`.** The statement assumes `M0` in Z, and RA's `A` uses `Poch_M0`, which is defined only there. The naive
analogue (elementary part of `Gamma(M0 + w)`, Bernoulli part of `Gamma(w)`) FAILS: pole at `z^1` with coefficient
`(beta_1(0) - beta_1(M0))/alpha0` (`= 1/(9 alpha0)` at `M0 = 1/3`; X6m, GEN H2). The correct analogue holds:

**Lemma F5.11.8 (PROVED here).** In `Q[h][[1/w]]`:
`(w + h - 1/2) log(1 + h/w) - h + Sum_k b_k (w + h)^(1-2k) = Sum_(n>=1) beta_n(h) w^(-n)`, `beta_n(h) = (-1)^(n+1) B_(n+1)(h)/(n(n+1))`.
*Proof.* Both sides have coefficients in `Q[h]`. At `h = 0` both are `Sum_k b_k w^(1-2k)` (`B_(n+1)(0) = B_(n+1)`, `B_odd>=3 = 0`). The
left side `Lhs(h)` is `log Gamma_St(w + h) - [(w + h - 1/2) log w - w + (1/2) log 2pi]` with `log(w + h) := log w + log(1 + h/w)`.
Then `Lhs(h + 1) - Lhs(h) = log(1 + h/w)` by FSL4 s.9.1(iv) (`rho = 0`, i.e. `log Gamma_St(Y+1) - log Gamma_St(Y) = log Y`) at
`Y = w + h`, a continuous substitution `1/Y = w^-1 (1 + h/w)^-1`. On the right, DLMF 24.4.1 gives
`Sum_n (-1)^(n+1) (n+1) h^n/(n(n+1)) w^-n = log(1 + h/w)`. By induction the two sides agree for `h = 0, 1, 2, ..`; the coefficients
are polynomials in `h`, so they agree identically. (24.4.1 itself follows from 24.2.3: the difference of the generating
functions is `t e^(xt)`.) QED

**Lemma L3-C (PROVED here).** For every `M` in C: `N[log t - t, t^(M-1); 1](z) = exp(Sum_(n>=1) beta_n(M) z^n)`; hence for `a` in `C^x`,
`N[a log x - x, x^(M-1); a](z) = a^(M-1) exp(Sum_n beta_n(M) (z/a)^n)`. For `M` in Z this is `a^(M-1) Poch_M(z/a) exp(Sum_k b_k (z/a)^(2k-1))`.
*Proof.* By (F1), `[z^n]` of the left side is a polynomial in the Taylor coefficients `binom(M-1, j)` (`j <= 2n`) of the amplitude
`(1 + x)^(M-1)`, hence a polynomial in `M`; the same holds on the right. For `M = m` in Z the left side is `St_m(z)` (FSL4 s.9.1(i),
proved for all `m` in Z), and `log St_m(w) = (1/w + m - 1/2) log(1 + mw) - m + Sum_k b_k w^(2k-1) (1 + mw)^(1-2k)` is Lemma
F5.11.8 at `w -> 1/w`, `h = m`. So the two polynomials agree on Z, hence everywhere. General `a`: CP2 s.2.8 step (1) (`x = a t`, (L1a),
(F3)) uses no integrality. The integer form: `St_(m+1) = (1 + mw) St_m` (FSL4 (iv)) and `St_0 = exp(Sum b_k w^(2k-1))`. QED

**Corollary (PROVED).** Theorem LMI holds for every `M0` in C, with `A := a0^(M0-1/2) exp(Sum_n beta_n(M0) (z/a0)^n) F`. *Proof:* RA's proof
verbatim, with (L3-C) in step 1; `Q^M0 := e^(-M0 lambda/alpha0)` is a constant; (L2') allows any `c` in C. VERIFIED: non-integer and
generic-`F_p` values of `M` and `M0` (REFEREE_2 s.4, E3, P, S1). Consequently RA's `Ginf` derivation (s.3.4) is entirely formal.

**2.4 The measures.** LMI encodes the measure in `M0` (insertion `x^M0`, measure `dx/x`). *Singular factor:* it matters. For
(Gamma) `M0 = 0` (`dx_0/x_0`); `Gamma(1 + (1 - sigma)/z)` (`dx_0`) gives a regular series whose value is the LMI value with `delta = z`,
which differs from `h` at `z^1` (Q3). *Twist factor `Gamma(1 + 5s)` (`dy`) versus `Gamma(5s)` (`dy/y`), regular at `L = 1`:* LMI
predicts the same value `h` (the extra factor `5 sigma` is `z`-free at `sigma = 1`, and `delta = 0`). VERIFIED (Q2); and `c_n(L)` does
change off `L = 1`, for `n >= 2` exactly (in agreement with CP s.1.7). So (Gamma) is blind to the twist measure. The correct `dy` enters
the chain through P(b) and the hatted/flat identity, not through (Gamma). At `L = 0` the twist factor is the singular one:
`M0 = 1` gives `Ginf`; `M0 = 0` gives `exp(5 Sum_(5 | 2k-1) b_k z^(2k-1))` != `Ginf` (Q4, Q5).

**2.5 "Holomorphic in `Q` at `Q = 0`" although `lambda` is singular at `L = s*`.** Legitimate. `lambda ~ -alpha0 log e` is singular,
but after `nu = mu - lambda/alpha0` it survives only in the constant `lambda s*` and in `Q`, which is holomorphic in `L`. Two facts
make this work, and both are essential: (i) `a0` is EXACTLY proportional to `sigma - s*`; (ii) `g`, `F` do not depend on `lambda`.
* Without (i), with `a = alpha0 x + alpha1 x^2`, the Taylor data contain `ell = log a(L)`. At fixed `ell` the series is regular,
  but `[e^1] c_n` depends on `ell` for `M0 = 0`: `c_n = value + O(e log e)`, continuous and not holomorphic (the observed limit is
  `F(s*, z)/F(s*, 0)`). For `M0 = 1` already `[e^0] c_n` depends on `ell`: the series diverges like `log(L - s*)` (X4).
* Without (ii), with `F = exp((2/3) z lambda (1 + x))`: `d c_1(s*)/d ell = 2/3` exactly, a logarithmic divergence (X5).
* With `F` depending holomorphically on `Q` (`F_Q = exp(z Q f_1(sigma))`): regular, with the value of `F_0` (X5+; `M0 = 0, 1, -1`).
  RA's proof covers this verbatim, since the chart data stay holomorphic in `(sigma, mu, Q)`.
In the quintic, `q` enters only through `e^(lambda sigma/z)` and every Gamma argument is linear in `sigma`.

**2.6 Branches.** If `log a0(L)` changes by `2 pi i k`, then `lambda` changes by `-2 pi i k alpha0` and `Q` is unchanged. The sign
of `q` (`lambda -> lambda + i pi`) multiplies `Q` by `e^(-i pi/alpha0)`; this changes only the dictionary `L <-> q` (CP s.4),
not `c` as a function of `L`. `beta(L) = 1`, and the normalised `a0^(M0 - 1/2)` is principal near `sigma = L`. No branch enters `c^F`.

**2.7 A series in `z` instead of a constant.** Impossible. (Fub) and (L3) are identities of series; the only prefactors, `beta`
and `N_0`, are `z`-free; (L2') is exact. VERIFIED: the whole value series agrees, the 2-variable chart equals the 1-variable series
coefficientwise in `z` AND `e` (S1), and three wrong-normalisation controls fail from `z^1` (X10, Q3, Q5).

## 3. Falsifier map (exact over Q, `n <= 12`; REFEREE_2 s.4.4 for the numbers)

| id | integrand (one change) | hypothesis violated | failing step | observed = predicted |
|---|---|---|---|---|
| X1 | `Gamma(a0/z)^2`, `alpha0 = -1, 2` | `g` holomorphic | 2 (rank-2 chart Hessian) | order `n`, `((2n-1)!!)^2/(n! 16^n) (-1/alpha0)^n` |
| X2, X3, X3b | `Gamma(a0/z)` times `Gamma(-a0/z)`, `Gamma(2a0/z)`, `Gamma(1+a0/z)` | `g` holomorphic | 2 | poles: orders `2n`, `n`, `n` (observed) |
| X4 | `a = alpha0 x + alpha1 x^2` | `a0` linear | 2 (data depend on `log Q`) | `e log e` (`M0 = 0`), `log e` (`M0 = 1`) |
| X5 | `F = exp(kappa z lambda (1+x))` | `F` free of `lambda` | 2 | `d c_1(s*)/d ell = -kappa alpha0` |
| X6 | singular `(z/a0)^n0` coefficient `+ eta`, `n0 = 1, 2, 3` | `A` exact Stirling form | 1 ((L3); S4: 1-var != chart) | first pole at `z^n0`, order `n0`, `eta/alpha0^n0` |
| X6m | elementary part `M0 = 1/3`, Bernoulli part `M0 = 0` | same | 1 | pole `z^1`, `1/(9 alpha0)` |
| X7 | exponent `a0^E`, `E != M0 - 1/2` | same | 1 | pole `z^1`, `(1/24 + beta_1(M0) - E^2/2)/alpha0` |
| X8 | `F = 1/(sigma - s*)` | `F` holomorphic | 2 (amplitude not holomorphic) | order `n`, `n! (-1/alpha0)^n` (observed; `E_1`, DLMF 6.12.1) |
| X9 | `F = (sigma - s*) F~`, `M0 != 0` | `F(s*, 0) != 0` | normalisation in 2 (`N_0(Q)` has a simple zero) | simple poles, `[e^-1] c_n = (-M0/alpha0) V_(n-1)` |

Here `V` is the LMI value of `(M0, F~)`. Positive controls (they hold): X5+ (`F` holomorphic in `Q`), X7+ (exponent `+1/2` with
`M0 = 0` is `Gamma(1 + w) = w Gamma(w)`), X9+ (`M0 = 0`: `x Gamma(a0/z) = (z/alpha0) Gamma(1 + a0/z)`, value of `M0 = 1`), Q2, Q5.
The X7 formula vanishes iff `E = +-(M0 - 1/2)` at order `z^1`. Every row fails where LMI says it must; none contradicts LMI.
