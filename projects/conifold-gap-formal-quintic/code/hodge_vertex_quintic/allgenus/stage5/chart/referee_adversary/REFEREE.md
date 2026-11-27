# Stage 5 / chart: adversarial referee report on the chart argument for (Gamma) (2026-11-14)

Folder `allgenus/stage5/chart/referee_adversary/` (task key `chart/referee_adversary`); new files only. AI-written by the
adversarial referee; no human has read this report. Labels as in the brief (PROVED, VERIFIED, CONJECTURAL, SKETCH, REFUTED).
CP, CP2, CP3, CP4 = `../CHART_PROOF.md`, `_2`, `_3`, `_4`; SYN4 = `../../../stage4/STAGE4_SYNTHESIS.md` (+ `_2`);
GEN = `../../general/GENERAL.md` (+ `GENERAL_2.md`). Continuation: `REFEREE_2.md` (s.4 my computations, s.5 comparison with
`stage5/general`, s.6 what I did not check, s.7 record).

**Read.** All of CP-CP4 and the worker's 16 logs; MECHANISM_STAGE4; SYN4 parts 1-2; POLY s.5-s.6 (for the normalisation of
`Phi^SP` that CP s.3.10 uses); GEN and GEN_2, read only after my own tests had been written and run (s.5). Web (public pages, no
personal data): DLMF s.5.11 and s.2.3(ii). The statements CP2 s.2.8 quotes as 5.11.1 (Stirling's series, `|ph z| <= pi - delta`)
and 2.3.8 (Watson's lemma, hypothesis 2.3.7 with `lambda, mu > 0`) are the ones there; 5.11.8 (shifted argument, Bernoulli
polynomials) is used in s.3.4 below. Not read: EVID, GAMMA, GAMMA_PROOFS, UR, GI and the Stage-4 referee reports (the inputs of
Theorem D, not audited here); arXiv:math/0108100 (CP s.3.8 is a remark used in no proof).

## 0. Verdict: **sound**

The chart argument proves Theorem Gamma (i)-(iv) of CP s.0.1 as stated: PROVED, with (L3) modulo DLMF 5.11.1 and 2.3.8, whose
statements I checked. Theorem D (CP s.0.3) is a correct deduction from its listed hypotheses; those hypotheses (R0, (R1) of UR,
R2, R3, R5, the residue part of (FL), (QRR)) are Stage-4 inputs, quoted by the worker and not audited by me. The precise list:

1. **No mathematical error** in CP s.1-s.4. I re-derived every definition, lemma, step, Hessian, prefactor and branch by hand
   (s.1). The numbers were recomputed with new code and a new method, a Laurent expansion at the colour point (s.4); everything the
   worker prints that my code can reach agrees (`c_1..c_5(2)`, `c_1..c_3(-2/3)`, `c_1, c_2` as polynomials, `h_1..h_5`).
2. **Every attack of the brief fails** (s.2): no hidden convergence; no step valid only for `L != 1` is used at `L = 1`; no
   normalisation can differ by a series in `z`; no branch survives; the 7-variable object and both iterated forms are well
   defined on all of `Lam`; order (a) contains all six Stirling series.
3. **Falsifiers** (REFEREE_2 s.4.4). Seven modifications make the analogue of (Gamma) false: the singular factor's Stirling series
   dropped; its `b_1`, `b_2` or `b_3` perturbed; `b_1` perturbed in all six factors; its exponent `-1/2` replaced by `0`; the
   singular factor squared. In every one, a hypothesis of the argument fails at a named step: (L3) / order (a) for the first six,
   step (c) for the square (the chart critical point runs into `x_0' = 0`). Conversely, every prediction holds exactly (REFEREE_2
   s.4.3): colours 0, 1, 2 of the quintic itself and 36 modified integrands. 17 of these are inside CP's argument as written
   (Gamma factors only, possibly after a Gamma identity). The other 19 need Theorem LMI of s.3 (non-Gamma regular data, or the
   point `q = oo`).
4. **Correct but overcomplete.** Only the factor that degenerates at `q = 0`, `Gamma((1 - sigma)/z)`, needs an Euler integral. A
   two-variable chart (s.3) proves (Gamma) from (L1a), (L1b), (L1c)(i), (L3) for ONE Euler integral, and the delta lemma with the
   roles of the variables exchanged. It needs no 7x7 or 5x5 Hessian, no Schur complement and no (L3) for the five regular factors.
   It proves more (Theorem LMI): the five regular factors may be replaced by ANY holomorphic phase and amplitude. This explains
   the falsifier pattern: perturbing a regular factor never creates a pole; it only changes the value (s.4.3, I1-I4).
5. **By-product** (s.3.4). The same chart at `L = 0` (`q = oo`), with a delta lemma with a source (L2', PROVED in s.3.2), gives
   `c(z; 0) = Ginf(z)`. The closed form `Sum_m f_m z^(5m) = Ginf(z)` of the free coefficients, CONJECTURAL in Stage 4, thus becomes
   PROVED given the inputs of Theorem D and DLMF 5.11.8; VERIFIED exactly through `z^18`. The parallel `stage5/general` agent found
   the same prediction independently (GEN_2 s.5.1 A), with a different, longer sketch (REFEREE_2 s.5).
6. **Remarks R1-R4** (s.1.4). None is needed for correctness.

## 1. What I checked by hand

### 1.1 Definitions (CP s.1)
* The Wick functional `N_n` of s.1.2 and (F1)-(F3): correct. (F3) is used with `mu = zeta` (s.3.9) and `mu = L^-4` (s.3.10); both
  uses are as stated.
* s.1.6: `S_0' = lambda + 5 log 5 sigma - Sum_j log(zeta^j - sigma)`, `S_0'' = 5/(sigma (1 - sigma^5))` (from
  `Sum_j 1/(zeta^j - sigma) = 5 sigma^4/(1 - sigma^5)`). The Taylor data of `S_0`, `S_1`, `S_2k` and the generating function of `P_k`
  are right. In the Stirling product the coefficient of `(log z)/z` is `-Sum_j zeta^j = 0`; the coefficient of `log z` is `-1/2 + 5/2 = 2`.
  So `c_n` lies in `Q[L, L^-1, (1 - L^5)^-1]` (Theorem Gamma (i)).
* s.1.7: the reflection identity `log Gamma_St(y) + log Gamma_St(1 - y) = log 2 pi + (y - 1/2)(log y - log(-y))` is right as formal
  Stirling expressions. `e = eps + Sum_j eps_j` is a sum of six odd integers, hence even, and criticality forces `e = 0`. POLY s.5
  (read) uses `G*` (all Bernoulli terms) for all six factors and divides by `sqrt(2 pi w/5)`, so `ph^SP_0 = 1`. CP's identification
  of `c` with GAMMA's and POLY's objects stands.

### 1.2 The lemmas (CP2)
* (S1): Gaussian integration by parts from `<eta_j P> = Sum_k C_jk <d_k P>` and `HC = -1`; the rescaling identity checked. Correct.
* (PG): both sides solve `d/ds X_s(P) = X_s(QP)`. I checked `dC/ds = C delta C`, the double Isserlis step
  `<QP> = (1/2) tr(delta C)<P> + (1/2) Sum (C delta C)_ij <d_i d_j P>`, and `d/ds log d = (1/2) tr(C(s) delta)`. I-adic convergence at
  `s = 1` holds. Correct. (L0): correct.
* (L1a): the homotopy `chi_t` keeps the Gaussian fixed. `I'(t)` is `N` of the amplitude `z^-1 (d_t Phi_t) A_t + d_t A_t`, which is
  legitimate because `N` is linear on `z^-k O[[z]]`. This amplitude is the `chi_t`-pull-back of an exact form, so (S1) kills it. Correct.
* (L1b): after the translation `v_c = 0` the Hessian is block diagonal, (PG) applies over `C[eta_w][[z^(1/2)]]`, and the
  `eta_v`-parity removes half-integer powers at fixed `w` before `w = z^(1/2) eta_w` is substituted. Correct.
* (L1c)(i)-(iii): correct. (iii) and step (e) are used for no theorem.
* (L2): from `C_(xi xi) = 0`, every monomial `s^a xi^b` of `Phi_3^V A_m` has `a <= b/2`, and a non-zero pairing needs `b <= a`. So
  only `V = 0` and the constant term survive. Correct, and in fact EXACT: no corrections in `z`.
* (L3): the reduction to `a > 0` ((L1a), (F3)) leaves two rational sequences. Each equals the asymptotic expansion of Euler's
  integral: Stirling (DLMF 5.11.1) on one side; on the other, the Morse substitution, Watson (DLMF 2.3.8, applied on a finite
  interval, so the error is exponentially small) and (L1a). Uniqueness of asymptotic expansions closes the proof. The analytic
  input is legitimate: it only identifies two sequences of rational numbers. Correct.

### 1.3 The steps (CP3)
* s.3.1: `E = -1 - (1/z)[5 sigma + Sum_j (zeta^j - sigma)] = -1` identically. s.3.2: I recomputed `H_7` (`H_sigma sigma = 0`,
  `H_sigma y = 1/L`, `H_sigma x_j = H_x_j x_j = -1/(zeta^j - L)`, `H_yy = -1/(5L)`), both Schur complements, `det H_7`,
  `det H_ch = -1/(L^2 Lam1^2)`, and its value `-1/25` at `L = 1`.
* s.3.3, order (a): `v_c(sigma) = (5 sigma, zeta^j - sigma)` and `phi = S_0`. (L3) applies with `M = 1` (`dy`) and `M = 0`
  (`dx_j/x_j`), `beta B = const z^-1 A_St`, and all six Stirling series are present.
* s.3.4, order (b): `X_0 = -q y^5/prod_(j>=1) x_j`, `sigma_c = 1 - X_0`, `beta = X_0/X_0*`, the chart phase `psi`, and `K = L - 1`
  at `p_ch` (`q (5L)^5 = L^5 - 1`) are all right.
* s.3.5: both orders give the same `N` by (Fub), an identity of series. `Nch_n` lies in
  `Q(zeta)[L, L^-1, (zeta^j - L)^-1]` because `1/det H_ch = -L^2 Lam1^2` is a polynomial. So two rational functions agree on the
  infinite set `Lam`. Correct.
* s.3.6: `T_k = -Sum_(w_0) a^(1-2k)` since `1 - 2k` is odd; the weight `-5` comes from `Gamma(1 + 5/z)`. s.3.7: I re-derived `I_d`
  (the signs `(-1)^(4d)` and the cancelling powers `z^(5d - d - 4d)`) and the prefactor `L`. s.3.9: `P_k(zeta L) = zeta^-k P_k(L)`,
  the substitution `t = zeta t'` and (F3) with `mu = zeta`. s.3.10: (F3) with `mu = L^-4`. All correct. s.3.11 is a correct
  deduction from R3, R5, s.3.10 and Theorem Gamma.
* s.4, the branches (B1)-(B10): correct. I add a stronger remark. `c(z; L)` depends only on the Taylor data of order `>= 2` of the
  phase at `L` and on the normalised amplitude, so `lambda` and the sign of `q` never enter it. They only fix the dictionary
  `L <-> q`, and `L = 1 <-> q = 0` in either sign convention.

### 1.4 Remarks (none needed for correctness)
* R1. CP s.3.8 (Givental, the Lagrange multiplier) was not checked by me. It is used in no proof.
* R2. CP s.1.2 (F4) cites Hoermander from memory and is used in no proof, as CP says.
* R3. CP2 s.2.8 "Not done" is accurate in substance. A formal proof of (L3) through `Gamma(w+1) = w Gamma(w)` must move the
  expansion point by `O(z)` (`x = a` to `x = a + z`), and this is not one of the operations of s.2. (L1c)(iii) would supply it.
* R4. `c_n in Q[L]` (not just `Q[L, L^-1]`) does not come from CP's chart. There it needs Theorem P(a) and the degree bound
  (SYN4 K2), as CP4 s.5 item 7 says correctly. The two-variable chart at `L = 0` (s.3.4) gives it directly.

## 2. The attacks of the brief, one by one

* **Hidden convergence.** None. Every step except (L3) is an identity of finite Wick sums over the rings of CP s.1.3. In (L3)
  convergent integrals only identify two rational sequences (s.1.2). No asymptotic expansion in `q` is used anywhere.
* **A step valid only for `L != 1` used at `L = 1`.** None. Orders (a) and (b) are used only on `Lam`. The value at `L = 1` is read
  off `Nch`, which is defined AT `L = 1` (non-degenerate chart Hessian, regular Taylor data), and the continuation is of rational
  functions. My Laurent computation (s.4) never uses order (b): it sees the poles cancel inside `c_n(1 + eps)` itself.
* **Normalisation differing by a series in `z`.** Impossible. (Fub) is an identity of the series `N` themselves, with no
  prefactors, and `Nhat = N/N_0` divides by a `z`-independent number. The symbol prefactors are constants times `z^k e^(c/z)`. They
  never enter (Gamma).
* **Branch choices.** See s.1.3, last item: the normalised series carries no branch.
* **Are the formal 7-variable integral and the iterated SPs well defined?** Yes. `det H_7 = 1/(L^2 (1 - L^5)^2)` vanishes nowhere
  on `Lam`. The inner Hessians are `diag(-1/(5 sigma), -1/(zeta^j - sigma))` (order (a), invertible near `L`) and determinant
  `-1/X_0^2` (order (b), `X_0 = 1 - L != 0`). The Taylor data are holomorphic in the outer variables.
* **Does order (a) contain all six Stirling corrections?** Yes (s.1.3). The falsifiers add something sharper. Only the
  SINGULAR factor's corrections are tested by the regularity: dropping or perturbing them gives poles, from `z^1` for `b_1` and from
  `z^(2k-1)` for `b_k`. The five regular factors' corrections enter only the VALUE, and the value test (C1, `c_n(1) = h_n`) checks
  them.
* **Would the same argument prove something false?** No (REFEREE_2 s.4.4, s.4.3). Every falsifier fails at a hypothesis of the
  argument, and every case that satisfies the hypotheses holds. Two of my first-draft "falsifiers" turned out to be covered cases
  in disguise: the exponent `+1/2` is `Gamma(1 + w) = w Gamma(w)`, and `1/Gamma` reduces to `Gamma(1 + (sigma - 1)/z)` by
  reflection. Both hold, with the values the argument predicts (P2b, P10).

## 3. A shorter proof: the two-variable chart (local Mellin inversion)

### 3.1 Statement
Fix `s*`, `alpha0 != 0`, `M0` in Z, `g` holomorphic near `s*`, `F` in `O_(s*)[[z]]` with `F(s*, 0) != 0`; put `a0 := alpha0 (sigma - s*)`.
For `L != s*` near `s*` choose `lambda = -g'(L) - alpha0 log a0(L)`, so that `L` is critical for
`S(sigma) := lambda sigma + g(sigma) + a0 log a0 - a0`. Let `A := a0^(M0 - 1/2) Poch_M0(z/a0) exp(Sum_k b_k (z/a0)^(2k-1)) F(sigma, z)`
(the Stirling form of `Gamma(M0 + a0/z)`, CP2 s.2.8 notation, times `F`), and `c^F(z; L) := Nhat[S, A; L]`.

**Theorem LMI (PROVED in s.3.3 from CP's (L1a), (L1b), (L1c)(i), (L3) and (L2') of s.3.2).** Every `c^F_n(L)` extends
holomorphically to `L = s*`, and
`c^F(z; s*) = exp(gt(delta)/z) F(s* + delta, z)/F(s*, 0)`, with `delta := -M0 z/alpha0` and `gt(x) := g(s* + x) - g(s*) - g'(s*) x`.
Both factors are re-expanded in `z`; for example `gt(delta)/z = Sum_(k>=2) g_k (-M0/alpha0)^k z^(k-1)`. In words: the saddle series,
continued to the colour point, is the normalised residue of the integrand at the first pole `sigma = s* + delta` of the singular
Gamma factor.

### 3.2 (L2'): the delta lemma with a source (PROVED)
*Statement.* Variables `(s, m)`. Phase `P(s) + alpha0 s m`, with `P = Sum_(k>=2) p_k s^k` and no pure `m` term. Amplitude
`e^(cm) A(s, z)`, with `c` in C and `A` in `O_0[[z]]`. Then `0` is a non-degenerate critical point, and
`N = exp(P(delta)/z) A(delta, z)` EXACTLY, with `delta = -cz/alpha0`. For `c = 0` this is CP's (L2) with the roles of multiplier and
constrained variable exchanged.

*Proof.* From `H = [[2p_2, alpha0], [alpha0, 0]]`: `C_ss = 0`, `C_sm = -1/alpha0`, `C_mm = 2 p_2/alpha0^2`. In the scaled variables the
integrand is `e^(J eta_m) Pi(eta_s)`, with `J = c z^(1/2)` and `Pi = exp(Sum_(k>=3) p_k z^(k/2-1) eta_s^k) A(z^(1/2) eta_s, z)`, a function
of `eta_s` alone. Since `C_ss = 0`, in `<e^(J eta_m) eta_s^a>` each `eta_s` pairs with one `eta_m` of the exponential, and the remaining
`eta_m` pair among themselves. Explicitly, `<(J eta_m)^r/r! eta_s^a> = J^r (r-a-1)!! C_sm^a C_mm^((r-a)/2)/(r-a)!` for even `r - a`, and
the sum over `r` is `(J C_sm)^a e^(J^2 C_mm/2)`. Hence `<e^(J eta_m) Pi> = e^(J^2 C_mm/2) Pi(J C_sm)`, termwise in `z^(1/2)` (each order
has finitely many terms). Finally `z^(1/2) J C_sm = delta` and `J^2 C_mm/2 = p_2 delta^2/z`. □
VERIFIED on seeded random data, with controls (`ra_l2src.py`, REFEREE_2 s.4.5).

### 3.3 Proof of Theorem LMI
1. *Two variables.* Write `Gamma(M0 + a0/z)` by its Euler integral, `x = e^nu`. The integrand is
   `Omega_2 := exp([lambda sigma + g(sigma) + a0(sigma) nu - e^nu]/z) e^(M0 nu) F(sigma, z) dsigma dnu`.
   * At fixed `sigma` the `nu`-part has the critical point `nu_c = log a0(sigma)` and Hessian `-a0(sigma)`. By (L1a) (`x = e^nu`) and
     (L3) (`a = a0`, `M = M0`), its `N` is `a0^M0 Poch_M0(z/a0) exp(Sum b_k (z/a0)^(2k-1)) F`.
   * The outer phase is `S`, and `beta = (a0(sigma)/a0(L))^(-1/2)`.
   * The total Hessian at `(L, nu_c(L))` has determinant `-g''(L) a0(L) - alpha0^2`, which is non-zero near `s*`.
   * So (Fub) gives `Nhat[Omega_2; (L, log a0(L))] = c^F(z; L)` for `L != s*`.
2. *Chart.* Put `nu = mu - lambda/alpha0` ((L1a), Jacobian 1) and `Q := e^(-lambda/alpha0)`; a change of branch of `log a0` leaves `Q`
   unchanged. The phase becomes `lambda s* + Psi_Q`, with `Psi_Q := g(sigma) + alpha0 (sigma - s*) mu - Q e^mu`; the amplitude becomes
   `Q^M0 e^(M0 mu) F`. The constants drop out of `Nhat` (F2).
   * `(Psi_Q, e^(M0 mu) F)` is holomorphic in `(sigma, mu, Q)`, `Q = 0` included. At `Q = 0` the critical point is `(s*, mu_*)`,
     `mu_* = -g'(s*)/alpha0`, and the Hessian `[[g''(s*), alpha0], [alpha0, 0]]` has determinant `-alpha0^2 != 0`.
   * So (CP s.1.4) the critical points `(sigma_c(Q), mu_c(Q))` form a holomorphic, non-degenerate family. Moreover
     `alpha0 (sigma_c - s*) = Q e^(mu_c)`, so `L := sigma_c(Q) = s* + Q e^(mu_*)/alpha0 + O(Q^2)` is a coordinate at `Q = 0`.
   * For `Q != 0` this is the critical point of step 1: `e^(nu_c) = Q e^(mu_c) = a0(L)`, and `S'(L) = g'(L) + alpha0 mu_c = 0`.
   * By (L1c)(i), every coefficient of `Nhat[Psi_Q, e^(M0 mu) F]` is holomorphic at `Q = 0`, hence at `L = s*`. By step 1 it equals
     `c^F_n(L)` for `L != s*`.
3. *Value.* At `Q = 0` the phase `Psi_0 = g(sigma) + alpha0 (sigma - s*) mu` is linear in `mu`. Apply (L2') with `s = sigma - s*`,
   `m = mu - mu_*`, `P = gt`, `c = M0`, `A = F(s* + s, z)`. It gives `N = e^(M0 mu_*) exp(gt(delta)/z) F(s* + delta, z)`; divide by
   `N_0 = e^(M0 mu_*) F(s*, 0)`. □

### 3.4 Consequences
* **(Gamma).** Take `s* = 1`, `alpha0 = -1`, `M0 = 0`. Let `g` be the phases of `Gamma(1 + 5s)` and of `Gamma((zeta^j - sigma)/z)`
  for `j >= 1`, and `F` their amplitudes (the `j >= 1` part of `e^(S_1)` and the five Stirling series). Then `c^F = c`, and
  `c(z; 1) = F(1, z)/F(1, 0) = exp(Sum_k b_k T_k z^(2k-1)) = h(z)` (`T_k` as in CP s.3.6). No Euler integral and no (L3) is needed for
  the five regular factors. Taking `s* = zeta^i` gives `c(z; zeta^i) = h(zeta^-i z)` directly.
* **The orbifold point `L = 0` (`q = oo`).** At `sigma = 0` the degenerate factor is `Gamma(1 + 5 sigma/z)`: `s* = 0`, `alpha0 = 5`,
  `M0 = 1`, `delta = -z/5`. The five factors `Gamma((zeta^j - sigma)/z)`, `j = 0..4`, are regular there and supply `g` and `F`.
  * LMI gives `c(z; 0)` = the normalised value of `prod_j Gamma_St((zeta^j - sigma)/z)` at `sigma = -z/5`, i.e. of
    `prod_j Gamma_St(zeta^j/z + 1/5)`.
  * DLMF 5.11.8 holds as an identity of formal series: 5.11.1 at the shifted argument, re-expanded, is also an asymptotic
    expansion of `Ln Gamma(w + h)`, and asymptotic expansions are unique. It turns this value into
    `exp(Sum_n beta_n(1/5) z^n Sum_j zeta^(-jn))`, `beta_n(h) = (-1)^(n+1) B_(n+1)(h)/(n(n+1))`.
  * Only `5 | n` survives: **`c(z; 0) = Ginf(z) = exp(Sum_l (-1)^(l+1) B_(5l+1)(1/5) z^(5l)/(l(5l+1)))`** (PROVED, as LMI).
  * Hence `c_n` has no pole at `L = 0`. With Theorem Gamma (iv), `c_n` lies in `Q[L]`, without Theorem P(a) and SYN4 K2 (cf. R4).
  * By CP s.3.10, `deg ph^SP_n <= floor(4n/5)` and `f^SP_m := [u^(4m)] ph^SP_(5m) = c_(5m)(0) = [z^(5m)] Ginf`.
  * `f^SP = f^can = f^QRR` given R5 and Theorem D. So the closed form, CONJECTURAL in SYN4 s.1.1, is **PROVED given the inputs of
    Theorem D**. VERIFIED exactly through `z^18` (O1, REFEREE_2 s.4.3).
* **What this says about the argument.** Regularity at `q = 0` is a property of the singular factor alone. That factor must be an
  exact `Gamma_St` (falsifiers X1-X3, X5, X7, X8), and it must be the only factor that degenerates there (X4). The regular factors and
  the twist enter only the value. So (Gamma) cannot hold by an accident of the regular data, and an error in them would show in the
  value, which C1 tests.
* **Status.** Theorem LMI and (L2') are PROVED here (AI-written, unrefereed); they use CP's lemmas as proved in CP2. This is a
  suggestion for the write-up: it would replace CP s.3.1-s.3.6 by about one page. CP's 7-variable argument remains correct.
