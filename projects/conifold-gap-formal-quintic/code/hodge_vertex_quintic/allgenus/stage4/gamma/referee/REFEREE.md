# Referee report: Stage 4, task "gamma" (Gamma / Mellin-Barnes route to D_0 = 1), 2026-11-14

Adversarial referee of `allgenus/stage4/gamma/` (GAMMA.md, GAMMA_PROOFS.md, gm_*.py, five logs). AI-written, AI-checked; no human
has read it. Labels: PROVED (argument written here or in the cited note), VERIFIED (exact computation, range and field stated),
CONJECTURAL, SKETCH, REFUTED. Notation of GAMMA.md: `u = 1 - 5^5 q`, `L = u^(-1/5)`, `s = 1/L`, `x = 1/sigma`, colour 0, `b_k = B_2k/(2k(2k-1))`.
My code: `rf_*.py` in this folder, stdlib only. It imports NO worker module (`gm_*`) and NO project module. Logs: `logs/rf_*.log`.
Primes: p127 = 2^127-1 and p107 = 2^107-1. The worker used p61 and p89, so every modular result below is at two NEW primes.

## 0. Verdict: SOUND WITH CORRECTIONS

* **Nothing is REFUTED.** Every VERIFIED claim was re-derived with independent code in a different formulation, exactly over Q
  and at two new primes, with failing controls. Every number the worker prints is reproduced digit for digit:
  `h_5`, `ph_5`, `x_1 = 4537/3937500 = B_6(1/5)/6`, `x_1 + m5`, `c_5(2) = 11495493197/31500`, the vanishing period `u + 7/10 u^2 + 41/75 u^3`,
  `phi_1`, `phi_2`.
* **Ranges extended** (all PASS):
  * identity (Gamma) `c = R`: exact over Q to `z^40` at six values of `L` (three new), and mod p127, p107 to `z^500`;
  * a seeded random-`L` test mod p127 to `z^80`, which makes `c_n(L) = R_n(L)` an identity of functions of `L` (s.3);
  * free constants `= [z^5m] Ginf` for `m <= 100`;
  * Theorem A(i) to `u^30` in three frames, and at `u = 0` for EVERY frame;
  * Lemma CFU for `k <= 16`; the single-Gamma lemma to `eps^40`.
* **The PROVED claims check out** given their stated inputs, with the corrections of s.5. None of the corrections changes a conclusion.
  * The reduction of s.1, the 1-form of s.2, Lemma CFU and the single-Gamma lemma are correct as written.
  * Theorem A(ii) is correct in structure, with four bookkeeping gaps (s.2.6). One sentence of its proof is wrong as worded:
    "absorbed by D(0) = 1". The factor in question is `(2 pi z/5)^(1/2)`: independent of `u`, but not of `z`. It cancels because it
    is common to both stationary phases.
* **Two referee results go beyond the worker.**
  * **(R-a) Theorem A(i) needs neither Theorem U nor CFU, Brieskorn or Theorem G** (s.2.5). Three inputs suffice: PF for the residue,
    uniqueness of regular solutions (GAUGE Lemma 1), and a FINITE exact computation at `u = 0`. I did that computation for all frames
    at once. So (i), and with it the vanishing-period closed form, is PROVED without Theorem U (given the formal PF derivation).
  * **(R-b) Theorem A(ii) has numerical support far beyond `z^24`.** Its only content is the scalar `c^fS(z) = 1`. Combine my (Gamma)
    data with the sibling task `stage4/evid` and its referee: `D_0 = 1` exact to `z^204` given Theorem U, and free of Lemma N mod two
    primes to `z^100` (canonical gauge computed order by order, not assumed). This gives `c^fS = 1` through `z^100` with the canonical
    gauge computed, and through `z^200` given Theorem U. The worker's "no direct test beyond E4's `z^24`" is superseded.
* **Status of the goal.** `D_0 = 1` is reduced to (Gamma) by Theorem A (PROVED given Theorem U, with the gaps of s.2.6). (Gamma) is
  VERIFIED to `z^500` and remains a SKETCH. The worker's labels are right on that point.

## 1. Independent machinery (what is shared with the worker, and what is not)

| piece | worker (`gm_*`) | referee (`rf_*`) |
|---|---|---|
| WKB series `R_n(L)` | hatted (HPF) recursion in `(zhat, u)` (as ref_rows) | the PF operator applied literally in the FLAT variable `L`: `D = ((L^6 - L)/5) d/dL`, conjugation `zD -> L + z delta`; transport operator `T_1 = (L^6 - L) d/dL` derived and asserted; `R_n` by exact division by `L^5 - 1` and integration |
| QRR constants `h_n` | cyclotomic power sums from a rational generating function | Galois trace of `(1 - zeta)^-m` in exact `Q(zeta_5)` arithmetic |
| saddle series `c(z;L)` | sigma variable; Taylor data of `S_0, S_1, S_2k`; cyclotomic sums as derivatives of `5 sigma^4/(1 - sigma^5)` | x variable; phase from `W' = x^-2 log((1-x^5)/(1-u))`, made exactly quadratic by Lagrange inversion (no composition); cyclotomic sums by a Newton recurrence for the roots of `(1-x^5)a^5 + 5a^4 + 10a^3 + 10a^2 + 5a + 1`; Kronecker products mod p |
| cross-check of the saddle code | - | second route with EXPLICIT 5th roots of unity in F_p61 and Horner composition: identical through `z^14` (rf_misc M4) |
| residue / Theorem A(i) | compared with the project's `Y^C_Q` (imports `ur_gauge`, `gi_core`) | own Frobenius periods (gamma FOUND as the log-partner condition: `7/10`); criterion = closing relation + no positive `z` + `Y(0) = -e_1`, then GAUGE Lemma 1 |
| Ginf | Bernoulli-polynomial closed form | closed form AND direct re-expansion of `log Gamma(y + 1/5)` in `1/y` |

## 2. The PROVED claims, line by line

### 2.1 The reduction (GAMMA s.1): correct
* `ph_n(0) = f_n` for every choice of free constants: the `u^0` part of the recursion only involves `ph_k(0)`; the free coefficient sits
  at `u^(4m)`, `m >= 1`. Correct. VERIFIED here for `n <= 40` against my own stationary phase of the model (`rf_run_gamma` A2).
  This also corrects the brief: the model fixes nothing at `u = 0`. The sibling tasks `evid` and `poly` reach the same conclusion.
* `R_n in Q[L]`: the quoted theorem was looked up. arXiv:0710.0889 Thm 4(ii) says that for fixed `s` and `n` the coefficient `Phi_s`
  lies in `L Q[L]`. The expansion (5) there is coefficientwise in `x` (their Lemma 1), and Thm 3 gives `Phi_1 = (3/20)(L - L^5)` at `n = 5`.
  Their `F(w,x)` is our colour-0 Frobenius series at `w = 1/z`, `x = q` (checked term by term).
  * **Imprecision (C1).** "our `psi_0/(L e^(u_0/z))` IS `h(z)` times the ZZ series" presupposes (H). The conclusion does not need it.
    The polynomial-recursion WKB equals `d(z) h(z) Sum_s Phi_s z^s/L` for SOME constant series `d`, hence lies in `Q[L]` anyway.
  * Minor: "the residue `[u^(4m)] rho_n = 0` holds" is asserted. It follows from the same citation: a log would contradict
    `Phi_s in L Q[L]`. The evid referee also proves it elementarily.
* `x_m = [L^0] R_5m = R_5m(0)`: correct. The `n = 5` polynomial, `h_5`, `x_1 = B_6(1/5)/6` and `x_1 + m5 = 4437977/4039875000` were all
  reproduced exactly (A3). The criterion "(N) iff `x_m = kappa_m`" is correct and tautological; the content is the identification of `kappa_m`.

### 2.2 The full-Stirling 1-form (GAMMA s.2): correct
Re-derived from the six Stirling series:
* `S_0' = 5 log(5 sigma) - log(1 - sigma^5) + log(-q)`, `S_1 = (1/2) log(5 sigma/(1 - sigma^5))`, `S_2k` as stated.
* The `log z` terms add up to `2 log z`, a constant: the elementary factor is `z^2`.
* With `5^5 q = 1 - u`: `dW/dx = -sigma^2 S_0' = x^-2 log((1-x^5)/(1-u))`.
* `(5 sigma/(1 - sigma^5))^(1/2) d sigma = const (1 - x^5)^(-1/2) dx`.
* `(zeta^j - sigma)^(1-2k) = x^(2k-1) (zeta^j x - 1)^(1-2k)`.
* Saddles `u sigma^5 = 1`; `W(zeta x) = zeta^-1 W(x)`; `A(zeta x; z) = A(x; zeta z)`; `D log(1-u) = 1`.
* The saddle prefactor is `(-W''(s))^(-1/2) A(s;0) = 5^(-1/2) L` (rf_misc M3), so `c` and `R` carry the same leading factor `L`.

**Overstatement (C6).** "the oscillatory-integral picture of UR s.8, now derived". The derived phase `log(1-u)/x + G(x)` is NOT the `W` of
UR s.8 (`(5/x)[1 - ((1-x^5)/(1-u))^(1/5)] - 5 Int t^3/(1-t^5)`). The two only share the critical data (points `x^5 = u`, values
`u_i - u_i^c`). Equivalent up to a CFU-type change of variable, not equal.

### 2.3 PF for both functionals; (H) for `Q^fS` (GAMMA_PROOFS s.1-s.3): correct
* Difference equation `(1 - sigma^5) a(s) = prod_(k=0..4)(5 sigma - kz) a(s-1)` re-derived. It holds for the Stirling integrand because
  `Gamma_St(y+1) = y Gamma_St(y)` in `Q[[1/y]]`, VERIFIED through `v^60` (rf_misc M1).
* **(S4)**: correct. Simpler proof: the `eps^n` coefficient of `F(sigma - eps)` is `(-1)^n F^(n)/n!`, and `SP[F^(n) d sigma] = SP[d F^(n-1)] = 0`
  by (S1). The valuation bound only legitimises the resummation.
* **(R3)**: invariance of the formal residue under `x -> x/(1 - zx)`: standard for Laurent series with a unit change of variable.
* PF for the residue is VERIFIED independently: the closing relation of the rows holds through `u^26` in three frames.
* (H) for `Q^fS`: correct. In hatted variables `A` is a power series in `u` with polynomial coefficients in `(xi, zhat)`. The `j`-sum
  kills the `u^(1/5)` powers: `Sum_j (zeta^j x - 1)^(-m) = (-1)^m 5 Sum_(5|r) C(m-1+r, r) x^r` (VERIFIED, M2). `What''(1; u) = -5/(1-u)` is a unit.
  * **Gap (C4).** "`Y^MB in Q[z^-1][[u]]^5` (only negative powers ...)" is justified only for row 0. For rows 2-4 the absence of positive powers
    of `z` is a CONSEQUENCE of (i), not an input; Theorem G only needs `L[[u]]`. I VERIFIED it for all rows to `u^30`.

### 2.4 Lemma CFU (s.4): PROVED; VERIFIED for `k <= 16` (mod `u^17`)
The induction was re-done:
* `[u^k](What(xi_<k(1); u)) = [u^k] What(1; u)` because `What'(1; u) = 0` for all `u` and `xi_<k(1) = 1` exactly.
* The derivative at `Xi = 1` gives `What''(1;0) phi_k(1) = -5 phi_k(1) = 0`.
* `What(1;u) = -Sum 5u^n/(5n+4) = -(5/4)/tau` is exact.
* `rf_cfu 16` checks: `deg phi_k = k`, `phi_k(1) = 0`, `phi_1`, `phi_2` as stated. The identity `What(Xi psi) = M(Xi)` is re-checked by
  DIRECT substitution at `Xi = 2, 3/2, -1/3, 5/7`. The control `tau = 1` is obstructed at `k = 1` (`r_1(1) = -5/9`).

### 2.5 Theorem A(i): correct, and provable without Theorem U (referee result R-a)
* **Proposition (referee).** For every conifold frame `(alpha, beta)`, `Res_{x=0} Omega_a = Y^C_Q,a`. The proof uses neither Theorem U nor
  Lemma CFU, the Brieskorn reduction or Theorem G.
* **Proof.**
  * `psi := Res Omega` lies in `Q[z^-1][[u]]` and solves PF (s.2.3). Its rows `Y_a` (`Y_0 = psi/C0`, `zD Y_a = c_a Y_(a+1)`) are therefore
    a solution in `L[[u]]^5`: `C0`, `C1` are units, and the pole of `c_1 = 1/(u C0^2 C1^2)` gives the factor `u` in `Y_2`.
  * By GAUGE Lemma 1 (re-checked: `det((n+1)z - B_(-1)) = ((n+1)z)^5` is a unit of `Q((z))`), `Y = Y^C_Q` iff `Y(0) = -e_1`.
  * `Y_a(0)` depends only on `[u^(<=3)] psi`, `C0 mod u^4`, `C1 mod u^3`, so it is a FINITE exact computation.
  * Put `w(f) := max_k (deg_(alpha,beta) [u^k] f - k)`. Then `w(fg) <= w(f) + w(g)`, `w(zD f) <= w(f) + 1`, `w(1/f) <= 0` if `f(0)` is a
    constant and `w(f) <= 0`. From `w(C0) = 0`, `w(t_c) = w(psi) = -1` one gets `w(Y_0..Y_4) <= -1, 0, 0, 1, 2`.
  * So `Y(0)` is a polynomial of degree `<= 2` in `(alpha, beta)`, and a `3 x 3` grid decides it.
  * VERIFIED exactly (`rf_run_res` R3, grid `{0, 1/3, -2} x {0, -2/7, 5}`, all powers of `z`). QED
* Consequently the vanishing-period closed form `varpi_0 t_c = -Sum_k alpha_(k-1) (log(1-u))^k/k!` is PROVED, as `[z^-1]` of (i).
* Independently VERIFIED to `u^30` against my own Frobenius `g` (R1, P1), with `gamma = 7/10` found, not assumed.
* Remark: the formula is also the transfer-theorem image of the Stirling expansion of `(5d)!/(d!)^5`.

### 2.6 Theorem A(ii) and Corollary A: correct in structure, four gaps
Steps (a) pull-back, (b) Brieskorn reduction and (c) one matrix for both functionals were re-derived. The relations
`Xi^n == Xi^(n-5) + (n-3) tau zhat Xi^(n-4)` and `Xi^(m-1) == Xi^(m+4) - (m+1) tau zhat Xi^m` are polynomial in `zhat` and terminate.
Step (d) is correct and is now independent of U (s.2.5). Gaps:
* **(G1) Normalisation factors (C2).** "the elementary normalisations of `Q^m` and `Q^m_nat` differ by colour-wise constants, absorbed by
  `D(0) = 1`" is wrong as worded: the factor is `u`-independent but `z`-dependent. Row 0 at colour 0 gives
  `SP[model row 0] = e^((u_0 - u_0^c)/z) (2 pi z/5)^(1/2) Q^m_0` and `SP[Omega_0] = e^((u_0 - u_0^c)/z) (2 pi z/5)^(1/2) L c(z;L)/C0`.
  The common factor `(2 pi z/5)^(1/2)` (a power of `z`) cancels in
  `P^MB` and in `D = c/R`. The colour phases `zeta^(+-i)` agree because both stationary phases are taken in the same variable `Xi` after
  (S2). The conclusion `D = c/R` stands; the sentence should be replaced by this computation.
* **(G2) Hatted/flat bookkeeping of Res** after the Brieskorn reduction. The worker flagged it. It is well defined: at fixed
  `(u^N, zhat^p)` only finitely many terms meet, because `Phat` at `u^N` is polynomial in `zhat` and `Y^C` is a series in `1/zhat`.
  This should be written out.
* **(G3)** `B_a` polynomial in `zhat`, INCLUDING the `z (D rho_a)` terms of the rows. `D` at fixed `x` produces `-xi/(5u) d/dxi`, which
  is compensated by `z = u^(4/5) zhat` and the hatted row weights. Asserted, not shown.
* **(G4)** Regularity of `P^MB` is taken from (H) for `Q^fS` plus Theorem U. This is fine, but it is the only place where U enters (ii).
* **Numerical status of (ii) (referee result R-b).** Given (H) for `Q^fS`, (ii) is equivalent to `c^fS(z) = c(z;L)/R^can(z;L) = 1`.
  * `c/R^QRR = 1` (here: two primes to `z^500`, exact to `z^40`);
  * `R^can = R^QRR` (sibling `evid`: exact to `z^204` given U, two primes to `z^504`; its referee, free of Lemma N with an own
    canonical-gauge solver: p89 and p61 to `z^100`, exact to `z^40`).
  * Hence `c^fS = 1` through `z^100` (gauge computed) and through `z^200` given U (the (Gamma) leg mod two primes beyond `z^40`).
  * A direct test of the mechanism (build `P^MB` from CFU + Brieskorn to `u^4`, compare `[u^4 zhat^5] Phat_11` with 0) was NOT done; it
    remains the best test of (b)-(c) as such.

### 2.7 The single-Gamma lemma and the two routes to (Gamma)
* **Single-Gamma lemma: PROVED.** The proof uses the functional equation and (S4) in `eta`. The constant `C = 1` comes from
  `(2 pi)^(1/2) Y^(1/2) (2 pi Y)^(-1/2)`. VERIFIED exactly for `r = 0..3` through `eps^40` (A6): the normalised Stirling saddle series
  of `Gamma_St(eta) Y^-eta eta^r` is `Sum_k S(r,k)(-eps)^(r-k)`, i.e. `1`, `1`, `1 - eps`, `1 - 3eps + eps^2`.
* **Formal route: SKETCH (agree).** Step (ii), the reorganisation lemma in `Q[[z/eps, eps^2/z]]`, is missing. The pole bound `p <= n`
  of step (i) is consistent with the diagram count (a phase-only diagram with `V` vertices, `E` propagators scales as `(z/eps)^(E-V)`);
  I did not check it with all amplitude insertions.
* **Analytic route: SKETCH (agree).**
  * (A1) is plausible as described. Two wording points: `E` is `e^(i pi/z)` times a positive function, not positive; and the
    Euler-Maclaurin remainder is small only because the summand varies on the scale `z^(1/2)` near the maximum.
  * (A2) is the real gap: the analytic validity, at fixed `q`, of the coefficientwise expansion (5) of arXiv:0710.0889. No citation
    was found or checked, by the worker or by me.

### 2.8 Remaining PROVED claims; literature
* `(H) + (N) <=> (Gamma)` (Corollary A): correct given Theorem A. Read `R` as the QRR-normalised WKB `h(z) x (Frobenius-normalised)`;
  then "`c(z;1) = h(z)`" needs no (H), and (H) for `R` FOLLOWS from (Gamma) and (H) for `Q^fS`.
* s.8 ("no ODE-formal argument can decide D"): correct for the properties listed (PF, symplecticity, Z/5, (H), polynomiality are all
  invariant under `D -> D exp(odd series in z^5)`). "No purely ODE-formal argument" is broader than what is shown. Theorem A is
  itself formal, with the residue as the exact input.
* Literature, checked:
  * arXiv:0710.0889 Thms 3, 4 and eq. (5): as quoted.
  * arXiv:1803.01409 and arXiv:1702.06096: abstracts only. LP18 s.2.4.4 / Prop. 5 was NOT checked.
  * **Correction (C7).** arXiv:1702.06096 is the local P^2 paper, but its abstract also defines the formal quintic (a twisted theory)
    and states its holomorphic anomaly equations. The brief's citation is therefore defensible; 1803.01409 is the dedicated paper.

## 3. The VERIFIED claims, re-derived (all PASS; my code; logs in `logs/`)

| claim (worker range) | referee range and field | log |
|---|---|---|
| (H) h_auto, `5 not| n` (`n <= 200`, p61/p89; `n <= 30`, Q) | `n <= 500` p127 and p107; `n <= 40` Q | rf_run_gamma_40_500 |
| `R_n in Q[L]`, degrees `0,0,0,2,3,4,4,5,..` (`n <= 30`) | `n <= 40` Q (degrees), `n <= 500` two primes (min exponent `>= 0`) | same |
| `ph_n(0) = f_n` (argument) | `n <= 40` Q, against my own model stationary phase | same, A2 |
| `ph_5`, `h_5`, `x_1 = B_6(1/5)/6`, `x_1 + m5`, `c_5(2)`; evid's `f_10`, `f_15` | exact | same, A3 |
| (Gamma) exact (`n <= 30`; `L = 2, 1/2, 7/3`) | `n <= 40`; `L = 2, 1/2, 7/3, 5/4, 2/3, -3` | same, A4 |
| (Gamma) mod p (`n <= 200`; `L = 2, 3/5`) | `n <= 500`: `L = 2, 3/5, 4/3` (p127), `2, 9/7` (p107); `n <= 80`: six seeded random `L` in F_p127 | same, B2-B3 |
| free constants `= [z^5m] Ginf` (`m <= 8` Q, `m <= 40` two primes) | `m <= 8` Q, `m <= 100` two primes; closed form = direct `Gamma(y + 1/5)` re-expansion to `z^500` | rf_ginf_500 |
| Theorem A(i) (`z^-1..z^-20`, `u^26`, two frames; project `Y^C_Q`) | all powers of `z`, `u^30` (closing `u^26`), three frames; `u = 0` for every frame | rf_run_res_30 |
| vanishing period formula (`u^26`) | `u^30`, own Frobenius, `gamma = 7/10` forced | same, P0-P1 |
| Lemma CFU (`k <= 12`) | `k <= 16`, plus direct substitution at four values of `Xi` | rf_cfu_16 |
| `Gamma_St(y+1) = y Gamma_St(y)`; cyclotomic closed form; saddle prefactor | `v^60`; `m < 30`, `x^60`; five `s` | rf_misc_60 |

* The random-`L` test: `c_n` and `R_n` are rational functions of `L` of degree `O(n)`, equal at six independent random points of
  F_p127 for every `n <= 80`. By Schwartz-Zippel they agree as functions with failure probability `< 10^-30`. So the `L`-independence of
  `c/R`, i.e. the WKB property of the full-Stirling series (Theorem A(a)), is tested and not assumed.
* Two internal cross-checks of my own code: the Newton-sum saddle series equals a second implementation with explicit roots of unity
  in F_p61 (through `z^14`, rf_misc M4); exact rationals reduced mod p agree with the mod-p runs.

## 4. Failing controls (each fails exactly where predicted)

| control | prediction | outcome |
|---|---|---|
| drop the `j = 0` Stirling term (the singular factor `Gamma((1 - sigma)/z)`) | (Gamma) fails | fails at `z^1` (Q) |
| `b_2 -> b_2 + 1` in every Stirling series | first failure at `z^3` | `z^3` |
| `x_1 -> x_1 + 1`; `x_16 -> x_16 + 3` in the QRR series | `D - 1` first non-zero at `z^5`, `z^80` | `z^5`, `z^80` |
| `Ginf + 10^-6` at every `z^5m` | matches no free constant | none matched (two primes) |
| residue without Stirling corrections | Theorem A(i) criterion fails | closing, positivity and `Y(0)` all fail |
| periods with a wrong vanishing period (`gamma + 1/7`) | criterion fails | all three fail |
| residue normalisation `K = 2` | `Y(0) = -2 e_1` | as predicted |
| CFU with `tau = 1` (`v = u`) | obstructed at `k = 1` | `r_1(1) = -5/9` |

## 5. Corrections (overclaims, gaps, errors); none changes a conclusion

* **C1** (s.1). The identification "our WKB IS `h` times the ZZ series" presupposes (H). The conclusion `R_n in Q[L]` holds without it
  (s.2.1). The unobstructed resonance should cite the same theorem.
* **C2** (GAMMA_PROOFS s.5(c)). "differ by colour-wise constants, absorbed by `D(0) = 1`" is wrong as worded. The factor is
  `(2 pi z/5)^(1/2)` (constant in `u` only), common to `Q^fS` and `Q^m_nat`, and it cancels (s.2.6 G1). Replace by that computation.
* **C3** (Theorem A(ii)). Gaps G2-G4 of s.2.6: hatted/flat bookkeeping of the residue, polynomiality of `B_a` including the `z D rho_a`
  terms, and the single use of Theorem U. Label: "PROVED given Theorem U, modulo the bookkeeping G2-G3" (not unconditional).
* **C4** (GAMMA_PROOFS s.2). "`Y^MB in Q[z^-1][[u]]^5` (only negative powers)" is argued only for row 0. It is true (VERIFIED to `u^30`),
  but as a consequence of (i).
* **C5** (GAMMA s.3, s.8, JSON open issue 3). "Theorem A(ii) has no direct numerical test beyond `z^24`" is superseded (s.2.6, R-b).
  Theorem A(i) and the vanishing-period formula are PROVED without Theorem U (s.2.5, R-a).
* **C6** (GAMMA s.0 item 2). "the oscillatory-integral picture of UR s.8, now derived": a different phase with the same critical data.
* **C7** (GAMMA s.6). The arXiv:1702.06096 remark (s.2.8).
* **C8** (scope of the numbers). "(F_g) for `g <= 67`" is correct modulo Theorem U and Theorem A. With (Gamma) to `z^500` it becomes
  `g <= 167` under the same hypotheses. The sibling evid referee has `D_0 = 1 mod z^1005` given Theorem U WITHOUT Theorem A. So the
  gamma route's value is the reduction of `D = 1` to a pure Gamma-function identity, not the numerical range.

## 6. What survives (for the synthesis)

* **PROVED.** Each item is checked here; AI-only throughout.
  * The reduction of s.1: the model fixes nothing at `u = 0`, and the free constants are the `q = infinity` values `R_5m(0)`.
  * The full-Stirling 1-form `Omega`.
  * PF for both formal functionals; (H) for `Q^fS`; Lemma CFU; the single-Gamma lemma.
  * Theorem A(i) for every frame, without Theorem U (referee proof, s.2.5), and with it the closed form of the vanishing period.
  * Theorem A(ii), given Theorem U and modulo G2-G3; hence `D = c/R`.
* **VERIFIED** (independent code, two new primes and exact):
  * (Gamma) through `z^500`, (H) through `n = 500`;
  * free constants `= Ginf` for `m <= 100`;
  * `c^fS = 1` through `z^100` (canonical gauge computed) and through `z^200` given Theorem U (combined with the evid data).
* **SKETCH:** (Gamma) for all orders (both routes). **CONJECTURAL:** `x(z) = Ginf(z)` for all `m`.
* **Net statement.** For all g, Conjecture N(g) follows from: Theorem U, the quoted theorems of STAGE2_SYNTHESIS s.2, Theorem A(ii) (as
  above), and (Gamma), a statement about Stirling series of six Gamma factors. The connection problem is thereby moved entirely to
  `q = 0` (a single-Gamma connection). The decisive missing step is the reorganisation lemma of s.6 (ii), or a citable form of (A2).

## 7. Files and runs (all `py -I -B`, stdlib, working directory outside the project)

| file | content |
|---|---|
| `rf_ring.py` | fields Q, F_p; truncated series; Kronecker products; Bernoulli numbers; `Q(zeta_5)` |
| `rf_wkb.py` | flat-`L` WKB recursion from the PF operator; QRR `h` via `Q(zeta_5)` |
| `rf_sp.py` | full-Stirling saddle series (Lagrange inversion, Newton sums); generic saddle series; explicit-root cross-check |
| `rf_res.py` | own Frobenius periods of L4; residue `psi`; rows and the Lemma-1 criterion |
| `rf_run_gamma.py` | (H), reduction, model, (Gamma) exact and mod p, random `L`, Ginf, single-Gamma lemma, controls |
| `rf_run_res.py` | Theorem A(i) in three frames and on the all-frames grid; vanishing period; controls |
| `rf_cfu.py`, `rf_ginf.py`, `rf_misc.py` | Lemma CFU; free constants vs Ginf to `z^500`; small identities |

* Runs: `rf_run_gamma.py 40 500 80` (869 s), `rf_run_res.py 30` (1 s), `rf_cfu.py 16` (1 s), `rf_ginf.py 500` (17 s), `rf_misc.py 60`.
  Every log ends with ALL CHECKS PASSED.
* Not done: a direct construction of `P^MB` by CFU plus Brieskorn reduction (it would test G2-G3 only). The worker's scripts were not
  re-run: my code reproduces every number in their logs.
