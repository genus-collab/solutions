# Stage 4 synthesis (critic): the connection problem `D_0 = 1` (2026-11-14)

Critic of the Stage-4 workflow `allgenus/stage4/`. Three workers, one referee each, all three referees "sound with corrections":
`evid` (falsifier beyond `z^24`, reduced criterion), `gamma` (Gamma / Mellin-Barnes route), `poly` (polynomial connection). The fourth
task key named in the brief, `osc`, has NO folder and no file: nothing was delivered. Its subject, an oscillatory-integral
representation, is covered in substance by gamma s.2-s.3 and poly s.5-s.6.
AI-written; no human has read any Stage-4 note. Labels: PROVED (argument written, or quoted theorem with arXiv number and statement
number), VERIFIED (exact computation; range, field or primes stated), CONJECTURAL, SKETCH, REFUTED. Notation of
`../stage2/STAGE2_SYNTHESIS.md`; UR = `natural_form/theory/UNIFORM_REDUCTION.md`, GI = `../gauge/GAUGE_IDENTITY.md`. As in the brief:
`Q D = P_can Q^m(v(u))`, (H) `D_0 in Q[[z^5]]`, (N) `D_0 = 1`, (F_g) `D_0 = 1 mod z^(3g-2)`. Colour 0, `u = 1 - 5^5 q`, `L = u^(-1/5)`,
`w = zhat`, `h = h_(w_0)` (QRR series). Part 2 (contradictions, overclaims, next steps, critic record): `STAGE4_SYNTHESIS_2.md`.

Read in full: `evid/EVID.md`, `gamma/GAMMA.md`, `gamma/GAMMA_PROOFS.md`, `poly/POLY.md`, the three `*/referee/REFEREE.md`;
`../MECHANISM_STAGE2.md`, `../stage2/STAGE2_SYNTHESIS.md` s.0-s.2, `STAGE2_SYNTHESIS_2.md` s.6; UR s.0-s.5, s.8; GI s.0-s.6; the
docstring of `natural_form/referee/ref_rows.py`. All 73 worker and referee logs were scanned: each ends in PASS, every FAIL line is a
deliberate control, and none prints an absolute path. Critic code (new, `critic/`, `py -I -B`, stdlib, imports no project and no
worker module): `cr4_core.py`, `cr4_sp.py`, `cr4_run.py`, `cr4_chart.py`; logs `critic/logs/` (s.3; record in part 2, s.7).

## 0. Bottom line

* **Nothing is REFUTED.** No order failed in any field, frame or route. Only a READING of the brief is refuted: "the model fixes the
  value of `ph_n` at `u = 0`". In fact `ph_n(0) = F_(0,n)` holds in every normalisation (PROVED three times, s.1.1).
* **(H) and (N) are not proved to all orders.** They are REDUCED, by PROVED steps (given Theorem U), to one identity at `q = 0`:
  **(Gamma)** `c(z; 1) = h(z)`. Here `c(z; L)` is the formal saddle-point series of the Mellin-Barnes summand of the I-function with all
  six Gamma factors replaced by their Stirling series; its `z^n` coefficient is a polynomial in `L`, so `L = 1` (`q = 0`) makes sense.
  * Lemma N + connection form (evid): `D = 1 <=> Phi^can(z, u=1) = h(z)`, `Phi^can` the canonically normalised WKB function.
  * Theorem P(b) (poly) = Theorem A(ii) (gamma): `Phi^can = Phi^SP` (the canonical gauge IS the full-Stirling normalisation). PROVED
    given (U), GI Theorem G and four standard formal lemmas; two independent write-ups, both refereed.
  * Hence `D = 1 <=> (Gamma)` (gamma Cor. A, poly Cor.).
* **(Gamma): CONJECTURAL; the critic gives a complete formal outline** (s.1.4, SKETCH, unrefereed). It is the poly referee's route (ii)
  (the chart argument of arXiv:math/0108100 s.10) made concrete in 7 variables. Formal Fubini identifies `c(z; L)` (`L != 1`) with the
  saddle series of a 5-variable chart integral whose critical point is analytic at `q = 0`. At `q = 0` the chart integral factorises
  into Euler integrals, and their Stirling factors are exactly `h`. No overlap-ring lemma (gamma's formal route) and no analytic
  uniformity (gamma's (A2)) are needed.
  * The step no worker computed, the chart series at `q != 0`, was tested by new critic code. It equals `R_n(L) = L^(4n) ph_n(L^-5)`
    through `z^5` at `L = 2, 3/2, -2/3`, mod `2^61-1` and `2^89-1`, and equals `h` at `L = 1`; `z^5` is the first order that sees (N).
  * If s.1.4 is confirmed, (H) and (N) hold to all orders modulo (U), Theorem G and standard formal lemmas.
* **Evidence, highest orders** (table s.1.5).
  * `D_0 = 1` exactly over Q through `z^204`, given (U): evid worker and evid referee, independent codes. So (F_g) holds for `g <= 69`.
  * `D_0 = 1` mod `2^89-1` and `2^61-1` through `z^1004`, given (U) (evid referee): (F_g) for `g <= 335`.
  * Without Lemma N (canonical gauge computed order by order, regularity asserted): `z^100` at two primes, `z^40` exact.
  * (Gamma), which needs no gauge: `z^500` mod `2^127-1` and `2^107-1` (gamma referee); exact `z^40` at six values of `L` (gamma
    referee) and `z^60` at `u = 1/2, -1/3, 3/7` (critic).
* **The routes agree** (s.3). The three formulas for the normalisation `[z^5q] c(z)` (evid Lemma N/P with `C0`, `C1`; the evid
  referee's period-free formula; poly's first-row formula) are the SAME affine functional of the free coefficients: critic, exact,
  `q <= 12` (`z^60`), two frames, three test functions. gamma's (Gamma), poly's "(A2)" and evid's connection form are one identity.
  The free coefficients agree four ways, exactly: `f^QRR = f^can` (fixed with no QRR input) `= f^SP = [z^5m] Ginf`.
* **Assembled statement** (s.2). Conjecture N(g), for all `g >= 2`, holds modulo (U), the quoted theorems (Q), the `q = 0`
  normalisation (QRR) (now listed explicitly), and either (F_g) as a finite identity (VERIFIED exactly for `g <= 69`, given (U))
  or, for all `g`, (Gamma) together with Theorem P (which adds only GI Theorem G and standard formal lemmas).

## 1. Status of (H) and (N)

### 1.1 The polynomial structure and what the model fixes

* (HPF), `w^(n+1)` coefficient: `(1-u)[(4n - 5 th_u) ph_n + r_n] = 0`. Re-derived independently by the workers, the referees (the
  evid referee on Laurent monomials, the gamma referee in the flat variable `L`, the poly referee at large radius without (HPF)) and
  the critic (`cr4_core.py`, from `M^5 - 1 - q prod_k (5M + kz)`).
* `ph_n(0) = F_(0,n)` (the model's coefficient) in EVERY normalisation (PROVED: evid s.2.1(a), gamma s.1, poly Lemma 1). The `u^0` part
  of (HPF) is the model's own recursion with pivot `4n`. The missing line (evid referee): the scaling limit `v th^5 + 3z th^4 + 1` has a
  unique homogeneous WKB solution.
* `deg ph_n <= floor(4n/5)` for all `n`, all colours and every choice of the free coefficients; every resonance is trivial. PROVED
  twice: elementarily at the orbifold point `x = 1/s` (evid referee s.4), and via arXiv:0710.0889 Thm 4(ii), whose statement and
  hypotheses the gamma and poly referees looked up. Equality holds for `n = 1` and `n >= 3` (`deg ph_2 = 0`). This is VERIFIED only:
  `n <= 1004` at two primes and `n <= 300` exactly (evid referee); `n <= 125` exactly (critic).
* So the free coefficient at `n = 5m` is the TOP coefficient, `f_m = [u^(4m)] ph_(5m) = R_(5m)(L = 0)`. That is the value of the
  colour-0 WKB at the orbifold point `q = oo` (gamma), and the `u^0` term at fixed `z` (evid).
* (H) iff `ph_n(1) = h_n` for every `n` with 5 not dividing `n`. This involves no free parameter.
* (N) iff one affine identity per `q`, with slope 1 in `f_q` (evid s.2.4, poly s.3).
* Closed form, VERIFIED and CONJECTURAL in general: `Sum_m f_m z^(5m) = Ginf(z) = exp(Sum_l (-1)^(l+1) B_(5l+1)(1/5) z^(5l)/(l(5l+1)))`,
  the Stirling factor of `prod_j Gamma(zeta^j/z + 1/5)`. Range: `m <= 100` at two primes (gamma referee), `m <= 25` exactly (critic).
  In particular `f_1 = B_6(1/5)/6 = 4537/3937500`.

### 1.2 The reduction chain (each step PROVED as stated; AI-written, AI-refereed)

| step | statement | hypotheses | source |
|---|---|---|---|
| R0 | `Q D = P_can Q^m(v(u))`, `D` constant diagonal, `D_i(z) = D_0(zeta^-i z)` | (U) | UR s.4 |
| R1 | Lemma P: `[zhat^N](Qhat Rhat(tau zhat)^-1)_ac = [N = a-c mod 5] Sum_n Qh^(0)_(a,n) (-tau)^(N-n) F_(3-c,N-n)` | model symplecticity, Prop. S | evid s.2.2 |
| R2 | Lemma N: every polynomial-type `Phi` is `Phi^can c(z)`, `c in 1 + z^5 Q[[z^5]]`, `[z^5q] c = [u^4q zhat^5q] Ghat_11` | (U), (R1) of UR s.5 | evid s.2.3 |
| R3 | connection form: `D = 1 <=> Phi^can(z, 1) = h(z)`; needs no separate (H) (evid referee) | R2 | evid s.2.4 |
| R4 | three formulas for `[z^5q] c`, and `m5 = [u^4] g_4` | R2, (R1) | evid s.2.3, evid referee s.3, poly s.3 |
| R5 | Theorem P(b) = Theorem A(ii): `Phi^can = Phi^SP` | (U), GI Thm G, (FL) | poly s.6 + poly referee s.3; gamma s.3, GAMMA_PROOFS s.5 |
| R6 | Theorem A(i): `Res_(x=0)` of the rows of the full-Stirling 1-form is `Y^C_Q`, in every frame | no (U) | gamma referee s.2.5 |
| R7 | `D = 1 <=> (Gamma)` | R5 | gamma Cor. A, poly Cor. |

* (FL) are the four standard formal lemmas listed by the poly referee (C6), plus Stirling's formal functional equation
  `Gamma_St(y+1) = y Gamma_St(y)`. The four: formal SP and formal residue kill exact forms; both are invariant under formal changes of
  variables, termwise in `u`; formal SP commutes with `D_q`.
* R5 in poly's version uses only the first row. With the poly referee's bookkeeping (`kappa = C0`, signs `eps_b`) it has no open
  step. gamma's version has the gaps G2-G3 of its referee; they are not needed.

### 1.3 Equivalent forms of (Gamma) (each equivalence PROVED, from R0-R7 or elementary)

* (Gamma1), gamma: `c(z; 1) = h(z)`, where `c` is the flat full-Stirling saddle series with `c_0 = 1`. Equivalently `c = R` for all `L`.
* (Gamma2), poly's "(A2)": `Phi^SP(z, u = 1) = h(z)`. This is the hatted form, `c_n(L) = L^(4n) ph^SP_n(L^-5)`; critic check E3.
* (Gamma3), poly referee: `Phi^SP = h(z) Phi^ZZ`, with `Phi^ZZ` the formal expansion (5) of arXiv:0710.0889.
* (Gamma4), evid connection form, given R5: `ph^can_n(1) = h_n` for all `n`.
* (Gamma5), given (H): `[z^5k] log Phi^can(z, 1) = B_(5k+1)/(5k(5k+1) 5^(5k))` for odd `k` (evid s.8.2). It rests on `p_(5k) = 0`
  for odd `k`, PROVED there.
* (Gamma6), modulo gamma's (A1) (Laplace and Euler-Maclaurin for the lattice sum): the expansion (5) of arXiv:0710.0889 is the
  analytic asymptotic expansion at fixed `0 < q < 5^-5` (gamma's (A2); poly referee's form (ii)).
* Naming hazard: poly's "(A2)" is the identity itself, gamma's "(A2)" is the analytic step of (Gamma6). This synthesis writes (Gamma).

### 1.4 The chart argument for (Gamma) (critic; SKETCH: complete formal outline, unrefereed)

*Integrand.* Write each Gamma factor of the summand `a(s)(-q)^s`, `a(s) = Gamma(1+5s) prod_j Gamma(zeta^j/z - s)`, `sigma = zs`, as an
Euler integral. The powers of `z` cancel because `Sum_j zeta^j = 0`, which leaves the formal 7-variable integrand
`Omega_7 = z^-1 exp(Phi/z) dsigma dy prod_(j=0..4) dx_j/x_j`, with
`Phi = sigma log(-q) - y + 5 sigma log y + Sum_j (-x_j + (zeta^j - sigma) log x_j)`.
Its colour-0 critical point is `x_j = zeta^j - L`, `y = 5L`, `sigma = L` (`u L^5 = 1`); it is non-degenerate for `L != 0, 1`.

* **(a) First `(x, y)` at fixed `sigma`, then `sigma`.** At fixed `sigma` the six integrals separate. Each is an Euler integral, whose
  formal SP is a Stirling series (L3). The SP in `sigma` of their product is the full-Stirling series `c(z; L)` (gamma s.2 = poly s.5).
  This is valid for `L != 1`.
* **(b) First `(sigma, x_0)`.** `Phi` is linear in `sigma`, with coefficient `G = log(-q y^5/prod_j x_j)`, and the amplitude does not
  involve `sigma`. In the coordinates `(sigma, G)` every Wick contraction `<G G>` vanishes. So the formal `(sigma, x_0)`-integral is the
  restriction to `x_0 = X_0 := -q y^5/prod_(j>=1) x_j`, times `2 pi z |x_0|` (L2). The Jacobian `|x_0|` cancels the `1/x_0` of the
  measure. What remains is the chart integrand
  `Omega_ch = 2 pi z (-q)^(1/z) exp([5 log y - y + Sum_(j=1..4) (-x_j + (zeta^j - 1) log x_j) + q y^5/prod_(j>=1) x_j]/z) dy prod_(j>=1) dx_j/x_j`.
  Its critical point `(zeta^j - L, 5L)` and its Hessian are analytic at `q = 0`, where the Hessian is diagonal
  (`-1/(zeta^j - 1)`, `-1/5`). So the normalised SP of `Omega_ch` has coefficients analytic at `q = 0`.
* **(c) Fubini (L1).** (a) and (b) compute the same 7-variable formal SP. So `c(z; L) = SP^norm(Omega_ch)` for `L != 1` near 1. Both
  sides are analytic at `L = 1` (`c_n in Q[L]`), so they also agree at `L = 1`.
* **(d) `q = 0`.** There `Omega_ch` is a product of the five Euler integrals for `Gamma(1 + 5/z)` and `Gamma((zeta^j - 1)/z)`,
  `j = 1..4`. Their Stirling factors multiply to `exp(Sum_k b_k z^(2k-1)[5^(1-2k) + Sum_(j>=1) (zeta^j - 1)^(1-2k)]) = h_(w_0)(z)`.
  Hence `c(z; 1) = h(z)`, which is (Gamma). □
* **(e) Consistency.** Expand `exp(q y^5/(z prod x))` in `q` and use (L3) together with `Gamma_St(y+1) = y Gamma_St(y)`. This gives
  `[q^d] SP(Omega_ch) / [q^0] = (1/d!) prod_(k=1..5d)(5/z + k) / prod_(j>=1) prod_(k=1..d)((1 - zeta^j)/z + k)` EXACTLY, which is the
  Frobenius coefficient of the I-function. So the chart series is `h(z)` times the I-function coefficientwise in `q`: this is (Gamma3).

*Lemmas.* All are standard, and none is written out in the project.
* (L1) Formal SP in several variables at a non-degenerate critical point: invariance under analytic changes of variables; iterated
  equals total (Fubini, by a Schur complement); analytic dependence on parameters, and Taylor expansion in them.
* (L2) The formal delta lemma above, proved by the Wick argument.
* (L3) Stirling's series is the formal SP of Euler's integral, also with `x^m` inserted. gamma's single-Gamma lemma is the
  one-factor case.

*VERIFIED by the critic* (`cr4_chart.py`, log `critic/logs/cr4_chart_5.log`).
* The normalised 5-variable SP of `Omega_ch`, computed by Wick's theorem mod `2^61-1` and `2^89-1`, equals
  `R_n(L) = L^(4n) ph_n(L^-5)` for `n <= 5` at `L = 2, 3/2, -2/3`, and equals `h_n` at `L = 1`.
* At `z^5` it differs from the stored normalisation by exactly `m5`. So this test sees (N).
* Controls: dropping the Jacobian `|x_0|` fails at `z^1`; the measure `dy/y` in place of `dy` fails at `z^2`.

This tests steps (b)-(c) and every normalisation; step (d) is classical.

*Status.* SKETCH (critic; complete outline, AI-written, unrefereed).
* It replaces gamma's two sketched routes (overlap-ring reorganisation; analytic uniformity (A2)).
* The continuation to `q = 0` acts on the coefficients (polynomial in `L`), not on asymptotic expansions, so no uniformity in `q` is
  needed.
* It does not contradict gamma s.8 ("no ODE-formal argument decides D"): its input is the integral representation, not the ODE.

### 1.5 Evidence (highest orders; all PASS, with failing controls)

| identity | highest order | field, frames or points | hypotheses | source |
|---|---|---|---|---|
| `D_0 = 1` ((H) + (N)) | `z^204` | exact Q, frame (0,0) | (U), via Lemma N | evid worker; evid referee (independent code) |
| `D_0 = 1` | `z^1004` | `2^89-1`, `2^61-1` | (U) | evid referee (worker: `z^504`) |
| `D_0 = 1`, Lemma-N-free | `z^100`; `z^40` | two primes; exact | gauge solved, regularity asserted | evid referee (`rf_direct`, `rf_conn`) |
| `D_0 = 1` by first row, Theorem G, closing | `z^150` | two primes, two frames | (U) | poly referee (`rf_gauge`) |
| (H) alone | `n <= 1004`; `n <= 300` | two primes; exact | none | evid referee |
| (Gamma), `c = R` | `z^500` | `2^127-1` (`L = 2, 3/5, 4/3`), `2^107-1` (`L = 2, 9/7`) | none | gamma referee |
| (Gamma) | `z^40`; `z^60` | exact, six `L`; exact, `u = 1/2, -1/3, 3/7` | none | gamma referee; critic `cr4_run` E2 |
| `f^can = f^QRR`, can side with no QRR input | `z^200`; `z^125` | exact | (U) | evid worker and referee; critic D1 |
| chart series `= R(z; L)` (s.1.4 (b)-(c)) | `z^5` | two primes, three `L`; `L = 1` gives `h` | none | critic `cr4_chart` |

The second and third conifold frames are no independent evidence for (N): `K_q` is frame independent given (U) (evid referee s.3;
critic probe C4).

## 2. Conjecture N(g) after Stage 4 (update of `../stage2/STAGE2_SYNTHESIS.md` s.2)

**Assembled theorem (AI-written and AI-refereed only).** For every `g >= 2`, Conjecture N(g) holds modulo exactly the following.
* **(U)** Theorem U (UR s.3). PROVED there, AI-audited only. The evid referee's solver asserts regularity for `k <= 100` at two primes
  and `k <= 40` exactly.
* **(Q)** The quoted theorems of STAGE2_SYNTHESIS s.2:
  * arXiv:1509.06954 Thm 8.1;
  * arXiv:2212.00320 Thm 7.1;
  * arXiv:math-ph/0702045 Thm 5.1, with its s.4;
  * arXiv:1303.1043 Def. 2.2, Prop. 2.4, Def. 2.5, Prop. 2.9;
  * optionally arXiv:2502.18115 Cor. 1.2.
* **(QRR)** The `q = 0` normalisation of the actual R-matrix: `diag(h_(w_i))`, `w_i = {zeta^i - zeta^j} U {-5 zeta^i}`.
  * It was implicit in Stage 2 ("the actual R-matrix"). It is listed now because (N) is exactly a statement about its `z^(5k)` part.
  * Sources as read by the referees:
    * the formal quintic's definition by localisation, arXiv:1702.06096 eq. (5) (evid referee);
    * the Bernoulli form of the `q = 0` R-matrix, arXiv:math/0108100 Thm 9.1 (poly referee).
  * The weight `-5 zeta^i` is pinned by (H) at `n = 1`.
  * The critic did not re-check statement numbers.
* **(F_g)** `D_0 = 1 mod z^(3g-2)`. Discharged in one of two ways:
  * for `g <= 69` (finite): VERIFIED exactly over Q given (U), by two independent codes; for `g <= 335`, mod two primes;
  * for all `g`: R5, Theorem P(b) (PROVED given (U), GI Theorem G and (FL)), plus (Gamma) (CONJECTURAL; complete outline in s.1.4).

Every other step is the Stage-2 chain (STAGE2_SYNTHESIS s.2, steps 1-8), unchanged: PROVED, AI-written, AI-refereed. Genus 2-5 are
also VERIFIED directly by graph sums.

| input | status after Stage 4 | what would close it |
|---|---|---|
| (U) | PROVED (UR s.3), AI-audited | human check |
| (Q) | quoted | hypotheses (STAGE2_SYNTHESIS_2 s.6.2) |
| (QRR) | quoted, standard; project convention since `FORMAL_QUINTIC_genus2.md` s.0 | human check of its `z^(5k)` part (part 2, s.6.2) |
| (FL), (L1)-(L3) | standard, not written out | write them out (Stage 5) |
| (Gamma) | CONJECTURAL; complete outline (s.1.4); VERIFIED to `z^500` (two primes), `z^60` (exact) | write s.1.4 out and referee it (Stage 5) |

**Consequence.** If s.1.4 is confirmed, Conjecture N(g) for all `g >= 2` follows from (U), (Q), (QRR) and standard formal lemmas alone.

## 3. Agreement between the routes (critic re-check; `critic/cr4_run.py`, log `critic/logs/cr4_run_12_25_60.log`, 180 s)

All exact over Q. `Q = 12` (`z^60`), `MCAN = 25` (`z^125`), `NSP = 60`. Smaller runs `4 6 20` and `8 12 40` give identical verdicts.
Every line PASS; ALL CHECKS PASSED.

| check | content | result |
|---|---|---|
| A1-A3 | own periods (`gamma = 7/10` found from the log condition), model rows, QRR constants, `p_m` | `t_c`, `f_0`, `F_0` through `t^7`, `h_1`, `h_5` as in UR/GI; `p_m = 0` for odd `m < 400` exactly at the 40 odd multiples of 5 |
| B1-B4 | QRR polynomials, `n <= 125` | (H) for all 100 `n` with 5 not dividing `n`; resonances trivial; `deg = floor(4n/5)` (`n >= 3`); `ph_n(0) = F_(0,n)`; `f_1, f_2, f_3` = workers' values |
| C1 | three formulas for `[z^5q] c`: `KA` (period free), `KE` (row 1, Lemma P, `C0`, `C1`), `KP` (first row, `C0`) | `KA = KE = KP` for `q <= 12`, frames (0,0) and (1/3,-2/7), on `Phi^QRR` (all 0) and on two perturbations (non-zero; one reproduces the evid referee's printed `rf_pfree` values) |
| C1 | `K_q = [z^5q] Phi(z,1)/h(z)` for all three test functions | PASS (part 2, K9) |
| C2 | `m5 = [u^4] g_4 = p'_(4,1)(0)` from the first row, no free coefficient involved | `-43397/807975000`; stored rows give `KA_1 = KE_1 = KP_1 = m5`, regularity blind |
| C3, C4 | control `tau -> 1`; probe: non-frame vanishing period `g + f_2/7` | `tau = 1` breaks `KA`/REG; the probe still gives `K = 0` (frame insensitivity) |
| D1 | `f^can` from `KA = 0` (no QRR input) | `= f^QRR` for `m <= 25`; `ph^can_n(1) = h_n` for ALL `n <= 125`: `D_0 = 1 mod z^126` |
| D2 | `Ginf` closed form | `f_m = [z^5m] Ginf` for `m <= 25`; `f_1 = B_6(1/5)/6` |
| E1-E3 | hatted full-Stirling SP at fixed `u` (`cr4_sp.py`) | `u = 0`: model `F_0` (`w^60`); `u = 1/2, -1/3, 3/7`: `Phi^SP = Phi^QRR` exactly through `w^60`, `Phi^SP/Phi^(0) in Q[[w^5]]`, `f^SP = f^QRR = f^can`; flat value `c_5(2) = 11495493197/31500` (gamma logs) from `u = 1/32` |
| E4 | controls | `b_2 + 1` fails at `w^3`; dropping the `Gamma(1+5s)` Stirling term fails at `w^1`; QRR weight `+5` breaks (H) at `n = 1` |

**Verdict.**
* The three reductions coincide. They are three readings of one scalar, `[z^5q] P_11(0, z)` of the gauge.
* The integrands of gamma (flat, `exp(W/z)`) and poly (hatted, `exp(-W/w)`) are the same 1-form, with `W_gamma = -W_poly` and equal
  amplitudes.
* The canonical normalisation, the full-Stirling normalisation, the QRR normalisation and the `Ginf` closed form give the same free
  coefficients.
