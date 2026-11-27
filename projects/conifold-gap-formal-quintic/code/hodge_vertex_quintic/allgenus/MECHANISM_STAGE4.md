# All-genus mechanism for Conjecture N: Stage 4, the connection problem (2026-11-14)

Entry point for Stage 4. It updates `MECHANISM_STAGE2.md` s.1 (assembled statement) and s.5 item 1; that note stays unchanged as the
record of Stage 2. Written from an 8-agent workflow (workers evid, gamma, poly, osc; one referee per delivered worker; a critic and
assembler) and from my own re-runs and reading (s.4). Everything is AI-written. Every PROVED step below is AI-refereed, except the chart
argument of s.3, which no one has refereed. No human has checked any proof. Labels as in `../NATURAL_FORM.md`.
Full record: `stage4/STAGE4_SYNTHESIS.md` (s.0-s.3) and `stage4/STAGE4_SYNTHESIS_2.md` (s.4-s.7: eleven contradictions adjudicated,
eighteen overclaims with corrections, next steps). Worker notes: `stage4/evid/EVID.md`, `stage4/gamma/GAMMA.md` and `GAMMA_PROOFS.md`,
`stage4/poly/POLY.md`; referee reports `stage4/*/referee/REFEREE.md`; critic code `stage4/critic/`. Written into the project by me:
this note and `stage4/independent_check/` (s.4.2); my re-run outputs stayed in the session scratch folder.
**Failed agent.** The fourth worker, `osc` (oscillatory-integral route), stopped with an API error: one reply exceeded the 128k
output-token limit. It wrote nothing. Its subject is covered by gamma s.2-s.3, poly s.5-s.6 and the critic's chart argument (s.3).

## 1. The problem and the answer

* The actual R-matrix is normalised at `q = 0` by quantum Riemann-Roch: `R(q=0) = diag(h_(w_i))`, with
  `h_w(z) = exp(-Sum_k B_2k/(2k(2k-1)) z^(2k-1) Sum_(a in w) a^(1-2k))` and `w_i = {zeta^i - zeta^j}_(j != i) U {-5 zeta^i}`.
* Theorem U gives the gauge `P_can`, regular at the conifold `u = 0` and normalised there by `P_11(0,z) = 1`.
* The two differ by a constant diagonal series: `Q D = P_can Q^m(v(u))`, `D_i(z) = D_0(zeta^-i z)`. Conjecture N(g) needs
  (F_g): `D_0 = 1 mod z^(3g-2)`.
* **Stage 4 result.**
  * `D_0 = 1` holds exactly over Q through `z^204`, given (U): (F_g) for every `g <= 69` (before: `g <= 9`).
  * `D_0 = 1` to ALL orders is equivalent, by steps PROVED given (U), GI Theorem G and standard formal lemmas, to one identity at
    `q = 0`, called (Gamma): the stationary-phase series of the Mellin-Barnes integral of the I-function, with all six Gamma factors
    replaced by their Stirling series, takes the value `h(z)` at `q = 0`.
  * The critic gives a complete formal outline of a proof of (Gamma): a change of chart in a 7-variable Euler-integral representation,
    the mechanism of Givental, arXiv:math/0108100 s.10. It is unrefereed.
* Nothing was refuted. No order failed in any field, frame or route; every FAIL line in the logs is a deliberate control.

## 2. The assembled statement after Stage 4

**Conjecture N(g) holds for every g >= 2, modulo exactly:**
* **(U)** Theorem U (`natural_form/theory/UNIFORM_REDUCTION.md` s.3): proved in the project, AI-audited only.
* **(Q)** the quoted theorems, unchanged from Stage 2: arXiv:1509.06954 Thm 8.1; arXiv:2212.00320 Thm 7.1; arXiv:math-ph/0702045
  Thm 5.1 with the linear loop equation of its s.4; arXiv:1303.1043 Def. 2.2, Prop. 2.4, Def. 2.5, Prop. 2.9; optionally
  arXiv:2502.18115 Cor. 1.2.
* **(QRR)** the `q = 0` normalisation `diag(h_(w_i))` of the actual R-matrix. It was implicit in Stage 2 ("the actual R-matrix"). It is
  listed now because (N) is exactly a statement about its `z^(5k)` part. Sources as read by the referees: arXiv:1702.06096 eq. (5)
  (definition of the formal quintic) and arXiv:math/0108100 Thm 9.1 (the Bernoulli form); statement numbers not re-checked by the
  critic or by me. The sign of the twist weight `-5 zeta^i` is pinned by (H) at `n = 1`.
* **(F_g)** `D_0 = 1 mod z^(3g-2)`, discharged in one of two ways:
  * for `g <= 69`, the finite identity: VERIFIED exactly over Q given (U), by two independent codes (evid worker, evid referee), and
    through `z^125` by the critic's third code; for `g <= 335`, VERIFIED mod `2^61-1` and `2^89-1`, given (U);
  * for all `g`: Theorem P(b) (PROVED given (U), GI Theorem G and four standard formal lemmas) together with (Gamma) (CONJECTURAL,
    complete formal outline, unrefereed; VERIFIED to `z^500` mod two primes and to `z^60` exactly).

Genus 2-5 remain VERIFIED directly by graph sums. Every other step is the Stage-2 chain (STAGE2_SYNTHESIS s.2), unchanged.
**If the chart argument survives refereeing, Conjecture N(g) holds for all g modulo (U), (Q), (QRR) and standard formal
stationary-phase lemmas.**

## 3. What Stage 4 established

* **What the model fixes (correction of my brief).** I had expected the model to fix the value of `ph_n` at `u = 0`. In fact
  `ph_n(0) = F_(0,n)` in every normalisation (PROVED three times), so the model fixes nothing there. The free coefficient at `n = 5m`
  is the TOP coefficient `f_m = [u^(4m)] ph_(5m)`, the value of the colour-0 WKB function at the orbifold point `q = oo`.
  `deg ph_n <= floor(4n/5)` with trivial resonances is PROVED twice (elementary; arXiv:0710.0889 Thm 4(ii)); equality is VERIFIED only.
* **Reduction chain** (each step PROVED as stated, AI-refereed).
  * Lemma P (colour reduction) and Lemma N (given U): `K_q := [u^(4q) zhat^(5q)] Ghat_11 = [z^(5q)] c`.
  * Connection form: `D = 1 <=> ph^can_n(1) = h_n` for all `n`. For `5 not| n` this is (H), with no free parameter; for `n = 5m` it is
    (N), with `f_m` fixed at the conifold by `K_m = 0`.
  * Theorem P(b): `Phi^can = Phi^SP`. Hence `D = 1 <=> (Gamma)`.
* **Theorem P(b), the mechanism** (poly s.6, with the poly referee's bookkeeping s.3).
  * The full-Stirling Mellin-Barnes 1-form, pulled back by the Chester-Friedman-Ursell map to the model phase `V/Y + Y^4/4`, reduces
    modulo exact forms to the model's five rows. Its stationary phase is `Phi^SP`.
  * Its residue at `x = 0`, the vanishing cycle, is `C0 (P^SP Y^C)_0`, which equals `C0 s(z) (Y^C_Q)_0` by Theorem G.
  * At order `u^1` both sides equal `-u/z`. So `s(z) = 1`: the full-Stirling normalisation IS the canonical one.
  * The residue also produces the conifold period: at `u^2`, `-1/(2z)` from the phase and `-1/(5z)` from the Stirling term `b_1 = 1/12`
    add up to `-7/(10z) = -gamma/z`.
* **(Gamma) in equivalent forms:** `c(z; L = 1) = h(z)`; `Phi^SP(z, u = 1) = h(z)`; `Phi^SP = h(z) Phi^ZZ`, with `Phi^ZZ` the
  coefficientwise expansion (5) of arXiv:0710.0889 (Zagier-Zinger).
* **The chart argument for (Gamma)** (critic, STAGE4_SYNTHESIS s.1.4; SKETCH, unrefereed).
  * Replace each Gamma factor of the summand by its Euler integral. This gives a formal 7-variable integral with phase
    `Phi = sigma log(-q) - y + 5 sigma log y + Sum_j (-x_j + (zeta^j - sigma) log x_j)`; the powers of `z` cancel because
    `Sum_j zeta^j = 0`.
  * Integrating `(x, y)` first and then `sigma` gives `c(z; L)` (for `L != 1`).
  * `Phi` is linear in `sigma`. Integrating `(sigma, x_0)` first therefore restricts to `x_0 = -q y^5/prod_(j>=1) x_j` and leaves a
    5-variable chart integral whose critical point and Hessian are analytic at `q = 0`.
  * At `q = 0` the chart integral is a product of Euler integrals for `Gamma(1 + 5/z)` and `Gamma((zeta^j - 1)/z)`, `j = 1..4`. Their
    Stirling factors multiply to exactly `h_(w_0)(z)`.
  * The poly referee located the mechanism in arXiv:math/0108100 s.10 (torus-equivariant CP^n, Thm 9.1, Lemma 10.4; an outline there).
    Here it is applied to the O(5)-twisted, colour-0 case.
* **Closed form of the free coefficients** (CONJECTURAL; VERIFIED for `m <= 100` at two primes, `m <= 25` exactly):
  `Sum_m f_m z^(5m) = exp(Sum_l (-1)^(l+1) B_(5l+1)(1/5) z^(5l)/(l(5l+1)))`, so `f_1 = B_6(1/5)/6 = 4537/3937500`.
* **The routes agree** (critic's own code): the three formulas for `[z^(5q)] c` are one affine functional of the free coefficients;
  the canonical, full-Stirling, QRR and closed-form free coefficients agree exactly.

## 4. My own checks (2026-11-14)

### 4.1 Re-runs of the agents' code

Run from a scratch copy of `stage4/` (the scripts write their logs next to themselves), `py -I -B`, with the stored arguments. Each new
log was compared line by line with the stored one, timings stripped.

| script and arguments | what it checks | result |
|---|---|---|
| `critic/cr4_run.py 12 25 60` | route agreement; `D_0 = 1 mod z^126` with `f^can` fixed at the conifold, no QRR input; (Gamma) exact through `w^60` at `u = 1/2, -1/3, 3/7`; controls | identical (196 s) |
| `critic/cr4_chart.py 5` | the chart integral of the outline equals `R_n(L)` through `z^5` at three `L`, and `h` at `L = 1`; `z^5` sees (N); controls | identical (197 s) |
| `evid/ev_direct.py 0 0 35 29` | E4 extended from the referee rows: `c^act = 1` through `z^35` (`a_1, a_3, a_5, a_7 = 0`) at two primes | identical (77 s) |
| `evid/ev_reduced.py Q can 40 0 0 12` | `D_0 = 1` exactly through `z^204`: `f^can` fixed by `K_q = 0`, then `ph^can_n(1) = h_n` | identical (152 s) |
| `evid/referee/rf_kq.py Q can 40` | the same, with the referee's independent code | identical (166 s) |
| `poly/referee/rf_gauge.py 30 10` | `c^act = 1` through `z^150`, two primes, two frames; Theorem-G identities; controls | identical (35 s) |
| `gamma/referee/rf_run_gamma.py 40 500 80` | (Gamma) as `c = R`: mod `2^127-1` and `2^107-1` through `z^500`, exact through `z^40`; controls | identical (870 s) |

Two of my first launches used default arguments that differ from the stored runs (`NR = 60` instead of 80; `QREG = 40` instead of 12).
I stopped those two processes and re-ran them with the stored arguments; only the re-runs with the stored arguments are reported.

### 4.2 My own check of (Gamma): `stage4/independent_check/gamma_interp.py` (logs in `stage4/independent_check/logs/`)

New code, stdlib only, sharing nothing with the agents' code. It uses no WKB recursion, no gauge and no Theorem U: only the
definition of `c(z; L)` and the definition of `h`.
* For up to 28 rational values `L != 0, 1` it computes `c_0(L), ..., c_N(L)` exactly by formal stationary phase in `sigma`.
* For each `n` it fits the predicted shape `c_n(L) = Sum_(k <= 4n/5) a_k L^(4n-5k)` on the first points and confirms the fit on three
  further points. A pole at `L = 1` would break the fit.
* It evaluates the fit at `L = 1` (`q = 0`) and compares with `h_n = [z^n] h_(w_0)(z)`.

| run | result |
|---|---|
| `gamma_interp.py 10 3` | polynomial shape and `c_n(1) = h_n` for all `n <= 10` (0.8 s); the three controls below fail as predicted |
| `gamma_interp.py 20 3` | the same for all `n <= 20` (15 s). Controls: Stirling term of `Gamma(1+5s)` dropped fails at `z^1`; `b_1 + 1/7` fails at `z^1`; `b_3 + 1/7` (visible only from `z^5` on) fails first at `z^5` |
| `gamma_interp.py 30 3 noctrl` | the same for all `n <= 30` (108 s); every fit confirmed on at least three further points |

So (Gamma), including (N) at `z^5, z^10, ..., z^30`, holds through `z^30` by a fourth, method-independent code. This is evidence, not
proof; the proof would be Stage 5.
*My own bug, found and fixed.* The first `z^30` run reported FAIL for `n >= 25`. The script had only 20 sample points: from `n = 25`
the fit had more unknowns than points (and from `n = 24` no point was left to check the fit). With a longer list of points and
assertions that every system is square and leaves three check points, `n <= 30` passes. The `z^10` and `z^20` runs needed at most 20
points and were not affected; all three logs are from the final script.

### 4.3 My reading of the chart argument

I went through STAGE4_SYNTHESIS s.1.4 by hand. Checked:
* The `z`-powers of the six Euler substitutions multiply to `z^-1`, because `Sum_j zeta^j = 0`.
* The 7-variable critical point is `(sigma, y, x_j) = (L, 5L, zeta^j - L)`: `G = 0` reads `-q (5L)^5 = prod_j (zeta^j - L) = 1 - L^5`,
  i.e. `u L^5 = 1`. The partial Hessians of both orders of integration are non-degenerate for `L != 0, 1`: diagonal entries
  `-1/(zeta^j - L)` and `-1/(5L)` in (a); determinant `-1/x_0^2` in `(sigma, x_0)` in (b).
* The formal delta lemma (L2). In the coordinates `(sigma, G)`, with `G = log(-q y^5/prod_j x_j)` and `dx_0/x_0 = -dG`, the phase is
  `sigma G + f(G; y, x_(j>=1))` and the amplitude does not involve `sigma`. The Gaussian propagator has `<G G> = 0`, and every vertex
  carries only `G`-legs, so every Feynman graph with a vertex vanishes. The formal `(sigma, x_0)` integral is therefore exactly the
  restriction to `G = 0` times its Gaussian factor, to all orders.
* At `q = 0` the twist and tangent weights match: since `1 - 2k` is odd,
  `Sum_(a in w_0) a^(1-2k) = -(5^(1-2k) + Sum_(j>=1) (zeta^j - 1)^(1-2k))`. So the Stirling factors of `Gamma(1 + 5/z)` and of
  `Gamma((zeta^j - 1)/z)` give exactly `h_(w_0)`, with `B_2k/(2k(2k-1))` as in QRR.
* The continuation to `L = 1` needs no a priori polynomiality. At each order `c_n(L)` is a rational function of `L`, equal on a
  punctured neighbourhood of `L = 1` to an analytic function (the chart series). It is therefore regular at `L = 1`, with value `h_n`.
Not checked by me: the full bookkeeping of the powers of `z` and of the Gaussian prefactors (the critic's code tests it through `z^5`),
and Theorem P(b) beyond reading its proof and the poly referee's s.3.

## 5. What is left

1. **Stage 5: write out and referee the chart argument.** State and prove, or cite with statement numbers, (L1) multivariable formal
   stationary phase (changes of variables, Fubini, parameters), (L2) the formal delta lemma, (L3) Stirling's series as the formal
   stationary phase of Euler's integral with `x^m` inserted. Write steps (a)-(e) with the `z`-power and prefactor bookkeeping. Optional:
   `critic/cr4_chart.py 6` (about 30 min). With Theorem P(b) this gives `D_0 = 1` to all orders, and Conjecture N(g) for all g modulo
   (U), (Q), (QRR) and standard formal lemmas.
2. **Human checks**, in this order: Theorem P(b) with the poly referee's bookkeeping (about two pages); the chart argument once
   written; (QRR) and its `z^(5k)` part; Lemma N (now needed only for the finite range `g <= 69`); then the Stage-2 list, starting with
   Theorem U. A cheap independent test of (QRR) at `z^5` (critic): the genus-3 degree-0 constant-map number from the project's R-matrix
   against the Hodge-integral formula; genus 2 is blind to a `z^5` change.
3. **Write-up**: a section "The connection problem" (10-15 pp) in the planned separate note `allgenus/note/`.
4. **Housekeeping, not done**: the critic's corrections O1-O18 (STAGE4_SYNTHESIS_2 s.5) and the Stage-2 corrections O1-O25 to the
   worker notes. The `osc` worker was not re-run; its natural brief is now Stage 5.
