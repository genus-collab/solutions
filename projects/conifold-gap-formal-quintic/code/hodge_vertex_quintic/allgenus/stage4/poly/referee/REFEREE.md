# Referee report: Stage 4 / poly, "polynomial connection: elementary route" (2026-11-14)

Adversarial referee of `allgenus/stage4/poly/` (POLY.md, `pp_*.py`, `logs/`). AI-written and AI-checked; no human has read this report.
Labels as in POLY.md: PROVED, VERIFIED, CONJECTURAL, SKETCH, REFUTED. Referee code `referee/rf_*.py` (stdlib, `py -I -B`, run from the
session scratch folder). No project module and no `pp_*` module is imported. The worker's logs are read as data only, to compare numbers.
Referee logs: `referee/logs/`. Primes p61 = 2^61-1, p89 = 2^89-1, both 1 mod 5, so zeta_5 lies in F_p. Conifold frames (0,0), (1/3,-2/7).

## 0. Verdict: sound, with corrections

* **No mathematical error found.** Every VERIFIED claim was reproduced by independent code, by a different method where one exists,
  with failing controls. Each claim was also pushed beyond the worker's range:
  * `D_0 == 1 mod z^151` (worker: `z^71`), mod two primes, in **both** frames.
  * (A2) directly through `w^50` (worker: `w^25`).
  * the residue identity through `u^30` in both frames (worker: `u^14`, one frame).
* **The reduction holds as stated.** `D_0 = 1 <=> (A2)` follows from Theorem P(b), given Theorem U, GI Theorem G and the standard
  formal lemmas the worker lists. The sign and normalisation bookkeeping that POLY.md leaves open is written in s.3 below: the
  conversion factor is `kappa(u) = C0` exactly, in every frame.
* **The corrections (s.4)** concern labels, independence counts, one wording overclaim and missing literature.
  * Lemma 3 can be upgraded to PROVED by quoting arXiv:0710.0889 Thm 4(ii).
  * (A2) is not a quotable theorem. The mechanism behind it is the q = 0 computation of arXiv:math/0108100 s.10 (an outline for
    torus-equivariant CP^n). The worker's route (ii) is exactly that chart argument.
* **Status.** (A2), hence `D_0 = 1` to all orders, stays CONJECTURAL.

## 1. Independent re-derivations (all PASS; every control fails as it should)

| claim (POLY) | referee method | range | result |
|---|---|---|---|
| (H) for `5 not| n`; `deg ph_n <= floor(4n/5)`; `ph_n(0)` = model; `f_m`; `ph_5` (s.1-s.2, s.4) | `rf_qexp.py`: the QRR WKB computed at LARGE RADIUS as q-series. `A_n(0) = h_n` is imposed for EVERY n, and HPF is not used. Transport operator `T(A) = 5 L^4 [(1-3125q) D A - 625 q A]` (derived, then checked numerically). Then `ph_n = (1-3125q)^((4n+1)/5) A_n`. | `n <= 150`, p61 and p89 (`rf_qexp_N70.log`, `_N150.log`) | Every q-series terminates at degree `floor(4n/5)`; this is (H) for all 120 values `5 not| n`. Degree is `= floor(4n/5)` except `n = 2`. `ph_n(0) = [t^n] F_0`. `f_1..f_14` equal the worker's exact values mod p. `ph_5` equals POLY s.1. Controls: shifting `A_3(0)`, `A_5(0)`, `A_10(0)` first breaks termination at `n = 3, 6, 11`, as predicted. |
| `c^act = 1` (s.3-s.4) | `rf_gauge.py`, three routes. (F) the worker's formula, re-implemented. (R) rows 1-4 rebuilt from row 0 by `P_(a+1,b) = (z P_ab' + bt_(b-1) P_(a,b-1))/b_a`, then closing of the 5-cycle, (R1), (R2), and `P_11(0,z)` read off directly. (G) row 0 of GI Theorem G, `Sum_b P_0b Y^C_b(v(u)) = (Y^C_Q)_0`, coefficientwise. `Y^C`, `Y^C_Q` are own code (`rf_vanish.py`). | `z^0..z^150`, p61 and p89, both frames (`rf_gauge_Q14_NU12.log`, `_Q30_NU10.log`) | (F): `c^act = 1`. (R): the cycle closes, (R1) and (R2) hold, and (R) agrees with (F). (G): 182 identities at Q=14 (window `u^0..u^12`, `z^-11..z^54`) and 330 at Q=30 (`u^0..u^10`, `z^-11..z^134`), counting only `z = 4 mod 5`; all other z-orders vanish on both sides by the grading. `p^act_(5q,0)(0)` equals the worker's logged values for `q <= 14`. Controls: the stored normalisation gives `[z^5] c = m5` and fails (G). No tau-weights and the wrong mu table (`4q + r`) each fail the rebuild (a required divisibility by u in a rebuilt row) and (G) (142-307 failures). `C0` of the other frame fails already at the regularity (divisibility) step of row 0. |
| (A2) directly (s.5, `pp_run_sp`) | `rf_sp.py`: formal SP by the Lagrange-Buermann (Morse lemma) formula `[y^2r] g(t(y)) t'(y) = [t^2r] g S^(-r-1/2)`, not by Wick contraction. The zeta-sums are done in F_p. | `w^0..w^50`, `u^0..u^42`, p61 and p89 (`rf_sp_N30.log`, `_N50.log`) | `Phi^SP` equals the QRR polynomials in every coefficient and is polynomial. Controls: no Stirling terms, and `b_k -> -b_k` (R versus R^-1), both fail from `w^1`. Model check (S2): `SP(Y^(c_b) e^(-f/w) dY) = eps_b F_b(w)`, `eps = (+,-,-,+,+)`, through `w^50`; the all-plus control fails. |
| residue identity (s.5) | `rf_res.py`, original variables. G is built from `Sum_j (1/x - zeta^j) log(1 - zeta^j x)` in F_p, not from its closed form. The Stirling sums are done in F_p. `Y^C_Q` is the referee's own Lemma-Q iteration. | `u^1..u^30`, both frames, p61 and p89 (`rf_res_NU30_K10.log`) | `V = C0 (Y^C_Q)_0`; V has no positive z-powers. Controls: dropping the Stirling terms, or using `gamma = 1/2`, fails at `u^2`. |
| CFU (s.5) | `rf_res.py`, exact over Q | `k < 10` | `Y_k in Q[T]`, `V = tau^(-5/4)`, `Y_1 = (1+T)/18`, and `Y(X=1)^5 = V` (critical point goes to critical point). |
| literature and closed forms | `rf_extra.py` (`rf_extra_M30.log`) | | `h_n` exact via Q(zeta_5) traces (`n <= 12`), `h_5 = -257071/8064000000`. arXiv:0710.0889 `Phi_0..Phi_3` (n = 5) equal `(1/h) x` (QRR WKB) as q-series. `F_inf = exp(Sum_l (-1)^(l+1) B_(5l+1)(1/5) z^(5l)/(l(5l+1)))` for `m <= 30`. |
| reproducibility | the worker's `pp_run_norm.py 7`, `pp_res.py 14` (scratch) | | Output identical to the stored logs, timing lines excluded. |

## 2. The PROVED claims, line by line

| POLY | claim | referee |
|---|---|---|
| s.1 | HPF and its `w^(n+1)` coefficient `(1-u)[(4n - 5 theta_u) ph_n + r_n]` | CONFIRMED by hand. `A_e` by the chain rule (`u d/du = theta_u - (4/5) theta_w` at fixed z); `z^m M^(5-m) psi = s^-1 u^(m-1) w^m Phi_(5-m)`; the coefficients `e_k(1/5..5/5)`; the `ph_n` part sums to `(4n + 3 - 5 theta) - 3`. Independently, the large-radius route (no HPF) gives the same polynomials. |
| s.2 L.1 | `ph_n(0)` = model for every normalisation (the naive criterion is vacuous) | CONFIRMED (also `rf_qexp` (M), `n <= 150`). |
| s.2 L.2 | power-series solutions are polynomial, equal to `Phi^(0) c(z^5)` | CONFIRMED. Multiplication by a function of z commutes with HPF because `((4/5) theta_w - theta_u) z^m = 0`. Theorem P(a) already supplies the power-series solution that the resonance argument needs, so Lemma 2 does not need Theorem U. |
| s.2 L.3 | `deg ph_n <= floor(4n/5)` | VERIFIED `n <= 150`. **Upgradable to PROVED (C1).** |
| s.2 L.4 | `F_inf^QRR = h/Phi^(0)(z,1)` | CONFIRMED. The log-coefficient remarks are also right: `z^10` has denominator `5^11`; `z^15` has primes `{2,3,5,17}`, as `B_16/240` does; `z^20` has primes `{2,5}`. Closed form in s.5. |
| s.3 | `c = P_00(0,z) + z d_u P_01(0,z)`, `[z^5q] c = [u^4q](g_5q + g_(5q-1))`, mu table | CONFIRMED. Row 0 of `z P' = A P - P At` gives `z P_01' = b_0 P_11 - bt_0 P_00`, with `b_0(0) = t_c'(0) = 1` and `bt_0(0) = v'(0) = 1`, and `P_01(0,z) = 0` by (R1). The formula is exactly the `u^1` coefficient of row 0 of Theorem G (`rf_gauge` (G)). See C7. |
| s.3 | `m5 = p'_(4,1)(0) = -p^can_(5,0)(0)` | CONFIRMED at both primes. |
| s.5 | `S(1/x) = -W(x;u)`; critical points `x^5 = u`; `W^(1) = (5/4)/tau`; `W^''(1) = 5/(1-u)`; amplitude `alpha = -(2pi)^-2 sqrt5 z at`, independent of u, `at(0,z) = 1` | CONFIRMED by hand. The `log z` terms cancel because `5 sigma = Sum_j (sigma - zeta^j)`. `sqrt(5 sigma/(sigma^5 - 1)) dsigma/z = -sqrt5 (1 - x^5)^(-1/2) dx/z`. The cyclotomic identity was also checked in F_p. |
| s.5 | CFU: `(T-1) Y_k = R_k(T) - V_k`, so `V = tau^(-5/4)` | CONFIRMED (exact, `k < 10`). |
| s.6 | Theorem P(a) | CONFIRMED modulo standard formal facts: formal SP commutes with `D_q`, kills total derivatives, and the formal Stirling series satisfies `St(y+1) = y St(y)`. `((z(s+1))^5 - 1) phi(s+1) = z^5 prod_k (5s+k) phi(s)` was checked by hand. |
| s.6 | Theorem P(b), `Phi^SP = Phi^can` | CONFIRMED as an outline; the referee's bookkeeping (s.3) closes the open step. Dependencies: Theorem U and GI Theorem G (AI-audited only), plus the standard lemmas of P(a) and of step 3. |
| s.6 | Corollary: `D = 1 <=> (A2)` | CONFIRMED given Theorem P. The QRR solution is the unique WKB column `Phi^can d(z)` with value h at u = 1. |

## 3. Theorem P, steps 4-6: the bookkeeping (referee; PROVED, VERIFIED below)

Notation: `f = V/Y + Y^4/4`, `V = tau^(-5/4)`, `Y_c = V^(1/5) = tau^(-1/4)`, `N_w = sqrt(2 pi w/5) exp(-f(Y_c)/w)`.
The model rows: `Y^C_b = eps_b I^C_(c_b)` with `eps = (+,-,-,+,+)` (UR s.1, GI Lemma V), and `F_b(0) = eps_b`.
* *SP of the model forms.* Substitute `Y = Y_c eta`; then `f(Y_c eta) = Y_c^4 (1/eta + eta^4/4)` and `w/Y_c^4 = tau w`. With (S2):
  `SP_(Y_c)(Y^m e^(-f/w) dY) = Y_c^(m+1) sqrt(2 pi tau w/5) e^(-f(Y_c)/w) eps_b F_b(tau w) = N_w eps_b tau^((1-c_b)/4) F_b(tau w)`
  for `m = c_b`. (S2) is VERIFIED through `w^50` and fixes the signs, which POLY left untracked.
* *Residue of the model forms.* Substitute `x = s Y`: `Res_(Y=0)(Y^m e^(-f/w) dY) = s^(-(m+1)) I^C_m(v(u), z)`.
* *Grading.* The reduction respects a Z/5 grading with `deg Y = 1` and `deg w = -1`, since `f(zeta Y) = f(Y)/zeta`. `B^ dX` has
  degree 1, so `pi_(c_b)` contains only `w^k` with `k = -b mod 5`: the grading of the hatted row-0 entry `G_b`.
* *Identification.* Write `Phi^SP/C0 = Sum_b G_b tau^((1-c_b)/4) F_b(tau w)`. The SP identity gives
  `Phi^SP = Sum_b eps_b pi_(c_b) tau^((1-c_b)/4) F_b(tau w)`, and the triangular expansion is unique, so `G_b = eps_b pi_(c_b)/C0`.
  Hatting with `c_0 = 0` gives `G_b = s^(c_b) P^SP_0b`.
* *Residue.* `V = Res_x(e^(-W/z) at dx) = s Res_X(e^(-W^/w) B^ dX) = s Sum_b pi_(c_b) s^(-(c_b+1)) eps_b Y^C_b
  = C0 Sum_b P^SP_0b Y^C_b = C0 (P^SP Y^C)_0`.
* **Conclusion: `kappa(u) = C0`, exactly, in every frame.** No power of s or tau survives. Step 7 then gives `s(z) = 1`. VERIFIED:
  `V = C0 (Y^C_Q)_0` through `u^30` in frames (0,0) and (1/3,-2/7) (`rf_res`).
* *Technical remark.* At each u-order the hatted integrand is `exp(-(1/X + X^4/4)/w)` times a Laurent polynomial in X with
  coefficients in `Q[w, 1/w]`. The residue functional is best defined as `(1/2 pi i)` times a contour integral around `X = 0` on C^*.
  Exact forms then integrate to zero, and the termwise-in-u change of variables of step 3 is immediate.
  The identity `F(X + uh)(1 + uh') = F + Sum_(j>=1) (u^j/j!) (F^(j-1) h^j)'` used there was checked by hand.

## 4. Corrections and overclaims

* **C1 (Lemma 3; label too weak).** arXiv:0710.0889 Thm 4(ii) (n = 5) gives `Phi_s in L Q[L]`. Their
  `F(w,x) = Sum_d x^d prod_(r=1..5d)(5w+r)/prod_(r=1..d)((w+r)^5 - w^5)` is the colour-0 series with `w = 1/z` and `x = q`. Their
  expansion (5) is coefficientwise in x, with `Phi_s(0) = delta_(s,0)`, so `Phi^ZZ = Phi^QRR/h(z)`; this is VERIFIED for `s <= 3`
  (`rf_extra` (E2)). Hence `A_n in L Q[L]`. Since `ph_n = u^((4n+1)/5) A_n` has no negative powers of u (every WKB column is `Phi^can d(z)`), `deg ph_n <= floor(4n/5)`:
  PROVED modulo the quoted theorem. The parallel Stage-4 "gamma" worker cites the same theorem.
* **C2 (independence counts).** POLY s.0 and s.4 count 14 identities for `c^act` and 56 for (H). Symplecticity (GI Prop. S) forces
  `c^act(z) c^act(-z) = 1` and `D_0(z) D_0(-z) = 1`. So only 7 of the c-identities are independent (`a_1, a_3, ..., a_13`) and 28 of
  the (H)-identities (odd n). The new independent content beyond E4 is `a_5..a_13` and (H) at odd n, `5 not| n`, `25 <= n <= 69`.
* **C3 (a "control" that is not one).** "c of `Phi^(0)` = `1/F_inf^QRR`" (`pp_run_norm`) is equivalent to the main result, because
  `c(Phi k(z)) = c(Phi) k(z)`. The genuine failing controls are `c^(0) != 1` and the stored normalisation. No check of the worker
  catches a convention error such as the wrong mu table, as POLY s.3 itself notes. The referee's closing and Theorem-G checks do
  (`rf_gauge` variants).
* **C4 (wording).** "(A2) involves only large-radius data" and "a statement at q = 0 alone" (s.0, s.6 Corollary, s.7) overstate.
  `Phi^SP` is produced as a power series around the conifold, `u = 0`; its coefficients happen to be polynomials. (A2) is the
  statement that their continuation to `u = 1` gives h, a connection statement. What is true is that the gauge, the periods and the
  model have disappeared from it.
* **C5 (labels in s.0 item 6 and s.7).** "(A2) VERIFIED ... via Theorem P ... through `w^70` (exact)" mixes proof and computation.
  Read instead: (A2) follows through `w^70` from Theorem P (PROVED mod U, G) and the VERIFIED `c^act = 1`.
  Referee: VERIFIED directly through `w^50`, and through `w^150` via Theorem P and `rf_gauge`.
* **C6 (Theorem P(b) label).** "PROVED modulo Theorem U and GI Theorem G" is acceptable now that s.3 is written down. The four
  standard formal lemmas should be listed explicitly in the statement: SP kills exact forms; the residue kills exact forms; both are
  invariant under formal changes of variables termwise in u; SP commutes with `D_q`. They are standard but not written in POLY.md.
* **C7 (s.3, "divisibility of `g_N` by `u^mu` is the regularity of P").** Divisibility is necessary but not sufficient. Regularity
  also needs the column-1 condition `p_(5q+4,1)(0) = 0` and regular rows 1-4. Both hold: they follow from Theorem U and (H), and
  `rf_gauge` (R) checks them through `z^150`.
* **C8 (missing literature).** The q = 0 mechanism of arXiv:math/0108100 s.10 is not cited; see s.5.
* **C9 (minor).** "(was `g <= 9`)" is out of date: the parallel `evid` worker reports `D_0 == 1 mod z^205` exact. This is not an
  error of this worker, and the referee has not checked that claim.

## 5. Is (A2) a known theorem?

Not as stated, as far as I could find. Searched: arXiv:0710.0889, arXiv:math/0108100 (read), and a web search on Givental's
oscillatory-integral R-matrices. No statement number applies to the twisted, one-dimensional, full-Stirling form.

* **The mechanism exists in the literature at outline level.** arXiv:math/0108100 s.10 covers torus-equivariant CP^n: Thm 10.1 (the
  mirror oscillating integral solves the equivariant QDE), Example 10.3 (SP), the paragraph before Lemma 10.4, Lemma 10.4 (Stirling),
  and "This coincides with `b_0(z)` in 9.1" (Thm 9.1: `R|_(Q=0) = exp diag(b_i)`, the Bernoulli/QRR normalisation). The argument:
  * in the chart `x_0 = Q/(x_1...x_n)` the phase stays analytic, with a nondegenerate critical point, up to Q = 0;
  * at Q = 0 the integral factorises into Gamma integrals, and Stirling gives `ln R^00 = Sum_k B_2k/(2k(2k-1)) z^(2k-1)
    Sum_(alpha != 0) (lambda_alpha - lambda_0)^(1-2k)`.

  The paper is an outline ("Details will be published elsewhere"), for the untwisted equivariant case, with a multidimensional
  integral. (A2) is its O(5)-twisted, one-dimensional analogue. POLY's route (ii), "keep the colour-0 Gamma exact", is this chart.
* **What a proof by this route needs** (SKETCH, the most promising route seen):
  * (a) the twisted mirror integral with the Laplace variable for `Gamma(1+5s)` (UR s.8, a SKETCH there). Its Stirling series is
    uniform near Q = 0, because its argument is at least `5/z`. Only the colour-0 factor `1/Gamma(1+d)` is non-uniform.
  * (b) formal Fubini and chart invariance, to identify the multidimensional SP with the full-Stirling one-dimensional SP for `Q != 0`.
  * (c) the Q = 0 factorisation, which yields exactly `h_(w_0)` (d = 0 prefactors: POLY s.5, re-derived by hand by the referee).
* **Equivalent forms found by the referee.**
  * (i) `(A2) <=> Phi^SP = h(z) Phi^ZZ`: the full-Stirling SP is h times the formal expansion (5) of arXiv:0710.0889.
  * (ii) Modulo the Laplace/Euler-Maclaurin step of POLY route (i), (A2) holds iff (5) is the true analytic asymptotic expansion of
    `F(w,x)` as `w -> +infinity`, at fixed `0 < x < 5^-5`. That analytic expansion is `e^(mu w) L Phi^SP/h`, and both are WKB columns.
    arXiv:0710.0889 proves (5) only coefficientwise in x (their Lemma 1). Status: SKETCH.
  * (iii) Closed form of the free constants (`rf_extra` (E3), VERIFIED `m <= 30`; first observed by the parallel "gamma" worker):
    `F_inf^QRR = exp(Sum_l (-1)^(l+1) B_(5l+1)(1/5) z^(5l)/(l(5l+1)))`, `f_1 = B_6(1/5)/6`. Combined with POLY's Lemma 4 form,
    (N) at order `z^(5m)` (`m <= 30`) holds iff `F_inf^can` equals this series.

## 6. What survives

| item | status after refereeing |
|---|---|
| HPF, recursion, Lemmas 1, 2, 4, the reduced criterion, the m5 explanation | PROVED (confirmed) |
| Lemma 3 (`deg ph_n <= floor(4n/5)`) | PROVED modulo arXiv:0710.0889 Thm 4(ii) (C1); VERIFIED `n <= 150` |
| phase, amplitude, CFU, Theorem P(a) | PROVED (confirmed; standard formal facts) |
| Theorem P(b) and the Corollary `D_0 = 1 <=> (A2)` | PROVED modulo Theorem U, GI Theorem G and four standard formal lemmas; bookkeeping written here (`kappa = C0`) |
| `D_0 == 1 mod z^151` (two primes, both frames), hence (F_g) for `g <= 51` | VERIFIED by the referee (independent code). Conditional, like E4, on Theorem U. The UR conventions are now tested by the closing and by Theorem G. The worker's `z^71` (exact, frame (0,0)) is CONFIRMED. |
| (A2) | VERIFIED directly `w <= 50`; CONJECTURAL in general. Mechanism: SKETCH (s.5). |
| `D_0 = 1` to all orders | CONJECTURAL, reduced to (A2) |

Not checked by the referee: Theorem U and GI Theorem G themselves (AI-audited dependencies); the 4th-order period operator of UR s.1,
taken as input but tested indirectly by closing, Theorem G and the residue identity; the parallel `evid` and `gamma` notes.

## 7. Files and runs (all `py -I -B`, working directory in the session scratch folder, at most one process at a time, none left running)

| file | content | log |
|---|---|---|
| `rf_core.py` | arithmetic mod p (Kronecker products), Bernoulli numbers, `h_(w_0)` (mod p and exact), periods, frames, model rows | (imported) |
| `rf_qexp.py` | the QRR WKB at large radius (q-expansion), hatted polynomials, (H), degrees, controls | `rf_qexp_N70.log`, `rf_qexp_N150.log` |
| `rf_vanish.py` | `Y^C` (closed form), `Y^C_Q` (Lemma Q) | (imported) |
| `rf_gauge.py` | row-0 extraction; (F), (R), (G); stored normalisation and convention controls | `rf_gauge_Q14_NU12.log` (4 s), `rf_gauge_Q30_NU10.log` (30 s) |
| `rf_sp.py` | Lagrange-Buermann SP; model rows with signs; controls | `rf_sp_N30.log` (9 s), `rf_sp_N50.log` (89 s) |
| `rf_res.py` | residue identity in both frames; CFU exact | `rf_res_NU30_K10.log` |
| `rf_extra.py` | exact `h_n`; arXiv:0710.0889 `Phi_0..Phi_3`; closed form of `F_inf` | `rf_extra_M30.log` |
