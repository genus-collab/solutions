# Referee report: Stage 4, task "evid" (falsifier for D = 1), 2026-11-14

Adversarial referee of `allgenus/stage4/evid/` (EVID.md, ev_*.py, 24 logs). AI-written, AI-checked; no human has read it.
Labels: PROVED (argument written here or in the cited note), VERIFIED (exact computation, range and field stated), CONJECTURAL,
SKETCH, REFUTED. Conventions of EVID.md, `natural_form/theory/UNIFORM_REDUCTION.md` (UR) and `allgenus/gauge/GAUGE_IDENTITY.md`
(gauge). p89 = 2^89-1, p61 = 2^61-1. My code: `rf_*.py` in this folder. It imports NO worker module and NO project module
(stdlib only). Logs: `logs/rf_*.log`.

## 0. Verdict: SOUND WITH CORRECTIONS

* **Nothing is REFUTED.** Every VERIFIED claim of the worker was re-derived with independent code, at two primes and/or exactly,
  with failing controls. All of them reproduce, digit for digit where numbers are printed: f_5 = 4537/3937500, f_10, f_15, the
  prefixes of f_20..f_200, m5, and K_19 = 5^-95/9120.
* **Every PROVED claim checks out given its stated inputs.** These are Theorem U (UR s.3, AI-audited only), (R1), formal WKB
  uniqueness, the model's symplecticity and Prop. S, covariance, and the q = 0 (QRR) convention (s.2).
* **Corrections (s.7).**
  * One labelling overclaim: "(F_g) VERIFIED for g <= 69 / 169" needs "given Theorem U".
  * One claim is not independent: the second and third frames are no independent test of (N), because by s.3 they cannot fail.
  * One small error in the exponent interval of s.2.1(b').
* **Two referee results go beyond the worker.**
  * **A (s.3).** The conifold periods C0, C1 DROP OUT of the normalisation functional K_q. (N) becomes a statement about the PF
    operator, the QRR constants, tau and the model row F_2 only.
  * **B (s.4).** deg ph_n <= floor(4n/5), and the triviality of every resonance, are PROVED for all n by an elementary argument at
    the orbifold point x = 1/s = 0. The worker had them as VERIFIED and "unexplained". Theorem U is not needed for the resonances.
* **Highest orders reached with my code** (no failure anywhere):
  * Exact over Q, frames (0,0) and (1/3,-2/7):
    * K_q = 0 for q <= 40, and (H) for n <= 204;
    * connection form ph^can_n(1) = h_n for all n <= 204 (frame (0,0));
    * (H) alone to n = 300.
    * Hence D_0 = 1 mod z^205, given Theorem U (Lemma N).
  * Mod p89 AND p61, frame (0,0): K_q = 0 for q <= 200 and (H) for n <= 1004. So **D_0 = 1 mod z^1005, given Theorem U**:
    (F_g) for g <= 335, a_1 = a_3 = ... = a_199 = 0. This goes beyond the worker's z^505.
  * **Without Lemma N, Lemma P or K_q:**
    * the actual gauge (five colours, matrix inverse) equals my own canonical gauge (own Theorem-U solver from the periods) in every
      slot through z^100 (u < 85) at p89 and p61. With the connection route D9 this gives D_0 = 1 mod z^101 (g <= 34). The worker
      reached z^55, p89 only.
    * Exactly over Q, the connection route gives the same through z^40.

## 1. Independent machinery (what is shared with the worker, and what is not)

| piece | worker (ev_core / project) | referee (rf_core, rf_gauge) |
|---|---|---|
| WKB polynomials ph_n | (HPF) of ref_rows, A_e operators in (zhat, u) | the 5th-order PF operator `M^5 - 1 - q prod(5M + kz)` applied literally to Laurent monomials `z^n s^e`; transport operator T_1 found and asserted |
| periods C0, C1 | dual numbers u^(n+eps) | Frobenius on `u^m (log u)^l`; the vanishing period FOUND as the analytic solution with a log partner (3x3 obstruction matrix, rank 2): gamma = 7/10 |
| QRR h_n | Newton sums for the quartic in y | Galois trace of (1-zeta)^-m in Q(zeta5) (exact), direct weights in F_p (all colours); Bernoulli numbers by the table algorithm |
| model rows F_b | affine trick on the closure | closure operator `-T_0 T_1 T_2 T_-2 T_-1` on monomials, pivot 4(n-1) zeta^(4i) asserted; colour i solved separately |
| actual gauge | colour reduction (Lemma P) / gi_actual_norm | five colours mod p, `R^` inverted as a 5x5 matrix series (no symplecticity, no covariance used), plus Lemma P exactly |
| canonical gauge | ur_gauge (project) | own Laurent-series solver of (E_k,b) with own ODE solver; regularity ASSERTED at every order; canonical and stored normalisations |

Checks of the shared conventions:
* The PF operator: re-derived (my derivation) from the I-function `Sum_d q^d prod_(k<=5d)(5H+kz) / prod_(k<=d)((H+kz)^5 - 1)`.
  * The I-function follows from the definition in arXiv:1702.06096 eq. (5): torus weights lambda_i = zeta^i lambda_0 and the twist
    e(R pi_*(S^-5)).
  * The pages I looked up (pp. 5 and 9) give eq. (5), L = (1 - 5^5 q)^(-1/5), C_0 = I_0 and C_1 = D(I_1/I_0). They do not print
    the I-function itself.
* The hatted rows (kh factors): checked by the CLOSING RELATION `Qh_0 = (1/C0) A_0 Qh_4` (D1). This is the PF equation as the step
  4 -> 0 of the cyclic system `c = (C1, 1/(u C0^2 C1^2), C1, C0, C0)`. It holds for all colours, n <= 100, u < 85, at two primes.
* The gauge comparison D6 checks the same conventions.
* Spot values reproduced:
  * my stored gauge gives P^s_11(0,z) = 1 - (43397/807975000) z^5 - (765368951400317/1305647201250000000) z^10 (= gauge s.1);
  * kappa residues (0, 1/5, -3/5, -2/5, -1/5);
  * p_(5,b)(0) = (43397/807975000, 0, 0, 283/2835000, -12659/201993750) (= UR s.4);
  * F_0 through t^7 (= UR s.2);
  * t_c = u + (7/10)u^2 + (41/75)u^3 (= gauge s.5).

## 2. The PROVED claims, line by line

* **s.2.1(a), ph_n(0) = F_(0,n) in every normalisation: CORRECT.**
  * The u^0 part of (HPF) closes on u^0 coefficients: theta_u kills u^0, and the bracket terms carry u^k, k >= 1.
  * The free coefficient sits at u^(4n/5), never at u^0.
  * Missing line: why the u^0 recursion is the model's. Fill it as follows. The scaling limit of the PF operator is
    `v th^5 + 3z th^4 + 1` (UR s.1). Its homogeneous WKB solution with leading coefficient 1 is unique (pivot 4n). It equals F_0
    (Prop. S, model side).
  * VERIFIED (mine): n <= 1000 at both primes, n <= 300 exact.
* **s.2.1(b), polynomiality: CORRECT.** Every non-identity word of `A_-5..A_-1` carries a LEFT factor (1-u). This holds for every
  colour, since zeta^(5i) = 1.
  * The resonance argument "via Theorem U" is also correct. It is superseded by s.4, which needs no Theorem U.
* **Lemma P: CORRECT.** Re-derived:
  `(Qhat R^-1)_ac [zhat^N] = Sum_n Qh^(0)_(a,n) (-tau)^(N-n) F_(3-c,N-n) (1/5) Sum_i zeta^(i(a+5-c-N))`.
  * It uses R^(t)^-1 = (1/5) R^(-t)^T eta'. The check: pinv^T pinv = eta'/5 gives pinv = Psi0^T eta'/5.
  * VERIFIED independently (D3, D4): model symplecticity to t^100; Lemma P = matrix inverse, z^30, u < 85, two primes, two frames.
* **Lemma N: CORRECT given Theorem U.**
  * Step 4 is formal WKB uniqueness in formal series in z with coefficients in Q((s)). It applies because the amplitude of
    P_can Q^m has such coefficients.
  * Step 5 is right (a ratio that is a power series in u at fixed zhat forces d in Q[[z^5]]).
  * "c the same for all colours" holds because d(zeta^-i z) = d(z) for d in Q[[z^5]].
  * K_q = [z^5q] c and the (2,2) slot follow from (R1).
  * Slope 1: VERIFIED (mine) at every q <= 40 exactly. In can mode K vanishes after the shift f -> f - K.
* **s.2.4 (reduced criterion), Parts 1 and 2: CORRECT.**
  * The **connection form** holds even WITHOUT assuming (H): if ph^can_n(1) = h_n for all n, then Phi^can satisfies the QRR
    condition. So Phi^can is the QRR solution, (H) holds and D = 1. Conversely D = 1 gives Phi^QRR = Phi^can.
  * The colour-i statements follow from h_(w_i)(z) = h_(w_0)(zeta^-i z), since w_i = zeta^i w_0. This is implicit in the note.
    VERIFIED: (H) for all five colours, D-runs.
* **s.8.2: CORRECT.**
  * p_(5k) = 0 for odd k: y_(-j) = 1 - y_j = -zeta^j y_j, so y_(-j)^5 = -y_j^5.
  * [z^5k] log h_(w_0) = B_(5k+1)/(5k(5k+1) 5^(5k)); the sign was re-derived. The a_k formula holds given (H).
  * VERIFIED: p_m exactly for m <= 300 (zero exactly at the odd multiples of 5, non-zero at every other odd m). The B1 control below
    depends on p_95 = 0.
* **REFUTED item (the naive reading): CORRECT as stated.**
  * At fixed z the model's homogeneous WKB function has u^0 part 1, while f_5 = 4537/3937500 != 0.
  * It refutes a reading of the brief, not D = 1. Suggested wording: "the model alone does not fix f_5m; the conifold condition
    involves tau (s.3)".

## 3. Referee result A: the conifold periods drop out of K_q (PROVED given Lemma N and (R1); VERIFIED)

**Claim.** Take any polynomial-type colour-0 WKB function Phi (any free coefficients). Put
`S_q(u) := Sum_(n<=5q) (A_-1 Phi)_n (-tau)^(5q-n) F_(2,5q-n)`, with `(A_-1 Phi)_n = ph_n + (1-u)((4n-3)/5 - th_u) ph_(n-1)`.
Then S_q = O(u^4q) and **K_q = -[u^4q] S_q**.

*Proof.*
1. `Qh_1 = (1/(C0 C1)) A_-1 Phi + (1/C1) zhat (1-u) (th C0/C0^2) Phi`. Use `th(1/C0) = -th C0/C0^2`.
2. The second part pairs to `(1/C1)(1-u)(th C0/C0) Ghat_01[zhat^(5q-1)]`. This is Lemma P with (a, c) = (0, 1).
3. By Lemma N, `Ghat_01[zhat^(5q-1)] = tau^(1/2) u^(4q-1) p_(5q-1,1)(u)`.
4. By (R1), `p_(5q-1,1)(0) = [z^(5q-1)] (P_can c)_01(0, z) = 0`, since column 1 of P(0, z) is c(z) e_1.
5. Also `th C0 = O(u)`. So the second part is O(u^(4q+1)).
6. REG (Lemma N) then forces S_q/(C0 C1) = O(u^4q). In every conifold frame C0(0) C1(0) = 1 (-1) = -1. QED

**VERIFIED** (`rf_pfree.py`):
* -[u^4q] S_q = K_q (rows WITH C0, C1) in frames (0,0) and (1/3,-2/7). Range: exact for q <= 20, p89 for q <= 60.
* Both the QRR coefficients (all values 0) and a perturbation f_5q -> f_5q + 1/(q+2) (all values non-zero) were tested.
* REG' holds throughout.

**Consequences.**
* (N) reads [u^4q] Sum_(n<=5q) (A_-1 Phi^QRR)_n (-tau)^(5q-n) F_(2,5q-n) = 0 for all q. The ingredients are the PF operator,
  the Bernoulli/QRR constants, the canonical coordinate tau and the model. EVID s.0 and s.2.3 list "the periods C0, C1" among the
  ingredients: correct, but not minimal.
* **The runs in frames (1/3,-2/7) and (7/11,-3/5) are NOT independent evidence for (N).** Given Theorem U, K_q is identically
  frame independent. Those runs test code paths only. The same holds for "f^can = f^QRR in both frames".
* Exploratory control (scratch only, not logged): replacing the vanishing period g by g + f_2/7, which is NOT a conifold frame,
  still gives K_q = 0. The reason is that K_q sees C1 only through C1(0).

## 4. Referee result B: deg ph_n <= floor(4n/5) and trivial resonances, for all n (PROVED; VERIFIED)

Write x = 1/s (= L; u = oo is the orbifold point q = oo). In x the operator is `M = zeta^i x - (z/5)(1 - x^5) th_x`, and
`(5M + kz) x^m = 5 zeta^i x^(m+1) + z (k - m) x^m + z m x^(m+5)`.

*Lemma.* `Pi := (1/3125) prod_(k=1..5) (5M + kz)` maps x^m (m >= 1) into x^6 Q[x][z].

*Proof.*
1. The five factors are polynomials in M, so they commute.
2. Apply k = m, m+1, ..., 5 first. Each kills the lowest power exactly when that power equals k, so the x-valuation rises to >= 6.
3. No factor lowers the valuation.

*Consequence.*
* `P = (x^-5 - 1) Pi + M^5 - 1` maps x Q[x][z] into itself, componentwise in z.
* T_1(x^m) = zeta^-i (1-m)(x^(m-1) - x^(m+4)). Induction on n then gives psi_n in x Q[x], i.e. s-exponents <= -1. Take the lowest
  x-power m_0 <= 0 of psi_n: its image term `(1 - m_0) c x^(m_0 - 1)` cannot be cancelled and is not in x Q[x], a contradiction.
* This is exactly deg ph_n <= floor(4n/5), for every colour and every choice of the free coefficients.
* The resonance is trivial, because the right-hand side has top s-exponent <= -1, so W stops at s^-6 and never reaches s^-5.
* So psi_n (the z^n coefficient of psi) is a polynomial in L = x, divisible by L, of degree <= 4n+1. The amplitude
  Phi = psi/L lies in Q[L][[z]], and its value at the orbifold point L = 0 is `1 + Sum_m f_5m z^(5m)`, which lies in Q[[z^5]].
  This makes EVID s.8.3 ("three points") rigorous.

**VERIFIED** (`rf_struct.py` T6, T7):
* T_j(s^e) has only s-exponents <= -1 for all -300 <= e <= -1 and j = 1..5.
* The valuation of Pi(x^m) is 6 for m = 1..5.
* Equality deg ph_n = floor(4n/5) for n >= 3 (non-vanishing of the top coefficients) remains VERIFIED only: n <= 1004 at two primes,
  n <= 300 exact.

## 5. VERIFIED claims re-derived (all logs end in "ALL CHECKS PASSED")

| worker claim | worker range | referee re-derivation (independent code) | log |
|---|---|---|---|
| (H): ph_n(1) = h_n, 5 not dividing n | exact n <= 300; p89, p61 n <= 1000 | same, and n <= 1004 inside the K_q runs; all five colours to n <= 100 (mod p) | rf_struct_*, rf_kq_* |
| K_q = 0 (qrr), frame (0,0) | exact q <= 40; mod p q <= 100 | exact q <= 40 (144 s); **p89, p61 q <= 200** (z^1000) | rf_kq_Q_qrr_0_0_200, rf_kq_{89,61}_qrr_0_0_1000 |
| frame (1/3,-2/7), qrr | exact q <= 40 | exact q <= 40 (not independent for (N), s.3) | rf_kq_Q_qrr_1o3_m2o7_200 |
| frame (7/11,-3/5) | p89 q <= 40 | p89 q <= 40 | rf_kq_89_qrr_7o11_m3o5_200 |
| can mode: ph^can_5q(1) = h_5q, f^can = f^QRR | exact q <= 40 | exact q <= 40; slope 1 at every q | rf_kq_Q_can_0_0_200 |
| K22 = K11 and REG | exact q <= 12, mod p q <= 40 | exact q <= 12; mod p K22 q <= 100, REG q <= 40 | same |
| ph_n(0) = F_(0,n); deg ph_n | n <= 1000 (mod p), n <= 300 (Q) | same ranges; upper bound now PROVED (s.4) | rf_struct_* |
| direct E4: actual = canonical in every slot | z^35, u < 29 (2 primes, 2 frames); z^55, u < 45 (p89) | z^35, u < 36 (2 primes, 2 frames); z^55, u < 50 (p89, p61 frame (0,0); p89 frame 2); **z^100, u < 85 (p89, p61; p89 frame 2)** | rf_direct_* |
| Lemma P = matrix inverse (V6) | z^35, u < 29 | z^30, u < 85, two primes, both frames (D4) | rf_direct_* |
| connection form, Lemma-N-free | - | exact N <= 40 (frame (0,0)), N <= 30 (frame 2); mod p N <= 100 (D9) | rf_conn_*, rf_direct_* |

## 6. Failing controls (each fails its test as predicted)

* **S.** m5 is taken from MY stored gauge, exactly (-43397/807975000). f_5 -> f_5 + m5 gives K11_1 = K22_1 = m5 exactly, and REG
  stays blind. Exact and p89.
* **C4.** The stored gauge violates the connection form by exactly m5 at N = 5.
* **D8.** The stored gauge differs from the actual gauge first at zhat^5, in 3017 to 20408 slots depending on the window.
* **T.** Setting tau = 1 breaks REG.
* **W.** Untwisted weights break (H) at n = 1.
* **B1.** B_96 + 1 leaves (H) and K_q (q <= 18) intact and gives **K_19 = 5^-95/9120 exactly**.
* **B2.** B_8 + 1 breaks (H) first at n = 7.
* **QRR sign test.** At n = 1, h_1 = -(1/12)(Sum_j (1-zeta^j)^-1 + w^-1) = -(1/12)(2 + w^-1) must equal ph_1 = -3/20. This forces
  w = -5. The alternatives +5, no twist, and the reversed tangent orientation (zeta^j - zeta^i) all fail at n = 1.

## 7. Overclaims, gaps, errors (quote, then correction)

1. **EVID s.0, s.4 and JSON:** "(F_g) is VERIFIED (exact) for every genus g <= 69" (and "g <= 169" mod p).
   * The computation (K_q = 0, h_auto) is exact. Its translation into D = 1 mod z^205 goes through Lemma N, i.e. Theorem U is used
     as an existence theorem.
   * Should read "VERIFIED given Theorem U". No new input is involved, since the conditional theorem already assumes (U).
   * Lemma-N-free coverage:
     * worker: z^35 (two primes), z^55 (p89);
     * referee: z^100 (two primes), z^40 exact.
2. **EVID s.4 / s.5 / JSON:** the second and third frames are offered as further confirmation of (N). By s.3 they cannot fail
   given Theorem U. Relabel them as code and convention checks.
3. **EVID s.2.1(b'), s.6 and s.8.3:** "deg ph_n = floor(4n/5) ... VERIFIED, not proved ... unexplained; a naive weight count gives
   only floor((4n+4)/5)". The upper bound and "deg rho_n < 4n/5 for 5 | n" are PROVED (s.4). Only the exactness of the degree stays
   VERIFIED.
4. **EVID s.2.1(b'):** "exponents in (-4n/5, 0]". The left end -4n/5 is attained, because ph_n(0) = F_(0,n) != 0. Read [-4n/5, 0].
   The same typo appears in the ev_struct.py docstring.
5. **EVID s.2.1(a):** the identification of the u^0 recursion with the model is asserted rather than argued. One line fixes it (s.2).
6. **QRR convention** (the worker's open issue). Looked up: arXiv:1702.06096 eq. (5) defines the formal quintic by e(R pi_*(S^-5))
   at lambda_i = zeta^i lambda_0. With the usual quantum Riemann-Roch the twist enters with the Euler sign, i.e. as the weight
   -5 zeta^i in the tangent-type formula. That is the project's w_i, and (H) pins the weight at every order with 5 not dividing n.
   * NOT re-derived: that the q = 0 normalisation carries no further factor in Q[[z^5]]. This matters, because such a factor is
     exactly what (N) measures.
7. **EVID s.8.1** "No analytic connection problem in hatted variables" (SKETCH): there is still a connection problem. By s.4 it is
   polynomial in L at every order and links three points: u = 0, u = 1 and u = oo.

No numerical error was found. Every number I checked agrees: f_5..f_200 (prefixes), m5, K_19, slot counts of the same order, and
the 259 s / 3 min timings.

## 8. What survives

* All worker results and the reduction (Lemma P, Lemma N, s.2.4 and the connection form) stand, with the corrections of s.7.
* **D = 1 is CONJECTURAL.** The evidence after this referee pass:
  * Given Theorem U: exact through z^204; mod p89 and p61 through z^1004.
  * Free of Lemma N: mod two primes through z^100; exact through z^40.
* By-product for the input (U) itself: my solver ASSERTS regularity of every p_(k,b) at every order. Ranges: k <= 100 at p89 and
  p61 (frame (0,0)) and at p89 (frame 2); k <= 40 exactly (frame (0,0)); k <= 30 exactly (frame 2). This is computational support,
  not a proof.
* The cleanest current statement of the open problem: with tau = 1/2F1(1,4/5;9/5;u) and the model row F_2, the free (top = orbifold)
  coefficients fixed by [u^4q] Sum_n (A_-1 Phi)_n (-tau)^(5q-n) F_(2,5q-n) = 0 make Phi(z, u = 1) = h_(w_0)(z).
  * At the independent orders z^5k (k odd) this says [z^5k] log Phi^can(z,1) = B_(5k+1)/(5k(5k+1) 5^(5k)).
  * SKETCH for a proof: an orbifold-point evaluation of Phi(z, L = 0), which lies in Q[[z^5]] by s.4, connected to both ends.

## 9. Files, runs, compliance

* Code (all < 250 lines):
  * `rf_core.py`: fields, periods, tau, QRR, model, PF-direct WKB, A_e;
  * `rf_gauge.py`: Theorem-U solver;
  * `rf_kq.py`: reduced test, modes, controls;
  * `rf_direct.py`: five colours, gauge comparison, D1-D9;
  * `rf_conn.py`: exact connection route;
  * `rf_struct.py`: (H), degrees, p_m, orbifold lemma;
  * `rf_pfree.py`: s.3.
* Runs: `py -I -B`, stdlib only, working directory in the session scratch folder, at most 4 concurrent processes, all finished.
  * Longest: about 3 min per prime for z^1000; 174 s exact for the can mode.
  * No existing file was modified. No `__pycache__` was written. Logs print relative paths only.
* Looked up: arXiv:1702.06096 (abstract page and PDF pp. 1-5, 7-10: eq. (5), C_0 = I_0, C_1 = D(I_1/I_0),
  L = (1 - 5^5 q)^(-1/5)). The q = 0 R-matrix normalisation itself is not stated on those pages.
