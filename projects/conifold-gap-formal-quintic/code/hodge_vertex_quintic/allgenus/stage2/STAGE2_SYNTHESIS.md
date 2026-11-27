# Stage 2 synthesis (critic and assembler): ML_K and Conjecture N(g) for all g, modulo exactly which inputs (2026-11-14)

Critic of the Stage-2 workflow `allgenus/stage2/`: four workers (`dict`, `logleg`, `flat`, `bvar`), four referees, each referee
"sound with corrections". AI-written; no human has read it. Labels: PROVED (argument written, or a quoted theorem with arXiv number
and statement number), VERIFIED (exact computation; range, field or primes stated), CONJECTURAL, SKETCH, REFUTED. Conventions of
`../SYNTHESIS_STAGE01.md` s.3: `H = d/dv`, `e = tau_1(1) + v tau_0(H)`, "(G)", "Gap Theorem", `K_vc` = the slice `K`
(`K_0 = K cap ker Pi_0`), "ML_K", "quintic Theorem R", `P_stored = P_can (1 + m5 z^5 + ..)`, `m5 = -43397/807975000`.
Part 2 (falsifiers, contradictions, overclaims, next steps, critic record): `STAGE2_SYNTHESIS_2.md`.

Read in full: `dict/DICT_PROOF.md`; `logleg/LOGLEG_PROOF.md`, `LOGLEG_PROOF_2.md`, `FALSIFIER_B.md`; `flat/S15_PROOF.md`; `bvar/BVAR.md`,
`BV_TABLES.md`; the four `*/referee/REFEREE.md` and `bvar/referee/REFEREE_S4.md`; Stage 0-1: `../MECHANISM.md`, `../SYNTHESIS_STAGE01.md`,
`../theory/REDUCTION_DRAFT.md`, `_2.md`, `../model/MODEL_SINGULAR.md`, `../tr/TR_MODEL.md`, `../xyswap/XYSWAP_MODEL.md` and its referee,
`../gauge/GAUGE_IDENTITY.md`, `natural_form/theory/UNIFORM_REDUCTION.md` (UR) s.2-s.4. All worker and referee logs were scanned: every
run ends in PASS; every FAIL line is a deliberate control. No worker folder is missing.
Critic code (new, `critic/`, `py -I -B`, stdlib): `cr2_legcore.py` + `cr2_legdict.py` (the leg dictionary LD, s.2 step 7.2),
`cr2_e2e.py` (LD on actual TR correlators), `cr2_count.py` (a disputed count). Logs `cr2_legdict_24.log`, `cr2_e2e.log`, `cr2_count.log`.

## 0. Bottom line

* **Nothing was REFUTED.** The three cheap falsifiers of SYNTHESIS_STAGE01 s.6 and the genus-3 s.1.5 test all PASSED (part 2, s.3).
* **ML_K is PROVED for every genus h >= 0 modulo quoted theorems only** (s.1). Quoted:
  * arXiv:1509.06954 Thm 8.1 (TR = the model CohFT);
  * arXiv:2212.00320 Thm 7.1 (L-regularity);
  * arXiv:math-ph/0702045 Thm 5.1 (variational formula), with the standard TR properties of its s.4;
  * for the constants `c_h`: arXiv:2502.18115 Cor. 1.2, or the Borel-sum route through arXiv:2212.00320 Thm 7.1 again.

  Every non-quoted step now has a written proof:
  * D(ii): dict referee s.3;
  * D(iii): DICT_PROOF;
  * (V): Theorem V of LOGLEG;
  * the **leg dictionary LD**: Stage 1 used it without proof. It is PROVED here (s.2, step 7.2) and VERIFIED by new critic code.
* **s.1.5 is no longer an input.** Theorem S15 (flat worker; referee CONFIRMED) derives the flat = hatted identification graph by graph
  from (F) = (H) + (N). The mechanism statement listed (F) already.
* **Assembled theorem (s.2).** Conjecture N(g) holds for every `g >= 2` modulo exactly the following:
  * (U) Theorem U;
  * (F_g) `D == 1 mod z^(3g-2)` in UR s.4's `Q D = P_can Q^m(v(u))`, i.e. (H) and `c^act = 1` through `z^(3g-3)`;
  * the quoted theorems above, plus arXiv:1303.1043 Def. 2.2, Prop. 2.4, Prop. 2.9 (R-matrix action).

  Every other step is PROVED, AI-written and AI-refereed; no human has checked any of it.
* **Genus 2-9.** (F_g) is the finite identity E4 (VERIFIED, not proved). So N(g) for `g <= 9` rests on (U), the quoted theorems and E4.
  Genus 2-5 are also VERIFIED directly. The only open mathematical input beyond quoted theorems and Theorem U is (F) for `g >= 10`:
  the connection problem, Stage 4. Its next independent coefficient is `a_5` at `z^25`.
* **B is a coordinate** (Lemma C of the bvar referee, PROVED). Falsifier (c) therefore tests `(F, Gm)`-universality on pulled-back
  systems plus the covariance of the pipeline, not a hidden use of `B`. With ML_K, flat universality (theory Cor. U) holds in all
  genera, modulo the same quoted theorems and Theorem U for the system.

## 1. ML_K after Stage 2, ingredient by ingredient (referees' corrections applied)

**ML_K** (MODEL_SINGULAR s.5; theory's ML_Lambda) is two statements.
* Insertions: let `h >= 0`, `n >= 1`, `2h - 2 + n > 0`, and `X_1..X_n in K` homogeneous with coefficients in `Q[v]`. Then
  `<X_1..X_n>_h(v) - F_h^(n)(v) prod_i Pi_0(X_i)` lies in `Q[v]`.
* Vertex graph: `F_h = c_h v^(2-2h)` for `h >= 2`.

Here the frozen free energies are `F_h = c_h v^(2-2h)`, `F_1 = -(1/12) log v`, `F_0 = -(5/2) v^2 log v`, and `c_h = |B_2h|/(2h(2h-2)5^(h-1))`.

| # | ingredient | status after Stage 2 | proof / evidence |
|---|---|---|---|
| M1 | **TR identification of the model CohFT**: `omega_{h,n}` of `x = W_v = p^4/4 + v/p`, `y = p` is `(-5)^(h-1) Sum <prod tau_{d_i}(H^{a_i})>^frozen_h prod d xif_{a_i,d_i}`; frozen = ancestor correlators of `M_0.omega^` at `v` | PROVED mod arXiv:1509.06954 Thm 8.1 (hypotheses checked, TR_MODEL s.1) + uniqueness of the homogeneous R-matrix (elementary recursion, TR_MODEL s.2.3); `F_h` by the dilaton equation on both sides (dilaton leg = LD at `(b,k) = (0,1)`) | VERIFIED: 72120 correlators, both primes (Stage 1); 52 frame one-point functions, `g = 2, 3`, exact over `Q(zeta5)` (`cr2_e2e`) |
| M2 | **LD, leg dictionary**: `q_{Phi_k^b} = -c_b tau_k(H^b)`, `c = (-1,1,-1,-3,-2)`; so the TR leg of `s in S_w` is `-Lambda(s)` | **PROVED here** (s.2, 7.2), given M1's functions `xif_a`. NEW: Stage 1 assumed it (critic's `cr_dictionary` docstring); the dict referee flagged it (corrections 4, 5) | VERIFIED `cr2_legdict 24`: `b <= 4`, `k <= 24`, `a <= 4`, `d <= 27`, at `v = 1` and `v = 3^5`, 35000 pairings per prime, 2^61-1 and 2^89-1; 4 controls fail |
| M3 | **D(ii)**: log coefficient of `Phi_k^b` = `c_b Sum_i Pi_i(tau_k H^b)(-x)^i/i!` | PROVED (dict referee s.3: Lemma G + residue identity), given Lemma V (gauge s.2) and `Pi_0 = PiC := (-1)^(n+1)[z^-n] Y^C_b`. The Stage-1 label "PROVED (s.6)" was sketch-level (part 2, O26) | VERIFIED: critic Test A `k <= 100` (Q); referee D1-D3 (`n <= 3000`, two primes); `cr2_legdict` L3 `k <= 24` (Q) |
| M4 | **D(iii)**: for `s in K_0 cap S_w`, `Lambda(s) in Q[p,q] = Q[v][p] + Q[v][x]` EXACTLY; `K_0 cap S_w -> (L/Q(v)[x])_w` bijective, `dim = floor(w/5)+1-[5|w]`; `N` free; Cor. 5: `K cap S_w = (K_0 cap S_w) (+) Q v^(w/5) tau_0(H)` | PROVED (DICT_PROOF s.2-s.5, s.9; the referee found no error or gap). (a), (d) need nothing; (b), (c), Cor. 4 as stated and Cor. 5 use M3. False over `Q(v)` (K7: 42 of 222) | VERIFIED exact `w <= 400`, two primes `w <= 1000` (referee); falsifier (a) |
| M5 | **Coefficient ring**: over `Q[[u]]`, `K_0` is spanned by the frames `frame(k,a) = z^k H^a + Sum_(j<k) f_j z^j H` (`a != 1`), with `f_j = const . v^((4(k-j)+c_a+1)/5)` and exponent `>= 1`. So every insertion of quintic Theorem R is a `Q[[u]]`-combination of elements of `K_0 cap S_w` and `tau_0(H)` | PROVED: theory Lemma 0 plus this exponent count. Closes dict referee correction 5(iii) | VERIFIED `cr2_legdict` L4 (`k <= 12`): legs log-free, in `Q[p,q]`, ancestor vector = frame |
| M6 | Functions of `x` pair to 0 in any slot | PROVED: linear loop equation (standard; arXiv:math-ph/0702045 s.4; statement number not re-checked) | consistency check only (logleg referee F4) |
| M7 | **L-regularity** (E9): all legs in `Q[v][p]` gives a pairing in `Q[v]` | PROVED mod arXiv:2212.00320 Thm 7.1 (statement and hypotheses checked by the Stage-1 xyswap referee, plus an off-line test) | VERIFIED in every falsifier-(b) pairing, `h <= 6` |
| M8 | **Gap Theorem**: `F_h^TR = B_2h/(2h(2h-2)) v^(2-2h)` (`h >= 2`) in the engines' convention (`F_2^TR = -1/240` on `x = W`; the kernel is MINUS the kernel printed in arXiv:math-ph/0702045 (4-3), part 2 K1); frozen `F_h = (-5)^(1-h) F_h^TR = c_h v^(2-2h)` | PROVED by two routes: arXiv:2502.18115 Cor. 1.2, formula (1.3) (its printed example s.3.2 has a sign slip, TR_MODEL s.4.1); or arXiv:2212.00320 Thm 7.1 + a Borel sum (XYSWAP s.2.3). Convention stated (logleg referee F5) | VERIFIED `h <= 6` (`h = 5, 6` by the referee's pure-log pairings) |
| M9 | **(V) = Theorem V**: (i) `a >= 1` legs in `Q[v][p]` and `k` legs `-log p` give a pairing in `Q[v]`; (ii) `a = 0` gives `F_h^(k)` (`h >= 2`; `h = 1`, `k >= 1`; `h = 0`, `k >= 3`); (iii) (V1) holds with "+" | PROVED (LOGLEG s.2-s.3; the referee re-derived every step) mod (Q1) = arXiv:math-ph/0702045 Thm 5.1, read through the referee's dictionary (`omega_{h,n} = (-1)^n W^print_{h,n}`, `F_h = -F^print_h`, `Lambda^print = -1`, "print" = the printed normalisation of arXiv:math-ph/0702045; temperature variation, its s.5.3 (5-27)-(5-30)), (Q2) = M7, (Q3) = M8. The base `h >= 2` of (ii) needs only M7 + M8 (F2). Base `h = 0` from the explicit (0,3) (LOGLEG s.2.4); base `h = 1` from homogeneity + one exact `Q(zeta5)` number, `-1/12` | VERIFIED: falsifier (b); variational formula as a tensor identity (14 cases worker, 7 referee, two primes); (V1) on 32559 (worker) and 21733 (referee) tuples; wrong-sign and no-xi controls fail |
| M10 | **Assembly** (LOGLEG_2 s.3 Corollary; F3 fixed by M4(b) and M2) | PROVED mod M1, M6, M7, M8, (Q1); steps 7.4-7.6 below | VERIFIED: genus 0-6 in the leg basis (falsifier (b) + extension); genus 0-4 in the K basis (Stage 1: 1395 singular correlators, all zero) |

**Status of ML_K: PROVED for every `h >= 0`, modulo the quoted theorems.**
* arXiv:1509.06954 Thm 8.1.
* arXiv:2212.00320 Thm 7.1.
* arXiv:math-ph/0702045 Thm 5.1, with the s.4 properties of TR.
* arXiv:2502.18115 Cor. 1.2. This one is optional: the Borel route needs only arXiv:2212.00320 Thm 7.1.

No non-quoted model-side input remains, and nothing about the quintic enters.

## 2. The assembled argument: Conjecture N(g), `g >= 2`

**Inputs that are not proved in the project.**
* **(U)** Theorem U (UR s.3). PROVED there; AI-audited only.
* **(F_g)** For the actual R-matrix of the formal quintic, the constant diagonal series `D` of UR s.4 (`Q D(z) = P_can Q^m(v(u))`) satisfies
  `D == 1 mod z^(3g-2)`. Equivalently, (H) (D is a scalar series in `z^5`) and (N) (`c^act = 1`) hold through `z^(3g-3)`.
  * At genus g only `P` through `z^(3g-3)` enters: the single-vertex leaf `z^(3g-2)` needs `P^-1_(3g-3)` (flat referee).
  * (F_g) implies (H) in the window, hence UR s.2's identification of the conifold symbol with the model's `M_0`.
* **(Q)** The quoted theorems:
  * arXiv:1509.06954 Thm 8.1;
  * arXiv:2212.00320 Thm 7.1;
  * arXiv:math-ph/0702045 Thm 5.1 and the linear loop equation of its s.4;
  * arXiv:1303.1043 Def. 2.2, Prop. 2.4, Def. 2.5, Prop. 2.9 (edge formula, composition of R-actions, translations; located in the PDF by the flat referee);
  * optionally arXiv:2502.18115 Cor. 1.2.

**Chain.** Each step is PROVED from the earlier ones and the inputs named.

1. **Vanishing solutions** (gauge s.2: Lemmas V, Q, 1).
   * `Y^C` (model) and `Y^C_Q` (quintic) lie in `Q[z^-1][[u]]^5` and equal `-e_1` at `u = 0`.
   * `[z^-1](Y^C_Q)_0 = -t_c`.
   * A regular solution is determined by its value at `u = 0`.
2. **(G)** (gauge s.3). Uses (U). Every regular gauge has `P Y^C(v(u)) = P_11(0,z) Y^C_Q`, so `P_can := P P_11(0,z)^-1` satisfies (A).
3. **Prop. S** (gauge s.4). Uses (U). `P_can(-z)^T eta P_can(z) = eta`.
4. **Transport** (theory s.3.1: Lemmas E, T, F, Cor. E' = model Lemma K). Uses 2, 3 and arXiv:1303.1043 Def. 2.2.
   * The edges of `P_can` lie in `K_0 (x) K_0`, as `Q[[u]]`-combinations of frames (M5).
   * The leaf point `d q_T` lies in `-z + K`, with `Pi_0 = d t_c = v + s` and `s in u^2 Q[[u]]`.
   * Every graph term is a finite sum (Lemma F).
5. **Quintic Theorem R** (theory s.3.2, model s.7; refereed twice).
   * If ML_K(h) holds for `h <= g-1`, every genus-g stable graph with an edge in the flat sum `P_can.Omega^mod_{v(u)}` lies in `Q[[u]]`.
   * If ML_K(g) holds as well, the vertex graph is `c_g t_c^(2-2g) + Q[[u]]`.
6. **Theorem S15** (flat s.2-s.4; referee CONFIRMED). Under (F_g), for every stable leg-free genus-g graph, the project's hatted outer term
   with `M = M_nat = M_0(tau zhat)` exactly equals the flat term of step 5. In the engine's normalisation, `G_Gamma = u^(2g-2) F_Gamma`.
   * Uses: (U) through UR s.4; Prop. S; arXiv:1303.1043 Prop. 2.4 and 2.9; the dilaton equation.
   * Also uses Lemma 3's description of the hatted variables. That description matches the `fq_sser`/`fq3_msum` docstrings, and the
     genus-2/3 matches pin it down.
7. **ML_K(h) for all h** (s.1).
   1. *M1.* `<q_1..q_n>_h = (-5)^(1-h) <Phi_1..Phi_n>^TR_h` whenever `q_i = q_{Phi_i}`, with `<Phi, w> := Sum_c Res_c Phi w`
      (`c` runs over the five branch points).
   2. *LD (M2).* The proof is in the box below.
   3. *D(ii) (M3).* The log coefficient of the leg `-Lambda(s)` is `-Sum_i Pi_i(s)(-x)^i/i!`.
   4. *M4, M5.* Take `X in K cap S_w`.
      * Write `X = Y + Pi_0(X) tau_0(H)` with `Y in K_0 cap S_w` (Cor. 5).
      * Then `leg(X) = -Lambda(Y) + Pi_0(X)(-log p)`.
      * By D(iii)(b) and Lemma 1(a), `-Lambda(Y) = Phi + f(x)` with `Phi in Q[v][p]` and `f in Q[v][x]`.
   5. *M6, M7, M9(i).* By multilinearity,
      `<X_1..X_n>_h = (-5)^(1-h) Sum_S prod_(i not in S) Pi_0(X_i) <(Phi_i)_(i in S), (-log p)^(n-|S|)>^TR_h`.
      Every term with `S` non-empty lies in `Q[v]`.
   6. *M8, M9(ii).* The term `S = {}` is `(-5)^(1-h) F_h^TR(n) prod Pi_0(X_i) = F_h^(n) prod Pi_0(X_i)`. For `n = 0` it is `F_h = c_h v^(2-2h)`.
8. **Conclusion.** Steps 5, 6 and 7 together give the following for every `g >= 2`, modulo (U), (F_g), (Q). Every stable graph with an edge
   in the hatted outer sum of the formal quintic is regular at `u = 0`, and the vertex graph is `c_g t_c^(2-2g)` + regular:
   **Conjecture N(g)**.

> **Step 7.2, LD (PROVED).**
>
> *Set-up.* Let `D = d/dx = p^2/(p^5 - v) d/dp`. Take `xif_a = eps_a v^alpha_a p^m_a/(p^5 - v)`, with `(m, eps, alpha) = (4,-1,0), (0,1,1),
> (1,1,0), (2,-1,0), (3,-1,0)`. These are TR_MODEL s.2.2 at `v = 1`; homogeneity forces the powers `alpha`. Put `xif_{a,d} = D^d xif_a`.
> The legs are `Phi_0^b = p, log p, v/p, p^3, p^2` and `D Phi_(k+1) = -Phi_k`.
>
> *Branch of `log p`.* Each leg is analytic near every branch point `c`, for any branch of `log p`. A change of branch adds
> `2 pi i f(x)`. Repeated integration by parts reduces `Res_c f(x) d xif_{a,d}` to `+-Res_c f^(d+1)(x) xif_a dx`. This vanishes because
> `xif_a dx = eps_a v^alpha_a p^(m_a - 2) dp` is holomorphic at `c`.
>
> *Four steps.*
> * (i) For `d >= 1`, integration by parts at each `c` gives `<Phi, d xif_{a,d}> = -<D Phi, d xif_{a,d-1}>`. Hence
>   `<Phi_(k+1), d xif_{a,d}> = <Phi_k, d xif_{a,d-1}>`.
> * (ii) For `k >= 1`, `<Phi_k, d xif_a> = Sum_c Res_c Phi_(k-1) xif_a dx = 0`.
> * (iii) For `d >= 1`, `<Phi_0^b, d xif_{a,d}> = -Sum_c Res_c D^d(xif_a) dPhi_0^b`. This form is meromorphic on `P^1` and holomorphic at
>   `p = 0`, because `D` maps functions holomorphic at 0 to `O(p^2)` while `dPhi_0^b = O(p^-2) dp`. At infinity it is `O(p^(1-4d)) dp`.
>   So the residues at the `c`'s add up to 0.
> * (iv) `<Phi_0^b, d xif_a> = -eps_a kappa_b v^(..) (1/5) Sum_c c^(m_a + n_b - 4)`, where `dPhi_0^b = kappa_b v^(..) p^(n_b) dp`. This is
>   non-zero only if `m_a + n_b = 4 mod 5`, i.e. only for `a = b`. There it equals `(1, -1, 1, 3, 2) = -c_b`, independent of `v` by homogeneity.
>
> Hence `<Phi_k^b, d xif_{a,d}> = -c_b delta_ab delta_dk`. □
>
> *Consequences.*
> * The TR leg of `v^m tau_k(H^b)` is `-c_b^-1 v^m Phi_k^b`.
> * `tau_0(1)` corresponds to `y = p`, the string pairing of TR_MODEL s.2.2.
> * The dilaton leg `Phi = Int y dx = -Phi_1^0` has `q_Phi = -tau_1(1)`.
> * The sign of SYNTHESIS_STAGE01 E11(ii) is confirmed.
>
> *VERIFIED.* `cr2_legdict` L1/L2: no mismatch in 2 x 35000 pairings; branch-independence `Res_c f(x) d xif_{a,d} = 0` at each `c`.
> `cr2_e2e`: `(-5)^(1-g) <leg, omega_{g,1}>` equals the frozen value from TR_MODEL's exact one-point table for all 52 frames at
> `g = 2, 3` (4 non-zero, all regular: `<frame(3,3)>_2 = 3/3200`, `<frame(4,2)>_2 = -41/9600`, `<frame(6,0)>_3 = 3/128000`,
> `<frame(7,4)>_3 = 223/2688000`).
>
> *Controls.* Step `+Int dx`; `Phi_0^2 := x`; `xif_1` without `v^1`; `log p` treated as a constant. All four fail, as they must.

**Theorem (assembled; PROVED modulo (U), (F), (Q); AI-written and AI-refereed only).** Assume three things: Theorem U; (H) and
`c^act = 1` through `z^(3g-3)`; and the quoted theorems of (Q). Then Conjecture N(g) holds for every `g >= 2`. The s.1.5 identification
is not an input: it is step 6. For `2 <= g <= 9` the input (F_g) is the finite identity E4. E4 is VERIFIED (two primes, two frames,
exact over Q on 397 Fourier components), not proved.

| input | used in step | status |
|---|---|---|
| Theorem U | 2, 3, 6 (through UR s.4) | PROVED in UR s.3; AI-audited only. A gauge-referee u-recursion confirms the solution space through `z^23` |
| (H) | 6 (in (F); the symbol `= M_0`) | CONJECTURAL; VERIFIED with (N) as E4 (`zhat <= 24`, `u < 18`), and mod `z^7`, `u^6` (`HatRowsU`) |
| (N) `c^act = 1` | 6 (in (F)) | CONJECTURAL; E4: `a_1 = a_3 = 0`. Independent coefficients are odd multiples of `z^5`; `a_5` (`z^25`) first enters at genus 10 |
| arXiv:1509.06954 Thm 8.1 | 7.1 (and M1's `F_h`) | quoted; hypotheses checked (TR_MODEL s.1) |
| arXiv:2212.00320 Thm 7.1 | 7.5 (M7); M8 route 2 | quoted; statement and hypotheses checked (Stage-1 xyswap referee s.2) |
| arXiv:math-ph/0702045 Thm 5.1 (+ s.4) | 7.5 (M9 (i)); M9 (ii) for `h = 0, 1` and the step `k -> k+1` | quoted; printed statement checked, dictionary written (logleg referee s.2) |
| arXiv:2502.18115 Cor. 1.2 | M8 route 1 (optional) | quoted; formula (1.3) re-evaluated (TR_MODEL s.4, xyswap referee) |
| arXiv:1303.1043 s.2 | 4, 5, 6 | quoted; Def. 2.2, Prop. 2.4, Def. 2.5, Prop. 2.9 located in the PDF (flat referee) |

**Corollary (universality, all genera).** Status: PROVED modulo (Q) and Theorem U for the system; flat form.
* Setting: a self-dual cyclic one-pole system with equal products, and its canonical gauge.
* Steps 1-5 and 7 give: every edge graph of the flat sum is regular, and the vertex graph is `c_g (T/lambda)^(2-2g)` + regular.
  This is theory Cor. U, with `|c0| = m^10/|f(0) g(0)|` for the bvar systems.
* B-perturbed systems are included by Lemma C.
* The hatted-engine form used in `univ/` and `bvar/` follows from step 6 with `D = 1` by construction. Matching the univ engine's
  conventions to Lemma 3 has not been written out: SKETCH.
