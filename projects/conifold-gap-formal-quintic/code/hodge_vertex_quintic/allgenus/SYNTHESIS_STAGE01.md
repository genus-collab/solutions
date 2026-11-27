# Stage 0-1 synthesis (critic): a mechanism for Conjecture N in all genera (2026-11-13)

Critic of the workflow `allgenus/`. AI-written; no human has read it. Labels: PROVED (argument written, or quoted theorem with arXiv number),
VERIFIED (exact computation; range stated), CONJECTURAL, REFUTED, SKETCH. Read: all six worker notes (`gauge/GAUGE_IDENTITY*.md`,
`univ/UNIVERSALITY*.md`, `model/MODEL_SINGULAR*.md`, `tr/TR_MODEL.md`, `xyswap/XYSWAP_MODEL.md`, `theory/REDUCTION_DRAFT*.md`) and all six
`*/referee/REFEREE.md`; no worker folder is missing. Critic code: `critic/cr_laurent.py`, `cr_conv.py`, `cr_dictionary.py` (py -I -B, stdlib,
exact over Q, under 2 s except the weight-240 run, 123 s), logs `cr_conv.log`, `cr_dictionary.log`, `cr_dictionary_60_240.log`: ALL CHECKS PASSED.

## 0. Bottom line

* The working hypothesis survives in a SHARPER form. (G) plus symplecticity puts every edge of the canonical gauge in `K_0 (x) K_0` and
  the leaf point in `-z + K` with `Pi_0 = t_c`. Here `K` is the slice where the shifted vanishing covectors `Pi_i`, `i >= 1`, all vanish,
  and `K_0 = K cap ker Pi_0`. Given the Model Lemma (ML) on `K`, Theorem R then gives Conjecture N in every genus.
  * (G) is PROVED given Theorem U.
  * The statement about edges and leaves is PROVED, in three independent write-ups.
  * Theorem R is PROVED as a conditional reduction, in two write-ups, both refereed.
  * The Model Lemma is: on `-z + K` the model potential is `F_h(Pi_0)` plus a regular part.
* The plan's non-linear `V(q)` is neither needed nor true (REFUTED for every genus-independent `V`). The LINEAR `Pi_0` is correct on `K`.
* NEW (critic): the model worker's `K_0` and the xyswap worker's TR leg class `L` are the same space.
  * The vanishing covector tower `Pi_i` consists of the Taylor coefficients of the monodromy of TR legs around the pole: PROVED, s.6.
  * `K_0` maps bijectively onto `L` modulo functions of `x`: VERIFIED exactly to weight 240.
  * Hence ML restricted to `K_0` IS the xyswap worker's L-regularity, which is PROVED modulo arXiv:2212.00320 Thm 7.1.
  * The only model-side step left is the `tau_0(H)` (log) legs. A short variational argument for it is sketched in s.6.
* The vertex constant `c_g` is PROVED for all `g >= 2`, modulo quoted theorems, by two independent routes.
* Still open: (H) and `c^act = 1` beyond `z^24` (the connection problem); the s.1.5 identification beyond genus 2; written proofs of
  D(iii) and (V) (s.6); a human check of Theorem U and Theorem R.

## 1. Established (status and scope)

| # | statement | status, scope, where |
|---|---|---|
| E1 | Theorem G: every regular gauge has `P Y^C(v(u)) = P_11(0,z) Y^C_Q`. So (G), meaning no positive `z`-power, holds iff the gauge is canonical | PROVED given Theorem U (Lemma 1 + (R1); no monodromy); VERIFIED Q, `z<=20, u<29`, 2 frames, canonical and stored (gauge + referee) |
| E2 | Prop. S: `P(-z)^T eta P(z) = c(z)c(-z) eta` | PROVED given Theorem U's uniqueness clause. The canonical gauge is symplectic; symplecticity fixes only the even part of `log c` |
| E3 | Residue relations to all orders: `(G)_n`, `(G*)_n` fix the H-column and the H^2-row of `[u^n]P`. (R1)-(R3) are orders 0-2 | PROVED given Theorem U; VERIFIED `n<=8` |
| E4 | The actual R-matrix's gauge equals the canonical one in the window `zhat<=24, u<18`: `c^act = exp(a_1 z^5 + a_3 z^15 + ..)` with `a_1 = a_3 = 0` | VERIFIED: 2 primes, 2 frames, and exact Q on 397 Fourier components. All orders: CONJECTURAL. The window covers genus `<= 9` |
| E5 | Edge/leaf transport (Lemma L = Lemma E; (E), (T); Lemma K = Cor. E' + Lemma T). Edges annihilate every `Pi_m`, `m >= 0`, so they lie in `K_0 (x) K_0`. The leaf point lies in `-z+K`, with `Pi_0 = t_c` and `s = dT_P - v = u^2/30 + (127/3150)u^3 + ..` (frame independent) | PROVED (symplecticity + (G) only); VERIFIED hatted mod p (`zhat<=24`) and flat over Q, with failing controls |
| E6 | Theorem R: flat-regular + symplectic + annihilating `P`, plus ML(`h <= g-1`) => every genus-g edge graph is regular. Adding ML(g) => the vertex graph is `c_g T_P^(2-2g)` + regular. The u-ODE is never used | PROVED, conditional (model s.7, theory s.3.2; both referees: no gap). VERIFIED in the flat frame at genus 2, per graph = hatted record, incl. a non-gauge `P` (theory referee) |
| E7 | Leaf bookkeeping `-(g-1)c_g/15` is derived, and is frame independent | PROVED given ML(g) (theory Cor. L); `cr_conv` (4) |
| E8 | Model gap: for `x = P(z) + v/(z-a)`, `y = z`, `F_g = B_2g/(2g(2g-2)) v^(2-2g)` for all `g >= 2` and any polynomial part. So the frozen `F_g = c_g = abs(B_2g)/(2g(2g-2)5^(g-1))`, and the primary-sector lemma holds for `g >= 2` | PROVED by two routes, each modulo a quoted theorem: arXiv:2502.18115 Cor. 1.2 (tr), or arXiv:2212.00320 Thm 7.1 plus a Borel sum (xyswap). Plus arXiv:1509.06954 Thm 8.1 with `Lambda = -5 kappa`; the Phi-pairing for all `d` was closed by the tr referee. VERIFIED `g <= 5` by TR |
| E9 | L-regularity: correlators whose legs all lie in `L = Q[params][z] + Q[params][x]` are polynomials on the Hurwitz space, so every `omega_(g,n)`, `n >= 1`, is regular at fixed `p` | PROVED mod arXiv:2212.00320 Thm 7.1 (xyswap; referee checked off the line). This upgrades tr item 4 from CONJECTURAL |
| E10 | ML on `K` (= theory's ML_Lambda) | VERIFIED genus 0-4: 1395 singular K-correlators, all zero; referee adds `(0,7)`, `(1,6)` and controls. PROVED on the `H^2` tower (arXiv:2211.12259 Thm 2). In general CONJECTURAL; reduced in s.4 |
| E11 | Dictionary D (critic). (ii) The TR leg of `tau_k(H^b)` is `phi(p) + log p . f(x)`, with `f(x) = -Sum_i Pi_i(tau_k H^b) (-x)^i/i!`. (iii) `K_0 -> L/(functions of x)` is bijective, weight by weight | (ii) PROVED (s.6); VERIFIED `k<=60`, all `b`, 61 non-zero + 244 zero slots. (iii) VERIFIED weights `<= 240`: 5880 basis vectors, dimensions equal (`cr_dictionary`) |
| E12 | Universality, genus 2-4 | VERIFIED (2^89-1; referee two primes at genus 4). Canonical gauges of perturbed systems: all 378 regular, vertex in the system's own `T`. Annihilating NON-gauge dressings are regular. Stored normalisation: exactly graphs 1, 2 (`[u^4]`) and 3, 4, 7, 8 (`[u^5]`) are singular, and the total is unchanged |
| E13 | The first attempt's (G*) failure was a float bug. (G*) is a theorem | VERIFIED reconstruction (gauge + referee) |

## 2. Refuted or weakened

* **Working hypothesis, the parts that fail.**
  * The plan's linear one-point law for every insertion is REFUTED off `K`: `<tau_5(H)>_3 = 53/1209600` against the law's `-1/604800`; the genus-4 remainders are `(19/100800000) Pi_5`. It holds exactly on `K`.
  * The plan's non-linear genus-independent `V(q)` is REFUTED: three pairs have genus-dependent Hessians. Genus-mixing forms (hbar-dependent, Toda) are not excluded, but they are not needed.
  * "The edge condition is (G) plus compatibility with the non-linear structure" is REFUTED as a requirement. (G) already annihilates the whole flag `Pi_m`. The string counterexample `<tau_1(H)tau_0(1)> = <H>` is off `K` and never occurs in a graph.
  * "The edges are built from symmetry directions" is REFUTED: random annihilating dressings are regular (univ E3).
  * The worry "(G) alone is not enough" survives only in this form: the model-side ML is needed. No further gauge-side condition is needed.
* **Gauge.**
  * The "frame dependence" of the leaf bookkeeping is REFUTED. It is `1/15` in every frame, since `gamma_t + w_1 = 7/10`.
  * "The mechanism needs the non-linear genus-0 descendant extension" (HAT s.5) is superseded by Lemma K plus ML on `K`, both linear.
* **Model.**
  * "Genus 4 is NOT a sharp test of the `z^5` normalisation", and its referee's "cannot see it at `u^-2`", are REFUTED for Conjecture N. `Delta_4(5) = 0` is the TOTAL; per graph the response is non-zero (s.3 C1).
  * "Rules out any jet expression" is WEAKENED to genus-independent `V` only.
  * Lemma K's hypothesis `phi_j -> 0` is false in the hatted frame. The conclusion holds, since edges are polynomial at each u-order.
* **Theory.**
  * "1013 singular lists" is WEAKENED: only 326 are genuine cancellations, and the deep frames are untested by th_mla.
  * Cor. N / N_9 omit the s.1.5 identification from their hypotheses.
  * "Each hypothesis is needed in examples" is WEAKENED: there is no example for (Sp).
* **Univ.**
  * "Conjecture A + (G) => Conjecture N" also needs (H) and `c^act = 1`.
  * "Flatness is not sufficient" holds by construction (Theorem U's freedom `c(z^5)`).
* **xyswap.**
  * "`Y^C` and the TR monodromy are different functionals; `Y^C` pairs non-trivially with L legs such as the unit" is REFUTED (critic, D(ii)). `PiC` equals minus the monodromy constant on every insertion and vanishes on `K_0 = L`. The non-zero pairing of `Y^C_0 = I^C_0` is with `tau_1(1), tau_6(1), ..`, which are log legs, not with `tau_0(1)`.
  * "(M1)/(M2) is the remaining untested step" is SETTLED conditionally: it is Lemma K plus D.
  * The referee's obstruction "every descendant acquires `log z`" dissolves: the logs are cancelled exactly by the calibration terms of `kappa(X,k)`.
  * "Mixed correlators are singular in general" is WEAKENED.
  * `F_1` is not computed from its definition.
  * Toda holds only up to a constant.
* **tr.**
  * "`x = -W` flips `psi -> -psi`" holds only up to `(-1)^(3g-3+n)`.
  * "any point" should read "any simply ramified point".
  * The metric-sign remark is incomplete. It is `(-1)^(g-1+D)`; harmless.
* **Record notes now outdated.** These are upgrades, not errors:
  * NATURAL_FORM s.2: the bookkeeping is "not derived". It is now derived.
  * NATURAL_FORM s.3: the "model gap known only to genus 3" is now known for all `g`, modulo quoted theorems.
  * UNIFORM_REDUCTION s.4: the canonical normalisation is "one number". It is now verified through `z^24`.

## 3. Contradictions between notes, and which side is right

* **C1. Is the `z^5` normalisation visible at genus 4?**
  * model, and its referee: no, `Delta_4(5) = 0`.
  * univ: yes, graphs 1 and 2 are singular at `[u^4]`.
  * The numbers agree exactly. Per unit `c_5` they are:
    * vertex `-1/8400000`;
    * bridge (1,3) `-1/1008000`;
    * bridge (2,2) `1/900000`;
    * loop 0.
  * The two bridges add to `1/8400000`, the model's "bridges" value, and the sum over graphs is 0 (`cr_conv` (3)).
  * The model's statement concerns the total. Conjecture N is per graph, so **univ is right**: genus 4 IS a sharp per-graph test.
* **C2. Vertex constant.** These four forms are equal for `g = 2..12` (`cr_conv` (1); `c_6 = 691/1023750000` is a prediction):
  * tr: `(-5)^(1-g) F_g(x=W)`;
  * brief, NATURAL_FORM: `abs(B_2g)/(2g(2g-2)5^(g-1))`;
  * UNIFORM_REDUCTION, model: `(-1)^(g-1) B_2g/(..)5^(1-g)`;
  * LG_MODEL s.5.4: `20^(1-g)(-4)^(3g-3) F^Hurw`.

  The tr worker's "`5^(1-g)` is wrong" concerns the TR-to-frozen factor for `x = W`, not `c_g`. On `x = p^4 + 4/p`, the tr and xyswap free energies agree: `-1/3840`, `1/258048`.
* **C3. `<H>_1 = -1/12` and `<e>_1 = +1/8`, against `+1/12` and `-1/8`.**
  * The first pair is from tr, model, xyswap and theory; the second from UNIFORM_REDUCTION s.6.3 and the plan.
  * The difference is convention, not substance: the project's `H` is `-d/dv` at its frozen point `t = -1 <-> v = -1`, and the sign of `e` is free.
  * Advice: fix `H = d/dv` and `e = tau_1(1) + v tau_0(H)` from now on.
* **C4. Signs of the vanishing covector.** All consistent (`cr_conv` (2), `m <= 60`; the tr worker's "signs not matched" is now matched):
  * model: `PiC(tau_k H^b) = (-1)^(k+1)[z^-k] Y^C_b`;
  * gauge: `phi = -PiC` times a positive factor;
  * theory: `V^(1) = PiC`, and its frames lie in `ker PiC`;
  * closed form: `1/(4^q q!(m-q)!)`, with support rule `b = 1-m`.
* **C5. Are `Y^C` and the TR monodromy the same functional?**
  * xyswap: different functionals.
  * gauge and model, implicitly: the same.
  * **Critic: the same.** Proof in s.6. VERIFIED `k <= 60`; the uniform sign `mu_0 = -PiC` holds for the natural legs `Int d_b W dp`.
* **C6. Leaf bookkeeping.**
  * gauge: frame dependent.
  * UNIFORM_REDUCTION s.6.2, gauge referee, theory: frame independent.
  * **Frame independent is right**: `2 gamma_t + 2 w_1 - 4/3 = 1/15` and `p_14(0)/2 = -1/30` in both frames (`cr_conv` (4)).
* **C7. Euler-dilaton edge leg.**
  * xyswap: `-(2/5) z^9 + (v/5) z^4`.
  * xyswap referee: `-(16/45) z^9 + (v/5) z^4`.
  * **The referee is right** (`cr_conv` (5)). The leg is still in `L`. The chain `psi^k e` is log-free for `k <= 4` and has a log at `k = 5`, consistent with `psi^k e in K_0` iff `k <= 4`.
* **C8. String-equation sign.** xyswap's prose has `-d/dx`; the referee has `+d/dx`. **The referee is right** (TR-verified there; not re-checked here).
* **C9. Name clashes, to rename.** "Theorem G" = gauge identity (gauge) and gap theorem (tr): use "(G)" and "Gap Theorem". "K", "Lemma K"
  (model) clash with UNIFORM_REDUCTION's frame `K = (1,uH,H^2,H^3,H^4)` and the plan's failed "Lemma K": use `K_vc`. "ML1/ML2" have swapped
  roles (model: pure-H sector / kappa-regularity; xyswap: L-regularity / gap). "Theorem R" also names the conifold theorem.
* **C10. Theory Cor. U's `lambda`.**
  * theory: `1/lambda = P_11(0,z)` with `lambda = f(0)`.
  * univ referee: `c0 = 1/(f(0)|g(0)|)`.
  * Define `lambda` by `r = -1/lambda^2`, i.e. `lambda = f(0)|g(0)|`. The two agree in every run made, since `g(0) = -1`.
* **C11. Sign of `c_5`.** GAUGE s.4 uses one sign, s.5-s.6 the other; the value is right. Unify the convention: `P_stored = P_can (1 + m5 z^5 + ..)`, `m5 = -43397/807975000`.

## 4. The sharpest surviving mechanism statement

**Statement.** For every `g >= 2`, the following implies Conjecture N(g):
* [Theorem U];
* [(H) and (N): `c^act = 1`, in the window of genus g];
* [s.1.5: the flat `P`-action graph sum equals the project's hatted outer sum, graph by graph];
* [ML(`h <= g`) on `K`].

The proof runs through (G) (E1), Prop. S (E2), Lemma K (E5) and Theorem R (E6).

Moreover, ML on `K` follows from:
* [gap, E8];
* [D(ii) PROVED];
* [D(iii)];
* [L-regularity, E9];
* [(V): adding `-log p` legs keeps polynomiality];
* [arXiv:1509.06954 Thm 8.1 with the tr normalisation].

| ingredient | status |
|---|---|
| Theorem U (existence; uniqueness up to `c(z^5)`) | PROVED, AI-audited only; the gauge referee's independent u-recursion confirms the solution space through `z^23` |
| (H) purity; (N) `c^act = 1` | CONJECTURAL; VERIFIED `zhat<=24, u<=17` (independent content `a_1 = a_3 = 0`). Enough for genus `<= 9` |
| s.1.5 identification | argued (frame changes, rescaling, dilaton factor); VERIFIED at genus 2 only (theory referee, flat sum) |
| (G), Prop. S, Lemma K / E', Lemma T | PROVED given Theorem U |
| Theorem R | PROVED, conditional; refereed twice |
| ML: gap and pure-H sector (`h >= 2`); `h = 0, 1` from `c_555 = 1/a_5` and `<H>_1 = -1/(12v)` | PROVED mod quoted theorems (`h >= 2`); `h = 0, 1` VERIFIED |
| ML: correlators with only `K_0` insertions | follows from D(iii) and L-regularity (PROVED mod arXiv:2212.00320 Thm 7.1). D(iii) VERIFIED weights `<= 240` |
| ML: mixed with `tau_0(H)` | (V) SKETCH (critic, s.6). VERIFIED genus 0-4 through the model data |

## 5. What universality implies for the shape of the proof

* **Only three properties of the left factor matter**: flat regularity, symplecticity, and annihilation (A).
  * Non-gauge dressings are regular at genus 2-4 (two constructions, two primes).
  * At first order, at genus 2, regularity is EQUIVALENT to annihilation.
  * So the proof must be edge-local, as in Theorem R. It must not use the u-ODE or any quintic-specific datum except `B = prod b_a`, through `v(u)` and `tau`.
* **The normalisation is invisible in the total.**
  * Genus 4, stored gauge: the total's polar part is unchanged through `u^5`, while six edge graphs and the vertex graph move.
  * So no argument about total potentials (KP/Toda for the full theory, total free energies, the gap of the total) can prove Conjecture N.
  * The proof must be graph-by-graph, and must use the ODD part of `log c` (fixed by (G)), which symplecticity does not see.
* **The vertex law is in the system's own coordinate.** It reads `c_g (T/lambda)^(2-2g)` with `T = Int b_0`, not in `t_c`.
  * `t_c` enters only as `Pi_0` of the leaf point (Lemma T).
  * The vertex constant is a pure model statement: the gap, now proved.
* **Flat-irregular dressings fail.**
  * First-order invisibility is not stable under products (C1c).
  * The genus-3 surplus directions die at genus 4.
  * So the proof must work in the flat frame with exact u-series (Lemma F finiteness), not in a truncated hatted window.
* **Not yet tested.** `B` was never varied. Theorem R predicts B-universality, which is a cheap extra check.

## 6. Stage 2: route, first lemma, falsifying computations

**Route (recommended): TR on the model curve `x = W = p^4/4 + v/p`, `y = p`, with legs.** This is not hypermap combinatorics, KP/Toda, or a codimension-2 extension: TR is allowed for the model, and both quoted inputs are in hand.

1. **D(ii) (PROVED here).**
   * Set-up:
     * Natural legs: `Phi_0 = Int d_b W dp`, with `d_b = (d_a4, -d_a5, -d_a1, d_a2, d_a3)` for `b = 0..4`.
     * Descendants: `Phi_(k+1) = -Int Phi_k dW`. This is the dual of the descendant functions of arXiv:1509.06954 (2-35); the tr worker's `y` and `-Int y dx` are the cases `b = 0`, `k = 0, 1`.
   * Argument:
     * The generating function `G = Sum_k z^k Phi_k` solves `dG/dW = dPhi_0/dW - zG`.
     * Hence `G = e^(-zW)(Int^p e^(zW) dPhi_0 + function of W)`.
     * Its monodromy around `p = 0` is `2 pi i e^(-zW) Res_(p=0) e^(zW) d_b W dp`.
     * By the LG dictionary (class of `H^b` = `d_b W dp`) and `Y^C` = residue period, this equals `2 pi i e^(-zW) Y^C_b(zeta = -1/z)`.
     * Taking `[z^k]` gives `f_k(x) = -Sum_i PiC(tau_(k-i) H^b)(-x)^i/i!`.
   * Consequences:
     * `K` corresponds to `L + C log p`, and `K_0` to the log-free legs.
     * The D-tower `tau_k(H)` corresponds to legs whose monodromy is not constant.
2. **First lemma to prove: D(iii) for all weights.** Statement: a log-free leg built from `v^m tau_k(H^b)` lies in `Q[v][p]` modulo functions of `x`.
   * Equivalently: `K_0 = L/Q(v)[x]`. Both sides have dimension `floor(w/5) + 1 - [5|w]`.
   * Route: injectivity of the leg map modulo `Q(v)[x]` (non-degeneracy of the residue pairing), plus the principal part at `p = 0` of the monodromy-free part of `G`.
3. **(V), log legs (SKETCH).**
   * Let `Phi_j` be in `Q[v][p]`, and `C(v) = (-1)^n Res_inf..Res_inf prod Phi_j omega_(h,n)`; it is a polynomial by L-regularity.
   * Use `d_v|_p = d_v|_x + Lie_xi`, with `xi = p/(p^5-v) d_p`.
   * The TR variational formula (arXiv:math-ph/0702045 s.5) applies with `d_v(y dx)|_x = -dp/p = Int_0^inf B(., p')`. So `d_v|_x omega_(h,n)` is the pairing with the extra leg `-log p'`.
   * The `Lie_xi` terms pair `m p^m/(p^5 - v)`, expanded at infinity. That gives finitely many polynomial legs; negative powers pair to zero at infinity.
   * So `<Phi_1..Phi_n, -log p>_h = C'(v) - (polynomial)`. Induct on the number of log legs.
     * Boundary terms of `Int_0^inf Lie_xi` vanish, since `xi(0) = 0` and `omega = O(dp/p^2)` at infinity.
     * Genus 0 starts from `omega_(0,2) = B`.
     * With no L-leg one gets `F_h^(m)`, the gap.
4. **Assemble.** ML on `K`, then Theorem R, then Conjecture N, modulo Theorem U, (H), (N) and s.1.5.
5. **In parallel.**
   * s.1.5 at genus 3 in the flat frame (extend `theory/referee/rf_g2flat.py`).
   * Stage 4: (H) and `c^act = 1`. The next independent coefficient is `a_5` at `z^25`, first used at genus 10.
   * A human check of Theorem U, Lemma E and Theorem R.

**Cheapest falsifying computations.**
* (a) D(iii) beyond weight 240: `py -I -B critic/cr_dictionary.py 80 320`, a few minutes. A single log-free leg outside `L` would break the route.
* (b) (V) at the first range the model data do not reach:
  * Compute TR `omega_(4,2)`, i.e. the model worker's missing `amp(4,2)`, and `omega_(3,3)`, mod two primes with the engines in `tr/` or `xyswap/`.
  * Pair them with every L-leg and with `-log p`.
  * ML predicts that every pairing with negative v-exponent vanishes. Expected cost: minutes per prime (the tr engine did `(4,1)`, `(3,2)` in
    under a minute; not timed for `(4,2)`), against about 4 CPU hours for the model engine.
* (c) B-universality: one univ run with `B` perturbed at genus 3.

## 7. Overclaiming sentences (quote, then correction)

* GAUGE s.0 item 5: "the production R' of genus 2-5 IS the canonical gauge at the orders used". This was tested on the referee rows. On the production rows it was tested only to `zhat <= 6`; beyond that it holds by transitivity, through NATURAL_FORM s.1 item 1.
* GAUGE_HAT s.5: "including its frame dependence -(g-1)(2 gamma_t + 2 w_1 - 4/3) c_g". The value is frame independent.
* GAUGE summary and HAT s.5: "needs the non-linear genus-0 descendant extension". Superseded: Lemma K plus ML on `K`.
* MODEL s.0 item 8 and s.9: "genus 4 is NOT a sharp test of it". Model referee item 10: "genus 4 cannot see the z^5 normalisation at u^-2". Both are true only for the total (C1).
* MODEL s.4: "any V, linear or not, jets of the genus-0 solution included". Read: any genus-independent `V`.
* MODEL_2 s.6, Lemma K: "phi_j -> 0". This fails in the hatted frame.
* THEORY s.0 item 3 and s.2.4: "VERIFIED on 1013 singular correlator lists". Of these, 326 test the frame corrections, and none tests the deep frames.
* THEORY Cor. N_9: "follows from ML_Lambda(h), h <= g". It also needs s.1.5 (genus 2 only) and Theorem U. Its "finite list" was never written down.
* UNIV s.5, Conjecture A: "With (G) ... this implies Conjecture N". It also needs (H) and `c^act = 1`.
* UNIV s.0 item 4: "flatness is neither sufficient nor necessary". "Not sufficient" holds by construction.
* XYSWAP s.3: "Y^C pairs non-trivially with L legs such as the unit, the monodromy does not". REFUTED (C5).
* XYSWAP s.0 item 4: "11 L-correlators that homogeneity would allow to be singular vanish". These zeros are forced by decay at infinity and are not evidence.
* XYSWAP s.0 item 5: "Mixed correlators ... are singular in general". Too strong: `(z^5, -log z)` gives 0.
* XYSWAP s.5: "Toda (PROVED given s.2.3)". It holds only up to a constant fixed by the genus-0 normalisation.
* TR s.0 item 4 understates: regularity at fixed `p` is PROVED mod arXiv:2212.00320 Thm 7.1 (E9).
* TR s.0 item 2: "x = -W, which flips psi -> -psi". This holds only up to `(-1)^(3g-3+n)`.

**Hard-rule deviations found by the referees, to be fixed by the workers.** Author names or initials in GAUGE s.6 (heading tag),
MODEL_2 s.8, TR_MODEL and one tr script name, XYSWAP (names, an acronym). Logs over 250 lines: `univ/e4.log` (444), `univ/un_e2.log` (347),
`model/ms_singlist_full.log` (987). Budgets: univ made 14 genus-4 runs (6 allowed); the theory referee briefly ran more than 2 processes.
