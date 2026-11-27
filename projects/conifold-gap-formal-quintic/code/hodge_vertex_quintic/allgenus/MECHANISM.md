# All-genus mechanism for Conjecture N: Stage 0-1 (2026-11-13)

Entry point of `allgenus/`. Written from the results of the 13-agent workflow (6 workers, 6 referees, 1 critic: all referees "sound
with corrections") and from my own re-runs (last section). AI-written; no human has checked the proofs. Labels as in `../NATURAL_FORM.md`:
PROVED (argument written, or a quoted theorem with arXiv number), VERIFIED (exact computation, range stated), CONJECTURAL, SKETCH, REFUTED.
Full synthesis: `SYNTHESIS_STAGE01.md` (critic). Worker notes: `gauge/GAUGE_IDENTITY*.md`, `univ/UNIVERSALITY*.md`, `model/MODEL_SINGULAR*.md`,
`tr/TR_MODEL.md`, `xyswap/XYSWAP_MODEL.md`, `theory/REDUCTION_DRAFT*.md`; referee reports in `*/referee/REFEREE.md`.

Notation. `Y^C` = the pure-`1/z` (residue, vanishing-cycle) solution of the model system, `Y^C_Q` the quintic's; `Pi_i` = the shifted
vanishing covectors; `K = {Pi_i = 0 for i >= 1}`, `K_0 = K cap ker Pi_0`; ML (Model Lemma) = on `-z + K` the model ancestor potential is
`F_h(Pi_0)` + regular; canonical gauge = the Theorem-U gauge with `P_11(0,z) = 1`.

## 1. The mechanism (sharpest surviving form)

For every `g >= 2`, Conjecture N(g) follows from four inputs: Theorem U; (H) together with `c^act = 1` in the window of genus g; the
identification of the flat `P`-action graph sum with the project's hatted outer sum (s.1.5 of the theory note); ML(`h <= g`) on `K`.
The chain:
1. **(G)**: every regular gauge has `P Y^C(v(u)) = P_11(0,z) Y^C_Q`, so "no positive powers of z" holds iff the gauge is canonical.
   PROVED given Theorem U (no monodromy needed). VERIFIED exactly over Q (`z <= 20`, `u < 29`, two frames).
2. **Symplecticity**: `P(-z)^T eta P(z) = c(z)c(-z) eta`; the canonical gauge is symplectic. It fixes only the even part of `log c`, so
   only (G) fixes the odd multiples of `z^5`, e.g. `c_5`. PROVED given Theorem U.
3. **Edges and leaves**: every edge annihilates every `Pi_m` (edges lie in `K_0 (x) K_0`); the leaf point lies in `-z + K` with
   `Pi_0 = t_c`. PROVED given Theorem U, three independent write-ups; VERIFIED on the actual rows.
4. **Theorem R** (conditional reduction): with ML(`h <= g-1`) every genus-g graph with an edge is regular; with ML(g) the vertex graph is
   `c_g T^(2-2g)` + regular, and the single-leaf bookkeeping `-(g-1)c_g/15` is derived (the same in every frame). PROVED, refereed twice.

ML on `K` itself reduces to: the gap `F_h = c_h v^(2-2h)` (PROVED for all `h >= 2` modulo arXiv:2502.18115 Cor. 1.2, or arXiv:2212.00320
Thm 7.1, together with arXiv:1509.06954 Thm 8.1); the dictionary D(ii) "the `Pi_i` tower = Taylor coefficients of the monodromy of the TR
legs at the pole" (PROVED, critic); D(iii) "`K_0` = the TR leg class `L` modulo functions of x" (VERIFIED for weights <= 240); L-regularity
(PROVED modulo arXiv:2212.00320 Thm 7.1); and (V), legs `tau_0(H)` = `-log p` (SKETCH). Independently, ML on `K` is VERIFIED at genus 0-4
(1395 singular `K`-correlators, all zero) and PROVED on the `H^2` tower (arXiv:2211.12259 Thm 2).

## 2. What the evidence now says

* **Normalisation matters graph by graph.** The actual R-matrix's gauge equals the canonical one through `zhat^24`, `u^17` (two primes,
  two frames; exact on 397 Fourier components), enough for genus <= 9; to all orders this is CONJECTURAL (a connection problem).
  In the stored normalisation (`c_5 = -43397/807975000` away from canonical) six genus-4 edge graphs become singular (graphs 1, 2 at
  `[u^4]`, 3, 4, 7, 8 at `[u^5]`) while the total polar part does not move. So the proof must be graph-local, and no argument about total
  potentials (KP/Toda of the full theory, the gap of the total) can prove Conjecture N.
* **Universality** (VERIFIED genus 2-4): canonical gauges of perturbed cyclic systems (not from any Picard-Fuchs equation) make every
  edge graph regular, with the vertex graph `c_g (T/lambda)^(2-2g)` in the system's own coordinate `T`; annihilating NON-gauge dressings
  are regular too (genus 4, two primes); violating annihilation or flat regularity breaks regularity. What matters is flat regularity,
  symplecticity and annihilation of the vanishing covector, not the u-ODE.
* **Refuted**: the plan's non-linear `V(q)` (no genus-independent `V` works, and none is needed: the linear `Pi_0` is right on `K`); the
  linear one-point law off `K` (`<tau_5(H)>_3 = 53/1209600` against `-1/604800`); "edges are built from symmetry directions".

## 3. Open inputs

Theorem U (AI-audited only); (H) and `c^act = 1` to all orders; s.1.5 beyond genus 2; written proofs of D(iii) and (V); a human check of
Theorem U, Lemma E and Theorem R.

## 4. Stage 2 (not started)

Recommended route (critic): TR on the model curve `x = p^4/4 + v/p`, `y = p`, with legs. First lemma: D(iii) for all weights; then (V).
Cheapest falsifiers first: `critic/cr_dictionary.py 80 320` (minutes); TR `omega_(4,2)` and `omega_(3,3)` paired with all `L`- and
`-log p` legs mod two primes; one genus-3 run with `B = prod b_a` perturbed. In parallel: s.1.5 at genus 3 in the flat frame; Stage 4,
the connection problem for (H) and `c^act = 1` (next independent coefficient `a_5` at `z^25`, first used at genus 10).

## 5. Verification done outside the workflow (2026-11-13)

* `critic/cr_conv.py` re-run: 29 PASS, 0 FAIL, check lines identical to `critic/cr_conv.log`.
* `critic/cr_dictionary.py 60 240` re-run: 5880 log-free basis vectors all in `L`, output identical to `critic/cr_dictionary_60_240.log` (146 s).
* `gauge/gi_relations.py 20 66 0 0` re-run from a scratch copy (so that the project log was not overwritten): PASS/FAIL lines identical to
  `gauge/gi_relations_ab_0_0.log` (35 s).
* Read `univ/e4.log` and `univ/referee/rf_e4.log`: the referee's module-disjoint engine reproduces the worker's genus-4 tables (0 of 379
  outer rows differ at 2^89-1; graphs with at most two edges also at 2^61-1), including the six singular graphs of the stored normalisation.
* Not re-run: the genus-4 graph sums (hours each), the TR and x-y swap computations, the model correlator tables.

## 6. Housekeeping flagged by the critic (not yet done)

Conventions: use `H = d/dv` and `e = tau_1(1) + v tau_0(H)`; unify the sign of `c_5` (`P_stored = P_can (1 + m5 z^5 + ...)`,
`m5 = -43397/807975000`); rename the clashing labels ("Theorem G", "Lemma K"/`K`, "ML1/ML2", "Theorem R"). Rule deviations in worker
output: three logs over 250 lines (`univ/e4.log`, `univ/un_e2.log`, `model/ms_singlist_full.log`); author names or initials in some notes.
