# Referee report: Stage 2, "s.1.5 identification at genus 3" (flat worker), 2026-11-14

Adversarial referee of `allgenus/stage2/flat/` (S15_PROOF.md, s15_*.py, logs/). AI-written; no human has read it. Labels: PROVED
(argument written, or quoted theorem with arXiv number and statement number), VERIFIED (exact computation; range and primes stated),
CONJECTURAL, SKETCH, REFUTED; for the worker's claims also CONFIRMED / CORRECTED. Conventions as in `../../../SYNTHESIS_STAGE01.md` s.3.
Read: S15_PROOF.md, all worker scripts and logs, MECHANISM.md, SYNTHESIS_STAGE01.md (s.1, s.3, s.4, s.6), REDUCTION_DRAFT s.1.3-s.1.5,
theory/referee/REFEREE.md and rf_g2flat.log, rf_flat_core.py, UNIFORM_REDUCTION.md (UR) s.1-s.4, ur_gauge.py, the docstrings of
fq_sser.py and fq3_msum.py, NATURAL_FORM.md s.5 (genus-3 line). Own code: `referee/rf_s15_*.py` (python -I -B, stdlib, mod 2^89-1 and
2^61-1, exact over Q where inputs are rational); logs in `referee/logs/`. No worker module is imported.

## 0. Verdict: sound with corrections

* **Theorem S15 (PROVED given (F)): CONFIRMED.** I checked every line (s.1 below). The argument is correct. In substance it is the
  composition of R-actions (arXiv:1303.1043 Prop. 2.4, Prop. 2.9), the dilaton equation and a degree count. Its only non-mathematical
  input is the description of the project's hatted code that Lemma 3 uses. I checked that description against the docstrings of
  fq_sser.py and fq3_msum.py (it matches), and the exact matches below pin it down. The consequence "s.1.5 is implied by Theorem U,
  (H), (N) and can be deleted from SYNTHESIS s.4 as a separate input" is CONFIRMED.
* **Every VERIFIED claim was re-derived with an independent engine** (s.2): own Witten-Kontsevich numbers, own graph enumeration, and
  an own frozen R-matrix derived from the Frobenius algebra along the v-line. Its unit row equals UR's F_0 through t^7, and it gives
  F_2 = 1/1200, F_3 = 1/25200, F_4 = 1/180000. The flat graph sum is computed as one red/blue Givental sum.
  * Results, frame (0,0), both primes: flat = production record A3 (k <= 3); flat = worker's hatted M_nat table and = worker's flat
    table (k <= 7, 336 coefficients); sum over graphs = total_id; all 41 edge graphs regular; vertex polar part = (1/25200) t_c^-4.
  * Frame (1/3,-2/7) at 2^89-1: the same.
  * Genus 2: = rf_g2flat.log in both frames.
  * All controls fail as they must.
* **New check R7 (stronger than F1).** I regrouped the same composite terms by the uncontracted graph, which gives the full R-matrix
  action graph by graph. This equals the stored point-frame production table A_id on all 252 coefficients (42 graphs, k <= 5), all of
  them non-zero, at both primes and in the second frame. The stored normalisation and a symplectic non-gauge dressing both break it.
* **Corrections** (s.3), none affecting the conclusion:
  * the count "112 non-zero regular coefficients of edge graphs" should be 108 (112 includes the vertex graph);
  * "independent confirmation of (N) at z^5" overstates what E4 does not already contain;
  * X4 is orientation-dependent;
  * F1's 168 coefficients are 164 zeros, so the non-zero per-graph comparison rested on the worker's own hatted driver (now
    supplemented by R7);
  * four cosmetic slips in s.1.3 and Lemma 3.

## 1. Line-by-line check of S15_PROOF.md

| item | verdict | remark |
|---|---|---|
| 1.1 definitions (E_O, q_O, Cont) | CONFIRMED | `E_O` is a power series by symplecticity (arXiv:1303.1043, Remark 2.1 and eq. (6)). For `O_0 = 1` this is Def. 2.2 (edges, legs) plus the translation of s.2.2-2.3 there. I located Def. 2.2, Prop. 2.4, Def. 2.5 and Prop. 2.9 in the PDF. The dilaton identity `d^n F_h(q/k) = k^(2h-2+n) d^n F_h(q)` holds for `2h-2+n > 0`, `n >= 1` at `h = 1`. Unstated: here `k = d(u)` has `d(0) = 1`, so all dilaton re-expansions converge u-adically |
| 1.2 flat sum | CONFIRMED | Observation: `prod_v d^(2h_v-2+n_v) = d^(2g-2)` for every graph (`Sum_v h_v = g-1-E+V`). The dilaton factor is therefore one global unit per graph, and the only per-vertex effect of `d` is the `d` in the leaves |
| 1.3 hatted sum | CONFIRMED, one slip | The leaf formula "`T'^(i)(z) = z - z Q'^(i)_0(-z)`" is wrong at `z^1`. The correct form is `T' = z y + q_{R'}`, where `y_i = Q'^(i)_0(0) = Delta_i^(-1/2)` is the unit of `omega^Q` (`Delta_i^(1/2) = zeta^(4i) x`, fq_sser). The code uses only `d >= 2` (`T'_d = (-1)^d Q'_{0,d-1}`, fq_sser docstring), so nothing downstream is affected. The edge formula equals `E_{R'}` because `R'` is symplectic. That holds for the ACTUAL `R' = R~ M_nat^-1` with or without (F), since `R~` and `M_nat` are symplectic |
| Lemma 1 | CONFIRMED | insertion of `A A^-1` |
| Lemma 2 | CONFIRMED | Correct when "q absolute" means: each CohFT's potential is evaluated at `t = q + z (its own unit)`. Then `F^(k.W)_h(q) = k^(2h-2) F^W_h(kq) = F^W_h(q)` by dilaton homogeneity of degree `2-2h` |
| Lemma 3 | CONFIRMED given the code description | Description: rows `Qhat = Q s^(4n+1-c_a)`, `xhat = x/s`, leaves from row 0. This matches fq3_msum.py (docstring line 5) and fq_sser.py. Degree count: `Sum_v (3h_v-3+n_v) + E = 3g-3`, so `s^(12(g-1))`; `c_a + c_b = 2` for `a+b = 3 mod 5`; TFT and leaves give `s^-(2g-2)`. Two slips: the text cites "the homogeneity of s.3", but s.3 of S15_PROOF is "Inputs" (homogeneity is I1), and `x` is never defined in S15_PROOF (it is fq_sser's `x = Cf(C_0)/L = s varpi_0`). Lemma 3 is a statement about the code, not about the mathematics; the genus-2/3 matches test it |
| Lemma 4 | CONFIRMED | `T` for `k.w` is `k^-1 T` (arXiv:1303.1043 Def. 2.5); the vertices give `k^(2h_v-2+n_v)`; the product is `k^(2g-2+n)` |
| (I1) | CONFIRMED (given (H)) | UR s.4: `Q^m = diag(sigma^(c-1)) R^(0) M_0(z/sigma^4)` with `z/sigma^4 = tau zhat`. UR s.2 Prop. S identifies the symbol with the model `M_0` given (H). My engine derives the model `R` independently and matches UR's `F_0` through `t^7` |
| (I2), (F) = (H) + (N) | CONFIRMED | UR s.4: `Q D = P Q^m`; (H) is equivalent to `D` being a scalar series in `z^5`; (N) says `P_can` is the actual gauge, so `D = 1`. Precision: at genus g, (F) is used only through `z^(3g-3)` (the single-vertex leaf `z^(3g-2)` needs `P^-1_(3g-3)`). So genus 3 sees `c_5` only, and `c_10` first enters at genus 5 |
| (I3), (I4) | CONFIRMED | In an orthonormal basis, `omega_{g,n}(e~_i^n) = y_i^-(2g-2+n)`. With `y^Q = d^-1 y^m` this gives `omega^Q = d.omega^m`, consistent with `d = tau^(1/4) varpi_0` (`d = 1 - u/9 - 23u^2/567 ..`, own run) |
| s.4 chain | CONFIRMED | Each step is one of 1.3, I4, Lemma 4, Lemma 2, Lemma 1, I1, as stated. Hypothesis (U) holds (`P_0` diagonal) |
| s.5 items 1-7 | CONFIRMED with s.3 notes | Item 1's "independent confirmation of (N)": see C2 |

Literature: arXiv:1303.1043 s.2 (Def. 2.2, Prop. 2.4, Def. 2.5, Prop. 2.9) read in the PDF. arXiv:0712.0160: title and abstract checked
(classification of semisimple CohFTs); it is used only as background, and the worker cites no statement number. Nothing else is quoted.

## 2. Independent engine (referee/rf_s15_*.py)

* **Frozen model** (`rf_s15_base.Frozen`). Built only from the Frobenius algebra along the v-line: `H * H^b = beta_b H^(b+1)`,
  `beta = (1, -1/v, 1, -1, -1)`, `eta = 5 x involution`.
  * Idempotents `eps_j` with `H^a = N_a v^(c_a/5) Sum_j zeta^(ja) eps_j`; `Delta_j = zeta^(3j)`; `U'_j = -zeta^j` at `v = 1`.
  * The R-matrix comes from `(A - s) R_k + R_k' = [U', R_(k+1)]` with homogeneity `R_k ~ v^(-4k/5)` (off-diagonal part from the
    commutator, diagonal part from the next order).
  * Checks (logs/rf_s15_frozen.log), all PASS at both primes:
    * symplectic for `diag(1/Delta)` through `z^10`;
    * unit row of colour 0 = UR s.2 `F_0` through `t^7` (8 coefficients);
    * `F_2(1) = 1/1200`, `F_3(1) = 1/25200`; `F_4(1) = 1/180000` (2^89-1 only, 299 s, logs/rf_s15_frozen_g4.log);
    * `R(-z)` gives `F_2 = -1/1200` (control; `F_3` happens to be unchanged).
  * WK by the DVV recursion (spot values, including `<tau_2^6>_3 = 1225/144`); graph counts 7, 42, 379.
* **Outer data** (`rf_s15_outer`). The only project input is `natural_form/theory/ur_gauge.gauge`, read-only: the stored Theorem-U
  gauge, built from the periods.
  * Own canonical normalisation `P c(z)^-1` with `c(z) = P_11(0,z)`; c_5 read off = -43397/807975000.
  * `A^-1 = J P^-1 J`; `E_A`, with divisibility by `psi + psi'` asserted (symplecticity).
  * Leaves `d T`; frozen scaling `v^(2-2h) prod v^w(k,a)`, `w = (c_a+4k-4)/5`, asserted integral on every non-zero entry; global
    `v^(2-2g) d^(2g-2)`.
* **Sum** (`rf_s15_sum`). One Givental sum over composite graphs whose edges are red (frozen `E_B`) or blue (`(L x L) v^2 E_A`,
  `L = B^-1 Psi^-1`), with translation `z(1 - B^-1 1) + L (d v^w T)`. This is arXiv:1303.1043 Prop. 2.4 (red/blue) with Prop. 2.9.
  * Grouping by contraction of the red edges gives the flat sum per outer graph.
  * Grouping by the uncontracted graph gives the full `R_tot`-action per graph, i.e. the point frame (R7).
  * No MVertex6, no `rf_*` code of other referees, no worker module. `fq3_graphs` is imported only for the stored index order and the
    `|Aut|` labels.
* **Cost.** About 20 s per genus-3 run, 6 s per genus-2 run; all runs together about 15 CPU minutes.

## 3. Results (logs/rf_s15_cmp_89.log, _61.log, _89_ab1o3_m2o7.log; logs/rf_s15_run_*.log). k = coefficient of u^(k-4)

| check | 2^89-1 (0,0) | 2^61-1 (0,0) | 2^89-1 (1/3,-2/7) |
|---|---|---|---|
| R1 own flat = production record A3 (`g3_nat_<p>_nat.json`), 42 graphs, k <= 3 (168) | PASS | PASS | PASS |
| R1' same at k = 4, 5 (record used `M mod u^4`, NATURAL_FORM s.5) | 10/84 differ (graphs 1,2,3,5,6,7,8): consistent with the worker's "valid only k <= 3" | same | same |
| R2 own flat = worker's hatted `M_nat` table, k = 0..7 (336; 108 non-zero edge-graph regular coefficients) | PASS | PASS | PASS |
| R3 own flat = worker's flat table, k = 0..7 | PASS | PASS | PASS |
| R4 sum over graphs = stored point-frame total `total_id`, k <= 5 | PASS | PASS | PASS |
| R5 41 edge graphs regular (k <= 3 vanish); vertex polar part = (1/25200) t_c^-4 (`t_c` from ur_gauge's Conifold) | PASS | PASS | PASS |
| R6 own genus 2 = rf_g2flat.log (7 graphs x 4 orders) | PASS | PASS | PASS |
| **R7** composite regrouped by uncontracted graph = stored point-frame table `A_id` (`g3_state_<p>.json`), 252, all non-zero | **PASS** | **PASS** | **PASS** |

Rationals (CRT of both primes, display, logs/rf_s15_rat.log) agree with the worker's s15_rational.log. Examples:
* vertex `u^-4 (1/25200 - u/9000 + 407u^2/3780000 - 151u^3/3780000 + 226213751u^4/62052480000000 + ..)`;
* graph 2 `[u^0] = -202813/38782800000000`; graph 4 `-113/35437500`; graph 39 `[u^2] = 512/7688671875`;
* graphs 37, 41, 42 are 0 through `u^3`.

**Controls** (vs the hatted M_nat table, 336 coefficients; my tables vs the worker's tables of the same name):

| control | result | agrees with the worker's table? |
|---|---|---|
| X1 constant-M_0 split `A_M0` | 94/252 differ, totals agree | (stored) |
| X5 `stored` = `P_can c(z)` | 21/336 differ, graphs 1, 2, 5, 6, 7, 8, regular orders only; graph 2 shifts by 43397/517104000000 at `u^0`, graphs 1, 5 by ±43397/37231488000000 | identical |
| X2 `nodil` | 77/336 | identical |
| X3 `noT` | 88/336 | identical |
| X4 `nosign` | 105/336 | NOT identical (C3) |
| own: `zminus` (frozen `R(-z)`) | 112/336 | - |
| own: `dress` `exp(u z N) P`, `N = E_04 + E_43` (symplectic, graded, flat-regular, not a gauge) | 74/336 differ; all 41 edge graphs still regular, the vertex polar part moves (`T_P != t_c`) | - |
| R7 `stored` vs `A_id` | 10/252 | - |
| R7 `dress` vs `A_id` | 168/252 | - |

## 4. Overclaims, errors, corrections (quote, then correction)

* **C1.** S15_PROOF s.0 item 2: "336 coefficients ..., 112 of them non-zero regular coefficients of edge graphs". Correction: 108 are
  edge-graph coefficients and 4 are the vertex graph's regular coefficients (k = 4..7). Same count in the worker's own hatted table.
* **C2.** s.5 item 1: "the per-graph genus-3 match F2 is an independent graph-sum-level confirmation of (N) at `z^5`". The match does
  depend on `c_5` (X5), so it is a valid consistency check. But E4 already VERIFIES (N) directly through `zhat^24`, `u^17`, so this adds
  no new information about (N). Read: "consistent with E4 at `z^5`, the only normalisation coefficient visible at genus 3".
* **C3.** Control X4 ("edge quotient without `(-1)^t`") is not a well-defined alternative convention. Without the sign the quotient is not
  symmetric under `(psi,x) <-> (psi',y)`, so its value depends on the orientation of each edge (loops included). My implementation gives
  105/336 mismatches against the worker's 103/336, and our tables differ already at graph 2, `k = 2`. It still fails as it must, but the
  count is implementation-specific and should not be quoted as a property of the mathematics.
* **C4.** F1 "flat = stored hatted record (168 coefficients)" is weaker than it looks. On `k <= 3`, 164 of the 168 are zeros (edge-graph
  regularity) and 4 are vertex-graph coefficients. Every NON-ZERO per-graph comparison (k = 4..7) was against the worker's own hatted
  driver, `s15_hat.py`. That driver is project code plus a worker-written `tau^D` vertex, checked only against stored tables at
  `k <= 3` and in M_0 mode. Two things now back the claim:
  * R2: two independent flat engines agree with that driver;
  * R7: the same composite terms, regrouped, reproduce the independent production table `A_id` on 252 non-zero coefficients.
  The hatted `M_nat` values at `k >= 4` themselves still come from one implementation only.
* **C5.** Cosmetic slips: s.1.3's leaf formula (see s.1); Lemma 3's reference "s.3" and the undefined `x`; the u-adic convergence of the
  dilaton re-expansions is left unstated (`d(0) = 1`).
* **C6.** Scope of the open issue "genus 4 out of reach". It is a limit of MVertex6 (psi-deficit <= 6), not of the method. My engine has
  no deficit limit (F_4(1) checked in 5 min). A genus-4 flat run with outer data, however, needs a contraction-based implementation:
  brute force over colourings and decorations is estimated at more than 10^5 CPU-s. NOT done.

## 5. What survives, and what remains open

* **Theorem S15**: PROVED given (F), and given the description of the hatted code in Lemma 3. That description is CONFIRMED against
  the docstrings, and the matches pin it down at genus 2 and 3.
* **Consequence**: SYNTHESIS s.4's bracket [s.1.5] may be removed; its content is [(H) and (N)], already listed. At genus g only
  (F) through `z^(3g-3)` is used.
* **VERIFIED, genus 3**: flat = hatted with `M_nat` exact, per graph on `u^-4..u^3`:
  * 2^89-1 and 2^61-1 in frame (0,0), and 2^89-1 in frame (1/3,-2/7);
  * two independent flat engines;
  * point-frame regrouping = production `A_id` (R7).
* **Edge graphs**: all 41 genus-3 edge graphs are regular in the flat frame, and the vertex graph is `(1/25200) t_c^-4` + regular.
  VERIFIED at the same primes and frames.
* **Still open**:
  * (F) = (H) + (N) to all orders: CONJECTURAL (E4 window `zhat <= 24`, `u <= 17`);
  * a human check of Lemmas 1-4;
  * the identity at genus >= 4 (no computation);
  * an independent hatted-side engine for `k >= 4` (C4).

## 6. Files (referee/) and reproduction

| file | content |
|---|---|
| `rf_s15_base.py` | mod-p utilities, series, WK (DVV), stable graphs, canonical forms and automorphisms, the frozen model (Frobenius recursion) |
| `rf_s15_outer.py` | outer data from ur_gauge: canonical normalisation, `A^-1`, flat edges and leaves, v-weights, map into the frozen idempotent frame |
| `rf_s15_sum.py` | red/blue composite Givental sum; grouping by contraction (flat) and by uncontracted graph (point frame) |
| `rf_s15_run.py` | driver: `py -I -B rf_s15_run.py <89/61> <g> <N> <can/stored/nodil/noT/nosign/zminus/dress> [a,b]` |
| `rf_s15_cmp.py`, `rf_s15_rat.py`, `rf_s15_frozen.py` | comparisons R1-R7 and controls; CRT display; calibration of the engine |
| `logs/` | `rf_s15_run_g*_*.log` (one per run), `rf_s15_cmp_*.log`, `rf_s15_rat.log`, `rf_s15_frozen.log`, `rf_s15_frozen_g4.log` |

JSON tables are in the scratch folder `allgenus_s2_flat_ref/`. Reproduce: run the genus-3 `can` runs at both primes (N = 10), the
variants, and the genus-2 runs (N = 6), then `rf_s15_cmp.py 89`, `61`, `89 ab`. No project file outside `referee/` was modified;
no wolframscript was run; no process was left running.
