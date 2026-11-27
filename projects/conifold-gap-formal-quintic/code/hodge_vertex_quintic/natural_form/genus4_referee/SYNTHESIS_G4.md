# Genus 4 of the formal quintic by edge splitting: synthesis of the referee round

Date 2026-11-11. Written by an AI agent from three AI reports (`review/REVIEW.md`, `indep/INDEP.md`, `genus5/GENUS5_ESTIMATE.md`),
the control logs `logs/genus4_split_g4_*_console.log` and a read-only cross-check of the stored tables (`syn_check.py`, `syn_check.log`,
this folder). No human has checked any of it. Labels: VERIFIED (exact finite computation; prime and orders stated), ARGUED, CONJECTURAL.
Finite computations are evidence, not proofs. All statements are modulo 2^89-1 and/or 2^61-1, `zhat <= 9`, `u^0..u^5` (K = 6).

## 1. What stands at genus 4

* **S1 (VERIFIED, both primes, frame (0,0)).** Point-frame total `= (1/180000) t_c^-6 + 0 t_c^-5 + ... + 0 t_c^-1`.
Support: project code; independent code (`indep/`, all 379 graph terms equal coefficient by coefficient); own `t_c` arithmetic in `review/rv_tc.py`.
* **S2 (VERIFIED, both primes, frame (0,0)).** With `M = M_0(tau(u) zhat)` mod `(zhat^10, u^6)`, `1/tau = 2F1(1,4/5;9/5;u)`, the regrouped split
sum has: vertex graph `= (1/180000)(u/t_c)^6 mod u^6`, i.e. `1/180000 - 7u/300000 + 701u^2/18000000 - 1459u^3/45000000 +
923639u^4/67500000000 - 171583u^5/67500000000`, and all 378 outer graphs with edges `O(u^6)`.
Support: (i) project code; (ii) independent re-implementation sharing no module with it (own tables from the earlier referee rows, own
Witten-Kontsevich numbers, Mathematica graph list with own automorphism counts, colour sum in momentum space checked against `5^V` brute
force): all 379 outer rows and all 7146 non-zero pieces (inner graph -> outer graph) equal; (iii) a full run of an unmodified copy in the
review reproduces the table of record (69 rows).
* **S3 (VERIFIED, both primes, project code only).** Same statement in the second conifold frame `(a,b) = (1/3, -2/7)`: all checks pass
(`genus4_split_g4_{89,61}_nat_ab_console.log`); total `= 1/180000 - 11u/900000 - 79u^2/378000000 + 20081u^3/2126250000 + ..`, the rationals
printed at 2^89-1 reduce to the residues at 2^61-1 (`syn_check.py` (5)). Not re-run by the independent code; the review checked only the
stored tables and the composition law at 2^89-1.
* **S4 (VERIFIED, 2^89-1 and 2^61-1 in frame (0,0), 2^89-1 in the second frame; `review/rv_compose.py`, own series arithmetic).** The tables
used are what they are said to be: `M` symplectic, `E_M = (1 - M^-1 M^-T)/(z+w)`, `E_out = M^-1 E_R' M^-T` with `R' = R M^-1`,
`T_R = T_M + M^-1 T_R'`, every `E_out` entry of `u`-valuation exactly 1, and `Psi^-1 R'` obeys the natural-form support rule through
`u^5`, `zhat^9` (which fails for a wrong `tau_4` or `tau_5`).
* **S5 (ARGUED; numerically supported at genus 2, 3 only).** The regrouped split sum IS the outer graph sum `R'.(M.omega^)`, graph by graph.
Written argument in `review/REVIEW.md` s.2 (composition law, orbit-stabiliser, exact sequence for `Aut(Gt,S)`, leg rule, `xhat` powers).
Its combinatorial core (weights `t^(inner edges)`, exact rationals) is VERIFIED for all 42 genus-3 and all 379 genus-4 outer graphs
(38811 pairs). The leg/leaf step is tested only by the genus-2/3 regression against the old vertex-tensor code. Both split codes compute
the same regrouped sum, so the independent run cannot detect an error in this identification at genus 4; the argument is genus independent.

## 2. How the claim should be worded (overclaims)

* "All 378 outer graphs with edges are regular" has less content than it sounds (VERIFIED, same counts at both primes, two agents):
  **229** graphs (6 to 9 edges) are `O(u^6)` trivially, since every outer edge is `O(u)`, and are never computed; **81** graphs
  (3 + 19 + 59 with 3, 4, 5 edges) receive only pieces that vanish identically (a selection rule of the true `E_out`, observed, not
  derived; a generic `E_out` gives non-zero pieces, and with the bare `M_0` the same graphs are non-zero); **68** graphs
  (3, 7, 18, 24, 16 with 1..5 edges) are regular by cancellation among 6767 non-zero pieces. Quote 68 (or 149), not 378.
* "REGULARITY" and "VERTEX GRAPH carries the polar part" are one check, not two (vertex = total - edge graphs).
* GAP, the per-inner-graph additivity, coverage and `frozen_ok` do NOT test the `u`-dependence of `E_M` (`E_out := E_R - E_M` telescopes).
  Mutants with a wrong `tau_4`, `tau_5` or a wrong tau power at `zhat >= 7` pass GAP and the whole genus-2/3 regression (VERIFIED,
  `review/out/mut_g3.log`); they are caught only by genus-4 REGULARITY and by S4. So "the regression validates the code" holds for the
  algorithm, the conventions, `tau_1..tau_3` and `zhat <= 6` only.
* Sensitivity is not uniqueness: regularity detects each of `tau_4`, `tau_5`; whether it determines them was not examined.

## 3. Bugs and weaknesses (none changed a stored number)

| id | where | defect | state |
|---|---|---|---|
| B1 | `genus4_split.py --collect` (08:44 version) | part files chosen by prefix `tag + "_"`: collecting `g4_89_nat` after the variants existed would have summed 2-3 runs and overwritten the table of record | fixed 09:40:48 (exact regex + tag assertion), re-tested |
| B2 | `merge` (08:44) | no overlap check between part files (matters for `--task` arrays) | fixed (identical pieces asserted, counted once) |
| W1 | `Split.__init__` (08:51) | `min valO >= 1`, needed for `popcount < K`, not asserted | fixed (assertion + check line) |
| W2 | `checks` (08:44) | GAP/VERTEX did not require six polar coefficients | fixed |
| W4 | regression | cannot see `tau_4`, `tau_5`, `zhat^7..9` | OPEN: add `rv_compose.py` (12 s) or an equivalent to the driver's checks |
| G1 | `checks` at `g = 5` | `c_g` known only for g = 2,3,4, so GAP and VERTEX are silently skipped | OPEN |
| G2 | `--task` mode | all array tasks append to one log; array index must be 0-based; mode read but never run | OPEN |

No bug was found in `fq4_split.py` (26 mutants: 23 caught by the genus-2/3 regression, the other 3 by the genus-4 checks; unit tests of
`contract`, `vertex_table` and the helpers against brute force). Process notes: the two project files were edited by their author during
the review (some mutant runs used the earlier source, identical apart from the new assertion); one agent deleted and regenerated three
of its own logs; one may have refreshed `__pycache__` of project modules. Longest referee run: 24 minutes (genus-5 generator timing).

## 4. What the controls show (project runs, 2^89-1 only; both reproduced by the independent code, all rows and pieces)

| control | GAP | REGULARITY | detail |
|---|---|---|---|
| bare symbol `M = M_0` | PASS | FAIL as expected | all 149 outer graphs with 1..5 edges non-regular, each with lowest order exactly `u^E` (3, 7, 21, 43, 75); vertex graph has `t_c^-5..t_c^-1` terms |
| `tau_5 -> tau_5 + 1` | PASS | FAIL as expected | exactly the three one-edge graphs, at `u^5` only; vertex graph acquires a `t_c^-1` term |
| second frame `(1/3, -2/7)`, both primes | PASS | PASS | S3 |

* The regularity check discriminates (no dead code path); the totals of nat, id, `M_0`, `nat_t5` are identical (VERIFIED, `syn_check.py` (4)).
* New consistency check between the two controls (VERIFIED at 2^89-1, `syn_check.py` (3)): for each of the three one-edge graphs
  `[u^5](tau_5 shifted by 1) = (9/4) [u^1](bare M_0)`. Both are the first-order response to a change of tau (`tau = 1 - 4u/9 + ..`, so
  bare `M_0` is `delta tau_1 = +4/9`): the same linear functional appears at `u^1` and `u^5` (ARGUED reading).

## 5. Not established

* S5 beyond the written argument; anything beyond `u^5`, two primes, two frames, genus 4. Conjecture N for `g >= 5`: CONJECTURAL.
* Shared inputs not rebuilt from scratch: rows (two derivations exist, neither by the independent agent), `t_c(u)` (not compared with
  the literature), symbol matrices, graph enumeration. The gap total with five vanishing coefficients is a strong joint test of them.
* Controls and all genus-4 mutant runs at 2^89-1 only; second frame without an independent re-run.

## 6. Long-run rule and SLURM, genus 4

No genus-4 run comes near 30 minutes: about 1450 CPU-seconds on an idle machine (5 minutes wall with 5 processes; 4.2-5.1 minutes wall
in the control logs); **no Stanage submission is needed at genus 4 and no SLURM script was prepared.** Further frames or controls cost the same.

## 7. Genus 5: cost estimate and recommendation (nothing about the genus-5 outcome was computed)

Facts (VERIFIED unless said): 4555 inner graphs; K = 8, `zhat <= 12`; set-up passes all assertions at both primes (11 s, 51 MB);
`E_out` has valuation exactly 1 on all 78 keys, so 2891 outer graphs (8..12 edges) are trivially regular and 1663 (1..7 edges) carry the
content (ARGUED from the valuation). Blocker: `fq3_graphs.legged_graphs(5,0,12)` needs about 30 minutes in EVERY process (canonical form
over `V!` orders); a refinement-based canonical form (`genus5/g5_lib.py`) does it in 10 s and reproduces the 379 genus-4 graphs.

| run (per prime, per variant) | CPU on this workstation | slowest inner graph | where |
|---|---|---|---|
| `id` (total only: gap test, value of `c_5`) | about 1.3 h (sampling estimate, 40 complete graphs) | 20 s | local, 15-20 min wall on 6 processes: close to the limit, use a resumable run with a deadline, or the cluster |
| natural `M`, `Split.graph` as it stands | about 277 h (sampling estimate, 30 complete graphs, ratio 0.90) | 4.1 h | cluster only |
| natural `M`, one elimination for all subsets (`genus5/g5_dp.py`, prototype; equals the stored genus-4 tables for all 379 graphs at 2^89-1) | about 46 h (17 complete graphs, ratio 0.92) | 0.7 h | cluster only |

Cluster factor f = 1.5-3 (CONJECTURAL, unmeasured). Expected `c_5 = |B_10|/(10 * 8 * 5^4) = 1/660000` (CONJECTURAL; pattern of g = 2, 3, 4).

Recommendation, in order:
1. Before any genus-5 run, in NEW files (a genus-5 driver): fast generator or cached graph list; subset-DP wired into `--task`/`--collect`
   and re-validated by `--regress` plus a genus-4 full run at both primes; task assignment longest-first from `genus5/g5_cost_table.json`;
   `c_5` and an unconditional test of the `t_c^-7..t_c^-1` coefficients in `checks` (G1); one log per task (G2); the composition/support
   test of S4 at `zhat <= 12`, `u <= 7` (W4), which is cheap and is the only check of `tau_6`, `tau_7` besides the claim itself.
2. Stage 1 (minutes): `id` variant at both primes: fixes `c_5` and tests the shared inputs at genus 5.
3. Stage 2 (Stanage; this is where the SLURM rule applies, each inner graph is independent): one array per (prime, variant), 50 one-core
   tasks, `--time=08:00:00 --mem=1G`, 0-based array index, Python >= 3.8 (static check only), resumable graph by graph. Natural `M` at both
   primes first; then the controls (bare `M_0`, a shifted `tau_k` for k = 6 or 7, second frame), each one more run of the same size:
   about 46 x f core hours each, 5 runs about 230 x f; a calibration job on three named graphs settles f first. **The SLURM package for
   genus 5 does not exist yet** (only a sketch in `genus5/GENUS5_ESTIMATE.md` s.6); prepare it (tarball, copy-paste commands) after step 1.
4. Cheaper first verdict (ARGUED, not timed): subsets with at most two outer edges only; exact for the vertex graph and the 14 outer graphs
   with one or two edges (the ones sensitive to the highest `tau_k`), with roughly ten times fewer colour sums on the expensive graphs.
