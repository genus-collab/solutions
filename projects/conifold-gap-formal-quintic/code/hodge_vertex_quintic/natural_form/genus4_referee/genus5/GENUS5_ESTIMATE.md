# Genus 5 by edge splitting: feasibility and cost estimate

Scope: `fq4_split.Split` / `genus4_split.py` (formal quintic, outer graph sum R'.(M.omega^), natural M) at genus 5,
K = 8 = 2g-2 coefficients in u, zhat-order KZ = 12.  **The genus-5 run itself was NOT attempted.**  Everything here is timing,
counting and small consistency checks; all scripts and logs are in this folder, nothing outside it was modified
(`fq4_split` is imported and two of its module names are rebound at run time, see `g5_lib.py`).
Labels: VERIFIED (computed, with prime and range), ARGUED, CONJECTURAL.
"Reference machine" = this workstation (20 cores / 28 threads, CPython 3.14.6), one core, with at most 6 of my processes running.

## 0. Bottom line

| quantity (per prime, natural M, K = 8) | value | status |
|---|---|---|
| leg-free genus-5 stable graphs | 4555 = (1, 3, 11, 34, 100, 239, 492, 784, 1002, 926, 632, 260, 71) for E = 0..12 | VERIFIED |
| original generator `fq3_graphs.legged_graphs(5, 0, 12)` | not finished in 24 min (E <= 11 done in 1379 s; E = 12 about 450 s more) | VERIFIED / extrapolated |
| fast canonical form (`g5_lib.py`) | 10 s for genus 5, same 379 genus-4 graphs and \|Aut\| | VERIFIED |
| `Split(p, 5, (0,0), 8, 'nat')` set-up | 11 s (2^89-1), 4 s (2^61-1), 51 MB, all asserts pass | VERIFIED |
| `Split.graph` as it stands, all 4555 inner graphs | **about 277 CPU hours** (250-290), slowest graph about 4.1 h, 296 graphs above 10 min | estimate, validated on 30 complete graphs |
| the same with one elimination for all subsets (prototype `g5_dp.py`) | **about 46 CPU hours**, slowest graph about 0.7 h, 24 graphs above 10 min | estimate, validated on 17 complete graphs |
| variant `id` (point frame, total only) | about 1.3 CPU hours, slowest graph 20 s | estimate, validated on 40 complete graphs |
| memory per process | 51 MB (`Split.graph`), below 100 MB (subset-DP) | VERIFIED on the samples |

Recommendation: **cluster** for the natural-M run (either code path is far beyond 30 minutes locally: 277 h or 46 h of CPU against
6 processes); only the `id` variant (the total gap, which fixes c_5) is small enough for a local run of about 15-20 minutes wall
with 6 processes, and only after the graph generator is replaced.  Details in section 6.

## 1. Graph generation

* `g5_gen_time.py` (the level loop of `fq3_graphs.legged_keys` with a 24-minute deadline, log `g5_gen_time.log`): levels E = 0..11
  complete after 1379 s (per level 0.4 s at E = 5, 41 s at E = 7, 136 s, 290 s, 338 s, 569 s at E = 8..11), stopped inside E = 12
  after 35 of 260 source graphs (61 s).  Counts of the complete levels: (1, 3, 11, 34, 100, 239, 492, 784, 1002, 926, 632, 260),
  sum 4484; extrapolated total about 30.5 minutes.  The cost is `canon` (minimum over all V! vertex orders; V = 8 gives 40320 orders
  per call).  **This is paid by every process** (`Split.__init__` calls `legged_graphs`), i.e. by every pool worker and every array task.
* `g5_lib.canon_count` (individualisation-refinement canonical form; the number of minimal leaves gives the vertex automorphisms):
  genus 5 in 10 s.  `g5_gen_check.log`, all PASS:
  - types (2,0), (3,0), (4,0), (2,2), (1,4), (3,1), (0,6): the fast list, pushed through the ORIGINAL `canon`, equals the original
    list with equal |Aut| graph by graph (genus 4: 379 graphs); invariance under random relabellings;
  - genus 5: 4555 graphs (the known number), per E as in the table, the E <= 11 counts equal those of the original generator;
    per V = 1..8: (6, 55, 262, 743, 1252, 1298, 726, 213); |Aut| of all 4555 graphs equals the brute-force `fq_graphs.aut_count`.
  The fast keys are a different representative of each class than `fq3_graphs.canon`; `Split` only needs consistency between the
  graph list and `contract`, so both names are rebound together (`g5_lib.make_split`).  `g5_setup.log` (a): with the fast list
  `Split.graph` reproduces the stored genus-4 per-inner-graph tables (2^89-1, natural M) for the 36 inner graphs tried.
* The other graph costs are harmless at genus 5 (`g5_misc.log`): `fq3_graphs._graph_dict` (`autos_of`) on all 4555 graphs 11 s;
  the contraction classes of the 3302 subsets of an E = 12, V = 8 graph 0.4 s with the original `canon` (0.1 s with the fast one).

## 2. Set-up at genus 5 (`g5_setup.log`, `g5_misc.log`)

`Split(p, 5, (0, 0), K = 8, 'nat')` with rows at Nz = 12: all assertions of `Split.__init__` pass at both primes (symbol valuations,
hatted-row purity `H.val_ok`, `H.min_prec >= 8`, M_0 M_0^-1 = 1 through zhat^12); `frozen_ok` = True ([u^0] E_R(a,b) = Eh_0(a,b) for
the 78 keys a + b <= 11).  Set-up 10.6 s / 4.3 s, peak memory 51 MB, all 29 vertex tables 0.8 s.  VERIFIED, 2^89-1 and 2^61-1.
* **u-valuation of E_out(a, b) is exactly 1 for all 78 keys and all 1950 colour entries** (both primes).  So the valuation pruning
  of `Split.graph` (`sum(vals) >= K`) never removes anything beyond the popcount restriction |S| < K.
* decorations per vertex (ordered half-edge exponents): (0,3): 1, (0,4): 5, (0,5): 21, (0,6): 84, (0,7): 330, (0,8): 1287,
  (0,9): 5005, (0,10): 19448, (1,4): 70, (2,3): 84, (3,4): 1001, (4,2): 78, ...
* the collect step works at g = 5: `tc_series` / `tc_expansion3` return (c, 0, ..., 0) on a synthetic total c t_c^-8.

## 3. Cost of `Split.graph` as it stands

Method.  For an inner graph the work is ncombo x nsub colour sums (ncombo = product of the decoration counts of its vertices,
nsub = number of subsets with |S| < 8).  `g5_sample.py` times random (decoration, subset) pairs and returns
ncombo x nsub x (mean time) + preparation; no full graph is computed.  `--timing` of `genus4_split.py` cannot be used at genus 5
because of section 1.
* Calibration at genus 4 (`g5_calib4.py`, K = 6): sampled total 1488 s; 15 graphs run completely under the same conditions:
  measured/estimate between 0.67 and 1.08.  The stored genus-4 run (13 processes) shows 3251 s CPU: a roughly uniform factor 2.2
  (1.7-2.6 for every E >= 3), i.e. machine load, not the estimator (ARGUED: the stored per-graph times also contain the set-up of
  each worker).  So genus 4 costs about 0.4 CPU hours on the reference machine.
* Variant `id` at genus 5 (`g5_sample_g5_89_id_K8.jsonl`): 1.34 CPU hours in total; 40 graphs run completely (the 10 heaviest and
  30 others): measured 241 s against 243 s estimated.
* Genus 5: all 4555 graphs sampled (510860 pairs), the 150 heaviest re-sampled independently (45141 pairs).  30 graphs, E = 0..12,
  each below 10 minutes, run completely (`g5_full_*.jsonl`): measured 4299 s against 4762 s estimated (ratio 0.90; per graph 0.60-1.25).
* Both primes cost the same (26 graphs, identical pairs: ratio 2^61-1 / 2^89-1 = 1.06).

| E | graphs | sum [CPU h] | mean [s] | max [h] |   | E | graphs | sum [CPU h] | mean [s] | max [h] |
|---|---|---|---|---|---|---|---|---|---|---|
| 0-3 | 49 | 0.05 | 4 | 0.01 |   | 8 | 1002 | 73.9 | 266 | 1.70 |
| 4 | 100 | 0.6 | 23 | 0.05 |   | 9 | 926 | 76.0 | 296 | 2.38 |
| 5 | 239 | 3.9 | 59 | 0.12 |   | 10 | 632 | 49.7 | 283 | 4.09 |
| 6 | 492 | 14.7 | 107 | 0.15 |   | 11 | 260 | 12.3 | 171 | 0.57 |
| 7 | 784 | 43.5 | 200 | 0.69 |   | 12 | 71 | 1.9 | 97 | 0.21 |

Total **277 CPU hours** (first pass alone 268).  296 graphs above 10 min (46 % of the total), 72 above 30 min, 15 above 1 h, 5 above 2 h.
Scaling from genus 4: 180 times more (decoration, subset) pairs (1.58e8 against 8.8e5), about 650 times more time.
Slowest graphs (all vertices of genus 0 unless said; dense graphs with one or more high-valence vertices, not the trivalent ones):
* E = 10, V = 6, valences (3,3,3,3,3,5), edges 01 02 05 13 15 24 25 34 35 45: 4.1 h (21 decorations x 968 subsets x 0.72 s);
* E = 9, V = 5, valences (3,3,4,4,4) (two graphs): 2.4 h and 2.3 h; (3,3,4,5,3): 2.2 h; (3,3,3,3,6): 2.0 h;
* E = 8, V = 5, g_v = (0,0,0,0,1), valences (3,3,4,4,2): 1.7 h; then eight E = 10, V = 6 graphs at 1.3-1.6 h.
The 71 trivalent graphs (E = 12, one decoration, 3302 subsets) are cheap: 1.9 h together.
Cost versus (E, V, valences): cost = ncombo x nsub x (one colour sum), the colour sum scaling like Sum over elimination steps of
5^(frontier + 1) x (1 + edges closed).  This structural proxy has log-log correlation 0.983 with the sampled estimates (slope 1.00,
scatter factor 1.26).  The proxy of `genus4_split.cost_order` has correlation -0.03 at genus 5 (its factor (12 - E)^2 was tuned to genus 4).

## 4. What has to change in the code for genus 5

Necessary:
1. **Graph generation** (section 1): replace `canon` by a refinement-based canonical form (or load a cached list); otherwise every
   process spends about half an hour (more on slower cores) before its first graph.
2. **Task balancing**: with the present proxy and round robin the heaviest of N = 100 tasks is 7.5 h against a mean of 2.8 h
   (N = 200: 5.4 h against 1.4 h).  Longest-first assignment with the estimates of `g5_cost_table.json` (per graph: gv, edges, |Aut|,
   estimated seconds for both code paths) gives 4.1 h at N = 100 (the single slowest graph) and a perfectly flat load at N <= 50.
3. **Collect step**: `checks()` knows c_g only for g = 2, 3, 4, so at genus 5 the GAP and VERTEX GRAPH checks are silently skipped
   (`cg is not None`).  The pattern c_g = |B_2g| / (2g (2g-2) 5^(g-1)) of the three known values gives c_5 = 1/660000 (CONJECTURAL;
   the run would decide it); at least the vanishing of the t_c^-7 .. t_c^-1 coefficients should be checked unconditionally.
4. **Log file**: in `--task` mode all array tasks append to the same `logs/genus4_split_<tag>.log`; use one log per task.

Worth doing (measured):
5. **One elimination for all outer-edge subsets** (`g5_dp.py`, prototype).  `Split.graph` repeats the colour sum for every subset
   (255 at E = 8, 502, 968, 1816, 3302 at E = 9..12).  Carrying the inner/outer marks of the closed edges in the elimination state
   costs Sum_steps 2^(closed edges) 5^(frontier+1) instead of 2^E Sum_steps 5^(frontier+1).
   - VERIFIED: `graph_dp` equals the stored genus-4 tables of `genus4_split` for ALL 379 inner graphs (2^89-1, natural M), and equals
     `Split.graph` on 18 genus-5 inner graphs with E = 0..12 (2^89-1, K = 8); `g5_dp.log`.
   - genus 4: 489 s CPU (against about 1490 s); genus 5 on the 18 graphs: 125 s against 628 s, per graph 1.5-5.4 times for 2 <= E <= 9 and
     4.4-17 times for E >= 10 (slower than `Split.graph` only for E <= 1, where both take under 0.1 s).  Sampled over all 4555 graphs (120 heaviest re-sampled, 17 run completely: ratio 0.92):
     **46 CPU hours, slowest graph 0.7 h (the same E = 10 graph), 24 graphs above 10 min, peak memory 99 MB**.
   - ARGUED, not measured: an elimination order chosen for this cost (close edges late, keep the frontier small) should gain more.
6. **Precision K - |S|**: every outer edge carries a factor u, so a term with |S| outer edges needs only K - |S| coefficients.
   Prototype `g5_trunc.py`: identical results on 360 random pairs, colour sums 1.6-1.7 times faster (cost-weighted).  The present
   `mul` already skips leading zeros, which is why the gain is modest; inside the subset-DP it is not additive with item 5.

Not a problem (measured):
7. **Subset restriction** |S| < K: keeps 3302 of 4096 subsets at E = 12, 1816 of 2048 at E = 11, 968 of 1024, 502 of 512, 255 of 256
   at E = 10, 9, 8, where 72 % of the time is spent; it cannot be tightened by valuations (section 2).
8. **`itertools.permutations`** for high-valence vertices: the worst case is the one-vertex graph with five loops (valence 10,
   3.6e6 orders for each of 45 exponent partitions): 12 s; all 4555 preparations together 128 s.
9. Memory, the part-file format, the collect arithmetic (`tc_expansion3` at g = 5) need no change.

## 5. Python version

Static check of the 23 project modules reachable from `fq4_split` / `genus4_split`: no syntax newer than Python 3.6 (no walrus,
match, positional-only arguments; `ast.parse` with feature_version 3.6-3.8 accepts all), library calls needing **Python 3.8**:
`math.prod` (fq4_split), `math.comb` (fq_rows), `math.isqrt` (fq_fp).  Nothing from 3.9 or later (`self.cache` are plain dicts,
`lru_cache(maxsize=None)` is old).  So 3.9-3.11 are fine; 3.6/3.7 are not.  The scripts of this folder use nothing newer either.
Not tested on an older interpreter (only 3.14 is installed here): ARGUED from the static check.  No data file is read on the run
path; the hard-coded Windows path in `fq_graphs.DATADIR` is used only by `--regress`.  The launcher `py -I` becomes `python3 -I`.
Speed: CPython 3.11+ is usually noticeably faster than 3.9/3.10 on this kind of pure-Python integer loop, and the cluster cores are
slower per thread than this workstation; neither factor was measured (CONJECTURAL: 1.5-3 together).  Calibrate with one short job:
inner graphs 808, 1109, 4511 of the fast list took 71 s, 96 s, 31 s here with `Split.graph` (11 s, 21 s, 3.5 s with `graph_dp`).

## 6. Recommendation

Stage 1, local (after change 1): variant `id`, both primes: about 1.3 CPU hours each, no graph above 20 s, about 15-20 minutes
wall with 6 processes.  It gives the total, i.e. the gap test and c_5.  (A natural-M run restricted to S = {} would give the vertex
graph alone, at a similar cost by the same count of one colour sum per decoration (ARGUED, not timed, and not an existing mode);
with the total it tests regularity of the SUM of the 4554 outer graphs with edges, not of each.)

Stage 2, cluster (per-outer-graph regularity; one array per prime; safety factor 3 on the reference times):

| code path | CPU on the reference machine | array | heaviest task (reference) | request per task |
|---|---|---|---|---|
| subset-DP (changes 1-5) | 46 h | 50 tasks, 1 core | 0.9 h longest-first; 1.7 h present round robin | `--time=08:00:00 --mem=1G` |
| as it stands + changes 1-4 | 277 h | 100 tasks, 1 core | 4.1 h longest-first | `--time=24:00:00 --mem=1G` |
| as it stands + change 1 only | 277 h | 100 tasks, 1 core | 7.5 h present round robin | `--time=48:00:00 --mem=1G` |

Core hours on the cluster: about 46 x f or 277 x f per (prime, variant) with f the unmeasured machine factor (1.5-3).  Each control
used at genus 4 (bare M_0, shifted tau, second frame) is one more run of the same size.  Tasks are resumable graph by graph
(part file rewritten after each inner graph), so a task that hits its limit is simply resubmitted.
I recommend the subset-DP path, after it is wired into the driver (`--task` / `--collect`) and re-validated by `--regress`.

```bash
#!/bin/bash
#SBATCH --job-name=fq5_89_nat
#SBATCH --array=0-49
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=08:00:00
#SBATCH --output=logs/slurm_%A_%a.out
# Python >= 3.8 required; <driver> = the genus-5 driver with changes 1-5
python3 -I <driver>.py --genus 5 --prime 89 --variant nat --task $SLURM_ARRAY_TASK_ID --ntasks 50
```

## 7. Caveats

* All times are estimates from sampling, checked against 30 (resp. 17) complete graphs; the per-graph values of the heaviest graphs
  rest on 22-115 sampled pairs (4-9 sampled decorations for the subset-DP) and can be off by 30 %.  The totals are much better determined.
* Times were taken with 5-6 of my processes on a machine that was also running other jobs (the genus-4 controls), on a CPU with
  fast and slow cores; the genus-4 comparison shows that 13 simultaneous processes cost a factor 2.2.
* The subset-DP is a prototype: functions only (`g5_dp.graph_dp`), natural-M type variants only, no driver, validated as stated in 4.5.
* Nothing here says anything about the outcome at genus 5: no genus-5 total, no regularity statement, no value of c_5 was computed.

## 8. Files (this folder)

`g5_lib.py` fast canonical form, generator, patched Split, cost sampler; `g5_gen_time.py`, `g5_gen_check.py` (+ `.log`) section 1;
`g5_setup.py`, `g5_misc.py` (+ `.log`) section 2; `g5_sample.py`, `g5_calib4.py`, `g5_analyse.py`, `g5_analyse.log`,
`g5_sample_*.jsonl`, `g5_resample_*.jsonl`, `g5_full_*.jsonl` section 3; `g5_dp.py`, `g5_dp_run.py`, `g5_dp.log`, `g5_dpsample_*`,
`g5_dpresample_*`, `g5_dpfull_*`, `g5_trunc.py`, `g5_trunc.log` section 4; `g5_cost_table.json` per-graph estimates;
`g5_graphs_4_0_9.json`, `g5_graphs_5_0_12.json` cached graph lists; `*_console.log`, `top_dp.txt`, `full_dp.txt` run by-products.
