# Genus-5 cluster evidence: full census and audit of the copied-back part files (2026-11-13)

What this is: the audit that could not be done while only `cluster_genus5_quintic/results/*.json` was on record. The cluster's
logs and per-graph pieces have now been copied back (`g5_evidence.tar.gz`, 46 MB, 117 MB extracted, `parts/` 113 MB), so the
result tables can be rebuilt from the pieces by independent code, the task decomposition can be recomputed, the cancellation
census can be taken from the cluster's own numbers rather than from a partial local re-run, and the per-inner-graph additivity
no longer has to be taken from the collect log.

Everything below is arithmetic on stored residues (plus one small new engine run, s.6). Wall time of the whole census: about
8 minutes. Code: `cen_lib.py` (readers, my own merge, my own longest-processing-time assignment), `cen_merge.py` (s.1),
`cen_prov.py` (s.2), `cen_census.py` (s.3-5), `cen_report.py` (assembles `census_g5.json`). None of the `g5_*`/`fq_*` modules
is imported, so the merge, the assignment and the additivity are not the package's own code; the *pieces* are of course the
package's output, so this is an audit of the stored evidence, not independence from `g5_engine`.

Nomenclature as in `NATURAL_FORM.md` s.5-6: a *piece* is a pair (inner graph `t`, outer graph `q`), `q` the graph obtained by
contracting the inner edges of `t`; `|S| = E(q)` is the number of outer edges; a genus-5 term is `u^-8 G(u)` and a row is the
eight polar coefficients `[u^0..u^7]`. A piece is stored only when some residue is non-zero.

## 1. The result tables rebuilt from the part files

For each of the seven runs I read every part file, took the union of the per-inner-graph pieces with my own code, summed them
into outer-graph rows modulo the run's prime and compared **every row and every residue** with `results/<tag>.json`.

| tag | part files | inner graphs | missing | duplicated | stored rows | rows non-zero | non-zero pieces | CPU | rows/residues equal |
|---|---|---|---|---|---|---|---|---|---|
| `g5_89_nat` | 100 | 4555 | 0 | 0 | 628 | 1 (the vertex row) | 340363 | 51.2 h | yes |
| `g5_61_nat` | 100 | 4555 | 0 | 0 | 628 | 1 | 340363 | 53.3 h | yes |
| `g5_89_nat_ab1o3_m2o7` | 100 | 4555 | 0 | 0 | 628 | 1 | 340369 | 52.6 h | yes |
| `g5_89_M0_out2` (control) | 20 | 4555 | 0 | 0 | 15 | 15 | 39670 | 14.7 h | yes |
| `g5_89_nat_t7_out2` (control) | 20 | 4555 | 0 | 0 | 15 | 4 | 39670 | 18.6 h | yes |
| `g5_89_id` (workstation) | 1 | 4555 | 0 | 0 | 4555 | 4555 | 4555 | 1.8 h | yes |
| `g5_61_id` (workstation) | 1 | 4555 | 0 | 0 | 4555 | 4555 | 4555 | 1.8 h | yes |

**No discrepancy of any kind.** 0 rows present on one side only, 0 differing residues, in all seven runs; the row sums also
reproduce the `info["total"]` stored in each result file, and the `outer` graph list (`E`, `gv`, `aut`, index by index) is
identical in all seven files. Of the 628 stored rows of a natural run, 627 are identically zero and one - the vertex row - is
non-zero: that is the regularity statement, and it is visible in the part files themselves.

Two cross-checks fall out: the two `--maxout 2` controls have exactly 39670 non-zero pieces each, which is precisely the
`|S| <= 2` part of the full census below (4555 + 10300 + 24815), and the `g5_89_nat_t7_out2` control has 4 non-zero rows (the
vertex row and the three one-edge rows), the `g5_89_M0_out2` control all 15 it computes, matching the `failed` lists.

## 2. Provenance and the task decomposition

The result files carry no host, date, CPU time or part-file inventory. From the logs:

| tag | SLURM array | collect | tasks | hosts | wall per task (min / mean / max) | CPU (sum of `times`) |
|---|---|---|---|---|---|---|
| `g5_89_nat` | 12149822 `[0..99]` | 12149823 | 100 | node003-009, 011, 013 | 1716 / 1849 / 3036 s | 51.2 h |
| `g5_61_nat` | 12149824 `[0..99]` | 12149825 | 100 | node003-009, 013 | 1780 / 1924 / 3147 s | 53.3 h |
| `g5_89_nat_ab1o3_m2o7` | 12156906 `[0..99]` | 12156907 | 100 | node003-006 | 1756 / 1900 / 3025 s | 52.6 h |
| `g5_89_M0_out2` | 12156902 `[0..19]` | 12156903 | 20 | node003-005 | 2427 / 2659 / 2888 s | 14.7 h |
| `g5_89_nat_t7_out2` | 12156904 `[0..19]` | 12156905 | 20 | node005 | 3108 / 3348 / 3570 s | 18.6 h |
| `g5_89_id`, `g5_61_id` | - (workstation, 13 processes) | - | 1 | - | wall 651 s / 653 s | 1.8 h each |

All cluster tasks ran under **Python 3.9.25** on `*.stanage.hpc.shef.ac.uk` (the workstation `id` runs under Python 3.14.6);
validation job 12149764 on node008, `machine factor f = 0.94` (cluster time / reference time), calibration
`{808: 9.8 s (ref 10.8 s), 1109: 19.6 s (ref 20.9 s), 4511: 3.6 s (ref 3.5 s)}`, 31 PASS / 0 FAIL, wall 377 s. (The
workstation's own `g5_validate`, in `logs/workstation/`, is 32 PASS / 0 FAIL, wall 490 s, estimate 39 CPU h per run.)

**The decomposition is exactly what the code says it is.** I recomputed the longest-processing-time assignment from
`data/cost_table_g5.json` (column `est_dp`, costs `max(est, 0.01)`, decreasing cost with ties by index, each graph to the
least-loaded task with ties by task index) with my own implementation, and compared it slice by slice with the inner graphs
actually stored in `task<I>of<N>.json`:

* 100-task runs: all 100 slices agree exactly (estimated heaviest task 0.704 h, mean 0.462 h on the reference machine;
  observed mean 0.514 h, heaviest 0.843 h, consistent with `f = 0.94` and a cold cache);
* 20-task runs: all 20 slices agree exactly (estimated 2.31 h per task);
* every inner graph lies in **exactly one** task file in every run (0 duplicates, so the disjointness question raised by the
  genus-4 referee is settled on the real files), the union is exactly `0..4554`, and no part file carries a foreign tag;
* **no task was resumed or restarted**: no `resuming: N done` line anywhere, exactly one header line per task log, and exactly
  4555 `inner graph ...` lines per run - one per graph, so nothing was computed twice or inherited from an earlier attempt;
* every task log ends with `TASK I OF N COMPLETE`, set-up 6-7 s, and the per-task `COMPLETE` wall times sum to 51.4 h / 53.4 h /
  52.8 h / 14.8 h / 18.6 h, matching the `times` fields (the difference is the set-up).

This is the provenance the result files lack; it is also recorded in `census_g5.json` under `step2_provenance`.

## 3. The cancellation census, from the cluster's own pieces

`g5_89_nat`: **340363 non-zero pieces**, 4555 of them on the vertex row and 335808 on the 627 cancelling rows.
`g5_61_nat`: the same 340363. `g5_89_nat_ab1o3_m2o7` (second frame): 340369.

By number of outer edges `|S| = E(q)` and by the edge number of the *inner* graph `t`:

| `\|S\|` | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 |
|---|---|---|---|---|---|---|---|---|
| pieces | 4555 | 10300 | 24815 | 49315 | 82393 | 91858 | 60835 | 16292 |

| `E(t)` | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| pieces | 1 | 6 | 39 | 209 | 1013 | 3959 | 12926 | 32503 | 63121 | 86067 | 81340 | 44843 | 14336 |

(the 4555 `|S| = 0` pieces are the vertex-row contributions, one per inner graph; the `|S| <= 2` block, 39670 pieces, is what the
two `--maxout 2` controls computed.)

**The 2891 / 1036 / 627 split re-derived** from the stored pieces, outer graph by outer graph (the edge histogram of the 4555
graphs is 1, 3, 11, 34, 100, 239, 492, 784, 1002, 926, 632, 260, 71 with 0..12 edges):

| `E` | graphs | rows that receive a non-zero piece | receive nothing | coefficients at stake `(8-E)` x graphs | in rows with pieces | pieces landing there |
|---|---|---|---|---|---|---|
| 1 | 3 | 3 | 0 | 21 | 21 | 10300 |
| 2 | 11 | 11 | 0 | 66 | 66 | 24815 |
| 3 | 34 | 34 | 0 | 170 | 170 | 49315 |
| 4 | 100 | 93 | 7 | 400 | 372 | 82393 |
| 5 | 239 | 168 | 71 | 717 | 504 | 91858 |
| 6 | 492 | 204 | 288 | 984 | 408 | 60835 |
| 7 | 784 | 114 | 670 | 784 | 114 | 16292 |
| 1-7 | **1663** | **627** | **1036** | **3142** | **1655** | 335808 |
| >= 8 | **2891** | - (not computed: regular a priori, `|S|` capped at 7) | | | | |

identical in both primes and - a fact worth recording - identical in the second frame: `g5_89_nat_ab1o3_m2o7` has the same 628
stored rows, the same 627/1036/2891 split and the same 3142/1655 counts, with only the individual pieces differing.

**How much of the 1655 is really cancellation.** For each of the 1655 (row, coefficient) cells at stake inside the 627 rows I
counted the pieces that are non-zero there:

* 517 cells receive **no** non-zero contribution at that order (the selection rule again, now coefficient by coefficient);
* **1138 cells are genuine cancellations**, 1075 of them sums of more than 100 non-zero residues, 63 of 3 to 100;
* **no cell has exactly one or exactly two** contributing pieces (one would be impossible, two would be a coincidence of sign).

So of the 3142 coefficients at stake, 1487 belong to graphs that receive nothing at all, 517 more vanish for lack of a
contribution at that order, and 1138 are cancellations among 3 or more non-zero terms. Pieces per cancelling row: 31 at the
minimum, median 347, mean 536, maximum 4525.

**The largest cancelling sets** (row, `E`, pieces, widest cell):

| row | `E` | pieces in the row | largest number of non-zero residues cancelling in one coefficient |
|---|---|---|---|
| 3 | 1 | 4525 | 4525 (at `[u^1]` through `[u^7]`) |
| 14 | 2 | 4360 | 4360 |
| 48 | 3 | 3812 | 3812 |
| 1 | 1 | 3479 | 3479 |
| 11 | 2 | 3447 | 3447 |
| 47 | 3 | 3267 | 3267 |
| 10 | 2 | 3095 | 3095 |
| 46 | 3 | 2971 | 2971 |
| 125 | 4 | 2654 | 2654 |
| 119 | 4 | 2618 | 2618 |

The three one-edge rows are sums of 3479, 2296 and 4525 non-zero pieces (rows 1, 2, 3), as `NATURAL_FORM.md` s.6 records; the
vertex row is a 4555-term sum with every inner graph contributing a residue non-zero at all eight coefficients.

**The two primes agree on which pieces are non-zero.** The piece sets of `g5_89_nat` and `g5_61_nat` are identical (340363 =
340363, 0 pieces on one side only), and so is the vanishing pattern inside each piece: over all 340363 x 8 = 2722904 residue
slots there is **not one** slot that vanishes mod one prime and not the other. There is therefore no place where the two-prime
evidence is thinner than elsewhere in the sense asked - no coefficient of no piece is zero by luck at one prime only.
(Between the *frames* the pattern does differ, as it must: 36 pieces occur only in frame `(0,0)`, 30 only in frame
`(1/3,-2/7)`, and 124 shared pieces differ in their vanishing pattern; the aggregate census is nevertheless the same.)

## 4. The valuation law on the real data

Claim: a stored piece landing on an outer graph with `|S|` edges has its first `|S|` coefficients zero (each outer edge carries
`E_out = O(u)`). Checked on every stored piece of all three natural runs: **0 violations** out of 340363 + 340363 + 340369
pieces. How often the floor `|S|` is attained, i.e. coefficient `[u^|S|]` of the piece is non-zero:

| `\|S\|` | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 |
|---|---|---|---|---|---|---|---|---|
| pieces | 4555 | 10300 | 24815 | 49315 | 82393 | 91858 | 60835 | 16292 |
| floor attained | 4555 | 4525 | 9350 | 16020 | 24162 | 27177 | 24803 | 16292 |
| share | 100% | 44% | 38% | 32% | 29% | 30% | 41% | 100% |

(`|S| = 7` is attained trivially - `[u^7]` is the only coefficient left, so a piece with seven outer edges is stored only when
it is non-zero there. `|S| = 0` is not trivial: every one of the 4555 inner graphs contributes a non-zero `[u^0]` to the vertex
row.) The non-attainment is the selection rule of the true `E_out` seen on individual pieces rather than on rows: well over
half of the pieces start strictly above their floor. The second frame gives the same picture (4525, 9350, 16004, 24149, 27181,
24811, 16303), with slightly different numbers from `|S| = 3` on.

## 5. Per-inner-graph additivity, redone here for all 4555 graphs at both primes

For every inner graph `t` the sum over the outer graphs of its pieces must equal its point-frame term, `pieces_id[t][t]`. I
recomputed all 4555 identities at both primes from `parts/g5_{89,61}_nat` against `parts/g5_{89,61}_id` with my own summation:

* `g5_89_nat` vs `g5_89_id`: 4555 inner graphs checked, **0 failures**;
* `g5_61_nat` vs `g5_61_id`: 4555 inner graphs checked, **0 failures**.

This replaces three weaker statements: the collect-log line (the logs do now confirm the check ran and passed in both runs,
which the result files could not), the 2586 identities redone locally on November 13 (`local/LOCAL.md`), and the single summed
identity (natural total = point-frame total). 36440 residues on each side at each prime.

## 6. The second frame: its additivity, which the cluster never checked, now holds

The collect step of `g5_89_nat_ab1o3_m2o7` skipped additivity and said so in its log
(`(no complete g5_89_id_ab1o3_m2o7 run: the additivity check is skipped)`); `parts/g5_89_id_ab1o3_m2o7/` in the package held
only a stub log, so the run had been started and abandoned. The frame matters: all 4555 piece sums of the second-frame run
differ from the frame-`(0,0)` point-frame terms, so the existing `id` runs cannot serve as the reference.

I made the missing run here - `genus5_run.py --genus 5 --prime 89 --variant id --ab 1/3,-2/7 --run --nproc 20`, 4555 inner
graphs, wall 650 s, CPU 2.0 h, Python 3.14.6 - and then checked the identities (`cen_ab_add.py`):

* **4555 inner graphs, 0 failures**; all 4555 point-frame terms in that frame are non-zero;
* the new run's own GAP check passes in the second frame with the same `c_5`:
  `total = (1/660000) t_c^-8 + 0 t_c^-7 + ... + 0 t_c^-1` (collect log `logs/g5_89_id_ab1o3_m2o7/collect.log`, result
  `results/g5_89_id_ab1o3_m2o7.json`, part file `parts/g5_89_id_ab1o3_m2o7/local.json` - three new files in the package, nothing
  overwritten).

This closes the third of the three open points of `NATURAL_FORM.md` s.6 ("the second frame is unchecked by additivity"). It is
the same engine (`g5_engine`) but the other code path: `graph_id` computes a point-frame term directly, `graph_dp` the split
pieces whose sum is now compared with it, for all 4555 graphs in a frame where nothing had been compared before.

## 7. Genus 4 against genus 5

| | genus 4 | genus 5 | factor |
|---|---|---|---|
| leg-free stable graphs | 379 | 4555 | 12.0 |
| polar coefficients per row (`K = 2g-2`) | 6 | 8 | |
| cap on `\|S\|` used by the engine (`K-1`) | 5 | 7 | |
| regular a priori (`E >= K`), not computed | 229 | 2891 | 12.6 |
| graphs with `1 <= E <= K-1` (the content) | 149 | 1663 | 11.2 |
| stored rows (vertex + rows receiving a piece) | 69 | 628 | 9.1 |
| receive no non-zero piece (selection rule) | 81 | 1036 | 12.8 |
| vanish by cancellation | 68 | 627 | 9.2 |
| polar coefficients at stake | 267 | 3142 | 11.8 |
| of those, in rows that receive pieces | 161 | 1655 | 10.3 |
| non-zero (inner graph, outer graph) pieces | 6767 | 340363 | 50.3 |
| largest cancelling set in one coefficient | (not recorded) | 4525 | |
| CPU per full natural run | 54 min | 51-53 h | 57 |

Genus-4 figures from `NATURAL_FORM.md` s.5 (whose `229 / 81 / 68` and `6767` I re-derived here from the genus-4 edge
histogram `1, 3, 7, 21, 43, 75, 89, 81, 42, 17`); genus-5 figures from this census. The phenomenon grows faster than the graph
count: the number of pieces that have to cancel grows by 50 while the number of graphs grows by 12, because a cancelling row
at genus 5 is a sum over hundreds or thousands of inner graphs rather than tens.

## 8. What this census does and does not establish

Establishes, from the stored evidence alone:

* the seven result tables are exactly the sum of the stored pieces - rebuilt row by row, residue by residue, by code that
  shares no module with the package;
* the task decomposition is the documented one, with no duplicate, no omission, no resumed or restarted task, and the pieces
  cover `0..4554` exactly once per run; provenance (hosts, Python, SLURM ids, wall and CPU times, machine factor) recovered;
* the cancellation census `2891 / 1036 / 627` and `3142 / 1655` from the cluster's own pieces, refined to 1138 genuine
  cancellations (1075 of them with more than 100 summands) and 517 coefficients that vanish for lack of any contribution;
* the valuation law `val >= |S|` on all 1.02 million stored pieces of the three natural runs, 0 violations;
* the per-inner-graph additivity for all 4555 graphs at both primes, and now also in the second frame;
* the two primes agree on the vanishing pattern of every one of the 2722904 residue slots, so no single coefficient of the
  evidence rests on one prime only.

Does not establish:

* **independence from `g5_engine`.** The pieces are the engine's output; this is an audit of the stored numbers, not a second
  computation of them. The independent re-implementation at genus 5 still covers only 17 graphs at one prime
  (`critic/CRITIC.md`), and `genus4_referee/indep/` can now be run against the stored pieces - that remains the missing check,
  and the part files are what makes it possible.
* **the missing controls.** `tau_6` is still untested at any genus, and both genus-5 controls were `--maxout 2`, so they
  discriminate on 14 of the 1663 content graphs while 613 of the 627 cancellations sit at `E >= 3`. This census says nothing
  about that.
* **the arithmetic margin.** Two primes agreeing on which residues vanish is evidence that the vanishing is not a coincidence
  of reduction, not a proof: heights up to `2^136` against `p89 . p61 = 2^150` stand as recorded in `NATURAL_FORM.md` s.6.

## 9. Files

* `cen_lib.py`, `cen_merge.py`, `cen_prov.py`, `cen_census.py`, `cen_extra.py`, `cen_ab_add.py`, `cen_report.py` (this folder);
  run as `py -I -B <script> <paths>`, each writing a JSON summary. Total wall time about 8 minutes plus the 11-minute
  point-frame run of s.6.
* `census_g5.json`: the machine-readable summary (steps 1-6 above, including the per-row piece and per-coefficient counts of
  all 628 rows of `g5_89_nat`, the provenance block and the genus-4/genus-5 comparison).
* Evidence read: `g5_evidence.tar.gz`, extracted to the session scratchpad (`g5ev/logs/`, `g5ev/parts/`);
  `../../../cluster_genus5_quintic/results/*.json`, `data/cost_table_g5.json`, `data/graphs_{4,5}_*.json`. New files written
  into the package by s.6: `parts/g5_89_id_ab1o3_m2o7/local.json`, `results/g5_89_id_ab1o3_m2o7.json`,
  `logs/g5_89_id_ab1o3_m2o7/collect.log`.
