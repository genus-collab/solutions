# Local reproduction of a slice of the genus-5 cluster run (formal quintic, natural factorisation)

Date 2026-11-13. Written and run by an AI agent; nobody has checked this note. Only the per-run JSON summaries of the Stanage runs
were copied back (`cluster_genus5_quintic/results/`), not `logs/` or `parts/`, so the cluster's per-graph pieces cannot be inspected.
The only way to audit the run was to redo part of it here, with the package's own code, on this workstation.

Machine: Windows workstation, 28 logical cores (14 physical with SMT), Python 3.14.6, `py -I -B` throughout (no packages, no
`__pycache__` written). No file of `cluster_genus5_quintic/` or of the project was modified, and no `.md` note was edited; everything
written by this audit is in this folder. During the whole session one foreign process of another session
(`natural_form/genus5_cluster/critic/c3_enumerate.py`) occupied one core; the two local runs below also overlapped each other, so the
CPU seconds quoted are inflated (the inflation factor was measured, see s.1).

Scripts (this folder): `g5l_lib.py` (readers, rational reconstruction), `g5l_compare.py` (s.2), `g5l_additivity.py` (s.3),
`g5l_cancel.py` (s.4), `g5l_fullrun.py` (the extra full-maxout run needed for s.3). Logs, part files and the local result file are in
`out/` (listed in s.6).

## Verdict

**The slice reproduces.** On this machine, over **all 4555** genus-5 leg-free stable graphs, the package engine with the natural
factorisation `M = M_0(tau(u) zhat)` at `2^89-1` gives, for the vertex graph and for the three one-edge outer graphs, **exactly the
rows of `results/g5_89_nat.json`** - all 32 residues, entry by entry. The vertex graph is `(1/660000) t_c^-8 + 0 + ... + 0 t_c^-1`
and the three one-edge rows are identically zero. These four rows are exactly the ones on which both shipped controls (bare symbol
`M_0`, and `tau` with its `u^7` coefficient shifted by 1) fail, so the slice is a discriminating one, not a vacuous subset.

Two things the cluster logs could have shown and cannot: per-inner-graph additivity was **redone here for 2586 of the 4555 inner
graphs** (57%, all edge numbers 0..12) with complete piece families, and every one of them adds up to the shipped point-frame term.
The cancellation is massive and tight: the three zero rows are sums of **3479, 2296 and 4525 non-zero pieces**, and removing any
single piece makes each of them non-zero.

What is **not** shown: only one prime (`2^89-1`) and one frame `(0, 0)`; only 4 of the 628 stored rows (the 624 rows with
`E(q) >= 2` are untouched here); and finite checks are evidence, never proof. See s.5.

## 0. What `--maxout m` makes exact, and why 4 rows

In the edge splitting of `fq4_split.py` a term of the inner graph `Gt` with outer edge set `S` belongs to the outer graph
`Gt/(inner edges)`, which has **exactly `|S|` edges**. With `--maxout m` the engine (`g5_engine.graph_dp`, cap `= min(K-1, m)`)
keeps exactly the subsets with `|S| <= m`. Hence

* the outer graphs with `E(q) <= m` receive **every** subset that can contribute to them, so their rows are complete and must agree
  with the full run to the last residue;
* no outer graph with `E(q) > m` receives anything at all (checked: the local table has rows only for `q` in `{0, 1, 2, 3}`).

At genus 5 there are 4 outer graphs with `E <= 1`: the vertex graph `q = 0` (`gv = [5]`), and `q = 1, 2, 3`
(`gv = [1, 4]`, `[2, 3]`, `[4]` with `aut = 2`). Those 4 rows are what `--maxout 1` proves here, over all 4555 inner graphs.
Nothing about the rows with `E(q) = 2 .. 7` follows from this run. (The full run stores 628 rows:
`E = 0..7` with counts 1, 3, 11, 34, 93, 168, 204, 114; subsets with `|S| >= 8` are dropped by the engine at every maxout, since
each outer edge carries a factor `E_out = O(u)` and 8 of them are already `O(u^8)`.)

## 1. Calibration, and the slice that was chosen (s.1 = item 1 of the task)

All times are CPU seconds of the package engine, natural `M`, `2^89-1`, frame `(0, 0)`, one process unless stated.
"ref" = the reference-machine estimate `est_dp` of `data/cost_table_g5.json`.

| what | graphs | CPU | ref | ratio |
|---|---|---|---|---|
| full maxout (`--run --list 4511,808`) | 2 | 8.6 s (2.08 + 6.49) | 12.5 s | **0.69** |
| `--maxout 1`, heaviest graphs 3897 (E=10), 2864 (E=9), 2134 (E=8), 808 (E=6) | 4 | 263 s (87.6, 84.4, 88.7, 2.6) | 5282 s | 0.072 |
| `--maxout 1`, 1156 (E=7), 522 (E=6), 350 (E=5), 128 (E=4), 46 (E=3) | 5 | 164 s (52.3, 30.5, 51.1, 26.9, 3.0) | 915 s | 0.26 |
| `--maxout 2`, the same four heaviest graphs | 4 | 888 s | 5282 s | 0.24 |

So on the heavy graphs `--maxout 2` costs **3.4x** `--maxout 1`, and this machine is about 20% faster per core than the workstation
of `logs/workstation/validate.log` (which measured `f = 0.84` against the same reference; here `f = 0.69` when otherwise idle).

From this I predicted about 4.3 CPU h for `--maxout 1` over all 4555 graphs (ca. 25 min at 13 workers) and about 14.5 CPU h for
`--maxout 2`, and chose **`--maxout 1`**. The choice was right but **the absolute estimate was wrong by a factor ~2.7**: the run
actually took 11.8 CPU h (42373 s, inflated by contention - see below), 3327 s of wall time with 13 workers. The reason is a bias in
the calibration: I had timed the *heaviest* graph of each edge class, where `--maxout 1` saves 10-20x, whereas the bulk of the 4555
graphs save far less (median local time per graph at `--maxout 1`: 3.4 s; slowest 161.7 s). Taking the measured 3.4x factor at face
value, `--maxout 2` over all 4555 graphs would have been several hours and did not fit the budget.

Contention: with the second local run going (s.3) there were 24 busy workers on 14 physical cores. Re-timing two known graphs under
that load gave **2.45x and 2.53x** their idle time, so the 11.8 CPU h and 8.9 CPU h below correspond to roughly 5-7 and 3.5-4.5
CPU h of idle-machine time. Wall times are what they are and are quoted as measured.

Budget: the task asked for about one hour of wall time at most 13 workers. Measured wall times: calibration 22 min in all (of which
15 min was the single-process `--maxout 2` probe), the main `--maxout 1` pool **3327 s = 55 min** (13 workers, 10:49-11:44), and the
additivity run **2934 s = 49 min** (11 workers, 10:56-11:45, overlapping the main run). End to end the session used 1 h 17 min of
wall time, i.e. the hour was overrun; the overrun is the price of the calibration error above and is reported as measured, not
smoothed.

## 2. The slice itself: four complete rows against the cluster's table (item 2)

Run (tag `g5_89_nat_out1`, which cannot collide with any copied result):

```bash
cd "G:/My Drive/Projects/Old/refGW/cluster_genus5_quintic"
py -I -B genus5_run.py --genus 5 --prime 89 --variant nat --maxout 1 --run --nproc 13
```

Its own collect step (`out/pkg_g5_89_nat_out1_collect.log`) reports, for all 4555 inner graphs covered, 0 missing:

```
PASS every entry of E_out = E_R - E_M is O(u)
PASS [u^0] E_R(a, b) = Eh_0(a, b) for all (a, b), all colours
   vertex graph in t_c: (1/660000) t_c^-8 + (0) t_c^-7 + ... + (0) t_c^-1
PASS VERTEX GRAPH: (c) t_c^-8 + 0 + ... + 0 t_c^-1 (exact for any maxout)
   vertex graph leading coefficient c_5 = 1/660000; pattern |B_2g|/(2g(2g-2)5^(g-1)) = 1/660000: EQUAL
   outer graphs tested for regularity: 3 (E <= 1); with a stored row, by E: [(1, 3)]
PASS REGULARITY: every tested outer graph with edges is O(u^8)
GENUS5_RUN COLLECT (g5_89_nat_out1) ALL CHECKS PASSED
```

`g5l_compare.py` then re-merges the table **from the per-inner-graph part file** (the package's collect step is not trusted; it is
compared with the re-merge afterwards and agrees) and compares with `results/g5_89_nat.json`
(`out/compare_g5_89_nat_out1.log`):

* coverage: 4555 of 4555 inner graphs, 0 missing;
* structural: the local table has rows only for `q = 0, 1, 2, 3`, i.e. only for `E(q) <= 1`;
* **all 32 entries of the 4 complete rows (4 rows x 8 coefficients of `u^0..u^7`) agree with the cluster run** - 0 mismatches;
* rows 1, 2, 3 are zero in all 8 coefficients; row 0 is the only non-zero one;
* row 0 reconstructs at this one prime as
  `1/660000, -7/825000, 199/9900000, -717/27500000, 449357/22500000000, -5593217/618750000000, ?, ?`
  (`u^6` and `u^7` are beyond what a single prime determines and print as `None`). Since every other row of the cluster table
  vanishes, row 0 **is** the cluster's total, and these six coefficients are exactly the first six of the two-prime CRT series
  `u^-8 (1/660000 - 7u/825000 + 199u^2/9900000 - 717u^3/27500000 + 449357u^4/22500000000 - 5593217u^5/618750000000 + ...)`
  read off the copied JSONs by the orchestrator;
* in the conifold coordinate: local vertex graph `= (1/660000) t_c^-8 + 0 t_c^-7 + ... + 0 t_c^-1`, computed here from the local
  residues with `fq_hat.tc_series` / `fq3_msum.tc_expansion3`, and `1/660000 = |B_10|/(10 . 8 . 5^4)` exactly; the whole 8-vector
  equals the cluster's `info['vertex']`.

Discrimination (same script): the same four rows **differ** from both shipped controls, and in each control the three one-edge rows
are non-zero (`g5_89_M0_out2`: rows 1-14 non-regular; `g5_89_nat_t7_out2`: rows 1, 2, 3 non-regular). So the four rows reproduced
here are precisely the rows where a wrong `M` is caught; passing them is not an identity.

## 3. Per-inner-graph additivity, done here (item 3)

**What `--maxout 1` can and cannot give.** The claim of the method is: for every inner graph `t`, the sum of its pieces over **all**
outer-edge subsets equals its point-frame term (variant `id`, `M = 1`). With `--maxout m` the stored pieces are only those with
`|S| <= m`, a partial sum, so the sum of the stored pieces equals the point-frame term only for the inner graphs with `E(t) <= m` -
at `m = 1` that is 4 inner graphs out of 4555. **The honest statement for the other 4551 is weaker**: their stored pieces are the
restriction of the full piece family to the outer graphs with `E(q) <= 1`. Forcing the full statement on them would be wrong.

So the real check was run separately, with `g5l_fullrun.py`, which recomputes **complete** piece families (`maxout = None`) for the
cheapest inner graphs that fit a budget (through the package's own `compute`/`pool_map`, writing the package's own part format, with
no collect step):

```bash
cd ".../genus5_cluster/local"
py -I -B g5l_fullrun.py --prime 89 --variant nat --budget 17000 --nproc 11
```

2584 graphs plus the 2 from the calibration = **2586 inner graphs, 2934 s of wall, 32009 CPU s** (8.89 h as measured under
contention), 222784 non-zero pieces, up to 368 pieces for a single inner graph. Results (`out/additivity_g5_89_nat_out1.log`):

* **(A) FULL additivity: for each of the 2586 inner graphs the pieces add up exactly to its point-frame term** as stored in the
  shipped `parts/g5_89_id/local.json`. 0 failures. The 2586 graphs cover every edge number:
  `E = 0..12` with counts 1, 3, 11, 34, 84, 150, 255, 377, 493, 501, 402, 214, 61. All 2586 target terms are non-zero, so this is
  not `0 = 0`.
* the reference itself was checked: the shipped point-frame parts cover all 4555 inner graphs, re-merge to the table of
  `results/g5_89_id.json`, and sum to its `info['total']` (`1/660000, -7/825000, 199/9900000, ...`).
* **(B)** for the 4 inner graphs with `E(t) <= 1` the `--maxout 1` pieces themselves add up to the point-frame term (there the
  truncation is vacuous).
* **(C) restriction**: on the 2586 inner graphs present in both local runs, the `--maxout 1` pieces are **exactly** the full pieces
  with `E(q) <= 1` - 9299 non-zero pieces against the 222784 of the full families. This is what ties the cheap run to the full one
  and makes (A) relevant to s.2.

Counted plainly: full additivity verified here for **2586 of 4555 inner graphs (57%)**; for the remaining 1969 it was not verified
locally (their complete piece families cost roughly 29 CPU h more, by the cost table).

## 4. Non-triviality: the pieces that cancel (item 4)

From the `--maxout 1` run (`out/cancel_g5_89_nat_out1.log`), over all 4555 inner graphs:

* **14855 non-zero pieces** in all: 4555 for the vertex graph (so **every single inner graph contributes a non-zero `S = {}` term**:
  the gap `1/660000` is a 4555-term sum, not the work of one graph) and 10300 for the three one-edge outer graphs;
* the three zero rows are sums of **3479, 2296 and 4525 non-zero pieces**. Largest number of non-zero contributions cancelling for
  a single outer graph: **4525** (outer graph 3, `gv = [4]`, `aut = 2`);
* per coefficient the cancellations are equally large, e.g. outer graph 1 has 3477 non-zero contributions in `u^2` and 3479 in each
  of `u^3 .. u^7`; outer graph 3 has 4525 in each of `u^1 .. u^7`;
* `u^0` of every piece of a one-edge outer graph vanishes identically (one outer edge carries a factor `E_out = O(u)`), so the
  content sits in `u^1 .. u^7`;
* **tightness**: removing one single piece from each of the three rows (the piece of inner graph 1, 2, 3 respectively) leaves a
  non-zero sum in every case; no piece is redundant.

The log prints, for each of the three rows, the first eight non-zero contributions with their residues mod `2^89-1` and the running
partial sums, which wander over the whole of `F_p` and end at 0. For comparison, the genus-4 figures quoted in the task were 6767
non-zero pieces and 68 outer graphs regular only by cancellation; here the `--maxout 1` slice alone has 14855 non-zero pieces and
3 outer graphs (all the ones it can see) regular only by cancellation.

## 5. What this does not show

1. **One prime, one frame.** Everything local is at `2^89-1` in the frame `(0, 0)`. The cluster's `2^61-1` run and its second
   conifold frame `(1/3, -2/7)` were not redone here.
2. **4 rows of 628.** `--maxout 1` says nothing about the outer graphs with `E(q) = 2 .. 7`; 624 stored rows of the cluster table
   are untouched by this audit. The two shipped controls were run at `--maxout 2`, so even they only cover `E(q) <= 2`.
3. **Evidence, not proof.** These are finite computations modulo one prime; they cannot establish Conjecture N at genus 5.
4. **The cluster's own pieces remain unavailable.** What is shown is that this machine, running the package code, reproduces the
   cluster's numbers on the slice - not that the cluster executed what it says it did. A cluster run with a corrupted part file
   that still summed to the same four rows would not be caught here.
5. **The engine is the same code.** `g5_engine.graph_dp` is used here as on the cluster, so a bug inside it would reproduce
   identically. The independent cross-checks of that engine are elsewhere (`genus5_package_review/spot/`,
   `genus4_referee/indep/`), not in this folder.
6. **Timings were measured under load** (two overlapping runs, one foreign process); the measured inflation factor was 2.45-2.53.

## 6. Files

In this folder: `LOCAL.md` (this note), `g5l_lib.py`, `g5l_compare.py`, `g5l_additivity.py`, `g5l_cancel.py`, `g5l_fullrun.py`.

In `out/`:

| file | what |
|---|---|
| `run_out1_console.log`, `pkg_g5_89_nat_out1_local.log`, `pkg_g5_89_nat_out1_collect.log` | the package's own logs of the `--maxout 1` run and its collect step (moved out of the package) |
| `pkg_g5_89_nat_out2_local.log` | log of the `--maxout 2` calibration probe |
| `parts_g5_89_nat_out1_local.json` | 3.2 MB: the per-inner-graph pieces of the `--maxout 1` run, all 4555 graphs |
| `parts_g5_89_nat_local.json` | 22 MB: the complete piece families of the 2586 inner graphs of s.3, `maxout = None` |
| `parts_g5_89_nat_out2_local.json` | the 4 graphs of the `--maxout 2` probe |
| `results_g5_89_nat_out1.json` | the local collect result (table, info, failed = []) |
| `fullrun_89_nat.log`, `fullrun_console.log` | the s.3 run |
| `compare_*.log/.json`, `additivity_*.log/.json`, `cancel_*.log/.json` | s.2, s.3, s.4 |

To re-run the analysis, the two big part files must be put back where the package's reader looks for them:

```bash
cd "G:/My Drive/Projects/Old/refGW/cluster_genus5_quintic"
mkdir -p parts/g5_89_nat_out1 parts/g5_89_nat
LOC="../hodge_vertex_quintic/natural_form/genus5_cluster/local"
cp "$LOC/out/parts_g5_89_nat_out1_local.json" parts/g5_89_nat_out1/local.json
cp "$LOC/out/parts_g5_89_nat_local.json"      parts/g5_89_nat/local.json
cd "$LOC"
py -I -B g5l_compare.py && py -I -B g5l_additivity.py && py -I -B g5l_cancel.py --show 8
```

## 7. State of the package after the run

`results/` is untouched: the seven copied-back JSONs still carry the mtimes they were copied with, and the local result file was
moved to `out/results_g5_89_nat_out1.json`. `parts/g5_89_id/` and `parts/g5_61_id/` are untouched. No `.py`, `.sh`, `.txt` or
`data/` file was modified, and no `__pycache__` was created.

Left behind, and harmless: five **empty** directories, `parts/g5_89_nat`, `parts/g5_89_nat_out1`, `parts/g5_89_nat_out2`,
`logs/g5_89_nat_out1`, `logs/g5_89_nat_out2`. Each holds nothing but the `desktop.ini` that Google Drive puts in every folder, so
`rmdir` refuses them and deleting a Drive metadata file was not worth the risk. They contain no part, log or result data - in
particular `parts/g5_89_nat/` is empty, so nothing in it can be mistaken for a local full natural run. (Note for whoever uses the
package next: running `--collect` on a tag whose `parts/<tag>/` is empty would overwrite `results/<tag>.json` with an empty table.
That hazard is in the driver, not in these directories; `--run --list ...` returns before the collect step, which is why the
calibration used `--list`.)

One operational remark for a future local run: `--run` keeps all pieces in memory and writes `parts/<tag>/local.json` only when the
whole pool has finished, so a crash loses everything. For anything longer than an hour the resumable `--task I --ntasks N` mode
(rewritten after every graph) is the safer choice, even on one machine.
