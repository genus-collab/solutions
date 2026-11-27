# Genus-5 spot checks of the cluster package (cluster_genus5_quintic), before the cluster run

Date 2026-11-11. Written and run by an AI agent; nobody has checked this note. Everything was run on a byte-identical copy of the
package in `spot/pkg/` (code, `data/`, `parts/`, `results/`). The original package folder and the project folders were not modified,
and no `__pycache__` was written (`py -I -B` throughout). At most 6 processes; the longest single run took 380 s of wall time.
Labels: BUG / RISK / MINOR / OK.

Scripts and logs (all in this folder):
`spot1.py` (step 1 and step 2), `spot3.py` (step 3), `spot4.py` (step 4); logs in `out/`: `spot1_89.log`, `spot1_61.log`,
`spot3_89.log`, `spot4.log`, and the console copies `*_console.txt`; per-graph pieces in `out/spot1_89.json`, `out/spot1_61.json`.

## Verdict

**No bug was found.** On 20 genus-5 inner graphs at 2^89-1 (15 of them also at 2^61-1), the fast engine `g5_engine.graph_dp` gives the
same result as the original slow path `fq4_split.Split.graph`, entry by entry. For each of these graphs, the pieces add up to the
shipped point-frame term. The independent genus-4 re-implementation runs at genus 5 without any code change and reproduces the pieces
of 17 of the 20 graphs.

The shipped point-frame totals agree at the two primes as rationals. The gap form `(1/660000) t_c^-8 + 0 + ... + 0` also holds with a
`t_c` computed by independent code (the referee module), at both primes.

There is one MINOR finding. The `g5_89_id` collect log prints a wrong "~" rational for the `u^6` coefficient of the total. Details
are in s.4.

## 1. graph_dp against Split.graph (natural M, K = 8, frame (0, 0))

`Split` is built once per process with `g5_engine.make_split`, so both paths use the same fast graph list and the same outer indices.
`A = graph_dp(sp, t)` and `B = sp.graph(t)` are compared as dictionaries {outer index: 8 residues}, after dropping rows that are zero
mod p.

The 20 inner graphs were chosen to cover the following (index: E, genera, feature):

| E | inner graphs | features covered |
|---|---|---|
| 1 | 1, 3 | bridge g1-g4; loop on genus 4 |
| 2 | 5 | chain g1-g1-g3 |
| 3 | 25, 28 | two loops + edge (valence 5); triple edge |
| 4 | 118, 119 | triple edge + loop; 4 loops on one genus-1 vertex (valence 8, aut 384) |
| 5 | 149, 245, 383 | 5 loops on one genus-0 vertex (valence 10, the maximum at genus 5, aut 3840); valence 9; star with 5 genus-1 leaves |
| 6 | 444 | genus-2 vertex, triple edge, loop |
| 7 | 1234 | triple edge, loop, two genus-1 vertices |
| 8 | 2318, 2539 | triple edge + loop; aut 4 (the most labelled subsets) |
| 9 | 2668 | 5 loops on genus-0 vertices of valence 6, aut 768 |
| 10 | 3638, 3847 | all genus 0, multi-edges and loops |
| 11 | 4356 | 3 loops + triple edge |
| 12 | 4484, 4552 | trivalent genus 0, V = 8 (the maximal frontier) |

Result (`out/spot1_89.log`): **`dp == split` for all 20 (OK)**. The pieces are generic non-zero residues: 2 to 223 non-zero outer rows
per inner graph, with outer-edge counts up to 7 = K - 1. The cap boundary of the subset restriction is therefore exercised: inner graph
4484 has 3 non-zero 7-edge pieces, 3847 has 19, and 4552 has 12.

The same comparison at 2^61-1 on 15 of the graphs (all except 119, 149, 245, 2668, 3847) is also **all equal** (`out/spot1_61.log`).

Timings at 2^89-1 with 6 concurrent processes:
- `graph_dp` took between 1.0x and 1.4x the `est_dp` of `data/cost_table_g5.json` on the heavy graphs (149: 59 s against 39 s;
  245: 97 s against 72 s). This is consistent with the machine factor that `--validate` measures.
- `Split.graph` was faster than `est_split`.

**maxout = 2** (inner graphs 383, 2318, 3638, 3847). `graph_dp(sp, t, 2)` was compared with my own copy of `Split.graph` restricted to
outer-edge sets of size <= 2 (`spot1.split_graph_restricted`). Results:
- all 4 are equal;
- the rows with <= 2 edges equal the rows of the full `graph_dp`, which is the claimed exactness of `--maxout`;
- no row with more than 2 edges is produced.

**OK.**

What this step does and does not test: the two paths share `vertex_table`, `plan`, `contract`, the edge tables `EM`/`EO` and the
graph list. It tests the one-elimination dynamic programme (subset bookkeeping, cap, loop and multi-edge handling, inner/outer
branching) against the per-subset colour sum. Independence of the shared parts comes from step 3.

## 2. Additivity against the shipped point-frame pieces

For every inner graph of step 1, sum_q A[q] and sum_q B[q] were compared with the shipped piece
`parts/g5_<p>_id/local.json: pieces[t][t]`. The result is **equal for all 20 at 2^89-1 and all 15 at 2^61-1 (OK)**. All of these id
pieces are non-zero, so the test is not vacuous.

## 3. Independent re-implementation at genus 5 (natural_form/genus4_referee/indep)

The independent code is genus-general apart from two things:
- its graph source, the Mathematica `StableGraph-g.bin`, does not exist for g = 5 (only genus 2-4 are in `cluster_genus4_hodge/`);
- its driver `ind_run.py` has the gap constants only up to genus 4.

`spot3.py` therefore imports its modules read-only (nothing in that folder was changed) and handles the graphs as follows:
- it takes only the genera and edge lists of the inner graphs from `data/graphs_5_0_12.json`;
- it rebuilds each graph with `ind_graphs.make`: own depth-first |Aut| and own brute-force canonical form;
- it contracts with `ind_graphs.contract`;
- it finds the package's outer index by matching its own canonical form against the package list (candidates with the same E and the
  same multiset of genera).

Tables come from `ind_tables.Tables(prime, 5, 8, 'nat')`. These are built from the referee hatted rows with their own E_R division,
WK numbers and Fourier transforms; all the Z/5 sparsity and symplecticity assertions passed at genus 5. The sum is
`ind_sum.Summer.graph` in momentum space.

Result (`out/spot3_89.log`, two runs): **17 of the 20 inner graphs are EQUAL to the package's `graph_dp` pieces, row by row, at
2^89-1:** 1, 3, 5, 25, 28, 118, 119, 383, 444, 1234, 2318, 2539, 3638, 3847, 4356, 4484, 4552. **|Aut| agrees on all 17. OK.**
- For 1, 3, 5, 25 and 28, the colour-space brute force (`graph_brute`) also equals the momentum sum.
- Not run independently: 149, 245, 2668. These are the loop-heavy graphs: 5 loops on high-valence vertices gives 5^5 momentum
  configurations times many decorations, which is too slow for the time budget. They are covered by step 1 only.
- Not run independently at 2^61-1.

This is the first genus-5 evidence that shares no module with the package for the tables, vertices, automorphisms and contractions.
It covers E = 1..12 except E = 9.

## 4. The shipped genus-5 point-frame totals at the two primes (`spot4.py`, `out/spot4.log`)

(a) **Totals [u^0..u^7]** (`results/g5_<p>_id.json: info.total`):

| u^k | rational | how it is established |
|---|---|---|
| 0 | 1/660000 | reconstructed at each prime alone; each reduces to the other prime's residue |
| 1 | -7/825000 | same |
| 2 | 199/9900000 | same |
| 3 | -717/27500000 | same |
| 4 | 449357/22500000000 | 2^89-1 alone; reduces to the 2^61-1 residue |
| 5 | -5593217/618750000000 | same |
| 6 | 14686042363/6496875000000000 | CRT only (den = 2^9 3^3 5^14 7 11) |
| 7 | -8370858511/32484375000000000 | CRT only (den = 2^9 3^3 5^15 7 11) |

For u^6 and u^7, the CRT rationals have height about 2^93 to 2^98, far below the bound of about 2^149, and denominators made of
2, 3, 5, 7 and 11 only. A random pair of residues would not do this. **The two residue lists come from the same rationals wherever
they can be reconstructed, and all 8 can. OK.**

(b) **Gap form with an independent t_c.** `t_c = g/f0` comes from the referee's `ref_rows.periods`, not from `fq_hat.tc_series`.
At both primes, [u^0..u^7] of the total equals (1/660000) [u^0..u^7] (u/t_c)^8, i.e. total = (1/660000) t_c^-8 + 0 t_c^-7 + ... +
0 t_c^-1. **OK.** 1/660000 = |B_10|/(10 * 8 * 5^4).

(c) **Per outer graph** (4555 x 8 entries of the two id tables):
- the outer lists are identical;
- no entry is zero at one prime only;
- 13129 entries have a small CRT rational (|num| den < 2^120) with a 1000-smooth denominator, so they are consistent;
- 23311 entries are undecidable with two primes (height > 2^120, which includes the u^0 term of the vertex graph);
- **no inconsistency.**

A single-prime reconstruction at the balanced bound succeeds for almost any residue, so it proves nothing by itself. One small "hit"
at 2^61-1 (inner graph 2419, u^0, 743929327/12818) was a coincidence: the CRT gives 379830627/20971520000000000000.

**MINOR: wrong "~" rational in `logs/g5_89_id/collect.log`.**
- The line "total over all graphs" prints the u^6 coefficient as `~-4668314908965/4939268731453`. This comes from the `fq3_rec.Rec.fmt`
  fallback at one prime. It does not reduce to the 2^61-1 residue, and the true value is 14686042363/6496875000000000.
- The "~" marks it as "confirm cross-prime", and the checks themselves work on residues, so nothing is falsified.
- However, the line would be wrong if quoted. The u^4 and u^5 "~" values in the same line are correct.

Suggestion: quote the totals from the table above, or let the collect step print a CRT reconstruction when both id runs are present.

## 5. Other observations

- **MINOR:** `g5_checks.costs_for` uses `est_split` as the LPT cost for the variant `id`. That is the cost of the natural
  all-subsets `Split.graph`, not of the single-subset id term, so the load balance of an id array run would be poor. This only matters
  if the id runs are ever redone on the cluster; they were done locally.
- **OK:** none of the 40 ingredient files differed between `spot/pkg` and the original package at the end. The `fq*.py` and
  `gsum_q.py` files are byte-identical to `hodge_vertex_quintic/`.

## 6. What is not covered here

- the full genus-5 natural run (that is the cluster job);
- the second conifold frame;
- the controls M0 and nat_t7;
- independent verification of inner graphs 149, 245 and 2668;
- independent verification at 2^61-1;
- per-graph rationals beyond 2^120 height (two primes are not enough for these).
