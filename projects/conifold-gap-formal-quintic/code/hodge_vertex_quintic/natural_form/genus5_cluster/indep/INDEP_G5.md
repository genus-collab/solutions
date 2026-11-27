# Genus 5: the independent re-implementation against the cluster's own stored pieces

Date 2026-11-13. Folder `hodge_vertex_quintic/natural_form/genus5_cluster/indep/`. Written and run by an AI agent; nobody has
checked this note. Exact arithmetic modulo `p89 = 2^89-1` and `p61 = 2^61-1`, pure Python, `py -I -B`, one workstation, at most 13
worker processes, no cluster. Finite computations are EVIDENCE, not proof. Status labels: VERIFIED (exact finite computation, prime
and scope stated), ARGUED, NOT COVERED.

This closes - as far as a sample can - the item that `NATURAL_FORM.md` s.6 lists first under "what is not established at genus 5":

> No independent re-implementation covers genus 5 beyond 17 graphs at one prime; all six runs share `g5_engine`, so a structural
> zero would be common-mode. [...] Copy back `logs/` and `parts/`; then `genus4_referee/indep/` can be run against a few hundred
> stored pieces, which is the independence that is missing.

The parts are now on this machine (tarball `g5_evidence.tar.gz`), and the second implementation has been run against them.

## 0. Verdict

**No discrepancy.** On a deliberately chosen sample of **1474 of the 4555** genus-5 inner graphs (32.4%, every edge number
`E = 0..12`), the second implementation reproduces the cluster's own stored pieces **entry by entry**: **93551 pieces**
(`inner graph -> outer graph` families), **748408 individual residues**, 0 mismatches at `2^89-1`; and on **460** of them
(E = 0..5 and E = 12 complete, plus graph 2668) **19591 pieces / 156728 residues** at `2^61-1`, also all equal (VERIFIED, natural
`M = M_0(tau(u) zhat)`, frame (0, 0), `K = 8` i.e. `[u^0..u^7]`).

Omissions count as claims and were checked both ways: the cluster stores a piece only when it is non-zero, and the independent
code produced **no non-zero piece anywhere the cluster stored none** (and vice versa) - 0 one-sided pieces out of the 93551 at
`p89` and the 19591 at `p61`.

Also VERIFIED:
* **point frame**: the independent point-frame term agrees with the shipped `parts/g5_<p>_id/local.json` for **1011** inner graphs
  at `2^89-1` and **460** at `2^61-1` (0 mismatches);
* **additivity inside the independent code**: for those same 1011 + 460 inner graphs the independent natural pieces sum over the
  outer graphs to the independent point-frame term (0 failures). This is the regrouping identity tested in code that cannot share
  the package's conventions - it is also the check the cluster may never have run (`NATURAL_FORM.md` s.6);
* **two complete outer-graph rows**: for the outer graphs `q = 839` (`E = 6`, 35 non-zero contributors) and `q = 1502` (`E = 7`,
  31 contributors) every inner graph that contributes a non-zero piece was run independently, and the independent sum of those
  35 resp. 31 pieces equals the cluster's table row, namely **zero**. These are two of the 627 outer graphs that are regular only
  by cancellation, re-derived end to end by the second code;
* **graphs 149, 245 and 2668** - the three loop-heavy graphs that no independent code had ever covered (`spot/SPOT_G5.md` s.3) -
  are covered at both primes, 149 (5 loops on one genus-0 vertex, `|Aut| = 3840`, 6272 CPU-s here) by splitting its decoration sum
  into 13 shards; the sharding was validated separately on graphs 119 and 444;
* **own `|Aut|` and own canonical form for all 4555 graphs** (`ig_canon.py`): 4555 pairwise distinct canonical forms, `|Aut|`
  equal to the package's for every graph, stability / genus 5 / connectivity re-checked, 9110 contraction sanity checks;
* the momentum-space sum agrees with the independent code's own `5^V` colour-space brute force on inner graphs 1, 3, 5, 25, 28.

What is still common-mode: the hatted rows and `tau` (a second derivation of the same analytic input, not a third), the genus-5
graph list itself, the frame (0, 0), and the method (the split sum and its identification with `R'.(M.omega^)`). See s.1 and s.6.

## 1. What is independent and what is still common-mode

The second implementation is `natural_form/genus4_referee/indep/` (`ind_wk`, `ind_graphs`, `ind_tables`, `ind_sum`), imported
read-only; nothing in that folder was modified and no `__pycache__` was written. It shares no module with `fq4_split` /
`g5_engine`; the only `fq_*` import anywhere in it is inside `ind_wk.selftest`, which is not called here.

| ingredient | cluster package (`g5_engine`, `fq4_split`) | this run |
|---|---|---|
| hatted rows `Qhat`, `xhat`, `t_c` | `fq_rows`, `fq_conifold`, `fq_hat`, `fq3_msum.HatRowsU` | `natural_form/referee/ref_{fp,rows,hat}`: WKB solution of the Picard-Fuchs equation directly in `(zhat, u)`, conifold periods by a dual-number Frobenius solution. A SECOND DERIVATION, not a second copy - but also AI-written, and it was refereed against the project rows in `REFEREE_NAT.md`, so an error common to both would pass unnoticed (see s.6). |
| `tau(u)` | `fq3_m`/`fq_mfit` tables | `ind_tables`: built from the closed series `tau = 1/Sum_n 4u^n/(5n+4)` with `Fraction`s, reduced mod `p`. Independent of the package, but the SAME analytic input. |
| `E_R(a,b)` | `fq_sser.EdgeLeaf`, closed alternating sum | `ind_tables.divide`: recursion `e(m,n) = d(m+1,n) - e(m+1,n-1)`, exactness of the division by `z+w` asserted |
| leaves `T_R` | `EdgeLeaf.T` | `(-1)^e Qhat_{0,e-1}` from `T(z) = z(1 - R^-1(z)1)` |
| `M_0`, `E_M`, `E_out` | colour matrices of `fq_symbol`, `fq3_mvertex.MVertex6.Eh` | Teleman form at `u = 0` from the rows; inverse by series inversion; symplecticity, `E_M = tau^(a+b+1) E_M0` and `[u^0] E_R = E_M0` all asserted at construction |
| Witten-Kontsevich numbers | `fq_hodge.wk` (DVV on the smallest index) | `ind_wk.wk` (DVV on the largest index) |
| graph list (genera + edge lists) | `data/graphs_5_0_12.json` | **SHARED.** At genus 4 the independent code read its own Mathematica export `StableGraph-4.bin`; there is no genus-5 export, so the 4555 pairs (genus vector, edge list) and their order - which defines the outer index - are taken from the package file. Everything computed from them is recomputed here (next row). Completeness and non-redundancy of that list is `critic/CRITIC.md`'s induction on the edge count, not re-done here; non-redundancy, stability, genus 5 and connectivity are re-checked here. |
| `\|Aut\|` | `fq_graphs.autos_of` (enumeration of vertex/edge permutations and loop flips) | `ind_graphs.n_aut`: depth-first count of half-edge bijections commuting with the edge involution. Recomputed here for all 4555 graphs. |
| piece -> outer graph | union-find + `fq3_graphs.canon` | `ind_graphs.contract` (own union-find) + `ind_graphs.canon` (minimum over vertex permutations). Recomputed at run time for every subset of every inner graph in the sample. |
| the colour sum | vertex elimination in colour space, then the `graph_dp` dynamic programme over subsets | momentum (Fourier) space with the `Z/5` selection rules asserted, a separate subset loop; on small graphs also the plain `5^V` colour-space brute force |
| frame | `(alpha, beta) = (0, 0)` | the same frame `(0, 0)` (`HatRows(alpha = 0, beta = 0)`); the second conifold frame is NOT covered here |

So: the rows and `tau` are a second derivation of the same analytic input, the graph list is shared, and everything between the
rows and the stored residues - edge tables, leaves, symbol, vertices, WK numbers, automorphisms, contractions, canonical forms and
the summation algorithm - is independent.

## 2. The sample and why it was chosen

The budget was about two hours of wall time with at most 13 worker processes. The independent code is the slower of the two: a
static cost proxy (`ig_cost.py`: decorations x subsets x `5^(E-V+1)`, calibrated on the 17 timings of `spot/out/spot3_89.log`,
2.5 to 6.8 s per `10^6` proxy units) predicts **311 CPU-h for all 4555 graphs against the cluster's 51 CPU-h**, so about 6% of the
full run was affordable. The proxy was measured before committing and is the basis of the plan (`ig_plan.py`, `ig_more.py`).

The sample was therefore chosen for structure, not for cheapness:

| rule | what it gives |
|---|---|
| **complete strata at the low end**: every graph with `E <= 5` (388 graphs, including the vertex graph `t = 0`) | nothing with few edges is left out, and these are the graphs whose pieces land on the low-`E` outer graphs - the ones that are regular only by cancellation |
| **complete stratum at the top**: every graph with `E = 12` (71 graphs, the trivalent frontier, `V = 8`) | the maximal-edge end, where the subset cap `\|S\| <= K-1 = 7` bites hardest (the cap boundary is exercised: 7-outer-edge pieces occur) |
| **the package-heaviest graphs** | 6 of the 20 largest cluster wall times, including the single heaviest, `t = 3897` (3030 s on the cluster, 506 s here), and `t = 2938`, `2864`, `2863`, `2930`, `2134` (1240-1868 s on the cluster) |
| **the three graphs no independent code had reached**: 149, 245, 2668 | 149 and 245 are in the complete `E = 5` stratum; 2668 was run separately (s.3) |
| **cost-spread samples of `E = 6..11`** | each of those strata sorted by predicted cost and sampled at evenly spaced quantiles, so cheap, median and expensive graphs all occur |
| **completion of two outer-graph rows** | the 56 extra inner graphs needed to make the non-zero contributor sets of `q = 839` and `q = 1502` complete (s.4) |

Resulting coverage at `2^89-1` (natural `M`):

| `E` | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | total |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| in the sample | 1 | 3 | 11 | 34 | 100 | 239 | 283 | 213 | 185 | 146 | 115 | 73 | 71 | **1474** |
| graphs of that `E` | 1 | 3 | 11 | 34 | 100 | 239 | 492 | 784 | 1002 | 926 | 632 | 260 | 71 | 4555 |
| pieces compared | 1 | 6 | 39 | 209 | 1013 | 3959 | 8113 | 9785 | 12754 | 14068 | 15649 | 13619 | 14336 | **93551** |

Largest `|Aut|` in the sample: 3840 (graph 149, the maximum at genus 5). The pieces of the sample land on **all 628** outer graphs
that the cluster records, i.e. on the vertex graph and on **every one of the 627 outer graphs that vanish only by cancellation**.

At `2^61-1` the sample is the complete strata `E <= 5` and `E = 12` (459 graphs, plus graph 2668, so 460 graphs and 19591 pieces, 10.1% of the list); the
cost-spread strata `E = 6..11` were not repeated at the second prime for lack of time. So: **p89 broad (all `E`, 32.4%), p61 on
the complete low and top strata only** - the asymmetry the task allows for, stated explicitly.

## 3. Results

Compared object, per inner graph `t`: the whole family `{outer graph q -> [8 residues]}` of
`parts/g5_<p>_nat/task<I>of100.json`, entry by entry, with a piece stored by only one of the two codes counting as a mismatch.

| check | prime | scope | outcome |
|---|---|---|---|
| pieces, natural `M` | `2^89-1` | 1474 inner graphs, 93551 pieces, 748408 residues | **all equal**; 0 pieces present in only one code |
| pieces, natural `M` | `2^61-1` | 460 inner graphs, 19591 pieces, 156728 residues | **all equal**; 0 one-sided |
| point-frame term vs `parts/g5_89_id/local.json` | `2^89-1` | 1011 inner graphs | **all equal** |
| point-frame term vs `parts/g5_61_id/local.json` | `2^61-1` | 460 inner graphs | **all equal** |
| additivity `Sum_q` (own nat pieces) `=` (own id term) | `2^89-1` / `2^61-1` | 1011 / 460 inner graphs | **0 failures** |
| complete outer rows (s.4) | `2^89-1` | `q = 839`, `q = 1502` | **both equal the cluster's zero row** |
| momentum sum vs own `5^V` colour-space brute force | `2^89-1` | inner graphs 1, 3, 5, 25, 28 | **equal** |
| shard consistency (s.5) | `2^89-1` | graph 119 in 6 shards, graph 444 in 4 shards | merged shards **equal the cluster family** |
| own `\|Aut\|`, own canonical form, stability / genus / connectivity | - | all 4555 graphs | **all equal**, 4555 distinct forms |

All the `ind_tables` construction assertions pass at genus 5 and `K = 8`: exact divisibility of `1 - N` by `z + w` at every degree,
symplecticity of `M_0`, `[u^0] E_R = E_M0`, `E_M = tau^(a+b+1) E_M0`, `E_R(a,b)^T = E_R(b,a)`, the `Z/5` selection rule on every
edge table (`k + l = -(a+b+1) mod 5`) and on every vertex vector (exactly one non-zero Fourier component). None of these is
assumed; each is computed and checked.

The pieces are generic non-zero residues, so the agreement is not vacuous: the sample's 93551 pieces spread over all 628 stored
outer rows, up to 223 non-zero rows for a single inner graph, with up to 7 outer edges (the cap).

The three graphs that `spot/SPOT_G5.md` could not reach independently are now all covered: **149** (`E = 5`, `V = 1`,
`|Aut| = 3840`, 6272 CPU-s at `p89` in 13 shards, 352 s of wall time; also at `p61`), **245** (`E = 5`, `V = 2`, `|Aut| = 384`) and
**2668** (`E = 9`, `V = 5`, `|Aut| = 768`, run on its own after the main passes), at both primes, natural and point frame.

## 4. Two complete outer-graph rows

A piece-by-piece agreement does not by itself re-derive a regularity statement: an outer row is a sum over hundreds or thousands
of inner graphs, and the content of "this edge graph is regular" is that the sum vanishes. The cheapest rows to complete were
sought (`ig_outer.py --mode extra`): the smallest non-zero contributor set among the 627 cancelling rows has **31** members
(`q = 1502`, `E = 7`), the next **35** (`q = 839`, `E = 6`); the median row has 349 contributors and the three one-edge rows have
2296 to 4525. The 56 inner graphs needed to complete those two sets were added to the sample.

Result (`ig_outer.py --mode rows`, `2^89-1`): for both rows the sum of the independent pieces equals the row of
`results/g5_89_nat.json`, i.e. **zero in all 8 coefficients**. So two of the 627 cancellations are now independent of
`g5_engine` from end to end.

Precisely what is complete about them: every inner graph that the cluster records with a **non-zero** piece on that row was run
independently, and for each of those inner graphs the independent run also confirms the absence of every piece the cluster does
not record. That no **further** inner graph contributes to the row is the cluster's claim; it is confirmed independently only for
the 1474 inner graphs of the sample, not for the other 3081. No outer row could be completed in the strict sense within this
budget.

## 5. Timings and reproduction

Machine: the same Windows workstation as `local/LOCAL.md` (28 logical cores, 14 physical with SMT), `py -I -B` throughout, no
packages, no `__pycache__` written, nothing outside this folder modified. During the whole session a foreign 13-worker job of
another session also occupied the machine, so the CPU seconds below are inflated and the wall times are not a clean benchmark.

| pass | graphs | wall | CPU |
|---|---|---|---|
| graph 149, natural, `p89`, 13 decoration shards | 1 | 485 s | 6272 s |
| main natural sample, `p89` (13 workers) | 984 | 1857 s | ~24000 s |
| the 25 graphs completing the two outer rows, `p89` | 25 | 477 s | ~2500 s |
| supplement `E = 6..11`, `p89` (deadline stop at 463 of 466) | 463 | ~600 s | ~4600 s |
| graph 149, natural, `p61`, 13 shards | 1 | 352 s | 4576 s |
| natural sample, `p61` | 458 | 1066 s | ~13800 s |
| point frame `p89` (incl. 149 in 13 shards) | 1011 | ~510 s | 6220 s |
| point frame `p61` | 460 | ~490 s | 6032 s |
| inner graph 2668 alone, both primes, nat and id (4 single-core jobs in parallel) | 1 | 163 s | 352 s |
| own `\|Aut\|` + canonical forms, all 4555 | - | 65 s | 65 s |
| **total** | | **~2 h 15 min** | **19.0 CPU-h** |

Every finished job is appended to `out/<tag>.jsonl` and `ig_run.py` skips what is already there, so a pass can be stopped at its
`--deadline` and continued by the same command; the 30-minute rule cannot be violated by it. Sharding (`ig_lib.ShardSummer`)
strides the decoration options of vertex 0; since `graph()` is a constant times a sum over decorations, the shards add up, which
was checked against the cluster families on graph 119 (6 shards, `V = 1`, the case 149 needs) and graph 444 (4 shards).

```bash
cd "G:/My Drive/Projects/Old/refGW/hodge_vertex_quintic/natural_form/genus5_cluster/indep"
py -I -B ig_index.py                      # index of the copied-back evidence; the census (628 / 627 / 1036 / 2891)
py -I -B ig_cost.py                       # cost proxy of the independent code, calibrated on the 17 spot timings
py -I -B ig_canon.py                      # own |Aut| and own canonical form for all 4555 graphs
py -I -B ig_plan.py --budget 32500 --n-heavy 6 --p61budget 6000
py -I -B ig_outer.py --mode extra         # the inner graphs completing the rows q = 1502, 839
py -I -B ig_run.py --prime 89 --variant nat --list 149 --shards 13 --nproc 13 --tag nat_89
py -I -B ig_run.py --prime 89 --variant nat --list @sample_nat89.txt  --nproc 13 --deadline 3900 --tag nat_89
py -I -B ig_run.py --prime 89 --variant nat --list @sample_extra89.txt --nproc 13 --tag nat_89
py -I -B ig_more.py --sup 11000 --p61 18000
py -I -B ig_run.py --prime 89 --variant nat --list @sample_sup89.txt  --nproc 13 --deadline 500 --tag nat_89
py -I -B ig_run.py --prime 61 --variant nat --list 149 --shards 13 --nproc 13 --tag nat_61
py -I -B ig_run.py --prime 61 --variant nat --list @sample_nat61c.txt --nproc 13 --deadline 1700 --tag nat_61
py -I -B ig_run.py --prime 89 --variant id  --list @sample_id89_all.txt --nproc 13 --tag id_89
py -I -B ig_run.py --prime 61 --variant id  --list @sample_id61b.txt  --nproc 13 --tag id_61
py -I -B ig_run.py --prime 89 --variant nat --list 2668 --nproc 1 --tag x2668_nat89   # and nat61, id89, id61; then
                                          # cat out/x2668_*.jsonl >> the matching out/<tag>.jsonl
py -I -B ig_cmp.py --prime 89 --save cmp_89.json
py -I -B ig_cmp.py --prime 61 --save cmp_61.json
py -I -B ig_outer.py --mode rows --prime 89
```

Files: `ig_index.py`, `ig_cost.py`, `ig_canon.py`, `ig_plan.py`, `ig_more.py`, `ig_lib.py`, `ig_run.py`, `ig_cmp.py`,
`ig_outer.py`; sample lists `sample_*.txt`; results, logs and per-graph pieces in `out/` (`nat_89.jsonl`, `nat_61.jsonl`,
`id_89.jsonl`, `id_61.jsonl`, `cmp_89.json`, `cmp_61.json`, `canon_map.json`, `index.json`, `cost_ind.json`, `sample.json`,
`cmp_*.log`, `run_*.txt`, `ig_*.log`). The cluster evidence itself was read from the extracted tarball
`g5_evidence.tar.gz` (path in `ig_lib.EV` / `ig_cmp.EV`, overridable with `G5EV`); it was only read, never written.

## 6. What this does not establish

1. **It is a sample.** 1474 of 4555 inner graphs at `p89` (460 at `p61`), 93551 of the 168475 non-zero pieces (56%). The strata
   `E = 6..11` are sampled, not complete (`E = 6..11`: 283/492, 213/784, 185/1002, 146/926, 115/632, 73/260). Nothing here rules
   out a wrong piece among the 3081 inner graphs not run.
2. **No outer row is complete in the strict sense** (s.4). Two rows are complete over their non-zero contributors; the regularity
   of the other 625 cancelling rows, and the vertex row `(1/660000) t_c^-8`, still rest on the cluster's engine for the pieces
   outside the sample. In particular the three one-edge rows (2296-4525 contributors each) are far out of reach.
3. **The rows and `tau` are a second derivation, not a third.** `ind_tables` builds everything on `natural_form/referee`'s hatted
   rows (own WKB solution of the Picard-Fuchs equation, own conifold periods) and on `tau = 1/Sum_n 4u^n/(5n+4)` written out
   directly. These were refereed against the project rows in `REFEREE_NAT.md` and agree, but an error common to both
   AI-written derivations of the *same* analytic input would pass unnoticed here. The strongest independent control on the rows
   remains the gap value `c_5 = 1/660000 = |B_10|/(10 . 8 . 5^4)` and the seven vanishing subleading coefficients, and that check
   is `M`-independent (`NATURAL_FORM.md` s.6).
4. **The genus-5 graph list is shared.** The 4555 (genus vector, edge list) pairs and their order come from
   `data/graphs_5_0_12.json`, because the independent Mathematica export exists only for genus 2-4. Non-redundancy, stability,
   genus 5, connectivity, `|Aut|`, the canonical forms and the piece-to-outer assignment are all recomputed here; completeness of
   the list is `critic/CRITIC.md`'s induction on the edge count and was not re-done.
5. **One frame.** Everything here is the frame `(0, 0)`. The second conifold frame `(alpha, beta) = (1/3, -2/7)` - which has no
   additivity check of any kind, and no point-frame run - is **not** covered; `HatRows` accepts `alpha`, `beta`, so the same
   scripts would do it, at the same cost again.
6. **One method.** Both codes compute the regrouped split sum. Its identification with the outer graph sum `R'.(M.omega^)` is the
   argued composition-law identity, tested against genuine `M`-vertices only at genus 2 and 3 (`INDEP.md` s.5 item 1). A second
   implementation of the split sum cannot detect an error there, at genus 5 or at any genus.
7. **Two primes, eight orders.** All statements are residues mod `2^89-1` and mod `2^61-1` for `[u^0..u^7]`. The two primes are
   not independent as a test of the inputs (the same rationals reduced twice). No third prime was used here.
8. The controls (`g5_89_M0_out2`, `g5_89_nat_t7_out2`) and the second-frame pieces are not touched; this note compares the
   `nat` and `id` parts only.
