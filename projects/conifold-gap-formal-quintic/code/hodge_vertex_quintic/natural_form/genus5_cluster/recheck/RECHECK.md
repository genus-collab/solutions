# Independent re-check of the copied genus-5 cluster results (formal quintic, natural form)

Date 2026-11-13. Scope: the seven JSON files in `cluster_genus5_quintic/results/` and the two part files
`cluster_genus5_quintic/parts/g5_{89,61}_id/local.json` that were shipped inside the package. The cluster's own
`logs/` and `parts/` were not copied back, so the collect logs and the per-inner-graph pieces of the natural runs
are unavailable; every statement below says explicitly what it rests on.

**Code written for this re-check** (`natural_form/genus5_cluster/recheck/`, pure Python, `py -I -B`):
`rc_lser.py` (exact Laurent series over `Fraction`), `rc_pf.py` (Frobenius solution of the quintic Picard-Fuchs
equation at the conifold point, `t_c(u)`, the change of variable), `rc_rec.py` (CRT, rational reconstruction,
Bernoulli numbers), and the five drivers `rc_t0_pf.py`, `rc_t1_recon.py`, `rc_t2_tc.py`, `rc_t3_census.py`,
`rc_t4_coherence.py`, `rc_t5_sound.py`.

**Project code relied on: none.** The re-check imports no project module at all (only `fractions`, `math`, `json`,
`os`, `sys`, `collections`): no `fq_fp`, no `fq_hat`, no `fq_ring`, no `fq_series`, no `fq3_msum`, no `g5_checks`.
The primes `2^89-1`, `2^61-1`, the Bernoulli numbers, the modular inverses, the Picard-Fuchs operator and `t_c` are
all recomputed here. I did **read** `fq_conifold.py` (the operator, the exponents `{0,1,1,2}`, the normalisations
`f_k = u^k + O(u^3)`, `g = f_1 + gamma f_2` fixed by the vanishing of the `u^-1` log obstruction,
`varpi_0 = f_0 + alpha g + beta f_2`, `t_c = g/varpi_0`), `fq_hat.tc_series`, `fq3_msum.tc_expansion3`,
`g5_checks.py` and `genus5_run.py` to learn the definitions and the control flow; nothing was executed or imported
from them.

**Provenance of the files checked.** All seven results were written at 2026-11-13 10:21:54 and are unchanged:

| file | SHA-256 (first 32 hex) | bytes |
|---|---|---|
| `g5_89_nat.json` | `1550838c4bf1f35be1fdab70fb35538f` | 555080 |
| `g5_61_nat.json` | `bb02101a7410b895e075650040d91366` | 554940 |
| `g5_89_nat_ab1o3_m2o7.json` | `36010076dcf9065835aefe7ea6541716` | 555091 |
| `g5_89_M0_out2.json` | `252e4b65fc709130c07df645ef4d2b0b` | 537442 |
| `g5_89_nat_t7_out2.json` | `0bbe027cc0df5d73cc33126976a67bdc` | 534785 |
| `g5_89_id.json` | `08d78c3684c09a3d94861e4c43e04368` | 1628630 |
| `g5_61_id.json` | `96cfd7a4cbfcfb26fff2365e6b701dae` | 1325988 |

All 5048 (resp. 36456, 128) stored residues of each run lie in `[0, p)` with the right prime.

---

## (1) Rational reconstruction: CONFIRMED, unambiguous with a very large margin

`M = (2^89-1)(2^61-1)` has 46 decimal digits; `sqrt(M/2)` has 23. My reconstruction is CRT followed by the
extended Euclidean algorithm on `(M, r)` stopped at the first remainder below `sqrt(M/2)`; a result `n/d` is
certified when `2 |n| d < M`, which is the standard uniqueness bound (at most one rational of that height is
congruent to `r`).

From `info['total']` of `g5_89_nat` and `g5_61_nat`:

| | value | margin `M/(2|n|d)` |
|---|---|---|
| `u^0` | `1/660000` | 1.08e+39 |
| `u^1` | `-7/825000` | 1.24e+38 |
| `u^2` | `199/9900000` | 3.62e+35 |
| `u^3` | `-717/27500000` | 3.62e+34 |
| `u^4` | `449357/22500000000` | 7.06e+28 |
| `u^5` | `-5593217/618750000000` | 2.06e+26 |
| `u^6` | `14686042363/6496875000000000` | 7.48e+18 |
| `u^7` | `-8370858511/32484375000000000` | 2.62e+18 |

All eight coefficients agree with the orchestrator's series exactly. The worst margin is `2.6e+18`: the true
rational would have to have height more than 2.6e18 times larger than the one found for the reconstruction to be
wrong, so it is unambiguous.

The same reconstruction applied to `table['0']` (the vertex graph row) gives the **identical** eight rationals;
indeed the stored residues of `table['0']` and `info['total']` are equal at both primes, because the other 627
stored rows are identically zero.

Caveat on what "certified" means. Reconstructing from `p89` alone also returns certified rationals, but at `u^6`
and `u^7` they are *wrong* (`-4668314908965/4939268731453` and `952715118971/17416291380050`, exactly the values
the project's own id-run log prints with a leading `~`). Certification is a statement about the modulus, not
evidence of correctness; the evidence here is that two independent prime computations CRT to a rational of height
`~1e26`, i.e. 19 orders of magnitude inside the bound.

## (2) The gap, re-derived from scratch: CONFIRMED, and it is not an artefact

My own Frobenius solution of
`P = u D^4 - 2(1-u) D^3 - (7/5)(1-u) D^2 - (2/5)(1-u) D - (24/625)(1-u)`, `D = -(1-u) d/du`,
reproduces the indicial polynomial `n(n-1)^2(n-2)`, gives `gamma = 7/10` and

`t_c = u + (7/10) u^2 + (41/75) u^3 + (9/20) u^4 + (35839/93750) u^5 + (310981/937500) u^6 + ...`

(frame `(alpha,beta) = (0,0)`; the second frame `(1/3,-2/7)` gives `t_c = u + (11/30) u^2 + (751/1575) u^3 + ...`).
The change of variable is implemented as: `u^-8 G(u) = t^-8 A(t)^-8 G(u(t))` with `u(t)` the compositional inverse
of `t_c` and `A(t) = u(t)/t`, `A(0)=1`, so the eight polar coefficients are `[t^0..t^7]` of `A^-8 (G o u)`. Only
`a_0..a_7` are needed, which is exactly what the data provide.

`B_10 = 5/66` (own Akiyama-Tanigawa), so `|B_10|/(10 . 8 . 5^4) = 1/660000`.

* **total**: `(1/660000) t_c^-8 + 0 t_c^-7 + 0 t_c^-6 + 0 t_c^-5 + 0 t_c^-4 + 0 t_c^-3 + 0 t_c^-2 + 0 t_c^-1`.
  All seven subleading coefficients vanish exactly over the rationals; the leading one equals `c_5` exactly.
* **vertex graph row**: identical (it is the same series).
* **second frame** (`p89` only, so the mod-`p` version of the same computation with my own `t_c(1/3,-2/7)`):
  subleading seven residues all `0`, leading `= c_5 mod p89`. Same for its vertex row.
* My change of variable reproduces the stored `info['total_tc']` and `info['vertex']` residue for residue in all
  three natural runs at both primes. So `fq_hat.tc_series` / `fq3_msum.tc_expansion3` are independently confirmed
  at the order used.

**Closed form.** The reconstructed `u`-series is *exactly* `c_5 (u/t_c)^8` through `u^7`, coefficient by
coefficient (`(u/t_c)^8 = 1 - 28u/5 + 199u^2/15 - 2151u^3/125 + 4942927u^4/375000 - 5593217u^5/937500 +
14686042363u^6/9843750000 - 8370858511u^7/49218750000`). This is the sharpest way to state the result: the whole
eight-term polar part of the genus-5 free energy is the single number `c_5` times a mirror-map factor.

**Is the vanishing strong, or an artefact of a triangular substitution?** It is strong. The Jacobian
`dT_k/da_j = [t^k](A^-8 u(t)^j)` is lower triangular with `1` on the diagonal:

```
 a_j \ T_k   t_c^-8 | t_c^-7 | t_c^-6  | t_c^-5   | t_c^-4       | t_c^-3           | t_c^-2                | t_c^-1
 a_0              1 |   28/5 | 1063/75 |     64/3 | 7917523/375000 | 20121749/1406250 | 16433495003/2460937500 | 2137560193/1025390625
 a_1              0 |      1 |   49/10 | 1603/150 |     13587/1000 |  1554322/140625  |         5547353/937500 |      59941643/29296875
 a_2              0 |      0 |       1 |     21/5 |        769/100 |        3989/500  |          952543/187500 |          310981/156250
 a_3              0 |      0 |       0 |        1 |            7/2 |          311/60  |              1249/300  |            35839/18750
 a_4              0 |      0 |       0 |        0 |              1 |            14/5  |                  19/6  |                   9/5
 a_5              0 |      0 |       0 |        0 |              0 |              1   |                 21/10  |                 41/25
 a_6              0 |      0 |       0 |        0 |              0 |              0   |                     1  |                   7/5
 a_7              0 |      0 |       0 |        0 |              0 |              0   |                     0  |                     1
```

Consequences, all verified numerically:

* Perturbing `a_j` by `+1` changes exactly `T_k` for `k >= j` and nothing below. So `T_1 = ... = T_7 = 0` is seven
  independent affine conditions on `(a_1,...,a_7)` and says nothing about `a_0`; a generic `u`-series with the same
  leading term satisfies none of them.
* Because the diagonal is `1`, for each `a_0` there is exactly **one** series with vanishing subleading polar part.
  Solving for it with `a_0 = c_5` reproduces the reconstructed total in all eight coefficients. Equivalently: the
  claim is `G(u) = c_5 (u/t_c)^8 mod u^8`, seven exact rational identities on top of the one value `a_0 = c_5`.
* Eight non-trivial facts in total, not one: `a_0 = |B_10|/(10 . 8 . 5^4)` plus seven vanishings.

## (3) Regularity as the data actually state it: CONFIRMED, with the content localised

Each natural run stores a table of `K = 8` residues (`[u^0..u^7]` of `G`) for 628 of the 4555 outer graphs.
**In all three natural runs: every stored row with `E >= 1` is identically zero (all eight residues zero), and the
only non-zero row is index 0, the vertex graph.** The three runs store rows for exactly the same 628 indices.

Census (identical in the three natural runs):

| `E` | outer graphs | stored rows | non-zero rows | `E <= 7` with no row |
|---|---|---|---|---|
| 0 | 1 | 1 | 1 | 0 |
| 1 | 3 | 3 | 0 | 0 |
| 2 | 11 | 11 | 0 | 0 |
| 3 | 34 | 34 | 0 | 0 |
| 4 | 100 | 93 | 0 | 7 |
| 5 | 239 | 168 | 0 | 71 |
| 6 | 492 | 204 | 0 | 288 |
| 7 | 784 | 114 | 0 | 670 |
| 8 | 1002 | 0 | 0 | - |
| 9 | 926 | 0 | 0 | - |
| 10 | 632 | 0 | 0 | - |
| 11 | 260 | 0 | 0 | - |
| 12 | 71 | 0 | 0 | - |

So of the 4554 edge graphs:

* **2891 (`E >= 8`) are `O(u^8)` for free**, not by computation: the engine's cap is `|S| <= K-1 = 7` outer edges,
  an outer graph with `E` edges only receives subsets with `|S| = E`, and every entry of `E_out` is `O(u)` (a check
  that the collect step runs and that passed). Their rows are absent by construction, and that is correct, not a
  gap in coverage.
* **1663 (`1 <= E <= 7`) carry the content.** Of these, 627 received at least one non-zero piece and cancel to zero
  (3, 11, 34, 93, 168, 204, 114 with 1..7 edges), and 1036 received no non-zero piece at all (7, 71, 288, 670 with
  4..7 edges) - a selection rule of the true `E_out`, not a cancellation. How many individual non-zero
  (inner graph, outer type) pieces enter the 627 cancellations cannot be counted here: the natural runs' `parts/`
  were not copied back (at genus 4 the figure of record is 6767).

Distribution of the 1664 graphs with `E <= 7` over the number of vertices `V` (graphs / with a stored row):
`V=1` 6/6, `V=2` 55/55, `V=3` 262/232, `V=4` 590/287, `V=5` 542/48, `V=6` 186/0, `V=7` 22/0, `V=8` 1/0.
Over the number of genus-0 vertices: 0 zeros 64/61, one 219/173, two 509/250, three 826/141, four 46/3.

**Genus-4 analogue** (`g4_split_g4_89_nat.json`, re-counted here): 379 graphs, cap 5, 69 stored rows, exactly one
non-zero (the vertex graph); stored rows 1, 3, 7, 18, 24, 16 with 0..5 edges; 229 graphs with `E >= 6` free; of the
149 edge graphs with `E <= 5`, 68 cancel and 81 get no piece. These are the same numbers as `NATURAL_FORM.md` s.5.
The genus-5 pattern is a plausible continuation: full stored coverage up to `E = 2` (genus 4) resp. `E = 3`
(genus 5), then a falling fraction; the share of `E <= cap` graphs that receive any piece goes 69/150 = 46 % at
genus 4 to 628/1664 = 38 % at genus 5. Nothing anomalous.

## (4) Soundness: what the copied files do and do not establish

How the checks work (`g5_checks.checks`, `genus5_run.main`): seven named checks exist; a name is appended to
`failed` **only when it fails**, and the JSON records `failed` but not which checks ran. For a full natural run
(`--maxout` absent, variant `nat`) the first six always run:

1. all 4555 inner graphs are covered by the part files
2. every entry of `E_out = E_R - E_M` is `O(u)`
3. `[u^0] E_R(a,b) = Eh_0(a,b)` for all `(a,b)`, all colours
4. VERTEX GRAPH: `(c) t_c^-8 + 0 + ... + 0 t_c^-1`
5. REGULARITY: every tested outer graph with edges is `O(u^8)`
6. GAP: total `= (c) t_c^-8 + 0 + ... + 0`

The seventh, ADDITIVITY (for every inner graph the pieces add up to its point-frame term), runs only if
`maxout is None`, `variant != 'id'` **and** a complete same-prime, same-frame `id` run is present in `parts/`; if
not, it is skipped with a log line and leaves no trace in `failed`.

Established by the copied files alone:

* **Coverage - YES, by the `failed` list.** Check 1 always runs; `failed == []` therefore asserts that part files
  with an entry for every one of the 4555 inner graphs were present at collect time. Entries are written only
  after `compute()` returns for that graph, so this is "processed", not "listed".
  There is a second, independent handle found here: each of the 4555 point-frame terms is non-zero at **both**
  primes (checked from the shipped `parts/g5_{89,61}_id/local.json`), and the natural total equals the sum of all
  4555 point-frame terms at both primes (recomputed here). An inner graph absent from the cluster part files would
  have to have its whole contribution cancel to zero mod both primes.
* **Regularity - YES for what the data cover, with the above reading.** A row absent from the table means no inner
  graph produced a non-zero piece for it (`compute()` drops all-zero results *after* computing them), not that it
  was skipped; combined with coverage, absence is a computed zero. For `E >= 8` the zero is a priori, from
  `E_out = O(u)` (check 2 passed) and the cap, so those 2891 graphs are not evidence.
* **The gap - YES**, and re-derived here from the residues with independent arithmetic and an independent `t_c`.
* **Per-inner-graph ADDITIVITY - NO.** It cannot be established from the copied files. For
  `g5_89_nat_ab1o3_m2o7` it was certainly skipped (the required `parts/g5_89_id_ab1o3_m2o7` does not exist - the
  tag includes the frame); for the two `--maxout 2` controls it is skipped by construction. For `g5_89_nat` and
  `g5_61_nat` it would have run if the shipped `parts/g5_{89,61}_id/` were present on Stanage, which the README
  says they were - but the collect log is not available, `failed == []` is equally consistent with "ran and
  passed" and with "silently skipped", and the natural runs' `parts/` were not copied back, so it cannot be redone
  here. What *is* established is the sum of the 4555 additivity identities: natural total `=` point-frame total at
  both primes (verified here from the shipped per-graph id parts, not from the id run's own `info`).
* Also **not** established from the copied files: that the cluster's `--validate` step passed (the only
  `validate.log` present is the workstation's, which does pass: genus 2/3 against the old vertex-tensor code,
  eight genus-4 runs against the tables of record, the genus-5 set-up at both primes, the `tau` identity); that
  the part files were produced by the shipped, unmodified code; the number of non-zero pieces behind the 627
  cancellations; and per-inner-graph timings or task bookkeeping.

## (5) Cross-run coherence: no contradiction found

* The `outer` array (4555 entries, fields `E`, `gv`, `ed`, `aut`) is **identical across all seven runs**, so all of
  them index the same graph list.
* Sum over the stored rows `= info['total']` in all five runs that store a total (`g5_89_nat`, `g5_61_nat`,
  `g5_89_nat_ab1o3_m2o7`, `g5_89_id`, `g5_61_id`). The two `--maxout 2` controls store no total, exactly as the
  code path dictates.
* `g5_89_nat` total `=` `g5_89_id` total; `g5_61_nat` total `=` `g5_61_id` total. The second frame's total differs
  from frame `(0,0)` at `p89`, as it must, while its `t_c` polar part is the same.
* The shipped id part files have the right tag, 4555 entries, no missing index and **no off-diagonal entry**
  (inner index `=` outer index, as the `id` variant requires); re-merging them here reproduces
  `results/g5_{89,61}_id.json`'s table exactly, and their sum reproduces the natural runs' totals.
* The `info` key sets match the code path in every run: `nat` full runs have
  `total, total_tc, vertex, nonregular, c_total_equals_pattern, c_vertex_equals_pattern`; the `--maxout` controls
  have only `vertex, nonregular, c_vertex_equals_pattern`; the `id` runs have only
  `total, total_tc, c_total_equals_pattern`.
* **Controls are a clean `--maxout 2` restriction.** Both store exactly 15 rows, which are exactly the 15 graphs
  with `E <= 2` (1 + 3 + 11), and both fail exactly the two expected checks (VERTEX GRAPH and REGULARITY).
  - `g5_89_M0_out2` (bare symbol, `tau = 1`): all 14 edge graphs non-regular, the three one-edge graphs at
    `[u^1..u^7]` and the eleven two-edge graphs at `[u^2..u^7]`; the vertex row differs from the natural run from
    `u^1` on and its `t_c` polar part is non-zero in all eight exponents. `c_vertex_equals_pattern` is still
    `True`, which is correct and instructive: `[u^0] E_M = Eh_0` whatever `tau` is, so the leading coefficient
    alone is not the test.
  - `g5_89_nat_t7_out2` (`[u^7] tau` shifted by 1): exactly the three one-edge graphs fail, exactly at `[u^7]`;
    the eleven two-edge rows are identical to the natural run's (identically zero); the vertex row agrees with the
    natural run through `u^6` and differs at `u^7`, and its polar part acquires exactly a `t_c^-1` term. This is
    precisely what a `u^7` perturbation of `tau` must do (`E_M(a,b) = tau^(a+b+1) Eh_0(a,b)`, so the change is
    `O(u^7)`, and an `E`-edge outer graph carries `u^E` on top), and it is the genus-5 copy of the genus-4
    `nat_t5` control. The test resolves the last coefficient of `tau` that genus 5 can see.

## Residual risks and one operational warning

* The two unverifiable items above (per-inner-graph additivity; cluster-side `--validate`) would both be settled
  by copying back `logs/g5_89_nat/collect.log`, `logs/g5_61_nat/collect.log` and `logs/validate.log` from Stanage.
  Copying back `parts/g5_{89,61}_nat/` as well would allow additivity to be redone locally graph by graph.
* **Warning.** During this re-check, concurrent local work in the same session created
  `cluster_genus5_quintic/parts/g5_89_nat/`, `parts/g5_89_nat_out1/`, `parts/g5_89_nat_out2/` with a handful of
  inner graphs. A `genus5_run.py --collect --prime 89 --variant nat` run from the package directory now would
  merge those few part files and **overwrite `results/g5_89_nat.json`** with a table built from them; the coverage
  check would fail loudly, but the cluster result would be gone. The hashes above pin the version this report
  checked. Those local part directories are not part of the copied cluster results and were used in no check here.

## Verdict

Everything the task asked to be tested is confirmed: the rational reconstruction (unambiguous, worst margin
2.6e18), the gap `(1/660000) t_c^-8` with seven exactly vanishing subleading coefficients and the stronger closed
form `G = c_5 (u/t_c)^8 mod u^8`, the vanishing being seven independent conditions rather than an artefact of the
substitution, the regularity of every stored edge-graph row in all three natural runs, the census and its
consistency with genus 4, and the cross-run coherence of all seven JSONs. The one claim that the copied files do
**not** support is the per-inner-graph additivity check: it is silently skipped in at least one of the five runs
that report `failed: []`, and for the two main natural runs the copied files cannot distinguish "ran and passed"
from "skipped". Only the aggregate additivity (total = sum of all 4555 point-frame terms) is established, and that
was re-verified here independently at both primes.
