# The two genus-5 controls: do they discriminate, and do they fail in the predicted way?

Formal quintic, natural factorisation `M = M_0(tau(u) zhat)`, `1/tau = 2F1(1, 4/5; 9/5; u)`; genus 5 near the conifold point
`u = 1 - 5^5 q = 0`. Source data: the six JSONs copied back from Stanage in
`cluster_genus5_quintic/results/` (`logs/` and `parts/` were not copied back), the stored genus-4 tables of record
`cluster_genus5_quintic/data/ref_g4.json`, and small local runs at genus 2 and 3 made for this note.
Code of this note: `controls/ctl_read.py`, `ctl_val.py`, `ctl_gsum.py`, `ctl_lin.py`, `ctl_floor.py`, `ctl_consist.py`,
`ctl_cross.py`, `ctl_g4.py`, `ctl_cost.py` (new files; nothing in the package or in any existing module was modified).

## 0. Verdict

Both controls fail, and they fail in exactly the pattern that edge splitting predicts, with no unexplained feature:

| control | predicted non-regular outer graphs | predicted orders | observed |
|---|---|---|---|
| `g5_89_M0_out2` (`tau = 1`) | every tested graph: all 3 with `E = 1`, all 11 with `E = 2` | from `u^E` to `u^7` | 14 graphs; `E = 1` at `u^1..u^7`, `E = 2` at `u^2..u^7` |
| `g5_89_nat_t7_out2` (`tau_7 + 1`) | only the 3 graphs with `E = 1` (two-edge response sits at `u^8`, outside the window) | `u^7` only | 3 graphs, each at `u^7` only |

and the vertex graph is spoilt at `t_c^-7..t_c^-1` in the first and at `t_c^-1` only in the second, again as predicted.

Three findings go beyond "the pattern matches":

* **(A) A sum rule.** The row sum over all outer graphs is the point-frame total, which does not contain `M` at all, so the
  perturbation of the rows must cancel exactly. Truncated by valuation this gives three identities computable from the
  copied-back data alone; all three hold mod `2^89-1` (`ctl_consist.py`).
* **(B) The two controls measure the same four numbers.** The floor-order response of the vertex row and of the three
  one-edge rows to `tau_k -> tau_k + c` is `c . A_Gamma` with `A_Gamma` **independent of `k`**. Since the bare symbol is
  `tau' = 1 = tau + (1 - tau)` with lowest term `c_1 = -tau_1 = 4/9`, this predicts
  `[u^1](row(M0) - row(nat)) = (4/9) . [u^7](row(nat_t7) - row(nat))` for those four rows. All four identities hold exactly
  (`ctl_cross.py`), and the same identity holds on the full genus-4 tables of record (`[u^1]` of `M0` against `[u^5]` of `nat_t5`).
* **(C) Hence the answer to the sharpness question is better than expected.** Because `A_Gamma` does not depend on `k` and
  the four measured `A_Gamma` are non-zero, the test is sensitive to `tau` at **every** order `u^1..u^7`, not only at `u^7`:
  any `tau' != tau` with `tau'(0) = 1` whose lowest deviation is at order `k0 <= 7` makes the three one-edge graphs fail at
  `[u^k0]`. Within the ansatz `M = M_0(tau(u) zhat)` the natural `tau` is the *unique* one mod `u^8` for which the one-edge
  graphs are regular. No further control run is needed for the `tau` direction.

The one thing to be suspicious of (`c_vertex_equals_pattern = True` in both controls) is not a problem: the leading polar
coefficient is fixed by the `u^0` data, which is identical in every variant with `tau(0) = 1`, so that flag cannot
discriminate and must not be counted as a partial success. See section 7.

## 1. What the JSONs actually say

`ctl_read.py` (reads `results/*.json` only). Both controls used `--maxout 2`, which is exact for the vertex graph and for the
outer graphs with at most two edges; their tables hold 15 rows (indices 0..14 = exactly those graphs). `info['nonregular']`
is indexed by **`u`-order** (the raw row `[u^0..u^7]`; regular means all eight vanish), whereas `info['vertex']` is the
**`t_c`-expansion** of the vertex row (coefficients of `t_c^-8 .. t_c^-1`). Both coordinates appear below; they are related by
the triangular substitution `t_c(u) = u + O(u^2)` (valuation 1, leading coefficient 1, `ctl_val.py`), so `[u^j]` of a row
feeds `[t_c^(j-8)]` and nothing below it.

### `g5_89_M0_out2` — bare symbol, `tau = 1`

`failed`: the VERTEX GRAPH check and the REGULARITY check. The coverage check, `frozen_ok`
(`[u^0] E_R = Eh_0`) and `E_out = O(u)` all passed, so all 4555 inner graphs were summed.
14 non-regular outer graphs = **every** graph tested:

| idx | E | `gv` | `ed` | \|Aut\| | failing `u`-orders |
|---|---|---|---|---|---|
| 1 | 1 | [1, 4] | [[0,1]] | 1 | 1,2,3,4,5,6,7 |
| 2 | 1 | [2, 3] | [[0,1]] | 1 | 1,2,3,4,5,6,7 |
| 3 | 1 | [4] | [[0,0]] | 2 | 1,2,3,4,5,6,7 |
| 4 | 2 | [0, 4] | [[0,0],[0,1]] | 2 | 2,3,4,5,6,7 |
| 5 | 2 | [1, 1, 3] | [[0,1],[1,2]] | 1 | 2,3,4,5,6,7 |
| 6 | 2 | [1, 1, 3] | [[0,2],[1,2]] | 2 | 2,3,4,5,6,7 |
| 7 | 2 | [1, 2, 2] | [[0,1],[0,2]] | 2 | 2,3,4,5,6,7 |
| 8 | 2 | [1, 2, 2] | [[0,1],[1,2]] | 1 | 2,3,4,5,6,7 |
| 9 | 2 | [1, 3] | [[0,0],[0,1]] | 2 | 2,3,4,5,6,7 |
| 10 | 2 | [1, 3] | [[0,1],[0,1]] | 2 | 2,3,4,5,6,7 |
| 11 | 2 | [1, 3] | [[0,1],[1,1]] | 2 | 2,3,4,5,6,7 |
| 12 | 2 | [2, 2] | [[0,1],[0,1]] | 4 | 2,3,4,5,6,7 |
| 13 | 2 | [2, 2] | [[0,1],[1,1]] | 2 | 2,3,4,5,6,7 |
| 14 | 2 | [3] | [[0,0],[0,0]] | 8 | 2,3,4,5,6,7 |

Vertex graph: `[t_c^-8]` correct, `t_c^-7 .. t_c^-1` all non-zero. In `u`-coordinates the vertex row agrees with the natural
run at `u^0` and differs at `u^1..u^7`.

### `g5_89_nat_t7_out2` — natural `tau` with `tau_7 -> tau_7 + 1`

`failed`: the same two checks. 3 non-regular outer graphs, the three with one edge:

| idx | E | `gv` | `ed` | \|Aut\| | failing `u`-orders |
|---|---|---|---|---|---|
| 1 | 1 | [1, 4] | [[0,1]] | 1 | 7 |
| 2 | 1 | [2, 3] | [[0,1]] | 1 | 7 |
| 3 | 1 | [4] | [[0,0]] | 2 | 7 |

The eleven two-edge graphs (indices 4..14) stay **exactly** regular. Vertex graph: `t_c^-8 .. t_c^-2` clean, `t_c^-1`
non-zero; in `u`-coordinates the vertex row differs from the natural run at `u^7` and nowhere else.

For contrast, in `g5_89_nat` all 628 stored rows vanish except the vertex row, so every failure above is a clean
`0 -> non-zero`, not a change of a non-zero number.

## 2. The prediction, and where it comes from

Write `E_R = E_M + E_out`, `E_M(a, b) = tau^(a+b+1) Eh_0(a, b)`, `E_out = E_R - E_M = M^-1 E_R' M^-T`. The row of an outer
graph `Gamma` is a sum over inner graphs `Gt` and over edge subsets `S` of `Gt` with `Gt/S^c = Gamma`, of a product over the
edges of `Gt`: a factor `E_out` on each of the `|E(Gamma)|` edges of `S` and a factor `E_M` on each contracted edge, times
leaves and trivial-theory vertices. Two facts drive everything.

**(i) `val_u E_out = 1` exactly.** `[u^0] E_R(a, b) = Eh_0(a, b)` (the check `frozen_ok`) and `tau(0) = 1`, so
`[u^0] E_out = 0`; measured with `ctl_val.py`, every one of the 1950 entries of `E_out` at genus 5 has valuation exactly 1,
for the natural `tau`, for `tau = 1` and for every `nat_tk`. Hence

> the row of an outer graph with `E` edges has `val_u >= E`,

which is why graphs with `E >= 2g-2 = 8` are regular for free (the engine's subset cap), and why in `g5_89_nat` no row with
`E >= 8` is stored at all. The stored rows run `E = 0..7` (1, 3, 11, 34, 93, 168, 204, 114 of them).

**(ii) A shift of `tau` at order `k` perturbs one edge factor at order exactly `k`.** For `tau' = tau + c u^k`,

```
E_M'(a,b) - E_M(a,b) = [(tau + c u^k)^m - tau^m] Eh_0(a,b)       (m = a+b+1)
                     = m tau^(m-1) c u^k Eh_0(a,b) + O(u^(2k)),
E_out' - E_out       = -(E_M' - E_M),
```

so both edge tables move at valuation exactly `k` (measured: all 1950 entries, every `k = 1..7`, `ctl_val.py`). Because the
row is multilinear in its edge factors, perturbing one edge and leaving the others alone gives

> `val_u (row'(Gamma) - row(Gamma)) >= k + (E - 1)` for `E >= 1`, and `>= k` for the vertex graph,

the `E - 1` because the remaining outer edges still carry a factor of valuation `>= 1`. Perturbing two edges costs `2k`.
Visible inside the window `u^0..u^(2g-3)` iff `k + E - 1 <= 2g-3`. At genus 5:

* `tau_7` shifted: `7 + E - 1 <= 7` forces `E <= 1`. Only the vertex graph (at `u^7`) and the three one-edge graphs (at
  `u^7`, the single order left in the window) can move. **Exactly what the run shows.**
* bare symbol: `tau' = 1` differs from `tau` first at order 1, so `k = 1` and `E <= 7` can fail; with `--maxout 2` the
  graphs tested are `E = 1, 2` and they fail from `u^1` and `u^2` respectively. **Exactly what the run shows**, including
  the fact that the two-edge graphs are clean at `u^1` (a two-edge row cannot have a `u^1` term at all, by (i)).

The bound is a lower bound on the first failing order; that the floor is **attained** (no accidental cancellation) is
an empirical fact, confirmed here in every one of the 173 rows that move at genus 4 and 5 (150 rows of the genus-4 `M0`
table, 4 of `nat_t5`, 15 of the genus-5 `M0` control, 4 of `nat_t7`) and violated in two known places (section 3).

### Genus 2, 3, 4 (same code, full sums)

`ctl_gsum.py` runs the complete outer graph sum for every `tau` control at a small genus (seconds). Writing `k = 1` for the
bare symbol, the prediction "fail set `= {E : k + E - 1 <= 2g-3}`, orders `u^(k+E-1)..u^(2g-3)`" gives:

| genus | control | predicted fail set | observed |
|---|---|---|---|
| 2 | `M0`, `nat_t1` | the 2 graphs with `E = 1` | 2 graphs at `u^1` — MATCHES |
| 3 | `M0`, `nat_t1` | `E <= 3`: 2 + 5 + 9 = 16 graphs | 16 graphs, first order `= E` — MATCHES |
| 3 | `nat_t2` | `E <= 2`: 7 graphs | 6 graphs (first order `k+E-1`); the chain `gv=[1,1,1]`, `ed=[[0,2],[1,2]]` stays regular |
| 3 | `nat_t3` | the 2 graphs with `E = 1`, at `u^3` | 2 graphs at `u^3` — MATCHES |
| 4 | `M0` (table of record) | `E <= 5`: 149 graphs | 149 graphs, failing orders exactly `(E, E+1, ..., 5)` — MATCHES |
| 4 | `nat_t5` (table of record) | the 3 graphs with `E = 1`, at `u^5` | 3 graphs at `u^5` — MATCHES |

`ctl_g4.py` turns the genus-4 record into a per-graph valuation statement on all 379 rows: for `M0` the perturbation has
valuation exactly `max(E, 1)` for every graph with `E <= 5` and vanishes mod `u^6` for `E >= 6`; for `nat_t5` it has
valuation exactly 5 on the vertex graph and the three one-edge graphs and vanishes mod `u^6` on all 375 others. No graph
violates the floor. The genus-4 record in `NATURAL_FORM.md` s.5 ("with `M = M_0` all 149 outer graphs with `E <= 5` are
singular; with `tau_5` shifted exactly the three one-edge graphs fail, at `[u^5]`, and the vertex graph acquires a `t_c^-1`
term") is reproduced here from the stored tables, and the genus-5 numbers are its exact analogue one genus up.

## 3. Where the valuation floor is not attained

The prediction is an inequality. Two places are known where the floor coefficient happens to vanish, so a graph that "could"
fail does not:

* genus 3, `nat_t2`, the two-edge chain `gv = [1, 1, 1]`, `ed = [[0,2],[1,2]]`: its only visible order `u^3` vanishes. The
  response is exactly linear in the shift `c` at that order (`ctl_lin.py`, `c = 1, 2, 3`), so it vanishes for every `c`:
  this row is exactly blind to `tau_2`, a selection rule of the true `E_out`, not a numerical accident.
* genus 4, mutant `m04_tau4` of `natural_form/genus4_referee/review/REVIEW.md` s.4.2 (`1/tau` wrong at `u^4`): five of the
  seven two-edge graphs fail at `u^5`, not seven.

Both are in the "other outer edges" part of the bound, i.e. for `E >= 2`. For `E <= 1`, where the floor response is the
`k`-independent number `A_Gamma` of section 5, no vanishing is known at any genus, and `A_Gamma != 0` is verified at genus 4
and 5. The lesson for reading a control: a graph that fails proves sensitivity; a graph that does not fail proves nothing.

Linearity in the shift (`ctl_lin.py`, genus 3): for `k >= 2` every visible perturbation is exactly linear in `c`, because the
second-order term sits at `u^(2k)`, outside the window. For `k = 1` the one-edge rows are still linear at their floor order
`u^1` but the two-edge rows are not at `u^2` (both edges can be perturbed: a `c^2` term). This is the only non-linearity in
the whole picture and it is predicted.

## 4. The numbers 3, 11, 14

Leg-free stable graphs of genus 5 by number of edges, from the shipped list (`data/graphs_5_0_12.json` through the `outer`
array of any result JSON; `ctl_read.py`):

```
E   = 0   1   2   3    4    5    6    7     8     9    10   11  12
#   = 1   3  11  34  100  239  492  784  1002   926   632  260  71      (4555 in all)
```

* **`E = 1`: 3 graphs.** A one-edge genus-`g` graph is either a loop at a vertex of genus `g-1` (1 graph) or a bridge
  splitting `g = g_1 + g_2` with `g_1, g_2 >= 1` (`floor(g/2)` graphs; a genus-0 end with one half-edge is unstable). At
  `g = 5`: the loop on `gv = [4]` and the bridges `gv = [1,4]`, `[2,3]` — the three graphs that fail in both controls. (The
  count is 3 at genus 4 too, `[3]` loop, `[1,3]`, `[2,2]`, which is why the genus-4 and genus-5 `nat_tk` controls both
  report "exactly three".)
* **`E = 2`: 11 graphs.** By topological type: two loops at one vertex, `gv = [3]` (1); loop + bridge, loop at the genus-`a`
  end with `a + b = 4`, `a >= 0`, `b >= 1`, giving `[0,4]`, `[1,3]`, `[3,1]`, `[2,2]` (4); double edge between two vertices of
  genera `(1,3)`, `(2,2)` — genus 0 with two half-edges is unstable (2); path of two bridges with all three genera `>= 1`
  summing to 5 (4). Total `1 + 4 + 2 + 4 = 11`, and the list of eleven in section 1 is exactly these.
* **14 = 3 + 11** is therefore precisely the number of outer graphs with `1 <= E <= 2`, i.e. *every* graph that `--maxout 2`
  computes exactly apart from the vertex graph. There is no discrepancy to explain: the bare symbol fails on everything the
  control can see, exactly as it failed on all 149 testable graphs at genus 4. The control's table has 15 rows = those 14
  plus the vertex graph, which is not an "outer graph with edges" and is judged by the separate VERTEX GRAPH check.

## 5. Two quantitative tests the copied-back data makes possible

### 5.1 The perturbation sums to zero (`ctl_consist.py`)

The sum of the rows over **all** outer graphs is the point-frame total, which contains no `M` (the copied-back `id` and `nat`
runs agree at both primes). So `sum_Gamma (row'(Gamma) - row(Gamma)) = 0` mod `u^8` identically. Truncating by the
valuations of section 2 leaves only the graphs the controls computed exactly:

| control | order | graphs that can contribute | sum of the 4 resp. 15 perturbations mod `2^89-1` |
|---|---|---|---|
| `nat_t7` | `u^7` | `E <= 1`: rows 0, 1, 2, 3 | 0 — PASS |
| `M0` | `u^1` | `E <= 1`: rows 0, 1, 2, 3 | 0 — PASS |
| `M0` | `u^2` | `E <= 2`: rows 0..14 | 0 — PASS |

Every one of the 23 individual perturbations entering these sums is non-zero, so these are not trivial identities; each is a
1-in-`2^89-1` coincidence if the runs were inconsistent. They tie the two control runs to the natural run arithmetically,
which is worth having since `parts/` and the collect logs were not copied back.

### 5.2 The floor-order response does not depend on `k`, and the two controls measure it twice (`ctl_floor.py`, `ctl_cross.py`)

At the floor order of a row with `E <= 1` the perturbed factor consumes the whole order `u^k`, so every other factor —
the inner `E_M`s, the leaves, the vertices — must be taken at `u^0`, and the perturbed factor contributes
`m c Eh_0(a,b)` with `[u^0] tau^(m-1) = 1`. Both are free of `k`. Hence

```
[u^k] (row(tau + c u^k) - row(tau)) = c . A_Gamma ,   A_Gamma independent of k,   for E = 0 and E = 1.
```

Verified directly at genus 3 for `k = 1, 2, 3` and `c = 1, 2`: the vertex row and the two one-edge rows give the same number
for all six (`k`, `c`) combinations, while the two-edge rows do not (their floor order forces one factor of `E_out` at `u^1`,
and for `k = 1` a second perturbation fits) — the predicted exception.

Now the bare symbol is itself a shift of `tau`: `tau' = 1 = tau + (1 - tau)`, lowest term `c_1 = -tau_1 = 4/9`. So the two
controls must satisfy, on the vertex row and the three one-edge rows,

```
[u^1] (row(M0) - row(nat))  =  (4/9) . [u^7] (row(nat_t7) - row(nat)).
```

All four identities hold exactly mod `2^89-1` (`ctl_cross.py`), with

```
A_0 = 220289747891123609804533358    (vertex graph, gv = [5])
A_1 = 515952397907607357026433825    (gv = [1,4])
A_2 = 594166047723601843970542697    (gv = [2,3])
A_3 = 526501865405737601547176453    (gv = [4], loop)
```

(residues only: single-prime reconstructions are not evidence). The same identity holds on the full genus-4 tables of
record, `[u^1]` of `g4_89_M0` against `(4/9) [u^5]` of `g4_89_nat_t5`, for the vertex graph and the three one-edge graphs.
Two runs made independently on the cluster, in different variants, agree on four numbers through a prediction derived from
the composition law — and those numbers are non-zero.

## 6. Sharpness: at which orders is the test really sensitive to `tau`?

**Read literally, the two runs resolve only two orders.** The `nat_t7` run isolates `tau_7`. The `M0` run changes all of
`tau_1..tau_7` at once; its failure at `u^1` on the one-edge rows does isolate `tau_1`, because `[u^1]` of a one-edge row
depends on `tau` through `tau_1` only (a change at order `k` cannot reach below `u^k`). Its failures at `u^2..u^7` are joint
functions of `tau_1..tau_j` and, with `tau_1` already wrong, carry no information about `tau_2..tau_6` individually. So from
the two JSONs alone: sensitivity proved at `u^1` and `u^7`, five orders unresolved.

**Read with section 5.2, all seven orders are covered, with no further run.** `A_Gamma` does not depend on `k`, and the
measured `A_1, A_2, A_3` are non-zero. Therefore for every `k = 1..7` the control `nat_tk` would make exactly the three
one-edge graphs fail, at `[u^k]`, with those same coefficients. More strongly, for any `tau2` with `tau2(0) = 1` whose first
deviation from the natural `tau` is at order `k0 <= 7`, with coefficient `c`, the one-edge rows have
`[u^k0] = c A_Gamma != 0`: the test detects **every** deviation of `tau` visible in the genus-5 window, and the natural `tau`
is the unique `tau` mod `u^8` making the one-edge graphs regular. (Both controls are instances: `M0` with `k0 = 1`,
`c = 4/9`, and `nat_t7` with `k0 = 7`, `c = 1`.) This is a consequence of the composition law plus four measured residues,
not a finite check; the measurement is at one prime only.

**If a direct run per order is wanted anyway** (the belt-and-braces version, which also re-tests the engine at each order):

* `--variant nat_tk --maxout 1` for `k = 2..6` is enough, since the one-edge graphs are the discriminating ones and
  `maxout 1` is exact for them. Timing three calibration graphs at `maxout 1, 2, None` (`ctl_cost.py`) gives
  `maxout1/full = 0.24`, `maxout2/full = 0.48`, so with the package's own cost table (46 CPU h per full natural run on the
  reference workstation) **one `maxout 1` control is about 11 CPU h, one `maxout 2` control about 22 CPU h**; five of them
  about 55 CPU h, i.e. a little over one full natural run, and with a 100-task array about 0.1-0.2 h per task — minutes of
  wall time on Stanage. (`--maxout 2` would also show the two-edge response for `k <= 6`, at twice the cost.)
* Much cheaper, and partly already done: the same controls at **genus 4** cost about 54 CPU min each (full, 379 graphs) and
  resolve `k = 1..5`. `nat_t5` is of record; `m02_tau1` and `m03_tau3` of the referee's mutant table cover `k = 1, 3` at
  genus 2/3, `m04_tau4` covers `k = 4` at genus 4, and the genus-3 `nat_t2` run of this note covers `k = 2`. **The only
  coefficient of `tau` never probed by a direct run at any genus is `tau_6`**, which needs genus 5 (the genus-4 window stops
  at `u^5`): one `--maxout 1` run, about 11 CPU h.
* A single run could in principle resolve all seven orders at once by carrying the shifts as seven square-zero parameters
  (`tau -> tau + sum_k eps_k u^k`, `eps_k eps_l = 0`), i.e. computing the linear responses `A_Gamma` directly. That is the
  same statement as section 5.2 and not worth coding.

**What the controls do not probe at all.** Only the `tau` direction is tested. The `zhat`-orders of the symbol `M_0` itself
are tested at genus 4 by the referee's mutant `m07_taupow_hi` (wrong `tau` power for `a+b+1 >= 7`: invisible at genus 3, 21
non-regular graphs at genus 4). Genus 5 uses symbol data that genus 4 never touches (`KZ = 12`, edge keys with `a + b` up to
11, leaves up to `psi^13`), and no control perturbs those orders. An `m07`-style control at genus 5 — a wrong `tau` power or
a shifted colour matrix at `zhat >= 10`, `--maxout 1` or 2 — is the missing control, not another `tau_k`.

## 7. The suspicious flag: `c_vertex_equals_pattern = True` in both controls

This is forced, and it is the right behaviour, not a sign that the controls half-worked.

1. `[u^0]` of every row is variant-independent whenever `tau(0) = 1`: then `[u^0] E_M = Eh_0 = [u^0] E_R` (the separately
   reported `frozen_ok`), and the leaves, vertices and `xhat` powers do not involve `M` at all. All three `2^89-1` runs do
   report the identical residue `460484621779905517725715216` for `[u^0]` of the vertex row, bit for bit, and the `M0` row
   first differs from the natural row at `u^1`, the `nat_t7` row at `u^7`.
2. `t_c(u)` has valuation 1 and leading coefficient 1, so `[t_c^-8]` of a row equals `[u^0]` of that row exactly; the whole
   polar part is triangular in the `u`-coefficients, which is why a defect at `u^j` shows up at `t_c^(j-8)` and in no higher
   pole: `M0` spoils `t_c^-7..t_c^-1`, `nat_t7` spoils `t_c^-1` alone.
3. Every graph with `E >= 1` has `val_u >= 1` in **every** variant (section 2 (i), checked in the controls' own collect step
   as `E_out = O(u)`), so `[u^0]` of the total equals `[u^0]` of the vertex row — and `[u^0]` of the total is
   `M`-independent, being the point-frame value (verified: the `id` and `nat` totals agree at both primes, and even the
   second conifold frame gives the same `[u^0]`).

So `c_vertex = [u^0](point-frame total) = c_5 = 1/660000` for *any* factorisation with `tau(0) = 1`, including deliberately
wrong ones. `c_vertex_equals_pattern` tests the gap value (a real result of the natural run and of the `id` runs) but has no
power to test `M`; it is not a discriminating check and should not be listed among the things a control confirms. The two
checks that do test `M` — "the vertex graph's polar part is `c t_c^-8` and nothing else" and "every tested outer graph with
edges is `O(u^8)`" — are precisely the two entries of `failed` in both controls. (Consistently, neither control reports
`total` or `c_total_equals_pattern`: with `--maxout 2` there is no total to report.)

A related non-check: the per-inner-graph ADDITIVITY check and the GAP check are also blind to `E_M` (they only use
`E_R = E_M + E_out`), as the genus-4 referee already noted; `m04_tau4` and `m07_taupow_hi` pass GAP at genus 4 while failing
regularity. The controls' failing exactly the two `E_M`-sensitive checks is the expected signature.

## 8. Caveats

* Both controls are single-prime (`2^89-1`) and `--maxout 2`. The sum rules and the cross-identities of section 5 are
  therefore verified mod `2^89-1` only. Rerunning either control at `2^61-1` costs the same as the original (about 22 CPU h)
  and would square the arithmetic confidence; it is not needed for the logic.
* `logs/` and `parts/` were not copied back, so coverage and per-inner-graph additivity could not be re-derived here. They
  are implied by the controls' own `failed` lists, which contain the two designed failures and **not** the coverage check,
  `frozen_ok`, or `E_out = O(u)` — all of which therefore passed inside the collect step on the cluster, over all 4555 inner
  graphs.
* `--maxout 2` means the rows with `E >= 3` in the control JSONs are incomplete and must not be read; only rows 0..14 are
  exact. The `nonregular` lists respect this (the engine only tests `E <= maxout`).
* Finite checks are evidence, never proof. What section 5.2 adds is an argument (the floor-order response is `k`-free) whose
  only numerical input is four non-zero residues; the argument itself uses only the composition law `E_R = E_M + E_out` with
  `E_M = tau^(a+b+1) Eh_0`, which is the definition in use.

## 9. Files

| file | what it does | runtime |
|---|---|---|
| `ctl_read.py` | tabulates the two controls from `results/*.json`; edge histograms; recomputes the genus-4 control record from `data/ref_g4.json` | instant |
| `ctl_val.py` | `val(E_out)` per entry for every variant at genus 2-5; `val(E_out(control) - E_out(nat)) = k`; `t_c` valuation | genus 5 set-up about 4 s |
| `ctl_gsum.py` | full outer graph sums at genus 2/3 for `nat`, `M0`, `nat_t1..nat_t(2g-3)`, with the fail-set prediction | 2-3 s per variant |
| `ctl_lin.py` | linearity of the response in the shift `c` (genus 3, `k = 1, 2, 3`, `c = 1, 2, 3`) | about 30 s |
| `ctl_floor.py` | the floor-order response is independent of `k` (genus 3) | about 20 s |
| `ctl_consist.py` | the three sum rules of 5.1, from the copied-back JSONs | instant |
| `ctl_cross.py` | the four cross-identities of 5.2 at genus 5 and their genus-4 analogue | instant |
| `ctl_g4.py` | per-graph valuation of the perturbation on all 379 genus-4 rows, both controls | instant |
| `ctl_cost.py` | times `graph_dp` at `--maxout 1, 2, None` on the package's three calibration graphs and scales the cost table | about 1 min |

All are run as `py -I -B <file>` from this directory; they read the package `cluster_genus5_quintic` and write nothing.
`ctl_lin.py` and `ctl_floor.py` patch `fq4_split.tau_pred` in memory only.
