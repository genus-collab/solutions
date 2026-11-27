# Genus-4 split sum of the formal quintic: independent re-implementation

Date 2026-11-11. Written and run by an AI agent; nobody has checked this note. Exact arithmetic modulo 2^89-1 and 2^61-1, pure Python,
`py -I`. Longest single run: 4.3 minutes of wall time with 5 processes (13.4 minutes for one single-process brute-force check); no cluster
job was needed and none is prepared. Status labels: VERIFIED (exact finite computation; prime and orders stated), ARGUED, CONJECTURAL.
Finite computations are evidence, not proofs.

## 0. Verdict

The numbers produced by `fq4_split.py` / `genus4_split.py` at genus 4 are reproduced by a second code that shares no module with them
(VERIFIED, frame (0, 0), `K = 6`, i.e. the coefficients `[u^0..u^5]` of `G` in `u^-6 G(u)`, at 2^89-1 and at 2^61-1):

* point frame (`M = 1`): all 379 graph terms agree coefficient by coefficient with `g4_split_g4_{89,61}_id.json`; the total is
  `(1/180000) t_c^-6 + 0 t_c^-5 + ... + 0 t_c^-1`;
* natural `M = M_0(tau(u) zhat)`: the vertex graph is
  `u^-6 (1/180000 - 7u/300000 + 701u^2/18000000 - 1459u^3/45000000 + 923639u^4/67500000000 - 171583u^5/67500000000)`
  `= (1/180000) t_c^-6 + 0 + 0 + 0 + 0 + 0` (the six rationals reduce to the computed residues at both primes), and all 378 outer graphs
  with edges are `O(u^6)`; the table agrees with `g4_split_g4_{89,61}_nat.json` for all 379 outer indices;
* finer than the tables: every one of the 7146 non-zero pieces (inner graph `Gt` -> outer graph `Gt`/inner edges) agrees with the
  project's part files `parts/g4_{89,61}_nat_run.json`. These pieces are generic non-zero residues; the zeros of the final table are
  cancellations among them (268, 136 and 368 non-zero pieces for the three one-edge outer graphs).
* controls (2^89-1 only): the bare symbol `M = M_0` (regularity fails for 149 outer graphs, as it should) and `tau` with its `u^5`
  coefficient shifted by 1 (only the three one-edge graphs acquire a `[u^5]` term) reproduce `g4_split_g4_89_M0.json` and
  `g4_split_g4_89_nat_t5.json` for all 379 outer indices and all pieces (9858 and 7146 non-zero pieces).

No discrepancy and no bug was found. What this does and does not establish is in s.5.

## 1. What is independent, what is shared

| ingredient | project code (`fq4_split.py`) | this folder |
|---|---|---|
| hatted rows `Qhat[a][i][n](u)`, `xhat`, `t_c` | `fq_rows`, `fq_conifold`, `fq_sser`, `fq3_msum.HatRowsU` (s-series, then `u`) | `natural_form/referee/ref_{fp,rows,hat}.py` (WKB polynomials directly in `(zhat, u)`; written earlier by another agent, not by me) |
| `E_R(a, b)[i][j]` | `fq_sser.EdgeLeaf` (closed alternating sum) | `ind_tables.Tables.divide`: recursion `e(m, n) = d(m+1, n) - e(m+1, n-1)` with the exactness of the division by `z + w` asserted |
| leaves `T_R(e)` | `EdgeLeaf.T` | `(-1)^e Qhat_{0, e-1}` from the definition `T(z) = z(1 - R^-1(z) 1)` |
| `M_0`, `E_M` | colour matrices of `fq_symbol`, `fq3_mvertex.MVertex6.Eh` | Teleman form at `u = 0` from the rows, inverse by series inversion (symplecticity asserted), same division routine; `E_M = tau^(a+b+1) E_M0` asserted against the division applied to `M_0^-1(tau z)` |
| Witten-Kontsevich numbers | `fq_hodge.wk` (string equation + DVV on the smallest index) | `ind_wk.wk` (DVV on the largest index for all `k >= -1`) |
| stable graphs | `fq3_graphs.legged_graphs` (one-edge degenerations) | Mathematica export `cluster_genus4_hodge/StableGraph-4.bin`, own parser with stability/genus/connectivity checks |
| automorphisms | `fq_graphs.autos_of`: enumeration of (vertex permutation, edge permutation, loop flips) | depth-first count of half-edge bijections commuting with the edge involution |
| vertex leaf sum | tuples of leaves with a dimension bound | partitions of the deficit `3g_v - 3 + n_v - |hp|` |
| colour sum | vertex elimination in colour space | momentum space (s.2); plus a plain `5^V` brute force as a check |
| regrouping `(Gt, S) -> Gt`/inner | union-find + `fq3_graphs.canon` | own union-find + own canonical form |

Shared by construction: the method itself (the split sum, s.5 item 1), the referee rows (a second derivation, but not mine), and the
reference tables that I compare with. The project's genus-4 code is imported in exactly one place, `ind_conv.py`, which is a convention
check and feeds nothing into the runs: my tables `E_R` (45 tables x 25 entries x 6 coefficients), `T_R`, `E_M`, `E_out`, the valuations of
`E_out`, the powers of `xhat` and all 231 vertex vectors (19 vertex types) are entrywise EQUAL to those inside `fq4_split.Split`
(genus 4 at both primes, genus 3 at 2^89-1). `ind_wk.selftest` compares 142 numbers with `fq_hodge.wk` (and checks closed
formulas, `<tau_10>_4 = 1/7962624`, and 13 classical genus-2/3 values). The Mathematica list equals the project's list as a set of 379
canonical forms and all `|Aut|` agree (also with my own evaluation of the product formula; genus 2: 7 graphs, genus 3: 42).

## 2. The colour sum in momentum space

For a 5-vector `V[c]` and a 5x5 matrix `E[i][j]` put `Vh[k] = (1/5) Sum_c zeta^(-kc) V[c]`, `Eh[k][l] = (1/25) Sum_ij zeta^(-ki-lj) E[i][j]`.
All components are computed and the following sparsity is ASSERTED on every table (it is the Z/5 covariance of the rows, not assumed):
a vertex vector has exactly one non-zero component `kappa`, and an edge table `E(a, b)` only components with `k + l = sigma = -(a+b+1) mod 5`.
Then for a decorated graph

    Sum_colours prod_v V_v[c_v] prod_e E_e[c_u][c_w] = 5^V Sum_{k_e in Z/5} prod_e Eh_e[k_e][sigma_e - k_e] prod_v Vh_v[kappa_v],

the sum running over the solutions of `kappa_v + Sum_{e leaves v} k_e + Sum_{e enters v} (sigma_e - k_e) = 0 mod 5` at every vertex:
`5^(E-V+1)` solutions (free momenta on the edges outside a spanning tree, at most 625 at genus 4). The vertex series factor out of the
sum. For an inner edge `Eh` is a number times `tau^(a+b+1)`, so for a given set `S` of outer edges the sum over the momenta of the inner
edges is a sum of numbers; only the outer edges carry `u`-series, each divided by its `u`-valuation (1 for every `(a, b)`).
This is a different algorithm from the project's, and it is about three times cheaper (1100 s of CPU against 3250 s for the natural run).

Checks of this algorithm against the plain colour-space sum over `5^V` colourings (`Summer.graph_brute`), all EQUAL (VERIFIED, 2^89-1):
genus 2, all 7 graphs, three variants; genus 3, all 42 inner graphs, variants nat and id, all subsets; genus 4, variant id, all 379 graphs;
genus 4, variant nat, all subsets, the first 26 inner graphs (all 11 with at most 2 edges and 15 of the 21 with 3 edges); genus 4, variant nat, subsets with at most one outer edge,
24 inner graphs with 5 or 6 vertices (every sixth one, up to `E = 9`).

## 3. Results

Genus 2 (`K = 4`; VERIFIED, natural `M` at both primes, `M_0` and id at 2^89-1): vertex `[1/1200, -7/6000, 287/864000]`, six edge
graphs `O(u^2)`, the `[u^2]` values `479/756000, 1/2250, -127/2016000, 0, 0, 0` of the notes; with `M_0`: `[u^1] = 1/750, -13/8000` on
graphs 2, 4 (as in the notes).

Genus 3 (`K = 6`, both primes; VERIFIED): against the stored tables of the OLD vertex-tensor code (a different formulation of the outer
graph sum): point frame, all 42 graphs, `u^0..u^5` (`g3_state_*.json: A_id`); `M = M_0`, all 42 graphs, `u^0..u^5` (`A_M0`, 37 non-zero
rows); natural `M`, all 42 graphs, `u^0..u^3` (`g3_nat_*_nat.json: A3`; the old code used `M` mod `u^4`, so `u^4, u^5` are not
comparable). Natural `M`: vertex `[1/25200, -1/9000, 407/3780000, -151/3780000]` `= (1/25200) t_c^-4 + 0 + 0 + 0`, all 41 edge graphs `O(u^4)`.

Genus 4 (`K = 6`; VERIFIED at 2^89-1 and 2^61-1 unless stated):

| run | compared with | outcome |
|---|---|---|
| point frame (id) | `g4_split_g4_*_id.json`, 379 graphs x 6 coefficients; 379 pieces of `parts/g4_*_id_run.json` | all equal; total `= (1/180000) t_c^-6 + 0 + ... + 0` |
| natural `M` | `g4_split_g4_*_nat.json`, 379 outer indices; 7146 non-zero pieces of `parts/g4_*_nat_run.json` | all equal; vertex graph as in s.0; 378 edge graphs `O(u^6)`; total = point-frame total; for each of the 379 inner graphs the pieces add up to my point-frame term |
| control `M = M_0` (2^89-1) | `g4_split_g4_89_M0.json`, 379; 9858 pieces | all equal; 149 edge graphs not `O(u^6)` (every outer graph with 1 to 5 edges), vertex graph not of gap form |
| control `tau_5 -> tau_5 + 1` (2^89-1) | `g4_split_g4_89_nat_t5.json`, 379; 7146 pieces | all equal; `[u^5]` of the three one-edge graphs `= -1297957089/160000000000, -115412058321/16000000000000, 409242577/25000000000` (rational reconstructions at one prime), vertex graph acquires a `t_c^-1` term |

Outer graphs verified independently: ALL 379 (the vertex graph, the three one-edge graphs: bridge g1-g3, bridge g2-g2, loop on genus 3,
and every other one), at both primes; the list with the number of non-zero pieces per outer graph is `out/g4_{89,61}_nat_verified.txt`.
Not verified by me: the second conifold frame `(1/3, -2/7)`; anything beyond `u^5`.

Structure of the genus-4 cancellation (natural `M`, the same at both primes):

| outer edges | outer graphs | with a non-zero piece | non-zero pieces | non-zero sums |
|---|---|---|---|---|
| 0 | 1 | 1 | 379 | 1 |
| 1 | 3 | 3 | 772 | 0 |
| 2 | 7 | 7 | 1363 | 0 |
| 3 | 21 | 18 | 2197 | 0 |
| 4 | 43 | 24 | 1839 | 0 |
| 5 | 75 | 16 | 596 | 0 |
| 6..9 | 229 | 0 | 0 | 0 |

So the content of "all 378 outer graphs with edges are regular" is: 229 graphs are `O(u^6)` trivially (every outer edge is `O(u)`:
`[u^0] E_R = E_M0`, asserted for all 45 tables), 81 graphs with at most 5 edges have all their pieces equal to zero through `u^5`, and
68 graphs vanish by cancellation among 6767 non-zero pieces. With the bare symbol the same 149 graphs with at most 5 edges have 9479
non-zero pieces and none of the sums vanishes.

## 4. Timings (this machine, at most 5 workers + the parent)

| run | wall | CPU |
|---|---|---|
| genus 4, id, 2^89-1 / 2^61-1 | 33 s / 56 s (4 workers) | 154 s / 211 s |
| genus 4, natural, 2^89-1 / 2^61-1 | 224 s / 213 s | 1103 s / 1043 s |
| genus 4, `M_0` control; `tau_5` control (4 workers) | 170 s; 259 s | 668 s; 1023 s |
| brute force, id, 379 graphs (1 process) | 13.4 min | |
| genus 3, any variant | 2 s | 3 s |

`ind_run.py` stores every finished inner graph, stops itself at `--deadline` (default 1500 s) and resumes on the next call, so the
30-minute rule cannot be violated by it. Nothing here needs a cluster.

## 5. What the agreement means, and what it does not

1. (ARGUED, not re-derived numerically at genus 4.) Both codes compute the same object: the point-frame graph sum of `R` with every edge
   split as `E_R = E_M + E_out` and regrouped by contracting the inner edges. That this regrouped sum IS the outer graph sum
   `R'.(M.omega^)` rests on the composition law of the Givental action (`E_out = M^-1 E_R' M^-T`, `T_R = T_M + M^-1 T_R'`, and the orbit
   count `Sum_S 1/|Aut Gt| = Sum 1/|Aut(Gt, S)|`). The numerical support for the identification is at genus 2 and 3 only: there my code
   (and the project's) reproduces the tables of the old vertex-tensor code, which evaluates `R'.(M.omega^)` with genuine `M`-vertices
   (`M_0`: 42 graphs through `u^5`; natural `M`: through `u^3`). An independent re-implementation of the split sum cannot detect an
   error in this identification at genus 4; the argument is genus-independent, which is the reason to believe it.
2. The rows are those of the earlier referee, not derived again here. An error common to the project rows and the referee rows would
   pass unnoticed; the total `= (1/180000) t_c^-6 + 0 + ... + 0` (the conifold gap with the Bernoulli coefficient, five vanishing
   subleading coefficients) is a strong joint test of rows, vertices, graphs and automorphisms, at both primes.
3. The statements are modulo two primes and through `u^5` only, frame (0, 0). The zeros are exact residues; a false zero would need a
   numerator divisible by both primes.
4. Regularity at genus 4 concerns the 149 outer graphs with at most 5 edges (s.3). The `tau_5` control shows that `[u^5]` of the
   one-edge graphs is sensitive to `tau_5`; the sensitivity of the multi-edge graphs to the higher `tau_k` was not probed here.

## 6. Files and commands

| file | content |
|---|---|
| `ind_wk.py` | Witten-Kontsevich numbers, self-test |
| `ind_graphs.py` | Mathematica graph list, canonical forms, automorphisms, contraction, self-test |
| `ind_tables.py` | rows -> `E_R`, `T_R`, `M_0`, `E_M`, `E_out`, Fourier transforms, vertex tables (all assertions run at construction) |
| `ind_sum.py` | the split sum of one inner graph: momentum space (`graph`) and colour-space brute force (`graph_brute`) |
| `ind_run.py` | driver (pool, deadline, resume, comparison with the project tables, gap and regularity checks) |
| `ind_post.py` | piecewise comparison with the project part files, cancellation statistics, `t_c` expansions, list of verified graphs |
| `ind_conv.py` | convention check against the tables inside `fq4_split.Split` |
| `out/<tag>.json`, `<tag>_parts.json`, `<tag>.log`, `<tag>_post.log`, `<tag>_verified.txt`, `*_console.txt`, `conv_*.log` | results and logs; tags `g4_{89,61}_{id,nat}`, `g4_89_M0`, `g4_89_nat_t5`, `g3_{89,61}_{id,nat,M0}_K6`, `g2_{89,61}_nat` |

```bash
cd "G:/My Drive/Projects/Old/refGW/hodge_vertex_quintic/natural_form/genus4_referee/indep"
py -I ind_wk.py
py -I ind_graphs.py
py -I ind_conv.py --genus 4 --prime 89
py -I ind_run.py --genus 3 --variant id --K 6 --nproc 5      # then --variant nat, --variant M0 (id first: nat compares with it)
py -I ind_run.py --genus 4 --variant id --nproc 5
py -I ind_run.py --genus 4 --variant nat --nproc 5           # add --prime 61 for the second prime
py -I ind_post.py --prime 89 --variant nat
py -I ind_run.py --genus 4 --variant nat --brute 0,1,2,3     # colour-space brute force on the listed inner graphs
```
