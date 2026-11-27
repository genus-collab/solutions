# Adversarial code review of `fq4_split.py` / `genus4_split.py` (genus-4 edge-splitting run)

Date 2026-11-11. Reviewer: an AI agent; no human has checked this note. All computations are exact, modulo 2^89-1 and/or 2^61-1 as stated.
Status labels: VERIFIED (finite exact computation; what, which prime, which orders), ARGUED (reasoning written here), CONJECTURAL.
Finite computations are evidence, not proofs. Nothing outside this folder was written; the project files were only read
(`sys.dont_write_bytecode` is set in every script, so not even byte-code caches were touched).

## 0. Verdict

**The genus-4 conclusion of the logs of record stands as a finite exact statement, and it is not vacuous.** No bug was found in
`fq4_split.py` that affects it. Two real defects (B1, B2) and two weaknesses (W1, W2) were found in the version that was given for
review (section 6); none of them affects a stored number.

Versions. The review started on `fq4_split.py` of 08:51 and `genus4_split.py` of 08:44 (2026-11-11). Both files were changed by
their author at 09:40:48, while the review was running; the change is exactly: an assertion `min valO >= 1` in `Split.__init__`
(W1), a check of the same in `checks`, `len(T) == npol` in the GAP / VERTEX tests (W2), an exact pattern for the part files (B1) and a
`merge` that counts an inner graph once and asserts identical content (B2). So B1, B2, W1, W2 below describe the earlier version and
are FIXED in the present one (read, and re-tested: the genus-2/3 mutation regression was repeated on copies of the present source,
`out/mut_g3.log` second block; the full run of m07 used the present source and the present `checks`; the light scripts were re-run).
The logs of record were re-collected by the author at 09:41 with the fixed driver ("part files: ['g4_89_nat_run.json']", all checks
pass); the stored table is identical to the full run of my copy (`rv_cmp.py`: 69 rows, all residues equal).

Statement that the review supports (VERIFIED, 2^89-1 and 2^61-1, frame (0, 0); second frame (1/3, -2/7) at both primes by the logs of
record, re-checked here only through the stored tables): with `M = M_0(tau(u) zhat)` mod `(zhat^10, u^6)`,
`1/tau = 2F1(1, 4/5; 9/5; u)`, the outer graph sum `R'.(M.omega^)` at genus 4 has all 378 outer graphs with edges equal to `O(u^6)` and
the vertex graph equal to `(1/180000) (u/t_c)^6 mod u^6`, i.e. `u^-6 G = (1/180000) t_c^-6 + O(t_c^0)`.

How much of this is content (section 3.6): 229 of the 378 graphs (those with 6 or more edges) are regular for the trivial reason that
every outer edge is `O(u)` and are never computed; 81 graphs with 3, 4 or 5 edges receive only pieces that vanish identically (a
selection rule of `E_out`; a generic `E_out` of the same valuation gives non-zero pieces); the remaining 68 graphs (3, 7, 18, 24, 16 with
1, .., 5 edges) receive 6767 non-zero pieces which cancel. So the non-trivial part of "378 graphs are regular" is 68 graphs, 149 if the
selection rule is counted.

| item | verdict |
|---|---|
| (a) regrouping identity (edge splitting = outer graph sum, graph by graph, with automorphism factors, leaves, xhat powers) | ARGUED (section 2), under hypotheses H1-H4; its combinatorial core VERIFIED with weights `t^(inner edges)` for all 42 graphs of genus 3 and all 379 of genus 4 |
| (a') `E_out = M^-1 E_R' M^-T`, `T_R = T_M + M^-1 T_R'`, `M` symplectic, for the tables actually used at genus 4 | VERIFIED independently (own series arithmetic): 2^89-1 and 2^61-1 in the frame (0, 0), 2^89-1 in the frame (1/3, -2/7); `zhat <= 9`, `u <= 5` |
| (b) `contract`, `plan`/`colour_sum`, `vertex_table`, pruning, popcount restriction, truncation, `checks` | correct (ARGUED by reading; VERIFIED where a test was possible); weaknesses W1, W2 (fixed since), W3, W4; none was active in the runs of record |
| (c) mutation tests | 26 mutants: 23 detected by the genus-2/3 regression; the 3 survivors (`tau_4`, `tau_5`, tau-power at `zhat >= 7`: invisible at genus 3 by construction) are all detected by the genus-4 checks |
| (d) independence of the internal checks | the GAP check, the per-inner-graph additivity and `frozen_ok` do NOT test `E_M`'s `u`-dependence; only REGULARITY / VERTEX GRAPH (one check, not two) and the genus-3 regression do; an independent test of `E_M` was added (section 5) |
| driver bugs (version of 08:44) | B1 (part-file prefix collision in `--collect`): real, latent, would have overwritten a table of record on a repeated collect; B2 (part files not checked for overlap): real, latent. Both FIXED in the version of 09:40:48 |

## 1. What was done (files, all in this folder)

| file | content | result |
|---|---|---|
| `rv_common.py` | mutant factory: copies of `fq4_split.py` in `mutants/` with one textual replacement each (asserted to match exactly once) | 27 copies (baseline + 26 mutants) |
| `rv_mut_g3.py` | the original `genus4_split.regress` run on every mutant (both primes) | `out/mut_g3.log`, `out/mut_g3.json` |
| `rv_g4mut.py` | genus-4 runs of a mutant, judged by the original `genus4_split.checks` (full) or by exact restricted checks (outer graphs with at most 2 edges) | `out/g4mut_*.log`, `out/g4mut_*.json` |
| `rv_compose.py` | independent check of the composition law on the genus-4 tables, and of the natural-form support rule | `out/compose_89.log`, `out/compose_61.log`, `out/compose_89_ab.log` |
| `rv_parts.py` | read-only analysis of the part files of record; independent contraction on all 38811 pairs (inner graph, subset) | `out/parts.log` |
| `rv_mass.py` | regrouping identity with weights `t^(inner edges)`, exact rationals | `out/mass_g3.log`, `out/mass_g4.log`, `out/mass_g4b.log` |
| `rv_tc.py` | independent form of the GAP / VERTEX checks; demonstration of bug B1 on the present content of `parts/` | `out/tc.log` |
| `rv_vertex.py` | unit tests of `partitions_min2`, `sorted_tuples`, `vertex_table` against brute force | `out/vertex.log` |
| `rv_absent.py`, `rv_probe.py`, `rv_cmp.py`, `rv_table.py` | probes: valuations of `E_out`, the outer graphs that receive no non-zero piece, comparison of the copy's full table with the table of record, the mutation table | console only (quoted in 3.6, 4) |

Run times: the longest single run was 9 minutes (`rv_mass.py 4`, one process); the full genus-4 runs took 5 minutes each with 5 worker
processes; never more than 5 workers plus the idle parent, plus at most one light single-process script.
`mutants/` now holds copies of the present source (09:40:48); the genus-4 mutant runs of sections 4.2, 4.3 other than the full run of m07
were made with copies of the 08:51 source, which differs only by the absence of the new assertion.

## 2. (a) The regrouping identity

### 2.1 Setting and definitions

Ground ring `A = F_p[u]/(u^K)` (the argument is the same over `Q[[u]]`). `V = A^5` with the orthonormal basis `e~_j` (normalised canonical
basis). All operator series `X(z)` below are in `End(V)[[z]]` with `X(0) = 1` (`z` stands for `zhat`), and are used modulo `z^(3g-2)`.

* Trivial theory with dilaton-shifted Kontsevich-Witten vertices: for a vertex of genus `h` and colour `j` with `n` insertions
  `e~_j psi^(a_1), .., e~_j psi^(a_n)`: `(zeta^(4j) xhat)^(2h-2+n) <tau_(a_1) .. tau_(a_n)>_h`; mixed colours give 0. Unit `1 = Sum_j 1_j e~_j`,
  `1_j = zeta^(-4j)/xhat` (in the code: `1_j = Qhat_(0,0)^(j)`, checked in `rv_compose.py` (4)).
* For symplectic `X` (`X^-1(z) = X(-z)^T`): edge `E_X(z, w) = (1 - X^-1(z) X^-1(w)^T)/(z + w)`, leaf `T_X(z) = z (1 - X^-1(z) 1)`, which is
  `O(z^2)`.
* Point-frame sum of `X` (the Givental action on the product of KW potentials, in Feynman form): for a leg-free stable graph `Gt`,
  `Cont_X(Gt) = (1/|Aut Gt|) Sum_(decorations) prod_v [ Sum_m (1/m!) KW vertex with the half-edge insertions and m leaves T_X ] prod_e E_X`,
  where a decoration is a colour per vertex and an exponent per half-edge, and `F_g(X) = Sum_(Gt) Cont_X(Gt)`.
  `Aut Gt` is the automorphism group of the half-edge structure (vertex permutations preserving genera, permutations of parallel edges,
  loop flips): this is what `fq3_graphs._graph_dict` counts.
* Vertex theory `Omega^M = M.omega^` with legs: `Omega^M_(h,n)(e~_(i_1) psi^(d_1), .., e~_(i_n) psi^(d_n))` := the sum over stable graphs
  `gamma` of type `(h, n)` with LABELLED legs of `(1/|Aut gamma|)` (automorphisms fixing the legs) times KW vertices, edges `E_M`, leaves
  `T_M`, and on leg `l` at the vertex `w` of colour `j` the insertion `Sum_k (M^-1_k)_(j, i_l) e~_j psi^(d_l + k)`.
  This is the definition used by the old vertex-tensor code (`fq3_mvertex.MVertex6.tensor`).
* Outer sum: `F_g(R'; Omega^M) = Sum_Gamma (1/|Aut Gamma|) < prod_v Sum_m (1/m!) Omega^M_(g_v, n_v + m)(.., T_R'^m), prod_e E_R' >`.

Hypotheses.
* H1 `R = R' M` with `R'`, `M` symplectic, all equal to 1 at `z = 0`. (For the tables used: VERIFIED, section 5.)
* H2 The vertex tensors of the outer sum are the ones defined above (correlators of the potential `M^ D_triv`, descendant insertions
  passing through `M^-1(psi)`), not some other normalisation of `M.omega^`. This is a convention; it is the one of the old code, so the
  regression of section 4 tests that the two codes mean the same object.
* H3 `|Aut|` and `m!` are invertible in `A` (true: `p` is a 61- or 89-bit prime), and all sums are finite (dimension constraints).
* H4 Leaves are `O(z^2)` (true by H1: `T_X(z) = z(1 - X^-1(z) 1)`).

### 2.2 Claim

For every `g >= 2`, with `E_R = E_M + E_out`, `E_out := E_R - E_M`:

    F_g(R) = Sum_(Gt) (1/|Aut Gt|) Sum_(S subset E(Gt)) W(Gt, S)        (W: the decoration sum with E_out on S, E_M off S, leaves T_R)
           = Sum_Gamma [ Sum over pairs (Gt, S) with Gt/(E(Gt) - S) isomorphic to Gamma of W(Gt, S)/|Aut Gt| ]
    and the bracket equals the outer term (1/|Aut Gamma|) < prod_v Sum_m (1/m!) Omega^M_(g_v, n_v+m)(.., T_R'^m), prod_e E_R' >  of Gamma.

### 2.3 Proof (ARGUED)

Step 0 (algebra, uses H1). `R^-1 = M^-1 R'^-1`, so
`1 - R^-1(z) R^-1(w)^T = (1 - M^-1(z) M^-1(w)^T) + M^-1(z) (1 - R'^-1(z) R'^-1(w)^T) M^-1(w)^T`, i.e.

    E_R(z, w) = E_M(z, w) + M^-1(z) E_R'(z, w) M^-1(w)^T,     T_R(z) = z(1 - M^-1 1) + M^-1(z) z (1 - R'^-1(z) 1) = T_M + M^-1 T_R'.

So `E_out = M^-1 E_R' M^-T`: coefficientwise
`E_out(a, b)[i][j] = Sum_(k+m=a, l+n=b) Sum_(i', j') (M^-1_k)_(i i') E_R'(m, n)[i'][j'] (M^-1_l)_(j j')`.

Step 1 (expansion). `prod_e E_R = Sum_S prod_(e in S) E_out prod_(e not in S) E_M` on a fixed representative of `Gt` (edges labelled). `W(Gt, S)`
is invariant under isomorphisms of marked graphs because all decorations are summed.

Step 2 (orbit-stabiliser). `Aut Gt` acts on the subsets of `E(Gt)`; the orbit of `S` has `|Aut Gt|/|Aut(Gt, S)|` elements, so
`(1/|Aut Gt|) Sum_S W(Gt, S) = Sum_(orbits [S]) W(Gt, S)/|Aut(Gt, S)|`, and
`F_g(R) = Sum over isomorphism classes of marked graphs (Gt, S) of W(Gt, S)/|Aut(Gt, S)|`.

Step 3 (nested graphs). A marked graph `(Gt, S)` is the same thing as: an outer graph `Gamma = Gt/(E - S)` (vertices = connected components
of the inner subgraph `(V, E - S)`, genus of a component = sum of its vertex genera + its first Betti number, half-edges = half-edges of
the edges in `S`) together with, for every vertex `v` of `Gamma`, a connected stable graph `gamma_v` of type `(g_v, n_v)` whose legs are
labelled by the set `H_v` of half-edges of `Gamma` at `v` (stability is inherited vertex by vertex). Let `D(Gamma)` be the set of families
`d = ([gamma_v])_v` of isomorphism classes of `H_v`-labelled graphs. `Aut Gamma` acts on `D(Gamma)` by transport of labels, and the
isomorphism classes of marked graphs with contraction `Gamma` are the orbits. Restriction to the outer half-edges gives an exact sequence

    1 -> prod_v Aut(gamma_v) -> Aut(Gt, S) -> Stab_(Aut Gamma)(d) -> 1

(kernel: automorphisms fixing every outer half-edge fix every component and its legs; when `S` is empty there is one component and the
statement is trivial; surjectivity: an element of the stabiliser lifts by choosing label-compatible isomorphisms
`gamma_v -> gamma_(sigma v)`). Hence `|Aut(Gt, S)| = |Stab(d)| prod_v |Aut gamma_v|`, and, the orbit of `d` having `|Aut Gamma|/|Stab(d)|`
elements,

    Sum over (Gt, S) with contraction Gamma of W/|Aut(Gt, S)|  =  (1/|Aut Gamma|) Sum_(d in D(Gamma)) W(Gamma, d)/prod_v |Aut gamma_v|.

Step 4 (the weight factorises into vertex tensors). Insert Step 0 on every edge of `S`: an outer half-edge `h` at the inner vertex `w`
(colour `j`, exponent `a`) becomes a leg of `gamma_v` with insertion `e~_(i') psi^m` and the leg rule `(M^-1_k)_(j i')`, `a = m + k`: exactly
the leg rule of `Omega^M`. Leaves: at an inner vertex, `Sum_m (1/m!) T_R^m = Sum_(m1, m2) (1/(m1! m2!)) T_M^(m1) (M^-1 T_R')^(m2)`; the
`m2` leaves of the second kind are extra legs carrying `T_R'` through the leg rule. For fixed `Gamma`, `v` and a total number `r` of such
extra legs on `gamma_v`: the graphs `gamma'` of type `(g_v, n_v + r)` all of whose vertices stay stable when the `r` extra legs are removed
correspond to pairs (`gamma_v`, a map `f` from the `r` extra legs to `V(gamma_v)`) modulo `Aut gamma_v`, with `Aut gamma' = Stab(f)`, so
`Sum_(gamma') (..)/|Aut gamma'| = Sum_(gamma_v) (1/|Aut gamma_v|) Sum_f (..)`, and `(1/r!) Sum_f = Sum_((r_w)) prod_w 1/r_w!`: the factors
`1/m2!` per inner vertex. The remaining graphs `gamma'` (a genus-0 vertex with `k <= 2` half-edges or legs that is stable only thanks to
`r_w >= 1` extra legs) contribute 0: the KW number lives on a space of dimension `k + r_w + m1 - 3` while the leaves alone have degree
`>= 2(r_w + m1)` (H4). Therefore

    Sum_(d) W(Gamma, d)/prod_v |Aut gamma_v| = < prod_v Sum_r (1/r!) Omega^M_(g_v, n_v + r)(.., T_R'^r), prod_e E_R' >,

which with Step 3 is the claim.

Step 5 (xhat powers and phases). The code puts `(zeta^(4j) xhat)^(2 g_w - 2 + n_w + m_w)` on every INNER vertex (`n_w` = all half-edges,
`m_w` = all leaves of `T_R`). Over a component this multiplies to `xhat^(2 g_v - 2 + n_v + m_total)`. The old code puts
`xhat^(2 g_v - 2 + n_v + r)` on the outer vertex (`r` = leaves of `R'`) and uses inner leaves `That` built with the unit of the theory with
`xhat = 1`. The two agree because `T_M = -M^-1_(d-1) 1` is linear in the unit and `1_j = zeta^(-4j)/xhat`: each inner leaf carries `1/xhat`,
which cancels its `xhat` from the vertex power. (VERIFIED: `rv_compose.py` (4) uses `1_j = Qhat_(0,0)^(j)`; the regression compares the
two normalisations on 42 graphs.)

Remarks. (i) The term `S = {}` has a single component: it is the vertex graph, with edges `E_M` and leaves `T_R = T_M + M^-1 T_R'`, as it
must be. (ii) The argument does not use any property of `tau`; it holds for every symplectic `M`. If the table `E_M` is NOT of the form
`(1 - M^-1 M^-T)/(z + w)` for a symplectic `M` (mutants m06, m07, m08), the code still computes a decomposition, but it is not an outer
graph sum of anything: that the table used has this form is a separate fact, VERIFIED in section 5.

### 2.4 Numerical test of Steps 2-3 (VERIFIED, exact rationals)

`rv_mass.py`: with every vertex, leaf and outer-edge weight replaced by 1 and every inner edge by `t`, the claim reads

    Sum_((Gt, S) -> Gamma) t^(|E(Gt)| - |S|)/|Aut Gt| = (1/|Aut Gamma|) prod_v N_(g_v, n_v)(t),   N_(h,n)(t) = Sum_(gamma) t^|E(gamma)|/|Aut gamma|.

Result: equality of polynomials for all 42 outer graphs of genus 3 (865 pairs) and for all 379 of genus 4 (38811 pairs; two passes,
370 graphs in `out/mass_g4.log` and 375 in `out/mass_g4b.log`, the untested sets being disjoint: those with a vertex of type (2,3), (1,6),
(3,1), (2,4), (3,2), (4,0) in the first pass and (0,7), (0,8) in the second). This tests the weight `1/|Aut Gt|` per labelled subset, the
contraction, the automorphism numbers and the completeness of the lists (e.g. `N_(4,0)` has top coefficient 1105/1152 and `N_(3,0)` 5/16,
the known weighted counts of trivalent graphs). It does not test Step 4 (legs and leaves), which is covered by the regression.

## 3. (b) Code review, item by item

### 3.1 `contract(G, S)` -- correct (VERIFIED)

Union-find over the edges NOT in `S`; component genus `Sum g_v + (#inner edges) - (#vertices) + 1` (an inner loop adds 1, an inner
multi-edge adds its extra edges); outer edges are mapped to component pairs (an outer edge inside one component becomes a loop); the key
is `fq3_graphs.canon` (minimum over all vertex permutations of (genera, sorted edge multiset)), the same function that produces the keys
of the list, so equal keys <=> isomorphic outer graphs (as genus-labelled multigraphs, which is what a leg-free stable graph is).
`rv_parts.py` (3): equal to an independent contraction (edges contracted one at a time on an explicit multigraph) on all 38811 pairs
(inner graph, subset) of genus 4; every key is one of the 379, has `|S|` edges and total genus 4; `S = {}` gives the vertex graph only and
`S` = all edges gives the graph itself; all 379 outer graphs are reached. A wrong key cannot pass silently: `outer_pos[...]` raises
`KeyError` (mutants m10, m12 crash at once).

### 3.2 `plan()` / `colour_sum()` -- correct (ARGUED; VERIFIED by regression and mutants)

Half-edge `2e` sits at the vertex `ed[e][0]`, `2e+1` at `ed[e][1]` (`fq3_graphs._graph_dict`: `H = [(e, s)]`, `hvtx = ed[e][s]`), so
`keys[e] = (hp[2e], hp[2e+1])` pairs the exponent `a` with the vertex `ed[e][0]`, and the tables are `E(a, b)[colour of ed[e][0]][colour of
ed[e][1]]` (`fq_sser.EdgeLeaf.E(m, n)[i][j]` is the coefficient of `psi^m` on the side of `i` and `psi^n` on the side of `j`;
`MVertex6.Eh[(a, b)][j][jj]` likewise). In `plan`, an edge between the eliminated vertex `v` and a frontier vertex is recorded with
`flip = True` iff `v = ed[e][0]`, and `colour_sum` then reads `mats[e][c_v][c_frontier]`, else `mats[e][c_frontier][c_v]`: consistent. A
vertex leaves the frontier only when all its neighbours are done, so every non-loop edge is contracted exactly once (when its second
endpoint is eliminated); loops use the diagonal entry; parallel edges are separate entries of `conn`. The elimination order affects only
the cost.
The symmetry `E(a, b)[i][j] = E(b, a)[j][i]` means that only the relative convention matters; mutants m13 (colours transposed), m20
(exponents swapped), m25 (loop off the diagonal) and m19 (only one assignment of exponents to the half-edges of a vertex) each break all
ten regression lines, including the `id` variant.

### 3.3 `vertex_table()` -- correct (VERIFIED by a unit test)

`rv_vertex.py` (2^89-1, the genus-4 tables): `partitions_min2` and `sorted_tuples` equal brute-force enumerations (no duplicates);
`vertex_table(g, hs)` equals a brute-force sum over ORDERED leaf tuples with `1/m!` (own series arithmetic) for all 19 vertex types that
occur at genus 4, 3496 ordered exponent tuples x 5 colours x 6 orders in `u`; this covers the multiplicity factor `1/prod(mult!)`, the phase
`zeta^(4 j n)` and `xhat^n`, `n = 2g - 2 + hs + m`. The cap `KZ + 1 = 10` on a leaf exponent is exactly the largest exponent that can occur
(the vertex graph with one leaf `psi^10`); the largest `xhat` power used is 15 (table up to 19). The range `m <= 3g - 3 + hs` is exact
(leaves have exponent `>= 2`). What the unit test cannot see (a wrong formula shared by both sides: the sign or definition of `T_R`, the
phase convention) is covered by the `id` variant of the regression and by the genus-4 GAP check (mutants m09, m14, m15, m18, m26).

### 3.4 Valuation pruning and the restriction `popcount < K` -- valid; its precondition was unchecked in the 08:51 version (W1, fixed since)

Every factor of a term is a power series in `u` (lists `u^0..u^(K-1)`, no negative powers), so a term with outer-edge set `S` has valuation
at least `Sum_(e in S) valO(key_e)`; skipping it when this is `>= K` is exact, and `valO` (minimum over the 25 colour entries) is
conservative. The restriction to `popcount(S) < K` is exact iff every `E_out` entry has valuation `>= 1`. The 08:51 version did not assert this (the present one does);
`frozen_ok` tests `[u^0] E_R = Eh_0`, which implies it only when `[u^0] E_M = Eh_0`, i.e. `tau_0 = 1` (true for `nat`, `M0`, `nat_t<k>` with
`k >= 1`; false for `nat_t0`, where the PASS line "so every outer edge is O(u)" would be printed although it is false -- seen with mutant m01:
`frozen_ok = True`, `min valO = 0`). In the runs of record: `valO = 1` for all 45 keys (`rv_probe.py`, `rv_compose.py` (5)), so the
restriction is valid, and the valuation pruning is inert (it never removes more than the popcount restriction).
That nothing non-zero below `u^6` was dropped is also seen directly: for each of the 379 inner graphs the retained pieces add up to the
point-frame term through `u^5` (`rv_parts.py` (1), both primes, 7162 pieces).

### 3.5 Truncation of series -- correct

`mul` truncates at `u^K`; `Split.__init__` asserts that every input series has precision `>= K` (`lst`), that the hatted rows are pure
`u`-series of valuation 0 (`H.val_ok`), and builds the rows with `prec_abs = 5K + 4`, `Nu = K + 4`. Accumulation in `graph()` is in
unreduced integers and reduced once. `E_R`, `T_R` recomputed with own arithmetic from the rows agree with the tables (`rv_compose.py` (2), (4)).

### 3.6 `checks()` -- correct for the runs of record; weaknesses W2 (fixed since), W3

* GAP: `tc_expansion3` composes `G(u)` with the inverse series of `t_c(u)` and multiplies by `u(t_c)^-6`. Independent form
  (`rv_tc.py`, own arithmetic): `u^-6 G = c t_c^-6 + O(1)` iff `G = c (u/t_c)^6 mod u^6`; the stored totals of all eight result files
  present (`nat`, `id`, both primes; `nat` in the second frame, both primes; the controls `M0`, `nat_t5` at 2^89-1) equal
  `(1/180000)(u/t_c)^6`, and the vertex rows of the four `nat` files equal the same series, with no other non-zero row.
  (`t_c(u) = u + 7u^2/10 + 41u^3/75 + 9u^4/20 + ..` in the frame (0, 0), from `fq_conifold`; not compared with the literature here.)
* REGULARITY tests `[u^0..u^5]` of every row with index `!= 0` present in the table; a row that is absent has no non-zero piece.
  Index 0 is the vertex graph (the list is sorted by number of edges; checked).
* What "378 graphs are O(u^6)" contains (2^89-1 and 2^61-1, identical counts; `out/parts.log`, `rv_absent.py`):

  | outer graphs | number | why regular |
  |---|---|---|
  | 6 or more edges | 229 (89 + 81 + 42 + 17) | not computed: every outer edge is `O(u)` (ARGUED from `valO >= 1`, which is VERIFIED) |
  | 3, 4, 5 edges, no non-zero piece at all | 81 (3 + 19 + 59) | every single colour sum vanishes identically; with a generic `E_out` of valuation 1 the same pieces are non-zero at `u^3, u^4, u^5` (3 graphs tried), so this is a selection rule of the true `E_out` (typically graphs with several one-valent vertices), not a dead code path |
  | 1..5 edges with non-zero pieces that cancel | 68 (3 + 7 + 18 + 24 + 16) | 6767 non-zero pieces from the 379 inner graphs; lowest non-zero order of a piece is `u^E` for 35 graphs, `u^(E+1)` for 26, `u^(E+2)` for 7 |

  The five-edge graphs test only the first-order data of `E_out` (their only possible coefficient is `u^5`, a product of five leading
  terms); `tau_4` is tested by the `u^4`, `u^5` coefficients of the three one-edge graphs and the `u^5` coefficient of five two-edge
  graphs, `tau_5` only by the `u^5` coefficient of the three one-edge graphs (section 4.2).
* W2 (08:44 version; fixed since): `checks` did not verify that `tc_expansion3` returned all six polar coefficients: `not any(T[1:npol])` passes on a shorter list. In
  the logs of record six coefficients are printed, so the check was complete.
* W3: VERTEX GRAPH and REGULARITY are one check, not two (vertex = total - edge graphs), and GAP is independent of `E_M` (section 5).

## 4. (c) Mutation tests

Method: `rv_common.py` writes copies of `fq4_split.py` to `mutants/`, each with ONE textual replacement (asserted to match exactly once in
the original source; m11 has two coupled replacements). Every copy, including `baseline`, also carries a non-mutating patch (a module global
`POPMAX`, default 99 = no effect) used for the restricted genus-4 runs. `baseline` passes everything, so the harness itself is sound.

### 4.1 Genus-2/3 regression (the original `genus4_split.regress`, both primes; `out/mut_g3.log`)

Ten lines per mutant (genus 3: `id`, `M0` through `u^5`, `nat` through `u^3`, 42 graphs each; genus 2: `nat`, `M0`; two primes). The two
primes always agree. "g3 id" etc. = the lines that FAIL.

| mutant | change | genus-2/3 regression | genus 4, 2^89-1 (4.2, 4.3) |
|---|---|---|---|
| baseline | exact copy | all 10 pass | `id` full: all pass; `nat` restricted: all pass; `nat` full: all pass, table = table of record |
| m01_EM_sign | `E_M -> -E_M` | FAIL: g3 M0, g3 nat, g2 nat, g2 M0 | restricted: REGULARITY, VERTEX fail |
| m02_tau1 | `1/tau`: `u^1` coefficient 4/9 -> 4/10 | FAIL: g3 nat, g2 nat | restricted: REGULARITY, VERTEX fail |
| m03_tau3 | `1/tau`: `u^3` coefficient 4/19 -> 4/20 | FAIL: g3 nat (graphs 1, 2, 5) | - |
| **m04_tau4** | `1/tau`: `u^4` coefficient 4/24 -> 4/25 | **all 10 pass** | restricted and full: REGULARITY, VERTEX fail; GAP passes |
| **m05_tau5** | `1/tau`: `u^5` coefficient 4/29 -> 4/30 | **all 10 pass** | restricted: REGULARITY, VERTEX fail |
| m06_taupow | `E_M(a,b) = tau^(a+b) Eh_0` | FAIL: g3 nat, g2 nat | restricted: REGULARITY, VERTEX fail |
| **m07_taupow_hi** | wrong tau power only for `a+b+1 >= 7` | **all 10 pass** | restricted and full: REGULARITY, VERTEX fail; GAP passes |
| m08_EM_transpose | `E_M(a,b)[i][j] -> [j][i]` | FAIL: g3 M0, g3 nat, g2 nat, g2 M0 | restricted: REGULARITY, VERTEX fail |
| m09_no_mult | drop `1/prod(mult!)` of equal leaves | FAIL: all | `id` full: GAP fails |
| m10_contract_loop | inner loop does not raise the component genus | crash `KeyError` | - |
| m11_contract_inverted | the OUTER edges are contracted (keys stay valid) | FAIL: g3 M0, g3 nat, g2 nat, g2 M0 | restricted: REGULARITY, VERTEX fail |
| m12_contract_gswap | component genera attached in reversed order | crash `KeyError` | - |
| m13_flip | `colour_sum`: edge matrix transposed | FAIL: all | - |
| m14_phase | vertex phase `zeta^(4j(n+1))` | FAIL: all | - |
| m15_xpow | `xhat^(n-1)` | FAIL: g3 id, g3 M0, g3 nat | `id` full: GAP fails |
| m16_popcount | subsets with `popcount < K-1` only | FAIL: g3 M0 (at `u^5`) | not run (see below) |
| m17_valprune | valuation pruning at `K-1` | FAIL: g3 M0 (at `u^5`) | not run (see below) |
| m18_leafmax | largest leaf exponent `3g-2` dropped | FAIL: g3 id, g3 M0, g3 nat, g2 nat | `id` full: GAP fails |
| m19_perms | only the sorted exponent assignment at a vertex | FAIL: all | - |
| m20_keys_swap | edge exponents swapped relative to colours | FAIL: all | - |
| m21_no_aut | drop the automorphism factor | FAIL: all | - |
| m22_EO_sign | `E_out = E_R + E_M` | FAIL: g3 M0, g3 nat, g2 nat, g2 M0 | - |
| m23_Minv_sign | `M_0^-1_k = C_k^T` without `(-1)^k` | crash: assertion `M_0 M_0^-1 != 1` | - |
| m24_tau_reciprocal | `tau = 2F1` instead of `1/2F1` | FAIL: g3 nat, g2 nat | - |
| m25_loopmat | loop edge off the diagonal | FAIL: all | - |
| m26_leaf_from3 | the `psi^2` leaf dropped | FAIL: all | `id` full: GAP fails |

The table is for copies of the 08:51 source. On copies of the present source (second block of `out/mut_g3.log`; first block saved as
`out/mut_g3_v0851.log`) the outcome is the same, except that m01, m08 and m22 now stop in the constructor at the new assertion
(`E_out is not O(u)`), which is a detection as well.

Result: 23 of 26 mutants are detected by the regression. The three survivors are exactly the ones that cannot be seen at genus 3 by
construction (`tau_4`, `tau_5` enter at `u^4`, `u^5`, and the stored `nat` table is valid through `u^3`; `zhat`-orders 7..9 do not exist at
genus 3). So the genus-3 regression validates the algorithm, the conventions and `tau_1..tau_3`, but NOT `tau_4`, `tau_5` nor the orders
7..9 of `M_0`; those rest on the genus-4 checks themselves and on section 5.

m16 / m17 (over-pruning by one order) are detected at genus 3 only by the `M0` variant at `u^5`. At genus 4 they would be invisible in
the final table if the result of record is true (the dropped five-edge subsets sum to zero per outer graph); they are visible in the
per-inner-graph additivity (pieces of an inner graph = its point-frame term, section 3.4), which the real code satisfies.

### 4.2 Genus 4, restricted runs (`rv_g4mut.py <mutant> nat 89 2`, about 2-3 minutes each with 5 processes)

Only the subsets with at most two outer edges are summed. This is exact for the vertex graph and the 10 outer graphs with one or two
edges; the control (baseline) reproduces these 11 rows of the table of record. Verdict = REGULARITY of the 10 graphs through `u^5`, vertex
graph in `t_c`, equality with the stored rows.

| mutant | non-regular outer graphs (of 10) | orders | vertex graph in `t_c` |
|---|---|---|---|
| baseline | 0 | - | `(1/180000) t_c^-6 + 0 + 0 + 0 + 0 + 0` |
| m04_tau4 | 8: the three one-edge graphs, five two-edge graphs | `u^4, u^5` (one edge), `u^5` (two edges) | `t_c^-2`, `t_c^-1` non-zero |
| m05_tau5 | 3: the three one-edge graphs | `u^5` | `t_c^-1` non-zero |
| m07_taupow_hi | 10 | from `u^1` (one edge), `u^2` or `u^3` (two edges) | `t_c^-5 .. t_c^-1` non-zero |
| m01, m02, m06, m08, m11 | non-regular | - | polar part wrong |

So the genus-4 checks are sensitive to `tau_4` and `tau_5` separately and to the orders `zhat^7..zhat^9` of the vertex symbol: the
genus-4 statement is a test of the natural form beyond what genus 3 could test. (The project's own controls `nat_t5` and `M0`, finished
during this review, agree: `nat_t5` has exactly 3 non-regular rows, `M0` has 149.)

### 4.3 Genus 4, full runs judged by the original `genus4_split.checks` (5 processes)

| run | wall, CPU | coverage | frozen | GAP | REGULARITY | VERTEX GRAPH |
|---|---|---|---|---|---|---|
| m04_tau4, `nat` | 299 s, 1454 s | PASS | PASS | **PASS** | FAIL: graphs 1, 2, 3 at `u^4, u^5`; 4, 7, 8, 9, 10 at `u^5` | FAIL |
| baseline, `nat` | 298 s, 1448 s | PASS | PASS | PASS | PASS | PASS |
| m07_taupow_hi, `nat` (present source and `checks`) | 287 s, 1398 s | PASS | PASS (and `min valO >= 1` PASS) | **PASS** | FAIL: 21 graphs (3 with one edge from `u^1`, 7 with two edges from `u^2`, 11 with three edges from `u^3`) | FAIL |
| baseline, `id` | 20 s, 65 s | PASS | - | PASS | - | - |
| m09, m15, m18, m26, `id` | about 20 s each | PASS | - | FAIL | - | - |

The full run of m04 finds the same eight non-regular graphs as the restricted run and no non-regular graph with three or more edges.
The table of the full baseline run of the copy equals the table of record (69 rows, all residues; `rv_cmp.py`).

## 5. (d) Independence of the internal checks: what still passes if `E_M` is wrong

`E_out` is DEFINED as `E_R - E_M`, so the sum over all subsets telescopes back to `E_R` for any table `E_M` whatsoever.

| check | still passes with a wrong `E_M`? | evidence |
|---|---|---|
| assertion `M_0 M_0^-1 = 1` | yes for any error in `tau` or in the assembly of `E_M`; it tests only that the colour matrices of the symbol are symplectic | m23 crashes, m01-m08 do not |
| `frozen_ok` (`[u^0] E_R = Eh_0`) | yes: it compares `E_R` with `Eh_0`, NOT with the table `E_M` that is used; it is blind to `tau`, to the sign and to a transposition of `E_M` | m01: `frozen_ok = True` while `min valO = 0`; m04 full run: PASS |
| coverage of the part files | yes | - |
| GAP (total over all outer graphs) | yes whenever `[u^0] E_M = [u^0] E_R` (then nothing is dropped by the popcount restriction); otherwise it fails, but for the wrong reason (truncation) | m04, m07 full runs: GAP PASS |
| per-inner-graph additivity (nat pieces = point-frame term) | yes, same as GAP | by construction |
| genus-3 regression, `id` | yes (`E_M` is not used) | m01-m08 pass `g3 id` |
| genus-3 regression, `M0` (through `u^5`) | no for errors in `Eh_0` or its conventions at `zhat <= 6`; yes for errors in `tau` | m01, m08 fail; m02-m07 pass |
| genus-3 regression, `nat` (through `u^3`) | no for `tau_1..tau_3` and `zhat <= 6`; yes for `tau_4`, `tau_5`, `zhat >= 7` | m02, m03, m06 fail; m04, m05, m07 pass |
| REGULARITY / VERTEX GRAPH at genus 4 | no (all `E_M` mutants fail) -- but these two are ONE check, and it is the claim itself | 4.2, 4.3 |

Consequence. Before this review, the only things tying the genus-4 table `E_M` to "the edge of `M_0(tau zhat)`" at the orders that genus 3
cannot see were the two source lines that build it and the outcome of the genus-4 run itself. A reader could object that a different
table might also give regular graphs. `rv_compose.py` closes this (VERIFIED, 2^89-1 and 2^61-1 in the frame (0, 0), 2^89-1 in the frame
(1/3, -2/7); `zhat <= 9`, `u <= 5`; own series arithmetic, own division by `psi + psi'`, `tau` from the Pochhammer ratio):

1. `M^-1_k := tau^k (-1)^k C_k^T` satisfies `M M^-1 = 1` and `M^-1(z) M^-1(-z)^T = 1`: `M` is symplectic.
2. `E_M` recomputed from this `M^-1` equals `Split.EM`; `E_R` recomputed from the rows equals `Split.ER` (45 keys x 25 colours x 6 orders).
3. With `R' := R M^-1` (rows `Q' = Q M^-1`), the edge `E_R'` recomputed from `Q'` and transported, `M^-1 E_R' M^-T`, equals `Split.EO`.
   So the pieces labelled "outer" really are the edges of `R' = R M^-1` and the pieces labelled "inner" the edges of `M`.
4. `T_R = T_M + M^-1 T_R'` for `d = 2..10`, with `T_M = -M^-1_(d-1) 1`, `1_j = Qhat_(0,0)^(j)`; `T_R'(1) = 0`.
5. `[u^0] E_R' = 0`, `[u^0] T_R' = 0`; every entry of `E_out` has valuation exactly 1.
6. The defining property of the natural factorisation (NATURAL_FORM.md s.1.1) holds on these tables: `G = Psi(u)^-1 R'` has
   `[u^k zhat^n] G = 0` for `n > floor((5k+2)/4)`, `k = 1..5`, `n <= 9` (degrees 1, 3, 4, 5, 6), and it FAILS for the mutated `tau` of m04
   (at `u^4 zhat^6..9`, `u^5 zhat^7..9`) and m05 (at `u^5 zhat^7..9`). This is an independent witness, not a graph sum, that `tau_4`, `tau_5`
   and `M_0` at `zhat^7..9` are the intended ones.

What remains shared between the code under review, the old code and these checks: the rows (`fq_rows`, `fq_conifold`, `fq_sser`,
`HatRowsU`), the symbol matrices `C_k` (`fq_symbol`), the KW numbers (`fq_hodge.wk`) and the graph generator (`fq3_graphs`). Evidence for
them inside this review: the total equals the gap at genus 4 (six conditions on the sum of 379 graphs; it would not survive an error in
rows, KW numbers, automorphism numbers or the list: mutants m09, m15, m18, m21, m26), the weighted counts of section 2.4, and item 2
(the tables are what the rows give). They were not rebuilt from scratch here.

## 6. Defects and weaknesses found

None of them changes a number in the logs or tables of record. B1, B2, W1, W2 were found in the version given for review
(`genus4_split.py` of 08:44, `fq4_split.py` of 08:51) and are FIXED in the version of 09:40:48 (changed by the author during the
review; the fixes were read and are adequate).

**B1 (BUG, driver; fixed).** Collect step of the 08:44 version:
`parts = [.. for f in sorted(os.listdir(parts)) if f.startswith(tag + "_") and f.endswith(".json")]`.
The tag of a base run is a prefix of the tags of its variants: `g4_89_nat_` is a prefix of `g4_89_nat_ab1o3_m2o7_run.json`,
`g4_89_nat_t5_run.json` and of any `g4_89_nat_K<k>_..`, and `merge` ADDED all files it was given. With the present content of `parts/`
(`out/tc.log`) a `--collect` or `--run` for the tag `g4_89_nat` would have merged three different runs, and for `g4_61_nat` two, and
would then have overwritten `g4_split_g4_<p>_nat.json` with the sum and appended FAIL lines to the log of record. It would not have
passed unnoticed (the totals add, so GAP reads 2/180000 or 3/180000; the coverage check passes), but the stored table would have been
replaced. The results of record were first collected before the other part files existed ("collect: 1 part files"), so they were never
affected. Present version: `re.compile("^" + re.escape(tag) + r"_(run|task\d+of\d+)\.json$")` plus an assertion on the tag stored in
each part file; with this pattern no tag picks up a foreign file (`rv_tc.py` (3)).

**B2 (BUG, driver; fixed; matters for the cluster mode).** `merge` of the 08:44 version did not check that the part files cover
disjoint sets of inner graphs (`covered` is a set): a `<tag>_run.json` next to `<tag>_task*of*.json`, or task files of two different
`--ntasks`, would have been double counted (GAP fails, REGULARITY still passes). Present version: an inner graph that occurs in several
part files must have identical pieces and is counted once.

**W1 (weakness; fixed).** The restriction `popcount(S) < K` needs `min valO >= 1`. The 08:51 version did not assert it, and
`frozen_ok` does not imply it (it compares `E_R` with `Eh_0`, not with the table `E_M` in use: mutant m01 had `frozen_ok = True` and
`min valO = 0`; the variant `nat_t0` would have been truncated silently). Present version: assertion in `Split.__init__` and a
separate line in `checks`. In all runs of record `valO = 1` for all 45 keys.

**W2 (weakness; fixed).** GAP / VERTEX tests of the 08:44 version did not require six polar coefficients from `tc_expansion3`
(`not any(T[1:npol])` passes on a shorter list). Present version: `len(T) == npol`. In the logs of record six coefficients are printed.

**W3 (reporting; open).** "REGULARITY" and "VERTEX GRAPH carries the polar part" are the same information given GAP (now said in the
text of the check). "All 378 outer graphs with edges are O(u^6)" should be quoted with the table of section 3.6: 229 trivially, 81 term
by term, 68 by cancellation.

**W4 (coverage of the regression; open).** The genus-3 regression cannot see `tau_4`, `tau_5` or the `zhat`-orders 7..9 (three
surviving mutants). `rv_compose.py` (items 1-6 of section 5; 12 s) is a test that does; it or an equivalent could be added to the
checks of the driver. The controls `nat_t4` and a `zhat >= 7` control are the genus-4 counterparts (done here as m04, m07).

Remarks on the cluster mode (`--task I --ntasks N`), read only, not run: the slices are `r % N == I` over the cost order, so the array
index must be 0-based (`0..N-1`); each task writes its own part file atomically (`.tmp` + `os.replace`) and resumes from it; all tasks
append to the same log file (interleaved lines, harmless). Cost for planning: one full genus-4 `nat` run is about 1450 CPU-seconds on
this machine when it is otherwise idle (298 s wall with 5 processes) and 3250 CPU-seconds in the logs of record (13 processes on a
loaded machine); the slowest single inner graph is below 80 s; the `id` run is 65 CPU-seconds. Nothing at genus 4 comes near the
30-minute limit, so no cluster submission is needed for this genus.

## 7. What is not established by this review

* Everything is a finite computation: genus 4, `zhat <= 9`, `u <= 5`, modulo 2^89-1 and 2^61-1 (frame (0, 0)); the second frame only
  through the stored tables (`rv_tc.py`) and `rv_compose.py` at 2^89-1. The genus-4 mutant runs of this review are at 2^89-1 only; the
  genus-2/3 mutation regression is at both primes.
* The identity of section 2 is ARGUED, not machine-checked; its hypothesis H2 (what "the vertex tensor of `M.omega^`" means) is the
  convention of the old code. Step 4 (legs, leaves) is supported numerically only through the genus-2/3 regression and through items
  3, 4 of section 5 (which check the algebra of Step 0 on the actual tables, not the graph combinatorics of Step 4).
* The shared inputs (rows, conifold frame and `t_c`, symbol matrices, KW numbers, graph generator) were not rebuilt independently here;
  `t_c(u)` was not compared with the literature.
* The selection rule behind the 81 identically vanishing graphs was observed (and shown not to be a dead code path), not derived.
* 26 mutants is a sample. A defect that changes nothing at genus 2, 3 and leaves the genus-4 total, the composition identities of
  section 5 and the unit tests of section 3 intact is not excluded.
* The controls with a shifted `tau_k` (project: `nat_t5`; here: m04, m05) show that the genus-4 regularity is sensitive to each of
  `tau_4`, `tau_5`; whether regularity DETERMINES them (no other value works) was not examined.
