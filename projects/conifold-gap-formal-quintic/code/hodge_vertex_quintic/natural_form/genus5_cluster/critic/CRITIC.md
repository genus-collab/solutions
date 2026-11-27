# Adversarial review of the genus-5 regularity run (formal quintic, natural factorisation)

Date 2026-11-13. Written and run by an AI agent in the role of a hostile referee; no human has checked it.
Material reviewed: `cluster_genus5_quintic/` (code and the six copied `results/*.json`), `hodge_vertex_quintic/NATURAL_FORM.md`
s.5-6, `natural_form/genus4_referee/`, `natural_form/genus5_package_review/`. The cluster's `logs/` and `parts/` were **not**
copied back, so everything below that concerns the run itself rests on the six JSON files plus the package's own code.
Nothing in `cluster_genus5_quintic/` or in any existing module or note was modified. The critic's own code is in this folder
(`c1`..`c13`), is new, and writes nothing outside the scratchpad.

## 0. Verdict

The claim survives the attacks on **execution** and on **arithmetic accident**: the copied files, checked with the critic's
own arithmetic, rule out a short or partial run, and the graph list (sound *and* complete), the automorphism counts, the
piece-to-outer-graph assignment, the Witten-Kontsevich numbers and the arithmetic at a third prime all hold up under
independent recomputation.

The claim does **not** survive the attack on **what is being claimed**. "All 4554 outer graphs with edges are regular at
genus 5" is, as a statement about the run, inflated by a factor of about ten:

* **2891 of the 4554 (63 %) were never computed.** The engine generates outer-edge sets of size at most `K-1 = 7`, so no
  outer graph with 8 or more edges receives a single piece. Their regularity is a one-line a-priori consequence of
  `E_out = O(u)`, not a computation. (`c5_contract_audit.py`.)
* Of the remaining 1663, **1036 receive only pieces that vanish identically** and have no stored row at all; their
  "regularity" is read off the *absence* of data, which the stored file cannot distinguish from "nothing was produced".
* **627 graphs vanish by actual cancellation.** That is the content. In coefficients: 3142 polar coefficients are at stake
  (not 4554 x 8 = 36432), of which 1655 are cancellations and 1487 are identically-vanishing pieces.

And the headline that reads best -- the predicted gap `c_5 = 1/660000` -- is **produced by the run regardless of whether the
factorisation is right**: both deliberately-wrong controls report `c_vertex_equals_pattern = True`. The Bernoulli
confirmation was already complete in the point-frame (`id`) run done on the workstation before the cluster job; the cluster
job adds nothing to it.

Four things must be fixed or copied back before this is written up as "VERIFIED at genus 5"; they are listed in s.8.

## 1. Attack (a): could the run have silently done less than it claims?

**Answer: no, not on the evidence of the copied files -- and this is the strongest part of the result.** Two independent
arguments:

1. `failed: []` in each natural run means the collect step's first check,
   `"all 4555 inner graphs are covered by the part files"`, passed. That check is
   `set(range(nin)) - set(pieces)`, and a task that died, was never submitted, or wrote nothing leaves its inner graphs out
   of `pieces` entirely. `genus5_run.py` writes `results/<tag>.json` only after `C.checks` returns, so a JSON with an empty
   `failed` list cannot come from a run with missing tasks. The genus-4 end-to-end test in
   `genus5_package_review/REVIEW_PKG.md` s.3.4 cut a task at 32 of 76 graphs and the coverage check duly failed.
2. Independently of that code: the critic re-summed the stored `table` of every run and compared with the point-frame
   totals (`c1_json_audit.py`). `sum(table) == info['total']` in all six runs; and
   `total(g5_89_nat) == total(g5_89_id)`, `total(g5_61_nat) == total(g5_61_id)` exactly. Every one of the 4555 rows of the
   `id` table is non-zero, so dropping any inner graph from the natural run would move its total by that graph's
   point-frame term. A false pass needs the omitted graphs' point-frame rows to sum to zero in all 8 coefficients at both
   primes.

So "half the tasks missing" is excluded. Three smaller execution holes remain, and they are real:

* **The additivity check leaves no trace.** `g5_checks.checks` runs the per-inner-graph additivity check only when
  `len(id_pieces) == len(sp.inner_graphs)`; otherwise it logs `"(no complete g5_89_id run: the additivity check is
  skipped)"` and *does not fail*. `info` has no entry for it, and `failed == []` is identical whether the check passed or
  was skipped. The collect log, the only record, was not copied back. For the second frame the check is skipped by design
  (no `id` run in that frame), so for `g5_89_nat_ab1o3_m2o7` the strongest per-graph consistency check certainly did not
  run. For the two frame-(0,0) runs it probably ran -- the `id` parts are shipped in the package -- but nothing in the
  copied files says so.
* **No provenance.** The results JSON holds only `tag`, `table`, `info`, `failed`, `outer`. No host, no date, no CPU hours,
  no part-file names, no code digest, no Python version. An auditor with `results/` alone cannot tell whether the file came
  from a 100-task SLURM array, from `--run` on a laptop, or from a different version of the code. The CPU total (39-46 h
  per prime) exists only in the collect log.
* **The cluster environment is uncertified.** `V_validate.sh` re-runs the entire genus-2/3/4 regression with the cluster's
  own interpreter (eight genus-4 tables of record, the old vertex-tensor tables at genus 2-3) and is the package's main
  guard against a platform or version problem. Its log, `logs/validate.log`, was not copied back either. Only the
  *workstation* copy (`logs/workstation/validate.log`, `G5_VALIDATE ALL CHECKS PASSED`) is available.

Two further hazards found in the code, neither affecting these results: `read_parts` merges **any** `task<I>of<N>.json`
regardless of `N`, and also the `local.json` that the package shipped (`parts/g5_89_nat/local.json`, 2 inner graphs computed
on the workstation), so the certified result silently mixes provenances; and `pieces[t]` may legitimately be `{}`, which
counts as covered. Both are caught in aggregate by the total-additivity argument above, not by the coverage check.

## 2. Attack (b): is 150 bits of arithmetic enough?

The naive form of this attack fails. The standard argument -- "a non-zero rational of the size seen elsewhere in the
computation cannot vanish mod a 61-bit prime, let alone mod both" -- is exactly the argument that worked at genus 4. At
genus 5 **it no longer works**, and that is a new fact nobody has written down.

The critic reconstructed the per-graph point-frame coefficients from three primes, adding 2^521-1 (also `= 1 mod 5`, the
prime an earlier referee used) with a fresh run of the `id` variant (`c8_thirdprime.py`, `c11_heights_by_E.py`):

* On a sample of 31 inner graphs, all 248 coefficients reconstruct from the three primes, every one is non-zero, and every
  reconstruction reduces correctly to each of the three stored residues: **0 inconsistencies**. That is a genuine
  third-prime check of the row machinery and of the shipped `id` table, at a modulus of a completely different size.
* **Only 11266 of the 36440 per-graph point-frame coefficients reconstruct at all from `p89 * p61`** (`c7_height.py`);
  the other 69 % have height above 2^70. The heights measured with the third prime are:

  | E of the graph | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 10 | 12 |
  |---|---|---|---|---|---|---|---|---|---|---|---|
  | max height (bits) | 136 | 130 | 126 | 123 | 122 | 117 | 107 | 107 | 104 | 95 | 89 |

  `p89 * p61` is `2^150`. So the margin between the observed size of these numbers and the threshold at which a non-zero
  rational *could* be zero mod both primes is **14 bits for the vertex graph and 20-50 bits for the others**. At genus 4
  the margin was of the order of a hundred bits.

What survives: the *heuristic* probability that a given non-zero coefficient is divisible by `p89 * p61` is about `2^-150`,
so across the 3142 polar coefficients at stake the chance of a false zero is about `2^-138`. That is negligible. But the
correct statement of the evidence is now **probabilistic, not deductive**, and the note should say so. A third prime (the
`id` run at 2^521-1 cost 6 s for 31 graphs here, so a few CPU hours for all 4555; the natural run would be dearer) would
restore a deductive margin and would also make the per-graph rationals visible, which two primes cannot do.

Two further points on this attack:

* **The two primes are not independent as a test of the code.** All rational inputs -- the Witten-Kontsevich numbers, the
  `tau` coefficients `4/(5n+4)`, the automorphism orders, the graph list, the frame `(alpha, beta)` -- are the same objects
  reduced modulo different primes. An error in any of them reproduces identically at every prime. Primes test the modular
  reduction only. Both primes also run the same `g5_engine` code, so every structural zero (a piece that is never produced,
  a subset that is never generated, an empty decoration table) is prime-independent. All of s.1 and s.3 below is therefore
  outside the reach of "two primes agree".
* **No division can silently wrap.** `FpRing.inv` raises `ZeroDivisionError` on 0 and is never called on anything larger
  than an automorphism order or a small-prime denominator (`|Aut| <= 2^12 * 12!`, denominators are products of 2, 3, 5, 7,
  11), none of which is divisible by `2^61-1` or `2^89-1`. `Split.__init__` asserts the symbol valuations, the hatted-row
  valuations and `min_prec >= K`. There is no silent wrap-around path.
* **`ratrec` honesty.** `FpRing.ratrec` caps the height at `min(isqrt(p/2)/100, 10^10)`, which is about `1.07e7` at
  2^61-1. Every rational printed in a single-prime log is therefore either small or flagged with `~`, and the
  genus-5 spot check already found one wrong `~` value in `logs/g5_89_id/collect.log` (the `u^6` coefficient). The
  rationals of record must come from the CRT of both primes, as they do in `NATURAL_FORM.md`.

## 3. Attack (c): is the right quantity being tested?

Yes as to the window, and yes as to the frame, but the *statement* is weaker than the genus-2/3 one.

* **The window is right.** `K = 2g-2 = 8`; a genus-5 term is `u^-8 G(u)`, and the regularity test is literally
  "the whole stored row of 8 coefficients of `G` is zero mod p" (`g5_checks.checks`: `tab.get(i)[:npol]` with
  `npol == K == 8`). The prefactor `u^-8` is uniform over graphs -- it comes from the `xhat` powers, whose total exponent
  is `sum_v (2g_v - 2 + n_v) = 2g-2` for every leg-free stable graph -- and the `id` table confirms this empirically:
  all 4555 point-frame rows are non-zero at `u^0`, i.e. in the point frame every graph has a full `u^-8` pole and nothing
  is regular. So the test has content and the normalisation is not graph-dependent.
* **The frame is right and was varied.** Two conifold frames were run; their `u`-series differ while their `t_c` polar
  parts agree (checked by the critic in `c1_json_audit.py`), which is the correct frame-independence statement.
* **But the quantity is the split sum, not a graph sum built from genuine vertex tensors.** The run never constructs `R'`.
  It computes the point-frame sum over inner graphs with each edge marked inner (`E_M`) or outer
  (`E_out = E_R - E_M = M^-1 E_R' M^-T`), regroups by contracting the inner edges, and calls the result the outer graph
  sum `R'.(M.omega^)`. That regrouping identity is an *argument* (`genus4_referee/review/REVIEW.md`), audited by an AI
  agent, and it has been checked against genuine vertex tensors **only at genus 2 and 3** (`data/ref_g23.json` is the old
  vertex-tensor code; `data/ref_g4.json` is `fq4_split`'s own output, so the genus-4 regression is self-referential on this
  point). The genus-5 run tests the same statement as genus 4 -- same code, same quantity -- and a strictly weaker one than
  genus 2-3, where the same numbers also came out of real vertex tensors.
* **A per-graph claim that no check constrains.** The total and the per-inner-graph additivity both sum over the outer-edge
  sets `S`, so neither constrains *which* outer graph a piece is added to. That assignment -- `Split.contract` composed
  with `Split.outer_pos` -- is the whole content of the per-graph claim, and the only genuinely independent validation of
  it is at genus 2-4. The critic closed this hole by hand (s.7, item 3): 168475 `(inner graph, S)` pairs recomputed with
  the critic's own contraction and canonical form, **0 disagreements**.

## 4. Attack (d): the predicted gap is not evidence for the factorisation

`c_5 = 1/660000 = |B_10|/(10 . 8 . 5^4)` was predicted before the run, and the run confirms it. **The confirmation is
worth nothing as a test of the natural form**, for a reason visible in the copied files themselves:

| run | verdict | `c_vertex_equals_pattern` | vertex `t_c` polar part non-zero at |
|---|---|---|---|
| `g5_89_nat` | passes | **True** | `t_c^-8` only |
| `g5_89_nat_ab1o3_m2o7` | passes | **True** | `t_c^-8` only |
| `g5_89_M0_out2` (bare symbol, must fail) | FAILS | **True** | all eight, `t_c^-8` .. `t_c^-1` |
| `g5_89_nat_t7_out2` (`tau_7` shifted, must fail) | FAILS | **True** | `t_c^-8` and `t_c^-1` |

Both deliberately broken controls reproduce the Bernoulli coefficient exactly. The reason is structural: `E_out` is
`O(u)` (valuation exactly 1, all 78 keys, checked in `c9_eout.py`), so at order `u^0` no outer edge can appear and the
`u^0` coefficient of the vertex graph equals the `u^0` coefficient of the *whole point-frame total*, for any `M` whose
`[u^0]` is the frozen symbol. `c_5` is therefore a restatement of the `id` run -- which was done on the workstation
*before* the cluster job and is already quoted in `NATURAL_FORM.md` s.6 -- and is insensitive to `tau`, to the frame, and
to everything the cluster run computed.

The same applies to two of the three checks the natural runs report:

* **`GAP: total = c t_c^-8 + 0 + ... + 0` is not a test of `M` either.** The `id` runs pass it (`c_total_equals_pattern =
  True`, `failed = []`). The total is `M`-independent: it is the point-frame total, as the critic verified
  (`total(nat) == total(id)` at both primes).
* **"The vertex graph carries the whole polar part" and "every edge graph is regular" are one statement, not two.** Given
  the gap, each follows from the other. The genus-4 referee said this already (`SYNTHESIS_G4.md` via
  `genus4_referee/review/REVIEW.md`); `info['vertex'] == info['total_tc']` in all three natural runs is an identity, not a
  check.

So: of the five checks reported by a natural run, exactly one (REGULARITY) is about the natural factorisation. What would
have to be wrong for the leading coefficient to come out right and regularity still to be false? Nothing exotic: any error
at orders `u^1..u^7` of `tau` or of `M_0` leaves `c_5` untouched by construction. The `nat_t7` control is the proof of
that in the copied files.

**And the discrimination controls are far too weak for genus 5.** Both were run with `--maxout 2`, which by construction
tests only the vertex graph and the 1- and 2-edge outer graphs: **14 of the 1663 graphs that carry content**, the ones
where the natural run has its *easiest* cancellations. The genus-5-specific content is elsewhere:

* the 3-to-7-edge sector: 613 of the 627 cancelling graphs lie there (34 + 93 + 168 + 204 + 114), against the 14 at `E <= 2` that the controls reach;
* `tau_6` and `tau_7`, which are invisible at genus 4 (`K = 6` there). **Only `tau_7` was control-shifted. There is no
  control for `tau_6` at all** -- nothing shows that the run would have failed had `tau_6` been wrong;
* the orders `zhat^10, zhat^11, zhat^12` of `M_0` (`KZ = 3g-3 = 12`), which enter for the first time at genus 5. The
  genus-4 referee already flagged that `tau_4`, `tau_5` and `zhat^7..zhat^9` were tested by genus-4 regularity alone; the
  genus-5 analogue is tested by genus-5 regularity alone, with no control, and the `maxout 2` controls cannot reach it.

A cheap fix: run `--variant nat_t6 --maxout 2` (same cost as the existing controls, about 10-15 CPU h), and one control
whose perturbation is invisible at `E <= 2` but visible at `E >= 4`, at full `maxout`.

## 5. Attack (e): the statement the files support, and the statement s.6 will be read as

### 5.1 The strongest statement the copied files support

> With `M = M_0(tau(u) zhat)`, `1/tau = 2F1(1, 4/5; 9/5; u)`, and `R' = Psi(u) G`, the genus-5 outer graph sum of the
> formal quintic produced by edge splitting (`fq4_split` / `g5_engine`) was computed over all 4555 leg-free stable graphs
> of genus 5, mod `2^89-1` in the conifold frames `(0,0)` and `(1/3,-2/7)` and mod `2^61-1` in the frame `(0,0)`. In each
> run every outer graph with at least one edge has `[u^0..u^7] G = 0`, and the vertex graph is exactly
> `(1/660000) t_c^-8 + 0 t_c^-7 + ... + 0 t_c^-1`, with `1/660000 = |B_10|/(10 . 8 . 5^4)`.
>
> Of the 4554 outer graphs with edges, **2891 (those with 8 or more edges) are regular a priori**, because every entry of
> `E_out = E_R - E_M` has `u`-valuation exactly 1, so a term with `s` outer edges is `O(u^s)`; the engine does not
> generate them. Of the remaining 1663, **1036 receive only pieces that vanish identically** and **627 vanish by
> cancellation** -- 3 with one edge, 11 with two, 34 with three, 93 of 100 with four, 168 of 239 with five, 204 of 492
> with six, 114 of 784 with seven. In coefficients: **3142 polar coefficients are tested, 1655 of them cancellations.**
> The leading coefficient `c_5` and the gap of the total are properties of the point frame, hold in every variant
> including the two broken controls, and were established by the `id` runs before the cluster job.
>
> Coverage of all 4555 inner graphs is certified by the collect step and, independently, by the natural total agreeing
> with the point-frame total at both primes. Everything is evidence at two primes of 150 bits in total, with a margin of
> 14-50 bits over the observed height of the numbers involved; it is not proof. All of it was produced and checked by AI
> agents; no human has refereed any of it.

The genus-4 wording in `NATURAL_FORM.md` s.5 already does this honestly ("the 229 outer graphs with at least six edges are
`O(u^6)` for free ... 81 ... receive only pieces that vanish identically ... 68 vanish by cancellation"). The critic
reproduced that exact genus-4 accounting (229 / 81 / 68) with independent code before applying the same method at genus 5,
so the genus-5 numbers above are the direct analogue. **s.6 must be written in the s.5 style, not in the s.0 style.**

### 5.2 What a reader of the current s.6 will take away, and the gaps

`NATURAL_FORM.md` s.6 still says "Regularity of the 4554 outer graphs with edges (natural `M`): **NOT YET RUN**". The
danger is the sentence that replaces it. A reader of s.0/s.1 item 7 updated in the obvious way ("every stable graph with
edges is regular at genus 2 (6 graphs), 3 (41), 4 (378) and 5 (4554); VERIFIED, two frames, both primes") would be
over-reading on every one of these points:

| what the reader takes away | what is true |
|---|---|
| 4554 graphs were checked regular | 1663 were computed; 627 of them cancel; 2891 are free |
| two frames at both primes | frame `(0,0)` at both primes; frame `(1/3,-2/7)` at `2^89-1` only |
| the second frame is as well checked as the first | in the second frame the per-inner-graph additivity check **cannot** run (no `id` run in that frame) and is silently skipped |
| the run confirms the Bernoulli pattern `c_5` | `c_5` comes from the point frame, is reproduced by both failing controls, and was already in s.6 before the run |
| the controls show the test discriminates | they discriminate on 14 of 1663 graphs (`--maxout 2`) and on `tau_7`; nothing tests `tau_6` or the `zhat^10..zhat^12` orders that are new at genus 5 |
| the per-graph values are known rationals | 69 % of the per-graph coefficients do not even reconstruct from the two primes |
| an independent implementation agrees | the independent genus-4 code was run at genus 5 on **17 of 4555** graphs (`genus5_package_review/spot`), none of the heavy ones, at one prime |
| the computation is re-auditable | `logs/` and `parts/` were not copied back; no single non-trivial outer graph can be re-checked without redoing most of the run |
| Conjecture N is verified | Conjecture N is a statement for all `g >= 2`; this is one more finite genus. Theorem U, the all-order identification of `M_0`, the regrouping identity and the purity hypothesis (H) are unchanged, and none has a human referee |

**Evidence versus proof.** Nothing in this run is proof of anything. What is *proved* (modulo the stated hypotheses) is
unchanged from s.1-s.3. What is *argued* and AI-audited: the regrouping of the split sum into the outer graph sum,
Theorem U, the symbol identification. What is *verified*: the finite table above. What is *AI-checked only*: all of it,
including this review.

**Not run.** The second frame at `2^61-1`; any third prime on the natural run; an order-by-order sensitivity scan (which
coefficient of `tau` or which order of `zhat` each outer-graph sector is actually sensitive to); a `tau_6` control; a
control that bites only at `E >= 4`; the independent genus-4 implementation on the full genus-5 set; left placement; the
non-formal quintic.

## 6. What the critic recomputed independently (all of it new code, in this folder)

| # | check | script | result |
|---|---|---|---|
| 1 | The six JSONs are internally consistent: `sum(table) == info['total']`; `total(nat) == total(id)` at each prime; the CRT of the two primes gives the stated `u`-series; `[u^0] == 1/660000` mod both primes | `c1_json_audit.py` | all pass |
| 2 | The shipped genus-5 graph list: 4555 entries, every one connected, stable, leg-free, of genus 5; **|Aut| recomputed by backtracking (not by refinement + individualisation) -- 0 disagreements**; **pairwise non-isomorphic under the critic's own brute-force canonical key -- 4555 distinct**; closed under one-edge contraction. Same checks at genus 4: 379, all pass | `c2_graphlib.py`, `c4_list_audit.py` | all pass |
| 3 | The piece-to-outer-graph assignment `outer_pos[contract(G, S)]`: **168475 `(inner graph, S)` pairs** (all 1663 graphs with `E <= 7`, plus 40 random graphs with `E >= 8`) recomputed with the critic's own contraction and canonical form; and the outer graph's edge count must equal `popcount(S)` | `c5_contract_audit.py` | 0 disagreements |
| 4 | Per-inner-graph **ADDITIVITY** (the check the copied files cannot report): the split pieces of the 250 cheapest inner graphs sum to their shipped point-frame terms. 19869 non-zero pieces seen; these touch 521 of the 627 cancelling outer graphs | `c6_resample.py` | 0 failures |
| 5 | `graph_dp` (the fast cluster engine) against `fq4_split.Split.graph` (the slow subset-by-subset path) on 80 genus-5 inner graphs, entry by entry | `c6_resample.py` | 0 disagreements |
| 6 | Valuation structure: every piece going to an outer graph with `s` edges is `O(u^s)` -- what makes the `|S| <= 7` cap sound and the 2891 many-edge graphs free | `c6_resample.py` | 0 violations |
| 7 | `E_out` is not degenerate: 1950 of 1950 colour entries non-zero, valuation exactly 1 for all 78 `(a,b)` keys | `c9_eout.py` | pass |
| 8 | **A third prime.** `2^521-1` on 31 inner graphs in the point frame; the three primes reconstruct one rational per coefficient that reduces correctly to each stored residue | `c8_thirdprime.py` | 248/248, 0 inconsistencies |
| 9 | Heights of the per-graph coefficients (the margin of the two-prime zero test) | `c7_height.py`, `c11_heights_by_E.py` | only 11266 of 36440 reconstruct from `p89 p61`; max height 2^136 against the threshold 2^150 |
| 10 | The **Witten-Kontsevich numbers** `fq_hodge.wk` -- rational inputs that no number of primes tests: `<tau_{3g-2}>_g = 1/(24^g g!)` for `g <= 5`; the genus-0 `n`-point formula for `n <= 7`; the dilaton equation on 2340 correlators with `g <= 5`; published low-genus values | `c10_wk.py` | all pass |
| 11 | Independent enumeration of the leg-free stable graphs from scratch (own route: vertex count, genus vector, degree sequence, multigraph fill, own canonical form), and the same count a second time by orbit counting with no canonical form at all (`sum |Aut_vertices| / V!`): genus 2, 3, 4 give **7, 42, 379** with the by-edge distributions `1,2,2,2` / `1,2,5,9,12,8,5` / `1,3,7,21,43,75,89,81,42,17` | `c3_enumerate.py`, `c12_count_orbits.py` | pass at `g <= 4` (genus 5 by this route is too slow in pure Python; superseded by item 12) |
| 12 | **COMPLETENESS of the 4555-graph list, independently proved.** Every stable graph with `E >= 1` edges contracts to a stable graph with `E-1` edges, so it is enough that each level is the set of one-edge expansions (vertex splittings and loop insertions) of the previous one. The critic generated the expansions with its own code and its own canonical form and compared level by level: `3, 11, 34, 100, 239, 492, 784, 1002, 926, 632, 260, 71` generated, **0 classes missing from the shipped list and 0 in the list that are not expansions**, from the base case `E = 0` (the single genus-5 vertex). Same at genus 2, 3, 4 | `c13_expand.py` | pass; the count 4555 no longer depends on the package's generator or on a literature value |

Two of these (items 2, 3) and item 10 close holes that nothing in the package covers at genus 5: `g5_validate.py` checks
`|Aut|` by brute force only for the 428 graphs of genus 2-4, and the assignment map and the Witten-Kontsevich table are
validated against genuinely independent code only at genus 2-4 as well.

## 7. Problems, in the order a referee would raise them

**P1 (wording; the one that matters).** "All 4554 outer graphs with edges are regular at genus 5" is true but
ten times weaker than it sounds: 2891 are regular by a valuation argument the engine never evaluates, 1036 have no stored
piece at all, 627 cancel. Write s.6 with the 2891 / 1036 / 627 split and the 3142 / 1655 coefficient count, exactly as s.5
does for genus 4.

**P2 (wording).** Do not present `c_5 = 1/660000` as a confirmation of the natural form. It is a property of the point
frame; it is reproduced by both failing controls; it was already established by the `id` runs quoted in the current s.6.
Likewise, do not report "the vertex graph carries the whole polar part" and "every edge graph is regular" as two
independent findings, and do not report the GAP check as evidence for `M` -- the `id` runs pass it.

**P3 (missing control).** `tau_6` enters for the first time at genus 5 and no control tests it; the orders
`zhat^10..zhat^12` of `M_0` likewise. The two controls that were run are `--maxout 2`, i.e. they discriminate on 14 of the
1663 content graphs. Run `--variant nat_t6 --maxout 2` (about 10-15 CPU h) and at least one control whose perturbation is
invisible at `E <= 2` and visible at `E >= 4`, at full `maxout`.

**P4 (missing run).** The second conifold frame was run at `2^89-1` only, and in that frame the per-inner-graph
additivity check cannot run at all (there is no `id` run in that frame, so `g5_checks` skips it silently). Either run
`--prime 61 --variant nat --ab 1/3,-2/7` and `--prime 89 --variant id --ab 1/3,-2/7`, or state in s.6 that the second
frame is single-prime and has no additivity check.

**P5 (auditability).** `logs/` and `parts/` were not copied back. Consequences: (i) it is not recorded whether the
additivity check ran or was skipped -- `failed: []` is identical either way and `info` has no entry for it; (ii) the CPU
time, host, interpreter and part-file inventory of the run are lost; (iii) the cluster's `logs/validate.log`, the only
evidence that the cluster environment reproduced the genus-2/3/4 tables, is lost; (iv) no single non-trivial outer graph
can be re-checked by anybody without redoing a large fraction of the 40-46 CPU h run, because each one collects pieces
from hundreds or thousands of inner graphs. Copy back `logs/` (a few MB) and `parts/` (tens of MB per run); they are what
makes the result re-auditable and what an independent re-implementation would be compared against piece by piece.

**P6 (code hardening, for the next run).** `g5_checks.checks` should put every check in `info` with one of
`pass / fail / skipped`, so a skipped check is visible in the result file; and the check
`"every entry of E_out = E_R - E_M is O(u)"` is `min(valO) >= 1`, which would also **pass if `E_out` were identically
zero** -- in which case every outer graph with edges would be trivially regular and the vertex graph would carry
everything, i.e. exactly the claimed output. Assert `min(valO) == 1` and that `E_out` has a non-zero entry (it does:
1950 of 1950, `c9_eout.py`). Add the host, the interpreter version, the CPU seconds, the part-file names and a digest of
the `.py` files to `info`.

**P7 (arithmetic, new at genus 5).** The argument "a non-zero rational of this size cannot vanish mod both primes" is no
longer available: the per-graph coefficients reach height `2^136` against the threshold `p89 p61 = 2^150`, a margin of
14-50 bits, and 69 % of them do not reconstruct from the two primes at all. The evidence is now probabilistic
(about `2^-138` over the 3142 coefficients at stake, still negligible) and s.6 should say so rather than repeat the
genus-4 deductive argument. A third prime -- `2^521-1` works with `FpRing` and the `id` variant costs a few CPU hours --
would restore a deductive margin and make the per-graph rationals visible.

**P8 (closed by this review).** The genus-5 graph list is sound and complete: 4555 pairwise non-isomorphic valid genus-5
stable graphs with independently recomputed automorphism orders (`c4_list_audit.py`), and the list is exactly the
one-edge-expansion closure of the single genus-5 vertex, checked level by level with the critic's own expansion code and
canonical form (`c13_expand.py`). So the count no longer rests on the package's generator, on `fq3_graphs.degenerations`,
or on a literature value. `g5_validate.py` should still be extended to check genus-5 `|Aut|` against a brute-force count
(it checks only the 428 graphs of genus 2-4).

**P9 (scope).** The genus-5 run tests the split-sum quantity. Its identification with the outer graph sum
`R'.(M.omega^)` is an AI-audited argument checked against genuine vertex tensors only at genus 2 and 3; at genus 4 the
"tables of record" are `fq4_split`'s own output. So genus 5 inherits exactly the genus-4 scope, not the genus-2/3 scope,
and `NATURAL_FORM.md` s.3 should keep saying so.

**P10 (independence).** No independent re-implementation covers genus 5 beyond the 17 graphs of
`genus5_package_review/spot` (one prime, none of the heavy graphs). Both primes, both frames and all six runs share
`g5_engine`, so every structural zero is common-mode. If `parts/` comes back, running `genus4_referee/indep/` against a
few hundred of the stored per-inner-graph pieces -- which it can already do, as the spot check showed -- is the cheapest
genuine independence available.

## 8. What would make this referee believe it

1. Copy back `logs/` and `parts/` of all five cluster runs (P5). Without `parts/` the per-graph claim is unfalsifiable by
   anyone who does not rerun the whole job.
2. Rewrite s.6 with the 2891 / 1036 / 627 accounting and the coefficient count, and demote `c_5` to a point-frame fact
   (P1, P2).
3. Run the `nat_t6` control and one control that bites only at `E >= 4` (P3); run the second frame at the second prime and
   its `id` partner (P4).
4. State the arithmetic evidence as probabilistic, or add a third prime (P7).
5. Have an independent implementation recompute a few hundred of the stored per-inner-graph pieces once `parts/` is back
   (P10). (The 4555-graph list no longer needs this: see s.6 item 12.)

None of these is expensive compared with the 80-90 CPU hours already spent, and the first two cost nothing.

## 9. Files of this review

All new, all in `hodge_vertex_quintic/natural_form/genus5_cluster/critic/`; run with `py -I -B <script>` from that folder.
Nothing outside the scratchpad is written and no existing module, note or part file was touched.

| script | what it does |
|---|---|
| `c1_json_audit.py` | audit of the six copied `results/*.json` alone: table sums, natural vs point-frame totals, CRT of the two primes, the by-edge accounting of the regularity claim, the controls' `c_vertex_equals_pattern` |
| `c2_graphlib.py` | the critic's own stable-graph library (multiplicity matrices, backtracking automorphism count, brute-force canonical key) |
| `c3_enumerate.py` | from-scratch enumeration of leg-free stable graphs (genus 2-4) |
| `c4_list_audit.py` | the shipped graph list: validity, `|Aut|`, pairwise non-isomorphism, closure under contraction |
| `c5_contract_audit.py` | the `(inner graph, S) -> outer graph` assignment, 168475 pairs |
| `c6_resample.py` | re-runs the engine: per-inner-graph additivity on 250 graphs, `graph_dp` vs `Split.graph` on 80, piece valuations |
| `c7_height.py` | how many per-graph coefficients reconstruct from the two primes |
| `c8_thirdprime.py` | `2^521-1` on a sample: three-prime consistency and reconstruction |
| `c9_eout.py` | `E_out` is not degenerate (the valuation check would pass if it were zero) |
| `c10_wk.py` | independent identities for the Witten-Kontsevich numbers (prime-independent inputs) |
| `c11_heights_by_E.py` | height of the per-graph coefficients by edge count, i.e. the margin of the two-prime zero test |
| `c12_count_orbits.py` | the graph count again, by orbit counting, with no canonical form |
| `c13_expand.py` | completeness of the shipped list by one-edge expansions, level by level |

Note for whoever edits the note: `NATURAL_FORM.md` s.6 and s.1 item 7 are stale (they still say "NOT YET RUN" and
"CONJECTURAL for g >= 5"). s.5.1 above is the sentence that the copied files support.
