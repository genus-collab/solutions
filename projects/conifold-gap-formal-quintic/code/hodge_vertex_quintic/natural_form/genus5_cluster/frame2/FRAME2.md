# Genus 5, the second conifold frame: the point-frame run, additivity, and whether the selection rule is frame dependent

Agent report, 2026-11-13.  Scope: the gap left open in `NATURAL_FORM.md` s.6, "**The second frame is single-prime (2^89-1) and
unchecked by additivity**".  The additivity half of that sentence is now closed, the frame's total is now two-prime, and a
structural question nobody had looked at (is the non-zero pattern of the pieces frame dependent?) is answered.
Everything below is a finite exact computation mod the two Mersenne primes `p89 = 2^89-1`, `p61 = 2^61-1`: **evidence, not proof**.

Code written for this report: `f2_lib.py`, `f2_add.py`, `f2_total.py`, `f2_sel.py`, `f2_seldetail.py` (this directory).
Logs and data: `out/`.  Nothing in `../../../../cluster_genus5_quintic/` was modified; the copied-back cluster evidence was read
only.  `g5_checks.checks` was **not** used: every sum in this report is computed by `f2_*.py`.

---

## 1. What was missing

The cluster ran the natural factorisation in the second conifold frame `(alpha, beta) = (1/3, -2/7)` at `p89`
(tag `g5_89_nat_ab1o3_m2o7`, full `--maxout`, 52.6 CPU h, 100 array tasks).  Its collect log says, verbatim:

        (no complete g5_89_id_ab1o3_m2o7 run: the additivity check is skipped)

so the one check that ties the 340369 stored pieces of that run to anything outside the run itself never ran in that frame.
`g5_checks.checks` skips it silently and still reports `failed = []`, which is why the note flagged it.

The missing ingredient is the point-frame (`--variant id`) run in the same frame: for each inner graph `t` it computes the single
number family `graph_id(sp, t)` = the point-frame term of `t`, with `M = 1` (`E_M = 0`), which is what the pieces of `t` must add
up to.

## 2. The two point-frame runs made here

        cd cluster_genus5_quintic
        py -I -B genus5_run.py --genus 5 --prime 89 --variant id --ab=1/3,-2/7 --run --nproc 13
        py -I -B genus5_run.py --genus 5 --prime 61 --variant id --ab=1/3,-2/7 --run --nproc 13

Both completed with `ALL CHECKS PASSED`; the full console logs are `out/run_g5_89_id_ab1o3_m2o7.log` and
`out/run_g5_61_id_ab1o3_m2o7.log`, the part files and result JSONs are in `out/parts/` and `out/results/`.

| run | wall | CPU (13 workers) | all 4555 graphs | GAP check |
|---|---|---|---|---|
| `g5_89_id_ab1o3_m2o7` | 650 s | 7193 s = 2.0 h | yes | PASS |
| `g5_61_id_ab1o3_m2o7` | 1003 s | 10159 s = 2.8 h | yes | PASS |

(the `p61` wall time is longer only because the audit scripts of s.3 and s.4 were running beside it on the same machine.)

**A concurrency event, for provenance.** `out/logs/g5_89_id_ab1o3_m2o7/local.log` shows that a *second* run of the same tag,
with `--nproc 20`, was started by someone else on this machine about five minutes into the run made here, finished 14 s after
it, and - `write_part` being an atomic `os.replace` - **overwrote `parts/g5_89_id_ab1o3_m2o7/local.json` with its own copy**.
So the `p89` part file now in `out/parts/` is the 20-worker run's, and the CPU seconds it records (10059 s = 2.8 h) are that
run's, not the 2.0 h of the 13-worker run.  Both runs then ran their own collect, and `out/logs/.../collect.log` holds the two
side by side: **all eight residues and the whole `t_c` polar part are identical**, so the duplication cost nothing and
incidentally confirms that the result does not depend on how the 4555 graphs are spread over workers.  Everything audited
below was audited against the file as it now stands, and it passes both the additivity check and the comparison with the
cluster's stored total.  The `p61` tag was run once only.

For comparison the frame-`(0, 0)` point-frame runs on this machine were 651 s wall / 1.8 CPU h, so the second frame costs about
10% more per prime (2.0 CPU h at `p89`; the 2.8 CPU h figures are inflated by other work sharing the 28 cores).  The package's
own estimate of 1.3 CPU h per prime was optimistic by 40-50%, in both frames.

**The GAP check holds in the second frame too.** At `p89` the point-frame total is

        residues  ['1/660000', '-1/225000', '323/207900000', '~66841/15592500000', '[...]', '~...', '~...', '[...]']
        in t_c    (1/660000) t_c^-8 + (0) t_c^-7 + ... + (0) t_c^-1,      c_5 = |B_10|/(10 . 8 . 5^4) = 1/660000.

Those eight residues are **bit-for-bit the eight residues of `info['total']` of the cluster's own
`results/g5_89_nat_ab1o3_m2o7.json`** (checked in `f2_add.py`, printed in `out/additivity_frame2_p89.log`): the natural run and
the point-frame run of the second frame produce the same total, as they must, the first as a sum of 340369 pieces and the second
as a sum of 4555 point-frame terms.

## 3. Per-inner-graph additivity in the second frame: 4555 / 4555

`f2_add.py` reads the cluster's `parts/g5_89_nat_ab1o3_m2o7/task{0..99}of100.json` (copied back, read-only) and the new
`parts/g5_89_id_ab1o3_m2o7/local.json`, merges each by its own `read_parts_dir` (asserting that a graph occurring in two part
files carries identical content, and that every part file carries the expected tag), and for every inner graph `t = 0..4554`
checks

        sum over the outer graphs q of piece(t -> q)   ==   point-frame term of t        (8 residues mod p89).

        natural run  g5_89_nat_ab1o3_m2o7: 100 part files, 4555 inner graphs, CPU 52.6 h
        point frame  g5_89_id_ab1o3_m2o7:    1 part file,  4555 inner graphs, CPU  2.0 h
        coverage: natural 4555/4555 (missing []), point frame 4555/4555 (missing [])
        point-frame part files storing a piece q != t: 0 []
        stored (t -> q) pieces in the natural run: 340369
        inner graphs whose point-frame term is identically zero: 0
        ADDITIVITY: 4555 / 4555 identities verified; 0 failures
        RESULT: PASS

**4555 identities, 0 failures, 36440 residues compared.**  The last line of the quoted block matters: *no* inner graph has a
vanishing point-frame term, so not one of the 4555 identities is the trivial `0 = 0`.  (`critic/CRITIC.md` had established the
same for the frame-`(0, 0)` runs.)

### 3a. The same code re-closes the frame-`(0, 0)` additivity

The note records the frame-`(0, 0)` additivity as "trusted from the cluster log" (point (ii) of the open list).  With the parts
copied back it costs nothing to redo it here with the same independent summation code:

| frame | prime | identities | failures | pieces summed |
|---|---|---|---|---|
| `(1/3, -2/7)` | 2^89-1 | 4555 | 0 | 340369 |
| `(0, 0)` | 2^89-1 | 4555 | 0 | 340363 |
| `(0, 0)` | 2^61-1 | 4555 | 0 | 340363 |

(`out/additivity_frame00_p89.log`, `out/additivity_frame00_p61.log`.)  **13665 additivity identities in all, 0 failures**, each
against a point-frame term that is non-zero, and in all three cases the summed total also equals `info['total']` of the
corresponding `results/*.json`.  So the cluster's two "ADDITIVITY PASS" log lines are now confirmed from the stored pieces
rather than trusted, and the third case, which never ran, passes as well.

---

## 4. Is the selection rule frame dependent?  Rows: no.  Individual pieces: yes, from `u^3` up.

Nobody had compared the two frames piece by piece.  `f2_sel.py` compares the two **full natural runs at the same prime**
`p89` -- the cluster's `g5_89_nat` and `g5_89_nat_ab1o3_m2o7` -- at three levels, and `f2_seldetail.py` localises the
differences.  All inner and outer indices refer to the same list of 4555 graphs, so the two runs are directly comparable.

| level | frame `(0, 0)` | frame `(1/3, -2/7)` | agreement |
|---|---|---|---|
| (a) stored pieces `(t -> q)` (somewhere non-zero) | 340363 | 340369 | 30 only in `(0, 0)`, 36 only in the second frame |
| (b) `u`-support of the 340333 common pieces | | | 124 differ |
| (c) outer graphs with a stored row | 628 | 628 | same 628 indices |
| (c) outer graphs whose row is non-zero | 1 (the vertex) | 1 (the vertex) | same |
| (c) outer graphs whose row cancels to 0 | 627 | 627 | **same 627 indices**, same `u`-support |

**(c) The row-level selection rule is exactly frame independent.**  The same 628 outer graphs receive a contribution (so the
same 1036 of the 1663 graphs with `E <= 7` receive none, with the same breakdown 7, 71, 288, 670 at `E = 4..7`), the same 627 of
them cancel to zero, the vertex graph is the only non-zero row in both, and no row differs in which of its eight coefficients
vanish.  The census of `NATURAL_FORM.md` s.6 ("1036 by selection, 627 by cancellation") is therefore a statement about the
theory, not about the frame `(0, 0)`.

**(a), (b) The piece-level pattern is *not* frame independent -- 190 of the 340429 pieces occurring in either frame differ.**
Before reading anything into that, the obvious alternative explanation (an accidental zero mod `p`) has to be excluded.
Running the same comparison between the **two primes in the same frame** (`g5_61_nat` vs `g5_89_nat`) gives

        (a) PIECES: 340363 in both; only in one: 0; (b) 0 of 340363 common pieces differ in u-support;
            u-support histogram identical: True; (c) rows identical.

(`out/selection_twoprime_frame00.log`.)  So in this computation there is **not a single** accidental mod-`p` vanishing among
340363 pieces x 8 coefficients at either prime, and the 190 frame differences are genuine.

Where they sit (`out/selection_detail.log`):

* the **66 pieces present in one frame only** (30 + 36) have `u`-support exactly `{u^7}` in the frame where they exist -- i.e.
  they live entirely in the **top coefficient of the `K = 8` window**, the one order the genus-5 run cannot see beyond;
* the lowest `u`-order at which any support differs at all is `u^3`; the distribution of that lowest order over the 190 pieces
  is `u^3`: 16, `u^4`: 13, `u^5`: 52, `u^6`: 43, `u^7`: 66.  **No piece differs at `u^0`, `u^1` or `u^2`**;
* all 190 have `E(t) >= 7` and `E(q) >= 3` -- the heavy inner graphs feeding the deep outer graphs.

Read together with (c): a change of conifold frame reshuffles individual pieces at high `u`-order (which is unsurprising, since
`M`, `E_out` and the vertex tensors all change with `(alpha, beta)`), but it moves **nothing** at low order and it does not
change which outer graph receives anything, nor which cancels.  The selection rule is a property of `E_out`, as the note says,
and it survives the frame change intact; the floor-order structure that the controls exploit is frame independent as far as this
comparison can see.

---

## 5. The second frame's total over the rationals, from the two primes

`f2_total.py` sums the 4555 point-frame terms of each prime with its own code (asserting completeness and that the `id` runs
store only the diagonal piece `q = t`), CRT-combines the eight coefficients against `p89 . p61 ~ 2^150`, reconstructs them by
continued fractions, re-reduces each reconstructed rational modulo **both** primes to confirm it, and then performs the change
of variable `u -> t_c` **over the rationals** (`fq_fp.QRing`, no prime at all; `fq_conifold.Conifold` is built over `Q`
internally, so `fq_hat.tc_series` / `fq3_msum.tc_expansion3` run exactly).

**Second frame, `(alpha, beta) = (1/3, -2/7)`:**

        total = u^-8 ( 1/660000 - u/225000 + 323 u^2/207900000 + 66841 u^3/15592500000 - 32135797 u^4/46777500000000
                       - 917102699 u^5/350831250000000 - 3515716092179 u^6/3683728125000000000
                       + 1892778904559 u^7/5023265625000000000 )

        t_c = u + (11/30) u^2 + (751/1575) u^3 + (454/945) u^4 + ...

        total = (1/660000) t_c^-8 + 0 t_c^-7 + 0 t_c^-6 + 0 t_c^-5 + 0 t_c^-4 + 0 t_c^-3 + 0 t_c^-2 + 0 t_c^-1   over Q.

The reconstruction is **certified, not probabilistic**: the largest height among the eight coefficients is `2^62.1`, and the
smallest margin `log2(p89 . p61 / (2 |n| d))` is **45.6 bits** (at `u^6`), so each rational is the unique one of that height
matching the residues.  The heights of the second frame are larger than those of frame `(0, 0)` (`2^62.1` against `2^54.9`,
margins 45.6 against 61.2 bits).

**The `u`-series totals of the two frames differ at every order `u^1 .. u^7`, and `t_c(u)` differs too** (`u + (7/10)u^2 + ...`
against `u + (11/30)u^2 + ...`), **yet the eight polar coefficients in `t_c` agree exactly.**  That is the content of the frame
check: seven vanishing conditions, computed from different numbers in different coordinates, come out zero both times.
As a by-product the same script reproduces the frame-`(0, 0)` total of `NATURAL_FORM.md` s.6 coefficient for coefficient
(`out/total_frame2.log`), from the copied-back `g5_{89,61}_id` parts and with an independent CRT, reconstruction and rational
`t_c` -- one more independent confirmation of that series, after `recheck/RECHECK.md`.

### 5a. Two more single-prime `~` values in the logs are spurious

`NATURAL_FORM.md` s.6 warns that the `~`-marked single-prime reconstructions in the November 11 point-frame logs are wrong.
The same holds in the second frame, and it is worth recording which, because the cluster's `g5_89_nat_ab1o3_m2o7` collect log
is the only record of that run's total:

| order | printed at `p89` alone | true value (two primes) |
|---|---|---|
| `u^3` | `~66841/15592500000` | `66841/15592500000` -- **confirmed** |
| `u^5` | `~-1764310128775/5569957551811` | `-917102699/350831250000000` -- **spurious** |
| `u^6` | `~1267427754587/8985989270387` | `-3515716092179/3683728125000000000` -- **spurious** |

(`u^4` and `u^7` printed as bare residues at `p89`, correctly declining to reconstruct.)  The `~` marker does its job; nothing
downstream used those two values.

---

## 6. What the second frame now has, and what it still lacks

**Has.**
* A complete point-frame run at **both** primes (`g5_89_id_ab1o3_m2o7`, `g5_61_id_ab1o3_m2o7`; 4555/4555 inner graphs each,
  2.0 and 2.8 CPU h, `ALL CHECKS PASSED` including GAP), hence a **two-prime, certified-rational total** whose `t_c` polar part is
  `(1/660000) t_c^-8 + 0 + ... + 0` although its `u`-series and its `t_c(u)` differ from frame `(0, 0)` at every order.
* **Per-inner-graph additivity, 4555/4555, 0 failures** against the cluster's 340369 stored pieces: the check `g5_checks` had
  skipped. This ties every one of the 4555 inner graphs of the regularity run to an independently computed point-frame term;
  it does **not** constrain how an inner graph's weight is distributed over the outer graphs (the split that regularity is
  about), which is what `critic/CRITIC.md`'s 168475 re-derived `(t, S)` assignments address.
* The row-level selection rule shown to be frame independent (s.4).

**Lacks.**
* **Its regularity run is still single-prime** (`p89` only).  A second-prime natural run in this frame is the expensive item
  (~52 CPU h, i.e. a cluster job, not a local one); nothing cheap closes it.  What the two-prime point-frame runs now give is
  the *total* at two primes in this frame, not the 628 rows.
* **Independence is unchanged.**  Both new runs use the package's own `fq4_split.Split` / `g5_engine`, so they are not an
  independent re-implementation; what they are is an independent *route* (`M = 1`, no edge splitting, `Split.graph` with the
  single subset `S = all`) to the number the 340369 pieces must add up to.  A structural zero inside `graph_dp` would still be
  common-mode, but a structural zero in the piece bookkeeping would not: it would break additivity, and it does not.
* Nothing here tests `tau_6`, or the `zhat^10..zhat^12` orders of `M_0`, which remain the missing controls of
  `NATURAL_FORM.md` s.6.

### 6a. One observation about the package directory

While this work ran, `cluster_genus5_quintic/logs/` acquired the directories `g5_89_nat_ab1o3_m2o7`, `g5_89_M0_out2` and
`g5_89_nat_t7_out2`, which were not there when the session started.  They were not created here (the only runs made here are the
two `id_ab1o3_m2o7` ones, whose parts, logs and results were moved into `out/` afterwards); they are the Google Drive sync of
the copied-back cluster evidence still landing.  The package's `results/` and `parts/g5_{89,61}_id/` are untouched.

## 7. Files

* `f2_lib.py` -- part-file merging (own `read_parts_dir`, with the duplicate-content and tag assertions), CRT, continued-fraction
  rational reconstruction.
* `f2_add.py` -- the additivity check; `--natTag/--idTag` let it run on any frame/prime pair.
* `f2_total.py` -- two-prime total, certified reconstruction, and the `u -> t_c` change of variable over `Q`.
* `f2_sel.py`, `f2_seldetail.py` -- the frame comparison of s.4 and its localisation.
* `out/run_g5_{89,61}_id_ab1o3_m2o7.log` -- the two run consoles; `out/{parts,logs,results}/` -- their part files, package logs
  and result JSONs (5.1 MB in all).
* `out/additivity_frame2_p89.{log,json}`, `out/additivity_frame00_p{89,61}.{log,json}`, `out/total_frame2.{log,json}`,
  `out/selection_rule.{log,json}`, `out/selection_twoprime_frame00.{log,json}`, `out/selection_detail.{log,json}`.

Reproduce (from this directory, `EV` = the extracted cluster evidence, `PK` = `../../../../cluster_genus5_quintic`):

        py -I -B f2_add.py   "$EV/parts/g5_89_nat_ab1o3_m2o7" out/parts/g5_89_id_ab1o3_m2o7 --prime 89
        py -I -B f2_total.py --dir89 out/parts/g5_89_id_ab1o3_m2o7 --dir61 out/parts/g5_61_id_ab1o3_m2o7 \
                             --dir89b "$EV/parts/g5_89_id" --dir61b "$EV/parts/g5_61_id"
        py -I -B f2_sel.py   "$EV/parts/g5_89_nat" "$EV/parts/g5_89_nat_ab1o3_m2o7"
