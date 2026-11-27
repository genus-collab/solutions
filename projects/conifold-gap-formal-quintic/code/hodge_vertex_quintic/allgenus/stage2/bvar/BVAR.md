# B-universality: Conjecture N for self-dual cyclic one-pole systems with `B = prod b_a` varied (Stage 2, falsifier (c))

Date 2026-11-14. Folder `allgenus/stage2/bvar/` (new files only; project and Stage-1 modules imported read-only). AI-run and AI-checked
only; no human has read it. Labels: PROVED (argument written here), VERIFIED (exact computation; range, primes stated), CONJECTURAL,
SKETCH, REFUTED. Finite checks are evidence, never proof. Conventions as in `../../SYNTHESIS_STAGE01.md` s.3 (`H = d/dv`, "(G)",
"quintic Theorem R", `K_vc`, `m5`); hatted theory `u = s^5`, `zhat = z/s^4`; a genus-g term is `u^(2-2g) G(u)` and a graph is regular
iff `[u^0..u^(2g-3)] G = 0`. Outer-graph indices = positions in `fq3_graphs.legged_graphs(g, 0, 3g-3)` (as in `univ/`).

## 0. Verdict

**Falsifier (c) PASSED (not refuted).** For five self-dual cyclic one-pole systems with `B` perturbed (three of them also with `f`, `g`,
frame or `lambda` perturbed, one with the residue of `B` changed), with `tau`, `v(u)`, `M = M_0(tau zhat)`, the Theorem-U gauge, its
canonical normalisation and the rows all RECOMPUTED from the system:
* genus 2 (6 edge graphs) and genus 3 (41 edge graphs): **every edge graph regular, in every system, at 2^89-1 and 2^61-1**;
* genus 4, outer graphs with at most two edges (10 edge graphs): **regular** for S1, S2, S3, S5; genus 4 **full**: **all 378 edge
  graphs regular** for S2 (VERIFIED, 2^89-1; s.5);
* the vertex graph is `c_g (T'')^(2-2g) + 0 + .. + 0` in every case, `T'' = -[z^-1 e_0] P Y^C_model = c0 T`, `T = Int_0^u f`,
  `c0 = P_11(0,0) = 1/lambda`, `c_2, c_3, c_4 = 1/1200, 1/25200, 1/180000`;
* the failing controls fail as predicted: stored normalisation (genus 4, the six graphs 1, 2, 3, 4, 7, 8), two non-annihilating
  dressings, and "M not recomputed" (stale `tau` in `E_M`).

This is the B-universality predicted by quintic Theorem R (E6 of the synthesis): the reduction never uses `B` except through `v(u)` and
`tau`, which enter `M` and the gauge, and here they are recomputed. VERIFIED in the stated range; all-genus status unchanged
(CONJECTURAL, conditional on the inputs listed in `../../MECHANISM.md` s.1).

## 1. What changes when B changes (PROVED, elementary)

System `th Y_a = b_a Y_{a+1}`, `th = z d/du`, `b = (f, B/(f^2 g^2), f, g, g)` (self-dual shape: `b_{a-1} = b_{3-a}`),
`f = t_c' F(u)`, `g = -(varpi_0/(1-u)) Gm(u)` in a conifold frame, `B = -m^20 (Pn/Pd)(u) / (u (1-u)^5)`, `Pn(0) = Pd(0) = 1`.
The model `z d/dv Y_b = bmod_b(v) Y_{b+1}`, `bmod = (1, -1/v, 1, -1, -1)`, pulled back by `v(u)` has `prod bt = -v'^5/v`. A
Theorem-U gauge between the two needs `prod b = prod bt` (the determinant of the u-equation of the gauge recursion), i.e.
`v'^5 = -B v`. Writing `beta = -u B` and `beta5 = beta^(1/5) = m^4 (1-u)^-1 (Pn/Pd)^(1/5)`:
`v^(4/5) = (4/5) Int_0^u s^(-1/5) beta5(s) ds`, so with `tau := (u/v)^(4/5)`:
    `1/tau = Sum_n 4 beta5_n/(5n+4) u^n`,   `v = u tau^(-5/4)`,   `tau(0) = m^-4`,   `v'(0) = m^5`.
(Quintic: `m = 1`, `Pn = Pd = 1`, `1/tau = 2F1(1,4/5;9/5;u)`.) The model's natural variable is `z v^(-4/5) = tau zhat`, so
`M = M_0(tau zhat)` with the SAME frozen symbol `M_0` (a model datum) and the NEW `tau`; the outer edges are `E_out = E_R - E_M`,
`E_M(a,b) = tau^(a+b+1) Eh_0(a,b)`. The hatted gauge is `Ghat_(k,b) = tau^((1-c_b)/4) u^mu(b,k) p_(k,b)`; for `m != 1` the prefactor
`tau(0)^((1-c_b)/4) = m^-(1-c_b)` is rational. So `B` enters the graph sum exactly through `v(u)` (the gauge), `tau` (rows and `E_M`).

**Order of the pole (PROVED; why no pole-order perturbation is run).** If `v` is a formal power series with `v(0) = 0` of order `n`,
then `ord_0 B = ord(v'^5/v) = 5(n-1) - n = 4n - 5`: `n = 1` gives a simple pole, `n >= 2` gives a ZERO of order `>= 3` (no pole; `v` is
a branched cover, the natural variable is `z/u^(4n/5)`, not `zhat = z/u^(4/5)`, so the hatted engine and Conjecture N as formulated do
not apply). Conversely, for `ord_0 B = -k`, `v'^5 = -Bv` forces `v ~ u^((5-k)/4)`: for `k = 2, 3, 4` (and every `k != 1 mod 4`) `v` is
not a power series in `u` (no single-valued pull-back, no Theorem-U gauge); `k = 5` gives `v^(4/5) ~ log u`; `k = 9, 13, ..` give a
pole of `v` at `u = 0`, i.e. the point `v = infinity` of the model, which is not its conifold point `v = 0`. Hence a self-dual cyclic
system to which the model pulls back regularly at a conifold point has `ord_0 B = -1` exactly: **changing the order of the pole is
excluded**, and the only remaining "pole" freedom is the residue `-m^20` (system S5).

## 2. Systems (`bv_sys.SYSTEMS`)

| name | F | Gm | frame | B multiplied by | m | lambda |
|---|---|---|---|---|---|---|
| R0 | 1 | 1 | (0,0) | 1 (regression: the formal quintic) | 1 | 1 |
| S1 | 1 | 1 | (0,0) | `1 + u/3 - 2u^2/7` | 1 | 1 |
| S2 | `1 + u/3 - 2u^2/7 + 5u^3/11` | `1 - 3u/5 + u^2/2 + 2u^3/3` | (0,0) | `1 - 5u/4 + 7u^2/3` | 1 | 1 |
| S3 | `3/2 + u/3` | `1 - 3u/5 + u^2/2` | (1/3,-2/7) | `(1-u)^2/(1+2u)^2` (pole at `u = 1` moved) | 1 | 3/2 |
| S4 | `1 - 7u/4 + 3u^2` | `1 + 11u/6` | (0,0) | `(1 + 3u + 5u^2 - 4u^3/3)/(1 + u/2)` | 1 | 1 |
| S5 | 1 | 1 | (0,0) | `2^20 (1 + u/3 - 2u^2/7)` (residue changed) | 2 | 1/1024 |

`tau` and `v` (from `prep_*.log`): S1 `tau = 1 - (64/135)u - ..`, `v = u + (16/27)u^2 + ..`; S2 `tau = 1 - u/3 - ..`; S3
`tau = 1 + (4/45)u - ..`; S4 `tau = 1 - (2/3)u - ..`; S5 `tau = 1/16 - (4/135)u - ..`, `v = 32u + (512/27)u^2 + ..`
(quintic: `tau = 1 - (4/9)u - ..`).

## 3. Inputs (VERIFIED; `logs/prep_R0_S1_S2_S3_S4_S5.log`, all PASS)

For every system, exactly over Q: every `p_(k,b)`, `k <= 12`, of the Theorem-U gauge is a power series in `u` (working precision
`O(u^24)`); in the canonical normalisation `P_11(0,z) = c0` (constant) and (G) `P Y^C_model(v(u)) = c0 Y^C_system` holds on
`z^-6..z^6` mod `u^6` (13 powers of z) with NO positive power of z; in the stored normalisation (G) holds with `P_11(0,z)` having
`z^5, z^10` terms (annihilation violated); the rows are symplectic mod `(zhat^(3g-2), u^(2g-2))`, `g = 2, 3, 4`, both primes.
Regressions: R0 rows = the project rows (25 x (3g-2) x (2g-2) residues, g = 2, 3, 4, both primes); R0's recomputed `E_M` = the
`fq4_split` base `E_M` (g = 2, 3, 4, both primes); R0 graph sums at genus 2, 3 (both primes): all edge graphs regular, vertex
`c_g t_c^(2-2g)` (rows identical to the project rows, so these are the natural runs; R0 was not re-run at genus 4).

*Symplectic scalar (S5).* Theorem U fixes `P` only up to `c(z^5)`; Prop. S gives `P(-z)^T eta P(z) = C eta`. For `m = 1` the gauge
code's normalisation has `C = 1`; for S5 it has `C = 2^-10` (`p_00 p_03 = p_01 p_02 = p_04^2 = 2^-10` at `u = 0`). The canonical gauge
is symplectic, so `P -> sigma P` with `sigma = 1/p_04(0) = 32` (`bv_run.normalise_sp`; `sigma = 1` for R0-S4). Without it the rows
fail the symplectic check (first run of prep, recorded in s.7). Then `c0 = P_11(0,0) = 1024`.

## 4. Genus 2 and 3 (VERIFIED, 2^89-1 and 2^61-1; window `u^0..u^(2g-1)`, i.e. two regular orders beyond the polar part)

| system | genus 2: edge graphs regular | vertex in `T''` | genus 3: edge graphs regular | vertex in `T''` | `c0` |
|---|---|---|---|---|---|
| R0 | 6/6 | `(1/1200) T''^-2 + 0` | 41/41 | `(1/25200) T''^-4 + 0 + 0 + 0` | 1 |
| S1 | 6/6 | same | 41/41 | same | 1 |
| S2 | 6/6 | same | 41/41 | same | 1 |
| S3 | 6/6 | same | 41/41 | same | 2/3 |
| S4 | 6/6 | same | 41/41 | same | 1 |
| S5 | 6/6 | same | 41/41 | same | 1024 |

Identical at both primes. Graph by graph: `BV_TABLES.md` s.2-3 (every entry `R`). Non-triviality: the rows differ from the quintic
rows in 95/100 (genus 2) and 170/175 (genus 3) residue slots, and in the window `u^0..u^(2g-1)` 6 of 7 (genus 2) and 25 of 42
(genus 3) outer graphs differ from R0 already for S1 (`compared with R0` lines of the run logs); i.e. B moves the regular parts of
the graphs but never creates a polar term.

Controls in the re-expansion (run logs): in the quintic `t_c` the vertex graph is NOT pure for S2-S4 (e.g. S2 genus 2:
`(1/1200) t_c^-2 - (1/3600) t_c^-1`); in the model coordinate `v(u)` it is not pure for any system (S1 genus 2: `- (29/162000) v^-1`).
For S1 `T = t_c` (same `f`), so its `t_c` line is not a control. The law is in the system's own `T'' = c0 Int f`, not in `v`.

## 5. Genus 4 (VERIFIED, 2^89-1; window `u^0..u^5`)

| system | run | edge graphs | vertex graph in `T''` | CPU s |
|---|---|---|---|---|
| S2 | maxout 2 | **10/10 regular** (graphs 1-10) | `(1/180000) T''^-6 + 0 + .. + 0` | 686 |
| S1 | maxout 2 | **10/10 regular** | same (`T'' = t_c`) | 633 |
| S3 | maxout 2 | **10/10 regular** | same (`c0 = 2/3`) | 635 |
| S5 | maxout 2 | **10/10 regular** | same (`c0 = 1024`) | 596 |
| S2 | **full** | **378/378 regular**; rows of graphs 0-10 identical to the maxout-2 run | `(1/180000) T''^-6 + 0 + .. + 0`; the total is the same | 1424 |
| S2 | maxout 2, STORED normalisation (control) | **6 of 10 singular**: graphs 1, 2 from `[u^4]`; 3, 4, 7, 8 at `[u^5]` | `T''^-2`, `T''^-1` terms | 704 |

maxout 2 is exact for the vertex graph and the outer graphs with at most two edges (`univ/UNIVERSALITY.md` s.1); graphs with three or
more edges are then incomplete and not reported. Genus-4 outer graphs 0-10: 0 vertex `(4)`; 1 bridge `(1,3)`; 2 bridge `(2,2)`;
3 loop on `(3)`; 4 `(0,3)`; 5, 6 `(1,1,2)` chains; 7, 8, 9 `(1,2)`; 10 two loops on `(2)`. In the quintic `t_c` the S2 and S3 vertex
graphs are not pure (S2: `(1/180000) t_c^-6 - (1/180000) t_c^-5 + ..`); in `v(u)` none is.

**Normalisation response of S2 (VERIFIED; `bv_norm.py`, `logs/norm_S2_g4_89.log`).** `m5 = [z^5] P^stored_11(0,z)/P_11(0,0) =
319129204271/48787200000` (quintic: `-43397/807975000`). `(stored - canonical)/m5` is non-zero exactly on graphs 0, 1, 2 (from `[u^4]`)
and 3, 4, 7, 8 (at `[u^5]`). Its `[u^4]` column is `-1/8400000, -1/1008000, 1/900000` on graphs 0, 1, 2: the SAME universal frozen
response as for the quintic and the univ system c1 (independent of `B`, as it must be: at `u^0` everything is the frozen model). The
`[u^5]` column (`563/3402000000, 947/272160000, -1181/324000000, 1/5040000, -1/1512000, 73/10800000, -17/2700000` on 0, 1, 2, 3, 4, 7, 8)
depends on the system. The sum over the graphs vanishes at `[u^4]` and `[u^5]` (the total is unchanged). So the stored normalisation
is again exactly the annihilation defect `m5 (z^5 Y^C)_+`, and the canonical normalisation (= (G)) makes the six graphs regular.

## 6. Failing controls at genus 2-3 (VERIFIED, 2^89-1; systems S1 and S2, same pattern)

| control | what is wrong | genus 2 | genus 3 | vertex law | total |
|---|---|---|---|---|---|
| `dress0`: rows `(1 + (1/7) zhat E_(2,1)) Qhat` | annihilation violated at `zhat^1`, `u^0` (flat-regular since `mu(1,1) = 0`; symplectic, `X^2 = 0`) | 6/6 singular at `[u^0]` | 41/41 singular at `[u^0]` | holds | changes (S1 g2: `(7223/2469600) T''^-2 + ..`) |
| `dress2`: `1 + (1/7) u^2 zhat E_(2,1)` | annihilation violated at `u^2` | regular (invisible) | 3 singular: 1, 2 at `[u^2]`, 3 at `[u^3]` | holds | changes at `T''^-2`, `T''^-1` |
| `staleM`: correct rows, `E_M` with the QUINTIC `tau` | `M` not recomputed from `B` | 2 singular: graphs 1, 2 at `[u^1]` | 16 singular: graphs 1-16, first orders 1, 2, 3 | FAILS (S1 g2: `+(7/360000) T''^-1`) | unchanged |

Readings. (i) `dress2` reproduces the univ control B2 (same graphs, orders and leading values `1/252000`, `1/14000`) on B-perturbed
systems. (ii) `staleM` is the control specific to this task: the rows, hence the TOTAL (which depends on `R` only), are right, but
splitting `R = R' M` with the wrong `tau` moves polar terms between the vertex graph and the edge graphs. So the graph-by-graph
statement genuinely depends on `B`, through `tau(u)` in `M`, and holds only with `tau` recomputed from `B`: the B-dependence is not
invisible, it is compensated exactly. (iii) No total-potential check can see `staleM` (cf. `../../MECHANISM.md` s.2).

## 7. Observations

* **B only moves regular parts.** For S1 (same `f`, `g` as the quintic) the polar parts of all graphs coincide with the quintic's
  (vertex `c_g t_c^(2-2g)` in the same `t_c`), while 6/7 (genus 2) and 25/42 (genus 3) graphs differ at `u^(2g-2)`, `u^(2g-1)`.
  In general the polar part of the vertex graph is fixed by `T'' = c0 Int f`; `B` enters `T''` only through the constant `c0`
  (S5: `c0 = 1024 = sigma P^code_11(0,0)`, `sigma = 32`).
* **Residue of B.** Changing the residue (S5, `v'(0) = 32`) needs the constant symplectic rescaling of s.3; with it every edge graph is
  regular and the vertex graph is `c_g (1024 T)^(2-2g)`. Without it the rows are not symplectic (observed in the first prep run; that
  log was overwritten when the rescaling was added and is not kept). Pitfall: S5's constants are powers of 2 and both test primes are
  Mersenne (`2^89 = 1` mod `2^89-1`), so rational reconstruction prints e.g. `16777216/5625` for `c_4 1024^-6` in `t_c`; the identity
  holds in F_p and that is what is checked, but powers of 2 should be avoided in future scale tests.
* **Order of the pole.** Excluded by s.1 (PROVED): no one-pole system with `ord_0 B != -1` admits a regular pull-back of the model at
  its conifold point.
* **What is tested.** The prediction of quintic Theorem R that only flat regularity, symplecticity and annihilation of the left factor
  (with ML on `K_vc` for the model) matter; `B` was the last quintic datum Stage 1 had not varied. A failure would have revealed a hidden
  use of the quintic's `B` (through `tau`, the hatted prefactors or the model coordinate).

## 8. Status summary

| statement | status |
|---|---|
| `ord_0 B = -1` is forced for a regular model pull-back at a conifold point; pole-order perturbations excluded | PROVED (s.1) |
| `tau`, `v(u)`, `M` recomputed from `B` reproduce the quintic data for R0 (rows, `E_M`, graph sums) | VERIFIED (g = 2, 3, 4; two primes) |
| Theorem-U gauges exist for S1-S5 (`k <= 12`); (G) holds canonically with no positive z-power | VERIFIED over Q (`z^-6..z^6`, mod `u^6`) |
| B-universality at genus 2, 3: every edge graph regular, vertex `c_g T''^(2-2g)` (S1-S5) | VERIFIED (2^89-1 and 2^61-1) |
| B-universality at genus 4, `<= 2` edges (S1, S2, S3, S5) and in full (S2) | VERIFIED (2^89-1) |
| stored normalisation, non-annihilating dressings, stale `tau`: singular as predicted | VERIFIED (failing controls) |
| B-universality in all genera | CONJECTURAL; implied by quintic Theorem R + ML on `K_vc` (+ Theorem U for the system) |

## 9. Files, code reuse, accounting

| file | content |
|---|---|
| `bv_sys.py` | systems with `B` varied; `tau`, `v(u)` from `B`; copy of `ur_gauge.gauge` without the quintic-`B` check; general `hatted`, rows; `E_M` with the new `tau` |
| `bv_engine.py` | univ engine `Var` on a base with the new `E_M`; small process pool |
| `bv_run.py` | prep checks and graph-sum runs (`--norm`, `--ctrl`, `--maxout`, `--kext`) |
| `bv_norm.py` | genus-4 normalisation response |
| `bv_tables.py`, `BV_TABLES.md` | run summary and graph-by-graph status tables |
| `logs/` | `prep_*.log`, `run_*.log` (one per run), `norm_S2_g4_89.log`, `g4_batch.sh` |
| `tables/` | `tab_*.json` (all residues mod p, per outer graph) |

Reused read-only by import: `natural_form/theory/ur_gauge.py` (`solve_ode`, `prod5`, `mu`), `univ/un_engine.py` (`Var`, `get_base`,
`get_nat`, `collect`), `univ/un_core.py` (rows/tables, symplectic check, `left_dress`), `univ/un_gauge.py` (`check_G`),
`univ/un_e4_cases.py` (`normalise`, `p11_series`), and the project's `fq_*`, `fq3_*`, `fq4_split` modules. Nothing outside
`stage2/bvar/` was written. CPU: genus 2-3 about 3 min in total; genus 4 4678 s (S2 full 1424; five maxout-2 runs of 596-704 each), about 82 CPU min in all. At most 4 concurrent processes (3 pool workers
plus the idle parent; a first launch with 4 workers was stopped within a minute and restarted with 3).

Reproduce: `py -I -B bv_run.py prep R0,S1,S2,S3,S4,S5 --out DIR`; `py -I -B bv_run.py run S2 3 89 --kext 2 --out DIR`;
`py -I -B bv_run.py run S2 4 89 --maxout 2 --nproc 3 --out DIR`; `py -I -B bv_norm.py DIR S2`; `py -I -B bv_tables.py DIR BV_TABLES.md`.
