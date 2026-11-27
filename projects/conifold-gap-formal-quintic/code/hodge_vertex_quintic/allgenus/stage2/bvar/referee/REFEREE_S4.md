# Referee report, part 2: further attempts to break falsifier (c) (`allgenus/stage2/bvar/`)

Continuation of `REFEREE.md` (s.4 there points here). Same labels and conventions. AI-written; no human has read it.

**4.1 A module-disjoint engine (`rf_indep2.py`).** My rows were summed on `natural_form/genus4_referee/indep`. That engine does
its own Fourier-space colour sum and its own division by `z + w`. It builds `E_M = (tau/tau0)^(a+b+1) E_M0` itself, from our `tau`
and from `M_0` read off our rows at `u = 0`, and it asserts `[u^0] E_R = E_M0`, the `E_M` division identity, the Z/5 charges and
`E_out = O(u)`. The genus-2/3 graph order differs from `un_engine`, so I compare multisets of per-graph series there; genus 4 is
compared graph by graph.

Results (VERIFIED; every engine assert passed):
* S5 (genus 3, 2^61-1), S6 (genus 2, 2^89-1), S7 (genus 3, 2^61-1): the multisets of per-graph series are IDENTICAL to the
  `un_engine` runs. No edge graph is singular, and the vertex law holds.
* S6, genus 4, maxout 2, 2^61-1: IDENTICAL graph by graph (11 graphs x 6 coefficients); 10/10 regular.
* S2, genus 4, FULL, 2^89-1: **378/378 edge graphs regular**, and the vertex law holds. The table is IDENTICAL to the worker's
  `tab_S2_g4_89_canon_none.json` on all 379 outer graphs x 6 coefficients. This re-derives the worker's full genus-4 claim with
  my rows on an engine that shares no module with theirs.

**4.2 New systems (VERIFIED, both primes at genus 2, 3; genus 4 in `REFEREE.md` s.2.1: S6 maxout 2, S7 in full, 2^61-1).**
* S6: `m = 3` (no powers of 2, which avoids the Mersenne pitfall the worker reports), `f(0) = 2`, frame `(-1/2, 1/5)`,
  `Pd = 1 + u/5`, a cubic `Pn`; `tau(0) = 1/81`, `v'(0) = 243`, `c0 = 59049/2`.
* S7: `f(0) = -5/7 < 0`, `Gm(0) = 4/3`, frame `(0, 3/11)`, `B` times `(1 + 6u^5)/(1-u)` (pole at `u = 1` of order 6).
In both: 6/6 and 41/41 edge graphs regular, and the vertex law `c_g T''^(2-2g) + 0` holds, at 2^89-1 and 2^61-1.

**4.3 The other real branch of the pull-back (S8; `rf_extra.py branch`).** `v'^5 = -B v` is invariant under `v -> -v`. S8 is the
formal quintic itself (B unperturbed) pulled back by `v = -u - (5/9)u^2 - ..` (`m = -1`, real fifth roots in the hatted
prefactor). The system's R-matrix is intrinsic, so nothing may change.

Result (VERIFIED, g = 2, 3, at 2^89-1 and 2^61-1):
* The S8 rows equal the project rows in every residue, and every S8 table equals the R0 table.
* Exactly over Q (`k <= 7`), the canonical gauge is `P_S8 = -P_q diag(1,-1,-1,1,-1)`.
* The real-root signs `(-1)^(1-c_b)` of the hatted prefactor cancel the diagonal signs.
So the branch `v -> -v` is a symmetry of the construction.

Disclosure: my first S8 run gave non-symplectic rows and singular graphs. The cause was a bug in my own sign patch:
`(-1)**(negative int)` is a float in Python, so rounding entered the exact series. After the fix the result is as stated. The bad
run is kept as `x_branch_floatbug.out` in the scratch folder.

**4.4 Further failing controls (2^61-1, where the worker used 2^89-1 only).**

| control (2^61-1) | genus 2 | genus 3 | vertex law | total |
|---|---|---|---|---|
| staleM, S1 (`E_M` with the quintic `tau`) | graphs 1, 2 singular | graphs 1-16 singular | FAILS (`+(7/360000) T''^-1`) | unchanged |
| staleM, S6 (`tau0 = 1/81`) | 6/6 singular | 41/41 singular | FAILS | unchanged |
| staleHat, S1 (only the hatted prefactor from the quintic `tau`) | regular | graphs 1, 2 singular | FAILS (`-(1/81000) T''^-1`) | changes |
| staleHat, S6 | regular | graphs 1, 2, 3 singular | FAILS (leading `1/10800`) | changes |
| dress0, S1 and S6 (`1 + (1/7) zhat E_21`) | 6/6 singular at `u^0` | 41/41 singular at `u^0` | holds | changes (S1 g2: `7223/2469600`, as the worker found at 2^89-1) |
| dress2, S1 and S6 (`1 + (1/7) u^2 zhat E_21`) | regular (invisible) | graphs 1, 2 singular at `u^2` (S1: `1/252000`, `1/14000`, as the worker found), graph 3 at `u^3` | holds | changes |
| staleboth, S1 (S1 gauge; quintic `tau` in the prefactor, in `tau^(n-k)` and in `E_M`; rows stay symplectic) | regular | graphs 1, 2 singular | FAILS (`+(1/16200) T''^-1` at g = 2) | changes |
| nosig, S5 (`P` off the symplectic scale by a constant: `C = 2^20`) | rows not symplectic; 6/6 singular at `u^0` | 41/41 singular | holds (in its own `c0`) | changes |

staleboth matters most here. A CONSISTENTLY wrong `tau` keeps the rows symplectic but still fails, so the pass needs the system's
own `tau`.

**4.5 Non-triviality.** Counted from my tables (2^89-1, window `u^0..u^(2g-1)`):
* genus 2: every outer graph that is non-zero in the window (6 of 7) differs from R0, for each of S1-S7;
* genus 3: 25 of the 26 non-zero graphs differ (the other 16 vanish in the window).
This agrees with the worker's "6 of 7" and "25 of 42". The differences sit only at the regular orders `u^(2g-2)`, `u^(2g-1)`.

