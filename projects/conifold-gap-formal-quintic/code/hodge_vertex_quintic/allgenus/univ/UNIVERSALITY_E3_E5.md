# Universality, part 2: E3 (annihilation of the vanishing covector) and E5 (characterisation)

Continuation of `UNIVERSALITY.md` (same folder, conventions, status labels and outer-graph numbering; split to keep pieces
below 250 lines). Sections keep their numbers.

## 3. E3: annihilation of the vanishing covector

### 3.1 First order (VERIFIED, genus 2 and 3, 2^89-1 and 2^61-1, frame (0, 0); `un_e3_lin.py`, `un_e3_lin.log`)

Variables: the Taylor coefficients `p_{k,b,j}` of a flat-regular, Z/5-graded, infinitesimally symplectic left dressing of `R'` visible
in the window (`j + mu(b,k) <= 2g - 3`). Response map = polar coefficients of all edge graphs (27 dual-number runs of `un_e3_resp.py`;
re-run today single-process, identical data); annihilation map = `C_m[u^i] = 0` (`(p Y^C_Q)_+ = 0` in the window). The response
to `p_{k,b}(u) H(k,b)` is `p_{k,b}(u) u^mu(b,k)` times the response to the constant `h(k,b)`: the `eps^1` part of the graph sum is
linear over u-series in the dressing (every ingredient is a polynomial in the row entries, `xhat` a rational function of the
`zhat^0` unit row; PROVED).

| genus | variables | rank Ann | dim ker Ann | rank Resp | dim ker Resp | rank(Resp on ker Ann) | rank(Ann on ker Resp) |
|---|---|---|---|---|---|---|---|
| 2 | 9 | 3 | 6 | 3 | 6 | **0** | **0** (ker Resp = ker Ann) |
| 3 | 30 | 10 | 20 | 8 | 22 | **0** | 2 |

* annihilation + flat regularity + symplecticity => every edge graph regular at first order (genus 2, 3);
* at genus 3 the two extra regular directions are supported on the top-order variables `(k, b, j) = (3,0,1), (3,1,1), (4,0,0), (4,1,0)`
  (hatted order `u^3`) and violate only the top-order conditions `(m, u^i) = (1,3), (2,2), (3,1), (4,0)`; they are singular at genus 4
  (s.4: a window artefact);
* on ker Ann the vertex graph moves exactly as `c_g T^(2-2g)` with `T -> t_c + delta T`, `delta T = -[z^-1 e_0] p Y^C_Q`;
* adding the vertex law to edge-graph regularity: genus 2: equivalent to annihilation; genus 3: dimension 21 against 20, the one extra
  combination being `5 w_1 + 2 w_2` of the two extra directions `w_1, w_2` (each alone violates the vertex law at `T^-1`, s.4).

*Dropping flat regularity* (all hatted variables `x_{k,b,j}`, `j <= 2g - 3`): genus 2: 20 variables (11 flat-irregular), rank(Resp on
ker Ann) = 5; genus 3: 68 (38), rank 23. So annihilation + symplecticity alone is far from sufficient. A few flat-irregular directions
are invisible at first order (genus 2: 3, pivots `(1,0,0), (2,2,0), (2,2,1)`; genus 3: 5, pivots `(1,0,0), (2,2,1), (3,2,2), (5,0,3),
(5,1,3)`; `(1,0,0)` is the string-type Pi direction `zhat (H^2(x)1 + 1(x)H^2)` at `u^0`).

### 3.2 Non-linear dressings at genus 2 and 3 (VERIFIED, 2^89-1; `un_e3_nl_89.log`, `un_e3_c_89.log`)

`D = exp(X)` built mod `(zhat^7, u^4)`, restricted to genus 2. Every `D` is symplectic (checked on the dressed rows). All rows of
the table were reproduced at the second prime 2^61-1 (`un_e3_nl_61.log`, `un_e3_c_61.log`; same graphs, same orders), and so was the
first-order analysis of s.3.1 (`un_e3_resp.py 61`, `un_e3_lin.log`).

| case | construction | `delta_0` (z^0 coeff. of `D Y^C_Q`, / `e_1`) | genus 2 | genus 3 |
|---|---|---|---|---|
| A1 | random flat-regular `X` (seed 20261113), column-1 orbits solved for annihilation | `-1 + 4u + (437/90)u^2 + ..` | regular, vertex `c_2 T_eff^-2` | **all 41 regular**, vertex `c_3 T_eff^-4 + 0` |
| A2 | as A1, seed 4242, with `u`-dependent `k = 0` diagonal terms | `-1 + u/3 - (11/360)u^2 + ..` | regular, `c_2 T_eff^-2` | **all 41 regular**, `c_3 T_eff^-4` |
| A3 | as A1, seed 1, `X` times 23/7 | `-1 + (115/14)u + ..` | regular, `c_2 T_eff^-2` | **all 41 regular**, `c_3 T_eff^-4` |
| A4 | seed 2, with `k = 0` terms, `X` times -17/5 | `-1 - (17/20)u - ..` | regular, `c_2 T_eff^-2` | **all 41 regular**, `c_3 T_eff^-4` |
| B1 | A1's free entries, random column-1 entries | | 1 singular (graph 2 at `[u^1]`) | 14 singular (E = 1, 2, 3 from `[u^1], [u^2], [u^3]`); vertex not `c_3 T_eff^-4` |
| B2 | A1 + `(1/7) u^2 zhat H(x)H` (flat-regular, violates `C_1` at `u^2`) | | regular (invisible) | 3 singular: one-edge graphs 1, 2 at `[u^2]` (`1/252000`, `1/14000`), graph 3 at `[u^3]`; vertex still `c_3 T_eff^-4` |
| C1a | A1 + `(1/3) u^0 zhat h(1,3)` (Pi-type `H^4(x)H^3 + H^3(x)H^4`, irregular by one order) | | 3 singular, from `[u^0]` | 28 singular; vertex not `c_3 T_eff^-4` |
| C1b | A1 + `(3/7) u^1 zhat^3 h(3,2)` (irregular by two orders) | | 1 singular (graph 1 at `[u^1]`) | 8 singular |
| C1c | A1 + `(3/7) u^2 zhat^3 h(3,2) + (2/5) u^1 zhat^2 h(2,2)` (irregular by one order each) | | regular | 4 singular (graphs 2, 4, 8, 12 at `[u^3]`); vertex `c_3 T_eff^-4` |

In all A cases `delta_0` is not constant, so `D` is not a gauge between the quintic system and any cyclic system (for such a gauge the
transported pure-`1/z` solution has constant `z^0` term), and `T_eff = -[z^-1 e_0] D Y^C_Q` differs from `t_c` already at `u^2`
(`T_eff = u + (49/45) u^2 + ..` for A1, `t_c = u + (7/10) u^2 + ..`): the vertex graph follows the transported vanishing coordinate.
C1c is regular at first order (both added directions are among the first-order invisible ones of s.3.1) and singular through the
products with A1's terms: flat irregularity adds up along products (`mu(b, k1) + mu(b + k1, k2) = mu(b, k1 + k2)`), so the
first-order invisibility of an irregular direction is not stable.
(The first non-linear run also contained a random flat-irregular case C1 that cannot be made annihilating inside the genus-3 window,
see the note in `un_e3_nl_89.log`, and a buggy first A3/A4 run with non-symplectic scaling, marked invalid in `un_e3_c_89.log`.)

### 3.3 Consistency of the annihilation system with symplecticity (VERIFIED; PROVED at first order)

`delta_m` for `m = 1 mod 5`, `m` even (first: `m = 6`) has as its column-1 entry the self-dual even entry `(6, 1)`, which vanishes for
symplectic `X`, and all its other entries in the window belong to column-1 orbits already fixed by `delta_7..delta_9`. At first order
`delta_6 = y_1 y_2 (x(9,4) - x(9,0)) = 0` because `(9,4)` and `(9,0)` are symplectic partners (`k` odd). Non-linearly, solving
`delta_m = 0` for `m != 6` gave `delta_6 = 0 mod u^6` for all four seeds tried at genus 4 (scratch probe, and `un_e3_g4.log`).

### 3.4 Case A at genus 4 (`un_e3_g4.py`, `un_e3_g4.log`)

Seed 20261113, `X` with 27 entries (`zhat^1..zhat^6`), flat-regular, symplectic, annihilating mod `u^6` for `m = 1..9`;
`delta_0 = -1 + 4u + (1003/210)u^2 + ..` (not a gauge); `T_eff = u + (89/70) u^2 + ..`.
* maxout 2 (one process, 537 s): the vertex graph and all 10 graphs with `<= 2` edges regular; vertex graph `(1/180000) T_eff^-6 + 0 +
  .. + 0`.
* **full run (4 workers, 413 s): all 378 edge graphs REGULAR**; vertex graph `(1/180000) T_eff^-6 + 0 + .. + 0`.
* **second dressing** (seed 2, with `u`-dependent `k = 0` diagonal terms; 29 entries `zhat^0..zhat^6`; `delta_0 = -1 - (4/7)u + ..`,
  `T_eff = u - (137/40) u^2 + ..`), **full run (433 s): all 378 edge graphs REGULAR**; vertex graph `(1/180000) T_eff^-6 + 0 + .. + 0`.

### 3.5 Dressings of the PURE MODEL (VERIFIED, genus 2 and 3, 2^89-1 and 2^61-1; `un_e3_model.py`, `un_e3_model_<p>.log`)

To see Conjecture A away from the quintic, the quintic gauge was dropped altogether: rows `Qhat = diag(tau^((1-c_b)/4)) Rhat(tau zhat)`
(= ur_gauge's hatted form of `P = 1`: the model `x^4 + v/x` at `v(u)` in the quintic's hatted frame; then `E_R = E_M`), vanishing
covector = the pure-`1/z` solution of the pulled-back model system (`y_1 = -v`, `y_2 = v^2/2`; checked), dressings as in s.3.2.

| case | `D` | genus 2 | genus 3 |
|---|---|---|---|
| M0 | 1 (control) | regular (trivially), vertex `c_2 v^-2` | regular, vertex `c_3 v^-4` (`T_eff = v`) |
| MA1 | random annihilating flat-regular symplectic (seed 7); `delta_0 = -1 - u/3 - ..` | **all 6 regular**, `c_2 T_eff^-2` | **all 41 regular**, `c_3 T_eff^-4` |
| MA2 | the same with `k = 0` terms (seed 8); `delta_0 = -1 + (7/3)u + ..` | **all 6 regular**, `c_2 T_eff^-2` | **all 41 regular**, `c_3 T_eff^-4` |
| MB | MA1 + `(1/7) u^2 zhat H(x)H` (annihilation violated) | regular (invisible) | **3 singular**: graphs 1, 2 at `[u^2]`, graph 3 at `[u^3]` (as B2) |

So the principle is a property of the model theory and of the class of dressings, not of the quintic. *Remark (reduction, PROVED,
elementary).* Around the pure model `R' = 1 + eps X`, so `E_out = O(eps)`: at first order only the vertex graph (one changed leaf) and
the one-edge graphs (one edge carrying the bivector of `X`) move. Hence first-order Conjecture A around the model is a finite statement
per genus about the polar parts (in `u`, through `v(u)`) of the model's one- and two-point descendant correlators `<tau_a(H^i)>_h`,
`<tau_a(H^i) tau_b(H^j)>_h` contracted with flat-regular annihilating bivectors; by homogeneity each correlator is a constant times a
power of `v`. This is the smallest model-only statement on which a proof could be attempted (TR is allowed for the model).

## 4. E5 (characterisation)

E3 answers the E5 question negatively: no property beyond (A) is needed through genus 4. In particular
* no flatness / gauge equation: the A dressings are not gauges; and the gauge equation is not sufficient (E4 b, d);
* no constancy of the transported `z^0` coefficient `delta_0`;
* no restriction to symmetry directions (string, Euler-dilaton): the A dressings are random.

**The genus-3 surplus is a window artefact (VERIFIED; `un_e5.py`, `un_e5.log`).** The two directions `w_1, w_2` of s.3.1 that are
regular at first order at genus 3 but violate annihilation (top order `u^3`) were dressed in as finite dressings `D = exp(X)` (`X` of
valuation 3, so `D = 1 + X` in the genus-3 window):
* genus 3: all 41 edge graphs regular for both (as predicted), but the vertex graph is `c_3 T_eff^-4 + 0 + 0 + (11/136080) T_eff^-1`
  (`w_1`) resp. `(-11/54432) T_eff^-1` (`w_2`): the vertex law fails, and `5 w_1 + 2 w_2` is the combination that keeps it;
* genus 4 (`X` completed by column-1 terms at hatted orders `u^4, u^5` so that only the genus-3 top-order conditions are violated;
  maxout 2): `w_1`: **all 10 edge graphs with at most two edges singular**, the one-edge graphs from `[u^3]`
  (`-1/378000`, `13/9720000`, `-1/283500` on graphs 1, 2, 3), and the vertex graph acquires a `T_eff^-3` term (`19/56700000`); `w_2`:
  **9 of 10 singular**, one-edge graphs from `[u^3]` (`1/151200`, `-541/97200000`, `-139/567000`), vertex `T_eff^-3` term `-19/22680000`.
  The vertex-law combination `5 w_1 + 2 w_2` (the response is linear: `X` has valuation 3) has `[u^3]` responses `0`, `-1/225000`,
  `-4/7875` on graphs 1, 2, 3 and `0` on the vertex `T_eff^-3` term: it keeps the vertex law but is singular on graphs 2, 3.
So what genus 3 leaves free at its top order is forced at genus 4: the first-order data are consistent with "regular <=> annihilating"
once the window is large enough.

What (A) does not cover: the first-order invisible flat-irregular directions (s.3.1); they are not stable under composition (C1c), so
(A) is the natural robust class.

