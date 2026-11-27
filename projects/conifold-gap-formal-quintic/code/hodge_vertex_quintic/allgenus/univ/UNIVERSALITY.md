# Universality of the regular-gauge factorisation: which property of R' makes the edge graphs regular (E1-E5)

Date 2026-11-13. Folder `hodge_vertex_quintic/allgenus/univ/` (new files only; project modules imported read-only). Exact computations
modulo p = 2^89-1 on exact rational inputs (gauges over Q); AI-run and AI-checked only, no human has read any of it. Status labels:
PROVED (argument written here), VERIFIED (exact computation; genus, orders, prime, frame stated), CONJECTURAL, REFUTED. Finite
computations are evidence, never proof. Conventions: NATURAL_FORM.md, UNIFORM_REDUCTION.md (Theorem U, (R1)-(R3)); hatted theory,
`u = 1 - 5^5 q`, `zhat = z/s^4`, `c = (0,-1,3,2,1)`, `mu(b,k) = (4k + c_b - c_{b+k})/5`; a genus-g term is `u^(2-2g) G(u)` and a graph is
regular iff `[u^0..u^(2g-3)] G = 0`. Outer-graph indices are positions in `fq3_graphs.legged_graphs(g, 0, 3g-3)` (the order of
`fq4_split.Split.inner_graphs`; listed in the tables below where needed), NOT Mathematica indices.

## 0. Summary

Question: which property of the left factor `R' = Psi G` (flat frame: the regular gauge `P(u, z)`) makes every edge graph of the outer
graph sum `R'.(M.omega^)`, `M = M_0(tau(u) zhat)`, regular? Candidates: (A) flat regularity + symplecticity + the annihilation
`(P Y^C)_+ = 0` of the vanishing covector (the content of (G)); (B) being a Theorem-U gauge between two cyclic systems (flatness),
in canonical normalisation; (C) something specific to the quintic. Results:

1. **Genus 4, canonical vs stored normalisation (E4/U4; VERIFIED, 2^89-1).** Two perturbed cyclic systems with their Theorem-U gauges
   in the CANONICAL normalisation `P_11(0, z) = const` (c1: generic `f`, `g`; c2: `f(0) = 3/2`): **all 378 edge graphs regular**,
   vertex graph = total = `(1/180000) (T/lambda)^-6 + 0 + ... + 0` in the system's own vanishing coordinate `T = Int_0^u f`. The quintic
   and the c1 system in the STORED normalisation (`P_00(0, z) = 1`, i.e. the canonical gauge times `1 + m5 z^5`): **singular, exactly
   as predicted**: the one-edge graphs 1, 2 from `[u^4]`, the one-edge graph 3 and the two-edge graphs 4, 7, 8 at `[u^5]`, the vertex
   graph acquires `t_c^-2`, `t_c^-1` terms; every other graph is unchanged (valuation argument, s.2.3, PROVED). The `[u^4]` responses
   per unit `m5` are universal (identical for the quintic and c1), and the change of the TOTAL vanishes through `u^5`: the normalisation
   only moves polar terms between the vertex graph and six edge graphs. The quintic in canonical normalisation reproduces the
   project rows exactly (25 x 10 x 6 residues) and the natural run (case a; s.2).
2. **Annihilation suffices; flatness is not needed (E3; VERIFIED genus 2, 3 at two primes, genus 4 at 2^89-1).** Random left
   dressings `D = exp(X)` of `R'` that are flat-regular, symplectic and annihilating, `(D Y^C_Q)_+ = 0`, but NOT gauges (the `z^0`
   coefficient `delta_0(u) e_1` of `D Y^C_Q` is not constant, so `D` cannot map the quintic system to any cyclic system): **all edge
   graphs regular at genus 2 and 3** (four dressings A1-A4) **and at genus 4** (two further dressings, all 378 edge graphs, s.3.4),
   and the vertex graph is `c_g T_eff^(2-2g)` + regular with `T_eff = -[z^-1 e_0] D Y^C_Q`. The same holds for annihilating dressings
   of the PURE MODEL, with no quintic input at all (genus 2-3, s.3.5). At first order this holds on the whole space: annihilation +
   flat regularity + symplecticity imply first-order regularity of every edge graph (genus 2: 9 variables, genus 3: 30; complete
   response map).
3. **Both conditions are needed.** Violating annihilation (B1 random; B2 = A1 plus the single flat-regular term `(1/7) u^2 zhat H(x)H`;
   E4 cases b, d: a scalar `z^5` factor, which is still a flat Theorem-U gauge) makes edge graphs singular. Violating flat regularity
   while keeping annihilation (C1a: `u^0 zhat` Pi-type; C1b: `u^1 zhat^3 h(3,2)`; C1c: two terms irregular by one order each) makes edge
   graphs singular at genus 3 (C1c only through non-linear terms: both its directions are invisible at first order). At first order,
   regularity is EQUIVALENT to annihilation on the flat-regular dressings at genus 2; at genus 3 regularity leaves free exactly two
   combinations of the top-order (`u^3`) annihilation conditions, and both are singular at genus 4 (s.4): a window artefact.
4. **Conclusions.** (A) regularity + (G) only: SUFFICIENT at genus 2-3 (non-linear, four random non-gauge dressings) and at genus 4
   (two random non-gauge dressings, s.3.4), also around the pure model (genus 2-3, s.3.5), and necessary in the sense of item 3.
   (B) any Theorem-U gauge: regular iff in the canonical normalisation (genus 4); a gauge in a non-canonical normalisation is still flat
   and singular, so flatness is neither sufficient nor (item 2) necessary. (C) nothing
   quintic-specific is used beyond the factorisation `R = R' M_0(tau zhat)` with `R'` in class (A); for the quintic that class membership
   is Theorem U + (G) + the canonical normalisation of the actual R-matrix. CONJECTURE A (s.5) is the all-genus statement.

## 1. Method

* **Engine.** The outer graph sum is a function of the hatted rows `Qhat` of `R = R' M` only: `E_R` (edge formula of `fq_sser`), leaves
  `T_R(d)_i = (-1)^d Qhat_0^(i)[zhat^(d-1)]`, vertex factor `xhat` from `Qhat_0^(i)[zhat^0]`, `E_out = E_R - E_M`, `M = M_0(tau zhat)`
  fixed.
  `un_engine.Var` (subclass of `fq4_split.Split`, imported read-only) rebuilds `E_R`, `T_R`, `xhat` from any table of rows (optionally
  with a dual number `eps`), allows `E_out` of valuation 0 in `eps`-parts, and has `maxout m` (only outer-edge subsets of size `<= m`:
  exact for the vertex graph and the outer graphs with `<= m` edges). Regression: `un_e1.log` (genus 2, 3; re-run today, PASS) and case a
  below (genus 4).
* **Rows from a cyclic system (E2, E4).** For `th Y_a = b_a Y_{a+1}`, `b = (f, B/(f^2 g^2), f, g, g)`, `B = -1/(u(1-u)^5)`, `f = t_c' F(u)`,
  `g = -(varpi_0/(1-u)) Gm(u)` (self-dual shape, `B` kept so that `v(u)`, `tau` are the quintic's), the Theorem-U gauge `P'` to the model
  is computed over Q by `natural_form/theory/ur_gauge.gauge` with its coefficient functions replaced in memory (`un_gauge.py`), and
  `Qhat' = Ghat'(u, zhat) Rhat(tau zhat)` with the frozen symbol rows `Rhat` of the project.
* **Normalisation (`un_e4_cases.py`).** Theorem U leaves `P -> P c(z^5)`. *Stored* = output of ur_gauge (free constant 0 at every
  resonance: `P_00(0, z) = 1`; for the quintic `P_11(0, z) = 1 - (43397/807975000) z^5 - ...`). *Canonical* = `P -> P c0/P_11(0, z)`,
  `c0 = P_11(0, 0)`: `P_11(0, z) = c0` constant (`c0 = 1/lambda`, `lambda = f(0)`), equivalently `(P Y^C_model)_+ = 0`. In hatted variables
  `z^5 E_{b+5,b} = u^4 zhat^5 E_{b+5,b}`: genus 2-3 cannot see the normalisation, genus 4 sees exactly its `z^5` coefficient, and the two
  normalisations differ by the scalar right factor `1 + m5 u^4 zhat^5` (checked on the rows, s.2.1).
* **Dressings (E3, `un_ann.py`).** `D = exp(X)`, `X = {(a, b, k): u-series}` hatted, Z/5-graded, infinitesimally symplectic
  (`x_{k,b} = -(-1)^k x_{k,b*}`, `b* = 3 - b - k`), flat-regular iff `x_{k,b}` is divisible by `u^mu(b,k)`; rows -> `D` rows. Hatted
  vanishing covector `Yhat^C = s^2 Sum_n zhat^-n yhat_n e_{1-n}`, `yhat_n = y_n/u^(n - [n/5])` (units; checked), annihilation
  `delta_m(D) = Sum_n D[(m+1, 1-n, m+n)] yhat_n = 0` for `m >= 1`; the column-1 orbits are solved for it by a u-adic Newton iteration.

## 2. E4/U4: genus 4 with and without the canonical normalisation

### 2.1 Inputs (VERIFIED, `e4.log`, `py -I -B un_e4.py prep a,b,c1,c2,d,c3,c4,e`)

For every case: every `p_{k,b}` (`k <= 12`) a power series in `u` (Theorem U); (G) `P Y^C_model(v(u)) = P_11(0, z) Y^C_system` on
`z^-6..z^6` mod `u^6` exactly over Q; rows symplectic mod `(zhat^10, u^6)`.

| case | system | normalisation | `P_11(0, z)` | extra check |
|---|---|---|---|---|
| a | quintic, frame (0,0) | canonical | `1` | rows = project rows (all 25 x 10 x 6 residues) |
| b | quintic, frame (0,0) | stored | `1 - (43397/807975000) z^5 - ...` | rows = project rows x `(1 + m5 u^4 zhat^5)`, `m5 = -43397/807975000` (125 slots differ, all at `zhat^n`, `n >= 5`) |
| c1 | `F = 1 + u/3 - 2u^2/7 + 5u^3/11`, `Gm = 1 - 3u/5 + u^2/2 + 2u^3/3` | canonical | `1` | |
| c2 | `F = 3/2 + u/3`, `Gm = 1 - 3u/5 + u^2/2` (`lambda = 3/2`) | canonical | `2/3` | |
| d | as c1 | stored | `1 + (31764344287/1163868750) z^5 + ...` | |
| e | quintic, frame (1/3,-2/7) | canonical | `1` | rows = project rows in that frame (25 x 10 x 6) |

(c3, c4 = un_e2 cases 6 and 8 in canonical normalisation passed the same input checks; not run at genus 4.)

### 2.2 Graph sums (VERIFIED, genus 4, 2^89-1; `e4_tab_*.json`)

| case | run | edge graphs | vertex graph | total |
|---|---|---|---|---|
| a | maxout 2 | 10/10 regular; the 11 rows equal `g4_split_g4_89_nat.json` | `(1/180000) t_c^-6 + 0 + .. + 0` | |
| a | **full** | **378/378 regular**; all 379 rows equal `g4_split_g4_89_nat.json` (the natural run of record) | `(1/180000) t_c^-6 + 0 + .. + 0` | same |
| c1 | maxout 2 and **full** | **378/378 regular** | `(1/180000) T^-6 + 0 + .. + 0`, `T = Int f` | `(1/180000) T^-6 + 0 + .. + 0` |
| c2 | maxout 2 and **full** | **378/378 regular** | `(1/180000) (T/lambda)^-6 + 0 + .. + 0`, `lambda = 3/2` | same |
| b | maxout 2 and **full** | **exactly 6 of 378 singular**: graphs 1, 2 from `[u^4]`; 3, 4, 7, 8 at `[u^5]`; the other 372 rows equal case a | `t_c^-2`, `t_c^-1` terms | `(1/180000) t_c^-6 + 0 + .. + 0` (unchanged) |
| d | maxout 2 | **6 of 10 singular**, same graphs and orders | `T^-2`, `T^-1` terms | unchanged (s.2.3) |

In the quintic coordinate `t_c` the c1/c2 vertex graphs are NOT pure (`(1/180000) t_c^-6 - (1/180000) t_c^-5 + ...` for c1): the
polar part is carried by the system's own vanishing coordinate. Genus-4 outer graphs with at most two edges: 0 vertex `(4)`;
1 bridge `(1,3)`; 2 bridge `(2,2)`; 3 loop on `(3)`; 4 `(0,3)` loop at the genus-0 vertex + bridge; 5, 6 `(1,1,2)` chains; 7 `(1,2)` loop
at the genus-1 vertex + bridge; 8 `(1,2)` double edge; 9 `(1,2)` bridge + loop at the genus-2 vertex; 10 two loops on `(2)`.

### 2.3 The normalisation response (PROVED + VERIFIED)

*Valuation argument (PROVED).* If the rows change by the right scalar factor `1 + m u^4 zhat^5` mod `(zhat^10, u^6)` and `E_out = O(u)`,
then `xhat` (from `zhat^0`) and `E_M` are unchanged, `E_R` changes by `O(u^4)` and the leaves `T_R(d)` by `O(u^4)` (only `d >= 6`). A piece
landing on an outer graph with `E` edges has valuation `>= E` (each outer edge `O(u)`, everything else `O(1)` relative to `u^(2-2g)`);
replacing one outer edge by its change raises the bound to `E + 3`, one leaf to `E + 4`; products of two changes are `O(u^8)`. Hence the
difference is linear in `m`, confined to the vertex graph and the graphs with `E <= 2`, at `[u^4], [u^5]` for `E <= 1` and `[u^5]` for
`E = 2`. (Checked on the data: `un_e4.py report`, PASS.)

*Differences divided by `m5` (VERIFIED; exact rationals):*

| outer graph | (b - a)/m5 `[u^4]`, `[u^5]` | (d - c1)/m5 `[u^4]`, `[u^5]` |
|---|---|---|
| 0 vertex | `-1/8400000`, `-187/680400000` | `-1/8400000`, `3673/567000000` |
| 1 bridge `(1,3)` | `-1/1008000`, `13/6804000` | `-1/1008000`, `1/280000` |
| 2 bridge `(2,2)` | `1/900000`, `-53/32400000` | `1/900000`, `-407/40500000` |
| 3 loop on `(3)` | `0`, `-41/68040000` | `0`, `169/17010000` |
| 4 `(0,3)` | `0`, `-11/34020000` | `0`, `1/850500` |
| 7 `(1,2)` | `0`, `49/16200000` | `0`, `-2/253125` |
| 8 `(1,2)` | `0`, `-17/8100000` | `0`, `-13/4050000` |
| **sum** | **`0`, `0`** | **`0`, `0`** |

The `[u^4]` column is the response of the frozen (`u^0`) theory and is the same for both systems; the sums vanish, so the total
potential's polar part does not depend on the normalisation: the stored normalisation keeps the gap of the total and only moves a
`t_c^-2`, `t_c^-1` part out of the vertex graph into six edge graphs. Vertex difference in `t_c`: `(-1/8400000) t_c^-2 -
(751/1701000000) t_c^-1`
(b - a), `(-1/8400000) t_c^-2 + (7157/1134000000) t_c^-1` (d - c1), per unit `m5`.

*Reading.* The stored gauge `P^stored = P^canon (1 + m5 z^5 + ..)` is still a Theorem-U gauge between the same two cyclic systems
(flatness holds), it is symplectic in the window, and its only defect is the failure of annihilation, the cleanest
possible annihilation violation: `P^stored Y^C = (1 + m5 z^5 + ..) Y^C_Q`, whose positive part
`m5 (z^5 Y^C_Q)_+ = -m5 (z^5 e_1 + t_c z^4 e_0 + ..)` is non-zero. So **flatness is not sufficient; the canonical normalisation (=
annihilation at `z^5`) is needed
from genus 4 on**, as predicted.

## 3-4. E3 and E5

See `UNIVERSALITY_E3_E5.md` (s.3: first order, non-linear dressings at genus 2-3, consistency, case A at genus 4, dressings of
the pure model; s.4: E5, the genus-3 surplus at genus 4).

## 5. Conclusions

* **(A) regularity + (G) only.** VERIFIED sufficient: genus 2, 3 non-linearly (A1-A4, two primes) and to first order on the whole
  class (two primes); genus 4 for two random dressings, all 378 edge graphs (s.3.4; 2^89-1). Necessary in the sense: violating
  annihilation (B1, B2, E4 b, d) or flat regularity (C1a-c) breaks
  regularity; at first order regularity = annihilation at genus 2, and the genus-3 surplus is singular at genus 4.
* **(B) any Theorem-U gauge in canonical normalisation.** VERIFIED: genus 2, 3 (six perturbed systems + two quintic frames, un_e2;
  there the normalisation is invisible) and
  genus 4 (c1, c2 full; a = the natural run). Non-canonical normalisation: singular at genus 4 exactly where predicted. (B) is the
  special case of (A) in which `D` is a gauge; the canonical normalisation is exactly (G).
* **(C) quintic-specific.** Nothing: the cyclic-system data `f`, `g` can be perturbed arbitrarily, and the gauge can be replaced by an
  arbitrary non-gauge dressing in class (A). The quintic enters only through `B = prod b_a` (fixing `tau`, `v(u)` and `M`; not varied
  here) and through the fact that its actual R-matrix factors with `R'` in class (A) (Theorem U + (G) + canonical normalisation (H)).

**Conjecture A (annihilation principle; CONJECTURAL, supported as above).** Let `P''(u, z)` be symplectic, Z/5-graded, with
coefficients regular at `u = 0` in the flat frame, such that `(P'' Y^C(v(u)))_+ = 0`. Then in the hatted outer graph sum of
`P''.Omega^LG_{v(u)}` every stable graph with an edge is regular at `u = 0`, and the vertex graph is `c_g T''^(2-2g)` + regular,
`T'' = -[z^-1 e_0] P'' Y^C`, `c_g = |B_2g|/(2g(2g-2) 5^(g-1))`. With (G) for the quintic's canonical gauge (gauge worker) this implies
Conjecture N, `T'' = t_c`. A proof would have to show that the polar parts of the model correlators pair to zero with every edge of
such a `P''`; the first-order computation of s.3.1 is the finite shadow of that statement.

## 6. Accounting of runs

Genus-4 runs, all at 2^89-1 on the frame-(0,0) base (CPU seconds as logged by the pool; the `fq4_split` reference run took 3251 CPU s):
* full (all 379 inner graphs, 4 workers): c1 1374, c2 1419, a 1470, b 1487, A seed 20261113 1576, A seed 2 1654 -- six full runs;
* maxout 2 (about a sixth of the work each): a 511, b 538, d 631, c1 630, c2 605 (5 workers); A 537, E5 `w_1` 615, `w_2` 577 (one process).
Total about 13600 CPU s = 227 CPU min, below the six full-size runs (6 x 54 CPU min) allowed; at most 6 processes at any time.

## 7. Files

| file | content |
|---|---|
| `un_core.py`, `un_engine.py` | rows -> edge/leaf tables; the engine `Var` (dual numbers, maxout; pool capped at 5 workers) |
| `un_gauge.py` | Theorem-U gauges of perturbed cyclic systems (in memory), rows from a gauge, (G) check |
| `un_e1.py/.log`, `un_e2.py/.log` | E1 regression, E2 genus 2-3 universality (first attempt; re-run today at 2^89-1) |
| `un_e4_cases.py`, `un_e4.py`, `e4.log`, `e4_tab_<case>_89[_out2].json` | E4/U4 genus 4 |
| `un_e3_resp.py/_89.log`, `un_resp_89.json`, `un_lin.py`, `un_e3_lin.py/.log/_89.json` | E3 first order |
| `un_ann.py`, `un_e3_nl.py/_89.log`, `un_e3_c.py/_89.log` | E3 non-linear, genus 2-3 (case C1 of `un_e3_nl.py` crashed by design, see s.3.2) |
| `un_e3_g4.py`, `un_e3_g4.log`, `e3g4_tab_*.json` | E3 case A at genus 4 |
| `un_e5.py`, `un_e5.log` | E5: the genus-3 surplus directions at genus 3 and 4 |
| `un_e3_model.py`, `un_e3_model_89.log`, `un_e3_model_61.log` | annihilating dressings of the pure model (s.3.5) |
| `un_e3_nl_61.log`, `un_e3_c_61.log`, `un_resp_61.json`, `un_e3_resp_61.log`, `un_e3_lin_61.json` | second prime 2^61-1, genus 2-3 |
| `UNIVERSALITY.md`, `UNIVERSALITY_E3_E5.md` | this note (two pieces) |
