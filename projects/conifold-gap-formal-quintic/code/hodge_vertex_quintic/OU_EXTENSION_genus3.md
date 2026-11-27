# Formal quintic at genus 3: the tautological factorisation and its u-corrections, graph by graph

Date 2026-11-09 (s.6d: 2026-11-10; s.6e: 2026-11-11). Exact finite computation (two primes 2^89-1, 2^61-1; rational reconstruction at 2^89-1 confirmed by the residues at
2^61-1 for steps 1-5b at frame (0, 0); the second frame of s.7 is single-prime), AI-checked only. Code: `fq3_graphs.py`, `fq3_mvertex.py`, `fq3_m.py`, `fq3_msum.py`, `fq3_par.py`, `fq3_steps.py`, `fq3_rec.py`,
`fq3_taut.py`, drivers `genus3_ou.py` (steps 1-4), `genus3_step5.py`, `genus3_step5b.py`, `genus3_step5b_post.py`, `genus3_cross.py`, and for s.6c
`fq3_quad.py`, `fq3_quadlib.py`, `genus3_quad.py`, `genus3_quadan.py`, `genus3_quadsolve.py`, `genus3_quadsolve2.py` (logs `logs/genus3_quad_*.log`, data `g3_quad_*.json`); logs `logs/genus3_ou_89.log`, `logs/genus3_ou_61.log`,
`logs/genus3_ou_89_ab1o3_m2o7.log`, `logs/genus3_cross.log`; states `g3_state_89.json`, `g3_state_61.json` (residues of every number
below). Referee reports: `g3_referee_modules.json`, `g3_referee_mvertex.json`, `g3_referee_results.json`, `g3_referee_final.json`,
`g3_referee_quad.json` (s.6c), `g3_referee_u3.json`, `g3_referee_u3b.json` (s.6d). This note continues OU_EXTENSION_genus2.md (genus 2) and uses
its conventions: formal quintic of LP18, `zeta` primitive 5th root of unity, `u = 1 - 5^5 q`, `s = u^{1/5} = 1/L`, `zhat = L^4 z`, conifold
frame `(alpha, beta) = (0, 0)` unless stated, `t_c = u + (7/10) u^2 + (41/75) u^3 + ...`.

## 0. Summary

1. **Machinery.** The hatted theory of OU_EXTENSION_genus2 s.1 extends verbatim: for ANY symplectic `M(zhat; u)` (now mod `zhat^7`, since
   the psi-deficit reaches 3g-3 = 6) every genus-3 graph term is `u^{-4} G_Gamma(u)` with `G_Gamma` a u-series; a graph is regular at the
   conifold point iff `[u^0] = [u^1] = [u^2] = [u^3]` of `G_Gamma` vanish. 42 stable graphs (1 vertex graph, 2 one-edge, 5 two-edge, 9
   three-edge, 25 with E >= 4). The vertex CohFT `M.omega^` is computed to M-order 6 (537 distinct vertex tensors per M); all 42 graphs with all
   psi-deficit decorations and colourings; one constant-M pass over the 42 graphs takes ~80 s, one u-dependent pass 19-28 min (17-42 graphs) on 24
   processes.
2. **Point frame (M = 1).** `F_3^cf u^4 = 1/25200 - u/9000 + 407 u^2/3780000 - 151 u^3/3780000 - 29663 u^4/9072000000 + 364981 u^5/75600000000`,
   i.e. in `t_c`: **`(1/25200) t_c^{-4} + 0 t_c^{-3} + 0 t_c^{-2} + 0 t_c^{-1} - 57737/8400000000 + ...`**: the genus-3 conifold gap holds
   for the formal quintic graph sum, with `1/25200 = (B_6/(6 . 4))/5^2` and `1/1200 = (-B_4/(4 . 2))/5`, i.e. `(-1)^{g-1} B_{2g}/(2g(2g-2))
   5^{-(g-1)}` at g = 2, 3: the gap coefficients in the coordinate `sqrt(5) t_c` (the factor `5^{-(g-1)}` is a rescaling of the string
   coupling by `1/5` in the formal-quintic normalisation). The total is independent of
   M (checked at M = 1, M_0 and a random order-6 M, through u^5).
3. **Tautological M_0 = R~ mod zhat^7** (the conifold symbol in Teleman's colour form, CONIFOLD_SYMBOL s.2, now to order 6): every edge
   graph loses its `u^0` term, an E-edge graph is `O(u^E)`, so the 25 graphs with E >= 4 are regular with M_0 alone; the vertex graph has
   `[u^0] = 1/25200` (the whole leading pole). The 16 graphs with E <= 3 are not regular (table in s.3).
4. **O(u): the covariant-scalar class FAILS.** With `M = M_0 exp(u N_1)` the `[u^1]` conditions on the two one-edge graphs are two
   independent linear equations on `N_1` (sp mod zhat^7, 75 coordinates). On the Z/5-covariant scalar class `N_1 = b_1 diag(zeta^{-i}) zhat
   + b_3 diag(zeta^{-3i}) zhat^3 + b_5 zhat^5` the genus-3 equations alone are consistent (a 1-parameter family), but TOGETHER with the two
   genus-2 equations of OU_EXTENSION_genus2 (which fix `(b_1, b_3) = (113/75, -21/125)`) the 4 x 3 system has rank 3 and is INCONSISTENT
   (right and left placement). The same holds for the 15-parameter diagonal class. The Z/5-covariant full-matrix (circulant) class is
   consistent (11-dimensional solution set) and so is the full class (71-dimensional); both contain the tautological first u-correction
   `N_1^taut = M_0^{-1} [u^1] R~(zhat; u)` of the symbol, which satisfies all four equations exactly (a check of the Jacobians and
   of g_1 at both genera through a different code path; it also makes the consistency of these two classes automatic, so the content
   of those rows is the rank 4). So **no q-dependent Z/5-covariant scalar (Hodge-type) factor placed on one side of M_0, `M = M_0
   exp(uS)` or `exp(uS) M_0`, regularises the edge graphs of genus 2 and genus 3 simultaneously, already at order u.** A scalar on both
   sides, `exp(u S_L) M_0 exp(u S_R)` (6 parameters), is not excluded at O(u): the 4 x 6 system has rank 4 and a 2-dimensional solution
   set, but it fails at O(u^2) at genus 3 (s.6b): the five two-edge `[u^2]` quadrics along the (effectively one-parameter) family
   have no common zero.
5. **O(u^2) at genus 3 alone.** Along the 1-parameter covariant-scalar family (right placement, frame (0, 0)) that solves the genus-3
   O(u) equations, the `[u^2]` coefficients of the five two-edge graphs are quadratics in the parameter with rational coefficients and
   no common zero over C (graph 6 has the two rational roots 2778636/528125, 2514114/1046875, at which the other four quadratics are
   non-zero; graphs 3, 8, 17 have irrational roots, graph 7 complex ones): the one-sided scalar class fails at O(u^2) even at genus 3 alone
   (s.6). Left placement and other frames were not run at this order.
6. **O(u^2) in the circulant and full classes (s.6c).** The five two-edge `[u^2]` coefficients are homogeneous quadrics in `Delta = N_1 -
   N_1^taut`; they vanish identically on the `zhat^{4,5,6}` components of `N_1` (a two-edge graph cannot carry two edges of psi-order >= 4),
   which absorb all but one condition in the circulant class and all of them in the full class (rational solution variety, covariance
   lost). The remaining circulant condition is a sextic on the six low-order directions compatible with genus 2 that FACTORS as
   `c_1(0) . m_1 . F_4` (the two rank-2 quadrics of graphs 6 and 17 split); `m_1 = 0` is a Q(zeta)-rational hyperplane of genuine solutions,
   so Q(zeta)-rational Z/5-covariant non-tautological `N_1` regularising all two-edge graphs at O(u^2) (and the one-edge graphs of both
   genera at O(u)) exist in a 4-parameter family; one was constructed explicitly and verified by an independent full pass. (At O(u^3)
   only `m_1 = 0` survives and the family is cut down, item 7.)
7. **O(u^3) (s.6d).** Beyond `N_1` the equations are linear (`N_2`: 7 equations, `N_3`: 2; where the `N_2` matrix drops rank its solvability
   is one more condition on `N_1`, satisfied wherever tested), so the new non-linear conditions are the nine three-edge cubics `c_Gamma(N_1 -
   N_1^taut)`. The star-graph cubic is a perfect cube, `c_32 = kappa m_1^3`, which forces `m_1 = 0`. All `N_1`-conditions are homogeneous in
   `Delta = N_1 - N_1^taut`; their Z/5-covariant solutions are, up to scaling and to affine fibres in the orders 4..6: the pure-high family
   (`zhat^{4,5,6}` parts subject to two linear conditions); a 2-dimensional linear space Pi, `Delta = M_0^{-1} (zhat D_1) M_0` with `D_1` a
   constant psi-free propagator `beta (1 (x) H^2 + H^2 (x) 1) + gamma (H^3 (x) H^4 + H^4 (x) H^3)`, which is kinematic (the conifold-limit
   vertex theory has higher-genus primary correlators only in the direction `H`; `beta` is a string-equation gauge direction, `gamma` is
   explained at leading order only), the third covariant propagator `H (x) H` being the one the conditions see; and 17 isolated rays (16
   with `c_1(0) = 0`: 13 at generic position, a conjugate pair and one rational point; and one rational point with `c_1(0) != 0`). For four
   representatives (pure-high, a point of Pi, one of the 13, the point with `c_1(0) != 0`) `N_2`, `N_3` were solved and a pass over all 42
   graphs confirms, at both primes, that ALL 41 edge graphs are regular, so that the vertex graph carries the whole polar part
   `(1/25200) t_c^-4`; the genus-2 edge graphs are regular as well. So the regularising covariant factorisation is not unique at genus 3,
   but the only positive-dimensional covariant families of first-order corrections are the kinematic plane and the pure-high family.
8. **Natural form (s.6e, NATURAL_FORM.md).** `M = M_0(tau(u) zhat)`, `1/tau = 2F1(1, 4/5; 9/5; u)` (the symbol in the variable `z` over the
   canonical coordinate measured from the conifold point), with no free parameter, makes all 41 genus-3 edge graphs and all six genus-2
   edge graphs regular (frame (0, 0), both primes; a second frame at one prime). `M_0` is the R-matrix of the Landau-Ginzburg model
   `x^4 + v/x` with primitive form `dx`; the other factor is a gauge regular at the conifold point. All-genus statement: conjectural.
9. Frame independence `(alpha, beta) = (1/3, -2/7)`: the O(u) system, its solutions and the gap statement are frame independent; the
   per-graph polar parts with `M_0` (the leading `[u^E]` of 30 of the 41 edge graphs) are not, so the O(u^2) and O(u^3) statements of s.6, s.6b, s.6c, s.6d are
   frame-(0, 0) statements (s.7). Cross-prime: all 412 reliably reconstructed values and all
   structural data agree at 2^61-1 (s.8).

## 1. Hatted theory at genus 3

`Hhat^a = s^{1-c_a} H^a`, `c = (0, -1, 3, 2, 1)`, `z = s^4 zhat`; `Qhat[a][i][n] = Q[a][i][n] s^{4n+1-c_a}` are pure u-series with
`Qhat|_{u=0} = R^` (the conifold symbol), `xhat = varpi_0^cf(u) = 1 + O(u^3)`. The s-exponent of a genus-g decoration is `-10(g-1)`
(OU_EXTENSION_genus2 s.1), so at genus 3 `T_Gamma = u^{-4} G_Gamma(u)`; in `t_c`: `t_c^{-4}: a_0`, `t_c^{-3}: a_1 + 4 c_1 a_0`, with
`c_1 = 7/10`, `t_c^{-2}: a_2 + 3 c_1 a_1 + (2 c_1^2 + 4 c_2) a_0`, ... (fq3_msum.tc_expansion3 composes the series; the displayed formulas
were checked exactly against the stored data, e.g. `407/3780000 + 3 (7/10)(-1/9000) + (2 (7/10)^2 + 4 (41/75))/25200 = 0`). The rows are needed to z-order 6 (`Nz = 6`: edges
`m + n + 1 <= 6`, leaves `d <= 7`); they are symplectic to order 6 as u-series. The code works in u directly (precision `O(u^6)`,
`prec_abs = 34` in s).

**Checks (all PASS at both primes):** hatted rows pure u-series with `Qhat|_{u=0} = R^`; `M_0` symplectic to order 6, `exp(log M_0) = M_0`;
direct s-series sums (SRows, M = 1) have valuation -20 and equal `s^{-20} G_Gamma(s^5)` for the graphs with E <= 1; frozen (u = 0) sums =
`[u^0]` of the hatted sums at M = 1 and at a random M; the total is M-independent through `u^5` (M = 1, random order-6 M, M_0); the
genus-2 graphs run through the genus-3 code reproduce the genus-2 tables (point frame, M_0, total through u^4); the genus-2 Jacobian
recomputed with the new code equals `fq_ou.jacobian_N` on the 40 genus-2 coordinates and vanishes on the 35 coordinates of orders 4-6.

## 2. The vertex CohFT to order 6 and the graph sum (new code, and a bug in the old one)

* `fq3_graphs.legged_graphs(g, n, Emax)`: stable graphs with legs by vertex splitting (non-separating / separating degenerations) with
  canonical forms; equals the Mathematica exports (7 at genus 2, 42 at genus 3) and brute force (fq_mvertex.legged_graphs) on all
  comparable cases; |Aut| checked by an independent half-edge brute force (referee). Sizes: (3,0): 1,2,5,9,12,8,5 per edge number;
  (2,2,E<=5): 75; (3,1,E<=5): 135; (1,4,E<=4): 163; (0,6,E<=3): 236.
* `fq3_mvertex.MVertex6`: `OmegaHat^M_{g,n}(i_1..i_n; d)` with psi-deficit `D <= 6`, inner edges `Ehat[(a,b)]` with `a + b + 1 <= 6`,
  leaves `That[e]` with `e <= 7`, legs `M^{-1}_k`, `k <= 6`; decorations grouped by (leg orders, half-edge exponents), legs contracted
  vertex by vertex. Checks: equality with the genus-2 `fq_mvertex.MVertex` on all tensors with `D <= 3`; diagonal
  `M = diag(exp(a_1 zeta^{-i} zhat + a_3 zeta^{-3i} zhat^3 + a_5 zeta^{-5i} zhat^5))` gives `zeta^{4i(2g-2+n)} Int exp(12 a_1 ch_1 - 720 a_3
  ch_3 + 30240 a_5 ch_5) prod psi^d` (23 tensors through genus 3, D = 6, off-diagonal entries zero); string and dilaton equations at a
  random full order-6 M at genus 3; a referee's plain nested-loop re-implementation agrees entry by entry on (2,(0,)), (2,(0,0)),
  (2,(1,0)), (3,()), (3,(1,)), (3,(2,)), (3,(2,1)), (1,(0,0,0)) and on slices of (1,(0,0,0,0)), (0,(0^6)).
* **Bug found and fixed.** `fq_msum.contract_pair`, `fq_msum.merge` and the leaf contraction in `fq_msum.GraphSummer.vertex_tensor`
  write their output tensors in `itertools.product` order (last index fastest) while every reader uses position `= Sum_l i_l 5^l` (first
  index fastest). At genus 2 no asymmetric intermediate tensor of rank >= 2 is ever read back (checked graph by graph by a referee; the
  genus-2 production values are unaffected), but at genus 3 they are (first at graph 17, the chain g1-g1-g1 with the middle vertex
  processed last): the total was M-dependent and the point-frame total had a non-zero `t_c^{-3}` coefficient. `fq3_msum.py` carries
  corrected versions (the genus-2 module is left untouched). Validation of the fix: the composition law `R'.(M_2.omega^) = (R' M_2).omega^`
  holds graph by graph for random full order-6 `R' = M_1 M_2 M_2^{-1}`, i.e. the outer sum with rows `M_1 M_2` and vertices `M_2.omega^`
  equals the inner tensor `MVertex6(M_1 M_2).tensor(3, ())` (42 graphs, two primes), and the M-independence checks of s.1 pass.
* `fq3_m.py`: sp mod zhat^7 (75 coordinates: `N_k` symmetric for odd k, antisymmetric for even k; the first 40 are the genus-2
  coordinates); classes `cov` (b_1, b_3, b_5), `diag` (15 Fourier modes), `circ` (Z/5-covariant: `(N_k)_{ij} = zeta^{-ki} c_k(j-i)`,
  `c_k(5-d) = (-1)^{k+1} zeta^{kd} c_k(d)`, 15 parameters, verified to span exactly the covariant subspace), `full`; `M(zhat; u) = M_0
  exp(Sum_k u^k N_k)` as a bivariate truncated exponential (right placement; left: `exp(..) M_0`); u-dependent vertex tensors by
  interpolation in eps over 21 nodes (eps-degree <= 18 asserted at run time).

## 3. Pass A: the tautological M_0 (reporting prime, confirmed at 2^61-1)

`M_0 = R~ mod zhat^7`, `(M_0)_k = C_k`, `(M_0^{-1})_k = (-1)^k C_k^T`; `R' = R^ M_0^{-1} = Psi + O(zhat^7) + O(u)` with `Psi` the constant
orthonormal change of basis flat -> colour (CONIFOLD_SYMBOL s.3), so every edge factor is `O(u)` (`Psi eta^{-1} Psi^T = eta^{-1}`). Hence every E-edge graph is
`O(u^E)` and `[u^0]` vanishes on all 41 edge graphs (checked). Vertex graph: `[u^0] = 1/25200`; in `t_c`: `(1/25200) t_c^{-4} + (20083/40000000)
t_c^{-3} + (129397/200000000) t_c^{-2} + (7813571/20250000000) t_c^{-1} + ...` (the subleading poles are cancelled by the edge graphs).
Non-regular graphs with M_0 alone (u^0..u^3; Mathematica indices, `gv` = vertex genera):

| graph | E | gv | u^0 | u^1 | u^2 | u^3 |
|---|---|---|---|---|---|---|
| 2 (loop on g2) | 1 | (2) | 0 | -1463/600000 | 4551/1562500 | -4649837/4500000000 |
| 5 (bridge g1-g2) | 1 | (1,2) | 0 | 232351/120000000 | -406121/240000000 | 1211567/2160000000 |
| 3 | 2 | (1) | 0 | 0 | -303/250000 | -2621/37500000 |
| 6 | 2 | (0,2) | 0 | 0 | -28329/50000000 | 141453/500000000 |
| 7 | 2 | (1,1) | 0 | 0 | 35267/22500000 | -115511/150000000 |
| 8 | 2 | (1,1) | 0 | 0 | -29153/45000000 | 197317/450000000 |
| 17 | 2 | (1,1,1) | 0 | 0 | 8263/180000000 | -18559/1080000000 |
| 4 | 3 | (0) | 0 | 0 | 0 | -27/125000 |
| 9 | 3 | (0,1) | 0 | 0 | 0 | 587/562500 |
| 10 | 3 | (0,1) | 0 | 0 | 0 | -2797/3750000 |
| 11 | 3 | (0,1) | 0 | 0 | 0 | 807/1250000 |
| 12 | 3 | (0,1) | 0 | 0 | 0 | -1319/7500000 |
| 18 | 3 | (0,1,1) | 0 | 0 | 0 | 127/11250000 |
| 19 | 3 | (0,1,1) | 0 | 0 | 0 | -59/1000000 |
| 20 | 3 | (0,1,1) | 0 | 0 | 0 | -79/2812500 |
| 32 | 3 | (0,1,1,1) | 0 | 0 | 0 | -1/4050000 |

All values are rational (colour sums of a Z/5-covariant M_0), with denominators of the form `2^a 3^b 5^c`. The 25 graphs with E >= 4 are
regular. Totals: `Sum_Gamma [u^1] = -1/9000 = -4 c_1/25200` (the `t_c^{-3}` cancellation) as at M = 1.

## 4. The O(u) system

`[u^1] G_Gamma[M_0 exp(u N_1)] = g1_Gamma + J[Gamma] . n_1` with `J[Gamma][v] = d/d eps|_0 G^0_Gamma[M_0 exp(eps N_v)]` (frozen theory,
75 directions, eps-polynomials of degree <= 6 (the bound is the weight count 3g-3 = 6; exact interpolation from 7 nodes, so the
logged degree PASS is tautological), right placement; `J` depends on the symbol only). `J` vanishes
identically on the E = 2 control graphs 3 and 6, `J[1] + J[2] + J[5] = 0`, rank of the 2 x 75 matrix (graphs 2, 5) = 2. `J[2]` has 65
non-zero entries, all Galois-reconstructible in Q(zeta) at the 1e10 bound (e.g. `n1_00: 11/400000`, `n1_01: -27501/256000000 (1 + zeta) -
1637723/4096000000 zeta^3`); `J[5]` has 75 non-zero entries, 55 reconstructible at that bound (`n1_00: -3971/80000000`, `n3_00:
59/1920000`; the covariant column entries are 5 x the diagonal ones, e.g. `5 . 59/1920000 = 59/384000`). Consistency check (same rows, different code path): the tautological
`N_1^taut = M_0^{-1} [u^1] R~(zhat; u)` (for which every edge graph is `O(u^2)` by construction) satisfies `J . n = -g1` for the genus-2
graphs 2, 4 (`-1/750`, `13/8000`) and the genus-3 graphs 2, 5 (`1463/600000`, `-232351/120000000`) exactly at both primes; its scalar
(diagonal covariant) part is `(b_1, b_3, b_5) = (-3/25, 3/625, 0)`.

Restricted to the covariant-scalar class (right placement; rows = equations, columns = b_1, b_3, b_5):

| equation | b_1 | b_3 | b_5 | rhs |
|---|---|---|---|---|
| genus 2, graph 2 (g1 loop) | -1/80 | -5/48 | 0 | -1/750 |
| genus 2, graph 4 (g1-g1) | 1/640 | 5/1152 | 0 | 13/8000 |
| genus 3, graph 2 (loop on g2) | 11/80000 | 1/9600 | 1/640 | 1463/600000 |
| genus 3, graph 5 (bridge g1-g2) | -3971/16000000 | 59/384000 | -1/46080 | -232351/120000000 |

The genus-2 rows give `(b_1, b_3) = (113/75, -21/125)` (OU_EXTENSION_genus2); the genus-3 rows then require `b_5 = 13492/9375` and
`b_5 = 5531436/78125` respectively: rank 3, inconsistent (left placement: the genus-2 rows give `(b_1, b_3) = (-113/75, 2/15)` (OU_EXTENSION_genus2) and the
genus-3 rows then require `b_5 = 674/625` and `4288718/46875`). Genus 3 alone: `(b_1, b_3, b_5) = (27664/2255, 36974/5125, 0) + t (-41875/8118, -2015/246,
1)`. Left placement (`exp(u N') M_0`): rows `(-1/120, -5/48, 0 | -1/750)`, `(-1/1440, 5/1152, 0 | 13/8000)`, `(-1/2000, 0, 1/640 |
1463/600000)`, `(-1/36000, 1/18000, -1/46080 | -232351/120000000)`: also rank 3, inconsistent; genus 3 alone: `(-1463/300,
-2237459/60000, 0) + t (25/8, 125/64, 1)`.

| class | parameters | genus 3 alone | genus 2 + 3 |
|---|---|---|---|
| cov (right / left) | 3 | rank 2, 1-dim family | rank 3, INCONSISTENT |
| diag Fourier modes (right / left) | 15 | rank 2, dim 13 | rank 3, INCONSISTENT |
| circulant (right / left) | 15 | rank 2, dim 13 | rank 4, dim 11 |
| full (right / left) | 75 | rank 2, dim 73 | rank 4, dim 71 |

(identical at both primes; over Q(zeta) the ranks are those over F_p unless p divides a denominator, excluded by the two primes). The
diagonal class fails exactly like the scalar one: at O(u) only the covariant Fourier modes `X_{1,1}, X_{3,3}, X_{5,0}` enter the equations
(their columns are `1/5` of the covariant ones), so the restricted system is the scalar one. This says nothing about the diagonal class
at O(u^2), where charge-neutral products of non-covariant modes enter.

## 5. Interpretation

At genus 2 the two O(u) equations had the unique covariant-scalar solution `(113/75, -21/125)` and it made all six edge graphs regular.
At genus 3 the corresponding equations for the two one-edge graphs are not compatible with that solution, and no `b_5` can repair this:
the "Hodge-type" q-dependent scalar factorisation of the symbol does not extend to genus 3 within the one-sided covariant-scalar class;
the genus-2 solution is not the first term of a scalar correction that regularises both genera.
What does persist is the tautological, colour-mixing correction `N_1^taut` (the actual first u-correction of the Teleman form of R),
which lies in the Z/5-covariant circulant class; the regularising factorisations of genus 2 + 3 at O(u) form an 11-dimensional affine
family inside that class (71-dimensional in the full class). Whether a non-tautological member of the circulant family also passes the
O(u^2), O(u^3) conditions of the two- and three-edge graphs is a system of quadratic and cubic equations in 11 unknowns (s.6 treats the
one-parameter scalar family at genus 3 alone; the question is answered in s.6c and s.6d: yes); the tautological member passes them
trivially since then `R' = Psi(u) + O(zhat^7) + O(u^2)`.

## 6. O(u^2) along the genus-3 covariant-scalar family

`genus3_step5.py` (log `logs/genus3_ou_89.log`, step (5); keys `A_cov5_right_t*`, `obs_u2_right_*` of the state): on the 1-parameter
family `(b_1, b_3, b_5)(t) = (27664/2255, 36974/5125, 0) + t (-41875/8118, -2015/246, 1)` solving the genus-3 O(u) equations (right
placement), four u-dependent passes `M = M_0 exp(u N_1(t))`, `t = 0, 1, 2, 3` (all 42 graphs at t = 0, the 17 graphs with E <= 3 otherwise;
eps-degree of the vertex tensors 6, bound 18). Checks at t = 0: M-independence of the total, the one-edge graphs vanish through `u^1`,
every graph is `O(u^E)`. The `[u^2]` coefficient of each two-edge graph is a polynomial of degree <= 2 in t (two affine edge factors; the
`t^3` coefficient of the 4-node interpolant vanishes: a genuine one-node redundancy check), the `[u^3]` coefficient of each three-edge
graph of degree <= 3 (three affine factors; not checkable with four nodes):

| graph (E = 2) | `[u^2](t)`: t^0, t^1, t^2 | roots mod 2^89-1 | roots mod 2^61-1 | over Q |
|---|---|---|---|---|
| 3 (g1, two loops) | -29316307753/23113750000, 3769939/1008600, -260979625/95857344 | none | two | D = 83488029443/527215392000 > 0, not a square: two real irrational roots |
| 6 (g0 loop + bridge to g2) | 32341702169/231137500000, -225954721/2662704000, 1415375/127809792 | 2778636/528125, 2514114/1046875 | the same | D a rational square: roots 2514114/1046875, 2778636/528125 |
| 7 (g1 = g1 double edge) | 840466283347/1830609000000, -5091835651/5857948800, 5788534625/12653169408 | two | none | D = -385384394149553/4555140986880000 < 0: two complex roots |
| 8 (g1 loop + bridge to g1) | -10982137777/13952812500, 54470011/35719200, -18991375/25717824 | two | none | D = 52926376057/94898770560000 > 0, not a square |
| 17 (chain g1-g1-g1) | 635251632541/27459135000000 (two-prime CRT), -16958212/411887025, 924120125/50612677632 | none | two | D = 1675487191/303676065792000 > 0, not a square |

(Coefficients: 1e13-bound reconstructions at 2^89-1 confirmed by the 2^61-1 residues of the step-5 run there (`logs/genus3_ou_61.log`,
which repeats all four passes without the control pass); the F_p root data are prime-specific, the invariant statements are the
discriminants.) Any common zero over C would be a root of graph 6, i.e. one of its two rational roots, and at both of them graphs 3, 7, 8,
17 are non-zero (e.g. graph 3 takes `-201348512665373/3534423750000` and `-560712332019/70140625000`; referee check, also all ten pairwise
resultants are non-zero at both primes). **Hence no one-sided covariant-scalar `N_1` that fixes the genus-3 one-edge graphs at O(u) makes
all two-edge graphs regular at O(u^2)**: this class fails at genus 3 even without the genus-2 constraint. A structural consistency check
of the step-5 data: at graph 6's root `t = 2514114/1046875` the family has `b_1 = -3/25 = c_1(0)^taut` exactly, as it must (the loop at
the trivalent genus-0 vertex has psi-deficit 0 and sees only the colour-diagonal of the `zhat^1` block, which vanishes at the tautological
value). The one FAIL line of the step-5 log ("KEY ... common roots: none") is this negative result, not an error. The `[u^3]` polynomials of the nine three-edge
graphs are also listed in the log (graph 32: `-1/4050000`, independent of t; graphs 18, 20: cubic coefficient 0).

Control: with the tautological `N_1^taut` (every edge graph `O(u^2)` by construction) the full pass gives exactly that: all 41 edge graphs
vanish through `u^1`, the 39 with E >= 2 are regular, the one-edge graphs have `[u^2] = 2623/1500000` (graph 2) and `-129629/72000000`
(graph 5) and would require `N_2` (and then `N_3`) - i.e. the next u-corrections of the symbol itself, which is the trivial factorisation
`M = R~ mod zhat^7`, `R' = Psi + O(zhat^7)`.

## 6b. The two-sided covariant-scalar dressing at O(u) and O(u^2)

(Prompted by the results referee.) `M = exp(u S_L) M_0 exp(u S_R)` with Z/5-covariant scalars `S_L = Sum_{k odd} b'_k diag(zeta^{-ki})
zhat^k`, `S_R = Sum_{k odd} b_k diag(zeta^{-ki}) zhat^k` (six parameters) has `N_1 = M_0^{-1} S_L M_0 + S_R` at O(u). The 4 x 6 system of the
genus-2 + genus-3 O(u) equations (`genus3_step5b.py`, log step (5b)) is

| equation | b_1 | b_3 | b_5 | b'_1 | b'_3 | b'_5 | rhs |
|---|---|---|---|---|---|---|---|
| genus 2, graph 2 | -1/80 | -5/48 | 0 | -1/120 | -5/48 | 0 | -1/750 |
| genus 2, graph 4 | 1/640 | 5/1152 | 0 | -1/1440 | 5/1152 | 0 | 13/8000 |
| genus 3, graph 2 | 11/80000 | 1/9600 | 1/640 | -1/2000 | 0 | 1/640 | 1463/600000 |
| genus 3, graph 5 | -3971/16000000 | 59/384000 | -1/46080 | -1/36000 | 1/18000 | -1/46080 | -232351/120000000 |

(the left columns are those of the left-placed scalar class of s.4): rank 4, consistent, solution set `(1600063/252735, -477191/421225,
5523764/2106125, 2032126/421225, 0, 0) + t_1 (-5375/16849, -15774/16849, -977/84245, -5375/16849, 1, 0) + t_2 (0, 0, -1, 0, 0, 1)`. The
`t_2` direction moves the central element `zhat^5 . 1` from the right factor to the left one and leaves `M` unchanged, so the family is
effectively one-dimensional. Seven u-dependent passes (`(t_1, t_2) = (0,0), (1,0), (2,0), (0,1), (0,2), (1,1)`, check point `(2,1)`) on the
graphs with E <= 2: the one-edge graphs vanish through `u^1` at every point (the O(u) solve is confirmed), every graph is `O(u^E)`, the
`[u^2]` coefficients of the five two-edge graphs are quadratics (degree <= 2 confirmed at the check point) with no `t_2` dependence.
Univariate analysis (`genus3_step5b_post.py`, log step (5b post)): the quadratics `Q_Gamma(t_1)` have pairwise `gcd(Q_3, Q_j) = 1` over F_p
for j = 6, 7, 8, 17, hence no common zero over the algebraic closure of F_p and, the coefficients being rational, none over C; mod 2^89-1 graphs 7 and 8 have no root at all, graph 6 has the
rational roots `-4710677392/18216875`, `7124167/403125`, graph 17 the rational roots `2200616/134375`, `-28895124/8111875`, graph 3 two
non-rational roots. The two residues determine all 15 coefficients by CRT (denominators `2^a 3^b 5^c 16849^2`, `16849 = 7 . 29 . 83` the
denominator of the `t_1` direction): over Q the discriminants of graphs 3, 7, 8 are positive non-squares and graphs 6, 17 have exactly
the rational roots above (all four also roots mod 2^61-1), so the gcd statement holds over C (final referee's check). **So a covariant scalar on both sides of `M_0` solves the four O(u) equations but fails at O(u^2) at genus 3 (frame
(0, 0)).** The 2^61-1 run of step 5b (`logs/genus3_ou_61.log`) gives the same structure (all four gcds constant; the rational roots
`7124167/403125` of graph 6 and `2200616/134375` of graph 17 reappear; the discriminants of graphs 3, 7, 8 are non-squares there, that of
graph 3 a square at 2^89-1, i.e. its roots are irrational). The step-5b log line "FAIL KEY ... gcd of resultants
= constant" is an artefact (the resultant in `t_2` of `t_2`-independent polynomials vanishes identically); the valid KEY line is the one of
step (5b post).

## 6c. The circulant and full classes at O(u^2)

(`fq3_quad.py`, `fq3_quadlib.py`, `genus3_quad.py`, `genus3_quadan.py`, `genus3_quadsolve.py`, `genus3_quadsolve2.py`, `genus3_quadsolve3.py`;
logs `logs/genus3_quad_89.log`, `logs/genus3_quad_61.log` - their first (D)/(E) blocks come from superseded runs (a random-combination bug and
a missing particular solution, both fixed; the stale keys `lambda`, `R_rank`, `R_matrix`, `quadsolve_examples` in `g3_quad_*.json` belong to
them); data `g3_quad_89.json`, `g3_quad_61.json`; referee report `g3_referee_quad.json`, whose two referees found the factorisation below.)

**Structure.** For `M = M_0 exp(u N_1)` every primed edge factor is `u E'_1 + O(u^2)` with `E'_1` affine in `N_1`, leaves are `O(u)`, the vertex
tensors are `Omega^{M_0} + O(u)` and `xhat = 1 + O(u^3)`. Hence the `u^2` coefficient of a two-edge graph is a quadratic polynomial
`Q_Gamma(N_1)` in the 75 sp coordinates (two edge factors, vertex tensors at `M_0`, no leaves), homogeneous in `Delta = N_1 - N_1^taut` since
`E'_1` vanishes at the tautological correction. The five polynomials were computed exactly by polarisation (2926 evaluations with the vertex
tensors frozen at `M_0` and the rows to `O(u^2)`, 3 s each) and checked against the stored full passes: `Q_Gamma(0)` = pass A, `Q_Gamma(N_1(t))`
= the step-5 passes at t = 0..3, `Q_Gamma(N_1^taut) = 0`, vanishing gradient at `N_1^taut`, a direct evaluation at a random `N_1`. The one-edge
graphs at `u^2` are affine in `N_2` through the same Jacobian `J` (rank 2 on the circulant and full classes), so they impose nothing on `N_1`
there; three-edge graphs are `O(u^3)`. The O(u^2) problem is therefore exactly: five homogeneous quadrics `q_Gamma(Delta)` on the kernel
`K_4` of the four O(u) equations (dim 71 full, 11 circulant). On `K_4` their ranks are 43, 10, 62, 34, 12 (graphs 3, 6, 7, 8, 17); on
`K_4^circ` 7, 2, 10, 6, 2; the five forms are linearly independent and the Hilbert function of their ideal agrees with that of a complete
intersection in degrees 3 (both classes) and 4 (circulant) - no linear syzygies, only the Koszul quadratic ones; this is NOT a proof of the
dimension, which is settled below by the explicit decomposition. The common zero locus has dimension >= 5 (circulant) / >= 65 (full) in any
case, so non-tautological complex solutions exist.

**Blindness to zhat^4, zhat^5, zhat^6.** A two-edge graph at genus 3 has `Sum hp <= 4` (graph 3: genus-1 vertex of dimension 4; 6: 0 + 4; 7:
1 + 1; 8: 1 + 1; 17: 1 + 0 + 1, with psi-deficit only lowering `Sum hp`), so it cannot carry two edge factors of psi-order >= 4 (each needs
`a + a' >= 3`): the forms vanish identically on the high-order coordinates of `N_1` (`q_HH = 0`, checked on the 7 high circulant and the 35
high full coordinates; the 10 coordinates of order 6 do not enter at all) and are at most linear in them through cross terms with the low
part. Writing `Delta = Delta_L + Delta_H` (orders 1-3 / 4-6), the genus-2 O(u) equations act on `Delta_L` only (orders <= 3 at genus 2) and
cut the 6-dimensional space `L_6` out of the 8 low circulant coordinates (`L_38` of the 40 low coordinates in the full class), while the
genus-3 O(u) equations are solved by `Delta_H` in an affine space of dimension 5 (circulant) / 33 (full); for fixed `Delta_L` the five O(u^2)
conditions are affine-linear in `Delta_H`. In particular `N_1^taut + Delta_H` is a solution for every admissible pure-high `Delta_H`: the
two-edge O(u^2) conditions do not see the `zhat^{>= 4}` part of `N_1` at all.

**Circulant class: the two factorisations.** For fixed `Delta_L` in `L_6` the 5 x 5 system for `Delta_H` has rank 4 off the hyperplane
`{c_1(0) = 0}` and rank 3 on it (both primes); its right kernel is the same vector throughout: the order-6 direction, invisible to all five
forms. Consistency is the vanishing of `F(Delta_L) = lambda(Delta_L) . q(Delta_L)` with `lambda` the cofactor left kernel (degree 4), a sextic
(degree exactly 6 along every random line), and the referees found that it factors: **`F = c_1(0) . m_1 . F_4`**. The two linear factors come
from the two rank-2 forms: `q_6 = c_1(0) . ell(Delta)` (the loop at the trivalent genus-0 vertex of graph 6 sees only the scalar zhat^1
coordinate, cf. the `b_1 = -3/25` remark of s.6) and `q_17 = m_1(Delta_L) . m_2(Delta)` (chain g1-g1-g1), where `m_1` is a linear form in the
low circulant coordinates `c_1(1), c_1(2), c_2(1), c_2(2)` only (no `c_1(0)`, no order 3). The factor `c_1(0)` is spurious: on `{c_1(0) = 0}`
the system drops to rank 3 and is generically inconsistent (its consistency locus is a codimension-2 piece with 2-dimensional fibres). The
factor `m_1` is genuine: on the hyperplane `{m_1 = 0}` of `L_6` the equation of graph 17 holds identically for every `Delta_H`, and the other
four have rank 4 on the four visible `zhat^{4,5}` directions and are consistent at every point with `c_1(0) != 0`, so **every point of
`{m_1 = 0}` extends to an O(u^2) solution**, uniquely up to the invisible order-6 direction where `c_1(0) != 0` and with a 2-dimensional
fibre on `c_1(0) = 0` (there `q_6` vanishes too; this is the codimension-2 piece just mentioned, over which 16 of the 17 isolated O(u^3)
rays of s.6d lie). `m_1` is defined over Q(zeta): it is the unique linear factor of the
Q(zeta)-rational form `q_17` on `K_4^circ` that vanishes on the Q(zeta)-rational pure-high subspace, hence Galois-stable up to scale. Therefore
**Q(zeta)-rational, Z/5-covariant, non-tautological solutions of the genus-2 + genus-3 O(u) equations and of all five genus-3 two-edge O(u^2)
conditions exist, in a 4-parameter projective family (the hyperplane `{m_1 = 0}` in `P(L_6) = P^5`) plus the order-6 line**. Off the two
hyperplanes the residual condition is the quartic `F_4 = 0` (irreducible over F_p at both primes on random lines), whose Q(zeta)-points are
not decided here. So the solution set in the circulant class is the union of the rational piece over `{m_1 = 0}` and the piece over
`{F_4 = 0}`, each of projective dimension 4 + 1 = 5 (plus the codimension-2 piece over `{c_1(0) = 0}`); equivalently, since `q_6` and `q_17`
split, it is the union over the choices of a factor of each of the loci `{factor of q_6 = 0, factor of q_17 = 0, q_3 = q_7 = q_8 = 0}` - three
quadrics in a P^8 each.

**An explicit rational solution** (`genus3_quadsolve3.py`, both primes): free low coordinates `(c_2(1), c_2(2), c_3(0), c_3(1), c_3(2)) = (1, 2,
-1, 3, 1)`, `c_1(0), c_1(1), c_1(2)` fixed by the two genus-2 rows and `m_1 = 0` (3 x 3 invertible), the three visible order-4, 5 coordinates
fixed by the equations of graphs 3, 6, 7, 8 (rank 4 on the five fibre directions, graph 17 automatic), order 6 normalised by `c_6(2) = 0`: all
five quadrics vanish and the four O(u) equations hold, and an independent full pass of the hatted theory (17 graphs with E <= 3, 21-node
vertex interpolation) confirms at both primes that the one-edge graphs vanish through `u^1` and all five two-edge graphs through `u^2`
(the three-edge graphs have non-zero `u^3` terms). Its coordinates are elements of Q(zeta) by construction (reduced mod both
primes in the state files; the two primes together do not suffice to reconstruct them, their denominators being large).

**Full class.** On `L_38` the 5 x 33 cross-term system has rank 5 for random `Delta_L` (both primes; effectively 5 x 24 since the order-6
block is invisible), so every generic low-order perturbation compatible with the genus-2 O(u) equations extends to an O(u^2) solution by a
choice of `Delta_H`: the main component of the solution variety is rational over Q(zeta) (an affine fibre over an open subset of the
38-dimensional `L_38`); one extended random `N_1` was verified by an independent full pass at both primes (one-edge graphs `O(u^2)`,
two-edge graphs `O(u^3)`). Components over the rank-drop locus and the pure-high family are not described.

**What this means.** At genus 3 the O(u^2) regularity of the two-edge graphs is cheap: the `zhat^{>= 4}` components of `N_1`, invisible to
every genus-3 condition quadratically, absorb all but one of the five conditions in the covariant class and all of them in the full class;
and the one remaining covariant condition is itself reducible, with a Q(zeta)-rational hyperplane of solutions coming from the factorisation
of the chain-graph quadric. So, unlike the scalar classes, the Z/5-covariant colour-mixing class does carry q-dependent corrections of
`M_0` that regularise the genus-2 edge graphs at O(u), the genus-3 one-edge graphs at O(u) and the two-edge graphs at O(u^2), with rational
coefficients. The O(u^3) conditions are treated in s.6d: the explicit solutions above do NOT satisfy the three-edge cubics, `m_1 = 0` becomes
necessary, and the 4-parameter family is cut down to the solution set described there. Genus 4, where the order-4, 5, 6 blocks first become
fully visible, is not treated.

## 6d. O(u^3): the three-edge cubics, N_2, N_3, and covariant factorisations that regularise every edge graph

(2026-11-10. Code: `fq3_cubic.py`, `genus3_cubic.py` (cubics); `fq3_circsys.py`, `genus3_cubicb.py`, `genus3_cubicc.py`, `genus3_cubicline.py`,
`genus3_pi.py`, `genus3_pi_meaning.py`, `genus3_qzeta_rec.py` (solution set of the N_1-conditions); `genus3_u3.py` (N_2, N_3, passes).
Data `g3_cubic_89.json`, `g3_cubic_61.json`, `g3_u3_<p>_<cls>.json`, `g3_qzeta_rec.json`; logs `logs/genus3_cubic_89.log`,
`logs/genus3_cubic_61.log`, `logs/genus3_u3_<p>_<cls>.log`, `logs/genus3_qzeta_rec.log`; superseded runs in `logs/old/`. Frame (0, 0), right
placement, both primes throughout. Two AI referee rounds with their own code: `g3_referee_u3.json` (direct evaluation, elimination, theory;
scripts in `referee_u3/`) and `g3_referee_u3b.json` (the write-up and the stored pass tables). The first round supplied the mechanism behind
the plane Pi, the two-dimensional elimination that completes the chart `c_1(0) != 0`, and the solution `Delta*`, which my line test had
missed; several analysis scripts had their docstrings and messages edited after the logs of record were written.)

**The complete system.** A genus-3 graph term is `u^-4 G_Gamma`, regular iff `[u^0..u^3] G_Gamma = 0`, and with `M = M_0 exp(u N_1 + u^2 N_2 +
u^3 N_3)` an E-edge graph is `O(u^E)`. Writing `Delta = N_1 - N_1^taut`, the conditions are

| unknown | conditions | type |
|---|---|---|
| `N_1` | (L2) one-edge graphs of genus 2 at `u^1` (2); (L3) one-edge graphs of genus 3 at `u^1` (2) | linear, homogeneous in `Delta` |
| `N_1` | (Q) the five two-edge graphs at `u^2` | quadrics `q_Gamma(Delta)` (s.6c) |
| `N_1` | (C) the nine three-edge graphs 4, 9, 10, 11, 12, 18, 19, 20, 32 at `u^3` | cubics `c_Gamma(Delta)` |
| `N_2` | one-edge graphs at `u^2`: `J . n_2 = -[u^2] G(N_1, 0)` (2); two-edge graphs at `u^3`: `grad Q_Gamma(N_1) . n_2 = -[u^3] G_Gamma(N_1, 0)` (5) | linear in `N_2`; solvability is a further condition on `N_1` where the matrix drops rank |
| `N_3` | one-edge graphs at `u^3`: `J . n_3 = -[u^3] G(N_1, N_2, 0)` (2) | linear, rank 2 |

and nothing else: graphs with four or more edges are regular for any `N_k`, and genus 2 needs only (L2) (it does not see `N_2`, `N_3`).
Derivation: the hatted rows are `Psi(u) R~(zhat; u)` with `Psi(u)` orthogonal for every u (it drops out of every edge), so `R' = Psi(u) (1 -
u D + O(u^2))` with `D = M_0 Delta M_0^{-1}`; edge factors are `u E'_1 + O(u^2)` with `E'_1(m, n) = (-1)^{n+1} D_{m+n+1}`, leaves are `O(u)`,
the vertex tensors are those of `M_0` up to `O(u)`, `xhat = 1 + O(u^3)`. Three edges already give `u^3`, so `c_Gamma` is the product of
three `E'_1` with the vertex tensors at `M_0`: homogeneous in `Delta`, independent of `N_2, N_3`. The `u^2` part of an edge factor is affine
in `N_2` through the same linear map as `E'_1` in `Delta`, so `[u^3]` of a two-edge graph is affine in `N_2` with linear part the polar form
of its quadric; the `N_2`-independent remainder (vertex u-dependence, a leaf, the quadratic part of `E'_2`) is not given by a formula but
read off a pass. (Checked by a referee at a generic non-solution point with true u-dependent vertices at 2^61-1, all graphs with at most
three edges except graph 4; the edge-table identities at both primes; graph 4 only through check (d) below and the passes; the affine
responses `A2 - A1 = J . n_2`, `grad Q . n_2` and `A3 - A2 = J . n_3` were confirmed on the stored pass tables by the second round.)

*The N_2 system.* Its 7 x 15 matrix (circulant class) has rank 7 at pure-high points and at the two special points with `b = zeta^2`
below, rank 6 at generic points of Pi, at the pass representative on `c_1(0) = 0`, at `Delta*` and at the special point with `b = -zeta^2`,
and rank 5 at the point of Pi with `c_1(0) = 0`. The drop is the vanishing of the whole row `grad q_17` (`q_17 = m_1 m_2`, so on `m_1 = 0`
the row is `m_2 grad m_1`), and of `grad q_6` as well at that point of Pi. Then `N_2` exists only if `[u^3] G_17(N_1, 0) = 0` (and `[u^3]
G_6(N_1, 0) = 0`). This held at every point tested, at both primes: the four pass representatives, random points of Pi, the special
points, and even random circulant `Delta` in the kernel of that row which solve nothing else (referee, graph 17 alone with true
u-dependent vertices). So consistency looks structural, but it is observed, not derived.

**Weights.** A stable graph of genus 3 with E edges has `Sum_v dim Mbar_{g_v, n_v} = 6 - E`, and the zhat-order k of `D` gives an edge of
psi-degree exactly `k - 1`. Hence a two-edge graph sees `Sum (k_e - 1) <= 4` (s.6c) and a three-edge graph `Sum (k_e - 1) <= 3`: an exact
grading in the left-placed generator `D`, a filtration in the right-placed coordinates of `Delta` used by the code (`D_k = Delta_k +`
lower orders conjugated by `M_0`). The cubics involve the orders 1..4 only, the order 4 linearly and only with two order-1 factors, the
order 3 linearly; the orders 5, 6 are invisible to them. (Numerically: `logs/genus3_cubic_*.log`, checks (b), (e).)

**The cubics** were computed exactly on the circulant class by polarisation on the ten coordinates `c_1(0..2), c_2(1,2), c_3(0..2),
c_4(1,2)` (220 points, vertex tensors frozen at `M_0`, rows to `O(u^3)`; 24 s per point, graph 4 dominating). Checks: `[u^k] = 0` for `k < 3`
at every point; `c_Gamma(0) = 0`; the fitted forms against direct evaluations at random points; homogeneity (one random point); and, as an
independent test of the frozen-vertex shortcut, check (d): the `[u^3]` coefficients of the three stored full passes of s.6c (21-node vertex
interpolation) at the three O(u^2) solutions constructed there; all nine graphs agree, at both primes. All eliminations below start from
these fitted forms and the quadrics of s.6c. Structure found (the perfect cube was confirmed as a polynomial identity by one referee and
re-checked by another from the stored forms; the vanishing on `m_1 = 0` and the `c_1(0)`-divisibility by two referees):

- `c_32 = kappa m_1^3` with `kappa != 0`, `m_1` the linear form of s.6c. Graph 32 is the star (a trivalent genus-0 centre with three genus-1
  leaves); `m_1` is the edge-dressed genus-1 one-point function, `m_1 ~ chat^D_1(3) - (5/2) chat^D_2(4)` in the Fourier coordinates of `D`
  (`<H>_1 = 1/12` with psi^0, `<tau_1(1)>_1 = 5/24` with psi^1), a single linear form by covariance. So **`m_1 = 0` is forced**: of the pieces
  of s.6c only those over `m_1 = 0` survive (the hyperplane and, inside it, the locus `c_1(0) = 0` with its larger fibres); the quartic
  piece `F_4 = 0` is eliminated. `c_18` and `c_20` are multiples of `m_1`; `c_11` and `c_19` are multiples of `c_1(0)` (loop on a trivalent
  genus-0 vertex, as for `q_6`).
- On `m_1 = 0` the three-edge conditions are the cubic of graph 4 (no order-4 coordinates) and those of graphs 9, 10, 11, 12, 19 (affine in
  the order-4 coordinates).

**The solution set of the N_1-conditions in the circulant class.** All conditions are homogeneous in `Delta`, so the solution set is a cone:
a "point" below is a ray `N_1^taut + s Delta`. It is described by the order-1 part `c_1 = (c_1(0), c_1(1), c_1(2))`; "low part" means the
orders 1..3. Each low part carries an affine fibre of orders 4..6: one combination of order-6 coordinates in the kernel of the genus-3 rows
is invisible to everything at genus 3; the fibre is that line at generic points of Pi and at `Delta*`, and 2-dimensional at the solutions
with `c_1(0) = 0`.

| stratum | low parts of the solutions, up to scale | how established |
|---|---|---|
| `Delta_L = 0` | the pure-high family: the orders 4, 5, 6 subject to the two genus-3 rows (5 parameters); all quadrics and cubics vanish identically there by weights | exact |
| `c_1 = 0` | nothing else: on `c_1 = 0` the cubics involve `c_2` only and those of graphs 9, 18, 32 survive; with `m_1 = 0` only `kappa' c_2(1)^3` remains, so `c_2 = 0`; then the only non-zero quadric is `q_7 ~ c_3(2)^2` (value `(5/192) zeta` on the line left by the genus-2 rows), so `c_3 = 0` | exact supports (two referees; the test in `genus3_cubicc.py` is single-point), both primes |
| `c_1(0) = 0`, `c_1 != 0` | 17 points. The cubics 11, 12, 18, 19, 20, 32 and the quadrics 6, 17 vanish identically; the cubics of graphs 4, 10 are conditions on the low part; the order-4 coordinates are constrained by `c_9, q_3, q_8` only (one 3 x 3 determinant); the genus-3 rows and `q_7` are always solvable in the orders 5, 6: three conditions (degrees 3, 3, 5) on the P^3 with coordinates `(c_1(1) : c_1(2) : c_2(1) : c_3(2))`. Result: 13 points at generic `b = c_1(2)/c_1(1)`, a pair at `b = zeta^2`, two at `b = -zeta^2`, of which one is the point of Pi with `c_1(0) = 0`. None on `a_4(b) = 0`, none at `c_1 = (0, 0, 1)` | resultant `U(b)` (degree 51 = `a_4^14` times degree 23; squarefree part 18; genuine part 13 with factor degrees (1, 2, 10) at 2^61-1 and (1, 4, 8) at 2^89-1); recomputed by a referee with three other eliminations; all 13 + 4 points lifted to full solutions over their root fields |
| `c_1(0) != 0` | the plane Pi (below) and ONE isolated point `Delta*` | line tests (only Pi off the conic `a_4 = 0`; they do not examine the fibre over the conic) and one referee's two-dimensional elimination of the unscaled system (resultants of degree 1586, common factor of degree 61, squarefree 41: the point over which the `c_1`-line of Pi touches the conic `a_4 = 0`, and 40 points at which the cubics and `q_3` are solvable but the full system is not; two seeds, both primes). `Delta*` sits over the tangency point |

So the circulant solutions with non-zero low part are the plane Pi and 17 isolated rays: 16 with `c_1(0) = 0` and `Delta*`. Arithmetic of the
isolated points, as far as two primes show: the 13 are either one Galois orbit or a Q(zeta)-rational point and an orbit of 12 (the factor
degrees do not decide; the F_p-rational one does not reconstruct in Q(zeta) with small height, which favours a single orbit); the pair at
`b = zeta^2` is conjugate over a quadratic extension of Q(zeta) that splits at both primes (sum and product of their coordinates
reconstruct, the points do not); the point at `b = -zeta^2` outside Pi is Q(zeta)-rational; `Delta*` is Q(zeta)-rational with small
coordinates (through order 4: denominators `2^a 5^b`, `c_1 = (1, (zeta + zeta^3)/2, (zeta + zeta^2)/2)`, `c_3(0) = -18/25` where the point of
Pi with the same `c_1` has 0). `a_4(c_1)` is the coefficient of `c_3(2)` in the cubic of graph 4 after `c_3(0), c_3(1)` are eliminated by the
genus-2 rows; all eliminations use only homogeneity (no denominators).

**The plane Pi.** In the chart `c_1(0) != 0` the solutions found on random lines span a 2-dimensional LINEAR space `Pi` of `Delta` on which the
four linear rows, the five quadrics and the nine cubics all vanish identically (also confirmed by direct graph evaluation at several points,
both primes, by two referees). Its coordinates of orders 1..4 are small elements of Q(zeta) with denominators `2^a 5^b`
(`genus3_qzeta_rec.py`: lattice reconstruction from the two primes, e.g. `A`: `c_1 = (1, 0, -zeta^4/2)`, `c_2(1) = (-1 - zeta + zeta^2 +
zeta^3)/10`, `c_3(0) = (18 + 6 zeta^2 + 6 zeta^3)/125`); the orders 5, 6 follow from the closed form (`genus3_pi_meaning.py`, both primes;
found independently by the theory referee):

    Pi = { Delta = M_0^{-1} (zhat D_1) M_0 :  D_1 constant, symmetric, Z/5-covariant, with no H (x) H component }   (+ the invisible direction),

i.e. `M = M_0 exp(u Delta + ...) = exp(u zhat D_1) M_0 (...)`: a constant zhat-linear dressing on the LEFT of `M_0`, `R' = Psi(u) (1 - u zhat
D_1) + O(u^2)`, and the first-order edge factor is the constant, psi-free propagator `D_1`. In the normalised canonical basis `phi_m =
(zeta^{m j})_j` is proportional to the flat class `H^{m-1}` (unit `phi_1`, pairing `m <-> -m`), and a covariant symmetric `D_1` is

    D_1 = alpha H (x) H + beta (1 (x) H^2 + H^2 (x) 1) + gamma (H^3 (x) H^4 + H^4 (x) H^3);     Pi is alpha = 0.

*Mechanism* (theory referee; the two inputs were checked numerically at both primes). For `D = zhat D_1` the `u^E` coefficient of an E-edge
graph is a contraction of PRIMARY (psi-free, top-degree) correlators of the frozen vertex theory `M_0 . omega^` with a propagator whose
ends lie in `W = span(1, H^2, H^3, H^4)`. At genus >= 1 the only non-zero primary correlators of that theory are `<H, ..., H>_g` (computed
for `(g, n) = (1, 1..4), (2, 1..3)`: 1/12, 1/12, 1/6, 1/2; 1/600, 1/200, 1/50, the derivatives of `-(1/12) log t` and `(1/1200) t^-2`: the
primary sector is the one-variable conifold free energy in the direction `H`), and `(0, 6)` has only `<H^6>` and `<H, (H^2)^5>`. Every
graph of genus 3 with at most three edges (and every one-edge graph of genus 2) has a vertex of genus >= 1 carrying a half-edge or is the
single `(0, 6)` vertex, so every term vanishes, linearly in `D_1`: this is why a linear space of solutions exists although 18 conditions
on 15 unknowns have negative expected dimension. Tested by the referee: in the FULL class every symmetric `D_1` with `D_1 phi_3 = 0` (10
parameters, propagator in `Sym^2 W`) solves all `N_1`-conditions (two random points per prime, direct evaluation); `H (x) x + x (x) H` with
`x` in `W` fails exactly the graphs 6, 9, 10, 18, 20, 32; a covariant `D_1` with `alpha != 0` fails all 18 conditions. The direction
`Delta = N_1^taut` (rescaling the tautological correction; `Delta = -N_1^taut` is `M_0` alone) is not in Pi: its `H (x) H` coefficient is
`-1/5`.

The two directions of Pi differ in nature. `beta` (string direction): `D_beta^2 = 0`, `D_beta 1 = 0`, and for the exact dressing `M = exp(u
zhat D_beta) R~(zhat; u)` one has `R' = Psi(u)(1 - u zhat D_beta)` exactly; every edge carries the unit at a psi-free end and there are no
leaves, so by the string equation, and `eta(x, y) = 0` on `{1, H^2}` at trivalent genus-0 vertices, every edge graph should vanish at
every genus: a second, trivially equivalent, tautological factorisation (the referee's argument; checked at 2^61-1 on the 16 graphs with at
most three edges through `u^3`). `gamma`: only the leading order is explained. The exact left dressing `exp(u zhat D_1) R~(zhat; u)` at a
random point of Pi does make all 16 graphs with at most three edges `O(u^4)` (2^61-1: all 16; 2^89-1: the 14 other than graphs 3, 4), but
it is not a symmetry: in the frozen theory the sums for `exp(eps zhat D_1) M_0` are non-zero at `eps^4..eps^6` for graphs 2, 3, 7, 10.
Nothing supports the `gamma` direction beyond the orders needed at genus 3.

**N_2, N_3 and the verification passes** (`genus3_u3.py`). Per case and prime, three successive passes with the 21-node vertex interpolation:
`(N_1)` and `(N_1, N_2)` on the 17 graphs with at most three edges, used to solve `N_2` and then `N_3`, and `(N_1, N_2, N_3)` on all 42
graphs. `N_2`, `N_3` are solved in the circulant class with a fixed combination of (the first eight) free directions added; they differ from
the tautological `N_2`, `N_3`, except that in the pure-high case the orders 1, 2 of `N_2` are forced to be the tautological ones; `N_3` was
given no component in the orders 5, 6. Four representatives were run:

| N_1 | prime | P1: cubics hold | N_2 system | P2 | P3: all 41 edge graphs O(u^4) | total | vertex polar part | genus 2 | script verdict |
|---|---|---|---|---|---|---|---|---|---|
| pure-high: `Delta` in the orders 4, 5, 6 only (coordinates (1, 2, -1, 3, 1) on the 5-dimensional kernel) | 2^89-1 | yes | rank 7, consistent, 8 free parameters | yes | yes | yes | yes | yes | all checks passed |
|  | 2^61-1 | yes | rank 7, consistent, 8 free parameters | yes | yes | yes | yes | yes | all checks passed |
| the point `A + 2 B` of Pi (`c_1 = (1, 2, -zeta^4/2 - 2 zeta^2)`) | 2^89-1 | yes | rank 6, consistent, 9 free parameters | yes | yes | yes | yes | yes | all checks passed |
|  | 2^61-1 | yes | rank 6, consistent, 9 free parameters | yes | yes | yes | yes | yes | all checks passed |
| the isolated solution on `c_1(0) = 0` that is F_p-rational (`c_1 = (0, 1, b*)`) | 2^89-1 | yes | rank 6, consistent, 9 free parameters | yes | yes | yes | yes | yes | all checks passed |
|  | 2^61-1 | yes | rank 6, consistent, 9 free parameters | yes | yes | yes | yes | yes | all checks passed |
| the isolated solution `Delta*` with `c_1(0) != 0` | 2^89-1 | yes | rank 6, consistent, 9 free parameters | yes | yes | yes | yes | yes | all checks passed |
|  | 2^61-1 | yes | rank 6, consistent, 9 free parameters | yes | yes | yes | yes | yes | all checks passed |

In all four cases, at both primes: after the first pass the one-edge graphs are `O(u^2)`, the two-edge graphs `O(u^3)` and the three-edge
graphs `O(u^4)` (the nine cubic conditions in the full hatted theory); in the last pass **all 41 edge graphs of genus 3 are `O(u^4)`, i.e.
regular at the conifold point**, and the total over the 42 graphs equals the point-frame total through `u^5`. Consequently the vertex graph
has `[u^0..u^3]` equal to those of the total and carries the whole polar part `(1/25200) t_c^-4 + 0 t_c^-3 + 0 t_c^-2 + 0 t_c^-1`. The six
genus-2 edge graphs are `O(u^2)` with the same `M` (a condition on `N_1` only). The second referee round reproduced every entry of the
table from the stored pass tables with its own code. Not run: other points of Pi, the three special points on `c_1(0) = 0`, and the other
twelve of the 13 (for these only the solvability of the `N_2` system was checked, by a referee, without a 42-graph pass). The pure-high,
Pi and `Delta*` representatives have Q(zeta)-rational `N_1`; the third is the F_p-rational root of the degree-13 factor at each prime, so
its two runs need not be reductions of the same point. `N_2`, `N_3` are prime-dependent choices: the two primes agree in structure (ranks,
zero patterns) and in the entries that reconstruct rationally (9 per class, 41 for `Delta*`), not coefficient by coefficient.

**Full class.** Not classified. Circulant solutions are full-class solutions and the 10-parameter family `D_1 phi_3 = 0` above solves the
`N_1`-conditions; no full-class pass with `N_2, N_3` was run. In the one direct construction tried (orders 1, 2 on the linear locus of
graph 32, order 3 from linear equations; `genus3_u3.py --cls full`, `fq3_fullsol.py`) the 12 linear equations for the orders 4..6 have rank
11 and are inconsistent at both primes; the attempt was abandoned and says nothing about other full-class points. The count of s.6c
("every generic low part extends") does not carry over, since the cubics of graphs 4, 18, 20, 32 do not see the order 4.

**What this means.** (i) At genus 3 the gap can be put entirely into the vertex graph by Z/5-covariant, q-dependent factorisations other than
the tautological one, with every one of the 41 edge graphs regular: verified by a 42-graph pass at both primes for four representatives
(pure-high, a point of Pi, one of the 13 points on `c_1(0) = 0`, `Delta*`). (ii) The `N_1`-conditions leave, up to the scaling `Delta -> s
Delta` and the affine fibres in the orders 4..6: the pure-high family (orders 4, 5, 6 subject to the two genus-3 rows); the plane Pi, which
is kinematic (left dressing by constant propagators avoiding `H`, the only direction in which the conifold-limit vertex theory has
higher-genus primary correlators; its string direction is, by the argument above, an exact gauge direction, its other direction is not
known to survive genus 4); and 17 isolated rays. The excluded propagator `H (x) H` is the one the conditions see. (iii) `N_3` is
constrained by 2 equations out of 15 covariant parameters; `N_2` by 7 equations of rank 7, 6 or 5 depending on the solution, with one
(or two) solvability conditions on `N_1` where the rank drops, satisfied wherever tested. (iv) So genus 3 does not single out the
tautological factorisation; but the only positive-dimensional covariant families of first-order corrections are the kinematic plane and
the pure-high family. The scalar (Hodge-type) class contains none of this (s.6, s.6b). Genus 4 constrains the orders 4, 5, 6 and `N_2`, and
would test the `gamma` direction and the isolated rays.

## 6e. The natural factorisation: `M = M_0(tau(u) zhat)` (2026-11-11; full account in NATURAL_FORM.md)

The solution set of s.6d contains a distinguished, parameter-free point. With `1/tau(u) = 2F1(1, 4/5; 9/5; u) = Sum_n 4 u^n/(5n+4)`, i.e.
`tau(u) zhat = -(5/4) zeta^i z/(u_i - u_i^c)` (`z` over the canonical coordinate measured from the conifold point), the Teleman form factorises as
`R~(zhat; u) = G(u, zhat) M_0(tau(u) zhat)` with `[u^n] G` a polynomial in `zhat` (degree 1, 3, 4, 5 for n = 1..4 in frame (0, 0); `G_1 = zhat X_1`
with `X_1` in the span of the two propagators of Pi). `genus3_natural.py` (logs `logs/genus3_natural_{89,61}_nat.log`,
`logs/genus3_natural_89_nat_ab1o3_m2o7.log`): with `M = M_0(tau zhat)` mod `(zhat^7, u^4)` and no adjustment, all 41 edge graphs of genus 3 are
`O(u^4)`, the vertex graph is `(1/25200) t_c^-4 + 0 + 0 + 0`, and the six genus-2 edge graphs are `O(u^2)` (frame (0, 0) at both primes; frame
`(1/3, -2/7)` at 2^89-1). In the terms of s.6d: `N_1 = N_1^taut - M_0^-1 (zhat X_1) M_0` (a point of Pi), and the canonical `N_2`, `N_3` satisfy
the seven and two linear conditions with no free parameter. In the family `M_0((1 + a u) zhat)` each of the four O(u) equations of s.4
separately gives `a = -4/9` (`natural_form/extra/nat_rescaling_forced.log`). `M_0` itself is the R-matrix of the Landau-Ginzburg model
`x^4 + v/x`, primitive form `dx` (`natural_form/lg/LG_MODEL.md`), and `G` is the hatted form of a gauge regular at the conifold point,
computable from the periods (`natural_form/theory/UNIFORM_REDUCTION.md`). Status, caveats and what is conjectural: NATURAL_FORM.md s.1, s.3.

## 7. Frame dependence

Steps 1-4 were repeated at the conifold frame `varpi_0^cf = f_0 + alpha g + beta f_2` with `(alpha, beta) = (1/3, -2/7)` (log
`logs/genus3_ou_89_ab1o3_m2o7.log`, state `g3_state_89_ab1o3_m2o7.json`; there `t_c = u + (11/30) u^2 + (751/1575) u^3 + ...`). Frame
independent: the symbol, `M_0`, both Jacobians `J` (genus 3) and `J2` (genus 2), `g1` of the one-edge graphs (`-1463/600000`,
`232351/120000000`), hence the whole O(u) system of s.4 (identical matrices, right-hand sides, ranks and solutions in every class), and
the leading coefficient `[u^E]` of the edge graphs 2, 5, 3, 17, 4, 32 at `M_0`. Frame dependent: the point-frame total (`1/25200 - 11
u/189000 - 197 u^2/8820000 + ...`), whose `t_c` expansion is nevertheless again `(1/25200) t_c^{-4} + 0 + 0 + 0 - 2567737/8400000000 +
...` (the gap is frame independent, as it must be), the `u^{>= 1}` coefficients of the vertex graph (`[u^1]`: `140747/360000000` vs
`3355687/7560000000`), and - unlike genus 2 - the leading coefficient `[u^E]` of the other 30 edge graphs with E <= 5 at `M_0` (e.g. graph 6:
`-28329/50000000` vs `9443/150000000`, graph 7: `35267/22500000` vs `4993/5625000`): the per-graph polar parts with `M_0` are frame dependent at
genus 3, so the O(u^2) obstruction of s.6 is a statement at frame (0, 0) (the O(u) system is not). The five six-edge graphs are `O(u^6)`
and carry no information at this precision. This frame was run at 2^89-1 only (its structural data coincide with those of frame (0, 0)). The two FAIL lines of
that log are the genus-2 reference checks, which compare against the `(0, 0)`-frame genus-2 tables and are not meaningful there.

## 8. Primes and reconstruction

All steps were run at 2^89-1 (values of record) and 2^61-1. `genus3_cross.py` (final run, after steps 5, 5b and the
tautological control pass at both primes): 494 values reconstructed at 2^89-1 with the default bound 1e10 agree with the 2^61-1 residues (0
mismatches, 0 missing); the 1e13 fallback reconstructions printed in the logs are reliable only when confirmed there (164 confirmed, 345 spurious: point-frame per-graph coefficients with denominators beyond 1e13 are given by their residues in
the state files); ranks, consistency, dimensions, regular sets and the Jacobian degree are identical at both primes. At 2^89-1 the 1e13
fallback has spurious rate `B^2/p ~ 0.16`, so an unconfirmed fallback value in a log (printed with a leading `~` by the final `fq3_rec.Rec`;
older log lines print them bare, e.g. the "roots" of graphs 7, 8 at 2^89-1) carries no evidence; the two residues together determine any
rational with `|num|, den < 1e22` by CRT (used for the t^0 coefficient of graph 17).

## 9. Caveats

Genus 3 only; exact but finite; step 5 (s.6) is right placement, frame (0, 0), two primes; s.6b is frame (0, 0), two primes; the second frame is
single-prime; the log lines `PASS g3:...` / `PASS g2+g3:...` of step 4 only record the ranks (they cannot fail), and the step-3 degree
PASS is tautological (7 nodes for degree <= 6). All other checks are internal consistency checks (M-independence, group law, frozen = `[u^0]`, genus-2
reproduction, Hodge integrals, string/dilaton) plus AI referee rounds with fresh code (two for s.1-7, one for s.6c, two for s.6d); no human has checked the code. The symplectic
conventions (`M(zhat)^T M(-zhat) = 1` in the normalised basis, `R^{-1}(z) = R(-z)^*`) are those of the genus-2 notes. The circulant and
full classes are analysed through O(u^2) (s.6c): the rational hyperplane of solutions is established at both primes with an explicit solution
verified by full passes; the exact dimension of all components (the Hilbert-function agreement is not a dimension proof) is open.
The O(u^3) analysis (s.6d) is for the circulant class, frame (0, 0), right placement, modulo the two primes. The classification of the
`N_1`-solutions rests on exact supports (`c_1 = 0`), on root counts of resultants with all points lifted over their root fields by a referee
(`c_1(0) = 0`), and, for `c_1(0) != 0` (completeness and `Delta*`), on a single referee's two-dimensional elimination with random
combinations and shears (two seeds); the project's own line tests discard the conic `a_4 = 0`, which is where `Delta*` sits. All
eliminations use the fitted cubic and quadric forms (validated at random points and, at the solution points, by direct graph
evaluation). None of it is a characteristic-0 proof and none was checked by a human. Whether the 13 points on `c_1(0) = 0` form one
Galois orbit is not decided by the two primes. The solvability of the `N_2` system where its rank drops is observed, not derived. The
mechanism behind the plane Pi uses primary correlators of the frozen vertex theory computed only for the (g, n) needed at genus 3; its
string direction rests on an argument checked at one prime through `u^3`, its `gamma` direction is explained at leading order only.
The full class was not classified at this order. Of the explicit solutions, those of s.6c are given by their F_p residues and their
Q(zeta)-rational defining data; Pi and `Delta*` are reconstructed in Q(zeta) through order 4; the third pass representative of s.6d is
known only by its residues.
