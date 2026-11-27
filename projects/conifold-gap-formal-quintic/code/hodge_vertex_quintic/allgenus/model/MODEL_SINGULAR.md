# The singular structure of the model vertex theory, and a reduction of Conjecture N (2026-11-13)

Worker key `model`, folder `allgenus/model/`. AI-derived and AI-checked only; no human has read it. Status labels: PROVED (complete
argument written here, or a quoted theorem with arXiv number), VERIFIED (exact computation; orders, field, frames stated),
CONJECTURAL, REFUTED. Conventions: NATURAL_FORM.md, UNIFORM_REDUCTION.md, LG_MODEL.md; `H^a` the flat classes (`H = d/dv`), `tau_k(H^a)`
an ancestor insertion `psi^k H^a`, `<..>_h` the genus-h ANCESTOR correlator (integral over Mbar_{h,n}) of the frozen theory
`V = M_0.omega^` = the LG model `x^4 + v/x` (form dx) at the point v; `F_h = c_h v^(2-2h)` (h >= 2), `F_1 = -(1/12) log v`,
`F_0 = -(5/2) v^2 log v` (n >= 3 derivatives), `c_h = (-1)^(h-1) B_2h/(2h(2h-2) 5^(h-1))`.

## 0. Summary

1. **Engine (VERIFIED).** `ms_engine.py` (written by the first attempt of this task; reviewed and re-run here): flat basis, Z/5 charges,
   exact over Q, graphs from `fq3_graphs` (or the faster `ms_graphs.py`, same graphs). It reproduces F_2, F_3, F_4, `<tau_1(1)>_1 = 5/24`,
   the primary H-sector, string and dilaton equations (`ms_regress.log`) and, mod 2^89-1, every flat correlator of the project's
   independent `MVertex6` that its symbol order reaches: 2512 + 3955 + 6125 + 23125 correlators, genus 0-4, n <= 4 (`ms_regress.log`,
   `ms_regress2.log`, `ms_regress3_33.log`, `ms_regress3_24.log`).
2. **The one-point law is REFUTED as stated; it holds on K.** `<tau_k(H^b)>_h = PiC(tau_k H^b) F_h'(v)` holds for every singular
   one-point function at h = 2, 3, 4 whose insertion lies in the subspace K of item 4, and fails for all those off K:
   `<tau_5(H)>_3` (e = -1): `53/1209600 = (-53/192) F_3'` against the law's `(1/96) F_3'`; at genus 4 `tau_5(H)`, `tau_6(1)`,
   `tau_7(H^4)`, whose three remainders are `Pi_5(X) x 19/100800000` (one number: the component along the descendant t^H_5).
   Regular one-point functions (v-exponent e >= 0) do not obey the law (e.g. `<tau_3(H^3)>_2`).
3. **No single non-linear V (REFUTED).** `Hess_h(X, Y) := (<XY>_h - F_h'' PiC(X) PiC(Y))/F_h'` depends on h for the pairs
   `t1(H) t2(H^3)`, `t2(H) t1(H^3)`, `t3(H) t0(H^3)` (singular at h = 2 and h = 3), so the singular part is not `F_h(V(q))` for any
   genus-independent V, linear or not, and not the string-only extension of PiC (42 singular pair instances differ from it). All
   genus dependence sits in insertions `tau_k(H)`, k >= 1.
4. **The Model Lemma (ML) on K (VERIFIED on all data; CONJECTURAL in general).** `K := {q : Pi_i(q) = 0 for i >= 1}`,
   `Pi_i(q) := (1/5)[z^i] eta(S(-z)H^2, q(z))` = the DESCENDANT time `t^H_i` of the ancestor vector q (S the homogeneous calibration,
   `S(z)H^2` = the vanishing flat section, `Pi_0 = PiC`). K is spanned by `tau_0(H)` and `kappa(X, k) = [S(z) X z^k]_+`,
   X in {1, H^2, H^3, H^4}. (ML): for `q in -z + K` the genus-h ancestor potential is `F_h(Pi_0(q)) + Reg_h(q)` with `Reg_h` regular:
   `<H^n>_h = F_h^(n)(v)` and every correlator of K-generators with at least one `kappa` has v-exponent >= 0 or vanishes.
   Data: genus 0-4 (see s.5): 5289 such correlators, 1395 of them singular, all 1395 vanish; `<H^n>_h = F_h^(n)` in all 17 cases.
   The genus-2 four-point table (376 singular instances, all zero) was computed after (ML) had been formulated: a prediction.
5. **Lemma K (PROVED) and its prediction (VERIFIED).** For a symplectic gauge with `P Y^C = Y^C_Q` (no positive power of z) every edge
   lies in `K_0 (x) K_0` (`K_0 = K cap ker PiC`) and the leaf point `-z + T` in `-z + K` with `Pi_0 = t_c`: ALL shifted vanishing
   covectors annihilate the edges, not only `Pi_0`. On the actual `R'` tables (referee rows zhat <= 24, u < 18, both primes, frames
   (0, 0) and (1/3, -2/7)): shifts k = 0..20 vanish (20790 non-trivial sums per run; random-covector controls non-zero), and
   `Pi_k(T) = 0` for k = 1..20 (controls non-zero). Nothing was fitted.
6. **Theorem R (PROVED, conditional).** (ML) at every genus + Lemma K's hypotheses + regularity of the gauge (Theorem U) + the
   factorisation (Corollary F) => Conjecture N: every stable graph with an edge is regular at u = 0, and the vertex graph is
   `F_g(t_c)` + regular. The D-directions `tau_k(H)` (where items 2-3 fail) never occur in the graph sum.
7. **Partial proof of (ML) from the literature.** On the tower of `H^2` (the pole position a_1 of the Hurwitz space H_{0;1,4}),
   arXiv:2211.12259 Thm 2 identifies the descendant potential with the generating function of rooted hypermaps; the v-exponent of a
   correlator of that tower equals the number of vertices of the hypermap, so it is >= 1 (part 2, s.8). The towers of 1, H^3, H^4
   are open (VERIFIED in the computed range only).
8. **Answers (part 2, s.9).** (b) on K the singular part is `F_h` of the LINEAR coordinate Pi_0 (no genus-0-point non-linearity,
   no jets; the family F_h itself is the Barnes/c = 1 free energy, which satisfies the second-difference (Toda) equation in Pi_0);
   (c) the primary-sector lemma follows from (ML) (every genus >= 1 primary correlator is singular by weight); (d) string, dilaton and
   Euler preserve the singular part on K. Model lemma in edge form: B in `K_0 (x) K_0` and leaves in K give regular graphs, and this is
   exactly what (G) produces. A non-canonical z^5 normalisation of the gauge breaks Lemma K (control VERIFIED), but its first-order
   effect on the genus-4 total vanishes at the leading order (`Delta_4(5) = 0`, exact), so genus 4 is NOT a sharp test of it.

## 1. Setting

**The model.** `V = M_0.omega^` is the Givental CohFT of the LG model `x^4 + v/x`, form dx, on its Z/5 line (LG_MODEL.md; the Hurwitz
Frobenius manifold H_{0;1,4} with metric x 20 and z -> -4z). Everything is homogeneous: with `w(a, k) := (4 - c_a - 4k)/5`,
`c = (0, -1, 3, 2, 1)`, a genus-h correlator is `const * v^e`, `e = 2 - 2h - Sum_i w(a_i, k_i)`; "singular" means e < 0 and value != 0.
The correlators computed are ANCESTOR correlators (CohFT integrals over Mbar_{h,n} with its own psi classes); the dilaton
equation reads `<tau_1(1) X>_h = (2h - 2 + n)<X>_h` and the string equation lowers a psi.

**The vanishing covector and its shifts.** S(z) = 1 + S_1/z + ..., `z d/dv S = (H *) S`, integration constants fixed by
homogeneity; the column `S(z)H^2` has no logarithm (only the H-column has: `S_1 H = -log(v) H^2`), and
`Y^C_b(z) = -(1/5) eta(S(z)H^2, H^b)` is the residue period at the pole (the vanishing cycle; ms_model.py). Define for `q in H[[z]]`

        Pi_i(q) := (1/5) [z^i] eta(S(-z)H^2, q(z)),        Pi_0 = PiC,   Pi_i(tau_k(H^b)) = PiC(tau_(k-i)(H^b)).

`Pi_i(q)` is the descendant time `t^H_i` of the ancestor vector q (descendant times = `[S^-1(z) q(z)]_+`, `S^-1(z) = S^*(-z)`).
`PiC = 1` on H, `-v` on tau_1(1), `-v^2/2` on tau_2(H^4), `-v^3/6` on tau_3(H^3), `v^4/24` on tau_4(H^2), `v^4/96` on tau_5(H),
`-v^5/480` on tau_6(1), ...; it vanishes unless `b + k = 1 mod 5` (the charge selection rule of one-point functions).

**K and its basis (PROVED, elementary).** `K := {q : Pi_i(q) = 0, i >= 1}` = descendant times with `t^H_k = 0` for k >= 1;
`K_0 := K cap ker PiC`. (i) `K = C tau_0(H) + span{kappa(X, k)}`, `kappa(X, k) := [S(z) X z^k]_+ = Sum_{j <= k} S_j X z^(k-j)`,
X in {1, H^2, H^3, H^4} (the ancestor images of the descendant insertions tau_k(X)); the kappa have coefficients that are
non-negative powers of v (no logarithm, since the columns X != H of S are log-free) and span K_0. (ii) `H[[z]] = K + span{tau_k(H), k >= 1}`
(direct): `Pi_i(tau_k(H)) = PiC(tau_(k-i)(H))` is triangular with diagonal 1. (iii) A basis insertion `tau_k(H^b)` lies in K iff
b = 0, k <= 1; b = 1, k = 0; b = 2, k <= 4; b = 3, k <= 3; b = 4, k <= 2. Its D-content is `{Pi_i(tau_k H^b), i >= 1}`.

## 2. The engine and its regressions (VERIFIED)

`ms_rows.py` (symbol rows over Q from `natural_form/lg/lg_symbol_rows_Q.json`, flat basis h^a = N_a H^a, edge/leaf/leg factors in Z/5
charge variables), `ms_engine.py` (amputated amplitudes summed over all legged stable graphs, cached in the scratch folder),
`ms_graphs.py` (a faster generator of the same graphs: equal counts and equal Sum 1/|Aut| per edge level as `fq3_graphs` for ten
types; `ms_chunk.py` reproduces `amp(3,2)` exactly). Both engine files were written by the first attempt of this task; they were
reviewed line by line here (edge formula, charges, leaves, TFT selection rule) and all amplitudes were recomputed (the old cache is kept
apart), except `amp(3,3)`, which a process of the first attempt finished at 18:57 with the same code and which is checked in full
against `MVertex6` below.
- `ms_regress.log`: F_2 = 1/1200, F_3 = 1/25200, F_4 = 1/180000, `<tau_1(1)>_1 = 5/24`, `<H>_1 = -1/12`, `<H^3>_0 = -5`, the primary
  H-sector `<H^n>_h = d_v^n F_h` for (h, n) = (1, 1..4), (2, 1..3), (3, 1..2), string and dilaton equations, permutation symmetry, and
  2512 flat correlators of `fq3_mvertex.MVertex6` on the project's own symbol (fq_rows -> fq_symbol), mod 2^89-1: ALL PASS.
- `ms_regress2.log`: every one- and two-point tensor at genus 1-4 and every three-point tensor at genus 2 with psi-deficit <= 6:
  3955 flat correlators = MVertex6 mod 2^89-1.
- `ms_regress3_33.log` and `ms_sd33.log` (the table not produced by this session's runs): all 49 tensors of (3,3) with psi-deficit <= 6,
  6125 flat correlators = MVertex6 mod 2^89-1; and all (3,3) correlators with a tau_0(1) or tau_1(1) insertion (140 string, 115 dilaton)
  agree with the string and dilaton equations applied to the (3,2) table computed here; permutation symmetry holds. (A full check with
  the symbol to zhat^9 was started and stopped at the 25-minute limit; MVertex6 is slow at psi-deficit 7-9.)
- `amp(2,4)` (5608 graphs, `ms_graphs` + `ms_chunk`, 627 s, `ms_chunk_2_4.log`): all 37 tensors with psi-deficit <= 6, 23125 flat
  correlators = MVertex6 mod 2^89-1 (`ms_regress3_24.log`); string (551), dilaton (391) and permutation checks against (2,3)
  (`ms_sd24.log`). `amp(4,2)` (28037 graphs) was not computed (s.10 of part 2).

## 3. One-point functions (`ms_onept.log`)

| h | singular insertions (e < 0) | law `PiC F_h'` holds | fails |
|---|---|---|---|
| 2 | t0(H), t1(1), t2(H^4) | all 3 | - |
| 3 | t0(H), t1(1), t2(H^4), t3(H^3), t4(H^2), t5(H) | 5 (all in K) | `t5(H)`: value 53/1209600, law -1/604800 |
| 4 | t0(H) .. t4(H^2), t5(H), t6(1), t7(H^4) | 5 (all in K) | t5(H), t6(1), t7(H^4) |

At genus 4 the three remainders are `19/100800000 x (1, -1, -1/2)` = `Pi_5(X) x 19/100800000`: the remainder of every singular one-point
function is carried by its component along the descendant `t^H_5`, the only D-direction allowed by the charge (i = 0 mod 5) in range.
So the one-point law is exactly the statement (ML) restricted to one point, and it fails exactly off K. Regular one-point functions
(e >= 0: t3(H^3), t4(H^2) at genus 2; t6(1), t7(H^4) at genus 3; t8, t9, t10 at genus 4) do not obey it. At genus 1 the law holds for H
and fails for the dilaton by the anomaly `5/24 - 1/12 = 1/8` (e = 0, regular).

## 4. Two points: no single V (REFUTED), and what the "genus-0 point" gives (`ms_twopt.log`)

For F_h(V(q)) with dV = PiC, `Hess_h(X, Y) := (<XY>_h - F_h'' PiC(X) PiC(Y))/F_h'` must be one bilinear form for all h on singular pairs.
- h-dependent (both genera singular): `t1(H) t2(H^3)`: -23/12 (h = 2), -13/48 (h = 3); `t2(H) t1(H^3)`: -11/6, -31/60;
  `t3(H) t0(H^3)`: -3/2, -11/60. So no genus-independent V exists (any V, linear or not, jets of the genus-0 solution included,
  since a jet-dependent V is still genus independent). Each of the three pairs has one insertion off K.
- Against the string-only extension `V_str(q) = PiC(q(s*))` (q flowed along the string vector field to q^1_0 = 0; the minimal
  X_{-1}-invariant, dilaton-degree-1 extension of PiC; Hessian = `[X = t0(1)] Pi_1(Y) + (X <-> Y)`): 42 singular instances differ.
- The ancestor genus-0 point `W = w^H / w^1_x` (string and dilaton covariant, degree 1) has `dW = 1` on H, `-v` on tau_1(1) and 0 on
  every other insertion; it is already refuted by `<tau_2(H^4)>_2 = c_2/v != 0` (s.3). The descendant genus-0 point `v^H(t)` has
  dilaton degree 0 and cannot carry `F_h` (degree 2 - 2h) at all.
- On K the correct "V" is the LINEAR function `Pi_0` (the descendant coordinate t^H_0, i.e. the vanishing period of q), with no jet and
  no genus correction (s.5). The genus dependence lives entirely in the D-directions tau_k(H), k >= 1 (complete list of singular
  correlators with remainders: `ms_singlist_full.log`).

## 5. The Model Lemma on K (`ms_ktest.log`, `ms_singlist.log`)

**(ML), statement.** For every h >= 0, n >= 1 (2h - 2 + n > 0) and `X_1, .., X_n in K`:
`<X_1 .. X_n>_h(v) - F_h^(n)(v) prod PiC(X_i)` is regular at v = 0 (`F_0 = -(5/2)v^2 log v`, n >= 3; `F_1 = -(1/12) log v`).
Equivalently, on the homogeneous basis of K: (ML1) `<H^n>_h = F_h^(n)(v)`; (ML2) a correlator of K-generators with at least one
`kappa(X, k)` has v-exponent >= 0 or vanishes. Equivalently: for `q in -z + K`, `Fbar_h(q) = F_h(Pi_0(q)) + Reg_h(q)` with all Taylor
coefficients of `Reg_h` at -z regular (`Pi_0(-z) = v`). In descendant language: on the slice `t^H_k = 0` (k >= 1) the singular part of
the genus-h descendant potential is `F_h(t^H_0)`, the c = 1 / Gaussian free energy of the single LINEAR coordinate t^H_0.
(ML) contains the primary-sector lemma (s.9 (c)) and the one-point law on K (s.3).

**Evidence (VERIFIED, exact over Q).** Every admissible correlator of K-generators (tau_0(H) and kappa(X, k)) for (h, n) = (0, 3..6),
(1, 1..5), (2, 1..3), (3, 1..3), (4, 1), with Sum k <= 3h - 3 + n, plus every SINGULAR one with larger k (kappa(X, k) has components of
lower psi-degree; none turned up): 3987 correlators with at least one kappa, 1019 of them with negative v-exponent, ALL 1019 equal to 0;
(ML1) in all 16 cases. Per type (h, n): (0,4) 6, (0,5) 37, (0,6) 152, (1,2) 3, (1,3) 19, (1,4) 83, (1,5) 313, (2,1) 2, (2,2) 17,
(2,3) 90, (3,1) 4, (3,2) 41, (3,3) 246, (4,1) 6 singular instances, all zero. The non-zero ones have v-exponents >= 0 (histograms in
the log; at genus 1 the exponent-0 values include the dilaton anomaly). Independently, on the basis insertions `tau_k(H^a)`:
of the 888 singular correlators computed (genus 0-4, n <= 3), the 361 whose insertions all lie in K have remainder exactly 0, while 489
of the 527 with an insertion off K have a non-zero remainder (`ms_singlist_full.log`, every singular correlator with its remainder).
**Prediction tested (genus 2, four insertions; `ms_ktest_24.log`).** (ML) was formulated from the data above; the table `amp(2,4)` was
computed afterwards: 1302 correlators of K-generators with a kappa, 376 of them singular, ALL 376 zero; `<H^4>_2 = 1/10 = F_2''''`;
on basis insertions 674 singular correlators, the 232 inside K with remainder 0, 398 of the 442 with D-content non-zero. Totals with
(2,4): 5289 / 1395 singular / 0 violations; 1562 singular basis correlators, 593 inside K, all with remainder 0.
**Status.** CONJECTURAL in general; PROVED on the `H^2`-tower (s.8); the primary part at genus >= 2 is proved by the TR worker
(allgenus/tr/TR_MODEL.md item 3, modulo arXiv:2502.18115: F_g depends only on the residue on the whole Hurwitz space), read for
comparison only.
