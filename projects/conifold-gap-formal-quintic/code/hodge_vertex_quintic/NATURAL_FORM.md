# Formal quintic at the conifold point: the natural factorisation (the analogue of the Hodge trick)

Date 2026-11-10/11 (s.5, genus 4: 2026-11-11; s.6, genus 5: point frame 2026-11-11, regularity collected from Stanage 2026-11-13). Exact finite computations modulo 2^89-1 and 2^61-1 (a referee also used 2^521-1), plus arguments written and audited
by AI agents only; no human has checked any of it. Conventions: OU_EXTENSION_genus3.md s.1 (hatted theory, `u = 1 - 5^5 q`, `s = u^(1/5)`,
`zhat = z/s^4`, `c = (0, -1, 3, 2, 1)`, Teleman form `R~(zhat; u)` of the hatted R-matrix, `M_0 = R~(zhat; 0)` the conifold symbol).
Detailed notes: `natural_form/lg/LG_MODEL.md` (the model), `natural_form/theory/UNIFORM_REDUCTION.md` (why), `natural_form/referee/REFEREE_NAT.md`
(independent check), `natural_form/SYNTHESIS.md` (critic: what stands, contradictions, overclaims in those three notes), verifiers' files in
`natural_form/verify_*/`. Where this note and those notes differ in the strength of a claim, SYNTHESIS.md s.1 and s.4 are the reference.
Status labels: VERIFIED (exact computation, orders and primes stated), ARGUED (a proof is written in the notes and was audited by a second
AI agent; it is not a refereed proof, and its hypotheses are named), CONJECTURAL.

## 0. Answer in one paragraph

There is a natural, parameter-free factorisation with the same shape as the Hodge one for the resolved conifold:

        M(u, zhat) = M_0( tau(u) zhat ),     1/tau(u) = 2F1(1, 4/5; 9/5; u) = Sum_n 4 u^n/(5n + 4),

i.e. the conifold symbol evaluated at `z` divided by the canonical coordinate measured from the conifold point,
`tau(u) zhat = -(5/4) zeta^i z/(u_i - u_i^c)`. `M_0` is not "tautological data of the quintic": it is the R-matrix of the Landau-Ginzburg model
`lambda = x^4 + v/x` with primitive form `dx` (five critical points coalescing onto a simple pole; the Hurwitz-space Frobenius manifold
`M_{0;3,0}`, `d = 1/2`, on its Z/5-symmetric line), so `M` is that model's R-matrix at the point `v(u) = u tau^(-5/4)` whose canonical coordinates
are those of the quintic. With this `M` and nothing else, every stable graph with edges is regular at the conifold point at genus 2
(6 graphs), genus 3 (41 graphs), genus 4 (378 graphs, of which the 149 with at most five edges carry the content, s.5) and genus 5
(4554 graphs, of which the 1663 with at most seven edges carry the content and 627 are regular only by cancellation, s.6), and the vertex
term carries the whole gap `(1/1200) t_c^-2`, `(1/25200) t_c^-4`, `(1/180000) t_c^-6`, `(1/660000) t_c^-8` (VERIFIED; two conifold frames,
both primes at genus 4; genus 3 and genus 5: second frame at one prime). The other factor `R' = Psi(u) G(u, zhat)` is the hatted form of a
gauge that is regular at the conifold point in the flat frame and is computable from the periods alone.

## 1. Statements and their status

1. **Factorisation (VERIFIED; all orders ARGUED modulo the purity hypothesis (H) below).**
   `R~(zhat; u) = G(u, zhat) M_0(tau(u) zhat)`, `G = 1 + Sum_{n >= 1} u^n G_n(zhat)`, `G_n` a polynomial in `zhat` of degree at most
   `floor((5n + 2)/4)` (attained in frame (0, 0): 1, 3, 4, 5, 6, 8, 9, ...; in special frames the degree drops). Checked: project code,
   `zhat <= 6`, `u <= 4`, both primes, two conifold frames (genus3_natural.py: 11 vanishing matrix coefficients against 0 free parameters);
   a referee's from-scratch rows, `zhat <= 16`, `u <= 11` at three primes and `zhat <= 24`, `u <= 17` at 2^89-1 and 2^61-1 (215 vanishing matrix
   coefficients); the project rows rebuilt at `zhat <= 16`, `u <= 11` (2^89-1). `M_0` must stand on the right: the same statement with `M_0` on
   the left is false.
2. **The rescaling is forced.** (a) ARGUED (modulo the hatted R-equation to all orders): at most one `tau` gives a polynomial `G`, and it is the
   canonical coordinate; neighbouring choices (`tau = 1`, `(1-u)^(4/9)`, the flat coordinate `(u/t_c)^(4/5)`) fail at the first order where
   they differ (VERIFIED, three primes). (b) VERIFIED on the graph sums (both primes): in the one-parameter family `M_0((1 + a u) zhat)` each of
   the four O(u) regularity equations (genus 2: graphs 2, 4; genus 3: graphs 2, 5) separately gives `a = -4/9 = tau_1`; the `[u^2]`, `[u^3]`
   conditions do not vanish on the rescaling direction (except `q_17`), so within rescalings `tau_2`, `tau_3` are fixed by regularity too.
3. **The model (VERIFIED through zhat^18 for all rows and zhat^40 for the identity row, both primes; all orders ARGUED modulo two standard
   inputs listed in SYNTHESIS.md item 4).** `R^_a^(i)(zhat) = N_a zeta^(i(a+1)) J_{c_a}(-4 zeta^-i zhat)`, `N = (1, -1, -1, 1, 1)`, where `J_k`
   is the normalised formal stationary-phase series of `Int exp(lambda/z) x^k dx`, `lambda = x^4 + u_LG/x`. The flat classes are
   `1, H, H^2, H^3, H^4 <-> dx, -dx/x, -x^3 dx, x^2 dx, x dx` (exponents `c_a`: the hatted weights are the Euler grading of the model; `H` is the
   derivative in the residue, the direction of the conifold logarithm). Gauss-Manin system: `z d/dv Y_b = beta_b Y_{b+1}` (`b` in Z/5),
   `beta = (1, -1/v, 1, -1, -1)`; scalar equation `(v (z d/dv)^5 + 3 z (z d/dv)^4 + 1) Y_0 = 0`, local exponents 0, 1, 1, 2, 3 (the conifold
   exponents of the quintic plus one). Closed recursion for the rows `g^a = J_{c_a}(-4 zhat)`:
   `g^(a+1)_n - g^a_n = ((4n - 3 - c_a)/5) g^a_(n-1)`, `Sum_a (4n + 1 - c_a) g^a_n = 0`: `M_0` to any order without the quintic (stored to
   `zhat^40` over Q). The primitive form `dx` is singled out by the data (`dx/x`, the mirror of P(1,4), has the same `M_0` but another unit).
   Not of Hodge/Gamma type; Mellin-Barnes kernel `Gamma(s) Gamma((s+1)/4)` in place of `Gamma(s)`.
4. **Regular gauge (ARGUED: Theorem U of UNIFORM_REDUCTION.md s.3, audited line by line by a second agent; VERIFIED).** Two cyclic systems
   `z d/du Y_a = b_a Y_(a+1)` with unit coefficients except for one simple pole and with equal products `prod b_a` are related by a formal gauge
   `P(u, z)` whose coefficients are regular at `u = 0`. Equality of the products is equality of the canonical one-forms, which forces
   `v = u tau^(-5/4)`. `G` is the hatted form of `P`. The gauge computed from the periods alone (`varpi_0^cf`, `t_c`, `tau`) reproduces the
   actual `G`: all 397 Fourier components for `zhat <= 14`, `u <= 11` as exact rationals, with no fitted constant once the normalisation
   `P_11(0, z) = 1` is imposed (that normalisation is itself only verified, at `z^5, z^10, z^15, z^20`). Support rule (VERIFIED as a bound,
   `zhat <= 24`, `u <= 17`, seven frames at lower order): `G - 1` is regular at `u = 0` in the frame `K = (1, uH, H^2, H^3, H^4)`. Consequences:
   the degree bound; `G_1 = zhat X_1`, `X_1 = -(8/225 - w_1/5)[1 (x) H^2 + H^2 (x) 1] - (1/75)[H^3 (x) H^4 + H^4 (x) H^3]` (Fourier normalisation of
   genus3_natural.py; `w_1 = [u^1] varpi_0^cf`; the second coefficient is frame independent), with no `H (x) H`: a point of the kinematic plane
   Pi of OU_EXTENSION_genus3.md s.6d, moved by the choice of conifold frame along the string direction only; `H (x) H` enters from `u^2` on and
   only inside the rank-one block `alpha zhat (H - zhat)(x)(H + zhat)` (two Euler-dilaton vectors `uH -+ z`), `alpha = -127/15750 - (2/5) beta` in
   the frame `(alpha, beta)`.
5. **Graph sum (VERIFIED).** `M = M_0(tau zhat)` mod `(zhat^7, u^4)`, right placement, `R' = Psi(u) G`:

   | run | genus-3 edge graphs | vertex term | genus 2 |
   |---|---|---|---|
   | frame (0, 0), 2^89-1 | all 41 are O(u^4) | `u^-4 (1/25200 - u/9000 + 407 u^2/3780000 - 151 u^3/3780000)` = `(1/25200) t_c^-4 + 0 + 0 + 0` | six edge graphs O(u^2); vertex `[u^0] = 1/1200` |
   | frame (0, 0), 2^61-1 | all 41 are O(u^4) | same | same |
   | frame (1/3, -2/7), 2^89-1 | all 41 are O(u^4) | `u^-4 (1/25200 - 11 u/189000 - ...)` = `(1/25200) t_c^-4 + 0 + 0 + 0` | same |

   No parameter was adjusted. In the solution set of s.6d this is one point: first order `N_1 = N_1^taut - M_0^-1 (zhat X_1) M_0` (on Pi; not the
   tautological point, not pure-high, not one of the 17 rays); the content beyond the kinematics of Pi is that the nine numbers which s.6d
   solved for with free `N_2`, `N_3` (`[u^2]`, `[u^3]` of the two one-edge graphs, `[u^3]` of the five two-edge graphs) vanish for the canonical
   `N_k = [u^k] log(M_0^-1 M_0(tau zhat))`, in each frame. (The frame-(0, 0) logs were produced before the `--ab` option was added to the
   script; the job specifications are identical. The second frame was not run at the second prime.)
6. **Regularity is not a property of the support rule alone (VERIFIED, both primes; scratch computation, linear responses only).** A left
   dressing of `M` by a single K-regular term is in general visible: at `u^2` the terms `zhat H (x) H`, `zhat^2 (H (x) 1 - 1 (x) H)`, `zhat^3 1 (x) 1`
   have non-zero responses on the one-edge graphs (`(-1/400, -1/7200)`, `(1/300, 7/14400)`, `(1/1200, 1/2880)` on graphs 2, 5) which cancel
   exactly in the Euler-dilaton combination that occurs in `G_2`; the Pi-type terms and `zhat^2 (H^4 (x) H^2 - H^2 (x) H^4)` are invisible. So the
   cancellation uses the residue relations of the gauge (UNIFORM_REDUCTION.md s.5), not just regularity in the frame `K`.
7. **Conjecture N (VERIFIED for g = 2, 3, 4, 5 - genus 5 in s.6, with the scope and the caveats stated there; CONJECTURAL for
   `g >= 6`).** For every `g >= 2`, with `M = M_0(tau(u) zhat)` and `R' = Psi(u) G`, every stable graph with at
   least one edge is regular at `u = 0`; equivalently `Omega^quintic = P . (Omega^LG at v(u))` with `P` regular and the vertex term equal to
   `c_g t_c^(2-2g)` + regular, `c_g = (-1)^(g-1) B_2g/(2g(2g-2)) 5^(1-g)`.

## 2. Conclusions from the genus-3 runs (s.6d of OU_EXTENSION_genus3.md re-read)

* The classification of s.6d said that genus 3 does not single out the tautological factorisation: besides it there are the kinematic plane
  Pi, the invisible orders 4-6 and 17 isolated rays. That remains true: through genus 3 the graph sums alone do not select a factorisation.
  What selects one is the reduction principle: the first u-correction of the R-matrix is the rescaling generator `tau_1 M_0^-1 zhat dM_0/dzhat`
  plus a constant `zhat`-linear dressing that lies exactly in Pi. So Pi is not an accident of genus 3: at order `u` the K-regular terms are
  exactly `zhat` times the two propagators spanning Pi, i.e. Pi is where the first-order regular gauge lives.
* The scalar (Hodge-type) corrections of the genus-2 note, `(b_1, b_3) = (113/75, -21/125)`, fitted two equations with two numbers; the
  object behind the O(u) equations is the rescaling, which is colour-mixing and needs one number, predicted in advance.
* The trick has the same shape as the Hodge one: *R = (gauge regular at the singular point) x (R-matrix of the universal local model in its
  natural variable)*. Resolved conifold: model `Int exp(-(x - w log x)/z) dx/x` (Gamma function, Stirling series, natural variable `z/w`,
  constant in `q`). Formal quintic: model `Int exp(-(x^4/4 + v/x)/z) x^m dx`, natural variable `z` over the canonical coordinate, which
  degenerates at the conifold point and therefore must depend on `u`. Side-by-side table: UNIFORM_REDUCTION.md s.7. The counterpart of the
  Hodge factor is `M_0(tau zhat)` together with the regular gauge, not the symbol `M_0` alone (which leaves simple poles in the one-edge
  graphs).
* The vertex theory is explicit: the Givental CohFT of the Hurwitz Frobenius manifold of `p^4 + a_2 p^2 + a_3 p + a_4 + a_5/(p - a_1)` (form
  `dp`) on its residue line. Its genus `>= 1` primary correlators, as far as they were computed, are those of the Gaussian/c = 1 free energy
  in the residue (`-(1/12) log v`, `B_2g/(2g(2g-2)) v^(2-2g)` up to the normalisation `5^(1-g)`, `g = 2, 3`). By a theorem quoted from the
  literature (arXiv:2211.12259, read only in part, not used) its total descendant potential is a tau function of the extended rationally
  4-constrained KP hierarchy.
* The vertex term is not the bare model potential: the leafless vertex is `c_g (u tau^(-3/2)/xhat)^(2-2g)`, and the dilaton leaves of `R'`
  convert the model coordinate into `t_c` (observed bookkeeping: `-(g-1) c_g/15` at order `u`; not derived).

## 3. What is not established

* Everything on the graph sum is genus 2 to 5, modulo two primes. At genus 3 the evidence beyond the kinematic first order was nine
  rational numbers per frame; genus 4 (s.5) adds the `[u^0..u^5]` coefficients of 149 outer graphs with at most five edges per frame, the
  first place where the diagonal `u^4 zhat^5` terms of the gauge and the genus-1 anomaly of the Euler-dilaton vector in multi-edge graphs
  enter; genus 5 (s.6) adds 3142 coefficients on 1663 outer graphs, of which 1655 vanish by cancellation. Genus 6 is out of reach with this
  code (the graph count and `K` both grow). At genus 5 the arithmetic evidence is probabilistic rather than deductive (s.6), the second
  frame's regularity was run at one prime (its total is two-prime, s.6), and `tau_6` and the orders `zhat^10..zhat^12` of the symbol are
  untested by any control.
* No mechanism is known for the cancellations beyond first order. The Euler-dilaton combination explains why the `H (x) H` block of `G_2` drops
  out of the one-edge graphs at genus `>= 3` (first order; checked at both primes at genus 3), but it is not an exact symmetry of the model
  theory (genus-1 anomaly `-1/8`), and the other `[u^2]`, `[u^3]` cancellations are observed, not derived.
* (H): the hatted rows are power series in `u` (purity). Verified to `zhat^16`, `u^11` on the project pipeline and to `zhat^24`, `u^17` in the
  referee's rows; not proved. The all-order factorisation, the uniqueness of `tau` and the normalisation of the gauge depend on it.
* The all-genus theory of the model: the Bernoulli gap `F_g = B_2g/(2g(2g-2)) a_5^(2-2g)` and the vanishing of the genus `>= 1` primary
  correlators away from the residue direction are known only for the `(g, n)` needed at genus 3. This is the missing analogue of the closed
  Hodge-integral evaluation of the vertices.
* No human has refereed Theorem U, the all-order identification of the symbol, or the uniqueness lemma.
* Left placement, the second frame at the second prime, and the non-formal quintic were not touched.

## 5. Genus 4 (2026-11-11)

**Method (new code `fq4_split.py`, driver `genus4_split.py`).** No vertex tensors: the edge bivector of `R = R' M` splits as
`E_R = E_M + E_out`, `E_out = M^-1 E_R' M^-T = E_R - E_M`, with `E_M(a, b) = tau^(a+b+1) Eh_0(a, b)` for `M = M_0(tau zhat)` and the leaves of
`R` unchanged, so the point-frame sum over the leg-free stable graphs with every edge marked inner or outer, regrouped by the contraction
of the inner edges, is the outer graph sum `R'.(M.omega^)`; all vertices are Witten-Kontsevich numbers with one colour, and each term is a
colour sum done by vertex elimination. `M` enters exactly in `u` (no truncation of `tau`), mod `zhat^(3g-2)`; series mod `u^(2g-2)`.
**Validation:** `py -I genus4_split.py --regress` reproduces, at both primes, the stored per-graph tables of the old vertex-tensor code at
genus 3 for `M = 1`, for `M = M_0` (all 42 graphs through `u^5`, non-zero edge graphs) and for the natural `M` (through `u^3`), and the
genus-2 values of the notes. For every inner graph the split pieces add up to the point-frame term. Cost: genus 3 in 7 s; genus 4 in
about 4.5 minutes of wall time with 13 processes (54 CPU minutes), so no cluster run was needed.

**Results (379 genus-4 graphs: 1, 3, 7, 21, 43, 75, 89, 81, 42, 17 with 0..9 edges; a term is `u^-6 G(u)`).**

| run | total | natural `M`: outer graphs with edges | vertex graph |
|---|---|---|---|
| frame (0, 0), 2^89-1 and 2^61-1 | `u^-6 (1/180000 - 7u/300000 + 701u^2/18000000 - 1459u^3/45000000 + 923639u^4/67500000000 - 171583u^5/67500000000)` = `(1/180000) t_c^-6 + 0 + 0 + 0 + 0 + 0` | all 378 are `O(u^6)`: regular | equals the total through `u^5`: `(1/180000) t_c^-6 + 0 + 0 + 0 + 0 + 0` |
| frame (1/3, -2/7), 2^89-1 and 2^61-1 | `u^-6 (1/180000 - 11u/900000 - 79u^2/378000000 + ...)` = `(1/180000) t_c^-6 + 0 + 0 + 0 + 0 + 0` | all 378 are `O(u^6)` | carries the whole polar part |

So the gap holds for the formal quintic at genus 4 (`1/180000 = (-B_8/48)/125`), and Conjecture N holds at genus 4 in both frames at
both primes. How much of "378 regular graphs" is content: the 229 outer graphs with at least six edges are `O(u^6)` for free (every
outer edge is `O(u)`); of the 149 with `E <= 5`, 81 (3, 19, 59 with 3, 4, 5 edges) receive only pieces that vanish identically (a
selection rule of the true `E_out`: a generic `E_out` of the same valuation gives non-zero pieces), and 68 (3, 7, 18, 24, 16 with 1..5
edges) vanish by cancellation among 6767 non-zero (inner graph, outer type) pieces.
**Controls (2^89-1, frame (0, 0)):** with the bare symbol `M = M_0` (`tau = 1`) all 149 outer graphs with `E <= 5` are singular; with the natural
`tau` and only its `u^5` coefficient shifted by 1, exactly the three one-edge graphs fail, at `[u^5]`, and the vertex graph acquires a
`t_c^-1` term: the test resolves the last coefficient of `tau` that genus 4 can see.
Files: `g4_split_g4_<p>_<variant>[_ab..].json` (per-outer-graph tables and the list of outer graphs), `parts/` (per-inner-graph pieces and
timings), `logs/genus4_split_*.log`.

**Refereeing (AI agents; `natural_form/genus4_referee/`; summary and wording advice in `SYNTHESIS_G4.md`).**
* `indep/INDEP.md`: a second code sharing no module with `fq4_split.py` (rows of the earlier from-scratch referee, own edge tables, own
  Witten-Kontsevich numbers, the Mathematica graph list with own automorphism counts, a different summation algorithm in Fourier space)
  reproduces, at both primes in frame (0, 0), all 379 outer-graph values of the natural run, all 379 point-frame terms, every one of the
  non-zero pieces of the part files, and both controls. It did not redo the second frame.
* `review/REVIEW.md`: a written argument that the regrouped split sum equals the outer graph sum graph by graph (automorphism factors by
  orbit-stabiliser and an exact sequence for nested graphs; leaves; `xhat` powers), with its combinatorial core checked for all 379
  genus-4 graphs; 26 code mutants, 23 caught by the genus-2/3 regression and the other three (`tau_4`, `tau_5`, powers of `tau` at
  `zhat >= 7`, invisible at genus 3 by construction) caught by the genus-4 checks. It found two latent defects in the driver (a
  prefix collision of part-file names in `--collect`, no disjointness check of part files), which did not affect the stored results and
  are fixed (the regression and all six genus-4 collections were repeated after the fix), and it points out that the GAP check and the
  per-inner-graph additivity do not test `E_M`; only the regularity check and the genus-3 regression do, and the regression cannot see
  `tau_4`, `tau_5` or the orders `zhat^7..zhat^9`: for those the evidence is the genus-4 regularity itself plus a direct check of the
  tables (`M` symplectic, `E_out = M^-1 E_R' M^-T`, the support rule of the gauge through `u^5`, `zhat^9`; `review/rv_compose.py`).
  "Vertex graph carries the polar part" and "all edge graphs regular" are one check, not two.
* `genus5/GENUS5_ESTIMATE.md`: cost of genus 5 (4555 graphs, `K = 8`, `zhat^12`; nothing at genus 5 was computed). Point-frame total
  alone: about 1.3 CPU hours per prime (local). Per-outer-graph regularity with the natural `M`: about 277 CPU hours with the present
  code, about 46 CPU hours with a prototype that does one elimination for all outer-edge subsets (validated on all 379 genus-4 graphs);
  slowest single graph 4.1 h, resp. 0.7 h. That is a cluster job (SLURM array, one core per task; Python >= 3.8); it also needs a faster
  graph generator (the present canonical form needs about 30 minutes for the genus-5 list; a refinement-based one takes 10 s) and
  cost-balanced tasks. (All of this was then built and run: s.6, s.6a. The estimate of 46 CPU hours and of the slowest graph held.)
* What the agreement cannot show: that the regrouped sum is the outer graph sum at genus 4 is an argument (genus independent), checked
  numerically against genuine vertex tensors only at genus 2 and 3; the rows come from two derivations, neither refereed by a human.

## 6. Genus 5 (point frame 2026-11-11; regularity on Stanage, collected and audited 2026-11-13)

**Conjecture N holds at genus 5** (VERIFIED in the sense of this note: exact finite computation, primes and scope as stated below;
evidence, not proof). With `M = M_0(tau(u) zhat)` and no parameter adjusted, every one of the 4554 leg-free genus-5 stable graphs with at
least one edge is regular at `u = 0`, and the vertex graph carries the whole polar part `(1/660000) t_c^-8`.

**The runs.** Package `../cluster_genus5_quintic/` (s.6a), run by the user on Stanage (SLURM, 100 one-core array tasks per prime). The
cluster's logs are now in the package (`../cluster_genus5_quintic/logs/`, 3.5 MB) and its per-graph pieces were copied back as
`../g5_evidence.tar.gz` (46 MB compressed, 113 MB of `parts/`, 1021095 stored pieces), so every number below has been re-derived from the
stored pieces rather than read off the collected tables. A genus-5 term is `u^-8 G(u)` and each table row is the eight polar coefficients
`[u^0..u^7]` of one outer graph.

| run | tag | outer graphs with edges | vertex graph | `failed` |
|---|---|---|---|---|
| natural `M`, frame (0, 0), 2^89-1 | `g5_89_nat` | all 4554 regular (628 rows stored, 627 of them zero) | `(1/660000) t_c^-8 + 0 + ... + 0` | `[]` |
| natural `M`, frame (0, 0), 2^61-1 | `g5_61_nat` | same, same 628 indices | same | `[]` |
| natural `M`, frame (1/3, -2/7), 2^89-1 | `g5_89_nat_ab1o3_m2o7` | same, same 628 indices | same | `[]` |
| control: bare symbol `M_0`, `--maxout 2` | `g5_89_M0_out2` | **14 singular** = every graph it computes | spoilt at `t_c^-7..t_c^-1` | 2 failures |
| control: `tau_7` shifted by 1, `--maxout 2` | `g5_89_nat_t7_out2` | **3 singular**, the one-edge graphs, at `[u^7]` only | spoilt at `t_c^-1` only | 2 failures |

**The total, and the gap.** Reconstructed from the two primes by CRT (the `u^6`, `u^7` coefficients were beyond one prime: the single-prime
values printed with a leading `~` in the November 11 point-frame logs are wrong, and the lattice certificate does not detect that):

        total = u^-8 (1/660000 - 7u/825000 + 199u^2/9900000 - 717u^3/27500000 + 449357u^4/22500000000
                      - 5593217u^5/618750000000 + 14686042363u^6/6496875000000000 - 8370858511u^7/32484375000000000)
              = (1/660000) t_c^-8 + 0 t_c^-7 + ... + 0 t_c^-1,

with `1/660000 = |B_10|/(10 . 8 . 5^4)`, the pattern `c_g = |B_2g|/(2g(2g-2) 5^(g-1))` of genus 2-4. Equivalently
`G = c_5 (u/t_c)^8 mod u^8` (the same statement packaged: `(u/t_c)^8 = 1 - 28u/5 + 199u^2/15 - 2151u^3/125 + 4942927u^4/375000 -
5593217u^5/937500 + 14686042363u^6/9843750000 - 8370858511u^7/49218750000`). The seven vanishing subleading coefficients are not an
artefact of the triangular change of variable: the Jacobian `dT_k/da_j = [t^k](A^-8 u(t)^j)` is lower triangular with 1 on the diagonal, so
a perturbation of the `u`-series at order `j` moves exactly `T_k` for `k >= j`; `T_1 = ... = T_7 = 0` are seven independent affine conditions,
and given `T_0 = c_5` there is exactly one `u`-series mod `u^8` satisfying them, namely the one computed.

**How much of "4554 regular graphs" is content** (the genus-4 accounting of s.5, continued; graphs by edge number
1, 3, 11, 34, 100, 239, 492, 784, 1002, 926, 632, 260, 71 for `E = 0..12`):
* **2891 graphs with `E >= 8`: regular a priori, not computed.** Every entry of `E_out = E_R - E_M` has `u`-valuation exactly 1 (all 1950
  entries, 78 edge keys; `[u^0] E_R = Eh_0` and `tau(0) = 1` force the vanishing, and non-vanishing at `u^1` was measured), so a row with `E`
  outer edges has valuation at least `E`. The engine therefore caps `|S|` at `K - 1 = 7` and stores no row for them. This is a one-line
  valuation argument, not evidence.
* **1036 of the remaining 1663 (7, 71, 288, 670 with 4, 5, 6, 7 edges) receive no non-zero piece at all** - the selection rule of the true
  `E_out` already seen at genus 4.
* **627 vanish by cancellation** (3, 11, 34, 93, 168, 204, 114 with 1..7 edges: every graph with `E <= 3` and a falling fraction above).
  1655 of the 3142 polar coefficients at stake are cancellations. The cancelling sums are large: the three one-edge rows are sums of 3479,
  2296 and 4525 non-zero (inner graph, outer type) pieces, and removing a single piece leaves a non-zero sum; the vertex graph is a
  4555-term sum, every inner graph contributing. (Genus 4: 229 free / 81 selection / 68 cancellation, 6767 pieces.)

**What the run does *not* show, although the checks pass.** Two of the five collect checks are independent of `M`:
* `c_5 = 1/660000` is **not** a confirmation of the factorisation. Since `E_out = O(u)`, `[u^0]` of the vertex graph equals `[u^0]` of the
  whole point-frame total for any `M` with the frozen symbol at `u^0`; both deliberately broken controls report
  `c_vertex_equals_pattern = True`, and the value was already established by the point-frame runs of November 11.
* the GAP check on the total is likewise `M`-independent (the `id` runs pass it), and "the vertex graph carries the polar part" and "all
  edge graphs are regular" are one check, not two (as `genus4_referee/review` already noted at genus 4).
So exactly one of the reported checks - the regularity of the edge graphs - is about the factorisation.

**Controls: what they establish.** A shift of `tau` at order `k` moves `E_M` and `E_out` at valuation exactly `k`
(`(tau + c u^k)^m - tau^m = m tau^(m-1) c u^k + O(u^(2k))`, verified on all 1950 entries for `k = 1..7`), so by multilinearity in the edge
factors a row with `E` edges moves at valuation at least `k + E - 1`, visible only if `k + E - 1 <= 2g - 3 = 7`. Hence `k = 7` can be seen on
the one-edge graphs only, and only at `[u^7]` - exactly the observed pattern - while the bare symbol is `k = 1` and spoils everything with
`E <= 7`, from `[u^E]`. The same prediction reproduces the genus-2, genus-3 and genus-4 control patterns (s.5: `tau_5`, 149 graphs).
Moreover the floor-order response of a one-edge row is `c . A_Gamma` with `A_Gamma` *independent of `k`*, which the two controls measure twice
over: `[u^1](M_0 - nat) = (4/9) [u^7](nat_t7 - nat)` on the vertex and the three one-edge rows (all four identities exact mod 2^89-1, all four
`A_Gamma` non-zero; the same identity holds on the genus-4 tables, with `tau_1 = -4/9`). Since `A_Gamma != 0` and `k` drops out, *any* rescaling
deviating from `tau` at any order `k <= 7` makes the one-edge graphs singular: within the ansatz `M = M_0(tau zhat)`, `tau` is the **unique**
rescaling mod `u^8` for which the one-edge graphs are regular. Also `Sum_Gamma (delta row) = 0 mod u^8` (the total is `M`-independent) holds in
all three instances where the controls computed every row that contributes.

**Independent verification of the collected result** (`natural_form/genus5_cluster/`; four AI agents, 2026-11-13).
* `recheck/RECHECK.md`: own CRT and lattice reconstruction (margins 10^18 to 10^39 against the modulus), own Frobenius solution of the
  Picard-Fuchs equation and own `t_c` (`t_c = u + (7/10)u^2 + (41/75)u^3 + ...`; second frame `u + (11/30)u^2 + ...`), which reproduces the
  stored `t_c` expansions residue for residue and so independently confirms `fq_hat.tc_series`/`fq3_msum.tc_expansion3` at the order used;
  own `B_10 = 5/66`; the census above; the seven JSONs mutually coherent (identical graph list, row sums equal to the stored totals,
  `nat` total = `id` total at both primes, the two frames sharing the `t_c` polar part).
* `local/LOCAL.md`: the package engine re-run **here** over all 4555 inner graphs at `--maxout 1` (55 minutes of wall time, 13 workers)
  reproduces the cluster's vertex row and its three one-edge rows exactly (32/32 residues), and the `--maxout 1` pieces are exactly the
  restriction of the full pieces where both runs overlap. Per-inner-graph additivity (the pieces of an inner graph add up to its
  point-frame term) was redone locally with complete piece families for **2586 of the 4555** inner graphs, 222784 non-zero pieces,
  0 failures. (Whether the cluster itself ran that check was open at the time; it is settled below.)
* `critic/CRITIC.md`: the graph list `data/graphs_5_0_12.json` proved **complete and non-redundant** by an independent induction on the edge
  count (one-edge expansions of each level, own canonical form; also at genus 2, 3, 4), `|Aut|` recomputed for all 4555 graphs by
  backtracking, and the piece-to-outer-graph assignment - which no genus-5 check constrains, since both the total and the additivity sum
  over `S` - recomputed for 168475 `(inner graph, S)` pairs with an own contraction and own canonical form, 0 disagreements; `graph_dp`
  against `Split.graph` on 80 genus-5 graphs; the Witten-Kontsevich inputs checked against closed formulas and the dilaton equation; a
  third prime 2^521-1 on 31 inner graphs consistent with both stored residues (248 coefficients). It also establishes that a short run
  could not have produced these files: all 4555 point-frame terms are non-zero at both primes and the natural total equals their sum, so a
  missing inner graph would move the total.
* `controls/CONTROLS.md`: the control analysis quoted above.

**The copied-back evidence** (the cluster logs and its 1021095 per-graph pieces; second wave, three AI agents, 2026-11-13). Everything in
the six paragraphs that follow was re-derived from the cluster'''s own part files rather than read off `results/*.json`.

**Provenance, recovered from the copied-back logs.** SLURM arrays 12149822 (`g5_89_nat`), 12149824 (`g5_61_nat`), 12156906 (second
frame), 12156902 and 12156904 (the two `--maxout 2` controls), with collect jobs 12149823/12149825/12156907/12156903/12156905 and validate
job 12149764; Python 3.9.25 on `node003-013.stanage.hpc.shef.ac.uk`; per-task wall time 1716-3570 s; CPU 51.2 / 53.3 / 52.6 h for the three
natural runs and 14.7 / 18.6 h for the controls, about 190 CPU hours in all, and 1.8 h for each of the two point-frame runs on the
workstation (Python 3.14.6, 13 processes). `V_validate.sh` ran on the cluster and passed 31 checks with 0 failures, so the genus 2-4 tables
of record were reproduced on a second machine with a second interpreter; it measured the machine factor `f = 0.94`, which makes the
predicted 0.66 h for the heaviest task right. No task was resumed or restarted (one header line and exactly 4555 `inner graph` lines per
run), and the task decomposition is the documented one: the longest-processing-time slices recomputed from `data/cost_table_g5.json` agree
task by task in all five cluster runs (100/100 and 20/20), every inner graph lies in exactly one task file, and the union is `0..4554`.

**The audit of the stored pieces** (`natural_form/genus5_cluster/census/CENSUS.md`). All seven result tables were rebuilt row by row and
residue by residue from the part files by code sharing no module with the package: 0 rows on one side only, 0 differing residues, row sums
equal to the stored totals, and the outer-graph list identical index by index in all seven files. There are 340363 non-zero
`(inner graph, outer graph)` pieces at `p89`, the same 340363 at `p61` and 340369 in the second frame; by `|S|` they are 4555, 10300, 24815,
49315, 82393, 91858, 60835, 16292 - so the `|S| <= 2` block, 39670 pieces, is exactly what the two `--maxout 2` controls computed - and by
inner-graph edge number 1, 6, 39, 209, 1013, 3959, 12926, 32503, 63121, 86067, 81340, 44843, 14336. The content split above (2891 free /
1036 selection / 627 cancellation, 3142 coefficients at stake, 1655 of them in rows that receive pieces) is re-derived graph by graph from
the pieces and comes out identical at both primes and in both frames. Refined: of the 1655, **517 receive no non-zero contribution at that
order and 1138 are genuine cancellations**, 1075 of them sums of more than 100 non-zero residues, and not one is a sum of just one or two
pieces; pieces per cancelling row min 31, median 347, max 4525. The valuation law (a piece landing on an `|S|`-edge outer graph has its
first `|S|` coefficients zero) holds on all 1021095 stored pieces with 0 violations, and the floor is attained for 44, 38, 32, 29, 30, 41%
of the pieces with `|S| = 1..6` (always for `|S| = 0` and `7`) - the selection rule visible on single pieces. The two primes agree on the
vanishing pattern of every one of the 2722904 residue slots, so no coefficient of the evidence rests on one prime alone and there is not one
accidental mod-`p` zero in the computation.

**Additivity: settled.** The copied-back collect logs say in words that the package ran and passed the per-inner-graph check for both
`g5_89_nat` and `g5_61_nat` (`PASS ADDITIVITY: for every inner graph the pieces add up to its point-frame term`) and that it was skipped in
the second frame for want of an `id` run there. It was then redone from the stored pieces by two separate audit codes, neither of them
calling `g5_checks.checks`: **all 4555 inner graphs at both primes, 0 failures** (340363 pieces summed per prime); and, after making the
missing point-frame run in the second frame here (`g5_89_id_ab1o3_m2o7`, `g5_61_id_ab1o3_m2o7`; 4555 graphs each, 650 s and 1003 s of wall
time with 13 workers, both passing their own GAP check with `c_5 = 1/660000`), **all 4555 in the second frame as well**, 0 failures over its
340369 pieces. 13665 identities in all. No inner graph has a vanishing point-frame term, so none of them is trivially `0 = 0`.

**The second frame now has a two-prime, certified total** (`natural_form/genus5_cluster/frame2/FRAME2.md`):

        total = u^-8 (1/660000 - u/225000 + 323u^2/207900000 + 66841u^3/15592500000 - 32135797u^4/46777500000000
                      - 917102699u^5/350831250000000 - 3515716092179u^6/3683728125000000000
                      + 1892778904559u^7/5023265625000000000)

(largest height `2^62.1`, smallest certificate margin 45.6 bits, so unique rather than probabilistic). It differs from the frame-`(0, 0)`
total at **every** order `u^1..u^7`, and `t_c = u + (11/30)u^2 + (751/1575)u^3 + ...` differs too, yet the polar part computed over `Q` is
`(1/660000) t_c^-8 + 0 + ... + 0` in both frames: the seven vanishing conditions come out zero from different numbers in different
coordinates. (Of the `~`-marked single-prime values in that frame's collect log, `u^3` is confirmed and those at `u^5` and `u^6` are
spurious, as warned above.)

**New structural fact: the selection rule is frame independent on rows and frame dependent on pieces.** Comparing the two full natural runs
at `p89`: the same 628 outer graphs have a stored row, the same 627 of them cancel to zero with the same `u`-support, and the vertex graph is
the only non-zero row in both - so "1036 by selection (7, 71, 288, 670 at `E = 4..7`), 627 by cancellation" is a statement about the theory,
not about the frame `(0, 0)`. Of the 340429 pieces occurring in either frame, however, 66 occur in one frame only - all with `u`-support
exactly `{u^7}`, the top of the `K = 8` window - and 124 more differ in `u`-support; the lowest order at which any support differs is `u^3`
(nothing at `u^0`, `u^1`, `u^2`), and all 190 differences have `E(t) >= 7`, `E(q) >= 3`. The same comparison between the two primes in the
same frame finds 0 differences among 340363 pieces, so the frame dependence is genuine and not a mod-`p` accident.

**The independent re-implementation, run against the cluster's own pieces** (`natural_form/genus5_cluster/indep/INDEP_G5.md`). This is what
`parts/` was copied back for. The second implementation of the graph sum (`genus4_referee/indep/`: own `E_R`, leaf, symbol and vertex tables
on the referee's rows, own Witten-Kontsevich numbers, own automorphism count, contraction and canonical form, momentum-space summation; no
module shared with `g5_engine`/`fq4_split`) was run on **1474 of the 4555** inner graphs - every edge number `E = 0..12`, the complete strata
`E <= 5` and `E = 12`, the vertex graph, the package's six heaviest graphs and the loop-heavy 149, 245, 2668 that no independent code had
reached - and **93551 pieces = 748408 residues agree entry by entry at 2^89-1**, with 19591 pieces = 156728 residues on 460 of them at
2^61-1: 0 mismatches, 0 pieces produced where the cluster stored none, 0 stored pieces missing (19 CPU hours here, no cluster). The
independent point-frame terms equal the shipped `parts/g5_<p>_id/local.json` for 1011 (`p89`) and 460 (`p61`) inner graphs, and for every one
of them the independent natural pieces add up to the independent point-frame term, so the regrouping identity is now tested in code that
cannot share the package's conventions. `|Aut|` and the canonical forms were recomputed independently for all 4555 graphs (4555 pairwise
distinct forms, every `|Aut|` equal), and for the outer graphs `q = 839` (`E = 6`, 35 contributors) and `q = 1502` (`E = 7`, 31 contributors)
every inner graph the cluster records there was re-run independently, so two of the 627 cancellations are now reproduced end to end.

**What is not established at genus 5.**
* **The second frame's regularity is single-prime** (2^89-1). Its total is now two-prime and certified and its 4555 additivity identities
  hold, but its 628 rows rest on one prime; a second-prime natural run in that frame is another ~52 CPU hours, i.e. a cluster job.
* **The independent coverage is a sample.** 1474 of 4555 inner graphs at `p89` (32%) and 460 at `p61`, 93551 of the 168475 non-zero pieces
  (56%); the strata `E = 6..11` are sampled, not complete (283/492, 213/784, 185/1002, 146/926, 115/632, 73/260). No outer row is complete in
  the strict sense: for `q = 839` and `q = 1502` every inner graph the cluster records there was re-run, but that no *further* inner graph
  contributes to those rows is still the cluster's claim, and the vertex row and the three one-edge rows (2296-4525 contributors) are out of
  reach. The independent code's rows and `tau` are a second derivation of the same analytic input, not a third; the genus-5 graph list is
  shared (`data/graphs_5_0_12.json` - there is no independent Mathematica export above genus 4, so completeness of the list remains the
  critic's induction, while non-redundancy, `|Aut|`, canonical forms and the piece assignment are recomputed); and the second frame and the
  two controls have no independent piece at all.
* **`tau_6` is untested by any control, at any genus** (the genus-4 window stops at `u^5`; only `tau_7` was shifted at genus 5), and so are
  the orders `zhat^10..zhat^12` of `M_0`, which genus 5 is the first to touch. Both controls were run at `--maxout 2`, so they discriminate
  on 14 of the 1663 content graphs, while 613 of the 627 cancellations sit at `E >= 3`. A `nat_t6 --maxout 2` run (about 10-15 CPU hours)
  and one control whose perturbation is invisible below `E = 4` at full `--maxout` are the controls actually missing.
* **The arithmetic evidence is now probabilistic.** The genus-4 argument ("a rational of this height cannot vanish mod both primes") no
  longer applies: per-graph point-frame coefficients reach height 2^136 against `p89 . p61 = 2^150`, a margin of 14-50 bits, and 69% of them
  do not reconstruct from the two primes at all. The failure probability is about 2^-138 over the 3142 coefficients at stake - negligible,
  but a bound, not a deduction. The two primes are also not independent as a test of the code (same engine, same rational inputs reduced
  mod two primes), which is why the inputs were checked separately and a third prime was spot-checked.
* **Scope is inherited from genus 4, not from genus 2-3.** The run tests the regrouped split sum; its identification with the outer graph
  sum `R'.(M.omega^)` is an argued, AI-audited identity checked against genuine vertex tensors only at genus 2 and 3 (`data/ref_g4.json` is
  `fq4_split`'s own output, so the genus-4 regression is self-referential on this point). The independent re-implementation shares the
  split-sum method, so a structural error in that identification is common-mode even where the two codes agree piece by piece; and all
  seven runs of the package share `g5_engine`.
* **The result files still carry no provenance of their own** - no host, date, CPU time, part-file inventory or code digest; what is
  recorded above was reconstructed from the copied-back logs, which is why `logs/` and `parts/` had to be kept (`../g5_evidence.tar.gz`).
  Without them no single non-trivial outer graph could be re-checked by anybody, each being a sum over hundreds or thousands of inner
  graphs.
* For the next run, `g5_checks.checks` should record every check as pass/fail/**skipped** in `info`, assert `min(valO) == 1` together with a
  non-zero entry of `E_out` (as written, the check would also pass if `E_out` vanished identically, which is exactly the claimed output),
  and stamp `info` with host, interpreter, CPU seconds and a digest of the `.py` files.

### 6a. The cluster package

`../cluster_genus5_quintic/` (tarball `../cluster_genus5_quintic.tar.gz`), pure Python >= 3.8, no packages, relative paths only; on Stanage it
lives in `~/formalquintic`. `genus5_run.py --validate | --task I --ntasks N | --run | --collect`; `g5_engine.graph_dp` does one vertex
elimination for all outer-edge subsets (about 5x faster than `fq4_split.Split.graph`, which it reproduces), `--maxout m` keeps only subsets
of size at most `m` and is exact for the vertex graph and the outer graphs with at most `m` edges; `g5_graphs.py` has a refinement-based
canonical form and ships the graph lists (indices are tied to `data/cost_table_g5.json`, which balances the array tasks longest-first);
`V_validate.sh` re-checks genus 2-4 against the stored tables and calibrates the machine; `S_array.sh`, `C_collect.sh`,
`submit_genus5.sh main|controls|all`. Cost as measured on the workstation: 39-46 CPU hours per prime at full `--maxout`, heaviest single
inner graph about 0.7 h, so with 100 tasks the wall time is set by that graph. Package review and a 20-graph spot check of the engine:
`natural_form/genus5_package_review/`.

## 4. Files

| file | content |
|---|---|
| `genus3_natural.py` | the factorisation checks and the 42-graph pass for `M = A M_0(tau zhat)` (`--variant nat`: `A = 1`; `--ab`: conifold frame) |
| `genus3_natfix.py` | linear algebra for a possible left dressing `A` (not needed and not used: no dressing is required through genus 3) |
| `g3_nat_{89,61}_nat.json`, `g3_nat_89_nat_ab1o3_m2o7.json` | states: `tau`, `G`, `n1..n3`, pass tables |
| `logs/genus3_natural_{89,61}_nat.log`, `logs/genus3_natural_89_nat_ab1o3_m2o7.log` | logs of record (ALL CHECKS PASSED) |
| `natural_form/lg/` | `LG_MODEL.md`, stationary-phase engine, Gauss-Manin recursion, symbol rows to `zhat^40` over Q |
| `natural_form/theory/` | `UNIFORM_REDUCTION.md`, the gauge over Q from the periods, residue relations, Euler-dilaton checks |
| `natural_form/referee/` | `REFEREE_NAT.md`, from-scratch rows, three primes, `G` over Q (`zhat <= 16`, `u <= 11`) |
| `natural_form/verify_*/`, `natural_form/critic/`, `natural_form/SYNTHESIS.md` | adversarial verifications and the critic's synthesis |
| `fq4_split.py`, `genus4_split.py`, `g4_split_*.json`, `parts/`, `logs/genus4_split_*.log` | genus 4 (and any genus) by edge splitting; regression, runs, controls (s.5) |
| `natural_form/genus4_referee/` | code review, independent re-implementation, genus-5 cost estimate, synthesis |
| `../cluster_genus5_quintic/` | the genus-5 package (s.6a) and the collected `results/*.json` of the Stanage runs; `../cluster_genus5_quintic.tar.gz` |
| `natural_form/genus5_package_review/` | review of the package before submission and a 20-graph spot check of `graph_dp` |
| `natural_form/genus5_cluster/` | verification of the collected genus-5 results: `recheck/RECHECK.md` (own CRT, own `t_c`, census), `local/LOCAL.md` (the `--maxout 1` re-run and 2586 additivity identities here), `controls/CONTROLS.md`, `critic/CRITIC.md` (graph-list completeness, piece assignment, what is overclaimed); then, on the copied-back evidence, `census/CENSUS.md` (the seven tables rebuilt from the part files, provenance, the piece census, additivity for all 4555 at both primes), `frame2/FRAME2.md` (the second frame's point-frame runs, its two-prime total, its additivity, the frame (in)dependence of the selection rule), `indep/INDEP_G5.md` (the module-disjoint re-implementation against 93551 stored pieces) |
| `../g5_evidence.tar.gz` | archival copy of the cluster's `logs/` and `parts/` (46 MB compressed; 1021095 stored pieces), the material every re-check of s.6 runs on |
