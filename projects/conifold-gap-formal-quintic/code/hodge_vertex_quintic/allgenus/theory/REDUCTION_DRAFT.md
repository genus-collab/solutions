# Reduction of Conjecture N to one lemma on the model: skeleton of an all-genus proof (part 1 of 3)

Date 2026-11-13. Folder `hodge_vertex_quintic/allgenus/theory/` (task key `theory`). AI-derived and AI-checked only; no human has read it.
Part 2 (`REDUCTION_DRAFT_2.md`): s.3, Lemma 1 and Theorem R with proofs, consequences, comparison with the experiments.
Part 3 (`REDUCTION_DRAFT_3.md`): s.4 deformation route, s.5 other candidate statements, s.6 ranked Stage-2 routes, s.7 gaps, s.8 files.
Labels: PROVED (complete argument written here, or a quoted theorem with hypotheses checked), VERIFIED (exact computation; orders,
primes, frames stated), CONJECTURAL, REFUTED. Inputs read today: `../gauge/GAUGE_IDENTITY.md` (Theorem G, Prop. S, actual normalisation),
`../univ/UNIVERSALITY.md` and its logs (E3, E4), `../model/` logs. Conventions: `NATURAL_FORM.md`, `natural_form/theory/UNIFORM_REDUCTION.md`.

## 0. Summary

1. **Theorem R (PROVED, s.3).** Let `P` be any flat-regular, Z/5-graded, symplectic dressing of the model at `v(u)` that **annihilates the
   vanishing period** (`P Y^C` has no positive power of `z`). If the model satisfies the **Model Lemma ML_Lambda** (s.2) for all genera
   `h <= g`, then in the genus-g graph sum of `P.Omega^mod_{v(u)}` every stable graph with an edge is a power series in `u`, and the vertex
   graph is `c_g T_P^(2-2g)` + power series, `T_P = -[z^-1](P Y^C)_0` (the transported vanishing coordinate; `= t_c` for the quintic).
   The gauge equation (flatness in `u`) is NOT used. Conversely each of the three hypotheses on `P` is needed in examples (s.3.4).
2. **Consequences.** Conjecture N for all `g` follows from Theorem U, (H), the normalisation (N) of the actual R-matrix (both VERIFIED
   through `z^24`, enough for every graph sum of genus `<= 9`, gauge worker) and ML_Lambda. So **Conjecture N for genus `g <= 9` is reduced
   to a finite statement about the model alone** (ML_Lambda(h), `h <= g`). The theorem also explains, with no further input, all of
   E2/E3/E4 of the univ worker (s.3.4) and derives the "dilaton-leaf bookkeeping" `-(g-1) c_g/15` that NATURAL_FORM s.2 lists as observed.
3. **The Model Lemma ML_Lambda (CONJECTURAL; VERIFIED on 1013 singular correlator lists at genus 0-3, two primes, s.2.4).** Let
   `V^(k)(q) = Omega(z^-k J^C, q)` (k >= 1) be the
   linear "c = 1 functionals" of the vanishing cycle, `T_v` their common kernel and `Lambda_v = {q : V^(k)(q) = 0, k >= 2}`, an affine space
   through the base point. Then on `Lambda_v` the genus-h ancestor potential of the model is `Phi_h(V^(1))` + regular, `Phi_h(x) = c_h x^(2-2h)`
   (`h >= 2`), `-(1/12) log x` (`h = 1`), `-(5/2) x^2 log x` (`h = 0`). Equivalently: every model correlator with at least one insertion
   from `T_v` (and any number of `tau_0(H)`) whose `v`-exponent is negative vanishes. The singular structure off `Lambda_v` (where the
   string equation forces Hessian terms, e.g. `<tau_1(H) tau_0(1)>_h = <H>_h`) is never used.
4. **Independent confirmation.** The model worker reached the same reduction in parallel (`../model/MODEL_SINGULAR*.md`: its `K` is
   `T_v (+) Q phi_1`, its Lemma K is Corollary E' + Lemma T, its Theorem R is Theorem R(a, b)); its exact-over-Q engine confirms ML_Lambda
   on 1019 singular correlators at genus 0-4, PROVES it on the `H^2`-tower for all genera (hypermaps, arXiv:2211.12259 Thm 2), and
   REFUTES every genus-independent singular structure on all of `Hc_+` (the Hessians off `Lambda` depend on `h`): the restriction to
   `Lambda` is essential, not a convenience.
5. **Answers to the task items.** (1) Lemma 1 = Lemmas E, T, F, C of s.3.1 (PROVED; the edge lemma is an exact formula, finiteness replaces
   u-adic convergence). (2) The deformation route reduces to the infinitesimal form of Theorem R; the Hessian term that worried the plan
   lives off `Lambda` and is killed on `Lambda` by annihilation of the WHOLE flag `V^(k)`, k >= 1, which the gauge equation provides but
   which any annihilating dressing also has; E3 (annihilation alone suffices) is exactly what Theorem R predicts (s.4). ML_Lambda is an
   invariant of the gauge class (Prop. I, s.4.3). (3) Toda / TR + x-y swap / Fock sheaf / Lie algebra: s.5. (4) s.3.
   (5) Ranked routes: s.6; the most promising is
   TR on the model curve `x = p^4 + v/p`, `y = p`, for which a closed formula for all `omega_{g,n}` is in the literature (y unramified,
   x rational: arXiv:2212.00320), to be combined with arXiv:1509.06954 (superpotential curve = CohFT).

## 1. Setting

**1.1 Spaces.** `H = Q^5`, flat basis `phi_a = H^a` (`a` in Z/5), `eta(phi_a, phi_b) = 5 delta_{a+b = 3}`, `c = (0, -1, 3, 2, 1)`. Loop space
`Hc = H((z^-1))` (coefficients in `Q[[u]]` where needed), `Omega(f, g) = Res_z eta(f(-z), g(z)) dz`, `Hc_+ = H[[z]]` (`tau_k(phi_a) = z^k phi_a`),
`Hc_- = z^-1 H[[z^-1]]`. Solutions of the Gauss-Manin systems are written with LOWERED components `Y_a = eta(phi_a, J)` as in the project;
`<Y, n> := Sum_a Y_a n^a = eta(J, n)` for `n = Sum n^a phi_a`. Since `eta^2 = 25`, statements for `P` acting on `Y` and for `eta^-1 P eta`
acting on `J` are equivalent; we do not distinguish them.

**1.2 The model.** `Omega^mod_v`: the CohFT of the LG model `x^4 + v/x`, form `dx` (LG_MODEL.md), on its residue line; its ancestor
correlators `<prod tau_ki(phi_ai)>_h(v)` are those of the frozen vertex theory `M_0.omega^` rescaled by homogeneity (PROVED, LG_MODEL s.3,
UNIFORM_REDUCTION s.1): `<..>_h(v) = C v^e`, `e = 2 - 2h + Sum_i w(k_i, a_i)`, `w(k, a) = (4k + c_a - 4)/5`, and `C = 0` unless `e` is an
integer (Z/5 selection). A correlator is **regular** if `e >= 0`. Ancestor potential: `Fb_h(v; q) = Sum_n <(q + z phi_0)^n>_h/n!`, base point
`q = -z phi_0`; dilaton: `Fb_h(lambda q) = lambda^(2-2h) Fb_h(q)` (`h != 1`), `Fb_1(lambda q) = Fb_1(q) - (5/24) log lambda`.
Values on record (model logs, VERIFIED): `<HHH>_0 = -5/v`, `F_1 = -(1/12) log v`, `F_h = c_h v^(2-2h)` with `c_h = |B_2h|/(2h(2h-2) 5^(h-1))`
for `h <= 5` (at `h = 4, 5` through the frozen point-frame totals, NATURAL_FORM s.5, s.6).

**1.3 The vanishing solutions and the c = 1 functionals.** `Y^C(v; z)` (model; residue at the pole) and `Y^C_A(u; z)` (any self-dual
cyclic one-pole system `A`; `Y^C_Q` for the quintic) as in GAUGE_IDENTITY s.2 (Lemmas V, Q, PROVED there): both lie in `Q[z^-1][[u]]^5`,
equal `-e_1` at `u = 0`, `Y^C_b` has only powers `z^-n` with `n = c_b + 1 mod 5`, and `[z^-1] Y^C_0 = -v`, `[z^-1] Y^C_{A,0} = -T_A`
(`T_Q = t_c`). For `k >= 1` define the linear functions on `Hc_+`

        V^(k)(q) := (-1)^k [z^(k-1)] <Y^C(-z), q(z)> = Omega(z^-k J^C, q),     J^C = eta^-1 Y^C,  z^-k J^C in Hc_-.

From the shape of `Y^C` (PROVED): `V^(1)(-z phi_0) = v`, `V^(k)(-z phi_0) = 0` (k >= 2), `V^(1)(phi_1) = 1`, `V^(k)(phi_1) = 0` (k >= 2).
`T_v := intersection of ker V^(k), k >= 1 = { n in Hc_+ : <Y^C(-z), n(z)> in z^-1 Q[[z^-1]] }`.

**Lemma 0 (frames; PROVED).** `Hc_+ = T_v (+) Q[[z]] phi_1`. For every `(k, a)` with `a != 1` there is a unique
`frame(k, a) = z^k phi_a + Sum_{j<k} f_j z^j phi_1` in `T_v` (`f_j` homogeneous monomials in `v`); the frame coordinates of `n` in `T_v`
are its flat `phi_a`-components, `a != 1`. *Proof.* `<Y^C(-z), z^j phi_1> = z^j Y_1(-z) = -z^j + O(z^(j-5))`: eliminate non-negative powers
from the top (`th_core.proj_T`); a combination `Sum s_j z^j phi_1` in `T_v` has `s_j = 0` (look at its top power). QED
Examples (`v = 1`, `th_core`): `frame(1,0) = z + H` (at general `v`: `z phi_0 + v phi_1`, the Euler-dilaton vector of UNIFORM_REDUCTION
s.6.4), `frame(2,4) = z^2 H^4 + H/2`, `frame(3,3) = z^3 H^3 + H/6`, `frame(4,2) = z^4 H^2 - H/24`, `frame(k,a) = z^k H^a` for `(0,0),
(0,2), (0,3), (0,4), (1,2), (1,3), (1,4), (2,3)`. There is no frame with `a = 1`: the descendants `tau_k(H)`, `k >= 1`, are transverse to `T_v`.
**The slice.** `Lambda_v := {q in -z phi_0 + Hc_+ : V^(k)(q) = 0 for all k >= 2} = -z phi_0 + (T_v (+) Q phi_1)` (PROVED from Lemma 0), and
`V^(1)(-z phi_0 + n + s phi_1) = v + s` for `n` in `T_v`: on `Lambda_v` the only coordinate seen by `V^(1)` is the residue direction `H`.

**1.4 Dressings.** `P(u, z) = Sum_k z^k P_k(u)`, `P_k` in `gl_5(Q[[u]])` (**flat-regular**), graded (`P_k = Sum_b p_{k,b} E_{b+k,b}`), `P_0`
invertible, acting on lowered components. Hypotheses: **(Sp)** `P(-z)^T eta P(z) = eta`; **(A)** `Y' := P Y^C(v(u))` has no positive power of
`z`. Put `T_P := -[z^-1] Y'_0` and `d := (P_0)_{33}` (the factor by which `P_0` rescales the unit in vector form; `d(0) = 1` for the gauges).
*Examples.* (i) The canonical Theorem-U gauge from the pulled-back model to any self-dual cyclic one-pole system `A`: flat-regular and
graded (Theorem U), (Sp) (Prop. S of GAUGE_IDENTITY), (A) with `Y' = Y^C_A`, `T_P = T_A` (Theorem G there); PROVED given Theorem U.
(ii) `D P_A` with `D` flat-regular, graded, symplectic and `(D Y^C_A)_{>0} = 0`: the univ worker's non-gauge cases A1-A4.

**1.5 The graph sum.** `Omega^P := P.Omega^mod_{v(u)}`: the R-matrix action (arXiv:1303.1043, s.2: edges `(eta^-1 - R^-1(psi) eta^-1
R^-1(psi')^T)/(psi + psi')`, legs `R^-1(psi)`, translation `z(1 - R^-1(z) 1)`) of `R = P_0^-1 P`, preceded by the constant frame change
`P_0`. For leg-free graphs, `F^P_g = Sum_Gamma |Aut|^-1 Cont_Gamma`, where `Cont_Gamma` contracts the edge bivectors `E_P` with the vertex
tensors `d^(n_v) Fb_(h_v)(q_T)`, `q_T := P^-1(-z phi_0)` (vector form; the dilaton part of the translation included). The project's outer
graph sum `R'.(M.omega^)` (NATURAL_FORM s.5-s.6) is this sum graph by graph in hatted variables (`R' = Psi G` is the hatted form of `P`;
hatting and `Psi` are frame changes, which do not change the graph decomposition; the dilaton factor `d^(2h-2+n)` of each vertex is
absorbed in the vertex TFT `omega^`, the `xhat` of the project).

## 2. The Model Lemma

**ML_Lambda(h) (CONJECTURAL).** For all `n, m >= 0` (`n + m >= 3` if `h = 0`, `>= 1` if `h = 1`) and homogeneous `e_1, .., e_n` in `T_v`:

        <e_1 .. e_n tau_0(H)^m>_h(v)  -  delta_{n,0} Phi_h^(m)(v)    is a regular monomial (zero if its exponent is negative),

`Phi_h(x) = c_h x^(2-2h)` (`h >= 2`), `-(1/12) log x` (`h = 1`), `-(5/2) x^2 log x` (`h = 0`; normalisation `<HHH>_0 = -5/v`).
**Function form (equivalent, by Taylor expansion along `Lambda_v`).** `Fb_h |_(Lambda_v) = Phi_h o V^(1) + R_h` with `R_h` regular (its Taylor
coefficients at the base point along `T_v (+) Q phi_1` are regular monomials; for `h = 0, 1` only orders `>= 3`, `>= 1` are compared).
The `n = 0` part is the gap: VERIFIED for `h <= 5` (s.1.2). The content is the `n >= 1` part.

**2.1 Consistency (PROVED).** (a) *String.* The string vector field `l_{-1} q = [q/z]_+` satisfies `l_{-1} V^(k) = -V^(k+1)`; so it preserves
`Lambda_v` and `V^(1)` on it, as it must (`Fb_h` is string invariant for `h >= 1`). Off `Lambda_v` no function of `V^(1)` alone can be the
singular part: `<tau_1(H) tau_0(1)>_h = <H>_h != 0` although both insertions are killed by `V^(1)`; `tau_1(H)` is not in `T_v (+) Q phi_1`
(`V^(2)(tau_1 H) = [z^0] Y_1(-z) = -1`). (b) *Dilaton.* `Lambda_v` is a cone, `V^(1)` is linear, so `Phi_h o V^(1)` has the dilaton weight of `Fb_h`.
(c) *Euler-dilaton.* `frame(1,0) = z phi_0 + v phi_1` in `T_v`; ML_Lambda gives `<frame(1,0)>_h = (2h-2) F_h + v F_h' = 0` (`h >= 2`); at `h = 1`
the exponent is `0`, so the anomaly `<frame(1,0)>_1 = 5/24 - 1/12 = 1/8` is allowed (it is regular): this is why the genus-1 anomaly of
Proposition E is harmless. (d) *Primary sector.* `phi_0, phi_2, phi_3, phi_4` lie in `T_v`; ML_Lambda contains the primary-sector lemma
(genus `>= 1` primary correlators with a non-`H` insertion vanish).

**2.2 Interpretation.** `Sigma_C := Q((z^-1)) J^C` is isotropic (`Y^C(-z)^T eta^-1 Y^C(z)` is constant in `u` by flatness and vanishes at
`u = 0` because `Y^C(0) = -e_1`, `(eta^-1)_11 = 0`; PROVED). The `V^(k)` are the coordinates dual to its negative half: `V^(1)` is the
residue (the coupling of the c = 1 theory), `V^(k)`, `k >= 2`, its descendants ("tachyon couplings"). ML_Lambda says: with all
tachyon couplings off, the singular part is the pure c = 1 free energy of `V^(1)` (cf. hep-th/9506122: c = 1 at self-dual radius = the
conifold). A natural strengthening, NOT needed: `Fb_h = F^(c=1)_h(V^(1), V^(2), ..; string coordinate) + regular` on all of `Hc_+`.

**2.3 What is known about the model (from LG_MODEL.md).** Frobenius manifold `H_{0;1,4}` (form `dp`), conformal, semisimple off `v = 0`;
`F_0` has the single logarithmic term in the residue coordinate (`c_555 = 1/a_5`); its descendant potential is a tau function of the
extended rationally 4-constrained KP hierarchy (arXiv:2211.12259, QUOTED there, not used).

**2.4 Tests (VERIFIED; `th_mla.py` on `th_core.Model` = MVertex6 on the project symbol, exact mod 2^89-1).** All lists of frames and `tau_0(H)`
with integer negative exponent, frames `frame(k, a)`, `k <= 3h + 3`, psi-deficit `<= 6` (lists outside: skipped and counted):

| genus | insertions | singular lists tested | vanish | control C1: same lists, frames -> leading monomials, non-zero | control C2: regular lists non-zero |
|---|---|---|---|---|---|
| 0 | <= 5 | 43 | 43 | 5 | 269 of 365 |
| 1 | <= 5 | 418 | 418 | 64 | 1264 of 1429 |
| 2 | <= 4 | 425 (102 skipped) | 425 | 165 | 740 of 878 |
| 3 | <= 3 | 127 (262 skipped) | 127 | 88 | 172 of 212 |

So ML_Lambda(0)-(3) hold on every tested list (1013 singular lists); the corrections `f_j z^j H` of the frames are essential (C1), and the
frames are not trivially annihilated (C2; e.g. the regular `<frame(1,0)>_1 = 1/8`, `<frame(3,3)>_2 = 3/3200`). Logs: `th_mla_89.log`;
`th_mla_61.log` (2^61-1: ALL PASS with identical counts). Genus 3 is limited by the psi-deficit `<= 6` of MVertex6.
Independent evidence: the model worker's exact engine (`../model/MODEL_SINGULAR.md` s.5: genus 0-4, 3987 correlators on its basis of `K`,
1019 of them singular, all vanish; `<H^n>_h = F_h^(n)` in 16 cases; the one-point law holds exactly on `K` and fails off it, first at
`<tau_5(H)>_3`), its proof on the `H^2`-tower (s.8 there), and E3 of the univ worker (s.3.4: annihilating non-gauge dressings are regular;
by Theorem R this tests ML_Lambda on the realised insertions). Both engines are AI-written; ours uses the production vertex code MVertex6,
theirs is exact over Q and was regressed against MVertex6 on 6182 flat correlators.
