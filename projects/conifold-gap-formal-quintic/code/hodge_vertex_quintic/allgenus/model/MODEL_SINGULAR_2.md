# MODEL_SINGULAR.md, part 2: Lemma K, the reduction, hypermaps, answers, open points, files

(Part 1, `MODEL_SINGULAR.md`, has the summary, the setting and the data; this file continues its numbering.)

## 6. Lemma K: the edges of a (G)-gauge avoid the whole vanishing tower

**Lemma K (PROVED).** Let `R(z) = Sum_k R_k z^k` (coefficients in End(H)[[u]]) be symplectic, `R^*(-z) R(z) = 1`, so that the edge
`V(w, z) := (eta^-1 - (R^-1(w) (x) R^-1(z)) eta^-1)/(w + z)` is a power series. Let `phi(w) = Sum_{j >= 0} phi_j w^-j` be a covector
series with no positive power of w such that `(R^-1(w))^T phi(w)` has no positive power of w either (u-adic convergence assumed:
`phi_j -> 0`). For m >= 0 put `Pi_m(f) := Sum_k phi_(k-m)(f_k)` on `f = Sum_k f_k z^k`. Then
 (a) `Pi_m` annihilates the first (and by symmetry the second) slot of V for EVERY m >= 0;
 (b) for the leaf point `-z R^-1(z) 1 = -z + T(z)`: `Pi_m(-z + T) = -[w^(m-1)] <(R^-1(w))^T phi(w), 1>`, which is 0 for m >= 2, equals
     `-<y_0, 1>` for m = 1 (`y_0` the w^0 coefficient of `(R^-1)^T phi`) and `-[w^-1]<(R^-1)^T phi, 1>` for m = 0.
*Proof.* Contract the first slot of `(w + z) V = eta^-1 - (R^-1(w) (x) R^-1(z)) eta^-1` with `phi(w)`: with `C(w, z) := <phi(w), V(w, z)>_1`,
`(w + z) C = eta^-1 phi(w) - R^-1(z) eta^-1 ((R^-1(w))^T phi(w))`. The right side has no positive power of w. Write
`C = Sum_m C_m(z) w^m` (all m in Z; `C_m = Pi_m(V)` for m >= 0, each a convergent sum). The coefficient of `w^(m+1)`, m >= 0, gives
`C_m + z C_(m+1) = 0`, so `C_m = (-z)^j C_(m+j)` for every j >= 0: the coefficient of z^r in `C_m` vanishes (take j > r). Hence
`C_m = 0` for all m >= 0, which is (a). (b) is the definition of `Pi_m` applied to `-w R^-1(w) 1`. QED

**Application.** Take `phi` = the vanishing covector of the model at v(u) (`phi(w) <-> Y^C(-w)`, so that `Pi_m` are the `Pi_i` of s.1)
and R = the gauge (hatted: `R'`). The hypothesis is exactly (G): the gauge carries the model's vanishing solution to a solution with no
positive power of z, `P Y^C = Y^C_Q` (the transpose form follows from symplecticity: `(R^-1(w))^T = eta R(-w) eta^-1`). Symplecticity is
automatic for the actual `R' = R M^-1` (both factors are symplectic R-matrices). With `Y^C_Q = -e_1 - t_c e_0/z + ...` (components
`Y_b` = pairing with `H^b`: the z^0 term has no `b = 0` component, the z^-1 term has `Y_0 = -t_c`) part (b) gives `Pi_1(-z + T) = 0`,
hence `Pi_1(T) = 0` (`Pi_1(-z) = -PiC(tau_0(1)) = 0`), and `Pi_0(-z + T) = t_c` (the quintic's vanishing period; sign as in the
gauge worker's leaf formula). Hence
**edges in `K_0 (x) K_0`, leaf point in `-z + K` with `Pi_0 = t_c`.** In the canonical normalisation only: for `P -> P c(z)`,
c scalar, `P Y^C = c(z) Y^C_Q` has positive powers and the hypothesis fails at the first non-zero `c_k`.

**Verification (prediction; nothing fitted).** `ms_edgeK.py`, actual rows, `R'` from `M = M_0(tau zhat)`, production EdgeLeaf tables:

| rows | prime | frame | shifts k checked (control non-zero) | coefficients per k | non-trivial vanishing sums | leaves Pi_k(T) = 0 |
|---|---|---|---|---|---|---|
| referee, zhat <= 24, u < 18 | 2^89-1, 2^61-1 | (0, 0), (1/3, -2/7) | 0..20 (21-23 vacuous: control zero) | 2160 | 20790 per run | k = 1..20, mod u^18 |
| production, zhat <= 6, u < 6 | 2^89-1, 2^61-1 | (0, 0); (1/3, -2/7) at 2^89-1 | 0..4 (5 vacuous) | 150 | 500 per run | k = 1..2, mod u^6 |

The gauge worker's check (E) is the case k = 0. Control (`--c5`): `R' -> R' exp(c5 z^5)` (still symplectic and a solution of the
gauge equation, Theorem U's freedom `c(z^5)`, but violating (G)), c5 = 1/7, referee rows, 2^89-1, frame (0, 0)
(`ms_edgeK_ref_89_0_0_24_18_c5ctrl.log`): EVERY shift k = 0..20 now fails (e.g. 750 of 2160 coefficients non-zero at k = 0, 180 at
k = 10), and the leaf acquires exactly a `Pi_5`-component (`T_5` fails, `T_1..T_4` still vanish: the extra leaf `c5 z^6 1` has
`Pi_i(tau_6(1)) = PiC(tau_(6-i)(1))`, non-zero for i = 5 only). So the test resolves the normalisation.

## 7. Theorem R: Conjecture N from the Model Lemma

**Theorem R (PROVED, conditional).** Assume (i) the factorisation `R~ = G M_0(tau zhat)` (Corollary F; needs (H)), i.e.
`Omega^Q_u = P(u, z).Omega^LG_(v(u))` with P regular at u = 0 (Theorem U) and symplectic; (ii) (G) for P (canonical normalisation);
(iii) (ML_h) for every h >= 0 (s.5; with `F_0` the n >= 3 Taylor part of `-(5/2)V^2 log V`, `F_1 = -(1/12) log V`). Then for every
g >= 2 every stable graph with at least one edge in the outer graph sum `R'.(M.omega^)` is regular at u = 0, and the vertex graph is
`F_g(t_c) + regular`, `F_g(t_c) = c_g t_c^(2-2g)`. That is Conjecture N.
*Proof.* Givental's formula for `P^ A^LG`: a graph contributes `1/|Aut|` times the contraction, along its edges, of the derivatives
`d^(n_v) Fbar_(h_v)(-z + T)` of the model's ancestor potentials at the leaf point, in the directions of the edge slots. By Lemma K the
slots lie in `K_0` and `-z + T` in `-z + K` with `Pi_0(-z + T) = t_c`. By (ML), for `q in -z + K`,
`Fbar_h(q) = F_h(Pi_0(q)) + Reg_h(q)` with every Taylor coefficient of `Reg_h` at -z (on the homogeneous basis tau_0(H), kappa(X, k)
of K) of the form `c v^e`, e >= 0. Since `Pi_0` is linear and vanishes on `K_0`, `F_h(Pi_0(-z + T + x)) = F_h(t_c)` for `x in K_0`:
every derivative of a vertex in an edge direction only sees `Reg_h`. In a graph with an edge every vertex carries an edge, so the graph
is a contraction of regular edges (power series in u) with derivatives of `Reg_h` at `-z + T`; these are u-adically convergent sums of
regular terms (the leaves are O(u) in hatted variables; in the flat frame each leaf, of psi-degree >= 2, raises the v-exponent by
>= 3/5), with coefficients in the kappa-basis that are polynomial in v (s.1 (i)-(ii)). The vertex graph is
`Fbar_g(-z + T) = F_g(t_c) + Reg_g(-z + T)`. QED
*Remarks.* (1) The argument never uses the D-directions `tau_k(H)`, k >= 1, where the singular structure is genus dependent (s.3, s.4):
by Lemma K they never occur. (2) It explains why "the vertex carries the whole gap" and why the dilaton leaves turn the model coordinate
into t_c (the gauge worker's leaf formula is `Pi_0(-z + T) = t_c`). (3) It explains the universality runs (univ/un_e2): for any self-dual
cyclic system whose Theorem-U gauge is normalised so that the analogue of (G) holds with the system's own vanishing solution, the same
proof gives regular edge graphs and the vertex `F_g(T_sys)`, `T_sys` the system's vanishing coordinate: the pattern observed in
univ/un_e2 (that the perturbed systems' canonical gauges satisfy (G) was not checked here).
(4) It reduces Conjecture N to statements about the LG model `x^4 + v/x` alone (ML) plus the gauge identity (G).

## 8. The H^2-tower is a theorem (hypermaps)

arXiv:2211.12259 (Carlet, van de Leur, Posthuma, Shadrin), Theorem 2 (read: pp. 1-11 of the arXiv PDF): for the Hurwitz Frobenius
manifold H_{0;1,N-1} with the calibration (1.24) of that paper, expanded at the special point `f = p^(N-1) + 1/p`, the descendant
potential restricted to the descendants `t^1_d` of the flat coordinate `t^1 = a_1` (all other descendant times zero) is the generating
function of rooted hypermaps: `log D = Sum_g eps^(2g-2) Sum_{n >= 1} (1/n!) Sum_k RHM_{g; k_1+1, .., k_n+1} prod t^1_(k_i)/(k_i+1)!`,
where `RHM_g` counts rooted genus-g hypermaps with n white polygons of sizes `k_i + 1` and `Sum (k_i + 1)/N` black N-gons.
*Dictionary (LG_MODEL.md s.5.2).* N = 5; `t^1 = a_1` (the position of the simple pole) is our `H^2` (up to a constant), its descendant
tower is `kappa(H^2, k)` (the column `S(z)H^2` of the homogeneous calibration is unique, no resonance; the paper's (1.24) at `beta = 1`
is homogeneous), the slice "all other times zero" is a sub-slice of K, and the special point is v = 1/4 (`a_5 = 4v`).
*Consequence (PROVED given the theorem and this dictionary).* A genus-h correlator of `kappa(H^2, k_i)`, i = 1..n, has v-exponent
`e = 2 - 2h - Sum (1 - 4k_i)/5 = 2 - 2h + 4B - n` with `B = Sum (k_i + 1)/5`, which by Euler's formula `V - 5B + (n + B) = 2 - 2h` is the
number V of vertices of the hypermaps it counts: `e = V >= 1`. Adding m insertions `tau_0(H) = d/dv` gives `v^(V-m) V!/(V-m)!`: regular.
So (ML) holds on `C tau_0(H) + span{kappa(H^2, k)}` for all genera. The engine agrees (pure kappa(H^2) correlators have e >= 1 in all
computed cases; e.g. `<kappa(H^2, 4)>_1 = v^3/32`, e = 3 = V for a white and a black pentagon glued into a torus).
*What remains:* the towers of the unit and of `H^3`, `H^4` (expansions at the pole of order 4 of the superpotential in the
spectral-curve picture of the paper's proof). For them (ML) is VERIFIED in the computed range only.

## 9. Answers to the task items (b)-(d) and the model lemma in edge form

(b) **Two/three-point singular parts.** On K: exactly `F_h^(n)(v) prod PiC(X_i)`, i.e. `F_h` of the LINEAR vanishing coordinate
`Pi_0` (the descendant time t^H_0), with no Hessian at all; this is neither `F_h` of the v-component of the genus-0 point (wrong
dilaton degree; its ancestor version W gives `<tau_2(H^4)>_2 = 0`, false) nor a jet expression (on K no jet and no genus mixing
occurs: `<H^n>_h = F_h^(n)` exactly and every K_0-correlator is regular). The second-difference (Toda) form is not an alternative but
a property of the one-variable function: `Sum_h hbar^(h-1) F_h(V)` (F_0, F_1 included) is the asymptotic expansion of
`log G(1 + N)`, G the Barnes function, at `N = V (-5/hbar)^(1/2)` (since `c_h hbar^(h-1) = B_2h/(2h(2h-2)) (-hbar/5)^(h-1)`), and
`G(N+2) G(N)/G(N+1)^2 = N`: the singular part on K solves the Toda equation in Pi_0 with step `(-hbar/5)^(1/2)` (formal; elementary). Off K (descendants `tau_k(H)`, k >= 1) the singular
parts are genus dependent (s.4); at one point they are `Sum_{i >= 1} Pi_i(X) D_{h,i}` with one number per (h, i) (s.3:
`D_{4,5} = 19/100800000`, `D_{3,5} = 11/241920`). Every singular correlator computed and its remainder: `ms_singlist_full.log`.
(c) **Primary-sector lemma** (genus >= 1 primary correlators vanish unless all insertions are H, and `<H^n>_h = F_h^(n)`): a
consequence of (ML), since `tau_0(H^a) = kappa(H^a, 0) in K_0` for a != 1 and every primary correlator at genus >= 1 has
`e = 2 - 2h - Sum w(a_i, 0) < 0` (all `w(a, 0) > 0`). VERIFIED wherever (ML) is.
(d) **Operators.** On `-z + K` the singular part `F_h(Pi_0(q))` is annihilated by the string vector field `q -> [q/z]_+` (it preserves
K and `Pi_0([q/z]_+) = Pi_1(q) = 0` there), is homogeneous of degree 2 - 2h for the dilaton vector field (Pi_0 is linear) and of the
right weight for the Euler field; it is even invariant under every `q -> [z^-m q]_+` (Pi_m = 0 on K), which are not symmetries of
the full potential for m >= 2. The higher Virasoro operators L_m (m >= 1) were not examined: their quadratic parts pair the
H-direction with the H^2-direction, so they do not act on the two factors separately (no claim).
**Model lemma in edge form (what the universality result needs).** Given (ML): for every bivector `B(psi, psi')` with coefficients
regular at u = 0 and `B in K_0 (x) K_0`, and every leaf vector T regular with `T in K`, all graphs built from model vertices at v(u)
with at least one B-edge are regular, and the leafy vertex is `F_h(v + PiC(T))` + regular. Every gauge in the image of the gauge equation
(Theorem U applied to a self-dual cyclic system with one simple pole) is of the form `P_can exp(Sum_j c_(5j) z^(5j))` (j odd) with
`P_can` the (G)-normalised one; Lemma K puts the edges of `P_can` in `K_0 (x) K_0`, and the vertex becomes `F_h(T_sys)`. For `c != 0` (G)
fails at the first non-zero `c_(5j)` and the edges leave `K_0 (x) K_0` (control in s.6).
**The genus-4 test of the normalisation.** If the gauge is normalised non-canonically, `P = P_can (1 + c_5 z^5 + ..)`, the genus-4
graph sum changes at first order, at its leading order `T^-2` (hatted u^4 x u^-6), by `c_5 Delta_4(5)` with
`Delta_g(m) = <tau_(m+1)(1)>_g + (1/2) Sum_j (-1)^j [<tau_(m-1-j)(e_a) tau_j(e^a)>_(g-1) + Sum_{h1+h2=g} <tau_(m-1-j)(e_a)>_h1 <tau_j(e^a)>_h2]`
(frozen values, `eta^-1 = (1/5) Sum_a H^a (x) H^(3-a)`; leaf `c z^(m+1) 1`, edge `c Sum_j (-1)^j w^(m-1-j) z^j eta^-1`).

**Delta_4(5) (`ms_delta4.log`, VERIFIED, exact).** `Delta_4(5) = -1/8400000 (leaf) + 0 (loop) + 1/8400000 (bridges) = 0`: at its
leading order (hatted `u^-2`) the z^5 normalisation of the gauge is INVISIBLE at genus 4, through two separate cancellations
(`Delta_3(5) = -1/640` has v-exponent 0 and is regular anyway; `Delta_2(5) = 0`). So genus 4 can see the normalisation at most at
`u^-1`, through products with the first u-corrections of the gauge, and the control `--c5` above (which breaks Lemma K at every
shift) is the sharper test. Whether `Delta_4(5) = 0` has a structural reason (e.g. a Mumford-type identity for the model) is open.

## 10. What is not established

* **(ML) itself** beyond the computed range and beyond the H^2-tower: it is the one model-side input of Theorem R. Data: genus 0-4,
  n <= 6/5/4/3/1 (s.5); the genus-2 four-point part was a prediction. A proof route: in the spectral-curve picture (TR allowed for the model), the K-towers are the expansions of
  the TR differentials at the two poles of `x = p^4 + v/p` (unit, H^3, H^4 at the pole of order 4; H^2 at the simple pole, where
  arXiv:2211.12259 gives hypermaps), while the D-tower is the logarithmic (extended) one; regularity of the expansions at the pole of
  order 4 would follow from the regularity of `omega_{g,n}(p)` for `|p|^5 > |v|` that the TR worker reports (VERIFIED there, 13 cases,
  not proved), and the mixed 0/infinity expansions are not covered by either.
* **Lemma K's hypothesis** (G) for the actual quintic is VERIFIED (zhat <= 24, u < 18, two primes, two frames) and argued by the gauge
  worker (monodromy + normalisation), not proved here; it rests on the canonical normalisation `P_11(0, z) = 1` and on (H).
* **Theorem R** was written and checked by one AI agent (this note); no referee has read it. The convergence statements in its proof
  (u-adic, via homogeneity in the flat frame) are argued, not formalised. Genus 0 and 1 vertices need (ML_0), (ML_1) with the
  stated normalisations of F_0, F_1 (data: genus 0 to n = 6, genus 1 to n = 5).
* **The D-directions** (descendants of H): their singular parts are genus dependent and were not classified (no Toda or jet form was
  established for them; at one point they are carried by one number per (h, i), s.3). They are irrelevant for Conjecture N by
  Lemma K, but a complete singular calculus of the model would need them.
* **(3,3)** was computed by a process of the first attempt (same code; finished 18:57); it passes string, dilaton and permutation
  checks against tables computed here (`ms_sd33.log`) and agrees with MVertex6 on all 6125 correlators of its 49 tensors with
  psi-deficit <= 6 (`ms_regress3_33.log`); the 7 <= deficit <= 9 part is covered only by string/dilaton. `amp(4,2)` (28037 graphs) was
  not completed: the first 2000 graphs alone exceeded 23 CPU minutes and the run was stopped (time limit); a profile of ten windows
  (0.02-2 s per graph, single graphs up to about 100 s) puts the whole table at roughly 4 CPU hours with this engine, so genus-4
  two-point data are missing.
* Higher Virasoro operators (task item 2(d)) were not examined.

## 11. Files (all in allgenus/model/; scratch: the amplitude caches, graph lists, part files)

| file | content |
|---|---|
| `MODEL_SINGULAR.md`, `MODEL_SINGULAR_2.md` | this note (two pieces) |
| `ms_rows.py`, `ms_engine.py`, `ms_model.py` | rows/edges/leaves in charge variables, the flat-basis correlator engine, S, PiC, kappa (first attempt; reviewed) |
| `ms_regress.py/.log`, `ms_regress2.py/.log`, `ms_regress3.py`, `ms_regress3_33.log`, `ms_sd33.py/.log` | regressions: record numbers, MVertex6 (mod 2^89-1), string/dilaton/permutations |
| `ms_graphs.py`, `ms_chunk.py`, `ms_chunk_2_4.log` | fast legged-graph generator (validated against fq3_graphs), time-bounded chunked amplitudes; the (2,4) run |
| `ms_regress3_24.log`, `ms_sd24.log`, `ms_ktest_24.log` | validation of amp(2,4) and the genus-2 four-point (ML) test (prediction) |
| `ms_corr.py` | helpers: F_h derivatives, PiC, admissible insertion lists |
| `ms_onept.py/.log`, `ms_twopt.py/.log` | s.3, s.4 |
| `ms_ktest.py/.log`, `ms_singlist.py/.log`, `ms_singlist_full.log` | s.5: (ML) on K; every singular correlator with its remainder |
| `ms_edgeK.py`, `ms_edgeK_*.log` | s.6: Lemma K on the actual R' tables (and the `--c5` control) |
| `ms_delta4.py/.log` | s.9: first-order effect of the z^5 normalisation |
