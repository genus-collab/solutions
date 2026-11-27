# Reduction of Conjecture N (part 2 of 3): Lemma 1 and Theorem R

Continuation of `REDUCTION_DRAFT.md` (s.0-s.2, notation s.1); part 3 is `REDUCTION_DRAFT_3.md`. Same labels and conventions.

## 3. The reduction

### 3.1 Lemma 1 (the transport lemma), in four parts

**Lemma E (edge contraction; PROVED).** Let `R(z) = 1 + O(z)` be symplectic (`R(-z)^T eta R(z) = eta`), `E_R(psi, psi') = (eta^-1 - R^-1(psi)
eta^-1 R^-1(psi')^T)/(psi + psi')` its edge, and `f` in `Hc_-` (vector form). Contracting the covector `Omega(f, .)` (which pairs `tau_m(e)`
with `eta([z^(-m-1)] f(-z), e)`) into the first slot of `E_R` gives the vector

        (Omega(f, .) (x) 1) E_R (psi')  =  - R(psi')^-1 (R f)_{>=0}(psi').

In particular `E_R` annihilates `Omega(f, .)` (in either slot, `E_R` being symmetric) iff `R f` lies in `Hc_-`.
*Proof.* Put `alpha(x) = eta f(-x)` (only negative powers of `x`), so the contraction is `C(psi') = Res_x alpha(x)^T E_R(x, psi') dx`. Expand
`1/(x + psi') = Sum_j (-psi')^j x^(-j-1)` in each of the two terms of the numerator (legitimate: their sum is the power series `E_R`). The
`eta^-1` term gives `Res alpha(x)^T x^(-j-1) = 0`. The other gives `-Sum_j (-psi')^j [x^j](alpha(x)^T R^-1(x)) eta^-1 R^-1(psi')^T`, i.e.
`C^T = -R^-1(psi') eta^-1 beta_{>=0}(-psi')` with `beta(x) = R^-1(x)^T alpha(x)`. Symplecticity gives `R^-1(x)^T = eta R(-x) eta^-1`, so
`beta(x) = eta (R f)(-x)` and `beta_{>=0}(-psi') = eta (R f)_{>=0}(psi')`. QED (Infinite sums converge u-adically in all applications
below: the coefficients of `z^-n` in `Y^C` are `O(v^(n - floor(n/5)))`.) This is Lemma (L) of the gauge worker, VERIFIED there with random
covectors; the proof here uses nothing but the edge formula of arXiv:1303.1043 and symplecticity.
**Corollary E'.** If `P = P_0 R` satisfies (Sp) and (A), then for every `k >= 1` the edges of `P` annihilate `V^(k)` in each slot, because
`(P z^-k J^C)_{>=0} = z^-k (P J^C)_{>=k} = 0`. Hence `E_P` lies in `T_v (x) T_v` (completed, coefficients in `Q[[u]]`, `v = v(u)`).

**Lemma T (translation; PROVED).** For every `f`: `Omega(f, q_T) = Omega(P f, -z phi_0)` (symplecticity). With `f = z^-k J^C`:

        V^(k)(q_T) = (-1)^(k+1) [z^(k-2)] Y'_0(-z) = T_P (k = 1),   0 (k >= 2),

by (A) for the positive powers and the grading for `z^0` (`Y'_0` has only powers `z^-1 mod 5`). Moreover `q_T = -(z/d) phi_0 + T_{>=2}`,
`T_{>=2}` in `z^2 Hc_+`, so `qt := d q_T = -z phi_0 + d T_{>=2}` lies in `Lambda_v`, with `V^(1)(qt) = d T_P = v + s`, and
**`s := d T_P - v` lies in `u^2 Q[[u]]`**: `T_P = v'(0) (P_0)_00(0) u + O(u^2)` (only `P_00(0, z) (-v'(0) e_0/z)` contributes at `u z^-1`)
and `(P_0)_00 (P_0)_33 = 1` by (Sp) at `z^0`. So the leaves move the vanishing coordinate from `v` to `T_P` and keep the vertex point on
`Lambda`: this DERIVES the "dilaton leaves convert the model coordinate into `t_c`" of NATURAL_FORM s.2.

**Lemma F (finiteness; PROVED).** For each stable graph `Gamma`, `Cont_Gamma` is a FINITE sum of products of model correlators at `v(u)`
(Laurent monomials in `u`) and entries of `P_k`, `(P^-1)_k`, `d^(+-1)` (power series in `u`). *Proof.* At a vertex, homogeneity gives
`d^n Fb_h(q_T) = d^(2h-2+n) d^n Fb_h(qt)`; the Taylor expansion of `Fb_h` at the base point in the direction `d T_{>=2}` has at most
`3h - 3 + n` terms (each insertion has psi-degree `>= 2`), and only finitely many psi-degrees of the edge factors survive the dimension
bound. QED This replaces the "u-adic convergence" of the task: in the flat frame there is nothing to converge (`P(0, z) != 1` is harmless).

**Lemma C (Fock-space form; PROVED given the formula `(R^)G(q) = [exp(hbar E_R(d, d)/2) G](R^-1 q)` of arXiv:math/0108100).** If `f` and
`P f` lie in `Hc_-`, then `P^ o (mult. by Omega(f, q)) = (mult. by Omega(P f, q)) o P^`; hence `P^ (Phi(V^(1), V^(2), ..) D) =
Phi(V'^(1), V'^(2), ..) P^ D` for every function `Phi` of the c = 1 functionals (`V'` built from `Y' = P Y^C`). *Proof.* `Omega(f, R^-1 q) =
Omega(R f, q)` and, by Lemma E, `exp(hbar E_R(d, d)/2)` commutes with multiplication by `Omega(f, .)`. QED Not used below; it is the
conceptual content: **the regular gauge transports the whole c = 1 sector of the model to that of the system.**

### 3.2 Theorem R

**Theorem R.** Fix `g >= 2`. Let `P` be flat-regular, graded, `P_0` invertible, with (Sp) and (A) (s.1.4), and assume ML_Lambda(h) for
`0 <= h <= g - 1`. Then (a) `Cont_Gamma` is in `Q[[u]]` for every stable graph `Gamma` of genus `g` with at least one edge.
If moreover ML_Lambda(g) holds, then (b) `Cont_vertex = c_g T_P^(2-2g) + Q[[u]]`, and (c) the ancestor potential of `Omega^P` restricted to its
own slice `Lambda'` (built from `Y' = P Y^C`) satisfies ML_Lambda(g) with `V'^(1)`: the Model Lemma transports.
*Proof.* (1) *Edges.* By Corollary E' write `E_P = Sum Et^(alpha beta) t_alpha (x) t_beta` with `t_alpha` the frames of `T_v` and `Et` the flat
components with `a, b != 1` (Lemma 0); the frame sums converge u-adically (the frame corrections `f_j(v) z^j H` of `frame(k, a)` carry
`v^(w(k,a) - w(j,1))`, while the leading terms die by dimension for large `k`). So every half-edge insertion is a frame.
(2) *Vertices of graphs with edges.* In a connected graph with an edge every vertex carries `n >= 1` half-edges. By Lemmas T, F its factor is
`d^(2h-2+n) Sum_m <t_alpha1 .. t_alphan (d T_{>=2})^m>_h / m!` (finite), and `d T_{>=2} = t' + s phi_1` with `t'` in `T_v` (frames with
`k >= 2`, coefficients in `Q[[u]]`) and `s` in `u^2 Q[[u]]`. Expanding multilinearly, every term is a model correlator with at least one
`T_v` insertion and all other insertions in `T_v (+) Q phi_1`, times power series; by ML_Lambda(h) each such correlator is a non-negative
power of `v(u)`. (Multilinear extension over frames: frames with `k >= 2` have weight `w >= 4/5`, so a correlator with many of them carries a
high power of `v`, and the frame expansion of `t'` converges u-adically inside every correlator; it agrees with the flat expansion because
only finitely many flat slots survive the dimension bound. The expansion never splits a frame into its flat components, so the
directions `tau_k(H)`, `k >= 1`, which occur inside frames, are never isolated.) So every vertex factor is in `Q[[u]]`; the edges are (flat regularity); the contraction is finite (Lemma F).
This proves (a), using ML_Lambda only for the genera `h <= g - 1` that occur at vertices of graphs with edges.
(3) *The vertex graph.* `Cont_vertex = d^(2g-2) Sum_m <(t' + s phi_1)^m>_g / m!`. Terms with a `t'` are regular by ML_Lambda(g); the rest is
`Sum_j s^j <phi_1^j>_g / j! = Sum_j s^j Phi_g^(j)(v)/j! + Q[[u]] = Phi_g(v + s) + Q[[u]]` (convergent: `s^j Phi_g^(j)(v) = O(u^(2-2g+j))`). With
`v + s = d T_P` and `d^(2g-2) Phi_g(d T_P) = c_g T_P^(2-2g)` this is (b). (4) For (c) add legs `w` in the slice directions of `Omega^P`:
the leg factor is `P^-1 w`, and `V^(k)(P^-1 w) = Omega(z^-k P J^C, w) = V'^(k)(w)`, so legs in `T'` become insertions in `T_v` and steps
(1)-(3) apply verbatim. QED

**Hypotheses needed on the Model Lemma (exactly what the proof uses).** (i) Only on `Lambda`: insertions from `T_v` and `tau_0(H)`; nothing
about the descendants `tau_k(H)`, `k >= 1`, alone. (ii) Only the vanishing direction: singular correlators with at least one `T_v` insertion
vanish; no formula for any non-zero singular part except the gap. (iii) Graphs with edges: genera `h <= g - 1` only, including `h = 0, 1`
(genus-0 and genus-1 vertices occur). (iv) Vertex graph: the gap `<tau_0(H)^m>_g = Phi_g^(m)(v)` + regular (equivalently `F_g = c_g v^(2-2g)`)
and the `n >= 1` statement at `h = g`. (v) For each `g`, only the finitely many correlators that reach hatted order `<= 2g - 3`.

**Comparison with the conifold Theorem R** (`../../../hodge_vertex_conifold/proof_status_mixed_scale.md` s.3, `leaf_proof/leaf_structure_proof.md`).
There: (i) homogeneity (dilaton) moves each vertex to the frozen CohFT at a shifted point; (ii) Lemma C expands every leg, loop and edge at
order `t^k` in geometric leaves `lam_j` of degree `D <= k` (proved order by order from the hypergeometric singular row); (iii) `M_irr(h)`:
the one-leaf functionals `Xi_h(R)` vanish for `D(R) <= 2h - 3`. Here: (i) the same (homogeneity of the model, exact dilaton factor `d`);
(ii) is replaced by Corollary E' and Lemma T, which are EXACT consequences of (G) - the edges lie in `T_v (x) T_v` identically in `u`, so
no degree bound and no order-by-order analysis is needed; (iii) is ML_Lambda, with "negative `v`-exponent" in the role of `D <= 2h - 3`.
As there, only the vanishing direction of the lemma is used, the single-vertex graph needs the gap separately (Theorem VI there;
ML_Lambda(g) at `n = 0` here, i.e. `F_g = c_g v^(2-2g)`, VERIFIED `g <= 5`), and graphs with an edge need the lemma only for `h <= g - 1`.
**Independent derivation.** The model worker (`../model/MODEL_SINGULAR.md`, `MODEL_SINGULAR_2.md`, s.6-s.7) reached the same reduction in
parallel: its `K` is `T_v (+) Q phi_1` (descendant times `t^H_k = 0`, `k >= 1`), its Lemma K is Corollary E' + Lemma T (with a different,
generating-function proof covering all shifts at once), its Theorem R is Theorem R(a), (b) for the canonical gauge. What this note adds:
general dressings (hypotheses (Sp), (A), flat regularity, and their necessity in s.3.4), part (c) and Prop. I (transport and gauge
invariance of the Model Lemma), Lemma F (finiteness in the flat frame), the exact normalisation (`s = d T_P - v = O(u^2)` from (Sp),
`d = tau^(1/4) xhat`, so that `v/d = vhat`, the project's `u tau^(-3/2)/xhat`), and the second test of ML_Lambda (s.2.4).

### 3.3 Consequences

**Corollary N (conditional all-genus theorem).** Conjecture N holds for every `g >= 2` if: Theorem U (PROVED in UNIFORM_REDUCTION s.3,
AI-audited), (H) and (N) [the actual LP R-matrix equals the canonical gauge applied to the model: `c^act = 1`] and ML_Lambda(h) for all `h`.
*Proof.* Under (H), (N) the outer graph sum of the project is the graph sum of `P_can.Omega^mod_{v(u)}` (s.1.5); `P_can` satisfies the
hypotheses of Theorem R (s.1.4 (i)) with `T_P = t_c`. QED
**Corollary N_9 (genus <= 9).** (H) and (N) are VERIFIED in the window `zhat <= 24`, `u < 18` (gauge worker, s.6 of GAUGE_IDENTITY), which
contains every order of `R'` used by a genus-g graph sum for `g <= 9`. Hence, for `g <= 9`, **Conjecture N(g) follows from ML_Lambda(h),
`h <= g`**, a finite list of model correlators for each `g` (only decorations of hatted order `<= 2g - 3` matter, NATURAL_FORM s.6.1).
Genus 6-9 (beyond any graph-sum computation in the project) are thereby reduced to computations in the model alone.
**Corollary U (universality; PROVED given ML_Lambda).** For every self-dual cyclic one-pole system `A` with equal products and its canonical
gauge, every edge graph is regular and the vertex graph is `c_g (T_A/lambda)^(2-2g)` + regular, `1/lambda = P_11(0, z)`. Here the residues
of the two poles may differ (`r = -1/lambda^2`, `rt = -1`): (R1) still gives column 1 of `P(0, z)` equal to `P_11(0, z) e_1`, so Theorem G
reads `P Y^C = P_11(0, z) Y^C_A`; Prop. S becomes `P(-z)^T eta P(z) = (rt/r) P_11(0, z) P_11(0, -z) eta` (the entry `(2, 1)` now carries
`P_22(0, z) = (rt/r) P_11(0, z)`), so the constant normalisation `P_11(0, z) = 1/lambda` is symplectic and annihilating, `Y' = Y^C_A/lambda`,
`T_P = T_A/lambda`. This is E2/E4 of `../univ/`, including the factor `lambda = 3/2` of the cases with `f(0) = 3/2`.
**Corollary L (leaf bookkeeping; PROVED given ML_Lambda(g)).** In the project's normalisation the leafless vertex is `c_g vhat^(2-2g)`; the
leaves turn it into `c_g t_c^(2-2g)` + regular, whose first term is `-(g-1) c_g (2/vhat^(2g-1))(t_c - vhat)`; with `t_c - vhat = u^2/30 + O(u^3)`
(VERIFIED, `gi_hat_table.log`) this is the observed `-(g-1) c_g/15 u^(3-2g)` (UNIFORM_REDUCTION s.6.2, listed there as not derived).

### 3.4 Comparison with the experiments (each hypothesis is needed in some example)

| experiment (`../univ/`, `../gauge/`) | (Sp) | (A) | flat-regular | prediction of Theorem R | observed |
|---|---|---|---|---|---|
| canonical gauges of perturbed systems (E2 genus 2, 3; E4 c1, c2 genus 4, all 378 graphs; c2 has `lambda = 3/2`) | yes | yes | yes | regular, vertex `c_g (T/lambda)^(2-2g)` | regular, vertex as predicted |
| A1-A4: annihilating NON-gauge dressings, genus 2, 3 (E3); one at genus 4, all 378 graphs (`un_e3_g4.log`) | yes | yes | yes | regular, vertex `c_g T_eff^(2-2g)` | regular; `T_eff = -[z^-1 e_0] D Y^C_Q` |
| first order, genus 2, 3 (E3 lin) | inf. | inf. | yes | response 0 on ker Ann | rank(Resp on ker Ann) = 0 |
| B1, B2: annihilation violated | yes | no | yes | no prediction | singular |
| C1a-C1c: flat-irregular | yes | yes | no | no prediction | singular (C1c regular at genus 2 only) |
| stored normalisation `P_11(0,z) = 1 - c_5 z^5 - ..` (E4 b) | mod `z^10` | no (`c_5 z^5 e_1` at `u^0`) | yes | no prediction | singular at genus 4, `u^4`, one-edge graphs |

Every observation is consistent with Theorem R. The failing families show that neither (A) (B1, B2, and the stored normalisation, which
violates (A) through `P Y^C = c(z) Y^C_Q`) nor flat regularity (C1a-C1c) can be dropped; (Sp) is needed for the edge formula itself (without it
the numerator of `E_R` does not vanish at `psi = -psi'` and the edge is not a power series). The flatness of the gauge (the u-ODE) is not a
hypothesis: A1-A4 are not gauges (the `z^0` coefficient of `D Y^C_Q` is not constant) and are regular, as predicted.
