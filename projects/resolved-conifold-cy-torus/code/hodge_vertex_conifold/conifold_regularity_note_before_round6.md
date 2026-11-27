# Why every stable graph with an edge is regular at the conifold point (Hodge-vertex frame B): structure, theorem, and what remains

Resolved conifold x C^2, fibrewise torus alpha0 = 0, alpha1 = eps1, alpha2 = eps2; Hodge-vertex frame B (weights w = (c, -eps1, -eps2), c = eps1 + eps2, R' = R_Z h_w^-1 with h_w = exp of the Stirling tail of log B(eps1/z, eps2/z)); conifold point Q = 1, t = log Q. Conventions of hvc_lib.m / rc_lib.m: a stable graph Gamma contributes

    (1/|Aut Gamma|) Sum_{colourings c: V -> {1,2}} Sum_{decorations} prod_v Delta_{c_v}^{g_v - 1} Int_{M_{g_v, n_v + m_v}} Hodge_w psi^a prod_{dilaton leaves} Dil_{c_v}(d) prod_e E_{c c'}(a_e, b_e),

idempotent 1 = the root U_1 -> 1 of U^2 - S U + P (the critical point that runs into the logarithmic singularity of the superpotential), idempotent 2 = the root U_2 -> -eps2/eps1. Below eps2 = 1, x = eps1, so c = 1 + x. Everything marked [exact] was computed exactly in Q and x from the stored R-matrix factors (z^9, hvc_asymptotics_g4.wls, hvc_normalised_R.wls, hvc_vertex_identity2.wls, hvc_percolouring.wls); [data] refers to the genus 2-4 graph sums.

## 1. Facts about the building blocks at Q = 1

F1 [exact; proved to all orders, section 7b] Delta_1 = -x c^3 / t^2 (1 + O(t)) (double pole), Delta_2 = x^3/c + O(t) (regular, nonzero).

F2 [exact to z^9; the first two statements proved to all orders in z, section 7b] The identity row of R' at the singular idempotent tends to the geometric series
        r_{0,1}(z; t) = c/(c + z) + O(t),
    and the H-row to x/(c+z) + (c/t)(1 + O(t)) (the pole is the canonical coordinate L_1 = -(x U_1 + 1) Q/(1 - Q)). In the idempotent basis:
        R'(z) e_1 = e_1 - z/(c + z) e_2 + O(t),     R'(z) e_2 = h_w(z)^-2 e_2 + O(t),
    and in the normalised basis (e~_i = Sqrt[Delta_i] e_i):
        R~'(z; t) = diag(1, h_w(z)^-2) + O(t),   with R~'_11 - 1 = O(t^2).
    In the point frame (R_Z) the same limit is diag(h_w, h_w^-1); in mode A it is diag(h_w^2, 1). So frame B is exactly the frame in which the singular idempotent is frozen (R-matrix 1), the Hodge class sitting there undressed.
    Derivation of the leading term (sketch; written out to all orders in z in leaf_proof/leaf_structure_proof.md s.3, see section 7b): with U = 1 - v and Q = e^t the superpotential is eps1 log(v - t) + eps2 log v + (analytic), so the oscillatory integral at the singular critical point is t^{c/z + 1} B(eps1/z + 1, eps2/z + 1) (1 + O(t)); Stirling gives the tail of log B(eps1/z, eps2/z), i.e. h_w(z), times the elementary factor z/(c + z) from B(a+1, b+1) = B(a, b) ab/((a+b)(a+b+1)); normalising the leading coefficient to 1 gives h_w(z) c/(c+z) in the point frame and c/(c+z) in frame B. At the regular critical point the integral is a Beta integral with the opposite sign of the tail, hence h_w^-1 (point frame) and h_w^-2 (frame B).

F3 [exact to z^9; proved to all orders, section 7b] In frame B all edge factors E_ij(m, n) and all dilaton factors Dil_i(d) are regular at t = 0, while in the point frame and in frame A the edge factor E_11 has a double pole (E_22 is regular in every frame: O(t^0) in the point frame, O(t^2) in frame A). The limits of E_11 and E_12 are of rank one, that of E_22 is not:
        E_11(z, w) = -(x^3/2) phi(z) phi(w),   E_12(z, w) = phi(z) chi(w),   phi(z) = 1/(c - z),   chi(w) = -x^3 h_w(w)^2/(c + w),
        E_22(z, w) = (x^3/c) [1 - h_w(z)^2 h_w(w)^2]/(z + w),
    [Erratum 2026-10-31: this line used to read "E_11 = E_22 = -(x^3/2) phi(z) phi(w)", and E_22 was said to have a double pole in the point frame
    and in frame A. Both came from reading the conjugate-encoded EdN[2,2] at U2 (rc_lib.m stores the idempotent-2 data as cconj of the idempotent-1
    data, so every stored factor must be read at U1; read at U2, EdN[2,2] is the E_11 table). Corrected in hvc_vertex_identity.wls line 14 and
    hvc_asymptotics*.wls line 24; leaf_proof/leaf_structure_proof.md Proposition 8.6(c), Remark 8.8; hvc_leaf_proof_check2.wls (M14, Q2). The graph
    sums are not affected (they are cconj-invariant), nor are Theorems A, B, C, which use E_22 only through its regularity.]
    and the dilaton-leaf series at the singular idempotent is T(psi) = -psi^2/(c - psi) = psi - c psi phi(psi).

F4 [exact, genus 2-4] Vertex functional. For the Hodge class Hodge_h = Lambda^vee_h(c) Lambda^vee_h(-eps1) Lambda^vee_h(-eps2) / (c (-eps1)(-eps2))^h put
        Z[h, n] = Sum_{m >= 0} (1/m!) Int_{M_{h, n+m}} Hodge_h prod_{i=1}^n phi(psi_i) prod_{j=1}^m T(psi_{n+j}).
    Then (-x c^3)^{h-1} Z[h, 0] = gap_h (the t^-(2h-2) coefficient of P_h Li_{3-2h}(Q)/(eps1 eps2)) for h = 2, 3, 4, and
        Z[2,1] = Z[2,2] = Z[3,1] = Z[3,2] = Z[3,3] = Z[4,1] = Z[4,2] = 0,
    whereas Z[1, n] != 0 (genus 1 carries no Delta pole, so nothing is needed there).

F5 [data, genus 3] Inside a stable graph with a genus-2 vertex the decoration sum is regular for every colouring separately, but neither the sum over dilaton patterns at fixed psi-assignment nor the sum over psi-assignments at fixed dilaton pattern is regular: the cancellation is a joint kappa/psi phenomenon (Mumford-type), localised at the singular vertex and its colouring.

## 2. Theorem A (decorated graphs)

Assume F1 and F3. Then every decorated graph of a stable graph all of whose vertices have genus <= 1 is regular at Q = 1, and the pole order of any decorated graph is at most Sum_v max(2 g_v - 2, 0).

Proof. The term is a product of factors that are regular at t = 0 (F3; the Hodge integrals are numbers) times prod_v Delta_{c_v}^{g_v - 1}. By F1, Delta_2^{g_v-1} is regular, Delta_1^{g_v - 1} has a pole of order 2(g_v - 1) for g_v >= 2, equals 1 for g_v = 1 and has a zero of order 2 for g_v = 0. Summing over colourings and decorations cannot raise the order. QED.

This is exactly the observed pattern at genus 2, 3 (local P2 and conifold) and 4: regular decorated graphs outside the stable graphs with a vertex of genus >= 2, and the bound Sum_v (2 g_v - 2)_+. It also explains the point frame: there E_11 ~ t^-2 adds 2 per edge, and Sum_v 2(g_v - 1) + 2E = 2(g - 1), so every decorated graph has order exactly 2g - 2, as observed.

## 3. Theorem B (leading pole of a stable graph with edges)

Assume F1, F3, F4 (F4 in the form Z[h, n] = 0 for all h >= 2, n >= 1). Let Gamma be a stable graph with at least one edge and fix a colouring with set S of singular vertices (colour 1). The coefficient of the maximal power t^{-2 Sum_{v in S}(g_v - 1)} of the decoration sum is, by F3 (rank-one edges at t = 0, singular half-edges carrying phi) and F1,
        (1/|Aut|) prod_{v in S} (-x c^3)^{g_v - 1} Z[g_v, n_v] x (factors of the regular vertices and of the chi-legs),
because at t = 0 the dilaton leaves and the singular half-edges decouple vertex by vertex. Since Gamma is connected with an edge, n_v >= 1 for every v, so every factor Z[g_v, n_v] with g_v >= 2 vanishes by F4 and the leading pole cancels, for every colouring. For the single vertex graph (n = 0) the same formula gives Delta_1^{g-1} Z[g, 0] = gap_g (F4), i.e. the vertex graph carries the whole gap at the leading order.

## 3a. Theorem C: the vertex identity (proof by the string equation)

For h >= 2 and all n >= 1, Z[h, n] = 0. Moreover Z[1, n] = (n-1)! A_1 / c^(2n-1) with A_1 = Int_{M_{1,1}} Hodge_1 = -(1/24) Sum_a 1/w_a (both statements checked against the exact values for n <= 3).

Proof. Put t_s(z) = c s phi(z) + T(z) = (c s - z^2)/(c - z) and consider the generating series of all correlators of the Hodge CohFT with the leaf t_s at every point,
        Phi_h(s) = Sum_{N >= 0} (1/N!) Int_{M_{h,N}} Hodge_h prod_{i=1}^N t_s(psi_i) = Sum_{n >= 0} (c s)^n Z[h, n] / n!
(expand (c s phi + T)^N and collect the n legs phi: the binomial N!/(n! m!) against 1/N! gives 1/(n! m!) and the factor (cs)^n). Differentiating in s inserts one extra point with the leaf d t_s/ds = c phi = 1 + z phi, so
        Phi_h'(s) = D_{1} Phi + D_{z phi} Phi,
where D_f denotes the insertion of a point carrying f. The insertion of the constant leaf 1 is a tau_0 insertion; the string equation (valid for Hodge integrals, since the lambda classes are pulled back along the forgetful map) replaces it by the string shift of every other leaf, (D t_s)(z) = (t_s(z) - t_s(0))/z = (s - z)/(c - z), plus the exceptional term delta_{h,1} A_1 (for h >= 2 the integral of a class pulled back from M_h over M_{h,1} vanishes; the genus-0 exceptional term does not concern us). Hence
        Phi_h'(s) = D_{D t_s + z phi} Phi + delta_{h,1} A_1,   and   D t_s + z phi = (s - z)/(c - z) + z/(c - z) = s phi = (s/c) d t_s/ds,
so Phi_h'(s) = (s/c) Phi_h'(s) + delta_{h,1} A_1, i.e. (c - s) Phi_h'(s) = c delta_{h,1} A_1. For h >= 2 this forces Phi_h' = 0 as a formal power series in s, i.e. Z[h, n] = 0 for every n >= 1; for h = 1, Phi_1'(s) = A_1 Sum_k (s/c)^k gives Z[1, n] = (n-1)! A_1/c^{2n-1}. QED

(Three independent attempts of a multi-agent proof search all arrived at this argument; the identities d t_s/ds = 1 + z phi and D t_s + z phi = s phi were also checked symbolically.) Nothing beyond the string equation is used, in particular no property of the Hodge class: the identity holds for any CohFT with a flat unit whose leaves are the geometric series phi = 1/(c - psi) and the dilaton leaves T = -psi^2/(c - psi) = psi - c psi phi. The gap Z[h, 0] = Phi_h(0) is the one quantity the argument leaves free.

Consequently Theorem B holds unconditionally (given F1, F3): the maximal pole of every stable graph with at least one edge cancels, colouring by colouring, at every genus.

## 4. What is proven, what is verified, what is open

- Status 2026-10-31, later (supersedes everything below where it differs): F1, F2, F3 and Lemma C (the leaf structure) are PROVED to all orders
  (section 7b, leaf_proof/leaf_structure_proof.md). With Theorem C and Theorem R (proof_status_mixed_scale.md s.3) the all-genus theorem is
  conditional only on Lemma A (M_irr(h) for h >= 5); it holds unconditionally through genus 4 (genus 2 without computer input).
- Status 2026-10-31 (supersedes the lines below where they differ): Conjecture M is established at genus 2, 3, 4 for every tau and every number of
  leaves, and the leaf route then gives Conjecture V and the theorem through genus 4 independently of the direct graph sums. The all-genus
  statement is reduced to Lemmas A-D of section 7a.

- Proven given F1, F3: Theorem A (all genera). F1 and F3 are properties of the explicit conifold R-matrix; they are verified exactly to z^9 (all orders relevant through genus 4), and the leading behaviour F2 has the Beta-integral derivation above. A complete proof of F3 to all orders in z needs the exact t -> 0 structure of the Wick expansion at the degenerate critical point (the pole of the H-row cancels precisely because of the geometric form c/(c + z) of the identity row: the residue of L_1 r_{0,1}[m] - Q r_{0,1}[m-1] dlogDelta_1 is (1/t)(-1/c)^{m-1}(c (-1/c) + 1) = 0).
- Proven (Theorem C, string equation): the vertex identity F4 for all h >= 2, n >= 1; hence Theorem B, the vanishing of the maximal pole of every stable graph with an edge, colouring by colouring, at every genus (given F1, F3).
- Observed structure of the corrections (genus-4 factors, z^9): the O(t) term of the identity row at the singular idempotent is t (1-x)^2 Sum_{m>=1} (-1)^{m-1}(2^m - 1) z^m / c^{m+2}, a rational function of z with poles at z = -c and z = -c/2; the O(t^2) term has poles at -c, -c/2, -c/3 (Pade). So the full t-dependence is organised by geometric series in z/c with scaled arguments, which suggests that the subleading cancellations are again string/dilaton-type identities for leaves of the form 1/(c - k psi); this is the natural next step towards the all-order statement.
- Open: the subleading poles t^{-(2 Sum (g_v - 1)) + k}, k >= 1, of a stable graph with edges. The exact-in-Q structures do not persist beyond t = 0 (r_{0,1} is geometric only at t = 0: r[2] - r[1]^2 = O(t); the rank-one form of E_11 fails at O(t^2)), so the cancellation of the lower orders mixes the O(t) corrections of the leaves, edges and Delta_1 with the t = 0 vertex data; F5 shows it is a per-colouring, joint kappa/psi cancellation. The data (genus 2, 3 for both theories, genus 4 for local P2 and the conifold) show full regularity, so the conjecture is:

  Conjecture. In frame B every stable graph with at least one edge is regular at the conifold point for every colouring, and the single vertex graph has polar part exactly the gap, at all genera.

## 5. Remarks

- The role of the triple class: with the geometric rank-4 class (eps1, eps2, -eps1, -eps2) one would have h = 1 and R' = R_Z, i.e. the point frame, where R~ -> diag(h_w, h_w^-1) and every graph is singular. The local-P2 triple class is the twist that removes exactly the Stirling tail produced by the degenerate critical point, R~' -> diag(1, h_w^-2).
- The cancellation of the maximal pole requires only the t = 0 data, which are universal (c/(c+z), phi, T, the Hodge class): this is why the leading polar coefficients of corresponding decorated graphs agree between local P2 and the conifold up to the universal factor (-3)^(g-1) (t_c^2 = -3 t^2), and the next order up to (-3)^(g-1)(-(eps1+eps2)^2/(9(eps1-eps2)^2)).

## 5a. Reduction Lemma (all genera): Conjecture V implies the theorem

Conjecture V (vertex form). For every h >= 2, n >= 1, every choice of n external legs, each a column of an edge term,
f_{j,b}(psi) = Sum_a E_{1j}(a, b) psi^a (j in {1, 2} the idempotent of the other end, b >= 0 its psi-power), and every number of loops
(pairings E_11(a, a') of two points of the vertex), the exact singular vertex functional
    V_h^(t) = Delta_1(t)^{h-1} Sum_{m >= 0} (1/m!) Int_{M_{h, n + 2 loops + m}} Hodge_w prod_i f_i(psi_i) prod_loops E_11(psi, psi') prod_j T_t(psi_j)
is O(t^0); for n = 0 and no loops it is gap_h t^{2-2h} + O(1). (Verified at genus 2 for all configurations tried, section 6.)

Lemma. Assume F1, F3 (Delta_2, all edge factors E^rc and all dilaton factors Dil^rc regular at t = 0) and Conjecture V. Then for every
stable graph Gamma with at least one edge and every colouring, the sum over all decorations is regular at t = 0; the single-vertex graph
has polar part exactly the gap.

Proof. Fix Gamma and the colouring; let v_1, ..., v_s be the singular vertices (colour 1) of genus >= 2 (singular vertices of genus 1
carry no pole, those of genus 0 a zero of order 2; regular vertices are regular). In a stable graph the psi-powers on the half-edges of a
vertex are bounded by its dimension, so for every edge e only a finite block of the matrix E_{c c'}(a, a') enters. For an edge between
two singular vertices the block of E_11 is a finite symmetric matrix and can be written as Sum_alpha f_alpha[a] f_alpha[a'] with vectors
f_alpha in its column space, i.e. in the span of the columns f_{1,b}; for an edge from a singular to a regular vertex the singular side is
the column f_{2,b} with b the power on the regular side; loops at a singular vertex are E_11-pairings; dilaton leaves and the Hodge
integral belong to the vertex. Hence the decoration sum of Gamma is a finite sum of products, over the singular vertices, of functionals
V_{g_v}^(t) with legs in the admissible family and loops, times regular factors (the regular vertices' sums, Delta_2 powers, E_22, the
regular sides of E_12, Dil_2). Each V_{g_v}^(t) is O(t^0) by Conjecture V (n_v >= 1 because Gamma is connected with an edge), so the
product is regular. For the single vertex graph the only term is V_g^(t) with n = 0, i.e. the gap. QED

Erratum (2026-10-31). The factorisation step fails for an E_11 edge between two singular vertices of genus >= 2. Conjecture V applies to legs
that are columns f_{1,b}, so the coefficients c_alpha in f_alpha = E c_alpha would have to be regular at t = 0. That is impossible as soon as the
block has two rows: E = E C E with C regular makes E C idempotent, so rank E(t) = trace(E(t) C(t)) would be constant near t = 0, whereas the block
has rank 1 at t = 0 (F3) and rank 2 for small t != 0 (det of the block over psi-orders {0, 1} = (1-x)^4 x^6 t^2/(18 c^10) + O(t^3); exact,
hvc_leaf_proof_check2.wls (Q1)). Such blocks first occur at genus 4 (two genus-2 vertices joined by an edge, both of colour 1). The valid
reduction is Theorem R (proof_status_mixed_scale.md s.3; leaf_proof/leaf_structure_proof.md s.13, Remark 13.1): it expands rho^2 E_11 in the
t-independent leaves lam_i(psi) lam_l(psi') with coefficients regular in t and splits the leaf degree between the two ends, so no factorisation
is needed. Conjecture V for a single singular vertex with legs f_{j,b} and loops is unaffected (it is steps (i)-(iii) of Theorem R).

So the all-genus statement is a property of the singular idempotent alone (Conjecture V, in the leaf form of Theorem R; see the erratum): its metric normalisation Delta_1(t), its
dilaton series T_t and the singular sides of its edge terms, all of which are determined by the thimble integral of the degenerate
critical point (section 6) and by the Hodge CohFT.

## 6. Pathway to the subleading poles (proposal, 2026-10-30)

Fact G [exact to z^9]: the t^k correction of the singular-idempotent column of R' (both flat components) is a rational function of z with
denominator exactly prod_{j=1}^{k+1} (c + j z) (hvc_tcorrections.wls; k = 1, 2, 3). [2026-10-31: the upper bound (poles only at -c/j, j <= k+1,
all simple) is proved for all k, section 7b; "exactly" holds for generic x and is checked for k <= 12, but fails at special x, e.g. the t^1
correction of rho r_{0,1} vanishes at x = 1.] This is the fingerprint of a Gauss hypergeometric function:
the thimble integral of the degenerate critical point is an Euler integral. With U = 1 - v, Q = e^t, tau = 1 - e^-t, a = eps1/z, b = eps2/z
the superpotential is eps1 log(v - tau) + eps2 log v - eps2 log(1 - v) and the measure dU/U = -dv/(1 - v), so
    I(z; t) = Int_0^tau (tau - v)^a v^b (1 - v)^{-b-1} dv = tau^{a+b+1} B(a+1, b+1) 2F1(b+1, b+1; a+b+2; tau),
and the k-th tau-coefficient (b+1)_k^2/(k! (a+b+2)_k) = prod_{j=1}^k (eps2 + j z)^2 / (k! z^k prod_{j=1}^k (c + (j+1) z)) has exactly the poles
of Fact G (the extra (c + z) comes from the Beta factor, B(a+1,b+1) = B(a,b) ab/((a+b)(a+b+1))). The R'-column is the z-asymptotics of I
divided by h_w(z) and the classical saddle factors; its regularity at tau = 0 (F2, F3) is the statement that the thimble is the solution of
exponent 0 of the Gauss equation tau(1-tau) F'' + [gamma - (alpha+beta+1) tau] F' - alpha beta F = 0, alpha = beta = eps2/z + 1, gamma = c/z + 2,
at the regular singular point tau = 0 (the other exponent, 1 - gamma = -c/z - 1, is the one that would carry a t^-1/z pole).

Proposed steps.
1. Make F1-F3 theorems: derive the singular column exactly from the Euler integral / Gauss equation (asymptotics in z of 2F1 with large
   parameters), identify h_w(z) c/(c+z) as the tau = 0 value, and read off the exact tau-expansion of r_{0,1}, r_{1,1}, hence of Delta_1,
   the dilaton series T_t and the singular legs of E_11, E_12: at order tau^k they are finite sums of products of scaled geometric series
   1/(c - j psi) (j <= k+1) with polynomial prefactors. Check against the stored factors (z^9) at every step.
2. Vertex formulation of the all-order statement: for a singular vertex of genus h with n half-edges, let V_h^(t)(f_1, ..., f_n) be
   Delta_1(t)^{h-1} Sum_m 1/m! < prod_i (L_t f_i)(psi_i) prod_j T_t(psi_j) >_h, where L_t is the exact singular-leg operator (the
   z-dependence of the singular side of E_11, E_12 applied to arbitrary external series f_i). Conjecture V: V_h^(t) = O(t^0) for n >= 1 and
   gap_h t^{2-2h} + O(1) for n = 0. By multilinearity and the finite rank of the tau^k corrections of the edges this implies the regularity
   of every stable graph with an edge for every colouring (F5 says the cancellation is local to the singular vertices and needs the
   psi-assignments and dilaton patterns jointly, exactly as V requires). Test V numerically first (genus 2, 3 with the exact leaves to
   O(t^{2h-2}), external series = monomials psi^p); the minimal one-idempotent cases are the genus-g vertex with loops (graph 2 at genus 3).
3. Mechanism: the string equation acts diagonally on scaled geometric leaves, D[1/(c - j psi)] = (j/c)/(c - j psi), and the dilaton equation
   handles the polynomial prefactors; all singular leaf data come from ONE solution of the Gauss equation, so the generating function
   Psi_h(s; tau) = Sum_N 1/N! < prod t_{s,tau}(psi_i) >_h (t_{s,tau} = the exact leaves with a leg-counting parameter s, as in Theorem C) should
   satisfy a closed first-order system in s coupled to the Gauss equation in tau. Theorem C is its tau = 0 slice (d/ds Psi_h = 0 for h >= 2);
   the aim is a homogeneous linear ODE in tau for d/ds Psi_h with zero initial condition, giving d/ds Psi_h = 0 to all orders. The Hodge class
   and the Calabi-Yau condition, unused at tau = 0, must enter here through Mumford's relation (the vertex graph's polar part being exactly
   the gap with no lower-order terms is itself an all-order statement about T_t alone and is the first thing to prove this way).
   TEST AT GENUS 2 (2026-10-30, hvc_vertex_conjecture_g2.wls, exact in Q and x, factors to z^9): Conjecture V holds in every case tried.
   Bare vertex: {t^-2, t^-1} = {gap_2, 0} exactly. One external leg through E_11 with external power b = 0..4, through E_12 with
   b = 0..4; two legs (E_11 E_11 with powers {0,0}, {0,1}, {1,2}, {2,3}; E_12 E_12 with {0,0}, {1,2}; E_11 E_12 with {1,0}); one loop;
   loop + leg (E_12, b = 0; E_11, b = 1); two loops; two loops + leg (E_12, b = 0): in all these the t^-2 and t^-1 coefficients vanish
   identically in x, i.e. the exact singular genus-2 vertex is regular as a multilinear functional of the external data. Two structural
   by-products: (i) the t^0 coefficient of a leg through E_11 depends on the external power b only through the factor c^-b (so the
   regular part of the E_11 leg operator is still rank one, with the geometric external leaf), while through E_12 it does not; (ii) the
   two-loop case shows the cancellation survives several singular half-edges on the same vertex, as needed for the genus-g vertex with
   loops. The genus-3 vertex (Delta_1^2, orders t^-4..t^-1) is the next test.
4. Frame refinement as an alternative: the frozen structure holds in the normalised basis only up to O(t) off-diagonal terms (R~'_12, R~'_21 =
   O(t)); a t-dependent symplectic frame change M(z; t) close to 1 that removes them exactly would freeze the singular idempotent to all
   orders and move the whole polar part into the vertex by construction, at the price of a t-dependent (but explicit, hypergeometric)
   vertex class on the singular idempotent; the group-action property guarantees the same totals.
5. Transfer to local P2: the leading Beta structure (h_w c/(c+z)) is the same, which is why the two theories agree at the two leading orders;
   the thimble of the local P2 conifold point obeys a third-order Picard-Fuchs equation, so its tau-corrections differ (hence the 609 odd
   leading orders there against none for the conifold); the same programme applies with that equation in place of Gauss's.

## 7. Mixed-scale vanishing (Conjecture M), 2026-10-30

Leaves lam_j(psi) = 1/(c - j psi) (j >= 0), dilaton leaves T = -psi^2/(c - psi), Hodge weights (c, -x, -1), c = 1 + x, and
Xi_h(j_1..j_r) := Sum_m 1/m! < lam_{j_1} ... lam_{j_r} T^m >_h = d_{lam_j1}...d_{lam_jr} F_h(q_0),  q_0 = c - c^2/(c - z).

**Conjecture M.** For h >= 2 and integers j_i >= 0: Xi_h(j) = 0 <=> 1 <= J = Sum j_i <= 2h - 3 + r (equivalently J >= 1 and D = Sum (j_i - 1) <= 2h - 3).

Evidence (exact; hvc_multiprobe.wls, hvc_multiprobe_g4.wls): both directions on all lattice points with J <= 2h - 1 + r, for genus 2 (r = 1, 2, 3),
genus 3 (r = 1, 2, 3) and genus 4 (r = 1, 2). For r = 1, Xi_h(j) = prod_{l=1}^{2h-2}(j - l) R_h(j) with deg_j R_h = h; a single leaf 1/(cp - psi) of
generic scale does not vanish (hvc_mixed_scale.wls). For r = 2 the polynomial does not factor through J.

Which classes (hvc_which_class.wls, single probe, genus 2 and 3): the zeros j = 1..2h-2 hold for the CY triple class with the leaf scale equal to a
weight (whole family (1, tau, -1-tau)) and for the single class Lambda^vee(1) (its tau -> infinity limit, extra zero at j = 2h). They fail for the trivial
class, for Lambda(1)Lambda(w2) with generic w2, for a non-CY triple and for a leaf scale that is not a weight. Only j = 1 (Theorem C) is universal.

Proved: (a) string and dilaton give Xi(1, j) = (J/c^2) Xi(j) and Xi(0, j) = ((D - 2h + 2)/c^2) Xi(j). (b) At x -> 0 (class lambda_h) the lambda_g formula
reduces Xi to genus-0 KW, Xi_lead = b_h N! [s^N] e^{Js}(s/(e^s - 1))^{N+1} = b_h (J-1)...(J-N), N = 2h - 3 + r (Norlund B_N^{(N+1)}), which is
Res_w (1+w)^{J-1}/w^{N+1} = C(J-1, N); checked against b_2, b_3, b_4.

Role: by Fact G the order-t^k corrections only involve leaves with D <= k, and Conjecture V at genus h needs total D <= 2h - 3, which is
exactly the range of Conjecture M. Together with a leaf-structure lemma for the corrections (no polynomial remainders, J >= 1), Conjecture M
gives the all-genus theorem by Theorem R (proof_status_mixed_scale.md s.3; the Reduction Lemma of s.5a has a gap, see its erratum). The
leaf-structure lemma is now proved (section 7b). Proof attempts for M: proof_status_mixed_scale.md.

## 7a. Status after the proof workflow (2026-10-31)

Details and proofs: proof_status_mixed_scale.md (synthesis), proof_M_hurwitz.md, proof_M_residue.md, referee_*.md, leaf_structure.md.
Normalisation c = 1: weights (1, tau, -1-tau), tau = -x/(1+x), sigma = tau(1+tau), N = 2h-3+r.

Proved for all genera.
- Kappa form: the dilaton leaves sum to K = exp(-Sum kappa_k/k), and Lambda^vee(1) K = c(T(-log)) (log tangent bundle, Chiodo/Mumford GRR).
- Components: Xi_h(R) = Sum_{k=0}^{h} (-1/sigma)^k Xi^(k)_h(R) with tau-independent Xi^(k). The lambda_g component is
  Xi^(h)_h(R) = b_h (J-1)(J-2)...(J-N). This gives the non-vanishing direction of Conjecture M for all h, r (also for non-integer scales),
  and shows that vanishing at h values of tau with distinct non-zero sigma implies vanishing for all tau.
- String/dilaton reduction: the vanishing direction of M at genus h is equivalent to M_irr(h): Xi_h(R) = 0 for all multisets R of
  integers >= 2 with D(R) <= 2h-3.
- Single-leaf answer: a leaf of scale c/j vanishes in the Theorem-C configuration iff j in {1, ..., 2h-2} (j = 0 vanishes trivially for n >= 2).
  The proved non-vanishing covers j >= 2h-1 and non-integer j; the vanishing is the statement M_irr(h) for r = 1.

M_irr(h) established.
- h = 2: proved without computer (no sections of omega_log^{-1}(2 x_1) on any stable curve, plus the lambda_g component).
- h = 3: exact lattice data (hvc_multiprobe.wls).
- h = 4: exact, all seven points with r = 3, 4, 5 that were missing (hvc_multiprobe_g4b.wls, 185 s); the six control points at D = 6 are non-zero.
- tau = 1 (and tau = -2), all genera and all r, modulo the dictionary P2: global-involution residue argument (proof_M_hurwitz.md s.3).
  tau = -1/2: the same for leaves with all parts odd.

Dictionary P2 (TR form). Xi_h(mu; tau) = sigma^(h-1) [x^mu] F_{h,r} / prod binom(mu_i (tau+1) - 1, mu_i - 1), where F_{h,r} are the
Eynard-Orantin correlators of C_tau: x = v^tau/(1+v)^(tau+1), y = -v (TR in X = log x; one branch point v = tau), expanded at v = oo.
Exact check with hvc_tr_tau.wls (general-tau exact TR): tau = 1, 2, 3, genus 2-3, r <= 3, 51 lattice points per tau, no mismatch;
tau = 2, 3 have no global involution. A proof route (BKMP for C^3 plus the variational formula) is sketched in referee_hurwitz_1.md s.4,
not written out.

Closed one-leaf formula (6.1): sigma^h Xi_h(j) = prod_{l=1}^{2h-2} (j-l) [lam^(2h)] S(e1 lam)^(1+j e2) S(e2 lam)^(1+j e1),
S(y) = (y/2)/sinh(y/2), e1 = -tau, e2 = 1+tau. Identical to the exact Xi_h(j) in j and x for h = 2, 3, 4 (hvc_check61.wls);
agrees with exact TR for h = 5, 6 at tau = 1 and -1/2. Unproved for h >= 5; it covers only r = 1.

Leaf structure (Lemma C): the order-t^k parts of rho q_t - q_0, of the legs and of the loops are combinations of lam_1..lam_{k+1}
(loops: pairs with D <= k), with no lam_0 and no polynomial remainder, verified exactly for k <= 5 (every coefficient that occurs
through genus 4). With M_irr(h) for h <= 4 this gives Conjecture V and the theorem through genus 4 by the leaf route. [Now proved to all
orders, section 7b.]

Open, as the remaining lemmas of the all-genus proof.
- A: M_irr(h) for h >= 5.
- B: P2 written out as a proof.
- C: the leaf structure to all orders (natural route: the exact hypergeometric thimble of the degenerate critical point). [PROVED, section 7b;
  so only A remains for the all-genus theorem, and B, D are one possible route to A.]
- D: the TR degree bound at integer tau >= 2 (total x-degree >= 2g-2+n at v = oo). Mod B, D at tau = 1..h gives A at genus h.
  Exact TR confirms D at tau = 2, 3 through genus 3.
Refuted heuristic: if F_g were invariant under polynomial shears X -> X + G(y), with the new branch points escaping to v = oo, the first
non-vanishing moment Res y^(2g-1) omega_{g,1} would be independent of tau. Exact values (hvc_moment_test.wls, genus 2-5) depend on tau,
so the logarithmic point v = oo must contribute and the shear route needs a correction term.

## 7b. Lemma C proved to all orders, without topological recursion (2026-10-31)

Full text: leaf_proof/leaf_structure_proof.md (assembled from the section files leaf_proof/S1_thimble.md, S2_metric_qde.md, S3_leaves_loops.md,
S4_legs_closed.md, each refereed with verdict "sound" and revised, then checked by a completeness critic; plan in leaf_proof/outline.md).
No topological recursion and no spectral curves are used.

Theorem LC (Lemma C with F1-F3, to all orders in t and in z). Frame B, rho = (Delta_1 t^2/(-x c^3))^(-1/2), lam_j(psi) = 1/(c - j psi).
 (F1) Delta_1 = -x c^3/(t^2 rho^2) with rho = tau/(t alg) analytic, rho(0) = 1; Delta_2 = x^3 delta s*^2 analytic, Delta_2(0) = x^3/c.
 (F2) r_{0,1}(z; 0) = c/(c + z), rho r_{1,1} = c/t + x/(c + z) + O(t).
 (F3) every Dil_i(d) and every E_ij(m, n) is analytic at t = 0.
 (L1) [t^0](rho q_t) = q_0 = c - c^2 lam_1 and [t^k](rho q_t) in span(lam_1, ..., lam_{k+1}) for k >= 1.
 (L2) [t^k](rho f_{j,b}) in span(lam_1, ..., lam_{k+1}) for j = 1, 2, every b and every k >= 0; E_21(b, a) = E_12(a, b).
 (L3) [t^k](rho^2 E_11(psi, psi')) in span{lam_i(psi) lam_l(psi') : i, l >= 1, (i-1) + (l-1) <= k}, for every k >= 0.
The bounds are sharp for generic x. Here tau = 1 - e^-t, delta = (c^2 - 4 x tau)^(1/2), s* = 2/(c + delta), alg = (delta/c)^(1/2)(c + delta)/(2c).

How it is proved.
1. The exact singular row. With U = 1 - tau s the singular thimble is Euler's integral Int s^(1/z) (1 - s)^(x/z) (1 - tau s)^(-1/z - 1) ds; its
   critical point s* stays nondegenerate at tau = 0, and the code's Wick row (in y = log U, degenerate at t = 0) equals the formal Laplace
   expansion in s (coordinate invariance, and homogeneity of the Wick monomials). The expansion satisfies Gauss's equation by a pointwise
   amplitude identity and the vanishing of total derivatives (no boundary terms), Gauss's recursion fixes it up to its tau = 0 value, and that
   value is the Beta-function constant h_w(z) c/(c + z) (Laplace's method plus Stirling's series, the one analytic input). Hence
       rho r_{0,1}(z; t) = (c/(c + z)) (tau/t) 2F1(1/z + 1, 1/z + 1; c/z + 2; tau) Exp[-dphi(tau)/z],   dphi = phi_tau(s*) - phi_0(1/c),
   phi_tau(s) = log s + x log(1 - s) - log(1 - tau s). The regularity at z = 0 of every tau-coefficient also has a purely formal proof (eikonal
   and transport equations).
2. Dilaton leaves (L1): each tau^n coefficient is a rational function of z with simple poles only at -c/j, j <= n + 1, no pole at z = 0 and the
   right limit at z = oo; partial fractions give L1 and the absence of lam_0.
3. Metric and QDE: Delta_1 = x^2 W_2(U_1) gives rho = tau/(t alg); the second QDE row holds at both idempotents (at U_2 by the same Laplace
   calculus, coefficientwise at t != 0); <v_i(w), v_j(-w)> = delta_ij Delta_i (K_1 = 1 is a Wronskian identity for Gauss's equation).
4. Loops (L3) by an exact flow, integrated from t = 0 (integration constant 0 by K_1 = 1):
       d/dt [t^2 rho^2 E_11(psi, psi')] = Theta(t) d/dt(t g(psi)) d/dt(t g(psi')),   Theta = -x^3 tau s*^2 = O(t),   g(psi) = rho r_{0,1}(-psi).
5. Legs (L2): the E_12 legs are rational in psi, vanish at oo and are power series at 0; E_21 = E_12^T by the orthogonality of the rows.

Closed forms (all k): rho q_t = c + Sum_j A_j(t) lam_j with
    A_j = (-1)^j C_j t^-1 zeta^j 2F1(j x/c, j x/c; j; tau),   C_j = P_{j-1}/(j! (j-1)! c^(2j-4)),   P_k = prod_{i=1}^k (k + 1 - i c)^2,   zeta = tau e^(dphi/c);
Sum_j A_j = -c^2 and Sum_j j A_j = -c^2 rho; loops e_il(t) = (i l/(c^2 t^2)) Int_0^t Theta (t A_i)' (t A_l)'; the top coefficient (-1)^(k+1) C_{k+1}
of leaf_structure.md (now proved for every k); E_12 at t = 0 is lam_1(psi) chi(w), chi(w) = -x^3 h_w(w)^2/(c + w).

Computer confirmations (none used in the proof): Mathematica hvc_leaf_proof_check.wls, hvc_pf_gauss.wls, hvc_leaf_closed_forms.wls (stored
z^9 factors; closed form to t^12), and the checks requested by the proof, hvc_leaf_proof_check2.wls, _check3.wls, _check4.wls (all passed;
among them the factors rebuilt from scratch to z^12 at x = 2, 3/5, E_11 built from the closed form alone for m + n <= 14 and t^8, and
Proposition 3.15 through z^15); Perl/BigRat section and referee scripts in leaf_proof/ and leaf_proof/referee_scripts/.

Errata found on the way (applied): F3 (E_22 is not of rank one, and has no double pole in any frame); hvc_vertex_identity.wls line 14 and
hvc_asymptotics*.wls line 24 (EdN[2,2] was read at U2); the Reduction Lemma of s.5a (erratum there: Theorem R is the valid reduction).

Consequence. With Theorem C, the string/dilaton reduction to M_irr and Theorem R (proof_status_mixed_scale.md s.3), every coloured stable
graph of genus g with at least one edge is regular at the conifold point, and the single-vertex graph carries exactly the gap, provided
M_irr(h) holds for 2 <= h <= g. M_irr(2) is proved and M_irr(3), M_irr(4) are exact computations, so the theorem holds through genus 4 (genus 2
without any computer input, genus 3 and 4 computer-assisted, independently of the direct graph sums); in every genus it is conditional only
on Lemma A, M_irr(h) for h >= 5. Lemmas B and D (the topological-recursion route to A) are not needed for the reduction.
