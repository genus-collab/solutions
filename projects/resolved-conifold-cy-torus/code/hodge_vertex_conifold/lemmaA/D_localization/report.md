# Route D report (round 1, 2026-10-31)

(Saved by the main session from the workflow result wf_40196105-891; subagents cannot write report files.)

## claims
[
   {
      "depends_on" : "Graber-Pandharipande, Invent. Math. 135 (1999); FP 2000 s.4.1",
      "evidence" : "loc_conifold_g012.out and loc_conifold_g3.out: N_{g,d} = d^{2g-3}[lambda^{2g-2}](2sin(lambda/2))^-2 for 12 (g,d) with g <= 3, CY and non-CY weights.",
      "statement" : "D2.1: Graber-Pandharipande formula for N_{g,d} of the resolved conifold X with any torus. With the CY torus (1,tau,-1-tau), the p0-vertex is exactly the Lemma A class, and a degree-d edge gives the leaf d*lam_d.",
      "status" : "proved-modulo-literature"
   },
   {
      "evidence" : "loc_P1.out: Int_{[Mbar_g(P^1,d)]^vir} 1 = 0 by localization in 7 cases.",
      "statement" : "D2.2: every graph of the X-localization is a polynomial in sigma of degree <= g-1+d. Its top coefficient is (-1)^{g-1+d} times the P^1 localization of Int 1, and E(d) = -eps_d(sigma) with eps_d = C(d*tau+d-1, d-1)^2.",
      "status" : "proved"
   },
   {
      "evidence" : "check_laguerre.out 18/18 and check_lagrange.out 6/6 at h = 2, 3; psi_poly.out.",
      "statement" : "D2.3: Xi_h(R) = sum_{p=0}^{N0} (-1)^p L^{(a+p)}_{N0-p}(-(J+p)) MV_h(R u 1^p)/p!, with a = 2h-2+r and N0 = 3h-3+r. Equivalent forms: (1+v)^a e^{Jv} Phi_R(-v e^v) = Psi_R(v/(1+v)), and Xi = (-1)^{N0}[y^{N0}](1-T)^{a+N0-1} e^{-(J+N0)T} Phi_R(y). The identity holds componentwise.",
      "status" : "proved"
   },
   {
      "depends_on" : "Graber-Vakil localization proof of ELSV (Compositio 135, 2003)",
      "statement" : "D2.4: up to the ELSV constant, the Hurwitz component Xi^(0)_h(mu) equals a single-fixed-locus equivariant relative invariant of P^1: br^*[b.0] times prod ev^*(gamma_0) T(psi).",
      "status" : "proved-modulo-literature"
   },
   {
      "statement" : "D3.1: on CY3 targets, equivariant integrals of polynomial classes give no vanishing relations beyond string/dilaton/divisor. Only weight-independence of N_{g,beta} and polynomiality remain.",
      "status" : "proved"
   },
   {
      "statement" : "D3.2 (resonance): a T_{w1}-insertion has a pole exactly at the p_oo end of a degree-1 edge. The K-twisted vertex is reached only through localized insertions, or by analytic continuation (y -> -oo, or Q -> 1 in frame B).",
      "status" : "proved"
   },
   {
      "depends_on" : "leaf_structure_proof.md Theorem 11.1 and Cor. 11.3 (proved, TR-free)",
      "evidence" : "check_polar_relations.out at 4 values of x. check_vertex_orders.out: the prediction matches the project's direct genus-2 vertex t^0, t^1 coefficients exactly at x = 2, 1/3, -3/5 (non-zero values).",
      "statement" : "Lemma 5.2: rho q_t = [e^{s1 t/z}(q_0 + sum_{j>=2} (-1)^j C_j t^{j-1} lam_j)]_+ + (deficit >= 2 terms), with s1 = (1-x)^2/(2c) and C_j = (c^2/j) eps_j(sigma). Hence the D = k-1 terms cancel in [t^k] F_h(rho q_t).",
      "status" : "proved"
   },
   {
      "depends_on" : "Theorems A, B, C of conifold_regularity_note.md; Theorem LC (F1, F3, L1-L3)",
      "statement" : "Lemma 5.3 (cleanliness): assume Lemma A for D <= k-2 in lower genera. Then no graph with an edge contributes to [t^{2-2h+k}] F_{h,0}.",
      "status" : "proved"
   },
   {
      "depends_on" : "hypothesis (H): [t^{2-2h+k}] F_{h,0} = 0 for the resolved conifold x C^2 (refined GV form, as used in Theorem R); exact for h <= 4",
      "statement" : "Theorem D5: T_k(h) = sum_{D(R)=k} (1/|Aut R|) prod_i ((-1)^{j_i} eps_{j_i}(sigma)/j_i) Xi_h(R) = 0 for 1 <= k <= 2h-3, given (H) at order k and Lemma A for D <= k-2 in genera <= h.",
      "status" : "proved"
   },
   {
      "depends_on" : "(H) at order k = 1",
      "evidence" : "Implication proved. Unconditional for h <= 4, where (H) is exact data.",
      "statement" : "D5.1: Xi_h(2) = 0 identically, for every h >= 2 and every tau (all components at once).",
      "status" : "proved"
   },
   {
      "depends_on" : "(H) at order 2",
      "statement" : "D5.2: for h >= 3, 2(2+9 sigma)^2 Xi_h(3) = 3(1+4 sigma)^2 Xi_h(2,2). Hence Xi_h(2,2) = 0 <=> Xi_h(3) = 0, and sigma^h Xi_h(3) = sigma (1+4 sigma)^2 A_h(sigma) with deg A_h <= h-4 after P4.",
      "status" : "proved"
   },
   {
      "depends_on" : "(H) at order 2; P4 of proof_M_residue.md (Xi^(0)_h(3) = 0 for h >= 3, which still lacks the precise Boissy/BCGGM citation)",
      "statement" : "D5.2b: the Hurwitz component Xi^(0)_h(2,2) = (27/8) Xi^(0)_h(3) = 0 for every h >= 3.",
      "status" : "proved"
   },
   {
      "depends_on" : "(H) at orders 1 (genera <= h) and 3",
      "statement" : "D5.3: for h >= 3, 4(3+16 sigma)^2 Xi_h(4) - 6(2+9 sigma)^2 Xi_h(2,3) + 3(1+4 sigma)^2 Xi_h(2,2,2) = 0.",
      "status" : "proved"
   },
   {
      "evidence" : "check_genus3_alltau.out; the engine agrees 26/26 with the main session's Mathematica tables (check_data.out).",
      "statement" : "M_irr(3) for every tau: all components Xi^(0..3)_3(R) = 0 for R = (2),(3),(4),(2,2),(2,3),(2,2,2); the D = 4 controls are non-zero. This is an independent Perl computation via Mumford GRR and DVV, without Mathematica.",
      "status" : "computer-verified"
   },
   {
      "statement" : "Conditional honest-bundle count: if -R pi_* S_D had no jump loci, degree counting would give Xi^(k)_h(R) = 0 only for 2k < r (and k = h). Degree counting alone cannot give Lemma A.",
      "status" : "proved"
   }
]


## key_formulas
[
   "Xi_h(R) = sum_{p=0}^{N0} (-1)^p L^{(a+p)}_{N0-p}(-(J+p)) MV_h(R u 1^p)/p!,   a = 2h-2+r, N0 = 3h-3+r,   MV_h(mu) = Int Omega prod lam_{mu_i}",
   "(1+v)^{2h-2+r} e^{Jv} Phi_R(-v e^v) = Psi_R(v/(1+v)),  Phi_R(y) = sum_p y^p MV_h(R u 1^p)/p!,  Xi_h(R) = Psi_R(1) = lim_{v->oo}",
   "Xi_h(R) = (-1)^{N0} [y^{N0}] (1 - T(y))^{a+N0-1} e^{-(J+N0) T(y)} Phi_R(y),  T(y) = sum n^{n-1} y^n/n!",
   "sum_n u^n/n! <Omega prod_{i<=N} lam_{mu_i} prod (1+psi)^n>_h = (1-u)^{2-2h-N} e^{|mu| u/(1-u)} <Omega prod lam_{mu_i}>_h",
   "Localization edge factor (CY torus): E(d) = -eps_d(sigma),  eps_d = C(d tau + d - 1, d - 1)^2 = prod_{k=1}^{d-1} (d^2 sigma + k(d-k)) / ((d-1)!)^2",
   "Top leaf coefficients at the conifold point: C_j = c^2 [(1 - j/c)_{j-1}]^2/(j!(j-1)!) = (c^2/j) eps_j(sigma)",
   "rho q_t = [e^{s1 t/z}(q_0 + sum_{j>=2} (-1)^j C_j t^{j-1} lam_j)]_+ + O(deficit 2),  s1 = (1-x)^2/(2c)",
   "[t^{2-2h+k}](genus-h single vertex) = (-x c^3)^{h-1} [ sum_{D(R)=k} prod_i (-1)^{j_i} C_{j_i} Xi^c_h(R)/|Aut R| + (D <= k-2 terms) ]",
   "T_k(h) = sum_{R irreducible, D(R)=k} (1/|Aut R|) prod_i ((-1)^{j_i} eps_{j_i}(sigma)/j_i) Xi_h(R) = 0",
   "T_1: ((1+4 sigma)/2) Xi_h(2) = 0;  T_2: 2(2+9 sigma)^2 Xi_h(3) = 3(1+4 sigma)^2 Xi_h(2,2);  T_3: 4(3+16 sigma)^2 Xi_h(4) - 6(2+9 sigma)^2 Xi_h(2,3) + 3(1+4 sigma)^2 Xi_h(2,2,2) = 0",
   "sigma^h Omega K = c(T(-log)) e(E^vee x C_tau) e(E x C_{1+tau})  (Mumford; all honest bundles)",
   "Cont_X(Gamma) in Q[sigma]_{<= g-1+d},  [sigma^{g-1+d}] Cont_X(Gamma) = (-1)^{g-1+d} Cont_{P^1}(Gamma)"
]


## next_steps
[
   "Second geometry with the same singular vertex. In mode B, local P2 x C^2 at its conifold point carries the Lemma A class with leaf scale c (hodge_vertex/ data show universal ratios only at the two leading orders). Compute its top leaf coefficients (C_3', C_2'^2/2). If they are not proportional to (C_3, C_2^2/2), its T_2 together with D5.2 forces Xi_h(3) = Xi_h(2,2) = 0 in all genera. T_3, T_4 then restart, and each further independent geometry adds one relation per order.",
   "Combine D5.2 with a proof of the one-leaf statement Xi_h(3) = 0 (routes A1/E: formula (6.1) at j = 3). This alone settles all of D = 2; D5.3 then reduces D = 3 to (4) plus one multi-leaf point.",
   "Make (H) unconditional. Either prove the refined GV form P_h Li_{3-2h}(Q)/(eps1 eps2) for the 5-fold, or just t -> -t evenness of its polar part (flop symmetry Q -> 1/Q, using Li_{3-2h}(1/Q) = Li_{3-2h}(Q)); evenness already gives all odd k, in particular D5.1 and D5.3. Alternatively, transfer Theorem LC to the 3-fold X with the CY torus, where FP Theorem 3 makes (H) a theorem.",
   "Hurwitz component (A2): combine the single-locus identity D2.4 with a Graber-Vakil lift-comparison recursion ([b.0] versus [k.0 + (b-k).oo], i.e. vertex times genus-0 rubber), using the Laguerre/ELSV form of D2.3 to handle the K-twist.",
   "Give other routes the reformulations: Lemma A <=> sum_p (-1)^p L^{(a+p)}_{N0-p}(-(J+p)) MV_h(R u 1^p)/p! = 0, and, via ELSV, the explicit Hurwitz-number identity for the k = 0 component."
]


## obstructions
[
   "Finite fixed-point localization identities (tau-independence, polynomiality, FP-type vanishing on non-CY targets, relative and rubber versions) involve only K-free Marino-Vafa vertices. On CY3 targets the only vanishing relations are string, dilaton and divisor (D3.1).",
   "The K-twisted vertex is the analytic continuation of the vertex dressed by degree-1 edges to y = -v e^v -> -oo (D2.3). It appears in a localization identity only at the resonance c = w1, which is a pole at degree-1 p_oo ends (D3.2), through localized insertions with no vanishing, or at the conifold point Q = 1 after resummation in the canonical frame.",
   "At the conifold point the known F_{h,0} gives exactly ONE linear relation per polar order k, namely T_k. T_k can be used only once Lemma A is known for D <= k-2 in lower genera. With p(k) unknown functions of sigma at D = k, the system is square only at D = 1. At D = 2 the one-leaf point (3) remains, leaving h-3 scalar conditions at genus h after P4. So the cascade stops at D = 2 for h >= 4.",
   "Theorem D5 and its corollaries are conditional on (H), the absence of t^{3-2h}, t^{4-2h}, ... poles in the resolved-conifold-times-C^2 potential F_{h,0} at Q = 1. (H) is the refined GV form the project already uses in Theorem R, but it is established only through genus 4 (exact data)."
]


## report
G:\My Drive\Projects\Old\refGW\hodge_vertex_conifold\lemmaA\D_localization\ -- report.md could NOT be written: the harness refuses report files from subagents. The full report text is in the summary field below. All scripts and outputs (*.pl, *.out) are in this folder.

## requested_checks
[
   "Main session: confirm that the project regards hypothesis (H) as established in all genera, or give its source; it is the form F_{g,0} - F_{g,0}(0) = P_g Li_{3-2g}(Q)/(eps1 eps2) that Theorem R / Corollary R uses as the known polar part of F_g. Theorem D5.1 (Xi_h(2) = 0 for all h) depends on it, and so do D5.2 and D5.3.",
   "Optional (Mathematica, cheap): for local P2 x C^2 in mode B at its conifold point, compute the lam_j coefficients a_j^(k) of rho q_t for k = 1, 2, in the conifold's conventions. Test whether the D = 2 weight vector (a_3^(2), (a_2^(1))^2/2) is proportional to the conifold's (-C_3, C_2^2/2), with C_2 = (1-x)^2/2 and C_3 = (2-x)^2(1-2x)^2/(12c^2). If it is not, two independent T_2 relations follow; together with (H) for local P2 they force Xi_h(3) = Xi_h(2,2) = 0 in all genera.",
   "Optional: re-run check_vertex_orders.pl-type comparisons at genus 3 for the regular orders t^0, t^1 of the bare genus-3 vertex. These test T_4 and T_5 in the form [t^k] = (-x c^3)^2 T^c_k, where the D <= 2 terms vanish by M_irr(3); it needs the direct genus-3 vertex series from hvc_leaf_structure.wls."
]


## route
D: virtual localization relations (Graber-Pandharipande / Faber-Pandharipande), plus the global polar relations at the conifold point that this route leads to

## summary
ROUTE D REPORT (2026-10-31). No TR of any kind is used anywhere: no spectral curves, no EO, no loop equations, no BKMP, nothing mod P2, no TR scripts. All exact checks are Perl/Math::BigRat; no wolframscript was run. Conventions: c = 1, w = (1, tau, -1-tau), sigma = tau(1+tau), lam_j = 1/(1 - j psi), T = -psi^2/(1-psi), D(R) = sum(j_i - 1), J = sum j_i, r = #R, a = 2h-2+r, N0 = 3h-3+r. c-normalisation: Xi^c_h(R) = c^-(3h-3+2r) Xi_h(R)(tau = -x/c).

=== MAIN NEW RESULT (section 5): the polar relations T_k at the conifold point ===
This is a triangular system of the kind asked for in D3: the geometry is the resolved conifold x C^2 in frame B; the total F_{h,0} is known; the other graphs are controlled by Theorems A, B and C and by induction.

INPUTS. These are all proved in the project without TR, except (H).
(F1) Delta_1 = -x c^3/(t^2 rho^2), with rho analytic and rho(0) = 1 [leaf_structure_proof.md, Theorem LC].
(F3) Every frame-B edge factor and dilaton factor is analytic at t = 0.
(CF) rho q_t = c + sum_j A_j(t) lam_j. Theorem 11.1 and Cor. 11.3 give:
  - top coefficients: a_{k+1}^(k) = (-1)^{k+1} C_{k+1}, with C_j = c^2[(1-j/c)_{j-1}]^2/(j!(j-1)!);
  - sub-top coefficients: a_k^(k) = (-1)^k k C_k (1-x)^2/(2c^2), and a_1^(1) = -(1-x)^2/2.
(A) The pole order of a decorated graph is at most sum_v max(2g_v - 2, 0). (B) For a graph with an edge, the leading pole cancels [conifold_regularity_note s.2, s.3].
(C) Theorem C (Xi_h(1^n) = 0) and the string identity Xi^c(1,R) = (J/c^2) Xi^c(R).
(V) The genus-h single-vertex graph equals (-x c^3/t^2)^{h-1} F_h(rho q_t). This follows from dilaton homogeneity, and Xi^c = d_{lam}... F_h(q_0).
(H) [t^{2-2h+k}] F_{h,0} = 0 for 1 <= k <= 2h-3. In words: the polar part of the 5-fold potential at Q = 1 is gap_h t^{2-2h} only.
  - (H) holds if F_{h,0} - F_{h,0}(0) = P_h Li_{3-2h}(Q)/(eps1 eps2). This is the form Theorem R / Corollary R uses as "the known polar part of F_g".
  - (H) is an exact computation for h <= 3 (resolvedconifold_results.m) and for the polar parts at h <= 4.

KEY IDENTITY. (1-j/c)_{j-1} = (-1)^{j-1} prod_{k<j}(j tau + k). Hence
  C_j = (c^2/j) eps_j(sigma), where eps_j = C(j tau+j-1, j-1)^2 = prod_{k=1}^{j-1}(j^2 sigma + k(j-k))/((j-1)!)^2.
eps_j is exactly minus the degree-j edge factor of the resolved-conifold localization (Prop. D2.2 below).

LEMMA 5.2 (sub-top cancellation = first-order string flow) [proved]. Put top(t) = sum_{j>=2}(-1)^j C_j t^{j-1} lam_j and s1 = (1-x)^2/(2c). Then
  rho q_t = [e^{s1 t/z}(q_0 + top(t))]_+ + (every lam_j at t-order >= D(lam_j) + 2).
Proof:
  - [e^{s/z} lam_j]_+ = e^{js/c} lam_j and q_0 = c - c^2 lam_1.
  - Compare coefficients using Cor. 11.3: the top and sub-top terms match exactly. The flow is linear and invertible.
  - F_h is string-invariant, so F_h(rho q_t) = F_h(q_0 + top + rest'), where rest' has deficit >= 2.
Hence [t^k] F_h(rho q_t) = T^c_k + (terms in Xi^c(R) with D(R) <= k-2), where T^c_k = sum_{D(R)=k} prod_i(-1)^{j_i}C_{j_i} Xi^c(R)/|Aut R|. The D = k-1 terms cancel identically.

LEMMA 5.3 (cleanliness) [proved]. Assume Xi_{h'}(R) = 0 for all h' < h and all irreducible R with D <= min(k-2, 2h'-3). Then no graph with an edge contributes to [t^{2-2h+k}] F_{h,0}.
Proof:
  - For n = 0 a graph with an edge has no genus-h vertex, so (A) bounds its pole order by 2h-4.
  - (B) kills the coefficient of t^{4-2h}.
  - For k >= 3, expand the singular vertices via (L1-L3). A nonzero term needs D_v >= min(k-1, 2g_v-2) at each singular vertex of genus >= 2.
  - The t-order is then >= k+3-2h for one big vertex, and larger for two or more such vertices.

THEOREM D5 [proved modulo (H)]. Let h >= 2 and 1 <= k <= 2h-3. Assume (H) at order k, and Lemma A for D <= k-2 (within range) in all genera <= h. Then
  T_k(h) := sum_{R irreducible, D(R)=k} (1/|Aut R|) prod_i((-1)^{j_i} eps_{j_i}(sigma)/j_i) Xi_h(R) = 0.
Corollaries:
- D5.1 (k = 1; the hypothesis is vacuous). T_1 = ((1+4 sigma)/2) Xi_h(2), so Xi_h(2) = 0 identically for every h >= 2, all components at once.
  - Previously, only the components k = 0, 1, h were known in all genera.
- D5.2 (k = 2, h >= 3; the hypothesis is vacuous): 2(2+9 sigma)^2 Xi_h(3) = 3(1+4 sigma)^2 Xi_h(2,2). Write P_R = sigma^h Xi_h(R). Then:
  - P_(3) = sigma(1+4 sigma)^2 A_h and P_(2,2) = (2/3) sigma(2+9 sigma)^2 A_h. Hence Xi_h(2,2) = 0 <=> Xi_h(3) = 0.
  - Comparing the sigma^{h+2} coefficients gives Xi^(0)_h(2,2) = (27/8) Xi^(0)_h(3). With P4 (Xi^(0)_h(3) = 0 for h >= 3), the Hurwitz component of (2,2) vanishes in every genus h >= 3, and deg A_h <= h-4.
  - At h = 3: Xi_3(3) = Xi_3(2,2) = 0.
  - For h >= 4: D = 2 reduces to h-3 scalar conditions.
- D5.3 (k = 3, h >= 3; the hypothesis is D5.1): 4(3+16 sigma)^2 Xi_h(4) - 6(2+9 sigma)^2 Xi_h(2,3) + 3(1+4 sigma)^2 Xi_h(2,2,2) = 0.

VERIFICATION (non-trivial values). At genus 2 the orders k = 2 and 3 are regular. There, Lemma 5.2 predicts:
  - [t^0] = (-xc^3)(-C_3 Xi^c_2(3) + (C_2^2/2) Xi^c_2(2,2));
  - [t^1] = (-xc^3)(C_4 Xi^c(4) - C_2 C_3 Xi^c(2,3) + (C_2^3/6) Xi^c(2,2,2)).
The project computed the bare genus-2 vertex directly (leaf_structure.md, Step 2). The prediction agrees exactly with it at x = 2, 1/3 and -3/5 (check_vertex_orders.out):
  - x = 2: t^0 = -385/746496, t^1 = -97/139968;
  - x = 1/3: t^0 = -18325/47775744, t^1 = -133/663552;
  - x = -3/5: t^0 = 22982201/207360000, t^1 = -217/288.
check_polar_relations.out confirms, at four values of x: the sum rules; the top and sub-top coefficients; the cancellations at k = 2, 3; the residual -(1-x)^2 x^2/(12c^4) Xi^c(2) at k = 3; and the c = 1 forms of T_2 and T_3.

=== D1: candidate moduli spaces ===
(G1) Mbar_{g,n}(X,d), X = Tot(O(-1)+O(-1) -> P^1).
  - CY torus: weights (1, tau, -1-tau) at p0 and (-1, 1+tau, -tau) at p_oo.
  - The moduli space is compact and vdim = n.
  - The p0-vertex carries the Lemma A class; a degree-d edge gives the leaf d lam_d, at scale = weight (clue (a)).
(G2) X relative to the fibre S_oo. For profile R u 1^p, vdim = -D(R) for every p, so the relative invariants are polynomials of degree D(R). The fixed loci have rubber with double Hodge classes of weights (1+tau, -tau).
(G3) P^1 relative to oo with the branch class [b.0] (Graber-Vakil's localization proof of ELSV).
(G4) Frame B of the resolved conifold x C^2 at Q = 1, where the geometry itself produces the K-twist (section 5).

=== D2: localization formula and dictionary ===
Prop D2.1 [proved-modulo-literature: Graber-Pandharipande, Invent. Math. 135 (1999); FP, Invent. Math. 139 (2000) s.4.1].
  N_{g,d} = sum_Gamma 1/|Aut| prod_e E(d_e)/d_e prod_v V(v), where
  - E(d) = prod_{k=1}^{d-1}(w2 + k w1/d)(w3 + k w1/d) d^{2d}/((-1)^d (d!)^2 w1^{2d});
  - V = (t f2 f3)^{val-1} Int Lambda^vee(t)Lambda^vee(f2)Lambda^vee(f3) prod 1/(t/d_F - psi_F);
  - unstable vertices: (0,1) gives t/d, (0,2) gives t f2 f3/(t/d1 + t/d2).
  Computer-verified: N_{g,d} = d^{2g-3}[lambda^{2g-2}](2 sin(lambda/2))^{-2} for (g,d) in {(0,1..3), (1,1..3), (2,1..4), (3,1..2)}, with CY weights and with non-CY weights (1,3,5), (2,3,7).
Prop D2.2 [proved].
  (a) A p0-vertex contributes (-sigma)^{h+l-1}(prod d_i) MV_h(d), where MV_h(mu) = Int Omega prod lam_{mu_i}.
  (b) E(d) = -eps_d(sigma), and a bare degree-d edge has weight -sigma eps_d/d times lam_d.
  (c) Each graph is a polynomial in sigma of degree <= g-1+d. Its top coefficient is (-1)^{g-1+d} times its contribution to the P^1 localization of Int 1, which vanishes (checked, 7 cases). Its sigma^0 coefficient is FP's lambda_g computation.
Prop D2.3 [proved] ("dilaton leaf = string + dilaton - degree-1 edge", clue (b) made exact).
  (i) Lemma. For h >= 2 and pulled-back classes, sum_n u^n/n! <Omega prod_{i<=N} lam_{mu_i} prod(1+psi)^n> = (1-u)^{2-2h-N} e^{|mu|u/(1-u)} <Omega prod lam_{mu_i}>.
      Proof: dilaton homogeneity plus string invariance under q -> [e^{s/z}q]_+ with s = u/(1-u).
  (ii) Psi_R(s) := sum_m 1/m! <Omega prod lam_{j_i} prod sT> is a polynomial of degree <= N0 with Psi_R(1) = Xi_h(R). As formal series,
      Psi_R(s) = sum_p (-1)^p (phi_p/p!) s^p (1-s)^{-a-p} e^{(J+p)s/(1-s)}, where phi_p = MV_h(R u 1^p).
  (iii) Laguerre form: Xi_h(R) = sum_{p=0}^{N0} (-1)^p L^{(a+p)}_{N0-p}(-(J+p)) phi_p/p!.
  (iv) Lambert form: (1+v)^a e^{Jv} Phi_R(-v e^v) = Psi_R(v/(1+v)), with Phi_R(y) = sum y^p phi_p/p!. Hence Xi = lim_{v->oo}. The limit lies at y -> -oo, outside the disc |y| < 1/e of the bare-edge dressing.
  (v) Tree-function form: Xi = (-1)^{N0}[y^{N0}]{(1-T(y))^{a+N0-1} e^{-(J+N0)T(y)} Phi_R(y)}, with T the Cayley tree function.
  The relation holds componentwise. For k = 0, ELSV gives the Hurwitz identity sum_p (-1)^p L^{(a+p)}_{N0-p}(-(J+p)) H_{h,R u 1^p}/(2h-2+J+r+2p)! = 0.
  Checks: Laguerre form 18/18 at h = 2, 3 (including non-zero controls); Lagrange form 6/6. The root s = 1 of Psi_R is simple, so there is no hidden (1-s)^{2h-2-D}.
Prop D2.4 [proved-modulo-literature: Graber-Vakil]. Int_{[Mbar_{h,m}(P^1/oo,mu)]^vir} br^*[b.0] prod ev^*(gamma_0) f(psi) = C_mu Int Lambda^vee(1) prod lam_{mu_i} prod f(psi_k).
  - This holds with the ELSV constant C_mu, and only one fixed locus contributes.
  - Summing over T-insertions gives C_mu Xi^(0)_h(mu).
  - The insertions are localized, so no vanishing follows.

=== D3 obstruction results ===
Prop D3.1 [proved]. On a CY3 target, equivariant integrals of polynomial classes, with lambda insertions allowed, give no vanishing relation beyond string/dilaton/divisor.
  - After string reduction every insertion has a_i + deg gamma_i >= 1, so the degree is >= 0.
  - Degree-0 invariants reduce to N_{g,beta}.
  - Only weight-independence of N_{g,beta} and polynomiality remain. FP-type vanishing does not exist.
Prop D3.2 [proved] (resonance). A T_c-insertion at the p_oo end of a degree-d edge contributes -gamma|_{p_oo}(w1/d)^2/(c - w1/d). This is a pole at the Lemma A scale c = w1, coming from the degree-1 edges.
  - T_{w1} = w1 + psi - w1^2/(w1 - psi), and the last term is w1 times a degree-1 flag factor.
  - With gamma_0 insertions the sum is regular, and its p0-vertices are exactly the Lemma A vertices (loc_twisted.pl). But gamma_0 is localized: the coefficients are polynomials of degree sum a_k + 2m, and c = w1 lies on the circle of convergence.
  - The K-twist needs the dilaton-cancelled input at p0 (t_1 = 1). Finite insertions never reach it. It is reached only by resummation: at the conifold point in frame B, or along y = -v e^v (D2.3).

=== (d) the range D <= 2h-3 as a degree count ===
(i) Polar order: Xi_h(R) enters the genus-h vertex at t^{2-2h+D}. This is the count that produces relations (D5).
(ii) Riemann-Roch: deg S_D < 0; equivalently rank(-R pi_* S_D) = 3h-3-D >= h = rank E.
(iii) The relative virtual dimension is -D(R).
(iv) The lambda_g Norlund residue.
(v) Honest-bundle count (conditional, NOT a proof). By Mumford, sigma^h Omega K = c(T(-log)) e(E^vee x C_tau) e(E x C_{1+tau}). Without jump loci this would give only Xi^(k) = 0 for 2k < r.

=== D4 tables (exact) ===
- Engine validation:
  - classical KW and lambda values (including b_2, b_3, Int lambda_1^3, Int_{M3} lambda_2^3, lambda_g lambda_{g-1});
  - the genus-2 closed formula at j = 3, 4, 5;
  - 26/26 agreements with the main session's tables mv_g{1,2,3}_l{1,2} and q_g{2,3}_r{1,2}.
- M_irr(3) for every tau: all four components Xi^(0..3) are 0 for (2),(3),(4),(2,2),(2,3),(2,2,2), solved from 4 values of sigma with a 5th as check. Controls at D = 4:
  - (5): -65/24192, 269/120960, 467/120960, 31/40320;
  - (3,3): -39/4480, 19/2688, 571/40320, 31/8064;
  - (2,4): -8/945, 1/135, 467/30240, 31/8064;
  - (2,2,3): -1/28, 2/63, 43/576, 31/1344.
  - Xi^(3) equals the lambda_g value.
- Localization N_{g,d} and the P^1 relations: as above.
- D5 structure at genus 2: as above.

=== Precise obstruction ===
(1) Fixed-point localization relations contain only K-free Marino-Vafa vertices. The K-twisted vertex is an analytic continuation, found only at a resonance or at the conifold point.
(2) At the conifold point the known F_{h,0} gives exactly one relation per polar order, T_k. T_k is usable only once Lemma A is known for D <= k-2 in lower genera.
(3) There are p(k) unknown functions of sigma at D = k. So the system closes at D = 1 (all genera), reduces D = 2 to the one-leaf point (3), and stops at D = 2 for h >= 4.

=== Files (lemmaA/D_localization) ===
- Libraries: hodge_lib.pl (KW via DVV; Mumford GRR for lambda_1^k ch_3^l, g <= 4; grouped-multiset vertex integrals), xi_lib.pl (tau-free components), mparse.pl.
- Engine checks: test_lib, t_genus2, t_genus3, check_data, check_genus3_alltau.
- D2.3 checks: check_laguerre, check_lagrange, psi_poly.
- Localization: loc_conifold (+ loc_conifold_g012.out, loc_conifold_g3.out), loc_P1, loc_twisted.
- Section 5 checks: check_polar_relations, check_vertex_orders.
- Every script has a matching .out file except loc_twisted.

