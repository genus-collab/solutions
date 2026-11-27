# Moved on 2026-11-27 from localp2/hodge_vertex/lp2_explore/round6_scripts/conifold_apply/apply_lsp.py (local P2 project) into the resolved-conifold
# project, because it edits only files of this project. One-off patch script, already applied: a record, not to be rerun.
# Targeted in-place edits of leaf_structure_proof.md; each old string must occur exactly once.
import os
P=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "hodge_vertex_conifold", "leaf_proof", "leaf_structure_proof.md")
s=open(P,encoding='ascii',newline='').read()
def rep(old,new):
    global s
    n=s.count(old)
    assert n==1,(n,old[:80])
    s=s.replace(old,new)

# header
rep("""confirmations only. Revised 2026-10-31 after the completeness critic and the Mathematica checks of s.12.2(e) (all passed).
""","""confirmations only. Revised 2026-10-31 after the completeness critic and the Mathematica checks of s.12.2(e) (all passed).
Revised 2026-11-02 (round 6: s.0.4, s.12.1, s.13, s.14 only; Theorem LC and sections 1-11 are unchanged).
""")
# contents
rep("""F1-F3. 11 Closed forms. 12 What is proved and what is computer-verified. 13 Consequence: the all-genus regularity theorem is conditional
only on M_irr. 14 Remaining gaps, errata, open points. Appendix A notation; Appendix B concordance with the section files.
""","""F1-F3. 11 Closed forms. 12 What is proved and what is computer-verified. 13 Consequence: the all-genus regularity theorem (M_irr = Lemma A,
proved modulo (MV)). 14 Remaining gaps, errata, open points. Appendix A notation; Appendix B concordance with the section files.
""")
# s.0.4
rep("""**Corollary R (section 13).** Fix g >= 2 and assume M_irr(h) for 2 <= h <= g (proof_status_mixed_scale.md s.1.2). Then every coloured
stable graph of genus g with at least one edge is regular at t = 0 after summing over decorations, and the single-vertex graph has the
single pole gap_g t^(2-2g). M_irr(2) is proved, and M_irr(3) and M_irr(4) are exact computations (proof_status_mixed_scale.md, addendum
(1) and s.1.9). So the theorem holds through genus 4 (computer-assisted, independently of the direct graph sums), and in every genus it
is conditional only on M_irr(h), h >= 5, i.e. on Lemma A of proof_status_mixed_scale.md s.2. Lemmas B (P2) and D, the
topological-recursion route to Lemma A, are not needed.
""","""**Corollary R (section 13).** Fix g >= 2. Every coloured stable graph of genus g with at least one edge is regular at t = 0 after summing over
decorations (this uses M_irr(h) only for 2 <= h <= g - 1, Lemma 13.0), and the single-vertex graph has the single pole gap_g t^(2-2g) (this
also uses M_irr(g), and the value comes from Theorem VI of the local-P2 draft, s.13.3 (5)). M_irr(h), for every h >= 2 and non-empty R, is
Theorem L of lemmaA/LEMMA_A_PROOF.md (s.7.4 there for the conifold convention), P modulo the Marino-Vafa formula (MV). So in every genus
Corollary R is unconditional modulo (MV), textbook facts and, only to identify the graph sum with F_g, the inputs of Proposition 13.2 (exact
list and refereeing status: s.13.2). Genus 3 and 4 are no longer computer-assisted: the exact computations there are now independent
confirmations. Lemmas B (P2) and D, the topological-recursion route to Lemma A, are not needed.
""")
# s.12.1
rep("""series, as Proposition 3.15].
Not proved here and not needed: see section 14.
""","""series, as Proposition 3.15].
 - s.13.0: the domain of Fpot_h, (H0), Lemma 13.0 [f + dilaton equation; port of the local-P2 draft, 2026-11-02]. Theorem R(a), (b) [given F1, F3, Lemma C, M_irr; the value of the
   pole by Theorem VI, s.13.3 (5)]. Corollary R [P modulo (MV) and the literature inputs listed in s.13.2].
Not proved here and not needed: see section 14.
""")
# s.13 heading and preamble
rep("""## 13. Consequence: the all-genus regularity theorem is conditional only on M_irr

Notation of proof_status_mixed_scale.md: Fpot_h denotes the genus-h potential of the Hodge CohFT Omega_w in dilaton-shifted coordinates (called F_h
there), Xi_h(j_1, ..., j_r) := d_{lam_{j_1}} ... d_{lam_{j_r}} Fpot_h(q_0) (a finite sum of Hodge integrals at each order, since the lam_j are power
series in psi), J_R := Sum j_i, D(R) := Sum (j_i - 1), and
    M_irr(h):  Xi_h(R) = 0 for every multiset R of integers >= 2 with D(R) <= 2h - 3.
In this section g denotes the genus.
""","""## 13. Consequence: the all-genus regularity theorem (M_irr = Lemma A, proved modulo (MV))

Notation of proof_status_mixed_scale.md: Fpot_h denotes the genus-h potential of the Hodge CohFT Omega_w (called F_h there; domain and definition in
s.13.0), Xi_h(j_1, ..., j_r) := d_{lam_{j_1}} ... d_{lam_{j_r}} Fpot_h(q_0) (a finite sum of Hodge integrals at each order, since the lam_j are power
series in psi), J_R := Sum j_i, D(R) := Sum (j_i - 1), and
    M_irr(h):  Xi_h(R) = 0 for every non-empty multiset R of integers >= 2 with D(R) <= 2h - 3.
The empty multiset is excluded: Xi_h(empty) = Z[h, 0] = Fpot_h(q_0) is the gap, non-zero in Q(x) (Theorem VI, s.13.3 (5); no proof step uses the non-vanishing).
In this section g denotes the genus.

### 13.0 The domain of Fpot_h, dilaton homogeneity (H0), and the range of genera
[Added 2026-11-02, round 6: port of LP2_REGULARITY_PROOF_draft.md (hodge_vertex/lp2_explore/) s.1.5, (H0) of s.7.1 and Lemma 7.1, with E -> c, t_c -> t, colour 0 ->
idempotent 1. Every line was re-derived for the port; nothing here depends on the rank. Not yet refereed in this form.]
Let h >= 2, Rr := C[[t]] (or Q(x)[[t]]), Omega_w as in Proposition 13.2 (a pull-back along forgetful maps, so the string and dilaton equations hold), and
T(psi) := -psi^2/(c - psi), so that q_0 = T - psi.
 - Base-point derivatives. For g_1, ..., g_s in Rr[[psi]],
       dd^s(g_1, ..., g_s) := Sum_m (1/m!) Int_{Mbar_{h,s+m}} Omega_w prod_i g_i(psi_i) prod_{k=1}^m T(psi_{s+k}),
   a finite sum, since T = O(psi^2) forces m <= 3h - 3 + s.
 - Fpot_h on the neighbourhood Dom of q_0. Dom := q_0 + t Rr[[psi]]. For q = q_0 + delta in Dom and f_1, ..., f_n in Rr[[psi]],
       d^n Fpot_h(q)[f_1, ..., f_n] := Sum_{r>=0} (1/r!) dd^(n+r)(f_1, ..., f_n, delta, ..., delta),     Fpot_h(q) := d^0 Fpot_h(q),
   convergent t-adically because delta = O(t). delta may have a psi-linear term: Dom contains rho q_t, whose psi-coefficient is -rho (s.11.1), so q need not be
   of the form -psi + O(psi^2). This is why the domain has to be stated.
 - The dilaton-shifted form. If q = -psi + O(psi^2), then Fpot_h(q) = Sum_m (1/m!) Int_{Mbar_{h,m}} Omega_w prod_{i=1}^m (q + psi)(psi_i), a finite sum, and similarly for
   d^n Fpot_h(q)[f_1, ..., f_n]. On the overlap with Dom the two definitions agree: expand prod (T + delta) binomially, and (1/m!) C(m, r) reindexes to (1/r!)(1/(m-r)!).
   This is the form of (C10): the decorations (dilaton degrees d >= 2, weight 1/m!) at a vertex of genus h >= 2, colour i, with half-edge series f_1, ..., f_n
   sum to Delta_i^(h-1) d^n Fpot_h(q_t^(i))[f_1, ..., f_n], multilinear in the slots (applied term by term to loops and to E_11 edges between rescaled vertices).
 - q_t lies in Dom at idempotent 1: Dil_1(d) is in Kx[[t]] with q_t|_{t=0} = q_0 (Theorem 5.4, remark (iv)). rho is in 1 + t Rr (F1), so rho q_t lies in Dom as well. At idempotent 2,
   q_t|_{t=0} = -psi h_w(-psi)^-2 is not q_0, so q_t is not in Dom; only the dilaton-shifted form is used there, and nothing is rescaled.
 - Xi_h(R) = dd^r(lam_{j_1}, ..., lam_{j_r}) for R = {j_1, ..., j_r}, and Z[h, n] = Xi_h(1^n) (F4), so Z[h, 0] = Fpot_h(q_0).
**(H0) Dilaton homogeneity on Dom.** For q in Dom and lambda in 1 + t Rr independent of psi, lambda q is in Dom and
       d^n Fpot_h(lambda q)[f_1, ..., f_n] = lambda^(2-2h-n) d^n Fpot_h(q)[f_1, ..., f_n],   equivalently
       d^n Fpot_h(q)[f_1, ..., f_n] = lambda^(2h-2) d^n Fpot_h(lambda q)[lambda f_1, ..., lambda f_n].                                     (13.0)
Proof. Step 1 (base point). Let A_m := (1/m!) Int Omega_w prod g_i prod_{k<=m} T. Since q_0 = T - psi, reindexing the T-sum and using the dilaton equation
Int Omega_w psi_0 prod g_i prod T = (2h - 2 + s + m) Int Omega_w prod g_i prod T give dd^(s+1)(g_1, ..., g_s, q_0) = Sum_m [m - (2h - 2 + s + m)] A_m = (2 - 2h - s) dd^s(g_1, ..., g_s).
Step 2 (Euler identity). For q = q_0 + delta, d^(n+1) Fpot_h(q)[f, q] = Sum_r (1/r!) [dd^(n+r+1)(f, delta^r, q_0) + dd^(n+r+1)(f, delta^(r+1))] = (2 - 2h - n) d^n Fpot_h(q)[f], because by Step 1
the terms with the factor r telescope against the last sum. Step 3 (integration). Put mu := lambda - 1 in t Rr, q_s := (1 + s mu) q with s a formal variable, and G(s) := d^n Fpot_h(q_s)[f].
Each t-coefficient of G is polynomial in s; the definition gives dG/ds = mu d^(n+1) Fpot_h(q_s)[f, q], and Step 2 at q_s gives (1 + s mu) dG/ds = mu (2 - 2h - n) G. Hence
G(s) = (1 + s mu)^(2-2h-n) G(0); put s = 1. QED
**Lemma 13.0 (range).** Let Gamma be a connected stable graph of genus g, without legs, with at least one edge. Then every g_v <= g - 1. So Theorem R(a) uses M_irr(h) only
for 2 <= h <= g - 1, and at g = 2 there is no vertex of genus >= 2 at all (Theorem A of the note). Proof. If h^1(Gamma) >= 1 this follows from g = h^1(Gamma) + Sum g_v. If Gamma is a
tree with >= 2 vertices it has at least two vertices of valence 1, each of genus >= 1 by stability, so no vertex carries all of g. QED
""")
# 13.1 statement
rep("""**Theorem R.** Fix g >= 2. Assume F1, F3, Lemma C through order t^(2h-3) for every 2 <= h <= g (at the psi-orders that occur in genus-g graphs), and
M_irr(h) for every 2 <= h <= g. Then every coloured stable graph of genus g with at least one edge is regular at t = 0 after summing over decorations,
and the single-vertex graph has the single pole (-x c^3)^(g-1) Fpot_g(q_0) t^(2-2g). Since the conifold F_g is P_g Li_{3-2g}(Q)/(eps1 eps2) up to
Q-independent terms, whose only singular term is (2g-3)! t^(2-2g), this pole is gap_g t^(2-2g), and F4 holds in genus g as a consequence.
""","""**Theorem R.** Fix g >= 2. Assume F1, F3, Lemma C through order t^(2h-3) for every 2 <= h <= g (at the psi-orders that occur in genus-g graphs).
 (a) If M_irr(h) holds for every 2 <= h <= g - 1, every coloured stable graph of genus g with at least one edge is regular at t = 0 after summing over decorations.
 (b) If M_irr(g) holds, the single-vertex graph has the single pole (-x c^3)^(g-1) Fpot_g(q_0) t^(2-2g). By Theorem VI (input (5), s.13.3) this is gap_g t^(2-2g), where
     gap_g := beta_{2g-2}(x, 1) = (2g-3)! P_g(x, 1)/x is the t^(2-2g) coefficient of P_g Li_{3-2g}(Q)/(eps1 eps2) (the only singular term of Li_{3-2g}(e^t) is
     (2g-3)! t^(2-2g)). So F4 at n = 0 holds in genus g.
""")
# 13.1 (i)
rep("""and all t-dependence sits in Delta, Dil and E. The dilaton equation gives Fpot_h(lambda q) = lambda^(2-2h) Fpot_h(q), so
d^n Fpot_h(q)[f_1, ..., f_n] = lambda^(2h-2) d^n Fpot_h(lambda q)[lambda f_1, ..., lambda f_n]; with lambda = rho and (4.2) the vertex factor is
""","""and all t-dependence sits in Delta, Dil and E. By (H0), (13.0), on the domain Dom of s.13.0 (q_t and rho q_t lie in Dom, although the psi-coefficient of
rho q_t is -rho), with lambda = rho and (4.2) the vertex factor is
""")
# 13.1 (ii)
rep("""D). Taylor-expanding Fpot_h(q_0 + Sum_k t^k q^(k)), the coefficient of t^(2-2h+k) of the vertex factor is a finite combination of Xi_h(R) with all parts
j >= 1 and D(R) <= k. (iii)""","""D). By the definition on Dom, Fpot_h(q_0 + Sum_k t^k q^(k)) is its Taylor series at q_0, so the coefficient of t^(2-2h+k) of the vertex factor is a finite combination of Xi_h(R) with all parts
j >= 1 and D(R) <= k. (iii)""")
# 13.1 (iii)
rep("""(Theorem C); if R' is non-empty, all its parts are >= 2 and D(R') = D(R) <= k, so Xi_h(R') = 0 by M_irr(h) whenever k <= 2h - 3. So a non-zero polar
""","""(Theorem C); if R' is non-empty, all its parts are >= 2 and D(R') = D(R) <= k, so Xi_h(R') = 0 by M_irr(h) whenever k <= 2h - 3 (h <= g - 1 in a graph with an edge, Lemma 13.0). So a non-zero polar
""")
# 13.1 (iv)-(v)
rep("""vertices) + 0 (genus-1 singular and all regular vertices, by F1 and F3) >= 0. (v) The single-vertex graph has no edges: its only term with R empty is
the order k = 0, (-x c^3)^(g-1) Fpot_g(q_0) t^(2-2g); all its other polar terms have R non-empty with D <= k <= 2g - 3 and vanish. Comparing the sum
over all graphs with the known polar part of F_g identifies this term with gap_g t^(2-2g). QED
""","""vertices) + 0 (genus-1 singular and all regular vertices, by F1 and F3) >= 0. This proves (a). (v) The single-vertex graph has no edges, and its colour-2 term is regular (F1, F3). At colour 1 its
only term with R empty is the order k = 0, (-x c^3)^(g-1) Fpot_g(q_0) t^(2-2g); all its other polar terms have R non-empty with D <= k <= 2g - 3 and vanish by (iii) with h = g. The value is
input (5): (-x c^3)^(g-1) Z[g, 0] = beta_{2g-2}(x, 1) = gap_g by Theorem VI. QED
Independent check of (5), not a step: by Proposition 13.2 the sum over all graphs is F_g, by (a) only the vertex graph is singular, and by [BS] Prop. 3.1 the polar part of F_g is gap_g t^(2-2g).
""")
# 13.2
rep("""**Corollary R.** Fix g >= 2 and assume M_irr(h) for every 2 <= h <= g. Then every coloured stable graph of genus g with at least one edge is regular at
t = 0 after summing over decorations, the single-vertex graph has the single pole gap_g t^(2-2g), and F4 holds in genus g. In particular:
 - genus 2: unconditional, and without any computer input (M_irr(2) is proved: proof_status_mixed_scale.md s.1.4);
 - genus 3 and 4: unconditional, computer-assisted (M_irr(3): the exact lattice data hvc_multiprobe_out.txt; M_irr(4): hvc_multiprobe_g4_out.txt for
   r <= 2 and hvc_multiprobe_g4b.wls for the seven points with r >= 3; proof_status_mixed_scale.md s.1.9 and addendum (1)); this is independent of
   the direct genus-4 graph sums, which already show regularity;
 - genus g >= 5: conditional only on M_irr(h) for 5 <= h <= g, i.e. on Lemma A of proof_status_mixed_scale.md s.2 (its case h = 4 is now settled).
Proof. Theorem LC gives F1, F3 and Lemma C (L1, L2 for j = 1, 2 and every b, L3) for every k and at every psi-order, so the hypotheses of Theorem R
other than M_irr hold in every genus. The orientation issue of a mixed edge (graphValueC attaches EdN[c(v1), c(v2)] in the order of the edge list,
""","""**Corollary R.** For every g >= 2 and every x with x c != 0: every coloured stable graph of genus g with at least one edge is regular at t = 0 after summing over decorations, the
single-vertex graph has the single pole gap_g t^(2-2g), and F4 at n = 0 holds in genus g.
Status (2026-11-02). P modulo: (MV) = LLZ eq. (1) (for M_irr, and for Theorem VI); the textbook facts of lemmaA/LEMMA_A_PROOF.md s.0; for the value of the pole only, the further inputs of Theorem VI (LLZ eq. (5), the
cut-and-join eigenvalue; the vertex-operator form of psi, psi^*). F1, F3 and Lemma C are Theorem LC. The inputs of Proposition 13.2 are needed only to identify the graph sum with F_g: regularity and the pole value are
statements about the sum as defined by (C6)-(C11). No TR is used. Refereeing: LEMMA_A_PROOF.md (frozen copy LEMMA_A_PROOF_round5.md) has had two round-6 referees (logic and literature; independent computation), both "sound with minor fixes", no major item, nothing that changes the mathematics of its s.2-s.6;
Theorem VI has had one referee pass (round 5); s.13.0, s.13.3 (5) and the field-of-definition paragraph of Proposition 13.2 have had none in this form. By genus:
 - genus 2: edge graphs use no M_irr (Lemma 13.0; all vertices have genus <= 1, Theorem A). The vertex graph uses M_irr(2), proved without (MV) and without a computer (proof_status_mixed_scale.md s.1.4);
 - genus 3: edge graphs use only M_irr(2). The vertex graph uses M_irr(3) = Theorem L;
 - genus g >= 4: edge graphs use M_irr(h) for 2 <= h <= g - 1, the vertex graph also M_irr(g), all Theorem L.
Genus 3 and 4 are no longer computer-assisted. The exact data (M_irr(3): hvc_multiprobe_out.txt; M_irr(4): hvc_multiprobe_g4_out.txt for r <= 2 and hvc_multiprobe_g4b.wls for the seven points with r >= 3;
proof_status_mixed_scale.md s.1.9 and addendum (1)) and the direct genus-4 graph sums are independent confirmations, not proof steps.
Proof. Theorem LC gives F1, F3 and Lemma C (L1, L2 for j = 1, 2 and every b, L3) for every k and at every psi-order, so the hypotheses of Theorem R
other than M_irr hold in every genus. M_irr(h) holds for every h >= 2, identically in x and at every x with x c != 0, x = 1 (tau = -1/2) included: lemmaA/LEMMA_A_PROOF.md Theorem L, Cor. L1, Cor. L3 and s.7.4 (E = c, tau = -x/c,
sigma = -x/c^2; parts 1 are stripped by the string identity, as in (iii)). The orientation issue of a mixed edge (graphValueC attaches EdN[c(v1), c(v2)] in the order of the edge list,
""")
# 13.3 inputs
rep("""Inputs (2)-(5) are proved in the cited notes, not here; input (1) is Proposition 13.2 below: (1) the frame-B Hodge-vertex graph-sum formula
for the conifold free energies (conifold_regularity_note.md,
preamble; hvc_lib.m), the basis of the whole project, with the vertex class Omega_w t-independent; (2) the dilaton equation for Fpot_h (homogeneity);
(3) Theorem C (conifold_regularity_note.md s.3a, by the string equation); (4) the string/dilaton reduction of Conjecture M to M_irr
(proof_status_mixed_scale.md s.1.2); (5) the known polar part of the conifold F_g (the definition of gap_g in F4).
""","""Input (1) is Proposition 13.2 below, input (2) is s.13.0, inputs (3), (4) are proved in the cited notes and input (5) in the local-P2 draft, not here: (1) the frame-B Hodge-vertex graph-sum formula
for the conifold free energies (conifold_regularity_note.md,
preamble; hvc_lib.m), the basis of the whole project, with the vertex class Omega_w t-independent; (2) the dilaton equation for Fpot_h (homogeneity), in the form (H0);
(3) Theorem C (conifold_regularity_note.md s.3a, by the string equation); (4) the string/dilaton reduction of Conjecture M to M_irr
(proof_status_mixed_scale.md s.1.2); (5) the value of the single pole, (-x c^3)^(g-1) Z[g, 0] = gap_g: Theorem VI of hodge_vertex/lp2_explore/LP2_REGULARITY_PROOF_draft.md (s.7.6, proof s.7.8), a statement about Hodge integrals
only, with no TR, no Teleman and nothing from Lemma A; [BS] Prop. 3.1 (BS v2: prop:ResConifold) is an independent check, no longer an input.

Specialisation of Theorem VI to the conifold (input (5)). Theorem VI_tau (draft s.7.8) reads sigma^(h-1) Z_h(tau) = beta_{2h-2}(a_1, a_2) in Q(tau), with a_1 = -tau, a_2 = 1 + tau, sigma = tau(1 + tau), Z_h the value at leaf scale 1 and
beta_n(a, b) := (n-1)! [t^n] 1/(4 sinh(a t/2) sinh(b t/2)). Frame B has w = (c, -x, -1) ((C8): the weights of h_w, with Sum_a w_a^(1-2k) = -St_k), and c - x - 1 = 0, so w = c (1, tau, -1-tau) with tau = -x/c, -1-tau = -1/c,
sigma = -x/c^2, a_1 = x/c, a_2 = 1/c; Q(tau) = Q(x), since tau = -x/(1+x) is a Moebius change of variable. Omega_w, T and lam_1 are those of F4, with leaf scale c = w_1. The weight count (lemmaA/LEMMA_A_PROOF.md Cor. L3 at r = 0; draft (7.2)) gives
Z[h, 0] = c^(-(3h-3)) Z_h(-x/c). So
    (-x c^3)^(h-1) Z[h, 0] = (-x)^(h-1) Z_h = sigma^(h-1) c^(2h-2) Z_h = c^(2h-2) beta_{2h-2}(x/c, 1/c) = beta_{2h-2}(x, 1),
using beta_n(l a, l b) = l^n beta_n(a, b). Z_h is a polynomial in 1/sigma (the expansion of lemmaA/LEMMA_A_PROOF.md s.7.1 at R empty), so this holds at every x with x c != 0, x = 1 included. On the other side, P_g := [t^(2g)] e1 e2 t^2/(4 sinh(e1 t/2) sinh(e2 t/2)) (rc_lib.m line 244), so
P_g/(eps1 eps2) = beta_{2g-2}(eps1, eps2)/(2g-3)!, and Li_{3-2g}(e^t) = (2g-3)! t^(2-2g) + (analytic), since Li_{-n}(e^t) = (d/dt)^n (e^t/(1 - e^t)) and the only singular term of e^t/(1 - e^t) is -1/t. Hence gap_g = beta_{2g-2}(x, 1), with no BS input. By hand at h = 2:
Z_2 = -(24 + 28 s + 7 s^2)/5760 with s = 1/sigma = -c^2/x gives -x Z_2 = (24 x^2 - 28 c^2 x + 7 c^4)/(5760 x) = (7 + 10 x^2 + 7 x^4)/(5760 x), the gap of README.txt. Non-vanishing (not used): with S(y) = (y/2)/sinh(y/2) = Sum c_i y^(2i), sign c_i = (-1)^i, one has beta_{2g-2}(a, b) = ((2g-3)!/(a b)) Sum_{i+j=g} c_i c_j a^(2i) b^(2j), all of whose terms have one sign at real a, b != 0.
Evidence, not a step: F4 [exact, genus 2-4, symbolic x] of the note; (H0) and this specialisation also re-checked 2026-11-02 in exact arithmetic from Mumford/DVV Hodge integrals at x = 2/3, 3, -5/7, 1 (h = 2, 3, 4; also h = 5 at x = 2/3).
""")
# Prop 13.2 (a) Teleman
rep("""Teleman's classification ([T12] s.1.9, Theorem (2); s.8.1, Proposition 8.1(ii); s.8.2; applied over the field Kw((Q))), Omega^T = R.omega for the
""","""Teleman's classification ([T12] (1.15), Theorem (2); Proposition 8.7(ii); applied over the field Kw((Q)), see the end of this step), Omega^T = R.omega for the
""")
rep("""QDE ([G01] s.1.3, (3)-(5)); two such R differ by right multiplication by Exp[diag Sum_k a_k z^(2k-1)] with constant a_k ([G01] s.1.3, Proposition
(c); its proof is formal""","""QDE ([T12] Propositions 7.7, 8.2(ii); [G01] s.1.3, (3)-(5)); two such R differ by right multiplication by Exp[diag Sum_k a_k z^(2k-1)] with constant a_k ([T12] Proposition 8.4; [G01] s.1.3, Proposition
(c); its proof is formal""")
rep("""Lemma 1.1 reads E_ij(z, w) = [delta_ij Delta_i - <v_i(-z), v_j(-w)>]/(z + w), the edge of [PPZ] Def. 2.2. So F^T_g := Int_{M_g} Omega^T_{g,0} is the
point-frame sum for the rows of R.
""","""Lemma 1.1 reads E_ij(z, w) = [delta_ij Delta_i - <v_i(-z), v_j(-w)>]/(z + w), the edge of [PPZ] Def. 2.2. So F^T_g := Int_{M_g} Omega^T_{g,0} is the
point-frame sum for the rows of R.
Field of definition. [T12] is stated over C (it allows the generic point of the Novikov variables, [T12] (1.6)). Kw((Q)) has characteristic 0 and cardinality 2^aleph_0, so it embeds in C as a field (send a transcendence basis over Q injectively into one of
C and extend, C being algebraically closed), and [T12] is applied to the base change of Omega^T along such an embedding. The embedding forgets d/dt, which comes back as follows [round 6, not refereed]: the datum R of a given CohFT is unique ([T12] (1.15), Proposition 8.7(ii)), so it is
fixed by Aut(C/Kw((Q))), has coefficients in Kw((Q)) and does not depend on the embedding; by the divisor equation (the twist is pulled back, [CG07] App. 1(iv)), Exp[s Q d/dQ] Omega^T is Teleman's family at u = s H ([T12] Definition 7.1, Example 7.17, the formal case included); and uniqueness
over C((s)), again of cardinality 2^aleph_0, gives R_{sH} = Exp[s Q d/dQ] R, so the ODE of [T12] Proposition 8.2(ii) along H is the QDE in d/dt = Q d/dQ.
""")
rep("""z^(2k-1)] acts on 1 by Exp[Sum_k b_k ((2k)!/B_2k) ch_{2k-1}(E)] ([T12] s.8.2; [G01] s.2.3, Lemma), so Hodge_v = h_v.1.""","""z^(2k-1)] acts on 1 by Exp[Sum_k b_k ((2k)!/B_2k) ch_{2k-1}(E)] ([T12] (8.6), Proposition 8.7(i); [G01] s.2.3, Lemma), so Hodge_v = h_v.1.""")
# refs
rep("""[T12] C. Teleman, The structure of 2D semi-simple field theories, Invent. Math. 188 (2012) 525-588 (0712.0160). [Z16]""","""[T12] C. Teleman, The structure of 2D semi-simple field theories, Invent. Math. 188 (2012) 525-588 (0712.0160; numbering as in arXiv v3, read 2026-11-02; journal numbering not compared). [Z16]""")
# G1
rep(""" (G1) The all-genus regularity theorem remains conditional on M_irr(h) for h >= 5 (Lemma A). Nothing in this text bears on it.
""",""" (G1) Closed modulo (MV), 2026-11-02: M_irr(h) for every h >= 2 and non-empty R is Theorem L of lemmaA/LEMMA_A_PROOF.md (s.7.4 for the conifold convention). Graphs with an edge need it only for h <= g - 1 (Lemma 13.0). Status of Corollary R: s.13.2.
""")
# G3
rep("""formula and the group property of the Givental action, which forces s = -1. Cited, not re-proved: Teleman (Invent. Math. 188, 2012: s.1.9
      Theorem (2), s.8.1-8.2); Givental""","""formula and the group property of the Givental action, which forces s = -1. Cited, not re-proved: Teleman (Invent. Math. 188, 2012; numbering of arXiv:0712.0160v3: (1.6), (1.15)
      Theorem (2), Definition 7.1, Example 7.17, Propositions 7.7, 8.2(ii), 8.4, (8.6), Proposition 8.7); Givental""")
rep("""      localisation theorem (compact targets) and Zong's preprint arXiv:1604.07270 are quoted only for comparison. Inputs (2)-(5) of s.13.3 are
      proved in the cited notes and are not re-examined here ((5): Brini-Schuler Prop. 3.1, proved by degeneration).
""","""      localisation theorem (compact targets) and Zong's preprint arXiv:1604.07270 are quoted only for comparison. Inputs (3), (4) of s.13.3 are
      proved in the cited notes and are not re-examined here; (2) is s.13.0; (5) is Theorem VI of the local-P2 draft, and Brini-Schuler Prop. 3.1 (proved by degeneration) is now an independent check, used by no proof step.
""")
# errata
rep("""Open points that are not gaps of the proof: the closed forms of section 11""",""" (E7) Round 6, 2026-11-02 (backup leaf_structure_proof_before_round6.md): s.13.0 added (domain of Fpot_h, (H0), Lemma 13.0); "non-empty" added to M_irr; Teleman citations corrected against arXiv v3 and the field of definition written out; input (5) and step (v) of Theorem R restated
      (Theorem VI, [BS] Prop. 3.1 only as a check); status updated (Lemma A = lemmaA/LEMMA_A_PROOF.md, s.13.2).

Open points that are not gaps of the proof: the closed forms of section 11""")
open(P,'w',encoding='ascii',newline='').write(s)
print("ok")
