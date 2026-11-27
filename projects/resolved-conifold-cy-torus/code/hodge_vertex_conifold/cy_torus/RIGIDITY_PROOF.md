# Corollaries of the generic-weight leaf theorem on the CY 4-torus: regularity (Corollary R_CY) and rigidity for the resolved conifold

File: hodge_vertex_conifold/cy_torus/RIGIDITY_PROOF.md. Written 2026-11-03. Companion of LC_CY_PROOF.md (same folder, written separately), whose main theorem is taken here as a hypothesis (s.2). Referee status: one hostile referee pass (AI), verdict "sound with minor fixes, no major items"; its items are applied (record in s.8.4). Nothing in this file has been read by a human.

Tags. **[P]** proved here. **[P mod X]** proved here modulo the named input X. **[C]** conditional on Hypothesis LC_CY of s.2. **[E]** evidence (finite-order or numerical), used by no proof step. **[cited]** a statement of another file of this project or of the literature, used as stated.

Rules respected. No topological recursion for the target theory, and none anywhere below. No finite-order computation fixes a formula or a constant: s.9 is evidence only. No wolframscript or Mathematica was run; the finite-order checks are exact rational arithmetic in Python (s.9).

Files referred to (read-only): LSP := hodge_vertex_conifold/leaf_proof/leaf_structure_proof.md; LAP := hodge_vertex_conifold/lemmaA/LEMMA_A_PROOF.md; LP2 := hodge_vertex/lp2_explore/LP2_REGULARITY_PROOF_draft.md; FB := hodge_vertex/lp2_explore/FRAMEB_PROOF.md; CRN := hodge_vertex_conifold/conifold_regularity_note.md; FEAS := hodge_vertex_conifold/cy_torus_explore/feasibility_result.json; [BS] := BS_v2_2026-11-02/main.tex, cited by label and line.

---

## 0. Statements and status

### 0.1 Setting (details in s.1)

Z = (resolved conifold) x C^2 with the 5-torus T acting with tangent weights wt_+ = (a0, a1, a2, -e1, -e2) at p_+ and wt_- = (-a0, a0+a1, a0+a2, -e1, -e2) at p_-, the two fixed points of the zero section P^1. The CY 4-torus is T_CY: E := a0 + a1 + a2 = e1 + e2. Put A := a0 + a2, B := a0 + a1, Kw := Q(a0, a1, a2, e1) (the function field of Lie(T_CY), with e2 := a0 + a1 + a2 - e1), K := Q(a1, a2, e1). Q is the degree variable, t := log Q, tau := 1 - 1/Q, so Q = 1 is t = tau = 0.

F_g(Q) := Sum_{d>=0} GW_{g,d} Q^d, where GW_{g,d} := GW^{Z, T_CY}_{g, d[P^1]} is the genus-g T_CY-equivariant GW invariant of Z in the localisation definition of [BS] (Def. defn: refined GW invariant, line ~1245), an element of Kw (the T_CY-fixed loci are the a0-fixed loci of the proper M_g(P^1, d), so the localisation definition applies). F_g(0) = GW_{g,0} is the constant-map term.

P_g := [t^{2g}] e1 e2 t^2/(4 sinh(e1 t/2) sinh(e2 t/2)) and beta_n(e1, e2) := (n-1)! [t^n] 1/(4 sinh(e1 t/2) sinh(e2 t/2)) ([BS] eq:betapol, line 5057), so that, for g >= 2,

    beta_{2g-2}(e1, e2) = (2g-3)! P_g/(e1 e2),                                                        (0.1)

both sides being (2g-3)!/(e1 e2) times [t^{2g}] S(e1 t) S(e2 t), S(y) := y/(2 sinh(y/2)). P_g is a symmetric polynomial in e1^2, e2^2 of degree 2g.

### 0.2 Hypothesis LC_CY (s.2)

The statement of LC_CY_PROOF.md as assumed here: on T_CY with a1 a2 A B e1 e2 E != 0, in frame B with the scalar class w = (E, -e1, -e2), Theorem LC of LSP s.0.2 holds with c -> E, leaves lam_j = 1/(E - j psi), Delta_1 = -e1 e2 E^3 t^{-2} times a unit of Kw[[t]] with value 1, everything to all orders in t and z. Exact wording in s.2.

### 0.3 Results

**Proposition 13.2_CY (frame-B formula on T_CY)** [P mod the inputs of LSP Proposition 13.2 (listed in s.3), applied with Kw = Q(a0,a1,a2,e1)]. For g >= 2, F_g lies in Kw(Q), is regular at Q = 0 with power series expansion Sum_d GW_{g,d} Q^d there, and equals the frame-B graph sum (13.1) of LSP with vertex class Omega_w, w = (E, -e1, -e2), and the frame-B factors R' = R_Z h_w^{-1}, where R_Z has the rows h_nw r^pt_{0,i}, nw = (a1, A, -e1, -e2). The graph sum is a finite sum of products of elements of the function field of the double cover C: y^2 = Disc(Q) (s.1.4), and each coloured term is regular over every Q not in {0, 1, oo} with Disc(Q) != 0.

**Corollary R_CY** [C; P mod (MV) and textbook facts given LC_CY; s.4]. Fix g >= 2. In frame B with w = (E, -e1, -e2), on T_CY with a1 a2 A B e1 e2 E != 0 (at the generic point, or at any such specialisation):
 (a) every coloured stable graph of genus g with at least one edge is regular at t = 0 after summing over decorations;
 (b) the single-vertex graph is regular at colour 2, and at colour 1 its polar part is exactly

        gap_g t^{2-2g},     gap_g := (-e1 e2 E^3)^{g-1} Z[g, 0] = beta_{2g-2}(e1, e2) = (2g-3)! P_g/(e1 e2),

     with no polar term of order t^{3-2g}, ..., t^{-1};
 (c) hence the expansion of the frame-B graph sum in Kw((t)) has polar part exactly gap_g t^{2-2g}; by Proposition 13.2_CY the same holds for F_g.

**Theorem Rig (rigidity for the resolved conifold)** [C; P mod (MV), [GP99] and the inputs of Proposition 13.2_CY, given LC_CY; s.5]. For g >= 2, in Kw(Q),

    F_g(Q) - F_g(0) = P_g Li_{3-2g}(Q)/(e1 e2),     i.e.   e1 e2 GW_{g,d} = P_g d^{2g-3}  for every d >= 1.           (0.2)

Since e1 e2 GW_{g,d} is a polynomial in (a0, a1, a2, e1, e2) (Lemma 5.2), (0.2) holds at every point of Lie(T_CY), including a1 = a2, E = 0, e1 = 0. In particular e1 e2 GW_{g,d} depends on (e1, e2) only, symmetrically and through e1^2, e2^2: it lies in Q[eps_-^2, eps_+^2]. Together with genus 0 (Proposition 13.2(a) of LSP: GW_{0,d} = 1/(d^3 e1 e2)) and genus 1 (Lemma 5.7, P mod [BS] prop:ResConifold at genus 1 for d >= 2; d = 1 by direct localisation), this is [BS] Conjecture conj: rigidity (line 1531) for X the resolved conifold, every d >= 1 and G = T_CY, hence for every torus acting through T_CY (s.5.6 says exactly what is and is not implied).

**Corollary 3T (three-torus a0 = 0)** [P mod the existing inputs: Theorem LC (LSP), Theorem L (LAP, mod (MV)), Theorem VI (LP2 s.7.8), Proposition 13.2 of LSP with (d) adapted; s.7]. On a0 = 0 with a1, a2, e1 free (E = a1 + a2 = e1 + e2), Corollary R_CY holds with no new row analysis: the frame-B rows are the fibrewise rows at x = a1/a2 with z -> z/a2, and Theorem L is used at tau_w = -e1/E.

**Remark (flop)** [C; s.6.1]. F_g(oo) = F_g(0): the constant-map term of Z at weights (wt_+, wt_-) equals that of the flopped geometry Z' x C^2 at weights (a1-a2, a2, A, -e1, -e2), (a2-a1, a1, B, -e1, -e2). F_g(0) itself depends on alpha (evidence only, FEAS crit.tested).

### 0.4 What changes relative to [BS]

[BS] proves (0.2) on the subtori a_i = e_j by degeneration (prop:ResConifold, line 1708) and, for general X, the case eps_+ = 0 (prop: rigidity unref limit, line 1543); it verifies (0.2) on T_CY by computer for d <= 3, g <= 5 (remark after line 1777) and states the general case as part of conj: rigidity. Theorem Rig proves (0.2) on all of T_CY for g >= 2 without degeneration, conditionally on LC_CY. s.6.3.

---

## 1. Conventions and the genus-0 data on T_CY

All statements of this section are [P] by hand; each formula was also checked in exact arithmetic (s.9, items B, C).

### 1.1 Weights and fields
As in s.0.1. On T_CY the elements a0, a1, a2, A, B, e1, e2, E, a1 - a2 are non-zero elements of Kw; this is what "a1 a2 A B e1 e2 E != 0" means at the generic point. The map (a0, a1, a2) -> (-a0, a0 + a1, a0 + a2) exchanges wt_+ and wt_-; it sends (a0, a1, a2, A, B, E) to (-a0, B, A, a2, a1, E), so it preserves E and exchanges the pairs (a1, A) <-> (B, a2). It is induced by the automorphism of Z lifting z -> 1/z on P^1 (the bundle O(-1)^2 + O^2 is SL(2)-homogeneous), so every GW_{g,d}, as an element of Kw, is invariant under it. Likewise a1 <-> a2 (exchange of the two O(-1) summands) and e1 <-> e2 are symmetries.

### 1.2 Quantum product, eigenvalues, metric
From LSP Proposition 13.2(b)(i) at indeterminate weights (the genus-0 potential PrepCl + Li_3(Q)/(e1 e2)): H * H = (L_1 + L_2) H - L_1 L_2 with

    L_1 + L_2 = -(a0 + (a1 + a2) Q)/(1 - Q) = E/tau - a0,       L_1 L_2 = -a1 a2 Q/(1 - Q) = a1 a2/tau,

so the L_i are the roots of tau L^2 - (E - a0 tau) L + a1 a2 = 0. The metric inverse is eta^{-1} = e1 e2 [[a0 a1 a2, a1 a2], [a1 a2, -E]] (FEAS (1); consistent with LSP (13.2)). The idempotent norms are

    Delta_i = e1 e2 (a0 a1 a2 + 2 a1 a2 L_i - E L_i^2) = e1 e2 (L_i - L_j)(a1 a2 - E L_i),     (a1 a2 - E L_1)(a1 a2 - E L_2) = a1 a2 A B,      (1.1)

the first equality being LSP (13.2) and the second following from it by the two relations above (checked by hand: expand (L_i - L_j)(a1 a2 - E L_i) with L_j = E/tau - a0 - L_i and L_i L_j = a1 a2/tau). The last identity shows a1 a2 - E L_i != 0 wherever L_i is finite.

### 1.3 The superpotential and the critical points
LSP Lemma 13.3 and FEAS (2) give the phase (in the variable U, theta := U d/dU)

    W(U) = -a2 log U + A log(1 - U) + a1 log(1 - Q U) + (terms independent of U),

so that W_1 := theta W = -a2 - A U/(1 - U) - a1 Q U/(1 - Q U), and the critical equation U W_U = 0 is, after multiplication by (1 - U)(1 - QU),

    P(U) := B Q U^2 - (a0 + Q (a1 - a2)) U - a2 = 0.                                                        (1.2)

Its discriminant is Disc(Q) of s.1.4. The two roots U_1, U_2 (colour 1 = the one with U_1 -> 1 as Q -> 1, FEAS (1)) have L_i = L(U_i) with L(U) := -a1 Q U/(1 - QU) (LSP (b)(i)); L(U) satisfies the quadratic of s.1.2 modulo P(U) (checked, s.9 C). Further

    P(1) = A (Q - 1),       P(1/Q) = a1 (1 - Q)/Q,                                                             (1.3)

so U_i = 1 or Q U_i = 1 happens only over Q = 1 (for Q finite and non-zero). Finally W_2 := theta^2 W = -A U/(1 - U)^2 - a1 Q U/(1 - QU)^2, and

    Delta_i = a1 A e1 e2 W_2(U_i)                                                                             (1.4)

(LSP (13.2), last equality; re-checked modulo P(U), s.9 C). At a0 = 0: P(U) = a2 (x Q U^2 + Q (1 - x) U - 1) with x := a1/a2, which is a2 times LeqU of LSP (C2); and W = a2 (-log U + log(1 - U) + x log(1 - QU)), which is a2 times the fibrewise phase LSP (C1). This is the rescaling used in s.7.

### 1.4 The discriminant
    Disc(Q) := a0^2 + (2 a0 (a1 + a2) + 4 a1 a2) Q + (a1 - a2)^2 Q^2.
 (i) Disc = Q^2 tau^2 (L_1 - L_2)^2 = (a0 + Q(a1 - a2))^2 + 4 a2 B Q, the discriminant of P(U) in U. Proof: tau^2 (L_1 - L_2)^2 = (E - a0 tau)^2 - 4 a1 a2 tau; multiply by Q^2 and expand with Q tau = Q - 1: Q^2 (E - a0)^2 - 4 a1 a2 Q^2 + Q (2 a0 (E - a0) + 4 a1 a2) + a0^2, and (E - a0)^2 - 4 a1 a2 = (a1 - a2)^2, 2 a0 (E - a0) = 2 a0 (a1 + a2). The second form is a direct expansion.
 (ii) disc_Q(Disc) = (2 a0 (a1 + a2) + 4 a1 a2)^2 - 4 (a1 - a2)^2 a0^2 = 16 a1 a2 (a0 + a1)(a0 + a2) = 16 a1 a2 A B.
 (iii) Disc(1) = E^2, Disc(0) = a0^2, and the leading coefficient is (a1 - a2)^2. So Disc has no root over Q in {0, 1} and, as a1 != a2 in Kw, two finite roots.
 (iv) Disc is irreducible over Kw = K(a0). Proof: a quadratic over a field is reducible iff its discriminant is a square; 16 a1 a2 (a0 + a1)(a0 + a2), as a polynomial in a0 over K, has the two simple roots -a1 != -a2, so it is not a square in K[a0], hence not in K(a0) (if f = (p/q)^2 with p, q coprime in the UFD K[a0], then q is a unit).
 (v) The double cover C -> P^1_Q defined by y^2 = Disc(Q) has function field Kw(C) = Kw(Q)[U]/(P(U)) (as 2 B Q U - (a0 + Q(a1 - a2)) = +-y), ramified exactly over the two roots of Disc, and unramified over Q = 0, 1, oo.

### 1.5 Behaviour at Q = 1, 0, oo
 - Q = 1 (tau -> 0): L_1 = E/tau + O(1) and L_2 -> a1 a2/E; Delta_1 ~ -e1 e2 E^3/t^2 and Delta_2(0) = e1 e2 a1 a2 A B/E != 0 (FEAS (1); the leading term of Delta_1 follows from (1.1) with L_1 ~ E/tau, L_2 -> a1 a2/E: e1 e2 (L_1 - L_2)(a1 a2 - E L_1) ~ e1 e2 (E/tau)(-E^2/tau) = -e1 e2 E^3/tau^2, and tau = t + O(t^2); the exact unit is (F1) of LC_CY).
 - Q = 0: the two roots are U_1 -> -a2/a0 and Q U_2 -> a0/B (LSP Lemma 13.3 with (a, b) = (-a2, A) and (a0, a1)); L_1 -> 0, L_2 -> -a0, Delta_1(0) = e(wt_+) = a0 a1 a2 e1 e2, Delta_2(0) = e(wt_-) = -a0 A B e1 e2 (LSP (F8)); all non-zero in Kw because a0 != 0. This is the semisimplicity of Q = 0 at a0 != 0.
 - Q = oo (s := 1/Q -> 0): dividing (1.2) by Q, B U^2 - (a0 s + a1 - a2) U - a2 s = 0, with simple roots U -> (a1 - a2)/B and U -> 0 at s = 0 (simple because a1 != a2); for the second put V := Q U, so B s V^2 - (a0 s + a1 - a2) V - a2 = 0 with simple root V -> -a2/(a1 - a2). Hence U_1(s) and V_2(s) are power series in s over Kw. L_i -> {a1, a2}: indeed L_1 + L_2 -> E - a0 = a1 + a2 and L_1 L_2 -> a1 a2 at tau = 1. By (1.1),

        Delta_1(oo) = e1 e2 (a1 - a2)(a1 a2 - E a1) = -a1 B (a1 - a2) e1 e2,     Delta_2(oo) = e1 e2 (a2 - a1)(a1 a2 - E a2) = a2 A (a1 - a2) e1 e2,   (1.5)

   using a2 - E = -B, a1 - E = -A; both non-zero in Kw. The same values follow from (1.4): W_2(U_1) -> -A U/(1-U)^2 at U = (a1-a2)/B, 1 - U = A/B, giving -(a1 - a2) B/A; and W_2(U_2) -> -a1 V/(1-V)^2 at V = -a2/(a1-a2), 1 - V = a1/(a1 - a2), giving a2 (a1 - a2)/a1. [Numerical confirmation s.9 C'.]

---

## 2. Hypothesis LC_CY (exact statement assumed)

**Hypothesis LC_CY.** Let the weights be generic on T_CY (the generic point of Lie(T_CY), or any specialisation with a1 a2 A B e1 e2 E != 0). Let the frame-B factors be those of LSP s.1.2 (C6)-(C11) for the rows R' = R_Z h_w^{-1} of Proposition 13.2_CY (s.3), with w = (E, -e1, -e2), and define, as in LSP (C13),

    rho := (Delta_1 t^2/(-e1 e2 E^3))^{-1/2},  the branch with rho(0) = 1;   lam_j(psi) := 1/(E - j psi) (j >= 1);   T(psi) := -psi^2/(E - psi);   q_0 := T - psi = -psi E/(E - psi).

Then Theorem LC of LSP s.0.2 holds with c -> E and x c^3 -> e1 e2 E^3 (and with its rational functions of x replaced by elements of Kw): (F1) Delta_1 = -e1 e2 E^3/(t^2 rho^2) with rho a unit of Kw[[t]] with rho(0) = 1, and Delta_2 analytic at t = 0 with Delta_2(0) != 0; (F2) r_{0,1}(z; 0) = E/(E + z); (F3) every Dil_i(d) and every E_ij(m, n) is analytic at t = 0; (L1) [t^0](rho q_t) = q_0 and [t^k](rho q_t) in Lam_{k+1} := span{lam_1, ..., lam_{k+1}} for k >= 1; (L2) for j = 1, 2 and all b, k: [t^k](rho f_{j,b}) in Lam_{k+1}, E_21(b, a) = E_12(a, b), and [t^k] f_{j,b} in Lam_{k+1}; (L3) [t^k](rho^e E_11(psi, psi')) in span{lam_i(psi) lam_l(psi') : (i-1) + (l-1) <= k} for every integer e, and E_11 symmetric; together with the large-z form (LSP Corollary 9.5) used by LSP s.9.6. All to all orders in t and z.

Remarks. (a) The hypothesis is about the decoration-sum factors as formal objects; it does not include Proposition 13.2_CY, which is an input of its own (s.3). (b) Only the following consequences of LC_CY are used below: F1-F3, L1-L3 and Corollary 9.5 (for Corollary R_CY, s.4); and (F1)-(F3) are also used, through Corollary R_CY, in Theorem Rig step (iv). The row analysis of s.5 steps (i)-(iii) does not use LC_CY at all: it uses Proposition 13.2_CY and s.1.

---

## 3. Proposition 13.2_CY: the frame-B formula on T_CY

**Statement.** As in s.0.3. In detail (notation of LSP Proposition 13.2 with Kw := Q(a0, a1, a2, e1) and e2 := E - e1): for g >= 2,

    F_g = Sum_Gamma |Aut Gamma|^{-1} Sum_{c: V -> {1,2}} Sum_dec prod_v Delta_{c_v}^{g_v - 1} Int Omega_w psi^a prod Dil_{c_v}(d) prod_e E_{c c'}(a_e, b_e)      (3.1)

in Kw(Q), with Omega_w := prod_{a=1}^3 Lambda^vee_{g_v}(w_a)/w_a^{g_v}, w = (E, -e1, -e2), and the frame-B factors computed from the rows v_i = (r_{0,i}, r_{1,i}), r_{0,i} := h_w^{-1} h_nw r^pt_{0,i}, r_{1,i} := (L_i + z D_i) r_{0,i}, D_i := d/dt - Delta_i'/(2 Delta_i), r^pt_{0,i} the normalised stationary-phase expansion (LSP (C7)) of e^{W/z} dU/U at U_i, and nw := (a1, A, -e1, -e2).

**Status.** [P mod: [T12] Theorem (2), Prop 8.7(ii), 2.18, 2.19; [GP99]; [CG07] App. 1; [PPZ] Defs 2.2, 2.13, Prop 2.14; [M83]; [FP00] Thm 4 and p. 3; the asymptotic series DLMF 5.11.1 through LSP Lemma 13.3; the further citations of LSP (G3). These are exactly the inputs of LSP Proposition 13.2; the proof below is that proof read at T_CY. This reading has had one referee pass (AI, s.8.4), which found minor items only (applied: the first sentence of the proof, the alternative route, the wording of the Q = 0 regularity); the parts (a)-(c), (e) of LSP had the passes recorded in LSP s.13.2.]

**Proof.** Run LSP Proposition 13.2 parts (a), (b), (c), (e) with the torus T_CY in place of T, i.e. with the field Kw = Q(a0, a1, a2, e1) and the T_CY-linearised bundle Nb = O(-1)_{a1} + O(-1)_{a2} + O_{-e1} + O_{-e2}, e2 = E - e1. Those parts treat the five weights as indeterminates and use no relation among them (FEAS (5); the critic's reading agrees); every hypothesis they need about the weights is a non-vanishing or non-constancy statement, listed below, and each holds in Kw = Q(a0, a1, a2, e1):
 (a) The CohFT: fixed loci of M_{g,n}(Z, d) are those of the proper M_{g,n}(P^1, d), unchanged for T_CY (a T_CY-fixed map has image in the zero section and is fixed by the a0-rotation, since a0 != 0 in Kw). Semisimplicity at Q = 0: L_1(0) - L_2(0) = a0 is a unit of Kw[[Q]] (s.1.5). The field-of-definition and descent paragraph uses only that Kw((Q)) is a field of characteristic 0 of cardinality 2^aleph_0. Flatness (F1)-(F9): (F4) is the divisor equation of P^1; (F8) needs Delta_+(0)/Delta_-(0) = -a1 a2/(A B) to be a non-constant element of Kw, which it is (a non-constant rational function of (a0, a1, a2), transcendental over Q, hence not a root of unity); Delta_+-(0) != 0 by s.1.5.
 (b) The rows: Lemma 13.3 (the Q = 0 constants) needs the limiting pairs (a, b) = (-a2, A) and (a0, a1) with a, b, a + b non-zero: a2, A, a0 and a0, a1, B are non-zero in Kw. Row 2 of the QDE and symplecticity use only the identities of s.1.2 and the invertibility of L_1 - L_2 in Kw[[Q]].
 (c) The constants: ch_{2k-1}(E) != 0 on M_{k,1} ([FP00] Thm 4, [M83]) and the weights wt_+-, all non-zero in Kw, so that prod_a Lambda^vee(w_a)/w_a^g = h_{wt_p}.1 makes sense.
 (e) The frame change R' = R_Z h_w^{-1}: w = (E, -e1, -e2) has non-zero entries, h_w is scalar, and Omega_w = h_w.1 by [M83] and [FP00] s.1.1 as in LSP (c).
 Part (d) of LSP (the passage to the non-semisimple T_F) is not needed: Q = 0 is semisimple on T_CY. (Alternative route, and a cross-check of the reading above: apply LSP Prop. 13.2 (a)-(c), (e) at the five indeterminate weights with the frame w = (a0 + a1 + a2, -e1, -e2), and specialise the resulting identity along e2 -> E - e1 exactly as in LSP (d): the colour-summed total is num/den with den|_{T_CY} != 0, because every per-term denominator is one of the factors listed in the description of the coloured terms below (W_2(U_i), 1 - U_i, 1 - Q U_i, y, B Q, and the constants a0, a1, A, e1, e2, E) and none of them involves E - e1 - e2; on the GW side GW^{T_CY}_{g,d} = GW^T_{g,d}|_{T_CY} by Lemma 5.2 for d >= 1 and, for d = 0, by the constant-map formula of the Remark after Lemma 5.2, whose denominators a0 a1 a2 e1 e2 and a0 A B e1 e2 do not vanish on T_CY. Either route gives (3.1).) What replaces it is the remark that the T_CY-invariants are the restrictions of the T-invariants (s.5.2, Lemma 5.2 and the paragraph after it), which is not used in the proof of (3.1) but identifies F_g with [BS]'s series.
 The description of the coloured terms: each z-coefficient of r^pt_{0,i} is a polynomial in W_3(U_i), W_4(U_i), ... and 1/W_2(U_i) (LSP (C7)); W_k(U) is a polynomial in U/(1 - U) and QU/(1 - QU) over Kw, except for the constant -a2 in W_1; D_i involves dU_i/dQ = -(partial_Q P)/(partial_U P)(U_i), whose denominator is +-y (s.1.4(v)), and Delta_i'/Delta_i; Delta_i = a1 A e1 e2 W_2(U_i) (1.4); h_nw, h_w^{-1}, eta^{-1} are constants of Kw[[z]] (denominators: a1, A, e1, e2, E). So every factor at colour i is an element of Kw(C) whose denominators are products of powers of W_2(U_i) (vanishing iff Delta_i = 0 iff, by (1.1), L_i = L_j or L_i = oo, i.e. Disc = 0 or Q = 1), 1 - U_i and 1 - Q U_i (only over Q = 1 by (1.3)), y (Disc = 0), and B Q (U_i finite iff Q != 0; at Q = oo see s.1.5). The two colours are exchanged by the deck involution of C (the factors at colour 2 are the same expressions evaluated at the other root, LSP (C3), and h_nw is the same scalar at both colours by Lemma 13.3). In genus g only finitely many z-coefficients enter (dilaton degrees and edge exponents are bounded by 3 g_v - 3 + n_v + m_v). So each coloured term is in Kw(C), regular over every Q not in {0, 1, oo} with Disc(Q) != 0, and the colour sum is invariant under the deck involution, hence in Kw(Q); by (b)(iii) the expansion of every coloured term at either of the two (unramified, s.1.4(v)) points of C over Q = 0 lies in Kw[[Q]] (the rows at both colours lie in Kw[[Q]][[z]], and the two points are exchanged by the deck involution together with the colours), so the total is regular at Q = 0 too, and its expansion there is Sum_d GW_{g,d} Q^d by (a), (c), (e). QED

---

## 4. Corollary R_CY

**Statement.** As in s.0.3 (a)-(c). Equivalently: (a) is Theorem R(a) of LSP s.13.1 with the inputs below, (b) is Theorem R(b) with the pole value (4.2), and (c) is the sum.

**Proof.** Follow LSP s.13 step by step. For each step I say whether it is verbatim or what the substitution is.

*Step 0 (LSP s.13.0: domain of Fpot_h, (H0), Lemma 13.0).* Verbatim with c -> E and Q(x) -> Kw: Rr := Kw[[t]], T(psi) = -psi^2/(E - psi), q_0 = T - psi, Dom := q_0 + t Rr[[psi]]. (H0) (dilaton homogeneity, d^n Fpot_h(lambda q)[f] = lambda^{2-2h-n} d^n Fpot_h(q)[f] for lambda in 1 + t Rr) uses only that Omega_w is pulled back along forgetful maps (string and dilaton equations), which holds for every triple w. That q_t and rho q_t lie in Dom is (L1) at k = 0 and (F1) of LC_CY. Lemma 13.0 (range: in a connected stable graph with an edge every g_v <= g - 1) is combinatorial, verbatim.

*Step 1 (Theorem R(i), homogeneity).* At a singular (colour-1) vertex of genus h >= 2 with n half-edges the factor is Delta_1^{h-1} d^n Fpot_h(q_t)[f_1, ..., f_n]. By (H0) with lambda = rho and by the definition of rho (s.2), Delta_1 rho^2 = -e1 e2 E^3/t^2 exactly, so the factor is (-e1 e2 E^3/t^2)^{h-1} d^n Fpot_h(rho q_t)[rho f_1, ..., rho f_n]. Substitution -x c^3 -> -e1 e2 E^3; nothing else changes.

*Step 2 (Theorem R(ii), leaf expansion).* Verbatim from (L1), (L2), (L3) and Corollary 9.5 of LC_CY, through LSP s.9.6: the order-t^k coefficient of every leg, loop and singular-singular edge at the vertex is a combination of leaves lam_j, j >= 1, of total degree D <= k, and Fpot_h(q_0 + Sum_k t^k q^{(k)}) is its power series at q_0 on Dom. So the coefficient of t^{2-2h+k} of the vertex factor is a finite Kw-combination of Xi_h(R) := d_{lam_{j_1}} ... d_{lam_{j_r}} Fpot_h(q_0) with parts j >= 1 and D(R) := Sum (j_i - 1) <= k. Here Xi_h is the leaf-scale-E integral of LAP Cor. L3 with w = (E, -e1, -e2): Xi_h(R) = Xi^E_h(R).

*Step 3 (Theorem R(iii), vanishing: string reduction and Lemma A).* Two inputs, both already in the local-P2 convention w = (E, -eps1, -eps2), E = eps1 + eps2, which is literally the present w with eps -> e:
 - String identity (S_P2) and Theorem C (C_P2), LP2 s.7.1 lines 775-779 [P]: Xi_h(1, R) = (J_R/E^2) Xi_h(R), so Xi_h(1^p, R') = prod_{i<p} (J_{R'} + i) E^{-2p} Xi_h(R'), and Xi_h(1^p) = 0 for p >= 1 (R' empty, J = 0). These use only the string equation and the leaf shapes lam_j = 1/(E - j psi), T = -psi^2/(E - psi) (CRN s.3a with c -> E); no weight enters. Verbatim.
 - M_irr(h) for w = (E, -e1, -e2): Xi_h(R) = 0 for every non-empty multiset R of parts >= 2 with D(R) <= 2h - 3. This is LAP Theorem L via Cor. L3 at tau_w = w_2/E = -e1/E, sigma = -w_2 w_3/E^2 = -e1 e2/E^2 != 0, exactly LAP s.7.5 ("Local P2 x C^2"), and by Cor. L1 it holds at every specialisation with e1 e2 E != 0, the diagonal e1 = e2 included. [P mod (MV)]. The leaf scale E is the first weight of w ("mode B"), as Cor. L3 requires; this is the one place where the CY condition E = e1 + e2 is used (FEAS (3), crit "The CY condition enters for one reason").
 Then, verbatim: strip the parts 1; if R' is empty the term dies by (C_P2); if R' is non-empty and D(R') = D(R) <= k <= 2h - 3 it dies by M_irr(h). A non-zero polar term needs D_v >= 2h_v - 2 at every singular vertex of genus >= 2 with at least one leaf.

*Step 4 (Theorem R(iv), counting).* Verbatim: every vertex of a connected graph with an edge has a half-edge, so r_v >= 1; the t-order of a non-zero term is >= Sum_{v sing, h_v >= 2} (2 - 2h_v + D_v) + 2 (number of singular genus-0 vertices) + 0 >= 0, using (F1) for genus-0 singular vertices (Delta_1^{-1} = O(t^2)) and (F1), (F3) for regular vertices (Delta_2^{h-1}, Dil_2, E_2j analytic, Delta_2(0) != 0). This proves (a). Graphs with an edge use M_irr(h) only for 2 <= h <= g - 1 (Lemma 13.0); at g = 2 no M_irr is used.

*Step 5 (Theorem R(v), the single vertex).* The colour-2 term is regular by (F1), (F3). At colour 1 the factor is (-e1 e2 E^3/t^2)^{g-1} Fpot_g(rho q_t); its order-k coefficient is a combination of Xi_g(R) with R empty only at k = 0, and for 1 <= k <= 2g - 3 of Xi_g(R) with R non-empty, parts >= 1, D(R) <= k, which vanish by Step 3 at h = g. So the polar part is exactly

    (-e1 e2 E^3)^{g-1} Fpot_g(q_0) t^{2-2g} = (-e1 e2 E^3)^{g-1} Z[g, 0] t^{2-2g},                               (4.1)

with no other negative powers. The value: Theorem VI of LP2 s.7.6, proof s.7.8 [P mod (MV) at mu = 1^a and textbook facts], states for w = (E, -eps1, -eps2), E = eps1 + eps2, and every g >= 2,

    (-eps1 eps2 E^3)^{g-1} Z[g, 0] = beta_{2g-2}(eps1, eps2).                                                 (4.2)

This is literally our situation with eps -> e: no homogeneity conversion is needed here (the conversion of LSP s.13.3 (5) is what specialises Theorem VI to the fibrewise convention (x, 1); the local-P2 form is the general one). Its own proof (LP2 s.7.8 Part A) uses homogeneity: Z[g, 0] = E^{-(3g-3)} Z_g(tau) at tau = -e1/E and (-e1 e2 E^3) = sigma E^5, with Z_g(tau) a polynomial in 1/sigma, so (4.2) holds identically in (e1, e2) and at every specialisation with e1 e2 E != 0. With (0.1), gap_g = beta_{2g-2}(e1, e2) = (2g-3)! P_g/(e1 e2). This proves (b).

*Hand check of the identification (0.1) and of gap_g at g = 2, 3.* S(y) = y/(2 sinh(y/2)) = 1 - y^2/24 + 7 y^4/5760 - 31 y^6/967680 + ..., so [t^{2g}] S(e1 t) S(e2 t) = Sum_{i+j=g} s_i s_j e1^{2i} e2^{2j}:
 g = 2: P_2 = (7 e1^4 + 10 e1^2 e2^2 + 7 e2^4)/5760 (s_0 s_2 = 7/5760, s_1^2 = 1/576 = 10/5760), and (2g-3)! = 1, so gap_2 = (7 e1^4 + 10 e1^2 e2^2 + 7 e2^4)/(5760 e1 e2) = beta_2 as displayed in [BS] line 5070; at (e1, e2) = (x, 1) it is (7 + 10 x^2 + 7 x^4)/(5760 x), the fibrewise gap of LSP s.13.3.
 g = 3: P_3 = -(31 e1^6 + 49 e1^4 e2^2 + 49 e1^2 e2^4 + 31 e2^6)/967680 (2 s_1 s_2 = -14/138240 = -98/967680 and 49 = 98/2; s_0 s_3 = -31/967680), and (2g-3)! = 6, so gap_3 = -(31 e1^6 + 49 e1^4 e2^2 + 49 e1^2 e2^4 + 31 e2^6)/(161280 e1 e2) = beta_4 of [BS] line 5070. [Exact re-computation s.9 A.]

*Step 6 (Corollary R of LSP s.13.2 -> (c)).* Each coloured graph term is an element of Kw(C) (s.3), so its expansion at the colour-1 or colour-2 point over Q = 1 lies in Kw((t)) (Q = e^t; Q = 1 is not a branch point, Disc(1) = E^2 != 0), and by LSP (C14) the expansions of the factors are the formal series to which LC_CY refers. Summing (a) and (b) over graphs and colourings gives (c) for the graph sum, and Proposition 13.2_CY transfers it to F_g: F_g in Kw(Q) maps to Kw((t)) under Q = e^t (Q - 1 = t(1 + t/2 + ...) is t times a unit; Q and Disc are units at Q = 1), and the image is the sum of the expansions of the terms. QED

**Specialisation.** At a point of Lie(T_CY) with a1 a2 A B e1 e2 E != 0 (complex or rational weights) the same proof applies, LC_CY being assumed there, with M_irr by LAP Cor. L1/L3 and (4.2) by its polynomiality in 1/sigma. The value gap_g is then a number, non-zero at real weights (all terms of beta_{2g-2}(a, b) have one sign at real a, b != 0, LSP s.13.3) but possibly zero at special complex weights; no proof step uses its non-vanishing.

**Status of Corollary R_CY.** [C] on LC_CY; given it, P modulo (MV) (LAP Theorem L for M_irr(h), 3 <= h <= g; M_irr(2) is proved without (MV), proof_status_mixed_scale.md s.1.4 as cited in LSP s.13.2) and, for the pole value, (MV) at mu = 1^a and the textbook facts listed in LP2 s.7.8 "Inputs". As a statement about the decoration sum it uses no TR, no [BS] input and no finite-order step. Its identification with F_g is Proposition 13.2_CY.

---

## 5. Theorem Rig: rigidity for the resolved conifold

### 5.1 Statement and reduction
For g >= 2, (0.2) in Kw(Q). The proof has four steps: (i) the poles of F_g lie over Q in {1, oo} or Disc(Q) = 0; (ii) no pole on Disc = 0; (iii) no pole at Q = oo; (iv) comparison with Li_{3-2g}. Steps (i)-(iii) use Proposition 13.2_CY and s.1 only; step (iv) uses Corollary R_CY. Throughout, G := F_g - F_g(0) = Sum_{d>=1} GW_{g,d} Q^d; subtracting the constant F_g(0) in Kw changes no pole.

### 5.2 Step (i): the possible poles [P mod Proposition 13.2_CY]

**Lemma 5.1.** F_g lies in Kw(Q), is regular at Q = 0, and its poles lie over Q = 1, Q = oo or the two roots of Disc.
Proof. This is the last paragraph of the proof of Proposition 13.2_CY: each coloured term is an element of Kw(C) regular over every Q not in {0, 1, oo} with Disc(Q) != 0, the colour sum is invariant under the deck involution hence lies in Kw(Q), and the expansion of every coloured term at either of the two (unramified) points of C over Q = 0 lies in Kw[[Q]] (LSP (b)(iii) at both colours, Q = 0 semisimple because a0 != 0; Q = 0 is unramified by s.1.4(v)). A rational function of Q over Kw which, as an element of Kw(C), is regular over a point of P^1 that is unramified, is regular there. QED

What is actually proved about the rows (to be read against the wording "algebraic functions of Q"): the rows are formal series in z whose coefficients are elements of Kw(C); the genus-g graph sum uses finitely many of them; the identity (3.1) is in Kw(Q), and the equality of its Q = 0 expansion with Sum GW_{g,d} Q^d is Proposition 13.2_CY. No analytic continuation in Q is used anywhere; everything is algebraic in Kw(C).

### 5.3 Step (ii): no pole on the discriminant

**Lemma 5.2 (polynomiality) [P mod [GP99]].** For d >= 1 and all g >= 0, e1 e2 GW^T_{g,d} lies in Q[a0, a1, a2, e1, e2], homogeneous of degree 2g, symmetric under e1 <-> e2, under a1 <-> a2 and under (a0, a1, a2) -> (-a0, a0 + a1, a0 + a2). The T_CY-invariant is its restriction e2 -> a0 + a1 + a2 - e1. Hence, on T_CY, GW_{g,d} lies in (1/(e1 e2)) Q[a0, a1, a2, e1], and therefore in K[[a0]], K = Q(a1, a2, e1): 1/e2 = 1/((a1 + a2 - e1) + a0) expands in K[[a0]] because a1 + a2 - e1 is a non-zero element of K.
Proof. This is the argument of FB B5 with O(-3) replaced by O(-1)_{a1} + O(-1)_{a2}, and it is [BS] lem: shape of refGW series (line 1472) for the proper M := Mbar_g(X, d), X the resolved conifold (equirigid for d >= 1, [BS] line 752-754): every compact curve in X lies in the zero section, so Mbar_g(X, d) = Mbar_g(P^1, d) as a space, proper. The T-fixed loci of Mbar_g(Z, d) are those of M (LSP Prop. 13.2(a)); for d >= 1 and C connected, H^0(C, f^*O(-1)) = 0 (a section vanishes on the non-contracted components, where f^*O(-1) has negative degree, and is constant on the contracted ones, which meet them), so R^1_a := R^1 pi_* f^* O(-1)_a is a vector bundle of rank d + g - 1 on M, and R pi_* f^* Nb = -R^1_{a1} - R^1_{a2} + (O - E^vee)_{-e1} + (O - E^vee)_{-e2}. Thus e_T(R pi_* f^* Nb)^{-1} = e_T(R^1_{a1}) e_T(R^1_{a2}) e_T(E^vee_{-e1}) e_T(E^vee_{-e2})/(e1 e2) is (1/(e1 e2)) times an honest T-equivariant class on M (e_T(R^1_a) depends on a0 through the base rotation); the localisation definition equals the T-equivariant integral over the proper M by [GP99] (applied to a generic one-parameter subgroup with the same fixed loci, both sides rational in the weights, as in FB B5), and the equivariant pushforward lands in H_T^*(pt) = Q[a0, a1, a2, e1, e2]. The degree is 2(d + g - 1) + 2g - 2(d + g - 1) = 2g (virtual dimension 2d + 2g - 2 of M). The symmetries are those of s.1.1. Restriction to T_CY: the fixed loci are the same, e_T of the fixed and moving parts restrict, and no denominator vanishes on T_CY (the weights of the virtual normal bundles are of the form (fibre weight at p_+-) + j a0/k and a0/k, none a multiple of a0 + a1 + a2 - e1 - e2); so GW^{T_CY}_{g,d} = GW^T_{g,d}|_{T_CY}, which is also [BS] lem: GW contravariance. QED

Remark. LSP Prop. 13.2(d) states the weaker fact GW^T_{g,d} in Q[a1, a2, e1, e2][1/(a1 a2 e1 e2)][[a0]], by expanding e_T^{-1} in a0. Either statement suffices for Lemma 5.3: what is used is only GW_{g,d} in K[[a0]] for d >= 1. For d = 0, GW_{g,0} = Sum_{p = +,-} Int_{M_g} prod_{a=1}^5 Lambda^vee(w^{(p)}_a)/w^{(p)}_a (localisation on M_g x {p}) is a rational function with denominators a0 a1 a2 e1 e2 and a0 A B e1 e2; it is a constant in Q and plays no role in this step.

**Lemma 5.3 (the small root).** Disc has a unique root Q_1 in K[[a0]] with Q_1 = O(a0); it has a0-adic valuation 2:

    Q_1 = -a0^2/(4 a1 a2) + (a1 + a2) a0^3/(8 a1^2 a2^2) - (5 a1^2 + 6 a1 a2 + 5 a2^2) a0^4/(64 a1^3 a2^3) + O(a0^5).

Proof. Modulo a0, Disc = 4 a1 a2 Q + (a1 - a2)^2 Q^2 has the simple root Q = 0 (the derivative there is 4 a1 a2 != 0 in K), which lifts uniquely to a root in K[[a0]] (simple-root lifting over the complete local ring K[[a0]], by successive approximation). Writing Q_1 = c_2 a0^2 + c_3 a0^3 + c_4 a0^4 + ... and expanding: order a0^2: 1 + 4 a1 a2 c_2 = 0; order a0^3: 4 a1 a2 c_3 + 2 (a1 + a2) c_2 = 0; order a0^4: 4 a1 a2 c_4 + 2 (a1 + a2) c_3 + (a1 - a2)^2 c_2^2 = 0, giving c_4 = -[4 (a1 + a2)^2 + (a1 - a2)^2]/(64 a1^3 a2^3) = -(5 a1^2 + 6 a1 a2 + 5 a2^2)/(64 a1^3 a2^3). (No a0^1 term is possible: the a0^1 coefficient of Disc(Q_1) is 4 a1 a2 c_1.) The other root is -4 a1 a2/(a1 - a2)^2 + O(a0), a unit. [Exact series check to O(a0^8), s.9 D.] QED

**Lemma 5.4 (no pole on Disc = 0) [P mod Lemma 5.1, Lemma 5.2].** G, and hence F_g, has no pole at the roots of Disc.
Proof. By Lemma 5.1 and the irreducibility of Disc (s.1.4(iv)), in the UFD Kw[Q] we can write

    G = N/((1 - Q)^n Disc^m),    N in Kw[Q], n, m >= 0, with Disc not dividing N,                                   (5.1)

(G is regular at Q = 0 and has no other finite poles; m is then the order of the pole of G along Disc, and we must show m = 0). Suppose m >= 1.
 1. Coefficients of N. The identity N = (1 - Q)^n Disc^m G holds in Kw[[Q]] (expand the rational functions at Q = 0; G is the series Sum_{d>=1} GW_{g,d} Q^d by Proposition 13.2_CY). Embed Kw = K(a0) into the field K((a0)) of formal series with finitely many negative powers of a0; the identity then holds in K((a0))[[Q]]. On the right, (1 - Q)^n Disc^m has coefficients in K[a0] and G has coefficients in K[[a0]] (Lemma 5.2), so the right side lies in K[[a0]][[Q]]. Hence every coefficient of N lies in K(a0) ∩ K[[a0]], and N is in K[[a0]][Q].
 2. Evaluation at Q_1. Let ev: K[[a0]][[Q]] -> K[[a0]] be substitution Q -> Q_1. It is well defined and a ring homomorphism because Q_1 has valuation 2 > 0: for a series Sum_d c_d Q^d with c_d in K[[a0]], the terms c_d Q_1^d have valuation >= 2d, so the sum converges a0-adically, and products are respected (K[[a0]] is complete and the substitution of a topologically nilpotent element is multiplicative). On the polynomial N, ev is the ordinary evaluation N(Q_1). Applying ev to the identity of step 1: N(Q_1) = (1 - Q_1)^n Disc(Q_1)^m ev(G) = 0, since Disc(Q_1) = 0 and m >= 1.
 This is the point flagged in the brief: the series G in K[[a0]][[Q]] can be evaluated at Q_1 only because val(Q_1) > 0, and the value of the rational function G at Q_1 is not used; what is used is that the polynomial identity N = (1 - Q)^n Disc^m G of K((a0))[[Q]] is preserved by ev, and N(Q_1) is the evaluation of a polynomial.
 3. Divisibility. N is in K(a0)[Q] and vanishes at Q_1 in the extension field K((a0)) of K(a0). The minimal polynomial of Q_1 over K(a0) is Disc/(a1 - a2)^2 (irreducible over K(a0) by s.1.4(iv), with Q_1 as a root), so Disc divides N in K(a0)[Q] = Kw[Q]. This contradicts (5.1). Hence m = 0. QED

### 5.4 Step (iii): regularity at Q = oo [P mod Proposition 13.2_CY]

**Lemma 5.5.** F_g is regular at Q = oo (as a1 != a2 in Kw).
Proof. Q = oo is unramified for C (s.1.4(v): the leading coefficient (a1 - a2)^2 of Disc is non-zero), so the two points of C over s := 1/Q = 0 are the two colours, and regularity has to be shown at both of them; this is what follows. In s, by s.1.5 the critical points U_1(s) and V_2(s) := Q U_2 are power series in s over Kw with U_1(0) = (a1 - a2)/B and V_2(0) = -a2/(a1 - a2). The ingredients of the factors at colour 1 are W_k(U_1), polynomials in U_1/(1 - U_1) and Q U_1/(1 - Q U_1) = 1/(s/U_1 - 1); both are power series in s (1 - U_1(0) = A/B != 0, and s/U_1 - 1 -> -1). At colour 2, U_2/(1 - U_2) = s V_2/(1 - s V_2) and Q U_2/(1 - Q U_2) = V_2/(1 - V_2), power series in s since 1 - V_2(0) = a1/(a1 - a2) != 0. Further W_2(U_i)(s = 0) != 0 by (1.5) and (1.4), so 1/W_2(U_i) is a power series in s; L_i = L(U_i) = -a1 Q U_i/(1 - Q U_i) is a power series in s with L_1(0) = a1, L_2(0) = a2; d/dt = Q d/dQ = -s d/ds maps power series in s to power series in s; Delta_i = a1 A e1 e2 W_2(U_i) is a power series in s with non-zero value (1.5), so Delta_i^{-1} and Delta_i'/Delta_i are power series in s; h_nw, h_w^{-1}, eta^{-1} are constants. Hence every z-coefficient of r_{0,i}, r_{1,i}, every Dil_i(d), E_ij(a, b) and Delta_i^{+-1} is a power series in s, so every coloured term of (3.1), and F_g, is a power series in s: F_g has no pole at Q = oo. (The normalised expansion (C7) has leading coefficient 1 and involves no square root, LSP s.1.2; so no branch issue arises at s = 0.) QED

The values at s = 0 (used only in s.6.1): L_i -> (a1, a2), Delta_1(oo) = -a1 B (a1 - a2) e1 e2, Delta_2(oo) = a2 A (a1 - a2) e1 e2, by (1.5).

### 5.5 Step (iv): comparison with the polylogarithm and conclusion

**Lemma 5.6 (Li_{-n}) [P].** For n >= 1, Li_{-n}(Q) := Sum_{d>=1} d^n Q^d = A_n(Q)/(1 - Q)^{n+1} with A_n in Z[Q] of degree n and A_n(0) = 0 (A_0 = Q, A_{m+1} = Q[(1 - Q) A_m' + (m + 1) A_m], from Li_{-m-1} = Q d/dQ Li_{-m}). Hence Li_{-n} is a rational function with its only pole at Q = 1, of order n + 1, vanishing at Q = 0 and at Q = oo (degree n < n + 1), and Li_{-n}(1/Q) = (-1)^{n+1} Li_{-n}(Q) (induction from n = 1, where Li_{-1}(1/Q) = Q/(Q - 1)^2 = Li_{-1}(Q); Q d/dQ anti-commutes with Q -> 1/Q, so Li_{-n-1}(1/Q) = -(Q d/dQ)[Li_{-n}(1/Q)]. For n = 0 instead Li_0(1/Q) = -1 - Li_0(Q), so Li_0 does not vanish at oo; this is why g >= 2 is needed here). Its expansion at t = 0, Q = e^t: Li_0(e^t) = e^t/(1 - e^t) = -1/t - 1/2 - t/12 + ..., so Li_{-n}(e^t) = (d/dt)^n Li_0(e^t) = (-1)^{n+1} n! t^{-n-1} + (power series in t), with no other negative power. For n = 2g - 3 (odd) the polar part is (2g-3)! t^{2-2g}. [Checks s.9 E, E'.]

**Proof of Theorem Rig.** Let g >= 2 and put

    D(Q) := F_g(Q) - F_g(0) - P_g Li_{3-2g}(Q)/(e1 e2)   in Kw(Q).

 - Poles of D at finite Q != 1: F_g has none by Lemmas 5.1 and 5.4 (regular at Q = 0, no pole on Disc = 0); Li_{3-2g} has none (Lemma 5.6).
 - Pole at Q = oo: none, by Lemma 5.5 for F_g and Lemma 5.6 for Li.
 - Pole at Q = 1: under Q = e^t, D maps to an element of Kw((t)) whose polar part is, by Corollary R_CY(c), gap_g t^{2-2g} minus, by Lemma 5.6 and (0.1), (P_g/(e1 e2)) (2g-3)! t^{2-2g} = beta_{2g-2}(e1, e2) t^{2-2g} = gap_g t^{2-2g}: zero. Here it matters that Corollary R_CY gives the polar part of F_g exactly, with no terms t^{3-2g}, ..., t^{-1}: these are what Theorem R(v) kills by M_irr(g) (s.4 Step 5), and Lemma 5.6 says Li has none either. A rational function of Q whose image in Kw((t)) has no negative powers is regular at Q = 1 (if it had a pole of order r >= 1 at Q = 1, its image would have the leading term (Q - 1)^{-r} -> t^{-r} times a unit).
 So D is a rational function of Q over Kw with no pole on P^1, hence a constant c in Kw, and D(0) = F_g(0) - F_g(0) - 0 = 0 gives c = 0. This is (0.2) in Kw(Q); comparing coefficients of Q^d, e1 e2 GW_{g,d} = P_g d^{2g-3} in Kw for every d >= 1.
 Extension to the excluded loci. Both sides of e1 e2 GW_{g,d} = P_g d^{2g-3} are polynomials in (a0, a1, a2, e1) on Lie(T_CY) (left: Lemma 5.2; right: P_g in Q[e1^2, e2^2] with e2 = a0 + a1 + a2 - e1). Two polynomials equal in Kw are equal, so the identity holds at every point of Lie(T_CY), including a1 = a2 (where Lemma 5.5's argument would fail), E = 0 (outside frame B), e1 = 0 or e2 = 0 (where GW_{g,d} itself is undefined but e1 e2 GW_{g,d} is not), and a0 = 0 (where Q = 0 is not semisimple). Nothing else is needed: the proof is an identity in the function field, and the excluded loci are reached by polynomial continuation, not by taking limits of the graph sum. QED

**Genus 0 and 1 (for the comparison with [BS]).**
 - g = 0: GW_{0,d} = 1/(d^3 e1 e2) at every weight (LSP Prop. 13.2(a), from [FP00] p. 3: M_{0,0}(Z, d) = M_{0,0}(P^1, d) x C^2 with obstruction bundle of rank equal to the virtual dimension). P_0 = 1, so e1 e2 GW_{0,d} = P_0 d^{-3}. [P mod [FP00].]
 - g = 1: **Lemma 5.7** [P mod [BS] prop:ResConifold at genus 1 for d >= 2, [GP99], [M83]; at d = 1 without [BS], see the direct localisation below]. e1 e2 GW_{1,d} = P_1/d = -(e1^2 + e2^2)/(24 d) on T_CY. Proof. In genus 1, Lambda^vee(-e) = -e - lambda_1 and lambda_1^2 = 0 ([M83]: c(E) c(E^vee) = 1), so e1 e2 GW_{1,d} = e1 e2 N_d + (e1 + e2) L_d (= e1 e2 N_d + E L_d on T_CY), an identity on the 5-torus, with N_d := Int_{[M_1(P^1,d)]^vir} e(R^1_{a1}) e(R^1_{a2}) in Q (degree 0) and L_d := Int_{[M_1(P^1,d)]^vir} e(R^1_{a1}) e(R^1_{a2}) lambda_1 in Q[a0, a1, a2]_1, a linear form (Lemma 5.2's argument: degrees 2d and 2d + 1 against virtual dimension 2d). L_d does not involve e1, e2. By the symmetries of s.1.1, L_d = c (a0 + a1 + a2) for some c in Q: a1 <-> a2 gives L_d = c_0 a0 + c_1 (a1 + a2), and invariance under (a0, a1, a2) -> (-a0, a0 + a1, a0 + a2) gives c_0 = c_1. On the subtorus a0 = 0, a_i = e_i of [BS] prop:ResConifold, e1 e2 GW_{1,d} = [t^2] S(e1 t) S(e2 t)/d = -(e1^2 + e2^2)/(24 d) = e1 e2/(12 d) - (e1 + e2)^2/(24 d); there (e1 + e2)(a0 + a1 + a2) = (e1 + e2)^2, so comparing with e1 e2 N_d + c (e1 + e2)^2 (the quadratic forms e1 e2 and (e1 + e2)^2 are linearly independent over Q) gives N_d = 1/(12 d) and c = -1/(24 d). Hence e1 e2 GW_{1,d} = e1 e2/(12 d) - (e1 + e2)(a0 + a1 + a2)/(24 d) on the 5-torus, which on T_CY is e1 e2/(12 d) - E^2/(24 d) = P_1/d, using E = e1 + e2 = a0 + a1 + a2. QED. (A direct computation of the single linear form L_d for all d, or a genus-1 analogue of Proposition 13.2_CY, would remove the [BS] input; for d = 1 it is done next, for d >= 2 it is not in the material.)

 *The case d = 1 by direct localisation* [P mod the localisation formula (L) below ([GP99] on the proper Mbar_g(P^1, 1), in the standard form; its ingredients are listed and its normalisation is confirmed by g = 0), [M83], and Int_{Mbar_{1,1}} lambda_1 = Int_{Mbar_{1,1}} psi_1 = 1/24]. For d = 1 the T-fixed loci of Mbar_g(P^1, 1) are indexed by g_+ + g_- = g: the zero section covered once, with a contracted component of genus g_p attached at p (absent if g_p = 0), fixed locus Mbar_{g_+,1} x Mbar_{g_-,1}. Collecting the factors (vertex at p: H^1 and H^0 of T_Z on the contracted component, prod_a Lambda^vee(w^{(p)}_a)/prod_a w^{(p)}_a; node: prod_a w^{(p)}_a from the normalisation sequence and 1/(w^{(p)}_1 - psi) from the smoothing, w^{(p)}_1 = +-a0 the tangent weight of P^1 at p; edge: H^0 = H^1 = 0 for O(-1), H^0(P^1, O_{-e_i}) giving 1/(e1 e2), and of H^0(P^1, T P^1) with weights (a0, 0, -a0) the part not cancelled by the automorphisms of the edge fixing its nodes, i.e. 1/w^{(p)}_1 for each node),

    e1 e2 GW_{g,1} = Sum_{g_+ + g_- = g} V_+(g_+) V_-(g_-),   V_p(0) := 1,   V_p(h) := Int_{Mbar_{h,1}} prod_{a=1}^5 Lambda^vee_h(w^{(p)}_a)/(w^{(p)}_1 (w^{(p)}_1 - psi)),            (L)

 with Lambda^vee_h(w) = Sum_i (-1)^i lambda_i w^{h-i}; at g = 0 it gives GW_{0,1} = 1/(e1 e2), the value of Prop. 13.2(a). At h = 1, Lambda^vee_1(w) = w - lambda_1 and lambda_1^2 = 0, so with Pi_p := prod_{a=1}^5 w^{(p)}_a and 1/(w_1 (w_1 - psi)) = w_1^{-2} (1 + psi/w_1 + ...) the degree-1 part of the integrand is (Pi_p/w_1^2) [psi/w_1 - lambda_1 Sum_a 1/w_a], and Int psi = Int lambda_1 = 1/24 gives

    V_p(1) = (Pi_p/(24 w_1^2)) [1/w_1 - Sum_{a=1}^5 1/w_a] = -(Pi_p/(24 w_1^2)) Sum_{a=2}^5 1/w_a.

 At p_+ (w = (a0, a1, a2, -e1, -e2), Pi_+ = a0 a1 a2 e1 e2): V_+(1) = -(1/(24 a0)) [(a1 + a2) e1 e2 - a1 a2 (e1 + e2)]. At p_- (w = (-a0, B, A, -e1, -e2), Pi_- = -a0 A B e1 e2, w_1^2 = a0^2): V_-(1) = (1/(24 a0)) [(A + B) e1 e2 - A B (e1 + e2)]. Since A + B - a1 - a2 = 2 a0 and A B - a1 a2 = a0 (a0 + a1 + a2),

    e1 e2 GW_{1,1} = V_+(1) + V_-(1) = e1 e2/12 - (a0 + a1 + a2)(e1 + e2)/24    on the 5-torus,

 i.e. N_1 = 1/12 and L_1 = -(a0 + a1 + a2)/24 without [BS]; on T_CY this is e1 e2 GW_{1,1} = P_1. So Lemma 5.7 holds at d = 1 modulo (L) and [M83] only, and the [BS] input is used for d >= 2 only. [E: the referee's independent direct-localisation test at d = 1 gave the same e1 e2 GW_{1,1} on T_max, and e1 e2 GW_{2,1} = P_2 on T_CY; both agree with the genus <= 2 part of [BS]'s display at line 1694 off T_CY.]

### 5.6 What exactly is implied for [BS] conj: rigidity

[BS] conj: rigidity (line 1531): if beta is a rigid curve class in X and the G-action on Z is CY, then eps1 eps2 GW(Z, G; beta) := eps1 eps2 Sum_{g>=0} GW_g(Z, G; beta) lies in Q[[eps_-^2, eps_+^2]], i.e. depends only on (eps1, eps2), symmetrically and evenly. The sum is over genera at fixed beta; beta = 0 is excluded (the series GW_g(Z, G) at line 1277 sums over non-zero classes, and "rigid" is a property of a non-zero class).

Implied by Theorem Rig with Lemma 5.7 and genus 0: for X the resolved conifold, every d >= 1, and G = T_CY (the maximal CY torus of Z = X x C^2 inside T_max, [BS] line 1703): e1 e2 GW_g(Z, T_CY; d[P^1]) = P_g d^{2g-3} for all g >= 0, so the conjecture holds, with the explicit value Sum_g P_g d^{2g-3}/(e1 e2) = (1/d)(2 sinh(d e1/2))^{-1} (2 sinh(d e2/2))^{-1}, i.e. [BS] eq:resconGWref on all of T_CY, as [BS]'s remark after line 1777 predicts. By contravariance ([BS] lem: GW contravariance, Lemma 5.2's last paragraph), the same holds for every subtorus of T_CY and, via R_G = R_T^W, for every connected G whose maximal torus acts through a subtorus of T_CY (in particular for every CY subtorus of T_max). Not implied: the conjecture for group actions on Z that are not conjugate into T_max, and the constant-map term d = 0, which is not part of the conjecture and does depend on alpha (s.6.2). Status: for g >= 2 [C; P mod (MV), [GP99] and the inputs of Proposition 13.2_CY]; g = 1 additionally mod [BS] prop:ResConifold for d >= 2 (d = 1 by direct localisation, Lemma 5.7); g = 0 mod [FP00].

---

## 6. Remarks

### 6.1 The flop identity F_g(oo) = F_g(0) [C; P mod Theorem Rig and the Q = oo values below]
By (0.2) and Lemma 5.6, F_g(oo) = F_g(0) + 0. The left side is the value at s = 0 of the graph sum (Lemma 5.5). At s = 0 all d/dt-derivatives vanish, so r_{1,i} = L_i r_{0,i}, and r^pt_{0,i}(z; oo) is the normalised stationary-phase expansion of the limiting phases: at colour 1, W -> (a1 - a2) log U + A log(1 - U) + const (as a1 log(1 - QU) = a1 log(-QU) + a1 log(1 - 1/(QU))), a Beta-type phase with (a, b) = (a1 - a2, A) in the notation of LSP Lemma 13.3, whose limit row is h_{(a2 - a1, -A, B)} (B = a1 - a2 + A); at colour 2, with V = QU, W -> -a2 log V + a1 log(1 - V) + const, (a, b) = (-a2, a1), limit row h_{(a2, -a1, a1 - a2)}. Multiplying by h_nw = h_{(a1, A, -e1, -e2)} (pairs v, -v cancel): the rows at s = 0 are diag(h_{wt'_-}, h_{wt'_+}) with

    wt'_- := (a2 - a1, a1, B, -e1, -e2),     wt'_+ := (a1 - a2, a2, A, -e1, -e2),

and Delta_1(oo) = e(wt'_-), Delta_2(oo) = e(wt'_+) by (1.5). These are the tangent weights at the two fixed points of the flop Z' = X' x C^2 (X' the other small resolution: the cone xy = zw has coordinate weights {a1, A} and {a2, B}, a1 + A = a2 + B = E; X pairs them as normal weights (a1, a2) at p_+ and (B, A) at p_-, X' pairs them as (a2, A) and (a1, B), with P^1' tangent weight +-(a1 - a2)); the CY condition is preserved. So the graph sum at s = 0 is the Q = 0 (constant-map) value of the frame-B sum of Z', i.e. F_g^{Z'}(0) = Sum_{p} e(wt'_p)^{g-1} Int_{M_g} prod_a Lambda^vee(w'^{(p)}_a)/(w'^{(p)}_a)^g by LSP Prop. 13.2(c) read for Z'. Hence

    Sum_{p = +,-} Int_{M_g} prod_{a=1}^5 Lambda^vee(w^{(p)}_a)/w^{(p)}_a   is the same for (wt_+, wt_-) and for (wt'_+, wt'_-),

an identity among constant-map lambda-class integrals, for which no direct proof is in the material. Status of the limiting rows: the Lemma 13.3 argument (stationary phase at a Beta-type phase, DLMF 5.11.1) applied at s = 0; derived here by hand, not re-written in full. [E: FEAS crit.tested, flop_const.log, g = 2..6 at two CY samples; fails off CY. The critical points at s = 0, the W_2 limits (1.5) and the flop weights wt'_+- were also checked exactly by the referee (evidence).]

### 6.2 F_g(0) depends on alpha [E]
GW_{g,0} = Sum_p Int_{M_g} prod_a Lambda^vee(w^{(p)}_a)/w^{(p)}_a is a rational function of degree 2g - 2 in the weights. On T_F it is -(e1^{2g-2} + e2^{2g-2}) |B_{2g} B_{2g-2}|/(4g (2g-2) (2g-2)!) (LSP s.13.3 evidence paragraph, proved there from [GP99] and [FP00] Thm 4). On T_CY it depends on alpha: FEAS crit.tested ("F_g(0) itself does depend on alpha, g = 2, 3, 4, three alpha samples"). This is consistent with [BS]: d = 0 is not in conj: rigidity. Theorem Rig says nothing about F_g(0) beyond s.6.1.

### 6.3 Relation to the existing partial results of [BS]
 - prop:ResConifold (line 1708): (0.2) for all d > 0 and all genera on subtori with a_i = e_j, by degeneration and the relative invariants of local curves. Theorem Rig contains it for g >= 2 (and reproves it without degeneration), uses it at g = 1 for d >= 2 (Lemma 5.7; d = 1 is done by direct localisation there), and does not touch g = 0.
 - prop: rigidity unref limit (line 1543): the case eps_+ = 0, i.e. e2 = -e1, E = 0, on which GW_g = (e1 e2)^{g-1} N_{g,d} by [M83] and the properness of M_g(X, d). On T_CY this is the hyperplane E = 0 (a0 + a1 + a2 = 0 too), outside frame B; Theorem Rig reaches it by the polynomial continuation of s.5.5, and there (0.2) reads GW_{g,d} = -P_g(e, -e) d^{2g-3}/e^2 = the e^{2g-2} term of -(1/d) (2 sinh(d e/2))^{-2} (e := e1 = -e2; P_g(e, -e) = [t^{2g}] S(e t)^2 = e^{2g} [t^{2g}] S(t)^2, and (2 sinh(d e/2))^{-2} = (d e)^{-2} S(d e)^2), the genus-g part of the unrefined conifold formula [BS] eq: res conifold unref (line 1780).
 - The computer verification of (0.2) on T_CY for d <= 3 up to genus 5 ([BS] remark after line 1777) becomes a theorem, conditionally on LC_CY.
 - The remark after conj: rigidity (line 1565): the conjecture implies that any admissible 2-torus recovers the full series. For the conifold this is now a theorem on T_CY.

---

## 7. Corollary 3T: the three-torus a0 = 0

**Setting.** T_3 := {a0 = 0} ∩ T_CY: weights (0, a1, a2, e1, e2) with a1 + a2 = e1 + e2 = E, a1, a2, e1 free. Function field K_3 := Q(a1, a2, e1). Non-degeneracy a1 a2 e1 e2 E != 0 (A = a2, B = a1). This contains the fibrewise torus T_F (a_i = e_i) of LSP.

**Statement.** Corollary R_CY holds on T_3 (at the generic point and at every specialisation with a1 a2 e1 e2 E != 0), and Proposition 13.2 of LSP holds on T_3 with h_nw = h_{(a1, a2, -e1, -e2)} kept.

**Proof.** (1) Rows. At a0 = 0, nw = (a1, a2, -e1, -e2) and w = (E, -e1, -e2), so h_nw h_w^{-1} = h_{(a1, a2, -e1, -e2, -E, e1, e2)} = h_{(a1, a2, -E)}, independent of e1, e2 (pairs v, -v cancel; h_{-v} = h_v^{-1}). The point-frame rows depend only on the phase W, which at a0 = 0 is a2 times the fibrewise phase at x = a1/a2 (s.1.3). Since e^{W/z} = e^{(W/a2)/(z/a2)} and the stationary-phase expansion (C7) is normalised to leading coefficient 1, r^pt_{0,i}(z) at (a1, a2) equals the fibrewise r^pt_{0,i}(z/a2) at x = a1/a2; and h_{(a1, a2, -E)}(z) = h_{(x, 1, -(1+x))}(z/a2), which is the fibrewise frame-B constant h_nw h_w^{-1}|_{T_F} = h_{(x, 1, -c)}, c = 1 + x (LSP (C8), with h_nw = 1 on T_F). The eigenvalues L_i = a2 L^F_i(x) have degree 1 in (a1, a2) (tau L^2 - E L + a1 a2 = a2^2 [tau (L/a2)^2 - c (L/a2) + x], c := 1 + x, the fibrewise quadratic of LSP (C5)). The metric inverse eta^{-1} = e1 e2 [[0, a1 a2], [a1 a2, -E]] is not a multiple of the fibrewise eta_F^{-1} = [[0, x^2], [x^2, -c x]] (LSP (C4)): its entries have different degrees, e1 e2 a1 a2 = (e1 e2/x) a2^2 x^2 (degree 2) and -e1 e2 E = (e1 e2/x) a2 (-c x) (degree 1). D_i = d/dt - Delta_i'/(2 Delta_i) is scale-invariant. Since r_{0,i}(z) = r^F_{0,i}(z/a2) and hence, by (C9), r_{1,i}(z) = (L_i + z D_i) r_{0,i}(z) = a2 [L^F_i + (z/a2) D_i] r^F_{0,i}(z/a2) = a2 r^F_{1,i}(z/a2), the pairing <v_i(z), v_j(w)> = eta^{01} [r_{0,i}(z) r_{1,j}(w) + r_{1,i}(z) r_{0,j}(w)] + eta^{11} r_{1,i}(z) r_{1,j}(w) equals (e1 e2/x) a2^3 times the fibrewise pairing at (z/a2, w/a2) (the eta^{01} term carries a2^2 x^2 times one factor a2, the eta^{11} term a2 (-c x) times a2^2), consistent with Delta_i = e1 e2 (2 a1 a2 L_i - E L_i^2) = (e1 e2/x) a2^3 Delta^F_i(x), Delta^F_i = x (2 x L^F_i - c (L^F_i)^2) (LSP (4.1)). Hence, by LSP (C10) (Dil_i(d) = (-1)^d r_{0,i}[d-1]), Dil_i(d) = a2^{1-d} Dil^F_i(d), and by LSP (C11) (E_ij(m, n) is a signed sum of <v_i[p], v_j[q]> over p + q = m + n + 1), E_ij(m, n) = (e1 e2/x) a2^{2-m-n} E^F_ij(m, n). So the frame-B rows on T_3 are the fibrewise frame-B rows at x = a1/a2 with z -> z/a2 (and r_{1,i} multiplied by a2), and every Delta_i, Dil_i(d), E_ij(m, n) is the fibrewise one times the explicit monomial in a2 and e1 e2/x just stated. [Checked: FEAS report 1, T0 "the a0 = 0 rescaling to fibrewise rows", and s.9 F for the constant.]
 (2) Theorem LC on T_3. Theorem LC of LSP holds at x = a1/a2 (x != 0, -1 since a1 a2 E != 0), for all orders in t and z. Under z -> z/a2 the leaves 1/(c - j psi) become a2/(E - j psi) = a2 lam^E_j (psi -> psi/a2, a2 c = E), Lam_{k+1} is preserved, rho = (Delta_1 t^2/(-e1 e2 E^3))^{-1/2} equals the fibrewise rho(x) (since Delta_1 t^2/(-e1 e2 E^3) = Delta^F_1 t^2/(-x c^3) by (1)), and the memberships F1-F3, L1-L3 and Corollary 9.5 hold with the same coefficients up to monomials in a2. So Hypothesis LC_CY holds on T_3, with LSP Theorem LC as its proof.
 (3) Corollary R_CY on T_3. Apply s.4 verbatim: the vertex class is Omega_w with w = (E, -e1, -e2), which is not the fibrewise class (that would be (E, -a1, -a2)/a2 rescaled) but is the class to which M_irr applies at tau_w = -e1/E (LAP Cor. L3, s.7.5) and for which Theorem VI gives (4.2). The decoupling makes this consistent: the rows and factors carry (a1, a2) only, the vertex class carries (e1, e2) only, and Theorem R needs no relation between them because each Xi_h(R) vanishes separately (M_irr is a vanishing statement, monomial by monomial in the leaf coefficients). The pole value is (-e1 e2 E^3)^{g-1} Z[g, 0] = beta_{2g-2}(e1, e2) = (2g-3)! P_g/(e1 e2).
 (4) Proposition 13.2 on T_3. Q = 0 is not semisimple at a0 = 0 (both roots of P degenerate, s.1.5), so (3.1) on T_3 is obtained as in LSP Prop. 13.2(d), by the ring map a0 -> 0 (then e2 -> a1 + a2 - e1) applied to the identity at indeterminate weights: the GW side lies in Q[a1, a2, e1, e2][1/(a1 a2 e1 e2)][[a0]] for d >= 1 (Lemma 5.2, or LSP (d)) and GW_{g,0} = Sum_p Int_{M_g} prod_{a=1}^5 Lambda^vee(w^{(p)}_a)/w^{(p)}_a is regular at a0 = 0: only the factors Lambda^vee(a0)/a0 (at p_+) and Lambda^vee(-a0)/(-a0) (at p_-) have a pole there, a simple one, and the two residues cancel. Indeed the a0^{-1} coefficient of the p_+ term is (-1)^g Int_{M_g} lambda_g Lambda^vee(a1) Lambda^vee(a2) Lambda^vee(-e1) Lambda^vee(-e2)/(a1 a2 e1 e2) (from Lambda^vee(a0)/a0 = Sum_i (-1)^i lambda_i a0^{g-i-1}, whose a0^{-1} term is (-1)^g lambda_g/a0; the other four factors at p_+ do not involve a0), and that of the p_- term is minus the same (from Lambda^vee(-a0)/(-a0) = Sum_i (-1)^i lambda_i (-a0)^{g-i-1}, whose a0^{-1} term is (-1)^g lambda_g/(-a0) = -(-1)^g lambda_g/a0, times the value at a0 = 0 of the factor Lambda^vee(B) Lambda^vee(A) Lambda^vee(-e1) Lambda^vee(-e2)/(A B e1 e2), which is regular there with A -> a2, B -> a1). Moreover GW^{T_3}_{g,0} in the [BS] localisation definition equals this value at a0 = 0, by the same contravariance/localisation argument as for d >= 1 (Lemma 5.2): the T_3-fixed locus of Mbar_g(Z, 0) = Mbar_g x Z is the proper Mbar_g x P^1 (the zero section is fixed pointwise by T_3, its normal weights a1, a2, -e1, -e2 are non-zero), with virtual class e(E^vee (x) T_{P^1}) ∩ [Mbar_g x P^1] (the T_3-fixed part of the obstruction bundle E^vee (x) T_Z) and virtual normal bundle Nb (x) (O - E^vee); the same integral taken T-equivariantly over the proper Mbar_g x P^1 lies in Q[a0, a1, a2, e1, e2][1/(a1 a2 A B e1 e2)] (its denominators are the T-weights of Nb at p_+-), restricts to T_3 by contravariance, and equals Sum_p prod_{a=1}^5 Lambda^vee(w^{(p)}_a)/w^{(p)}_a by localisation at the two T-fixed points of P^1 ([GP99]; at p the factor Lambda^vee(+-a0)/(+-a0) is e(E^vee (x) T_pP^1)/e(T_pP^1)); on the graph side the factors are regular along a0 = 0 at generic Q (roots distinct: Disc|_{a0=0} = 4 a1 a2 Q + (a1 - a2)^2 Q^2 != 0 for generic Q; U_i not in {0, 1, 1/Q}; W_2, Delta finite non-zero; eta^{-1} polynomial; h_nw -> h_{(a1, a2, -e1, -e2)} finite), so the colour-summed total num/den has den|_{a0 = 0} != 0 and the specialisation is the point-frame sum with the rows of (1). This is LSP (d) with "h_nw -> 1" replaced by "h_nw -> h_{(a1, a2, -e1, -e2)}"; nothing else in (d) changes. Then (e) verbatim. QED

**Status.** P modulo the existing inputs: Theorem LC (LSP, with its inputs DLMF 5.11.1 and the citations of LSP s.0.5), Theorem L (LAP, mod (MV)), Theorem VI (LP2 s.7.8, mod (MV) at mu = 1^a), Proposition 13.2 of LSP with (d) adapted as in (4) (mod its literature inputs). No new row analysis and no new vanishing. Rigidity on T_3 is not a separate statement: it is the restriction of Theorem Rig by polynomiality (s.5.5), and the a0-adic step of s.5.3 is not available on T_3 itself (a0 = 0 there), so Corollary 3T alone does not give (0.2).

---

## 8. Status table, gaps found, open points

### 8.1 Status
    Statement                              Status
    s.1 genus-0 data, Disc, limits          P (hand; exact checks s.9 B, C, C', D)
    Proposition 13.2_CY                     P mod the inputs of LSP Prop. 13.2 (list in s.3), read at T_CY; one AI referee pass of this reading (s.8.4), minor items applied
    Corollary R_CY (a), (b), (c)            C on LC_CY; then P mod (MV) (Theorem L; Theorem VI at mu = 1^a) and textbook facts
    Lemma 5.1 (poles)                       P mod Proposition 13.2_CY
    Lemma 5.2 (polynomiality)               P mod [GP99] (= [BS] lem: shape of refGW series for the proper M_g(X, d))
    Lemma 5.3 (small root)                  P
    Lemma 5.4 (no pole on Disc)             P mod Lemmas 5.1, 5.2
    Lemma 5.5 (Q = oo)                      P mod Proposition 13.2_CY; needs a1 != a2 (true in Kw)
    Lemma 5.6 (Li)                          P
    Theorem Rig, g >= 2                     C on LC_CY; then P mod (MV), [GP99], inputs of Proposition 13.2_CY
    Lemma 5.7 (g = 1)                       d = 1: P mod the degree-1 localisation formula (L), [M83]; d >= 2: P mod [BS] prop:ResConifold at g = 1, [GP99], [M83]
    g = 0                                   P mod [FP00] p. 3 (LSP Prop. 13.2(a))
    Remark 6.1 (flop)                       C; limiting rows derived by hand, Lemma 13.3 argument at s = 0 not re-written; critical points, W_2 limits, flop weights checked exactly by the referee (E)
    Remark 6.2 (F_g(0) alpha-dependent)     E only
    Corollary 3T                            P mod existing inputs (Theorem LC, Theorem L, Theorem VI, LSP Prop. 13.2 with (d) adapted)

### 8.2 Gaps found in the proposed route, and how they are handled
 1. Genus 1 is not covered by the graph argument (Proposition 13.2 and Theorem R are for g >= 2), while [BS] conj: rigidity concerns the full genus sum. Handled by Lemma 5.7, at the price of one [BS] input (prop:ResConifold at genus 1). A direct evaluation of the linear form L_d = -(a0 + a1 + a2)/(24 d) would remove it; this is done at d = 1 (Lemma 5.7, direct localisation), so the [BS] input is used for d >= 2 only.
 2. The brief suggested K[[a0]] with K = Q(a1, a2, e1) as the coefficient ring, and asked whether F_g Disc^m (1 - Q)^n Q^p has coefficients there when m, n, p are chosen over Kw. Answer (Lemma 5.4 step 1): p = 0 since G is a power series; m, n are integers and (1 - Q)^n Disc^m has coefficients in K[a0]; so N has coefficients in K(a0) ∩ K[[a0]], computed inside K((a0))[[Q]]. The evaluation at Q_1 is legitimate because val(Q_1) = 2 > 0 and because it is applied to the polynomial identity in K[[a0]][[Q]], not to the rational function. No gap remains in this step.
 3. The polynomiality input: the brief asked whether only GW_{g,d} in K[[a0]] is available. The stronger statement e1 e2 GW_{g,d} in Q[a0, a1, a2, e1, e2] (d >= 1) is available (Lemma 5.2 = FB B5 argument = [BS] lem: shape of refGW series), and it is what makes the extension to the excluded loci a one-line polynomial identity. LSP Prop. 13.2(d) states only the K[[a0]] form; both suffice for Lemma 5.4.
 4. Lemma 5.5 needs a1 != a2; this is automatic in Kw, and the locus a1 = a2 is recovered by polynomiality. Similarly E = 0, e_i = 0, a0 = 0.
 5. The brief's formula "QU = 1 forces a1(1 - Q) = 0" (FEAS crit) is correct: P(1/Q) = a1 (1 - Q)/Q (a first hand computation of the writer gave B(1 - Q)/Q, which the exact check s.9 C rejected; (1.3) is the corrected form). The brief's "(a1 a2 - E L_i) = 0 only at tau = 0" is correct by (1.1).
 6. Nothing in s.5 steps (i)-(iii) uses LC_CY; the only use of LC_CY in Theorem Rig is the polar part at Q = 1 through Corollary R_CY. Conversely Corollary R_CY does not use steps (i)-(iii). So a failure of LC_CY would leave Lemmas 5.1-5.6 intact and would remove only the identification of the polar part.
 7. Proposition 13.2_CY is "LSP Prop. 13.2 read at T_CY". The reading in s.3 lists the points where weights enter; it is not a line-by-line re-derivation of LSP (a)-(c), (e), which were refereed in their original form only.

### 8.3 Open points
 - A human reading of this file (one AI referee pass is recorded in s.8.4), in particular of Lemma 5.4 (the a0-adic step), of the claim of Lemma 5.1 that the only possible poles are over {0, 1, oo} and Disc = 0, and of s.3.
 - A conceptual proof of the flop identity of s.6.1, and a direct proof of Lemma 5.7 without [BS] for d >= 2 (d = 1 is done).
 - Whether the same route (graph regularity + Disc-type lemma + regularity at the other boundary points of the Q-line) applies to other equirigid targets; not examined.

### 8.4 Referee record
One hostile referee pass (AI), verdict: sound with minor fixes, no major items. The eight minor items (s.7 proof (1): the metric is not a rescaled fibrewise metric entrywise, replaced by the entrywise and pairing scalings; s.7 proof (4): the a0 = 0 regularity of GW_{g,0} argued by the residue cancellation and the proper T_3-fixed locus; Lemma 5.7: the 5-torus factor is (e1 + e2), and the d = 1 case is proved by direct localisation; s.6.3: "the e^{2g-2} term of"; Lemma 5.2: the coefficient clause deleted; s.3: the first sentence of the proof and the alternative route through LSP (d); Lemma 5.1, s.3 and Lemma 5.5: regularity at both (unramified) points of C over Q = 0 and Q = oo; Remark 6.1: the referee's exact checks) are applied; each item that changes mathematics was re-derived by hand before being applied, and none changed a statement of s.0.3. A fresh direct-localisation test at d = 1, g <= 2 by that referee confirmed e1 e2 GW_{1,1} = P_1 and e1 e2 GW_{2,1} = P_2 on T_CY and the [BS] display (line 1694) off T_CY (evidence only; the g = 1 part is now also a hand derivation, Lemma 5.7). Nothing in this file has been read by a human.

---

## 9. Checks run (evidence only; no proof step depends on them)

Script: <user home>/AppData/Local/Temp/claude/G--My-Drive-Projects-Old-refGW/37f24460-3601-4c85-b123-5156e19eb475/scratchpad/cy_torus/rig_writer/checks.py, log checks.log in the same folder. Python 3.14, exact rational arithmetic (fractions), random rational weights, no symbolic package. All items report True.
 A. S(y) = y/(2 sinh(y/2)) = 1 - y^2/24 + 7 y^4/5760 - 31 y^6/967680 + 127 y^8/154828800 - ...; P_2, P_3 as polynomials; (2g-3)! P_g/(e1 e2) = beta_2, beta_4 as displayed in [BS] line 5070; gap_2 at (x, 1) = (7 + 10 x^2 + 7 x^4)/(5760 x).
 B. 50 random points: Disc = Q^2 [(E - a0 tau)^2 - 4 a1 a2 tau]; disc_Q(Disc) = 16 a1 a2 A B; Disc(1) = E^2; (a1 a2 - E L_1)(a1 a2 - E L_2) = a1 a2 A B.
 C. 30 random points, exact polynomial arithmetic modulo P(U): U W_U (1 - U)(1 - QU) = P(U); disc_U(P) = Disc; L(U) = -a1 Q U/(1 - QU) satisfies tau L^2 - (E - a0 tau) L + a1 a2 = 0; e1 e2 (2L - E/tau + a0)(a1 a2 - E L) = a1 A e1 e2 W_2(U) = e1 e2 (a0 a1 a2 + 2 a1 a2 L - E L^2); P(1) = A (Q - 1); P(1/Q) = a1 (1 - Q)/Q; P|_{a0 = 0} = a2 (x Q U^2 + Q(1 - x) U - 1).
 C'. Floating point at Q = 10^6, 10^9: U_1 -> (a1 - a2)/B, Q U_2 -> -a2/(a1 - a2), L_i -> a1, a2, Delta_1 -> -a1 B (a1 - a2) e1 e2, Delta_2 -> a2 A (a1 - a2) e1 e2.
 D. Simple-root lifting iteration in K[[a0]] to O(a0^8) at 5 random (a1, a2): Disc(Q_1) = 0 to that order, Q_1 = -a0^2/(4 a1 a2) + (a1 + a2) a0^3/(8 a1^2 a2^2) - (5 a1^2 + 6 a1 a2 + 5 a2^2) a0^4/(64 a1^3 a2^3) + ... (the a0^4 coefficient corrected a first hand value).
 E. Li_{-n} = A_n/(1 - Q)^{n+1}, A_n(0) = 0, deg A_n = n (A_1 = Q, A_2 = Q + Q^2, A_3 = Q + 4Q^2 + Q^3, A_4 = Q + 11 Q^2 + 11 Q^3 + Q^4); Li_{-n}(1/Q) = (-1)^{n+1} Li_{-n}(Q) for n = 1..7 at 10 random Q each; agreement with the series at Q = 1/3 to 10^{-30}; Li_0(1/Q) = -1 - Li_0(Q).
 E'. Expansion of Li_{-n}(e^t) at t = 0 from e^t/(1 - e^t) = -1/t - 1/2 - t/12 + ...: polar part exactly (-1)^{n+1} n! t^{-n-1} for n = 1..6 (so (2g-3)! t^{2-2g} for n = 2g - 3), no other negative powers; Li_{-1}(e^t) = t^{-2} - 1/12 + O(t^2), Li_{-3}(e^t) = 6 t^{-4} + 1/120 + O(t^2).
 F. Power sums Sum v^{1-2k}, k = 1, 2, 3, at 20 random points: nw + (E, -a1, -A) = (E, -e1, -e2) (the constant cancellation of FEAS (2)); at a0 = 0, nw - w = (a1, a2, -E); at Q = oo, (a2 - a1, -A, B) + nw = (a2 - a1, a1, B, -e1, -e2) and (a2, -a1, a1 - a2) + nw = (a1 - a2, a2, A, -e1, -e2).
 Not re-run here, cited as evidence from FEAS: graph regularity on T_CY at g = 2..5; direct localisation GW_{g,d} = P_g d^{2g-3}/(e1 e2) for g <= 7 at d = 1, g <= 5 at d <= 3 (crit.tested); rational reconstruction in Q at g = 2, 3 (edge graphs have Disc poles, orders 4 and 8, cancelling in the total; totals (1 - Q)^{-2}, (1 - Q)^{-4} with numerator degrees 2, 4, so regular at oo); flop identity g = 2..6; alpha-dependence of F_g(0).
