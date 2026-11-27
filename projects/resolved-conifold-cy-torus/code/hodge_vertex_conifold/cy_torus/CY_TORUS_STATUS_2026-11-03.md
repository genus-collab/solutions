# Resolved conifold x C^2 on the CY 4-torus: status, 2026-11-03

Folder: hodge_vertex_conifold/cy_torus/. Files: LC_CY_PROOF.md (Theorem LC_CY, generic-weight leaf structure; Theorems Q_CY, CF_CY; Corollary 12.1), RIGIDITY_PROOF.md (Proposition 13.2_CY, Corollary R_CY, Theorem Rig, Lemma 5.7, Corollary 3T, flop remark), scripts/ (writers', referees' and fixers' code and logs). Feasibility study: ../cy_torus_explore/ (STATUS_2026-11-02.md, feasibility_result.json).

Notation: E = a0+a1+a2 = e1+e2, A = a0+a2, B = a0+a1, u = (E-2a1)(E-2a2), tau = 1-1/Q, t = log Q; (ND): a1 a2 A B e1 e2 E != 0.

## Result

Conditional chain, with every input named in the files:

1. Theorem LC_CY [P mod the analytic inputs of the fibrewise Theorem LC: LSP Prop 3.15 and Lemma 13.3 (Laplace + Stirling)]. At every weight point with (ND), in frame B with class w = (E,-e1,-e2) (rows h_{(a1,A,-E)} r^pt), the leaf structure holds with leaves 1/(E - j psi): F1-F3, L1-L3 as in Theorem LC with c -> E, x c^3 -> e1 e2 E^3. The CY condition E = e1 + e2 is NOT used (Remark 0.1).
2. Corollary R_CY [P mod (MV) and textbook facts, given 1; the CY condition enters here, so that w is a sum-zero triple with first weight E]: on T_CY every coloured stable graph with an edge is regular at Q = 1; the vertex graph has polar part exactly beta_{2g-2}(e1,e2) t^{2-2g} = (2g-3)! P_g/(e1 e2) t^{2-2g}.
3. Proposition 13.2_CY [P mod the inputs of LSP Prop 13.2, read at indeterminate weights]: identification of the graph sum with F_g on T_CY; easier than on T_F since Q = 0 is semisimple (a0 != 0).
4. Theorem Rig [P mod (MV), GP, the inputs of 3, given 1]: for g >= 2, F_g - F_g(0) = P_g Li_{3-2g}(Q)/(e1 e2) on T_CY, i.e. e1 e2 GW_{g,d} = P_g d^{2g-3} for d >= 1. Together with genus 0 and genus 1 (Lemma 5.7: d = 1 by hand, d >= 2 mod BS prop:ResConifold at genus 1), this is BS conj: rigidity for the resolved conifold (d >= 1; the d = 0 term depends on alpha and is excluded, as in BS). Steps: poles of F_g only in {0, 1, oo, Disc = 0}; no pole on Disc by the a0-adic evaluation at the small root Q_1 (valuation 2); regular at Q = oo (a1 != a2); exact polar part at Q = 1 from 2; Li_{3-2g} has the same polar part and vanishes at 0 and oo; the difference is constant and zero.
5. Corollary 3T: on a0 = 0 (a1, a2, e1 free) the frame-B rows are the fibrewise rows at x = a1/a2, z -> z/a2, so Corollary R holds there from the existing Theorem LC by rescaling.
6. Remark: F_g(oo) = F_g(0) (flop identity of constant-map Hodge integrals); conditional.

## Refereeing (all AI; no human has read these texts)

- LC_CY_PROOF.md: hand referee line by line against LSP (sound with minor fixes, no major; six items, two of them typos in intermediate displays, applied) and independent numerical referee (own F_p code; KEY FORMULA, F1-F3, L1 k <= 13, L2, L3 k <= 9, flow, QDE, symplecticity, regular colour, a0 = 0 rescaling, non-CY points; nine weight points; four controls fail as required; no refutation).
- RIGIDITY_PROOF.md: hostile hand referee (sound with minor fixes, no major; eight items applied), with a fresh direct localisation at d = 1 giving e1 e2 GW_{1,1} = P_1, e1 e2 GW_{2,1} = P_2 on T_CY and the BS display off T_CY.
- Main-session hand checks: constant cancellation h_nw h_{(E,-a1,-A)} = h_w; the a0-adic lemma; P_2/(e1 e2) = beta_2.
- Earlier evidence (cy_torus_explore): graph regularity exact through genus 4 at generic and special CY weights; direct localisation e1 e2 GW_{g,d} = P_g d^{2g-3} through g = 7 at d = 1; rational reconstruction in Q at g = 2, 3 (edge graphs have Disc poles that cancel in the total); flop identity g = 2..6.

## Not done / open

- No human referee. Teleman, PPZ, GP, Mumford/FP, (MV), Laplace/Stirling inputs cited, not re-proved.
- Barnes-type analytic inputs: Prop 3.15 and Lemma 13.3 of LSP are used through rescaling and a multiset identity only.
- Off T_CY the regularity fails (polar terms appear from genus 2 in the vertex graph and genus 3 in edge graphs), so the extension stops at T_CY.
- Lemma 5.7 for d >= 2 at genus 1 uses one BS input; a direct computation of L_d = -E/(24 d) would remove it.
- Whether the same method applies to local P2 on its CY subtorus was not examined.
