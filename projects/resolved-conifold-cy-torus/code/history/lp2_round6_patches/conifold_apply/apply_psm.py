# Moved on 2026-11-27 from localp2/hodge_vertex/lp2_explore/round6_scripts/conifold_apply/apply_psm.py (local P2 project) into the resolved-conifold
# project, because it edits only files of this project. One-off patch script, already applied: a record, not to be rerun.
# Targeted in-place edits of proof_status_mixed_scale.md, conifold_regularity_note.md, lemmaA/RESUME.md.
import os
B=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "hodge_vertex_conifold", "")
def load(p): return open(B+p,encoding='ascii',newline='').read()
def save(p,s): open(B+p,'w',encoding='ascii',newline='').write(s)
def rep(s,old,new):
    n=s.count(old)
    assert n==1,(n,old[:80])
    return s.replace(old,new)

# ---------- proof_status_mixed_scale.md
p="proof_status_mixed_scale.md"; s=load(p)
s=rep(s,"""> mixed edges in graphValueC); both are part of what is now proved.

""","""> mixed edges in graphValueC); both are part of what is now proved.
>
> **Addendum 3, 2026-11-02 (round 6; supersedes the status statements below where they differ).** Lemma A = M_irr(h), for every h >= 2 and non-empty R, is now
> lemmaA/LEMMA_A_PROOF.md, Theorem L (s.7.4 there for the conifold convention), P modulo the Marino-Vafa formula (MV) and textbook facts, with no TR and no computer step; two round-6 referees, "sound with minor
> fixes", no major item. So Corollary R (leaf_proof/leaf_structure_proof.md s.13.2) is unconditional modulo (MV), textbook facts and, only to identify the graph sum with F_g, the inputs of Proposition 13.2 there; genus 3 and 4 are no longer
> computer-assisted (their exact data are independent confirmations). Graphs with an edge need Lemma A only for h <= g - 1; the single-vertex graph also needs h = g. The value of the single pole is Theorem VI of the local-P2 draft (draft s.7.8;
> further inputs LLZ eq. (5) and the vertex-operator form of psi, psi^*; one referee pass; specialisation in leaf_structure_proof.md s.13.3 (5)), and [BS] Prop. 3.1 is now an independent check. Lemmas B and D, and "mod P2" below, are no longer needed for anything.

""")
s=rep(s,"""   all genera, general tau ............ OPEN; reduced (1.7) to the degree bound for C_tau at h integer values of tau.
""","""   all genera, general tau ............ OPEN; reduced (1.7) to the degree bound for C_tau at h integer values of tau. [Superseded by Addendum 3: proved modulo (MV), all h, all r, all tau with sigma != 0.]
""")
s=rep(s,"""   M_irr(h):  Xi_h(R) = 0 for every multiset R of integers >= 2 with D(R) <= 2h-3   (so r <= 2h-3, 2r <= J <= 2h-3+r).
This is exactly the input of the leaf-structure proposition (section 3).
""","""   M_irr(h):  Xi_h(R) = 0 for every non-empty multiset R of integers >= 2 with D(R) <= 2h-3   (so 1 <= r <= 2h-3, 2r <= J <= 2h-3+r).
The empty multiset is excluded: Xi_h(empty) = Z[h,0] is the gap.
This is exactly the input of the leaf-structure proposition (section 3).
""")
s=rep(s,"""Lemma A (= M_irr(h), the vanishing direction of Conjecture M after the string/dilaton reduction). For h >= 4 and every tau:
Xi_h(R) = 0 for every multiset R of integers >= 2 with D(R) = Sum (j_i - 1) <= 2h-3.
""","""[Status 2026-11-02: Lemma A is proved modulo (MV), lemmaA/LEMMA_A_PROOF.md Theorem L; Addendum 3. The route lists below are history.]
Lemma A (= M_irr(h), the vanishing direction of Conjecture M after the string/dilaton reduction). For h >= 4 and every tau:
Xi_h(R) = 0 for every non-empty multiset R of integers >= 2 with D(R) = Sum (j_i - 1) <= 2h-3.
""")
s=rep(s,"""Theorem R (conditional). Fix g >= 2. Assume F1, F3, Lemma C through order t^(2h-3) for every 2 <= h <= g (only at the
psi-orders that occur in genus-g graphs), and M_irr(h) for every 2 <= h <= g. Then every coloured stable graph of genus g with
at least one edge is regular at t = 0 after summing over decorations, and the single-vertex graph has the single pole
(-x c^3)^(g-1) F_g(q_0) t^(2-2g). Since the conifold F_g is P_g Li_{3-2g}(Q)/(eps1 eps2) up to Q-independent terms (the definition of gap_g in F4)
and the only singular term of Li_{3-2g}(e^t) is (2g-3)! t^(2-2g), its polar part is gap_g t^(2-2g), so this pole is gap_g t^(2-2g); so F4 holds in genus g as a consequence
(F4 is verified independently only for g <= 4), which determines F_g(q_0) = gap_g/(-x c^3)^(g-1). Assuming the hypotheses
for all g gives the all-genus theorem.
""","""Theorem R (conditional). Fix g >= 2. Assume F1, F3, Lemma C through order t^(2h-3) for every 2 <= h <= g (only at the
psi-orders that occur in genus-g graphs). (a) If M_irr(h) holds for every 2 <= h <= g - 1, every coloured stable graph of genus g with
at least one edge is regular at t = 0 after summing over decorations (no vertex of such a graph has genus g: leaf_proof/leaf_structure_proof.md Lemma 13.0). (b) If M_irr(g) holds, the single-vertex graph has the single pole
(-x c^3)^(g-1) F_g(q_0) t^(2-2g). By Theorem VI of the local-P2 draft (draft s.7.8; specialisation h_w, w = (c, -x, -1) written out in leaf_structure_proof.md s.13.3 (5)), (-x c^3)^(g-1) F_g(q_0) = beta_{2g-2}(x, 1) = gap_g, the t^(2-2g) coefficient of P_g Li_{3-2g}(Q)/(eps1 eps2);
so F4 at n = 0 holds in genus g, with no input from F_g itself. Independent check: the conifold F_g is P_g Li_{3-2g}(Q)/(eps1 eps2) up to Q-independent terms ([BS] Prop. 3.1) and the only singular term of Li_{3-2g}(e^t) is (2g-3)! t^(2-2g), so by (a) and Proposition 13.2 the pole must be gap_g t^(2-2g). Assuming the hypotheses
for all g gives the all-genus theorem.
""")
s=rep(s,"""Delta_1^(h-1) d^(n_v) F_h(q_t)[f_1, ..., f_{n_v}], F_h the genus-h potential of the CohFT Omega_w in dilaton-shifted coordinates
(the vertex class is the t-independent Omega_w in frame B, all t-dependence sits in Delta, Dil, E). The dilaton equation gives
F_h(lambda q) = lambda^(2-2h) F_h(q), so this equals""","""Delta_1^(h-1) d^(n_v) F_h(q_t)[f_1, ..., f_{n_v}], F_h the genus-h potential of the CohFT Omega_w
(the vertex class is the t-independent Omega_w in frame B, all t-dependence sits in Delta, Dil, E). The dilaton equation, in the form (H0) on the domain q_0 + t C[[t]][[psi]] of leaf_structure_proof.md s.13.0 (needed because rho q_t has psi-coefficient -rho), gives
F_h(lambda q) = lambda^(2-2h) F_h(q), so this equals""")
s=rep(s,"""F3) >= 0. (v) The single-vertex graph has no edges and no loops: its only term with R empty is the order k = 0, which is
(-x c^3)^(g-1) F_g(q_0); all its other polar terms have R non-empty with D <= k <= 2g-3 and vanish. Comparing the sum over
all graphs with the known polar part of F_g identifies this term with gap_g t^(2-2g).  QED
""","""F3) >= 0. This proves (a). (v) The single-vertex graph has no edges and no loops: its only term with R empty is the order k = 0, which is
(-x c^3)^(g-1) F_g(q_0); all its other polar terms have R non-empty with D <= k <= 2g-3 and vanish by M_irr(g). The value is Theorem VI, as in (b).  QED
""")
s=rep(s,"""Obstruction: none found. What is unconditional:
""","""Obstruction: none found. [Status 2026-11-02: see Addendum 3; the three items below are superseded.] What is unconditional:
""")
save(p,s)

# ---------- conifold_regularity_note.md
p="conifold_regularity_note.md"; s=load(p)
s=rep(s,"""F5 [data, genus 3] Inside a stable graph""","""    [2026-11-02: n >= 1 in all genera is Theorem C below. n = 0 in all genera is Theorem VI of the local-P2 draft (hodge_vertex/lp2_explore/LP2_REGULARITY_PROOF_draft.md s.7.8), (-x c^3)^{h-1} Z[h, 0] = beta_{2h-2}(x, 1) = gap_h; specialisation in leaf_proof/leaf_structure_proof.md s.13.3 (5).]

F5 [data, genus 3] Inside a stable graph""")
s=rep(s,"""## 4. What is proven, what is verified, what is open

""","""## 4. What is proven, what is verified, what is open

- Status 2026-11-02 (round 6; supersedes everything below where it differs): Lemma A = M_irr(h), for every h >= 2 and non-empty R, is lemmaA/LEMMA_A_PROOF.md, Theorem L (s.7.4 for this convention), P modulo the Marino-Vafa formula (MV) and textbook facts; two round-6 referees, sound with minor fixes, no major item.
  Corollary R (leaf_proof/leaf_structure_proof.md s.13.2) is therefore unconditional in every genus modulo (MV), textbook facts and, only to identify the graph sum with F_g, the inputs of Proposition 13.2 there. Genus 3 and 4 are no longer computer-assisted. Graphs with an edge need Lemma A only for h <= g - 1;
  the value of the single pole is Theorem VI of the local-P2 draft (one referee pass), with Brini-Schuler Prop. 3.1 as an independent check.
""")
s=rep(s,"""- String/dilaton reduction: the vanishing direction of M at genus h is equivalent to M_irr(h): Xi_h(R) = 0 for all multisets R of
  integers >= 2 with D(R) <= 2h-3.
""","""- String/dilaton reduction: the vanishing direction of M at genus h is equivalent to M_irr(h): Xi_h(R) = 0 for all non-empty multisets R of
  integers >= 2 with D(R) <= 2h-3 (at R empty, Xi_h is the gap).
""")
s=rep(s,"""M_irr(h) established.
- h = 2:""","""M_irr(h) established. [2026-11-02: every h >= 2, lemmaA/LEMMA_A_PROOF.md Theorem L, modulo (MV); the items below are now independent confirmations, h = 2 being also MV-free.]
- h = 2:""")
s=rep(s,"""Open, as the remaining lemmas of the all-genus proof.
- A: M_irr(h) for h >= 5.
""","""Open, as the remaining lemmas of the all-genus proof. [2026-11-02: A is proved modulo (MV), section 4; B and D are not needed.]
- A: M_irr(h) for h >= 5.
""")
s=rep(s,"""graph of genus g with at least one edge is regular at the conifold point, and the single-vertex graph carries exactly the gap, provided
M_irr(h) holds for 2 <= h <= g. M_irr(2) is proved and M_irr(3), M_irr(4) are exact computations, so the theorem holds through genus 4 (genus 2
without any computer input, genus 3 and 4 computer-assisted, independently of the direct graph sums); in every genus it is conditional only
on Lemma A, M_irr(h) for h >= 5. Lemmas B and D (the topological-recursion route to A) are not needed for the reduction.
""","""graph of genus g with at least one edge is regular at the conifold point (using M_irr(h) for 2 <= h <= g - 1 only), and the single-vertex graph carries exactly the gap (using M_irr(g) too, with the value from Theorem VI of the local-P2 draft).
[2026-11-02: M_irr(h) holds for every h >= 2 by lemmaA/LEMMA_A_PROOF.md Theorem L, modulo (MV); see section 4 and leaf_proof/leaf_structure_proof.md s.13.2. M_irr(2) is also proved without (MV), and the exact computations of M_irr(3), M_irr(4) are independent confirmations, not proof steps.]
Lemmas B and D (the topological-recursion route to A) are not needed for the reduction.
""")
save(p,s)

# ---------- lemmaA/RESUME.md
p="lemmaA/RESUME.md"; s=load(p)
s=rep(s,"""# RESUME NOTE: Lemma A / all-genus regularity, state at wrap-up (2026-10-31, about 21:40)

""","""# RESUME NOTE: Lemma A / all-genus regularity, state at wrap-up (2026-10-31, about 21:40); resume point updated 2026-11-02 (round 6)

## RESUME POINT, 2026-11-02 (round 6; supersedes everything below where it differs)
- Lemma A is now lemmaA/LEMMA_A_PROOF.md (Theorem L; s.7.4 for the conifold convention), P modulo (MV) and textbook facts. Frozen copy refereed in round 6: LEMMA_A_PROOF_round5.md (two referees, sound with minor fixes, no major item; their minor items are about labels, stale pointers and citations in that file).
- Corollary R (leaf_proof/leaf_structure_proof.md s.13.2) is unconditional in every genus modulo (MV), textbook facts and, only to identify the graph sum with F_g, the inputs of Proposition 13.2; for the value of the pole also the further inputs of Theorem VI (LP2 draft s.7.8: LLZ eq. (5), vertex-operator form of psi, psi^*; one referee pass). Genus 3 and 4 are no longer computer-assisted. Graphs with an edge need Lemma A only for h <= g - 1.
- Round-5 critic item S3 applied to the conifold notes (backups *_before_round6.md): s.13.0 of leaf_structure_proof.md ports Fpot_h on a domain and (H0) (S3(a)); "non-empty" in M_irr (S3(b)); Teleman citations checked against the arXiv:0712.0160v3 text, field embedding added (S3(c)); input (5), step (v) of Theorem R restated via Theorem VI, [BS] Prop. 3.1 as a check (S3(d)).
- Not refereed yet: s.13.0, s.13.3 (5), and the field-of-definition paragraph of Proposition 13.2 (the part recovering d/dt through uniqueness and the divisor equation goes beyond the critic's item).
- Next: a referee of those three passages; proof_M_residue (round2_status s.5.5), README.txt, outline.md and the memory notes still carry the old status; hvc_lib.m line ~181 still calls the closed form a prediction; Teleman's journal numbering not compared with arXiv v3.
- Local-P2 status and open items (Prop 13.2_P2, gap corollary): hodge_vertex/lp2_explore/ROUND5_STATUS_2026-11-02.md and later round-6 files.

""")
s=rep(s,"""   Unverified points: Teleman was read through tool summaries; the CohFT axioms for non-compact Z go through the twisted P^1 theory.
""","""   Unverified points: Teleman was read through tool summaries [2026-11-02: arXiv v3 text read, citations corrected]; the CohFT axioms for non-compact Z go through the twisted P^1 theory.
""")
s=rep(s,"""   - M_irr only below the genus: leaf_structure_proof.md s.0.4, s.13.1, s.13.2 and G1; proof_status_mixed_scale.md s.3 and Addendum 2;
     conifold_regularity_note.md s.4 and s.7b.
   - Cite Brini-Schuler Prop. 3.1 for (H) at s.13.3(5).
""","""   - DONE 2026-11-02: M_irr only below the genus: leaf_structure_proof.md s.0.4, s.13.1, s.13.2 and G1; proof_status_mixed_scale.md s.3 and Addendum 2 (now Addendum 3);
     conifold_regularity_note.md s.4 and s.7b.
   - DONE 2026-11-02, differently: s.13.3(5) now rests on Theorem VI; Brini-Schuler Prop. 3.1 is cited as an independent check.
""")
save(p,s)
print("ok")
