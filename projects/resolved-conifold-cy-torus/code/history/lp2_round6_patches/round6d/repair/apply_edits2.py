# Moved on 2026-11-27 from localp2/hodge_vertex/lp2_explore/round6_scripts/round6d/repair/apply_edits2.py (local P2 project) into the resolved-conifold
# project, because it edits only files of this project. One-off patch script, already applied: a record, not to be rerun.
import os
LEAF = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "..", "hodge_vertex_conifold", "leaf_proof", "leaf_structure_proof.md")
s = open(LEAF, encoding="utf-8", newline="").read()
eds = [
 ("f_i^*(H) f_i^*(((B(psi_l) - 1)/psi_l) x_l). The identities", "f_i^*(H) f_i^*(((B(psi_l) - 1)/psi_l) x_l), times N_i and the other legs. The identities"),
 ("the i-th column of z S' = H * S is z D_i v_i", "the i-th column of z S' = H * S, lowered by eta, is z D_i v_i"),
 ("Rank 2 therefore needs exactly three facts about the target: L_+ - L_- is a unit (for the e_i, Psi, and the ambiguity statement), (F4), and (F8).",
  "The facts about the target used in rank 2 are: L_+ - L_- is a unit (for the e_i, Psi and the ambiguity statement); Delta_i(0) != 0 (for f_i and the invertibility of D); (F4); and (F8)."),
]
for a, b in eds:
    assert s.count(a) == 1, a
    s = s.replace(a, b)
open(LEAF, "w", encoding="utf-8", newline="").write(s)
print("ok")
