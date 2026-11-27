# Moved on 2026-11-27 from localp2/hodge_vertex/lp2_explore/round6_scripts/round6d/repair/apply_edits3.py (local P2 project) into the resolved-conifold
# project, because it edits only files of this project. One-off patch script, already applied: a record, not to be rerun.
import os
LEAF = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "..", "hodge_vertex_conifold", "leaf_proof", "leaf_structure_proof.md")
s = open(LEAF, encoding="utf-8", newline="").read()
eds = [
 ("Sum_i D_i N_i(kappa) prod_l a^i(psi_l; x_l),   D_i := Delta_i^((2g-2+n)/2),   N_i(0) = 1.",
  "Sum_i P_i N_i(kappa) prod_l a^i(psi_l; x_l),   P_i := Delta_i^((2g-2+n)/2),   N_i(0) = 1"),
 ("Only omega(f_i^{(x)n}) = D_i and the vanishing", "(P_i is the D_i of A2.8; D_i below is the operator of (b)). Only omega(f_i^{(x)n}) = P_i and the vanishing"),
 ("= Sum_i D_i N_i Sum_l f_i^*(", "= Sum_i P_i N_i Sum_l f_i^*("),
 ("(F3) and (F4) at n = 2, where D_i = Delta_i^g, read", "(F3) and (F4) at n = 2, where P_i = Delta_i^g, read"),
 ("      Sum_i D_i N_i [(M(psi_1) a_1)^i a_2^i", "      Sum_i P_i N_i [(M(psi_1) a_1)^i a_2^i"),
 ("Gtilde_i := G_i - theta(D_i N_i)", "Gtilde_i := G_i - theta(P_i N_i)"),
 ("with x_1, x_2 arbitrary and D := diag(Delta_+^g, Delta_-^g),", "with x_1, x_2 arbitrary and P := diag(Delta_+^g, Delta_-^g),"),
 ("      D M(psi) + M_0^T D = Gtilde(0), diagonal, modulo that degree.", "      P M(psi) + M_0^T P = Gtilde(0), diagonal, modulo that degree."),
 ("(for f_i and the invertibility of D)", "(for f_i and the invertibility of P)"),
 ("is (Delta_+^g - Delta_-^g)(M_0)_{+-} = 0.", "is (Delta_+^g - Delta_-^g)(M_0)_{+-} = 0, for any g >= 2."),
 ("since Delta_i lies in Kw[[Q]] with Delta_i(0) = e(wt_i) != 0 by (13.2).", "since Delta_i lies in Kw[[Q]] with Delta_i(0) = e(wt_i) != 0 by (13.2) in (b)(i), which uses only the Frobenius data, not this step."),
 ("(F8) Non-degeneracy (the second target input). By (13.2), with L_+(0) = 0", "(F8) Non-degeneracy (the second target input). By (13.2), with sig := alpha0 + alpha1 + alpha2, L_+(0) = 0"),
]
for a, b in eds:
    assert s.count(a) == 1, a
    s = s.replace(a, b)
s.encode("ascii")
open(LEAF, "w", encoding="utf-8", newline="").write(s)
print("ok")
