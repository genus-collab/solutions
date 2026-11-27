# Moved on 2026-11-27 from localp2/hodge_vertex/lp2_explore/round6_scripts/round6d/fix/apply_conifold.py (local P2 project) into the resolved-conifold
# project, because it edits only files of this project. One-off patch script, already applied: a record, not to be rerun.
# Exact-match edits of the conifold files, round 6d minor items. All-or-nothing across files.
import os
import io, sys
R = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "..", "hodge_vertex_conifold", "")
files = {k: R + v for k, v in dict(leaf="leaf_proof/leaf_structure_proof.md", ps="proof_status_mixed_scale.md",
         note="conifold_regularity_note.md", res="lemmaA/RESUME.md", lem="lemmaA/LEMMA_A_PROOF.md").items()}
txt = {k: io.open(p, encoding="utf-8", newline="").read() for k, p in files.items()}
def rep(k, old, new, n=1):
    c = txt[k].count(old)
    if c != n:
        print("COUNT", c, "expected", n, k, ":", old[:120]); sys.exit(1)
    txt[k] = txt[k].replace(old, new)

MARK_OLD = "[repaired in round 6d; the repaired paragraph has had the passes recorded in the round 6 status file]"
MARK_NEW = "[repaired in round 6d; passes on the repaired paragraph: hodge_vertex/lp2_explore/ROUND6_STATUS_2026-11-02.md]"
rep("leaf", MARK_OLD, MARK_NEW, 8)
rep("ps", MARK_OLD, MARK_NEW, 1)
rep("res", MARK_OLD, MARK_NEW, 2)
rep("note", MARK_OLD, MARK_NEW, 1)

# leaf: status line
rep("leaf", "Mumford, Faber-Pandharipande Thm 4 and Stirling; Teleman Prop 7.7",
    "Mumford, Faber-Pandharipande Thm 4, Stirling and the further citations of (G3); Teleman Prop 7.7")
# leaf: S^i
rep("leaf", "(S^i the i-th solution, in the flat basis {1, H})",
    "(S^i the i-th column of S, paired by eta with the flat basis {1, H}: v_i(x) = Delta_i^(1/2) e^{-u_i/z} eta(S^i, x))")
# leaf: F6
rep("leaf", "statements about moduli of curves, monomials in the kappa_j (j >= 1), psi_1, psi_2 are independent in H^*(M^o_{g,2}; Q) in total degree < (g - 1)/2, and over F_0 by base change.",
    "statements about moduli of curves, monomials in psi_1, psi_2 and the classes kappa_j^M pulled back from M^o_g (j >= 1) are independent in H^*(M^o_{g,2}; Q) in total cohomological degree < (g - 1)/2 (deg psi_l = 2, deg kappa_j = 2j). The kappa_j here are the Arbarello-Cornalba classes, and on M^o_{g,2} they satisfy kappa_j = kappa_j^M + psi_1^j + psi_2^j ([T12] (2.15)), a triangular degree-preserving change of generators. So monomials in the kappa_j (j >= 1), psi_1, psi_2 are independent in the same range, and over F_0 by base change. Only independence in some range growing with g is used.")
rep("leaf", "Since M does not depend on g, M_n = 0 for n >= 1 (take g > 2n + 1).",
    "Since M does not depend on g, M_n = 0 for n >= 1 (take g > 4n + 1: psi_1^n has cohomological degree 2n, and the range of [T12] 2.18 is in cohomological degree).")
# leaf: target facts
rep("leaf", "The facts about the target used in rank 2 are: L_+ - L_- is a unit (for the e_i, Psi and the ambiguity statement); Delta_i(0) != 0 (for f_i and the invertibility of P); (F4); and (F8).",
    "Beyond the semisimple Frobenius data (L_+ - L_- a unit, for the e_i, Psi and the ambiguity statement; Delta_i(0) != 0, for f_i and the invertibility of P), the target enters only through (F4) and (F8).")
# leaf: evidence
rep("leaf", "(rank2_check.py, round 6d repair folder)", "(hodge_vertex/lp2_explore/round6_scripts/round6d/repair/rank2_check.py)")
rep("leaf", "It does not test (F1)-(F7).",
    "It does not test (F1)-(F7). Evidence, not a step: g1_check.py and g0_check.py (hodge_vertex/lp2_explore/round6_scripts/round6d/ref_con/) compare the rows with the genus-0 potential and genus-1 localisation; the first separates R(z) from R(-z), the second R from R^*, at the stated orders only. Round 6d evidence, not a step (SP/lp2_round6/round6d/classlevel/part3.py, part2b.py, SP as in the Appendix of the LP2 draft; own Wick expansion (C7) at generic weights): the defect of (F5) built from the rows of (b) vanishes exactly over Q at three samples with both e(wt) squares (through z^5, Q^4, and z^7, Q^6 at one of them) and in F_p with the square roots as elements of F_p at four samples, one through z^9, Q^6. Controls R^-1(-z), R(-z), H * -> -H * fail; dropping h_nw gives M = 0, as the ambiguity statement says. GP localisation on P^1 with psi classes pulled back from M_{g,n} equals the graph sum with legs for (g,n) = (1,1), (1,2), (2,0), (2,1), (2,2) through d = 3 at two weight points; legs at +psi, R^-1(-z), R(-z) fail. This observes M = 0 at finite order; it does not test the stability argument of (F6).")
# leaf: (c)
rep("leaf", "pi^*: A^*(M_k) -> A^*(M_{k,1})", "pi^*: H^*(M_k; Q) -> H^*(M_{k,1}; Q)")
# leaf: G3
rep("leaf", "Graber-Pandharipande (virtual localisation);",
    "Graber-Pandharipande (virtual localisation); the fundamental-class and splitting axioms of the virtual class of P^1 (Behrend-Manin, not opened);")

# LEMMA_A
rep("lem", "P modulo literature as stated there, assembled text refereed only as recorded there",
    "P modulo the literature and the draft statements listed there, conditional on OPEN item O1 of FRAMEB_PROOF.md (source of eq:J=I at T_max), assembled text refereed only as recorded there")
rep("lem", "are now in FRAMEB_PROOF.md (P modulo literature)",
    "are now in FRAMEB_PROOF.md (P modulo the literature and the draft statements listed there, conditional on its OPEN item O1; assembled text refereed only as recorded there)")

# RESUME
rep("res", "hodge_vertex/lp2_explore/ROUND5_STATUS_2026-11-02.md and later round-6 files.",
    "hodge_vertex/lp2_explore/ROUND5_STATUS_2026-11-02.md and later round-6 files. Round 6d: Prop 13.2_P2 and the gap corollary are conditional on OPEN item O1 of FRAMEB_PROOF.md (citation level: the source of eq:J=I at generic alpha on T_max).")

for k, p in files.items():
    io.open(p, "w", encoding="utf-8", newline="").write(txt[k])
print("conifold edits applied")
