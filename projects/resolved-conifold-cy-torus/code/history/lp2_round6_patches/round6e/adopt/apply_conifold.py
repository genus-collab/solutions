# Moved on 2026-11-27 from localp2/hodge_vertex/lp2_explore/round6_scripts/round6e/adopt/apply_conifold.py (local P2 project) into the resolved-conifold
# project, because it edits only files of this project. One-off patch script, already applied: a record, not to be rerun.
# Round 6e: exact-match edits of LEMMA_A_PROOF.md (lines 427, 552 only), RESUME.md, leaf_structure_proof.md.
import os
import sys
R = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "..", "hodge_vertex_conifold", "")
fails = []


def load(p):
    return open(p, encoding='utf-8', newline='').read()


def rep(s, old, new, count=1):
    c = s.count(old)
    if c != count:
        fails.append((c, count, old[:120]))
        return s
    return s.replace(old, new)


# LEMMA_A: only the two lines; line 545 (that file's own O1) is left alone.
pL = R + 'lemmaA/LEMMA_A_PROOF.md'
L = load(pL)
lines = L.split('\n')
assert lines[426].startswith('- The frame-B formula identifying the graph sum')
assert lines[551].startswith('  - lem:modrel (critic 1 S4)')
lines[426] = rep(lines[426],
                 "P modulo the literature and the draft statements listed there, conditional on OPEN item O1 of FRAMEB_PROOF.md (source of eq:J=I at T_max), assembled text refereed only as recorded there)",
                 "P modulo the literature and the draft statements listed there, among them eq:J=I (CCIT15 Thm 31), subject to the Q1 proviso of item O1 of FRAMEB_PROOF.md (the T_max use awaits the user's confirmation), assembled text refereed only as recorded there)")
lines[551] = rep(lines[551],
                 "P modulo the literature and the draft statements listed there, conditional on its OPEN item O1; assembled text refereed only as recorded there)",
                 "P modulo the literature and the draft statements listed there, among them eq:J=I (CCIT15 Thm 31), subject to the Q1 proviso of the item O1 of that file, which is not the O1 of s.10.2 here; assembled text refereed only as recorded there)")
L = '\n'.join(lines)

# RESUME
pRes = R + 'lemmaA/RESUME.md'
Res = load(pRes)
Res = rep(Res,
          "- Next: the passes on the repaired paragraph of Proposition 13.2; proof_M_residue (round2_status s.5.5),",
          "- The repaired paragraph of Proposition 13.2 has had its round-6d referee pass (refCon2; hodge_vertex/lp2_explore/ROUND6_STATUS_2026-11-02.md). AI pass only; no human referee has read it.\n- Next: proof_M_residue (round2_status s.5.5),")
Res = rep(Res,
          "- Local-P2 status and open items (Prop 13.2_P2, gap corollary): hodge_vertex/lp2_explore/ROUND5_STATUS_2026-11-02.md and later round-6 files. Round 6d: Prop 13.2_P2 and the gap corollary are conditional on OPEN item O1 of FRAMEB_PROOF.md (citation level: the source of eq:J=I at generic alpha on T_max).",
          "- Local-P2 status and open items (Prop 13.2_P2, gap corollary): hodge_vertex/lp2_explore/ROUND6_STATUS_2026-11-02.md (full round-6 status; ROUND5_STATUS_2026-11-02.md for the rounds before). Round 6e: eq:J=I at generic alpha on T_max is now cited to CCIT15 Thm 31 (FRAMEB_PROOF.md s.1.2, item O1), so Prop 13.2_P2 and the gap corollary are P mod eq:J=I (CCIT15 Thm 31) and the other inputs listed there. The T_max use is outside the literal wording of ruling Q1, which names statements on T_F, and awaits the user's confirmation.")

# leaf: (F1), (F3) to match A2.8; evidence pointer
pLeaf = R + 'leaf_proof/leaf_structure_proof.md'
Leaf = load(pLeaf)
Leaf = rep(Leaf,
           "so the translation legs push forward to a polynomial in the kappa_j, j >= 1, with no boundary correction:",
           "so the translation legs push forward to a series in the kappa_j, j >= 1, finite in each degree, with no boundary correction:")
Leaf = rep(Leaf,
           "      with G_i a polynomial in the kappa_j, j >= 1.",
           "      with G_i a series in the kappa_j, j >= 1, finite in each degree.")
Leaf = rep(Leaf,
           "(SP/lp2_round6/round6d/classlevel/part3.py, part2b.py, SP as in the Appendix of the LP2 draft; own Wick expansion (C7) at generic weights)",
           "(hodge_vertex/lp2_explore/round6_scripts/round6d/classlevel/part3.py, part2b.py; own Wick expansion (C7) at generic weights)")

if fails:
    for f in fails:
        print('FAIL', f)
    sys.exit(1)
open(pL, 'w', encoding='utf-8', newline='').write(L)
open(pRes, 'w', encoding='utf-8', newline='').write(Res)
open(pLeaf, 'w', encoding='utf-8', newline='').write(Leaf)
print('ok')
