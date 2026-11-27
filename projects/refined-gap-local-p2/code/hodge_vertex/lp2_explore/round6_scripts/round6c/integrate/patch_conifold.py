# 2026-11-27: this one-off patch script (already applied; a record, not to be rerun) edits files of both the local P2 project and the
# resolved-conifold project. To keep the two projects standalone, the resolved-conifold root is no longer reached by a relative path:
# it is read from the environment variable GENUS_RESOLVEDCONIFOLD, and the script stops if that is not set.
import os as _os_genus, sys as _sys_genus
def _genus_rc_root():
    r = _os_genus.environ.get("GENUS_RESOLVEDCONIFOLD")
    if not r: _sys_genus.exit("set GENUS_RESOLVEDCONIFOLD to the resolved-conifold project folder (this script edits files there)")
    return r
# ASCII edits only; aligns the inputs of the vertex identity with draft s.7.8. Aborts if a target is not unique.
import os
R=os.path.join(_genus_rc_root(), "hodge_vertex_conifold", "")
D=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "LP2_REGULARITY_PROOF_draft.md")
def ed(path,old,new):
    s=open(path,encoding='utf-8',newline='').read()
    assert s.count(old)==1,(path,old[:60],s.count(old))
    s=s.replace(old,new)
    open(path,'w',encoding='utf-8',newline='').write(s)
ed(R+"leaf_proof/leaf_structure_proof.md",
"for the value of the pole only, the further inputs of Theorem VI (LLZ eq. (5), the\ncut-and-join eigenvalue; the vertex-operator form of psi, psi^*). F1, F3",
"for the value of the pole only, the further inputs of Theorem VI (draft s.7.8: (MV) only at mu = 1^a; lemmaA/LEMMA_A_PROOF.md Prop. 2.1, Lemma 3.2 and eq. (4.2) at k = 2, the last derived there from Murnaghan-Nakayama and the Frobenius formula; Frobenius' transposition identity for the\ncut-and-join eigenvalue on the p_1-slice; the principal specialisation; draft (7.2). LLZ eq. (5) and the vertex-operator form of psi, psi^* are no longer inputs). F1, F3")
ed(R+"proof_status_mixed_scale.md",
"> further inputs LLZ eq. (5) and the vertex-operator form of psi, psi^*; one referee pass;",
"> further inputs, beyond (MV) at mu = 1^a and textbook facts: lemmaA/LEMMA_A_PROOF.md Prop. 2.1, Lemma 3.2 and eq. (4.2) at k = 2, with the Frobenius formula; LLZ eq. (5) and the vertex-operator form of psi, psi^* are no longer inputs; one referee pass;")
ed(R+"lemmaA/RESUME.md",
"(LP2 draft s.7.8: LLZ eq. (5), vertex-operator form of psi, psi^*; one referee pass)",
"(LP2 draft s.7.8: lemmaA/LEMMA_A_PROOF.md Prop. 2.1, Lemma 3.2 and eq. (4.2) at k = 2, with the Frobenius formula, and (MV) only at mu = 1^a; LLZ eq. (5) and the vertex-operator form of psi, psi^* are no longer inputs; one referee pass)")
ed(D,"the parenthesised labels (Cn) are never to be read as C1cf or C1i.\n",
"the parenthesised labels (Cn) are never to be read as C1cf or C1i, nor as the parts (C1), (C2) of Theorem FB_C.\n")
print('ok')
