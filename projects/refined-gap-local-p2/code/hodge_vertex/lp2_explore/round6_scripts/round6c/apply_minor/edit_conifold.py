# 2026-11-27: this one-off patch script (already applied; a record, not to be rerun) edits files of both the local P2 project and the
# resolved-conifold project. To keep the two projects standalone, the resolved-conifold root is no longer reached by a relative path:
# it is read from the environment variable GENUS_RESOLVEDCONIFOLD, and the script stops if that is not set.
import os as _os_genus, sys as _sys_genus
def _genus_rc_root():
    r = _os_genus.environ.get("GENUS_RESOLVEDCONIFOLD")
    if not r: _sys_genus.exit("set GENUS_RESOLVEDCONIFOLD to the resolved-conifold project folder (this script edits files there)")
    return r
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib import Doc

OPEN1 = "OPEN (round 6c referee: sign u = s H in the field-of-definition paragraph of (a); by [T12] Def. 7.1 the family is at u = -s H, and as written the step reaches the QDE with z -> -z)"

# ------------------------------------------------------------------ LSP
d = Doc(os.path.join(_genus_rc_root(), "hodge_vertex_conifold", "leaf_proof", "leaf_structure_proof.md"))
r = d.rep

r("""Revised 2026-11-02 (round 6: s.0.4, s.12.1, s.13, s.14 only; Theorem LC and sections 1-11 are unchanged).""",
  """Revised 2026-11-02 (round 6: s.0.4, s.12.1, s.13, s.14 only; Theorem LC and sections 1-11 are unchanged; round 6c: minor referee items in s.13, s.14 and one OPEN mark in Proposition 13.2).""")

r("""Corollary R is unconditional modulo (MV), textbook facts and, only to identify the graph sum with F_g, the inputs of Proposition 13.2 (exact
list and refereeing status: s.13.2).""",
  """Corollary R, a statement about the decoration sum, is proved modulo (MV) and textbook facts. Its identification with F_g is Proposition 13.2, which carries an OPEN mark from round 6c (exact
list and refereeing status: s.13.2).""")

r(""" - Proposition 13.2 (the frame-B formula, input (1) of Theorem R) [cited: Teleman,""",
  """ - Proposition 13.2 (the frame-B formula, input (1) of Theorem R) [""" + OPEN1 + """; otherwise cited: Teleman,""")

r("""The empty multiset is excluded: Xi_h(empty) = Z[h, 0] = Fpot_h(q_0) is the gap, non-zero in Q(x)""",
  """The empty multiset is excluded: (-x c^3)^(h-1) Xi_h(empty) = (-x c^3)^(h-1) Fpot_h(q_0) is the gap, non-zero in Q(x)""")

r("""Every line was re-derived for the port; nothing here depends on the rank. Not yet refereed in this form.]""",
  """Every line was re-derived for the port; nothing here depends on the rank. One referee pass (round 6c).]""")

r("""Hence
G(s) = (1 + s mu)^(2-2h-n) G(0); put s = 1. QED""",
  """The solution with value G(0) at s = 0 is unique in Rr[[s]] (compare s-coefficients; characteristic 0), so
G(s) = (1 + s mu)^(2-2h-n) G(0); put s = 1. QED""")

# Corollary R statement
r("""single-vertex graph has the single pole gap_g t^(2-2g), and F4 at n = 0 holds in genus g.
Status (2026-11-02). P modulo: (MV) = LLZ eq. (1) (for M_irr, and for Theorem VI); the textbook facts of lemmaA/LEMMA_A_PROOF.md s.0; for the value of the pole only, the further inputs of Theorem VI (draft s.7.8:""",
  """polar part of the single-vertex graph is gap_g t^(2-2g) (non-zero in Q(x) and at real x; it vanishes at the finitely many complex x with beta_{2g-2}(x, 1) = 0), and F4 at n = 0 holds in genus g.
Status (2026-11-02, round 6c). [P modulo: (MV) = LLZ eq. (1), for M_irr(h) with h >= 3, and, only at mu = 1^a, for the pole value; the textbook facts of lemmaA/LEMMA_A_PROOF.md s.0, Frobenius' transposition identity, the principal specialisation; Stirling's series with remainder, through Theorem LC. It is a statement about the decoration sum (C6)-(C11) as defined. No TR, no BS input, no finite-order step. Its identification with the GW potential F_g is Proposition 13.2, which is not claimed as proved: """ + OPEN1 + """. Apart from that step, Proposition 13.2 rests on Teleman Theorem (2) and Prop 7.7 (proof read, not re-derived), Graber-Pandharipande, Coates-Givental App. 1, PPZ Defs 2.2, 2.13, Prop 2.14, Mumford, Faber-Pandharipande Thm 4 and Stirling. Refereeing: Lemma A two passes (round 6); Theorem VI two passes (rounds 5, 6); s.13.0, s.13.3 (5), and the field-of-definition paragraph one pass each (round 6c). Section 13 as a whole has had no pass of its own. "Unconditional" is not to be used without this list.] Further inputs of Theorem VI, for the value of the pole only (draft s.7.8:""")

r("""F1, F3 and Lemma C are Theorem LC. The inputs of Proposition 13.2""",
  """F1, F3 and Lemma C are Theorem LC, which cites Stirling's series with remainder (DLMF 5.11.1, s.5.11(ii); (G2)). The inputs of Proposition 13.2""")

r("""Theorem VI has had one referee pass (round 5); s.13.0, s.13.3 (5) and the field-of-definition paragraph of Proposition 13.2 have had none in this form. By genus:""",
  """Theorem VI has had two referee passes (round 5; round 6, D_VI); s.13.0, s.13.3 (5) and the field-of-definition paragraph of Proposition 13.2 have had one pass (round 6c), and the pass on that paragraph found the major item marked OPEN there. By genus:""")

# s.13.3 (5)
r("""a statement about Hodge integrals
only, with no TR, no Teleman and nothing from Lemma A;""",
  """a statement about Hodge integrals
only, with no TR, no Teleman and nothing from Theorem L (it quotes only lemmaA/LEMMA_A_PROOF.md Prop. 2.1, Lemma 3.2 and eq. (4.2), whose proofs do not involve Theorem L);""")

# Prop 13.2 heading
r("""**Proposition 13.2 (the frame-B formula).** Let T""",
  """**Proposition 13.2 (the frame-B formula).** [Status lowered in round 6c: not claimed as proved. """ + OPEN1 + """. Regularity and the pole value in Corollary R do not use this paragraph. The other steps rest on the citations of (G3).] Let T""")

# ambiguity
r("""with constant a_k ([T12] Proposition 8.4; [G01] s.1.3, Proposition
(c); its proof is formal and holds over Kw((Q)) with d/dt = Q d/dQ).""",
  """with constant a_k (FRAMEB_PROOF.md A2.9, with no citation: X := Psi^-1 R_1^-1 R_2 Psi satisfies z X' = [Lambda, X], so (L_i - L_j)(X_n)_ij = (X_{n-1}')_ij with L_+ - L_- a unit; hence X_n is diagonal with Q d/dQ-constant entries, the constants of Q d/dQ on Kw((Q)) are Kw, and symplecticity makes the exponent odd; compare [T12] Proposition 8.4, [G01] s.1.3, Proposition (c)).""")
r("""and R = R_gen Exp[diag Sum_k a_k z^(2k-1)]
by [G01] Proposition (c).""",
  """and R = R_gen Exp[diag Sum_k a_k z^(2k-1)]
by the ambiguity statement in (a).""")

# OPEN mark at the paragraph
r("""which comes back as follows [round 6, not refereed]: the datum""",
  """which comes back as follows [""" + OPEN1 + """; repair not applied]: the datum""")

# colour isolation
r("""for k = 1, Int_{M_{1,1}} lambda_1 = 1/24. Induction on k gives a^p = 0""",
  """for k = 1, Int_{M_{1,1}} lambda_1 = 1/24. Apply this on M_{k,1} with the insertion e_p, since on M_k, which has no marking, Omega_{k,0} is a sum over the colours and one colour cannot be separated. By the flat unit, (s_p.1)_{k,1} = pi^*(s_p.1)_{k,0}, and pi^*: A^*(M_k) -> A^*(M_{k,1}) is injective for k >= 2 because pi_*(psi_1 pi^* y) = (2k - 2) y, so ch_{2k-1}(E) != 0 on M_{k,1}. For k = 1, M_{1,1} is used directly. Induction on k gives a^p = 0""")

# (G3)
r(""" (G3) Closed as a write-up item. The premise that""",
  """ (G3) Written up, not closed: """ + OPEN1 + """. The premise that""")
r("""is now Proposition 13.2, proved from: the CohFT""",
  """is now Proposition 13.2, whose proof, apart from that step, rests on: the CohFT""")
r("""numbering of arXiv:0712.0160v3: (1.6), (1.15)
      Theorem (2), Definition 7.1, Example 7.17, Propositions 7.7, 8.2(ii), 8.4, (8.6), Proposition 8.7);""",
  """numbering of arXiv:0712.0160v3: (1.6), (1.15)
      Theorem (2), Definition 7.1, (7.6), Propositions 7.7, 8.2(ii), 8.7; Example 7.17 for the formal case (the sentence marked OPEN in Proposition 13.2(a) follows its printed sign, which is the defect named there); Proposition 8.4 and (8.6) for comparison; the proof of Proposition 7.7 was read, not re-derived);""")

d.save()

# ------------------------------------------------------------------ status note
p = Doc(os.path.join(_genus_rc_root(), "hodge_vertex_conifold", "proof_status_mixed_scale.md"))
r = p.rep
r("""So Corollary R (leaf_proof/leaf_structure_proof.md s.13.2) is unconditional modulo (MV), textbook facts and, only to identify the graph sum with F_g, the inputs of Proposition 13.2 there; genus 3 and 4""",
  """So Corollary R (leaf_proof/leaf_structure_proof.md s.13.2), a statement about the decoration sum, is proved modulo (MV) and textbook facts. Its identification with F_g is Proposition 13.2 there, which is not claimed as proved: """ + OPEN1.replace('[T12]', 'Teleman') + """; genus 3 and 4""")
r("""LLZ eq. (5) and the vertex-operator form of psi, psi^* are no longer inputs; one referee pass; specialisation""",
  """LLZ eq. (5) and the vertex-operator form of psi, psi^* are no longer inputs; two referee passes (round 5; round 6, D_VI); specialisation""")
p.save()

# ------------------------------------------------------------------ draft and FRAMEB small fixes
q = Doc(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "LP2_REGULARITY_PROOF_draft.md"))
q.rep("""with the printed Asym) , item 9 (lines""", """with the printed Asym), item 9 (lines""")
q.rep("""P mod C     Stages A, B; Prop 2.10; Mumford through A3.2; PPZ Prop 2.14;""",
      """P mod C     Stages A, B; Prop 2.10; Mumford through A3.2; PPZ Defs 2.2, 2.5, 2.13, Prop 2.14;""")
q.save()
