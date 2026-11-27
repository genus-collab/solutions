# 2026-11-27: this one-off patch script (already applied; a record, not to be rerun) edits files of both the local P2 project and the
# resolved-conifold project. To keep the two projects standalone, the resolved-conifold root is no longer reached by a relative path:
# it is read from the environment variable GENUS_RESOLVEDCONIFOLD, and the script stops if that is not set.
import os as _os_genus, sys as _sys_genus
def _genus_rc_root():
    r = _os_genus.environ.get("GENUS_RESOLVEDCONIFOLD")
    if not r: _sys_genus.exit("set GENUS_RESOLVEDCONIFOLD to the resolved-conifold project folder (this script edits files there)")
    return r
# Exact-match edits of round 6d (repair). Every old string must occur exactly once.
import os
import sys

ROOT = os.path.join(_genus_rc_root(), "")
LEAF = ROOT + "hodge_vertex_conifold/leaf_proof/leaf_structure_proof.md"
STAT = ROOT + "hodge_vertex_conifold/proof_status_mixed_scale.md"
NOTE = ROOT + "hodge_vertex_conifold/conifold_regularity_note.md"
RESU = ROOT + "hodge_vertex_conifold/lemmaA/RESUME.md"
LEMA = ROOT + "hodge_vertex_conifold/lemmaA/LEMMA_A_PROOF.md"
FRAM = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "FRAMEB_PROOF.md")
DRAF = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "LP2_REGULARITY_PROOF_draft.md")

MARK = "[repaired in round 6d; the repaired paragraph has had the passes recorded in the round 6 status file]"
OPENT = "OPEN (round 6c referee: sign u = s H in the field-of-definition paragraph of (a); by [T12] Def. 7.1 the family is at u = -s H, and as written the step reaches the QDE with z -> -z)"


def load(p):
    with open(p, encoding="utf-8", newline="") as f:
        return f.read()


def save(p, s):
    with open(p, "w", encoding="utf-8", newline="") as f:
        f.write(s)


PARAGRAPH = r"""The embedding forgets d/dt, which comes back as follows """ + MARK + r""".
Descent. Freeness: if R.omega = omega over a field of characteristic 0 then R = 1 (FRAMEB_PROOF.md A1.7; [T12] (1.15), Proposition 8.7(ii) state the uniqueness). A1.7 uses no rank. For R = 1 + R_m z^m + ... and i != j, the degree-m part of (R.omega)_{0,m+3}(f_i, f_j, ..., f_j) is the single leg term -(R_m)_{ji} Delta_j^((m+1)/2) psi_1^m: the other legs and the translations give mixed insertions, and a graph with an edge would need a vertex carrying f_i and the edge alone, which is unstable. So R_m is diagonal, symplecticity makes m = 2k - 1 odd, and its entries vanish by Mumford's formula and ch_{2k-1}(E) != 0 on M_{k,1}, as in (c) below. Read R in End(V (x) C), where V has the basis {1, H} over Kw((Q)). The e_i = (H - L_j)/(L_i - L_j) lie in V, because L_+ - L_- is a unit of Kw[[Q]] (Hensel; its value at Q = 0 is alpha0), so omega and Omega^T are defined over Kw((Q)) and no square root of Delta_i is needed for this step. Let phi be in Aut(C/Kw((Q))). The action is given by polynomial formulas over Q, so phi(R).omega = phi(R.omega) = Omega^T = R.omega; hence R^-1 phi(R) fixes omega ([PPZ] Prop. 2.14) and phi(R) = R. The fixed field of Aut(C/F) is F for every subfield F of C (an algebraic element outside F has a different conjugate, a transcendental y can be moved to y + 1, and these extend to C). So R has coefficients in Kw((Q)), and by freeness over Kw((Q)) it does not depend on the embedding.
Flatness. [T12] Proposition 7.7 is not used. The proof is FRAMEB_PROOF.md A2.8, from "Set-up" to "The equation", in rank 2 with the following data. F_0 := K~w((Q)), K~w := Kw(e(wt_+)^(1/2), e(wt_-)^(1/2)), a finite extension of Kw((Q)) containing Delta_i^(1/2), since Delta_i lies in Kw[[Q]] with Delta_i(0) = e(wt_i) != 0 by (13.2). ' := theta := Q d/dQ extends coefficientwise to F_0, with constants K~w. On V_0 := Kw 1 + Kw H, theta vanishes and eta has values in Kw. f_i := Delta_i^(1/2) e_i, and Psi is the matrix of (f_+, f_-) in the basis {1, H}, so Psi^T eta Psi = 1 and H * Psi = Psi Lambda with Lambda := diag(L_+, L_-). B(psi) := R^-1(psi), a^i(psi; x) := f_i^*(B(psi) x). M^o_{g,n} is the moduli space of smooth curves (M_{g,n} in this note is the compactification). Each line of A2.8 was re-derived by hand for these data. The reason why each line does not see the rank is as follows.
 (F1) Edgeless form. On M^o_{g,n} graphs with an edge vanish, being pushed forward from the boundary, and T(z) = z (1 - R^-1(z) 1) = O(z^2), so the translation legs push forward to a polynomial in the kappa_j, j >= 1, with no boundary correction:
      Omega^T_{g,n}(x_1, ..., x_n) = Sum_i D_i N_i(kappa) prod_l a^i(psi_l; x_l),   D_i := Delta_i^((2g-2+n)/2),   N_i(0) = 1.
      Only omega(f_i^{(x)n}) = D_i and the vanishing of mixed insertions are used.
 (F2) Pushforward. Let p: C -> M^o_{g,n} be the restriction of M_{g,n+1} -> M_{g,n}, with sections delta_l. On C: psi_l = p^*psi_l + delta_l, psi_l delta_l = psi_{n+1} delta_l = 0, delta_l delta_l' = 0 for l != l', kappa_j = p^*kappa_j + psi_{n+1}^j, p_* psi_{n+1}^k = kappa_{k-1} (kappa_0 = 2g - 2 + n is a number). These are statements about curves, with no target.
 (F3) Bubble. Fibres over M^o_{g,n} are smooth, so the only graphs with an edge that meet C are the bubbles (0,3) carrying the legs l, n+1. A translation leg on a bubble vanishes (degree >= 2m on M_{0,3+m}, of dimension m). The bubble vertex is the product and eta^-1, and the edge factor at psi_h = 0 is (eta^-1 - (B(psi_h') (x) 1) eta^-1)/psi_h', so the bubble puts ((1 - B(psi_l))/psi_l)(H * x_l) on leg l of the genus-g vertex. The edgeless graph contributes on delta_l the term Delta_i^((2g-1+n)/2) f_i^*(H) f_i^*(((B(psi_l) - 1)/psi_l) x_l). The identities H = Sum_i L_i e_i, i.e. f_i^*(H) = L_i Delta_i^(-1/2), and f_i^*(a * b) = Delta_i^(1/2) a^i b^i hold in any semisimple Frobenius algebra, and the two terms add up to (1/psi_l)[H *, B(psi_l)] x_l on leg l. Hence
      p_* Omega^T_{g,n+1}(x, H) = Sum_i D_i N_i Sum_l f_i^*((1/psi_l)[H *, B(psi_l)] x_l) prod_{l' != l} a^i(psi_l'; x_l') + Sum_i G_i(kappa) prod_l a^i(psi_l; x_l),
      with G_i a polynomial in the kappa_j, j >= 1.
 (F4) Divisor equation (A1.5 for P^1; the first of the two target inputs). For x in V_0: ft_*(ev_{n+1}^* H) = d, the fibre degree of f^*O(1), for the equivariant lift of H used here as for any other; the virtual class (fundamental-class axiom) and the twist ([CG07] App. 1(iv)) are pulled back. So p_* Omega^T_{g,n+1,d}(x, H) = d Omega^T_{g,n,d}(x), also at d = 0, and the left side of (F3) is theta Omega^T_{g,n}(x).
 (F5) Defect. Define M(psi) in Mat_2(F_0[[psi]]) by theta a^i(psi; x) - f_i^*((1/psi)[H *, B(psi)] x) = Sum_k M(psi)^i_k a^k(psi; x), which is possible because x -> a(psi; x) is invertible. M is built from R, Psi and H * alone, so it does not depend on g. (F3) and (F4) at n = 2, where D_i = Delta_i^g, read
      Sum_i D_i N_i [(M(psi_1) a_1)^i a_2^i + a_1^i (M(psi_2) a_2)^i] = Sum_i Gtilde_i(kappa) a_1^i a_2^i,   Gtilde_i := G_i - theta(D_i N_i),   a_l := a(psi_l; x_l).
 (F6) Stability. By [T12] 2.18 and 2.19 (Harer stability with Looijenga; Madsen-Weiss), statements about moduli of curves, monomials in the kappa_j (j >= 1), psi_1, psi_2 are independent in H^*(M^o_{g,2}; Q) in total degree < (g - 1)/2, and over F_0 by base change. The coefficient of kappa^0 psi_2^0 gives, with x_1, x_2 arbitrary and D := diag(Delta_+^g, Delta_-^g),
      D M(psi) + M_0^T D = Gtilde(0), diagonal, modulo that degree.
      Since M does not depend on g, M_n = 0 for n >= 1 (take g > 2n + 1).
 (F7) Skewness of M_0. M_0 = -W_Psi + [Lambda, R^f_1], with W_Psi := Psi^-1 theta Psi and R^f := Psi^-1 R Psi. W_Psi is skew because Psi^T eta Psi = 1 and theta eta = 0. R^f_1 is symmetric by symplecticity, so its commutator with the diagonal Lambda is skew. No rank enters. The off-diagonal part of (F6) at psi^0 is (Delta_+^g - Delta_-^g)(M_0)_{+-} = 0.
 (F8) Non-degeneracy (the second target input). By (13.2), with L_+(0) = 0 and L_-(0) = -alpha0: Delta_+(0) = alpha0 alpha1 alpha2 eps1 eps2, and Delta_-(0) = eps1 eps2 (alpha0 alpha1 alpha2 - 2 alpha0 alpha1 alpha2 - sig alpha0^2) = -alpha0 (alpha0 + alpha1)(alpha0 + alpha2) eps1 eps2. So Delta_+(0)/Delta_-(0) = -alpha1 alpha2/((alpha0 + alpha1)(alpha0 + alpha2)). This element of Kw is non-constant, and Q is algebraically closed in the purely transcendental extension Kw, so it is transcendental over Q and is not a root of unity. Hence Delta_+^g - Delta_-^g has a non-zero value at Q = 0, for every g >= 1, and M_0 = 0. The weights are indeterminate here. On T_F, alpha0 = 0 and Q = 0 is not a semisimple point; the passage to T_F is (d), not this step.
 (F9) The equation. So M = 0, that is, theta(Psi^-1 B(psi)) = Psi^-1 (1/psi)[H *, B(psi)]. Taking eta-adjoints, with B(psi)^dagger = R(-psi), Psi^dagger = Psi^-1 and theta eta = 0, and putting z = -psi gives
      z (R(z) Psi)' = (H * R(z) - R(z) H *) Psi.
      With H * Psi = Psi Lambda and U' = Lambda = diag(L_i), this says that S := R Psi e^{U/z} satisfies z S' = H * S, the QDE. The R here is the one of the legs R^-1(psi_l) of [PPZ] Def. 2.2, so no comparison of conventions is needed. (With R written in the frame f_i, as in the display above, S = Psi (Psi^-1 R Psi) e^{U/z}.) In terms of the rows of (b): by symplecticity and e_i^* = Delta_i eta(e_i, .), v_i(x) = e_i^*(R^-1(-z) x) = Delta_i^(1/2) eta(R(z) f_i, x), so with D_i := d/dt - Delta_i'/(2 Delta_i) the i-th column of z S' = H * S is z D_i v_i = (A - L_i) v_i, the equation of (b)(ii), with the same z.
Rank 2 therefore needs exactly three facts about the target: L_+ - L_- is a unit (for the e_i, Psi, and the ambiguity statement), (F4), and (F8). Evidence, not a step: exact arithmetic at three rational weight points through Q^7, z^3 (rank2_check.py, round 6d repair folder) confirms (13.2), the ratio in (F8), that the solution of the row equation with constants h_{wt_i} equals the Laplace row h_nw Lnorm[W] of (b) and fails at -z, and v_i = e_i^* o R^-1(-z). It does not test (F1)-(F7).
Remark (not a proof step): the route through [T12] Proposition 7.7. Let sigma := Exp[s Q d/dQ]: Kw((Q)) -> Kw((Q))[[s]], Q -> Q e^s. By (F4), Sum_m (s^m/m!) p^m_* Omega^T_{g,n+m}(x, H^m) = sigma(Omega^T_{g,n}(x)). [T12] Definition 7.1 carries the factor (-1)^m, so, applying the embedding coefficientwise in s, sigma(Omega^T) is Teleman's family at u = -s H, not at u = s H. (Against the true divisor equation in (F4), the divisor equation printed in [T12] Example 7.17 has the opposite sign; that printed sign is not relied on.) Applying sigma to coefficients gives sigma(R).sigma(omega) = sigma(Omega^T), so E_s := sigma(R) is the datum of this family, unique by freeness. (7.7.a) with v = -H gives at s = 0 z (R Psi)' = (H * R - R H *) Psi, the equation of (F9), and [T12] Proposition 8.2(ii) along v = -H is z S' = H * S. The legs of (7.6) carry E^-1(psi_i) and those of [PPZ] Def. 2.2 carry R^-1(psi_i), so in this comparison Teleman's E is the R of [PPZ], with the same z; this is a reading of (7.6) against Def. 2.2. The proof of Proposition 7.7 was read, not re-derived. No step of the proof above uses Proposition 7.7, Proposition 8.2(ii) or Example 7.17."""

EDITS = {}

EDITS[LEAF] = [
    ("round 6c: minor referee items in s.13, s.14 and one OPEN mark in Proposition 13.2).",
     "round 6c: minor referee items in s.13, s.14; round 6d: the field-of-definition paragraph of Proposition 13.2 rewritten " + MARK + ")."),
    ("Its identification with F_g is Proposition 13.2, which carries an OPEN mark from round 6c (exact",
     "Its identification with F_g is Proposition 13.2 " + MARK + " (exact"),
    ("[" + OPENT + "; otherwise cited: Teleman, Givental,",
     MARK + " [cited: Teleman, Givental,"),
    ("Status (2026-11-02, round 6c). [P modulo:", "Status (2026-11-02, round 6d). [P modulo:"),
    ("Its identification with the GW potential F_g is Proposition 13.2, which is not claimed as proved: " + OPENT + ". Apart from that step, Proposition 13.2 rests on Teleman Theorem (2) and Prop 7.7 (proof read, not re-derived), Graber-Pandharipande, Coates-Givental App. 1, PPZ Defs 2.2, 2.13, Prop 2.14, Mumford, Faber-Pandharipande Thm 4 and Stirling. Refereeing: Lemma A two passes (round 6); Theorem VI two passes (rounds 5, 6); s.13.0, s.13.3 (5), and the field-of-definition paragraph one pass each (round 6c).",
     "Its identification with the GW potential F_g is Proposition 13.2 " + MARK + ", P modulo Teleman Theorem (2), Prop 8.7(ii) and 2.18, 2.19 (Harer stability with Looijenga; Madsen-Weiss), Graber-Pandharipande, Coates-Givental App. 1, PPZ Defs 2.2, 2.13, Prop 2.14, Mumford, Faber-Pandharipande Thm 4 and Stirling; Teleman Prop 7.7 appears only in a remark, as no proof step. Refereeing: Lemma A two passes (round 6); Theorem VI two passes (rounds 5, 6); s.13.0, s.13.3 (5) one pass each (round 6c); the field-of-definition paragraph one pass in its round-6c form, which found the sign error now removed, and its repaired form as in the mark above."),
    ("s.13.0, s.13.3 (5) and the field-of-definition paragraph of Proposition 13.2 have had one pass (round 6c), and the pass on that paragraph found the major item marked OPEN there. By genus:",
     "s.13.0 and s.13.3 (5) have had one pass (round 6c); the field-of-definition paragraph of Proposition 13.2 had one pass (round 6c) in its earlier form, which found the sign error u = s H for u = -s H, " + MARK + ". By genus:"),
    ("**Proposition 13.2 (the frame-B formula).** [Status lowered in round 6c: not claimed as proved. " + OPENT + ". Regularity and the pole value in Corollary R do not use this paragraph. The other steps rest on the citations of (G3).]",
     "**Proposition 13.2 (the frame-B formula).** [P modulo Teleman Theorem (2), Prop 8.7(ii) and 2.18, 2.19 (Harer stability with Looijenga; Madsen-Weiss), Graber-Pandharipande, Coates-Givental App. 1, PPZ Defs 2.2, 2.13, Prop 2.14, Mumford, Faber-Pandharipande Thm 4, Stirling (through Lemma 13.3), and the further citations of (G3). Teleman Prop 7.7 appears only in a remark, as no proof step. Field-of-definition paragraph of (a): " + MARK + ". Regularity and the pole value in Corollary R do not use this paragraph.]"),
    ("QDE ([T12] Propositions 7.7, 8.2(ii); [G01] s.1.3, (3)-(5)); two such R",
     "QDE (proved under \"Flatness\" at the end of this step, with no use of [T12] Proposition 7.7; [G01] s.1.3, (3)-(5) for the form of the system); two such R"),
    (" (G3) Written up, not closed: " + OPENT + ". The premise",
     " (G3) Written up " + MARK + ". The premise"),
    ("whose proof, apart from that step, rests on: the CohFT of Z at indeterminate weights (virtual localisation, as a twisted theory",
     "whose proof rests on: the flatness of R by the direct derivation of Proposition 13.2(a), (F1)-(F9); the CohFT of Z at indeterminate weights (virtual localisation, as a twisted theory"),
    ("Theorem (2), Definition 7.1, (7.6), Propositions 7.7, 8.2(ii), 8.7; Example 7.17 for the formal case (the sentence marked OPEN in Proposition 13.2(a) follows its printed sign, which is the defect named there); Proposition 8.4 and (8.6) for comparison; the proof of Proposition 7.7 was read, not re-derived);",
     "Theorem (2), 2.18, 2.19 (Harer stability with Looijenga; Madsen-Weiss; originals not opened), Proposition 8.7; Definition 7.1, (7.6), Propositions 7.7, 8.2(ii) and Example 7.17 only in the remark closing Proposition 13.2(a), as no proof step, the printed sign of Example 7.17 not being relied on; Proposition 8.4 and (8.6) for comparison; the proof of Proposition 7.7 was read, not re-derived);"),
]

EDITS[STAT] = [
    ("Its identification with F_g is Proposition 13.2 there, which is not claimed as proved: OPEN (round 6c referee: sign u = s H in the field-of-definition paragraph of (a); by Teleman Def. 7.1 the family is at u = -s H, and as written the step reaches the QDE with z -> -z); genus 3 and 4",
     "Its identification with F_g is Proposition 13.2 there " + MARK + ", P modulo the literature inputs listed in its heading and in (G3); genus 3 and 4"),
]

EDITS[DRAF] = [
    ("its field-of-definition paragraph carries an OPEN mark from round 6c (LSP Prop 13.2(a), sign of u), and it keeps",
     "its field-of-definition paragraph (LSP Prop 13.2(a)) " + MARK + ", and it keeps"),
]

EDITS[FRAM] = [
    ("conifold Prop 13.2's proof should use the Teleman numbers Theorem (2), Prop 7.7, Prop 8.4, (8.6), Prop 8.7, and the steps",
     "conifold Prop 13.2's proof uses the Teleman numbers Theorem (2), 2.18, 2.19, Prop 8.7, with Prop 8.4, (8.6) for comparison and Prop 7.7 only in a remark, and the steps"),
    ("The field-of-definition paragraph of the conifold Prop 13.2 carries an OPEN mark from round 6c (sign of u); no statement of this file rests on that paragraph.",
     "The field-of-definition paragraph of the conifold Prop 13.2 " + MARK + "; no statement of this file rests on that paragraph."),
]

EDITS[RESU] = [
    ("is unconditional in every genus modulo (MV), textbook facts and, only to identify the graph sum with F_g, the inputs of Proposition 13.2; for the value",
     "is proved in every genus, as a statement about the decoration sum, modulo (MV) and textbook facts; its identification with F_g is Proposition 13.2 " + MARK + ", P modulo the literature inputs listed there; for the value"),
    ("no longer inputs; one referee pass).", "no longer inputs; two referee passes (round 5; round 6, D_VI))."),
    ("- Not refereed yet: s.13.0, s.13.3 (5), and the field-of-definition paragraph of Proposition 13.2 (the part recovering d/dt through uniqueness and the divisor equation goes beyond the critic's item).",
     "- One referee pass (round 6c): s.13.0, s.13.3 (5). The field-of-definition paragraph of Proposition 13.2 had one pass in its round-6c form, which found the sign error (the family is at u = -s H). It was then rewritten with the direct derivation of FRAMEB_PROOF.md A2.8, with Teleman's Prop 7.7 only in a remark: " + MARK + "."),
    ("- Next: a referee of those three passages; proof_M_residue",
     "- Next: the passes on the repaired paragraph of Proposition 13.2; proof_M_residue"),
]

EDITS[NOTE] = [
    ("is therefore unconditional in every genus modulo (MV), textbook facts and, only to identify the graph sum with F_g, the inputs of Proposition 13.2 there.",
     "is therefore proved in every genus, as a statement about the decoration sum, modulo (MV) and textbook facts; its identification with F_g is Proposition 13.2 there " + MARK + ", P modulo the literature inputs listed there."),
    ("Theorem VI of the local-P2 draft (one referee pass), with",
     "Theorem VI of the local-P2 draft (two referee passes (round 5; round 6, D_VI)), with"),
]

EDITS[LEMA] = [
    ("the s.1.5 extension and (H0) are marked unrefereed in the round-5 draft. For the conifold it is leaf_structure_proof s.13.1(i), which still needs the domain fix of round-5 critic item S3(a) (the form of draft s.1.5 and (H0), under the dictionary c -> E).",
     "the s.1.5 extension and (H0) had a pass in round 6 (D_rows, D_hodge). For the conifold it is leaf_structure_proof s.13.1(i), with the domain and (H0) in leaf_structure_proof s.13.0."),
    ("Proposition 13.2_P2 (LP2 draft s.8.3, open) and conifold Proposition 13.2 (citation fixes of round-5 critic item S3(c)).",
     "Proposition 13.2_P2 (FRAMEB_PROOF.md; LP2 draft s.8.3, P modulo literature as stated there, assembled text refereed only as recorded there) and conifold Proposition 13.2 (leaf_structure_proof s.13.3)."),
    ("The draft's s.7.8 is the record of which applies.",
     "The draft's s.7.8 is the record: it uses LLZ eq. (1) only, at mu = 1^a, and neither eq. (5) nor the vertex-operator form."),
    ("(local P2 written in the draft, conifold pending S3(a)).",
     "(local P2 written in the draft, conifold: leaf_structure_proof s.13.0, s.13.1)."),
    ("Still open: Prop. 13.2_P2, the conifold Prop. 13.2 citation fixes (S3(c)), and hence the gap corollary (S5).",
     "Prop. 13.2_P2 and the local P2 gap corollary are now in FRAMEB_PROOF.md (P modulo literature); the conifold Prop. 13.2 is leaf_structure_proof s.13.3."),
]


def main():
    out = {}
    for p, eds in EDITS.items():
        s = load(p)
        for old, new in eds:
            n = s.count(old)
            if n != 1:
                print("MISMATCH", n, p, repr(old[:90]))
                sys.exit(1)
            s = s.replace(old, new)
        out[p] = s
    # the paragraph of Proposition 13.2(a)
    s = out[LEAF]
    a = s.index("The embedding forgets d/dt")
    end = "so the ODE of [T12] Proposition 8.2(ii) along H is the QDE in d/dt = Q d/dQ."
    b = s.index(end, a) + len(end)
    assert s.count("The embedding forgets d/dt") == 1
    assert b - a < 1500
    out[LEAF] = s[:a] + PARAGRAPH + s[b:]
    for p, s in out.items():
        if p == LEAF:
            s.encode("ascii")
        save(p, s)
    print("applied")


main()
