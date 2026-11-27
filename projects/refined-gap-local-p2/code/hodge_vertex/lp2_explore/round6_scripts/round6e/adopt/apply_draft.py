# Round 6e: exact-match edits of LP2_REGULARITY_PROOF_draft.md.
import os
import sys
P = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "LP2_REGULARITY_PROOF_draft.md")
s = open(P, encoding='utf-8', newline='').read()
fails = []


def rep(old, new, count=1):
    global s
    c = s.count(old)
    if c != count:
        fails.append((c, count, old[:120]))
        return
    s = s.replace(old, new)


# line 1 comment
rep("the one major item (report lit, source of eq:J=I at T_max) is OPEN item O1 of FRAMEB_PROOF.md, not fixed here. -->",
    "the one major item (report lit, source of eq:J=I at T_max) is item O1 of FRAMEB_PROOF.md. Round 6e: that citation repaired (CCIT15 Thm 31, report verify:ccit), and the marker ROUND6C-REFEREES of s.0.4 filled in. -->")

# s.0.2 labels
rep("inputs in s.8.2a; conditional on OPEN item O1 (FRAMEB_PROOF.md, head) through Stages A, B]",
    "inputs in s.8.2a; through Stages A, B, also P mod eq:J=I (CCIT15 Thm 31, FB s.1.2), subject to the Q1 proviso of item O1 (FRAMEB_PROOF.md, head)]")
rep("listed in s.8.2a. Conditional on OPEN item O1 (FRAMEB_PROOF.md, head) through Proposition 13.2_P2. Proof in s.7.6.]",
    "listed in s.8.2a. Among those inputs is eq:J=I (CCIT15 Thm 31, FB s.1.2), through Proposition 13.2_P2, subject to the Q1 proviso of item O1 (FRAMEB_PROOF.md, head). Proof in s.7.6.]")
rep("[Same status, including OPEN item O1, plus the polynomiality of Theorem FB_B(b).",
    "[Same status, including eq:J=I (CCIT15 Thm 31) and the Q1 proviso of item O1, plus the polynomiality of Theorem FB_B(b).",
    count=2)

# s.0.2 referee status (coldB 4)
rep("FRAMEB_PROOF.md as assembled has had no referee pass of its own, and neither has the chain of Corollary G as assembled.",
    "The passes of FRAMEB_PROOF.md as assembled, and of the chain of Corollary G, are recorded at the marker ROUND6C-REFEREES in s.0.4.")

# s.0.3
rep("Ruling Q1 names statements on T_F. This use is outside the literal wording of Q1 and is recorded here pending the user's confirmation. **OPEN item O1 (round 6d, report lit; citation level; stated at the head of FRAMEB_PROOF.md).** BS cites eq:J=I to MR2486673 (Coates, arXiv:0810.2200v1), which covers T_F only (Cor 4.2), so the printed source does not supply the input at generic α. A candidate source is Coates-Corti-Iritani-Tseng, Compos. Math. 151 (2015) 1878-1912, Theorem 31, applied to K_P2 x C^2 (FB s.1.2). It has not been adopted. Proposition 13.2_P2, Corollaries G, G_{k,h}, Theorem FB_A (A2), (A3), Theorem FB_B, Proposition 2.6, Theorem 2.7 and Proposition 2.11 are conditional on O1; Theorems LC_P2, R_P2, Corollary SV_P2 and Theorem VI, as statements about the decoration sum, are not.",
    "**Item O1 (citation repaired in round 6d-6e; stated at the head of FRAMEB_PROOF.md).** BS cites eq:J=I to MR2486673 (Coates, arXiv:0810.2200v1), which covers T_F only (Cor 4.2). eq:J=I at T_max is CCIT15 = Coates-Corti-Iritani-Tseng, Compos. Math. 151 (2015) 1878-1912, Theorem 31, applied to K_P2 x C^2 (FB s.1.2), a published all-orders genus-0 theorem whose proof uses localisation and the Lagrangian cone, with no Eynard-Orantin TR, remodelling or other forbidden input. The genus-0 QDE at T_max is derived in FB s.1.2. The T_max use is outside the literal wording of ruling Q1, which names statements on T_F, and awaits the user's confirmation. Until it is given, Proposition 13.2_P2, Corollaries G, G_{k,h}, Theorem FB_A (A2), (A3), Theorem FB_B except (b), Proposition 2.6, Theorem 2.7 and Proposition 2.11 carry this one proviso (the Q1 proviso of item O1); Theorems LC_P2, R_P2, Corollary SV_P2 and Theorem VI, as statements about the decoration sum, do not.")
rep("which covers T_F only (Cor 4.2); at generic α the source is OPEN item O1, with the candidate above; it is not re-proved.",
    "which covers T_F only (Cor 4.2); at generic α the source is Coates-Corti-Iritani-Tseng, Compos. Math. 151 (2015) 1878-1912, Theorem 31, applied to K_P2 x C^2 (FB s.1.2); it is not re-proved.")
rep("eq:J=I at T_max, candidate only: CCIT15, Compos. Math. 151 (2015) 1878-1912, arXiv:1310.4163v5, Theorem 31, Definition 28, eq. (19), eq. (4) (OPEN item O1).",
    "eq:J=I at T_max: CCIT15, Compos. Math. 151 (2015) 1878-1912, arXiv:1310.4163v5, Theorem 31, Definition 28, eq. (19), eq. (4), Remark 3 (item O1; pinned in round 6d, adopted in round 6e). Genus-0 QDE at T_max: derived in FB s.1.2 from Pandharipande, math/9806133v1, Lemma 1, and FB:A1.4.")

# s.0.4
rep("- A valid source for eq:J=I at generic α on T_max: OPEN item O1 (s.0.3), with a candidate pinned in FB s.1.2.",
    "- That the T_max use of eq:J=I (CCIT15 Thm 31, FB s.1.2) falls under ruling Q1: it is outside the literal wording of Q1 and awaits the user's confirmation (item O1, s.0.3).")
rep("(MV), Teleman's Theorem (2) with 2.18, 2.19, and eq:J=I stay cited,",
    "(MV), Teleman's Theorem (2) with 2.18, 2.19, and eq:J=I (CCIT15 Thm 31) stay cited,")
rep("**Not refereed:** FRAMEB_PROOF.md as assembled, and the chain of Corollary G as assembled; the round-6c edits of s.0, s.1.2-s.1.4, s.2.4, s.2.6, s.2.7, s.7.5, s.7.6, s.8, s.9 and the Appendix of this text; the items of FRAMEB_PROOF.md s.1.4, in particular SP1-SP4 line by line and the proof of Lemma 2.4 (s.2.2-2.5 as a whole have not had a pass in their present form; the part of Theorem 2.7(e) that the proof of record uses is FB:B4, which the referees of the reports re-derived by hand);",
    "**Not refereed, beyond the reads recorded at the marker ROUND6C-REFEREES below:** the round-6c edits of s.0, s.1.2-s.1.4, s.2.4, s.2.6, s.2.7, s.7.5, s.7.6, s.8, s.9 and the Appendix of this text; the items of FRAMEB_PROOF.md s.1.4 (refS2 re-derived SP1-SP5 and Lemma 2.4 as used by FB:A3.5, B2, B4 in round 6c; the rest of s.2.2-2.5 has not had a pass in its present form; the part of Theorem 2.7(e) that the proof of record uses is FB:B4, which the referees of the reports re-derived by hand);")

marker_old = "- ROUND6C-REFEREES: placeholder for the main session. The round-6c referee passes on FRAMEB_PROOF.md and on this draft had not run when this line was written, and none is claimed here."
marker_new = ("- ROUND6C-REFEREES: record of the round-6 referee passes on FRAMEB_PROOF.md, on this draft and on the chain of Corollary G. Reports: round6_scripts/round6c_result.json, round6_scripts/round6d_result.json; each report's own not_checked list states its limits. All these passes are by AI referees, not human ones, and the proof has not been read by a human referee. "
              "*Round 6c.* refA (Stage A, with A2.8 still in its earlier form through Teleman's Prop 7.7), refBC (Stages B, C), refS2 (s.2.1-2.5 as used by FB:A3.5, B2, B4), refGap (Corollary G, restricted as its report says: Stage A not re-derived and A2.8-A2.9 not read; B4, SP1-SP5, Lemmas 2.3-2.5, Theorem A and Prop 2.10 S1-S6 not re-derived; Theorems LC_P2 and VI at statement level or low genus), critic (assembled text; draft s.3-s.6, s.7.8, LEMMA_A_PROOF.md s.2-s.6 and the SP1-SP4 calculus used through their statements): all sound with minor fixes, no major item. refCon: one major item, only in the conifold Prop 13.2(a) paragraph (the sign of the family, in the comparison with Teleman's Prop 8.2(ii)), now repaired; no statement of FRAMEB_PROOF.md rests on that paragraph. "
              "*Round 6d.* ref28hand and coldF re-derived the rewritten A2.8 by hand, with ref28num as evidence only (finite order; it tests neither the stability step nor the derivation of M = 0); refNew and coldB also re-derived it once. refNew: one pass, by hand re-derivation, on the text written after round 6c (the A2.1 row-2 paragraph, draft Lemmas 2.3, 2.5 as rewritten, A3.5, the g = 3 constant, the Pi_{gamma_cf} paragraph, the C Step 7 Taylor paragraph and Definition (Cf)); not covered by it: SP1-SP4 line by line, Lemma 2.4, the Barnes constant, Theorems R_P2, SV_P2, VI. refCon2: the repaired rank-2 paragraph (ROUND6_STATUS_2026-11-02.md); this does not raise the status of the rank-3 A2.8. coldF: forward cold read of FRAMEB_PROOF.md. coldB: backward read of the chain of Corollary G. lit: citations. All sound with minor fixes, except lit's citation-level major item (item O1 of FRAMEB_PROOF.md). "
              "*Round 6e.* verify:ccit: one independent reading of CCIT15 Theorem 31 and of its application to Z (statement, definitions, substitution, mirror map, comparison with FB:A1.2, QDE derivation, all by hand; exact check through d = 4 as evidence only). Sound with minor fixes, no major item. It did not re-derive the proof of Theorem 31 and did not read the journal print. "
              "*Not covered by any of these passes:* the proofs of the cited theorems, namely Teleman's Theorem (2), Prop 8.7 and 2.17-2.19; CCIT15 Theorem 31; Pandharipande's Lemma 1; Graber-Pandharipande; (MV); Mumford and FP Theorem 4; PPZ Prop 2.14; Barnes's second lemma (primary source not opened). Finite-order computations in these reports are evidence only. Whether the T_max use falls under ruling Q1 is the user's decision (item O1, s.0.3).")
rep(marker_old, marker_new)

# Prop 2.6, Thm 2.7
rep("[P mod eq:J=I and the genus-0 QDE, the reduced (C1) of s.0.3, with the source of eq:J=I at generic α OPEN (item O1, s.0.3);",
    "[P mod eq:J=I (CCIT15 Thm 31, FB s.1.2; T_max use under Q1 pending the user, s.0.3) and the genus-0 QDE (derived in FB s.1.2), the reduced (C1) of s.0.3;",
    count=2)
rep("eq:J=I (source at generic α: OPEN item O1, s.0.3)", "eq:J=I (CCIT15 Thm 31, FB s.1.2)", count=2)

# Prop 2.11
rep("[P mod the inputs of Theorem FB_B, conditional on OPEN item O1 (s.0.3): Theorem FB_A;",
    "[P mod the inputs of Theorem FB_B, among them eq:J=I (CCIT15 Thm 31), subject to the Q1 proviso of item O1 (s.0.3): Theorem FB_A;")

# s.7.6
rep("Theorem VI, and conditional on OPEN item O1 through Proposition 13.2_P2; end-to-end inputs",
    "Theorem VI, with eq:J=I (CCIT15 Thm 31) through Proposition 13.2_P2, subject to the Q1 proviso of item O1; end-to-end inputs")

# s.8.1 rows
rep("eq:J=I (OPEN item O1, s.0.3; candidate CCIT15 Thm 31), genus-0 QDE; SP1-SP5;",
    "eq:J=I (CCIT15 Thm 31, FB s.1.2; T_max use under Q1 pending the user), genus-0 QDE (derived in FB s.1.2); SP1-SP5;",
    count=2)
rep("= Theorem FB_B(d), hence OPEN item O1; step",
    "= Theorem FB_B(d), hence eq:J=I (CCIT15 Thm 31) and the Q1 proviso of item O1; step")
rep("eq:J=I (OPEN item O1, s.0.3; candidate CCIT15 Thm 31), genus-0 QDE, eq:SHb;",
    "eq:J=I (CCIT15 Thm 31, FB s.1.2; T_max use under Q1 pending the user), genus-0 QDE (derived in FB s.1.2), eq:SHb;")
rep("Theorem FB_A (hence OPEN item O1);",
    "Theorem FB_A (hence eq:J=I, CCIT15 Thm 31, and the Q1 proviso of item O1);")
rep("Stages A, B (hence OPEN item O1);",
    "Stages A, B (hence eq:J=I, CCIT15 Thm 31, and the Q1 proviso of item O1);")
rep("and the inputs of FRAMEB_PROOF.md s.1 (s.8.2a), conditional on OPEN item O1; Prop 13.2_P2,",
    "and the inputs of FRAMEB_PROOF.md s.1 (s.8.2a), among them eq:J=I (CCIT15 Thm 31), Q1 proviso of item O1; Prop 13.2_P2,")
rep("assembled text not refereed (s.0.4)", "assembled text: passes as in s.0.4")
rep("P mod (MV)  same, including OPEN item O1, plus polynomiality",
    "P mod (MV)  same, including eq:J=I (CCIT15 Thm 31) and the Q1 proviso of item O1, plus polynomiality")

# s.8.2a item 1
rep("which covers T_F only (Cor 4.2). At generic α the source is OPEN item O1 (s.0.3; candidate CCIT15 Theorem 31, FB s.1.2, not adopted).",
    "which covers T_F only (Cor 4.2). At generic α the source is CCIT15 Theorem 31 (Coates-Corti-Iritani-Tseng, Compos. Math. 151 (2015) 1878-1912), applied to K_P2 x C^2 (FB s.1.2; item O1, s.0.3); it is not re-proved. The genus-0 QDE at T_max is derived in FB s.1.2 from Pandharipande, math/9806133v1, Lemma 1, and FB:A1.4.")

# s.8.3 Verification
rep("5. *Verification.* The reports had the referee passes listed in s.0.4. The assembled text has not; see the marker ROUND6C-REFEREES in s.0.4.",
    "5. *Verification.* The reports had the referee passes listed in s.0.4. The assembled text had the passes recorded at the marker ROUND6C-REFEREES in s.0.4. These are AI passes; no human referee has read the proof.")

# s.9 item 21
rep("and item 12 (lines 2633, 2658: MR2486673 gives eq:J=I on T_F only; OPEN item O1).",
    "and item 12 (lines 2633, 2658: MR2486673 gives eq:J=I on T_F only; the correct citation at T_max is MR3414388, Theorem 31; item O1).")

if fails:
    for f in fails:
        print('FAIL', f)
    sys.exit(1)
open(P, 'w', encoding='utf-8', newline='').write(s)
print('ok, chars', len(s))
