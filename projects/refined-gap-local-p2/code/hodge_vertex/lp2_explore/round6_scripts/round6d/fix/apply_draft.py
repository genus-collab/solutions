# Exact-match edits of LP2_REGULARITY_PROOF_draft.md, round 6d minor items. All-or-nothing.
import os
import io, sys
P = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "LP2_REGULARITY_PROOF_draft.md")
s = io.open(P, encoding="utf-8", newline="").read()
log = []
def rep(old, new, n=1):
    global s
    c = s.count(old)
    if c != n:
        print("COUNT", c, "expected", n, "for:", old[:120]); sys.exit(1)
    s = s.replace(old, new)
    log.append(old[:60])

O1 = "OPEN item O1 (FRAMEB_PROOF.md, head)"

# head comment
rep("""Later in round 6c: minor items of the round-6c referee reports applied in place, where they could be re-derived. -->""",
"""Later in round 6c: minor items of the round-6c referee reports applied in place, where they could be re-derived. Round 6d: minor items of the round-6d reports applied in place after checking; the one major item (report lit, source of eq:J=I at T_max) is OPEN item O1 of FRAMEB_PROOF.md, not fixed here. -->""")

# s.0.2 labels
rep("""the logarithmic form in Lemma 1.2(c) (DLMF 15.6.7, Gauss), not Barnes's second lemma; inputs in s.8.2a]""",
"""the logarithmic form in Lemma 1.2(c) (DLMF 15.6.7, Gauss), not Barnes's second lemma; inputs in s.8.2a; conditional on """ + O1 + """ through Stages A, B]""")
rep("""and the inputs of FRAMEB_PROOF.md s.1, listed in s.8.2a. Proof in s.7.6.]""",
"""and the inputs of FRAMEB_PROOF.md s.1, listed in s.8.2a. Conditional on """ + O1 + """ through Proposition 13.2_P2. Proof in s.7.6.]""")
rep("""**Corollary G_{k,h}.** [Same status, plus the polynomiality of Theorem FB_B(b). Proof in s.7.6.]""",
"""**Corollary G_{k,h}.** [Same status, including OPEN item O1, plus the polynomiality of Theorem FB_B(b). Proof in s.7.6.]""")
rep("""vanishes linearly at Δ_c = 0 with leading coefficient 1 (BS 929, 937, 955-957).""",
"""vanishes linearly at Δ_c = 0 with leading coefficient 1 (L2b∘θ is eq:pfinf, whose non-constant solutions BS identifies with periods, BS 2755-2760; the expansion is BS 929; A_cf = θ Π_{γ_cf} is BS 937; the asymptotics are BS 955-957).""")
rep("""the integrality wording of BS 925-927 is not used. F -> F_cf in def:dirint is ev_cf on the polynomial of Theorem FB_B(e), whose coefficients lie in""",
"""the integrality wording of BS 927-929 is not used. F -> F_cf in def:dirint is ev_cf on the polynomial f_j of FB:C Step 7 (Theorem FB_B(e) after extracting monomials in ε), whose coefficients lie in""")

# s.0.3
rep("""Ruling Q1 names statements on T_F. eq:J=I at T_max is a genus-0 mirror theorem, with no TR, cited by BS to Coates, Commun. Math. Phys. 287 (2009) 1071-1108 (MR2486673), not opened here. This use is outside the literal wording of Q1 and is recorded here pending the user's confirmation.""",
"""Ruling Q1 names statements on T_F. This use is outside the literal wording of Q1 and is recorded here pending the user's confirmation. **OPEN item O1 (round 6d, report lit; citation level; stated at the head of FRAMEB_PROOF.md).** BS cites eq:J=I to MR2486673 (Coates, arXiv:0810.2200v1), which covers T_F only (Cor 4.2), so the printed source does not supply the input at generic α. A candidate source is Coates-Corti-Iritani-Tseng, Compos. Math. 151 (2015) 1878-1912, Theorem 31, applied to K_P2 x C^2 (FB s.1.2). It has not been adopted. Proposition 13.2_P2, Corollaries G, G_{k,h}, Theorem FB_A (A2), (A3), Theorem FB_B, Proposition 2.6, Theorem 2.7 and Proposition 2.11 are conditional on O1; Theorems LC_P2, R_P2, Corollary SV_P2 and Theorem VI, as statements about the decoration sum, are not.""")
rep("""eq:asymp 2943 and eq:canR 3566 are no longer cited for it. eq:J=I is cited by BS to MR2486673 (Coates, Commun. Math. Phys. 287 (2009) 1071-1108), which was not opened; it is not re-proved.""",
"""eq:asymp 2943 and eq:canR 3566 are no longer cited for it. BS cites eq:J=I to MR2486673 (Coates, arXiv:0810.2200v1), which covers T_F only (Cor 4.2); at generic α the source is OPEN item O1, with the candidate above; it is not re-proved.""")
rep("""  - Teleman, arXiv:0712.0160: Theorem (2), Prop 8.7(ii) (FB:A1.8), and 2.18, 2.19, that is Harer stability with Looijenga, and Madsen-Weiss (FB:A2.8).""",
"""  - Teleman, arXiv:0712.0160: Theorem (2), Prop 8.7(ii) (FB:A1.8), and 2.18, 2.19, that is Harer stability with Looijenga Props (2.1), (2.2), and Madsen-Weiss Theorem 1.1.1 (FB:A2.8; pinned in FB s.1.2).""")
rep("""  - Pandharipande-Pixton-Zvonkine (PPZ), arXiv:1303.1043: Defs""",
"""  - Pandharipande-Pixton-Zvonkine (PPZ), J. Amer. Math. Soc. 28 (2015) 279-309, arXiv:1303.1043v3: Defs""")
rep("""  - Graber-Pandharipande (GP), arXiv:alg-geom/9708001: formula (1), Apps. A, C""",
"""  - Graber-Pandharipande (GP), Invent. Math. 135 (1999) 487-518 (journal data recalled, not verified), arXiv:alg-geom/9708001v2: formula (1), Apps. A, C""")
rep("""(FB:A3.2; the original paper was not opened), and FP Theorem 4 (FB:A3.3).""",
"""(FB:A3.2; original read in round 6d: (5.2), (5.3), (5.4)), and FP Theorem 4 (FB:A3.3).""")
rep("""the primary source was not opened; the form used was confirmed through a secondary statement. DLMF 15.6.7 and Gauss's values of ψ at 1/3, 2/3 (same step).""",
"""the primary source and the books (Bailey s.6.2, AAR Thm 2.4.3) were not opened; numbers and form come from secondary statements listed in FB s.1.2. DLMF 15.6.7 and Gauss's values of ψ at 1/3, 2/3 (same step).
  - Pinned in round 6d, details in FB s.1.2. eq:J=I at T_max, candidate only: CCIT15, Compos. Math. 151 (2015) 1878-1912, arXiv:1310.4163v5, Theorem 31, Definition 28, eq. (19), eq. (4) (OPEN item O1). Stability: Looijenga, J. Algebraic Geom. 5 (1996) 135-150, arXiv:alg-geom/9401005v2; Madsen-Weiss, Ann. of Math. 165 (2007) 843-941, arXiv:math/0212321v3; Harer, Ann. of Math. 121 (1985) 215-249, and Ivanov, Leningrad Math. J. 1 (1990) 675-691, not opened. Mumford, Progr. Math. 36 (1983) 271-328. FP, journal data recalled as Invent. Math. 139 (2000) 173-199, not verified; c(𝔼) c(𝔼^vee) = 1 is in s.2 before (22). Comparison only, not an input: Lho-Pandharipande, Adv. Math. 332 (2018) 349-402, arXiv:1702.06096v4. (MV): LLZ print read through the EMIS mirror, eq. (1) on p. 291; Okounkov-Pandharipande, Geom. Topol. 8 (2004) 675-699, Theorem 1.""")

# s.0.4
rep("""- A literature-free proof of Lemma A,""",
"""- A valid source for eq:J=I at generic α on T_max: OPEN item O1 (s.0.3), with a candidate pinned in FB s.1.2.
- A literature-free proof of Lemma A,""")
rep("""Remark 2.13(a)-(c); Lemma 1.4 at genus 4.""",
"""Remark 2.13(b), (c) and the audit text of (a) (the mathematics of (a) is the row-2 paragraph of FB:A2.1, which is self-contained and has the round-6d pass recorded at ROUND6C-REFEREES); Lemma 1.4 at genus 4.""")

# Lemma 2.3
rep("""The constant c_0 comes from dΩ = (1/3) d^3x and det Hess 𝒲 = W/3 (Lemma 2.1). ∎""",
"""This finite sum, divided by its z^0 term, defines Φ; it is the series Phi of FB, Theorem FB_B(a). The constant c_0 comes from dΩ = (1/3) d^3x, det Hess(-𝒲) = -W/3 and det Hess 𝒲 = W/3 (Lemma 2.1). ∎""")

# Lemma 2.5
rep("""The θ-constants are then F, and f, g are series in z with coefficients in F[[y]].
""",
"""The θ-constants are then F, and f, g are series in z with coefficients in F[[y]]. In (i), W^{-1/2} need not lie in F[[y]]. The statement there is that every non-zero element of ker 𝒫_1 on F[[y]] is a unit (compare y-orders, using θ log W in y F[[y]]) and that two such elements have ratio in F, which is all (ii) uses. With F = K^ of FB:A3.5, ker 𝒫_1 = F (3W)^{-1/2}.
""")

# Prop 2.6, Thm 2.7, Prop 2.11
rep("""**Proposition 2.6 (R2:Prop 5).** [P mod eq:J=I and the genus-0 QDE, the reduced (C1) of s.0.3;""",
"""**Proposition 2.6 (R2:Prop 5).** [P mod eq:J=I and the genus-0 QDE, the reduced (C1) of s.0.3, with the source of eq:J=I at generic α OPEN (item O1, s.0.3);""")
rep("""4. By Lemma 2.3 at y = 0 and SP5:""",
"""4. By the first sentence of Lemma 2.3, applied to 𝒲_y of step 2, and SP5:""")
rep("""[P mod eq:J=I and the genus-0 QDE, the reduced (C1) of s.0.3; (C3) is discharged:""",
"""[P mod eq:J=I and the genus-0 QDE, the reduced (C1) of s.0.3, with the source of eq:J=I at generic α OPEN (item O1, s.0.3); (C3) is discharged:""")
rep("""**Proposition 2.11 (Givental-Teleman on T_F; replaced by Theorem FB_B of FRAMEB_PROOF.md).** [P mod the inputs of Theorem FB_B: Theorem FB_A;""",
"""**Proposition 2.11 (Givental-Teleman on T_F; replaced by Theorem FB_B of FRAMEB_PROOF.md).** [P mod the inputs of Theorem FB_B, conditional on OPEN item O1 (s.0.3): Theorem FB_A;""")

# s.7.6
rep("""its field-of-definition paragraph (LSP Prop 13.2(a)) [repaired in round 6d; the repaired paragraph has had the passes recorded in the round 6 status file], and it keeps a BS citation.""",
"""its field-of-definition paragraph (LSP Prop 13.2(a)) was repaired in round 6d (passes on the repaired paragraph: ROUND6_STATUS_2026-11-02.md, same folder), and it keeps a BS citation.""")
rep("""**Corollary G (BS rmk:cfgap).** [Status as in s.0.2: P mod Proposition 13.2_P2, Lemma 1.2(c), Theorem R_P2, Corollary SV_P2, Theorem VI; end-to-end inputs""",
"""**Corollary G (BS rmk:cfgap).** [Status as in s.0.2: P mod Proposition 13.2_P2, Lemma 1.2(c), Theorem R_P2, Corollary SV_P2, Theorem VI, and conditional on OPEN item O1 through Proposition 13.2_P2; end-to-end inputs""")
rep("""at c1 = C1cf, up to the constant GW^{T_F}_{g,0} on the d > 0 convention. Each T(Γ, c)""",
"""at c1 = C1cf. Each T(Γ, c)""")
rep("""GW^cf_g is such an element, up to the constant, and t_c^{2-2g}""",
"""GW^cf_g is such an element, and t_c^{2-2g}""")
rep("""**Corollary G_{k,h}.** [Same status, plus the polynomiality of Theorem FB_B(b).] F^cf_{k,h}""",
"""**Corollary G_{k,h}.** [Same status, including OPEN item O1, plus the polynomiality of Theorem FB_B(b).] F^cf_{k,h}""")

# s.8.1 rows
rep("""    Prop 2.6 (normalisation; eq:normNLR derived)     P mod C     eq:J=I, genus-0 QDE; SP1-SP5;""",
"""    Prop 2.6 (normalisation; eq:normNLR derived)     P mod C     eq:J=I (OPEN item O1, s.0.3; candidate CCIT15 Thm 31), genus-0 QDE; SP1-SP5;""")
rep("""    Thm 2.7 (continuity, closed form, Rec = Wick)    P mod C     eq:J=I, genus-0 QDE; SP1-SP5;""",
"""    Thm 2.7 (continuity, closed form, Rec = Wick)    P mod C     eq:J=I (OPEN item O1, s.0.3; candidate CCIT15 Thm 31), genus-0 QDE; SP1-SP5;""")
rep("""    Theorem FB_A (CohFT, R, constants, generic α)    P mod C     eq:J=I, genus-0 QDE, eq:SHb; Teleman""",
"""    Theorem FB_A (CohFT, R, constants, generic α)    P mod C     eq:J=I (OPEN item O1, s.0.3; candidate CCIT15 Thm 31), genus-0 QDE, eq:SHb; Teleman""")
rep("""A2.8 as rewritten in round 6c (no Prop 7.7): no further pass""",
"""A2.8 as rewritten in round 6c (no Prop 7.7): two independent hand re-derivations, recorded in FB:A2.8""")
rep("""    Theorem FB_B (specialisation to T_F)             P mod C     Theorem FB_A; SP1-SP5,""",
"""    Theorem FB_B (specialisation to T_F)             P mod C     Theorem FB_A (hence OPEN item O1); SP1-SP5,""")
rep("""    Theorem FB_C (frame change, dictionary)          P mod C     Stages A, B; Prop 2.10; Mumford through A3.2; PPZ Defs 2.2, 2.5, 2.13, Prop 2.14;""",
"""    Theorem FB_C (frame change, dictionary)          P mod C     Prop 2.10; Mumford through A3.2; PPZ Defs 2.2, 2.5, 2.13, Prop 2.14; Stages A, B only to identify Tbar_g with the reduction of 𝒯_g;""")
rep("""    Proposition 13.2_P2 (G_g = GW^cf_g, d >= 0)      P mod C     Stages A, B; Props 2.10, 2.12;""",
"""    Proposition 13.2_P2 (G_g = GW^cf_g, d >= 0)      P mod C     Stages A, B (hence OPEN item O1); Props 2.10, 2.12;""")
rep("""    Corollary G (BS rmk:cfgap)                       P mod (MV)  and the inputs of FRAMEB_PROOF.md s.1 (s.8.2a); Prop 13.2_P2, Lemma 1.2(c),""",
"""    Corollary G (BS rmk:cfgap)                       P mod (MV)  and the inputs of FRAMEB_PROOF.md s.1 (s.8.2a), conditional on OPEN item O1; Prop 13.2_P2, Lemma 1.2(c),""")
rep("""    Corollary G_{k,h} (BS def:dirint, conifold)      P mod (MV)  same, plus polynomiality of Theorem FB_B(b)""",
"""    Corollary G_{k,h} (BS def:dirint, conifold)      P mod (MV)  same, including OPEN item O1, plus polynomiality of Theorem FB_B(b)""")

# s.8.2a
rep("""Corollary G = Proposition 13.2_P2 + Lemma 1.2(c) + Corollary 2.8 + Theorem R_P2 + Corollary SV_P2 + Theorem VI. Details: FRAMEB_PROOF.md s.0.3, s.1.""",
"""Corollary G = Proposition 13.2_P2 + Lemma 1.2(c) + Corollary 2.8 + Theorem R_P2 + Corollary SV_P2 + Theorem VI, with Theorem FB_B(e) for the clause on Taylor coefficients. Details: FRAMEB_PROOF.md s.0.3, s.1.""")
rep("""eq:J=I is cited by BS to MR2486673 (Coates, Commun. Math. Phys. 287 (2009) 1071-1108), not opened. Definitions only:""",
"""BS cites eq:J=I to MR2486673 (Coates, arXiv:0810.2200v1), which covers T_F only (Cor 4.2). At generic α the source is OPEN item O1 (s.0.3; candidate CCIT15 Theorem 31, FB s.1.2, not adopted). Definitions only:""")

# s.9
rep("""and item 11 (lines 925-927: the expansions""",
"""item 11 (lines 927-929: the expansions""")
rep("""and hold for an integral cycle only up to a constant prefactor).
""",
"""and hold for an integral cycle only up to a constant prefactor), and item 12 (lines 2633, 2658: MR2486673 gives eq:J=I on T_F only; OPEN item O1).
""")

io.open(P, "w", encoding="utf-8", newline="").write(s)
print("draft edits applied:", len(log))
