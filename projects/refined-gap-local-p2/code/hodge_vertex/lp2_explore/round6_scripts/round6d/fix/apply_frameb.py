# Exact-match edits of FRAMEB_PROOF.md, round 6d minor items. All-or-nothing: writes only if every old string has the stated count.
import os
import io, sys
P = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "FRAMEB_PROOF.md")
s = io.open(P, encoding="utf-8", newline="").read()
log = []
def rep(old, new, n=1):
    global s
    c = s.count(old)
    if c != n:
        print("COUNT", c, "expected", n, "for:", old[:120]); sys.exit(1)
    s = s.replace(old, new)
    log.append(old[:60])

# ---------- Head
rep("""Minor items of the round-6c referee reports are applied in place where they could be re-derived; A2.8 is the one place where a proof route changed.""",
"""Minor items of the round-6c referee reports are applied in place where they could be re-derived; A2.8 is the one place where a proof route changed. Minor items of the round-6d reports ref28hand, ref28num, refNew, refCon2, lit, coldF, coldB are applied in place after checking against this text. The one major item of round 6d (report lit) is not fixed here: it is OPEN item O1 below.""")

rep("""Ruling Q1 names statements on T_F. eq:J=I at T_max is a genus-0 mirror theorem, with no TR, cited by BS to Coates, Commun. Math. Phys. 287 (2009) 1071-1108 (MR2486673), not opened here. This use is outside the literal wording of Q1 and is recorded here pending the user's confirmation.
""",
"""Ruling Q1 names statements on T_F. This use is outside the literal wording of Q1 and is recorded here pending the user's confirmation.

**OPEN item O1 (round 6d, report lit; citation level; no mathematics of this file changes).** The chain needs eq:J=I at generic alpha on T_max, and the source printed in BS does not give it. BS 2658 cites MR2486673 = Coates, Commun. Math. Phys. 287 (2009) 1071-1108, arXiv:0810.2200v1, opened in round 6d: its torus covers the trivial action on P2, and its Prop 4.1, Cor 4.2 are eq:J=I on T_F only (lambda = E, p = H + E/3). Stage A cannot be moved to T_F, because the Q = 0 product is not semisimple there. A candidate source is pinned in s.1.2 (CCIT15 Theorem 31, a published genus-0 theorem whose proof uses localisation and the Lagrangian cone, with no TR, remodelling or other forbidden input). It is recorded as a candidate and has not been adopted. The genus-0 QDE at T_max is also not pinned to a numbered theorem (CCIT15 Theorem 2 with eq. (6) was read as a source; the derivation through genus-0 TRR and WDVV was not written out). Until O1 is closed, and the T_max use is confirmed under Q1, every statement below that rests on eq:J=I is conditional on O1: A2.1, A2.4-A2.6, A2.9, A3, parts (A2), (A3) and the row form of the graph sum of Theorem FB_A; Theorem FB_B; Proposition 13.2_P2; Corollaries G and G_{k,h}; the discharge paragraph at the end of s.2. Not affected: A1, A2.3, A2.8, A3.2, the identity Tbar_g[𝖥] = G_g[𝖥] of Theorem FB_C (C1), (C2), Lemma (t_cf), and draft Theorems LC_P2, R_P2, Corollary SV_P2, Theorem VI as statements about the decoration sum.
""")

# ---------- s.0.2 labels
rep("""**Theorem FB_A (generic weights).** [P mod BS eq:J=I and the genus-0 QDE, with eq:SHb; Teleman Theorem (2), Prop 8.7(ii), and 2.18, 2.19 (Harer stability with Looijenga; Madsen-Weiss); PPZ Defs 2.2, 2.13, Prop 2.14; Graber-Pandharipande;""",
"""**Theorem FB_A (generic weights).** [P mod BS eq:J=I and the genus-0 QDE (eq:SHb is derived from them in A2.1), with the source of eq:J=I at generic alpha OPEN (item O1; candidate CCIT15 Thm 31, s.1.2), so that (A2), (A3) and the row form of the graph sum are conditional on O1; Teleman Theorem (2), Prop 8.7(ii), and 2.18, 2.19 (Harer stability with Looijenga; Madsen-Weiss); PPZ Defs 2.2, 2.5, 2.13, Prop 2.14; Graber-Pandharipande;""")
rep("""that derivation has had one hand re-derivation and no further pass.]""",
"""that derivation has had two independent hand re-derivations, recorded in A2.8.]""")
rep("""c(𝔼)c(𝔼^vee) = 1 (FP s.1.1), for the constant term. One referee pass, conditional on Stage A.]""",
"""c(𝔼)c(𝔼^vee) = 1 (Mumford (5.4), (5.5); FP s.2 before (22)), for the constant term. One referee pass, conditional on Stage A. Conditional on OPEN item O1 through Theorem FB_A, except (b).]""")
rep("""**Theorem FB_C (C1): frame change, for all g >= 2.** [P mod Stages A, B and draft Prop 2.10, Mumford (through A3.2), PPZ Defs 2.2, 2.5, 2.13, Prop 2.14.]""",
"""**Theorem FB_C (C1): frame change, for all g >= 2.** [P mod draft Prop 2.10, A3.2 (Mumford as in FP), PPZ Defs 2.2, 2.5, 2.13, Prop 2.14. Stages A, B are needed only to identify Tbar_g[𝖥] with the reduction of 𝒯_g, B6, not for this identity.]""")
rep("""T'(z) = z(1 - R'^{-1}(z)1) = Sum_i e_i z(1 - S(v_i,-z)) = Sum_i e_i Sum_{d>=2} Dil_i(d) z^d, with vanishing z^0 and z^1 terms. The code's q_t^{(i)} = -psi S(v_i,-psi) equals T'_i - psi, so -psi is exactly the dilaton shift (Dil_i(1) = -1). The edge is (1.1) of the draft, and E_ij[𝖥](psi,w) = E_ij[0](psi,w) + 3 eps1 eps2 Delta_c 𝖥 Q_i(-psi) Q_j(-w).""",
"""T'(z) = z(1 - R'^{-1}(z)1) = Sum_i e_i z(1 - S(v_i,-z)) = Sum_i e_i Sum_{d>=2} Dil^_i(d) z^d, with Dil^_i(d) := (-1)^d [z^{d-1}] S(v_i,z) as in draft s.1.4, and with vanishing z^0 and z^1 terms. The code's q_t^{(i)} = -psi S(v_i,-psi) equals T'_i - psi, so -psi is exactly the dilaton shift (Dil^_i(1) = -1). The edge is (1.1) of the draft, and Ehat_ij[𝖥](psi,w) = Ehat_ij[0](psi,w) + 3 eps1 eps2 Delta_c 𝖥 Q_i(-psi) Q_j(-w), with Ehat as in (1.1) of the draft; in the point frame the 𝖥-term of Ebar carries the additional factor h_w(-psi) h_w(-w).""")
rep("""for the clause F_cf = XX only, the logarithmic form in Lemma (t_cf) (DLMF 15.6.7, Gauss), not Barnes's second lemma.] For g >= 2, Cf(F_g^{T_F})""",
"""for the clause F_cf = XX only, the logarithmic form in Lemma (t_cf) (DLMF 15.6.7, Gauss), not Barnes's second lemma. Conditional on OPEN item O1 through Stages A, B.] For g >= 2, Cf(F_g^{T_F})""")
rep("""The latter three are P mod (MV) = LLZ eq. (1) and textbook facts, s.1.3.] For g >= 2, with GW^cf_g""",
"""The latter three are P mod (MV) = LLZ eq. (1) and textbook facts, s.1.3. Conditional on OPEN item O1 through Proposition 13.2_P2.] For g >= 2, with GW^cf_g""")
rep("""**Corollary G_{k,h}.** [Same status, plus the polynomiality of Theorem FB_B(b).] For k+h >= 2,""",
"""**Corollary G_{k,h}.** [Same status, including OPEN item O1, plus the polynomiality of Theorem FB_B(b).] For k+h >= 2,""")
rep("""BS's Pi_{gamma_cf} solves L2b o theta = 0 and vanishes linearly at xi = 0 with leading coefficient 1 (BS 929, 937, 955-957).""",
"""BS's Pi_{gamma_cf} solves L2b o theta = 0 and vanishes linearly at xi = 0 with leading coefficient 1 (L2b o theta is eq:pfinf, whose non-constant solutions BS identifies with periods, BS 2755-2760; the expansion is BS 929; A_cf = theta Pi_{gamma_cf} is BS 937; the asymptotics are BS 955-957).""")
rep("""the integrality wording of BS 925-927 is not used. F -> F_cf in def:dirint is ev_cf on the polynomial of Theorem FB_B(e), whose coefficients lie in""",
"""the integrality wording of BS 927-929 is not used. F -> F_cf in def:dirint is ev_cf on the polynomial f_j of C Step 7 (Theorem FB_B(e) after extracting monomials in eps), whose coefficients lie in""")

# ---------- s.0.3
rep("""Corollary G = Proposition 13.2_P2 + Lemma (t_cf) + Corollary 2.8 + Theorem R_P2 + Corollary SV_P2 + Theorem VI.

End-to-end inputs of Corollary G: eq:J=I and the genus-0 QDE, with eq:SHb; Teleman""",
"""Corollary G = Proposition 13.2_P2 + Lemma (t_cf) + Corollary 2.8 + Theorem R_P2 + Corollary SV_P2 + Theorem VI, with Theorem FB_B(e) for the clause on Taylor coefficients.

End-to-end inputs of Corollary G: eq:J=I (at generic alpha on T_max: OPEN item O1; candidate CCIT15 Thm 31, s.1.2) and the genus-0 QDE (eq:SHb is derived from them in A2.1); Teleman""")

# ---------- s.1.1
rep("""read and used. eq:J=I is cited by BS to MR2486673 (Coates, Commun. Math. Phys. 287 (2009) 1071-1108), which was **not opened**. It is not re-proved.""",
"""read and used. BS cites eq:J=I to MR2486673 (Coates, Commun. Math. Phys. 287 (2009) 1071-1108, arXiv:0810.2200v1, opened in round 6d). That paper proves it only on T_F (Prop 4.1, Cor 4.2). **OPEN item O1:** at generic alpha this citation does not supply the input; the candidate source CCIT15 Theorem 31 (s.1.2) has not been adopted. It is not re-proved.""")

# ---------- s.1.2
rep("""### 1.2 Literature

- **Teleman, arXiv:0712.0160.**""",
"""### 1.2 Literature

- **eq:J=I at generic alpha on T_max: candidate source (OPEN item O1; not adopted).** Coates-Corti-Iritani-Tseng (CCIT15), A mirror theorem for toric stacks, Compos. Math. 151 (2015) 1878-1912 (MR3414388), arXiv:1310.4163v5, Theorem 31 with Definition 28, eq. (19), and eq. (4) for valued points of the cone. Opened at arXiv v5 in round 6d; the journal print was not read. Application: the smooth toric variety Z = K_P2 x C^2 with divisor classes u_0, ..., u_5 of charges (1,1,1,-3,0,0), S empty, T = (C^*)^5. Z has a fixed point and is projective over the affine C^3/mu_3 x C^2, so it is semi-projective. Put H := -u_3/3, alpha_i := u_i - H (i = 0,1,2), eps_j := -u_{3+j}; these are the weights of BS at about line 2383, and alpha_i, eps_j form a basis of H^2_T(pt;Q). In (19) take t_3 = -t/3 and the other t_i = 0, with y = Q_nov e^t. The modification factor is Prod_{m=0}^{3d-1}(-3H-mz) for i = 3, 1/Prod_{m=1}^d(H+alpha_i+mz) for i = 0,1,2, and 1 for i = 4,5. After removing Q_nov^{H/z} from both sides, (19) is eq:Ifun (checked by hand against BS 2621 in round 6d). By eq:Ian, I(y,-z) = -z + I_0^{[1]} H + O(1/z), with I_0^{[1]} vanishing at t = Q_nov = 0, so (4) with t(z) = I_0^{[1]} H gives I(y,z) = J(H I_0^{[1]}, z). The divisor equation then gives eq:J=I with Q = exp I_0^{[1]}. No reduction from K_P2 is needed. If K_P2 with T = (C^*)^3 is used instead, the reduction is immediate in genus 0: R pi_* f^* O_w = O_w on genus-0 trees, so each invariant is multiplied by 1/w and eta^{ab} by w, and J does not change. Other sources read in round 6d and found not to cover T_max: MR2486673 = Coates, arXiv:0810.2200v1, Prop 4.1 and Cor 4.2, prove eq:J=I on T_F only (lambda = E, p = H + E/3), since its torus covers the trivial action on P2. Coates-Givental (Ann. of Math. 165 (2007) 15-53, math/0110142) Theorem 2 and CCIT (Duke Math. J. 147 (2009) 377-438, math/0702234) Theorem 4.6 have T acting on fibres only. Givental alg-geom/9701016 Theorem 0.2 has compact X and non-negative bundles. Givental alg-geom/9603021 Theorem 9.5 gives the T-equivariant small J-function of P^n only. The proofs of these theorems were not re-checked. The T-equivariant setting of CCIT15 was accepted as matching BS's localisation definition by reading its s.2.2-2.3 and A1.1, A1.2 here, with no comparison of the two definitions on Z written out.
- **Teleman, arXiv:0712.0160.**""")
rep("""2.17 (Harer stability), 2.18 (Looijenga), 2.19 (Madsen-Weiss), with the originals not opened.""",
"""2.17 (Harer stability), 2.18 (Looijenga), 2.19 (Madsen-Weiss). Originals pinned in round 6d. Looijenga, J. Algebraic Geom. 5 (1996) 135-150, arXiv:alg-geom/9401005v2 (opened), Props (2.1), (2.2): H^*(M_g;Q)[u_i] -> H^*(M_g^X;Q) is an isomorphism in degree <= N(g). M_g^X is the open moduli space of curves with distinct marked points, u_i = c_1 of the relative tangent sheaf = -psi_i, and N(g) >= g/2 - 1 (Ivanov, as quoted in Looijenga s.1). Madsen-Weiss, Ann. of Math. 165 (2007) 843-941, arXiv:math/0212321v3 (opened), Theorem 1.1.1 and the deduction of Mumford's conjecture after it; Teleman's item [MW] is their ECM 2005 survey, not this paper. Harer, Ann. of Math. 121 (1985) 215-249, and Ivanov, Leningrad Math. J. 1 (1990) 675-691: not opened; their ranges are as reported by Looijenga s.1 and Wahl, arXiv:1006.4476v2, s.1. Range: for integer k, k < (g-1)/2 is equivalent to k <= g/2 - 1, which is Ivanov's bound. It is conservative against Wahl Theorems 1.1, 1.2 (Randal-Williams, arXiv:0909.4278v4, s.1.3; 3k <= 2g - 2, closed case 3k <= 2g). A2.8 needs only some range growing with g, so the choice of bound plays no role. Madsen-Weiss cannot be removed from the chain: Theorem (2) of Teleman uses it (his abstract and (1.5)), although A2.8 itself uses only injectivity.""")
rep("""- **Pandharipande-Pixton-Zvonkine, arXiv:1303.1043.** Opened (provers, referees A2, C, C2). Numbers confirmed:""",
"""- **Pandharipande-Pixton-Zvonkine, J. Amer. Math. Soc. 28 (2015) 279-309, arXiv:1303.1043v3.** Opened (provers, referees A2, C, C2; numbers as in the arXiv version). Numbers confirmed:""")
rep("""Definitions 2.2, 2.5, 2.13; Propositions 2.8, 2.9, 2.14; the line on "the same convention" as Givental and Teleman;""",
"""Definitions 2.2, 2.5, 2.13; Propositions 2.3, 2.8, 2.9, 2.14 (left group action of the unit-preserving action); the line on "the same convention" as Givental and Teleman;""")
rep("""- **Graber-Pandharipande, arXiv:alg-geom/9708001v2.** Provers did not open it.""",
"""- **Graber-Pandharipande, Invent. Math. 135 (1999) 487-518 (journal data recalled, not verified), arXiv:alg-geom/9708001v2.** Numbers as in the arXiv version: formula (1) in the Introduction; Appendix A (global embeddings), B (obstruction theory), C (localization for DM stacks). Provers did not open it.""")
rep("""- **Mumford.** The original paper was **not opened**. The theorem is used""",
"""- **Mumford.** Original opened in round 6d (OCR scan via INSPIRE): D. Mumford, Towards an enumerative geometry of the moduli space of curves, Arithmetic and Geometry II, Progr. Math. 36 (1983) 271-328, formula (5.2), Corollary (5.3) (ch_{2k} 𝔼 = 0), (5.4) with (5.5) (c(𝔼) c(𝔼^vee) = 1 on M_g, extended to Mbar_g), pp. 304-307. It treats unmarked curves, uses another Bernoulli convention, and FP s.1.1 records a typographic error in its boundary term, so the form of record stays FP's. The theorem is used""")
rep("""- **Faber-Pandharipande, math/9810173.** Opened: Theorem (Mumford) in s.1.1; Theorem 4; c(𝔼) c(𝔼^vee) = 1.""",
"""- **Faber-Pandharipande, Hodge integrals and Gromov-Witten theory, math/9810173 (journal data recalled as Invent. Math. 139 (2000) 173-199, not verified).** Opened: Theorem (Mumford) in s.1.1; Theorem 4, which reads Int_{Mbar_g} lambda_{g-1}^3 = (|B_{2g}|/2g)(|B_{2g-2}|/(2g-2))/(2g-2)!, g >= 2; c(𝔼) c(𝔼^vee) = 1, quoted in s.2 just before eq. (22). With lambda_{g-1}^2 = 2 lambda_g lambda_{g-2} this gives Int lambda_{g-2} lambda_{g-1} lambda_g = half that value.""")
rep("""- **Barnes's second lemma.** Barnes, Quart. J. Math. 41 (1910), 136-140: **primary source not opened**. It is not in DLMF s.5.13 (which has only the beta integral 5.13.3). Confirmed through a secondary statement (Wikipedia, "Barnes integral") and numerically (s.A.3). The form used is in C Step 6.""",
"""- **Barnes's second lemma.** E. W. Barnes, Quart. J. Math. 41 (1910) 136-140: **not opened**. It is not in DLMF s.5.13 (which has only the beta integral 5.13.3). Bailey, Generalized Hypergeometric Series, Cambridge Tract 32 (1935), s.6.2, pp. 42-43: section and pages taken from MathWorld (BarnesLemma), book not opened. Andrews-Askey-Roy, Special Functions (1999), Theorem 2.4.3: number taken from arXiv:2107.13871v2 App. A (A.3), book not opened. The form in C Step 6 is the one displayed on MathWorld (beta_1 = d, beta_2 = e, beta_1 + beta_2 = alpha_1 + alpha_2 + alpha_3 + 1) and on Wikipedia (Barnes integral). The symmetric form is Jantzen, arXiv:1211.2637, eq. (2) (opened), with lambda = (1/3, 2/3), rho = (0,0,0), z = -s, 0 < Re z < 1/3. Substitution checked by hand: under s = -(sigma + 2/3) the orientation change is absorbed by ds = -dsigma, so the upward contour in s becomes the upward contour in sigma with sign +, and -2/3 < Re sigma < -1/3 corresponds to -1/3 < Re s < 0. Numerical confirmation is in s.A.3 [E]. The constant t_c = t_cf rests on this lemma and numerics are not a proof step, so before submission the form should be compared with a printed source.""")
rep("""Givental's Theorem 1 does not apply directly, since Z is not compact. Behrend-Manin axioms: not opened.
""",
"""Givental's Theorem 1 does not apply directly, since Z is not compact. Behrend-Manin axioms: not opened.
- **Lho-Pandharipande, Stable quotients and the holomorphic anomaly equation, Adv. Math. 332 (2018) 349-402 (BS LP18:HAE), arXiv:1702.06096v4. Opened; comparison only.** They work at (19): lambda_i = zeta^i lambda_0, so e_1(alpha) = e_2(alpha) = 0, with alpha_i = -lambda_i, q = y, and the K_P2 factor only. Their (24) is eq:Ifun. Their rows (27), (30) have the form of A2.1, A2.5 at these weights: Y_i = lambda_i L, L = (1+27y)^{-1/3}; W(Y_i) = 3 e_3(alpha) is constant, so Psi_{0i} is constant, matching R_0 = 1; S(H), S(H^2) carry L lambda_i/C_1 and L^2 lambda_i^2/(C_1 C_2), with C_1 = A, C_2 = D. Their (29), C_1^2 C_2 = X, is the relation A^2 D = X of Prop 2.12. In Lemma 8 the coefficient of X_LP := D C_1/C_1 = F in R_{2k} is -R_{1,k-1}/L, which is the term -z F S(H) from S(H^2) = z theta_y S(H)/D at e_1 = 0 (A2.1), checked by hand. Their R_k vanish at q = 0 for k >= 1, with the Hodge class (36) kept at the vertex; this is R~ in R = R~ diag(h_wt) of A3.6, and (36) matches A3.1 with the C^2 weights removed. Their route is stable quotients, wall-crossing and Givental, not Teleman, and only at this subtorus; it is a cross-check of form and constants, not of the C^2 factor or of T_F. R_1 = (1-L^2)/18, R_2 = (1-24L-2L^2+25L^4)/648 were compared exactly with eq:Ifun at H = lambda_0 = 1 through y^8 [E].
""")

# ---------- s.1.3
rep("""p. 291, numbered as in arXiv:math/0306434v3), with Okounkov-Pandharipande Theorem 1 as a second published proof.""",
"""p. 291, numbered as in arXiv:math/0306434v3; JDG print opened in round 6d through the EMIS mirror: (1)-(3) on printed p. 291, and the name Marino-Vafa formula for (1) on p. 292. LLZ Theorems 1, 2 are the cut-and-join equations, not (MV). Okounkov-Pandharipande: Geom. Topol. 8 (2004) 675-699, Theorem 1, opened), with Okounkov-Pandharipande Theorem 1 as a second published proof.""")

# ---------- s.1.4
rep("""The direct derivation now in A2.8 has had no further pass.""",
"""The direct derivation now in A2.8 has had two independent hand re-derivations, recorded in A2.8.""")
rep("""- The class-level identity Omega = R.omega beyond integrated numbers (n <= 1 in the round-6b runs, no psi insertions at d > 0) was not tested numerically and rests on Teleman. The R versus R^* test of s.A.4 reaches n = 5 at genus 0 and n = 2 at genus 1, still with no psi insertions.
""",
"""- The class-level identity Omega = R.omega rests on Teleman. Evidence only: the round-6d computation of s.A.1 pairs it with psi classes pulled back from Mbar_{g,n} and with powers of lambda_1 at d <= 3 (g <= 2) and d <= 2 (g = 3, n = 1, psi only). It was not paired with boundary strata or with general tautological classes, and it is not a proof step.
""")

# ---------- A2.1
rep("""**A2.1 Rows are formal [P mod eq:J=I and the genus-0 QDE].**""",
"""**A2.1 Rows are formal [P mod eq:J=I and the genus-0 QDE; source of eq:J=I at generic alpha: OPEN item O1].**""")
rep("""This is the content of draft Remark 2.13(a), used here as a proof step.""",
"""Draft Remark 2.13(a) records the same values as an audit of BS; this paragraph is the proof and does not use the remark.""")

# ---------- A2.8
rep("""**A2.8 Flatness of R [P mod Teleman Theorem (2) through A1.8, and Teleman 2.18, 2.19 (Harer stability with Looijenga; Madsen-Weiss)].**""",
"""**A2.8 Flatness of R [P mod Teleman Theorem (2) through A1.8, A1.5 (fundamental-class axiom), PPZ Defs 2.2, 2.13 and the symplectic group (5), Mumford and FP Thm 4 through A1.7 in A1.8, and Teleman 2.17-2.19 (Harer stability with Looijenga; Madsen-Weiss; originals as in s.1.2)].**""")
rep("""It has had one hand re-derivation and no further pass. In A2.8 only,""",
"""It has had two independent hand re-derivations (the proposer's and the round-6d referee's, with PPZ Defs 2.2, 2.13 and Teleman 2.15-2.19, 7.1, 7.6-7.8, 8.2 compared against the arXiv text). The round-6d reads refNew, coldF and coldB also re-derived it by hand and found only the notation corrected here. The finite-order tests of s.A.1 are evidence, not passes. In A2.8 only,""")
rep("""with D_i := Delta_i^{(2g-2+n)/2} and N_i(kappa) := Sum_m (1/m!) p_{m*} prod_k (Delta_i^{1/2} f_i^* T(psi_{n+k})), a polynomial in the kappa_j with constant term 1.""",
"""with D_i := Delta_i^{(2g-2+n)/2} and N_i(kappa) := Sum_m (1/m!) p_{m*} prod_k (Delta_i^{1/2} f_i^* T(psi_{n+k})), a power series in the kappa_j, j >= 1 (finite in each degree, and given by exp(Sum_j s^i_j kappa_j) with exp(-Sum_j s^i_j psi^j) = b_i(psi) := Delta_i^{1/2} a^i(psi;1), by T(z) = z(1 - R^{-1}(z)1)), with constant term 1, the same series for every (g,n).""")
rep("""On C: psi_l = p^*psi_l + delta_l, psi_{n+1} delta_l = 0, delta_l delta_l' = 0 for l != l', kappa_j = p^*kappa_j + psi_{n+1}^j, and the only graphs with an edge are the bubbles (0,3) carrying legs l, n+1. So B(psi_l) = B(p^*psi_l) + delta_l (B(psi_l) - 1)/psi_l, and the bubble puts (1 - B(psi_l))/psi_l applied to H*x_l""",
"""On C: psi_l = p^*psi_l + delta_l, psi_{n+1} delta_l = 0, psi_l delta_l = 0 (hence psi_l^k = (p^*psi_l)^k + delta_l (p^*psi_l)^{k-1}, as in Teleman (2.15)), delta_l delta_l' = 0 for l != l', kappa_j = p^*kappa_j + psi_{n+1}^j, and the only graphs with an edge are the bubbles (0,3) carrying legs l, n+1 (a point of C has at most one node, so graphs with two or more edges miss C, and the other one-edge graphs lie over the boundary of Mbar_{g,n}). The bubble vertex is omega_{0,3} with no translation leg, since p_{m*} of a class of degree >= 2m on Mbar_{0,3+m} vanishes on Mbar_{0,3} for m >= 1. Its leg factors are R^{-1}(0) = 1, since psi vanishes on Mbar_{0,3}. On delta_l the classes kappa_j restrict to p^*kappa_j and a^i(psi_{n+1};H) restricts to f_i^*H, since psi_{n+1} delta_l = 0. So B(psi_l) = B(p^*psi_l) + delta_l (B(p^*psi_l) - 1)/p^*psi_l. After p_*, using p_*(delta_l p^*c) = c, the factor becomes (B(psi_l) - 1)/psi_l with psi_l on M_{g,n}. The bubble puts (1 - B(psi_l))/psi_l applied to H*x_l""")
rep("""with G_i a polynomial in kappa.
""",
"""with G_i = D_i N_i(kappa) p_*[Delta_i^{1/2} a^i(psi_{n+1};H)/b_i(psi_{n+1})], a series in the kappa_j, j >= 1, with coefficients in F_0 and with the number kappa_0 = 2g-2+n among them. G_i contains no psi_l. Here N_i(p^*kappa + psi^j) = N_i(p^*kappa)/b_i(psi) is used.
""")
rep("""*Stability.* By Teleman 2.18, 2.19, monomials in kappa, psi_1, psi_2 are independent in H^*(M_{g,2};F_0) in total degree < (g-1)/2. Take the coefficient""",
"""*Stability.* By Teleman 2.17-2.19, in cohomological degree < (g-1)/2 (deg psi_l = 2, deg kappa_j = 2j), H^*(M_{g,2};Q), the rational cohomology of the open moduli space, is the polynomial ring on psi_1, psi_2 and the pull-backs of the kappa_j of M_g. On M_{g,2}, where delta_{12} = 0, the classes used here, which are the Arbarello-Cornalba classes, satisfy kappa_j = pi^*kappa_j + psi_1^j + psi_2^j (Teleman (2.15)). This is a triangular degree-preserving change of generators, so monomials in kappa_j (j >= 1), psi_1, psi_2 form a basis in that range. Over F_0 they form a basis because H^*(M_{g,2};F_0) = H^*(M_{g,2};Q) (x)_Q F_0. theta acts as 1 (x) theta, killing the monomials. All coefficient extractions below are in this basis, and kappa^0 means no kappa_j with j >= 1. Only independence in some range growing with g is used. Take the coefficient""")
rep("""Since g is arbitrary and M does not depend on g, M_n = 0 for n >= 1. Further""",
"""Since g is arbitrary and M does not depend on g, M_n = 0 for n >= 1 (take g > 4n + 1, so that psi_1^n has degree 2n < (g-1)/2). Further""")
rep("""and Delta_i^g != Delta_k^g because W_i^0/W_k^0 is a non-constant element of Kw, not a root of unity. Hence M = 0""",
"""and Delta_i^g != Delta_k^g because their values at Q = 0 are (W_i^0)^g, (W_k^0)^g, and, with {i,k,m} = {0,1,2}, W_i^0/W_k^0 = -alpha_i (alpha_m - alpha_i)/(alpha_k (alpha_m - alpha_k)) is a non-constant element of Q(alpha), hence transcendental over Q and not a root of unity. Any g >= 1 serves here, since only degree 0 is used. Hence M = 0""")
rep("""*The equation.* Taking eta-adjoints, with B(psi)^dagger = R(-psi) and Psi^dagger = Psi^{-1}, and putting z = -psi gives""",
"""*The equation.* Taking eta-adjoints, with B(psi)^dagger = R(-psi) (this is R^{-1}(z)^dagger = R(-z), the symplectic condition of PPZ (5)), Psi^dagger = Psi^{-1}, (H*)^dagger = H*, and theta eta = 0, so that theta commutes with dagger, and putting z = -psi gives""")
rep("""it is equivalent to z theta_Q S = H* S for S = R Psi e^{U/z}, which is Givental (3).
""",
"""it is equivalent to z theta_Q S = H* S for S = R Psi e^{U/z}, which is Givental (3).

*Remark: the kappa terms carry no further constraint [P; a consistency statement, not a repair; re-derived by hand in round 6d].* With M = 0 and c_i(psi) := Delta_i^{1/2} a^i(psi;H)/b_i(psi) = lambda_i + c_{i,1} psi + ..., one has c_{i,1} = Delta_i^{1/2} Sum_k (W_Psi)_ik Delta_k^{-1/2} = theta log Delta_i^{1/2}, which matches the kappa_0 term of theta log D_i. The flatness equation applied to x = 1 gives -psi theta b_i = c_i b_i - (lambda_i + c_{i,1} psi) b_i, which is theta s^i_j = c_{i,j+1}. So Gtilde_i(kappa) = 0 identically: M = 0 implies A1.5 on M_{g,n} in full. Compare (7.7.c) of Teleman, which was not re-derived in his notation.
""")
rep("""By A1.7 E_t is the unique such data.""",
"""By A1.7, applied over the field ℂ((t)), E_t is the unique such data.""")
rep("""they are two separate statements. In this proof the sign comes from Def 7.1.""",
"""they are two separate statements. In the proof the sign comes directly from PPZ Def 2.2, A1.5 and the adjoint step. Def 7.1 enters only in the comparison remark, where it carries the sign of that route.""")

# ---------- A2.9
rep("""**A2.9 Ambiguity [P, no citation].** X_0 := Psi^{-1} R~^{-1} R Psi satisfies z theta_Q X_0 = [Lambda, X_0]. With""",
"""**A2.9 Ambiguity [P, no citation].** X_0 := Psi^{-1} R~^{-1} R Psi satisfies z theta_Q X_0 = [Lambda, X_0]. Here R~ is the operator with matrix R~_{aj} = eta(H^a, R~ f_j), that is, X_0 is the matrix R~_{aj}^{-1} R_{aj} with R_{aj} = eta(H^a, R(z) f_j) (A2.7). Both matrices solve z theta_Q Y = C Y - Y Lambda, by A2.4 and A2.8. With""")
rep("""So R = R~ diag(exp Sum a_k^{(j)} z^{2k-1}) with a in K~w, and R has coefficients""",
"""So R_{aj} = R~_{aj} exp Sum_k a_k^{(j)} z^{2k-1}, with a in K~w, and R has coefficients""")

# ---------- A3.5
rep("""Prop 2.6 is therefore P mod eq:J=I, the genus-0 QDE and SP1-SP5.""",
"""Prop 2.6 is therefore P mod eq:J=I (source at T_max: OPEN item O1), the genus-0 QDE and SP1-SP5.""")
rep("""The argument below runs over K^ := K~w(i, sqrt(eps1 eps2),""",
"""The argument below runs over K^ := K~w(sqrt(-1), sqrt(eps1 eps2),""")
rep("""In the linear coordinates (s, xi_j, xi_k) the phase splits, so SP1 (linear case) and SP4 factorise Phi into three one-dimensional expansions. Each factor is SP5 with z -> -z,""",
"""In the linear coordinates (s, xi_j, xi_k) the phase splits, so SP1 (linear case) and SP4 factorise Phi into three one-dimensional expansions, with phases a(e^xi - xi), b(e^xi - xi), 3 alpha_i (e^s - s). With q = a e^xi (SP1 in dimension one), each factor is SP5 with beta = 0 and z -> -z,""")

# ---------- B3
rep("""is injective, since 𝒞 (x)_O 𝒞 is free over O and O embeds in Kw((y)).""",
"""is injective, since 𝒞 (x)_O 𝒞 is free over O with basis Y^a (x) Y^b, 0 <= a,b <= 2, O embeds in Kw((y)), and the matrix (Y_i^a Y_j^b) is the Kronecker square of a Vandermonde matrix, invertible because the Y_i are distinct.""")

# ---------- Stage C
rep("""so q_t = T_i - psi exactly.""",
"""so q_t = T^{(i)} - psi exactly, where T^{(i)} is the coefficient of e_i in T(psi), not the operator T_i of the Notation.""")
rep("""**C Step 2. Frame change (C1) [P mod Stages A, B: A3.2 (hence Mumford), B4].**""",
"""**C Step 2. Frame change (C1) [P mod draft Prop 2.10, A3.2 (Mumford as in FP), PPZ Defs 2.2, 2.5, 2.13, Prop 2.14. Stages A, B are needed only to identify Tbar_g[𝖥] with the reduction of 𝒯_g, B6, not for this identity.]**""")
rep("""- h_w contains only z and eps, and log h_w is odd in z. So h_w is t-independent, colour-independent, and h_w(z) h_w(-z) = 1. By B4, pi(P_0^{(i)})""",
"""- h_w contains only z and eps, and log h_w is odd in z. So h_w is t-independent, colour-independent, and h_w(z) h_w(-z) = 1. For the identification with the reduction of 𝒯_g, which is the only use of Stage B in this step: by B4, pi(P_0^{(i)})""")
rep("""at c1 = C1cf (up to the constant). Each T(Gamma,c)""",
"""at c1 = C1cf. Each T(Gamma,c)""")
rep("""GW^cf_g is such an element, up to the constant, and t_c^{2-2g}""",
"""GW^cf_g is such an element, and t_c^{2-2g}""")

# ---------- s.5.1
rep("""11. **Lines 925-927.** The expansions""",
"""11. **Lines 927-929.** The expansions""")
rep("""at the conifold the prefactor is fixed by line 5145.
""",
"""at the conifold the prefactor is fixed by line 5145.
12. **Lines 2633, 2658.** MR2486673 treats only the fibre-rotating C^* (T_F); it does not give eq:J=I at generic alpha. CG07 line 2633 also has the torus on fibres only. The commented pointer Theorems 3.5, 3.7 at 2659 does not match the paper (concave case: Theorem 3.4, Prop 4.1, Cor 4.2). A candidate citation at T_max is MR3414388, Theorem 31 (s.1.2; OPEN item O1 until adopted).
""")

# ---------- s.5.2
rep("""A1.5, A1.7, A1.8, A2.8, A2.9 are rank-free. The field-of-definition paragraph of the conifold Prop 13.2 [repaired in round 6d; the repaired paragraph has had the passes recorded in the round 6 status file]; no statement of this file rests on that paragraph.""",
"""A1.5, A1.7, A1.8, A2.9 are rank-free. A2.8 is rank-free except for the non-degeneracy Delta_i^g != Delta_k^g for some g, which is checked per target (conifold: (F8)). M_0 = 0 also follows with no non-degeneracy from A1.5 at (g,n) = (0,3): the coefficient of a_j b_i c_i in the degree-1 part of (R.omega)_{0,4}(a,b,c,H) is -Delta_i^{1/2}(lambda_i - lambda_j) r_ij, and the same coefficient in theta C_abc is Delta_i^{1/2} (W_Psi)_ji (the computation of the Sign check in A2.8, re-derived by hand in round 6d). The field-of-definition paragraph of the conifold Prop 13.2 was repaired in round 6d (passes on the repaired paragraph: hodge_vertex/lp2_explore/ROUND6_STATUS_2026-11-02.md); no statement of this file rests on that paragraph.""")

# ---------- Appendix A.1
rep("""Controls (sign flip of c_j, dropping the -eps1, -eps2 weights, no h) fail, already at y^1. At g = 3, y^0 the sign flip is invisible, because dim Mbar_3 is even.
""",
"""Controls (sign flip of c_j, dropping the -eps1, -eps2 weights, no h) fail, already at y^1. At g = 3, y^0 the sign flip is invisible, because dim Mbar_3 is even.
- **round6d/classlevel (own code, no project import; SP/lp2_round6/round6d/classlevel/, SP as in the draft's Appendix).** Defect M of A2.8 computed directly from the rows of A2.4-A2.6, with C^c_{1a} solved from the z = infinity expansion of the I-function. M = D M_rat D^{-1} with D = diag((W_j^0)^{1/2}) constant. M_rat = 0 over Q at 6 rational samples through z^5, y^3, and through z^7, y^5 at (1,2,5;3,7). M = 0 in F_p with the square roots as elements of F_p at 5 samples, p about 2^89, and through z^9, y^6 at (1,2,5;3,7), p about 2^127. Controls R^*, R(-z), H* -> -H*, theta -> -theta fail, R^* first at psi^1, y^1. h -> h^{-1} and no h give M = 0, as A2.9 says. Class level: GP localisation with psi classes pulled back from Mbar_{g,n} against the R.omega graph sum with legs, exact over Q: 150 comparisons per weight set at (g,n) = (1,1), (1,2), (2,1), (2,2) through d = 3 at (1,2,5;3,7), (1,3,8;5/2,-3/7), (3,7,-11/2;5/3,-2/7); n = 3 at g = 1, 2 through d = 2 at the first set; g = 3, n = 1 through d = 2 at the first set and d = 1 at the second; no mismatch. Also through d = 2 at two weight sets, no mismatch: lambda_1 insertions; Hodge-vertex form with legs; descendants at k <= 6 from ancestors through S_z of genus-0 GP. Divisor equation with the bubble term on both sides, n = 1 -> 2 through d = 3, n = 2 -> 3 through d = 2, two weight sets, no failure; the same equation without the bubble term fails. Controls that fail: legs at +psi, R^*, R(-z), h^{-1}, no h.
- **round6d/cite (SP/lp2_round6/round6d/cite/scripts/lp_rows_check.py).** Comparison with Lho-Pandharipande, s.1.2: R_1 = (1-L^2)/18, R_2 = (1-24L-2L^2+25L^4)/648 compared exactly with eq:Ifun at H = lambda_0 = 1 through y^8; a control fails.
""")

# ---------- Appendix A.3
rep("""### A.4 Limits and warnings
""",
"""- **round6d/cite (SP/lp2_round6/round6d/cite/scripts/barnes_check.py, barnes_run2.py).** Barnes, parameters (2/3,2/3,2/3,4/3,5/3), trapezoid on Re sigma = -1/2, h = 1/120, 60 digits: LHS - (2 pi/sqrt 3)^3 = 8.1e-47; generic control (0.3,0.45,0.7,0.2): 8.6e-32; wrong e' = 2 fails (47.98 against 67.04). Stdlib decimal, independent of earlier checks [E].

### A.4 Limits and warnings
""")

# ---------- Appendix A.4
rep("""SP/lp2_round6/round6c/refA/rstar.py, rstar2.py, rstar3.py (SP as in the draft's Appendix; scratchpad)""",
"""LP/round6_scripts/round6c/refA/rstar.py, rstar2.py, rstar3.py (LP as in the draft's Appendix)""")
rep("""-406623543620925998705062409/62912004762894336000000. Evidence only.
""",
"""-406623543620925998705062409/62912004762894336000000. Evidence only. R^* differs from R only in class degree >= 2 of Omega. In round 6d, R^* passes every test at (g,n) = (1,1) and every test with total psi power k >= dim Mbar_{g,n} - 1 at (1,2), (2,1), (2,2). Those tests cannot separate R from R^*. Below that power it fails at (1,2) k = 0 (3 of 6), (2,1) k = 0, 1, 2 (2/3, 3/3, 3/3), (2,2) k = 0, 1, 2, 3 (3/6, 8/9, 15/15, 18/18), first at y^1.
""")
rep("""g = 4 with one or two eps pairs and t points; not beyond.""",
"""g = 4 with one or two eps pairs and t points; g = 5 per colouring for graphs with a vertex of genus >= 2 at three eps pairs (round 6d); g = 6 for about 1000 of 5226 such graphs at one eps pair, with an unfinished log (round 6d); single-vertex polar part through g = 11 at (3,7) (round 6d). Nothing beyond these.""")
rep("""The class-level identity Omega = R.omega was not tested beyond integrated numbers.
""",
"""Before round 6d the class-level identity Omega = R.omega was not tested beyond integrated numbers. Round 6d extends the rows at generic alpha to z^9, y^6 (one sample, mod p), and adds the class-level comparison with psi classes of s.A.1. It is still finite-order evidence.
- **Sensitivity to F_cf (round 6d, fresh code; SP/lp2_round6/round6d/backward/, F_p with p = 2^127-1).** Replacing F_cf by F_cf + (1/9) xi^k in the edge graphs changes the sum at valuation k - (2g-3) at g = 2, 3, for k = 0..5. So the polar part sees F_cf through xi^{2g-4} at these genera, namely xi^0 at g = 2 and xi^2 at g = 3. At g = 5, with 356 edge graphs and 9484 coloured terms, xi^6 gives poles (3 colourings) and xi^7 gives none. Mode A gives poles in 2926 of 9484 terms. Evidence only.
""")

io.open(P, "w", encoding="utf-8", newline="").write(s)
print("FRAMEB edits applied:", len(log))
