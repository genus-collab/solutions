# Round 6e: exact-match edits of FRAMEB_PROOF.md. Every old string must occur the stated number of times.
import os
import sys
P = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "FRAMEB_PROOF.md")
s = open(P, encoding='utf-8', newline='').read()
fails = []


def rep(old, new, count=1):
    global s
    c = s.count(old)
    if c != count:
        fails.append((c, count, old[:120]))
        return
    s = s.replace(old, new)


# ---- line 3, 5: provenance (fixer flag; coldB 4; coldF 6)
rep("File: hodge_vertex/lp2_explore/FRAMEB_PROOF.md. Written 2026-11-02 (round 6c).",
    "File: hodge_vertex/lp2_explore/FRAMEB_PROOF.md. Written 2026-11-02 (round 6c); edited in place in rounds 6c, 6d and 6e.")
rep("**This assembled text had no referee pass of its own when it was written; the record of the round-6c passes is the marker ROUND6C-REFEREES in s.0.4 of the draft.**",
    "**The referee passes of this assembled text are recorded at the marker ROUND6C-REFEREES in s.0.4 of the draft. They are AI passes; no human referee has read this proof.**")
rep("The one major item of round 6d (report lit) is not fixed here: it is OPEN item O1 below.",
    "The one major item of round 6d (report lit, citation level) is item O1 below. Its citation was repaired in round 6e, after the verification report verify:ccit, whose corrections to lit's text are applied.")

# ---- line 9: rulings
rep("Stage A uses eq:J=I and the genus-0 QDE at generic alpha on T_max, not on T_F. Ruling Q1 names statements on T_F. This use is outside the literal wording of Q1 and is recorded here pending the user's confirmation.",
    "Stage A uses eq:J=I and the genus-0 QDE at generic alpha on T_max, not on T_F. eq:J=I at T_max is CCIT15 Theorem 31 (s.1.2), a published all-orders genus-0 theorem whose proof uses localisation and the Lagrangian cone, with no Eynard-Orantin TR, remodelling or other forbidden input. BS's own citation MR2486673 covers T_F only. The genus-0 QDE at T_max is derived in s.1.2. Ruling Q1 names statements on T_F. This use of CCIT15 Theorem 31 is outside the literal wording of Q1 and is recorded here pending the user's confirmation.")

# ---- line 11: the O1 paragraph
a = s.find("**OPEN item O1 (round 6d, report lit; citation level; no mathematics of this file changes).**")
b = s.find("\n", a)
if a < 0:
    fails.append((0, 1, 'O1 paragraph'))
else:
    old = s[a:b]
    assert old.endswith("as statements about the decoration sum.")
    new = ("**Item O1 (citation repaired in round 6d-6e; citation level; no mathematics of this file changes).** "
           "The chain needs eq:J=I at generic alpha on T_max. "
           "*What BS cites.* BS 2658 cites MR2486673 = Coates, Commun. Math. Phys. 287 (2009) 1071-1108, arXiv:0810.2200v1. "
           "*What it covers.* Its torus covers the trivial action on P2, so that its Prop 4.1, Cor 4.2 are eq:J=I on T_F only (lambda = E, p = H + E/3). Stage A cannot be moved to T_F, because the Q = 0 product is not semisimple there. "
           "*What is now cited.* The source of record is CCIT15 Theorem 31 (s.1.2), a published genus-0 theorem whose proof uses localisation and the Lagrangian cone, with no Eynard-Orantin TR, remodelling or other forbidden input. "
           "*How it was verified.* Report lit (round 6d) found the defect and proposed the source. An independent reader (report verify:ccit, round 6e, an AI pass) re-read the statement and the definitions, did the substitution, the mirror-map step and the comparison with A1.2 by hand, and ran an exact finite-order check, which is evidence only (s.A.1). The proof of Theorem 31 was read for structure and not re-derived, and the journal print was not read. "
           "*QDE source.* No numbered theorem at T_max was found. The genus-0 QDE at T_max is derived in s.1.2 from Pandharipande's Lemma 1 and A1.4. "
           "The T_max use is outside the literal wording of ruling Q1, which names statements on T_F, and awaits the user's confirmation. "
           "Until it is given, every statement below that rests on eq:J=I carries this one proviso, called the Q1 proviso of item O1: A2.1, A2.4-A2.6, A2.9, A3, parts (A2), (A3) and the row form of the graph sum of Theorem FB_A; Theorem FB_B except (b); Proposition 13.2_P2; Corollaries G and G_{k,h}; the discharge paragraph at the end of s.2. "
           "Not affected: A1, A2.3, A2.8, A3.2, the identity Tbar_g[𝖥] = G_g[𝖥] of Theorem FB_C (C1), (C2), Lemma (t_cf), and draft Theorems LC_P2, R_P2, Corollary SV_P2, Theorem VI as statements about the decoration sum. "
           "This partition is a dependency reading of the current text (round-6d fix step); no referee has checked it.")
    s = s[:a] + new + s[b:]

# ---- line 49: Theorem FB_A label
rep("[P mod BS eq:J=I and the genus-0 QDE (eq:SHb is derived from them in A2.1), with the source of eq:J=I at generic alpha OPEN (item O1; candidate CCIT15 Thm 31, s.1.2), so that (A2), (A3) and the row form of the graph sum are conditional on O1; Teleman",
    "[P mod eq:J=I (CCIT15 Thm 31, s.1.2; T_max use under Q1 pending the user, item O1) and the genus-0 QDE (derived in s.1.2; eq:SHb is derived from them in A2.1), on which (A2), (A3) and the row form of the graph sum rest; Teleman")
rep("Two referee passes on the report form. A2.8 was rewritten in round 6c",
    "Two referee passes on the report form; passes on the assembled text as recorded at the marker ROUND6C-REFEREES in s.0.4 of the draft. A2.8 was rewritten in round 6c")

# ---- line 66: FB_B
rep("One referee pass, conditional on Stage A. Conditional on OPEN item O1 through Theorem FB_A, except (b).]",
    "One referee pass in report form, conditional on Stage A; later passes as recorded at the marker ROUND6C-REFEREES. Through Theorem FB_A, also P mod eq:J=I (CCIT15 Thm 31, s.1.2), except (b); subject to the Q1 proviso of item O1.]")

# ---- line 112: 13.2_P2
rep("Conditional on OPEN item O1 through Stages A, B.]",
    "Through Stages A, B, also P mod eq:J=I (CCIT15 Thm 31, s.1.2); subject to the Q1 proviso of item O1.]")

# ---- line 116: G
rep("Conditional on OPEN item O1 through Proposition 13.2_P2.]",
    "Through Proposition 13.2_P2, also P mod eq:J=I (CCIT15 Thm 31, s.1.2); subject to the Q1 proviso of item O1.]")

# ---- line 122: G_{k,h}
rep("[Same status, including OPEN item O1, plus the polynomiality of Theorem FB_B(b).]",
    "[Same status, including eq:J=I (CCIT15 Thm 31) and the Q1 proviso of item O1, plus the polynomiality of Theorem FB_B(b).]")

# ---- line 141
rep("and of eq:J=I, which stay cited.", "and of eq:J=I (CCIT15 Thm 31), which stay cited.")

# ---- line 145
rep("End-to-end inputs of Corollary G: eq:J=I (at generic alpha on T_max: OPEN item O1; candidate CCIT15 Thm 31, s.1.2) and the genus-0 QDE (eq:SHb is derived from them in A2.1);",
    "End-to-end inputs of Corollary G: eq:J=I (CCIT15 Thm 31, s.1.2; T_max use under Q1 pending the user) and the genus-0 QDE (derived in s.1.2 from Pandharipande math/9806133v1 Lemma 1 and A1.4; eq:SHb is derived from them in A2.1);")

# ---- line 147 (coldB 4)
rep("The chain as assembled here has not.",
    "The chain as assembled here had the passes recorded at the marker ROUND6C-REFEREES in s.0.4 of the draft.")

# ---- line 159
rep("**OPEN item O1:** at generic alpha this citation does not supply the input; the candidate source CCIT15 Theorem 31 (s.1.2) has not been adopted. It is not re-proved.",
    "At generic alpha the source of record is CCIT15 Theorem 31, s.1.2 (item O1; T_max use under Q1 pending the user), and the genus-0 QDE at T_max is derived in s.1.2. eq:J=I is not re-proved.")

# ---- line 168: the literature bullet, and the new QDE bullet
a = s.find("- **eq:J=I at generic alpha on T_max: candidate source (OPEN item O1; not adopted).**")
b = s.find("\n", a)
if a < 0:
    fails.append((0, 1, 's.1.2 bullet'))
else:
    bullet1 = (
        "- **eq:J=I at generic alpha on T_max: source of record (item O1; source settled in round 6e; the T_max use under Q1 is pending the user's confirmation).** "
        "Coates-Corti-Iritani-Tseng (CCIT15), A mirror theorem for toric stacks, Compos. Math. 151 (2015) 1878-1912 (MR3414388), arXiv:1310.4163v5, Theorem 31 with Definition 28, eq. (19), eq. (4) and Remark 3. Read at arXiv v5 in rounds 6d and 6e. The journal print was not read, so the number Theorem 31 in print rests on a search snippet of the Imperial copy. No erratum was found on the Cambridge page or in one search; this does not exclude one. The structure of the proof (Theorem 41, Props 45-47) was read; the proof was not re-derived. "
        "*Scope.* Z = K_P2 x C^2 is the smooth toric variety with N = Z^5, rays b_0 = (1,0,1,0,0), b_1 = (0,1,1,0,0), b_2 = (-1,-1,1,0,0), b_3 = (0,0,1,0,0), b_4 = (0,0,0,1,0), b_5 = (0,0,0,0,1), and maximal cones spanned by b_3, b_4, b_5 and two of b_0, b_1, b_2. These cones are unimodular, so Box = {0}, and L = Z(1,1,1,-3,0,0). Z is semi-projective in the sense of CCIT15 s.3.1: it has three fixed points, and Z -> Spec H^0(O_Z) = C^3/mu_3 x C^2 is projective, with -D_3 = 3H relatively ample. In fan terms the support Cone(b_0,b_1,b_2) x R_{>=0}^2 is convex and full-dimensional, and the piecewise linear function with value -1 at b_3 and 0 at the other rays is the maximum of its linear pieces, strictly. NE(Z) is generated by the line. "
        "*Torus.* CCIT15's torus is N x C^* = (C^*)^6/C^*. BS's T_max (weights (alpha_0,alpha_1,alpha_2,0,-eps1,-eps2) on C^6, BS line 2349) maps onto it with kernel mu_3. The character lattices differ by index 3, so S_T = C(alpha_0,alpha_1,alpha_2,eps1,eps2) = Q_{T_max}, and the localisation formulas agree. "
        "*Dictionary.* For the T_max action the toric divisor classes are u_i = H + alpha_i (i = 0,1,2), u_3 = -3H, u_4 = -eps1, u_5 = -eps2, with H = c_1 of O(1) with its BS lift, H|P_i = -alpha_i. So H = -u_3/3, alpha_i = u_i - H, eps_j = -u_{3+j}. Then u_j|P_i are BS's tangent weights w^{(j)}_i (BS line 2383), and the relation u_0 u_1 u_2 = 0 is BS's prod(H + alpha_i) = 0. "
        "*Substitution.* Take S empty and all t_i = 0 in (19). Then lambda = d(1,1,1,-3,0,0), d >= 0, and the factors are 1/prod_{m=1}^d (H+alpha_i+mz) for i = 0,1,2, prod_{m=0}^{3d-1}(-3H-mz) for i = 3, and 1 for i = 4,5. So I_CCIT(Q,z) = z Sum_d Q^d I_d(z), where z y^{H/z} Sum_d y^d I_d(z) is eq:Ifun. "
        "*Mirror map.* Each summand with d >= 1 is -3H times a function homogeneous of degree -1 in (H,alpha,z), with expansion at z = infinity beginning 3H(-1)^d (3d-1)!/(d!)^3 z^{-1}. Reducing H^3 by the relation keeps the power of z. Hence I_CCIT(Q,-z) = -z + g(Q)H + O(1/z), with g(Q) = 3 Sum_{d>=1} (3d-1)!/(d!)^3 (-Q)^d independent of alpha and g(0) = 0; g = I_0^{[1]} - log y in the notation of eq:Ian. Theorem 31 says this is a point (4). The H_+ part of (4) is -z + t(z), so t(z) = g(Q)H, and by Remark 3 I_CCIT(Q,z) = J(g(Q)H,z). "
        "*Divisor equation.* J(tau H,z) = z e^{tau H/z}(1 + Sum_d (Q e^{tau})^d <<phi_a/(z(z-psi))>>_{0,1,d} phi^a), including the d = 0 classical terms z e^{tau H/z}. Hence Sum_d Q^d I_d = e^{g(Q)H/z}(1 + Sum_d (Q e^{g(Q)})^d <phi_a/(z(z-psi))>_{0,1,d} phi^a). Rename Q as y and multiply by z y^{H/z}. This is eq:J=I with Q = y e^{g(y)} = exp I_0^{[1]}, by eq:Jred. "
        "*Comparison with A1.2.* CCIT15 s.2.2 defines the invariants on Z by virtual localisation over the proper fixed locus (via Liu, arXiv:1107.4712v3, eq. (51)), which is BS's Notation at line 1233 for this Z and T. A1.1 and A1.2 identify this with the twisted integral on Mbar_{0,n}(P2,d) in genus 0. For d = 0, Mbar_{0,n}(Z,0) = Mbar_{0,n} x Z, and both sides are Sum_i prod a_j|P_i /W_i^0 times a psi integral on Mbar_{0,n}; the pairing is BS's eta. On genus-0 trees R pi_* f^* O_w = O_w, so the C^2 factor multiplies each invariant by 1/(eps1 eps2) and eta^{ab} by eps1 eps2, and J is unchanged. "
        "*Inputs of CCIT15's proof.* Brown's recursive characterisation of the cone in the form of Theorem 41; Liu's localisation formula; the Coates-Givental/CCIT hypergeometric modification with Tseng's quantum Riemann-Roch at the fixed points (here BG_sigma is a point, so that theory is on Mbar_{0,n}); and Theorem 2, whose TRR is the genus-0 psi-class relation. Eynard-Orantin topological recursion, remodelling and mirror curves do not occur; a text search finds no such term, and [FLT12] is cited only as an application upstream of nothing used in the proof. "
        "*Second source, read only at the level of the statements, not adopted.* Ciocan-Fontanine-Kim, arXiv:1304.7056v3, Theorem 1.2.2 and Corollaries 1.4.1, 1.4.2 (targets W//G with isolated torus-fixed points, semi-positive case). "
        "*Other sources read in round 6d and found not to cover T_max.* MR2486673 = Coates, arXiv:0810.2200v1, Prop 4.1 and Cor 4.2, prove eq:J=I on T_F only (lambda = E, p = H + E/3), since its torus covers the trivial action on P2. Coates-Givental (Ann. of Math. 165 (2007) 15-53, math/0110142) Theorem 2 and CCIT (Duke Math. J. 147 (2009) 377-438, math/0702234) Theorem 4.6 have T acting on fibres only. Givental alg-geom/9701016 Theorem 0.2 has compact X and non-negative bundles. Givental alg-geom/9603021 Theorem 9.5 gives the T-equivariant small J-function of P^n only. "
        "*Evidence [E], not a proof step.* s.A.1 records the exact check of eq:J=I through y^4 at generic alpha.")
    bullet2 = (
        "- **Genus-0 QDE at T_max [P mod the following two statements].** No numbered theorem at T_max was found: Iritani arXiv:0903.1463v3 Prop 2.4 assumes X proper, and Pandharipande math/9806133v1 s.1 assumes X convex. The derivation is short. "
        "Input (a): Pandharipande, Rational curves on hypersurfaces (after A. Givental), math/9806133v1, Lemma 1: for n >= 3, psi_1 on Mbar_{0,n}(P2,d) is the sum of the boundary divisors separating marking 1 from markings 2, 3 (P2 is convex; opened in round 6e). "
        "Input (b): the boundary structure with the normalisation sequence, as in A1.4 (restriction of e^{-1}(R pi_* f^* N) to a boundary divisor is the product of the two factors times e(ev^* N), and e(N)[Diagonal of P2] = Sum eta^{ef} phi_e (x) phi_f). "
        "This gives the genus-0 topological recursion relation (TRR) of the theory of Z, namely <<c, a, psi^{k+1} b>> = Sum_{e,f} <<psi^k b, e>> eta^{ef} <<f, c, a>>, with only stable splittings on the right. This is the psi-class relation on Mbar_{0,n}; it has nothing to do with Eynard-Orantin TR. "
        "Put V(a,b) := eta(a,b) + <<a, b/(z-psi)>>, so that S(a) = Sum H^c eta^{cb} V(a,H^b) by eq:fundsol. Then z d_c V(a,b) = <<c,a,b>> + Sum_{k>=0} z^{-k-1} <<c,a,psi^{k+1} b>> = Sum_{e,f} V(e,b) eta^{ef} <<f,c,a>> = V(c*a, b). So z d_c S(a) = S(c*a). Restrict to t_0 = t_2 = 0, t_1 = log Q, where d_1 = theta_Q. WDVV is not used; commutativity of the operators c* follows by comparing z d_c z d_e with z d_e z d_c, since S(a) = a + O(1/z) is invertible. "
        "Same form of (6) and of the fundamental solution: Pandharipande s.1.2 (6), s.1.3 Prop 2. Consistent with CCIT15 Theorem 2 with (6) and Givental math/0305409v1 Theorem 1 (TRR, string and dilaton equivalent to the cone property); these were read and are not used as inputs. Evidence [E]: s.A.1.")
    s = s[:a] + bullet1 + "\n" + bullet2 + s[b:]

# ---- line 188: SP row (consistency with the marker)
rep("No pass in the present text at the time of writing. |",
    "No pass in the present text in round 6b. In round 6c, refS2 re-derived the calculus as used by FB:A3.5, B2, B4 (marker ROUND6C-REFEREES in s.0.4 of the draft). |")

# ---- line 193
rep("P mod eq:J=I (source at generic alpha: OPEN item O1), genus-0 QDE as a whole theorem;",
    "P mod eq:J=I (CCIT15 Thm 31, s.1.2), genus-0 QDE (derived in s.1.2) as a whole theorem;")

# ---- line 213: s.1.4 last bullet
rep("- The assembled text.\n",
    "- The assembled text. (Later passes: marker ROUND6C-REFEREES in s.0.4 of the draft.)\n")

# ---- A1.7, A3.3: H^* as in leaf (c)
rep("p^* is injective because p_*(psi_1 p^* x) = (2k-2) x.",
    "p^*: H^*(Mbar_k;Q) -> H^*(Mbar_{k,1};Q) is injective because p_*(psi_1 p^* x) = (2k-2) x.")
rep("p^*: A^*(Mbar_k) -> A^*(Mbar_{k,1}) is injective",
    "p^*: H^*(Mbar_k;Q) -> H^*(Mbar_{k,1};Q) is injective")

# ---- line 253: A2.1
rep("**A2.1 Rows are formal [P mod eq:J=I and the genus-0 QDE; source of eq:J=I at generic alpha: OPEN item O1].**",
    "**A2.1 Rows are formal [P mod eq:J=I (CCIT15 Thm 31, s.1.2; T_max use under Q1 pending the user) and the genus-0 QDE (derived in s.1.2)].**")

# ---- A2.8 status (ref28hand 9, with the record)
rep("(the proposer's and the round-6d referee's, with PPZ Defs 2.2, 2.13 and Teleman 2.15-2.19, 7.1, 7.6-7.8, 8.2 compared against the arXiv text). The round-6d reads refNew, coldF and coldB also re-derived it by hand and found only the notation corrected here.",
    "(the proposer's and that of the round-6d referee ref28hand, with PPZ Defs 2.2, 2.13 and Teleman 2.15-2.19, 7.1, 7.6-7.8, 8.2 compared against the arXiv text). The round-6d reads refNew, coldF and coldB also re-derived it by hand and found only the notation corrected here. The record of these AI passes is the marker ROUND6C-REFEREES in s.0.4 of the draft.")

# ---- line 382
rep("Prop 2.6 is therefore P mod eq:J=I (source at T_max: OPEN item O1), the genus-0 QDE and SP1-SP5.",
    "Prop 2.6 is therefore P mod eq:J=I (CCIT15 Thm 31, s.1.2), the genus-0 QDE (derived in s.1.2) and SP1-SP5.")

# ---- s.5.1 item 12
rep("A candidate citation at T_max is MR3414388, Theorem 31 (s.1.2; OPEN item O1 until adopted).",
    "The correct citation at T_max is MR3414388, Theorem 31 (s.1.2).")

# ---- s.5.2 last paragraph (coldB 5, merged with the status-file pointer)
rep("The field-of-definition paragraph of the conifold Prop 13.2 was repaired in round 6d (passes on the repaired paragraph: hodge_vertex/lp2_explore/ROUND6_STATUS_2026-11-02.md); no statement of this file rests on that paragraph.",
    "The field-of-definition paragraph of the conifold Prop 13.2 was repaired in round 6d, after the major item of refCon in round6_scripts/round6c_result.json (passes on the repaired paragraph: hodge_vertex/lp2_explore/ROUND6_STATUS_2026-11-02.md). No statement of this file rests on that paragraph.")

# ---- evidence pointers
rep("(own code, no project import; SP/lp2_round6/round6d/classlevel/, SP as in the draft's Appendix)",
    "(own code, no project import; round6_scripts/round6d/classlevel/)")
rep("**round6d/cite (SP/lp2_round6/round6d/cite/scripts/lp_rows_check.py).**",
    "**round6d/cite (round6_scripts/round6d/cite/scripts/lp_rows_check.py).**")
rep("**round6d/cite (SP/lp2_round6/round6d/cite/scripts/barnes_check.py, barnes_run2.py).**",
    "**round6d/cite (round6_scripts/round6d/cite/scripts/barnes_check.py, barnes_run2.py).**")
rep("(round 6d, fresh code; SP/lp2_round6/round6d/backward/, F_p with p = 2^127-1)",
    "(round 6d, fresh code; round6_scripts/round6d/backward/, F_p with p = 2^127-1)")

# ---- s.A.1: new item, after the Lho-Pandharipande comparison
anchor = "compared exactly with eq:Ifun at H = lambda_0 = 1 through y^8; a control fails.\n"
new_item = ("- **round6e/citeO1 (round6_scripts/round6e/citeO1/scripts/jI_check.py, output out_d4.txt).** Round 6e, exact fractions, finite order, evidence only. "
            "Direct graph sum of virtual localisation on Mbar_{0,n}(Z,d), n = 1, 2, d <= 4, independent of the I-function, at three weight points: alpha = (2/7,5/11,-3/13), (1/3,-4/5,11/2), and (1/7,3/5,11/13) with e_1(alpha) = eps1 + eps2. "
            "Results at all three: eq:J=I at P_0, P_1, P_2 through y^4 at a rational z; eq:Ian 1/z coefficients; eq:etaflat; S(1) = J/z and S(H) = theta_Q J from two-point sums; the QDE z theta_Q S(H^a) = Sum C^c_{1a} S(H^c), a = 0,1,2, through Q^4, with C from three-point invariants. "
            "Control: dropping the mirror map fails. GW_{0,d} eps1 eps2 = 3, -45/8, 244/9, -12333/64 at every point, and <H^2>_{0,1,d} eps1 eps2 / e_1(alpha) = -5, 87/4, -1544/9, 27671/16.\n")
rep(anchor, anchor + new_item)

if fails:
    for f in fails:
        print('FAIL', f)
    sys.exit(1)
open(P, 'w', encoding='utf-8', newline='').write(s)
print('ok, chars', len(s))
