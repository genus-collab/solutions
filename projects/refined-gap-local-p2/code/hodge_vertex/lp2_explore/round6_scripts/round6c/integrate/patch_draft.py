# Integrates FRAMEB_PROOF.md into LP2_REGULARITY_PROOF_draft.md. Exact replacements; aborts with no write if any target is not unique.
import os
import sys
P = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "LP2_REGULARITY_PROOF_draft.md")
s = open(P, encoding="utf-8", newline="").read()
assert "\r\n" not in s
errors = []

def sub(old, new):
    global s
    n = s.count(old)
    if n != 1:
        errors.append((n, old[:90]))
        return
    s = s.replace(old, new)

def line(prefix, new):
    """Replace the unique line starting with prefix."""
    global s
    lines = s.split("\n")
    hits = [i for i, l in enumerate(lines) if l.startswith(prefix)]
    if len(hits) != 1:
        errors.append((len(hits), "LINE " + prefix[:90]))
        return
    lines[hits[0]] = new
    s = "\n".join(lines)

# ---------------------------------------------------------------- header, labels
sub("Not refereed as a whole: see s.0.4. -->",
    "Not refereed as a whole: see s.0.4. Round 6c (2026-11-02): edited in place from LP2_REGULARITY_PROOF_draft_round6a.md, integrating FRAMEB_PROOF.md (Theorems FB_A, FB_B, FB_C, Proposition 13.2_P2, Lemma (t_cf), Corollaries G and G_{k,h}; report journal: round6_scripts/round6b_result.json). Corollary G now carries its status as stated in s.0.2. FRAMEB_PROOF.md as assembled and these edits are not refereed at the time of this edit: see s.0.4. -->")
sub("the round-5 reports modrel, VI, ThmL, MVmatch are in journal wf_9d7fd94e-8d7).",
    "the round-5 reports modrel, VI, ThmL, MVmatch are in journal wf_9d7fd94e-8d7; FB:x = FRAMEB_PROOF.md, same folder, round 6c, item x).")

# ---------------------------------------------------------------- s.0.2
line("**Corollary G (refined gap, BS rmk:cfgap, line 5422).**",
"""**Proposition 13.2_P2 (frame-B formula on T_F, d >= 0 reading).** [P mod Stages A, B of FRAMEB_PROOF.md and Propositions 2.10, 2.12; inputs in s.8.2a]

Theorems FB_A, FB_B, FB_C are stated in FRAMEB_PROOF.md s.0.2. For g >= 2 put F_g^{T_F}(y) := Sum_{d>=0} 𝒬(y)^d GW^{T_F}_{g,d}, in Q(ε1, ε2)[[y]], with 𝒬(y) the mirror map (BS's cQ) and F := θ log A. By Theorem FB_B(d), (e), F_g^{T_F} is the point-frame graph sum, and it lies in Q(ε1, ε2)[y^{±1}, Δ_c^{-1}][F], with deg_F <= 3g-3 and a unique polynomial in F. With Cf as in s.1.3,

    Cf(F_g^{T_F}) = G_g of (0.1) at c1 = C1cf.

With BS's printed d > 0 convention (lines 1277, 2446), Cf(GW_g) = G_g - GW^{T_F}_{g,0}, where GW^{T_F}_{g,0} in Q(ε1, ε2) is the constant of Theorem FB_B(b) (E^2/640 at g = 2). BS's eq:F2ref includes the d = 0 term (s.9 item 17).

**Corollary G (refined gap, BS rmk:cfgap, line 5422).** [P mod Proposition 13.2_P2, Lemma 1.2(c), Theorem R_P2, Corollary SV_P2 and Theorem VI. End to end: P mod (MV) with its textbook facts, and the inputs of FRAMEB_PROOF.md s.1, listed in s.8.2a. Proof in s.7.6.]

For g >= 2, with GW^cf_g := Cf(F_g^{T_F}),

    GW^cf_g = (-3)^{g-1} β_{2g-2}(ε1, ε2)/t_cf^{2g-2} + O(1).

O(1) means: the difference is holomorphic in a neighbourhood of t_cf = 0, equivalently of y = -1/27, with no negative powers of t_cf, and its Taylor coefficients lie in Q(ε1, ε2). The value of the O(1) term is not determined. The statement holds with either convention for d = 0, since the constant is O(1).

**Corollary G_{k,h}.** [Same status, plus the polynomiality of Theorem FB_B(b). Proof in s.7.6.] Write GW_g = Sum_{k+h=g} E^{2k} (ε1 ε2)^{h-1} F_{k,h}. For k + h >= 2,

    F^cf_{k,h} = 𝔟_{k,h} t_cf^{2-2h-2k} + O(1),   𝔟_{k,h} = (-3)^{k+h-1} [E^{2k} (ε1 ε2)^{h-1}] β_{2k+2h-2}.

This is (Conifold asymptotics) of BS def:dirint. In particular F^cf_{0,h} and F^cf_{h,0} have the leading terms of thm:cfasym with remainder O(1).

**Referee status of these three statements.** They are not labelled refereed beyond the passes that happened. Each link had at least one referee pass in report form (s.0.4). FRAMEB_PROOF.md as assembled has had no referee pass of its own, and neither has the chain of Corollary G as assembled. FRAMEB_PROOF.md s.1.4 lists what the passes did not cover.""")
sub("If in addition G_g = GW^cf_g (Proposition 13.2_P2, open), then GW^cf_g = (-3)^{g-1} β_{2g-2}/t_cf^{2g-2} + O(1). The identification t_c = t_cf is Lemma 1.2(c). This would prove BS's conjecture rmk:cfgap. **This corollary is NOT claimed** while Proposition 13.2_P2 is open.\n\n", "")

# ---------------------------------------------------------------- s.0.3
line("- **BS genus-0 inputs (allowed).**",
"""- **BS genus-0 inputs (allowed).** Since round 5 these are **not** inputs of Theorems LC_P2 and R_P2, or of Corollary SV_P2 on the recorded route for (VI): Proposition 2.10 is proved from Theorem A alone, and s.3-s.6 contain no other BS link (dependency reading by the round-5 modrel author and referee). (η^{ab} and L_pf enter those theorems as definitions.) They enter only Proposition 13.2_P2, Corollaries G, G_{k,h}, and the interface statements Proposition 2.6, Theorem 2.7, Proposition 2.11. The labels (C1)-(C5) of R2:s.0 and later rounds are kept for traceability, with their present role as follows; the parenthesised labels (Cn) are never to be read as C1cf or C1i.""")
line("  - (C1) eq:Ifun 2621,",
"""  - (C1) reduced to eq:J=I 2663 and the genus-0 QDE, with eq:SHb 3305 and the definition eq:Ifun 2621. That e^{u_i/z} R^α_{0i} formally solves eq:pf (2628) is derived in FB:A2.1-A2.5, with no analytic asymptotics; eq:asymp 2943 and eq:canR 3566 are no longer cited for it. eq:J=I is cited by BS to MR2486673, which was not opened; it is not re-proved.""")
line("  - (C2) (Normalisation) 3497,",
"""  - (C2) discharged. (Normalisation) 3497 is a theorem, FB:A3.3, modulo the literature below; Psi_{0i} = Delta_i^{-1/2} and Delta_i = 3 ε1 ε2 W_i are FB:A2.5; the weights at about 2383 are FB:A1.2. eq:phi 3031 is not used.""")
line("  - (C3) eq:RHb 3573,",
"""  - (C3) discharged. (Symplecticity) 3496 is FB:A2.6, and eq:RHb 3573 is derived there from the QDE. eq:ADEFXdef 2681 and eq:etaflat 2419 are definitions. lem:modrel 2695 is no longer cited: both of its relations are proved in Proposition 2.12.""")
line("  - Also used: (C4) prop:algindABC 2854,",
"""  - (C4) prop:algindABC 2854, with eq:XFABC 2860 (primary form BFGW Lem 4.2, s.9 item 20), only for statements involving A: ev_A on the ring with generator A (FRAMEB_PROOF.md, Definition (Cf)), Cf(R_{ai}), and the legged form of eq:gtcf. It is not used in Proposition 13.2_P2 or Corollaries G, G_{k,h}: there A has cancelled, and only the transcendence of F over C(y) is needed, which is proved by a monodromy argument (FB:s.3, proof of (e)). Its former use in Proposition 2.11 is deleted (s.2.6).
  - t_cf at line 5145 with eq:pfinf 2757 is now only a definition; its value is Lemma 1.2(c).""")
line("  - For Proposition 2.11 and Proposition 13.2_P2 only: (C5)",
"""  - (C5) no longer used: thm:haegen 3902, lem:lho 3247, eq:gtH 3799 at rational tori, and the regularity on the diagonal at about line 4056. Stages A, B of FRAMEB_PROOF.md replace them; BS 4056 is proved in a stronger form in Theorem FB_B(b), FB:B5 [P mod Graber-Pandharipande].
- **Literature for Proposition 13.2_P2 and Corollaries G, G_{k,h}.** FRAMEB_PROOF.md s.1.2 says what was opened and by whom.
  - Teleman, arXiv:0712.0160: Theorem (2), Prop 7.7, Prop 8.7(ii) (FB:A1.8, A2.8). The proof of Prop 7.7 was read, not re-derived; it is the one non-classical citation on which the graph sum rests.
  - Pandharipande-Pixton-Zvonkine (PPZ), arXiv:1303.1043: Defs 2.2, 2.5, 2.13, Prop 2.14 (FB:A2.7, A3.2, C Steps 1, 2); Props 2.8, 2.9 only for the reconciliation remark (FB:C Step 9).
  - Graber-Pandharipande (GP), arXiv:alg-geom/9708001: formula (1), Apps. A, C (FB:A1.2, B5).
  - Mumford's theorem in the form of Faber-Pandharipande (FP), math/9810173, s.1.1 (FB:A3.2; the original paper was not opened), and FP Theorem 4 (FB:A3.3).
  - Barnes's second lemma, Quart. J. Math. 41 (1910) 136-140 (FB:C Step 6): the primary source was not opened; the form used was confirmed through a secondary statement. DLMF 15.6.7 and Gauss's values of ψ at 1/3, 2/3 (same step).""")

# ---------------------------------------------------------------- s.0.4
line("- The identification of G_g with a GW potential. This is Proposition 13.2_P2, s.8.3.",
"""- Anything off T_F. All statements are on the fibrewise torus, per colouring, with no extension to general CY-torus polar parts.
- The value of the regular part of GW^cf_g.
- BS's F_{k,h} in Q[X^{±1}, F]. Only F_g^{T_F} in Q(ε1, ε2)[y^{±1}, Δ_c^{-1}][F] is proved (Theorem FB_B(e)). Regularity at y = 0 is established for the total, and BS 4071-4076 argues from the total, not coefficient by coefficient. Stored data at genus 2-4 have BS's form [E]. Open for all g, independent of Corollary G.
- (Orbifold regularity) of BS def:dirint.
- The holomorphic anomaly equations.
- Uniqueness in BS def:dirint.""")
line("- Corollary G.", "@@DELETE@@")
s = s.replace("@@DELETE@@\n", "")
line("- A literature-free proof of Lemma A or of Theorem VI:",
"- A literature-free proof of Lemma A, of Theorem VI, of Proposition 13.2_P2 or of Corollary G: they are proved modulo the literature listed in s.8.2 and s.8.2a. (MV), Teleman's Prop 7.7 and eq:J=I stay cited, and their proofs were not re-derived.")
line("- That the whole text has been refereed.",
"""- That the whole text has been refereed. Round 6 gave each assembled part one referee pass of its own, each with verdict sound with minor fixes and no major item; the fixes are applied here. Coverage: s.0-s.2 apart from the items below, including the extension of Fpot_h in s.1.5 and (H0) of s.7.1 (D_rows); s.3-s.6 (D_core); s.7.0-s.7.7, with (H0) re-derived again (D_hodge); s.7.8 (D_VI); LEMMA_A_PROOF.md (L_logic). Proposition 2.10, Proposition 2.12 and Theorem VI also had a pass in round 5. The frame-B reports of round 6b had passes in report form: Stage A (refFB_A, refFB_A2), Stage B (refFB_B, conditional on Stage A), Stage C (refFB_C, and refFB_C2, conditional on Stages A, B), all sound with minor fixes and no major item. **Not refereed:** FRAMEB_PROOF.md as assembled, and the chain of Corollary G as assembled; the round-6c edits of s.0, s.1.2-s.1.4, s.2.4, s.2.6, s.2.7, s.7.5, s.7.6, s.8, s.9 and the Appendix of this text; the items of FRAMEB_PROOF.md s.1.4, in particular SP1-SP4 line by line and the proof of Lemma 2.4 (s.2.2-2.5 as a whole have not had a pass in their present form; the part of Theorem 2.7(e) that the proof of record uses is FB:B4, which the referees of the reports re-derived by hand); Remark 2.13(a)-(c); Lemma 1.4 at genus 4. D_rows and D_hodge judged by reading that s.2.2-2.5 and s.2.7 are not used in the proofs of Theorems LC_P2, R_P2 and Corollary SV_P2. All finite-order computations of those passes are evidence only.
- ROUND6C-REFEREES: placeholder for the main session. The round-6c referee passes on FRAMEB_PROOF.md and on this draft had not run when this line was written, and none is claimed here.""")

# ---------------------------------------------------------------- s.1.2
line("**Lemma 1.2.** [P]",
"**Lemma 1.2.** [(a), (b): P. (c): P mod DLMF 15.6.7, Barnes's second lemma and Gauss's values of ψ at 1/3, 2/3; this is FRAMEB_PROOF.md Lemma (t_cf)]")
line("- (c) t_c = t_cf and C1cf = A_cf, with t_cf, A_cf as in BS sec:cfgen, line 5145.",
"- (c) On y in (-1/27, 0), with the expression at BS 5145 taken on the branch log y = log(-y) + iπ: t_c = t_cf, C1cf = A_cf and XX|_{c1=C1cf} = F_cf, with t_cf, A_cf, F_cf as in BS sec:cfgen, line 5145.")
line("- (c) t_cf is a combination of the periods of eq:pfinf (2757),",
"- (c) This is FB:C Step 6. Put l := log(-y), u := -27y, g := I_0^{[1]} - log y, and let g_2 be the series part of I_1^{[2]}. On log y = l + iπ the expression at BS 5145 reads t_cf = (sqrt3/(2π))(l^2/2 + g l + g_2 - π^2/2). DLMF 15.6.7 and Gauss's values give the logarithmic form wc = -(sqrt3/(2π))(A l + h), with h as in Proposition 2.12 Step 0, and eq:Ian gives θ t_cf = (sqrt3/(2π))(A l + h). So A_cf = -wc = C1cf and F_cf = XX. Thus t_cf and t_c have the same θ, and it remains to compare constants. The Barnes representation of t_c, with Barnes's second lemma for its value at u = 1 and the triple pole at s = 0 for u -> 0, gives t_c = (sqrt3/(2π))(l^2/2 - π^2/2) + o(1) as u -> 0, the same as t_cf, since g, g_2 -> 0. So t_cf = t_c, and the expansion printed at BS 5145 follows from (a). On the other branch, log y = l - iπ, direct substitution shows that the printed expression equals t_cf - sqrt3 π - i sqrt3 (l + g). Since l, g are real on (-1/27, 0), its real part at Δ_c = 0 is -sqrt3 π != 0, and it does not vanish there. ∎")
sub("Remark: (c) was new in round 4. Critics 1 and 2 of round 4 read eq:pfinf (it is L2b∘θ) and BS 5145 and agree with the proof. BS asserts the value of the combination at 5145 without proof; round 6 (D_rows) tested it only by float [E]. It is used only for Corollary G. The normalisation",
    "Remark: the round-4 proof of (c) relied on the expansion printed at BS 5145, which BS assert without proof and which round 6 (D_rows) tested only by float [E]. The proof above replaces that reliance; in report form it had the passes refFB_C, refFB_C2. (c) is used only for Corollary G. The normalisation")

# ---------------------------------------------------------------- s.1.3
sub("With c1 = C1cf = A_cf they are Cf(R_{ai}). Here Cf is BS's ring homomorphism G_{T_F}[A,F] -> G_{T_F}[A_cf,F_cf] induced by Cf(A) = A_cf, Cf(F) = F_cf (BS lines 5148-5160; well defined by (C4)), so it is the identity on G_{T_F}, and D_cf := X/A_cf^2 treats X as fixed. BS's further gloss Cf(f)(τ) = f(-1/(3τ)) cannot be literal on G_{T_F}, since it does not fix the Hauptmodul. A general solution",
    "With c1 = C1cf = A_cf they are Cf(R_{ai}), for Cf as follows (this statement involves A, so it is P mod (C4)). **Cf is BS's definition read through generators (FRAMEB_PROOF.md, Definition (Cf)).** BS 5148-5160 define the ring homomorphism G_{T_F}[A,F] -> G_{T_F}[A_cf,F_cf] on generators, by Cf(A) = A_cf, Cf(F) = F_cf. The lines fixing X and Y_i (BS 5159, 5174) are commented out, and BS 5163-5167 contradicts them (s.9 item 8). The reading adopted as the definition is: Cf is the identity on coefficients, Cf(X) = X and Cf(Y_i) = Y_i. It is forced by BS 5172 (D_cf = X/A_cf^2) and BS 5340 (the formula for Y_0^cf), and Proposition 13.2_P2 is a theorem about this Cf. On leg-free terms A has cancelled and Cf is F -> F_cf, well defined on functions by the transcendence of F (FB:C Step 4(i)). For a single coloured function Cf depends on the branch of t (t -> ωt relabels colours); only the colour-summed total is independent of it (FB:C Step 4(ii)). A general solution")
sub("everything in s.2.6 and s.3-s.6 is algebra on the closed form, including symplecticity (Proposition 2.10). The identification is needed for Proposition 2.11 and Proposition 13.2_P2.",
    "everything in s.2.6 except Proposition 2.11, and everything in s.3-s.6, is algebra on the closed form, including symplecticity (Proposition 2.10). The identification is needed for Proposition 2.11 and Proposition 13.2_P2, where it is carried by FB:B2, B4.")

# ---------------------------------------------------------------- s.1.4
sub("This lemma relates (0.1) to the code. The GW meaning of Dil^ (Givental's T(z) = z(1 - R^{-1}(z) 1)) is **not** derived here; it belongs to Proposition 13.2_P2.",
    "This lemma relates (0.1) to the code. The GW meaning of Dil^ (Givental's T(z) = z(1 - R^{-1}(z) 1)) is **not** derived here; it is FB:C Step 1 and Theorem FB_C (C2), for Proposition 13.2_P2.")

# ---------------------------------------------------------------- s.2.4
sub("**Proposition 2.6 (R2:Prop 5).** [P mod (C1), (C2); uses SP1-SP5] Fix generic α, and let R^α_{0i} be BS's canonical row (eq:canR; (C1): e^{u_i/z} R^α_{0i} formally solves eq:pf; the citation is eq:asymp with eq:canR, not (Flatness)).",
    "**Proposition 2.6 (R2:Prop 5).** [P mod eq:J=I and the genus-0 QDE, the reduced (C1) of s.0.3; uses SP1-SP5. (C2) is discharged: FB:A2.5, A3.3.] Fix generic α, and let R^α_{0i} be BS's canonical row (eq:canR; e^{u_i/z} R^α_{0i} formally solves eq:pf by FB:A2.4, A2.5, which replace the former citation of eq:asymp).")
sub("by (C1), Lemma 2.4 and Lemma 2.5.", "by FB:A2.4, A2.5, Lemma 2.4 and Lemma 2.5.")
sub("5. By (C2) (e^{-c_i} = prod_{k=0}^{4} Θ(-z/w_i^(k)) over the five tangent weights {3α_i, α_j-α_i, α_k-α_i, -ε1, -ε2}, Psi_{0i} = Delta_i^{-1/2} from Psi = Phi Delta^{-1/2}, BS 3014, eq:phi 3031), R^α_{0i}(0) = prod_5 w^{-1/2} Θ(-z/w).",
    "5. By eq:canR, with (Normalisation) a theorem (FB:A3.3: e^{-c_i} = prod_{k=0}^{4} Θ(-z/w_i^(k)) over the five tangent weights {3α_i, α_j-α_i, α_k-α_i, -ε1, -ε2}), and Psi_{0i} = Delta_i^{-1/2} (FB:A2.5): R^α_{0i}(0) = prod_5 w^{-1/2} Θ(-z/w).")
sub("So eq:normNLR (4327) is derived rather than cited.",
    "So eq:normNLR (4327) is derived rather than cited. FB:A3.5 writes this argument in the formal ring K^[[y]][[z]], including the transfer of Lemma 2.4 from an analytic domain to that ring, and finds c^2 = -1 with BS's printed Asym (c = 1 with (det Hess 𝒲)^{-1/2} in place of (det Hess(-𝒲))^{-1/2}). As a statement about the row eq:canR, (Normalisation) is its definition; that this row is Teleman's R, with these constants, is FB:A2.9, A3.3, modulo the literature of s.0.3.")
sub("**Theorem 2.7 (R2:Theorem 6, Prop 8, R3:Cor B; continuity and identification).** [P mod (C1)-(C3); uses SP1-SP5]",
    "**Theorem 2.7 (R2:Theorem 6, Prop 8, R3:Cor B; continuity and identification).** [P mod eq:J=I and the genus-0 QDE, the reduced (C1) of s.0.3; (C2), (C3) are discharged: FB:A2.5, A2.6, A3.3; uses SP1-SP5. For the role of this theorem after round 6c, see the status paragraph below.]")
sub("The dependence chain is SP1-SP5 (written) + Lemma 2.4 (written above) + (C1)-(C3).",
    "The dependence chain is SP1-SP5 (written) + Lemma 2.4 (written above) + eq:J=I and the genus-0 QDE, through Proposition 2.6.\n\nStatus after round 6c. FB:B4 writes step 2 in full, including the power count showing that the s_3-dependence is a pure power, which step 2 above does not show. FB:B4 also replaces step 5 by its identity π(Φ) = Θ(-z/E) S(v, z) in the single rational variable Y, which holds at all three roots with no chart or monodromy discussion. Parts (a)-(c) stay correct and give regularity locally near (y0, α*), but they are stated for germs continued along paths from y near 0. What the identification with Stage A's series needs is FB:B2: rationality in (α, Y) with denominator a power of W. So the route to Proposition 13.2_P2 is FB:B2, B4, B6, and Theorem 2.7 is kept as the analytic statement.")

# ---------------------------------------------------------------- s.2.6: Prop 2.11
line("**Proposition 2.11 (R2:Prop 11; Givental-Teleman on T_F).**",
"""**Proposition 2.11 (Givental-Teleman on T_F; replaced by Theorem FB_B of FRAMEB_PROOF.md).** [P mod the inputs of Theorem FB_B: Theorem FB_A; SP1-SP5, Lemmas 2.1, 2.3-2.5, Theorem A, Corollary 2.2, Propositions 2.10, 2.12; Graber-Pandharipande] Not used in Theorems LC_P2, R_P2. For g >= 2, in Q(ε1, ε2)((y)), F_g^{T_F}(y) equals the point-frame graph sum built from h_w times the rows of Definition 1.3 at c1 = A (Theorem FB_B(d)). *Proof.* FRAMEB_PROOF.md s.3. At generic α the graph sum is a multi-trace over a local ring 𝒞, and it is a polynomial in an indeterminate for F with coefficients in the localisation O of Q[α, ε, y] at the prime 𝔮 of T_F (Theorem FB_B(a), (c)). The identity is carried to T_F by the ring map α_i -> E/3 on a subring of power series in y, with a denominator outside 𝔮 (FB:B6). ∎

Continuation to the universal cover of C \\ {0, -1/27} means continuation of the right side; it is not a statement about GW series at large y. The earlier proof along eq:adelta is withdrawn. Its step from "δ divides Q" to "some φ_{k,m} blows up" is deleted, not repaired: at fixed y0 the graph sum is one number, so a pole of some φ_{k,m} need not show (round-5 critic item S1(c)). Non-divisibility is now the fact that the denominator is not in 𝔮, obtained a priori from Theorem FB_B(a). Algebraic independence, eq:adelta and rational tori are not used. BS line ~4056 is Theorem FB_B(b). BS's prop:fingenR (3611) is stated for rational tori only; at T_F, R_0's independence of A follows from Theorem 2.7(e), since the row is rational in t.""")

# ---------------------------------------------------------------- s.2.7 scope
sub("identifying row 0 of Definition 1.3 with BS's R_0 still needs Theorem 2.7(e), hence SP1-SP5, (C1) and (C2).",
    "identifying row 0 of Definition 1.3 with BS's R_0 still needs Theorem 2.7(e), now in the form FB:B2, B4, hence SP1-SP5, eq:J=I, the genus-0 QDE and the constants of FB:A3.")

# ---------------------------------------------------------------- s.7.5
sub("Identification with GW needs (C1)-(C5) and Proposition 13.2_P2, which remains open.",
    "Identification with GW is Proposition 13.2_P2 (s.0.2, s.8.3), which has its own inputs (s.8.2a).")

# ---------------------------------------------------------------- s.7.6
sub("For general g this reduction is a derivation sketch; it was tested exactly only at g = 2,",
    "For general g this reduction is now written in FB:C Step 9 [P mod PPZ Props 2.8, 2.9; not an input], so eq:betagiv is equivalent to Theorem VI. Before that it was tested exactly only at g = 2,")
line("**Corollary G.** [P mod Corollary SV_P2 and Proposition 13.2_P2 [O]]",
"""**Corollary G (BS rmk:cfgap).** [Status as in s.0.2: P mod Proposition 13.2_P2, Lemma 1.2(c), Theorem R_P2, Corollary SV_P2, Theorem VI; end-to-end inputs in s.8.2a; referee status in s.0.2, s.0.4.] For g >= 2, with GW^cf_g := Cf(F_g^{T_F}),

    GW^cf_g = (-3)^{g-1} β_{2g-2}(ε1, ε2)/t_cf^{2g-2} + O(1),

with O(1) as in s.0.2: holomorphic near t_cf = 0, Taylor coefficients in Q(ε1, ε2), value not determined, and either convention for d = 0.

*Proof (FB:C Step 7).* By Proposition 13.2_P2, GW^cf_g = Sum_Γ |Aut Γ|^-1 Sum_c T(Γ, c) at c1 = C1cf, up to the constant GW^{T_F}_{g,0} on the d > 0 convention. Each T(Γ, c) is a meromorphic germ at t_c = 0 (Corollary 2.8). Theorem R_P2 removes all polar terms of graphs with an edge, colouring by colouring. By Corollary SV_P2, the single vertex is analytic at colours 1, 2, and at colour 0 its polar part is (-3)^{g-1} β_{2g-2} t_c^{2-2g}. By Lemma 1.2(c), t_c = t_cf. For the Taylor coefficients: wc, t_c and F_cf lie in Q[[Δ_c]] and y^{-1} is regular at Δ_c = 0, so an element of Q(ε1, ε2)[y^{±1}, Δ_c^{-1}][F_cf] has a Laurent expansion in Δ_c over Q(ε1, ε2), and t_c is a coordinate over Q. The difference is such an element, up to the constant, and it has no negative powers. ∎

**Corollary G_{k,h}.** [Same status, plus the polynomiality of Theorem FB_B(b).] F^cf_{k,h} = 𝔟_{k,h} t_cf^{2-2h-2k} + O(1), as in s.0.2. *Proof (FB:C Step 7).* By Theorem FB_B(b), ε1 ε2 F_g^{T_F} = Sum_j E^{2g-2j} (ε1 ε2)^j f_j(y), with f_j in Q[[y]]. Solving at finitely many rational ε (Vandermonde) puts f_j in Q[y^{±1}, Δ_c^{-1}][F], and uniqueness of the polynomial in F gives Cf(GW_g) = Sum_j E^{2g-2j} (ε1 ε2)^{j-1} Cf(f_j). The Laurent coefficients of Cf(f_j) in t_c lie in Q and the monomials E^{2k} (ε1 ε2)^{h-1}, k + h = g, are linearly independent over Q. Comparing at each polar order in Corollary G gives the statement. The limits F^cf_{0,h}, F^cf_{h,0} are coefficient extractions, so E = 0 and ε1 = 0 need not be points of T_F. ∎

Comparison at g = 2, not a step: 𝔟_{2,0} = -7/1920, 𝔟_{1,1} = 7/480, 𝔟_{0,2} = -1/80 are the coefficients of X^2 in the three lines of eq:F2ref (checked by hand).""")

# ---------------------------------------------------------------- s.8.1
sub("""    Lemma 1.2(c) (t_c = t_cf)                        P mod C     BS 5145, eq:pfinf (read by round-4 critics);
                                                                 used only for Corollary G
""", """    Lemma 1.2(c) (t_c = t_cf, A_cf = C1cf, branch)   P mod C     DLMF 15.6.7, Barnes's second lemma, Gauss's ψ values; FB:C Step 6;
                                                                 used only for Corollary G
""")
sub("""    Prop 2.6 (normalisation; eq:normNLR derived)     P mod C     (C1), (C2); R2:Prop 5; interface only
    Thm 2.7 (continuity, closed form, Rec = Wick)    P mod C     (C1)-(C3); R2:Thm 6, Prop 8; R3:Cor B; interface only;
                                                                 reading of BS 5257-5267 re-checked by referee
""", """    Prop 2.6 (normalisation; eq:normNLR derived)     P mod C     eq:J=I, genus-0 QDE; SP1-SP5; (C2) discharged (FB:A2.5, A3.3);
                                                                 R2:Prop 5; formal form FB:A3.5; interface only
    Thm 2.7 (continuity, closed form, Rec = Wick)    P mod C     eq:J=I, genus-0 QDE; SP1-SP5; (C2), (C3) discharged (FB:A2.6);
                                                                 R2:Thm 6, Prop 8; R3:Cor B; steps 2, 5 of (e) as in FB:B4; interface only;
                                                                 reading of BS 5257-5267 re-checked by referee
""")
sub("    Prop 2.11 (Givental-Teleman on T_F)              P mod C     (C1)-(C5); not used for Theorems LC_P2, R_P2; not refereed in round 6\n",
    "    Prop 2.11 (Givental-Teleman on T_F)              P mod C     = Theorem FB_B(d); step \"δ divides Q\" deleted; not used for Theorems LC_P2, R_P2\n")
sub("""    Proposition 13.2_P2 (G_g = GW^cf_g)              O           s.8.3
    Corollary G (BS rmk:cfgap)                       Conditional on Proposition 13.2_P2. NOT CLAIMED.
""", """    Theorem FB_A (CohFT, R, constants, generic α)    P mod C     eq:J=I, genus-0 QDE; Teleman Thm (2), Prop 7.7, Prop 8.7(ii); PPZ; GP; Mumford, FP Thm 4;
                                                                 A3.5 also SP1-SP5, Lemmas 2.1, 2.3-2.5; two referee passes on the report form
    Theorem FB_B (specialisation to T_F)             P mod C     Theorem FB_A; SP1-SP5, Lemmas 2.1, 2.3-2.5, Thm A, Cor 2.2, Props 2.10, 2.12; GP;
                                                                 one referee pass on the report form, conditional on Stage A
    Theorem FB_C (frame change, dictionary)          P mod C     Stages A, B; Prop 2.10; Mumford through A3.2; PPZ Prop 2.14;
                                                                 two referee passes on the report form
    Definition (Cf)                                  definition  BS 5148-5160 read through generators, forced by BS 5172, 5340; s.1.3
    Proposition 13.2_P2 (G_g = GW^cf_g, d >= 0)      P mod C     Stages A, B; Props 2.10, 2.12; FB:C Step 5; s.8.2a
    Corollary G (BS rmk:cfgap)                       P mod (MV)  and the inputs of FRAMEB_PROOF.md s.1 (s.8.2a); Prop 13.2_P2, Lemma 1.2(c),
                                                                 Theorem R_P2, Cor SV_P2, Theorem VI; assembled text not refereed (s.0.4)
    Corollary G_{k,h} (BS def:dirint, conifold)      P mod (MV)  same, plus polynomiality of Theorem FB_B(b)
    Remark (eq:betagiv equivalent to (VI))           P mod C     PPZ Props 2.8, 2.9; FB:C Step 9; not an input
    F_{k,h} in Q[X^{±1}, F] for all g                not proved  s.0.4; independent of Corollary G
""")

# ---------------------------------------------------------------- s.8.2
sub("Identification with GW still needs (C1)-(C5) and Proposition 13.2_P2, s.8.3.",
    "Identification with GW is Proposition 13.2_P2, whose inputs are in s.8.2a.")
sub("F2 has **not** been taken from thm:cfasym: since round 2 it follows from the closed form, and eq:Rcf is then a consistency check.\n",
    """F2 has **not** been taken from thm:cfasym: since round 2 it follows from the closed form, and eq:Rcf is then a consistency check.

### 8.2a Inputs of Proposition 13.2_P2 and Corollaries G, G_{k,h}

Dependency: Stage A -> Stage B -> Stage C (Proposition 13.2_P2). Lemma (t_cf) is independent of Stages A, B. Corollary G = Proposition 13.2_P2 + Lemma 1.2(c) + Theorem R_P2 + Corollary SV_P2 + Theorem VI. Details: FRAMEB_PROOF.md s.0.3, s.1.

Cited, end to end:
1. **BS (allowed).** eq:J=I 2663 and the genus-0 QDE, with eq:SHb 3305. eq:J=I is cited by BS to MR2486673, not opened. Definitions only: eq:Ifun 2621, eq:pf 2628, eq:Ian 2647, eq:ADEFXdef 2681, eq:etaflat 2419, eq:spotP2 4301, Asym at about 4245, t_cf at 5145.
2. **Teleman**, Theorem (2), Prop 7.7, Prop 8.7(ii). **PPZ**, Defs 2.2, 2.5, 2.13, Prop 2.14. **GP**, formula (1), Apps. A, C. **Mumford** as in **FP** s.1.1, and FP Theorem 4.
3. **DLMF 15.6.7, Barnes's second lemma, Gauss's values of ψ at 1/3, 2/3.** Used only for Lemma 1.2(c).
4. **(MV)** with the textbook facts of s.8.2 items 1, 2. Used only through Theorems R_P2, VI and Corollary SV_P2.

Proved in this text and used: SP1-SP5, Lemmas 2.1, 2.3-2.5, Theorem A, Corollaries 2.2, 2.8, Propositions 2.10, 2.12, Theorem LC_P2. BFGW Lem 4.2 is not among the inputs: it enters only for statements involving A ((C4), s.0.3).

Not inputs: BS rmk:refmm, BKMP/FLZ remodelling, HK10, Coates-Iritani 2018, BFGW Conj 1.9, BS prop:cfgraphsum, thm:cfasym, lem:cfpole, TR on the local-P2 mirror curve, and any finite-order computation.

What the referee passes did not cover is FRAMEB_PROOF.md s.1.4: the proofs of Teleman's Prop 7.7 and Theorem (2) (read only); SP1-SP4 line by line and the proof of Lemma 2.4; eq:J=I, GP's proof, (MV); Prop 2.10 steps S1-S6 and Theorem A in that round; the class-level identity Omega = R.omega beyond integrated numbers (it rests on Teleman); the assembled text.
""")

# ---------------------------------------------------------------- s.8.3
i0 = s.index("### 8.3 Frame-B formula on T_F (Proposition 13.2_P2) and the gap corollary as the final step")
i1 = s.index("## 9. Errata and imprecisions found in BS v2")
s = s[:i0] + """### 8.3 Frame-B formula on T_F (Proposition 13.2_P2) and the gap corollary

**Proposition 13.2_P2** is stated in s.0.2, with status [P mod Stages A, B of FRAMEB_PROOF.md and Propositions 2.10, 2.12]. Its proof of record is FRAMEB_PROOF.md: Stage A s.2, Stage B s.3, Stage C s.4 (C Steps 0-5). Lemma (t_cf) is C Step 6, here Lemma 1.2(c). The gap corollaries are s.0.2, s.7.6 (C Step 7). No TR for the target theory enters. Inputs: s.8.2a. Referee status: s.0.4.

The list of what remained in the earlier versions of this section is closed, as follows.
1. *Frame B.* Omega_w = h_w·1 with markings is FB:A3.2. h_w.ω = Omega_w (x) ω and R' = R_{T_F} h_w^{-1} are FB:C Step 2. The dilaton leaf is q_t = T' - ψ, with T' the Givental translation of R' and -ψ the dilaton shift: FB:C Step 1, Theorem FB_C (C2). The worry about conifold Lemma 13.3 is void: h_{wt_i} is singular on the diagonal, but Stage B's Φ absorbs it, and the constants are fixed at generic α and carried to T_F by a ring map (FB:A3, B2, B4). The route through a family and pointwise continuation is not used.
2. *N4.* lem:modrel: Proposition 2.12. N1, N2, N5: Proposition 2.6, Theorem 2.7, Definition 1.3, now with FB:A2, B4.
3. *BS line ~4056.* Theorem FB_B(b), FB:B5 [P mod GP].
4. *N3.* FB:C Step 4 and Definition (Cf). Independence of A and F is needed only for statements involving A; on leg-free terms the transcendence of F suffices.
5. *Verification.* The reports had the referee passes listed in s.0.4. The assembled text has not; see the marker ROUND6C-REFEREES in s.0.4.

**Remark on N4 and Proposition 2.10.** In round 4, Proposition 2.10 step 3 compared with BS's rows at c1 = A, which made lem:modrel a load-bearing input of Theorem LC_P2, and the draft wrongly attributed c2 = D to eq:ADEFXdef (it is the first relation of lem:modrel). Both issues are gone: Proposition 2.10 is now proved from Theorem A, and lem:modrel is Proposition 2.12. The audit of BS's own use of lem:modrel is Remark 2.13.

Round-1 items N1-N8 are recorded in journal wf_e928f771-9e3 (only N4 appears in ROUND1_STATUS_2026-11-02.md).


""" + s[i1:]

# ---------------------------------------------------------------- s.9
line("5. **BS ~4344-4346**",
"5. **BS 4345-4347** assert regularity at α_i = E/3 and eq:oscintP \"by continuity\". Fix: Theorem 2.7(a) near (y0, α*). At every z-order, in the form used by the proof of record: FB:B2 and the analytic corollary of FB:B4, where all factors are rational in (α, ε, y, Y) with no pole. No limit is needed for Proposition 13.2_P2 (item 19).")
line("6. **BS 3497-3510 (Normalisation)**",
"6. **BS 3497-3510 (Normalisation)** is asserted with five-weight constants without proof. Fix: FB:A3.3, a theorem at generic α with BS's sign, modulo the literature of s.0.3. Proposition 2.6 is its consequence for the oscillating integral.")
line("8. **sec:cfgen (5140).**",
"8. **sec:cfgen (5140), lines 5163-5167.** Cf is defined on generators, Cf(A) = A_cf, Cf(F) = F_cf (5148-5160), and the lines fixing X and Y_i (5159, 5174) are commented out. BS add that in terms of quasi-modular forms the map is induced by Cf(f)(τ) = f(-1/(3τ)). This is false as a description of Cf. By eq:tau (2794) and the transformation of η, τ -> -1/(3τ) sends u := -27y to 1 - u, so it sends X = 1/(1-u) to 1/u, which contradicts D_cf = X/A_cf^2 at 5172. GW_g(-1/(3τ)) is regular at the conifold, since it is GW_g near y = 0, which contradicts thm:cfasym. F fails too, since θ does not commute with this map. Corrected statement: Cf is the substitution on generators, extended to A^{-1}, with identity on coefficients, Cf(X) = X and Cf(Y_i) = Y_i (s.1.3; forced by 5172 and 5340). What survives of the display concerns the generator A only: A(-1/(3τ')) = -A_cf, with τ' = τ - 1 on the branch of item 16.")
line("12. **Line ~4056:**",
"12. **Lines 4055-4057.** Regularity of GW^{T_max} on the diagonal is argued in one clause, by properness of the T_F-fixed locus and lem: GW contravariance, followed by \"by continuity\" from rational weights. The T_max-fixed loci change on the diagonal, and individual T_max-localisation terms are singular there, so the clause needs a localisation statement over the prime of T_F, which is not written. Fix: Theorem FB_B(b), FB:B5 [P mod GP], which rests on H^0(C, f^*O(-3)) = 0 on the proper Mbar_g(P2, d) for d >= 1, and on e_3 for d = 0. It gives more: ε1 ε2 GW^T_{g,d} lies in Q[α, ε] for d >= 1. The extension of thm:haegen(i) by continuity is replaced by Theorem FB_B(e).")
sub("Earlier round notes quote slightly different numbers.\n",
    """Earlier round notes quote slightly different numbers.
15. **Line 3537.** f_j^* is printed as the dual basis "at t". It must be the constant covector at Q = 0: with a moving covector the (Flatness) claimed at 3543 fails, and 3559 uses the value at Q = 0. Fix: FB:A2.2.
16. **Line 5145.** On y in (-1/27, 0), t_cf needs the branch log y = log(-y) + iπ. On the other branch the printed expression equals t_cf - sqrt3 π - i sqrt3 (l + g), with l = log(-y), g = I_0^{[1]} - log y, and it does not vanish at y = -1/27 (Lemma 1.2(c)). On the correct branch eq:tau places the conifold at τ = 1, with A_cf = i sqrt3 (τ - 1) A, which is T-equivalent to the τ = 0 stated at 2798.
17. **Lines 1277, 2446 versus eq:gtH (3799) and eq:F2ref (4217).** GW_g and the correlators of eq:gwprim are printed with β != 0, resp. d > 0, while the graph sum includes d = 0. The y = 0 values (X = 1, F = 0) of eq:F2ref are 0, 1/640, 0 in F_{0,2}, F_{1,1}, F_{2,0}, that is GW^{T_F}_{2,0} = E^2/640, which is the constant of Theorem FB_B(b). Fix: read GW_g with d >= 0, or read the identities up to the constant (s.0.2).
18. **Lines 4196-4198.** The ansatz for ρ_{i,n}^{[0]} is said to give a full-rank system by direct inspection "for the first few values of n". Not needed: by FB:B2, Theorem FB_B(a), the z-coefficients of the normalised row lie in the ring 𝒞 at every order, with denominators only powers of W, ε1, ε2 and no factor of the discriminant.
19. **Line 3983.** T_F as a "limiting version" of a rational torus, and the limit δ -> 0 in eq:adelta. Not needed: the specialisation is a ring map α_i -> E/3 on power series in y with coefficients in the local ring at the prime of T_F, with a denominator outside that prime (FB:B6). The normalised rows contain no square roots, so the reason for passing through rational tori disappears. Coefficientwise in y a coloured factor blows up on the diagonal, so a limit of the y-series would not be valid in any case.
20. **prop:algindABC (2847-2854), attribution.** BS cite BFGW Lem 4.2 for the independence of A, B, C. According to referee refFB_C, who opened arXiv:2001.05347, that independence is in BFGW s.4.1 (Lem 4.1 plus Zagier Prop 20), and Lem 4.2 is its corollary, reading "S, X, I11 are algebraically independent", with I11 = A and S = F - (X-1)/3. eq:ABCmod is called well known at 2832 and is cited to other papers in BFGW. Used here only for statements involving A ((C4), s.0.3).
21. **Further items, not repeated here:** FRAMEB_PROOF.md s.5.1 item 3 (the constant c in eq:normNLR with the printed Asym) and item 9 (lines 4071-4076: regularity at y = 0 is argued for the coefficients of the polynomial in F, while what is given is regularity of the total, and F = O(y); no counterexample is known, s.0.4).
""")

# ---------------------------------------------------------------- Appendix
sub("Copies of the round 1-3 scripts are in LP/round1_scripts, LP/round2_scripts, LP/round3_scripts. Scratchpad paths may be purged.\n",
    """Copies of the round 1-3 scripts are in LP/round1_scripts, LP/round2_scripts, LP/round3_scripts. Scratchpad paths may be purged.

**Round 6b (frame-B; LP/round6_scripts/<key>/, with key FB_A, FB_B, FB_C for the authors and refFB_A, refFB_A2, refFB_B, refFB_C, refFB_C2 for the referees).** Listed with outcomes in FRAMEB_PROOF.md Appendix s.A.1-s.A.3; limits and warnings in s.A.4, including that genus 2 sees only F_cf(0) at the O(1) level. Coverage is finite: genus <= 3 for most items, genus 4 with one or two weight pairs, and nothing to all orders. Evidence only; no proof step of FRAMEB_PROOF.md rests on it.
""")

if errors:
    for e in errors:
        print("ERR", e)
    sys.exit(1)
open(P, "w", encoding="utf-8", newline="").write(s)
print("ok", len(s))
