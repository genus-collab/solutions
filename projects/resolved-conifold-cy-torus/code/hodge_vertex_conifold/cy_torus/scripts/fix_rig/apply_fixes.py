import sys, io
path = sys.argv[1]
raw = io.open(path, encoding="utf-8", newline="").read()
crlf = raw.count("\r\n")
lf_only = raw.count("\n") - crlf
print("CRLF lines:", crlf, " bare-LF lines:", lf_only)
if crlf and lf_only:
    raise SystemExit("mixed line endings; refusing")
src = raw.replace("\r\n", "\n")
orig = src

def rep(old, new):
    global src
    n = src.count(old)
    if n != 1:
        raise SystemExit("MATCH COUNT %d for:\n%s" % (n, old[:200]))
    src = src.replace(old, new)

# --- header status line
rep("Nothing in this file has been refereed.",
    "Referee status: one hostile referee pass (AI), verdict \"sound with minor fixes, no major items\"; its items are applied (record in s.8.4). Nothing in this file has been read by a human.")

# --- s.0.3 Theorem Rig: genus-1 input is [BS] for d >= 2 only
rep("and genus 1 (Lemma 5.7, P mod [BS] prop:ResConifold at genus 1), this is [BS] Conjecture",
    "and genus 1 (Lemma 5.7, P mod [BS] prop:ResConifold at genus 1 for d >= 2; d = 1 by direct localisation), this is [BS] Conjecture")

# --- item 6: s.3 first sentence
rep("Those parts are written at indeterminate weights (FEAS (5); the critic's reading agrees), and the places where the weights enter are the following, each of which holds on T_CY:",
    "Those parts treat the five weights as indeterminates and use no relation among them (FEAS (5); the critic's reading agrees); every hypothesis they need about the weights is a non-vanishing or non-constancy statement, listed below, and each holds in Kw = Q(a0, a1, a2, e1):")

# --- item 6: alternative route
rep(" Part (d) of LSP (the passage to the non-semisimple T_F) is not needed: Q = 0 is semisimple on T_CY. What replaces it",
    " Part (d) of LSP (the passage to the non-semisimple T_F) is not needed: Q = 0 is semisimple on T_CY. (Alternative route, and a cross-check of the reading above: apply LSP Prop. 13.2 (a)-(c), (e) at the five indeterminate weights with the frame w = (a0 + a1 + a2, -e1, -e2), and specialise the resulting identity along e2 -> E - e1 exactly as in LSP (d): the colour-summed total is num/den with den|_{T_CY} != 0, because every per-term denominator is one of the factors listed in the description of the coloured terms below (W_2(U_i), 1 - U_i, 1 - Q U_i, y, B Q, and the constants a0, a1, A, e1, e2, E) and none of them involves E - e1 - e2; on the GW side GW^{T_CY}_{g,d} = GW^T_{g,d}|_{T_CY} by Lemma 5.2 for d >= 1 and, for d = 0, by the constant-map formula of the Remark after Lemma 5.2, whose denominators a0 a1 a2 e1 e2 and a0 A B e1 e2 do not vanish on T_CY. Either route gives (3.1).) What replaces it")

# --- item 7a: s.3 last paragraph
rep("by (b)(iii) every coloured term is in Kw[[Q]], so the total is regular at Q = 0 too,",
    "by (b)(iii) the expansion of every coloured term at either of the two (unramified, s.1.4(v)) points of C over Q = 0 lies in Kw[[Q]] (the rows at both colours lie in Kw[[Q]][[z]], and the two points are exchanged by the deck involution together with the colours), so the total is regular at Q = 0 too,")

# --- s.3 status line: referee pass
rep("It has had no referee pass in this form; the parts (a)-(c), (e) of LSP had the passes recorded in LSP s.13.2.]",
    "This reading has had one referee pass (AI, s.8.4), which found minor items only (applied: the first sentence of the proof, the alternative route, the wording of the Q = 0 regularity); the parts (a)-(c), (e) of LSP had the passes recorded in LSP s.13.2.]")

# --- item 7b: Lemma 5.1
rep("and every coloured term lies in Kw[[Q]] (LSP (b)(iii), Q = 0 semisimple because a0 != 0).",
    "and the expansion of every coloured term at either of the two (unramified) points of C over Q = 0 lies in Kw[[Q]] (LSP (b)(iii) at both colours, Q = 0 semisimple because a0 != 0; Q = 0 is unramified by s.1.4(v)).")

# --- item 5: Lemma 5.2
rep("is (1/(e1 e2)) times an honest T-equivariant class on M, with coefficients in Q[a1, a2, e1, e2]; the localisation definition",
    "is (1/(e1 e2)) times an honest T-equivariant class on M (e_T(R^1_a) depends on a0 through the base rotation); the localisation definition")

# --- item 7c: Lemma 5.5
rep("**Lemma 5.5.** F_g is regular at Q = oo (as a1 != a2 in Kw).\nProof. In s := 1/Q, by s.1.5",
    "**Lemma 5.5.** F_g is regular at Q = oo (as a1 != a2 in Kw).\nProof. Q = oo is unramified for C (s.1.4(v): the leading coefficient (a1 - a2)^2 of Disc is non-zero), so the two points of C over s := 1/Q = 0 are the two colours, and regularity has to be shown at both of them; this is what follows. In s, by s.1.5")

# --- item 3: Lemma 5.7
rep("so e1 e2 GW_{1,d} = e1 e2 N_d + E L_d with N_d :=",
    "so e1 e2 GW_{1,d} = e1 e2 N_d + (e1 + e2) L_d (= e1 e2 N_d + E L_d on T_CY), an identity on the 5-torus, with N_d :=")
rep("e1 e2 GW_{1,d} = [t^2] S(e1 t) S(e2 t)/d = -(e1^2 + e2^2)/(24 d) = e1 e2/(12 d) - E^2/(24 d); comparing with e1 e2 N_d + c E^2 (the quadratic forms e1 e2 and E^2 = (e1 + e2)^2 are linearly independent over Q) gives N_d = 1/(12 d) and c = -1/(24 d). Hence e1 e2 GW_{1,d} = e1 e2/(12 d) - E^2/(24 d) = P_1/d on T_CY, using E = e1 + e2. QED. (A direct computation of the single linear form L_d, or a genus-1 analogue of Proposition 13.2_CY, would remove the [BS] input; neither is in the material. At d = 1 the value L_1 = -E/24 agrees with the genus-1 part of [BS]'s display at line 1694, which is evidence.)",
    "e1 e2 GW_{1,d} = [t^2] S(e1 t) S(e2 t)/d = -(e1^2 + e2^2)/(24 d) = e1 e2/(12 d) - (e1 + e2)^2/(24 d); there (e1 + e2)(a0 + a1 + a2) = (e1 + e2)^2, so comparing with e1 e2 N_d + c (e1 + e2)^2 (the quadratic forms e1 e2 and (e1 + e2)^2 are linearly independent over Q) gives N_d = 1/(12 d) and c = -1/(24 d). Hence e1 e2 GW_{1,d} = e1 e2/(12 d) - (e1 + e2)(a0 + a1 + a2)/(24 d) on the 5-torus, which on T_CY is e1 e2/(12 d) - E^2/(24 d) = P_1/d, using E = e1 + e2 = a0 + a1 + a2. QED. (A direct computation of the single linear form L_d for all d, or a genus-1 analogue of Proposition 13.2_CY, would remove the [BS] input; for d = 1 it is done next, for d >= 2 it is not in the material.)\n\n *The case d = 1 by direct localisation* [P mod the localisation formula (L) below ([GP99] on the proper Mbar_g(P^1, 1), in the standard form; its ingredients are listed and its normalisation is confirmed by g = 0), [M83], and Int_{Mbar_{1,1}} lambda_1 = Int_{Mbar_{1,1}} psi_1 = 1/24]. For d = 1 the T-fixed loci of Mbar_g(P^1, 1) are indexed by g_+ + g_- = g: the zero section covered once, with a contracted component of genus g_p attached at p (absent if g_p = 0), fixed locus Mbar_{g_+,1} x Mbar_{g_-,1}. Collecting the factors (vertex at p: H^1 and H^0 of T_Z on the contracted component, prod_a Lambda^vee(w^{(p)}_a)/prod_a w^{(p)}_a; node: prod_a w^{(p)}_a from the normalisation sequence and 1/(w^{(p)}_1 - psi) from the smoothing, w^{(p)}_1 = +-a0 the tangent weight of P^1 at p; edge: H^0 = H^1 = 0 for O(-1), H^0(P^1, O_{-e_i}) giving 1/(e1 e2), and of H^0(P^1, T P^1) with weights (a0, 0, -a0) the part not cancelled by the automorphisms of the edge fixing its nodes, i.e. 1/w^{(p)}_1 for each node),\n\n    e1 e2 GW_{g,1} = Sum_{g_+ + g_- = g} V_+(g_+) V_-(g_-),   V_p(0) := 1,   V_p(h) := Int_{Mbar_{h,1}} prod_{a=1}^5 Lambda^vee_h(w^{(p)}_a)/(w^{(p)}_1 (w^{(p)}_1 - psi)),            (L)\n\n with Lambda^vee_h(w) = Sum_i (-1)^i lambda_i w^{h-i}; at g = 0 it gives GW_{0,1} = 1/(e1 e2), the value of Prop. 13.2(a). At h = 1, Lambda^vee_1(w) = w - lambda_1 and lambda_1^2 = 0, so with Pi_p := prod_{a=1}^5 w^{(p)}_a and 1/(w_1 (w_1 - psi)) = w_1^{-2} (1 + psi/w_1 + ...) the degree-1 part of the integrand is (Pi_p/w_1^2) [psi/w_1 - lambda_1 Sum_a 1/w_a], and Int psi = Int lambda_1 = 1/24 gives\n\n    V_p(1) = (Pi_p/(24 w_1^2)) [1/w_1 - Sum_{a=1}^5 1/w_a] = -(Pi_p/(24 w_1^2)) Sum_{a=2}^5 1/w_a.\n\n At p_+ (w = (a0, a1, a2, -e1, -e2), Pi_+ = a0 a1 a2 e1 e2): V_+(1) = -(1/(24 a0)) [(a1 + a2) e1 e2 - a1 a2 (e1 + e2)]. At p_- (w = (-a0, B, A, -e1, -e2), Pi_- = -a0 A B e1 e2, w_1^2 = a0^2): V_-(1) = (1/(24 a0)) [(A + B) e1 e2 - A B (e1 + e2)]. Since A + B - a1 - a2 = 2 a0 and A B - a1 a2 = a0 (a0 + a1 + a2),\n\n    e1 e2 GW_{1,1} = V_+(1) + V_-(1) = e1 e2/12 - (a0 + a1 + a2)(e1 + e2)/24    on the 5-torus,\n\n i.e. N_1 = 1/12 and L_1 = -(a0 + a1 + a2)/24 without [BS]; on T_CY this is e1 e2 GW_{1,1} = P_1. So Lemma 5.7 holds at d = 1 modulo (L) and [M83] only, and the [BS] input is used for d >= 2 only. [E: the referee's independent direct-localisation test at d = 1 gave the same e1 e2 GW_{1,1} on T_max, and e1 e2 GW_{2,1} = P_2 on T_CY; both agree with the genus <= 2 part of [BS]'s display at line 1694 off T_CY.]")

# --- s.5.6 status
rep("g = 1 additionally mod [BS] prop:ResConifold; g = 0 mod [FP00].",
    "g = 1 additionally mod [BS] prop:ResConifold for d >= 2 (d = 1 by direct localisation, Lemma 5.7); g = 0 mod [FP00].")

# --- item 8: Remark 6.1 status
rep("[E: FEAS crit.tested, flop_const.log, g = 2..6 at two CY samples; fails off CY.]",
    "[E: FEAS crit.tested, flop_const.log, g = 2..6 at two CY samples; fails off CY. The critical points at s = 0, the W_2 limits (1.5) and the flop weights wt'_+- were also checked exactly by the referee (evidence).]")

# --- s.6.3 first bullet
rep("uses it at g = 1 (Lemma 5.7), and does not touch g = 0.",
    "uses it at g = 1 for d >= 2 (Lemma 5.7; d = 1 is done by direct localisation there), and does not touch g = 0.")

# --- item 4: s.6.3 second bullet
rep("GW_{g,d} = -P_g(e, -e) d^{2g-3}/e^2 = -(1/d) [e^{2g-2}] (2 sinh(d e/2))^{-2} (e := e1 = -e2; P_g(e, -e) = [t^{2g}] S(e t)^2),",
    "GW_{g,d} = -P_g(e, -e) d^{2g-3}/e^2 = the e^{2g-2} term of -(1/d) (2 sinh(d e/2))^{-2} (e := e1 = -e2; P_g(e, -e) = [t^{2g}] S(e t)^2 = e^{2g} [t^{2g}] S(t)^2, and (2 sinh(d e/2))^{-2} = (d e)^{-2} S(d e)^2),")

# --- item 1: s.7 proof (1)
rep("The eigenvalues L_i and the metric eta^{-1} = e1 e2 [[0, a1 a2], [a1 a2, -E]] = (e1 e2/x) a2^2 eta_F^{-1}(x) are homogeneous of degree 1 and 2 in (a1, a2), and D_i = d/dt - Delta_i'/(2 Delta_i) is scale-invariant. So the frame-B rows on T_3 are the fibrewise frame-B rows at x = a1/a2 with z -> z/a2 (and r_{1,i} multiplied by a2), Delta_i = (e1 e2/x) a2^3 Delta^F_i(x), and every Dil_i(d), E_ij(a, b) is the fibrewise one times an explicit monomial in a2 and e1 e2/x.",
    "The eigenvalues L_i = a2 L^F_i(x) have degree 1 in (a1, a2) (tau L^2 - E L + a1 a2 = a2^2 [tau (L/a2)^2 - c (L/a2) + x], c := 1 + x, the fibrewise quadratic of LSP (C5)). The metric inverse eta^{-1} = e1 e2 [[0, a1 a2], [a1 a2, -E]] is not a multiple of the fibrewise eta_F^{-1} = [[0, x^2], [x^2, -c x]] (LSP (C4)): its entries have different degrees, e1 e2 a1 a2 = (e1 e2/x) a2^2 x^2 (degree 2) and -e1 e2 E = (e1 e2/x) a2 (-c x) (degree 1). D_i = d/dt - Delta_i'/(2 Delta_i) is scale-invariant. Since r_{0,i}(z) = r^F_{0,i}(z/a2) and hence, by (C9), r_{1,i}(z) = (L_i + z D_i) r_{0,i}(z) = a2 [L^F_i + (z/a2) D_i] r^F_{0,i}(z/a2) = a2 r^F_{1,i}(z/a2), the pairing <v_i(z), v_j(w)> = eta^{01} [r_{0,i}(z) r_{1,j}(w) + r_{1,i}(z) r_{0,j}(w)] + eta^{11} r_{1,i}(z) r_{1,j}(w) equals (e1 e2/x) a2^3 times the fibrewise pairing at (z/a2, w/a2) (the eta^{01} term carries a2^2 x^2 times one factor a2, the eta^{11} term a2 (-c x) times a2^2), consistent with Delta_i = e1 e2 (2 a1 a2 L_i - E L_i^2) = (e1 e2/x) a2^3 Delta^F_i(x), Delta^F_i = x (2 x L^F_i - c (L^F_i)^2) (LSP (4.1)). Hence, by LSP (C10) (Dil_i(d) = (-1)^d r_{0,i}[d-1]), Dil_i(d) = a2^{1-d} Dil^F_i(d), and by LSP (C11) (E_ij(m, n) is a signed sum of <v_i[p], v_j[q]> over p + q = m + n + 1), E_ij(m, n) = (e1 e2/x) a2^{2-m-n} E^F_ij(m, n). So the frame-B rows on T_3 are the fibrewise frame-B rows at x = a1/a2 with z -> z/a2 (and r_{1,i} multiplied by a2), and every Delta_i, Dil_i(d), E_ij(m, n) is the fibrewise one times the explicit monomial in a2 and e1 e2/x just stated.")

# --- item 2: s.7 proof (4)
rep("and GW_{g,0} is regular at a0 = 0 (the 1/a0 terms of the two fixed-point contributions cancel, as on T_F);",
    "and GW_{g,0} = Sum_p Int_{M_g} prod_{a=1}^5 Lambda^vee(w^{(p)}_a)/w^{(p)}_a is regular at a0 = 0: only the factors Lambda^vee(a0)/a0 (at p_+) and Lambda^vee(-a0)/(-a0) (at p_-) have a pole there, a simple one, and the two residues cancel. Indeed the a0^{-1} coefficient of the p_+ term is (-1)^g Int_{M_g} lambda_g Lambda^vee(a1) Lambda^vee(a2) Lambda^vee(-e1) Lambda^vee(-e2)/(a1 a2 e1 e2) (from Lambda^vee(a0)/a0 = Sum_i (-1)^i lambda_i a0^{g-i-1}, whose a0^{-1} term is (-1)^g lambda_g/a0; the other four factors at p_+ do not involve a0), and that of the p_- term is minus the same (from Lambda^vee(-a0)/(-a0) = Sum_i (-1)^i lambda_i (-a0)^{g-i-1}, whose a0^{-1} term is (-1)^g lambda_g/(-a0) = -(-1)^g lambda_g/a0, times the value at a0 = 0 of the factor Lambda^vee(B) Lambda^vee(A) Lambda^vee(-e1) Lambda^vee(-e2)/(A B e1 e2), which is regular there with A -> a2, B -> a1). Moreover GW^{T_3}_{g,0} in the [BS] localisation definition equals this value at a0 = 0, by the same contravariance/localisation argument as for d >= 1 (Lemma 5.2): the T_3-fixed locus of Mbar_g(Z, 0) = Mbar_g x Z is the proper Mbar_g x P^1 (the zero section is fixed pointwise by T_3, its normal weights a1, a2, -e1, -e2 are non-zero), with virtual class e(E^vee (x) T_{P^1}) ∩ [Mbar_g x P^1] (the T_3-fixed part of the obstruction bundle E^vee (x) T_Z) and virtual normal bundle Nb (x) (O - E^vee); the same integral taken T-equivariantly over the proper Mbar_g x P^1 lies in Q[a0, a1, a2, e1, e2][1/(a1 a2 A B e1 e2)] (its denominators are the T-weights of Nb at p_+-), restricts to T_3 by contravariance, and equals Sum_p prod_{a=1}^5 Lambda^vee(w^{(p)}_a)/w^{(p)}_a by localisation at the two T-fixed points of P^1 ([GP99]; at p the factor Lambda^vee(+-a0)/(+-a0) is e(E^vee (x) T_pP^1)/e(T_pP^1));")

# --- s.8.1 table rows
rep("    Proposition 13.2_CY                     P mod the inputs of LSP Prop. 13.2 (list in s.3), read at T_CY; no referee pass of this reading",
    "    Proposition 13.2_CY                     P mod the inputs of LSP Prop. 13.2 (list in s.3), read at T_CY; one AI referee pass of this reading (s.8.4), minor items applied")
rep("    Lemma 5.7 (g = 1)                       P mod [BS] prop:ResConifold at g = 1, [GP99], [M83]",
    "    Lemma 5.7 (g = 1)                       d = 1: P mod the degree-1 localisation formula (L), [M83]; d >= 2: P mod [BS] prop:ResConifold at g = 1, [GP99], [M83]")
rep("    Remark 6.1 (flop)                       C; limiting rows derived by hand, Lemma 13.3 argument at s = 0 not re-written",
    "    Remark 6.1 (flop)                       C; limiting rows derived by hand, Lemma 13.3 argument at s = 0 not re-written; critical points, W_2 limits, flop weights checked exactly by the referee (E)")

# --- s.8.2 item 1
rep("A direct evaluation of the linear form L_d = -(a0 + a1 + a2)/(24 d) would remove it.",
    "A direct evaluation of the linear form L_d = -(a0 + a1 + a2)/(24 d) would remove it; this is done at d = 1 (Lemma 5.7, direct localisation), so the [BS] input is used for d >= 2 only.")

# --- s.8.3 open points
rep(" - A referee pass on this file, in particular on Lemma 5.4 (the a0-adic step), on the claim of Lemma 5.1 that the only possible poles are over {0, 1, oo} and Disc = 0, and on s.3.\n - A conceptual proof of the flop identity of s.6.1, and a direct proof of Lemma 5.7 without [BS].",
    " - A human reading of this file (one AI referee pass is recorded in s.8.4), in particular of Lemma 5.4 (the a0-adic step), of the claim of Lemma 5.1 that the only possible poles are over {0, 1, oo} and Disc = 0, and of s.3.\n - A conceptual proof of the flop identity of s.6.1, and a direct proof of Lemma 5.7 without [BS] for d >= 2 (d = 1 is done).")

# --- s.8.4 referee record (new subsection, before s.9)
rep("\n---\n\n## 9. Checks run (evidence only; no proof step depends on them)",
    "\n### 8.4 Referee record\nOne hostile referee pass (AI), verdict: sound with minor fixes, no major items. The eight minor items (s.7 proof (1): the metric is not a rescaled fibrewise metric entrywise, replaced by the entrywise and pairing scalings; s.7 proof (4): the a0 = 0 regularity of GW_{g,0} argued by the residue cancellation and the proper T_3-fixed locus; Lemma 5.7: the 5-torus factor is (e1 + e2), and the d = 1 case is proved by direct localisation; s.6.3: \"the e^{2g-2} term of\"; Lemma 5.2: the coefficient clause deleted; s.3: the first sentence of the proof and the alternative route through LSP (d); Lemma 5.1, s.3 and Lemma 5.5: regularity at both (unramified) points of C over Q = 0 and Q = oo; Remark 6.1: the referee's exact checks) are applied; each item that changes mathematics was re-derived by hand before being applied, and none changed a statement of s.0.3. A fresh direct-localisation test at d = 1, g <= 2 by that referee confirmed e1 e2 GW_{1,1} = P_1 and e1 e2 GW_{2,1} = P_2 on T_CY and the [BS] display (line 1694) off T_CY (evidence only; the g = 1 part is now also a hand derivation, Lemma 5.7). Nothing in this file has been read by a human.\n\n---\n\n## 9. Checks run (evidence only; no proof step depends on them)")

assert src != orig
out = src.replace("\n", "\r\n") if crlf else src
io.open(path, "w", encoding="utf-8", newline="").write(out)
print("applied; chars before/after:", len(orig), len(src))
