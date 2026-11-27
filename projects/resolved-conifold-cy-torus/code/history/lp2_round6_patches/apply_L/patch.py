# Moved on 2026-11-27 from localp2/hodge_vertex/lp2_explore/round6_scripts/apply_L/patch.py (local P2 project) into the resolved-conifold
# project, because it edits only files of this project. One-off patch script, already applied: a record, not to be rerun.
import os
P = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "hodge_vertex_conifold", "lemmaA", "LEMMA_A_PROOF.md")
t = open(P, encoding="ascii", newline="").read()
E = []
def r(a, b): E.append((a, b))

# ---- header, provenance
r("(round 5, plan item P1 of lemmaA/round2/round2_status.md s.6).",
  "(round 5, plan item P1 of lemmaA/round2/round2_status.md s.6). Edited in round 6 with the referee fixes of s.11.4 (round-5 text: LEMMA_A_PROOF_round5.md).")
r("The assembled text has not yet had its own referee (open point O1, s.10).",
  "The assembled text was refereed in round 6 under two lenses, logic and literature (L_logic) and independent computation (L_comp): both \"sound with minor fixes\", no major item; the fixes are applied (s.11.4, and O1 in s.10.2). (MV) itself stays cited and not re-proved (O3).")

# ---- s.0
r("**Corollary L1 (specialisation).** [P] Write", "**Corollary L1 (specialisation).** [P modulo (MV)] Write")
r("Hence Xi_h(R; tau_0) = 0 at every complex tau_0 with sigma(tau_0) != 0, including the diagonal",
  "Hence Xi_h(R; tau_0) = 0 at every tau_0 in a field of characteristic 0 (in particular at every complex tau_0) with sigma(tau_0) != 0, including the diagonal")
r("**Corollary L2 (parts >= 1).** [P] The same vanishing holds for every non-empty",
  "**Corollary L2 (parts >= 1).** [P modulo (MV)] The same vanishing, in Q(tau) and at every such tau_0, holds for every non-empty")
r("**Corollary L3 (general scale).** [P] Let E,", "**Corollary L3 (general scale).** [scale law P; vanishing P modulo (MV)] Let E,")

# ---- s.3
r("expand the k-th term of Q_g:", "expand the m-th term of Q_g:")

# ---- s.6: multiset T -> V, level D -> d, |alpha| -> l(alpha)
r("on D = 1, ..., 2h-3. Fix G = h and D, and assume:", "on d = 1, ..., 2h-3. Fix G = h and d, and assume:")
r("parts >= 2 and D(R') < D.", "parts >= 2 and D(R') < d.")
r("U non-empty, D(U) < D.", "U non-empty, D(U) < d.")
r("Let T be a non-empty multiset with parts >= 2 and D(T) = D. Pick **any** part k_0 of T and put R = T - k_0, which may be empty; then D(R) = D - k_0 + 1.",
  "Let V be a non-empty multiset with parts >= 2 and D(V) = d. Pick **any** part k_0 of V and put R = V - k_0, which may be empty; then D(R) = d - k_0 + 1.")
r("Since Sum_B D(U_B) = D(R) + k_0 - |alpha|, where |alpha| = number of parts of alpha = Sum_B |B|,",
  "Since Sum_B D(U_B) = D(R) + k_0 - l(alpha), where l(alpha) = number of parts of alpha = Sum_B |B|,")
r("= D(R) + k_0 - Sum_B n_B = D + 1 - 2G <= -2 < 0.", "= D(R) + k_0 - Sum_B n_B = d + 1 - 2G <= -2 < 0.")
r("<= D(R) + k_0 - |alpha| = D + 1 - |alpha|.", "<= D(R) + k_0 - l(alpha) = d + 1 - l(alpha).")
r("- If |alpha| >= 2: D(U) <= D - 1, so", "- If l(alpha) >= 2: D(U) <= d - 1, so")
r("- If |alpha| = 1: alpha = (k_0), the only block is B = {k_0}, k = 1, R_B = R, so U = T.",
  "- If l(alpha) = 1: alpha = (k_0), the only block is B = {k_0}, k = 1, R_B = R, so U = V.")
r("tau delta_G(T) X_G(T) = 0,      delta_G(T) = D + 2 - 2G <= -1.", "tau delta_G(V) X_G(V) = 0,      delta_G(V) = d + 2 - 2G <= -1.")
r("as k_0 >= 2. So X_G(T) = 0.", "as k_0 >= 2. So X_G(V) = 0.")
r("h = 2, D = 1, T = (2), R = empty, k_0 = 2, n = 4.", "h = 2, d = 1, V = (2), R = empty, k_0 = 2, n = 4.")

# ---- s.7.1 - 7.3
r("### 7.1 Corollary L1 [P]", "### 7.1 Corollary L1 [P modulo (MV)]")
r("So the vanishing holds at every tau_0 with sigma(tau_0) != 0, also where the proof's Pi_{k_0} vanishes. QED",
  "So the vanishing holds at every tau_0 in any field of characteristic 0 with sigma(tau_0) != 0, also where the proof's Pi_{k_0} vanishes. QED")
r("### 7.2 Corollary L2 and the parts 1 and 0 [P]", "### 7.2 Corollary L2 [P modulo (MV)] and the parts 1 and 0 [P]")
r("and (S) with Theorem L gives the claim.",
  "and (S) with Theorem L gives the claim. The argument of s.7.1 applies to U: each Xi^(k)_h(U) = 0, so L2 also holds at every such tau_0.")
r("### 7.3 Corollary L3: general scale [P]", "### 7.3 Corollary L3: general scale [scale law P; vanishing P modulo (MV)]")

# ---- s.7.5
r("- The identification of the Taylor coefficients of the graph sum with Xi_h(R). For local P2 this is draft s.7.5(ii) with interface item (h), and it has rank-3 features with no conifold counterpart.",
  "- The identification of the Taylor coefficients of the graph sum with Xi_h(R). For local P2 this is written in LP2 draft s.1.5 (Fpot_h on the neighbourhood of q_0, dilaton-shifted form), s.7.1 (H0) and s.7.5(i)-(ii), tagged [P] there; the s.1.5 extension and (H0) are marked unrefereed in the round-5 draft. For the conifold it is leaf_structure_proof s.13.1(i), which still needs the domain fix of round-5 critic item S3(a) (the form of draft s.1.5 and (H0), under the dictionary c -> E).")
r("- The vertex identity for R empty (the gap): (VI) of the LP2 draft, F4 at n = 0 for the conifold.",
  "- The vertex identity for R empty (the gap): Theorem VI of the LP2 draft (s.7.6, proof s.7.8), equal to F4 at n = 0 for the conifold. It is a separate result with its own inputs (O5), and it uses nothing from Theorem L. (Referee D_VI's simplification of its Step 2 would use (4.2) at k = 2 of s.4 here, which does not depend on Theorem L, so no circularity arises.)")
r("- The frame-B formula identifying the graph sum with the GW potential.",
  "- The frame-B formula identifying the graph sum with the GW potential: Proposition 13.2_P2 (LP2 draft s.8.3, open) and conifold Proposition 13.2 (citation fixes of round-5 critic item S3(c)).")

# ---- s.8.2
r("For R empty and h >= 2 the same holds for a >= 1, plus the a = 0 term L_N^{(2h-2)}(0) <Omega_h>, with <Omega_h> := Int_{Mbar_h} Omega_h (s.8.5).",
  "For R empty and h >= 2, N = 3h-3: Xi_h(empty) = L_N^{(2h-2)}(0) <Omega_h> + sigma^{1-h} Sum_{a=1}^{N} (-1)^a sigma^{-a} L_{N-a}^{(2h-2+a)}(-a) [x^{2h-2+a}] [p_1^a] F, with <Omega_h> := Int_{Mbar_h} Omega_h (s.8.5); the a = 0 term does NOT carry the prefactor sigma^{1-h}.")

# ---- s.8.5
r("The a = 0 coefficient is L_{N+1}^{(2h-3)}(0) = binom(N+2h-2, N+1) != 0. So M_0 = <Omega_h> is a finite combination of MV coefficients M_a, a >= 1.",
  "The a = 0 coefficient is (-1)^{N+1} L_{N+1}^{(2h-3)}(0) = (-1)^{N+1} binom(N+2h-2, N+1) != 0. So <Omega_h>, the a = 0 value of Lemma 3.2 at U = empty, is a finite combination of MV coefficients M_a, a >= 1.")
r("Non-vanishing as an element of Q[tau, 1/sigma]: [C: R2 Cor. 7, sigma^h Xi_h(empty) = -(2h-3)! [l^{2h}] S(e1 l) S(e2 l), e1 = -tau, e2 = 1+tau, S(y) = (y/2)/sinh(y/2); Hodge-side TR, allowed]. This non-vanishing,",
  "Closed form: [C: LP2 draft Theorem VI_tau, s.7.8, TR-free] sigma^{h-1} Xi_h(empty) = beta_{2h-2}(-tau, 1+tau), beta_n(a,b) = (n-1)! [t^n] 1/(4 sinh(at/2) sinh(bt/2)); this is R2 Cor. 7 (Hodge-side TR, cross-check only). Non-vanishing: with S(y) = (y/2)/sinh(y/2) = Sum c_i y^{2i}, sign c_i = (-1)^i, c_i != 0, one has beta_{2h-2}(a,b) = ((2h-3)!/(ab)) Sum_{i+j=h} c_i c_j a^{2i} b^{2j}, whose terms all have sign (-1)^h at real a, b != 0. So Xi_h(empty) != 0 at every real tau with sigma != 0, hence in Q(tau). This non-vanishing,")
r("which has zeros at sigma = -0.804..., -0.363... (complex tau).",
  "which has zeros at sigma = -0.804..., -0.363...; both satisfy sigma < -1/4, so tau is complex there, in agreement with the real non-vanishing.")

# ---- s.9
r("Introduction read directly (author and referees of round 5): C_{g,mu} with -sqrt(-1)^{|mu|+l}/|Aut mu|; the l < 3 extension; eqs. (1), (2), (3); kappa_nu. Theorems 1, 2 are the cut-and-join for R and C, so \"LLZ Thm 1\" for the formula is wrong by arXiv numbering. Prop. 2.2 (tau = 0) agrees. Journal data confirmed on Project Euclid. | Journal text not opened, so \"eq. (1)\" is confirmed for arXiv v3 only (dated 31 Dec 2004, after publication). LLZ's proof not re-checked. |",
  "Introduction read directly in arXiv v3 (author and referees of round 5) and in the JDG print (round 6, EMIS mirror, https://www.maths.tcd.ie/EMIS/journals/NYJM/jdg/p/2003/65-2-4.pdf), same numbering: C_{g,mu} with -sqrt(-1)^{|mu|+l}/|Aut mu| and the l < 3 extension (JDG pp. 290-291); eqs. (1), (2), (3) (p. 291); \"we will call (1) the Marino-Vafa formula\", Theorem 1 with (4), (5) (p. 292); Theorem 2 with (6) (p. 293); (8), (9), Prop. 2.1 (p. 294); kappa_nu. Theorems 1, 2 are the cut-and-join for R and C, so \"LLZ Thm 1\" for the formula is wrong; \"eq. (1)\" is confirmed in JDG and arXiv v3. Prop. 2.2 (tau = 0) agrees. Journal data confirmed on Project Euclid. | LLZ's proof not re-checked. |")
r("No theorem number pinned (e.g. Faber-Pandharipande, Invent. Math. 139 (2000), not opened for one) |", "No theorem number pinned |")
r("| Exponent computation done in place | Macdonald I.3 Ex. 2",
  "| Exponent computation done in place; LLZ eq. (12) (JDG p. 297) states s_nu(1,q,q^2,...) = q^{n(nu)}/prod(1-q^{h}) | Macdonald I.3 Ex. 2")
r("| A_g closed form as quoted in OP s.0.6 | FP not opened |",
  "| A_g closed form as quoted in OP s.0.6. FP, arXiv:math/9810173v1, read: Theorem 4 and s.4.3 give Int lambda_{g-1}^3 = (|B_{2g}|/2g)(|B_{2g-2}|/(2g-2))/(2g-2)! for g >= 2, and that this equals 2 Int lambda_g lambda_{g-1} lambda_{g-2} from c(E)c(E^vee) = 1. Agrees with OP s.0.6 and s.8.5 | Journal numbering of FP not compared |")
r("Web sources read in round 5: arXiv:math/0306434, projecteuclid.org/euclid.jdg/1090511689, arXiv:math/0307209, arXiv:1407.7800.",
  "Web sources read in round 5: arXiv:math/0306434, projecteuclid.org/euclid.jdg/1090511689, arXiv:math/0307209, arXiv:1407.7800. In round 6: the JDG print of LLZ (EMIS mirror above), arXiv:math/9810173v1.")

# ---- s.10.1
r("| (2.2)-(2.4) | P | Prop 2.1, Pi_k != 0 in Q(tau) |",
  "| (2.2) | P mod (MV) | restates (2.1); Pi_k != 0 in Q(tau) |\n| (2.3), (2.4) | P | definition and formal |")
r("| Cor L1, L2, L3 | P mod (MV) | Theorem L; L3 by weight count |",
  "| Cor L1, L2, L3 | P mod (MV) | Theorem L; the scale law of L3 and the identities for parts 1, 0 of s.7.2 are P |")
r("| s.8.5, non-vanishing of Xi_h(empty) | C | R2 Cor. 7, Hodge-side TR; outside Lemma A |",
  "| s.8.5, closed form and non-vanishing of Xi_h(empty) | C | LP2 draft Theorem VI_tau (s.7.8); outside Lemma A |")

# ---- s.10.2
r("- **O1 (refereeing).** This assembled document has not had its own referee. Its two source reports had two referees each, and s.2-s.6 were re-derived by the writer. Notation was renamed during assembly (dictionary in s.1), and doc_check.py (s.12) was written against the final formulas to guard against transcription errors, but that is evidence only.",
  "- **O1 (refereeing): closed in round 6.** The assembled text (frozen copy LEMMA_A_PROOF_round5.md) was refereed under two lenses. L_logic: s.2-s.7 line by line by hand, including every sign, (F1)-(F6) and the induction of s.6, a search for hidden finite-order inputs and circularity, and the literature of s.9 read directly. L_comp: new code written from the document's formulas, plus hand checks of the same steps: (4.2), (5.1), the logic of s.6 term by term, Theorem L and its polynomial form in 1/sigma; these runs are evidence only. Both verdicts: sound with minor fixes, no major item. The confirmed fixes are in s.11.4. Their stated limits: the proof of (MV) was not re-checked (O3); s.8.6 and the runs of s.12 were not repeated; the interfaces listed in s.7.5 were not assessed. The referees had not seen this edited text, and the edits change no mathematics of s.2-s.6.")
r("- **O2 (citation details; none is mathematical).** Journal numbering of LLZ. Macdonald and Stanley",
  "- **O2 (citation details; none is mathematical).** Macdonald and Stanley")
r("- **O5 (outside Lemma A, unchanged by this document).**\n  - Empty profile: the gap formula and F4 / (VI) still rest on R2 Cor. 7 (Hodge-side TR) or on the frame-B route. No TR-free proof of the gap formula exists yet (plan P5). Critic 1 S3 of round 4 is not addressed here.\n  - The identification of Theorem R's Taylor coefficients with Xi_h(R), conifold (Theorem R) and local P2 (interface item (h), restated Dil_0).\n  - Prop. 13.2 and Prop. 13.2_P2 (frame-B), lem:modrel (critic 1 S4), the gap corollary (S5).",
  "- **O5 (outside Lemma A, unchanged by this document).**\n  - Empty profile: the gap formula and (VI) are Theorem VI of the LP2 draft (s.7.8), TR-free, P modulo literature. In the round-5 draft that literature is LLZ eqs. (1), (5) and the vertex-operator form of psi, psi*; referee D_VI of round 6 (sound with minor fixes, no major item; second pass on Theorem VI) gives replacements that reduce it to LLZ eq. (1) at mu = 1^a plus textbook facts and (4.2) at k = 2 here. The draft's s.7.8 is the record of which applies. This closes plan P5 for j = 0, and critic 1 S3 of round 4 is answered there. R2 Cor. 7 is a cross-check only (draft s.7.6 (B-a)).\n  - The identification of Theorem R's Taylor coefficients with Xi_h(R): see s.7.5 (local P2 written in the draft, conifold pending S3(a)).\n  - lem:modrel (critic 1 S4) is closed by Proposition 2.12 of the LP2 draft. Still open: Prop. 13.2_P2, the conifold Prop. 13.2 citation fixes (S3(c)), and hence the gap corollary (S5).")

# ---- s.11
r("Item 2 of LP2 draft s.8.2: its first bullet (2D Toda input) is replaced by s.4, its second bullet (LLZ/OP numbers) is settled against arXiv v3 only (O2), and its third bullet (matching) is closed by s.2.",
  "Item 2 of the round-4 draft's s.8.2 (backup LP2_REGULARITY_PROOF_draft_round4.md): its first bullet (2D Toda input) is replaced by s.4, its second bullet (LLZ/OP numbers) is settled against the JDG print and arXiv v3 (s.9), and its third bullet (matching) is closed by s.2.")
r("and not extended to item (h): s.7.5, O5.", "and not extended to the Taylor identification: s.7.5, O5.")
r("A_mv/report.md line 8, R3 report line 250, round2_status lines 26 and 460, LP2 draft s.7.4 and s.8.",
  "A_mv/report.md line 8, R3 report line 250, round2_status lines 26, 42 and 460.")
r("as printed has 2^l (the text extraction shows it as 2l) where", "as printed has 2^l (JDG p. 294) where")
r("- **LP2 draft.** s.0.3, s.7.4(b) with its input list, s.8.1 label of Lemma A_P2(b) (to \"P mod (MV), LEMMA_A_PROOF.md\", after O1), s.8.2 item 2 (replace by a pointer here), s.8.2 \"Nothing else\" and Appendix.",
  "- **LP2 draft.** s.0.3, s.7.4(b), s.8.2 were edited in round 5, and \"Nothing else\" is gone. Remaining: s.8.1 label of Lemma A_P2(b) to \"P mod (MV), LEMMA_A_PROOF.md\", and the statements in s.0.3, s.7.4 and s.8.2 that this file has had no referee, now that O1 is closed; the last sentence of the s.7.4 bullet \"Not provided by that file\" (see O5); \"JDG numbering not confirmed\" in s.7.4, s.7.8, s.8.2 (see s.9); \"LLZ eq. (13) states it\" in s.7.8 should read eq. (12).")
r("---\n\n## 12. Evidence [E]",
  "### 11.4 Round-6 referee fixes applied here\n\nFrom L_logic and L_comp (scripts: scratchpad lp2_round6/L_logic/, lp2_round6/L_comp/; no mathematics of s.2-s.6 changed):\n- Stale passages synced with the LP2 draft: s.7.5, s.8.5, s.10.1, O5, s.11.1-s.11.3.\n- Xi_h(empty): closed form quoted from Theorem VI_tau, and a written proof of non-vanishing in Q(tau) from the signs of the coefficients of S: s.8.5.\n- Letters in s.6: multiset V, level d, l(alpha), so that T, D and |.| keep one meaning each (s.11.2 claimed this before it was true); \"m-th term\" in Lemma 3.2.\n- Status labels: L1, L2, the vanishing in L3 and (2.2) are P mod (MV), not P: s.0, s.7.1-7.3, s.10.1.\n- Scope of specialisation: any field of characteristic 0, and parts >= 1: s.0, s.7.1, s.7.2.\n- Sign (-1)^{N+1} of the a = 0 coefficient, and <Omega_h> in place of M_0: s.8.5. Placement of the prefactor sigma^{1-h} for R empty: s.8.2.\n- Citations: JDG print of LLZ read, LLZ eq. (12), FP Theorem 4 and s.4.3, E5: s.9, O2, s.11.3.\n\n---\n\n## 12. Evidence [E]")

# ---- s.12 letter consistency
r("| Mechanical run of the s.6 logic over all terms | all T for G <= 5 |", "| Mechanical run of the s.6 logic over all terms | all V for G <= 5 |")
r("h = 4 (r <= 5, all in-range T),", "h = 4 (r <= 5, all in-range V),")

for a, b in E:
    n = t.count(a)
    assert n == 1, (n, a[:70])
    t = t.replace(a, b)
t.encode("ascii")
open(P, "w", encoding="ascii", newline="").write(t)
print("applied", len(E))
