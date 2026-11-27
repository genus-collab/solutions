# Referee report on route D (round 1, 2026-10-31)

(Saved by the main session from the workflow result wf_40196105-891; subagents cannot write report files.)

## confirmed
[
   "Key identity (1-j/c)_{j-1} = (-1)^{j-1} prod_{k<j}(j tau + k), since k - j/c = -(j tau + j - k) and 1/c = 1 + tau. Also (j tau + k)(j tau + j - k) = j^2 sigma + k(j-k). Hence eps_j = C(j tau+j-1, j-1)^2, C_j = (c^2/j) eps_j(sigma) (checked numerically for j <= 9 at 6 values of x), and E(d) = -eps_d(sigma), the localization edge factor.",
   "Cor. 11.3 (top a_{k+1}^(k) = (-1)^{k+1} C_{k+1}; sub-top a_k^(k) = (-1)^k k C_k (1-x)^2/(2c^2)), re-derived from the closed form of Theorem 11.1. ref_closedform.pl rebuilds A_j and rho from scratch as exact t-series and checks them: sum rule; Sum j A_j = -c^2 rho to t^8; rho and every a_j^(k), k <= 5, equal the independent Mathematica fits hvc_leaf_structure_data_kt5.m; top and sub-top to k = 8 (2418 checks, 6 values of x).",
   "Lemma 5.2 (string-flow cancellation). In resummed form, the coefficient of Xi^c(R) in F_h(rho q_t) is G_R = (1 - Y_1/c^2)^{-J_R} prod_{j>=2} Y_j^{m_j}/m_j!. Its t^{D+1} coefficient vanishes identically, its t^D coefficient is prod (-1)^j C_j/|Aut R|, and its t^{D+2} coefficient is generically non-zero. Verified for all 44 irreducible R with D <= 7 at 6 values of x.",
   "Lemma 5.3 (cleanliness), re-derived: a non-zero term needs D_v >= min(k-1, 2h_v-2), so its order is >= k+3-2h (one vertex) or >= 2k+2-2h (two or more vertices), both > k+2-2h for k >= 1.",
   "Theorem D5 and the c = 1 form of T_k: the c-power c^{3-3h} is independent of R. T_1, T_2, T_3 re-derived: eps_2 = 1+4s, eps_3 = (2+9s)^2/4, eps_4 = (3+16s)^2(1+4s)/9.",
   "D5.1 (Xi_h(2) = 0 for all h, all tau), D5.2 (2(2+9s)^2 Xi_h(3) = 3(1+4s)^2 Xi_h(2,2); divisibility, P_(3) = s(1+4s)^2 A_h, the factor 27/8, deg A_h <= h-4 after P4) and D5.3 are correct. Since (H) holds in all genera by Brini-Schuler Prop. 3.1, all three are unconditional.",
   "(H) itself: F_g - F_g(0) = P_g Li_{3-2g}(Q)/(e1 e2) on the project torus for every g (Brini-Schuler Prop. 3.1). Confirmed independently by my GP localization on Z with alpha0 != 0 (ref_loc5.pl, ref_loc5_more.pl): exact agreement up to (g,d) = (4,1), (3,2), (2,3) on three tori of the Proposition's class and on one generic T_CY torus.",
   "Genus-2 test of section 5, reproduced and extended with an evaluator independent of hodge_lib.pl, built from the Mathematica psi_g2 table with string/dilaton reduction. It agrees with the q-tables (8 multisets R, 4 values of x). It reproduces the direct bare-vertex t^0 and t^1 coefficients at x = 2, 1/3, -3/5, 5/2 and -7/3, both with the full closed-form leaf expansion and with the truncated T_2/T_3 forms (the latter includes Xi_2(2,2,2)).",
   "D2.1: the GP formula for the resolved conifold and the GV checks (classical; FP 2000 s.4).",
   "D2.2(a)-(c), re-derived: each graph is a polynomial in sigma of degree g-1+d, and its top coefficient is (-1)^{g-1+d} times the P^1 graph contribution (vertex sign (-1)^{g_v+val_v-1}, edge sign (-1)^{d_e+1}).",
   "D2.3(i)-(v), re-derived: the dilaton-homogeneity plus string-flow lemma, the Laguerre form (Laguerre generating function and partial sums), the Lambert form, and the tree-function form by Lagrange inversion with v = -T(y). The Laguerre form was verified independently at genus 2 for (4), (5), (2,3), (3,3), (2,2,2) at two values of x.",
   "D3.2: the resonance factor -(w1/d)^2/(c - w1/d) at a p_oo end, and T_{w1} = w1 + psi - w1^2/(w1 - psi).",
   "The conditional honest-bundle count (Xi^(k) = 0 only for 2k < r) and sigma^h Omega K = c(T(-log)) e(E^vee x C_tau) e(E x C_{1+tau}).",
   "M_irr(3) for every tau (D4 table). The engine's non-zero genus-3 values agree with the Mathematica q-tables (68/68, ref_engine_g3.pl). The tables give M_irr(3) and M_irr(4) at every point with r <= 2, and the D = 2h-2 controls are non-zero.",
   "check_polar_relations.pl and check_vertex_orders.pl reproduce their recorded outputs exactly.",
   "TR audit: Route D, its scripts and its dependencies (Theorem LC; Theorems A, B, C; proof_status s.1.1-1.5 including P4; Brini-Schuler Prop. 3.1, which uses degeneration and Bryan-Pandharipande) use no topological recursion and no mod-P2 result."
]


## issues
[
   {
      "description" : "This error works in the author's favour. The report says (H) ([t^{2-2h+k}] F_{h,0} = 0 for 1 <= k <= 2h-3) is 'established only through genus 4', so that Theorem D5 and its corollaries are conditional. That is wrong. (H) is a theorem in every genus on exactly the project torus: Brini-Schuler, 'Refined Gromov-Witten invariants' (project file arXiv/main_DMJ.tex, s.3.1, Proposition 3.1 = Theorem A1). For T inside T_CY with alpha_i = eps_j for some i, j, it proves GW^{Z,T}_{d[P1]} = (1/d)(2 sinh(d e1/2))^-1 (2 sinh(d e2/2))^-1 for every d > 0 and every genus. The proof uses Li's degeneration formula, Mumford's relation, the dimension count of BP08 Lem. 6.3 and BP05 Thm 5.1; no TR is involved. The project torus (alpha0 = 0, alpha1 = e1, alpha2 = e2) satisfies alpha1 = e1. Its weight conventions ({alpha0, alpha1, alpha2, -e1, -e2}, CY alpha0+alpha1+alpha2 = e1+e2) are the paper's. Hence F_g - F_g(0) = P_g Li_{3-2g}(Q)/(e1 e2). Since Li_{-n}(e^t) = n!(-t)^{-n-1} + regular, (H) holds for all h. I confirmed this independently with ref_loc5.pl and ref_loc5_more.pl, using my own GP localization on Z with alpha0 != 0. The prediction P_g d^{2g-3}/(e1 e2) is reproduced exactly at (g,d) up to (4,1), (3,2) and (2,3), on three tori of the Proposition's class and on one generic T_CY torus. Independence of alpha0 together with polynomiality carries the result to alpha0 = 0. So D5, D5.1, D5.2 and D5.3 are unconditional, modulo the project framework.",
      "location" : "Summary, 'Precise obstruction' item 4; hypothesis (H) of Theorem D5; claims D5, D5.1, D5.2, D5.3 (depends_on)",
      "severity" : "major",
      "suggested_fix" : "Cite Brini-Schuler Prop. 3.1 for (H). Change the status of Theorem D5, D5.1, D5.2 and D5.3 to proved-modulo-literature (Prop. 3.1, Theorem LC, Theorems A-C, the frame-B Givental-Teleman premise). Delete obstruction 4. Add the same citation to Theorem R / Corollary R (leaf_structure_proof.md s.13.3(5): 'the known polar part of F_g')."
   },
   {
      "description" : "(V) says the genus-h single-vertex graph equals (-x c^3/t^2)^{h-1} F_h(rho q_t). This is only its colour-1 (singular-idempotent) term. The colour-2 term Delta_2^{h-1}(...) is also present; it is regular at t = 0 by F1 and F3 of Theorem LC, so it does not affect the polar orders. The proof of D5 needs this sentence.",
      "location" : "Section 5, input (V) and proof of Theorem D5",
      "severity" : "minor",
      "suggested_fix" : "State (V) for the colour-1 term. Add that the colour-2 single-vertex term is analytic at t = 0 (F1: Delta_2 analytic; F3: Dil_2 analytic)."
   },
   {
      "description" : "D5.1 has a much shorter proof that needs neither Lemma 5.2 nor Lemma 5.3. Graphs with an edge have pole order <= 2h-4 (Theorem A), and the colour-2 vertex is regular. Also [t^1](rho q_t) = ((1-x)^2/2)(lam_2 - lam_1) (L1 with Cor. 11.3, or s.11.4 by hand), so [t^1] F_h(rho q_t) = ((1-x)^2/2) Xi^c_h(2) by Theorem C. Hence [t^{3-2h}] F_{h,0} = (-x c^3)^{h-1} ((1-x)^2/2) Xi^c_h(2) exactly. Inside the framework, '(H) at order 1' and 'Xi_h(2) = 0' are therefore equivalent: the whole content of D5.1 is the refined GV closed form.",
      "location" : "Corollary D5.1 (k = 1)",
      "severity" : "minor",
      "suggested_fix" : "Give the three-line proof and state the equivalence explicitly. Keep Lemmas 5.2 and 5.3 for k >= 3, where they are really needed: there the D = k-1 symbols Xi(3) and Xi(2,2) are not known to vanish for h >= 4."
   },
   {
      "description" : "Precisely, D3.1 is a degree count: after string reduction, no class of degree < vdim is integrated over compact moduli. Its stated conclusion ('no vanishing relation beyond string/dilaton/divisor') is informal. It does not exclude relations coming from weight-independence and polynomiality, which the report itself uses in D2.2(c). It also does not cover degree-0 (non-compact) contributions, where equivariant residues of low degree need not vanish.",
      "location" : "Proposition D3.1 (status 'proved')",
      "severity" : "minor",
      "suggested_fix" : "Restate D3.1 as the precise degree-count statement for d >= 1, or relabel it as a heuristic obstruction."
   },
   {
      "description" : "D2.4 is labelled proved-modulo-literature, but the exact statement and normalisation are not given (the constant C_mu, the weight at the fixed point, the ELSV #Aut(mu) convention). It is also not checked against the data, as the brief requires for that status. Its use is only as an obstruction remark, so nothing downstream depends on it.",
      "location" : "Proposition D2.4 (Graber-Vakil)",
      "severity" : "minor",
      "suggested_fix" : "Quote Graber-Vakil (Compositio 135 (2003) 25-36) with its normalisation and check one case against psi_g2/mv tables, or downgrade it to a remark."
   },
   {
      "description" : "Three qualifications. (1) 'The root s = 1 of Psi_R is simple' is a computer observation (psi_poly.out, h = 2, 3, the listed R), not a proved statement. (2) In (iv), 'Xi = lim_{v->oo}' must be read through the rational function Psi_R(v/(1+v)): the identity holds only as a formal series, and Phi_R has a finite radius of convergence. (3) The Hurwitz identity with H/(2h-2+J+r+2p)! is correct only with the ELSV convention that includes 1/#Aut(mu), #Aut(R u 1^p) = #Aut(R) p!; otherwise a 1/p! is missing.",
      "location" : "Proposition D2.3, remarks after (ii) and (iv); the ELSV identity",
      "severity" : "minor",
      "suggested_fix" : "Label (1) computer-verified, with its range. Phrase (2) via the rational function. State the ELSV convention used in (3)."
   },
   {
      "description" : "In check_data.pl, the Perl engine's dilaton-leaf functional was compared with the Mathematica q-tables at genus 3 only at points where Xi = 0. The genus-3 M_irr table rests on it. This gap is now closed: ref_engine_g3.pl compares every non-zero genus-3 value printed in t_genus3.out and check_genus3_alltau.out (r <= 2, 7 values of tau) with the q-tables, and all 68 checks agree. ref_g2.pl corroborates the engine at genus 2 for r = 3.",
      "location" : "check_data.pl / D4 engine validation",
      "severity" : "minor",
      "suggested_fix" : "Add a non-zero comparison point to check_data.pl, for example Xi_3(5), Xi_3(2,4) or Xi_3(3,3)."
   },
   {
      "description" : "'With the CY torus the p0-vertex is exactly the Lemma A class' is too strong. The localization vertex carries the Hodge class Omega and the leaves d lam_d, but no K-twist (no dilaton leaves). The report itself says so in D3.",
      "location" : "Proposition D2.1 wording",
      "severity" : "cosmetic",
      "suggested_fix" : "Write 'carries the Hodge class Omega of Lemma A (without the K-twist)'."
   },
   {
      "description" : "loc_twisted.pl has no .out file. D3.2's statement that the sum with gamma_0 insertions is regular, with p0-vertices equal to Lemma A vertices, is therefore unrecorded. By the code's construction it is regular at c = w1, since T-points at p0 ends sit at psi = -w1/d.",
      "location" : "loc_twisted.pl / D3.2",
      "severity" : "cosmetic",
      "suggested_fix" : "Run loc_twisted.pl on one case and store the output, or mark the remark as unverified."
   },
   {
      "description" : "The claims list gives D5, D5.1, D5.2 and D5.3 the status 'proved', while the summary says 'proved modulo (H)'. D5.2b also rests on P4, whose Boissy/BCGGM citation is still incomplete; the author already flags this.",
      "location" : "Claims list versus summary text",
      "severity" : "cosmetic",
      "suggested_fix" : "After the first fix, use proved-modulo-literature consistently, and keep the P4 dependency explicit for D5.2b."
   }
]


## notes
REFEREE REPORT, Route D (virtual localization) of the Lemma A attack, 2026-10-31.
The harness forbids report files from subagents, so referee.md was not written; this text is the report. Scripts and outputs are in G:\My Drive\Projects\Old\refGW\hodge_vertex_conifold\lemmaA\D_localization\referee\ : ref_lib.pl, ref_leafseries.pl, ref_closedform.pl/.out, ref_g2.pl/.out, ref_engine_g3.pl/.out, ref_loc5.pl/.out, ref_loc5_core.pl, ref_loc5_more.pl/.out. All are Perl + Math::BigRat; no wolframscript and no TR of any kind were used.

VERDICT: sound. The section-5 result is stronger than the author claims: its only open input (H) is a theorem in every genus. The output of Route D is:
- Xi_h(2) = 0 for every h >= 2 and every tau, with no TR. Previously only the components k = 0, 1 and h were known.
- In every genus, one relation at D = 2: 2(2+9s)^2 Xi_h(3) = 3(1+4s)^2 Xi_h(2,2). This closes genus 3 without a computer, modulo P4; for h >= 4 it leaves h-3 scalars after P4.
- In every genus, one relation at D = 3.
Lemma A remains open from D = 2 on for h >= 5 (D = 2 at h = 4 is covered by the exact data).

1. (H) IS PROVED FOR ALL h. This is the main correction, in the author's favour.
- Source: Brini-Schuler, 'Refined Gromov-Witten invariants' (project file arXiv/main_DMJ.tex), s.3.1, Prop. 3.1 = Theorem A1. For T inside T_CY with alpha_i = eps_j for some i, j, and every d > 0: GW^{Z,T}_{d[P1]} = (1/d)(2 sinh(d e1/2))^-1 (2 sinh(d e2/2))^-1, in all genera.
- Its proof: Li's degeneration formula; Mumford's relation; a dimension count (BP08 Lem. 6.3); Bryan-Pandharipande BP05 Thm 5.1. No TR.
- Applicability: the project torus (README: alpha0 = 0, alpha1 = e1, alpha2 = e2; weights {alpha0, alpha1, alpha2, -e1, -e2}; CY alpha0+alpha1+alpha2 = e1+e2) uses the paper's conventions and satisfies alpha1 = e1.
- Consequence: F_g - F_g(0) = P_g Li_{3-2g}(Q)/(e1 e2), with P_g = [u^{2g}] S(e1 u) S(e2 u) and S(y) = (y/2)/sinh(y/2). Since Li_{-n}(e^t) = n!(-t)^{-n-1} + regular with n = 2h-3, (H) holds for every h.
- Independent check (ref_loc5*.pl). I wrote my own GP localization on Z with a base-acting torus:
  - five weights per fixed point;
  - edge factor = the two O(-1) H^1 terms, divided by the moving H^0(TP1) and by (-e1)(-e2);
  - vertex = e(T_pZ)^{val-1} Int prod_5 Lambda^vee(w) prod 1/(t/d - psi), with the usual unstable-vertex conventions.
  It reproduces P_g d^{2g-3}/(e1 e2) exactly:
  - at (g,d) = (0,1), (0,2), (1,1), (1,2), (2,1), (2,2), (3,1), on the tori a1 = e1 with alpha0 = 1 and alpha0 = -3/7, on a2 = e2 with alpha0 = 2, and on a generic T_CY torus;
  - at (2,3), (3,2), (4,1) on two of them.
  Independence of alpha0 together with polynomiality (compact moduli for d >= 1) carries this to alpha0 = 0.
- So obstruction 4 is refuted, and D5 and its corollaries are proved-modulo-literature. The same citation should be added to Theorem R / Corollary R ('known polar part', leaf_structure_proof.md s.13.3(5)).

2. SECTION 5, STEP BY STEP.
(a) KEY IDENTITY and C_j = (c^2/j) eps_j: correct (re-derived). Also E(d) = -eps_d.
(b) Cor. 11.3: re-derived. From A_j = (-1)^j C_j t^{j-1} (zeta/t)^j 2F1(jx/c, jx/c; j; tau), the [t^1] coefficient is j(1/c^2 - 1/2) + j x^2/c^2 = j(1-x)^2/(2c^2). ref_closedform.pl builds A_j and rho from scratch at x = 2, 1/3, -5/7, 3/2, 7/4, -2/9 and checks:
- the sum rule and Sum j A_j = -c^2 rho, to t^8;
- rho to t^5 and every a_j^(k), k <= 5, against the independent Mathematica fits hvc_leaf_structure_data_kt5.m;
- top and sub-top coefficients to k = 8.
All agree, which also confirms the closed form independently.
(c) Lemma 5.2: correct. Cleanest form: the coefficient of Xi^c(R) in F_h(rho q_t) is G_R = (1 - Y_1/c^2)^{-J} prod_{j>=2} Y_j^{m_j}/m_j!, where Y_1 = A_1 + c^2 and Y_j = A_j. The first factor is the genus-independent string resummation of the lam_1 parts. With Y_j = (-1)^j C_j t^{j-1}(1 + j s1 t/c + ...) and Y_1 = -c s1 t + ..., the t^{D+1} coefficient is (J - J)(s1/c) = 0. Verified for all 44 irreducible R with D <= 7 at 6 values of x. The t^{D+2} coefficient is generically non-zero, so the induction in D5 is needed.
(d) Lemma 5.3: correct (re-derived). A non-zero term needs D_v >= min(k-1, 2h_v-2) at each singular vertex of genus >= 2. Its order is then >= k+3-2h for one such vertex and >= 2k+2-2h for two or more, both > k+2-2h.
(e) D5: correct. The c-power c^{3-3h} is independent of R. The proof must also mention that the colour-2 single-vertex term is regular (F1, F3).
(f) Independent genus-2 test (ref_g2.pl). The evaluator uses only the Mathematica psi_g2 table with string/dilaton reduction, not hodge_lib.pl.
- It agrees with the q-tables for (2), (3), (4), (5), (2,2), (2,3), (3,3), (2,4) at 4 values of x.
- It reproduces the direct genus-2 bare-vertex t^0 and t^1 coefficients (leaf_structure.md) at x = 2, 1/3, -3/5, 5/2, -7/3, both with the full closed-form leaf expansion and with the truncated T_2/T_3 forms.

3. SIMPLER PROOF OF D5.1. Graphs with an edge have pole order <= 2h-4 (Theorem A), and the colour-2 vertex is regular. Also [t^1] F_h(rho q_t) = ((1-x)^2/2) Xi^c_h(2) by L1 and Theorem C. Hence [t^{3-2h}] F_{h,0} = (-x c^3)^{h-1} ((1-x)^2/2) Xi^c_h(2), and Prop. 3.1 gives Xi_h(2) = 0. Within the framework the two statements are equivalent.

4. D2.
- D2.1: GP/FP formulas correct. Wording: the p0-vertex carries Omega but no K-twist.
- D2.2: re-derived, including the sigma-degree g-1+d and the sign (-1)^{g-1+d}.
- D2.3(i)-(v): re-derived (homogeneity plus string flow; T = 1 + psi - lam_1; Laguerre generating function; Lagrange inversion with v = -T(y)). The Laguerre form was verified independently at genus 2 for (4), (5), (2,3), (3,3), (2,2,2). Qualifications:
  - 'root s = 1 simple' is a computer observation, not a proof;
  - the 'lim' in (iv) must be read via the rational function Psi_R(v/(1+v));
  - the ELSV identity needs the #Aut(mu) convention.
- D2.4: no exact statement or normalisation, and no check against data.

5. D3 AND THE OBSTRUCTIONS.
- D3.1: really a degree count. Its 'no relation' conclusion is informal: it does not exclude weight-independence/polynomiality relations (which D2.2(c) uses) or degree-0 residues. Downgrade or restate.
- D3.2: computations correct. loc_twisted.pl has no .out file.
- Obstructions (2)-(3): a correct count (one scalar-function relation per (h,k) against p(k) unknowns).
- Obstruction (1): a fair heuristic.
- Obstruction (4): refuted (point 1).

6. ENGINE.
- The author's section-5 scripts reproduce their outputs exactly.
- check_data.pl compared the genus-3 dilaton functional only at zero values. ref_engine_g3.pl compares the non-zero genus-3 values with the q-tables: 68/68 agree.
- The tables also confirm M_irr(3) and M_irr(4) at every point with r <= 2, with non-zero controls at D = 2h-2.

7. TR AUDIT: clean. Dependencies:
- Theorem LC (one analytic input: Stirling);
- Theorems A, B, C;
- the frame-B Givental-Teleman premise;
- proof_status s.1.1-1.5, including P4, whose citation is incomplete (used only in D5.2b);
- Brini-Schuler Prop. 3.1.

8. SUGGESTION (speculative). Other tori in the class of Prop. 3.1 (alpha0 != 0) also satisfy (H). If the frame-B conifold analysis can be redone there and gives leaf coefficients not proportional to eps_j(sigma), it would give independent extra relations at each D. That could close D = 2.

## refuted
[
   "Obstruction 4 / the hypothesis status of Theorem D5. The claim that (H) ('the polar part of F_{h,0} at Q = 1 is gap_h t^{2-2h} only') is established only through genus 4, so that D5 and its corollaries are conditional, is false. (H) holds for every h by Brini-Schuler, 'Refined Gromov-Witten invariants', Prop. 3.1 (= Theorem A1), whose hypothesis alpha_i = eps_j is met by the project torus (alpha1 = eps1). This was confirmed independently by 5-fold localization through genus 4."
]


## route
D_localization (Route D: virtual localization relations / conifold-point polar relations T_k)

## verdict
sound

