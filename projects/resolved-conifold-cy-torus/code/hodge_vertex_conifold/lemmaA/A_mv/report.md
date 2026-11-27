# Route A report (round 1, 2026-10-31)

(Saved by the main session from the workflow result wf_40196105-891; subagents cannot write report files.)

## claims
[
   {
      "depends_on" : "Liu-Liu-Zhou JDG 65 (2003) Thm 1; Okounkov-Pandharipande GT 8 (2004)",
      "evidence" : "a1_check.pl, runs/a1new_*.out: genus 1 (l<=3), genus 2 (l<=3), genus 3 (l<=2), genus 4 (l=1) as full identities in m and sigma; genus 3 l=3 (parts<=6) and genus 4 l=2 (parts<=6) partially; 0 mismatches in 2070 comparisons",
      "statement" : "(A1) The LLZ Marino-Vafa formula, in the real form [x^(2g-2+l)][p_mu] log sum_nu e^{(tau+1/2)kappa_nu x/2} prod_box (2sinh(hx/2))^-1 s_nu(p) = (-1)^g sigma^(l-1) Pi_mu(tau) I_g(mu)/|Aut mu| (x = sqrt(-1) lambda), with data P_g = c^-(3g-3+2l)(-sigma)^-g I_g, matches lemmaA/data/mv_g<g>_l<l> exactly; no correction of the published formula is needed.",
      "status" : "computer-verified"
   },
   {
      "evidence" : "string+dilaton expansion of s(lam_1-1-psi) and Laguerre resummation (report section 2.1); checked against q_g2/q_g3/q_g4 tables, hvc_multiprobe (r=2), (6.1) (r=1) and the degree bound: 287 checks, 0 failures (runs/a2_*.out)",
      "statement" : "(A2) For h>=2: [s^m]Q_h(s;R) = sum_{a<=m} (M_a/a!)(-1)^(m-a) L_{m-a}^{(2h-3+r+a)}(-(J+a)) and Xi_h(R) = sum_{a=0}^{N} ((-1)^a/a!) L_{N-a}^{(2h-2+r+a)}(-J-a) M_a with M_a = <Omega lam_R lam_1^a>, N = 3h-3+r; equivalently Q_h(s;R) = (1-u)^(2h-2+r) e^{-uJ} G_R(u e^{-u}), u = s/(1+s).",
      "status" : "proved"
   },
   {
      "evidence" : "substitution of (1.1) into Cor. 2.2 (report 2.3)",
      "statement" : "Generating-function form: sum_h hbar^(h-1) Q_h(s;P) = log sum_nu e^{-u|nu|} q'^{tau kappa_nu/2} s_nu(q'^rho) s_nu(P~) + const, q' = e^{x/(1+s)}, P~_k = P_k/(sigma x Pi_k), P_1 = s, hbar = sigma x^2; Xi is its genus-by-genus value at s = -1. Lemma A <=> each sum_h hbar^(h-1) Xi_h(R) is a Laurent polynomial in hbar of degree <= floor(D(R)/2).",
      "status" : "proved"
   },
   {
      "evidence" : "lambda_g formula + Lagrange resummation; matches q_g2_r1 (7/5760)(j^2+s+3sj+3s^2); Laguerre transform of b_h(J+a)^(N'+a) reproduces b_h(J-1)...(J-N') (216 cases)",
      "statement" : "lambda_g component for all s: <lambda_h K_s lam_R> = b_h N'! Res_z e^{Jz} dz/(z - s(e^z-1-z))^(N'+1), N' = 2h-3+r.",
      "status" : "proved"
   },
   {
      "depends_on" : "MV cut-and-join (LLZ)",
      "evidence" : "MV cut-and-join restricted to p = Z delta_1 gives (g-1+a+sigma d/dsigma)M_g(1^a) = a(a-1)M_g(2,1^(a-2)) (checked exactly g = 2,3,4); with Cor. 2.2 this is s^2 Q_h(s;(2)) = (1+s)[(h-1)(1+2s+2s^2)Q_h(s;0) + s(1+s)^2 Q_h' + sigma dQ_h/dsigma], evaluate at s = -1; the a = 0 term uses <Omega>_g ~ sigma^(1-g) (CY degree argument)",
      "statement" : "Theorem A: Xi_h(2) = 0 for every h >= 2 and every tau (all Hodge components, not only k = 0, 1, h).",
      "status" : "proved-modulo-literature"
   },
   {
      "depends_on" : "MV cut-and-join; genus-0 Gamma formula and genus-1 value (proved, report 3.1)",
      "evidence" : "labelled cut-and-join at the dilaton point (tau-derivative term drops with a factor (1+s)); checked exactly outside the range: g=2 R=(2),(3),(2,2); g=3 R=(4),(5),(3,2),(3,3),(4,2) at tau=2/7, 3/2 (runs/relation_*.out)",
      "statement" : "Theorem B (2-insertion relation): for all g >= 2 and R, Xi_g(R+2) - sum gam_j Xi_g(R_i^+) + sum gam_jk Xi_g(R_ik) + (sigma/(2 sigma')) sum al be (Pi_al Pi_be/Pi_j)[Xi_{g-1}(al,be,R-j) + sum Xi_{g1}(al,R_I) Xi_{g2}(be,R_J)] = 0, gam_j = (j+1)Pi_{j+1}/(Pi_j sigma'), gam_jk = (j+k)Pi_{j+k}/(Pi_jPi_k sigma').",
      "status" : "proved-modulo-literature"
   },
   {
      "depends_on" : "Theorem B",
      "evidence" : "residues at sigma = -1/4 of Theorem B for R = T - t + (t-1); induction on (h, D, r)",
      "statement" : "Theorem C: at tau = -1/2, Xi_h(T) = 0 for every h and every multiset T of odd parts >= 3 with D(T) <= 2h-3 (TR-free version of the 'odd parts' statement of proof_status 1.6).",
      "status" : "proved-modulo-literature"
   },
   {
      "depends_on" : "Theorem D, Theorem A",
      "evidence" : "two-leaf KP relation (3, c-2) at the simple zero of Pi_3, induction on c from Xi(1) = Xi(2) = 0",
      "statement" : "Theorem C': assuming Lemma A in genera < G, Xi_G(c) = 0 at tau = -1/3 for every 3 <= c <= 2G-2 not divisible by 3.",
      "status" : "proved-modulo-literature"
   },
   {
      "depends_on" : "Sato / Date-Jimbo-Kashiwara-Miwa; Orlov-Scherbin content-product tau functions. NOT in the brief's explicit list (flagged): integrable hierarchy, not TR",
      "evidence" : "fermionic form <0|exp(sum t_k alpha_k) e^{tau x F_2} Gamma_-(q^rho)|0>; checked: first KP equation exactly at p = 0 and p_1-coefficient, x^-5..x^9 (t_kp.pl)",
      "statement" : "Fact: the MV generating function R = sum_nu e^{tau x kappa_nu/2} s_nu(e^{x rho}) s_nu(p) is a KP tau function in t_k = p_k/k.",
      "status" : "proved-modulo-literature"
   },
   {
      "depends_on" : "Fact 3.5; Takasaki-Takebe dispersionless Hirota (Rev. Math. Phys. 7 (1995)); MV",
      "evidence" : "order-x^(2G) part of the KP hierarchy = linearised dispersionless Hirota + intermediate-genus terms; uniform order rule and excess homogeneity at the dilaton point. Validated with ALL terms: full DJKM Hirota equations at the dilaton point vanish exactly for genus 1 (40 equations, weight <= 8) and genus 2 (14 equations, weight <= 6); perturbing one convention breaks all 14; dispersionless Hirota checked on genus-0 MV data (405 relations); closed two-leaf form checked algebraically",
      "statement" : "Theorem D (KP at the dilaton point): if Lemma A holds in genera < G, then for every R and al, be >= 2 with excess 2G-D(R)-al-be >= 1: sum_{a<=al,b<=be} sum_{B subset R} Pi_aPi_b Xi_G(R_B+{a,b}) E_0(al-a,be-b;R-R_B) = (J+al+be-1) Pi_{al+be-1} Xi_G(R+(al+be-1)), E_0 from exp(sum z^-a w^-b sigma Pi_aPi_b Xi_0(.+{a,b})), Xi_0(T) = Gamma(J)/Gamma(J-r+3). For R = empty: sum Pi_aPi_b Xi_G(a,b) z^-a w^-b = (V(z)-V(w))/(Y(w)-Y(z)), V = sum c Pi_c Xi_G(c) z^-c, Y = z - sigma sum Pi_k z^-k/(k+1).",
      "status" : "proved-modulo-literature"
   },
   {
      "depends_on" : "Theorem D",
      "evidence" : "induction on r using Theorem D with the genus-0 corrections vanishing at lower levels; at sigma = 0 it reproduces the lambda_g formula",
      "statement" : "Corollary E (reduction to one leaf): if Lemma A holds in genera < G and at genus G for levels < D, then Xi_G(T) = (J-1)...(J-r+1) Pi_{D+1}/prod_{t in T} Pi_t Xi_G(D+1) for T at level D. Hence Lemma A(G) follows from Lemma A in lower genera and the one-leaf case r = 1 at genus G.",
      "status" : "proved-modulo-literature"
   },
   {
      "depends_on" : "Theorem D (MV + KP + dispersionless Hirota)",
      "evidence" : "full system nullity 0 for G = 5..10 (a3_kp_nullity.pl G<=8, a3_kp_nullity2.pl G<=10; exact over Q at G = 5 and 6); two-leaf subsystem nullity 0 for every G = 7..20 (mod p; exact over Q at G = 7); rigorous lifting argument for the modular rank; independent implementations agree for G = 3..8",
      "statement" : "Theorem F: for every 5 <= G <= 20, Lemma A in all genera < G implies Lemma A in genus G (full column rank of the KP system of Theorem D + polynomiality in sigma + lambda_g).",
      "status" : "computer-verified"
   },
   {
      "depends_on" : "Theorem F; genus-3/4 exact data; KP input (flagged)",
      "evidence" : "h = 2 proved (proof_M_residue T7); h = 3, 4 exact Hodge data of the main session (proof_status 1.9 and Addendum); h = 5..20 by Theorem F inductively",
      "statement" : "Lemma A (M_irr(h)) holds for every h <= 20, for all tau and all multisets R.",
      "status" : "computer-verified"
   },
   {
      "evidence" : "verified for G = 7..20 (nullity 0); fails for G = 4, 5, 6 (nullity 11, 11, 7)",
      "statement" : "Conjecture TLR(G): the two-leaf KP system at genus G (unknown one-leaf polynomials sigma^G Xi_G(c), c = 3..2G-2, of degree <= G with zero constant term; constraints: every coefficient [z^-a w^-b](V(z)-V(w))/(Y(w)-Y(z)), a,b >= 2, a+b <= 2G-1, equals Pi_aPi_b times such a polynomial over sigma^G) has only the zero solution for all G >= 7. Together with Theorem F for G = 5, 6 and the data for h <= 4 it implies Lemma A for all h.",
      "status" : "conjectural"
   },
   {
      "evidence" : "runs/nullity_cj_only.out, runs/nullity_kp_variants.out",
      "statement" : "The cut-and-join relations alone (Theorem B + lambda_g + polynomiality) do not determine Lemma A: nullity 3, 9, 24 at G = 3, 4, 5; with KP + cut-and-join + the proved Hurwitz j=3 component, genus 3 and 4 still have nullity 2 and 1 (they rest on the exact data).",
      "status" : "computer-verified"
   },
   {
      "evidence" : "a3_sorder.pl, runs/sorder_g2.out, runs/a3_sorder_g3.out",
      "statement" : "Q_h(s;R) has a simple zero at s = -1 at every irreducible point of genus 2 and 3 (no higher-order vanishing).",
      "status" : "computer-verified"
   },
   {
      "evidence" : "a3_xisigma.pl on q_g2_r2 (runs/xisigma_g2.out)",
      "statement" : "The genus-2 two-leaf component Xi^(1)_2(j1,j2) vanishes on the lattice range but not on the lines J = 1, 2, 3 (e.g. -1/3072 at (1/2,5/2)): the vanishing is a lattice phenomenon.",
      "status" : "computer-verified"
   }
]


## key_formulas
[
   "MV (real form): [x^(2g-2+l)][p_mu] log sum_nu e^{(tau+1/2)kappa_nu x/2} prod_box (2 sinh(h x/2))^-1 s_nu(p) = (-1)^g sigma^(l-1) Pi_mu(tau) I_g(mu)/|Aut mu|, Pi_mu = prod_i prod_{a=1}^{mu_i-1}(mu_i tau+a)/(mu_i-1)!, x = sqrt(-1) lambda",
   "data: P_g(m;x) = c^-(3g-3+2l) (-sigma)^-g I_g(m), c = 1+x = 1/(1+tau), sigma = -x/c^2",
   "[s^m] Q_h(s;R) = sum_{a=0}^m (M_a/a!) (-1)^(m-a) L_{m-a}^{(2h-3+r+a)}(-(J+a)),  M_a = <Omega lam_R lam_1^a>",
   "Xi_h(R) = sum_{a=0}^{N} ((-1)^a/a!) L_{N-a}^{(2h-2+r+a)}(-J-a) M_a,  N = 3h-3+r",
   "Q_h(s;R) = (1-u)^(2h-2+r) e^{-uJ} G_R(u e^{-u}),  u = s/(1+s),  G_R(y) = sum_a y^a M_a/a!",
   "sum_h hbar^(h-1) Q_h(s;P) = log sum_nu e^{-u|nu|} q'^{tau kappa_nu/2} s_nu(q'^rho) s_nu(P~),  q' = e^{x/(1+s)}, P~_k = P_k/(sigma x Pi_k), P_1 = s",
   "<lambda_h K_s lam_R> = b_h N'! Res_z e^{Jz} dz/(z - s(e^z-1-z))^(N'+1),  N' = 2h-3+r",
   "(g-1+a+sigma d/dsigma) M_g(1^a) = a(a-1) M_g(2,1^(a-2))   (MV cut-and-join on the unit-leg line)",
   "s^2 Q_h(s;(2)) = (1+s)[(h-1)(1+2s+2s^2)Q_h(s;0) + s(1+s)^2 dQ_h(s;0)/ds + sigma dQ_h(s;0)/dsigma]  =>  Xi_h(2) = 0",
   "genus-0 dilaton values: Xi_0(T) = Gamma(J)/Gamma(J-r+3);  string: Xi_g(R+1^k) = J(J+1)...(J+k-1) Xi_g(R)",
   "KP at the dilaton point: sum_{a<=al,b<=be} sum_{B subset R} Pi_a Pi_b Xi_G(R_B+{a,b}) E_0(al-a,be-b;R-R_B) = (J+al+be-1) Pi_{al+be-1} Xi_G(R+(al+be-1))",
   "two-leaf closed form: sum_{a,b>=1} Pi_a Pi_b Xi_G(a,b) z^-a w^-b = (V(z)-V(w))/(Y(w)-Y(z)),  V = sum_c c Pi_c Xi_G(c) z^-c,  Y = z - sigma sum_k Pi_k z^-k/(k+1)",
   "merge formula (in the range): Xi_G(T) = (J-1)(J-2)...(J-r+1) Pi_{D+1}/prod_{t in T} Pi_t  Xi_G(D+1)"
]


## next_steps
[
   "Prove TLR(G) for all G >= 7 (or full-system rigidity for G >= 5): start from the closed form A = (V(z)-V(w))/(Y(w)-Y(z)); at tau = -k/a the conditions say the a-section of A vanishes; combine the pointwise propagation (values for a not dividing c determined by Xi(3..a-1)) with interpolation of the top polynomial v_{2G-2} at the ~G^2 points sigma = -k(a-k)/a^2; an induction in G (nullity drops to 0 exactly at G = 7) is another possibility.",
   "Extend exact (over Q) certificates: full system G = 7, 8 and two-leaf G = 8..12 with a3_kp_nullity2_exact.pl (G = 6 full took 18 min).",
   "Feed 'Lemma A for h <= 20' into leaf_proof/leaf_structure_proof.md s.13: the all-genus regularity theorem then holds through genus 20 (conditional on accepting the KP input).",
   "Decide with the user whether the KP hierarchy of the MV vertex is an admissible input (it is not TR); if yes, record Theorems D-F as part of the no-TR proof.",
   "Optional: use the two-leaf closed form to attack (6.1) outside the range (r = 1 closed formula) and the first non-vanishing layer (excess 0), where the KP relations become inhomogeneous with explicit lower-genus sources."
]


## obstructions
[
   "All genera: the only unproved step for Lemma A at every h is Conjecture TLR(G) (G >= 7): an explicit finite linear rigidity statement for the two-leaf KP system (one-leaf polynomials of degree <= G in sigma; poles of Pi_a Pi_b at the fractional framings tau = -k/a must cancel). Verified for 7 <= G <= 20 and false for G <= 6; no structural proof found: at a fixed framing the relations only propagate values upward in the leaf size (they force vanishing only at tau = -1/2 and -1/3, Theorems C, C'); the rigidity is global (polynomiality across ~G^2 special framings), and subsystems with only small smallest parts are not rigid.",
   "Genus 3 and 4 are not closed by the MV/KP relations (nullity 3 and 2 with KP + cut-and-join, 2 and 1 after adding the proved Hurwitz j = 3 component); the induction from genus 5 on relies on the main session's exact genus-3/4 data.",
   "The cut-and-join equation alone cannot prove Lemma A: at the dilaton point its tau-derivative term is multiplied by (1+s) and drops, and the only join creating a new part from unit legs is (1,1) -> 2, so no relation produces parts >= 3 (nullity 3, 9, 24 at G = 3, 4, 5). The KP hierarchy supplies the missing merge relations.",
   "Level-by-level induction fails: the merge formula forces Z(D) zeros on the one-leaf polynomial at level D, Z(D) = 0,1,2,5,6,11,... < G for small D; only the coupled system (genus-0 terms E_0 linking levels) is rigid.",
   "Input outside the brief's explicit list: the KP hierarchy of the MV vertex (Sato/DJKM) and its dispersionless limit (Takasaki-Takebe). Not topological recursion, but flagged; if excluded, Theorems D, E, F and the h <= 20 statement become conditional on it.",
   "Formula (6.1) (r = 1, all tau, outside the range) is not proved; it is no longer needed for Lemma A up to genus 20."
]


## report
NOT WRITTEN AS A FILE: the subagent harness refused to create A_mv\report.md ("Subagents should return findings as text, not write report files"). The complete report is the 'summary' field of this output. All scripts and raw outputs are in G:\My Drive\Projects\Old\refGW\hodge_vertex_conifold\lemmaA\A_mv\ (outputs in A_mv\runs\).

## requested_checks
[
   "Mathematica, exact over Q: build the KP equations of Theorem D at G = 5 symbolically in tau (clear denominators) and compute NullSpace of the polynomial ansatz sigma^5 Xi_5(T) = sum_{k=1}^5 c_{T,k} sigma^k; expected: only the zero solution (independent implementation of a3_kp_nullity.pl / a3_kp_nullity2_exact.pl).",
   "Exact Hodge computation (hvc_lib.m, Mumford's formula) of a few genus-5 irreducible points at one or two rational values of tau, e.g. Xi_5(3), Xi_5(4), Xi_5(3,3,3), Xi_5(2,2,2,2,2,2,2): Theorem F predicts 0 (a test independent of KP).",
   "Generate the KP Hirota equations to weight 8-10 (bilinear identity) and check them at the dilaton point at genus 2 and 3 with all terms, outside the range (A_mv/a3_hirota_check.pl does genus 1 to weight 8 and genus 2 to weight 6).",
   "Referee review of the inputs of Theorem D: Fact 3.5 (MV vertex = KP tau function), the dispersionless Hirota form exp(D(z)D(w)F) = 1 - (d_1D(z)F - d_1D(w)F)/(z-w), the uniform-order rule (3.1) with the genus-0 Gamma formula, and excess homogeneity."
]


## route
A: Marino-Vafa formula (Schur side, cut-and-join, and the KP hierarchy of the MV vertex evaluated at the dilaton point)

## summary
ROUTE A REPORT (Marino-Vafa), 2026-10-31. Folder HVC\lemmaA\A_mv (HVC = G:\My Drive\Projects\Old\refGW\hodge_vertex_conifold); outputs in A_mv\runs. Exact arithmetic throughout (Perl Math::BigRat/BigInt; ranks also mod p = 2^31-1 with a rigorous lifting argument, and exactly over Q at G = 5, 6, 7). No wolframscript.

THE RULE. No topological recursion of any kind was used: no spectral curve, no Eynard-Orantin invariant, no TR recursion kernel, no BKMP/remodelling, no TR dictionary, no result proved only mod P2, not even as a cross-check. Inputs: the MV formula (Liu-Liu-Zhou, J. Diff. Geom. 65 (2003) 289-340, Thm 1; independent proof Okounkov-Pandharipande, Geom. Topol. 8 (2004) 675-699) and its cut-and-join equation; string/dilaton equations; lambda_g formula; Mumford relation; and (sections 3.5-3.6) the KP hierarchy satisfied by the MV generating function. FLAG FOR THE USER: KP is an integrable hierarchy (Sato / Date-Jimbo-Kashiwara-Miwa), not topological recursion, and it is a consequence of the Schur expansion of the MV formula; but it is not in the brief's explicit list. If it is not admitted, Theorems D, E, F below become conditional on it; everything else is unaffected.

Notation (c = 1): sigma = tau(1+tau), Omega = Lambda(1)Lambda(tau)Lambda(-1-tau)/(-sigma)^h, lam_j = 1/(1-j psi), Tt = psi^2/(1-psi), K_s = exp(sum((1+s)^k-1)kappa_k/k). For R = (j_1..j_r): J = sum j_i, D = J - r, N = 3h-3+r, N' = 2h-3+r, excess e(h,R) = 2h-2-D(R); Lemma A: Xi_h(R) = 0 if all j_i >= 2 and e >= 1. Pi_j(tau) = prod_{a=1}^{j-1}(j tau+a)/(j-1)! (Pi_1 = 1, Pi_2 = 2tau+1 = sigma'). M_h(mu) = <Omega prod lam_{mu_i}>_h, I_h = (-sigma)^h M_h. Q_h(s;R) = <Omega K_s lam_R> = sum_k s^k/k! <Omega lam_R Tt^k> (Q(0) = MV integral, Q(-1) = Xi).

=== 0. RESULTS ===
A1 MV formula fixed exactly; it matches every data table (genus 1-4) with no correction.
A2 finite Laguerre formula for Q and Xi in terms of MV integrals with unit legs; Lambert form; generating-function form (Kaehler-twisted framed vertex).
Thm A: Xi_h(2) = 0 for EVERY h >= 2 and EVERY tau (new; before only components k = 0, 1, h).
Thm B: "2-insertion" relation for all (g,R) from the cut-and-join at the dilaton point.
Thm C: at tau = -1/2, Xi_h(T) = 0 for every multiset T of odd parts in the range, all h (TR-free version of proof_status 1.6). Thm C': at tau = -1/3 (sigma = -2/9) one-leaf Xi_h(c) = 0 for 3 not dividing c (new).
Thm D: KP relations at the dilaton point (all genera). Cor E: Lemma A(G) follows from Lemma A in genera < G plus the ONE-LEAF case at genus G (merge formula); the multi-leaf problem is reduced to r = 1.
Thm F (computer-assisted, rigorous): for 5 <= G <= 20, Lemma A in genera < G implies Lemma A in genus G. With Lemma A for h <= 4 (h = 2 proved, h = 3, 4 exact Hodge data), LEMMA A HOLDS FOR ALL h <= 20.
Remaining gap (all h): one explicit algebraic rigidity statement (Conjecture TLR(G), section 3.7), verified for 7 <= G <= 20.

=== 1. (A1) THE MV FORMULA EXACTLY ===
LLZ Thm 1: C_{g,mu}(tau) = -(sqrt(-1))^l/|Aut mu| [tau(tau+1)]^(l-1) prod_i [prod_{a=1}^{mu_i-1}(mu_i tau+a)/(mu_i-1)!] Int_{Mbar_{g,l}} Lambda(1)Lambda(-tau-1)Lambda(tau)/prod(1-mu_i psi_i), with Int_{Mbar_{0,1}} 1/(1-mu psi) = 1/mu^2, Int_{Mbar_{0,2}} = 1/(mu_1+mu_2); C(lambda;tau;p) = sum_mu sum_g lambda^(2g-2+l) C_{g,mu} p_mu = log R, R = sum_nu e^{sqrt(-1)(tau+1/2) kappa_nu lambda/2} V_nu(lambda) s_nu(p) = sum_mu sum_{|nu|=|mu|} chi_nu(mu)/z_mu e^{...} V_nu p_mu, V_nu = prod_{a<b}[nu_a-nu_b+b-a]/[b-a] prod_i prod_{v=1}^{nu_i} 1/[v-i+l(nu)] = prod_boxes 1/[h], [m] = q^(m/2)-q^(-m/2), q = e^{sqrt(-1)lambda}, kappa_nu = sum nu_i(nu_i-2i+1). Equivalently R = sum_nu q^{tau kappa_nu/2} s_nu(q^rho) s_nu(p).
Real form (x = sqrt(-1) lambda, all sqrt(-1) powers fixed: lambda^(2g-2+l)(-(sqrt(-1))^l) = (-1)^g x^(2g-2+l)):
  [x^(2g-2+l)][p_mu] log sum_nu e^{(tau+1/2)kappa_nu x/2} prod_box (2sinh(hx/2))^-1 s_nu(p) = (-1)^g sigma^(l-1) Pi_mu I_g(mu)/|Aut mu| = sigma^(g+l-1) Pi_mu M_g(mu)/|Aut mu|. (1.1)
Data conversion (w = (c,-x,-1) = c(1,tau,-1-tau), c = 1+x = 1/(1+tau), sigma = -x/c^2): P_g(m;x) = c^-(3g-3+2l) (-sigma)^-g I_g(m). (1.2) No correction of the published formula is needed. Genus 0: I_0(mu) = |mu|^(l-3) for all l >= 1. One-part closed form (checked): sum_g (-1)^g I_g((d)) x^2g = d^-2 S(dx) prod_{a=1}^{d-1} S(ax)/S((d tau+a)x), S(y) = (y/2)/sinh(y/2).
Implementation (mvlib.pm, a1_check.pl): Murnaghan-Nakayama on beta-numbers, hook lengths, W_nu as integer series over a common denominator, [p_mu] log R by the cumulant formula over labelled parts. 7 values of tau outside [-1,0]; I_g is a polynomial of degree <= g in sigma, so the comparison is an identity in sigma; when all parts up to 3g-3+l+1 are covered the data polynomial is determined (identity in m). RESULTS (runs/a1new_*.out), all exact, 0 mismatches: g1 l1 (parts<=5, 35/35), g1 l2 (105/105), g1 l3 (245/245), g2 l1 (<=8, 56/56), g2 l2 (<=6, 147/147), g2 l3 (<=7, 588/588), g3 l1 (<=9, 63/63), g3 l2 (<=9, 315/315) [all of these identities in m], g3 l3 (<=6, 392/392, partial), g4 l1 (<=11, 77/77, identity), g4 l2 (<=6, 147/147, partial). (The first data versions had unevaluated WK symbols; I evaluated them by string/dilaton from the classical base table and got the same matches; all checks were rerun on the 12:01 tables.)

=== 2. (A2) DILATON LEAVES AS UNIT LEGS ===
Prop 2.1 [proved]. h >= 2 (or h = 1, r >= 1), M_a := <Omega lam_R lam_1^a>_{h,r+a}. Then [s^m]Q_h(s;R) = sum_{a=0}^m (M_a/a!)(-1)^(m-a) L_{m-a}^{(2h-3+r+a)}(-(J+a)); Q_h has degree <= N in s; Xi_h(R) = Q_h(-1;R) = sum_{a=0}^N ((-1)^a/a!) L_{N-a}^{(2h-2+r+a)}(-J-a) M_a (Laguerre L_n^(al)(x) = sum_i (-1)^i C(n+al,n-i)x^i/i!). Same for any class pulled back from Mbar_h. Proof: s Tt = s(lam_1 - 1 - psi); expand; remove psi-leaves by dilaton (factor (2h-2+r+a+b)...(2h-3+r+a+b+c)), then constants by string (S lam_j = j lam_j, factor (J+a)^b); resum: Q = sum_a s^a/a! (1+s)^-(2h-2+r+a) e^{-s(J+a)/(1+s)} M_a; Laguerre generating function at t = -s; sum_{n<=K} L_n^(al) = L_K^(al+1); degree bound from 2k <= 3h-3+r+k. QED.
Cor 2.2 [proved]. u = s/(1+s), G_R(y) = sum_a y^a M_a/a!: Q_h(s;R) = (1-u)^(2h-2+r) e^{-uJ} G_R(u e^{-u}); so G_R(ue^{-u}) = e^{uJ}P_R(u)/(1-u)^(5h-5+2r), deg P_R <= N, and Xi_h(R) = (-1)^N [u^N]P_R: Lemma A <=> deg P_R <= N-1 (s = -1 is u = infinity).
Cor 2.3 [proved]. G^(h)_R(y) = |Aut R|/(Pi_R sigma^(h+r-1)) [x^(2h-2+r)][p_R] C(x;tau;p)|_{p_1 = y/(sigma x)}; all R and genera at once (hbar = sigma x^2): sum_h hbar^(h-1) Q_h(s;P) = log sum_nu e^{-u|nu|} q'^{tau kappa_nu/2} s_nu(q'^rho) s_nu(P~) + (unstable constant), q' = e^{x/(1+s)}, P~_k = P_k/(sigma x Pi_k), P_1 := s. Xi = this vertex at s = -1, genus by genus. Reformulation: Lemma A <=> for every R, sum_h hbar^(h-1) Xi_h(R) is a Laurent polynomial in hbar of degree <= floor(D(R)/2).
2.4 [proved] lambda_g part for all s: <lambda_h K_s lam_R> = b_h N'! Res_z e^{Jz} dz/(z - s(e^z-1-z))^(N'+1); s = -1 gives b_h(J-1)...(J-N'). The Laguerre transform of the lambda_g input b_h(J+a)^(N'+a) reproduces it (t_lambdag_laguerre.pl, 216 cases).
2.5 Verification (a2_check.pl; runs/a2_*.out): degree bound (m = N+1, N+2), (2.1a) vs (2.1b), q-tables (all s-coefficients), (6.1) for r = 1, hvc_multiprobe for r = 2: genus 2 r=1 (32 checks), r=2 (48), genus 3 r=1 (100), r=2 (80), genus 4 r=1 (24, two tau), r=2 (3): 0 failures.

=== 3. (A3) ATTACK ===
3.1 Dilaton-point calculus [proved]. For (g,R') with g >= 2, or g = 1 and r' >= 1: G^(g)_{R'}(ue^{-u}) = (1+s)^(2g-2+r') e^{uJ'} Q_g(s;R'); d/dy = (1+s)^3 e^u d/ds gives d^k G/dy^k = e^{u(J'+k)}(1+s)^(2g-2+r'+k)[Xi_g(R'+1^k) + O(1+s)], Xi_g(R'+1^k) = J'(J'+1)...(J'+k-1)Xi_g(R') (string). Uniform order rule: the order is 2g-2+(number of parts incl. unit parts) in all cases, also: genus 0, Lemma 3.1: leading value Xi_0(T) = Gamma(J)/Gamma(J-r+3) for every r >= 1 (unit parts included) [proof: A_m(x) = sum_a (x+a)^(a+m)y^a/a! = e^{xu}(1+s)^(m+1)Q_m, Q_m = xQ_{m-1} + msQ_{m-1} + s(1+s)Q'_{m-1}, so Q_m(-1) = Gamma(x)/Gamma(x-m)]; genus 1 without parts: dG^(1)_empty/dy = e^u(1+s)(s-1/sigma)/24, value Xi_1(1) = -(1+1/sigma)/24. Genus 0/1 values have excess <= -2 / <= 0.
3.2 Cut-and-join on the unit-leg line. MV cut-and-join: dC/dtau = x[(1/2)sum(i+j)p_ip_j d_{i+j}C + (1/2)sum ij p_{i+j}(d_id_jC + d_iC d_jC)]. Prop 3.2 [proved mod MV]: (g-1+a+sigma d/dsigma) M_g(1^a) = a(a-1) M_g(2,1^(a-2)), g >= 2 (on p = Z delta_1 only p_1^2 d_2 survives; a = 0 by the CY degree argument <Omega>_g ~ (w1w2w3)^(1-g) ~ sigma^(1-g)). Checked exactly g = 2 (a<=6), 3 (a<=7), 4 (a<=6).
THEOREM A [proved mod MV]: Xi_h(2) = 0 for all h >= 2, all tau. Proof: insert Cor 2.2 into Prop 3.2: s^2 Q_h(s;(2)) = (1+s)[(h-1)(1+2s+2s^2)Q_h(s;empty) + s(1+s)^2 dQ_h(s;empty)/ds + sigma dQ_h(s;empty)/dsigma], an identity of polynomials; set s = -1. QED.
3.3 THEOREM B [proved mod MV; checked]. For all g >= 2, R (parts >= 2), gam_j = (j+1)Pi_{j+1}/(Pi_j sigma'), gam_jk = (j+k)Pi_{j+k}/(Pi_jPi_k sigma'):
Xi_g(R+2) - sum_i gam_{j_i} Xi_g(R-j_i+(j_i+1)) + sum_{i<k} gam_{j_ij_k} Xi_g(R-j_i-j_k+(j_i+j_k)) + (sigma/(2sigma')) sum_i sum_{al+be=j_i} al be (Pi_al Pi_be/Pi_{j_i})[Xi_{g-1}(al,be,R-j_i) + sum_{g1+g2=g} sum_{I+J=R-j_i} Xi_{g1}(al,R_I)Xi_{g2}(be,R_J)] = 0 (conventions of 3.1). Derivation: labelled cut-and-join for sigma^(g+l-1)Pi_mu M_g(mu) with mu = R+1^a, summed with y^a/a!, at the dilaton point the tau-derivative term acquires a factor (1+s) and drops out, all other terms have the same order. Checked exactly, outside the range (all terms non-zero): g=2: R=(2),(3),(2,2); g=3: R=(4),(5),(3,2),(3,3),(4,2) at tau=2/7 and 3/2 (runs/relation_*.out). Inside the range it determines every Xi with a part 2 from Xi's with fewer 2's. Alone it is insufficient: with lambda_g it leaves nullity 3, 9, 24 at G = 3, 4, 5 (runs/nullity_cj_only.out).
3.4 THEOREM C [proved mod MV]: for all h and every T of odd parts >= 3 with D(T) <= 2h-3: Xi_h(T)|_{tau=-1/2} = 0. Proof: induction on (h, D, r); for t in T apply (3.3) to R = T - t + (t-1) and take the residue at sigma = -1/4 (Pi_m vanishes simply iff m even, sigma' = 0): gam_{t-1}Xi(T) has a simple pole with non-zero residue, the other polar terms are all-odd with fewer parts or lower genus/products, all regular terms drop. THEOREM C' [proved mod MV+KP, given Lemma A in lower genera]: at tau = -1/3, Xi_G(c) = 0 for all 3 <= c <= 2G-2 with 3 not dividing c (two-leaf KP relation (3,c-2) at the simple zero of Pi_3, induction on c from Xi(1) = Xi(2) = 0). Generally at tau = -k/a the two-leaf relations determine the values Xi_G(c)(sigma_0) (a not dividing c) from Xi_G(3..a-1)(sigma_0).
3.5 KP AT THE DILATON POINT. Fact 3.5: R = <0|exp(sum t_k alpha_k) e^{tau x F_2} Gamma_-(q^rho)|0> (t_k = p_k/k) with e^{tau x F_2}, Gamma_- exponentials of fermion bilinears, so R is a KP tau function (Sato; Date-Jimbo-Kashiwara-Miwa 1983; = Orlov-Scherbin content-product tau function, TMP 128 (2001)). Checked: first KP equation F_1111 + 6F_11^2 + 3F_22 - 4F_13 = 0 for C at p = 0 and its p_1-coefficient, x^-5..x^9 (t_kp.pl). Dispersionless Hirota (Takasaki-Takebe, Rev. Math. Phys. 7 (1995) 743): exp(D(z)D(w)F) = 1 - (d_1D(z)F - d_1D(w)F)/(z-w), D(z) = sum z^-n d_{t_n}/n; checked on the genus-0 MV potential at P = 0 (a+b <= 10) and on the genus-0 dilaton values (360 relations). With P_k = xp_k, C = sum x^(2g-2)F_g(P): at order x^(2G) the F_G-linear terms of the KP hierarchy are exactly the linearisation of the dispersionless Hirota equation at F_0; every other term contains a factor of genus 1..G-1. At the dilaton point every term has the same (1+s)-order and THE SAME EXCESS e = 2G - D(R) - (al+be) (excess homogeneity); since genus-0/1 values have e <= 0, in a relation with e >= 1 every intermediate-genus factor that occurs is itself in the range.
THEOREM D [proved mod MV + Fact 3.5 + TT]: G >= 2, Lemma A in all genera < G. For every R (parts >= 2) and al, be >= 2 with e >= 1:
sum_{a<=al,b<=be} sum_{B subset R} Pi_aPi_b Xi_G(R_B+{a,b}) E_0(al-a,be-b;R-R_B) = (J+al+be-1) Pi_{al+be-1} Xi_G(R+(al+be-1)),
E_0(m,n;S) = coefficient of z^-m w^-n p_S in exp(sum_{a,b>=1} z^-a w^-b sigma Pi_aPi_b Xi_0(.+{a,b})), unit parts of genus-G arguments by string, Xi_G(1^k) = 0. For R = empty this is the closed form (verified, t_twoleaf_closed.pl): sum_{a,b>=1} Pi_aPi_b Xi_G(a,b) z^-a w^-b = (V(z)-V(w))/(Y(w)-Y(z)), V(z) = sum_c c Pi_c Xi_G(c) z^-c, Y(z) = z - sigma sum_k Pi_k z^-k/(k+1). Validation with ALL terms (independent of the dispersionless formula): the full DJKM Hirota equations at the dilaton point vanish exactly for genus 1 (all 40 equations up to weight 8, R = empty,(2),(3)) and genus 2 (all 14 up to weight 6, R = empty,(2)) (a3_hirota_check.pl); deliberately perturbing one convention makes all 14 fail (sensitivity); KP1, KP2 checked separately at genus 2.
COROLLARY E [proved mod D]: if moreover Lemma A holds at genus G for all levels < D, then for every T at level D: Xi_G(T) = (J-1)(J-2)...(J-r+1) Pi_{D+1}/prod_{t in T}Pi_t Xi_G(D+1). Hence LEMMA A(G) <= LEMMA A(<G) + LEMMA A(G, r=1). Poles of the prefactor force zeros of sigma^G Xi_G(D+1): Z(D) = 0,1,2,5,6,11,13,19,22,29 for D = 1..10 (~D^2); a level closes alone iff Z(D) >= G, never for small D.
3.6 THEOREM F (computer-assisted, rigorous). Unknowns: sigma^G Xi_G(T) = sum_{k=1}^G c_{T,k} sigma^k (degree <= G by proof_status 1.1(e); no sigma^0 term by lambda_g 1.3). Equations: Theorem D for all (R,al,be) in the range, at many sample tau. Full column rank implies no non-zero rational solution (mod-p argument: a coprime integral solution would reduce to a non-zero solution mod p), hence Lemma A(G). Results (nullity): full system G=2:2(KP)/0(+CJ); G=3: 6 (KP), 3 (KP+CJ), 2 (+P4 Hurwitz j=3), 0 only with Xi_3(3)=Xi_3(4)=0 added; G=4: 6, 2, 1, 0 with Hurwitz(4) or Xi_4(3)=0 added; G = 5, 6, 7, 8 (all T as unknowns: 220, 576, 1358, 2976 unknowns) and G = 9, 10 (reduced: one-leaf unknowns 135, 170): NULLITY 0 with KP ALONE. Exact over Q (no modular step): G = 5 full (35/35), G = 6 full (54/54), G = 7 two-leaf (77/77). Without the lambda_g input: nullity 4 at G=5, 0 at G=6. TWO-LEAF SUBSYSTEM (R = empty only + polynomiality of two-leaf T): nullity 11, 11, 7 at G = 4, 5, 6 and 0 for EVERY G = 7, ..., 20 (unknowns 77..740); three-leaf: 4 at G=5, 0 at G=6. Two independent implementations (a3_kp_nullity.pl, a3_kp_nullity2.pl) agree for G = 3..8; second prime/seed agrees. CONCLUSION: Lemma A(<G) => Lemma A(G) for 5 <= G <= 20; with h <= 4 known, Lemma A holds for every h <= 20.
3.7 WHERE IT STOPS. (1) All genera: Conjecture TLR(G): if v_c in sigma^-G sigma Q[sigma]_{<=G-1} (c = 3..2G-2; v_1 = v_2 = 0), V = sum c Pi_c v_c z^-c, and every coefficient [z^-a w^-b](V(z)-V(w))/(Y(w)-Y(z)) with a,b >= 2, a+b <= 2G-1 equals Pi_aPi_b times an element of sigma^-G sigma Q[sigma]_{<=G-1}, then V = 0. TLR(G) for all G >= 7 (with the computed G = 5, 6 and the data for G <= 4) implies Lemma A for all h. Verified 7 <= G <= 20; false for G <= 6. No structural proof found: pointwise at tau = -k/a the relations only propagate values upward in c (vanishing at a = 2, 3: Thms C, C'), the rigidity is global (polynomiality across ~G^2 special framings against 2G^2 unknowns); restricting to leaves with smallest part 2 or 3 leaves large nullity (so all pairs (a,b) are needed). (2) Genus 3, 4 are not closed by these relations; they rest on the exact data. (3) Cut-and-join alone is insufficient: at the dilaton point it only produces the (1,1)->2 join and its tau-derivative term drops; KP supplies the merges. (4) (6.1) not proved (no longer needed for Lemma A up to genus 20).
3.8 Other observations (computer-verified): Q_h(s;R) has a simple zero at s = -1 at every irreducible point of genus 2, 3; Xi^(1)_2(j1,j2) vanishes on the lattice range but not on the lines J = const (e.g. -1/3072 at (1/2,5/2)): lattice phenomenon, consistent with the integer (Pi_j, merge) structure of MV/KP.

=== 4. FILES (A_mv) ===
mvlib.pm (WK numbers, characters, hooks, MV Schur side mv_conn/mv_unit, Laguerre, InputForm parser); t_basic.pl, t_unit.pl, t_kp.pl, t_lambdag_laguerre.pl, t_twoleaf_closed.pl; a1_check.pl; a2_check.pl; a3_sorder.pl, a3_components.pl, a3_qsigma.pl, a3_xisigma.pl, a3_cj_line.pl, a3_relation.pl, a3_nullity.pl (cut-and-join system), a3_kp_check.pl, a3_kp_genus0_test.pl, a3_hirota_check.pl (+ a3_hirota_sens1/2.pl sensitivity), a3_kp_nullity.pl (all T), a3_kp_nullity_nolg.pl, a3_kp_nullity2.pl (reduced; args G [samples] [rmax]), a3_kp_nullity2_exact.pl (over Q), a3_kp_nullity2_sub.pl, a3_zcount.pl. Outputs: runs\ (a1new_*, a2_*, cj_line_*, relation_*, hirota_*, kp_genus0_test, nullity_*, kpn_*, kpn2_*, kpn2_twoleaf_g7..g20, exact_g5_full, exact_g6_full, exact_g7_twoleaf, zcount, sorder_*, components_*, xisigma_*, qsigma_*). Examples: perl a3_kp_nullity2.pl 9; perl a3_kp_nullity2.pl 12 0 2; perl a3_kp_nullity.pl 4 C H=3; perl a1_check.pl 3 2 9; perl a3_hirota_check.pl 1 2/7 8 ";2;3".

