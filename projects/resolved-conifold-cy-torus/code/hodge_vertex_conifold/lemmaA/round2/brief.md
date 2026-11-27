# Lemma A, round 2: the all-genus step (shared brief, 2026-10-31)

Project folder: G:\My Drive\Projects\Old\refGW\hodge_vertex_conifold (HVC below). Round-2 work goes in HVC\lemmaA\round2\<route>\.
Round-1 brief (background, definitions, data): HVC\lemmaA\brief.md. Round-1 reports: HVC\lemmaA\<route>\report.md and referee.md
for <route> = A_mv, B_hurwitz, C_intersection, D_localization, E_formulas (Route E's referee is still running).

## 0. Rules (they override everything)
 - TOPOLOGICAL RECURSION. User clarification, 2026-10-31, replacing the earlier blanket ban: "I don't want you to use TR applied
   to the cohft of the resolved conifold times a plane. But it's ok to use TR considerations applied to the Hodge CohFT."
   FORBIDDEN: any use of TR (Eynard-Orantin recursion, BKMP/remodelling, loop equations, mirror-curve methods) for the resolved
   conifold x C^2 theory itself (its free energies, its polar part, its graph sum). The input (H), the conifold closed form, is the
   Brini-Schuler degeneration theorem and must stay that.
   ALLOWED: TR considerations for the Hodge CohFT, i.e. the triple Hodge class Omega with CY weights (1, tau, -1-tau) = the framed
   C^3 vertex. This includes:
    - TR on the curve C_tau: x = v^tau/(1+v)^(tau+1), y = -v;
    - BKMP for C^3 (L. Chen; J. Zhou; Eynard-Orantin);
    - the DOSS identification of TR with Givental's formula;
    - the dictionary P2 (Lemma B) and Theorem 3.1 of proof_M_hurwitz.md (tau = 1);
    - Lemma D, and the "mod P2" results of proof_status_mixed_scale.md s.1.6-1.8 (conditional on P2 until P2 is written out);
    - the TR scripts as cross-checks: review_checks/status_tr*.pl, referee1_tr*.pl, proof_M_scripts/trtau.pl, hvc_tr_tau.wls.
   Lemma A is a statement about the Hodge CohFT, so TR on C_tau may be used for it (see s.5). Cite Faber-Pandharipande, Ann. of
   Math. 157 (2003) (localization), for the lambda_g formula.
 - NEW, user ruling 2026-10-31: the KP hierarchy of the Marino-Vafa generating function IS an allowed input. So are: KP and 2D Toda
   for content-product (hypergeometric) tau functions in the GL(infinity) orbit of the vacuum (Takasaki 1984; Orlov-Scherbin, Theor.
   Math. Phys. 128 (2001)); Okounkov's Toda equation for Hurwitz numbers (math/0004128); the dispersionless Hirota equations
   (Takasaki-Takebe, Rev. Math. Phys. 7 (1995)); the fermionic / E-operator formalism of Okounkov-Pandharipande; the MV formula and
   its cut-and-join equation (Liu-Liu-Zhou, J. Diff. Geom. 65 (2003); Okounkov-Pandharipande, Geom. Topol. 8 (2004)); ELSV;
   character formulas; virtual localization (Graber-Pandharipande; Faber-Pandharipande, Invent. Math. 139 (2000)); GRR, Mumford's
   relations, Chiodo's formula; Farkas-Pandharipande strata dimensions; Lagrange inversion; exact algebra. Genus-0 (dispersionless)
   data and changes of variables are fine as algebra; setting up a recursion on a spectral curve is not.
 - Working rules. No wolframscript (licence: about two kernels, used by the main session): write any Mathematica script into your
   folder and list it under requested_checks. Exact checks in Perl with Math::BigRat (no python on this machine). Write files with
   the Write tool (long shell heredocs have failed here); ASCII only. You cannot kill processes: test every script on small cases
   first, wrap long runs in "timeout <seconds> perl ...", and run at most 4 heavy processes at a time (28 cores shared with other
   agents). Do not edit files outside your own folder (copy what you need). Subagents cannot write report files (report.md and the
   like are refused): put your complete report, in markdown, in the report_text field of your structured answer; the main session
   saves it. Scripts and outputs may be written to your folder.
 - Every claim carries a status: proved (complete argument written out), proved-modulo-literature (exact statement and
   normalisation of the cited result checked against the data), computer-verified (state the range), conjectural, or refuted.

## 1. Statement and conventions (c = 1)
Weights w = (1, tau, -1-tau), sigma = tau(1+tau); Omega = Lambda^vee(1) Lambda^vee(tau) Lambda^vee(-1-tau)/(-sigma)^h;
lam_j(psi) = 1/(1 - j psi); dilaton leaves T(psi) = -psi^2/(1 - psi). For a multiset R = {j_1..j_r}:
    Xi_h(R) = Sum_{m>=0} (1/m!) Int_{Mbar_{h,r+m}} Omega prod_i lam_{j_i}(psi_i) prod_{k=1}^m T(psi_{r+k}).
Level D(R) = Sum (j_i - 1), J = Sum j_i, r = #R, Pi_j(tau) = prod_{a=1}^{j-1} (j tau + a)/(j-1)! (Pi_1 = 1, Pi_2 = 2 tau + 1).
Xi depends on tau only through sigma; Pi_j(-1-tau) = (-1)^(j-1) Pi_j(tau).
**Lemma A = M_irr(h): Xi_h(R) = 0 identically in tau for every R with all parts >= 2 and D(R) <= 2h - 3.**
Why it matters: Corollary R (HVC\leaf_proof\leaf_structure_proof.md s.13, the all-genus regularity theorem for the conifold
graph sum) follows from Lemma A in all genera. (Main-session observation, to be checked by the critic: at genus g it needs Lemma A
only up to genus g - 1, because a stable graph of genus g with an edge has no genus-g vertex and the known polar part of F_g then
fixes the single-vertex graph.)

## 2. Known for all genera (cite these; do not re-derive them)
 (a) lambda_g component (k = h in Xi = Sum_k (-1/sigma)^k Xi^(k)): b_h (J-1)...(J-N), N = 2h-3+r; it vanishes in the range.
 (b) Hurwitz component (k = 0): vanishes for every h and R in the range (Route C, Theorem H; Route B, Theorem B, also on the
     extra layer D = 2h - 1). With (a): in the range sigma^h Xi_h(R) lies in span(sigma, ..., sigma^(h-1)).
 (c) Xi_h(2) = 0 for all h (Route A, Theorem A, by cut-and-join; Route D, Corollary D5.1, by the closed form).
 (d) Route A (A_mv\report.md; referee A_mv\referee.md: sound). Prop. 2.1 (Laguerre formula: Xi as a finite combination of MV
     integrals with unit legs); Cor. 2.2 (u = s/(1+s): G_R(u e^-u) = e^{uJ} P_R(u)/(1-u)^(5h-5+2r), deg P_R <= N = 3h-3+r, and
     Xi_h(R) = (-1)^N [u^N] P_R, so Lemma A <=> deg P_R <= N - 1; the dilaton point s = -1 is u = infinity, the other sheet of
     y = u e^-u); Cor. 2.3 (all genera at once: a Kaehler-twisted framed vertex; Lemma A <=> Sum_h hbar^(h-1) Xi_h(R) is a Laurent
     polynomial in hbar of degree <= floor(D(R)/2)); Theorem B (2-insertion relation from cut-and-join); Theorems C, C' (vanishing at
     tau = -1/2 and tau = -1/3); Fact 3.5 (the MV series is a KP tau function); Theorem D (KP at the dilaton point; it holds given
     Lemma A in genera < G); Corollary E (one-leaf reduction: if Lemma A holds in genera < G and at genus G for all levels < D, then
     for every T at level D: Xi_G(T) = (J-1)(J-2)...(J-r+1) Pi_{D+1} / prod_{t in T} Pi_t * Xi_G(D+1)); Theorem F (exact ranks:
     Lemma A in genera < G implies Lemma A(G) for 5 <= G <= 20, hence Lemma A for h <= 20); Conjecture TLR(G) (s.3.7: the two-leaf
     KP system is rigid), verified for 7 <= G <= 20, false for G <= 6. Scripts: A_mv\mvlib.pm, a3_kp_nullity*.pl (usage in s.4).
 (e) Route D (D_localization\report.md; referee: sound). Theorem D5: the closed form F_g - F_g(0) = P_g Li_{3-2g}(Q)/(eps1 eps2) of
     the conifold (a theorem: Brini-Schuler, "Refined Gromov-Witten invariants", Prop. 3.1 = Thm A1, file G:\My Drive\Projects\Old\
     refGW\arXiv\main_DMJ.tex s.3.1; degeneration + Bryan-Pandharipande, no TR) gives, for 1 <= k <= 2G-3,
         T_k(G) = Sum_{D(R)=k} (1/|Aut R|) prod_i ((-1)^{j_i} Pi_{j_i}^2 / j_i) Xi_G(R) = 0
     given Lemma A for levels <= k-2 in all genera <= G. Lemma 5.2 (the level k-1 terms cancel), Lemma 5.3 (graphs with an edge do
     not contribute given Lemma A in lower genera). Without the reduction hypothesis, the single-vertex-graph coefficient of every
     polar order k <= 2G-3 vanishes given only Lemma A in genera < G; its coefficients come from the leaf closed forms (Theorem CF,
     leaf_structure_proof.md s.0.3 and s.11).
 (f) Main session (lemmaA\parity\README.txt, exact check): substituting Corollary E into T_D gives T_D(G) = Pi_{D+1} kappa_D Xi_G(D+1)
     with kappa_D = [y^D] (tau y + log((e^y - 1)/y)): kappa_1 = tau + 1/2, kappa_D = B_D/(D D!) != 0 for even D, and kappa_D = 0
     for every odd D >= 3. So the closed form settles every even level and gives nothing at odd levels.
 (g) Route E (E_formulas\report.md): exact Perl Hodge engine (WK.pm, Hodge.pm, HodgePt.pm); M_irr(5), M_irr(6) complete (every
     component); M_irr(7) complete too (Route E referee engine, 194 points, all components zero; update of 2026-10-31);
     genus-5 Mathematica cross-checks at the non-special framings tau = 2/7 and 3/2 confirm 16 engine values; the one-leaf formula (6.1) of proof_M_residue.md s.6,
         sigma^h Xi_h(j) = prod_{l=1}^{2h-2} (j - l) [lambda^(2h)] S(e1 lambda)^(1 + j e2) S(e2 lambda)^(1 + j e1),
         e1 = -tau, e2 = 1 + tau, S(y) = (y/2)/sinh(y/2),
     holds identically for h <= 7 (conjectural beyond); conjectural one-leaf universality (SN) and two-leaf form (NR) (report s.0(c),(d)),
     generating-function form s.1.3 (Lemma A <=> Zhat_R(u) = Sum_h sigma^h Xi_h(R) u^{N_h}/N_h! is a polynomial of degree <= J - 1).
     Exact component tables: E_formulas\engine_h<h>_r<r>.txt ((h, r) up to (5, 7) and (6, 3)), mirr_g5.txt, mirr_g6.txt, mirr_g7.txt.
 (h) Route B (B_hurwitz\report.md): ELSV + Okounkov's Toda equation + a T -> infinity leading-order analysis give the recursion (REC)
     for the Hurwitz component and Theorem B: Xi^(0)_h(j) = (j-1)...(j-2h+2) [t^(2h)] S^2 H^j. Referee B (B_hurwitz\referee.md)
     discusses the CY analogue and why it is non-linear.
 (i) Route C (C_intersection\report.md): Theorem H (localized top Chern class on P(ker phi), ker phi = H^0(S_D), dimension bound
     from Farkas-Pandharipande); partial analysis of the k = 1 component (conjectural modulo a plumbing model).

## 3. What remains
By (c), Corollary E, Theorem D5 and (f), with induction on G and, inside genus G, on the level D (D = 1 by (c); even D by (f);
Corollary E then gives every T at level D), Lemma A for all h is equivalent to the following.
 (OL) Odd-level one-leaf vanishing. For every G >= 3 and every odd D with 3 <= D <= 2G - 3: if Lemma A holds in all genera < G and
      Xi_G(T) = 0 for every T with D(T) < D, then Xi_G(D + 1) = 0.
(OL) holds for G <= 20 (Theorem F, exact data for G <= 4, Route E's engine for G = 5, 6), so a proof for all G >= G0 with G0 <= 21
suffices. Each of the following implies Lemma A for all h: (OL) for all large G; TLR(G) for all G >= 7; the one-leaf formula (6.1)
for all h (only its vanishing for even j, under the hypotheses of (OL), is needed); or any direct all-genus argument. You may use the
inductive hypotheses (Lemma A in lower genera, lower levels at genus G) wherever they are stated.

## 4. Code and data you may reuse (read-only; copy into your folder to modify)
 lemmaA\data (exact Hodge tables, genus <= 4; read WARNING_READ_FIRST.txt); E_formulas (engine and component tables); A_mv (MV Schur
 side, characters, Laguerre transform, InputForm parser, KP relation builders and nullity scripts, runs\); B_hurwitz (HB.pm,
 HBconn.pm, Rec.pm); C_intersection (polylib.pm, components.pl); D_localization (hodge_lib.pl, xi_lib.pl, check_polar_relations.pl);
 lemmaA\parity\kappa_parity.pl.

## 5. The TR route for the Hodge CohFT (allowed since the clarification of s.0)
 - proof_M_hurwitz.md s.2: the base curve C_tau (x = w/(1+w)^(tau+1), y = -1/w; with v = 1/w: x = v^tau/(1+v)^(tau+1), y = -v), its
   Laplace leaves xi_0 = Sum_mu C(mu(tau+1)-1, mu-1) x^mu, and the dictionary P2:
       Xi_h(mu; tau) = sigma^(h-1) [x^mu] F^{C_tau}_{h,r} / prod_i C(mu_i(tau+1)-1, mu_i-1), expanded at w = 0 (v = infinity).
 - proof_M_hurwitz.md s.3, Theorem 3.1: at tau = 1 the global involution v -> 1/v (X o sigma = X) and the residue theorem give the
   degree bound for all g, n, hence Lemma A at tau = 1 and tau = -2 for all h and r; at tau = -1/2 (involution v -> -1-v) only odd
   parts. Referee: referee_hurwitz_1.md s.1-3 (line-by-line check), s.4 (a rigorous route to P2: BKMP for C^3, the EO variational
   formula along the REGULAR family Y_e = Y_B + e Y_0, rationality in e, homogeneity and the dilaton equation as e -> infinity).
 - proof_status_mixed_scale.md s.1.6-1.8 and Lemma D (s.2): for integer tau >= 2 there is no global involution; the global linear
   loop equation holds (omega_{g,n} at the sheet near v = infinity = minus the sum over the tau sheets near v = 0); what is missing
   is the analogue of the residue step (the quadratic loop equation summed over pairs of sheets, including mixed pairs, with the
   kernel expanded at v_0 = infinity). At tau = 2, 3 the correlators have low degree at the other zero v = 0 of x.
 - Combination with s.3: at genus G and an odd level D, under the inductive hypotheses, (OL) needs G - 1 independent conditions on
   the one-leaf polynomial sigma^G Xi_G(D+1), which lies in span(sigma, ..., sigma^(G-1)) (proof_status s.1.7(a)). Theorem 3.1
   gives one (sigma = 2, given P2); Route A's divisibility argument gives Z(D) more (A_mv\report.md s.3.5); Lemma D at further values
   of tau, or any TR argument valid for all tau, would give the rest.
