# Route C report (round 1, 2026-10-31)

(Saved by the main session from the workflow result wf_40196105-891; subagents cannot write report files.)

## claims
[
   {
      "depends_on" : "Lemma 2 of proof_M_residue.md (c(R^1 pi_* omega_log^{-1}) = Lambda^vee(1) K, Mumford GRR; Gauss-Bonnet checked h<=4); pure dimension 2g-3+n of H_g(mu) for strictly meromorphic mu (Farkas-Pandharipande, J. Inst. Math. Jussieu 17 (2018); Boissy, Comment. Math. Helv. 90 (2015)); Fulton's intersection theory (projective bundles, excision), Vistoli for DM stacks",
      "evidence" : "Complete proof in the summary: principal-parts map phi: P -> F = R^1 pi_* omega_log^{-1}; projective-bundle formula on P(P) localizes the integral to a class beta in A_{J-1}(P(ker phi)); dimension bound J_Y - e - 1 for (support, section up to scale); Lemma V. Computer-verified at all lattice points h = 2,3,4, r = 1,2.",
      "statement" : "Theorem H: for all h >= 0, r >= 1 (2h-2+r>0), integers j_i >= 0 with J >= 1 and D = Sum(j_i - 1) <= 2h-3, the Hurwitz component Xi^(0)_h(R) = Int_{Mbar_{h,r}} c(T(-log)) prod 1/(1 - j_i psi_i) vanishes. In particular sub-lemma A2 (Hurwitz component of Lemma A) holds for every h >= 2 and every irreducible R, any r.",
      "status" : "proved-modulo-literature"
   },
   {
      "depends_on" : "Theorem H proof",
      "evidence" : "Same proof (deg alpha >= J_Y - 1 > dim T'); partial fractions identify the leaf spans. Computer-verified: 286 (j,b) pairs, h = 2..4, r = 1,2, all zero (check_Hprime_out.txt).",
      "statement" : "Strong form H': Int c_m(-R pi_* S_D) prod psi_i^{b_i} = 0 for all m >= N-J+1 and 0 <= b_i <= j_i (Sum b_i = N - m); equivalently Xi^(0) vanishes for every choice of leaves in span(lam_0..lam_{j_i}) at x_i, not all lam_0.",
      "status" : "proved-modulo-literature"
   },
   {
      "evidence" : "check_H_out.txt, h = 2,3,4, r = 1,2, all boundary points D = 2h-2",
      "statement" : "The range is sharp for this mechanism: at D = 2h-2 the support can be the whole curve (holomorphic strata, one extra dimension) and Xi^(0) is non-zero at the tested boundary points (e.g. Xi^(0)_2(3) = -1/320, Xi^(0)_3(5) = -65/24192, Xi^(0)_4(7) = -1967/276480).",
      "status" : "computer-verified"
   },
   {
      "depends_on" : "lambda_g formula <lambda_g prod tau> = b_g <tau_0^{2g} prod tau>_0 (Faber-Pandharipande / Getzler-Pandharipande)",
      "evidence" : "Genus-0 psi integrals + dilaton sum as in proof_status 1.3; D(R u 0^{2h}) <= -3 iff D <= 2h-3",
      "statement" : "Genus-0 form: Xi^(0)_0(R'') = (J-1)(J-2)...(J-n+3) on Mbar_{0,n}; with the lambda_g formula, Xi^(h)_h(R) = b_h Xi^(0)_0(R u 0^{2h}), so the lambda_g component vanishing in the range is Theorem H at genus 0 with 2h zero legs.",
      "status" : "proved-modulo-literature"
   },
   {
      "evidence" : "c_1(T_{2,1}(-log)) = pi^*c_1(T_2) - psi, Int_{Mbar_2} lambda_1 lambda_2 = 1/5760, Xi^(2)_2(2) = 0 (lambda_g); cross-checked with comp_g2_r1.txt",
      "statement" : "Individual Hodge insertions do not vanish: Int_{Mbar_{2,1}} c(T(-log)) lambda_1 lambda_2 lam_2 = 1/2880 and Int c(T(-log)) lambda_2 lam_2 = -1/2880; only e_2(alpha+alpha^2) = lambda_2 + lambda_1 lambda_2 pairs to zero. Any proof of the components 1 <= k <= h-1 must use the CY relation of the weights.",
      "status" : "proved"
   },
   {
      "evidence" : "Degree count in (2.3) with the extra Hodge factor; string equation Xi(1,R) = J_R Xi(R).",
      "statement" : "Proposition L: since e_k(alpha^2) = 0 (Mumford), deg e_k(alpha+alpha^2) <= 2k-1, and for r >= 2k every component is fully localized: Xi^(k)_h(R) = Int_beta [e_k(alpha+alpha^2) prod p_{j_i}/(1-zeta)]_{J-1}. Any R reaches r >= 2k by adding legs j = 1 (string equation, component-wise).",
      "status" : "proved"
   },
   {
      "evidence" : "lambda_1|_Gamma = Sum_v lambda_1(v) = c_1(det E_v); Lemma V kills Y-terms and all Z-terms with r_{Z,+} + e >= 2.",
      "statement" : "Proposition L1 (lambda_1-component, every R): only components of P(ker phi) whose section is supported on a subcurve containing all markings, attached to the rest by exactly one separating node, with saturated Y-part (dimension J-2, exact double zero at the node) can contribute to Xi^(1)_h(R).",
      "status" : "proved"
   },
   {
      "depends_on" : "local model (T1); Fulton excess intersection formula",
      "evidence" : "Complete except for writing out (T1) (first-order plumbing: sigma = (sigma_Y^{(-1)} + O(t), -t c delta_{q'} + O(t^2))). Computer-verified for h <= 4, r <= 2 (and r <= 5 at genus 4 via the existing multiprobe data).",
      "statement" : "The lambda_1-component Xi^(1)_h(R) vanishes for all h and all R in the range. Argument: by L + L1 and the local plumbing model (T1) near a separating node, the zero scheme is Z_Y x {t=0} x Mbar_{g_Z,1} with excess bundle pi^*T_{Mbar_{g_Z}}(-log) (x) O(1); the contribution is a multiple of Int pi^*(c_top(T_{g_Z}) lambda_1) = 0.",
      "status" : "conjectural"
   },
   {
      "evidence" : "Local computation at attaching nodes plus degree/range bookkeeping (D_Y >= 2g(Y) - 2 + Sum m_e).",
      "statement" : "For k >= 2 a product formula on support types would give contributions (Y-factor) x Xi^(k_Z)_{g_Z}(R_Z u {m_e - 1}), with the Z-factor always in the range (D_Z <= 2g_Z - 3 for at least one component of Z). Lemma A would then follow by induction, except for self-referential rational tails (s = z^{1-j} d/dz, m = j+1) when j <= 2k-1, whose coefficient is not computed.",
      "status" : "conjectural"
   },
   {
      "evidence" : "fit_decomp_out.txt: NOT in span at genus 2 (k=1), 3 (k=1,2), 4 (k=1,2,3), r = 2",
      "statement" : "The decomposition Xi^(k)_h(R) = Sum_{g=0}^{h-k} c_{h,k,g} Xi^(0)_g(R u 0^{2(h-g)}) with R-independent constants fails for r = 2 (it is automatic for r = 1).",
      "status" : "refuted"
   },
   {
      "evidence" : "gauss_bonnet_out.txt (9/9 OK)",
      "statement" : "Gauss-Bonnet check of the log-tangent identification: Xi^(0)_h(0^r) = chi(M_{h,r}) for h = 1..4, r = 1, 2.",
      "status" : "computer-verified"
   }
]


## key_formulas
[
   "c(T(-log)) prod_i 1/(1 - j_i psi_i) = c(F - P) prod_i p_{j_i}(psi_i),  F = R^1 pi_* omega_log^{-1} (rank N = 3h-3+r),  P = pi_*(S_D|_{Sum j_i x_i}) (rank J),  c(P) = prod_i prod_{k=1}^{j_i} (1 - k psi_i)",
   "ker(phi_C: P_C -> F_C) = H^0(C, S_D),  S_D = omega_log^{-1}(Sum j_i x_i) = omega^{-1}(Sum (j_i - 1) x_i)",
   "on p: P(P) -> X:  p_*( c_N(p^*F (x) O(1)) zeta^j ) = c_{N-J+1+j}(F - P),  zero locus of the section = P(ker phi)",
   "Xi^(0)_h(R) = Int_beta [ prod_i p_{j_i}(psi_i) / (1 - zeta) ]_{J-1},  beta in A_{J-1}(P(ker phi))",
   "dim{(support Y, [s])} <= Sum_{v in Y}(2g_v - 3 + n_v) + |V(Y)| - 1 <= J_Y - e - 1  (dim H_g(mu) = 2g-3+n for strictly meromorphic mu)",
   "Y = C impossible iff Sum_v d_v = D - (2h-2) < 0, i.e. D <= 2h-3",
   "Xi^(0)_0(R) = (J-1)(J-2)...(J-n+3) on Mbar_{0,n};  Xi^(h)_h(R) = b_h Xi^(0)_0(R u 0^{2h})",
   "Mumford: c_u(E) c_u(E^vee) = 1, e_k(alpha^2) = 0 (k>=1), so deg e_k(alpha + alpha^2) <= 2k-1; full localization iff r >= 2k",
   "Int_{Mbar_{2,1}} c(T(-log)) lambda_2 lam_2 = -1/2880,  Int c(T(-log)) lambda_1 lambda_2 lam_2 = +1/2880"
]


## next_steps
[
   "Write out (T1) in plumbing coordinates (zw = t at the separating node; first-order obstruction -t c delta_{q'}; identification of F near the boundary with F_Y^{(-1)} (+) F_Z). Together with Props L and L1 this proves the lambda_1-component Xi^(1)_h(R) = 0 for all h and R. Lemma A at genus 3 would then need only k = 2.",
   "Formulate the moduli of non-zero sections of S_D up to scale (= P(ker phi), perfect obstruction theory (R pi_* S_D)^vee, virtual class = localized top Chern class) and prove a splitting formula along the boundary strata (support types), analogous to the splitting axiom of GW virtual classes.",
   "Compute the coefficient of the self-referential rational-tail contributions (Y = P^1 with x_1 at 0, attaching point at infinity, s = z^{1-j} d/dz, p_Y legs j = 1) and compare it with the string factor (j)_p. If they differ, induction on (h, r, j) closes Lemma A for all k.",
   "Alternatively, a cosection-localization (Kiem-Li) treatment of the CY triple S_D (+) O (+) O with weights (1, tau, -1-tau), using the CY relation that statement C3-a shows to be essential.",
   "Have Theorem H checked by an independent referee (Steps 2-6; the use of Farkas-Pandharipande / period coordinates; Lemma V on DM stacks via coarse spaces). Then record A2 as proved in proof_status_mixed_scale.md s.2."
]


## obstructions
[
   "Components 1 <= k <= h-1: after full localization (Prop. L) the integrand contains Hodge classes of the ZERO part Z of the section; the dimension argument cannot control them. One needs the localized (excess) class of the jump locus on the support types with deg(Hodge_Z) >= r_{Z,+} + e + (deficiency of the Y-part).",
   "k = 1: the only remaining input is the local plumbing model (T1) at a separating node (first-order term -t c delta_{q'} of the obstruction, regular embedding of the zero scheme, excess bundle pi^*T_{Mbar_{g_Z}}(-log) (x) O(1)). It is sketched, not written out.",
   "k >= 2: non-saturated support types with higher-order zeros m_e >= 3 at the attaching points, and Z-parts with their own jumps, need a general product (splitting) formula for the localized class (the virtual class of the moduli of sections of S_D up to scale). With it, rational tails with s = z^{1-j} d/dz reproduce the unknown Xi^(k)_h(j, 1^p) when j <= 2k-1; their coefficient must be shown to differ from the string factor (j)_p, otherwise the recursion is vacuous.",
   "Route C2 cannot remove the jump loci: roots of S_D^rho coincide with S_D on rational tails (no torsion; zero twist at separating nodes when s, a_i are divisible by rho); other choices change the Chiodo kappa/psi coefficients, and JPPZ polynomiality only holds for rho >> 0."
]


## report
G:\My Drive\Projects\Old\refGW\hodge_vertex_conifold\lemmaA\C_intersection\ (report.md was NOT written: the harness refuses report files from subagents ('return findings as text'); the full report, including the complete proof of Theorem H, is in the summary field below. Scripts and outputs are in this folder: polylib.pm, components.pl, comp_g{2,3,4}_r{1,2}.txt, comp_g4_r2.log, dilsum.pl, comp_g{1,2}_r*_psi.txt, check_H.pl + check_H_out.txt, check_Hprime.pl + check_Hprime_out.txt, gauss_bonnet.pl + gauss_bonnet_out.txt, fit_decomp.pl + fit_decomp_out.txt)

## requested_checks
[
   "Referee check of the proof of Theorem H (summary, steps (0)-(6)), in particular the dimension bound J_Y - e - 1 (per component J_v - n_v^bd - 1 from dim H_g(mu) = 2g-3+n) and the factorisation of the integrand on the lifted cycle W' in P(P_Y) x M_Z.",
   "Optional Mathematica check (main session, once hvc_lib's hodgeClass4 has a ch_5 branch): genus 5 Hurwitz component Xi^(0)_5(j) = Int Lambda^vee(1) K lam_j for j = 1..9, and a few r = 2 points (2,2), (2,3), (3,4). Theorem H predicts zeros for D <= 7 (j <= 8; (2,2), (2,3), (3,4)) and a non-zero value at j = 9. Only the single class Lambda^vee(1) is needed, which is cheaper than the full CY class.",
   "Optional: re-run C_intersection/check_H.pl, check_Hprime.pl and gauss_bonnet.pl after any regeneration of the lemmaA/data tables (perl components.pl h r first)."
]


## route
C (intersection theory on Mbar_{h,r}: jump loci, excess intersection; r-th roots; intermediate Hodge components)

## summary
MAIN RESULT (new): the Hurwitz component (k = 0) of Lemma A vanishes for every h and every R in the range, any r. This settles sub-lemma A2 of proof_status_mixed_scale.md completely (previously only j = 2 and j = 3 were known). No TR of any kind is used anywhere.

THEOREM H. Let h >= 0, r >= 1, 2h-2+r > 0, j_i >= 0 integers with J >= 1 and D = Sum(j_i - 1) <= 2h-3. Then
Xi^(0)_h(R) = Int_{Mbar_{h,r}} c(T(-log)) prod_i 1/(1 - j_i psi_i) = 0.
Strong form (H'): Int c_m(-R pi_* S_D) prod_i psi_i^{b_i} = 0 for every m >= N-J+1 and all 0 <= b_i <= j_i with m + Sum b_i = N (N = 3h-3+r). By partial fractions this is the same as the vanishing for all leaves in span(lam_0..lam_{j_i}) at x_i (not all lam_0).

PROOF.
(0) Setting. X = Mbar_{h,r} (Q-coefficients, Vistoli), S_D = omega_log^{-1}(Sum j_i x_i) = omega^{-1}(Sum(j_i - 1)x_i). F := R^1 pi_* omega_log^{-1} is locally free of rank N (omega_log^{-1} has negative degree on every component), and c(F) = Lambda^vee(1) K (Lemma 2 of proof_M_residue.md, GRR). P := pi_*(S_D restricted to Sum j_i x_i) is locally free of rank J with graded pieces L_i^{-k}, k = 1..j_i, so c(P) = prod_i prod_{k=1}^{j_i}(1 - k psi_i). The sequence 0 -> omega_log^{-1} -> S_D -> S_D|_Q -> 0 gives a bundle map phi: P -> F with ker(phi_C) = H^0(C, S_D) at every point (base change), and -R pi_* S_D = F - P in K-theory, i.e. identity (S): c(T(-log)) prod 1/(1 - j_i psi_i) = c(F - P) prod p_{j_i}(psi_i).
(1) Projective bundle. On p: P(P) -> X (lines), zeta = c_1(O(1)), O(-1) -> p^*P -> p^*F is a section sigma of p^*F (x) O(1) (rank N). From p_*(zeta^{J-1+k}) = [1/c(P)]_k: p_*(c_N(p^*F (x) O(1)) zeta^j) = c_{N-J+1+j}(F - P). sigma vanishes exactly on Zt = P(ker phi) = {(C,[s]) : 0 != s in H^0(C, S_D)}; off Zt the top Chern class restricts to 0, so by excision c_N(...) cap [P(P)] = i_* beta with beta in A_{J-1}(Zt) (Fulton Prop. 1.8 / localized top Chern class, Prop. 14.1). Since deg prod p_{j_i} <= J - r_+ <= J - 1, only the localized degrees occur and
Xi^(0)_h(R) = Int_beta [prod_i p_{j_i}(psi_i)/(1 - zeta)]_{J-1}.
It suffices that every monomial zeta^a prod psi_i^{b_i} (b_i <= max(j_i - 1, 0)) integrates to 0 on every irreducible (J-1)-dimensional W in Zt.
(2) Sections on stable curves. For a component C_v with normalisation Ct_v: nu^* S_D = omega_{Ct_v}^{-1}(Sum_{i in v}(j_i - 1)x_i - Sum branches). Let Y be the support of s (components where s != 0), Z the rest, e = number of Y-Z nodes; s vanishes to order >= 2 at the Y-branch of each Y-Z node. s_v is a non-zero section of a bundle of degree d_v = 2 - 2g_v + D_v - n_v' - n_v^bd >= 0. If Y = C then Sum d_v = D - (2h-2) < 0: impossible in the range. Hence Y != C, e >= 1, and every component of Y has a branch (where s vanishes).
(3) Dimension bound. eta_v = 1/s_v is a meromorphic differential with at least one pole; (Ct_v, markings, branches, free poles (at most d_v)) lies in H_{g_v}(mu_v) = {omega = O(Sum mu_k p_k)}, of pure dimension 2g - 3 + n when mu has a negative entry (Farkas-Pandharipande, J. Inst. Math. Jussieu 17 (2018): twisted canonical divisors of strictly meromorphic type have pure codimension g in Mbar_{g,n}; equivalently period coordinates, Boissy, Comment. Math. Helv. 90 (2015) s.2.3; BCGGM Duke 167 (2018)). This gives dim <= J_v - n_v^bd - 1 per component, and with the relative scalars of the eta_v modulo one overall scalar: the data (Y, [s]) vary in dimension <= J_Y - e - 1 (J_Y = Sum of j_i over markings on Y).
(4) Factorisation. Lift W to W' over the gluing map xi_Gamma of the generic dual graph; xi^*P = P_Y (+) P_Z, and W' lies in P(P_Y) x M_Z (principal parts at markings on Z vanish). On W': zeta^a prod psi^b = f^*alpha . g^*beta', alpha = zeta^a prod_{i on Y} psi_i^{b_i} (first Chern classes of line bundles on P(P_Y)), beta' = prod_{i on Z} psi_i^{b_i}; dim f(W') <= J_Y - e - 1 by (3).
(5) Lemma V: on a projective T of dimension t, any product of p > t first Chern classes of line bundles, pulled back to any V -> T, annihilates A_*(V) (induction on t, cutting by very ample divisors; DM stacks via coarse spaces).
(6) deg alpha = J - 1 - Sum_{Z} b_i >= J - 1 - Sum_{Z, j_i >= 1}(j_i - 1) = J_Y + r_{Z,+} - 1 > J_Y - e - 1 >= dim f(W'), since e >= 1. Lemma V gives Int_{W'} = 0. QED. (H': with b_i <= j_i, deg alpha >= J_Y - 1, same conclusion.)
Remarks: the range enters only through Y != C. At D = 2h-2 the support can be all of C (holomorphic strata have one more dimension) and the bound fails by exactly one: Weierstrass divisor, Xi^(0)_2(3) = -1/320 != 0. Props 3, 4 of proof_M_residue.md are the cases j = 2, 3; the new point is to work on P(P), where the localized class has the single dimension J - 1 and zeta remembers the section, so the generalisation of bound (4.2) suffices for every R.

COROLLARY H0 (lambda_g component re-derived). Theorem H at genus 0 is the closed formula Xi^(0)_0(R'') = (J-1)(J-2)...(J-n+3) (zero exactly for 1 <= J <= n-3, i.e. D'' <= -3). By the lambda_g formula Xi^(h)_h(R) = b_h Xi^(0)_0(R u 0^{2h}), and D(R u 0^{2h}) <= -3 iff D <= 2h-3: the k = 0 and k = h components are instances of the same vanishing (genus h, resp. genus 0 with 2h zero legs).

ROUTE C2 (r-th roots). Not needed for k = 0. Negative observation: for roots of S_D^rho (s = -rho, a_i = -rho j_i) all twists at separating nodes are 0 and rational curves carry no torsion, so the root equals S_D on every rational tail and the jump loci there persist for all rho; any other (rho, s, a_i) changes the kappa/psi coefficients B_{m+1}(s/rho), B_{m+1}(a_i/rho) of Chiodo's formula, and JPPZ polynomiality holds only for rho >> 0, so rho = 1 cannot be reached by evaluation. The only GRR/Chiodo-type input used is Mumford's GRR through Lemma 2, checked by Gauss-Bonnet (below).

ROUTE C3 (components 1 <= k <= h-1): partial results and the precise obstruction.
(a) The CY structure is essential: Int_{Mbar_{2,1}} c(T(-log)) lambda_1 lambda_2 lam_2 = 2 Int_{Mbar_2} lambda_1 lambda_2 = 1/2880 and hence (Xi^(2)_2(2) = 0, e_2 = lambda_2 + lambda_1 lambda_2) Int c(T(-log)) lambda_2 lam_2 = -1/2880. Individual Hodge insertions do not vanish.
(b) Mumford gives e_k(alpha^2) = 0, so deg e_k(alpha+alpha^2) <= 2k-1. Proposition L (proved): if r >= 2k, Xi^(k)_h(R) = Int_beta [e_k(alpha+alpha^2) prod p_{j_i}/(1-zeta)]_{J-1} (fully localized). Any R can be brought there by adding j = 1 legs, using the string equation Xi(1, R) = J_R Xi(R) (valid component-wise).
(c) Proposition L1 (proved): for k = 1 and every R, only components of Zt whose section is supported on a subcurve Y carrying ALL markings and attached to the rest by ONE (separating) node can contribute, with Y-part saturated (dim = J - 2, which forces an exact double zero at the node). Proof: lambda_1 restricts to Sum_v lambda_1(v); Y-terms die by Lemma V; Z-terms need r_{Z,+} + e <= 1.
(d) Local model (T1), sketched but not written out: with plumbing zw = t at the separating node, sigma = (sigma_Y^{(-1)} + O(t), -t c delta_{q'} + O(t^2)), where c != 0 is the second-order coefficient of s at q and delta_{q'} spans L_{q'}^vee in T_{Mbar_{g_Z,1}}(-log). Then the zero scheme is Z_Y x {t=0} x Mbar_{g_Z,1} with excess bundle pi^*T_{Mbar_{g_Z}}(-log) (x) O(1), whose contribution is a multiple of Int pi^*(c_{3g_Z-3}(T_{g_Z}) lambda_1) = 0. So the k = 1 component vanishes for all h, R MODULO (T1) (conjectural; consistent with data for h <= 4).
(e) For k >= 2 the same scheme suggests a product formula in which each support type contributes (Y-factor) x Xi^(k_Z)_{g_Z}(R_Z u {m_e - 1}). The Z-factor is always in the range (D_Z <= 2g_Z - 3 for at least one component of Z; computed). So Lemma A would follow by induction, EXCEPT that rational tails with s = z^{1-j} d/dz (m = j + 1) reproduce the unknown Xi^(k)_h(j, 1^p) itself when j <= 2k-1. Their coefficient is not computed. This is the precise point where I stop.
(f) The naive ansatz Xi^(k)_h(R) = Sum_g c_g Xi^(0)_g(R u 0^{2(h-g)}) with R-independent constants is refuted at r = 2 (genus 2, 3, 4).

EXACT CHECKS (Perl + Math::BigRat, tables of lemmaA\data; component extraction by solving the polynomial in -1/sigma at h+1 values of x, consistency verified at an extra x):
- Gauss-Bonnet check of Lemma 2: Xi^(0)_h(0^r) = chi(M_{h,r}) for h = 1..4, r = 1, 2 (-1/12, 1/12, 1/120, -1/40, -1/252, 5/252, 1/240, -7/240): 9/9 OK.
- Theorem H and all components at every lattice point j_i >= 1, D <= 2h-3, r = 1, 2, h = 2, 3, 4: all zero. At D = 2h-2 the Hurwitz component is non-zero, e.g. h=3 j=5: -65/24192; h=4 j=7: -1967/276480; h=4 (4,4): -44/1575.
- H' (b_i <= j_i, Sum b <= J-1, j_i >= 0 allowed): 286 (j, b) pairs over h = 2..4, r = 1, 2, all zero.
- The psi-table route (dilsum.pl) reproduces the q-table components at genus 2 (independent check of the data pipeline).
- r >= 3 points (genus 3: {2,2,2}; genus 4: the seven missing points) are covered by the existing hvc_multiprobe outputs (full Xi = 0 identically in x, hence every component).

