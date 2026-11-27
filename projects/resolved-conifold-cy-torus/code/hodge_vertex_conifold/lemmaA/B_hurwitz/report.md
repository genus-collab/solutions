# Route B report (round 1, 2026-10-31)

(Saved by the main session from the workflow result wf_40196105-891; subagents cannot write report files.)

## claims
[
   {
      "depends_on" : "ELSV (Invent. Math. 146, 2001); Okounkov, Toda equations for Hurwitz numbers (MRL 7, 2000), eq. (10), specialised at P' = (1,0,..), q = 1 and checked on 50 coefficients; Hurwitz genus-0 formula (Goulden-Jackson, PAMS 1997); string/dilaton.",
      "evidence" : "Proof in the summary, sections 2-5: Lambert structure, normalised Toda equation, T -> oo leading-order identity, pure-unit cancellation, and a counting argument showing that a non-zero term needs D >= 2h. Independent exact checks: (REC) against Hurwitz numbers at genus 1-3 (80 values) and at genus 4-5 (20 points, every r <= 5 irreducible genus-4 point included); (REC) against the Hodge data at genus 2-4 (110 values).",
      "statement" : "Theorem A: for every h >= 1 and every multiset R of integers >= 2, Xi^(0)_h(R) = Int c(T(-log)) prod 1/(1-j_i psi_i) vanishes whenever D(R) <= 2h-3 (this is sub-lemma A2 of Lemma A, for all h) or D(R) = 2h-1 (a new extra zero, for every r).",
      "status" : "proved-modulo-literature"
   },
   {
      "depends_on" : "Theorem C (same literature inputs).",
      "evidence" : "For r = 1, (REC) becomes 2h Phi_h = (j-2h) sum_{g'<h} (2/(2(h-g'+1))!) Phi_{g'}. Its generating function solves t((1+C)F)' = jCF with 1+C = S^{-2}, giving F = S^2 H^j. Checked against (REC) up to genus 8 (109 values), against the stored Hodge polynomials for h <= 4, and against a Hurwitz-number value at genus 5 (Xi^(0)_5(9) = -4401/112640).",
      "statement" : "Theorem B: Xi^(0)_h(j) = (j-1)(j-2)...(j-2h+2) [t^{2h}] S(t)^2 H(t)^j = (j-1)...(j-2h+2) [y^j] x(y)^{j-2h}, with x = y H(x), for all h, as polynomials in j. This is the Hurwitz limit of (6.1), proof_M_residue s.6a(iii), including the extra zero j = 2h.",
      "status" : "proved-modulo-literature"
   },
   {
      "depends_on" : "ELSV, Okounkov's Toda eq. (10), Hurwitz genus-0 formula.",
      "evidence" : "Lemma 6.2 combined with exp(Lambda_empty) = 1. Matches independent Hurwitz-number values at genus 1-5 and the exact Hodge data at genus 2-4, with no mismatch.",
      "statement" : "Theorem C (Toda recursion): for delta = D-2h+2 != 0, (e+1) delta Xi^(0)_h(R) = sum over set partitions pi of R and (h_B,k_B) with sum(h_B+k_B-1) = h, excluding the single term (h,1,R), of prod_B (2/(2k_B)!) (delta_B)_{(2k_B)} Xi^(0)_{h_B}(R_B). Genus-0 input: Xi_0(R) = (J-1)...(J-r+3), with 1/J for r = 2 and 1/(j(j+1)) for r = 1. It determines all values with D != 2h-2 from genus 0.",
      "status" : "proved-modulo-literature"
   },
   {
      "depends_on" : "String and dilaton equations; Lemma 1 of proof_M_residue (kappa form).",
      "evidence" : "Proof uses string and dilaton (binomial expansion of psi^2 lam_1 = lam_1 - tau_0 - tau_1) and the substitution u = s/(1+s). The whole s-polynomial was checked against the q tables (429 comparisons).",
      "statement" : "Lemma 4.2 (Lambert structure): e^{-uJ} G_{R,h}(u e^{-u}) = sum_{k=0}^{3h-3+r} p_k (1-u)^{-(2h-2+r+k)}, where P(eps) = sum_k p_k eps^k = <Lambda^vee(1) K_{eps-1} lam_R> and p_0 = Xi^(0)_h(R). It holds for any class pulled back from Mbar_h, in particular the CY class. Hence Lemma A at genus h is equivalent to the decay e^{-uJ} G^Omega_{R,h}(u e^{-u}) = O(u^{-(2h-1+r)}) at u = oo for R in the range.",
      "status" : "proved"
   },
   {
      "depends_on" : "Lemma 4.1; ELSV for the Hurwitz form.",
      "evidence" : "Laguerre summation of Lemma 4.1. Agrees with the Hodge data at genus 2-4 (429 + 8 comparisons, including r = 3..5 at genus 4).",
      "statement" : "Proposition D (finite formula): Xi^(0)_h(R) = sum_{c=0}^{n0} (-1)^c (A_h(R,1^c)/c!) sum_{i=0}^{n0-c} C(n0+e, n0-c-i) (J+c)^i/i!, with n0 = 3h-3+r, e = 2h-2+r, and A_h(R,1^c) the ELSV-normalised connected Hurwitz numbers H_{h,R+1^c}.",
      "status" : "proved"
   },
   {
      "depends_on" : "Okounkov's Toda equation; ELSV.",
      "evidence" : "u_1 = 0 explicitly; for m >= 2 it follows from Lemma 6.2 at (m, empty). The Bernoulli identity (x f' - f)(cosh x - 1) = -x^2/2 identifies the solution.",
      "statement" : "Pure-unit cancellation: u_m = sum_{g'+k=m+1} (2/(2k)!) L_{g',k,empty} = 0 for all m >= 1. As a byproduct this gives Harer-Zagier, chi(M_g) = B_{2g}/(2g(2g-2)) (Dubrovin-Yang-Zagier's argument).",
      "status" : "proved-modulo-literature"
   },
   {
      "depends_on" : "Okounkov, MRL 7 (2000) 447-453, eq. (10).",
      "evidence" : "Homogeneity tau(P,(p'_1,0,..),beta,q) = Z((q p'_1) o P). Checked exactly on all 50 non-zero coefficients (degree <= 6, beta^10) computed from the character formula (t_toda.out).",
      "statement" : "Okounkov's eq. (10), specialised at P' = (1,0,0,..) and q = 1, gives d/dp_1(E H) = exp(H(e^beta o p) - 2H(p) + H(e^{-beta} o p)) for the connected single-Hurwitz series.",
      "status" : "proved-modulo-literature"
   },
   {
      "depends_on" : "HB.pm, HBconn.pm, Data.pm.",
      "evidence" : "81 comparisons on the first tables (wk symbols evaluated with my own DVV code) and 236 on the regenerated tables (genus 1-4, l <= 3; genus-4 full unisolvent sets), 0 mismatches.",
      "statement" : "ELSV normalisation: Int Lambda^vee(1) prod 1/(1-m_i psi_i) = H_{g,m} |Aut m| prod(m_i!/m_i^{m_i}) / b!, with H_{g,m} = (1/|m|!) #{transitive (t_1..t_b, sigma)}. This matches the exact Hodge tables (sigma -> oo component).",
      "status" : "computer-verified"
   },
   {
      "depends_on" : "Theorem A.",
      "evidence" : "Part of Theorem A; every listed case also confirmed exactly by Hurwitz numbers and/or the Hodge data.",
      "statement" : "Extra zero for every r: Xi^(0)_h(R) = 0 at D = 2h-1, e.g. (4), (2,3), (2,2,2) at h = 2; (6), (2,5), (3,4), (2,2,4), (2,3,3), (2,2,2,3) at h = 3; (8), (2,7), (3,6), (4,5), (2,2,2,5), (2,2,2,2,4) at h = 4; (10) at h = 5. Values at D = 2h-2 are non-zero.",
      "status" : "proved-modulo-literature"
   },
   {
      "depends_on" : "Okounkov/Ueno-Takasaki.",
      "evidence" : "t_mkp.out: 52 coefficients (degree <= 4, beta^8) vanish; with the opposite sign 29 fail. The transfer to Marino-Vafa is via the Ueno-Takasaki bilinear identity at equal second times; it is standard but not written out here.",
      "statement" : "The modified-KP relation 2 d_{p_2}(F - F_1) + d_{p_1}^2 (F + F_1) + (d_{p_1}(F - F_1))^2 = 0, with F_1(p) = F(e^beta o p), holds for the single-Hurwitz series (and, by 2D Toda theory, for any P', in particular for the Marino-Vafa framed vertex).",
      "status" : "computer-verified"
   }
]


## key_formulas
[
   "Xi^(0)_h(R) = Int_{Mbar_{h,r}} Lambda^vee(1) K prod 1/(1 - j_i psi_i) = 0 for D(R) <= 2h-3 or D(R) = 2h-1 (parts >= 2, all h, all r)",
   "Lambert structure: e^{-uJ} G_{R,h}(u e^{-u}) = sum_{k=0}^{3h-3+r} p_k (1-u)^{-(2h-2+r+k)}, G_{R,h}(y) = sum_c y^c/c! Int Lambda^vee(1) lam_R lam_1^c, p_0 = Xi^(0)_h(R)",
   "Finite formula: Xi^(0)_h(R) = sum_{c=0}^{n0} (-1)^c (A_h(R,1^c)/c!) sum_{i=0}^{n0-c} C(n0+e, n0-c-i)(J+c)^i/i!, n0 = 3h-3+r, e = 2h-2+r, A_h(R,1^c) = H_{h,R+1^c} |Aut| prod(j!/j^j)/b!",
   "Toda (Okounkov eq. (10) at P' = (1,0,..)): d/dp_1 (E H) = exp(H(e^beta o p) - 2H(p) + H(e^{-beta} o p))",
   "Normalised Toda: theta D_J g_{g,R} = T S'_{g,R}, D_J = (theta+J)/(1-T), theta = T d/dT, z = T e^{-T}, g_{g,R} = e^{-JT} [q_R]H_g",
   "Leading order: (e+1)(J-e) a_{g,R} = [beta^{2g} q_R] exp(Lambda) - L_{g,1,R}, L_{g',k,R'} = (J'-e')_{(2k)} a_{g',R'}",
   "(REC): (e+1) delta Xi_h(R) = sum' prod_B (2/(2k_B)!) (delta_B)_{(2k_B)} Xi_{h_B}(R_B), delta = D-2h+2, sum_B (h_B+k_B-1) = h; Xi_0(R) = (J-1)...(J-r+3)",
   "Counting: sum_i delta_i - sum_i 2k_i = D - 2h, so every term vanishes for D <= 2h-1",
   "Theorem B: Xi^(0)_h(j) = (j-1)...(j-2h+2) [t^{2h}] S^2 H^j = (j-1)...(j-2h+2) [y^j] x^{j-2h}, x = y H(x), S = (t/2)/sinh(t/2), H = S e^{-t S'/S}",
   "Pure-unit cancellation: sum_{g+k=m+1} (2/(2k)!) (2g-1) B_{2g} (2g+2k-3)!/(2g)! = 0 for m >= 1 (equivalently (x f' - f)(cosh x - 1) = -x^2/2, f = (x/2)coth(x/2))"
]


## next_steps
[
   "Main session: save the 'summary' field as lemmaA\\B_hurwitz\\report.md (the subagent harness blocked writing report files).",
   "Use Remark 4.4: Lemma A in full is equivalent to the decay at u = oo of the CY unit-part generating function. Look for a charge-direction equation for the framed one-leg vertex at fixed P' = q^rho; integer tau suffices by proof_status s.1.7(a). The q-shift property s_nu(q^rho) e^{n beta |nu|} = s_nu(q^{rho + n tau}) combined with Okounkov's bilinear identity (9) is a candidate.",
   "Alternative for the CY components: carry out the T -> oo leading-order analysis of the Marino-Vafa cut-and-join flow in tau (all orders give algebraic constraints on the leading coefficients), with the k = 0 values proved here as boundary data at tau = oo.",
   "Referee pass on Lemmas 5.1, 6.1-6.3 and the counting in section 5 (the algebra is short; the numerical checks already exercise every step).",
   "Rec.pm computes the Hurwitz component for any (h,R) with D != 2h-2 in seconds. It can supply the main session with sigma -> oo boundary values for any further work on the remaining components."
]


## obstructions
[
   "The components 1 <= k <= h-1 of the CY class (full Lemma A) are not reached. The proof needs a Toda-type equation (first order in the unit-part derivative, with the degree operator, exponential of second differences in the charge direction) for the generating function of one-partition triple Hodge integrals. Marino-Vafa identifies it with the 2D Toda tau function tau(P,P',beta,1) at the principal specialisation P' = q^rho. There the d/dp'_1 of Okounkov's eq. (10) is no longer the degree operator: it is the vertex-operator expression <0|Gamma_+ E_{-1}(beta) e^{beta F_2} Gamma_-|0>, involving all p_k. This destroys the triangular induction on (g,r).",
   "The P-only modified-KP relation holds for the CY vertex, but it adds 2-parts, so it is an upward recursion. Its T -> oo analysis requires explicit cancellations, not the counting argument.",
   "The r = 1 CY formula (6.1) is not generated by the same type of ODE as the Hurwitz case, checked at order lambda^2. So a naive deformation of (REC) to the CY class fails.",
   "Literature dependence: Theorems A-C rest on ELSV and on Okounkov's Toda equation (MRL 2000, eq. (10)). Both are cited with exact statements, and their normalisations are checked numerically (236+81 and 50 exact comparisons). They are not re-proved here."
]


## report
NOT WRITTEN AS A FILE. The subagent harness refused report.md ('Subagents should return findings as text, not write report files'), so the complete report is the 'summary' field below. Scripts, libraries and all outputs are in G:\My Drive\Projects\Old\refGW\hodge_vertex_conifold\lemmaA\B_hurwitz\ (HB.pm, HBconn.pm, Data.pm, Rec.pm, *.pl, *.out, runs\).

## requested_checks
[
   "Save the summary text of this result as HVC\\lemmaA\\B_hurwitz\\report.md; the subagent harness refused to write report files.",
   "Optional independent Hodge-side confirmation at genus 5, via a cheap Mathematica run at h+1 rational x values: the sigma -> oo component (t^0 coefficient in t = c^2/x of c^{3h-3+2r} Xi) of Xi_5(R) should be 0 for R = (10), (2,2,2), (3,3,3), (2,2,2,2); 1327385/202752 for (11); 80 for (12); 102723/256 for (2,2,9); and -4401/112640 for (9). These are already confirmed by Hurwitz numbers.",
   "Pending result: runs\\b2pt_5_2_2_2_2_2_2_2.out, Xi^(0)_5(2^7) from Hurwitz numbers, a degree-33 computation still running; it should print 0. Do not kill it; it ends on its own.",
   "Referee re-derivation of Lemma 5.1 (normalised Toda equation), Lemma 6.1 (orders at T = oo) and the counting step of section 5."
]


## route
B: Hurwitz component of Lemma A, via ELSV, the character formula and Okounkov's Toda equation for Hurwitz numbers (no topological recursion of any kind)

## summary
ROUTE B REPORT (inline; to be saved by the main session as lemmaA\B_hurwitz\report.md).

I did not use topological recursion in any form: no spectral curve, no EO invariants, no loop equations, no BKMP/remodelling, no 'mod P2' result, and not as a cross-check either.

Tools used:
- ELSV;
- the Frobenius-Burnside character formula and the exponential formula;
- Hurwitz's genus-0 formula;
- Okounkov's Toda equation for Hurwitz numbers (a theorem about the character generating function, proved in the infinite wedge);
- string/dilaton;
- Lagrange inversion.

All checks are exact (Perl Math::BigInt/BigRat). I did not run wolframscript.

=== 0. CONVENTIONS AND RESULTS ===
- c = 1; R = (j_1..j_r), J = sum j_i, D = J - r, lam_m(psi) = 1/(1 - m psi).
- Xi^(0)_h(R) := Int_{Mbar_{h,r}} Lambda^vee(1) K prod lam_{j_i} = Int c(T(-log)) prod 1/(1 - j_i psi_i), with K = exp(-sum kappa_k/k). This is the sigma^0 (k = 0, Hurwitz) component of Xi_h(R).
- e := 2h-2+r and delta := J - e = D - 2h + 2. So D <= 2h-3 means delta <= -1, and D = 2h-1 means delta = 1.

THEOREM A (Hurwitz component of Lemma A, all genera, all r, plus an extra zero). For h >= 1 and every multiset R of integers >= 2:
  Xi^(0)_h(R) = 0 if D(R) <= 2h-3 or D(R) = 2h-1.
Status: proved modulo literature (ELSV; Okounkov MRL 2000 eq. (10); Hurwitz's genus-0 formula). Every other step is proved below.
Consequences:
- Sub-lemma A2 of proof_status s.2 holds for all h.
- With string/dilaton (s.1.2, valid componentwise), the Hurwitz component of the vanishing direction of Conjecture M holds for all h and r.
- The single class has an extra zero at D = 2h-1 for every r; for r = 1 this is j = 2h.

THEOREM B (r = 1: the Hurwitz limit of (6.1), proof_M_residue s.6a(iii), for all h). As polynomials in j,
  Xi^(0)_h(j) = (j-1)(j-2)...(j-2h+2) [t^{2h}] S(t)^2 H(t)^j = (j-1)...(j-2h+2) [y^j] x(y)^{j-2h}, with x = y H(x),
  S(t) = (t/2)/sinh(t/2), H = S exp(-t S'/S).
This includes the extra zero at j = 2h. Status: proved modulo literature, same inputs.

THEOREM C (Toda recursion). For h >= 1, parts >= 2 and delta != 0:
  (e+1) delta Xi^(0)_h(R) = sum' prod_{B in pi} (2/(2k_B)!) (delta_B)_{(2k_B)} Xi^(0)_{h_B}(R_B)        (REC)
- pi runs over set partitions of the r labelled parts.
- Each block gets h_B >= 0 and k_B >= 1 with sum_B (h_B + k_B - 1) = h.
- delta_B = J_B - (2h_B - 2 + r_B), and (x)_(n) = x(x-1)...(x-n+1).
- The single term pi = {R}, (h_B, k_B) = (h, 1) is omitted.
- Genus-0 input: Xi_0(R) = (J-1)(J-2)...(J-r+3) for r >= 3, 1/J for r = 2, 1/(j(j+1)) for r = 1.
(REC) determines every value with D != 2h-2 from genus 0 alone. Status: proved modulo literature.

PROP D (finite formula, B2). With n0 = 3h-3+r,
  Xi^(0)_h(R) = sum_{c=0}^{n0} (-1)^c (A_h(R,1^c)/c!) sum_{i=0}^{n0-c} C(n0+e, n0-c-i) (J+c)^i/i!,
where A_h(R,1^c) = Int Lambda^vee(1) lam_R lam_1^c = H_{h,R+1^c} |Aut(R+1^c)| prod_i (j_i!/j_i^{j_i}) / b! and b = e + J + 2c.
Status: proved (string/dilaton); the Hurwitz form uses ELSV.

LEMMA 4.2 (Lambert structure; it also holds verbatim for the CY class).
  e^{-uJ} G_{R,h}(u e^{-u}) = (1-u)^{-e} P(1/(1-u)) = sum_{k=0}^{n0} p_k (1-u)^{-(e+k)},  with p_0 = Xi^(0)_h(R),
where G_{R,h}(y) = sum_c A_h(R,1^c) y^c/c! and P(eps) = <Lambda^vee(1) K_{eps-1} lam_R>. Status: proved.

Byproduct: the case R = empty gives Harer-Zagier, chi(M_g) = B_{2g}/(2g(2g-2)). This is exactly Dubrovin-Yang-Zagier's argument. Our proof extends their T -> oo analysis from the empty profile to an arbitrary non-unit profile R.

=== 1. LITERATURE (exact statements, normalisations checked) ===

(L1) ELSV (Ekedahl-Lando-Shapiro-Vainshtein, Invent. Math. 146 (2001) 297-327, Thm 1.1; Graber-Vakil, Compositio 135 (2003) 25-36).
- Statement, for 2g-2+l > 0: H_{g,mu} = (b!/|Aut mu|) prod_i (mu_i^{mu_i}/mu_i!) Int_{Mbar_{g,l}} (1 - lambda_1 + ... + (-1)^g lambda_g)/prod_i(1 - mu_i psi_i), with b = 2g-2+l+|mu|.
- Normalisation of H_{g,mu}: H_{g,mu} = (1/|mu|!) #{(t_1..t_b, sigma): t_i transpositions, sigma in C_mu, product 1, transitive}, i.e. covers weighted by 1/|Aut|.
- Checked against the mv tables: 81 comparisons on the first tables (wk symbols evaluated by my own DVV routine), then 236 comparisons on the regenerated tables (genus 1-4, l <= 3; at genus 4 the full unisolvent simplex). 0 mismatches.

(L2) Frobenius-Burnside: #{g_i in C_i : g_1...g_k = 1} = (prod|C_i|/|G|) sum_chi prod chi(C_i)/chi(1)^{k-2}.
- Hence T(nu;n,b) = sum_{lambda |- n} dim(lambda) chi^lambda(nu u 1^{n-|nu|}) f_2(lambda)^b / (prod nu_i (n-|nu|)!), where T(nu;n,b) counts (t_1..t_b; labelled disjoint cycles of lengths nu) with product 1.
- Connected counts come from the exponential formula.
- Checked against brute force in S_3 and S_4 (17 cases).

(L3) Hurwitz genus 0 (Hurwitz 1891; Goulden-Jackson, Proc. AMS 125 (1997), Thm 1.1): H_{0,mu} = (b!/|Aut mu|) |mu|^{l-3} prod mu_i^{mu_i}/mu_i! for all mu.

(L4) Okounkov, 'Toda equations for Hurwitz numbers', Math. Res. Lett. 7 (2000) 447-453, arXiv:math/0004128, eq. (10). I read the PDF.
- Let tau(P,P',beta,q) = sum_lambda q^{|lambda|} e^{beta f_2(lambda)} s_lambda(P) s_lambda(P'), with f_2 = sum of contents. Then
  d^2 log tau/dp_1 dp'_1 = q tau(..,e^beta q) tau(..,e^{-beta} q)/tau^2.
- Specialise to P' = (1,0,0,..) and q = 1. By homogeneity, d/dp'_1 becomes E = sum k p_k d/dp_k. This gives, with H = connected single-Hurwitz series sum H_{g,mu} p_mu beta^b/b!:
  (T)  d/dp_1(E H) = exp(H(e^beta o p) - 2H(p) + H(e^{-beta} o p)),  where (e^beta o p)_k = e^{k beta} p_k.
- Checked exactly on all 50 non-zero coefficients (degree <= 6, beta^10; t_toda.out).

(L5) String and dilaton equations for classes pulled back from Mbar_g.

Context: Dubrovin-Yang-Zagier, Moscow Math. J. 17 (2017) 601-633, s.1.3 (empty profile); Giacchetto-Lewanski-Norbury arXiv:2112.11137, Thm 1.2 (our dilaton-leaf formula at R = 0^n).

=== 2. UNIT PARTS: LAMBERT STRUCTURE AND FINITE FORMULA ===

Lemma 4.1. Let Q(s) := sum_m (s^m/m!) <Lambda^vee(1) lam_R (psi^2 lam_1)^m>_{h,r+m}. Then Q is a polynomial of degree <= n0, Q(-1) = Xi^(0)_h(R), and, as formal series in s,
  Q(s) = sum_c (s^c/c!) (1+s)^{-(e+c)} exp(-s(J+c)/(1+s)) A_h(R,1^c).
Proof.
1. psi^2 lam_1 = lam_1 - tau_0 - tau_1 as leaves; expand binomially.
2. Remove the tau_1's by dilaton, with factors (e+c+e0)^(e1) (rising factorial).
3. Remove the tau_0's by string, each with factor J+c (S lam_m = m lam_m, S(1) = 0).
4. Sum the series: sum_{e1} (-s)^{e1} x^(e1)/e1! = (1+s)^{-x}.
5. -psi^2 lam_1 = -psi^2/(1-psi) is the dilaton leaf T, hence Q(-1) = Xi^(0) by Lemma 1 of proof_M_residue. QED

Lemma 4.2 follows by substituting u = s/(1+s).

Prop D. Write Q(-1) = sum_{m<=n0} (-1)^m q_m and exchange the sums. Then sum_{k<=N} (-1)^k [s^k] f(s) = [t^N] f(-t)/(1-t), which is a Laguerre value. QED

Remark 4.4. Only string and dilaton were used, so Lemma 4.2 holds for the CY class Omega. Hence Lemma A at genus h is equivalent to: e^{-uJ} G^Omega_{R,h}(u e^{-u}) = O(u^{-e-1}) at u = oo for every R in the range.

Checks of Lemma 4.1 and Prop D: the whole k = 0 polynomial Q(s) against q_g2_r1, q_g2_r2, q_g3_r1, q_g3_r2 (at s = 0..n0, and the value at s = -1), plus the stored r = 1 polynomials for h = 2, 3, 4. 429 comparisons, 0 mismatches.

=== 3. TODA IN GENUS-GRADED, NORMALISED FORM ===

Genus expansion. Put q_k = beta^{k+1} p_k, so that H = sum_g beta^{2g-2} H_g(q). Then (T) is equivalent to
  (T_g)  d/dq_1 E H_g = [beta^{2g}] exp( sum_{g'>=0, k>=1} beta^{2g'+2k-2} (2/(2k)!) E^{2k} H_{g'} ).

Fixed non-unit profile. Write H_g = sum_R q_R G_{g,R}(z), with z = q_1 and R a multiset of integers >= 2 (R = empty allowed). By ELSV, G_{g,R} = kappa_R^{-1} G_{R,g}, where kappa_R = |Aut R| prod_{j in R} j!/j^j.

Lambert variable. z = T e^{-T}, D := z d/dz = theta/(1-T), theta := T d/dT, and D_J := (theta+J)/(1-T).

Normalised functions. g_{g,R}(T) := e^{-JT} G_{g,R}. By Lemma 4.2, for (g,R) stable:
  (5.1) g_{g,R} = a T^{-e} + O(T^{-e-1}), with a = a_{g,R} = kappa_R^{-1} (-1)^e Xi^(0)_g(R).

Unstable pieces (explicit, from (L3), Lagrange inversion and Lemma 4.2):
- (U1) G_{0,empty} = H_0 = sum_d d^{d-3} z^d/d!, with D^2 H_0 = T and D^3 H_0 = T/(1-T).
- (U2) D H_1 = T^2/(24(1-T)^2), from Lemma 4.2 with R = (1): Q(s) = <(1 - lambda_1)(1 + s kappa_1) lam_1>_{1,1} = s/24.
- (U3) g_{0,(j)} = kappa_j^{-1} (1/j^2 - T/(j(j+1))).
- (U4) g_{0,(j1,j2)} = kappa^{-1}/(j1+j2).
(U3) and (U4) fit (5.1) with e = r-2 and Xi_0(j) = 1/(j(j+1)), Xi_0(j1,j2) = 1/J. For stable genus 0, Xi_0(R) = (J-1)...(J-r+3) (proof_status s.1.3 at h = 0).

Lemma 5.1. For (g,R) != (0,empty),
  theta D_J g_{g,R} = T S'_{g,R},
  S'_{g,R} = [beta^{2g} q_R] exp( sum' beta^{2(g'+k-1)} (2/(2k)!) q_{R'} Y_{g',k,R'} ) - Y_{g,1,R},
where Y_{g',k,R'} := D_{J'}^{2k} g_{g',R'} and sum' runs over (g',k,R') != (0,1,empty) with k >= 1.
Proof.
1. Identities: (J + D)(e^{JT} h) = e^{JT} D_J h; D(J + D)(e^{JT} h) = e^{JT} (D_J - J) D_J h; (D_J - J - T D_J) h = theta h.
2. The term (0,1,empty) of the exponent is D^2 H_0 = T, which gives e^T = T/z.
3. The only appearance of g_{g,R} on the right is the linear term Y_{g,1,R}, with coefficient 1.
4. Both sides are rational in T: H_1 (which has a log) and H_0 enter only through D^{2k}, k >= 1. Equality as power series at T = 0 therefore gives equality as rational functions. QED

=== 4. ASYMPTOTICS AT T = oo ===

Lemma 6.1. Every Y_{g',k,R'} equals L T^{-o} + O(T^{-o-1}), with o = 2(g'+k-1) + r'. The coefficient L is:
- L = (delta')_{(2k)} a_{g',R'} if R' != empty or g' >= 2, where delta' = J' - (2g'-2+r');
- L_{0,k,empty} = -(2k-3)! (k >= 2);
- L_{1,k,empty} = (2k-1)!/12 (k >= 1).
Proof: D_J(c T^{-n} + ...) = -(J-n) c T^{-n-1} + ..., iterated from (5.1) and (U1)-(U4).

Lemma 6.2 (leading-order identity). For (g,R) stable,
  (e+1) delta a_{g,R} = [beta^{2g} q_R] exp(Lambda) - L_{g,1,R},  where Lambda := sum' beta^{2(g'+k-1)} (2/(2k)!) q_{R'} L_{g',k,R'}.
Proof.
1. The left side of Lemma 5.1 is (e+1)(J-e) a T^{-e-1} + O(T^{-e-2}).
2. Every monomial of T S' has order >= sum o_i - 1 = 2g + r - 1 = e + 1.
3. Compare the T^{-e-1} coefficients.

Lemma 6.3 (pure-unit cancellation). u_m := sum_{g'+k=m+1} (2/(2k)!) L_{g',k,empty} = 0 for all m >= 1, i.e. exp(Lambda_empty) = 1.
Proof.
1. u_1 = (1/12)(-1) + 1/12 = 0.
2. For m >= 2, Lemma 6.2 at (m, empty) (J = 0) gives [beta^{2m}] exp(Lambda_empty) = (e+1)(-e)a + e(e+1)a = 0.
3. Induction on m.

=== 5. PROOF OF THEOREM A ===

Claim (V), proved by induction on (g,r) in lexicographic order: if g >= 1, R != empty has parts >= 2, and delta <= -1 or delta = 1, then a_{g,R} = 0.

Setting up. By Lemma 6.3, Lemma 6.2 becomes (e+1) delta a_{g,R} = a sum over factor multisets {(g_i,k_i,R_i)}, where:
- the R_i are non-empty, their union is R, and sum (g_i + k_i - 1) = g;
- the single term (g,1,R) is excluded;
- each term is (symmetry factor) prod_i (2/(2k_i)!) (delta_i)_{(2k_i)} a_i.
Every factor is smaller than (g,R): if s = 1 then k_1 >= 2 and g_1 < g; if s >= 2 then r_i < r.

When a factor can be non-zero:
- if 0 <= delta_i <= 2k_i - 1, the falling factorial vanishes;
- if delta_i <= -1, then g_i >= 1, because genus-0 factors have delta_i = D_i + 2 >= 3; so a_i = 0 by the induction hypothesis.
Hence a non-zero term needs delta_i >= 2k_i for every i.

Counting. Sum_i delta_i = D - 2 sum g_i + 2s and sum_i 2k_i = 2g - 2 sum g_i + 2s. So a non-zero term needs D >= 2g.

Conclusion. For D <= 2g-1 every term vanishes, so (e+1) delta a = 0 with e+1 >= 2. Hence a = 0 when delta = -(something) <= -1 (D <= 2g-3), and when delta = 1 (D = 2g-1).

Base. Genus 0 carries no claim. At genus 1, R = (2): the only other term is (0,2,(2)), with (3)_(4) = 0. QED

Remarks.
- Values with D = 2g-2 (delta = 0) never enter, since they come with the factor (0)_(2k) = 0.
- Those values are non-zero in general: Xi^(0)_2(3) = -1/320, Xi^(0)_3(2,4) = -8/945, Xi^(0)_4(2,2,2,2,3) = -8.

Theorem C is Lemma 6.2 rewritten.
- For multisets: [q_R] exp(sum_{R'} q_{R'} F_{R'}) = (1/|Aut R|) sum_pi prod_B |Aut R_B| F_{R_B}.
- |Aut R_B| a_{g',R_B} = (-1)^{e_B} prod_{j in R_B} (j^j/j!) Xi^(0)_{g'}(R_B).
- The signs cancel because sum_B e_B = e (mod 2).

=== 6. THEOREM B AND HARER-ZAGIER ===

Theorem B, proof.
1. For r = 1, (REC) reads 2h(j-2h+1) Xi_h(j) = sum_{g'<h} (2/(2(h+1-g'))!) (j-2g'+1)_{(2(h+1-g'))} Xi_{g'}(j). This is a polynomial identity in j (it holds for infinitely many j).
2. Put Xi_g = (j-1)_{(2g-2)} Phi_g, so Phi_0 = 1. Then 2h Phi_h = (j-2h) sum_{g'<h} (2/(2(h-g'+1))!) Phi_{g'}.
3. With F = sum Phi_g t^{2g} and C = (2cosh t - 2 - t^2)/t^2, this is t F' = j C F - t(CF)'. Hence W = (1+C) F satisfies t W' = j (C/(1+C)) W.
4. 1 + C = S^{-2}, and 1 - S^2 = -t^2 (log S)''. So W = (S e^{-t S'/S})^j = H^j, and F = S^2 H^j.
5. Since S^2 H^j = H^j - (t/j)(H^j)', Lagrange inversion gives [y^j] x^{j-2h}. QED

Harer-Zagier.
- At R = empty, a_{g,empty} = Int_{Mbar_g} c(T(-log)) = chi(M_g).
- Lemma 6.3 fixes it recursively, and the solution is B_{2g}/(2g(2g-2)).
- Uniformly, L_{g,k,empty} = (2g-1) B_{2g} (2g+2k-3)!/(2g)!.
- u_m = 0 is the identity (x f' - f)(cosh x - 1) = -x^2/2 for f = (x/2)coth(x/2). This is DYZ, Thm 3 and its Corollary.

=== 7. EXACT CHECKS (all passed) ===
- ELSV against mv tables: 81 comparisons (old tables) + 236 (new, genus <= 4, l <= 3).
- B2 against q tables: 429.
- Toda (T): 50 coefficients.
- (REC) against independent Hurwitz-number values (B2), all R with parts >= 2, J <= 9, delta != 0:
  - genus 1: 29; genus 2: 27; genus 3: 24.
  - Included zeros:
    - genus 2: (2), plus D = 3: (4), (2,3), (2,2,2);
    - genus 3: (2), (3), (4), (2,2), (2,3), (2,2,2), plus D = 5: (6), (2,5), (3,4), (2,2,4), (2,3,3), (2,2,2,3).
- Genus-4 points, all verified 0 by Hurwitz numbers:
  - every irreducible point with r >= 3: (2,2,2), (2,2,3), (2,2,4), (2,3,3), (2,2,2,2), (2,2,2,3), (2,2,2,2,2);
  - D = 7 zeros: (2,2,2,5), (2,2,2,2,4).
- Non-zero genus-4 checks: (2,2,7) = 239927/8640 and (2,3,6) = 1707/70, matching (REC).
- Genus-5 points by Hurwitz numbers:
  - in range: (2,2,2) = 0, (3,3,3) = 0, (2,2,2,2) = 0;
  - extra zero: (10) = 0;
  - non-zero, matching (REC) or Theorem B: (11) = 1327385/202752, (12) = 80, (2,2,9) = 102723/256, (9) = -4401/112640.
- (REC) against the Hodge data (q_g2..g4, r <= 2, s = -1): 110 comparisons. These include every in-range r <= 2 genus-4 point and the D = 7 zeros (8), (2,7), (3,6), (4,5).
- Hurwitz numbers against the k = 0 part of hvc_multiprobe_g4b (genus 4, r = 3..5): 8 comparisons, including the delta = 0 values -25/144, -4/25, -1003/6400, -736/675, -417/400, -8.
- Theorem B against (REC), genus <= 8: 109 comparisons.
- Pending: Xi^(0)_5(2^7) by Hurwitz numbers. This run (degree up to 33) was still running when I returned; its output will appear in runs\b2pt_5_2_2_2_2_2_2_2.out. It is covered by Theorem A.

=== 8. NOT PROVED / OBSTRUCTION ===
Not proved: the components k = 1..h-1 of the CY class (full Lemma A).

The method needs, for the CY generating function, an equation of type (T): first order in the unit derivative, with the degree operator, in the charge direction.
- By Marino-Vafa, exp(C) = sum_nu s_nu(i o p) s_nu(q^rho) q^{tau kappa_nu/2}. This is the same 2D Toda tau function tau(P,P',beta,1) with beta = i lambda tau, but at P' = q^rho.
- So Okounkov's (10) holds, but d/dp'_1 at q^rho is not E. It is the vertex-operator expression <0|Gamma_+ E_{-1}(beta) e^{beta F_2} Gamma_-|0>, which involves all p_k and destroys the triangularity used in section 5.
- The P-only modified-KP relation 2 d_{p_2}(F - F_1) + d_{p_1}^2 (F + F_1) + (d_{p_1}(F - F_1))^2 = 0, with F_1 = F(e^beta o p), holds for any P'. It is checked for Hurwitz numbers in t_mkp.out (the opposite sign fails), so it also holds for Marino-Vafa. But it raises the profile (it adds a 2-part), and its leading-order analysis needs explicit cancellations rather than counting.
- The r = 1 CY generating function S(e1 l)^{1+j e2} S(e2 l)^{1+j e1} does not satisfy the analogue of the r = 1 ODE (checked at order lambda^2). So the CY recursion is not a simple deformation of (REC).

=== 9. FILES (lemmaA\B_hurwitz\) ===
Libraries:
- HB.pm: DVV/KW numbers, InputForm evaluator, partitions, Murnaghan-Nakayama characters, T(nu;n,b).
- HBconn.pm: exponential formula, A_h(R,1^c), finite formula (xi0), Q(s).
- Data.pm: k = 0 extraction from the tables by exact interpolation in t = c^2/x.
- Rec.pm: (REC).

Scripts and outputs:

| Script | Output |
|---|---|
| t_basic.pl | t_basic.out |
| t_toda.pl | t_toda.out |
| t_mkp.pl | t_mkp.out |
| b1_check.pl | b1_check.out, b1_check_oldtables.out |
| b2_check.pl | b2_check.out |
| rec_check.pl | rec_check_g1/2/3.out |
| rec_vs_data.pl | rec_vs_data.out |
| closed_r1_check.pl | closed_r1_check.out |
| b2_point.pl | runs\b2pt_*.out |
| g4b_check.pl | g4b_check.out |
| explore_P.pl | explore_P_h2/3.out |

Usage examples: perl rec_check.pl 3 9; perl b2_point.pl 4 2 2 7; perl closed_r1_check.pl 8.

