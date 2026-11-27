# Referee report on route B (round 1, 2026-10-31)

(Saved by the main session from the workflow result wf_40196105-891; subagents cannot write report files.)

## confirmed
[
   "Theorem A (proved modulo ELSV, Okounkov's eq. (10) and Hurwitz's genus-0 formula). For h >= 1 and any multiset R of parts >= 2, Xi^(0)_h(R) = 0 when D <= 2h-3 or D = 2h-1. Every step was re-derived by hand, including the count Sum delta_i - Sum 2k_i = D - 2g, the lexicographic induction (genus-0 blocks have delta = D_i + 2 >= 3) and e+1 >= 2. 27 further zeros were computed from independent Hurwitz numbers, none of them computed by the author: genus 3 (2,2,2,2,2); genus 4, D = 7: (2,2,6), (2,3,5), (2,4,4), (3,3,4), (2,2,3,4), (2,3,3,3), (2,2,2,3,3), (2,2,2,2,2,3), (2^7); genus 5 in range: (2,3), (3,4), (4,4), (2,6), (3,5), (2,2,4), (2,3,4), (2,2,2,2,2), (2,3,3,3); genus 5, D = 9: (2,9), (5,6), (2,2,8); genus 6: (2,2), (5,5), (7), (10), (12). The 7 in-range genus-4 Hodge-data zeros also agree.",
   "Theorem C, the recursion (REC). The referee's implementation (a subset recursion, different code) agrees with the author's Rec.pm on 464 points (g <= 5, J <= 13, r <= 7); all 111 Theorem-A points among them are 0. It agrees with independent Hurwitz numbers on 93 new points with delta != 0 (genus 1-6, r <= 7), and with the Mathematica Hodge data on 126 values (q tables, r <= 2, genus 2-4, extracted at x = -2, -3, ...) plus the 7 genus-4 zeros for r = 3..5.",
   "Lemma 5.1 (normalised Toda equation). It holds as an EXACT identity of Laurent polynomials in w = 1/(1-T), to all orders, in 12 cases (g <= 3, r <= 3). Deliberate perturbations of the inputs make it fail.",
   "Lemma 6.1 and the unstable inputs (U1)-(U4) were re-derived: g_{0,(j)} = kappa^{-1}(1/j^2 - T/(j(j+1))), g_{0,(j1,j2)} = kappa^{-1}/J, D H_1 = T^2/(24(1-T)^2), L_{0,k,0} = -(2k-3)!, L_{1,k,0} = (2k-1)!/12. The genus-0 stable input Xi_0(R) = (J-1)...(J-r+3) agrees in 8 cases.",
   "Lemma 6.3 (pure-unit cancellation): u_m = 0 for m <= 25. L_{g,k,0} = (2g-1) B_{2g} (2g+2k-3)!/(2g)! agrees with DYZ eq. (38). From the referee's Hurwitz engine, chi(M_2) = -1/240 and chi(M_3) = 1/1008, and DYZ's closed forms H_2 and H_3 are reproduced (16 coefficients).",
   "Theorem B. Two independently coded closed forms, [t^{2h}] S^2 H^j and the Lagrange form [y^j] x^{j-2h}, agree with each other and with REC for g <= 10 (341 checks), with the Hodge data for r = 1, g = 2-4 (45 checks), and with Hurwitz numbers at g = 5 (j = 13) and g = 6 (j = 7, 10, 12, 13). The derivation (C/(1+C) = 1 - S^2 = -t^2 (log S)'') was re-done by hand.",
   "Lemma 4.2 (Lambert structure) and Proposition D. Both were re-derived by hand (string/dilaton, rising factorials, Laguerre sum). The structure was checked in an overdetermined way, with 2 extra orders, at all 96 referee points. Prop D and the Lambert extraction agree on every point. The lemma holds for any class pulled back from Mbar_h, hence for the CY class (Remark 4.4).",
   "The specialised Toda equation (T). Okounkov's statement was read in arXiv:math/0004128v1, eqs. (4), (5) and (10) and s.2.7 (no stray constants at n = +-1). The specialisation to P' = (1,0,...), q = 1 was re-derived. Referee re-check with a cut-and-join engine that uses no characters: monomials of degree <= 9 (after d/dp_1), beta^<=18, 1843 coefficients (424 non-zero), 0 mismatches. A sign-flipped right side fails.",
   "ELSV normalisation H_{g,mu} = (b!/|Aut mu|) prod mu^mu/mu! Int Lambda^vee(1)/prod(1 - mu psi), with covers weighted 1/|Aut|: 19 independent comparisons against the Hodge mv tables (extracted at negative x), 0 mismatches. The author's 236 also stand.",
   "The mKP relation for the single-Hurwitz series: referee re-check at degree <= 6, beta^<=12, 390 coefficients, 0 failures; the opposite sign fails (119).",
   "No topological recursion is used: there is no spectral curve, EO recursion, loop equation, BKMP or mod-P2 input, and no TR script. A grep of B_hurwitz finds none. The Lambert substitution z = T e^{-T} is Lagrange inversion, as in DYZ's elementary s.1.3, and the Toda equation is a Hurwitz-theory and infinite-wedge theorem."
]


## issues
[
   {
      "description" : "The status 'computer-verified' covers only the single-Hurwitz series. The author checked degree <= 4 and beta^8. The referee re-checked degree <= 6 and beta^<=12 with independent code: 390 coefficients, 0 failures, and the opposite sign gives 119 failures. The extension in the same statement, '(and, by 2D Toda theory, for any P', in particular for the Marino-Vafa framed vertex)', is an argument, not a computation. That argument is correct: by Okounkov s.2.7, for every P', tau_1(P,P') = q^{1/2} e^{beta/8} tau(e^beta o P, P'), and the first-time-set (mKP) equations of the 2D Toda hierarchy hold for each fixed P'.",
      "location" : "Claims list, item 'modified-KP relation' (status computer-verified); report section 8",
      "severity" : "minor",
      "suggested_fix" : "Split the claim. Mark the Hurwitz-series part computer-verified, giving the range. Mark the Marino-Vafa part proved-modulo-literature, citing Ueno-Takasaki and Okounkov s.2.7 and writing out the exact Hirota equation used."
   },
   {
      "description" : "81 of the ELSV comparisons were made on the first data tables. The main session declared those tables wrong and asked that anything derived from them be discarded (lemmaA\\data\\WARNING_READ_FIRST.txt). The author evaluated the unevaluated wk symbols with his own DVV code. That is plausible, but it means these 81 comparisons are not independent evidence.",
      "location" : "Report s.1 (L1), s.7 and the claims list ('81 comparisons on the first tables')",
      "severity" : "minor",
      "suggested_fix" : "Drop the 81. Quote only the 236 comparisons on the regenerated tables. The referee adds 19 more, independent comparisons (own parser, extraction at negative x, including l = 3 at genus 3 and 4 with parts up to 6), with 0 mismatches."
   },
   {
      "description" : "The 'analogue of the r = 1 ODE' is never written down, and no script or output documents the check. As it stands the remark cannot be verified.",
      "location" : "Section 8 (obstructions), third bullet: r = 1 CY ODE 'checked at order lambda^2'",
      "severity" : "minor",
      "suggested_fix" : "State the ODE that was tested and add the script and output, or remove the remark."
   },
   {
      "description" : "The bullet asserts the obstruction without showing it, and it describes the derivative as less explicit than it is. At P' = q^rho, d/dp'_1 acts diagonally on Schur functions: sum_{mu = nu - box} s_mu(q^rho) = (q^{1/2} - q^{-1/2}) sum_{box in nu} q^{-c(box)} s_nu(q^rho). This is the principal-specialisation branching rule; ref_branch.pl confirms it exactly in 192 cases. The CY analogue of (T) is therefore d/dp_1 (e^{-H} Eq e^{H}) = exp(H(e^beta o p) - 2H + H(e^{-beta} o p)). Here Eq is an explicit diagonal operator, sum_n (-i lambda)^n/n! F_n, where F_n has eigenvalue sum_box c(box)^n (F_0 = E, F_1 = cut-and-join, and so on). The real obstruction is that e^{-H} F_n e^{H} is nonlinear for n >= 1. The loss of an explicit operator is not the obstruction.",
      "location" : "Section 8, first bullet ('d/dp'_1 at q^rho ... destroys the triangularity')",
      "severity" : "minor",
      "suggested_fix" : "Rephrase the obstruction accordingly. Before declaring the route closed, try the T -> oo leading-order analysis component by component (expanding in 1/sigma or 1/tau) on this deformed equation. That is allowed: Marino-Vafa plus Hurwitz/completed-cycle operators."
   },
   {
      "description" : "Theorem A is stated for h >= 1. The string/dilaton relations it relies on (proof_status_mixed_scale.md s.1.2) are stated only for h >= 2.",
      "location" : "Theorem A, 'Consequences' (Conjecture M for all h and r)",
      "severity" : "cosmetic",
      "suggested_fix" : "State the Conjecture-M consequence for h >= 2."
   },
   {
      "description" : "Xi_0(j) = 1/(j(j+1)) and Xi_0(j1,j2) = 1/J are not intersection numbers, since Mbar_{0,1} and Mbar_{0,2} do not exist. They are conventions fixed by (5.1) through (U3) and (U4). The theorem should say so.",
      "location" : "Theorem C (statement) and Rec.pm header",
      "severity" : "cosmetic",
      "suggested_fix" : "Add one sentence to Theorem C defining the r = 1, 2 genus-0 inputs through (5.1), (U3) and (U4)."
   },
   {
      "description" : "In Dubrovin-Yang-Zagier (Moscow Math. J. 17 (2017) 601-633), the relevant text is s.1.3: Theorem 3, eqs. (35) and (38), and the Bernoulli-recursion step in its proof. The Corollary in DYZ s.1.2 concerns fixed degree. Okounkov's statement was verified against arXiv:math/0004128v1 (eqs. (4), (5), (10), s.2.7). The arXiv page lists no journal reference, so 'MRL 7 (2000) 447-453' was not checked by the referee.",
      "location" : "Citations: 'DYZ, Thm 3 and its Corollary'; Okounkov journal data",
      "severity" : "cosmetic",
      "suggested_fix" : "Make the DYZ pointer precise and confirm the MRL reference."
   },
   {
      "description" : "The run has finished. runs\\b2pt_5_2_2_2_2_2_2_2.out reads 'Hurwitz/B2 0 | recursion 0 ok [1454 s]'.",
      "location" : "Section 7, 'Pending: Xi^(0)_5(2^7)'",
      "severity" : "cosmetic",
      "suggested_fix" : "Update the text."
   }
]


## notes
REFEREE REPORT: ROUTE B (Hurwitz component via ELSV, Hurwitz theory and Toda), attack on Lemma A, 2026-10-31.
Not written to referee.md, because the harness rule forbids subagent report .md files. This text is the report.
Scripts and outputs are in G:\My Drive\Projects\Old\refGW\hodge_vertex_conifold\lemmaA\B_hurwitz\referee\ :
- RefHur.pm, RefRec.pm, RefData.pm;
- ref_selftest, ref_toda, ref_mkp, ref_lemma51 (+ _perturb), ref_rec, ref_data, ref_genus0, ref_branch, ref_point (.pl/.out);
- runs\pt_*.out (96 points); make_jobs.pl, run_batch.sh, jobs.txt.
None of this code is shared with HB.pm, HBconn.pm, Data.pm or Rec.pm. No wolframscript was run. Nothing outside lemmaA\B_hurwitz\ was written.

0. VERDICT: SOUND.
- Theorems A, B and C, Lemma 4.2, Proposition D, Lemma 6.3 and the specialised Toda equation (T) all survived every attempt to refute them.
- Theorem A is a complete proof modulo three published theorems: ELSV; Okounkov, Toda equations for Hurwitz numbers, eq. (10); and Hurwitz's genus-0 formula. Their statements and normalisations are confirmed below.
- The report claims nothing false. The issues are minor or cosmetic: status labels, superseded evidence, one undocumented remark, an overstated obstruction, and citation precision.
- Full Lemma A (the components k = 1..h-1 of the CY class) is NOT proved, as the author says himself.
- The brief's TR ban is respected.

1. LITERATURE, CHECKED AT THE SOURCE.
(L4) Okounkov, arXiv:math/0004128v1 (PDF read).
- eq. (4): tau(P,P',beta,q) = Sum_lambda q^{|lambda|} e^{beta f_2(lambda)} s_lambda(P) s_lambda(P').
- eq. (5): f_2 is the content sum.
- eq. (10): d^2 log tau/dp_1 dp'_1 = q tau(..,e^beta q) tau(..,e^{-beta} q)/tau^2.
- After eq. (9): tau(P,P',beta,q) = tau(q p_1, q^2 p_2, .., P', beta, 1).
- s.2.7: tau_n = q^{n^2/2} e^{n(4n^2-1)beta/24} tau(.., e^{n beta} q). At n = +-1 the constant factors multiply to q, so there are no stray constants.
- The report's specialisation is correct. With P' = (p'_1,0,..), log tau = H((q p'_1) o P), and d_{p_1} d_{p'_1} log tau at p'_1 = q = 1 equals d_1 H + E d_1 H = d_1 (E H).
(L1) ELSV.
- Normalisation: H_{g,mu} = (b!/|Aut mu|) prod mu^mu/mu! Int Lambda^vee(1)/prod(1 - mu psi), where H counts covers weighted by 1/|Aut|, i.e. (1/d!) times the number of transitive tuples.
- Re-checked by hand on H_{0,(1,1,1)} = 4.
- Re-checked against the Mathematica mv tables with the referee's own parser, extracting the k = 0 part at x = -2, -3, ... (the author used x = 1, 2, ...): 19 comparisons, including l = 3 at genus 3 and 4, e.g. (6,6,6) at genus 4. 0 mismatches.
(L3) Hurwitz genus 0: A_0(R,1^c) = (J+c)^{r+c-3} (22 values), and Xi_0(R) = (J-1)...(J-r+3) (8 values).
DYZ (Moscow Math. J. 17 (2017) 601-633), s.1.3.
- Theorem 3, eq. (35): kappa_{g,2g-2} = B_2g/(2g(2g-2)).
- eq. (38): H_g^{(n)} ~ (2g-1) B_2g (2g+n-3)!/((2g)! T^{2g+n-2}).
- This is exactly the report's R = empty argument and its L_{g,k,0}. The phrase 'extends DYZ to arbitrary profile R' is accurate.
GLN (arXiv:2112.11137): the abstract is consistent with the context citation. It is not used in any proof.

2. INDEPENDENT RE-DERIVATIONS (by hand). All correct.
- Lemma 4.1/4.2:
  - psi^2 lam_1 = lam_1 - tau_0 - tau_1;
  - dilaton first, giving (e+c+e0)^(e1);
  - then string: S lam_m = m lam_m and S(1) = 0, giving (J+c)^{e0};
  - Sum (-s)^n x^(n)/n! = (1+s)^{-x}, then u = s/(1+s);
  - degree bound n0 = 3h-3+r.
  Only string and dilaton are used, so the lemma holds for any class pulled back from Mbar_h.
- Prop D: [t^N](1-t)^{-(a+1)} exp(xt/(1-t)) = Sum_i C(N+a, N-i) x^i/i!.
- Lemma 5.1:
  - (J+D)(e^{JT}h) = e^{JT} D_J h;
  - D(J+D) = e^{JT}(D_J - J)D_J;
  - D_J - J - T D_J = theta;
  - z e^T = T;
  - g_{g,R} appears only through the linear term Y_{g,1,R};
  - e^{J'T} can be factored out by the J-grading.
- Lemma 6.1, with the unstable inputs:
  - g_{0,(j)} solves (theta+j)g = (1-T)/j;
  - P(eps) = (eps-1)/24 at (1,(1)) gives D H_1;
  - L_{0,k,0} = -(2k-3)! and L_{1,k,0} = (2k-1)!/12.
- Lemmas 6.2 and 6.3, the counting step and the lexicographic induction on (g,r): no gap.
- (REC) from 6.2: the multiset-to-set-partition factor |Aut R_B|, the cancellation of kappa, and the sign parity are all right.
- Theorem B: the substitution Phi = Xi/(j-1)_(2h-2) holds as rational functions, so no divisibility is assumed. Then 1+C = S^-2, 1 - S^2 = -t^2 (log S)'', S^2 H^j = H^j - (t/j)(H^j)', and Lagrange inversion.

3. INDEPENDENT COMPUTATIONS. All exact (Math::BigRat). 0 failures anywhere.
Engine:
- disconnected counts by cut-and-join on class types (no characters) and by the referee's own Murnaghan-Nakayama character sums;
- connected counts by the integer E-logarithm Hhat(mu) = Zhat(mu) - Sum C(d-1,d1-1) C(b,b1) Hhat(m1) Zhat(mu/m1). This is a different algorithm from the author's moment-cumulant recursion;
- p_0 extracted by an overdetermined Lambert fit (2 extra orders) AND by Prop D.
Results:
(a) Self-tests, 165 checks:
- cut-and-join = characters for all mu |- d <= 9 and b <= 12;
- brute force in S_4 and S_5;
- genus 0 and genus 1 values;
- DYZ H_2 and H_3 (16 coefficients);
- chi(M_2) = -1/240 and chi(M_3) = 1/1008 obtained as p_0 at R = empty.
(b) Toda (T): degree <= 9 after d/dp_1 and beta^<=18. 1843 coefficients (424 non-zero), 0 mismatches. The author checked 50. A sign-flipped right side fails.
(c) Lemma 5.1 exactly in w = 1/(1-T), to all orders: 12 cases, namely (1,2), (1,3), (1,5), (2,2), (2,3), (2,4), (3,2), (1,(2,3)), (2,(2,2)), (2,(2,3)), (1,(2,2,2)), (3,(3,4)). The identity holds in all of them. Perturbing the H_1 input or g_{0,(j)} breaks it.
(d) Lemma 6.3: u_m = 0 for m <= 25, plus L-formula consistency for g <= 8 (60 checks).
(e) REC: the referee's implementation agrees with the author's on 464 points (g <= 5, J <= 13, r <= 7). The 111 Theorem-A points among them are all 0.
(f) Theorem B: both forms agree with REC for g <= 10 (341 checks).
(g) 96 new Hurwitz-number points. 93 agree with REC, and the 3 with delta = 0 match the Hodge data. 27 Theorem-A zeros are listed in 'confirmed'. Sample non-zero values: g4 (3,3,5) = 195303/8960; g4 (2,2,2,6) = 2997/14; g5 (13) = 3463937/6144; g6 (13) = 9908066657/100638720. Every r >= 3 point with J = 10-12 (genus 1) or J = 10-11 (genus 2, 3) is included.
(h) Hodge data (Mathematica, independent of Hurwitz theory), 213 checks:
- 19 ELSV;
- REC = k = 0 part of q_g2..g4 at s = -1: 126 values, 26 of them zeros, including all in-range r <= 2 genus-4 points and the D = 7 zeros;
- Theorem B = data: 45;
- Hurwitz/Lambert = data: 8, including the delta = 0 values -1/320, -1/180, -65/24192, -8/945, -39/4480 and -9/280;
- genus-4 r = 3..5 zeros: 7;
- in addition, genus-4 delta = 0 Hurwitz values -25/144, -1003/6400 and -4/25 equal the g4b data.
The author's own outputs were re-read and agree with his claims. The pending genus-5 (2^7) run has finished and gives 0, matching REC.

4. ISSUES: see the issues list. Nothing is fatal or major.

5. OBSTRUCTIONS.
- The obstruction is real for a verbatim transfer: at P' = q^rho, d/dp'_1 is not the degree operator.
- It is, however, explicit and diagonal: eigenvalue (q^{1/2}-q^{-1/2}) Sum_{box in nu} q^{-c(box)}, by the principal-specialisation branching rule (ref_branch.pl, 192 exact checks).
- So the CY analogue of (T) is d_{p_1}(e^{-H} Eq e^{H}) = exp(Delta H), with Eq = Sum_n (-i lambda)^n/n! F_n built from completed-cycle operators, which the brief allows.
- The genuine difficulty is that e^{-H} F_n e^{H} is nonlinear for n >= 1. Whether the T -> oo counting survives component by component (in 1/sigma) was not tried. That is the natural next step for the full Lemma A.
- The mKP remark is correct.
- The r = 1 CY ODE remark is undocumented.
- For Theorem A itself the Toda route is already short. I see no simpler path.

6. NO-TR COMPLIANCE.
- No TR, spectral curve, loop equation, BKMP or mod-P2 result appears in the arguments or the code (grep finds nothing).
- Cited project results are the TR-free ones: the kappa form (Lemma 1), the log-tangent identity (L), string/dilaton s.1.2, and the genus-0/lambda_g formula s.1.3.
- The Lambert variable is used only through Lagrange inversion.

## refuted
[]


## report
NOT WRITTEN. The subagent harness rule forbids report .md files, so the complete referee report is the 'notes' field. Referee scripts and outputs: G:\My Drive\Projects\Old\refGW\hodge_vertex_conifold\lemmaA\B_hurwitz\referee\

## route
B: Hurwitz component of Lemma A via ELSV, Hurwitz theory and Okounkov's Toda equation (lemmaA\B_hurwitz)

## verdict
sound

