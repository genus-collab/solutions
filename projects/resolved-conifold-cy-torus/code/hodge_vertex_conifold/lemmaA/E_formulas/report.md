# Route E report (round 1, 2026-10-31)

(Saved by the main session from the workflow result wf_40196105-891; subagents cannot write report files.)

## claims
[
   {
      "depends_on" : "Mumford's GRR formula (Mumford 1983; Faber-Pandharipande 2000 s.2), DVV recursion",
      "evidence" : "compare.pl (0 mismatches in 221 coefficients), check_g4b.pl (26/26), hodge_known.pl, consistency.pl, E_crosscheck_g5_x2_out.txt",
      "statement" : "The Perl engine (Mumford GRR recursion + DVV + AC kappa calculus + forgetful-map recursion with K) computes Xi_h(R) and all components exactly; it reproduces all stored Mathematica data (genus 2,3 r=2; genus 4 r=1; 13 genus-4 points r=3..5), Faber's lambda_g lambda_{g-1} lambda_{g-2}, Getzler-Pandharipande's psi^{g-1} lambda_g lambda_{g-1} and the lambda_g formula for g=2..6, and an independent Mathematica implementation at genus 5 (main session run of E_crosscheck_g5.wls at x=2: 11/11 confirmed).",
      "status" : "computer-verified"
   },
   {
      "evidence" : "proof in summary s.2.2 (kappa_k = pi^*kappa_k + psi^k, f(psi_i) = f(pi^*psi_i) + D_in (Sf)(pi^*psi_i), psi_n D_in = 0, pi_*D_in = 1)",
      "statement" : "Forgetful recursion: W(beta;a) = Sum_{A in beta}[W(beta-A+{e-1};a') - W(beta-A+{e};a')] + [a_n=0] Sum_{i<n} W(beta;a'-e_i) (e = a_n+|A|, kappa_{-1}=0, kappa_0 = 2h-3+n), and its point version with leaves lam_j.",
      "status" : "proved"
   },
   {
      "depends_on" : "engine",
      "evidence" : "check61_engine.pl: 112 (h=5) and 133 (h=6) exact points, >= deg+1 in j and >= h+1 values of sigma",
      "statement" : "Formula (6.1) holds identically in j and tau at genus 5 and genus 6 (Hodge side, no TR).",
      "status" : "computer-verified"
   },
   {
      "depends_on" : "engine; Theorem R and Lemma C of proof_status_mixed_scale.md Addendum 2 for the regularity consequence",
      "evidence" : "mirr_g5.txt, mirr_g6.txt (point evaluation, exact); genus 5 also via full symbolic polynomials (check_stream_out*.txt) and partially by the independent Mathematica run",
      "statement" : "M_irr(5) and M_irr(6) hold: Xi^(k)_h(R) = 0 for every k and every irreducible R (44 points at h=5, 96 at h=6, r up to 9); all control points D = 2h-2 are non-zero. Hence the all-genus regularity theorem holds through genus 6 (with Theorem R and Lemma C).",
      "status" : "computer-verified"
   },
   {
      "depends_on" : "engine",
      "evidence" : "mirr_g7.txt (HODGE_NCH=4, ch_7 branch)",
      "statement" : "Genus 7 (partial): 161 of the 194 irreducible points of M_irr(7) (all with r <= 5, and 12 of 19 with r = 6) vanish in every component; 50 control points are non-zero; no violation. 33 points with r >= 6 were not finished.",
      "status" : "computer-verified"
   },
   {
      "evidence" : "Moebius inversion; Newton differences of 1/(1-j psi); (3.1) of proof_M_residue.md",
      "statement" : "Binomial support: a polynomial vanishes on {j_i >= 1, Sum(j_i-1) <= M} iff its coefficients in the basis prod C(j_i-1,a_i) vanish for |a| <= M; the coefficients equal Xi_h(phi_{a_1},...,phi_{a_r}) with phi_a = a! psi^a/prod_{i=1}^{a+1}(1-i psi), and for the Hurwitz component they are prod a_i! Int c(-R pi_* S_{D'}) prod psi^{a_i}, D' = Sum (a_i+1)x_i (excess Chern numbers of excess r).",
      "status" : "proved"
   },
   {
      "evidence" : "summary s.1.4",
      "statement" : "Kappa reduction for two leaves: Xi_h(j1,j2) = j1(j1-1)P(j1) + j2(j2-1)P(j2) + (j1-1)(j2-1)Y(j1,j2) with P(j) = <<kappa(j)>> and Y = <<kappa(j1)kappa(j2) + Sum_b h_b(j1,j2)kappa_b>>.",
      "status" : "proved"
   },
   {
      "evidence" : "summary s.3 and s.0(f)(2)",
      "statement" : "No formula depending on R only through J (any r) is compatible with both string and dilaton except the lambda_g shape; the top j-degree part of Xi_h(R) is the Witten-Kontsevich r-point function (Dijkgraaf's function for r = 2).",
      "status" : "proved"
   },
   {
      "depends_on" : "engine data",
      "evidence" : "exact overdetermined tests: r=2 m=0..16 (h<=5), m=2..5 (h=6); r=3 eight prefixes (h<=4; (2,2),(3,4) h<=5); r=4 (2,2,2) h<=5; r=5 (2,2,2,2) h<=5; null tests detect Lemma-A-preserving perturbations",
      "statement" : "Conjecture SN (one-leaf universality): for fixed other leaves m (mass M), sigma^h Xi_h(m,j) = Sum_{s=0}^{M-1} prod_{l=s+1}^{N}(M+j-l) [l^{2h}](B^j C^{(m)}_s(l)) with genus-independent series C_s; r = 1 is (6.1).",
      "status" : "conjectural"
   },
   {
      "depends_on" : "SN for the genus in question",
      "evidence" : "two-line argument in summary 0(c)",
      "statement" : "SN implies Lemma A (all h, all r): for R = (m,j) with parts >= 1 and D <= 2h-3, every term of (SN) contains the factor (J-J).",
      "status" : "proved"
   },
   {
      "depends_on" : "SN slice data, HFun.pm",
      "evidence" : "PsiNewton.pl/fitNR.pl: consistent with a unique solution at every order through l^10, i.e. an exact identity for h <= 5; leading A_kk forced by Dijkgraaf (proved given the form); A_11 = e^{-beta}-1 forced by string and checked to l^10",
      "statement" : "Newton-regular form for two leaves: sigma^h Xi_h(j1,j2) = Sum_{k} Sum_{n'=k}^{max(2k-1,0)} (j1)_k(j2)_k prod_{l=n'+1}^{2h-1}(J-l)[l^{2h}](S1S2 B^J A_{k,n'}) with universal series A_{k,n'} = O(l^{2k}sigma^{2k-n'}), A_{00}=1, A_{11} = e^{-beta}-1, leading A_{kk} = (-sigma l^2/8)^k/(k!(2k+1)); it makes Lemma A for r = 2 evident.",
      "status" : "computer-verified"
   },
   {
      "evidence" : "comp_hm1_closed.pl: exact agreement for h = 2..6; found from h <= 5 data, genus 6 is a parameter-free prediction",
      "statement" : "Closed formula for the component next to lambda_g with two leaves: (-1)^{h-1}Xi^(h-1)_h(j1,j2) = prod_{l=1}^{2h-1}(J-l)[l^{2h}](l S' - J S log S) + Sum_{k=1}^{h}(j1)_k(j2)_k prod_{l=2k}^{2h-1}(J-l)[l^{2h}] S G_k, G_k = (-1)^{k+1}(1/k)Sum_n C(2n-k,k-1)c_n l^{2n} = [y^k] of -Sum_n c_n l^{2n}(a^{2n+1}+b^{2n+1}-1)/(2n+1), a+b=1, ab=y. Lemma A for this component (r=2) follows.",
      "status" : "conjectural"
   },
   {
      "evidence" : "identity via prod(1-alpha_i^2) = 1 (Mumford); cancellation pattern computer-verified at h = 3 (comp_hm1.pl)",
      "statement" : "Mumford decomposition of the lambda_{h-1}-type component: c(E^vee) e_{h-1}(alpha(1+alpha)) = lambda_{h-1} - h lambda_h + lambda_h Sum_{n odd} n! ch_n; in the Lemma-A range the two non-trivial pieces Int K lambda_{h-1} prod lam and Int K lambda_h Sum n! ch_n prod lam cancel without vanishing separately.",
      "status" : "proved"
   },
   {
      "depends_on" : "engine",
      "evidence" : "hurwitz_extra_g2..5.txt: h = 2,3,4,5 exhaustive (3,7,15,30 points, all r); zeros_scan.pl: h = 6, r = 2",
      "statement" : "The Hurwitz component Xi^(0)_h vanishes at every R with parts >= 2 and D = 2h-1, and at no R with D = 2h-2 or D = 2h (r >= 2 version of the extra zero j = 2h).",
      "status" : "computer-verified"
   },
   {
      "evidence" : "explore_lines.pl",
      "statement" : "For r = 2 and h <= 6 only the lambda_g component is divisible by a linear form J - l; all other components vanish only at lattice points.",
      "status" : "computer-verified"
   },
   {
      "evidence" : "candK1.pl, candY1.pl, residual_rank.pl, fitD2.pl, fitBB.pl, logfit.pl, comp_hm1_closedR.pl (exact counterexamples listed in summary s.3)",
      "statement" : "Candidates K1 (kappa shift), Y1=Y2 (shifted pull-back residue), D2 (diagonal binomial x Norlund), symmetric bi-binomial with B^J, separable corrections, exponent forms of the slice series, and the pairwise-only extension of the lambda_{h-1} formula to r = 3 are false.",
      "status" : "refuted"
   }
]


## key_formulas
[
   "(SN) sigma^h Xi_h(m_1..m_{r-1}, j) = Sum_{s=0}^{M-1} prod_{l=s+1}^{2h-3+r} (M + j - l) [lambda^{2h}] ( B(lambda)^j C^{(m)}_s(lambda) ),  M = Sum m_i, C_s independent of h, [lambda^0]C_s = delta_{s0}; r = 1: C_0 = S(e1 lambda)S(e2 lambda) (= (6.1))",
   "Borel form of (SN): sigma^h Xi_h(m, f) = [lambda^{2h}] Res_v f^(v + beta(lambda)) Sum_s C_s (N-s)! e^{(M-s)v} (e^v - 1)^{s-N-1} dv",
   "(SN) => Lemma A: R = (m,j), parts >= 1, J <= N, J >= M+1 > s+1 for s <= M-1, so prod_{l=s+1}^{N}(J-l) = 0",
   "(NR, r=2) sigma^h Xi_h(j1,j2) = Sum_{k>=0} Sum_{n'=k}^{max(2k-1,0)} (j1)_k (j2)_k prod_{l=n'+1}^{2h-1}(J-l) [lambda^{2h}] ( S1 S2 B^J A_{k,n'} ),  A_00 = 1, A_11 = e^{-beta} - 1, A_kk = (-sigma lambda^2/8)^k/(k!(2k+1)) + O(lambda^{2k+2})",
   "(-1)^{h-1} Xi^(h-1)_h(j1,j2) = prod_{l=1}^{2h-1}(J-l)[l^{2h}](l S' - J S log S) + Sum_{k=1}^{h}(j1)_k(j2)_k prod_{l=2k}^{2h-1}(J-l)[l^{2h}] S(l) G_k(l),  Sum_k G_k y^k = -Sum_n c_n l^{2n}(a^{2n+1}+b^{2n+1}-1)/(2n+1), a+b = 1, ab = y, c_n = -B_{2n}/(2n(2n)!)",
   "Zhat_j(u) = Sum_h sigma^h Xi_h(j) u^{2h-2}/(2h-2)! = Res_{l=0} l^{-3} (1 + u/l)^{j-1} S1 S2 B^j dl; Lemma A <=> Zhat_R(u) polynomial of degree <= J-1",
   "alpha^(0)_a = prod a_i! Int_{Mbar_{h,r}} c(-R pi_* S_{D'}) prod psi_i^{a_i},  D' = Sum (a_i+1) x_i  (binomial coefficients of the Hurwitz component)",
   "Xi_h(j1,j2) = j1(j1-1)P(j1) + j2(j2-1)P(j2) + (j1-1)(j2-1) Y(j1,j2),  Y = <<kappa(j1)kappa(j2) + Sum_b h_b(j1,j2) kappa_b>>",
   "c(E^vee) e_{h-1}(alpha(1+alpha)) = lambda_{h-1} - h lambda_h + lambda_h Sum_{n odd} n! ch_n(E)",
   "W(beta;a) = Sum_{A in beta}[W(beta-A+{e-1};a') - W(beta-A+{e};a')] + [a_n=0] Sum_{i<n} W(beta;a'-e_i)"
]


## next_steps
[
   "Prove (SN) for r = 2 (or directly the Newton-regular form): show that in the Givental R-matrix graph sum for Omega K the leg dressing R^{-1}(psi) acts on the Borel transform of a leaf as translation by beta(lambda) uniformly in the rest of the graph, with the pairwise contractions (j1)_k(j2)_k arising from edges; start with the sigma-linear order, where the closed formula (e) with the Lucas power sums a^{2n+1}+b^{2n+1} (a+b=1, ab=y) is now known.",
   "Try a Faber-Pandharipande-type localization relation (one free marked point alone at a torus-fixed point, background leaves elsewhere) to derive (6.1) and the free-leaf structure of (SN); (6.1) is a two-weight analogue of FP's ((t/2)/sin(t/2))^{k+1}.",
   "Determine the A_{k,n'} at order sigma^2 and beyond using the second pair of roots (a,b) with a+b=1, ab=y suggested by (e); a closed form for Psi = H/(S1S2) would give a complete r = 2 formula.",
   "Extend the Newton-regular form to r = 3 with a minimal basis (string/dilaton force A_{(k,k,0)} = A^{(2)}_k, A_{(1,1,1),1} = (e^{-beta}-1)^2, A_{(1,1,1),2} = e^{-beta}(e^{-beta}-1)) and test uniqueness.",
   "Finish M_irr(7) (33 points with r >= 6) on a machine without the 2 GB per-process limit; run slice tests at genus 6-7 for more m and prefixes.",
   "Prove the Hurwitz extra layer D = 2h-1 (all r) from the twisted Chern number form of Proposition 1.2 (deg S_{D'} = 1 there)."
]


## obstructions
[
   "No mechanism was found that produces genus-independent coefficient series C_s (resp. A_{k,n'}); this universality is the whole content of (SN)/(NR) and is verified only in the stated ranges (h <= 5, some slices h = 6). A genus-by-genus argument (excess intersection on the jump loci of S_D) can at best give each h separately: its obstruction is the vanishing of excess Chern numbers (excess r) of the K-class -R pi_* S_{D'} in negative degree.",
   "Formula (6.1) itself (the r = 1 case of (SN)) remains unproved for h >= 7 (verified h <= 6 here without TR).",
   "For r >= 3 the symmetric Newton-regular basis is not unique at fixed genus with the general falling-factorial monomials, so only the slice form (SN) was tested there; the triple-contraction coefficients have no closed form yet.",
   "The closed form of A_{k,n'} beyond sigma-linear order is open; at sigma^2 the shortest-shortening coefficients follow binomial patterns only up to lambda^8 and deviate at lambda^10 (A_{4,6}).",
   "The dilaton dressing is a p_1-shift of the Marino-Vafa generating function by a constant ~ sinh(1)/(sigma lambda) outside its radius of convergence, so a formal cut-and-join/character argument does not apply directly; the 33 remaining genus-7 points (r >= 6) exceeded the local per-process memory/time."
]


## report
G:\My Drive\Projects\Old\refGW\hodge_vertex_conifold\lemmaA\E_formulas  (report.md was NOT written: the harness refuses report files from subagents; the full report is the 'summary' field below; all scripts, data and outputs are in this folder)

## requested_checks
[
   "E_formulas/E_crosscheck_g5.wls with xval = None (symbolic x) or another rational x (e.g. 1/3): independent Mathematica genus-5 check of 11 engine predictions (8 Lemma-A zeros, 3 controls); the x = 2 run by the main session already confirmed all 11 (E_crosscheck_g5_x2_out.txt).",
   "E_formulas/E_crosscheck_g4r2.wls (uses hvc_mp_header.m): symbolic Xi[4,2] and Xi[3,3] in Mathematica, then `perl compare_g4r2.pl` compares with engine_h4_r2.txt and engine_h3_r3.txt (the data behind the slice/Newton-regular structures; comparator tested on the genus 2,3 r=2 data: 0 mismatches).",
   "On Stanage (Linux perl, no 2 GB limit): `HODGE_NCH=4 perl mirr_points.pl 7` to finish M_irr(7) (33 points with r >= 6 remaining), and `ORD=2 perl slice_fitJ.pl <m> 7` / `perl slice_fitR.pl 2,2 6` to extend the (SN) tests to genus 6-7.",
   "Housekeeping (needs the user or main session): stop my two runaway early runs of check_lemmaA.pl (msys PIDs 3721 and 3931, Windows PIDs 11816 and 24484, about 2 GB each); killing them was not permitted to me."
]


## route
E: closed formulas and structure of Xi_h(R; tau) for r >= 2 (no topological recursion anywhere)

## summary
ROUTE E REPORT (Lemma A = M_irr(h)). No TR of any kind was used (no spectral curve, EO invariants, loop equations, BKMP, nothing proved mod P2, no TR script, not even as a cross-check). Inputs: Mumford's GRR formula for ch(E), Witten-Kontsevich numbers (DVV), Arbarello-Cornalba kappa calculus, string/dilaton, exact rationals (Perl Math::BigRat). Conventions as in the brief (c = 1, w = (1,tau,-1-tau), sigma = tau(1+tau), v = -1/sigma, Xi = Sum_k v^k Xi^(k), J = Sum j_i, D = J - r, N = 2h-3+r, S(x) = (x/2)/sinh(x/2), e1 = -tau, e2 = 1+tau, S1 = S(e1 l), S2 = S(e2 l), beta = e2 log S1 + e1 log S2, B = e^beta, log S = Sum_n c_n l^{2n}, c_n = -B_{2n}/(2n(2n)!), (j)_k = falling factorial).

0. MAIN RESULTS
(a) Independent exact engine (Hodge.pm, HodgePt.pm, WK.pm): reproduces every stored Mathematica datum (genus 2,3 r=2 symbolic: 46+130 coefficients; genus 4 r=1 symbolic: 45; the 13 genus-4 points r=3..5 of hvc_multiprobe_g4b at x = 2, 1/3: 26/26), Faber's lambda_g lambda_{g-1} lambda_{g-2}, Getzler-Pandharipande's psi^{g-1} lambda_g lambda_{g-1} and the lambda_g formula for g = 2..6, top j-degree = Witten-Kontsevich numbers, string/dilaton across r. The main session independently re-ran my Mathematica script E_crosscheck_g5.wls at x = 2 (E_crosscheck_g5_x2_out.txt): all 11 genus-5 engine predictions (8 Lemma-A zeros, 3 non-zero controls) CONFIRMED, Faber and GP genus-5 values confirmed.
(b) New exact results without TR: (6.1) holds identically in (j, tau) for h = 5 and 6 (112 and 133 points, enough for identities). M_irr(5): all 44 irreducible points, M_irr(6): all 96 irreducible points (r up to 9), every component Xi^(k) zero; all control points D = 2h-2 non-zero (22 resp. 42). Genus 7 (ch_7 branch): 161 of 194 irreducible points verified zero in every component (all points with r <= 5, and 12 of 19 with r = 6), 50 control points non-zero, no violation; the remaining 33 points (r >= 6) were still running at submission. With Theorem R and Lemma C (proof_status_mixed_scale.md, Addendum 2) the all-genus regularity theorem now holds through genus 6 (computer-assisted).
(c) THE STRUCTURE THAT MAKES LEMMA A EVIDENT ("one-leaf universality", SN): for any fixed values m = (m_1..m_{r-1}) >= 0 of all leaves but one (M = Sum m_i) there are power series C^{(m)}_s(l) in Q[sigma][[l^2]], 0 <= s <= M-1, INDEPENDENT OF h, with [l^0]C_s = delta_{s0}, such that for every h:
   sigma^h Xi_h(m, j) = Sum_{s=0}^{M-1} prod_{l=s+1}^{N} (M + j - l) [l^{2h}] ( B(l)^j C^{(m)}_s(l) ).   (SN)
r = 1 is exactly (6.1) (C_0 = S1 S2); m = (0), (1) reproduce dilaton/string (C_0 = S1S2). Equivalent Borel form: sigma^h Xi_h(m, f) = [l^{2h}] Res_v f^(v+beta) K_m(v;l) dv with K_m = Sum_s C_s (N-s)! e^{(M-s)v}(e^v-1)^{s-N-1}: the free leaf enters only through the beta-shift of its Borel transform (the e^{beta S} of (6a)(ii)), against a kernel whose numerator is a polynomial of degree <= M in e^v.
PROPOSITION (SN => Lemma A) [proved]: write R = (m, j) with all parts >= 1. D <= 2h-3 means J <= N; j >= 1 gives J >= M+1; for s <= M-1, s+1 <= M < J <= N, so the factor (J-J) occurs in prod_{l=s+1}^{N}(J-l): every term vanishes. With string/dilaton this is Lemma A at genus h, for all r.
Evidence (each genus an overdetermined test: the residual quotient must drop from degree h+M-1 to M-1, h conditions per power of sigma; null tests with Lemma-A-preserving perturbations are detected): r = 2, m = 0..16, h = 1..5 and m = 2,3,4,5 at h = 6; r = 3 prefixes (2,2),(3,4) h <= 5 and (2,3),(3,3),(1,2),(0,3),(2,4),(1,3) h <= 4; r = 4 prefix (2,2,2) h <= 5; r = 5 prefix (2,2,2,2) h <= 5. No failure.
(d) TWO LEAVES IN CLOSED SYMMETRIC FORM ("Newton-regular", NR) [computer-verified h <= 5, unique solution at every order]:
   sigma^h Xi_h(j1,j2) = Sum_{k>=0} Sum_{n'=k}^{max(2k-1,0)} (j1)_k (j2)_k prod_{l=n'+1}^{2h-1}(J-l) [l^{2h}] ( S1 S2 B^J A_{k,n'}(l) ),
with universal series A_{k,n'} = O(l^{2k} sigma^{2k-n'}), A_{0,0} = 1, A_{1,1} = e^{-beta} - 1 (forced by string; checked to l^10), leading A_{k,k} = (-sigma l^2/8)^k/(k!(2k+1)) + O(l^{2k+2}) (these Gaussian moments are forced by Dijkgraaf's two-point function: proved given NR). The k = 0 term is the J-only candidate (e^{beta S} on both leaves); the corrections are pairwise contractions (j1)_k(j2)_k with shortened Norlund products. Lemma A for r = 2 is evident: at a lattice point with j_i >= 1, J <= 2h-1, either J >= n'+1 (product vanishes) or J <= n' <= 2k-1, so min(j1,j2) <= k-1 and (j1)_k(j2)_k = 0. Equivalent forms checked: H(J, P = j1 j2; l) = S1S2 Psi symmetric with [l^{2n}]H = N_n(J,P)/prod_{l=1}^{2n-1}(J-l) through l^10 (H2_numerators.txt, verifyH.pl: exact identity h <= 5); slice series C_s(m) are polynomials in m divisible by C(m, s+1) (17 slices, degrees <= 10).
(e) CLOSED FORMULA FOR THE COMPONENT NEXT TO lambda_g, two leaves [conjectural; verified h = 2..6, genus 6 being a prediction]:
   (-1)^{h-1} Xi^(h-1)_h(j1,j2) = prod_{l=1}^{2h-1}(J-l)[l^{2h}](l S'(l) - J S(l) log S(l)) + Sum_{k=1}^{h} (j1)_k(j2)_k prod_{l=2k}^{2h-1}(J-l) [l^{2h}] S(l) G_k(l),
   G_k = (-1)^{k+1}(1/k) Sum_n C(2n-k,k-1) c_n l^{2n} = [y^k] G,  G(y;l) = -Sum_n c_n l^{2n}(a^{2n+1}+b^{2n+1}-1)/(2n+1) = -(1/l)[Int_0^{al}+Int_0^{bl}-Int_0^{l}] log S,  a + b = 1, ab = y.
The contraction variable y enters exactly like -sigma enters the CY weights (e1+e2 = 1, e1e2 = -sigma): Lucas power sums of a second pair of roots. Lemma A for this component (r = 2) follows as in (d). For r = 3 the pairwise extension fails: triple contractions are needed (they are consistent but not unique at sigma^1).
(f) Structural facts: (1) for r = 2 and h <= 6 only the lambda_g component is divisible by any (J-l); all other components vanish only at lattice points; (2) the top j-degree of Xi_h(R) is the Witten-Kontsevich r-point function (proved), which rules out every J-type formula and forces Dijkgraaf structure for r = 2; (3) the Hurwitz component Xi^(0) vanishes on the whole extra layer D = 2h-1 (all R with parts >= 2: h = 2,3,4,5 exhaustively (3,7,15,30 points), h = 6 for r = 2) and at no point with D = 2h-2 or 2h (r >= 2 version of the extra zero j = 2h); (4) by Mumford (prod(1-alpha_i^2) = 1) c(E^vee)e_{h-1}(alpha(1+alpha)) = lambda_{h-1} - h lambda_h + lambda_h Sum_{n odd} n! ch_n, so Xi^(h-1) = X_a - h Xi^(h) + X_b with X_a = Int K lambda_{h-1} prod lam, X_b = Int K lambda_h Sum n! ch_n prod lam; in the Lemma-A range X_a = -X_b != 0 (h = 3): neither piece vanishes alone.

1. EQUIVALENT FORMS OF LEMMA A [proved]
1.1 Binomial support: f = Sum_a alpha_a prod_i C(j_i-1,a_i) vanishes on {j_i >= 1, Sum(j_i-1) <= M} iff alpha_a = 0 for |a| <= M (unitriangular Moebius inversion alpha_b = Sum_{a<=b}(-1)^{|b-a|}prod C(b_i,a_i) f(a+1)). So Lemma A <=> every component has binomial support |a| >= 2h-2.
1.2 alpha_a(Xi_h) = Xi_h(phi_{a_1},..,phi_{a_r}) with phi_a = a! psi^a/prod_{i=1}^{a+1}(1-i psi) (Newton differences: Delta^a[1/(1-j psi)] = a!psi^a/prod_{i=0}^a(1-(j+i)psi), induction), and by (3.1) of proof_M_residue.md alpha^(0)_a = prod a_i! Int c(-R pi_* S_{D'}) prod psi_i^{a_i}, D' = Sum(a_i+1)x_i, rank 3h-3-|a|, required Chern degree 3h-3+r-|a|: Lemma A for the Hurwitz component is the vanishing of excess Chern numbers (excess r) of the K-class -R pi_* S_{D'} when deg S_{D'} = |a|-2h+2 < 0.
1.3 Finite genus expansion: Lemma A for all h <=> Zhat_R(u) = Sum_h sigma^h Xi_h(R) u^{N_h}/N_h! is a polynomial of degree <= J-1 for every R with parts >= 1. For r = 1, (6.1) gives Zhat_j(u) = Res_{l=0} l^{-3}(1+u/l)^{j-1} S1S2 B^j dl (genus 1 included: sigma Xi_1(j) = (sigma(j-2)-1)/24).
1.4 Kappa reduction (r = 2): Xi_h(j1,j2) = j1(j1-1)P(j1) + j2(j2-1)P(j2) + (j1-1)(j2-1)Y(j1,j2), P(j) = <<kappa(j)>>, Y = Int_{Mbar_{h,2}} pi^*(Omega K_h) psi_1psi_2 lam_{j1}lam_{j2} = <<kappa(j1)kappa(j2) + Sum_b h_b(j1,j2)kappa_b>> (<<X>> = Int_{Mbar_h} Omega K_h X). Proof: K_{h,2} = pi^*K_h(1-psi_2)(1-psi_1+D_12) (kappa_k = pi^*kappa_k + psi^k), psi_i D_12 = 0, pi_*D_12 = 0, (1-psi)lam_j = 1+(j-1)psi lam_j, pi_*(psi_1^{a+1}) = kappa_{a-1}, pi_*(psi_1^{a+1}psi_2^{b+1}) = kappa_a kappa_b + kappa_{a+b}.

2. ENGINE [recursions proved]
2.1 Mumford (Mumford 1983 s.5; with markings as in Faber-Pandharipande, Invent. Math. 139 (2000) s.2): ch_l(E) = B_{l+1}/(l+1)! [kappa_l - Sum_i psi_i^l + (1/2) xi_*(Sum_{p+q=l-1}(-1)^p psi^p psi'^q)], l odd; restrictions to the boundary additive; kappa_b = pi_*(psi^{b+1}), pi^*kappa_b = kappa_b - psi^b; lambda_k = [deg k]exp(Sum_{l odd}(l-1)!ch_l); ch_l with l <= h suffice. Omega = Lambda^vee(1) prod_i(1 + v alpha_i(1+alpha_i)) (from (1-a/tau)(1+a/(1+tau)) = 1 - a/sigma - a^2/sigma).
2.2 Forgetting points with K (Prop. W): W(beta;a) = Int Omega K kappa_beta prod psi^{a_i} = Sum_{A in beta}[W(beta-A+{e-1};a') - W(beta-A+{e};a')] + [a_n=0]Sum_{i<n} W(beta;a'-e_i), e = a_n+|A|, kappa_{-1} = 0, kappa_0 = 2h-3+n. Proof: K kappa_beta = pi^*(K kappa_{beta-A})(1-psi_n)psi_n^{|A|}; f(psi_i) = f(pi^*psi_i) + D_in (Sf)(pi^*psi_i); psi_n D_in = 0, D_in D_kn = 0, pi_*D_in = 1, pi_*psi_n^c = kappa_{c-1}. Point version with leaves lam_j (S lam_j = j lam_j): HodgePt.pm.
DVV recursion (Dijkgraaf-Verlinde-Verlinde 1991; Witten-Kontsevich).

3. CANDIDATES TESTED AND REFUTED (exact)
- J-only sigma^h Xi = prod_{l=1}^N(J-l)[l^{2h}]S1S2B^{J-c}: impossible [proved: string+dilaton force f_r(J+1)/f_r(J) = J/(J+3-r-2h), i.e. the lambda_g shape only] and fails in the KW top degree.
- K1 kappa shift sigma^h<<kappa_{c_1}..kappa_{c_s}>> = [l^{2h}]S1S2 L_0(prod(e^{beta d}kappa)_{c_i}), L_0 = Int_{Mbar_{0,2h}}K: exact for s = 1 (=(6.1)); s = 2: 2/2, 9/9, 20/20 monomials wrong (h = 2,3,4); s = 3 all wrong.
- Y1 (= Y2) shifted pull-back residue N![l^{2h}]S1S2B^{N+1}Res(1+w)^{r-1}prod((1+w)^{a_i}-1)/a_i (w-(B-1))^{-N-1}: exact r = 1; r = 2: 15/15 (h=2), 19/21 (h=3) points wrong, wrong in the KW top degree.
- separable (rank-one) correction to Y1: residual matrix full rank (4, 7).
- D2 Sum_m C(j1-1,m)C(j2-1,m)prod_{l=2m+1}^{2h-1}(J-l)g_m(J[,P]): inconsistent for all sigma^{>=1} coefficients (h = 2,3).
- symmetric bi-binomial with B^J and shortened products: fails at genus 1 (the -P/24 of Dijkgraaf).
- genus-0 reduction with deformed leaves at one genus-0 vertex: always J-only (Borel form), excluded.
- exponent form S1^{x(e2)}S2^{x(e1)} for the slice series (deg x <= 4): inconsistent.
- Hessian (proof_M_residue s.8): H = j1(j1-1)j2(j2-1)Q, deg Q = 3h-5-k, no line factors except lambda_g; no closed form.
- pairwise-only extension of (e) to r = 3: fails (triple contractions needed).

4. WHAT A PROOF NEEDS
(SN) (or NR for r = 2) is the precise identity; its content is GENUS-INDEPENDENCE of the coefficient series (a genus-by-genus excess-intersection argument can at best give each h separately; obstruction there: excess Chern numbers of -R pi_* S_{D'}). Routes: (i) localization: (6.1) is a two-weight analogue of Faber-Pandharipande's one-point formula Sum t^{2g}k^i Int psi^{2g-2+i}lambda_{g-i} = ((t/2)/sin(t/2))^{k+1} (FP 2000, Thm 2); a relation in which the free leaf is alone at one fixed point would explain "free leaf only through B^j"; (ii) Givental R-matrix: Omega K = R.(trivial), vertex kappa-exponent log R(z) + log(1-z), legs R^{-1}(psi)lam_j; one must show the leg dressing acts on the leaf's Borel transform as translation by beta(l) uniformly in the graph, the contractions (j1)_k(j2)_k coming from edges; (iii) Marino-Vafa: the dressing is a p_1-shift outside the radius of convergence, so it needs the analytic (Lambert-type) continuation, not formal cut-and-join. Obstruction where I stop: no mechanism producing the genus-independent coefficients was found; all statements of 0(c)-(e) are computer-verified in the stated ranges only.

5. FILES (all in E_formulas): engine WK.pm, Hodge.pm, HodgePt.pm, Poly.pm, MParse.pm, Comp.pm, Comps.pm, CPoly.pm, SerS.pm, YFun.pm, HFun.pm; data engine_h<h>_r<r>.txt (components as polynomials, (h,r) up to (2,4),(3,4),(4,5),(5,7),(6,3),(1,3)), mirr_g2..7.txt, hurwitz_extra_g2..5.txt, slice_m0..16.txt, slice6_m2..5.txt, sliceR_*.txt, sliceR5_*.txt, mpoly_out.txt, Cs_polynomials_r2.txt, H2_numerators.txt, PsiNewton_out.txt, expected_g5.m, E_crosscheck_g5_x2_out.txt (main session run); checks check61_engine.pl, check_g4b.pl, consistency.pl, hodge_known.pl, mirr_points.pl, check_stream.pl, explore_lines.pl, zeros_scan.pl, hurwitz_extra.pl, binom_basis.pl, candK1.pl, candY1.pl, residual_rank.pl, fitD2.pl, fitBB.pl, logfit.pl, hessian.pl, slice_fitJ.pl, slice_fitR.pl, slice_fitJ_pert*.pl (null tests), mpoly.pl, mpoly_show.pl, Hsym.pl, Hsym3.pl, Psi.pl, PsiNewton.pl, verifyH.pl, fitNR.pl, comp_hm1.pl, comp_hm1_closed.pl, comp_hm1_closedR.pl, comp_hm1_r3fit.pl, E_crosscheck_g5.wls, E_crosscheck_g4r2.wls, compare_g4r2.pl. Housekeeping: two early runaway runs of my check_lemmaA.pl (msys PIDs 3721, 3931; Windows 11816, 24484, about 2 GB each) could not be stopped by me (kill not permitted); the current check_lemmaA.pl is fixed.

