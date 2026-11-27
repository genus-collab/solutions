# Referee report on route C (round 1, 2026-10-31)

(Saved by the main session from the workflow result wf_40196105-891; subagents cannot write report files.)

## confirmed
[
   "Theorem H holds: the k = 0 component vanishes for all h >= 0, r >= 1, j_i >= 0, J >= 1, D <= 2h-3. The referee re-derived every step of the proof: F, P and phi; Lemma 2 via the Chiodo r=1 GRR (AC kappas); the projective-bundle formula with Fulton's conventions; the localized top Chern class; Y != C exactly when D <= 2h-3; the strata bound J_Y - e - 1; the degree count. It is correct modulo FP Theorem 3 / s.0.1, which states exactly what is used.",
   "Theorem H settles sub-lemma A2 for every h and every r, including legs j_i = 0. Nothing in it uses TR or a mod-P2 result.",
   "An independent, TR-free computation of Xi^(0) from Hurwitz numbers (Burnside + Murnaghan-Nakayama + exponential formula + ELSV + string/dilaton) agrees exactly with the author's tables at all 56 lattice points tested (h = 2,3,4; r = 1,2; parts >= 1; D <= 2h-1).",
   "Theorem H holds at 29 in-range points outside the author's component tables, all computed by Hurwitz numbers and all zero: genus 3 {2,2,2}; all 7 genus-4 irreducible points with r = 3,4,5; genus 5 j = 2..8 and (2,2), (3,3), (2,7), (3,6), (4,5); genus 6 j = 2..10.",
   "The Hurwitz component is non-zero at D = 2h-2 at every tested point. Six genus-4 r>=3 boundary values (-25/144, -4/25, -1003/6400, -736/675, -417/400, -8) agree exactly between the Hurwitz-number computation and the Mathematica multiprobe data.",
   "H' holds at all 266 (j,b) pairs (referee's own code). It is sharp in b: with one b_i = j_i + 1, 70 of 84 pairs are non-zero.",
   "Corollary H0 holds: the genus-0 formula (J-1)...(J-n+3) checks out for n = 3..10, J = -4..16 (168 pairs, 0 mismatches), and the lambda_g component b_h (J-1)...(J-N) matches the tables at h = 2,3,4, r = 1,2.",
   "C3(a) holds: Int_{Mbar_{2,1}} c(T(-log)) lambda_1 lambda_2 lam_2 = 1/2880 and Int c(T(-log)) lambda_2 lam_2 = -1/2880 (analytic, with Int_{Mbar_2} lambda_1 lambda_2 = 1/5760 read off psi_g2).",
   "Prop. L and the logic of Prop. L1 are correct. The component-wise string and dilaton relations hold as polynomial identities for h = 2,3,4 and every k.",
   "The C3(f) refutation is reproduced by an independent test that evaluates at lattice points: inconsistent for r = 2 at (h,k) = (2,1), (3,1), (3,2); consistent for r = 1 with the author's constants.",
   "The Gauss-Bonnet check (9/9) and the normalisation of Lemma 2 are confirmed.",
   "All of the author's scripts were re-run and give byte-identical outputs. A separate exact decomposition in u = -1/sigma of all six q-tables and the three symbolic multiprobe outputs gives the same components.",
   "Byproduct: the conjectural Hurwitz limit of (6.1) (proof_M_residue.md 6a(iii)) matches the TR-free values at genus 5 (j = 9, 10, 11) and genus 6 (j = 11, 12)."
]


## issues
[
   {
      "description" : "H' is presented as a strengthening, but it is equivalent to Theorem H on the box. By partial fractions, span{psi^b/prod_{k<=j}(1-k psi) : 0<=b<=j} = span(lam_0..lam_j), and the excluded monomial b = j is the only one with a lam_0 x...x lam_0 component. So H' at j holds exactly when Theorem H holds at every m with 0 <= m_i <= j_i, m != 0. Checked numerically: the number of H' pairs equals the number of box points at every (h,r) (3, 14, 10, 58, 21, 160). H' therefore adds nothing beyond Theorem H with legs j_i >= 0 allowed.",
      "location" : "Theorem H, 'Strong form (H')' (statement and claims list)",
      "severity" : "minor",
      "suggested_fix" : "Call it the 'box form' or an equivalent reformulation (a corollary), not a strong form."
   },
   {
      "description" : "Boissy is cited at 's.2.3'. The local-parameter count 2g+r+s-2 (r zeros, s poles) is Lemma 3.6 / Remark 3.8 (s.3.3.4 of arXiv:1211.4951). Farkas-Pandharipande is cited without a theorem number. What is used is FP Theorem 3 (s.0.4): for mu with a negative part, H~_g(mu) in Mbar_{g,n} has pure dimension 2g-3+n. Also used is the s.0.1 statement that the open locus H_g(mu) in M_{g,n} is nonsingular of dimension 2g-3+n, extending Polishchuk; parts m_i = 0 are allowed. Only this open-locus statement is needed.",
      "location" : "Proof step (3) and the depends_on literature of Theorem H",
      "severity" : "minor",
      "suggested_fix" : "Cite FP Thm 3 and s.0.1 (and Polishchuk) explicitly, and correct the Boissy location."
   },
   {
      "description" : "Getzler-Pandharipande derived the lambda_g formula from the Virasoro constraints, a loop-equation-type input that the 'no TR at all' rule excludes. Faber-Pandharipande (Ann. of Math. 157 (2003)) prove it by virtual localization, which is allowed.",
      "location" : "Corollary H0, depends_on ('Faber-Pandharipande / Getzler-Pandharipande')",
      "severity" : "minor",
      "suggested_fix" : "Cite Faber-Pandharipande only."
   },
   {
      "description" : "Prop. L assumes r >= 2k, but full localization needs r_+ >= 2k when legs j_i = 0 are present. Prop. L1 is labelled 'proved', but it rests on the same strata-dimension bound as Theorem H. Its phrase 'Y carrying ALL markings' should read 'all markings with j_i >= 1'.",
      "location" : "Route C3(b) Proposition L and C3(c) Proposition L1",
      "severity" : "minor",
      "suggested_fix" : "Use the hypothesis r_+ >= 2k. Relabel L1 as proved-modulo-literature and fix its wording."
   },
   {
      "description" : "The summary says 286 (j,b) pairs. check_Hprime_out.txt gives 3+14+10+58+21+160 = 266, which the referee re-runs confirm.",
      "location" : "Summary, EXACT CHECKS / evidence for H'",
      "severity" : "cosmetic",
      "suggested_fix" : "Correct the count to 266."
   },
   {
      "description" : "Lemma V is correct; for DM stacks the cutting argument needs powers of line bundles that descend to the projective coarse space. It is not needed: W' is proper, so Int_{W'} f^*alpha . gamma = Int_T alpha cap f_*gamma, and this is 0 because A_{deg alpha}(T)_Q = 0 when deg alpha > dim T (projection formula, Fulton Prop. 2.5(c) / Thm 3.2(c)).",
      "location" : "Proof step (5), Lemma V",
      "severity" : "cosmetic",
      "suggested_fix" : "Replace Lemma V by this one-line projection-formula argument."
   },
   {
      "description" : "'At D = 2h-2 ... the bound fails by exactly one' is imprecise. When Y = C, the margin e + r_{Z,+} >= 1 drops to 0, and the holomorphic stratum adds one dimension, so deg alpha - dim T >= -1 instead of >= 1. 'Theorem H at genus 0 is the closed formula' is also loose: the closed formula is a genus-0 computation, and Theorem H gives only its vanishing part.",
      "location" : "Remarks after Theorem H; wording of Corollary H0",
      "severity" : "cosmetic",
      "suggested_fix" : "Reword both sentences."
   },
   {
      "description" : "Several points the argument relies on are left implicit. W is chosen with constant generic dual graph Gamma and support Y. The generic point of W' lies over the open stratum (smooth components), so only open-strata dimensions are needed. Free poles that collide with each other or with special points only lower the count. Zero-order marked points are allowed in the strata.",
      "location" : "Proof steps (3)-(4), presentation",
      "severity" : "cosmetic",
      "suggested_fix" : "Add two sentences making these explicit."
   }
]


## notes
REFEREE REPORT - Route C (intersection theory on Mbar_{h,r}), attack on Lemma A. 2026-10-31. The 'no TR at all' rule was followed throughout: no topological recursion, nothing proved only mod P2, no wolframscript. referee.md was NOT written, because the harness refuses report files from subagents; this field is the report. Scripts and outputs are in G:\My Drive\Projects\Old\refGW\hodge_vertex_conifold\lemmaA\C_intersection\referee\: rpoly.pm, ref_ucomp.pl, ref_elsv.pl, ref_elsv_tests.pl, ref_Hprime.pl, ref_genus0.pl, ref_fit.pl, ref_stringdil.pl, ref_hurwitz_limit61.pl. Outputs are in out\, out\r3plus\ and out\g4b\; re-runs of the author's scripts are in rerun\. out\elsv_new3r.txt may be empty: that sequential run was replaced by the per-point runs in out\r3plus\.

VERDICT: SOUND.
- Theorem H is correct, and its proof is complete modulo literature that is used correctly. Theorem H says the Hurwitz component (k = 0) of Lemma A vanishes for every h, every R in the range and any r.
- I re-derived every step and could not break it.
- A separate, TR-free computation from Hurwitz numbers reproduces the author's values exactly. It also confirms Theorem H at 29 points outside the author's component tables (genus 3 to 6, r up to 5).
- All other claims with status proved, proved-modulo-literature or computer-verified hold. The defects are minor or cosmetic.
- Sub-lemma A2 is settled. Lemma A still needs the components 1 <= k <= h-1.

1. THEOREM H, STEP BY STEP
(0) F, P and phi.
- On every stable curve C, omega_log has degree 2g_v-2+n_v > 0 on each component, so H^0(C, omega_log^{-1}) = 0.
- Hence F = R^1 pi_* omega_log^{-1} is locally free of rank N. The fibre of phi: P -> F is the connecting map, and its kernel is H^0(C, S_D).
- c(F) = Lambda^vee(1) K follows from Chiodo's r=1 GRR. Going from s=0 to s=-1 changes only the kappa term: B_{m+1}(0) - B_{m+1}(-1) = (m+1)(-1)^m gives -Sum kappa_m/m.
- This matches K = exp(-Sum kappa_k/k) of Lemma 1, with Arbarello-Cornalba kappas.
- Toy checks: on Mbar_{0,4}, c_1 = -1 = chi(M_{0,4}); on Mbar_{1,1}, -lambda_1 - kappa_1 = -1/12.
- P has graded pieces L_i^{-m}, m = 1..j_i, so identity (S) holds, also for j_i = 0.

(1) Projective bundle and localization.
- With Fulton's conventions (P = lines, s = 1/c), p_*(c_N(p^*F (x) O(1)) zeta^j) = c_{N-J+1+j}(F-P).
- deg prod p_{j_i} = J - r_+, so only c_m with m >= N-J+1 occur.
- Excision (Fulton Prop. 1.8) and the localized top Chern class (Prop. 14.1), with Vistoli's rational Chow groups of DM stacks, give Xi^(0) = Int_beta [prod p_j/(1-zeta)]_{J-1}.

(2) Support of a section.
- A section of the line bundle S_D on a nodal curve is a family of sections of nu^*S_D = omega^{-1}(Sum (j_i-1)x_i - Sum branches) that agree in the node fibres.
- If Y = C, then Sum d_v = D - (2h-2). So in the range D <= 2h-3 we get Y != C and e >= 1.
- Every component of Y then carries a branch, so each eta_v = 1/s_v has a pole there and only strictly meromorphic strata occur.

(3) Dimension bound.
- div eta_v = Sum(j_i-1)x_i - Sum q - E.
- The pointed curve together with at most d_v free poles lies in a strictly meromorphic stratum, of projectivized dimension 2g_v - 3 + (number of points).
- This gives dimension at most J_v - n_v^bd - 1 per component and J_Y - e - 1 in total. Correct.
- In fact each Y-Y node costs one more dimension; this is not needed, but it confirms the 'saturated' picture in L1.
- Sanity example, h=3, R={2,2}: the jump locus contains the whole rational-tail divisor, T is a surface in P^3, and zeta^3 vanishes on T.

(4)-(6) Factorisation and degree count.
- The generic point of W' lies over the open stratum, so dim f(W') <= J_Y - e - 1 < J_Y + r_{Z,+} - 1 <= deg alpha. Correct.
- Lemma V is correct, but the projection formula is enough.

Literature (checked online):
- Farkas-Pandharipande, J. Inst. Math. Jussieu 17 (2018) 615-672, Theorem 3 (s.0.4): pure dimension 2g-3+n when mu has a negative part.
- FP s.0.1: the open locus H_g(mu) is nonsingular of dimension 2g-3+n, and parts m_i = 0 are allowed. This is all that is used.
- Boissy, Comment. Math. Helv. 90 (2015): the dimension count is Lemma 3.6 / Remark 3.8 in s.3.3.4, not s.2.3.
- Signs and normalisations cannot affect a vanishing statement. Xi^(0) is the sigma -> infinity coefficient, symmetric under tau <-> -1-tau, and the c-normalisation only rescales.

2. H'
- H' is correct but equivalent to Theorem H on the box {0 <= m_i <= j_i, m != 0}.
- The number of pairs is 266, not 286. All 266 vanish.
- H' is sharp in b: with one b_i = j_i+1, 70 of 84 pairs are non-zero (for example h=2, j=(0,3), b=(1,0) gives -1/32).

3. COROLLARY H0
- Confirmed: 168 (n, J) pairs with 0 mismatches.
- Cite Faber-Pandharipande (virtual localization) for the lambda_g formula, not Getzler-Pandharipande (whose route is via the Virasoro constraints).

4. COMPUTER CHECKS
(a) Re-runs: all of the author's outputs are byte-identical.
(b) Exact decomposition in u = -1/sigma, without sampling in x: all six q-tables and the three symbolic multiprobe outputs give the same components.
(c) Hurwitz-number computation (ref_elsv.pl):
- It uses Burnside/Frobenius with Murnaghan-Nakayama characters, the exponential formula and ELSV; the dilaton leaves are handled by writing T = (1+psi) - lam_1 and applying the string and dilaton equations.
- It uses no Mumford formula, no Witten-Kontsevich numbers, no tables and no TR.
- Against the author's tables: 56 of 56 points agree exactly.
- New in-range points, all zero (29 in total):
  * genus 3: {2,2,2};
  * genus 4: {2,2,2}, {2,2,3}, {2,2,4}, {2,3,3}, {2,2,2,2}, {2,2,2,3}, {2,2,2,2,2};
  * genus 5: j = 2..8, and (2,2), (3,3), (2,7), (3,6), (4,5);
  * genus 6: j = 2..10.
- Boundary values D = 2h-2, all non-zero:
  * genus 3 {2,2,3}: -1/28;
  * genus 4: {2,2,5} -25/144, {2,3,4} -4/25, {3,3,3} -1003/6400, {2,2,2,4} -736/675, {2,2,3,3} -417/400, {2,2,2,2,3} -8. These six equal the k=0 parts of the Mathematica values in hvc_multiprobe_g4b_out.txt.
  * genus 5: j=9 -4401/112640, (5,5) -261355/1419264, (2,8) -1024/4455, (4,6) -72/385;
  * genus 6: j=11 -221948177/603832320.
- Byproduct: the conjectural Hurwitz limit of (6.1), proof_M_residue.md 6a(iii), matches genus 5 (j = 9, 10, 11) and genus 6 (j = 11, 12) exactly.
- Observation: the Hurwitz component also vanishes at D = 2h-1 at every tested point. This covers r = 1, 2 at h = 2..4, r = 3 at h = 3 and 4, and the extra zeros j = 2h at h = 5, 6. It is non-zero at D = 2h-2 and at D >= 2h. So the extra zero j = 2h is really a whole hyperplane D = 2h-1, and Theorem H does not explain it.

5. ROUTE C3
- (a) Confirmed.
- (b) Prop. L is correct, but it needs r_+ >= 2k.
- (c) Prop. L1 is correct; its status should be proved-modulo-literature.
- (d) (T1) is plausible and correctly labelled conjectural.
  * The obstruction space at (Z, q') is T_{Mbar_{g_Z,1}}(-log), and the plumbing term lands on the vertical line L_{q'}^vee.
  * That leaves pi^*T_{Mbar_{g_Z}}(-log) (x) O(1), whose rank matches the excess dimension, and the lambda_1-term then vanishes.
  * To be proved: the Y-part is regular, and the excess bundle has no psi_{q'} part.
- (e) Conjectural. The vacuous-recursion risk is real.
- (f) The refutation is reproduced independently.
- C2 (r-th roots): the negative observation is correct.

6. NO-TR AUDIT
Clean. Nothing reads the TR files, and nothing uses a mod-P2 statement.

7. OBSTRUCTIONS AND SIMPLER PATHS
- Theorem H needs no simpler proof.
- For 1 <= k <= h-1 the obstruction is genuine. The Z-part of the jump locus is a whole moduli space, so a dimension count cannot kill the Hodge classes of Z. A per-support excess (virtual) class is needed; (T1) is its first case.
- Mumford's formula for lambda_1 does not give a shortcut: the term Int c(T(-log)) kappa_1 prod lam is a new functional.
- Theorem H supplies one of the h linear conditions per lattice point (proof_status_mixed_scale.md s.1.7(a)).

Sources: https://arxiv.org/abs/1508.07940 (Farkas-Pandharipande, Theorem 3 and s.0.1, read via https://ar5iv.arxiv.org/html/1508.07940); https://ar5iv.arxiv.org/html/1211.4951 (Boissy, Lemma 3.6 / Remark 3.8); https://ems.press/journals/cmh/articles/13114 (Boissy, CMH 90 (2015) 255-286).

## refuted
[
   "The summary's '286 (j,b) pairs' for the H' check is wrong: the author's own output check_Hprime_out.txt gives 266.",
   "The Boissy reference 's.2.3' for the dimension of meromorphic strata points to the wrong place: the count is Lemma 3.6 / Remark 3.8 in s.3.3.4 of arXiv:1211.4951."
]


## report
referee.md not written: the harness refuses report files from subagents. The full report is in the notes field. Referee scripts and outputs: G:\My Drive\Projects\Old\refGW\hodge_vertex_conifold\lemmaA\C_intersection\referee\

## route
C: intersection theory on Mbar_{h,r} (jump loci, r-th roots, Chiodo classes)

## verdict
sound

