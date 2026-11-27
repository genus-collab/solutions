# Referee report on route A (round 1, 2026-10-31)

(Saved by the main session from the workflow result wf_40196105-891; subagents cannot write report files.)

## confirmed
[
   "A1 normalisation. The referee's own MV engine (rmv.pm, mod p) reproduces the exact Hodge tables mv_g1..4_l1..3 (with unit legs) and q_g2..4_r1,r2 at s = -1, 2 and 5, at random tau: 996 comparisons, 0 mismatches. This confirms the factor (-1)^g sigma^(g+l-1) Pi_mu / |Aut mu|, the weights (1, tau, -1-tau), the leaf weight c and the c-normalisation c^-(3h-3+2r). The LLZ factor V_nu equals the product over boxes of 1/[h], checked by hand.",
   "A2 (Prop 2.1, Laguerre formula) re-derived by hand: s Tt = s(lam_1 - 1 - psi); dilaton gives (1+s)^-(2h-2+r+a+b); string gives (J+a)^b; then the Laguerre generating function at t = -s and the partial-sum identity. Confirmed numerically for every s-coefficient, not only at s = -1.",
   "Cor 2.2 (u = s/(1+s), Lambert form) and Cor 2.3 (Kaehler-twisted framed vertex) re-derived: substituting x' = x/(1+s), p'_k = e^{-uk} P_k/(sigma x Pi_k), p'_1 = e^{-u} s/(sigma x) into (1.1) gives exactly sum_h hbar^(h-1) Q_h(s;R)/|Aut R|. The reformulation of Lemma A as a Laurent-degree bound follows.",
   "Section 2.4 (lambda_g residue formula) at genus 2: the residue gives (7/5760)(j^2 + s + 3sj + 3s^2), which equals the x -> 0 limit of sigma^2 q_g2_r1.",
   "Theorem A (Xi_h(2) = 0 for all h and all tau) re-derived. On p = p_1 only, d_tau C = x p_1^2 d_2 C, which gives (g-1+a+sigma d_sigma) M_g(1^a) = a(a-1) M_g(2,1^(a-2)). The a = 0 case follows from the CY degree argument (the only symmetric cubic surviving on w1+w2+w3 = 0 is e_3). Then y d_y = s(1+s)^2 d_s and Cor 2.2 give the polynomial identity in s; set s = -1. Direct MV computation also gives Xi_h(2) = 0 mod p for h = 2..7.",
   "Theorem B re-derived: d_tau C is O((1+s)^m) while every cut and join term is O((1+s)^(m-1)), m = 2g-2+r. Checked exactly from the tables alone at genus 2 for R = (3) and R = (2) at four values of tau, with non-zero terms (referee\\out\\misc_g2.out). The R = (2) relation at genus 2 coincides with KP1 at genus 2.",
   "Theorem C (tau = -1/2, odd parts): every coefficient is invariant under tau -> -1-tau, so the residue in sigma is legitimate. Pi_m(-1/2) = 0 exactly when m is even, and the zero is simple. All pole terms are odd-part terms in the range. Base case h = 3, T = (3) checked by hand.",
   "Theorem C' (tau = -1/3): the logic via the closed two-leaf form is correct, given Theorem D and Lemma A in lower genera.",
   "Fact 3.5 (the MV generating function is a KP tau function): standard GL(infinity)-orbit theorem. Independently confirmed at the dilaton point: KP1 and KP2 at genus 2 with all terms, using only the Hodge tables and hand-derived genus-0/1 values, vanish exactly at six values of tau. Perturbing Xi_1(1) breaks KP1 (referee\\out\\kp12_g2.out).",
   "Theorem D re-derived from the hbar-expansion of the differential Fay identity (leading order exp(D(z)D(w)F_0) = 1 - d_1(D(z)-D(w))F_0/(z-w); the [z^-2 w^-2] coefficient gives dispersionless KP). For each factor the x-order equals the (1+s)-order. The excess identity (sum of factor excesses = x-order - J(R) - w) holds term by term, with genus-0 excess <= -2 and genus-1 excess <= 0, so at e >= 1 only the genus-G-linear terms survive.",
   "Cor E (merge formula, reduction of Lemma A(G) to the one-leaf case) follows from Theorem D. Z(3) = 2 checked by hand.",
   "Theorem F reproduced by a referee implementation written from scratch (labelled-subset Leibniz rule, w-degree exp recursion, all-T unknowns or Lagrange interpolation, primes 1e9+9 and 998244353). Full-system nullity is 6, 6, 0, 0, 0, 0 at G = 3..8. Two-leaf subsystem nullity is 11, 11, 7 at G = 4, 5, 6 and 0 for every G = 7..20. All identical to the route. Genus-0 self-test: 690 relations, 0 mismatches.",
   "The mod-p lifting argument of Theorem F is correct: integer tau give relations whose denominators are coprime to p, and the mod-p Vandermonde left kernel is the reduction of the saturated rational kernel.",
   "The route's own scripts reproduce when re-run: a3_kp_nullity.pl 4 and 5, a3_kp_nullity2.pl 6 and 7 (two-leaf), t_twoleaf_closed.pl, and a3_kp_genus0_test.pl at a new tau = 3/11 (420 relations).",
   "Induction base, checked independently and without KP by direct MV + Prop 2.1 (mod p): all 6 genus-3 and all 18 genus-4 points of the range vanish (up to five parts). The controls at D = 2h-2 are non-zero.",
   "Conclusion 'Lemma A for h <= 20' corroborated without KP by direct MV computation (mod p, random tau). Genus 5: all 44 points of the range vanish and all 14 controls at D = 8 or 9 are non-zero. Genus 6: one-leaf c <= 10, all 20 two-leaf points and 8 three-leaf points vanish. Genus 7: one-leaf c <= 12 and 7 two-leaf points vanish. All tested controls just outside the range are non-zero.",
   "Section 3.8 observations confirmed exactly from the tables: Q_h(s;R) has a simple zero at s = -1 at every irreducible point of genus 2 and 3; the lambda_1 component Xi^(1)_2(1/2, 5/2) = -1/3072 (and Xi^(0) = -1/6144, Xi^(2) = 0).",
   "Obstructions confirmed as real: genus 3 and 4 are not closed by KP (nullity 6 reproduced); TLR(G) fails for G <= 6 (two-leaf nullities 11, 11, 7) and holds for 7 <= G <= 20; the all-genus gap is TLR(G) for G >= 21.",
   "TR audit: no topological recursion, spectral-curve input, loop equations, BKMP, dictionary P2 or mod-P2 result anywhere in A_mv or in the referee's checks. Y(z) is the genus-0 dispersionless Lax inverse computed from MV genus-0 data, not a TR input."
]


## issues
[
   {
      "description" : "The KP hierarchy is not on the brief's explicit list of allowed inputs. It is an integrable-hierarchy theorem, not topological recursion. The content-product series sum_nu r_nu s_nu(t) s_nu(t*) lies in the GL(infinity) orbit of the vacuum (Takasaki 1984; Orlov-Scherbin TMP 128 (2001); same mechanism as Okounkov's Toda equations for Hurwitz numbers). So it is mathematically solid and TR-free. Still, every result after section 3.4 depends on it, and its admissibility is a ruling only the user can make. The author flagged this correctly.",
      "location" : "Fact 3.5 (KP for the MV generating function), used by Theorems C', D, Cor E, Theorem F and the claim 'Lemma A for all h <= 20'",
      "severity" : "minor",
      "suggested_fix" : "Ask the user to rule explicitly on admitting KP and its dispersionless limit. Cite the GL(infinity)-orbit statement precisely. Independently of KP, genus 5 can be proved directly (see the last issue)."
   },
   {
      "description" : "As stated this is false. At order x^(2G) the full hierarchy also contains terms made only of genus-0 factors with more derivatives, for example F_{0,1111} in KP1 at G = 1. Generally these are terms with sum over factors of (d_i - 2) = 2G. They are harmless: by excess homogeneity every term has total excess e, genus-0 factors have excess <= -2, so they cannot occur when e >= 1. The conclusion of Theorem D is unaffected.",
      "location" : "Report section 3.5 (derivation of Theorem D): 'every other term contains a factor of genus 1..G-1'",
      "severity" : "minor",
      "suggested_fix" : "Replace the sentence with: 'every other term either contains a factor of genus 1..G-1 or consists of genus-0 factors only; the latter have total excess <= -2 and cannot occur for e >= 1'."
   },
   {
      "description" : "Genus-0 factors whose legs are all unit legs have several exponential sectors. For example d^2G0_empty/dy^2 = (e^{2u} - 1)/2 and dG0_empty/dy = e^u(1 - u/2) - 1 - y/2. Here u = s/(1+s) tends to +infinity as s -> -1 from below and to -infinity from above, so the two one-sided limits are dominated by different sectors, and 'leading order as s -> -1' is ambiguous. The correct and rigorous procedure is: write each relation coefficient as a finite sum over W of e^{uW} R_W(u) with R_W rational; the e^{uW} are independent over C(u); project onto the top sector e^{u(J+al+be)}; then take the leading coefficient at u = infinity. The checks (route's a3_hirota_check.pl and the referee's ref_kp12_g2.pl) implicitly use exactly this projection.",
      "location" : "Report sections 3.1 and 3.5 (uniform order rule, 'leading order at s = -1')",
      "severity" : "minor",
      "suggested_fix" : "State the sector projection explicitly. Write out Lemma 3.1 including the r = 1, 2 cases and the genus-1 factor without parts: Q_1(s;(j)) = (j - 1 - 1/sigma + s)/24, which gives Xi_1(1) = -(1+1/sigma)/24 and Xi_1(2) = -1/(24 sigma)."
   },
   {
      "description" : "The two scripts are by the same author and share the same E-recursion and multiset-binomial design, so they are not independent in the sense that would catch a shared conceptual error. This is now remedied by the referee's re-implementation, which reproduces every nullity.",
      "location" : "Theorem F evidence: 'two independent implementations (a3_kp_nullity.pl, a3_kp_nullity2.pl)'",
      "severity" : "minor",
      "suggested_fix" : "Cite the referee reproduction: ref_kp_null.pl and outputs referee\\out\\full_g5.out, full_g6.out, full_g7g8.out, twoleaf_g4to8.out, twoleaf_g9to12.out, twoleaf_g13to20.out."
   },
   {
      "description" : "'0 mismatches in 2070 comparisons': the run summaries add up to 2170 (35+105+245+56+147+588+63+315+392+77+147).",
      "location" : "Claim A1, evidence field",
      "severity" : "cosmetic",
      "suggested_fix" : "Change 2070 to 2170."
   },
   {
      "description" : "Prop 3.2 is said to have been checked exactly at g = 2 (a <= 6), but runs\\ holds only cj_line_g3.out and cj_line_g4.out. report.md is also absent from A_mv because the harness refused to write it, while the scripts refer to 'report section ...'.",
      "location" : "Theorem A / Prop 3.2 evidence; A_mv\\runs",
      "severity" : "cosmetic",
      "suggested_fix" : "Re-run 'perl a3_cj_line.pl 2 6 > runs/cj_line_g2.out'. Have the main session save the route summary as A_mv\\report.md."
   },
   {
      "description" : "Not a defect but a gap that is easy to close. For any fixed genus h there is a KP-free route: MV Schur side, then Prop 2.1, then exact arithmetic at h values of sigma (proof_status 1.7(a)). The referee's mod-p version covers all 44 genus-5 points in about 4 minutes. That is evidence, not proof; done in exact arithmetic it would make Lemma A(5) independent of the KP ruling and cross-validate Theorem F.",
      "location" : "Scope of Theorem F / Lemma A for h <= 20",
      "severity" : "minor",
      "suggested_fix" : "Run the direct MV computation in exact rational arithmetic (mvlib's integer W_nu series, or a BigRat port of referee\\rmv.pm) at 5 values of sigma, for the 44 genus-5 points."
   }
]


## notes
REFEREE REPORT, Route A (Marino-Vafa, full tau-dependence), 2026-10-31.

I did not write referee.md: the subagent rules forbid report files, the same rule that blocked the author's report.md. This text is the report; the main session can save it as A_mv\referee.md.

Scripts (all in G:\My Drive\Projects\Old\refGW\hodge_vertex_conifold\lemmaA\A_mv\referee\):
- rparse.pm: own parser for the Mathematica tables.
- rmv.pm: own MV Schur-side engine mod p (beta-set Murnaghan-Nakayama, W_nu series, cumulants by set partitions and by a multiset recursion, Laguerre transform).
- ref_xi.pl: table checks and direct Xi_h(B).
- ref_kp_null.pl: own Theorem D relations and nullities.
- ref_kp12_g2.pl: KP1 and KP2 at the dilaton point from the tables.
- ref_misc_g2.pl: Theorem B, simple-zero and off-lattice checks.

Raw outputs are in ...\A_mv\referee\out\*.out. No wolframscript was run. No TR file or TR-derived result was read or used.

VERDICT: sound. Nothing mathematical was refuted. Every proved / proved-modulo-literature / computer-verified claim survived re-derivation or independent recomputation.

The main result, Lemma A for every h <= 20, is correct with one dependency: the KP hierarchy of the MV generating function (Fact 3.5). It is a TR-free theorem but is not on the brief's explicit list, so the user should rule on it. Independently of KP, direct MV + Prop 2.1 computations (mod p, so evidence rather than proof) confirm Lemma A at every point of the range in genus 3, 4 and 5, and at large samples in genus 6 and 7.

1. NORMALISATIONS (A1, A2)
- My engine reproduces the mv and q tables for genus 2-4, one to three leaves, with unit legs, at s = -1, 2 and 5: 996 comparisons, 0 mismatches. This fixes (1.1) with (-1)^g sigma^(g+l-1) Pi_mu/|Aut mu|, the weights (1, tau, -1-tau) and the factor c^-(3h-3+2r).
- Re-derived by hand: Prop 2.1, Cor 2.2, Cor 2.3 (the Kaehler-twisted vertex substitution is exact) and the lambda_g residue formula at genus 2.
- LLZ's V_nu equals the product over boxes of 1/[h], checked on (1), (2), (1,1), (2,1).
- Count typo: 2070 should read 2170.

2. THEOREM A, proof re-derived
- On the p_1-line the cut-and-join gives d_tau C = x p_1^2 d_2 C, hence (g-1+a+sigma d_sigma) M(1^a) = a(a-1) M(2,1^(a-2)).
- The a = 0 case holds because the integral of Lambda(w1)Lambda(w2)Lambda(w3) is C_g w1 w2 w3 on w1+w2+w3 = 0.
- With y d_y = s(1+s)^2 d_s this gives s^2 Q(s;2) = (1+s)[...], an identity of polynomials; s = -1 gives Xi_h(2) = 0. Correct.
- This is the full-tau version of the main session's Hurwitz-component argument.

3. THEOREMS B, C, C'
- B: in the (1+s)-counting the tau-derivative term (including the y d_y piece from y = sigma x p_1) is one order more regular, so it drops. Exact check from the tables at genus 2, R = (3) and R = (2), four values of tau, all with non-zero terms.
- C: all coefficients are functions of sigma (Pi_m(-1-tau) = (-1)^(m-1) Pi_m(tau)); the pole bookkeeping at sigma = -1/4 and the induction are correct.
- C': correct given Theorem D and Lemma A in lower genera.

4. KP AND THEOREM D
- Re-derived from the hbar-expansion of the differential Fay identity. The leading order is the Takasaki-Takebe form (the [z^-2 w^-2] coefficient is dispersionless KP), and the F_G-linear part is its linearisation.
- At the dilaton point, each factor's x-order equals its (1+s)-order.
- Excess identity: summing over factors, the total excess equals (x-order) - J(R) - w, the same for every term. Genus-0 factors have excess <= -2 and genus-1 factors <= 0. So for e >= 1 only genus-G-linear terms with genus-0 co-factors remain.
- Write-up fixes needed (minor):
  (a) Pure genus-0 dispersive terms do occur at order x^(2G), contrary to the text; they are excluded by excess homogeneity.
  (b) State the projection onto the top exponential sector. All-unit-leg genus-0 factors have several sectors, e.g. (e^{2u}-1)/2, and the one-sided limits at s = -1 differ; the projection is rigorous because the e^{uW} are independent over C(u).
  (c) Write out Lemma 3.1 and the genus-1 values. Mine: Q_1(s;(j)) = (j - 1 - 1/sigma + s)/24.
- Independent all-term test from the tables only: at genus 2,
  KP1: 6 sigma Xi_1(1)(1+Xi_1(1)) + 12 Pi_2^2 Xi_2(2,2) - 36 Pi_3 Xi_2(3) = 0
  KP2: sigma Pi_2 [24 Xi_1(2) + 6(Xi_2(2) + 2 Xi_1(1) Xi_1(2))] + 6 Pi_2 Pi_3 Xi_2(2,3) - 24 Pi_4 Xi_2(4) = 0
  Both hold exactly at six values of tau, and perturbing Xi_1(1) breaks KP1. Their Xi_2-linear parts coincide with Theorem D's operator.
- The author's full DJKM check is a separate code path and a genuine validation.

5. THEOREM F, independent reproduction
- Full system: nullity 6, 6, 0, 0, 0, 0 at G = 3..8 (2976 unknowns at G = 8).
- Two-leaf subsystem: nullity 11, 11, 7 at G = 4, 5, 6 and 0 for every G = 7..20.
- These match the route exactly, using different primes, seeds, sample counts and formulation.
- The mod-p lifting argument is correct, so the result is rigorous given KP and the h <= 4 base.
- The base (main-session data) is confirmed by direct MV at genus 3 and 4.

6. KP-FREE DIRECT TEST
- Genus 5: all 44 points in the range vanish.
- Genus 6: one-leaf c <= 10, all 20 two-leaf points and 8 three-leaf points vanish.
- Genus 7: one-leaf c <= 12 and 7 two-leaf points vanish.
- Every control just outside the range is non-zero.
- Section 3.8 observations confirmed: simple zeros at s = -1, and -1/3072 at (1/2, 5/2).

7. TR AUDIT: clean.
- Inputs are MV (LLZ/OP), cut-and-join, string/dilaton, lambda_g, Mumford-side tables, KP and its dispersionless limit.
- Y(z) is a genus-0 dispersionless Lax function computed from MV data. The user may want to know that in TR language it is the framed C^3 curve, though it is not used as a TR input.

8. OBSTRUCTIONS AND SIMPLER PATH
- The obstructions are real. The only all-genus gap is TLR(G) for G >= 21; genus 3 and 4 are not reached by KP.
- Simpler path for fixed genus: direct MV + Laguerre in exact arithmetic at h values of sigma is KP-free and cheap at genus 5.
- For all genera nothing simpler was found. One suggestion: write Y through Lagrange inversion (Pi_k = binom(k(tau+1) - 1, k - 1), an allowed tool) to study pole cancellation at tau = -k/a structurally.

Housekeeping: runs\cj_line_g2.out is missing; report.md is missing.

## refuted
[]


## route
A_mv: Marino-Vafa formula, full tau-dependence (G:\My Drive\Projects\Old\refGW\hodge_vertex_conifold\lemmaA\A_mv)

## verdict
sound

