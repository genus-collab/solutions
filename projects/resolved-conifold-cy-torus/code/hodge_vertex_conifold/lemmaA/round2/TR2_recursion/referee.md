# Lemma A round 2: referee TR2 (saved by the main session from the workflow journal, 2026-10-31)

Subagents cannot write report files, so this file holds the agent's structured return value: the fields first, then the full report text (report_text) verbatim.

## confirmed

- Claim 1: ABDKS formula (17) (Prop. 3.1 / Thm 2.8, arXiv v3) applies to C_tau for tau not in {0,-1}, and the transcription (2.1)-(2.2) is exact. Checked against v3 for: S(u) = u^-1(e^{u/2}-e^{-u/2}), (15)-(17), Remark 3.2, Def. 2.4, Remark 2.6, the standing assumptions, the Prop. 3.1 setting, and Lemma 5.5, whose common-pole clause covers v = oo (dX simple pole, dv double pole). The vital points are v = 0 and v = -1; v = oo is not vital.
- Claim 2 (Theorem 1, Lemma D for all tau, g, n): the scaling argument was re-derived with ord as the m-adic valuation, and no gap was found. Numerically, the referee's independent TR satisfies the bound at 7 new tau values (494 zeros, 0 violations), and the formula gives 246 more zeros.
- Claim 3 (Corollary 2): the logic is correct mod P2 (positive binomials at integer tau, polynomiality in 1/sigma, D(R) <= 2h-3 iff |R| <= 2h-3+r).
- Claim 4 (Theorem 3, (5.1)-(5.3)): integration by parts, the sinh form, the product form and the residue at oo were re-derived. The referee's own implementation of (5.1) agrees with independent TR up to genus 6.
- Claim 5 (Corollary 4): the term-by-term identity was re-derived and the g = 1, 2 hand checks agree with Route E tables. It is valid mod P2 and is equivalent to P2 for r = 1 (mod ABDKS).
- Claim 6 (Norlund identities): re-derived (Norlund polynomial residue, Appell property, exact Gamma-ratio expansion, Laurent-polynomial extension in tau).
- Claim 7 (two-point formula (5.4)), together with the general-n integrated formula, which the referee verified against independent TR for n = 3 (g <= 4), n = 4 (g <= 3) and n = 5 (g <= 3).
- Claim 8 (Proposition 5): re-derived for integer tau; the local form was verified independently at tau = 3, 5/2 and -7/3, with a failing negative control.
- Claim 9: the reduction statements are correct; the tau = 1 global identity omega_{g,1} = Brk(v,1/v)/((y(v)-y(1/v))dX) is correct.
- Claim 10: P2 for r = 1 holds for h <= 7 mod ABDKS. The referee's direct TR also confirms P2 at h = 6 (r = 1, 2) without ABDKS.
- Claim 11: Lemma D holds directly in exact TR; reproduced with the independent engine at tau = 5/2, -7/3, 2/9, 4, 1/2, 3/7 and -2/5, including the tau = 3/7 run whose output was missing.
- P2 itself (not a TR2 claim) agrees with exact TR-free Hodge tables in 1397 comparisons at 7 tau values, (h,r) up to (6,2), (5,2), (4,4) and (3,4); this tests the sign, the factor sigma^(h-1), the binomial leaves and repeated parts.
- Rules respected: TR is used only on the Hodge curve C_tau, never on the conifold x C^2; KP is not used; no wolframscript was run.

## issues

### [1]

#### description

Both corollaries are labeled 'proved-modulo-literature', but they depend on P2 (Lemma B). P2 is a project lemma that has not been written out (route TR1); it is not literature. The text correctly says 'conditional only on P2 and [ABDKS]', but the status labels do not match it.

#### location

Claims 3 and 5 (Corollary 2, Corollary 4): status field

#### severity

minor

#### suggested_fix

Relabel both as 'conditional on P2 (mod P2)'. Record Lemma A for all h, (OL) for all G, Conjecture M (vanishing) for all h and r, Corollary R in all genera, and (6.1) for all h as mod P2. Note that (6.1) for all h is equivalent to P2 for r = 1 (mod ABDKS).

### [2]

#### description

The reduction is proved: Lemma D(g,1) is equivalent to the vanishing, up to order (2g-3) tau, of the tau-section of Phi_g at v = 0. The assertion that this section is 'not controlled by induction' and is 'the precise obstruction' is an informal remark, yet it carries the status 'proved'.

#### location

Claim 9 / report s.6 line (1)

#### severity

minor

#### suggested_fix

Split the claim: mark the reduction 'proved' and the obstruction statement 'heuristic/remark'.

### [3]

#### description

The report says the identity holds 'for all tau in its equivalent local form' (the sigma-even polar part at v = tau vanishes), but it does not write the argument for non-integer tau. That argument is a one-line rationality step. The statement is true: the referee checked it independently at tau = 5/2 and -7/3, with a failing negative control.

#### location

Claim 8 (Proposition 5), extension to non-integer tau

#### severity

minor

#### suggested_fix

Add: the polar coefficients of Phi_g at v = tau are rational functions of tau and vanish at every integer tau >= 2, so they vanish identically.

### [4]

#### description

Several cited results have no stored output. The Lemma D run at tau = 3/7 (g <= 4, n <= 3, 82 coefficients) has no output file. The output of tau1_values.pl (the genus-5 and genus-6 tau = 1 comparison) was not stored. The summary says '21 points' where claim 1 says '22 points'.

#### location

report s.7 checks / out folder

#### severity

cosmetic

#### suggested_fix

Store the outputs. The referee's reruns are in referee\out\hodge_tau3_7_G4N3.txt and referee\out\author_tau1_values_rerun.txt; both agree with the cited values.

### [5]

#### description

For h = 6 and 7 the TR side comes from Theorem 3, i.e. modulo ABDKS, not from a TR computation. The label 'computer-verified' hides that dependency.

#### location

Claim 10 (one-leaf P2 cross-check, h <= 7)

#### severity

cosmetic

#### suggested_fix

Label the claim 'mod ABDKS, plus Route E computer verification'. The referee's direct TR runs confirm P2 at h = 6, r <= 2 (tau = 5/2, -7/3) without ABDKS.

### [6]

#### description

arXiv v3 (20 Jan 2025), the version read, lists 'several corrections and clarifications' and postdates the IMRN 2024 publication. Equation numbers (15)-(17), (92) and Lemma 5.5 are those of v3.

#### location

report s.2.1 citation

#### severity

cosmetic

#### suggested_fix

Cite arXiv:2312.16950v3 explicitly next to the IMRN reference.

### [7]

#### description

The 'one-line proof of Theorem 3.1 for n = 1 given the bound in lower genera' needs the (g-1, 2) bound, i.e. the n = 2 case in lower genus.

#### location

Claim 9, tau = 1 consequence

#### severity

cosmetic

#### suggested_fix

State that the induction hypothesis includes n = 2 at genus g-1.


## notes

I read ABDKS (arXiv:2312.16950v3) through a web-fetch summarising tool, which refuses verbatim quotation. It returned formula (17), (15), (16), Remark 3.2, the definition of S, Def. 2.4, Remark 2.6, the standing assumptions, the Prop. 3.1 setting, and the statement and proof structure of Lemma 5.5, including the (92) exponent 1-s+2(m+s-2)Sum r_i. All of these match the author's use. I did not re-derive the proof of Lemma 5.5 line by line. The numerical agreement of (17) with independent EO TR in the x-expansion at v = oo, for n <= 5 and g <= 6, is strong evidence that the theorem applies at v = oo.

All referee scripts are in G:\\My Drive\\Projects\\Old\\refGW\\hodge_vertex_conifold\\lemmaA\\round2\\TR2_recursion\\referee\\ (RefMod.pm, RefTR.pm, RefXY.pm, RefHodge.pm, cmp_hodge.pl, cmp_xy.pl, qle_check.pl, spot_values.pl, t_ref.pl). Outputs are in referee\\out\\. No files outside the referee folder were modified; the author's tau1_values.pl was only re-run read-only.

The main remaining dependency for Lemma A in all genera is P2 (Lemma B, route TR1). P2 is now strongly supported numerically but still unproven.

## refuted


## route

TR2 (round 2): topological recursion on C_tau through the x-y swap. Lemma D for every tau, g, n; Lemma A, (OL) and (6.1) in all genera (mod P2)

## verdict

sound

## report_text

# Referee report: route TR2, round 2 (x-y swap on C_tau)

Folder: HVC\lemmaA\round2\TR2_recursion (HVC = G:\My Drive\Projects\Old\refGW\hodge_vertex_conifold). My scripts and outputs are in TR2_recursion\referee\ and referee\out\.

## 0. Verdict

**Sound.** I tried to refute every claim and could not.

- **Claims that are fully established.** On the TR side of C_tau, the following hold with exactly the dependencies the author states:
  - Theorem 1 (Lemma D for every tau, g, n) and the closed formulas (Theorem 3, (5.4)) are correct modulo one published theorem, ABDKS Thm 2.8 / Prop. 3.1.
  - The Norlund identities are proved.
  - Proposition 5 is proved.
- **Hypotheses and transcription of ABDKS.** I checked both against the paper (arXiv:2312.16950v3). They are correct.
- **The Hodge-integral consequences are conditional on P2 (Lemma B).** These are Corollary 2 (Lemma A for all h, (OL), Conjecture M vanishing, Corollary R) and Corollary 4 ((6.1) for all h).
  - P2 is a lemma of this project that has not yet been written out. It is not literature.
  - So these two corollaries should carry the status "conditional on P2", not "proved-modulo-literature" (issue 1).
  - The author's text says this correctly; only the JSON status labels are off.
- **Independent numerics.** I wrote my own exact (mod p) EO recursion and my own implementation of the integrated n-point formula. Compared with each other and with Route E's TR-free Hodge tables through P2, they give:
  - 1397 TR-vs-Hodge comparisons, 660 formula-vs-TR comparisons and 616 formula-vs-Hodge comparisons, all with 0 mismatches;
  - Lemma D with 0 violations;
  - 7 values of tau that the author did not use, up to genus 6, and up to n = 5 points.

  In particular, P2 now has direct exact confirmation (no ABDKS involved) at:
  - (h, r) = (6, 1) and (6, 2);
  - (5, <= 2);
  - (4, <= 4);
  - (3, <= 4).

  These confirmations include the sign, the factor sigma^(h-1) that does not depend on n, the binomial leaves, and repeated parts (no automorphism factors needed).

## 1. The cited theorem (ABDKS, arXiv:2312.16950, v3)

I read v3 (HTML) through the web-fetch tool and checked the following.

- **Definition of S.** Notation s.1.5 defines S(u) = u^(-1)(e^(u/2) - e^(-u/2)). This is the author's S (s_k = 1/(4^k (2k+1)!) in the code).
- **Standing assumptions (s.2.2-2.3).**
  - Sigma is compact.
  - dx and dy are meromorphic and may have residues.
  - The zeros of dx and of dy are simple and their zero sets are disjoint.
  - dy is regular at the zeros of dx, and vice versa.
- **Definition 2.4 (LogTR-vital point).** A point is LogTR-vital if it is a simple pole of dy with non-zero residue at which dx has no pole.
  - Remark 2.6: a simple pole of dy that is also a pole of dx is not vital, and it contributes nothing for g > 0.
  - Remarks 2.12 and 3.4: LogTR = TR when y = z is affine.
- **Prop. 3.1 setting.**
  - Sigma = CP^1.
  - y = z, with z = oo a pole of dy.
  - x = Sum_i alpha_i^(-1) log(z - a_i) + R(z), where the a_i are the vital points of the dual, 1/alpha_i are the residues of dx there, and R is "not assumed rational".
  - dy has no zeros.
- **Formula (17) for 2g-2+n > 0:**

      omega^(g)_n = [hbar^(2g-2+n)] (-1)^n prod_i ( Sum_r (d_i 1/dx_i)^r [u_i^r] e^{-u_i(S(u_i hbar d_{y_i}) xhat_i - x_i)} ) (-1)^(n-1) Sum_{sigma n-cycle} prod_i sqrt(dw_i dwbar_sigma(i))/(w_i - wbar_sigma(i))

  - (16): w = e^{u hbar d_y/2} z.
  - (15): xhat = Sum_i (1/S(alpha_i hbar d_y)) log(z - a_i)/alpha_i + R.
  - Remark 3.2: for y = z, w_i - wbar_j = z_i - z_j + (u_i+u_j) hbar/2, and the square roots give prod dz_i.
- **Conclusion on the transcription.** The author's (2.1)-(2.2) is an exact transcription, including:
  - the overall sign (-1)^n (-1)^(n-1) = -1;
  - the operator d(1/dx) without a minus sign;
  - 1/S (not S) in xhat.

  My re-derivation of (2.2) from (15) agrees, using d^(2k) log(v-a) = -(2k-1)!/(v-a)^(2k).

**Hypotheses for C_tau.** Take (x, y) = (X, v), with -y dx = -v dX = omega_{0,1} of P2.
- dX = (tau - v) dv/(v(1+v)) has a single simple zero, at v = tau (tau not in {0, -1}).
- dv has no zeros.
- Vital points of the dual: v = 0 (residue tau) and v = -1 (residue -(tau+1)).
- v = oo is a simple pole of dX with residue +1, but dv has a pole there, so it is not vital (Def. 2.4 and Remark 2.6).
- R = 0. Non-integer residues are allowed.
- The conventions agree: ABDKS's TR (blobbed TR plus projection, omega^(0)_1 = -y dx) is EO TR with omega_{0,1} = -v dX. The numerics fix the remaining sign (y -> -y would flip omega_{g,n} by (-1)^n, which the n = 1 and n = 3 comparisons with the Hodge tables exclude).

**Lemma 5.5.**
- It treats a simple pole of dx with dy of any order.
- Its last clause says: if the point is also a pole of dy, the right-hand side of (12) is holomorphic for all (g,n) != (0,1).
- The proof uses (92), with exponent 1 - s + 2(m+s-2) Sum r_i. For s = 1 and m = 2 this is 2 Sum r_i >= 2.
- This is exactly the situation at v = oo for C_tau: dX has a simple pole and dv a double pole. Lemma 5.6 would cover the same point with the roles of x and y swapped.

**Limitation of my reading.** I checked the statement and the structure of the proof of Lemma 5.5, not every line of it. My reading went through a summarising fetch tool, which refused verbatim quotation. However, formula (17) agrees with TR in the x-expansion at v = oo for n <= 5 (s.5 below). Any failure of the theorem at v = oo would show up exactly there.

**Version.** v3 (20 Jan 2025) has the comment "several corrections and clarifications" and postdates the IMRN 2024 version. Equation numbers (15)-(17) and (92), and Lemma 5.5, should be cited as those of v3 (issue 6).

## 2. Theorem 1 (the scaling argument): re-derived, no gap

- **ord.** ord is the order valuation on C(t_1..t_n): the lowest homogeneous degree, with ord(P/Q) = ord P - ord Q. So ord(fg) = ord f + ord g and ord(f+g) >= min(ord f, ord g).
- **The operator d(1/dX).** dX = [(1 - tau t)/(1+t)] dt/t. So d_i(1/dX_i) acts on the coefficient psi of dt_i/t_i as t_i d/dt_i (psi/q). This does not lower ord, because t_i d/dt_i preserves homogeneous degrees.
- **Orders of the ingredients.**
  - hbar^(2k) in E_i carries (v_i - a)^(-2k) = t_i^(2k)(1 - a t_i)^(-2k), which has ord 2k.
  - hbar^m in a cycle factor carries (v_i - v_j)^(-m-1), which has ord m+1, since (v_i - v_j)^(-1) = t_i t_j/(t_j - t_i).
  - For n = 1, 1/(u hbar) contributes ord 0 at hbar^(-1).
  - Each dv_i contributes ord -1.
- **Conclusion of the count.** Every term has ord >= 2g-2+n. omega_{g,n} = W prod dt_i/t_i with W rational and holomorphic at 0 (poles only at v_i = tau). So the Taylor components of W below degree 2g-2+n vanish. x_i = t_i (unit) is tangent to the identity, which gives [x^mu]F_{g,n} = 0 for |mu| <= 2g-3+n.
- **Uniformity.** The argument is uniform in g and n and uses no induction.
- **Consistency checks of the argument.**
  - Tau -> oo: the base curve C, X = -log v - 1/v, with R = -1/v and S(u hbar d)R - R of ord >= 2k+1. The same count gives the vanishing of the Hurwitz component, which is known independently (Route B, Thm B; Route C, Thm H).
  - At v = 0 the count gives nothing, because of the factors v^(-2k). This is consistent with the degree-1 behaviour there.

## 3. Theorem 3, the Norlund identities, Corollary 4: re-derived

- **Integration by parts.** Res x^(-j) d D^(r-1) psi = j^r Res x^(-j) eta. Then Sum_r j^r [u^(r+1)] E = (E(j) - 1)/j and x^(-j) E(j) = Psi_j. This gives (5.1).
- **Sinh form.** (j/alpha) S(j t)/S(alpha t) = sinh(jt/2)/sinh(alpha t/2).
- **Product form.** For integer N:
  - sinh(N u)/sinh(u) = Sum_{k<N} e^((2k-N+1)u), and e^(c d) log v = log(v + c);
  - Psi_j = prod (v+1+hbar s'_k)/prod (v + hbar s_k);
  - Res_{v=oo} f dv = -[w^1] f, which gives (5.3).
- **Expansion at oo.** The expansion of Psi_j at oo gives (5.2). The definitions of A_k, B_k and the sign of the a = -1 term are correct.
- **Norlund identities.**
  - B_N^(N+1)(x) = (x-1)...(x-N), from the residue N! Res (1+u)^(x-1) du/u^(N+1).
  - Appell property, and B_n^(a)(a/2) = n! [y^n] S1^a.
  - The inverse follows from Gamma(Z+a)/Gamma(Z+b) with a - b = -N. That ratio is a rational function, so its expansion is exact and not merely asymptotic.
  - Extension to all tau: both sides are Laurent polynomials in tau at fixed j.
  - All of this is correct.
- **Corollary 4.** I re-derived the term-by-term identity (-N)_{2g-2p} (M)_{2p} C(M-2p, j+1-2g)/C(M-1, j-1) = M N (j-1)!/(j+1-2g)!. The product equals M (M-1)_j, using M - j = N. MN/j^2 = sigma.
  - The h = 1 case is consistent: both p-terms equal MN.
  - The g = 2 hand formula agrees with every component of E_formulas\engine_h2_r1.txt (Xi^(0), Xi^(1), Xi^(2)).
- **What Corollary 4 does and does not add.** Without P2, the TR side equals (6.1) for all h. So (6.1) for all h is equivalent to P2 for r = 1 (mod ABDKS). It is a reformulation, not independent evidence for P2. Route E's h <= 7 verification of (6.1) therefore gives P2(r = 1) for h <= 7, mod ABDKS.

## 4. Proposition 5 and line (1)

- **Integer tau: correct.**
  - Q(x) is holomorphic at x(tau): by the quadratic loop equation for the branch pair and the linear loop equation for mixed pairs.
  - It is bounded at x = 0, using omega_{g,1}(v_oo)/dX = O(x) and v_oo ~ 1/x.
  - It vanishes at x = oo, where omega/dX -> 0 at v = -1. So Q = 0.
  - The global linear loop equation gives Sum_{i<j} Brk = -(1/2) Sum_i Brk(v_i, v_i). This needs g >= 2, so that omega_{0,2} is not inside Brk.
  - The omega_{0,1} terms give Sum_i v_i omega_{g,1}(v_i) dX. Hence Tr Phi_g = 0.
  - The equivalence with "the sigma-even polar part of Phi_g at v = tau vanishes" holds for integer tau.
- **Non-integer tau.** This needs one extra line: the polar coefficients are rational in tau and vanish at every integer tau >= 2. The author did not write this line (issue 3).
- **tau = 1.** The identity omega_{g,1} = Brk(v, 1/v)/((y(v) - y(1/v)) dX) is correct. The resulting one-line re-proof needs the (g-1, 2) bound, i.e. the n = 2 case in lower genus (issue 7).
- **The tau >= 2 "obstruction".** The reduction to the tau-section at v = 0 (in the coordinate s with s^tau = x), up to order (2g-3) tau, is correct. The statement that "induction gives no control there" is an explanatory remark, not a theorem (issue 2).

## 5. Independent computations (all exact mod p, p near 2^31; 1 or 2 primes)

Scripts are in TR2_recursion\referee. Every output cited is in referee\out\.

- **RefTR.pm: an independent EO recursion on C_tau.**
  - It works in the coordinate e = v - tau, with the global involution sigma(v) = tau + s(e). s(e) is solved order by order from X(tau+s) = X(tau+e) and verified.
  - Kernel: (1/2)[1/(v0 - v) - 1/(v0 - sigma v)] dv0/((y(v) - y(sigma v)) dX).
  - Correlators are stored in the pole basis prod dv_i/(v_i - tau)^(k_i).
  - [x^mu] is computed by Lagrange inversion: [x^mu] Int_oo dv/(v-tau)^k = -(1/mu)[t^(mu+1-k)](1 - tau t)^(-k)(1+t)^(mu(tau+1)).
  - It shares no code or coordinates with TRtau.pm or TRmod.pm.
  - It found a bug in my own first draft (the loop range), which the Hodge comparison caught. The fixed engine reproduces omega_{0,3} = -sigma, Xi_1, and Xi_2(3) = -73/20736 at tau = 2.
- **RefXY.pm: my own implementation of the integrated formula for general n.**
  - Formula: [x^mu]F = -(prod 1/mu_i) [hbar^(2g-2+n)] Res..Res prod Psi_{mu_i} Sum_{n-cycles} prod dv_i/(v_i - v_sigma(i) + (mu_i+mu_sigma(i)) hbar/2).
  - The residues are computed in C((t_1))...((t_n)); n = 1 is treated separately.
  - I checked that the integration by parts is legitimate in this iterated Laurent ring.
- **RefHodge.pm.** It reads Route E's tables E_formulas\engine_h<h>_r<r>.txt, which were computed without any TR, and evaluates Xi_h(R; tau) = Sum_k (-1/sigma)^k Xi^(k).

| check | tau | range | comparisons | mismatches |
|---|---|---|---|---|
| TR vs Hodge (P2) | 2 | g<=2, r<=3 | 51 | 0 |
| TR vs Hodge (P2) | 5/2 | g<=3 r<=4; g<=4 r<=3 (2 primes); g<=5 r<=2 (2 primes) | 121+157+150 | 0 |
| TR vs Hodge (P2) | -7/3 | g<=4 r<=3; g<=5 r<=2 (2 primes) | 157+150 | 0 |
| TR vs Hodge (P2) | 2/9 | g<=4, r<=2 (2 primes) | 96 | 0 |
| TR vs Hodge (P2) | 4 | g=3,4, r<=4 | 201 | 0 |
| TR vs Hodge (P2) | 1/2 | g<=4, r<=3 | 157 | 0 |
| TR vs Hodge (P2) | 3/7 | g<=4, r<=3 (2 primes) | 157 | 0 |
| formula vs TR / vs Hodge | 2 | g<=2, n<=3 | 36 / 36 | 0 |
| formula vs TR / vs Hodge | -7/3 | g=0..4, n<=3 | 158 / 157 | 0 |
| formula vs TR / vs Hodge | 4 | g=0..3, n<=4 | 127 / 121 | 0 |
| formula vs TR / vs Hodge | -2/5 | g=0..3, n<=5 | 112 / 75 | 0 |
| formula vs TR / vs Hodge | 2/9 | g<=5, n<=2 | 150 / 150 | 0 |
| formula vs TR / vs Hodge | 5/2 | g=6, n=1 | 15 / 15 | 0 |
| formula vs TR / vs Hodge | -7/3 | g=6, n<=2 | 62 / 62 | 0 |

- **Lemma D in TR.** 494 zeros in range across the TR runs, 0 violations. This includes tau = 3/7, g <= 4, n <= 3, which replaces the author's missing output. The formula runs add 246 more zeros.
- **tau = 1/2.** TR equals Hodge at every mu, including the vanishing all-even first-layer entries. So that curiosity holds on the Hodge side too (sigma = 3/4).
- **New range for the general-n integrated formula.** It is now verified for n = 3 (g <= 4), n = 4 (g <= 3) and n = 5 (g <= 3). The author had checked only n <= 2.
- **Author's values.** spot_values.pl confirms, with both my TR and Hodge and modulo 2 primes:
  - Xi_5(9) at tau = 3, 3/7, -5/3;
  - Xi_4(7) and Xi_3(5) at tau = 2;
  - Xi_2(3) at tau = 3/7.

  I also re-ran the author's out\tau1_values.pl, whose output had not been stored. It gives -2205/1048576 (g = 5) and -297675/33554432 (g = 6), the referee values in referee_hurwitz_1.md s.2(d) and proof_status 1.8.
- **Proposition 5, local form.** qle_check.pl finds that the polar part of Phi_g(v) + Phi_g(sigma v) vanishes at:
  - tau = 3 (g <= 4);
  - tau = 5/2 (g <= 4; non-integer);
  - tau = -7/3 (g <= 3).

  Phi_g alone has 10 to 18 polar coefficients. A negative control (1 in place of 1/2 in front of Brk) fails, as it should.

**Normalisation items in the brief's checklist.**
- **Not used by TR2.** BKMP/Chen/Zhou framing and signs, the EO variational formula and DOSS enter only through P2. TR2 treats P2 as a black box.
- **Tested directly against exact data.** The Hodge comparisons above test P2's own normalisation:
  - c = 1;
  - sigma^(h-1) for every r;
  - the leaf binomial C(mu(tau+1)-1, mu-1) = Pi_mu(tau);
  - the sign at odd n;
  - repeated parts, with no 1/|Aut| factor needed.

## 6. Rules

- TR and the x-y swap are applied only to the Hodge base curve C_tau. Nothing is applied to the resolved conifold x C^2: grep finds no conifold, Q or polylog use in the TR2 scripts, and (H) remains Brini-Schuler.
- KP is not used.
- No wolframscript was run. tr2_closed_check.wls is written but not run, and concerns only C_tau.
- I ran no wolframscript either.

## 7. Issues

1. **(minor) Status labels for claims 3 and 5 (Corollaries 2 and 4).** They depend on P2 (Lemma B). P2 is not literature and has not been written out (route TR1). The correct status is "conditional on P2". The report text already says "conditional only on P2 and [ABDKS]". Lemma A, (OL) for all G, Conjecture M (vanishing) and Corollary R in all genera must be recorded as **mod P2**.
2. **(minor) Claim 9.** The reduction (Lemma D(g,1) is equivalent to the vanishing of the tau-section of Phi_g at v = 0 up to order (2g-3) tau) is proved. "Not controlled by induction / the precise obstruction" is heuristic. Split the status.
3. **(minor) Claim 8, "for all tau in its equivalent local form".** Add the one-line rationality argument: the polar coefficients at v = tau are rational in tau and vanish at all integers tau >= 2.
4. **(cosmetic) Missing outputs.**
   - There is no file for the tau = 3/7 Lemma D run (82 coefficients cited).
   - The output of tau1_values.pl was not stored.
   - cmp_n1_tau1_g6.txt is empty; this was acknowledged.
   - The summary says "21 points", claim 1 says "22".

   My runs cover the first three.
5. **(cosmetic) Claim 10's label.** For h = 6 and 7 the TR side comes from Theorem 3 (mod ABDKS), not from a TR run. Label it "mod ABDKS, plus Route E's computer verification". My direct runs confirm P2 at h = 6, r <= 2 without ABDKS.
6. **(cosmetic) Citation.** Cite arXiv v3 for the numbering of (15)-(17), Lemma 5.5 and (92). v3 lists "several corrections and clarifications" relative to earlier versions.
7. **(cosmetic) The tau = 1 re-proof of Theorem 3.1 (n = 1).** State that it uses the (g-1, 2) bound.
8. **(upgrade, not an error) The general-n integrated formula** is now verified for n <= 5 (s.5).

No claim is refuted.

## 8. Files (TR2_recursion\referee)

- **Modules:** RefMod.pm (mod-p series), RefTR.pm (independent EO), RefXY.pm (integrated formula (17), general n), RefHodge.pm (Route E tables).
- **Scripts:**
  - cmp_hodge.pl (TR vs Hodge through P2, plus the Lemma D check): `perl cmp_hodge.pl tn td G N [X] [NP] [GMIN]`
  - cmp_xy.pl (formula vs TR vs Hodge): `perl cmp_xy.pl tn td G N [X] [NP] [GMIN]`
  - qle_check.pl: `perl qle_check.pl tn td G`
  - spot_values.pl, t_ref.pl
- **Outputs:** out\hodge_*.txt, out\xy_*.txt, out\qle_check.txt, out\spot_values.txt, out\author_tau1_values_rerun.txt.

