# Lemma A round 2: referee2 TR2 (saved by the main session from the workflow journal, 2026-10-31)

Subagents cannot write report files, so this file holds the agent's structured return value: the fields first, then the full report text (report_text) verbatim.

## confirmed

- Claim 1: [ABDKS] Prop. 3.1, formula (17), applies to C_tau for every tau not in {0,-1}. Hypotheses checked against the arXiv v2 and v3 text: standing assumptions of s.2.3; Def. 2.4 (no vital point on (X,v); v=0 and v=-1 vital on the dual; v=oo not vital); Rem. 2.6; Rem. 2.12 (LogTR = TR); Lemma 5.5's last clause covers v=oo, where dX has a simple pole and dy a double pole. The transcription (2.1)-(2.2) agrees with (17), (15), (16), (20) of v3, sign included.
- Claim 2 (Theorem 1, Lemma D for all tau, g, n): scaling argument checked line by line. Direct exact-TR test with an independent engine: all (g,n) with g<=8, n<=9-g at tau=3, 5/2, -3/4, 1/3 and generic tau=123456789/987654321; g<=9, n<=10-g at tau=7, -7/3. 3259 in-range coefficients, each zero mod two primes; 0 violations.
- Claim 3 (Corollary 2): the logic is correct given P2. P2 is used only at integer tau>=1, then polynomiality in 1/sigma. (OL), the vanishing direction of Conjecture M (via the string/dilaton reduction 1.2) and Corollary R (via Lemma C, proved) follow. No circularity, no induction.
- P2 (assumed) agrees numerically with Route E's exact Hodge tables, multi-leaf, with independent TR: (h,r) up to (6,3), (5,4), (4,5) at tau=2/5 (523 values) and generic tau (499 values), including two layers beyond the range; 0 mismatches.
- Claim 4 (Theorem 3, (5.1)-(5.3)): re-derived. (5.2) matches TR for g<=8 (tau=3), g<=7 (generic tau), g<=6 (tau=-7/3), g<=4 (tau=2).
- Claim 5 (Corollary 4): the termwise identity was re-derived by hand. (6.1) matches TR in the same ranges.
- Claim 6 (Norlund identities): re-derived (Norlund polynomial, Appell property, DLMF 5.11.17, extension in tau). 1488 numerical checks against exp-coefficients computed from first principles.
- Claim 7 (two-point formula (5.4)): matches TR in both residue orders: g<=8 (tau=3, 288 coefficients), g<=7 (generic, 252), g<=6 (tau=-7/3, 126), g<=4 (tau=2, 60).
- General-n integrated formula (author: checked only for n<=2): n=3 to g=5 (tau=-7/3) and g<=4 (generic tau); n=4 to g=2 (generic tau); 0 mismatches.
- Claim 8 (Proposition 5): proof re-derived (behaviour of Q at x(tau), 0, oo; global linear loop equation).
- Claim 9: the tau=1 identity x F'_{g,1} = t Brk_g(v,v)/((1-t^2) dX^2) holds exactly for g=2..7 through t^34.
- Claim 11: consistent; the extra first-layer zeros at tau=-3/4 are explained by the P2 normalisation zeros.

## issues

### [1]

#### description

Labelled 'proved-modulo-literature', but they depend on P2 (Lemma B), which is not a published result (sketched route via BKMP for C^3 and the EO variational formula, being written out by TR1). The statements themselves say 'Assuming P2'.

#### location

Structured claims 3 (Corollary 2) and 5 (Corollary 4); claim 10

#### severity

minor

#### suggested_fix

Relabel as 'conditional on Lemma B (P2); otherwise proved modulo [ABDKS]'. Claim 10 holds mod [ABDKS].

### [2]

#### description

The equation and lemma numbering used is that of arXiv:2312.16950v3 (Jan 2025; arXiv comment: 'several corrections and clarifications'). The explicit formula is (17) in v3 but (15) in v2. The IMRN 2024 numbering was not verified. Prop. 3.1 and Lemma 5.5, including the common-pole clause, read the same in v2 and v3.

#### location

Report s.2.1 and all citations of [ABDKS]

#### severity

minor

#### suggested_fix

Cite arXiv v3 explicitly with its numbering (Thm 2.8, Prop. 3.1, (17), Lemma 5.5), note the version, or check the IMRN numbering.

### [3]

#### description

'For integer j tau' should be 'for N = j tau a positive integer'. The product formula (5.3) needs N >= 1; the report body says this correctly.

#### location

Structured claim 4 (Theorem 3)

#### severity

cosmetic

#### suggested_fix

Replace with 'when N = j tau is a positive integer'.

### [4]

#### description

'The first layer is non-zero' holds only at the tested tau. At tau = -3/4 the first layer vanishes exactly on the mu with a part divisible by 4, where the P2 normalisation C(mu_i(tau+1)-1, mu_i-1) vanishes.

#### location

Structured claim 11

#### severity

cosmetic

#### suggested_fix

Say 'generically non-zero (zero where the P2 normalisation vanishes, e.g. tau = -3/4)'.

### [5]

#### description

For non-integer tau there is no global fibre sum, so the 'equivalent local form' is really a local version. It follows from the linear and quadratic loop equations alone.

#### location

Structured claim 8 (Proposition 5) and report s.6(1)

#### severity

cosmetic

#### suggested_fix

Write 'local version (sigma-even polar part of Phi_g at v = tau vanishes), valid for all tau'.

### [6]

#### description

The pole cancellation at the vital point v = 0 is Lemma 5.2 of [ABDKS] applied to the dual pair (v, X), whose y is X. The remark does not say so.

#### location

Report s.3, Remark (b)

#### severity

cosmetic

#### suggested_fix

Add 'applied to the dual pair'.


## notes

No failure found. The all-genus claim stands as stated. Theorem 1 holds modulo [ABDKS] (published; hypotheses verified for C_tau, Lemma 5.5 covering v = oo). Lemma A, (OL), Conjecture M (vanishing) and Corollary R in all genera hold given P2 (Lemma B), which TR1 must still write out; P2 now also has independent multi-leaf numerical support through genus 6. All issues are labels or wording. Independent engine and scripts: G:\My Drive\Projects\Old\refGW\hodge_vertex_conifold\lemmaA\round2\TR2_recursion\referee2\ (RefTR.pm, RefXYc.pm, lemmaD_ref.pl, closed_ref.pl, closed_n34.pl, p2_engine_ref.pl, tau1_identity.pl, rerun_ref2.sh; outputs in out\).

## refuted


## route

TR2 (round 2): x-y swap for log TR on C_tau, giving Lemma D for every tau, g, n, and Lemma A, (OL) and (6.1) in all genera (mod P2). Second, independent referee.

## verdict

sound

## report_text

# Referee 2 report: route TR2, round 2 (x-y swap route to Lemma D and Lemma A)

Referee folder: G:\My Drive\Projects\Old\refGW\hodge_vertex_conifold\lemmaA\round2\TR2_recursion\referee2\ (scripts, plus outputs in out\). No wolframscript was run. No TR was applied to the resolved conifold x C^2. KP was not used.

## Verdict: sound. I found no failure.

I set out to break the argument and could not.
- **Theorem 1 holds.** Lemma D at every tau not in {0,-1} and every (g,n) follows correctly from one published theorem: [ABDKS] = Alexandrov-Bychkov-Dunin-Barkowski-Kazarian-Shadrin, arXiv:2312.16950 (IMRN 2024), Thm 2.8 and Prop. 3.1. I checked its hypotheses for C_tau against the text of arXiv v2 and v3. The step from the formula to the degree bound, a termwise scaling count, is correct as written.
- **Corollary 2 holds, given P2.** Lemma A for all h, then (OL), the vanishing direction of Conjecture M, and Corollary R follow from Theorem 1 and P2 (Lemma B). P2 is the only unproved input and is used exactly where the author says. P2 is not a literature result: it is still to be written out (route TR1). So "conditional on P2" is the correct status, not "proved modulo literature" (issue 1).
- **The closed formulas are right.** The one-point formulas (5.1)-(5.3), the identity behind (6.1), the two-point formula (5.4), the Norlund identities and Proposition 5 are correct.
- **The general-n formula, which the author checked only for n <= 2, also holds for n = 3 and 4.** I checked it against exact TR.

Remaining issues are labels, citations and wording only (s.4).

## 1. Steps asserted for all genera, and tests beyond the author's range

I wrote a new Eynard-Orantin engine, RefTR.pm, from scratch, sharing no code with TRmod.pm or TRtau.pm. It works mod p.
- **Setup.** The local coordinate z is defined by X = X(tau) - z^2/(2 sigma). The basis is xi_k = (1/(2k+1)) Res B z^(-2k-1). The kernel is reduced to omega_{g,n+1} = Sum_k0 xi_k0 (sigma/2) [z^(-2k0)] (r/Y), with Y = (u(z)-u(-z))/(2z). B-terms are kept only through their even part; the odd part provably drops out. The expansion at v = oo is done directly in t = 1/v.
- **Validation (t_engine.pl, out\t_engine.txt).**
  - omega_{0,3} = -sigma, and omega_{1,1} = (163/1728) xi_0 - (3/4) xi_1 at tau = 2.
  - The author's tau = 2 values Xi_2(3), Xi_3(5), Xi_4(7).
  - The Hodge formulas for Xi_1(j) and Xi_2(j) at tau = 2, 3/7, -5/3, 5/2.
  - The earlier referee's tau = 1 values -2205/1048576 (g = 5) and -297675/33554432 (g = 6, zeros for j <= 10).
- **Speed.** The engine computes only the needed (g',n'), those with n' <= n + g - g'. Genus 8 takes about 20 s per prime; genus 9 about 3 min.

"Generic tau" below means tau = 123456789/987654321. Every coefficient is a rational function of tau over Q, so agreement or vanishing at this point, mod two primes, certifies it identically in tau with negligible error probability.

| Step asserted for all genera | Author's range | Referee test (independent code), result |
|---|---|---|
| (2.1) = TR(C_tau), all (g,n) (ABDKS) | n=1 exact to g=5; n=2,3,4 at 22 points (n=2 to g=5, n=3 to g=3, n=4 to g=1) | Via its integrated consequences, coefficientwise at v = oo: n = 1 to g = 8; n = 2 to g = 8; n = 3 to g = 5; n = 4 to g = 2. 0 mismatches (rows below). |
| Theorem 1 (Lemma D), all tau, g, n | TR to g=5 (n<=2) and g=4 (n<=3), at 3 values of tau | Direct TR, all (g,n) with g <= 8 and n <= 9-g at tau = 3, 5/2, -3/4, 1/3 and generic tau (381 in-range coefficients each); g <= 9 and n <= 10-g at tau = 7, -7/3 (677 each). Total 3259 coefficients, each zero mod two primes, 0 violations. Largest n per genus: (2,8), (3,7), (4,6), (5,5), (6,4), (7,3), (8,2), (9,1). |
| Theorem 3, formula (5.2) | g<=5 (tau=2), g<=4 (3/7, -5/3) | Against TR: g <= 8 at tau = 3 (104 values); g <= 7 at generic tau (84); g <= 6 at tau = -7/3 (54); g <= 4 at tau = 2 (28). 0 mismatches. |
| Norlund identities (claim 6) | 224 checks | Exp-coefficients from first principles against the closed forms: 540 + 432 + 336 + 180 checks, 0 mismatches. Hand check at N = 1, 2. |
| Corollary 4: (5.2) = (6.1) termwise | algebraic check to g=9 | Identity re-derived by hand. (6.1) against TR in the same ranges as (5.2), 0 mismatches. |
| Two-point formula (5.4) | g<=4, 148 coefficients | Against TR, both residue orders (which tests the symmetry): g <= 8 at tau = 3 (288); g <= 7 generic (252); g <= 6 at tau = -7/3 (126); g <= 4 at tau = 2 (60). 0 mismatches. |
| General-n integrated formula (report s.5) | n <= 2 only | n = 2 matches (5.4). n = 3: g <= 5 at tau = -7/3 (120 values, mu_i <= 4); g <= 4 generic (100, of which 78 non-zero); g <= 3 at tau = 2 (40). n = 4: g <= 2 generic (45, of which 43 non-zero); g <= 1 at tau = 2 (10). 0 mismatches. |
| Proposition 5, claim 9 (tau = 1 identity) | QLE local form to g=4 | Proof re-derived (s.2.6). The identity x F'_{g,1} = t Brk_g(v,v)/((1-t^2) dX^2) holds exactly through t^34 for g = 2..7. |
| P2, used by Corollaries 2 and 4 | r = 1, h <= 7 (via (6.1) and Route E) | Independent TR against Route E's exact Hodge tables (engine_h*_r*.txt). All mu with |mu| <= 2h-2+r+2, which is the range plus two further layers. (h,r) = (2,1), (2,2), (3,3), (4,2), (4,4), (4,5), (5,2), (5,3), (5,4), (6,2), (6,3) at tau = 2/5 (523 values). (2,3), (3,4), (4,4), (4,5), (5,2), (5,3), (5,4), (6,2), (6,3) at generic tau (499 values). 0 mismatches. |

Side observation. At tau = -3/4 some first-layer coefficients vanish: exactly those mu with a part divisible by 4. These are the zeros of the P2 normalisation C(mu_i(tau+1)-1, mu_i-1), so they are consistent with P2 and are not a problem.

I did not use Route A's KP code. All claims under test are TR statements or P2, and the Hodge side was covered by Route E's tables.

## 2. Logic

### 2.1 Applying [ABDKS] to C_tau, hypothesis by hypothesis

Their convention is omega^(0)_1 = -y dx. Take x = X and y = v; then -v dX is our omega_{0,1}, and B agrees.

- **Standing assumptions (s.2.2-2.3):**
  - Sigma = P^1 is compact.
  - dX = (tau-v) dv/(v(1+v)) and dv are meromorphic for every complex tau, integer or not.
  - The zeros of dX are {tau}: simple, since tau is not 0 or -1. dv has no zeros, so the zero sets are disjoint.
  - dv is regular and non-zero at v = tau.
- **Def. 2.4 (vital points).** A point is LogTR-vital if it is a simple pole of dy with non-zero residue where dx has no pole.
  - On (X, v), dv has no simple pole, so there is no vital point. By Rem. 2.12, LogTR(X,v) = TR(X,v).
  - On the dual (v, X), v = 0 and v = -1 are vital: dv is regular there, with residues tau and -(tau+1).
  - v = oo is not vital, because dv has a pole there. By Rem. 2.6 its principal part vanishes for g > 0.
- **Prop. 3.1.** Its hypotheses: y = z is an affine coordinate with z = oo a pole of dy, and x = Sum alpha_i^(-1) log(z - a_i) + R.
  - Here R = 0. The log singularity of X at oo (residue +1) is produced automatically by the two logs.
  - Its S-deformation terms d^(2k) log(v-a) ~ v^(-2k) vanish at oo, consistent with oo being non-vital.
- **Proof structure of Thm 2.8 (Lemmas 5.1-5.6), point by point on P^1:**
  - v = tau (zero of dX) and the diagonals: Lemma 5.1.
  - v = 0, -1 (simple poles of dX, dy regular): Lemma 5.5, first two clauses. They give the vital behaviour on the dual side.
  - v = oo (simple pole of dX, double pole of dy): Lemma 5.5, last clause: "if z = 0 is also a pole of dy, the right side is holomorphic for all (g,n) != (0,1)".
  - Lemmas 5.2-5.4 are not needed: no pole of dy is free of a pole of dx. Lemma 5.6 is not needed either: dX has no higher-order poles.
  - I read Lemma 5.5's statement in both v2 and v3, and they agree. The point the author asked a referee to re-read is covered.
  - The dual LogTR has no ramification points, so it consists of the vital terms only. Uniqueness then gives Prop. 3.1.

### 2.2 Transcription

The author's (2.1)-(2.2) agree with [ABDKS] v3 (17), (15), (16), (20):
- The operator is (d_i (1/dx_i))^r, with no sign.
- E_i = exp(-u_i (S(u_i hbar d) xhat_i - x_i)).
- The signs are (-1)^n (-1)^(n-1) = -1.
- For y = z: w_i - wbar_j = z_i - z_j + (u_i + u_j) hbar/2.
- The expansion (2.2) follows from d^(2k) log(v-a) = -(2k-1)!/(v-a)^(2k).

The overall sign of (2.1) is correct as printed. (5.1), (5.4) and the n = 3, 4 formulas, all derived from it with its leading minus sign, match my TR coefficientwise. The internal sign flip in XY.pm is an internal convention only.

### 2.3 Theorem 1 (scaling argument), line by line

All steps are correct.
1. For 2g-2+n > 0, omega_{g,n} has poles only at v = tau, so W is holomorphic at t = 0 in the frame prod dt_i/t_i.
2. In this frame, d(1/dX) acts by psi -> t d/dt (psi/q(t)) with q = (1 - tau t)/(1+t) and q(0) = 1. So it does not lower ord. The reason: (t_i d_i f)(lambda theta) = theta_i d/dtheta_i [f(lambda theta)], applied termwise to the lambda-expansion.
3. [hbar^(2m)] E_i has ord >= 2m in every u-degree.
4. [hbar^m] of the cycle product has ord exactly m+n for generic theta with distinct entries. For n = 1 the factor 1/(u hbar) gives ord 0 at hbar^(-1).
5. Each dv_i contributes ord -1.
6. At fixed (g,n) the u-degrees are bounded, so the sums are finite.
7. Hence ord W >= 2g-2+n. All homogeneous Taylor components of degree below 2g-2+n vanish. The change of variables t -> x is tangent to the identity, so [x^mu] F_{g,n} = 0 for |mu| <= 2g-3+n.

There is no induction, so no base case is involved.

The argument is specific in the right way:
- It needs R = 0. A polynomial part of x would put non-decaying derivatives into E_i.
- It works only at v = oo. At the vital points the factors (v-a)^(-2k) blow up, consistent with the degree-1 behaviour at v = 0.
- So it does not over-reach to the counterexample classes.

### 2.4 Corollary 2 and the final statements

- **Corollary 2.** P2 gives [x^R] F_{h,r} = sigma^(1-h) prod C(R_i(tau+1)-1, R_i-1) Xi_h(R).
  - It is needed only at integer tau >= 1, where the binomials are positive.
  - Xi_h(R) = Sum_{k<=h} (-1/sigma)^k Xi^(k)_h(R) with Xi^(k) independent of tau (status file 1.1(e)). So vanishing at infinitely many sigma means vanishing identically.
  - Parts >= 2 with D(R) <= 2h-3 are exactly the case |R| <= 2h-3+r. This is Lemma A for every h and every r <= 2h-3.
  - Theorem 1 is in fact stronger, since it also covers parts equal to 1. This is consistent with the string relation Xi(1,R) = J_R Xi(R).
- **(OL).** R = {D+1} gives (OL) for every G and D, with no inductive hypothesis. This is more than brief s.3 asks for: brief s.3 says any direct all-genus argument implies Lemma A.
- **Conjecture M.** The vanishing direction follows from M_irr through the proved string/dilaton reduction (status 1.2).
- **Corollary R.** It is conditional on Lemma A only, since Lemma C is proved (status s.2, Addendum 2), together with input (H) of Brini-Schuler, which is a theorem.
- **No circularity.** Theorem 1 is a pure TR statement on C_tau. P2's proof route (BKMP for C^3, EO variational formula) does not use Lemma A.

### 2.5 Theorem 3 and Corollary 4

I re-derived each step.
- **Integration by parts.** The residues of exact forms at oo vanish, because x^(-j) = t^(-j)(1+t)^M is single-valued at t = 0 even for non-integer tau. Then Sum_r j^r [u^(r+1)] E = (E(j) - 1)/j, and x^(-j) E(u = j) = Psi_j.
- **Expansion at oo.** It gives (5.2).
- **Norlund identities.**
  - B_N^(N+1)(x) = (x-1)...(x-N), via N! Res (1+u)^(x-1) du/u^(N+1).
  - B_n^(a)(a/2) = n! [y^n] S1^a.
  - The inverse via DLMF 5.11.17 (Tricomi-Erdelyi).
  - Extension to all tau: both sides are Laurent polynomials in tau and agree at the points j tau in Z_{>0}.
- **Termwise identity.** (-N)_{2g-2p} (M)_{2p} C(M-2p, j+1-2g) / C(M-1, j-1) = M N (j-1)!/(j+1-2g)! holds, using (M)_{2p}(M-2p)_{j+1-2g} = (M)_{j+1-2g+2p} and M - j = N. With MN/j^2 = sigma this gives (6.1) for all h (mod P2).
- **Hand check.** h = 1 gives ((j-2) sigma - 1)/24.

### 2.6 Proposition 5 and claim 9

- **Proposition 5.** The rational function Q(x) is holomorphic at x(tau): by the quadratic loop equation for the ramified pair and the linear one for mixed pairs.
  - At x = 0 and x = oo, Q = O(|x|^(-1+1/tau)) and O(|x|^(-1/(tau+1))) respectively. A rational function with these bounds is holomorphic at 0 and vanishes at oo, so Q = 0.
  - The global linear loop equation turns this into -Tr Phi_g.
  - The local version for all tau follows directly from the loop equations.
- **Claim 9 at tau = 1.** It follows from anti-invariance, and it is the exact identity omega_{g,1} = Brk(v,1/v)/((y(v) - y(1/v)) dX).
- **The obstruction for tau >= 2.** The description is accurate. It is an equivalence given the lower-genus bound, plus an observation.

### 2.7 Rules

- TR is used only on C_tau.
- The x-y swap is a TR result for the Hodge CohFT, so it is allowed.
- KP, Corollary E, Theorem D5 and the parity argument are not used.
- tr2_closed_check.wls is only requested, not run.
- No report file is in the author's folder.

## 3. Literature

- **[ABDKS], arXiv:2312.16950.** Read from the arXiv HTML of v2 (Jan 2024) and v3 (Jan 2025). v3's arXiv comment reads "several corrections and clarifications".
  - Def. 2.4, Rem. 2.6, Rem. 2.12, the standing assumptions of s.2.3, Prop. 3.1 and Lemma 5.5 are as the author states.
  - The explicit formula is (17) in v3 but (15) in v2, with xhat = (13) and w = (14) there. I did not check the IMRN numbering.
  - The proof of Thm 2.8 relies on [BDBKS21] (divisibility lemmas) and on the universal x-y swap (JEMS 2025). I did not re-prove it; the author does not claim to have either.
- **DLMF 5.11.17** (Gamma(z+a)/Gamma(z+b) ~ z^(a-b) Sum C(a-b,k) B_k^(a-b+1)(a) z^(-k)) and **Norlund's** B_n^(n+1)(x) = (x-1)...(x-n): standard; I re-derived both above.
- **P2 (Lemma B)** is not literature. Status file 1.6(b): it is a theorem modulo BKMP for C^3 (Chen; Zhou) and the EO variational formula, along a route that is sketched but not written out.

## 4. Issues (none affects correctness)

1. **(minor) Status labels.** Claims 3 (Corollary 2) and 5 (Corollary 4) are labelled "proved-modulo-literature", but they rest on P2, which is not a published result. The label should read "conditional on Lemma B (P2); otherwise proved modulo [ABDKS]". The statements themselves already say "Assuming P2". In the same spirit, claim 10 holds mod [ABDKS].
2. **(minor) Citation and version.** Cite arXiv:2312.16950v3 explicitly, with its numbering: Thm 2.8, Prop. 3.1, formula (17), Lemma 5.5. Also note that v3 corrects the earlier versions; for example the formula is (15) in v2. Either check the IMRN numbering or cite by arXiv version.
3. **(cosmetic) Claim 4.** "For integer j tau" should read "for N = j tau a positive integer". The report body says this correctly.
4. **(cosmetic) Claim 11.** "Non-zero first layer" holds only at the tested tau. At tau = -3/4 the first layer vanishes exactly on the mu with a part divisible by 4, the zeros of the P2 normalisation. Say "generically non-zero".
5. **(cosmetic) Claim 8.** For non-integer tau there is no global fibre sum, so "equivalent local form" should be "local version". That version follows from the loop equations alone.
6. **(cosmetic) Remark (b) in s.3.** The pole cancellation at v = 0 is Lemma 5.2 applied to the dual pair (v, X), whose y is X. The text could say so.

## 5. Files

All in G:\My Drive\Projects\Old\refGW\hodge_vertex_conifold\lemmaA\round2\TR2_recursion\referee2\:
- **Engine and formulas:**
  - RefTR.pm: independent EO engine.
  - RefXYc.pm: closed formulas implemented from first principles: Psi_j, (5.2), the Norlund forms, (5.4) and the general-n formula.
- **Check scripts:**
  - t_engine.pl
  - lemmaD_ref.pl
  - closed_ref.pl
  - closed_n34.pl
  - p2_engine_ref.pl
  - tau1_identity.pl
  - rerun_ref2.sh (all commands).
- **Outputs in out\:**
  - t_engine.txt
  - lemmaD_tau{3_1, 5_2, -3_4, 1_3}_G8.txt, lemmaD_tauRandom_G8.txt
  - lemmaD_tau{7_1, -7_3}_G9.txt
  - closed_tau3_1_G8.txt, closed_tau-7_3_G6_n3.txt, closed_tauRandom.txt
  - p2_engine_tau2_5.txt, p2_engine_tauRandom.txt
  - tau1_identity.txt
