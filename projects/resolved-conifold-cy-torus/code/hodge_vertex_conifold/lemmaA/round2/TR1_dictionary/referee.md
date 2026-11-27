# Lemma A round 2: referee TR1 (saved by the main session from the workflow journal, 2026-10-31)

Subagents cannot write report files, so this file holds the agent's structured return value: the fields first, then the full report text (report_text) verbatim.

## confirmed

- Theorem P2 (Lemma B): omega^{C_tau}_{g,n} = sigma^(1-g) Sum_k X_{g,n}[k] prod dxi_{k_i} for every tau not in {0,-1} and every stable (g,n), i.e. [x^mu]F = sigma^(1-g) prod Pi_{mu_i} Xi_g(mu), with the brief's normalisation; proved modulo BKMP for framed C^3 and the EO symmetry theorem (status proved-modulo-literature justified, given the base-case remark).
- BKMP in the author's conventions, omega^{log(1+w)} = sigma^(g-1+n) Sum M[k] prod dxi_k (Zhou's y = -log(1+w) gives (-1)^n): Zhou's definitions (3)-(7), curve (8)-(11) and recursion (19)-(23) checked against the paper, including |Aut| factors, kernel orientation, the B vs B - S issue and framing a = tau. Holds for 2g-2+n >= 2 from Zhou; (0,3), (1,1) from s.5. Independently verified at tau = 4 (g <= 3) and Zhou-literally at tau = 5/2 (g <= 2, n <= 3).
- Lemma 1 (loop equations), Lemma 2 (variational formula, including the implicit bijection of E-terms, the identity (**) and all order estimates), Lemma 3 (dual cycle; halves irregular, sum -1/w regular), Lemma 4 (rationality in e, Yh(0) = 2(tau/(1+tau))(1+e sigma), e -> oo limit), Lemma 5 (dilaton resummation): all re-derived, correct.
- Steps 1-2 of s.4: Taylor coefficients at e = 0, identification of rational functions, the limit giving the n-independent factor sigma^(1-g), extension from Z>=1 to all tau via Q[tau,1/tau,1/(1+tau)].
- Base cases re-derived by hand for all tau and e: omega^{Y0}_{0,3} = sigma dxi_0^3, omega^{yB}_{0,3} = sigma^2 dxi_0^3, omega^{Y0}_{1,1} = -((2+1/sigma)/24)dxi_0 + (1/24)dxi_1, omega^{yB}_{1,1} = -((1+sigma)/24)dxi_0 + (sigma/24)dxi_1, and Theorem F at (1,1); numerically at 24 tau (384 coefficients, 0 mismatches).
- Unstable cases: [x^mu] Int Y_0 dX = sigma Pi_mu/(mu(mu+1)) (Lagrange inversion) and [x1^m1 x2^m2] log((w1-w2)/(x1-x2)) = sigma Pi Pi/(m1+m2) (Zhou (67), (74)).
- Theorem F (interpolating family), confirmed independently at (tau,e) = (3/5,-2) and (4,-1/21), g = 2,3, n <= 2; its x-form equals Route A's Cor. 2.2 with u = -e sigma (for g >= 2 or g = 1, n >= 1).
- Moebius classification: a non-identity s with s^*dX = dX exists iff tau is in {1,-2,-1/2} (residues 1, -(1+tau), tau of dX).
- P2 independently verified (referee code TRw.pm, unisolvent lattices, Hodge data from Mathematica q/mv tables and the Route E engines) at tau = 1, 4, 3/5, -7/3, -2/7, 5/2. At tau = 4: all (g,n) with g <= 5, n <= 3, g+n <= 6, including (3,3), (4,2), (5,1).
- P2 directly verified at genus 7 (n = 1) at seven tau: the author's saved TR coefficients equal the Route E referee engine values (539 coefficients, 0 mismatches); the TR-free tables reproduce the referee_hurwitz_1 tau = 1 values and -15436575/268435456 at h = 7, j = 13.
- Lemma A at sigma = 2 (tau = 1, -2) in every genus, via Theorem 3.1 (TR statement re-checked) plus P2; Lemma A at sigma = -1/4 for odd parts via the tau = -1/2 variant plus P2 (already Route A's Theorem C).
- proof_status 1.7(b),(c) and (6.1) at sigma = 2, -1/4 (h <= 6) are now unconditional modulo P2's literature inputs; (6.1) holds identically for h <= 7 (the sigma-degree argument is correct).
- TR is applied only to the Hodge CohFT (C_tau, framed C^3 vertex), never to the resolved conifold x C^2.
- Author's unstored checks reproduced with the author's scripts: Zhou's (6) literally at tau = 2/3 (139/139), genus 0 up to n = 8 (1728/1728), tau = -5/7 g <= 1 (29/29).

## issues

### [1]

#### description

Zhou's Theorem 1 (arXiv:0910.4320v2) does not cover the base cases (0,3) and (1,1). At (1,1) his recursion (21) contains W_0(q,q~) = B(q,q~) - dx(q)dx(q~)/(x(q)-x(q~))^2, whose second term is infinite since x(q~) = x(q); his derivation of (65) from the cut-and-join (68) uses the stable expansions (77)/(79) of Phi_{g-1,n+1} and Phi_{g,n-1}, valid only for 2g-2+n >= 2; his s.3.3 checks only (19)-(20). The proof of P2 is nevertheless complete because report s.5 computes (0,3) and (1,1) directly for all tau (re-derived by the referee: omega_{0,3} = e_0^3/(2 c2 y1), omega_{1,1} = (1/(4 c2 y1))[(beta0 - y3/(4y1)) e_0 + e_1/12]).

#### location

report s.2.1 (transcription of Zhou) and s.4 'Uniformity in the genus' ('Zhou's theorem is an all-(g,n) statement')

#### severity

minor

#### suggested_fix

State BKMP as: (0,3), (1,1) by the direct computation of s.5; for 2g-2+n >= 2 by Zhou's recursion (which then never involves W_0(q,q~)) and induction. Note that the (1,1) term must be read with B(q,q~).

### [2]

#### description

(6.1) at h = 7 is not new: Route E's referee verified it identically without TR (E_formulas/referee.md; brief s.2(g)). The genus-7 Hodge tables E_formulas/referee/ref_engine_h7_r1..3.txt exist; comparing the author's saved one-leaf TR coefficients with them (7 tau, g = 1..7, 539 coefficients, 140 at genus 7) gives 0 mismatches, so P2 is directly verified at genus 7.

#### location

summary ('New result ... (new at h = 7)'), claims list ((6.1) for h <= 7), s.6 table ('(7,1): no Hodge table'), s.7(e)

#### severity

minor

#### suggested_fix

Describe the genus-7 TR computation as an independent confirmation, add the direct genus-7 P2 comparison with ref_engine_h7_r1.txt, and remove '(no Hodge table)'.

### [3]

#### description

Zhou (0910.4320) and Chen (0910.3739) have no journal reference on arXiv (unrefereed preprints). The published proofs (Eynard-Orantin CMP 337 (2015); Fang-Liu-Zong JAMS 33 (2020)) are cited, but their conventions were not transcribed. The normalisation is fixed by exact data, so only the provenance is affected.

#### location

s.2.1 and literature list

#### severity

minor

#### suggested_fix

Cite EO 2015 / FLZ 2020 as the refereed source for BKMP, and say that the normalisation was fixed on Zhou's explicit statement plus the exact genus <= 4 checks.

### [4]

#### description

'About 25,000 comparisons' double-counts nested runs (the g <= 4 coefficients inside the g <= 5 and g <= 6 runs at the same tau, about 4,800). Several listed checks have no stored output (Zhou's (6) literally, genus 0 to n = 8, tau = -5/7, closed forms at five tau, the TRtau runs). The referee re-ran three of them with the author's scripts: 139/139, 1728/1728, 29/29.

#### location

summary and s.0 item 3, s.6 results table

#### severity

cosmetic

#### suggested_fix

Report distinct comparisons and store all outputs in runs\.

### [5]

#### description

TRtau.pm and TRxi.pm share the zeta-coordinate ODE, the kernel expansion and the residue tables: two bases for one method, not independent implementations.

#### location

s.0 item 3 ('Two independent implementations agree'), s.6

#### severity

cosmetic

#### suggested_fix

Rephrase; the referee's TRw.pm (w-basis, direct involution series) provides a methodologically independent confirmation.

### [6]

#### description

Both values are correct for their own y (Zhou's y_Z = -log(1+w) versus y_B = log(1+w)), but the text does not say which y is meant.

#### location

summary ('B = -sigma') vs claims list ('B = sigma')

#### severity

cosmetic

#### suggested_fix

State which y each value refers to.

### [7]

#### description

The x-form uses Route A's Cor. 2.2, which is stated for h >= 2 (or h = 1, r >= 1). The x-form should carry the same restriction (g = 0 is not covered by the cited corollary).

#### location

Theorem F, x-form (s.1.3)

#### severity

cosmetic

#### suggested_fix

Add the range g >= 2, or g = 1 with n >= 1.

### [8]

#### description

The text calls the sigma = -1/4 odd-part statement new in every genus. It is Route A's TR-free Theorem C, as s.7(b) itself notes. Only the sigma = 2 statement is new.

#### location

s.7, closing sentence

#### severity

cosmetic

#### suggested_fix

Restrict the novelty claim to sigma = 2.


## notes

I read Zhou's arXiv:0910.4320v2 in full, through WebFetch of the arXiv PDF. arXiv shows no journal reference for it or for Chen's 0910.3739. All referee scripts and outputs are in G:\My Drive\Projects\Old\refGW\hodge_vertex_conifold\lemmaA\round2\TR1_dictionary\referee\ (TRw.pm, HodgeRef.pm, run_w.pl, check_11_03.pl, cmp_saved_vs_refengine.pl; outputs in runs\). The independent TR code passed 2511 lattice comparisons with 0 failures. The tau = 4 genus-5 run took 21 minutes. No author file was modified. The only new mathematical facts are: Zhou's theorem does not cover (0,3) and (1,1), which s.5 already handles; and P2 holds exactly at genus 7 against Route E's referee engine tables. All mathematical claims stand. The required changes are minor attribution and presentation fixes.

## refuted

- The novelty claim that (6.1) holding identically in (j,tau) at h = 7 is 'new': Route E's referee had already established it without TR (E_formulas/referee.md; brief s.2(g)). The statement itself is true.
- The report's '(7,1): no Hodge table': E_formulas/referee/ref_engine_h7_r1..3.txt exist, and P2 checks exactly against them at genus 7.
- The citation claim that Zhou's Theorem 1 covers all stable (g,n) ('Zhou's theorem is an all-(g,n) statement'): it does not cover (0,3) and (1,1). The proof survives because s.5 computes these directly.

## route

TR1_dictionary (round 2): the dictionary P2 (Lemma B) between the Hodge CohFT and TR on C_tau

## verdict

sound

## report_text

# Referee report: TR1, round 2. The dictionary P2 (Lemma B)

**Referee folder.** `HVC\lemmaA\round2\TR1_dictionary\referee\`, where HVC = `G:\My Drive\Projects\Old\refGW\hodge_vertex_conifold`.
- Scripts: `TRw.pm`, `HodgeRef.pm`, `run_w.pl`, `check_11_03.pl`, `cmp_saved_vs_refengine.pl`.
- Outputs: `referee\runs\*.out`.

**Ground rules followed.**
- All arithmetic is exact (Perl 5.38, Math::BigRat). No wolframscript, no python.
- TR was applied only to C_tau and the framed C^3 vertex curve, i.e. the Hodge CohFT. It was never applied to the resolved conifold x C^2. The same holds for the author's work.
- No file outside the referee folder was modified.

## 0. Verdict: sound

P2 (Lemma B) is correctly proved, modulo the literature, for every tau not in {0,-1} and every stable (g,n), with the brief's normalisation:

    [x^mu]F_{g,n} = sigma^(1-g) prod_i Pi_{mu_i}(tau) Xi_g(mu; tau).

What I did:
- re-derived every step of report s.3-4;
- read Zhou's paper (arXiv:0910.4320v2) in full and checked the transcription equation by equation;
- re-derived the base cases (0,3) and (1,1) by hand, for every tau and also for every e in the family y_e;
- confirmed the dictionary with an independent TR implementation at six new values of tau. It differs from the author's code in coordinate, basis, involution and residue method, and it was compared against Hodge data the author did not use for P2;
- verified P2 directly at genus 7.

No mathematical error was found. The issues are minor:
1. Zhou's theorem does not cover (0,3) and (1,1). The report's own s.5 supplies both cases, so the proof is complete, but s.2.1 and s.4 misstate what Zhou covers.
2. The formula (6.1) at h = 7 is not new. Route E's referee already verified it without TR. The genus-7 Hodge table that the report says is missing exists, and with it P2 holds exactly at genus 7.
3. Presentation: Zhou and Chen are unrefereed preprints; the comparison count is inflated; some checks have no stored output; there are small inconsistencies.

## 1. The literature input: J. Zhou, arXiv:0910.4320v2 (read in full)

I checked the transcription against the paper item by item.

- **W_{g,mu}(a) (p.3).**
  - It carries the sign (-1)^{g+l(mu)}, the factor (a(a+1))^{l(mu)-1}/|Aut mu|, the leaf factors prod_{j=1}^{mu_i-1}(mu_i a + j)/(mu_i - 1)! = Pi_{mu_i}(a), and the class Lambda^vee_g(1) Lambda^vee_g(-a-1) Lambda^vee_g(a).
  - Unstable conventions: (1) gives m^{-2}, (2) gives 1/(m1+m2).
  - The weights are therefore (1, tau, -1-tau) with a = tau, and the leaf normalisation is the binomial C(mu(tau+1)-1, mu-1).
- **(3)-(7) and automorphism factors.**
  - With z_mu = |Aut mu| prod mu_i and m_mu = (1/|Aut mu|) Sum_{S_n} x^{mu_sigma}, the coefficient of an ordered monomial x^m is prod m_i (-1)^{g+n} sigma^{n-1} prod Pi Int(.)/prod(1 - m_i psi_i).
  - This is Zhou's (6), with dphi_b = Sum Pi_m m^{b+1} x^m dx/x.
  - Since <prod tau_b T_g(a)> = (-sigma)^g M[b], one gets omega^{y_Z} = (-1)^n sigma^{g-1+n} Sum_k M[k] prod dxi_k.
  - The transcription is correct, including all |Aut| factors.
- **Curve (8)-(11).**
  - x = y^a - y^{a+1}, y = 1-u, x = u(1-u)^a.
  - With y_Zhou = 1/(1+w): x = w(1+w)^{-a-1}, and the ramification point y = a/(a+1) is w = 1/a.
  - Remark 2.1 ((13)-(14)) only converts to BKMP's curve and is not used. Correct.
- **Recursion (19)-(23).**
  - (19) is W_0 = ln y dx/x. So the TR "y" is y_Z = ln y_Zhou = -log(1+w), and dX = dx/x.
  - (23) is dE_z(y_1) = (1/2)(1/(y_1 - y(q)) - 1/(y_1 - y(q~))) dy_1, which equals (1/2) Int_{q~}^{q} B because B is Moebius invariant.
  - This is the same kernel orientation as the author's code and as review_checks/referee1_trC1.pl. The sign-sensitive n = 1 values -1/16, -1/512, -3/4096 agree.
  - Homogeneity then gives (BKMP): omega^{log(1+w)}_{g,n} = sigma^{g-1+n} Sum M[k] prod dxi_k.
- **The (0,2) term.** Zhou uses W_0 = B - S with S = dx dx/(x-x)^2. The report's (LL) argument removing the S-terms is correct; Zhou's (65) uses the same oddness of S.
- **Range of a.** Zhou assumes a in Z_{>0} (after his (31)). The report uses tau in Z_{>=1} only, then extends by rationality. Correct.

### Coverage gap (minor; repaired by the report's own s.5)

Zhou's Theorem 1 states only that the Bouchard-Marino recursion can be derived from the cut-and-join equation. It does not cover (0,3) and (1,1):

- **At (1,1).** (21) contains W_0(q, q~) = B(q,q~) - dx(q)dx(q~)/(x(q) - x(q~))^2. The second term is infinite because x(q~) = x(q). Zhou's remark after (23) mentions only the W_0(q, y_i) terms.
- **The derivation of (65).** It starts from the cut-and-join equation (68) and expands Phi_{g-1,n+1} and Phi_{g,n-1} by the stable formulas (77) and (79). Those formulas hold only for stable indices, which requires 2g-2+n >= 2 for the target.
- **Initial values.** Zhou's s.3.3 verifies only (19) and (20).

So BKMP at (0,3) and (1,1) has to be computed directly. The report does this in s.5 for all tau, and I re-derived both cases (s.2 below). For 2g-2+n >= 2, Zhou's recursion never contains W_0(q, q~), so BKMP for every stable (g,n) follows by induction on 2g-2+n. Only the sentence "Zhou's theorem is an all-(g,n) statement" (s.4, "Uniformity") needs correcting.

### Literature status

- arXiv lists no journal reference for Zhou (0910.4320) or for Chen (0910.3739).
- Published proofs of BKMP for framed C^3 exist: Eynard-Orantin, CMP 337 (2015), and Fang-Liu-Zong, JAMS 33 (2020). Both are cited, but their conventions were not transcribed.
- The normalisation is pinned down by exact data. The author checked BKMP to genus 4 at three values of tau. I checked it at tau = 4 for g <= 3, and Zhou's (6) literally at tau = 5/2 for g <= 2, n <= 3.
- So only the provenance of the theorem is affected, not the normalisation.

## 2. Re-derivation of the proof (report s.3-5)

**Lemma 1 (loop equations): correct.**
- Since Res_a sigma^*(.) = Res_a(.), one gets omega = [Phi]_-.
- Phi is even, so it has no residue at a.
- (LL) holds because an invariant form f dzeta has f odd.
- (QL) holds because E = (omega - h)(Y - Ybar).

**Lemma 2 (variational formula for a regular deformation): correct.**
- The proof leaves one combinatorial step implicit, and it holds. Apart from the two new terms B(q,r) omega(sigma q, J) and omega(q, J) B(sigma q, r), the terms of E_{g,n+2}(q; J, r) are in bijection with pairs (a term of E_{g,n+1}(q; J), one factor of that term). The bijection appends r to the factor.
- The bijection includes omega_{0,2}(q, p_j) -> omega_{0,3}(q, p_j, r), and at (1,1) it includes omega_{0,2}(q, sigma q) -> omega_{0,3}(q, sigma q, r).
- Those terms are killed by Int Lambda omega_{0,3} = 0, because the bracket eta(-zeta) b(zeta) - eta(zeta) b(-zeta) is odd.
- The identity A Bbar + Abar B = (1/2)(A + Abar)(B + Bbar) - (1/2)(A - Abar)(B - Bbar) is right, and so are the two decompositions built from it.
- The order counts are right:
  - K = O(zeta^-1) dzeta^-1;
  - Omega~/omega~ = O(1);
  - Q, (Y + Ybar)(omega + sigma^*omega) and (Omega + sigma^*Omega)(omega + sigma^*omega) are each O(zeta^2) dzeta^2.

**Lemma 3 (dual cycle): correct.**
- Partial fractions: -(1 - tau w)/(w^2(1+w)) = -1/w^2 + (1+tau)/w - (1+tau)/(1+w).
- Res_0 dxi_k/w = [x^1] xi_k = 1.
- xi_k(-1) = 0 for k >= 1, because D carries the factor 1+w.
- Hence v = -psi/(1-psi) = lam_0 - lam_1.
- The two halves eta_1 = -(1+w)/(w(1 - tau w)) and eta_0 = (1+tau)/(1 - tau w) each have a pole at a. Only their sum, -1/w, is regular.
- For tau < -1 the segment [-1, 0] passes through a. This is harmless: the path is deformed, and only tau in Z_{>=1} is used.

**Lemma 4 (rationality in e, and the limit): correct.**
- Yh(0) = 2 (tau/(1+tau)) (1 + e sigma).
- Homogeneity: omega^{lambda y} = lambda^{2-2g-n} omega^y.
- C^{Y_0 + delta y_B} is regular at delta = 0.

**Lemma 5 (dilaton resummation): correct.**
- Removing the -psi leaves gives the factor (-1)^{K-j} (N_j)_{K-j}, with N_j = 2g-2+n+j > 0.
- Summing over K is then the binomial series.

**Steps 1-2: correct.**
- The limit: e^{2g-2+n} sigma^{g-1+n} (1 + e sigma)^{2-2g-n} H(e sigma/(1+e sigma)) tends to sigma^{1-g} H(1) = sigma^{1-g} X.
- The extension to all tau:
  - C^{Y_0}[k] lies in Q[tau, 1/tau, 1/(1+tau)]. The coefficients u_m come from the ODE with kappa = tau^2/(1+tau), and the only divisions are by c2, a, 1+a, 2 tau^2 and T_kk.
  - The right side lies in Q[1/sigma].
  - Agreement on Z_{>=1} therefore implies an identity in tau.

### Base cases, by hand for every tau and e

For a regular y on a curve with one branch point, I derived two general formulas.

Local data:
- X - X(a) = -c2 zeta^2, and y = y(a) + y1 zeta + y2 zeta^2 + y3 zeta^3 + ... .
- beta0 = b(0,0) = u3 - u2^2 (the Schwarzian).
- e_0 = Res B/zeta and e_1 = Res B/zeta^3.

Formulas:

    omega_{0,3} = e_0^3 / (2 c2 y1),
    omega_{1,1} = (1/(4 c2 y1)) [ (beta0 - y3/(4 y1)) e_0 + e_1/12 ].

Change to the Laplace basis:
- e_0 = tau^2 dxi_0;
- e_1 = 2 c2 tau^2 dxi_1 + (u2^2 - u3) tau^2 dxi_0;
- u2 = tau(tau+2)/(3(1+tau)) and u3 = tau^2 (tau^2 + 13 tau + 13)/(36 (1+tau)^2).

Results:
- **Y_0 = -1/w:** (C_{1,1}[0], C_{1,1}[1]) = (-(2 + 1/sigma)/24, 1/24), and omega_{0,3} = sigma dxi_0^3.
- **y_B = log(1+w):** (C_{1,1}[0], C_{1,1}[1]) = (-(1 + sigma)/24, sigma/24), and omega_{0,3} = sigma^2 dxi_0^3.
- **y_e = y_B + e Y_0:** (C_{1,1}[0], C_{1,1}[1]) = (-[(1+sigma) + e sigma (1 + 2 sigma)]/(24 (1 + e sigma)^2), sigma/(24 (1 + e sigma))). This is exactly Theorem F at (1,1).

All of these agree with the report's s.5. As a numerical cross-check, `check_11_03.pl` compares 384 x-coefficients at 24 values of tau, for both curves, with 0 mismatches.

### Other items checked

- **(0,1):** Lagrange inversion gives [x^{mu+1}] (1+w)^{-tau-1} = -sigma Pi_mu/(mu+1), hence [x^mu] Int Y_0 dX = sigma Pi_mu/(mu(mu+1)).
- **(0,2):** from Zhou's (67) and (74), Phi_{0,2} = -log((w1 - w2)/(x1 - x2)) plus terms in one variable only.
- **Moebius classification:** s^*dX = dX with s not the identity happens iff tau is in {1, -2, -1/2}, because the residues of dX at 0, -1, oo are the CY weights 1, -(1+tau), tau.
- **Theorem F's x-form:** it is Route A's Cor. 2.2 with s = -e sigma/(1+e sigma), u = s/(1+s) = -e sigma and 1 + e sigma = 1 - u. It inherits the range of Cor. 2.2: h >= 2, or h = 1 with r >= 1.

## 3. Independent computations

### 3.1 Independent TR code (`TRw.pm`, `run_w.pl`, `HodgeRef.pm`)

The code shares nothing with the author's TRxi.pm or TRtau.pm. Its method:
- **Coordinate and basis.** It works in the global coordinate w, with basis dw/(w - a)^m.
- **Involution.** s(u) is computed directly from the ODE s s' (a+u)(1+a+u) = u (a+s)(1+a+s), and checked against X(a+s) = X(a+u).
- **Kernel.** It is expanded directly in u = w - a.
- **x-expansion.** By Lagrange inversion, at every lattice point mu >= 1 with |mu| <= 3g-3+2n+2. That set is unisolvent for total degree 3g-3+n, with two redundant layers.
- **Hodge side.** Read with my own parser from:
  - the Mathematica q-tables, evaluated at s = -1 for P2 (the author did not use this data path for P2);
  - the Mathematica mv-tables;
  - the Route E author and referee engines.

All results below are exact, with 0 failures and 0 symmetry or pole-order violations:

| identity | y | tau (sigma) | (g,n) compared | Hodge data | lattice points |
|---|---|---|---|---|---|
| P2 | -1/w | 1 (2) | g <= 1, n <= 3; (2,1) | genus-0 formula, engine, q-table | 96 |
| P2 | -1/w | 4 (20) | all g <= 5 with n <= 3, g + n <= 6: (0,3), (1,1)-(1,3), (2,1)-(2,3), (3,1)-(3,3), (4,1), (4,2), (5,1) | q-table at s = -1 (g = 2-4, r <= 2); engine (g = 1, r = 3, g = 5) | 867 |
| P2 | -1/w | 3/5 (24/25) | g <= 3, n <= 2; and (4,1) | engine (g = 1), q-table | 138 + 34 |
| P2 | -1/w | -7/3 (28/9) | g <= 3, n <= 2 | engine (g = 1), q-table | 138 |
| P2 | -1/w | -2/7 (-10/49) | g <= 3, n <= 2 | engine (g = 1), q-table; Pi_7 = 0, so 122 of the 138 points are informative | 138 |
| P2 | -1/w | 5/2 (35/4) | g <= 2, n <= 3 | Route E referee engine | 293 |
| BKMP | log(1+w) | 4 | g <= 3, n <= 2 | mv-table | 138 |
| Zhou's (6) literally | -log(1+w) | 5/2 | g <= 2, n <= 3 | mv-table, with sign (-1)^n | 293 |
| Theorem F | log(1+w) - e/w | 3/5, e = -2 | g = 2, 3; n <= 2 | q-table at s = -e sigma/(1+e sigma) | 119 |
| Theorem F | log(1+w) - e/w | 4, e = -1/21 | g = 2, 3; n <= 2 | q-table | 119 |
| base cases | -1/w and log(1+w) | 24 values | (0,3), (1,1) | closed forms | 384 |

The earlier tau = 4 run for g <= 3, n <= 2 (138 points) is contained in the 867-point row.

### 3.2 Genus 7

The report says "(7,1): no Hodge table". In fact `E_formulas\referee\ref_engine_h7_r1..3.txt` exist. They are Route E's referee engine tables and use no TR.

- **Comparison.** `cmp_saved_vs_refengine.pl` compares the author's saved one-leaf TR coefficients (`runs\save_tau*_g7n1*.txt`) with sigma^{1-g} X_{g,1}[k] from these tables, at tau = 1, 2, 3, 1/2, -1/2, -1/3, 2/3 and g = 1..7. Result: 539 coefficients, 140 of them at genus 7, with 0 mismatches.
- **Consequence.** P2 is verified directly at genus 7, n = 1.
- **TR-free values.** The same tables reproduce, without any TR, the tau = 1 values of referee_hurwitz_1: -1/512, -3/4096, -225/262144, -2205/1048576 and -297675/33554432. They also reproduce the author's new value -15436575/268435456 at h = 7, j = 13.

### 3.3 Author checks with no stored output, re-run with the author's own scripts

| check | result |
|---|---|
| Zhou's (6) literally, tau = 2/3, g <= 2, n <= 3 | 139/139 |
| genus 0 up to n = 8, tau = 2/3 | 1728/1728 |
| tau = -5/7, g <= 1 | 29/29 |

## 4. Consequences (report s.7-8)

**Theorem 3.1 (the TR statement) re-checked.**
- The kernel has the product form, and (*) holds.
- The residue theorem applies to f_m, with sigma^* f_m = f_m.
- For g >= 2 there is no pole at v = -1.
- The degree count is right.
- The tau = -1/2 variant also holds: s(v) = -1-v fixes oo, x o s = -x, and degree 0 for B(sigma v, v_j) is enough.

With P2 this gives:
- **Lemma A at sigma = 2, every genus:** correct. It is new as an all-genus statement.
- **Lemma A at sigma = -1/4, R with odd parts:** correct, but it is already Route A's TR-free Theorem C. The report says so in s.7(b).

**Other consequences.**
- **proof_status 1.7(b) and (c):** now unconditional, modulo P2's literature and the other inputs that 1.7(b) cites.
- **(6.1) at h <= 7:**
  - The argument that both sides have sigma-degree at most h is correct: the power sum p_k has sigma-degree floor(k/2).
  - However, h = 7 is not new. Route E's referee verified (6.1) identically at h = 7 without TR (E_formulas\referee.md; brief s.2(g)).
  - After 3.2, the genus-7 input is in any case available from Hodge-side data alone.
- **s.8 (what P2 does not give):** accurate.

## 5. Issues

1. **(Minor) Zhou's coverage.** Zhou's Theorem 1 does not cover the base cases (0,3) and (1,1), as explained in s.1. The proof is complete only because s.5 computes these cases directly. The report should say this explicitly.
2. **(Minor) Genus 7.**
   - "(6.1) ... new at h = 7" is not new; Route E's referee established it without TR.
   - "(7,1): no Hodge table" is wrong; `ref_engine_h7_r*.txt` exist.
   - With those tables, P2 holds exactly at genus 7 (3.2).
3. **(Minor) Literature status.** Zhou and Chen are arXiv preprints with no journal reference. The report should state that the refereed sources are Eynard-Orantin (2015) and Fang-Liu-Zong (2020), and that the normalisation was fixed on Zhou's statement together with exact data.
4. **(Cosmetic) Bookkeeping.**
   - "About 25,000 comparisons" counts nested runs twice. The g <= 4 runs sit inside the g <= 5 and g <= 6 runs at the same tau, about 4,800 comparisons.
   - Several of the listed checks have no stored output: Zhou's (6) literally, genus 0 to n = 8, the tau = -5/7 run, the closed forms at five tau, and the TRtau runs.
   - I re-ran three of them (3.3); all pass.
5. **(Cosmetic) "Two independent implementations".** TRtau.pm and TRxi.pm share the zeta-coordinate ODE, the kernel expansion and the residue tables. They are two bases for one method, so the phrase overstates their independence. My TRw.pm is methodologically independent.
6. **(Cosmetic) Value of B.** The summary says B = -sigma and the claims list says B = sigma. Each is right for its own y (y_Z versus y_B); the report should say which y is meant.
7. **(Cosmetic) Range of the x-form.** Theorem F's x-form should carry the range of Route A's Cor. 2.2: g >= 2, or g = 1 with n >= 1.
8. **(Cosmetic) Novelty wording in s.7.** The closing sentence of s.7 presents the sigma = -1/4 odd-part statement as new. It is Route A's TR-free Theorem C, as s.7(b) itself says. Only the sigma = 2 statement is new.

## 6. Files in `referee\`

- `TRw.pm`: independent exact TR on C_tau (w-basis, involution series, residues in u, x-expansion by Lagrange inversion).
- `HodgeRef.pm`: independent readers for the q- and mv-tables and the engines.
- `run_w.pl TAU G N [base|bkmp|zhou|mix] [e] [q|eng|ref]`, with env `NCMP` for the maximum n compared.
- `check_11_03.pl`: base cases at 24 values of tau.
- `cmp_saved_vs_refengine.pl TAU SAVEFILE`: the genus-7 comparison.
- `runs\`: all outputs listed above.

