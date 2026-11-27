# Lemma A round 2: referee2 TR1 (saved by the main session from the workflow journal, 2026-10-31)

Subagents cannot write report files, so this file holds the agent's structured return value: the fields first, then the full report text (report_text) verbatim.

## confirmed

- THEOREM P2 (Lemma B) for every tau not in {0,-1} and every stable (g,n), with the brief's normalisation [x^mu]F = sigma^(1-g) prod Pi_{mu_i} Xi_g(mu). Status proved-modulo-literature is justified (Zhou/Chen/EO-2015 BKMP for framed C^3, local TR symmetry, dilaton equation).
- Logic of s.4 re-derived. (i) Iterating (V) at e = 0 gives the K-th Taylor coefficient sigma^(g-1+n) sigma^K <Omega prod psi^k prod v(psi)> via BKMP. (ii) Lemma 5 resums this to sigma^(g-1+n)(1+e sigma)^(2-2g-n) H(e sigma/(1+e sigma)). (iii) Both sides are rational in e with the single pole e = -1/sigma (Yh(0) = 2(tau/(1+tau))(1+e sigma)), so they coincide. (iv) The homogeneity limit e^(2g-2+n) C^{y_e} -> C^{Y_0} gives sigma^(1-g) H(1;k) = sigma^(1-g) X[k]. (v) Agreement at tau in Z_{>=1} plus C^{Y_0} in Q[tau,1/tau,1/(1+tau)] gives every tau. No circularity: BKMP is used only at e = 0, and P2 only at e = oo.
- Lemma 3 re-derived by hand: Y_0 dX = -dw/w^2 + (1+tau)/w dw - (1+tau)/(1+w) dw; pairing v_0 = 0, v_k = -1 (k >= 1), i.e. v = -psi/(1-psi); the two halves are irregular, their sum -1/w is regular.
- Variational formula (Lemma 2 + Lemma 3) tested exactly with the referee's independent TR code R2TR.pm (dual numbers): tau = 2, e0 = 1/3, g <= 5 (1277 coefficients); tau = -2/5, e0 = 2 and tau = -2, e0 = -3/7, g <= 4 (299 each); homogeneity/dilaton identity at the base curve, tau = 1/2 (598). 0 failures.
- BKMP normalisation (sigma^(g-1+n) M[k] for y = log(1+w)) tested against the referee's TR-free Hodge engine RefHodge at tau = 2 up to genus 5, all (g,n) with g+n <= 6: 1277/1277.
- Lemma 5 (dilaton resummation) checked as a power series in e, with v-leaves as explicit marked points (RefHodge): (2,1) at tau = 2/3 to order e^4, (2,2) at tau = 3 to e^3, (3,1) at tau = -2/5 to e^2: 133/133.
- Theorem F (e-family) at genus 5, 6, 7 (one leaf), identity in j: (tau,e) = (3/5,-2/3) g = 5, (-2/5,5/2) g = 6, (4/3,-1/7) g = 7; BKMP and P2 at the same points; 0 mismatches.
- P2 far beyond the author's range with the author's TRxi.pm (md5-identical copy) against TR-free tables (Route E engine + Route E referee engine; genus 0 formula): tau = 1, 2, 3/4, all (g,n) with g+n <= 9, g <= 8, n >= 1, including (6,3), (7,2) and (8,1): 38952 Hodge comparisons per tau, 0 failures, 0 symmetry violations. (8,1) agrees with (6.1) at 26 j-values, where (6.1) is still conjectural.
- P2 end-to-end in x-form ([x^mu]F by Lagrange inversion of global forms, no xi-basis) with the independent R2TR.pm at tau = -2 and -3/2, g <= 4: 1935/1935 each.
- Unstable cases: (0,1) [x^mu] Int Y_0 dX = sigma Pi_mu/(mu(mu+1)) and (0,2) sigma Pi_m1 Pi_m2/(m1+m2), checked at tau = 1, 2/3, -1/2, 3, -5/7, 7/3 with m1+m2 <= 10: 330/330.
- Moebius classification (s^*dX = dX, s != id iff tau in {1,-2,-1/2}) re-derived: residues 1, -(1+tau), tau sum to 0, so they cannot all coincide; with two equal, the swap is unique and correct.
- Consequences: Lemma A at sigma = 2 in every genus (Theorem 3.1 + P2; Pi_mu(1) = C(2mu-1,mu-1) != 0); the odd-part statement at sigma = -1/4 (= Route A Theorem C); proof_status 1.7(b),(c) and 1.8 now unconditional; (6.1) identically for h <= 7 (8 sigma values, sigma-degree <= h). The degree bound |mu| <= 2g-3+n holds in the TR data at tau = 1 and also at tau = 2 and 3/4 through genus 8, one leaf, which is consistent with Lemma A (known for h <= 20).
- The final statements are not overclaimed: the report states correctly that P2 gives neither (OL) nor Lemma A for all tau, and that Lemma D for tau >= 2 remains open.
- Rules respected: TR is used only for the Hodge CohFT (C_tau); none for the resolved conifold x C^2; no wolframscript.

## issues

### [1]

#### description

EO07 s.4 proves symmetry for a compact spectral curve with x and y meromorphic. Here the recursion is taken with respect to dX = dx/x, and y_e contains log(1+w), so (x, y_e) is not an algebraic curve for any e != oo, even at integer tau. Lemma 1 and Lemma 2 need symmetry at every e near 0. The statement that EO07 'applies verbatim' is plausible: with a single branch point, EO's exchange-of-residues argument is local at a. But that is not what EO07 states. The correct sources are the local-curve versions: Dunin-Barkowski-Orantin-Shadrin-Spitz CMP 328 (2014), which identifies TR on a local curve with Givental's formula and so gives symmetry, or Borot-Eynard-Orantin CNTP 9 (2015) (abstract loop equations). All correlators computed by both referees are exactly symmetric, roughly 120,000 coefficients up to genus 8, so the gap concerns citation only.

#### location

report s.2.2 (symmetry) and Lemma 1 / Lemma 2, used along the whole family y_e = log(1+w) + e(-1/w)

#### severity

minor

#### suggested_fix

Cite symmetry for local spectral curves with one simple branch point (DOSS 2014, or BEO 2015 with the precise theorem number), and say that this is the form of the theorem being used.

### [2]

#### description

I confirm referee 1's point. Zhou's Theorem 1 is derived for 2g-2+n >= 2; the (1,1) term of his (21) involves W_0(q,q~), whose dx dx/(x-x)^2 part is singular. The proof is still complete, because s.5 computes (0,3) and (1,1) directly for all tau and all y_e. The sentence 'Zhou's theorem is an all-(g,n) statement' is nevertheless false as written.

#### location

report s.2.1 and s.4 'Uniformity in the genus' (already raised by referee 1, issue 1)

#### severity

minor

#### suggested_fix

Say: base cases (0,3) and (1,1) by s.5; Zhou for 2g-2+n >= 2.

### [3]

#### description

For real tau <= -1 the branch point a = 1/tau lies on [-1, 0]. For example tau = -2 gives a = -1/2. The cycle then has to be deformed around U. This is harmless: the correlators have zero residue at a by Lemma 1(a), and Step 1 uses only tau in Z_{>=1}, where a lies in (0,1]. Still, the text states Lemma 3 for every tau without the caveat.

#### location

Lemma 3 / Lemma 2 hypothesis 'path from -1 to 0 avoiding Ubar', and the claim that Theorem F holds for all tau

#### severity

cosmetic

#### suggested_fix

State Lemma 3 for tau in Z_{>=1}, which is all that is used, or add 'path homotopic in P^1 minus Ubar; independent of the path because Res_a = 0'.

### [4]

#### description

The report does not show that the sigma = 2 condition is linearly independent of Route A's Z(D) divisibility conditions (brief s.5), so 'one of the G-1' may overcount when the two are combined. This affects bookkeeping for (OL), not any claim about P2.

#### location

s.7(f) and summary ('sigma = 2 is one of the G-1 needed conditions')

#### severity

cosmetic

#### suggested_fix

Say 'one condition (sigma = 2)'; its independence from Route A's Z(D) conditions is not checked.

### [5]

#### description

I confirm referee 1's issues 2, 4 and 5. (6.1) at h = 7 was already established without TR by Route E's referee. The comparison count double-counts nested runs. TRtau.pm and TRxi.pm share the method. R2TR.pm, below, is a genuinely independent implementation.

#### location

summary and claims ((6.1) at h = 7 'new'; 'about 25,000 comparisons'; TRtau.pm vs TRxi.pm called 'two independent implementations')

#### severity

cosmetic

#### suggested_fix

As in referee 1.


## notes

I am the second, adversarial referee. This session resumes earlier referee-2 work. All scripts and outputs are in G:\My Drive\Projects\Old\refGW\hodge_vertex_conifold\lemmaA\round2\TR1_dictionary\referee2\ (R2TR.pm is an independent TR implementation written from scratch; RefHodge.pm/RefWK.pm are TR-free Hodge engines; r2_big.pl, r2_post.pl, r2_check.pl, r2_thF.pl, r2_lemma5.pl, r2_unstable.pl; outputs in runs\). Runs added in this session: runs\var_taum2_5_e2_S6G4.out, runs\var_taum2_em3_7_S6G4.out, runs\vbase_tau1_2_S6G4.out (0 failures); r2_lemma5.pl and r2_unstable.pl results are quoted in the report. No author file was modified. Despite assuming the argument was wrong, I found no mathematical failure. The issues concern the precision of two literature citations (EO symmetry outside its stated hypotheses; Zhou's base cases) and bookkeeping.

## refuted

- 'Zhou's theorem is an all-(g,n) statement': it covers 2g-2+n >= 2 only (as referee 1 found); the proof survives through s.5.
- 'EO07's symmetry theorem applies verbatim': EO07's hypotheses (meromorphic x, y on a compact curve) fail for y_e = log(1+w) + e Y_0; a local-curve symmetry theorem is needed (it exists: DOSS 2014 / BEO 2015).
- Novelty of (6.1) at h = 7 (referee 1's refutation confirmed).

## route

TR1_dictionary (round 2): the dictionary P2 (Lemma B) between the Hodge CohFT and TR on C_tau, second independent referee

## verdict

sound

## report_text

# Referee 2 report: TR1 round 2, the dictionary P2 (Lemma B)

**Referee folder.** `HVC\lemmaA\round2\TR1_dictionary\referee2\`, where HVC = `G:\My Drive\Projects\Old\refGW\hodge_vertex_conifold`.

**Code.**
- `R2TR.pm`: an independent exact TR implementation, written from scratch with no code shared with TRxi or TRtau. It uses zeta by Lagrange inversion from the Taylor series of X, a global basis E_l, the output read off as the principal part of E/omega~, and optional dual numbers.
- `RefHodge.pm` / `RefWK.pm`: a TR-free Hodge engine with no kappa classes, dilaton and v-leaves as explicit marked points, and Mumford's formula.
- `TRxi.pm`: the author's module, copied unchanged (md5 7de3a1e9...).
- Drivers: `r2_big.pl`, `r2_post.pl`, `r2_check.pl` (modes p2, bkmp, var, vbase), `r2_thF.pl`, `r2_lemma5.pl`, `r2_unstable.pl`.
- Outputs: `runs\`.

**Ground rules followed.** All arithmetic is exact (Perl, Math::BigRat). There was no wolframscript. TR was used only on C_tau (the Hodge CohFT); nothing was done for the resolved conifold x C^2.

## 0. Verdict: sound

I assumed the argument was wrong and attacked every all-genus step, numerically and logically. No mathematical failure was found. P2 holds as stated: for every tau not in {0,-1} and every stable (g,n), [x^mu]F_{g,n} = sigma^(1-g) prod Pi_{mu_i}(tau) Xi_g(mu;tau). It is proved modulo the literature:
- BKMP for the framed one-leg vertex (Zhou; Chen; Eynard-Orantin 2015);
- symmetry of TR correlators;
- the dilaton equation.

The defects:
- **Two literature statements are used outside their stated hypotheses (minor).**
  - The EO07 symmetry theorem is cited for non-algebraic y.
  - Zhou's theorem is cited for all (g,n), including (0,3) and (1,1).
  - Correct replacements exist (see s.3), so the proof is unaffected.
- **Bookkeeping and novelty (cosmetic).** These agree with referee 1.

## 1. All-genus steps, and numerical tests beyond the author's range

Each step asserted for all (g,n) is listed with the independent test I ran. All tests are exact.

| step | test (referee 2) | result |
|---|---|---|
| P2, xi-form, all coefficients | author's TRxi (unchanged) against TR-free tables (Route E engine and Route E referee engine; genus 0: [j^k] prod (J-i)); all (g,n) with g+n <= 9, g <= 8, n >= 1, at tau = 1, 2, 3/4 | 38952 Hodge comparisons per tau, 0 failures, 0 symmetry violations, nothing outside \|k\| <= 3g-3+n. Covers (6,3), (7,2), (8,1), beyond the author's g <= 6, n <= 3 |
| P2, x-form end-to-end | independent R2TR, Lagrange inversion of global forms, no xi-basis; tau = -2, -3/2; g <= 4 | 1935/1935 each |
| BKMP normalisation sigma^(g-1+n) M[k] | R2TR, y = log(1+w), against RefHodge; tau = 2; g+n <= 6, g <= 5 | 1277/1277 |
| Lemma 2 + Lemma 3 (variational formula, pairing v = -psi/(1-psi)) | R2TR with dual numbers: d/de C_{g,n}[k] = Sum_l C_{g,n+1}[k,l] v_l | tau = 2, e0 = 1/3, g <= 5: 1277/1277. tau = -2/5, e0 = 2: 299/299. tau = -2, e0 = -3/7: 299/299 |
| homogeneity / dilaton identity at Y_0 | R2TR vbase; tau = 1/2; g <= 4 | 598/598 |
| Lemma 5 (dilaton resummation) | RefHodge with explicit v-leaves, power series in e | (2,1) at tau = 2/3 to e^4; (2,2) at tau = 3 to e^3; (3,1) at tau = -2/5 to e^2: 133/133 |
| Theorem F (e-family), one leaf, identity in j | TRxi y = mix, bkmp, base against RefHodge | (3/5, -2/3) at g = 5; (-2/5, 5/2) at g = 6; (4/3, -1/7) at g = 7: 0 mismatches. BKMP and P2 at the same tau also 0 mismatches |
| unstable (0,1), (0,2) | Lagrange inversion | tau = 1, 2/3, -1/2, 3, -5/7, 7/3; m1+m2 <= 10: 330/330 |

By-products:
- **(6.1) at h = 8.** (6.1) is still conjectural at h = 8. The genus-8 one-leaf TR polynomial agrees with it at tau = 1, 2, 3/4, at 26 values of j each.
- **Degree bound.** At tau = 1, 2 and 3/4, the TR correlators satisfy the degree bound [x^mu]F = 0 for \|mu\| <= 2g-3+n in every computed case, through genus 8 for one leaf. Every coefficient of the first layer \|mu\| = 2g-2+n is non-zero. This agrees with Lemma A, which is known for h <= 20, and gives no counterexample.

## 2. Logic

**Lemma 1 (loop equations).** Re-derived.
- Phi = E/omega~ is odd under the local involution.
- omega = [Phi]_- in the global coordinate w. The pole of dw_{p0}/(w_{p0} - w_q) at p0 = oo has residue -Res_a Phi = 0, so omega is holomorphic at oo.
- (LL) and (QL) follow as stated.
- Lemma 1 uses symmetry; see s.3.

**Lemma 2 (variational formula).**
- The bijection between the e-derivative terms of E and the terms of E_{g,n+2}(q;J,r) is correct. The only new terms are B(q,r) omega(sigma q) and omega(q) B(sigma q,r).
- The identity (**) reduces, through the A Bbar + Abar B decomposition, to integrands that are O(zeta) dzeta.
- The base case Int Lambda omega_{0,3} = 0 is correct: the bracket eta(-zeta) b(zeta) - eta(zeta) b(-zeta) is odd. It extends to q near a by meromorphy in (z1, z2).
- The exact dual-number tests above, which include negative and non-integer tau, confirm the formula.

**Lemma 3 (dual cycle).** Re-derived by partial fractions:

    -(1 - tau w)/(w^2 (1+w)) = -1/w^2 + (1+tau)/w - (1+tau)/(1+w).

- The pairing is Res_{r->0} dxi_k/w = [x^1] xi_k = Pi_1 = 1.
- Also xi_0(-1) = -1/(1+tau), and xi_k(-1) = 0 for k >= 1, because D = w(1+w)/(1 - tau w) d/dw vanishes at -1.
- Caveat (cosmetic): for real tau <= -1 the point a = 1/tau lies on [-1, 0]. The path must then avoid U, which is possible because Res_a = 0. Step 1 uses only tau >= 1, where a lies in (0, 1].

**Lemma 4 (rationality in e).** Correct. 1/Yh has coefficients polynomial in e over powers of Yh(0) = 2(tau/(1+tau))(1 + e sigma). Every other ingredient (B, the basis, T) is independent of e. Homogeneity holds: omega^{lambda y} = lambda^(2-2g-n) omega^y, checked at (0,3) and (1,1).

**Lemma 5 (dilaton resummation).** Correct: Sum_l ((-e sigma)^l/l!) (N)_l = (1 + e sigma)^(-N). Checked exactly (s.1).

**Step 1.**
- Iterating (V) at e = 0 needs (V) for every (g, n+K). That is a single induction, so it is available.
- At e = 0, BKMP applies to (g, n+K). This is the only place BKMP enters. P2 enters only as the e -> oo limit, so there is no circularity.
- Identifying the two rational functions is legitimate: both are regular at e = 0, and their only pole is at e = -1/sigma.
- The limit produces sigma^(g-1+n) sigma^(2-2g-n) H(1;k) = sigma^(1-g) X[k]. The exponents match, so the factor is independent of n.

**Step 2 (extension to all tau).**
- The ODE coefficients for u, c2 = tau^3/(2(1+tau)), a = 1/tau and Yh_0(0) = 2 tau^2 all lie in Q[tau, 1/tau, 1/(1+tau)].
- On the Hodge side, the elementary symmetric functions of the weights (1, tau, -1-tau) are s1 = 0, s2 = -(1+sigma), s3 = -sigma. So Omega is a Laurent polynomial in sigma.
- Agreement on Z_{>=1} therefore gives an identity.

**Base cases.** Zhou does not give (0,3) and (1,1). Report s.5 supplies them for all tau and e; referee 1 re-derived them by hand.

**Consequences against brief s.3 / s.5.**
- **Lemma A at sigma = 2 in every genus.** This follows from Theorem 3.1 (degree bound at tau = 1) plus P2, using Pi_mu(1) != 0. It is a specialisation only, and the report says so.
- **sigma = -1/4.** The odd-part statement equals Route A's TR-free Theorem C.
- **proof_status 1.7(b): Xi_3(3).** The conditions at sigma = 2 and sigma = -1/4 have determinant -9.
- **proof_status 1.7(b): Xi_3(2).** This is already Route A's Theorem A.
- **proof_status 1.7(c).** Distinct integer tau give distinct non-zero sigma, as required.
- **(6.1) for h <= 7.** Eight sigma values against sigma-degree <= h is enough.
- **(OL).** The report correctly does not claim it.
  - P2 supplies one condition (sigma = 2) on sigma^G Xi_G(D+1).
  - Whether that condition is independent of Route A's Z(D) conditions is not shown. So "one of the G-1" should read "one condition".
  - tau = -1/2 contributes nothing to (OL), since D+1 is even.

The all-genus claim is exactly P2, plus Lemma A at sigma = 2 (and at sigma = -1/4 for odd parts). It is correct.

## 3. Literature statements and their hypotheses

**(a) BKMP / Zhou (arXiv:0910.4320).**
- The normalisation and transcription were checked equation by equation by referee 1. I checked them independently against RefHodge at tau = 2 up to genus 5.
- Hypotheses:
  - a in Z_{>0}: respected, since Step 1 uses only these values.
  - 2g-2+n >= 2: not stated by the report, but the gap is covered by s.5.
- Zhou's and Chen's papers are arXiv preprints. The refereed proofs are Eynard-Orantin, CMP 337 (2015), and Fang-Liu-Zong, JAMS 33 (2020); cite those as well.

**(b) Symmetry of correlators.**
- EO07, s.4, is stated for a compact curve with meromorphic x and y.
- Lemmas 1-2 need symmetry for y_e = log(1+w) + e(-1/w), with the recursion taken with respect to dX = dx/x. For e != oo this is not an algebraic curve, even at integer tau.
- With one simple branch point the symmetry argument is local, so the claim is true. The correct citation is:
  - the local-spectral-curve theory of DOSS, CMP 328 (2014), which identifies TR with Givental's formula and so gives symmetric correlators; or
  - Borot-Eynard-Orantin, CNTP 9 (2015).
- The roughly 120,000 exactly symmetric correlator coefficients computed by both referees support this.
- Severity: minor (citation).

**(c) Dilaton equation for pulled-back classes.** Omega is pulled back along the forgetful maps (pi^* lambda_i = lambda_i), so the equation applies; standard.

**(d) Theorem 3.1 (proof_M_hurwitz s.3).** It was checked line by line by referee_hurwitz_1. My genus-8 one-leaf data at tau = 1 confirm its degree bound.

## 4. Issues

1. **(minor)** EO07 symmetry is cited outside its hypotheses. Cite DOSS 2014 or BEO 2015, the local-curve versions.
2. **(minor, as referee 1)** Zhou covers 2g-2+n >= 2 only. Say that (0,3) and (1,1) come from s.5.
3. **(cosmetic)** The Lemma 3 path for real tau <= -1 needs the caveat above.
4. **(cosmetic)** "One of the G-1 conditions" for (OL): its independence from Route A's Z(D) conditions is not shown.
5. **(cosmetic, as referee 1)**
   - (6.1) at h = 7 is not new.
   - The comparison count double-counts nested runs.
   - TRtau and TRxi are not independent implementations. R2TR.pm is.

## 5. Conclusion

P2 (Lemma B) is a theorem modulo the literature in every genus and for every tau not in {0,-1}. Theorem F and the low and unstable cases hold. The all-genus consequences, Lemma A at sigma = 2 and the odd-part statement at sigma = -1/4, follow. Numerically, P2 is now confirmed exactly through genus 8 (one leaf), (7,2) and (6,3) at three values of tau. The variational formula and Lemma 5, the two steps specific to this proof, have also been confirmed exactly at new and negative tau.

None of this proves (OL) or Lemma A for all tau. The open input is still Lemma D for tau >= 2, as the report itself says.
