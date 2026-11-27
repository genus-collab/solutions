<!-- workflow wf_e928f771-9e3 (7 tasks, 5 independent verifiers, synthesis). Scripts/outputs: <user home>/AppData/Local/Temp/claude/G--My-Drive-Projects-Old-refGW/37f24460-3601-4c85-b123-5156e19eb475/scratchpad/lp2_round1/. Main session spot-checked BS v2 eq:Rcf (line 5326) and Upsilon (line 4282). -->
# Local P2 regularity: round-1 status (2 Nov 2026)

Seven tasks ran. Five of them had an independent verifier, who wrote fresh exact Perl/BigRat scripts and reread BS v2. Two had no verifier: the frame-B skeleton and the v1/v2 diff. Status labels below are **proof** (proved by hand), **exact-to-N** (exact computer check to order N) and **open**.

## Stage by stage

**Stage 1: F2 (cfasym).** Reporter's claims confirmed by the verifier.
- **What cfasym gives.** The theorem itself states only the leading polar coefficient, (-3)^{g-1} β_{2g-2} / t_cf^{2g-2}. The full gap is still conjectural in v2 (rmk:cfgap); BS check it only at g=2.
- F2 comes from eq:Rcf, inside the proof: R^cf_00 = -ξ Υ(z/e1, z/e2) / (√(3e1e2E) (E+z)) + O(ξ²). This is a full z-series.
- Two facts finish the step:
  - Υ(z/e1, z/e2) = h_w for w = (E, -e1, -e2). This is an identity of the definitions.
  - The identity row does not depend on the polarisation (prop:fingenR).
- Together they make the frame-B identity row at t_c=0 equal to E/(E+z), to all orders in z.
- BS's derivation (5249–5273) actually gives an exact closed form in ξ: Σ_l r_l (vz)^l / Π_{i=1}^{l+1}(1+iz/E). This is our KEY FORMULA_P2, with ρ_l = (-1)^l r_l. Both agents computed r_l exactly: 1, 1/9, 2/81, 14/2187, …
- Verifier's stronger check: at finite rational L, the stored cPP_z9 row equals this closed form, with v = 3/(E+3L), exact to z⁹ at 4 weight pairs. The z → -z control fails.
- No TR anywhere in the chain. TR appears only in sec:refTR and rmk:refmm.
- **Status:** proof, modulo the genus-0 inputs (Givental's mirror theorem, the CCIT normalisation N_i) and a continuity lemma for passing to T_F.
- Verifier's correction: the "3d vs iterated Laplace expansion" issue is not a real gap. The Hessian and the non-Gaussian part split, so the operators commute. N_i does not depend on α, which also simplifies the continuity lemma.

**Stage 2: closure lemma.** Confirmed by the verifier.
- Hand identity, valid for any invertible c1: C1i Th v = (A_0 − z² c1 Rd E_21) v.
  - A_0 is z-free, η-self-adjoint and independent of XX.
  - Rd is the Riccati defect.
- On shell, Rd = 0, because C1cf = −₂F₁(1/3, 2/3; 1; Del) solves L2b. So the conifold-frame QDE holds to all orders.
- This explains the earlier XXSHIFT observations, including the off-shell pole coefficient e1e2E²(XX(0) + 2/9).
- It does not need algebraic independence.
- **Status:** proof, conditional on H1 (the code rows are frame-B BS R_0i satisfying eq:pf). H1 is exact to z⁹ only.

**Stage 3: Lemma H2.** Confirmed by the verifier.
- On shell, σ₂ has no ψ-linear term through t_c⁷ at 3 weight pairs. It satisfies the sharper bound σ₂^(k) ∈ Λ_k, and in the true gauge it has no constants.
- The ψ-coefficient vanishes exactly because XX = θ log C1cf, i.e. by the definition of the conifold polarisation. The reporter derives this from (I1) and (I2).
- **Disagreement (minor).** The verifier notes that the reporter's off-shell controls (XXSHIFT) do not give a consistent buildRmatrix section. With a consistent deformation c1 = −wc(1 + Q·Del), the conclusion still holds, with a1 = E·Del·wc²·θf / t_c.
- The verifier proposed a shorter all-orders route that needs only Z∞ (lim z ρ P̂0 = E) and L1. It did not recheck (I1) and (I2).
- **Status:** exact to t_c⁷. The all-orders argument depends on Z∞, L1 and F2; Z∞ is unproved. Removing the constants also needs symplecticity.

**Residue lemma (E_00).** Confirmed by the verifier.
- Pure algebra on the basis {1, λ_j}. It holds for any rank and any constant symmetric pairing.
- Under H1–H3 it gives regularity of E_00, L3 support and rank one at s⁰, with A_11 = −e1e2E⁴/27.
- It reproduces the code's E_00 exactly, 45 entries, s⁻² to s⁶, at 3 weight pairs. A conifold control gives 644/0. The off-shell controls fail, as they should.
- **Status:** proof. It makes the rank-3 flow argument unnecessary for E_00. Its inputs H1–H3 are open. It does not cover the mixed edges E_0j.

**Stage 5: regular colours.** Confirmed by the verifier.
- Exact recursion in s = 1/(L + α), with closed forms for Op1, Op2 and Op3. Induction gives S_k ∈ sQ[s] of degree ≤ k+1, unique up to c(z).
- Consequences: analyticity (R) at colours 1 and 2, and the stored cPP form for symbolic ε.
- C1i and XX drop out of the orthogonality relation; the cross-colour terms vanish via a Lagrange identity.
- **Status:**
  - The closed forms and the induction are proved.
  - The Lagrange identity was not rechecked by the verifier; it is consistent with the exact-to-z⁹ orthogonality checks.
  - K(z) = 1 reduces to stage 2.
  - The WKB identification on T_F is open; it is the same issue as H1.

**Frame-B skeleton for the gap corollary.** No verifier, so this is unchecked.
- Proposed proof of Prop 13.2_P2: Cf(GW_g) equals the frame-B sum, in 8 steps.
  - Givental–Teleman on T_max.
  - Rows from eq:canR and eq:RHb.
  - The normalisation constants, by porting conifold 13.2(c).
  - Passage T_max → T_F by an implicit-function-theorem argument (W_i ≠ 0).
  - Frame B, by porting 13.2(e).
  - The polarisation change Cf is a ring homomorphism, so it applies graph by graph and needs no on-shell input.
  - F2 from eq:Rcf.
  - The gap: the single-vertex coefficient is (-3)^{g-1} β_{2g-2}.
- Gaps it flags:
  - N1: BS assert the normalisation constants without proof.
  - N2: analyticity on the diagonal α_i = E/3.
  - N4: lem:modrel E = (D−X)/2 is only a sketch in BS. It is checked to y⁴⁰.
  - N5: a dictionary from the code to BS.
  - N6: the gap needs Lemma A up to h = g, not g−1.
- One tension: frame-B uses algebraic independence to lift the identities to G_{T_F}[A^{±1}, F]. Closure shows this is not needed if eq:pf is cited.

## Convention changes in v2 (diff task, no verifier)

The changes found all agree with our code. lp2_lib.m, the README and the (-3)^{g-1} ratios need no change.
- A = y d_y Π, with no log.
- F_{k,g} is the [(e1+e2)^{2k}] coefficient, with pole Π^{2−2g−2k}.
- The gap coefficient is β_{2g−2}, not β_{2g}.
- D_cf = X/A_cf².
- Measure dΩ = 9 ds / (s3 √Disc). Our "typo" suspicion is resolved.
- The r_l signs are corrected.
- Θ(−z/w) appears in N_i, and eq:Rcf has E+z.
- New: d_F R_2 = −z R_1 / D.
- New HAE prefactor 2X/(3A²).
- F^cf_{0,g} carries 3^{g−1}.
- New TR-free lem:cfpole: pole order ≤ 2g−2, and only all-colour-0 graphs reach it.
- **This affects our notes.** BS state the s3 and s3² insertion identities only at leading order in ξ, and the s3² one changed. The pathway's Route B must not call them exact.
- Text to update:
  - Pathway: stage 1 citation, the "z⁹ only" line for F2, and Route B.
  - Scoping note: index slip.
  - leaf_structure_proof cites v1.
  - Gloss our "F_{g,0}" as GW_g.

## Common bottleneck

H1 (the code's rows are, up to z-only factors, BS's frame-B R_0i on T_F and formal WKB solutions of eq:pf) appears in closure, H2, the residue lemma, the regular colours and the frame-B skeleton (N5). So does the T_F continuity lemma. Proving these turns most of the "exact-to-z⁹" results into proofs.

## Ranked next steps

1. **H1 plus the T_F continuity lemma.** BS's R_0i is defined on rational tori; show its limit to T_F is a WKB solution of eq:pf with z-coefficients analytic at α_i = E/3 for y ≠ 0, −1/27. Then define the theorem's rows as BS's R, not the stored truncation.
2. **Stage 1, written from the exact closed form:** F2, L1 (simple poles) and Z∞. Z∞ is the only extra input that H2 needs.
3. **Symplecticity, K(z) = 1, to all orders.** Do it via the closure lemma, the orthogonality argument and F2.
4. **Typeset the proved pieces:** closure, (I1)/(I2), the residue lemma and the s-recursion. Have (I1)/(I2) and the Lagrange identity independently re-derived.
5. **Mixed edges E_0j:** port Prop 9.2 using (R).
6. **Verify the frame-B skeleton and the diff independently.** Then do N1 (normalisation port) and N4 (a written proof of E = (D−X)/2).
7. **Gap corollary, last:** Lemma A up to h = g, plus lem:cfpole and the single-vertex coefficient.

Optional: generate cPP_z12 rows to get spare orders at the top of the fits.

=== VERDICTS ===
cfasym: outcome partial; verifier agrees -- Minor understatement: the reporter says the identification of the code row (cPP_z9.m) with BS's R_00 is only float-supported, at L = 1e8, to about 1e-5. My exact BigRat test is stronger. At finite rational L the normalised frame-B code row equals BS's exact closed form Sum_l r_l (vz)^l / prod(1 + iz/E), with v = 3/(E+3L), through z^9 for four weight pairs, and the z -> -z control fails. This is still finite-order evidence, not a proof of the identification, but it now tests the closed form at generic xi, not only its xi -> 0 limit. | Open problem 3 (equality of the 3d formal Laplace expansion with the iterated (s1,s2)-then-s3 one) is weaker than stated. In eq:asympRp2 the Hessian h is block-diagonal and the non-Gaussian exponent splits as W1_ng(s1) + W2_ng(s3). So exp(z/2 h^{mn} d_m d_n) factorises into commuting operators, and evaluating at s1 = s2 = 0 commutes with the s3 operator. The coupling through g(s) is then handled exactly by doing the (s1,s2) contraction first. This is a one-line formal identity, not a real gap. | Easing of open problem 2: N_i in eq:normNLR does not depend on alpha at all. So the continuity step to T_F needs only regularity of the left-hand saddle coefficients (equivalently R_0i) at alpha_i = E/3 for y != 0, -1/27. N_i itself needs no limit. The reporter's lemma is still needed but is simpler than suggested. | Line ranges are slightly off, with no substantive effect. Thm cfasym label is at 5301; eq:Rcf label at 5334; lem:cfpole label at 5407 (lemma about 5400-5407); rmk:cfgap label at 5422 (remark about 5409-5422). The quoted contents are verbatim-accurate.
frameB: outcome partial; verifier none
diff: outcome confirmed; verifier none
closure: outcome confirmed; verifier agrees -- Precision: in my independent run the closure and PF residuals are exact only to s-order 2 overall, because negative s-powers in the high z-orders cost precision. The report's 's^4/s^5' figures refer to its own extraction. Both are finite-order evidence only, and the all-orders status rests on the hand identity. | Minor: the report says the z^1 cancellation is 'because theta d2 = -E XX C1i^2/Del - (27y/Del) d2'. In my derivation it cancels directly between the term -(E/2)(1+C1i^2) c1 XX from re-expressing Th^2 P0 and the term from theta of g = Del c1 d2 = (E/2)(C1i - c1). The two routes are equivalent and the identity quoted is correct, so this is only a matter of wording. | Minor: 'hypergeometric equation multiplied by (1-Del)': in my convention L2b is (Del-1) times minus the hypergeometric equation. This is only a sign or normalisation difference and is harmless. | Not rerun by me: the off-shell flow identity, the claim that [s^-2] Flow_0 = 0, and the XXSHIFT pole values 210 and 49 and the XX=0 value 1400/3. I checked only that the quoted numbers fit the closed form e1 e2 E^2 (XX(0)+2/9). The report itself marks this part as diagnostic only. | The report frames hypothesis H1 (rows satisfy PF, sourced from BS eq:oscintP and eq:pf) correctly as an open identification. The sqD conjugation (Dop contains -1/3 + t/(1+3t)) changes the explicit form of Th. PF in the conjugated form is checked numerically, not derived from BS, so the open-problem item is real and needed.
h2: outcome confirmed; verifier agrees -- Minor, does not affect the conclusion. h2's off-shell controls (XXSHIFT / XXSHIFT1) shift XX while the gauge prefactors stay at c1 = C1cf, so they are not a consistent buildRmatrix section. In buildRmatrix, XX = y c1'/c1 is tied to the section. My consistent control varies the section itself: c1 = -wc(1+Q Del). It also turns on psi-linear terms and constants (e.g. (7,3), Q=1/10: 't_c^0 sigma_2: a1=-1', 't_c^1 sigma_2: a1=7/6 a0=11/18'). So h2's conclusion that the psi-linear vanishing selects the polarisation stands. Its closed-form prefactor C1cf^2 corresponds to its own scheme. In the consistent scheme, the exact psi-coefficient I find is a1 = E Del wc^2 theta(f)/t_c. | Correction to my own script's printout: the 'my_pred_a1' column has the wrong overall sign. I had printed the z-coefficient, not the psi-coefficient. It equals minus the fitted a1 at every order, k = 0..7, in both off-shell runs, e.g. '-7/6' vs 7/6 and '403/486' vs -403/486. So up to that sign the prediction matches exactly. | The all-orders psi-linear argument is simpler than h2's (I2). It needs only Z-infinity, lim_{z->oo} z sigma_0 = E, together with L1 so that the large-z limit makes sense on rational rows. The derivation: Q1 tends to C1i theta A with A = E/(rho u) = -t_c/alpha^2. Then a1 is proportional to Del c1^2 theta(C1i theta A) = theta(wc/c1)-type, and this vanishes because theta t_c = -wc = c1. F2 is not needed for this step. I did not re-check (I1)/(I2) independently; I relied on this route instead. | Not refuted, but unproved, as h2 itself states: Z-infinity to all orders (exact here only through t_c^6 for sigma_0, from the fits); the absence of constants in sigma_2 to all orders (needs symplecticity); and t_c^8 for sigma_2 and t_c^7 for sigma_0 and sigma_1, which can be tested only under the no-linear hypothesis.
residue: outcome confirmed; verifier agrees -- None that change the conclusion. | Minor point on the conifold control: the sigma rows there are constructed from the conifold proof's (6.4) leaf formulas rather than read from raw rows. It is a consistency check of (6.4) plus the residue formula against Mathematica, not a fully row-level test. | Minor point on the off-shell control: under XXSHIFT=1/10 in the C1cf gauge, sigma_2 at s^2 happens to pass the leaf fit while every other order from s^0 on fails. This does not affect the conclusion that (H2) fails off shell.
regular: outcome confirmed; verifier agrees -- None material. Normalisation: I normalise S_k by Sqrt3^-1 sq^-1 and get S_0 = s/100 at (7,3). The report's S_0 = s/300 comes from normalising by Sqrt3^+1 sq^-1. This is a convention-only factor of 3 and has no effect on PF or the form. | Control subtlety, not an error. Perturbing S_3 along s^1 does not show up at z^4, because Op1(c s) = 0 and this is exactly the c(z) freedom; my PF failures start at z^5. The report's S_3 control (PERT=3,2,...) perturbs the s^2 coefficient, which is why it fails from z^4. So its controls are real tests, but one has to perturb off the s^1 direction (or perturb the relative normalisation) to see a failure at the first order. | Not independently checked by me: the hand-derived Lagrange identity th B = (eta^00 - alpha^3/k)(f'g + fg'). Only its consequence, orthogonality to z^9, was checked; the consistency eta^00 = 3 e1 e2 alpha^3 = alpha^3/k does hold. Also unchecked: the large-radius values rho_l = 1, -1/9, 2/81, and the analyticity of C1i/XX on the conifold section with C1i(0) = -1.
