# Lemma A round 2: referee R4 (saved by the main session from the workflow journal, 2026-10-31)

Subagents cannot write report files, so this file holds the agent's structured return value: the fields first, then the full report text (report_text) verbatim.

## confirmed

- Theorem K1 (Xi^(1)_h(j) = 0 for all h >= 2, 2 <= j <= 2h-2). Every step re-derived by hand. Correct modulo Route B Theorems B, C (ELSV, Okounkov eq. (10), Hurwitz genus 0) and Mumford's ch_1 formula. No KP, (6.1) or TR. The induction is genus-uniform; the only base case is H(1) at j = 1. Numerically K1 holds for h <= 22 (ref_rec_A.pl) and h <= 10 (full bivariate REC, ref_k1.pl).
- K1' (12 Xi^(1)_h(2h-1) = (2h-1-1/(2h-1)) X_h(2h-1) + X_{h-1}(2h-1)): follows from H(h) at j = 2h-1; checked for h <= 22; hand value 1/30 at h = 2.
- Prop. M1 (Mumford reduction): re-derived, including the normalisations of the AC kappa_1, delta, string/dilaton for arbitrary classes and chi_m = -B_2m/(2m). With REC-derived B_h it equals the k = 1 component of the Mathematica tables (genus 2-4), Route E's tables (genus 1-6) and the genus 7-9 tables, coefficient by coefficient.
- Prop. M2 and recursion (A): the REC holds as a polynomial identity in (j, y). The referee's bivariate Y_h (h <= 10) reproduces string/dilaton at y = 0 and y = 1 (not inputs), equals the Mathematica r = 2 data (genus 2-4) and Route E's two-leaf tables (genus <= 7), and satisfies (A) as polynomials in j.
- Prop. M3 (boundary-sum identity for 1 <= j <= 2g-1): induction checked; numerically h <= 22.
- Claim E: all algebraic steps re-derived. (V2) verified symbolically in a (author: a = 1..12 only). P_j = 0 verified symbolically in j (N <= 20). deg R_j <= j+1 for j <= 36. E(h,j) = 0 directly for h <= 16.
- Corollary: sigma^h Xi_h(j) lies in span(sigma..sigma^{h-2}) in the range.
- C1 closed form: equals REC + Mumford for h <= 22, and is exactly minus the sigma^(h-1) coefficient of (6.1) for h <= 14.
- M3' at j = 2g (equivalently 12 Xi^(1)_g(2g) = 3 X_{g-1}(2g)): h <= 22.
- KP nullities reproduced: tlr_k01.pl two-leaf 7/4/1 (G=3), 11/6/1 (G=4), 0 (G=5, K01=2), 0 (G=6, K01=1); full 1 (G=3), 0 (G=4). Independently, Route A's all-T a3_kp_nullity.pl gives the same nullities with H/L options, also with another seed and prime.
- (6.1) as an exact polynomial identity in (j, sigma) at genus 2-9, with an independent implementation. The genus-8/9 tables match Theorem B (k = 0), REC + Mumford (k = 1) and lambda_g (k = h).
- N2: all 52 values of irr_test_out.txt reproduced exactly from REC data alone (no Hodge engine).
- N3 (the CY analogue of M3 fails for k >= 1): reproduced on Mathematica data at genus 2, 3, 4.
- N4 (higher Taylor coefficients are not genus-only boundary sums): reproduced from the referee REC, c = 2, 3, 4, g <= 6.
- K2 reduction of the k = 2 component: reproduced exactly on Mathematica one- and two-leaf data at genus 2, 3, 4; the H_2 coefficients 5, -7, -12 re-derived.
- Reformulation W_h(z) = (2h-2)! [l^2h] S1 S2 (1 + (zB/(1-zB))^(2h-1)) and Lemma A(r = 1) <=> W_h(z) - W_h(0) = O(z^(2h-1)).
- No topological recursion and no wolframscript anywhere in R4_geometry.

## issues

### [1]

#### description

The residual null vector of the KP + Theorem H + K1 system at G = 3 is not supported on Xi^(2)_3(4) alone. Route A's all-T script (H=2 H=3 H=4 L=2 L=3 L=4) gives nullity 1 and lists (2), (4), (3,2), (2,2,2) as not forced, and adding Xi_3(3) = 0 does not remove it. The reduced system in tlr_k01.pl leaves the t = 2 unknowns free, although Xi_G(2) = 0 is proved for all G and all tau (Route A, Theorem A; brief s.2(c)). With it imposed (referee/tlr_k01_z2.pl, Z2=1, or a3_kp_nullity.pl 3 Z=2 H=3 H=4 L=3 L=4), G = 3 has nullity 0, both full and two-leaf, and the G = 4 two-leaf system also has nullity 0. So Lemma A for every h <= 20 needs no exact Hodge data.

#### location

Report s.3.3; summary ('Lemma A for h <= 20 then needs exact data only for one genus-3 number, Xi^(2)_3(4)'); claims list, KP item ('At G = 3 one unknown remains: Xi^(2)_3(4)')

#### severity

minor

#### suggested_fix

Restate: Lemma A(h) for h <= 20 follows from Lemma A(2), KP (Theorem D), Route A Theorem A, Theorem B/H, K1 and Theorem F, with no exact data. The constrained two-leaf system (H + K1 + Xi(2) = 0) is rigid for 3 <= G <= 6.

### [2]

#### description

True, but not new as a consequence. Route A's round-1 run (KP + cut-and-join + H=3 + H=4, runs/nullity_kp_variants.out) already gives nullity 0 at G = 4, and the Hurwitz components are proved for all j. K1 is also not needed at G = 4 once Xi(2) = 0 is imposed (tlr_k01_z2.pl 4 full, K01=1 Z2=1: nullity 0). Round-1 results also close G = 3 without data (a3_kp_nullity.pl 3 C Z=2 Z=3 H=4: nullity 0).

#### location

Report s.3.3 ('Lemma A(4) follows from Lemma A(<=3), KP, Theorem H and K1, without exact genus-4 Hodge data')

#### severity

minor

#### suggested_fix

Say where K1 is genuinely needed: at G = 3 on the KP-only route. It gives a route independent of cut-and-join and of the conifold closed form.

### [3]

#### description

lam1_rec_h8.txt shows 'genus 1 DIFFERS' and '1 MISMATCHES'. It was written by a version without the delta_irr = 1/2 term; the current lam1_rec.pl 8 prints 'equal' throughout (re-run by the referee). lam1_rec_h16.txt was written by a version without the comparison blocks.

#### location

lam1_rec_h8.txt, lam1_rec_h16.txt (outputs cited as evidence)

#### severity

minor

#### suggested_fix

Regenerate both outputs with the current script.

### [4]

#### description

delta_bdry.pl evaluates Xi^(1) through the closed form C1, which is only computer-verified, so it is not evidence for K1'. The REC-only check is hbdry_rec.pl (M3 at j = 2g-1). The proof covers K1' anyway.

#### location

Claims list, K1' ('evidence: ... delta_bdry.pl for g <= 9')

#### severity

minor

#### suggested_fix

Cite hbdry_rec.pl, or the referee's ref_rec_A.pl (h <= 22), as the evidence.

### [5]

#### description

This is a heuristic. No cosection-localised virtual class on the log moduli space with leaves is constructed, and nothing is derived from it.

#### location

Claims list, 'Cosection observation' (status proved); report s.5

#### severity

minor

#### suggested_fix

Change the status to 'observation (heuristic)'.

### [6]

#### description

N2-N4 are computer refutations of specific ansaetze in finite ranges, not proofs that no extension exists. For example, N4 tests only boundary sums whose constants depend only on the genus, not sums with psi-insertions at the node, which the author's own next steps propose.

#### location

Summary ('Several mechanisms provably do not extend'); report s.4, s.8

#### severity

minor

#### suggested_fix

Write 'the following ansaetze fail (computer-verified, ranges ...)'.

### [7]

#### description

phi_a is undefined, and the inference 'a pure localization argument at r = 1 can only yield the even j' is not explained. The outputs show E^(1) = NE^(1) = 0 throughout the range at genus 3 and 4, which is the relevant fact.

#### location

Report s.4 item 2; claims list, N1 (k = 1 formula)

#### severity

cosmetic

#### suggested_fix

Define phi_a or drop the formula, and state the k = 1 vanishing of both parts.

### [8]

#### description

The right side of H_2 is not written out. It also contains genus-h known data (Xi^(1)_h(0), Taylor coefficients of X_h, B_h = BS_h), not only lower-genus data. The coefficients 5, -7, -12 re-derive correctly.

#### location

Report s.4 item 6 (identity H_2)

#### severity

cosmetic

#### suggested_fix

Write H_2 out explicitly; k2_reduce.pl effectively contains it.

### [9]

#### description

'Q_h(j+1) - Q_h(j) = Sum_{m>=1} b_m Q_{h-m}(j), where Sum_m b_m mu^{2m} = B' should read B - 1 = Sum_{m>=1} b_m mu^{2m}.

#### location

Report s.6, Sheffer form

#### severity

cosmetic

#### suggested_fix

Replace B by B - 1.

### [10]

#### description

The test (irr_test.pl) subtracts only the irreducible-boundary part. Equivalently, the l >= 3 analogue of M3 (kappa_l - psi^l + separating part = 0) fails.

#### location

Summary ('Single ch_l insertions (l >= 3) do not localize to the boundary')

#### severity

cosmetic

#### suggested_fix

Say 'differ from their irreducible-boundary part (the l >= 3 analogue of M3 fails)'.

### [11]

#### description

Three small points. (1) The restriction j + y != 2h is unnecessary, since Lemma 6.2 holds at delta = 0. (2) The headers of claimE.pl, delta_bdry.pl and hbdry_rec.pl cite 'report s.4' for the K1 induction, which is in s.2. (3) 'TLR holds for 5 <= G <= 20' redefines Route A's unconstrained TLR(G).

#### location

Prop. M2 step 1; script headers; report s.3.3

#### severity

cosmetic

#### suggested_fix

Drop the restriction, fix the section references, and give the constrained rigidity statement its own name.


## notes

Headline: Theorem K1 is correct. It is the first all-genus vanishing of a CY component beyond Hurwitz and lambda_g. It leaves the components 2 <= k <= h-1 open, as the author says.

Main improvement for the main session (report s.3): with the proved Xi_G(2) = 0 added, genus 3 closes, so Lemma A for every h <= 20 needs no exact Hodge data. The route is: Lemma A(2), KP Theorem D, Route A Theorem A, Theorem B/H, K1 at G = 3, then Theorem F. Round-1 results alone also close genus 3 if cut-and-join and the even-level conifold argument are used (a3_kp_nullity.pl 3 C Z=2 Z=3 H=4: nullity 0).

Not re-checked independently:
- N1, the excess/non-excess split, which needs a Hodge engine for K-twisted moments; its internal E + NE = table check passes.
- The unit-leg dressing reformulation.
- The genus-8/9 middle components, except through (6.1) itself.

Runtimes: data_components.pl about 29 min (genus-4 r = 2 Mathematica table); ref_k1.pl 10 about 6 min; ref_rec_A.pl 22 about 19 min; ref_claimE.pl 40 20 36 16 about 1.5 h. Referee files are in G:\My Drive\Projects\Old\refGW\hodge_vertex_conifold\lemmaA\round2\R4_geometry\referee\ (outputs in out\). Nothing outside that folder was written.

## refuted

- 'At G = 3 one unknown remains: Xi^(2)_3(4) ... the only exact-data input left for Lemma A with h <= 20.' The residual null vector of the KP + H + K1 system at G = 3 also involves T = (2); Xi_G(2) = 0 is already proved (Route A, Theorem A). Imposing it gives nullity 0 at G = 3 (full and two-leaf; tlr_k01_z2.pl with Z2=1, and a3_kp_nullity.pl 3 Z=2 H=3 H=4 L=3 L=4). So Lemma A for every h <= 20 needs no exact Hodge data.

## route

R4 (geometry): intersection theory, localization, Givental. Main claim: Theorem K1, the lambda_1-component of the one-leaf case in all genera. Referee pass, round 2.

## verdict

sound

## report_text

# Referee report on R4 (geometry), round 2 of the attack on Lemma A

**Folders.**
- Referee scripts: `G:\My Drive\Projects\Old\refGW\hodge_vertex_conifold\lemmaA\round2\R4_geometry\referee\`
- Referee outputs: the `out\` subfolder of that folder.

**How the checks were run.**
- All arithmetic is exact (Perl Math::BigRat).
- The referee code was written from scratch (RL.pm, RParse.pm). It shares nothing with PJ.pm, RefHodge.pm or the author's scripts.
- No wolframscript was run. No topological recursion was used in any form. At most 4 heavy processes ran at a time.
- Files were written only in `referee\`. The author's scripts and Route A's scripts were only run; they write nothing.

## 0. Verdict: SOUND

**Theorem K1 is correct.** I re-derived every step by hand:
- Prop. M1;
- the polynomial REC and recursion (A);
- the induction for M3;
- Claim E, including its generating-function proof;
- the conclusion, and K1'.

I found no gap. K1 holds modulo Route B's Theorems B and C (ELSV, Okounkov's Toda eq. (10), Hurwitz's genus-0 formula; the Route B referee found these sound) and Mumford's ch_1 formula. It uses no KP, no (6.1) and no TR. The induction is uniform in the genus, and its only base case is H(1) at j = 1.

**Every positive computational claim I tested survived.** Most were re-checked with independent code and independent data:
- the Mathematica q tables in `lemmaA\data` (genus 2-4, r = 1, 2), extracted with my own parser;
- Route E's engine tables (genus up to 6, and 7);
- the genus-7 tables;
- the author's new genus-8 and genus-9 tables.

**One consequence is stated too weakly.** The report says a single genus-3 number, Xi^(2)_3(4), is the only exact-data input left for Lemma A with h <= 20. That is not so. The residual null vector at G = 3 also has a component on T = (2), and Xi_G(2) = 0 is already proved for all G (Route A, Theorem A). Once it is imposed, genus 3 is rigid. So **Lemma A for every h <= 20 needs no exact Hodge data at all** (s.3).

**The other issues are minor or cosmetic.** They concern stale output files, an undefined symbol, status labels and wording.

**Not proved, as the author says:** (OL), (6.1) and Lemma A in all genera. The components 2 <= k <= h-1 of the one-leaf case remain open.

## 1. Theorem K1, line by line

### 1.1 Normalisations

- **Components.** (tau - a)(-1 - tau - a) = -sigma + a + a^2. So the CY class splits as Lambda^vee(1) Sum_k (-1/sigma)^k e_k(alpha + alpha^2).
- **The lambda_1-component.** e_1(alpha + alpha^2) = lambda_1 + p_2(alpha) = lambda_1 + 2 ch_2(E) = lambda_1, because ch_2(E) = 0 (Mumford). So Xi^(1) is indeed Int c(T(-log)) lambda_1 lam_j.
- **Forgetting a point.** Arbarello-Cornalba give kappa_k = pi^* kappa_k + psi_{n+1}^k. Hence K_{n+1} = pi^* K_n (1 - psi_{n+1}), and c(T_{h,n+1}(-log)) = pi^* c(T_{h,n}(-log)) (1 - psi_{n+1}).
- **Gluing maps.** Lambda^vee(1) and K restrict to products of vertex classes, with the node branches marked. This is consistent with the residue sequence for T(-log).
- **String and dilaton.** Both hold for an arbitrary class on Mbar_{h,n}, not only for classes pulled back from Mbar_h. They give:
  - Y_h(j,1) = j X_h;
  - Y_h(j,0) = (j - 2h + 1) X_h;
  - Xi^(0)_{h-1}(j,0,0) = (j - 2h + 2)(j - 2h + 3) X_{h-1}(j).
- **Mumford on Mbar_{h,1}.**
  - ch_1 = (B_2/2!) [kappa_1 - psi + (1/2) Sum_xi xi_* 1] = (kappa_1 - psi + delta)/12, where kappa_1 is the Arbarello-Cornalba class pi_*(psi_2^2).
  - The ordered separating splittings count each divisor Delta_{g1,{1}} twice, so with the 1/2 each appears once. Delta_irr is (1/2) xi_* 1 from Mbar_{h-1,3}.
- **chi_m.** chi_m = X_m(0) = (2m-2)! [t^2m] S^2 = -B_2m/(2m) = chi(M_{m,1}) (Harer-Zagier). Check: X_2(0) = 1/120 agrees with the Mathematica data.
- **zeta_m.** zeta_m = T_m'(0) = (2m)! [t^2m] log H = (2m-1) B_2m/(2m). This uses log S = -Sum_m B_2m t^2m/(2m (2m)!).

### 1.2 Prop. M1 (Mumford reduction)

**Steps 1-6 are correct.** The kappa_1 term is the only delicate one:
- pi^* psi_1^a times psi_2^2 equals psi_1^a psi_2^2, because psi_2 D_12 = 0;
- pi^* c(T_{h,1}(-log)) = c(T_{h,2}(-log)) lam_1(psi_2).

Together these give p_1 = Sum_{c>=2} [y^c] Y_h = (2h - 1) X_h - B_h. For h = 1 the formula gives 1/24 - 1/24 + 1/2 = 1/2 = 12 Int_{Mbar_{1,1}} lambda_1.

**Independent confirmation.** I took B_h from my own REC (s.1.3) and X_h from Theorem B. The right side of M1 then equals the k = 1 component, coefficient by coefficient, in three sets of tables:
- Mathematica, genus 2-4;
- Route E's engine, genus 1-6;
- the genus 7-9 tables.

**Hand values.**
- Xi^(1)_2(2) = 0.
- Xi^(1)_2(3) = 1/360, so 12 Xi^(1)_2(3) = 1/30. This matches K1': (8/3)(-1/320) + 1/24.

### 1.3 Prop. M2 and recursion (A)

**The two-leaf REC is reproduced.** Route B's Theorem C at r = 2 has:
- e = 2h and delta = j + y - 2h;
- one-block terms with h'' = h + 1 - k;
- two-block terms with h1 + k1 + h2 + k2 = h + 2.

**Polynomiality.**
- The genus-0 terms become polynomials: (j + y)_(2h+2) Y_0 = (j + y - 1)_(2h+1), and (x + 1)_(2k) X_0 = (x - 1)_(2k-2).
- The lattice {j, y >= 2} is Zariski dense, so the identity holds as polynomials.
- The author's restriction j + y != 2h is not needed: Lemma 6.2 also holds at delta = 0. It is harmless.

**Two-block simplification.** (j - 2h1 + 1)_(2k1) (j - 1)_(2h1-2) = (j - 1)_(2(h1+k1-1)), and Sum_k c_k t^(2k-2) = S^-2. Together these give Sum T_m1(j) T_m2(y).

**Derivation of (A).** Differentiate at y = 0, using:
- B_0 = -1/j^2;
- Y_g(j,0) = T_g(j)/j, which follows from S^2 H^j = H^j - (t/j)(H^j)' and t (log H)' = 1 - S^2.

**Referee computation (`ref_k1.pl`).** I computed the full bivariate polynomial Y_h(j,y) from the REC by exact division by (j + y - 2h), for h <= 10. Findings:
1. The two-block double sum equals Sum T_m1(j) T_m2(y) exactly.
2. The REC polynomial reproduces the geometric values at y = 0 and y = 1 (dilaton and string). These values are not inputs: REC is derived only for parts >= 2. So this confirms the polynomial continuation that Theorem K1 relies on.
3. (A) holds as a polynomial identity in j for h <= 10.
4. Y_h equals the Mathematica r = 2 data (k = 0) at genus 2, 3, 4 (21, 45 and 78 coefficients) and Route E's two-leaf tables for genus <= 7.

### 1.4 Prop. M3: the induction

In (A), at an integer 1 <= j <= 2h-1:
- (j - 2h + 1) X_h(j) = 0, by the zeros j = 1..2h-2 of Theorem B and the factor at j = 2h-1;
- (j - 2h'')_(2h+2-2h'') vanishes for 2h'' <= j <= 2h+1. So only the h'' with j <= 2h''-1 survive, and for those H(h'') applies (h'' < h, same integer j);
- T_{h''}(j) = 0 for 1 <= j <= 2h'';
- j - 2h != 0.

So H(h) at j is equivalent to E(h,j) = 0. The base case is H(1) at j = 1, with B_1(1) = -1/12. No hypothesis on other components, and none on Lemma A in lower genera, is used.

Remark: at j = 2h the coefficient (j - 2h) vanishes, so (A) does not determine B_h(2h). M3' is therefore a genuinely separate, computer-verified observation, not a consequence of this induction.

### 1.5 Claim E

**All steps re-derived:**
- The chi_M coefficient: (2M + 2n + 1)(a - 2M) + (a + 2n)(2M - 1) = 2M(2a - 2M - 1) + 2n(a - 1).
- Extending the m-sum to m >= 1. This needs 2M >= a + 1. Every n with T_n(j) != 0 has M = h - n and 2M >= a + 1, which is equivalent to j <= 2h-1.
- Extending the n-sum. This uses T_n(j) = 0 for 1 <= j <= 2n.
- The generating function. It uses (2a-1)D - D^2 + (a-D)(a-1-D) = a(a-1), and the log terms cancel by (V1).
- P_j. Its w^2N coefficient is 2 (j-1)_(2N-1) [...], by (j-1)_(2n) (j-2n)_(2N-2n) = (j-1)_(2N-1)(j - 2n). Its vanishing follows from (cosh t - 1)(t G' - j(1 - S^2) G) and t G'/G = j(1 - S^2).
- R_j has degree <= j + 1 for integer j, so [w^(2h+2)] R_j = 0 for j <= 2h.

(V1) is the asymptotic expansion of psi(x+1) - psi(x) = 1/x at x = 1/w. Equivalently, it is a formal Bernoulli identity.

**Referee checks (`ref_claimE.pl`):**
- (V1) to order w^40.
- (V2) as an identity of polynomials in a, to w^40. The author checked only a = 1..12.
- P_j = 0 as polynomials in j, for N <= 20.
- deg R_j <= j + 1 for j <= 36.
- E(h,j) = 0 directly for h <= 16 (256 pairs).

### 1.6 Conclusion of K1, and K1'

**The final identity checks out.** B_h - BS_h = (2h - 1 - 1/j) X_h - 12 Xi^(1)_h + (1/2)(j - 2h + 2)(j - 2h + 3) X_{h-1}, with BS_h = chi_h/j + (the separating sum). In the range:
- X_h(j) = 0;
- X_{h-1}(j) = 0 for j <= 2h-4, and the prefactor kills j = 2h-3 and j = 2h-2.

So K1 follows. At j = 2h-1 the same identity gives K1'.

**Independent numerical confirmation:**
- `ref_rec_A.pl` uses recursion (A) in the T/zeta form, a different code path from the author's explicit double sum. It confirms K1, K1', M3 at j = 1..2h-1 and M3' at j = 2h for every h <= 22.
- `ref_k1.pl` confirms the same from the full bivariate REC for h <= 10.

## 2. The other claims

**C1 (closed form of the lambda_1-component).**
- It equals REC + Mumford for h <= 22; the author went to 16.
- It equals, exactly, minus the sigma^(h-1) coefficient of the right side of (6.1) for h <= 14 (`ref_c61_c1.pl`, symbolic in j and sigma). So the eps-expansion of s.3.1 is right.
- For j >= 2h it remains conjectural beyond these ranges. That is the correct status.

**M3' (j = 2g), equivalently 12 Xi^(1)_g(2g) = 3 X_{g-1}(2g).** Holds for h <= 22.

**(6.1) at genus 7, 8, 9.**
- Confirmed with a separate implementation (`ref_c61.pl`), as exact polynomial identities in (j, sigma) for genus 2-9. It uses p_n = e1^n + e2^n, with p_n = p_{n-1} + sigma p_{n-2}.
- The new genus-8 and genus-9 tables come from one engine only. Three of their components are checked independently: k = 0 equals Theorem B, k = 1 equals REC + Mumford, and k = h equals b_h (j-1)...(j-2h+2).
- The middle components are cross-validated only through (6.1) itself, which an engine error would not reproduce.
- The author's comparison at sample points already proves the identity: it uses 3h+2 values of j for degree <= 3h-2, and h+3 distinct values of sigma for degree <= h.

**Corollary.** In the range, sigma^h Xi_h(j) lies in span(sigma, ..., sigma^(h-2)). Correct, given K1, Theorem B or H, and the lambda_g formula.

**KP claims, reproduced with the author's `tlr_k01.pl`:**

| G | two-leaf nullity: none / +H / +H+K1 | full system nullity |
|---|---|---|
| 3 | 7 / 4 / 1 | 1 |
| 4 | 11 / 6 / 1 | 0 |
| 5 | (+H+K1) 0 | – |
| 6 | (+H) 0 | – |

- **Independent cross-check.** Route A's all-T script `a3_kp_nullity.pl` uses a different formulation, imposing the constraints through its H= and L= options. It gives:
  - G = 3, KP + H + K1: nullity 1;
  - G = 4, KP + H: nullity 1;
  - G = 4, KP + H + K1: nullity 0.

  The same results come out with another seed and the prime 10^9+7 (`kp_null_seedprime.pl`).
- **The constraints are legitimate.** They act only on the one-leaf unknowns c_{t,G} and c_{t,G-1} for t = 2..2G-2, which is exactly K1's range.
- **Lifting.** The modular lifting argument is the one the Route A referee accepted.

**K2 (exact reduction of the k = 2 component).**
- Reproduced exactly at every j on the Mathematica one- and two-leaf data for genus 2, 3, 4. The three-leaf genus-(h-1) Hurwitz input, used only for the irreducible boundary of ch_3, was taken from Route E's tables.
- The coefficients 5, -7, -12 in H_2 re-derive as 1440 x (1/288, 1/288 - 1/120, -1/120).

**N2 (single ch_l insertions).** All 52 values of `irr_test_out.txt` (l = 3 at genus 3-6, l = 5 at genus 5 and 6) were reproduced exactly from REC data alone: my bivariate Y_h plus Theorem B, with no Hodge engine. At genus 2, the l = 3 quantity vanishes on j = 2..4.

**N3 (CY analogue of M3).** Reproduced on the Mathematica data at genus 2, 3, 4:
- for k >= 1, equality holds only at j = 1;
- for k = 0 it holds at j = 1..2g, which is M3 together with M3'.

**N4 (higher Taylor coefficients).** Reproduced from my REC for c = 2, 3, 4 and g <= 6. For c = 1 the fitted constants are chi_m, consistent with M3.

**N1 (excess part at r = 1).** Not recomputed independently: it needs K-twisted moments of Hodge monomials, i.e. a Hodge engine. The script's internal check that E + NE equals the table value passes at every point. The outputs also show E^(1) = NE^(1) = 0 throughout the range at genus 3 and 4. The text should say this plainly; the formulas with phi_a are not interpretable (phi_a is undefined).

**Reformulations (s.6).**
- The W_h(z) form, and Lemma A(r = 1) <=> W_h(z) - W_h(0) = O(z^(2h-1)), re-derive correctly. They use Xi_h(1) = 0, which follows from the string equation.
- The Sheffer recursion needs B - 1 = Sum_{m>=1} b_m mu^2m.
- The unit-leg dressing was not checked.

**TR.** A grep of the R4 folder finds no spectral curve, Eynard-Orantin, loop equation, BKMP, remodelling or TR script, and no wolframscript.

**Re-runs of the author's own scripts.** Each reproduced its stored output: `claimE.pl 12`, `rec2_test.pl 4`, `p1moments.pl 3`, `lam1_rec.pl 8`, `tlr_k01.pl` (G = 3, 4, 5, 6). The one exception is the stale `lam1_rec_h8.txt` (issue 3).

## 3. Genus 3 and Lemma A for h <= 20: no exact data needed

**The author's G = 3 statement mis-describes the null vector.** The report says that, after Theorem H and K1, the genus-3 KP system has one unknown left, Xi^(2)_3(4), and that it is the only exact-data input for h <= 20.
- `a3_kp_nullity.pl 3 H=2 H=3 H=4 L=2 L=3 L=4` gives nullity 1, with (2), (4), (3,2), (2,2,2) not forced. So the null vector is not supported on T = (4) alone.
- Adding Xi_3(3) = 0 does not remove it either (`a3_kp_nullity.pl 3 Z=3 H=2 H=4 L=2 L=4`: nullity 1).
- The reduced system in `tlr_k01.pl` leaves the t = 2 unknowns free, although Xi_G(2) = 0 is proved for every G and every tau (Route A, Theorem A; brief s.2(c)).

**With Xi_G(2) = 0 added the systems become rigid.** I added it to a copy, `referee\tlr_k01_z2.pl`, run with environment variable Z2=1:

| G | two-leaf, +H+K1+Xi(2) | full, +H+K1+Xi(2) | full, +H+Xi(2) (no K1) |
|---|---|---|---|
| 3 | 0 | 0 | 1 |
| 4 | 0 | 0 | 0 |

Route A's script agrees: `a3_kp_nullity.pl 3 Z=2 H=3 H=4 L=3 L=4` gives nullity 0.

**Consequences:**
- Lemma A for every h <= 20 follows from:
  - Lemma A(2) (by hand, or from Route A's Theorem A, since (2) is the only genus-2 point);
  - KP (Theorem D);
  - Route A's Theorem A;
  - Theorem B or H;
  - K1;
  - Theorem F.

  No exact Hodge data is needed.
- K1 is needed only at G = 3 along this route. At G = 4, H + Xi(2) = 0 already gives nullity 0.
- With these proved constraints, the two-leaf system is rigid for every 3 <= G <= 6. Route A's unconstrained TLR covers 7 <= G <= 20.
- Round-1 results alone also close genus 3 and 4 without data:
  - `a3_kp_nullity.pl 3 C Z=2 Z=3 H=4` gives nullity 0, using the cut-and-join relations (Route A, Theorem B) and Xi_3(3) = 0 from the even-level conifold argument (brief s.2(f));
  - Route A's round-1 run `KP + C + H=3 + H=4` at G = 4 also gives nullity 0.

  So "Lemma A(4) without genus-4 data" is not new as a consequence. K1 supplies an independent route that needs neither cut-and-join nor the conifold closed form.

## 4. Issues

1. **(minor)** The G = 3 "one unknown Xi^(2)_3(4)" and "only exact-data input" statements are inaccurate (s.3). Restate as above.
2. **(minor)** "Lemma A(4) without exact genus-4 data" is already available from round 1, and K1 is not needed at G = 4 once Xi(2) = 0 is imposed. Say where K1 is actually needed, namely G = 3 on the KP-only route.
3. **(minor)** Stale outputs.
   - `lam1_rec_h8.txt` shows genus-1 "DIFFERS" and "1 MISMATCHES". It was written by a version without the delta_irr = 1/2 term; the current script prints "equal".
   - `lam1_rec_h16.txt` was written by a version without the comparison blocks.

   Regenerate both.
4. **(minor)** The evidence quoted for K1' includes `delta_bdry.pl`, which uses the closed form C1 (only computer-verified). The REC-only check is `hbdry_rec.pl` at j = 2g-1. The proof covers K1' in any case.
5. **(minor)** The cosection observation is labelled "proved". It is a heuristic: no cosection-localised virtual class is constructed (with the log structure and the leaves), and nothing is derived from it. Label it as an observation.
6. **(minor)** "Several mechanisms provably do not extend" overstates N2-N4. They are computer refutations of specific ansätze in finite ranges. For example, N4 tests only boundary sums whose constants depend only on the genus, not sums with psi-insertions at the node, which the author's own next steps propose.
7. **(cosmetic)** N1: phi_a is undefined, and "can only yield the even j" is not explained. State that E^(1) = NE^(1) = 0 in the range at genus 3 and 4.
8. **(cosmetic)** The right side of H_2 is not written out. It also contains genus-h known data (Xi^(1)_h(0), Taylor coefficients of X_h, B_h = BS_h), not only lower-genus data.
9. **(cosmetic)** In the Sheffer form, Sum_{m>=1} b_m mu^2m equals B - 1, not B.
10. **(cosmetic)** The summary's "do not localize to the boundary" should say "differ from their irreducible-boundary part", i.e. the l >= 3 analogue of M3 fails.
11. **(cosmetic)** Housekeeping:
    - M2's restriction j + y != 2h is unnecessary.
    - The script headers cite "report s.4" for the K1 induction, which is in s.2.
    - "TLR(G) holds for 5 <= G <= 20" redefines Route A's unconstrained TLR(G). Give the constrained statement its own name.

## 5. Obstructions and simpler paths

**The obstructions are real as stated.**
- (OL) needs every component 2 <= k <= G-1. After K1 that is G-2 numbers per (G, D).
- N2-N4 show that the mechanism of K1 does not extend term by term:
  - single ch_3 and ch_5 insertions do not localise;
  - [y^c] Y for c >= 2 are not genus-only boundary sums;
  - the CY analogue of M3 fails.
- They do not exclude:
  - boundary sums with psi-insertions at the node;
  - or an argument that uses the CY combination e_k(alpha + alpha^2) as a whole.

  N2 indicates that only such a combination can work.

**No simpler proof of K1 is visible.** M3 is equivalent to "p_1 (the first Lambert coefficient of the Hurwitz class) is a boundary sum in the low range". That is the next-to-leading order of Route B's T -> oo analysis, which the author in effect carries out through the two-leaf REC at y -> 0.

**The simplification found is the one in s.3:** impose the proved Xi_G(2) = 0.

**For an all-genus attack, two suggestions.**
- Push the constrained two-leaf rigidity (H + K1 + Xi(2) = 0, rigid for 3 <= G <= 6) toward a structural proof of TLR.
- Extend K1 to r >= 2: REC with a zero leaf, plus Mumford.

## 6. Referee files (all in `referee\`)

| File | What it does | Output (`out\`) |
|---|---|---|
| `RL.pm` | Library: univariate and bivariate polynomials, power series, Bernoulli numbers | – |
| `RParse.pm` | Evaluator for Mathematica InputForm expressions | – |
| `data_components.pl` | Mathematica q tables to sigma-components | `data_h{2,3,4}_r{1,2}.txt` |
| `cmp_tables.pl` | Coefficient-by-coefficient table comparison | – |
| `ref_k1.pl` | Full bivariate REC, M1, C1, K1, M3 | `ref_k1_G10.txt`, `xi1_ref.txt`, `Y_ref_h*.txt`, `X_ref_h*.txt` |
| `ref_rec_A.pl` | Recursion (A) in T/zeta form, h <= 22 | `ref_rec_A_22.txt` |
| `ref_claimE.pl` | Symbolic checks of Claim E | `ref_claimE_big.txt` |
| `ref_c61.pl` | (6.1) against the tables | `ref_c61_out.txt` |
| `ref_c61_c1.pl` | C1 as the 1/sigma part of (6.1) | `ref_c61_c1_out.txt` |
| `ref_n2_n4.pl` | N2 and N4 from REC data | `ref_n2_n4_out.txt` |
| `ref_n3_k2.pl` | N3 and K2 on Mathematica data | `ref_n3_k2_out.txt` |
| `kp_null_seedprime.pl` | Route A's script with another seed and prime | `kp_runs_out.txt` |
| `tlr_k01_z2.pl` | Author's KP script plus Xi(2) = 0 | `kp_runs_out.txt` |
| (author's `tlr_k01.pl`) | Reproduction runs | `tlr_g5_two_K2.txt`, `tlr_g6_two_K1.txt` |

