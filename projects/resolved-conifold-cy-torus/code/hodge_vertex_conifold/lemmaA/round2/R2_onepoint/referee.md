# Lemma A round 2: referee R2 (saved by the main session from the workflow journal, 2026-10-31)

Subagents cannot write report files, so this file holds the agent's structured return value: the fields first, then the full report text (report_text) verbatim.

## confirmed

- Theorem 1 (closed n-point formula on C_tau): proved-modulo-literature is justified. The statement and hypotheses of ABDKS arXiv:2312.16950v3 Prop. 3.1 / eq. (17), with Rem. 2.6, 3.2, 3.3, 3.4 and the standing assumptions of Thm 2.8, were checked against the paper, and the derivation steps (i)-(iv) were re-done. The referee's formal implementation (any rational tau) matches the Route E engine at 3058 values (n = 1..7, g = 1..7, integer and non-integer tau, 0 mismatches mod p) and the exact EO recursion without P2: 403 values at tau = 1, 316 at tau = -1/2.
- Theorem 2 (degree bound = Lemma D): the homogeneity argument is correct, and in the formal form it holds for every tau != 0, -1, not only integers. Range zeros were checked to genus 12 (n = 1) and genus 9 (n = 2).
- Theorem 3 (dictionary P2 = Lemma B) for all g, n and tau with sigma != 0. ABDKS (40)-(50) were checked against the source. The re-derivation covers p = 1/tau, q = -1/e2, r = 1, the gamma -> 0 limit, q_k = (-1)^k kappa, s_i = (-1)^{i+1}/i, the degree involution and the factor sigma^{1-g}.
- Moment lemma: Res_{v=tau} xi_k dv = 1 for all k (proof re-done; moments_check.pl re-run).
- Corollary 4 (Lemma A = M_irr(h) for all h, so (OL) for all G without inductive hypotheses), for non-empty R: correct given Theorems 2 and 3.
- Theorem 5 (divided-difference evaluation): the elementary proof was re-done. Formal Theorem 1 at n = 1 agrees with (6.1) at 3028 values (h <= 12, 13 framings).
- Corollary 6 ((6.1) for all h >= 1, all j, all tau): correct. (6.1) agrees with the engine at 1274 values (h <= 7, including the new genus-7 table). Route E's Zhat_j(u) form follows.
- Corollary 7 (gap formula and the string/dilaton consequences): correct. 154 checks against (6.1) at j = 0 up to h = 12; the sigma -> oo limit gives the Harer-Zagier chi(M_h).
- Two-point difference form: re-derived independently from Theorem 1 at n = 2.
- F1 (Murnaghan-Nakayama, coefficient of p_j): correct.
- Computer verification: confirmed, with the corrections in issues 1 and 2. Independent extra checks also agree: 358 values against the Mathematica q-tables, 57 against Route A's Marino-Vafa (Schur-side) values, 24 first-layer values against the main session's Mathematica EO at tau = 2, 3, and 840 genus-0 values.

## issues

### [1]

#### description

V11 is reported as 21 engine comparisons at (g,n) = (2,5), tau = 1, with 0 mismatches. No engine table E_formulas/engine_h2_r5.txt exists, and all 21 lines of runs_np_g2_n5_t1.out read '(no data)'. So the engine total is 1051, not 1072, and the author compared n = 5 only with the exact EO recursion at tau = 1 (26 values, V13).

#### location

R2 report s.10 row V11; R2_onepoint/runs_np_g2_n5_t1.out; the totals 1072 engine values and n <= 5

#### severity

minor

#### suggested_fix

Correct the V11 row and the totals. The referee's runs now give engine comparisons at n = 5 (g = 4, 5; tau = 2, 1/2, 3) and at n = 6, 7 (g = 5), all with 0 mismatches.

### [2]

#### description

Tuples with no table data ('(no data)') count towards the total but never as mismatches. A run with a missing table therefore prints 'total N, mismatches 0'. Reproduced: perl xy_npoint.pl 1 7 1 21 without the table argument gives 21 '(no data)' lines and 'total 21, mismatches 0'. The bare summary lines (e.g. the genus-7 rows) do not prove that a comparison took place.

#### location

R2_onepoint/xy_npoint.pl and the scripts built the same way; the summary lines in runs_all_checks.txt

#### severity

minor

#### suggested_fix

Print the compared and no-data counts separately, and fail when the table is missing. With the table, the genus-7 run at tau = 1 gives 21 genuine OK lines, and the referee's sweep covers genus 7 at 9 framings.

### [3]

#### description

At h = 0 the coefficient of lambda^0 on the right side is (j-1)!/(j+1)! = 1/(j(j+1)). It is not polynomial in j and is undefined at j = 0. With the brief's definition Xi_0(j) = 0, because every term vanishes by degree. The identity is right for h >= 1 only.

#### location

Corollary 6, generating-function form (sum over h >= 0), and the matching key formula

#### severity

cosmetic

#### suggested_fix

Sum over h >= 1 (subtract the 1/(j(j+1)) term), or state that the h = 0 term is a convention.

### [4]

#### description

The wording includes R = empty, but Xi_h(empty) = Int K Omega is non-zero by Corollary 7. The proof uses Theorems 2 and 3 with n = r >= 1, so it covers non-empty R only.

#### location

Corollary 4, the summary and all_genus_statement: 'every multiset R of positive integers with D(R) <= 2h-3'

#### severity

cosmetic

#### suggested_fix

Write 'every non-empty multiset R'.

### [5]

#### description

Integrality is not needed. In the formal form Q_mu(v;c)/v^mu = (1+cw)^{e2 mu} exp(Sum_k (2k-1)! hb^{2k}[phi_k(tau mu) tau^{-2k} w^{2k} - phi_k(e2 mu) e2^{-2k} w^{2k}(1+cw)^{-2k}]), every power of w = 1/v carries one factor of c or hb. So R_mu(c,hb) is a homogeneous polynomial of degree |mu| for every tau != 0, -1. The referee's formal implementation confirms the formula at 7 non-integer framings, against the engine and against the exact EO recursion at tau = -1/2.

#### location

Theorem 2 (hypothesis 'integer tau >= 1') and the polynomial-in-tau step of Corollary 4

#### severity

cosmetic

#### suggested_fix

State Theorems 1 and 2 for every tau with sigma != 0, in the formal form. Keep the finite-product form as the integer-tau special case. The polynomiality step in Corollary 4 then becomes optional.

### [6]

#### description

The range zeros and the absence of hb^p with p > |mu| follow from the shape of (17) alone. The referee's negative control (S in place of Sc in the dressing) fails 14 of 24 out-of-range values but still gives every range zero. V20 therefore says nothing about whether Theorem 1 is correct; the out-of-range agreements are what validate it.

#### location

R2 report s.10, row V20 (the homogeneity check), counted as evidence

#### severity

cosmetic

#### suggested_fix

Describe V20 as a code-consistency test, and point to the out-of-range comparisons as the validation of Theorem 1.

### [7]

#### description

The report labels Eynard's CNTP 8 (2014) paper [Eyn14a]. In ABDKS's reference list that paper is [Eyn14], and [Eyn14a] is Eynard's ICM 2014 overview. ABDKS attribute the kappa-insertion formula (50) to [Eyn14], and the plain triple-Hodge formula (44) to Bouchard-Marino, Zhou and Chen.

#### location

R2 report s.12 (literature list)

#### severity

cosmetic

#### suggested_fix

Align the labels with ABDKS so that a reader following ABDKS's numbering finds the right paper.


## notes

The all-genus result depends on TR for the Hodge CohFT curve C_tau (ABDKS log-TR x-y swap, Prop. 3.1; Eynard's kappa-class formula as in ABDKS (50)). Brief s.0 allows this and the user confirmed it as an allowed input. Nothing applies TR to the resolved conifold x C^2 theory. All referee comparisons are modulo p = 2^31-1, with 2^31-19 on subsets; a mod-p mismatch is certain, an accidental match has probability about 1/p. One bug in the referee's own script (sigma^{-1} at g = 0) was found and fixed; it affected only the genus-0 run. The proofs of ABDKS Thm 2.8 and of Eynard 2014 were not re-checked, only their statements and hypotheses. Scripts and outputs: G:\My Drive\Projects\Old\refGW\hodge_vertex_conifold\lemmaA\round2\R2_onepoint\referee\ (main implementation ref_xy_mod.pl; sweeps sweep1.sh and sweep1b.sh).

## refuted


## route

R2_onepoint (round 2): one-leaf function, (6.1) and Lemma A in all genera via the x-y swap on the Hodge curve C_tau

## verdict

sound

## report_text

# Referee report: R2, the one-leaf function, (6.1) and Lemma A via the x-y swap on C_tau (round 2)

Referee scripts and outputs are in `G:\My Drive\Projects\Old\refGW\hodge_vertex_conifold\lemmaA\round2\R2_onepoint\referee\`. HVC below means `G:\My Drive\Projects\Old\refGW\hodge_vertex_conifold`. Date: 2026-10-31.

## 0. Verdict

**Sound.** I found no mathematical error. The chain that proves Lemma A for every h is:

- ABDKS (17) gives Thm 1.
- Homogeneity gives Thm 2 (the degree bound).
- Eynard's formula (50) gives Thm 3 (the dictionary P2).
- Thms 2 and 3 together give Cor 4.

Each step is correct and uniform in the genus. Nothing uses induction on the genus or the level, and every hypothesis holds on C_tau.

I checked the two literature inputs against the source, arXiv:2312.16950v3. They are stated and normalised correctly, and their hypotheses hold on C_tau. Their proofs were not re-checked. So the status proved-modulo-literature is justified.

I re-derived every proof step by hand. I also wrote my own implementation of Theorem 1, in a form that works for every rational tau, not only integers. With it I made about 9,500 new comparisons against four kinds of independent data, mostly at framings and (g, n) the author never used:
- Route E's TR-free engine;
- the Mathematica tables in lemmaA\data;
- Route A's Marino-Vafa (Schur-side) values;
- three exact EO recursions.

There were no mismatches.

The issues I found are all minor or cosmetic. None affects a theorem:
- one verification row (V11) reports comparisons that never took place;
- the comparison script counts missing data as passes;
- Cor 6 is misstated at h = 0;
- Cor 4 should exclude the empty multiset R;
- Thm 2 assumes integer tau, which it does not need.

## 1. Literature check (arXiv:2312.16950v3, read on arxiv.org/html)

### 1a. ABDKS Proposition 3.1 and eq. (17)

What the paper states:
- **Setting.** Sigma = P^1; y (or e^y) is an affine coordinate z; x = Sum_i log(z - a_i)/alpha_i + R(z).
- **The points a_i.** They are the LogTR-vital points, with 1/alpha_i the residues of dx there. By Def. 2.4, a point is LogTR-vital (here, for the dual curve) when it is a simple pole of the logarithmic differential and the other differential is regular there.
- **(15).** xhat = Sum_i Sc(alpha_i hb d_y)^{-1} log(z - a_i)/alpha_i + R. ABDKS's calligraphic S is the author's Sc, so this is the author's xhat.
- **(16).** w_i = exp(u_i hb d_y/2) z_i and wbar_i = exp(-u_i hb d_y/2) z_i.
- **(17).**

      ^log omega^(g)_n = [hb^{2g-2+n}] (-1)^n prod_i ( Sum_r (d_i 1/dx_i)^r [u_i^r] e^{-u_i (Sc(u_i hb d_{y_i}) xhat_i - x_i)} ) (-1)^{n-1} Sum_{n-cycles} prod_i sqrt(dw_i dwbar_sigma(i)) / (w_i - wbar_sigma(i)).

- **Rem. 3.2.** For y = z, w = z + u hb/2 and sqrt(dw dwbar) = dz.
- **Rem. 3.3.** log(z - a_i) is to be read as log((z - a_i)/(z - p_0)), with p_0 a pole of dy. Here p_0 = oo, so R = 0.
- **Rem. 3.4.** LogTR = TR when y is affine.
- **Rem. 2.6.** A simple pole of the logarithmic differential at which the other differential also has a pole is not LogTR-vital and does not contribute. For the dual of C_tau this is v = oo.
- **Status of the result.** Prop. 3.1 is stated as a direct corollary of eq. (13), i.e. Theorem 2.8, via the Cauchy determinant (30). Theorem 2.8 is proved in s.5 under standing assumptions: simple, disjoint zeros of dx and dy; dy regular and non-zero at the zeros of dx, and vice versa. There is no condition on the alpha_i.

Check on C_tau:
- dX = (tau - v) dv/(v(1+v)) has one simple zero, v = tau, and dv has no zeros.
- a_1 = 0 with 1/alpha_1 = tau; a_2 = -1 with 1/alpha_2 = -e2; R = 0.
- v = oo is non-vital.

All hypotheses hold. The author's transcription of (15)-(17) and of the remarks is exact.

### 1b. ABDKS s.3.4-3.5, eqs. (40)-(50)

The setting:
- (40): dx = z dz/((1+alpha z)(1+beta z)(1+gamma z)), with omega_{0,1} = -y dx and B the Bergman kernel.
- (41)-(42): x and the parameters p, q, r as quoted by the author.
- (43): xi_k = (-d_x)^k(-1/z - b). The constant b drops out of every d xi_k and out of q_k, because q_k uses xi_{k+1}.

The Eynard formula, s.3.5, (49)-(50):
- **Hypothesis on y.** Any power series in z with y'(0) != 0.
- **The q_k and s_i.** q_k = -res_{z=0} xi_{k+1} y dx, and 1 + Sum_{k>=1} (q_k/q_0) t^k = exp(-Sum_i s_i t^i).
- **(50).** omega^(g)_n = q_0^{2-2g-n} Sum_k Int exp(Sum s_j kappa_j) Lambda(p) Lambda(q) Lambda(r) prod psi_i^{k_i} d xi_{k_i}.

All of this matches the author.

Citation details:
- ABDKS attribute (50) to [Eyn14] = Eynard, CNTP 8 (2014) 541-588. The author's own label [Eyn14a] belongs, in ABDKS's list, to Eynard's ICM 2014 overview. This is a cosmetic mismatch.
- ABDKS attribute the plain triple-Hodge formula (44) to Bouchard-Marino, Zhou and Chen.

## 2. Re-derivations (all correct)

**Thm 1.**
- (i) The coefficient extraction prod mu_i [x^mu]F = iterated Res at v = oo is valid because omega_{g,n} is regular at (oo, ..., oo).
- (ii) Integration by parts: Res x^{-mu} d(eta/dX) = mu Res x^{-mu} eta. It is valid over the field of the later variables. The operators in distinct variables commute with the residues.
- At n = 1 the cycle factor 1/(u hb) has no u^0 term. Its -1/(mu hb) piece sits at hb^{-1} and does not contribute.
- (iii) The operator identity -mu Sc(mu hb d) xhat = log Q_mu follows from sinh(N t/2)/sinh(t/2) = Sum_{m in M(N)} e^{m t}.
- (iv) The signs combine to (-1)^{n-1}.

**Two-point difference form.** It follows from [v^-1] G(v)/((v-a)(v-b)) = G_+[a, b], the divided difference of the polynomial part. I re-derived it independently.

**Thm 2 (homogeneity).** The argument is correct. It also holds without integrality: see issue 5.

**Thm 3, steps (a)-(f).**
- With v = tau + kappa z: alpha = kappa/tau, beta = kappa/e2, gamma = 0, so p = 1/tau, q = -1/e2, r = 1.
- The limit gamma -> 0 is harmless: the third pole of dx goes to z = oo, and r -> -alpha beta stays finite.
- Lambda(1/t) = (-t)^{-g} Lambda^vee(-t), which gives the factor sigma^{-g}.
- q_0 = kappa and q_1 = kappa alpha beta = -kappa, by direct residues.
- The moment lemma gives q_k = (-1)^k kappa, hence s_i = (-1)^{i+1}/i.
- The degree involution then gives K = exp(-Sum kappa_j/j) and the overall factor sigma^{1-g}. The sign bookkeeping, (-1)^{6g-3+n} = (-1)^{n+1}, is right.
- The dilaton-leaf pushforward gives K (as in brief (K)).

**Moment lemma.** Proved; xi_0 = 1/(v - tau) checked by Lagrange inversion.

**Thm 5.**
- The sum of finite residues is a centred divided difference: (2 sinh(b d/2)/b)^{N-1} P(0)/(N-1)!.
- Duality: [G(d_v) P](0) = [P(d_phi) G](0).
- Conjugation by e^phi turns d + 1 into d.
- The centred product identity: prod_{m in M(M)} (s + m a) = M! [phi^M] e^{s phi} S(a phi)^{M+1}, via 1 + w = e^{a phi}.

**Cor 6.** The prefactor collapses: n1! n2! (j-1)!/(sigma j^2 (n1-1)! (n2-1)!) = (j-1)!. The h = 1 case reproduces Route E's sigma Xi_1(j) = (sigma(j-2) - 1)/24.

**Cor 7.**
- The string and dilaton steps are correct; they are brief (SD).
- In the limit sigma -> oo the formula gives B_{2h}/(2h(2h-2)) = chi(M_h).

**F1.** Murnaghan-Nakayama: correct.

**Cor 4.** Correct: Thm 2 and Thm 3 at tau = 1, 2, ..., where every Pi_j > 0. Then sigma^h Xi_h(R) is a polynomial in sigma, so it vanishes identically.

## 3. Independent computations

All comparisons below are modulo p = 2^31 - 1, with a second prime, 2^31 - 19, on subsets. The arithmetic is exact native integers, and engine rationals are reduced mod p. A mismatch mod p is a certain mismatch; an accidental match has probability about 1/p.

My main implementation, `ref_xy_mod.pl`, evaluates Theorem 1 in **formal form**: Q_mu/v^mu = exp(e2 mu log(1+w) + Sum_k (2k-1)! hb^{2k} [phi_k(tau mu) tau^{-2k} w^{2k} - phi_k(e2 mu) e2^{-2k} w^{2k}(1+w)^{-2k}]), with w = 1/v. So it works at every rational tau, and it does not reuse the author's finite-product code. The cycle factors are handled with complete homogeneous polynomials of the nodes, and the result is passed through P2.

| # | Check | Script | Range | Values | Mismatches |
|---|---|---|---|---|---|
| R1 | Thm 1 + P2 vs Route E engine (TR-free) | sweep1.sh, sweep1b.sh | 142 configurations (n = 1..7, g = 1..7, framings listed below) | 3058 (1250 in range, all zero) | 0 |
| R2 | Thm 1 alone (no P2) vs exact EO, tau = 1 (status_trC1) | ref_cmp_tr.pl | g <= 7, n <= 6 (includes n = 6 at g = 2) | 403 | 0 |
| R3 | Thm 1 alone vs exact EO, tau = -1/2 (status_trhalf) | ref_cmp_half.pl | n = 1 up to g = 8; n = 2 up to g = 7 (odd parts) | 316 | 0 |
| R4 | Thm 1 + P2 vs main session's Mathematica EO, first layer, tau = 2, 3 | ref_mma_first.pl | g <= 3, n <= 3 (9 values the author did not list) | 24 | 0 |
| R5 | vs lemmaA\data q-tables, evaluated directly with my own parser (no interpolation) | ref_qtable.pl, ref_cmp_q.sh | g = 2, 3, 4; r = 1, 2; tau in {1/2, 3, -2/7, 5/3, 4} | 358 | 0 |
| R6 | vs Route A's Marino-Vafa Schur-side values (A_mv\a2_g*_r*.out) | ref_cmp_mv.pl | g = 2, 3; r = 1, 2; tau in {2/7, 3/2, 5, -7/3, 7/11} | 57 | 0 |
| R7 | Thm 1 at n = 1 (formal) vs (6.1), i.e. Theorem 5 | ref_61.pl | h <= 12, j <= 3h, 13 framings, one repeated with the second prime | 3028 | 0 |
| R8 | (6.1) vs engine, including the author's new genus-7 table | ref_61.pl | h <= 7, j = 0..3h, 14 runs | 1274 | 0 |
| R9 | gap formula vs (6.1) at j = 0 | ref_61.pl | h = 2..12, 14 runs | 154 | 0 |
| R10 | genus 0 vs (J-1)...(J-r+3) | ref_genus0.pl | r = 3..6, 4 framings | 840 | 0 |
| R11 | degree bound of the formula beyond the engine's range | ref_highgenus.out | n = 1, g = 12 and n = 2, g = 9 (each at tau in {1/2, 3, -5/3}); n = 3, g = 7; n = 4, g = 5 | 344 range zeros; 30 sampled values outside the range, all non-zero | -- |

Framings used in R1:
- n = 1, g = 2..6: tau in {1/2, -2/7, 3/2, 3, 5, 7, -5/3, 1/9, 11/4}.
- n = 1, g = 7: tau in {1/2, -2/7, 3/2, 5, 6, 7, 8, -5/3}. This closes the author's unfinished tau = 5..8 runs.
- n = 2, g = 1..6: tau in {1/2, -2/7, 3/2, 3, 4, 5, -5/3}.
- n = 3, g = 1..6: tau in {1/2, 3, -2/7, 5/3, 4}.
- n = 4, g = 2..5: tau in {2, 1/2, 3, -2/7}.
- n = 5, g = 4, 5: tau in {2, 1/2, 3}.
- n = 6 and n = 7 at g = 5.

Notes on these runs:
- **Pi_j = 0.** At negative tau some Pi_j vanish (for example Pi_7(-2/7)). In every such case I checked that [x^mu]F = 0, as P2 requires, and found no inconsistency.
- **My own script bug.** The first genus-0 run failed everywhere. The cause was my sigma^{-1} handling at g = 0 (the modular power routine does not take negative exponents). After the fix all 840 values match. The fix does not affect R1-R9, which have g >= 1.

**Negative controls** (`ref_xy_mod_variant.pl`):
- Replacing Sc by S in the dressing makes 14 of 24 values fail.
- Reversing the sign of the cycle shift at n = 3 makes 8 of 10 fail. (At n = 2 the product is invariant under this flip, so it is not a test there.)
- Both wrong variants still produce **all** the range zeros.

So the range zeros, and the author's homogeneity test V20, are structural consequences of the shape of (17). They do not validate Theorem 1. The validation comes from the out-of-range agreement: about 1800 engine values in R1, plus R2-R6.

**Re-runs of the author's scripts:**
- `xy_npoint.pl 1 7 1 21 engine_h7_r1.txt`: 21 real OK lines (3.5 min).
- The same command without the table prints "(no data)" on every line and still ends "total 21, mismatches 0" (see issue 2).
- `moments_check.pl`: Res xi_k dv = 1 for k = 0..6.
- The h7 table copies are identical.

## 4. Issues

1. **(minor) V11 never compared anything.**
   - The report lists V11 as 21 engine comparisons at (g, n) = (2, 5), tau = 1, 0 mismatches. There is no `E_formulas\engine_h2_r5.txt`, and every line of `runs_np_g2_n5_t1.out` reads "(no data)".
   - So the engine total is 1051, not 1072. The claim "n <= 5 vs engine" is unsupported by the author's runs. The only n = 5 data are 26 values against the exact EO at tau = 1 (V13).
   - Fix: correct s.10 and the totals. R1 now supplies engine comparisons at n = 5 (g = 4, 5) and at n = 6, 7 (g = 5).
2. **(minor) The comparison scripts count missing data as passes.**
   - `xy_npoint.pl` (and the family of scripts built the same way) counts "(no data)" lines in the total but not as mismatches. A run with no table therefore ends "total N, mismatches 0", so the summary lines in `runs_all_checks.txt` do not show that a comparison took place.
   - Fix: print the compared and no-data counts separately, and fail when no table is found.
   - My re-run at genus 7 and the R1 runs confirm the genus-7 claims anyway.
3. **(cosmetic) Cor 6 is wrong at h = 0.**
   - The generating-function form sums over h >= 0. Its lambda^0 coefficient is (j-1)!/(j+1)! = 1/(j(j+1)), which is not polynomial in j and is undefined at j = 0.
   - By the brief's definition Xi_0(j) = 0: every term vanishes by degree, since 2m > m - 2.
   - Fix: sum over h >= 1, or declare the h = 0 term a convention.
4. **(cosmetic) Cor 4 should exclude R = empty.**
   - Cor 4 and the summary say "every multiset R of positive integers with D(R) <= 2h-3". This includes R = empty, but Xi_h(empty) = Int K Omega != 0 (that is Cor 7).
   - The proof needs n = r >= 1, so it covers only non-empty R.
   - Fix: add "non-empty".
5. **(cosmetic, simplification) Thm 2 does not need integer tau.**
   - In the formal form, Q_mu(v; c)/v^mu = (1+cw)^{e2 mu} exp(...), every power of w carries one factor of c or hb. So R_mu(c, hb) is a homogeneous polynomial of degree |mu| for every tau != 0, -1.
   - Hence Lemma D holds at every such tau directly, and the polynomial-in-tau step of Cor 4 can be dropped. R1 and R3 confirm the formal form at seven non-integer framings.
6. **(cosmetic) V20 is not evidence for Thm 1.** See the negative controls in section 3; it should not be counted among the checks that validate the formula.
7. **(cosmetic) Citation label.** [Eyn14a] versus ABDKS's [Eyn14] (section 1b).

## 5. Claim by claim

| Claim | Assessment |
|---|---|
| Thm 1 | Confirmed: literature statement checked; derivation re-done; independent formal implementation matches the engine, three exact EO codes, the q-tables and the MV values, including non-integer tau and n up to 7. |
| Thm 2 | Confirmed (Lemma D for every tau != 0, -1, not only integers). |
| Thm 3 (P2 = Lemma B) | Confirmed (re-derived; every engine comparison goes through P2). |
| Moment lemma | Confirmed. |
| Cor 4 (Lemma A for all h) | Confirmed for R non-empty. No induction is used, so the brief's (OL) holds outright. |
| Thm 5 | Confirmed (elementary proof re-done; R7). |
| Cor 6 ((6.1) for all h) | Confirmed for h >= 1 (issue 3). |
| Cor 7 (gap formula) | Confirmed (R9; the engine at j = 0 in R8). |
| Two-point difference form | Confirmed. |
| F1 | Confirmed. |
| Computer verification | Confirmed, with the corrections in issues 1 and 2. |

**Uniformity in the genus.** There are no base cases and no induction. The hypotheses are sigma != 0, 2g - 2 + n > 0 and n >= 1; integer tau >= 1 is used only for the finite-product form, and Pi_j != 0 holds for tau >= 1.

## 6. Dependence on TR and the rules

- **What depends on TR.** Thms 1 and 3, and therefore Cor 4, Cor 6, Cor 7 and the two-point form, rest on TR results for the Hodge curve C_tau: the log-TR x-y swap, and Eynard's / BKMP formula for C^3. This is TR applied to the Hodge CohFT. Brief s.0 allows it, and the user's relayed message confirms it is an allowed input.
- **What does not.** Nothing in the R2 work applies TR to the resolved conifold x C^2 theory. Input (H) is untouched. Thm 5, the moment lemma and F1 are TR-free.
- **Cross-checks.** My comparisons against the engine, the q-tables and the MV values are TR-free. R2-R4 are TR cross-checks, which brief s.0 also allows.
- **Downstream.** The step "Lemma A => Corollary R (all-genus conifold regularity)" belongs to the main session and is not part of R2. The author flags it correctly.

## 7. Obstructions and simpler paths

- **The obstructions the author lists are real.** They are:
  - the dependence on literature;
  - the use of TR, which is allowed;
  - no all-order control on the Marino-Vafa side;
  - multi-leaf closed forms and the Hurwitz-limit extra zero, not obtained.
  The tau = 5..8 genus-7 gap is now closed by R1.
- **What Lemma A actually needs.** Only (17) (giving Thms 1 and 2) and (50) (giving Thm 3). Thm 5 and Cors 6-7 are by-products.
- **One simplification.** Drop integrality (issue 5).
- **A possible TR-free route.** ABDKS prove (17) through a free-fermion Cauchy determinant. The Bychkov-Dunin-Barkowski-Kazarian-Shadrin closed formulas for Orlov-Scherbin n-point functions are TR-free and use inputs brief s.0 allows. They could give a TR-free proof of Thm 1, but only if the K-twisted series Sum [x^mu]F_{g,n} at v = oo is first identified with such a tau function. Not done here.
- **No simpler argument for the degree bound.** Rescaling (v, c) -> (lambda v, lambda c) shows [x^mu]F_{g,n}(c) = const * c^{|mu| - (2g-2+n)}. Lemma D is therefore equivalent to regularity at c = 0, where the branch point collides with the log points. The x-y swap is what makes that regularity visible.

## 8. Files (referee folder)

| File | Contents | Usage |
|---|---|---|
| ref_xy_mod.pl | formal Theorem 1 + P2, mod p | perl ref_xy_mod.pl TAU G N JMAX [TABLE or -] [PIDX] |
| ref_xy_mod_variant.pl | negative controls | VARIANT=1 or VARIANT=2, same arguments |
| sweep1.sh, sweep1b.sh, sweep1.out, sweep1b.out | R1 | |
| ref_cmp_tr.pl, ref_cmp_one.out | R2 | |
| ref_cmp_half.pl, ref_cmp_half.out | R3 | |
| ref_mma_first.pl, ref_mma_first.out | R4 | |
| ref_qtable.pl, ref_cmp_q.sh, ref_cmp_q.out | R5 | |
| ref_cmp_mv.pl, ref_cmp_mv.out | R6 | |
| ref_61.pl, ref61_runs.out | R7-R9 | |
| ref_genus0.pl, ref_genus0.out | R10 | |
| ref_highgenus.out | R11 | |
| ref_cmp_one.pl | placeholder only (an aborted edit) | |
| tmp*.txt | scratch | |

