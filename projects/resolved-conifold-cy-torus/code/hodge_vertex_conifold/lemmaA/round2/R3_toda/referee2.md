# Lemma A round 2: referee2 R3 (saved by the main session from the workflow journal, 2026-10-31)

Subagents cannot write report files, so this file holds the agent's structured return value: the fields first, then the full report text (report_text) verbatim.

## confirmed

- The first mKP equation P_2 holds for Z_MV with tau_1 = Z(e^{tau x} o p). Verified in bilinear (Hirota) form on Z itself (no logarithm, no cumulants) by an independent Schur-side implementation: through x^120 for all |lambda| <= 8 and through x^80 for |lambda| <= 6. Both controls (partner tau_{-1}; reversed order) fail in every case.
- P_k = D_{t_k} + p_k(-D~) is the u_k-coefficient of the Ueno-Takasaki 2D Toda bilinear identity at equal second times with n' = n+1 (re-derived). P_3, P_4, P_5, P_6 verified in bilinear form through x^60 for |lambda| <= 5. The Jimbo-Miwa family (n' = n-1) gives the same P_2 and P_3.
- Lemmas 2.1-2.3 re-derived. Every term at (x^n, q_R) lies in the single exponential sector e^{u(J+2)}, because log Z_MV contains all unstable genus-0 terms. The empty-profile table (L0) was recomputed exactly, including that no term lies below the nominal order.
- Route A Prop. 2.1 / Cor. 2.2 degree bound deg_s Q_g <= N confirmed from MV data for genus 1-8 (56 cases, B = (1), (2), (3), (5), (2,3), (4,4), (2,2,3)). Q_g(-1) agrees with exact Hodge values for genus <= 7.
- Theorem M re-derived and re-implemented from scratch as a general Hirota leading-order engine. It equals the author's formula (times 2) in 592 exact evaluations at G <= 6. It holds exactly on Route E tables for P_2 (G <= 6, r <= 2), P_3 (G <= 3), P_4 (G <= 2) and P_5 (G <= 2), and mod p at G = 3, 4 with three- and four-part R (16 relations).
- Theorem M at G = 7 from Hodge data alone (exact genus-7 engine points plus tables): R = (13) and (14), both parities, are exact identities in tau (22 rational points, above the degree bound).
- Theorem M (P_2) and Theorems M_3, M_4 hold at G = 7, 8, 9, 10 outside the range (excess -1 to -3), mod p, with all genus >= 1 values from the referee's own MV engine (validated against more than 2,200 exact Hodge-table values and the 4 exact genus-7 points). 44 relations (33 non-trivial), including R = (6,8), (5,9), (2,2,11), (7,9), (6,10), (8,10), (9,11). All 22 perturbation controls give nonzero residuals.
- Lemma A holds directly at genus 7-10 (MV, mod p) for 16 in-range profiles at the top level D = 2G-3: Xi_7 at (12), (2,11), (6,7), (2,2,10), (3,10); Xi_8 at (14), (2,13), (7,8), (3,12); Xi_9 at (16), (2,15), (8,9); Xi_10 at (18), (2,17), (9,10); plus Xi_10(17).
- Route A Theorem B (used in proof (i)) holds at g = 7 (R = (12), (13), (14), (15)), g = 8 ((14), (15), (16)), g = 9 ((17), (18)) and g = 10 ((19)), mod p with own values.
- Corollary I_2 counting checked mechanically: every term of every relation enumerated for G = 2..12 and all R (6259 relations), and for G = 13..21 with r <= 3. No surviving non-target term, no genus-G term above level D, exactly one target per relation. Corollary I_P checked for P_3 (G <= 12) and P_4 (G <= 6).
- Theorem-B induction step of proof (i) checked mechanically for G = 2..12 (1957 steps): every profile at every level D <= 2G-3 is reached, and every cut and join term is killed by (H1) or (H2) by excess counting with string reduction.
- Base case h = 2 needs no hypothesis. The hypotheses at (G, D) are exactly (H1) and (H2), and there is no circularity: Theorems M, M_P and B are hypothesis-free identities. The final statement implies (OL) for every G >= 3 and every odd and even D <= 2G-3, hence Lemma A for all h (directly, or with brief s.3).
- The lambda_g limit of Theorem M holds for every G and R (algebra redone).
- Literature: the content-product weight r(c) = e^{tau x c} satisfies the hypothesis of the 2D Toda theorem. Specialising the second times to q^rho is legitimate coefficientwise. tau_n(p, q^rho) = Z(e^{n tau x} o p).

## issues

### [1]

#### description

The w^{n+r} leading coefficient is well defined only because every term at (x^n, q_R) lies in the single exponential sector e^{u(J+2)}. That in turn holds only because F = log Z_MV includes all unstable genus-0 terms (l = 1, 2), giving for example d_y^2 G^(0)_0 = e^{2u}/2. With stable-only genus-0 series one gets two sectors, (e^{2u}-1)/2, as Route A's referee found. The report says the unstable terms are included but never draws this consequence.

#### location

Report s.2 Lemma 2.1 and s.4 proof, step 3 (passage to the dilaton point)

#### severity

minor

#### suggested_fix

Add: each factor d_y^v (J_U + D)^k G^(g)_U lies in e^{u(J_U+v)} C[w, 1/w], so every term at (x^n, q_R) lies in e^{u(J+2)} C[w, 1/w] and no sector projection is needed.

### [2]

#### description

Proof (ii) needs P_k for every k up to 2G-2, yet the form of P_k is said to be 'fixed by direct checks for k <= 4'.

#### location

Report s.3 'Members used' and obstructions (proof (ii))

#### severity

minor

#### suggested_fix

Derive P_k from the Ueno-Takasaki bilinear identity with n' = n+1 at equal second times: the right side is tau_{n+1}(t) tau_n(t') and the u^0 term fixes the orientation; the u_k-coefficient is (D_{t_k} + p_k(-D~)) tau_0.tau_1 = 0. The referee re-derived this and checked P_2..P_6 in bilinear form.

### [3]

#### description

Proof (i) via Cor. E rests on Theorem D, i.e. on KP for the MV series and Takasaki-Takebe's dispersionless Hirota equation. Proof (i) via Theorem B needs only the MV formula, its cut-and-join equation, the single 2D Toda equation P_2 and Route A's dilaton calculus. The summary does not separate the two.

#### location

Summary / s.6 (dependencies of proof (i))

#### severity

minor

#### suggested_fix

State the two dependency sets explicitly.

### [4]

#### description

log Z_MV has no constant term, so the G^(g)_0 obtained from F lacks <Omega>_g. It lies in w^(2g-2) Q[w] only up to an additive constant; harmless because only derivatives occur.

#### location

Report s.2 Lemma 2.1, empty profile, g >= 2

#### severity

cosmetic

#### suggested_fix

Add 'up to an additive constant (only derivatives occur)'.

### [5]

#### description

'Every multiset R of integers >= 2 with D(R) <= 2h-3' includes R = empty (D = 0). The induction (D = 1, ..., 2h-3) does not cover it, and it is not part of Lemma A (round-1 brief: J >= 1).

#### location

Theorem L statement (claim 1, report s.0)

#### severity

cosmetic

#### suggested_fix

Say 'every non-empty multiset'.

### [6]

#### description

The report says the exact genus-7 Hodge-engine evaluation (g7_points.pl) did not finish within 2 h. R3_toda\out\g7_points_A.out contains all four points, Xi_7 at (13), (14), (2,13), (2,14); the run finished at 19:45, after the report was written.

#### location

Report s.7 last paragraph; obstructions

#### severity

cosmetic

#### suggested_fix

Update the text. Add the referee's exact genus-7 check: with these values Theorem M at G = 7, R = (13) and (14), both parities, holds as exact identities in tau (referee2\out\g7_exact.out).

### [7]

#### description

The referee could not confirm the section reference 'Okounkov math/0004128 s.2.7'. The content-product 2D Toda statement should be cited precisely.

#### location

Report s.0 (I2) citations

#### severity

cosmetic

#### suggested_fix

Cite Orlov-Scherbin (math-ph/0003011, Physica D 152-153 (2001); Theor. Math. Phys. 128 (2001)) and Harnad-Orlov, Comm. Math. Phys. 338 (2015), for the general content-product statement, and Okounkov for double Hurwitz numbers.

### [8]

#### description

The claim is correct but the reason is not given. The relation has degree at most 2G + eps + 2 in tau and every term contains a factor tau, so 2G+3 non-zero points suffice for both parities.

#### location

Report s.7 table, 'polynomial identity in tau (2G+3 points)'

#### severity

cosmetic

#### suggested_fix

State the degree bound and the factor tau.


## notes

Verdict: sound. I found no mathematical gap. The proof of Lemma A for all h is correct modulo the MV formula and the 2D Toda theorem for content-product tau functions.

All referee scripts and outputs are in G:\My Drive\Projects\Old\refGW\hodge_vertex_conifold\lemmaA\round2\R3_toda\referee2 (outputs in out\). The core libraries are R2mv.pm (independent MV Schur engine), R2rel.pm (independent Hirota leading-order engine) and R2tab.pm.

Two long logic runs were still going at submission. logic_P2_G13to30_r3.out had reached G = 21 with 0 survivors. logic_P2_G40_r2.out had not finished. Both are wrapped in timeout 9000 and stop on their own.

No wolframscript was run and no topological recursion was used.

## refuted

- Statement in report s.7 and the obstructions that 'the exact Hodge-engine evaluation of genus-7 points (g7_points.pl, HODGE_NCH = 4) did not finish within 2 h locally': R3_toda\out\g7_points_A.out contains all four genus-7 points (finished at 19:45). This is a factual statement about the runs, not a mathematical one; with these points Theorem M at genus 7 holds exactly.

## route

R3 (Toda/KP), round 2: Lemma A for all genera from the first-time-set (mKP) equations of the 2D Toda hierarchy of the Marino-Vafa series; second, independent referee

## verdict

sound

## report_text

# Second referee report on R3 (Toda/KP), round 2: Lemma A for all genera

Referee 2 (adversarial, independent). My scripts are in HVC\lemmaA\round2\R3_toda\referee2 and my outputs in referee2\out (HVC = G:\My Drive\Projects\Old\refGW\hodge_vertex_conifold).

- **Independence.** I wrote all of my code from scratch. It shares nothing with the author's mkp_*.pl, mkp_eqs.pm or XiData.pm, nor with mvlib.pm or rmv.pm.
- **Arithmetic.** Low-genus checks are exact (Perl Math::BigRat). High-genus checks are mod p, with p = 1000000007 or 998244353 and tau chosen at random mod p.
- **Rules.** No wolframscript was run. No topological recursion of any kind was used: no spectral curve, loop equation, EO or BKMP, dictionary P2, mod-P2 result or TR script. The only Hodge data used are Route E's exact tables and the exact genus-7 engine points.

## Verdict: SOUND

I assumed the argument was wrong. I then tried to break every step that is asserted for all G, numerically at G = 7, 8, 9 and 10 and logically. I found no failure.

Theorem L (Lemma A = M_irr(h) for every h) is correctly proved modulo the cited literature: the Marino-Vafa formula, and the 2D Toda hierarchy for content-product tau functions. This is exactly the status the author claims.

The new mechanism is genuine:
- The Toda shift tau_1 = Z(e^{tau x} o p) enters the first-time-set equations only through A = (e^{beta E} - 1)F.
- At the dilaton point, E^k produces the falling factorial (delta)_k.
- The first mKP equation therefore gives a relation that is inhomogeneous in the excess. Its diagonal coefficient -2 Pi_2 sigma^G tau delta_G(R+2) is not zero.

That is precisely what the homogeneous KP relations of round 1 lacked.

All issues I found are minor or cosmetic (Section 5).

## 1. Independent tools

| file | what it does |
|---|---|
| R2mv.pm | Marino-Vafa Schur side mod p, from scratch: Murnaghan-Nakayama on beta-sets; w_nu = (1/H) exp((tau+1/2) cs x + sum c_k P_2k x^2k) with conjugate pairing; unit legs; labelled set-partition cumulants; Laguerre transform (Route A Prop. 2.1). One run gives Xi_g(B) for every g <= gmax. |
| R2tab.pm | Reader for Route E's exact tables. |
| R2rel.pm | Dilaton-point leading coefficient at (x^n, q_R) of ANY Hirota member P(D) tau_0.tau_1 = 0, derived from first principles (details below). Exact or mod p. Includes the generator of the members P_k = D_{t_k} + p_k(-D~). |
| r2_mkp_bilinear.pl | Hirota equations checked in bilinear form on Z itself, with no logarithm and no cumulants. |
| r2_logic.pl | Mechanical check of Corollaries I_2 / I_P and of the Theorem-B induction step. |
| r2_highg.pl, r2_g7_exact.pl, r2_thmB_highg.pl, r2_degree.pl, r2_validate.pl, r2_rel_tables.pl | Tests below. |

R2rel.pm works as follows:
- It expands Faa di Bruno over set partitions.
- It uses d_{t_i} = i x sigma Pi_i d_{q_i}, block weight 2g-2+|B|+k_B, and factor (J_U)^(v) (delta_g(U))_k X_g(U).
- It computes the empty-profile coefficients L0(g,k,v) by exact Laurent-polynomial arithmetic in w, starting from G^(0)_0 = u - 3u^2/4 + u^3/6 and y d_y G^(1)_0 = (w-1)(w-1-1/sigma)/24.
- It asserts in every case that no term lies below the nominal w-order.

**Validation of R2mv.** It reproduces:
- Route E's exact Hodge tables at genus 1-6, more than 2,200 value comparisons over all runs, 0 mismatches;
- the four exact genus-7 Hodge-engine points of R3_toda\out\g7_points_A.out (4/4);
- the author's mod-p genus-7 values at the same tau, for example Xi_7(2,13) = 1413483.

## 2. Every step asserted for all G, tested at G = 7, 8, 9, 10

| # | step | why it holds for all G | referee's test | result |
|---|---|---|---|---|
| S1 | First mKP equation P_2 for Z_MV (s.3) | 2D Toda for content-product tau functions | Bilinear form on Z, all p_lambda with \|lambda\| <= 8 through x^120 (67 monomials), and \|lambda\| <= 6 through x^80 (30). Connected contributions up to genus about 55. | 0 failures; both controls (partner tau_{-1} = Z(e^{-tau x} o p), reversed order) nonzero in 97/97 |
| S2 | P_k for k >= 3 (proof (ii)) | u_k-coefficient of the 2D Toda bilinear identity (s.4) | Bilinear form, P_3, P_4, P_5, P_6 through x^60, \|lambda\| <= 5 (19 monomials each) | 0 failures; controls nonzero in all |
| S3 | Lemma 2.1 = Route A Prop. 2.1 / Cor. 2.2 (deg_s Q_g <= N and Q_g(-1) = Xi_g) | string and dilaton equations (proved) | [s^(N+1)]Q = [s^(N+2)]Q = 0 from MV data, genus 1..8, B = (1), (2), (3), (5), (2,3), (4,4), (2,2,3): 56 cases. Q_g(-1) agrees with the exact Hodge values (genus <= 7). | all pass |
| S4 | Lemmas 2.2 and 2.3, including the empty-profile table | elementary | Recomputed exactly for v = 1..4 (table below). No term below the nominal order. | equal to the report's table |
| S5 | Theorem M (s.4) and Theorem M_P (s.5) | formal leading order of an exact identity | Own engine with own MV values, outside the range (excess -1 to -3), both parities (cases below) | 44 relations, 33 of them non-trivial, all exactly 0 mod p; each of the 22 perturbation controls (genus-G target + 1) gives a nonzero residual |
| S5' | Theorem M at G = 7, exactly | | Hodge data only (exact genus-7 engine points and tables): R = (13) and (14), both parities, at 22 rational tau, more than the degree bound | exact identities in tau, so these cases are proved |
| S6 | Corollaries I_2 / I_P (counting) | combinatorics | Enumeration of every term. P_2: all R at all levels for G = 2..12 (6259 relations); G = 13..21 with r <= 3 (one run, still going). P_3: G = 2..7 for all R, G = 8..12 with r <= 2 (540). P_4: G = 2..6 (47). | no surviving non-target term; no genus-G term above level D; exactly one target per relation |
| S7 | Route A Theorem B (proof (i)) | MV cut-and-join (LLZ) | Own values: g = 7, R = (12), (13), (14), (15); g = 8, (14), (15), (16); g = 9, (17), (18); g = 10, (19) | all 0 (10 relations) |
| S7' | Theorem-B induction step of proof (i) | combinatorics | Mechanical, G = 2..12, every T at every level: 1957 steps | every T reached; every cut and join term killed by (H1) or (H2) |
| S8 | lambda_g limit (s.4) | algebra | Redone by hand | holds for all G |
| S9 | The conclusion itself | | Direct MV values in the range, at the top odd level D = 2G-3 (Xi_10(17) is at D = 16). G = 7: (12), (2,11), (6,7), (2,2,10), (3,10). G = 8: (14), (2,13), (7,8), (3,12). G = 9: (16), (2,15), (8,9). G = 10: (18), (2,17), (9,10), (17). | all 16 vanish |

The empty-profile values recomputed in S4 are:

| (g, k) | L0(g, k, v) |
|---|---|
| (0, 0) | (v-1)!/2 |
| (0, 1), (0, 2) | (v-1)! |
| (0, k >= 3) | 0 |
| (1, 0) | (v-1)! X_1(1), with X_1(1) = -(1 + 1/sigma)/24 |
| (1, k >= 1) | 0 |
| (g >= 2, any k) | 0 |

The relations tested in S5 are:

| G | prime | P_2 at R = | P_3 at R = | P_4 at R = |
|---|---|---|---|---|
| 7 | 998244353 | (13), (14), (15), (6,8), (5,9), (2,2,11) | (13), (6,8) | (13) |
| 8 | 1000000007 | (15), (16), (7,9), (6,10) | (15) | |
| 9 | 998244353 | (17), (18), (8,10) | (17) | |
| 10 | 1000000007 | (19), (20), (9,11) | (19) | |

These relations are genuinely non-trivial. At G = 10 and R = (19), for example:
- the (a) terms carry X_g(2,19) for every g = 0..10, with the nonzero falling factorials (21-2g)_(21-2g);
- the (c) terms carry X_g(19) for every g < 10;
- (b) carries X_10(19), which has excess 0 and is nonzero.

The two- and three-part R exercise splittings into two non-empty parts in the product terms.

Further exact checks at low genus, on the tables:
- P_2 is zero for G = 0..3, r <= 2, J <= 7 at two values of tau (104 relations); G = 1..3, J <= 9 (126); and G = 4..6, r <= 2, J <= 13 (258).
- P_3 is zero for G <= 3 (104), P_4 for G <= 2 (44), P_5 for G <= 2 (36).
- The author's printed Theorem M formula equals twice my engine's value for P_2 = (D_1^2 + D_2)/2 in all 592 evaluations.
- Mod p at G = 3 and 4 with three- and four-part R: 16 relations, all 0.

## 3. Logic

### 3.1 Theorem M
I re-derived it. These are the points a referee must check, and all of them hold.

**(a) Single exponential sector.** At fixed (x^n, q_R) every term is e^{u(J+2)} times a Laurent polynomial in w. So the identity between power series in y continues to an identity between rational functions, and the w^{n+r} coefficient is well defined.

This holds because F = log Z_MV contains every unstable genus-0 term, M_0(mu) = |mu|^{l-3} for l = 1, 2. The following are each single-sector (the last two re-derived from the Lambert recursion (D + x) A_{-2} = A_{-1} = e^{xu}/x):
- G^(0)_0 = u - 3u^2/4 + u^3/6, d_y G^(0)_0 = e^u (1 - u/2), d_y^2 G^(0)_0 = e^{2u}/2;
- G^(0)_(j) = e^{ju} (1/j^2 - u/(j(j+1)));
- G^(0)_(j1,j2) = e^{Ju}/J.

With stable-only genus-0 series one gets two sectors instead, for example (e^{2u} - 1)/2 (Route A's referee raised this). That does not occur here.

The genus-1 empty function contains a logarithm, and the constant <Omega>_g is absent from log Z. Neither matters: every Hirota block has |B| >= 1, so only derivatives occur.

**(b) No term below its nominal order.** A block has nominal order 2g-2+r(T_B)+|B|+k_B, and these orders add up to n + r. This follows from Lemma 2.1 (a polynomial in w, by Prop. 2.1) and Lemma 2.2. For empty U it follows from the L0 table.

**(c) Leading coefficients.** (J + D) raises the order by one with coefficient J - m = delta, so E^k gives (delta)_k, and d_y gives (J)^(v). Leading coefficients multiply in products, and q_R is extracted over labelled splittings. My independent engine reproduces the author's formula.

So Theorem M follows from (I1)-(I4) uniformly in G, with no inductive hypothesis.

### 3.2 Corollary I_2
I re-derived it by hand:
- **Falling-factorial dichotomy.** For 0 <= delta < k the factor (delta)_k vanishes. delta <= -1 forces g >= 2, which is (H1).
- **Genus-G product factor.** It has k = 1, its partner is (0,1), and T_i is contained in R, so (H2) applies.
- **Count.** Sum (delta_i - k_i) = delta_G(R+2) - 1 <= -2, while empty factors have delta - k in {0, 1}.

The target coefficient -2 Pi_2 sigma^G tau delta, with delta <= -1, is not identically zero. Check S6 confirms all of this mechanically.

### 3.3 From profiles with a part 2 to all profiles
**Route via Theorem B at R (level D-1).**
- R+2, R_i^+ and R_ik all lie at level D.
- The cut term Xi_{h-1}(al, be, R - j_i) has level D-2 and excess 2h-2-D >= 1 at genus h-1 >= 2. For h = 2 we have D = 1 and R empty, so there are no cut terms.
- The join terms have total excess 2h-2-D >= 1. Genus-0 factors have excess <= -2 and genus-1 factors <= 0, so the factor in the range has genus in [2, h-1] (then (H1) applies) or genus h with level <= D-2 (then (H2) applies).
- Unit parts are removed by the string equation, and Xi_g(1) = 0 for g >= 2.
- The induction on (r, min T) closes (S7').

**(OL) itself needs only one step.** Corollary I_2 at R = (D) gives Xi_G(2,D) = 0. Theorem B at R = (D) then gives gamma_D Xi_G(D+1) = 0.

**Alternative route.** Cor. E, or the KP merge relation Pi_2 Pi_D Xi(2,D) = (D+1) Pi_{D+1} Xi(D+1). The latter is Theorem D at R = empty and (al, be) = (2, D), with excess 2G-2-D >= 1; its a = 1 terms vanish by (H2).

### 3.4 Base case, hypotheses, circularity
- **Base case** h = 2, D = 1, R = empty: Theorem M gives -2 Pi_2 sigma^2 tau (-1) X_2(2) = 0, with no hypothesis needed. This agrees with Theorem A.
- **Hypotheses at (G, D).** (H1): Lemma A in genera 2..G-1, supplied by the outer induction. (H2): all levels below D at genus G, supplied by the inner induction, which settles every T at each level.
- **No circularity.** Theorems M, M_P and B are hypothesis-free identities. Cor. E uses exactly (H1) and (H2).

### 3.5 Consequences
(OL) holds for every G >= 3 and every D <= 2G-3, odd and even. Lemma A therefore holds for all h, either through brief s.3 or directly, since every level is reached. Corollary R then follows, given Lemma C and Theorem R; I did not re-examine those.

### 3.6 Proof (ii) is also correct
- **Budget.** Every block has weight >= 0. A genus-G block uses the whole budget 2G, so (|B|, k_B) is (1,1) or (2,0).
- **Level of the genus-G block.** It is at most D_R + k - |alpha| <= D. Equality holds only for alpha = {k} and T_B = R, which is the target, with coefficient -(k-1) sigma^G tau Pi_k delta. Every other genus-G block falls under (H2).
- **Terms without a genus-G block.** Sum (delta_B - k_B) = delta_G(R+k) - 1 <= -2, and empty blocks have delta - k >= 0. So some non-empty block has delta < k, and it vanishes by the falling factorial or by (H1).

Check S6 confirms this mechanically for P_3 and P_4.

## 4. Literature: exact statements and hypotheses

**(I1) Marino-Vafa formula** (Liu-Liu-Zhou, J. Diff. Geom. 65 (2003); Okounkov-Pandharipande, Geom. Topol. 8 (2004)), in Route A's real form (1.1):

[x^(2g-2+l)][p_mu] log Z = sigma^(g+l-1) Pi_mu M_g(mu)/|Aut mu|, with the LLZ unstable conventions.

The normalisation is independently confirmed: my engine, built on (1.1), reproduces Route E's Hodge-integral tables and the exact genus-7 points.

**(I2) Content-product 2D Toda tau functions.** For any content-product weight, tau_n(t, tbar) = sum_lambda r_lambda(n) s_lambda(t) s_lambda(tbar), with r_lambda(n) = prod over boxes of r(n + c), is a tau function of the Ueno-Takasaki 2D Toda hierarchy. References:
- Takasaki 1984;
- Orlov-Scherbin, math-ph/0003011 (Physica D 152-153 (2001)) and Theor. Math. Phys. 128 (2001);
- the general form also appears in Harnad-Orlov, Comm. Math. Phys. 338 (2015);
- Okounkov, math/0004128, proves the Toda property for double Hurwitz numbers. I confirmed this from the arXiv abstract.

Here r(c) = e^{tau x c}, and s_nu(q^rho) = q^{kappa/4}/prod[h], so tau_n(p, q^rho) = Z(e^{n tau x} o p). Specialising tbar = q^rho is legitimate coefficientwise: each p-monomial coefficient of any Hirota equation is a finite polynomial in tbar and in e^{tau x}.

**(I3) First-time-set equations.** Start from the 2D Toda bilinear identity:

oint_inf z^(n-n') e^{xi(t-t',z)} tau_n(t - [1/z], tbar) tau_{n'}(t' + [1/z], tbar') dz = oint_0 z^(n-n') e^{xi(tbar-tbar',1/z)} tau_{n+1}(t, tbar - [z]) tau_{n'-1}(t', tbar' + [z]) dz.

- **Family n' = n+1.** At tbar = tbar' the right side is tau_{n+1}(t) tau_n(t'). Putting t -> t+u and t' -> t-u gives the author's identity. Its u_k-coefficient is (D_{t_k} + p_k(-D~)) tau_0.tau_1 = 0, and the u^0 term fixes the orientation.
- **Family n' = n-1** (Jimbo-Miwa). Its lowest members, (D_1^2 - D_2) tau_1.tau_0 = 0 and (D_1^3 + 3 D_1 D_2 - 4 D_3) tau_1.tau_0 = 0, coincide with P_2 and P_3.

I re-derived both families and checked P_2 to P_6 numerically.

**(I4) Route A Prop. 2.1, Cor. 2.2 and Lemma 3.1.** Proved by the string and dilaton equations and the Lambert recursion. They hold for h >= 2, or h = 1 with r >= 1, and also for B = (1). Confirmed numerically in S3.

**(I5) Route A Theorem B.** Proved modulo the LLZ cut-and-join equation, re-derived by Route A's referee, and confirmed by me at g = 7, 8, 9, 10 (S7).

**Cor. E.** It needs Theorem D, i.e. KP for the MV series plus Takasaki-Takebe's dispersionless Hirota equation. The user ruling allows this input, but it is not needed when Theorem B is used.

## 5. Issues (none affects the conclusion)

1. **(minor) The single-sector point should be stated.** The report says "All unstable terms are included" but does not draw the consequence. Add one sentence: every factor d_y^v (J_U + D)^k G^(g)_U lies in e^{u(J_U+v)} C[w, 1/w], so every term at (x^n, q_R) lies in e^{u(J+2)} C[w, 1/w]. No sector projection is then needed, unlike the stable-only setting, where d_y^2 G_0 = (e^{2u} - 1)/2.
2. **(cosmetic) Lemma 2.1, empty profile, g >= 2.** log Z_MV has no constant term, so the G^(g)_0 obtained from F lacks <Omega>_g. It lies in w^(2g-2) Q[w] only up to an additive constant; this is harmless because only derivatives occur.
3. **(cosmetic) Statement of Theorem L.** "Every multiset R" includes R = empty (D = 0), which the induction (D = 1, ..., 2h-3) does not cover. The empty profile is also not part of Lemma A (round-1 brief: J >= 1). Say "every non-empty multiset".
4. **(cosmetic) Genus-7 engine run.** The report says the exact genus-7 Hodge-engine evaluation did not finish. R3_toda\out\g7_points_A.out contains all four points (Xi_7 at (13), (14), (2,13), (2,14); finished 19:45). With them, Theorem M at G = 7, R = (13) and (14), both parities, holds exactly as identities in tau (referee2\out\g7_exact.out). Update s.7 and the obstructions.
5. **(minor) Proof (ii) needs P_k for every k <= 2G-2.** The caveat "form fixed by direct checks for k <= 4" can go. Cite the Ueno-Takasaki identity with n' = n+1 at tbar = tbar' and give the derivation of s.4 (right side tau_{n+1}(t) tau_n(t'), orientation from the u^0 term). I checked k <= 6 in bilinear form.
6. **(minor) Dependency bookkeeping.**
   - Proof (i) via Theorem B uses only the MV formula, its cut-and-join equation, the single 2D Toda equation P_2, and Route A's dilaton calculus.
   - Proof (i) via Cor. E additionally uses KP for the MV series (Theorem D) and Takasaki-Takebe's dispersionless Hirota equation.
   
   Record this in the summary. TLR, Theorem F, the closed form and the parity argument are indeed not needed.
7. **(cosmetic) Citations.** I could not confirm "Okounkov math/0004128 s.2.7". Cite the content-product theorem precisely (Orlov-Scherbin; Harnad-Orlov, Comm. Math. Phys. 338 (2015), general form).
8. **(cosmetic) Polynomial-identity checks.** "Polynomial identity in tau, checked at 2G+3 points" is correct, but say why. The relation has degree at most 2G + eps + 2 and an overall factor tau, so 2G+3 non-zero points suffice for both parities.

## 6. Status recommendations
I agree with the author's statuses:
- Theorem L: proved-modulo-literature.
- Theorem M: proved-modulo-literature.
- Corollary I_2: proved.
- Theorem M_P / Corollary I_P: proved-modulo-literature, given the explicit bilinear identity.

For the mKP equation, "computer-verified" understates its standing. It is a literature theorem, and its normalisation is now verified in bilinear form through x^120.

## 7. Files
All in HVC\lemmaA\round2\R3_toda\referee2.

| script | output (in out\) |
|---|---|
| r2_validate.pl | validate_g7_p1.out |
| r2_rel_tables.pl | rel_tables_P2_g0to3.out, rel_tables_P2_g1to3_r3.out, rel_tables_P2_g4to6.out, rel_tables_P3_g0to3.out, rel_tables_P4_g0to2.out, rel_tables_P5_g0to2.out |
| r2_mkp_bilinear.pl | bilinear_P2_d6_x80.out, bilinear_P2_d8_x120.out, bilinear_P3..P6_d5_x60.out |
| r2_highg.pl | highg_G7_p2.out, highg_G8_p1.out, highg_G9_p2.out, highg_G10_p1.out, highg_G3G4_r3.out |
| r2_g7_exact.pl | g7_exact.out |
| r2_thmB_highg.pl | thmB_G7.out, thmB_G8.out, thmB_G9.out, thmB_G10.out |
| r2_degree.pl | degree_bound_g1to8.out |
| r2_logic.pl | logic_P2_G2to9.out, logic_P2_G10to12.out, logic_P3_G2to7.out, logic_P3_G8to12_r2.out, logic_P4_G2to6.out |

Two logic runs were still going when I stopped: logic_P2_G13to30_r3.out (reached G = 21, 0 survivors) and logic_P2_G40_r2.out.

Runtimes: the G = 10 MV engine run took 2070 s (909470 conjugate pairs of partitions, Smax = 52, x-order 103).

Usage examples:
```bash
perl r2_highg.pl 1000000007 1010 10 "2:10:19;2:10:20;3:10:19;2:10:9,11" "18;2,17;9,10;17"
perl r2_mkp_bilinear.pl 998244353 123456789 2 8 120
perl r2_logic.pl 2 2 9
perl r2_thmB_highg.pl 1000000007 1001 10 "19"
```
