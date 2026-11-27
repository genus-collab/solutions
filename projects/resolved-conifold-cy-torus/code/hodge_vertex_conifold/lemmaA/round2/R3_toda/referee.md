# Lemma A round 2: referee R3 (saved by the main session from the workflow journal, 2026-10-31)

Subagents cannot write report files, so this file holds the agent's structured return value: the fields first, then the full report text (report_text) verbatim.

## confirmed

- Literature input (I2)-(I3), re-derived from the fermionic bilinear identity. Charges (1,0) give the DJKM/Jimbo-Miwa mKP equations and charges (0,1) give the author's z^{-1} identity, whose u_k-part is P_k. The orientation is tau_1 = Z(e^{+tau x} o p).
- Generic exact test with random B and random second times, degree <= 9: (D_1^2+D_2) tau_0.tau_1 = 0, JM2 tau_1.tau_0 = 0 and P_2..P_6 tau_0.tau_1 = 0. Controls with the pair swapped or with tau_{-1} fail. KP holds for each tau_n separately.
- The author's Schur-side P_2 check reproduces at tau = 5/11.
- Theorem M re-derived by hand: block leading-coefficient rule, single exponential sector, w-order >= n+r. The unstable genus-0 terms are included, so the multi-sector problem does not arise. Hand expansion at G=2, R=(3) vanishes.
- Empty-profile table (Lemma 2.3): all 50 entries with v >= 1 (g = 0, 1; k <= 4; v <= 5) agree with operator-derived values. The v = 0 entries are more singular but never occur.
- The report's Theorem M formula equals the general leading-order expansion of (D_1^2+D_2) tau_0.tau_1 identically on random inputs: 104 relations.
- Theorem M holds as exact identities in Q[tau] on Route E data: G=0-2 (78), G=3-4 with r<=4 (164), G=5 with r<=3 (94), G=6 with r<=2, J<=13 (86), and G=7 exact with R=(13),(14) (4, with a perturbation control). Total 426, 0 failures.
- The genus-7 engine values agree mod 998244353 with the MV-Schur-side values of rmv.pm.
- Theorem M_P leading orders are exact in Q[tau]: P_3 at G=2-5 (96), P_4 at G=2-4 (60), P_5 at G=2 (16). 0 failures.
- Corollary I_2 re-derived by hand and checked mechanically: every term classified for G=2..11 (3545 cases, 37.7 million terms). 0 survivors and exactly one target, with coefficient -2 Pi_2 sigma^G tau delta.
- Corollary I_P checked mechanically: P_3 at G=2..8 (365 cases), P_4 at G=2..6 (47), P_5 at G=2..5 (9). 0 survivors; the target factor is -(k-1).
- Theorem L induction: base case h=2 needs no hypothesis, (H1)/(H2) are used correctly, there is no circularity, and only n=2G is needed.
- Proof (i) via Theorem B (double induction on r_T and min T, with unit parts treated by the string equation) and via Cor. E at T=(D,2), whose factor (D+1)Pi_{D+1}/(Pi_D Pi_2) is nonzero. Proof (ii) via P_k with c_k = 1-1/k.
- Route A's Theorem B re-implemented: 105 relations (genus 2-5, r <= 3, three new tau), 0 failures.
- lambda_g limit identity re-derived for all G.
- Author scripts reproduce at new parameters: Hurwitz limit G=1-3 (48), line (2) checks (Q=7/4, tau=3/5), and genus-7 mod-p with a new prime/seed and new targets.
- TR audit clean: no topological recursion anywhere in R3_toda or in the referee checks.

## issues

### [1]

#### description

For odd parity n = 2G+1 the relation has degree <= n+2 = 2G+3 in tau, not 2G+2. Example: G = 2, R = (9), ODD: the (a)-term g = 0, k = 6 equals -84(2 tau+1) tau^6, of degree 7. With 2G+3 sample points (e.g. the usage example 'perl mkp_lead.pl 5 5 1 14 1,2,...,13') the odd half of the 182 relations is not proved by the sampling. The even relations are proved, and only those enter the induction.

#### location

report s.7 (row 'Theorem M as a polynomial identity in tau') and summary ('2G+3 points, so these cases are proofs')

#### severity

minor

#### suggested_fix

State the degree bound n+2 and use n+3 points, or use exact polynomial arithmetic. The referee's ref_hirota_lead.pl checks 426 relations of both parities exactly in Q[tau] (genus 0-7).

### [2]

#### description

(a) G^(g)_0 as defined from log Z has no y^0 term. It equals Route A's G_empty minus the constant <Omega>_g, so it lies in w^{2g-2}Q[w] + Q, not w^{2g-2}Q[w]. (b) Undifferentiated empty-profile functions are more singular than their x-order: D G^(0)_0 = 1/2 - w^{-2}/2, D^2 G^(0)_0 = 1 - 1/w, D^3 G^(0)_0 = w - 1. Both are harmless, because a block with U empty consists only of d_{p_1}'s, so v = |B| >= 1. But this is a needed step of the proofs of Theorems M and M_P, and it is not stated.

#### location

s.2, Lemma 2.1 (empty profile, g >= 2) and Lemma 2.3; s.4/s.5 (blocks with U empty)

#### severity

minor

#### suggested_fix

Write 'up to an additive constant' in Lemma 2.1. Add 'a block with U empty has v = |B| >= 1, so only the v >= 1 rows of Lemma 2.3 occur'. Cite Route A Prop. 2.1 at r = 0 explicitly.

### [3]

#### description

'its partner is Phi^0_1(0)': the partner of a genus-G factor Phi^G_1(T_i) is Phi^0_1(R - T_i), and its argument may be nonempty (Phi^0_1(T) != 0 because delta_0(T) >= 3). The conclusion is unaffected, since T_i is a subset of R with level <= D-1, so (H2) applies.

#### location

s.5, proof of Corollary I_2, (d) step 4

#### severity

cosmetic

#### suggested_fix

Replace by 'its partner is Phi^0_1(R - T_i), with k_i = 1 and T_i a subset of R'.

### [4]

#### description

The products (d) also contain genus-G terms Phi^G_1(T) Phi^0_1(R - T), with T a subset of R and level <= D-1. They vanish by (H2), as s.5 says, but the summary sentence omits them.

#### location

summary ('Its only genus-G terms are ... and ...')

#### severity

cosmetic

#### suggested_fix

Add 'and the products Phi^G_1(T) Phi^0_1(R-T), T a subset of R, all of level <= D-1'.

### [5]

#### description

The genus-7 engine run finished. R3_toda\out\g7_points_A.out holds exact v-polynomials of Xi_7(13), Xi_7(14), Xi_7(2,13) and Xi_7(2,14). With them Theorem M holds exactly at G = 7, R = (13) and (14), both parities, as identities in Q[tau]. A perturbation control fails as it should. Xi_7(14) and Xi_7(2,14) agree mod 998244353 with rmv's MV-side values.

#### location

s.7 last paragraph; obstructions list (genus 7-9 checks only mod p)

#### severity

cosmetic

#### suggested_fix

Update s.7 and the obstruction list. Cite referee\out\thmM_poly_g7_exact.out.

### [6]

#### description

What the checks verify is the bosonised Eq (29 checks) and the genus-0 empty-profile equation through y^8. The conclusion that the p'-equation 'admits no termwise leading-order analysis' is a structural judgement that no check verifies.

#### location

claim 'Line (2)' (status computer-verified) and s.8(2)

#### severity

cosmetic

#### suggested_fix

Mark the formulas computer-verified and the no-go statement as a remark. It is not needed.

### [7]

#### description

The form of the bilinear identity follows in three lines from the fermionic bilinear identity on |0>(x)|1>, because contents shift by +n in charge n. That fixes it for every k with tau_1 = Z(e^{+tau x} o p). On |1>(x)|0> the same identity gives Jimbo-Miwa's mKP equations. The referee verified the members P_2..P_6 and JM2 exactly for generic content-product tau functions.

#### location

s.3 'Members used' / obstructions ('form fixed by direct checks for k <= 4')

#### severity

cosmetic

#### suggested_fix

Add the derivation (referee report s.1.1) to s.3.

### [8]

#### description

This step also uses Lemma C and Theorem R of leaf_structure_proof.md s.13. They were not part of R3 and were not re-examined.

#### location

Theorem L consequences ('Corollary R ... follows in every genus')

#### severity

cosmetic

#### suggested_fix

State it as 'with Lemma C and Theorem R (leaf_structure_proof.md s.13), as in the brief'.


## notes

The verdict is sound. Lemma A for all h >= 2 is correctly proved, modulo the MV formula and the first-time-set (mKP) equations of the 2D Toda hierarchy for content-product tau functions. The brief allows both, and I verified their exact form generically.

All referee files are in G:\My Drive\Projects\Old\refGW\hodge_vertex_conifold\lemmaA\round2\R3_toda\referee\, with outputs in out\:
- ref_toda_generic.pl: generic mKP members;
- ref_hirota_lead.pl: exact Q[tau] leading orders of any member, and the random-input transcription test;
- ref_count.pl: mechanical vanishing argument;
- ref_thmB.pl: Theorem B;
- ref_lib.pm and ref_data.pm: libraries.

I used the author's finished genus-7 engine values (R3_toda\out\g7_points_A.out) for an exact genus-7 check. My own runs have all finished; the perl processes still running on the machine belong to other agents.

Mathematica: none needed and none run.

## refuted


## route

R3 (round 2): all-genus Lemma A from the first-time-set (modified-KP) equations of the 2D Toda hierarchy of the Marino-Vafa tau function at the dilaton point (G:\My Drive\Projects\Old\refGW\hodge_vertex_conifold\lemmaA\round2\R3_toda)

## verdict

sound

## report_text

# Referee report on R3, round 2: Lemma A for all genera from the modified-KP equations of the 2D Toda hierarchy

**Where things are**
- Referee scripts: `G:\My Drive\Projects\Old\refGW\hodge_vertex_conifold\lemmaA\round2\R3_toda\referee\`. Outputs are in `referee\out\`. Below, HVC = `G:\My Drive\Projects\Old\refGW\hodge_vertex_conifold`.
- Arithmetic is exact (Perl Math::BigRat) unless marked mod p.
- No wolframscript was run.
- I read and used no TR file, TR script or TR result: no spectral curve, loop equation, BKMP, dictionary P2 or mod-P2 result.

## Verdict: sound

I tried to refute every claim marked proved, proved-modulo-literature or computer-verified. None was refuted.

Theorem L holds as stated: Lemma A = M_irr(h) for every h >= 2. Hence (OL) holds for every G >= 3 and every D <= 2G-3, odd and even. The proof depends on:
- the MV formula in Route A's normalisation (1.1), refereed with Route A;
- the first-time-set (modified-KP) equations of the 2D Toda lattice for content-product tau functions at fixed second times. The brief allows this input explicitly. I re-derived the exact form from the fermionic bilinear identity and checked it exactly for generic content-product tau functions (s.1);
- Route A's proved Prop. 2.1, Cor. 2.2 and Lemma 3.1 (the dilaton-point calculus);
- for proof (i) only, Route A's Theorem B (cut-and-join) or Cor. E (KP). Proof (ii) needs neither.

All issues found are minor or cosmetic (s.6), and none affects the proof:
- two imprecise sentences in the vanishing argument;
- a point count that does not prove the odd-parity sampling checks;
- an unstated, but true, step: empty-profile blocks always carry d_y's;
- a stale remark about the genus-7 run.

## 1. The literature input (I2)-(I3): confirmed

### 1.1 Form and orientation, re-derived

**Setup.**
- Take G = e^{beta F_2} Gamma_-(p'), with beta = tau x. Here F_2 is the diagonal fermion bilinear with eigenvalue the content sum, and Gamma_- is a vertex operator. So G lies in the GL(infinity) group.
- In the charge-n sector every content is shifted by +n. Hence tau_n(p) = <n|e^{H(t)} G|n> = c_n Sum_nu e^{beta(f_2(nu) + n|nu|)} s_nu(p') s_nu(p) = c_n Z_MV(e^{n tau x} o p) at p' = q^rho.
- The constant c_n (the vacuum energy) drops out of every bilinear equation.

**The bilinear identity.** Apply Sum_i psi_i G (x) psi*_i G = Sum_i G psi_i (x) G psi*_i to |n> (x) |m> and pair with <n+1|e^{H(t)} (x) <m-1|e^{H(t')}. The right side is a sum over n <= i < m.
- **(n,m) = (1,0).** The right side is 0. This gives the DJKM 1-step mKP identity Res_z z tau_1(t-[1/z]) tau_0(t'+[1/z]) e^{xi(t-t',z)} dz = 0. Put t = T+u, t' = T-u. The u^0 part is Jimbo-Miwa's (D_1^2 - D_2) tau_1.tau_0 = 0. The u_1 part is (D_1^3 + 3D_1D_2 - 4D_3) tau_1.tau_0 = 0.
- **(n,m) = (0,1).** The right side is tau_1(t) tau_0(t'). This is exactly the author's identity Res_z z^{-1} tau_0(T+u-[1/z]) tau_1(T-u+[1/z]) e^{2xi(u,z)} dz = tau_1(T+u) tau_0(T-u). Its u^0 part is trivial, and its u_k-linear part is P_k = D_{t_k} + p_k(-D~).

So 2P_2 = D_1^2 + D_2 on tau_0.tau_1 is Jimbo-Miwa's first mKP equation, and P_3 is equivalent to Jimbo-Miwa's second. The two forms the author uses are consistent with each other, with tau_1 = Z(e^{+tau x} o p). The derivation fixes the form for every k, not only for k <= 4.

### 1.2 Exact generic test (`ref_toda_generic.pl`)

The script uses its own Murnaghan-Nakayama characters and its own Hirota operators. It takes tau_n = Sum_{|nu| <= d} B^{f_2 + n|nu|} s_nu(p') s_nu(p) with random rational B and random rational second times p'_1..p'_d, and tests every monomial of degree <= d minus the weight.

| run | members tested | result |
|---|---|---|
| d = 6, B = 3/2 | (D_1^2+D_2) tau_0.tau_1, (D_1^3+3D_1D_2-4D_3) tau_1.tau_0, P_2..P_4 | all exactly 0 |
| d = 9, B = -5/7 | the same, with P_2..P_6 | all exactly 0 |

Controls:
- Swapping tau_0 and tau_1 gives nonzero results in every case, as does replacing tau_1 by tau_{-1}.
- (D_1^2+D_2) tau_{-1}.tau_0 = 0, as it must, since that is the next adjacent pair.
- KP holds for tau_0 and for tau_1 separately (a sanity check of the code).

Every coefficient is a polynomial in B^{+-1} and p'. So the MV case B = e^{tau x}, p' = q^rho follows as formal Laurent series in x.

The author's Schur-side script also reproduces at a new value, tau = 5/11: 7 monomials through x^14 vanish, and the wrong sign fails in all 7 (`author_mv_mkp_direct_rerun_t511.out`).

## 2. Theorem M: confirmed

### 2.1 Re-derivation by hand

**Operators.** On functions of w:
- D = y d_y = (w-1) w^2 d_w, and d_y = e^u w^3 d_w.
- Hence e^{-uJ}(J+D)e^{uJ} w^m = (J-m) w^{m+1} + m w^{m+2}.
- And e^{-u(J+1)} d_y e^{uJ} w^m = J w^{m+1} + m w^{m+2}.

Both operators raise the order by at least 1, so leading coefficients multiply.

**Block rule.** A Hirota block B has U = T_B + {i in B : i >= 2} and v = #{i in B : i = 1}. At x-power 2g-2+|B|+k it contributes:
- the factor sign x prod_{i in B} i x [2 if B is even and k = 0] x prod_{i >= 2} Pi_i x sigma^{|B|+g-1} tau^k/k!;
- times the leading coefficient (J_U)^(v) (delta_g(U))_k X_g(U);
- at w-order (x-power) + r_{T_B}.

So every term of the q_R coefficient at x^n lies in the single sector e^{u(J_R + weight)} and has w-order >= n + r. This is the rule of the report's s.5; for P = D_1^2 + D_2 it gives Theorem M.

**Points checked explicitly.**
- F = log Z contains the unstable genus-0 terms. Therefore G^(0)_0 = u - 3u^2/4 + u^3/6 and d_y^2 G^(0)_0 = e^{2u}/2, a single sector. The multi-sector problem noted by Route A's referee arises only when these terms are dropped.
- The log in G^(1)_0 never occurs undifferentiated.
- At fixed (x^n, q_R) the equation is an identity in Q[[y]]. With u = T(y) it says that e^{u(J+2)} times a rational function of u vanishes, so every w-coefficient vanishes. The argument is purely formal.
- Hand expansion at G = 2, R = (3), even parity: eight nonzero terms, sum exactly 0 at tau = 1, 2/7 and -5/3.

### 2.2 Empty-profile table (Lemma 2.3)

I recomputed the table by applying D and d_y to the closed forms of G^(0)_0, D G^(1)_0 and d_y G^(1)_0 (`lead0_table_check.out`).
- All 50 entries with g = 0, 1, k <= 4 and 1 <= v <= 5 agree with the report's table.
- For v = 0 these functions are more singular than their x-order: D G^(0)_0 = 1/2 - w^{-2}/2, D^2 G^(0)_0 = 1 - 1/w, D^3 G^(0)_0 = w - 1. So the restriction v >= 1 is essential.
- It always holds, because a block with U empty consists only of d_{p_1}'s, so v = |B| >= 1 (issue 4).

### 2.3 Transcription

`ref_hirota_lead.pl` is a from-scratch general leading-order implementation: set partitions, sub-multiset assignments, the (g,k) budget, and operator-derived empty-profile values. I replaced every sigma^g X_g(T) by a random rational. The general expression for D_1^2 + D_2 then equals the report's Theorem M formula identically: 104 relations at G = 0..3, of which 98 are nonzero, with 0 differences (`thmM_random_transcription_g0to3.out`).

### 2.4 Exact checks as identities in Q[tau]

These are not sampled at points of tau, so each checked case is a proof of that case. Data come from Route E's exact tables (read with my own parser) and Route E's point engine.

| genus | range | relations | new? | output file |
|---|---|---|---|---|
| G = 0-2 | r <= 2, J <= 7 | 78 | | `thmM_poly_g0to2_r2_J7.out` |
| G = 3-4 | r <= 4, J <= 10 (R+2 with up to 5 parts) | 164 | yes | `thmM_poly_g3g4_r4_J10.out` |
| G = 5 | r <= 3, J <= 11 (R+2 with 4 parts) | 94 | yes | `thmM_poly_g5_r3_J11.out` |
| G = 6 | r <= 2, J <= 13 | 86 | | `thmM_poly_g6_r2_J13.out` |
| G = 7 | R = (13), (14), both parities, exact | 4 | yes | `thmM_poly_g7_exact.out` |

Total: 426 exact polynomial identities, 0 failures.

**Genus 7.** The genus-7 values come from the author's finished engine run `R3_toda\out\g7_points_A.out` (see issue 5).
- Perturbing the constant term of Xi_7(2,13) breaks the even relation (`thmM_poly_g7_perturbed_control.out`).
- The engine values Xi_7(14) and Xi_7(2,14) reduce mod 998244353 to exactly the MV-Schur-side values computed by `rmv.pm`.

### 2.5 General members, exact in Q[tau]

These runs use my implementation with operator-derived empty-profile values for v up to 5.

| member | genus | range | relations | failures |
|---|---|---|---|---|
| P_3 | G = 2-5 | r <= 1, J <= 12 | 96 | 0 |
| P_4 | G = 2-4 | r <= 1, J <= 10 | 60 | 0 |
| P_5 | G = 2 | r <= 1, J <= 8 | 16 | 0 |

The author checked M_3 only up to genus 4 and M_4 only up to genus 3, so P_3 at genus 5, P_4 at genus 4 and all of P_5 are new.

## 3. Corollaries I_2 and I_P (the vanishing argument): confirmed

**By hand.** I re-derived:
- the falling-factorial dichotomy;
- the count Sum(delta_i - k_i) = delta_G(R) <= -2, where an empty argument has delta := 2 and nonzero empty factors have delta - k in {0, 1};
- the budget argument for general members: every block contributes 2g-2+|B|+k >= 0, so a genus-G block is (1,1) or (2,0) and has level <= D, with equality only for the linear term;
- the target coefficients: -2 Pi_2 sigma^G tau delta, and -(k-1) sigma^G tau Pi_k delta for P_k (c_k = 1 - 1/k).

**Mechanically (`ref_count.pl`).** I enumerated every term of the leading-order expansion, for all D <= 2G-3 and all R at the required level. Every non-target term is killed by one of: a vanishing falling factorial, (H1), (H2), or a zero empty-profile coefficient. The target occurs exactly once, with a nonzero coefficient.

| equation | genus range | (G,D,R) cases | terms | survivors | target factor -ck |
|---|---|---|---|---|---|
| D_1^2+D_2 (Theorem M) | G = 2..11 | 3545 | 37.7 million | 0 | -2 |
| P_3 | G = 2..8 | 365 | 49.0 million | 0 | -2 |
| P_4 | G = 2..6 | 47 | 11.4 million | 0 | -3 |
| P_5 | G = 2..5 | 9 | 3.6 million | 0 | -4 |

Two wording slips (issues 2 and 3) do not affect the argument.

## 4. Theorem L, the induction: confirmed

**Structure.** The induction runs on h, and then on D = 1, ..., 2h-3.
- Only the even relations n = 2G are used.
- No circularity: Theorem M and Theorem B are derived without Lemma A, Cor. I_2 uses only (H1) and (H2), and Cor. E uses exactly (H1) and (H2).

**Base case.** h = 2 needs no hypothesis. At n = 4 and R = empty, Theorem M reduces to -2 Pi_2 sigma^2 tau (-1) X_2(2) = 0:
- the lower-genus terms die by (1)_3 = (3)_5 = 0;
- B_2(0) = 0;
- products of empty Phi's have k_1 + k_2 <= 4 < 6.

At level D = 1 in general, only (H1) for Xi_g(2), g < G, is used.

**Proof (i), via Theorem B.** Under (H1), (H2) and I_2, Theorem B at R (level D-1) reduces to Sum_i gamma_{j_i} X(R_i^+) = Sum_{i<k} gamma_{j_i j_k} X(R_ik). I checked every lower-order term, including unit parts:
- the string equation applies, and Xi_g(1) = 0 for g >= 2;
- the genus-(h-1) three-leaf term has level D-2;
- Xi_1(1) != 0 occurs only with a genus-(h-1) partner of level D-2 <= 2(h-1)-3;
- genus-0 cofactors force the genus-h partner to level <= D-2.

The double induction on r_T and min T is correct. The case r_T = 1 uses R = (D) and needs D >= 2; D = 1 is covered by I_2.

Alternatively, Cor. E at T = (D,2) has the nonzero factor (D+1) Pi_{D+1}/(Pi_D Pi_2). So I_2 at R = (D) kills Xi_G(D+1), and Cor. E then kills all of level D.

**Proof (ii), via P_k.** Any member with nonzero D_{t_k} coefficient works.

**Theorem B, re-implemented (`ref_thmB.pl`)**, each relation at tau = 5/11, -7/3 and 4/9:
- genus 2-3, r <= 2, J <= 7: 24 relations;
- genus 3-5, r <= 3, J <= 9: 81 relations. The three-part R are new.

0 failures in both runs.

## 5. Other claims

- **lambda_g limit:** re-derived by hand. Both terms equal 2 b_G J(J+1) prod_{m=1}^{2G-4+r}(J-m)(J+3-2G-r), and the other terms are O(tau^2).
- **Author scripts reproduce at new parameters:**
  - `mkp_hurwitz_limit.pl`, G = 1-3: 48 relations;
  - `eq_bosonization_check.pl` at Q = 7/4 and `toda_eq_genus0_check.pl` at tau = 3/5;
  - `g7_mkp_modp.pl` with prime 998244353, seed 31 and targets (14), (16): 4 relations, and the control fails.
- **Counts in R3_toda\out:** every count quoted in the report matches its output file (374, 182, 60, 14, 314, 162, 136, 122, 29).
- **Line (2):** the formulas reproduce. The no-termwise-analysis statement is a remark, not a verified fact (issue 6).
- **Line (3):** conjectural, and not needed.

## 6. Issues (none fatal or major)

1. **(minor) Odd-parity point count.** The report says "degree <= 2G+2, checked at 2G+3 points, so these cases are proved". For odd parity n = 2G+1 the degree bound is n+2 = 2G+3. Example: at G = 2, R = (9), odd parity, the (a)-term with g = 0, k = 6 equals -84(2 tau+1) tau^6, of degree 7. So the odd half of the 182 "proved" relations is not proved by the sampling. The even relations are proved, and only those enter the induction. My exact Q[tau] checks (s.2.4) cover both parities.
2. **(cosmetic) Cor. I_2, step (d)4.** The partner of a genus-G factor Phi^G_1(T_i) is Phi^0_1(R - T_i), and its argument may be nonempty: Phi^0_1(T) != 0 because delta_0 >= 3. The conclusion still holds, since T_i is a subset of R and so has level <= D-1, and (H2) applies.
3. **(cosmetic) Summary sentence.** "Its only genus-G terms are ..." omits the products Phi^G_1(T) Phi^0_1(R-T), with T a subset of R and level <= D-1. These vanish by (H2), as the proof itself says.
4. **(minor) Lemmas 2.1 and 2.3.**
   - G^(g)_0 for g >= 2, as defined from log Z, equals Route A's G_empty minus <Omega>_g. So it lies in w^{2g-2} Q[w] + Q, not in w^{2g-2} Q[w].
   - With v = 0 the empty-profile functions are more singular than their x-order (s.2.2).
   - Both points are harmless because v = |B| >= 1 for every empty-U block. This step should be stated, and Route A's Prop. 2.1 at r = 0 should be cited.
5. **(cosmetic) The genus-7 run did finish.** `R3_toda\out\g7_points_A.out` holds exact Xi_7(13), Xi_7(14), Xi_7(2,13) and Xi_7(2,14). The first point took 3674 s; the other three took seconds. These values now give an exact genus-7 check (s.2.4).
6. **(cosmetic) Status of the line (2) claim.** What is computer-verified is the bosonised Eq and the genus-0 equation. That the p'-equation "admits no termwise leading-order analysis" is an interpretation; mark it as a remark.
7. **(cosmetic) Corollary R.** It also uses Lemma C and Theorem R (leaf_structure_proof.md s.13). These were not part of R3 and were not re-examined here.
8. **(cosmetic) Form of the bilinear identity.** It is fixed "by direct checks for k <= 4". The three-line fermionic derivation of s.1.1 fixes it for all k and should be added to the report's s.3.

## 7. Obstructions and a simpler path

**Obstructions.**
- The report's statement "no obstruction for Lemma A, modulo the cited literature" is correct.
- The external inputs are standard, and their exact form was verified generically.
- Lines (2) and (3) are side remarks, not obstructions to the result.

**Simplest complete proof.** Theorem M (only n = 2G is needed), then Cor. I_2, then Cor. E at T = (D,2).

**Most self-contained proof.** Proof (ii): only the 2D Toda members P_k, MV, and Route A's dilaton calculus.

**No longer needed.** TLR(G), Theorem F, the closed form D5 and the parity argument; also the odd relations and the Hurwitz and lambda_g limits.

## 8. TR audit: clean

- The author's code imports only Route A (`mvlib`, and `rmv.pm` from Route A's referee), Route B (`Rec.pm`), Route E (`Hodge.pm`, `HodgePt.pm`, `WK.pm`) and the author's own modules.
- A grep finds no TR-related import or script.
- `R3_thmM_check.wls` was not run.

## 9. Referee files

All in `HVC\lemmaA\round2\R3_toda\referee\`.

**Scripts:**
- `ref_lib.pm`: characters, power-sum polynomials, Hirota operators.
- `ref_toda_generic.pl`: the generic test of s.1.2.
- `ref_hirota_lead.pl`: modes `data` and `random`; members M, P2-P6 and JM2.
- `ref_count.pl`: the mechanical vanishing check of s.3.
- `ref_data.pm`: table parser and engine access.
- `ref_thmB.pl`: Theorem B.

**Outputs, in `out\`:**
- generic test: `toda_generic_d6_B32_s1.out`, `toda_generic_d9_Bm57_s2.out`;
- Theorem M: `thmM_random_transcription_g0to3.out`, `thmM_poly_g0to2_r2_J7.out`, `thmM_poly_g3g4_r4_J10.out`, `thmM_poly_g5_r3_J11.out`, `thmM_poly_g6_r2_J13.out`, `thmM_poly_g7_exact.out`, `thmM_poly_g7_perturbed_control.out`, `g7_points_perturbed.txt`;
- empty-profile table: `lead0_table_check.out`;
- general members: `P3_poly_g2to5_r1_J12.out`, `P4_poly_g2to4_r1_J10.out`, `P5_poly_g2_r1_J8.out`;
- vanishing check: `count_M_g2to11.out`, `count_P3_g2to8.out`, `count_P4_g2to6.out`, `count_P5_g2to5.out`;
- Theorem B: `thmB_ref_g2g3_r2_J7.out`, `thmB_ref_g3to5_r3_J9.out`;
- reruns of the author's scripts: `author_mv_mkp_direct_rerun_t511.out`, `author_hurwitz_limit_rerun_g1to3.out`, `author_line2_rerun.out`, `author_g7_modp_rerun_seed31.out`.

