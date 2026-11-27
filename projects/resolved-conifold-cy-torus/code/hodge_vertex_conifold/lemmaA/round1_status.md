> **Main-session note (2026-10-31, added when this file was saved).** The synthesis below could not see two rulings that the user has since made explicit:
> 1. The KP hierarchy of the Mariño–Vafa generating function is an allowed input.
> 2. The TR restriction applies only to the resolved conifold × C² CohFT. TR considerations applied to the Hodge CohFT are allowed: the framed C³ curve C_τ, BKMP for C³, the dictionary P2, Theorem 3.1 of proof_M_hurwitz.md, Lemma D and the TR scripts.
>
> So the conflict recorded in §5.1(b) is resolved, and round 2b (round2\TR1_dictionary, round2\TR2_recursion) is within the rules. No TR is used, or allowed, for the conifold × C² theory itself. Everything else below is the synthesis output, saved verbatim.

# Lemma A (M_irr(h) for all h): round-1 status

Synthesis of round 1 (routes A-E and their referees), 2026-10-31. Intended file: HVC\lemmaA\round1_status.md, where HVC = G:\My Drive\Projects\Old\refGW\hodge_vertex_conifold. The main session saves this text, because subagents cannot write report files.

**Governing rule.** This workflow was started by the request "no TR at all please". Nothing below uses topological recursion: no spectral curve, Eynard-Orantin invariant, loop equation, BKMP/remodelling, dictionary P2, Lemma D, "mod P2" result or TR script, not even as a cross-check. Section 5.1(b) records a conflict with round2\brief.md.

**What I read.**
- lemmaA\brief.md.
- report.md and referee.md of A_mv, B_hurwitz, C_intersection and D_localization.
- E_formulas\report.md and the Route E referee result. E_formulas\referee.md has not been saved yet.
- lemmaA\parity\README.txt.
- leaf_proof\leaf_structure_proof.md s.0, 13 and 14.
- conifold_regularity_note.md (preamble, s.1-3a).
- proof_status_mixed_scale.md s.1-2.
- Brini-Schuler, arXiv\main.tex lines 1660-1790.
- round2\brief.md, for context only.

No wolframscript was run. The new computations are in s.1b, with the commands that reproduce them.

**Conventions (c = 1).**
- Weights w = (1, tau, -1-tau), and sigma = tau(1+tau).
- Pi_j = prod_{a=1}^{j-1} (j tau + a)/(j-1)!. In particular Pi_2^2 = 1+4sigma, Pi_3 = (9sigma+2)/2, and Pi_4 = (16sigma+3)(2tau+1)/3.
- For a multiset R: level D(R) = Sum (j_i - 1), J = Sum j_i, r = #R, and P_R := sigma^G Xi_G(R).
- Components: Xi = Sum_{k=0}^{h} (-1/sigma)^k Xi^(k). The k = 0 component is the Hurwitz component; k = h is the lambda_g component.
- Z(D) is the number of zeros (at sigma != 0, counted with multiplicity) that Corollary E forces on P_(D+1). It equals the maximal pole order of Pi_{D+1}/prod_t Pi_t over the T at level D, summed over the special framings tau = -k/a. The first values are Z(1..13) = 0, 1, 2, 5, 6, 11, 13, 19, 22, 29, 31, 44, 47 (A_mv\a3_zcount.pl, runs\zcount.out; re-derived here).

## 0. Bottom line

1. Lemma A is still not proved for all h.

2. The following hold for all h. They are TR-free and their referees judged them sound.
   - The Hurwitz component vanishes on the range. There are two independent proofs, so sub-lemma A2 is settled.
   - Xi_h(2) = 0 identically, again with two independent proofs.
   - Route D's conifold polar relations T_k hold. Their input (H) is a theorem (Brini-Schuler Prop. 3.1).
   - Several exact reformulations of Lemma A.

3. Conditional on the KP hierarchy of the Marino-Vafa series, Lemma A holds for every h <= 20 (Route A, Theorem F, computer-assisted). KP is TR-free. round2\brief.md s.0 records it as admitted by the user; please confirm, since this synthesis cannot see that ruling. Without KP, exact computation with two independent engines gives M_irr(h) for every h <= 7.

4. New in this synthesis (s.1b; not yet refereed). Assume Lemma A in genera < G. Then:
   - Lemma A(G) is equivalent to Xi_G(D+1) = 0 at the odd levels D >= 3 with Z(D) <= G-2 only.
   - That is level D = 3 for 4 <= G <= 7; levels 3 and 5 for 8 <= G <= 14; levels 3, 5 and 7 for 15 <= G <= 23. These are the one-leaf values at j = 4, 6, 8, ...
   - At genus G these are polynomials with G-3, G-7, G-14, ... unknown coefficients.
   - The base cases h = 3 and h = 4 no longer need the Mathematica tables.
   - The regularity theorem holds through genus 8 without KP, and through genus 21 with KP.

5. The decisive open target is Xi_G(4) = 0 for all G >= 4, and then the same for j = 6, 8, .... Two routes lead there: KP rigidity at special framings, or the one-leaf formula (6.1).

## 1. Rigorously proved (only claims that survived the referees)

**1.1 Hurwitz component (sub-lemma A2), all h.** For h >= 1 and every multiset R of integers >= 2 with D(R) <= 2h-3, and also with D(R) = 2h-1:
    Xi^(0)_h(R) = Int_{Mbar_{h,r}} c(T(-log)) prod_i 1/(1 - j_i psi_i) = 0.
- **Route B, Theorem A.**
  - Status: proved modulo ELSV (Invent. Math. 146 (2001)); Okounkov, arXiv:math/0004128 eq. (10); and Hurwitz's genus-0 formula (Goulden-Jackson, PAMS 125 (1997)).
  - Written in B_hurwitz\report.md, summary s.2-5 (Lemma 4.1/4.2, Lemma 5.1, Lemmas 6.1-6.3, counting in s.5).
- **Route C, Theorem H**, an independent proof. It also allows j_i >= 0, with J >= 1 and D <= 2h-3.
  - Status: proved modulo Farkas-Pandharipande (J. Inst. Math. Jussieu 17 (2018), Thm 3 and s.0.1); Fulton (Prop. 1.8 and 14.1, projection formula); Vistoli; and Lemma 2 of proof_M_residue.md.
  - Written in C_intersection\report.md, steps (0)-(6).
- Both referees re-derived every step and computed 27 and 29 further zeros from independent Hurwitz numbers.
- Consequence, together with (G): in the range, P_R lies in span(sigma, ..., sigma^(h-1)).

**1.2 Hurwitz closed forms (Route B, Theorems B and C).**
- Theorem B: Xi^(0)_h(j) = (j-1)...(j-2h+2) [t^(2h)] S^2 H^j, with H = S exp(-tS'/S). This is the Hurwitz limit of (6.1) (proof_M_residue.md s.6a(iii)), now proved for all h.
- Theorem C: the recursion (REC), which determines every Xi^(0) with D != 2h-2 from genus 0.

**1.3 Xi_h(2) = 0 for every h >= 2 and every tau.**
- Route A, Theorem A: proved modulo the Marino-Vafa cut-and-join (Liu-Liu-Zhou, JDG 65 (2003)). Written in A_mv\report.md s.3.2.
- Route D, Corollary D5.1: proved modulo Brini-Schuler and the conifold framework (s.5.2). Written in D_localization\report.md s.5; the referee's three-line proof is in D_localization\referee.md s.3.

**1.4 Reformulations.** These are proved using string and dilaton only, and were derived three times independently (A: Prop. 2.1, Cor. 2.2-2.3; B: Lemma 4.1-4.2, Prop. D; D: Prop. D2.3).
- Laguerre form: Xi_h(R) = Sum_{a=0}^{N} ((-1)^a/a!) L^{(2h-2+r+a)}_{N-a}(-J-a) M_a, with N = 3h-3+r and M_a = <Omega lam_R lam_1^a> (Marino-Vafa integrals with a unit legs).
- Lemma A <=> deg P_R(u) <= N-1, where Q_h(s;R) = (1-u)^(2h-2+r) e^{-uJ} G_R(u e^{-u}) and u = s/(1+s).
- Lemma A <=> e^{-uJ} G_R(u e^{-u}) = O(u^-(2h-1+r)) at u = oo.
- Lemma A <=> Sum_h hbar^(h-1) Xi_h(R) is a Laurent polynomial of degree <= floor(D/2). This uses Route A's Kaehler-twisted framed vertex, Cor. 2.3.
- Also proved: the lambda_g residue formula for every s (A s.2.4).

**1.5 Cut-and-join consequences** (proved modulo the Marino-Vafa cut-and-join).
- A Theorem B: the 2-insertion relation (A s.3.3).
- A Theorem C: at tau = -1/2, Xi_h(T) = 0 for every T of odd parts in the range (A s.3.4). This is a TR-free replacement of the tau = -1/2 statement of proof_status s.1.6.

**1.6 Conifold polar relations (Route D, D_localization\report.md s.5).**
- **Input (H) is a theorem in every genus.** Brini-Schuler, "Refined Gromov-Witten invariants", Prop. 3.1 = Theorem A1 (arXiv\main.tex lines 1704-1772; main_DMJ.tex).
  - Statement: for T inside T_CY with alpha_i = eps_j, GW_{d[P1]} = (1/d)(2 sinh(d eps1/2))^-1 (2 sinh(d eps2/2))^-1.
  - Proof: Li's degeneration formula, Mumford's relation, BP08 Lem. 6.3 and BP05 Thm 5.1; no TR.
  - The project torus (alpha0 = 0, alpha1 = eps1, alpha2 = eps2) lies in its class. I checked the statement and its hypotheses in the tex. The referee confirmed it independently by localization up to (g,d) = (4,1), (3,2), (2,3).
- **The relations.** For 1 <= k <= 2G-3, given Lemma A for levels <= k-2 in all genera <= G:
    T_k(G) = Sum_{D(R)=k} (1/|Aut R|) prod_i ((-1)^{j_i} Pi_{j_i}^2/j_i) Xi_G(R) = 0.
  - Lemma 5.2: the level k-1 terms cancel identically.
  - Lemma 5.3: graphs with an edge do not contribute.
- Special cases:
  - D5.2: 2(2+9sigma)^2 Xi_h(3) = 3(1+4sigma)^2 Xi_h(2,2), for h >= 3.
  - D5.3: 4(3+16sigma)^2 Xi_h(4) - 6(2+9sigma)^2 Xi_h(2,3) + 3(1+4sigma)^2 Xi_h(2,2,2) = 0.
- Status: proved modulo literature. It depends on Brini-Schuler, Theorem LC, Theorems A-C of conifold_regularity_note.md, and the frame-B premise (s.5.2).

**1.7 Route C structure.**
- Prop. L. Its hypothesis must read r_+ >= 2k.
- Prop. L1, proved modulo the strata bound.
- C3(a): single Hodge insertions do not vanish, so the Calabi-Yau relation of the weights is essential.

**1.8 Route E structure** (E_formulas\report.md s.1-2).
- The forgetful recursion with K (Prop. W).
- The binomial-support equivalence, and the excess-Chern form of the Hurwitz coefficients.
- The kappa reduction for r = 2.
- No formula in J alone is compatible with both string and dilaton.
- The top j-degree part is the Witten-Kontsevich r-point function.
- The Mumford decomposition of the k = h-1 class.

**1.9 Route D localization facts.**
- D2.1: the Graber-Pandharipande formula for X (classical).
- D2.2: the sigma-degree and top coefficient of each graph.
- D3.2: the resonance at c = w1.

**1.10 Conditional on KP** (TR-free; A_mv\report.md s.3.4-3.6). The referee judged this sound and reproduced Theorem F with independent code.
- Fact 3.5: the Marino-Vafa series is a KP tau function (GL(oo) orbit: Takasaki 1984; Orlov-Scherbin, TMP 128 (2001)).
- Theorem D: KP at the dilaton point. It needs Lemma A in genera < G.
- Corollary E (merge formula): Xi_G(T) = (J-1)...(J-r+1) Pi_{D+1}/prod_t Pi_t * Xi_G(D+1).
- Theorem C': vanishing at tau = -1/3.
- Theorem F: for 5 <= G <= 20, Lemma A(<G) implies Lemma A(G). The certificates are exact over Q at G = 5, 6 (full system) and G = 7 (two-leaf system), and mod p otherwise, with a lifting argument the referee confirmed.
- Hence Lemma A for every h <= 20.

**1.11 Computer-verified without KP** (exact rationals).
- M_irr(5) (44 points), M_irr(6) (96 points) and M_irr(7) (194 points): every component is zero, and every control is non-zero in every component.
  - Engines: Route E's engine and the referee's independent RefWK/RefHodge (E_formulas\referee\tally_final.txt).
  - The author's own genus-7 run has 191 of 194 points so far, with no violation.
- (6.1) holds identically in (j, tau) for h <= 7.
- Normalisations against the tables: the Marino-Vafa formula (2170 + 996 comparisons) and ELSV (236 + 19).

**1.12 Refuted.**
- Route E's candidate closed forms (K1, Y1, D2, bi-binomial, separable, exponent forms, and the pairwise r = 3 extension).
- Route C's R-independent decomposition (r = 2).
- Route D's claim that (H) is known only through genus 4.

## 1b. Combinations checked in this synthesis (not yet refereed)

**S1. Parity** (main session, lemmaA\parity\README.txt; re-derived here). Substituting Corollary E into T_D gives
    T_D(G) = Pi_{D+1} kappa_D Xi_G(D+1),   kappa_D = [y^D](tau y + log((e^y - 1)/y)).
- So kappa_1 = tau + 1/2; kappa_D = B_D/(D D!) != 0 for even D; and kappa_D = 0 for odd D >= 3.
- Proof checked:
  1. kappa_D = (1/D)[y^D](1-alpha)^-D.
  2. Apply Lagrange inversion.
  3. Use Sum_t (Pi_t/t) z^t = log B_{tau+1}(z), the generalised binomial series.
  4. Then w(y) = e^{tau y}(e^y - 1).
- Exact check for D <= 10 at four values of tau.

**S2. Level shift.**
- In Theorem D, E_0(m,n;S) != 0 forces m, n >= 1 (or m = n = 0 with S empty). So a relation whose main term is at level D contains genus-G values only at level D and at levels <= D-2. Unit parts are included via string.
- Hence Corollary E at level D needs only levels <= D-2 to vanish (A states "< D"), and so does T_D.
- Consequence: given Lemma A in genera < G, all points of levels 1, 2 and 4 vanish at genus G, whatever happens at level 3.
  - Level 1 is 1.3.
  - Levels 2 and 4 follow from S1, with kappa_2 = 1/24 and kappa_4 = -1/2880.
- Level 3 is the first open level, and every level >= 5 depends on it.

**S3. Counting.**
- At an odd level D whose lower levels <= D-2 vanish:
  - P_(D+1) lies in span(sigma, ..., sigma^(G-1)), by 1.1 and (G).
  - P_(D+1) is divisible by a polynomial of degree Z(D) with non-zero constant term, by Corollary E.
- So level D closes by counting alone once Z(D) >= G-1. The top level 2G-3 always closes this way.
- The remaining odd levels are those with Z(D) <= G-2:

| level D | j = D+1 | open from genus |
|---|---|---|
| 3 | 4 | G >= 4 |
| 5 | 6 | G >= 8 |
| 7 | 8 | G >= 15 |
| 9 | 10 | G >= 24 |
| 11 | 12 | G >= 33 |
| 13 | 14 | G >= 49 |

- At an open level, G-1-Z(D) coefficients are free:

| G | free coefficients |
|---|---|
| 4 | 1 |
| 5 | 2 |
| 6 | 3 |
| 7 | 4 |
| 8 | 5+1 |
| 10 | 7+3 |
| 15 | 12+8+1 |
| 20 | 17+13+6 |
| 21 | 18+14+7 |

- So, given Lemma A(<G): Lemma A(G) <=> Xi_G(D+1) = 0 for these D, taken in increasing order. This uses KP and the conifold. Without the conifold, Corollary E alone reduces Lemma A(G) to all one-leaf values Xi_G(j), 2 <= j <= 2G-2 (A s.3.5).

**S4. Base cases without the Hodge tables.**
- **Genus 3, by hand.**
  - Level 1 is 1.3.
  - Level 2: Corollary E gives Xi_3(2,2) = (3/2)(9sigma+2)/(1+4sigma) Xi_3(3). Then D5.2 becomes -(1/2)(9sigma+2) Xi_3(3) = 0.
  - Level 3: Xi_3(2,3) = (8/3)(16sigma+3)/(9sigma+2) Xi_3(4), and Xi_3(2,2,2) = (20/3)(16sigma+3)/(1+4sigma) Xi_3(4). So P_(4) is divisible by (9sigma+2)(4sigma+1) and lies in span(sigma, sigma^2), which forces P_(4) = 0.
  - Inputs: KP, T_2, 1.1, (G) and M_irr(2).
- **Genus 4, by rank.** The KP relations of Theorem D, plus Xi_4(2) = 0, plus the Hurwitz theorem (the sigma^4 coefficient of all 18 multisets T), have full column rank 72: 498 equations, mod p = 2^31-1, at 28 random integer tau.
  - KP + cut-and-join + Hurwitz also has rank 72.
  - Without the Hurwitz theorem, Route A had nullity 1.
  - At genus 3, KP + cut-and-join + Hurwitz leaves nullity 1, in (3) and (2,2); T_2 removes it.
- So Lemma A for h <= 20 follows from TR-free theorems plus exact rank computations, with no Mathematica data. The genus-3 and genus-4 tables remain an independent check.
- Commands:
```bash
cd "/g/My Drive/Projects/Old/refGW/hodge_vertex_conifold/lemmaA/A_mv"
H4="H=2 H=3 H=2.2 H=4 H=3.2 H=2.2.2 H=5 H=4.2 H=3.3 H=3.2.2 H=2.2.2.2 H=6 H=5.2 H=4.3 H=4.2.2 H=3.3.2 H=3.2.2.2 H=2.2.2.2.2"
perl a3_kp_nullity.pl 4 Z=2 $H4    # rank 72, nullity 0
perl a3_kp_nullity.pl 4 C $H4      # rank 72, nullity 0
perl a3_kp_nullity.pl 3 C H=2 H=3 H=4 H=2.2 H=3.2 H=2.2.2   # nullity 1: (3) (2,2)
```

**S5. Corollary R needs M_irr only below the genus.**
- A stable graph of genus g without legs and with at least one edge has no vertex of genus g. So Theorem R (ii)-(iv) uses M_irr(h) only for h <= g-1. The polar part of the single-vertex graph is then fixed by (H) rather than by M_irr(g). This is the main-session observation in round2\brief.md s.1; I checked it.
- Consequences: the all-genus regularity theorem holds through genus 8 without KP (from M_irr(h) for h <= 7), and through genus 21 with KP.
- leaf_structure_proof.md s.13.2 still assumes M_irr(h <= g); update it after a referee check.

**S6. The odd-order conifold relations add nothing to KP.**
- Setup: a patched copy of round2\R1_rigidity\global_suite.pl and kpmod.pm. This code is unrefereed; the copy is in scratchpad\synth\r1copy, and the patch adds a variant with odd orders only (cf=3) and a mode with no regularity conditions (none).
- With levels >= 3 unknown and propagated by KP, the closed-form rows of odd order have rank 0 at G = 4, 5, 6. The even orders give rank 3, 11 and 22.
- So new information at odd levels can only come from KP regularity at the special framings, from the one-leaf formula, or from geometry beyond the conifold polar part.

**S7. Two-leaf rigidity with the Hurwitz window.** Single-level problem at D = 3, with levels 1 and 2 imposed and window sigma..sigma^(G-1):
- Two-leaf regularity gives nullity 0 at G = 5 and 6, and also at G = 9 and 10 in R1's outputs.
- At G = 4 it gives nullity 2; full regularity gives 0.
- Compare Route A's TLR (no Hurwitz window): nullity 11, 11 and 7 at G = 4, 5 and 6.

## 2. Routes: state, obstruction, best continuation

**A (Marino-Vafa + KP).**
- State: 1.3-1.5 and 1.10.
- Obstructions:
  - TLR(G) for G >= 21 has no structural proof. At a fixed framing the relations only propagate upward; the rigidity is global over about G^2 special framings.
  - Cut-and-join alone creates no parts >= 3.
  - A level-by-level argument stops at the low odd levels (S3).
- Continuation:
  - P3/P4, now a small residual problem (S3).
  - Exact certificates.
  - The write-up fixes of Theorem D (s.5.3).

**B (Hurwitz, ELSV, Toda).**
- State: 1.1-1.2.
- Obstruction: the CY analogue of Okounkov's equation at P' = q^rho is d_{p_1}(e^{-H} Eq e^{H}) = exp(H(e^beta o p) - 2H + H(e^{-beta} o p)).
  - Eq = Sum_n (-i lambda)^n/n! F_n is diagonal: d/dp'_1 has eigenvalue (q^{1/2} - q^{-1/2}) Sum_box q^{-c(box)}.
  - But e^{-H} F_n e^{H} is non-linear for n >= 1, which breaks the counting argument.
- Continuation: P6, starting with r = 1.

**C (intersection theory).**
- State: Theorem H, Props. L and L1, C3(a).
- Obstruction: for 1 <= k <= h-1, the Hodge classes of the zero part Z of the section are not controlled by a dimension count. A per-support excess (virtual) class is needed. There is also a risk of a vacuous recursion through rational tails.
- Continuation: P5.

**D (localization, conifold).**
- State: 1.6 and 1.9; (H) is a theorem.
- Obstruction: the conifold gives one relation per polar order, against p(k) unknown functions at level k. After S1 and S6 its polar part is used up: it fixes the even levels and is implied by KP at the odd ones.
- Continuation: write-up only, plus P8.

**E (formulas, Perl engine).**
- State: 1.8 and 1.11. SN, NR and the k = h-1 closed form are conjectural, verified through genus 7.
- Obstruction: at fixed genus, SN is equivalent to Lemma A; its only extra content is that the coefficient series C_s are independent of the genus, which is unexplained. (6.1) is unproved for h >= 8.
- Continuation: P2. Together with Corollary E, (6.1) alone settles Lemma A for all h.

## 3. Ranked plan for round 2 (TR-free)

**P1. Reduction theorem, written and refereed.**
- Difficulty: low. Value: essential.
- Write S1-S5 as one proposition with exact hypotheses.
- Give exact-over-Q certificates for the genus-4 rank, and repeat it with an independent implementation (the A referee's ref_kp_null.pl with Hurwitz rows).
- Update Corollary R and Theorem F's base.

**P2. One-leaf vanishing: Xi_h(j) = 0 for 2 <= j <= 2h-2, all h (the vanishing part of (6.1)).**
- Difficulty: medium-hard. Payoff: decisive.
- With Corollary E this gives Lemma A for all h, using KP only. With S1-S3, only the even j = 4, 6, 8, ... at the open levels are needed.
- Methods:
  1. Schur side: in sum_a y^a/a! M_a the unit legs resum to a Gamma_-(y) insertion (A Cor. 2.3). With the Okounkov-Pandharipande E-operator formalism this should give a closed contour formula, like the one-part Marino-Vafa formula of A s.1. Then take the u = oo limit genus by genus (A Cor. 2.2).
  2. Route E's residue form: Zhat_j(u) = Res_{l=0} l^-3 (1+u/l)^{j-1} S1 S2 B^j dl must be a polynomial of degree <= j-1 (E s.1.3).
  3. P6.
- Checks:
  - (6.1) is verified for h <= 7. Round 2's R4 is running h = 8.
  - Every derivation must reproduce the two proved components, k = 0 (B Theorem B) and k = h ((G)).
- Main risk: the u = oo continuation lies outside the radius of convergence.

**P3. Level-3 rigidity: Xi_G(4) = 0 for every G >= 4, given lower genera (levels 1 and 2 are automatic by S2).**
- Difficulty: medium-hard.
- Unknown: P_(4) = sigma(1+4sigma)(2+9sigma) q(sigma) with deg q <= G-4.
- Method:
  1. Write the KP relations at levels 5, 6, 7 with level-3 sources (for example E_0(1,1; empty) = sigma/2).
  2. Impose regularity at tau = -k/a, a <= 2G-2. This gives an explicit, heavily overdetermined linear system on q (and on P_(6), P_(7)).
  3. Prove it is injective uniformly in G, using a framing triangularity or an interpolation in the special points.
- Evidence: R1's single-level runs, with nullity 0 for 4 <= G <= 10.
- For 4 <= G <= 7, P3 is all of Lemma A(G).

**P4. All odd levels with Z(D) <= G-2 (equivalently TLR-type rigidity), all G.**
- Difficulty: hard.
- The same method as P3, with the coupling D -> D+2, D+3. Alternatively, a structural proof of TLR from (V(z)-V(w))/(Y(w)-Y(z)), using Lagrange inversion for Y with Pi_k = C(k(tau+1)-1, k-1).
- The residual counts are small (G-3, G-7, G-14, ...).

**P5. lambda_1 component for every R (Route C's (T1)).**
- Difficulty: medium.
- Method: a plumbing model at a separating node, Fulton's excess formula, and Props. L and L1.
- Payoff:
  - It removes the sigma^(G-1) coefficient of every P_R.
  - Genus 4 then closes on its own.
  - Each open level loses one unknown.
- Data: Route E's component tables.

**P6. CY Toda for the framed vertex (Route B continuation).**
- Difficulty: hard.
- Do the T -> oo leading-order analysis of the deformed equation (s.2B), component by component in 1/sigma, starting with r = 1. This targets (6.1) at k = 1, 2, .... The boundary data k = 0 and k = h are already proved.

**P7. Exactness and bookkeeping.**
- Difficulty: low.
  - Exact certificates for Theorem F at G = 8..20.
  - DVV-free exact Marino-Vafa confirmation of M_irr(5) (and samples at genus 6 and 7).
  - Finish the author's genus-7 run (3 points left).
  - The write-up fixes of s.5.3 and s.5.4.

**P8. Further TR-free geometric relations.**
- Status: speculative.
- Candidates: other tori in the class of Brini-Schuler Prop. 3.1 (alpha0 != 0), or the shifted local P1 (Brini-Schuler s.3.2).
- Such a geometry helps only if it has a frame whose singular vertex carries the Lemma A class, and if its odd-order relations are not implied by KP. By S6, the conifold's odd-order relations are.
- Local P2 is excluded unless its all-genus conifold polar part has a TR-free proof, which the project does not have.

**Mapping to the existing round-2 folders.**
- R1_rigidity corresponds to P3 and P4.
- R2_onepoint corresponds to P2. Audit where its x-y swap formula comes from (s.5.1), and leave out its trcheck\ material.
- R3_toda corresponds to P6.
- R4_geometry corresponds to P5, and also runs (6.1) at genus 8.
- TR1_dictionary and TR2_recursion are excluded under "no TR at all".

## 4. Consolidated Mathematica checks (for the main session)

**M1 (medium).** Run E_formulas\E_crosscheck_g5.wls with xval = 1/3 (or symbolic x).
- Points: 8 zeros, (2), (3), (8), (2,2), (2,7), (3,6), (2,2,2), (2,2,6); and 3 controls, (9), (2,8), (2,2,7).
- The x = 2 run confirmed all 11 (E_crosscheck_g5_x2_out.txt).
- Expected: the same result, with control values equal to expected_g5.m. For example, the k = 0 part of (9) is -4401/112640, which is Route B's Hurwitz value.

**M2 (medium).** Copy E_crosscheck_g5.wls and add the points Route A and Route C asked for: (4), (3,3,3), (2,2,2,2,2,2,2), (2,3) and (3,4), at one rational x.
- Expected: all 0.
- Why it matters: this is the only Hodge-side genus-5 check that is independent of the Perl engines (which use DVV) and of KP.

**M3 (low).** Run E_formulas\E_crosscheck_g4r2.wls (it loads hvc_mp_header.m), then `perl compare_g4r2.pl`. This compares symbolic Xi[4,2] and Xi[3,3] with engine_h4_r2.txt and engine_h3_r3.txt, the data behind SN and NR. Expected: 0 mismatches.

**M4 (optional).** Route A's exact NullSpace of the Theorem D system at G = 5, symbolic in tau. No script exists yet, and exact Perl certificates already exist (A_mv\runs\exact_g5_full.out). An exact certificate of S4 is more useful.

**No longer needed:**
- The optional genus-5 Hurwitz-component runs of B and C. The referees' independent Hurwitz numbers, Route E's engine and the x = 2 Mathematica run already agree, e.g. Xi^(0)_5(9) = -4401/112640 and Xi^(0)_5(2,8) = -1024/4455.
- Route D's local-P2 leaf run (see P8 and S6).
- Route C's re-runs, unless lemmaA\data is regenerated.

**Perl checks (no wolframscript):**
1. An exact-over-Q certificate of S4.
2. `perl a3_cj_line.pl 2 6 > runs/cj_line_g2.out` (missing in A_mv\runs).
3. DVV-free exact Marino-Vafa evaluation of the 44 genus-5 points at 5 values of sigma.
4. KP Hirota equations to weight 8-10 at genus 2 and 3, with all terms.
5. R1-type single-level runs at D = 3, 5, 7 for G = 21..30 with window G-1, as evidence beyond Theorem F.

## 5. Completeness critique

### 5.1 Topological recursion

**(a) Round 1 is clean.** All five referees audited their routes and found no TR input of any kind (no TR, spectral curve, loop equation, BKMP, P2, Lemma D or mod-P2 result), and this synthesis uses none. Route A's Y(z) = z - sigma Sum Pi_k z^-k/(k+1) is the genus-0 dispersionless Lax inverse built from Marino-Vafa genus-0 data. In TR language it is the framed C^3 curve, but here it enters only as genus-0 algebra.

**(b) Conflict with round 2.** round2\brief.md s.0 records a later user clarification. It allows TR for the Hodge CohFT, which covers the curve C_tau, BKMP for C^3, DOSS, P2, Theorem 3.1 of proof_M_hurwitz.md, Lemma D, the mod-P2 results and the TR scripts. The request governing this synthesis says "no TR at all". The round-2 folders that use TR are:
- TR1_dictionary. TRtau.pm is an exact Eynard-Orantin recursion on the framed-vertex curve C_tau.
- TR2_recursion.
- R2_onepoint\trcheck (status_trC1.pl, status_trhalf.pl, status_cmp61.pl).
- R2_onepoint\cmp_tr_xy.pl, which compares against EO output.

Also check where R2_onepoint's x-y swap closed formula comes from. It is TR-derived if it rests on TR x-y duality, and TR-free if it comes from the Marino-Vafa formula or Lagrange inversion. Until the user confirms which rule applies, none of this should enter the proof record or the status files, not even as cross-checks.

**(c) DVV.** The Perl engines of D and E, and of E's referee, compute Witten-Kontsevich numbers with the DVV (Virasoro) recursion for the point, which is formally equivalent to TR on the Airy curve.
- It is used only to evaluate classical intersection numbers. These are covered by the Witten-Kontsevich theorem, with TR-free proofs by Kontsevich, Okounkov-Pandharipande and Kazarian-Lando.
- DVV-free engines cross-validate the tables: the Marino-Vafa Schur side (A and A's referee) and Hurwitz numbers (B and C).
- M_irr(5..7) rests on the DVV engines, plus partial DVV-free mod-p Marino-Vafa checks. P7 would remove this.

**(d) Inputs that are near TR but are not TR.**
- KP and 2D Toda, and the dispersionless Hirota equations.
- Givental-Teleman reconstruction and the genus-0 equivariant mirror theorem behind the frame-B premise. Theorem LC's Beta integral is a genus-0 oscillatory integral, not a mirror curve, and DOSS is never invoked.
- Getzler-Pandharipande's psi^{g-1} lambda_g lambda_{g-1}, which they derived from Virasoro. It appears only as a check value in E.
- The lambda_g formula. Cite Faber-Pandharipande, Ann. of Math. 157 (2003) (localization), not Getzler-Pandharipande (C's H0; also check A s.2.4 and proof_status s.1.3).

### 5.2 Dependence on literature and normalisations

- **ELSV**: the normalisation is checked against the data (236 + 19 comparisons). Drop the 81 comparisons made on the superseded tables.
- **Okounkov eq. (10)**: checked against arXiv:math/0004128v1, including s.2.7. The MRL 7 (2000) journal reference is unverified.
- **Farkas-Pandharipande Thm 3 / s.0.1**: checked; cite them precisely.
  - Boissy's dimension count is Lemma 3.6 / Remark 3.8, in s.3.3.4 of arXiv:1211.4951.
  - P4 (proof_status s.1.5, whose citation is incomplete) is now subsumed by Theorem H and by B's Theorem A, so D5.2b no longer depends on it.
- **Marino-Vafa formula (Liu-Liu-Zhou Thm 1; Okounkov-Pandharipande)**: the normalisation is fixed exactly.
- **KP**: cite the GL(oo)-orbit theorem precisely.
  - The genus expansion (differential Fay identity) was re-derived by A's referee.
  - The full Hirota equations were checked only at genus 1 (to weight 8) and genus 2 (to weight 6).
- **Brini-Schuler Prop. 3.1**: the statement and hypotheses are checked here. Its own inputs (Li's degeneration formula, BP08 Lem. 6.3, BP05 Thm 5.1) are cited, not re-proved.
- **Frame-B premise** (Theorem R input (1); leaf_structure_proof.md gap G3).
  - The premise: the frame-B graph sum (rc_lib.m/hvc_lib.m, with R' = R_Z h_w^-1 and R_Z from the stationary-phase expansion of the equivariant superpotential) is the GW potential of Z.
  - It should follow from Givental-Teleman reconstruction, quantum Riemann-Roch and Mumford/Faber-Pandharipande. However, the project has no written proof with citations; it is only validated by the exact genus 2-4 graph sums.
  - What depends on it: 1.6, the second proof in 1.3, S1, the genus-3 case of S4, S5, and Corollary R.
  - What does not: the statement of Lemma A, and Routes A, B, C and E.
- **Theorem LC's analytic input**: Stirling's series with remainder (DLMF 5.11), cited.
- **Graber-Vakil (D2.4)**: no exact statement or normalisation is given, but nothing depends on it.

### 5.3 Logical gaps

1. **The main gap: Lemma A for all h.** It reduces to any one of the following:
   - the one-leaf vanishing at the odd levels with Z(D) <= G-2 (S3), using KP and the conifold;
   - TLR(G) for G >= 21, using KP;
   - P2, using KP.
2. **KP dependence.** Everything after A s.3.4, and all of S1-S4, uses KP. The KP-free results are 1.1-1.9 and 1.11.
3. **Theorem D write-up.** It must state:
   - the projection onto the top exponential sector at u = oo;
   - the corrected sentence on pure genus-0 terms;
   - Lemma 3.1 with the r = 1, 2 genus-0 values and the genus-1 values, e.g. Q_1(s;(j)) = (j-1-1/sigma+s)/24.
4. **Route D's (V)** is the colour-1 term only. Add that the colour-2 single-vertex term is analytic (F1, F3).
5. **New and unrefereed material.**
   - S2, S4 and S5 have not been refereed.
   - S4's genus-4 rank uses one prime, one seed and the author's code.
   - S6 and S7 use unrefereed round-2 code.
6. **Theorem F.** It is exact over Q only at G = 5 and 6 (full system) and G = 7 (two-leaf system); elsewhere it is mod p with lifting. It is rigorous up to the correctness of three implementations: two by the author, sharing one design, and one by the referee.
7. **Route E.**
   - SN and NR are equivalent to Lemma A at fixed genus.
   - "A_kk proved given NR" is unsupported; it is computer-verified for k <= 7.
   - SN needs its ORD=2 normalisation stated.
8. **Route C.** (T1) and the product formula for k >= 2 are conjectural.
9. **Route D.**
   - D3.1 is a degree count, not a no-go theorem.
   - "The root s = 1 is simple" in D2.3 is a computer observation.

### 5.4 Bookkeeping

- Save the Route E referee result as E_formulas\referee.md.
- **A**: change 2070 to 2170; add the missing runs\cj_line_g2.out; cite the referee's reproduction of Theorem F.
- **B**:
  - The Hurwitz vanishing theorem is B's Theorem A; round2\brief.md s.2(b) calls it Theorem B.
  - Split the mKP claim (Hurwitz series computer-verified; Marino-Vafa case proved modulo literature).
  - Document or delete the r = 1 CY ODE remark.
  - State the Conjecture-M consequence for h >= 2.
  - Say that the genus-0 inputs for r = 1, 2 are conventions.
- **C**:
  - Change 286 to 266.
  - H' is the box form of Theorem H, not a strengthening.
  - Correct the Boissy location.
  - Prop. L needs r_+ >= 2k.
  - Relabel L1 as proved modulo literature.
  - Replace Lemma V by the projection formula.
- **D**:
  - D5-D5.3 are proved modulo literature, no longer conditional.
  - Delete obstruction 4.
  - loc_twisted.pl has no output file.
- **E**:
  - M_irr(7) is complete.
  - Check (c) of consistency.pl is not implemented.
  - The explore_lines claim holds for every complex l (gcd over Q).
  - comp_hm1.pl needs HODGE_NCH >= h-1.
- **round2\brief.md s.2(g)** still calls genus 7 partial.

## 6. Main-session results of the requested checks (2026-10-31)

All these runs avoid the special framings tau = -k/a. Every rational x > 0 gives tau = -x/(1+x) in (-1, 0), which is special. The earlier x = 2 run was at tau = -2/3, where Theorem C' already forces one-leaf zeros with 3 not dividing j. The runs below use x < 0 instead.

- **M1: done, passed.** E_formulas\E_crosscheck_g5_M1.wls at x = -2/9 (tau = 2/7, sigma = 18/49). All 11 engine predictions are confirmed exactly: 8 zeros, and the 3 controls (9), (2,8), (2,2,7) at their non-zero values. Faber and Getzler-Pandharipande sanity values are correct. Run time 31 min. Output: E_crosscheck_g5_M1_out.txt.
- **M2: done, passed.** E_formulas\E_crosscheck_g5_M2.wls at x = -3/5 (tau = 3/2, sigma = 15/4). The predictions come from expected_g5_M2.pl, the engine run on the new points. (4), (2,3), (3,4) and (3,3,3) vanish, and the control (4,6), at D = 2h-2, equals its non-zero engine value. Run time 31 min. Output: E_crosscheck_g5_M2_out.txt.
  - The point (2,2,2,2,2,2,2) was left out. This implementation runs genus-5 integrals with up to 26 marked points, which is out of reach at r = 7.
  - Like the Perl engines, this Mathematica code computes Witten-Kontsevich numbers by the DVV recursion. So it is an independent implementation, but not a DVV-free one.
- **M3: done, passed.** E_formulas\E_crosscheck_g4r2.wls (hvc_mp_header.m + hvc_lib.m, symbolic in x and in the parts jj[i]), then compare_g4r2.pl. Symbolic Xi[4,2] matches engine_h4_r2.txt and Xi[3,3] matches engine_h3_r3.txt: 280 and 589 component coefficients, 0 mismatches, 1/sigma-structure ok in both. These are the data behind Route E's SN and NR fits. Run time 18 s. Output: compare_g4r2_out.txt.
