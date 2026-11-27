# Lemma A round 2: route R4 (saved by the main session from the workflow journal, 2026-10-31)

Subagents cannot write report files, so this file holds the agent's structured return value: the fields first, then the full report text (report_text) verbatim.

## all_genus_claim

0

## all_genus_statement

No proved claim implies (OL), (6.1) or Lemma A for all genera. Proved for all genera here: the lambda_1-component (k = 1) of the one-leaf case (Theorem K1) and the boundary value K1'. Together with Theorem H (k = 0) and the lambda_g formula (k = h), this leaves the components 2 <= k <= h-1.

## claims

### [1]

#### depends_on

Route B Theorem C (REC) and Theorem B (ELSV; Okounkov math/0004128 eq. (10); Hurwitz genus-0 formula; referee: sound); Mumford's ch_1 formula; boundary restriction of the log tangent bundle; K-twisted string/dilaton. Not KP, not (6.1), no TR.

#### evidence

Complete proof in report s.2 (Props M1, M2, M3, Claim E). Every step checked exactly: M1 at genus 2-5 (p1moments.pl), M2 against the two-leaf tables for h <= 7 (lam1_rec.pl), Claim E's identities (V1)-(V3) (claimE_proof_checks.pl) and directly for h <= 12 (claimE.pl); the resulting Xi^(1) equals the Hodge-engine tables for h <= 9.

#### statement

Theorem K1: for every h >= 2 and every integer 2 <= j <= 2h-2, Xi^(1)_h(j) = Int_{Mbar_{h,1}} c(T(-log)) lambda_1/(1 - j psi) = 0 (the lambda_1-component of the one-leaf case of Lemma A, all genera).

#### status

proved-modulo-literature

### [2]

#### depends_on

Same as Theorem K1 (it is H(h) at j = 2h-1).

#### evidence

Hand check at h = 2 (1/30 on both sides); delta_bdry.pl for g <= 9.

#### statement

Boundary value: 12 Xi^(1)_h(2h-1) = (2h-1 - 1/(2h-1)) X_h(2h-1) + X_{h-1}(2h-1), X_g = Hurwitz one-leaf function of Theorem B.

#### status

proved-modulo-literature

### [3]

#### evidence

Mumford's 12 lambda_1 = kappa_1 - psi + delta; xi^* T(-log) = T_Gamma(-log) + O; dilaton; kappa_1 = pi_* psi^2. Exact check at genus 2-5 for all j, with p_1 computed independently (p1moments.pl).

#### statement

Prop. M1 (Mumford reduction): 12 Xi^(1)_h(j) = p_1 - (X_h(j) - X_h(0))/j + (1/2)(j-2h+2)(j-2h+3) X_{h-1}(j) + Sum_{g1=1}^{h-1} (j-2g1+1) X_{h-g1}(0) X_{g1}(j), with p_1 = Int c(T(-log)) kappa_1 lam_j = (2h-1) X_h - B_h, B_h = d/dy Xi^(0)_h(j,y)|_{y=0}.

#### status

proved

### [4]

#### depends_on

Route B Theorem C and Theorem B.

#### evidence

rec2_test.pl (REC at random rational points, genus <= 5); B_h from (A) equals the two-leaf engine tables for h <= 7 (lam1_rec.pl).

#### statement

Prop. M2: Route B's two-leaf REC holds as a polynomial identity in (j,y); differentiated at y = 0 it gives the recursion (A) for B_h in terms of B_{h''} (h'' < h), Theorem B data T_n(j) = (j-1)_(2n)[t^{2n}]H^j and zeta_m = (2m-1)B_{2m}/(2m).

#### status

proved-modulo-literature

### [5]

#### depends_on

Props M2 and Claim E.

#### evidence

Induction in report s.2; hbdry_rec.pl confirms it from REC data for g <= 12.

#### statement

Prop. M3 (boundary-sum identity): d/dy Xi^(0)_g(j,y)|_{y=0} = Sum_{g1=0}^{g-1} chi(M_{g-g1,1}) T_{g1}(j)/j for every integer 1 <= j <= 2g-1 (separating-boundary localization of the kappa_1 insertion of the Hurwitz class).

#### status

proved-modulo-literature

### [6]

#### evidence

Generating-function proof in report s.2: reduces to X(w/(1+w)) - X(w) = -log(1+w) + w/2 + w/(2(1+w)) (expansion of psi(x+1) - psi(x) = 1/x) and t(log H)' = 1 - S^2. Checks: claimE_proof_checks.pl ((V1) to w^30, (V2) 840 coefficients, (V3) N <= 16) and claimE.pl (144 pairs, h <= 12).

#### statement

Claim E (Bernoulli/Lambert identity closing the induction) holds for all h >= 1 and integers 1 <= j <= 2h-1.

#### status

proved

### [7]

#### depends_on

Theorem K1, Theorem H (Route C) or Route B Theorem A, lambda_g formula (Faber-Pandharipande 2003).

#### statement

Corollary: in the Lemma-A range for one leaf, sigma^h Xi_h(j) lies in span(sigma, ..., sigma^{h-2}); only the components 2 <= k <= h-1 remain open.

#### status

proved-modulo-literature

### [8]

#### evidence

Equal as polynomials in j to the REC+Mumford values for h = 1..16 (lam1_rec_h16.txt) and to the Hodge-engine tables for h <= 9.

#### statement

Closed form of the lambda_1-component for all j: Xi^(1)_h(j) = -(1/4)(j-1)_(2h-2) [t^{2h}] S^2 H^j [(1+j/3)(S^2-1) - (j/3) t S S' + t S'/S] (the order-1/sigma part of (6.1)).

#### status

computer-verified

### [9]

#### evidence

hbdry_rec.pl, g <= 12 (REC data only); delta_bdry.pl, g <= 9.

#### statement

The boundary-sum identity M3 also holds at j = 2g (equivalently 12 Xi^(1)_g(2g) = 3 X_{g-1}(2g)).

#### status

computer-verified

### [10]

#### depends_on

Route A Theorem D (KP, allowed); Theorem H; Theorem K1.

#### evidence

tlr_k01.pl (Route A's a3_kp_nullity2.pl plus constraints), mod p = 2^31-1; tlr_k01_twoleaf.txt, tlr_k01_full.txt. Rank deficiency 0 is rigorous by Route A's lifting argument.

#### statement

With Theorem H and K1 added as linear constraints, Route A's KP system (Theorem D + polynomiality + lambda_g) has nullity 0 at G = 4 (full system; it was 6) and G = 5; the two-leaf subsystem has nullity 1, 1, 0, 0 at G = 3, 4, 5, 6 (it was 7, 11, 11, 7). At G = 3 one unknown remains: Xi^(2)_3(4). Hence Lemma A(4) follows from Lemma A(<=3) without exact genus-4 data, and TLR holds for 5 <= G <= 20.

#### status

computer-verified

### [11]

#### evidence

c61.pl: genus 7 (referee table) 230, genus 8 (new table) 286, genus 9 (new table) 348 comparisons, 0 mismatches. The new genus-8/9 one-leaf tables were computed with the referee engine (t_moments.pl; 1320 s and 4927 s).

#### statement

Formula (6.1) holds identically in (j, tau) at genus 7, 8 and 9 (Hodge side, no TR).

#### status

computer-verified

### [12]

#### evidence

excess_split.pl at genus 3 and 4 (excess_split_h3.txt, excess_split_h4.txt); the k = 1 parity identity is a direct degree count.

#### statement

For r = 1 and k >= 2, the excess (jump-locus-localized) part E^(k)_h(j) of Xi^(k)_h(j) does not vanish in the range (e.g. genus 3: E^(2)(4) = 1/360; genus 4: E^(2)(5) = 1/1728). So a localization-only proof at r = 1 cannot work for k >= 2. For k = 1, E^(1) = Sum_{a<=j-2} binom(j-1,a) Xi^(1)(phi_a) + (1+(-1)^j) Xi^(1)(phi_{j-1}) and NE^(1) = (-1)^{j-1} Xi^(1)(phi_{j-1}).

#### status

computer-verified

### [13]

#### evidence

irr_test.pl, genus 3-6 (irr_test_out.txt).

#### statement

Single Chern-character insertions ch_l (l = 3, 5) into the Hurwitz one-leaf function do not localize to the irreducible boundary in the Lemma-A range; only the CY combination e_k(alpha+alpha^2) vanishes.

#### status

computer-verified

### [14]

#### evidence

cyH_test.pl: fails for every k >= 1 and g <= 7 at all j >= 2.

#### statement

The CY analogue of the boundary-sum identity M3 (component k >= 1, with the Hodge class distributed over the separating boundary) holds.

#### status

refuted

### [15]

#### evidence

bc_fit.pl, g <= 7, c = 2, 3, 4 (bc_fit_out.txt); the same script validates the y-expanded REC against the tables.

#### statement

The higher Taylor coefficients [y^c] Xi^(0)_g(j,y), c >= 2, are boundary sums Sum_{g1} A^(c)_{g-g1} T_{g1}(j)/j with universal constants in the low range.

#### status

refuted

### [16]

#### evidence

k2_reduce.pl: exact agreement with the component tables at genus 2-6 for all j (and the two-leaf lambda_1 reduction at genus 2-5).

#### statement

Exact reduction of the k = 2 component: Xi^(2) = (1/24) <lambda_1(kappa_1 - psi + delta)> + (1/120)[<kappa_3> - <psi^3> + (1/2)<iota_*(psi'^2 - psi'psi'' + psi''^2)>]. Each term reduces to Hurwitz zero-leaf derivative data with at most three leaves; in the range, the k = 2 statement is equivalent to the single identity H_2: 5 D_h - 7 [y^2]Y_h - 12 [y^3]Y_h = explicit lower-genus data.

#### status

proved

### [17]

#### evidence

Elementary manipulations in report s.6.

#### statement

Reformulations: Lemma A(r=1) <=> W_h(z) - W_h(0) = O(z^{2h-1}), W_h = Sum_j sigma^h Xi_h(j) z^j; (6.1) <=> W_h = (2h-2)! [l^{2h}] S1 S2 (1 + (zB/(1-zB))^{2h-1}) <=> Sheffer form Sum_h Q_h(j) mu^{2h} = S1 S2 B^j <=> unit-leg dressing sigma^h <Omega K f> = [l^{2h}] S1 S2 b_h^{-1} Sum_m ((1-B^{-1})^m/m!) <lambda_h K f lam_1^m>.

#### status

proved

### [18]

#### evidence

Mumford's relation Lambda^vee(w) = w^{2h}/e(E(x)C_w); Serre duality. No vanishing statement is derived from it (the degeneracy locus contains the zero section).

#### statement

Cosection observation: the CY one-leaf integrand is the zero-section localization of Tot(E(x)C_tau + E(x)C_{-1-tau}) with obstruction T(-log)(x)C_1, and the product of abelian differentials gives a cosection of weight 1 + tau + (-1-tau) = 0, invariant exactly by the CY condition.

#### status

proved


## key_formulas

- 12 Xi^(1)_h(j) = p_1 - (X_h(j) - X_h(0))/j + (1/2)(j-2h+2)(j-2h+3) X_{h-1}(j) + Sum_{g1=1}^{h-1} (j-2g1+1) X_{h-g1}(0) X_{g1}(j),   p_1 = Int c(T(-log)) kappa_1 lam_j = (2h-1) X_h - B_h
- X_g(j) = (j-1)_(2g-2) [t^{2g}] S^2 H^j,   T_n(j) = (j-1)_(2n) [t^{2n}] H^j,   Y_g(j,0) = T_g(j)/j,   chi_m = -B_{2m}/(2m),   zeta_m = (2m-1)B_{2m}/(2m)
- (2h+1)[(j-2h+1)X_h + (j-2h)B_h] = Sum_{h''<h} c_{h+1-h''}[ d/dx (x)_(2h+2-2h'')|_{j-2h''} T_{h''}/j + (j-2h'')_(2h+2-2h'') B_{h''} ] + Sum_{m=1}^h zeta_m T_{h-m},   B_0 = -1/j^2
- Boundary-sum identity: B_g(j) = Sum_{g1=0}^{g-1} chi_{g-g1} T_{g1}(j)/j for integers 1 <= j <= 2g-1 (observed also at j = 2g)
- X(w/(1+w)) - X(w) = -log(1+w) + w/2 + w/(2(1+w)),   X(w) = -Sum_{m>=1} B_{2m} w^{2m}/(2m)
- (N-1) Gamma_{N-1} = Sum_{n<=N-2} (j-2n) Gamma_n/(2N-2n)!   <=>   t (log H)' = 1 - S^2,   Gamma_n = [t^{2n}] H^j
- Xi^(1)_h(j) = -(1/4)(j-1)_(2h-2) [t^{2h}] S^2 H^j [ (1+j/3)(S^2-1) - (j/3) t S S' + t S'/S ]   (verified h <= 16)
- 12 Xi^(1)_h(2h-1) = (2h-1 - 1/(2h-1)) X_h(2h-1) + X_{h-1}(2h-1)
- S(e1 l)^{1+j e2} S(e2 l)^{1+j e1} = S^{2+j} e^{-j t S'/S} [1 + eps^2((1+j/2) t^2 f'' - (j/6) t^3 f''') + O(eps^4)],  w = tau+1/2, eps = 1/(2w), l = t/w, f = log S
- Xi^(2) = (1/24) <lambda_1(kappa_1 - psi + delta)> + (1/120)[<kappa_3> - <psi^3> + (1/2)<iota_*(psi'^2 - psi'psi'' + psi''^2)>],  <kappa_3> = (2h-1)X_h - Sum_{c=1}^3 [y^c] Y_h
- W_h(z) = Sum_j sigma^h Xi_h(j) z^j = (2h-2)! [l^{2h}] S1 S2 (1 + (zB/(1-zB))^{2h-1})  (form of (6.1));  Lemma A(r=1) <=> W_h(z) - W_h(0) = O(z^{2h-1})

## next_steps

- Prove H_2 (the k = 2 component): derive the REC recursions for [y^2]Y_h, [y^3]Y_h (two leaves, y-expanded REC as in bc_fit.pl) and for D_h = [yz] Xi^(0)_h(j,y,z) (three-leaf REC), find their low-range formulas (boundary sums with psi-insertions at the node), then prove an analogue of Claim E.
- Extend Theorem K1 to all R (r >= 2) with REC for R plus a zero leaf; this would replace Route C's plumbing model (T1) and give the lambda_1-component of Lemma A for all h and R.
- For R1: attempt a structural proof of the two-leaf KP rigidity TLR(G) using the two extra vanishing coefficients per one-leaf polynomial (unknowns sigma^1..sigma^{G-2}); rigidity now starts at G = 5.
- Next-to-lambda_g component (k = h-1): reduce Int K lambda_{h-1} lam_j via c(E^vee) e_{h-1}(alpha(1+alpha)) = lambda_{h-1} - h lambda_h + lambda_h Sum n! ch_n; together with K1 it would make genus 3 data-free and leave only k = 2..h-2.
- Look for a uniform framework: Omega K = R_sigma.(Hurwitz K-twisted CohFT); the order-sigma^{-k} term involves Hurwitz polynomials with <= k zero leaves; a general low-range formula for their zero-leaf Taylor coefficients, proved by REC, would give all components.

## obstructions

- (OL) needs Xi^(k)_G(D+1) = 0 for every 2 <= k <= G-1 (G-2 unknown coefficients per (G,D) after K1); only k = 0 (Theorem H), k = 1 (K1) and k = G (lambda_g) are proved for all genera.
- K1's mechanism uses one Mumford relation whose kappa-part is a single zero-leaf derivative that equals a boundary sum in the low range and closes under REC. For k >= 2 the required Hurwitz derivative data ([y^2]Y, [y^3]Y, three-leaf [yz]) are not boundary sums (refuted, g <= 7), and their REC recursions are not closed on the combination H_2, so the induction needs separate low-range closed forms that are not known.
- Jump-locus localization at r = 1: for k >= 2 the localized part does not vanish in the range (genus 3, 4), so it must be combined with non-localized terms; with unit legs (Route C, Prop. L) the self-referential rational tails s = z^{1-j} d/dz remain, with uncomputed coefficients.
- Single Hodge insertions (ch_3, ch_5) do not localize; the vanishing is a genuine CY cancellation (consistent with Route C 3(a)).
- Finite fixed-point localization on CY3 targets produces only K-free vertices; the K-twist appears only after resummation at the conifold point, so no FP-type finite relation for (6.1) was found. The invariant cosection gives localization to {ab = 0}, not vanishing.
- Genus independence of the SN/(6.1) series is explained only at order 1/sigma (via REC + Mumford); higher orders need the higher Mumford relations and multi-leaf REC.

## requested_checks

- Referee check of Theorem K1: Prop. M1 (Mumford's lambda_1 relation, boundary restriction xi^* c(T(-log)) = c(T_Gamma(-log)), dilaton factors), Prop. M2 (Route B's REC as a polynomial identity and its derivative at y = 0, including the genus-0 conventions), the induction M3 (which B_{h''} enter, j <= 2h''-1), and the generating-function proof of Claim E, in particular the formal identity X(w/(1+w)) - X(w) = -log(1+w) + w/2 + w/(2(1+w)).
- Optional independent confirmation (Mathematica, main session): Xi^(1)_h(j) = Int c(T(-log)) lambda_1/(1 - j psi) at h = 10 for two values of j outside the range, compared with the closed form C1 (report s.3.1); script not written (no wolframscript allowed here).
- Re-run tlr_k01.pl (K01=2, args G 0 2 for the two-leaf subsystem; G alone for the full system) at G = 3..8 with a second prime/seed, and the exact-over-Q variant of Route A for G = 4, to certify the rank statements of report s.3.3.
- If wanted: finish the NR two-leaf test at genus 6-7 on a machine with more time (in R4/nr_test run: perl fitNR.pl 2 7 1 and perl fitNR.pl 2 7; both timed out locally after genus 5 and genus 4).

## route

R4 (geometry): Mumford's GRR relations and boundary restrictions on Mbar_{h,1} (intersection theory), Route B's Hurwitz REC (ELSV + Okounkov's Toda) as the only recursion input, jump-locus/excess analysis of Route C, Givental structure; no TR of any kind

## summary

Main new result: Theorem K1. The lambda_1-component of the one-leaf case vanishes in every genus: Xi^(1)_h(j) = Int_{Mbar_{h,1}} c(T(-log)) lambda_1/(1 - j psi) = 0 for all h >= 2 and 2 <= j <= 2h-2. Status: proved modulo literature. The inputs are ELSV, Okounkov's Toda equation and Hurwitz's genus-0 formula (through Route B's REC and Theorem B), plus Mumford's lambda_1 formula. It does not use KP, (6.1) or any TR.

The proof has four steps.
- Mumford's relation 12 lambda_1 = kappa_1 - psi + delta, with boundary restrictions of the log tangent bundle, rewrites Xi^(1) in Hurwitz terms: a kappa_1 insertion, a psi term, and boundary terms (Prop. M1).
- The kappa_1 insertion is a zero-leaf derivative of the two-leaf Hurwitz polynomial. Route B's two-leaf REC, read as a polynomial identity and differentiated at the second leaf, gives a recursion for it (Prop. M2).
- An induction on h shows that this derivative equals the separating boundary sum for 1 <= j <= 2g-1 (Prop. M3). The induction closes by one explicit identity (Claim E).
- Claim E is proved with generating functions. It reduces to the formal functional equation of the Bernoulli generating function (from psi(x+1) - psi(x) = 1/x) and to t(log H)' = 1 - S^2, which is the ODE behind Theorem B.

Every step is also checked exactly against the Hodge engines (genus 2-9). Consequences:
- In the range, sigma^h Xi_h(j) lies in span(sigma, ..., sigma^{h-2}), so only the components k = 2..h-1 remain open.
- The boundary value 12 Xi^(1)_h(2h-1) = (2h-1-1/(2h-1)) X_h(2h-1) + X_{h-1}(2h-1) holds (proved).
- The closed form of the lambda_1-component predicted by (6.1) is confirmed through genus 16 by this engine-free route.
- Added to Route A's KP system, Theorem H and K1 make the full system rigid already at G = 4 (it was nullity 6) and the two-leaf subsystem rigid at G = 5, 6. Lemma A for h <= 20 then needs exact data only for one genus-3 number, Xi^(2)_3(4).
- (6.1) holds identically at genus 7, 8 and 9 (new one-leaf tables for genus 8 and 9, computed with the referee's engine).

Where it stops: (OL) needs every component 2 <= k <= G-1. Several mechanisms provably do not extend.
- At r = 1 the jump-locus (excess) part does not vanish for k >= 2.
- Single ch_l insertions (l >= 3) do not localize to the boundary.
- The CY analogue of the boundary-sum identity fails for every k >= 1.
- The higher Taylor coefficients of the two-leaf Hurwitz polynomial are not boundary sums.

The k = 2 component is reduced exactly (verified genus 2-6) to Hurwitz derivative data, all determined by REC with at most three leaves. What remains is the single identity H_2 (stated in the report), with no inductive closure yet. No all-genus proof of (OL), (6.1) or Lemma A is claimed.

## report_text

# Round 2, R4 (geometry): report

Folder: `G:\My Drive\Projects\Old\refGW\hodge_vertex_conifold\lemmaA\round2\R4_geometry\` (called R4 below). All arithmetic is exact (Perl Math::BigRat). Route A's KP ranks are taken mod p = 2^31-1, and a rank deficiency of 0 is rigorous by Route A's lifting argument.

No wolframscript was run, and at most 3 heavy processes ran at a time.

**No topological recursion of any kind was used:**
- no spectral curve, Eynard-Orantin invariants or loop equations;
- no Virasoro constraints used as loop equations;
- no BKMP or dictionary P2, and nothing proved only mod P2;
- no TR scripts, not even as a cross-check.

**Inputs:**
- Mumford's formula for ch_l(E) (Mumford 1983; Faber-Pandharipande, Invent. Math. 139 (2000), s.2).
- The restriction of the log tangent bundle to boundary divisors.
- String and dilaton equations for pulled-back classes with the K-twist.
- Route B (`B_hurwitz/report.md`, referee: sound): Theorem C (the recursion REC), Theorem B (closed one-leaf Hurwitz form) and Theorem A. Their inputs are ELSV, Okounkov's Toda equation (math/0004128, eq. (10)) and Hurwitz's genus-0 formula. Route C's Theorem H is the geometric alternative for the Hurwitz-component vanishing.
- KP (Route A, Theorem D) only in s.3.3, as allowed by today's ruling.

## Conventions
- c = 1; weights w = (1, tau, -1-tau); sigma = tau(1+tau).
- lam_j = 1/(1 - j psi); K = exp(-Sum kappa_k/k); c(T(-log)) = Lambda^vee(1) K (Lemma 2).
- Xi_h(R) = Sum_{k=0}^h (-1/sigma)^k Xi^(k)_h(R), with Xi^(k)_h(R) = Int c(T(-log)) e_k(alpha + alpha^2) prod lam_{j_i}.
- S(t) = (t/2)/sinh(t/2); H = S exp(-t S'/S); log H = log S - t (log S)'.
- (x)_(n) is the falling factorial; c_k = 2/(2k)!.
- X_g(j) := Xi^(0)_g(j) = (j-1)_(2g-2) [t^{2g}] S^2 H^j (Route B, Theorem B).
- T_n(x) := (x-1)_(2n) [t^{2n}] H(t)^x.
- chi_m := X_m(0) = chi(M_{m,1}) = -B_{2m}/(2m); zeta_m := (2m-1) B_{2m}/(2m).
- Y_g(j,y) := Xi^(0)_g(j,y), the two-leaf Hurwitz polynomial.
- Genus-0 conventions of Route B (fixed by (U3), (U4)): X_0 = 1/(j(j+1)), Y_0 = 1/(j+y).
- B_g(j) := d/dy Y_g(j,y) at y = 0.
- Useful identities: Y_g(j,0) = (j-2g+1) X_g(j) = T_g(j)/j; [t^{2g}] S^2 H^j = ((j-2g)/j) [t^{2g}] H^j, from S^2 H^j = H^j - (t/j)(H^j)'.

## 1. Results at a glance

| # | Statement | Status |
|---|---|---|
| K1 | Xi^(1)_h(j) = 0 for all h >= 2, 2 <= j <= 2h-2 (lambda_1-component of the one-leaf case) | proved-modulo-literature (REC/Theorem B of Route B; Mumford) |
| K1' | 12 Xi^(1)_h(2h-1) = (2h-1 - 1/(2h-1)) X_h(2h-1) + X_{h-1}(2h-1) (boundary value) | proved-modulo-literature (same) |
| M1 | Mumford reduction of Xi^(1)_h(j) to Hurwitz data | proved; checked exactly at genus 2-5, all j |
| M2 | Recursion for B_h from Route B's two-leaf REC | proved-modulo-literature; checked: B_h equals the two-leaf engine tables, h <= 7 |
| M3 | Boundary-sum identity B_g(j) = Sum_{g1=0}^{g-1} chi_{g-g1} T_{g1}(j)/j for integers 1 <= j <= 2g-1 | proved-modulo-literature (part of K1) |
| E | The Bernoulli/Lambert identity closing the induction | proved; identities (V1)-(V3) checked; direct check for h <= 12 |
| C1 | Closed form Xi^(1)_h(j) = -(1/4)(j-1)_(2h-2) [t^{2h}] S^2 H^j [(1+j/3)(S^2-1) - (j/3) t S S' + t S'/S] for all j | computer-verified h <= 16 (REC route) and h <= 9 (Hodge engines); conjectural beyond |
| M3' | M3 also holds at j = 2g, i.e. 12 Xi^(1)_g(2g) = 3 X_{g-1}(2g) | computer-verified g <= 12 |
| KP | Route A's KP system plus Theorem H and K1: full system rigid for G = 4, 5; two-leaf subsystem rigid for G = 5, 6; at G = 3 one unknown, Xi^(2)_3(4) | computer-verified (mod p) |
| 6.1 | (6.1) identically in (j, tau) at genus 7, 8, 9 | computer-verified (Hodge side, referee engine) |
| N1 | r = 1, k >= 2: the jump-locus (excess) part does not vanish in the range | computer-verified at genus 3, 4 (refutes a localization-only proof at r = 1) |
| N2 | Single ch_l insertions (l = 3, 5) do not localize to the irreducible boundary in the range | computer-verified at genus 3-6 |
| N3 | The CY analogue of M3, component by component (k >= 1), fails | computer-verified g <= 7 |
| N4 | The naive analogue of M3 for the higher Taylor coefficients [y^c] Y_g (c >= 2) fails | computer-verified g <= 7 |
| K2 | Exact reduction of Xi^(2)_h(j) to Hurwitz derivative data | proved (Mumford + boundary); checked at genus 2-6, all j |
| R | Reformulations (W_h(z), Sheffer form, unit-leg dressing) and the cosection observation | proved (equivalences / equivariance); no vanishing derived from them |

## 2. Theorem K1 and its proof

**Theorem K1.** For every h >= 2 and every integer 2 <= j <= 2h-2:
Xi^(1)_h(j) = Int_{Mbar_{h,1}} c(T_{Mbar_{h,1}}(-log)) lambda_1 / (1 - j psi) = 0.

### Prop. M1 (Mumford reduction)
For h >= 2 and all j:

    12 Xi^(1)_h(j) = p_1(h,j) - (X_h(j) - X_h(0))/j + (1/2)(j-2h+2)(j-2h+3) X_{h-1}(j)
                     + Sum_{g1=1}^{h-1} (j-2g1+1) X_{h-g1}(0) X_{g1}(j),

where p_1(h,j) := Int c(T(-log)) kappa_1 lam_j = (2h-1) X_h(j) - B_h(j). For h = 1 the irreducible term is 1/2.

*Proof.*
1. On Mbar_{h,1}, Mumford gives ch_1 = lambda_1 = (kappa_1 - psi + delta)/12. Here kappa_1 is the Arbarello-Cornalba class and delta = (1/2) iota_* 1 is the total boundary class.
2. psi/(1 - j psi) = (lam_j - 1)/j.
3. For a boundary gluing map xi, the residue sequence gives xi^* T(-log) = T_{Mbar_Gamma}(-log) (+) O. Hence xi^* c(T(-log)) is the product of the vertex classes c(T_v(-log)), and the node branches carry the trivial leaf lam_0 = 1.
4. These leaves are removed by the K-twisted dilaton equation Xi(0,R) = (D_R - 2g + 2) Xi(R), which holds componentwise. On Delta_irr this gives (1/2)(j-2h+2)(j-2h+3) X_{h-1}(j). On Delta_{g1} it gives (j-2g1+1) X_{g1}(j) chi_{h-g1}.
5. For the kappa_1 term, kappa_1 = pi_* psi_2^2 and c(T_{h,2}(-log)) = pi^* c(T_{h,1}(-log)) (1 - psi_2) give Int c(T) kappa_1 lam_j = Sum_{c>=2} [y^c] Y_h(j,y).
6. That sum equals Y_h(j,1) - Y_h(j,0) - B_h = j X_h - (j-2h+1) X_h - B_h, by string and dilaton.

Check: `p1moments.pl` computes p_1 independently from explicit dilaton points with the referee engine RefHodge.pm. Result: genus 2-5, all j = 1..3h+1, 0 mismatches.

### Prop. M2 (recursion for B)
Route B's REC for R = (j, y) reads

    (2h+1)(j+y-2h) Y_h = Sum_{k=2}^{h+1} c_k (j+y-2(h+1-k))_(2k) Y_{h+1-k}
                         + Sum_{h1+k1+h2+k2 = h+2} c_{k1} c_{k2} (j-2h1+1)_(2k1) (y-2h2+1)_(2k2) X_{h1}(j) X_{h2}(y).

1. It holds at all lattice points j, y >= 2 with j + y != 2h.
2. Both sides are polynomials in (j, y), because the genus-0 terms cancel against the falling factorials: (j+y)_(2k) Y_0 = (j+y-1)_(2k-1) and (x+1)_(2k) X_0(x) = (x-1)_(2k-2).
3. By Zariski density it is an identity of polynomials. Checked at random rational points, genus <= 5, in `rec2_test.pl`.
4. By Theorem B and Sum_{k>=1} c_k t^{2k-2} = S^{-2}, the two-block part equals Sum_{m1+m2=h} T_{m1}(j) T_{m2}(y).
5. Use T_m(0) = delta_{m0} and T'_m(0) = (2m)! [t^{2m}] log H = zeta_m (m >= 1), and differentiate at y = 0. For h >= 1 this gives (A):

    (2h+1)[(j-2h+1) X_h + (j-2h) B_h]
        = Sum_{h''=0}^{h-1} c_{h+1-h''} [ d/dx (x)_(2h+2-2h'')|_{x=j-2h''} T_{h''}(j)/j + (j-2h'')_(2h+2-2h'') B_{h''}(j) ]
          + Sum_{m=1}^h zeta_m T_{h-m}(j),          with B_0 = -1/j^2.

Check: `lam1_rec.pl` computes B_h from (A). It equals the two-leaf engine tables (E_formulas, and the referee's tables at h = 6, 7) for h <= 7.

### Prop. M3 (boundary sum) by induction
**Claim H(g).** B_g(j) = BS_g(j) := Sum_{g1=0}^{g-1} chi_{g-g1} T_{g1}(j)/j for all integers 1 <= j <= 2g-1.

Fix h and an integer j in [1, 2h-1], and use (A).
1. The term (2h+1)(j-2h+1) X_h vanishes. Indeed X_h(j) = 0 for j <= 2h-2 by Theorem B, and the factor (j-2h+1) is zero at j = 2h-1.
2. The coefficient (j-2h'')_(2h+2-2h'') = prod_{m=2h''}^{2h+1} (j - m) vanishes unless j <= 2h''-1. In that case H(h'') applies, with h'' < h.
3. T_{h''}(j) = 0 for 1 <= j <= 2h'', so only h'' < j/2 survive in the derivative terms.
4. Hence (2h+1)(j-2h) B_h(j) is explicit. Since j - 2h != 0, H(h) at j is equivalent to E(h,j) = 0 (Claim E below).
5. The base case h = 1, j = 1 is E(1,1) = 0, which gives B_1(1) = -1/12 = chi_1.

### Claim E
For all h >= 1 and integers 1 <= j <= 2h-1, E(h,j) = 0, where

    E(h,j) = Sum_{0 <= n < j/2} T_n(j) [ (2h+1)(j-2h) chi_{h-n} - c_{h+1-n} d/dx (x)_(2h+2-2n)|_{x=j-2n}
             - Sum_{h''=ceil((j+1)/2)}^{h-1} c_{h+1-h''} (j-2h'')_(2h+2-2h'') chi_{h''-n} - j zeta_{h-n} ].

*Proof.*
1. Put a = j - 2n and M = h - n. With zeta_M = -(2M-1) chi_M, the bracket equals

       C(M; a, n) = chi_M [2M(2a-2M-1) + 2n(a-1)] - 2 d/da binom(a, 2M+2) - 2 Sum_m binom(a-2m, 2M+2-2m) chi_m.

2. For integer a >= 1 and 2M >= a+1 (which is the range h >= (j+1)/2), the m-sum may start at m = 1. The added terms vanish because 0 <= a-2m < 2M+2-2m.
3. The n-sum may run over all n >= 0, because T_n(j) = 0 for 1 <= j <= 2n.
4. Hence E(h,j) = [w^{2h+2}] E_j(w), where E_j(w) := Sum_n T_n(j) w^{2n} Sum_{M>=1} C(M; j-2n, n) w^{2M+2}.
5. With the formal series X(w) := Sum_{m>=1} chi_m w^{2m}:

       Sum_{M>=1} C(M;a,n) w^{2M+2} = X(w) [(a-1)(a+2n) w^2 + 2 - (1+w)^a - (1-w)^a] + (2a-1) w^2
                                       - (w/2) [(1+w)^a - (1-w)^a + (1+w)^{a-1} - (1-w)^{a-1}].

   Two facts are used. First, (2a-1)D - D^2 + (a-D)(a-1-D) = a(a-1), with D = w d/dw. Second, the formal identity (V1): X(w/(1+w)) - X(w) = -log(1+w) + w/2 + w/(2(1+w)). This is the expansion at x = 1/w of psi(x+1) - psi(x) = 1/x, using psi(x) ~ log x - 1/(2x) + X(1/x).
6. So E_j = X P_j + R_j with deg R_j <= j+1, and (since a + 2n = j)

       P_j(w) = Sum_n T_n(j) w^{2n} [2n(j-2n-1) w^2 - 2 Sum_{k>=2} binom(j-2n, 2k) w^{2k}].

7. P_j = 0. Its w^{2N} coefficient is 2 (j-1)_(2N-1) [ (N-1) Gamma_{N-1} - Sum_{n<=N-2} (j-2n) Gamma_n/(2N-2n)! ], where Gamma_n = [t^{2n}] H^j.
   - The generating function of the bracket is (t^3/2) Gamma' - (j Gamma - t Gamma')(cosh t - 1 - t^2/2) = (cosh t - 1)(t Gamma' - j(1 - S^2) Gamma).
   - This vanishes because t Gamma'/Gamma = j t (log H)' = j(1 - S^2), and 2(cosh t - 1) = t^2/S^2. This is Route B's identity 1 - S^2 = -t^2 (log S)'', i.e. the ODE behind Theorem B.
8. Therefore E(h,j) = [w^{2h+2}] R_j = 0 whenever 2h+2 >= j+2. ∎

Checks:
- `claimE_proof_checks.pl`: (V1) to order w^30; (V2) the generating function, 840 coefficients; (V3) the bracket identity as polynomials in j, N <= 16. All pass.
- `claimE.pl`: E(h,j) = 0 directly, 144 pairs, h <= 12.

### Conclusion of the proof of K1
Prop. M1 gives B_h - BS_h = (2h-1 - 1/j) X_h - 12 Xi^(1)_h + (1/2)(j-2h+2)(j-2h+3) X_{h-1}.
1. For 2 <= j <= 2h-2, X_h(j) = 0 (Theorem B or Theorem H).
2. X_{h-1}(j) = 0 for j <= 2h-4, and the prefactor kills j = 2h-3 and j = 2h-2. So the third term vanishes.
3. Therefore H(h) gives Xi^(1)_h(j) = 0. ∎
4. At j = 2h-1, H(h) gives K1'.

The proof is an induction on h. Each step uses REC at genus h, Theorem B at all genera <= h, Mumford at genus h, and Claim E, which is uniform in (h, j). It needs no hypothesis on the other components or on lower genera.

## 3. Consequences, checks beyond the known range

### 3.1 Closed form C1 (conjectural for j >= 2h)
- Expand (6.1) at large sigma, with w = tau + 1/2, eps = 1/(2w), lambda = t/w, sigma = w^2 - 1/4.
- One gets S(e1 l)^{1+j e2} S(e2 l)^{1+j e1} = S^{2+j} e^{-j t S'/S} [1 + eps^2 Phi + O(eps^4)], with Phi = (1 + j/2) t^2 f'' - (j/6) t^3 f''' and f = log S.
- Order eps^0 is Theorem B. Order eps^2 gives C1, after simplifying with S^2 = 1 + t^2 f''.
- `lam1_rec.pl 16` computes Xi^(1) from REC + Mumford alone (no Hodge engine, no KP, no (6.1)). It equals C1 as a polynomial in j for h = 1..16 and equals the Hodge tables for h = 1..9. (Genus 1 includes the delta_irr term 1/2.)

### 3.2 Unknowns left
In the range, sigma^h Xi_h(j) lies in span(sigma, ..., sigma^{h-2}). These are the h-2 components k = 2..h-1; before this round there were h-1.

### 3.3 KP plus Theorem H plus K1
`tlr_k01.pl` is Route A's `a3_kp_nullity2.pl` with the constraints c_{t,G} = 0 (Theorem H) and c_{t,G-1} = 0 (K1) added; set the environment variable K01 = 0, 1 or 2. Nullities mod p:

| G | two-leaf subsystem: none / +H / +H+K1 | full system: none / +H / +H+K1 |
|---|---|---|
| 3 | 7 / 4 / 1 | 6 / 3 / 1 |
| 4 | 11 / 6 / 1 | 6 / 1 / **0** |
| 5 | 11 / 4 / **0** | 0 / 0 / 0 |
| 6 | 7 / **0** / 0 | – |

Consequences:
- Lemma A(4) follows from Lemma A(<= 3), KP, Theorem H and K1, without exact genus-4 Hodge data.
- The two-leaf rigidity TLR(G) with these constraints holds for 5 <= G <= 20; without them it held for 7 <= G <= 20.
- At G = 3 one unknown remains, the k = 2 component of Xi_3(4). It is the only exact-data input left for Lemma A with h <= 20 (h = 2 is proved by hand).

### 3.4 (6.1) at genus 7, 8, 9
- `c61.pl` compares both sides as polynomials in (j, v).
- Genus 7: referee table, 230 comparisons. Genus 8: new table, 286 comparisons. Genus 9: new table, 348 comparisons. 0 mismatches.
- The genus-8 and genus-9 one-leaf tables were computed here with the referee engine (`t_moments.pl`; ch_9 enabled at genus 9 with REF_NCH=5; 1320 s and 4927 s).
- (6.1) is now verified for h <= 9.

### 3.5 Other checks

| Check | Script | Range | Result |
|---|---|---|---|
| M3 and its extension to j = 2g, from REC data only | `hbdry_rec.pl` | g <= 12 | holds at j = 1..2g |
| y-Taylor coefficients of the two-leaf REC vs the engine | `bc_fit.pl` | c <= 4, g <= 7 | equal |
| Moments framework (one-leaf K-twisted moments for every ch-monomial) vs the component tables | `moments.pl`, `chk_moments.pl`, `Pair.pm` | g <= 6 | all equal |

## 4. Line (i): what jump loci and excess intersection can and cannot do

1. **Split at r = 1.**
   - Write Xi^(k)_h(j) = Sum_m Int c_m(F-P) [p_j e_k]_{N-m}, with F - P = -R pi_* S_{jx} of rank N - j and N = 3h-2.
   - Call the part with m > N - j "excess" (E^(k)). It is localized on P(ker phi), as in Theorem H.
   - Call the rest "non-excess" (NE^(k)).
2. **Parity at k = 1.** One has NE^(1)_h(j) = (-1)^{j-1} Xi^(1)_h(phi_{j-1}) and E^(1)_h(j) = Sum_{a<=j-2} binom(j-1,a) Xi^(1)(phi_a) + (1 + (-1)^j) Xi^(1)(phi_{j-1}). So a pure localization argument at r = 1 can only yield the even j, i.e. the odd levels D, which are exactly the (OL) levels.
3. **N1, k >= 2.** The excess part does not vanish in the range.
   - Genus 3: E^(2)(4) = 1/360 = -NE; E^(3)(4) = -1/1728.
   - Genus 4: E^(2)(4,5,6) = -1/1260, 1/1728, 3/128; E^(3)(4,5,6) = -61/272160, 677/2177280, -27617/1088640; E^(4)(5,6) = 127/7257600, 12749/3628800.
   - Hence a jump-locus induction with lattice Z-factors cannot work at r = 1. The Z-factors carry psi-powers at the node from the normal bundle (-psi - psi'), which are not lattice leaves.
   - With unit legs (Route C, Prop. L) the integral is fully localized, but the self-referential rational tails (s = z^{1-j} d/dz) remain. This is Route C's obstruction.
4. **N2.** Single ch_3 and ch_5 insertions do not localize to the irreducible boundary in the range (genus 3-6). Only the CY combination e_k(alpha + alpha^2) vanishes.
5. **N3, N4.** The CY analogue of M3, i.e. [y] Xi^(k)(j,y) equal to the separating sum with the Hodge class distributed, fails for every k >= 1. The naive analogue of M3 for [y^c] Y_g, c >= 2, fails too. So K1's mechanism does not extend term by term.
6. **K2, the exact reduction of k = 2 (`k2_reduce.pl`, exact at genus 2-6, all j).** Mumford gives e_2(alpha + alpha^2) = lambda_1^2/2 - 6 ch_3. Then

       Xi^(2) = (1/24) <lambda_1 (kappa_1 - psi + delta)> + (1/120)[<kappa_3> - <psi^3> + (1/2) <iota_*(psi'^2 - psi' psi'' + psi''^2)>],

   where <.> = Int c(T(-log)) (.) lam_j. The pieces are:
   - <kappa_3> = (2h-1) X_h - Sum_{c<=3} [y^c] Y_h;
   - <lambda_1 kappa_1> = (2h-1) Xi^(1) - [y] Xi^(1)(j,y);
   - the irreducible boundary of ch_3 is (1/2)(2 [y^2 z^0] - [yz]) Xi^(0)_{h-1}(j,y,z);
   - the separating boundary is Sum_{g1} ([y^2] Y_{g1} X_{g2}[0] - [y] Y_{g1} X_{g2}[1] + Y_{g1}(j,0) X_{g2}[2]).

   The two-leaf lambda_1 term reduces in the same way to Hurwitz data, through 12 [y] Xi^(1)_h(j,y) = [y] of <kappa_1 - psi_1 - psi_2 + delta> on Mbar_{h,2}. That term needs D_h := [yz] Xi^(0)_h(j,y,z) and [y^2] Y_h (checked at genus 2-5).

   With K1 and M3, the k = 2 statement in the range becomes a single identity:

   **H_2:** 5 D_h - 7 [y^2] Y_h - 12 [y^3] Y_h = an explicit expression in lower-genus data, for 2 <= j <= 2h-2.

   All quantities in H_2 are determined by Route B's REC with at most three leaves.

## 5. Line (ii): localization and the Faber-Pandharipande analogy
- In our notation FP's Theorem 2 is the K-free single-Hodge one-leaf function: Int_{Mbar_{g,1}} Lambda^vee(1)/(1 - k psi) = k^{2g-2} [t^{2g}] S(t)^{1-k}.
- Its K-twisted version is Theorem B (S^2 H^j). The K-free CY version is the one-leg Marino-Vafa formula. (6.1) is the K-twisted CY version, and its large-sigma expansion reproduces Theorem B and C1.
- On CY3 targets, finite fixed-point localization produces only K-free vertices (Route D, D3.1, D3.2). The K-twist appears only after resumming degree-1 edges at the conifold point (analytic continuation). So no FP-type finite relation exists for (6.1) itself.
- **Equivariant reformulation (proved).** By Mumford, Lambda^vee(w) = w^{2h}/e(E (x) C_w). The CY integrand c(T(-log)) Lambda^vee(tau) Lambda^vee(-1-tau) is therefore sigma^{2h} times the zero-section localization of Tot(E (x) C_tau (+) E (x) C_{-1-tau}) with obstruction bundle T(-log) (x) C_1.
- **Cosection observation.** Multiplying abelian differentials gives a cosection v -> <v, ab> (Serre duality H^1(T_C(-x)) x H^0(omega^2(x)) -> C). Its weight is 1 + tau - 1 - tau = 0, so it is torus-invariant exactly because of the CY condition.
  - With the leaf, use V_D and E(-alpha x), E(-beta x) with alpha + beta = j - 1.
  - Kiem-Li then moves the virtual class to {ab = 0}. This locus still contains the zero section and the coordinate subbundles, so no vanishing follows.
  - Recorded as an open route.

## 6. Line (iii): Givental structure and the genus independence in SN
**Proved equivalent forms of (6.1) and of Lemma A for r = 1:**
1. With W_h(z) := Sum_j sigma^h Xi_h(j) z^j, (6.1) is equivalent to W_h(z) = (2h-2)! [l^{2h}] S1 S2 (1 + (zB/(1-zB))^{2h-1}). Lemma A for r = 1 is equivalent to W_h(z) - W_h(0) = O(z^{2h-1}).
2. Sheffer form: Sum_h Q_h(j) mu^{2h} = S1 S2 B^j, with Q_h = sigma^h Xi_h/(j-1)_(2h-2). Equivalently Q_h(j+1) - Q_h(j) = Sum_{m>=1} b_m Q_{h-m}(j), where Sum_m b_m mu^{2m} = B. The j = 0 case is the gap formula, fixed by the conifold polar part.
3. Unit-leg dressing: e^{beta S} = (1 - y)^{-S} with y = 1 - B^{-1}. Hence sigma^h <Omega K f>_{h,1} = [l^{2h}] S1 S2 b_h^{-1} Sum_m ((1 - B^{-1})^m/m!) <lambda_h K f lam_1^m>_{h,1+m}. The CY one-leaf vertex is the lambda_g vertex dressed by unit legs, each weighted by 1 - e^{-beta(l)}.

**Genus independence.** At order 1/sigma, the genus independence of the series in (6.1)/(SN) is now explained in the range.
- The first-order R-matrix action is Mumford's lambda_1 relation. Acting on the Hurwitz theory, it is controlled by the genus-uniform Toda recursion (REC), whose generating-function solutions are S, H and the Bernoulli series X.
- The vanishing at this order is K1. The closed form C1 is verified to genus 16.

**Beyond first order:**
- Not done. The higher Mumford relations do not localize term by term (N2-N4).
- Route E's two-leaf NR form stays consistent through genus 5. Re-running it at genus 6-7 on the referee's two-leaf tables (`nr_test/`) did not finish in 50 minutes; the full fit reached genus 4 and is consistent there.

## 7. Uniformity in genus
- **K1, K1', M1-M3, E:** proved for all genera. The induction on h uses only REC (all genera, Route B), Theorem B (all genera) and Mumford's lambda_1 relation (all genera). Claim E is proved by generating functions uniformly in (h, j). No hypothesis on lower genera or other components is used.
- **Range-limited statements (computational):** C1 (h <= 16), M3' (g <= 12), (6.1) (h <= 9), KP tables (G <= 6), N1-N4 and K2 checks (the ranges listed).

## 8. Exact obstruction
1. (OL) at genus G and odd level D needs Xi^(k)_G(D+1) = 0 for 2 <= k <= G-1. After K1 these are G-2 numbers per (G, D).
2. K1's mechanism uses one Mumford relation, whose kappa-part is a single zero-leaf derivative. That derivative is a boundary sum in the low range (M3), and the sum closes under REC.
3. For k >= 2, the needed Hurwitz derivative data ([y^c] Y for c = 2, 3, three-leaf [yz]) are not boundary sums (N4). Their REC recursions are not closed on the combination in H_2. So an induction needs low-range closed forms for each of them separately, and none is known.
4. Geometric localization fails at r = 1 for k >= 2 (N1). With unit legs it runs into Route C's self-referential tails.

## 9. Next steps
1. Prove H_2. Derive the REC recursions for [y^2] Y, [y^3] Y (two leaves) and [yz] Xi^(0)(j,y,z) (three leaves). Look for low-range formulas mixing boundary sums with psi-insertions at the node, then prove an analogue of Claim E.
2. Extend K1 to r >= 2 with REC for R and a zero leaf. This would replace Route C's (T1) for the whole lambda_1-component.
3. For R1: a structural proof of TLR with the two extra zero coefficients per one-leaf polynomial (unknowns sigma, ..., sigma^{G-2}). Rigidity now starts at G = 5.
4. The next-to-lambda_g component, through Int K lambda_{h-1} lam_j. Proving it would remove the last genus-3 data point and leave k = 2..h-2.

## 10. Files (all in R4)
- **Libraries:** `RefHodge.pm`, `RefWK.pm`, `RefParse.pm` (copies of the referee engine); `Pair.pm` (Hodge-class pairing); `PJ.pm` (polynomial and series algebra).
- **New tables:** `refdata/ref_engine_h8_r1.txt`, `refdata/ref_engine_h9_r1.txt` (one-leaf components, genus 8 and 9); `refdata/moments_h2..7.txt`.
- **Scripts:**
  - (6.1) and tables: `c61.pl`, `t_moments.pl`, `moments.pl`, `chk_moments.pl`.
  - K1 and its proof: `p1moments.pl`, `rec2_test.pl`, `lam1_rec.pl`, `delta_bdry.pl`, `hbdry_rec.pl`, `claimE.pl`, `claimE_proof_checks.pl`.
  - Line (i) diagnostics and k = 2: `excess_split.pl`, `irr_test.pl`, `cyH_test.pl`, `bc_fit.pl`, `k2_reduce.pl`.
  - KP with constraints: `tlr_k01.pl`.
  - NR tests: `nr_test/` (fitNR.pl and copies of the two-leaf tables, genus 1-7).
- **Outputs:** `c61_out.txt`, `c61_h8_out.txt`, `c61_h9_out.txt`, `g8_r1.log`, `g9_r1.log`, `p1moments_out.txt`, `rec2_test_out.txt`, `lam1_rec_h8.txt`, `lam1_rec_h16.txt`, `lam1_rec_h9_vs_engine.txt`, `delta_bdry_out.txt`, `hbdry_rec_out.txt`, `claimE_h12.txt`, `claimE_proof_checks_out.txt`, `excess_split_h3.txt`, `excess_split_h4.txt`, `irr_test_out.txt`, `cyH_test_out.txt`, `bc_fit_out.txt`, `k2_reduce_out.txt`, `tlr_k01_twoleaf.txt`, `tlr_k01_full.txt`, `nr_test/nr_diag_h7.txt`, `nr_test/nr_full_h7.txt`.
