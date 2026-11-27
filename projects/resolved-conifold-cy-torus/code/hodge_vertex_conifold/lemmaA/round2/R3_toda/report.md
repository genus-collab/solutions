# Lemma A round 2: route R3 (saved by the main session from the workflow journal, 2026-10-31)

Subagents cannot write report files, so this file holds the agent's structured return value: the fields first, then the full report text (report_text) verbatim.

## all_genus_claim

1

## all_genus_statement

Lemma A (M_irr(h)) holds for every h >= 2: Xi_h(R) = 0 identically in tau for every multiset R of integers >= 2 with D(R) <= 2h-3. Consequently (OL) holds for every G >= 3 and every odd or even D <= 2G-3, and so does Corollary R (all-genus regularity of the conifold Hodge-vertex graph sum).

The proof is by induction on h and on the level D:
- The first modified-KP equation of the 2D Toda hierarchy of the Marino-Vafa content-product tau function, at the dilaton point, gives Theorem M.
- Under the inductive hypotheses every term of Theorem M except -2 Pi_2 sigma^h tau delta Xi_h(R+2) vanishes, by falling factorials or by induction. So every level-D value with a part 2 vanishes.
- Route A's Theorem B (cut-and-join) or Cor. E gives the remaining profiles. Alternatively, the u_k-linear members of the same hierarchy give every part k directly.

Status: proved modulo the MV formula and the 2D Toda / modified-KP equations for content-product tau functions; both were verified exactly on the MV Schur side.

## claims

### [1]

#### depends_on

MV formula (LLZ/OP); first-time-set (modified-KP) equations of the 2D Toda hierarchy for content-product tau functions (Takasaki 1984, Ueno-Takasaki 1984, Jimbo-Miwa 1983, Orlov-Scherbin 2001, Okounkov math/0004128 s.2.7); Route A Prop. 2.1/Cor. 2.2/Lemma 3.1 (proved); proof (i) also Route A Theorem B (cut-and-join) or Cor. E

#### evidence

complete argument, report s.2-5; induction on h and level D uses only (H1) lower genera and (H2) lower levels; all identities used verified (report s.7)

#### statement

Theorem L: Lemma A (M_irr(h)) holds for every h >= 2: Xi_h(R) = 0 identically in tau for every multiset R of integers >= 2 with D(R) <= 2h-3. In particular (OL) holds for every G >= 3 and every D <= 2G-3 (odd and even), and Corollary R (all-genus regularity) follows.

#### status

proved-modulo-literature

### [2]

#### depends_on

I1-I4 (MV, 2D Toda mKP equation, dilaton calculus)

#### evidence

374 relations genus 0-6 at 3 tau; 182 exact polynomial identities in tau (2G+3 points) genus 2-6; 60 on Mathematica q-tables genus 2-4; mod p with MV-side values at genus 7 (8), 8 (4), 9 (2); Hurwitz limit 314 relations to genus 10; lambda_g limit proved for all genera; perturbation controls fail

#### statement

Theorem M (2-insertion relation): for all G >= 0, eps in {0,1}, n = 2G+eps and all R, the w^{n+r} coefficient (dilaton point) of the first mKP equation is: -2 Pi_2 Sum_{2g+k=n+1} sigma^g tau^k (delta_g(R+2))_k X_g(R+2)/k! + [eps=0] 2 sigma^{G+1} J(J+1) X_G(R) + Sum_{2g+k=n,k>=1} sigma^{g+1} tau^k J(J+1)(delta_g(R))_k X_g(R)/k! + Sum_{S subset R} Sum_{Sum(2g_i+k_i)=n+2} Phi^{g1}_{k1}(R_S) Phi^{g2}_{k2}(R-R_S) = 0, with Phi^g_k(T) = sigma^g tau^k J_T (delta_g(T))_k X_g(T)/k! and the explicit empty-profile values (report s.4).

#### status

proved-modulo-literature

### [3]

#### depends_on

Theorem M

#### evidence

report s.5: falling-factorial dichotomy and count Sum_i(delta_i - k_i) = delta_G(R) <= -2

#### statement

Corollary I_2: under (H1) Lemma A in genera 2..G-1 and (H2) Xi_G(T) = 0 for D(T) < D (1 <= D <= 2G-3), Xi_G(R+2) = 0 for every R at level D-1; the only surviving term of Theorem M is -2 Pi_2 sigma^G tau delta_G(R+2) Xi_G(R+2), delta <= -1.

#### status

proved

### [4]

#### depends_on

2D Toda bilinear identity at equal second times (form fixed by direct checks); I1, I4

#### evidence

budget argument (every block contributes >= 0) + count Sum_B(delta_B - k_B) = D(R)+k-n; M_3 verified genus 0-4 (162 relations), M_4 genus 0-3 (136); P_3, P_4 verified on the Schur side

#### statement

Theorem M_P and Corollary I_P: for any member P(D) tau_0.tau_1 = 0 of the hierarchy of weight k with non-zero coefficient c of D_{t_k} (e.g. P_k = D_{t_k} + p_k(-D~)), the dilaton-point leading order obeys the factor rule of report s.5, and under (H1),(H2) Xi_G(R+k) = 0 for every R at level D-k+1 (target coefficient -c k sigma^G tau Pi_k delta). This gives a second proof of Theorem L using only the 2D Toda hierarchy.

#### status

proved-modulo-literature

### [5]

#### depends_on

2D Toda (Jimbo-Miwa mKP) for content-product tau functions

#### evidence

exact Schur-side check (mvlib): tau = 2/7, all monomials of degree <= 5 through x^30; tau = -5/3, degree <= 4 through x^34; P_3 through x^14 (degree <= 3), P_4 through x^12 (degree <= 2)

#### statement

The first mKP equation 2 d_{p2}(F - F1) + d_{p1}^2(F + F1) + (d_{p1}(F - F1))^2 = 0, F = log Z_MV, F1(p) = F(e^{tau x} o p), holds for the Marino-Vafa generating function; the opposite sign fails.

#### status

computer-verified

### [6]

#### evidence

both terms equal 2 b_G J(J+1) prod_{m=1}^{2G-4+r}(J-m) (J+3-2G-r)

#### statement

lambda_g limit of Theorem M: the O(tau) part is -2 delta_G(R+2) Lambda_G(R+2) + 2J(J+1) Lambda_G(R) = 0 with Lambda_G(T) = b_G (J_T-1)...(J_T-N_T); it holds for every G and R.

#### status

proved

### [7]

#### evidence

eq_bosonization_check.pl: all s_nu with |nu| <= 5 at two values of q, branching rule included (29 checks); toda_eq_genus0_check.pl through y^8 at tau = 2/7, -5/3

#### statement

Line (2): referee B's operator is explicit: Eq = (1/vs_1)[w^0](exp(-Sum vs_k w^k p_k/k) exp(-Sum vs_k w^{-k} d_{p_k}) - 1), vs_k = q^{k/2} - q^{-k/2}, so e^{-H}Eq e^H is a q-Miwa-shift residue; at genus 0 and empty profile the CY Toda equation reads sigma d_y [w^0] e^{-Y(w)} = e^{tau T(y)/(1+tau)}, so at the dilaton point it carries the non-integer exponential sector e^{tau u/(1+tau)} and has no termwise leading-order analysis (all completed cycles contribute at genus 0).

#### status

computer-verified

### [8]

#### evidence

thmB_check.pl: 32 relations (g = 2,3) + 90 (g = 4-6), r <= 2

#### statement

Route A's Theorem B (2-insertion relation from cut-and-join) holds exactly on the engine tables at genus 2-6.

#### status

computer-verified

### [9]

#### evidence

48 values at two tau, 0 mismatches (mkp_lead_qdata.pl)

#### statement

The Mathematica q-tables (lemmaA\data) and Route E's engine agree on Xi_g, g = 2-4.

#### status

computer-verified

### [10]

#### evidence

structural remark (Route A s.3.8: simple zeros at s = -1); not pursued since Theorem M suffices

#### statement

Line (3): the cut-and-join (framing) flow at the next order in (1+s) does not close on Xi: it introduces Q'_h(-1;R) = <Omega K kappa_1 lam_R>, which do not vanish; it is not needed.

#### status

conjectural


## key_formulas

- First mKP equation (2D Toda, fixed p' = q^rho): 2 d_{p2}(F - F1) + d_{p1}^2 (F + F1) + (d_{p1}(F - F1))^2 = 0, F = log Z_MV, F1(p) = F(e^{tau x} o p)
- Toda shift: A = F1 - F = Sum_{k>=1} (tau x)^k E^k F/k!, E = Sum_k k p_k d_{p_k};  S = 2F + A;  Hirota blocks: odd -> -d^B A, even -> d^B S
- Dilaton point: y = u e^{-u}, w = 1/(1-u) = 1+s;  e^{-uJ}(J+D)e^{uJ} w^n = (J-n) w^{n+1} + n w^{n+2},  e^{-u(J+1)} d_y e^{uJ} w^n = J w^{n+1} + n w^{n+2}
- Leading coefficient: d_y^v (J_T+D)^k G^(g)_T ~ e^{u(J_T+v)} (J_T)^(v) (delta_g(T))_k Xi_g(T) w^{2g-2+r_T+k+v},  delta_g(T) = D(T) - 2g + 2
- Theorem M: -2 Pi_2 Sum_{2g+k=n+1} sigma^g tau^k (delta_g(R+2))_k X_g(R+2)/k! + [n even] 2 sigma^{G+1} J(J+1) X_G(R) + Sum_{2g+k=n,k>=1} sigma^{g+1} tau^k J(J+1)(delta_g(R))_k X_g(R)/k! + Sum_{S subset R} Sum_{Sum(2g_i+k_i)=n+2} Phi^{g1}_{k1}(R_S) Phi^{g2}_{k2}(R-R_S) = 0,  Phi^g_k(T) = sigma^g tau^k J_T (delta_g(T))_k X_g(T)/k!,  Phi^0_1(empty) = tau, Phi^0_2(empty) = tau^2/2
- Induction step (n = 2G, R at level D-1): only -2 Pi_2 sigma^G tau delta_G(R+2) Xi_G(R+2) survives; general member: -(k-1) sigma^G tau Pi_k delta_G(R+k) Xi_G(R+k)
- Count: Sum_B (delta_B - k_B) = D(R) + k - n = delta_G(R+k) - 1 <= -2;  budget: every Hirota block contributes 2g - 2 + |B| + k_B >= 0
- lambda_g limit: delta_G(R+2) Lambda_G(R+2) = J(J+1) Lambda_G(R),  Lambda_G(T) = b_G (J_T-1)...(J_T-N_T)
- Eq = (1/vs_1)[w^0](exp(-Sum vs_k w^k p_k/k) exp(-Sum vs_k w^{-k} d_{p_k}) - 1),  vs_k = q^{k/2} - q^{-k/2};  genus 0: sigma d_y [w^0] e^{-Y(w)} = e^{tau T(y)/(1+tau)}

## next_steps

- Referee pass on the new step: derive Theorem M from P_2 and the dilaton calculus (Lemmas 2.1-2.3, especially the empty-profile values), then check Corollary I_2 term by term ((a), (b), (c), (d)) and the budget and count argument of Corollary I_P.
- Main session: run R3_thmM_check.wls, a symbolic-in-tau Mathematica check of Theorem M at genus 2-4.
- Main session: record Lemma A for all h, and update leaf_structure_proof.md s.13, conifold_regularity_note.md and proof_status_mixed_scale.md. Corollary R (all-genus regularity) is now proved modulo the MV formula and the 2D Toda hierarchy; Route A's TLR/Theorem F, Route D's closed form and the parity argument become unnecessary.
- Optional: an exact genus-7 check of Theorem M with Hodge-side values. Compute Xi_7(13), Xi_7(2,13), Xi_7(14), Xi_7(2,14) with the Route E engine on Stanage (g7_points.pl, HODGE_NCH = 4) or in Mathematica (hvc_lib), then evaluate M_{7,eps}((13)) with mkp_lead.pl.
- Optional: use Theorem M_P as a recursion. Away from delta = 0 it determines every Xi_G(T) from lower genera and the excess-0 values. This could lead to closed formulas off the range, for example (6.1) or the NR/SN forms of Route E.

## obstructions

- None for Lemma A, modulo the cited literature. The status is proved-modulo-literature because the MV formula and the 2D Toda / modified-KP equations for content-product tau functions are cited. Their normalisations were verified directly on the MV Schur side: P_2 through x^30 / x^34, P_3 through x^14, P_4 through x^12.
- Proof (ii) uses the u_k-linear members for every k, taken from the 2D Toda bilinear identity at equal second times. Conventions differ between references, so the form was fixed by direct checks for k <= 4. The same argument works with the standard DJKM members that have a non-zero linear coefficient. Proof (i) needs only P_2 (Jimbo-Miwa's first mKP equation) plus Route A's Theorem B.
- Theorem M relies on Route A's Prop. 2.1 / Cor. 2.2 and on the explicit unstable functions: the empty-profile table of Lemma 2.3. A referee should re-derive that table (in particular that no factor is more singular than its x-order) and the counting in Corollaries I_2 and I_P.
- The genus 7-9 checks of Theorem M are mod p (MV Schur side). Exact identity checks reach genus 6. The exact Hodge-engine evaluation of genus-7 points did not finish locally in 2 h.
- Line (2), the p'-direction Toda equation with Eq, admits no termwise dilaton-point analysis: the right side has the non-integer exponential sector e^{tau u/(1+tau)}, and all completed cycles contribute at genus 0. Line (3), the next order of the cut-and-join, brings in the new unknowns Q'(-1). Neither is needed.

## requested_checks

- Mathematica (main session): run 'wolframscript -file lemmaA/round2/R3_toda/R3_thmM_check.wls' from the project folder. This is a symbolic-in-tau check of Theorem M on the q-tables for G = 2-4 and R = (j), j = 2..12, both parities. Expected: every line 0, 'ALL ZERO: True', and a non-zero control.
- Referee: re-derive Theorem M from the first mKP equation, with Lemma 2.3's empty-profile values and the claim that no term is more singular than its x-order; re-check Corollary I_2 term by term and the budget and count argument of Corollary I_P.
- Optional exact genus-7 test: on Stanage run 'HODGE_NCH=4 perl g7_points.pl 7 "13;14;2,13;2,14"', or compute these values with hvc_lib in Mathematica; then evaluate M_{7,eps}((13)) and M_{7,eps}((14)) exactly. The mod-p versions already pass for two primes.

## route

R3: the 2D Toda hierarchy of the Marino-Vafa content-product tau function, used through its first-time-set (modified-KP) equations linking tau_n(p) = Z(e^{n tau x} o p) and tau_{n+1}, analysed at the dilaton point (Route B's Toda + T -> oo method generalised to the full CY vertex)

## summary

Lemma A, M_irr(h), is proved for every h >= 2, modulo the cited literature. The new input is the first modified-KP equation of the 2D Toda hierarchy for the pair tau_0 = Z_MV(p), tau_1 = Z_MV(e^{tau x} o p) at fixed second times p' = q^rho: 2 d_{p2}(F - F1) + d_{p1}^2 (F + F1) + (d_{p1}(F - F1))^2 = 0. The Toda shift enters only through A = F1 - F = Sum_k (tau x)^k E^k F/k!. At the dilaton point (w = 1+s -> 0) the Euler operator E raises the w-order by one with coefficient delta = D - 2g + 2, so E^k produces Route B's falling factorials (delta)_k. The leading order of the equation is an exact all-genus relation (Theorem M, the "2-insertion relation") among the values Xi_g. Its only genus-G terms are -2 Pi_2 sigma^G tau delta_G(R+2) Xi_G(R+2) and 2 sigma^{G+1} J(J+1) Xi_G(R), one level lower. Every other term carries a falling factorial or is a product, and all of them vanish in the range by Route B's count. So, under the hypotheses of (OL), every level-D value with a part 2 vanishes (Corollary I). Route A's Theorem B (cut-and-join), or Cor. E, then gives every T at level D, so (OL) holds for all G >= 3 and every odd or even D <= 2G-3, and with it Lemma A. A second proof uses only the 2D Toda hierarchy: the u_k-linear members P_k give k-insertion relations with target coefficient -(k-1) sigma^G tau Pi_k delta. TLR, Theorem F, the closed form and the parity argument are not needed. With Lemma C and Theorem R, the all-genus regularity theorem (Corollary R) follows in every genus.

All checks are exact unless marked mod p, and no topological recursion was used anywhere.
- Theorem M against Route E's engine tables: genus 0-6, 374 relations at 3 values of tau.
- Theorem M as a polynomial identity in tau (2G+3 points, so these cases are proofs): genus 2-6, 182 relations.
- Theorem M on the Mathematica q-tables (lemmaA\data), genus 2-4: 60 relations.
- Theorem M with genus-7, 8 and 9 values from the MV Schur side, mod two primes: 14 relations.
- Hurwitz limit to genus 10 (314 relations); the lambda_g limit is proved for all genera.
- The mKP equation itself, directly on the Schur side, through x^30 and x^34; the generated 2D-Toda members P_3 and P_4 also checked.
- Theorem M_3 to genus 4 (162 relations) and M_4 to genus 3 (136); Theorem B re-checked at genus 2-6 (122 relations); all perturbation controls fail as they must.

Lines (2) and (3) are explained but not needed. Referee B's p'-direction Toda equation was made explicit (Eq as a q-Miwa-shift residue, verified). At the dilaton point it carries a non-integer exponential sector e^{tau u/(1+tau)}, so it admits no termwise leading-order analysis. The next order of the cut-and-join brings in new unknowns.

## report_text

# Round 2, R3 (Toda/KP): Lemma A for all genera from the modified-KP equations of the 2D Toda hierarchy at the dilaton point

Folder: HVC\lemmaA\round2\R3_toda (HVC = G:\My Drive\Projects\Old\refGW\hodge_vertex_conifold); outputs in out\. Arithmetic is exact rational (Perl Math::BigRat), except the genus 7-9 tests, which are mod p. No wolframscript was run.

No topological recursion of any kind was used: no spectral curve, loop equation, Eynard-Orantin invariant, BKMP, dictionary P2, mod-P2 result or TR script, not even as a cross-check.

## 0. Results

**Theorem L (Lemma A = M_irr(h) for every h).** For every h >= 2 and every multiset R of integers >= 2 with D(R) <= 2h-3, Xi_h(R) = 0 identically in tau.
Status: proved modulo literature.

Inputs:
- (I1) The Marino-Vafa formula (Liu-Liu-Zhou, JDG 65 (2003) Thm 1; Okounkov-Pandharipande, GT 8 (2004)), in Route A's normalisation (1.1).
- (I2) The content-product (hypergeometric) 2D Toda tau functions tau_n(p,p') = Sum_nu e^{tau x (f_2(nu) + n|nu|)} s_nu(p) s_nu(p') (Takasaki 1984; Orlov-Scherbin, TMP 128 (2001); Okounkov math/0004128 s.2.7). Hence tau_n(p, q^rho) = Z(e^{n tau x} o p).
- (I3) The first-time-set equations of this 2D Toda lattice at fixed second times p' = q^rho, i.e. the modified-KP hierarchy for the pair (tau_0, tau_1) (Jimbo-Miwa, Publ. RIMS 19 (1983); Ueno-Takasaki, ASPM 4 (1984)).
- (I4) Route A's Prop. 2.1 / Cor. 2.2 and Lemma 3.1 (dilaton-point calculus; proved).
- (I5) Route A's Theorem B (MV cut-and-join at the dilaton point; proved modulo LLZ). Used only in proof (i).

Two proofs:
- (i) The first mKP equation (Theorem M, 2-insertion) together with Theorem B.
- (ii) All u_k-linear members of the same hierarchy (Theorem M_P, k-insertion), using (I1)-(I4) only: no Theorem B, no KP of Route A, no closed form.

Consequences:
- (OL) holds for every G >= 3 and every D <= 2G-3, odd and even.
- TLR(G), Theorem F, Route D's closed form and the parity reduction (brief (f)) are no longer needed.
- With Lemma C and Theorem R (leaf_structure_proof.md s.13), Corollary R (all-genus regularity of the conifold Hodge-vertex graph sum) holds in every genus.

**Mechanism.** The shift tau_1 = tau_0(e^{tau x} o p) enters only through A := log tau_1 - log tau_0 = (e^{beta E} - 1) F = Sum_{k>=1} beta^k E^k F/k!, with beta = tau x and E = Sum k p_k d_{p_k}. At the dilaton point E raises the w-order by one with coefficient delta = D - 2g + 2, so E^k produces the falling factorial (delta)_k. This is exactly Route B's mechanism, now for the full CY class. The leading order of the equation contains the genus-G unknown at level D only through -2 Pi_2 sigma^G tau delta_G(R+2) Xi_G(R+2). Every other term is one level lower, has lower genus with a falling factorial, or is a product killed by Route B's count.

## 1. Conventions

**Basic quantities.** c = 1, sigma = tau(1+tau), Pi_k = prod_{a=1}^{k-1} (k tau + a)/(k-1)!.
- For a multiset T of integers >= 2: r_T, J_T, D_T = J_T - r_T.
- X_g(T) := Xi_g(T) (brief s.1) for g >= 1.
- Genus 0: X_0(T) = (J-1)(J-2)...(J-r+3) for r >= 3, 1/J for r = 2, 1/(j(j+1)) for r = 1 (Route A Lemma 3.1; Route B (U3), (U4)).
- delta_g(T) := D_T - 2g + 2. Then delta_g(T) <= -1 iff T is in the Lemma-A range at genus g. For T nonempty, delta_0(T) >= 3 and delta_1(T) >= 1.
- (a)_k = a(a-1)...(a-k+1) (falling factorial); (a)^(u) = a(a+1)...(a+u-1) (rising factorial).

**Generating function.**
- Z = Sum_nu e^{tau x f_2(nu)} s_nu(q^rho) s_nu(p), with q = e^x and p_k(q^rho) = 1/(q^{k/2} - q^{-k/2}). This is Route A's real form.
- F = log Z = Sum_g hbar^{g-1} F_g(q), with hbar = sigma x^2 and q_k = x sigma Pi_k p_k.
- F_g = Sum_T q_T G^(g)_T(y)/|Aut T|, with y = q_1 and G^(g)_T(y) = Sum_a M_g(T + 1^a) y^a/a!. All unstable terms are included, with M_0(mu) = |mu|^{l-3}.
- Dilaton point: y = u e^{-u} on the second sheet, w := 1/(1-u) (= 1+s of Route A), and u -> oo corresponds to w -> 0.

## 2. Dilaton-point calculus

**Lemma 2.1 (expansions; this is I4).**
- For T nonempty and every g >= 0: e^{-u J_T} G^(g)_T = w^{2g-2+r_T} * (a polynomial in w), with constant term X_g(T).
  - g >= 1: this is Route A Cor. 2.2, e^{-uJ}G = w^e Q_g(w-1;T) with Q_g(-1;T) = Xi_g(T).
  - g = 0: Route A Lemma 3.1, together with G^(0)_(j) = e^{ju}(1/j^2 - u/(j(j+1))) and G^(0)_(j1,j2) = e^{Ju}/J.
- Empty profile:
  - G^(0)_0 = u - 3u^2/4 + u^3/6 (tree function).
  - y d_y G^(1)_0 = (w-1)(w-1-1/sigma)/24 (Route A s.3.1).
  - For g >= 2, G^(g)_0 lies in w^{2g-2} Q[w], d_y G^(g)_0 = G^(g)_(1), and X_g(1) = 0 by the string equation.

**Lemma 2.2 (operators).** With D = y d_y = (w-1) w^2 d_w and 1/y = e^u w/(w-1):
- e^{-uJ}(J + D) e^{uJ} w^n = (J-n) w^{n+1} + n w^{n+2};
- e^{-u(J+1)} d_y e^{uJ} w^n = J w^{n+1} + n w^{n+2}.

Both operators raise the order by at least 1. When J = 0, d_y raises it by 2.

**Lemma 2.3 (leading coefficients).**
- For T nonempty: e^{-u(J_T+v)} d_y^v (J_T + D)^k G^(g)_T = (J_T)^(v) (delta_g(T))_k X_g(T) w^{2g-2+r_T+k+v} + higher powers of w.
- For T empty and v >= 1: d_y^v D^k G^(g)_0 = e^{vu} (L0 w^{2g-2+k+v} + higher powers of w), where L0 takes the values below. A value 0 means that the term is of strictly higher order.

| (g, k) | L0 |
|---|---|
| (0, 0) | (v-1)!/2 |
| (0, 1) | (v-1)! |
| (0, 2) | (v-1)! |
| (0, k >= 3) | 0 |
| (1, 0) | (v-1)! X_1(1), with X_1(1) = -(1+sigma)/(24 sigma) |
| (1, k >= 1) | 0 |
| (g >= 2, any k) | 0 |

These follow from Lemmas 2.1 and 2.2, for example d_y D G^(0)_0 = e^u, d_y D^2 G^(0)_0 = e^u w and D^3 G^(0)_0 = w - 1. Hence no function in the expansion is more singular than its x-order predicts.

## 3. The equations

**F/A/S form.** Write tau_0 = Z, tau_1 = Z(e^beta o p) with beta = tau x, F = log tau_0, A = log tau_1 - log tau_0 and S = 2F + A.
- For a Hirota polynomial P: P(D) tau_0.tau_1/(tau_0 tau_1) = P(d_s) e^{Phi}|_{s=0}, with Phi(s) = F(t+s) + F1(t-s) - F - F1.
- Hence d_s^B Phi|_0 = d^B S if |B| is even and -d^B A if |B| is odd (Faa di Bruno over set partitions).

**Members used.**
- P_k = D_{t_k} + p_k(-D~), with D~ = (D_{t1}, D_{t2}/2, ...) and t_k = p_k/k. P_k is the u_k-coefficient of the 2D Toda bilinear identity at equal second times, in the form Res_z z^{-1} tau_0(t+u-[1/z]) tau_1(t-u+[1/z]) e^{2 xi(u,z)} dz = tau_1(t+u) tau_0(t-u).
- P_2 and P_3 are Jimbo-Miwa's first two mKP equations: (D_1^2 - D_2) tau_1.tau_0 = 0 and (D_1^3 + 3 D_1 D_2 - 4 D_3) tau_1.tau_0 = 0.
- Any member with a nonzero coefficient of D_{t_k} works in Section 5. One example is the u_{k-2}-linear member (D_{t_{k-2}} p_2(-D~) + 2 p_k(-D~)) tau_1.tau_0 = 0 of the standard DJKM identity Res z tau_1(t-[1/z]) tau_0(t'+[1/z]) e^{xi(t-t',z)} dz = 0.

**Explicit forms** (generator mkp_eqs.pm):
- P_2: -2 d_{p2} A + d_{p1}^2 S + (d_{p1} A)^2 = 0, i.e. 2 d_{p2}(F - F1) + d_{p1}^2 (F + F1) + (d_{p1}(F - F1))^2 = 0.
- P_3: -12 d_{p3} A + 6 d_{p1} d_{p2} S + 6 d_{p1}A d_{p2}A + d_{p1}^3 A + 3 d_{p1}^2 S d_{p1}A + (d_{p1}A)^3 = 0.

**Direct verification on the Schur side of (I1)** (exact; the sign flip of A fails in every case):
- P_2 through x^30 for all monomials of degree <= 5 (tau = 2/7), and through x^34 for degree <= 4 (tau = -5/3).
- P_3 through x^14 (degree <= 3).
- P_4 through x^12 (degree <= 2).

Because conventions for the bilinear identity vary between references, these checks, not a transcription, fix the form used. Proof (i) needs only P_2, the most standard of these equations; Route B's referee had also confirmed it for the Hurwitz series.

## 4. Theorem M (2-insertion relation)

**Statement.** For all G >= 0, eps in {0,1}, n = 2G + eps, and every multiset R of integers >= 2 (J = J_R, r = r_R):

    -2 Pi_2 Sum_{g>=0, k>=1, 2g+k=n+1} sigma^g tau^k (delta_g(R+2))_k X_g(R+2)/k!
    + [eps = 0] 2 sigma^{G+1} B_G(R)
    + Sum_{g>=0, k>=1, 2g+k=n} sigma^{g+1} tau^k C^g_k(R)/k!
    + Sum_{S subset R labelled} Sum_{g_i>=0, k_i>=1, Sum (2g_i+k_i) = n+2} Phi^{g1}_{k1}(R_S) Phi^{g2}_{k2}(R - R_S) = 0.

The symbols are defined as follows.
- For nonempty arguments: B_G(R) = J(J+1) X_G(R), C^g_k(R) = J(J+1) (delta_g(R))_k X_g(R), and Phi^g_k(T) = sigma^g tau^k J_T (delta_g(T))_k X_g(T)/k!.
- B for the empty profile: B_0(0) = 1/2, B_1(0) = X_1(1), and B_G(0) = 0 for G >= 2.
- C for the empty profile: C^0_1(0) = C^0_2(0) = 1, and C^g_k(0) = 0 otherwise.
- Phi for the empty profile: Phi^0_1(0) = tau, Phi^0_2(0) = tau^2/2, and Phi^g_k(0) = 0 otherwise.

**Proof.**
1. Expand in x. In the q-variables, d_{p_k} = x sigma Pi_k d_{q_k}, beta^k = tau^k x^k and hbar^{g-1} = sigma^{g-1} x^{2g-2}. Therefore:
   - -2 d_{p2}A = -2 Pi_2 Sum x^{2g+k-1} sigma^g tau^k d_{q2}E^kF_g/k!;
   - d_{p1}^2 S = 2 Sum x^{2g} sigma^{g+1} d_y^2 F_g + Sum x^{2g+k} sigma^{g+1} tau^k d_y^2 E^k F_g/k!;
   - d_{p1}A = Sum x^{2g+k-1} sigma^g tau^k d_y E^k F_g/k!.
2. Extract the q_R/|Aut R| coefficients. For d_{q2} E^k F_g it is (J + 2 + D)^k G^(g)_{R+2}. For d_y^2 E^k F_g it is d_y^2 (J + D)^k G^(g)_R. For a product it is the sum over labelled splittings of R.
3. Pass to the dilaton point. At fixed (x^n, q_R) the equation is a finite identity between power series in y. By Lemma 2.1 every term is e^{u(J+2)} times a Laurent polynomial in w. The log in G^(1)_0 never occurs undifferentiated. Analytic continuation in u = T(y) therefore turns the identity into one between Laurent polynomials.
4. By Lemmas 2.2 and 2.3, every term has w-order at least n + r, and its w^{n+r} coefficient is the one written. This proves the theorem.

**Remarks.**
- The relation holds for every G, inside and outside the range. Outside the range its terms are nonzero, which is what the checks below test.
- tau -> 0 limit. The O(tau) part is -2 delta_G(R+2) Lambda_G(R+2) + 2 J(J+1) Lambda_G(R) = 0, with Lambda_G(T) = b_G (J_T - 1)...(J_T - N_T) (the lambda_g formula). This is an identity for every G, because both terms equal 2 b_G J(J+1) prod_{m=1}^{2G-4+r} (J - m) (J + 3 - 2G - r).
- tau -> oo limit. The top coefficient (Pi_2 ~ 2 tau, sigma^g ~ tau^{2g}) is the same relation for the Hurwitz components.

## 5. Insertion vanishing and the all-genus theorem

**Hypotheses at (G, D)**, with G >= 2 and 1 <= D <= 2G-3:
- (H1) Lemma A holds in genera 2, ..., G-1.
- (H2) X_G(T) = 0 whenever D_T < D.

**Corollary I_2.** Under (H1) and (H2), X_G(R+2) = 0 for every R with D_R = D - 1 (R = 0 if D = 1).

**Proof.** In M_{G,0}(R), the term (a) with g = G, k = 1 equals -2 Pi_2 sigma^G tau delta_G(R+2) X_G(R+2), and delta_G(R+2) = D - 2G + 2 <= -1. All other terms vanish:
- (a) with g < G, k = 2(G-g) + 1. Here delta_g(R+2) = delta_G(R+2) + 2(G-g) <= k - 2. If delta_g >= 0, the falling factorial (delta_g)_k is 0. If delta_g <= -1, then g >= 2 and R + 2 lies in the range at genus g, so the term vanishes by (H1).
- (b). X_G(R) is at level D - 1, so it vanishes by (H2). For R = 0 and G >= 2, B_G(0) = 0.
- (c), with g = G - m, k = 2m >= 2. Here delta_g(R) = delta_G(R) + 2m <= 2m - 2 < k, and the same dichotomy applies. For R = 0 the only possible term is (g,k) = (0,2), which occurs only at G = 1, so there is none.
- (d), the products.
  1. Set delta = 2 for an empty argument. Then Sum_i (delta_{g_i}(T_i) - k_i) = D_R + 4 - Sum (2g_i + k_i) = delta_G(R) <= -2.
  2. Empty factors are nonzero only with delta - k in {1, 0}.
  3. Hence some factor with T_i nonempty has delta < k, and it vanishes by the dichotomy.
  4. If that factor has genus G, then its partner is Phi^0_1(0), k_i = 1, and T_i is a subset of R. Its level is then <= D - 1, so it vanishes by (H2).

Since Pi_2 sigma tau is not identically zero, X_G(R+2) = 0.

**Theorem M_P (general member).** Let P be a Hirota polynomial of weight k with P(D) tau_0.tau_1 = 0. The w^{n+r} coefficient at (x^n, q_R) is

    Sum_{monomials c_alpha D^alpha} c_alpha Sum_{set partitions pi of alpha} prod_{B in pi} (-1)^{|B|} (prod_{i in B} i) Sum_{labelled assignments of R} Sum_{(g_B,k_B): Sum_B (2g_B - 2 + |B| + k_B) = n} prod_B f_B,

and it vanishes. The block factors are as follows.
- A block B of odd size is an A-factor, with k_B >= 1.
- A block of even size is an S-factor, with k_B >= 0 and an extra factor 2 when k_B = 0.
- f_B = sigma^{g-1+|B|} tau^k (prod_{i in B, i>=2} Pi_i) L/k!.
- L = (J_U)^(v) (delta_g(U))_k X_g(U), where U = T_B + {i in B : i >= 2} and v = #{i in B : i = 1}. If U is empty, L is the value L0 of Lemma 2.3.

The proof is the same as for Theorem M; for P_2 the formula reduces to Theorem M.

**Corollary I_P.** If the coefficient c of D_{t_k} in P is nonzero, then under (H1) and (H2), X_G(R+k) = 0 for every R with D_R = D - k + 1. The target coefficient is -c k sigma^G tau Pi_k delta_G(R+k); for P_k it is -(k-1) sigma^G tau Pi_k delta.

**Proof.**
- (alpha) Budget. Every block contributes 2g - 2 + |B| + k_B >= 0, because odd blocks have k_B >= 1 and even blocks have |B| >= 2; the total is n = 2G.
  - So a genus-G block has (|B|, k_B) in {(1,1), (2,0)}, and every other block contributes 0.
  - Its profile is at level <= D_R + k - 1 = D. Equality holds only for the linear term D_{t_k} with T_B = R, which is the target. Otherwise the block vanishes by (H2).
- (beta) Count. Sum_B (delta_B - k_B) = D_R + k - n = delta_G(R+k) - 1 <= -2, with delta_B := 2 - 2g_B for an empty U. Nonzero empty blocks have delta - k >= 0. So some nonempty block has delta_B < k_B and vanishes, as in the proof of Corollary I_2.

**Proof of Theorem L.** Induction on h, and for fixed h on D = 1, ..., 2h-3. The base case h = 2 (D = 1, R = 0) needs no hypothesis.
- (ii) Every T at level D has a part k >= 2. Corollary I_{P_k}, applied with R = T - k (level D - k + 1 < D), gives X_h(T) = 0.
- (i) Corollary I_2 gives every T that contains a part 2.
  - Under (H1), (H2) and I_2, Theorem B at R (level D - 1) reduces to Sum_i gamma_{j_i} X_h(R_i^+) = Sum_{i<k} gamma_{j_i j_k} X_h(R_ik). Its genus-(h-1) terms and product terms have the same excess >= 1 and vanish by (H1) or (H2).
  - Now induct on r_T and then on min T.
    - r_T = 1: taking R = (D) gives gamma_D X_h(D+1) = 0.
    - min T = m >= 3: take R = T - m + (m-1). Raising the unique part m-1 of R gives T. Raising any other part gives a profile containing m-1, which has a smaller minimum and is already known to vanish. The R_ik have fewer parts.
    - Hence gamma_{m-1} X_h(T) = 0, and gamma_j = (j+1) Pi_{j+1}/(Pi_j sigma') is not identically zero.
  - Alternatively one can use Route A's Cor. E, or the KP merge Pi_2 Pi_D X_h(2,D) = (D+1) Pi_{D+1} X_h(D+1).

## 6. Uniformity in the genus, and which step uses what

- Theorems M and M_P are identities for every G and every R, inside and outside the range. They are derived from (I1)-(I4) with no inductive hypothesis.
- Lemmas 2.1-2.3 are genus-independent, because Prop. 2.1 holds for all h.
- Corollaries I_2 and I_P use only (H1), for lower-genus blocks with delta <= -1, and (H2), for genus-G blocks of lower level. The counting identities hold uniformly in G.
- Theorem L at genus G needs only Lemma A in genera < G and the lower levels at genus G.
- No data enter any proof.

## 7. Exact checks (all zero; outputs in out\)

| identity | data | range | count |
|---|---|---|---|
| P_2 (first mKP equation) directly on the MV Schur side | mvlib (exact) | tau = 2/7: degree <= 5 through x^30; tau = -5/3: degree <= 4 through x^34 | 19 + 12 monomials; sign flip of A fails in all |
| P_3 and P_4 on the Schur side | mvlib | P_3: degree <= 3 through x^14; P_4: degree <= 2 through x^12 | 7 + 4 monomials; controls fail |
| Theorem M at 3 values of tau (2/7, 3/2, -5/3) | Route E engine tables (XiData.pm) | G = 0-2 (r <= 2, J <= 7), G = 3, 4 (r <= 3, J <= 10), G = 5, 6 (r <= 2, J <= 12) | 78 + 74x4 = 374 |
| Theorem M as a polynomial identity in tau (degree <= 2G+2, checked at 2G+3 points, so these cases are proved) | engine tables | G = 2, 3, 4 (r <= 2, J <= 9); G = 5, 6 (r <= 1, J <= 14) | 42x3 + 28x2 = 182 |
| Theorem M | lemmaA\data q-tables (Mathematica, Mumford formula) | G = 2-4, R = (j), j <= 11, 3 values of tau | 60; q-tables vs engine: 48 values, 0 mismatches |
| Theorem M at genus 7-9 | genus >= 7 values from the MV Schur side + Laguerre transform (referee's rmv.pm, mod p); genus <= 6 from the tables | G = 7: R = (13), (14) (p = 10^9+7) and (13), (15) (p = 998244353); G = 8: (15), (16); G = 9: (17); both parities | 14; perturbing Xi_G(R+2) by 1 gives a nonzero residual |
| Hurwitz limit of Theorem M | Route B's (REC) (Rec.pm) | G = 1-5 (r <= 2), 6-8 (r <= 2), 9-10 (r = 1) | 160 + 138 + 16 = 314 |
| lambda_g limit | lambda_g formula | all G | proved (Section 4) |
| Theorem M_3 | engine tables | G = 0-4 (r <= 2) | 78 + 84 = 162 |
| Theorem M_4 | engine tables | G = 0-3 (r <= 2, J <= 8) | 136 |
| Route A's Theorem B (re-check) | engine tables | g = 2-6 (r <= 2) | 32 + 90 = 122 |
| perturbation controls | | Xi_2(2,3) and Xi_1(2,3) in M; Xi_2(3,3) in M_3 | exactly the affected relations fail |
| Line (2): branching rule and bosonised Eq | exact symmetric functions | all s_nu with |nu| <= 4 (Q = 3/2) and <= 5 (Q = 5/7) | 29, 0 failures |
| Line (2): genus-0 empty-profile CY Toda equation | exact series | through y^8, tau = 2/7 and -5/3 | OK |

The exact Hodge-engine evaluation of genus-7 points (g7_points.pl, HODGE_NCH = 4) did not finish within 2 h locally; this is why the genus 7-9 tests use the MV Schur side mod p.

## 8. The three lines of the task, and comparison with Theorem D

**(1) 2D Toda for the family tau_n = Z(e^{n tau x} o p).** What works is the first-time-set (mKP) part at fixed p' = q^rho.
- Its dispersionless limit is the G = 0 part of Theorem M.
- In the empty profile it reads -2 tau Pi_2 d_{q2} E F_0 + 2 sigma d_y^2 F_0 + tau^2 (d_y E F_0)^2 = 0. Exactly: -tau Pi_2 e^{2u} + sigma e^{2u} + tau^2 e^{2u} = 0.

**(2) Referee B's CY Toda equation d_{p1}(e^{-H} Eq e^H) = exp(Delta H).**
- Branching rule: d_{p'1} s_nu(p')|_{q^rho} = eps(nu) s_nu(q^rho), with eps(nu) = (q^{1/2} - q^{-1/2}) Sum_box q^{-c(box)}.
- Explicit form: Eq = -E_0(-x) - 1/vs(x) = (1/vs_1) [w^0](exp(-Sum vs_k w^k p_k/k) exp(-Sum vs_k w^{-k} d_{p_k}) - 1), with vs_k = q^{k/2} - q^{-k/2} (verified).
- Hence e^{-H} Eq e^H = (1/vs_1) [w^0](exp(-Sum vs_k w^k p_k/k + H(p - vs w^{-.}) - H(p)) - 1), a q-Miwa-shift two-point function.
- Why it is non-linear: x^{n+1} e^{-H} F_n e^H ~ x^{-1} for every completed cycle F_n, so all of them contribute already at genus 0.
- Why Route B's method fails on it. At genus 0, in the empty profile, the equation is sigma d_y [w^0] e^{-Y(w)} = e^{tau T(y)/(1+tau)}, with Y = w y/sigma + Sum_k k Pi_k w^{-k} G^(0)_(k)(y) (verified).
  - At u -> oo the right side has the non-integer exponential sector e^{tau u/(1+tau)}. In the Hurwitz case the corresponding factor is e^u, which is absorbed by z = T e^{-T}.
  - The left side produces e^{tau u/(1+tau)} only after resumming the w-residue (a saddle point). So no termwise leading-order analysis exists.
- The Hirota form Z d_1(Eq Z) - d_1 Z Eq Z = Z(e^beta o p) Z(e^{-beta} o p) is bilinear, but Eq has infinite order, so nothing is gained.
- The cure is the members of the same hierarchy that do not involve p'. There the shift enters only as A = (e^{beta E} - 1) F, which is polynomial in F, and E raises the order.

**(3) Cut-and-join beyond leading order.**
- At the next order in (1+s), sigma d_sigma Xi enters together with the new unknowns Q'_h(-1;R) = <Omega K kappa_1 lam_R>. Q has a simple zero at s = -1 (Route A s.3.8), so the flow does not close.
- It is not needed: the useful flow is the discrete Toda shift e^{tau x E}, not the framing derivative, which is subleading at the dilaton point.

**Comparison with Theorem D.**
- KP at fixed n gives excess-homogeneous merge relations. At the induction step these relations, Theorem B, and T_D at odd D are all satisfied by Cor. E's one-parameter family for every value of X = Xi_G(D+1) (checked algebraically). This is why round 1 needed global rigidity (TLR).
- The mKP equations are excess-inhomogeneous. Theorem M relates level D to level D - 1 with the diagonal coefficient -2 Pi_2 sigma^G tau delta, which is nonzero.
- After Cor. E, at R = (D) this coefficient becomes -2 sigma^G tau (D+1) delta Pi_{D+1}/Pi_D. It has no definite parity under tau <-> -1-tau, because the shift e^{tau x} breaks that symmetry. So the parity obstruction of T_D does not apply.

## 9. Files (HVC\lemmaA\round2\R3_toda)

**Libraries.**
- XiData.pm: exact Xi_g(T) from Route E's tables (lazy loading) or the Route E point engine.
- Copies of other routes' code: Hodge.pm, HodgePt.pm and WK.pm (Route E), mvlib.pm (Route A), Rec.pm (Route B).
- mkp_eqs.pm: generator of the members P_k in F/A/S form.

**Checks of Theorem M and its variants.**
- mkp_lead.pl: Theorem M; set PERTURB for perturbation tests.
- mkp_lead_qdata.pl: Theorem M on the q-tables.
- mkp_general_lead.pl: Theorem M_P for any k.
- mkp_hurwitz_limit.pl: the Hurwitz limit.
- g7_mkp_modp.pl and g8_mkp_modp.pl: mod-p tests at genus 7-9; they use A_mv\referee\rmv.pm in place, read-only.
- g7_points.pl: exact engine evaluation of genus-7 points (did not finish).

**Direct Schur-side checks of the equations.**
- mv_mkp_direct.pl: P_2.
- mv_hirota_direct.pl: P_k for k = 2, 3, 4.

**Other checks.**
- thmB_check.pl: Theorem B.
- eq_bosonization_check.pl and toda_eq_genus0_check.pl: line (2).
- R3_thmM_check.wls: Mathematica symbolic check, not run.

**Outputs:** out\*.out, one file per check above.

**Usage examples:**
- perl mkp_lead.pl 3 3 3 10
- perl mkp_lead.pl 5 5 1 14 1,2,...,13
- perl mkp_general_lead.pl 3 0 2 2 7
- perl mv_mkp_direct.pl 2/7 5 30
- perl g8_mkp_modp.pl 1000000007 41 9 17
