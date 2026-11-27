# Lemma A round 2: route R2 (saved by the main session from the workflow journal, 2026-10-31)

Subagents cannot write report files, so this file holds the agent's structured return value: the fields first, then the full report text (report_text) verbatim.

## all_genus_claim

1

## all_genus_statement

Proved modulo literature, uniformly in the genus with no induction on the genus or the level:
- For every h >= 2 and every multiset R of positive integers with D(R) = Sum(j_i - 1) <= 2h-3 (in particular every R with parts >= 2), Xi_h(R) = 0 identically in tau. This is Lemma A = M_irr(h) for all h, and hence (OL) for all G.
- The one-leaf formula (6.1) holds for all h, j and tau.
- The gap formula sigma^h Int_{Mbar_h} K Omega = -(2h-3)! [lambda^{2h}] S(e1 lambda)S(e2 lambda) holds for all h >= 2.
The literature inputs are ABDKS arXiv:2312.16950 Prop. 3.1 with Rem. 3.4 (x-y swap with trivial dual; LogTR = TR for affine y) and Eynard's kappa-class formula (CNTP 8 (2014)) as stated in ABDKS eqs. (49)-(50). Their hypotheses are checked and their normalisations are confirmed by exact data. TR is used only for the Hodge CohFT, on the curve C_tau, as allowed by brief s.0.

## claims

### [1]

#### depends_on

ABDKS arXiv:2312.16950v3 Prop. 3.1 (eq. 17) with Remarks 2.6, 3.2, 3.4 (log TR x-y swap with trivial dual; LogTR = TR for affine y)

#### evidence

Proof in report s.3: hypotheses of ABDKS Prop. 3.1 checked for C_tau, integration by parts, operator identity, order of residues, signs. Exact checks with 0 mismatches: 396 values vs the exact EO recursion at tau = 1 (g <= 7, n <= 5); 15 first-layer values vs the main session's Mathematica EO at tau = 2, 3; through P2, 1072 values vs Route E's TR-free engine (g <= 7, n <= 5, tau <= 4) and 216 values vs the lemmaA/data tables.

#### statement

Theorem 1 (closed n-point formula on the Hodge curve C_tau): for integer tau >= 1 and 2g-2+n > 0, [x^mu]F_{g,n} = ((-1)^{n-1}/prod mu_i) [hb^{2g-2+n}] [v_1^-1]...[v_n^-1] ( prod_i Q_{mu_i}(v_i) Sum_{n-cycles sigma} prod_i 1/(v_i - v_sigma(i) + (mu_i+mu_sigma(i)) hb/2) ), with Q_mu(v) = prod_{m in M(e2 mu)} (v+1+m hb/e2) / prod_{m in M(tau mu)} (v + m hb/tau), M(N) = {k-(N+1)/2 : k=1..N}, Laurent coefficients at v_i = oo extracted in the order v_1,...,v_n.

#### status

proved-modulo-literature

### [2]

#### depends_on

Theorem 1

#### evidence

Homogeneity proof (report s.4). Replacing the shift 1 by c in Q_mu, the iterated residue R_mu(c,hb) is a homogeneous polynomial of degree |mu| in (c, hb), so its hb-degree is at most |mu|. Numerically: 0 homogeneity violations in every xy_npoint.pl run, and every Lemma-A range point vanishes.

#### statement

Theorem 2 (degree bound = Lemma D of proof_status_mixed_scale.md for every integer tau >= 1): the TR correlators of C_tau satisfy [x^mu]F_{g,n} = 0 whenever |mu| <= 2g-3+n, for all g, n with 2g-2+n > 0.

#### status

proved-modulo-literature

### [3]

#### depends_on

Eynard, Commun. Number Theory Phys. 8 (2014) 541-588, as stated in ABDKS arXiv:2312.16950 eqs. (49)-(50); Kaufmann-Manin-Zagier dilaton pushforward

#### evidence

Proof in report s.5. The reparametrisation v = tau + kappa z, kappa^2 = -sigma, gives (p,q,r) = (1/tau, -1/e2, 1); the moment lemma gives q_k = (-1)^k kappa and s_i = (-1)^{i+1}/i; the degree involution gives K and sigma^{1-g}. Exact checks: Theorem 1 + P2 vs engine (1072 values); (6.1) vs exact EO at tau = 1 (h <= 7) and at tau = -1/2 (h <= 8).

#### statement

Theorem 3 (dictionary P2 = Lemma B, all g, n, tau with sigma != 0): omega_{g,n} = sigma^{1-g} Sum_k Int_{Mbar_{g,n}} K Omega prod_i psi_i^{k_i} d xi_{k_i}, with K = exp(-Sum kappa_j/j) and xi_0 = 1/(v-tau) = Sum Pi_mu x^mu. Hence Xi_g(mu) = Int K Omega prod 1/(1-mu_i psi_i) = sigma^{g-1} [x^mu]F_{g,n}/prod Pi_{mu_i}.

#### status

proved-modulo-literature

### [4]

#### depends_on

Lagrange inversion: xi_0 = Sum_mu Pi_mu x^mu

#### evidence

xi_k has poles only at v = tau and xi_k = 1/v + O(v^-2) at v = oo (report s.5(c)). moments_check.pl confirms 17 values at tau = 2, 5, 1/3.

#### statement

Moment lemma: Res_{v=tau} xi_k dv = 1 for all k >= 0 and all tau, where xi_k = (x d/dx)^k (1/(v-tau)) on C_tau.

#### status

proved

### [5]

#### depends_on

Theorems 2 and 3 (ABDKS Prop. 3.1; Eynard 2014)

#### evidence

Theorem 2 and Theorem 3 at every integer tau >= 1, then polynomiality of sigma^h Xi_h(R) in tau (report s.6). Consistent with Route E's complete M_irr(5), M_irr(6) and partial M_irr(7), and with every range point checked in report s.10.

#### statement

Corollary 4 (Lemma A = M_irr(h) for all h): for every h >= 2 and every multiset R of positive integers with D(R) <= 2h-3, Xi_h(R) = 0 identically in tau. Consequently (OL) holds for every G, with no inductive hypothesis needed.

#### status

proved-modulo-literature

### [6]

#### depends_on

elementary

#### evidence

Divided differences on a centred grid plus the generating function of centred factorial products (report s.7). xy_onepoint.pl: 120 exact values at tau = 1, 2, 3; the generic-tau version gives 42 more values.

#### statement

Theorem 5: for integers N >= 1 and M >= 0, [v^-1]_{v=oo} prod_{m in M(M)} (v+1+m a) / prod_{m in M(N)} (v + m b) = (M!/(N-1)!) [phi^{M-N+1}] e^phi S(a phi)^{M+1} S(b phi)^{1-N}, where S(y) = (y/2)/sinh(y/2).

#### status

proved

### [7]

#### depends_on

Theorems 1, 3, 5

#### evidence

Theorem 1 (n = 1) + Theorem 3 + Theorem 5 (report s.7). Exact checks with 0 mismatches: engine h <= 7 (180 + 84 + 154 values), exact EO at tau = 1 (h <= 7, 90 values) and at tau = -1/2 (h <= 8, 287 values).

#### statement

Corollary 6: the one-leaf formula (6.1), sigma^h Xi_h(j) = prod_{l=1}^{2h-2}(j-l) [lambda^{2h}] S(e1 lambda)^{1+j e2} S(e2 lambda)^{1+j e1}, holds for all h, j, tau. Equivalently Sum_h sigma^h Xi_h(j) lambda^{2h} = (j-1)! [phi^{j+1}] e^phi S(e1 lambda phi)^{1+j e2} S(e2 lambda phi)^{1+j e1}. In particular Xi_h(j) = 0 for 1 <= j <= 2h-2, and Route E's Zhat_j(u) = Res_{l=0} l^-3 (1+u/l)^{j-1} S1 S2 B^j dl holds for all h.

#### status

proved-modulo-literature

### [8]

#### depends_on

Corollary 6

#### evidence

(6.1) at j = 0 plus the string and dilaton equations (report s.7). gap_check.pl: h = 2..7 at 4 values of tau, 24 values, 0 mismatches.

#### statement

Corollary 7 (gap, r = 0): for h >= 2, sigma^h Int_{Mbar_h} K Omega = -(2h-3)! [lambda^{2h}] S(e1 lambda) S(e2 lambda). Also Xi_h(1 cup R') = |R'| Xi_h(R') and Xi_h(0) = (2-2h) Int_{Mbar_h} K Omega.

#### status

proved-modulo-literature

### [9]

#### depends_on

Theorem 1

#### evidence

Partial fractions in Theorem 1 at n = 2. xy_twopoint_diff.pl: 36 values vs the engine, 0 mismatches.

#### statement

Two-point difference form: [x^{m1}x^{m2}]F_{g,2} = (1/(m1 m2)) [hb^{2g}] (1/(2k hb)) [v^-1] Q_{m2}(v) (Q_{m1,+}(v + k hb) - Q_{m1,+}(v - k hb)), where k = (m1+m2)/2 and Q_+ is the polynomial part at v = oo.

#### status

proved-modulo-literature

### [10]

#### depends_on

character theory of the symmetric group

#### evidence

Murnaghan-Nakayama rule

#### statement

F1: for the MV series R = Sum_lambda W_lambda s_lambda(p), [p_j]R at p = 0 equals (1/j) Sum over hooks lambda of size j of (-1)^{leg(lambda)} W_lambda. With a p_1 background the hooks become j-border strips.

#### status

proved

### [11]

#### depends_on

Route E engine (copied to R2_onepoint/hodge_engine/); referee EO codes (copied to R2_onepoint/trcheck/)

#### evidence

report s.10, checks V1-V20, 0 mismatches; raw outputs in R2_onepoint/runs_all_checks.txt

#### statement

Computer verification of the chain beyond the previously known range: genus-7 one-leaf values (new engine table; tau = 1..4 via Theorem 1, 7 values of tau via (6.1)), genus-8 one-leaf values at tau = -1/2 (exact EO), and multi-leaf values up to n = 5. About 2400 exact comparisons in total.

#### status

computer-verified


## key_formulas

- Q_mu(v) = prod_{m in M(e2 mu)} (v + 1 + m hb/e2) / prod_{m in M(tau mu)} (v + m hb/tau), M(N) = {k - (N+1)/2 : k = 1..N}, e1 = -tau, e2 = 1+tau, sigma = tau(1+tau)
- Theorem 1: [x^mu]F_{g,n} = ((-1)^{n-1}/prod mu_i) [hb^{2g-2+n}] [v_1^-1]...[v_n^-1] ( prod_i Q_{mu_i}(v_i) Sum_{n-cycles sigma} prod_i 1/(v_i - v_sigma(i) + (mu_i + mu_sigma(i)) hb/2) ); for n = 1: [x^j]F_{g,1} = j^-2 [hb^{2g}] [v^-1] Q_j(v)
- xhat = tau S(hb d_v/tau) log v - e2 S(hb d_v/e2) log(1+v), and -mu Sc(mu hb d) xhat = log Q_mu(v), where S(y) = (y/2)/sinh(y/2) and Sc = 1/S
- Theorem 2 (homogeneity): R_mu(c, hb) is homogeneous of degree |mu|, so [x^mu]F_{g,n} = 0 for |mu| <= 2g-3+n (Lemma D, every integer tau)
- Theorem 3 (P2): omega_{g,n} = sigma^{1-g} Sum_k Int_{Mbar_{g,n}} K Omega prod psi_i^{k_i} d xi_{k_i}, K = exp(-Sum kappa_j/j); Xi_g(mu) = sigma^{g-1} [x^mu]F_{g,n} / prod Pi_{mu_i}
- v = tau + kappa z, kappa^2 = -sigma: (p,q,r) = (1/tau, -1/e2, 1), q_k = (-1)^k kappa, s_i = (-1)^{i+1}/i; Res_{v=tau} xi_k dv = 1
- Theorem 5: [v^-1] prod_{m in M(M)} (v+1+m a) / prod_{m in M(N)} (v+m b) = (M!/(N-1)!) [phi^{M-N+1}] e^phi S(a phi)^{M+1} S(b phi)^{1-N}
- (6.1) for all h: Sum_h sigma^h Xi_h(j) lambda^{2h} = (j-1)! [phi^{j+1}] e^phi S(e1 lambda phi)^{1+j e2} S(e2 lambda phi)^{1+j e1}, i.e. sigma^h Xi_h(j) = prod_{l=1}^{2h-2}(j-l) [lambda^{2h}] S(e1 lambda)^{1+je2} S(e2 lambda)^{1+je1}
- Gap: sigma^h Int_{Mbar_h} K Omega = -(2h-3)! [lambda^{2h}] S(e1 lambda) S(e2 lambda); Xi_h(0) = (2-2h) Int_{Mbar_h} K Omega; Xi_h(1 cup R') = |R'| Xi_h(R')
- Zhat_j(u) = Sum_h sigma^h Xi_h(j) u^{2h-2}/(2h-2)! = Res_{l=0} l^-3 (1+u/l)^{j-1} S(e1 l) S(e2 l) B(l)^j dl, B = S(e1 l)^{e2} S(e2 l)^{e1}
- Two-point: [x^{m1}x^{m2}]F_{g,2} = (1/(m1 m2)) [hb^{2g}] (1/(2k hb)) [v^-1] Q_{m2}(v) (Q_{m1,+}(v+k hb) - Q_{m1,+}(v-k hb)), k = (m1+m2)/2

## next_steps

- Referee/critic pass on report s.3-5. Points to check: the application of ABDKS Prop. 3.1 to C_tau; the integration-by-parts step and the order of residues; Eynard's formula (50) against the original CNTP paper; the gamma = 0 parametrisation; the degree-involution bookkeeping.
- Main session: record Lemma A (all h), Lemma B (P2) and Lemma D as proved modulo literature in proof_status_mixed_scale.md and leaf_structure_proof.md, and feed Lemma A into Corollary R (all-genus regularity).
- TR-free proof of Theorem 1: derive the dilaton-point n-point formula directly from the Okounkov-Pandharipande E-operators or Orlov-Scherbin closed formulas (both allowed inputs), using that the dual side is trivial. The n = 1 case is a one-part MV-type evaluation.
- Extract Route E's multi-leaf structures (SN) and (NR) from Theorem 1 (divided differences with shifted nodes), and the boundary layer D = 2h-2 from the c = 0 products.
- Optional: finish the genus-7 one-leaf framings tau = 5..8 (perl xy_npoint.pl TAU 7 1 21 engine_h7_r1.txt in R2_onepoint/).

## obstructions

- Literature dependence. Theorem 1 rests on ABDKS arXiv:2312.16950 Prop. 3.1 with Rem. 3.4 (log-TR x-y swap with trivial dual). Theorem 3 rests on Eynard's kappa-class formula (CNTP 8 (2014)) in the form of ABDKS eqs. (49)-(50). Both are applied with their hypotheses checked, and their normalisations are confirmed by exact data (genus <= 7, and genus 8 for n = 1), but they are not re-proved here.
- TR is used, but only for the Hodge CohFT on its curve C_tau, which brief s.0 allows (user clarification of 2026-10-31). The computed task text's older blanket TR ban is superseded by the brief. If a TR-free proof is wanted, Theorems 1 and 3 need independent (fermionic or localization) proofs.
- The fermionic MV route as briefed did not close by itself. On the MV side alone I found no all-order control of the semiclassical expansion of the one-leaf bilocal at the dilaton point u = oo, on the other sheet (F2).
- The genus-7 one-leaf checks at framings tau = 5..8 were still running at submission. tau = 1..4 are done, and (6.1) was checked against the genus-7 engine at 7 values of tau.
- Not obtained: closed multi-leaf forms (Route E's SN and NR) and the extra zero at j = 2h in the Hurwitz limit.

## requested_checks

- Referee: check the hypotheses and conventions in applying ABDKS arXiv:2312.16950 Prop. 3.1 + Rem. 3.4 to C_tau: x_A = log x, y_A = v; LogTR-vital points v = 0, -1; v = oo not vital; xhat = tau S(hb d/tau) log v - e2 S(hb d/e2) log(1+v), with ABDKS's calligraphic S = 1/S. Also check the integration-by-parts step and the order of residues (report s.3).
- Referee: compare ABDKS eqs. (49)-(50) with Eynard, CNTP 8 (2014) 541-588 (statement and normalisation of the kappa-class formula), and re-check report s.5 (a)-(f): the gamma = 0 parametrisation, the q_k, and the degree involution.
- Critic: confirm that the x-y swap on C_tau and Eynard's formula for C_tau fall under brief s.0's ALLOWED list (TR considerations for the Hodge CohFT only), and that nothing in the argument touches the resolved conifold x C^2 theory.
- Optional (main session, Mathematica): extend the hvc_tr_tau.wls runs at tau = 2, 3 to genus 4-5, n <= 3, and compare full coefficient tables with perl xy_npoint.pl TAU G N JMAX. This checks Theorem 1 at tau != 1 independently of P2.

## route

R2 (round 2): one-leaf function -> x-y swap with trivial dual on the Hodge curve C_tau (log TR of ABDKS) + Eynard's kappa-class dictionary; Lemma A and (6.1) for all h

## summary

Lemma A (M_irr(h)) is proved for all h, together with the one-leaf formula (6.1) for all h and the r = 0 gap formula, modulo two literature results about TR for the Hodge CohFT, which brief s.0 allows. The key observation: on the Hodge spectral curve C_tau (x = v^tau/(1+v)^(tau+1), y = -v), X = log x = -e1 log y - e2 log(y-1) + const has no critical points as a function of y, so the x-y dual curve is unramified. The log-TR x-y swap formula of Alexandrov-Bychkov-Dunin-Barkowski-Kazarian-Shadrin (arXiv:2312.16950, Prop. 3.1 with Rem. 3.4) then gives every omega_{g,n} in closed form. Integration by parts at x = 0 turns this into a finite residue of centred products Q_mu(v) (Theorem 1: all g, n, integer tau). A homogeneity argument gives [x^mu]F_{g,n} = 0 for |mu| <= 2g-3+n (Theorem 2 = Lemma D for every integer tau). Eynard's kappa-class formula (CNTP 8 (2014), as stated in ABDKS (49)-(50)) gives the dictionary P2 for all g, n, tau (Theorem 3 = Lemma B). Together: Xi_h(R) = 0 identically in tau for every R with parts >= 1 and D(R) <= 2h-3, every h, with no induction on the genus or the level. For one leaf the residue is a divided difference that evaluates in closed form (Theorem 5) to exactly (6.1), for all h; at j = 0 it gives the gap formula sigma^h Int_{Mbar_h} K Omega = -(2h-3)! [l^{2h}] S(e1 l)S(e2 l). All of this was checked exactly: about 2400 comparisons with 0 mismatches. They are against Route E's TR-free Hodge engine (1072 values, genus <= 7, n <= 5, tau <= 4; the genus-7 one-leaf table was computed here), the exact EO recursion at tau = 1 (396 values, g <= 7, n <= 5) and at tau = -1/2 ((6.1), up to genus 8), the main session's Mathematica EO recursion at tau = 2, 3, and the lemmaA/data tables.

Where it stops: the all-genus statements are proved modulo ABDKS Prop. 3.1 and Eynard's formula. Their hypotheses are checked and their normalisations confirmed by the data, but they are not re-proved. The fermionic Marino-Vafa route as briefed did not close by itself: no all-order control of the semiclassics at the dilaton point on the MV side alone. A TR-free (E-operator) proof of Theorem 1 remains to be written. TR is used only on C_tau (the Hodge CohFT), never for the resolved conifold x C^2 theory. The computed task text's older blanket TR ban is superseded by brief s.0. The genus-7 one-leaf runs at framings tau = 5..8 were still running at submission.

## report_text

# Route R2 (round 2): the one-leaf function, formula (6.1) and Lemma A in all genera, via the x-y swap on the Hodge curve C_tau

Folder: HVC\lemmaA\round2\R2_onepoint\ (HVC = G:\My Drive\Projects\Old\refGW\hodge_vertex_conifold). 2026-10-31.

## 0. Results

The route was briefed as "the one-leaf function from the fermionic Marino-Vafa formula". The fermionic computation by itself did not close (section 9). What closed it is the x-y swap for the Hodge spectral curve C_tau. After the swap the dual curve is unramified, and the log-TR x-y swap theorem of Alexandrov-Bychkov-Dunin-Barkowski-Kazarian-Shadrin (ABDKS) gives every correlator in closed form.

| | Statement | Status | Depends on |
|---|---|---|---|
| Thm 1 | closed n-point formula for the TR correlators of C_tau at x = 0, all g, n, integer tau >= 1 (s.3) | proved-modulo-literature | ABDKS arXiv:2312.16950v3, Prop. 3.1 (eq. (17)), Rem. 3.4 |
| Thm 2 | degree bound: [x^mu]F_{g,n} = 0 for abs(mu) <= 2g-3+n, every integer tau >= 1 (= Lemma D of proof_status_mixed_scale.md) (s.4) | proved-modulo-literature | Thm 1; homogeneity argument proved here |
| Thm 3 | dictionary P2 (= Lemma B) for all g, n, tau (s.5) | proved-modulo-literature | Eynard, CNTP 8 (2014), in the form of ABDKS (49)-(50); moment lemma proved here |
| Cor 4 | **Lemma A = M_irr(h) for every h** (indeed Xi_h(R) = 0 for every R with parts >= 1 and D(R) <= 2h-3); hence (OL) for every G (s.6) | proved-modulo-literature | Thm 2, Thm 3 |
| Thm 5 | closed evaluation of the one-leaf residue (divided differences) (s.7) | proved | elementary |
| Cor 6 | **the one-leaf formula (6.1) for all h, j, tau**; Route E's generating function Zhat_j(u) (s.7) | proved-modulo-literature | Thm 1 (n = 1), Thm 3, Thm 5 |
| Cor 7 | gap formula sigma^h Int_{Mbar_h} K Omega = -(2h-3)! [l^{2h}] S(e1 l) S(e2 l), all h >= 2 (s.7) | proved-modulo-literature | Cor 6, string and dilaton equations |
| F1-F3 | fermionic MV route (s.9) | F1 proved; F2 is the obstruction; F3 a remark | |

No step uses induction on the genus or on the level. No inductive hypothesis of (OL) is used. Thm 1 and Thm 3 are identities of all-genus generating series; Thm 2 is a scaling statement (s.8).

Verification (s.10): exact Perl/BigRat comparisons, about 2400 in all, with no mismatch anywhere.
- 1072 values against Route E's TR-free Hodge engine (genus <= 7, n <= 5, tau <= 4). The genus-7 one-leaf table was computed here.
- 396 values against the exact EO recursion at tau = 1 (g <= 7, n <= 5).
- 287 values of (6.1) against the exact EO recursion at tau = -1/2, up to genus 8.
- 15 values against the main session's independent Mathematica EO recursion at tau = 2, 3.
- 216 values against the converted lemmaA/data tables.

**Use of TR (brief s.0).** TR enters only for the Hodge CohFT Omega = Lambda^vee(1)Lambda^vee(tau)Lambda^vee(-1-tau)/(-sigma)^h, i.e. the framed C^3 vertex. It is used in three ways, all for that CohFT:
- TR on its curve C_tau;
- the x-y swap for that curve;
- Eynard's intersection-theoretic formula for that curve.

This is exactly what the user's clarification of 2026-10-31 allows: "it's ok to use TR considerations applied to the Hodge CohFT". The computed task text still contained the older blanket "no TR". Brief s.0, which the task text itself declares overriding, supersedes it. Nothing in this report uses TR, or anything else, for the resolved conifold x C^2 theory, its free energies, its polar part or its graph sum. Input (H) is not used.

## 1. Conventions

- c = 1; e1 = -tau, e2 = 1 + tau, sigma = tau(1 + tau) = -e1 e2, e1 + e2 = 1.
- S(y) = (y/2)/sinh(y/2) (the brief's S, even in y); Sc(y) = 1/S(y) = 2 sinh(y/2)/y. ABDKS write calligraphic S for Sc.
- For an integer N >= 0: M(N) = {k - (N+1)/2 : k = 1, ..., N}, a centred set symmetric about 0 with spacing 1.
- Pi_mu = C(mu(tau+1) - 1, mu - 1) = prod_{a=1}^{mu-1} (mu tau + a)/(mu-1)!.
- Xi_h(R) is as in brief s.1. Pushing forward the dilaton leaves (Kaufmann-Manin-Zagier; proof_M_hurwitz.md s.1a) gives Xi_h(R) = Int_{Mbar_{h,r}} K Omega prod_i 1/(1 - j_i psi_i), with K = exp(-Sum_{k>=1} kappa_k/k).
- The curve C_tau is P^1 with coordinate v, and:
  - x = v^tau/(1+v)^(tau+1), X = log x = tau log v - e2 log(1+v);
  - y = -v, omega_{0,1} = y dX, B = dv1 dv2/(v1 - v2)^2;
  - TR is taken with respect to X. Then dX = (tau - v) dv/(v(1+v)) has one simple zero, v = tau, which is the only branch point;
  - omega_{g,n} = d_1 ... d_n F_{g,n}, with F_{g,n} expanded in the x_i at x_i = 0, i.e. at v_i = oo (x ~ 1/v there);
  - these are the conventions of proof_M_hurwitz.md s.2 and of review_checks/status_tr*.pl.
- Laplace leaves: xi_0 = 1/(v - tau), xi_k = D^k xi_0 with D = d/dX = x d/dx.
  - Lagrange inversion with w = 1/v, x = w/(1+w)^(tau+1) gives Sum_{mu>=0} C((tau+1)mu, mu) x^mu = (1+w)/(1 - tau w) = 1 + (tau+1) xi_0.
  - Hence xi_0 = Sum_{mu>=1} Pi_mu x^mu and xi_k = Sum_{mu>=1} mu^k Pi_mu x^mu.

## 2. The idea

Write X as a function of y = -v. Then X = -e1 log y - e2 log(y - 1) + const, so X has no critical points in y at all.

- After the x-y swap (x^vee = y, y^vee = X), the dual curve is unramified.
- Its only features are two logarithmic points: y = 0 with weight -e1 = tau, and y = 1 with weight -e2.
- For such curves, Prop. 3.1 of ABDKS (log TR via the x-y swap) gives omega_{g,n} in closed form, for all g and n at once.
- The entire all-genus content sits in an hbar-deformation xhat of X. In xhat each log point is dressed by S(hbar d/weight).
- Evaluated at x = 0 (v = oo), the closed formula becomes a finite residue of products of shifted linear factors.
  - The vanishing range then follows from scaling (Thm 2).
  - The one-leaf case evaluates to (6.1) (Thm 5).
- The two factors S(e1 lambda), S(e2 lambda) of (6.1) are the two log points.

In the language of Route A (Cor. 2.2), v = oo is the dilaton point (u = oo, on the other sheet). The x-y swap is what makes all orders explicit there. This is precisely what the fermionic computation lacked (section 9).

## 3. Theorem 1 (closed n-point formula)

**Theorem 1.** Let tau >= 1 be an integer, g >= 0, n >= 1 with 2g - 2 + n > 0, and mu = (mu_1, ..., mu_n) positive integers. Put

    Q_mu(v) = prod_{m in M(e2 mu)} (v + 1 + m hb/e2) / prod_{m in M(tau mu)} (v + m hb/tau)

(so Q_mu = x^{-mu} + O(hb^2)). Then

    [x_1^{mu_1} ... x_n^{mu_n}] F_{g,n} = ((-1)^{n-1}/(mu_1 ... mu_n)) [hb^{2g-2+n}] [v_1^{-1}] ... [v_n^{-1}] ( prod_i Q_{mu_i}(v_i) * Sum_{sigma n-cycle} prod_{i=1}^n 1/(v_i - v_{sigma(i)} + (mu_i + mu_{sigma(i)}) hb/2) ).

Here [v_i^{-1}] is the coefficient of v_i^{-1} in the Laurent expansion at v_i = oo. The coefficients are extracted in the order v_1, ..., v_n, each over the field of rational functions of the remaining variables; equivalently, in the region |v_1| >> ... >> |v_n|. For n = 1 the cycle sum is 1/(mu hb), so

    [x^j] F_{g,1} = j^{-2} [hb^{2g}] [v^{-1}] Q_j(v).

**Literature input: ABDKS, arXiv:2312.16950v3, Prop. 3.1, eq. (17), with Rem. 3.2 (eqs. (18), (20)) and Rem. 3.4.**

Setting.
- y = z is an affine coordinate on P^1.
- x(z) = Sum_i log(z - a_i)/alpha_i + R(z). Here the a_i are the LogTR-vital singular points of x: simple poles of dx, with residue 1/alpha_i, at which dy is regular.

Then the LogTR differentials with omega_1^(0) = -y dx are

    omega^(g)_n = [hb^{2g-2+n}] (-1)^n prod_i ( Sum_{r>=0} (d_i (1/dx_i))^r [u_i^r] exp(-u_i (Sc(u_i hb d_{y_i}) xhat_i - x_i)) ) (-1)^{n-1} Sum_{n-cycles} prod_i dz_i/(w_i - wbar_{sigma(i)}),

with xhat = Sum_i (1/alpha_i) S(alpha_i hb d_y) log(z - a_i) + R(z), w_i = z_i + u_i hb/2 and wbar_i = z_i - u_i hb/2.

- d(1/dx) maps a differential eta to d(eta/dx).
- The operators act on everything to their right.
- By Rem. 3.4, LogTR coincides with ordinary TR when y = z is affine, because dy then has no simple poles.

**Hypotheses on C_tau.** Take x_A := X and y_A := v, so that -y_A dx_A = -v dX = y dX = omega_{0,1}.
- dX = (tau - v) dv/(v(1+v)) has the single simple zero v = tau, and dv is regular there.
- dv has no zeros, so the dual TR is trivial.
- X = tau log v - e2 log(1+v), so a_1 = 0 with 1/alpha_1 = tau, a_2 = -1 with 1/alpha_2 = -e2, and R = 0.
- Both of these points are LogTR-vital: dX has simple poles there and dv is regular.
- v = oo is a simple pole of dX (residue 1), but dv has a double pole there. So it is not LogTR-vital and does not contribute (ABDKS Rem. 2.6).

Therefore xhat = tau S(hb d_v/tau) log v - e2 S(hb d_v/e2) log(1+v). (S is even, so the sign of alpha_2 does not matter.)

**Proof of Theorem 1.**

(i) Coefficient extraction.
- For 2g-2+n > 0, omega_{g,n} has poles only at v_i = tau and none on the diagonals. Hence prod_i mu_i [x^mu]F_{g,n} = Res_{v_1=oo} ... Res_{v_n=oo} prod_i x_i^{-mu_i} omega_{g,n}, taken in the stated order. Each residue is over the field of rational functions of the later variables, and the result equals the joint Taylor coefficient.
- At each order in hb, the right side of (17) is a finite sum of rational differentials. Indeed, the exponent is O(hb^2), with hb^{2k} accompanied by at most u^{2k+1}, and the cycle factor is polynomial in u at each order. So the residues may be taken term by term.

(ii) Integration by parts.
- Let eta be a rational differential in v_1 whose coefficients are rational in the other variables. Then Res_{v_1=oo} x_1^{-mu_1} d_1(eta/dX_1) = -Res d_1(x_1^{-mu_1}) eta/dX_1 = mu_1 Res x_1^{-mu_1} eta. This uses d(x^{-mu}) = -mu x^{-mu} dX and the fact that an exact differential has zero residue.
- The operators in the other variables commute with Res_{v_1}.
- Hence, for every i, Sum_r (d_i (1/dX_i))^r [u_i^r] becomes evaluation at u_i = mu_i.

(iii) Operator identity.
- For an integer N >= 1, Sum_{m in M(N)} e^{m t} = sinh(N t/2)/sinh(t/2). Hence:
  - mu tau Sc(mu hb d) S(hb d/tau) = Sum_{m in M(tau mu)} e^{m hb d/tau};
  - mu e2 Sc(mu hb d) S(hb d/e2) = Sum_{m in M(e2 mu)} e^{m hb d/e2}.
- Combining these with Taylor's formula log(v + c) = e^{c d} log v (formal in hb) gives
  -mu Sc(mu hb d) xhat = Sum_{m in M(e2 mu)} log(1 + v + m hb/e2) - Sum_{m in M(tau mu)} log(v + m hb/tau).
- Therefore x^{-mu} exp(-mu (Sc(mu hb d) xhat - X)) = exp(-mu Sc(mu hb d) xhat) = Q_mu(v), as formal series in hb with rational coefficients.

(iv) Signs.
- y_A = z gives sqrt(dw dwbar) = dz. At u = mu, w_i - wbar_{sigma(i)} = v_i - v_{sigma(i)} + (mu_i + mu_{sigma(i)}) hb/2.
- The prefactors of (17) give (-1)^n (-1)^{n-1} = -1.
- Each Res_{v=oo} f dv = -[v^{-1}] f contributes a further -1.
- In total the sign is (-1)^{n-1}. QED

**Two-point difference form.** For n = 2, put k = (m1+m2)/2 and let Q_+ be the polynomial part at v = oo. Partial fractions give

    [x^{m1} x^{m2}] F_{g,2} = (1/(m1 m2)) [hb^{2g}] (1/(2 k hb)) [v^{-1}] Q_{m2}(v) (Q_{m1,+}(v + k hb) - Q_{m1,+}(v - k hb)).

## 4. Theorem 2 (degree bound; Lemma D for every integer tau)

**Theorem 2.** For every integer tau >= 1 and all g, n with 2g-2+n > 0, [x^mu] F_{g,n} = 0 whenever |mu| <= 2g - 3 + n.

*Proof.*
- Replace the 1 in the numerator factors by a parameter c: Q_mu(v; c) = prod_{m in M(e2 mu)} (v + c + m hb/e2) / prod_{m in M(tau mu)} (v + m hb/tau).
- Let R_mu(c, hb) be the iterated coefficient of Theorem 1, before taking [hb^{2g-2+n}], built from the Q_{mu_i}(v_i; c).
- Give v_i, c and hb degree 1.
  - Q_mu(v; c) is homogeneous of degree e2 mu - tau mu = mu.
  - Each cycle factor 1/(v_i - v_j + k hb) is homogeneous of degree -1.
  - So the integrand is homogeneous of degree |mu| - n. Taking the coefficient of prod_i v_i^{-1} (degree -n) leaves something homogeneous of degree |mu| in (c, hb).
- R_mu(c, hb) is a polynomial:
  - At v = oo, Q_mu(v; c) = v^mu Sum_a q_a(c, hb) v^{-a}, with q_a a homogeneous polynomial of degree a.
  - The cycle factors expand as Sum_l (v_j - k hb)^l v_i^{-l-1} (or with the roles of i and j exchanged), which is polynomial in the later variables and in hb.
  - Every extraction [v_i^{-1}] is a finite sum, because the exponent of v_i is bounded above at each step.
- Hence R_mu(c, hb) = Sum over b + p = |mu|, b, p >= 0, of r_{b,p} c^b hb^p, and [hb^p] R_mu(1, hb) = r_{|mu|-p, p} = 0 for p > |mu|.
- Theorem 1 with p = 2g-2+n gives the claim. QED

*Remarks.*
- (a) This is Lemma D of proof_status_mixed_scale.md (s.1.7(c), s.2), now for every integer tau >= 1 and all g, n.
  - At tau = 1 it re-proves Theorem 3.1 of proof_M_hurwitz.md, which was proved there with the involution v -> 1/v.
  - The "missing analogue of the residue step" of proof_status s.2 is not needed.
- (b) The first layer that can be non-zero, |mu| = 2g-2+n, is the c^0 part. It is computed from the scale-free products prod (v + m hb/e2)/prod (v + m hb/tau).
- (c) The bound is sharp in general: Route E's control points with D = 2h-2 are non-zero.
- (d) The scaling is invisible in the Hodge integrals themselves. It lives in the shift 1 of v + 1, i.e. in the log point v = -1, and it is uniform in g.

## 5. Theorem 3 (dictionary P2 = Lemma B)

**Theorem 3.** For every tau with sigma != 0 and all g, n with 2g-2+n > 0,

    omega_{g,n} = sigma^{1-g} Sum_{k_1..k_n >= 0} ( Int_{Mbar_{g,n}} K Omega prod_i psi_i^{k_i} ) prod_i d xi_{k_i}(v_i).

Hence

    Xi_g(mu) = Int_{Mbar_{g,n}} K Omega prod_i 1/(1 - mu_i psi_i) = sigma^{g-1} [x^mu] F_{g,n} / prod_i Pi_{mu_i}.

**Literature input.** B. Eynard, "Invariants of spectral curves and intersection theory of moduli spaces of complex curves", Commun. Number Theory Phys. 8 (2014) 541-588, in the form stated by ABDKS, arXiv:2312.16950 s.3.4-3.5, eqs. (40)-(50).

Setting.
- dx = z dz/((1+alpha z)(1+beta z)(1+gamma z)) on P^1, with affine coordinate z and z = 0 the zero of dx.
- p = -(alpha-beta)(alpha-gamma), q = -(beta-gamma)(beta-alpha), r = -(gamma-alpha)(gamma-beta), so that 1/p + 1/q + 1/r = 0.
- xi^A_k = (-d_x)^k (-1/z).
- y is any power series at z = 0 with y'(0) != 0.
- Define q_k = -Res_{z=0} xi^A_{k+1} y dx, and define s_i by 1 + (q_1/q_0) t + (q_2/q_0) t^2 + ... = exp(-Sum_i s_i t^i).

Then the TR differentials with omega_{0,1} = -y dx are

    omega^(g)_n = q_0^{2-2g-n} Sum_k Int_{Mbar_{g,n}} exp(Sum_j s_j kappa_j) Lambda(p) Lambda(q) Lambda(r) prod_i psi_i^{k_i} d xi^A_{k_i}(z_i),

where Lambda(t) = 1 + t lambda_1 + ... + t^g lambda_g.

For the special y of (40) this is BKMP for C^3 (Bouchard-Marino; L. Chen; J. Zhou). DOSS (CMP 328 (2014)) identify TR with Givental's formula.

**Proof of Theorem 3.**

(a) Coordinates.
- Put v = tau + kappa z with kappa^2 = -sigma.
- Then tau - v = -kappa z, v = tau (1 + (kappa/tau) z) and 1 + v = e2 (1 + (kappa/e2) z). So dX = (-kappa^2/sigma) z dz/((1+alpha z)(1+beta z)) = z dz/((1+alpha z)(1+beta z)(1+gamma z)), with alpha = kappa/tau, beta = kappa/e2, gamma = 0.
- Hence p = -(alpha-beta) alpha = 1/tau, q = -beta(beta-alpha) = -1/e2, r = -alpha beta = 1, and X = (1/p) log(1+alpha z) + (1/q) log(1+beta z) + const.
- gamma = 0 only places the third pole of dX at z = oo; the three poles stay distinct. Both sides of Eynard's formula are continuous in gamma, or one may use a Moebius change of the affine coordinate.

(b) Leaves. xi^A_0 = -1/z = -kappa/(v - tau) = -kappa xi_0 and -d_x = -D, so xi^A_k = (-1)^k (-kappa) xi_k.

(c) Moment lemma. Res_{v=tau} xi_k dv = 1 for all k >= 0 and all tau.
- *Proof.* xi_k is rational with poles only at v = tau, because D = (v(1+v)/(tau - v)) d/dv preserves this property.
- xi_k vanishes at v = oo, since it equals Sum_mu mu^k Pi_mu x^mu = O(x).
- Hence Res_{v=tau} xi_k dv = -Res_{v=oo} xi_k dv = [v^{-1}]_{v=oo} xi_k = 1, because x = 1/v + O(v^{-2}) gives xi_k = Pi_1 x + O(x^2) = 1/v + O(v^{-2}). QED
- Consequently, with y = y_A = v: q_k = -Res xi^A_{k+1} v dX = Res v d xi^A_k = -Res xi^A_k dv = (-1)^k kappa.
- So q_0 = kappa and Sum (q_k/q_0) t^k = 1/(1+t), which gives s_i = (-1)^{i+1}/i.
- Direct check: q_0 = Res (tau + kappa z) z^{-2} dz = kappa, and q_1 = kappa alpha beta = -kappa.

(d) Hodge factor. Lambda(1/t) = (-t)^{-g} Lambda^vee(-t), where Lambda^vee(s) = Sum_i (-1)^i lambda_i s^{g-i}. So Lambda(1/tau) Lambda(-1/e2) Lambda(1) = sigma^{-g} Lambda^vee(-1) Lambda^vee(-tau) Lambda^vee(1+tau).

(e) Assembly.
- So far, omega_{g,n} = kappa^{2-2g-n} (-kappa)^n sigma^{-g} Sum_k Int C_k prod_i d xi_{k_i}, with C_k = exp(Sum_j (-1)^{j+1} kappa_j/j) Lambda^vee(-1) Lambda^vee(-tau) Lambda^vee(1+tau) prod_i (-psi_i)^{k_i}.
- Let iota be the degree involution (multiplication by (-1)^d on H^{2d}(Mbar_{g,n})). It is a ring automorphism, and Int alpha = (-1)^{3g-3+n} Int iota(alpha).
  - iota(kappa_j) = (-1)^j kappa_j turns the exponent into -Sum kappa_j/j, which is the class K.
  - iota(psi) = -psi.
  - iota(Lambda^vee(s)) = (-1)^g Lambda^vee(-s), so iota of the Hodge product is (-1)^{3g} Lambda^vee(1) Lambda^vee(tau) Lambda^vee(-1-tau) = (-1)^{3g} (-sigma)^g Omega.
- Therefore Int C_k = (-1)^{n+1} (-sigma)^g Int K Omega prod_i psi_i^{k_i}.
- Using kappa^{2-2g} = (-sigma)^{1-g}: omega_{g,n} = (-1)^n (-sigma)^{1-g} sigma^{-g} (-1)^{n+1} (-sigma)^g Sum_k (...) = sigma^{1-g} Sum_k Int K Omega prod_i psi_i^{k_i} d xi_{k_i}.

(f) Coefficients.
- F_{g,n} = sigma^{1-g} Sum_k Int K Omega prod psi^{k_i} xi_{k_i}(x_i), since all xi_k vanish at x = 0.
- [x^mu] xi_k = mu^k Pi_mu.
- Xi_g(mu) equals Int K Omega prod_i 1/(1 - mu_i psi_i) by the dilaton-leaf pushforward:
  - Sum_m (1/m!) pi_{m*} prod_k T(psi_{n+k}) = exp(Sum_j p_j kappa_j) with Sum_j p_j z^j = log(1 - z), i.e. p_j = -1/j;
  - the psi of the remaining points may be replaced by pull-backs, because every forgotten point carries psi^{>=2}. QED

This is P2 exactly as written in proof_M_hurwitz.md s.2 and brief s.5, now for all g, n, tau. It also explains the n-independent factor sigma^{1-g}, and it replaces the sketched route of referee_hurwitz_1.md s.4.

## 6. Corollary 4: Lemma A for all h

**Corollary 4.** For every h >= 2 and every multiset R of positive integers with D(R) <= 2h - 3, Xi_h(R) = 0 identically in tau. In particular:
- Lemma A (= M_irr(h), parts >= 2) holds for every h;
- (OL) holds for every G, with its conclusion Xi_G(D+1) = 0 true even without its hypotheses.

*Proof.*
- Fix an integer tau >= 1.
- Theorem 3 gives Xi_h(R) = sigma^{h-1} [x^R] F_{h,r} / prod Pi_{j_i}, with all Pi_j > 0.
- Theorem 2 gives [x^R] F_{h,r} = 0, because |R| = D(R) + r <= 2h - 3 + r.
- So Xi_h(R; tau) = 0 for tau = 1, 2, 3, ....
- sigma^h Xi_h(R) is a polynomial in tau (the numerator of Omega is), so it vanishes identically. QED

The reductions of brief s.3 (Corollary E, Theorem D5, the parity argument, Theorem F) are not needed.

## 7. The one-leaf function: Theorem 5, (6.1), the gap, Zhat

**Theorem 5.** For integers N >= 1, M >= 0 and parameters a, b,

    [v^{-1}]_{v=oo} prod_{m in M(M)} (v + 1 + m a) / prod_{m in M(N)} (v + m b) = (M!/(N-1)!) [phi^{M-N+1}] e^phi S(a phi)^{M+1} S(b phi)^{1-N}.

*Proof.*
1. The left side is the sum of the residues at the finite poles. That sum is the divided difference P[x_1..x_N] of P(v) = prod_{m in M(M)} (v + 1 + m a) at the nodes -m b, m in M(N). The nodes are equally spaced and symmetric about 0, so P[x_1..x_N] = (2 sinh(b d/2)/b)^{N-1} P(0)/(N-1)! = (Sc(b d) d)^{N-1} P |_{v=0} / (N-1)!.
2. For a polynomial P and a power series G, [G(d_v) P](0) = [P(d_phi) G](0).
3. P(d) = e^{-phi} prod_{m in M(M)} (d + m a) e^{phi}.
4. For every power series H, [prod_{m in M(M)} (d + m a) H](0) = M! [phi^M] H(phi) S(a phi)^{M+1}.
   - It suffices to take H = e^{s phi}, i.e. to show prod_{m in M(M)} (s + m a) = M! [phi^M] e^{s phi} S(a phi)^{M+1}.
   - This follows from prod_{m in M(M)} (s + m a) = a^M M! [w^M] (1+w)^{s/a + (M-1)/2} and the substitution 1 + w = e^{a phi}.

Steps 1-4 together give (M!/(N-1)!) [phi^M] e^phi (Sc(b phi) phi)^{N-1} S(a phi)^{M+1}, which is the right side. QED

**Corollary 6 ((6.1) for all h).** For all h >= 0, all j >= 1 (and, as a polynomial identity, all j), and all tau:

    Sum_{h>=0} sigma^h Xi_h(j) lambda^{2h} = (j-1)! [phi^{j+1}] e^phi S(e1 lambda phi)^{1+j e2} S(e2 lambda phi)^{1+j e1},

that is,

    sigma^h Xi_h(j) = prod_{l=1}^{2h-2} (j - l) [lambda^{2h}] S(e1 lambda)^{1+j e2} S(e2 lambda)^{1+j e1}.   (6.1)

*Proof (integer tau >= 1 first).*
- Apply Theorem 1 (n = 1) with hb = sigma lambda, so that hb/e2 = tau lambda and hb/tau = e2 lambda.
- Apply Theorem 5 with N = n1 = tau j, M = n2 = e2 j, a = tau lambda, b = e2 lambda. This gives [x^j] F_{g,1} = j^{-2} sigma^{-2g} [lambda^{2g}] (n2!/(n1-1)!) [phi^{j+1}] e^phi S(tau lambda phi)^{1+n2} S(e2 lambda phi)^{1-n1}.
- By Theorem 3, sigma^g Xi_g(j) = sigma^{2g-1} [x^j] F_{g,1}/Pi_j.
- Pi_j = (n2-1)!/(n1! (j-1)!), so the prefactor is n1! n2! (j-1)!/(sigma j^2 (n1-1)! (n2-1)!) = (j-1)!.
- S is even, so S(tau lambda phi) = S(e1 lambda phi). Also 1 - n1 = 1 + j e1 and 1 + n2 = 1 + j e2.
- Expand S(e1 l phi)^{1+je2} S(e2 l phi)^{1+je1} = Sum_h c_h(j) (l phi)^{2h}. Then [phi^{j+1}] gives Sum_h c_h lambda^{2h}/(j+1-2h)!, and (j-1)!/(j+1-2h)! = prod_{l=1}^{2h-2} (j - l), which is 0 for 1 <= j <= 2h-2.
- Both sides of (6.1) are polynomials in tau, and in j of degree <= 3h-2, so the identity holds for all tau and all j. QED

*Consequences.*
- Xi_h(j) = 0 for 1 <= j <= 2h-2, for all tau. This is the one-leaf part of Lemma A and (OL), and contains the even-j statement asked for in the route brief.
- Route E's generating function (E_formulas report s.1.3) holds for all h:
  Zhat_j(u) = Sum_{h>=1} sigma^h Xi_h(j) u^{2h-2}/(2h-2)! = Res_{l=0} l^{-3} (1 + u/l)^{j-1} S(e1 l) S(e2 l) B(l)^j dl, with B = S(e1 l)^{e2} S(e2 l)^{e1}.
  (The coefficient of u^{2h-2} is C(j-1, 2h-2) [l^{2h}] S1 S2 B^j = sigma^h Xi_h(j)/(2h-2)!.)
- Theorem 1 at n = 1 is a dilaton-point analogue of the one-part Marino-Vafa formula: a residue of a ratio of centred products.

**Corollary 7 (gap, r = 0).** For h >= 2, F_h := Int_{Mbar_h} K Omega = Xi_h(empty) satisfies

    sigma^h F_h = -(2h-3)! [lambda^{2h}] S(e1 lambda) S(e2 lambda).

*Proof.*
- String equation: forget a point carrying lam_0 = 1.
  - A leaf lam_j contributes D lam_j = j lam_j.
  - A dilaton leaf contributes DT = (T(psi) - T(0))/psi = -psi/(1-psi) = lam_0 - lam_1.
  - Hence Xi_h(0 cup R') = |R'| Xi_h(R') + Xi_h(0 cup R') - Xi_h(1 cup R'), i.e. Xi_h(1 cup R') = |R'| Xi_h(R'). In particular Xi_h(1) = 0.
- Dilaton equation: Xi_h[psi] - Xi_h[T] = (2h-2) F_h, where Xi_h[f] is the value with f inserted at one point. Since psi - T = psi/(1-psi) = lam_1 - lam_0, (2h-2) F_h = Xi_h(1) - Xi_h(0) = -Xi_h(0).
- (6.1) at j = 0 gives sigma^h Xi_h(0) = (2h-2)! [lambda^{2h}] S(e1 lambda) S(e2 lambda). QED

This is the formula of proof_M_residue.md s.5 (F_h(q_0), checked there for h <= 4), now for all h. Its limits are chi(M_h) = B_{2h}/(2h(2h-2)) as sigma -> oo (Harer-Zagier) and the lambda_g value as sigma -> 0.

## 8. Uniformity in the genus

| Step | Why it is uniform in g |
|---|---|
| Theorem 1 | ABDKS (17) is an identity of formal series in hb, valid for all (g, n) at once. Steps (i)-(iv) are done at each fixed order in hb and do not depend on g. |
| Theorem 2 | A homogeneity property of the all-genus polynomial R_mu(c, hb). |
| Theorem 3 | Eynard's formula is an all-genus statement. The moment lemma, the degree involution and the pushforward are uniform in g. |
| Corollary 4 | Uses only Theorems 2 and 3 plus polynomiality in tau. No inductive hypothesis (Lemma A in lower genera, lower levels at genus G) is used anywhere. |
| Theorem 5, Corollaries 6 and 7 | Uniform in h and j. |

Hypotheses used: integer tau >= 1 (so that the products are finite; the statements then extend because they are polynomial in tau), 2g-2+n > 0, and sigma != 0.

## 9. The fermionic MV route (as briefed)

**F1 (proved; Murnaghan-Nakayama).** Let R(p) = Sum_lambda W_lambda s_lambda(p), with W_lambda = e^{tau x f_2(lambda)} s_lambda(q^rho), the framed one-leg vertex weight of Route A's mvlib.
- The coefficient of p_j at p = 0 is (1/j) Sum over hooks lambda of size j of (-1)^{leg(lambda)} W_lambda.
- In the background p_1 = y/(sigma x), the hooks are replaced by j-border strips: [p_j] R = (1/j) Sum_a (p_1^a/a!) Sum_{|lambda| = j+a} W_lambda chi_lambda(j, 1^a).
- This is the generating function of the MV integrals with one j-leg and unit legs that enters Route A's Laguerre formula (Prop. 2.1).

**F2 (the obstruction).**
- At fixed y, the genus expansion of [p_j] log R is a semiclassical expansion. Its leading order produces the Lambert-type relation y = u e^{-u}.
- Evaluating at the dilaton point u = oo, on the other sheet, needs every order of that expansion.
- I could organise the expansion order by order, but found no closed all-order control on the MV side alone. In effect this is the quantum-curve problem.
- This is why the route was redirected to the x-y swap.

**F3 (remark).**
- ABDKS's proof of (17) is the Cauchy determinant (their (30)), i.e. a free-fermion correlator.
- The trivial dual means that the hbar-dressing of the two log points is all there is.
- Theorem 1 is therefore exactly the closed all-genus semiclassics that F2 lacked. A TR-free derivation from E-operators or Orlov-Scherbin formulas looks possible but was not completed (next steps).

## 10. Verification (exact, Perl Math::BigRat)

There are no mismatches anywhere. Raw outputs are in runs_all_checks.txt.

| # | Check | Script | Range | Values | Mismatches |
|---|---|---|---|---|---|
| V1 | Theorem 5 (finite products) = (6.1) | xy_onepoint.pl | tau=1 (j<=6, g<=3); tau=2, 3 (j<=8, g<=5) | 120 | 0 |
| V2 | n=1 formula for rational tau (formal xhat, no finite products) = (6.1) | xy_onepoint_generic.pl | tau = 1, 1/2, -2/7; g = 1..3 | 42 | 0 |
| V3 | Thm 1 + P2 vs Route E engine, n=1 | xy_npoint.pl | h = 2..6 at tau = 1, 2, 4 (j <= 3h) | 180 | 0 |
| V4 | same, genus 7 (new table hodge_engine/engine_h7_r1.txt, 3779 s) | xy_npoint.pl | tau = 1, 2, 3, 4; j <= 21 | 84 | 0 |
| V5 | (6.1) vs genus-7 engine | check61_engine.pl | 7 values of tau, j = 0..21 | 154 | 0 |
| V6 | n=1, genus 1 | xy_npoint.pl | tau = 2 | 6 | 0 |
| V7 | n=2 vs engine | xy_twopoint.pl, xy_npoint.pl | g=1 (tau 1,3): 30; g=2,3,4 (tau 1,2,3): 252; g=5,6 (tau 1,2): 144 | 426 | 0 |
| V8 | two-point difference form | xy_twopoint_diff.pl | g=2 tau=1; g=3 tau=3 | 36 | 0 |
| V9 | n=3 vs engine | xy_npoint.pl | g=1 (tau 2): 20; g=2 (tau 1): 10; g=3 (tau 1,2): 40; g=4 (tau 1,2): 30; g=5 (tau 1): 20 | 120 | 0 |
| V10 | n=4 vs engine | xy_npoint.pl | g = 2, 3, 4 at tau = 1 | 45 | 0 |
| V11 | n=5 vs engine | xy_npoint.pl | g = 2 at tau = 1 | 21 | 0 |
| V12 | genus 0 vs Gamma(J)/Gamma(J-r+3) | genus0_check.sh | r = 3, 4, 5; tau = 1, 2, 3 | 132 | 0 |
| V13 | Thm 1 vs the exact EO recursion (referee code status_trC1.pl), tau = 1 | cmp_tr_xy.pl | g = 2..7, n <= 3: 307; g <= 4, n = 4, 5: 89 | 396 | 0 |
| V14 | Thm 1 vs the main session's Mathematica EO recursion (hvc_tr_tau2_g3n3.txt, hvc_tr_tau3_g3n3.txt) | xy_npoint.pl | tau = 2, 3; g <= 3, n <= 3, first layer | 15 | 0 |
| V15 | (6.1) vs exact EO recursion, tau = 1 | status_cmp61.pl | h = 2..7 | 90 | 0 |
| V16 | (6.1) vs exact EO recursion, tau = -1/2 (odd j <= 45) | status_trhalf.pl, status_cmp61.pl | h = 2..7 (126); h = 2..8 (161) | 287 | 0 |
| V17 | gap formula vs engine | gap_check.pl | h = 2..7 at tau = 1, 2, -3/7, 5/2 | 24 | 0 |
| V18 | moment lemma (now proved) | moments_check.pl | tau = 2 (k <= 6); tau = 5, 1/3 (k <= 4) | 17 | 0 |
| V19 | Thm 1 vs lemmaA/data q-tables (converted by q_to_engine.pl) | xy_npoint.pl | r=1: h = 2..4 at tau = 1..3 (81); r=2: h = 2..4 at tau = 1..3, parts <= 5 (135) | 216 | 0 |
| V20 | homogeneity (no hb^p with p > abs(mu) in R_mu(1, hb)) | xy_npoint.pl | every run above | -- | 0 violations |

Notes on the table.
- Engine total: 1072 values; the engine is TR-free (Mumford GRR plus Witten-Kontsevich numbers).
- The genus-7 framings tau = 5..8 (V4) were still running at submission. They are redundant with V5 and Theorem 5.
- An early version of xy_threepoint.pl had a sign error in two cycle factors. This was a script bug, caught by the engine comparison and fixed; the superseded runs are marked as such in runs_all_checks.txt.
- Route A's MV tables are not compared directly. The Hodge engine was validated against the one-part MV formula by Route E, and the lemmaA/data tables are compared in V19.

## 11. Consequences for the project

- Lemma A (M_irr(h)) holds for all h, modulo the two literature inputs.
  - Combine with Lemma C (leaf_structure_proof.md, proved) and Corollary R (s.13 there; the deduction "Lemma A => Corollary R" is the main session's and was not re-checked here).
  - Then the all-genus regularity theorem for the conifold graph sum follows modulo the same inputs. TR enters only on the Hodge side.
- Status of the project lemmas:
  - Lemma B (P2) is Theorem 3; Lemma D is Theorem 2.
  - Theorem 3.1 of proof_M_hurwitz.md is the tau = 1 case of Theorem 2.
  - (6.1), the gap formula and Route E's Zhat_j are now theorems, modulo literature, in all genera.
- Not obtained here:
  - closed multi-leaf forms ((SN), (NR) of Route E). They are contained in Theorem 1, but are not extracted;
  - the extra zero j = 2h in the Hurwitz limit;
  - a TR-free proof.

## 12. Literature

- [ABDKS-log] A. Alexandrov, B. Bychkov, P. Dunin-Barkowski, M. Kazarian, S. Shadrin, "Log topological recursion through the prism of x-y swap", arXiv:2312.16950v3 (Jan 2025). Used: Prop. 3.1 (eq. (17)), Rem. 2.6, 3.2, 3.4, and s.3.4-3.5 (eqs. (40)-(50)).
- [ABDKS-swap] same authors, "A universal formula for the x-y swap in topological recursion", arXiv:2212.00320. This is the meromorphic x-y swap; Theorem 2.8 of [ABDKS-log] extends it.
- [Eyn14a] B. Eynard, "Invariants of spectral curves and intersection theory of moduli spaces of complex curves", Commun. Number Theory Phys. 8 (2014) 541-588.
- [DOSS] P. Dunin-Barkowski, N. Orantin, S. Shadrin, L. Spitz, Comm. Math. Phys. 328 (2014) 669-700.
- BKMP for C^3: V. Bouchard, M. Marino; L. Chen; J. Zhou.
- R. Kaufmann, Yu. Manin, D. Zagier (kappa classes and dilaton pushforward).
- A. Okounkov, R. Pandharipande, Geom. Topol. 8 (2004) (for section 9).

Normalisations of the two cited results were fixed against the data in section 10: Theorem 1 against the exact EO recursion, and Theorem 1 + Theorem 3 against the TR-free engine.

## 13. Files (all in HVC\lemmaA\round2\R2_onepoint\)

| File | Contents | Usage |
|---|---|---|
| xy_npoint.pl | Theorem 1 for general n, with the homogeneity check; compares with engine-format tables | perl xy_npoint.pl TAU G N JMAX [TABLE] |
| xy_onepoint.pl | Theorem 5 algebra vs (6.1) | perl xy_onepoint.pl TAU JMAX G |
| xy_onepoint_generic.pl | n = 1 for rational tau | perl xy_onepoint_generic.pl TAU G JMAX |
| xy_twopoint.pl, xy_twopoint_diff.pl, xy_threepoint.pl | n = 2, 3 versions | |
| cmp_tr_xy.pl | Theorem 1 vs status_trC1.pl output | perl cmp_tr_xy.pl TRFILE NMAX |
| moments_check.pl | moment lemma | |
| gap_check.pl | Corollary 7 vs engine | perl gap_check.pl TABLE h |
| genus0_check.sh | genus-0 check | |
| q_to_engine.pl, qeng_h<h>_r<r>.txt | data q-tables in engine format | perl q_to_engine.pl h r |
| runs_qtables_r1.out, runs_qtables_r2.out | V19 outputs | |
| queue_multileaf.sh, runs_np_g*_n*_t*.out, runs_cmp_tr_n45.out | multi-leaf runs | |
| runs_all_checks.txt | raw outputs of every check | |
| hodge_engine/ | copy of Route E's engine; engine_h7_r1.txt (new), run_h7_r1.log | perl xi_engine.pl 7 1 |
| engine_h7_r1.txt | copy of the new genus-7 table | |
| trcheck/ | copies of status_trC1.pl, status_trhalf.pl, status_cmp61.pl; outputs trC1_g3.out, trC1_g7.out, trC1_g7_n45.out, trhalf_g7.out, trhalf_g8.out | |
| lib/MP.pm | multivariate polynomials over BigRat | |
| lib/mvlib.pm | copy of Route A's library | |
| q_extract.pl, q_explore.pl | exploration of the data q-tables | |

