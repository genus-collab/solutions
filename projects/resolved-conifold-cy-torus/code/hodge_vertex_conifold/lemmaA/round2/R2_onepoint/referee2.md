# Lemma A round 2: referee2 R2 (saved by the main session from the workflow journal, 2026-10-31)

Subagents cannot write report files, so this file holds the agent's structured return value: the fields first, then the full report text (report_text) verbatim.

## confirmed

- Theorem 1 (closed n-point formula on C_tau): the derivation from ABDKS Prop. 3.1 (eq. 17; Rem. 2.6, 3.2, 3.4) was checked line by line. The ABDKS statement and its hypotheses were checked against arXiv:2312.16950v3 (IMRN 2024), and all hypotheses hold for C_tau when tau is not in {0, -1}. Own implementations (n = 1, 2, 3) agree with TR-free engines through genus 9 (n = 1), 8 (n = 2) and 8 (n = 3), and with exact EO at tau = -1/2 through genus 9.
- Theorem 2 (degree bound = Lemma D for every integer tau >= 1): the homogeneity argument is correct, including the n = 1 hb^-1 term.
- Theorem 3 (dictionary P2 for all g, n, tau): the reparametrisation, (p,q,r) = (1/tau, -1/(1+tau), 1), q_k = (-1)^k kappa, s_i = (-1)^(i+1)/i, the Hodge factor, the degree involution and sigma^(1-g) were all re-derived. The normalisation is confirmed TR-free at genus 2..9.
- Moment lemma: correct and immediate.
- Corollary 4: Lemma A = M_irr(h) for every h, and the conclusion of (OL) for every G without induction. Consistent with the TR-free engine at genus 8 (r = 1, 2, 3, every component: 1737 zero values, with the first layer non-zero) and genus 9 (r = 1), and with Route A's KP ranks (nullity 0, rerun at G = 7, 8, 9).
- Theorem 5 (divided-difference evaluation): the proof is correct. Theorem 1 (n = 1) equals (6.1) for tau = 1..5, g <= 12, j <= 40 (2400 values).
- Corollary 6 ((6.1) for all h, j, tau): exact identity in (tau, j) against the TR-free engine at genus 7, 8 and 9 (518 coefficients). Its Hurwitz and lambda_g components match Route B's Theorem B and the Faber-Pandharipande lambda_g formula for all h (analytically) and h <= 12 (exactly).
- Corollary 7 (gap formula; the string and dilaton relations): correct. The engine agrees at genus 2..9, and the sigma -> oo limit equals chi(M_h).
- Two-point difference form: re-derived, and the independent implementation agrees with the engine at genus 6, 7 and 8.
- Zhat_j generating function: algebraically equivalent to (6.1).
- The R2 formulas satisfy the TR-free Marino-Vafa cut-and-join relations (R = (j): G = 7..10; R = (j1,j2): G = 8, 9, 10) and the full DJKM/KP Hirota equations at the dilaton point (G = 7..10), with no failure.
- F1 (Murnaghan-Nakayama): correct; not used in the main line.

## issues

### [1]

#### description

The proof uses TR on the Hodge side: the LogTR x-y swap of ABDKS, EO on C_tau and Eynard's kappa-class formula. Brief s.0 (the user's clarification of 2026-10-31, also in the user's memory) allows exactly this. The workflow's computed task text still says 'no topological recursion of any kind'. Under the stale wording Theorems 1 and 3, and with them Cor 4, 6 and 7, would be inadmissible.

#### location

R2 summary and report s.0 ('Use of TR'); computed task text

#### severity

minor

#### suggested_fix

Keep R2's statement. In the main notes record Lemma A as 'proved modulo TR for the Hodge CohFT (ABDKS 2024 Prop. 3.1; Eynard 2014)', not as TR-free. A TR-free derivation of Theorem 1 via the closed n-point formulas for Orlov-Scherbin (hypergeometric) KP tau functions, an allowed input, is the natural follow-up.

### [2]

#### description

ABDKS state (50) without proof, quoting it from [Eyn14]. Its Hodge-class form relies on the C^3 case (44), which ABDKS attribute to BMn08, Zho09 and Che18 (BKMP for C^3). The first expression of the chain (49) in ABDKS has an index slip: it is printed with xi_k, while the consistent form, and q_0 = y'(0), require xi_{k+2}. R2 uses the consistent form q_k = -res xi_{k+1} y dx, and the data confirm it.

#### location

report s.5 (Theorem 3, literature input) and s.12

#### severity

minor

#### suggested_fix

Cite Eynard, CNTP 8 (2014), for the one-branch-point kappa-class theorem, and BKMP for C^3 (L. Chen; J. Zhou), with theorem numbers. Add a footnote on the (49) index slip.

### [3]

#### description

(6.1) vs the genus-7 engine is called 'identically in j and tau' on the basis of 7 framings. sigma^7 Xi_7(j) has degree 7 in sigma, so 8 sigma-values are needed, or the known lambda_g constant term. The check is therefore one condition short of an identity.

#### location

report s.10, V5 (check61_engine.pl at genus 7)

#### severity

minor

#### suggested_fix

Compare componentwise. Done here: at genus 7, 8 and 9 every (k, j^b) coefficient agrees (132, 171 and 215 coefficients).

### [4]

#### description

The h = 0 term Xi_0(j) = 1/(j(j+1)) is the dilaton-point convention (Route A, Lemma 3.1). It is not an intersection number, and Theorem 1 does not cover the unstable (0,1) case. Genus 1 is covered, since (1,1) is stable.

#### location

report s.7, Corollary 6 ('for all h >= 0')

#### severity

minor

#### suggested_fix

State the generating-function identity for h >= 1, with the h = 0 term as the normalisation convention.

### [5]

#### description

The formula Xi_g(mu) = sigma^{g-1}[x^mu]F_{g,n}/prod Pi_{mu_i} 'for every tau with sigma != 0' divides by Pi_{mu_i}, which vanishes at some framings (tau = -1/2 with even parts). The proof of Cor 4 uses only integer tau >= 1, where Pi > 0.

#### location

report s.5, statement of Theorem 3

#### severity

cosmetic

#### suggested_fix

Write [x^mu]F_{g,n} = sigma^{1-g} prod Pi_{mu_i} Xi_g(mu).

### [6]

#### description

R_mu(c, hb) is called a homogeneous polynomial. For n = 1 it contains the term c^{|mu|+1} hb^{-1}, the unstable (0,1) part. The degree argument is unaffected.

#### location

report s.4, proof of Theorem 2

#### severity

cosmetic

#### suggested_fix

Say 'polynomial in c, Laurent in hb with hb-exponent >= -1, homogeneous of degree |mu|'.

### [7]

#### description

The remark that both sides are continuous in gamma, or a Moebius change, is unnecessary. ABDKS (40)-(42) only require alpha, beta, gamma to be distinct, so gamma = 0 is admissible.

#### location

report s.5, proof step (a)

#### severity

cosmetic

#### suggested_fix

Drop the remark.

### [8]

#### description

engine_h7_r1.txt is textually identical, after sorting lines, to the pre-existing E-referee table E_formulas/referee/ref_engine_h7_r1.txt from the independent RefHodge engine. So it is not new data; its value is that two engines agree.

#### location

summary, V4 and s.13 ('the genus-7 one-leaf table was computed here')

#### severity

cosmetic

#### suggested_fix

Say so.

### [9]

#### description

Brief s.2(g) records M_irr(7) as complete (Route E referee engine, 194 points).

#### location

summary and s.6 ('partial M_irr(7)')

#### severity

cosmetic

#### suggested_fix

Update the wording.

### [10]

#### description

R2 labels Eynard's CNTP paper [Eyn14a]. In ABDKS that paper is [Eyn14], and [Eyn14a] is Eynard's ICM 2014 survey.

#### location

report s.12 (literature keys)

#### severity

cosmetic

#### suggested_fix

Use ABDKS's keys, or full references, to avoid confusion.


## notes

I found no failure. The all-genus argument (ABDKS x-y swap, then scaling, plus Eynard's formula, plus polynomiality in tau) is complete and non-inductive, and both cited theorems apply to C_tau.

New data written to referee2/eng: TR-free engine tables at genus 8 (r = 1, 2, 3) and genus 9 (r = 1). Structural TR-free tests (MV cut-and-join, full KP Hirota) reach genus 10. No genus-10 engine table was computed (estimated 6-7 h).

Admissibility depends on the user's TR ruling for the Hodge CohFT (brief s.0). The computed task text's 'no TR of any kind' is stale. All issues are minor or cosmetic.

No Mathematica checks are requested.

## refuted


## route

R2 (round 2): the one-leaf function / x-y swap on the Hodge curve C_tau -- second, independent referee

## verdict

sound

## report_text

# Second referee report on Route R2 (round 2): Lemma A, (6.1) and the gap formula in all genera via the x-y swap on C_tau

Referee 2 (independent), 2026-10-31 to 11-01. HVC = G:\My Drive\Projects\Old\refGW\hodge_vertex_conifold.
- My scripts and outputs are in HVC\lemmaA\round2\R2_onepoint\referee2\.
- Everything is exact (Perl Math::BigRat). No wolframscript was used.
- Nothing outside referee2\ was modified. Route A, Route E and author code were copied in before being reused.

## 0. Verdict: SOUND (modulo the two cited theorems, which apply)

I tried to break the all-genus claim at every step and could not. The argument is:

    Thm 1 (ABDKS x-y swap, Case I, trivial dual)  =>  Thm 2 (scaling / homogeneity)
    Thm 2 + Thm 3 (Eynard's kappa-class formula = P2) + polynomiality in tau  =>  Cor 4: Xi_h(R) = 0 for every R (parts >= 1), D(R) <= 2h-3, every h.

- **No induction.** There are no base cases or inductive hypotheses to check, and no circularity: the data are used only to confirm normalisations, never as input.
- **Each step is correct as written.** Only cosmetic points remain (s.5).
- **Both literature inputs are stated correctly, and their hypotheses hold for C_tau.** I checked them against the source, arXiv:2312.16950v3 = IMRN 2024 no. 21, 13461-13487 (s.3).
- **Lemma A for every h.** Cor 4 contains Lemma A = M_irr(h) for every h. It also contains the conclusion of (OL) for every G, without the hypotheses of (OL). The reductions of brief s.3 (Cor E, Thm D5, parity, Thm F) are not needed.
- **No discrepancy anywhere in the tests at G = 7, 8, 9, 10** (s.1). The tests use:
  - four NEW TR-free Hodge-engine tables: genus 8 with r = 1, 2, 3, and genus 9 with r = 1;
  - the Marino-Vafa cut-and-join, through genus 10;
  - the full DJKM/KP Hirota equations of the MV tau function at the dilaton point, through genus 10;
  - Route A's KP ranks, rerun at G = 7, 8, 9;
  - the exact EO recursion at tau = -1/2, through genus 9.

  This is about 6800 exact comparisons and relations.
- **Not TR-free.** The result rests on ABDKS Prop. 3.1 and on Eynard's formula, i.e. on TR for the Hodge CohFT. Brief s.0 allows this (s.4).

## 1. Every step asserted for all G, and its tests at G = 7, 8, 9, 10

- "engine" means a TR-free Hodge engine: Mumford GRR plus Witten-Kontsevich, with the dilaton leaves kept as explicit points. This is RefHodge, copied from E_formulas/referee into referee2/eng.
- "identity" means that every (component k, power j^b) coefficient is compared, so the check is an identity in tau and j.

| # | Step asserted for all G | Test (script) | G = 7 | G = 8 | G = 9 | G = 10 |
|---|---|---|---|---|---|---|
| 1 | Thm 1 (n = 1) + Thm 5: [x^j]F_{g,1} = j^-2[hb^2g][v^-1]Q_j equals (6.1), integer tau | ref_xy1.pl (own code: [w^{j+1}] of the product series); tau = 1..5, j <= 40, every g <= 12; 2400 values (1740 non-zero) | 0 mism. | 0 | 0 | 0 |
| 2 | Thm 1 + Thm 3 = (6.1) vs engine, identity | ref_cmp61_engine.pl | 132/132 coeff. (E-referee table) | 171/171, NEW table h8 r1 (899 s) | 215/215, NEW table h9 r1 (4159 s) | not run (est. 6-7 h) |
| 3 | Thm 1 + Thm 3 at n = 2 vs engine | ref_xy2.pl (own code for the 2-cycle residue) | tau=1: 66 pairs (25 non-zero); tau=2: 55 (16) | tau=1: 105 (49); tau=2: 78 (25); NEW table h8 r2 | - | - |
| 4 | Thm 1 + Thm 3 at n = 3 vs engine | ref_xy3.pl (own divided-difference code) | 120 triples (50 non-zero) | 120 (31); NEW table h8 r3 | - | - |
| 5 | String equation across n (R2 values on both sides) | ref_string.pl | - | Xi(m1,m2,1) = (m1+m2)Xi(m1,m2): 36 | Xi(m,1) = m Xi(m): 26 (tau=2) | 30 (tau=1) |
| 6 | MV cut-and-join, R = (j) (Route A Thm B; TR-free; validated on engine data at G = 2..7, 90 relations) | ref_cj.pl, with (6.1) + Thm 1 (n = 2) | 23 relations | 26 (tau=1) + 26 (tau=2) | 29 | 32 (tau=1) + 32 (tau=2) |
| 7 | MV cut-and-join, R = (j1,j2) (uses Thm 1 at n = 3 in genera G and G-1; validated on engine data at G = 2, 3, 4) | ref_cj2.pl | - | 45 | 66 | 33 (j1 = 10..12) |
| 8 | Full DJKM Hirota (KP) equations at the dilaton point, all genus splittings (Route A's generator; validated on engine data at G = 2, 3, 6, 7) | kp/ref_hirota.pl, W = 4, R = (j), j = 2G-3..2G+3 | 21 eqs | 21 | 21 | 21 |
| 9 | Thm 2 / Cor 4: vanishing in the range, TR-free, every component | engine tables (ref_mirr_engine*.pl) | r = 1 | r = 1, 2, 3: 126 + 504 + 1107 values all zero; first layer D = 14 non-zero in every component (9/9, 72/72, 216/216) | r = 1: zeros at j = 1..16 in every component (item 2) | - |
| 9' | The same, via Route A's KP ranks (TR-free) | kp/a3_kp_nullity(2).pl, rerun | rank 1358/1358 | 104/104 | 135/135 | 170/170 (stored run) |
| 10 | Cor 7, gap formula | engine, 6 framings each | OK | OK | OK | - |
| 11 | End components of (6.1): Hurwitz (k = 0) and lambda_g (k = h) | vs Route B Thm B (ELSV + Toda) and the Faber-Pandharipande lambda_g formula; ref_limits.pl, h = 2..12, plus an analytic proof for all h (s.2.6) | OK | OK | OK | OK |
| 12 | TR side alone: Thm 1 (n = 1) = (6.1) vs the exact EO recursion on C_{-1/2} | tr/status_trhalf.pl (rerun, 3286 s), compared with my own (6.1) code | 0 mism. | 0 | 0 (21 values, 13 non-zero) | - |

Two things were not extended:
- **n >= 4.** Thm 1 at n >= 4 was not tested beyond the author's range (g <= 4). This is not needed for the proof, since Thm 1 is uniform in n.
- **Engine at genus 10.** No genus-10 engine table was computed. At genus 10 the R2 formulas are tested only structurally:
  - items 1, 5, 6, 7, 8 and 11;
  - the Route A KP rank, from a stored run (item 9').

## 2. Logic, step by step

### 2.1 Theorem 1 (re-derived)

(a) R2's restatement of ABDKS (17) is exact:
- the signs (-1)^n and (-1)^(n-1);
- the operator Sum_r (d (1/dx))^r [u^r];
- the exponent -u(S_A(u hb d_y) xhat - x), where S_A(u) = (e^{u/2} - e^{-u/2})/u is R2's Sc;
- xhat built with 1/S_A, which is R2's S;
- w, wbar = z +- u hb/2 (ABDKS Rem. 3.2).

(b) [x^mu] may be extracted as successive residues at v_i = oo, each taken over the field of rational functions of the later variables:
- omega_{g,n} has no poles on the diagonals or at v_i = oo;
- each hb-coefficient of (17) is a finite sum of rational differentials, so the residues can be applied term by term.

(c) Integration by parts: Res x^-mu d(eta/dX) = mu Res x^-mu eta.
- The operators in the other variables commute with Res_{v_1}.
- At each order in hb the u-dependence is polynomial, so u_i -> mu_i is legitimate.

(d) The identity mu tau Sc(mu hb d) S(hb d/tau) = Sum_{M(tau mu)} e^{m hb d/tau}, together with the analogous e2 identity, gives Q_mu. The overall sign is (-1)^(n-1).

(e) For n = 1, the u^-1 of the cycle factor contributes only an hb-independent term (E(0) = 1), which matters only for g = 0.

(f) Independent implementations:
- My three implementations (n = 1, 2, 3) of the resulting formula agree with TR-free data through genus 9 (n = 1), 8 (n = 2) and 8 (n = 3).
- The iterated residues are symmetric in the mu_i, as they must be. Example: xi_3(2,2,21) = xi_3(21,2,2) at genus 10.

### 2.2 Theorem 2 (correct)

- Q_mu(v; c) is homogeneous of degree mu in (v, c, hb), and each cycle factor is homogeneous of degree -1.
- The expansion region |v_1| >> ... >> |v_n| >> |c|, |hb| is invariant under scaling. So the coefficient of prod v_i^-1 is homogeneous of degree |mu| in (c, hb).
- That coefficient is a polynomial in c. Its hb-exponents are >= -1, and hb^-1 occurs only for n = 1 (the unstable (0,1) part).
- Hence [hb^p] R_mu(1, hb) = 0 for p > |mu|, i.e. [x^mu]F_{g,n} = 0 for |mu| <= 2g-3+n.
- This is exactly Lemma D of proof_status_mixed_scale.md s.2, at every integer tau >= 1. It uses nothing but Theorem 1.

### 2.3 Theorem 3 (re-derived; every intermediate value agrees with R2)

- **Coordinates.** Put v = tau + kappa z with kappa^2 = -sigma. Then dX = z dz/((1+alpha z)(1+beta z)) with alpha = kappa/tau, beta = kappa/(1+tau), gamma = 0. These are distinct, and gamma = 0 is admissible in ABDKS (40)-(42).
- **Parameters.** (p, q, r) = (1/tau, -1/(1+tau), 1), and 1/p + 1/q + 1/r = 0.
- **Leaves.** xi^A_k = (-1)^k (-kappa) xi_k.
- **Moments.** q_0 = y'(0) = kappa and q_1 = kappa alpha beta = -kappa, both checked by direct residue computation.
  - In general q_k = (-1)^k kappa, by the moment lemma.
  - The moment lemma is immediate: xi_k = Sum mu^k Pi_mu x^mu, and xi_k is rational with poles only at v = tau.
  - Hence s_i = (-1)^(i+1)/i.
- **Hodge factor.** Lambda(1/tau) Lambda(-1/(1+tau)) Lambda(1) = sigma^-g Lambda^v(-1) Lambda^v(-tau) Lambda^v(1+tau).
- **Degree involution.** It turns exp(Sum s_j kappa_j) into K = exp(-Sum kappa_j/j), and the Hodge product into (-1)^(3g) (-sigma)^g Omega. The powers of kappa and the signs then combine to sigma^(1-g).
- **KMZ pushforward.** With T(psi) = -psi^2/(1-psi), the KMZ pushforward gives s_j = -1/j with Arbarello-Cornalba kappa classes. The engine keeps the dilaton leaves as explicit points, so its agreement through genus 9 independently confirms the identification Xi = Int K Omega prod lam.

### 2.4 Corollary 4 (correct)

- (-sigma)^h Xi_h(R) is a polynomial in tau of degree <= 2h.
- It vanishes at every integer tau >= 1, where Pi_j > 0. Hence it vanishes identically.
- Lemma A and (OL) follow at once.
- No step uses Lemma A in a lower genus or at a lower level.
- Multisets R with unit parts are consistent with the string equation.

### 2.5 Theorem 5, Cor 6, Cor 7, the two-point form and Zhat (re-derived)

- **Theorem 5.** The proof uses a divided difference on a centred grid, the duality [G(d)P](0) = [P(d)G](0), and centred factorial products.
- **Cor 6.** The prefactor (M!/(N-1)!) / (sigma j^2 Pi_j) equals (j-1)!.
- **Cor 7.**
  - String: f -> (f - f(0))/psi sends lam_j to j lam_j and T to lam_0 - lam_1.
  - Dilaton: psi - T = lam_1 - lam_0.
- **Two-point form.** It follows from the 2-cycle residue. My implementation uses the equivalent expansion
  Sum_{n,l} (K hb)^(2n) C(2n+1+l, l) q_{m1-2n-1-l} q'_{m2+l+1}.
- **Zhat.** The u^(2h-2) coefficient is C(j-1, 2h-2) [l^(2h)] S1 S2 B^j, and the odd powers of u vanish.

### 2.6 All-genus analytic cross-checks added here

- **Hurwitz limit.** The leading sigma-coefficients of p_{2k} and p_{2k-1} (where p_m = e1^m + e2^m) are 2 and 2k-1. So the top component of (6.1) is
  (j-1)...(j-2h+2) [t^(2h)] S^2 H^j,  with H = S e^{-tS'/S}.
  This is Route B's Theorem B (ELSV + Toda, TR-free), for every h.
- **lambda_g limit.** At sigma = 0, (6.1) gives [l^(2h)] S(l) prod(j-l) = (-1)^h b_h prod(j-l), with b_h = [t^(2h)] (t/2)/sin(t/2). This is the Faber-Pandharipande lambda_g formula.
- **Gap as sigma -> oo.** [y^(2h)] S(y)^2 = (1-2h) B_{2h}/(2h)!, so F^(0)_h = B_{2h}/(2h(2h-2)) = chi(M_h) (Harer-Zagier). R2's remark is right.

## 3. Literature: exact statements and hypotheses

### 3.1 ABDKS, arXiv:2312.16950v3 (= IMRN 2024 no. 21, 13461-13487)

Read from the arXiv HTML. The parts R2 uses:
- **Standing assumptions** (s.2.2-2.3): Sigma is a compact Riemann surface; dx and dy are meromorphic, and residues are allowed; the zeros of dx and of dy are simple and disjoint; dy is regular at the zeros of dx, and dx at the zeros of dy.
- **Def. 2.4:** a LogTR-vital point is a simple pole of dy, with residue 1/alpha, at which dx has no pole.
- **Def. 2.5:** omega_1^(0) = -y dx, omega_2^(0) = B, blobbed TR plus the log projection property.
- **Rem. 2.6:** a simple pole of dy at which dx also has a pole is not vital and does not contribute.
- **Thm 2.8:** the x-y swap for LogTR, proved in s.5.
- **s.3.1:**
  - Sigma = CP^1;
  - the dual TR is trivial iff dy has no zeros;
  - Case I is y = z;
  - (14): x = Sum log(z-a_i)/alpha_i + R(z), where the a_i are the vital points and R is unrestricted;
  - (15)-(16): definitions of xhat and of w, wbar.
- **Prop. 3.1 (17):** stated as "a direct corollary of Equation (13)", together with the Cauchy determinant (30).
- **Rem. 3.2 (18)-(20)** and **Rem. 3.4:** LogTR = TR when y = z.

Application to C_tau (z = v, x = X = tau log v - (1+tau) log(1+v), y = v):
- dX has a single simple zero, at v = tau. dv has no zeros and is regular at v = tau.
- The vital points are a_1 = 0 (1/alpha_1 = tau) and a_2 = -1 (1/alpha_2 = -(1+tau)).
- v = oo is not vital, since dv has a pole there (Rem. 2.6). R = 0.
- Every hypothesis holds for every tau not in {0, -1}. Nothing requires the alpha_i to be integral.
- omega_{0,1} = -v dX = y dX with y = -v is the project's convention. Exact-EO agreement confirms that the TR conventions match: at tau = 1 through genus 7 (author), and at tau = -1/2 through genus 9 (here).

### 3.2 Eynard's formula

- ABDKS (50) is quoted there from [Eyn14] = B. Eynard, CNTP 8 (2014) 541-588, without proof.
- Its hypotheses: x as in (40)-(41) with alpha, beta, gamma distinct, and y any power series at z = 0 with y'(0) != 0.
- ABDKS (44) is the case of a special y, giving triple Hodge integrals with no kappa classes. ABDKS attribute it to [BMn08], [Zho09] and [Che18], i.e. BKMP for C^3 (L. Chen; J. Zhou).
- So Theorem 3 rests on Eynard's general one-branch-point kappa-class theorem, combined with the C^3 identification of the Hodge class.
- Index slip in ABDKS (49): the first expression of the chain is printed with xi_k. Consistency with the other two expressions, and with q_0 = y'(0), requires xi_{k+2}. R2 uses the consistent form, q_k = -res xi_{k+1} y dx.
- The normalisation of the composite (Thm 1 + Thm 3) is confirmed TR-free at every genus from 2 to 9.

### 3.3 Other inputs

- **Kaufmann-Manin-Zagier / Arbarello-Cornalba** dilaton pushforward: correct. Replacing psi by its pull-back is valid because the forgotten points carry psi^(>=2).
- **Lagrange inversion:** xi_0 = Sum Pi_mu x^mu, checked.
- **Murnaghan-Nakayama** (F1): correct, but not used in the main line.
- **Harer-Zagier:** used only in a remark.

## 4. Rule compliance (TR)

- **Where TR is used.** Only for the Hodge CohFT on C_tau: the LogTR x-y swap, the EO recursion and Eynard's formula.
- **Allowed by brief s.0.** Brief s.0 (the user's clarification of 2026-10-31, also recorded in the user's project memory) allows exactly this.
- **Nothing touches the conifold.** No step involves the resolved conifold x C^2 theory, its free energies, its polar part, its graph sum, or input (H).
- **Stale wording.** The workflow's computed task text still says "no topological recursion of any kind". That is a stale summary of brief s.0.
  - Under the user's ruling R2 is admissible.
  - Under the stale wording, Theorems 1 and 3 would be inadmissible.
  - The main session should record Lemma A as proved modulo TR for the Hodge CohFT (ABDKS 2024; Eynard 2014), not as TR-free.
- **My own checks** are TR-free, except the EO cross-check at tau = -1/2. Brief s.0 lists that one as an allowed cross-check.

## 5. Issues (none affects the theorem)

1. (minor) Admissibility documentation; see s.4.
2. (minor) Citation precision for Theorem 3.
   - Cite Eynard's 2014 theorem together with BKMP for C^3 (Chen; Zhou), not only "ABDKS (50)".
   - Note the index slip in ABDKS (49).
3. (minor) V5 claims "identically in j and tau" at genus 7 from 7 framings.
   - sigma^7 Xi_7(j) has degree 7 in sigma, so 8 values (or the known lambda_g constant term) are needed.
   - The componentwise comparisons made here, at genus 7, 8 and 9 with every coefficient, close the gap.
4. (minor) Cor 6 says "for all h >= 0".
   - The h = 0 term Xi_0(j) = 1/(j(j+1)) is the dilaton-point convention (Route A, Lemma 3.1).
   - Theorem 1 does not cover (g, n) = (0, 1). Genus 1 is covered, since (1,1) is stable.
5. (cosmetic) Theorem 3 divides by Pi_{mu_i}.
   - Pi_{mu_i} vanishes at some framings (tau = -1/2 with even parts).
   - The statement is better written as [x^mu]F = sigma^(1-g) prod Pi Xi.
   - The proof of Cor 4 uses only integer tau >= 1, where Pi > 0, so nothing is affected.
6. (cosmetic) The proof of Theorem 2 calls R_mu a "homogeneous polynomial". For n = 1 it contains a c^(|mu|+1) hb^-1 term. The degree argument is unaffected.
7. (cosmetic) Theorem 3, step (a): the continuity remark about gamma = 0 is unnecessary, since gamma = 0 is an admissible value.
8. (cosmetic) "Genus-7 table computed here": engine_h7_r1.txt is textually identical to the pre-existing E-referee table ref_engine_h7_r1.txt. What it actually shows is that two independent engines agree.
9. (cosmetic) The summary calls Route E's M_irr(7) "partial". Brief s.2(g) records it as complete (194 points).
10. (cosmetic) R2's literature key [Eyn14a] for Eynard's CNTP paper differs from ABDKS's key [Eyn14]. In ABDKS, [Eyn14a] is Eynard's ICM 2014 survey.

## 6. What remains / not covered

- **A TR-free proof of Theorems 1 and 3.** A possible route is the closed n-point formulas for hypergeometric (Orlov-Scherbin) KP tau functions, which are an allowed input. This was not attempted here.
- **Closed multi-leaf forms** (SN, NR) were not extracted.
- **The extra zero at j = 2h in the Hurwitz limit** was not addressed by R2. It is contained in Route B's Theorem B, which is consistent with the top component of (6.1) (s.2.6).
- **Lemma A => Corollary R** (conifold all-genus regularity) is the main session's deduction and was not checked here.

## 7. Files (HVC\lemmaA\round2\R2_onepoint\referee2\)

Libraries:
- RefR2.pm: Bernoulli numbers; (6.1) as an exact polynomial in (sigma, j) and its components; engine-table parser.
- RefXY.pm: Theorem 1 at n = 1, 2 (own code).
- RefXY3.pm: Theorem 1 at n = 3 (own divided-difference code).

Scripts:
- ref_cmp61_engine.pl TABLE h: (6.1) and the gap formula vs an engine table, componentwise.
- ref_mirr_engine.pl and ref_mirr_engine_fast.pl TABLE h r: Lemma A, every component.
- ref_xy1.pl TAU GMAX JMAX; ref_xy2.pl TAU G JMAX [TABLE]; ref_xy3.pl TAU G JMAX [TABLE] [PERM].
- ref_string.pl TAU G JMAX [J3].
- ref_cj.pl MODE G TAU JMIN JMAX and ref_cj2.pl MODE G TAU JMAX [JMIN]: MV cut-and-join, with MODE = engine or formula.
- ref_limits.pl HMIN HMAX.
- ref_cmp_trhalf.pl FILE.

Subfolders:
- eng/: copies of RefHodge.pm, RefWK.pm, RefParse.pm and t_moments.pl. NEW tables ref_engine_h8_r1.txt, ref_engine_h8_r2.txt, ref_engine_h8_r3.txt and ref_engine_h9_r1.txt, with their logs.
- kp/: copies of Route A's a3_kp_nullity.pl, a3_kp_nullity2.pl and mvlib.pm. ref_hirota.pl is Route A's Hirota generator (sections 1-2 verbatim) with the values taken from the engine or from the R2 formulas. Outputs: kpn_g7_rerun.out, kpn2_g8_rerun.out, kpn2_g9_rerun.out, run_hirota_engine.out, run_hirota_formula.out.
- tr/: status_trhalf.pl, trhalf_g7_ref.out, trhalf_g9_ref.out.

Outputs: run_*.out, runs_referee2_foreground.txt, runs_referee2_background.txt.

Reproduce, for example:
```
REF_NCH=5 perl eng/t_moments.pl 9 1 eng/ref_engine_h9_r1.txt
perl ref_cmp61_engine.pl eng/ref_engine_h9_r1.txt 9
perl ref_cj.pl formula 10 1 1 32
MODE=formula perl kp/ref_hirota.pl 10 1 4 "17;18;19;20;21;22;23"
```

