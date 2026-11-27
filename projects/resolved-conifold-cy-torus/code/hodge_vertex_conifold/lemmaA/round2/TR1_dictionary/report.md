# Lemma A round 2: route TR1 (saved by the main session from the workflow journal, 2026-10-31)

Subagents cannot write report files, so this file holds the agent's structured return value: the fields first, then the full report text (report_text) verbatim.

## all_genus_claim

1

## all_genus_statement

P2 (Lemma B) holds for every genus g, every n >= 1 and every tau not in {0,-1}, with the brief's normalisation (proved modulo Zhou's BKMP theorem and the Eynard-Orantin symmetry theorem). Consequently, in every genus: Lemma A holds at sigma = 2 (tau = 1, -2) for all R, and at sigma = -1/4 for all R with odd parts. This is not Lemma A for all tau, and it does not prove (OL) or (6.1) in all genera.

## claims

### [1]

#### depends_on

J. Zhou arXiv:0910.4320 Thm 1 (BKMP for the framed one-leg vertex; also Chen arXiv:0910.3739, Eynard-Orantin CMP 337 (2015)); Eynard-Orantin CNTP 1 (2007) symmetry theorem; dilaton equation

#### evidence

complete proof in report s.3-4 (Lemmas 1-5, rationality in e and tau, limit e->oo); exact coefficientwise agreement with Route E's Hodge engine: g<=4,n<=3 at tau=1,2,3,-3,1/2,-1/2,2/3,-1/3,5/7,3/2; g<=5,n<=3 at tau=1,-1/2,2/3; g<=6,n<=3 at tau=-1/2; one leaf g<=6 at 7 tau; 0 failures

#### statement

THEOREM P2 (Lemma B): for every tau in C\{0,-1} and every (g,n) with 2g-2+n>0, omega^{C_tau}_{g,n} = sigma^(1-g) Sum_k <Omega prod psi^{k_i}>^{q_0}_{g,n} prod dxi_{k_i} (xi_0 = w/(1-tau w), xi_k = D^k xi_0); equivalently [x^mu]F_{g,n} = sigma^(1-g) prod Pi_{mu_i}(tau) Xi_g(mu;tau), i.e. Xi_h(mu;tau) = sigma^(h-1)[x^mu]F_{h,r}/prod C(mu_i(tau+1)-1,mu_i-1) whenever the denominator is non-zero; normalisation of [x^mu]F and of the TR kernel as in review_checks/referee1_trC1.pl

#### status

proved-modulo-literature

### [2]

#### depends_on

Zhou arXiv:0910.4320; W_0 = B - dx dx/(x-x)^2 vs B shown irrelevant via (LL)

#### evidence

(0,3),(1,1) by hand for all tau; exact for all (g,n) with g<=4, n<=3 at tau=2/3, 3, -2/5 vs Mathematica Hodge integrals (956 each); Zhou's (6) literally at tau=2/3, g<=2 (139)

#### statement

BKMP in our conventions (transcription of Zhou's Theorem 1, eqs. (3),(6),(8),(19)-(23)): omega^{y=log(1+w)}_{g,n} = sigma^(g-1+n) Sum_k <Omega prod psi^{k_i}>_{g,n} prod dxi_{k_i} (equivalently (-1)^n sigma^(g-1+n) for Zhou's y = -log(1+w)); i.e. A_g = sigma^(g-1), B = sigma, q_B = -psi; framing a = tau

#### status

proved-modulo-literature

### [3]

#### depends_on

same as P2; Route A Prop 2.1/Cor 2.2 for the x-form

#### evidence

exact vs Route A's q-tables at (tau,e) = (2/3,3/5),(2,-7/3),(1/2,4),(-1/2,1), g=2..4, n<=2 (168 each)

#### statement

Theorem F (interpolating family): for e != -1/sigma, omega^{log(1+w)-e/w}_{g,n} = sigma^(g-1+n)(1+e sigma)^(2-2g-n) Sum_k H_{g,n}(e sigma/(1+e sigma);k) prod dxi_{k_i}, H(s;k) = Sum_m (s^m/m!)<Omega prod psi^{k_i} T^m>; in x-form [x^mu]F = sigma^(g-1+n) prod Pi e^{-u|mu|} G_mu(u e^{-u}) at u = -e sigma (Route A's Lambert variable)

#### status

proved-modulo-literature

### [4]

#### depends_on

symmetry of TR correlators (EO 2007)

#### evidence

proof in report s.3

#### statement

Lemma 1 (loop equations): for a regular y on the one-branch-point curve (P^1, dX, B) the TR correlators satisfy omega = [E/omega~]_- (principal part), omega + sigma^*omega = O(zeta)dzeta, and E + Y sigma^*omega + Ybar omega = O(zeta^2)dzeta^2

#### status

proved

### [5]

#### depends_on

Lemma 1 (hence EO symmetry)

#### evidence

proof by induction via the identity (**) in report s.3; consistent with the exact family checks

#### statement

Lemma 2 (variational formula for a regular deformation y -> y + e eta at fixed X and B, eta holomorphic at the branch point, eta dX = Int_gamma Lambda B with gamma away from a): d/de omega_{g,n}(J) = Int_gamma Lambda(r) omega_{g,n+1}(J,r) for all stable (g,n)

#### status

proved

### [6]

#### evidence

direct computation (report s.3)

#### statement

Lemma 3: Y_0 dX = -Res_{r->0}B(.,r)/w_r + (1+tau) Int_{-1}^{0} B(.,r), and Int_gamma Lambda dxi_k = delta_{k0} - 1, i.e. the deformation direction is v = -psi/(1-psi) = lam_0 - lam_1; each half separately is an irregular deformation (pole at the branch point), only the sum -1/w is regular

#### status

proved

### [7]

#### evidence

dilaton equation; report s.3

#### statement

Lemma 5: Sum_K ((e sigma)^K/K!) <Omega prod psi^{k_i} v^K>_{g,n+K} = (1+e sigma)^(2-2g-n) H_{g,n}(e sigma/(1+e sigma);k) (dilaton resummation); both sides of Theorem F have their only pole at e = -1/sigma

#### status

proved

### [8]

#### evidence

Lagrange inversion; Zhou (71)-(76); hand computation with u_2 = tau(tau+2)/(3(tau+1)), u_3 = (u_2^2+kappa)/4; exact checks (0,2) m1+m2<=9 at 5 tau, (0,3),(1,1) at 5 more tau, genus 0 up to n=8

#### statement

Unstable and low cases: [x^mu] Int Y_0 dX = sigma Pi_mu/(mu(mu+1)); [x1^m1 x2^m2] log((w1-w2)/(x1-x2)) = sigma Pi_m1 Pi_m2/(m1+m2) (both = sigma^(1-0) prod Pi times Route A's Gamma(J)/Gamma(J-r+3)); omega^{Y_0}_{0,3} = sigma dxi_0^3, omega^{Y_0}_{1,1} = -((2+1/sigma)/24) dxi_0 + (1/24) dxi_1, omega^{log(1+w)}_{1,1} = -((1+sigma)/24) dxi_0 + (sigma/24) dxi_1

#### status

proved

### [9]

#### evidence

residue argument, report s.5

#### statement

A non-identity Moebius map s with s^*dX = dX exists iff tau in {1,-2,-1/2} (the residues 1, -(1+tau), tau of dX at 0, -1, oo are the CY weights); then s is w->1/w, w->-1-w, w->-w/(1+w) respectively

#### status

proved

### [10]

#### depends_on

Theorem 3.1 (TR statement, proof_M_hurwitz s.3, referee_hurwitz_1 s.1) + P2

#### evidence

uniform induction; exact TR zeros confirmed to genus 7 (one leaf) and genus 5 (n<=3)

#### statement

Lemma A at sigma = 2 (tau = 1, -2) in every genus: Xi_h(mu;1) = 0 for all h>=2, r>=1, mu in Z_{>=1}^r with |mu| <= 2h-3+r; in particular Sum_k (-1/2)^k Xi^(k)_h(R) = 0 on the irreducible range (proof_status 1.6(c), now unconditional)

#### status

proved-modulo-literature

### [11]

#### depends_on

tau = -1/2 variant of Theorem 3.1 + P2

#### statement

Lemma A at sigma = -1/4 (tau = -1/2) in every genus for every R with all parts odd, D(R) <= 2h-3 (proof_status 1.6(c), now unconditional; same as Route A's Theorem C)

#### status

proved-modulo-literature

### [12]

#### depends_on

P2

#### statement

proof_status 1.7(b) (Xi_3(2) = Xi_3(3) = 0 for all tau) and 1.7(c) (Lemma D at h distinct integer tau => M_irr(h)) are unconditional

#### status

proved-modulo-literature

### [13]

#### depends_on

P2

#### evidence

status_cmp61 runs (main session) and this report's cmp61_h2to7.txt

#### statement

(6.1) at sigma = 2 and sigma = -1/4 for h <= 6 (main session's exact TR runs, proof_status 1.8) is now an exact verification independent of the P2 hypothesis

#### status

computer-verified

### [14]

#### depends_on

P2 (exact TR = Hodge), lambda_g formula at sigma = 0, sigma-degree <= h of both sides

#### evidence

exact genus-7 one-leaf TR at tau = 1,2,3,1/2,-1/2,-1/3,2/3 (sigma = 2,6,12,3/4,-1/4,-2/9,10/9), 23 values of j each, 0 mismatches (runs/cmp61_h2to7.txt); Xi_7(j) = 0 exactly for j = 1..12

#### statement

(6.1) holds identically in j and tau for every h <= 7 (new at h = 7): sigma^h Xi_h(j) = prod_{l=1}^{2h-2}(j-l) [lambda^(2h)] S(e1 lambda)^(1+j e2) S(e2 lambda)^(1+j e1)

#### status

computer-verified


## key_formulas

- P2: omega^{C_tau}_{g,n} = sigma^(1-g) Sum_k X_{g,n}[k] prod dxi_{k_i}, with X[k] = <Omega prod psi^k>^{q_0}; [x^mu]F_{g,n} = sigma^(1-g) prod Pi_{mu_i}(tau) Xi_g(mu;tau), where Pi_mu = C(mu(tau+1)-1, mu-1)
- Laplace basis: xi_0 = w/(1-tau w) = D log(1+w) = Sum_mu Pi_mu x^mu, xi_k = D^k xi_0, D = x d/dx, [x^mu] xi_k = mu^k Pi_mu
- BKMP (Zhou) in our conventions: omega^{log(1+w)}_{g,n} = sigma^(g-1+n) Sum_k M_{g,n}[k] prod dxi_{k_i}, with M[k] = Int Omega prod psi^{k_i} (Zhou's y = -log(1+w) gives the factor (-1)^n)
- Theorem F: omega^{log(1+w) - e/w}_{g,n} = sigma^(g-1+n) (1+e sigma)^(2-2g-n) Sum_k H_{g,n}(e sigma/(1+e sigma); k) prod dxi_{k_i}, with H(s;k) = Sum_m (s^m/m!) <Omega prod psi^{k_i} T^m>, T = -psi^2/(1-psi)
- Dual cycle: Y_0 dX = -Res_{r->0} B(.,r)/w_r + (1+tau) Int_{-1}^{0} B(.,r); Int Lambda dxi_k = delta_{k0} - 1, so v(psi) = -psi/(1-psi)
- Lemma 5: Sum_K ((e sigma)^K/K!) <Omega prod psi^{k_i} v^K> = (1+e sigma)^(2-2g-n) H(e sigma/(1+e sigma); k)
- Variational identity (**): -Res K (Omega~/omega~) E = Res K [Omega sigma^*omega + omega sigma^*Omega], proved from Q = E + Y sigma^*omega + Ybar omega = O(zeta^2)
- x-form of F: [x^mu]F^{y_e} = sigma^(g-1+n) prod Pi_{mu_i} e^{-u|mu|} G_mu(u e^{-u}) at u = -e sigma, with G_mu(Y) = Sum_a (Y^a/a!) <Omega prod lam_{mu_i} lam_1^a>
- Unstable cases: [x^mu] Int Y_0 dX = sigma Pi_mu/(mu(mu+1)); [x1^m1 x2^m2] log((w1-w2)/(x1-x2)) = sigma Pi_m1 Pi_m2/(m1+m2)
- Low cases: omega^{Y_0}_{0,3} = sigma dxi_0^3; omega^{Y_0}_{1,1} = -((2+1/sigma)/24) dxi_0 + (1/24) dxi_1; omega^{log(1+w)}_{1,1} = -((1+sigma)/24) dxi_0 + (sigma/24) dxi_1
- Local coordinate: X(a+u) - X(a) = -c2 zeta^2, with c2 = tau^3/(2(1+tau)) and u u' = (tau^2/(1+tau)) zeta (a+u)(1+a+u), a = 1/tau
- Residues of dX at 0, -1, oo are 1, -(1+tau), tau (the CY weights); a Moebius s with s^*dX = dX and s != id exists iff tau is in {1, -2, -1/2}

## next_steps

- Lemma D for integer tau >= 2: the TR2 route. With P2 proved, Lemma D at h distinct values of tau implies Lemma A(h) unconditionally.
- Try to read Route A's criterion deg P_R <= N-1 (equivalent to Lemma A) as a statement about the e -> oo behaviour of TR on the family y = log(1+w) - e/w (Theorem F, u = -e sigma). This could give a TR proof of the degree bound at fixed tau.
- Extend (6.1) to h = 8: genus-8 one-leaf TR (run_xi.pl t 8 1 with SAVE) at 8 values of tau, compared with cmp61_xi.pl. Expect about 1 h per run; Stanage is suitable.
- Write out the DOSS/R-matrix route (Laplace transform of B on C_tau, Stirling expansion versus Mumford's formula for the triple Hodge class) as a BKMP-free second proof of P2.

## obstructions

- None for P2 itself: the proof is complete modulo two literature results, BKMP for the framed one-leg vertex (Zhou; Chen; Eynard-Orantin 2015) and the EO symmetry theorem. Their normalisations and symmetry were checked exactly to genus 4-6.
- P2 is a dictionary. It yields Lemma A only at values of tau where TR has a degree bound at the zero v = oo of x: tau = 1 for all parts, tau = -1/2 for odd parts. For (OL) at genus G it gives exactly one of the G-1 needed conditions (sigma = 2). tau = -1/2 gives nothing, because j = D+1 is even.
- The remaining TR-side input for Lemma A in general is Lemma D (degree bound for tau >= 2). The obstruction is unchanged: no Moebius symmetry with s^*dX = dX exists for tau not in {1,-2,-1/2} (now proved), and for integer tau >= 2 the correlators have low degree at the other zero v = 0 of x.
- Theorem F (the TR family y = log(1+w) - e/w) works at fixed tau and supplies no new values of sigma.

## requested_checks

- Referee: check Lemma 1 (loop equations from the recursion plus symmetry), Lemma 2 (the variational formula via the identity (**)), and Step 2 of s.4 (rationality in tau).
- Referee: check the transcription of Zhou's Theorem 1 (arXiv:0910.4320v2, eqs. (3), (6)-(8), (19)-(23)) into our conventions: y_Zhou = 1/(1+w), TR y = -log(1+w), a = tau, phi_b = xi_b, <prod tau T_g> = (-sigma)^g M.
- Main session (Mathematica, optional): run hvc_tr_tau.wls at a non-integer tau (e.g. 2/3) with 'hodge' through genus 3-4. This gives an independent implementation to cross-check TRxi.pm.
- Optional (Stanage, long): `SAVE=s_t.txt perl run_xi.pl t 8 1 base` for 8 values of t (e.g. 1, 2, 3, 4, 1/2, -1/2, 2/3, -1/3), then `perl cmp61_xi.pl t s_t.txt`. Together with sigma = 0 this would prove (6.1) at h = 8.

## route

TR1: the dictionary P2 (Lemma B) via BKMP for the framed one-leg C^3 vertex (Zhou's Theorem 1, transcribed and checked) + the Eynard-Orantin variational formula along the REGULAR family y_e = log(1+w) - e/w (proved here from the loop equations) + rationality in e + homogeneity/dilaton as e -> oo + rationality in tau; exact Perl/BigRat TR on C_tau for rational tau (two independent implementations) compared coefficient by coefficient with Route E's Hodge engine and the Mathematica Hodge tables.

## summary

P2 (Lemma B) is now a theorem for every tau in C\{0,-1} and every stable (g,n), in particular every h >= 2 and r >= 1, with the normalisation of the brief: omega^{C_tau}_{g,n} = sigma^(1-g) Sum_k <Omega prod psi^{k_i}>^{q_0}_{g,n} prod dxi_{k_i}, i.e. [x^mu]F_{g,n} = sigma^(1-g) prod Pi_{mu_i}(tau) Xi_g(mu;tau), and Xi_h(mu;tau) = sigma^(h-1)[x^mu]F_{h,r}/prod C(mu_i(tau+1)-1,mu_i-1) whenever the denominator is non-zero. The proof follows referee_hurwitz_1.md s.4 and is written out completely. Its only external inputs are J. Zhou's BKMP theorem for the framed one-leg vertex (arXiv:0910.4320; also proved by Chen and by Eynard-Orantin, CMP 337), the Eynard-Orantin symmetry theorem, and the dilaton equation. Zhou's statement is transcribed into our conventions (TR y = -log(1+w), A_g = sigma^(g-1), B = -sigma) and checked exactly. Everything else is proved: the loop equations, the variational formula for the regular deformation (the separate lam_0 and lam_1 deformations are irregular, as the referee warned, but their sum is -1/w), the dual cycle with pairing v = -psi/(1-psi), rationality in e with the single pole e = -1/sigma on both sides, the limit e -> oo that produces the n-independent factor sigma^(1-g), and the extension from integer to all tau. A byproduct is Theorem F: TR on y = log(1+w) - e/w gives the whole family between the Marino-Vafa point (e = 0) and the dilaton point (e = oo). Its parameter is Route A's Lambert variable u = -e sigma. The unstable cases (0,1) and (0,2) fit the same formula with Route A's genus-0 values Gamma(J)/Gamma(J-r+3).

About 25,000 exact coefficient comparisons were made, with no failure. They cover P2 for g <= 4, n <= 3 at ten values of tau including non-integer ones; g <= 5 at three values; g <= 6, n <= 3 at tau = -1/2; Zhou's theorem to genus 4 at three values; Theorem F at four (tau, e) pairs; and genus 0 up to n = 8. With P2 proved, every "mod P2" result of proof_status s.1.6-1.8 becomes unconditional, uniformly in the genus. Lemma A holds at sigma = 2 (tau = 1, -2) for every h and every R. It holds at sigma = -1/4 for every R with odd parts. Also unconditional: Xi_3(2) = Xi_3(3) = 0 for all tau, the implication "Lemma D at h integer tau => Lemma A(h)", and the TR checks of (6.1). New result (computer-assisted, via P2): (6.1) holds identically in (j, tau) for every h <= 7, from exact genus-7 TR at seven non-zero sigma plus the lambda_g limit. Where it stops: P2 is a dictionary and gives Lemma A only where TR has a degree bound (tau = 1; tau = -1/2 for odd parts). A non-trivial Moebius symmetry with s^*dX = dX exists iff tau is in {1,-2,-1/2}; this is now proved. For (OL), sigma = 2 is one of the G-1 needed conditions, and the remaining input is Lemma D for tau >= 2, which this work does not touch.

## report_text

# Route TR1 (round 2): the dictionary P2 (Lemma B): proof, exact verification, consequences

Folder: HVC\lemmaA\round2\TR1_dictionary (HVC = G:\My Drive\Projects\Old\refGW\hodge_vertex_conifold). All arithmetic is exact (Perl 5.38, Math::BigRat); no wolframscript, no python. TR is used only for the Hodge CohFT (the framed C^3 vertex, i.e. the curve C_tau), as allowed by brief s.0. Nothing below concerns the resolved conifold x C^2.

## 0. Outcome

1. **P2 is a theorem.** It holds for every tau in C \ {0,-1} and every (g,n) with 2g-2+n > 0, so for every h >= 2 and r >= 1. It also holds for g = 0, 1 and for the unstable (0,1), (0,2) with the standard conventions. The normalisation is that of the brief (s.1.3).
   - Status: proved modulo two published results. (i) BKMP for the framed one-leg C^3 vertex in J. Zhou's form (arXiv:0910.4320v2, Theorem 1 with eqs. (3), (6), (8), (19)-(23)); it is also proved by L. Chen (arXiv:0910.3739) and by Eynard-Orantin (CMP 337 (2015)). It is transcribed into our conventions in s.2.1 and checked exactly. (ii) The Eynard-Orantin symmetry theorem (CNTP 1 (2007), s.4). Plus the dilaton equation.
   - Everything else is proved here: the loop equations (Lemma 1), the variational formula for the regular deformation (Lemma 2), the dual cycle (Lemma 3), rationality in e and the limit e -> oo (Lemma 4), the Hodge side of the family (Lemma 5), and the extension to all tau (s.4).
   - This is the route of referee_hurwitz_1.md s.4, written out.
2. **Theorem F (byproduct).** TR on y_e = log(1+w) - e/w interpolates between the Marino-Vafa point (e = 0, Zhou) and the dilaton point (e -> oo, P2):
       omega^{y_e}_{g,n} = sigma^(g-1+n) (1+e sigma)^(2-2g-n) Sum_k H_{g,n}(e sigma/(1+e sigma); k) prod dxi_{k_i}.
   Its parameter is Route A's Lambert variable, u = -e sigma.
3. **Exact verification (s.6).** About 25,000 coefficient comparisons, 0 failures, 0 symmetry violations:
   - P2 against Route E's Hodge engine (no TR) for all (g,n) with g <= 4, n <= 3 at tau = 1, 2, 3, -3, 1/2, -1/2, 2/3, -1/3, 5/7, 3/2;
   - g <= 5, n <= 3 at tau = 1, -1/2, 2/3; g <= 6, n <= 3 at tau = -1/2;
   - Zhou's theorem against Hodge integrals at tau = 2/3, 3, -2/5 (g <= 4, n <= 3);
   - Theorem F at four (tau, e);
   - genus 0 up to n = 8, and the (0,2) case.
   Two independent implementations agree, and they reproduce every number of the referee's tau = 1 code.
4. **Consequences (s.7).** Every "mod P2" item of proof_status_mixed_scale.md s.1.6-1.8 is now unconditional:
   - Lemma A at sigma = 2 (tau = 1, -2) for every h and every R (Theorem 3.1);
   - Lemma A at sigma = -1/4 for every R with odd parts;
   - Xi_3(2) = Xi_3(3) = 0 for all tau;
   - "Lemma D at h integer tau => Lemma A(h)";
   - (6.1) at sigma = 2 and sigma = -1/4 for h <= 6.
   For (OL) this gives, in every genus, the sigma = 2 condition: one of the G-1 needed conditions.
5. **New (computer-assisted, via P2):** (6.1) holds identically in (j, tau) for every h <= 7. The evidence is exact genus-7 TR at seven non-zero sigma, plus sigma = 0.
6. **Remark (proved):** a non-trivial Moebius s with s^*dX = dX exists iff tau is in {1,-2,-1/2}, because the residues of dX are the CY weights. So the involution of Theorem 3.1 has no analogue in Aut(P^1) at other tau. P2 does not touch the remaining obstruction (Lemma D, s.8).

## 1. Conventions and statements

### 1.1 Hodge side (c = 1)
- Weights (1, tau, -1-tau), sigma = tau(1+tau) != 0.
- Omega = Lambda^vee(1) Lambda^vee(tau) Lambda^vee(-1-tau)/(-sigma)^g, with Lambda^vee(u) = Sum_i (-1)^i lambda_i u^(g-i). It is pulled back along forgetful maps.
- Dilaton leaf: T(psi) = -psi^2/(1-psi).

For k in Z_{>=0}^n:

    H_{g,n}(s;k) := Sum_{m>=0} (s^m/m!) Int_{Mbar_{g,n+m}} Omega prod_i psi_i^{k_i} prod_{j<=m} T(psi_{n+j})

This is a polynomial in s of degree <= 3g-3+n. Define:
- M_{g,n}[k] := H(0;k), the Marino-Vafa integrals in descendant form (no dilaton leaves);
- X_{g,n}[k] := H(1;k), the dilaton point q_0 = lam_0 - lam_1.

Then Xi_g(mu) = Sum_k X[k] prod mu_i^{k_i}, and Route A's Q_g(s;R) = Sum_k H(-s;k) prod j_i^{k_i}.

Data sources:
- X: Route E's engine tables E_formulas\engine_h<h>_r<r>.txt (Mumford GRR + DVV, no TR), via X[k] = Sum_c (-1/sigma)^c [j^k] Xi^(c).
- M and H: the Mathematica tables lemmaA\data\mv_g<g>_l<l>.txt and q_g<g>_r<r>.txt. These use weights (c,-x,-1) and are multiplied by c^(3g-3+2r), with c = 1/(1+tau), x = -tau/(1+tau).
- Genus 0: X_{0,n}[k] = [j^k] prod_{i=1}^{n-3} (J-i) (Route A's Gamma(J)/Gamma(J-r+3)), and M_{0,n}[k] = (n-3)!/prod k_i! for |k| = n-3.

### 1.2 TR side

**Curve.** P^1 with coordinate w.
- x = w (1+w)^(-tau-1), a power series in w for any tau; X = log x.
- dX = (1 - tau w) dw/(w(1+w)). It has simple poles at 0, -1, oo with residues 1, -(1+tau), tau (the CY weights), and one simple zero a = 1/tau.
- B(p,q) = dw_p dw_q/(w_p - w_q)^2.

**Local coordinate at a.** On a disc U around a define zeta by X(a+u) - X(a) = -c2 zeta^2, with c2 = tau^3/(2(1+tau)) and u = zeta + O(zeta^2). Equivalently u u' = kappa zeta (a+u)(1+a+u), kappa = tau^2/(1+tau). The local involution is sigma: zeta -> -zeta.

**Correlators.** Call y holomorphic on U with dy(a) != 0 "regular". Its correlators (2g-2+n > 0) are given by omega_{0,2} = B and

    omega_{g,n+1}(p0,J) = Res_{q->a} K^y(p0,q) E_{g,n+1}(q;J),
    K^y(p0,q) = (1/2) Int_{sigma q}^{q} B(p0,.) / ((y(q) - y(sigma q)) dX(q)),
    E_{g,n+1}(q;J) = omega_{g-1,n+2}(q,sigma q,J) + Sum'_{h,I} omega_{h,1+|I|}(q,I) omega_{g-h,1+|J\I|}(sigma q,J\I).

Sum' omits (h,I) = (0,empty) and (g,J); omega_{0,2} factors are included. This is the convention of review_checks/referee1_trC1.pl and of Zhou's (21)-(23), and it gives omega_{0,3} = -Res BBB/(dX dy). Three choices of y appear:
- Y_0 := -1/w, the base curve C_tau of P2;
- y_B := log(1+w);
- y_e := y_B + e Y_0.

**Laplace basis.** xi_0 := w/(1 - tau w), xi_k := D^k xi_0, with D = d/dX = x d/dx.
- By Lagrange inversion, xi_0 = D log(1+w) = Sum_{mu>=1} Pi_mu(tau) x^mu with Pi_mu = C(mu(tau+1)-1, mu-1). Hence [x^mu] xi_k = mu^k Pi_mu(tau).
- The xi_k are rational and vanish at w = 0. xi_0(-1) = -1/(1+tau), and xi_k(-1) = 0 for k >= 1 (D vanishes at -1).
- xi_k has one pole, at a, of order 2k+1, with odd principal part in zeta.
- Let V be the space of meromorphic differentials on P^1 with poles only at a and even principal part in zeta. Then {dxi_k} is a basis of V: P^1 has no holomorphic differentials, and the principal parts are triangular.

**Expansion at x = 0.** [x^mu]F_{g,n} (mu_i >= 1) is the x^mu coefficient of a primitive of omega_{g,n} expanded at w_i = 0, as in referee1_trC1.pl. If omega_{g,n} = Sum_k C[k] prod dxi_{k_i}, then [x^mu]F_{g,n} = Sum_k C[k] prod mu_i^{k_i} Pi_{mu_i}.

### 1.3 Statements

**THEOREM P2.** For tau in C\{0,-1} and 2g-2+n > 0:

    omega^{Y_0}_{g,n} = sigma^(1-g) Sum_k X_{g,n}[k] prod_i dxi_{k_i}(w_i).

Equivalently, [x^mu]F_{g,n} = sigma^(1-g) prod_i Pi_{mu_i}(tau) Xi_g(mu;tau) for all mu in Z_{>=1}^n. Hence

    Xi_h(mu;tau) = sigma^(h-1) [x^mu]F^{C_tau}_{h,r} / prod_i C(mu_i(tau+1)-1, mu_i-1)

whenever the denominator is non-zero. The denominator vanishes iff tau = -b/mu_i with 1 <= b < mu_i (for example tau = -1/2 with mu_i even); the coefficient form still determines Xi there.

**THEOREM F.** For tau in C\{0,-1}, e != -1/sigma and 2g-2+n > 0:

    omega^{y_B+eY_0}_{g,n} = sigma^(g-1+n) (1+e sigma)^(2-2g-n) Sum_k H_{g,n}(e sigma/(1+e sigma); k) prod_i dxi_{k_i}.

- At e = 0 this is BKMP (s.2.1); as e -> oo it gives P2 (s.4).
- With Route A's Prop. 2.1 / Cor. 2.2 (u = s/(1+s), here s = -e sigma/(1+e sigma), so u = -e sigma):
      [x^mu]F^{y_e}_{g,n} = sigma^(g-1+n) prod Pi_{mu_i} e^{-u|mu|} G_mu(u e^{-u}) at u = -e sigma, with G_mu(Y) = Sum_a (Y^a/a!) <Omega prod lam_{mu_i} lam_1^a>_{g,n+a}.
- So the deformation by -e/w adds unit legs (a p_1-shift of the MV vertex), and Route A's dilaton point u = oo is e = oo.

## 2. Inputs from the literature

### 2.1 BKMP for the framed one-leg vertex (J. Zhou, arXiv:0910.4320v2), transcribed

**Zhou's definitions (s.2-3).** For a partition mu of length l and framing a:

    W_{g,mu}(a) = ((-1)^(g+l)/|Aut mu|) (a(a+1))^(l-1) prod_i Pi_{mu_i}(a) Int Lambda^vee_g(1) Lambda^vee_g(-a-1) Lambda^vee_g(a) / prod(1 - mu_i psi_i).

By eqs. (3)-(7), for 2g-2+n > 0:

    W_g(x_1..x_n;a) := Sum_{l(mu)=n} z_mu W_{g,mu} m_mu prod dx_i/x_i = (-1)^(g+n) (a(a+1))^(n-1) Sum_b <prod tau_{b_i} T_g(a)>_g prod_i dphi_{b_i}(x_i;a),

where phi_b(x;a) = Sum_m Pi_m(a) m^b x^m and <prod tau_b T_g(a)>_g = Int prod psi_i^{b_i} Lambda^vee_g(1) Lambda^vee_g(a) Lambda^vee_g(-1-a).
- Curve (8): x = y^a - y^(a+1); near (0,1), y = 1-u and x = u(1-u)^a.
- Initial data: (19) W_0(y;a) = ln y dx/x; (20) W_0(y1,y2) = dy1 dy2/(y1-y2)^2 - dx1 dx2/(x1-x2)^2.
- Recursion (21): W_g(y_[n]) = Res_{z=0} (dE_z(y_1)/omega(z)) (W_{g-1}(q, q~, y_[n]_1) + Sum W_{g1}(q,y_A) W_{g2}(q~,y_B)), with the (0,1) terms set to 0.
- Kernel: (22) omega(z) = (ln y(q) - ln y(q~)) dx(q)/x(q); (23) dE_z(y_1) = (1/2)(1/(y_1 - y(q)) - 1/(y_1 - y(q~))) dy_1.
- Theorem 1 (the Bouchard-Marino conjecture for the framed one-leg vertex): the W_g of (3) satisfy (19)-(21).

**Transcription.**
- Set y_Zhou = 1/(1+w) (so u = w/(1+w)) and a = tau. Then x = u(1-u)^a = w/(1+w)^(tau+1) is our x, and the ramification point y = a/(a+1) is w = a.
- Zhou's TR "y" is y_Z := ln y_Zhou = -log(1+w), and his dx/x is dX.
- dE_z = (1/2) Int_{q~}^{q} B, because B is Moebius invariant. omega(z) = (y_Z(q) - y_Z(q~)) dX(q). So Zhou's recursion is ours with y = y_Z (same kernel, same order of q and q~).
- phi_b = xi_b.
- Lambda^vee(1) Lambda^vee(a) Lambda^vee(-1-a) = (-sigma)^g Omega, so <prod tau_b T_g(a)> = (-sigma)^g M_{g,n}[b].
- Zhou's (0,2) term W_0 = B - S, with S(q,p) = dx(q)dx(p)/(x(q)-x(p))^2, gives the same stable correlators as B. Indeed S is sigma-invariant in q and O(zeta) dzeta (dx(a) = 0), so the extra terms are S(omega + sigma^*omega) = O(zeta^2) dzeta^2 by Lemma 1(b) (respectively S(B + sigma^*B), and S S at (0,3)), and K = O(zeta^-1) dzeta^-1 kills them.

Hence Zhou's Theorem 1 reads

    omega^{y_Z}_{g,n} = (-1)^(g+n) sigma^(n-1) (-sigma)^g Sum_k M_{g,n}[k] prod dxi_{k_i} = (-1)^n sigma^(g-1+n) Sum_k M_{g,n}[k] prod dxi_{k_i},

and homogeneity (omega^{lambda y} = lambda^(2-2g-n) omega^y with lambda = -1) gives

    (BKMP)  omega^{y_B}_{g,n} = sigma^(g-1+n) Sum_k M_{g,n}[k] prod_i dxi_{k_i}   for 2g-2+n > 0.

**The factors A_g B^n.** In the referee's notation q_B = -psi (kappa = 1), A_g = sigma^(g-1), and B = sigma (B = -sigma for Zhou's y_Z). Either way A_g B^(2-2g) = sigma^(1-g).

**Framing and signs.** Zhou's framing a is our tau: the CY weights are (1, tau, -1-tau) and the leaves are Pi_mu(tau). BKMP's own curve x + y^a + y^(a+1) = 0 (Zhou's (14)) differs by x -> -(-1)^a x, y -> -y and is not used.

**Range of a.** We use Theorem 1 only for a in Z_{>=1}, where (8) is a plane algebraic curve. s.4 extends P2 to all tau. Zhou's argument is in fact local at the ramification point, with coefficients rational in a, and our exact checks at a = 2/3 and -2/5 agree with that.

**Checks.**
- (0,3) and (1,1) by hand for all tau (s.5).
- (BKMP) exactly for all (g,n) with g <= 4, n <= 3 at tau = 2/3, 3, -2/5, against Hodge integrals computed without TR or MV: 956 coefficients each, 0 failures.
- Zhou's (6) literally (y = y_Z, sign (-1)^n) at tau = 2/3 for g <= 2, n <= 3: 139/139.

### 2.2 Symmetry
omega^y_{g,n} is symmetric in all its arguments [EO07, s.4]. The recursion is symmetric in J by construction; the content of the theorem is the symmetry between p0 and J.
- EO's proof uses only the recursion near the branch points and the properties of B (symmetric, double pole without residue on the diagonal). It applies verbatim to TR with respect to the meromorphic form dX on P^1.
- Symmetry also follows from DOSS's identification of TR with Givental's formula for local spectral curves (CMP 328 (2014)).
- Every correlator computed in s.6 is exactly symmetric.

### 2.3 Dilaton equation
For a class pulled back along the forgetful maps (Omega is, since pi^* lambda_i = lambda_i) and 2g-2+N > 0:

    Int_{Mbar_{g,N+1}} Omega psi_{N+1} prod_{i<=N} psi_i^{k_i} = (2g-2+N) Int_{Mbar_{g,N}} Omega prod psi_i^{k_i}.

## 3. Lemmas (complete proofs)

**Notation.** Fix a regular y and J away from a. Write
- omega := omega^y_{g,n+1}(.,J) and E := E_{g,n+1}(.;J);
- Y := y dX, Ybar := sigma^*Y, omega~ := Y - Ybar = (y - y o sigma) dX;
- Phi := E/omega~.

Orders at a: dX = -2c2 zeta dzeta, omega~ = c zeta^2 (1 + O(zeta)) dzeta with c != 0, and K = (1/2)(Int_{-zeta}^{zeta} B(p0,.))/omega~ = O(zeta^-1) dzeta^-1.

**Lemma 1 (loop equations).** For 2g-2+n+1 > 0:
- (a) omega = [Phi]_-, the principal part of Phi at a. Phi is even in zeta, so Res_a Phi = 0 and omega is holomorphic at w = oo.
- (b) (LL) omega + sigma^*omega = O(zeta) dzeta.
- (c) (QL) Q := E + Y sigma^*omega + Ybar omega = O(zeta^2) dzeta^2.

*Proof.*
1. By 2.2 and the structure of the sum, sigma^*E = E (swap (h,I) <-> (g-h, J\I)). Also sigma^*omega~ = -omega~. Hence sigma^*Phi = -Phi, i.e. Phi = phi(zeta) dzeta with phi even.
2. Int_{sigma q}^{q} B(p0,.) = dw_{p0}(1/(w_{p0} - w_q) - 1/(w_{p0} - w_{sigma q})). Substitute q -> sigma q in the second term (residues are invariant under sigma, which fixes a). This gives omega(p0) = Res_{q->a} dw_{p0} Phi(q)/(w_{p0} - w_q) = [Phi]_-(p0). So h := omega - Phi is holomorphic near a.
3. (b): omega + sigma^*omega = h + sigma^*h is holomorphic and sigma-invariant. An invariant f(zeta) dzeta has f odd, so it is O(zeta) dzeta.
4. (c): E = Phi omega~ = (omega - h)(Y - Ybar), so Q = Y(omega + sigma^*omega) - h omega~ = O(zeta)O(zeta) + O(1)O(zeta^2). QED

**Lemma 2 (variational formula for a regular deformation at fixed X and B).**
Hypotheses:
- eta is holomorphic on U, and Omega := eta dX is meromorphic on P^1 with no pole at a.
- Int_gamma Lambda(r)(.) is a finite combination of residues Res_{r->b} Lambda_b(r)(.) (b not in Ubar, Lambda_b meromorphic) and path integrals c Int_{gamma_j}(.) (gamma_j in P^1\Ubar), with Int_gamma Lambda(r) B(p,r) = Omega(p) for p in U.
- y_e := y + e eta is regular except at one point e_*.

Then for every stable (g,n) and e != e_*:

    (V)  d/de omega^{y_e}_{g,n}(J) = Int_gamma Lambda(r) omega^{y_e}_{g,n+1}(J,r).

*Proof (induction on 2g-2+n).*

(i) Genus 0. d/de omega_{0,2} = 0, since B is fixed. Also Int Lambda(r) omega_{0,3}(z1,z2,r) = 0: by the recursion and Fubini (gamma is away from U) it equals

    Res_q K(z1,q) dX(q) [eta(sigma q) B(q,z2) + eta(q) sigma^*B(.,z2)(q)].

Writing B(q,z2) = b(zeta) dzeta, the bracket is (eta(-zeta) b(zeta) - eta(zeta) b(-zeta)) dzeta = O(zeta) dzeta, so the residue is 0. Both identities extend to z near a by analytic continuation.

(ii) Inductive step. Differentiate the recursion:

    d/de omega_{g,n+1} = Res K' E + Res K E',   with K' = -K Omega~/omega~ and Omega~ := Omega - sigma^*Omega.

Every factor of E is stable with smaller 2h-2+m, or is omega_{0,2}. By the induction hypothesis and (i), its e-derivative is Int Lambda(r) of the same factor with r appended (slot order is irrelevant by symmetry).

E_{g,n+2}(q;J,r) consists of exactly these terms plus two new ones, B(q,r) omega(sigma q) and omega(q) B(sigma q,r), where r sits alone with q in an omega_{0,2}. Using Int Lambda(r) B(q,r) = Omega(q) for q in U:

    Int Lambda(r) omega_{g,n+2}(p0,J,r) = Res K E' + Res K [Omega sigma^*omega + omega sigma^*Omega].

So (V) for (g,n+1) is equivalent to

    (**)  -Res K (Omega~/omega~) E = Res K [Omega sigma^*omega + omega sigma^*Omega].

Use the identity A Bbar + Abar B = (1/2)(A+Abar)(B+Bbar) - (1/2)(A-Abar)(B-Bbar) and Lemma 1(c):
- E = Q - (1/2)(Y+Ybar)(omega + sigma^*omega) + (1/2) omega~ (omega - sigma^*omega);
- Omega sigma^*omega + omega sigma^*Omega = (1/2)(Omega + sigma^*Omega)(omega + sigma^*omega) - (1/2) Omega~ (omega - sigma^*omega).

The terms -(1/2) Res K Omega~ (omega - sigma^*omega) agree on both sides. All remaining terms vanish:
- Omega~/omega~ = (eta(zeta) - eta(-zeta))/(y_e(zeta) - y_e(-zeta)) = O(1);
- Q = O(zeta^2) dzeta^2;
- Y + Ybar and Omega + sigma^*Omega are O(zeta) dzeta;
- omega + sigma^*omega = O(zeta) dzeta by Lemma 1(b).

So each remaining integrand is O(zeta) dzeta. The base cases (0,3) and (1,1) are covered by the same computation: their E consists of omega_{0,2} factors only, and Lemma 1 holds for them. QED

This is EO07 s.5 ("variations of the curve") in the regular case. A proof is given because TR with respect to log x with non-algebraic y is not literally the setting of EO07.

**Lemma 3 (dual cycle of the base curve).**

    Y_0 dX = -dw/w^2 + (1+tau)(dw/w - dw/(1+w)) = Int_gamma Lambda B,
    Int_gamma Lambda(r)(.) := -Res_{r->0} (1/w_r)(.) + (1+tau) Int_{r=-1}^{0} (.)   (path avoiding Ubar).

For every k, Int_gamma Lambda dxi_k = delta_{k0} - 1 =: v_k. So v(psi) := Sum_k v_k psi^k = -psi/(1-psi) = lam_0 - lam_1, the dilaton direction.

*Proof.* Res_{r=0} dr/(r(p-r)^2) = 1/p^2, and Int_{-1}^{0} dr/(p-r)^2 = 1/p - 1/(p+1). Next, Res_{r->0} dxi_k(r)/w_r = [w^1] xi_k = [x^1] xi_k = 1, because xi_k(0) = 0 and x = w + O(w^2). Finally Int_{-1}^{0} dxi_k = xi_k(0) - xi_k(-1) = delta_{k0}/(1+tau). QED

The two halves separately are irregular. They are eta_1 dX with eta_1 = -(1+w)/(w(1 - tau w)), and eta_0 dX with eta_0 = (1+tau)/(1 - tau w), and both have a pole at a. Only their sum eta_0 + eta_1 = -1/w is regular. This is the referee's point; the combined deformation is the one used.

**Lemma 4 (rationality in e, and the limit).** By Lemma 1(a) and symmetry, omega^y = Sum_k C^y[k] prod dxi_{k_i}. For y_e = y_B + e Y_0:
- C^{y_e}[k] lies in Q(tau)[e, 1/(1+e sigma)];
- lim_{e->oo} e^(2g-2+n) C^{y_e}[k] = C^{Y_0}[k].

*Proof.* The recursion uses y only through 1/(y(zeta) - y(-zeta)) = 1/(zeta Yh(zeta)), where Yh = Yh_B + e Yh_0 is even. Its constant term is Yh(0) = 2(y_B'(0) + e Y_0'(0)) = 2(tau/(1+tau))(1 + e sigma), since y_B'(0) = 1/(1+a) and Y_0'(0) = 1/a^2 = tau^2. Each residue needs only finitely many coefficients of 1/Yh, and these are polynomials in e over powers of (1+e sigma). For the limit, homogeneity gives C^{y_e} = e^(2-2g-n) C^{Y_0 + y_B/e}, and C^{Y_0 + delta y_B} is regular at delta = 0 because Y_0'(0) != 0. QED

So the only pole of the TR side is at e = -1/sigma, where dy_e(a) = 0.

**Lemma 5 (Hodge side of the family).** As formal power series in e, for 2g-2+n > 0:

    Sum_{K>=0} ((e sigma)^K/K!) <Omega prod psi_i^{k_i} prod_{j<=K} v(psi_{n+j})>_{g,n+K} = (1+e sigma)^(2-2g-n) H_{g,n}(e sigma/(1+e sigma); k).

*Proof.* Write v = -psi + T. For j leaves T and K-j leaves -psi, the dilaton equation removes the -psi's with the factor (-1)^(K-j) (N_j)(N_j+1)...(N_j+K-j-1), where N_j = 2g-2+n+j. Summing over K at fixed j gives Sum_l ((-e sigma)^l/l!) (N_j)_l = (1+e sigma)^(-N_j). QED

The right side is rational in e, and its only pole is also at e = -1/sigma: there the dilaton coefficient of q_B + e sigma v = -(1+e sigma) psi - e sigma psi^2/(1-psi) vanishes.

## 4. Proofs of Theorems F and P2

**Step 1 (tau in Z_{>=1}).**
1. Apply Lemma 2 with y = y_B, eta = Y_0 (regular at a, since a != 0) and the cycle of Lemma 3. Here e_* = -1/sigma.
2. Iterate (V); this is allowed near e = 0, and d/de commutes with the residues and integrals because everything is rational in e. The K-th e-derivative at 0 of C^{y_e}[k] is then Sum_l C^{y_B}_{g,n+K}[k,l] prod v_{l_j}.
3. By (BKMP) this equals sigma^(g-1+n) sigma^K <Omega prod psi^{k_i} prod v(psi_{n+j})>_{g,n+K}.
4. By Lemma 5 these are the Taylor coefficients at 0 of R(e) := sigma^(g-1+n) (1+e sigma)^(2-2g-n) H(e sigma/(1+e sigma); k).
5. Both C^{y_e}[k] and R(e) are rational in e and regular at 0 (Lemma 4; H is a polynomial). Two rational functions with the same Taylor series at a regular point coincide. This is Theorem F at tau.
6. By Lemma 4, C^{Y_0}[k] = lim_{e->oo} e^(2g-2+n) R(e) = sigma^(g-1+n) sigma^(2-2g-n) H(1;k) = sigma^(1-g) X[k]. This is P2 at tau.

The factor sigma^(1-g) is independent of n because the TR scaling exponent 2-2g-n equals the dilaton exponent: B^n B^(2-2g-n) = B^(2-2g).

**Step 2 (all tau).**
- C^{Y_0}_{g,n}[k](tau) is produced by finitely many field operations from Taylor coefficients in Q[tau, 1/tau, 1/(1+tau)]. The divisions are only by integers, c2, a, 1+a, Yh(0) = 2 tau^2, and the diagonal entries T_kk of the basis change, which are monomials in tau, (1+tau) and c2.
- So C^{Y_0}[k] lies in Q[tau, 1/tau, 1/(1+tau)], and specialising to any tau not in {0,-1} commutes with the recursion.
- The right side sigma^(1-g) X[k] is a polynomial in 1/sigma.
- The two sides agree on Z_{>=1} (Step 1), hence identically, hence at every tau not in {0,-1}.
- Theorem F extends in the same way (rational in (tau, e)). QED

**Uniformity in the genus.** No step uses a genus-dependent input:
- Lemma 1 holds for each (g,n);
- Lemma 2 is a single induction over all (g,n);
- Zhou's theorem is an all-(g,n) statement;
- Lemma 5 and s.4 are per-(g,n) arguments.

## 5. Low and unstable cases; remarks

**(0,3) and (1,1), by hand for all tau.** From the ODE, u_2 = tau(tau+2)/(3(tau+1)) and u_3 = (u_2^2 + kappa)/4.
- omega_{0,3} = (1/(2 c2 y'(0))) e_0^3 with e_0 = tau^2 dxi_0. So omega^{Y_0}_{0,3} = sigma dxi_0^3 and omega^{y_B}_{0,3} = sigma^2 dxi_0^3, matching X_{0,3} = M_{0,3} = 1.
- For (1,1): B(q, sigma q) = -(1/4) zeta^-2 (1 + (kappa - 3 u_2^2) zeta^2 + ...) dzeta^2. The basis change is T_00 = 1/tau^2, T_11 = 1/(2 c2 tau^2), T_10 = (u_3 - u_2^2)/(2 c2 tau^2). These give
      omega^{y_B}_{1,1} = -((1+sigma)/24) dxi_0 + (sigma/24) dxi_1 = sigma (M_{1,1}[0] dxi_0 + M_{1,1}[1] dxi_1),
      omega^{Y_0}_{1,1} = -((2+1/sigma)/24) dxi_0 + (1/24) dxi_1 = X_{1,1}[0] dxi_0 + X_{1,1}[1] dxi_1.
  Here M_{1,1}[0] = -(1+1/sigma)/24, X_{1,1}[0] = -(2+1/sigma)/24 (one dilaton leaf, <tau_0 tau_2>_1 = 1/24), and M[1] = X[1] = 1/24.
- These closed forms were also checked exactly at tau = 7/3, -5/2, 11, -1/7, 13/5.

**Genus 0.** X_{0,n}[k] = [j^k] prod_{i=1}^{n-3}(J-i) follows from P2; it was checked for n <= 8 (1728 coefficients).

**Unstable cases.** Both fit the same formula with Route A's genus-0 values Xi_0(T) = Gamma(J)/Gamma(J-r+3) for r = 1, 2:
- (0,1): [x^mu] Int Y_0 dX = sigma Pi_mu/(mu(mu+1)). *Proof:* D F_{0,1} = -1/w and 1/w = (1+w)^(-tau-1)/x. Lagrange gives [x^(mu+1)](1+w)^(-tau-1) = -((tau+1)/(mu+1)) C(mu(tau+1)-1, mu) = -sigma Pi_mu/(mu+1).
- (0,2): [x1^m1 x2^m2] log((w1-w2)/(x1-x2)) = sigma Pi_m1 Pi_m2/(m1+m2) (Zhou (71)-(76)). Note d1 d2 log((w1-w2)/(x1-x2)) = B - dx1dx2/(x1-x2)^2. Checked for m1+m2 <= 9 at 5 values of tau.
- At the BKMP point the same cases give Pi_mu/mu^2 and sigma Pi Pi/(m1+m2), Zhou's unstable conventions (1), (2).

**Vanishing leaf normalisations.** If Pi_mu(tau) = 0 then [x^mu]F = 0 whatever Xi is. This is why TR at tau = -1/2 sees only odd parts through x. The coefficient form of P2 still holds: all coefficients were compared at tau = -1/2, -1/3, 5/7.

**tau <-> -1-tau.** sigma is invariant and Pi_mu(-1-tau) = (-1)^(mu-1) Pi_mu(tau). So C_tau and C_{-1-tau} are different curves with the same Hodge side, and both satisfy P2: tau = 2 and tau = -3 give the same sigma = 6 values.

**Global symmetries (proved; settles referee_hurwitz_1 s.5 item 4 within Aut(P^1)).**
- A Moebius s with s^*dX = dX permutes the poles 0, -1, oo of dX and preserves their residues 1, -(1+tau), tau.
- If the residues are distinct, s fixes three points, so s = id.
- Two residues coincide iff tau = 1 (s(w) = 1/w), tau = -2 (s(w) = -1-w) or tau = -1/2 (s(w) = -w/(1+w)). The involution swapping the two poles with equal residue and fixing the third satisfies s^*dX = dX, because a form on P^1 with only simple poles is determined by its residues.
- So the global involution of Theorem 3.1 exists exactly when two CY weights coincide.

**Alternative route (not needed).** By DOSS, TR on a one-branch-point curve is Givental's formula: the R-matrix comes from the Laplace transform of B (for C_tau a Beta integral, whose Stirling expansion is Mumford's formula for the triple Hodge class), and the translation is linear in y. That gives a second proof without BKMP; its normalisations were not written out.

## 6. Exact verification

**Code.**
- **TRxi.pm** (production): exact EO recursion in the Laplace basis.
  - u(zeta) comes from the ODE in O(L^2) steps and is checked against the Taylor series of X.
  - xi_0 = -1/tau - (1/tau^2)/u, and xi_{k+1} = -(1/(2 c2 zeta)) d xi_k/dzeta.
  - The local basis e_l = (2l+1) zeta^(-2l-2) dzeta + hol. enters only through dxi_k = Sum_l T_kl e_l (T_kl = [zeta^(-2l-2)] xi_k'/(2l+1)) and K = Sum_l e_l(p0) zeta^(2l+1)/((2l+1) omega~).
  - Residue tables: rawAB(k0,a,b) = [zeta^(-2k0)](X_a(zeta) X_b(-zeta)/Yh)/(2 c2 (2k0+1)), and rawC(k0,l,b) = [zeta^(-2k0)](zeta^(2l)(X_b(zeta)+X_b(-zeta))/Yh)/(2 c2 (2k0+1)). In rawC, B(q,p) = Sum_l e_l(p) zeta^(2l) dzeta plus odd terms, which drop out by (LL).
  - The output is converted to the xi-basis by T^-1.
- **TRtau.pm** (independent check): recursion in the local basis with global forms e_k = (2k+1) Sum_j [zeta^(2k+1)]u^j dw/(w-a)^(j+1), and x-expansion by Lagrange inversion.
- **HodgeSide.pm** reads the Hodge tables.
- **run_xi.pl / run_p2.pl** compare every coefficient: all ordered k with |k| <= 3g-3+n, plus a check that nothing lies outside that range.

**Results.** Everything is exact, with 0 failures and 0 symmetry violations everywhere. "coeff/run" is the number of xi-basis coefficients compared per run.

| identity | y | tau | range | coeff/run |
|---|---|---|---|---|
| P2 vs Route E engine | -1/w | 1, 2, 3, -3, 1/2, -1/2, 2/3, -1/3, 5/7, 3/2 | g<=4, n<=3 | 956 |
| P2 | -1/w | 1, -1/2, 2/3 | g<=5, n<=3 | 1906 |
| P2 | -1/w | -1/2 | g<=6, n<=3 | 3424 |
| P2 (one leaf) | -1/w | 1, 2, 3, 1/2, -1/2, -1/3, 2/3 | g<=6, n=1 | 57 |
| P2, genus 0 | -1/w | 2/3 | (0,3)..(0,8) | 1728 |
| P2, g <= 1 | -1/w | -5/7 | (0,3..6), (1,1..3) | 139 |
| (0,3), (1,1) closed forms | -1/w, log(1+w) | 7/3, -5/2, 11, -1/7, 13/5 | | 3 |
| (0,2) | | 1, 2/3, -1/2, 3, -5/7 | m1+m2<=9 | 36 |
| BKMP (Zhou) vs MV integrals | log(1+w) | 2/3, 3, -2/5 | g<=4, n<=3 | 956 |
| Zhou's (6) literally | -log(1+w) | 2/3 | g<=2, n<=3 | 139 |
| Theorem F vs Route A's Q_g(s;R) | log(1+w)-e/w | (2/3,3/5), (2,-7/3), (1/2,4), (-1/2,1) | g=2..4, n<=2 | 168 |
| TRtau.pm (e-basis, Lagrange) | -1/w | 1, 2/3 | g<=2, n<=3 | 139 + 36 lattice points |
| (6.1) | -1/w | 1, 2, 3, 1/2, -1/2, -1/3, 2/3 | h = 2..7, j = 0..3h+1 | identity in j |

**Agreement with the referee's independent tau = 1 code.** Sum_k C_{h,1}[k] j^k (= [x^j]F/Pi_j) reproduces:

| quantity | value |
|---|---|
| h=1, j=1 | -1/16 |
| Xi_2(3) | -1/512 |
| Xi_3(5) | -3/4096 |
| genus 4, j = 7, 8, 9 | -225/262144, -1823/262144, 669/327680 |
| Xi_5(9) | -2205/1048576 |
| Xi_6(11) | -297675/33554432 |

New: at h = 7 and tau = 1 the first non-zero value is at j = 13, equal to -15436575/268435456.

**Run times (one core).**

| range | time |
|---|---|
| genus 4, n <= 3 | 1-4 min |
| genus 5, n <= 3 | 6-15 min |
| genus 6, n <= 3 | 39 min |
| genus 7, n = 1 | 6-22 min |

**(6.1) at h = 7.**
- At sigma in {2, -1/4, 6, 3/4, 12, -2/9, 10/9} (seven values), the exact genus-7 one-leaf TR polynomial in j agrees with (6.1) at 23 values of j. Both sides have j-degree <= 3h-2 = 19.
- With sigma = 0 (the lambda_g limit, proved for all h) this is 8 values of sigma.
- Both sides of (6.1) have sigma-degree <= h:
  - left side: sigma^h Xi_h = Sum_k (-1)^k sigma^(h-k) Xi^(k);
  - right side: with e1+e2 = 1 and e1 e2 = -sigma, the lambda^(2m) coefficient of the exponent is c_m (p_{2m} - j sigma p_{2m-1}). The power sums p_k have sigma-degree floor(k/2), so that coefficient has degree m.
- Hence (6.1) holds identically at h = 7, and by the same data at every h <= 7. For h <= 6 this independently confirms Route E's TR-free verification.
- Zeros: Xi_7(j) = 0 exactly for j = 1..12 at all seven values.

## 7. What becomes unconditional

P2 was the only conditional input of the items below (proof_status_mixed_scale.md). They are now proved modulo the literature inputs of P2, uniformly in the genus:

- **(a) s.1.6(c) / Theorem 3.1.** For every h >= 2, every r and every mu in Z_{>=1}^r with |mu| <= 2h-3+r: Xi_h(mu;tau) = 0 at tau = 1 and tau = -2 (sigma = 2).
  - In particular Lemma A holds at sigma = 2 in every genus: Sum_{k=0}^h (-1/2)^k Xi^(k)_h(R) = 0 for every R with parts >= 2 and D(R) <= 2h-3.
  - TR input: Theorem 3.1 as a TR statement, proved in proof_M_hurwitz.md s.3 and checked line by line in referee_hurwitz_1.md s.1. Its consequence (zeros of the TR polynomial for j <= 2h-2) is confirmed here exactly up to genus 7 for one leaf and genus 5 for n <= 3.
- **(b) s.1.6(c) at tau = -1/2.** For every h and every R with all parts odd and D(R) <= 2h-3: Xi_h(R; -1/2) = 0 (sigma = -1/4). This is the same statement as Route A's TR-free Theorem C.
- **(c) s.1.7(b).** Xi_3(2) = Xi_3(3) = 0 for all tau, without computer data.
- **(d) s.1.7(c).** Lemma D at h distinct integer values of tau >= 1 implies M_irr(h). The implication is unconditional; Lemma D itself is still open for tau >= 2.
- **(e) s.1.8.** The exact TR checks of (6.1) at sigma = 2 (h = 2..6, all j) and sigma = -1/4 (h = 2..6, odd j, enough for an identity in j) are now exact verifications of (6.1) at those sigma. With this report, (6.1) holds identically in (j, tau) for h <= 7 (s.6).
- **(f) brief s.5 / (OL).**
  - At genus G and odd D, the sigma = 2 condition on the one-leaf polynomial sigma^G Xi_G(D+1) is now unconditional. It is one of the G-1 conditions needed.
  - tau = -1/2 gives nothing for (OL): D+1 is even, so Pi_{D+1}(-1/2) = 0 in the x-expansion, and the TR degree bound there only concerns odd parts.
- **(g) Lemma B of proof_status s.2:** proved.

All of (a)-(f) were already known for h <= 20 through Route A (KP) and for h <= 6 through Route E. What is new is that the sigma = 2 statement, and the odd-part statement at sigma = -1/4, hold in every genus by a proof independent of KP.

## 8. What P2 does not give; obstruction

P2 is a dictionary. It transports a statement about the TR correlators of C_tau at one value of tau to Xi at the same tau.
- Lemma A at genus h needs h independent conditions on the tau-independent components Xi^(0..h-1)_h(R) (proof_status 1.7(a)).
- TR supplies them only where a degree bound at the zero v = oo of x is known: tau = 1 (all parts) and tau = -1/2 (odd parts).
- For other tau the missing statement is Lemma D, and the obstruction is unchanged. There is no Moebius symmetry (s.5), and for integer tau >= 2 the correlators have low degree at the other zero v = 0 of x (proof_M_hurwitz s.3), so the residue argument at v = oo does not close with the induction hypothesis alone.
- Theorem F is a family at fixed tau and produces no new values of sigma.

The exact checks of s.6 test the identities the argument relies on beyond the earlier ranges:
- Zhou's normalisation, to genus 4;
- the interpolating family (variational formula + rationality), to genus 4;
- P2, to genus 6 for all coefficients with n <= 3, and to genus 6/7 for one leaf.

## 9. Files (HVC\lemmaA\round2\TR1_dictionary)

**Modules and scripts.**

| file | role |
|---|---|
| TRxi.pm | exact TR on C_tau in the Laplace basis (production) |
| TRtau.pm | exact TR in the local basis with global forms and Lagrange x-expansion (independent check) |
| HodgeSide.pm | Hodge side: Route E engine tables (X), Mathematica tables mv_*/q_* (M, H), genus-0 formulas |
| run_xi.pl tau G N [base\|bkmp\|zhou\|mix] [e] [nmin] | runs TRxi and compares every coefficient. Env SAVE=file stores the xi-coefficients; TRVERB=1 prints timings |
| run_p2.pl tau G N [ycase] [e] | the same with TRtau.pm, plus x-expansion spot checks |
| check_02.pl tau [D] | the (0,2) case |
| cmp61_xi.pl tau savefile | compares saved one-leaf TR values with (6.1) |
| prof_setup.pl | timing diagnostic |
| MParse.pm, Poly.pm, WK.pm | copies from E_formulas |

**Outputs (runs\).**
- p2_tau<t>_g<G>n<N>.out and .err
- bkmp_*.out, mix_*.out
- save_tau<t>_g7n1*.txt
- cmp61_h2to7.txt, check02.txt

**Examples.**

    perl run_xi.pl 2/3 4 3 base
    perl run_xi.pl 2/3 4 3 bkmp
    perl run_xi.pl 2/3 4 2 mix 3/5
    SAVE=s.txt perl run_xi.pl 1 7 1
    perl cmp61_xi.pl 1 s.txt
    perl check_02.pl 2/3 9

**Literature.**
- J. Zhou, arXiv:0910.4320v2.
- L. Chen, arXiv:0910.3739.
- B. Eynard, N. Orantin, CNTP 1 (2007) 347-452 [EO07].
- B. Eynard, N. Orantin, CMP 337 (2015) 483-567.
- B. Fang, C.-C. M. Liu, Z. Zong, JAMS 33 (2020) 135-222.
- P. Dunin-Barkowski, N. Orantin, S. Shadrin, L. Spitz, CMP 328 (2014) 669-700.
- G. Borot, B. Eynard, N. Orantin, CNTP 9 (2015) 51-187.
- C.-C. M. Liu, K. Liu, J. Zhou, JDG 65 (2003).
