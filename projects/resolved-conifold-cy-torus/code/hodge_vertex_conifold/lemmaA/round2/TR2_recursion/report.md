# Lemma A round 2: route TR2 (saved by the main session from the workflow journal, 2026-10-31)

Subagents cannot write report files, so this file holds the agent's structured return value: the fields first, then the full report text (report_text) verbatim.

## all_genus_claim

1

## all_genus_statement

Two statements, with different dependencies. First, a pure TR statement on C_tau that needs only the published x-y swap theorem for log TR (ABDKS, IMRN 2024, Thm 2.8/Prop. 3.1): for every tau not in {0,-1} and all g, n, the EO correlators of C_tau have x-degree >= 2g-2+n at v = oo (Lemma D), with a uniform proof and no induction; the closed one-point formula (5.1)-(5.2) also holds. Second, adding P2: Lemma A = M_irr(h) holds for all h >= 2, all R and all tau, so (OL) holds for every G and every level, Conjecture M (vanishing) holds for all h and r, Corollary R holds in all genera, and formula (6.1) holds for all h.

## claims

### [1]

#### depends_on

A. Alexandrov, B. Bychkov, P. Dunin-Barkowski, M. Kazarian, S. Shadrin, IMRN 2024 no. 21, 13461-13487 (arXiv:2312.16950), Thm 2.8 and Prop. 3.1; universal x-y swap, JEMS 2025 (arXiv:2212.00320), Thm 1.14

#### evidence

Hypotheses checked: dX meromorphic with one simple zero; dy = dv without zeros and non-zero at v = tau; y meromorphic, so LogTR = TR; v = oo is covered by Lemma 5.5. Normalisation and transcription checked against exact EO recursion: n = 1 exact rationals (132 values, g <= 5, tau = 2 (g <= 4), 3, 3/7, -5/3, 0 mismatches); the referee's tau = 1 genus-5 and genus-6 values; pointwise mod p at 22 points (n = 1..4, g <= 5, 9 values of tau).

#### statement

Formula (17) of ABDKS (Log TR through the prism of x-y swap, IMRN 2024, Thm 2.8 + Prop. 3.1) applies to C_tau for every tau not in {0,-1}. In their convention (omega^(0)_1 = -y dx) the curve is z = v, x = X = tau log v - (tau+1) log(1+v), y = v, with dual vital points a = 0 (1/alpha = tau) and a = -1 (1/alpha = -(tau+1)) and R = 0. All EO correlators of C_tau are then given in closed form by (2.1)-(2.2) of the report.

#### status

proved-modulo-literature

### [2]

#### depends_on

claim 1

#### evidence

Scaling argument (report s.3): restricted to t = lambda theta, every term of (2.1) is O(lambda^(2g-2+n)), because each hbar comes with (v-a)^(-1) or (v_i-v_j)^(-1). Direct TR checks: tau = 3/7 (g <= 4, n <= 3); tau = 2 and 1/2 (g <= 5, n <= 2): 238 coefficients, 0 violations.

#### statement

Theorem 1 (Lemma D for every tau, g, n). For tau not in {0,-1} and 2g-2+n > 0, omega_{g,n} of C_tau has total x-degree >= 2g-2+n at v = oo: [x^mu]F_{g,n} = 0 for all mu_i >= 1 with |mu| <= 2g-3+n. The proof is uniform in g and n and uses no induction.

#### status

proved-modulo-literature

### [3]

#### depends_on

claim 2; P2 = Lemma B (route TR1); proof_status 1.1(e)

#### evidence

Report s.4: by P2 and Theorem 1, Xi_h(R; tau) = 0 for every integer tau >= 1, where the binomials are positive; polynomiality in 1/sigma then gives vanishing identically. P2 is used only here and in Corollary 4.

#### statement

Corollary 2. Assuming P2: Lemma A = M_irr(h) holds for all h >= 2, all multisets R and all tau. Hence (OL) holds for every G and every level D (odd or even, with no inductive hypothesis), Conjecture M (vanishing direction) holds for all h and r, and Corollary R holds in all genera, conditional only on P2 and ABDKS.

#### status

proved-modulo-literature

### [4]

#### depends_on

claim 1; claim 6

#### evidence

Integration by parts in (2.1) and the Norlund identities (claim 6). Checks: closed form = full formula at tau = 2 (g <= 5), 3/7 and -5/3 (g <= 4); closed form = (6.1) up to genus 9.

#### statement

Theorem 3 (closed one-point function of C_tau, TR side, no P2). [x^j]F_{g,1} = -(1/j^2)[hbar^(2g)] Res_{v=oo} exp(-Sum_{a=0,-1} [sinh(j hbar d/2)/sinh(alpha_a hbar d/2)] log(v-a)) dv = (1/j^2) Sum_{p=0}^g A_{g-p} B_p C(j(tau+1)-2p, j+1-2g), with A_m = (-N)_{2m}[y^(2m)]S1^(1-N) tau^(-2m), B_p = (M)_{2p}[y^(2p)]S1^(1+M)(tau+1)^(-2p), N = j tau, M = j(tau+1). For integer j tau it is [hbar^(2g)][w^(j+1)] prod(1+w+hbar s'_k w)/prod(1+hbar s_k w) divided by j^2.

#### status

proved-modulo-literature

### [5]

#### depends_on

claim 4; P2

#### evidence

Term-by-term identity (-N)_{2g-2p}(M)_{2p} C(M-2p, j+1-2g)/C(M-1,j-1) = M N (j-1)!/(j+1-2g)! applied to (5.2) (report s.5). Numerical checks: 310 values, g <= 9 at tau = 2, g <= 8 at tau = 5/3, g <= 6 at tau = -2/7, 0 mismatches. Hand checks at g = 1, 2 against known data.

#### statement

Corollary 4: the one-leaf formula (6.1), sigma^h Xi_h(j) = prod_{l=1}^{2h-2}(j-l) [lambda^(2h)] S(e1 lambda)^(1+j e2) S(e2 lambda)^(1+j e1), holds for all h >= 1, all j >= 1 and all tau.

#### status

proved-modulo-literature

### [6]

#### evidence

For integer N: the geometric-sum form of sinh(Nu)/sinh(u) plus B_N^(N+1)(x) = (x-1)...(x-N) (a residue computation), or Tricomi-Erdelyi (DLMF 5.11.17) for the inverse. All tau then follow from Laurent-polynomial identities in tau at fixed j. norlund_check.pl: 224 checks at tau = 3/7, -13/5, 2, 0 mismatches.

#### statement

Norlund identities: exp(Sum_k A_k t^(2k)) = Sum_n (-N)_n [y^n] S1(y)^(1-N) (t/tau)^n and exp(Sum_k B_k t^(2k)) = Sum_n (M)_n [y^n] S1(y)^(1+M) (t/(tau+1))^n, for all tau. Here A_k = (2k-1)! j tau c_k(j,1/tau), B_k = -(2k-1)! j(tau+1) c_k(j,-1/(tau+1)), S1(y) = (y/2)/sinh(y/2).

#### status

proved

### [7]

#### depends_on

claim 1

#### evidence

Integration by parts in (2.1). Exact check against TR mod each of 3 primes: 148 coefficients (tau = 2: g <= 3; tau = 3/7: g <= 3; tau = -5/3: g <= 4), 0 mismatches.

#### statement

Two-point closed formula: [x1^j1 x2^j2]F_{g,2} = (1/(j1 j2))[hbar^(2g)] Res_{v1=oo}Res_{v2=oo} Psi_{j1}(v1)Psi_{j2}(v2) dv1 dv2/((v1-v2)^2 - (j1+j2)^2 hbar^2/4), with |v2| > |v1| in the inner residue and Psi_j = exp(-j S(j hbar d) xhat).

#### status

proved-modulo-literature

### [8]

#### evidence

Standard linear and quadratic loop equations of TR together with the global linear loop equation (report s.6(1)). trace_qle.pl: holds at tau = 2 (g <= 4), 3/7 (g <= 3), 1 (g <= 4).

#### statement

Proposition 5 (global quadratic loop equation in trace form): Sum over the fibre of x of [Brk_g(v,v)/(2dX^2) - v omega_{g,1}(v)/dX] = 0 on C_tau. This holds for integer tau, and for all tau in its equivalent local form (the sigma-even polar part at v = tau vanishes).

#### status

proved

### [9]

#### depends_on

claim 8

#### evidence

Report s.6(1): anti-invariance at tau = 1; global linear loop equation for tau >= 2.

#### statement

At tau = 1, Proposition 5 gives x dF_{g,1}/dx = t Brk_g(v,v)/((1-t^2) dX^2) at v = 1/t, a one-line proof of Theorem 3.1 for n = 1 given the bound in lower genera. For integer tau >= 2 it reduces Lemma D(g,1) to the vanishing of the tau-section of the v = 0 expansion of Phi_g up to order (2g-3)tau, which involves data at v = 0 not controlled by induction: this is the precise obstruction of the loop-equation route.

#### status

proved

### [10]

#### depends_on

claims 4, 5 (TR side, mod ABDKS); Route E engine

#### evidence

Route E verified (6.1) identically against its Hodge engine for h <= 7 (brief s.2(g)); Theorem 3 with the term-by-term identity shows that the TR side equals (6.1) for all h, without P2.

#### statement

One-leaf P2 cross-check: TR on C_tau and the Hodge integrals agree for r = 1 and h <= 7, identically in tau.

#### status

computer-verified

### [11]

#### evidence

lemmaD_check.pl; out/lemmaD_tau2_g5.txt, out/lemmaD_tau1_2_g5.txt; 238 coefficients, 0 violations.

#### statement

Exact TR on C_tau satisfies the degree bound directly: tau = 3/7 for g <= 4, n <= 3; tau = 2 and tau = 1/2 for g <= 5, n <= 2; the first layer |mu| = 2g-2+n is non-zero, except that at tau = 1/2 the all-even first-layer coefficients vanish.

#### status

computer-verified


## key_formulas

- C_tau in the ABDKS convention (omega^(0)_1 = -y dx): x = X = tau log v - (tau+1) log(1+v), y = v; dual vital points (a, 1/alpha) = (0, tau), (-1, -(tau+1)); R = 0
- omega_{g,n} = -[hbar^(2g-2+n)] prod_i ( Sum_{r_i} (d_i(1/dX_i))^(r_i) [u_i^(r_i)] E_i ) Sum_{sigma n-cycle} prod_i dv_i/(v_i - v_sigma(i) + (u_i+u_sigma(i)) hbar/2)   (n = 1: cycle factor dv/(u hbar))
- E_i = exp( Sum_{k>=1} hbar^(2k) u_i (2k-1)! Sum_{a in {0,-1}} c_k(u_i,alpha_a)/(alpha_a (v_i-a)^(2k)) ),  c_k(u,al) = [t^(2k)] S(ut)/S(al t),  S(t) = 2 sinh(t/2)/t
- Theorem 1: ord_{v=oo} omega_{g,n} >= 2g-2+n, i.e. [x^mu]F_{g,n} = 0 for |mu| <= 2g-3+n (all tau not in {0,-1}, all g, n)
- [x^j]F_{g,1} = -(1/j^2) [hbar^(2g)] Res_{v=oo} exp(-Sum_a [sinh(j hbar d/2)/sinh(alpha_a hbar d/2)] log(v-a)) dv
- [x^j]F_{g,1} = (1/j^2) Sum_{p=0}^g A_{g-p} B_p C(j(tau+1)-2p, j+1-2g),  A_m = (-N)_{2m}[y^(2m)]S1(y)^(1-N) tau^(-2m),  B_p = (M)_{2p}[y^(2p)]S1(y)^(1+M)(tau+1)^(-2p),  N = j tau, M = j(tau+1), S1(y) = (y/2)/sinh(y/2)
- integer j tau, g >= 1: [x^j]F_{g,1} = (1/j^2)[hbar^(2g)][w^(j+1)] prod_{k<M}(1+w+hbar s'_k w)/prod_{k<N}(1+hbar s_k w),  s_k = (k-(N-1)/2)/tau,  s'_k = (k-(M-1)/2)/(tau+1)
- term-by-term identity giving (6.1): (-N)_{2g-2p} (M)_{2p} C(M-2p, j+1-2g)/C(M-1, j-1) = M N (j-1)!/(j+1-2g)!,  MN/j^2 = sigma
- Norlund: prod_{k=0}^{N-1}(1+s_k t) = Sum_n (N)_n [y^n] S1(y)^(N+1) (t/tau)^n   (from B_N^(N+1)(x) = (x-1)...(x-N))
- [x1^j1 x2^j2]F_{g,2} = (1/(j1 j2))[hbar^(2g)] Res_{v1=oo}Res_{v2=oo} Psi_{j1}(v1)Psi_{j2}(v2) dv1 dv2/((v1-v2)^2 - (j1+j2)^2 hbar^2/4)   (|v2| > |v1|)
- trace-form QLE: Tr_x [Brk_g(v,v)/(2 dX^2) - v omega_{g,1}(v)/dX] = 0;  at tau = 1: x F'_{g,1}(x) = t Brk_g(v,v)/((1-t^2) dX^2), v = 1/t

## next_steps

- Combine with route TR1's write-up of P2. Lemma A, (OL), Conjecture M (vanishing) and (6.1) then hold in all genera, modulo ABDKS, and Corollary R holds in all genera. Update proof_status_mixed_scale.md s.2 and the regularity note accordingly.
- Referee: check the hypotheses of ABDKS Thm 2.8 / Prop. 3.1 for C_tau, in particular Lemma 5.5 at v = oo and the sign convention omega^(0)_1 = -y dx (their y = v). Check the degree argument of s.3 and the integration by parts and Norlund identities of s.5.
- Optional, to remove the ABDKS dependence for one leaf: derive the closed one-point formula (5.1)-(5.2), equivalently (6.1), directly on the Hodge side. Possible routes are MV with unit legs plus the Laguerre transform (Route A, Prop. 2.1), or Faber-Pandharipande localization. With Corollary E (KP) this would give Lemma A without P2 or ABDKS.
- Use the general-n closed formula (s.5) to prove Route E's conjectural one-leaf universality (SN) and two-leaf form (NR). Also extract the first non-vanishing layer (degree exactly 2g-2+n) explicitly.
- Run tr2_closed_check.wls (symbolic in tau and j) as an independent Mathematica check of the Norlund identities and of the term-by-term (6.1) identity at fixed g.

## obstructions

- Conditional on P2 (Lemma B, the dictionary between TR on C_tau and the Hodge integrals Xi_h), which route TR1 is writing out. P2 is used only to translate Theorem 1 and Theorem 3 (pure TR statements on C_tau) into Lemma A and (6.1).
- Conditional on the published x-y swap theorem for logarithmic TR (ABDKS, IMRN 2024, Thm 2.8 and Prop. 3.1), which was not re-proved here. Its hypotheses were verified for C_tau, and the transcribed formula was checked against exact EO recursion for n = 1..4 up to genus 5. A referee should re-read its Lemma 5.5 (the point v = oo, where dX has a simple pole and dy a double pole) and the transcription.
- Without the x-y swap, the global-loop-equation route of line (1) does not close for integer tau >= 2. The trace-form quadratic loop equation equates the v = oo expansion of x F'_{g,1} with the tau-section of the v = 0 expansion of Phi_g. That section involves omega_{g,1} and lower correlators at v = 0 (where they have degree 1), and the induction gives no control there.
- The general-n integrated closed formula (iterated residues of prod Psi_{mu_i} times cycle factors) was derived but checked only for n <= 2.

## requested_checks

- Mathematica: run TR2_recursion\tr2_closed_check.wls <gmax> (e.g. gmax = 6). Symbolic in tau and j, it checks the two Norlund identities and the term-by-term identity that turns (5.2) into (6.1); expected: all differences 0.
- Referee: read ABDKS (IMRN 2024, arXiv:2312.16950) Thm 2.8, Prop. 3.1 and Lemma 5.5. Confirm that C_tau (z = v, x = X, y = v, vital points 0 and -1, R = 0) satisfies the hypotheses, including the point v = oo (simple pole of dX at a double pole of dy), and that formula (17) as transcribed in report s.2.1 is correct. The transcription is validated numerically in report s.7.
- Optional independent check: implement (2.1) for n = 2, 3 symbolically in Mathematica and compare with hvc_tr_tau.wls at tau = 2, 3 (g <= 3, n <= 3).

## route

TR2 (round 2): topological recursion on C_tau, the Hodge CohFT / framed C^3 vertex base curve, through the x-y swap theorem for logarithmic TR (Alexandrov-Bychkov-Dunin-Barkowski-Kazarian-Shadrin, IMRN 2024). The theorem gives closed formulas for all correlators. A scaling argument then gives Lemma D for every tau, g and n. With P2 this gives Lemma A, (OL) and (6.1) in all genera.

## summary

The degree bound behind Lemma A now holds for all genera, all numbers of points and every tau, not only integer tau. It rests on one published theorem: the x-y swap theorem for logarithmic topological recursion, A. Alexandrov, B. Bychkov, P. Dunin-Barkowski, M. Kazarian, S. Shadrin, IMRN 2024, Thm 2.8 and Prop. 3.1 [ABDKS]. On C_tau, y = -v is an affine coordinate and X = tau log v - (tau+1) log(1+v) is a sum of two logarithms. The dual curve therefore has no branch points, only two log-vital points (v = 0 and v = -1), and the theorem gives every correlator of C_tau in closed form (formula (2.1) of the report). In that formula each power of hbar comes with a factor (v-a)^(-1), a in {0,-1}, or (v_i-v_j)^(-1). A scaling argument then shows that omega_{g,n} has total x-degree at least 2g-2+n at v = oo, i.e. Lemma D, uniformly in g and n with no induction. I checked the hypotheses of [ABDKS] for C_tau. I also checked that my transcription of its formula reproduces exact Eynard-Orantin recursion on C_tau: n = 1 as exact rationals through genus 5 at tau = 2 (genus 4), 3, 3/7 and -5/3; the referee's independent genus-5 and genus-6 values at tau = 1; and n = 2, 3, 4 pointwise mod p at 21 points up to genus 5, plus genus 5 for n = 2 at tau = 3/7. Lemma D was also checked directly in TR through genus 5 (n <= 2, tau = 2 and 1/2) and genus 4 (n <= 3, tau = 3/7).

P2 is used once, to translate this into Hodge integrals. It gives Lemma A for all h, all R and all tau. So (OL) holds for every G and every level, odd or even, with no inductive hypothesis. Conjecture M (vanishing direction) holds for all h and r. Corollary R holds in all genera, conditional only on P2 and [ABDKS]. KP, Corollary E, Theorem D5 and the parity argument are not needed. The same formula gives the one-point function of C_tau in closed form, [x^j]F_{g,1} = -(1/j^2)[hbar^(2g)] Res_{v=oo} exp(-j S(j hbar d) xhat) dv, as a finite Norlund-type sum. Matching that sum term by term proves the conjectural one-leaf formula (6.1) for all h (mod P2); it was also checked against (6.1) through genus 9. A closed two-point formula was checked exactly against TR. For line (1) I proved the global quadratic loop equation of C_tau in trace form. At tau = 1 it re-proves Theorem 3.1 for n = 1 in one line. For tau >= 2 it isolates the exact obstruction of the loop-equation route: the tau-section of the v = 0 expansion, about which induction says nothing. What remains: P2 (route TR1), and an independent referee reading of [ABDKS] Lemma 5.5 (the case at v = oo, where dX has a simple pole and dy a double pole) and of my transcription of formula (17). The transcription is validated numerically here but the theorem is not re-proved.

## report_text

# Route TR2 (round 2): topological recursion on C_tau through the x-y swap. Lemma D for every tau, g, n; Lemma A, (OL) and (6.1) in all genera (mod P2)

Folder: HVC\lemmaA\round2\TR2_recursion (HVC = G:\My Drive\Projects\Old\refGW\hodge_vertex_conifold).
- Arithmetic is exact throughout. Perl with Math::BigRat is used directly. The large EO runs use exact arithmetic modulo three primes near 2^31, combined by CRT and rational reconstruction, and are also compared separately modulo each prime.
- No wolframscript was run.
- TR and the x-y swap are applied only to the base curve C_tau of the Hodge CohFT (the triple Hodge class with CY weights = the framed C^3 vertex; brief s.0 and s.5). Nothing is applied to the resolved conifold x C^2, whose input (H) remains Brini-Schuler.
- KP, Corollary E, Theorem D5 and the parity argument are not used.

## 0. Results

- **Theorem 1 (Lemma D for every tau, g, n).** Let tau be in C minus {0,-1}.
  - Statement: for all g >= 0, n >= 1 with 2g-2+n > 0, the Eynard-Orantin correlators of C_tau, expanded at v = oo in x as in P2, satisfy [x_1^mu_1 ... x_n^mu_n] F_{g,n} = 0 whenever all mu_i >= 1 and |mu| <= 2g-3+n.
  - Proof: one closed formula valid for all (g,n), followed by a termwise degree count. It is uniform in g and n and uses no induction.
  - Status: proved modulo one published theorem, the x-y swap theorem for logarithmic TR [ABDKS] (Thm 2.8 and Prop. 3.1). Its hypotheses hold for C_tau (s.2.2). My transcription of its formula reproduces exact EO recursion on C_tau for n = 1, 2, 3, 4 (s.7).
- **Corollary 2 (uses P2, exactly once).**
  - Lemma A = M_irr(h) holds for all h >= 2, all R and all tau.
  - Hence (OL) holds for every G and every D, odd or even, with no inductive hypothesis.
  - Conjecture M (vanishing direction) holds for all h and r.
  - Corollary R (brief s.1) holds in all genera, conditional only on P2 and [ABDKS].
- **Theorem 3 (closed one-point function of C_tau; TR side, no P2).**
  - [x^j] F_{g,1} = -(1/j^2) [hbar^(2g)] Res_{v=oo} exp(-Sum_{a=0,-1} [sinh(j hbar d_v/2)/sinh(alpha_a hbar d_v/2)] log(v-a)) dv, with alpha_0 = 1/tau and alpha_{-1} = -1/(tau+1).
  - It has the explicit finite form (5.2).
  - There is also a closed two-point formula (5.4), checked exactly against TR.
- **Corollary 4 (uses P2).** Formula (6.1) of proof_M_residue.md s.6 holds for all h >= 1, all j >= 1 and all tau. The identity holds term by term in a single summation index of (5.2).
- **Proposition 5 (line (1)).** The global quadratic loop equation of C_tau holds in trace form. At tau = 1 it gives Theorem 3.1 for n = 1 in one line. For tau >= 2 it isolates the precise obstruction of the loop-equation route.

## 1. Setting

C_tau has Sigma = P^1 with coordinate v and X = log x = tau log v - (tau+1) log(1+v).
- dX = (tau - v) dv/(v(1+v)) has exactly one zero, a simple zero at v = tau, which is the only branch point. It has simple poles at v = 0, -1, oo with residues tau, -(tau+1), +1.
- omega_{0,1} = y dX with y = -v, and omega_{0,2} = B = dv1 dv2/(v1-v2)^2. These are the conventions of P2, hvc_tr_tau.wls, TRtau.pm and TRmod.pm.
- At v = oo put t = 1/v; then x = t(1+t)^(-tau-1). F_{g,n} = int_oo ... int_oo omega_{g,n}.
- P2: Xi_h(mu; tau) = sigma^(h-1) [x^mu] F_{h,r} / prod_i C(mu_i(tau+1)-1, mu_i-1), with sigma = tau(1+tau).
- Notation:
  - S(u) = (e^(u/2) - e^(-u/2))/u is the function of [ABDKS].
  - S1(y) = 1/S(y) = (y/2)/sinh(y/2) is the S of (6.1).
  - (a)_n = a(a-1)...(a-n+1) is the falling factorial.
  - N = j tau, M = j(tau+1).

## 2. The x-y swap theorem and C_tau

### 2.1 The cited result

- [ABDKS] A. Alexandrov, B. Bychkov, P. Dunin-Barkowski, M. Kazarian, S. Shadrin, "Log topological recursion through the prism of x-y swap", Int. Math. Res. Not. IMRN 2024, no. 21, 13461-13487, doi:10.1093/imrn/rnae213, arXiv:2312.16950.
- It extends [ABDBKS]: same authors, "A universal formula for the x-y swap in topological recursion", J. Eur. Math. Soc. (2025), doi:10.4171/JEMS/1615, arXiv:2212.00320 (Thm 1.14 there).

I read it from the arXiv HTML (v3). The formula below is my transcription; s.7 validates it numerically.

- **Standing assumptions (s.2.3 of [ABDKS]).**
  - Sigma is compact.
  - dx and dy are meromorphic and may have residues, so x and y may be logarithmic.
  - The zeros of dx and of dy are simple and disjoint.
  - dy is regular and non-zero at the zeros of dx, and dx is regular and non-zero at the zeros of dy.
  - TR conventions there: omega^(0)_1 = -y dx, omega^(0)_2 = B, and the ramification points are the zeros of dx.
- **LogTR (Def. 2.4-2.5).** LogTR is TR plus extra terms at the "LogTR-vital" points: the simple poles of dy with non-zero residue at which dx is regular. These terms enter only omega^(g)_1. If y is meromorphic there is no vital point, and LogTR(x,y) = TR(x,y) (Remark 2.12).
- **Theorem 2.8.** LogTR(x,y) and LogTR(y,x) are related by the universal x-y swap relation of [ABDBKS].
  - The proof is a local analysis at the zeros of dx and dy, at the diagonals and at the poles of dx and dy (Lemmas 5.1-5.6).
  - Lemma 5.5 covers a simple pole of dx, including the case where that point is also a pole of dy.
- **Proposition 3.1 (explicit case).**
  - Hypotheses: Sigma = P^1; y = z is an affine coordinate, with z = oo a pole of dy; x = Sum_{i=1}^M alpha_i^(-1) log(z - a_i) + R(z).
  - Here the a_i are the simple poles of dx at which dy is regular, and 1/alpha_i = Res_{a_i} dx.
  - Then dx^vee = dy has no zeros, and the dual LogTR consists of the vital terms only.
  - For 2g-2+n > 0 the swap relation becomes

        omega_{g,n} = [hbar^(2g-2+n)] (-1)^n prod_{i=1}^n ( Sum_{r_i>=0} (d_i(1/dx_i))^(r_i) [u_i^(r_i)] E_i ) (-1)^(n-1) Sum_{sigma n-cycle} prod_i sqrt(dw_i dwbar_sigma(i)) / (w_i - wbar_sigma(i)),   (17)

    where
    - E_i = exp(-u_i (S(u_i hbar d_y) xhat(z_i) - x(z_i))),
    - xhat = Sum_i alpha_i^(-1) S(alpha_i hbar d_y)^(-1) log(z - a_i) + R(z),
    - w_i = exp(u_i hbar d_{y_i}/2) z_i and wbar_i = exp(-u_i hbar d_{y_i}/2) z_i,
    - (d(1/dx)) eta := d(eta/dx).
  - With y = z: w_i - wbar_j = z_i - z_j + (u_i+u_j) hbar/2, and the square roots give prod dz_i. For n = 1 the cycle factor is dz/(u hbar).

### 2.2 Application to C_tau

Take z = v, x = X and, in the convention of [ABDKS], y = v. Then -y dx = -v dX is our omega_{0,1}, and B agrees.

- **Assumptions.**
  - dX is meromorphic with one zero, at v = tau. That zero is simple for tau not in {0,-1}.
  - dy = dv has no zeros and is regular and non-zero at v = tau.
- **No correction terms on our side.** y is meromorphic, so LogTR(X,v) = TR(X,v) = the correlators of P2.
- **Dual vital points.**
  - a_0 = 0 with 1/alpha_0 = tau.
  - a_{-1} = -1 with 1/alpha_{-1} = -(tau+1).
  - R = 0.
  - v = oo is a log point of X at the pole of dy, so it is not vital (it is covered by Lemma 5.5).
- **Non-integer tau** is allowed, since only dX enters.

Hence, for 2g-2+n > 0,

    omega_{g,n} = - [hbar^(2g-2+n)] prod_i ( Sum_{r_i} (d_i(1/dX_i))^(r_i) [u_i^(r_i)] E_i ) Sum_{sigma n-cycle} prod_i dv_i / (v_i - v_sigma(i) + (u_i+u_sigma(i)) hbar/2),   (2.1)

with the cycle factor dv/(u hbar) for n = 1, and

    E_i = exp( Sum_{k>=1} hbar^(2k) u_i (2k-1)! Sum_{a in {0,-1}} c_k(u_i,alpha_a) / (alpha_a (v_i - a)^(2k)) ),   c_k(u,al) = [t^(2k)] S(u t)/S(al t).   (2.2)

(2.2) is the expansion of -u(S(u hbar d)/S(alpha hbar d) - 1) log(v-a)/alpha, using d^(2k) log(v-a) = -(2k-1)!/(v-a)^(2k).

## 3. Proof of Theorem 1

Fix (g,n) with 2g-2+n > 0 and put t_i = 1/v_i.

- omega_{g,n} has poles only at v_i = tau, so omega_{g,n} = W(t) prod dt_i/t_i with W rational and holomorphic near t = 0.
- Fix theta in (C^*)^n with distinct entries. For a rational function f(t), let ord f be the order at lambda = 0 of f(lambda theta), for generic theta.
- Properties of ord:
  - ord(fg) = ord f + ord g, and ord(f+g) >= min(ord f, ord g).
  - Factors (1 + O(t_i)) do not change ord.
  - ord(t_i df/dt_i) >= ord f, because (t_i df/dt_i)(lambda theta) = theta_i d/dtheta_i [f(lambda theta)].

Work in the frame prod dt_i/t_i. The ingredients of (2.1) behave as follows.

- dv_i = -t_i^(-1) (dt_i/t_i), so each dv_i contributes ord -1.
- dX_i = q(t_i) dt_i/t_i with q(t) = (1 - tau t)/(1+t). The operator d_i(1/dX_i) acts on the coefficient psi of dt_i/t_i by psi -> t_i d/dt_i (psi/q(t_i)), so it does not lower ord.
- (v_i - a)^(-2k) = t_i^(2k)(1 - a t_i)^(-2k) has ord 2k. By (2.2), every u-coefficient of [hbar^(2m)] E_i has ord >= 2m.
- For n >= 2: 1/(v_i - v_j + c hbar) = Sum_m (-c hbar)^m (v_i - v_j)^(-m-1), and (v_i - v_j)^(-1) = t_i t_j/(t_j - t_i) has ord 1 because theta_i != theta_j. So [hbar^m] of the product over a cycle has ord >= m + n.
- For n = 1 the cycle factor 1/(u hbar) has ord 0 at hbar^(-1).

In both cases [hbar^K](E_1...E_n x cycle sum) has ord >= K + n, for every monomial in the u_i.

- With K = 2g-2+n and the n factors dv_i, every term of (2.1) has ord >= 2g-2+n. Hence ord W >= 2g-2+n.
- W is holomorphic at 0, so W(lambda theta) = Sum_d lambda^d W_d(theta), where the W_d are the homogeneous Taylor components.
- Therefore W_d(theta) = 0 for d < 2g-2+n and generic theta, so W_d = 0.
- The change of variables x_i = t_i(1+t_i)^(-tau-1) is tangent to the identity in each variable. So no x-monomial of degree < 2g-2+n occurs either.
- Finally [x^mu] F_{g,n} = [x^mu](W in the x-frame)/prod mu_i. QED.

Remarks.
- (a) Only the structure of (2.1)-(2.2) is used: every power of hbar carries a factor (v-a)^(-1), a in {0,-1}, or (v_i-v_j)^(-1).
- (b) At the other zero v = 0 of x the factors v^(-2k) blow up. The poles cancel only in the sum (the vital-point analysis, [ABDKS] Lemma 5.2), so no bound follows there. This is consistent with the degree-1 behaviour at v = 0 found in proof_M_hurwitz.md s.3.
- (c) At tau = 1 this is an independent proof of Theorem 3.1.
- (d) At tau = -1/2 both vital points have residue -1/2, and v -> -1-v exchanges them. This is the TR explanation of the tau = -1/2 involution.

## 4. Corollary 2: the only use of P2 for Lemma A

Let R have r parts >= 1 with |R| <= 2h-3+r.
- By P2, [x^R] F_{h,r} = sigma^(1-h) prod_i C(R_i(tau+1)-1, R_i-1) Xi_h(R; tau).
- For integer tau >= 1 the binomials are positive. Theorem 1 gives [x^R] F_{h,r} = 0, so Xi_h(R; tau) = 0 for every integer tau >= 1.
- Xi_h(R; tau) is a polynomial in 1/sigma (proof_status 1.1(e)), so it vanishes identically.
- For parts >= 2, D(R) <= 2h-3 is the same as |R| <= 2h-3+r. This is Lemma A.
- R = {D+1} is (OL), for every G and D.

P2 is used here and in Corollary 4, nowhere else. Theorems 1 and 3 and Proposition 5 are statements about TR on C_tau and do not depend on P2.

## 5. Closed one- and two-point functions; proof of (6.1)

**Theorem 3.** Let tau not be in {0,-1}, g >= 1 and j >= 1, with N = j tau and M = j(tau+1). Then

    [x^j] F_{g,1} = -(1/j^2) [hbar^(2g)] Res_{v=oo} Psi_j(v) dv,   Psi_j := exp(-j S(j hbar d_v) xhat) = exp(-Sum_{a in {0,-1}} [sinh(j hbar d_v/2)/sinh(alpha_a hbar d_v/2)] log(v-a))   (5.1)
                  = (1/j^2) Sum_{p=0}^{g} A_{g-p} B_p C(M - 2p, j + 1 - 2g),                                                                      (5.2)

where A_m = (-N)_{2m} [y^(2m)] S1(y)^(1-N) tau^(-2m) and B_p = (M)_{2p} [y^(2p)] S1(y)^(1+M) (tau+1)^(-2p).

When N is a positive integer, Psi_j = prod_{k=0}^{M-1}(v + 1 + hbar s'_k) / prod_{k=0}^{N-1}(v + hbar s_k), with s_k = (k-(N-1)/2)/tau and s'_k = (k-(M-1)/2)/(tau+1). Then, for g >= 1,

    [x^j] F_{g,1} = (1/j^2) [hbar^(2g)] [w^(j+1)] prod_{k=0}^{M-1}(1 + w + hbar s'_k w) / prod_{k=0}^{N-1}(1 + hbar s_k w).   (5.3)

Every hbar comes with a w. So (5.3) vanishes for j <= 2g-2, and at fixed j the hbar-series is a polynomial in hbar^2 of degree <= (j+1)/2. This is Route A's Cor. 2.3 form of Lemma A for one leaf.

Proof.
- (i) The case n = 1 of (2.1) reads omega_{g,1} = -[hbar^(2g)] Sum_r (d(1/dX))^r ([u^(r+1)]E) dv.
- (ii) Start from [x^j] F = (1/j) Res_{v=oo} x^(-j) omega_{g,1}.
  - For r >= 1 the term is d D^(r-1) psi, with psi = [u^(r+1)]E dv/dX and D = d/dX.
  - Since x^(-j) = e^(-jX), r integrations by parts give Res x^(-j) d D^(r-1) psi = j^r Res x^(-j) [u^(r+1)]E dv. The r = 0 term is the same with j^0.
  - Sum_r j^r [u^(r+1)]E = (E(u=j) - 1)/j, and x^(-j) E(u=j) = Psi_j. This gives (5.1).
  - The identity (j/alpha) S(j hbar d)/S(alpha hbar d) = sinh(j hbar d/2)/sinh(alpha hbar d/2) gives the second form of Psi_j.
- (iii) Expand at v = oo, with w = 1/v.
  - The exponent of Psi_j is -Sum_a (j/alpha_a) log(v-a) + Sum_{k>=1} (hbar w)^(2k) [A_k + B_k (1+w)^(-2k)].
  - Here A_k = (2k-1)! j tau c_k(j, 1/tau) and B_k = -(2k-1)! j(tau+1) c_k(j, -1/(tau+1)).
  - So [hbar^(2g)] Psi_j = w^(2g-j)(1+w)^M Sum_p A_{g-p} B_p (1+w)^(-2p), with A_m := [t^(2m)] exp(Sum A_k t^(2k)) and B_p := [t^(2p)] exp(Sum B_k t^(2k)).
  - Taking the residue gives (5.2) with these A and B.
- (iv) Norlund identities, which put A and B in the stated form.
  - For integer N, sinh(Nu)/sinh(u) = Sum_{k<N} e^((2k-N+1)u) and e^(c d) log v = log(v+c). These give exp(Sum A_k t^(2k)) = prod_k (1+s_k t)^(-1), exp(Sum B_k t^(2k)) = prod_k (1+s'_k t), and the product form of Psi_j.
  - prod_{k<N}(1+s_k t) = Sum_n (N)_n [y^n] S1(y)^(N+1) (t/tau)^n. This is B_N^(N+1)(x) = (x-1)...(x-N) (proof: N! Res_t e^(xt) dt/(e^t-1)^(N+1) = N! Res_u (1+u)^(x-1) du/u^(N+1), u = e^t - 1), together with the Appell property and B_n^(N+1)((N+1)/2) = n! [y^n] S1(y)^(N+1).
  - The inverse is prod(1+s_k t)^(-1) = Sum_n (-N)_n [y^n] S1(y)^(1-N) (t/tau)^n. This is the Tricomi-Erdelyi expansion Gamma(Z+a)/Gamma(Z+b) ~ Z^(a-b) Sum_n C(a-b,n) B_n^(a-b+1)(a) Z^(-n) (DLMF 5.11.17). Equivalently: both sides have coefficients polynomial in the exponent and agree at positive integers.
  - For fixed j, both sides of each identity are Laurent polynomials in tau. They agree for the infinitely many tau with j tau in Z_{>0}, hence for all tau. QED (uniform in g).

**Corollary 4 ((6.1) for all h; uses P2).**
- By P2, sigma^g Xi_g(j) = sigma^(2g-1) [x^j] F_{g,1}/C(M-1, j-1).
- Three elementary facts: (M)_{2p}(M-2p)_{j+1-2g} = (M)_{j+1-2g+2p}; (-N)_{2g-2p} = N(N+1)...(N+2g-2p-1); M - j = N. They give
  (-N)_{2g-2p} (M)_{2p} C(M-2p, j+1-2g)/C(M-1, j-1) = M N (j-1)!/(j+1-2g)!.
- Since MN/j^2 = sigma, the p-th term of (5.2) becomes
  prod_{l=1}^{2g-2}(j-l) tau^(2p) [y^(2p)] S1(y)^(1+M) (1+tau)^(2g-2p) [y^(2g-2p)] S1(y)^(1-N).
- Summing over p:
  sigma^g Xi_g(j) = prod_{l=1}^{2g-2}(j-l) [lambda^(2g)] S1(tau lambda)^(1+M) S1((1+tau) lambda)^(1-N).
- This is (6.1): e1 = -tau, e2 = 1+tau, 1 + j e2 = 1+M, 1 + j e1 = 1-N, and S1 is even. It holds for all g >= 1 and j >= 1, as rational functions of tau.
- Hand checks:
  - g = 1 gives sigma Xi_1(j) = ((j-2)sigma - 1)/24 (referee A).
  - g = 2 gives Xi_2(j) = (j-1)(j-2)[(5j^2-26j+24) + (28-12j)/sigma + 7/sigma^2]/5760 (proof_M_hurwitz.md s.4).
- P2 cross-check for one leaf: Route E verified (6.1) identically against its Hodge engine for h <= 7 (brief s.2(g)). The TR side equals (6.1) for all h, without P2. So P2 for r = 1 holds identically in tau for h <= 7, mod [ABDKS].

**Two-point function.** The same integration by parts in both variables applies; for n >= 2, Sum_r j^r [u^r] f = f(j). It gives

    [x_1^j1 x_2^j2] F_{g,2} = (1/(j1 j2)) [hbar^(2g)] Res_{v1=oo} Res_{v2=oo} Psi_{j1}(v1) Psi_{j2}(v2) dv1 dv2/((v1-v2)^2 - (j1+j2)^2 hbar^2/4),   (5.4)

where the inner residue is taken with |v2| > |v1|, so it encloses the diagonal poles. (5.4) agrees exactly with TR (s.7).

The general-n analogue is
[x^mu] F_{g,n} = -(prod 1/mu_i) [hbar^(2g-2+n)] Res...Res prod_i Psi_{mu_i}(v_i) Sum_{sigma n-cycle} prod_i dv_i/(v_i - v_sigma(i) + (mu_i + mu_sigma(i)) hbar/2),
with iterated residues taken in the order |v_n| > ... > |v_1|. It is derived the same way and checked here only for n <= 2.

## 6. The four lines of the brief

### Line (1): Lemma D at integer tau >= 2 from the global loop equations

Proposition 5 (proved). For g >= 2 let
- Brk_g(v,v') = omega_{g-1,2}(v,v') + Sum_{h=1}^{g-1} omega_{h,1}(v) omega_{g-h,1}(v'), and
- Phi_g = Brk_g(v,v)/(2 dX^2) - v omega_{g,1}(v)/dX.

Then Sum_{x(v')=x(v)} Phi_g(v') = 0.

Proof, for integer tau, where x has degree tau+1.
- Let Q(x) = Sum_{i<j} [Brk_g(v_i,v_j) + omega_{0,1}(v_i) omega_{g,1}(v_j) + omega_{0,1}(v_j) omega_{g,1}(v_i)]/dX^2. It is a rational function of x.
- Q is holomorphic at x(tau): for the pair of sheets meeting at the branch point by the quadratic loop equation (with its double zero), and for the mixed pairs by the linear loop equation.
- Q is holomorphic at x = 0 and vanishes at x = oo. Hence Q = 0.
- The global linear loop equation Sum_i omega(v_i) = 0, in each slot, gives Sum_{i<j} Brk(v_i,v_j) = -(1/2) Sum_i Brk(v_i,v_i). The omega_{0,1} terms become Sum_i v_i omega_{g,1}(v_i) dX. So Q = -Tr Phi_g.
- For any tau the identity is equivalent to a local statement: the sigma-even polar part of Phi_g at v = tau vanishes (checked in trace_qle.pl).

Consequence at tau = 1.
- The two sheets are v and 1/v, and anti-invariance gives Tr Phi_g = Brk(v,v)/dX^2 - (v - 1/v) omega_{g,1}/dX.
- Hence x dF_{g,1}/dx = t Brk_g(v,v)/((1-t^2) dX^2) at v = 1/t. This is O(x^(2g-1)) given the bound in lower genera.
- So Theorem 3.1 for n = 1 follows in one line. At tau = 1, TR is the exact global identity omega_{g,1} = Brk(v,1/v)/((y(v) - y(1/v)) dX).

The obstruction for tau >= 2, stated precisely.
- Phi_g(v_oo(x)) = -Sum_{l=1}^{tau} Phi_g(v_l(x)). So Lemma D(g,1) is equivalent to the vanishing, up to order (2g-3)tau, of the tau-section of the v = 0 expansion of Phi_g.
- That section involves v omega_{g,1}/dX at v = 0, where multiplication by v mixes all residue classes mod tau. It also involves products of lower correlators at v = 0.
- At v = 0 the inductive hypothesis gives nothing: the correlators have degree 1 there.
- The higher global loop equations (orders <= tau+1) bring in the same v = 0 data.
- The x-y swap removes this obstruction: it writes every correlator through the log points v = 0, -1 themselves.

### Line (2): other tau with extra symmetry

- In (2.1), tau enters only through the vital data (0, 1/tau) and (-1, -1/(tau+1)).
- At tau = -1/2 the two residues coincide (alpha = -2). The involution v -> -1-v exchanges the two vital points and preserves (2.1). This is the TR explanation of the symmetry behind Route A's Theorem C.
- Even parts are covered as well. P2's normalisation C(mu/2 - 1, mu - 1) vanishes at tau = -1/2 for even mu, so TR at tau = -1/2 alone cannot see them. But Theorem 1 at generic tau, together with polynomiality in 1/sigma, gives Xi_h(R; -1/2) = 0 for every R in the range, even parts included.
- The same holds at tau = -1/3 (Theorem C') and at every special framing.
- Route A's divisibility zeros are the zeros of the factor C(j(tau+1)-1, j-1) relating [x^j]F to Xi.

### Line (3): a degree argument valid for all tau at once

This is Theorem 1. The correlators are written through the other two log points v = 0, -1 of X. Each hbar-order brings a factor (v-a)^(-1), which vanishes at v = oo and blows up at v = a.

### Line (4): the x-y swap

It gives closed forms at the dilaton point. C_tau is the base-point curve q_0, and (5.1)-(5.4) are its closed forms. Lemma D and (6.1) follow.

## 7. Checks

All checks are exact. Scripts are in TR2_recursion and outputs in TR2_recursion\out.

- **(2.1), n = 1, against EO recursion, as exact rationals.** TR is computed mod 3 primes with CRT and rational reconstruction, and also compared modulo each prime.
  - cmp_n1.pl; out/cmp_n1_tau2_g4.txt, cmp_n1_tau3_g5.txt, cmp_n1_tau3_7_g5.txt, cmp_n1_taum5_3_g5.txt.
  - tau = 2: g = 2..4, mu <= 2g+2, 24 values.
  - tau = 3, 3/7, -5/3: g = 2..5, mu <= 2g+2, 36 values each.
  - 0 mismatches.
  - Example values: Xi_5(9) = -558409145/13759414272 (tau = 3) and 7833182869/1343692800000 (tau = 3/7).
  - The tau = 1, g = 6 comparison hit its 50-minute timeout and is not included.
- **(2.1) against the referee's independent exact TR at tau = 1.**
  - out/tau1_values.pl.
  - Xi_5(8) = 0, Xi_5(9) = -2205/1048576, Xi_6(10) = 0, Xi_6(11) = -297675/33554432, in the TR normalisation sigma^(1-g)Xi.
  - These equal the values in referee_hurwitz_1.md s.2(d) and proof_status 1.8.
- **(2.1) pointwise at random points mod p.** The formula is evaluated as Taylor series and compared with TR in the e_k basis; 21 + 1 points, all agree.
  - xy_point.pl; out/xy_point_all.txt, out/xy_point_n2_g5.txt.
  - n = 1: tau = 2, g = 1..3; tau = -1/2, -2/7, 7, 1/5, g = 4.
  - n = 2: tau = 3/7, g = 1..3 and g = 5; tau = -5/3 and 2, g = 4; tau = 1/2, g = 2; tau = -1/2, g = 3.
  - n = 3: tau = -5/3, g = 0..2; tau = 3/7, g = 3.
  - n = 4: tau = 2, g = 0, 1; tau = 3/7, g = 1.
- **Lemma D directly in TR**, for all mu with |mu| <= 2g-2+n.
  - lemmaD_check.pl; out/lemmaD_tau2_g5.txt, lemmaD_tau1_2_g5.txt.
  - tau = 3/7: g = 2..4, n = 1..3, 82 coefficients.
  - tau = 2 and tau = 1/2: g = 2..5, n = 1, 2, 78 coefficients each.
  - 0 violations. The first layer |mu| = 2g-2+n is non-zero, except that at tau = 1/2 the all-even first-layer entries vanish.
- **(5.2) against the full formula (2.1).**
  - onepoint.pl; out/onepoint_*.txt.
  - tau = 2: g <= 5; tau = 3/7 and -5/3: g <= 4.
  - 0 mismatches.
- **(5.2) against (6.1).**
  - closed_vs_61.pl; out/closed_vs_61_*.txt.
  - tau = 2: g <= 9, j <= 3g+1, 140 values.
  - tau = 5/3: g <= 8, 112 values.
  - tau = -2/7: g <= 6, 58 values. The 7 points where C(j(tau+1)-1, j-1) = 0 are skipped.
  - 0 mismatches.
- **Two-point (5.4) against TR**, compared modulo each of 3 primes.
  - two_point.pl; out/two_point_*.txt.
  - tau = 2: g <= 3, j <= 6, 63 coefficients.
  - tau = 3/7: g <= 3, j <= 5, 45 coefficients.
  - tau = -5/3: g <= 4, j <= 4, 40 coefficients.
  - 0 mismatches.
- **Norlund identities at non-integer tau.**
  - norlund_check.pl; out/norlund_check.txt.
  - tau = 3/7 (84 checks), -13/5 (70), 2 (70).
  - 0 mismatches.
- **Trace-form quadratic loop equation**, i.e. the sigma-even polar part of Phi_g vanishes.
  - trace_qle.pl; out/trace_qle.txt.
  - tau = 2: g = 2..4; tau = 3/7: g = 2, 3; tau = 1: g = 2..4.
  - Holds in every case.
- **Engine validation.**
  - t_basic.pl, t_mod.pl.
  - tau = 2: omega_{0,3} = -6, omega_{1,1} = (163/1728)e_0 - (1/4)e_1, Xi_2(3) = -73/20736 and Xi_3(5) = -245/82944, the values in hvc_tr_tau2_g3n3.txt.
  - tau = 1: Xi_2(3) = -1/256, Xi_3(5) = -3/1024, and Xi_4(7..9) = 8 x the referee's values.

How this goes beyond the earlier range:
- Before this round, exact TR had been compared with Lemma D at tau = 1 and -1/2 up to genus 6 (n = 1), and at tau = 2, 3 up to genus 3 (n <= 3).
- New here: genus 5 at tau = 2, 3, 1/2, 3/7, -5/3 (n <= 2); genus 4 at seven further tau; n up to 4; and the closed formulas, (5.2) = (6.1), up to genus 9.

## 8. Files (TR2_recursion)

Core libraries:
- TRtau.pm: exact (Math::BigRat) EO recursion on C_tau for rational tau, in the local coordinate z at v = tau with basis e_k. It is slow and was used for validation.
- TRmod.pm: the same mod p, fast. Functions: new_curve(tn, td, kmax, p), run_tr, xcoef, xi_tr.
- ratrec.pm: primes, CRT, rational reconstruction.
- XY.pm: formula (2.1) for n = 1, expanded at v = oo, exact.

Check scripts:
- cmp_n1.pl: (2.1) against TR for n = 1.
- xy_point.pl: pointwise mod-p test of (2.1) for any n, over all n-cycles.
- onepoint.pl and closed_vs_61.pl: (5.2) and (6.1).
- two_point.pl: (5.4) against TR.
- norlund_check.pl, lemmaD_check.pl, trace_qle.pl.
- t_basic.pl, t_mod.pl: engine tests.

Other:
- xy_n1.pl: first version, written with y = -v and so carrying the opposite overall sign. Superseded by XY.pm.
- tr2_closed_check.wls: requested Mathematica check, not run.
- out/: all outputs. out/rerun_quick.sh re-runs the quick checks.

Usage examples:
- perl cmp_n1.pl 3 7 5 4 3
- perl xy_point.pl -5 3 2 3 7
- perl two_point.pl 2 1 3 6
- perl closed_vs_61.pl 2 9
- perl lemmaD_check.pl 3 7 4 3
- perl onepoint.pl 2 5 3

## 9. Dependencies, and what is not done

- **P2 (Lemma B).** Used only in Corollaries 2 and 4. Route TR1 is writing it out.
- **[ABDKS] Thm 2.8 and Prop. 3.1.** Used for Theorems 1 and 3. This is a published theorem. I checked its hypotheses for C_tau and checked my transcription of its formula against exact TR for n <= 4, but I did not re-prove it. A referee should re-read its Lemma 5.5 at v = oo, where dX has a simple pole and dy a double pole.
- **Not used:** KP, Corollary E, Theorem D5, the parity argument, and TR on the conifold.
- **Line (1) alone does not close.** The loop-equation route does not close for tau >= 2 without the x-y swap (s.6(1)).
