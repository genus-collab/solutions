# Lemma A round 2: route R1 (saved by the main session from the workflow journal, 2026-10-31)

Subagents cannot write report files, so this file holds the agent's structured return value: the fields first, then the full report text (report_text) verbatim.

## all_genus_claim

1

## all_genus_statement

Lemma A (M_irr(h): Xi_h(R) = 0 for all R with parts >= 2 and D(R) <= 2h-3) holds for every h, hence Corollary R holds in all genera. Proof: T1 (proved from Route A Theorem D via the flow theorem and the leaf-curve identity), the special-framing lemma and the forcing theorem give simple zeros of sigma^G Xi_G(D+1) at sigma0 = -k0(q-k0)/q^2. Together with Route A's merge zeros Z(D) there are at least G-1 of them for every G >= 8 and odd D, so (OL) holds for G >= 8. Theorem F covers G = 5,6,7 and the exact data cover G <= 4. The statement is proved modulo the project's inputs listed under obstructions; no TR of any kind was used.

## claims

### [1]

#### depends_on

Route A Theorem D (MV + KP Fact 3.5 + Takasaki-Takebe) at every point of the slice; genus-0 Xi_0 = Gamma(J)/Gamma(J-r+3); leaf curve b_j (Theorem 11.1 + string)

#### evidence

proof s.2-s.7; exact symbolic in tau k <= 20 (2876 evaluations), full KP mod p k <= 23 at 4 framings, ODE to order 100 at 4 framings, independent z-expansion implementation (t_zexp.pl) at 3 framings

#### statement

T1 (closed form of the polar-order relations): if Lemma A holds in genera < G then for 1 <= k <= 2G-3, relation_k := [t^k] Sum_T prod b_t Xi_G(T)/|Aut T| = kappa_k Sum_{L<=k} y_L (-1)^(k-L)[(1+tau)binom((L+1)tau,k-L) - tau binom((L+1)(1+tau)+k-L-1,k-L)], y_L = Pi_{L+1} Xi_G(L+1); equivalently Sum_k t^k relation_k = Res Psi(lambda)[(1+tau)Kappa(t lambda) - tau Kappa(t(lambda-1))].

#### status

proved-modulo-literature

### [2]

#### depends_on

Route A Theorem D; branch-point identity; genus-0 self-similarity

#### evidence

proof s.5-s.6; ode_path.pl/ode_path2.pl along random paths (exact, 0 mismatches; perturbation detected)

#### statement

Flow theorem: along any path of the slice P_1 = 0 the genus-G one-point function satisfies hatU_P = [Psi o phi_P]_- with phi_P(lambda) = W lambda - [T_P(W/u(lambda))]_+ + [T_P(W/u(lambda))]_+(0); and F_G(P) - F_G(0) = Res_mu Psi(mu) log(W u(mu)/u(phi_P^-1(mu))) for every P (in range).

#### status

proved-modulo-literature

### [3]

#### evidence

residue proof s.4; t_generator.pl part 1

#### statement

Branch-point identity: d/dlambda [u^-j]_+ at lambda = 1+tau equals j Pi_j for all j; hence N(1+tau) = 0 on every path.

#### status

proved

### [4]

#### evidence

Lagrange inversion + Pfaff + polynomiality argument (s.7); t_phi_explicit.pl (5 framings, order 6), t_res_ut.pl, t_crit.pl

#### statement

Curve identity on the leaf curve: C u_t(phi_t(lambda)) = u(lambda) with u_t(mu) = t(e^(t(mu-1))-1)^tau (e^(t mu)-1)^(-1-tau) e^(-tau t(mu+tau)) and C = alpha e^(tau^2 t)/H_1(1-e^-t); via Res_mu u_t^n = (-1)^(n-1) t^(n-1) e^(-n tau^2 t) Pi_n 2F1(-n tau,-n tau;n;1-e^-t).

#### status

proved

### [5]

#### depends_on

Route A Theorem D with R empty (two-leaf closed form)

#### evidence

proof s.8; t_special.pl all framings G = 3..9

#### statement

Special-framing lemma: at tau0 = -k0/q (gcd 1, k0 < q) first-order two-leaf regularity (m+n <= 2G-1) is equivalent to (1-y)^(k0-1)(1-(1+tau0)y)F(y) having no coefficients y^N, q <= N <= 2G-2.

#### status

proved-modulo-literature

### [6]

#### depends_on

T1, special-framing lemma, Route D polar relations (I1)

#### evidence

odd functional equation argument s.9; t_forcing.pl: 1131 tuples G <= 20, 0 failures

#### statement

Forcing theorem: under the hypotheses of (OL), Xi_G(D+1) vanishes at tau0 = -k0/q whenever gcd(k0,q) = 1, k0 <= D+1, 2k0 < q, q >= D+2, 3(q-k0) < 2G.

#### status

proved-modulo-literature

### [7]

#### evidence

exact integer count G <= 1000 (count_zeros2.pl, cross-checked with count_zeros.pl); analytic bounds (B1), (B2) for G > 300 (s.10)

#### statement

Zero count: Z1(D) + Z2(G,D) >= G-1 for all G >= 8 and odd 3 <= D <= 2G-3 (fails only at G = 4..7, D = 3).

#### status

proved

### [8]

#### depends_on

T1, forcing theorem, zero count, Route A Corollary E (Z1 zeros), brief 2(a),(b)

#### evidence

s.11

#### statement

(OL) holds for every G >= 8.

#### status

proved-modulo-literature

### [9]

#### depends_on

(OL) for G >= 8; Route A Theorem D, Corollary E, Theorem F (G = 5,6,7); Route D (H) = Brini-Schuler, Lemma 5.3, (V); leaf_structure_proof Theorem 11.1; Theorem A; brief 2(a),(b); exact data h <= 4

#### evidence

s.11

#### statement

Lemma A (M_irr(h)) holds for every h; hence Corollary R (all-genus regularity of the refined resolved-conifold Hodge-vertex graph sum).

#### status

proved-modulo-literature

### [10]

#### evidence

global_suite.pl, mod p = 2^31-1, G = 3..10

#### statement

Global nullities (window G-1): TLR 2,3 at G = 3,4 and 0 for G = 5..10; TLR+CF 1,2 then 0; full 1,0 then 0; full+CF 0 for all G = 3..10. Single-level problems at odd D have nullity 0 for G = 5..10 with TLR alone.

#### status

computer-verified

### [11]

#### evidence

scan_points.pl, out/scan_G12_two_cf.txt

#### statement

Pointwise scan (two-leaf + closed-form rows from the full KP at all special framings) closes (OL) for 5 <= G <= 12 (with Route A Z(5) = 6 at (5,5)).

#### status

computer-verified

### [12]

#### evidence

t_flow_curve.pl (fails at t^2), t_ansatz2.pl (two-leaf test fails at J = 5), t_zeta.pl

#### statement

u_t(phi_t) = u/alpha (constant alpha); deformed-curve ansatz Y_P(z) - Sum P_j z^j = lambda + tau; zeta = (e^t-1)e^(tau t)((1+tau)e^t-tau)^tau.

#### status

refuted


## key_formulas

- u(lambda) = (lambda-1)^tau lambda^(-1-tau),  h = u/u' = -lambda(lambda-1)/(lambda-1-tau),  Psi(lambda) = Sum_c Pi_c Xi_G(c) u(lambda)^c
- genus 0 on the slice: Y_P(z) = z - sigma Sum_n Pi_n W^(n+1) z^-n/(n+1) = W Y_0(z/W),  W = 1 + Sum_j Pi_j P_j W^j,  Y_0 = lambda + tau at 1/z = u(lambda)
- flow: phi_P(lambda) = W lambda - [T_P(W/u(lambda))]_+(lambda) + [T_P(W/u(lambda))]_+(0),  hatU_P = [Psi o phi_P]_-,  F_G(P) - F_G(0) = Res_mu Psi(mu) log(W u(mu)/u(phi_P^-1(mu)))
- branch point: d/dlambda [u^-j]_+ (1+tau) = j Pi_j
- leaf curve: P_j(t) = (-1)^j (Pi_j/j) t^(j-1) H_j(T)/H_1(T)^j,  H_j = 2F1(-j tau,-j tau;j;T),  T = 1-e^-t
- curve identity: alpha e^(tau^2 t) H_1^-1 u_t(phi_t(lambda)) = u(lambda),  u_t(mu) = t (e^(t(mu-1))-1)^tau (e^(t mu)-1)^(-1-tau) e^(-tau t(mu+tau))
- residue: Res_mu u_t(mu)^n dmu = (-1)^(n-1) t^(n-1) e^(-n tau^2 t) Pi_n 2F1(-n tau,-n tau;n;1-e^-t)
- T1: Sum_k t^k relation_k = Res_lambda Psi(lambda)[(1+tau)Kappa(t lambda) - tau Kappa(t(lambda-1))],  Kappa(x) = tau x + log((e^x-1)/x);  beta_{k,L} = kappa_k (-1)^(k-L)[(1+tau)binom((L+1)tau,k-L) - tau binom((L+1)(1+tau)+k-L-1,k-L)]
- s-form: Sum_k (k+1) R_k s^(k+1) = (1+tau)(1-(1+tau)s)F(s)/(1-s) - tau(1-tau s)F(s/(1+s))/(1+s),  F(y) = Sum_c c y_{c-1} u(y)^c
- special framing tau0 = -k0/q: two-leaf regularity <=> deg[(1-y)^(k0-1)(1-(1+tau0)y)F(y)] <= q-1 mod y^(2G-1)
- forcing: X(lambda) = A r(lambda) + B r(lambda+1) odd, r = pi/(lambda^A (lambda-1)^B), A = q-k0, B = k0 => pi = 0 when k0 <= D+1, 3(q-k0) < 2G
- zero count: Z1(D) + Z2(G,D) >= G-1 for G >= 8; Z2 >= (16/9)G - (8/3)(D+1) - 5/3; Z1(D) >= 0.0376 D^2 + 0.5D - 2.75

## next_steps

- Critic/referee pass on s.2-s.11 of the report, in particular the level bookkeeping of Lemma 2.1, the Toeplitz lemma, the flow theorem (s.5), step (d) of the curve identity, Lemma 8.1 and the degree count of the forcing theorem.
- Main session: run R1_checks.wls (independent Mathematica checks C1-C6) and record the output.
- Write the T1 proof, the flow theorem and the forcing/zero-count closing of (OL) into the project notes (leaf_proof / conifold note) as the all-genus completion of Corollary R.
- Optional: use the closed formula F_G(P) - F_G(0) = Res Psi log(W u/u o phi_P^-1) to give explicit multi-leaf merge formulas (generalising Corollary E) and a proof of (OL) for G = 5,6,7 without Theorem F (the pointwise scan suggests full regularity at all special framings suffices).

## obstructions

- No open mathematical obstruction on the main line. The result is conditional on the project's inputs: Route A Theorem D (KP of the MV generating function + Takasaki-Takebe, used in generating form at every point of the slice and in range), Corollary E, Theorem F (only for G = 5,6,7); Route D's single-vertex polar relation ((H) = Brini-Schuler Prop. 3.1, Lemma 5.3, (V)); leaf_structure_proof Theorem 11.1; brief 2(a),(b); Theorem A; exact data h <= 4.
- The uniform zero count fails for G = 4..7 at D = 3 (Z1 + Z2 = 2, 3, 3, 5 < G-1), so the argument needs Theorem F / the data below genus 8; TLR(G) itself is still unproved for all G (no longer needed).

## requested_checks

- Run: wolframscript -file "G:/My Drive/Projects/Old/refGW/hodge_vertex_conifold/lemmaA/round2/R1_rigidity/R1_checks.wls" (expected: (C1) True, (C2) True, (C3) {True, True} three times, (C4) lists of True at tau = 2/7, -3/5, 3, (C5) {tuples, 0}, (C6) failures = 0 for 8 <= G <= 60 and 4 for G = 4..7). (C3) symbolic tau and (C4) are the heavy parts.
- Critic: referee the proof chain s.2-s.11 (Lemma 2.1 level bookkeeping with Route A Theorem D at every point of the slice; Toeplitz Lemma 3.1; flow theorem 5.1; curve identity 7.1 step (d); Lemma 8.1; Theorem 9.1; Proposition 10.1 and Route A's Z(D) multiplicities).
- Main session: confirm that relation_k as defined here (b_j = A_j alpha^-j/c^2, string elimination of lambda_1 leaves; t_leafcurve.pl part 2) is, up to a unit prefactor, the single-vertex polar coefficient whose vanishing Route D proves from (H) and Lemma 5.3.

## route

R1 (rigidity of the KP system), round 2: exact integration of the genus-G linearised KP system along the leaf curve (flow theorem), closed form T1 of the polar-order relations, special-framing forcing, zero count

## summary

The closed form T1 of the polar-order relations is now proved, not just observed (earlier this round it had only been checked, up to order 100). The proof uses only the KP input (Route A Theorem D at every point of the slice P_1 = 0), Lagrange inversion, formal residues and one Pfaff transformation; no TR of any kind. Two structural results on the whole slice carry it. First, along any path the genus-G one-point function in the uniformising coordinate lambda (u = (lambda-1)^tau lambda^(-1-tau)) moves by composition with the explicit polynomial map phi_P(lambda) = W lambda - [T_P(W/u(lambda))]_+ + [T_P(W/u(lambda))]_+(0), where T_P(z) = Sum_j P_j z^j and W = 1 + Sum_j Pi_j P_j W^j (the flow theorem). Second, this gives the closed formula F_G(P) - F_G(0) = Res Psi(mu) log(W u(mu)/u(phi_P^-1(mu))). On the leaf curve the map phi_t links u to the explicit curve u_t(mu) = t(e^(t(mu-1))-1)^tau (e^(t mu)-1)^(-1-tau) e^(-tau t(mu+tau)) up to the constant C = alpha e^(tau^2 t)/H_1. The proof of that step is Lagrange inversion plus the residue identity Res u_t^n = (-1)^(n-1) t^(n-1) e^(-n tau^2 t) Pi_n 2F1(-n tau,-n tau;n;1-e^-t), obtained by a Pfaff transformation. Putting these together gives Sum_k t^k relation_k = Res Psi(lambda)[(1+tau)Kappa(t lambda) - tau Kappa(t(lambda-1))], which is T1.

From T1, a special-framing lemma (two-leaf regularity at tau0 = -k0/q is equivalent to (1-y)^(k0-1)(1-(1+tau0)y)F(y) being a polynomial of degree <= q-1 mod y^(2G-1)) and an exact odd-functional-equation argument give the forcing theorem. Under the hypotheses of (OL), Xi_G(D+1) vanishes at sigma0 = -k0(q-k0)/q^2 whenever gcd(k0,q) = 1, k0 <= D+1, 2k0 < q, q >= D+2 and 3(q-k0) < 2G. These simple zeros, added to Route A's merge-formula zeros Z(D), give at least G-1 zeros of sigma^G Xi_G(D+1) (which lies in span(sigma..sigma^(G-1))) for every G >= 8 and every odd D. This is an exact integer count for G <= 1000 plus an explicit analytic bound for G > 300. So (OL) holds for all G >= 8. With Theorem F for G = 5,6,7 and the exact data for G <= 4, Lemma A holds in all genera, and with it Corollary R. This is conditional only on the project's established inputs: Route A Theorem D, Corollary E and Theorem F; Route D's single-vertex polar relation, i.e. (H) = Brini-Schuler plus Lemma 5.3 and (V); the leaf closed form (Theorem 11.1); the lambda_g and Hurwitz components; Theorem A. It stops short of an unconditional statement only because those inputs, and this proof, still need a referee. Every new lemma also has an independent exact check, and a Mathematica script (R1_checks.wls) is ready for the main session to run.

## report_text

# Route R1 (rigidity of the KP system), round 2: report

Work folder: G:/My Drive/Projects/Old/refGW/hodge_vertex_conifold/lemmaA/round2/R1_rigidity/ (all scripts and outputs; out/ holds outputs).
NO TOPOLOGICAL RECURSION of any kind was used: no spectral-curve recursion, no loop equations, no Virasoro constraints as loop equations, no BKMP/remodelling, no dictionary P2, no TR script, not even as a cross-check. Inputs are the KP hierarchy of the Marino-Vafa generating function through Route A's Theorem D (allowed by the user ruling of 2026-10-31), Lagrange inversion, formal residues, one Pfaff transformation of 2F1, exact algebra, and the project inputs listed in s.1. Exact checks were done in Perl (Math::BigRat) or with exact modular arithmetic; no wolframscript was run.

## 0. Main result and chain of the proof

THEOREM (all genera; proved modulo the inputs of s.1). Lemma A = M_irr(h) holds for every h >= 1. Hence Corollary R (the all-genus regularity of the refined resolved-conifold Hodge-vertex graph sum) holds.

1. T1, the closed form of the polar-order relations (s.2-s.7; proved modulo Route A Thm D). On the leaf curve P(t),
   Sum_k t^k relation_k = Res_lambda Psi(lambda) [ (1+tau) Kappa(t lambda) - tau Kappa(t(lambda-1)) ],   Kappa(x) = tau x + log((e^x-1)/x),
   equivalently relation_k = kappa_k Sum_{L<=k} y_L (-1)^(k-L) [ (1+tau) binom((L+1)tau, k-L) - tau binom((L+1)(1+tau)+k-L-1, k-L) ].
   It is obtained by an exact integration of the genus-G linearised KP system, using a FLOW THEOREM (s.5) that holds on the whole slice and a curve identity (s.7) specific to the leaf curve.
2. Special-framing lemma (s.8): at tau0 = -k0/q, two-leaf regularity is equivalent to (1-y)^(k0-1)(1-(1+tau0)y)F(y) being of degree <= q-1 mod y^(2G-1).
3. Forcing theorem (s.9): under the hypotheses of (OL), Xi_G(D+1)(tau0) = 0 for tau0 = -k0/q with gcd(k0,q) = 1, k0 <= D+1, 2k0 < q, q >= D+2, 3(q-k0) < 2G.
4. Zero count (s.10): Z1(D) + Z2(G,D) >= G-1 for every G >= 8 and every odd 3 <= D <= 2G-3.
5. So (OL) holds for all G >= 8 (s.11). Theorem F covers G = 5,6,7 and the exact data cover G <= 4, which gives Lemma A for all h.

Side result (s.6): a closed formula for the genus-G free energy on the whole slice P_1 = 0 in terms of the one-leaf values,
   F_G(P) - F_G(0) = Res_mu Psi(mu) log( W u(mu) / u(phi_P^-1(mu)) ),   phi_P(lambda) = W lambda - [T_P(W/u(lambda))]_+(lambda) + [T_P(W/u(lambda))]_+(0).

## 1. Conventions and inputs

- c = 1, framing tau, sigma = tau(1+tau), Pi_j = prod_{a=1}^{j-1}(j tau+a)/(j-1)! = binom(j(1+tau)-1, j-1), kappa_k = [x^k](tau x + log((e^x-1)/x)): kappa_1 = tau+1/2, kappa_k = B_k/(k k!) for k >= 2, and kappa_k = 0 for odd k >= 3.
- MV normalisation: Y(T) = Pi_T Xi_G(T), y_L := Y(L+1) (L = level). Unit parts are handled by string, Y(T+{1}) = |T| Y(T). On the slice (coordinates P_j, j >= 2), F_G(P) = Sum_T Y(T) P^T/|Aut T| and d_1 F_G = E F_G with E = Sum j P_j d_j.
- Leaf curve: P_j(t) = b_j/Pi_j = (-1)^j (Pi_j/j) t^(j-1) H_j(T)/H_1(T)^j, with H_j = 2F1(-j tau,-j tau; j; T) and T = 1 - e^-t. Here b_j = A_j alpha^-j/c^2, where A_j is given by leaf_structure_proof.md Theorem 11.1 and the lambda_1-leaves are eliminated by string; this was checked in t_leafcurve.pl (216 comparisons, 0 mismatches).
- Polar-order relations: relation_k := [t^k] Sum_T (1/|Aut T|) prod_{t in T} b_t(t) Xi_G(T) = [t^k] F_G(P(t)).
- INPUT (I1) (Route D: (H) = Brini-Schuler Prop. 3.1, Lemma 5.3, (V); leaf closed forms): relation_k = 0 for 1 <= k <= 2G-3, given Lemma A in genera < G. Any unit prefactor of the single-vertex graph is irrelevant.
- INPUT (I3) (Route A Theorem D for all R = Taylor coefficients at P of the linearised dispersionless Hirota identity; Fact 3.5 + Takasaki-Takebe; given Lemma A in genera < G; valid for coefficients of level <= 2G-3, i.e. excess >= 1): at every point P of the slice,
  Sum_{m,n>=1} G_mn z^-m w^-n = (V_P(z)-V_P(w))/(Y_P(w)-Y_P(z)),   G_mn = d_m d_n F_G,   V_P(z) = Sum_n d_1 d_n F_G z^-n.
- INPUT (I4) (genus 0; consequence of Xi_0(T) = Gamma(J)/Gamma(J-r+3), used by Route A): Y_P(z) = z - sigma Sum_n Pi_n W^(n+1) z^-n/(n+1) = W Y_0(z/W), where W = W(P) solves W = 1 + Sum_{j>=2} Pi_j P_j W^j. Proof: d_1 d_n F_0 = sigma Pi_n Sum_S Xi_0(S+{1,n}) prod_s Pi_s P_s/|Aut S|, with Xi_0(S+{1,n}) = (s+n)!/(s+n+1-m) (s = |S|, m = #S). The plane-forest count [prod q_j^(m_j)] W^a = a (s+a-1)!/((s-m+a)! prod m_j!) gives W^(n+1)/(n+1). Also Y_0(z) = lambda + tau when 1/z = u(lambda), by Lagrange: [u^n](1/y) = -sigma Pi_n/(n+1).
- Other inputs: the lambda_g and Hurwitz components (brief 2(a),(b): sigma^G Xi_G(R) is in span(sigma..sigma^(G-1)) in the range); Route A Corollary E (merge formula) and its zero count Z(D); Theorem A (Xi_h(2) = 0); Theorem F (5 <= G <= 20); exact data for h <= 4.
- Uniformising coordinate: u(lambda) = (lambda-1)^tau lambda^(-1-tau) = lambda^-1 (1-1/lambda)^tau, y = 1/lambda, h := u/u' = -lambda(lambda-1)/(lambda-1-tau) = -(lambda+tau) - sigma/(lambda-1-tau). [f]_+ and [f]_- are the polynomial part (lambda^(>=0)) and the negative part of a Laurent series at lambda = infinity; Res f := [lambda^-1] f. The one-leaf generating function is Psi(lambda) := Sum_{c>=2} y_{c-1} u(lambda)^c, which is O(lambda^-2).

## 2. The linear system along a path

Let P(t) be any path in the slice with P_j(t) = O(t^(j-1)), T_t(z) := Sum_j P_j(t) z^j and W(t) := W(P(t)) (= alpha(t) on the leaf curve). The formal system, in power series of t, is:
 (S) V_n = n U_n + Sum_{j>=2} j P_j G_jn   [string d_1 F = E F, differentiated];
 (H) Sum G_mn z^-m w^-n = (V(z)-V(w))/(Y_{P(t)}(w) - Y_{P(t)}(z))   [I3 at P(t)];
 (E) U_n' = Sum_{i>=2} P_i' G_in,  Phi' = Sum_{i>=2} P_i' U_i,  U_n(0) = y_{n-1} (n>=2),  U_1(0) = 0,  Phi(0) = 0   [chain rule].
LEMMA 2.1. The system has a unique formal solution. If Lemma A holds in genera < G, then U_n(t) = d_n F_G(P(t)), G_mn(t) = d_m d_n F_G(P(t)) and Phi(t) = F_G(P(t)) - F_G(0) agree with that solution on every t-coefficient of level <= 2G-3, where [t^m]U_n has level m+n-1, [t^r]G_ab has level r+a+b-2 and [t^k]Phi has level k. In particular relation_k = [t^k]Phi for k <= 2G-3.
Proof. Every equation relates coefficients of equal or lower level, since P_j has t-order >= j-1. (H) holds for every coefficient of level <= 2G-3 (excess homogeneity of Theorem D), and (S), (E) hold exactly. The claim follows by induction on the t-order. QED
[ode_leaf.pl integrates exactly this system and reproduces the full KP propagation of kp_point.pl/kpmod.pm, which implements Theorem D with all R: t_T1_modp.pl, order 23, 4 framings.]

## 3. Toeplitz form in lambda

Put w = W/u(lambda), hatU(lambda) := U(w) and hatV(lambda) := V(w). Then Y_{P(t)} = W(lambda + tau) is affine in lambda.
LEMMA 3.1. For a polynomial f, Res_z f(z) G(z,w) dz = (1/W) [ (f(z(lambda)) z'(lambda))_+ hatV ]_-, evaluated at lambda_w.
Proof. Formal residues are invariant under z = W(lambda + tau + O(1/lambda)). The divided difference (hatV(l)-hatV(l_w))/(l_w-l) may be expanded in |l| > |l_w|, and Res_l A(l)/(l_w-l) = -A_+(l_w). So the result is -[A hatV]_+ + hatV A_+ = [A_+ hatV]_-. QED
With Pcal_t(lambda) := T_t(W/u(lambda)) = Sum_j P_j W^j u^-j, a := W - [d_lambda Pcal]_+, N := [d_t Pcal|_lambda / h]_+ + W' and b := [d_t T|_z / h]_+ = N - (W'/W) a, the system becomes
 (S') [a hatV]_- = W h d_lambda hatU;   (E') d_t hatU|_lambda = -(1/W)[N hatV]_-;   (P') Phi' = -Res_lambda b hatU.
This uses theta U = -w d_w U = h d_lambda hatU, d_t hatU|_lambda = d_t U|_w - (W'/W) h d_lambda hatU, z'/z = -1/h and d_t Pcal|_lambda = d_t T|_z - (W'/W) h d_lambda Pcal. Since a = W + O(t), the Toeplitz operator [a.]_- is invertible, so (S')(E') determine hatU uniquely from hatU_0 = Psi.

## 4. Branch-point identity

LEMMA 4.1. For every j >= 1, d_lambda [u^-j]_+ (1+tau) = j Pi_j. Consequently N(1+tau) = 0 on every path of the slice.
Proof. We have d[f]_+ = [f']_+ and, formally, [g]_+(l0) = Res g(l)/(l-l0). For g = (u^-j)' = j u^-j (l-1-tau)/(l(l-1)) and l0 = 1+tau the factor (l - l0) cancels, so [g]_+(1+tau) = j Res u^-j/(l(l-1)) = j [l^-(j-1)](1-1/l)^(-j tau-1) = j binom(j tau+j-1, j-1) = j Pi_j. Then [d_t Pcal/h]_+(1+tau) = -Sum_j (d_t P^_j)(1/j) d[u^-j]_+(1+tau) = -d_t Sum_j Pi_j P^_j = -W', by the genus-0 equation, where P^_j := P_j W^j. QED
(Checked in t_generator.pl part 1.)

## 5. Flow theorem

THEOREM 5.1. Let phi_t(lambda) := W lambda - [Pcal_t]_+(lambda) + [Pcal_t]_+(0). Then:
 (i) phi_t = lambda + O(t) is a t-series of polynomials with phi_t(0) = 0 and phi_t(1) = 1;
 (ii) d_lambda phi_t = a, and d_t phi_t = -h N = v_t d_lambda phi_t with v_t := -hN/a;
 (iii) the solution of (S')(E') is hatU_t = [ Psi o phi_t ]_-.
Proof. (i) [u^-j]_+(1) - [u^-j]_+(0) = binom(j tau+j, j) - binom(j tau+j-1, j) = Pi_j, so phi_t(1) - 1 = W - 1 - Sum_j Pi_j P^_j = 0.
(ii) Write h = -(lambda+tau) - sigma/(lambda-1-tau). Since N(1+tau) = 0, hN is a polynomial. Then [d_t Pcal]_+ = [h (d_t Pcal/h)]_+ = hN + W'(lambda+tau) - nu, where nu = [lambda^-1](d_t Pcal/h). Hence d_t phi = W' lambda - [d_t Pcal]_+ + d_t[Pcal]_+(0) = -hN + const, and the constant is 0 because both sides vanish at lambda = 0.
(iii) Put X := Psi o phi_t, R := [X]_+, hatU := [X]_-. Then d_t X = v_t d_lambda X, so d_t hatU = [v_t d_lambda X]_-. Let Vfull := (W h/a) d_lambda X and delta := W sigma c_R/(a(1+tau)(lambda-1-tau)) with c_R := d_lambda R(1+tau), and set hatV := [Vfull]_- + delta.
 - (S') holds: [a hatV]_- = [W h d X]_- + [a delta]_- = W h d hatU + W [h dR]_- + W sigma c_R/(lambda-1-tau) = W h d hatU, because [h poly]_- = -sigma poly(1+tau)/(lambda-1-tau) and [a/(lambda-1-tau)]_- = a(1+tau)/(lambda-1-tau).
 - (E') holds: -(1/W)[N hatV]_- = -(1/W)[N Vfull]_- - (1/W)[N delta]_- = [v_t d X]_- - 0. Here -(1/W) N Vfull = v_t dX exactly, and N delta is a polynomial because N(1+tau) = 0.
 Uniqueness then gives hatU_t = [Psi o phi_t]_-. QED
Nothing specific to the leaf curve is used. Checks: ode_path.pl compares the exact state with [Psi o phi]_- along random polynomial paths (tau = 2/7, -3/5, 3, -1/5; 15-28 coefficients each, 0 mismatches), while a deliberately perturbed phi gives 6 mismatches. t_generator.pl confirms that the explicit generator equals the generator of the fitted flow (tau = 1 to order 4, tau = 2 and 1/2 to order 3). The flow had first been found by fitting (fit_flow.pl; at tau = 0 phi_t = log(1+(e^t-1)lambda)/t).

## 6. Integration: the free energy on the slice

THEOREM 6.1. On any path, Phi(t) = Res_mu Psi(mu) Kflow_t(mu) dmu with Kflow_t(mu) := log W + log u(mu) - log u(psi_t(mu)) and psi_t := phi_t^-1. Kflow_t is a t-series of polynomials because psi_t fixes 0 and 1. Equivalently, the formula of s.0 holds for every P in the slice (path independence is automatic).
Proof. Substituting mu = phi_t(lambda) (formal residues are invariant under lambda -> lambda + t poly),
d/dt Res Psi Kflow = Res_lambda X (W'/W + v_t/h) phi_t' = -Res_lambda X b phi_t'/a = -Res_lambda X b = -Res_lambda b hatU = Phi'.
This uses v_t/h = -N/a, N = b + (W'/W)a, phi_t' = a and Res(bR) = 0. Both sides vanish at t = 0. QED
(ode_path2.pl: Phi from the system equals the closed formula along random paths, exact, tau = 2/7 to order 4 and tau = -3/5 to order 5. Values are nontrivial, e.g. -7340807/1372.)

## 7. The leaf curve: proof of T1

THEOREM 7.1 (curve identity FC). On the leaf curve put zeta~ := t alpha/H_1 and C := alpha e^(tau^2 t)/H_1(T). Then
 C u_t(phi_t(lambda)) = u(lambda),   u_t(mu) := t (e^(t(mu-1))-1)^tau (e^(t mu)-1)^(-1-tau) e^(-tau t(mu+tau)) = u(mu) e^(-K_t(mu)),
 K_t(mu) := (1+tau)Kappa(t mu) - tau Kappa(t(mu-1)) = Sum_k kappa_k t^k ((1+tau)mu^k - tau(mu-1)^k).
COROLLARY 7.2 (T1). Kflow_t - K_t = log C does not depend on mu. Since Psi = O(mu^-2), Phi(t) = Res Psi K_t, i.e. relation_k = kappa_k Res Psi ((1+tau)lambda^k - tau(lambda-1)^k) for k <= 2G-3, which is the beta-form of s.0. To see this, [lambda^-1] u^N lambda^k = (-1)^d binom(N tau, d), and [lambda^-1] u^N (lambda-1)^k = (-1)^d binom(N(1+tau)+d-1, d) with N = L+1, d = k-L.
Proof of FC.
(a) P^_j = P_j alpha^j = (-1)^j (Pi_j/j) t^-1 zeta~^j H_j for j >= 2, and -alpha is the j = 1 term of the same expression. Hence, with Q_t(w) := (1/t) Sum_{j>=1} (-1)^j (Pi_j/j) (zeta~ w)^j H_j(T) and w(lambda) := 1/u(lambda) = lambda + tau + O(1/lambda), we get Pcal_t = Q_t(w) + alpha w and phi_t = -[Q_t(w(lambda))]_+ + [Q_t(w(lambda))]_+(0).
(b) Let g be the compositional inverse of f(mu) := 1/(C u_t(mu)), a t-series of Laurent series. Lagrange inversion (with the formal change of variables) gives [w^n] g = (C^n/n) Res_mu u_t(mu)^n dmu for n >= 1.
(c) Residue identity: Res_mu u_t^n dmu = (-1)^(n-1) t^(n-1) e^(-n tau^2 t) Pi_n H_n(T). Both sides are t-series with coefficients polynomial in tau, so it suffices to take real tau in (0,1) and small t, where u_t^n = u^n e^(-nK_t) holds analytically for |t mu| < 2pi and the formal residue is the contour integral over |mu| = r. Substitute Lambda = e^(t mu) and then Lambda = 1 + (e^t-1)s. Expanding (1+(e^t-1)s)^(-n tau-1) binomially gives (-1)^(n-1) t^(n-1) e^(-n sigma t) Pi_n 2F1(n+n tau, -n tau; n; 1-e^t). Pfaff, 2F1(a,b;c;z) = (1-z)^-b 2F1(c-a,b;c;z/(z-1)), turns this into e^(n tau t) H_n(T). Consequently [w^n](g + Q_t) = 0 for every n >= 1 exactly when C^n t^n e^(-n tau^2 t) = zeta~^n, which is the stated C. So g(w) = -Q_t(w) + c0 + O(1/w), and lambda~ := g(w(lambda)) = phi_t(lambda) + eps(lambda) with eps = O(lambda^0).
(d) lambda~ solves C u_t(lambda~) = u(lambda). The quantity Z(lambda) := log C - K_t(phi_t) + tau log((phi_t-1)/(lambda-1)) - (1+tau) log(phi_t/lambda) is a t-series of polynomials, since phi_t fixes 0 and 1, and it equals log u_t(phi_t) - log u_t(lambda~). If eps != 0 had lowest order t^m, then [t^m] Z = eps_m (lambda-1-tau)/(lambda(lambda-1)) = O(1/lambda) would be a polynomial, hence 0, hence eps_m = 0, a contradiction. So eps = 0, which is FC. QED
Remarks. C also equals u(1+tau)/u_t(mu_c), where mu_c is the critical point of u_t, (1+tau)Lambda^2 - (1+2tau)e^t Lambda + tau e^t = 0 with Lambda = e^(t mu); this involves sqrt(1+4 sigma T) (t_crit.pl, 4 framings, order 6). With the constant alpha instead of C the identity fails at t^2 (refuted, t_flow_curve.pl). Equivalent forms of T1 used below: with F(y) := Sum_c c y_{c-1} u(y)^c = h d_lambda Psi and relation_k = kappa_k R_k,
 Sum_k (k+1) R_k s^(k+1) = (1+tau)(1-(1+tau)s)F(s)/(1-s) - tau(1-tau s)F(s/(1+s))/(1+s),
which follows from Res Psi s/(1-s lambda)^2 = -(1/s) Psi'(1/s) and the analogous identity at lambda-1.
Checks of T1 and its ingredients:

| check | range | result |
|---|---|---|
| t_beta_closed.pl: exact polynomials beta_{k,L}(tau) from full KP (interp_beta.pl) vs closed form | k <= 20 | 2876 evaluations, 0 mismatches |
| t_T1_modp.pl: full KP mod p | k <= 23, tau = 2/7, -1/5, -3/7, 5/3 | 1104 coefficients, 0 mismatches |
| ode_leaf.pl | order 100, tau = 2/7, -1/5, -3/7, 3/11 | confirmed |
| t_zexp.pl: independent z-expansion implementation, no y-Hankel form | tau = 2/7 (order 4), -3/5 and 3 (order 6) | confirmed |
| t_phi_explicit.pl: FC with log C = log C_a | 5 framings, order 6 | holds |
| t_res_ut.pl: residue identity | n <= 5, 3 framings | all agree |
| t_generator.pl part 1: branch-point identity on the leaf curve | order 6 | holds |

## 8. Special framings

LEMMA 8.1. Let tau0 = -k0/q with gcd(k0,q) = 1 and 1 <= k0 < q, and put F(y) := Sum_c c y_{c-1}(tau0) u(y)^c with u = y(1-y)^tau0. The two-leaf conditions Y(m,n)(tau0) = 0 for all m >= 1, q|n, m+n <= 2G-1 are equivalent to [y^N] p = 0 for q <= N <= 2G-2, where p(y) := (1-y)^(k0-1)(1-(1+tau0)y)F(y) and Y(m,n) = [u1^m u2^n] y1 y2 (F(y1)-F(y2))/(y1-y2) (Theorem D, R empty). These conditions hold at every tau0 because Xi_G(m,n) is finite at sigma0 != 0 and Pi_n(tau0) = 0 when q|n. They include the one-leaf conditions (m = 1).
Proof. Lagrange inversion (y = u(1-y)^-tau0) gives [u^(qm)] g = (1/(qm)) [y^(qm-1)] g'(y)(1-y)^(m k0). The y1-coefficients of A = y1 y2 D_F are y2 S_i(y2), with S_i(y) = Sum_j f_{i+j+1} y^j, so the condition (i,m) reads Sum_r (-1)^r binom(m k0, r)(qm-r) f_{i+qm-r} = 0. Since Sum_r (-1)^r binom(mk0,r)(qm-r) y^r = mq (1-y)^((m-1)k0) (1-y)^(k0-1)(1-(1+tau0)y), the condition (i,m) is [y^(i+qm)]((1-y)^((m-1)k0) p) = 0 for i+qm <= 2G-2. The conditions with m = 1 say that p mod y^(2G-1) has degree <= q-1. Those with m >= 2 then follow, because (m-1)k0 + q - 1 < qm. QED
(t_special.pl: equality of row spaces at every framing with q <= 2G-2 and k0 < q, for G = 3..9, 0 differences; with a wrong exponent 20 of 21 framings differ.)

## 9. Forcing theorem

THEOREM 9.1. Let G >= 3, D odd with 3 <= D <= 2G-3, Lemma A in genera < G, and y_L = 0 for L < D. Let tau0 = -k0/q with gcd(k0,q) = 1, 1 <= k0 <= D+1, 2k0 < q, q >= D+2 and 3(q-k0) < 2G. Then Xi_G(D+1)(tau0) = 0, so P_{D+1} := sigma^G Xi_G(D+1) vanishes at sigma0 = -k0(q-k0)/q^2.
Proof. By Lemma 8.1, F = p/((1-y)^(k0-1)(1-(1+tau0)y)) mod y^(2G-1), with p a polynomial of degree <= q-1 (q < 4G/3 <= 2G-2), p = O(y^(D+1)) and [y^(D+1)]p = (D+1) y_D(tau0).
By (I1) and T1, R_k = 0 for even k <= 2G-4 (kappa_k != 0). In the s-form, at tau0, with A := q-k0 and B := k0:
 q R(s) = A(1-s)^-B p(s) + B(1+s)^(B-1) p(s/(1+s))  mod s^(2G-1), and its odd coefficients up to s^(2G-3) vanish.
With s = 1/lambda, pi(lambda) := lambda^(q-1) p(1/lambda) (deg pi <= q-D-2) and r := pi/(lambda^A (lambda-1)^B), one has q R(1/lambda) = lambda X(lambda), where X := A r(lambda) + B r(lambda+1) = O(lambda^(-D-2)).
So E := X(lambda) + X(-lambda) = O(lambda^-2G). E is a rational function whose denominator divides (lambda(lambda^2-1))^A, so its numerator has degree <= 3A - 2G < 0. Hence E = 0 and X is odd.
X has a pole of order <= B at lambda = 1 and of order A - ord_0 pi at lambda = -1. Oddness forces ord_0 pi >= A - B, so pi = lambda^(A-B) pi_1 with deg pi_1 <= 2B-D-2.
Then X = Nn/(lambda(lambda^2-1))^B with Nn = A pi_1(lambda)(lambda+1)^B + B pi_1(lambda+1)(lambda-1)^B. Oddness is equivalent to (lambda+1)^B s1(lambda) + (lambda-1)^B s1(-lambda) = 0, where s1(lambda) := A pi_1(lambda) + B pi_1(1-lambda). This makes (lambda-1)^B divide s1. Since deg s1 <= 2B-D-2 < B (because B <= D+1), s1 = 0.
Then pi_1(1-lambda) = -(A/B) pi_1(lambda), applied twice, gives pi_1 = (A/B)^2 pi_1, so pi_1 = 0 (A != B). Hence p = 0 and y_D(tau0) = 0. Finally Pi_{D+1}(tau0) != 0 because q does not divide D+1. QED
(t_forcing.pl: all 1131 tuples in the hypothesis range for G <= 20 are forced, 0 failures; outside the range for G <= 9, 126 of 320 are still forced. eq_kernel.pl handles single cases.)

## 10. Zero count

Z1(D) = floor(D/2) + Sum_{a=3}^{D+1} (phi(a)/2) max(0, floor(D/(a-1)) - [a|D+1]) counts the zeros, with multiplicity, of P_{D+1} at the framings -c/a with a <= D+1 (Route A, Corollary E; re-derived: T = (a^m, r+1) with m = floor(D/(a-1)) has m parts divisible by a; at a = 2 the sigma-order is half the tau-order).
Z2(G,D) = #{reduced k0/q : k0 <= D+1, 2k0 < q, q >= D+2, 3(q-k0) < 2G} counts the simple zeros of Theorem 9.1. All these sigma0 are nonzero and distinct, and the two sets have disjoint denominators.
PROPOSITION 10.1. Z1(D) + Z2(G,D) >= G-1 for all G >= 8 and all odd 3 <= D <= 2G-3.
Proof.
- Exact integer count for 3 <= G <= 1000 (count_zeros2.pl: Moebius counting, cross-checked with the direct enumeration of count_zeros.pl). The only failures are G = 4..7 at D = 3; the minimal margin is at D = 3, e.g. 618 at G = 800 and 771 at G = 1000.
- For G > 300 two bounds are used:
  (B1) Using k0 = 1..4 and g = ceil(2G/3) - 1: Z2 >= (g-D) + (g-D)/2 + (2(g-D)-2)/3 + (g-D-2)/2 >= (16/9)G - (8/3)(D+1) - 5/3.
  (B2) From 2 Phi(n) - 1 = Sum_d mu(d) floor(n/d)^2 >= (2 - pi^2/6) n^2 and Z1 >= floor(D/2) + (1/2) Sum_{j=2..4} (Phi(floor(D/j)+1) - 2): Z1(D) >= 0.0376 D^2 + 0.5D - 2.75 for D >= 8.
  If D+1 <= (7G-6)/24, then (B1) gives Z2 >= G-1. Otherwise D >= 87 and G-1 < (24D+23)/7 <= Z1(D) by (B2), which holds for D >= 80.
- (B1) and (B2) were also checked numerically: 0 violations for G <= 1000 and D <= 2000. QED

## 11. Conclusion: (OL) for G >= 8 and Lemma A in all genera

Under the hypotheses of (OL), P_{D+1} = sigma Q(sigma) with deg Q <= G-2 by (a),(b). Q vanishes at the Z1(D) + Z2(G,D) >= G-1 nonzero points (with multiplicity), so Q = 0. Hence (OL) holds for every G >= 8.
Induction on h: Lemma A holds for h <= 4 (exact data) and for 5 <= h <= 7 by Theorem F. For G >= 8, inside genus G proceed by induction on the level D:
- D = 1 is Theorem A;
- even D: T1 at k = D gives relation_D = kappa_D y_D, equivalently the parity README (f);
- odd D >= 3: (OL);
- Corollary E then gives every T at level D.
This proves Lemma A for all h, and Corollary R follows by brief s.1.
The proof depends on: Route A Theorem D (MV + Fact 3.5 + Takasaki-Takebe), Corollary E and Theorem F; Route D's single-vertex polar relation ((H) = Brini-Schuler Prop. 3.1, Lemma 5.3, (V)); leaf_structure_proof.md Theorem 11.1; brief 2(a),(b); Theorem A; the exact data for h <= 4.

## 12. R1a computations (computer-verified, exact modular arithmetic; nullity 0 mod p implies nullity 0 over Q)

Global nullities of the whole system (one-leaf unknowns c = 3..2G-2) in the correct window W = G-1 (span sigma..sigma^(G-1)):

| G | TLR | TLR+CF | full | full+CF |
|---|---|---|---|---|
| 3 | 2 | 1 | 1 | 0 |
| 4 | 3 | 2 | 0 | 0 |
| 5..10 | 0 | 0 | 0 | 0 |

With W = G: TLR has nullity 4,7,6,1,0 at G = 3..7 (Route A's 11,11,7 included c = 2). With the closed-form relations (CF) added the nullities are 2,4,0,0 at G = 3..6, and full+CF has 1,0,0 at G = 3..5.
Single-level problems at odd D: nullity 0 for every G >= 5 (G <= 10), even with TLR alone. The exceptions are G = 3, D = 3 (TLR: 1) and G = 4, D = 3 (TLR: 2); full regularity gives 0 there too.
Pointwise scan (scan_points.pl: two-leaf + closed-form rows from the full KP, all special framings with a <= 2G-2) closes (OL) for every (G,D) with 5 <= G <= 12 except (5,5), which Route A's Z(5) = 6 closes (out/scan_G12_two_cf.txt, 43 min).
Smallest rigid subsystem: two-leaf regularity at the framings -k0/q plus the closed-form (polar) relations. With T1 this is exactly what s.9 uses: one relation per even order, and pointwise oddness of X(lambda).

## 13. Refuted or abandoned claims

- (FC1) u_t(phi_t) = u/alpha: refuted (fails at t^2; the correct constant is C = alpha e^(tau^2 t)/H_1).
- Deformed-curve ansatz Y_P(z) - Sum P_j z^j = lambda + tau: refuted (the two-leaf test fails at J = 5; t_ansatz2.pl).
- A guessed closed form zeta = (e^t-1)e^(tau t)((1+tau)e^t - tau)^tau: refuted (t_zeta.pl); zeta involves sqrt(1+4 sigma T).

## 14. Points for the referee

1. The use of Route A Theorem D at every point of the slice, in generating form along a curve, together with the level bookkeeping of Lemma 2.1.
2. The string normalisation Y(T+{1}) = |T| Y(T) (consistent with Route A and kp_point.pl).
3. The identification relation_k = [t^k] F_G(P(t)) via b_j = A_j alpha^-j/c^2.
4. Step (d) of s.7 (polynomiality forces eps = 0).
5. The Lemma 8.1 bookkeeping of the range m+n <= 2G-1.
6. The degree count in s.9.
7. Route A's Z(D) multiplicities (re-derived in s.10).

## 15. Files (all in the work folder)

Library and engines:
- r1lib.pm: BigRat helpers.
- kpmod.pm: modular Theorem-D propagation.
- kp_point.pl, kp_single.pl, kp_dump.pl: exact single-level KP problem and closed-form rows.
- t_kpmod.pl: modular engine vs exact.

Leaf curve:
- t_leafcurve.pl: Theorem 11.1, sum rule, b_j.

T1 evidence:
- interp_beta.pl (out/beta_poly_K14.txt, out/beta_poly_K20.txt) and t_beta_closed.pl.
- t_T1_modp.pl (out/T1_modp_K23.txt).
- ode_leaf.pl (out/ode_leaf_K100.txt) and ode_state.pl.
- t_zexp.pl (out/t_zexp.txt).
- t_Rformula.pl, t_depth2.pl.
- beta_sigma.pl, beta_kappa.pl, R_bases.pl, pk_fit.pl, beta_table.pl, beta_to_F.pl: discovery.

Flow theorem and its ingredients:
- fit_multiplier.pl, fit_flow.pl (out/flow_tau*.txt): discovery.
- t_flow_K.pl, t_flow_curve.pl, t_generator.pl.
- ode_path.pl (out/ode_path.txt) and ode_path2.pl (out/ode_path2.txt).

Curve identity:
- t_phi_explicit.pl (out/t_phi_explicit.txt), t_res_ut.pl (out/t_res_ut.txt), t_crit.pl (out/t_crit.txt), t_zeta.pl.

Special framings and forcing:
- t_special.pl (out/t_special_G8_9.txt).
- eq_kernel.pl, t_forcing.pl (out/t_forcing_G20.txt).

Zero count:
- count_zeros.pl, count_zeros2.pl (out/count_zeros2_G3to1000_raw.txt).

Global nullities and scans:
- global_null.pl, global_suite.pl (out/global_suite_G3to8.txt, out/global_suite_G5.txt, out/global_suite_G9.txt, out/global_suite_G10.txt).
- scan_points.pl (out/scan_G7_two_cf.txt, out/scan_G12_two_cf.txt, out/scan_G9_*.txt).

Exploratory (not used in the proof):
- t_ansatz.pl, t_ansatz2.pl, ycoord3.pl, t_curve_deform*.pl, t_genus0_curve.pl, eq_scan.pl, t_time.pl.

Mathematica:
- R1_checks.wls: checks (C1)-(C6) for the main session.

Copies:
- routeA_copy/: read-only copies of Route A scripts.
