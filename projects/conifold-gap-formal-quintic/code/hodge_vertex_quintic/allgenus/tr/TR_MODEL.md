# Topological recursion for the LG model x^4 + v/x: dictionary, frozen theory, all-genus gap, v -> 0 (allgenus/tr, 2026-11-13)

Labels: PROVED (argument here, or quoted theorem with arXiv number, statement, checked hypotheses), VERIFIED (exact computation; field
or primes and orders stated), CONJECTURAL, SKETCH, REFUTED. AI-derived and AI-checked only.
TR is used only for the MODEL curve (user ruling 2026-11-13), never for the quintic's mirror curve. All code is new (this folder);
the unreviewed scripts of the first attempt (`tr_core.py`, `tr_rec.py`, `tr_xi.py`) were replaced. Project modules imported read-only:
`fq_fp`, `fq3_mvertex` (MVertex6, MConst); data read: `natural_form/lg/lg_symbol_rows_Q.json`, `lg_symbol_89.json`, `lg_symbol_61.json`.

## 0. Results in one page

1. **Literature (s.1).** arXiv:1509.06954 (DNOPS) Theorem 8.1 is exactly our case: for f(p,a) = p^n + a_2 p^(n-2) + .. + a_n +
   a_(n+1)/(p - a_1) the Givental-Teleman CohFT of the Hurwitz Frobenius manifold (at any semisimple point, DNOPS sign of the metric)
   equals the CohFT reconstructed from TR on (P^1, x = f, y = p, B = dp1 dp2/(p1-p2)^2) via their (2-30)-(2-35). Hypotheses hold on the
   line a = (0,0,0,0,v), v != 0 (five simple critical points p^5 = v/4, distinct critical values 5 p_c^4). QUOTED, hypotheses checked.
2. **Dictionary, derived and VERIFIED (s.2, s.3).** Use x = kappa W, W = p^4/4 + u/p (u = 1: c_j = zeta^-j, colour j of the project).
   The Laplace transform (2-31) of B gives R_TR(z) = M_0(z/kappa) EXACTLY through z^8 (all 225 entries, Q(zeta), kappa = 1, 4, -1);
   (2-32) (Laplace transform of dy, i.e. omega_{0,1}) gives the frozen unit M_0^-1 1^ and hence the dilaton leaf, through z^8.
   With Lambda = -5 kappa:  omega_{g,n} = Lambda^(g-1) Sum kappa^-(3g-3+n-Sum d) <prod tau_{d_l}(v_l)>^frozen_g prod d xi_{v_l,d_l}.
   For kappa = 1 (x = W) the TR CohFT IS the project's frozen theory M_0.omega^ with metric times -1/5: the normalisation factor
   is (-5)^(1-g), not 5^(1-g) (5^(1-g) belongs to x = -W, which flips psi -> -psi: observed first, s.2.4).
   **Check against the project's frozen theory:** TR x (-5)^(1-g) == MVertex6(M_0) on all 72120 colour-basis correlators with
   psi-deficit <= 6 of (0,3), (0,4), (0,5), (1,1), (1,2), (1,3), (2,1), (2,2), (3,1) (4665 non-zero) at 2^89-1 AND 2^61-1, plus F_2, F_3.
   Exact flat values: <H>_1 = -1/12, <tau_1(1)>_1 = 5/24, <H,H>_1 = 1/12, F_2 = 1/1200, <H>_2 = -1/600 (H = N_1 Sum zeta^2i et_i;
   the project's quoted list 1/12, 1/12, 1/6; 1/600, 1/200 uses phi = -H and is reproduced). Full tables: tr_run_exact.log.
3. **All-genus model gap (s.4): F_g = B_2g/(2g(2g-2)) v^(2-2g) for x = p^4 + v/p (DNOPS conventions), i.e. the frozen c_g =
   |B_2g|/(2g(2g-2) 5^(g-1)) for every g >= 2.** PROVED modulo one quoted theorem (Hock, arXiv:2502.18115, Cor. 1.2: the x-y swap
   formula for F_g when y is unramified); the local evaluation (residue at the simple pole of x = Bernoulli by Lagrange inversion,
   residue at infinity = 0 by a degree count) is proved here. VERIFIED independently by TR: g = 2, 3 exact over Q(zeta), g = 2..5 at
   both primes; the formula itself evaluated exactly for g <= 12 and for several polynomial parts.
   Stronger (same proof): on the WHOLE Hurwitz space H_{0;1,r} (any r >= 1, any point) F_g depends only on the residue v of x at its
   simple pole. Consequence: every genus-g >= 2 primary correlator of the model with an insertion other than H vanishes (the
   "primary-sector lemma" of UNIFORM_REDUCTION s.6.5 item 3, for g >= 2). VERIFIED by TR at random points of H_{0;1,4} (2 points x 2
   primes, g = 2, 3) and for r = 1..6 at the special point.
4. **v -> 0 (s.5).** At fixed p != 0 every computed omega_{g,n} (n >= 1; 13 cases up to (4,1), (3,2), (2,3), (1,4), (0,5), both
   primes: VERIFIED; all (g,n): CONJECTURAL) is a convergent power series in u with a non-zero limit: no omega_{g,n} is singular, there
   is no leading singular part. The singular F_g (n = 0) come from the residues at the branch points, which collide at p = 0
   (p ~ v^(1/5)). The limits obey the selection rule of the 4-spin theory of x = p^4/4 (genus-1 value <tau_1(e_0)>^(4-spin)_1 = 1/8).

## 1. The theorem from the literature (task item 1)

**DNOPS (arXiv:1509.06954v2, J. Inst. Math. Jussieu 18 (2019)); read: s.1-8 from the PDF.**
* TR conventions (2-26)-(2-29): omega_{0,1} = y dx, omega_{0,2} = B,
  omega_{g,k+1}(z_0, J) = Sum_i Res_{z->c_i} [Int_z^{sigma_i(z)} omega_{0,2}(., z_0)] / [2(omega_{0,1}(sigma_i z) - omega_{0,1}(z))]
  x [omega_{g-1,k+2}(z, sigma_i z, J) + Sum' omega omega] (Sum' includes omega_{0,2}-factors). This is the Eynard-Orantin kernel (both
  numerator and denominator flip relative to (1/2) Int_{sigma z}^{z} B/((y(z) - y(sigma z)) dx)). Same kernel in Hock (2.1).
* Dictionary (2-30)-(2-35): local coordinate x = x(c_i) - w_i^2/2; Delta_i^(-1/2) = dy/dw_i(0) (2-30); R^-1(zeta^-1)^j_i =
  -(2 pi zeta)^(-1/2) Int B(w_i, w_j)/dw_i|_{w_i=0} e^((x(w_j) - x(c_j)) zeta) (2-31); Sum_k R^-1(zeta^-1)^i_k Delta_k^(-1/2) =
  (zeta/2pi)^(1/2) Int dy e^((x - x(c_i)) zeta) (2-32); xi^i(z) = Int^z B(w_i, .)/dw_i|_0 (2-34);
  omega_{g,k} = Sum Int alpha_{g,k}(e_i1, .., e_ik) prod psi_j^(d_j) d((d/dx)^(d_j) xi^(i_j)) (2-35), e_i the normalised canonical basis.
* Theorem 5.1: if Dubrovin's construction gives a compact genus-0 curve with global coordinate p, lambda(p = infinity) = infinity,
  d_lambda p without zeros off the cuts, and exactly one critical point in each singular fibre, then the CohFT and TR on
  (x = lambda, y = p, B = dp1 dp2/(p1-p2)^2) are related by (2-35). **Theorem 8.1** (Hurwitz case, superpotential (8-1)
  f = p^n + a_2 p^(n-2) + .. + a_n + a_(n+1)/(p - a_1)): the CohFT of the Frobenius manifold with superpotential f coincides with
  the one reconstructed from (CP^1, f, p, B); for n > 1 the proof uses p_Dub = p - a_1/(n+1) (a constant shift of y). Metric
  (2-7): eta = -Sum_crit Res d'(lambda dp) d''(lambda dp)/d_p lambda, i.e. MINUS the metric (1.6) of arXiv:2211.12259 and of
  LG_MODEL.md s.5.1 (= +Sum_crit Res); the two CohFTs differ by Omega_g -> (-1)^(g-1) Omega_g.
* **Application to lambda = p^4 + v/p (n = 4, a = (0,0,0,0,v)), v != 0:** critical points p^5 = v/4, non-degenerate (lambda'' = 20 p_c^2),
  critical values 5 p_c^4 pairwise distinct, so the point is semisimple and each singular fibre has one critical point: Theorem 8.1
  applies at every point of the line with v != 0 (and on an open dense set of the whole Hurwitz space). Conformal: d = (n-2)/n = 1/2,
  unit d/da_4, Euler sum (i/4) a_i d/da_i, all as in LG_MODEL.md s.5.
* arXiv:1211.4021 (DOSS, CMP 328 (2014)): the local version (TR on a union of discs <-> Givental formula), the source of (2-30)-(2-35).
* arXiv:2211.12259 (Carlet-van de Leur-Posthuma-Shadrin, CNTP 17 (2023)); read s.1-2.1: same manifold H_{0;1,N-1}, special point
  f = p^(N-1) + 1/p, flat coordinates (1.7), metric (1.8), d = (N-3)/(N-1), R from (1.18); main theorem: the total descendant
  potential is a tau function of the extended rationally (N-1)-constrained KP hierarchy. Not used for any result here.
* arXiv:2502.18115v2 (Hock, 2025; read s.1-3.2 and Prop. 3.4-3.5): Theorem 1.1/3.2 expresses (2-2g)(F_g - F_g^vee) (F_g^vee the free
  energy of the x-y swapped curve) as a sum of local residues at the points where y dx or x dy has a non-zero residue; Corollary 1.2:
  if y = z is unramified, F_g^vee = 0 and F_g is given by the explicit residue formula (1.3) (quoted in s.4). Hypotheses: compact
  Sigma, x, y meromorphic (log singularities allowed) with simple ramification points, B the Bergman kernel; TR kernel and F_g
  definition (1.1) identical to ours. Used in s.4.

## 2. The dictionary (task item 1, derived)

### 2.1 Set-up
x = kappa W(p), W = p^4/4 + u/p, y = p, B = dp1 dp2/(p1-p2)^2; u = 1 (s = 1) in all computations (general u by homogeneity, s.5).
kappa = 1: x = W (the project's W-form, LG_MODEL.md s.3); kappa = 4: x = lambda = p^4 + v/p with v = 4u (the task's curve);
kappa = -1: x = -W. Branch points c_j = zeta^-j (colour j: the project's colour i sits at the critical point zeta^-i s),
x(c_j) = (5/4) kappa zeta^j, x''(c_j) = 5 kappa c_j^2. Local coordinates: x - x(c_j) = (5/2) kappa wt_j^2, wt_j = c_j (p - c_j) + ..;
DNOPS w_j^2 = Lambda wt_j^2 with **Lambda = -5 kappa**.

### 2.2 TFT and xi-functions (PROVED, elementary)
(2-30): Delta_j^TR = Lambda c_j^2; the frozen theory has Deltahat_j^(1/2) = zeta^4j = c_j, so e_j = sqrt(Lambda) et_j and the TR TFT on
the vectors et_j is Lambda^(g-1) times the frozen TFT. (2-34): xi^j = xit^j/sqrt(Lambda), xit^j(p) = -1/(c_j (p - c_j)).
Flat basis H^a = N_a Sum_i zeta^(i(a+1)) et_i (N = (1,-1,-1,1,1), eta(H^a,H^b) = 5 delta_{a+b=3}); the function multiplying a
<.. tau_d(H^a) ..>-insertion is xif_{a,d} = (d/dx)^d xif_a with
   xif_0 = -p^4/(p^5-1),  xif_1 = 1/(p^5-1),  xif_2 = p/(p^5-1),  xif_3 = -p^2/(p^5-1),  xif_4 = -p^3/(p^5-1)
(xif_a is the xi-function of the dual vector H^(3-a)/5; xi of a flat vector d/dt_alpha is -d_alpha x/x' = d_alpha y|_x: 1 <-> d_4,
H <-> d_5 (residue), H^2 <-> d_1, H^3 <-> d_2, H^4 <-> d_3, as in LG_MODEL.md s.5.2). VERIFIED (tr_run_exact.log): string pairing
Sum Res y d xif_{a,d} = delta_{(a,d),(0,0)} and dilaton pairing Sum Res Phi d xif_{a,d} = -delta_{(a,d),(0,1)} (Phi = Int y dx), a <= 4,
d <= 3, exact; i.e. the unit and the dilaton vector are where they must be with these normalisations.

### 2.3 The R-matrix (task item 3) -- VERIFIED exactly, and why
`tr_run.py rmat` (log tr_run_rmat.log), exact over Q(zeta): with X(z)[j][i] := (2-31) (row j = branch point where the Laplace transform
is taken, column i = the xi-function),
   **X(z) = M_0(z/kappa)^-1 = M_0(-z/kappa)^T  for z^0..z^8, all 25 entries, kappa = 1, 4, -1**,
M_0 = R~ from lg_symbol_rows_Q.json (C_n[j][i] = (1/5) Sum_a zeta^(j(4-a)+i(a+1-n)) g^a_n), itself equal to the project's symbol dump
(lg_symbol_89/61.json, 175 entries at each prime). The transposed reading X = (M_0^-1)^T fails from z^2 (z^0, z^1 are symmetric).
So R_TR(z) = M_0(z/kappa): for x = W exactly M_0; for the task's lambda = p^4 + 4/p exactly M_0(z/4) = M_0(zhat) with z = 4 zhat.
Why the sign of kappa matters (ARGUED): integrating by parts, the Laplace transform of d xi_{(a)} is z^-1 Int e^(x/z) d_a x dp, the
oscillating integral of e^(x/z); its stationary phase in the normalised canonical frame is R_osc(z)^T = R_osc(-z)^-1, so (2-31) returns
R_DNOPS(z) = R_osc(-z), the R-matrix of Int e^(-x/z). The project's M_0 is the stationary phase of Int e^(-W/z) (LG_MODEL.md s.0,
lambda/z_LG = -W/z), hence kappa = +1. All orders: both R_TR and M_0 are homogeneous R-matrices of the same (U, V)
([R_(m+1), U] = (m + V) R_m determines R from R_0 = 1 for semisimple conformal Frobenius manifolds); U agrees (critical values) and
V = [R_1, U] agrees because R_1 agrees -- so the z^8 check plus uniqueness gives all orders (PROVED modulo that standard uniqueness).
**Dilaton leaf (omega_{0,1}):** (2-32) holds exactly through z^8 for all three kappa, and its right-hand side equals M_0(z/kappa)^-1 1^,
1^ = Sum_k et_k/c_k (the frozen unit); so omega_{0,1} = y dx produces exactly the frozen dilaton leaf T(z) = z (1^ - M_0^-1(z) 1^)
(kappa = 1); colour-0 coefficients z^2..z^5: -3/20, -9/800, 11571/80000, 2415609/6400000.

### 2.4 The normalisation and a negative result
From 2.2-2.3: omega_{g,n} = Lambda^(g-1) Sum kappa^-(3g-3+n-Sum d) <prod tau_{d_l}(v_l)>^frozen_g prod d xi_{v_l,d_l} (xi descendants with
d/dx, x = kappa W). For kappa = 1: factor (-5)^(g-1), no psi-rescaling. **The factor 5^(g-1) suggested in the task is wrong in sign
for x = W**; it is right for x = -W, which however gives the frozen theory with psi -> -psi: the first exploratory run (x = -W,
factor 5^(g-1)), compared by hand with MVertex6 on 12 correlators, agreed exactly on those of even class degree D = 3g-3+n-Sum d and
had the opposite sign on all with odd D (<H>_1, F_2, <tau_1(1)>_2, <tau_3(H^3)>_2, <t0(H^3) t1(H^3)>_1). This is the z -> -z of 2.3.
Task normalisation x = lambda = p^4 + v/p, v = 4: omega = (-20)^(g-1) Sum 4^-D <..>^frozen prod d xi (d/d lambda); F_g = (-5/16)^(g-1) c_g.
General u (s^5 = u): p -> s p maps W_1 to s^-4 W_u, so c_j = zeta^-j s, R_TR(z) = M_0(z/s^4) = M_0(zhat) and Delta_j = -5 zeta^-2j s^2:
the hatted frozen theory (zhat = z/s^4, Hhat^a = s^(1-c_a) H^a, colour basis) is TR on W_u read in these variables (homogeneity, s.5).

## 3. TR computations and comparison with the frozen theory (task item 2)

Engine `tr_rec.py`: recursion (2-28) by exact residue calculus in the basis e_{j,k} = dp/(p - c_j)^(k+1); residues at c_0 only and the
Z/5 rotation C[((r,m),(j+r,k)..)] = zeta^(-r(m + Sum k)) C[..] (TR is equivariant: y dx, B invariant under p -> zeta^-1 p); full mode
(all branch points, no symmetry) for generic curves. Fields: exact Q(zeta) (`tr_field.Qz`) and F_p (2^89-1, 2^61-1).
* **Regression (tr_run_airy.log):** Airy curve x = -p^2/2, y = p: omega = Sum <tau..>_g prod (2d+1)!! dp/p^(2d+2) (Witten-Kontsevich)
  for (0,3), (1,1), (0,4), (1,2), (2,1), (0,5), (1,3), (2,2), (3,1), and the dilaton equation; Q(zeta) and both primes: PASS.
* **Exact run on x = W (tr_run_exact.log, 37 s):** omega_{0,3}, (1,1), (0,4), (1,2), (2,1), (1,3), (2,2), (3,1): symmetric in all arguments;
  in the span of the d xi with zero remainder (colour and flat bases); all flat correlators rational; TR dilaton equation
  Sum Res Phi omega_{g,n}(., J) = (2-2g-|J|) omega_{g,n-1}(J); string and dilaton equations of the flat correlators: all PASS.
  F_2 = -1/240, F_3 = 1/1008 (x = W), i.e. frozen F_2 = 1/1200, F_3 = 1/25200.
* **Against the project's code (tr_run_frozen.log):** see s.0 item 2: 72120 correlators (4665 non-zero) per prime, 0 mismatches.
* **Selected exact flat correlators of the frozen theory** (H = H^1; `t_d(H^a)` = tau_d(H^a); complete lists in tr_run_exact.log):
  - genus 0: <1,1,H^3> = <1,H,H^2> = <1,H^4,H^4> = <H^2,H^2,H^4> = 5, <H,H,H> = <H,H^3,H^4> = <H^2,H^3,H^3> = -5  (= 5 N_a N_b N_c delta_{a+b+c=3});
    <H,H,H,H>_0 = 5, <H^2,H^2,H^2,H^3>_0 = 10.
  - genus 1: <t1(1)> = 5/24, <H> = -1/12; <H,H> = 1/12, <t0(H) t2(H^4)> = <t1(H) t1(H^4)> = <t2(H) t0(H^4)> = -5/24,
    <t0(H^2) t2(H^3)> = <t1(H^2) t1(H^3)> = <t2(H^2) t0(H^3)> = -5/24, <t0(H^3) t1(H^3)> = 1/24; <H,H,H>_1 = -1/6.
  - genus 2 one-point: <t1(1)> = 1/600, <H> = -1/600, <t2(H^4)> = 1/1200, <t3(H^3)> = 7/5760, <t4(H^2)> = -5/1152 (nothing else).
    genus 2, <H,H> = 1/200 and 44 two-point values in the log.
  - genus 3 one-point: <t1(1)> = 1/6300, <H> = -1/6300, <t2(H^4)> = 1/12600, <t3(H^3)> = 1/37800, <t4(H^2)> = -1/151200,
    <t5(H)> = 53/1209600, <t6(1)> = -1/46080, <t7(H^4)> = 5/82944.
  Remarks (observations, data for Stage 1 of the plan, no claim): <H>_g + <tau_1(1)>_g = 0 for g = 2, 3 (the Euler-dilaton vector
  e = z1 + uH of UNIFORM_REDUCTION s.6.4 annihilates the one-point functions; genus 1: 5/24 - 1/12 = 1/8, the anomaly).
  Exponents: <tau_d(H^a)>_g(u) ~ u^e, e = 2 - 2g + (4d - beta)/5, beta = 4, 5, 1, 2, 3 for a = 0..4 (homogeneity); singular means e < 0.
  The linear "one-point law" |<xi>_g| = |Pi_C(xi) d_v F_g| with the plan's Y^C = (-v/z, -1, -v^4/(24z^4), -v^3/(6z^3), v^2/(2z^2)) + ..
  holds in absolute value for every SINGULAR genus-2 one-point function (t1(1), H, t2(H^4)) and for the genus-3 ones t1(1), H, t2(H^4),
  t3(H^3), t4(H^2) (e.g. |<t3(H^3)>_3| = (1/6) 4 c_3 = 1/37800, |<t4(H^2)>_3| = (1/24) 4 c_3 = 1/151200), but NOT for the singular
  <t5(H)>_3 = 53/1209600 (the next term +v^4/(96 z^5) of Y^C_1 predicts 1/604800); the regular genus-2 values t3(H^3), t4(H^2) are not
  of that form either (they need not be). Signs were not matched (H, z conventions).

## 4. The free energies for all genera (the model gap)

### 4.1 Statement and proof
**Theorem G.** Let x = P(z) + v/(z - a), P a polynomial of degree r >= 1, v != 0, with simple ramification; y = z; B = dz1 dz2/(z1-z2)^2.
Then for every g >= 2:   F_g = B_2g/(2g(2g-2)) v^(2-2g)   (EO/DNOPS conventions, F_g = (1/(2-2g)) Sum Res Phi omega_{g,1}).
In particular for lambda = p^4 + v/p: F_g = B_2g/(2g(2g-2)) v^(2-2g); for x = W: B_2g/(2g(2g-2)) u^(2-2g); and the frozen theory has
F_g = (-5)^(1-g) B_2g/(2g(2g-2)) = |B_2g|/(2g(2g-2) 5^(g-1)) = c_g for every g (Conjecture N's vertex constant).
Status: **PROVED modulo the quoted Corollary 1.2 of arXiv:2502.18115** (and, for the frozen statement, modulo DNOPS Thm 8.1 + s.2.3).

*Proof.* (i) Quoted (Hock, Cor. 1.2, eq. (1.3)): for y = z and x without logarithmic singularities,
F_g = (1/(2-2g)) Sum_b Res_{z=b} [hbar^(2g-1)] dz Sum_{m>=2} (d^(m-1) z/dx^(m-1)) [u^m] exp[(S(u hbar d_z) - 1) x u]/(u hbar),
S(t) = sinh(t/2)/(t/2), b running over the singular points of x, y where y dx or x dy has a non-zero residue: here b = a and b = infinity.
Hypotheses: Sigma = P^1 compact, x, y meromorphic, x with simple ramification (assumed; true at all semisimple points), y unramified
(F_g^vee = 0). (ii) Extraction (elementary): (S(t) - 1) = Sum_{k>=1} t^(2k)/(4^k (2k+1)!), so [hbar^(2g-1)][u^m] keeps the products of
j factors with Sum k_i = g and m = 2g-1+j:
   F_g = (1/(2-2g)) Sum_b Res_{z=b} I_g dz,   I_g = Sum_{j=1}^{g} (1/j!) D^(2g-2+j)(z) Sum_{k_1+..+k_j=g} prod_i x^(2k_i)(z)/(4^k_i (2k_i+1)!),
D = (1/x') d/dz. (iii) z = infinity: deg D^n(z) <= 1 - n r and deg x^(2k) <= r - 2k, so every term of I_g has degree
<= 1 - (2g-2) r - 2g <= -5: Res = 0. (iv) z = a (put a = 0; y -> y - a does not change D^n(y) for n >= 1): x = v/z + regular, so
D^n(z) = (-1)^n n! z^(n+1) v^-n (1 + O(z^2)) and x^(2k) = (2k)! v z^-(2k+1) (1 + O(z^(2k+1))); every term of I_g has lowest order
exactly z^-1, so only the leading terms contribute and the regular part of x drops out:
   Res_{z=0} I_g dz = v^(2-2g) A_g,   A_g = Sum_{j>=1} (-1)^j (2g-2+j)!/j! [t^g] phi(t)^j,   phi = Sum_{k>=1} (t/4)^k/(2k+1).
(v) A_g = (2g-2)! [t^g] (1 + phi)^(1-2g) (binomial series Sum_j (2g-1)_j (-phi)^j/j!), 1 + phi = artanh(s)/s with t = 4 s^2, so
A_g = (2g-2)! 4^-g [s^1] w(s)^(1-2g), w = artanh s; Lagrange inversion (s = tanh w): [s^1] w^(1-2g) = (1-2g) [w^2g] (w coth w)
= (1-2g) 4^g B_2g/(2g)!. Hence A_g = -B_2g/(2g) and F_g = A_g v^(2-2g)/(2-2g) = B_2g/(2g(2g-2)) v^(2-2g). QED
*Sign remark.* Hock's s.3.2 prints B_2g/(2g(2-2g)) for the Harer-Zagier curve; evaluating his own formula (1.3) exactly gives
B_2g/(2g(2g-2)) (= chi(M_g), and = our TR values); the slip is the expansion "d^(m-1)z/dx^(m-1) = (-1)^m (m-1)! z^m" (correct:
(-1)^(m-1) (m-1)! z^m for x = 1/z + ..). Our sign is fixed by TR (EO conventions) and agrees with chi(M_2) = -1/240.

### 4.2 Verification (independent of the quoted theorem)
* TR, x = W: F_2 = -1/240, F_3 = 1/1008 exact (Q(zeta)); F_2..F_5 = -1/240, 1/1008, -1/1440, 1/1056 at 2^89-1 and 2^61-1
  (tr_run_fg.log, 137 s); these are c_2..c_5 = 1/1200, 1/25200, 1/180000, 1/660000 of the frozen theory, the vertex values of the
  project's genus 2-5 graph sums (NATURAL_FORM.md s.5-6), here by a different method.
* The formula of 4.1 evaluated exactly (tr_hock.py, tr_hock.log): Res_inf = 0 and Res_0 = A_g v^(2-2g) for P = z^4/4 (v = 1), z^4 (v = 4),
  z (Harer-Zagier), z^2/2 (v = 3), two random rational quartics, g = 2..6; A_g = -B_2g/(2g) for g = 2..12.
* Whole Hurwitz space (tr_generic.log): TR at all five branch points of x = z^4/4 + al z^3 + be z^2 + ga z + v/z at 2 random points
  over each prime: F_2 v^2 = -1/240, F_3 v^4 = 1/1008; family x = p^r/r + 1/p, r = 1..6: F_2 = -1/240, F_3 = 1/1008 for every r.

### 4.3 Consequences
* **Primary-sector lemma for g >= 2 (PROVED modulo the same quoted theorem + DNOPS Thm 8.1).** F_g is constant along every direction
  of H_{0;1,4} except the residue t_5 = a_5; since the family of CohFTs over the Frobenius manifold is the shifted family,
  <tau_0(d_alpha1) .. tau_0(d_alphan)>_g = d_alpha1 .. d_alphan F_g, which vanishes unless all alpha_i = 5 (H). This is item 3 of
  UNIFORM_REDUCTION s.6.5 for genus >= 2 (genus 1 not covered). Consistent with the exact tables of s.3.
* The model gap of UNIFORM_REDUCTION s.6.5 item 2 is settled for g >= 2 (modulo the quoted theorem): it is the Harer-Zagier /
  Gaussian (c = 1) free energy of the simple pole of x, independent of the polynomial part (the 4-spin part contributes nothing:
  Prop. 3.4 of Hock, and the degree count at infinity).

## 5. Behaviour as v -> 0 (task item 4)

**Homogeneity (PROVED).** p -> k p, u -> k^5 u gives x -> k^4 x, y -> k y, hence omega -> k^(5(2-2g-n)) omega; with s = u^(1/5):
omega^(u)_{g,n}(p) = s^(5(2-2g-n)) f(p/s) prod dp_l/s if omega^(1) = f(p) prod dp. Expanding omega^(1) at infinity (|p| > 1),
omega^(1) = Sum_m c_m prod dp_l/p_l^(m_l+1) (m_l >= 1; |m| = 0 mod 5 by Z/5), gives for |p_l| > |u|^(1/5) the convergent expansion
   omega^(u)_{g,n}(p) = Sum_{K>=0} u^(2-2g-n+K) Omega_K(p),   Omega_K = Sum_{|m|=5K} c_m prod dp_l/p_l^(m_l+1).
Inner region: omega^(u)(s q) = u^(2-2g-n) omega^(1)(q) exactly (pull-back), i.e. at p ~ v^(1/5) everything scales like the free energy.

**Result (VERIFIED, tr_run_v0.log, both primes):** c_m = 0 for all |m| < 5(2g-2+n), and the K = 2g-2+n term is non-zero, for
(g,n) = (0,3), (0,4), (0,5), (1,1), (1,2), (1,3), (1,4), (2,1), (2,2), (2,3), (3,1), (3,2), (4,1). Hence at fixed p != 0 every computed
omega_{g,n} is a power series in u with finite non-zero limit:
   omega_{1,1} -> -(5/8) dp/p^6,  omega_{2,1} -> (2079/128) dp/p^16,  omega_{3,1} -> (3132675/1024) dp/p^26,  omega_{4,1} -> -(45541532421/32768) dp/p^36,
   omega_{0,3} -> -(3 [1,1,3] + 4 [1,2,2] + 3 [1,3,1] + 4 [2,1,2] + 4 [2,2,1] + 3 [3,1,1]),  [m1,m2,m3] = prod dp_l/p_l^(m_l+1),
   omega_{1,2} -> (45/8)[1,9] - (21/8)[3,7] + (25/8)[5,5] - (21/8)[7,3] + (45/8)[9,1]  (x = W; (4,1) not reconstructible at 2^61-1).
   **Answer to item 4: none of the omega_{g,n} with n >= 1 is singular at fixed p; there is no leading singular part.** Singular are
   the F_g (n = 0, order v^(2-2g); their residues sit at the branch points, which collide at p = 0) and the flat correlators: e.g.
   <H>_g ~ u^(1-2g) multiplies xif_1-type functions that vanish at fixed p like u (xi of d/dt_1 = -u/(p^5 - u)), and the remaining
   singular terms cancel; at the first order this is the Euler-dilaton identity <H>_g + <tau_1(1)>_g = 0 (the dp/p^6 coefficient of
   omega_{g,1} at u = 1 is -5 (-5)^(g-1) (<H>_g + <tau_1(1)>_g): d xif_{1,0}, d xif_{0,1} both start with -5 dp/p^6). Regularity of
   omega_{g,1} at fixed p is thus a compact form of 2g-2 linear identities among genus-g one-point functions. The u^0 limits obey the
   4-spin selection rule of x = p^4/4 (n = 1: a single monomial dp/p^(5(2g-1)+1), i.e. a/4 + d = 3g - 2 - (g-1)/2), and the genus-1 limit is
   4 <tau_1(e_0)>^(4-spin)_1 d((d/dx)(-1/p))|_{x = p^4} = 4 (1/8)(-5/4) dp/p^6 = -(5/8) dp/p^6. CONJECTURAL: lim_{v->0} omega_{g,n} =
   omega_{g,n} of the Bouchard-Eynard recursion on x = p^4/4, y = p (the general collision theorem does not cover a merging pole).

## 6. What this suggests for Conjecture N (SKETCH, not established)

The frozen/vertex theory is TR on x = W_v. In the local TR of DOSS a left Givental action by the regular gauge R' replaces B by
B + Delta B (whether Delta B is global is part of the question), and R'.(M.omega^) is the Bergman-kernel deformation expansion: vertices =
model omega_{g,n}, edges = pairings against Delta B. Since the model omega_{g,n} (n >= 1) are regular at fixed p and all their
singular behaviour lives at p ~ v^(1/5), an edge pairing that can be written as a contour integral over a FIXED contour |p| = R
is automatically regular; the singular part of an edge pairing is then the residue at p = 0, the pole of x, i.e. the vanishing cycle.
This is the TR form of the plan's hypothesis "the regular gauge preserves the vanishing period": (G) would say that the dual functions
of the gauge have no residue against the omega's at p = 0. Neither the identification of the quintic with a Delta B-deformation of
this curve nor the contour argument has been checked; the vertex graph (n = 0) is covered by Theorem G.

## 7. Files (this folder) and runs

`tr_field.py` (exact Q(zeta5), F_p with zeta of any order), `tr_core.py` (series; local data tau, sigma, x' at a branch point; ModelCurve
x = kappa W; AiryCurve), `tr_family.py` (x = p^r/r + 1/p; prime search), `tr_rec.py` (recursion (2-28), symmetric and full modes; dilaton,
string, F_g), `tr_xi.py` (xi-functions, d/dx, conversion to the colour and flat d xi bases, Laplace transforms (2-31), (2-32), pairings),
`tr_frozen.py` (M_0 from the model rows; dump check; frozen correlators via the project's MVertex6), `tr_compare.py` (normalisation of
s.2.4, comparisons, Bernoulli numbers), `tr_v0.py` (expansion at infinity), `tr_hock.py` (Hock's formula, A_g), `tr_generic.py`
(generic points of H_{0;1,4}, r-family), `tr_test_airy.py` (Witten-Kontsevich), `tr_run.py` (driver: airy, rmat, exact, frozen, fg, v0).
Logs, all ending "ALL CHECKS PASSED": tr_run_{airy,rmat,exact,frozen,fg,v0}.log, tr_hock.log, tr_generic.log. Run times (one process
each, at most two at a time): 1, 22, 38, 9, 137, 25, 2, 58 s.

## 8. Caveats

* The identification TR = CohFT (DNOPS Thm 8.1, (2-35)) and Hock's Cor. 1.2 are quoted; Hock's paper is a 2025 preprint whose closed
  formulas carry a sign slip (4.1); its general theorem was checked here only through the consequences listed in 4.2.
* "R_TR = M_0 to all orders" rests on the z^8 check plus the uniqueness of the homogeneous R-matrix (standard, not re-proved).
* The comparison with MVertex6 is mod two primes; it covers psi-deficit <= 6 and the (g,n) listed (no (3,1) with d = 0, which needs
  deficit 7). Exact rational values are from the Q(zeta) run (through (3,1)) and rational reconstruction (genus 4-5 F_g).
* Item 4 is VERIFIED for 13 (g,n) only; regularity for all (g,n) is CONJECTURAL, as is the 4-spin identification of the limits.
* 4.3 uses that the CohFTs over the Frobenius manifold form the shifted family (flat derivatives of F_g = primary correlators):
  standard for Givental-Teleman CohFTs, not re-proved; genus 1 is not covered.
