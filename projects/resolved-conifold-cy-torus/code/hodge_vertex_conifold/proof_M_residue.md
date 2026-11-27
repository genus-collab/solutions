# Conjecture M via a closed residue formula (proof attempt, 2026-10-30)

Notation as in mixed_scale_findings.md and conifold_regularity_note.md section 7.
Weights w = (c, -x, -1), c = 1 + x; leaves lam_j(psi) = 1/(c - j psi); dilaton leaves T(psi) = -psi^2/(c - psi);
Xi_h(j_1..j_r) = Sum_m 1/m! < Omega_w prod lam_{j_i} prod T >_h.
Tools: no Mathematica (licence reserved); checks done in Perl with Math::BigRat (exact rationals), scripts in the scratchpad.

## 0. Plan
1. Reduce the dilaton-leaf sum to a shift of the KdV/Givental times (dilaton-shifted coordinate q_0), and rewrite Xi_h as
   a derivative of the Hodge potential F_h at a point of the "geometric" locus.
2. Use the Faber-Pandharipande / Mumford R-matrix action (Hodge CohFT = R . KW with R the Stirling series) to get a
   closed formula for the single-leaf case r = 1 at all genera, which is the user's question
   ("does a single geometric leaf of a different scale in the Theorem-C configuration vanish").
3. Try to extend to general r.

## 1. Removing the dilaton leaves: a kappa class  [proved]

Normalise c = 1 (homogeneity: Xi_h is homogeneous of degree -(3h-3+2r) in (c, w)), so w = (1, tau, -1-tau),
lam_j = 1/(1 - j psi), T = -psi^2/(1-psi).

Lemma 1.  Xi_h(j_1..j_r) = Int_{Mbar_{h,r}} Omega_w  K  prod_i 1/(1 - j_i psi_i),   K := exp(- Sum_{k>=1} kappa_k / k)
(Arbarello-Cornalba kappa classes on Mbar_{h,r}).

Proof. Every dilaton leaf is divisible by psi^2, so on Mbar_{h,r+m} one may replace psi_i (i <= r) by pi^* psi_i (the
correction divisors put i on a genus-0 bubble with a set S of forgotten points, where prod_{k in S} psi_k^{>=2} has degree
2|S| > |S| - 1 = dim). Omega_w is pulled back. Then pi_*(prod psi_{r+k}^{b_k+1}) = Sum_{sigma in S_m} prod_{cycles} kappa_{b(c)}
gives Sum_m 1/m! pi_*(prod T(psi)) = exp( Sum_k kappa_k [z^k] (-log(1 - T(z)/z)) ), and 1 - T(z)/z = 1/(1-z), so the
exponent is Sum_k kappa_k [z^k] log(1-z) = - Sum kappa_k/k.  QED

Corollary 1 (r = 1).  On Mbar_{h,1}, kappa_k = pi^*kappa_k + psi^k (pi: Mbar_{h,1} -> Mbar_h), so K = pi^*K_h (1 - psi) and
    Xi_h(j) = Int_{Mbar_{h,1}} pi^*(Omega_w K_h) (1 - psi)/(1 - j psi) = (j - 1) Sum_{b>=0} j^b Int_{Mbar_h} Omega_w K_h kappa_b
(kappa_0 = 2h - 2). So Xi_h(j)/(j-1) = P_h(j) := Int_{Mbar_h} Omega_w K_h pi_*(psi/(1 - j psi)), a polynomial of degree 3h-3.

## 2. K is the log tangent bundle  [proved; GRR]

Lemma 2.  On Mbar_{g,n}:  Lambda^vee(1) * K = c( T_{Mbar_{g,n}}(-log boundary) ),  where T(-log) = R^1 pi_* (omega_log^{-1})
(= dual of pi_*(omega (x) omega_log), the bundle of quadratic differentials with simple poles at the markings).

Proof. Mumford/Chiodo GRR for r = 1: for S = omega_log^s(-Sum a_i x_i),
  ch_m(R pi_* S) = B_{m+1}(s)/(m+1)! kappa_m - Sum_i B_{m+1}(a_i)/(m+1)! psi_i^m + (1/2) B_{m+1}/(m+1)! j_*(gamma_{m-1}),
and c(-R pi_* S) = exp(Sum_m (-1)^m (m-1)! ch_m). s = 0, a = 0 gives c(E^vee) = Lambda^vee(1). Passing to s = -1 changes only the
kappa part, by Sum_m (-1)^m [B_{m+1}(-1) - B_{m+1}(0)]/(m(m+1)) kappa_m = - Sum_m kappa_m/m  (B_{m+1}(0) - B_{m+1}(-1) = (m+1)(-1)^m).
omega_log^{-1} has negative degree on every component of a stable curve, so R^0 = 0 and R^1 pi_* omega_log^{-1} is the log
tangent bundle (equisingular deformations plus one trivial "log" direction per node).  QED

Checks against the stored exact data (hvc_which_class_out.txt, single class Lambda^vee(1), leaf scale 1):
  Xi^H_2(0) = (-4)(-2)(-1)(-6)/5760 = 1/120 = zeta(-3) = chi(M_{2,1});   Xi^H_3(0) = -1/252 = zeta(-5) = chi(M_{3,1})
(Gauss-Bonnet: Int c_top(T(-log)) = orbifold Euler characteristic of the open part, Harer-Zagier). Both agree.

Consequence (the CY family).  With u_2 = 1/tau, u_3 = -1/(1+tau), c_u(E^vee) := Sum_i c_i(E^vee) u^i = Lambda^vee(w)/w^h:
    Xi_h(j_1..j_r) = Int_{Mbar_{h,r}} c(T(-log)) c_{u_2}(E^vee) c_{u_3}(E^vee) prod_i 1/(1 - j_i psi_i).        (2.1)
The leaf scale = the weight w_1 = 1 is exactly what lets Lambda^vee(w_1)/w_1^h absorb K: the dilaton leaves at scale c turn
Lambda^vee(c) into the log tangent bundle.  The single class Lambda^vee(1) (tau -> infinity) gives the "Hurwitz" functional
    Xi^H_h(j_1..j_r) = Int_{Mbar_{h,r}} c(T(-log)) prod_i 1/(1 - j_i psi_i),   Xi^H_h(0^r) = chi(M_{h,r}).
Writing sigma = tau(1+tau): Lambda^vee(tau)Lambda^vee(-1-tau) is a polynomial in sigma of degree h; sigma^h coefficient gives
Xi^H, sigma^0 gives the lambda_g case. Genus 2 data in this form:
    5760 sigma^2 Xi_2(j) = (j-1)(j-2) [ 7 - (12 j - 28) sigma + (5 j^2 - 26 j + 24) sigma^2 ]
so each sigma-coefficient vanishes separately at j = 1, 2.

## 3. The geometric mechanism behind the range J <= 2h-3+r: a line bundle of negative degree

Put D = Sum j_i x_i and S_D := omega_log^{-1}(D) on the universal curve. From 0 -> S(-x_i) -> S -> S|_{x_i} -> 0 and
omega_log|_{x_i} = O (residue), O(k x_i)|_{x_i} = L_i^{-k} (L_i = cotangent line, c_1 = psi_i):
    c(-R pi_* S_D) = c(T(-log)) / prod_i prod_{k=1}^{j_i} (1 - k psi_i),
hence (p_0 := 1, p_j(psi) := prod_{k=1}^{j-1}(1 - k psi) for j >= 1, a polynomial of degree j - 1)
    c(T(-log)) prod_i 1/(1 - j_i psi_i) = c(-R pi_* S_D) prod_i p_{j_i}(psi_i).                                  (3.1)
rank(-R pi_* S_D) = -chi(S_D) = (3h-3+r) - J =: N - J;  deg S_D = J - (2h-2+r) on a smooth curve.
So "deg S_D < 0 on smooth curves" is EXACTLY the range J <= 2h-3+r of Conjecture M.

Lemma 3 (degree vanishing).  If J >= 1 and H^0(C, S_D|_C) = 0 for every stable curve C in Mbar_{h,r}, then
Xi^H_h(j_1..j_r) = 0.
Proof. Then R^1 pi_* S_D is locally free of rank N - J, so c_k(-R pi_* S_D) = 0 for k > N - J. In (3.1) the psi-polynomial has
degree J - r_+ (r_+ = #{i : j_i >= 1} >= 1 because J >= 1), so the degree-N part needs c_k with k >= N - J + r_+ > N - J.  QED

Jump locus.  On a smooth curve of genus h, deg S_D < 0 gives H^0 = 0 in the whole range. On a stable curve a section of S_D
vanishes on every component without markings (degree -(2g_v - 2 + n_v) < 0) and must vanish at the nodes joining the
support to the rest. For r = 1 a section can live only on the component C_0 containing x_1 (normalisation genus g0, l0 self-nodes,
e0 >= 1 nodes to other components), where the relevant degree is  j + 1 - 2(g0 + l0 + e0).  Stability of C_0 forces
g0 + l0 + e0 >= 2, hence:

Proposition 3 [proved].  For the single class Lambda^vee(1) (Hurwitz case) and every h >= 2:  Xi^H_h(1) = Xi^H_h(2) = 0.
(j = 2 is new; j = 1 is Theorem C.)  For j >= 3 the jump locus is non-empty for h >= 3: e.g. x_1 on a rational component
with exactly two nodes (degree j - 3 >= 0, a codimension-2 stratum, kernel of rank j - 2), elliptic tails with
2x_1 ~ 2p, nodal rational tails with a cross-ratio condition. For r >= 2 already rational tails carrying x_1, x_2 with
j_1 = j_2 = 2 give a jump on a whole divisor. The data (exact, genus 3: Xi^H_3(3) = Xi^H_3(4) = 0) show that the excess
contributions of the jump locus nevertheless integrate to zero; Lemma 3 does not see this.

Remark (why the Hodge class matters).  For the trivial (KW) class there is no K-absorption: Lambda^vee(1) is what turns the
dilaton leaves into c(T(-log)) and makes the Riemann-Roch count of S_D visible; this is why the zeros 2..2h-2 need the
Hodge class with leaf scale = weight (hvc_which_class_out.txt).  For the CY triple, by Mumford Lambda^vee(-1-tau) =
(-1)^h (1+tau)^{2h}/Lambda^vee(1+tau), so the extra factor is the ratio e(E^vee (x) C_tau)/e(E^vee (x) C_{tau+1}) of virtual
rank 0: the integrand is the equivariant Euler class of a K-class of virtual rank N - r_+, the same as in the Hurwitz case,
but not of an honest bundle (the denominator E^vee (x) C_{tau+1} has no partner), so Lemma 3 does not apply directly.

## 4. Beyond Lemma 3: a localized top Chern class, and Xi^H_h(3) = 0

V_2 := -R pi_* S_{2x_1} = R^1 pi_* omega_log^{-1}(2 x_1) is LOCALLY FREE of rank N - 2 on Mbar_{h,1} for every h >= 2
(no jumps for j <= 2, section 3). Polynomial division of (1 - psi)(1 - 2 psi) by (1 - j psi) and degree vanishing give,
for EVERY j,
    Xi^H_h(j) = (j-1)(j-2) E_h(j),     E_h(j) := Int_{Mbar_{h,1}} e( V_2 (x) L_1^{(x) j} ) psi_1^2 ,                    (4.1)
a polynomial identity in j (deg E_h = 3h - 4).  Check, genus 2 (V_2 of rank 2, c(V_2) = c(T(-log))/((1-psi)(1-2psi))):
E_2(j) = Int c_2(V_2)psi^2 + j Int c_1(V_2)psi^3 + j^2 Int psi^4 with Int psi^4 = 1/1152 = 5/5760 and
Int c_1(V_2) psi^3 = -<lambda_1 tau_3> - <tau_3 tau_2> + 3<tau_4> = (-12 - 29 + 15)/5760 = -26/5760, as in the data
E_2(j) = (5j^2 - 26j + 24)/5760.

Proposition 4 [proved]. Xi^H_h(3) = 0 for all h >= 3.
Proof. 0 -> S_2 -> S_3 -> L^{-3} -> 0 gives a sheaf injection L^{-3} -> V_2 (pi_* S_3 = 0: torsion free and generically 0),
i.e. a section s of V_2 (x) L^3 vanishing exactly on the jump locus Z_3 = {H^0(C, S_3) != 0} (base change for the top
cohomology). Hence e(V_2 (x) L^3) = [Z_3]^loc in A_2(Z_3) (Fulton, localized top Chern class) and Xi^H(3) = 2 Int_{[Z_3]^loc} psi_1^2.
Dimension estimate. Let V be an irreducible surface in Z_3, Gamma the dual graph of its generic point, v0 the vertex carrying
x_1. Then psi_1 is pulled back from Mbar_{g_{v0}, n_{v0}} (gluing map), and the image of V there lies in the closure of
Y_0 = {(C_0, x_1, attachment points a, self-node points q) : H^0(C_0, omega_{C_0}^{-1}((j-1)x_1 - 2a - q)) != 0}.
A non-zero section with divisor E gives a meromorphic differential eta with div(eta) = (j-1)x_1 - 2a - q - E, of degree
d_0 = j + 1 - 2(g0 + l0 + e0) >= 0, with e0 >= 1 double poles (e0 = 0 is impossible for j <= 2h-2). Meromorphic strata have
dimension 2g - 2 + n (period coordinates), so dim Y_0 <= 2 g0 - 3 + (1 + e0 + 2 l0 + d_0) = j - 1 - e0 <= j - 2.        (4.2)
For j = 3 this is <= 1 < 2, so the push-forward of [V] to Mbar_{g_{v0},n_{v0}} is zero and Int_V psi_1^2 = 0. Hence
psi_1^2 kills A_2(Z_3) and Xi^H(3) = 0.  QED
(For h = 2, j = 3 the jump locus contains the Weierstrass divisor (e0 = 0) and indeed Xi^H_2(3) != 0.)

Why it stops at j = 3: for j >= 4, S_{j-1} has jumps, R^1 pi_* S_{j-1} is not locally free, and the excess Porteous classes
of J_j^vee -> T(-log) on Z_j occur in all dimensions d <= j - 1, while (4.2) only kills d = j - 1.  General statement left:
    E_h(j) = 0 for 4 <= j <= 2h - 2   (Hurwitz case), and the whole CY family.

## 5. The gap in closed form (r = 0), interpolating Harer-Zagier and lambda_g  [data h = 2,3; identity]

With c = 1, w = (1, tau, -1-tau), sigma = tau(1+tau), eps = (-tau, 1+tau):
    F_h(q_0) = Int_{Mbar_h} c(T(-log)) c_{u2}(E^vee) c_{u3}(E^vee) = -(2h-3)! P_h(-tau, 1+tau) / sigma^h,
    P_h(e1, e2) = Sum_{a+b=h} B_{2a}(1/2) B_{2b}(1/2) e1^{2a} e2^{2b} / ((2a)!(2b)!)  = [lambda^{2h}] (e1 l/2)/sinh(e1 l/2) (e2 l/2)/sinh(e2 l/2).
Checked: h = 2: 5760 P_2 = 7 + 28 sigma + 24 sigma^2;  h = 3: 967680 P_3 = -(31 + 186 sigma + 328 sigma^2 + 160 sigma^3);
genus 4: P_4(x,1) is proportional to the stored gap numerator 381x^8+620x^6+686x^4+620x^2+381.
sigma -> infinity: [sigma^h]P_h = -(2h-1)B_{2h}/(2h)! gives chi(M_h) = B_{2h}/(2h(2h-2)) (Harer-Zagier);
sigma -> 0: -(2h-3)! b_h (-1)^h, the lambda_g value. Consistent with Xi(0) = (2-2h) F_h(q_0).

## 6. CLOSED FORMULA for one leaf (r = 1), all tau  [found 2026-10-31; verified exactly h = 1..4]

Notation: c = 1, w = (1, tau, -1-tau) = (e1 + e2, -e1, -e2) with e1 = -tau, e2 = 1 + tau (so e1 + e2 = 1 = leaf scale),
sigma = tau(1+tau) = -e1 e2,  S(x) := (x/2)/sinh(x/2) = Sum_g (-1)^g b_g x^{2g}.

    sigma^h Xi_h(j)  =  prod_{l=1}^{2h-2} (j - l)  *  [lambda^{2h}]  S(e1 lambda)^{1 + j e2}  S(e2 lambda)^{1 + j e1} .        (6.1)

Equivalently, with (2h-2)! C(j-1, 2h-2) = (2h-2)! Res_w (1+w)^{j-1} w^{-(2h-1)}:
    sigma^h Xi_h(j) = (2h-2)!  Res_{w=0} Res_{lambda=0}  (1+w)^{j-1} w^{1-2h} lambda^{-2h-1} S(e1 l)^{1+j e2} S(e2 l)^{1+j e1} dw dl.
Evidence: exact rational comparison with the stored symbolic data (hvc_which_class_out.txt, hvc_multiprobe_g4_out.txt
converted to sigma; genus 1 from Int_{Mbar_{1,1}}) at 180 points (5 values of tau, 9 of j, h = 1..4): no mismatch.
Genus 2 and the j^2 coefficient of genus 3 fixed nothing: the exponent E = e2 log S(e1 l) + e1 log S(e2 l) was READ OFF from
h = 2 and the j-linear part of h = 3 (coefficients -s/24, s(p-s)/2880, -s(p^2-ps-s^2)/181440, s = e1e2, p = e1^2+e2^2,
i.e. s(e1^{2k-1} + e2^{2k-1}) times the coefficients of log S); the j^2 term of h = 3 and ALL of genus 4 are predictions.
(Script: check_closed_formula_r1.pl in this folder, Perl Math::BigRat; run: perl check_closed_formula_r1.pl.)
Special cases: j = 0 gives the gap formula of section 5; e1 -> 0 (lambda_g point) gives b_h prod (j - l) (Norlund);
sigma -> infinity (Hurwitz) gives [lambda^{2h}] of the e1^{2h}-part: Xi^H_h(j) = prod(j-l) * coefficient, including the extra
zero j = 2h. In the c = 1 normalisation the polynomial R_h of mixed_scale_findings.md is R_h(j) = sigma^{-h} [lambda^{2h}] S(e1 l)^{1+je2} S(e2 l)^{1+je1}.

Consequence for the user's question (CONDITIONAL on (6.1); unconditional for h = 2, section 7): a single leaf lam_j (scale c/j) in the Theorem-C configuration
(Xi(j, 1, ..., 1) = j(j+1)...(j+r-2) c^{-2(r-1)} Xi(j) by the string equation) vanishes for ALL tau exactly when
j in {1, ..., 2h-2} (the bracket in (6.1) is a non-zero polynomial in tau at every integer j, its tau -> 0 value being
(-1)^h b_h != 0), and for special tau at the further real roots of the degree-h polynomial in (6.1).

### 6a. Equivalent forms of (6.1)
(i) Shifted Norlund residue (generalises Xi_lead = b_h Res_w (1+w)^{J-1}/w^{N+1}; put 1+w = e^v):
    sigma^h Xi_h(j) = (2h-2)!  Res_{lambda=0} Res_{v=0}  e^{j v} (e^{v - beta(lambda)} - 1)^{1-2h} S(e1 l) S(e2 l) lambda^{-2h-1} dv dl,
    beta(lambda) = e2 log S(e1 lambda) + e1 log S(e2 lambda)  (= log B, B = S(e1 l)^{e2} S(e2 l)^{e1}).
    At e1 = 0: beta = 0, S(e1 l)S(e2 l) = S(l), [l^{2h}]S(l) = (-1)^h b_h: the lambda_g formula.
(ii) "String-shifted lambda_g": for every leaf f(psi) (Borel transform f^(v) = Sum f_a v^a/a!),
    sigma^h <Omega_w K f>_{h,1} = [lambda^{2h}] S(e1 l)S(e2 l) * (1/b_h) <lambda_h K (e^{beta(lambda) S} f)>_{h,1},
    where S f = (f(psi) - f(0))/psi is the string derivation on leaves (S lam_j = j lam_j, f^(v) -> f^(v + beta)).
    Coefficientwise: <Omega_w K psi^a> = sigma^{-h} Sum_m [l^{2h}](S1 S2 beta^m/m!) * <lambda_h K psi^{a-m}>/b_h.
(iii) Hurwitz limit (sigma -> infinity, the single class Lambda^vee(1), i.e. c(T(-log))): with
    log H(x) := Sum_{k>=1} (2k-1) B_{2k} x^{2k} / (2k (2k)!)   (H = S e^{-x S'/S}),  and  1 + x^2 (log S)'' = S^2,
    Xi^H_h(j) = prod_{l=1}^{2h-2}(j-l) [x^{2h}] S^{j+2} e^{-j x S'/S} = prod_{l=1}^{2h-2}(j-l) * ((j - 2h)/j) [x^{2h}] H(x)^j
             = prod_{l=1}^{2h-2}(j-l) * [y^j] x(y)^{j-2h},   x = y H(x)  (Lagrange),
    and y = x/H(x) = 2 e sinh(x/2) exp(-(x/2)coth(x/2)) = e w (1+w)^{-(1+w)/w}, w = e^x - 1.
    The factor (j - 2h) is the EXTRA ZERO j = 2h of the single class (hvc_which_class_out.txt), now explained by the
    identity S^2 H^j = H^j - (x/j)(H^j)' (from 1 + x^2(log S)'' = S^2). Check h = 2: (1/j)[x^4]H^j = (5j-6)/5760.  [identity proved]
(iv) In the binomial basis Xi_h(j) = Sum_{k=0}^{h} alpha_{h,k} C(j-1, 2h-2+k): Conjecture M for r = 1 is exactly the
    absence of C(j-1, m) with m < 2h-2; (6.1) gives alpha_{h,k} = (2h-2+k)! [lambda^{2h}] S1 S2 (B-1)^k ... (expansion
    of e^{j beta} = B^j in the Newton basis of j around 2h-1 is not used below).

### 6b. What (6.1) would prove
If (6.1) holds for all h (conjecture, verified h <= 4), then for r = 1 and integers j >= 0:
Xi_h(j) = 0 iff j in {1..2h-2}  (for generic tau), i.e. Conjecture M for r = 1, all genera, all tau, with the explicit
non-vanishing factor [lambda^{2h}] S(e1 l)^{1+je2} S(e2 l)^{1+je1} (value (-1)^h b_h != 0 at tau = 0, so not identically 0 in tau).
Via the string equation this covers every Theorem-C configuration with one leaf of another integer scale:
Xi(j, 1^{r-1}) = j (j+1) ... (j+r-2) Xi(j) (c = 1).

## 7. Splitting the CY class into Hodge components; Conjecture M at genus 2 is PROVED (all r, all tau)

(a) Component decomposition [proved]. With alpha_i the Chern roots of E and u2 + u3 = 1/sigma = -u2 u3,
    c_{u2}(E^vee) c_{u3}(E^vee) = prod_i (1 - alpha_i(1 + alpha_i)/sigma) = Sum_{k=0}^h (-1/sigma)^k e_k(alpha + alpha^2),
    Xi_h(j..) = Sum_k (-1/sigma)^k Xi^{(k)}_h(j..),   Xi^{(k)} := Int c(T(-log)) e_k(alpha+alpha^2) prod 1/(1 - j_i psi_i).
    e_0 = 1 (Hurwitz), e_1 = lambda_1 (because ch_2(E) = 0), e_2 = lambda_2 - 6 ch_3(E), ..., e_h = lambda_h c(E), and
    c(T(-log)) e_h = lambda_h c(E) c(E^vee) K = lambda_h K (Mumford): the k = h component is the lambda_g case (Norlund, proved).
    Conjecture M (r = 1) <=> every component vanishes at j = 1..2h-2 (necessary since Xi is a polynomial in 1/sigma).

(b) Lemma 5 [proved]. On Mbar_{h,1}:  psi * e(V_2) = -(1/2) pi^* c_{3h-3}(T_{Mbar_h}(-log)).
    Proof: V_2 = pi^*T_h(-log) / L^{-2} (the Schiffer line; no jumps for S_{2x_1}), so c(V_2)(1 - 2psi) = pi^*c(T_h); take the
    degree 3h-3 part and use c_{3h-3}(V_2) = 0 (rank 3h - 4).  QED
    Consequence: Theta_2 := pi_*( c(V_2)(1 - psi) ) = c(T_h(-log)) Sum_b 2^b kappa_b  has NO components in degrees 3h-3 and 3h-4:
    [Theta_2]_{3h-3} = pi_*(c_{3h-2}(V_2) - psi c_{3h-3}(V_2)) = 0,  [Theta_2]_{3h-4} = -pi_*(psi e(V_2)) = (1/2) pi_* pi^* c_top = 0.
    (Tautological relations Sum_b 2^b c_{3h-3-b}(T_h) kappa_b = 0 and Sum_b 2^b c_{3h-4-b}(T_h) kappa_b = 0 on Mbar_h.)

(c) Proposition 6 [proved]. For all h >= 2:  Xi^{(1)}_h(2) = Int_{Mbar_h} lambda_1 Theta_2 = 0  (the lambda_1-component at j = 2).
    (Xi^{(k)}(j) = (j - 1) Int_{Mbar_h} e_k Theta_j by Corollary 1; at j = 2 only [Theta_2]_{3h-4} pairs with the degree-1 class lambda_1.)

(d) Theorem 7 [proved]. Genus 2, every tau, every r >= 1: Xi_2(j_1..j_r) = 0 whenever the j_i >= 0 are integers with
    1 <= J <= 1 + r.
    Proof. r = 1: j = 1 is Theorem C; j = 2: Xi_2(2) = Xi^{(0)}(2) - Xi^{(1)}(2)/sigma + Xi^{(2)}(2)/sigma^2 with Xi^{(0)}(2) = 0
    (Proposition 3), Xi^{(1)}(2) = 0 (Proposition 6), Xi^{(2)}(2) = b_2 (2-1)(2-2) = 0 (lambda_g).  r >= 2: any lattice point in the
    range with all j_i >= 2 has 2r <= J <= 2h-3+r, i.e. r <= 2h-3 = 1; so every point with r >= 2 has some j_i in {0, 1}, and
    the string/dilaton relations Xi(1, j') = J' Xi(j'), Xi(0, j') = (D' - 2h + 2) Xi(j') reduce it to a point in the range with
    r - 1 leaves, or produce the vanishing factor J' = 0 / D' = 2h - 2 at the boundary of the range.  QED
    In general genus the same reduction shows: Conjecture M (vanishing direction) follows from its "irreducible" points
    (all j_i >= 2, so r <= 2h-3 and 2r <= J <= 2h-3+r).

(e) Status of the irreducible points for h >= 3 (r = 1, j in 2..2h-2; components k = 0..h):
    proved: (j = 2, k in {0, 1, h}), (j = 3, k in {0, h}) for h >= 3, (all j, k = h);  open: all other (j, k), and r >= 2.
    Genus 3 open list: r = 1: (2, k=2), (3, k=1), (3, k=2), (4, k=0), (4, k=1), (4, k=2); r = 2: (2,2), (2,3); r = 3: (2,2,2).

## 8. Two leaves (r = 2): what the data say about a closed formula  [exploratory]

Data: genus 2, 3 symbolic r = 2 (hvc_multiprobe_out.txt), converted to sigma (scratchpad polylib.pl, g2r2*.pl, r2guess.pl).
* The J-only candidate  sigma^h Xi(j1,j2) =? prod_{l=1}^{2h-1}(J - l) [l^{2h}] S1 S2 B^J  is exact on j1 = 0 and j2 = 0 (an
  integration by parts in v, (e^u-1)^{-2h} = -(1/(2h-1)) d/du (e^u-1)^{1-2h} - (e^u-1)^{1-2h}), but on j1 = 1 it is off by
  exactly one factor B (string equation): no formula depending on J alone can satisfy both string and dilaton unless the
  r = 1 bracket is constant in j, which happens only at sigma = 0 (lambda_g). Consistent with "does not factor through J".
* "Flow only the deviation from q_0" (Psi_h[rho] = sigma^{-h}[l^{2h}] S1S2 Psi^{lambda_g}_h[rho_0 + e^{beta S}(rho - rho_0)]/b_h,
  which reproduces (6.1) at first order) FAILS at second order: the residual is j1(j1-1)j2(j2-1) * (non-zero, vanishing at sigma = 0).
* The two-leaf functional is equivalent to the Hessian H(j1,j2) = <Omega K Pt_{j1} Pt_{j2}>, Pt_j = lam_j - j lam_1 + (j-1) (zero for
  j in {0,1}), of the IZ potential Psi_h(x_2, x_3, ...) at x = (-1,-1,...):
  H = Xi(j1,j2) - [j1 + (2h-1)(j2-1)] Xi(j1) - [j2 + (2h-1)(j1-1)] Xi(j2) + (2h-1)(2h-2)(j1-1)(j2-1) F_h(q_0).
  A closed formula for r >= 2 needs second-order information not contained in (6.1). Not found.

## 9. Status (2026-10-31)

PROVED (all h >= 2 unless stated, all tau):
  L1  dilaton leaves = K = exp(-Sum kappa_k/k);  r = 1: Xi_h(j) = (j-1) Int_{Mbar_h} Omega_w K kappa(j), kappa(j) = Sum j^b kappa_b.
  L2  Lambda^vee(w1) K = c(T(-log)) when the dilaton/leaf scale equals w1; Xi = Int c(T(-log)) H_tau prod 1/(1 - j_i psi_i).
  L3  line-bundle mechanism: the range J <= 2h-3+r is exactly deg omega_log^{-1}(Sum j_i x_i) < 0; honest-bundle vanishing.
  P3  Hurwitz component Xi^{(0)}_h(2) = 0;  P4  Xi^{(0)}_h(3) = 0 (h >= 3, localized top Chern class + strata dimension bound).
  P6  lambda_1-component Xi^{(1)}_h(2) = 0;  lambda_g component Xi^{(h)}_h(j) = b_h prod(j - l) (Norlund).
  T7  Conjecture M (vanishing direction) at genus 2 for all r; reduction of general genus to irreducible points.
  identity: Hurwitz extra zero j = 2h follows from (6.1)-Hurwitz via 1 + x^2 (log S)'' = S^2.
CONJECTURE (exact evidence h = 1..4, genus 4 a parameter-free prediction): the closed formula (6.1); it implies Conjecture M
  for r = 1 in both directions (generic tau) and gives the gap formula F_h(q_0) = -(2h-3)! P_h / sigma^h.
OPEN: proof of (6.1) (natural routes: Marino-Vafa with one j-leg and winding-1 legs + the IZ/Lambert evaluation at s = -1,
  where Q_h(s,j) = <Omega K_s lam_j>, K_s = exp(Sum ((1+s)^k - 1) kappa_k/k), interpolates between the MV one-leg formula (s = 0)
  and (6.1) (s = -1); or an inductive excess-intersection computation on the jump loci Z_j); the components 1 <= k <= h-1 for
  j >= 2 (except k = 1, j = 2) and k = 0 for 4 <= j <= 2h-2; any closed formula for r >= 2.

### 9a. All-genus consistency checks of (6.1) [proved identities]
1. tau -> 0: lambda_g/Norlund, all h.   2. sigma -> infinity, j = 0: (2h-2)! [x^{2h}] S(x)^2 = -B_{2h}/(2h) = chi(M_{h,1})
   (Harer-Zagier; S^2 = Sum (1-2n) B_{2n} x^{2n}/(2n)!), all h.
3. Top j-coefficient j^{3h-2}: sigma^{-h} [l^{2h}] beta^h/h! = 1/(24^h h!) = <tau_{3h-2}>_h, all h.
4. Next coefficient j^{3h-3}: geometric side -(1 + 1/sigma)<lambda_1 tau_{3h-3}> - <tau_2 tau_{3h-3}> with
   <tau_2 tau_{3g-3}>_g = (5 + 12 g(g-1))/(5 24^g g!) (DVV with tau_2 removed: ((6g-3)(6g-5) + 12g)/(15 24^g g!)) and, from Mumford's lambda_1 = (kappa_1 - psi + delta)/12,
   <lambda_1 tau_{3g-3}>_g = g(g+4)/(5 24^g g!);  (6.1) gives (with s = -sigma)
   5 24^h h! sigma^h [j^{3h-3}]Xi / (-s)^h = h(h+4)/s - (13h^2 - 8h + 5)  on BOTH sides, for all h.
