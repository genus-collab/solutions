# Conjecture M via the single-Hodge (Hurwitz) endpoint and topological recursion

Status file of the proof attempt (2026-10-30). Appended as the work proceeds; see the final section "Summary of what is proved".

Notation as in mixed_scale_findings.md: Omega_w = prod_a Lambda^vee(w_a)/w_a^h, leaves lam_j(psi) = 1/(c - j psi), dilaton leaves
T(psi) = -psi^2/(c - psi), Xi_h(j_1..j_r) = Sum_m 1/m! < Omega_w prod lam_{j_i}(psi_i) prod T(psi_{r+k}) >_h.
Single class: Omega^(1) = Lambda^vee(1) = Sum_i (-1)^i lambda_i (u = 1), leaf scale c = 1.

## 1. Three reformulations (single class Omega^(1) = Lambda^vee(1), c = 1; the same works for Omega_w with general c)

**1a. kappa form (proved).** The dilaton leaves carry psi^(>=2), T(psi) = Sum_{a>=1} s_a psi^(a+1) with s_a = -c^(-a). For forgotten points carrying
psi^(>=2) the psi classes of the remaining points may be replaced by pull-backs (the correction divisors carry psi_{forgotten} restricted to an
M_{0,3} bubble, which is 0), and pi_*(prod psi_{n+i}^(a_i+1)) = Sum_{sigma in S_m} prod_{cycles} kappa_{sum of a over the cycle}. The exponential
formula gives Sum_m 1/m! pi_{m*} prod T(psi_{n+i}) = exp(Sum_k p_k kappa_k) with Sum_k p_k u^k = -log(1 - S(u)), S(u) = Sum s_a u^a = -(u/c)/(1-u/c),
so p_k = -c^(-k)/k. Lambda classes pull back, hence
    Xi_h(j_1..j_r) = < Omega_w exp(-Sum_k kappa_k/(k c^k)) prod_i 1/(c - j_i psi_i) >_{h,r}.

**1b. Lambert / one-point-function form (proved, tautological but useful).** Put x = z e^(-z) (Lambert curve), eta = 1/(1-z),
xi_0 = z/(1-z) = Sum_{mu>=1} mu^mu/mu! x^mu, xi_k = D^k xi_0 with D = x d/dx = (z/(1-z)) d/dz (xi_k is a polynomial of degree 2k+1 in eta, vanishing at
z = 0). Since [x^mu] xi_k = mu^(mu+k)/mu!, the "Lambert one-point function at the base point"
    Omega_h(z) := Sum_k Xi_{h,k} xi_k(z),   Xi_{h,k} := Sum_m 1/m! < Omega^(1) psi_1^k T^m >_h,
satisfies [x^j] Omega_h = (j^j/j!) Xi_h(j). Hence (r = 1)
    Xi_h(j) = 0 for j = 1..2h-2   <=>   Omega_h(z) = O(x^(2h-1)) = O(z^(2h-1)) at z = 0,
and the extra zero j = 2h says that the next coefficient [x^(2h)] also vanishes. Example: h = 1, Xi_1(j) = (j-2)/24 (computed by hand:
<(1-lambda_1)(1 + j psi)>_{1,1} = (j-1)/24 and the single T-term gives -1/24), Omega_1 = (xi_1 - 2 xi_0)/24 = (eta-1)(eta^2-2)/24 = -x/24 + 0 x^2 + O(x^3).
For multi-probes: [x_1^mu_1 ... x_r^mu_r] of Omega_{h,r}(z_1..z_r) := Sum <Omega psi^k ... T^m> prod xi_{k_i}(z_i) equals prod(mu_i^mu_i/mu_i!) Xi_h(mu),
so Conjecture M (=> direction) says: Omega_{h,r} has no monomial of total x-degree 1..2h-3+r.

**1c. Hurwitz form (proved, formal).** By string + dilaton (F_h(t + a) = F_h(t with leaves c_mu -> c_mu e^(a mu)); F^s(q) = lim_{L->oo} L^(2h-2) F_h(L q)),
    Xi_h(j) = lim_{L -> oo} L^(2h-1) e^(jL) W_j(-L e^L),   W_j(X) := Sum_l X^l/l! < Lambda^vee(1) (1-j psi_0)^(-1) prod_{i=1}^l (1-psi_i)^(-1) >_{h,1+l},
and by ELSV (j^j/j!) W_j(X) = d/dp_j of the genus-h simple Hurwitz potential at p = X delta_1 (covers with one j-cycle and l = d - j unramified
points over infinity). The base point is X -> -oo along the Lambert curve: X = T e^(-T) with T = T(X) = -L. So the vanishing is a statement about
the behaviour of W_j at the far end T -> -oo of the Lambert curve.

## 2. The base point is a spectral curve: the dictionary (derived, checked at genus 1)

The Hurwitz curve (Bouchard-Marino / Eynard-Mulase-Safnuk): X = log x = log z - z, y = z, omega_{0,1} = y dX, B = dz dz'/(z-z')^2; its
Eynard-Orantin correlators are omega_{g,n} = d_1..d_n < Lambda^vee(1) prod xi_{k_i}(z_i) psi_i^{k_i} > at the background t = 0 (q = -psi).
The background is a LINEAR function of omega_{0,1} (Laplace transform of dy at the branch point z = 1; the Hodge class, i.e. the R-matrix,
depends only on X and B). The leaf directions are dual, by the Eynard-Orantin variational formula (delta omega_{0,1} = Int Lambda B <=>
delta F_g = Int Lambda omega_{g,1}), to the following Lambda's:
  lam_mu (mu >= 1):  Lambda = (mu!/mu^(mu+1)) x^(-mu) on a small circle around z = 0 (because (mu!/mu^(mu+1)) Res_{z=0} x^(-mu) d xi_k = mu^k);
                     for mu = 1: delta omega_{0,1} = Res_{z'=0} (e^z'/z') B(z, z') = dz/z^2;
  lam_0 = tau_0:     xi_k(z = oo) = -delta_{k,0}, so <tau_0> = -Int_0^oo omega_{g,1}: Lambda = -1 on [0, oo], delta omega_{0,1} = dz/z.
Hence the base point q_0 = lam_0 - lam_1 corresponds to omega_{0,1} = dz/z - dz/z^2 = (z-1) dz/z^2, i.e. to
    THE BASE CURVE  C:  x = z e^(-z) (same X as Hurwitz),  y = -1/z.
(Separately lam_0 and lam_1 give y with a pole at the branch point z = 1 -- matching the fact that q = lam_0 and q = -lam_1 are singular points of
F^s, see 1c -- and the poles cancel exactly in the combination q_0.) Equivalently, by y -> y + f(x) (f = 1/x - 1): y ~ (e^z - 1 - z)/z.

**Check (genus 1, exact, Perl/BigRat EO recursion, scratch tr11.pl).** omega_{1,1} = Sum_k c_k du/u^(k+1), u = z - 1:
  Hurwitz curve: 24 (c_1, c_2, c_3) = (-1, 2, 3)  = d[(xi_1 - xi_0)/24]   (ELSV: <Lambda^vee tau_0>_1 = -1/24, <tau_1>_1 = 1/24);
  base curve C:  24 (c_1, c_2, c_3) = (-2, 2, 3)  = d[(xi_1 - 2 xi_0)/24], i.e. Xi_1(j) = (j-2)/24 exactly as computed in 1b.

**Consequences.** (i) [x_1^mu_1..x_n^mu_n] of the C-correlators F^C_{g,n} = Int omega^C_{g,n} equals prod(mu_i^mu_i/mu_i!) Xi_g(mu). So Conjecture M
(vanishing direction, single class) <=> F^C_{g,n} has no monomial of total degree < 2g-2+n at z = 0.
(ii) r = 1 in residue form: with v = 1/z (= -y), omega^C_{g,1} has poles only at v = 1, so its expansion at v = oo (z = 0) starts at v^(-2g) dv iff
    Res_{v=1} v^N omega^C_{g,1} = 0 for N = 0, 1, ..., 2g-2,  i.e.  Res_{branch point} y^N omega^C_{g,1} = 0  (N <= 2g-2).
N = 0 is automatic (residue-free) and N = 1 is the Eynard-Orantin "string" identity Sum Res y omega_{g,1} = 0 (g >= 2): this is Theorem C.
The content of the conjecture is that the higher y-moments of omega_{g,1} vanish up to N = 2g-2 on this particular curve.

**Check (genus 2 and 3, Eynard-Orantin recursion in floating point, scratch trf.pl).** TR on C gives
  g = 2: Xi_2(j) = 0, 0, -0.003125, 0, 0.0395833, 0.1666667, 0.453125 (j = 1..7)  = (j-1)(j-2)(j-4)(5j-6)/5760 to 1e-13;
  g = 3: Xi_3(j) = 0 (j = 1..4), -0.0026868, 0, 0.1132192, 0.7592593, 3.092014 (j = 5..9) = the hvc_which_class value to 1e-7.
So the dictionary "single-Hodge base point = EO recursion on C: x = z e^(-z), y = -1/z" is confirmed through genus 3 (the extra zero j = 2h included).
(It is a theorem modulo the standard facts used: ELSV/EMS for the Hurwitz curve, linearity of the background in y (DOSS), and the
variational formula; the genus 1-3 agreement is a check of the normalisations.)

**The CY family.** The same derivation with the framed-vertex (Marino-Vafa / BKMP, proved by L. Chen and J. Zhou) curve
X = log w - (tau+1) log(1+w), whose Laplace leaves are xi_0 = D log(1+w) = w/(1 - tau w) = Sum_mu binom(mu(tau+1)-1, mu-1) x^mu (D = d/dX,
branch point w = 1/tau, xi_k(oo) = -delta_{k0}/tau), gives L(lam_1) = dw/w^2, L(lam_0) = tau dw/w, hence the base curve
    C_tau:  x = w/(1+w)^(tau+1),  y = -1/w   (equivalently, with v = 1/w:  x = v^tau/(1+v)^(tau+1),  y = -v),
and C_tau -> C as tau -> oo after w = z/tau. For tau = 1 (weights (1,1,-2)) y = -1/w ~ w mod functions of x, so C_1 is the curve
1/x = w + 2 + 1/w, y = w of the Orlov-Scherbin model with r(c) = (1 + hbar c)^2, i.e. bipartite maps / dessins (Kazarian-Zograf type curve).

**Check (CY family, genus 2, scratch trtau.pl).** EO recursion on C_tau at tau = 1 and tau = 2 reproduces the hvc_which_class "CY triple, leaf
scale = w1" polynomial up to the overall normalisation (tau(tau+1))^(1-g): e.g. tau = 1: Xi_2(3..7) = -0.001953125, -0.004296875, 0.0049479,
0.0481771, 0.1582031 = (1/2) x [(j-1)(j-2)(20j^2-128j+159)/23040]; tau = 2: ratio 1/6. Zeros at j = 1, 2 only (the extra zero j = 2h is a
tau = oo phenomenon). [Side remark: for tau = 1 the base curve coincides with the Orlov-Scherbin curve of r(c) = (1+hbar c)^2, p~ = delta_1
(bipartite maps), but the OS hook formula for that tau function does NOT reproduce the TR correlators (it vanishes for mu <= 2g, while
Xi_2(3) != 0), so the naive OS/TR identification fails here; not used below.]

## 3. THEOREM (tau = 1): Conjecture M, vanishing direction, all genera and all r, via a global involution

Work on C_1 in the coordinate v = 1/w:  X = log v - 2 log(1+v) (x = v/(1+v)^2),  y = -v,  B = dv dv'/(v-v')^2.
Facts: dX = (1-v) dv/(v(1+v)); the only zero is the branch point v = 1; x(1/v) = x(v), so the local involution near v = 1 is the GLOBAL
Moebius involution sigma(v) = 1/v (fixed points v = 1 and v = -1; at v = -1, x has a double pole and dX a simple pole).
The two zeros of x, v = oo and v = 0, are exchanged by sigma; t = 1/v is a local coordinate at oo, and at 0 the coordinate is v = sigma(1/t).
Since y = -v, y(v) - y(sigma v) = 1/v - v and the EO kernel simplifies to
    K(v_0, v) = - dv_0 / ( 2 (v_0 - v)(v_0 - sigma(v)) dX(v) ) = -(dv_0/2) Sum_{m>=0} v_0^(-m-2) h_m(v, sigma v)/dX(v),
h_m(a, b) = Sum_{i+k=m} a^i b^k. With t_0 = 1/v_0, v_0^(-m-2) dv_0 = -t_0^(m+1) dt_0/t_0, so
  (*)  [t_0^(m+1) dt_0/t_0] omega_{g,n+1}(v_0, J) = (1/2) Res_{v=1} f_m(v, J),   f_m := h_m(v, sigma v) Brk_{g,n+1}(v, J) / dX(v),
where Brk = omega_{g-1,n+2}(v, sigma v, J) + Sum' omega_{h,1+|I|}(v, I) omega_{g-h,1+|J-I|}(sigma v, J-I) (omega_{0,2} terms included, omega_{0,1} excluded).

**Anti-invariance.** omega_{g,n} = d_1..d_n Sum c_k prod xi_{k_i}(v_i) with xi_0 = 1/(v - 1) (here tau = 1) and xi_k = D^k xi_0, D = d/dX.
Since x o sigma = x, D commutes with sigma^*, and xi_0 o sigma = -1 - xi_0, so sigma^* d xi_k = -d xi_k for every k. Hence omega_{g,n}
changes sign under sigma applied to ANY ONE argument. Consequently the expansion of omega_{g,n} with some arguments near v = 0 equals (up to
sign) the expansion with all arguments near v = oo, in the coordinates t_i (t_i = 1/v_i, resp. t_i = v_i).

**Degree.** Say a multidifferential near the points {0, oo} has degree >= d if, written as G(t_1..t_n) prod dt_i/t_i, every monomial of the
Laurent expansion of G has total degree >= d. (x ~ t at both points, so "degree" = total x-degree.) P(g,n): omega_{g,n} has degree
>= 2g-2+n. P(g,n) is automatic for g <= 1 (every xi_k vanishes at t = 0, so each variable contributes >= 1).

**Theorem 3.1 (tau = 1).** Assume the dictionary of section 2 for tau = 1 (base point of the triple (1,1,-2) <=> EO recursion on C_1). Then for
all g >= 2, n >= 1, omega^{C_1}_{g,n} has degree >= 2g-2+n at the zeros {0, oo} of x. Equivalently Xi_g(mu_1..mu_n) = 0 whenever all
mu_i >= 1 and |mu| <= 2g-3+n; with the string/dilaton relations (a) this is the full vanishing direction of Conjecture M for the triple (1,1,-2)
(and, by w2 <-> w3, for tau = -2).

*Proof (induction on 2g-2+n; P(g,n) is automatic for g <= 1).* Let g >= 2 and assume P for all smaller 2g'-2+n'. By anti-invariance it
suffices to expand with all of J = (v_1..v_n) near oo (t_j = 1/v_j small). Fix m >= 0 and consider f_m(v, J) as a meromorphic 1-form in v on
P^1. Its poles: v = 1 (branch point), v = 0 and v = oo (from the Laurent polynomial h_m(v, 1/v) = Sum_{i+k=m} v^(i-k)), v = v_j and
v = sigma(v_j) (from the omega_{0,2} terms). It has NO pole at the second fixed point v = -1 of sigma: for g >= 2 every term of Brk is regular
there (omega_{0,2}(v, sigma v), the only term singular at v = -1, occurs only for (g,n) = (1,0)), h_m is regular and 1/dX vanishes there.
Moreover sigma^* f_m = f_m (h_m symmetric, Brk symmetric, X o sigma = X), so Res_0 f_m = Res_oo f_m and Res_{sigma(v_j)} f_m = Res_{v_j} f_m. Hence
    Res_{v=1} f_m = -2 Res_{v=oo} f_m - 2 Sum_j Res_{v=v_j} f_m.
Degrees: at v = oo put t = 1/v; then sigma(v) = t is the local coordinate of the other zero 0 of x, h_m has degree >= -m, and
1/dX = (t/dt)(1+t)/(1-t). Every term of Brk has degree >= 2g-2+n in (t, t, t_J) (the stable factors by P and anti-invariance -- genus-0
factors even have degree >= their number of points -- and omega_{0,2}(v, v_j) = (t t_j/(t-t_j)^2)(dt/t)(dt_j/t_j) is a degree-0 Laurent series
for |t| < |t_j|, omega_{0,2}(sigma v, v_j) has degree >= 2). Thus Res_{t=0} f_m, the t^0 coefficient of h_m Brk (1+t)/(1-t)/(dt/t)^2, has
t_J-degree >= 2g-2+n-m. At v = v_j: Res_{v=v_j} B(v,v_j) phi(v) = d_{v_j} phi(v_j) with phi = h_m(v,sigma v) omega_{g,n}(sigma v, J-j)/dX(v),
which again has degree >= 2g-2+n-m (d_{v_j} preserves degree). By (*), the coefficient of t_0^(m+1) dt_0/t_0 in omega_{g,n+1} has t_J-degree
>= 2g-2+n-m for every m, i.e. omega_{g,n+1} has joint degree >= 2g-1+n. QED

Remarks. (1) For (g,n) = (1,0) the residue at v = -1 is nonzero (Res_{-1} of h_0 B(v,1/v)/dX = v dv/((v-1)^3 (v+1)) is 1/8): this is the
genus-1 exceptional term (Xi_1(1) != 0, cf. A_1 in Theorem C). (2) The proof gives exactly the range |mu| <= 2g-3+n, i.e. it explains the
number 2h-3+r: it is "number of points of x = 0 seen by Brk" bookkeeping, deg Brk >= 2g-2+n against h_m of degree -m.
(3) The proof uses only: y linear in the global coordinate v (so that K has the product form), a global involution sigma with X o sigma = X
exchanging two simple zeros of x, and no other poles of f_m. It is NOT specific to Hodge integrals once the dictionary is granted.

**Variant tau = -1/2 (triple (1,-1/2,-1/2) ~ (-2,1,1), leaf at the distinct weight).** x^2 = 1/(v(1+v)), sigma(v) = -1-v is again a global
involution with sigma^* dX = dX, fixed points v = -1/2 (branch point) and v = oo (the zero of x), and x o sigma = -x near oo. The same residue
argument (now both arguments of Brk sit at the same point oo, and h_m(v,-1-v) ~ v^m) proves degree >= 2g-2+n for omega^{C_{-1/2}}_{g,n}. But
here omega is odd in every x_i and the leaf normalisation c_mu = binom(mu(tau+1)-1, mu-1) vanishes for even mu at tau = -1/2, so this only
proves Xi_g(mu; -1/2) = 0 in the range when all mu_i are odd (no information on even parts).

**Why tau = 1 (and -2, -1/2) only.** A global involution with X o sigma = X exists iff x has degree 2, i.e. iff two of the three weights
coincide. For other tau the argument fails for a real reason, not a technical one: sigma (continued from the branch point) maps a neighbourhood of
the good zero v = oo of x to the OTHER zero v = 0, and there the correlators do NOT vanish to high order for tau != 1. Numerically (trtau.pl,
genus 2): omega_{2,1}/d nu at w = oo (nu = 1/w = v -> 0) starts 0, 0, 0.0586 for tau = 1 (degree 3, as the theorem says) but -8.7e-5, ... for
tau = 2 and -3.2e-5, ... for tau = 3 (degree 1). Local reason: at v = oo, x ~ 1/v and y = -v has a pole ("good" zero); at v = 0, x ~ v^tau
and y is regular (a Hurwitz-type zero, like z = 0 of the ordinary Lambert curve, where one-part Hurwitz numbers do not vanish). For tau = 1 the
second zero is good as well, because y = -v ~ 1/v mod functions of x (v + 1/v = 1/x - 2).

## 4. Structure of the CY family in s = tau(1+tau): reduction to h+1 functionals (proved)

With Chern roots x_i of E and w = (1, tau, -1-tau): (tau - x)(-(1+tau) - x) = -(s - x - x^2), so
    Omega_w = Lambda^vee(1) prod_{i=1}^h (1 - (x_i + x_i^2)/s) = Sum_{k=0}^{h} (-1/s)^k Lambda^vee(1) e_k(x + x^2),   s = tau(1+tau) = -w2 w3.
Hence Xi_h(mu; s) = Sum_{k=0}^h (-s)^(-k) Xi^(k)_h(mu), Xi^(k) := < Lambda^vee(1) e_k(x+x^2) exp(-Sum kappa_m/m) prod lam_{mu_i} > (exactly a
polynomial of degree <= h in 1/s. Genus 2 data: Xi_2 = (j-1)(j-2)[(5j^2-26j+24) + (28-12j)/s + 7/s^2]/5760, so each Xi^(k)_2 vanishes at j = 1, 2 separately).
  k = 0: the single class (s = oo, tau = oo);  k = h: e_h(x+x^2) = lambda_h c(E), and lambda_h c(E) Lambda^vee(1) = lambda_h (Mumford), i.e. the
  lambda_g point, PROVED (Norlund formula, mixed_scale_findings.md);  s = 2 (tau = 1): Sum_k (-1/2)^k Xi^(k) = 0 in the range, PROVED above
  (given the dictionary);  s = -1/4 (tau = -1/2): Sum_k 4^k Xi^(k) = 0 in the range for all-odd mu.
So Conjecture M for the whole CY family at genus h is equivalent to the vanishing (in the range) of the h+1 functionals Xi^(0..h); for h = 2
the three constraints are: Xi^(2) = 0 (proved), Xi^(0) - Xi^(1)/2 = 0 (proved, tau = 1), and the single-class statement Xi^(0) = 0 (open),
so the single class is the only missing input at genus 2 (which is of course also verified directly).

## 5. The single class (tau = oo): what is established, and what is missing

* Dictionary (section 2): Xi^(0)_h = single class <=> EO recursion on C: x = z e^(-z), y = -1/z, expanded at z = 0 (checked through genus 3,
  including the extra zero j = 2h). In v = 1/z: X = -log v - 1/v, y = -v, branch point v = 1; the kernel has the product form of section 3
  and (*) holds verbatim with t = z: [z_0^(m+1) dz_0/z_0] omega_{g,n+1} = (1/2) Res_{z=1} h_m(1/z, 1/sigma(z)) Brk(z, sigma z)/dX(z).
  So r = 1 is the statement Res_{z=1} y^N omega_{g,1} = 0 for N <= 2g-2 (N = 1 = Theorem C = EO string identity).
* Why section 3 does not apply: the Lambert involution sigma continued from z = 1 sends z -> 0 (the good zero of x) to z' -> +oo (the essential
  singularity of x, where xi_k(z') ~ polynomial in 1/(1-z'), i.e. only 1/log(1/x) decay); sigma has a logarithmic branch point at z = 0 and
  square-root branch points at the complex preimages of x = 1/e. The "residue at the other zero" of section 3 becomes an integral along a cut
  with no degree control.
* Proved by cut-and-join (formal, via 1c): Xi_g(2) = 0 for every g >= 1 (single class). On the line p = X delta_1 the genus-g CJ equation
  reduces exactly to (2g-2 + 2X d/dX) H_g = X^2 dH_g/dp_2 (only the join of two unramified sheets survives); with X = -L e^L, the left side is
  O(L^(2-2g)) while the right side is 2 L^(3-2g) Xi_g(2) + O(L^(2-2g)), so Xi_g(2) = 0 (for g = 1 the left side is O(1)). For g = 1 this is the
  extra zero j = 2 = 2h of Xi_1(j) = (j-2)/24. (Rigour: needs the analytic continuation of the p_1-specialised Hurwitz series to X -> -oo along
  the Lambert parametrisation, which holds because they are polynomial in 1/(1-T); I have not written this out.) The same CJ equation with
  probes gives, at leading order in L, relations inserting a leaf of size 2 (Xi_g(mu, 2) in terms of Xi_g(mu + e_k), merged probes and lower genus);
  they respect the degree pattern of Conjecture M but do not close (they never produce probes of size >= 3 from smaller ones), so they do not
  prove the general case.
* Scaling heuristic (not a proof; it predicts exactly the range). Res_{y=oo} y^N omega_{g,1} is the response of F_g to X -> X + c y^(N-1);
  the new branch points come in from y = oo where X ~ -log y + c y^(N-1); rescaling y gives a contribution c^((2g-2)/(N-1)) (homogeneity
  omega_{g,n}[x, lambda y] = lambda^(2-2g-n) omega_{g,n}), which is o(c) iff N <= 2g-2. The same heuristic applies to every C_tau because only the
  local behaviour X ~ -log y, y -> oo at the good zero enters; Theorem 3.1 is the case where it can be made rigorous by a global involution.

## 6. The first non-vanishing layer (observation, for the <= direction)

In the proof of Theorem 3.1 the minimal-degree part of omega_{g,n+1} is produced only by the minimal-degree parts of Brk (only t^(-m) of h_m
and degree exactly 2g-2+n of Brk contribute), so the leading homogeneous parts Q_{g,n} (degree 2g-2+n, i.e. Xi_g(mu) with |mu| = 2g-2+n)
satisfy a CLOSED homogeneous recursion; e.g. for r = 1, L_g = Q_{g-1,2}(1,1) + Sum_{h=1}^{g-1} L_h L_{g-h} (signs from anti-invariance).
The degree 2g-2+n equals that of the Bessel (BGW, Norbury Theta-class) correlators, but the layer is not literally Bessel: at tau = 1,
omega_{2,2} has both the (1,3) and the (2,2) coefficient non-zero (-0.176, -0.117; Bessel would have odd exponents only). Solving this
homogeneous recursion would give non-vanishing at J = 2h-2+r for tau = 1. Values (tau = 1, TR normalisation):
[x^1]F_{1,1} = -1/16, Xi_2(3) = -1/512, Xi_3(5) = -3/4096 (data: -1/256 and -3/1024 before the factor 2^(g-1)). Not proved.

## 7. Summary of what is proved (precise)

1. (Proved) kappa form Xi_h = < Omega_w exp(-Sum kappa_k/(k c^k)) prod 1/(c - j_i psi_i) >; Lambert/one-point form (r = 1 <=> the Laplace one-point
   function at the base point vanishes to order 2h-1 at x = 0); Hurwitz form Xi_h(j) = lim L^(2h-1) e^(jL) W_j(-L e^L) (T -> -oo end of the
   Lambert curve); exact decomposition Xi_h(mu; s) = Sum_{k=0}^h (-s)^(-k) Xi^(k)_h(mu), s = tau(1+tau).
2. (Derived from standard theorems -- ELSV/EMS, DOSS/Eynard linearity of the background in y, EO variational formula, BKMP for C^3 (Chen, Zhou) --
   and checked: genus 1 exactly, genus 2-3 numerically for the single class, genus 2 for tau = 1, 2) the base point q_0 = lam_0 - lam_1 is the
   spectral curve C_tau: x = w/(1+w)^(tau+1), y = -1/w (single class: x = z e^(-z), y = -1/z); the probes lam_j are the x^j-coefficients at
   the zero w = 0 of x, normalised by binom(j(tau+1)-1, j-1) (single class: j^j/j!). Conjecture M (vanishing) <=> the TR correlators have
   total x-degree >= 2h-2+n there; r = 1 <=> Res_{branch} y^N omega_{h,1} = 0 for N <= 2h-2 (N = 1 is the EO string identity = Theorem C).
3. (Theorem 3.1, proved given 2) For the triple (1,1,-2) (tau = 1 or -2) the vanishing direction of Conjecture M holds for ALL h >= 2 and ALL r:
   Xi_h(j_1..j_r) = 0 whenever 1 <= J <= 2h-3+r (j_i >= 1 by the theorem, j_i = 0 by the string/dilaton relations). Mechanism: global involution
   v -> 1/v of x = v/(1+v)^2, residue theorem, degree count deg Brk >= 2h-2+n against deg h_m = -m. For tau = -1/2 the same argument proves it
   for all-odd j_i.
4. (Proved, formal) single class: Xi_g(2) = 0 for all g >= 1 by cut-and-join on the line p = X delta_1 (gives the extra zero j = 2h at h = 1).
5. NOT proved: the single class for h >= 3 beyond j = 1, 2 (and the extra zero j = 2h for h >= 2); the CY family for tau not in {1, -2} (and
   tau = -1/2 with even parts); the non-vanishing direction. The obstruction is identified: no global involution, and the Lambert/framed
   involution maps the good zero of x to a point (essential singularity, resp. the Hurwitz-type zero v = 0) without vanishing.
   Most promising next steps: (a) Liu-Liu-Zhou framing cut-and-join (d/dtau of the Marino-Vafa generating function = CJ operator) at the base
   point with Theorem 3.1 as initial value at tau = 1; (b) the explicit x-y swap formula (dual curve (y, X) has no branch points; log-TR of
   Alexandrov-Bychkov-Dunin-Barkowski-Kazarian-Shadrin) evaluated at the pole of y; (c) for integer tau, global loop equations over the tau+1
   sheets of x = v^tau/(1+v)^(tau+1).
Scripts (copied to hodge_vertex_conifold/proof_M_scripts/, Perl, no Mathematica; run e.g. "perl trtau.pl 1 3 30" from that folder): tr11.pl (exact genus 1), trf.pl (Lambert curves, float, genus <= 3), trtau.pl (C_tau),
hookC*.pl / hookD.pl (Orlov-Scherbin hook formula tests; negative results recorded in section 2).
Additional checks at tau = 1 (trtau.pl, float): genus 3, Xi_3(1..4) = 6e-7 (= the uniform float error; Xi_3(5) = -7.318e-4 vs exact
-3/4096 = -7.324e-4); omega_{2,2} w-coefficients (1,1), (1,2) = 2e-10, 4e-10 (vanish, |mu| <= 3), (1,3), (2,2) non-zero, as Theorem 3.1 predicts.
