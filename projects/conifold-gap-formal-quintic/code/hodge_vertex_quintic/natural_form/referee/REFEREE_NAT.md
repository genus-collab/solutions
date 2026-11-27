# REFEREE_NAT -- independent adversarial check of the natural-variable factorisation

Date: 2026-11-10.  Folder: `hodge_vertex_quintic/natural_form/referee/`.  All arithmetic exact over F_p,
p in {2^521-1, 2^89-1, 2^61-1} (all = 1 mod 5; 2^521-1 is the referee's third prime), single process, every run < 1 minute.

**Claim refereed.**  With tau(u) = 1/2F1(1, 4/5; 9/5; u) = 1/Sum_n 4u^n/(5n+4), the Teleman form of the hatted R-matrix of the
formal quintic factorises as  R~(zhat; u) = G(u, zhat) M_0(tau(u) zhat),  G = 1 + Sum_{n>=1} u^n G_n(zhat),  G_n a polynomial in
zhat of degree 1, 3, 4, 5 (n = 1..4);  G_1 = zhat X_1 with no H (x) H component.

Status labels: PROVED (argument given), VERIFIED (exact computation; orders and primes stated), CONJECTURAL, REFUTED.

## 0. Verdict in one paragraph

The claim is **confirmed and survives every attempt to break it**.  A from-scratch re-derivation of the hatted rows (not using the
ring N, the s-series or any `fq_*` module) reproduces the project rows coefficient by coefficient and the factorisation holds
far beyond the orders originally tested: through u^11 and zhat^16 at three primes, through u^16 (tau fitted) and u^17 (tau
predicted) at 2^89-1 with zhat <= 22, 24.  The statement is not vacuous: tau is **unique** (PROVED) and is forced to be the
canonical coordinate; any other rescaling, and the placement of M_0 on the left, fail at the first order.  Three refinements /
corrections (referee's findings, all VERIFIED, all-order statements CONJECTURAL):
 1. deg_zhat G_n = floor((5n+2)/4) = 1, 3, 4, 5, 6, **8**, 9, 10, 11, 13, ... (not n+1: the pattern 1,3,4,5 does not continue
    as "n+1 for n >= 2"; G_5 has degree 6, G_6 has degree 8).
 2. Exact support rule ("rule W"): the component H^a (x) H^b of [u^n zhat^k] G is nonzero iff
    5n - 4k + 2 - w_a - w_b >= 0, w = (0, 4, 3, 2, 1).  Equivalently G - 1 is a formal power series in z whose coefficients are
    REGULAR at the conifold (power series in s = u^(1/5)) in the frame (1, uH, H^2, H^3, H^4).  The absence of H (x) H in G_1 is
    the first instance of this rule (H (x) H needs u^2).
 3. The leading (weight-0) part of G contains the classes H and 1 only through uH + z (left leg) and uH - z (right leg).

## 1. Item 1 -- independent recomputation  (verdict: CONFIRMED, VERIFIED)

**Method (ref_rows.py, ref_hat.py; nothing imported from the project).**  With D = q d/dq = -(1-u) d/du, L = 1/s, z = s^4 zhat,
the operator zD + zeta^i L acts on s^e Fh(zhat, u) as s^(e-1) A_e,  A_e = zeta^i + zhat (1-u) ((4/5) theta_zhat - e/5 - theta_u).
Writing the WKB function as psi_i = L Phi_i(zhat, u), the Picard-Fuchs equation of the formal quintic becomes

    (HPF)  Phi_5 - Phi_0 - (1-u) [ 3 zhat Phi_4 + (17/5) u zhat^2 Phi_3 + (9/5) u^2 zhat^3 Phi_2 + (274/625) u^3 zhat^4 Phi_1
                                   + (24/625) u^4 zhat^5 Phi_0 ] = 0,      Phi_j = A_{-j} Phi_{j-1},

whose zhat^(n+1) coefficient is zeta^(-i) (1-u)(4n - 5 theta_u) ph_n + (lower) = 0.  Hence the coefficients ph_n(u) of Phi are
POLYNOMIALS in u, unique unless 5 | n, where the coefficient of u^(4n/5) is free and is fixed by ph_n(1) = [z^n] h_{w_i}(z)
(QRR constant at q = 0).  Rows: Qh_0 = (zeta^i/C0) Phi, Qh_(a+1) = kh_a A_(e_a) Qh_a, e = (-1,-2,2,1,0),
kh = (1/C1, C0^2 C1^2, 1/C1, 1/C0), C0 = varpi_0, C1 = D t_c, with the conifold periods recomputed by a dual-number Frobenius
method (g = u + (7/10) u^2 + ..., f_0 = 1 + (2/625) u^3 + ...).  Teleman form Rt_ji = (1/5) Sum_{a+b=3} Qh_{b,0}^(j) Qh_a^(i).
Cost: 0.2 s for (zhat^6, u^5), 4 s for (zhat^24, u^17).

**Internal checks (all True at the three primes, all orders computed):** ph_n(1) = h_n automatically for every n not divisible
by 5 (tests the WKB and the QRR convention at 13 orders up to n = 16); vanishing resonance obstruction at n = 5, 10, 15, 20;
Rt(0; u) = 1 for all u; Rt(zhat)^T Rt(-zhat) = 1; covariance Qh_{a,n}^(i) = zeta^(i(a+1-n)) Qh_{a,n}^(0); circulant property of Rt.

**Agreement with the project rows** (`ref_run4.py --project`: fq_rows.Rows -> fq_conifold -> fq_sser.SRows -> fq3_msum.HatRowsU,
Nz = 6, prec_abs = 34, K = 6): all 1050 coefficients Qh[a][i][n][u^k] agree, mismatches 0, at 2^89-1 and 2^61-1.

**The claim at the original orders** (`ref_run1.py`, zhat <= 6, u <= 5, three primes): identity row
f = 1 - 3/20 zhat + 9/800 zhat^2 + 11571/80000 zhat^3 - 2415609/6400000 zhat^4 + 47735919/128000000 zhat^5 ...;
degrees of G_n = 0, 1, 3, 4, 5, (6 = truncation, untested there); tau fitted from the single top order zhat^6 =
(-4/9, -50/567, -4300/96957, -344425/12216582) = prediction; every Fourier component quoted in the claim is reproduced exactly:
u^1 zhat^1: (0,1),(1,0): -1/75, (2,4),(4,2): -8/225;  u^2 zhat^1: (0,1): -1/378, (2,4): 1/1134, (3,3): -127/15750;
u^2 zhat^2: (0,2): -17/23625, (1,1): 1/2250, (2,0): 73/23625, (3,4): -127/15750, (4,3): 127/15750;  u^2 zhat^3: (4,4): 127/15750
(normalisation Xhat = (1/25) Sum_ij zeta^(mi) X_ij zeta^(m'j)).

## 2. Item 2 -- higher orders  (verdict: CONFIRMED; degree pattern CORRECTED; VERIFIED)

| run | primes | zhat <= | u <= | result |
|---|---|---|---|---|
| ref_cross.py | 2^521-1, 2^89-1, 2^61-1 | 16 | 11 | degrees 1,3,4,5,6,8,9,10,11,13,14; 92 vanishing matrix coefficients (460 Z/5-allowed Fourier components) per prime; tau_1..tau_11 fitted from zhat^16 alone = prediction |
| ref_run1.py --NZ 22 --K 17 | 2^89-1 | 22 | 16 | tau_1..tau_16 fitted = prediction; degrees ..., 15, 16, 18, 19, 20 |
| ref_run2.py --NZ 24 --K 18 | 2^89-1 | 24 | 17 | degrees = floor((5n+2)/4) for all n <= 17; 215 vanishing matrix coefficients (1075 Fourier components); rule W exact |

* tau_5 = -32029400/1594263951, tau_6 = -497204834600/32441677138899 (= prediction; fitted values agree as residues at all primes).
* All 397 nonzero Fourier components of G (zhat <= 16, u <= 11) are RATIONAL: reconstructed at 2^521-1 (heights <= 40 digits,
  bound 78 digits) and their reductions agree with the residues at 2^89-1 and 2^61-1 (each prime has its own zeta).
  Stored in `G_fourier_NZ16_K12.json`.
* **Degree.**  deg G_n = floor((5n+2)/4) for 1 <= n <= 17 (VERIFIED).  The original "1, 3, 4, 5" is correct, but its natural
  reading "n+1" is wrong from n = 6 on (deg G_6 = 8, deg G_10 = 13, deg G_14 = 18).
* **Rule W (sharp support rule, VERIFIED in all 120 (component, k) pairs with bound < 18 at zhat <= 24; also at 2^61-1 and
  2^521-1 for zhat <= 16; all-order statement CONJECTURAL).**  Label the Fourier component (m, m') by the flat classes
  (a, b) = (-m-1, -m'-1) (item 4).  Then the minimal u-order at which the component H^a (x) H^b of zhat^k occurs in G is EXACTLY
      n_min(k; a, b) = max(1, ceil((4k - 2 + w_a + w_b)/5)),   w = (0, 4, 3, 2, 1).
  With the hatted weights c = (0,-1,3,2,1) the same formula is only a lower bound (25 non-sharp cases of 80): the class H needs
  one extra power of u per leg.  Since u^n zhat^k Hh^a (x) Hh^b = z^k s^(5n-4k+2-w_a-w_b) K^a (x) K^b with
  K = (1, uH, H^2, H^3, H^4) (hatted normalisations), rule W says: **G - 1 is a formal power series in z whose coefficients are
  power series in s = u^(1/5), i.e. regular at the conifold, in the frame K (x) K.**  The set {1 + regular K-bivectors} is a
  group (the weight is super-additive under composition: PROVED, elementary), so G^(-1) and log G obey the same rule
  (log G: 0 violations, 8 non-sharp cases of 67).  The polynomiality in zhat at fixed u-order and the degree formula are
  corollaries: the top coefficient sits at weight j := 5n - 4k + 2 - w_a - w_b = 0.
* **Leading coefficients.**  The top-degree coefficient of G_n is supported on weight 0 (Fourier-normalised values, labels (a,b)):
  n=1 (zhat^1): (0,2),(2,0): -8/225, (3,4),(4,3): -1/75;   n=2 (zhat^3): (0,0): 127/15750;
  n=3 (zhat^4): (0,4): -1958/4809375, (4,0): -8783/67331250;
  n=4 (zhat^5): (0,3): 283/14175000, (3,0): 43397/4039875000, (4,4): -12659/1009968750;
  n=5 (zhat^6): (0,2): 60247216/659004609375, (2,0): -150732703/1318009218750, (3,4): -34035397/439336406250, (4,3): 2401069/30830625000;
  n=6 (zhat^8): (0,0): 1416472237/3626797781250000.   (Full weight-0 and weight-5 tables: `logs/ref_run2b_521.log`.)
  No closed form was found for the numbers, but there is an exact structural pattern (**chain relation**, VERIFIED in all 13 + 13
  available cases, zhat <= 16, three primes): for weight-0 components, in flat components X^{ab} (of Hh^a (x) Hh^b),
      X^{0,b}[u^n zhat^(k+1)] = + X^{1,b}[u^n zhat^k],        X^{a,0}[u^n zhat^(k+1)] = - X^{a,1}[u^n zhat^k]      (k >= 1),
  i.e. the weight-0 part of G contains Hh^1 and Hh^0 only in the combinations Hh^1 + zhat Hh^0 = s^(-3)(uH + z) (left leg) and
  Hh^1 - zhat Hh^0 = s^(-3)(uH - z) (right leg).  Example: -(127/15750) [ zhat H(x)H-type, zhat^2 (H(x)1 - 1(x)H)-type,
  zhat^3 1(x)1-type ] at u^2, and again at u^6 (zhat^6, zhat^7, zhat^8) and u^10.  In particular the top coefficient of G_n for
  n = 2 mod 4 is a pure 1 (x) 1 term equal to minus the H (x) H coefficient two zhat-orders below.  The relation fails at
  weight 5 (0 of 12) and for the k = 0 -> 1 step in the flat normalisation.  It is the lowest-weight content of the G-equation
  of s.3 once rule W is granted (sketch only; the relation itself is VERIFIED).  REMARK (CONJECTURAL interpretation): uH - z is
  the Euler-minus-dilaton combination, which annihilates the homogeneous correlators t^(2-2g) of the frozen vertex theory; this
  is a plausible mechanism for the regularity of the edge graphs with M = M_0(tau zhat).

## 3. Item 3 -- attempts to break the claim  (verdict: NOT BROKEN; one variant REFUTED)

**(a) Is it trivially true for any reasonable tau?  NO.**
* Counting.  At u^n the claim is that [u^n zhat^k] G = 0 for all k > deg G_n.  Each such matrix coefficient has 5 Fourier
  components allowed by the Z/5 selection rule (m + m' = k mod 5), 2 or 3 of them independent after the symplectic symmetry.
  Original test (zhat <= 6): 5+3+2+1 = 11 matrix coefficients (55 components, >= 22 independent rational identities) against 4
  parameters tau_1..tau_4 - which are moreover PREDICTED, so against 0 parameters.  Referee test: 92 matrix coefficients (460
  components) per prime against 0 parameters; 215 (1075 components) at zhat <= 24.
* Perturbation (`ref_run3.py`, three primes): tau_1 -> tau_1 + 1 makes 67 of the 75 allowed components of [u^1] G above zhat^1
  nonzero (degree = truncation 16); tau_2 -> tau_2 + 1: 58 of 65; tau_5 -> tau_5 + 1: 45 of 50.
* Neighbouring "natural" choices all fail at the first order where they differ from tau: tau = 1; tau = (1-u)^(4/9) (same tau_1,
  fails at u^2); tau = 1/Sum 9u^n/(5n+9); tau = (u/t_c)^(4/5) (flat coordinate in place of the canonical one).
* **Uniqueness of tau: PROVED.**  If G = Rt M_0(c zhat)^(-1) and G' = Rt M_0(c' zhat)^(-1) were both polynomial at every u-order,
  then (G^(-1) = G(-zhat)^T is polynomial too) M_0(c zhat) M_0(c' zhat)^(-1) would be; at the first order where c, c' differ this
  says that zhat M_0' M_0^(-1) = (5/4) [ (Z - M_0 Z M_0^(-1))/zhat + B_0 ] is a polynomial (M-eq below), i.e. that A := M_0 Z M_0^(-1)
  is.  But (4/5) w^2 A' = [Z + w B_0, A] forces deg A <= 1 (ad B_0 has eigenvalues (mb - mb')/5, |.| <= 4/5) and then
  A_1 = B_0 + (diagonal), which is impossible because B_0 has four distinct nonzero Fourier-diagonal eigenvalues.  Hence at most
  one rescaling works; the fit from a single top-order entry (20 entries used, all consistent) finds it and it is the prediction.
* **tau is forced to be the canonical coordinate: PROVED** (given existence).  u_i/z = -(5/4) zeta^i/(tau(u) zhat) exactly, because
  Int_0^u t^(-1/5) (1-t)^(-1) dt = (5/4) u^(4/5) 2F1(1, 4/5; 9/5; u).  Equivalently (ref_run5.py):
  nabla f(c(u) zhat) = (4/5) w^2 f'(w) for nabla = zhat (1-u) ((4/5) theta_zhat - theta_u) iff (1-u)((4/5) F + theta_u F) = 4/5,
  F = 1/c, iff F = Sum 4u^n/(5n+4).  Sign remark: with D = q d/dq one has Int_{q_c}^{q} L dq/q = -(5/4) u^(4/5) 2F1, so
  tau zhat = (5/4) z / Int_{q}^{q_c} L dq/q (the claim's formula has the opposite orientation; immaterial).
* The three differential equations (referee's derivation; VERIFIED coefficient by coefficient for zhat <= 16, u <= 11, 3 primes):
      (R-eq)  nabla Rt = [Z, Rt] + zhat (1-u) B(u) Rt,      Z = diag(zeta^i),  B(u) = Psi^(-1) (E/5 + theta_u) Psi,  E = diag(-1,-2,2,1,0)
      (M-eq)  (4/5) w^2 dM_0/dw = [Z, M_0] + w B_0 M_0,      B_0 = B(0)
      (G-eq)  nabla G = [Z, G] + zhat ( (1-u) B(u) G - tau(u) G B_0 ).
  B(u) is antisymmetric, circulant, and diagonal in the Fourier (= flat) basis; B_0 f_m = (mb/5) f_m with mb in {0,1,2,-2,-1}
  the symmetric residue of m.  G is the formal solution of (G-eq) with G(0; u) = 1, unique up to the known z^(5k) diagonal
  constants (resonances 5n = 4k, fixed by the q = 0 normalisation; standard order-by-order recursion, sketch only).  What is NOT
  proved is that this solution terminates: the all-order existence of a polynomial G is CONJECTURAL (it is a uniform
  turning-point reduction statement).
* PROVED (conditional): if G_n is a polynomial of degree d and d > deg G_m for all m < n, then the Fourier component (m, m')
  of its top coefficient vanishes unless 5n - 4d = mb + mb' (read off (G-eq) at zhat^(d+1): ((4/5)d - n) g = [B_0, g]).  Hence
  d <= floor((5n+4)/4); with rule W (no index m = 3, i.e. no H leg, at weight 0) d <= floor((5n+2)/4), which is what is observed.

**(b) Is the Teleman form (left factor Psi(u)) essential?  NO for the polynomial property, YES for the normalisation.**
For any zhat-independent invertible A(u), A(u) G is polynomial of the same degrees (PROVED, trivial).  VERIFIED: the flat-colour
rows Qh(zhat; u) M_0(tau zhat)^(-1) = Psi(u) G have degrees 0,1,3,4,5,6,8,9,10,11,13,14.  The Teleman normalisation is what
makes G(0; u) = 1, G symplectic and G_1 = zhat X_1 without a zhat^0 term.  The statement is also independent of the conifold
frame: for (alpha, beta) = (1/3, -2/7) and (0, 5) the same degrees and rule W (0 violations, 0 non-sharp of 80) hold (2^89-1).

**(c) M_0 on the left: REFUTED.**  G' = M_0(tau zhat)^(-1) Rt is not polynomial at any order u^1..u^11 (degree = truncation 16),
and no other rescaling repairs it: at u^1 the 25 entries of the zhat^16 (and zhat^15) coefficient require 2 distinct values of c_1.
The factorisation is one-sided: R~ = (polynomial) . M_0(tau zhat), M_0 on the RIGHT (the side of the canonical exponentials).

**(d) Is G symplectic?  YES: PROVED and VERIFIED.**  Rt and M_0(tau zhat) are symplectic, so G(zhat)^T G(-zhat) = 1
automatically; checked mod (zhat^17, u^12).  Consequently G^(-1) = G(-zhat)^T is polynomial as well.

**(e) Is log G simpler than G?  Only marginally.**  Degrees of [u^n] log G: 1,3,4,5,6,7,9,10,11,13,14 (one lower at n = 6);
nonzero Fourier components per u-order 4,10,16,21,27,30,38,43,49,55,61 against 4,11,17,23,29,36,42,48,54,61,67 for G; log G is
in sp (zhat^k coefficient symmetric for k odd, antisymmetric for k even); rule W holds (not sharp in 8 of 67 cases).
The zhat-linear part of log G EQUALS that of G at every u-order (PROVED: G - 1 = O(zhat)).  Fourier components (m,m'):

    u^1 zhat^1: (0,1),(1,0): -1/75;  (2,4),(4,2): -8/225                                  [(3,3) = 0]
    u^2 zhat^1: (0,1),(1,0): -1/378; (2,4),(4,2): 1/1134;  (3,3): -127/15750
    u^2 zhat^2: (0,2): -1/525, (2,0): 1/525;  (3,4): -127/15750, (4,3): 127/15750           [(1,1) = 0 by antisymmetry]
    u^2 zhat^3: (4,4): 127/15750
    u^3 zhat^1: (0,1),(1,0): -43/32319; (2,4),(4,2): 716117/302990625;  (3,3): -167269/33665625
    u^3 zhat^2: (0,2): 1/71820, (2,0): -1/71820;  (3,4): -9128/4809375, (4,3): 9128/4809375
    u^3 zhat^3: (0,3),(3,0): 18629/134662500;  (1,2),(2,1): -12281/22443750;  (4,4): -13159/11221875
    u^3 zhat^4: (0,4): 18629/134662500, (4,0): -18629/134662500

The (3,3) = H (x) H component of the zhat-linear part (same for G and log G): u^1: 0;  u^2: -127/15750;  u^3: -167269/33665625;
u^4: -316834391/84837375000;  u^5: -8458079614/2767819359375;  u^6: -92268153113065624/35201472589951171875;
u^7: -47706385187645528392/20592861465121435546875 (all negative; rationals from 2^521-1, consistent with the other primes).
The structure of log G is the weight structure of s.2 (the -(127/15750) chain at u^2, the +-(18629/134662500) chain at u^3), not
a simplification beyond it.

## 4. Item 4 -- the dictionary and G_1  (verdict: CONFIRMED; PROVED + VERIFIED)

e~_j = Sum_a v_j^a Hh^a with v_j^a = (1/5) Qh_{3-a,0}^(j), and by covariance Qh_{b,0}^(j) = zeta^(j(b+1)) psi_b.  Hence (algebra)
Sum_j zeta^(mj) e~_j is proportional to H^(m-1), and a colour matrix X is the bivector Sum_ab X^{ab} Hh^a (x) Hh^b with
X^{ab} = psi_{3-a} psi_{3-b} Xhat(-a-1, -b-1):  (2,4) <-> 1 (x) H^2 (and (4,2) <-> H^2 (x) 1), (0,1) <-> H^4 (x) H^3, (3,3) <-> H (x) H.
VERIFIED numerically (`ref_run4.py`): psi(u=0) = (1,-1,-1,1,1); single flat component of every Fourier vector; the formula on 5
random matrices (125 components).  Flat components of G_1 computed DIRECTLY (no Fourier transform):
    [u^1 zhat^1] G = (8/225)(Hh^0 (x) Hh^2 + Hh^2 (x) Hh^0) - (1/75)(Hh^3 (x) Hh^4 + Hh^4 (x) Hh^3),   H (x) H component = 0 exactly;
    [u^2 zhat^1] G: H (x) H component -127/15750;   [u^3 zhat^1] G: -167269/33665625.
(The sign of the 1 (x) H^2 coefficient differs from the Fourier value -8/225 by psi_1 psi_3 = -1; a normalisation, not an error.)
G_1 has no H (x) H component at the three primes; by rule W this is structural (H (x) H at zhat^1 needs u^2).

## 5. Side verifications (not requested, cheap)

* **LG hand check extended: VERIFIED through zhat^10.**  The formal stationary-phase expansion of Int exp((x^4 + u/x)/z) dx at
  x_c (phase y^4 + 4/y after x = x_c y; its cubic Taylor coefficient at y = 1 vanishes) is 1 + (3/80) eps + (9/12800) eps^2
  - (11571/5120000) eps^3 + ..., and with eps = -4 zhat it equals the identity row f of the symbol for all orders 0..10
  (exact rationals; f from 2^521-1).  This is as expected: the leading part of the Picard-Fuchs equation at u -> 0, zhat fixed,
  is u (z d/du)^5 + 3 z (z d/du)^4 + 1 = 0 (D -> -d/du), which is the LG equation quoted in the task after u_LG = c u,
  c^4 = -1/4 (the same c that gives eps = -4 zhat), and f is its unique formal WKB solution.  The LG equation itself was taken
  from the task, not re-derived.  Only the identity row was compared; the rows a >= 1 were not checked separately.
* **The symbol is characterised by (M-eq)**: M_0 is the unique formal solution with M_0(0) = 1 of
  w^2 M' = (5/4) [Z, M] + (5/4) w B_0 M.  In Dubrovin's language U is proportional to diag(zeta^i) and V = -+(5/4) B_0 has
  spectrum mu(H^a) = +-(c_a - 1)/4 = +-(-1/4, -1/2, 1/2, 1/4, 0): the unit has |mu| = 1/4, i.e. d = 1/2, and the spectrum is
  {(b+1)/4 : b = -3..1}, the degrees of the forms x^b dx, b = a mod 5 (x^5 = u/4).  This is consistent with the conjectured
  dictionary H^a <-> x^a dx on M_{0;3,0} with d = 1/2 (up to the overall sign convention of mu).  The local exponents
  0,1,1,2,3 of the LG equation were NOT checked.

## 6. Remarks for the authors

1. State the result as: R~ = G . M_0((5/4) z / T), T the canonical coordinate from the conifold, with **G - 1 regular at the
   conifold in the frame (1, uH, H^2, H^3, H^4)** (rule W).  "G_n polynomial of degree 1,3,4,5" is the shadow of this at fixed
   u-order; the degree is floor((5n+2)/4), not n+1.
2. The absence of H (x) H in G_1 should not be presented as a separate observation: it is rule W.  H (x) H is present from u^2 on.
3. The original test left G_5 untested (degree 6 = truncation); it is now tested with 10 vanishing orders (deg G_5 = 6).
4. Nothing here proves the factorisation to all orders.  A proof would be a formal uniform reduction of (R-eq) to (M-eq) near
   the point where the five canonical coordinates coalesce; uniqueness, symplecticity and the conditional degree bound are
   proved above, existence is not.
5. Not refereed here: the graph-sum consequences (regularity of the edge graphs with M = M_0(tau zhat)), the dressing variants,
   and the identification of the rows a >= 1 of the symbol with LG integrals.

## 7. Files

Scripts (all new, referee's own; run as `py -I -B <script>`):
`ref_fp.py` (F_p helpers), `ref_rows.py` (periods, QRR constants, hatted WKB), `ref_hat.py` (rows, Teleman form, tau, G, fit),
`ref_run1.py` (items 1-2: fit and degrees), `ref_cross.py` (three primes, rational reconstruction), `ref_run2.py` (rule W),
`ref_run2b.py` (weight-0/5 parts of G and log G), `ref_run3.py` (item 3), `ref_run4.py` (item 4, LG check, project cross-check),
`ref_run5.py` (the three ODEs, chain relation).
Data: `G_fourier_NZ16_K12.json` (the 397 rational Fourier components of G, zhat <= 16, u <= 11; identity row of the symbol).
Logs (`logs/`): `ref_run1_{89,61,521}_NZ6_K6.log`, `ref_run1_89_NZ22_K17.log`, `ref_cross_NZ16_K12.log`,
`ref_run2_{89,61,521}.log`, `ref_run2_89_NZ24_K18.log`, `ref_run2_89_frame_1o3_m2o7.log`, `ref_run2_89_frame_0_5.log`,
`ref_run2b_521.log`, `ref_run3_{89,61,521}.log`, `ref_run4_89_project.log`, `ref_run4_61_project.log`, `ref_run4_521.log`,
`ref_run5_{89,61,521}.log`.
Project modules imported (read-only, only in `ref_run4.py --project`): fq_fp, fq_rows, fq_conifold, fq_sser, fq3_msum.
No existing project file was modified; the production runs were not touched.
