# The conifold symbol of the formal quintic is the R-matrix of the Landau-Ginzburg model x^4 + u/x (Hurwitz space H_{0;1,4} = M_{0;3,0}, primitive form dx) (2026-11-10)

Folder: hodge_vertex_quintic/natural_form/lg/. Everything here is NEW code and NEW files in this folder; no existing project file was changed (project modules are imported read-only; the symbol is re-derived in memory or read from the dumps lg_symbol_89.json / lg_symbol_61.json written here by lg_dump_symbol.py). Exact arithmetic only: rationals (fractions.Fraction) for the model, residues mod p = 2^89-1 and 2^61-1 for the quintic. AI-derived and AI-checked only; no human has read this. Status labels: PROVED (argument given here, unrefereed), VERIFIED (exact computation, order and primes stated), QUOTED (literature, what was actually read is stated), CONJECTURAL.

## 0. Summary

1. **Result (VERIFIED through zhat^18 for all five rows at both primes, through zhat^40 for the identity row; PROVED to all orders in s.4).** Let s = u^{1/5}, u = 1 - 5^5 q, zhat = z/s^4. The conifold symbol R^_a^{(i)}(zhat) (flat index a, colour i; CONIFOLD_SYMBOL_genus2.md s.1) is

        R^_a^{(i)}(zhat) = N_a zeta^{i(a+1)} J_{c_a}(-4 zeta^{-i} zhat),      c = (0, -1, 3, 2, 1),   N = (1, -1, -1, 1, 1),

   where J_k(eps) = 1 + O(eps) is the normalised formal stationary-phase series of Int exp(lambda/z_LG) x^k dx for lambda = x^4 + u_LG/x at a critical point x_c, eps = z_LG/x_c^4. The exact dictionary is **u_LG = 4u, z_LG = -4z**: then the critical points are x_c^5 = u, i.e. x_c = zeta^{-i} s (colour i), the critical values are lambda(x_c)/z_LG = -(5/4) zeta^i s^4/z = u_i/z (the canonical coordinate of the quintic measured from the conifold point, to leading order), and eps = -4 zeta^{-i} zhat. In particular R^_{a,0}^{(i)} = N_a zeta^{i(a+1)} and the identity row is f = J_0(-4 zhat).
2. **Correction of the conjectured dictionary (VERIFIED, negative result for the first guess).** The flat class H^a corresponds to x^{-a} dx (mod 5), NOT to x^a dx: 1 <-> dx, H <-> -dx/x, H^2 <-> -x^3 dx = -u x^{-2} dx, H^3 <-> x^2 dx, H^4 <-> x dx (exponents k = c_a; these are the five classes d_alpha lambda dx of the flat coordinates, s.5). The guess "row a = J_a" holds only for a = 0. The phase zeta^{i(a+1)} is right because colour i is the critical point zeta^{-i} s.
3. **Primitive form (VERIFIED).** The identity row equals J_k(-4 zhat) only for k = 0 among k = -3..6 (J_1 agrees through zhat^1 and fails at zhat^2; all others fail at zhat^1): the primitive form is dx, charge d = 1/2. The alternatives dx/x (k = -1, d = 1: the mirror of the weighted projective line P(1,4)) and x dx (k = 1) have the SAME five classes and the SAME Teleman matrix M_0 (s.5.3: Legendre-type transformations; proved directly for dx/x) but a different unit, hence a different CohFT; the quintic picks dx.
4. **Model connection and closed recursion (PROVED from the Gauss-Manin relations, VERIFIED on the stationary-phase series through zhat^10).** With g^a(zhat) := J_{c_a}(-4 zhat) (a in Z/5):

        (4/5) zhat (zhat d/dzhat) g^a = g^{a+1} - g^a + (zhat/5)(c_a - 1) g^a,        i.e.   g^{a+1}_n - g^a_n = ((4n - 3 - c_a)/5) g^a_{n-1},     Sum_a (4n + 1 - c_a) g^a_n = 0  (n >= 1).

   This gives the symbol, hence M_0, to ANY order without the quintic (lg_gm.symbol_rows; lg_symbol_rows_Q.json holds the rows through zhat^40). In the canonical frame: M_0 = Sum C_n zhat^n is the unique solution of C_0 = 1, (5/4)[C_{n+1}, Z] = (V - n) C_n with Z = diag(zeta^i), V_{ji} = (1/5) Sum_a zeta^{(i-j)(a+1)} mu_a, mu = (-1/4, -1/2, 1/2, 1/4, 0): the R-matrix of a conformal Frobenius manifold with Euler operator U = -(5/4) s^4 Z and spectrum mu, d = 1/2.
5. **Frobenius manifold (VERIFIED locally, QUOTED globally).** It is the Hurwitz space of rational functions f(p) = p^4 + a_2 p^2 + a_3 p + a_4 + a_5/(p - a_1) with primitive differential dp (Dubrovin's M_{0;3,0}; H_{0;1,4} in the notation of Carlet-van de Leur-Posthuma-Shadrin, arXiv:2211.12259, N = 5), on the line a = (0,0,0,0,4u), which is a rescaling of their "special point" p^4 + 1/p. Flat coordinates (a_1, a_2, a_3, a_4 - a_2^2/8, a_5) of degrees 1/4, 1/2, 3/4, 1, 5/4; unit d/da_4 <-> 1; the residue coordinate a_5 = 4u <-> H carries the conifold logarithm (c_555 = 1/a_5 exactly). By the quoted theorem its total descendant potential is a tau function of the extended rationally 4-constrained KP hierarchy.
6. **What fails / what is not shown.** The symbol is NOT of Hodge/Gamma type (known); nothing here says how the u-corrections G(u, zhat) or tau(u) of genus3_natural.py arise from the model (the formal quintic is not this Frobenius manifold away from u = 0); the identification of the all-genus frozen theory with the KP tau function is by the quoted theorem, not checked here beyond F_1 (s.5.4).

## 1. The stationary-phase engine (task item 1)

lg_sp.py. For Phi(y) = (a/2) y^2 + Sum_{m>=3} phi_m y^m and g(y) = Sum g_j y^j,

    Int exp(-Phi(y)/h) g(y) dy = sqrt(2 pi h/a) Sum_n A_n h^n      (formal),

computed by y = w t, w^2 = h/a, exp(-t^2/2) exp(-Sum_m (phi_m/a) w^{m-2} t^m) g(w t), Gaussian moments <t^{2j}> = (2j-1)!!. Exact rationals.

Validation (all PASS, lg_sp.py run as a script):
- Stirling: Gamma(n+1) = sqrt(2 pi n)(n/e)^n (1 + 1/(12n) + 1/(288 n^2) - 139/(51840 n^3) - 571/(2488320 n^4) + 163879/(209018880 n^5) + 5246819/(75246796800 n^6) - 534703531/(902961561600 n^7) - 4483131259/(86684309913600 n^8)): all eight coefficients reproduced.
- Airy saddle Phi = y^2 + y^3/3: A_k = (6k-1)!!/((2k)! 72^k) through k = 8.
- an exact Gaussian with polynomial amplitude, and reparametrisation invariance (Phi = (y + y^2)^2/2 with amplitude 1 + 2y gives exactly 1 through h^6).

The model: lambda = x^4 + u/x, x = x_c(1 + y), x_c^5 = u/4, eps = z/x_c^4, lambda/z = phi(y)/eps with phi = (1+y)^4 + 4/(1+y) = 5 + 10 y^2 + 0 y^3 + 5 y^4 - 4 y^5 + ..., so

    Int exp(lambda/z) x^k dx = x_c^{k+1} exp(5/eps) sqrt(2 pi (-eps)/20) J_k(eps).

J_k(eps) through eps^6 (lg_match.log; lg_sp prints eps^7, lg_gm uses eps^10):

| k | eps^1 | eps^2 | eps^3 | eps^4 | eps^5 | eps^6 |
|---|---|---|---|---|---|---|
| -2 | -9/80 | -63/12800 | 9009/5120000 | 2945943/1638400000 | 23432409/26214400000 | -9454287843/104857600000000 |
| -1 | -1/80 | -111/12800 | -13191/5120000 | -8841/1638400000 | 116525493/131072000000 | 93088145997/104857600000000 |
| 0 | 3/80 | 9/12800 | -11571/5120000 | -2415609/1638400000 | -47735919/131072000000 | 52989974037/104857600000000 |
| 1 | 3/80 | 21/2560 | 5229/5120000 | -1411641/1638400000 | -138080943/131072000000 | -57474780363/104857600000000 |
| 2 | -1/80 | 81/12800 | 16569/5120000 | 1504503/1638400000 | -47810763/131072000000 | -93810960243/104857600000000 |
| 3 | = J_{-2} identically | | | | | |
| 4 | -21/80 | -231/12800 | -14091/5120000 | 916839/1638400000 | 164837673/131072000000 | 102635329797/104857600000000 |
| 5 | -37/80 | -231/12800 | -32571/5120000 | -3252249/1638400000 | 8729721/131072000000 | 108222351237/104857600000000 |

The hand check of the task is confirmed and extended: J_0 = 1 + (3/80) eps + (9/12800) eps^2 - ..., and J_0(-4 zhat) = 1 - (3/20) zhat + (9/800) zhat^2 + (11571/80000) zhat^3 - (2415609/6400000) zhat^4 + (47735919/128000000) zhat^5 + (52989974037/25600000000) zhat^6 = f through zhat^6 (all seven stored coefficients).

## 2. The exact relation between the symbol and the J_k (task item 2)

lg_match.py (log lg_match.log; the first run with the naive expectation N_a = 1 is kept as lg_match_run1.log), both primes.

**(A) Agnostic fit.** For every flat index a and every k in -3..6: N_a := R^_{a,0}^{(0)}, c from the zhat^1 coefficient, then the orders zhat^2..zhat^6 are tests (5 per row). Result, identical at both primes:

| a | matching k | N_a | c (eps = c zhat) | class |
|---|---|---|---|---|
| 0 | 0 | 1 | -4 | dx (unit) |
| 1 | -1 | -1 | -4 | dx/x |
| 2 | 3 and -2 (J_3 = J_{-2}) | -1 | -4 | x^3 dx = u x^{-2} dx |
| 3 | 2 | 1 | -4 | x^2 dx |
| 4 | 1 | 1 | -4 | x dx |

No other (a, k) passes. So T is a signed permutation in the window k in [-1, 3]: R^_a^{(0)}(zhat) = N_a J_{c_a}(-4 zhat). In the window k in [0, 4] of the task (using J_{k+5} = J_k - ((k+2) eps/4) J_{k+1}, s.3): R^_0 = J_0, R^_1 = -(J_4 - zhat J_0), R^_2 = -J_3, R^_3 = J_2, R^_4 = J_1 (arguments -4 zhat): monomial except for a = 1, where the polynomial correction -zhat J_0 appears.

**(B) All colours.** All 175 stored residues satisfy R^_{a,n}^{(i)} = N_a zeta^{i(a+1-n)} (-4)^n j_{c_a,n} (0 mismatches, both primes). Negative: "row a = N_a J_a" (H^a <-> x^a dx) holds only for a = 0. Unit test: f = J_k(-4 zhat) only for k = 0; first failing order for k = -3..6: 1, 1, 1, none, 2, 1, 1, 1, 1, 1.

**(C) Controls.** The identity row against N J_k(c zhat), k = -6..6, for the other rank-5 potentials x^P + u/x^Q, P + Q = 5: (3,2) and (2,3): no match; (1,4): only k = -2, c = -1, which is the same model under x -> 1/x (x + u/x^4, x^{-2} dx = -d(1/x)); (4,1): only k = 0, c = -4.

**(D) Teleman form.** With N_a N_{3-a} = 1 and eta^{ab} = (1/5) delta_{a+b = 3 mod 5} (checked on the dump):

    R~_{ji}(zhat) = (1/5) Sum_a zeta^{(i-j)(a+1)} J_{c_a}(-4 zeta^{-i} zhat),     C_n[j][i] = (1/5) Sum_a zeta^{j(4-a) + i(a+1-n)} g^a_n:

all 175 entries of C_0..C_6 rebuilt from the dumped symbol agree at both primes.

**Higher orders (beyond the stored zhat^6).** lg_rows_highorder.py re-runs the project pipeline (Rows -> Conifold -> SRows -> Symbol) in memory with larger Nz: all 25 (Nz+1) residues agree with the recursion of s.3 for Nz = 12 (325 residues) and Nz = 18 (475 residues) at both primes, with symplecticity and covariance of the symbol holding at that order. lg_wkb_check.py takes the identity row straight from the WKB recursion of the Picard-Fuchs operator (fq_rows.WKB): deg_L phi_n = 4n and [L^{4n}] phi_n = zeta^{-in} g^0_n for all n <= 40, colours 0 and 1, both primes (a wrong row, g^4 in place of g^0, is rejected).

**Independence of the conifold frame.** With (alpha, beta) = (1/3, -2/7) instead of (0, 0) the symbol is unchanged through zhat^8 (225 residues, both primes): this settles the caveat of CONIFOLD_SYMBOL_genus2.md s.5 ("plausible but not tested") to that order; s.4 proves it.

**Normalisation at zhat = 0.** R^_{a,0}^{(i)} = N_a zeta^{i(a+1)}; equivalently Q_a^{(i)} = N_a x_c^{c_a - 1}(1 + O(zhat)) with x_c = zeta^{-i} s: the leading term of d_alpha lambda(x_c)/sqrt(lambda''(x_c)) times sqrt(20) (s.5.2).

The rows g^a_n = [zhat^n] J_{c_a}(-4 zhat) (the symbol is N_a zeta^{i(a+1-n)} g^a_n):

| n | a = 0 (k = 0) | a = 1 (k = -1) | a = 2 (k = 3) | a = 3 (k = 2) | a = 4 (k = 1) |
|---|---|---|---|---|---|
| 1 | -3/20 | 1/20 | 9/20 | 1/20 | -3/20 |
| 2 | 9/800 | -111/800 | -63/800 | 81/800 | 21/160 |
| 3 | 11571/80000 | 13191/80000 | -9009/80000 | -16569/80000 | -5229/80000 |
| 4 | -2415609/6400000 | -8841/6400000 | 2945943/6400000 | 1504503/6400000 | -1411641/6400000 |
| 5 | 47735919/128000000 | -116525493/128000000 | -23432409/25600000 | 47810763/128000000 | 138080943/128000000 |
| 6 | 52989974037/25600000000 | 93088145997/25600000000 | -9454287843/25600000000 | -93810960243/25600000000 | -57474780363/25600000000 |
| 7 | -7749624763881/512000000000 | -2450627360181/512000000000 | 7230539823507/512000000000 | 6398562493323/512000000000 | -2232045849033/512000000000 |
| 8 | 720285444605727/16384000000000 | -3590224557852933/81920000000000 | -5942826823626693/81920000000000 | 72982309531131/81920000000000 | 5601340303762203/81920000000000 |

(orders 7, 8 are predictions of the model relative to the stored symbol; they are confirmed mod both primes by lg_rows_highorder.py. lg_symbol_rows_Q.json has n <= 40.)

## 3. Gauss-Manin system, ODEs, model connection, closed recursion (task item 3)

Use the equivalent W-form, which removes all constants: with x = xi, u_LG = 4u, z_LG = -4z,  lambda/z_LG = -W/z,

    W(xi) = xi^4/4 + u/xi,     W' = (xi^5 - u)/xi^2,     I_k(u, z) := Int exp(-W/z) xi^k dxi,
    I_k ~ exp(-5/(4 zhat)) (2 pi z/5)^{1/2} s^{k-1} g_k(zhat),   g_k(zhat) := J_k(-4 zhat)     (colour 0: xi_c = s, W(xi_c) = (5/4) s^4, W''(xi_c) = 5 s^2).

**Relations (PROVED for the integrals by differentiation under the integral and by d(exp(-W/z) xi^m) = exp(-W/z)(-(xi^5 - u) xi^{m-2}/z + m xi^{m-1}) dxi; VERIFIED on the formal series for k = -4..8 through zhat^10, lg_gm.py):**

    (R1)  z d/du I_k = -I_{k-1}                 <=>   g_{k-1} = g_k - (zhat/5)(k - 1 - 4 theta) g_k,       theta = zhat d/dzhat,
    (R2)  I_{m+3} = u I_{m-2} + m z I_{m-1}      <=>   g_{k+5} = g_k + (k+2) zhat g_{k+1}     (k = m - 2; in particular I_3 = u I_{-2}, g_3 = g_{-2}),
    (R3)  [u (z d/du)^5 - (k-3) z (z d/du)^4 + 1] I_k = 0
                                                <=>   zhat^5 T(T-1)(T-2)(T-3)(T-k-1) g_k + g_k = 0,   T = -1/zhat + (k-1)/5 - (4/5) theta.

In the task's variables (lambda = x^4 + u_LG/x, D = z_LG d/du_LG) (R3) reads u_LG D^5 I_k - (k-3) z_LG D^4 I_k - 4 I_k = 0; for k = 0 this is the hand-derived u D^5 I + 3 z D^4 I - 4 I = 0. Local exponents at u = 0: 0, 1, 2, 3, k+1 (for k = 0: 0, 1, 1, 2, 3, logarithmic; for k = -1: 0, 0, 1, 2, 3). Note the factorisation u (z d/du)^5 + 3 z (z d/du)^4 = (z d/du)^3 u (z d/du)^2.

**The model connection in u (the 5 x 5 system).** Put Psi_a := N_a I_{c_a}, i.e. (Psi_0, ..., Psi_4) = (I_0, -I_{-1}, -I_3, I_2, I_1). Then (R1), (R2) give the closed cyclic system

    z d/du Psi_0 = Psi_1,    z d/du Psi_1 = -Psi_2/u,    z d/du Psi_2 = Psi_3,    z d/du Psi_3 = -Psi_4,    z d/du Psi_4 = -Psi_0,

i.e. z d/du Psi = A(u) Psi with A = E_{01} - u^{-1} E_{12} + E_{23} - E_{34} - E_{40}: regular singular at u = 0 (the conifold point, with the 1/u in the H o H = u^{-1} H^2 entry), irregular at infinity. Eliminating gives (R3) for Psi_0: [(z d/du)^3 u (z d/du)^2 + 1] Psi_0 = 0.

**The model connection in zhat.** By homogeneity (xi -> c xi, u -> c^5 u, z -> c^4 z) everything depends on zhat only. For g^a := g_{c_a}, a in Z/5 (so g^{a+1} = g_{c_a - 1}, using g_{-2} = g_3 at a = 1):

    (MC)   (4/5) zhat theta g^a = g^{a+1} - g^a + (zhat/5)(c_a - 1) g^a.

For the full matrix of rows Rhat_{ai}(zhat) := R^_a^{(i)}(zhat) and Y := Rhat exp(Uhat/zhat), Uhat = -(5/4) diag(zeta^i):

    zhat d/dzhat Y = ( (5/4) P / zhat + mu ) Y,      (P Y)_a = sigma_a Y_{a+1},  sigma = N_{a+1}/N_a = (-1, 1, -1, 1, 1),     mu = diag((c_a - 1)/4) = diag(-1/4, -1/2, 1/2, 1/4, 0).

This is Dubrovin's z-connection of a conformal Frobenius manifold: (5/4) P = -(E o) in the flat frame (eigenvalues -Uhat_i), mu_a = q_a - d/2 with q_a = c_a/4 and d = 1/2 (mu_0 = -1/4 = -d/2 for the unit; mu_a + mu_{3-a} = 0).

**Closed recursion (flat frame).** From (MC):  g^a_0 = 1;  for n >= 1, d^0 := 0, d^{a+1} := d^a + ((4n - 3 - c_a)/5) g^a_{n-1} (a = 0..3), t := -Sum_a (4n + 1 - c_a) d^a/(20 n), g^a_n = t + d^a. (The cycle closes: d^5 = 0 by the level-(n-1) constraint; asserted in the code.) lg_gm.symbol_rows(N): VERIFIED against the stationary-phase series for all five rows through zhat^10 (50 non-trivial coefficients) and against the quintic as in s.2. Its output is symplectic through zhat^40 (exact rationals): for every r mod 5, Sum_{a, m : a+1-m = r} (-1)^{n-m} g^a_m g^{3-a}_{n-m} = delta_{n0}.

**Scalar recursion for the identity row f = g^0** (from (R3), k = 0; P_4(m) = -4m for every k, so the formal solution with f_0 = 1 is unique):

    4 n f_n = -(1/5)(32 m^2 + 32 m + 3)|_{m = n-1} f_{n-1} - (4/25)(32 m^3 + 96 m^2 + 73 m + 12)|_{m = n-2} f_{n-2}
              - (4/625)(320 m^4 + 1920 m^3 + 3700 m^2 + 2580 m + 459)|_{m = n-3} f_{n-3} - (1/3125)(4m+1)(4m+6)^2(4m+11)(4m+16)|_{m = n-4} f_{n-4}

(the polynomials P_3, ..., P_0 printed in lg_gm.log; e.g. f_1 = -3/20, f_2 = (-(67/5)(-3/20) - 48/25)/8 = 9/800).

**Closed recursion (canonical frame).** With Z = diag(zeta^i) and V := Psi^{-1} mu Psi, V_{ji} = (1/5) Sum_a zeta^{(i-j)(a+1)} mu_a = (1/20)(zeta^{4d} - zeta^d) + (1/10)(zeta^{3d} - zeta^{2d}), d = i - j (antisymmetric):

    C_0 = 1,     (5/4) [C_{n+1}, Z] = (V - n) C_n     (off-diagonal part of C_{n+1}),     (n+1) (C_{n+1})_{ii} = Sum_{k != i} V_{ik} (C_{n+1})_{ki}   (diagonal part, from the next order).

VERIFIED mod both primes for the dumped symbol (n <= 5) and for the recursion output (n <= 30), lg_frob.py (e). This is equation (1.18) of arXiv:2211.12259 ([Psi^{-1} R_{m+1} Psi, U] = (m + mu) Psi^{-1} R_m Psi) up to the conventions for rows versus columns.

## 4. Why it is true: the scaling limit of the Picard-Fuchs operator (PROVED, unrefereed)

Notation of fq_rows.py: D = q d/dq, L = (1 - 5^5 q)^{-1/5} = 1/s, l = zeta^i L, M = zD, and the PF operator conjugated by exp(u_i/z) (D u_i = l) is

    Pcal = (M + l)^5 - 1 - q [3125 (M+l)^5 + 9375 z (M+l)^4 + 10625 z^2 (M+l)^3 + 5625 z^3 (M+l)^2 + 1370 z^4 (M+l) + 120 z^5],   3125 q = 1 - L^{-5}.

**Grading.** Give z^m L^e D^r the degree deg* := e + 5r - 4m (s has weight 1, z weight 4, D L^k = (k/5)(L^{k+5} - L^k) raises the L-degree by at most 5). The filtration by deg* <= d is multiplicative under composition. M + l has deg* <= 1; 1 - 3125 q = L^{-5} has deg* = -5; so Pcal = L^{-5}(M+l)^5 - 1 - (1 - L^{-5})[3 z (M+l)^4 + (17/5) z^2 (M+l)^3 + ...] has deg* <= 0, the terms z^j (M+l)^{5-j} having deg* <= 5 - 5j <= -5 for j >= 2, and

    gr_0(Pcal) = u (z Dbar + l)^5 - 1 - 3 z (z Dbar + l)^4,      Dbar L^k := (k/5) L^{k+5},  i.e.  Dbar = -d/du  (u = L^{-5} = s^5).

**Lemma 1 (degree bound).** psi_n = L phi_n with phi_n a Laurent polynomial in L of degree <= 4n. Proof: induction. z^m Pcal_m has deg* <= 0, so Pcal_m raises the L-degree by at most 4m; with deg psi_j <= 4j + 1 for j < n, r_n = -(zeta^i/5) Sum_{m=2}^{n+1} Pcal_m psi_{n+1-m} has degree <= 4n + 5, and D phi_n = r_n with [L^{k+5}] D phi = (k/5)[L^k] phi (k > 0 top) gives deg phi_n <= 4n. (The integration constant of the WKB recursion, fixed at L = 1, sits in degree 0 < 4n and does not matter.) VERIFIED with equality for n <= 40 (lg_wkb_check.py).

**Lemma 2 (top parts).** The top parts psi_n^top := ([L^{4n+1}] psi_n) L^{4n+1} satisfy gr_0(Pcal) Sum_n z^n psi_n^top = 0 (the deg* = 1 component of Pcal psi = 0 involves top parts only). Undoing the conjugation with uhat_i := -(5/4) zeta^i L^{-4} (Dbar uhat_i = l), the series Psi_0 := exp(uhat_i/z) Sum_n z^n psi_n^top = exp(uhat_i/z) s^{-1} f^{(i)}(zhat), f^{(i)}_n = [L^{4n}] phi_n, is a formal solution of

    [u (z d/du)^5 + 3 z (z d/du)^4 + 1] Psi_0 = 0,      equivalently  [(z d/du)^3 u (z d/du)^2 + 1] Psi_0 = 0,

which is (R3) with k = 0, the Gauss-Manin equation of I_0 = Int exp(-W/z) dxi, W = xi^4/4 + u/xi.

**Lemma 3 (uniqueness).** (R3), k = 0, has exactly one formal solution of the form exp(-5/(4 zhat)) s^{-1} (1 + Sum_{n>=1} f_n zhat^n): its recursion is 4 n f_n = (lower terms) (s.3; P_4(m) = -4m). The other colours follow by xi -> zeta xi. Hence f^{(i)}(zhat) = g^0(zeta^{-i} zhat) to all orders, with g^0 defined by the recursion.

**Rows a >= 1.** fq_rows: X_0 = psi_i/C0, X_{a+1} = (C1 kappa_a)^{-1} (zD + l) X_a with (C1 kappa_a)^{-1} = (1/C1, C0^2 C1^2 L^{-5}, 1/C1, 1/C0) for a = 0..3, P_a = (C0/L) X_a, Q_a = P_a/Delta_i^{1/2}, Delta_i^{1/2} = zeta^{4i} C0/L. In the conifold frame C0 -> varpi_0 = 1 + O(u) and C1 -> D t_c = -1 + O(u) (checked for (alpha, beta) = (0,0) and (1/3,-2/7): lg_cf_leading.log), L^{-5} = u, D = -(1-u) d/du. All corrections have relative weight 5 (a factor u), so the leading (hatted, u^0) parts are obtained with the leading operators: with Psi_a := exp(uhat_i/z) X_a^{top},

    Psi_1 = z d/du Psi_0,    Psi_2 = -u z d/du Psi_1,    Psi_3 = z d/du Psi_2,    Psi_4 = -z d/du Psi_3,

and the closure -z d/du Psi_4 = Psi_0 is exactly Lemma 2 ((z d/du)^3 u (z d/du)^2 Psi_0 = -Psi_0). This is the cyclic system of s.3, whose solution is Psi_a = N_a I_{c_a} with N = (1, -1, -1, 1, 1): the signs found numerically are the signs of this chain (H acts as -d/du, and the H o H structure constant kappa_1 carries u). Therefore:

**Theorem A (PROVED; all orders; any conifold frame (alpha, beta)).** The symbol exists to all orders in zhat and R^_a^{(i)}(zhat) = N_a zeta^{i(a+1)} g^a(zeta^{-i} zhat), where the g^a are the unique solution of (MC) with g^a(0) = 1 (equivalently of the cyclic u-system). Inputs: the definitions above, Lemmas 1-3, and the passage from the u-system to (MC) by homogeneity (a substitution; checked as (R1)-(R2) of lg_gm.py).

**Theorem B (PROVED modulo one standard fact; independently VERIFIED through zhat^10 in exact rationals).** g^a = J_{c_a}(-4 zhat), the stationary-phase data of lambda = x^4 + 4u/x, z_LG = -4z, form x^{c_a} dx. The standard fact: the formal stationary-phase series of an oscillating integral satisfies the differential relations (R1), (R2) that hold for the integral (compatibility of formal Gaussian integration with differentiation under the integral sign and integration by parts). It is not re-proved here; the 50 coefficients checked by lg_gm.py (MC) do not depend on it.

Remark. The rows are u-series because all corrections come with u = s^5; the symbol is the gr_0 part. The first corrections come from (i) the factor (1 - u) in 3125 q and in D, (ii) the term (17/5) z^2 (M+l)^3 of Pcal (deg* = -5), (iii) the O(u) parts of C0, C1. They are not part of the model.

## 5. The Frobenius structure (task item 4)

### 5.1 The manifold, flat coordinates, charge (VERIFIED by lg_frob.py at random rational points; QUOTED formulas agree)

M = {f(p) = p^4 + a_2 p^2 + a_3 p + a_4 + a_5/(p - a_1)}: rational functions on P^1 with a pole of order 4 at infinity and a simple pole; Dubrovin's Hurwitz space M_{0;3,0} (genus 0, n_0 = 3, n_1 = 0), rank 5 = number of critical points. With the primitive differential dp:

- metric (d', d'') = -(Res_{p=inf} + Res_{p=a_1}) d'f d''f/f' dp = sum over the critical points of d'f d''f/f''. In t = (t1..t5) = (a_1, a_2, a_3, a_4 - a_2^2/8, a_5) it is constant: eta_15 = 1, eta_24 = eta_33 = 1/4 (4 random rational points; control: with a_4 in place of a_4 - a_2^2/8 the (2,2) entry is -a_2/8). These are (1.7)-(1.8) of arXiv:2211.12259 for N = 5 (t^alpha = -(4/(alpha-1)) Res_inf f^{(alpha-1)/4} dp gives a_2, a_3, a_4 - a_2^2/8).
- unit e = d/dt4 (f -> f + c; c_{4 alpha beta} = eta_{alpha beta} checked): the class 1.dp, i.e. the quintic's H^0. Under dx/x the unit would be d/dt5.
- Euler field E = Sum (alpha/4) t_alpha d_alpha (homogeneity p -> c p, f -> c^4 f): deg t = (1/4, 1/2, 3/4, 1, 5/4); the pairs of eta have total degree 3/2 = 2 - d, so **d = 1/2** (= 1 - 2/(n_0 + 1), = (N-3)/(N-1) of the quoted (1.10)). mu_alpha = (2-d)/2 - deg t_alpha = (3 - alpha)/4 = (1/2, 1/4, 0, -1/4, -1/2).
- **the conifold logarithm**: c_555 = 1/a_5 exactly at every point (4 random points), so the prepotential contains (1/2) t5^2 log t5: the logarithm sits in the RESIDUE coordinate t5 = a_5, of degree 5/4 (so q = 1 - 5/4 = -1/4 < 0: t5 is the one coordinate below the unit in the spectrum).

### 5.2 The line and the dictionary (VERIFIED)

On the line a = (0,0,0,0,r): critical points p^5 = r/4, critical values 5 p_c^4, E = (5/4) r d_5, and multiplication by d_5 is cyclic: d_4 -> d_5 -> (1/r) d_1, d_1 -> 4 d_2, d_2 -> d_3, d_3 -> d_4, (d_5 o)^5 = 4/r (lg_frob.py (c)). The tangent classes are d_alpha f dp = (r p^{-2}, p^2, p, 1, p^{-1}) dp for alpha = 1..5, with r p^{-2} dp = 4 p^3 dp modulo exact forms. Dictionary (r = 4u, z_LG = -4z):

| quintic flat class | 1 | H | H^2 | H^3 | H^4 |
|---|---|---|---|---|---|
| Hurwitz flat field (alpha = a - 1 mod 5) | sqrt(20) d_4 | -sqrt(20) d_5 | -(sqrt(20)/4) d_1 | sqrt(20) d_2 | sqrt(20) d_3 |
| class x^k dx, k = c_a | dx | dx/x | x^3 dx = u x^{-2} dx | x^2 dx | x dx |
| mu_a = (c_a - 1)/4 | -1/4 | -1/2 | 1/2 | 1/4 | 0 |
| weight of t (= 4 - c_a) / 4 | 1 | 5/4 | 1/4 | 1/2 | 3/4 |

This map is an isometry eta_quintic = 5 delta_{a+b=3} <-> eta_Hurwitz (e.g. eta_q(H, H^2) = 20 (1/4) eta_15 = 5). Unit-preserving version: 1 -> d_4 etc. with eta_q = 20 eta_Hurwitz. The exponents c = (0,-1,3,2,1) of the hatted theory are 4 - (weight of the flat coordinate): the "hat" is the Euler grading of this Frobenius manifold. H is (minus) the derivative along the residue coordinate t5 = 4u: the conifold coordinate t_c = u + ... IS the residue of the superpotential (up to 4), and the only direction with a logarithm.

### 5.3 Primitive forms and Legendre-type transformations (PROVED for dx versus dx/x; QUOTED in general)

The primitive differentials of the Hurwitz space are the d_kappa f dp, kappa = 1..5; on the line these are the classes x^k dx with k = -2, 2, 1, 0, -1. Direct argument for dx versus dx/x: for dx/x the natural parametrisation is f = x^4 + b_3 x^3 + b_2 x^2 + b_1 x + b_0 + q/x (x = p - a_1 fixed), its tangent functions on the line are x^3, x^2, x, 1, 1/x of five distinct weights, so the flat frame is, up to scalars, the set of classes x^k dx/x, i.e. again x^k dx with k in {2, 1, 0, -1, -2 = 3}; the normalised canonical frame is the same because the Hessian is taken in the coordinate adapted to the form (d^2 f/d(log x)^2 = x^2 f''), which compensates the change of form. Hence Psi R (the matrix of normalised oscillating integrals) and the Teleman matrix M_0 are common to the two structures. In general (QUOTED from memory: Dubrovin, App. B) the structures for different primitive differentials are related by Legendre-type transformations S_kappa, which preserve the canonical coordinates and the rotation coefficients, hence (U, V) and the R-matrix in the normalised canonical frame. They differ in the unit (the class of the primitive form itself), i.e. in the TFT Delta_i and in d = -2 mu_unit: u x^{-2} dx: d = -1; x^2 dx: d = -1/2; x dx: d = 0; dx: d = 1/2; dx/x: d = 1. The frozen conifold theory has unit 1 <-> J_0 and Delta^_i^{1/2} = zeta^{4i} = 1/x_c up to scale: primitive form dx (s.2 (B)). The d = 1 member (dx/x) is the Landau-Ginzburg mirror of the weighted projective line P(1,4) (x^4 + Q/x, form dx/x; Milanov-Tseng's spaces of Laurent polynomials; from memory, exact form not re-checked): M_0 is also the R-matrix of its small quantum orbifold cohomology in the normalised canonical frame, but that CohFT is a different one.

### 5.4 Genus one and the gap (VERIFIED at genus 1 from the model; consistency at genus 2, 3)

Givental's formula dF_1 = (1/48) Sum_i dlog Delta_i + (1/2) Sum_i (R_1)_{ii} du_i with Delta_i ~ s^2, (R_1)_{ii} = zeta^{-i} C_1[0][0]/s^4, C_1[0][0] = (1/5) Sum_a g^a_1 = 1/20, du_i = -zeta^i s^4 du/u gives dF_1 = (1/24 - 1/8) du/u = -(1/12) du/u, the conifold value quoted in OU_EXTENSION_genus3.md s.6d. Scaling: the frozen quintic theory is the Hurwitz CohFT with eta -> 20 eta and z -> -4z, so F_g^quintic(u) = 20^{1-g} (-4)^{3g-3} F_g^Hurwitz(a_5 = 4u); if F_g^Hurwitz = c_g a_5^{2-2g} then F_g^quintic = (-1)^{g-1} 5^{1-g} c_g u^{2-2g}. The project values 1/1200 and 1/25200 correspond to c_2 = -1/240, c_3 = 1/1008, i.e. c_g = B_{2g}/(2g(2g-2)): the c = 1 / Gaussian-matrix-model free energy that goes with the prepotential term (1/2) a_5^2 log a_5. (That the Hurwitz theory has exactly these F_g on the line is here an inference from the project's frozen sums plus Theorem A, not an independent computation.)

### 5.5 Literature (QUOTED; what was read is stated)

- B. Dubrovin, Geometry of 2D topological field theories (1996), Lecture 5: Frobenius structures on Hurwitz spaces M_{g;n_0,...,n_m}, flat coordinates, primitive differentials, d = 1 - 2/(n_0+1) for the differential dp; Appendix B: Legendre-type transformations. (From memory; not re-read.)
- G. Carlet, J. van de Leur, H. Posthuma, S. Shadrin, "Enumeration of hypermaps and Hirota equations for extended rationally constrained KP", arXiv:2211.12259 (read: abstract, introduction, s.1 pp.1-5). Same manifold H_{0;1,N-1}, f(p) = p^{N-1} + a_2 p^{N-3} + ... + a_{N-1} + a_N/(p - a_1); e = d/da_{N-1}; flat coordinates (1.7), eta (1.8), d = (N-3)/(N-1), mu (1.11); special point f = p^{N-1} + 1/p; R-matrix defined by (1.18). Main result as stated there: the total descendant potential of H_{0;1,N-1} is a tau function of the extended rationally (N-1)-constrained KP hierarchy (a version of a conjecture of Liu-Zhang-Zhou), and its restriction to one set of times enumerates rooted hypermaps. For N = 5 this is the manifold of the conifold symbol; the special point is a_5 = 1, our line is a_5 = 4u.
- S.-Q. Liu, Y. Zhang, X. Zhou, "Central invariants of the constrained KP hierarchies", arXiv:1502.01774 (abstract only): the central invariants are 1/24, the constrained KP hierarchies are topological deformations of their dispersionless limits.
- S. Aoyama, Y. Kodama, "Topological Landau-Ginzburg theory with a rational potential and the dispersionless KP hierarchy", Comm. Math. Phys. 182 (1996) 185-220, hep-th/9505122 (abstract only): topological LG theory with rational superpotential, residue formulas, descendants, relation to Dubrovin's flat solutions.
- arXiv:1401.5725, "Rational reductions of the 2D-Toda hierarchy and mirror symmetry", JEMS 19 (2017) 835-880 (abstract only): rational reductions of 2D Toda; their dispersionless limits are NON-conformal Frobenius manifolds with flat identity, mirror to equivariant local orbifold lines. This is the Toda (logarithmic) relative of the present structure, not the same manifold: the structure found here is the conformal KP-type one (primitive form dp). The Ablowitz-Ladik / resolved-conifold case sits in that Toda family.
- T. Milanov, H.-H. Tseng, "The spaces of Laurent polynomials, Gromov-Witten theory of P^1-orbifolds, and integrable hierarchies", J. reine angew. Math. (2008), arXiv:math/0607012 (search summary only): Laurent polynomials with form dx/x, orbifold lines C_{k,m}, extended bigraded Toda.

## 6. Consequences for the programme

1. M_0 is no longer "tautological data of the quintic": it is the R-matrix of the conformal Frobenius manifold H_{0;1,4} (d = 1/2) on its residue line, computable to any order from (U, mu) alone. The frozen vertex theory M_0.omega^ is the Givental CohFT of that manifold (metric x 20, z -> -4z), whose all-genus potential is, by the quoted theorem, a tau function of the extended rationally 4-constrained KP hierarchy. This is the formal-quintic counterpart of the Hodge (C^3-vertex, Gamma-type) factor of the resolved conifold: a universal, explicitly solvable vertex theory carrying the c = 1 gap in its residue direction.
2. The natural variable of genus3_natural.py: in the model the canonical coordinate is exactly -(5/4) zeta^i s^4, so zhat = -(5/4) zeta^i z/u_i. Replacing the model's canonical coordinate by the true one, M_0(tau(u) zhat) is the Hurwitz R-matrix at the point a_5 = 4 u_eff with (5/4) u_eff^{4/5} := Int_0^u t^{-1/5} (1-t)^{-1} dt = |Int_{q_c}^{q} L dq/q|, i.e. u_eff = u tau(u)^{-5/4} = u (2F1(1, 4/5; 9/5; u))^{5/4} (a restatement of the finding in the model's terms, not a new check).
3. The kinematic plane Pi of OU_EXTENSION_genus3.md s.6d in Hurwitz terms: the allowed constant propagators 1 (x) H^2 + H^2 (x) 1 and H^3 (x) H^4 + H^4 (x) H^3 are d_4 (x) d_1 + d_1 (x) d_4 (unit with pole position) and d_2 (x) d_3 + d_3 (x) d_2; the forbidden H (x) H is d_5 (x) d_5, the residue direction, the only one in which the frozen theory has genus >= 1 primary correlators.
4. OPEN (CONJECTURAL / not attempted): whether the u-corrections (the polynomial left factor G(u, zhat) of genus3_natural.py, the O(u) parts of the TFT) have a description inside H_{0;1,4} (a path leaving the residue line, or a non-conformal deformation of the superpotential); whether F_g of H_{0;1,4} on the line equals B_{2g}/(2g(2g-2)) a_5^{2-2g} for all g (expected; here only g <= 3 via the project).

## 7. Caveats

- AI-derived, AI-checked; unrefereed. The proof of s.4 uses the row definitions as documented in fq_rows.py / fq_sser.py / fq_conifold.py and that the conifold-frame map is a differential homomorphism; the leading terms of Cf(C0), Cf(C1) were checked numerically for two frames only (they are immediate from varpi_0 = 1 + O(u), t_c = u + O(u^2)).
- Finite-order evidence: stored symbol zhat^6 (175 residues per prime); re-derived symbol zhat^18 (475 residues per prime); identity row zhat^40 (two colours, both primes); stationary-phase series versus recursion zhat^10 (exact). Residue comparisons are mod p (2^89-1 and 2^61-1) against exact rationals reduced mod p, so no rational reconstruction is involved in any PASS.
- The stationary-phase series are formal; no analytic (Borel/thimble) statement is made, and the overall constants of the integrals (sqrt(2 pi z ...), the orientation sign of the thimble) are normalised away.
- The first version of lg_match.py expected N_a = 1 for all rows and reported FAIL for a = 1, 2 (lg_match_run1.log): the signs N = (1,-1,-1,1,1) were found from the data and then explained in s.4. The task's guess H^a <-> x^a dx is wrong for a != 0 (correct: x^{-a} mod 5).
- Literature: only the items marked as read were read; statements about Dubrovin's notes and Milanov-Tseng are from memory.
- Nothing here was run with more than one process; the two background production runs and their files were not touched.

## 8. Files in this folder

| file | content |
|---|---|
| lg_sp.py | stationary-phase engine (Fractions), validation, J_series, J_general |
| lg_dump_symbol.py -> lg_symbol_89.json, lg_symbol_61.json | dump of sym.Rh and the hatted rows (u-series) from fq3_par.setup |
| lg_match.py -> lg_match.log (lg_match_run1.log: first run) | s.2 blocks (A)-(D) |
| lg_gm.py -> lg_gm.log, lg_symbol_rows_Q.json | (R1)-(R3), (MC), symplecticity to zhat^40, data comparison, F_1; rows g^a_n, n <= 40 |
| lg_wkb_check.py -> lg_wkb_check.log | identity row from the PF WKB recursion to zhat^40 |
| lg_rows_highorder.py -> lg_rows_highorder.log | all rows to zhat^12, zhat^18; frame (1/3, -2/7) to zhat^8 |
| lg_frob.py -> lg_frob.log | s.5 checks (a)-(e) |
| lg_cf_leading.log | leading terms of Cf(C0), Cf(C1) |
| RUNLOG.md | chronological log of the runs |
