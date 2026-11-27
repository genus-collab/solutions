# Stage 4, Gamma route: proofs (companion of GAMMA.md, 2026-11-14)

AI-written, unrefereed. Notation of `GAMMA.md` s.1-s.3. `Omega = e^(W/z) A dx` in `x = 1/sigma`, or `e^(Phi(sigma)/z) H(sigma) dsigma` in `sigma`,
`Phi = S_0` (the phase, containing `sigma log(-q)`), `H = exp(S_1 + Sum_k z^(2k-1) S_2k)`. Theorem U, Theorem G, Lemma 1 are those of UR s.3 and
GAUGE s.2-s.3 (quoted, AI-audited).

## 1. The two formal functionals

**Stationary phase `SP_c`.** For `f = e^(phi/z) g`, `phi` analytic near `c`, `phi'(c) = 0 != phi''(c)`, `g` a formal power series in `z` with
coefficients analytic near `c`: Taylor-expand at `c`, put `x = c + z^(1/2) y`, integrate monomials against `e^(phi''(c) y^2/2)` (Gaussian
moments); the result is `e^(phi(c)/z) (2 pi z/(-phi''(c)))^(1/2) (g_0(c) + O(z))`, a series in integer powers of `z`. Standard properties
(formal Laplace method): (S1) `SP_c[d(e^(phi/z) h)] = 0` for `h` with analytic coefficients (formal integration by parts); (S2)
`SP_c[f dx] = SP_(c')[chi^* (f dx)]` for `x = chi(y)` analytic near `c'`, `chi(c') = c`, `chi'(c') != 0`, also for `chi` a formal power series
in a parameter `u` whose `u^0` term is such a map; (S3) smooth dependence on parameters (`u`, `q`).
*(S4) translation invariance:* `SP_mu[F(sigma - z) dsigma] = SP_mu[F(sigma) dsigma]` for `F = e^(Phi/z) g`. *Proof.* `F(sigma - eps) =
e^(Phi(sigma)/z) g_eps(sigma)`, `g_eps = e^((Phi(sigma - eps) - Phi(sigma))/z) g(sigma - eps)`, a formal series in `eps` whose `eps^n` coefficient has
`z`-valuation `>= -n/2` at the saddle (`Phi'(mu + t) = O(t)`, `t ~ z^(1/2)`). `G(eps) := SP_mu[F(sigma - eps) dsigma]` satisfies `dG/deps = -SP_mu[d(F(sigma - eps))] = 0`
by (S1), so `G(eps) = G(0)` coefficientwise; the substitution `eps -> z` is legitimate on series with that valuation bound. QED

**Residue `Res`.** For `omega = e^(W/z) f dx`, `W = l/x + G(x)`, `l = log(1-u) in u Q[[u]]`, `f` a formal series in `z` with coefficients in
`Q[[u]]((x))`: expand `e^(l/(xz)) = Sum_k l^k x^-k z^-k/k!`; at each order `u^N` only `k <= N` occur and `[x^m] e^(G/z)` is a polynomial in `1/z`,
so `Res_(x=0)` is a finite computation per order. (R1) it kills `d(e^(W/z) h)`, `h` Laurent in `x` (residue of a derivative of a function
single-valued at each order); (R2) it is invariant under `x = chi(y)`, `chi(y) = y(1 + O(y))` analytic, or a formal perturbation of the
identity in `u`; (R3) in particular under `sigma -> sigma - z`, i.e. `x -> x/(1 - z x)`.

## 2. Both functionals produce solutions of the quintic system

*Difference equation.* For the exact integrand `a(s) (-q)^s`: `(1 - sigma^5) a(s) = prod_(k=0..4)(5 sigma - k z) a(s-1)` (from `Gamma(y+1) = y Gamma(y)`).
The full-Stirling integrand satisfies the same identity as formal series, because each `Gamma_St` satisfies `Gamma_St(y+1) = y Gamma_St(y)`
in `Q[[1/y]]` (classical: `(y + 1/2) log(1 + 1/y) - 1 + Sum_k b_k ((y+1)^(1-2k) - y^(1-2k)) = 0`).
*PF.* Write the integrand as `a(s)(-q)^s ds` (`a` exact or `a_St`), `psi = F[.]` (`F = SP_c` or `Res`). Since `zD (-q)^s = sigma (-q)^s`:
`prod_j (zD - zeta^j) psi = F[(sigma^5 - 1) a(s)(-q)^s]` and `q prod_(k=1..5)(5zD + kz) psi = -F[prod_k (5 sigma + kz) a(s) (-q)^(s+1)]`, which by the
shift `s -> s - 1` ((S4), resp. (R3)) is `-F[prod_(k=0..4)(5 sigma - kz) a(s-1)(-q)^s]`. So `PF psi = F[{(sigma^5 - 1) a(s) + prod_k (5 sigma - kz) a(s-1)}(-q)^s] = 0`.
*Rows.* `zD` acts on `Omega` as multiplication by `x^-1` (`A` does not depend on `q`). With the row recursion of GAMMA s.3 the vector
`(F[Omega_a])_a` solves `zD Y_a = c_a Y_(a+1)` for `a = 0..3` by construction and for `a = 4` by PF (UR s.1: the closing relation is PF).
So `Q^fS := (SP_(colour i)[Omega_a])` is a formal WKB fundamental matrix with exponentials `e^((u_i - u_i^c)/z)`: the critical values of
the single-valued phase are exactly `W(zeta^-i s) = zeta^i s^4 What(1; u) = -zeta^i (5/4) u^(4/5)/tau = u_i - u_i^c` (s.4), as for the model at
`v(u)`. `Y^MB := (Res[Omega_a])` lies in `Q[z^-1][[u]]^5` (only negative powers: `[x^(k-1)] A` has `z`-degree `<= k - 1`).

## 3. Hatted structure; (H) for the full-Stirling WKB

`x = s xi`, `z = s^4 zhat`: `W/z = What(xi; u)/zhat` with `What = (l/u) xi^-1 - Sum_l u^l xi^(5l+4)/((l+1)(5l+4))` and
`A = (1 - u xi^5)^(-1/2) exp(Sum_k b_k (u zhat xi)^(2k-1) [5^(1-2k) - 5 Sum_(r>=0) binom(2k-2+5r, 5r) u^r xi^(5r)])`
(`Sum_j (zeta^j x - 1)^(1-2k) = -5 Sum_(5|r) binom(2k-2+r, r) x^r`). Both are power series in `u` with Laurent-polynomial coefficients in `xi`
(polynomial in `zhat`), and `dWhat/dxi = 0` at `xi = 1` for all `u` (`xi^5 = 1`). So the colour-0 saddle series is, at each order in `zhat`, a
power series in `u`; the row operators preserve this (UR s.4). This is (H) for `Q^fS`. With UR s.4: `Q^fS D^fS = P_can Q^m` for a constant
diagonal `D^fS`, which (H) forces into `Q[[z^5]]` and Z/5 makes scalar: `P^fS := Q^fS (Q^m)^-1 = P_can / c^fS(z)`, regular (Theorem U).
Consequence: `c_n(L) = L^(4n) ph^fS_n(L^-5)` with `ph^fS_n` a power series, hence (WKB recursion) a polynomial: `c_n in Q[L^(+-1)]`, exponents `= 4n mod 5`.

## 4. Lemma CFU (hatted Chester-Friedman-Ursell map)

*Statement.* There is a unique `xi = Xi psi(Xi^5; u)`, `psi = 1 + Sum_(k>=1) u^k phi_k(y)`, `phi_k in Q[y]`, with
`What(Xi psi; u) = M(Xi; u) := -(Xi^4/4 + Xi^-1)/tau(u)`, `1/tau = 2F1(1, 4/5; 9/5; u) = Sum_n 4u^n/(5n+4)`. Moreover `phi_k(1) = 0`.
*Proof.* Multiply by `Xi`: `F(psi) := (l/u)/psi + Sum_l g_l u^l y^(l+1) psi^(5l+4) + (1 + y/4)/tau = 0`, `g_l = -1/((l+1)(5l+4))`, `y = Xi^5`. At `u = 0`,
`psi = 1` solves it (`-1/psi - (y/4) psi^4 + 1 + y/4`). Given `psi` mod `u^k`, the `u^k` coefficient reads `(1 - y) phi_k = r_k(y)` with `r_k` a
polynomial (all terms are `Xi^-1` times polynomials in `Xi^5`), since `dF/dpsi|_(u=0,psi=1) = 1 - y`. Solvable in `Q[y]` iff `r_k(1) = 0`.
Induction hypothesis `phi_j(1) = 0` (`j < k`): then `xi(1; u) = 1 + O(u^k)` and, `xi = 1` being a critical point of `What` for all `u`,
`[u^k] What(xi(1; u); u) = [u^k] What(1; u)`; and `What(1; u) = l/u - Sum_l u^l/((l+1)(5l+4)) = -Sum_n u^n (1/(n+1) + 1/((n+1)(5n+4))) =`
`-Sum_n 5u^n/(5n+4) = -(5/4)/tau(u) = M(1; u)`, so `r_k(1) = 0`. Differentiating the order-`k` identity in `Xi` at `Xi = 1` gives
`What''(1; 0) phi_k(1) = 0` with `What''(1; 0) = -5`, so `phi_k(1) = 0`. Uniqueness: the homogeneous equation `(1-y) phi = 0`. QED
*VERIFIED* `gm_cfu.py 12`: `phi_1 = (1-y)/18`, `phi_2 = 13/567 - (101/4536) y - y^2/1512`, `deg phi_k = k` (`k <= 12`), identity mod `u^13`.

## 5. Proof of Theorem A

(a) *Pull-back.* By Lemma CFU, `Omega_a = s e^(M(Xi)/zhat) B_a(Xi; u, zhat) dXi` with `B_a` a Laurent polynomial in `Xi` at each order in `u`.
(b) *Brieskorn reduction* for `M`: `d(e^(M/zhat) Xi^(m+1)) = e^(M/zhat)[(m+1) Xi^m - (Xi^(m+4) - Xi^(m-1))/(tau zhat)] dXi`, so modulo exact forms
`Xi^n == Xi^(n-5) + (n-3) tau zhat Xi^(n-4)` (`n >= 4`) and `Xi^(m-1) == Xi^(m+4) - (m+1) tau zhat Xi^m` (`m <= -1`); every Laurent monomial
reduces to `{Xi^-1, .., Xi^3}` with coefficients in `Q[tau, zhat]`. Hence `Omega_a == Sum_b Phat_ab omega^m_b` mod `d(e^(M/zhat) Q[Xi^(+-1)])`, `omega^m_b`
the model's row forms (`Y = (I_0, -I_-1, -I_3, I_2, I_1)`, rescaled to `Xi`), `Phat_ab` power series in `u` with polynomial coefficients in `zhat`.
(c) *Intertwining.* By (S1), (S2) at the five points `Xi^5 = 1` (the colour saddles; the map is Z/5-equivariant and fixes them, Lemma CFU)
and (R1), (R2): `Q^fS = P^MB Q^m_nat` and
`Y^MB = P^MB Y^C`, with ONE matrix `P^MB` (flat form of `Phat`). The first identity gives `P^MB = Q^fS (Q^m_nat)^-1 = P_can/c^fS` (s.3; the
elementary normalisations of `Q^m` and `Q^m_nat` differ by colour-wise constants, absorbed by `D(0) = 1`).
(d) *Normalisation.* Theorem G: `P^MB Y^C = (1/c^fS) Y^C_Q`. Row 1 at `u = 0`: `Y_1 = zD(psi/C_0)/C_1`, `C_0(0) = 1`, `C_1(0) = -1` in every
conifold frame, and `zD psi|_(u=0) = Res_(x=0)[x^-1 e^(G/z) A dx] = A(0; z) = 1`. So `Y^MB_1(0) = -1 = Y^C_Q,1(0)` and `c^fS = 1`.
Therefore `Y^MB = Y^C_Q` (i) and `Q^fS = P_can Q^m_nat` (ii). Corollary A: the canonical colour-0 identity row is `c(z; L)` (normalised
`c_0 = 1`), so `D = c/R` and (H) + (N) <=> (Gamma). QED
*Remark.* (d) contains the identity behind the numbers: at `u^2 z^-1`, `-1/2` from `l = log(1-u)` and `-1/5 = b_1 (1/5 - 5)/2` from the
Stirling term give `-7/10`, the `u^2` coefficient of the vanishing period; generally `[z^-1] Res Omega = Sum_k alpha_(k-1) l^k/k!`
(`alpha` = coefficients of the `z`-top part `A_top(y) = exp(Sum_k b_k (5^(1-2k) - 5) y^(2k-1))` of `A` at `zx = y`), which is `-varpi_0 t_c`.

## 6. Identity (Gamma), formal route (single-Gamma lemma PROVED; matching SKETCH)

*Single-Gamma lemma (PROVED).* Let `F_r(Y) := SP[Gamma_St(eta) Y^-eta eta^r d eta]` at the saddle `eta = Y` (large `Y`, expansion in `1/Y`), normalised
by `1/(2 pi i)` as for the Cahen-Mellin integral. Then `F_r = (-Y d/dY)^r e^(-Y)`. *Proof.* `-Y F_0' = SP[eta Gamma_St(eta) Y^-eta] =`
`SP[Gamma_St(eta + 1) Y^-eta] = Y SP[Gamma_St(eta) Y^-eta] = Y F_0` by the functional equation and translation invariance (S4) in `eta`; so
`F_0 = C e^(-Y)` with `C` a constant (no `1/Y` corrections survive), `C = 1` from the Gaussian term; `F_r = (-Y d/dY)^r F_0`. QED
So for every polynomial weight the Stirling saddle series equals the exact residue sum `Sum_d Res_(eta = -d) Gamma(eta) Y^-eta eta^r`.
*Matching (SKETCH).* Near `sigma = 1` put `sigma = 1 - z eta`; the full-Stirling integrand is `const * Gamma_St(eta) Y^-eta e^(V(eta, z))`,
`V = Sum_(m>=1) z^m P_m(eta)`, `deg P_m <= m + 1`, `Y = exp(dPsi_0/dsigma(1))` (`Psi` = the other five factors), constant prefactor `h(z) *`
elementary. In the ring `Q[[alpha, beta]]`, `alpha = z/eps = 1/Y`, `beta = eps^2/z` (`eps = L - 1`, so `z = alpha^2 beta`, `eps = alpha beta`):
(i) the outer series `c(z; 1 + eps) = Sum c_(n,k) z^n eps^k` lies in it (`alpha^(2n+k) beta^(n+k)`), and so do the individual outer Feynman
terms (poles `(z/eps)^p`, `p <= n`); (ii) at the saddle `z^m P_m(eta) = O(alpha^(m-1) beta^m)`, so `e^V` is a convergent perturbation and the
outer saddle series equals `Sum_n SP[Gamma_St Y^-eta V^n/n!]` (reorganisation lemma: both are the same convergent sum of Gaussian moments in
the ring); (iii) by the single-Gamma lemma this is `h(z) e^(-Y) Sum_m z^m Qt_m(Y)` with polynomials `Qt_m`, the inner (Frobenius) expansion of
the QRR WKB; comparing the `eps^0` parts gives `c(z; 1) = h(z)`. Step (ii) is the only part not written out.

## 7. Identity (Gamma), analytic route (SKETCH)

*(A1) PROVED modulo routine estimates.* `Z_0 := Gamma(1+5/z) prod_(j!=0) Gamma((zeta^j - 1)/z) q^(1/z) f_0(q, z)` (exact). With the reflection
formula `Gamma((1-sigma)/z) = -2 pi i e^(i pi (1-sigma)/z)/((1 - e^(2 pi i (1-sigma)/z)) Gamma((sigma-1)/z + 1))`, `Z_0 = const * Sum_(d>=0) E(1 + d z)`,
`E(sigma) = e^(i pi(1-sigma)/z) Gamma((sigma-1)/z + 1)^-1 Gamma(5 sigma/z + 1) prod_(j!=0) Gamma((zeta^j - sigma)/z) (-q)^(sigma/z)`, and `E(1 + dz) = 0` for
`d <= -1`. For real `z > 0`, `0 < q < q_c`: `E > 0` on `[1, inf)` (conjugate pairs `j, 5-j`), its phase has `Phi_E'' = 5/(sigma(1 - sigma^5)) < 0`
for `sigma > 1` and a unique maximum at `sigma = L`. Euler-Maclaurin in `d` (summand analytic, all endpoint terms exponentially small,
remainder `O(z^K)` relative) and Laplace's method with Stirling uniform on `[1 + delta, M]` give `Z_0 ~ (elementary) e^(u_0/z) L c(z; L)`.
*(A2) SKETCH (classical input).* The Frobenius-normalised WKB is uniformly valid down to the regular singular point `q = 0`, whose
leading matrix has the distinct eigenvalues `zeta^j`; then `Z_0 ~ (elementary) e^(u_0/z) L R(z; L)` with the QRR series `R` (the Gamma
prefactor of `Z_0` contributes `h(z)`, and the Frobenius-normalised WKB `R/h` has the value 1 at `q = 0`).
(A1) + (A2) at one `q` give `c = R`. A concrete form of (A2) is the uniform representation `q^d A_d(z) = (X^d/d!) exp(Theta(dz, z))`,
`X = 5^4 q/z`, `Theta` smooth in `t = dz` on `[0, T]` uniformly in `z` (Euler-Maclaurin in the product index `r`), which reduces (A2) to a
matched-asymptotics statement for `Sum_d X^d e^(Theta(dz,z))/d!`; the error control in the overlap `1 << X << z^(-1/2)` is not written out.

## 8. Why no ODE-formal argument can replace s.6-s.7

Any two colour-0 WKB series differ by a constant factor `D(z)`; the PF, symplecticity (`D(z) D(-z) = 1` only), Z/5, (H) and the polynomial
structure in `L` are all invariant under `D -> D * exp(odd series in z^5)`. Theorem A removes the conifold end by an exact identity (the
residue at the pole); (Gamma) needs one more exact input at `q = 0`: the residue sum of the single singular Gamma (s.6) or analysis (s.7).
