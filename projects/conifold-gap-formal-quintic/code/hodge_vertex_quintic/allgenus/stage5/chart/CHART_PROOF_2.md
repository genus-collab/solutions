# Stage 5 / chart, part 2: the formal lemmas (s.2)

Continuation of `CHART_PROOF.md` (s.0-s.1; notation and the definition of `N`, `Nhat`, `SP` in s.1.2). Every lemma holds over the
rings of s.1.3 and with holomorphic parameters (s.1.4), with the same proof. The proofs are algebraic (Wick's theorem), except (L3),
which uses two classical analytic theorems. A parallel agent (`stage5/lemmas/`) writes independent proofs; this part stands alone.
Exact tests of every lemma on seeded random data, with failing controls: `ch_lemmas.py` (s.6).

## 2.1 (S1) Formal integration by parts (PROVED)

**Statement.** For `h` in `O_p[[z]]` and `i = 1..n`: `N[Phi, (d_i Phi) h + z d_i h; p] = 0`, i.e. `SP_p[d_i(e^(Phi/z) h) d^n x] = 0`.
Hence `SP_p` kills every exact form `d(e^(Phi/z) g)` (s.1.5).

*Proof.* Put `T(eta) := exp(Phi_3(z^(1/2) eta)/z) h(p + z^(1/2) eta, z)`, a series in `z^(1/2)` with polynomial coefficients.
Since `(d_i Phi)(p + xi) = (H xi)_i + d_i Phi_3(xi)` and `d/d eta_i = z^(1/2) d/d xi_i`,
`exp(Phi_3/z) [(d_i Phi) h/z + d_i h] = z^(-1/2) [(H eta)_i T + d T/d eta_i]`. Gaussian integration by parts,
`<(H eta)_i P + d_i P>_C = 0` for every polynomial `P`, follows from `<eta_j P>_C = Sum_k C_jk <d_k P>_C` and `H C = -1`. Apply it to
each coefficient of `z^(r/2)`. QED

## 2.2 (PG) The perturbed Gaussian (PROVED)

**Statement.** Let `R` be a commutative Q-algebra, complete for an ideal `I`, `H` an invertible symmetric matrix over Q (or C),
`C = -H^-1`, and `delta` a symmetric matrix with entries in `I`. Put `C_delta := -(H + delta)^-1` (Neumann series). Then for every
polynomial `P` over `R`: `< exp((1/2) eta^T delta eta) P >_C = det(1 + H^-1 delta)^(-1/2) < P >_(C_delta)`, both sides `I`-adically
convergent, the determinant power defined by the binomial series.

*Proof.* Replace `delta` by `s delta` (`s` a formal variable); put `Q := (1/2) eta^T delta eta`, `Lam_s(P) :=` left side,
`rho_s(P) :=` right side, with `C(s) := -(H + s delta)^-1` and `d(s) := det(1 + s H^-1 delta)^(-1/2)`. Then `Lam_0 = rho_0 = < >_C`
and `d/ds Lam_s(P) = Lam_s(Q P)`. For `rho`: `dC/ds = C delta C`; since `<P>_C = [exp((1/2) Sum C_ij d_i d_j) P](0)`,
`d/ds <P>_(C(s)) = (1/2) Sum_ij (C delta C)_ij <d_i d_j P>_(C(s))`; applying the Isserlis recursion twice,
`<Q P>_(C(s)) = (1/2) tr(delta C) <P> + (1/2) Sum_ij (C delta C)_ij <d_i d_j P>`; and `d/ds log d(s) = (1/2) tr(C(s) delta)`.
Together `d/ds rho_s(P) = rho_s(Q P)`. By induction both sides have the same `s`-Taylor coefficients `<Q^k P>_C`; set `s = 1`.
Both sides are `I`-adically continuous and linear in `P`, so the identity extends to the `I`-adic completion of `R[eta]`. QED
The infinitesimal form, used in (L1c): `d/dt [det(-H_t)^(-1/2) <P>_(C_t)] = det(-H_t)^(-1/2) <(1/2) eta^T (dH_t/dt) eta P>_(C_t)`.

## 2.3 (L0) Products (PROVED)

If `x = (x', x'')`, `p = (p', p'')`, `Phi = Phi'(x') + Phi''(x'')` and `A = A'(x', z) A''(x'', z)`, then
`N[Phi, A; p] = N[Phi', A'; p'] N[Phi'', A''; p'']` and `SP` factorises (with `det(-H)^(-1/2) = det(-H')^(-1/2) det(-H'')^(-1/2)`).
*Proof.* The covariance is block diagonal, `exp(Phi_3/z) A` is a product of functions of `eta'` and of `eta''`, and Gaussian moments
of products of polynomials in disjoint blocks factorise. QED

## 2.4 (L1a) Invariance under analytic changes of variables (PROVED)

**Statement.** Let `chi: (C^n, p') -> (C^n, p)` be a biholomorphic germ, `chi(p') = p`, `J := det D chi`. Then `p'` is a
non-degenerate critical point of `Phi o chi`, with Hessian `H' = D chi(p')^T H D chi(p')`, and
`N[Phi o chi, (A o chi) J; p'] = J(p') N[Phi, A; p]`. Equivalently `SP_(p')[chi^*(e^(Phi/z) A d^n x)] = SP_p[e^(Phi/z) A d^n x]` with
the branch `det(-H')^(-1/2) := det(-H)^(-1/2)/J(p')`. In particular `Nhat` is invariant.

*Proof.* *Translations:* `N` is defined by the Taylor data at the critical point; nothing changes. *Linear maps* `chi = M`: the
Taylor data transform tensorially, `C' = M^-1 C M^-T`, and `<P(M eta')>_(C')= <P(eta)>_C` (a linear image of a Gaussian, Isserlis);
the Jacobian is the constant `det M`, so `N' = det M N`. *General `chi`:* write `chi = M o chit` with `chit` tangent to the identity;
it remains to treat `chit`, with `p = p' = 0`. Put `chi_t(x) := x + t (chit(x) - x)`, `t` in C. Then `D chi_t(0) = 1` for all `t`, so
`chi_t` is a biholomorphic germ, `Phi_t := Phi o chi_t` has `Phi_t(0) = Phi(0)`, `grad Phi_t(0) = 0`, `Hess Phi_t(0) = H`, and the
Taylor coefficients of `Phi_t` and of `A_t := (A o chi_t) J_t` at 0 are polynomials in `t`. So `I(t) := N[Phi_t, A_t; 0]` has
polynomial coefficients; since the Gaussian does not move and `d_t Phi_t = O(|x|^3)`, differentiating the finite sums of s.1.2 gives
`I'(t) = N[Phi_t, z^-1 (d_t Phi_t) A_t + d_t A_t; 0]`. Let `Y_t := (chit - id) o chi_t^-1` (a holomorphic vector field germ), so
`d_t chi_t = Y_t o chi_t`. For `f := e^(Phi/z) A`, the chain rule and Jacobi's formula `d_t J_t = J_t (div Y_t) o chi_t` give
`d_t [(f o chi_t) J_t] = [div(f Y_t)] o chi_t . J_t`, and naturality of the Lie derivative
(`chi_t^*(L_Y omega) = L_(chi_t^* Y) chi_t^* omega`, `L_Y (g d^n x) = div(g Y) d^n x`) turns the right side into
`Sum_i d_i(e^(Phi_t/z) g_i)` with `g := A_t . D chi_t^-1 (Y_t o chi_t)` in `O_0[[z]]^n`. By (S1), `I'(t) = 0`. So `I(1) = I(0)`. QED

## 2.5 (L1b) Fubini for iterated stationary phase (PROVED)

**Statement.** Let `x = (v, w)` in `C^(n1) x C^(n2)`, `p = (v_0, w_0)` a non-degenerate critical point of `Phi` with `H_vv := d_v^2 Phi(p)`
invertible.
* (i) There is a unique holomorphic germ `v_c(w)`, `v_c(w_0) = v_0`, with `d_v Phi(v_c(w), w) = 0`; `phi(w) := Phi(v_c(w), w)` has a
  critical point at `w_0` with `Hess phi(w_0) = S := H_ww - H_wv H_vv^-1 H_vw` (Schur complement); `det H = det H_vv det S`, so `S` is
  invertible; `H_vv(w) := d_v^2 Phi(v_c(w), w)` is invertible near `w_0`.
* (ii) `B(w, z) := N[Phi(., w), A(., w, .); v_c(w)](z)` lies in `O_(w_0)[[z]]` (each coefficient holomorphic in `w`).
* (iii) With `beta(w) := (det H_vv(w)/det H_vv(w_0))^(-1/2)` (the branch with `beta(w_0) = 1`):

      N[Phi, A; p] = N[phi, beta B; w_0].                                                         (Fub)

  Equivalently `SP_p = SP_(w_0)[ SP_(v_c(w))[ . ] ]`, with the branch `det(-H)^(1/2) = det(-H_vv)^(1/2) det(-S)^(1/2)`.

*Proof.* (i) Implicit function theorem; differentiating `d_v Phi(v_c(w), w) = 0` gives `v_c' = -H_vv^-1 H_vw` at `w_0`; the chain rule
gives `phi'(w) = d_w Phi(v_c(w), w)` and then `phi''(w_0) = H_ww + H_wv v_c'(w_0) = S`; the block determinant formula gives
`det H = det H_vv det S`. (ii) is (L1c)(i) below.
(iii) By (L1a) with `chi(v', w) = (v_c(w) + v', w)` (Jacobian 1) we may assume `v_c(w) = 0` for all `w`, `w_0 = 0`; the inner data at
fixed `w` are merely translated, so `B` is unchanged. Now `d_v Phi(0, w) = 0` for all `w`, hence `H = diag(H_vv, S)`. Write
`Phi(v, w) = phi(w) + (1/2) v^T H_vv(w) v + Phi_(3,v)(v; w)` with `Phi_(3,v) = O(|v|^3)` at fixed `w`, and
`delta(w) := H_vv(w) - H_vv(0)`. In the scaled variables `(v, w) = z^(1/2)(eta_v, eta_w)` the non-Gaussian part of the total
integrand is `exp([phi(w) - phi(0) - (1/2) w^T S w]/z) . exp((1/2) eta_v^T delta(w) eta_v) . exp(Phi_(3,v)(z^(1/2) eta_v; w)/z) A`, with
`w = z^(1/2) eta_w` everywhere (each term of the three exponents has degree >= 3 in `(v, w)`, as the definition requires).
The total covariance is `diag(C_v, C_w)`, so `< >_C = < < >_(C_v) >_(C_w)` (moments of a block-diagonal Gaussian factorise).
Inner step: apply (PG) in `R = C[eta_w][[z^(1/2)]]`, `I = (z^(1/2))`; the entries of `delta(z^(1/2) eta_w)` lie in `I`. Since `det(1 + H_vv(0)^-1 delta(w)) = det H_vv(w)/det H_vv(0)`, (PG) produces the factor `beta(w)` and replaces `C_v` by
`-H_vv(w)^-1`; what is left inside is exactly the inner Wick expansion at fixed `w`, i.e. `B(w, z)` with `w = z^(1/2) eta_w`
(the inner expansion is a polynomial in the Taylor data at `(0, w)`, which are Taylor-expanded in `w`; both operations commute).
What remains is `< exp([phi(w) - phi(0) - (1/2) w^T S w]/z) beta(w) B(w, z) >_(C_w)`, which is `N[phi, beta B; 0]` by definition. All
rearrangements take place in `R`, with finitely many terms at each order of `z^(1/2)`. QED

## 2.6 (L1c) Parameters (PROVED)

Let `Phi(x, t)`, `A(x, t, z)` be holomorphic near `(p_0, t_0)`, `p_0` a non-degenerate critical point at `t_0`, `p(t)` as in s.1.4.
* **(i) Holomorphic dependence.** `t -> N_n[Phi(., t), A(., t); p(t)]` is holomorphic near `t_0`, and at `t_0` it is `N_n` of the
  `t_0`-data. If the data and `p(t)` are rational in `t` over a field `k`, then `N_n(t)` is in `k(t)`, with poles only where the
  Taylor data at `(p(t), t)` have poles or `det H(t) = 0`. *Proof.* (F1): `N_n(t)` is a polynomial in the Taylor coefficients at
  `(p(t), t)` and in the entries of `-H(t)^-1 = -adj H(t)/det H(t)`.
* **(ii) Taylor expansion in `t`.** The Taylor coefficients of `N_n(t)` at `t_0` are obtained by running the Wick formula over the ring
  `C[[t - t_0]]` with the Taylor-expanded data. *Proof.* (i) and (F1).
* **(iii') Differentiation in a parameter** ("SP commutes with `d_t`"; with `t = log q` this is "SP commutes with `D_q`"):

      d/dt SP_(p(t))[e^(Phi_t/z) A_t] = SP_(p(t))[e^(Phi_t/z) (z^-1 (d_t Phi_t) A_t + d_t A_t)],

  the left side being the `t`-derivative of the symbol `(2 pi z)^(n/2) det(-H_t)^(-1/2) e^(Phi_t(p(t))/z) N_t` (product rule).
  *Proof.* Translate by `tau_t(x) = x + p(t) - p(t_0)` ((L1a), Jacobian 1): `Psi_t := Phi_t o tau_t`, `B_t := A_t o tau_t` have the
  fixed critical point `p_0`. At a fixed critical point the formula holds term by term: the constant `Psi_t(p_0)` gives the
  derivative of `e^(Psi_t(p_0)/z)`; the change of the Hessian is the infinitesimal (PG) (the linear term of `d_t Psi_t` at `p_0`
  vanishes because `grad Psi_t(p_0) = 0` for all `t`); the cubic and higher terms and the amplitude are differentiated directly.
  Then `d_t Psi_t = (d_t Phi_t) o tau_t + p'(t) . (grad Phi_t) o tau_t` and `d_t B_t = (d_t A_t) o tau_t + p'(t) . (grad A_t) o tau_t`;
  the `p'`-terms add up to `SP_(p_0)[Sum_i p'_i d_i(e^(Psi_t/z) B_t)] = 0` by (S1). Translating back gives the formula. QED
* **(iii) A parameter-dependent part of the phase as a perturbation.** If `Phi_t = Phi_0 + t K` and `A_t = A`, iterating (iii') gives
  `d^d/dt^d SP_(p(t))[e^(Phi_t/z) A] = z^-d SP_(p(t))[e^(Phi_t/z) K^d A]`. At `t = 0`: the `d`-th Taylor coefficient in `t` of
  `e^(-Phi_0(p_0)/z) SP_(p(t))[e^(Phi_t/z) A]` (an element of `z^(n/2) C((z))`: each derivative of `e^(Phi_t(p(t))/z)` produces
  finitely many negative powers of `z`) is `(1/d!) z^-d e^(-Phi_0(p_0)/z) SP_(p_0)[e^(Phi_0/z) K^d A]`. So the `t`-expansion of the SP
  at the moving critical point equals, coefficientwise in `t`, the SP at the fixed point `p_0` with `e^(tK/z)` expanded in the
  amplitude. QED (VERIFIED in a one-variable instance with the `q`-dependent critical point, `ch_lemmas.py`.)

## 2.7 (L2) The formal delta lemma (PROVED)

**Statement.** Variables `(sigma, x_0)` in `C^2` and parameters `w` in `C^r`. Let `Phi = sigma G(x_0, w) + f(x_0, w)` (EXACTLY linear in
`sigma`) and let the amplitude `A(x_0, w, z)` not depend on `sigma`. Suppose `G(X_*, w_*) = 0` and `d_(x_0) G(X_*, w_*) != 0`. For `w`
near `w_*` let `X_0(w)` be the solution of `G(X_0(w), w) = 0` near `X_*`, and `sigma_c(w) := -[d_(x_0) f/d_(x_0) G](X_0(w), w)`. Then
`(sigma_c(w), X_0(w))` is a non-degenerate critical point of `Phi(., ., w)`, `det Hess = -(d_(x_0) G)^2`, critical value `f(X_0(w), w)`, and

    N[Phi(., ., w), A(., w, .); (sigma_c(w), X_0(w))](z) = A(X_0(w), w, z)        EXACTLY (no corrections in z),

so `SP = 2 pi z (-(d_(x_0) G)^2)^(-1/2) e^(f(X_0, w)/z) A(X_0, w, z)`: the restriction to `G = 0` times the Gaussian factor.

*Proof.* Fix `w`. Put `s = sigma - sigma_c`, `xi = x_0 - X_0`, `g_k := [d^k G/dx_0^k](X_0)/k!`, `fh_k := [d^k (sigma_c G + f)/dx_0^k](X_0)/k!`.
Since `G(X_0) = 0` and `fh_1 = 0` (definition of `sigma_c`): `Phi - Phi(p) = s Sum_(k>=1) g_k xi^k + Sum_(k>=2) fh_k xi^k`. So
`H = [[0, g_1], [g_1, 2 fh_2]]`, `det H = -g_1^2 != 0`, and `C = -H^-1 = g_1^-2 [[2 fh_2, -g_1], [-g_1, 0]]`: **`C_(xi xi) = 0`.**
Every monomial of `Phi_3 = s Sum_(k>=2) g_k xi^k + Sum_(k>=3) fh_k xi^k` has `s`-degree at most 1, and if it is 1 then its
`xi`-degree is at least 2; `A` has `s`-degree 0. Hence every monomial `s^a xi^b` of `Phi_3^V A_m` has `a <= b/2`. In a Wick pairing
of `s^a xi^b` each `xi` must be paired with an `s` (because `C_(xi xi) = 0`), so `b <= a`. Thus `b <= b/2`, `a = b = 0`: only `V = 0` and
the constant term of `A_m` survive, `N_m = A_m(X_0)`. QED
*Remark.* The same proof in the coordinates `(sigma, G)` (critic, SYN4 s.1.4) needs (L1a); the count above does not.

## 2.8 (L3) Stirling's series is the formal stationary phase of Euler's integral (PROVED modulo two quoted theorems)

**Statement.** For `a` in `C^x` and `M` in Z, the phase `a log x - x` (any branch near `x = a`) has the non-degenerate critical point
`x = a`, Hessian `-1/a`, critical value `a log a - a`, and

    N[a log x - x, x^(M-1); a](z) = a^(M-1) Poch_M(z/a) exp( Sum_(k>=1) b_k (z/a)^(2k-1) ),
    Poch_M(e) := prod_(i=0..M-1) (1 + i e)  (M >= 0),      Poch_M(e) := prod_(i=1..-M) (1 - i e)^-1  (M < 0).

In symbol form: `SP_a[x^(a/z + M) e^(-x/z) dx/x] = z^(a/z + M) Gamma_St(a/z + M)`, where the right side denotes
`(2 pi z)^(1/2) a^(M - 1/2) e^((a log a - a)/z) Poch_M(z/a) exp(Sum_k b_k (z/a)^(2k-1))` (Stirling's formula with `log(a/z) = log a - log z`
and `Gamma_St(w + M) = (w)_M Gamma_St(w)`). Insertion of `x^m` with the measure `dx/x` is `M = m`; with the measure `dx` it is `M = m + 1`.

*Proof.* (1) *Scaling.* By (L1a) with `x = a t` and (F3) with `mu = a`:
`N[a log x - x, x^(M-1); a](z) = a^(M-1) N[log t - t, t^(M-1); 1](z/a)`; the coefficients of the last series are rational numbers. The
right side of the claim is also `a^(M-1)` times a series in `z/a` with rational coefficients. So it suffices to prove the identity of
these two rational sequences, and we may take `a > 0` real.
(2) *The integral.* For real `a > 0`, real `z > 0` with `a/z + M > 0`: `I(z) := Int_0^oo x^(a/z + M - 1) e^(-x/z) dx = z^(a/z+M) Gamma(a/z+M)`
(substitute `x = z t` in Euler's integral).
(3) *Stirling.* DLMF 5.11.1 (read 2026-11-14): `ln Gamma(w) ~ (w - 1/2) ln w - w + (1/2) ln(2 pi) + Sum_(k>=1) B_2k/(2k(2k-1) w^(2k-1))`
as `w -> oo` in `|ph w| <= pi - delta`. With `Gamma(w + M) = (w)_M Gamma(w)` (resp. `Gamma(w)/prod_(i=1..-M)(w - i)`) and `w = a/z` this gives
`I(z) ~ (2 pi z)^(1/2) a^(M - 1/2) e^((a log a - a)/z) Poch_M(z/a) exp(Sum_k b_k (z/a)^(2k-1))` as `z -> 0+`.
(4) *Laplace.* `phi(x) = a log x - x` has on `(0, oo)` the unique maximum `x = a`, `phi''(a) = -1/a`. Outside `U = [a - d, a + d]`,
`phi <= phi(a) - eta` with `eta > 0` and the tails are integrable, so that part of `I(z)` is `O(e^((phi(a) - eta/2)/z))`. On `U` use the
Morse substitution `x = a + mu(u)` (`mu` analytic, `mu(0) = 0`, `phi(a + mu(u)) = phi(a) - u^2/2`); the integrand becomes
`e^(phi(a)/z) e^(-u^2/(2z)) g(u) du`, `g = (a + mu)^(M-1) mu'` analytic. Splitting `u > 0`, `u < 0` and putting `t = u^2/2`, Watson's
lemma (DLMF 2.3.8, read 2026-11-14: `Int_0^oo e^(-x t) q(t) dt ~ Sum_s Gamma((s + lam)/mu) a_s x^(-(s + lam)/mu)` for
`q(t) ~ Sum_s a_s t^((s + lam - mu)/mu)`; here `x = 1/z`, `mu = 2`, `lam = 1`, finite upper limit changing only exponentially small
terms) gives `e^(phi(a)/z) (2 pi z)^(1/2) Sum_k g_2k (2k-1)!! z^k`, which is the Wick series of the Morse normal form; by (L1a) it
equals `e^(phi(a)/z) (2 pi z)^(1/2) a^(1/2) N[phi, x^(M-1); a](z)`.
(5) *Uniqueness* of asymptotic expansions in powers of `z` (after dividing by `(2 pi z)^(1/2) e^(phi(a)/z)`) identifies (3) and (4). QED
*VERIFIED* (`ch_euler.py 40 40`): the identity holds exactly over Q through `z^40` for `m = 0, 1, 2, 3` with both measures and for
`m = -1, -2, -3` with both measures (the negative `m` are used in step (e)), at `a = 1, -2/7, 3/2`; the first coefficients
`1, 1/12, 1/288, -139/51840, -571/2488320, 163879/209018880` are the `g_k` of DLMF 5.11.4. Controls fail: `b_k -> B_2k/(2k)` (first
at `z^3`), the Pochhammer factor dropped (`z^1`), `exp(-Sum ...)` (`z^1`).
*Not done:* a purely formal proof of (L3). (The natural one, the difference equation `Gamma_St(w+1) = w Gamma_St(w)` transported to
`N`, needs a shift of the parameter `a` by the formal variable `z`, which is not a formal operation; I did not pursue it.)
