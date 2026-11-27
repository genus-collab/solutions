# Stage 5 / lemmas, part 4: Stirling (L3) and the formal residue (FL2, FL3)

Continuation of `FORMAL_SP_LEMMAS.md`, `_2.md`, `_3.md`. Same notation and labels.

## 9. Lemma L3: Stirling's series is the formal SP of Euler's integral (PROVED)

### 9.1 Statement
Let `a in A^x` and `m in Z`. The phase `F(x) = a log x - x` is taken as its Taylor series at `x = a`,
`F(a + t) - F(a) = a log(1 + t/a) - t = -t^2/(2a) + t^3/(3a^2) - ..`; the constant `F(a)` enters only through the symbol
`[e^(F(a)/z)] = [a^(a/z) e^(-a/z)]`. So `x = a` is a nondegenerate critical point, `F''(a) = -1/a`, `C = a`. Put (`b_k = B_2k/(2k(2k-1))`)
`St_m(w) := exp( (1/w + m - 1/2) log(1 + m w) - m + Sum_(k>=1) b_k w^(2k-1) (1 + m w)^(1-2k) ) in Q[[w]]`
(the exponent has no constant term, since `(1/w) log(1 + m w) = m + O(w)`). Then:
* **(i)** `SP^norm_a[F; (x/a)^(m-1)] = St_m(z/a)`.
* **(ii)** `SP_a[x^(a/z+m) e^(-x/z) dx/x] = [a^(a/z) e^(-a/z)] (2 pi z)^(1/2) a^(m-1/2) St_m(z/a)`; and `SP_a[x^(a/z+m) e^(-x/z) dx]` is the same with
  `m + 1` in place of `m` (the measure `dx` is `x dx/x`).
* **(iii)** With `Gamma_St(Y) := exp((Y - 1/2) log Y - Y + (1/2) log 2pi + Sum_k b_k Y^(1-2k))` and the conventions
  `log(a/z + m) := log a - log z + log(1 + mz/a)`, `z^(a/z+m) := exp((a/z + m) log z)`, the expansion of `z^(a/z+m) Gamma_St(a/z + m)` is the
  right side of (ii): the `log z` terms combine to `(1/2) log z`, the `log a` terms to `(a/z + m - 1/2) log a`, and the rest is
  `-a/z + (1/2) log 2pi + log St_m(z/a)` (use `(a/z + m)^(1-2k) = (z/a)^(2k-1) (1 + mz/a)^(1-2k)`). So the formal SP of
  `x^(a/z+m) e^(-x/z) dx/x` at `x = a` is `z^(a/z+m) Gamma_St(a/z+m)`, and that of `x^(a/z+m) e^(-x/z) dx` is `z^(a/z+m+1) Gamma_St(a/z+m+1)`.
* **(iv)** *(Formal functional equation.)* `rho(Y) := (Y + 1/2) log(1 + 1/Y) - 1 + mu(Y+1) - mu(Y) = 0` in `Q[[1/Y]]`, `mu(Y) := Sum_k b_k Y^(1-2k)`;
  this is `Gamma_St(Y+1) = Y Gamma_St(Y)`. Equivalently `St_(m+1)(w) = (1 + m w) St_m(w)` for every `m in Z`.
* **(v)** *(Recursion.)* `SP_a[x^(a/z+m+1) e^(-x/z) dx/x] = (a + m z) SP_a[x^(a/z+m) e^(-x/z) dx/x]` for every `m in Z`.

*Bookkeeping for (iv).* Substituting `Y = 1/w + m` (so `1/Y = w/(1 + m w) in w Q[[w]]`, a continuous substitution),
`log(St_(m+1)(w)/St_m(w)) = log(w Y) + rho(Y) = log(1 + m w) + rho(1/w + m)`; so `rho = 0` gives the recursion for all `m`, and the case
`m = 0` (`St_1 = St_0`) alone gives `rho = 0`.

### 9.2 Proof of (iv) (PROVED, formal)
`D := d/dY` acts on `Q[[1/Y]]` raising the order in `1/Y`, so power series in `D` act, compose as power series, and `e^D f(Y) = f(Y+1)`
(compare coefficients of `(Y+1)^-k = Y^-k (1 + 1/Y)^-k`). `D^2` is injective on `Y^-1 Q[[1/Y]]`, and `rho in Y^-2 Q[[1/Y]]` (each part is
`O(Y^-2)`). So it suffices to show `D^2 rho = 0`.
* `D^2 mu = Sum_k B_2k Y^(-2k-1) = T - 1/Y - 1/(2Y^2)`, where `T := (D/(e^D - 1))(1/Y) = Sum_n (B_n/n!) D^n (1/Y)`
  (`D^n(1/Y) = (-1)^n n! Y^(-n-1)`, `B_1 = -1/2`, `B_odd>1 = 0`).
* `(1 - e^D) T = -D(1/Y) = 1/Y^2`; `(1 - e^D)(1/Y) = 1/(Y(Y+1))`; `(1 - e^D)(1/Y^2) = 1/Y^2 - 1/(Y+1)^2`.
  Hence `(e^D - 1) D^2 mu = -(1/Y^2 - 1/(Y(Y+1)) - (1/2)(1/Y^2 - 1/(Y+1)^2)) = -(1/2)(1/Y - 1/(Y+1))^2 = -1/(2 Y^2 (Y+1)^2)`.
* With `l := log(1 + 1/Y)`, `l' = -1/(Y(Y+1))`, `l'' = 1/Y^2 - 1/(Y+1)^2`:
  `D^2[(Y + 1/2) l] = 2 l' + (Y + 1/2) l'' = (-2Y(Y+1) + 2(Y + 1/2)^2)/(Y^2 (Y+1)^2) = 1/(2 Y^2 (Y+1)^2)`.
So `D^2 rho = 1/(2Y^2(Y+1)^2) - 1/(2Y^2(Y+1)^2) = 0`. QED (Classically: `Gamma_St` is the asymptotic expansion of `Gamma`, DLMF 5.11.1, and
`Gamma(Y+1) = Y Gamma(Y)`; the proof above needs no analysis.)

### 9.3 Proof of (i)-(iii) (PROVED, formal; Proof A)
*(H) Homogeneity.* `[z^k] SP^norm_a[F; (x/a)^e] = c_(k,e) a^-k` with `c_(k,e) in Q` independent of `A` and `a`. *Proof:* in Lemma U each
term is `(1/r!)[t^alpha](F_3^r (1 + t/a)^e) <t^alpha>_C` with `|alpha| = 2(k + r)`; a coefficient of `t^l` in `F_3` carries `a^(1-l)`, one in
`(1 + t/a)^e` carries `a^-l`, and `<t^alpha>_C = (|alpha| - 1)!! a^(|alpha|/2)`; the total power of `a` is `r - |alpha| + |alpha|/2 = -k`. QED So
`N_m := SP^norm_a[F; (x/a)^(m-1)] = psi_m(z/a)` for universal `psi_m in 1 + w Q[[w]]`.
*Step A1 ((v); FL1).* FL1 in one variable with `h = x^m = a^m (1 + t/a)^m`: `h' + z^-1 F' h = m x^(m-1) + z^-1 (a x^(m-1) - x^m)`, so
`SP^norm_a[F; x^m] = (a + m z) SP^norm_a[F; x^(m-1)]`. This is (v), and in normalised form `N_(m+1) = (1 + m w) N_m`, `w = z/a`.
*Step A2 (the a-shift; L1c(iv)).* Apply L1c(iv) with `g(x) := log(x/a)` (`g(a) = 0`) and the amplitude `(x/a)^-1`. Left side:
`SP^norm_a[F; e^g (x/a)^-1] = N_1 = N_0` (A1 with `m = 0`). Right side: `F_eps = (a + eps) log x - x - eps log a` has the critical point
`a + eps`; `F_eps(a+eps) - F(a) = (a + eps) log(1 + eps/a) - eps`; `Delta_eps/Delta = a/(a + eps)`; and by (H) over `A[[eps]]` (the unit
`a + eps`; a constant in the phase does not matter) `SP^norm_(a+eps)[F_eps; (x/a)^-1] = (1 + eps/a)^-1 psi_0(z/(a + eps))`. So
`S(eps; z) = exp(((a + eps) log(1 + eps/a) - eps)/z) (1 + eps/a)^(-1/2) psi_0(z/(a + eps))`, and at `eps = z` (`w = z/a`):
`psi_0(w) = (1 + w)^(1/w + 1/2) e^-1 psi_0(w/(1 + w))`,  with `(1 + w)^(1/w+1/2) e^-1 := exp((1/w + 1/2) log(1 + w) - 1) in 1 + w Q[[w]]`.
*Step A3 (uniqueness).* `Lam(Y) := log psi_0(1/Y) in Y^-1 Q[[1/Y]]`. With `w = 1/Y`, `w/(1+w) = 1/(Y+1)`, Step A2 reads
`Lam(Y) - Lam(Y + 1) = (Y + 1/2) log(1 + 1/Y) - 1`. The operator `Lam -> Lam(Y) - Lam(Y+1)` sends `Y^-k` to `k Y^(-k-1) + O(Y^(-k-2))`, so it is
injective on `Y^-1 Q[[1/Y]]`; by (iv) `mu` solves the same equation. Hence `Lam = mu` and `psi_0(w) = exp(Sum_k b_k w^(2k-1)) = St_0(w)`.
*Step A4 (all m).* `N_0 = St_0`, `N_(m+1) = (1 + mw) N_m` (A1) and `St_(m+1) = (1 + mw) St_m` ((iv)); `1 + mw` is a unit, so `N_m = St_m` for
all `m in Z`. This is (i); (ii) is Definition 1.3 with `Delta^(-1/2) = a^(1/2)` and the amplitude value `a^(m-1)`; (iii) is the bookkeeping
stated in s.9.1. QED

### 9.4 The classical derivation by Laplace's method (Proof B of N_0 = St_0; PROVED modulo the cited statements)
Sources read for this note (DLMF, NIST Digital Library of Mathematical Functions, dlmf.nist.gov, accessed 2026-11-14):
**DLMF s.2.3(iii)**, Laplace's method, eqs. (2.3.13)-(2.3.16): for `I(x) = int_a^b e^(-x p(t)) q(t) dt` with the minimum of `p` at the left
endpoint, hypotheses (a) `p'`, `q` continuous near `a` and the minimum approached only at `a`, (b) `p(t) ~ p(a) + Sum p_s (t-a)^(s+mu)`,
`q(t) ~ Sum q_s (t-a)^(s+lam-1)` with `p_0 > 0` (the `p`-expansion differentiable), (c) absolute convergence for large `x`; conclusion
(2.3.15) `I(x) ~ e^(-x p(a)) Sum_s Gamma((s+lam)/mu) b_s x^(-(s+lam)/mu)`, the `b_s` defined by (2.3.16) `q(t)/p'(t) ~ Sum b_s v^((s+lam-mu)/mu)`,
`v = p(t) - p(a)`. **DLMF 5.11.1** (Stirling's series): `Ln Gamma(y) ~ (y - 1/2) ln y - y + (1/2) ln 2pi + Sum_k B_2k/(2k(2k-1) y^(2k-1))` as
`y -> oo` in `|ph y| <= pi - delta`. **DLMF s.2.1(iii)**: asymptotic expansions (2.1.13)-(2.1.14) are unique; addition, multiplication,
"substitution, logarithms, and powers" are permitted.
*Derivation.* By (H) take `a = 1`, `z = 1/y`, `y > 0`. Exactly, `J(y) := int_0^oo e^(-y(x - log x)) x^-1 dx = y^-y Gamma(y)`. Put `p(x) = x - log x`
(minimum 1 at `x = 1`, `p'' = 1`). By the analytic Morse lemma (`p - 1 = (x-1)^2 u(x)`, `u(1) = 1/2`) there is a real-analytic
`x = X(w)`, `X(0) = 1`, `X'(0) = sqrt 2`, with `p(X(w)) = 1 + w^2`; write `X'(w)/X(w) = Sum g_n w^n`. Split `J` at `x = 1`. On `[1, oo)` (2.3.15)
applies with `mu = 2`, `lam = 1`, `p_0 = 1/2`: (a) `p' > 0` on `(1, oo)`, (b) analytic, (c) `int x^(y-1) e^(-yx) dx < oo`. By (2.3.16),
`q/p' = q dx/dv = g(w)/(2w)` with `v = w^2`, so `b_s = g_s/2`. On `(0, 1]` (reflect `x - 1 -> 1 - x`) the same holds with `b_s = (-1)^s g_s/2`.
Adding: `J(y) ~ e^(-y) Sum_k Gamma(k + 1/2) g_(2k) y^(-k-1/2) = e^(-y) (pi/y)^(1/2) Sum_k g_(2k) (2k-1)!! (2y)^(-k)`.
The right side is the full formal SP of `e^(-(1 + w^2)/z) g(w) dw` at `w = 0` (Hessian `-2`, `<w^(2k)>_(1/2) = (2k-1)!! 2^-k`), which by L1a
(the formal map `w -> X(w)`, `X'(0) = sqrt 2 > 0`, positive roots) is the full formal SP `e^(-y) (2 pi/y)^(1/2) N_0(1/y)` of
`e^((log x - x)/z) x^-1 dx` at `x = 1`. On the other hand DLMF 5.11.1 with the permitted operations gives
`J(y) e^y (2 pi/y)^(-1/2) = exp(Ln Gamma(y) - (y - 1/2) ln y + y - (1/2) ln 2pi) ~ exp(Sum_k b_k y^(1-2k))`. By uniqueness,
`N_0(1/y) = exp(Sum_k b_k y^(1-2k))`, i.e. `N_0 = St_0`. QED
This route uses the analytic inputs DLMF 2.3.15 and 5.11.1; Proof A (s.9.3) uses none and in addition proves (iv).

### 9.5 Comparison with gamma's single-Gamma lemma (GAMMA_PROOFS s.6)
* gamma's lemma: the formal saddle series of the inverse Mellin (Cahen-Mellin) integrand `Gamma_St(eta) Y^-eta eta^r` at `eta = Y` is
  `(-Y d/dY)^r e^(-Y)`. Its inputs are `Gamma_St(eta+1) = eta Gamma_St(eta)` (assumed there as classical; proved here as (iv)) and the
  translation invariance (S4) (proved there by a valuation argument; here it is Lemma T(b) of s.3.2, with the weight bound explicit).
* L3 is the dual statement: `Gamma_St` is the OUTPUT of the formal SP of Euler's x-integral, not the input of an eta-integral. Its
  proof uses FL1 in `x` (the analogue of the functional equation) and the shift `a -> a + z` (L1c(iv), the analogue of (S4)).
* Neither lemma implies the other by itself. The chart argument needs L3 (the Euler integrals of steps (a), (d), (e)); gamma's lemma
  belongs to gamma's overlap-ring route, which the chart argument replaces.

## 10. The formal residue: FL2 and FL3 (PROVED)

### 10.1 The functional
Setting of Theorem P(b) (POLY s.6 steps 2-3; the poly referee's technical remark, REFEREE s.3): `K` a field of characteristic 0
(`Q` or `Q(zeta)`), the residue variable `X`, the parameter `u`, `v := 1/w`. Base phase `f_0(X) = P_+(X) + P_-(1/X)`, `P_+-` polynomials
without constant term (here `f_0 = 1/X + X^4/4`); `d := max(deg P_+, deg P_-)`.
**The ring L.** `L := { Sum_(j in Z) c_j X^j : c_j in K((v)), ord_v(c_j) >= eps|j| - M for some eps > 0, M }`.
**Lemma 10.1 (PROVED).** `L` is a `K((v))`-algebra containing `K((v))[X, 1/X]` and `e^(-v f_0)`; `d/dX` is a derivation of `L`.
*Proof.* Product: the coefficient of `X^j` is `Sum_(i) c'_i c''_(j-i)`; `ord_v >= eps'|i| + eps''|j - i| - M' - M''`, which tends to infinity with
`|i|` and is `>= min(eps', eps'')|j| - M' - M''`, so the sum converges in `K((v))` and the bound is of the required form.
`e^(-v P_-(1/X)) = Sum_k (-v P_-(1/X))^k/k!`: the coefficient of `X^-j` has `ord_v >= j/d`; likewise `e^(-v P_+(X))`; multiply. QED
**Definition.** `Res : L -> K((v))`, `Sum c_j X^j -> c_(-1)`; on `L[[u]]` coefficientwise in `u`. An integrand of Theorem P(b) is, termwise
in `u`, `e^(-v f_0) A(X)` with `A` a Laurent polynomial in `X` with coefficients in `K[w, 1/w]`: an element of `L[[u]]`.
*Analytic description (the referee's definition).* For `K = Q`, complex `w != 0` and a Laurent polynomial `A`,
`(1/2 pi i) oint_(|X|=1) e^(-f_0(X)/w) A(X) dX` equals `Res[e^(-v f_0) A]` evaluated at `v = 1/w` (the Laurent expansion on `C^*` of the
holomorphic integrand; each `c_j(v)` is an entire series in `v`). Not needed for the proofs.

### 10.2 FL2: the residue kills exact forms
For every `Om in L`: `Res(d Om/dX) = 0`, since `d/dX Sum c_j X^j = Sum j c_j X^(j-1)` has `X^-1`-coefficient `0 . c_0`. In particular, termwise
in `u`, for `f = f_0 + f_1`, `f_1 in u K[X, 1/X][[u]]`, and `Q in K[w, 1/w][X, 1/X][[u]]`: `e^(-f/w) Q = e^(-v f_0) (e^(-v f_1) Q) in L[[u]]` and
`Res[ e^(-f/w) (Q' - (f'/w) Q) ] = Res[ d/dX (e^(-f/w) Q) ] = 0`. This covers the exact forms `d(Y^m e^(-f/w))` of Theorem P(b), step 2. QED

### 10.3 FL3: invariance under X -> X(1 + O(u)), termwise in u
**Lemma.** Let `eps(X; u) in u K[X, 1/X][[u]]` (e.g. the CFU map `Y = X + eps = X(1 + Sum_k u^k Y_k(X^5))`). For `Om in L[[u]]` define the
pull-back `Om(X + eps) := Sum_(j>=0) Om^(j)(X) eps^j/j!` (u-adically convergent). Then `Res_X[ Om(X + eps)(1 + eps') ] = Res[Om]`.
*Proof.* `Om(X + eps)(1 + eps') = Om(X) + d/dX [ Sum_(j>=1) Om^(j-1)(X) eps^j/j! ]`: expanding the derivative gives
`Sum_(j>=1) Om^(j) eps^j/j! + Sum_(j>=1) Om^(j-1) eps^(j-1) eps'/(j-1)!`, which is the left side minus `Om`. Apply FL2 termwise in `u`. QED
*Consistency with the naive pull-back.* `Om -> Om(X + eps)` is a continuous ring homomorphism (Leibniz for higher derivatives), and for
Laurent monomials it is substitution (`(X + eps)^-1 = Sum_j (-1)^j X^(-1-j) eps^j`). So for `Om = e^(-f(Y)/w) beta(Y)` it is
`e^(-f(X+eps)/w) beta(X+eps)` with `f(X+eps)`, `beta(X+eps)` expanded termwise in `u`: the pull-back used in POLY s.6 step 1.

### 10.4 The residue in the original variables (SKETCH; outside FL2-FL3)
Theorem P(b) step 6 and the referee's s.3 also use `V = Res_(x=0)(e^(-W/z) at dx) = s Res_X(e^(-What/w) Bhat dX)` (`x = sX`, `z = s^4 w`,
`s^5 = u`). Here `s = u^(1/5)` is not in the coefficient ring, so this is not an instance of FL3. In the original variables the
`u^N`-coefficient of `V` is a finite sum (`Lam = -log(1-u) = O(u)`, `G = O(x^4)`); the monomial `u^N z^r x^al dx` becomes
`s^(5N + 4r + al + 1) w^r X^al dX`, and both organisations are expansions of the same product of the same elementary factors. Both
residues are then the coefficient of one monomial of one element of the ring `B` of formal sums `Sum c s^nu w^r X^al` with `nu >= 0`,
`r <= r_0(nu)` and finitely many `al` per `(nu, r)` (a ring by the argument of Lemma 10.1). What is not written out: that each
exponential identity used to regroup the factors (e.g. `e^(-Lam/(xz)) = e^(-v/X) e^(-v(Lam/u - 1)/X)`) holds in `B`; this is routine.
The poly referee checked the identity numerically through `u^30` in two frames (`rf_res`).
