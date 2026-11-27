# GAPS_G2 part 2: the commutant and the projection (G1, G2, step 4), the lifting and scalar lemmas (G3)

Continuation of `GAPS_G2.md` (setting and local notation in its s.2).

## 3. `lem:gau-comm` (items 1a-1d)

**Statement (proposed final form).** Let `F` be a field of characteristic 0 (in the application `F = KK = Q((u))`, or `k((u))`),
`b~_0, .., b~_4` in `F^x`, `A~ = Sum_a b~_a E_(a,a+1)`, `calB_a^(m)`, `Bprod = calB_a^(5)` as in s.2, `gl^(m) = (+)_a F E_(a,a+m)`.
* (a) `A~^m = Sum_a calB_a^(m) E_(a,a+m)` for `0 <= m <= 5`, so `A~^5 = Bprod 1`; `A~` is invertible, `calY := A~^-1 = Sum_a
  b~_a^-1 E_(a+1,a)`, `calY^j = Sum_a (calB_a^(j))^-1 E_(a+j,a)` for `0 <= j <= 4`, `calY^5 = Bprod^-1 1`; `det(x 1 - A~) = x^5 - Bprod`.
* (b) `ad A~ = [A~, .]` maps `gl^(m)` into `gl^(m+1)`: for `X = Sum_a xi_a E_(a,a+m)`, `[A~, X] = Sum_a (b~_a xi_(a+1) - xi_a b~_(a+m))
  E_(a,a+m+1)`.
* (c) (commutant) `ker(ad A~) cap gl^(m) = F A~^m`; hence `ker(ad A~) = F[A~] = F[calY] = (+)_(j=0..4) F calY^j`, of dimension 5.
* (d) (image) For `X = Sum_a xi_a E_(a,a+m)` put `lambda_m(X) := Sum_a xi_a/calB_a^(m)` (for `m = 0` read `calB^(5) = Bprod`, i.e.
  `lambda_0(X) = tr(X)/Bprod`, and read `A~^m` as `A~^5 = Bprod 1`, so that `lambda_m(A~^m) = 5` for every `m`). Then
  `im(ad A~) cap gl^(m) = ad A~(gl^(m-1)) = ker lambda_m`, `gl^(m) = F A~^m (+) ker lambda_m`, and `gl_5(F) = F[calY] (+) im(ad A~)`.
* (e) (projection) The projection `pr` onto `F[calY]` along `im(ad A~)` is `F`-linear, preserves the shift, and on `gl^(m)`
  `pr(X) = (1/5) lambda_m(X) A~^m = (1/5) [Sum_a xi_a calB^(5-m)_(a+m)] calY^(5-m)` (`m = 1..4`), `pr(X) = (1/5) tr(X) 1` (`m = 0`).
  Equivalently: `pr(X)` is the unique element of `F[calY]` with `tr(pr(X) Z) = tr(X Z)` for every `Z` in `F[calY]` (trace form).
* (f) (integrality; `F = k((u))` and (H) of s.2) `calY` lies in `gl_5(k[[u]])`. An element `c = Sum_(j=0..4) c_j calY^j` (`c_j` in
  `k((u))`) lies in `gl_5(k[[u]])` iff every `c_j` lies in `k[[u]]`; so the commutant of `A~` in `gl_5(k[[u]])` is
  `k[[u]][calY] := (+)_(j=0..4) k[[u]] calY^j`, a ring. For `X` in `gl_5(k[[u]])` the coefficient of `calY^j` in `pr(X)` lies in
  `u^-1 k[[u]]` for `j = 1..4` and in `k[[u]]` for `j = 0`.
* (g) (remark, not used) `x^5 - Bprod` is separable, so `A~` is regular semisimple over `F` (five distinct eigenvalues over an
  algebraic closure); under (H) it is irreducible over `k((u))`, so `F[calY]` is a field, totally ramified of degree 5 (UR s.3 Rmk (b)).

**Proof.** (a) Induction: `A~^(m+1) = (Sum_a calB_a^(m) E_(a,a+m))(Sum_c b~_c E_(c,c+1)) = Sum_a calB_a^(m) b~_(a+m) E_(a,a+m+1)` and
`calB_a^(m) b~_(a+m) = calB_a^(m+1)`. For `m = 5`, `calB_a^(5) = Bprod` for every `a` and `E_(a,a+5) = E_(aa)`. The product
`A~ (Sum_c b~_c^-1 E_(c+1,c)) = Sum_a b~_a b~_a^-1 E_(aa) = 1` gives `calY`, and the formula for `calY^j` follows as for `A~^m`. In the
Leibniz expansion of `det(x 1 - A~)` a permutation `pi` contributes only if `pi(a)` is in `{a, a+1}` for every `a`; if `pi(a) = a + 1`
for one `a`, then `pi(a+1) != a+1`, so `pi(a+1) = a+2`, and so on: `pi` is the 5-cycle `a -> a+1` (an even permutation), contributing
`prod_a (-b~_a) = -Bprod`; otherwise `pi = id`, contributing `x^5`.
(b) `A~ X = Sum_a b~_a xi_(a+1) E_(a,a+1+m)` and `X A~ = Sum_a xi_a b~_(a+m) E_(a,a+m+1)`.
(c) By (b), `[A~, Sum_m X^(m)] = Sum_m [A~, X^(m)]` with `[A~, X^(m)]` in `gl^(m+1)`; the sum `(+)_m gl^(m)` is direct, so `X` commutes
with `A~` iff every shift component does. For `X` in `gl^(m)` this means `xi_(a+1) = xi_a b~_(a+m)/b~_a` for all `a` in Z/5: a solution
is determined by `xi_0` (the `b~` are invertible), and `xi_a = calB_a^(m)` is a non-zero solution, since `calB_(a+1)^(m) b~_a =
calB_a^(m) b~_(a+m)`. So `ker(ad A~) cap gl^(m) = F A~^m` and `ker(ad A~) = (+)_m F A~^m = F[A~]`. Finally `A~^m = Bprod calY^(5-m)`
(`m = 1..4`, by (a)) and `calY = Bprod^-1 A~^4`, so `F[A~] = F[calY] = (+)_(j=0..4) F calY^j`.
(d) For `Z` in `ker(ad A~)` and any `Y`: `tr(Z [A~, Y]) = tr(Z A~ Y) - tr(Z Y A~) = tr(A~ Z Y) - tr(A~ Z Y) = 0` (use `Z A~ = A~ Z` and
`tr(Z Y A~) = tr(A~ Z Y)`). Take `Z = A~^(5-m)` (`Z = 1` if `m = 0`) and `Y` in `gl^(m-1)`. For `X` in `gl^(m)`, `X A~^(5-m) =
Sum_a xi_a calB^(5-m)_(a+m) E_(aa)` by (a), and `calB_a^(m) calB^(5-m)_(a+m) = calB_a^(5) = Bprod`, so `tr(X A~^(5-m)) = Bprod
lambda_m(X)` (for `m = 0`: `tr(X) = Bprod lambda_0(X)` by definition). Hence `lambda_m` vanishes on `ad A~(gl^(m-1))`. By (c) this
image has dimension `5 - 1 = 4 = dim ker lambda_m` (`lambda_m != 0`), so `ad A~(gl^(m-1)) = ker lambda_m`. Since `lambda_m(A~^m) =
Sum_a calB_a^(m)/calB_a^(m) = 5 != 0` (characteristic 0; for `m = 0` with `A~^5 = Bprod 1`: `Sum_a Bprod/Bprod`), `gl^(m) = F A~^m
(+) ker lambda_m`. By (b) `im(ad A~) = (+)_m ad A~(gl^(m-1))`, whose
intersection with `gl^(m)` is `ad A~(gl^(m-1))`. Summing over `m`: `gl_5(F) = F[A~] (+) im(ad A~)`.
(e) By (d), `X = (lambda_m(X)/5) A~^m + (X - (lambda_m(X)/5) A~^m)` with the second term in `ker lambda_m = im(ad A~) cap gl^(m)`.
The `calY`-form uses `A~^m = Bprod calY^(5-m)` and `Bprod/calB_a^(m) = calB^(5-m)_(a+m)`. Trace form: `X - pr(X)` lies in `im(ad A~)`,
which is orthogonal to `F[calY] = ker(ad A~)` by the first line of (d); uniqueness because the trace form on `F[A~]` is non-degenerate:
`A~^n` has shift `n`, so `tr(A~^n) = 0` unless `5 | n`, `tr(A~^0) = 5`, `tr(A~^5) = 5 Bprod`; in the basis `1, A~, .., A~^4` the Gram
matrix has `5` at `(0,0)`, `5 Bprod` at `(i, 5-i)` (`i = 1..4`), zeros elsewhere, determinant `5^5 Bprod^4 != 0`.
(f) The entries of `calY` are `b~_a^-1`, in `k[[u]]` (for `a = 1`: `b~_1^-1 = u (u b~_1)^-1`). By (a) `calY^j` is supported on the five
positions `(a+j, a)`, those of shift `-j`; for `j = 0..4` these sets partition the 25 positions, so `(Sum_j c_j calY^j)_(a+j,a) =
c_j (calB_a^(j))^-1`. For `j <= 4` the window `{a, .., a+j-1}` has at most 4 elements and contains the pole index 1 at most once, so
`calB_a^(j)` is a unit of `k[[u]]` or `u^-1` times a unit, and `(calB_a^(j))^-1` lies in `k[[u]]`; for `a = 2` the window
`{2, .., j+1}` lies in `{2, 3, 4, 0}`, so `(calB_2^(j))^-1` is a unit. Hence `c` is integral iff every `c_j` is. `k[[u]][calY]` is a
ring because `calY^5 = Bprod^-1 = u beta^-1` is integral; it is the integral commutant by (c). For `X` in `gl^(m) cap gl_5(k[[u]])`,
`m = 1..4`, the `calY^(5-m)`-coefficient `(1/5) Sum_a xi_a calB^(5-m)_(a+m)` is a sum of integral `xi_a` times products of `5 - m <= 4`
consecutive weights (at most one of them `b~_1`), so it lies in `u^-1 k[[u]]`; for `m = 0` it is `(1/5) tr X`, integral.
(g) `gcd(x^5 - Bprod, 5x^4) = 1` as `Bprod != 0`. With `Bprod = beta/u`, `x = 1/X`: `x^5 - Bprod = -(beta/u) x^5 (X^5 - u beta^-1)`,
and `X^5 - u beta^-1` is Eisenstein at the prime `u` of `k[[u]]`. QED

*Source comparison.* UR s.3 (ii) states (c), (d), (e) and the pole bound of (f) without proof ("regular semisimple over K, so ...";
"Then pr(X) = ..."); UR's "for each `j <= 4` some entry of `y^j` is a unit" is (f) with `a = 2`. Nothing else is needed: the proof
above uses only that the five weights are invertible and that the characteristic is 0 (`lambda_m(A~^m) = 5`).

## 4. The lifting lemma and the scalar lemma

### 4.1 `lem:gau-lift` (UR "Lemma 1"; [W] in UR s.3 (iii), re-derived; stated here because s.5 uses its kernel clause)
Assume (H). Let `X` in `gl_5(k[[u]])` have shift `m + 1` (`0 <= m <= 4`) and `pr(X) = 0`. Then there is `calZ` in `gl_5(k[[u]])` of shift
`m` with `[A~, calZ] = X`; the holomorphic solutions of shift `m` are `calZ + k[[u]] calY^(5-m)` (`m >= 1`) resp. `calZ + k[[u]] 1`
(`m = 0`); every holomorphic solution of `[A~, calZ'] = X` (any shift) is `calZ + c` with `c` in `k[[u]][calY]`. Data over `Q` give
`calZ` over `Q`.
*Proof.* Write `X = Sum_a xi_a E_(a,a+m+1)` and look for `calZ = Sum_a rho_a calB_a^(m) E_(a,a+m)`. By (b),
`[A~, calZ] = Sum_a (rho_(a+1) - rho_a) calB_a^(m+1) E_(a,a+m+1)` (use `b~_a calB_(a+1)^(m) = calB_a^(m) b~_(a+m) = calB_a^(m+1)`). So
the equation is `rho_(a+1) - rho_a = sg_a := xi_a/calB_a^(m+1)` (`a` in Z/5), solvable iff `Sum_a sg_a = lambda_(m+1)(X) = 0`, i.e.
iff `pr(X) = 0` (s.3 (e)). Each `sg_a` is integral, and `sg_a(0) = 0` when the window `{a, .., a+m}` contains 1 (then `calB_a^(m+1)` has
the pole; for `m + 1 = 5` every window contains 1). Put `rho_1 := 0` and `rho_(a+1) := rho_a + sg_a` for `a = 1, 2, 3, 4`. The entry
`rho_a calB_a^(m)` is integral iff `rho_a` is integral and, when the window `{a, .., a+m-1}` contains 1 (i.e. `a` in `{2-m, .., 1}`),
`rho_a(0) = 0`. For `a` in `{2-m, .., 0}` one has `rho_a = rho_1 - Sum_(c=a..0) sg_c` and each window `{c, .., c+m}` (`2-m <= c <= 0`)
contains 1, so `rho_a(0) = rho_1(0) = 0`: `calZ` is integral. The difference of two holomorphic solutions is in `ker(ad A~) cap
gl_5(k[[u]]) = k[[u]][calY]` (s.3 (c), (f)). QED (UR's text, with the choice `rho_1 = 0` made explicit.)

### 4.2 `lem:gau-scal` (UR "Lemma 2"), with the two projection facts (item 1e, G3)
**Statement.** Assume (H), (EP); `Theta_U`, `dlog Theta_U`, `frak d` as in s.2 (the `d^U_a` are units of `k[[u]]`: `b~_a/b_a` is a unit
for `a != 1`, and `b~_1/b_1 = (u b~_1)/(u b_1)`; the cyclic definition closes up because `prod_a b~_a/b_a = 1` by (EP)).
* (a) For `c = Sum_(j=0..4) c_j calY^j`, `c_j` in `k((u))`: `pr(c' + dlog Theta_U c) = Sum_j (c_j' + kappa^U_j c_j) calY^j` with
  `kappa^U_j := -(j/5) Bprod'/Bprod + (1/5) tr(dlog Theta_U)`. In particular
  `pr((calY^j)') = -(j/5) (Bprod'/Bprod) calY^j` and `pr(dlog Theta_U calY^j) = (1/5) tr(dlog Theta_U) calY^j`.
  Moreover `kappa^U_j = j/(5u) + h_j` with `h_j := -(j/5) beta'/beta + (1/5) tr(dlog Theta_U)` in `k[[u]]`.
* (b) For `Z` in `gl_5(k[[u]])` the equation `pr(c' + dlog Theta_U c) = -pr(Z)` has solutions `c` in `k[[u]][calY]`; the components
  `c_1, .., c_4` are unique and `c_0` is unique up to `c_0 -> c_0 + gamma frak d`, `gamma` in `k`.
* (c) For `Z = 0` the solutions in `k[[u]][calY]` are exactly `gamma frak d 1`, `gamma` in `k`.
* (d) If the data are over `Q`, so are the solutions, with `gamma` in `Q` (`frak d` in `1 + u Q[[u]]`).
(`kappa^U_j` is a local letter: `kappa_a` of NOTATION table 4 are the cyclic coefficients; the editor may prefer another letter.)

**Proof.** (a) `(calY^j)' = Sum_a ((calB_a^(j))^-1)' E_(a+j,a)` has shift `5 - j`. In the form `Sum_(a') xi_(a') E_(a', a'+5-j)` one
has `a' = a + j` and `xi_(a') = ((calB_(a'-j)^(j))^-1)'`. By s.3 (e) with `m = 5 - j` (and `a' + 5 - j = a' - j = a` in Z/5),
`pr((calY^j)') = (1/5)[Sum_a ((calB_a^(j))^-1)' calB_a^(j)] calY^j = -(1/5)[Sum_a (calB_a^(j))'/calB_a^(j)] calY^j`. Now
`(calB_a^(j))'/calB_a^(j) = Sum_(i=0..j-1) b~_(a+i)'/b~_(a+i)`, and as `a` runs over Z/5 each index `c` occurs for exactly the `j` pairs
`(a, i) = (c - i, i)`: the sum is `j Sum_c b~_c'/b~_c = j Bprod'/Bprod`. ("Every weight occurs in `j` of the five windows", UR.)
Next `dlog Theta_U calY^j = Sum_a e_(a+j) (calB_a^(j))^-1 E_(a+j,a)`, `e_c := (d^U_c)'/d^U_c`; here `xi_(a') = e_(a') (calB^(j)_(a'-j))^-1`
and s.3 (e) gives `pr = (1/5)[Sum_(a') e_(a')] calY^j = (1/5) tr(dlog Theta_U) calY^j` (for `j = 0` it is the diagonal case `m = 0`).
`pr` is `k((u))`-linear, so `pr(c' + dlog Theta_U c) = Sum_j pr(c_j' calY^j + c_j (calY^j)' + c_j dlog Theta_U calY^j) =
Sum_j (c_j' + kappa^U_j c_j) calY^j`. Finally `Bprod'/Bprod = -1/u + beta'/beta` and `tr(dlog Theta_U) = Sum_c (d^U_c)'/d^U_c` is
integral (units), which gives `kappa^U_j = j/(5u) + h_j`.
(b) By s.3 (f), `pr(Z) = Sum_j phi_j calY^j` with `u phi_j` integral (`j >= 1`) and `phi_0` integral. By (a) the equation splits into
`c_j' + kappa^U_j c_j = -phi_j`. For `j = 1..4` multiply by `u`: `u c_j' + (j/5 + u h_j) c_j = -u phi_j`; at `u^n` this reads
`(n + j/5) c_(j,n) = -[u^n](u phi_j) - Sum_(l=0..n-1) [u^(n-1-l)] h_j c_(j,l)`, and `n + j/5 != 0` for every `n >= 0` because `j/5` is not
an integer: a unique solution in `k[[u]]`. For `j = 0`: `c_0 = frak d (gamma - Int_0^u phi_0/frak d du)` solves it
(`frak d' = -(1/5) tr(dlog Theta_U) frak d`), and the difference `e` of two solutions satisfies `(e/frak d)' = 0`, so `e` is in `k frak d`.
(c) is (b) with `phi = 0`. (d) Every pivot is a non-zero rational and `frak d = exp(-(1/5) Int_0^u tr(dlog Theta_U) du)` has rational
coefficients (exponential of a rational series without constant term). QED

*Remark (why the two facts are needed).* `pr` is `k((u))`-linear but does not commute with `d/du`; the facts compute `pr` on the two
terms where `d/du` and the diagonal matrix `dlog Theta_U` meet the non-scalar elements `calY^j`. UR s.3 (iv) states both facts; the
first with the one-line count, the second without derivation (G3). Evidence: s.10, `g2_lemma2.log`.
