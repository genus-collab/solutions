# GAPS_G2 part 3: Theorem U (existence, the leading-term lemma, the uniqueness clause G4, the field G7)

Continuation of `GAPS_G2_2.md`. Setting and local notation: `GAPS_G2.md` s.2.

## 5. `thm:gau-U`

### 5.1 Statement (proposed final form)
Let `k` be a field of characteristic 0 (`k = Q` in the note) and `A`, `A~` cyclic systems with one simple pole over `k[[u]]` ((H) of
s.2) with equal products ((EP)). Let `V(A, A~)` be the set of all solutions `P` in `gl_5(k[[u]])[[z]]` of `z P' = A P - P A~`.
* (a) (existence) There is `P^nat` in `V(A, A~)` with `P^nat_0 = frak d Theta_U`, diagonal with unit entries.
* (b) (uniqueness) The map `k[[z]] -> V(A, A~)`, `c -> P^nat c(z)`, is bijective. Consequently every solution `P` of the gauge
  equation with coefficients in `gl_5(k[[u]])` is `P^nat c(z)`, its `z^0`-term is `P_0 = c(0) frak d Theta_U` (diagonal), and `P` is a
  regular gauge iff `c(0) != 0`. If `P` is one regular gauge, the regular gauges are exactly `P c(z)`, `c` in `k[[z]]^x`.
* (c) (graded gauges) Graded regular gauges exist; if `P` is one, the graded regular gauges are exactly `P c(z)`, `c` in `k[[z^5]]^x`
  (s.6.2).
* (d) (field) If all `b_a`, `b~_a` lie in `Q((u))`, then `P^nat` can be chosen in `gl_5(Q[[u]])[[z]]`; for every field `k` containing
  `Q`, `V_k(A, A~) = P^nat k[[z]]`, and `P^nat c` has rational coefficients iff `c` is in `Q[[z]]`.
* (e) (remark: necessity of (EP); `prop:gau-equalprod`) If `GE(A, A~)` has a solution with `P_0 != 0`, then `prod b_a = prod b~_a`:
  order `z^0` gives `A P_0 = P_0 A~`, hence `A^5 P_0 = P_0 A~^5`, i.e. `(prod b_a - prod b~_a) P_0 = 0` by s.3 (a).

The statement of UR s.3 is (a), (c) and the last sentence of (b) over `C[[u]]`; GI s.1 quotes it over `Q[[u]]`. The stronger form
(b) costs nothing extra and is what `prop:gau-symp` step 3 actually uses (s.5.4, consequence 2).

### 5.2 Reduction (UR s.3 (i), [W]; re-derived)
`Theta_U` (s.2, s.4.2) has unit entries and `A Theta_U = Theta_U A~`, because `(A Theta_U)_(a,a+1) = b_a d^U_(a+1) = d^U_a b~_a =
(Theta_U A~)_(a,a+1)`. Write `P = Theta_U Phat` (`Phat` in `gl_5(k[[u]])[[z]]` iff `P` is). Then `z P' = z Theta_U' Phat + z Theta_U
Phat'` and `A P - P A~ = Theta_U (A~ Phat - Phat A~)`, so
    `GE(A, A~)`  <=>  `z Phat' = [A~, Phat] - z dlog Theta_U Phat`  <=>  `[A~, Phat_0] = 0`  and  `[A~, Phat_(k+1)] = Phat_k' + dlog Theta_U Phat_k`  (`k >= 0`).

### 5.3 Existence (UR s.3 (v), [W]; written out with the order of choices)
Construct `Phat_k = calZ_k + c_k` (`calZ_k` in `gl_5(k[[u]])`, `c_k` in `k[[u]][calY]`) such that for every `k >= 0`
(I_k) `[A~, Phat_k] = Phat_(k-1)' + dlog Theta_U Phat_(k-1)` (with `Phat_(-1) := 0`), and (S_k) `pr(Phat_k' + dlog Theta_U Phat_k) = 0`.
* `k = 0`: `calZ_0 := 0`, `c_0 := frak d 1`. (I_0) holds since `c_0` commutes with `A~`; (S_0) is `lem:gau-scal` (c) with `gamma = 1`.
* `k -> k+1`: `X_(k+1) := Phat_k' + dlog Theta_U Phat_k` is integral and `pr(X_(k+1)) = 0` by (S_k). `lem:gau-lift`, applied to each shift
  component of `X_(k+1)` (`pr` preserves the shift), gives an integral `calZ_(k+1)` with `[A~, calZ_(k+1)] = X_(k+1)`. Choose
  `c_(k+1)` in `k[[u]][calY]` by `lem:gau-scal` (b) with `Z := calZ_(k+1)' + dlog Theta_U calZ_(k+1)` (free constant: any `gamma`, e.g.
  0). Then (I_(k+1)) holds because `[A~, c_(k+1)] = 0`, and (S_(k+1)) holds because `pr(Phat_(k+1)' + dlog Theta_U Phat_(k+1)) =
  pr(Z) + pr(c_(k+1)' + dlog Theta_U c_(k+1)) = 0`.
`P^nat := Theta_U Sum_k z^k Phat_k` solves `GE(A, A~)` by 5.2, and `P^nat_0 = frak d Theta_U`. The only choices are the constant
`gamma` of `lem:gau-scal` (b) at each order (and the normalisation `rho_1 = 0` inside `lem:gau-lift`, which is absorbed by `c_(k+1)`).

### 5.4 The leading-term lemma and the uniqueness clause (item 1f, G4)
**`lem:gau-lead` (new; proposed label).** Assume (H), (EP). Let `P` be in `V(A, A~)` and `n >= 0` with `P_0 = .. = P_(n-1) = 0`. Then
`P_n = gamma frak d Theta_U` for a unique `gamma` in `k`, namely `gamma = (P_n)_(00)(0)` (as `d^U_0 = 1`, `frak d(0) = 1`).
*Proof.* Put `Phat = Theta_U^-1 P`; then `Phat_l = 0` for `l < n`. The equation of 5.2 at order `z^n` reads `[A~, Phat_n] =
Phat_(n-1)' + dlog Theta_U Phat_(n-1) = 0` (for `n = 0` it is the order-`z^0` equation). So `Phat_n` lies in `ker(ad A~) cap
gl_5(k[[u]]) = k[[u]][calY]` (s.3 (c), (f)). At order `z^(n+1)`: `Phat_n' + dlog Theta_U Phat_n = [A~, Phat_(n+1)]` lies in
`im(ad A~)`, so `pr(Phat_n' + dlog Theta_U Phat_n) = 0`. By `lem:gau-scal` (c), `Phat_n = gamma frak d 1`. QED

**Proof of 5.1 (b).** *Injective:* if `c != 0` has order `n` in `z`, then `[z^n](P^nat c) = c_n P^nat_0 = c_n frak d Theta_U != 0`.
*Surjective:* let `P` be in `V(A, A~)`. Define `c_0, c_1, ..` in `k` recursively. Suppose `R^(n) := P - P^nat (c_0 + c_1 z + .. +
c_(n-1) z^(n-1))` is `0 mod z^n` (true for `n = 0`); `R^(n)` lies in `V(A, A~)` (a `k[[z]]`-module). By `lem:gau-lead`,
`[z^n] R^(n) = gamma_n frak d Theta_U`; put `c_n := gamma_n`. Then `R^(n+1) = R^(n) - c_n z^n P^nat` lies in `V(A, A~)` and
`[z^n] R^(n+1) = gamma_n frak d Theta_U - c_n P^nat_0 = 0`, so `R^(n+1) = 0 mod z^(n+1)`. With `c := Sum_n c_n z^n`, for every `l` and
every `n > l`: `[z^l](P - P^nat c) = [z^l] R^(n) - [z^l] P^nat (c - Sum_(m<n) c_m z^m) = 0 - 0`. Hence `P = P^nat c`.
The `z^0`-term is `P_0 = c_0 P^nat_0 = c(0) frak d Theta_U`, invertible iff `c(0) != 0`. If `P = P^nat c_P` is a regular gauge, then
every regular gauge `P' = P^nat c'` equals `P (c_P^-1 c')` with `c_P^-1 c'` in `k[[z]]^x`, and conversely. QED

*Reading.* UR's sentence "the only freedom is one constant per order, which is `P -> P c(z)`" is the surjectivity: the constant
`gamma` chosen at order `k` in 5.3 is the coefficient `c_k`, and every solution arises from the construction with some choice of the
constants. The point that UR does not write is `lem:gau-lead`: a solution that vanishes to order `n` in `z` has a leading term
that is forced (by the order-`z^(n+1)` solvability condition and the scalar lemma) to be a constant multiple of `P^nat_0`.

*Consequences used elsewhere.*
1. Every regular gauge has `P_0 = c(0) frak d Theta_U`: "`P_0` diagonal with unit entries" holds for EVERY regular gauge, not only
   for the constructed one (sharpening of the statement in `def:gau-onepole`/`thm:gau-U`).
2. `prop:gau-symp` step 3 (GI s.4; G4 dependency): for the pair `(A~, A~)` one has `Theta_U = 1`, `frak d = 1`, and `1` is in
   `V(A~, A~)`, so `V(A~, A~) = k[[z]] 1`: every solution `T` of `z T' = A~ T - T A~` with integral coefficients is `t(z) 1`. (Hypotheses
   (H), (EP) hold trivially for the pair. Invertibility of `T_0` is not needed.)
3. `thm:gau-G` Cor. (a) (s.8) and `prop:gau-graded` (s.6) use (b).

### 5.5 The field (item 1h, G7)
UR s.3 works over `O = C[[u]]` ("or convergent series"); GI quotes the theorem over `Q[[u]]` without comment. Proof of 5.1 (d):
* Every object of s.3-s.5 is built from the `b_a`, `b~_a` by field operations of `Q((u))`, `d/du`, and integration of power series
  without `u^-1` term: `Theta_U` (ratios of the weights), `dlog Theta_U`, `calY`, `calB_a^(m)`, `Bprod`, `beta = u Bprod`; `pr`
  (s.3 (e): the factor `1/5` and the `calB`); `frak d = exp(-(1/5) Int_0^u tr dlog Theta_U)` has rational coefficients; the lift of
  `lem:gau-lift` with `rho_1 = 0` (sums of products of the given series); the solutions of `lem:gau-scal` (b) (rational pivots
  `n + j/5`, the rational constant `gamma`). So with every `gamma` chosen in `Q` the construction of 5.3 runs inside `gl_5(Q[[u]])`.
* For `k` containing `Q` the proofs of s.3-s.5 hold verbatim over `k`, so `V_k(A, A~) = P^nat k[[z]]` by (b). If `P = P^nat c` has
  rational coefficients, the recursion of the surjectivity proof computes `c_n = gamma_n = ([z^n] R^(n))_(00)(0)` from rational data,
  so `c` is in `Q[[z]]`.
* Convergence is never used: all statements of s.3-s.5 are identities of formal series. (The analytic statement over convergent series
  is not needed by the chain; it is not claimed.)

**Field of the application.** In the programme family the coefficients of `A` and `A~` lie in `Q[[u]]` (s.2, last paragraph); hence
`P_can` (s.6.3) has coefficients in `Q[[u]]`, and `P_11(0, z)` of a rational graded gauge lies in `Q[[z^5]]`. This settles D15
("Theorem U and `P_can` over `Q`") on the gauge side. `Q(zeta)` enters only through the colour-`i` WKB columns (s.9).

### 5.6 Remark on what the proof uses (for `rmk:gau-Uremarks` (a))
Only (H) and (EP): invertibility of the five weights with the stated orders at `u = 0`, and characteristic 0 (in `lambda_m(A~^m) = 5`,
in the pivots `n + j/5 != 0`, in `frak d`). The pole of `b~_1` enters exactly through the integrality statements: s.3 (f) and the
window argument of `lem:gau-lift`. This is UR's "versality" remark (a), now with its proof.
