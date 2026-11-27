# GAPS_G4, part 4: item 3, `thm:mod-Lreg` with `lem:mod-abdks`, route 2 `prop:mod-period`, `lem:mod-vanper` (W5) -- PROVED mod (Q)

Conventions: GAPS_G4.md s.1. Sources: XY s.0-s.3, XYR s.2-s.3, DC1 TR-7..TR-9, DC2 `thm:mlk-lreg` steps 1-5 (step 3 [W-], step 5's transfer
[M]). ABDKS was read in the full text (current arXiv PDF, s.1.1 and s.7-7.1, Thm 7.1 with (7.5)-(7.12), Remarks 7.2, 7.3). Status: every
step below is written; the only input beyond the definitions is ABDKS Thm 7.1 (already in (Q)) and EO Thm 4.2 (poles only at branch points).

## 1. The quoted theorem (pinned; read in the full text)

ABDKS s.7.1, **Theorem 7.1**: for an arbitrary rational function `x(z)`, all `omega^(g)_(m,0)` of topological recursion with the spectral
curve `(x = x(z), y = z)` and initial form `omega^(0)_(1,0) = -y dx` are given by the closed formula (7.9):
`omega^(g)_(m,0) = [hbar^(2g)] Sum_(r_1..r_m >= 0) prod_i D_(x_i)^(r_i) [w_i^(r_i)] Sum_Gamma hbar^(2 g(Gamma)) prod_((i,j) in E) W_(i,j) prod_(i in V) W_i`,
`Gamma` connected simple graphs on `m` numbered vertices, `D_x omega := d(omega/dx)`, vertex weight (7.6) `W_i = -exp(-w_i (S(w_i hbar
d_(z_i)) - 1) x(z_i)) dz_i/w_i`, edge weight (7.7) `W_(i,j) = exp(hbar^2 w_i w_j S(w_i hbar d_(z_i)) S(w_j hbar d_(z_j)) (z_i - z_j)^-2) - 1`,
`S(u) = sinh(u/2)/(u/2)` ((1.12)). One point (7.10): `Sum_(g >= 1) hbar^(2g) omega^(g)_(1,0) = Sum_(r >= 0) D_x^r [w^r] W_1`. Standing
assumptions (s.1.1): compact `S`, `x`, `y` meromorphic, critical points of `x` simple, `y` regular there with `dy != 0`, `B` the Bergman
kernel; the s.7 introduction restricts to "rational x such that dx has only simple zeros". The kernel (1.1) of ABDKS with `omega^(0)_(1,0)
= -y dx` is `(1/2) Int_z^(sigma z) B/((y(z) - y(sigma z)) dx(z))`, EO's printed kernel; hence `omega^(g)_(m,0) = (-1)^m omega_(g,m)` (engine).
**Hypotheses for `x = W_v`, `y = p` (`v != 0`)** -- the theorem applies to the note's curve DIRECTLY (W5 dissolves): `W_v` is rational;
`dW_v = (p^5 - v) dp/p^2` has exactly five zeros, simple (`p^5 = v`; at `p = 0, oo` it has poles); `y = p` is regular with `dy = dp != 0`
there; `S = P^1`, `B = dp_1 dp_2/(p_1 - p_2)^2` (no A-cycles). Example (7.2) is not needed. Remark 7.3 (limits to the discriminant) is not
used (XYR: `v = 0` is a boundary point, not the discriminant).

## 2. `lem:mod-abdks` (structure; PROVED from Thm 7.1)

**Lemma.** For `2h - 2 + n > 0` there are `K >= 0` and `N_(h,n)` in `Q[v][p_1, .., p_n]` with
`omega_(h,n) = N_(h,n)(p; v) prod_i dp_i / prod_i (p_i^5 - v)^K` (one may take `K = 6h - 4 + 2n`).
*Proof.* (1) *Each term of (7.9) lies in `Q[v][p_i, p_i^-1, (p_i^5 - v)^-1, (p_i - p_j)^-1] prod dp_i`, and for fixed `(h, n)` there are finitely
many.* With `S(u) - 1 = Sum_(k>=1) u^(2k)/(4^k (2k+1)!)`, `-w (S(w hbar d_p) - 1) W_v(p) = Sum_(k>=1) w^(2k+1) hbar^(2k) eps_k(p)`,
`eps_k := -W_v^((2k))/(4^k (2k+1)!)`, and `W_v^((2k)) = (p^4/4)^((2k)) + (2k)! v p^(-2k-1)` lies in `Q[v][p, p^-1]`. Each power `hbar^(2k)` in
`W_i` comes with `w_i^(2k+1)`, each `hbar^2` in `W_(i,j)` with `w_i w_j` and a derivative of `(p_i - p_j)^-2` (here `z_i = p_i`), and each
`hbar` from an `S`-operator with one `w`; so at fixed `hbar^(2h)` the
`w_i`-degrees are bounded, only finitely many `r_i` and graphs contribute, and the coefficients are polynomials in the `eps_k(p_i)` and in
`(p_i - p_j)^-1` with rational coefficients. `D_x(f dp) = d(p^2 f/(p^5 - v))` keeps the ring and adds only powers of `(p^5 - v)^-1`.
(2) *Clearing.* Write the finite sum as `omega_(h,n) = P prod dp_i / (prod p_i^a prod (p_i^5 - v)^K' prod_(i<j) (p_i - p_j)^e)`, `P` in
`Q[v][p]`. For `2h - 2 + n > 0`, `omega_(h,n)` has poles only at branch points (EO Thm 4.2; `p = 0` and the diagonals contain none, `v != 0`).
So for every `v != 0`, `P` vanishes to order `>= e` on `p_i = p_j` and to order `>= a` on `p_i = 0`. Divide `P` by `(p_i - p_j)^e`
(resp. `p_i^a`) in the variable `p_i` over the ring `Q[v][p_k : k != i]` (monic divisors, so the division algorithm works over that
ring): the remainder vanishes for all `v != 0` and all `p`, hence is `0`. So these factors cancel in `Q[v][p]`. (If one wants `K = 6h-4+2n`:
the pole order at the branch points is `<= 6h - 4 + 2n`, e.g. from (2-35), `(d/dx)^d xi` having poles of order `2d + 2` and `d <= 3h-3+n`;
`p_i^5 - v` is monic in `p_i`, so the same division argument removes the surplus powers.) QED

## 3. `thm:mod-Lreg` (PROVED mod Thm 7.1)

**Theorem.** For `2h - 2 + n > 0` and `phi_1, .., phi_n` in `Q[v][p]`, `<phi_1 .. phi_n>_h := Sum_(b_1..b_n) Res_(varsigma_(b_1)) ..
Res_(varsigma_(b_n)) prod_l phi_l(p_l) omega_(h,n)` lies in `Q[v]`. For legs in `LL = Q[v][p] + Q(v)[x]` the residue pairing lies in `Q[v]`.
*Proof.* (1) Fix all variables but `p_l` off the branch points. `phi_l(p_l) omega_(h,n)` is a rational 1-form in `p_l` whose poles are the
branch points (from `(p_l^5 - v)^-K`) and possibly `p_l = oo`; it is holomorphic at `p_l = 0` (`phi_l` polynomial, `N` polynomial). By the
residue theorem `Sum_b Res_(varsigma_b) = -Res_(p_l = oo)`. (2) The result is again of the form (polynomial)/`prod_(m != l) (p_m^5 - v)^K`,
because for one variable, with `(p^5 - v)^-K = Sum_(j>=0) C(K+j-1, j) v^j p^(-5K-5j)` (`|p|^5 > |v|`),
  `Res_(p = oo) p^a dp/(p^5 - v)^K = -C(K+j-1, j) v^j` if `a + 1 = 5(K + j)` with `j >= 0` an integer, and `0` otherwise,
an element of `Q[v]`; iterated residues commute (no diagonal poles). So `<phi_1 .. phi_n>_h = (-1)^n Res_(p_1 = oo) .. Res_(p_n = oo)
[prod_l phi_l N_(h,n) prod dp_l/prod (p_l^5 - v)^K]` is a finite `Q[v]`-combination: in `Q[v]`. (3) Legs in `Q(v)[x]` pair to zero under the
residue pairing (`lem:mod-xfun`), so for legs in `LL` only the `Q[v][p]` components matter. QED
*Remarks.* (i) The circle form of `def:mod-pair` equals the residue pairing for `Q[v][p]` legs only (G4_3, F2). (ii) By homogeneity
(TRRF s.3.4) a monomial pairing is `C v^e`, `e = 2 - 2h - n + Sum m_l/5`; the theorem says `C = 0` whenever `e < 0`; XYR: XY's "11 forced
zeros" are not evidence. (iii) **W5, if XY's normalisation is wanted** (only to reuse its numbers): with `z = p`, `V = 4v`, `x_XY = z^4 +
V/z = 4 W_v`, `lem:mod-scale` gives `omega_(h,n)[W_v] = 4^(2h-2+n) omega_(h,n)[x_XY]|_(V = 4v)`; legs `Q[v][p] = Q[V][z]`, `Q[v][x] = Q[V][x_XY]`,
`-log p = -log z`, the same branch points (`4 z^5 = V` iff `p^5 = v`), and every pairing scales by the same `4^(2h-2+n)`. E.g. XY s.2.2's
`omega_(1,1) = -5 z^4 (3V^2 + 36 V z^5 + 8 z^10)/(4z^5 - V)^4 dz` becomes `omega_(1,1) = -(5/8) p^4 (p^10 + 18 v p^5 + 6 v^2)/(p^5 - v)^4 dp`.

## 4. Route 2 of the Gap Theorem: `prop:mod-period`, `lem:mod-vanper` (PROVED mod Thm 7.1)

**Proposition.** For `h >= 1`: (a) `omega_(h,1) = d Phi_h` with `Phi_h` rational, holomorphic at `p = 0` and `p = oo`; (b)
`Int_(oo)^0 omega_(h,1) = (B_2h/(2h)) v^(1-2h)` (any path avoiding the branch points), equivalently `<-log p>_h := Int_(frak r) omega_(h,1) =
Sum_b Res_(varsigma_b) (-log p) omega_(h,1) = -(B_2h/(2h)) v^(1-2h)` (`lem:mod-vanper`); (c) for `h >= 2`, `F^TR_h = B_2h v^(2-2h)/(2h(2h-2))`.
*Step 1 (exactness).* By (7.10), `omega^A_(h,1) = [hbar^(2h)] Sum_(r>=1) D_x^r [w^r] W_1` (the `r = 0` term is `-[w^1] e^E dp = 0`, `E :=
Sum_k w^(2k+1) hbar^(2k) eps_k = O(w^3)`). With `delta := (1/x') d/dp`, `D_x^r(a dp) = d(delta^(r-1)(a/x'))`. Hence `omega^A_(h,1) = d Phi^A_h`,
`Phi^A_h := Sum_(r>=1) delta^(r-1)(a_r/x')`, `a_r := -[hbar^(2h)][w^(r+1)] e^E = -Sum (1/j!) prod_(i=1..j) eps_(k_i)` over `j >= 1` and
`k_1 + .. + k_j = h` with `r = 2h - 1 + j` (finite: `1 <= j <= h`); in the engine convention `omega_(h,1) = -d Phi^A_h`.
*Step 2 (`p = oo`).* `deg eps_k <= 4 - 2k`, `1/x' = p^2/(p^5 - v)` has degree `-3`, `delta` lowers the degree by 4: every term has degree
`<= (4j - 2h) - 3 - 4(2h - 2 + j) = 5 - 10h <= -5`. So each term vanishes at `oo`, `Phi^A_h(oo) = 0`, and `Phi^A_h` is holomorphic there.
*Step 3 (`p = 0`, the monomial count).* `eps_k = eps^0_k + (terms of order >= 0)`, `eps^0_k := -v p^(-2k-1)/(4^k (2k+1))`;
`1/x' = -(p^2/v)(1 - p^5/v)^-1`; `delta = delta_0 + (an operator raising the order by >= 6)`, `delta_0 := -(p^2/v) d/dp`. The *pure part*
(all `eps_k`, `1/x'`, `delta` replaced by `eps^0_k`, `-p^2/v`, `delta_0`) of a term: in `xi := v/p` one has `delta_0 = d/dxi`, `-p^2/v =
dp/dxi = -v/xi^2`, `eps^0_k = -xi^(2k+1) v^(-2k)/(4^k (2k+1))`; with `t := w xi`, `epsl := hbar/(2v)` (local letters),
`a^0_r = -(2v)^(-2h) xi^(r+1) c_(h,r)`, `c_(h,r) := [epsl^(2h)][t^(r+1)] exp(E_0)`, `E_0 := -Sum_(k>=1) epsl^(2k) t^(2k+1)/(2k+1)`, and the pure
term is `(d/dxi)^(r-1)[a^0_r (-v/xi^2)] = v (2v)^(-2h) (r-1)! c_(h,r)`, a constant. Each replacement of a pure factor by its correction raises
the order in `p` by `>= 3` (`>= 2k + 1` for `eps_k`, `>= 5` for `1/x'` and `delta`; `delta` never lowers the order). So every term is a
constant plus `O(p^3)`: `Phi^A_h` is holomorphic at 0 and `Phi^A_h(0) = v (2v)^(-2h) S_h`, `S_h := Sum_r (r-1)! c_(h,r)`, independent of
the polynomial part of `x` (the family version is a remark).
*Step 4 (the identity, exact; replaces the "Borel identity").* With `r + 1 = 2h + j`: `c_(h,r) = ((-1)^j/j!) [T^h] phit(T)^j`, `phit(T) :=
Sum_(k>=1) T^k/(2k+1)` (the `epsl^(2h) t^(2h+j)` coefficient of `E_0^j/j!`, after `t^j` is factored out). Hence
  `S_h = Sum_(j=1..h) (-1)^j ((2h-2+j)!/j!) [T^h] phit^j = (2h-2)! [T^h] (1 + phit(T))^(1-2h)`
(binomial series `(1 + phit)^(1-2h) = Sum_j (-1)^j C(2h-2+j, j) phit^j`; the `j = 0` term has no `T^h`, `h >= 1`). With `T = s^2`,
`1 + phit(s^2) = artanh(s)/s`, so `[T^h](1 + phit)^(1-2h) = [s^1] artanh(s)^(1-2h)`; Lagrange-Buermann (`s = tanh w`; for `n != 0`,
`[s^n] H(w(s)) = (1/n) Res_(w=0) H'(w) tanh(w)^-n dw`) gives `[s^1] w^(1-2h) = (1-2h) [w^(2h)] (w coth w) = (1-2h) 4^h B_2h/(2h)!`. So
**`S_h = -4^h B_2h/(2h)`** and `Phi^A_h(0) = -(B_2h/(2h)) v^(1-2h)`. (This is route 1's `A_h` of TRM s.4.1 (v): `S_h = 4^h A_h`; XY's Borel
integral `Int_0^oo e^-u epsl^2/sinh^2(epsl u) du` is a heuristic for the same identity and need not be printed, finding F5.)
*Step 5.* (b): `Int_oo^0 omega_(h,1) = -(Phi^A_h(0) - Phi^A_h(oo)) = (B_2h/2h) v^(1-2h)`, path-independent (`Phi^A_h` rational). Contour lemma:
`Res_(varsigma_b)(log p dPhi) = -Res_(varsigma_b) Phi dp/p` (`log p` single-valued near `varsigma_b`), and the residue theorem for the rational
form `Phi dp/p` (poles at the `varsigma_b`, `0`, `oo`): `Sum_b Res_(varsigma_b) log p . omega_(h,1) = Phi(0) - Phi(oo) = Int_oo^0 omega_(h,1)`
(`Phi := -Phi^A_h`); with `<-log p>_h = Int_(frak r) omega_(h,1) = -Int_oo^0 omega_(h,1)` this is `lem:mod-vanper` (orientation as XYR s.2).
(c): `Phi_ydx = p^5/5 - v log p`, so `(2-2h) F^TR_h = (1/5) <p^5>_h + v <-log p>_h`; `<p^5>_h` lies in `Q[v]` (s.3) and is homogeneous of
exponent `2 - 2h < 0`, hence `0`; so `F^TR_h = -v (B_2h/2h) v^(1-2h)/(2-2h) = B_2h v^(2-2h)/(2h(2h-2))`. QED

## 5. Evidence (VERIFIED; `g4_period.py 5` -> `g4_period_5.log`, exact over Q, 0.2 s; `g4_w11.log` (B))

* `S_h = -4^h B_2h/(2h)` from the finite sum for `h = 1..15`; control (coefficients `1/(2k+2)` in `phit`) agrees for no `h`.
* From (7.10) for `x = A z^4 + V/z`, `A = 1, 1/4`, `V = 1, 3`, `h = 1..5`: every term `delta^(r-1)(a_r/x')` has no negative power at 0,
  the constant terms sum to the pure-pole value `= v(2v)^(-2h) S_h = -(B_2h/2h) V^(1-2h)` (e.g. `-1/12, 1/120, -1/252, 1/240, -1/132` at
  `V = 1`), every term vanishes at `oo` (exponents `> 0` in `1/z`).
* `omega_(1,1)` of the closed form = `omega_(1,1)` of the recursion (local formula of G4_5, summed over the five branch points) at 8 points.
