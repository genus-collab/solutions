# Gap cluster G3, part 3 of 4: item G3-2, `lem:tra-conv` (a)-(d) with proofs, and the added parts (e), (f)

Continues `GAPS_G3.md` (conventions s.2) and `GAPS_G3_2.md` (the frame lemma; cited as "FL Step n"). u-adic setting throughout:
`v = v(u)` in `u Q[[u]]`, `v != 0`; `(e)` assumes moreover `v'(0) != 0` (O-2). Model correlators: for flat insertions,
`<z^(k_1) H^(a_1) .. z^(k_n) H^(a_n)>_h(v)` is the integral over `Mbar_(h,n)` of the model class against `prod psi_i^(k_i)`
(`def:set-modcohft`); it vanishes unless `Sum k_i <= D := 3h - 3 + n` (the class has non-negative degree and `dim Mbar_(h,n) =
3h - 3 + n`), and it is `C v^e`, `e = 2 - 2h + Sum wt(k_i, a_i)`, `wt(k,a) = (4k + nu_a - 4)/5`, `C = 0` unless `e` is in `Z`.

## 1. Proposed statement of `lem:tra-conv`

**Lemma (u-adic convergence).**
* (a) *Frames.* For every family `(c_(k,a))_(k >= 0, a != 1)` in `Q[[u]]` the family `(c_(k,a) frame(k,a))` is summable; the
  `z^j H`-coefficient of `frame(k,a)` has u-order `>= (4(k-j) + nu_a + 1)/5`. The decomposition `q = Sum q_(k,a) frame(k,a) +
  sigma(q)(z) H` of `lem:tra-frames` (c) holds with these sums; the same holds in two slots for edges (`GAPS_G3_2.md` s.3).
* (b) *Frame expansions inside correlators.* For `h >= 0`, `n >= 1`, `2h - 2 + n > 0`, the correlator extends from flat insertions
  to a symmetric `Q[[u]]`-multilinear map `<.>_h : (HH^u_+)^n -> Q((u))`,
  `<X_1 .. X_n>_h := Sum_(Sum k_i <= D) prod_i (X_i)_(k_i, a_i) <z^(k_1) H^(a_1) .. z^(k_n) H^(a_n)>_h(v(u))` (a finite sum),
  of finite type in each slot. If `X_i = Sum_(beta in B) c_(i,beta) t_beta` with `B = {(k,a) : a != 1} u {H}`, `t_(k,a) = frame(k,a)`,
  `t_H = H`, `c_(i,beta)` in `Q[[u]]` (every element of `K^u` is of this form, FL Step 4 and Cor. (ii)), then the family
  `(prod_i c_(i,beta_i) <t_(beta_1) .. t_(beta_n)>_h)_(beta in B^n)` is summable in `Q((u))` and its sum is `<X_1 .. X_n>_h`. Every
  term is a correlator of homogeneous elements of `K` over `Q[v]` evaluated at `v(u)`; no flat component `z^j H`, `j >= 1`, of a
  frame occurs as an insertion of its own. The same holds for the vertex factors with leaves and for `Cont_Gamma` with every
  edge expanded in frames (s.2, (b3)-(b4)).
* (c) *Vanishing series; pairings.* `[z^-n] Y^C_b(v(u), z) = eps_b (-1)^n v(u)^(n-q)/(4^q q! (n-q)!)`, `q = floor(n/5)`, for `n = nu_b + 1
  (mod 5)`, else `0`; its u-order is `>= n - floor(n/5) >= 4n/5`. Hence `Y^C(v(u))`, `J^C`, `z^-k J^C` (`k >= 1`) lie in `HH^u`
  (the last in `HH^u_-`), and every pairing and product of s04 (`Pi_i(q)` for `q` in `HH^u_+`; `P Y^C`, `P J^C`, `R f`; `Omega_s(f, g)`
  for `f, g` in `HH^u`; the contraction of `lem:tra-E`) is a ring operation in `Q((z))[[u]]`, resp. `Q((x))[[u, psi']]` (part 4),
  i.e. a sum that is finite modulo every power of `u`.
* (d) *Leaves.* A flat insertion `z^k H^b` with `k >= 2` has `wt(k,b) >= 3/5`. At a vertex of genus `h` with `n` slots and `m`
  insertions of psi-degree `>= 2`, a non-zero correlator has slot degrees `Sum_i k_i <= 3h - 3 + n - m`; in particular `m <= 3h - 3 + n`.
* (e) *Vertex resummation* (added; for `thm:red-main` steps 8-9). See s.4.
* (f) *Dilaton re-expansion* (added; for `lem:tra-F` step 1 and `def:red-flatsum`). See s.5.

## 2. Proofs of (a)-(c)

**(c) first** (the other parts use it). The closed form is `lem:gau-vanshape` (GI s.2, Lemma V): `Y^C_b = eps_b I^C_(nu_b)`,
`I^C_m = Sum_(i >= 0) (-1/(4z))^i/i! (-v/z)^(m+4i+1)/(m+4i+1)!`; the term with `z^-n` has `n = 5i + m + 1` and `v`-power `m + 4i + 1 =
n - i`, and since `0 <= nu_b + 1 <= 4`, `i = floor(n/5)`; the sign is `(-1)^(i + m + 4i + 1) = (-1)^n`. Since `ord_u v(u) >= 1`, the
u-order is `>= n - floor(n/5) >= 4n/5`, so modulo `u^N` only `n < 5N/4 + 1` occur: `Y^C(v(u), z)` lies in `Q[z^-1][[u]]^5` (as
`lem:gau-vanshape` states), and `J^C = eta^-1 Y^C` and `z^-k J^C` (`k >= 1`) follow. The rest is `GAPS_G3.md` s.2.1: `Q((z))[[u]]`
is a ring containing `Q[[z]][[u]]` (the coefficients of `P`, `P^-1`, `R`, `q_T`) and `Q[z^-1][[u]]` (the vanishing series); products,
`[z^i]`, `( )_(>=0)`, `( )_(<0)`, `z -> -z` and `Res_z` are defined coefficientwise in `u` with finite sums. In particular:
`Pi_i(q) = -[z^i] <Y^C(-z), q(z)>` is defined for every `q` in `HH^u_+` and equals the u-adically convergent sum of FL Step 1;
`(P z^-k J^C)_(>=0) = z^-k (P J^C)_(>=k)` (`cor:tra-edges` step 3) is an identity in `Q((z))[[u]]` (multiplication by `z^-k` commutes with
`P` and shifts the splitting); the identity `Omega_s(f, P^-1 g) = Omega_s(P f, g)` of `lem:tra-T` (a) is associativity in
`gl_5(Q((z))[[u]])`: `f(-z)^T eta P^-1(z) g(z) = f(-z)^T P(-z)^T eta g(z)` by (Sp) in vector form, `eta P(z)^-1 = P(-z)^T eta`, then
`Res_z`; and `T_P = -[z^-1] Y'_0` lies in `Q[[u]]`. The flat-frame hypothesis of `rmk:tra-K` ("the covector series decays") is
the bound above: `phi_j = O(u^(j - floor(j/5)))`.

**(a)** is FL Step 4 (first bullet) and FL s.3; the u-order of the `z^j H`-coefficient is FL Step 3 (exponent `(4(k-j)+nu_a+1)/5`)
with `ord_u v >= 1`.

**(b).** (b1) *Definition, finiteness, finite type.* Only the finitely many flat slots with `Sum k_i <= D` contribute; each flat
correlator is a Laurent monomial `C v^e`, `e >= 2 - 2h - n` (`wt >= -1`), so at `v(u)` it lies in `u^(e_0) Q[[u]]` with `e_0 :=
(2 - 2h - n) ord_u v` if `2 - 2h - n < 0` (else `e_0 = 0`). Hence `<X_1 .. X_n>_h` lies in `u^(e_0) Q[[u]]` for all `X_i` in `HH^u_+`,
it is symmetric and `Q[[u]]`-multilinear, and with `X_2 .. X_n` fixed it is of finite type in `X_1`: modulo `u^N` it depends only on
the `(X_1)_(k,a)`, `k <= D`, taken modulo `u^(N - e_0)`. (b2) *The frame expansion.* Let `beta = (beta_1, .., beta_n)` in `B^n` and put
`kappa(beta_i) := k` if `beta_i = (k, a)`, `kappa(H) := 0`. By (b1) the term `prod_i c_(i,beta_i) <t_(beta_1) .. t_(beta_n)>_h` is the
finite sum over flat slots of `prod_i c_(i,beta_i) (t_(beta_i))_(k_i,a_i) <slot>_h(v(u))`. A frame `t_(beta_i)` with `kappa(beta_i) > D`
enters a slot with `k_i <= D` only through its `H`-coefficient at `z^(k_i)`, of u-order `>= (4(kappa(beta_i) - D) + 1)/5` (by (a)).
So the term has u-order `>= e_0 + Sum_i max(0, (4(kappa(beta_i) - D) + 1)/5)`, and for each `N` only finitely many `beta` (all
`kappa(beta_i) <= D + 5(N - e_0)/4`) give a term outside `u^N Q[[u]]`: the family is summable. By the summation lemma applied in
slot 1 (finite type, (b1)), `<X_1 .. X_n>_h = Sum_(beta_1) c_(1,beta_1) <t_(beta_1), X_2, .., X_n>_h`; then in slot 2 for each
`beta_1`; and so on. This expresses `<X_1 .. X_n>_h` as an iterated sum of the family, which equals its unordered sum
(regrouping, `GAPS_G3.md` s.2.2). Each `t_beta` is a homogeneous element of `K` over `Q[v]` (a frame, or `H`), evaluated at `v(u)`;
the expansion is indexed by `B`, so a flat `z^j H`, `j >= 1`, never appears as an insertion: ML_K(h) applies to each term. QED (b)
(b3) *Vertex factors with leaves.* For `L` in `z^2 HH^u_+` (in `thm:red-main`: `L = d T_(>=2)`, `lem:tra-T` (c)) the vertex factor
`Sum_(m=0)^(D) (1/m!) <X_1 .. X_n L^m>_h` (finite by (d)) is, for each `m`, of finite type in each slot. If moreover `L` lies in `K^u`
(in `thm:red-main` by `lem:tra-T` (c): `L = d q_T + z 1`), write `L = t' + delta H` with `delta := Pi_0(L)` (`= delta_P`) and `t' := L -
delta H` in `K_0^u`, whose frame coordinates `t'_(k,a) = L_(k,a)` vanish for `k <= 1` (FL Step 4); the multinomial expansion
for fixed `m` is finite, and each term is handled by (b2) with `n + m` slots. (b4) *Graph terms.* For a stable graph `Gamma`,
`Cont_Gamma` is a finite sum (`lem:tra-F`): for an edge from a vertex `w` to a vertex `w'` only its `psi^k psi'^l`-coefficients with
`k <= D_w`, `l <= D_(w')` matter (`D_w = 3h_w - 3 + n_w`). As a function of one edge bivector (the others fixed) it is linear and
depends, modulo `u^N`, on finitely many coefficients of that edge modulo a fixed power of `u`. Expanding every edge by
`GAPS_G3_2.md` s.3 gives a family indexed by frame labels of all half-edges; by the bound of (b2) at each vertex it is summable,
and by the summation lemma, edge by edge, its sum is `Cont_Gamma`. This is the precise form of `thm:red-main` steps 1 and 6.

## 3. Proof of (d)

`wt(k,b) = (4k + nu_b - 4)/5 >= (8 - 1 - 4)/5 = 3/5` for `k >= 2`, since `nu_b >= -1`. A correlator with `n` slots of degrees `k_i` and
`m` insertions of degrees `l_j >= 2` at genus `h` vanishes unless `Sum k_i + Sum l_j <= 3h - 3 + n + m` (degree bound with `n + m`
points), so `Sum k_i <= 3h - 3 + n + m - 2m = 3h - 3 + n - m`, and `m <= 3h - 3 + n` because `Sum k_i >= 0`. QED
*Uses.* `lem:tra-F`: the leaf sum at a vertex stops at `m = 3h - 3 + n`, and only slot degrees `<= 3h - 3 + n` occur. In the
hatted second proof (`rmk:tra-K`, MS2 s.7): together with the first statement, each u-order receives finitely many terms. The
source's "each leaf raises the exponent by `>= 3/5`" is the first statement; no proof in the flat frame needs more.

## 4. (e) The vertex resummation (`thm:red-main` (b), steps 8-9)

**Statement.** Let `g >= 2`, `v'(0) != 0`, `delta` in `u^2 Q[[u]]`, `t'` in `K_0^u` with `t'_(k,a) = 0` for `k <= 1`, and `L := t' +
delta H` in `z^2 HH^u_+`. (e1) `<L^m>_g = 0` for `m > 3g - 3`. (e2) The family `(delta^j/(l! j!) <t'^l H^j>_g)_((l,j) in N^2)` is
summable in `Q((u))`, and `Sum_m (1/m!) <L^m>_g = Sum_(l,j) delta^j/(l! j!) <t'^l H^j>_g`. (e3) If ML_K(g) holds, then
`Sum_j (delta^j/j!) <H^j>_g = c_g (v + delta)^(2-2g)` and `Sum_j (delta^j/j!) Sum_(l >= 1) (1/l!) <t'^l H^j>_g` lies in `Q[[u]]`; hence
`Sum_m (1/m!) <L^m>_g = c_g (v + delta)^(2-2g) + Q[[u]]`. (Here `<H^0>_g := F^m_g(v)`, the vertex graph without insertions.)
*Proof.* (e1) by (d) with `n = 0`. (e2) For fixed `m`, `<L^m>_g = Sum_(l+j=m) m!/(l! j!) delta^j <t'^l H^j>_g` (multilinearity and
symmetry). Expand `t'` in frames with `k >= 2`; these have `wt(k,a) >= (8 + 0 - 4)/5 = 4/5` (`a != 1`, so `nu_a >= 0`), and by
homogeneity (`def:set-modcohft`) a correlator of `l` such frames and `j` insertions `H` (`wt(0,1) = -1`) is `C v^e` with
`e >= 2 - 2g + 4l/5 - j`; with `ord_u delta^j >= 2j` and `ord_u v = 1`, every term of the expanded `(l, j)` family has u-order
`>= 2 - 2g + j + 4l/5`. So the family over `(l, j, frame labels)` is summable (for each `N` finitely many `(l, j)`, and for each `(l, j)`
the frame family is summable by (b2)); grouping by `m = l + j` gives the left side, grouping by `(l, j)` the right side. (e3) For `j >= 1`,
ML_K(g)(i) with `X_1 = .. = X_j = H` (`Pi_0(H) = 1`) says that `<H^j>_g - F^m_g^((j))(v)` lies in `Q[v]`; both terms are homogeneous of
exponent `2 - 2g - j < 0`, so the difference is `0`; for `j = 0`, ML_K(g)(ii) gives `F^m_g = c_g v^(2-2g)`. Hence `Sum_j (delta^j/j!) <H^j>_g =
c_g v^(2-2g) Sum_j binom(2-2g, j) (delta/v)^j = c_g v^(2-2g) (1 + delta/v)^(2-2g) = c_g (v + delta)^(2-2g)`: `delta/v` lies in `u Q[[u]]`
(`v = u` times a unit), and the binomial identity `(1 + x)^r = Sum_j binom(r, j) x^j` of `Q[[x]]` (`r = 2 - 2g`) is preserved by the
continuous substitution `x -> delta/v`. For `l >= 1`, each term of the frame expansion is a correlator of elements of `K` with at least
one frame (`Pi_0 = 0`), which lies in `Q[v]` by ML_K(g)(i), times an element of `Q[[u]]`; the sum of a summable family in `Q[[u]]` is in
`Q[[u]]`. QED
*Remark (a shorter route, no infinite family).* By (e1), `Sum_m (1/m!) <L^m>_g = Sum_(l + j <= 3g-3) delta^j/(l! j!) <t'^l H^j>_g`, a
finite sum. Its `l >= 1` terms are regular as above. Its `l = 0` terms are `Sum_(j <= 3g-3) delta^j F^m_g^((j))(v)/j! = c_g (v +
delta)^(2-2g) - Sum_(j > 3g-3) delta^j F^m_g^((j))(v)/j!`, and every term of the tail has u-order `>= 2j + 2 - 2g - j >= g >= 2`, so the
tail lies in `Q[[u]]`. Either route can replace the sentence "the sum over leaves may be extended to all `m` and re-ordered" (TRF s.2).

## 5. (f) Dilaton re-expansion and `d(0) = 1`

**Statement.** Let `2h - 2 + n > 0`, `X_1, .., X_n` in `HH^u_+`, `kappa` in `Q[[u]]^x`, `T` in `z^2 HH^u_+`, and `c := 1 - 1/kappa`.
(f1) If `kappa(0) = 1`, the family `(c^(m_1)/(m_1! m_2!) <X_1 .. X_n (z 1)^(m_1) T^(m_2)>_h)_(m_1, m_2)` is summable and
`Sum_m (1/m!) <X_1 .. X_n (c z 1 + T)^m>_h = kappa^(2h-2+n) Sum_m (1/m!) <X_1 .. X_n (kappa T)^m>_h`,
the right side a finite sum (`m <= 3h - 3 + n`). (f2) If `kappa(0) != 1` and `<X_1 .. X_n T^(m_2)>_h != 0` for some `m_2`, the series in
`m_1` has no u-adic limit.
*Proof.* The dilaton equation (`lem:set-dilaton`, correlator form `<z 1 Y_1 .. Y_r>_h = (2h - 2 + r) <Y_1 .. Y_r>_h` for `2h - 2 + r > 0`),
applied `m_1` times, gives `<X (z 1)^(m_1) T^(m_2)>_h = (N)_(m_1) <X T^(m_2)>_h`, `N := 2h - 2 + n + m_2 >= 1`, `(N)_(m_1) := N (N + 1) ..
(N + m_1 - 1)`. So the `m_1`-series is `<X T^(m_2)>_h Sum_(m_1) binom(N + m_1 - 1, m_1) c^(m_1)`. If `kappa(0) = 1`, `ord_u c >= 1` and the
binomial series converges to `(1 - c)^(-N) = kappa^N`; the family is summable (`m_2 <= 3h - 3 + n` by (d), `ord c^(m_1) >= m_1`, and the
correlators lie in one `u^(e_0) Q[[u]]`), and regrouping gives `Sum_(m_2) kappa^(2h-2+n+m_2) <X T^(m_2)>_h/m_2!`, the right side. If
`kappa(0) != 1`, then `c(0) != 0` and the lowest u-coefficient of the `m_1`-th term is `binom(N + m_1 - 1, m_1) c(0)^(m_1)` times that of
`<X T^(m_2)>_h`, non-zero for every `m_1`: no limit. QED
**Consequence for the note.** With `kappa = d` and `T = T_(>=2)` (`lem:tra-T` (c): `q_T = -(z/d) 1 + T_(>=2)`, so `q_T + z 1 = c z 1 + T_(>=2)`),
(f1) is `lem:tra-F` step 1, `d^n Fbar_h(q_T) = d^(2h-2+n) d^n Fbar_h(d q_T)`, PROVED when `d(0) = 1`, which holds for `P_can`
(`lem:gau-diag`: `p_(0,0)(0) = 1`, `d = 1/p_(0,0)` by (Sp) at `z^0`). For a general dressing (`def:tra-dressing` allows `d(0) != 1`, e.g.
`D P_A` with `D_0 = diag(lambda, mu, 1/mu, 1/lambda, 1)`) the left side is undefined by (f2), and `def:red-flatsum` should DEFINE the vertex
tensor at `q_T` as `d^(2h-2+n) d^n Fbar_h(d q_T)` (this is the form in which S15 s.1.1 and the computations use it; it is also what the
R-matrix action of `R = P_0^-1 P`, `R(0) = 1`, produces after the constant frame change, cf. S15 Lemmas 1-2). Theorem R is unaffected.
