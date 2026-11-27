# Gap cluster G3, part 2 of 4: item G3-1, the frame lemma with step 4 written out; G3-2g, frame expansion of edges

Continues `GAPS_G3.md` (conventions s.2 there: the graded setting over `Q[v]`, the u-adic setting `HH^u`, summable families, maps of
finite type, the summation lemma). Notation of NOTATION.md / NOTATION_2.md: `H = Q^5` with basis `H^a`, `a` in `Z/5`; `nu = (0, -1, 3,
2, 1)`, `eps = (1, -1, -1, 1, 1)`; lowered components, `<Y, n> = Sum_a Y_a n^a`; `Pi_i(q) = -[z^i] <Y^C(-z), q(z)>` (`def:tra-tower`);
`K`, `K_0` (`def:tra-slice`); `frame(k,a)`; a flat insertion `psi^k H^b` is the vector `z^k H^b`; `q_(k,b) := [z^k] q^b`.

## 1. Proposed statement of `lem:tra-frames`

Let `A` be `Q[v]` (graded setting; `HH^A_+ := H (x) Q[v][z]`, `A{z} := Q[v][z]`) or `Q[[u]]` with `v = v(u)` in `u Q[[u]]` (u-adic
setting; `HH^A_+ := HH^u_+`, `A{z} := Q[[u, z]]`); write `K^A`, `K_0^A` for the slice and its kernel part in `HH^A_+`.
* (a) `HH^A_+ = K_0^A (+) A{z} H`. In the graded setting `V_w = (K_0 cap V_w) (+) (A{z} H cap V_w)` for every `w`.
* (b) For every `(k, a)` with `a != 1` there is exactly one element of `K_0^A` of the form `z^k H^a + s(z) H`, `s` in `A{z}`: it is
  `frame(k,a) = z^k H^a + Sum_(j<k) f_j z^j H`, `f_j = c_(k,a,j) v^((4(k-j) + nu_a + 1)/5)`, `c_(k,a,j)` in `Q`, and `c_(k,a,j) = 0`
  unless the exponent is an integer; it is then `>= 1`. No element of `K_0^A` has the form `z^k H + Sum_(j<k) s_j z^j H`.
* (c) (frame coordinates) The map `K_0^A -> prod_(k >= 0, a != 1) A`, `n -> (n_(k,a))`, is bijective (graded setting: onto the
  families with finitely many non-zero members, each homogeneous); the inverse sends `(c_(k,a))` to the sum of the summable family
  `(c_(k,a) frame(k,a))`. In particular `pr(q) = Sum_(k, a != 1) q_(k,a) frame(k,a)` for the projection `pr` of (a).
* (d) The frames of the u-adic setting are the graded frames with `v = v(u)`.

## 2. Proof (Steps 1-5 correspond to DC1 A1 steps 1-5; Step 4 is the step marked [S])

**Step 0 (two values of `Pi_0`).** By `lem:gau-vanshape`, `[z^-n] Y^C_b = eps_b (-1)^n v^(n-q) / (4^q q! (n-q)!)`, `q = floor(n/5)`, if
`n = nu_b + 1 (mod 5)`, and `0` otherwise. With `def:tra-tower`, `Pi_0(z^m H^b) = (-1)^(m+1) [z^-m] Y^C_b`, i.e.
`Pi_0(z^m H^b) = -eps_b v^(m-q)/(4^q q! (m-q)!)` if `b + m = 1 (mod 5)`, and `0` otherwise. In particular
* (0.1) `Pi_0(z^m H) = pi_q v^(4q)` if `m = 5q`, and `0` if `5` does not divide `m`, where `pi_q := 1/(4^q q! (4q)!)`, `pi_0 = 1`;
* (0.2) `Pi_0(H^a) = 0` for `a != 1` (the condition `b + m = 1` fails at `m = 0`);
* (0.3) the `v`-exponent of `Pi_0(z^m H^b)` is `m - floor(m/5) >= 4m/5`; so in the u-adic setting `ord_u Pi_0(z^m H^b) >= 4m/5`.
Examples: `Pi_0(H) = 1`, `Pi_0(z 1) = -v`, `Pi_0(z^2 H^4) = -v^2/2`, `Pi_0(z^5 H) = v^4/96`, `Pi_0(z^6 1) = -v^5/480`.

**Step 1 (`Pi_i` on `HH^A_+`; the operator `Theta`).** For `q` in `HH^A_+`, the shift rule gives
`Pi_i(q) = Sum_(k >= i, b) q_(k,b) Pi_0(z^(k-i) H^b)`.
In the graded setting this is a finite sum. In the u-adic setting it is the product of `Y^C(v(u), -z)` (in `Q[z^-1][[u]]^5`) and `q` (in
`Q[[z]][[u]]^5`) in the ring `Q((z))[[u]]` (`GAPS_G3.md` s.2.1); by (0.3) the coefficient of `q_(k,b)` has u-order `>= 4(k-i)/5`, so
`Pi_i(q)` modulo `u^N` depends only on the `q_(k,b)` with `k < i + 5N/4 + 1`, taken modulo `u^N`: `Pi_i` is of finite type
(`GAPS_G3.md` s.2.2). For a sequence `s = (s_j)_(j >= 0)` in `A` (finitely many non-zero terms in the graded setting) put `s(z) :=
Sum_j s_j z^j` and `Theta s := (Pi_i(s(z) H))_(i >= 0)`. By the shift rule and (0.1),
        `(Theta s)_i = s_i + (N s)_i`,   `(N s)_i := Sum_(q >= 1) pi_q v^(4q) s_(i+5q)`.
*Claim: `Theta` is bijective, with inverse `Sum_(r >= 0) (-N)^r`.* Graded setting: if `s` vanishes outside `[0, L)`, then `N s`
vanishes outside `[0, L - 5)`, so `N^r s = 0` for `5r >= L` and the sum is finite. u-adic setting: `N` is well defined (the sum over `q`
converges, `ord_u v^(4q) >= 4q`), and if `ord_u s_j >= m` for all `j` then `ord_u (N s)_i >= m + 4` for all `i`; hence `(N^r s)_i` has
u-order `>= 4r + min_j ord_u s_j` uniformly in `i`, the family `((-N)^r s)_r` is summable coefficientwise, `N` commutes with its sum
(the double family `pi_q v^(4q) ((-N)^r s)_(i+5q)` over `(q, r)` is summable at each `i`), and
`(1 + N) Sum_(r <= R) (-N)^r s = s - (-N)^(R+1) s` tends to `s`. So `Theta Sum_r (-N)^r = 1`, and in the same way `Sum_r (-N)^r Theta
= 1`. In both settings `N` maps sequences vanishing outside `[0, L)` to sequences vanishing outside `[0, L)`, hence so does
`Theta^-1`. In the graded setting `N` maps a sequence with `s_j` of weight `w - 4j` (for all `j`) to one of the same kind, since
`v^(4q) s_(i+5q)` has weight `20q + w - 4i - 20q`; so `Theta^-1` preserves this homogeneity.

**Step 2 (existence of the decomposition: (a)).** For `q` in `HH^A_+` put `sigma(q) := Theta^-1 ((Pi_i(q))_(i >= 0))` and
`pr(q) := q - sigma(q)(z) H`. Then `Pi_i(pr(q)) = Pi_i(q) - (Theta sigma(q))_i = 0` for every `i >= 0`, so `pr(q)` lies in `K_0^A`, and
`q = pr(q) + sigma(q)(z) H`. If `n = s(z) H` lies in `K_0^A`, then `Theta s = 0`, so `s = 0`: the sum `K_0^A + A{z} H` is direct.
Graded setting: for `q` in `V_w` the sequence `(Pi_i(q))_i` vanishes for `i` larger than the `z`-degree of `q` and `Pi_i(q)` has weight
`w - 4i` (`GAPS_G3.md` s.2.1); by Step 1 `sigma(q)_j` has weight `w - 4j`, so `sigma(q)(z) H` and `pr(q)` lie in `V_w`. This is (a).
*Remark (source step 3).* The source proves uniqueness by "a combination `Sum s_j z^j H` in `K_0` has all `s_j = 0` (look at its top
power)". That argument covers polynomial combinations only, i.e. the graded setting. In the u-adic setting the elements of
`A{z} H` are power series in `z` (edges and leaf points have infinitely many `z`-powers), there is no top power, and the
injectivity of `Theta` (Step 1, u-adic) is what proves uniqueness. The source's elimination "from the top" (step 2) is
`Theta^-1` applied to a finitely supported vector.

**Step 3 (frames: (b), (d)).** For `a != 1` put `frame(k,a) := pr(z^k H^a)`. The sequence `(Pi_i(z^k H^a))_i = (Pi_0(z^(k-i) H^a))_i`
vanishes for `i > k` (shift rule) and at `i = k` by (0.2); so it vanishes outside `[0, k)`, hence so does `sigma(z^k H^a)` (Step 1),
and `frame(k,a) = z^k H^a + Sum_(j<k) f_j z^j H` with `f_j := -sigma(z^k H^a)_j`. *Uniqueness:* if `z^k H^a + s(z) H` lies in
`K_0^A` for some `s` in `A{z}`, its difference with `frame(k,a)` lies in `K_0^A cap A{z} H = 0`. *No frame with `a = 1`:* an element
`z^k H + Sum_(j<k) s_j z^j H` lies in `A{z} H`; if it lay in `K_0^A` it would be `0` by (a), but its `z^k H`-coefficient is `1`.
*Exponents* (graded setting): `z^k H^a` lies in `V_w`, `w = 4k + nu_a + 1`, so by Step 2 `f_j` has weight `w - 4j`, i.e. `f_j =
c_(k,a,j) v^((4(k-j) + nu_a + 1)/5)` with `c_(k,a,j)` in `Q`, and `c_(k,a,j) = 0` unless the exponent is a non-negative integer;
for `j <= k - 1` and `a != 1` (`nu_a + 1 >= 1`) it is `>= (4 + 1)/5 = 1`. Integrality of the exponent is the support rule
`j = k - nu_a - 1 (mod 5)`. *(d):* a graded frame is a polynomial in `z` with coefficients in `Q[v]`; `v -> v(u)` is a ring map
`Q[v] -> Q[[u]]` that commutes with each `Pi_i` on polynomials (finite sums of products), so the substituted frame lies in the
u-adic `K_0`, has the form `z^k H^a + (polynomial) H`, and equals the u-adic frame by uniqueness.
Examples (these are the outline's, computed here from Step 1 and (0.1); appendix C may print them): `frame(1,0) = z 1 + v H = e`,
`frame(2,4) = z^2 H^4 + (v^2/2) H`, `frame(3,3) = z^3 H^3 + (v^3/6) H`, `frame(4,2) = z^4 H^2 - (v^4/24) H`, `frame(6,0) = z^6 1 +
v z^5 H - (v^5/120) H` (here `f_5 = -Pi_0(z 1) = v`, `f_0 = -Pi_0(z^6 1) - pi_1 v^4 f_5 = v^5/480 - v^5/96`).
*Trivial frames.* `frame(k,a) = z^k H^a` iff `k < nu_a + 1`: if `k >= nu_a + 1`, the index `j_0 = k - nu_a - 1` is admissible and the
recursion `f_i = -Pi_0(z^(k-i) H^a) - Sum_(q >= 1) pi_q v^(4q) f_(i+5q)` (Step 1 read backwards) has no correction at `j_0` (`j_0 + 5q > k`),
so `f_(j_0) = -Pi_0(z^(nu_a+1) H^a) = eps_a v^(nu_a+1)/(nu_a+1)! != 0`; if `k < nu_a + 1` no `j < k` is admissible. So the trivial frames
are exactly the ten `(0,0)`, `(0,2), (1,2), (2,2), (3,2)`, `(0,3), (1,3), (2,3)`, `(0,4), (1,4)`. The list of eight in RD1 s.1.3 / DC1 D9
omits `(2,2)` and `(3,2)` (an example list, not a classification; do not print it as complete).

**Step 4 (the step marked [S]: frame coordinates (c), the decomposition of infinite elements).** Let `c = (c_(k,a))_(k >= 0, a != 1)`
be a family in `A` (finitely many non-zero, homogeneous of the right weights, in the graded setting).
* *Summability.* The coefficient `[z^j]^b` of `c_(k,a) frame(k,a)` is: for `b != 1`, `c_(j,b)` if `(k,a) = (j,b)` and `0` otherwise; for
  `b = 1`, `c_(k,a) f_j^((k,a))(v(u))`, of u-order `>= (4(k-j) + nu_a + 1)/5 >= (4(k-j) + 1)/5` (Step 3, and `ord_u v >= 1`). For fixed
  `(j, b, N)` only finitely many `(k, a)` give a coefficient outside `u^N Q[[u]]`. So `(c_(k,a) frame(k,a))` is summable; put
  `n(c) := Sum c_(k,a) frame(k,a)`. (The family is also `(u, z)`-adically summable: `frame(k,a)` lies in `(u, z)^(floor(4k/5)) (x) H`,
  since `j + (4(k-j) + 1)/5 >= (4k + 1)/5` and `k >= 4k/5`.) In the graded setting the family is finite.
* *`n(c)` lies in `K_0^A`.* Each `Pi_i` is of finite type (Step 1), so by the summation lemma `Pi_i(n(c)) = Sum c_(k,a) Pi_i(frame(k,a)) = 0`.
* *Its flat components.* By the first bullet, `n(c)_(k,a) = c_(k,a)` for `a != 1`.
* *Injectivity and surjectivity.* For `n` in `K_0^A` put `c_(k,a) := n_(k,a)`. Then `n - n(c)` lies in `K_0^A` and has no flat component
  with `a != 1`, i.e. lies in `A{z} H`; by (a) it is `0`. So `n = n(c)`, and `n -> (n_(k,a))` is bijective with inverse `c -> n(c)`.
* *The projection.* For any `q`, `pr(q)` and `Sum_(k, a != 1) q_(k,a) frame(k,a)` lie in `K_0^A` and have the same components with
  `a != 1` (those of `q`, since `pr` changes only `H`-components); by injectivity they are equal. Hence
  `q = Sum_(k, a != 1) q_(k,a) frame(k,a) + sigma(q)(z) H`,   `sigma(q)_j = q_(j,1) - Sum_(k > j, a != 1) q_(k,a) f_j^((k,a))`,
  the last sum being u-adically convergent (finite in the graded setting). This is the sentence "decompose each flat component with
  `a != 1` by step 2, keep the `H` components" of the source, with its convergence. QED (c).

**Step 5 (homogeneity; `lem:tra-coeff`).** Step 3 gives the exponent count `wt(k,a) - wt(j,1) = (4(k-j) + nu_a + 1)/5 >= 1` of TRF
s.2 / ST2 M5. With Step 4: every element of `K_0` (u-adic) is the sum of the summable family `(n_(k,a) frame(k,a)(v(u)))`, a
`Q[[u]]`-combination of homogeneous elements of `K_0` (graded); every element of `K^u` adds `Pi_0(q) H` (`cor:tra-slice` below).

**Corollaries.** (i) (`cor:mod-basis`, see `GAPS_G3.md` s.3 item 1) The `v^m frame(k,a)` (`a != 1`) of weight `w` are a basis of
`K_0 cap V_w`. (ii) (`cor:tra-slice`) `K^A = K_0^A (+) A H`: `Pi_i(H) = 0` for `i >= 1` (shift rule) and `Pi_0(H) = 1`, so `q` in
`K^A` gives `q - Pi_0(q) H` in `K_0^A`. With `Pi_0(-z 1) = v`, `Pi_i(-z 1) = 0` for `i >= 1` (Step 0 and (0.2)):
`-z 1 + K^A = {q in -z 1 + HH^A_+ : Pi_i(q) = 0, i >= 1}`, and `Pi_0(-z 1 + n + c H) = v + c` for `n` in `K_0^A`, `c` in `A`
(write a local letter for `c`, R8 (a)).

## 3. G3-2g: frame expansion of bivectors (for `cor:tra-edges` step 5; u-adic setting)

For `E` in `H (x) H (x) Q[[psi, psi']][[u]]` write `E_(kl)^(ab) := [psi^k psi'^l] E^(ab)`. For fixed `(l, b)` let `E_(.,(l,b)) :=
Sum_(k,a) E_(kl)^(ab) psi^k H^a` (an element of `HH^u_+` in the variable `psi`), and define `(Pi_i (x) 1) E := Sum_(l,b) Pi_i(E_(.,(l,b)))
psi'^l H^b`; `(1 (x) Pi_i) E` likewise with the slots exchanged.
**Lemma.** The following are equivalent: (i) `(Pi_i (x) 1) E = 0` and `(1 (x) Pi_i) E = 0` for all `i >= 0`; (ii) `E` is the sum of the
summable family `(E^(alpha beta) frame_alpha(psi) (x) frame_beta(psi'))`, `alpha = (k,a)`, `beta = (l,b)`, `a, b != 1`, `E^(alpha beta) :=
E_(kl)^(ab)`. Moreover the frame coefficients in (ii) are unique.
*Proof.* Summability of the family in (ii), for any coefficients in `Q[[u]]`: the coefficient of `psi^(k') psi'^(l') H^(a') (x) H^(b')`
of the `(alpha, beta)` term is a product of a coefficient of `frame_alpha` and one of `frame_beta`, and Step 4 (first bullet) applies in
each variable. (ii) => (i): `Pi_i (x) 1` is of finite type coefficientwise in the second slot, and kills each term. (i) => (ii): for each
`(l, b)`, `E_(.,(l,b))` lies in `K_0` (first half of (i)), so by Step 4 `E_(.,(l,b)) = Sum_alpha E_(kl)^(ab) frame_alpha(psi)`, the sum
over `alpha = (k, a)`, `a != 1`. The double family `(E_(kl)^(ab) frame_alpha(psi) psi'^l H^b)` over `(alpha, (l, b))` (all `b`) is
summable (Step 4 in `psi`; the `psi'^l H^b` are distinct), so `E = Sum_alpha frame_alpha(psi) (x) G_alpha(psi')` with
`G_alpha := Sum_(l,b) E_(kl)^(ab) psi'^l H^b`. Now `1 (x) Pi_i` is of finite type coefficientwise in the first slot, so by the
summation lemma `0 = (1 (x) Pi_i) E = Sum_alpha Pi_i(G_alpha) frame_alpha(psi)`, and by the uniqueness of frame coordinates (Step 4,
in `psi`) `Pi_i(G_alpha) = 0` for all `alpha` and `i`. So `G_alpha` lies in `K_0` (in `psi'`), and Step 4 gives `G_alpha = Sum_beta
E^(alpha beta) frame_beta(psi')`. The resulting double family is summable and its sum is `E` (regrouping). Uniqueness: the
coefficients are the flat components with `a, b != 1`. QED
**Use.** With `lem:tra-E` (part 4) and (A), `E_P` satisfies (i) (`cor:tra-edges` steps 1-4), so `E_P = Sum E_P^(alpha beta) frame_alpha
(x) frame_beta`: "edges in `K_0 (x) K_0`, as `Q[[u]]`-combinations of frames", with the completed tensor product made precise as
the summable families above. This is the "duality step" that DC1 B2 step 5 asks to be written: `Pi_i` restricted to `A{z} H` is
`Theta`, which is invertible (Step 1), so the tower detects the `H`-part.
