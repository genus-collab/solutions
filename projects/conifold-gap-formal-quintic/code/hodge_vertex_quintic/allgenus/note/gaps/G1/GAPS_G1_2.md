# Gap cluster G1, part 2 of 5: item 2, `prop:set-cyclicLR` (inputs and lemmas)

Continuation of `GAPS_G1.md`. Resolution: **PROVED** (all orders in `q`), given (M) of s.1 and the genus-0 axioms (G0) of s.2.1.
Notation: `D = q d/dq`; `C_0, I_1 = C_0 log q + Itilde_1, T = I_1/C_0, C_1 = DT, L = (1 - 5^5 q)^(-1/5)` (`def:set-genus0`); the flat
basis `1, H, .., H^4` with `H^5 = 1`; `eta(H^a, H^b) = 5 delta_(a+b = 3 mod 5)`, `eta^(ab) = (1/5) delta_(a+b = 3 mod 5)`; indices `a` in Z/5.
`d_T := C_1^(-1) D` (derivative in `T`; NOT the note's `d_t` of s03, which is defined in a conifold frame; the writer should pick a
local letter, e.g. `\partial_T`). `ell` is a formal `log q` with `D ell = 1` (local to this proof; the note's `ell` of s03 is `log u`).

## 2.0 Statement (proposed text of `prop:set-cyclicLR`)

**Proposition.** Assume (M). At `t = T(q) H` the quantum product of the formal quintic satisfies, in the flat basis,
`H * H^a = kappa_a H^(a+1)` (`a` in Z/5, `H^5 = 1`) with `kappa = (kappa_0, .., kappa_4) = (1, L^5/(C_0^2 C_1^3), 1, C_0/C_1, C_0/C_1)`.
Equivalently, for every solution `v` of the quantum differential equation `z d_T v = H * v`, the lowered components `v_a := eta(H^a, v)`
satisfy `z D v_a = C_1 kappa_a v_(a+1)` for all `a` in Z/5; for the colour-`i` column `v = R(z) e~_i e^(u_i/z)` of `def:set-rows` this is
`X_(a+1) = (C_1 kappa_a)^(-1) (z D + zeta^i L) X_a` (`a = 0..3`) for `X_a := zeta^(-i) Q^(i)_a`, with `X_0 = L Phi_i/C_0` (= `psi_i/C_0` in the
sources' letters). Moreover `C_0 v_0` solves PF5 for every such `v` (Corollary 2.10 (b)).

## 2.1 Inputs

* **(M)** as in s.1.1.
* **(G0) genus-0 axioms of the twisted theory** (the genus-0 classes `Omega_(0,n)(g_1, .., g_n) := Sum_d Q^d rho_*(e_T(E_(0,n,d)) prod
  ev_j^* g_j)` form a CohFT with unit `1` for the pairing `eta`, in the sense of PPZ15 Def. 0.3 (read); GJR18 s.4.1 defines the
  `lambda`-twisted CohFT with this pairing and unit (read)). Concretely used:
  * (G0-a) *polynomiality and symmetry:* `<H^(a_1), .., H^(a_n)>_(0,n,d)` (`d >= 1`) lies in `Q[lambda_0, .., lambda_4]`, is symmetric
    under permutations of the `lambda_j`, and is homogeneous of degree `Sum a_j - n` (s.1.4: `E_(0,n,d)` is a vector bundle of rank `5d+1`,
    `dim Mbar_(0,n)(P^4,d) = 5d + n + 1`; the moduli space is proper; a permutation of the coordinates of `C^5` is an automorphism of `P^4`
    that permutes the weights and preserves `H` and `O(5)` with their linearisations).
  * (G0-b) *string and divisor equations* (`E_(0,n+1,d) = pi^* E_(0,n,d)` for the map forgetting a point; `pi_* ev^* H = d`).
  * (G0-c) *splitting, hence WDVV*: the normalisation sequence along a boundary divisor gives `e(E)|_(boundary) e(ev_node^* O(5)) =
    e(E_1) e(E_2)`, i.e. the splitting axiom for the pairing `eta`; and the genus-0 TRR (`psi_1 = ` boundary on `Mbar_(0,n)(P^4,d)`).
  * (G0-d) *QDE:* with `J` of s.1.1 and `s_a := d J/d t^a` (derivative along `H^a`), `z d_(t^b) s_a = Sum_c (H^b * H^a)^c s_c`, where
    `x * y := Sum_(c,c') <<x, y, H^c>> eta^(c c') H^(c')`; `s_a = H^a + O(1/z)`; `s_0 = J/z` (string equation). This follows from (G0-c) by
    the standard computation (TRR applied to `z <<H^b, H^a, H^c/(z - psi)>>`); stated for the formal quintic in GJR18 s.5.2 eq. (14)
    ("Recall the S-matrix is the solution of the quantum differential equation (QDE) [13, 25] `z dS_t(z) = dt *_t S_t(z)`", read).
* Elementary: `lem:set-pf5` (`PF5 = z^5 D PF4 - 1`, i.e. `D^5 - q prod_(k=1..5)(5D + k) = D (D^4 - 5q prod_(k=1..4)(5D + k))` because
  `D q = q(D + 1)`).

## 2.2 Lemma A (Z/5 grading)

*At the zeta point, `<<H^a, H^b, H^c>>(TH) = 0` unless `a + b + c = 3 mod 5`; hence `H^a * H^b` is a multiple of `H^(a+b mod 5)`, and
`H * H^a = kappa_a(q) H^(a+1)` with `kappa_a` in `Q[[q]]`, `kappa_a = 1 + O(q)`.*

*Proof.* By (G0-b), `<<H^a,H^b,H^c>>(TH) = eta(H^a H^b, H^c) + Sum_(d>=1) e^(dT) <H^a,H^b,H^c>_(0,3,d)` (with `Q = 1`; for `d = 0` the
insertions of `TH` vanish because `Mbar_(0,3+n)(P^4,0) = Mbar_(0,3+n) x P^4`). By (G0-a) the degree-`d` term is a symmetric polynomial
of degree `m = a + b + c - 3`, i.e. a polynomial in the elementary symmetric functions `e_1, .., e_5` of weighted degree `m`. At
`lambda_j = zeta^j`, `prod_j (x - zeta^j) = x^5 - 1` gives `e_1 = .. = e_4 = 0`, `e_5 = 1`; so the value is 0 unless `5 | m`. The classical
term is `5 delta_(a+b+c = 3 mod 5)`. With `eta^(cc')` supported on `c + c' = 3 mod 5`, `H^a * H^b` is supported on `c' = a + b mod 5`.
`kappa_a` is a power series in `e^T` (structure constants), hence in `q` by the mirror map `e^T = q exp(Itilde_1/C_0)`; at `q = 0` the
product is the classical one (`H * H^a = H^(a+1)`), so `kappa_a(0) = 1`. QED

## 2.3 Lemma B (unit and Frobenius)

*`kappa_0 = kappa_2 = 1` and `kappa_3 = kappa_4`.*

*Proof.* `1` is the unit of `*` (string equation), so `H * 1 = H`: `kappa_0 = 1`. The Frobenius property `eta(x * y, w) = eta(x, y * w)`
(both equal `<<x, y, w>>`) with `x = H`, `y = H^a`, `w = H^(2-a)` gives `5 kappa_a = eta(kappa_a H^(a+1), H^(2-a)) = eta(H^a, kappa_(2-a) H^(3-a)) =
5 kappa_(2-a)`; so `kappa_a = kappa_(2-a)` (indices mod 5): `kappa_0 = kappa_2`, `kappa_3 = kappa_4`. QED

## 2.4 Lemma C (Birkhoff recursion)

*Write, by (M), `s_0 = J/z = I/(z C_0) = Sum_(k>=0) phi_k H^k z^(-k)` (`H^k := H^(k mod 5)`; distinct `k` carry distinct powers of `z`,
so the coefficients are well defined), `phi_k := c_k/C_0` (Lemma D below). Put `phi^(0)_k := phi_k` and, recursively,
`kappa_a = d_T phi^(a)_(a+1)`, `phi^(a+1)_k := d_T phi^(a)_k / kappa_a`. Then `s_a = Sum_(k>=a) phi^(a)_k H^k z^(a-k)` with `phi^(a)_a = 1`
for `a = 0..4`, and these `kappa_a` are the structure constants of Lemma A.*

*Proof.* Induction on `a`; `a = 0` is (M) with `phi_0 = c_0/C_0 = 1`. By (G0-d) with `b = 1` and Lemma A, `z d_T s_a = kappa_a s_(a+1)`
(`a = 0..3`). Differentiating the expansion, `z d_T s_a = Sum_(k>=a+1) (d_T phi^(a)_k) H^k z^(a+1-k)` (the `k = a` term is `z d_T 1 = 0`).
Since `s_(a+1) = H^(a+1) + O(1/z)`, the `z^0` coefficient gives `kappa_a = d_T phi^(a)_(a+1)`, a unit of `Q[[q]]` (Lemma A), and
`s_(a+1) = z d_T s_a / kappa_a` gives the expansion of `s_(a+1)` with `phi^(a+1)_k = d_T phi^(a)_k/kappa_a`, `phi^(a+1)_(a+1) = 1`. QED

In closed form (used below; `d_T = C_1^(-1) D`): `kappa_1 = d_T^2 phi_2`, `kappa_2 = d_T(d_T^2 phi_3 / kappa_1)`, and, when `kappa_2 = 1`,
`kappa_3 = d_T^2 (d_T^2 phi_4 / kappa_1)`. (Also `phi_1 = T`, so `d_T phi_1 = 1`, consistent with `kappa_0 = 1`.)

## 2.5 Lemma D (the periods `c_k`)

*Let `w` be a formal variable and `Fhat(w, q) := Sum_(d>=0) q^d prod_(k=1..5d)(5w + k) / prod_(k=1..d)((w + k)^5 - w^5)` (in `Q[[q]][[w]]`;
the denominators have constant term `(d!)^5`). The map `Q[[w]] -> (Q[H]/(H^5 - 1))[[z^(-1)]]`, `w^k |-> H^(k mod 5) z^(-k)`, is an injective
ring homomorphism, and `I/z` is the image of `e^(ell w) Fhat(w, q)`. Hence `c_k = Sum_(j=0..k) (ell^j/j!) ctilde_(k-j)` with `ctilde_m :=
[w^m] Fhat` in `Q[[q]]`, `ctilde_m(0) = delta_(m0)`; `c_0 = C_0`, `c_1 = I_1`. Moreover `PF4 c_k = 0` for `k = 0..3` and `PF4 c_4 = 1`.*

*Proof.* Multiplicativity: `H^5 = 1`. At the zeta point `prod_j (H - zeta^j + kz) = (H + kz)^5 - 1`, and under the map `(w + k)^5 - w^5 |->
z^(-5)((H + kz)^5 - 1)`, `5w + k |-> z^(-1)(5H + kz)`; the powers of `z` cancel in each `q^d` term, and `q^(H/z) = e^(ell H/z)`.
*PF5.* For each fixed point, `prod_j (zD - zeta^j)` acting on the `q^d` term of `I_i` multiplies it by `prod_j (zeta^i - zeta^j + dz)`, which
cancels the `k = d` factors of the denominator and leaves `q prod_(k=1..5)(5 zeta^i + (5d - 5 + k) z)` times the `q^(d-1)` term, i.e.
`q prod_(k=1..5)(5 zD + kz)` applied to it; and `prod_j (zD - zeta^j) = (zD)^5 - 1`. So PF5 `I_i = 0`. With `lem:set-pf5` and
`I_i/z = Sum_k zeta^(ik) c_k z^(-k)`: `Sum_k zeta^(ik) (z^(5-k) D PF4 c_k - z^(-k) c_k) = 0`; the coefficients of `z^1, .., z^5` give
`D PF4 c_k = 0` for `k = 0..4`. *Constants.* `Q[[q]][ell]` with `D ell = 1` has constant ring `Q` (for `f = Sum_(j<=N) a_j ell^j`, `Df = 0`
forces `D a_N = 0`, so `a_N` in `Q`, then `D a_(N-1) = -N a_N`; but `D(Q[[q]])` has no constant term, so `N a_N = 0`). Hence `PF4 c_k`
is a constant, equal to its `q^0` part. Taking the `q^0` part commutes with `D` (it acts as `d/d ell`) and kills `q`; `PF4 = (1 - 3125q)
D^4 - q(...)`, so the `q^0` part of `PF4 c_k` is `(d/d ell)^4 (ell^k/k!)` = 0 for `k <= 3`, = 1 for `k = 4`. QED

## 2.6 Lemma E (Wronskian)

*`W := det(D^j c_k)_(j,k = 0..3) = L^10`.*

*Proof.* Write `PF4 = a_4 D^4 + a_3 D^3 + ..`, `a_4 = 1 - 3125 q` (a unit), `a_3 = -6250 q`. Differentiating the determinant row by row,
only the last row changes to `D^4 c_k = -a_4^(-1)(a_3 D^3 c_k + ..)` (Lemma D), so `D W = -(a_3/a_4) W = 6250 q (1 - 3125q)^(-1) W`. Also
`D L^10 = 6250 q (1 - 3125 q)^(-1) L^10` (`L^10 = (1 - 3125 q)^(-2)`), so `D(W L^(-10)) = 0` and `W L^(-10)` is a constant (Lemma D), equal to
the `q^0` part of `W`, which is `det((d/d ell)^j ell^k/k!) = 1` (upper unitriangular). QED

## 2.7 Lemma F (operators of order at most 3)

*If `Theta = Sum_(j=0..3) r_j D^j` with `r_j` in `Q((q))` satisfies `Theta c_k = 0` for `k = 0..3`, then `Theta = 0`.*

*Proof.* Since `r_j` do not depend on `ell`, `Theta(e^(ell w) G) = e^(ell w) Sum_j r_j (D + w)^j G` for `G` in `Q((q))[[w]]` (Leibniz,
`D e^(ell w) = w e^(ell w)`). With `G = Fhat`: `Sum_k (Theta c_k) w^k = e^(ell w) Ghat(w)`, `Ghat := Sum_j r_j (D + w)^j Fhat`. The hypothesis
says `e^(ell w) Ghat = O(w^4)`, hence `Ghat = O(w^4)` (`e^(ell w)` is invertible with constant term 1). The coefficients `[w^n] Ghat =
Sum_j m_(nj) r_j`, `m_(nj) = Sum_(i <= min(j,n)) binom(j,i) D^(j-i) ctilde_(n-i)`, `n, j = 0..3`. Modulo `q`: `D^e ctilde_m = 0` unless
`e = 0` and `m = 0` (`ctilde_m(0) = delta_(m0)`, and `D` kills `q^0`), so `m_(nj) = delta_(nj) mod q`; `det(m_(nj))` is a unit of
`Q[[q]]`, and `r_j = 0`. QED
