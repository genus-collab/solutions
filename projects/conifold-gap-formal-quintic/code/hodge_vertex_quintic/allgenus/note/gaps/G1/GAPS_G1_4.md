# Gap cluster G1, part 4 of 5: item 3 (`thm:set-qrr`) and item 4 (`lem:set-dilaton`)

Continuation of `GAPS_G1_3.md`. Notation as in `GAPS_G1_2.md`; `b_k := B_2k/(2k(2k-1))` (the note's sans-serif `b_k`), `h_w(z) :=
exp(-Sum_k b_k z^(2k-1) Sum_(a in w) a^(1-2k))`, `w_i = {zeta^i - zeta^j}_(j != i) U {-5 zeta^i}`, `h = h_(w_0)`, `h_n = [z^n] h`.

## 3. Item 3: `thm:set-qrr`

(QRR) as printed: in Givental's convention (`def:set-rows`: the columns `R(z) e~_i e^(u_i/z)` solve `z d_T v = H * v`, `R = 1 + O(z)`,
symplectic), `R(q = 0) = diag(h_(w_0), .., h_(w_4))` in the normalised canonical basis `e~_i(q=0)`. Here `R` is a formal series in `z` whose
coefficients are power series in `q` along the line `t = T(q) H` (GJR18 s.4.1: `R_tau(z)` in `End(H) (x) Q(lambda)[[q, z]]`, read), and
`R(q = 0)` is the constant term in `q`.

### 3.1 (a) The passage to `A^QRR_n(q = 0) = h_n`: PROVED

*Definitions (recommended for `def:set-rows`).* `Q^(i)_a(z) := eta(H^a, R(z) e~_i)`; `Phi_i(zhat, u) := C_0 Q^(i)_0(z)/(zeta^i L)` with
`zhat = L^4 z`; flat amplitude `A_n := [z^n](L Phi_0) = [z^n](C_0 Q^(0)_0)`, so `A_n = L^(4n+1) ph_n(L^(-5))` when `Phi_0 = Sum_n zhat^n ph_n(u)`.
(This is the normalisation of the sources: DA D6, "`Q = P/Delta^(1/2) = zeta^i X_a`", `X_0 = psi_i/C_0`, `psi_i = L Phi_i`.) By Corollary 2.10
(b), `e^(u_i/z) z L Phi_i` solves PF5: the "WKB form of colour `i`" is a consequence, not a definition.

*Claim.* `Phi_i(z, 1) = h_(w_i)(z)` for every `i`; in particular `Sum_n A^QRR_n(q=0) z^n = Phi_0(z, 1) = h(z)`, i.e. `A^QRR_n(q = 0) = h_n`.

*Proof.* `q = 0` is `u = 1`, `L = 1`, `zhat = z`, `C_0 = 1`. For every `q`, `eta(1, e~_i) = Delta_i^(1/2) eta(1, e_i) = Delta_i^(1/2) eta(e_i, e_i)
= Delta_i^(-1/2) = zeta^(-4i) L/C_0 = zeta^i L/C_0` (Corollary 2.10 (a) and the convention `Delta_i^(1/2) = zeta^(4i) C_0/L`, NOTATION_2 s.10
item 7). By (QRR), `R(0, z) e~_i(0) = h_(w_i)(z) e~_i(0)`, so `Q^(i)_0(q = 0) = h_(w_i)(z) eta(1, e~_i)|_(q=0) = zeta^i h_(w_i)(z)`, and
`Phi_i(z, 1) = Q^(i)_0(q=0)/zeta^i = h_(w_i)(z)`. QED (Consistency: `Phi_i(zhat, u) = Phi_0(zeta^(-i) zhat, u)` and `h_(w_i)(z) = h(zeta^(-i) z)`,
`lem:set-hrot`.) The `\draftnote` (a) of OUTLINE_3 s02 item 6 can be deleted once the definitions above are adopted.

### 3.2 (b) Uniqueness: `R(q = 0)` is the vertex R-matrix: CITED (Teleman), with the argument written

*Step 1 (PROVED here; QDE only).* `R(0, z)` is diagonal in `e~_i(0)` with entries `c_i(z)`, `c_i(z) c_i(-z) = 1`, `c_i(0) = 1`. Proof: the
columns satisfy `(z D + D u_i) R e~_i = C_1 H * (R e~_i)` (`D = C_1 d_T`). Take the `q^0` part: `D` kills it, `D u_i(0) = zeta^i`, `C_1(0) = 1`,
and `H*` at `q = 0` is the classical product. So `R(0, z) e~_i(0)` lies, coefficientwise in `z`, in the `zeta^i`-eigenline of classical
multiplication by `H`, which is spanned by `e~_i(0)` (the eigenvalues `zeta^0..zeta^4` are distinct). Symplecticity (PPZ15 (5),
`R(z) R^*(-z) = 1`; `eta(e~_i, e~_j) = delta_ij`) gives `c_i(z) c_i(-z) = 1`. Hence `c_i = exp(odd series)`.

*Step 2 (the CohFT at `q = 0`).* The `q^0` part of the formal quintic CohFT comes from constant maps: `Mbar_(g,n)(P^4, 0) = Mbar_(g,n) x P^4`,
virtual class `e(E^vee [x] T_(P^4))`, `R pi_* f^* O(5) = O(5) - E^vee (x) O(5)`; localisation on `P^4` gives `Omega_(g,n)(g_1..g_n)|_(q=0) =
Sum_j prod_k g_k(p_j) H_g^(p_j)` with LP18 eq. (26) (read, s.3.2 p.23): `H_g^(p_i) = [e(E_g^* (x) T_(p_i)P^4)/e(T_(p_i)P^4)] [(5 lambda_i)/
e(E_g^* (x) (5 lambda_i))]`. With Mumford's relation in the form `Lv(w) Lv(-w) = (-1)^g w^(2g)`, `Lv(w) := e(E^* (x) L_w)` (LIT s.1.4 (a), whose
exact check `N_(2,0) = -19/9600`, `N_(3,0) = 1/224000` is evidence only): `H_g^(p_i) = W_i^(g-1) prod_(a in w_i) Lv(a)/a^g`, `W_i = prod_(j != i)
(lambda_i - lambda_j)/(5 lambda_i) = Delta_i(0)`, and `Lv(a)/a^g = exp(-Sum_k (2k-2)! ch_(2k-1)(E) a^(1-2k))` (`ch_(2k)(E) = 0`). In the basis `e~_i(0)`:
`Omega|_(q=0) = (+)_i Omega^(i)`, `Omega^(i)_(g,n)(e~_i, .., e~_i) = Delta_i(0)^((2g-2+n)/2) exp(Sum_k s^(i)_k ch_(2k-1)(E))`, `s^(i)_k = -(2k-2)!
Sum_(a in w_i) a^(1-2k)`; mixed colours give 0.

*Step 3 (Teleman v3, (8.6) and 8.7 Proposition (i), read verbatim, p.40).* For an odd series `h(z) = Sum_j h_(2j-1) z^(2j-1)`, the Hodge
CohFT `Zbar[h] := (n-th co-power of) alpha^g exp{Sum_j h_(2j-1) ((2j)!/B_2j) mu_(2j-1)}`, `mu_j = ch_j` of the Hodge bundle, satisfies
"(i) The theory `Zbar[h]` is the transform `T_z^(-1) o exp h(z) o T_z (I_A)` of the trivial A-theory", i.e. (PPZ15 Def. 2.13: `R.Omega = R T Omega`,
`T(z) = z 1 - z R^(-1)(z) 1`) `Zbar[h] = exp(h).omega`. Matching `h_(2k-1) (2k)!/B_2k = s^(i)_k` gives `h_(2k-1) = -b_k Sum_(a in w_i) a^(1-2k)`,
i.e. `exp h = h_(w_i)`. So `Omega^(i) = h_(w_i).omega^(i)` (rank one; `alpha^g` and the co-powers are the TFT `omega^(i)`). Sign cross-check: for
`c(E) = exp(Sum_(k odd) (k-1)! ch_k)` the same matching gives `exp(+Sum b_k z^(2k-1)) = h_({-1})`, the value LIT s.1.3 records from
Pandharipande 1712.02528 (not re-read here).

*Step 4 (Givental-Teleman for the formal quintic; GJR18 Theorem 4.1, read p.19).* `Omega = R.omega` (legs `R^(-1)(psi)`, edges `(sum e_a (x)
e^a - R^(-1) e_a (x) R^(-1) e^a)/(z_1 + z_2)`, vertices `T omega` with `T(z) = z(1 - R^(-1)(z) 1)`). For fixed `(g, n)` and cohomological degree
the graph sum is finite and polynomial in the coefficients of `R` and `omega`, so taking the `q^0` part commutes with it:
`Omega|_(q=0) = R(0).omega(0) = (+)_i c_i.omega^(i)` (Step 1; a diagonal `R` acting on a direct sum of rank-one TFTs keeps the colours: every
vertex term vanishes unless its half-edges have one colour, and the edge bivector is diagonal).

*Step 5 (uniqueness in rank one).* From Steps 2-4, `c_i.omega^(i) = h_(w_i).omega^(i)`; with `delta_i := c_i/h_(w_i) = exp(Sum_k d_k z^(2k-1))`
and the left group action (PPZ15 Prop. 2.14, read), `Zbar[log delta_i] = Zbar[0]`, i.e. `exp(Sum_k d_k ((2k)!/B_2k) mu_(2k-1)) = 1` on every
`Mbar_(g,n)`; comparing degrees, `d_k mu_(2k-1) = 0` for every `k`. Cited: Teleman 8.7 Proposition (ii) ("All cohomological Field theories
with flat vacuum based on a fixed semi-simple, pointed Frobenius manifold are classified by matrices `E o exp h(z)` in `Sp_L^+`, with
arbitrary `h` but the same `E`"), with Theorem (2) in (1.15) (p.11: "The CohFT's with underlying semi-simple Frobenius algebra A
constitute the `Sp_L^+ |x H^++`-orbit of the trivial theory `I_A`. The theories with flat vacuum form the orbit of the subgroup
`T_z o Sp_L^+ o T_z^(-1)`") and 8.10 Proposition (rank one, p.41): the classification datum is an invariant of the theory, so
`d_k = 0`. Honest note: "classified" (bijective) is how Teleman's 8.7 (ii) is read; the injectivity is not a separately numbered
statement. Elementary route for the only case needed (rank one): over `M_g`, `mu_(2k-1) = (B_2k/(2k)!) kappa_(2k-1)` (the Riemann-Roch
identity Teleman quotes on p.40), and `kappa_(2k-1) != 0` in `H^*(M_g)` for `g` large (Madsen-Weiss; quoted by Teleman, NOT read here); hence
`mu_(2k-1) != 0` and `d_k = 0`. Therefore `c_i = h_(w_i)`, i.e. `R(q = 0) = diag(h_(w_i))`. QED given the citations.

*Relation to the sources' quotes.* GJR18 Lemma 5.2 (read, p.27) asserts the Birkhoff factorisation `S_(i alpha) C_alpha^(-1) =
e^(u_alpha/z) R_(i alpha)` with (17) `C_alpha = h_(w_alpha)^(-1)` "by the results of [16, 12]" (Givental math/0008067, Coates-Givental), without
proof; Steps 1-5 prove `R(0) = diag(h_(w_i))` from the CohFT, given Teleman. *Numbering:* in the v3 PDF the paragraphs are "(1.15)" with
"Theorem (2)", "8.2 Proposition", "8.4 Proposition ([D, G1])" (solutions `R = Id mod z`, symplectic, unique up to right multiplication by
`exp(H_1 z + H_3 z^3 + ..)`, `H_(2i+1)` constant diagonal), "8.5 Corollary", "(8.6) Ambiguity for inhomogeneous theories", "8.7 Proposition",
"8.10 Proposition"; v3 has no heading "s.8.2". The project's "Prop. 8.5(ii) (s.8.2)" (MS5 s.2) is the ar5iv label of v3's 8.7 (ii).

## 4. Item 4: `lem:set-dilaton`: PROVED

*Statement.* Let `Omega` be a CohFT with unit `1` on `(V, eta)` (PPZ15 Def. 0.3 (i)-(iii)), `<prod_j tau_(k_j)(a_j)>_h := Int_(Mbar_(h,n))
Omega_(h,n)(a_1..a_n) prod_j psi_j^(k_j)` (`2h - 2 + n > 0`). (i) `<tau_1(1) prod_j tau_(k_j)(a_j)>_h = (2h - 2 + n) <prod_j tau_(k_j)(a_j)>_h`.
(ii) The model CohFT `Omega^mod_v` is a CohFT with unit `phi_0` (the flat unit `H^0`), and `<tau_1(phi_0)>_1 = 5/24`. (iii) Let `Lam` be a
complete local `Q`-algebra with maximal ideal `m` (e.g. `Q[[u]]`) and `kappa` in `1 + m`. With `Fbar_h(q) := Sum_n <(q + z 1)^n>_h/n!` (a
formal function of `y = q + z 1` in `V[[z]]`) and `d^n Fbar_h(q)[x_1..x_n] := Sum_m <x_1..x_n (q + z 1)^m>_h/m!`: for `2h - 2 + n > 0`,
`d^n Fbar_h(q/kappa) = kappa^(2h-2+n) d^n Fbar_h(q)`; for `n = 0`, `Fbar_h(kappa q) = kappa^(2-2h) Fbar_h(q)` (`h != 1`), `Fbar_1(kappa q) =
Fbar_1(q) - (5/24) log kappa` (model).

*Proof.* (i) PPZ15 Def. 0.3 (iii) (read): `Omega_(h,n+1)(a_1..a_n, 1) = p^* Omega_(h,n)(a_1..a_n)`. Arbarello-Cornalba alg-geom/9406008 (1.7)
(read): `pi_*(psi_1^(k_1)..psi_n^(k_n) psi_(n+1)^(a+1)) = psi_1^(k_1)..psi_n^(k_n) kappa_a` for the map forgetting the last point, the `psi_j`
on the left being those of `Mbar_(h,n+1)`, and `kappa_0 = 2h - 2 + n` (their p.5: "Notice that `kappa_0 = 2g - 2 + n`"). With `a = 0` and the
projection formula, `Int p^*Omega psi_(n+1) prod psi_j^(k_j) = Int Omega pi_*(psi_(n+1) prod psi_j^(k_j)) = (2h - 2 + n) Int Omega prod psi_j^(k_j)`.
(ii) PPZ15 Prop. 2.12 with Def. 2.13 (read, p.27-30): `R.Omega = R T Omega` is a CohFT with the same unit; the model TFT `omega^m_v`
(normalised canonical basis, `Deltahat_i^(1/2) = zeta^(4i) varsigma`) has unit `Sum_i Deltahat_i^(-1/2) e~_i`, which the isometry `Psi^m_v`
maps to `phi_0` (GAPS_G1_5 s.5.2, item (U)). `<tau_1(phi_0)>_1 = omega_(1,1)(phi_0) Int_(Mbar_(1,1)) psi_1`: only the degree-0 part of
`Omega_(1,1)` pairs with `psi_1`, and it is `omega_(1,1)` (`R = 1 + O(z)`, `T = O(z^2)`); `omega_(1,1)(1) = Sum_(a,b) eta^(ab) omega_(0,3)(1,
phi_a, phi_b) = dim V = 5` (gluing at the boundary point); `Int psi_1 = 1/24` (standard; not re-read). Only the genus-1 log term uses it.
(iii) Write `q/kappa + z 1 = y/kappa + eps z 1`, `eps := 1 - 1/kappa` in `m`. Expanding and applying (i) `b` times (allowed: `2h - 2 + n + a >
0` at every step), `<x, y^a, (z1)^b>_h = (N)_b <x, y^a>_h`, `N := 2h - 2 + n + a`, `(N)_b := N(N+1)..(N+b-1)`. So `d^n Fbar_h(q/kappa)[x] =
Sum_a (kappa^(-a)/a!) Sum_b (eps^b/b!) (N)_b <x, y^a>_h = Sum_a (kappa^(-a)/a!) (1 - eps)^(-N) <x, y^a>_h = kappa^(2h-2+n) d^n Fbar_h(q)[x]`
(binomial series, `m`-adically convergent; `(1 - eps)^(-1) = kappa`). The `n = 0` statements follow the same way; for `h = 0` the terms with
fewer than three `y` vanish (`psi = 0` on `Mbar_(0,3)`); for `h = 1`, `<(z1)^m>_1 = (m-1)! <tau_1(1)>_1` gives `Sum_(m>=1) (-eps')^m (m-1)!/m!
(5/24) = -(5/24) log kappa` (`kappa = 1 + eps'`). QED

*Hypothesis check for the chain.* The lemma is applied with `kappa = d = (P_0)_33` (S15 s.1.1-1.2; `lem:tra-F`, `lem:flh-dil`). The
formal statement needs `d` in `1 + uQ[[u]]`; the sources give `d = tau^(1/4) varpi_0` (S15 s.2 remark, line 90), so `d(0) = 1` in the
programme family (`varpi_0(0) = 1`, `tau(0) = 1`). The writers of s04/s07 should state `d(0) = 1` where the lemma is invoked (for a
constant `kappa != 1` the expansion at `q/kappa` is not a formal identity).
