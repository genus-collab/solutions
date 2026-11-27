# Gap cluster G1, part 3 of 5: item 2, `prop:set-cyclicLR` (proof, corollaries, literature, evidence)

Continuation of `GAPS_G1_2.md` (Lemmas A-F, notation there).

## 2.8 Proof of the Proposition

1. `kappa_0 = kappa_2 = 1`, `kappa_3 = kappa_4` (Lemma B).
2. *`kappa_1`.* By Lemma D, `c_k = C_0 phi_k` with `phi_0 = 1`, `phi_1 = T`. For any functions, `det(D^j (C_0 f_k)) = C_0^4 det(D^j f_k)`
   (`n = 4` functions), and `D^j = C_1^j d_T^j + (terms with fewer d_T)` (since `D = C_1 d_T`), so `det(D^j f_k) = C_1^(0+1+2+3)
   det(d_T^j f_k)`. With `f = (1, T, phi_2, phi_3)` and `d_T T = 1`:
   `det(d_T^j f_k) = d_T^2 phi_2 d_T^3 phi_3 - d_T^3 phi_2 d_T^2 phi_3 = (d_T^2 phi_2)^2 d_T(d_T^2 phi_3 / d_T^2 phi_2) = kappa_1^2 kappa_2`
   (Lemma C). Lemma E gives `C_0^4 C_1^6 kappa_1^2 kappa_2 = L^10`; with `kappa_2 = 1`, `(kappa_1 - L^5/(C_0^2 C_1^3))(kappa_1 + L^5/(C_0^2 C_1^3))
   = 0` in the domain `Q[[q]]`; the second factor has constant term 2, so `kappa_1 = L^5/(C_0^2 C_1^3)`.
3. *`kappa_3`.* Let `Mop := d_T^2 o kappa_1^(-1) o d_T^2 o C_0^(-1)` (an operator with coefficients in `Q[[q]]`; `C_0, C_1, kappa_1` are
   units). By Lemma C: `Mop c_0 = d_T^2 kappa_1^(-1) d_T^2 1 = 0`; `Mop c_1 = 0` (`d_T^2 T = 0`); `Mop c_2 = d_T^2(kappa_1^(-1) kappa_1) = 0`;
   `Mop c_3 = d_T(kappa_2) = 0`; `Mop c_4 = kappa_3`. In powers of `D`, `Mop = (kappa_1 C_0 C_1^4)^(-1) D^4 + (lower)`, and `PF4 = L^(-5) D^4 + ..`.
   So `Mop - g PF4`, `g := L^5/(kappa_1 C_0 C_1^4)`, has order `<= 3` and kills `c_0..c_3` (Lemma D); by Lemma F it is 0. Hence
   `kappa_3 = Mop c_4 = g PF4 c_4 = g` (Lemma D), and by step 2 `g = C_0^2 C_1^3 / (C_0 C_1^4) = C_0/C_1`.
4. `kappa_4 = kappa_3 = C_0/C_1` (Lemma B).
5. *Equivalent forms.* For a solution `v` of `z d_T v = H * v`, `z d_T eta(H^a, v) = eta(H^a, H * v) = eta(H * H^a, v) = kappa_a eta(H^(a+1), v)`
   (Frobenius), i.e. `z D v_a = C_1 kappa_a v_(a+1)`, `a` in Z/5. For `v = R(z) e~_i e^(u_i/z)` (`def:set-rows`: Givental's convention, the
   columns `R e~_i e^(u_i/z)` solve the QDE) and `D u_i = zeta^i L` (Corollary 2.10 (a)), `v_a = e^(u_i/z) Q^(i)_a` and `z D (e^(u_i/z) f) =
   e^(u_i/z)(z D + zeta^i L) f`, which is the row recursion with `X_a := zeta^(-i) Q^(i)_a` (the sources' `Q = zeta^i X`). QED

Everything is an identity in `Q[[q]]` or `Q[[q]][ell]`: no truncation in `q`, no analytic continuation.

## 2.9 What the proof uses, and what it does not

* Uses: (M); (G0-a)-(G0-d); `lem:set-pf5`; the Frobenius property. Only the fixed-point component structure of `I` at one colour is
  used in Lemma D (any `i` works).
* Does not use: Zagier-Zinger (the identities `I_(0,0)..I_(4,4) = L^5`, `I_(p,p) = I_(4-p,4-p)` of the literature are re-proved: s.2.11);
  the closing relation `4 -> 0` (it follows; Corollary 2.10 (b)); (CYC) or any conifold-frame statement; any finite computation.
* Remark (why `kappa_3` is "genuinely equivariant", R1): `kappa_3 = Mop c_4` involves the period `c_4`, which satisfies the INHOMOGENEOUS
  equation `PF4 c_4 = 1`; in the non-equivariant quintic `H^4 = 0` and `c_4` does not occur. The coefficients `c_0..c_4` coincide with the
  non-equivariant ones (`(w + k)^5 - w^5 = (w + k)^5 mod w^5`), so `kappa_0..kappa_3` are the same series as in the non-equivariant
  Birkhoff recursion; only the closing (`H^5 = 1`) is new (s.2.12, control c1).

## 2.10 Corollaries (all PROVED from the Proposition; they close further implicit steps of `def:set-genus0`, `def:set-rows`)

(a) *Eigenvalues, idempotents, `Delta_i`, canonical coordinates.* `(H*)^5 = (prod_a kappa_a) Id = (L/C_1)^5 Id`, so the eigenvalues of `H*`
are `mu_i = zeta^i L/C_1` (the branch with `mu_i(q=0) = zeta^i`, the classical eigenvalue at `p_i`). The idempotents are `e_i =
(1/5) Sum_a (prod_(b<a) kappa_b) mu_i^(-a) H^a` (eigenvector equation `x_a = x_(a-1) kappa_(a-1)/mu_i`; `Sum_i e_i = 1` forces `x_0 = 1/5`).
`eta(e_i, e_i) = eta(1, e_i) = (1/5) kappa_0 kappa_1 kappa_2 mu_i^(-3) eta(1, H^3) = kappa_1 mu_i^(-3) = L^2/(C_0^2 zeta^(3i))`, so `Delta_i
= C_0^2 zeta^(3i)/L^2` (`def:set-frob`). The canonical coordinates satisfy `d_T u_i = mu_i` (eigenvalue of `H*` on `e_i`), i.e. `D u_i =
zeta^i L` (`def:set-genus0`; additive constants as fixed in `lem:set-canon`).

(b) *Large-radius cyclic system closes on PF5.* For every solution `v` of `z d_T v = H * v` (in any differential `Q[[q]]`-algebra, e.g.
formal WKB columns `e^(u_i/z) Q[[q]]((z))`), `C_0 v_0` solves PF5. Proof: from `z d_T v_a = kappa_a v_(a+1)` (`a = 0..3`), `v_4 = z^4 kappa_3^(-1)
Mop(C_0 v_0)`; by step 3, `kappa_3^(-1) Mop = kappa_3^(-1) g PF4 = PF4` (`g = kappa_3`); the closing `z d_T v_4 = kappa_4 v_0` becomes `z^5 C_1^(-1)
D PF4(C_0 v_0) = (C_0/C_1) v_0`, i.e. `(z^5 D PF4 - 1)(C_0 v_0) = 0`, which is PF5 (`lem:set-pf5`). Applied to `v = R e~_i e^(u_i/z)`: the
function `e^(u_i/z) z C_0 Q^(i)_0 = zeta^i e^(u_i/z) z L Phi_i(zhat, u)` solves PF5. This is the "WKB form of colour `i`" of `def:set-rows`,
now derived from the definition `Phi_i := C_0 Q^(i)_0/(zeta^i L)` (recommended as the definition in the note; GAPS_G1_4 s.3.1).

(c) *Membership in the differential ring (R1 fix item 3).* By the Proposition the recursion of the rows has coefficients `(C_1 kappa_a)^(-1)
= (1/C_1, C_0^2 C_1^2 L^(-5), 1/C_1, 1/C_0)` (`a = 0..3`) and `D L = (L^6 - L)/5`; so once `A^act_n` is known to lie in `Q[L, 1/L]`
(`thm:gau-R0` (iii), s03) all rows lie in the ring generated by `D^j C_0`, `D^j Itilde_1` (`j <= 3`, `D^4` reduced by PF4), `L^(+-1)`,
`C_0^(-1)`, `C_1^(-1)`. This is the forward dependency R1 asks to make explicit (`cor:gau-rowsNN`); nothing more is needed from G1.

## 2.11 The statement in the literature (FOUND; corroboration, not the proof of record)

* **GJR18 (1812.11908v1) s.5.1 eq. (13)** (read): "The matrix of quantum multiplication by `tau-dot = H + q d/dq(I_1/I_0)` in the flat basis
  `{H^i}` is `A` = the transpose of the cyclic matrix with entries `I_(1,1), I_(2,2), I_(3,3), I_(4,4)` and `I_(5,5) lambda^5`", citing
  [21, s.6] (Guo-Ross 1703.06955) and [20, s.7] (Guo-Ross 1611.08876, not read). (`tau = (I_1/I_0) H` there, so `tau-dot` must be read
  `q d/dq(I_1/I_0) H = C_1 H`; the printed "`H + ..`" is a misprint.) With `I_(a+1,a+1) = C_1 kappa_a` this is the Proposition at `lambda = 1`.
* **Guo-Ross 1703.06955v1 s.6.2** (read): (16) quantum product; (17) `I_(0,q) := I_q`, `I_(p,q) := D(I_(p-1,q)/I_(p-1,p-1))`, so that
  `M^p(I/z) = Sum_q I_(p,p+q) z^(-q) Phi_(p+q)`, `M F := z D (F/F(x, oo))`; **Prop. 6.1**: `S(x,z)^*(Phi_m) = M^m(I/z)/I_(m,m)` ("This
  follows from standard properties of Givental's Lagrangian cone"); **Prop. 6.2** ("[ZZ08] Theorem 2, [GR16] Lemma 7.6"): (1) `I_(0,0)
  .. I_(4,4) = L^5`; (2) `I_(5+p,5+p) = lambda^5 I_(p,p)`; (3) `I_(p,p) = I_(4-p,4-p)` (`0 <= p <= 4`). Dictionary: `I_(m,m) = C_1 kappa_(m-1)`
  (`m >= 1`), `I_(0,0) = C_0`: (3) gives `kappa_0 = kappa_2`, `kappa_3 = C_0/C_1`; (1) gives `kappa_1`; (2) gives `kappa_4 = C_0/C_1`.
* **LP18 (1803.01409v2) s.2.4.3-2.4.4** (read): (17) `S(H^j) = z(d/dt)S(H^(j-1)) / [z (d/dt) S(H^(j-1))]|_(t=0,H=1,z=oo)` (stable-quotient side,
  "Birkhoff factorization, see [18]"), (19)-(20) `C_0 := I|_(z=oo)`, `C_i := z(d/dt) S(H^(i-1))|_(z=oo,t=0,H=1)`, and "The following relations
  were proven in [27] [= ZZ08]: `C_0 C_1 C_2 C_3 C_4 = L^5`, `C_i = C_(4-i)`". With `C_(a+1) = C_1 kappa_a` (a = 0..3) this is
  `kappa_0..kappa_3`; the closing is not stated there.
* Zagier-Zinger arXiv:0710.0889 itself was NOT read; its relations concern the non-equivariant series, whose first five Birkhoff
  coefficients coincide with ours (s.2.9). The proof of s.2.8 makes the citation unnecessary.

## 2.12 Evidence (VERIFIED; never a proof step)

`g1_cyclic_check.py` (this folder; exact rationals; log `g1_cyclic_check.log`): with `Fhat` of Lemma D and the recursion of Lemma C in
the `w`-form (`g_0 = Fhat/C_0`, `I_(p+1) = 1 + D [w^1] g_p`, `g_(p+1) = (1 + D/w) g_p / I_(p+1)`, `I_(p+1) = C_1 kappa_p`), to `O(q^31)`:
`I_1 = C_1`, `I_2 C_0^2 C_1^2 = L^5`, `I_3 = C_1`, `I_4 = C_0`, `I_5 = C_0` and the closing `g_5 = g_0` through `w^7`; the ingredients of
Lemmas D and E in `Q[[q]][ell]`: `PF4 c_k = 0` (`k <= 3`), `PF4 c_4 = 1`, `det(D^j c_k) = L^10` (no `ell`): all PASS. Failing
controls, all fail as required: (c1) the non-equivariant series (`(w + k)^5` in the denominator): `I_1..I_4` unchanged, closing fails;
(c2) twist `O(4)` or `O(6)` instead of `O(5)`: the `kappa_1` identity fails; (c3) the guess `kappa_3 = C_1/C_0` fails. The project's earlier
check (FQ2 s.1, s.5: `q^12`) is superseded by the proof.
