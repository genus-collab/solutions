# Gap closing, cluster G5, part 3: the natural factorisation (`cor:mai-natfact`) and the degree bound (`prop:con-deg`)

Continuation of `GAPS_G5.md`, `GAPS_G5_2.md`. Same rules.

## 3. Item 3: `cor:mai-natfact` (the natural factorisation and `Psi(u)`)

### 3.1 What is stated without proof
In `OUTLINE_8.md`: the definition of `Psi(u)` (`\draftnote{the definition of Psi(u) (DA s.7 item 7) to be checked against UR s.4 and OU3 s.3}`);
UR s.4's Corollary F is one line ("PROVED given (H)") with `Psi(u)` undefined; OU3 s.3 calls `Psi` "the constant orthonormal change of
basis flat -> colour" with the relation `Psi eta^-1 Psi^T = eta^-1`; OU3 s.6d says `Psi(u)` is "orthogonal for every `u`"; RN s.3(d)
states the symplecticity of `G` in one line; (b) uses `p_(0,a)(0) = 1` (the one line G8 of `lem:gau-diag` step 4) and (R1) without saying
so. Code: `natural_form/referee/ref_rows.py` (docstring) defines the Teleman form `Rt_ji = (1/5) Sum_(a+b = 3 mod 5) Qh_(b,0)^(j) Qh_a^(i)`,
i.e. `R~ = Qhat(u,0)^T eta^-1 Qhat`: this fixes the definition (DA s.7 item 7 read it so).

### 3.2 Definitions and two lemmas (all PROVED here)
* `Qhat(u, zhat)`: the hatted lowered rows of `def:set-hat` (flat index `a`, colour `i`), `Qhat_a = s^(1-nu_a) Q_a`, `Q` expanded in the frame.
* **`Psi(u) := Qhat(u, zhat = 0)`**, the `zhat^0` coefficient (it does not depend on the normalisation series `D`, since `D(0) = 1`).
* **`R~(zhat; u) := Psi(u)^-1 Qhat(u, zhat)`**, the Teleman (colour-colour) form.
* *Lemma 3.A (explicit `Psi`; orthonormality).* `Psi(u)_(a j) = zeta^(j(a+1)) w_a(u)` with `w = (1/varpi_0, 1/(varpi_0 D t_c), varpi_0 D t_c, varpi_0, 1)`;
  hence `Psi(u)` lies in `GL_5(Q(zeta)[[u]])`, `Psi(0) = Psi_0` (entries `eps_a zeta^(j(a+1))`, `def:set-M0`), and `Psi(u)^T eta^-1 Psi(u) = 1`,
  so `Psi^-1 = Psi^T eta^-1` and `R~ = Psi^T eta^-1 Qhat` (the code's formula).
  *Proof.* At `z^0` the WKB amplitude of every colour is `L` (`ph_0 = 1`, `def:set-rows`), and the forward relations of `prop:set-cyclicLR`
  reduce to `X_(a+1)(0) = zeta^j L X_a(0)/(C_1 kappa_a)`, `X_0(0) = L/C_0`, `Q^(j)_a = zeta^j X_a`. With `C_1 kappa = (C_1, L^5/(C_0^2 C_1^2), C_1, C_0, C_0)`:
  `Q^(j)_a(0) = zeta^(j(a+1)) m_a`, `m = (L/C_0, L^2/(C_0 C_1), C_0 C_1/L^2, C_0/L, 1)`. Expanding in the frame (`C_0 -> varpi_0`, `C_1 -> D t_c`, `L -> 1/s`)
  and hatting (`s^(1-nu_a)`, `nu = (0,-1,3,2,1)`) gives `w_a`; at `u = 0`, `varpi_0 = 1`, `D t_c = -1` give `w(0) = eps`. Orthonormality: `eta^-1_ab = (1/5)
  delta_(a+b = 3)`, and `m_0 m_3 = m_1 m_2 = m_4^2 = 1`, so `(Psi^T eta^-1 Psi)_ij = (1/5) Sum_a zeta^(i(a+1) + j(4-a)) = (1/5) zeta^(i-j) Sum_a
  zeta^((i-j)a) = delta_ij` (equivalently: `e~_i` is orthonormal, `def:set-frob`, and hatting preserves `eta`, `nu_a + nu_b = 2` on `eta`-pairs). QED
* *Lemma 3.B (the hatted gauge).* For the graded canonical gauge `P_can = Sum_k z^k Sum_b p_(k,b)(u) E_(b+k, b)` (`p_(k,b)` in `Q[[u]]`, `thm:gau-U`,
  `prop:gau-graded`, `def:gau-can`) put `Ghat := Lambda P_can(u, s^4 zhat) Lambda'^-1`, `Lambda = diag(s^(1-nu_a))`, `Lambda' = diag(varsigma^(1-nu_b))`.
  Then `Ghat_ab = tau^((1-nu_b)/4) Sum_k zhat^k u^(mu(b,k)) p_(k,b)(u) delta_(a, b+k)`, `mu(b,k) = (4k + nu_b - nu_(b+k))/5`, an integer with
  `mu(b,k) >= max(0, (4k-4)/5)`; `mu(b,k) = 0` only for `k = 0` and for `(b,k) = (1,1)` (entry `(2,1)`). Hence `[u^n] Ghat` is a polynomial in `zhat`
  of degree `<= floor((5n+4)/4)`. *Proof.* `varsigma = s tau^(-1/4)` (`def:set-tau`); the `s`-exponent of the `(b+k, b)` entry is
  `1 - nu_(b+k) + 4k - (1 - nu_b) = 4k + nu_b - nu_(b+k)`; `nu_a = -a mod 5` makes it divisible by 5; `nu` takes values in `{-1,..,3}`, so it is
  `>= 4k - 4`, and `>= 0` (`k = 0`: 0; `k >= 1`: `>= 0`); it vanishes for `k >= 1` only if `nu_(b+k) - nu_b = 4`, i.e. `k = 1`, `b = 1`. QED
  (VERIFIED for `k <= 40` by `g5_natfact_checks.py`; not needed.)

### 3.3 Statement and proof (proposed text of `cor:mai-natfact`)
Assume `D = 1` (`thm:gam-D`). In every frame of `def:set-frames`:
* (a) `Qhat(u, zhat) = Ghat(u, zhat) Rhat^mod(tau(u) zhat)`; `Ghat` is a power series in `u` whose coefficients are polynomials in `zhat` (3.B).
* (b) `Qhat(0, zhat) = Rhat^mod(zhat)`: the conifold symbol is the model's; so `M_0 = R~(zhat; 0)` (NF's definition), and `Psi(0) = Psi_0`.
* (c) `R~(zhat; u) = G(u, zhat) M_0(tau(u) zhat)` with `G := Psi(u)^-1 Ghat Psi_0 = Psi_0^-1 Ghat(u,0)^-1 Ghat(u,zhat) Psi_0`; `G(u, 0) = G(0, zhat) = 1`;
  `[u^n] G` is a polynomial in `zhat`; `G(-zhat)^T G(zhat) = 1`; `M_0(tau zhat)` is the same in every frame, `G` depends on the frame.
*Proof.* (a) `thm:gau-R0` (iv) with `D = 1` gives `Q = P_can Q^m(v(u))` (expanded rows) with `Q^m_b = varsigma^(nu_b - 1) Rhat^mod_b(z/varsigma^4)`,
`z/varsigma^4 = tau zhat` (`prop:set-modrows`, `def:set-tau`, `lem:set-canon`); multiply by `Lambda` and insert `Lambda'^-1 Lambda'`: `Qhat = Ghat Rhat^mod(tau zhat)`.
Polynomiality is Lemma 3.B.
(b) At `zhat = 0`, (a) gives `Psi(u) = Ghat(u, 0) Psi_0` (`Rhat^mod(0) = Psi_0`, `prop:set-modrows`: `Fs_b(0) = eps_b`). With Lemma 3.A, `Ghat(0,0) =
Psi(0) Psi_0^-1 = 1`. By Lemma 3.B the only `u^0` entries of `Ghat` beyond the diagonal sit at `(2,1)` with coefficient `zhat p_(1,1)(0)`, and
`p_(1,1)(0) = 0` by (R1) (`lem:gau-R1`: column 1 of `P(0, z)` is `P_11(0, z) e_1`). So `Ghat(0, zhat) = 1`, and `tau(0) = 1` gives `Qhat(0, zhat) =
Rhat^mod(zhat)`. Then `R~(zhat; 0) = Psi_0^-1 Rhat^mod(zhat) = M_0(zhat)` (`def:set-M0`).
(c) `R~ = Psi(u)^-1 Ghat Rhat^mod(tau zhat) = [Psi(u)^-1 Ghat Psi_0][Psi_0^-1 Rhat^mod(tau zhat)] = G M_0(tau zhat)`; the second form of `G` uses
`Psi(u) = Ghat(u,0) Psi_0`; `Ghat(u,0) = diag(tau^((1-nu_b)/4) p_(0,b))` is invertible over `Q[[u]]` (`p_(0,b)` units), so `[u^n] G` is a
polynomial; `G(u,0) = 1` is clear and `G(0, zhat) = 1` is (b). Symplecticity: `Cf(Q) = P_can Q^m` with `P_can(-z)^T eta P_can(z) = eta`
(`prop:gau-symp`; `eta^-1 = eta/25`) and the model rows symplectic (`M_0(-t)^T M_0(t) = 1`, `def:set-M0`; `Psi_0` orthonormal) give
`Qhat(-zhat)^T eta^-1 Qhat(zhat) = 1` (hatting preserves `eta^-1`), hence `R~(-zhat)^T R~(zhat) = 1` (Lemma 3.A: `Psi^-T Psi^-1 = eta^-1`), and
`G(-zhat)^T G(zhat) = M_0(-tau zhat)^-T M_0(tau zhat)^-1 = 1`. Frame independence of `M_0(tau zhat)`: `M_0` and `tau` do not involve the frame. QED
*Remarks.* (i) The line G8 of `lem:gau-diag` step 4 (`p_(0,a)(0) = 1` for `P_can`) follows in s03 directly: `z^0` of `z P' = A P - P A~` gives
`A P_0 = P_0 A~`, so `p_(0,a+1)/p_(0,a) = b~_a/b_a`; at `u = 0` all these ratios are 1 (`b_a(0) = b~_a(0)`, namely `1, 1, -1, -1` for
`a = 0, 2, 3, 4`, by `t_c'(0) = varpi_0(0) = v'(0) = 1`; `lim u b_1 = -1/(varpi_0^2 (D t_c)^2) = -1 = lim u b~_1`), so all `p_(0,a)(0)` equal
`p_(0,1)(0) = [z^0] P_11(0,z) = 1`. After `D = 1` it
also follows from Lemma 3.A: `Ghat(u,0) = Psi(u) Psi_0^-1 = diag(1/varpi_0, -1/(varpi_0 D t_c), -varpi_0 D t_c, varpi_0, 1)` (consistent with
`Ghat(u,0)_00 = 1/varpi_0`, `g_0 g_3 = 1`, `prod = 1` of `lem:gau-diag`). (ii) (H) and the canonical normalisation are not hypotheses:
they are `D = 1`. (iii) "`M_0` on the left fails" is computed only (RN s.3(c)): appendix C. (iv) Not part of the corollary: the uniqueness
of `tau` (`rmk:mai-tau`) and the degree bound of `[u^n] G` (`rmk:mai-degree`, only the a priori bound of 3.B is proved).

### 3.4 Erratum (for `design/SOURCE_ERRATA.md`; not applied by me)
OU3 s.3, parenthesis "`Psi eta^-1 Psi^T = eta^-1`": with `Psi = Qhat(u,0)` (flat index first) this is false (`Psi_0 eta^-1 Psi_0^T =
diag(zeta^(3(a+1)))`, `g5_natfact_checks.py`); the correct relation is `Psi^T eta^-1 Psi = 1` (equivalently `Psi Psi^T = eta`). The
consequence drawn in OU3 ("it drops out of every edge") is right with the correct relation: `R'(-psi)^T eta^-1 R'(-psi') = G(-psi)^T G(-psi')`.

## 4. Item 4: `prop:con-deg` step 5 (the triviality of the resonance)

### 4.1 Source and what is missing
ER s.4 (stage4/evid/referee/REFEREE.md, "Referee result B"): the lemma on `Pi`, the map `P = (x^-5 - 1) Pi + M^5 - 1` of `x Q[x][z]`, the
induction on the lowest power of `psi_n`; step 5 is the single sentence "the right-hand side has top `s`-exponent `<= -1`, so W stops at
`s^-6` and never reaches `s^-5`". No other HVQ source expands it (searched all `allgenus/*.md` for "resonan": EVID s.2.1 (b') has
"`deg rho_n < 4n/5` when `5 | n`, so every resonance is unobstructed for a trivial reason" as VERIFIED only, `n <= 1000`; EVID (b) and
POLY s.2 obtain the resonances from Theorem U or from P(a); PR s.2 confirms; the Stage-4 syntheses quote ER). Below, step 5 is written
via a factorisation that ER does not state; it proves EVID's (b') and makes step 4 independent of `lem:gau-poly`.

### 4.2 Set-up (colour 0; colour `i` is colour 0 in `x' = zeta^i x`, since `M_i(x) = M_0(zeta^i x)`, `x'^5 = x^5`, `theta_(x') = theta_x`)
`x := L = 1/s`; `D = -(x/5)(1 - x^5) d/dx`, so the conjugated operator `M := zD + x = x - (z/5)(1 - x^5) theta_x` (`theta_x = x d/dx`);
`q = (1 - x^-5)/5^5`; `P := e^(-u_0/z) PF5 e^(u_0/z) = M^5 - 1 + (x^-5 - 1) Pi`, `Pi := 5^-5 prod_(k=1..5)(5M + kz)`; `P = Sum_(j=0..5) z^j P_j`.
`psi = Sum_n A_n z^n`, `A_n = x^(4n+1) ph_n(x^-5)` (`def:set-rows`; `A_0 = x`). `[zhat^(n+1)] HPF = 0` is `[z^(n+1)] P psi = 0` (`lem:gau-HPF`).
* (L0) `P_0 = 0`: at `z^0`, `M = x`, `Pi = x^5`, and `x^5 - 1 + (x^-5 - 1) x^5 = 0`.
* (L1) For `j >= 1`, `P_j = (1 - x^5) o P~_j` with `P~_j` an operator on `Q[x, x^-1]`. Indeed `M = x + z N`, `N := -(1/5)(1 - x^5) theta_x`, and
  `N o (1 - x^5) = (1 - x^5) o (N + x^5)` (as `theta_x(1 - x^5) = -5 x^5`), so `M o (1 - x^5) = (1 - x^5) o (M + z x^5)`: every `z^j`-term (`j >= 1`)
  of `M^5` has a factor `N`, hence a left factor `1 - x^5`; and `(x^-5 - 1) Pi = x^-5 (1 - x^5) Pi`; the `z^0` terms cancel by (L0). (This is the
  flat form of "every raising carries a LEFT factor `(1-u)`" of POLY s.1 and EVID s.2.1 (b): `1 - x^5 = -x^5 (1 - u)`.)
* (L2) (ER's lemma) `Pi(x^m)` lies in `x^6 Q[x][z]` for `m >= 1`: `(5M + kz) x^m = 5 x^(m+1) + z(k - m) x^m + z m x^(m+5)`; the factors commute;
  apply `k = m, m+1, .., 5` first (each kills the lowest power when it equals `k`); no factor lowers the valuation. Hence `P_j(x Q[x]) ⊂ x Q[x]`.
* (L3) `P_1 x^m = (1 - m)(x^(m-1) - x^(m+4))`, i.e. `P_1 = (1 - x^5) o E`, `E x^m := (1 - m) x^(m-1)` (ER: `T_1`; colour `i`: a factor `zeta^-i`).
  Direct: the `z^1` part of `M^5 x^m` is `-(2 + m)(1 - x^5) x^(m+4)`, that of `Pi x^m` is the same plus `3 x^(m+4)` (`Sum k/5 = 3`), and
  `-(2+m)(1 - x^5) x^(m+4) + x^-5 (1 - x^5)[-(2+m)(1 - x^5) x^(m+4) + 3 x^(m+4)]` collapses to `(1 - x^5) x^(m-1) [3 - (2 + m)]`.

### 4.3 Proof of `prop:con-deg` with step 5 written out (PROVED)
Fix free coefficients `frak f_m` (arbitrary). *Claim:* for every `n` there is a unique `A_n` in `x^(4n+1) Q[[x^-5]]` with `[z^(n+1)] P psi = 0` and
`[x^1] A_(5m) = frak f_m`, and it lies in `x Q[x]`. Induction on `n`; `A_0 = x`. Assume `A_0, .., A_(n-1)` in `x Q[x]`. The equation is
`P_1 A_n = RHS_n := -Sum_(j=2..5) P_j A_(n+1-j)`.
* (Step 4, ER) `RHS_n` lies in `x Q[x]` by (L2).
* (**Step 5**) `RHS_n` lies in `(1 - x^5) Q[x, x^-1]` by (L1). As `1 - x^5` is a unit of `Q[[x]]`, `rho_n := RHS_n/(1 - x^5)` lies in
  `x Q[[x]] ∩ Q[x, x^-1] = x Q[x]`. By (L3) and since `1 - x^5` is not a zero divisor, the equation is `E A_n = rho_n`, i.e.
  `(1 - m)[x^m] A_n = [x^(m-1)] rho_n` for every `m`. At `m = 1` it reads `0 = [x^0] rho_n`, which HOLDS because `rho_n` lies in `x Q[x]`: the
  resonance is trivial. (Under the dictionary of `lem:gau-HPF`, `rho_n = x^(4n) r_n(x^-5)`, so `[x^0] rho_(5m) = [u^(4m)] r_(5m)` is exactly the
  resonance coefficient of `lem:gau-HPF`, and `rho_n` in `x Q[x]` is `deg r_n < 4n/5`.) For `m != 1`, `[x^m] A_n = [x^(m-1)] rho_n/(1 - m)`: zero
  for `m <= 0`, determined for `m >= 2`; `[x^1] A_n` is free.
* (Grading) `P` maps `x^m z^j` to terms `x^(m') z^(j')` with `m' + j' = m + j mod 5`; `psi` has grade 1, so `A_n` has only exponents `= 4n + 1 mod 5`
  (consistent with `x^(4n+1) Q[[x^-5]]`): `[x^1] A_n` is allowed only for `5 | n`, where it is `frak f_(n/5)` (`x^1 = x^(4n+1) u^(4n/5)`).
* Hence `A_n` lies in `x Q[x]`, and it is the unique solution in `x^(4n+1) Q[[x^-5]]` (the equation `E A_n = rho_n` there has kernel `Q x`).
* (Step 6) `A_n` in `x Q[x]` iff every `u^k` in `ph_n` has `4n + 1 - 5k >= 1`, i.e. `deg ph_n <= floor(4n/5)`; `A_n/x` at `x = 0` is `frak f_(n/5)`
  (`5 | n`) or 0: the value at `L = 0` (the orbifold point).
Consequences: every power-series solution of HPF with `ph_0 = 1` is this one for its free coefficients (uniqueness of the `u`-recursion), so
the degree bound holds for every colour and every choice; and, independently of `thm:con-Pa` and of Theorem U, HPF has a polynomial
solution for every choice of the free coefficients (the hypothesis of `lem:gau-poly` (6), "the resonances hold for every choice"). QED

### 4.4 For the writer of s08
Replace the `\draftnote` of step 5 by the bullet "Step 5" (with (L1)); state (L1) as part of step 3; step 4 no longer needs `lem:gau-poly`
(the induction constructs `A_n` in `x Q[x]`). Keep: equality `deg ph_n = floor(4n/5)` is not claimed. Evidence (appendix C only):
`g5_degree_bound.py`, s.5.
