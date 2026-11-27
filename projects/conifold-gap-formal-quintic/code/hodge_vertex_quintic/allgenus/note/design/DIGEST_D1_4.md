# Digest D1, part 4: s.3 proofs

**Piece map (final; it supersedes the list in the header of `DIGEST_D1.md`, written before the split):** `DIGEST_D1.md` s.0-s.1;
`_2.md` s.2.0-s.2.3; `_3.md` s.2.4-s.2.5; this file s.3; `_5.md` s.4-s.7.

For each statement of s.2: the proof as numbered steps, each with where it is written. Markers: **[W]** written out in the cited
place; **[T]** written tersely (one sentence or less), the writer expands it from the same place; **[V]** only VERIFIED; **[M]** missing
(needs a `\draftnote`). No step below is [V] or [M] unless marked.

## 3.1 Polynomial structure

* **`lem:pol-HPF`.** (1) Chain rule at fixed `z` (`u d/du = th_u - (4/5) th_w`): on `F = s^e Fh(w, u)`, `(zD + L) F = s^(e-1) A_e Fh` [W:
  POLY s.1; PR s.2; PB s.2.2]. (2) `z^m M^(5-m) psi = s^-1 u^(m-1) w^m Phi_(5-m)` for `psi = s^-1 Phi` [W: POLY s.1]. (3) Insert into the
  `u`-form of PF5 (D1.2) and divide by `s^-1`: HPF [W: POLY s.1]. (4) `[w^m] Phi_j = [w^m] Phi_(j-1) + (1-u)((4(m-1)+j)/5 - th_u)[w^(m-1)] Phi_(j-1)`,
  so every raising carries a left factor `(1-u)` [W: POLY s.1; EVID s.2.1 (b)]. (5) The `ph_n`-part of `[w^(n+1)] HPF` sums to
  `(4n + 3 - 5 th_u) ph_n - 3 ph_n`; the rest is `r_n` (displayed in POLY s.1, with `Phi^b_(j,n)` computed with `ph_n := 0`) [W: POLY s.1].
* **`lem:pol-u0`.** (1) The `u^0` part of HPF involves only `ph_k(0)`: `th_u` kills constants and the bracket terms carry `u^k`, `k >= 1`
  [W: POLY s.2 Lemma 1; GR s.2.1]. (2) The constant term of `(4n - 5 th_u) ph_n` is `4n ph_n(0)`, pivot `4n != 0` [W: POLY s.2]. (3) This
  `u^0` recursion is that of the scaling limit `v th^5 + 3z th^4 + 1` (UR s.1), whose homogeneous WKB solution with leading term 1 is
  unique and is `F_0` (Prop. S) [W: ER s.2, the line ER s.7 item 5 asks to add]. (4) The free coefficient sits at `u^(4m)`, `m >= 1` [W: EVID
  s.2.1 (a)].
* **`lem:pol-poly`.** (1) `4n - 5 th_u` is diagonal on `u^k` with eigenvalue `4n - 5k` [W: POLY s.2; PB s.2.2]. (2) `r_n` is a polynomial if
  `ph_0..ph_(n-1)` are [W: POLY s.2]. (3) Hence `[u^k] ph_n = [u^k] r_n/(5k - 4n)` vanishes for large `k`, except the free coefficient at
  `5k = 4n` [W: POLY s.1-s.2]. (4) Multiplication by a function of `z` commutes with HPF, since `((4/5) th_w - th_u) z^m = 0` [W: PR s.2].
  (5) `z^n = u^(4n/5) w^n` is a power series in `u` iff `5 | n`, so the power-series solutions are `Phi^(0) c(z)`, `c` in `1 + z^5 Q[[z^5]]`
  [W: POLY s.2]. (6) Resonances: they hold for every choice of the `f_m` once one power-series solution exists [T: POLY s.2, one
  sentence; expand with (4): every other choice is that solution times some `c(z)`]; P(a) supplies one without [U] [W: PR s.2].
* **`prop:pol-deg`** (ER's elementary proof; recommended, it keeps the note free of a quotation). (1) In `x = 1/s = L` the conjugated
  operator is `M = zeta^i x - (z/5)(1 - x^5) th_x`, and `(5M + kz) x^m = 5 zeta^i x^(m+1) + z(k - m) x^m + z m x^(m+5)` [W: ER s.4]. (2) Lemma:
  `Pi := (1/3125) prod_(k=1..5)(5M + kz)` maps `x^m` (`m >= 1`) into `x^6 Q[x][z]` (the factors commute; apply `k = m, .., 5` first, each
  kills the lowest power when it equals `k`; no factor lowers the valuation) [W: ER s.4]. (3) `P = (x^-5 - 1) Pi + M^5 - 1` maps
  `x Q[x][z]` into itself, componentwise in `z` [W: ER s.4]. (4) `T_1(x^m) = zeta^-i (1 - m)(x^(m-1) - x^(m+4))`; by induction on `n`, `psi_n`
  lies in `x Q[x]`: a lowest power `x^(m_0)`, `m_0 <= 0`, would give the uncancellable term `(1 - m_0) c x^(m_0 - 1)` [W: ER s.4]. The
  induction presupposes that `psi_n` is a Laurent polynomial in `x` (true by `lem:pol-poly`); ER does not say so [T]. (5) Triviality of the
  resonance: "the right-hand side has top `s`-exponent `<= -1`, so W stops at `s^-6` and never reaches `s^-5`" [T: ER s.4, one sentence; the
  writer must expand it from ER s.4, else `\draftnote`]. (6) Translation: `psi_n = A_n` in `L Q[L]` of degree `<= 4n + 1` iff
  `deg ph_n <= floor(4n/5)`; the value at `L = 0` [W: ER s.4].
  *Second proof* (for a remark): `F(w, x)` of arXiv:0710.0889 is the colour-0 Frobenius series (`w = 1/z`, `x = q`), its expansion (5)
  is coefficientwise in `x` (their Lemma 1), and Thm 4(ii) gives `Phi_s` in `L Q[L]`; the recursion WKB is `d(z) h(z) Sum_s Phi_s z^s/L` for
  some `d`, so `A_n` lies in `L Q[L]`; no negative powers of `u` because every WKB column is `Phi^can d(z)` [W: GAM s.1; GR s.2.1, C1; PR C1].
* **`lem:pol-Finf`.** Evaluate `Phi^QRR = Phi^(0) F(z)` (`lem:pol-poly`) at `u = 1`, where `z = w` [W: POLY s.2 Lemma 4].
* **`lem:pol-p5k`.** (1) `y_j := 1/(1 - zeta^j)` satisfy `y_(-j) = 1 - y_j = -zeta^j y_j`, so `y_(-j)^5 = -y_j^5` and the pairs `j, -j` cancel in
  `Sum_j y_j^(5k)` for odd `k` [W: EVID s.8.2; ER s.2]. (2) Insert in `log h = -Sum_k b_k z^(2k-1) Sum_(a in w_0) a^(1-2k)` at `2k - 1 = 5k'`;
  only `a = -5` survives [W: EVID s.8.2; sign re-derived ER s.2].
* **`prop:pol-Ginf`.** (1) [LMI0] `c(z; 0) = Ginf(z)` [W: RA s.3.4; LF s.4.2, other digest]. (2) No pole at `L = 0`; with Theorem Gamma
  (i), (iv): `c_n` in `Q[L]` [W: LF s.4.3]. (3) [HF]: `f^SP_m = c_(5m)(0) = [z^(5m)] Ginf` [W: RA s.3.4; LF s.4.3]. (4) `thm:thp-b`: `f^can = f^SP`
  [W: SYN5b K7]. (5) `f^QRR = f^can` from `D = 1` [W: SYN5b K7].

## 3.2 The finite route

* **`lem:fin-P`.** (1) `M_0(t)^-1 = M_0(-t)^T`, `Psi0^-1 = pinv` with `pinv^T pinv = (1/5) eta'`, so `pinv = Psi0^T eta'/5` and
  `R^(t)^-1 = (1/5) R^(-t)^T eta'` [W: EVID s.2.2; ER s.2]. (2) Prop. S form of `R^`, covariance of `Qhat`, `Sum_i zeta^(im) = 5 [5 | m]`;
  the re-derived sum `(1/5) Sum_i zeta^(i(a+5-c-N))` gives the selection rule [W: EVID s.2.2; ER s.2].
* **`lem:fin-N`.** (1) By [U], `P_can Q^m` is a formal WKB solution of the quintic with the same exponentials [W: EVID s.2.3]. (2) Hatted it
  is `Ghat^can R^(tau w)`, a power series in `u` at fixed `w` [W]. (3) Its colour-0 row 0 is `C0^-1 Phi'`, and HPF forces `ph'_n` to solve the
  same recursion (rows by the cyclic system: [CYC], PS-C9) [W: EVID s.2.3]. (4) Two WKB solutions with the same exponential differ by a
  constant scalar series `d(z)` (`lem:glu-wkb`, in formal series in `z` over `Q((s))`) [W: EVID s.2.3; ER s.2]. (5) Both sides are power series
  in `u` at fixed `w` and `z^j = w^j u^(4j/5)`, so `d` lies in `Q[[z^5]]`, the same for all colours (`d(zeta^-i z) = d(z)`) [W: EVID; ER s.2]. (6)
  (R1) gives `P_can,11(0,z) = P_can,22(0,z) = 1`; the hatted `(1,1)` entry at `w^(5q)` is `tau^(1/2) u^(4q) p_(5q,1)(u)`; hence the claims; slope
  `F_2(0)/(C1(0) C0(0)) = 1` [W: EVID s.2.3; ER s.2].
* **`thm:fin-conn`.** Part 1: the general solution is `Phi^pol d(z)` with value `Phi^pol(w, 1) d(w)` at `u = 1`; [QRR] demands `h`; if
  `ph_n(1) = h_n` for all `5 ∤ n` then `d = 1`, otherwise `d` has a `z^n` term with `5 ∤ n` [W: EVID s.2.4]. Part 2: `Q = P_can c Q^m` (Lemma N),
  so `D = c^-1`, and `D = 1` iff `K_q = 0` for all `q`; slope 1 fixes `f^can_q` [W: EVID s.2.4]. (iii) If `ph^can_n(1) = h_n` for all `n`, then
  `Phi^can` satisfies the QRR condition, so it is the QRR solution, (H) holds and `D = 1`; conversely [W: ER s.2].
* **`prop:fin-pfree`.** ER s.3 steps 1-6: `Qh_1 = (1/(C0 C1)) A_(-1) Phi + (1/C1) w (1-u)(th C0/C0^2) Phi`; the second part pairs, by Lemma P
  with `(a, c) = (0, 1)` and Lemma N, to a multiple of `p_(5q-1,1)` which vanishes at `u = 0` by (R1), and `th C0 = O(u)`; then REG gives
  `S_q/(C0 C1) = O(u^(4q))` with `C0(0) C1(0) = -1` [W: ER s.3].
* **`prop:fin-row`.** (1) Row 0 of `z P' = A P - P At`: `z P_01' = b_0 P_11 - bt_0 P_00`, `b_0(0) = bt_0(0) = 1`, `P_01(0, z) = 0` by (R1) [W: POLY
  s.3; PR s.2]. (2) Read `[z^(5q)]` through the first-row expansion D1.16 with the `mu` table [W: POLY s.3]. (3) `m5`: `p_(4,1)` does not depend on
  the normalisation; `p'_(4,1)(0) = m5` [W: POLY s.3; numbers V].

## 3.3 Theorem P

* **`lem:thp-form`.** [W: POLY s.5 (phase, amplitude, cyclotomic identities); PB s.2.1 (re-derived); PR s.2 (by hand); GR s.2.2 (from the
  six Stirling series)]. (d): `W_GAM = -W_POLY` and the amplitude identity `Sum_j (zeta^j x - 1)^(1-2k) = -p_(2k-1)(x)` [W: Stage-4 K6]; the two
  summands: `d^2/d sigma^2` of `(c - sigma) log(c - sigma)` and of `-(sigma - c) log(sigma - c)` agree, `b_k((c-sigma)/z)^(1-2k) = -b_k((sigma-c)/z)^(1-2k)`,
  branches only move the critical point by `2 pi i k` in `S'` [W: PB P9]; the reflection version [W: CP s.1.7].
* **`thm:thp-a`.** (1) `PF5 = z^5 D PF4 - 1` [W: PB s.2.2]. (2) `PF5[q^s phi] = g(s) - g(s+1)`, `g = q^s((zs)^5 - 1) phi`, from
  `((z(s+1))^5 - 1) phi(s+1) = z^5 prod_k (5s + k) phi(s)` [W: POLY s.6 (a); checked by hand PR s.2, PB s.2.2]. (3) The Stirling form obeys the
  same recursion by L3(iv) (`Gamma_St(y+1) = y Gamma_St(y)`), the shift realised as `sigma -> sigma + z` [W: PB s.2.2; FSL4 s.9.2]. (4) FL4
  (`P(zD_q) SP = SP[P(sigma) ..]`) and the shift lemma: `PF5[SP] = SP[g(sigma) - g(sigma + z)] = 0`; `g(sigma + z) - g(sigma) =
  d/dsigma Sum_(k>=1) z^k g^(k-1)/k!` is at each `z`-order a convergent sum of exact forms [W: PB s.2.0 FL7, s.2.2; FSL2 s.3.1 Remark,
  s.3.2 (Lemma T(b)), s.4.5; RL s.3.4]. (5) Flat SP = hatted SP: L1a for `sigma = L/X` and Lemma H for `z = L^-4 w` [W: RL s.3.4 table; CP3
  s.3.10 with RF2 C2]; PB states it "on `0 < u < 1`" [T: PB s.2.2; the comparison ring should be named, s.7]. (6) The hatted Taylor data lie in
  `Q[[u]]` (`[t^k] What(1+t) = Sum_m u^(m-1)[(-1)^k + C(5m-1, k)/(5m-1)]/m`, `1/What''(1) = (1-u)/5`, `B` polynomial in `X`, `w` per `u`-order),
  so `[w^n] Phi^SP` lies in `Q[[u]]` and solves HPF [W: POLY s.5-s.6 (a); PB s.2.2; GP s.3]. (7) `lem:pol-poly` [W: PB s.2.2].
* **`lem:thp-cfu`.** Write `Y = X y(T)`; then `V/Y + Y^4/4 = X^-1 (V/y + T y^4/4)`; at order `u^k`: `(T - 1) Y_k + V_k = R_k(T)`, `R_k` a polynomial
  from lower orders; so `V_k = R_k(1)` and `Y_k` is a polynomial; uniqueness; `X = 1` goes to a critical point since `Y'(1)` is in
  `1 + uQ[[u]]`; `What(1) = (5/4) V^(4/5)` gives `V = tau^(-5/4)`; inverse of the same form [W: PB s.2.4; POLY s.5]; `deg_T Y_k <= k` from
  `deg R_k <= k + 1` [W: PS s.2.3]. GP s.4 proves the normalised variant by the same induction [W: GP s.4; GR s.2.4].
* **`lem:thp-red`.** (1) The exact-form identity and the relation `Y^n == (n-3) w Y^(n-4) + V Y^(n-5)` [W: PB s.2.5]. (2) Termination: each step
  lowers the degree by at least 4, the lowest exponent produced is `-1` [W: PB s.2.5]. (3) `Gamma` is one amplitude in `A[[t]][[w]]`, so FL1
  over `A = Q[[u]]` applies once [W: RL s.3.4]. (4) Z/5 grading [W: PB s.2.5; PR s.3].
* **`lem:thp-S2`.** PS s.1 (i)-(vi) [W]: (i) data over `A = Q[[v - 1]]`, `J_m := SP_(x = sigma)[x^m e^(lambda/z) dx]`, Hessian `5 sigma^2` a unit;
  (ii) FL4 with a moving critical point gives `z dJ_m/dv = -J_(m-1)`, FL1 for `d(x^0 e^(lambda/z))` gives `J_3 = v J_(-2)`, so
  `(J_0, -J_(-1), -J_3, J_2, J_1)` solves the model system (row `b = 1` from `J_(-2) = J_3/v`); (iii) L1a with `x = sigma x'` and Lemma H with
  `mu = sigma^4`: `J_(c_b) = sqrt(2 pi z/5) e^(U/z) sigma^(c_b - 1) G_b(z/sigma^4)`; (iv) the chain rule turns the system into Prop. S's T-system;
  (v) uniqueness: `F_0` enters the closing relation at `t^(n+1)` with coefficient `4n` (`Sum k_i = 0`); (vi) `G_0(0) = 1`, so `eps_b G_b = F_b`.
* **`lem:thp-id`.** (a) `Y = Y_c eta` (L1a, linear over `Q[[u]]`), `f(Y_c eta) = Y_c^4 (1/eta + eta^4/4)`, `w/Y_c^4 = tau w` (Lemma H with
  `mu = Y_c^4`), `Y_c^(c_b + 1) tau^(1/2) = tau^((1-c_b)/4)` [W: PB s.2.6; PS s.1 "Weights"]. (b) SP at `X = 1` equals SP at `Y_c` (FSL s.5.4); FL1
  kills `d(Gamma e^(-f/w))`; `A[[w]]`-linearity [W: PS s.1]. (c) `Q^m_b = sigma^(c_b - 1) F_b(z/sigma^4)`, `sigma = s tau^(-1/4)`, `Qhat_a = s^(1-c_a) Q_a`,
  so `Ghat_ab = s^(c_b - c_a) tau^((1-c_b)/4) P_ab` [W: PB s.2.6; PS s.1 "Gauge side"]. (d) Least `w`-degree argument (`F_b(0) = eps_b != 0`) [W: PB
  s.2.6; PS s.1]. No circularity: s.2.8 of PB uses neither s.2.6 nor s.2.7 [W: PS s.1].
* **`lem:thp-relab`.** PS s.2.2 [W]: (1) labels and `iota`, injective; (2) the family: `(k, r, multiset of atoms of l = log(e^(-G/z) at))`, atoms
  `-x^(5m-1)/(m(5m-1) z)`, `x^(5n)/(2n)`, `b_k c_(k,P) (xz)^(2k-1) x^(5P)`; (3) original fibres finite (`k <= E`; atoms of positive `x`-degree);
  (4) hatted fibres finite (atoms of hatted `u`-order `>= 0`, only `-X^4/(4w)` of order 0; the `X`- and `w`-balances fix `k` and its
  multiplicity, determinant 5); (5) `H` computed in `L[[u]]` is the hatted-order sum: `e^(-v(Lam/u)/X) = e^(-v/X) e^(-v(Lam/u - 1)/X)` holds
  coefficient by coefficient (also in FSL's ring `B`); (6) one family summed in two orders with finite fibres. PB s.2.7 has the coarse
  version (index `(k, r, b)`), legitimate by PS s.2.2's last remark.
* **`lem:thp-key`.** The exponent count `beta = j - 5i' - m - 1`, `E = i + l - (4beta - 1)/5` [W: PB s.2.7; PS s.2.3]. Attainment [V: PS T3].
* **`lem:thp-weight`.** PS s.2.3 (a)-(b) [W]: monomials of `B` have weight 0; CFU and its inverse are `X(1 + Sum c u^k T^d)`, `d <= k`; the
  two reduction steps keep `omega` or raise it by `5(l + 1)`.
* **`prop:thp-res`.** (1) FL2-FL3 for the residue in `L[[u]]` (PB's ring with PS-C2's condition, or FSL's `L`) [W: FSL4 s.10.1-s.10.3; PB s.2.7].
  (2) Model residues: `x = sY` gives `f/w = (x^4/4 + v/x)/z` with `v = uV` [W: PB s.2.7]; termwise identity with the `i`-th term of
  `s^(-(m+1)) I^C_m` [W: PS s.2.1]. (3) `lem:thp-relab`, `lem:thp-key`. (4) The hatted product `pi_m . Res_Y[..]` has finite fibres and `iota` is a
  monoid map [W: PS s.2.3, last bullet]. (5) `s^(-c_b) G_b = P_0b`, `Y^C_b = eps_b I^C_(c_b)`; signs cancel (`eps_b^2 = 1`), powers of `s` cancel
  (`s . s^(-(c_b+1)) . s^(c_b) = 1`), no `tau` survives [W: PS s.1 "`kappa = C0`"; PR s.3].
* **`thm:thp-b`.** Step 1 CFU [W: PB s.2.4]. Step 2 reduction [W: PB s.2.5]. Step 3 invariance of both functionals and vanishing on exact forms
  [W: PB s.2.0, s.2.7; FSL s.5.4, s.10]. Step 4 identification [W: PB s.2.6; PS s.1]. Step 5 [W: PB s.2.8; PS s.3]: both `Phi^SP` (P(a)) and `Phi^can`
  ([CYC]) solve HPF with `ph_0 = 1`; `lem:pol-poly` gives `Phi^SP = Phi^can s(z)`, `s` in `1 + z^5 Q[[z^5]]`; rows and colours commute with `s(z)`
  (`s(zeta z) = s(z)`), so `P^SP = s(z) P_can`, regular and graded; Theorem G for `P_can` times `s(z)`: `(P^SP Y^C)_0 = s(z)(Y^C_Q)_0` (PS-C6: PB's
  four-item re-check is redundant). Step 6 `prop:thp-res`: `V = C0 s(z) (Y^C_Q)_0` [W: PB s.2.7]. Step 7 [W: PB s.2.9; PS s.3]: `V_1 = -1/z` (only
  `k = 1`; `at(0, z) = 1`; `G = O(x^4)`); `[u^1](Y^C_Q)_0 = -1/z` (Lemma Q(a): only `n = 1`, `y^(1)_0 = -t_c`); `[u^0](Y^C_Q)_0 = 0` (PS-C7);
  `C0(0) = 1`; so `s = 1`. Remark: the `z^(5q-1)` coefficient at `u^1` is `-s_q` and uses the reduction up to `u`-order `4q` [W: PB s.2.7].
* **`cor:thp-A2`.** The QRR solution is the unique WKB column `Phi^can d(z)` with value `h` at `u = 1`, and `Phi^SP` is a polynomial solution
  [W: POLY s.6 Corollary; PR s.2]; then [HF] [W: CP3 s.3.10].
* **`rmk:thp-Ai`.** GR s.2.5 [W]: PF for the residue (GP s.2), rows in `L[[u]]^5` (needs [CYC], PB P1), GI Lemma 1, `Y(0)` depends on
  finitely many coefficients, a weight function `w(f)` shows `Y(0)` has degree `<= 2` in `(alpha, beta)`; the value on a 3x3 grid [V: GR,
  `rf_run_res` R3]. The last step is a finite computation: see s.7.

## 3.4 The glue

* **`lem:glu-wkb`.** [T: PS-C10 and PS s.4.2 (b), one line]: the exponent is a simple root, the transport operator is first order with kernel
  `Q L` (ER s.4's `T_1`); the induction on the first differing coefficient (subtract `c_n z^n A`) is implicit and is expanded from the same line.
* **`thm:glu-main`.** PS s.4.2 [W]: (a) `A^can_n = L^(4n+1) ph^can_n(L^-5)` lies in `Q[L, 1/L]` (`lem:pol-poly`); (b) in `Q[[q]]` both `A^can` and
  `A^QRR` are WKB columns of PF5 with exponential `e^(u_0/z)`, so `A^can = A^QRR d(z)`, `d` in `1 + zQ[[z]]` (`lem:glu-wkb`); (c) hence `A^QRR_n` lies
  in `Q[L, 1/L]`; (d) `Q[L, 1/L][[z]]` embeds in `Q((s))[[w]]`, where [R0] holds; `Phi^QRR` has constant term 1, not a zero divisor; row 0,
  colour 0, hatted, of [R0] reads `Phi^QRR(w, u) D_0(s^4 w) = Phi^can(w, u)` (the factor `1/C0` is common) [W: PB s.3], so `D_0 = d`; (e) at
  `L = 1` (`s = 1`, `z = w`), `Phi^QRR(z, 1) = Sum_n A^QRR_n(0) z^n = h(z)` by H9; so `h D_0 = Phi^can(z, 1)`. Then `thm:thp-b` and [HF] [W: PB s.3].
  PS-C8's sentence (the ring of comparison; `D` after stripping the same exponential) must be in the text.
* **`cor:glu-Laurent`.** Step (c) above [W: PS s.4.2 (c)].
* **`prop:glu-same`.** Lemma N with `K_q = [z^(5q)] c` and `P_11(0, z) = c(z)` for `P_can c`: `K_q = 0` for all `q` iff `c = 1` iff `Phi` is
  `C0 [hatted P_can Q^m]_00` [W: PS s.4.4].
* **`cor:glu-colours`.** [W: PS s.4.5; UR s.4 (graded `P`, covariance under `z -> zeta z` with `diag(zeta^-a)`)].
* **`thm:glu-D`.** Assembly [W: SYN5 s.2 S1-S7; RS2 s.6.2]: [GAM] gives `c(z; 1) = h(z)`; `thm:glu-main` gives `D_0 = c(z; 1)/h(z) = 1`;
  `cor:glu-colours` gives `D = 1`; (F_g) for every `g`. (CP3 s.3.11 uses Lemma N instead of the glue; both valid, K3.)
