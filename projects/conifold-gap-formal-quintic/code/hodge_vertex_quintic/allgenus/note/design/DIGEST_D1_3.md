# Digest D1, part 3: s.2 statements (continued): Theorem P and the glue

Continuation of `DIGEST_D1_2.md`. Same format and abbreviations.

## 2.4 Theorem P: the canonical normalisation is the full-Stirling normalisation

* **`lem:thp-form`** (lem) *one 1-form.* (a) With Stirling for `Gamma(1 + 5s)` and `1/Gamma(1 + s - zeta^j/z)`, `sigma = zs`, `x = 1/sigma`: the
  exponent is `S(sigma)/z`, `S = sigma log q + 5 sigma log 5 sigma - Sum_j (sigma - zeta^j) log(sigma - zeta^j)` (the `log z` terms and the
  linear terms cancel, `Sum_j zeta^j = 0`); `S(1/x) = -W`, `W = Lam/x + G(x)`, by the cyclotomic identity for `G` (only `5 | n` and
  `5 | n+1` survive the sum over `j`); `dW/dx = -(Lam + log(1 - x^5))/x^2`, critical points `x^5 = u`. (b) `sqrt(5 sigma/(sigma^5 - 1)) d sigma
  = -sqrt5 (1 - x^5)^(-1/2) dx` and `log at` as in D1.13 (from `Sum_j (1 - zeta^j x)^(-n) = 5 Sum_(5|i) C(n-1+i, i) x^i`); `at` does not
  depend on `u`, `at(0, z) = 1`. (c) The hatted data of D1.13. (d) GAM's flat form (`e^(W/z)`) and POLY's hatted form (`e^(-W/w)`) are the
  same 1-form (`W_GAM = -W_POLY`, equal amplitudes); GAM's summand `(-q)^s Gamma(1+5s) prod_j Gamma(zeta^j/z - s)` and POLY's
  `q^s Gamma(1+5s)/prod_j Gamma(1 + s - zeta^j/z)` give the same normalised series. Src. POLY s.5; PB s.2.1; PR s.2; GAM s.2; GR s.2.2;
  CP s.1.7; Stage-4 K6; PB P9. St. PROVED. Corr. Stage-4 O7 / GR C6: this phase is NOT UR s.8's `W`; they share only the critical data.
* **`thm:thp-a`** (thm; Theorem P(a)) `Phi^SP` solves HPF, and every `[w^n] Phi^SP = ph^SP_n(u)` is a polynomial in `u`. Hyp. [FSP]: L1c(iii)
  and FL4 (SP commutes with `D_q`; `P(zD_q) SP[..] = SP[P(sigma) ..]` when `D_q F = sigma`), L3(iv), Lemma T(b) (or FL1 with a convergent
  sum, FSL2 s.3.1 Remark), L1a with Lemma H (flat to hatted); `lem:thp-form`; `lem:pol-poly`. No [U] (PR s.2). Src. POLY s.6 (a);
  PB s.2.2; PR s.2; RL s.3.4. St. PROVED (confirmed by PR; RL table "PROVED"). Not re-checked by PS (RS s.3 row P(b)). Corr. PB P5 (the
  total derivative is, at each `z`-order, an infinite convergent sum: PB's FL7); O12 (FL3, FL6 are PROVED in FSL, not "only sketched").
* **`lem:thp-cfu`** (lem; Lemma CFU) There are unique `Y = X(1 + Sum_(k>=1) u^k Y_k(X^5))`, `Y_k` in `Q[T]`, and `V(u)` with
  `What(X; u) = V/Y + Y^4/4`. Then `V = tau^(-5/4)` (so `v = u tau^(-5/4)` is recovered from the integral), `Y_1 = (1 + T)/18`, `deg_T Y_k <= k`;
  `X = 1` goes to the critical point `Y_c = V^(1/5)` of `f = V/Y + Y^4/4`, `Y'(1)` in `1 + uQ[[u]]`; the inverse `X(Y) = Y x(Y^5)` is
  polynomial at each `u`-order; the map commutes with `X -> zeta X`. Src. PB s.2.4; POLY s.5; PS s.2.3 (degree); GP s.4 (normalised
  variant, s.4 below). St. PROVED; VERIFIED (PR exact `k < 10`; GR `k <= 16` for GP's variant; PB `u^20`).
* **`lem:thp-red`** (lem; reduction modulo exact forms) Put `beta(Y) dY := B(X(Y)) X'(Y) dY`, `beta` in `Q[Y, w][[u]]`. From
  `w d(Y^m e^(-f/w)) = (m w Y^(m-1) + V Y^(m-2) - Y^(m+3)) e^(-f/w) dY`: `Y^n == (n-3) w Y^(n-4) + V Y^(n-5)` (`n >= 4`); the reduction
  terminates, `beta e^(-f/w) dY = Sum_(m=-1..3) pi_m Y^m e^(-f/w) dY + d(Gamma e^(-f/w))`, `Gamma` in `Y Q[Y, w][[u]]`, `pi_m` in `Q[w][[u]]`.
  Z/5 grading (`deg Y = 1`, `deg w = -1`; `f(zeta Y) = f(Y)/zeta`, `B(zeta X; w) = B(X; zeta w)`): `pi_m` contains only `w^k`, `k = m mod 5`;
  for `m = c_b`, `k = -b mod 5`. Src. PB s.2.5; POLY s.6 step 2; PR s.3 (grading); RL s.3.4 (`Gamma` is one amplitude, no `u`-adic sum of
  exact forms). St. PROVED. Uniqueness of the `pi_m` is not needed (`lem:thp-id` pins them).
* **`lem:thp-S2`** (lem; (S2)) For each `b`: `SP_(eta=1)[eta^(c_b) e^(-(1/eta + eta^4/4)/t) d eta] / (sqrt(2 pi t/5) e^(-5/(4t))) = eps_b F_b(t)`.
  Hyp. [FSP] FL4 with a moving critical point (PS-C1), FL1 once (at `m = 0`), L1a with Lemma H (PS-C3); [PropS]. Src. PB s.2.6; PS s.1.
  St. PROVED, refereed once (PS); VERIFIED. Corr. K12: PR s.3 had it VERIFIED only, so Stage 4's "R5 PROVED" rested on one VERIFIED step.
* **`lem:thp-id`** (lem; signs, weights, identification) (a) `SP_(Y_c)[Y^(c_b) e^(-f/w) dY] = N_w eps_b tau^((1-c_b)/4) F_b(tau w)`,
  `N_w := sqrt(2 pi w/5) e^(-f(Y_c)/w)`, `f(Y_c) = What(1)`. (b) `Phi^SP = Sum_b eps_b pi_(c_b) tau^((1-c_b)/4) F_b(tau w)`. (c) For every
  regular gauge, row 0 reads `Phi/C0 = Sum_b tau^((1-c_b)/4) G_b F_b(tau w)`, `G_b = s^(c_b) P_0b` graded (only `w^k`, `k = -b mod 5`; `[w^k] G_b`
  divisible by `u^((c_b + 4k)/5)`). (d) Triangular uniqueness: a graded family `Delta_b` with `Sum_b tau^((1-c_b)/4) Delta_b F_b(tau w) = 0` is
  zero. Hence `G^SP_b = eps_b pi_(c_b)/C0`. Hyp. `lem:thp-S2`, `lem:thp-cfu`, `lem:thp-red`; [FSP] s.5.4, FL1 over `Q[[u]]`, `A[[w]]`-linearity;
  `P^SP` graded (step 5 of `thm:thp-b`; no circularity, PS s.1). Src. PB s.2.6; PS s.1; PR s.3. St. PROVED, refereed once; VERIFIED.
* **`lem:thp-relab`** (lem; relabelling) `V = s Res_X[e^(-What/w) B dX]`, coefficientwise: `[u^N w^j] H = [u^(N + (1-4j)/5) z^j] V` with
  `H := Res_X[e^(-What/w) B dX]` computed in `L[[u]]`. Both sides are sums of ONE family (D1.13 atoms) over two labellings related by the
  injective map `iota(E, b, al) = (E + (4b + al)/5, b, al)`, with finite fibres in both orders. Src. PB s.2.7; PS s.2.2 (complete argument,
  PS-C4). St. PROVED, refereed once. History: FSL4 s.10.4 SKETCH; RL s.3.4 "open"; K5 (critic: complete).
* **`lem:thp-key`** (lem; key lemma) In `s Sum_m pi_m Res_Y[Y^m e^(-f/w) dY] = Sum_b pi_(c_b) s^(-c_b) I^C_(c_b)(v(u), z)`, the key `(u^E, z^beta)`
  receives `pi` only from `u`-orders `i <= E + (4 beta - 1)/5`, through finitely many terms; only `beta = 4 mod 5` occurs. The bound is
  attained. Src. PB s.2.7; PS s.2.3. St. PROVED, refereed once; VERIFIED. Corr. PS-C5 / A6: finite fibres give a formal sum only;
  membership in `Q((z))[[u]]` needs `lem:thp-weight`, or follows from the coefficientwise identity with `V` in `Q[z^-1][[u]]` (the chain).
* **`lem:thp-weight`** (lem; optional) With `omega(u^N Y^n w^j) := 5N - n - 4j`: `omega >= 0` on `beta`, `X'(Y)` and every reduction step,
  so `5i - m - 4j >= 0` for every `u^i w^j` in `pi_m`, and `E >= 0`; at `i = 4q + 3` the top coefficient is `-C0(0) p_(5q+4,1)(0) = 0`, so the
  mechanism reproduces (R1). Src. PS s.2.3; RS A6. St. PROVED (PS; checked by the critic); not used by the chain.
* **`prop:thp-res`** (prop; `kappa = C0`) `V = Sum_b eps_b pi_(c_b) s^(-c_b) Y^C_b = C0 Sum_b P^SP_0b Y^C_b = C0 (P^SP Y^C)_0` in
  `Q((z))[[u]]`, in every frame. Hyp. `lem:thp-relab`, `lem:thp-key`, [FSP] FL2-FL3 (residue) in the ring `L` (PS-C2), the model residues
  `Res_(Y=0)[Y^m e^(-f/w) dY] = s^(-(m+1)) I^C_m(v(u), z)` (put `x = sY`, `v = uV`), [GI] Lemma V, `lem:thp-id`. Src. PB s.2.7; PR s.3; PS s.1
  ("`kappa = C0`: confirmed"), s.2. St. PROVED; refereed (PR, PS); VERIFIED.
* **`thm:thp-b`** (thm; Theorem P(b)) `Phi^SP = Phi^can`. **Complete hypothesis list** (union of SYN5 s.1 table row P(b), SYN5b K11 and
  O10, RS s.3 row P(b), PS2 s.7; PB s.1 H1-H7):
  * [U] for the pair `(A, At)`, whose hypotheses hold (`b_0 = b_2 = t_c'`, `b_3 = b_4 = -varpi_0/(1-u)` units, `u b_1` a unit,
    `bt = (v', -v'/v, v', -v', -v')`, `Prod b = Prod bt` iff `v^(4/5) = u^(4/5)/tau`) (PB s.2.8; PS s.3);
  * [GI] Theorem G with Lemma 1 and (R1) (applied to `P_can`, then multiplied by `s(z)`, PS-C6); [GI] Lemma V, Lemma Q(a), `y^(0) = -e_1`;
  * [CYC] (step 5: `Phi^can` solves HPF) (PB P1; K4);
  * [PropS] uniqueness half (step 4, (S2));
  * `lem:pol-poly`; `thm:thp-a`; `lem:thp-form`, `lem:thp-cfu`, `lem:thp-red`, `lem:thp-S2`, `lem:thp-id`, `lem:thp-relab`,
    `lem:thp-key`, `prop:thp-res`;
  * [FSP]: FL1, L1a and s.5.4, L1c(iii)/FL4, FL2-FL3 (residue) with the ring `L`, L3(iv), Lemma T(b), Lemma H;
  * the frame lies in the family `varpi_0 = f_0 + alpha g + beta f_2`, `t_c = g/varpi_0` (`C0(0) = 1`, `t_c'(0) = 1`) (RS2 items 27-28; A1);
  * `ph^can_0 = 1`: from `p_(0,0) = tau^(-1/4)/varpi_0` (GI s.5, not read by me; PB s.2.8; PS s.3).
  * Not used: [QRR], Lemma N, Lemma P, the degree bound (only polynomiality), (Gamma), the chart argument (PB s.1 H9 as rephrased by O13).
  Src. POLY s.6 (b); PR s.3; PB s.2; PS s.1-s.3. St. PROVED: Stage 4 (POLY + PR, refereed), re-derived in PB, PB s.2.6-s.2.9 refereed once
  (PS); VERIFIED. Corr. O10 (Stage 4's "(U), GI Thm G and four formal lemmas" is incomplete), PR C6, PB P1, PS-C1..C7, O8 (CP s.0.3's list).
* **`cor:thp-A2`** (cor) `D = 1` iff `Phi^QRR = Phi^SP` iff `Phi^SP(z, 1) = h(z)` (POLY's "(A2)") iff, by [HF], (Gamma) `c(z; 1) = h(z)`. This
  is a connection statement (a series built near `u = 0` compared at `u = 1`) from which the gauge, the periods and the model have dropped
  out. Src. POLY s.6 Corollary; PR s.2; GAM Cor. A; SYN4 s.1.2 (R7), s.1.3. St. PROVED given `thm:thp-b`. Corr. Stage-4 K5, O5, O12 (PR C4:
  not "a statement at `q = 0` alone"); Stage-4 K4 (naming: POLY's "(A2)" is this identity, GAM's "(A2)" an analytic step; write (Gamma)).
* **`cor:thp-frame`** (cor) `Phi^can` does not depend on the conifold frame. Src. PB s.0 (it follows from P(b), `Phi^SP` being frame-free).
  St. stated in PB; VERIFIED in two frames (PB2 s.4).
* **`rmk:thp-Ai`** (rmk; GAM's Theorem A(i)) `Res_(x=0)` of the rows of the full-Stirling 1-form is `Y^C_Q` in every conifold frame, without
  [U]; hence `varpi_0 t_c = -Sum_(k>=1) alpha_(k-1) (log(1-u))^k/k!`, `alpha_m = [y^m] exp(Sum_k b_k (5^(1-2k) - 5) y^(2k-1))`. An alternative
  ending of step 7 (PB s.2.9). Src. GAM s.3; GP s.5; GR s.2.5 (R-a). St. PROVED (GR), but its proof decides `Y(0)` by a finite exact
  computation on a 3x3 grid in `(alpha, beta)` (see s.7); its rows need [CYC] (PB P1). Not needed by the chain.
* **`rmk:thp-710`** (rmk) At `u^2`: `-1/(2z)` from `Lam` and `-1/(5z)` from the Stirling term `b_1 = 1/12` add up to `-7/(10z) = -gamma/z`, the
  `u^2` coefficient of the vanishing period. Src. POLY s.5; GP s.5 Remark; MS4 s.3. St. part of `prop:thp-res` (with `s = 1`); VERIFIED.

## 2.5 The glue and the proof of `D_0 = 1`

* **`lem:glu-wkb`** (lem; formal WKB uniqueness) Two colour-0 WKB columns of PF5 with the same exponential `e^(u_0/z)` and the same leading
  term `L` differ by a constant series: `A' = A d(z)`, `d` in `1 + zQ[[z]]`; this holds in `Q[[q]]` (`L = (1 - 3125q)^(-1/5)`) and, for
  Lemma N, in formal series in `z` with coefficients in `Q((s))`. Proof in one line: the exponent is a simple root of the symbol (the roots
  `zeta^i L` are distinct), so the transport operator is first order with kernel `Q L`. Src. PS s.4.2 (b) and PS-C10; ER s.2 (Lemma N step 4);
  ER s.4 (`T_1(x^m) = zeta^-i (1-m)(x^(m-1) - x^(m+4))`); PB s.1 H8. St. PROVED (PS-C10; PB had "standard").
* **`thm:glu-main`** (thm; the glue) `h(z) D_0(z) = Phi^can(z, 1)` in `Q[[z]]`. With `thm:thp-b` and [HF]: **`D_0(z) = Phi^SP(z, 1)/h(z) =
  c(z; 1)/h(z)` exactly.** Hyp. [R0] (hence [U], [CYC]); polynomiality of `ph^can_n` (`lem:pol-poly` with [U], [CYC]; or P(a) with P(b));
  `lem:glu-wkb`; [QRR] as H9. Not used: Lemma N's formula for `K_q`, Lemma P, the degree bound, (Gamma), the chart argument. Src. PB s.3;
  PS s.4.1-s.4.5; RS A4; SYN5 S5; MS5 s.2 item 5. St. PROVED, refereed once (PS); the refereed alternative is EVID s.2.3-s.2.4 with ER.
  Corr. PS-C8 (one sentence: the ring of comparison `Q[L, 1/L][[z]]` inside `Q((s))[[w]]`; `D` after stripping the same exponential); K3, P7
  (Lemma N not needed); O13.
* **`cor:glu-Laurent`** (cor) Every `A^QRR_n` lies in `Q[L, 1/L]`: the QRR amplitudes are Laurent polynomials in `L` (no degree bound used),
  so the actual rows exist as `s`-series near `u = 0`; this is not an extra input. Src. PS s.4.2 (c); RS A4; RS2 s.6.3. St. PROVED (by-product;
  checked by the RS critic).
* **`prop:glu-same`** (prop) EVID's `Phi^can` (`K_q = 0` for all `q`) and POLY's (`P_11(0, z) = 1`) are the same element of `Q[u][[w]]`.
  Src. PB s.3; PS s.4.4. St. PROVED given Lemma N (hence [U], [CYC]); VERIFIED (PB2 s.4).
* **`cor:glu-colours`** (cor) `D_i(z) = D_0(zeta^-i z)` and `h_(w_i)(z) = h_(w_0)(zeta^-i z)` (`Sum_(a in w_i) a^(1-2k) = zeta^(i(1-2k))
  Sum_(a in w_0) a^(1-2k)`); hence `D_0 = 1` gives `D = 1`. Src. UR s.4; PB s.3; PS s.4.5; CP s.0.2, CP3 s.3.9. St. PROVED.
* **`thm:glu-D`** (thm; CP's "Theorem D") `D_0 = 1` in `Q[[z]]` to all orders; hence `D = 1` and (F_g) for every `g >= 2`. Hyp. external
  inputs exactly [U] and [QRR]; everything else is proved in the note: [CYC], [FSP], [GAM], [HF], `thm:thp-a`, `thm:thp-b`, `thm:glu-main`,
  [GI], [PropS], `lem:pol-poly`, `lem:glu-wkb`. Src. SYN5 s.2 (S1-S7), s.3; RS2 s.6.2; MS5 s.1-s.2; CP s.0.3, CP3 s.3.11. St. PROVED modulo
  [U] and [QRR], every step AI-refereed. Corr. O8 (CP s.0.3's "complete list" omits (CYC), (S2), the relabelling, and includes R2, which
  is not needed); CP3 s.3.11's route through Lemma N is replaced by `thm:glu-main` (K3).
