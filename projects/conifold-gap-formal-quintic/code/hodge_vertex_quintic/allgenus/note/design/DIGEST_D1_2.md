# Digest D1, part 2: s.2 statements

Continuation of `DIGEST_D1.md` (abbreviations, s.0, s.1). Proposed label prefixes: `nrm` normalisation problem, `pol` polynomial
structure, `fin` finite route, `thp` Theorem P, `glu` the glue. Format: **`label`** (type) *name.* Statement. Then: Hyp. (hypotheses
and scope), Src. (source), St. (status in the sources), Corr. (corrections applied), Dep. (dependencies).

## 2.0 Inputs from other digests (statements as the chain uses them; their labels are set by those digests)

* **[U]** Theorem U (UR s.3): two cyclic systems with one simple pole and equal products `Prod b = Prod bt` have a regular gauge
  (D1.7), unique up to `P -> P c(z)` (graded: `c(z^5)`). PROVED in the project, AI-audited only (MS5 s.2).
* **[R0]** UR s.4: `Q D = P_can Q^m(v(u))` (D1.8). PROVED given [U]; (CYC) supplies that the actual rows solve the system Theorem U is
  applied to (RS2 item 13).
* **[GI]** GI s.2-s.3: Lemma V, Lemma Q, Lemma 1 (a regular solution of `z dY/du = B Y` in `L[[u]]^5`, `L = Q((z))`, with `Y(0) = 0` is 0),
  (R1), Theorem G: for every regular gauge `P`, `P(u,z) Y^C(v(u),z) = P_11(0,z) Y^C_Q(u,z)` in `L[[u]]^5`. PROVED (given [U]), AI-refereed.
* **[PropS]** UR s.2 Prop. S, uniqueness half: the model T-system has a unique solution with `F_0(0) = 1` (pivot `4n`); this half does not
  use UR's (H) (PS s.1 (v)).
* **[CYC]** In every conifold frame, `PF4 = (1/varpi_0) D (1/Dt) D (u varpi_0^2 (Dt)^2) D (1/Dt) D (1/varpi_0)`, so the cyclic system with
  its closing relation `4 -> 0` is equivalent to PF5. PROVED, refereed twice; scope widened (A1): any power-series `varpi_0` not
  proportional to `g`; `varpi_0(0) != 0` is used only downstream (unit shape for [U]; `C0(0) = 1`, `t_c'(0) = 1`). (MS5 s.3.4; RS s.3.)
* **[FSP]** formal SP lemmas in FSL numbering: FL1 (FSL2 s.3.1), Lemma T(b) (FSL2 s.3.2), L1c(iii) and FL4 (FSL2 s.4, s.4.5), L1a
  (FSL2 s.5) with s.5.4 (CFU termwise in `u`), L3(iv) (FSL4 s.9.1-s.9.2), the ring `L` and FL2, FL3 for the residue (FSL4 s.10.1-s.10.3);
  Lemma H, z-rescaling (RL s.3.2). PROVED, refereed (RL s.1; SYN5 S1).
* **[HF]** hatted/flat identity: `c_n(L) = L^(4n) ph^SP_n(L^-5)` as rational functions, so `c_n(1) = ph^SP_n(1)` and `c(z; 1) = Phi^SP(z, 1)`.
  PROVED (CP3 s.3.10 with RF2 C2: the Jacobian `-L/X^2` belongs to the hatted amplitude; inputs L1a for `sigma = L/X`, Lemma H with
  `mu = L^-4`, L1c(i), and P(a)). SYN5 S3.
* **[GAM]** Theorem Gamma (iii): every `c_n` is regular at `L = 1` and `c(z; 1) = h(z)`. PROVED twice, each refereed twice (MS5 s.1).
* **[LMI0]** LMI at `s* = 0`, `alpha0 = 5`, `M0 = 1`: `c(z; 0) = Ginf(z)` (RA s.3.4; LF s.4.2, with the formal DLMF 5.11.8). PROVED,
  refereed twice (RS s.4).
* **[QRR]** used only as H9: `A^QRR_n(q = 0) = h_n` for every `n` (PB s.1; PS s.4.3; MS5 s.2).

## 2.1 The normalisation problem

* **`defn:nrm-D`** (defn) *the normalisation constant.* D1.8: R0 defines `D`; (F_g), (H), (N). Src. UR s.4; SYN4 header; MS4 s.1.
  St. R0 PROVED given [U] (SYN4 s.1.2 R0). Corr. PS-C8 (strip the same exponential). Dep. [U], [R0], [CYC].
* **`rmk:nrm-genus`** (rmk) *genus bookkeeping.* `D_0 = 1 mod z^(N+1)` iff (H) for `5 ∤ n <= N` and `K_q = 0` for `5q <= N`; (F_g) needs
  `N = 3g - 3`. With `log c = Sum_j a_j z^(5j)` (odd `j` independent), `a_j` first matters at genus `ceil((5j+3)/3)` (`a_5` at 10, `a_7`
  at 13, `a_19` at 33). Symplecticity forces `D_0(z) D_0(-z) = 1`: only odd `a_j` are independent. Src. EVID s.2.4; PR C2 (from GI
  Prop. S); GP s.8. St. PROVED (EVID, given Lemma N); the symplectic relation is stated by PR C2 (I did not read GI s.4).
* **`rmk:nrm-model`** (rmk) *the conifold model does not fix the normalisation:* `ph_n(0) = F_(0,n)` in every normalisation
  (`lem:pol-u0`); the free coefficient is the top coefficient `f_m`, the value at `q = oo`. Src. SYN4 s.0, s.1.1; MS4 s.3. St. PROVED.
  Corr. Stage-4 O17 (the brief's reading is refuted, not `D = 1`).

## 2.2 The polynomial structure

* **`lem:pol-HPF`** (lem) *HPF.* Put `A_e := 1 + w(1-u)((4/5) th_w - e/5 - th_u)`, `Phi_0 := Phi`, `Phi_j := A_(-j) Phi_(j-1)`. For
  `psi = L e^(u_0/z) Phi(w, u)`, `PF5 psi = 0` iff
  `(HPF) Phi_5 - Phi_0 - (1-u)[3w Phi_4 + (17/5) u w^2 Phi_3 + (9/5) u^2 w^3 Phi_2 + (274/625) u^3 w^4 Phi_1 + (24/625) u^4 w^5 Phi_0] = 0`.
  Every raising carries a left factor `(1-u)`, and `[w^(n+1)] HPF = (1-u)[(4n - 5 th_u) ph_n + r_n]` with `r_n` built from
  `ph_0..ph_(n-1)` (explicit formula POLY s.1). Hence `[u^k] ph_n = [u^k] r_n/(5k - 4n)` for `5k != 4n`; for `n = 5m` the resonance
  `[u^(4m)] r_n = 0` is required and `f_m` is free. Src. POLY s.1 (via `(zD + L)(s^e Fh) = s^(e-1) A_e Fh`); EVID s.2.1; PB s.2.2. St. PROVED;
  re-derived by PR (by hand), ER (Laurent monomials), GR (flat `L`), the Stage-4 critic (SYN4 s.1.1).
* **`lem:pol-u0`** (lem; POLY Lemma 1) *the conifold value.* For every `n >= 1` and every choice of the free coefficients,
  `ph_n(0) = F_(0,n)`. Src. POLY s.2; EVID s.2.1 (a); GAM s.1; GR s.2.1. St. PROVED three times. Corr. ER s.7 item 5: add the line
  identifying the `u^0` recursion with the model (scaling limit `v th^5 + 3z th^4 + 1`, unique homogeneous WKB solution, Prop. S).
* **`lem:pol-poly`** (lem; POLY Lemma 2) *power-series solutions are polynomial.* If `Phi` solves HPF with `ph_0 = 1` and every
  `[w^n] Phi` in `Q[[u]]`, then every `ph_n` is a polynomial in `u`, and the set of such solutions is `{Phi^(0) c(z) : c in 1 + z^5 Q[[z^5]]}`
  (`z^(5m) = u^(4m) w^(5m)`); multiplication by a function of `z` commutes with HPF (`((4/5) th_w - th_u) z^m = 0`). Hyp. the resonance
  conditions hold for every choice of the `f_m` once one power-series solution exists: P(a) supplies one without [U] (PR s.2), or
  [U]; independently `prop:pol-deg`. Src. POLY s.2; PR s.2; EVID s.2.1 (b); PB s.2.2. St. PROVED.
* **`prop:pol-deg`** (prop; POLY Lemma 3 = ER result B) *degree bound.* For every `n`, every colour and every choice of the free
  coefficients, `deg ph_n <= floor(4n/5)`, and every resonance is trivial (`deg rho_n < 4n/5` for `5 | n`). Equivalently `A_n` lies in
  `L Q[L]` with degree `<= 4n + 1`; the flat amplitude lies in `Q[L][[z]]` and its value at `L = 0` is `1 + Sum_m f_m z^(5m)`; `f_m` is the
  top coefficient of `ph_(5m)`, `f_m = R_(5m)(0)`. Equality `deg ph_n = floor(4n/5)` for `n = 1` and `n >= 3` (`deg ph_2 = 0`) is VERIFIED
  only and must not be claimed. Src. ER s.4 (elementary proof); GAM s.1, GR s.2.1, PR C1 (via arXiv:0710.0889 Thm 4(ii)). St. PROVED twice
  (Stage-4 K2). Corr. Stage-4 O3 (EVID, POLY called it VERIFIED: underclaim), O4 (EVID's interval `(-4n/5, 0]` is `[-4n/5, 0]`), GR C1
  (the 0710.0889 route must not presuppose (H)). Dep. `lem:pol-HPF`; for the 0710.0889 route also `lem:glu-wkb` (PR C1).
* **`lem:pol-Finf`** (lem; POLY Lemma 4) `Phi^QRR = Phi^(0) F_inf^QRR(z)` with `F_inf^QRR(z) = h(z)/Phi^(0)(z, 1)`; so (N) iff
  `F_inf^QRR = F_inf^can := 1/c^(0)(z)` (`c^(0)` = the `c` of `Phi^(0)`). Src. POLY s.2-s.3; PR s.2. St. PROVED. Optional in the note.
* **`lem:pol-p5k`** (lem) `p_m := Sum_(j=1..4) (1 - zeta^j)^(-m) = 0` for every odd multiple `m` of 5. Hence, for odd `k`,
  `[z^(5k)] log h = B_(5k+1)/(5k(5k+1) 5^(5k))` (only the twist weight contributes), and `log h` has no `z^(5k)` term for even `k`.
  Src. EVID s.8.2; ER s.2 (`y_(-j) = 1 - y_j = -zeta^j y_j`, `y_j = 1/(1 - zeta^j)`). St. PROVED. Optional (it explains control B1).
* **`prop:pol-Ginf`** (prop) *free coefficients in closed form.* `Sum_m f^can_m z^(5m) = Ginf(z)`; in particular
  `f_1 = B_6(1/5)/6 = 4537/3937500`; `f^QRR = f^can` follows from `D = 1`. Also `c_n` lies in `Q[L]` (exponents in `4n + 5Z`) without the
  degree bound. Src. RA s.3.4; LF s.4.2-s.4.3; SYN5 s.3; SYN5b K7, O14; RS s.4. St. `c(z; 0) = Ginf` PROVED, refereed twice;
  `f^can_m = [z^(5m)] Ginf` PROVED given P(b), refereed once (LF s.4.3). Corr. CONJECTURAL in SYN4 s.1.1 and MS4 s.3: superseded; RA's
  "given the inputs of Theorem D" is more than needed (O14). Dep. [LMI0], [HF], P(a) (no negative powers of `u`), `thm:thp-b`.
  Prior art: GJR arXiv:1812.11908 Lemma 7.7, eq. (37), for their LG integral (MS5 s.3.6).
* **`rmk:pol-points`** (rmk) *three points.* `q = 0` (`L = 1`; value `h`, Gamma factors `Gamma(1 + 5/z) prod_(j>=1) Gamma((zeta^j - 1)/z)`),
  `q = oo` (`L = 0`; value `Ginf`, factors `prod_j Gamma(zeta^j/z + 1/5)`), conifold (`L = oo`; top coefficients, canonical gauge).
  Src. GAM s.0 table; EVID s.8.3 (SKETCH there), made rigorous by ER s.4; ER s.7 item 7.

## 2.3 The finite route (not used by the proof of the main theorem)

* **`lem:fin-P`** (lem; Lemma P) For every covariant `Qhat` (`Qh^(i)_(a,n) = zeta^(i(a+1-n)) Qh^(0)_(a,n)`):
  `[w^N](Qhat R^(tau w)^-1)_ac = [N = a - c mod 5] Sum_(n<=N) Qh^(0)_(a,n)(u) (-tau)^(N-n) F_(3-c, N-n)`. Hyp. model symplecticity
  `M_0(t)^-1 = M_0(-t)^T`, `pinv^T pinv = (1/5) eta'` (`eta'_(bb') = [b + b' = 3 mod 5]`), so `R^(t)^-1 = (1/5) R^(-t)^T eta'`; Prop. S;
  `Sum_i zeta^(im) = 5 [5 | m]`. Src. EVID s.2.2; ER s.2. St. PROVED; VERIFIED. Not used by the all-orders chain (RF2 s.5).
* **`lem:fin-N`** (lem; Lemma N) Take any choice of the free coefficients, `Phi` the resulting polynomial-type colour-0 solution, `Qhat`
  its hatted rows (rows by the cyclic system, colours by covariance), `Ghat := Qhat R^(tau w)^-1`. Then `Ghat` is the hatted form of
  `P_can c(z)` for a unique `c` in `1 + z^5 Q[[z^5]]`; for that gauge `P_11(0, z) = c(z)`; `[u^j w^(5q)] Ghat_11 = 0` for `j < 4q`; and
  `K_q = [z^(5q)] c = [u^(4q) w^(5q)] Ghat_22`. Slope `dK_q/df_q = F_2(0)/(C1(0) C0(0)) = 1`. Hyp. [U], (R1), `lem:glu-wkb`, [CYC] (PS-C9).
  Src. EVID s.2.3; ER s.2; PS s.4.4. St. PROVED given [U]. Corr. PS-C9 (reading; add (CYC)); K3 (only steps 1-5 are needed for the
  all-orders statement; the formula for `K_q` serves the finite route).
* **`thm:fin-conn`** (thm; connection form) Given [U] and [QRR]: (i) (H) holds iff `ph_n(1) = h_n` for every `n` with `5 ∤ n` (the unique
  polynomial solution; no free parameter); (ii) given (H), (N) holds iff `K_q = 0` for every `q >= 1`, and `K_q = 0` fixes `f^can_q`; (iii)
  `D = 1` iff `ph^can_n(1) = h_n` for every `n`, without assuming (H). Src. EVID s.2.4; ER s.2 (iii without (H)). St. PROVED given [U] and
  [QRR]. Corr. RF2 C4 ("the u = 1 convention" is [QRR] in its flat colour-0 form); K3 (superseded by `thm:glu-main` for all orders).
* **`prop:fin-pfree`** (prop; ER result A) For every polynomial-type colour-0 `Phi`: `S_q = O(u^(4q))` and `K_q = -[u^(4q)] S_q`; the periods
  `C0`, `C1` drop out. Consequence: given [U], `K_q` is frame independent, so further conifold frames give no independent evidence for
  (N). Src. ER s.3. St. PROVED given Lemma N and (R1); VERIFIED. Corr. Stage-4 K7, O2.
* **`prop:fin-row`** (prop) Row 0 of `z P' = A P - P At` at `u = 0`: `P_11(0,z) = P_00(0,z) + z (d/du) P_01(0,z)`, so
  `[z^(5q)] P_11(0,z) = p_(5q,0)(0) + p'_(5q-1,1)(0) = [u^(4q)](g_(5q) + g_(5q-1))`; (N) iff this vanishes for all `q >= 1` on the QRR
  polynomials (slope 1 in `f_q`). `m5 = p'_(4,1)(0) = -43397/807975000 = -p^can_(5,0)(0)`; the stored normalisation has
  `P_11(0,z) = 1 + m5 z^5 + ..`. Hyp. `b_0(0) = t_c'(0) = 1`, `bt_0(0) = v'(0) = 1` (fq_conifold family, RS2 item 28), (R1). Src. POLY s.3;
  PR s.2. St. PROVED. Corr. PR C7 (Stage-4 O15): divisibility of `g_N` by `u^mu` is necessary, not sufficient, for regularity.
* **`rmk:fin-routes`** (rmk) The three formulas for `[z^(5q)] c` are three readings of the one scalar `[z^(5q)] P_11(0, z)`. Src. SYN4
  s.3 verdict (agreement VERIFIED `q <= 12`, check C1); Stage-4 K9 (two perturbation rules, both predicted by `K_q = [z^(5q)] Phi(z,1)/h(z)`).
