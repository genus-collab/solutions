# DIGEST_B_5 (key read_B): s.3 continued (S20-S38), s.4 conventions, s.5 evidence

Continuation of `DIGEST_B_4.md` (marks [W], [W*], [S], [V], [M] defined there). Note: this piece also carries the end of s.3, which
did not fit in `DIGEST_B_4.md`; the piece table in `DIGEST_B.md` should be read with that shift.

## s.3 Proofs (continued)

**S20** (R0). 1. `Q^m(v(u))` solves the pulled-back model system (chain rule) and a regular gauge maps solutions to solutions: `P_can
Q^m(v(u))` is a formal solution of the quintic system [W] UR s.4; EV s.2.3 step 1. 2. Its exponentials: the model exponent at `v(u)` is
`-(5/4) zeta^i v^(4/5) = -(5/4) zeta^i u^(4/5)/tau = u_i - u_i^c` (D6, D7) [W] UR s.4. 3. The actual columns are WKB solutions of PF5 whose
rows obey the cyclic system with its closing relation (S10) [W] UR s.1; CYC Cor. 4. Formal WKB uniqueness (distinct `u_i`; transport
operator first order with kernel `Q L`) [W*] EV s.2.3 step 4; PS s.4.2 (b), C10. 5. The same exponential on both sides, so D is constant,
diagonal, in `Q[[z]]` [W] PS s.4.1. 6. Hatted: power series in u at fixed zhat since `mu(b,k) >= 0` [W] UR s.4; EV s.2.3 step 2.
7. Colour 0, row 0 [W] PB s.3. 8. The ring of comparison [W] PS s.4.2 (d). 9. `D_i(z) = D_0(zeta^-i z)`: S17 step 3 [S].

**S21** ((R1)). 1. `z P'` has no `u^-1` term; order `u^-1` of `A P - P At` is `r E_12 P(0) - rt P(0) E_12` [W] UR s.5; GI s.3. 2. Reading off
row 2 and column 1 [W] GI s.3.

**S22** (Lemma V). 1. Expand `e^(-x^4/(4z)) e^(-v/(xz))`; the residue picks `j = m + 4i + 1` [W] GI s.2; GR s.2. 2. (i) differentiation under
the residue; 3. (ii) residue of an exact form [W]. 4. The five rows of the model system (only (i), (ii) are used) [W] GI s.2; GR s.2.
5. `Y^C(0) = -e_1` [W].

**S23** (shape). 1. Monomials `v^j z^-n`: `n = 5i + c_b + 1`, `j = c_b + 4i + 1` [W] GI s.2. 2. `Q[z^-1][[u]]` since `v = O(u)` [W]. 3. Homogeneity
and the closed form with signs `(-1)^(m+i+j) = 1`, phase `zeta^(i(m-2))` [W] GH s.1-s.2; GR s.2.

**S24** (Lemma Q). 1. Grading by induction [W] GI s.2. 2. `ord y^(n+1)_a = ord b_a + ord y^(n)_(a+1) + 1 = (4(n+1) + c_a + 1)/5` using
`c_(a+1) = c_a - 1 + 5 delta_(a,1)` [W]. 3. Exactness of leading terms (integrands of valuation `>= 0`) [W] GI s.2; GR s.2. 4. Same leading
coefficients as `[z^-n] Y^C_a(v(u))` (`b_a`, `bt_a` both `bmod_a(u)(1 + O(u))`; `Y^C(v(u))` is the iterated-integral solution of the pulled-back
system by S25) [W*] GI s.2. 5. (b), (c): `z dY/du = A Y` coefficientwise [W]. 6. (d): `y^(1)_0 = -Int_0^u t_c' = -t_c` [W].
7. `[u^0](Y^C_Q)_0 = 0` [W*] PS s.3 (C7).

**S25** (Lemma 1). 1. Order `u^n`: `((n+1)z - B_(-1)) Y^(n+1) = Sum_(j<=n) B_j Y^(n-j)` [W] GI s.2. 2. The determinant is a unit of L,
`((n+1)z)^5` here [W] GI s.2; GR s.2. 3. Induction [W].

**S26** [W*] GI s.2 Remark; `W_2(0) = z` [W*] GR s.2.

**S27** (diagonal). 1. Order `z^0`: `b_a p_(0,a+1) = p_(0,a) bt_a` [W*] (UR s.3 (i) for `D_0`). 2. Order `z^1`, diagonal [W] GI s.5. 3. The
weighted sum cancels, `prod_a p_(0,a)` constant [W] GI s.5. 4. "`= 1` for `P_can`" needs `p_(0,a)(0) = 1` for all a, from `p_(0,1)(0) =
P_11(0,0) = 1` and `bt_a(0)/b_a(0) = 1` [S] (one line, not written). 5. `p_(0,0)^5 = prod_(j<4)(b_j/bt_j)^(4-j) = v/(u varpi_0^5)` using
`v'^5 = v/(u(1-u)^5)` [W] GR s.2 (GI states the result). 6. `p_(0,0) = tau^(-1/4)/varpi_0`, `Ghat(u,0)_00 = tau^(1/4) p_(0,0) = 1/varpi_0` [W*] GI s.5.

**S28** ((G)). 1. `Yt := Y^C(v(u),z)` solves `th Yt = At Yt` [W] GI s.3. 2. `P Yt` lies in `L[[u]]^5` [W] GI s.3; GR s.2. 3. `th(P Yt) =
(A P - P At) Yt + P At Yt = A(P Yt)` [W]. 4. `u^0` value `-P_11(0,z) e_1` by S21 [W]. 5. `P_11(0,z) Y^C_Q` has the same value [W]. 6. S25 with
`B = A` [W]. Cor. (a): `P_11(0,z) = Sum_q p_(5q,1)(0) z^(5q)` with `p_(0,1)(0) != 0` [W]; uniqueness of `P_can` among regular gauges with
`P Y^C = Y^C_Q` (every regular gauge is `P_can c(z)` by S13, and then `P Y^C = c(z) Y^C_Q`) [S] (implicit in GI s.3). (b), (c) [W*] GI s.3.

**S29** (symplecticity). 1. `S := P(-z)^T eta P(z)`: `z S' = At^T S - S At`, using the gauge equation at `z` and `-z` and `A^T eta = eta A`
(entry `(c,d)` of `eta A` is `5 b_(d-1) delta_(c+d=4)`, so this is `b_0 = b_2`, `b_3 = b_4`) [W] GI s.4; GR s.2. 2. `S = eta T`: `z T' = At T -
T At`, T regular, `T_0 = eta^-1 P_0^T eta P_0` invertible [W]. 3. Uniqueness clause of S13 for `(At, At)`: `T = t(z) 1` [W] (rests on S13
step 10 [S]). 4. Entry (2,1) at `u = 0` by S21: `5 P_22(0,-z) P_11(0,z) = 5 t(z)` [W]. 5. `c(z)c(-z) = 1`, `c` in `Q[[z^5]]`: `log c` odd,
`c = +-exp(a_1 z^5 + a_3 z^15 + ..)` [W*] GI s.4.

**S30** ((G*)) [W] GI s.4: (G) at `-z`, multiplied by `eta P(z)`, with S29. **S31** [W] GI s.5 (`u^n`-coefficients; explicit `n <= 4`);
GR s.2 (re-derived, incl. `(Ychi_4)_2 = (theta_4)_2 = -1/(24 z^4)`); values [V].
**S32** (R2) order `u^0`, entry (2,1) [W] UR s.5. (R3) "the same orders on (2,0), (3,1) and order `u^1` on (2,1)" [W*] UR s.5; a second
derivation for `P_can` from `(G)_1, (G)_2, (G*)_1` [W] GI s.5. The hatted sign [W*] NFS s.2 item 8; the explicit Fourier dictionary is code
(VT `v3_fourier.py`) [V].
**S33** (Theorem D). 1. Entry `(b+k,b)` of `[zhat^k] Ghat` is `u^(mu(b,k)+eps)` times a power series [W] UR s.5. 2. `eps >= 2` for (2,1),
`eps >= 1` if `a = 2` or `b = 1` (`5 !| k`), by (R1), (R2) [W]. 3. Table: first order `n(k) = 4q + (0,1,2,2,3)_r` [W]. 4. `Psi(u)` is a power
series [W*] [conv: OU3 s.1]. Sharpness [V].
**S34** derivation [M]; values [V] `ur_chain_check.log` (X1), VT `v2_run.log`.
**S35** (a) [W] UR s.4 (given (H)). (b) [W] RN s.3(a), but it rests on RN's (R-eq), (M-eq), (G-eq), whose derivation is not written
("referee's derivation; VERIFIED coefficient by coefficient") [V]. (c) [V] RN s.3(c). (d) [W*] RN s.3(d).
**S36** (a) [W] GI s.6 (given (H)); not circular (GR s.2: the rows fix their free coefficients at `u = 1`, the canonical gauge is
normalised at `u = 0`). (b) [W*] ST2b s.6.1.
**S37** (i) Apply `Lambda` to (G); the regrouping `(a, z^m u^n) -> (a, zhat^m, u^(n + (4m+4-c_a)/5))` is bijective (needs `m = a - 1 mod 5`)
[W] GH s.1; GR s.2. (ii) [W] GH s.2. (iii) difference quotient and the sign `(-1)^m` from `Q'(-psi)` [W] GH s.3. (iv) `h = 0` [W] GH s.4.
(v) `[x^-1](Q' f)_0 = -t_c` (S24 (d)), unit row `zeta^i/xhat` (S27), `[x^-1] f_i = -u tau^(-3/2) zeta^(-i)/5` [W] GH s.4; GR s.2.
**S38** (b) [W] UR s.6.1; (c) [W] UR s.6.3 given primary-sector inputs that are [V] for listed `(h, n)` only; (e) [W] UR s.8.

## s.4 Conventions, sign dictionaries, normalisations

**s.4.1 Fixed by the sources.** `D = q d/dq = -(1-u) d/du` (also in CYC); `th = z d/du`; `zD = -(1-u) th`; `b_a = -c_a/(1-u)` (UR s.1).
Index `a` in Z/5 = {0..4}; `A = Sum_a b_a E_(a,a+1)`; graded `P_k = Sum_b p_(k,b) E_(b+k,b)`; entry `(a,b)` of P is the bivector
`H^(3-a) (x) H^b` (UR s.5); components lowered, `Y_a = eta(phi_a, J)`, and P on Y is equivalent to `eta^-1 P eta` on J (TH s.1.1).
`eta_ab = 5 delta_(a+b=3)`, `eta^-1 = eta/25`; symplecticity written `P(-z)^T eta P(z) = eta` (GI) or `P(z) eta^-1 P(-z)^T = eta^-1` (VT v4).
Residues `r = rt = -1`, `B_(-1) = r E_12`. Normalisations: `varpi_0(0) = 1`, `t_c = u + O(u^2)` (programme family); `v = u + O(u^2)`;
`F_0(0) = 1`; `P_11(0,z) = 1` (canonical); `c(z) := P_11(0,z)` in (G) and S29; `M_cyc` is independent of the scales of `y_0` and g.
**s.4.2 Places where the sources differ.**
* Model normalisation: UR, GI, TH `lambda = -(x^4/4 + v/x)`; LG `lambda = x^4 + u_LG/x`, `u_LG = 4u`, `z_LG = -4z`, `eps = -4 zeta^-i zhat`;
  RN s.5 "`u_LG = c u`, `c^4 = -1/4`": the same statement, only `eps = -4 zeta^-i zhat` is fixed by the data (NFS s.2 item 7).
* Orientation: `tau zhat = -(5/4) zeta^i z/(u_i - u_i^c)` (UR s.4) vs `(5/4) z/Int_q^(q_c) L dq/q` (RN s.3(a)); immaterial (RN).
* Sign of `c_5`: unified as `P_stored = P_can(1 + m5 z^5 + ..)`, `m5 = -43397/807975000` (S01 C11; MEC s.6). UR s.4's `c_5 = +43397/807975000`
  is the constant with `P_actual = P_stored(1 + c_5 z^5)`; GI s.4 writes `2c_10 - c_5^2` with `c = P^s_11(0,z)` (its `c_5 = m5`), GI
  s.5-s.6 use `+43397/807975000` (GR s.4 item 6). The "stored" normalisation (`p_(5q,0)(0) = 0`) is code bookkeeping: keep it out of the
  note, or name it by its defining condition in the evidence appendix.
* Fourier labels: RN s.2 `(a,b) = (-m-1, -m'-1)`; GR s.3 (rf_actual) `Xhat(m,m') = (1/5) eps_a eps_b [Ghat(u,0)^-1 Ghat]_ab`, `a = m-1`,
  `b = -m'-1`. The two printed dictionaries differ in the first index; the code dictionary is VT `v3_fourier.py` (not read). Only the
  evidence appendix needs it.
* Euler-dilaton signs: flat `-pi z (uH + z)(x)(uH - z)` vs hatted `alpha zhat (H - zhat)(x)(H + zhat)`: `N_1 = -1` (NFS s.2 item 8), not an
  ordering convention as UR s.5 says.
* `H` and `e`: of record `H = d/dv`, `e = tau_1(1) + v tau_0(H)`, giving `<H>_1 = -1/12`, `<e>_1 = +1/8`; UR s.6.3 has `+1/12`, `-1/8`
  (the frozen-point H is `-d/dv` at `t = -1 <-> v = -1`; S01 C3).
* `c_g`: four forms, equal for `g = 2..12` (S01 C2); UR s.6.5 uses `(-1)^(g-1) B_2g/(2g(2g-2)) 5^(1-g)`.
* Pivot sign of S18: `-4n` (UR) vs `4n` (PS s.1 (v)); immaterial.
* Frame labels: programme `(alpha, beta)` vs the critic's `(a, b)` (D4; RS1 A2).
**s.4.3 Notation clashes to resolve in NOTATION.md** (source letters; proposals in brackets).
`c`: hatted exponents `(0,-1,3,2,1)` vs cyclic coefficients `c_a` vs `c(z)` vs `c_g` vs LD's `c = (-1,1,-1,-3,-2)` (ST2 s.2 step 7.2) vs `c_s(n)`.
`D`: `q d/dq` vs `D(z)`, `D_i` vs UR's `D_0 = diag(d_a)` vs "Theorem D" vs D(ii), D(iii) of Stage 2. `tau`: `tau(u)` vs CYC's automorphism
vs CYC's `tau_k` vs descendants `tau_k(.)` [`sigma` for the automorphism, CF3 E6]. `M`: `zD` vs CYC's factorised operator [`M_cyc`] vs
`M = M_0(tau zhat)`. `P`: the gauge vs `P(D)` vs the jet matrix [`Pi_ab`] vs UR's PF operator `P`. "Prop. S": UR (model rows) vs GI
(symplecticity). "Theorem G": GI vs the Gap Theorem [write "(G)"]. "Lemma 1": UR s.3 (iii) vs GI s.2; "Lemma 2": UR s.3 (iv) vs POLY.
`(R1)`: UR residue relation vs LG s.3 model relation vs R0-R5 of Stage 4. `eta`: pairing vs GI's `eta_n(z)` [`Ychi_n`]. `B`: `prod bt_a`
vs the concomitant vs RN's `B(u)`. `S`: solution space vs GI's `S = P(-z)^T eta P(z)`. `K`: `O[1/u]` / `Q((u))` vs the slice K vs RN's frame
K vs `K_q`. `L`: `1/s` vs `Q((z))` vs the TR leg class. `y`: `At^-1` vs solutions `y_i` vs `y^(n)` vs `y_b(t)`. `f`: Frobenius basis vs
CYC's `y_i/y_0` [`phi_i`] vs GH's `f(x)` vs the identity row `f(zhat)`. `g`: vanishing period vs genus. `h`: QRR series vs the power
series in `y_3` vs UR's `h_j` vs genus h. `e`: basis vectors vs `e = (0,0,1,1,1)` vs `e(z)` vs `e = Dt` (CF) vs UR's e, f. `Pi`: UR's
lift vs the kinematic plane vs `Pi_i` vs GH's `Pi(T')`. `w`: zhat (PB) vs weights.

## s.5 Evidence for the appendix (finite checks; none is used by a proof)

* **E1 (CYC).** pb_cyc: exact Q, every coefficient below `u^57`, programme frames (0,0), (1/3,-2/7), (7/11,-3/5) (PB2 s.4). Critic
  cr5_cyc_run: exact Q to `u^70`, each lemma separately, frames programme (0,0), (1/3,-109/210), (7/11,-23/22) (CYC s.8; RS2 s.5 item 3).
  CF: `M = PF4` below `u^77` in five frames, identically in `(alpha, beta)` below `u^27`, all 16 entries of C3 below `u^75`, B constant
  below `u^78`, `det Wr = -2/u^2` (CF3 s.10). CA: 17 frames to `u^74` (incl. `y_0(0) = 0` to `u^68`); both directions of S10 on actual
  objects (WKB rows `z^0..z^12` to `u^52`; `Y^C_Q` solves PF5 for `n <= 16`); `Cf(kinv[a]) = 1/c_a` mod `2^89-1` through `s^127`, frames
  (0,0), (1/3,-2/7) (CA s.2.2, CA2 s.6). Critic rc_cyc: `u^70`, 10 frames incl. `y_0(0) = 0`, identities on `u^x` at six non-integer x
  (RS1 s.0). independent_check/cyc_check: `u^60`, 4 frames (MS5 s.4.2). ur_model_check (A): closing relation on the s-series rows, all
  colours, `z^1..z^6`, `>= 32` s-orders, p89, p61 (UR s.1). Failing controls in every code (non-isotropic pairs, wrong gamma, `PF4 + uD`).
* **E2 (Theorem U).** ur_gauge_run: exact Q, `k <= 16`, every `p_(k,b)` a power series, resonances only at `k = 0 mod 5`, residues
  `(0,1/5,-3/5,-2/5,-1/5)` (UR s.3 Rmk (d)). VT v2_run: proof-route re-implementation, K = 8, `u^0..u^9`, equal to ur_gauge (log; NFS
  item 2 says `k <= 10`). Critic: K = 14, periods only, all holomorphy assertions passing (NFS s.5). GR rf_frob_run: an independent route
  (u-recursion with unknown `P(0,z)`, 1574 linear equations over Q): solution space exactly `c(z^5)` below `z^24`, canonical `P(0,z)`
  equal to ur_gauge's, both frames (GR s.3). PB2: the construction through `z^80` at p61, p89, every resonance unobstructed (PB2 s.4).
* **E3 (R1)-(R3).** ur_chain_check: `k <= 14`, exact Q (UR s.5); GI gi_relations: `(G)_n`, `(G*)_n`, `n <= 8`, `z <= 20`, two frames (GI s.5).
* **E4 ((G)).** gi_check: exact Q, 176 positive-z coefficients vanish (152 genuine cancellations, up to 142 terms), 1038 coefficients equal
  `Y^C_Q`, `z^-38..z^20`, `u < 29`, two frames; the stored normalisation fails in 127 of 176 (GI s.0, s.7). GR rf_gi: direct substitution
  (2940 coefficients), 1214 slots for both normalisations, controls (GR s.3).
* **E5 (symplecticity).** gi_relations: canonical gauge symplectic `z^0..z^20` mod `u^29`; stored defect exactly `(c(z)c(-z) - 1) eta^-1`,
  first at `z^10`, `2c_10 - c_5^2 = -18223115111903/15543419062500000` (GI s.4). VT v4_symplectic: `z^1..z^9` vanish identically in u,
  symplectic `c_10 = 153074543599907/261129440250000000` (log). Critic: canonical `c_10` = symplectic `c_10` (NFS s.5).
* **E6 (Lemmas V, Q; diagonal).** gi_formulas: valuations `n <= 45`, `p_(0,0)` to `u^37`, closed form of phi `m <= 80` (GI s.7); GR rf_misc,
  rf_gi (`n <= 35`); CA cy_rows (Lemma Q (a), (d), three frames).
* **E7 (actual = canonical).** gi_actual_norm: 1058 slots `zhat <= 24`, `u < 18`, p89, p61, two frames; `c^act = 1 + O(z^25)` (GI s.6). GR
  rf_actual: exact over Q on 397 Fourier components; 2250 graded slots mod p with `M_0` from the model recursion (GR s.3). UR: one
  constant `c_5` at both primes (UR s.4). Later windows: `c^act = 1` through `z^35` (ev_direct), `z^150` (poly referee rf_gauge), `D_0 = 1`
  exact through `z^204` given (U) (MS4 s.1, s.4.1). Superseded by the proof of `D = 1`; keep as evidence.
* **E8 (model rows).** ur_model_check (B): symbol = model, all rows, colours, `n <= 6`, both primes; stationary phase = `F_0` through
  `t^10` exact (UR s.2). VT v1_model: 175 entries, both primes. LG: all rows `zhat^18`, identity row `zhat^40` (NFS item 4). pb_s2 /
  PS T4: (S2) through `t^60` / `t^80`.
* **E9 (Corollary F, degrees, rule W).** RN: `zhat <= 16`, `u <= 11`, three primes; `zhat <= 24`, `u <= 17` at p89; degrees
  `floor((5n+2)/4)` for `n <= 17`; rule W in 120 pairs (RN s.2); 7 frames at `zhat <= 12`, `u <= 8` (NFS item 13).
* **E10 (genus 2, 3 with M_nat).** 41 genus-3 edge graphs `O(u^4)`, 6 genus-2 edge graphs `O(u^2)`, frames (0,0) at two primes and
  (1/3,-2/7) at p89, `M` mod `(u^4, zhat^7)` (UR s.6.2; NFS s.4 wording). Genus 4-5 evidence: other digests.
* **E11 (GH).** 2160 edge coefficients (965 with `>= 2` non-zero terms), control without `(-1)^m` 880 non-zero, `zhat <= 24`, `u < 18`,
  two primes, two frames; production rows `zhat <= 6`, `u < 6` (GH s.6). GR rf_hat: flat frame, exact Q, (T) through `u^11`, leading `-1/30`
  in both frames (GR s.3).
* **E12 (Euler-dilaton).** ur_euler_linear (genus 2, 3, both primes); the negative result ur_euler_check (p89) (UR s.6.3).
