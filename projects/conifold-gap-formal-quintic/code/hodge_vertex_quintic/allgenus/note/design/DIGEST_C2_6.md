# DIGEST_C2, part 6: s.4 Conventions, sign dictionaries, normalisations; s.5 Evidence for the appendix

Continues `DIGEST_C2_5.md`.

## s.4 Conventions (as the sources fix them; every place where two sources differ is flagged)

**C1. Orientation of the TR kernel (SYN2b K1; LR s.2; XYR s.2; TRM s.1).**
* Engine convention (TRM engines, LL, FB, LR, XY "code"): kernel `(1/2) Int_{sigma z}^{z} B/((y(z) - y(sigma z)) dx)`; DNOPS (2-28) is the same
  (numerator and denominator both flipped). SYN2b K1 calls this "the common convention".
* Printed EO (4-3) (its footnote 3) and ABDKS (1.1) with `omega_{1,0} = -y dx` are one and the same expression, opposite to the engines.
* Dictionary: `omega^engine_{g,n} = (-1)^n W^print_{g,n} = (-1)^n omega^ABDKS_{g,n}`; `F^engine_g = -F^print_g`; engine
  `omega_{0,3} = -Sum Res BBB/(dx dy)` against EO Thm 4.1 (4-20) `+Res BBB/(dx dy)`; for the temperature variation `Lambda^print = -1`.
* Conflict settled (K1): LL s.2.3 ("our kernel is EO's") and TRM s.1 / XYR s.2 ("minus the kernel of arXiv:math-ph/0702045") use the phrase
  for different objects; the logleg referee is right about the printed paper; no result is affected (signs fixed by regressions: Airy
  intersection numbers, `F_2^TR = -1/240`, the (0,2) -> (0,3) case, tensor-level variational checks). The critic checked this "only on
  the formulas as quoted in the notes, not against the PDF". [DW] Printed EO Thm 4.7 (4-30) has `(2g + k - 2) W_k`; with
  `omega_n = (-1)^n W_n` it becomes `(2 - 2g - n) omega_n`, the engine dilaton equation of TRM s.3: consistent.

**C2. Free energies.** `F_h^TR` (engine) `= B_2h/(2h(2h-2)) v^(2-2h)`, `F_2^TR = -1/240` (= `chi(M_2)`); TRM s.4.1 writes "EO/DNOPS conventions"
meaning C1's common convention; TRM s.1: Hock's kernel and `F_g` (1.1) are "identical to ours" (LR F5: not checked by the referee).
Frozen: `F_h = (-5)^(1-h) F_h^TR = c_h v^(2-2h)`; `F_0 = (-5)(1/2) v^2 log v = -(5/2) v^2 log v`; `F_1 = -(1/12) log v` in both.

**C3. TR versus frozen normalisation (TRM s.2.2-s.2.4).** `Lambda = -5 kappa` for `x = kappa W`; for `x = W`: TR pairing
`= (-5)^(h-1) x` frozen ancestor correlator, no psi rescaling. The factor `5^(g-1)` belongs to `x = -W` (psi -> -psi, up to `(-1)^(3g-3+n)`:
S01 s.7). DNOPS form `lambda = p^4 + v/p` at `v = 4`: `omega = (-20)^(g-1) Sum 4^-D <..>^frozen prod d xi`, `F_g = (-5/16)^(g-1) c_g`.
TRM s.0: for `kappa = 1` the TR CohFT is the frozen theory with metric times `-1/5`.

**C4. Metric signs (TRM s.1; MSG s.1).** DNOPS (2-7) `eta = -Sum_crit Res d'(lambda dp) d''(lambda dp)/d_p lambda` is MINUS the metric (1.6) of
arXiv:2211.12259 and of LG_MODEL s.5.1; the CohFTs differ by `Omega_g -> (-1)^(g-1) Omega_g`. MSG s.1: the frozen theory is the Hurwitz
Frobenius manifold `H_{0;1,4}` "with metric x 20 and z -> -4z". Flat pairing `eta(H^a, H^b) = 5 delta_{a+b=3}` (RD s.1.1; `eta^2 = 25`).

**C5. Leg dictionaries (K10; LD box; XYR s.4; S01 C7, C8).** Of record: LD (`tau_0(1) <-> +p`, step `-Int dx`, `q_{Phi_k^b} = -c_b tau_k(H^b)`,
`c = (-1,1,-1,-3,-2)`; `-log p <-> tau_0(H)`; `Int y dx <-> -tau_1(1)`; `e <-> -p^5/5`). XY s.0: `tau_0(1) <-> -z`, step `+Int dx` (both pairs
satisfy string and dilaton; only LD is dual to the `xif_{a,d}`; no XY result depends on it). XY prose "the unit acts by `-d/dx`" is wrong,
`+d/dx` (S01 C8). On `x = z^4 + v/z`: `e <-> -(4/5) z^5`; `psi e` leg `-(16/45) z^9 + (v/5) z^4` (S01 C7; XY's `-(2/5)` is a slip).

**C6. The letter c (clash).** LD constants `c_b = (-1,1,-1,-3,-2)`; grading `c_a = (0,-1,3,2,1)` (RD, S15, BV, GI); branch points `c_b = zeta^-b v^(1/5)`
(LL, FB); genus constants `c_h`; normalisation series `c(z)`, `c^act`, `c_5`, `c_10`; `c0 = P_11(0,0)`. Needs distinct symbols in NOTATION.md.

**C7. Weights and exponents (four conventions).** RD s.1.2 and S15 s.1.2: `e = 2 - 2h + Sum w(k,a)`, `w(k,a) = (4k + c_a - 4)/5`; FR: frozen
scaling `v^(2-2h) prod v^(w(k,a))`, same `w`; MSG s.1: `w(a,k) = (4 - c_a - 4k)/5`, `e = 2 - 2h - Sum w` (opposite sign of `w`, same `e`);
TRM s.3: `e = 2 - 2g + (4d - beta)/5`, `beta = (4,5,1,2,3) = 4 - c_a`; LL (F2), FB s.1: legs `p^(m_i)` and `k` log legs give
`e = 2 - 2h - n + Sum m_i/5`, Z/5 forces `Sum m_i = 0 mod 5`. DP: leg weights `w_b = (1,0,4,3,2) = c_b + 1`; `S_w` weight `w`. "Singular" =
`e < 0` and value non-zero; "regular" = `e >= 0` (RD s.1.2).

**C8. `Pi_i` and the vanishing period (MSG s.1; RD s.1.3; GI s.2).** `PiC(tau_n H^b) = (-1)^(n+1) [z^-n] Y^C_b`; `PiC = 1` on `H`, `-v` on
`tau_1(1)`, `-v^2/2` on `tau_2(H^4)`, `-v^3/6` on `tau_3(H^3)`, `v^4/24` on `tau_4(H^2)`, `v^4/96` on `tau_5(H)`, `-v^5/480` on `tau_6(1)`; zero unless
`b + k = 1 mod 5`. RD: `V^(k)(q) = (-1)^k [z^(k-1)] <Y^C(-z), q(z)> = Omega(z^-k J^C, q)`, `T_v = cap_k ker V^(k)`. [DW] On `tau_n(H^b)`,
`V^(k) = (-1)^(k-1) Pi_(k-1)` (from the two closed forms), so `T_v = K_0` and `V^(1) = Pi_0`. `Y^C = (I_0, -I_-1, -I_3, I_2, I_1)`,
`I^C_m = Res p^m exp(-(p^4/4 + v/p)/z) dp` (GI writes the variable as `x`); `Y^C(0,z) = -e_1`; leading terms
`(-v/z, -1, -v^4/(24 z^4), -v^3/(6 z^3), v^2/(2 z^2))`. Lowered components `Y_a = eta(phi_a, J)` (RD s.1.1).

**C9. Gauge normalisations (S01 C11; GI s.1).** Canonical: `P_11(0,z) = 1`; stored: `p_{5q,0}(0) = 0`; `P_stored = P_can (1 + m5 z^5 + ..)`,
`m5 = -43397/807975000` (GI s.4 used the opposite sign convention in one place; unified in C11). `d = (P_0)_33`, `d(0) = 1`.

**C10. Hatting and engine normalisation (S15 Lemma 3; FR s.1; BV header).** `u = s^5`, `zhat = z/s^4`; rows `Qhat_a = s^(1-c_a) Q_a`
(FR: `[zhat^n]` coefficient `= s^(4n+1-c_a)` times that of `Q`); `xhat = x/s`, `x = s varpi_0`; `d = tau^(1/4) varpi_0`, `xhat = varpi_0`;
`G_Gamma = u^(2g-2) F_Gamma`; regular iff `[u^0..u^(2g-3)] G_Gamma = 0`. Genus-3 tables index `k` = coefficient of `u^(k-4)`.

**C11. Ancestor potential and dilaton (RD s.1.2).** `Fb_h(v; q) = Sum_n <(q + z phi_0)^n>_h/n!`, base point `-z phi_0`;
`Fb_h(lambda q) = lambda^(2-2h) Fb_h(q)` (`h != 1`), `Fb_1(lambda q) = Fb_1(q) - (5/24) log lambda`. S15: `d^n F_h(q/kappa) = kappa^(2h-2+n) d^n F_h(q)`.

**C12. Homogeneity of TR (LL (F2); TRM s.5).** `p -> s p`, `v -> s^5 v`: `y dx -> s^5 y dx`, `omega_{h,n} -> s^(5(2-2h-n)) omega_{h,n}`;
with `s^5 = u`: `R_TR(z) = M_0(z/s^4) = M_0(zhat)`.

**C13. Other name clashes to resolve in NOTATION.md** (S01 C9 and this digest): `Lambda` (TR normalisation `-5`; DP leg map; EO's `Lambda(q)`;
RD's slice `Lambda_v`, `ML_Lambda`); `K` (slice; TR kernel; Stage 4 `K_q`); `M` (DP module; split matrix `M_nat`); `N` (DP module; condition (N);
BR's nilpotent dressing `N`); `q` (`v/p`; ancestor vector `q_Phi`; dilaton-shifted point; Kahler parameter); `x` (TR coordinate `W`; LG variable in
GI; `x = s varpi_0`); `D` (`d/dx`; the series `D(z)`; labels D(ii), D(iii)); `L` (class; Stage-5 parameter; number of leaves); `d` (`d(w)`;
`(P_0)_33`; descendant index); `e` (Euler-dilaton vector; `e(w)`; basis `e_{b,k}`; `e_1`); `S` (`S_w`; `S(z)`; `diag(s^(1-c_a))`; subsets);
`s` (`u^(1/5)`; elements of `S_w`; `s = d T_P - v`); "Theorem G" (gauge identity; gap theorem in TRM s.4.1); "Lemma C" (BR; RD_3 s.7);
"Lemma 1" (DP, S15, GI, RD_2); "(U)" (Theorem U; S15's unit hypothesis); "(ML1)/(ML2)" (MSG vs XY); (F1)/(F2) (LL properties) vs (F), (F_g).

## s.5 Evidence for the appendix (finite checks; none may be cited by a proof)

| # | what | range, fields, counts | source |
|---|---|---|---|
| E1 | D(iii) falsifier (a): log-free legs outside `L` | none: 10400 vectors (`k <= 80`, `w <= 320`), 16200 (`k <= 100`, `w <= 400`), exact over Q; re-run byte-identical (DR, and the main session in 12 min 28 s) | DP s.0, s.8; DR s.0, s.5; SYN2b s.3; MS2 s.4 |
| E2 | D(iii) steps T1-T7, X1-X4 | T1 `w <= 400`; T2 861 in / 859 out (`w <= 80`); T3 1042; T4 `k <= 101`; T5-T7 `w <= 400` (16200 raw-leg kernel vectors, O1); X1 505/505, X2 505, X3 231, X4 101; controls C1-C4 fail | DP s.7-s.8 |
| E3 | independent re-derivation of D(iii) | R0 `w <= 4000`; R1 `w <= 100` (1326 in, 1374 out, 606/606 random); R2 1602 elements; R3 `k <= 200` (Q), `k <= 1000` (2^61-1, 2^89-1); R5 `w <= 400` (Q), `<= 300` (primes); R6 16200 / 36300 / 100500 vectors (`w <= 400` Q; `<= 600`, `<= 1000` primes); R7; 13 failing controls; detection power: bad `H^2` leg gives 40 FAIL lines | DR s.4, s.5 |
| E4 | D(ii) | critic Test A `k <= 100` (Q); DR D1 `k <= 100`, D2 `n <= 2000` (Q) and `<= 3000` (primes), D3 `w <= 400` (Q), `<= 800` (primes); 4 controls (2402, 480, 80 weights, 390 of 401) | DR s.3; DP s.0 |
| E5 | LD | `cr2_legdict 24`: `b, a <= 4`, `k <= 24`, `d <= 27`, `v = 1, 3^5`, 35000 pairings per prime, 0 mismatches; D(ii) 125/125 (Q); 52/52 frames; controls C1 60, C2 25, C3 25, C4 1385 mismatches | SYN2 s.2 box, s.7; MS2 s.4 (re-run 74 s) |
| E6 | LD on actual TR correlators | 52 frame one-point functions, `g = 2, 3`, exact over `Q(zeta5)`: `<frame(3,3)>_2 = 3/3200`, `<frame(4,2)>_2 = -41/9600`, `<frame(6,0)>_3 = 3/128000`, `<frame(7,4)>_3 = 223/2688000`; all singular ones 0 | SYN2 s.2 box |
| E7 | TR = frozen model | 72120 colour-basis correlators (psi-deficit `<= 6`; (0,3)..(3,1)), 4665 non-zero, two primes; `F_2`, `F_3`; R-matrix `z^0..z^8` (225 entries, `Q(zeta)`) | TRM s.0 item 2, s.2.3 |
| E8 | falsifier (b) (worker) | (4,2), (3,3), (2,4), (1,5), (0,6) and all `2h - 2 + n <= 6`; 17416 singular ordered tuples per prime (4678 mixed), all 0; 18 pure-log values = Gap derivatives (new `<H^2>_4 = 7/30000`); exact `Q(zeta5)` for (0,3), (0,4), (0,5), (1,1), (1,2), (1,3), (2,1), (2,2), (3,1); vacuous at (0,3), (1,1) (F4) | FB s.0-s.3; SYN2b s.3 |
| E9 | falsifier (b) (referee, independent engine) | `g + n <= 7`, `g <= 6`: 40739 singular multisets per prime (19194 mixed), 0 non-zero; 26 pure-log values incl. `<H^3>_4 = -7/3750`, `<H>_5 = -1/82500`, `<H^2>_5 = 3/27500`, `<H>_6 = -691/102375000` (re-derived by hand from `c_h`, SYN2b s.7); `F_2, F_3, F_4 = -1/240, 1/1008, -1/1440` from the definition; control: log leg based at `p = 2` | LR s.1 |
| E10 | variational formula and (V1) | tensor identity 14 cases (0,3)..(4,1) (worker), 7 cases (referee); (V1) on 32559 tuples (32 triples, `h <= 4`) and 21733 (referee); wrong-sign / no-xi controls fail (467 / 1920); (D) dilaton on the same tuples | LL2 s.4; LR s.1 |
| E11 | Gap Theorem | TR `F_2..F_5` at two primes; formula (1.3) exactly `g <= 12`, several polynomial parts; generic points of `H_{0;1,4}`; family `x = p^r/r + 1/p`, `r = 1..6`; Borel identity `g = 1..10` | TRM s.4.2; XYR s.2 |
| E12 | S15 at genus 3 (worker) | 42 graphs x `u^-4..u^3` = 336 coefficients per prime and frame; 108 non-zero edge-graph + 4 vertex coefficients (K4); frames (0,0) at 2^89-1 and 2^61-1, (1/3,-2/7) at 2^89-1; all 41 edge graphs regular in the flat frame; vertex `(1/25200) t_c^-4` + regular; genus 2: 7 graphs x 4 orders | S15 s.0, s.6; SYN2b s.3 |
| E13 | S15 at genus 3 (referee) | independent red/blue engine; R1-R6 pass; R7: point-frame regrouping = production table on 252/252 non-zero coefficients, both primes and the second frame; own frozen `F_2 = 1/1200`, `F_3 = 1/25200`, `F_4 = 1/180000` | FR s.2-s.3 |
| E14 | S15 controls | constant-`M_0` split 94/252; dilaton dropped 77/336; leaves removed 88/336; edge sign dropped 103 (worker) / 105 (referee) of 336, orientation-dependent (K5); stored normalisation 21/336 (graphs 1, 2, 5-8); `R(-z)` 112/336; dressing `exp(u z N)` 74/336 | S15 s.6; FR s.3 |
| E15 | S15 rational values | vertex `u^-4(1/25200 - u/9000 + 407u^2/3780000 - 151u^3/3780000 + 226213751 u^4/62052480000000 + ..)`; graph 2 `[u^0] = -202813/38782800000000`; graph 4 `-113/35437500`; graph 9 `127/425250000`; graph 39 `[u^2] = 512/7688671875`; graphs 37, 41, 42 zero through `u^3`; stored shifts `43397/517104000000` (graph 2), `+-43397/37231488000000` (graphs 1, 5) | S15 s.5-s.6; FR s.3 |
| E16 | B-variation | systems R0, S1-S5 (BV s.2) and S6, S7, S8 (BR); genus 2 (6 edge graphs), genus 3 (41): all regular, two primes; genus 4 maxout 2 (10 edge graphs): S1, S2, S3, S5 (2^89-1), S2, S4, S6 (2^61-1); genus 4 full 378/378: S2 (2^89-1, worker + disjoint engine), S7 (2^61-1); vertex `c_g T''^(2-2g)` with `c_2, c_3, c_4 = 1/1200, 1/25200, 1/180000` | BV s.4-s.5; BVT; BR s.2; BR4 s.4.1-s.4.2 |
| E17 | B-variation controls and responses | stored normalisation: genus 4 graphs 1, 2 from `[u^4]`, 3, 4, 7, 8 at `[u^5]`; S2 `m5 = 319129204271/48787200000`; `[u^4]` response `-1/8400000, -1/1008000, 1/900000` (graphs 0, 1, 2; universal); dress0, dress2, staleM, staleHat, staleboth, nosig fail as predicted; non-triviality 6 of 7 (g = 2), 25 of 42 (g = 3) graphs differ from R0 | BV s.5-s.6; BR s.2.1; BR4 s.4.4-s.4.5 |
| E18 | Lemma C numerically | `phi'^5 B_X(phi) = B_q` exactly; per-graph equality with the univ pipeline: S1 (g = 2, 3), S4, S5, S6 (g = 3, 2^61-1), S7 (g = 3, 2^89-1); S8 branch: `P_S8 = -P_q diag(1,-1,-1,1,-1)` exactly (`k <= 7`) | BR s.3; BR4 s.4.3 |
| E19 | model data (Stage 1) | 1395 singular K-correlators, genus 0-4, all zero; 72120 TR correlators | MSG s.5; TRM s.0 |
| E20 | genus 2-5 graph sums | Conjecture N VERIFIED directly at genus 2-5 | SYN2 s.0; MS2 s.1; MS5 s.2 |
