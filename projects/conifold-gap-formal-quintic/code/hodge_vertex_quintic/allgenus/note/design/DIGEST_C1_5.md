# DIGEST_C1 (piece 5 of 6): s.4 Conventions, sign dictionaries, normalisations; s.5 Evidence for the appendix

## s.4 Conventions (as the sources fix them; every place where two sources differ is flagged, with the adjudication)

**CV1. Tower: theory's V^(k) versus the model's Pi_i.** RD1 s.1.3: `V^(k)(q) := (-1)^k [z^(k-1)] <Y^C(-z), q(z)>`, `k >= 1`. MS1 s.1:
`Pi_i(q) := (1/5) [z^i] eta(S(-z) H^2, q(z))` with `Y^C_b(z) = -(1/5) eta(S(z) H^2, H^b)`, i.e. `Pi_i(q) = -[z^i] <Y^C(-z), q(z)>`. Hence
`V^(k) = (-1)^(k+1) Pi_(k-1)` (my reading of the two definitions; the sources state only `V^(1) = PiC = Pi_0`, SYN01 C4, and that the
kernels agree, RD2 s.3.2 "Independent derivation"). Consistent with both string statements: `l_(-1) V^(k) = -V^(k+1)` (RD1 s.2.1) iff
`Pi_0([q/z]_+) = Pi_1(q)` (MS2 s.9(d)). Writers: use one family (the Stage-2 convention is `Pi_i`, `K`, `K_0`) and check the sign once.

**CV2. Names of record** (SYN01 C9; S2A preamble): "(G)" for the gauge identity (not "Theorem G"); "Gap Theorem" for TRM's Theorem G;
`K` (alias `K_vc`) and `K_0`; "ML_K"; "quintic Theorem R". Old names: theory `T_v = K_0`, `Lambda_v = -z + K`, ML_Lambda = ML_K; model's
"Lemma K" = Cor. E' + Lemma T. MS1's (ML1)/(ML2) (pure-H sector / kappa-regularity) are NOT XY's ML1/ML2 (L-regularity / gap).

**CV3. Weights and regularity.** RD1: `w(k,a) = (4k + c_a - 4)/5`, `e = 2 - 2h + Sum w`; MS1: `w(a,k) = (4 - c_a - 4k)/5`,
`e = 2 - 2h - Sum w`; TRM s.3: `e = 2 - 2g + (4d - beta_a)/5`, `beta = (4, 5, 1, 2, 3) = 4 - c_a`. Same exponent. "Regular" = `e >= 0`
(model correlators; "lies in `Q[v]`", S2A s.1); = `Cont_Gamma in Q[[u]]` (flat graph terms, RD2 s.3.2); = `[u^0..u^(2g-3)] G = 0` for a
hatted term `u^(2-2g) G(u)` (UN1 header). "Regular" never means "zero" (TRF s.4).

**CV4. Vertex constant.** Four equal forms (SYN01 C2, `g = 2..12`): `(-5)^(1-g) F_g(x = W)` (tr); `|B_2g|/(2g(2g-2) 5^(g-1))` (brief, NF);
`(-1)^(g-1) B_2g/(2g(2g-2)) 5^(1-g)` (UR, model); `20^(1-g) (-4)^(3g-3) F^Hurw` (LG s.5.4). TRM's "the factor `5^(1-g)` is wrong" concerns
the TR-to-frozen factor for `x = W`, not `c_g` (SYN01 C2).

**CV5. H, the unit, e** (SYN01 C3): `H = d/dv`; `<H>_1 = -1/12`, `<tau_1(1)>_1 = 5/24`, `<e>_1 = +1/8` (tr, model, xyswap, theory). UR s.6.3
and the plan have `+1/12`, `-1/8`: their `H` is `-d/dv` at the frozen point `t = -1 <-> v = -1` (XYRF s.4 confirms `t = -1 <-> v = -1`);
the sign of `e` is free. Adopt `H = d/dv`, `e = tau_1(1) + v tau_0(H)`. TRM s.0 item 2: the project's quoted list `1/12, 1/12, 1/6;
1/600, 1/200` uses `phi = -H`.

**CV6. Lowered versus vector form.** Solutions: lowered components `Y_a = eta(phi_a, J)`; `eta = 5 J` (`J`: `a <-> 3 - a`), `P^vec = J P J`,
`d = (P_0)_33 = (P_0^vec)_00` (TRF s.2). RD1 s.1.1 does not distinguish the two forms ("statements ... are equivalent"); Lemma E and Lemma T
are in vector form, (A) in lowered form. The note must fix one form and state the conversion once (s.7 O-7).

**CV7. The z^5 normalisation** (SYN01 C11): use `P_stored = P_can (1 + m5 z^5 + ..)`, `m5 = -43397/807975000`. Conflicting writings:
RD2 s.3.4 table "`P_11(0,z) = 1 - c_5 z^5 - ..`" (so `c_5 = -m5`); MECH1 s.2 "`c_5 = -43397/807975000` away from canonical" and MS2 s.9
"`P = P_can (1 + c_5 z^5 + ..)`" (so `c_5 = m5`); GI s.4 versus s.5-s.6 (opposite signs; SYN01 C11). MS2 s.6's control uses
`R' -> R' exp(c5 z^5)` with `c5 = 1/7` (a test value).

**CV8. Normalisations of the model curve.** (i) LG model "`x^4 + v/x`, form `dx`" (RD1 s.1.2, MS1 s.0, LG); relative to `H_(0;1,4)`:
"metric x 20 and `z -> -4z`" (MS1 s.1; TRRF correction 6 citing LG s.5.4). (ii) TR: `x = kappa W`, `W = p^4/4 + u/p`, final `kappa = 1`
(`x = W_v`, S2A M1); in TRM `u` denotes the residue, not the conifold coordinate. (iii) `lambda = p^4 + v/p` (`kappa = 4`, `v = 4u`;
arXiv:1509.06954's form; XY's `x = z^4 + v/z`, with `z` the curve coordinate). Frozen `= Lambda^(1-g) kappa^D x` (coefficient of
`prod d xit`), `Lambda = -5 kappa`: `kappa = 1`: `(-5)^(1-g)`; `kappa = 4`, `v = 4`: `omega = (-20)^(g-1) Sum 4^-D <..>^frozen prod d xi`,
`F_g^TR = (-5/16)^(g-1) c_g` (`F_2^TR = -1/3840`, `F_3^TR = 1/258048`); `kappa = -1`: `5^(1-g) (-1)^D`. General residue (`p -> s p`,
`s^5 = u`): `R_TR(z) = M_0(z/s^4) = M_0(zhat)`, `Delta_j = -5 zeta^(-2j) s^2` (TRM s.2.4).

**CV9. Orientation of the TR kernel** (S2B K1, the critic's resolution). Engines (XY conventions): `K(z0, z) = (1/2) Int_(sigma z)^z B(z0, .)/
((y(z) - y(sigma z)) dx(z))`, `omega_(0,1) = y dx`; equal to arXiv:1509.06954 (2-28) (numerator and denominator both flip, TRM s.1) and to
arXiv:2502.18115 (2.1), (1.1) = (3.1) (TRRF s.1). The printed kernel (4-3)/(4-9) of arXiv:math-ph/0702045 is the OPPOSITE (logleg referee,
via S2B K1); it equals arXiv:2212.00320 (1.1) with `omega_(1,0) = -y dx`. So `omega^code_(g,n) = (-1)^n omega^ABDKS_(g,n) =
(-1)^n omega^print_(g,n)` (XY, XYRF s.2; S2B K1), `F_h = -F^print_h`, `Lambda^print = -1` (S2A M9). "The kernel of arXiv:math-ph/0702045"
in TRM s.1 and XYRF s.2 means the common convention, not the printed (4-3). Sign anchors: psi-class intersection numbers on the Airy
curve; `F_2^TR = -1/240` for `x = z + 1/z` (TRRF s.4; XYRF s.2: no tuning; flipping `y -> -y` in kernel and `Phi` leaves `F_g` unchanged).
The critic checked K1 on the formulas quoted in the notes, not against the PDF (S2B K1).

**CV10. Leg dictionaries** (S2B K10). XY: `tau_0(1) <-> -z`, `tau_0(H) <-> -log z`, `tau_k(a) <-> +Int tau_(k-1)(a) dx` (so
`tau_1(1) <-> -(4/5) z^5 + v log z`, `tau_1(H) <-> -(x log z - z^4/4 + v/z)`; XY's log label `tau1H = x log z - z^4/4 + v/z` is minus it,
XYRF s.4). LD, the dictionary of record (S2A step 7.2): `q_(Phi_k^b) = -c_b tau_k(H^b)`, `c = (-1, 1, -1, -3, -2)`,
`Phi_0^b = p, log p, v/p, p^3, p^2`, `D Phi_(k+1) = -Phi_k` (`D = d/dx`); the TR leg of `v^m tau_k(H^b)` is `-c_b^-1 v^m Phi_k^b`;
`tau_0(1) <-> y = p`; dilaton leg `Int y dx = -Phi_1^0`. Both pass the string and dilaton checks; only LD is dual to the `xif_(a,d)`
(TRM s.2.2: `Sum Res y d xif_(a,d) = delta_((a,d),(0,0))`, `Sum Res Phi d xif_(a,d) = -delta_((a,d),(0,1))`). No XY result depends on the
difference. String sign in XY's dictionary: the unit acts by `+d/dx` (XYRF s.4; SYN01 C8). LD's `c_b` clash with the charges `c_a` (s.7).

**CV11. Euler-dilaton edge leg**: `psi A1 <-> -(16/45) z^9 + (v/5) z^4`, not `-(2/5) z^9 + ..` (XYRF s.4; SYN01 C7; `cr_conv` (5)).

**CV12. TR metric and psi signs.** arXiv:1509.06954 (2-7) is MINUS the metric (1.6) of arXiv:2211.12259 and of LG s.5.1 (TRM s.1);
correct combined statement: `TR(lambda) = (-1)^(g-1+D) x` (Hurwitz CohFT in the LG convention), the two signs cancel for `F_g`
(`D = 3g - 3`) (TRRF correction 6). "`x = -W` flips `psi -> -psi`" holds only up to `(-1)^(3g-3+n)` (TRRF s.3.1, correction 2).

**CV13. MRF's flat normalisation**: "`Pi_0(-z R^-1 1) = -t_c` (normalisation `Pi_0(-z) = -v`)" is the sources' `+t_c` with the overall
sign of `Pi_0` flipped (MS1, RD1: `Pi_0(-z phi_0) = +v`). Same statement.

**CV14. Three leaf differences** (not a conflict): `s = d T_P - v = u^2/30 + (127/3150) u^3 + ..` (frame independent; TRF s.2-s.3);
`t_c - vhat = s/d = u^2/30 + O(u^3)` with `vhat = v/d` (RD2 s.3.3), `d = tau^(1/4) xhat` (RD2 s.3.2) `= tau^(1/4) varpi_0` (TRF s.2);
`t_c - v = (13/90) u^2 + ..` (MRF item 8(b), no factor `d`).

**CV15. lambda and c0** (SYN01 C10; S2B K6): RD2 s.3.3 `1/lambda = P_11(0, z)`; UN1 s.1 `c0 = 1/lambda`, `lambda = f(0)`; URF s.2.1
`c0 = 1/(f(0) |g(0)|)`; BVRF `|c0| = m^10/|f(0) g(0)|` (general; sign by `sigma > 0`); define `lambda` by `r = -1/lambda^2`. The univ
runs have `g(0) = -1`, where all agree.

**CV16. "Frame" has four meanings**: `frame(k, a)` in `K_0` (D9); the conifold frames `(alpha, beta)`, labels `(0, 0)` and `(1/3, -2/7)`
(normalisation `y_0 = f_0 + alpha g + beta f_2` of the conifold period; MECH5 s.3.7 records harmlessly mislabelled frames
`y_0 = f_0 + a f_1 + b f_2` in the Stage-5 critic's checks); the flat versus hatted frame (variables); UR's frame `K = (1, uH, H^2, H^3, H^4)`.

**CV17. Genus one**: `Fb_1(lambda q) = Fb_1(q) - (5/24) log lambda`; `<tau_1(1)>_1 = 5/24 = rank/24` (XYRF s.4).

**CV18. Fields and primes**: `2^89 - 1`, `2^61 - 1` (both `= 1 mod 5`; MRF s.2), `p3 = 1208925819614629174706261` (least prime `= 1 mod 10`
above `2^80`) and a 128-bit prime (TRF preamble), exact `Q`, `Q(zeta_5)`. Outer-graph indices of univ are positions in the project's
genus-g graph list (not Mathematica indices); genus-4 graphs with at most two edges: 0 vertex (4); 1 bridge (1,3); 2 bridge (2,2); 3 loop
on (3); 4 (0,3) loop at the genus-0 vertex + bridge; 5, 6 (1,1,2) chains; 7 (1,2) loop at the genus-1 vertex + bridge; 8 (1,2) double
edge; 9 (1,2) bridge + loop at the genus-2 vertex; 10 two loops on (2) (UN1 header, s.2.2). Appendix only.

## s.5 Evidence for the appendix (finite checks; none may be cited by a proof)

**Transport.**
* EV-T1 Lemma E: exact Q, random symplectic `R = exp(X)` (generic and Z/5-graded, `z^8..z^10`), random `f` (3-6 terms), both slots,
  6 cases; `E_R` divisible by `psi + psi1` iff (Sp); non-symplectic and wrong-sign controls fail. TRF s.3 (`rf_lemE`).
* EV-T2 flat Cor. E' and Lemma T on `P_can`: frames (0,0), (1/3,-2/7), `z^0..z^16`, `u < 16`, exact Q: (Sp) exact; (G) on 575 complete
  coefficients; `T_P = t_c` through `u^14`; `V^(1)(q_T) = t_c`, `V^(2..8)(q_T) = 0`; `d (P_0)_00 = 1`; the series `s`; edge divisible by
  `psi + psi1`, contraction with `V^(1..8)` zero on 2985 coefficients; control without `(-1)^n`: 92 slots non-zero; stored: (G) fails (77),
  `V^(6)(q_T) != 0`, edge contraction fails (92), symplectic defect at `z^10`. TRF s.3 (`rf_flat`).
* EV-T3 model Lemma K on the actual hatted rows `R'` (`M = M_0(tau zhat)`): referee rows `zhat <= 24`, `u < 18`, `2^89-1` and `2^61-1`,
  frames (0,0), (1/3,-2/7): shifts `k = 0..20` vanish (21-23 vacuous), 2160 coefficients per `k`, 20790 non-trivial sums per run,
  `Pi_k(T) = 0` for `k = 1..20` mod `u^18`; production rows `zhat <= 6`, `u < 6`: `k = 0..4` (150 per `k`, 500 sums), leaves `k = 1..2`;
  control `c5 = 1/7`: every shift fails (750 of 2160 at `k = 0`, 180 at `k = 10`), the leaf acquires exactly a `Pi_5` component. MS2 s.6;
  MRF item 7 (re-run identical).
* EV-T4 flat Lemma K (referee): `P` to `z^22`, `u^9`, exact Q, `P_11(0,z) = 1`: symplectic through `z^22`; 780 positive-power coefficients
  of `P Y^C` vanish; `Pi_m` of the first edge slot zero on 4515 coefficients (`m <= 6`; 742 sums of `>= 2` terms; a random covector is
  non-zero on 868); leaf `Pi_1..Pi_6 = 0`, `Pi_0 = t_c` (CV13) through `u^8`; controls: stored (not symplectic from `z^10`, (G) fails on 102,
  edge statement 627 non-zero, leaf `Pi_0` moves at `u^5`); `P exp(z^5/7)` (leaf `Pi_5 = (1/7) t_c`). MRF item 6 (`rf_lemmaK`).
* EV-T5 (G) itself: exact over Q, `z <= 20`, `u < 29`, two frames, canonical and stored; residue relations `n <= 8`. SYN01 s.1 E1, E3 (GI).

**Theorem R.**
* EV-R1 genus 2 in the flat frame (never computed before TRF): frame (0,0): vertex `u^-2 (1/1200 - 7u/6000 + 287 u^2/864000 + ..)`; edge
  graphs `[u^0] = [u^1] = 0`, `[u^2] = 479/756000` (loop on g = 1), `-127/2016000` (bridge g1-g1), `1/2250` (two loops on g = 0), `0, 0, 0`,
  equal to the hatted values on record (UR s.6.2); frame (1/3,-2/7): edge graphs regular, vertex `(1/1200) t_c^-2` + regular; controls X1
  (leaves removed: `-1/900 u^-1`), X3/X4 (`(1 + u^m z E_21) P`, `m = 0, 1`: singular from `[u^m]`); non-gauge `exp(u^m z N) P`,
  `N = E_04 + E_43`: `m = 0` regular with `T_P = u + u^2/5 + ..`, `m = 1` regular with `T_P = t_c + O(u^3)`. Model correlators from MVertex6 with
  the referee's own LG `M_0` (128-bit prime, rationally reconstructed, re-checked mod `2^89-1`). TRF s.3 (`rf_g2flat`).
* EV-R2 (cross-reference, other readers): genus 2-5 by hatted graph sums (MECH5 s.2; NF); genus 3 flat = hatted, 42 graphs x `u^-4..u^3`,
  41 edge graphs regular, vertex `(1/25200) t_c^-4`, 108 + 4 non-zero regular coefficients (S2B s.3 (d), K4).

**Model Lemma and singular structure.**
* EV-M1 ML_Lambda lists (theory): genus 0/1/2/3 with `<= 5/5/4/3` insertions, frames `k <= 3h + 3`, psi-deficit `<= 6`: 43/418/425/127 singular
  lists, all vanish (`2^89-1`, `2^61-1`; referee at `p3` with own frames at general `v` and own `M_0`); controls C1 5/64/165/88 non-zero,
  C2 269/365, 1264/1429, 740/878, 172/212 non-zero. Refined (TRF s.3a): genuine cancellations 5/65/168/88 (326), every flat term zero although
  allowed 22/328/253/39 (642), trivial 16/25/4/0 (45); singular lists skipped 60 (g = 2), 164 (g = 3) (RD1's "102", "262" include regular
  lists); C3 non-zero 4/10, 46/212, 123/274, 76/95; C5 4/9, 50/170, 120/209, 62/69. Deep part (frames using the `z^-5` terms): 15 lists
  (1, 3, 6, 5 at genus 1-4) vanish at `p3`, C4 non-zero on all; e.g. `<frame(6,0) H^5>_1 = -5/24 + 7/24 - 1/12 = 0` and `<frame(6,0)>_4 = 0`.
  RD1 s.2.4; TRF s.0 item 3, s.3a.
* EV-M2 ML on K (model engine, exact over Q): regressions 2512 + 3955 + 6125 + 23125 flat correlators = MVertex6 mod `2^89-1`; 3987
  correlators with a `kappa` at `(0,3..6), (1,1..5), (2,1..3), (3,1..3), (4,1)`, 1019 singular, all 0; ML1 in 16 cases; singular per type
  (0,4) 6, (0,5) 37, (0,6) 152, (1,2) 3, (1,3) 19, (1,4) 83, (1,5) 313, (2,1) 2, (2,2) 17, (2,3) 90, (3,1) 4, (3,2) 41, (3,3) 246, (4,1) 6;
  prediction `(2,4)`: 1302 correlators, 376 singular, all 0, `<H^4>_2 = 1/10`; totals 5289 / 1395 singular / 0 violations; basis insertions:
  1562 singular, the 593 inside K with remainder 0, off K mostly non-zero. Referee: same counts at both primes; controls `kappa^-` 476/1395,
  bare `tau_k(X)` 456 non-zero; new `(0,7)` 191 singular, `(1,6)` 376 singular, 0 non-zero (controls 51, 131); 986 singular-list entries
  agree. `amp(4,2)` not computed. MS1 s.2, s.5; MRF items 1, 4.
* EV-M3 one-point and two-point tables, Hessians, H^2-tower data (181 pure `kappa(H^2)` correlators, `e >= 1`), `Delta_g(5)` and the
  per-graph responses: MS1 s.3-s.4; MS2 s.9; MRF items 2, 3, 5, 10; TRF s.3a (`rf_deep`); SYN01 C1.
* EV-M4 (cross-reference) Stage-2 falsifier (b): TR `omega_(4,2)`, `omega_(3,3)` with every L-leg and `-log p`: 17416 singular tuples per prime,
  all 0; referee to `g + n <= 7`, `g <= 6`: 40739 multisets, 0 non-zero; pure-log values = Gap derivatives. S2B s.3 (b).

**Universality and the failing controls** (hatted engine, `2^89-1` unless stated; gauges exact over Q).
* EV-U1 genus 4, canonical versus stored (E4/U4): inputs (every `p_(k,b)`, `k <= 12`, a power series; (G) on `z^-6..z^6` mod `u^6`; rows
  symplectic mod `(zhat^10, u^6)`); case a (quintic, canonical): 378/378 edge graphs regular, all 379 rows equal the natural run of record,
  vertex `(1/180000) t_c^-6 + 0 + .. + 0`; c1, c2 (perturbed `F`, `Gm`; c2 has `lambda = 3/2`): 378/378 regular, vertex
  `(1/180000) (T/lambda)^-6 + 0`; b (stored): exactly graphs 1, 2 (from `[u^4]`) and 3, 4, 7, 8 (`[u^5]`) singular, the other 372 rows equal
  case a, total unchanged; d (c1 stored, maxout 2): 6 of 10 singular; responses per unit `m5` (vertex `-1/8400000, -187/680400000`; graph 1
  `-1/1008000, 13/6804000`; graph 2 `1/900000, -53/32400000`; graphs 3, 4, 7, 8 at `[u^5]`: `-41/68040000, -11/34020000, 49/16200000,
  -17/8100000`), sums `0, 0`; the `[u^4]` column is universal (quintic = c1). Referee: two new systems r1 (`f(0) = 2`, `g(0) = -3/2`), r2
  (frame (1/3,-2/7)) full, 378 regular; all responses at both primes; module-disjoint engine reproduces c1 (379 rows) and case a (0 of 379
  differ). UN1 s.2; URF s.2; MECH1 s.5.
* EV-U2 genus 2-3 universality: six perturbed systems and two quintic frames in canonical normalisation, all edge graphs regular, vertex in
  the system's own `T`. UN1 s.5 (B); SYN01 E12.
* EV-U3 annihilating non-gauge dressings: A1-A4 (genus 2: 6 regular; genus 3: 41 regular; both primes), `delta_0` non-constant,
  `T_eff = -[z^-1 e_0] D Y^C_Q` (A1: `u + (49/45) u^2 + ..` versus `t_c = u + (7/10) u^2 + ..`); genus 4: seed 20261113 (27 entries,
  maxout 2 and full) and seed 2 (29 entries, full): 378/378 regular, vertex `(1/180000) T_eff^-6 + 0`; referee: Cayley dressings, A1 full at
  genus 4 at both primes; pure-model dressings MA1, MA2 (genus 2-3; referee mA7 at genus 4, E <= 2). UN2 s.3.2-s.3.5; URF s.3, s.5.
* EV-U4 controls: B1 (genus 2: graph 2 at `[u^1]`; genus 3: 14 singular), B2 (genus 3: graphs 1, 2 at `[u^2]` = `1/252000`, `1/14000`, graph 3
  at `[u^3]`), C1a (3; 28 singular), C1b (1; 8), C1c (regular at genus 2; genus 3: graphs 2, 4, 8, 12 at `[u^3]`), MB (pure model, as B2),
  xA1 (wrong covector; genus 3: two graphs at `[u^3]`, `91/810000`, `-91/6480000`), R1 (random column-1 orbits: 16 singular). UN2 s.3.2,
  s.3.5; URF s.3.1, s.5.
* EV-U5 first order (genus 2: 9 variables, rank Ann 3, rank Resp 3, rank(Resp on ker Ann) 0, ker Resp = ker Ann; genus 3: 30 variables,
  rank Ann 10, rank Resp 8, rank(Resp on ker Ann) 0, rank(Ann on ker Resp) 2; without flat regularity: 20 resp. 68 variables, rank(Resp
  on ker Ann) 5 resp. 23); E5: the surplus `w_1`, `w_2` keep the genus-3 edge graphs regular but break the vertex law (`T_eff^-1`:
  `11/136080`, `-11/54432`; `5 w_1 + 2 w_2` keeps it) and are singular at genus 4 (all 10, resp. 9 of 10, graphs with `<= 2` edges, from
  `[u^3]`). UN2 s.3.1, s.4; URF s.3.2-s.3.3 (inclusion confirmed, ranks not re-derived).
* EV-U6 (cross-reference) B-universality (Stage-2 falsifier (c)): S1-S5 at genus 2, 3 (both primes), genus 4 maxout 2 for S1, S2, S3, S5 and
  S2 in full (378/378); referee S4, S6, S7 (S7 full at `2^61-1`); controls stored, dress0, dress2, staleM, staleHat, staleboth, nosig fail.
  S2B s.3 (c).

**TR and the x-y swap on the model.**
* EV-TR1 TR versus the frozen theory: `TR x (-5)^(1-g) = MVertex6(M_0)` on 72120 colour-basis correlators (psi-deficit `<= 6`; `(0,3)`, `(0,4)`,
  `(0,5)`, `(1,1)`, `(1,2)`, `(1,3)`, `(2,1)`, `(2,2)`, `(3,1)`; 4665 non-zero) at both primes; referee: `kappa = 1, 4, -1`, 6495 correlators each
  (1660 non-zero), 0 mismatches; `R_TR = M_0(z/kappa)` exactly through `z^8` (worker), `z^10` (referee), all 25 entries over `Q(zeta_5)`;
  dilaton leaf through `z^8`; exact flat values (genus-2 one-point `<t1(1)> = 1/600`, `<H> = -1/600`, `<t2(H^4)> = 1/1200`,
  `<t3(H^3)> = 7/5760`, `<t4(H^2)> = -5/1152`; genus 3 `<t5(H)> = 53/1209600`, `<t6(1)> = -1/46080`, `<t7(H^4)> = 5/82944`, ..). TRM s.0,
  s.2.3, s.3; TRRF s.2, s.4.
* EV-TR2 gap: `F_2..F_5 = -1/240, 1/1008, -1/1440, 1/1056` on `x = W` (exact over `Q(zeta_5)` for `g = 2, 3`; both primes `g <= 5`) =
  `c_2..c_5`; formula (1.3) evaluated for several polynomial parts, `g = 2..6`; `A_g = -B_2g/(2g)` for `g = 2..12` (referee `2..14`); TR at
  generic points of `H_(0;1,4)` and on `x = p^r/r + 1/p`, `r = 1..6`; referee exact over Q at rational branch points, `r = 1..5`, pole moved,
  and a curve with TWO simple poles (`F_2`, `F_3` = the sum of the two pole terms; one-pole value rejected). TRM s.4.2; TRRF s.4.
* EV-TR3 `v -> 0`: 13 `(g, n)`, both primes (referee 8), limits such as `omega_(1,1) -> -(5/8) dp/p^6`, `omega_(2,1) -> (2079/128) dp/p^16`.
  TRM s.5; TRRF claim 11.
* EV-XY1 closed formulas = own TR (exact `Q(zeta_5)`: `(0,3), (1,1), (0,4), (1,2), (2,1)`; mod `2^89-1`: `(3,1), (2,2), (1,3)`; referee:
  exact Q at generic rational points); `F_2 = -1/3840`, `F_3 = 1/258048` at `v = 4`; generic points of the family, both primes; the
  vanishing period for `g = 1..10`; referee: (1.3) of arXiv:2502.18115 for `g = 2..6` at five points, double-pole control gives `F_2 = 0`.
  XY s.2; XYRF s.2.
* EV-XY2 legs: 28 non-zero L-correlators at `(1,1), (2,1), (1,2), (2,2)`, all with exponent `>= 0`; bad legs singular (`<omega_(g,1), -log z>`
  `= -(1/12) v^-1, (1/120) v^-3`; `<omega_(2,1), tau_1(1)> = 2 F_2`; `<omega_(1,1), z^-5> = 15 v^-2`; `<omega_(2,2), -log z, -log z> = F_2''`;
  mixed `<omega_(1,2), -(4/5) z^5, z^-5> = 60 v^-2`); referee run B off the line (12 points, 36 L-correlators polynomial in `S`, value at
  `S = 0` = TR of the limiting curve). The "11 forced zeros" are not evidence. XY s.3; XYRF s.3.
* EV-O1 critic checks: `cr_conv` 29 PASS ((1) the four vertex forms `g = 2..12`; (2) covector signs `m <= 60`; (3) per-graph `c_5` responses
  sum to 0; (4) frame independence `1/15`; (5) the Euler-dilaton leg `-16/45`); `cr_dictionary` to weight 240 / 320 / 400 (5880 / 10400 /
  16200 log-free vectors, all in L); `cr2_e2e`: 52 frame one-point functions at `g = 2, 3` (`<frame(3,3)>_2 = 3/3200`, `<frame(4,2)>_2 =
  -41/9600`, `<frame(6,0)>_3 = 3/128000`, `<frame(7,4)>_3 = 223/2688000`). SYN01 preamble, s.3; MECH1 s.5; MECH2 s.4; S2A step 7.2 box.
