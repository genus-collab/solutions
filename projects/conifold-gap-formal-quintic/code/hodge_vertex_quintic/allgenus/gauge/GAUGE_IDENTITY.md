# The gauge identity (G): the regular gauge preserves the vanishing period

Date 2026-11-13. Folder `allgenus/gauge/` (task key `gauge`). AI-derived and AI-checked only; no human has read it. Status labels: PROVED (complete
argument written here, or a quoted theorem with its hypotheses checked), VERIFIED (exact computation; orders, field/primes and frames stated),
CONJECTURAL, REFUTED. Conventions: `natural_form/theory/UNIFORM_REDUCTION.md` s.1-s.5 (flat frame, index `a` in Z/5, `th = z d/du`, `c = (0,-1,3,2,1)`,
`eta_ab = 5 delta_{a+b = 3 mod 5}`, `E_ab` matrix units, conifold frame `(alpha, beta)`), `natural_form/lg/LG_MODEL.md` (the model), hatted
conventions of `OU_EXTENSION_genus3.md` s.1. Primes p89 = 2^89-1, p61 = 2^61-1. The hatted translation (covector, (L), (E), (T)) is in
`GAUGE_IDENTITY_HAT.md`. Every log quoted below was produced (or re-produced) by this task today.

## 0. Summary

1. **Theorem G (PROVED given Theorem U of UNIFORM_REDUCTION s.3).** For every regular gauge `P(u, z)` from the pulled-back model to the
   quintic system, `P(u, z) Y^C(v(u), z) = P_11(0, z) Y^C_Q(u, z)` in `Q((z))[[u]]^5`. Hence the canonical gauge (`P_11(0, z) = 1`) maps the
   vanishing solution of the model to that of the quintic, with no positive power of z: **(G)**. The proof needs no monodromy: a solution
   regular at `u = 0` is determined by its value at `u = 0` (Lemma 1), and both sides have the value `-P_11(0, z) e_1` by (R1). VERIFIED
   (inherited from the first attempt, re-run): 176 positive-z coefficients vanish (152 genuine cancellations, up to 142 terms) and 1038
   coefficients equal `Y^C_Q`, `z^-38..z^20`, `u < 29`, frames (0,0) and (1/3,-2/7), exact over Q (`gi_check_*.log`).
2. **Symplecticity (PROVED given Theorem U).** For every regular gauge `P(-z)^T eta P(z) = c(z) c(-z) eta`, `c(z) = P_11(0, z)`. The canonical
   gauge is symplectic; the symplectic regular gauges are `P^can c(z)` with `c(z) c(-z) = 1` (graded: `c = +-exp(odd series in z^5)`), so
   symplecticity fixes only the EVEN part of the normalisation and cannot see `c_5`; only (G) singles out the canonical gauge. VERIFIED:
   canonical gauge symplectic `z^0..z^20` mod `u^29`; the stored `ur_gauge` normalisation has EXACTLY the defect `(c(z)c(-z) - 1) eta^-1`
   (first at `z^10`, constant in u), both frames (`gi_relations_*.log`).
3. **(G*) and the diagnosis of the first attempt's failure.** `(G*)`: `Y^C_Q(-z)^T eta P(z) = Y^C(-z)^T eta` follows from (G) + 2 (PROVED). The
   reported failure was **a floating-point bug in the test, not truncation and not a wrong statement**: the covector map used `(-1) ** m` with
   `m < 0`, which in Python is the float `-1.0`; 73 of the 119 tested slots then failed by rounding residues `5.1e-21..8.9e-16` (the old log
   printed only the first four, which are the four slots quoted in the brief). With an integer sign (G*) holds on all 119 (old bound) resp.
   165 (sharp bound) complete coefficients `u^0..u^8`, both frames, exact (`gi_diag_old.log`; `gi_relations_*.log`, section DIAGNOSTIC; s.8).
4. **Residue relations to all orders (PROVED).** Expanding (G) and (G*) in u expresses the H-column and the H^2-row of every `[u^n] P` through
   lower orders. Order 0: (R1); order 1: (R2), (R3a), (R3b); order 2: (R3c); the u^3, u^4 relations are written out in s.5, the u^4 one
   containing `(1 - P_22(0, z))/(24 z^4)`, the first place where the normalisation enters. VERIFIED exactly, `n <= 8`, `z <= 20`, two frames.
5. **The actual R-matrix (VERIFIED; all orders CONJECTURAL).** From the referee's from-scratch rows (QRR normalisation, no gauge input) the
   hatted gauge `Ghat^act = Qhat M_0(tau zhat)^-1 Psi0^-1` equals the canonical gauge computed from the periods in ALL 1058 slots
   `zhat^k u^N` (`k <= 24`, `N < 18`, about 995 of them non-zero), at p89 and p61, in frames (0,0) and (1/3,-2/7). So `c^act(z) := P^act_11(0, z)
   = 1 + O(z^25)`. Since both gauges are symplectic, `c^act = exp(a_1 z^5 + a_3 z^15 + ...)`: the `z^10`, `z^20` coefficients are forced, and the
   independent tests are `z^5` (the constant `c_5 = 43397/807975000` of UNIFORM_REDUCTION s.4) and **`z^15`** (read off rule W before, SYNTHESIS
   item 8; here tested directly, with every other slot of the window) (`gi_actual_norm_*.log`). The brief's `z^10` test is not independent.
   This covers every order of R' used by a graph sum of genus `g <= 9` (`zhat^(3g-3)`, `u^(2g-3)`): the production `R'` of genus 2-5 IS the
   canonical gauge at the orders used.
6. **Hatted translation (GAUGE_IDENTITY_HAT.md).** The vanishing covector in closed form `phi_{i,m} = (u tau^(-1/2)/5) tau^(-m) zeta^(i(m-2)) /
   (4^q q! (m-q)!)`, `q = floor(m/5)`; the contraction Lemma (L) (PROVED from symplecticity alone); (E) every edge contracts phi to zero and
   (T) the leaves give `Pi(T') = u tau^(-3/2)/xhat - t_c` (both PROVED for the canonical gauge given Theorem U, all genera, as identities of
   u-series). VERIFIED on the actual rows (2160 edge coefficients, `zhat <= 24`, `u < 18`, two primes, two frames).
7. **Not established.** (H) and `c^act = 1` to all orders (a connection problem between `q = 0` and the conifold, s.6); and how (E), (T) turn into
   regularity of the graph sum: that needs the model's singular structure beyond the vanishing covector (HAT s.5). Theorem U is AI-audited only.

## 1. Setting

Quintic: `th Y = A Y`, `A = Sum_a b_a E_{a,a+1}`, `b = (t_c', -1/(u varpi_0^2 (1-u)^3 t_c'^2), t_c', -varpi_0/(1-u), -varpi_0/(1-u))` (frame
`(alpha, beta)`: `varpi_0 = varpi_0^cf`, `t_c = u + O(u^2)`), so `A = u^-1 r E_12 + A_0 + u A_1 + ...`, `r = -1`, `A_j` rational constant matrices.
Model: `z d/dv Y_b = bmod_b Y_{b+1}`, `bmod = (1, -1/v, 1, -1, -1)`, pulled back by `v(u) = u tau^(-5/4)` (`1/tau = 2F1(1,4/5;9/5;u)`):
`th Yt = At Yt`, `At = v'(u) Amod(v(u)) = u^-1 rt E_12 + ...`, `rt = -1`. Both have `b_0 = b_2`, `b_3 = b_4`, which is exactly `eta A = A^T eta`
(entry `(c, d)` of `eta A` is `5 b_{d-1} delta_{c+d = 4}`). Rings: `Q[[u]]`, `L = Q((z))` (finitely many negative powers), `L[[u]]`.
**Theorem U** (quoted, UNIFORM_REDUCTION s.3, PROVED there, AI-audited): there is `P = Sum_k z^k P_k(u)`, `P_k` in `gl_5(Q[[u]])`, `P_0` diagonal
with unit entries, `z P' = A P - P At`; it can be taken graded, `P_k = Sum_b p_{k,b} E_{b+k,b}`; every such `P` with `P_0` invertible is `P c(z)`,
`c` in `Q[[z]]^x` (graded: `c(z^5)`). `P` lies in `gl_5(Q[[z]][[u]])` inside `gl_5(L[[u]])`. Normalisations: **stored** (`ur_gauge`: `p_{5q,0}(0) = 0`,
`q >= 1`) and **canonical** (`P_11(0, z) = 1`). The stored gauge has `P^s_11(0,z) = 1 - (43397/807975000) z^5 - (765368951400317/1305647201250000000) z^10
- ...` in both frames.

## 2. The two vanishing solutions

**Lemma V (model; PROVED).** `I^C_m(v, z) := Sum_{i >= 0} (-1/(4z))^i/i! (-v/z)^(m+4i+1)/(m+4i+1)!` (terms with `m+4i+1 >= 0`) is the formal residue
`Res_{x=0} x^m exp(-(x^4/4 + v/x)/z) dx`. (i) `z d/dv I^C_m = -I^C_{m-1}` (differentiate the residue); (ii) `m z I^C_{m-1} - I^C_{m+3} + v I^C_{m-2} = 0`
(the residue of the exact form `d(x^m e^(lambda/z))`). The three-line argument of UNIFORM_REDUCTION s.1 uses only (i), (ii), so
`Y^C := (I_0, -I_{-1}, -I_3, I_2, I_1)^C` solves the model system. `Y^C(0, z) = -e_1` (only `I_{-1}` has a `v^0` term). `Y^C_b` is a sum of monomials
`v^j z^-n` with `n = 5i + c_b + 1`, `5j = 4n + c_b + 1`: at order `u^N` of `Y^C(v(u), z)` only `n <= (5N - c_b - 1)/4`, so `Y^C(v(u), z)` is in
`Q[z^-1][[u]]^5`. Leading terms `Y^C = (-v/z, -1, -v^4/(24 z^4), -v^3/(6 z^3), v^2/(2 z^2)) + ...`.

**Lemma Q (quintic; PROVED).** Put `y^(0) = -e_1`, `y^(n+1)_a = Int_0^u b_a y^(n)_{a+1} du`. (a) `y^(n)_a = 0` unless `n = 1 - a mod 5`, and then
`ord_u y^(n)_a = (4n + c_a + 1)/5` EXACTLY, with the same leading coefficient as `[z^-n] Y^C_a(v(u))`. (b) Every integrand is a power series:
the only pole is that of `b_1`, which multiplies `y^(n)_2` of order `(4n + 4)/5 >= 1`; so no logarithm arises. (c) `Y^C_Q := Sum_n z^-n y^(n)` lies in
`Q[z^-1][[u]]^5` (order `>= 4n/5` at `z^-n`), solves the quintic system, and `Y^C_Q(0) = -e_1`. (d) `[z^-1] Y^C_Q = -t_c e_0`.
*Proof.* (a) by induction: grading, since `y^(n+1)_a` is built from `y^(n)_{a+1}`; valuation: `ord y^(n+1)_a = ord b_a + ord y^(n)_{a+1} + 1 =
(4n + c_{a+1} + 1)/5 + 1 - delta_{a,1} = (4(n+1) + c_a + 1)/5`, using `c_{a+1} = c_a - 1 + 5 delta_{a,1}`; integration of `u^k`, `k >= 0`, is exact on
leading terms, and `b_a = bmod_a(u)(1 + O(u))`, `bt_a = bmod_a(u)(1 + O(u))`, while `Y^C(v(u), z)` is the iterated-integral solution of the
pulled-back system (Lemma 1 below: same value `-e_1`), so the leading coefficients obey the same recursion. (b), (c)
follow; the equation is `y^(n+1)' = (b_a y^(n)_{a+1})_a`, i.e. `z dY/du = A Y` coefficientwise. (d) `y^(1)_0 = -Int_0 t_c' = -t_c`. QED
VERIFIED (`gi_formulas.log` (V)): exact valuations and equal leading coefficients for all graded slots `n <= 45`, both frames.

**Lemma 1 (uniqueness of regular solutions; PROVED).** Let `B = u^-1 B_{-1} + Sum_{j>=0} u^j B_j`, `B_j` constant matrices (no z). If
`Y` in `L[[u]]^5` satisfies `z dY/du = B Y` and `Y(0) = 0`, then `Y = 0`. *Proof.* Order `u^n`: `((n+1) z - B_{-1}) Y^(n+1) = Sum_{j=0}^n B_j Y^(n-j)`.
`det((n+1) z - B_{-1})` is a non-zero polynomial in z, hence a unit of `L = Q((z))`; induction. QED
*Remark (monodromy).* Over `L` the regular solutions are parametrised by `Y(0)` with `Y(0)_2 = 0` (order `u^-1`); the fifth solution is
`log(u) Y^C_Q + W` with `W` regular, `W_2(0) = z`, so the logarithm of the monodromy maps it to `Y^C_Q` and kills the regular ones: `Y^C_Q`
spans the image of the log-monodromy (the vanishing cycle), and the same for `Y^C`. This is the picture of the plan; the proof below needs
only Lemma 1.

## 3. Theorem G

**(R1) (PROVED, UNIFORM_REDUCTION s.5).** Order `u^-1` of `z P' = A P - P At`: `r E_12 P(0) = rt P(0) E_12`, i.e. row 2 of `P(0, z)` is
`P_22(0, z) e_2^T`, column 1 is `P_11(0, z) e_1`, and `r P_22(0, z) = rt P_11(0, z)`; here `r = rt`, so `P_22(0, z) = P_11(0, z)`.

**Theorem G (PROVED given Theorem U).** For every regular gauge `P`: `P(u, z) Y^C(v(u), z) = P_11(0, z) Y^C_Q(u, z)` in `L[[u]]^5`.
*Proof.* `Yt := Y^C(v(u), z)` solves `th Yt = At Yt` (chain rule, `v` a power series with `v(0) = 0`); `P Yt` is in `L[[u]]^5` (products of
`Q[[z]][[u]]` and `Q[z^-1][[u]]`) and `th(P Yt) = (A P - P At) Yt + P At Yt = A (P Yt)`. Its value at `u = 0` is `P(0, z)(-e_1) = -P_11(0, z) e_1`
by (R1). `P_11(0, z) Y^C_Q` is a solution (the factor is constant in u) with the same value. Lemma 1 with `B = A`. QED
**Corollaries.** (a) Graded `P`: `P_11(0, z) = Sum_q p_{5q,1}(0) z^(5q)` is a unit of `Q[[z^5]]` (`p_{0,1}(0) != 0`), so `P^can = P P_11(0, z)^-1` is a
graded regular gauge (inside the freedom `c(z^5)` of Theorem U), and it is the unique regular gauge with `P Y^C = Y^C_Q`. (b) For any regular
gauge, `P Y^C` has no positive power of z iff `P_11(0, z)` is constant (the `u^0` term of `P Y^C` is `-P_11(0, z) e_1`). (c) The statement
uses nothing but the shape of the two systems: for any pair of cyclic systems with one simple pole and equal products, the regular gauge
maps the vanishing solution of the second (iterated integrals, value `-e_1`) to `P_11(0, z)` times that of the first, whose `[z^-1]` term is
`-Int_0^u b_0 e_0` (the system's own vanishing coordinate in place of `t_c`): the universality seen in `univ/`.

## 4. Symplecticity and (G*)

**Proposition S (PROVED given Theorem U).** For every regular gauge: `P(u,-z)^T eta P(u,z) = c(z) c(-z) eta`, `c(z) = P_11(0, z)`.
*Proof.* `S := P(-z)^T eta P(z)`. From `z P'(z) = A P(z) - P(z) At` and `z P'(-z) = -A P(-z) + P(-z) At` (replace z by -z), `z S' = -P(-z)^T A^T eta P(z)
+ At^T S + P(-z)^T eta A P(z) - S At = At^T S - S At` since `A^T eta = eta A`. With `S = eta T` and `At^T eta = eta At`: `z T' = At T - T At`, so `T` is a
regular gauge from the model to itself (coefficients in `Q[[u]]`, `T_0 = eta^-1 P_0^T eta P_0` invertible). By the uniqueness clause of
Theorem U applied to the pair `(At, At)`, whose gauges include the identity, `T = t(z) 1`. Entry `(2, 1)` of `S` at `u = 0`: by (R1) only `a' = 1`
survives in `Sum_{a,a'} P_{a2}(0,-z) eta_{aa'} P_{a'1}(0,z)`, giving `5 P_22(0,-z) P_11(0,z)`; `eta_21 = 5`; so `t(z) = P_22(0,-z) P_11(0,z) = c(-z) c(z)`. QED
**Consequences.** Canonical implies symplectic. Symplectic regular gauges: `P^can c(z)` with `c(z) c(-z) = 1`; for graded gauges
`c = +-exp(a_1 z^5 + a_3 z^15 + a_5 z^25 + ...)`. So symplecticity leaves every ODD multiple of `z^5` free: it does not determine `c_5`.
The stored gauge: defect `(c(z) c(-z) - 1) eta^-1 = ((2 c_10 - c_5^2) z^10 + ...) eta^-1`, with `2 c_10 - c_5^2 = -18223115111903/15543419062500000`
(`c = P^s_11(0, z)`). The first attempt's control expected the stored gauge to fail at `z^5`; that expectation was wrong (`c(z)c(-z)` is even);
the corrected check compares the whole defect with this prediction. VERIFIED exactly, `z^0..z^20`, mod `u^29`, both frames.
**Corollary (G*) (PROVED).** `Y^C_Q(-z)^T eta P(z) = c(z) Y^C(-z)^T eta` for every regular gauge; `= Y^C(-z)^T eta` for the canonical one.
*Proof.* (G) at `-z`: `Y^C_Q(-z) = c(-z)^-1 P(-z) Y^C(-z)`; multiply by `eta P(z)` and use Proposition S. QED

## 5. The u-expansion: residue relations to all orders

Write `P = Sum_n u^n P^(n)(z)` (canonical), `Y^C(v(u), z) = Sum_n u^n eta_n(z)`, `Y^C_Q = Sum_n u^n theta_n(z)`, `eta_0 = theta_0 = -e_1`. Then
`(G)_n`: `P^(n) e_1 = Sum_{j<n} P^(j) eta_{n-j} - theta_n` and `(G*)_n`: `5 e_2^T P^(n) = Sum_{j<n} theta_{n-j}(-z)^T eta P^(j) - eta_n(-z)^T eta`
(PROVED: coefficients of `u^n`). So **the H-column and the H^2-row of every `[u^n] P` are fixed by the lower orders**: this is the all-order
form of the residue relations. Data (frame (0,0); `v = u + (5/9)u^2 + (220/567)u^3 + ...`, `t_c = u + (7/10)u^2 + (41/75)u^3 + ...`):
`eta_1 = theta_1 = -e_0/z`, `eta_2 = -(5/9) e_0/z + e_4/(2z^2)`, `theta_2 = -(7/10) e_0/z + e_4/(2z^2)`; `Y_2` starts at `u^4`: `eta_4, theta_4` have
`-(1/24) z^-4` in slot 2.
* `n = 0`: `P(0,z) e_1 = e_1`, `e_2^T P(0,z) = e_2^T`: (R1) in the canonical normalisation.
* `n = 1`: `P^(1)_{a1} = (delta_{a0} - P^(0)_{a0})/z` and `e_2^T P^(1) = (e_3^T P^(0) - e_3^T)/z`. Row 2: `P^(1)_21 = -P^(0)_20/z = 0`, so `P_21 = O(u^2)`
  (**R2**); row 3: `P^(1)_31 = -pi/z` (**R3b**); column 0 of the dual: `P^(1)_20 = pi/z` (**R3a**); `pi(z) := P_30(0, z)`.
* `n = 2`, row 2: `P^(2)_21 = -P^(1)_20/z = -pi/z^2` (**R3c**).
* `n = 3`, row 2: `P^(3)_21 = -P^(2)_20/z - v_2 P^(1)_20/z + P^(1)_24/(2z^2)` (value `-(336179/6733125) z + ...` in frame (0,0)).
* `n = 4`, row 2: `P^(4)_21 = -P^(3)_20/z - v_2 P^(2)_20/z - v_3 P^(1)_20/z + P^(2)_24/(2z^2) + v_2 P^(1)_24/z^2 - P^(1)_23/(6z^3) + (1 - P_22(0,z))/(24 z^4)`;
  the last term is the normalisation: zero for the canonical gauge, `(43397/19391400000) z + ... = (c_5/24) z + ...` for the stored one.
Canonical-gauge values (frame (0,0)): `pi(z) = (127/3150) z^3 + (1416472237/725359556250000) z^8 - ...` (the `z^8` coefficient differs from
UNIFORM_REDUCTION's `-77141089/362679778125000`, which is the stored normalisation; the two differ by `c_5 . 127/3150`, checked); frame
(1/3,-2/7): `pi = -(239/450) z^3 + ...`. VERIFIED (`gi_relations_*.log`): (G)_n, (G*)_n for `n <= 8` on all complete coefficients (122 and 165
identities per frame), (R1)-(R3), the `n = 3, 4` relations, both frames, exact.
**The diagonal (PROVED; VERIFIED `u^37`, both frames, `gi_formulas.log` (D)).** Order `z^1` of the gauge equation, diagonal: `p_{0,a}' = b_a p_{1,a}
- bt_{a-1} p_{1,a-1}`; with weights `1/p_{0,a}` (which satisfy `lambda_{a+1}/lambda_a = b_a/bt_a`) the right sides cancel, so `prod_a p_{0,a}` is
constant, `= 1` canonically. With `p_{0,a+1} = p_{0,a} bt_a/b_a` and `v'^5 = v/(u (1-u)^5)` (equal products) this gives
`p_{0,0}(u) = tau^(-1/4)/varpi_0`, hence `Ghat(u, 0)_00 = 1/xhat`, `xhat = varpi_0^cf`.

## 6. The actual R-matrix of the formal quintic (LP18)

**What connects the canonical gauge to the actual R-matrix.** UNIFORM_REDUCTION s.4: `Q D(z) = P Q^m(v(u))` with `D` a constant diagonal series;
under **(H)** (purity of the hatted rows) `D` is a scalar series in `z^5` and is absorbed into the gauge: the actual R-matrix defines
`P^act = P^can c^act(z^5)` with `Qhat = Ghat^act R^(tau zhat)`. Since `R~` and `M_0` are symplectic, `G^act` and hence `P^act` are symplectic, so by
Proposition S `c^act(z) c^act(-z) = 1` and (`Ghat^act(0, 0) = 1`) `c^act = exp(a_1 z^5 + a_3 z^15 + a_5 z^25 + ...)`. **"The actual R-matrix is in the
canonical normalisation" means `a_1 = a_3 = ... = 0`.** The `z^5` value was the constant `c_5` of UNIFORM_REDUCTION s.4 (one number); SYNTHESIS
item 8 read `z^5..z^20` off rule W. Here the test is direct and complete in the window.
**Test (VERIFIED, `gi_actual_norm_ab_0_0.log`, `gi_actual_norm_ab_1o3_m2o7.log`).** Rows of `natural_form/referee/ref_hat.HatRows` (`zhat <= 24`,
`u < 18`, own checks h_auto, psi_orth, sympl, covar, circ all pass); `Ghat^act = Qhat M_0(tau zhat)^-1 Psi0^-1`, `M_0` from the same rows;
canonical side `ur_gauge.gauge(24, 72)` over Q, normalised, hatted, reduced mod p. Results, at p89 and p61, frames (0,0) and (1/3,-2/7):
(0) `Ghat^act` is Z/5-graded with `[zhat^k] Ghat^act_{b+k,b} = O(u^mu(b,k))` (this is (H) + Theorem U in the window); (1) `Ghat^act = Ghat^can` in all
1058 slots (993 resp. 997 non-zero); (C) control: the stored normalisation differs in 603 resp. 607 slots, first at `zhat^5 u^4`, entry (0,0);
(2) `[u^4q zhat^5q] Ghat^act_11 = [u^4q zhat^5q] Ghat^act_22 = 0` for `q = 1..4`, while the stored values are non-zero. So `c^act = 1 + O(z^25)`;
`z^5` and `z^15` are the independent tests, `z^10` and `z^20` follow from them by symplecticity.
**Status.** (H) and `c^act = 1` to all orders: CONJECTURAL. Both compare two normalisations of the formal WKB columns: the QRR constants at
`q = 0` (only those at `z^(5q)` are free; the others are automatic, SYNTHESIS item 7) against the conifold normalisation (model stationary phase
for `R^`, the vanishing cycle for `P`). A proof needs a connection formula (e.g. the thimble decomposition of the vanishing cycle in an
oscillatory-integral representation of the formal quintic, UNIFORM_REDUCTION s.8, a SKETCH there); nothing formal fixes `a_1, a_3, ...`.

## 7. Verification table (all logs re-produced today; scripts in this folder)

| check | script / log | result |
|---|---|---|
| (G): positive powers vanish, `P Y^C = Y^C_Q`, quintic system, `[z^-1] = -t_c` | `gi_check.py` / `gi_check_ab_0_0.log`, `gi_check_ab_1o3_m2o7.log` | ALL PASS, Q, `z <= 20`, `u < 29`, 176 + 1038 + 384 identities per frame; stored normalisation fails (127 of 176) |
| (S), (S') exact defect, (G)_n, (G*)_n, (R1)-(R3), u^3, u^4 relations, float diagnosis | `gi_relations.py` / `gi_relations_ab_*.log` | ALL PASS, Q, `n <= 8`, `z <= 20`, both frames |
| the first version's (G*) failure reproduced with its own bound and sign | `gi_diag_old.py` / `gi_diag_old.log` | 73 of 119 fail as floats (rounding only); 0 with an integer sign |
| valuations of `Y^C_Q`, `Y^C`; `p_{0,0}`; closed form of phi | `gi_formulas.py` / `gi_formulas.log` | ALL PASS, Q, `n <= 45`, `u^37`, `m <= 80`, both frames |
| actual gauge = canonical gauge, `c^act = 1 + O(z^25)` | `gi_actual_norm.py` / `gi_actual_norm_ab_*.log` | ALL PASS, p89 and p61, both frames, 1058 slots each |
| hatted (G) on the referee's G table over Q; leaf; dictionary regression | `gi_hat_table.py` / `gi_hat_table.log` | ALL PASS (frame (0,0), `zhat <= 16`, `u <= 11`) |
| (L), (E), control, (T) on edge/leaf tables of the actual R' | `gi_hat_edge.py` / `gi_hat_edge_{ref,prod}_{89,61}_ab_*.log` | ALL PASS, 8 runs (ref: `zhat <= 24`, `u < 18`; prod: `zhat <= 6`, `u < 6`) |

## 8. Diagnosis of the (G*) failure, in detail

The first version (`gi_relations.py` before today's edit) built `Y(-z)^T eta` with `{m: 5 c (-1) ** m}`; for the negative z-powers `m < 0` of the
vanishing solutions `(-1) ** m` is a Python float, so every product and sum involving it became a float, and float-versus-float or
float-versus-Fraction comparisons failed where rounding did not cancel. Reproduced with the first version's completeness bound (`2n + 2`):
119 identities, **73** failing, all floats, residues `|lhs - rhs|` in `[5.1e-21, 8.9e-16]` (first four: `(1,0,2)`, `(1,2,0)`, `(1,2,10)`, `(1,2,15)`, the
slots in the brief), and 0 failing with an integer sign (`gi_diag_old.py`, `gi_diag_old.log`; the same comparison is built into the new log as DIAGNOSTIC, with the
sharper bound `floor(5n/4)`: 165 identities, 112 resp. 101 float failures, all rounding, 0 exact failures). Neither truncation nor the statement
was at fault; the statement is a theorem (s.4). Separately, the old control "(S') first failure at z^5" was a wrong prediction (s.4).

## 9. Files

`gi_core.py` (systems, vanishing solutions, gauge application; unchanged), `gi_check.py` (unchanged), `gi_relations.py` (rewritten: integer sign,
exact (S') prediction, diagnostic), `gi_diag_old.py` (new: the old failure reproduced), `gi_formulas.py` (new), `gi_actual_norm.py` (new), `gi_hat_core.py`, `gi_hat_table.py`, `gi_hat_edge.py`
(unchanged; reviewed), logs `*.log` as in s.7, `GAUGE_IDENTITY_HAT.md`. Project modules are imported read-only (`ur_gauge`, `fq_conifold`, `fq_series`,
`fq_fp`, `fq_sser.EdgeLeaf`, `fq3_par`, `natural_form/referee/ref_hat`).
