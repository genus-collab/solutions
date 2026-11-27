# SYNTHESIS -- natural form of the conifold factorisation of the formal quintic (completeness critic, 2026-11-11)

Sources read: `lg/LG_MODEL.md`, `theory/UNIFORM_REDUCTION.md`, `referee/REFEREE_NAT.md`, the three verifier folders (`verify_lg`, `verify_theory`,
`verify_referee`), the three production logs `logs/genus3_natural_*_nat*.log` (read-only; all three runs had FINISHED: 23:28 and 23:50 on
2026-11-10), `OU_EXTENSION_genus3.md` s.6d. Primes: p89 = 2^89-1, p61 = 2^61-1, p521 = 2^521-1 (referee only). Labels: PROVED-mod(X) = argument
written, AI-audited only, depends on X; VERIFIED = exact computation (orders, primes stated); QUOTED = literature, not checked; CONJECTURAL.
Nothing here has been read by a human. Two light checks of my own are in `critic/` (s.5); they are marked [critic].

## 0. Short answer to the question

- **What the genus-3 run shows (VERIFIED).** With `M = M_0(tau(u) zhat)`, `tau = 1/2F1(1,4/5;9/5;u)`, and nothing else in the vertex factor,
  all 41 edge graphs of genus 3 are `O(u^4)` and all 6 edge graphs of genus 2 are `O(u^2)`; the vertex graph carries `(1/25200) t_c^-4 + 0 + 0 + 0`
  resp. `(1/1200) t_c^-2`. Frame (0,0): p89 and p61. Frame (1/3,-2/7): p89 only. `M` taken mod `(u^4, zhat^7)`.
- **Natural form (yes, as a conjecture with a precise shape).** `R = [regular gauge Psi(u) G, computable from the periods] x [R-matrix of the
  Landau-Ginzburg model x^4 + v/x, primitive form dx, evaluated at z over the canonical coordinate measured from the conifold point]`. This is
  the same uniform-asymptotics pattern as the Hodge factor (one saddle next to a logarithm, variable z/w) with five saddles coalescing on a
  simple pole and variable `tau zhat = -(5/4) zeta^i z/(u_i - u_i^c)`.
- **What is NOT yet the Hodge trick.** (i) Regularity is known at genus 2, 3 only, with no mechanism beyond first order. (ii) The vertices have
  no closed evaluation: the model's `F_g` is known only through the project's own frozen sums (g <= 3). (iii) The graph sums through genus 3
  do not single out this `M` (s.6d has a large solution family); it is singled out by the reduction principle (uniqueness of tau, regular gauge).

## 1. Statements that stand

| # | statement | status |
|---|---|---|
| 1 | `R~(zhat;u) = G(u,zhat) M_0(tau zhat)`, `[u^n]G` polynomial in zhat, `M_0` on the RIGHT | VERIFIED zhat<=16, u<=11 at p521, p89, p61 (referee's from-scratch rows); same on the PROJECT rows at p89 (zhat<=16, u<=11) and p61 (zhat<=12, u<=8); zhat<=24, u<=17 at p89 and, [critic], p61. All orders: PROVED-mod((H), Theorem U) |
| 2 | Theorem U: two cyclic systems with one simple pole and equal products are related by a gauge regular at u = 0, unique up to `c(z^5)` | PROVED (written proof, audited line by line by a second AI and re-implemented along the proof route, k <= 10; [critic] k <= 14 with all holomorphy assertions passing). Unrefereed by a human |
| 3 | tau is unique and equals the canonical coordinate; any other rescaling, and `M_0` on the left (scalar rescalings, zhat<=16), fail | uniqueness PROVED-mod(hatted R-equation to all orders, i.e. Picard-Fuchs + (H)); failures VERIFIED three primes |
| 4 | `M_0` = stationary-phase data of `x^4 + v/x`: `R^_a^(i)(zhat) = N_a zeta^(i(a+1)) J_{c_a}(-4 zeta^-i zhat)`, `c = (0,-1,3,2,1)`, `N = (1,-1,-1,1,1)`; closed recursion to any order | VERIFIED all rows zhat<=18 (p89, p61), identity row zhat<=40; two independent stationary-phase engines. All orders: PROVED-mod((a) WKB coefficients are Laurent polynomials in L, (b) the conifold-frame map respects D on the row subring, a standard stationary-phase lemma). Independent of (H) |
| 5 | Dictionary: `H^a <-> x^(c_a) dx` (1, H, H^2, H^3, H^4 <-> dx, -dx/x, -x^3 dx, x^2 dx, x dx); primitive form dx, d = 1/2 | VERIFIED (only k = c_a matches, both primes; the task's guess `x^a dx` fails for a != 0) |
| 6 | The model is the Hurwitz Frobenius manifold of `p^4 + a_2 p^2 + a_3 p + a_4 + a_5/(p - a_1)`, form dp, on the residue line; flat metric, unit `d/da_4`, `c_555 = 1/a_5`; H is the residue direction | VERIFIED locally (exact at 4 random points; floating point at 6). Identification with H_{0;1,4} of arXiv:2211.12259 and its KP tau-function theorem: QUOTED |
| 7 | (H) purity: the hatted rows are power series in u | VERIFIED mod (z^17, u^12) at p89 and (z^13, u^9) at p61 on the project pipeline (purity is asserted in `HatRowsU`); identity row: QRR normalisation holds automatically for all n <= 24 not divisible by 5 (referee code, p89; [critic] p61). All orders: CONJECTURAL |
| 8 | Canonical normalisation: the gauge of the actual R-matrix has `P_11(0,z) = P_22(0,z) = 1` | VERIFIED at z^5, z^10 (three primes), z^15, z^20 (p89; [critic] p61): read off rule W, see s.2 item 4. [critic]: the canonical `c_10` equals the symplectic value `153074543599907/261129440250000000`. All orders: CONJECTURAL |
| 9 | Periods-only gauge + model series reproduce the actual G | VERIFIED: 25 entries, zhat<=6, u<=5, two frames, both primes (theory + verifier); [critic] ALL 397 Fourier components (392 non-trivial) for zhat<=14, u<=11, exact rationals, NO fitted constant (c_5, c_10 from item 8) |
| 10 | Rule W as a BOUND: `G - 1` is regular at u = 0 in the frame (1, uH, H^2, H^3, H^4); hence `deg_zhat [u^n]G <= floor((5n+2)/4)` | VERIFIED zhat<=24, u<=17 (p89, [critic] p61); 7 frames at zhat<=12, u<=8 (p89). [critic] bookkeeping: rule W is EQUIVALENT to Theorem U + (R1) + (R2) + item 8, so PROVED-mod((H), item 8) |
| 11 | (R1)-(R3): `P(0,z)` commutes with the residue `E_12`; no `H(x)H` in `G_0, G_1`; the H-block starts as `alpha zhat (H - zhat)(x)(H + zhat)` at u^2 (two Euler-dilaton vectors); chain relation at weight 0 | (R1)-(R3) PROVED-mod(Theorem U), VERIFIED k <= 14 over Q; general chain relation VERIFIED 13/13 per leg (three primes), 9/9 in 7 frames |
| 12 | `G_1 = zhat X_1`, components `(alpha - 8/45)/5` (string type, frame dependent) and `-1/75` (frame independent), no `H(x)H` | PROVED by hand from the G-equation mod the hatted R-equation; VERIFIED 7 frames (p89), 2 frames (both primes) |
| 13 | Frame dependence: `M_0`, tau, polynomiality, rule W (bound), chain relation are frame independent; all values, sharpness and degrees are not; `H(x)H` coefficient of `[u^2 zhat]G = -127/15750 - (2/5) beta` | VERIFIED 7 frames, p89 (zhat<=12, u<=8). The beta-formula also matches the production log of the second frame (239/2250) |
| 14 | Graph sum, genus 2 and 3, regular with `M_0(tau zhat)`; first-order correction is the point `-X_1` of the plane Pi | VERIFIED (s.0). Content beyond kinematics: 9 rational numbers per frame (`[u^2],[u^3]` of the 2 one-edge graphs, `[u^3]` of the 5 two-edge graphs), with no free parameter |
| 15 | Euler-dilaton element: zero first-order response on one-edge graphs with a vertex of genus >= 2; NOT an exact symmetry (genus-1 anomaly -1/8) | VERIFIED genus 2, 3, both primes; general statement conditional on the primary sector at h >= 3. Consistency seen in the logs: genus-2 bridge graph `[u^2] = (H(x)H coefficient)/128` in both frames |
| 16 | Genus one from the model alone: `dF_1 = -(1/12) du/u` | VERIFIED (arithmetic, both primes for `C_1`) |
| 17 | Conjecture N (all g), model `F_g = B_2g/(2g(2g-2)) a_5^(2-2g)` for all g, primary-sector lemma, oscillatory-integral (mirror) picture | CONJECTURAL |

Remark [critic] (elementary): at order u the K-regular terms of rule W are exactly `zhat x (1(x)H^2 + H^2(x)1, H^3(x)H^4 + H^4(x)H^3)`, i.e. the
kinematic plane Pi of s.6d found from the graph sums. This suggests the working hypothesis "K-regular left factors are harmless" (Lemma K, s.3).

## 2. Contradictions between the reports and their resolution

1. **Dictionary.** Task and referee side remark: `H^a <-> x^a dx` ("consistent with the spectrum"). LG and theory: refuted, `k = c_a`. Resolution:
   the spectrum is symmetric under mu -> -mu and cannot decide; the row-by-row test does. The referee's item is withdrawn. Also "x^(-a) dx mod 5"
   is not an identity of classes (only `x^3 dx = v x^-2 dx`); state `k = c_a` in the window -1..3.
2. **Does a polynomial G exist to all orders?** Referee: CONJECTURAL ("a uniform turning-point reduction, not proved"). Theory: PROVED given (H).
   Resolution: Theorem U is exactly the missing existence statement; the referee wrote before seeing it. Status: PROVED-mod((H), Theorem U).
3. **Frame independence.** Referee: "independent of the conifold frame, degrees floor((5n+2)/4), rule W sharp". Verifier: false in special
   frames (beta = -127/6300: deg G_2 = 2 and the u^2 chain vanishes; alpha = 8/45: the string part of G_1 vanishes). Theory had "frame-dependent
   G" correctly. Resolution: table item 13. floor((5n+2)/4) is a bound attained in generic frames.
4. **How far are (H) and the normalisation tested?** Theory and its verifier: "(H) only mod z^7, u^6; normalisation = one number at z^5;
   z^7..z^9 untested, so the genus-4 shortcut is unsupported". This is too pessimistic: the referee's rows are pure by construction and meet
   the QRR normalisation automatically to z^24, the referee's verifier rebuilt the project rows (purity asserted) at zhat<=16, and rule W
   forbids the `H(x)H^2`, `H^2(x)H` components at `u^(4q) zhat^(5q)`, which IS `P_11(0,z) = P_22(0,z) = 1` at z^(5q), q = 1..4. Table items 7-9.
5. **"The u-corrections are not explained by the model" (LG) versus theory.** Resolved at the level of the Gauss-Manin systems: tau is forced
   by equality of the canonical one-forms and G is the regular gauge (`u_eff` of LG = `v(u)` of theory). Still open: a description inside the
   Frobenius manifold (the formal quintic is not a path in it).
6. **Rule W "iff" (referee) versus bound.** The iff holds in frame (0,0) (880 slots checked); in general only the bound. The two a priori
   bounds agree: floor((5n+4)/4) from the grading, improved to floor((5n+2)/4) by (R1), (R2) (theory) = "no H leg at weight 0" (referee).
7. **Normalisation of the model.** "u_LG = 4u, z_LG = -4z exactly" (LG), "c^4 = -1/4" (theory, referee): the same statement; only
   `eps = -4 zeta^-i zhat` is fixed by the data, the rest is the scaling freedom of the model. The orientation of the integral in the task's
   formula for tau is opposite to the sign used in the notes; immaterial.
8. **Sign in the Euler-dilaton block.** The flat form `-pi z (uH + z)(x)(uH - z)` and the hatted `alpha zhat (H - zhat)(x)(H + zhat)` differ
   because the Fourier vector of H is minus the flat class (`N_1 = -1`), not by an "ordering convention" (theory note); explicit dictionary in
   `verify_theory/v3_fourier.py`.

## 3. What is missing or unverified (ranked) and the next computation

1. **Conjecture N beyond genus 3, and a mechanism.** Evidence is genus 2, 3; beyond the kinematic first order it is 9 numbers per frame. Nothing
   explains the vanishing of the non-H parts of `G_2, G_3`, the `G_1 G_1` terms, the leaves and the tau-dependence of the vertices.
2. **The all-genus theory of the model** (the analogue of the closed Hodge-integral vertex): `F_g` on the residue line (Bernoulli constants)
   and the primary-sector lemma (only H-insertions survive at genus >= 1). Known: `F_1` from the model; g = 2, 3 only as the project's
   frozen sums; primary sector for (1, <=4), (2, <=3), (3, 0). The KP statement is quoted and has not been used for anything.
3. **Is `M_0(tau zhat)` singled out by the graph sum?** Not through genus 3. Whether K-regular parts of G may be moved into the vertex factor
   without spoiling regularity is untested (Lemma K).
4. **(H) and the canonical normalisation to all orders** (now well tested, items 7-8; no proof; equivalent to a connection formula between
   the q = 0 and the conifold normalisations of the WKB series).
5. **Human refereeing** of Theorem U (s.3 of the theory note), of Theorem A of the LG note with its inputs (a), (b), and of the uniqueness
   lemma for tau (last step to be written out).
6. **Second frame at the second prime; provenance.** Frame (1/3,-2/7) graph-sum pass exists at p89 only. The frame-(0,0) logs (23:28) predate
   the last edit of `genus3_natural.py` (23:30). One re-run per case (about 20 min, 13 processes) closes both.
7. Minor: the single-leaf bookkeeping `-(g-1) c_g/15` is observed, not derived; Proposition E for h >= 3; sharpness of the degree; the general
   chain relation; dx versus dx/x (a sketch); orders u = 12..17 rest on the referee's code alone (now two primes).
8. Speculative: a preferred frame `(alpha, beta) = (8/45, -127/6300)` in which `G_1` is the single frame-independent term and the u^2
   Euler-dilaton chain vanishes; untested beyond u^4, no interpretation.
9. Untouched: the mirror derivation (sketch in the theory note s.8) and any statement about the true, non-formal quintic.

**Most informative next computation (cheap): the Lemma-K test at genus 3.** Use the existing `--variant` mechanism of `genus3_natural.py`
(`M = A M_0(tau zhat)`, `A = exp(Sum u^n a_n)`) with a RANDOM symplectic, Z/5-covariant, K-regular dressing (`a_1` a random element of Pi; `a_2`,
`a_3` random terms of rule-W weight >= 0; at u^2 these are 6 parameters). One pass per prime. If all 41 edge graphs stay regular, regularity is
a property of the K-regular class: Conjecture N then reduces to items 7-8 plus a lemma about the model theory alone, and the evidence grows from
9 numbers to a family. If not, the failing components show what the mechanism uses beyond rule W (the G-equation). **The decisive but expensive
computation is genus 4** (379 graphs, `[u^0..u^5]`): first test of the `u^4 zhat^5` diagonal terms (where the z^5 normalisation enters) and of
the genus-1 anomaly in multi-edge graphs. Its inputs now exist without new WKB runs: the gauge and the model series reproduce the actual G to
zhat^14 (item 9), and `M_0` is stored to zhat^40. A cheap prelude is the frozen genus-4 vertex number, predicted `1/180000`.

## 4. Overclaims to remove or reword

- LG note s.6 and summary: "the formal-quintic counterpart of the Hodge factor, carrying the gap". True only for the symbol; `M_0` alone leaves
  simple poles in the one-edge graphs (OU_EXTENSION_genus2). The counterpart is `M_0(tau zhat)` together with the regular gauge.
- LG note: Theorem A "PROVED" -> PROVED-mod((a), (b)); s.5.3 "PROVED" -> sketch; "under dx/x the unit would be d/dt5" is wrong as worded;
  "u-corrections not explained" is superseded (s.2 item 5); gap constants of the Hurwitz theory are inferred, not computed.
- Theory note: Theorem D equality and Proposition E are labelled PROVED; only the bound, resp. the genus-2/3 cases, are. "(H) verified only mod
  z^7" and "normalisation: one number" understate the evidence (items 7-9). "v is unique given v(0) = 0": only up to a fourth root of unity.
- Referee note: "independent of the conifold frame", "deg = floor((5n+2)/4)", "nonzero iff", "H(x)H present from u^2 on, all negative": frame
  (0,0) statements only. All-order existence "CONJECTURAL" -> PROVED-mod((H), Theorem U). Counts: 104 of 120 and 67 of 80 testable pairs.
- Everywhere: "all 41 edge graphs are regular" should carry "M mod (u^4, zhat^7); frame (0,0) two primes, second frame one prime; non-kinematic
  content 9 numbers per frame". "Uniqueness of tau PROVED" should carry "modulo the hatted R-equation to all orders".
- "Natural" must not be read as "forced by the graph sum": through genus 3 a whole family regularises every edge graph (s.6d).
- The KP tau-function statement, the Bernoulli formula for all g, and the Mellin-Barnes comparison are context (quoted, inferred, analogy).
- "Genus 4 needs no new R-matrix computation" is now supported to zhat^14, but still assumes (H) and item 8, which are unproved.

## 5. Critic's own checks (new files in `critic/`, single process, seconds)

- `critic_rerun_ref_run2_61_NZ24_K18.log`: the referee's `ref_run2.py` run unchanged (`py -I -B`, output redirected) at p61 with zhat<=24,
  u<=17: degrees floor((5n+2)/4) for n <= 17, rule W 0 violations and 0 non-sharp of the testable pairs, `H(x)H^2`, `H^2(x)H` absent at
  `u^(4q) zhat^(5q)`, q = 1..4. Second prime for the orders u = 12..17 and for the normalisation at z^15, z^20.
- `critic_gauge_vs_table.py`, `.log`: the verifier's proof-route gauge (`verify_theory/v2_gauge.py`, periods only, K = 14) with the normalisation
  `P_11(0,z) = 1` against the referee's table `referee/G_fourier_NZ16_K12.json`: 397 of 397 components agree as exact rationals, none missing
  on either side; `c_5 = 43397/807975000`; canonical `c_10` = symplectic `c_10`. ALL PASS. (The table itself was confirmed on the project
  rows by the referee's verifier at p89, and for zhat<=12, u<=8 at p61.)
- The equivalence in table item 10 is index bookkeeping: with `w_(3-a) = 2 - c_a + 5 delta_(a,2)` and `w_b = c_b + 5 delta_(b,1)`, rule W for the
  entry `(a,b) = (b+k, b)` reads `n >= mu(b,k) + delta_(a,2) + delta_(b,1)`. Unrefereed.

## 6. Housekeeping

- A verifier's re-run of `lg/lg_gm.py` rewrote `lg/lg_symbol_rows_Q.json` (deterministic, same size, content re-checked; disclosed by that
  verifier). No project file outside `natural_form/` was modified by anyone; I wrote only `SYNTHESIS.md` and `critic/`.
- The LG and theory notes cite papers with author surnames; if rule (3) is meant strictly, replace them by arXiv numbers.
- The production runs were described as still running; all three had finished and passed before the reports were completed.
