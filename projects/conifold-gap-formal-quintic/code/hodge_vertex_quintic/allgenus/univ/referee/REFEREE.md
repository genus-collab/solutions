# Referee report on the universality experiment (allgenus/univ): E4/U4, E3, E5

Date 2026-11-13. Adversarial referee, AI agent; nobody else has read this. Status labels as in the project: PROVED, VERIFIED (exact
computation; genus, orders, prime, frame stated), CONJECTURAL, REFUTED; verdicts on the worker's claims: confirmed / weakened /
refuted / unverifiable. Finite computations are evidence, not proof. Own code `rf_*.py` in this folder (each < 250 lines, ASCII,
`py -I -B`, stdlib only). It imports NO module of the worker (`un_*.py`). Project modules imported read-only: `fq_fp` (primes and
the project's zeta), the row pipeline `fq_rows/fq_conifold/fq_sser/fq3_msum` ONLY to compare rows, and the module-disjoint engine
`natural_form/genus4_referee/indep/` (`ind_tables`, `ind_graphs`, `ind_sum`; it shares no module with `fq4_split` or `un_engine`).
Logs: `rf_g_run.log`, `rf_e4.log`, `rf_e3.log` + `rf_e3_b.log`, `rf_e5.log`, `rf_p61chk.log` (appended; background runs interleave). Rows, parts and tables were
written to the scratch folder only. At most two of my processes ran at any time. Prime 2^89-1 unless 2^61-1 is stated; gauges
exact over Q.

## 0. Verdict: sound with corrections

Every central number of the worker that I re-derived by an independent route came out identical, residue by residue; the
controls fail where they should; I found no tuned convention and no test that cannot fail. Corrections (wording/scope, not results):
1. "Conjecture A + (G) for the canonical gauge implies Conjecture N" omits two inputs: purity (H) of the hatted rows and the
   identification of the ACTUAL R-matrix's gauge with the CANONICAL one (UNIFORM_REDUCTION s.4: verified at the single order z^5,
   which is all that genus 4 and 5 can see). Weakened, not wrong; the worker's s.5 item (C) states it correctly, the claim list
   and Conjecture A's last sentence do not.
2. "Flatness (the Theorem-U gauge equation) is not sufficient" is true but holds BY CONSTRUCTION: Theorem U leaves P -> P c(z^5)
   free and every c gives a flat gauge; the stored normalisation is one such c. The experimental content of b, d is the size and
   pattern of the failure, not the failure itself.
3. The delta_6 consistency is labelled "PROVED (first order)" and "VERIFIED non-linearly"; it is in fact PROVED non-linearly by a
   two-line argument (s.3.4); the worker's one-line first-order formula is that identity truncated at k <= 9 (right in its window).
4. Hard rules: the task allowed at most 6 genus-4 runs; 14 were made (6 full, 8 with maxout 2; the worker discloses this and
   argues by CPU time). e4.log (444 lines) and un_e2.log (347) exceed the 250-line limit if logs count. Nothing else found.

## 1. Independent route

* `rf_series.py`: exact Laurent series over Q with tracked absolute precision. `rf_pf.py`: own Frobenius solution of the conifold
  PF operator (only the operator is taken from the docstring of fq_conifold): t_c = u + (7/10)u^2 + (41/75)u^3 + (9/20)u^4 + ..
  (frame (0,0)), u + (11/30)u^2 + (751/1575)u^3 + .. (frame (1/3,-2/7)): the values on record.
* `rf_gauge.py`: the regular gauge by DIRECT undetermined coefficients in Laurent series (order z^k: the singular cyclic operator L_k
  with kernel w^(k) and cokernel lam^(k); particular solution plus kernel coefficient kappa_k, fixed one order later by the
  regular-singular ODE kappa' (lam.w) + kappa (lam.w') + lam.q' = 0). The ODE residues come out 0, 1/5, -3/5, -2/5, -1/5 (k = 0..4
  mod 5) as in Theorem U; free constants only at k = 0, 5, 10, 15. Not the K[y]-projection of ur_gauge. Normalisations imposed inside
  the recursion: stored p_{k,0}(0) = 0, canonical p_{k,1}(0) = 0; overall scale from d_0(0) d_3(0) = 1.
* `rf_gchk.py`: model vanishing solution from the closed formula, system vanishing solution by integration, (G) coefficient by
  coefficient in z (with an asserted bound on the z-truncation), symplecticity P(-z)^T eta P(z) = eta.
* `rf_rows.py`: symbol rows from the LG closed recursion (LG_MODEL s.3), hatted rows Ghat Rhat(tau zhat), left dressings.
* `rf_indep.py`: graph sums of arbitrary rows on the module-disjoint engine (its row class replaced in memory; the engine asserts the
  Z/5 selection rules, exact division by z + w (symplecticity), [u^0] E_R = E_M0 and E_out = O(u)); own re-expansion in T.
* `rf_e3core.py`: own dressings D = Cayley(X) = (1 - X/2)^-1 (1 + X/2) (the worker used exp), column-1 orbits fixed by a Jacobi
  iteration, other seeds; symplecticity, flat regularity, annihilation and non-constancy of delta_0 are CHECKED on D itself.

## 2. E4/U4 and E2: gauges, normalisation, genus-4 universality (all VERIFIED here unless stated)

**2.1 Gauges (exact over Q; `rf_g_run.py`, `rf_g_run.log`, all PASS).** For the quintic in both frames, the worker's c1, c2, c3, c4
and two NEW systems (r1: F = (2, -1/4, 5/3, 1/6), Gm = (3/2, 2/9, -1, 1/5), frame (0,0), so f(0) = 2 AND g(0) = -3/2, which the worker
never varied; r2: F = (1/2, 3, -2), Gm = (1, -1/3, 4), frame (1/3,-2/7)), with k <= 16: the gauge equation holds to the tracked
precision; every p_{k,b} is a power series (Theorem U); the STORED constants of P_11(0, z) printed in the worker's e4.log are reproduced
exactly (quintic: -43397/807975000 at z^5 and -765368951400317/1305647201250000000 at z^10, identical in both frames; c1, c2 (incl.
P_11(0,0) = 2/3), c3, c4); canonical normalisation imposed inside the recursion = stored times c0/P^stored_11(0, z) (post hoc, the
worker's definition); (G) P Y^C_model = c0 Y^C_system holds on z^-8..z^8 mod u^8 in canonical normalisation for all eight systems,
and fails at z^1..z^8 in the stored one; the canonical gauge is symplectic through z^16, the stored one only through z^9 (first
defect at z^10, where c(z)c(-z) = 1 - (m5^2 - 2 m10) z^10 + ..). So symplecticity fixes only the even part of log c(z^5); (G) fixes the
odd part. Remark: c0 = 1/lambda holds only for g(0) = -1; in general c0 = 1/(f(0)|g(0)|) (r1: 1/3).

**2.2 Rows.** My canonical quintic rows (own PF, own gauge, own symbol) equal the project rows in all 25 x 10 x 6 residues at genus 4,
frames (0,0) and (1/3,-2/7) (2^89-1): the worker's case a/e regression is confirmed by a route that shares nothing with ur_gauge.

**2.3 Graph sums on the module-disjoint engine (`rf_e4.py`, `rf_e4.log`).**

| system (canonical) | genus 2, 3 (full) | genus 4 | vertex graph |
|---|---|---|---|
| quintic (0,0) / (1/3,-2/7) | regular | (0,0), E <= 2: = worker's e4_tab_a_89_out2.json (11 rows); also regular at 2^61-1 | c_g t_c^(2-2g) + 0 |
| c1 | regular | **full: all 378 edge graphs regular; all 379 rows = worker's e4_tab_c1_89.json residue by residue**; E <= 2 also at 2^61-1 | c_4 (T/lambda)^-6 + 0 |
| c2 (lambda = 3/2) | regular | E <= 2: = worker's e4_tab_c2_89_out2.json | c_4 (T/lambda)^-6 + 0 |
| r1 (NEW: f(0) = 2, g(0) = -3/2) | regular | **full: all 378 regular** (vertex [u^0] = 81/1280000 = c_4 (3/2)^6) | c_4 T''^-6 + 0, T'' = c0 Int f |
| r2 (NEW: frame (1/3,-2/7), f(0) = 1/2) | regular | **full: all 378 regular** | c_4 T''^-6 + 0 |

**2.4 Normalisation response.** P_canon x (1 + m z^5) with m = 1 (my choice; the window sees only z^5): the differences to the canonical
run, per unit m, equal the worker's (b - a)/m5 and (d - c1)/m5 tables EXACTLY on all seven graphs that move (vertex; graphs 1, 2 from
[u^4]; 3, 4, 7, 8 at [u^5]), at 2^89-1 and 2^61-1 (`rf_p61chk.py`, `rf_p61chk.log`; so the worker's rational reconstructions are right).
The sums over graphs vanish at [u^4] and [u^5] (exact rationals). The valuation argument of the worker's s.2.3 was re-derived: rows
change by m u^4 Q[n-5]; xhat (zhat^0) and E_M do not move; an outer graph with E edges gets O(u^(E+3)) from one changed edge and
O(u^(E+4)) from one changed leaf, products of changes are O(u^8): CONFIRMED (PROVED, elementary).

## 3. E3 and E5 (own dressings D = Cayley(X); `rf_e3.py`, `rf_e3.log`, `rf_e5.py`, `rf_e5.log`)

**3.1 Non-linear (VERIFIED, 2^89-1; genus 3 also 2^61-1).** Every D below was checked to be symplectic and flat-regular; annihilation
was imposed for m != 6 and CHECKED for m = 6; delta_0 is never constant (no z-independent target connection).

| D | genus 2 | genus 3 | genus 4 |
|---|---|---|---|
| A1 (k >= 1) | 6 regular | 41 regular (both primes) | **full: 378 regular, at 2^89-1 and 2^61-1** |
| A3; K2 (with u-dependent k = 0 terms) | A3: 6 regular | 41 regular | K2, E <= 2: 10 regular |
| mA7, mK8: PURE MODEL base (P = 1, model covector) | 6 regular | 41 regular | mA7, E <= 2: 10 regular (beyond the worker: genus 4) |
| B1 = A1 + (1/7) u^2 zhat H(x)H (control) | regular (u^2 outside the window) | 3 singular, [u^2] = 1/14000, 1/252000 (= worker's B2 values) | E <= 2: all 10 singular |
| mB7 (pure model, same violation) | | 3 singular, same [u^2] values | |
| R1: random column-1 orbits (control) | | 16 singular, vertex law broken | |
| C1c-type: A + two terms irregular by one order, re-annihilated (control) | | 4 singular at [u^3]; annihilating | |
| xA1: annihilates the MODEL covector, not the quintic's (control) | regular (no column-1 entry visible) | 2 singular at [u^3] | |

In every regular case the vertex graph is c_g T_eff^(2-2g) + 0, T_eff = -[z^-1 e_0] D Y^C_Q. The [u^2] response of the H(x)H violation
does not depend on the random part of D, which is why my B1 reproduces the worker's B2 numbers.

**3.2 First order (VERIFIED, genus 2 and 3, both primes).** For random X in the kernel of the LINEARISED annihilation map (k >= 1, and
with k = 0 terms), D_t = Cayley(t X) at t = 0..4: since X = O(u) the polar coefficients are polynomials of degree <= 2g-3 in t (the
t^4 coefficient vanishes, checked); their t^1 coefficients vanish on all 41 (genus 3) and 6 (genus 2) edge graphs (cases L1, L2, N4,
L3). Control P5 (random column-1 orbits): non-zero on 9 graphs. One random kernel vector with zero response shows ker Ann <= ker Resp
except with probability about 1/p. Not re-derived: the worker's ranks (genus 2: rank Ann = rank Resp = 3; genus 3: rank Ann = 10,
rank Resp = 8, dim ker Resp = 22) and the directions with a u^0
hatted part (the engine requires E_out = O(u)); for these: annihilation forces the u^0 part of zhat H(x)H to vanish (Delta_1 at u^0),
and constant symplectic diagonal D leave E_R invariant and rescale the genus-g sum by d_0^(2-2g) (T_eff -> d_0 t_c): trivial.

**3.3 E5 (VERIFIED).** The worker's surplus directions w1, w2 (coefficients copied from un_e5.log as data): genus 3: all 41 regular,
vertex T^-1 coefficients 11/136080 and -11/54432, and 5 w1 + 2 w2 keeps the vertex law: identical. Genus 4 (E <= 1, no completion
needed: the [u^3] values of one-edge graphs see only the u^3 part of X): w1: -1/378000, 13/9720000, -1/283500, vertex T^-3 19/56700000;
w2: 1/151200, -541/97200000, -139/567000, -19/22680000: identical to the worker's.

**3.4 delta_6 is implied NON-LINEARLY (PROVED).** Let D be symplectic, D(-z)^T eta D(z) = eta, and Y = Y^C_Q (only z^(<= 0), z^0 term
-e_1). Put W = D Y = Sum_m delta_m z^m e_{m+1}. Then W(-z)^T eta W(z) = Y(-z)^T eta Y(z) has no positive power of z. If delta_m = 0 for
all m >= 1 except m = 6, the z^6 coefficient of the left side is ((-1)^6 + 1) delta_6 delta_0 eta_{21} = 2 delta_6 delta_0, and delta_0 =
-1 + O(u) is a unit: delta_6 = 0. By induction the same holds for every self-dual even m (16, 26, ..). In a truncated window the identity
fails by terms of order zhat^(>= KZ+1) paired with negative powers of Y, harmless for flat-regular D (my A1: delta_6 = 0 mod u^6) and
not for flat-irregular D (the worker's crashed case C1). This upgrades the worker's "PROVED (first order) / VERIFIED (non-linear)".

## 4. Hard rules and process

* Files: the worker's 48 files (plus desktop.ini) are inside allgenus/univ, and by inspection its scripts write only there; no file of
  the project outside allgenus changed during its run (modification times, 17:40-19:40); all ASCII; no e-mail address or person's name (only standard terminology). py -I -B throughout.
* Budget: "at most 6 genus-4 runs in all" was exceeded: 6 full + 8 maxout-2 runs (disclosed; about 227 CPU min, i.e. within six
  full-size runs by CPU). The "at most 6 processes" claim could not be checked from the logs (pool of 4-5 workers plus single runs).
* e4.log (444 lines) and un_e2.log (347) exceed 250 lines (if the rule covers logs). un_e3_nl_89.log and un_e3_c_89.log carry
  hand-added NOTE lines in otherwise program-written logs (they are marked as such and explain crashed/invalid blocks; acceptable).
* My own process: two development bugs (PF monomial call, cokernel index) were fixed before any logged run; the first B1 genus-2
  attempt in rf_e3.log stops after its header (index error: the u^2 term lies outside the genus-2 window), rerun after a guard.
* Inherited, common-mode caveat (both engines): the regrouped split sum is identified with the outer graph sum by an argued,
  genus-independent identity, numerically checked against genuine vertex tensors only at genus 2-3 (INDEP.md s.5).

## 5. Second prime, genus 4, and an extra control

* Own non-gauge dressing A1 at genus 4, FULL run, also at 2^61-1: all 378 edge graphs regular, vertex c_4 T_eff^-6 + 0 + .. + 0,
  delta_6 = 0 mod u^6 without being imposed. So the genus-4 statement of E3 now rests on two primes and on dressings of two
  different constructions (worker: exp, Newton; referee: Cayley, Jacobi).
* Wrong-covector control (not in the worker's set): D built to annihilate the MODEL's vanishing covector instead of the quintic's
  (`rf_e3.py xA1`): genus 3: annihilation of the quintic covector fails at m = 1, 2 and two one-edge graphs are singular at [u^3]
  (91/810000, -91/6480000). The principle is about the covector of the theory actually dressed, as the worker states.

## 6. What remains open (not defects of the worker's report)

Conjecture A is unproved and tested only at genus 2-4 (one or two primes) on random dressings and gauges; B = prod b_a (hence tau,
v(u), M) was never varied by anyone; the first-order rank statements and the flat-irregular first-order ranks (worker's s.3.1) were
not re-derived here (the module-disjoint engine needs E_out = O(u)); the identification of the actual R-matrix's normalisation with
the canonical one is a single verified number (z^5).

## 7. The worker's claims, one by one

| # | claim (worker's list) | verdict | basis here |
|---|---|---|---|
| 1 | canonical quintic gauge = project rows / natural run, both frames | confirmed | s.2.1, s.2.2 (own route), E <= 2 table |
| 2 | genus-4 universality c1, c2 (378 regular, vertex c_4 (T/lambda)^-6) | confirmed | c1 full: 379 rows identical; c2 E <= 2; NEW r1, r2 full |
| 3 | stored normalisation: exactly graphs 1, 2 [u^4], 3, 4, 7, 8 [u^5]; total keeps the gap | confirmed | s.2.4, both primes, valuation argument |
| 4 | valuation argument (PROVED) | confirmed | re-derived, s.2.4 |
| 5 | per-unit-m responses, [u^4] universal, sums vanish | confirmed | identical rationals at two primes |
| 6 | first order: ker Ann in ker Resp (genus 2, 3) | confirmed (inclusion) | s.3.2; ranks not re-derived |
| 7 | first order without flat regularity: ranks 5, 23 | unverifiable here | engine limit; qualitative point supported by C1c-type |
| 8 | non-linear A1-A4 regular at genus 2, 3 | confirmed | own dressings, other construction, two primes |
| 9 | two random dressings regular at genus 4 | confirmed | own A1 full at two primes; K2, mA7 at E <= 2 |
| 10 | violations (B1, B2, b, d, C1a-c) singular | confirmed | B1, R1, C1c-type, B1 at genus 4; plus wrong-covector control |
| 11 | gauge equation neither necessary nor sufficient | weakened | "not sufficient" is by construction (s.0 item 2) |
| 12 | genus-3 surplus singular at genus 4 | confirmed | s.3.3, identical numbers |
| 13 | pure-model dressings regular (genus 2, 3) | confirmed | mA7, mK8; extended to genus 4 at E <= 2 |
| 14 | delta_6 implied (PROVED first order) | confirmed, strengthened | s.3.4: true non-linearly |
| 15 | Conjecture A, and A + (G) => Conjecture N | weakened | needs (H) and the actual-normalisation identification (s.0 item 1) |
