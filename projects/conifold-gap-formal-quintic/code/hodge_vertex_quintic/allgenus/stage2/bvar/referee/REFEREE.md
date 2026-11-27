# Referee report: Stage 2, falsifier (c) "B-universality" (`allgenus/stage2/bvar/`)

Date 2026-11-14. Adversarial referee (AI agent). No human has read this. Labels: PROVED (argument written here, or a quoted
theorem), VERIFIED (exact computation; range and primes stated), CONJECTURAL, SKETCH, REFUTED. Verdicts on the worker's claims:
confirmed, corrected, weakened, refuted. Finite checks are evidence, never proof. Conventions as in `../BVAR.md` and
`../../../SYNTHESIS_STAGE01.md` s.3. A genus-g term is `u^(2-2g) G(u)`, and a graph is regular iff `[u^0..u^(2g-3)] G = 0`.

## 0. Verdict: sound with corrections

The falsifier outcome is CONFIRMED and survives every attack I tried. What has to change is how much it is said to prove.

* **Numbers.** I re-derived the worker's VERIFIED numbers with my own construction of the system, of `v(u)`, `tau`, the
  Theorem-U gauge, the normalisations, (G), the hatting, the rows and `E_M` (s.2). This covers:
  * every canonical genus-2/3 run, at both primes;
  * at genus 4, S3 (maxout 2) and S2 in full, plus S2 at the second prime;
  * the controls and the normalisation response.
  My per-graph tables agree with the worker's tables residue by residue, for every graph and every coefficient in the window.
  A second, module-disjoint graph-sum engine gives the same per-graph series (`REFEREE_S4.md` s.4.1). S1 and S5 at genus 4 were
  not re-run.
* **Correction 1 (main).** `B` is not an invariant of the system. Lemma C (s.3, PROVED; the per-graph identity is also VERIFIED by
  a separate pipeline) works as follows. Any one-pole system with `B` perturbed becomes, after the coordinate change
  `u = phi(w)`, `phi = v_X^-1 o v_q`, a system with the QUINTIC `B`, i.e. a univ `(F, Gm)` system. Every outer-graph amplitude is
  carried along: `V~_Gamma(w) = V_Gamma(phi(w))`. So B-universality is EQUIVALENT to `(F, Gm)`-universality of the pulled-back
  systems. A failure could not have revealed "a hidden use of the quintic's B" in the mechanism (BVAR s.7, "What is tested").
  It could only have revealed one of two things: a non-covariant implementation (a wrong `tau`, `v` or hatted prefactor, which is
  what the staleM and staleHat controls show), or a failure of `(F, Gm)`-universality.
  The test is still worth having. The pulled-back `F~, Gm~` are full power series (univ used cubic polynomials), and the run checks
  the coordinate covariance of the whole pipeline.
* **Correction 2.** The row "B-universality in all genera: CONJECTURAL; implied by quintic Theorem R + ML on K_vc (+ Theorem U)"
  is right, but it is not a new conjecture. By Lemma C it is equivalent to `(F, Gm)`-universality in all genera, i.e. to the
  all-genus extension of E12.
* **Correction 3 (minor).** "c0 = P_11(0,0) = 1/lambda" holds only for `m = 1`, `g(0) = -1`. In general
  `abs(c0) = m^10/abs(f(0) g(0))`, with the sign fixed by the choice `sigma > 0`; the vertex law is even in `c0`. VERIFIED for
  S1-S7 here: S5 1024, S6 59049/2, S7 21/20 with `f(0) = -5/7`. The table entry "lambda = 1/1024" for S5 is a relabelling, not `f(0)`.
* **Minor (rules, disclosure).**
  * A 4-worker pool plus its parent made 5 processes for under a minute (disclosed).
  * `logs/g4_batch.sh` contains an absolute scratch path with the local login name. It is not a personal name, but it is best
    avoided.
  * The S5 run without the symplectic scalar is not logged; I reproduced it (`REFEREE_S4.md` s.4.4).
* Nothing is REFUTED.

## 1. The PROVED claim (BVAR s.1, pole order): confirmed, with a scope remark

* (i) A gauge needs `prod b = prod bt`. Correct. The `z^0` equation is `b_(c-1) q^0_c = bt_(c-1) q^0_(c-1)` (`P_(-1) = 0`). An
  invertible `P_0 = diag(q^0)` forces `prod b = prod bt = -v'^5/v`.
* (ii) Let `v` be a power series with `v(0) = 0` of order `n`. Then `ord_0(v'^5/v) = 5(n-1) - n = 4n - 5`, so `ord_0 B` lies in
  `{-1, 3, 7, ..}`. A pole forces `n = 1`. Correct and complete.
* The converse paragraph ("`ord B = -k` forces `v ~ u^((5-k)/4)`; `k = 5` gives a log; `k = 9, 13, ..` give a pole") assumes a
  leading power behaviour of a non-power-series `v`. That part is SKETCH. It is also unnecessary, because (ii) already excludes
  everything.
* Scope.
  * The exclusion is relative to the fixed model `bmod = (1, -1/v, 1, -1, -1)`. A system with a pole of order 2 might be the pull-back
    of a different model. That is outside Conjecture N, so "excluded" is fair.
  * From (ii) with `n = 1`: `v'(0)^4 = -res_0 B`. A rational `v'(0) = a` gives residue `-a^4`. The worker's `-m^20` also needs
    `a = m^5`, so that `tau(0)^(1/4) = 1/m` is rational. This restricts the residues sampled, but by Lemma C (with `phi = c w`)
    the residue is a coordinate artefact anyway.
* The formula `1/tau = Sum 4 beta5_n/(5n+4) u^n`: re-derived. I computed `v` instead by a fixed-point recursion `(n + 4/5) w_n = ..`
  for `w = v/u` (`rf_bsys.build`); `tau` and `v` agree with `prep_*.log` for R0-S5.

## 2. VERIFIED claims, re-derived independently

My construction (`rf_bsys.py`, no `bv_*` import) has these features.
* Own Laurent series over Q.
* `v` by a fixed point, not by the worker's closed form.
* Own gauge solver: particular solution started at component 2 (the worker starts at 0), left null vector, own ODE solver.
* Normalisations imposed intrinsically, `P_11(0,z)` constant or `P_00(0,z)` constant.
* `sigma = C^(-1/2)` from `P_0(0)`. My raw scale differs from the worker's: `C = 2^20` for S5, where the worker has `2^-10`. The
  final `c0` agrees.
* Own pure-1/z solutions and (G) check.
* Hatted prefactor computed as `(u/v)^((1-c_b)/5)`.

`E_M = tau^(a+b+1) Eh_0` takes `Eh_0` from a `fq4_split.Split` of variant `'M0'` (`tau = 1`), not from the `u^0` part of the `'nat'`
base. Graph sums run on `un_engine.Var`, which is shared with the worker; `REFEREE_S4.md` s.4.1 repeats them on a disjoint engine. Run logs are in
the scratch folder and are summarised in `rf_runs.log`.

| worker claim | my check | result |
|---|---|---|
| R0 regression: rows = project rows, recomputed `E_M` = base `E_M` | own rows; `E_M` from `'M0'` | rows identical (g = 2, 3, 4; 2^89-1, 2^61-1); R0 tables identical to the worker's |
| Theorem-U gauge holomorphic, `k <= 12`, S1-S5 | own solver, exact over Q | holomorphic; resonances only at `k = 0, 5, 10` (ODE residue 0), residues `k/5` mod 1 otherwise, as Theorem U predicts; also S6, S7 |
| (G) canonical, `z^-6..z^6` mod `u^6`, no positive `z`-power | own `check_G` | PASS for R0, S1-S7; the stored normalisation fails it, as it must |
| stored `m5` | `normalise(.., 'stored')` | identical: S1 `-146707861988473/35275380525000`, S2 `319129204271/48787200000`, S3, S4, S5 (= S1/1024^2) |
| rows symplectic g = 2, 3, 4, both primes | `un_core.symplectic_check` on my rows | PASS (g = 2, 3, 4, both primes, R0 and S1-S7; `rf_extra.py rows4`) |
| genus 2, 3: every edge graph regular, vertex `c_g T''^(2-2g) + 0`, S1-S5, both primes | own rows, own `E_M`, window `u^0..u^(2g-1)` | 6/6 and 41/41 regular, vertex law holds; tables IDENTICAL to the worker's json in all 7 x 4 (g = 2) and 42 x 6 (g = 3) residues, both primes |
| genus 4, maxout 2 (S1, S2, S3, S5) and S2 full (2^89-1) | S3 maxout 2 (`un_engine`); S2 FULL (disjoint engine) | both IDENTICAL to the worker's tables (11 resp. 379 graphs x 6 coefficients); 10/10 and 378/378 regular |
| controls staleM, dress0 | S1 at 2^89-1 | tables IDENTICAL to the worker's (`tab_S1_g2_89_canon_staleM_k4`, `tab_S1_g3_89_canon_dress0_k6`) |
| stored-normalisation response (S2, genus 4) | S2 at 2^61-1 (s.2.1) | see s.2.1 |
| "6/7 (g = 2), 25/42 (g = 3) graphs differ from R0" | counts from my tables | `REFEREE_S4.md` s.4.5 |

### 2.1 Genus 4 at the second prime and on new systems (VERIFIED, `rf_run.py`, `--kext 0`, window `u^0..u^5`)

| system | prime | run | edge graphs | vertex graph in `T''` | remark |
|---|---|---|---|---|---|
| S3 | 2^89-1 | maxout 2 | 10/10 regular | `(1/180000) T''^-6 + 0 + .. + 0` | table IDENTICAL to the worker's |
| S4 | 2^61-1 | maxout 2 | 10/10 regular | same | new (the worker ran S4 only at genus 2, 3) |
| S6 | 2^61-1 | maxout 2 | 10/10 regular | same, `c0 = 59049/2` | new system, `m = 3` |
| S2 | 2^61-1 | maxout 2 | 10/10 regular | same | second prime |
| S2 | 2^61-1 | maxout 2, STORED (control) | 6/10 singular: graphs 1, 2 from `[u^4]`, graphs 3, 4, 7, 8 at `[u^5]` | `T''^-2`, `T''^-1` terms: law FAILS | the worker's pattern, at the second prime |
| S7 | 2^61-1 | FULL | **378/378 regular** | `(1/180000) T''^-6 + 0 + .. + 0`; total the same | new system, all 379 inner graphs |

Normalisation response of S2 at 2^61-1 (`rf_extra.py norm`; `m5 = 319129204271/48787200000`, as the worker has).
* `(stored - canon)/m5` is non-zero exactly on graphs 0, 1, 2 (from `[u^4]`) and 3, 4, 7, 8 (at `[u^5]`).
* The `[u^4]` column is `-1/8400000, -1/1008000, 1/900000`.
* All seven `[u^5]` rationals of the worker (`563/3402000000`, .., `-17/2700000`) reduce to my residues mod 2^61-1.
* The sum over the graphs vanishes through `[u^5]`.
CONFIRMED at the second prime.

## 3. Lemma C: a change of B is a change of coordinate (PROVED; per-graph identity VERIFIED)

**Setting.** `X = (f, g, B)` is a self-dual cyclic system, `b^X = (f, B/(f^2 g^2), f, g, g)`, with `ord_0 B = -1`. Its model
pull-back is `v_X` (`v_X'^5 = -B v_X`, `v_X(0) = 0`, `a := v_X'(0)`, so `a^4 = -res_0 B`). The quintic has `B_q`, `v_q`
(`v_q'(0) = 1`). Put `phi := v_X^-1 o v_q`. It lies in `w Q[[w]]` with `phi'(0) = 1/a`.

**Lemma C.**
* (a) `X~ := phi^* X` (coefficients `b~ = phi' . b^X o phi`) has the self-dual shape and `prod b~ = B_q`. So it is a univ system
  `f~ = t_c' F~`, `g~ = -(varpi_0/(1-w)) Gm~` with
  `F~ = phi' f(phi)/t_c'`, `Gm~ = phi' g(phi)/g_q`, `F~(0) = f(0)/a`, `Gm~(0) = -g(0)/a`.
  Its model pull-back is `v_X o phi = v_q`, so its `tau` is the quintic `tau`.
* (b) If `P(u, z)` is a Theorem-U gauge of `X`, then `P(phi(w), z)` is one of `X~`, with the same `P_11(0,z)`, `P_00(0,z)` and
  `P_0(0)`. Canonical goes to canonical and stored to stored, with the same `C`, `sigma`, `c0`, `m5`.
* (c) The pure-1/z solutions pull back, so (G) for `X` iff (G) for `X~`, and `T''_X~ = T''_X o phi`.
* (d) `R_X~ = R_X o phi`, `M_X~ = M_X o phi`, `R'_X~ = R'_X o phi`, and every outer-graph amplitude satisfies
  `V^X~_Gamma(w) = V^X_Gamma(phi(w))`.
* (e) Hence a graph `Gamma` is regular for `X` iff it is regular for `X~`, and the vertex law transports.

*Proof.*
* (a) `prod b~ = phi'^5 B(phi) = -phi'^5 v_X'(phi)^5/v_X(phi) = -(v_X o phi)'^5/(v_X o phi) = -v_q'^5/v_q = B_q`. Each `b~_a` is
  `phi' b_a(phi)`, so the shape `(f, ., f, g, g)` is kept.
* (b) Chain rule: `z d/dw P(phi) = phi' (N P - P Nt)(phi) = N~ P(phi) - P(phi) Nt~`. Holomorphy is clear, and the normalisations are
  conditions at the origin, where `phi(0) = 0`.
* (c) `y_(n+1) = Int_0 N y_n`, and `Int_0^phi(w) h(u) du = Int_0^w h(phi) phi' dw`.
* (d) `M = M_0(z v^(-4/5))` depends only on the model coordinate, and `v_X~ = v_X o phi`. `R = Psi^-1 P Psi_mod R_mod(v)` is built
  pointwise from the same data. The amplitude of `Gamma` is a polynomial in the coefficients of `R'`, `M`, the leaves and the TFT at
  the same base point.
* (e) `phi(0) = 0` and `phi'(0) != 0` preserve `ord_0`. `QED`

*The hatted encoding.* With `lam := (w/phi(w))^(1/5)`, the hatted gauge formula gives `[zhat^n]` of the rows of `X~` at `w`
equal to `lam^(1-c_a+4n)` times those of `X` at `phi(w)`. The same holds for `E_M`, since `tau_X~(w)/tau_X(phi(w)) = lam^4`. In
every amplitude the powers of `lam` add up as follows:
* an edge `(m, n)` gives `lam^(4(m+n+1))`, since `c_a + c_(3-a) = 2`;
* a leaf `T(e)` gives `lam^(4e-3)`;
* a vertex gives `lam^-(2g_v-2+n_v)`.
With `Sum psi-exponents = 3g-3-E+L`, the total is `lam^(10(g-1))` for every graph. So `G~(w) = (w/phi)^(2g-2) G(phi(w))`, i.e.
`V~ = V o phi` in the engine's own normalisation.

Status: (a)-(e) are PROVED (elementary; the engine-level statement by the count above). The identity is also VERIFIED numerically.

Conversely, every univ system `(F~, Gm~)` and every `phi` with `phi(0) = 0`, `phi'(0) != 0` give back a B-perturbed system. So
one-pole systems modulo coordinate changes are exactly the quintic-B family, and `B` carries no invariant information. The residue
is a coordinate artefact too (`phi = c w`).

**Check (VERIFIED, `rf_reparam.py`).** I built `phi` and `F~`, `Gm~` over Q and checked `phi'^5 B_X(phi) = B_q` exactly. I then
ran `X~` through the Stage-1 univ pipeline only. That pipeline uses the `ur_gauge` solver with the quintic `tau`,
`un_e4_cases.normalise`, `un_gauge.rows_from_gauge`, and the `'nat'` base with the quintic `E_M`, plus the constant `sigma`. It
shares neither the `tau` recomputation nor the gauge solver of either `bv_sys` or `rf_bsys`. I compared each outer graph with
`(phi/w)^(2-2g) G^X_Gamma(phi(w))`, taken from my direct run of `X`, mod `w^(2g)`.

| X | genus, prime | `phi'(0)` | univ pipeline: `C`, `sigma`, `c0` | `w^(2g-2) V~ = (phi/w)^(2-2g) G(phi)` | `X~` itself |
|---|---|---|---|---|---|
| S1 | 2, 2^61-1 | 1 | 1, 1, 1 | ALL 7 graphs equal (6 non-zero) | regular; `(1/1200) T~''^-2 + 0` |
| S1 | 3, 2^61-1 | 1 | 1, 1, 1 | ALL 42 equal (25 non-zero) | regular; `c_3 T~''^-4 + 0` |
| S4 | 3, 2^61-1 | 1 | 1, 1, 1 | ALL 42 equal | same |
| S5 | 3, 2^61-1 | 1/32 | 1/1024, 32, 1024 | ALL 42 equal | same |
| S6 | 3, 2^61-1 | 1/243 | 1/59049, 243, 59049/2 | ALL 42 equal | same |
| S7 | 3, 2^89-1 | 1 | 16/9, 3/4, -21/20 | ALL 42 equal | same |

So the B-perturbed runs, done with `tau` recomputed, coincide graph by graph with univ-type runs that use only the quintic `tau`.
The univ route itself has no `m`-branch, so this also checks the worker's `m^-(1-c_b)` prefactors for `m = 2, 3`. The sign of `c0`
differs between raw conventions, and the amplitudes do not see it.

Consequences for BVAR.
* s.7 "What is tested ... A failure would have revealed a hidden use of the quintic's B": WEAKENED, as in s.0.
* s.6(ii), "So the graph-by-graph statement genuinely depends on B": correct, but read it as "on the model coordinate". staleM
  uses the `M` of another system. Its failure is the expected non-covariance, of the same kind as the Stage-1 control `'nat_t<k>'`
  of `fq4_split`.
* The status row "B-universality in all genera" is equivalent to `(F, Gm)`-universality in all genera.

## 4. Further attempts to break the outcome (details in `REFEREE_S4.md`)

* **Disjoint engine.** My rows summed on `natural_form/genus4_referee/indep` agree with the `un_engine` runs. The genus-2/3
  multisets of per-graph series are identical for S5, S6, S7, and genus 4 (S6, maxout 2) agrees graph by graph. The genus-4 FULL
  S2 run (378/378 regular) is IDENTICAL to the worker's full table on all 379 graphs.
* **New systems** S6 (`m = 3`, frame `(-1/2, 1/5)`, `f(0) = 2`) and S7 (`f(0) < 0`, a `u^5` term in `B`): regular, vertex law holds.
* **Other real branch** `v -> -v` of the pull-back (S8): the rows and tables are those of the quintic exactly. The branch is a symmetry.
  Disclosure: a float bug in my own patch first produced a spurious failure; it is fixed and documented.
* **Failing controls at 2^61-1** (staleM, staleHat, dress0, dress2, staleboth, nosig) all fail as they should. A CONSISTENTLY wrong
  `tau` (staleboth) keeps the rows symplectic and still breaks regularity.
* **Non-triviality** counts agree with the worker's.

## 5. What survives

| worker claim (BVAR s.8) | referee verdict |
|---|---|
| `ord_0 B = -1` is forced; pole-order perturbations excluded | CONFIRMED (PROVED). The converse paragraph is SKETCH and not needed. The scope is relative to the fixed model |
| `tau`, `v(u)`, `M` recomputed from `B` reproduce the quintic data for R0 | CONFIRMED (VERIFIED independently, g = 2, 3, 4, two primes) |
| Theorem-U gauges for S1-S5 (`k <= 12`); (G) canonical with no positive `z`-power | CONFIRMED (own solver and own (G) check, exact over Q; also S6, S7) |
| B-universality at genus 2, 3, S1-S5, two primes | CONFIRMED, residue by residue (also S6, S7; and through a disjoint engine) |
| B-universality at genus 4: maxout 2 for S1, S2, S3, S5, and in full for S2 (2^89-1) | CONFIRMED for S3 (identical) and S2 full (re-derived on the disjoint engine, `REFEREE_S4.md` s.4.1); S1 and S5 at genus 4 were not re-run. EXTENDED to 2^61-1 (S2, S4, S6; S7 in full; s.2.1) |
| stored normalisation and non-annihilating dressings singular as predicted; staleM | CONFIRMED (also at 2^61-1; plus staleHat, staleboth, nosig) |
| "A failure would have revealed a hidden use of the quintic's B" | WEAKENED: by Lemma C, B is a coordinate. The test checks `(F, Gm)`-universality on pulled-back systems plus the coordinate covariance of the pipeline |
| B-universality in all genera: CONJECTURAL | CONFIRMED as a label, but it is EQUIVALENT to `(F, Gm)`-universality in all genera (Lemma C), so it is not a new conjecture |
| `c0 = 1/lambda` | CORRECTED: `abs(c0) = m^10/abs(f(0) g(0))` |

## 6. Files

All referee files are in this folder. Code: `py -I -B`, stdlib only, each file under 250 lines. No `bv_*` module is imported.

| file | content |
|---|---|
| `REFEREE_S4.md` | part 2 of this report: disjoint engine, new systems, branch test, controls at 2^61-1, counts |
| `rf_ser.py` | own truncated Laurent series over Q; residues mod p |
| `rf_bsys.py` | systems R0, S1-S5 (from BVAR s.2) and S6, S7 (new); `v` by fixed point; own Theorem-U solver; normalisations; (G); hatting |
| `rf_run.py` | rows, `E_M` from the `'M0'` base, graph sum (`un_engine.Var`), vertex law, comparison with the worker's json; controls staleM, staleHat, dress0, dress2 |
| `rf_indep2.py` | the same rows on the module-disjoint engine `natural_form/genus4_referee/indep` (`E_M` rebuilt by that engine) |
| `rf_reparam.py` | Lemma C check through the Stage-1 univ pipeline |
| `rf_extra.py` | S8 branch test, nosig, staleboth, genus-4 normalisation response, genus-4 row checks, non-triviality counts |
| `rf_runs.log` | summary of every referee run (the full logs and json tables stay in the scratch folder `allgenus_s2_bvar_ref`) |

Project and Stage-1 modules are used read-only:
* `fq_conifold` (periods `t_c`, `varpi_0` only), `fq_fp`, `fq3_rec` (printing), `fq4_split.Split` (`'M0'` base);
* `un_engine` (`Var`, `collect`, `get_nat`, `get_base`) and `un_core` (`Arith`, `lift_rows`, `symplectic_check`, `left_dress`);
* in `rf_reparam` only, as the deliberately different route: `un_gauge`, `un_e4_cases`, `ur_gauge`;
* `natural_form/genus4_referee/indep` and `natural_form/referee`, in `rf_indep2`, with in-memory patches only.

Processes: at most 4 of mine at a time (a pool of 3 plus its parent, or three single runs). There was one exception of about a
second: a text edit run with `py` while a pool was busy. No process was left running. CPU: about 1.7 h in total, with genus 4
about 1.4 h of it.
