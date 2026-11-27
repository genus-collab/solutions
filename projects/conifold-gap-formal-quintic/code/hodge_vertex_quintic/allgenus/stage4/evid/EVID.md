# Stage 4, evidence task (falsifier): the connection problem `D = 1` (2026-11-14)

Folder `allgenus/stage4/evid/` (worker "evid" of Stage 4). AI-written and AI-checked only; no human has read it. Labels: PROVED
(argument written here or in the cited project note), VERIFIED (exact computation; field/primes, frames and range stated),
CONJECTURAL, SKETCH, REFUTED. Conventions of `natural_form/theory/UNIFORM_REDUCTION.md` (UR) s.1-s.4, `allgenus/gauge/GAUGE_IDENTITY.md`
(gauge) s.1-s.6 and the docstring of `natural_form/referee/ref_rows.py`. Primes p89 = 2^89-1, p61 = 2^61-1; frames (0,0), (1/3,-2/7).

## 0. Bottom line

* **No failure anywhere.** Nothing is REFUTED. Every test PASSED in every field, frame and window. The deliberate controls failed
  as predicted; their log lines read PASS because they test that the control fails.
* **Highest order reached** (s.4). `D_0` is the constant series of UR s.4, `Q D = P_can Q^m(v(u))`.
  * Exact over Q, frame (0,0): `D_0 == 1 mod z^205`, i.e. (H) at every `n <= 204` and `a_1 = a_3 = ... = a_39 = 0`. So (F_g) is
    VERIFIED (exact) for every genus `g <= 69`. The exact z^200 runs (qrr in both frames, can in frame (0,0)) are listed in s.4.
  * Mod p89 and p61, frame (0,0): `D_0 == 1 mod z^505` (q <= 100), so (F_g) is VERIFIED mod two primes for `g <= 169`. (H) alone: two primes to
    `n = 1000`, exact to `n = 300`.
  * Direct test of the brief, part (1): the referee rows from scratch against the canonical gauge from the periods, in every slot.
    Through `z^35` (`a_5`, `a_7`) in both frames at both primes. Through `z^55` (`a_9`, `a_11`) at p89, frame (0,0).
* **The reduction (s.2), PROVED given Theorem U** and the standard facts listed there. `D = 1` holds iff two things hold.
  * (H): `ph_n(1) = h_n` for every `n` with `5 ∤ n`. This is a parameter-free identity: the WKB polynomial of colour 0 at `q = 0`.
  * (N): `K_q = 0` for every `q >= 1`. Here `K_q := [u^4q zhat^5q] Ghat_11` is a linear functional of the colour-0 WKB polynomials
    (their jets to order `u^4q` at the conifold), the periods `C0`, `C1`, `tau` and the model row `F_2`; and `K_q = f_5q + (lower data)`.

  **Connection form:** fix the free coefficients at the conifold by `K_q = 0`, with no QRR input. Then `D = 1` iff the resulting
  polynomials satisfy `ph^can_n(1) = h_n` for EVERY `n`. One polynomial identity per `n`, as the brief anticipated.
* **What the model fixes.** It is NOT `ph_n(0)`: `ph_n(0) = F_{0,n}` holds in every normalisation (PROVED, s.2.1).
  * `deg ph_n = floor(4n/5)` (VERIFIED `n <= 1000`). So the free coefficient `f_5m = [u^4m] ph_5m` is the TOP coefficient.
  * Equivalently, `f_5m` is the `u^0` term at fixed `z` of the WKB function at the conifold.
  * The conifold side fixes `f_5m` through `K_q`, which involves the gauge data. The naive guess "`f_5m` = the model's value 0" is
    false: `f_5 = 4537/3937500`.
* **Validation against E4 and `m5`** (s.3, s.5).
  * `K_1 = .. = K_4 = 0`, which reproduces `a_1 = a_3 = 0`.
  * In the window `z^35`, `u^28`, the reduced formula equals the direct E4 machinery in all 25 entries, slot by slot.
  * The stored normalisation gives `K_1 = m5 = -43397/807975000` exactly (control S). The regularity checks cannot see it.

## 1. Part (1): E4 extended directly (`ev_direct.py`)

Same construction as `gauge/gi_actual_norm.py`, whose functions are imported read-only.
* Actual side: `ref_hat.HatRows` (QRR normalisation, all colours, own checks `h_auto`, `psi_orth`, `sympl`, `covar`, `circ`), the symbol
  `M_0` from the same rows, and `Ghat^act = Qhat M_0(tau zhat)^-1 Psi0^-1`.
* Canonical side: `ur_gauge.gauge` over Q from `(varpi_0, t_c, tau)`, normalised by `P_11(0,z) = 1`, then hatted.
* Measured cost (probe, then extrapolated). Rows: 2 s at `(z^24, u^17)`, 5 s at `(z^30, u^23)`. Gauge over Q: 6.7 s at `k = 20`,
  13 s at `k = 26`; 42 s at `k = 24` with `NU = 72` (E4 log). Predicted about 1 min per frame at `(z^35, u^28)`; observed 70 s.

| window | frame | primes | (0) graded + valuations | (1) slots = canonical | stored differs | normalisation slots | time |
|---|---|---|---|---|---|---|---|
| `zhat <= 35, u < 29` | (0,0) | p89, p61 | PASS | 2700/2700 (2598 non-zero) | 1939 slots, first `(5,0,0,4)` | `a_1..a_7 = 0` | 70 s |
| `zhat <= 35, u < 29` | (1/3,-2/7) | p89, p61 | PASS | 2700/2700 (2602 non-zero) | 1943 slots, first `(5,0,0,4)` | `a_1..a_7 = 0` | 73 s |
| `zhat <= 55, u < 45` | (0,0) | p89 | PASS | 6440/6440 (6282 non-zero) | 5231 slots, first `(5,0,0,4)` | `a_1..a_11 = 0` | 259 s |

* (H) holds in the window: `h_auto` passes for every `5 ∤ n <= NZ`.
* The odd `a_j` are the independent tests; the even ones are forced by symplecticity and come out 0 as well. The exact stored-gauge
  values `log P^s_11(0,z)` are printed for contrast (`a_1 = m5`; `a_2 = -18223115111903/31086838125000000`, etc.).
* Status: VERIFIED. The gauge's `c^act` equals 1 through `z^55` (p89) and `z^35` (two primes, two frames).

## 2. The reduced criterion

### 2.1 The colour-0 WKB polynomials (ref_rows docstring; `ev_core.Wkb`)

The hatted WKB function of colour `i` is `Phi_i = Sum_n zhat^n ph_n(u)`. For colour 0 everything is over Q. Write
`Phi_j = A_{-j} Phi_{j-1}`, `A_e = 1 + zhat (1-u)((4/5) th_zhat - e/5 - th_u)`. Then the `zhat^(n+1)` coefficient of (HPF) reads
`(1-u)(4n - 5 th_u) ph_n = -R_n(ph_0..ph_(n-1))`. With `rho = -R_n/(1-u)` this gives `(4n - 5k) [u^k] ph_n = rho_k`.

* **(a) PROVED.** `ph_n(0) = F_{0,n}`, the model's homogeneous coefficient, in every normalisation. The `u^0` part of (HPF) is the
  model's limit equation, the pivot `4n` never vanishes, and the free coefficient (`k = 4n/5`, `5 | n`) is never at `k = 0`.
  * So "ph_n at the conifold equals the model's" holds automatically and cannot be the criterion.
  * VERIFIED: S2, `n <= 1000` at two primes and `n <= 300` exact; V5 (all rows, all colours, `n <= 35`).
* **(b) PROVED: `ph_n` is a polynomial.** `R_n(1) = 0` always: every word of `A_{-5}..A_{-1} - 1`, and every other term of (HPF), has a
  LEFT factor `(1-u)`.
  * The resonance (`rho_{4n/5} = 0` for `5 | n`) is unobstructed given Theorem U: the lower data extend to the regular-type solution
    `P_can Q^m c`, so a power-series `ph_n` exists.
* **(b') VERIFIED, not proved.** These hold for `n <= 1000` at two primes and `n <= 300` exactly (`ev_struct`):
  * `deg ph_n = floor(4n/5)` for `n >= 3`;
  * `ph_1 = -3/20` and `ph_2 = 9/800`;
  * `deg rho_n < 4n/5` when `5 | n`, so every resonance is unobstructed for a trivial reason;
  * for `n = 5m` the free coefficient `f_5m` is the TOP coefficient of `ph_5m`.

  At fixed `z` the term `z^n u^(-4n/5) ph_n(u)` therefore has exponents in `(-4n/5, 0]`. The `u^0` part of `Phi` at fixed `z` is
  `1 + Sum_m f_5m z^(5m)`: the regular part at the conifold, which is also the termwise limit `u -> oo` (the orbifold point).
* **(c) The two normalisations.**
  * QRR: `f_5m := h_5m - (ph_5m(1) - f_5m)`, with `h_n = [z^n] h_{w_0}`, `w_0 = {1 - zeta^j} U {-5}`, exact via Newton sums for the
    roots `y_j = 1/(1 - zeta^j)` of `5y^4 - 10y^3 + 10y^2 - 5y + 1`.
  * Canonical: s.2.3.

### 2.2 Lemma P: the colour-reduced pairing (PROVED; VERIFIED V6)

Let `R^(t) = Psi0 M_0(t)` be the model's hatted WKB matrix (rows `b`, colours `i`).
* Symplecticity `M_0(t)^-1 = M_0(-t)^T` and `Psi0^-1 = pinv` with `pinv^T pinv = (1/5) eta'` (`eta'_bb' = [b+b' = 3 mod 5]`;
  `eps_b eps_b' = 1` on those pairs) give `R^(t)^-1 = (1/5) R^(-t)^T eta'`.
* With Prop. S (`R^_b^(i)(t) = zeta^(i(b+1)) F_b(zeta^-i t)`), the covariance `Qh^(i)_{a,n} = zeta^(i(a+1-n)) Qh^(0)_{a,n}` and
  `Sum_i zeta^(im) = 5[5 | m]`, the following holds for every covariant `Qhat`:

      [zhat^N] (Qhat R^(tau zhat)^-1)_ac = [N = a - c mod 5] Sum_{n<=N} Qh^(0)_{a,n}(u) (-tau)^(N-n) F_{3-c, N-n}.

* Only colour 0 enters, so the whole actual gauge is computable exactly over Q (`ev_core.pair`).

### 2.3 Lemma N: reading the normalisation (PROVED given Theorem U)

**Statement.** Take any choice of the free coefficients. Let `Phi` be the resulting polynomial-type colour-0 solution, `Qhat` its hatted
rows (all colours by covariance) and `Ghat := Qhat R^(tau zhat)^-1`. Then `Ghat` is the hatted form of `P_can c(z)` for a unique
`c in 1 + z^5 Q[[z^5]]`. Moreover `[u^j zhat^5q] Ghat_11 = 0` for `j < 4q`, and `K_q := [u^4q zhat^5q] Ghat_11 = [z^5q] c`
`= [u^4q zhat^5q] Ghat_22`.

*Proof.*
1. By Theorem U, `P_can Q^m` is a formal WKB solution of the quintic with the same exponentials (choice of `v(u)`).
2. In hatted form it is `Ghat^can R^(tau zhat)`: a power series in `u` at fixed `zhat`.
3. Its colour-0 row 0 is `C0^-1 Phi'`. (HPF) forces each `ph'_n` to solve `(4n-5k) c_k = rho_k`, so it is the same kind of object as
   `ph_n`.
4. Two WKB solutions with the same exponential differ by a constant scalar series `d(z)` (formal WKB uniqueness; the `u_i` are distinct).
5. Both are power series in `u` at fixed `zhat`, and `z^j = zhat^j u^(4j/5)`. Hence `d in Q[[z^5]]`, the same for all colours
   (covariance fixes `z^5`). Put `c = d`.
6. (R1) for `P_can` gives `P_can,11(0,z) = P_can,22(0,z) = 1` (gauge s.3). The hatted `(1,1)` entry at `zhat^5q` is
   `tau^(1/2) u^4q p_{5q,1}(u)`, and the claims follow. QED

**Slope.** `dK_q/df_5q = F_2(0)/(C1(0) C0(0)) = (-1)/((-1)(1)) = 1`, printed in every log.

### 2.4 Theorem (reduced criterion; PROVED given Theorem U and the u = 1 convention of ref_rows)

**Part 1: (H) iff `h_auto` for all `5 ∤ n`.** Let `Phi^pol` carry the QRR rule. The general solution is `Phi^pol d(z)`, whose hatted
value at `u = 1` (`s = 1`, `z = zhat`) is `Phi^pol(zhat,1) d(zhat)`. QRR demands that this equal `h(zhat)`.
* If `ph_n(1) = h_n` for every `5 ∤ n`, then `d = 1`, and the QRR rows are `Phi^pol`, a power series in `u`.
* Otherwise `d` has a `z^n` term with `5 ∤ n` (the first failing `n`), so (H) fails.

**Part 2: (N) given (H).** `Q = P_can c Q^m` by Lemma N, so `D = c^-1`, and `D = 1` iff `K_q = 0` for all `q`. Since `K_q` is affine
in `f_5q` with slope 1, `K_q = 0` determines `f^can_5q`.

**Connection form.** `D = 1` iff `Phi^can = Phi^QRR`, iff `ph^can_n(1) = h_n` for all `n`. By Theorem G (gauge s.3), `c(z)` is also the
factor in `P^act Y^C = c(z) Y^C_Q`, so this is the vanishing-solution form (N) of the brief, read off one entry.

**Genus bookkeeping.**
* `D_0 == 1 mod z^(N+1)` iff (H) holds for `5 ∤ n <= N` and `K_q = 0` for `5q <= N`.
* (F_g) needs `N = 3g - 3`.
* `a_j` (`log c = Sum a_j z^5j`, odd `j` independent) first matters at genus `ceil((5j+3)/3)`: `a_5` at 10, `a_7` at 13, `a_19` at 33 (`z^95`).

## 3. Validation of the reduced machinery (`ev_validate.py`; window `zhat <= 35`, `u < 29`; p89 and p61; both frames)

All PASS (`logs/ev_validate_*.log`, 53 s per frame). `ev_core` shares no code with the project; the comparisons are with the referee's
from-scratch modules (module-disjoint from `fq_*`) and with the direct E4 machinery.

| check | content | result |
|---|---|---|
| V1 | `C0 = varpi_0`, `C1 = D t_c` from my dual-number periods = `ref_hat.HatRows(..).C0, .C1` | equal, `u < 29` |
| V2 | `h_n` by Newton power sums (exact) = `ref_rows.hcoef` (colour 0) | equal, `n <= 39` |
| V3 | incremental `ph_n` (QRR rule) = `ref_rows.wkb` | equal, `n <= 39` |
| V4 | colour-0 hatted rows `Qh_a` = `HatRows.Q[a][0]` | all `a`, `n <= 35`, `u < 29` |
| V5 | model rows `F_b` (model recursion only) = quintic symbol `Psi0 M_0`, i.e. `Qhat(u=0)` (Prop. S) | all `b`, all 5 colours, `n <= 35` |
| V6 | Lemma P: colour-reduced pairing = `gi_actual_norm.actual`'s `Ghat^act` (five colours, matrix inverse) | 25 entries x 36 orders x 29 `u`-terms (2598 / 2602 non-zero) |
| V7 | the exact-Q computation reduced mod p = the mod-p computation | `ph_n`, rows, `F_b` |

## 4. Results of the reduced test (`ev_reduced.py`, `ev_struct.py`; all ALL CHECKS PASSED)

| field | frame | mode | (N) or connection identity | (H): `ph_n(1) = h_n`, `5 ∤ n` | REG + (R1) | time |
|---|---|---|---|---|---|---|
| Q | (0,0) | qrr | `K_q = 0`, `q <= 40` (`z^200`) | `n <= 204` (164 values) | `q <= 12` | 131 s |
| Q | (0,0) | can | `ph^can_5q(1) = h_5q`, `q <= 40` | `n <= 204` | `q <= 12` | 151 s |
| Q | (1/3,-2/7) | qrr | `K_q = 0`, `q <= 40` | `n <= 204` | `q <= 12` | 135 s |
| Q | (0,0) and (1/3,-2/7) | qrr, can | `q <= 20` (`z^100`) | `n <= 104` | `q <= 20` | 45 s each |
| p89, p61 | (0,0) | qrr | `K_q = 0`, `q <= 100` (`z^500`) | `n <= 504` | `q <= 20` | 18 s each |
| p89, p61 | (0,0) | qrr | `K_q = 0`, `q <= 40` | `n <= 204` | `q <= 40` | 8-9 s |
| p89 | (7/11,-3/5) | qrr | `K_q = 0`, `q <= 40` | `n <= 204` | `q <= 40` | 8 s |
| p89, p61 / Q | - | (H) alone (`ev_struct`) | - | `n <= 1000` / `n <= 300` | - | 13 s / 16 s |

* In every run the slope `dK_q/df_5q = 1`, every resonance is unobstructed, and `[u^4q] Ghat_22 = [u^4q] Ghat_11`.
* In mode `can`, `f^can_5q` agrees with the frame-independent `f^QRR_5q` in both frames. Sample values: `f_5 = 4537/3937500`,
  `f_10 = -177640111151/31007812500000`, `f_15 = -284328315069009259339/6226756347656250000000` (more in `logs/ev_struct_Q_300.log`).
* Status: VERIFIED. Exact over Q, frame (0,0): `D_0 == 1 mod z^205` (`g <= 69`). Mod p89 and p61: `mod z^505` (`g <= 169`).
* Not run: `z^1000` mod p. Estimate: about 3 min per prime. Launching it as an extra background job was declined by the session's
  permission check, and it was not needed for the brief.

## 5. Failing controls (each must fail the test it targets; all did)

| control | where | prediction | outcome |
|---|---|---|---|
| (C) stored gauge `P^s = P_can(1 + m5 z^5 + ..)` | `ev_direct` | differs from the actual gauge | 1939 / 1943 / 5231 slots, first at `zhat^5 u^4`, entry (0,0) |
| S: `f_5 -> f_5 + m5` (stored normalisation of the rows) | `ev_reduced .. 1`, Q and p89 | `K_1 = m5` exactly; regularity blind | `K_1 = -43397/807975000`; REG passes |
| T: `tau -> 1` (`v = u`, wrong natural variable) | same | regularity fails | fails at 2 of 2 orders |
| W: QRR weights without `-5 zeta^i` | same | (H) fails at `n = 1` | fails at `n = 1, 2, 3, 4` |
| B1: `B_96 -> B_96 + 1` inside `h` | `ev_ctrl.py`, Q and p89 | (H) unchanged; `K_q = 0` for `q <= 18`; `K_19 = -S_95/(96*95)` | exactly as predicted; `K_19 = 1/(9120 * 5^95)` |
| B2: `B_8 -> B_8 + 1` | same | (H) fails first at `n = 7` | fails at `n = 7, 8, 9` |

* B1 shows that the identity at `z^95` (`a_19`, genus 33) genuinely involves `B_96`, and that nothing at lower order or in (H) sees it.
  The high-order test is therefore not vacuous.
* The B1 value `K_19 = 5^-95/9120` uses `p_95 = 0` (s.8).

## 6. Status and caveats

* **PROVED here:**
  * s.2.1(a) and (b);
  * Lemma P;
  * Lemma N, given Theorem U;
  * the theorem of s.2.4, given Theorem U and the `u = 1` (q = 0) normalisation convention of `ref_rows.py` (the QRR constants
    `h_{w_i}`, used throughout the project and not re-derived from LP here);
  * `p_{5k} = 0` for odd `k` (s.8).
* **Quoted, PROVED elsewhere in the project (AI-audited only):** Theorem U (UR s.3); (R1) and Theorem G (gauge s.3); symplecticity and
  Prop. S of the model (UR s.2); covariance. Standard: uniqueness of formal WKB solutions up to a constant scalar series.
* **CONJECTURAL:** `D = 1` to all orders, i.e. (H) and (N) for all `n`. **VERIFIED** in the ranges of s.1 and s.4.
* **Not independent of the inputs, by design.** The direct and reduced tests share the PF operator, the QRR constants and the model.
  They differ in everything else (rows, gauge, all colours versus colour 0, mod p versus Q). The regularity checks are automatic
  consequences of Theorem U, so they test the conventions and the code, not (N).
* **Finite checks are evidence, not proof.** In particular the degree pattern `deg ph_n = floor(4n/5)` is VERIFIED only.

## 7. Files (all new; `py -I -B`, stdlib only; run with the working directory in the scratch folder; logs use relative paths)

* **`ev_core.py`** (exact colour-0 machinery, generic field). Imports nothing from the project.
* **`ev_reduced.py`** (the reduced criterion, modes `qrr`/`can`, controls S/T/W), **`ev_struct.py`** (s.2.1 facts, (H) alone) and
  **`ev_ctrl.py`** (B1/B2). These import `ev_core` only.
* **`ev_direct.py`** (part 1). Read-only imports: `gauge/gi_core.py` (`gauge_frame`, `p11_at_zero`, `inv_poly`, `normalise`),
  `gauge/gi_actual_norm.py` (`actual`, `hatted_can`, `mu`), and through them `natural_form/theory/ur_gauge.py`, `fq_conifold`,
  `fq_series`, `fq_fp` and `natural_form/referee/ref_hat.py`, `ref_rows.py`, `ref_fp.py`.
* **`ev_validate.py`**. Read-only imports: `gi_actual_norm`, `ref_rows` (`wkb`, `hcoef`), `ref_hat` (`HatRows`, `M0_of`), `fq_fp`.
* **`logs/`**: 24 logs, each under 250 lines, each ending in `ALL CHECKS PASSED`.
* No existing file was modified. No `__pycache__` was created. Every process I started has finished. The maximum was 4 concurrent
  processes.
* `ur_model_check.py` was deliberately NOT imported: it overwrites its own log at import time.

## 8. Observations for the proof of `D = 1` (SKETCH unless labelled)

1. **No analytic connection problem in hatted variables.** Both halves of (F) are identities between a polynomial fixed at the
   conifold and its value at `u = 1`.
   * (H): `ph_n(1) = h_n` for `5 ∤ n`, with no parameter.
   * (N): `ph^can_5m(1) = h_5m`, with `f^can` fixed by the linear functional `K_q`.

   A closed formula for `Phi^can(zhat, 1)` would settle both. A candidate is an integral representation valid for all `u` (UR s.8)
   evaluated at `q = 0` by stationary phase, i.e. Stirling.
2. **PROVED: `p_m := Sum_{j=1..4} (1 - zeta^j)^-m = 0` for every odd multiple `m` of 5.** With `y_j = 1/(1 - zeta^j)` one has
   `y_j^5 = (y_j - 1)^5`, so `y_j^5 = -(1 - y_j)^5 = -y_{-j}^5`, and the pairs `j, -j` cancel in `y^(5k)` for odd `k`.
   VERIFIED for `m <= 200`.
   * Hence at the independent orders `z^5k` (`k` odd) the QRR exponent is `[z^5k] log h_{w_0} = B_{5k+1}/(5k(5k+1) 5^(5k))`: only the
     twisting weight `-5` contributes.
   * At even `k`, `log h` has no `z^5k` term at all (only odd powers occur), matching the symplectic forcing of the conifold side.
   * PROVED from Lemma N: `Phi^QRR = c Phi^can`, so `h = c Phi^can(., 1)` at `u = 1`. Given (H) this yields
     `a_k = B_{5k+1}/(5k(5k+1) 5^(5k)) - [z^5k] log Phi^can(z, 1)` for odd `k`.
   * So (N) asks the conifold-normalised WKB function at `q = 0` to reproduce, at each independent order, exactly the Stirling
     coefficient of the twisting weight.
3. **Three points, not two.** The free coefficients `f_5m` are the top coefficients, i.e. the termwise limit `u -> oo` at fixed `z`
   (`q = oo`). So (N) ties together `q = 0` (QRR), the conifold (model plus gauge) and the `u^0` part at fixed `z`.
   * The degree bound `deg ph_n <= 4n/5` (VERIFIED `n <= 1000`) is equivalent to the absence of positive `u`-powers at fixed `z`.
   * A naive weight count gives only `deg ph_n <= floor((4n+4)/5)`, so extra top-weight cancellations occur. They are unexplained.
