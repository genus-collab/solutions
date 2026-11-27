# Referee cyc_adversary, part 2: frames, deformations, code, corrections, compliance (2026-11-14)

Continuation of `REFEREE.md`: s.0 verdict, s.1 notation, s.2 task (a), s.3 task (b). Same labels and notation.

## 4. Task (c): does "every conifold frame" cover every frame the programme uses?

**The programme's frames.** Every code I inspected builds the same family. It sets `C0 = varpi_0 = f_0 + alpha g + beta f_2`,
`y_1 = g`, `t_c = g/C0` and `C1 = D t_c`. The codes are:
* `fq_conifold` and `gi_core` (a wrapper of it);
* `ur_gauge` (frame (0,0)) and `pp_gauge` (frame (0,0));
* the poly referee's `rf_core` and the gamma referee's `rf_res`;
* `cr4_core` (it also has a probe `dg`, a deliberately wrong `g`), `ev_core`, and `pb_core.set_frame`.

The `(alpha, beta)` values that occur, with where they occur:

| (alpha, beta) | where |
|---|---|
| (0,0) | everywhere |
| (1/3,-2/7) | UR s.6.2, GI s.5-s.7, poly and evid referees, pb |
| (7/11,-3/5) | gamma referee, Stage-4 critic, pb_cyc |
| (1,0), (0,1), (2/3,-1/7) | FORMAL_QUINTIC_genus2 s.4 |
| alpha in {0, 1/3, -2} | gamma referee s.3.6 |

GI s.5 has no frame of its own: it uses (0,0) and (1/3,-2/7). All of these frames have `varpi_0(0) = 1` and `g = u + O(u^2)`, so all
are CYC frames.

**VERIFIED**, the operator identity `M = PF4` coefficientwise, exact over Q, by `cy_frames` (precision `u^74` at `N = 80`):
* the frames (0,0), (1/3,-2/7), (7/11,-3/5), (1,0), (0,1), (2/3,-1/7), (-2,0);
* the closing relation and the PF5 property on actual objects, in (0,0), (1/3,-2/7), (7/11,-3/5) (`cy_rows`).

pb_cyc had checked the programme's frames (0,0), (1/3,-2/7), (7/11,-3/5) through `u^56`.

**Label clash (K1).** `cr5_cyc_run` uses `y_0 = f_0 + a f_1 + b f_2`; CYC s.8 row (e) states this.
* Its "(1/3,-2/7)" is `f_0 + f_1/3 - 2f_2/7`, with u^2 coefficient `-2/7`.
* The programme's frame of that name is `f_0 + g/3 - 2f_2/7 = f_0 + f_1/3 - (11/210) f_2`.
* These are different frames (`cy_frames` (E0)). The same holds for "(7/11,-3/5)".
* This is harmless: CYC covers both families, and both pass in my run.

**Frames in the other sense.** FORMAL_QUINTIC_genus2 s.3 calls the R-matrix normalisations F0-F5 (Hodge classes, `rows_times(h)`)
"frames". Each one multiplies every row by a z-series constant in q, possibly depending on the colour. That commutes with `zD` and
with the `c_a`, so the cyclic system, including the closing relation, is unchanged.

**Scope (K2).** CYC's proof never uses `y_0(0) != 0`. In `K = Q((u))` every non-zero element is a unit.
* The same proof gives `M(y_0, y_1) = PF4` in more cases:
  * for every non-zero power-series solution `y_0` not proportional to `g` (so `Dt != 0`), with `y_1 = g`;
  * for every pair `(y_0, y_1)` in `S cap K` that spans a B-Lagrangian plane.
* B has rank 2 on `S cap K = span(f_0, f_1, f_2)`, with radical `Q g` (C4; `B(f_0, f_1) = -7/5`, `B(f_0, f_2) = 2`). So these
  Lagrangian planes are exactly the planes that contain `g`. VERIFIED to `u^68`:
  * `y_0 = f_1`, `f_2`, `f_1 - 2 f_2`, each with `y_1 = g`;
  * the pairs `(f_1, f_2)` and `(f_2, f_1)`.
* What `y_0(0) != 0` buys is downstream.
  * (i) Theorem U (UR s.3) needs the unit shape `val(c) = (0, -1, 0, 0, 0)`, and that shape holds iff `y_0(0) != 0`. For example,
    `val(c) = (0,-3,0,1,1)` for `y_0 = f_1` and `(-2,-1,-2,2,2)` for `y_0 = f_2`.
  * (ii) The normalisations `C0(0) = 1`, `t_c'(0) = 1` and `t_c = u + O(u^2)` hold in the fq_conifold family only, not in every CYC
    frame (e.g. not in `(2 f_0, 3 g)`). They are used in PB s.2.8 ("in every frame"), POLY s.3 (`b_0(0) = t_c'(0) = 1`), GI
    Lemma Q(d) and P(b) step 7 (the comparison at `u^1`).
* The programme uses only the fq_conifold family, so nothing breaks; but the two statements should be kept apart.

**The unit shape does not see the frame (K3).**
* The non-isotropic pair `(f_0, f_1)` and the wrong-gamma pairs `y_1 = f_1 + f_2/2` and `f_1 + (701/1000) f_2` all give
  coefficients of exactly the Theorem-U shape. Theorem U would accept them.
* But the system is then not the quintic system (`cy_rows`):
  * the closing relation of the WKB rows fails from `z^2` on;
  * the vanishing solution fails PF5 at `n = 6`.
* So the identification of the Theorem-U system with the quintic system rests on `gamma = 7/10` exactly. That is CYC s.1 together
  with C4, and it is PROVED.

## 5. Task (d): deformations, and which lemma fails

| deformation | outcome (VERIFIED unless stated) | lemma that fails |
|---|---|---|
| `PF4 + (1/7)(1-u) D` (not self-adjoint) | `PF4* != PF4`. Every `M(y_0, y_1)` is formally self-adjoint for ARBITRARY units `y_0`, `Dt` (it is palindromic; checked on random non-solutions, `cy_run` (B4)), so no frame can work | C1; then C2(ii) |
| `A_2 + u^2/3` alone (not self-adjoint) | The Frobenius structure survives (`gamma = 7/10`). `B(f_0, f_1)` is not constant. CYC fails in the operator's own frame and in the old one | C1, C2(ii) |
| `A_0 + (1/7)(1-u)` (self-adjoint) | CYC holds for the perturbed operator in its own Frobenius frame (`gamma = 7/10`); it fails with the unperturbed periods | none. CYC is structural; the periods must belong to the operator |
| `A_2 + u^2/3`, `A_1 + D(u^2)/3` (self-adjoint) | as in the previous row | none |
| `A_4 = u + u^2/4`, `A_3 = 2 D A_4`, `A_1 = D A_2 - D^3 A_4` (self-adjoint) | Exponents unchanged; `gamma = 19/20`. CYC holds with middle factor `A_4 y_0^2 Dt^2` and fails with the literal `u` | none. CYC's `u` is `A_4 = L^-5`, which makes `prod c_a = L^5` match the symbol of PF5 |
| `A_3 + 1/2` (exponents `0, 1, 3/2, 2`) | Not self-adjoint, by the criterion of REFEREE s.3. The obstruction map has rank 3 on `span(f_0, f_1, f_2)` (VERIFIED), so no solution `g log u + h` lies in R. Argument: the top log-coefficient of a solution in R is a power-series solution, so `S cap R` is the 3-dimensional power-series space, and an alternating form on it is degenerate | s.1 (Frobenius structure), C1, C2(iv) |
| `A_4(0) != 0`, e.g. `PF4 + D^4/10` (self-adjoint) | Argument, not run: `u = 0` is an ordinary point (the pivot is `(1/10) n(n-1)(n-2)(n-3)`). There are four power-series solutions and no log solution, so `tau` acts trivially, C4 says nothing, there is no `g`, and "conifold frame" is undefined | s.1 |
| non-isotropic pairs `(f_0, f_1)`, `(f_0, f_2)`, `(f_0 + f_1, f_2)` | CYC fails. `(f_0, f_1)` has the Theorem-U shape | the hypothesis of C5: the span is not Lagrangian |
| wrong `gamma'`: `y_1 = f_1 + f_2/2`, `f_1 + (701/1000) f_2` | CYC fails. The WKB closing relation fails from `z^2`, and `Y^C_Q` fails PF5 at `n = 6`. The Theorem-U shape is intact | C4 does not apply, since `y_1` is not the log coefficient; `B(f_0, y_1) = -7/5 + 2 gamma' != 0` |
| `varpi_0(0) = 0`: `y_0 = f_1, f_2, f_1 - 2 f_2`; pairs `(f_1, f_2)`, `(f_2, f_1)` | CYC HOLDS (to `u^68`); the `c_a` do not have the Theorem-U shape | no lemma of CYC; Theorem U's hypotheses fail |

**Remark (SKETCH; VERIFIED on the three self-adjoint deformations above).** CYC's proof uses about PF4 only three things: (i)
self-adjointness, (ii) `A_k` in `Q[[u]]` with `A_4 = u` times a unit, (iii) `det Bmat != 0`.
* Then `A_3(0) = -2 A_4'(0)`, so the indicial polynomial is `A_4'(0) n(n-1)^2(n-2)` and the exponents are `{0, 1, 1, 2}`.
* The obstruction map always has a one-dimensional kernel, because `P_-3'(0) != 0` and `P_-3'(2) != 0`.
* So CYC, with middle factor `A_4 y_0^2 Dt^2`, holds for every such operator. Point (iii) is the Wronskian lemma for a basis of
  solutions; I did not prove it in general.
* Consequence: CYC is not a numerical accident of the quintic coefficients, and no finite coincidence is hidden in it.

## 6. Own code and logs (`py -I -B`, working directory the session scratch folder, exact over Q unless stated)

| file | content | recorded runs (`logs/`) | result |
|---|---|---|---|
| `cy_core.py` | truncated Laurent series over Q; `K[log u]`; operators (compose, adjoint, apply) | imported | - |
| `cy_frob.py` | PF4 (perturbable); closed form `P_j(x)`; Frobenius basis and log obstruction by recurrence | imported | - |
| `cy_log.py` | check bookkeeping (a control MUST fail) | - | - |
| `cy_run.py N` | (A) Frobenius basis and `gamma` (found); (B) C1, own criterion, M always self-adjoint; (C) Bmat, C4, obstruction = `B(f_0, .)`; (D) C3; (G) PF5 identity, HPF u-form | `24`, `40`, `64` | ALL PASS (precision `u^59..u^62`); controls B3, B5, C5 fail |
| `cy_frames.py N` | (E) `M = PF4` coefficientwise: 17 frames and 5 controls, plus the Theorem-U shape; (F) deformations | `24`, `48`, `80` (13 s) | ALL PASS (`u^74`, or `u^68` when `y_0(0) = 0`); every control fails |
| `cy_rows.py N NZ` | (G) cyclic system => PF5 on `Y^C_Q`, `n <= 16`; (H) WKB rows of PF5 => closing relation | `30 5`, `30 6`, `44 9`, `64 12` (9 s) | ALL PASS ((G) to `u^57`; (H) `z^0..z^12` to `u^52`); 8 controls fail |
| `cy_proj.py NU` | (P1) project periods = own; (P2) `Cf(kinv[a]) = 1/c_a` mod 2^89-1; (P3) UR's `b` | `20`, `40` | ALL PASS; the factor-5 control fails |
| `cy_extra.py N` | by-product: `Psi = varpi_0 (Y^C_Q)_0` is the same in five frames, and scales with `g` (`Psi -> 3 Psi` for `3g`), `n <= 16`. This agrees with PB s.2.7, `V = C0 (Y^C_Q)_0` "in every frame" | `40` | ALL PASS; the wrong-gamma control fails |

Methods differ from the critic's and pb's in three ways:
* The Frobenius basis comes from the closed form `PF4 u^x = Sum_j P_j(x) u^(x+j)`. The log solution and `gamma` come from `d/dx`,
  not from applying the operator to series.
* The operator identity is compared coefficient by coefficient, not on solutions (critic) or on monomials (pb).
* Both directions of the Corollary are tested on actual objects: WKB rows, and `Y^C_Q`. Neither earlier check did this.

The critic's `cr5_cyc_run.py 24` was re-run from a scratch copy (together with `cr5_cyc.py`); its log is identical to the stored one,
timings stripped.

## 7. Precise list of corrections (none changes a result)

* **K1** (label).
  * Where: CYC s.8, rows (d) and (e); `cr5_cyc_run`; SYN s.1 table, row (CYC), "3 frames"; SYN s.4.1 (e).
  * The problem: the critic's frames "(1/3,-2/7)" and "(7/11,-3/5)" are `f_0 + a f_1 + b f_2`, not the programme's `(alpha, beta)`
    frames.
  * Correction: rename them. The programme's frames are VERIFIED by pb_cyc (to `u^56`) and by `cy_frames` (to `u^74`).
* **K2** (scope).
  * Where: CYC s.0, "any power-series solution with `y_0(0) != 0`"; SYN2 s.9.
  * The theorem needs only `y_0 != 0` not proportional to `g`, or more generally a Lagrangian pair in `S cap K`.
  * `y_0(0) != 0` is the condition for the Theorem-U shape of the `c_a`.
  * PB s.2.8 ("`t_c'(0) = varpi_0(0) = 1` in every frame") should read "in every frame `(alpha, beta)` of fq_conifold". So should the
    normalisations of POLY s.3 and of P(b) step 7.
* **K3** (robustness; not an error). The Theorem-U shape does not detect a wrong `gamma` or a non-isotropic pair (s.4).
  Recommendation: state in the write-up that `gamma = 7/10` is the unique `gamma` with `B(f_0, f_1 + gamma f_2) = 0`, and the unique
  unobstructed one; the proof is the two-line closed form of REFEREE s.3.
* **K4** (statement).
  * The Corollary should say where the `Y_a` live: any differential extension of `(K, D)` in which the `c_a` are invertible. In the
    programme these are three rings.
    * For the WKB rows, `Q(zeta)((s))((z))`. The conjugated form `zD + zeta^i L` is the image of the identity in `K[z][D]` under the
      ring automorphism `D -> D + zeta^i L/z`.
    * For `P_can Q^m` and `Y^C_Q`, `L((u))` with `L = Q((z))`.
    * The hatted `Q[[u]][[w]]`, via `z = s^4 w`.
  * Remark: CYC s.0's "the closing relation of the quintic rows is a theorem" presupposes the rows, i.e. WKB solutions of PF5 with
    Laurent-polynomial coefficients. At all orders these exist by P(a) and POLY Lemma 2, which give a polynomial solution of HPF,
    hence such a psi. This is a separate input, not CYC's.
* **K5** (missing use).
  * fq_ring eliminates `i_3` with `B(c_0, c_0 Lam + i_0) = 0`.
  * Identities proved in N modulo this relation (e.g. the exact symplecticity of FORMAL_QUINTIC_genus2 s.1) reach the conifold
    frame only because `Cf` kills it, i.e. because `B(varpi_0, g) = 0`. That is C4.
  * C4 thereby PROVES FORMAL_QUINTIC_genus2 s.1 "Cf preserves the Wronskian relation" and fq_conifold's `pairing_vacuous`.
  * Correction: list C4 next to (CYC) in the input table.
* **K6** (wording).
  * Where: CYC s.7, "the Yukawa coupling ... `Y = 1/(u y_0^2 (Dt)^3)` ... normalisation `1/u`".
  * Add: this is relative to the Darboux-normalised B. The GW Yukawa coupling of FORMAL_QUINTIC_genus2 s.1 is `5 Y` (because
    `eta = 5 delta`), and the row coefficient `c_1 = Y Dt` has no 5 (VERIFIED against fq_rows).
  * Also, the `u` of the formula is `A_4`.
* **K7** (wording).
  * Where: SYN s.2 S9, "so Theorem U applies to the actual quintic".
  * Theorem U needs no CYC. CYC supplies two things: (i) the actual rows solve the system to which Theorem U is applied (UR s.4);
    (ii) every solution of that system gives a solution of PF5 (P(b) step 5).

## 8. Compliance and record

* **Files.** New files only, all in this folder: the two notes, eight scripts and thirteen logs. The `desktop.ini` files are created
  by the Drive client, not by me. `cy_extra.py` was written after the split and run once, on the final layout.
* **Project code.** Project modules were imported read-only, by `cy_proj` only, with `-B`. A `find -newer` showed that no file
  outside this folder changed during my work, apart from other agents' files in their own folders. No existing file was modified.
* **Own-file history.**
  * (1) My first `cy_run.py 24` hit a bug in my `cy_log.check` signature: the info string landed in the control flag, so the log
    showed spurious "CONTROL PASSED (ERROR)" lines. I fixed the helper and re-ran, overwriting my own `logs/cy_run_24.log`. This is
    the only overwrite.
  * (2) `cy_core.py` was first written as one 285-line piece, against the 250-line rule. I split it mechanically into `cy_core.py`
    (178 lines) and `cy_frob.py` (112 lines) and changed four import lines.
  * After the split, all twelve recorded runs were reproduced from copies of the final scripts in the scratch folder and compared,
    timings stripped. Eleven logs are identical. `logs/cy_rows_30_5.log` comes from the previous `cy_rows.py`: its four control
    lines lack the "fails at z^.." suffix, with the same checks and verdicts. `logs/cy_rows_30_6.log` is its final-code
    counterpart.
* **Environment.** Python 3.14 stdlib, one process at a time, longest run 13 s, nothing left running, nobody else's process touched.
  No wolframscript, no Mathematica, no cluster, no web access. The logs contain relative paths only.
* **Not done.** I did not referee PB s.2.5-s.2.7, Theorem U or GI Theorem G. I did not test colours `i != 0`: by Z/5 covariance,
  `Y_a(z) -> zeta^-a Y_a(zeta z)` preserves the system. I did not test the all-orders existence of the WKB rows. I did not read the
  special-geometry literature; CYC cites none and needs none.
