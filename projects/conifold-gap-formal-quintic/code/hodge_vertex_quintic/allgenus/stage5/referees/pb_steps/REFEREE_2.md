# Referee report, part 2: own code and results, corrections, what survives, record (2026-11-14)

Continuation of `REFEREE.md` (s.0 verdict; s.1 PB s.2.6; s.2 PB s.2.7; s.3 PB s.2.8-s.2.9; s.4 PB s.3). Same labels and notation.

## 5. Own code and results (VERIFIED: exact over Q, `fractions.Fraction`; every test with failing controls)

**Independence.** I wrote the code from the formulas re-derived in part 1. It imports no project module and opened no pb, poly,
evid, gamma, lemma or critic script.
* Only the closed forms of `G` and `log at` (PB s.2.1, POLY s.5) are transcribed, and I re-derived both by hand.
  * `Sum_(j=0..4) (1/x - zeta^j) log(1 - zeta^j x) = Sum_m x^(5m-1)/(m(5m-1))`.
  * `Sum_(j=0..4) (1 - zeta^j x)^(-n) = 5 Sum_(5|i) C(n-1+i, i) x^i`. The `j = 0` factor is included: the colour-0 summand contains
    `1/Gamma(1+d)`.
* Methods differ from pb's (as PB_AUDIT_2 s.4 describes them).
  * The hatted residue is computed DIRECTLY; pb compared only the key-lemma sum with `V`.
  * `Phi^SP` is the hatted SP over `Q[[u]]`; pb used the flat SP at numerical saddles in `F_p`, interpolated in `u`.
  * Everything is exact over Q; pb's mechanism and (S2) tests are mod `p`.
* *Limitation.* T1, T3 and T5 use the same amplitude formula on both sides. They test the relabelling, the mechanism and the
  identification, not the amplitude against the I-function, which is outside this task (POLY s.5, poly referee, PB/RF residue
  identity `V = C0 (Y^C_Q)_0`).

**Driver** `rps_run.py E NH NPI M` (with `rps_run_b.py`, `rps_run_c.py`). Final runs: `120 28 24 80` (132 s), `80 20 16 60` (20 s),
`20 6 5 8`; all ALL CHECKS PASSED. The numbers below are those of `logs/rps_run_v3_120_28_24_80.log`.

| test | what is compared | range | result |
|---|---|---|---|
| T1 relabelling (PB s.2.7) | `H = Res_X[e^(-What/w) B dX]` computed DIRECTLY at each hatted order (`e^(-f0/w)` whole: `Res_X[X^p w^q e^(-f0/w)] = w^q Sum_l (-1)^(a+l) w^(-a-l)/(a! 4^l l!)`), versus `V` (original variables, `Lam`-expansion), via `[u^N w^j] H = [u^(N+(1-4j)/5) z^j] V` | hatted `u^0..u^28`, every `w^j` that `V` through `u^120` reaches; 965 coefficients, 779 non-zero | PASS. Also: `H` has no positive w-power and none outside `4 mod 5`; `H_0 = Sum_l (-1)^(5l+1) w^-(5l+1)/((4l+1)! 4^l l!)` (RL's example); hatted `u^0` collects original `u^1, u^5, u^9, ...` |
| T2 mechanism, hatted (PB s.2.4, s.2.5, s.2.7) | CFU (`V^4 = tau^-5`, `y_1 = (1+T)/18`), inverse (`x(T) y(T x^5) = 1`), `beta(Y)`, reduction, grading; then `H = Sum_m pi_m Res_Y[Y^m e^(-f/w) dY]` | `u^0..u^25`; 4283 coefficients | PASS. Weight checks: every monomial of `beta` has `n + 4j <= 5i`, every `u^i w^j` of `pi_m` has `m + 4j <= 5i`; the a-priori degree bound equals `floor((5i+1)/4)` and is attained except at `i = 4q+3` (the (R1) zero) |
| T3 key lemma (PB s.2.7) | `Sum_m pi_m s^(-m) I^C_m(v(u), z)` in ORIGINAL variables, with `pi` truncated at `u^24`, versus `V` | 884 admissible keys `E + (4b-1)/5 <= 24` (684 non-zero, 200 with `b > 0`, 134 with `E <= 0`) | PASS. Top order contributes on 744 keys. Layer `= 25`: 30 of 32 keys wrong with `pi` to `u^24` (predicted), all right with `u^25` |
| T4 (S2), Prop. S (PB s.2.6) | pivots of `F_0[n]` in the closing relation; `F_0` (UR s.2 values); `F_b(0) = eps_b`; `eps_b SP^norm[eta^(c_b) e^(-(1/eta + eta^4/4)/t)] = F_b` (Lagrange-Buermann); the SP rows satisfy the T-system | `t^0..t^80`, `n <= 80` | PASS; pivot `= 4n` |
| T5 s.2.6 without the gauge | `Phi^SP` (hatted SP at `X = 1` over `Q[[u]]`, Lagrange-Buermann) versus `Sum_b eps_b pi_(c_b) tau^((1-c_b)/4) F_b(tau w)` | `w^0..w^25`, `u^0..u^25` | PASS. Also `ph_0 = 1`, `ph_1 = -3/20`, `ph_2 = 9/800`, `deg ph_n <= floor(4n/5)` where visible (degrees `0,0,0,2,3,4,4,5,..`) |

**Failing controls (each failed exactly where it must).**

| control | target | observed |
|---|---|---|
| `Lam/u -> 1` in the hatted phase (no regrouping factor `e^(-(Lam/u-1)/(Xw))`) | T1 | fails at the first compared coefficient of hatted `u^1` |
| hatted amplitude without `(1-uX^5)^(-1/2)`; with `b_k -> -b_k`; with an extra `u` per Stirling term | T1 | each fails at hatted `u^1` |
| original `G` cut to `x^4/4` | T1 | fails at hatted `u^1` |
| pull-back without `X'(Y)`; no CFU (`X = Y`) | T2 | fail at hatted `u^1` |
| reduction with `V -> 1` | T2 / T3 | fails at hatted `u^2` / on 684 of 884 keys |
| truncation of `pi` below the key-lemma bound | T3 | 30 of 32 keys of the first inadmissible layer |
| sign of `pi_-1` flipped; `pi_2` and `pi_3` exchanged | T3 | 630, resp. 679, of 884 keys |
| all signs `+` in (S2); `c_b -> c_b + 1`; system `T_k -> T_(k+1)` | T4 | all fail |
| without `eps_b`; without `tau^((1-c_b)/4)`; `F_b(w)` for `F_b(tau w)` | T5 | first at `(w^3,u^3)`, `(w^0,u^1)`, `(w^1,u^1)` |

The sign control of T5 first bites at `w^3 u^3`, which is predictable from the weight bound.
* `pi_3` (`b = 2`, `eps = -`) can first occur at `w^3 u^3`, and its coefficient there is non-zero.
* `pi_-1` (`b = 1`, `eps = -`) can first occur at `w^4 u^3`, but its coefficient there is the (R1) zero.

## 6. Corrections (the precise list; none changes a statement or a label of the chain)

* **C1 (PB s.2.0 versus s.2.6).** "(FL3) ... used only in P(a)": it is also used in (S2), as the v-derivative at the moving critical
  point. Cite FSL FL4 (s.4.5, full-SP form).
* **C2 (s.2.7, FL4).** The ring `R` is described imprecisely. Either require `ord_(1/w) c_n -> +oo` as `|n| -> oo`, which makes it
  a ring, or use FSL's ring `L` (s.10.1, refereed).
* **C3 (s.2.6, (S2) and Weights).** Cite L1a (linear maps over `A`) and Lemma H (RL s.3.2) for `x = sigma x'` with
  `z -> sigma^4 z'`, and for `Y = Y_c eta` with `w -> tau w`. This is RL's C1, applied here.
* **C4 (s.2.7, two expansions).** Name the family precisely: `(k, r)` together with the multisets of atoms of `l = log(e^(-G/z) at)`.
  Add two sentences:
  * the hatted fibre is finite, because `X^4/(4w)` is the only atom of hatted order 0 and the X- and w-balances fix `k` and its
    multiplicity;
  * `H` computed in `L[[u]]` is the hatted-order sum of the family, since the regrouping
    `e^(-v(Lam/u)/X) = e^(-v/X) e^(-v(Lam/u-1)/X)` holds coefficient by coefficient (REFEREE s.2.2).
* **C5 (s.2.7, key lemma; and K5 (iv)).** "Therefore a well-defined element of `Q((z))[[u]]`": finite fibres give a formal sum only.
  Add the weight lemma (REFEREE s.2.3): `5i - m - 4j >= 0` on `pi_m`, hence `E >= 0`. The chain does not need it.
* **C6 (s.2.8).** The four-item re-check of Theorem G for `P^SP` is redundant. `P^SP = s(z) P_can` is a regular gauge, and GI
  Theorem G holds for every regular gauge.
* **C7 (s.2.9).** State `[u^0](Y^C_Q)_0 = 0` (`y^(0) = -e_1`). The comparison at `u^1` uses it.
* **C8 (s.3).** Say where the two identities are compared (REFEREE s.4.2):
  * `A^can` and `A^QRR = A^can/d` lie in `Q[L, 1/L][[z]]`, inside `Q((s))[[w]]`, where UR s.4 holds. This also shows that the QRR
    amplitudes are Laurent polynomials in `L`, with no degree bound needed.
  * `D` is defined after stripping the same exponential on both sides. UR's `e^((u_i - u_i^c)/z)` differs from `e^(u_i/z)` by a
    constant that is not a power series in `z`.
* **C9 (s.3, "Same object").** Read: "Lemma N says `K_q = [z^(5q)] c`, where `Ghat(Phi)` is the hatted form of `P_can c`; for that
  gauge `P_11(0, z) = c(z)`". Add (CYC) to Lemma N's inputs: the rows are built by the cyclic system.
* **C10 (labels, s.1 H8).** Formal WKB uniqueness is elementary and can be labelled PROVED with one line. The exponent is a simple
  root of the symbol (`zeta^i L` distinct), so the transport operator is first order with kernel `Q L`. PB's label "standard" is
  weaker than necessary.

## 7. What survives

| item (PB section) | PB's label | after this referee |
|---|---|---|
| (S2) (s.2.6) | PROVED | PROVED (confirmed; C1, C3); VERIFIED `t^80` |
| weights, signs `eps_b`, triangular uniqueness, `kappa = C0` (s.2.6) | PROVED | PROVED (confirmed); the gauge-free identity VERIFIED `w^25 u^25` |
| FL4, FL5, model residues (s.2.7) | PROVED | PROVED (FSL s.10.1-10.3; C2) |
| relabelling `V = s Res_X[...]` (s.2.7) | PROVED | PROVED, complete with C4; VERIFIED (hatted `u^28`, `V` to `u^120`) |
| key lemma: bound and finiteness | PROVED | PROVED; bound attained; VERIFIED (884 keys, truncation layer) |
| key lemma: membership in `Q((z))[[u]]` | PROVED | PROVED with the weight lemma (C5); not needed by the chain |
| `V = C0 (P^SP Y^C)_0` (s.2.7) | PROVED | PROVED given s.2.6 and GI Lemma V; VERIFIED (T2 hatted, T3 original) |
| regularity, Theorem U, Theorem G for `P^SP` (s.2.8) | PROVED given H1-H3 | PROVED given (U), (CYC), GI Theorem G, POLY Lemma 2, P(a) (C6) |
| step 7, `s = 1` (s.2.9) | PROVED | PROVED (C7) |
| glue `D_0 = Phi^can(z,1)/h(z) = c(z;1)/h(z)` (s.3) | PROVED given P(b), H8, H9 | PROVED given the polynomiality of `Phi^can` (P(b), or Lemma 2 with (CYC)), H8 (C10) and H9 (QRR); C8 |
| EVID's `Phi^can` = POLY's `Phi^can` (s.3) | PROVED via Lemma N; VERIFIED (pb_glue) | PROVED given Lemma N (Theorem U, (CYC)); literally the same element (C9) |

* **On the critic's adjudications.**
  * K5: I agree that the argument is complete. Point (iv) needs the weight lemma (C5), and point (iii) needs C4's two sentences.
  * K12: I agree that (S2) is PROVED.
  * The synthesis lists PB s.2.6-s.2.7 and the glue among the (u) items. After this report they are refereed: sound, with the
    corrections C1-C10.
* **Not checked here (outside the task).**
  * P(a) (PB s.2.2) and the hatted/flat identity (S3);
  * (CYC) (CYC_PROOF; refereed by other agents);
  * Theorem U and GI Theorem G themselves;
  * the passage from `R(q=0) = diag(h_(w_i))` to H9;
  * pb's gauge-side check `eps_b pi_(c_b) = C0 G^can_b`, which needs a Theorem-U gauge that I did not build.

## 8. Files, runs, my own failures, compliance

| file | content |
|---|---|
| `REFEREE.md`, `REFEREE_2.md` | this report |
| `rps_core.py` | Bernoulli numbers; `V` in original variables; the hatted integrand, its direct residue, `B` |
| `rps_mech.py` | CFU (`(T-1) y_k + V_k = R_k`), inverse map, `beta(Y)`, reduction, the residue sides (hatted, and original for the key lemma) |
| `rps_s2.py` | Prop. S (T-system, pivots), (S2) by Lagrange-Buermann |
| `rps_sp.py` | hatted `Phi^SP` over `Q[[u]]`; the right side of the s.2.6 identity with controls |
| `rps_run.py`, `rps_run_b.py`, `rps_run_c.py` | driver (T1-T2; T3-T4 and summary; T5); the second and third parts are executed inside the first |
| `logs/` | 8 logs, relative paths only (list below) |

**Runs** (all `py -I -B`, Python 3.14, stdlib, working directory in the session scratch folder; one process at a time).
* `v3` logs come from the final code: `rps_run_v3_120_28_24_80` (132 s), `rps_run_v3_80_20_16_60` (20 s), `rps_run_v3_20_6_5_8`.
  `rps_run_v3_40_10_8_24` predates T5.
* `rps_run_16_5_4_8`, `rps_run_40_10_8_24` and `rps_run_80_20_16_60` (no tag) come from the first version, before the T2 degree
  diagnostic and before T5. Their checks are a subset of v3's, and all passed.
* The driver raises `NH` to `NPI + 1` before naming the log, so the first log's name and title say `16 5 4 8` for the command
  `16 4 4 8`.

**My own failures (reported plainly).**
* The first command (`16 4 4 8`) crashed with an IndexError in T2: the direct `H` was computed only to `NH < NPI + 1`. No log was
  written. Fix: `NH = max(NH, NPI + 1)`.
* **A wrong prediction of mine.** In version 2 I added the check "the regularity bound is NOT implied by the reduction alone". It
  FAILED (`logs/rps_run_v2_40_10_8_24.log`, kept as written).
  * The a-priori bound of the reduction equals `floor((5i+1)/4)` exactly, because every object of the mechanism has weight
    `omega >= 0` (REFEREE s.2.3).
  * I replaced the check by the three weight checks (v3) and turned the observation into the weight lemma.
  * This was an error in my expectation, not in PB.
* No other failure. No log was overwritten: each version writes under a new name.

**Compliance.**
* New files only, all inside `stage5/referees/pb_steps/` (2 notes, 7 scripts, 8 logs). I created the folders `referees/pb_steps/`
  and `logs/`; the parent `referees/` already existed (other referees' folders).
* No existing file was modified, moved or deleted; no project script was run.
* Every file is under 250 lines (largest: REFEREE.md, 225; `rps_mech.py`, 207).
* `py -I -B` (no `__pycache__`). The scripts add their own folder to `sys.path`, because `-I` drops it.
* At most one of my processes ran at a time, and all have finished; no one else's process was touched. No wolframscript,
  Mathematica or cluster job.
* No web access was used: the steps refereed here cite no external statement.
* No absolute path and no user name appears in any note or log. The scratch folder holds only a copy of one run's stdout.
