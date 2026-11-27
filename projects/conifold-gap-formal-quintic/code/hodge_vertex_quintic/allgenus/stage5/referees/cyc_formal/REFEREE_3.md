# Referee report `cyc_formal`, part 3: own code, errata, record (2026-11-14)

Continuation of `REFEREE.md` (s.0-s.5) and `REFEREE_2.md` (s.6-s.9). Same labels and abbreviations.

## 10. Own code (exact over Q; Python 3.14 stdlib; `py -I -B`, working directory the session scratch folder)

**Modules.** I wrote these from scratch for this report:
* `cyc_ref_core.py`: truncated Laurent series with absolute precision, and polynomials in `(alpha, beta)`;
* `cyc_ref_ops.py`: the ring `R = K[l]`, the map `l -> l + c`, and operators with left coefficients (composition, application, formal
  adjoint);
* `cyc_ref_jet.py`: bilinear jets, the symbolic-n calculus, and PF4 built from its large-radius form;
* `cyc_ref_frob.py`: the Frobenius basis, the obstruction matrix, and B from the exact jet matrix.

**Independence.** I did not import or copy `critic/cr5_cyc*.py`, `pb/pb_cyc.py`, `pb/pb_core.py` or any other agent module. The only
project code I ran is in `cyc_ref_xproj.py`: scratch copies of `fq_conifold.py`, `fq_fp.py` and `fq_series.py`, imported read-only, for
item (1).

**Methods differ from the critic's and from pb's:**

| step | mine | the critic's / pb's |
|---|---|---|
| the recursion | symbolic in `n` | the critic: residuals of PF4 on monomials |
| B | the exact 4x4 jet matrix | the critic: the nested definition |
| (CYC) | the operator-coefficient expansion of `M` | the critic: `M` applied to solutions and to `u^0..u^4`; pb: to `u^0..u^6` |

I also test things neither of them tests:
* a generic frame with indeterminate coefficients;
* the Wronskian determinant;
* all 16 entries of C3;
* the Leibniz identities on random data.

**Runs.** All of these passed.

| script (log) | checks | range | controls (all fail as they must) | time |
|---|---|---|---|---|
| `cyc_ref_alg.py` (`logs/cyc_ref_alg.log`) | A1 PF4, UR u-form; A2 `D^5 - q prod(5D+k) = D PF4`; A3 `c_s(n)`; A4 C1 and its two conditions; A5 C2(i) and alternation (jets); A6 `P_ab`, `det P = u^4`, `P^-1` | exact (no truncation) | `q` on the right; `k = 0..4`; `PF4 + uD`; `A_3 -> -3(1-u)`; `q` on the right in C1; no signs `(-1)^j`; alternation for `PF4 + uD` | 0.2 s |
| `cyc_ref_sol.py 80` (`logs/cyc_ref_sol_80.log`) | S1 Frobenius; S2 obstruction matrix, kernel `(0,1,7/10)`, `[u^-3] r = -2 f(0)`, `y_3`; S3 B log-free, constant, alternating, values, `det = 4`, `Pf = 2`; S4 `det Wr = -2/u^2`; S5 C3, 16 entries; S6 C4 | PF4 kills `f_i`, `y_3` below `u^77`; B constant below `u^78`; Wr below `u^76`; C3 below `u^75` | `gamma = 1/2`; `g + f_0`; concomitant without `k = 1`; perturbed `Binv`; `l -> l + u` (twice); Lemma 0 with `D l' = 1 - u` | 15 s |
| `cyc_ref_frame.py 80` (`logs/cyc_ref_frame_80.log`) | F1 C5: adapted basis, `tau_k`, `Y`; F2 operator expansion `M = PF4`, `c_2 - c_3 t`; F3 Corollary operator form, `prod c = 1/u`, UR/GI `b_a`, `t_c` | five frames (s.8.1); below `u^76` / `u^77` | `t = f_1/y_0`; `t = (f_1 + f_2/2)/y_0`; middle factor `Dt^3`; `c_1` with `Dt^3`; `c_2 <-> c_3`; C5's span for `t = f_1/y_0` | 19 s |
| `cyc_ref_sym.py 30` (`logs/cyc_ref_sym_30.log`) | Y0 the symbols are live; Y1 `B(y_0, g) = 0`; Y2 `M = PF4`; Y3 Corollary; all identically in `(alpha, beta)` | below `u^27` (Y2), `u^26` (Y3) | `y_1 = f_1` | 47 s |
| `cyc_ref_leib.py 40` (`logs/cyc_ref_leib_40.log`) | L1-L2 Leibniz identities; L3 adapted-basis expansion; L4 `D(1/e)D = e th^2`, the th-form of `M`, leading coefficient `u` (random `y_0`, `t`, `f_i`, antisymmetric matrix) | below `u^34` | a wrong binomial; a missing factor 3; the th-form with `Dt^2` | 0.5 s |
| `cyc_ref_xproj.py` (`logs/cyc_ref_xproj.log`) | X1-X3: project code = mine (`f_i`, `gamma`, `h`, pairing matrix and flags, `varpi0`, `tc` in two frames) | below `u^28` | (comparison) | 2 s |

**Re-run of the critic's code.** I ran scratch copies of `cr5_cyc.py` and `cr5_cyc_run.py` with the stored arguments 24, 48 and 80.
All three logs are identical to the stored `critic/logs/cr5_cyc_run_{24,48,80}.log` once timings are stripped. One caveat: the
critic's gamma check searches only `span(f_1, f_2)` (E1).

## 11. Errata (precise; none affects validity)

**E1. CYC s.1** says the condition "holds exactly for `g = f_1 + (7/10) f_2` (computed, `cr5_cyc_run`)".
* The cited code, like `fq_conifold.solve_log`, searches only `span(f_1, f_2)`.
* Replace with: "on `span(f_0, f_1, f_2)` the obstruction matrix (rows `u^-3, u^-2, u^-1`) is
  `[[-2, 0, 0], [17/5, 0, 0], [-7/5, -7/5, 2]]`. Its kernel is `Q (0, 1, 7/10)`, so `g` is unique up to scale."
* This is VERIFIED (`cyc_ref_sol` S2) and checked by hand (REFEREE s.2.3).
* The uniqueness is needed only for the words "the vanishing period", not by the proof.

**E2. CYC s.3 (iv)** says "computed exactly from the solutions".
* Add: "`Pf(Bmat) = -B(f_0, f_2) B(f_1, y_3) = -(2)(-1) = 2`, since `B(f_1, f_2) = B(f_2, y_3) = 0`. These four values are `u^0`
  coefficients (REFEREE s.4.4)."
* Optionally add: "`det Wr = -2/u^2` (Abel's identity and Lemma 0)."
* C2 (iv) is then proved without reference to code.

**E3. CYC s.8 (rows (d), (e)) and STAGE5_SYNTHESIS (s.1 table "3 frames", s.4.1).**
* The critic's frames are `y_0 = f_0 + a f_1 + b f_2`.
* In the project's parametrisation `f_0 + alpha g + beta f_2` they are `(0,0)`, `(1/3, -109/210)` and `(7/11, -23/22)` (VERIFIED).
* They are not pb's `(1/3,-2/7)` and `(7/11,-3/5)`. Relabel them.

**E4. CYC s.0 and the brief.** Read `fq_conifold.py` (package root) for `natural_form/fq_conifold.py`. The latter does not exist;
byte-identical copies exist only under `natural_form/genus5_package_review/`.

**E5. CYC s.0, last line.** Read "in every conifold frame (every power-series solution `y_0` with `y_0(0) != 0`, every non-zero scale
of `g`)" for "in every frame". The large-radius frame and the ring N are not covered, and they are not needed (REFEREE_2 s.9).

**E6. Notation clashes in CYC.**
* `M` is both the factorised operator (s.0, s.7) and `zD` (the Corollary).
* `tau` is both the automorphism (s.5) and `tau_k` (s.6); elsewhere in the project `tau(u) = 1/2F1(1, 4/5; 9/5; u)`.
* `f_i` is both the Frobenius basis and `y_i/y_0` (s.6).
* `P` is both `prod(D + k/5)` (s.2) and the concomitant / jet matrix (s.3).
* `y_3` (the log solution) sits next to `y_3'` (the adapted basis).
* Suggestion: write `M_cyc`, `sigma` for the automorphism, `phi_i = y_i/y_0`, and `Pi_ab` for the jet matrix.

**E7. CYC s.6.** Add: "`tau_1 = tau_2 = 0` and `tau_3 = 1/(u y_0^2 Dt^3)` hold for every unit `y_0` and every `t` with `Dt` a unit,
since C3 is basis-free. Isotropy enters only through the adapted basis, which exists iff `B(y_0, y_1) = 0`." (REFEREE_2 s.7.)

**E8. CYC s.0, Corollary.** Add: "The identities lie in `Q[[u]][D]`. They hold on every differential `Q[[u]]`-module in which `z` is
constant. That includes the s-series rows (`D s = -(1-s^5)/(5 s^4)`) and, after conjugation by `e^(u_i/z)`, the rows `X_a` of
`fq_rows`." (REFEREE_2 s.8.2.)

**Status lines.**
* CYC s.9 "PROVED (critic; AI-written; unrefereed)" now reads "refereed once (`referees/cyc_formal`)".
* "VERIFIED twice (pb, critic)" now reads three times.
* STAGE5_SYNTHESIS s.1 table, "(CYC) ... what would close it: referee CYC": done (AI); a human check remains.

**Aside, not CYC.** `fq_conifold`'s docstring says "elements that are Lam-free in N give Lam-independent images". Read
"Lam-independent when written in the `D^j I_1`". For example `i_0` is Lam-free in N but its image depends on the chosen value of
`Lam`. The rows are of the good kind (REFEREE_2 s.9).

## 12. Record and compliance

**Files** (new, in `allgenus/stage5/referees/cyc_formal/`):

| file | content |
|---|---|
| `REFEREE.md` | s.0 verdict; s.1-s.5 (items 1-5) |
| `REFEREE_2.md` | s.6-s.9 (items 6-9) |
| `REFEREE_3.md` | s.10-s.12 (code, errata, record) |
| `cyc_ref_core.py`, `cyc_ref_ops.py`, `cyc_ref_jet.py`, `cyc_ref_frob.py` | the modules of s.10 |
| `cyc_ref_alg.py`, `cyc_ref_sol.py`, `cyc_ref_frame.py`, `cyc_ref_sym.py`, `cyc_ref_leib.py`, `cyc_ref_xproj.py` | the checks of s.10 |
| `logs/cyc_ref_alg.log`, `logs/cyc_ref_sol_80.log`, `logs/cyc_ref_frame_80.log`, `logs/cyc_ref_sym_30.log`, `logs/cyc_ref_leib_40.log`, `logs/cyc_ref_xproj.log` | ALL CHECKS PASSED; relative paths only |

`cyc_ref_xproj.py` expects scratch copies of the three project modules in a subfolder `xproj/` next to it. Those copies are not
stored here.

**My own failures.** Each was fixed before any recorded run.
* The first `cyc_ref_core.py` had 290 lines. I split it into `cyc_ref_core.py` and `cyc_ref_ops.py`.
* The first `cyc_ref_frame.py` tested "`u y_0^2 Dt^3 th^2(y/y_0) = c_2 - c_3 t`" for l-freeness using the key set of the element. The
  elements that depend on `y_3` carry an `l^1` key whose series is zero, so the check FAILED in all five frames at `N = 20`. I changed
  it to test that the `l^k` parts (`k >= 1`) vanish. A debugging script confirmed the cause.
* I reworded one control label (Lemma 0); the computation is unchanged.

No other failure.

**What I did not do.**
* I did not read the special-geometry literature (CYC s.9 cites hep-th/9204035). CYC is self-contained, and I checked it as such.
* I did not re-run `pb_cyc`.
* I did not audit the project's `Cf` map beyond its docstring and REFEREE_2 s.9.
* I did not referee the uses of (CYC) beyond the direction and the frame (PB s.2.1, s.2.3, s.2.8, s.3; UR s.1; GI s.1).

**Compliance.**
* Only new files, all in my output folder. I created the parent folder `stage5/referees/`, which did not exist; the `desktop.ini`
  files in it are the file-sync client's.
* No existing file was modified, moved or deleted. Project and critic scripts were run only from scratch copies.
* Python stdlib only. One process at a time; every process I started has finished. No wolframscript, no Mathematica, no cluster job,
  no web access.
* Logs and notes contain no absolute path and no user name. Every file is under 250 lines.
