# Stage 5 referees: synthesis, part 3 — the critic's checks, human checks, next steps, record (2026-11-14)

Continuation of `REFEREES_SYNTHESIS.md` (s.0-s.4) and `REFEREES_SYNTHESIS_2.md` (s.5 corrections, s.6 final status). Same labels
and abbreviations. AI-written; no human has read it.

## 7. The critic's own checks (Python 3.14 stdlib, `py -I -B`, working directory the session scratch folder)

### 7.1 Re-runs of the referees' code from scratch copies
Each referee's scripts were copied to the scratch folder and run there with the stored arguments, at most three of my processes at a
time. For `cyc_ref_xproj` (CF) and `cy_proj` (CA), byte-identical copies of `fq_conifold`, `fq_ring`, `fq_rows`, `fq_fp`, `fq_series`,
`fq_hodge`, `fq_ifun` were placed so that the scripts' relative imports resolve (`cmp`-checked against the package root). The comparison
is `critic/rc_rerun_compare.py` (timings `d.d s` replaced on both sides; CF's stored logs carry one extra header line `# logs/...`,
which CF added when saving stdout, dropped); its output is `critic/logs/rc_rerun_compare.log`.

| referee | scripts and stored arguments | result |
|---|---|---|
| cyc_formal | `cyc_ref_alg`; `cyc_ref_sol 80`; `cyc_ref_frame 80`; `cyc_ref_sym 30`; `cyc_ref_leib 40`; `cyc_ref_xproj` | 6 logs identical |
| cyc_adversary | `cy_run 64`; `cy_frames 80`; `cy_rows 64 12` and `30 6`; `cy_proj 40`; `cy_extra 40` | 6 logs identical |
| pb_steps | `rps_run 80 20 16 60` and `120 28 24 80` (v3 logs) | 2 logs identical |
| lmi_formal | `lf_dlmf 60 80`; `lf_quintic 14 60`; `lf_lmi 12 12`; `lf_lmi2 14`; `lf_fit 30` | 5 logs identical |
| lmi_adversary | `lm_run_eps 10 Q`, `12 p1`, `14 p2`; `lm_run_q 12 Q`, `16 p1`, `16 p2`; `lm_run_x 12` (the 27 falsifiers) | 7 logs identical |

All 26 re-runs end in ALL CHECKS PASSED. Not re-run: the lower-range duplicates (CA `cy_run 24, 40`, `cy_frames 24, 48`,
`cy_proj 20`, `cy_rows 44 9`; PS v3 `20 6 5 8`, `40 10 8 24`; LF's lower ranges), the superseded logs (PS's first version; CA
`cy_rows_30_5`), the two drafts, and the Stage-5 critic's `cr5_cyc_run` (re-run by CF and by CA, identical). About 24 CPU-minutes.

### 7.2 New code (`critic/`, exact over Q)
* **`rc_cyc.py N`** (with `rc_core.py`, `rc_cyc_frames.py`): (CYC) and every point the two CYC referees raise, 62 checks.
  * *Method, new relative to the five earlier CYC checks* (cr5 on solutions and on `u^0..u^4`, pb on `u^0..u^6`, the main session's
    `cyc_check` on `u^0..u^5`, CF and CA coefficientwise on the operator): every operator identity is tested on `u^x` at six NON-integer
    rational `x` through the twisted derivation `Dx g = x (1 - 1/u) g + D g`. For an operator `E = Sum_(k<=4) e_k D^k` over K,
    `u^-x E u^x` is a polynomial in `x` of degree `<= 4` whose `x^k` coefficient is `e_k (1 - 1/u)^k` plus terms in `e_(k+1..4)`; so five
    distinct `x` determine E. The concomitant is computed from its definition in `K[l]` (not from the jet matrix), `det Wr` by the
    Leibniz expansion over `K[l]`, the Frobenius basis from the pivots `[u^(n-3)] PF4 u^n`.
  * *Results* (`logs/rc_cyc_40.log`, `logs/rc_cyc_60.log`: ALL CHECKS PASSED). Exact: PF4 from its large-radius form; `D^5 - q
    prod_(k<=5)(5D+k) = D o PF4`; C1; PF4 maps `Q[[u]]` to itself; the full obstruction matrix `[[-2,0,0],[17/5,0,0],[-7/5,-7/5,2]]`
    (CF-E1) and `gamma = 7/10`; `Bmat` with the hand values `b01 = -7/5, b02 = 2, b03 = -9/10, b12 = 0, b13 = -1, b23 = 0`, `Pf = 2`
    (CF-E2); the `u^-1` obstruction row equals `B(f_0, .)` (CA s.3); C4; the frame relabelling (A2). In range (`N = 60`): `f_i` solve
    PF4 below `u^70`, `y_3` below `u^68`; B log-free and constant below `u^64`; `det Wr = -2/u^2` below `u^62`; `M = PF4` below `u^70` in
    six frames with `y_0(0) != 0` (programme `(0,0)`, `(1/3,-2/7)`, `(7/11,-3/5)`; critic `(1/3,-2/7)`, `(7/11,-3/5)`; `(2 f_0, 3 g)`) and
    below `u^69` in four with `y_0(0) = 0` (`f_1`, `f_2`, `f_1 - 2 f_2`, pair `(f_1, f_2)`; A1); the Corollary as the operator identity
    `(1/c_4)D(1/c_3)D(1/c_2)D(1/c_1)D(1/c_0)D o (1/y_0) = (1/y_0) D PF4` below `u^66`-`u^69` in the same ten frames; UR s.0 item 1's
    `b = -c/(1-u)`; `t_c = u + (7/10)u^2 + (41/75)u^3 + (9/20)u^4`; unit shapes (A1).
  * *Controls, each fails:* `PF4 + uD` self-adjoint; `q` on the right of the product; `g + f_0` and `gamma = 1/2` unobstructed;
    `M = PF4` for the non-isotropic `(f_0, f_1)` and for `gamma' = 1/2`, `701/1000` (whose coefficients nevertheless keep the
    Theorem-U shape: CA-K3 reproduced); middle factor without `u`; the Corollary with `c_1`, `c_2` swapped and with `5 c_1` (CA-K6).
* **`rc_l3c.py NW NL`** (`logs/rc_l3c_60_12.log`: ALL CHECKS PASSED, 8 checks): formal DLMF 5.11.8 as an identity of polynomials in
  `h` for `w^-1..w^-60`, and the two inputs of its interpolation proof (`B_m(h+1) - B_m(h) = m h^(m-1)`, `m <= 61`; `beta_n(0)` is the
  Stirling series); Lemma L3-C as an identity of polynomials in a SYMBOLIC `M` through `z^12` (Wick sum of
  `N[log t - t, t^(M-1); 1]`; `[z^n]` has degree exactly `2n`), and for `M = -3..4` the Pochhammer form of (L3). Controls fail
  (`M -> M + 1`, i.e. the wrong measure; `beta_n` at `1/5` in place of `M = 0`).
* **`rc_rerun_compare.py`**: s.7.1.

### 7.3 By hand
* CYC s.0-s.7 line by line: Lemma 0; C1 (`P(-D-1) = P(D)`); C2 (i) telescoping, (ii) alternation, (iii) the anti-diagonal
  `(u, -u, u, -u)`, (iv) `det P = u^4`; C3 (`P^-1 = Wr Bmat^-1 Wr^T`, anti-triangular); C4 (`tau` commutes with D because `D l in K`);
  C5 (adapted basis, `Bmat^-1 = -Bmat`, the Leibniz expansions, `tau_2 = -th^2 f_3 - tY`, `tau_3 = Y`); the Theorem; the Corollary
  (`q(5D + 5) = 5Dq`). Also `O_-3 = -2 b_0`, `O_-2 = (17/5) b_0` (`P_-3'(1) = 0`), `O_-1 = 2 b_2 - (7/5) b_1 - (7/5) b_0`, and CF's
  `Y_1 = -1 + (3/10)u`, `Y_2 = u^-1 - 9/10`. No error; I agree with CF and CA.
* PS: (S2) row by row (`bmod = (1, -1/v, 1, -1, -1)`), the weight lemma (A6), the glue's steps (a)-(e) (A4).
* LF, LA: the (L2') moments and resummation; RA s.3.3 steps 1-3 against LF s.3 and LA s.1; Lemma F5.11.8 and Lemma L3-C (A3, A7).

## 8. What a human must check (SYN2 s.7, updated), and next steps

1. **Theorem Gamma:** RA s.3 (one page) with LF s.3 (the identification of the critical points, LF-C1), or CP s.3; with (L1a), (L1b),
   (L1c)(i), (L2'), (L3) (FSL s.5, s.6, s.4, s.7, s.9.3).
2. **(CYC):** CYC_PROOF.md (two pages) with CF s.2.3, s.4.4 (the two hand computations) and the scope of correction 5 (part 2).
3. **Theorem P(b):** POLY s.6 with PB s.2.5-s.2.9 and PS s.1-s.3: (S2), the relabelling with PS-C4, the key lemma and the weight lemma.
4. **The glue:** PB s.3 with PS s.4.2 and UR s.4.
5. **(QRR):** as SYN2 s.7 item 5 (GJR arXiv:1812.11908 eq. (17) typeset, its sources, Teleman s.8.2 verbatim).
6. **Theorem U** (UR s.3) and **GI Theorem G** (GI s.3). 7. **The Stage-2 list** (SYN2b s.6.2).

*Next steps.* (1) The write-up of SYN2 s.8 item 2, with the corrections of part 2, s.5: (CYC) with the wider scope and the hand
computations; the glue with PS-C8; Theorem Gamma by RA s.3 with LF-C1-C5 and LA-C1-C6; `Ginf` by the formal 5.11.8. (2) The owners of
CYC_PROOF, SYN, PB_AUDIT, RA, GEN, POLY, GI and of the project files apply part 2, s.5 (no agent of this round edited any of them).
(3) Optional: a second AI referee of PB s.2.6-s.3, the only chain steps refereed once in this round; a human check is worth more.
(4) Optional, unchanged: the genus-3 code test `F_3(0) = 1/224000` of SYN2 s.7.

## 9. Record

**Files** (all new; the three synthesis files at the place the brief names, everything else in `critic/`):

| file | content |
|---|---|
| `REFEREES_SYNTHESIS.md`, `_2.md`, `_3.md` | this synthesis |
| `critic/rc_core.py` | Laurent series in `u` with absolute precision (inverse by the coefficient recursion), `K[l]` with `D l = 1 - 1/u`, the twisted derivation `Dx`, operator chains |
| `critic/rc_cyc.py`, `critic/rc_cyc_frames.py` | the (CYC) re-check of s.7.2 (the second file is executed by the first) |
| `critic/rc_l3c.py` | formal 5.11.8 and Lemma L3-C (s.7.2) |
| `critic/rc_rerun_compare.py` | the log comparison of s.7.1 |
| `critic/logs/rc_cyc_40.log`, `rc_cyc_60.log`, `rc_l3c_60_12.log`, `rc_rerun_compare.log` | ALL CHECKS PASSED (62, 62, 8 checks); 26 of 26 logs identical; relative paths only; the timing is left out of the final line |

**Runs.** `rc_cyc 12` (development), `40`, `60`; `rc_l3c 20 6` (development), `60 12`; `rc_rerun_compare` twice (25 logs, then 26 after
the falsifier re-run); the re-runs of s.7.1. The recorded logs are from the final code.

**My own failures (all before any recorded run).**
* `rc_cyc.py`, found by reading before the first run: the formal adjoint truncated its accumulator; the control "q on the right" was
  built wrongly; the solver lacked a precision guard for its right-hand side. `rc_core.LS.inv` would have looped (almost) for ever on
  an exact series; it was replaced by the coefficient recursion with a cap. All fixed before any run.
* `rc_cyc 12`, first development run: two FAIL lines, both mine. Two TRUE statements ("`g + f_0` is obstructed", "`gamma = 1/2` is
  obstructed") had been wrapped as controls, so the harness reported failed controls. They were rewritten as controls of the false
  statements ("... is unobstructed"); the computation was unchanged. That development log, in the scratch folder, was overwritten.
* `rc_l3c 20 6`, first development run: a crash (Bernoulli table too short); no log; fixed. One control label was reworded later.

**Compliance.**
* New files only. The brief asks for files only inside `referees/critic/` and also names `referees/REFEREES_SYNTHESIS.md` (+ `_2`);
  I wrote the synthesis there (three pieces, each under 250 lines) and nothing else outside `critic/`, which I created.
* No existing file was modified, moved or deleted. Referee and project scripts were run only from scratch copies (they write their logs
  next to themselves); the main session's `independent_check/` was read (README, log tails) as data and not run.
* Python 3.14 stdlib, `py -I -B`, working directory the session scratch folder; at most three of my processes at a time; every
  process I started has finished; none of anyone else's was touched. No wolframscript, Mathematica, cluster job or web access.
* Notes and logs contain no absolute path and no user name. Every file is under 250 lines.

**What I did not do.** I did not referee Theorem U (UR s.3), GI s.3, the Stage-2 chain, (QRR), or the proofs of the formal lemmas
(FSL, CP2); PB s.2.2-s.2.5 were read as context only; I did not re-run the lower-range logs, and I read no literature myself.
