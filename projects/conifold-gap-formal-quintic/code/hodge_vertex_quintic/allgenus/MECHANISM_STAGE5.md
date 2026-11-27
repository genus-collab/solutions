# All-genus mechanism for Conjecture N: Stage 5, the proof of (Gamma) (2026-11-14)

Entry point for Stage 5. It supersedes `MECHANISM_STAGE4.md` s.2 (assembled statement) and s.5 item 1; that note stays unchanged as
the record of Stage 4. Written from two workflows and my own checks (s.4). The first workflow had 9 agents: workers chart, lemmas,
general, pb and lit; two referees on chart and one on lemmas; a critic. The second had 6: referees of the three pieces the first left
unrefereed, and a critic. Everything is AI-written. Every step of the chain is AI-refereed; no human has checked any proof. Labels as
in `../NATURAL_FORM.md`.
Full record: `stage5/STAGE5_SYNTHESIS.md` and `_2.md` (first workflow: bottom line, assembled theorem, proof chain S1-S9,
contradictions K1-K14, overclaims O1-O19, human checks); `stage5/referees/REFEREES_SYNTHESIS.md`, `_2.md`, `_3.md` (second workflow:
verdicts, four disagreements settled, all 38 corrections, final status, the critic's checks). Notes: `stage5/chart/CHART_PROOF.md`
(+`_2`..`_4`), `stage5/lemmas/FORMAL_SP_LEMMAS.md` (+`_2`..`_5`), `stage5/general/GENERAL.md` (+`_2`), `stage5/pb/PB_AUDIT.md` (+`_2`),
`stage5/lit/LIT.md` (+`_2`, `_3`), `stage5/critic/CYC_PROOF.md`. Referee reports: `stage5/chart/referee_formal/`,
`stage5/chart/referee_adversary/` (which also contains Theorem LMI, s.3), `stage5/lemmas/referee/`, `stage5/referees/*/`. Written into
the project by me: this note and `stage5/independent_check/` (s.4.2); my re-run outputs stayed in the session scratch folder.

## 1. The result

* **(Gamma) is PROVED**: the full-Stirling saddle series `c(z; L)` of the Mellin-Barnes summand has no pole at `L = 1` (`q = 0`),
  and `c(z; 1) = h_(w_0)(z)`. There are two independent proofs, each refereed twice:
  * the 7-variable chart (`chart/CHART_PROOF_3.md` s.3), the Stage-4 outline written out with every power of `z` and every
    prefactor: `referee_formal` "sound with corrections" (six minor ones), `referee_adversary` "sound";
  * a one-page two-variable chart, Theorem LMI ("local Mellin inversion", `chart/referee_adversary/REFEREE.md` s.3):
    `referees/lmi_formal` "sound with corrections" (presentational), `referees/lmi_adversary` "sound" (27 falsifiers fail as LMI
    predicts).
* Its only inputs are formal stationary-phase lemmas, each proved twice independently (CHART_PROOF_2 s.2; FORMAL_SP_LEMMAS s.1-s.9)
  and refereed. Stirling as the stationary phase of Euler's integral (L3) has a purely formal proof (FORMAL_SP_LEMMAS s.9.3), so
  no analytic theorem is needed.
* **Hence `D_0 = 1` to all orders, (F_g) for every `g`, and Conjecture N(g) for every `g >= 2`, modulo (U), (Q) and (QRR)** (s.2).
* **A missing hypothesis, found and proved.** The second audit of Theorem P(b) found that (CYC), the cyclic form of the quintic
  system, is used in P(b) step 5 and, through UR s.1, by the whole programme, but appeared in no hypothesis list. The Stage-5 critic
  proved it (`critic/CYC_PROOF.md`); both referees found it sound (s.3.4).
* **Free coefficients, by-product.** The two-variable chart at the orbifold point gives `c(z; 0) = Ginf(z)` (refereed twice). So the
  closed form `Sum_m f_m z^(5m) = exp(Sum_l (-1)^(l+1) B_(5l+1)(1/5) z^(5l)/(l(5l+1)))`, CONJECTURAL in Stage 4, is PROVED given P(b).
  It also gives `c_n in Q[L]` directly.
* **(QRR) now has a literature statement for the formal quintic itself**: Guo-Janda-Ruan arXiv:1812.11908v1 eq. (17), proof of
  Lemma 5.2. Stage 4 cited arXiv:math/0108100 Thm 9.1, which is the untwisted statement.
* **Nothing was refuted.** The generalisation falsifier found the mechanism correct in every theory where its hypotheses hold, and
  failing exactly where one fails (s.3.5). No order failed in any field, frame or route. Every FAIL line in the logs is a deliberate
  control, except two draft logs kept and explained by their referees (wrong predictions in the referees' own tests).

## 2. The assembled statement after Stage 5

**Theorem (AI-written and AI-refereed; no human check).** For every `g >= 2`, Conjecture N(g) holds modulo exactly:
* **(U)** Theorem U (`natural_form/theory/UNIFORM_REDUCTION.md` s.3): proved in the project, AI-audited only.
* **(Q)** the quoted theorems, unchanged since Stage 2: arXiv:1509.06954 Thm 8.1; arXiv:2212.00320 Thm 7.1; arXiv:math-ph/0702045
  Thm 5.1 with the linear loop equation of its s.4 (on the model curve `x^4 + v/x` only, as allowed); arXiv:1303.1043 Def. 2.2,
  Prop. 2.4, Def. 2.5, Prop. 2.9; optionally arXiv:2502.18115 Cor. 1.2.
* **(QRR)** `R(q = 0) = diag(h_(w_i))`, `w_i = {zeta^i - zeta^j}_(j != i) U {-5 zeta^i}`, used only in the form `A^QRR_n(q = 0) = h_n`.
  * Sources: GJR arXiv:1812.11908v1 eq. (17); GJR arXiv:1709.07392v1 eqs. (19)-(20) at `t = 0`. Both justify it by citing
    arXiv:math/0008067 Thm 1 with its Remark and arXiv:math/0110142 Thm 1, Cor. 1.
  * A vertex-level derivation (lit, PROVED given its quotes): arXiv:1803.01409v2 eq. (26), Mumford's relation, arXiv:math/0108100v2
    Thm 4.1.
  * Uniqueness of R for the CohFT: arXiv:0712.0160v3 Theorem (2) (s.1.9) and Prop. 8.5(ii) (s.8.2), read only through a summariser.
  * The twist sign `-5` is settled three ways and pinned by `h_1 = ph_1 = -3/20`.

The chain (STAGE5_SYNTHESIS s.2, S1-S9):
1. the formal lemmas;
2. Theorem Gamma;
3. flat to hatted, `c_n(L) = L^(4n) ph^SP_n(L^-5)`, so `Phi^SP(w; 1) = h(w)`;
4. Theorem P(b), `Phi^SP = Phi^can`;
5. the glue, `D_0(z) = Phi^can(z, 1)/h(z)`, an exact identity (pb s.3; it uses UR s.4, formal WKB uniqueness and (QRR); Lemma N's
   formula for `K_q` is no longer needed). So `D_0 = c(z; 1)/h(z) = 1`, and all colours follow from `D_i(z) = D_0(zeta^-i z)`;
6. (F_g) for all `g`, then the unchanged Stage-2 chain.

Genus 2-5 remain VERIFIED directly by graph sums.

| input or step | status after Stage 5 | what would close it |
|---|---|---|
| (U) | PROVED (UR s.3), AI-audited; VERIFIED in range (independent construction regular through `z^80`, pb) | human check of UR s.3 |
| (Q) | quoted with statement numbers; hypotheses checked in Stages 1-2 | human check (STAGE2_SYNTHESIS_2 s.6.2) |
| (QRR) | quoted; sign pinned | human reading of GJR (17) in the typeset PDF and of its sources; Teleman s.8.2 verbatim |
| (CYC) | PROVED, refereed twice; VERIFIED by six codes | human check (two pages) |
| formal SP lemmas | PROVED twice, refereed | human check of FORMAL_SP_LEMMAS s.4-s.7, s.9.3 |
| Theorem Gamma | PROVED twice, each proof refereed twice | human check of LMI (one page) |
| Theorem P(b) | PROVED (Stage 4, refereed; re-derived by pb, whose new steps s.2.6-s.2.9 are refereed once) | human check |
| the glue | PROVED, refereed once; refereed alternative EVID s.2.3-s.2.4 | - |

## 3. What Stage 5 established

### 3.1 The proof of (Gamma)
* **The 7-variable chart** (CHART_PROOF_3 s.3.1-s.3.6).
  * Replace the six Gamma factors by Euler integrals. This gives a formal integrand `Omega_7 = z^-1 e^(Phi/z) prod_j x_j^-1 d sigma dy
    dx_0..dx_4`; the powers of `z` cancel identically in `sigma`. Its critical point `p(L) = (L, 5L, zeta^j - L)` has
    `det H_7 = 1/(L^2 (1 - L^5)^2)`.
  * Fubini in the order `(y, x_j)` then `sigma` gives `c(z; L)`.
  * Fubini in the order `(sigma, x_0)` then the rest, with the exact delta lemma, gives a 5-variable chart with phase
    `psi = 5 log y - y + Sum_(j>=1) (-x_j + (zeta^j - 1) log x_j) + q y^5/prod x_j`. Its data are rational in `L` and regular at
    `L = 1`, because `1/det H_ch = -L^2 (1 + L + .. + L^4)^2` is a polynomial.
  * Two rational functions agreeing off `{0, 1, zeta^j}` are equal, so `c_n` has no pole at `L = 1`.
  * At `q = 0` the chart is a product of the Euler integrals of `Gamma(1 + 5/z)` and `Gamma((zeta^j - 1)/z)`. Their Stirling factors
    give exactly `h_(w_0)`; the weight `-5` comes from `Gamma(1 + 5/z)`.
* **The two-variable chart, Theorem LMI** (referee_adversary s.3).
  * Only the factor that degenerates at `L = 1`, `Gamma((1 - sigma)/z)`, needs an Euler integral. With `x = e^nu`,
    `nu = mu - lambda/alpha0` and `Q = e^(-lambda/alpha0)`, the phase becomes `g(sigma) + (1 - sigma) mu - Q e^mu`, where `g` is the sum
    of the phases of the five regular factors.
  * This phase is holomorphic in `Q` at `Q = 0`, with Hessian determinant `-1` there, and `Q = (1 - L) e^(-g'(L))` is a coordinate at
    `L = 1`.
  * At `Q = 0` the phase is linear in `mu` and the covariance has `C_ss = 0`. So the series is exactly the normalised product of the
    five regular Stirling factors at `sigma = 1`, which is `h(z)`.
  * In general (any regular phase and amplitude, and a singular factor `Gamma(M0 + alpha0 (sigma - s*)/z)` with `M0` in Z, or in C by
    the adversarial referee's extension): the saddle series, continued to the colour point, is the normalised residue of the integrand
    at the first pole `sigma = s* - M0 z/alpha0`.

### 3.2 Formal stationary phase
* Two independent written sets of lemmas: change of variables; Fubini, with the factor `beta = (det H_vv(w)/det H_vv(w_0))^(-1/2)`;
  parameters; the exact delta lemma; Stirling as the SP of Euler's integral with `x^m`, for every `m` in Z and both measures.
* Refereed: `referee_formal` s.2 (CHART_PROOF_2) and `lemmas/referee` s.1 (FORMAL_SP_LEMMAS). The lemmas referee's corrections C1-C8
  include the homogeneity Lemma H, which P(b) and the hatted form use; it is now stated and proved.
* No source read by lit states these formal lemmas for complex data with statement numbers. The project's proofs are self-contained.

### 3.3 Theorem P(b) and the glue, re-audited (pb; refereed by `referees/pb_steps`)
* Theorem P(b) was re-derived from scratch, with no mathematical gap. Three steps are written out for the first time:
  * (S2), the signs `eps_b`, previously VERIFIED only;
  * the relabelling of the residue, as one summable family summed in two orders;
  * a key lemma: the coefficient of `u^E z^b` involves only finitely many reduction terms. The referee added a weight lemma that puts
    the sum in `Q((z))[[u]]`; the chain does not need it.
* The glue is an exact identity, `D_0(z) = c(z; 1)/h(z)`. The referee proved formal WKB uniqueness in one line and asked for one
  sentence naming the ring of comparison. Lemma N's `K_q` formula is needed only for the finite route `g <= 69`.
* pb's own code: `Phi^SP = Phi^can` in every coefficient through `w^80`, at two primes and in two frames; `V = C0 (Y^C_Q)_0` exactly
  through `u^40`; the mechanism itself (CFU plus reduction reproduces the canonical gauge's first row) through `u^20`.

### 3.4 (CYC)
* **Statement.** In every conifold frame `PF4 = (1/varpi_0) D (1/Dt) D (u varpi_0^2 (Dt)^2) D (1/Dt) D (1/varpi_0)`, so the cyclic
  system with its closing relation `4 -> 0` is PF5. The referees widened the scope: any power-series solution `varpi_0` not
  proportional to the vanishing period `g` will do; `varpi_0(0) != 0` is used only downstream, for the unit shape that Theorem U needs.
* **Proof** (`critic/CYC_PROOF.md`):
  * PF4 is formally self-adjoint for `D = q d/dq`;
  * the Lagrange concomitant is a constant, non-degenerate, alternating form on the solutions;
  * transversality: `T_k = 0` for `k <= 2` and `T_3 = 1/u`;
  * the automorphism `log u -> log u + 1` of `Q((u))[log u]` makes every power-series solution orthogonal to `g`, so every conifold
    frame is isotropic;
  * special geometry then gives the factorisation, with Yukawa coupling `1/(u varpi_0^2 (Dt)^3)`.
* The two computations the proof leaned on are now derived by hand: `gamma = 7/10`, and the Pfaffian 2 of the form.
* The adversarial referee re-derived the cyclic system from the project's `fq_rows` and `fq_conifold`. It is exactly the statement
  the programme uses, in both directions of the Corollary.

### 3.5 What the mechanism needs (general: generalisation falsifier, 12 theories, 40 colours)
* The predictions held in every case tested, with the value `h` carrying the QRR weights (an Euler-class twist contributes minus its
  fibre weight, an inverse Euler class plus its fibre weight):
  * the O(5)-twisted `P^4` with four sets of distinct rational weights, every colour, exact to `z^24`;
  * `P^1`, `P^2`, local `P^2` (two lifts), `O(-1) + O(-1)` over `P^1` (two lifts), and `O(3)`, `O(4)`, `O(6)` over `P^4`, exact to
    `z^20`.
* Regularity is a property of the degenerating factor alone: it must be an exact Stirling form, and the only factor degenerating
  there.
  * Coincident weights give poles of order exactly `n` at `z^n`.
  * Two singular factors of different slope give poles.
  * The regular factors can be arbitrary; they enter only the value. This is what LMI says.

### 3.6 Literature (lit) and corrections to Stage 4
* Prior art:
  * Givental arXiv:math/0108100 s.10: untwisted `CP^n`, the chart `x_0 = Q/(x_1..x_n)`, Lemma 10.4. It is an outline, and continuity
    in `Q` is asserted.
  * GJR arXiv:1812.11908: Remark 7.2 states the analogue of (Gamma3) for their LG integral without proof, and Lemma 7.7, eq. (37)
    proves the orbifold-point formula (`Ginf`) for it.
  * Iritani arXiv:0903.1463v3 Lemma 4.20: an exact `q`-expansion at `z = 1`.
  * No paper read proves (Gamma). The write-up should cite these and claim no more novelty.
* Corrections to `MECHANISM_STAGE4.md` (left unchanged):
  * s.2: the (QRR) source (O11);
  * s.2: the P(b) hypothesis list is incomplete; add (CYC), (S2), the relabelling, L3(iv), Lemma T(b) and Lemma H (O10);
  * s.1-s.2: the status of (Gamma) (O18).
  * Full list: STAGE5_SYNTHESIS_2 s.6.

### 3.7 The referee round (second workflow)
* Verdicts: `cyc_formal` and `cyc_adversary` sound with corrections; `pb_steps` sound with corrections; `lmi_formal` sound with
  corrections; `lmi_adversary` sound. None found a mathematical error or a gap. The 38 numbered corrections concern justification,
  wording, labels and bookkeeping.
* The critic settled four disagreements. A1: the wider scope of (CYC). A3: LMI holds for non-integer `M0`. A4: who supplies log-free
  rows; the glue derives them, so they are not a new input. A6: finite fibres alone give only a formal sum; the weight lemma completes
  the step.
* Two referees independently found the same mislabelled frames in the critic's checks (`y_0 = f_0 + a f_1 + b f_2`, not the
  programme's `f_0 + alpha g + beta f_2`; harmless) and the same formal proof of DLMF 5.11.8.
* New statements made by the referees, not used by the chain: the weight lemma; a formal DLMF 5.11.8, which makes `Ginf` purely
  formal; LMI for every `M0` in C; (CYC) for every Lagrangian pair.

## 4. My own checks (2026-11-14)

### 4.1 Re-runs of the agents' code
Every script was run from a scratch copy (the scripts write their logs next to themselves), with `py -I -B` and the stored
arguments. Each new log was compared line by line with the stored one, timings stripped.

**First workflow: all 42 identical.** That is 41 logs, plus the stdout of `lit_vertex_check.py` (which writes no log of its own)
compared with its stored log.

| agent | scripts and arguments | what the decisive ones check |
|---|---|---|
| chart (8) | `ch_chart 3 4 6 5 7` (39 min), `ch_hess`, `ch_zpow`, `ch_stirling 30`, `ch_sigma 10`, `ch_eexp 6 10`, `ch_lemmas 6 4`, `ch_euler 40 40` | chart = `h` at `L = 1` exactly to `z^7`; chart and 7-variable series = `c(z; L)` at several `L`; Hessians in both orders; (L3) to `z^40` |
| referee_formal (6) | `rf_sigma 40 10`, `rf_chart 3 2 4 3`, `rf_hat 30`, `rf_hess`, `rf_eexp 6 12`, `rf_lemmas 7 6` | (Gamma) exactly to `z^40`, by Laurent expansion at `L = 1 + eps` |
| referee_adversary (4) | `ra_run 18 24`, `ra_2var 9`, `ra_l2src 8 8`, `ra_poly 8 6` | LMI's predictions and falsifiers to `z^18`; the two-variable chart to `z^9` |
| lemmas (6) | `chk_fubini 8 3`, `chk_delta 8 2`, `chk_change 8 2`, `chk_param 8`, `chk_res 6 14`, `chk_stirling 80 12` | the formal lemmas on random data |
| lemmas/referee (3) | `rl_toy 6 14 52`, `rl_l3 40 16 8`, `rl_fl4 8 8` | a 4-variable toy chart; (L3) and its formal proof |
| general (6) | `g5_t1 60 60 Q 6`, `g5_t2 24 10`, `g5_t3 20`, `g5_t6 20`, `g5_t6c 20`, `g5_valid 6 24` | (Gamma) exactly to `z^60`; twisted `P^4`, all colours; 12 theories, 40 colours |
| pb (7) | `pb_run_phi 80 72`, `pb_res 40`, `pb_cyc 60`, `pb_glue`, `pb_s2 60`, `pb_exact 15`, `pb_mech 20` | `Phi^SP = Phi^can` to `w^80` at two primes; `V = C0 (Y^C_Q)_0` exactly to `u^40`; (CYC) to `u^56` |
| critic (1) | `cr5_cyc_run 80` | (CYC), lemma by lemma |
| lit (1) | `lit_vertex_check` | `N_(2,0) = -19/9600`, `N_(3,0) = 1/224000` from LP (26); the weight `+5` control |

My first launches of `g5_t1 60 60 Q` and `g5_t2 24` used the default control ranges (`z^12`, `z^8`) instead of the stored ones
(`z^6`, `z^10`). Every main-result line was identical; the control lines differed. Re-run with the stored ranges, both are identical.

**Referee round: 46 final logs re-run, 44 identical, and the other 2 have the same verdicts.**
* The runs covered all five referees and their critic:
  * `cyc_formal`: 6 scripts, to `u^80`;
  * `cyc_adversary`: 13 logs, including both directions of the Corollary on the project's own rows;
  * `pb_steps`: the four final `rps_run` logs, to `E = 120`;
  * `lmi_formal`: 13 logs, including the chart at `L = 1` against `h` to `z^60`;
  * `lmi_adversary`: 7 logs, 219 checks;
  * the critic: `rc_cyc 40`, `rc_cyc 60`, `rc_l3c 60 12`.
* Two scripts import package-root modules by a relative path (`cy_proj.py`, `cyc_ref_xproj.py`). I ran them in a scratch mirror of the
  package layout.
* The six `cyc_formal` scripts print to stdout. I compared their output with the stored logs, without the logs' header line.
* The two remaining stored logs come from earlier versions of their own scripts:
  * in `cyc_adversary/logs/cy_rows_30_5.log`, four control lines lack the "fails at `z^[..]`" detail that the final script prints;
    the verdicts are the same;
  * `pb_steps/logs/rps_run_v3_40_10_8_24.log` predates the script's T5 section; its first 56 lines are identical, and the three other
    `v3` logs, which include T5, are identical.
* Not re-run: the two kept draft logs (`rps_run_v2_*`, `DRAFT1_*`) and the superseded first-version `rps_run` logs.

### 4.2 Two independent checks (`stage5/independent_check/`, new code sharing nothing with the agents')
* `cyc_check.py 60` (11 s):
  * builds the Frobenius basis of PF4 at `u = 0` from the recursion;
  * finds the vanishing period from the `log u` obstruction, `g = f_1 + (7/10) f_2`, with nothing assumed;
  * checks that the factorised operator equals PF4 on `u^0..u^5` and kills `f_0, f_1, f_2`, through `u^60`, in four frames
    `y_0 = f_0 + a f_1 + b f_2` and with `g` rescaled;
  * controls (the non-isotropic pairs `t = f_1/y_0` and `t = (f_1 + f_2/2)/y_0`, and `g + f_2/7`) fail at `u^0`.
* `lmi_check.py 16` (201 s), exact over Q (all symmetric functions of `zeta^1..zeta^4` come from `1 + x + .. + x^4`):
  * the two-variable chart series equals the one-variable full-Stirling `c(z; L)` at `L = 2, -2/3, 1/2, 3` through `z^16`;
  * at `L = 1` it equals `h_(w_0)`, with `h` computed from its definition by Newton sums;
  * `c_1..c_3(2) = -12/5, 72/25, 354278/625`, as in the chart worker's log;
  * controls fail: the singular factor's Stirling series dropped; twist `+5`; `Q e^(2 mu)` for `Q e^mu`.

### 4.3 My reading
* CHART_PROOF s.0-s.3 in full, with CHART_PROOF_2 s.2 (all lemmas). Checked by hand:
  * the critical point and `det H_7`;
  * the order-(a) amplitude: `beta` times the six inner Stirling series is `A_St` up to a constant;
  * the delta lemma's exactness: every vertex has `xi`-degree at least twice its `s`-degree, and `C_xixi = 0`;
  * the chart phase `psi`;
  * the rationality of `Nch_n` and the continuation;
  * the value at `q = 0` and the weight `-5`.
* referee_adversary s.3 (LMI and (L2')). Checked by hand:
  * the chart phase after `nu = mu - lambda/alpha0`;
  * `Q(L)` is holomorphic with a simple zero at `s*`;
  * the Wick sum of (L2');
  * the quintic case `M0 = 0`, `delta = 0`, with value `F(1, z)/F(1, 0) = h`.
* critic/CYC_PROOF.md in full: the adjoint computation, the concomitant, isotropy via `log u -> log u + 1`, special geometry, the
  final uniqueness step.
* FORMAL_SP_LEMMAS s.9.2-s.9.3 (the formal proof of (L3)): the functional equation via `D^2 rho = 0`, and the `a`-shift with
  `eps = z`.
* lit s.1 ((QRR)).
* pb s.0-s.3: s.2.6-s.2.9 read but not re-derived line by line; the glue checked.
* The referee round: the five verdicts and the critic's synthesis, s.0-s.2.
* Not read by me: the rest of FORMAL_SP_LEMMAS (only through its referee), GENERAL, LIT_2 and LIT_3, and the cited papers.

## 5. What is left
1. **Human checks**, in this order (REFEREES_SYNTHESIS_3 s.8):
   1. Theorem Gamma by LMI (one page), with the lemmas it uses;
   2. (CYC) (two pages, with the hand computations of `gamma = 7/10` and the Pfaffian);
   3. Theorem P(b), above all pb s.2.6-s.2.7;
   4. the glue;
   5. (QRR): GJR (17) in the typeset PDF, its sources, and the uniqueness step (Teleman s.8.2);
   6. Theorem U and GI Theorem G;
   7. the Stage-2 list.
2. **Write-up** in a separate new draft (`allgenus/note/`, reusing `refs.bib` if wanted; the draft with Yannik stays untouched):
   * the polynomial structure;
   * Theorem P and the glue;
   * Theorem Gamma by LMI, with the 7-variable chart as a remark;
   * `Ginf` via the formal DLMF 5.11.8;
   * (CYC) as a set-up lemma;
   * an appendix on formal stationary phase.
3. **Optional**:
   * a second AI referee for pb s.2.6-s.3, the only chain steps refereed once (a human check is worth more);
   * the genus-3 code test of (QRR): `F_3(0) = 1/224000` from the project's R-matrix graph sum at `u = 1`;
   * Teleman s.8.2 verbatim, and GJR (17) in the typeset PDF.
4. **Housekeeping, not done**: corrections to the worker notes.
   * Stage 5: RF C1-C6, RL C1-C8, O1-O19, and the 38 corrections of REFEREES_SYNTHESIS_2 s.5.
   * Stage 4: O1-O18. Stage 2: O1-O25.
   * Compliance points reported by the agents: pb and the adversary overwrote or deleted logs of their own; `pb/pb_core.py`
     (286 lines) and `lemmas/referee/rl_core.py` (307 lines) exceed 250 lines; `referees/lmi_adversary` once created, and then
     deleted, an empty file in `/tmp`.
