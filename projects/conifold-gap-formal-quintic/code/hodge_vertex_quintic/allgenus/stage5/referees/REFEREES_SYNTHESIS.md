# Stage 5 referees: synthesis (critic of the five referees), 2026-11-14

Critic of the referee round `allgenus/stage5/referees/`. Five referees: `cyc_formal` and `cyc_adversary` (object: the Stage-5 critic's
`critic/CYC_PROOF.md`), `pb_steps` (object: `pb/PB_AUDIT.md` s.2.6-s.2.9 and s.3), `lmi_formal` and `lmi_adversary` (object:
`chart/referee_adversary/REFEREE.md` s.3, Theorem LMI and (L2')). AI-written; no human has read any of the referee reports or this
synthesis. Labels as in STAGE5_SYNTHESIS: PROVED, VERIFIED (range, field or primes stated), CONJECTURAL, SKETCH, REFUTED.
Abbreviations: SYN, SYN2 = `../STAGE5_SYNTHESIS.md`, `_2.md`; CYC, PB, RA, CP, FSL, RL, RF, GEN, UR, GI, POLY, EVID as in SYN;
CF = `cyc_formal/REFEREE.md` (+`_2`, `_3`), CA = `cyc_adversary/REFEREE.md` (+`_2`), PS = `pb_steps/REFEREE.md` (+`_2`),
LF = `lmi_formal/REFEREE.md` (+`_2`), LA = `lmi_adversary/REFEREE.md` (+`_2`). Corrections keep their authors' labels, prefixed:
CF-E1..E8, CA-K1..K7, PS-C1..C10, LF-C1..C6, LA-C1..C7 (CA-K2 is not SYN2's K2). Part 2 (`REFEREES_SYNTHESIS_2.md`): s.5 all
corrections, s.6 final status. Part 3 (`REFEREES_SYNTHESIS_3.md`): s.7 the critic's checks, s.8 human checks and next steps, s.9
record. My own code and logs: `critic/`.

**Read in full:** the eleven referee notes; SYN and SYN2; the three objects (CYC_PROOF.md, PB_AUDIT.md, RA s.3); UR s.0-s.1 and
s.3-s.4; GI s.4; STAGE4_SYNTHESIS s.1.2; the chain of STAGE2_SYNTHESIS s.2; FORMAL_QUINTIC_genus2 s.1; the docstring of `fq_ring.py`
(the `i_3` relation); by search, the import and log-writing code of the 44 referee scripts, and the headers of the two that import
project code. All 48 referee logs were scanned: 46 end in ALL CHECKS PASSED; the other two are kept drafts, each reported by its author
(PS `rps_run_v2_40_10_8_24.log`, a wrong prediction of the referee; LA `DRAFT1_lm_run_x_4.log`, a removed sub-claim); every other
FAIL-type line is a labelled control or falsifier; no log prints an absolute path.

## 0. Bottom line

* **No referee found a mathematical error or a gap in any of the three objects.** Verdicts: CF sound-with-corrections, CA
  sound-with-corrections, PS sound with corrections, LF sound-with-corrections, LA sound. All 38 numbered corrections concern
  justification, wording, labels or bookkeeping; none changes a statement, a hypothesis list or a label of the chain (part 2, s.5).
* **Every "(u)" item of SYN s.0 is now refereed (AI):** (CYC) twice (CF, CA); PB s.2.6-s.2.9 and the glue PB s.3 once (PS);
  Theorem LMI with (L2') twice (LF, LA); the by-product `c(z; 0) = Ginf(z)` twice (LF, LA); `c_n in Q[L]` once (LF). s.4 lists every
  occurrence of "(u)" in SYN with its new status.
* **No step of the chain S1-S9 of SYN s.2 is unrefereed.** Steps refereed by ONE agent only: PB s.2.6-s.2.9 and s.3 (PS), on top of
  their Stage-4 counterparts (POLY s.6 with the poly referee; EVID s.2.3-2.4 with its referee). No step has been checked by a human.
* **(Gamma): PROVED**, by two independent proofs, each refereed twice: CP s.3 (seven-variable chart; RF, RA) and RA s.3 (Theorem LMI;
  LF, LA). Inputs: the formal stationary-phase lemmas (S1) only.
* **`D_0 = 1` in `Q[[z]]` to all orders: PROVED modulo (U) and (QRR).** **Conjecture N(g) for every `g >= 2`: PROVED modulo (U), (Q),
  (QRR).** The inputs are those of SYN s.1, unchanged; this round adds no hypothesis. Two bookkeeping items must enter the (CYC) row
  (C4 also makes the project's conifold map `Cf` well defined, CA-K5; EVID's Lemma N builds rows by the cyclic system, PS-C9), and one
  sentence must enter the glue (the ring of comparison, PS-C8). Neither is a new input (s.2, A4-A5).
* **Disagreements (s.2): four, all settled.** A1, scope of (CYC): CA is right, `y_0(0) != 0` is never used (the narrower statement of
  CF and SYN2 s.9 is correct but weaker). A3, non-integer `M0`: LA is right, LF's SKETCH label is superseded by LA's interpolation
  proof. A4, the role of log-free rows: CF and CA are both right, about different statements. A6, SYN2's K5 (iv): PS is right that
  finite fibres give only a formal sum; K5's conclusion stands. Two referees independently found the same mislabelled frames in SYN
  and CYC (CF-E3 = CA-K1; A2) and the same formal proof of DLMF 5.11.8 (LF s.4.2 = LA Lemma F5.11.8; A7).
* **New statements made by the referees (not used by the chain):** the weight lemma (PS s.2.3; checked here by hand), the formal DLMF
  5.11.8 (LF, LA), Lemma L3-C and LMI for every `M0` in C (LA s.2.3), (CYC) for every Lagrangian pair in `S cap K` (CA-K2). I
  re-proved L3-C and 5.11.8 by hand and VERIFIED them with new code (`critic/rc_l3c.py`); status in part 2, s.6.3.
* **Re-runs (part 3, s.7).** 26 referee logs (for every recorded check the stored log of highest range; lower-range duplicates,
  superseded first versions and the two drafts were not re-run) re-run from scratch copies: all identical to the stored logs once
  timings are stripped. New exact code `critic/rc_cyc.py` re-verifies (CYC) by a
  method none of the five CYC checks used (operator identities tested on `u^x` at six non-integer `x`), the full 3x3 obstruction
  matrix, the hand values of `Bmat`, `det Wr = -2/u^2`, the frame relabelling, the frames with `y_0(0) = 0`, and the Corollary in
  operator form; every control fails.

## 1. The five reports at a glance

| referee | object | verdict | numbered corrections | own exact checks (each with failing controls) |
|---|---|---|---|---|
| cyc_formal | CYC s.0-s.9 | sound-with-corrections | E1-E8 | `M = PF4` below `u^77` (5 frames), identically in `(alpha, beta)` below `u^27`; all 16 entries of C3; `det Wr = -2/u^2`; project code cross-check |
| cyc_adversary | CYC and its uses | sound-with-corrections | K1-K7 | 17 frames to `u^74`; both directions of the Corollary on actual objects (WKB rows, `z^0..z^12`; `Y^C_Q`, `n <= 16`); `Cf(kinv) = 1/c_a` mod `2^89-1`; deformations |
| pb_steps | PB s.2.6-s.2.9, s.3 | sound, with corrections | C1-C10 | relabelling (965 hatted coefficients), mechanism (4283), key lemma (884 keys), (S2) to `t^80`, the s.2.6 identity gauge-free (`w^25 u^25`) |
| lmi_formal | RA s.3 (LMI, (L2'), consequences) | sound-with-corrections | C1-C6 | chart at `L = 1` = `h` to `z^60`; global fit `n <= 30`; Method A `z^14`; formal 5.11.8 to `w^-60`; `Ginf` to `z^80` |
| lmi_adversary | RA s.3, adversarially | sound | C1-C7 | 219 checks: exact `z^10`/`z^12`, mod two new primes `z^12`-`z^16`; 27 falsifiers fail as predicted |

All five worked from the definitions with new code; none imported the code of the object it refereed (CF and CA imported read-only
scratch copies or the package copies of `fq_conifold`, `fq_ring`, `fq_rows`, `fq_fp`, `fq_series` for the cross-checks only).

## 2. Disagreements and points needing a ruling (which side is right, and why)

* **A1. Scope of (CYC).** CF-E5 and SYN2 s.9: "every conifold frame (power-series `y_0` with `y_0(0) != 0`, every scale of `g`)"; CF
  REFEREE_2 s.7 calls `y_0` and `Dt` "units of `Q[[u]]`". CA-K2: `y_0(0) != 0` is never used; the Theorem holds for every non-zero
  power-series solution `y_0` not proportional to `g`, and for every pair spanning a B-Lagrangian plane of `S cap K` (the planes that
  contain `g`). **Verdict: CA.** Every step of CYC s.1-s.7 takes place in `K = Q((u))` or `K[l]`, where every non-zero element is a
  unit; C5 needs `y_0`, `Dt` invertible in K only, and s.7 needs only `det Wr != 0`. VERIFIED here (`rc_cyc 60`, E1-E2):
  `M = PF4` below `u^69` for `y_0 = f_1`, `f_2`, `f_1 - 2 f_2` and the pair `(f_1, f_2)`. Both agree that the large-radius frame and the
  ring N are not covered and not needed (CF REFEREE_2 s.9). What `y_0(0) != 0` buys is downstream: the Theorem-U unit shape
  `val(c) = (0,-1,0,0,0)` (`(0,-3,0,1,1)` for `y_0 = f_1`, `(-2,-1,-2,2,2)` for `f_2`: CA, and `rc_cyc` E3), and the normalisations
  `C0(0) = 1`, `t_c'(0) = 1` of PB s.2.8-s.2.9, POLY s.3 and GI Lemma Q(d), which hold in the fq_conifold family only. Nothing breaks:
  the programme uses only that family.
* **A2. The critic's frame labels.** CF-E3 and CA-K1 independently: the frames `(1/3,-2/7)`, `(7/11,-3/5)` of `cr5_cyc_run`, CYC s.8,
  SYN s.1 (row (CYC), "3 frames") and SYN s.4.1 are `y_0 = f_0 + a f_1 + b f_2`, not the programme's `f_0 + alpha g + beta f_2`. **Both
  right, and they agree:** critic `(1/3,-2/7)` = programme `(1/3,-109/210)`, critic `(7/11,-3/5)` = programme `(7/11,-23/22)`, programme
  `(1/3,-2/7)` = `f_0 + f_1/3 - (11/210) f_2` (VERIFIED, `rc_cyc` E0). Harmless (the Theorem covers both families); relabel.
* **A3. Non-integer `M0` in Theorem LMI.** LF-C6: RA_2 s.5's "for non-integer `M0` ... same proof" is a SKETCH (VERIFIED for
  `h = 1/5, 1/3, -1/2, 7/4` to `z^14`). LA s.2.3: PROVED, by Lemma F5.11.8 (formal DLMF 5.11.8) and Lemma L3-C ((L3) for every `M` in C).
  **Verdict: LA.** LF's label was right for RA_2's sentence; LA supplies the missing lemma. The proof of L3-C is an interpolation: by
  (F1) each `[z^n]` of `N[log t - t, t^(M-1); 1]` is a polynomial in the Taylor coefficients `binom(M-1, j)` (`j <= 2n`), hence in `M`;
  so is `[z^n] exp(Sum beta_n(M) z^n)`; they agree on Z (FSL4 s.9.1(i), (iv) and F5.11.8 at `h = m`). I checked each step by hand and
  VERIFIED the identity as polynomials in `M` through `z^12` with a symbolic-`M` Wick sum (`critic/rc_l3c.py`, L1-L2). Not used by
  the chain (the quintic needs `M0 = 0` at `L = 1` and `M0 = 1` at `L = 0`). Status: PROVED (AI; written by LA, checked by the critic
  only).
* **A4. Is the Laurent-polynomiality of the actual rows an input?** CF REFEREE_2 s.8.2: for the closing relation, "whether the `psi_n`
  are Laurent polynomials in L plays no role". CA-K4: "the closing relation of the quintic rows is a theorem" presupposes WKB solutions
  of PF5 with Laurent-polynomial coefficients, supplied at all orders by P(a) and POLY Lemma 2, not by (CYC). **Verdict: both right,
  about different statements.** The Corollary needs only a PF5 solution in a differential module (CF-E8, CA-K4 name the rings). That
  the actual rows exist as `s`-series near `u = 0` (continuation from `q = 0`) is a separate fact. In the all-orders chain it is not an
  extra input: PS s.4.2 (a)-(d) derives it inside the glue. `Phi^can` has polynomial coefficients (Theorem U, (CYC), Lemma 2); formal WKB
  uniqueness in `Q[[q]]` gives `A^QRR = A^can/d`, so `A^QRR_n in Q[L, 1/L]`; then UR s.4 is compared in `Q[L,1/L][[z]] c Q((s))[[w]]`.
  I re-checked (a)-(e) of PS s.4.2 line by line. The write-up must state it (PS-C8).
* **A5. C4 and the conifold map `Cf` (CA-K5).** `fq_ring` eliminates `i_3` by `B(c_0, c_0 Lam + i_0) = 0`; its `Cf` image is
  `B(varpi_0, g) = 0`, i.e. C4. CF REFEREE_2 s.9 says the same ("`B(varpi_0, g) = 0` holds by C4"). **Verdict: correct, no conflict.**
  C4 PROVES what FORMAL_QUINTIC_genus2 s.1 recorded as checked ("Cf preserves the Wronskian relation"; `pairing_vacuous`). It is used by
  every computation that maps N-rows to the conifold frame (all finite-genus evidence, `fq_rows` via `Cf`). The abstract chain S1-S9
  works with the WKB columns directly (A4) and does not pass through N. C4 is part of the twice-refereed CYC: record it, no new input.
* **A6. SYN2 K5 (iv)** ("finitely many `pi`-terms per target monomial puts the hatted residue sum in `Q((z))[[u]]`"). PS-C5: finite
  fibres give only a formal sum; `E >= 0` needs a bound on the w-degree of `pi_m`, the weight lemma `5i - m - 4j >= 0` (PS s.2.3).
  **Verdict: PS.** I checked the weight lemma by hand: monomials of `B` have weight 0 (`x^a z^b -> u^((a+4b)/5) X^a w^b`); the CFU map
  and its inverse are `X(1 + Sum c u^k T^d)` with `d <= k`, weight `>= 0`; the two reduction steps keep or raise the weight by `5(l+1)`.
  K5's conclusion is unaffected and the chain does not need the lemma: the identity is coefficientwise and `V in Q[z^-1][[u]]`
  (PB s.2.1), so membership follows from the identity. K5 (iii) needs PS-C4's two sentences.
* **A7. Formal DLMF 5.11.8.** LF s.4.2 and LA Lemma F5.11.8 give the same proof independently (difference `D_n(h)` of the coefficients
  is a polynomial with `D_n(h + 1) = D_n(h)`, `D_n(0) = 0`); RA's route (uniqueness of asymptotic expansions) is also valid (LF).
  **Agreement.** Inputs: Stirling's formal functional equation (FSL4 s.9.1(iv)) and `B_k(h+1) - B_k(h) = k h^(k-1)` (DLMF 24.4.1, which
  LA derives from 24.2.3). VERIFIED here as polynomials in `h` through `w^-60` (`rc_l3c`, D1-D3). So `c(z; 0) = Ginf(z)` is purely formal.
* **A8. LMI's input list and the (Gamma) bullet.** LF-C1-C4 and LA-C1-C4 overlap: identification of the critical points written out
  (LF-C1; LA s.1.3); `Q^M0` drops out by `Nhat` (LF-C2); the punctured disc (LF-C3 = LA-C3); `F` must contain `(5 sigma)^(1/2)` (LF-C4 =
  LA-C4); `Omega_2` is a formal integrand (LA-C1; LF s.3.1); inputs (F1)-(F2), the translation case of (L1a), CP s.1.4 and the
  `C[[z]]`-linearity of N (LA-C2). **No conflict**; merged in part 2, s.5.
* **A9. Labels of LMI.** GEN s.3.3: CONJECTURAL, SKETCH; SYN2 K6: PROVED (u); LF, LA: PROVED. **Verdict: PROVED**, refereed twice; GEN's
  statement is RA's once `Lam(0)` is read as the constant term of `Lam(z)` (LA s.5). LF and LA differ only in the verdict word.
* **A10. Smaller points, no conflict.** PS-C1 (PB's FL3 is also used in (S2)) and PS-C6 (Theorem G for `P^SP` is redundant) are right.
  LF's clarification (RA's "without Theorem P(a)" covers `c_n in Q[L]` only; no negative powers of `u` in `ph^SP_n` still needs P(a))
  agrees with SYN s.3. LA's remark that (Gamma) is blind to the twist measure (`Gamma(1+5s)` vs `Gamma(5s)`) agrees with RF's remark
  (SYN2 s.8 item 3). CA-K3 (a wrong `gamma` or a non-isotropic pair still has the Theorem-U shape) is a robustness remark: the
  identification rests on `gamma = 7/10` exactly, which is PROVED (CYC s.1 with CF-E1; `rc_cyc` B3, E4).

## 3. The input table of SYN s.1, row by row, after the referees

| input (SYN s.1) | status in SYN | status now | what would close it |
|---|---|---|---|
| (U) | PROVED (UR s.3), AI-audited; VERIFIED `z^80` (PB) | unchanged (no referee of UR s.3 in this round). PS s.3 re-derived U's hypotheses for `(A, At)` (`Prod b = Prod bt` iff `v^(4/5) = u^(4/5)/tau`); CA-K2: the unit shape holds iff `y_0(0) != 0`, true in every fq_conifold frame | human check of UR s.3 |
| (Q) | quoted with statement numbers | unchanged (not in this round) | human check (SYN2b s.6.2) |
| (QRR) | quoted; sign pinned | unchanged. PS s.4.3: the glue uses exactly H9, `A^QRR_n(q = 0) = h_n`; the passage `R(q = 0) = diag(h_(w_i))` to H9 belongs to (QRR) and was not re-checked | human reading of GJR (17), its sources, Teleman s.8.2 |
| (CYC) | PROVED (CYC, (u)); VERIFIED (pb_cyc `u^56`, cr5 `u^70`, "3 frames") | **PROVED, AI-refereed twice** (CF, CA: no error, no gap; corrections CF-E1..E8, CA-K1..K7). The two computed inputs are now hand-derived: `gamma = 7/10` (CF s.2.3; CA s.3) and `Pf(Bmat) = -B(f_0,f_2) B(f_1,y_3) = 2` (CF s.4.4). VERIFIED five times by different methods: pb_cyc (`u^56`, programme frames), cr5_cyc_run (`u^70`; its frames are programme `(0,0)`, `(1/3,-109/210)`, `(7/11,-23/22)`), CF (`u^77`, 5 frames; generic `(alpha,beta)` to `u^27`), CA (17 frames, `u^74`; both directions on actual objects), critic `rc_cyc` (`u^70`, 10 frames incl. `y_0(0) = 0`); also the main session's `independent_check/cyc_check` (`u^60`, 4 frames). Used in: P(b) step 5; the glue; POLY Lemma 2 for `Phi^can`; UR s.1/s.4 (closing relation of the actual rows), hence R0 and the Stage-2 chain; PB's `Q^SP`; EVID Lemma N (PS-C9); C4 also makes `Cf` well defined on N (CA-K5) | human check of CYC_PROOF (two pages) with CF s.2.3, s.4.4 |
| formal SP lemmas (S1), (PG), (L0), (L1a)-(L1c), (L2), (L3) | PROVED twice, refereed | unchanged. LF and LA checked that every use in LMI meets the hypotheses of (L1a), (L1b), (L1c)(i), (L3), (F1)-(F2), CP s.1.4 (statements, not proofs) | human check of FSL s.5-s.7, s.9.3 |
| homogeneity (F3) = Lemma H; reflection | PROVED | unchanged; PS-C3: PB s.2.6 should cite L1a and Lemma H (two places) | - |
| residue FL2, FL3; relabelling `V = s Res_X` | PROVED (FSL refereed; PB s.2.7 (u)) | **PROVED, PB s.2.7 refereed (PS)**: complete with PS-C4's two sentences; key-lemma bound proved and attained; membership in `Q((z))[[u]]` by the weight lemma (PS-C5; not needed). VERIFIED exactly (PS T1-T3: 965 hatted coefficients, 4283, 884 keys) | human check of PB s.2.7 with PS s.2.2-s.2.3 |
| Theorem Gamma (i)-(iv) | PROVED (CP s.3; RF, RA) | **PROVED; second proof RA s.3 (LMI, (L2')) refereed twice** (LF, LA) | human check of RA s.3 (one page) with LF s.3 |
| hatted/flat identity | PROVED (CP s.3.10, RF C2) | unchanged. LF s.4.3: with `c_n in Q[L]` it bounds the top `u`-degree of `ph^SP_n` by `floor(4n/5)`; "no negative powers of `u`" still needs P(a) | - |
| Theorem P(b) | PROVED (POLY s.6 + poly referee; PB s.2 (u)) | **PROVED; PB s.2.6-s.2.9 refereed (PS, corrections C1-C7)**: (S2) PROVED (uniqueness half of UR Prop. S, pivot `4n`, independent of (H)); signs `eps_b`, weights, triangular uniqueness, `kappa = C0` confirmed; gauge-free identity VERIFIED `w^25 u^25`. Inputs as SYN, plus Lemma H; GI Theorem G for `P^SP` is redundant (PS-C6). PB s.2.2 (P(a)), s.2.4, s.2.5 not re-checked by PS: they re-derive POLY s.6 (poly referee) and P(a) (RL s.3.4) | human check of PB s.2.5-s.2.9 |
| glue `D_0 = c(z;1)/h(z)` | PROVED (PB s.3 (u); refereed alternative EVID) | **PROVED, refereed (PS, C8-C10)**: ring of comparison `Q[L,1/L][[z]] c Q((s))[[w]]`; formal WKB uniqueness PROVED in one line (simple root, transport kernel `Q L`); `A^QRR_n in Q[L,1/L]` as a by-product; EVID's and POLY's `Phi^can` are literally one element (given Lemma N, which needs (U), (CYC)) | human check of PB s.3 with PS s.4.2 |
| GI Theorem G, Lemmas 1, V, Q | PROVED, AI-refereed | unchanged; CA re-VERIFIED Lemma Q(a), (d) in three frames (`cy_rows`); PS-C7 states the used fact `[u^0](Y^C_Q)_0 = 0` | human check of GI s.3 |
| Stage-2 chain, steps 1-8 | PROVED, AI-refereed | unchanged (not in this round) | SYN2b s.6.2 |

## 4. The "(u)" items of SYN, every occurrence

| item | where in SYN | referee(s) now | new status |
|---|---|---|---|
| (CYC) | s.0 bullets 3-4; s.1 table; s.2 S9; s.3; SYN2 K4, s.9 | CF, CA | PROVED, refereed twice; scope wider than stated (A1) |
| (S2), signs `eps_b`, weights, `kappa = C0` (PB s.2.6) | s.0 bullet 4; s.2 S4 step 4; s.3; SYN2 K12 | PS | PROVED, refereed once |
| relabelling `V = s Res_X` and key lemma (PB s.2.7) | s.0 bullet 4; s.1 table; s.2 S4 step 6; s.3; SYN2 K5 | PS | PROVED, refereed once (K5 (iv) corrected, A6) |
| PB s.2.8-s.2.9 (regularity, Theorem G for `P^SP`, step 7) | s.1 table (P(b), "re-derived PB s.2 (u)") | PS | PROVED, refereed once (C6, C7) |
| the glue PB s.3 | s.0 bullet 4; s.1 table; s.2 S5 | PS | PROVED, refereed once (C8-C10) |
| Theorem LMI and (L2') (RA s.3) | s.0 bullets 1, 4; s.2 S2; SYN2 K6 | LF, LA | PROVED, refereed twice |
| `c(z; 0) = Ginf(z)` (RA s.3.4) | s.0 last bullet; s.3; SYN2 K7 | LF, LA | PROVED, refereed twice; purely formal (A7) |
| `f^can_m = [z^(5m)] Ginf` given P(b) | s.0 last bullet; s.3 | LF s.4.3 | PROVED given P(b), refereed once |
| `c_n in Q[L]` without the degree bound (RA R4) | s.3 | LF s.4.3 | PROVED (LMI at `s* = 0` with Theorem Gamma (i), (iv)), refereed once |
