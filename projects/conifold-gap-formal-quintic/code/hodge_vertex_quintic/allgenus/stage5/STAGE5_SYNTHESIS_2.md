# Stage 5 synthesis, part 2: contradictions, overclaims, human checks, next steps, critic record (2026-11-14)

Continuation of `STAGE5_SYNTHESIS.md` (s.0 bottom line, s.1 the assembled theorem and the input table, s.2 the proof chain S1-S9,
s.3 status verdicts, s.4 the critic's checks). Same labels, notation and abbreviations. AI-written; no human has read it.

## 5. Contradictions between notes, and which side is right

* **K1. Status of (Gamma).** CP s.0.5, RF s.0 and RA s.0: PROVED. GEN s.0 and s.7 (S1): "CONJECTURAL with a SKETCH proof (SYN4 s.1.4),
  unrefereed"; MECHANISM_STAGE4 s.2: CONJECTURAL. **Verdict:** PROVED (AI-written, refereed twice). GEN worked in parallel from the
  Stage-4 outline and never saw CP; its label describes its input. Its computations support the proof, they do not contradict it.
* **K2. Can (L3) be proved formally?** CP s.2.8 and s.5 item 4: "not done; the shift of `a` by `z` is not a formal operation"; RA R3:
  "(L1c)(iii) would supply it". FSL s.9.3 (Proof A): FL1 gives the recursion in `m`, L1c(iv) realises `a -> a + z` as the perturbation
  `eps log(x/a)` with `eps = z` legitimate by the weight bound, and a difference equation in `1/Y` closes it; RL s.2: "complete and purely
  formal"; RF C6: the Pochhammer part is formal by (S1). **Verdict:** FSL and RL. (L3) needs no analytic input; CP's remark undersells.
* **K3. Is Lemma N needed for `D_0 = 1` to all orders?** CP s.0.3 (R2 + R3 among the hypotheses of Theorem D) and RF s.5 ("R3 with R2"):
  yes. PB s.3 and P7: no; MECHANISM_STAGE4 s.5 item 2 already said "now needed only for g <= 69". **Verdict:** both chains are valid. The
  all-orders statement uses Lemma N's steps 1-5 (EVID s.2.3: `P_can Q^m` is a WKB solution; formal WKB uniqueness), not its formula for
  `K_q` (step 6). PB's glue is the shorter statement; CP's list is correct but longer than needed.
* **K4. Is (CYC) an input?** PB s.2.3 and P1: yes, in P(b) step 5 and through UR s.1 everywhere. Not listed in POLY s.6, the poly referee's
  C6, SYN4 R5, MECHANISM_STAGE4 s.2 or CP s.0.3; UR s.0 item 1 labels the cyclic form VERIFIED; `fq_conifold` records the isotropy
  `B(f_a, g) = 0` as "checked". **Verdict:** PB is right. (CYC) is now PROVED (CYC, critic, (u)) and VERIFIED twice.
* **K5. The relabelling `V = s Res_X(hatted integrand)`** (P(b) step 6). FSL s.10.4: SKETCH; RL s.0, s.3.4: "a genuine summability
  issue ... until refereed, this step of R5 stays open". PB s.2.7: PROVED. **Verdict (critic's reading of PB s.2.7):** the argument is
  complete. (i) In the original variables each `u^N`-coefficient of `V` is a finite sum (`Lam = O(u)`, `G = O(x^4)`, finitely many
  `z`-powers per `x^a`). (ii) Each hatted coefficient of `u^E w^beta` is a finite sum. (iii) The original term `(k, r, a, b)` becomes
  the hatted monomial `u^(r + (a + 4b)/5) w^(b - k) X^(a - k)`; only finitely many terms land on a given hatted monomial, and the
  monomial map `u^N w^beta -> u^(N - 4beta/5) z^beta` is injective. So the two expansions are one family summed in two orders, with
  finite fibres both ways (the regrouping `e^(-Lam/(xz)) = e^(-1/(Xw)) e^(-(Lam/u - 1)/(Xw))` holds in FSL's ring B). (iv) The key lemma (finitely many
  `pi`-terms per target monomial) puts the hatted residue sum in `Q((z))[[u]]`. Status: PROVED (u); a referee is recommended. Also
  VERIFIED: `pb_mech` (`u^20`, mod p, 185 keys), `pb_res` (exact, `u^40`), the poly referee's `rf_res` (`u^30`).
* **K6. Theorem LMI.** GEN s.3.3: CONJECTURAL, SKETCH (formal Legendre dual and a Zariski-density step). RA s.3.3: PROVED (two-variable
  chart and the delta lemma with a source (L2')). **Verdict:** RA's proof is correct as far as I checked: at `L != s*` the total Hessian
  determinant `-g''(L) a0(L) - alpha0^2` is non-zero and (Fub) applies; after `nu = mu - lambda/alpha0` the data are holomorphic in
  `Q = e^(-lambda/alpha0)` at `Q = 0` with determinant `-alpha0^2`; (L1c)(i) gives holomorphy in `L`; (L2') is a 2x2 Wick computation
  with `C_ss = 0`. PROVED (u). LMI is not needed for (Gamma), but it is the shortest proof.
* **K7. The free coefficients and `Ginf`.** SYN4 s.1.1: CONJECTURAL. GEN s.5.1 A: VERIFIED plus SKETCH. RA s.3.4: "PROVED given the
  inputs of Theorem D". LIT_3 s.4 item 2: GJR arXiv:1812.11908 Lemma 7.7, eq. (37), proves the formula for GJR's LG series; a transfer
  needs (L1) + (L2). **Verdict:** `c(z; 0) = Ginf(z)` is PROVED (u) directly (RA). `f^can_m = [z^(5m)] Ginf` needs only S3 and P(b) (not
  (QRR), not Theorem D); `f^QRR = f^can` needs `D = 1`. The LIT transfer route is superseded; GJR is prior art and should be cited.
* **K8. Regularity at `L = 1`.** The brief makes regularity part of (Gamma); PB P8: regularity is automatic, since by P(a)
  `c_n(L) = L^(4n) ph^SP_n(L^-5)` is a Laurent polynomial. **Verdict:** no conflict. CP proves regularity without P(a); given P(a) the
  content of (Gamma) is the value.
* **K9. The source of (QRR).** MECHANISM_STAGE4 s.2 and SYN4 s.2: arXiv:1702.06096 eq. (5) and arXiv:math/0108100 Thm 9.1 ("the Bernoulli
  form"). CP s.5 item 2, RF s.5: Thm 9.1 is untwisted, no twisted statement located. LIT s.1.2: GJR arXiv:1812.11908v1 eq. (17) and
  arXiv:1709.07392v1 eqs. (19)-(20). **Verdict:** LIT. Thm 9.1 (cotangent characters) and math/0008067 Thm 1 (tangent characters) are the
  untwisted analogue; 1702.06096 eq. (5) and 1803.01409 eq. (5) define the theory but state no `R(q = 0)` (their `R_k` carry the
  Zagier-Zinger normalisation, LIT s.1.1).
* **K10. "The 7-variable SP is not defined at `L = 1`"** (FSL s.8, warning 1) versus CP's conclusion at `L = 1`. **Verdict:** no conflict.
  CP uses the two orders only on `Lam` and reads the value off the chart, which is defined at `L = 1`.
* **K11. Which formal lemmas does R5 use?** The poly referee's C6 (four) plus Stirling's functional equation (SYN4); RL s.3.4 adds L3(iv),
  Lemma T(b), Lemma H and the relabelling; PB s.2.0 has FL1-FL7 in its own numbering. **Verdict:** the union, with this dictionary:
  PB FL1 = FSL FL1 (exact forms); PB FL2 (rank-one quotient, invariance) = FSL L1a and FL3-SP (s.5.4); PB FL3 (parameter derivatives)
  = FSL L1c(iii)/FL4; PB FL4-FL5 (residue) = FSL FL2-FL3 (s.10.2-10.3); PB FL6 = FSL L3(iv); PB FL7 (shift) = FSL Lemma T(b);
  plus Lemma H (RL s.3.2) and the relabelling (PB s.2.7). Numbering hazard: PB's FL2, FL3 are not FSL's FL2, FL3.
* **K12. (S2)** (P(b) step 4). Poly referee s.3: "VERIFIED through `w^50` and fixes the signs". PB s.2.6: PROVED (the model system,
  homogeneity, Prop. S). **Verdict:** PB (an upgrade). So SYN4's "R5 PROVED" rested on one VERIFIED-only step; it is now PROVED (u).
* **K13. Ranges.** (Gamma) exact `z^60` (GEN), `z^40` (RF), `z^18` (RA); mod p `z^80` (PB), `z^24` (RA), `z^500` (Stage 4). Methods differ;
  no conflict.
* **K14. The twist-weight sign** (`-5`): CP s.3.6, LIT s.1.4 (three ways), Stage-4 (H) at `n = 1`, GEN T2-T3 controls. No conflict.

## 6. Overclaims and underclaims (quote or point, then correction)

* **O1** GEN s.0, s.7 S1: (Gamma) "CONJECTURAL with a SKETCH proof". Superseded: PROVED (K1).
* **O2** GEN s.3.3, s.7 S4: LMI "CONJECTURAL", "SKETCH". Superseded: PROVED (u) (K6).
* **O3** GEN s.2 Lemma 2 "PROVED": the residue identity is proved; its identification with the twisted GW I-functions rests on formulas
  quoted without statement numbers (GEN says so in s.2 and s.7). Read "PROVED as a residue identity". Not used by the chain.
* **O4** CP s.3.2, last bullet: "`det H_ch` ... vanishes only at `L = 0, zeta^j`". False: it has double poles there (RF C1).
* **O5** CP s.3.10: "its Jacobian is removed by normalisation". The Jacobian `-L/X^2` belongs to the hatted amplitude (RF C2).
* **O6** CP s.5 item 7: polynomiality "PROVED elsewhere per SYN4 s.1.1". It needs P(a), s.3.10 and the degree bound (RF C3), or RA's
  `L = 0` chart (RA R4).
* **O7** CP s.2.8 and s.5 item 4: "a purely formal proof of (L3)" not available. Underclaim (K2).
* **O8** CP s.0.3, "this is the complete list" (and RF s.0, "the list of hypotheses in s.0.3 is complete"). It omits (CYC) (K4), (S2)
  (K12) and the relabelling (K5), and includes R2, which is not needed (K3). Correction: the list of part 1, s.1.
* **O9** CP s.0.3 and s.5 item 2: the twisted (QRR) "not a statement I have located". Superseded (K9). CP s.0.4 "No finite identity
  remains among the hypotheses" is true; read it with (QRR) (RF s.5).
* **O10** MECHANISM_STAGE4 s.2 and SYN4 s.1.2 (R5), s.2: "Theorem P(b) (PROVED given (U), GI Theorem G and four standard formal lemmas)".
  It also needs (CYC), (S2) (then VERIFIED only), the relabelling, L3(iv), Lemma T(b) and Lemma H. Correction: part 1, s.1.
* **O11** MECHANISM_STAGE4 s.2 (QRR): "arXiv:math/0108100 Thm 9.1 (the Bernoulli form)" as a source. Thm 9.1 is untwisted; cite GJR (K9).
* **O12** PB s.2.0: FL3 and FL6 "standard; only sketched here". Underclaim: both are PROVED in FSL (L1c(iii), FL4; L3(iv)).
* **O13** PB s.1, H9: "Not used anywhere: Lemma N, Lemma P, the degree bound ..., (Gamma), the chart argument". Garbled: it means that
  these items are not used by P(b) and the glue. Rephrase.
* **O14** RA s.0 item 5, s.3.4: `Ginf` "PROVED given the inputs of Theorem D": for `f^can` less is needed (K7).
* **O15** LIT s.1.4 "PROVED, given the statements of s.1.1-1.3": the GJR and Coates-Givental signs were read in the worker's own PDF
  text extraction, which drops minus signs (LIT_3 s.5); signs were reconstructed from ar5iv or the formula's structure. Read
  "PROVED given the quoted statements, whose signs need a human reading". The sign itself is independently pinned by `h_1 = ph_1`.
* **O16** FSL s.0 "Who uses what" and s.12: the four lemmas of C6 as the formal input of R5. Incomplete (RL C6; K11).
* **O17** LIT_2 s.2.3, (L3) "PROVED modulo writing it out": superseded by FSL s.9.3-s.9.4 (written, refereed).
* **O18** MECHANISM_STAGE4 s.1 and s.2: "(Gamma) (CONJECTURAL, complete formal outline, unrefereed)". Status update: PROVED (s.3).
* **O19** Compliance, not mathematics: pb overwrote its own first log `pb_mech_6.log` and deleted its own `pb_cyc_20.log` (reported, PB_AUDIT_2
  s.5); RA overwrote or deleted first-draft logs of its own (reported, RA_2 s.7); `pb/pb_core.py` (286 lines, reported as written in
  pieces) and `lemmas/referee/rl_core.py` (307 lines, not mentioned) exceed 250 lines.

## 7. What a human must check, in this order

1. **Theorem Gamma.** Either RA s.3 (one page: Theorem LMI with (L2')) or CP s.3 (three pages: steps (a)-(d), s.3.1-s.3.6), together
   with the lemmas they use: (L1a), (L1b), (L1c)(i), (L2) or (L2'), (L3), i.e. FSL s.5, s.6, s.4 (i)-(ii), s.7, s.9.3 (about six pages).
   Points to watch: the measures (`dy` for `Gamma(1 + 5s)`, `dx_j/x_j` for the others), the factor `beta` of (Fub), and the continuation
   (two rational functions agreeing on `Lam`, CP s.3.5).
2. **(CYC)** (`critic/CYC_PROOF.md`, two pages of elementary algebra; it is new and unrefereed).
3. **Theorem P(b)**: POLY s.6 with the poly referee's s.3 and PB s.2.5-s.2.9, above all PB s.2.6 ((S2), signs `eps_b`, weights,
   `kappa = C0`) and PB s.2.7 (relabelling, key lemma) (about four pages).
4. **The glue** PB s.3 with UR s.4 (one page): `Phi^QRR D_0 = Phi^can`, evaluated at `u = 1`.
5. **(QRR)**: GJR arXiv:1812.11908 eq. (17) in the typeset PDF (the sign of `-5 h_alpha` and `R = e^(-u/z) S C_alpha^-1`), and its sources
   arXiv:math/0008067 s.2.2-s.2.3 (Thm 1 and the Remark after it) and arXiv:math/0110142 Thm 1, Cor. 1; the vertex derivation of LIT
   s.1.4(a) (arXiv:1803.01409 eq. (26), Mumford's relation, arXiv:math/0108100 Thm 4.1); and the uniqueness step: does the full CohFT fix
   the Hodge twist `h` of Teleman's s.8.2 (Prop. 8.5(ii)), so that the `q = 0` R-matrix is the vertex R-matrix?
6. **The deepest project inputs**: Theorem U (UR s.3) and GI Theorem G (GI s.3).
7. **The Stage-2 list** (SYN2b s.6.2): LD, D(ii), D(iii), Theorem V with the sign dictionary, quintic Theorem R, Theorem S15, and the
   hypotheses of the quoted theorems (Q).

Optional, a code test rather than a proof check: the genus-3 degree-0 number `F_3(0) = 1/224000` from the project's R-matrix graph sum
at `u = 1` with `R = diag(h_(w_i))`. It would test the project's code convention at `z^5` (`b_3`), which the theorems of item 5 fix; the
genus-2 calibration (`F_2(0) = -19/9600`, FORMAL_QUINTIC_genus2 s.2) tests it at `z^1`, `z^3` only; LIT's `N_(3,0)` is a localisation number.

## 8. Next steps

1. **Referee the two new written pieces that the chain needs:** CYC_PROOF.md and PB s.2.6-s.2.7 (one referee each, about an hour each).
   Nothing else in the chain is unrefereed.
2. **Write-up** (`allgenus/note/`, SYN4b s.6.3, revised): (i) the polynomial structure; (ii) the glue (PB s.3) and, for the finite range,
   Lemma N; (iii) Theorem P (POLY s.6 with PB s.2); (iv) Theorem Gamma by the two-variable chart (RA s.3), with CP's seven-variable chart
   as a remark; (v) `Ginf` and the free coefficients (RA s.3.4; cite GJR Lemma 7.7); (vi) (CYC) as a lemma of the set-up; (vii) an
   appendix on formal stationary phase (FSL, with RL's Lemma H). Cite GJR Remark 7.2, which states the analogue of (Gamma3) for their LG
   integral without proof, and Givental arXiv:math/0108100 s.10 for the mechanism (untwisted, outline). Do not claim more novelty.
3. **Housekeeping** for the owners of the notes (the critic edited no existing file): RF C1-C6 and the remark on the `dy` versus `dy/y` blindness at `L = 1`; RL
   C1-C8; O1-O18 above; MECHANISM_STAGE4 s.2 (the input list) to be superseded by part 1, s.1.
4. **Optional computations:** the genus-3 code test of s.7; LMI for non-integer `M0` (GEN's `e_*` cases) through DLMF 5.11.8; GEN's
   twisted-`P^4` results (T2, all colours) as a corollary of LMI.
5. **Literature:** read Teleman s.8.2 verbatim; read GJR (17) and s.7 in the typeset PDF (the three glyph issues of LIT_3 s.5 concern s.7
   only); optionally a classical reference for special geometry (arXiv:hep-th/9204035), not needed since (CYC) is proved.

## 9. (CYC): what the critic proved (`critic/CYC_PROOF.md`; AI-written, unrefereed)

* *Statement.* In every conifold frame (`varpi_0` any power-series solution of PF4 with `varpi_0(0) != 0`, `t = g/varpi_0`, `g` the
  vanishing period), `PF4 = (1/varpi_0) D (1/Dt) D (u varpi_0^2 (Dt)^2) D (1/Dt) D (1/varpi_0)`; hence the cyclic system of UR s.1, with its
  closing relation `4 -> 0`, is equivalent to PF5 (`PF5 = z^5 D PF4 - 1`).
* *Proof, in five lemmas.* C1: PF4 is formally self-adjoint for `D = q d/dq` (`(q P(D))* = q P(-D - 1) = q P(D)`). C2: the Lagrange
  concomitant is a constant, alternating, non-degenerate form `B` on the solutions; its jet matrix is anti-triangular with anti-diagonal
  `+-u`. C3: hence `Sum Binv_ij y_i D^k y_j = 0` for `k <= 2` and `= 1/u` for `k = 3` (Griffiths transversality). C4: the automorphism
  `log u -> log u + 1` of `Q((u))[log u]` commutes with `D` and preserves `B`, so every power-series solution is B-orthogonal to `g`
  (formal Picard-Lefschetz). C5: in a symplectic basis adapted to the Lagrangian `span(varpi_0, g)`, C3 gives
  `th^2 (y/varpi_0) in Y span(1, t)` with `Y = 1/(u varpi_0^2 (Dt)^3)`, `th = (Dt)^-1 D`. Then the factorised operator kills all four
  solutions and has leading coefficient `u`; an operator of order `<= 3` killing a basis with invertible Wronskian matrix is zero.
* *Checks.* Part 1, s.4.1; each lemma separately, exact over Q to `u^70`, three frames, failing controls.

## 10. Critic record

**Files** (new; the two synthesis files at the place the brief names, everything else in `critic/`):

| file | content |
|---|---|
| `STAGE5_SYNTHESIS.md`, `STAGE5_SYNTHESIS_2.md` | this synthesis (163 and 175 lines) |
| `critic/CYC_PROOF.md` | the proof of (CYC) (s.9) |
| `critic/cr5_cyc.py` | Laurent series in `u` with absolute precision and powers of `log u`; `D`; PF4; the formal adjoint; the concomitant |
| `critic/cr5_cyc_run.py` | the checks of s.4.1 with controls |
| `critic/logs/cr5_cyc_run_{24,48,80}.log` | ALL CHECKS PASSED (0 s, 2 s, 8 s); relative paths only |

**Runs.** `cr5_cyc_run.py 16` (development) crashed in a control: for `t = f_2/varpi0` the series `Dt` is not a unit (`f_2 = O(u^2)`), so
my inverse refused it. I replaced that control by `t = f_1/varpi0` before any recorded run; no log was written by the crashed run. Then
`24`, `48`, `80`: all PASS. A left-over log line still printed `B(varpi0, f_2)`; I replaced it by check (d) and re-ran the three runs from
the final code (overwriting my own logs of minutes earlier; same verdicts). Re-runs of 20 scripts of pb, lit, general and the adversarial referee from scratch copies (part 1, s.4.2):
all outputs identical to the stored logs, timings stripped. No other failure.

**What I did not do.** I did not re-run the long computations (part 1, s.4.2), the Stage-4 codes, or the chart and lemmas workers' scripts
(their referees did). I did not read GJR, LP, Givental or Coates-Givental myself (I rely on LIT for them) and read Teleman only through a
summariser. I did not attempt the genus-3 code test of s.7.

**Compliance.**
* No existing file was modified, moved or deleted. The agents' scripts were run only from copies in the session scratch folder (they
  write their logs next to themselves).
* The brief asks for new files only inside `critic/` and also names `stage5/STAGE5_SYNTHESIS.md` and `_2.md`. I wrote those two new files
  at the named place (as Stage 4 did), and nothing else outside `critic/`.
* Python 3.14 stdlib only, `py -I -B`, working directory the session scratch folder; at most two of my processes at a time; every process
  I started has finished; none of anyone else's was touched. No wolframscript, no Mathematica, no cluster job.
* Web: the arXiv abstract page of 0712.0160 and three reads of its ar5iv rendering (fetch tool); no personal data in any request.
* Logs and notes contain no absolute path and no user name.
