# Stage 5 synthesis (critic): (Gamma), `D_0 = 1` to all orders, Conjecture N(g) for all g (2026-11-14)

Critic of the Stage-5 workflow `allgenus/stage5/`. Five workers: `chart` (two referees, `referee_formal` and `referee_adversary`),
`lemmas` (one referee), `general`, `pb`, `lit` (no referee). AI-written; no human has read any Stage-5 note. Labels: PROVED (argument
written, or quoted theorem with arXiv number and statement number), VERIFIED (exact computation; range, field or primes stated),
CONJECTURAL, SKETCH, REFUTED; (u) marks a PROVED item that no referee has checked. Notation of the brief and of
`../stage4/STAGE4_SYNTHESIS.md` (SYN4; part 2 SYN4b). Abbreviations: CP = `chart/CHART_PROOF.md` (+`_2`,`_3`,`_4`), RF =
`chart/referee_formal/REFEREE.md` (+`_2`), RA = `chart/referee_adversary/REFEREE.md` (+`_2`), FSL = `lemmas/FORMAL_SP_LEMMAS.md`
(+`_2`..`_5`), RL = `lemmas/referee/REFEREE.md` (+`_2`), GEN = `general/GENERAL.md` (+`_2`), PB = `pb/PB_AUDIT.md` (+`_2`), LIT =
`lit/LIT.md` (+`_2`,`_3`), CYC = `critic/CYC_PROOF.md` (this critic), UR, GI, EVID, POLY as in SYN4.
Part 2 (`STAGE5_SYNTHESIS_2.md`): contradictions, overclaims, human checks, next steps, the critic's record.

Read in full: the 22 Stage-5 notes and reports above; `../MECHANISM_STAGE4.md`; SYN4 and SYN4b; `../MECHANISM_STAGE2.md`;
STAGE2_SYNTHESIS s.2; EVID s.2.3-2.4; UR s.0-s.1; the docstrings of all 65 Stage-5 scripts and of `fq_conifold.py`. All 77 Stage-5 logs
were scanned: each contains ALL CHECKS PASSED (except `lemmas/referee/logs/rl_blind_9.log`, a recorded observation with no checks),
every FAIL-type line is a labelled control or falsifier, and none prints an absolute path.

## 0. Bottom line

* **(Gamma) is PROVED** (Theorem Gamma, CP s.0.1 and s.3; AI-written; refereed twice, both "sound"). Its only inputs are formal
  stationary-phase lemmas, each PROVED twice by independent agents (CP s.2; FSL s.1-s.9) and refereed (RF s.2; RL s.1). With FSL's
  purely formal proof of (L3) (s.9.3, refereed) the argument uses no analytic theorem at all; CP's own (L3) quotes DLMF 5.11.1 and 2.3.8.
  RA s.3 gives a shorter proof (two-variable chart, Theorem LMI): PROVED (u).
* **`D_0 = 1` in `Q[[z]]` is PROVED** given (U), (QRR) and the AI-refereed Stage-4 chain, through the exact identity
  `D_0(z) = c(z; 1)/h(z)` (PB s.3) and Theorem Gamma. Hence `D = 1` (all colours) and (F_g) for every `g >= 2`. Theorem P(b), the main
  Stage-4 input, was audited a second time from scratch (PB s.2): no mathematical gap; three implicit steps are now written out.
* **(CYC) is PROVED here** (CYC; critic; (u)). pb found that the cyclic form of the quintic system, used by P(b) step 5 and implicitly
  by the whole programme (UR s.0 item 1: "VERIFIED"), is missing from every hypothesis list. CYC proves it in a page of algebra
  (self-adjointness of PF4, the Lagrange concomitant, a formal monodromy `log u -> log u + 1`, special geometry), and the critic's
  new code VERIFIES each step exactly (s.4).
* **Conjecture N(g) holds for every `g >= 2` modulo exactly (U), (Q), (QRR)** (s.1). No finite identity is left among the hypotheses.
  Every other step is PROVED in the project, AI-written and AI-refereed, except four (u) items: CYC; PB s.2.6-s.2.7 (signs, weights,
  relabelling of the residue, the key summability lemma); the glue PB s.3 (an alternative to the refereed EVID s.2.3-2.4, so not
  needed); RA's Theorem LMI (an alternative proof, not needed). No human has checked any step.
* **(QRR) now has a literature statement for the formal quintic itself:** Guo-Janda-Ruan arXiv:1812.11908v1 eq. (17), in the proof of
  Lemma 5.2, gives `R(q=0) = C_alpha^-1 = h_(w_alpha)` to all orders, justified there only by citation of arXiv:math/0008067 and
  arXiv:math/0110142 (LIT s.1.2). The twist-weight sign is settled three ways (LIT s.1.4). The step "R of the CohFT at `q = 0` is the
  vertex R-matrix" rests on Teleman's classification arXiv:0712.0160; the critic located Theorem (2) of s.1.9 and Prop. 8.5(ii) of s.8.2
  through the ar5iv rendering (fetch-tool summary, not verbatim; s.4.3).
* **Nothing is REFUTED.** No order failed in any field, frame or route. By-products: `c(z; 0) = Ginf(z)` (RA s.3.4, PROVED (u);
  GEN VERIFIED to `z^20`), hence the free coefficients `f^can_m = [z^(5m)] Ginf` given P(b): the closed form that SYN4 s.1.1 listed as
  CONJECTURAL is PROVED (u); GJR arXiv:1812.11908 Lemma 7.7 is prior art for the same formula for their LG integral (LIT_2 s.3.6).

## 1. The assembled theorem after Stage 5

**Theorem (assembled; AI-written; AI-refereed except the (u) items; no human check).** For every `g >= 2`, Conjecture N(g) holds
modulo exactly:
* **(U)** Theorem U (UR s.3). PROVED in the project, AI-audited only. Evidence: PB's independent Theorem-U construction has every
  `p_(K,b)` regular and every resonance unobstructed through `z^80` at `2^61-1` and `2^89-1` (PB_AUDIT_2 s.4).
* **(Q)** The quoted theorems, unchanged since Stage 2: arXiv:1509.06954 Thm 8.1; arXiv:2212.00320 Thm 7.1; arXiv:math-ph/0702045
  Thm 5.1 with the linear loop equation of its s.4 (applied to the MODEL curve `x^4 + v/x`, as the user's ruling allows; topological
  recursion on the quintic's own mirror curve is used nowhere); arXiv:1303.1043 Def. 2.2, Prop. 2.4, Def. 2.5, Prop. 2.9; optionally
  arXiv:2502.18115 Cor. 1.2.
* **(QRR)** The `q = 0` normalisation of the actual R-matrix, `R(q = 0) = diag(h_(w_i))`, `w_i = {zeta^i - zeta^j}_(j != i) U {-5 zeta^i}`,
  used in the flat colour-0 form `Phi^QRR(z, u = 1) = h(z)` (EVID's "u = 1 convention of ref_rows", RF C4; PB's H9). Sources: GJR
  arXiv:1812.11908v1 eq. (17) (proof of Lemma 5.2), and GJR arXiv:1709.07392v1 eqs. (19)-(20) at `t = 0`; both by citation of
  arXiv:math/0008067 Thm 1 with its Remark and arXiv:math/0110142 Thm 1, Cor. 1. Vertex-level derivation (LIT s.1.4(a), PROVED given
  the quotes): arXiv:1803.01409v2 eq. (26), Mumford's relation `c(E) c(E^*) = 1` (arXiv:math/9810173 s.2.1, unnumbered) and
  arXiv:math/0108100v2 Thm 4.1. Uniqueness of R for the CohFT: arXiv:0712.0160v3 Theorem (2) (s.1.9) and Prop. 8.5(ii) (s.8.2), [H].

| input | used in | status after Stage 5 | what would close it |
|---|---|---|---|
| (U) | R0 (UR s.4); P(b) step 5; Stage-2 steps 2, 3, 6 | PROVED (UR s.3), AI-audited; VERIFIED in range (`z^80`, PB) | human check of UR s.3 |
| (Q) | Stage-2 steps 4-7 | quoted with statement numbers; hypotheses checked in Stages 1-2 | human check (SYN2b s.6.2) |
| (QRR) | the glue (PB s.3) or R3 (EVID s.2.4) | quoted (sources above); sign pinned (LIT s.1.4; `h_1 = ph_1 = -3/20`) | human reading of GJR (17) and its sources; Teleman s.8.2 verbatim |
| (CYC) | P(b) step 5; POLY Lemma 2 for `Phi^can`; UR s.1 hence R0 and Stage-2 steps 1-6 | PROVED (CYC, (u)); VERIFIED exactly (pb_cyc to `u^56`; cr5_cyc_run to `u^70`; 3 frames) | referee CYC |
| formal SP lemmas (S1), (PG), (L0), (L1a)-(L1c), (L2), (L3) | Theorem Gamma | PROVED twice (CP s.2; FSL s.1-s.9), refereed (RF s.2; RL s.1); VERIFIED on random data | human check of FSL s.5-s.7, s.9.3 |
| homogeneity (F3) = Lemma H; reflection identity | CP s.1.7, s.3.9, s.3.10; P(b) step 4 | PROVED (CP s.1.2; RL s.3.2-s.3.3) | - |
| residue FL2, FL3; relabelling `V = s Res_X` | P(b) steps 3, 6 | PROVED (FSL s.10.1-s.10.3, refereed; PB s.2.7, (u), re-checked by the critic, part 2 K5) | referee PB s.2.7 |
| Theorem Gamma (i)-(iv) | (Gamma) | PROVED (CP s.3; RF, RA) | human check of CP s.3 or RA s.3 |
| hatted/flat identity `c_n(L) = L^(4n) ph^SP_n(L^-5)` | (Gamma) => poly's (A2) | PROVED (CP s.3.10 with RF C2); VERIFIED (`rf_hat`, `n <= 12`) | - |
| Theorem P(b) `Phi^SP = Phi^can` | R5 | PROVED (POLY s.6 + poly referee s.3; re-derived PB s.2 (u)) given (U), (CYC), GI Thm G + Lemmas 1, V, Q, UR Prop. S, POLY Lemma 2, the formal lemmas | human check of PB s.2.5-s.2.9 |
| glue `D_0 = c(z;1)/h(z)` | `D_0 = 1` | PROVED (PB s.3 (u); refereed alternative: EVID s.2.3-2.4, R2 + R3) given R0, POLY Lemma 2, formal WKB uniqueness, (QRR) | - |
| GI Theorem G, Lemmas 1, V, Q | P(b) steps 5-7; Stage-2 steps 1-2 | PROVED (GI s.2-s.3), AI-refereed | human check of GI s.3 |
| Stage-2 chain, steps 1-8 | N(g) from (F_g) | PROVED, AI-refereed (STAGE2_SYNTHESIS s.2) | human checks of SYN2b s.6.2 |

## 2. The proof chain, step by step (with the exact sections)

* **S1. Formal stationary phase.** The functional: Wick contraction of the Taylor data at a non-degenerate critical point (CP s.1.2;
  FSL s.1.4 Def. 1.2; the same functional, RL s.3.1). Lemmas, each with a written proof: (S1) = FL1 integration by parts (CP s.2.1;
  FSL s.3.1); (PG) = (G5)-(G5') reweighting (CP s.2.2; FSL s.2); (L0) = (G4) products; (L1a) changes of variables (CP s.2.4; FSL s.5,
  Piola); (L1b) Fubini with `beta = (det H_vv(w)/det H_vv(w_0))^(-1/2)` (CP s.2.5; FSL s.6); (L1c) parameters (CP s.2.6; FSL s.4, incl.
  the device (iv) `eps -> z`); (L2) the delta lemma, EXACT (CP s.2.7 by a degree count, FSL s.7 by a change of variables); (L3) Stirling
  = SP of Euler's integral with `x^m`, every `m` in Z, both measures (CP s.2.8 mod DLMF 5.11.1, 2.3.8; FSL s.9.3 Proof A, formal);
  (F3) homogeneity = RL Lemma H (s.3.2). Refereed: RF s.2 (CP's), RL s.1-s.2 (FSL's); no error found.
* **S2. Theorem Gamma** (CP s.3; RF s.3-s.5; RA s.1-s.2). For `L` in `Lam = C \ {0, 1, zeta^j}`:
  * (a) `Omega_7 = z^-1 e^(Phi/z) d sigma dy prod dx_j/x_j` (CP s.3.1); the six Euler substitutions leave `z^-1` identically in `sigma`;
    critical point `p(L) = (L, 5L, zeta^j - L)`, `det H_7 = 1/(L^2 (1 - L^5)^2)` (CP s.3.2);
  * order (a): (L1b) with inner `(y, x_j)`, then (L0) and (L3), gives `Nhat_7 = c(z; L)` (CP s.3.3);
  * order (b): (L2) and (L1b) with inner `(sigma, x_0)` give `Nhat_7 = Nch(z; L)`, the chart series of
    `psi = 5 log y - y + Sum_(j>=1)(-x_j + (zeta^j - 1) log x_j) + q y^5/prod x_j` (CP s.3.4);
  * (c) `c_n` and `Nch_n` are rational in `L` and equal on `Lam`, and `1/det H_ch = -L^2 (1 + L + .. + L^4)^2` is a polynomial; so `c_n`
    is regular at `L = 1` with `c_n(1) = Nch_n(1)` (CP s.3.5);
  * (d) at `L = 1` (`q = 0`) the chart is the product of the Euler integrals of `Gamma(1 + 5/z)` and `Gamma((zeta^j - 1)/z)`; by (L3)
    `Nch(z; 1) = h_(w_0)(z)`, twist weight `-5` (CP s.3.6). **This is (Gamma).** Z/5: `c_n(zeta L) = zeta^-n c_n(L)` (CP s.3.9).
  * Shorter route (RA s.3, (u)): only the singular factor `Gamma((1 - sigma)/z)` becomes an Euler integral; the two-variable chart
    `g(sigma) + alpha0 (sigma - s*) mu - Q e^mu` is holomorphic at `Q = 0`, and the delta lemma with a source (L2') gives the value.
* **S3. Flat to hatted** (CP s.3.10 with RF C2; (L1a) for `sigma = L/X`, Lemma H for `z = L^-4 w`): `c_n(L) = L^(4n) ph^SP_n(L^-5)`.
  With S2: `ph^SP_n(1) = h_n` for all `n`, i.e. `Phi^SP(w; 1) = h(w)` (poly's "(A2)").
* **S4. Theorem P(b)** `Phi^SP = Phi^can` (POLY s.6; poly referee s.3; PB s.2, re-derived step by step):
  * P(a): `Phi^SP` solves HPF (FL4 with a `q`-independent amplitude, L3(iv), Lemma T(b) or FL7), so its coefficients are polynomials
    (POLY Lemma 2) (PB s.2.2; RL s.3.4);
  * step 1, the CFU map `Y = X psi(X^5; u)`, unique and polynomial at each `u`-order (PB s.2.4); step 2, reduction modulo exact forms to
    `Y^-1 .. Y^3`, Z/5 grading (PB s.2.5); step 4, signs `eps_b`, weights, `kappa = C0`, with (S2) PROVED (PB s.2.6, (u));
  * step 5: both functions solve HPF (for `Phi^can` this is (CYC)) and have `ph_0 = 1`, so `Phi^SP = Phi^can s(z)` and `P^SP = s P_can`
    is regular and graded; Theorem G gives `(P^SP Y^C)_0 = s(z) (Y^C_Q)_0` (PB s.2.8);
  * step 6: the residue `V = s Res_X(..)` (relabelling and key lemma, PB s.2.7 (u)) equals `C0 s(z) (Y^C_Q)_0`; step 7: at `u^1` both
    sides are `-1/z`, so `s = 1` (PB s.2.9).
* **S5. The glue** (PB s.3 (u); refereed alternative EVID s.2.3-2.4): UR s.4 (`Q D = P_can Q^m(v(u))`, given (U)), polynomiality of
  `ph^can_n` (POLY Lemma 2, via (CYC)) and formal WKB uniqueness give `Phi^QRR D_0 = Phi^can`; at `u = 1`, (QRR) gives
  `h(z) D_0(z) = Phi^can(z, 1)`. With S4 and S3: `D_0(z) = c(z; 1)/h(z)`, EXACTLY. Lemma N's formula for `K_q` is not used.
* **S6. `D_0 = 1`** by S2. All colours: `D_i(z) = D_0(zeta^-i z)` (UR s.4), so `D = 1` (equivalently CP s.3.9 with `w_i = zeta^i w_0`).
* **S7. (F_g) for every `g`**: `D_0 = 1 mod z^(3g-2)` (S6).
* **S8. Conjecture N(g)** from (F_g), (U), (Q) by the Stage-2 chain (STAGE2_SYNTHESIS s.2, steps 1-8: vanishing solutions, (G),
  symplecticity, transport, quintic Theorem R, Theorem S15, the Model Lemma on K, conclusion).
* **S9. (CYC)** (CYC s.0-s.7) feeds S4 (step 5), S5 (polynomiality of `Phi^can`) and S8: by UR s.1 the quintic rows of `fq_rows` then
  satisfy the closing relation `4 -> 0` of the cyclic system, so Theorem U applies to the actual quintic.

## 3. Status of the three questions

* **(Gamma): PROVED**, modulo no external input. Proof: CP s.3 (or RA s.3) from the formal lemmas of S1; the DLMF statements of
  CP s.2.8 are avoidable (FSL s.9.3). AI-written; refereed twice (RF: sound, six minor corrections; RA: sound); not human-checked.
  Evidence, independent of the proof: exact over Q through `z^60` (GEN T1, Laurent expansion at `L = 1`), `z^40` (RF `rf_sigma`),
  `z^30` (Stage-4 `gamma_interp`), `z^18` (RA); mod two primes through `z^80` (PB) and `z^500` (Stage-4 gamma referee).
* **`D_0 = 1` to all orders: PROVED** modulo (U) and (QRR), using the Stage-4 chain R0 + R5 (AI-refereed; re-audited by PB), (CYC) (u),
  PB s.2.6-s.2.7 (u: (S2), the relabelling, the key lemma) and the formal lemmas. Exact evidence: `D_0 = 1 mod z^205` (Stage 4, given (U)).
* **Conjecture N(g), every `g >= 2`: PROVED** modulo (U), (Q), (QRR) (s.1): every step written; AI-refereed except the (u) items of
  s.0; none human-checked. Genus 2-5 remain VERIFIED directly by graph sums.
* **Free coefficients** (by-product): RA s.3.4 proves `c(z; 0) = Ginf(z)` (LMI at `s* = 0`, `alpha0 = 5`, `M0 = 1`; DLMF 5.11.8 as a
  formal identity). With S3 and S4: `f^can_m = [u^(4m)] ph^can_(5m) = [z^(5m)] Ginf`, given P(b) alone (no (QRR)); `f^QRR = f^can` needs
  `D = 1`. Also `c_n in Q[L]` without the degree bound (RA R4). PROVED (u); VERIFIED (GEN `z^20`, RA `z^18`, Stage 4 `m <= 100` mod p).

## 4. The critic's own checks (new code or scratch copies; `py -I -B`, working directory in the session scratch folder)

### 4.1 (CYC): `critic/cr5_cyc.py` + `critic/cr5_cyc_run.py` (exact over Q; logs `critic/logs/cr5_cyc_run_{24,48,80}.log`)

* Frobenius basis at `u = 0` from the operator alone; `gamma = 7/10` FOUND from the `u^-1` obstruction; `f_0, f_1, f_2` agree with UR.
* (a) `PF4* = PF4` exactly; (b) the concomitant `B` is log-free and constant on all pairs of solutions (to `u^72`), alternating,
  `det = 4`, `B(f_i, g) = 0` for `i = 0, 1, 2`; (c) `T_0 = T_1 = T_2 = 0`, `T_3 = 1/u` (to `u^71`); (d) `B(varpi0, g) = 0` in a frame, while
  the two control pairs have `B = -7/5, -2/5`; (e) `M y = 0` for all four solutions and
  `M u^n = PF4 u^n` (`n <= 4`) in three frames (to `u^70`). Controls fail as they must: `PF4 + uD` is not self-adjoint; the top-order part
  of `B` alone is not constant; the non-isotropic pairs `t = f_1/varpi0` and `t = (f_1 + f_2/2)/varpi0` violate (CYC). 8 s at `N = 80`.

### 4.2 Re-runs of the unrefereed agents and of the adversarial referee (scratch copies; logs compared with timings stripped)

| agent | scripts and arguments (as stored) | result |
|---|---|---|
| pb | `pb_cyc 60`, `pb_glue`, `pb_res 10` and `40`, `pb_s2 60`, `pb_exact 15`, `pb_mech 6`, `12`, `20`, `pb_run_phi 10 14` and `40 40` | 11 logs identical |
| lit | `lit_vertex_check` (`N_(2,0) = -19/9600`, `N_(3,0) = 1/224000`; the weight `+5` control) | output identical |
| general | `g5_valid 6 24`, `g5_t6c 20`, `g5_t1 40 40 Q 12` ((Gamma) exact to `z^40`), `g5_t3 20` (12 theories, 40 colours) | 4 logs identical |
| chart/referee_adversary | `ra_l2src 8 8`, `ra_poly 8 6`, `ra_run 8 10`, `ra_2var 6` | 4 logs identical |

Not re-run: the long runs (`ch_chart 3 4 6 5 7`, 2137 s; `g5_t1 60`; `g5_t2 24`; `pb_run_phi 80 72`; `rf_sigma 40 10`; `ra_run 18 24`).
The chart and lemmas workers' logs were re-run by their referees (RF s.7, RL s.6): all identical.

### 4.3 Literature (web; public pages only, no personal data in any request)

* arXiv:0712.0160 (Teleman, "The structure of 2D semi-simple field theories"; v3 of 2012 agrees with the published version): the
  abstract page, and the ar5iv rendering through the fetch tool's summariser (NOT verbatim). Located: Theorem (1) in s.1.7 (homogeneous
  theories are reconstructible from genus-zero data); Theorem (2) in s.1.9 (the CohFTs with a given semi-simple Frobenius algebra form
  the orbit of the trivial theory); s.8.2 "ambiguity for inhomogeneous theories", Prop. 8.5(ii) (over a given Frobenius manifold the
  theories differ by a Hodge-type twist `h`, the same `E`); Prop. 8.3, Cor. 8.4 (homogeneous case). The step "the R-matrix is fixed by
  the full CohFT, so `R(q = 0)` is the vertex R-matrix" is not literally stated in what I saw: a human check (part 2, s.7 item 5).
