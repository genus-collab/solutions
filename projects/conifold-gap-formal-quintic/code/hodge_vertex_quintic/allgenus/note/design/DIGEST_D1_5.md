# Digest D1, part 5: s.4 conventions, s.5 evidence, s.6 literature, s.7 open issues

Continuation of `DIGEST_D1_4.md` (piece map there). Same abbreviations.

## 4. Conventions, sign dictionaries, normalisations (as the sources fix them; differences flagged)

* **4.1 Variables.** D1.1. Flag: UR writes `th := z d/du`, EVID/POLY `th_u := u d/du`; ER uses `x = 1/s = L`; GR's flat derivation
  `D = ((L^6 - L)/5) d/dL`. Hatted variable: `zhat` (UR, EVID, GAM), `w` (POLY, PB, PS); `xi`/`Xi` (GAM) for POLY's `X`/`Y`.
* **4.2 Exponents and signs.** `c_b = (0, -1, 3, 2, 1)`, `eps_b = (+, -, -, +, +)`, `F_b(0) = eps_b`, `Y_b = eps_b I_(c_b)`, `Y^C_b = eps_b I^C_(c_b)`
  (UR s.1-s.2; GI s.2; PR s.3; PB s.2.6). The all-plus sign choice fails (PR s.1; PB2 s.4; PS2 s.5). In `kappa = C0` the signs cancel by
  `eps_b^2 = 1` (PS s.1).
* **4.3 Phases.** GAM: flat `e^(W/z)`; POLY: hatted `e^(-W/w)`; `W_GAM = -W_POLY`, equal amplitudes (Stage-4 K6; SYN4 s.3). FSL's phase `F`
  is `Phi` (GAM) or `-W`, `z = w` (POLY) (FSL s.1.4). Summands: GAM `(-q)^s Gamma(1+5s) prod_j Gamma(zeta^j/z - s)`, POLY
  `q^s Gamma(1+5s)/prod_j Gamma(1 + s - zeta^j/z)`: same normalised series (PB P9; CP s.1.7, reflection, branch integer `e = 0` forced by
  criticality). Critical equation in every form `-q (5 sigma)^5 = prod_j (zeta^j - sigma)`, i.e. `u sigma^5 = 1`; with `+q` it would be
  `sigma^5 (1 + 5^5 q) = 1` (CP s.1.7).
* **4.4 Measures.** `Gamma(1 + 5s) = Int y^(5s) e^-y dy` carries `dy`, `Gamma(w) = Int x^w e^-x dx/x` carries `dx/x` (CP s.1.7). (Gamma) is
  blind to the twist measure; the measure `dy` enters through P(b) and [HF] (RS2 item 40).
* **4.5 Normalisations of series.** `ph_0 = 1`; `c_0 = 1`; `Phi^SP` is divided by `sqrt(2 pi w/5) e^(-What(1)/w)`, the Gaussian factor times
  the amplitude value `(1-u)^(-1/2)` (CP s.1.7); `A_n = L R_n` (D1.9). GP s.5 (c) "absorbed by `D(0) = 1`" is wrong as worded: the factor
  `(2 pi z/5)^(1/2)` depends on `z` and cancels because it is common to both stationary phases (Stage-4 K10; GR C2).
* **4.6 Exponentials.** `D` is defined after stripping the same exponential on both sides; UR's `e^((u_i - u_i^c)/z)` and `e^(u_i/z)` differ
  by a constant that is not a power series in `z` (PS-C8). `u_i - u_i^c = -zeta^i (5/4) u^(4/5)/tau` (UR s.1). PB s.2.1 writes
  `e^(-What(1)/w) = e^(u_0/z)`; with D1.13 and UR s.1 one has `-What(1)/w = (u_0 - u_0^c)/z` (my arithmetic, not in the sources): the note
  should use one convention throughout.
* **4.7 Gauge normalisations.** Canonical `P_11(0, z) = 1`; stored `p_(5q,0)(0) = 0`, `P^s_11(0, z) = 1 - (43397/807975000) z^5 -
  (765368951400317/1305647201250000000) z^10 - ..` (GI s.1); `m5 = -43397/807975000`; canonical `p_(5,0)(0) = -p'_(4,1)(0) = 43397/807975000`
  (POLY s.3, = UR's `c_5`); the ODE residues `(0, 1/5, -3/5, -2/5, -1/5)` for `k = 0..4 mod 5` (UR s.3 Remark (d); ER s.1).
* **4.8 CFU, two normalisations.** PB/POLY: `Y(X)`, `What = V/Y + Y^4/4`, `V = tau^(-5/4)`, `Y_c = tau^(-1/4)`, `Y_1 = (1+T)/18`. GP s.4:
  `xi = Xi psi(Xi^5; u)`, `phi_k(1) = 0`, `What_GAM(Xi psi) = -(Xi^4/4 + Xi^-1)/tau`, `phi_1 = (1-y)/18`, `phi_2 = 13/567 - (101/4536) y - y^2/1512`.
  They agree under `Y = Y_c Xi` (my check at order `u`: `(1 + u/9)(1 - u(1-T)/18) = 1 + u(1+T)/18`; not in the sources). The Brieskorn
  reduction of GP s.5 (b) is in `Xi`; PB s.2.5's in `Y`. Use PB's versions (the chain is PB's).
* **4.9 Frames.** Programme frames `(alpha, beta)`: `varpi_0 = f_0 + alpha g + beta f_2`; used `(0,0)`, `(1/3, -2/7)`, `(7/11, -3/5)`. The Stage-5
  critic's (CYC) frames are `f_0 + a f_1 + b f_2`: its `(1/3, -2/7)`, `(7/11, -3/5)` are programme `(1/3, -109/210)`, `(7/11, -23/22)` (RS A2).
  Every frame named in s.5 is a programme frame.
* **4.10 QRR sign.** Twist weight `-5 zeta^i` (Euler sign); pinned by `h_1 = ph_1 = -3/20`: weight `+5`, no twist and the reversed tangent
  orientation all fail at `n = 1` (ER s.6); settled three ways (SYN5b K14). The `z^(5k)` parts are pinned by no such check (RF2 s.5): they
  are exactly what (N) measures and are fixed only by [QRR].
* **4.11 Free-coefficient indices.** `f_m` (POLY, the chain) = `f_(5m)` (EVID) = `x_m` (GAM); `kappa_m` (GAM) the canonical values. Values:
  `f_1 = 4537/3937500 = B_6(1/5)/6`, `f_2 = -177640111151/31007812500000`, `f_3 = -284328315069009259339/6226756347656250000000` (EVID s.4); the
  stored normalisation has `x_1 + m5 = 4437977/4039875000` (GAM s.1).
* **4.12 Numbering and label hazards.** PB's FL1-FL7 vs FSL (K11): PB FL1 = FSL FL1; PB FL2 = FSL L1a and s.5.4; PB FL3 = FSL L1c(iii)/FL4;
  PB FL4-FL5 = FSL FL2-FL3 (s.10.2-s.10.3); PB FL6 = FSL L3(iv); PB FL7 = FSL Lemma T(b); plus Lemma H and the relabelling. "(A2)": POLY's
  identity vs GAM's analytic step (Stage-4 K4). "(S2)": the model identity (PR, PB) vs GP s.1's change of variables; SYN5's "(S1)" (the
  lemma set) vs GP's (S1). "(H)": `D_0` in `Q[[z^5]]` vs UR s.4 purity vs FSL s.9.3 homogeneity vs RL's Lemma H. "(R1)": the relation of
  UR s.5 vs SYN4's step R1 (Lemma P). "(U)": the input vs FSL's Lemma U. POLY's Lemmas 1-4 vs GI Lemma 1. GAM's "Theorem A" = P(b) plus
  A(i); CP's "Theorem D" = `thm:glu-D`. Project tags (R0-R7, S1-S9, K, O, A, H1-H9, PS-C, E4) must not appear in the note (rule 5).
* **4.13 Symbol clashes** (for NOTATION.md): `c` (exponents `c_b`; cyclic coefficients `c_a`; the constant series `c(z)` of Lemma N and
  POLY s.3; `c(z; L)`, `c_n(L)`); `D` (`q d/dq`; `diag(D_i)`; `D_q`); `L` (`u^(-1/5)`; GI's field `Q((z))`; FSL's ring `L`); `Lam` (`-log(1-u)`;
  CP's `C \ {0, 1, zeta^j}`); `V` (`tau^(-5/4)`; the residue); `s` (`u^(1/5)`; the Mellin variable; `s(z)`; LMI's `s*`); `h` (QRR series; CFU
  `Y = X + uh`; DLMF shift); `tau` (`tau(u)`; CYC's automorphism); `w` (hatted variable; multisets `w_i`; GR's weight `w(f)`); `q`
  (parameter; index in `z^(5q)`); `kappa` (fq_rows `kappa_a`; `kappa(u) = C0`; GAM's `kappa_m`; `kappa = +-i` in (L2)); `eps` (`eps_b`;
  `L - 1`; L1c(iv)); `G` (phase part `G(x)`; Theorem G; `Ghat`; `G_b`; chart `G`); `P` (gauges; CP's power sums `P_k(L)`; ER's operator);
  `Q` (rationals; actual rows; `Q^m`, `Q^SP`; LMI's `Q`); `M` (`zD`; `M_0`; LMI's `M0`); `A` (cyclic matrix; amplitudes `A_n`; `A_e`; FSL's
  ring `A`; `A_St`); `B` (amplitude; CYC's concomitant; FSL's ring `B`); `F` (`F_b`; `F_inf`; FSL phase; ZZ's `F(w,x)`); `R` (R-matrix; `R^`;
  `R_n(L)`; PB's ring); `T` (`X^5`; `T_k`; transport `T_1`); `Y`, `y` (rows; CFU variable; `Y^C`; `y^(n)`; `y_j`; Morse `y`); `u_i` vs `u`.

## 5. Evidence for the appendix (finite checks; never cited by a proof; all logs end in PASS, every FAIL line a control)

| # | identity | range | field, frames, points | hypotheses | source |
|---|---|---|---|---|---|
| E1 | `D_0 = 1` ((H) + (N)) | `z^204` (`K_q = 0`, `q <= 40`; (H) `n <= 204`) | exact Q, frame (0,0); two independent codes | (U), via Lemma N | EVID s.4; ER s.5; MS4 s.4.1 |
| E2 | `D_0 = 1` | `z^1004` (`q <= 200`) | `2^89-1`, `2^61-1` | (U) | ER s.0, s.5 (worker `z^504`, EVID s.4) |
| E3 | gauge of the actual rows = independent canonical gauge, every slot | `z^100`, `u < 85`; connection route exact `z^40` | two primes; exact | gauge computed, regularity asserted | ER s.0, s.5 |
| E4 | `c^act = 1` (direct test) | `z^35`; `z^55` | two primes, two frames; `2^89-1` | (U) | EVID s.1 |
| E5 | `D = 1` by the first row, closing, Theorem G | `z^70` exact; `z^150` | frame (0,0); two primes, two frames | (U) | POLY s.4; PR s.1, s.6 |
| E6 | (H): `ph_n(1) = h_n`, `5 ∤ n` | `n <= 1004` / `300`; `n <= 500` | two primes / exact; `2^127-1`, `2^107-1` | none | ER s.5; GR s.3 |
| E7 | (Gamma) `c = R`, `c(z;1) = h` | `z^500`; `z^40`; `z^80` random `L`; `z^60` | `2^127-1` (`L = 2, 3/5, 4/3`), `2^107-1` (`L = 2, 9/7`); exact, six `L`; `F_p127`; exact `u = 1/2, -1/3, 3/7` | none | GR s.3; SYN4 s.3 E2 |
| E7' | (Gamma), further codes | `z^30`; `z^60`, `z^40`, `z^18`; `z^80`; `z^16` | exact, fourth code; exact (Stage 5); two primes; exact at `L = 2, -2/3, 1/2, 3` | none | MS4 s.4.2; SYN5 s.3; RS2 s.6.1; MS5 s.4.2 |
| E8 | `Phi^SP = Phi^can`, every coefficient | `w^80`, `u^75` (6156 per frame); `w^15` | `2^61-1` frames (0,0), (1/3,-2/7); `2^89-1` (0,0); exact | gauge by the U-construction | PB2 s.4 |
| E8' | `Phi^SP = Phi^QRR = Phi^can`; (A2) | `w^60`; `w^50`; `w^25` | exact at `u = 1/2, -1/3, 3/7`, `f^can` without QRR input; two primes | (U) | SYN4 s.3 D1, E2; PR s.1; POLY s.5 |
| E9 | `V = C0 (Y^C_Q)_0` | `u^40` (220 keys); `u^30`; `u^14` | exact, two frames; two primes, two frames; exact | none | PB2 s.4; PR s.1; POLY s.5 |
| E10 | mechanism of P(b): `eps_b pi_(c_b) = C0 G^can_b`, residue side | `u^20`, `w^26`, 185 keys (45 with `z^(>0)`) | `2^61-1` | (U) for `G^can` | PB2 s.4 |
| E10' | relabelling; mechanism; key lemma; (S2); gauge-free identity of `lem:thp-id` | 965 hatted coefficients; 4283; 884 keys; `t^80`; `w^25 u^25` | exact | none | PS2 s.5 (T1-T5) |
| E11 | (S2) | `t^60`; `t^80`; `w^50` | two primes; exact | none | PB2; PS2; PR s.1 |
| E12 | CFU, `V = tau^(-5/4)` | `u^20`; `k < 10`; GP variant `k <= 16` | | none | PB2; PR s.1; GR s.2.4 |
| E13 | `f^can = f^QRR`; `= [z^(5m)] Ginf` | `m <= 25` (no QRR input; then `D_0 = 1 mod z^126`); `m <= 100`; `m <= 30` | exact; two primes; exact | (U) for `f^can` | SYN4 s.3 D1, D2; GR s.3; PR s.1 |
| E14 | `deg ph_n <= floor(4n/5)`, equality `n != 2` | `n <= 1004` / `300`; `n <= 150`; `n <= 80` | two primes / exact | none | ER s.4-s.5; PR s.1; PB2 |
| E15 | Theorem U construction regular, resonances unobstructed | `z^80`; `k <= 100` / `40` | two primes; two primes / exact | none (evidence for (U)) | PB2 s.4; ER s.8 |
| E16 | three formulas for `[z^(5q)] c` agree | `q <= 12` | exact, frames (0,0), (1/3,-2/7) | (U) | SYN4 s.3 C1 |

**Failing controls** (each fails where predicted): stored normalisation differs first at `w^5 u^4` by `m5` (EVID s.5; ER s.6; PB2 s.4);
`f_1 -> f_1 + m5` gives `K_1 = m5`, regularity blind (EVID s.5); `tau -> 1` (`v = u`) breaks regularity, obstructs CFU at `k = 1`
(`r_1(1) = -5/9`) and leaves the gauge equation without kernel (EVID s.5; GR s.4; PB2 s.4); QRR weights without `-5` break (H) at `n = 1`
(EVID s.5); `B_96 -> B_96 + 1` leaves (H) and `K_q`, `q <= 18`, intact and gives `K_19 = 5^-95/9120` exactly (non-vacuity at `z^95`; uses
`lem:pol-p5k`) (EVID s.5; ER s.6); `B_8 -> B_8 + 1` breaks (H) first at `n = 7`; dropping the Stirling term of `Gamma(1+5s)`, or `b_k -> -b_k`,
fails at `w^1` (PB2; PR); `b_2 + 1` fails at `z^3` (GR s.4); dropping the `j = 0` Stirling term fails at `z^1` (GR s.4); residue without
Stirling terms, or `gamma = 1/2`, fails at `u^2 z^-1` (PB2; PR); all signs `+` in (S2) fails (PB2; PS2); PS2 s.5's controls (`Lam/u -> 1`,
no `X'(Y)`, no CFU, `V -> 1`, truncation below the key-lemma bound, `pi_-1` sign, `pi_2 <-> pi_3`, no `eps_b`, no `tau`-weights).
**Independence caveats** (state them in the appendix): further conifold frames are no independent evidence for (N) (ER s.3; Stage-4 K7);
regularity checks are automatic given (U) (EVID s.6); only odd `a_j` are independent (PR C2); the direct and reduced tests share the PF
operator, the QRR constants and the model (EVID s.6); finite checks are evidence, not proof (EVID s.6). Rule 5: no script names in the note.

## 6. Literature used in this scope

| ref | statement used | use | read in full text by the sources? | root `refs.bib` key |
|---|---|---|---|---|
| arXiv:0710.0889 (Zagier-Zinger) | Thm 4(ii) (`Phi_s` in `L Q[L]`, `n = 5`); Thm 3 (`Phi_1 = (3/20)(L - L^5)`); eq. (5), Lemma 1 (coefficientwise in `x`) | second proof of `prop:pol-deg`; (Gamma3) form | yes: GR s.2.1, s.2.8 ("looked up", "as quoted"); PR C1, s.5 | none (new key) |
| arXiv:1702.06096 (Lho-Pandharipande) | eq. (5) (definition of the formal quintic), pp. 5, 9 (`C_0 = I_0`, `C_1 = D(I_1/I_0)`, `L`) | set-up; states no `R(q = 0)` (K9) | yes, PDF pp. 1-5, 7-10 (ER s.9) | none |
| arXiv:1803.01409 (Lho-Pandharipande) | the dedicated paper; v2 eq. (5), eq. (26) via LIT | cite with 1702.06096 (Stage-4 K1); GAM's s.2.4.4/Prop. 5 unconfirmed (Stage-4 O11) | abstract only in Stage 4; LIT in Stage 5 (not read by me) | none |
| arXiv:math/0108100 (Givental) | s.10: Thm 10.1, Example 10.3, Lemma 10.4 (outline, untwisted, chart `x_0 = Q/(x_1..x_n)`); Thm 9.1 (untwisted `R`, NOT a source of the twisted (QRR): K9, O11) | prior art for the mechanism | yes: PR s.5; RF2 s.5 (Thm 9.1, p. 12) | MR1901075 |
| arXiv:1812.11908v1 (Guo-Janda-Ruan) | eq. (17) in the proof of Lemma 5.2 ((QRR)); Remark 7.2; Lemma 7.7, eq. (37) | (QRR) source; prior art for (Gamma3) and `Ginf` | by LIT only; signs need a human reading (O15) | none |
| arXiv:1709.07392v1; arXiv:math/0008067; arXiv:math/0110142; arXiv:0712.0160v3; arXiv:0903.1463v3 | eqs. (19)-(20); Thm 1 + Remark; Thm 1, Cor. 1; Theorem (2), Prop. 8.5(ii); Lemma 4.20 | (QRR) sources and uniqueness; prior art | via LIT / a summariser (SYN5 s.4.3) | none; none; MR2276766; MR2917177; MR2553377 (matched by title; verify its eprint) |

Chester-Friedman-Ursell is named for the CFU map (POLY s.5) with no reference in any source and no root key: name it only as Lemma CFU,
or add a reference by hand. The (QRR) and (Q) items belong to other digests; listed here only as the glue uses them.

## 7. Open issues for the writers

1. **Notation** (s.4.13): settle `c`, `D`, `V`, `s`, `h`, `L`, `tau`, `kappa` first; give `Q` (actual rows) a non-blackboard letter.
2. **(QRR) in the H9 form.** The glue uses only `A^QRR_n(q = 0) = h_n`; the passage from `R(q = 0) = diag(h_(w_i))` to H9 was not
   re-checked (PS s.4.3; RF2 C4): `\draftnote` at the statement of (QRR) (QRR digest). Never write "the u = 1 convention of ref_rows" (rule 5).
3. **Definition of the actual rows** in [R0] and of `A^QRR` as a WKB column of PF5 near `q = 0`: the glue's step (b) presupposes it; the
   set-up digest must state it with its source (UR s.4; GI s.6, unread by me). `cor:glu-Laurent` then removes the continuation issue (A4).
4. **P(b)'s list** must name (CYC), Lemma H, L3(iv), Lemma T(b) and the relabelling (O10, O8); FSL numbering only (K11).
5. **Relabelling and key lemma**: write PS s.2.2 with PS-C4's two sentences; replace "therefore a well-defined element of `Q((z))[[u]]`"
   by the A6 sentence; the weight lemma is an optional remark.
6. **`prop:pol-deg`**: ER s.4 step (5) is one terse sentence; expand from ER s.4 or `\draftnote`; say first that `psi_n` is a Laurent
   polynomial. Never claim the equality `deg ph_n = floor(4n/5)` (VERIFIED only; POLY's Lemma 3 states it).
7. **P(a)**: name the ring/domain where the flat and hatted SPs are compared (PB "on `0 < u < 1`"; CP3 s.3.10 "as rational functions" by
   L1c(i)); state that P(a) does not use [HF] while [HF] uses P(a) (no circularity). P(a) was re-checked by PR and RL, not by PS (author-facing).
8. **`rmk:thp-Ai`**: its proof ends with a finite exact computation (GR's 3x3 grid), forbidden as a proof step by rule 3; the grid values
   are not in the notes I read. Recommendation: omit, or state it as a remark without proof; it is not needed.
9. **Finite route**: present Lemma N and the connection form as "reading the normalisation at the conifold", flagged as not used by the
   main proof; (F_g) for `g <= 69` is evidence only ("VERIFIED given Theorem U", Stage-4 O1).
10. **Glue wording**: PS-C8's ring sentence; `D` after stripping the same exponential (s.4.6); `[u^0](Y^C_Q)_0 = 0` in step 7 (PS-C7).
11. **Unread dependencies** (other digests must cover): GI s.5 (`p_(0,0) = tau^(-1/4)/varpi_0`, i.e. `ph^can_0 = 1`); GI Prop. S
    (`D_0(z) D_0(-z) = 1`); UR s.4 covariance (`D_i(z) = D_0(zeta^-i z)`); (CYC) and its proof.
12. **Corrected wordings to carry**: "a statement at `q = 0` alone" -> a connection statement (PR C4); "the picture of UR s.8, now derived"
    -> a different phase with the same critical data (GR C6); "absorbed by `D(0) = 1`" -> a common factor (GR C2).
13. **Author-facing (DRAFT_NOTES.md, not the note)**: PB s.2.6-s.2.9, s.3 and `f^can_m = [z^(5m)] Ginf` are refereed once; no human check.
14. **Novelty**: cite math/0108100 s.10, GJR Remark 7.2 and Lemma 7.7, Iritani Lemma 4.20; claim no more (MS5 s.3.6; SYN5b s.8 item 2).
15. **Conflicts settled** (record only): Stage-4 K1 (cite both LP papers), K2 (degree bound PROVED), K3 (Theorem A(i) needs no U), K4,
    K5 (still a connection problem), K6, K7, K9, K10, K11 (ranges differ by method); Stage-5 K3 (both chains valid; glue shorter), K4
    ((CYC) is an input, now proved), K5 + A6 (relabelling complete; formal sum vs `Q((z))[[u]]`), K8 (regularity automatic given P(a)), K9
    (QRR source is GJR, not math/0108100 Thm 9.1), K11 (FL dictionary), K12 ((S2) PROVED); RS A1 ((CYC) scope), A4 (log-free rows derived).
16. **Placement** (aligning SYN2b s.6.3 s.8, SYN4b s.6.3, SYN5b s.8 item 2, MS5 s.5 item 2): a section "The normalisation problem and the
    polynomial structure" (`nrm`, `pol`; put `prop:pol-Ginf` after Theorem Gamma, it needs [LMI0], [HF], P(b)); a section "Theorem P"
    (`thp`); a section "The glue and `D_0 = 1`" (`glu`; `thm:glu-D` needs Theorem Gamma, so Gamma comes first or `thm:glu-D` closes the
    Gamma section); a subsection or appendix "The finite route" (`fin`) feeding the evidence appendix; s.4.2-s.4.10 into the sign-convention
    appendix; s.5 into the evidence appendix.
