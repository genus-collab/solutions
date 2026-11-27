# Stage 5 referees: synthesis, part 2 — every correction; final status (2026-11-14)

Continuation of `REFEREES_SYNTHESIS.md` (s.0 bottom line, s.1 the reports, s.2 adjudications A1-A10, s.3 the input table, s.4 the
"(u)" items). Same labels and abbreviations. Part 3 (`REFEREES_SYNTHESIS_3.md`): the critic's checks, human checks, next steps,
record. AI-written; no human has read it. No existing file was edited: each correction is addressed to the owner of the file.

## 5. Every correction, by the file it concerns (duplicates merged; none changes a statement, a hypothesis or a chain label)

**`critic/CYC_PROOF.md`** (Stage-5 critic)
1. s.1 (CF-E1 with CA-K3): "holds exactly for `g = f_1 + (7/10) f_2` (computed)" claims more than the code shows (it searches only
   `span(f_1, f_2)`). Write the 3x3 obstruction matrix on `span(f_0, f_1, f_2)`, rows `u^-3, u^-2, u^-1`:
   `[[-2,0,0],[17/5,0,0],[-7/5,-7/5,2]]`, kernel `Q(0, 1, 7/10)`; equivalently `gamma = -P_-2'(1)/P_-3'(2) = 7/10`, the unique `gamma`
   with `B(f_0, f_1 + gamma f_2) = -7/5 + 2 gamma = 0`. (Hand: CF s.2.3, CA s.3; VERIFIED `rc_cyc` B3, C4.)
2. s.3 (iv) (CF-E2): replace "computed exactly" by `Pf(Bmat) = -B(f_0, f_2) B(f_1, y_3) = 2` (`B(f_1, f_2) = B(f_2, y_3) = 0`, all
   `u^0` coefficients); optionally `det Wr = -2/u^2` (Abel's identity and Lemma 0). (VERIFIED `rc_cyc` C2, C3, C6.)
3. s.8 rows (d), (e) (CF-E3 = CA-K1): relabel the frames (s.2 A2).
4. s.0 (CF-E4): the file is `fq_conifold.py` at the package root; `natural_form/fq_conifold.py` does not exist (confirmed here;
   copies only under `natural_form/genus5_package_review/`). The same slip is in the Stage-5 brief.
5. s.0 statement (CF-E5 with CA-K2): the Theorem holds for every non-zero power-series solution `y_0` with `t = g/y_0` non-constant,
   and for every B-Lagrangian pair of `S cap K`; `y_0(0) != 0` is not used, it is the condition for the Theorem-U unit shape. Read
   "in every frame" (last line of s.0) as "in every conifold frame"; the large-radius frame and the ring N are not covered (not needed).
6. Notation (CF-E6): `M` (operator and `zD`), `tau` (automorphism, `tau_k`, the project's `tau(u)`), `f_i` (basis and `y_i/y_0`), `P`
   (`prod(D + k/5)` and the jet matrix), `y_3` next to `y_3'`. CF suggests `M_cyc`, `sigma`, `phi_i`, `Pi_ab`.
7. s.6 (CF-E7): add that `tau_1 = tau_2 = 0`, `tau_3 = 1/(u y_0^2 Dt^3)` hold for every unit `y_0` and `t` with `Dt` a unit (C3 is
   basis-free); isotropy enters only through the adapted basis, which exists iff `B(y_0, y_1) = 0`.
8. Corollary (CF-E8 with CA-K4): name the rings. CF: the identities lie in `Q[[u]][D]` (frames with `y_0(0) != 0`), so they hold on
   every differential `Q[[u]]`-module in which `z` is constant (the `s`-series rows; the rows `X_a` of `fq_rows` after conjugation by
   `e^(u_i/z)`). CA: the `Y_a` may live in any differential extension of `(K, D)` in which the `c_a` are invertible; the programme uses
   `Q(zeta)((s))((z))` for the WKB rows (`zD + zeta^i L` is the image of `zD` under `D -> D + zeta^i L/z`), `L((u))` with `L = Q((z))`
   for `P_can Q^m` and `Y^C_Q`, and the hatted series via `z = s^4 w`. Add CA's remark: the existence of log-free rows at all orders
   is a separate fact, supplied in the chain by PS s.4.2 (s.2 A4).
9. s.7 (CA-K6): `Y = 1/(u y_0^2 Dt^3)` is normalised by the Darboux form of B; the GW Yukawa coupling of FORMAL_QUINTIC_genus2 s.1 is
   `5 Y` (`eta = 5 delta`); the row coefficient `c_1 = Y Dt` has no 5 (a 5 there breaks the Corollary: `rc_cyc` E4); the `u` is `A_4`.
10. s.9 (CA-K5; status): list C4's further use (it proves that `Cf` respects `fq_ring`'s `i_3` relation). Status: "refereed twice
    (`referees/cyc_formal`, `referees/cyc_adversary`)"; "VERIFIED twice" becomes five times (s.3 row (CYC)).

**`STAGE5_SYNTHESIS.md`, `_2.md`** (Stage-5 critic)
11. s.1 table row (CYC) and s.4.1 (d), (e): the "3 frames" of `cr5_cyc_run` are programme `(0,0)`, `(1/3,-109/210)`, `(7/11,-23/22)`.
12. s.1 table row (CYC), column "used in": add the glue, PB's `Q^SP`, EVID Lemma N (PS-C9) and `Cf` on N (CA-K5).
13. s.2 S9 (CA-K7): "so Theorem U applies to the actual quintic" should read: Theorem U needs no (CYC); (CYC) supplies (i) that the
    actual rows solve the system to which Theorem U is applied (UR s.4) and (ii) that every solution of that system gives a PF5
    solution (P(b) step 5).
14. SYN2 s.9 statement: the scope of item 5 (CA-K2).
15. SYN2 K5: point (iii) needs PS-C4's two sentences; point (iv) needs the weight lemma, or the remark that the identity is
    coefficientwise with `V in Q[z^-1][[u]]` (s.2 A6).
16. Every "(u)": as in part 1, s.4. SYN s.0 bullet 4 ("except four (u) items") and SYN2 s.8 item 1 are now discharged (AI).

**`pb/PB_AUDIT.md`** (pb)
17. s.2.0 (PS-C1): FL3 is used in (S2) as well (the `v`-derivative at the moving critical point): cite FSL FL4 (s.4.5).
18. s.2.7, FL4 (PS-C2): the ring R needs `ord_(1/w) c_n -> +oo` as `|n| -> oo`; or use FSL's ring L (s.10.1, refereed).
19. s.2.6 (PS-C3): cite L1a and Lemma H (RL s.3.2) for `x = sigma x'`, `z -> sigma^4 z'` and for `Y = Y_c eta`, `w -> tau w`.
20. s.2.7 (PS-C4): name the family (`k`, `r`, multisets of atoms of `log(e^(-G/z) at)`) and add: the hatted fibre is finite (the only
    hatted-order-0 atom is `X^4/(4w)`; the X- and w-balances fix `k` and its multiplicity); `H` computed in `L[[u]]` is the family's
    hatted-order sum (the regrouping holds coefficient by coefficient).
21. s.2.7 key lemma (PS-C5): "therefore a well-defined element of `Q((z))[[u]]`" needs the weight lemma `5i - m - 4j >= 0` (PS s.2.3).
22. s.2.8 (PS-C6): the four-item re-check of Theorem G for `P^SP` is redundant (`P^SP = s P_can` is a regular gauge).
23. s.2.9 (PS-C7): state `[u^0](Y^C_Q)_0 = 0` (`y^(0) = -e_1`), used in the comparison at `u^1`.
24. s.3 (PS-C8): state the ring of comparison `Q[L, 1/L][[z]] c Q((s))[[w]]`; `D` is defined after stripping the same exponential on
    both sides (UR's `e^((u_i - u_i^c)/z)` differs from `e^(u_i/z)` by a constant that is not a power series in `z`).
25. s.3 "Same object" (PS-C9): read Lemma N as "`K_q = [z^(5q)] c` where `Ghat(Phi)` is the hatted form of `P_can c`, and
    `P_11(0, z) = c(z)` for that gauge"; add (CYC) to Lemma N's inputs.
26. s.1 H8 (PS-C10): formal WKB uniqueness is PROVED in one line (simple root; transport operator first order with kernel `Q L`).
27. s.2.8, s.2.9 (CA-K2): "`t_c'(0) = varpi_0(0) = 1` in every frame" holds in every frame `(alpha, beta)` of fq_conifold.

**Other Stage-4 notes** (their owners)
28. POLY s.3 (`b_0(0) = t_c'(0) = 1`), GI Lemma Q(d), P(b) step 7 (CA-K2): state that these normalisations refer to the fq_conifold
    family `varpi_0 = f_0 + alpha g + beta f_2`, `t_c = g/varpi_0`.

**`chart/referee_adversary/REFEREE.md` s.3 and `REFEREE_2.md` s.5** (RA)
29. s.3.3 step 2 (LF-C1; LA s.1.3): write out the identification of the critical points (for `0 < |Q| < r`, `L := sigma_c(Q)`,
    `lambda := -alpha0 log Q`, `nu_c := mu_c + log Q`); no uniqueness statement is needed.
30. s.3.3 step 2 (LF-C2): (F2) removes `lambda s*`; `Q^M0` drops out because `Nhat = N/N_0`.
31. s.3.1 (LF-C3 = LA-C3): define `g_k`; "`L != s*` near `s*`" means the punctured disc where `g`, `F_m` are holomorphic and
    `alpha0 + (L - s*) g''(L) != 0`.
32. s.3.4, (Gamma) bullet (LF-C4 = LA-C4): `F` is `e^(S_1)` without its `j = 0` factor `a0^(-1/2)`, i.e. it contains `(5 sigma)^(1/2)`.
33. s.3.4, orbifold bullet (LF-C5): state the normalisation of the value (division by `(2 pi)^(5/2) z^(3/2) e^(g(0)/z) e^(-g'(0)/5) F(0,0)`).
34. s.3.3 step 1 (LA-C1): introduce `Omega_2` as a formal integrand (the Euler integral's `z^(a0/z)` would shift `lambda`).
35. s.3.1 input list (LA-C2): add (F1)-(F2), the translation case of (L1a), CP s.1.4, the `C[[z]]`-linearity of N.
36. s.3.4 `Ginf` (LA-C5): DLMF 5.11.8 can be replaced by the formal Lemma F5.11.8 (s.2 A7).
37. s.0 item 4, s.3.4 "What this says" (LA-C6): "ANY holomorphic phase and amplitude" needs: holomorphic at `s*` and free of `lambda`
    (holomorphic dependence on `Q` allowed), `a0` exactly linear, `F(s*, 0) != 0` (a genuine hypothesis: LA's X9 has poles).
38. s.3.4 (LF clarification): "without Theorem P(a)" applies to `c_n in Q[L]` only; no negative powers of `u` in `ph^SP_n` needs P(a).
39. RA_2 s.5 (LF-C6; LA s.2.3): non-integer `M0` is now PROVED by LA (s.2 A3), no longer "same proof".
40. Write-up note (LA): (Gamma) is blind to the twist measure; the measure `dy` enters through P(b) and the hatted/flat identity.
41. Status (LA-C7): LMI and (L2') "PROVED, unrefereed" becomes "refereed twice (`referees/lmi_formal`, `referees/lmi_adversary`)".

**`general/GENERAL.md` s.3.3** (general)
42. LA s.5: read `Lam(0)` as the constant term of `Lam(z)`; the label CONJECTURAL is obsolete (as SYN2 O2 said).

**Outside Stage 5** (owners of the files; NOT edited)
43. `fq_conifold.py` docstring (CF aside): "Lam-free elements give Lam-independent images" should read "Lam-independent when written in
    the `D^j I_1`" (`i_0` is a counterexample); the rows are of the good kind.
44. FORMAL_QUINTIC_genus2 s.1: "Cf preserves the Wronskian relation" and "`B(g, analytic) = 0` automatically" are now PROVED (C4).
45. UR s.0 item 1, "Cyclic form (VERIFIED)": PROVED by CYC (as SYN K4 said).

**The referees' own reports** (my corrections; none affects a verdict)
46. CF REFEREE_2 s.7: "`y_0` and `e` are units of `Q[[u]]`": units of K suffice (A1).
47. LF REFEREE s.0 table, row "non-integer `M0`": SKETCH superseded by LA (A3).
48. LA REFEREE s.0 item 5 and CA REFEREE_2 s.5 Remark: their new lemmas are unrefereed by any other agent; status in s.6.3.

## 6. Final status after the referees

### 6.1 (Gamma): PROVED
*Statement.* The colour-0 full-Stirling saddle series `c(z; L)` of the Mellin-Barnes summand (`c_0 = 1`) has every `c_n` regular at
`L = 1` (`q = 0`), and `c(z; 1) = h_(w_0)(z)`; by Z/5, `c(z; zeta^i) = h(zeta^-i z)`.
* *Two independent proofs, each AI-refereed twice:* CP s.3 (seven-variable chart; RF, RA) and RA s.3 (Theorem LMI at `s* = 1`,
  `alpha0 = -1`, `M0 = 0`, with (L2'); LF, LA). AI-written; no human check.
* *Inputs:* only the formal stationary-phase lemmas of S1, PROVED twice and refereed ((S1), (PG), (L0), (L1a)-(L1c), (L2), (L3), with
  (L3) formal in FSL s.9.3 or modulo DLMF 5.11.1, 2.3.8 in CP s.2.8; (F1)-(F3); LMI adds (L2'), PROVED in RA s.3.2 and refereed by
  LF, LA), the holomorphic implicit function theorem (CP s.1.4), and CP s.1.6 (cancellation of the `sigma`-dependent powers of `z`,
  which identifies `c` with the summand's saddle series).
* *Evidence (independent of the proofs):* exact over Q at `L = 1` through `z^60` (GEN T1; LF Q3), `z^40` (RF), global Laurent fit
  `n <= 30` (LF); mod primes `z^80` (PB), `z^500` (Stage-4 gamma referee).

### 6.2 `D_0 = 1` to all orders, and Conjecture N(g) for every `g >= 2`
* **`D_0 = 1` in `Q[[z]]`: PROVED modulo (U) and (QRR).** Chain S1-S6 and S9 of SYN s.2, every step AI-refereed: (Gamma) (6.1); the
  hatted/flat identity (CP s.3.10, RF); Theorem P(b) (POLY s.6 with the poly referee; PB s.2 with PS); the glue `D_0 = c(z;1)/h(z)`,
  EXACT (PB s.3 with PS; alternative EVID s.2.3-2.4 with its referee); (CYC) (CYC with CF, CA); GI Theorem G and Lemmas 1, V, Q; UR
  Prop. S; POLY Lemma 2; formal WKB uniqueness (PROVED, PS-C10). (U) enters through R0 (UR s.4) and P(b) step 5; (QRR) enters only as
  H9, `A^QRR_n(q = 0) = h_n`. All colours: `D_i(z) = D_0(zeta^-i z)`, so `D = 1`. Exact evidence: `D_0 = 1 mod z^205` (Stage 4, given (U)).
* **Conjecture N(g), every `g >= 2`: PROVED modulo exactly (U), (Q), (QRR)**: (F_g) for every `g` from `D = 1` (S7), then the Stage-2
  chain (S8; STAGE2_SYNTHESIS s.2, AI-refereed). Genus 2-5 remain VERIFIED directly by graph sums.
* **The inputs, exactly** (unchanged from SYN s.1; this round adds none and removes none):
  * (U) Theorem U (UR s.3): PROVED in the project, AI-audited only;
  * (Q) arXiv:1509.06954 Thm 8.1; arXiv:2212.00320 Thm 7.1; arXiv:math-ph/0702045 Thm 5.1 with the linear loop equation of its s.4
    (applied to the model curve `x^4 + v/x` only; no topological recursion on the quintic's mirror curve); arXiv:1303.1043 Def. 2.2,
    Prop. 2.4, Def. 2.5, Prop. 2.9; optionally arXiv:2502.18115 Cor. 1.2;
  * (QRR) `R(q = 0) = diag(h_(w_i))`, used as H9; sources as in SYN s.1 (GJR arXiv:1812.11908v1 eq. (17); arXiv:1709.07392v1
    eqs. (19)-(20); Teleman arXiv:0712.0160v3 Theorem (2), Prop. 8.5(ii) for uniqueness, read only through a summariser).
  * No finite identity remains among the hypotheses: the two finite computations inside (CYC) (`gamma = 7/10`, `Pf(Bmat) = 2`) are
    now written out by hand (CF s.2.3, s.4.4) and re-derived here (`rc_cyc` B3, C2-C3).

### 6.3 What is still unrefereed, or refereed by one agent only
* **Steps of the chain with no AI referee: none.**
* **Refereed by one agent:** PB s.2.6-s.2.9 and s.3 (PS); `c_n in Q[L]` and `f^can_m = [z^(5m)] Ginf` (LF s.4.3). (Their Stage-4
  counterparts, POLY s.6 and EVID s.2.3-2.4, have their own referees.)
* **New statements of this round, checked by no second agent except the critic** (none is used by the chain):
  * PS s.2.3, the weight lemma: PROVED (PS); checked by hand here (s.2 A6).
  * Formal DLMF 5.11.8: PROVED twice independently (LF s.4.2, LA Lemma F5.11.8, the same argument); checked by hand and VERIFIED here
    through `w^-60` (`rc_l3c` D1-D3).
  * LA Lemma L3-C and LMI for every `M0` in C: PROVED (LA); checked by hand and VERIFIED here (symbolic `M`, `z^12`; `rc_l3c` L1-L3).
  * CA-K2, (CYC) for every Lagrangian pair of `S cap K`: the proof is CYC's verbatim; VERIFIED (CA to `u^68`; `rc_cyc` to `u^69`).
  * CA REFEREE_2 s.5 Remark, (CYC) for every self-adjoint operator of the same shape: SKETCH (the Wronskian step is not proved in
    general); VERIFIED on three deformations (CA).
  * PS s.4.2 (c), `A^QRR_n in Q[L, 1/L]`: a by-product of the glue; checked here (s.2 A4).
* **Human checks:** none of any step (part 3, s.8).

### 6.4 Nothing is REFUTED
No referee found a counterexample, no re-run differed, no label is lowered. (CYC), PB s.2.6-s.3 and LMI move from "PROVED (u)" to
"PROVED, AI-refereed"; the statement of (CYC) is strengthened (A1); LMI extends to `M0` in C (A3); SYN2 K5's argument is repaired (A6).
