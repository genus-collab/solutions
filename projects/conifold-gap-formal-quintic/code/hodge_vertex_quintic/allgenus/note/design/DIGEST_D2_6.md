# DIGEST_D2, part 6: s.5 Evidence, s.6 Literature, s.7 Open issues

Continuation of `DIGEST_D2_5.md`. Script and log names below are for traceability only; the note's evidence appendix states content,
ranges, fields and primes, never file or script names (rule 5). All numbers as stated by the cited reports (I re-ran nothing; I only
re-computed `h_1..h_5`, the `+5` control `-11/60` and `[z^5] Ginf = 4537/3937500` from the definitions).

## 5. Evidence for the appendix (finite checks; every one with failing controls in its source)

| # | what | range; field or primes | source |
|---|---|---|---|
| E1a | (Gamma) directly: Laurent expansion of `c_n` at `L = 1 + eps`, no negative power of `eps`, `[eps^0] c_n = h_n` | `n <= 60`, exact Q | GEN s.4 T1 (`g5_t1 60 60 Q`, 750 s); re-run M5 s.4.1 |
| E1b | the same with rigorous precision tracking; control `b_1 + 1/7` in the `j = 0` term leaves `eps^-1` | `n <= 40`, exact Q | RF2 s.7 (`rf_sigma 40 10`) |
| E1c | global fit of `c_n` at 32 rational points (Laurent window with slots for poles), `c_n(1) = h_n` | `n <= 30`, exact Q | LF2 s.5.2 F |
| E1d | two-variable chart at `L = 1` equals `h` computed in `Q(zeta_5)` from its definition; `det H = -1` | `z^60`, exact | LF2 Q3 |
| E1e | (Gamma) by a scaled Wick sum in `eps` | `z^18` exact Q; `z^24` mod `999999999999999631`, `999999999999999999661` | RA2 s.4.3 C1-C2 |
| E1f | unscaled Laurent Wick sum with tracked precision | `z^12` exact; `z^16` mod `4611686018427387761`, `1180591620717411303161` | LA2 s.4.4 Q1 |
| E1g | two-variable chart = one-variable `c(z; L)` at `L = 2, -2/3, 1/2, 3`; `= h` at `L = 1`; `c_1..c_3(2)` | `z^16`, exact Q | IC (`lmi_check 16`, 201 s); M5 s.4.2 |
| E1h | Stage 4: shape fit, `c_n(1) = h_n`; `c = R` | `z^30` exact; `z^500` mod `2^127-1`, `2^107-1`, `z^40` exact | M4 s.4.1-s.4.2 |
| E2 | 5-variable chart = `c(z; L)`: exact `Q(zeta)` `z^3` (`L = 2, -2/3`), mod `2^62+177`, `2^80+85` `z^6` (`L = 3, 5/4, -1/2`); chart at `L = 1` = `h` (exact `z^7`, mod p `z^6`); 7-variable SP = `c` (exact `z^4` at `L = 2`, mod p `z^5` at `L = 3`); controls: no Jacobian (`z^1`), `dy/y` (`z^2`) | as stated | CP4 s.6 (`ch_chart 3 4 6 5 7`, 2137 s) |
| E2b | the same with new Gaussian-expectation algorithms: chart mod `7000000000000000061`, `2900000000000000000000291` through `z^5` (`L = -3, 5/3`; `= h` at `L = 1`), exact `z^3`; 7-variable `z^3` mod p, `z^2` exact | as stated | RF2 s.7 (`rf_chart`) |
| E2c | the Stage-4 critic's 5-variable test: `L = 1` chart = `h`; `L = 2, 3/2, -2/3` chart = `R_n(L) = L^(4n) ph_n(L^-5)`; `z^5` sees (N) | `z^6` mod `2^61-1`, `2^89-1` | GEN2 s.6 (T5) |
| E3 | two-variable chart = one-variable series | exact `z^9` (RA2 s.4.3, `ra_2var 9`); `z^14` at six `L` (LF2 Q2); coefficientwise in `z` and `e` (LA2 S1) | as stated |
| E4 | values: `c_1 = -(3/20) L^4`, `c_2 = (9/800) L^8`; `c_1..c_5(2) = -12/5, 72/25, 354278/625, -73040834/3125, 11495493197/31500`; `c_1..c_3(-2/3) = -4/135, 8/18225, 4067314/110716875`; `h_1..h_5` (D3) | exact | RA2 s.4.2 E1-E2, C0; LA2 E1; CP4 s.6 (`ch_sigma`) |
| E5 | Z/5: `c_n(zeta^i L) = zeta^(-in) c_n(L)`, `i = 1..4`; colours 1, 2 with an explicit `zeta` in `F_p` | `n <= 10` exact `Q(zeta)` (CP4 `ch_sigma`; RF2 `rf_sigma` (D)); `z^24` mod p (RA2 C3); `z^16` mod p (LA2 Q6) | as stated |
| E6 | orbifold point: `c(z; 0) = Ginf` (chart at `L = 0`, `det H = -25`: `z^60`; value step `z^80`; Method A `z^14`; fit `m <= 6`); `[z^5] Ginf = 4537/3937500`; controls `theta = 0, -1/5` fail at `z^5`, `theta = 4/5` only at `z^10`, `M0 = 0` at `z^5` | exact | LF2 Q4, D3, M2, F |
| E6b | `c(z; 0) = Ginf`: `z^18` exact (RA2 O1); `n <= 20` exact, unshifted tail fails at `z^5` (GEN2 s.5.1 A); `z^12` exact, `z^16` mod p (LA2 Q4); Stage 4 `m <= 100` at two primes, `m <= 25` exact (SYN4 s.1.1) | as stated | as stated |
| E7 | shape: `c_n in Q[L]` with exponents `4n - 5k in [0, 4n]` (`n <= 30`, LF2 F; `n <= 8`, RA2 E4; `n <= 16`, CP4 `ch_eexp`); `deg ph^SP_n = floor(4n/5)` for `3 <= n <= 30`, 0 for `n <= 2`, `ph^SP_n(1) = h_n` (`n <= 30`), `c_n(L) = L^(4n) ph_n(L^-5)` (`n <= 12`) | exact | RF2 s.7 (`rf_hat 30`), RF-C3 |
| E8 | LMI on general data: `g = 0, F = 1` gives `c^F = 1` (`M0 = -2..3`, `z^12`); ten random cases (`alpha0 in {-1, 2, -1/3, 5/2, 5, -7/4}`, `M0 in {-1,0,1,2}`), four entire cases (`z^14`), non-integer shifts `1/5, 1/3, -1/2, 7/4` (`z^14`) | exact | LF2 L2, L3, M3, M4 |
| E8b | twelve random cases incl. `M0 = 1/3, -2/5, 7/2`; four mod p cases (`alpha0 = sqrt(-1)`, generic `M0`); `g''(s*) = 0`; RA's 36 modified integrands (P1-P10, I1-I4) | `z^10` exact, `z^12`, `z^14` mod p; `z^18` exact | LA2 s.4.2 P; RA2 s.4.3 |
| E9 | (L2'): random `P`, `A`, `alpha0`, `c in {0, +-1, 2, -3, 1/5}`; moment table vs Wick recursion (degree `<= 24`); controls: wrong sign of `delta`, pure `m^3` term, `s^2 m` term | `z^8` (RA2 s.4.5); `z^12` (LF2 L4) | as stated |
| E10 | falsifiers (GAM-15 (c)): RA2 s.4.4 X1-X8 (`z^18`); LA s.3 with LA2 s.4.3 X1-X10 (`n <= 12`, 37 checks; square: `1/16, 9/512, 75/8192, 3675/524288`; X8: `1, 2, 6, 24, 120, 720` at `alpha0 = -1`); GEN2 s.5.3 G, H1, H2 (`n <= 20`; cube: `1/9, 2/81, 14/2187, 8/19683, -440/177147`) | exact | as stated |
| E11 | generalisation (GAM-16): twisted `P^4`, weights `W1 = (1,3,-2,5/2,-7/3)`, `W2 = (2,-1,3,-5,1)`, `W3 = (1,1/2,1/3,1/4,1/5)`, `W4 = (7,-3,4/3,-11/2,2/5)`, all 20 colours (`n <= 24`); 12 theories, 40 colours (`n <= 20`); GEN2 s.5.1 B-E, s.5.2 F-F3 | exact Q | GEN s.4, GEN2 s.5 |
| E12 | formal lemmas on random data: CP4 s.6 (`ch_lemmas 6 4`: (S1), (L1a), (L1b) inner dim 1-2, (L2), (L1c)(iii); `ch_euler 40 40`: (L3), `m = -3..3`, both measures, `a = 1, -2/7, 3/2`, `z^40`); RF2 s.7 (`rf_lemmas 7 6`); FSL5 s.11.2 (Fubini to `z^8`/`z^6`, changes of variables, delta lemma, T-q, T-z, T(b), FL4, (L3) to `z^80` at `a = 1`, residue); RL2 s.7 (a 4-variable toy chart through `z^6`, toy (Gamma) through `z^14`; (L3) by Lagrange-Buermann, `m = -6..6`, `z^40`; FL4 and T(b)); LF2 L1, D1-D2 (formal 5.11.8 through `w^-60`); RS3 s.7.2 (L3-C with symbolic `M`, `z^12`) | exact Q | as stated |
| E13 | step (e): all 98 non-zero coefficients `q^d z^n`, `d <= 6`, `-d <= n <= 10` (CP4 `ch_eexp 6 10`); 112 coefficients to `z^12` with `c_n` from the hatted polynomials (RF2 `rf_eexp 6 12`) | exact Q | as stated |
| E14 | QRR sign: from LP (26) mod Mumford's relation, `N_(2,0) = -19/9600`, `N_(3,0) = 1/224000`; weight `+5`: `-97/28800`, `1/81000` | exact `Q(zeta_5)` | LIT s.1.5 |

*Caveats for the appendix text.* A control can be blind at one parameter value (RL2 s.7.4: the toy chart at `L = 2` through `z^9`);
the reported controls were run at generic points. RL-C8: one residue check is close to a tautology; the lemma workers' data are
polynomial phases (log-type data in RL2 s.7). Two kept draft logs (LA `DRAFT1_lm_run_x_4`, with a removed sub-claim; PS `rps_run_v2`) are
not evidence (RS s.0). Re-runs: M5 s.4.1 (42 first-workflow logs identical; 46 referee logs, 44 identical, 2 with the same verdicts);
RS3 s.7.1 (26 identical); RF2 s.7, RL2 s.6. Independent methods: Wick/Isserlis (CP), Lagrange-Buermann (RF, GEN), scaled Wick in `eps`
(RA), unscaled Laurent Wick (LA), separated two-variable engine (LF), new main-session code (IC).

## 6. Literature (statement numbers as read by the sources; mode of reading; key in the root `refs.bib` if present)

| item | statements | use | read by, how | root key |
|---|---|---|---|---|
| Givental, arXiv:math/0108100v2 | Thm 9.1 (untwisted; `R|_(Q=0) = exp diag(b_i)`, cotangent characters; proved by citation, Rem. 9.2(a)); Thm 4.1 (Hodge potential), Rem. 4.2; s.10: Thm 10.1, Rem. 10.2, Ex. 10.3, Lemma 10.4, chart `x_0 = Q/(x_1..x_n)`, Thm 10.5, Cor. 10.6 | prior art (GAM-17); QRR-3 (Thm 4.1); not a (QRR) source | CP (pp. 11-15), RF2 s.8 (pp. 12-15), LIT2 s.3.1 (own PDF extraction) | `MR1901075` (Mosc. Math. J. 1 (2001)) |
| Givental, arXiv:math/0008067v4 | Thm 1 (s.2.2) and the Remark after it; s.2.3; Example 6 | the source GJR cite for (QRR) | LIT s.1.3 (own PDF extraction) | none (new key) |
| Coates-Givental, arXiv:math/0110142v2 | Thm 1 (s.3, eq. (7)), Remark (2), s.4 eq. (9), Cor. 1, Cor. 3, s.7 eq. (11), Thm 2, eqs. (12)-(13), App. 1 | (QRR) as cited by GJR; QRR-2 (b) sign; prior art (twist half, genus 0) | LIT s.1.3, LIT2 s.3.2 (PDF and ar5iv) | `MR2276766` (Ann. of Math. 165 (2007)) |
| Guo-Janda-Ruan, arXiv:1812.11908v1 | s.3.1 eq. (9); s.4.1; Lemma 5.2 proof, eq. (17); s.7: Thm 7.1, Prop. 7.5, Cor. 7.6, Lemma 7.7 eq. (37), Lemma 7.9, Rem. 7.2, Cor. 7.3-7.4; Lemma 5.5 | (QRR) (QRR-1); prior art: Rem. 7.2 (unproved analogue of (Gamma3)), Lemma 7.7 (`Ginf` for their LG integral) | LIT s.1.2, LIT2 s.3.6 (own PDF extraction; ar5iv truncated) | none (new key); "only version, no journal reference" |
| Guo-Janda-Ruan, arXiv:1709.07392v1 | eq. (7); Lemma 4.6 proof, eqs. (19)-(20) | (QRR) at `t = 0` | LIT s.1.2 (own PDF extraction) | none |
| Lho-Pandharipande, arXiv:1803.01409v2 | s.0.2 eq. (5) (definition; footnotes 3-4); s.1.1; s.2.4.4 eq. (18), Prop. 5, eq. (21); s.3.2 eq. (26); Props. 16, 21-22 | definition of the theory; QRR-3 (vertex class) | LIT s.1.1 (PDF, in full) | none |
| Lho-Pandharipande, arXiv:1702.06096 | s.0.2 eq. (5); s.0.3 eq. (6); s.0.6 Thms 3-4 | definition (same formula) | LIT s.1.1 (ar5iv, summariser) | none |
| Teleman, arXiv:0712.0160v3 | Thm (1) s.1.7; Thm (2) s.1.9; Prop. 8.3, Cor. 8.4, Prop. 8.5 (ii) s.8.2 | uniqueness of R (QRR-4) | SYN5 s.4.3 only, ar5iv through a summariser, NOT verbatim; LIT did not open it | `MR2917177` (Invent. Math. 188 (2012); v3 agrees with it, SYN5 s.4.3) |
| Faber-Pandharipande, arXiv:math/9810173 | s.1.1 Thm (Mumford, unnumbered), Props. 1-2, Thm 4; s.2.1 Mumford's relation (unnumbered), eq. (22) | QRR-3; E14 | LIT s.1.3 (ar5iv, summariser) | `MR1728879` (Invent. Math. 139 (2000)) |
| Mumford 1983 | (as quoted by FP, Givental, CG) | QRR-3 | not read | `MR717614` |
| Pandharipande, arXiv:1712.02528 | Def. 5; example after Thm 7; footnote 11 | convention `R` vs `R^-1` | LIT s.1.3 (ar5iv, summariser) | none |
| Zagier-Zinger, arXiv:0710.0889v1 | eq. (1); Lemma 1; eq. (5); Thm 3 eq. (6); Thm 4 (i)-(iii) | GAM-12 remark; normalisation conventions; (Thm 4(ii): Stage-4 degree bound) | LIT2 s.3.5 (own PDF extraction); CP, RF did not read eq. (5) | none |
| Iritani, arXiv:0903.1463v3 | Thm 4.14 (= 1.2), Conj. 4.3, eqs. (70)-(71), Rem. 4.15, Thm 4.17, Rem. 4.18, Lemma 4.20, Lemma 4.21 | prior art (Lemma 4.20) | LIT2 s.3.3 (own PDF extraction) | `MR2553377` (Adv. Math. 222 (2009)) |
| Galkin-Golyshev-Iritani, arXiv:1404.6407v4 | Thm 1.3.1 (= 5.0.1), eq. (5.2.3), Prop. 5.2.4, Cor. 5.2.6 | set aside (leading order, large `t`) | LIT2 s.3.4 | none |
| Coates-Iritani, arXiv:1203.4193 | - | set aside | LIT2 s.3.7 [H] | `Coates:2012vs` |
| DLMF | 5.11.1, 5.11.3-5.11.4, 5.11.8; 2.3.8 (with 2.3.7); s.2.3(iii) (2.3.13)-(2.3.18); s.2.1(iii) (2.1.13)-(2.1.16); 24.2.3, 24.2.5, 24.4.1; 6.12.1; 5.2.1 | Proof B of (L3) (remark); 24.4.1 in SP-5118; 6.12.1 only as an explanation (E10) | CP, RF, RA, FSL, RL on the site; LF, LA via fetch summaries/renderings (not verbatim) | none |
| Etingof arXiv:2409.03117v2; Kirwin arXiv:0810.1700; Polyak arXiv:math/0406251; Bar-Natan-Garoufalidis-Rozansky-Thurston arXiv:math/9801049 | Etingof Thms 2.16-2.18, Cor. 2.19, Thms 3.1, 3.5, 3.8, 3.11, Rem. 3.6; Kirwin Thm 1.1, Cor. 3.6; Polyak Props. 2.2, 2.4, 2.17, 2.20, Rem. 2.3, Thms 2.7-2.8; Aarhus II Defs. 2.7-2.8, Ex. 2.11, Props. 2.13, 2.15 | optional "related formal Gaussian integration" citations; none states the formal lemmas for complex data with numbers (LIT s.0) | LIT2 s.2 [H] (summariser; numbers confirmed by a second query) | none |
| Hoermander, ALPDO I, Thm 7.7.5 | - | CP's remark (F4) only | NOT read (cited from memory, CP s.1.2; LIT3 s.5) | - |

Root keys found by searching `refGW/refs.bib` (read-only) for titles and authors: MR1901075, MR2276766, MR1728879, MR2917177,
MR2553377, MR717614, Coates:2012vs; also MR1408320 (Givental, "Equivariant GW invariants", IMRN 1996), which LIT3 s.5 lists among the
papers cited by GJR and Givental that were not opened (the sources do not record whether it is the "[15]" of GJR s.7 to which LIT2
s.3.6 says GJR attribute the mirror integrals). The root keys are journal (MR) entries; the statement numbers above
were read in the arXiv versions named (check the published numbering, `\draftnote`).

## 7. Open issues for the writers

1. **Main proof.** Theorem Gamma by LMI (GAM-4, GAM-3) with LF-C1..C5 and LA-C1..C6 applied; the 7-variable chart as a remark (GAM-11,
   with RF-C1) (M5 s.5 item 2; SYN5b s.8 item 2; RS3 s.8). The note's numbering of Theorem Gamma's parts differs from CP's (GAM-5).
2. **Appendix base.** FSL with RL's Lemma H and RL-C1..C8 (SYN5b s.8 (vii)); shorter CP proofs where available: (L2) by the degree count,
   (PG) algebraically, the sharp (F1) bound. Use descriptive names (no "FL2", "(U)").
3. **Fully formal route.** (L3) by Proof A and `Ginf` by SP-5118: Theorem Gamma then uses no analytic theorem; DLMF enters only as
   24.4.1 (one line from 24.2.3) and in the remark on Proof B.
4. **GAM-2 / the definition of `Phi^SP`.** RF-C5's direct route lacks one sentence (the `(log z)/z` cancellation for POLY's form); the
   reflection route supplies it. Decide with the Theorem-P digest which integrand defines `Phi^SP` in the note; if it is the hatted
   pull-back of D6, GAM-2 is not needed.
5. **Dependencies on Theorem P(a).** GAM-7 (polynomiality of `ph^SP_n`) and GAM-8 (`deg c_n <= 4n`) cite P(a): coordinate the order of sections.
6. **(QRR)** (`\draftnote`s): GJR (17) read only in a text extraction, the sign of `-5 h_alpha` and `R = e^(-u/z) S C_alpha^-1` to be read in
   the typeset PDF (SYN5b O15, s.7 item 5); the uniqueness step via Teleman (QRR-4) is not literally stated in what was read; the
   `z^(5k)` parts are not pinned by any check and the genus-3 code test was not done (QRR-5). State (QRR) exactly as in QRR-1, used only as H9.
7. **Prior-art wording** (`\draftnote`s): Givental's "outline" (CP3 s.3.8 says the text calls itself one; RF2 s.8 did not find the word on
   pp. 12-15; LIT2 s.3.1: the abstract presents the Virasoro result as an outline); GJR (37) factor `m`, Rem. 7.2 sign, display of `W`
   (LIT3 s.5); "no paper read proves (Gamma)" stated modestly (GJR may regard Rem. 7.2 as standard via Givental's mirror theorem,
   LIT3 s.4 item 3). Do not claim that GJR's LG series equals `c(z; L)` (LIT2 s.3.6).
8. **Step (e) and (Gamma3).** The coefficientwise formula (GAM-12) is proved; its identification with ZZ eq. (5) is written nowhere
   (CP3 s.3.7): omit or `\draftnote`.
9. **Unrefereed items used only in remarks:** GEN Lemmas 2-4 (GEN had no referee); LMI for `M0 in C` and L3-C (checked only by the RS
   critic). Transcribe their written proofs if used; mark with `\draftnote`.
10. **Measure of the singular factor.** It matters (LA s.2.4: `Gamma(1 + (1 - sigma)/z)` gives the value with `delta = z`, which differs
    from `h` at `z^1`); the summand has `Gamma(zeta^j/z - s)` (measure `dx_j/x_j`, `M0 = 0`) and `Gamma(1 + 5s)` (`dy`). Say this in the
    set-up, together with GAM-14.
11. **Placement.** s.1: (QRR) as input (QRR-1, QRR-5) and the conditionality; set-up: D1-D3, D6 (or in the Gamma section); the section
    "the connection problem" / "the value at `q = 0`": GAM-1..GAM-17; appendix "formal stationary phase": SP-*; evidence appendix: s.5.
    Statuses, referee counts, correction tags and file names go to DRAFT_NOTES.md only.
12. **Naming clashes** (s.4, last item) for NOTATION.md; in particular rename the Wick functional `N` and the DLMF shift `h`.
13. **Bibliography.** New keys for: math/0008067, 1812.11908, 1709.07392, 1803.01409, 1702.06096, 0710.0889, 1712.02528, 1404.6407
    (if cited), DLMF, and the optional formal-integration papers; reuse MR1901075, MR2276766, MR1728879, MR2917177, MR2553377, MR717614.
