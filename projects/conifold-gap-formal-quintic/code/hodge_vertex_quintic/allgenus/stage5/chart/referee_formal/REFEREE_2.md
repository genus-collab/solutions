# Referee report, part 2: colours, Theorem D, corrections, code and runs (2026-11-14)

Continuation of `REFEREE.md` (s.0 verdict, s.1-s.3 line-by-line check of CHART_PROOF s.1-s.3.7). Same labels; AI-written.

## 4. Colours and the Z/5 covariance (CHART_PROOF s.0.2, s.3.9): CONFIRMED

* `P_k(zeta L) = zeta^-k P_k(L)`. The Taylor coefficients of `S_0`, `S_1`, `S_2k` at `zeta L` are `zeta^-k`, `zeta^-k`, `zeta^(1-2k-r)` times
  those at `L`. With `t = zeta t'` ((L1a), constant Jacobian) the phase becomes `zeta` times the phase at `L`, and the normalised
  `e^(S_1)` part is unchanged. `z^(2k-1) S_2k` becomes `(zeta^-1 z)^(2k-1) S_2k(L + t')`. (F3) with `mu = zeta` then gives
  `c(z; zeta L) = c(zeta^-1 z; L)`. Every step checked.
* *Colour `i`.* `q(zeta^i L) = q(L)`. The colour-`i` point `sigma = zeta^i L` of the parameter `L` is the colour-0 point of the
  parameter `zeta^i L`, the logarithm of `-q` being chosen to make it critical (s.4 (B1)). So
  `c^(i)(z; L) = c(z; zeta^i L) = c(zeta^-i z; L)`.
  Moreover `w_i = zeta^i w_0` gives `h_(w_i)(z) = h(zeta^-i z)`. So (Gamma) at colour `i` is equivalent to (Gamma) at colour 0.
* Consequences: `c_n` is regular at the fifth roots of unity, `c_n in Q[L, L^-1]`, and its exponents are `= 4n (mod 5)`. CONFIRMED.
* Theorem D does not need the covariance: R0 already gives `D_i(z) = D_0(zeta^-i z)`. The note says so.
* VERIFIED (`rf_sigma.py 40 10`, (D)): `c_n(zeta^i L) = zeta^(-in) c_n(L)` exactly in `Q(zeta)` for `i = 1..4`, `n <= 10`,
  `L = 2, -1/3`. The SP code is my Lagrange-Buermann code, not the worker's Wick routine. Control `zeta^(+in)` fails.

## 5. Theorem D: the assembly from R0-R7 (s.0.3, s.3.10, s.3.11): CONFIRMED given its inputs

**The chain.** Each implication was checked against the cited statement.
1. Theorem Gamma (iii): `c_n(1) = h_n`. Proved in the note (REFEREE s.3).
2. s.3.10: `c_n(L) = L^(4n) ph^SP_n(L^-5)` as rational functions, by (L1a) for `sigma = L/X` and (F3) for `z = L^-4 w`. So `ph^SP_n` is
   regular at `u = 1`, by 1 (or by Theorem P(a)), and `ph^SP_n(1) = h_n`. CONFIRMED, with correction C2.
3. R5, Theorem P(b) (POLY s.6 with the poly referee's REFEREE s.3): `Phi^SP = Phi^can`. Hence `ph^can_n(1) = h_n` for all `n`.
4. R3 with R2 (EVID s.2.3-2.4): `Phi^can` and `Phi^QRR` are polynomial-type WKB columns, so `Phi^can = Phi^QRR d(z)`.
   * At `u = 1` (`s = 1`, `z = w`): `Phi^can(z, 1) = h(z) d(z)`.
   * So `ph^can_n(1) = h_n` for all `n` is equivalent to `d = 1`, to `Phi^can = Phi^QRR`, and to `D_0 = 1`.
5. R0: `D_i(z) = D_0(zeta^-i z)`, so `D = 1`, and (F_g) holds for every `g`.

**Inputs that are actually used, with their status (none is re-proved in the note or here; all are AI-written, none human-checked).**
* (U) Theorem U (UR s.3), AI-audited only. It enters through R0, R2 (Lemma N) and R5 (step 5 of P(b)).
* (R1) of UR s.5, used by Lemma N and by P(b).
* R3 needs (QRR) and, in EVID's formulation, the "u = 1 convention of ref_rows" (C4).
* R5 needs:
  * (U), GI Theorem G with GI Lemma 1, Lemma V, Lemma Q;
  * the model rows, the CFU map and the Brieskorn reduction;
  * (FL).
  The SP items of (FL) are now written out in s.2 of the note. Two items are not written out, neither in the note nor in the
  Stage-4 notes I read (SYN4 s.2 lists them as "standard, not written out"): the residue items of (FL), and Stirling's functional
  equation `Gamma_St(y+1) = y Gamma_St(y)` in `Q[[1/y]]`. Both are standard; the poly referee
  gives the contour-integral definition of the residue that makes them immediate. The note is honest about this (s.5 item 3).
* Not used, confirmed: R1 (Lemma P; Lemma N's proof uses only (U), formal WKB uniqueness and (R1)), R4, and R6 (Theorem A(i);
  step 7 of P(b) uses only the `u^1` term of the residue, `-u/z`).

**What Theorem D does NOT settle.** The external normalisation is (QRR): the actual R-matrix at `q = 0` is `diag(h_(w_i))`, including
its `z^(5k)` part. That part is exactly what (N) tests. Theorem Gamma computes `c(z; 1) = h(z)` with no reference to the actual
R-matrix. So "D = 1" is proved for the R-matrix normalised by (QRR).
* The overall sign and the twist weight `-5` are pinned: `h_1 = ph_1 = -3/20`, and `ph_1` is forced by HPF.
* The `z^(5k)` parts are not pinned by any such check.
* Nobody in Stage 4 or Stage 5 has located the twisted statement with a statement number. Givental's Thm 9.1, re-read by me on p. 12,
  is the untwisted equivariant statement. The worker says this plainly (s.5 item 2).
* Before (F_g) for all `g` is advertised as a theorem, the source of (QRR) should be pinned down: a twisted QRR statement with a
  number (e.g. in arXiv:math/0110142, or LP arXiv:1702.06096 / 1803.01409), or the genus-3 constant-map test of SYN4b s.6.2.

**Labels.** Theorem Gamma: PROVED (AI-written, now AI-refereed; modulo DLMF 5.11.1 and 2.3.8). Theorem D: PROVED given (U), (R1), R2,
R3, R5 with its inputs, R0 and (QRR). The note's own labels (s.0.5, s.0.3) say exactly this and do not overclaim. One sentence of
s.0.4 reads more broadly than it should: "No finite identity remains among the hypotheses". It is true, but should be read together
with the (QRR) remark above.

## 6. Corrections (all minor; none changes a conclusion)

* **C1 (s.3.2, last bullet; false as written).** "`det H_ch = -1/(L^2 Lam1(L)^2)` vanishes only at `L = 0, zeta, .., zeta^4`". It never
  vanishes; it has double poles there. Read: "`1/det H_ch = -L^2 Lam1(L)^2` vanishes only at `L = 0, zeta, .., zeta^4`. These are
  the `L` at which `p_ch(L)` leaves the chart (`y = 5L` or some `x_j = zeta^j - L` is 0). For every other `L`, `H_ch` is invertible,
  and `det H_ch(1) = -1/25`." s.3.4 and s.3.5 use the right statement, so nothing downstream changes.
* **C2 (s.3.10).** "The map `sigma = L/X` is analytic near `X = 1` ((L1a); its Jacobian is removed by normalisation)". The Jacobian
  `d sigma/dX = -L/X^2` is not constant. It belongs to the hatted amplitude: with the factor `X^2` of `e^(S_1(L/X))` it gives POLY's
  `(1 - u X^5)^(-1/2)`. Only its value at `X = 1` is removed by normalisation. The identity is unaffected (VERIFIED, `rf_hat` (H3)).
* **C3 (s.5 item 7; JSON "PROVED elsewhere per SYN4 s.1.1").** Polynomiality `c_n in Q[L]` follows from three facts together:
  Theorem P(a) (`Phi^SP` solves HPF), s.3.10, and the degree bound valid for every HPF solution (SYN4 s.1.1 and K2:
  arXiv:0710.0889 Thm 4(ii), or the evid referee's elementary proof). SYN4 s.1.1 alone does not give it. It is not needed.
  VERIFIED here (`rf_hat`): `deg ph^SP_n = floor(4n/5)` for `3 <= n <= 30`, and 0 for `n <= 2`.
* **C4 (s.0.3, R2 + R3).** EVID's Theorem 2.4 is stated "given Theorem U and the u = 1 convention of ref_rows". List that
  convention, or say that it is part of (QRR) as formulated in s.0.3 (the flat colour-0 value at `L = 1`). Bookkeeping only.
* **C5 (s.1.7, POLY bullet; optional).** The reflection argument is correct. The direct argument of REFEREE s.1 is shorter: the
  same `S''`, the same normalised amplitude, and `L` critical by the choice of logarithm.
* **C6 (s.2.8; optional).** The Pochhammer part of (L3) is formal, by (S1) with `h = x^M` (REFEREE s.2). The DLMF input is then
  needed only for the bare Stirling series.
* **Remark (s.1.7, `dx/x` versus `dx`).** At `L = 1` the series is blind to `dy` versus `dy/y` (REFEREE s.1). The worker's controls
  run at `L != 1`, so they are valid. Worth a sentence in the note.

## 7. My own code and runs (all new; `py -I -B`, stdlib, working directory in the session scratch folder; logs in `logs/`)

No import of `ch_*.py` or of any Stage-4 module. The algorithms differ from the worker's.
* One-variable SP by the Morse lemma and Lagrange-Buermann inversion (`rf_core.sp1`):
  `N_n = Sum (2k-1)!! C^k [t^(2k)](A_m r^(-(2k+1)/2))`. The worker uses Wick/Isserlis.
* Multivariable Gaussian expectations, two ways (`rf_wick`), neither of them the worker's Isserlis recursion:
  * a closed form for `C = diag + rank one` (`eta = eta_0 + u gamma`);
  * the differential-operator form `(1/k!) Delta_C^k`.
  The two agree (`rf_chart`, first line), and the second agrees with Lagrange-Buermann on random data (`rf_lemmas` (T0)).
* `c_n` at `L = 1` from a Laurent expansion with rigorous precision tracking, instead of interpolation.
* POLY's hatted integrand with `u` formal, instead of the sigma-form.

Primes: `7000000000000000061` and `2900000000000000000000291`. Both are `= 1 mod 10`, and Miller-Rabin with the bases 2..41 is a
proof below `3.317e24`. They differ from the worker's `2^62 + 177` and `2^80 + 85` (both checked prime here), from the critic's
`2^61-1` and `2^89-1`, and from the gamma referee's `2^127-1` and `2^107-1`.

| file (run, log) | what is checked | field, range | result |
|---|---|---|---|
| `rf_core.py`, `rf_rings.py` | modules: `Q(zeta)` (inverse via Galois norm), `F_p`, Bernoulli, Newton sums, `h`, `sp1`; Laurent series `LS` with absolute precision; `Q[[u]]/(u^n)` | - | - |
| `rf_sigma.py 40 10` (`rf_sigma_40_10.log`, 662 s; also `16 6`) | (A) `c_1 = -(3/20) L^4`, `c_2 = (9/800) L^8` at four `L`; stored `c_5(2) = 11495493197/31500`. (B) `c_n(1 + eps)`: no pole and `[eps^0] = h_n`. `h` from the five Stirling factors in `Q(zeta)` = Newton-sum `h`. (D) Z/5 covariance | Q; Q((eps)), `n <= 40`; `Q(zeta)`, `n <= 10` | ALL PASS. Controls: twist `+5` (`z^1`); `Gamma(1+5s)` Stirling term dropped (`z^1`); `b_1 + 1/7` in `j = 0` only leaves a pole `eps^-1`; `Gamma(5s)` changes `c(z; 2)` from `z^2` (blind at `L = 1`, PASS as predicted); `zeta^(+in)` |
| `rf_hat.py 30` (`rf_hat_30.log`, 49 s) | POLY's hatted SP in `Q[[u]]/(u^29)`: (H1) `deg ph_n <= floor(4n/5)`; (H2) `ph_n(1) = h_n`; (H3) `L^(4n) ph_n(L^-5) = c_n(L)` (`n <= 12`, `L = 2, -2/3, 3`); (H4) `ph_n(0)` = model SP; `ph_5` = GAMMA s.1's polynomial | Q, `n <= 30` | ALL PASS. Controls: twist `+5` and `(1 - uX^5)^(-1/2)` dropped both break (H1)+(H2) at `w^1` |
| `rf_hess.py` (`rf_hess.log`, 1 s) | s.3.2 in full (REFEREE s.3.2) | `Q(zeta)`, 8 `L`; polynomial identities by 40 points | ALL PASS. Controls: `+q`; wrong `det H_7`; `q -> -q` |
| `rf_chart.py 3 2 4 3` (`rf_chart_3_2_4_3.log`, 264 s) | chart SP = `c(z; L)`; chart SP at `L = 1` = `h`; `Omega_7` SP = `c(z; L)`; `K = L - 1` from `q`; covariance closed form asserted | exact `Q(zeta)`: chart `z^3` (`L = 1/2`, `1`), `Omega_7` `z^2` (`L = 1/2`); mod both primes: chart `z^4` (`L = -3, 5/3, 1`), `Omega_7` `z^3` (`L = -3`) | ALL PASS. Controls: no Jacobian (`z^1`); `dy/y` (`z^2`); `q -> -q` not critical; `Omega_7` with `dy/y` differs from `c` at `z^2` and equals the `Gamma(5s)` sigma-form |
| `rf_chart.py 2 1 5 3` (`rf_chart_2_1_5_3.log`, 606 s) | as above, with the chart through `z^5`, the first order that sees (N) | mod both primes: chart `z^5` (`L = -3, 5/3`: `= c`; `L = 1`: `= h`), `Omega_7` `z^3` (`L = -3`); exact `Q(zeta)` chart `z^2`, `Omega_7` `z^1` | ALL PASS; the same controls fail |
| `rf_lemmas.py 7 6` (`rf_lemmas_7_6.log`, 55 s; also `5 4`) | (T0) `sp1` = Wick; (T1) (L3); (T2) (L2); (T3) (L1b) with `beta`; (T4) (S1); (T5) (L1a) | Q (and `Q(zeta)` for (T1)), 6 seeds, `z^7` | ALL PASS; nine controls fail |
| `rf_eexp.py 6 12` (`rf_eexp_6_12.log`, 4 s) | step (e), 112 non-zero coefficients | Q, `q^6`, `z^-d..z^12` | ALL PASS. Controls: `L` dropped; `e^(Dpsi/z)` dropped; twist `+5` |

**The worker's code.**
* I read the Taylor-data code (`ch_phases.py`), which agrees with my construction.
* I re-ran eight scripts from a scratch COPY, because the scripts write their logs next to themselves: `ch_hess`, `ch_zpow`,
  `ch_sigma 10`, `ch_stirling 30`, `ch_euler 8 6`, `ch_lemmas 3 1`, `ch_eexp 2 3`, `ch_chart 1 1 2 2 2`. Every new log is identical to
  the stored one, timings stripped.
* The long runs (`ch_chart 3 4 6 5 7`, `ch_euler 40 40`, `ch_eexp 6 10`, `ch_lemmas 6 4`) were not re-run.
* File sizes (all under 250 lines) and the logs (no absolute paths) comply with the brief.

**My development record.** These bugs were found and fixed before the reported runs.
1. `geom()` in `rf_chart` used `comb(e+k-1, k)` for negative exponents. It crashed at the first smoke test.
2. I first expected the `Gamma(5s)` control to fail at `L = 1`. It does not (blindness, REFEREE s.1). It now runs at `L = 2` (fails
   at `z^2`), and the blindness at `L = 1` is recorded as a PASS.
3. The `Omega_7` `dy/y` control was first run at `z^1`, where it cannot fail. It now runs at `z^2`.
4. Two syntax slips: a dict with tuple keys passed as `**`, and quoting in one message.

The logs of the smoke runs (`rf_sigma 6 3`, `rf_hat 10`, `rf_chart 2 1 2 1`, `rf_lemmas 2 2`, `rf_eexp 2 3`) were my own throw-away
files and were deleted. No check other than a deliberate control ever failed.

## 8. Literature, re-checked by me

* DLMF 5.11.1: Stirling's series for `Ln Gamma(z)`, valid as `z -> oo` in `|ph z| <= pi - delta`. DLMF 5.11.3-5.11.4: the `g_k`, with
  `g_1..g_5 = 1/12, 1/288, -139/51840, -571/2488320, 163879/209018880`, as quoted in the note.
* DLMF 2.3.8 is Watson's lemma. Its hypothesis on `q(t)` is (2.3.7), with `lambda, mu > 0`, and its conclusion is
  `Int_0^oo e^(-xt) q(t) dt ~ Sum Gamma((s+lambda)/mu) a_s x^(-(s+lambda)/mu)`. Correctly quoted and correctly used (`mu = 2`,
  `lambda = 1`).
* arXiv:math/0108100v2, pp. 12-15: Thm 9.1's `b_i(z) = Sum_k N^(i)_(2k-1) (B_2k/2k) z^(2k-1)/(2k-1)`; s.10, Thm 10.1, Example 10.3;
  the chart `x_0 = Q/(x_1..x_n)` with the Lagrange multiplier `P`; Lemma 10.4; and "This coincides with `b_0(z)` in 9.1 since
  `lambda_alpha - lambda_0` are exactly the infinitesimal characters of the torus `G` action on `T*_(w_0) CP^n`". The note's s.3.8
  reports all of this accurately. On these pages I did not see the word "outline" that the note quotes; it may be on pp. 1-11, which
  I did not read. It is irrelevant to the proof.
* Not read: arXiv:0710.0889 eq. (5); arXiv:1702.06096 eq. (5); arXiv:1803.01409; Hoermander ALPDO I Thm 7.7.5 (remark (F4), used
  in no proof).

## 9. What I did not check; compliance

* Not re-checked (inputs of Theorem D only): Theorem U, GI Theorem G and its lemmas, Lemma N, the CFU map, the Brieskorn reduction,
  the residue items of (FL), Stirling's functional equation, (QRR). I read Theorem P(b) and the connection form only to the extent
  needed for the chain of s.5.
* Compliance:
  * New files only, all in `stage5/chart/referee_formal/`; every file is under 250 lines.
  * No existing file was modified, moved or deleted. The worker's scripts were run only from a copy in the scratch folder.
  * At most three of my processes ran at a time, and none is left running.
  * Logs print relative paths only.
  * Web: the DLMF pages s.5.11 and s.2.3, and the arXiv PDF of math/0108100v2. No personal data in any request.
