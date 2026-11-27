# Stage 4 synthesis, part 2: contradictions, overclaims, next steps, critic record (2026-11-14)

Continuation of `STAGE4_SYNTHESIS.md` (s.0-s.3: bottom line, status of (H) and (N), the assembled statement, agreement of the routes).
Same labels and notation. AI-written; no human has read it.

## 4. Contradictions between notes, and which side is right

* **K1. Which paper defines the formal quintic.**
  * gamma s.6: "arXiv:1702.06096, cited in the brief, is the local P^2 paper; the formal quintic paper is arXiv:1803.01409".
  * evid referee: the formal-quintic twist is eq. (5) of 1702.06096. gamma referee (C7): its abstract defines the formal quintic.
  * Critic, from the arXiv abstract pages (2026-11-14):
    * 1702.06096 (Lho-Pandharipande, "Stable quotients and the holomorphic anomaly equation") has local CP^2 as its main subject. It
      also defines the formal quintic, via the Euler class of a complex, and proves its holomorphic anomaly equations.
    * 1803.01409 ("Holomorphic anomaly equations for the formal quintic", same authors) is the dedicated paper.
  * **Verdict:** the brief's citation is correct, and the project's `FORMAL_QUINTIC_genus2.md` cites 1702.06096 s.0.2. gamma's remark
    is an overstatement. Cite both.
* **K2. Is the degree bound `deg ph_n <= floor(4n/5)` proved?**
  * VERIFIED and "unexplained": evid s.2.1(b'), s.8.3; poly Lemma 3.
  * PROVED: gamma s.1 (arXiv:0710.0889 Thm 4(ii)); evid referee s.4 (elementary); poly referee C1.
  * **Verdict:** PROVED, by two independent proofs. Equality (`n = 1` and `n >= 3`) is VERIFIED only.
* **K3. Does Theorem A(i) need Theorem U?** gamma: "given Theorem U". gamma referee (R-a): no. A finite exact computation of `Y(0)`
  (degree `<= 2` in the frame parameters, decided on a 3 x 3 grid) together with GI Lemma 1 suffices.
  **Verdict:** the referee; the stronger statement holds.
* **K4. The name "(A2)".** poly's "(A2)" is the identity (Gamma); gamma's "(A2)" is the analytic uniformity step. This is not a
  mathematical conflict: the poly referee showed (A2)_poly <=> (A2)_gamma modulo gamma's (A1). This synthesis writes (Gamma).
* **K5. Is there still a connection problem?**
  * Workers: evid s.8.1, "no analytic connection problem in hatted variables"; poly s.0 and s.6, "(A2) ... a statement at q = 0 alone,
    no longer involves the conifold".
  * Referees: evid referee item 7 and poly referee C4 say it is still a connection statement, linking `u = 0, 1, oo`.
  * **Verdict:** the referees are right in substance. (Gamma) compares a series built at `L != 1` with its continuation to `L = 1`. The
    workers are right that the gauge, the periods and the model have dropped out.
  * The chart argument (part 1, s.1.4) settles the connection by a change of chart in which the coefficients are analytic at `q = 0`.
* **K6. Is the derived phase UR s.8's `W`?**
  * gamma s.0 item 2: the oscillatory-integral picture of UR s.8 is "now derived". gamma referee (C6): it is a different phase with the
    same critical data.
  * **Verdict:** the referee.
  * Critic: gamma's 1-form and poly's 1-form ARE equal. `W_gamma = -W_poly`, and the amplitudes agree by
    `Sum_j (zeta^j x - 1)^(1-2k) = -p_(2k-1)(x)`.
  * Numerical confirmation: the hatted SP at `u = 1/32` reproduces gamma's flat value `c_5(2) = 11495493197/31500` (`cr4_run` E3).
* **K7. Are other frames independent evidence?**
  * evid s.4-s.5 count frames (1/3,-2/7) and (7/11,-3/5) as further evidence for (N). evid referee (result A): given (U), `K_q` is
    frame independent, so those runs cannot fail.
  * **Verdict:** the referee.
  * Critic probe C4: even a non-frame vanishing period `g + f_2/7` gives `K = 0`. `KE` sees `C1` only through `C1(0)`.
* **K8. Does Theorem A(ii) have a direct test beyond `z^24`?**
  * gamma: none. gamma referee (R-b): `c^fS = 1` through `z^100`, with the gauge computed.
  * **Verdict:** the referee.
  * Critic, directly: `Phi^SP = Phi^QRR = Phi^can` at three values of `u`, exact through `z^60`, with `f^can` computed using no QRR
    input (`cr4_run` D1, E2).
* **K9. Perturbed values of `K_q`.**
  * The evid referee's `rf_pfree` prints `K_2 = 6048257071/24192000000` for "`f_5q + 1/(q+2)`". The critic's first run gave
    `737147/2953125`.
  * This is not a contradiction: there are two perturbation rules.
    * Referee: `ph_5q(1) = h_5q + 1/(q+2)`, with later `f` recomputed.
    * Critic's first run: `f_q = f_q^QRR + 1/(q+2)`, fixed numbers.
  * Both values are predicted exactly by `K_q = [z^5q] Phi(z,1)/h(z)`: `K_2 = 1/4 - h_5/3`, respectively `1/4 - f_1/3`. All three
    formulas agree in both cases.
  * The referee's wording should name the rule.
* **K10. GAMMA_PROOFS s.5(c), "absorbed by D(0) = 1".** gamma referee (G1): the factor `(2 pi z/5)^(1/2)` depends on `z`. It cancels
  because it is common to both stationary phases. **Verdict:** the referee. The conclusion `D = c/R` stands.
* **K11. Ranges for (F_g).** The stated ranges differ:
  * `g <= 24` (poly) and `51` (poly referee);
  * `67` and `167` (gamma and its referee, modulo Theorem A);
  * `69` and `169` (evid, exact and two primes);
  * `335` (evid referee, two primes).

  They do not conflict: the methods and hypotheses differ (part 1, s.1.5). Headline: exact `g <= 69`, given (U).

## 5. Overclaims (quote, then correction)

* **O1** EVID s.0, s.4 and its JSON: "(F_g) is VERIFIED (exact) for every genus g <= 69" (and "g <= 169" mod p).
  Correction: "VERIFIED given Theorem U". Lemma N uses (U) as an existence theorem (evid referee 7.1).
* **O2** EVID s.4-s.5: the second and third frames as further confirmation of (N). Correction: they are code and convention checks
  (K7).
* **O3** EVID s.2.1(b'), s.6, s.8.3 and POLY Lemma 3: the degree bound is "VERIFIED, not proved", "unexplained". This is an
  underclaim: it is PROVED (K2).
* **O4** EVID s.2.1(b') and the `ev_struct.py` docstring: "exponents in (-4n/5, 0]". Correction: `[-4n/5, 0]`.
* **O5** EVID s.8.1: "No analytic connection problem in hatted variables". Correction: there is one (K5).
* **O6** GAMMA s.0 item 3: "Theorem A (PROVED given Theorem U)". Correction, by part:
  * (i) is PROVED without (U).
  * (ii) is PROVED given (U), GI Theorem G and (FL). gamma's version has the bookkeeping gaps G2-G3. poly's Theorem P(b), with the poly
    referee's s.3, closes them.
* **O7** GAMMA s.0 item 2: the picture of UR s.8 is "now derived". Correction: a different phase with the same critical data (K6).
* **O8** GAMMA s.6: "arXiv:1702.06096, cited in the brief, is the local P^2 paper". Correction: it also defines and treats the formal
  quintic (K1).
* **O9** GAMMA s.3, s.8 and JSON open issue 3: "Theorem A(ii) has no direct numerical test ... beyond z^24". Superseded (K8).
* **O10** GAMMA s.8: "No purely ODE-formal argument can replace s.6-s.7". Correction: true for the properties listed. Theorem A and
  the chart argument are formal, but their input is the integral.
* **O11** GAMMA s.6: "arXiv:1803.01409 ... s.2.4.4 eq. (18) and Prop. 5" are listed as "read for this note". No one has confirmed
  this: the gamma referee and the critic saw the abstract only. Treat the details as unverified.
* **O12** POLY s.0 item 6, s.6 and s.7: "(A2) ... no longer involves the conifold", "a statement at q = 0 alone". Correction: it is a
  connection statement; the gauge, the periods and the model have dropped out (K5; poly referee C4).
* **O13** POLY s.0 item 6, s.7: "(A2) VERIFIED ... via Theorem P and item 3, through w^70 (exact)". This mixes proof and computation.
  Correction: (A2) follows through `w^70` from Theorem P together with the VERIFIED `c^act = 1` (poly referee C5).
* **O14** POLY s.0 item 3, s.4: "56 identities", "14 identities". Correction: 28 and 7 are independent; symplecticity forces the
  even ones (poly referee C2).
* **O15** POLY s.3: "divisibility of `g_N` by `u^mu` is the regularity of P". Correction: necessary, not sufficient (poly referee C7).
* **O16** POLY s.4: "(was g <= 9)". Out of date: evid reaches `g <= 69` (poly referee C9).
* **O17** The brief (KEY STRUCTURE): "D = 1 should be equivalent to ... the value at u = 0 (or the coefficient that the model fixes)
  required by the model". Correction, by all three workers: the model fixes nothing at `u = 0`.
  * The free coefficient is the top coefficient, i.e. the value at `q = oo`.
  * The conifold fixes it through the gauge (Lemma N), and the full-Stirling integral fixes it through Theorem P.
* **O18** A general caution: every "PROVED" of Stage 4 is AI-written and AI-refereed only. The chart argument of part 1, s.1.4 is not
  refereed at all.

## 6. Next steps

**6.1 Stage 5: prove (Gamma) by the chart argument (part 1, s.1.4).**
* Write out (L1)-(L3), with statements and proofs, or with citations giving statement numbers:
  * multivariable formal SP: changes of variables, Fubini, parameters;
  * the formal delta lemma;
  * Stirling's series as the formal SP of Euler's integral, with `x^m` inserted.
* Write steps (a)-(e) in full. Include the bookkeeping of powers of `z` and of the Gaussian prefactors: both orders give
  `(2 pi z)^(7/2)`, and the `z`-powers of the Euler substitutions cancel because `Sum zeta^j = 0`.
* Referee it. Optionally extend `cr4_chart.py` to `z^6` (estimated 20-30 min). The identity itself is already VERIFIED to `z^500` by
  the one-variable route.
* Combined with R5 this gives `D_0 = 1` to all orders. Conjecture N(g) for all `g` then rests on (U), (Q), (QRR) and standard formal
  lemmas only.

**6.2 Human checks, in priority order.**
1. Theorem P(b) (POLY s.6) together with the poly referee's bookkeeping (REFEREE s.3, `kappa = C0`): about two pages.
2. The chart argument, once written (Stage 5).
3. (QRR): LP's definition of the formal quintic and the `q = 0` R-matrix `diag(h_(wt_i))`, in particular its `z^(5k)` part
   (Mumford's GRR for the Hodge bundle).
   * A cheap independent test: compute the genus-3, degree-0 constant-map number of the formal quintic from the project's R-matrix and
     compare it with the Hodge-integral formula. A `z^5` change of the `q = 0` normalisation would generically alter it; the
     unpointed genus 2 (dimension 3) is blind to it.
4. Lemma N and the connection form (EVID s.2.3-2.4). These are now needed only for the finite statement `g <= 69`.
5. The Stage-2 list (STAGE2_SYNTHESIS_2 s.6.2), beginning with Theorem U.

**6.3 Write-up.** In the planned `allgenus/note/`, a section "The connection problem" with these parts:
* (i) the polynomial structure (part 1, s.1.1);
* (ii) Lemma N and the connection form;
* (iii) Theorem P, canonical = full-Stirling;
* (iv) (Gamma) and the chart argument;
* (v) the evidence tables and the failing controls;
* (vi) `Ginf` as a remark, labelled CONJECTURAL.

Estimated 10-15 pages, on top of the 45-60 planned in STAGE2_SYNTHESIS_2 s.6.3.

**6.4 Housekeeping**, for the owners of the notes (the critic edited no existing file):
* evid: O1-O5.
* gamma: O6-O11 and the sentence of K10. `gm_compare_20_200.log` was produced by a copy in the session scratch folder (disclosed in
  GAMMA s.7).
* poly: O12-O16. `pp_run_sp_N25_K22.log` ends with shell `time` lines, which is harmless.
* evid referee: name the perturbation rule of `rf_pfree` (K9).
* The `osc` task was not delivered. If it is re-run, its natural brief is now Stage 5 (6.1).

## 7. Critic record

New files in `critic/`, all under 250 lines. Each was run with `py -I -B` and stdlib only, with the working directory in the session
scratch folder. No project or worker module is imported. Logs print relative paths only.

| file | content |
|---|---|
| `cr4_core.py` | field Q or F_p; truncated series. Bernoulli numbers; QRR `h_(w_0)` by Newton sums of `y^4 - 2y^3 + 2y^2 - y + 1/5`; `tau`. Conifold periods: PF4 with dual exponents, so `gamma = 7/10` is found from the log condition, not assumed. Model rows `F_b`, closure asserted. Colour-0 (HPF) recursion with modes qrr / zero / list / can. The three functionals `KA`, `KE`, `KP` |
| `cr4_sp.py` | hatted full-Stirling SP at fixed rational `u` by Morse / Lagrange-Buermann. Phase from `dW/dX = -(1/(uX^2)) log((1-uX^5)/(1-u))`; amplitude by Newton sums of `(1-x^5)Y^5 - 5Y^4 + ...`; controls |
| `cr4_run.py` | the checks of part 1, s.3 |
| `cr4_chart.py` | the 5-variable chart integral of part 1, s.1.4: Taylor data at the colour-0 critical point, Wick's theorem with covariance `-H^-1`, mod `2^61-1` and `2^89-1` |

Runs:
* `cr4_run.py 4 6 20`: 2 s. `8 12 40`: 25 s. `12 25 60`: 180 s. Each ends ALL CHECKS PASSED.
* `cr4_chart.py 2`, `3`, `4`, `5`: under 1 s, 5 s, 36 s and 179 s. Each ends ALL CHECKS PASSED. Logs `logs/cr4_chart_<N>.log`; at
  `N = 5` the stored-normalisation line shows that `z^5` sees (N).
* Throw-away smoke tests stayed in the session scratch folder; nothing else was written into the project.

One wrong expectation was corrected during development: the control "measure `dy/y`" first fails at `z^2`, not `z^1`. The `z^1`
coefficient is blind to it, as is the whole series at `q = 0`, where `Gamma(1+y) = y Gamma(y)` changes only the elementary factor.

Not re-run: the workers' and referees' long computations, such as `z^1000` mod p and `z^500` (Gamma). The critic relied on their logs
(73, all PASS) and on the referees' independent re-derivations, and reproduced every number that is printed in both a note and the
critic's own range: `f_1, f_2, f_3`, `h_5`, `m5`, `c_5(2)`, `t_c` through `u^4`, and the referee-A perturbation values.

Compliance:
* No existing file was modified; no `__pycache__` was created.
* At most 2 processes ran at once; every process the critic started has finished.
* The only web access was the two arXiv abstract pages of K1.
