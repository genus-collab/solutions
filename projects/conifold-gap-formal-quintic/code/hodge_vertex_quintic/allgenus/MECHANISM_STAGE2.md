# All-genus mechanism for Conjecture N: Stage 2 (2026-11-14)

Entry point for Stage 2. It supersedes the statuses in `MECHANISM.md` s.1, s.3 and s.4 (Stage 0-1), which are kept unchanged as the
record of 2026-11-13. Written from a 9-agent workflow (4 workers, 4 referees, 1 critic and assembler) and from my own re-runs (s.4).
Everything is AI-written and AI-refereed: no human has checked any proof. Labels as in `../NATURAL_FORM.md`.
Full record: `stage2/STAGE2_SYNTHESIS.md` (s.0-s.2: bottom line, the Model Lemma ingredient by ingredient, the assembled argument) and
`stage2/STAGE2_SYNTHESIS_2.md` (s.3-s.7: falsifiers, 11 adjudicated contradictions, 25+ overclaims with corrections, next steps).
Worker notes: `stage2/dict/DICT_PROOF.md`, `stage2/logleg/LOGLEG_PROOF*.md` and `FALSIFIER_B.md`, `stage2/flat/S15_PROOF.md`,
`stage2/bvar/BVAR.md` and `BV_TABLES.md`; referee reports in `stage2/*/referee/REFEREE.md`.

## 1. The assembled statement

**Conjecture N(g) holds for every g >= 2, modulo exactly:**
* **(U)** Theorem U (`natural_form/theory/UNIFORM_REDUCTION.md` s.3): proved in the project, AI-audited only;
* **(F_g)** the normalisation of the actual R-matrix through `z^(3g-3)`: in `Q D = P_can Q^m(v(u))` the constant diagonal series `D` is
  `1 mod z^(3g-2)`; equivalently (H) and `c^act = 1` through `z^(3g-3)`. For `2 <= g <= 9` this is the finite identity E4, VERIFIED
  (two primes, two frames, exact over Q on 397 Fourier components) but not proved; for `g >= 10` it is OPEN (the next independent
  coefficient `a_5`, at `z^25`, first matters at genus 10);
* **quoted theorems**: arXiv:1509.06954 Thm 8.1; arXiv:2212.00320 Thm 7.1; arXiv:math-ph/0702045 Thm 5.1 with the linear loop equation of
  its s.4; arXiv:1303.1043 Def. 2.2, Prop. 2.4, Def. 2.5, Prop. 2.9; optionally arXiv:2502.18115 Cor. 1.2.

The chain (STAGE2_SYNTHESIS s.2): vanishing solutions -> (G) -> symplecticity -> transport (edges in `K_0 x K_0`, leaf point in `-z + K`
with `Pi_0 = t_c`) -> quintic Theorem R -> Theorem S15 (flat = hatted, graph by graph; derived from (F), so the s.1.5 identification
is no longer a separate input) -> the Model Lemma on K for all genera -> Conjecture N(g). Genus 2-5 remain VERIFIED directly by graph sums.

## 2. What Stage 2 proved (each refereed; corrections applied in the synthesis)

* **D(iii)** for all weights: the legs of `K_0` lie in `Q[v][p] + Q[v][x]`, and `K_0` maps bijectively onto the leg class L modulo
  functions of x (DICT_PROOF). VERIFIED exactly for weight <= 400, at two primes for weight <= 1000.
* **D(ii)** for all k, now with a complete argument (dict referee s.3); the Stage-1 "PROVED" label was a sketch.
* **LD**, the leg dictionary `q_{Phi_k^b} = -c_b tau_k(H^b)`, `c = (-1,1,-1,-3,-2)`: Stage 1 used it without proof; PROVED by the
  critic (STAGE2_SYNTHESIS s.2 step 7.2).
* **Theorem V**, the log legs: pairings with L-legs and `-log p` legs lie in `Q[v]`, pure log pairings equal `F_h^(k)`. PROVED modulo
  arXiv:math-ph/0702045 Thm 5.1 (via the referee's sign dictionary), arXiv:2212.00320 Thm 7.1 and the Gap Theorem.
* **The Model Lemma on K** for every genus h >= 0: PROVED modulo the quoted theorems; nothing about the quintic enters.
* **Theorem S15**: the hatted outer term with `M = M_nat` equals the flat `P_can` term graph by graph in every genus, given (F).
* **Lemma C**: a change of `B = prod b_a` is a change of coordinate `u = phi(w)`, so B-universality is coordinate covariance.

## 3. The four falsifiers (all passed; nothing REFUTED)

* (a) D(iii) beyond the Stage-1 range: `critic/cr_dictionary.py 80 320` and `100 400`: 10400 and 16200 log-free vectors, all in L.
* (b) TR correlators `omega_(4,2)`, `omega_(3,3)` and the mixed L / log pairings: no non-zero pairing with negative v-exponent among 17416
  singular tuples per prime (two primes; exact over `Q(zeta5)` in the low cases); pure log pairings equal the Gap derivatives.
* (c) B varied: canonically normalised B-perturbed systems give regular edge graphs and the vertex law at genus 2-3 (two primes) and
  genus 4 (at most two edges for four systems; all 378 for two of them); failing controls fail.
* s.1.5 at genus 3: flat and hatted sums agree graph by graph (336 coefficients per prime and frame; 108 non-zero regular coefficients
  of edge graphs + 4 of the vertex graph); five deliberately wrong conventions all fail.

## 4. My own re-runs (2026-11-14)

* `stage2/critic/cr2_e2e.py` (52 frame one-point functions at genus 2-3, exact over `Q(zeta5)`): identical to `cr2_e2e.log`.
* `stage2/critic/cr2_legdict.py 24` (the leg dictionary, 35000 pairings per prime, two values of v, failing controls): identical to
  `cr2_legdict_24.log` (74 s).
* `critic/cr_dictionary.py 80 320` (D(iii) falsifier to weight 320): 10400 log-free basis vectors, all in L, controls fail as they
  should; identical to `stage2/dict/cr_dictionary_80_320.log` (12 min 28 s).
* Not re-run: the workers' long computations (TR correlators, genus-3 flat sums, B-perturbed genus-4 runs); the referees re-derived those
  with independent code.

## 5. What is left

1. **Stage 4, the connection problem**: prove `D = 1` (equivalently: the actual gauge maps the model's vanishing solution to the
   quintic's). This is the only open mathematical input besides Theorem U. Evidence task first: extend E4 to `z^34` from rows computed
   from scratch, which would make N(g) for g <= 12 conditional only on Theorem U and the quoted theorems. No total-potential argument
   can decide it (the stored normalisation leaves the total unchanged). TR on the quintic's own mirror curve stays excluded.
2. **Human checks**, in this order: LD, D(ii), D(iii) (elementary algebra); Theorem V with the sign dictionary; quintic Theorem R
   (Lemmas 0, E, T, F); Theorem S15; Theorem U, (G), symplecticity; the hypotheses of the quoted theorems.
3. **Write-up**: a separate LaTeX note (the critic suggests `allgenus/note/`, 45-60 pp) built from the corrected statements.
4. Housekeeping: the critic's corrections O1-O25 (STAGE2_SYNTHESIS_2 s.5) are not yet applied to the worker notes; 23 Stage-2 logs and
   one batch script echo absolute scratch paths that contain the local login name.

## 6. Re-run log

All four re-runs of s.4 finished on 2026-11-14 and agree with the stored logs. The re-run outputs stayed in the session scratch
folder; nothing was written into the project except this note.
