# Stage 2 synthesis, part 2: falsifiers, contradictions, overclaims, next steps, critic record (2026-11-14)

Continues `STAGE2_SYNTHESIS.md` (s.0-s.2; same labels and conventions). AI-written; no human has read it.

## 3. Falsifier outcomes (SYNTHESIS_STAGE01 s.6, plus the parallel s.1.5 test)

| falsifier | what would have broken the route | outcome | independent re-derivation (referee) |
|---|---|---|---|
| (a) D(iii) beyond weight 240: `critic/cr_dictionary.py 80 320` and `100 400` | one log-free leg outside `L` | **PASSED**: all 10400, resp. 16200, log-free basis vectors lie in `L` (13 min, 77 min) | <ul><li>re-run of `80 320` is byte-identical;</li><li>own code finds none outside `L`: `w <= 400` exact, `w <= 1000` mod both primes (100500 vectors per prime);</li><li>detection power: a bad `H^2` leg gives 40 FAIL lines;</li><li>the test cannot fail once Lemma 2 (PROVED) holds.</li></ul> |
| (b) TR `omega_(4,2)`, `omega_(3,3)` paired with every L-leg and `-log p` | a non-zero pairing with an L-leg and negative v-exponent; a pure log pairing `!= F_h^(n)` | **PASSED**: 17416 singular ordered tuples per prime (4678 mixed), all 0; 18 pure-log values equal the Gap derivatives, including the new frozen `<H^2>_4 = 7/30000` (the model engine's missing amp(4,2)) | <ul><li>a different engine and basis, extended to `g + n <= 7`, `g <= 6`;</li><li>40739 singular multisets per prime (19194 mixed), 0 non-zero;</li><li>26 pure-log values, including `<H^3>_4 = -7/3750`, `<H>_5 = -1/82500`, `<H^2>_5 = 3/27500`, `<H>_6 = -691/102375000`;</li><li>all four new values re-derived here from `c_h` (s.7).</li></ul> |
| (c) B-universality | an edge graph singular, or the vertex law off, for a B-perturbed system | **PASSED**: S1-S5 at genus 2, 3 (both primes), S1, S2, S3, S5 at genus 4 (maxout 2), and S2 in full (378/378, 2^89-1); controls `stored`, `dress0`, `dress2`, `staleM` fail as predicted | <ul><li>the tables agree residue by residue;</li><li>S2 in full re-run on a module-disjoint engine;</li><li>S4, S6 (maxout 2) and S7 in full (378/378) at 2^61-1;</li><li>controls `staleHat`, `staleboth`, `nosig`;</li><li>**Lemma C**: B is a coordinate, so this tests `(F, Gm)`-universality plus covariance (K7).</li></ul> |
| (d) s.1.5 at genus 3 (flat worker) | flat `!=` hatted for some graph | **PASSED**: 42 graphs x `u^-4..u^3` agree (frame (0,0) at both primes, frame (1/3,-2/7) at 2^89-1); all 41 edge graphs regular in the flat frame; vertex `(1/25200) t_c^-4` | <ul><li>an independent red/blue engine;</li><li>R7: the composite terms regrouped by uncontracted graph equal the production table `A_id` on 252/252 non-zero coefficients;</li><li>controls fail.</li></ul> |

Critic checks of this stage (s.7): LD (M2) VERIFIED; LD on actual TR correlators VERIFIED; the count of K4 re-counted.

## 4. Contradictions between notes, and which side is right

* **K1. Orientation of the TR kernel against arXiv:math-ph/0702045.**
  * The claims that disagree:
    * LOGLEG s.2.3 says its kernel is that of arXiv:math-ph/0702045; TR_MODEL s.1 says the same of the kernel of arXiv:1509.06954;
    * the xyswap referee says arXiv:2212.00320 (1.1) is "minus the kernel of arXiv:math-ph/0702045";
    * the logleg referee read the printed paper: its (4-3)/(4-9) kernel is the OPPOSITE of the engines' kernel (its footnote 3), and
      its Thm 4.1 has `W_3 = +Res BBB/(dx dy)`, whereas the engines have `-`.
  * **Resolution (critic).**
    * The printed kernel `(1/2) Int_q^{q-bar} B/((y(q) - y(q-bar)) dx(q))` and arXiv:2212.00320 (1.1) with `omega_{1,0} = -y dx` are the same expression.
    * Hence "code = `(-1)^n` x arXiv:2212.00320" (VERIFIED, xyswap) and "code = `(-1)^n` x printed" (logleg referee) agree.
    * The phrase "the kernel of arXiv:math-ph/0702045" in TR_MODEL and in the xyswap referee therefore means the common convention
      `(1/2) Int_{sigma z}^{z} B/((y(z) - y(sigma z)) dx)`, not the printed (4-3).
    * The logleg referee is right about the printed paper. I checked this only on the formulas as quoted in the notes, not against
      the PDF.
  * **No result is affected.** Every sign is fixed by regressions: the psi-class intersection numbers on the Airy curve, `F_2^TR = -1/240`, the case
    (0,2) -> (0,3) re-proved in LOGLEG s.2.4, and tensor-level checks of the variational formula.
* **K2. The sign of the xi-term in (V1).**
  * LOGLEG s.0 says the critic's sketch had `-`; its referee (F6) says that is not established.
  * **Both are right.** The sketch reads `C'(v) - (polynomial)`, with the polynomial coming from the `Lie_xi` terms, and
    `<Phi, Lie_xi w> = -<(xi dPhi)_+, w>` (LOGLEG s.2.2(a)). So with polynomial `:= Sum <.., Phi_i Lie_xi, ..>` the sketch equals (V1).
  * Worker and referee VERIFIED the "+" form with `(xi dPhi)_+`; it is the precise statement.
* **K3. Does Cor. 4 use D(ii)?**
  * DICT_PROOF s.6 says no; its referee says yes.
  * **The referee is right**: with `K_0` in the statement, Cor. 4 needs `ker(log map) = K_0`. (d) is D(ii)-free.
* **K4. "112 non-zero regular coefficients of edge graphs"** (S15 s.0) against the flat referee's 108.
  * **The referee is right.** Re-counted from the worker's own tables (`cr2_count.log`): 108 edge-graph + 4 vertex-graph coefficients,
    identical in five tables (flat at both primes, flat in the second frame, hatted K8 at both primes).
* **K5. Control X4: 103/336 (worker) against 105/336 (referee).**
  * Without the sign, the edge quotient is not symmetric, so the count depends on edge orientation.
  * Neither number is a property of the mathematics; both runs fail, as required. **The referee is right.**
* **K6. `c0 = 1/lambda`** (BVAR) against `|c0| = m^10/|f(0) g(0)|` (bvar referee).
  * The two agree for `m = 1`, `g(0) = -1`, which covers S1-S4 and SYNTHESIS_STAGE01 C10's `lambda = f(0)|g(0)|`.
  * **The referee's formula is the general one**: S5 gives `2^10 = 1024`, and S6 gives `59049/2`.
* **K7. What falsifier (c) tests.**
  * BVAR s.7: "a failure would have revealed a hidden use of the quintic's B". The referee's Lemma C says `B` is a coordinate.
  * **The referee is right.** Re-derived here: `prod b~ = phi'^5 B_X(phi) = -(v_X o phi)'^5/(v_X o phi) = -v_q'^5/v_q = B_q`, and
    each `b~_a = phi' b_a(phi)` keeps the self-dual shape.
  * So `(F, Gm)`-universality with full power series `F~`, `Gm~` is what was tested, together with covariance (`staleM`, `staleHat`,
    `staleboth`).
* **K8. Is s.1.5 an input?**
  * SYNTHESIS_STAGE01 s.4 lists it separately; S15 derives it from (F).
  * **S15 and its referee are right**: the bracket [s.1.5] is deleted.
* **K9. The window of (F).**
  * SYNTHESIS_STAGE01 s.4 says "in the window of genus g"; the flat referee says "through `z^(3g-3)`". Consistent; the second is sharper.
  * Genus 3 sees only `c_5`. Symplecticity makes `log c` odd, so `c_10 = c_5^2/2` is forced; the next independent content is `a_3`
    (`z^15`, genus `>= 6`) and then `a_5` (`z^25`, genus `>= 10`).
* **K10. Sign conventions of the TR legs.**
  * XYSWAP s.0 uses `tau_0(1) <-> -z` with step `tau_k = +Int tau_(k-1) dx`. LD uses `tau_0(1) <-> +p` with step `-Int dx`.
  * Both pairs satisfy the string and dilaton equations, so the TR string and dilaton checks cannot tell them apart.
  * Only LD is the dictionary dual to the functions `xif_{a,d}` (TR_MODEL's string pairing `<y, d xif_{a,d}> = delta_{(a,d),(0,0)}`;
    s.2 step 7.2).
  * No XYSWAP result depends on these signs, since they concern membership in `L` and the log legs. **LD is the dictionary of record.**
* **K11. Status of D(iii).**
  * LOGLEG (written in parallel) says "VERIFIED to weight 240"; DICT_PROOF says PROVED. This is a matter of timing.
  * **D(iii) is PROVED** (refereed). The logleg Corollary needs it in the form "`Q[v][p] + Q(v)[x]`", which is exactly D(iii)(b).

## 5. Overclaims (quote, then correction)

* **dict.**
  * O1. DICT_PROOF s.1, "Rescaling by the `c_b` changes no statement below". This holds for the SET of legs, but not for `K_0` as a
    subspace (`e = tau_1(1) + v tau_0(H)` lies in `K_0`, while its raw leg `-p^5/5 + 2v log p` is not log-free). The "16200 `K_0` vectors"
    of T7 are vectors of the raw-leg kernel.
  * O2. "`dim K_0 cap S_w = d(w)` by Lemma 1(d) and (a)". It follows from (d) or Cor. 4, together with D(ii).
  * O3. s.6, "(d) and Cor. 4 do not use D(ii)". This is true for (d) only (K3).
  * O4. s.6, "Exactly what is used". The list omits LD (now proved, s.2 step 7.2) and the closed form of `PiC`. The theorem is sensitive
    to the `H^2` leg: `Phi_0^2 = a p^4 + b q` works iff `a != b/4`.
  * O5. s.9, "this closes ML_K restricted to `K_0` insertions modulo E9". It also needs M1, LD, M6 and the `Q[v]` coefficients (M5).
    All of these are now in place.
  * O6. s.1, "D(ii) (PROVED, SYNTHESIS s.6)". The Stage-1 argument was a sketch; the proof of record is dict referee s.3.
* **logleg.**
  * O7. LOGLEG s.2.3, "our kernel is" that of arXiv:math-ph/0702045, cited as its "s.5.1". The kernel is opposite to the printed (4-3), and the temperature
    variation is in s.5.3. The "form used" is a correct translation only through the referee's dictionary (K1).
  * O8. (E1) is listed as evidence for (Q1). It is equivalent to `<p^5>_h = 0`, an instance of L-regularity (F2). (F2 works in the
    worker's favour: the base `h >= 2` of (ii) needs no (Q1).)
  * O9. LOGLEG_2 s.3, "by homogeneity ... `Phi_i in Q[v][p]`". This is a gap, because homogeneity allows `v^-j p^(w+5j)`. It is now
    closed by D(iii)(b), which is PROVED.
  * O10. "The critic's sketch had `-`" (K2).
  * O11. "This ... leaves D(iii) as the only unproved model-side input of ML_K". The same holds for LD and M5 (both proved now), and the
    status line should list arXiv:1509.06954 Thm 8.1.
  * O12. FALSIFIER_B lists "functions of x pair to exactly zero" among the controls. It is a consistency check, not a failing control.
    "Every (h,n)" hides that the test is vacuous at (0,3) and (1,1) (F4).
  * O13. The Gap Theorem (Q3) is quoted without its convention (F5). LOGLEG s.2.4's sign of `omega_{0,3}` is opposite to Thm 4.1 of arXiv:math-ph/0702045
    and should say so.
* **flat.**
  * O14. S15 s.0 item 2, "112 non-zero regular coefficients of edge graphs". The correct figure is 108 + 4 (K4).
  * O15. s.5 item 1, "an independent graph-sum-level confirmation of (N) at `z^5`". It is consistent with E4 at `z^5` and adds nothing new.
  * O16. s.6, X4 "103/336" is quoted as a fixed count. It depends on orientation (K5).
  * O17. F1, "flat = stored hatted record (168 coefficients)". 164 of these are zeros. The non-zero comparison now rests on R2 and R7.
  * O18. s.1.3, the leaf formula `T' = z - z Q'_0(-z)`, is wrong at `z^1` (it should be `T' = z y + q_{R'}`). Lemma 3 cites "s.3" for
    homogeneity (should be I1) and never defines `x`. The u-adic convergence (`d(0) = 1`) is not stated.
  * O19. "Genus 4 out of reach". This is a limit of MVertex6 only; the referee's engine has no deficit limit.
  * O20. s.0 item 1 suggests (F) "to all orders". At genus g only `z^(3g-3)` is used (K9).
* **bvar.**
  * O21. s.7, "a failure would have revealed a hidden use of the quintic's B". Weakened by Lemma C (K7).
  * O22. s.8, "B-universality in all genera: CONJECTURAL" as a separate conjecture. It is equivalent to `(F, Gm)`-universality and is now
    PROVED in flat form modulo (Q) and Theorem U (s.2 Corollary).
  * O23. s.0, "`c0 = P_11(0,0) = 1/lambda`" (K6). The S5 table entry "`lambda = 1/1024`" is a relabelling, not `f(0)`.
  * O24. s.1, the converse paragraph (pole orders `k = 2..13`) is labelled PROVED. It is SKETCH, and not needed: the exclusion of
    `ord_0 B != -1` already follows from `ord_0(v'^5/v) = 4n - 5`.
  * O25. s.6(ii), "genuinely depends on B". It depends on the model coordinate.
* **Stage 1 (the critic's own; now superseded).**
  * O26. SYNTHESIS_STAGE01 E11(ii) and s.6 item 1, and MECHANISM s.1 item "D(ii) (PROVED, critic)". The "LG dictionary" step was one
    sentence, and the leg dictionary itself was assumed. Proofs of record: dict referee s.3 (D(ii)) and s.2 step 7.2 here (LD).
  * O27. SYNTHESIS_STAGE01 s.0, "Hence ML restricted to `K_0` IS the xyswap worker's L-regularity". This also needed LD, D(iii) and M5.

## 6. Next steps

**6.1 Stage 4: the connection problem.** This is the only open mathematical input besides Theorem U's human check.
* *Statement.* Prove `D(z) = 1`, where `Q D(z) = P_can Q^m(v(u))` (UR s.4).
  * Here `Q` are the WKB columns of the formal quintic at the conifold, normalised by the QRR constants at `q = 0`.
  * Under (H), `D = c(z) 1` with `c in Q[[z^5]]`. Symplecticity gives `c(z) c(-z) = 1`, so `c = exp(a_1 z^5 + a_3 z^15 + a_5 z^25 + ..)`.
  * E4: (H) holds, and `a_1 = a_3 = 0`, through `z^24`.
  * `a_(2j+1)` first matters at the genus `g` with `3g - 3 >= 5(2j+1)`: `a_5` at genus 10, `a_7` at genus 13.
* *Reformulation* (PROVED from gauge s.3, Cor. (b), given (H)). (N) holds iff the actual gauge `P^act = Q (Q^m)^-1` satisfies (G),
  i.e. `P^act Y^C(v(u)) = Y^C_Q`: a statement about ONE solution, the conifold vanishing period.
  * Heuristic (SKETCH; only meaningful after Borel summation, not formally): the vanishing cycle has the same thimble
    coordinates in the quintic's QRR-normalised WKB frame as in the model's stationary-phase frame. For the model, the residue loop
    around the pole is a signed sum of the five thimbles.
  * Suggested route (SKETCH):
    * compute both coordinate vectors in closed form, the model one from the residue loop at the pole, the quintic one from the
      QRR constants and the conifold period;
    * an oscillatory-integral representation (UR s.8) or a Gamma-type integral-structure argument is the natural tool;
    * allowed: quintic Picard-Fuchs data and KP/Toda; NOT TR on the quintic's mirror curve.
  * Total potentials cannot decide (N): the stored normalisation leaves the total unchanged (MECHANISM s.2).
* *Evidence task.* Extend E4 to `z^34` (`a_5`) from from-scratch rows. This would make N(g) for `g <= 12` conditional only on (U) and (Q).

**6.2 Human checks** (priority order; all AI-written so far).
1. LD (s.2 step 7.2, 15 lines), D(ii) (dict referee s.3, one page), D(iii) (DICT_PROOF s.2-s.5): elementary algebra.
2. Theorem V (LOGLEG s.2-s.3) with the sign dictionary of the logleg referee (s.2): the analytic core. Check Thm 5.1 of arXiv:math-ph/0702045 (temperature variation)
   and the differentiation under the contours.
3. Quintic Theorem R and theory Lemmas 0, E, T, F (theory s.1.3, s.3.1-s.3.2).
4. Theorem S15 (flat s.2-s.4), in particular Lemma 3 against the definition of the hatted outer sum (NATURAL_FORM, UR).
5. Theorem U (UR s.3), Theorem G and Prop. S (gauge s.2-s.4).
6. Hypotheses of the quoted theorems:
   * arXiv:1509.06954 Thm 8.1 (case (8-1), `n = 4`, at `a = (0,0,0,0,v)`);
   * arXiv:2212.00320 Thm 7.1 (`y = z`, rational `x` with simple critical points);
   * arXiv:math-ph/0702045 Thm 5.1;
   * arXiv:2502.18115 Cor. 1.2 (optional);
   * arXiv:1303.1043 s.2.

**6.3 Write-up plan.** A LaTeX note in a new folder `allgenus/note/`, built from the corrected statements of s.1, s.2 and s.5:
* §1 Statement of Conjecture N and of the conditional theorem, with the input table of s.2.
* §2 The two cyclic systems, Theorem U, the canonical gauge, (G), symplecticity.
* §3 The vanishing tower, the slice `K`, frames; transport Lemmas E, T, F.
* §4 Quintic Theorem R.
* §5 The model as TR: arXiv:1509.06954, LD, D(ii); the lattice `K_0 = L` (D(iii)).
* §6 L-regularity, Theorem V, the Gap Theorem; ML_K.
* §7 Flat versus hatted (S15).
* §8 The main theorem; the connection problem (H) + (N) and E4; genus `<= 9`.
* §9 Universality and Lemma C; the failing controls as necessity examples.
* Appendices: (A) sign conventions (the kernel-orientation dictionary, `H = d/dv`, `c_b`, `m5`); (B) the computations (ranges, primes, logs);
  (C) errata to the working notes (s.5).

Estimated 45-60 pages; the style of `resolvedconifold/hodge_vertex_conifold/note/` can be reused.

**6.4 Further evidence (optional).**
* An independent hatted engine for `k >= 4` at genus 3 (flat referee C4).
* A genus-4 flat run with a contraction-based engine (C6).
* A (Q1)-free proof of (V). The dilaton route bounds the pole by `v^-1`, so only exponent `-1` remains to be excluded; log-TR swap
  (arXiv:2312.16950) or hypermaps/KP are candidates.
* Falsifier (b) at `g + n = 8`.

**6.5 Errata and housekeeping.** These are for the owners of the notes; the critic edits no existing file.
* O1-O25.
* MECHANISM.md s.1, s.3, s.4: D(ii), D(iii) and (V) are PROVED; the bracket [s.1.5] is removed; ML_K is PROVED modulo the quoted theorems.
* Absolute scratch paths containing the local login name: in nine `flat/logs/*.log` and in `bvar/logs/g4_batch.sh`.
* The bvar pool briefly ran 5 processes (disclosed).

## 7. Critic record (new files in `critic/`; `py -I -B`, stdlib; working directory in the scratch folder; nothing else written)

| file | content | run | result |
|---|---|---|---|
| `cr2_legcore.py` | exact legs (`Phi_(k+1) = -Int Phi_k dx`), closed-form `PiC` (Lemma V), `x`-polynomial reduction; residues at the five branch points mod `P` via principal parts of `D^d xif_a` | imported by `cr2_legdict` and `cr2_e2e` | - |
| `cr2_legdict.py` | L1 LD; L2 branch independence; L3 D(ii) over Q; L4 frames end to end; controls C1-C4 | `cr2_legdict.py 24`, 75 s, log `cr2_legdict_24.log` | ALL CHECKS PASSED: 0 mismatches in 35000 pairings per prime (2^61-1, 2^89-1; `v = 1, 3^5`); D(ii) 125/125; 52/52 frames per prime; C1 60, C2 25, C3 25, C4 1385 mismatches |
| `cr2_e2e.py` | LD on actual `omega_{g,1}` (allgenus/tr engine via `logleg/ll_core.py`, read-only), exact over `Q(zeta5)` | 15 s, log `cr2_e2e.log` | 52 frames, `g = 2, 3`: 0 mismatches against TR_MODEL's exact one-point table; every singular one is 0 |
| `cr2_count.py` | recount for K4 from the flat worker's scratch tables | `cr2_count.log` | 108 + 4, in five tables |

* Recomputed by hand from `c_h = |B_2h|/(2h(2h-2)5^(h-1))`:
  * `F_4' = -1/30000`, `F_4'' = 7/30000`, `F_4''' = -7/3750`;
  * `F_5' = -1/82500`, `F_5'' = 3/27500`;
  * `F_6' = -691/102375000`.

  These are the referee's new pure-log values.
* Not re-run: the workers' and referees' long computations. I relied on their logs (all PASS) and on the referees' independent re-derivations.
* No process is left running, and no `__pycache__` was created.
