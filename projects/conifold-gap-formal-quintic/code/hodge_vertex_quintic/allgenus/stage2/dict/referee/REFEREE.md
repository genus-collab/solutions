# Referee report: Stage 2, "D(iii) for all weights + falsifier (a)" (2026-11-14)

Adversarial AI referee of `../DICT_PROOF.md` and of the scripts and logs in `..`. No human has checked this report or the proof.
Labels as in `../../../MECHANISM.md`: PROVED, VERIFIED, CONJECTURAL, SKETCH, REFUTED. Conventions of `../../../SYNTHESIS_STAGE01.md` s.3.
My code: `rf_*.py` in this folder (`py -I -B`, stdlib, exact over Q or mod 2^61-1 and 2^89-1). It imports nothing from `../dict_*.py`
or `critic/`, except `rf_crctrl.py`, which imports `critic/cr_dictionary.py` read-only in order to test the critic's own falsifier.

## 0. Verdict: SOUND WITH CORRECTIONS

* **Theorem D(iii) (a)-(d) is PROVED as written, given D(ii).** I checked every line of s.1-s.5 and s.9.
  * Lemmas 0-3, the graded Nakayama step and the proofs of (b), (c), (d), Cor. 4 and Cor. 5 have no gap.
  * The corrections (s.2, listed in s.6) concern attributions, one misleading sentence about the constants `c_b`, and the
    hypotheses of the s.9 consequence. None of them changes a statement of the theorem.
* **D(ii) in the form the proof uses: PROVED here for ALL k (s.3).** The proof is a generating-function lemma plus one residue
  identity. It needs only Lemma V of `gauge/GAUGE_IDENTITY.md` (PROVED there) and the Stage-1 definition of `PiC` through `Y^C`.
  So D(iii) no longer leans on the sketch-level "LG dictionary" sentence of SYNTHESIS s.6 item 1.
* **Falsifier (a): PASSES, and I could not break it.**
  * My re-run of `critic/cr_dictionary.py 80 320` (unchanged copy, SHA-256 equal to the original) is byte-identical to
    `../cr_dictionary_80_320.log`. My re-runs of `../dict_verify.py 160 60` and `../dict_xcheck.py 100` are byte-identical to the logs.
  * The critic's test has detection power. With the H^2 leg replaced in memory by `v^2/p^6`, it prints 40 FAIL lines, 39 of them
    "log-free combination NOT in L" (`rf_crctrl.py`).
  * Independent code finds no log-free leg outside `L`: exact over Q for `w <= 400`; mod both primes for `w <= 1000` (s.4).
  * Structurally the falsifier cannot fail. It follows from Lemma 2 (every single leg lies in `M`), which is PROVED.
    I VERIFIED the per-leg statement for `k <= 200` over Q and `k <= 1000` mod both primes.
* **All VERIFIED claims of the worker were re-derived with my own code** (s.4), with larger ranges and 13 failing controls.

## 1. Files (this folder)

| file | content |
|---|---|
| `rf_core.py` | homogeneous Laurent-log calculus; legs by iterated integration AND by a closed generating-function formula; normal forms; `in_Qpq`, `in_L`, `in_M`; rank, kernel |
| `rf_check.py` | R0-R7: re-derivation of T1-T7, X1-X4 and of Theorem D(iii) (a)-(d), Cor. 4, Cor. 5 |
| `rf_dii.py` | D1-D3: D(ii) in the used form (s.3); K_0 from the closed-form `Pi_i` against the log-free kernel of the legs |
| `rf_controls.py` | K1-K7: failing controls (bad legs, log-dropping integration, L-test anchors, perturbations, `Q(v)` sharpness) |
| `rf_crctrl.py` | the critic's falsifier with a bad H^2 leg (detection power) |
| `rf_counts.py` | reproduces every count quoted in DICT_PROOF s.8: 861/859 (T2, worker's range), 1042 (T3), 231 (X3), 2640/10400/16200 |
| `rf_h2leg.py` | sensitivity of (d) to the H^2 leg: `Phi_0^2 = a p^4 + b q` gives a bijection iff `a != b/4` |
| `rf_nvn.py` | the target of Lemma 3: `dim Q[p,q]_w/(Q[v,x] + pq Q[p,q])_w = [w >= 1]`, `q^n = -p^(4n)/4^n` mod it (linear algebra in `Q[p,q]`) |
| `rf_*.log` | run logs (all under 250 lines) |

## 2. Line-by-line check of DICT_PROOF.md

**s.1 (set-up).** Correct.
* Legs and weights are as stated.
* "A different constant convention changes legs only by elements of `Q(v)[x]`" is right. In fact the termwise antiderivative of a
  Laurent-log polynomial never produces a `p^0` term (`Int dp/p = log p`; `Int p^a log p` has no `p^0` term for `a != -1`). So
  "p-constant dropped" simply means "integration constant 0".
* The D(ii) formula `log coefficient of Phi_k^b = c_b Sum_i Pi_i(tau_k H^b)(-x)^i/i!` with `c = (-1, 1, -1, -3, -2)` is right
  for every k (s.3).
* **Correction 1 (the `c_b` are not cosmetic).** The sentence "Rescaling by the `c_b` changes no statement below" is true for
  statements about the SET of legs `Lambda(K_0 cap S_w)`. It is false for the identification of `K_0` inside `S_w`.
  * Example: `e = tau_1(1) + v tau_0(H)` lies in `K_0`, since `Pi_0(e) = -v + v = 0` and `Pi_i(e) = 0` for `i >= 1`.
  * Its raw leg `Phi_1^0 + v Phi_0^1 = -p^5/5 + 2 v log p` is NOT log-free, whereas `Lambda(e) = p^5/5` is.
  * With raw legs, the log-free kernel differs from `K_0` (defined by the `Pi_i`) at 390 of the 401 weights `w <= 400`
    (`rf_dii.py` control).
  * `../dict_verify.py` T7 and the critic's Test B both compute the raw-leg kernel, i.e. the image of `K_0` under the invertible
    diagonal map `a_t -> a_t/c_b`. The set of legs tested is therefore the same, and no numerical claim is affected. But the
    "16200 K_0 basis vectors" are vectors of the rescaled kernel, and the sentence should say so.

**s.2, Lemma 0.** PROVED, correct. The b-weights `1, 2, 3, 4` (for `b = 0, 4, 3, 2`) meet each class mod 4 once. The class-0
candidate fails exactly when `w - 5m = 0`.

**s.3, Lemma 1.** PROVED, correct.
* (a): the injectivity of `Q[x,T] -> Q[p,x]` needs `p, x` algebraically independent. This holds because `v = px - p^5/4` and `p, v`
  are independent. It is implicit in the text; one line suffices.
* (b) and (c) are right as written. In (b) the top-coefficient induction works for arbitrary `a_n` in `Q(v)`.
* (d) and the byproduct `floor(w/4) + 1 = d(w) + e(w)` are right.

**s.4, Lemma 2.** PROVED, correct.
* (i): I re-derived the formula, including both exceptional cases.
* (ii): the termwise antiderivative and the formula `-log p . G + Int G dp/p` agree exactly. Write `G(x(p)) = Sum_j g_j p^j`. Then
  both sides equal `log p (G - g_0) - Sum_(j != 0) g_j p^j/j`.
* The remark "no p-constant appears because `G(0) = 0`" is not needed: the identity holds for every primitive. It is harmless.
  Note that `g_0 != 0` does occur when 5 divides the x-degree. For example, `G = x^5/5` has `g_0 = v^4/4`, which produces the
  exceptional term `(v^4/4) log p`.
* (iii) is right. `M` is moreover a `Q[v,x]`-module, which s.3 below uses.

**s.5, Lemma 3 and the Theorem.** PROVED, correct.
* Lemma 3: the source and target bases mod `v`, the relation `q = -p^4/4` mod `Q[v,x] + pq Q[p,q]`, the closed images and the
  triangularity are all right. `v` is a non-zero-divisor on `N`, so `dim (N/vN)_w = dim N_w - dim N_(w-5)`, consistent with the
  stated bases.
* (d): right exactness, graded Nakayama (`C_w = 0` for `w < 0`, `v` of positive weight) and the dimension count are right.
* (b) is right. In (c), injectivity, surjectivity and the use of D(ii) in the direction "log-free implies all `Pi_i = 0`" are right.
* **Correction 2 (citation).** The last line of the proof of (c), and s.0, say "`dim K_0 cap S_w = d(w)` by Lemma 1(d) and (a)".
  But (a)/Lemma 0 count triples. `dim K_0 cap S_w = d(w)` follows from (d) (or Cor. 4) together with D(ii).
* **Correction 3.** s.6 and the JSON say "(d) and Cor. 4 do not use D(ii)". As stated, with `K_0`, Cor. 4 does use D(ii): it needs
  `ker(log map) = K_0`. Cor. 4 is D(ii)-free with "log-free combinations" in place of `K_0`.
* Cor. 5: PROVED given D(ii). It agrees with the model note's list of basis insertions in `K`. For instance, at `w = 5`,
  `K cap S_5 = S_5` contains `tau_1(1)`; at `w = 9`, `K_0 cap S_9 = <v tau_0(H^2), tau_2(1) + v tau_1(H)>`.

**s.6 (what is used).** **Correction 4.** The list should add two Stage-1 inputs.
* (i) The leg dictionary (SYNTHESIS s.6 item 1): the TR legs of the model's ancestor insertions are `Phi_k^b` up to `c_b` and
  `Q[v][x]`. The theorem is sensitive to the primary legs:
  * Replacing `Phi_0^2 = q` by `x` kills injectivity at `w = 4` (control K1).
  * More precisely, for `Phi_0^2 = a p^4 + b q` Lemma 3 predicts that (d) holds iff `a != b/4`, i.e. iff `Phi_0^2` is not a
    multiple of `x`. `rf_h2leg.py` confirms this for 7 choices of `(a, b)`, `w <= 40`.
  * By weight, `Phi_0^2` is the only primary leg not fixed up to scale: `M_4 = <p^4, q> (+) log p . <x>`. In weights 0-3 the legs
    `log p`, `p`, `p^2`, `p^3` are forced up to scale and constants.
* (ii) The closed form of `PiC` (s.3).

**s.9 (consequences).** **Correction 5.** "With D(ii), this closes ML_K restricted to `K_0` insertions modulo E9" needs three more
hypotheses:
* (i) the Stage-1 identification of the model's ancestor correlators with the TR residue pairings against these legs
  (arXiv:1509.06954 Thm 8.1, together with the dualised descendant functions (2-35));
* (ii) the vanishing of pairings with functions of `x` (standard, from the linear loop equations);
* (iii) the model note's statement that the generators `kappa(X, k) = tau_k(X) + (lower k)` of `K_0` have coefficients in `Q[v]`.
  Given that, D(iii) itself shows that the `v^m kappa(X, k)` of weight `w` form a basis of `K_0 cap S_w`: they are independent,
  lie in `K_0 cap S_w`, and number `d(w)` (one for each triple with `b != 1`). The coefficient ring matters, because the theorem
  is false over `Q(v)`: control K7 finds 42 of 222 log-free combinations outside `L` once `m = -1` is allowed.

With these listed, the s.9 claim is fine. I confirmed that arXiv:2212.00320 is the x-y swap paper, whose abstract announces explicit
formulas for unramified `y` and arbitrary rational `x`. I did not re-read its Theorem 7.1. That theorem is not used by D(iii) itself,
and the Stage-1 xyswap referee checked its hypotheses.

## 3. D(ii) in the form used: a complete proof (PROVED here)

**Lemma G (generating function).** For all `b, k`: the log coefficient of `Phi_k^b` is `Sum_(i=0..k) R^b_(k-i) (-x)^i/i!`, where
`R^b_n := Res_(p=0) W^n dPhi_0^b / n!`, a monomial in `v`.

*Proof.*
* Put `B_n := Sum_(j=0..n) x^j Phi_(n-j) / j!`, the coefficient `[z^n]` of `e^(zx) Sum_k z^k Phi_k`.
* Since `dPhi_m/dp = -x' Phi_(m-1)` for `m >= 1`, the sum telescopes: `dB_n/dp = x^n Phi_0'/n!`, a Laurent polynomial.
* `B_n` lies in `M`, because `M` is a `Q[v,x]`-module (Lemma 2). Write `B_n = A + log p . B` with `B` in `Q[v,x]`.
* Then `dB_n/dp = A' + B/p + log p . B'` has no log term. So `B' = 0`: `B` is a `p`-constant, and it is the `p^-1` coefficient of
  `x^n Phi_0'/n!`, i.e. `B = R_n`. (`A'` has no `p^-1` term.)
* Inverting, `Phi_k = Sum_j (-x)^j/j! B_(k-j)`. Take log coefficients. □

**Residue identity.**
* Lemma V of `gauge/GAUGE_IDENTITY.md` gives `I^C_m = Res_(p=0) p^m e^(-W/z) dp`. Hence `[z^-n] I^C_m = (-1)^n Res p^m W^n dp/n!`.
* With `Y^C = (I_0, -I_-1, -I_3, I_2, I_1)` and `PiC(tau_n H^b) = (-1)^(n+1) [z^-n] Y^C_b`, this gives
  `PiC(tau_n H^b) = -s_b Res p^(m_b) W^n dp / n!`, with `(m_b, s_b) = (0,1), (-1,-1), (3,-1), (2,1), (1,1)`.
* Compare with `dPhi_0^b = dp, dp/p, -v p^-2 dp, 3p^2 dp, 2p dp`. For `b = 2` use
  `Res W^n (p^3 - v p^-2) dp = Res d(W^(n+1))/(n+1) = 0`.
* Result: `R^b_n = c_b PiC(tau_n(H^b))` with `c = (-1, 1, -1, -3, -2)`, for EVERY `n`.
* Together with `Pi_i(tau_k H^b) = PiC(tau_(k-i) H^b)` this is exactly the D(ii) statement of DICT_PROOF s.1. □

**Status and checks.**
* PROVED, given Lemma V (PROVED in the gauge note) and the Stage-1 definition `Pi_0 = PiC` through `Y^C`.
  * The identification of `Y^C` with `-(1/5) eta(S(z)H^2, .)` is Stage 1 (model note s.1, `ms_model.py`).
  * Its sign conventions were cross-checked in SYNTHESIS s.3 C4 (`cr_conv` (2), `m <= 60`).
* VERIFIED by `rf_dii.py`:
  * D1, Lemma G on the iterated legs: `k <= 100`, over Q and mod both primes.
  * D2, the residue identity: W^n by repeated products for `n <= 400`; binomial form for `n <= 2000` over Q and `n <= 3000` mod
    both primes; the 7 `PiC` values listed in the model note.
  * D3: the kernel of the closed-form `Pi_i` equals the `c_b`-rescaled log-free kernel of the legs, with dimension `d(w)`, for
    `w <= 400` over Q and `w <= 800` mod both primes.
  * Controls, all failing as they should: `c_3 -> +3` (480 mismatches); sign `(-1)^n` (2402); `Pi_i` shifted by one (80 weights);
    no `c_b` (390 of 401 weights).

## 4. Re-derivation of the VERIFIED claims (own code; logs in this folder)

| worker claim | my test | range | result |
|---|---|---|---|
| T1: Lemma 0 and Lemma 1(d) counts, `w <= 400` | R0: brute-force counts of `S_w`, `Q[v,x]_w`, `Q[p,q]_w`, `Q[v][p]/Q[v]` | `w <= 4000` (`S_w` listed for `w <= 1000`) | pass |
| T2: Lemma 1(c) rule; 861 in, 859 out, `w <= 80` | R1: my elimination L-test against the support rule `j >= -m`, on every monomial with `m <= w/5 + w/4 + 4`, plus 606 random elements (half with one bad monomial) | `w <= 100`; Q, 2^61-1, 2^89-1 | 1326 in (all `j >= -m`), 1374 out (all `j < -m`); 606/606 agree |
| T3: Lemma 2 (i), (ii); 1042 elements, `w <= 80` | R2: closed formulas (i) and (ii), and image in `M`, on the monomial basis of `M_w` | `w <= 100` | 1602 basis elements, pass |
| T4: every leg in `M`, `k <= 101` | R3: `in_M` for each leg; the reduction of Lemma 1(a) never needs a `Q(v)[x]` multiplier | `k <= 200` over Q; `k <= 1000` mod both primes | pass |
| X1: legs equal the critic's, `k <= 100` | R3: iterated legs against the closed generating-function legs, a different algorithm; the difference must lie in `Q[v][x]` | `k <= 100` over Q; `k <= 120` mod both primes | pass |
| X2 and T5: Lemma 3 mod-`v` images | R4: closed images mod `v` and the 2x2 determinants in weights `4n` | `k <= 200` over Q; `k <= 1000` mod both primes | pass |
| Lemma 3, the basis of `N/vN` | `rf_nvn.py`: quotient dimension `[w >= 1]`, `p^w` not in the subspace, `q^n = -p^(4n)/4^n`; wrong-sign control | `w <= 320`, over Q | pass; control fails at all 80 weights |
| T6: Theorem (d), rank = `|S_w|`, `w <= 400` | R5: generic elimination of the full matrix `S_w -> N_w` | `w <= 400` over Q; `w <= 300` mod both primes | pass |
| T7: Theorem (b), (c); 16200 vectors, `w <= 400` | R6: log-free kernel by generic elimination. Each leg is log-free, in `Q[p,q]` (support rule) and in `L` (elimination test); the normal form is `Q[v]`-linear; rank mod `Q(v)[x]` = `d(w)` | `w <= 400` over Q (16200 vectors); `w <= 600` (36300) and `w <= 1000` (100500) mod both primes | pass |
| X4: Cor. 4; Cor. 5 | R7: the `b = 1` log rows are triangular with diagonal `(-1)^k/k!` and rank `e(w)`; `dim K cap S_w = d(w) + [5|w]` | as R6 | pass |
| critic Test A, D(ii), `k <= 100` | `rf_dii.py` D1-D3 | see s.3 | pass |
| controls C1-C4 | `rf_controls.py` K1-K7, the 4 `rf_dii.py` controls, `rf_crctrl.py`, the `rf_nvn.py` control | `w <= 40` / `w <= 400` | all 13 fail as they should |

Notes on the table:
* `rf_counts.py` reproduces the worker's quoted counts exactly with my code: 861/859, 1042, 231, and 2640/10400/16200.
* Mod-prime runs use 2^61-1 and 2^89-1. A rational coefficient at a forbidden key could vanish modulo one prime; the two primes and
  the exact range `w <= 400` rule this out in practice. A rank drop mod `p` would show up as a dimension mismatch; none occurred.
* The controls are listed in `rf_controls.py`. K1 (`Phi_0^2 := x`) first fails at `w = 4`. K2 (`v^2/p^6`) and K3 (`p + v/p^4`)
  each push 13 legs out of `M`, and 31 resp. 61 of 180 log-free combinations out of `Q[p,q]` and out of `L`. K4 (a log-dropping
  integration) fails from `w = 0`. K5 checks 8 anchors of the L-test. K6: 156/156 perturbations `v^m p^j` with `j < -m` are
  rejected and 170/170 with `j >= -m` are accepted. K7 is the `Q(v)` case (s.2, Correction 5).

## 5. Attempts to break the falsifier

* **The logs are genuine.** `critic/cr_dictionary.py 80 320`, re-run from my own unchanged copy (12 min 56 s), is byte-identical to
  `../cr_dictionary_80_320.log`. The worker's scratch copies have the same SHA-256 as `critic/`. I did not re-run `100 400`: it is
  the same code (77 min), and my independent code covers `w <= 400` exactly.
* **The test can fail.** `rf_crctrl.py` imports the critic's script read-only and replaces the H^2 leg by `v^2/p^6` in memory.
  At `KMAX = 12`, `WMAX = 48` it prints 40 FAIL lines: 39 "log-free combination NOT in L" and one Test A ratio. `critic/` is
  unchanged afterwards.
* **The critic's code.** I read `cr_laurent.l_reduce` and `cr_dictionary.rank_kernel` line by line.
  * `l_reduce` is a correct L-test for `w >= 0`: once the negative powers of `p` are removed by `Q(v)[x]`-multiples, the element is in
    `L` iff the remainder lies in `Q[v][p]` up to a constant.
  * `rank_kernel` returns a correct left kernel: the unused rows are zero after the elimination.
  * `KMAX = WMAX/4` reaches every triple, since `b = 1` needs `k <= w/4`.
* **Larger range, independent code.** No log-free leg outside `L` in any weight `w <= 1000`, modulo both primes. The dimensions and
  the bijectivity also hold in every weight.
* **Why it cannot fail.** Every check of the falsifier is implied by the PROVED Lemmas 1-3. The "in `L`" part needs only Lemma 2,
  a closure property of single legs, which is VERIFIED per leg for `k <= 1000`.
* **Weak points probed.** None breaks the theorem.
  * D(ii) beyond Test A's range: closed for all `k` by s.3.
  * Coordinates of `K_0` against leg sets: Correction 1.
  * The `Q[v]` scope: K7.
  * The primary H^2 leg: K1.
  * Normalisation constants of integration: they change legs by `Q[v][x]`, which lies in `Q[p,q]`. Membership "exactly in `Q[p,q]`"
    therefore holds for every homogeneous convention.

## 6. Summary

**Errors in the proof: none. Gaps: none in D(iii) itself.** Corrections and imprecisions (none changes a theorem statement):
1. s.1: "Rescaling by the `c_b` changes no statement below" holds for leg sets only, not for `K_0` as a subspace of `S_w` (example
   `e = tau_1(1) + v tau_0(H)`). The "16200 `K_0` vectors" of T7 are raw-leg kernel vectors.
2. Proof of (c), and s.0: `dim K_0 cap S_w = d(w)` comes from (d) (or Cor. 4) plus D(ii), not from (a).
3. s.6 and the JSON: Cor. 4 uses D(ii) as stated with `K_0`; it is D(ii)-free only for "log-free combinations".
4. s.6: list the Stage-1 inputs, namely the leg dictionary (with `Phi_0^2 = v/p`, to which the theorem is sensitive) and the closed
   form of `PiC`. Cite D(ii) through s.3 of this report rather than the "LG dictionary" sentence of SYNTHESIS s.6.
5. s.9: "closes ML_K on `K_0` modulo E9" also needs the TR identification of the model's ancestor correlators (arXiv:1509.06954
   Thm 8.1 and the dualised (2-35)), the vanishing of pairings with functions of `x`, and `kappa(X, k)` having `Q[v]` coefficients.
   Given the last, the `v^m kappa(X, k)` form a basis of `K_0 cap S_w`, by D(iii)'s count.
6. Lemma 1(a) uses the algebraic independence of `p` and `x`, which is not stated. Lemma 2(ii)'s "because `G(0) = 0`" is unnecessary.

**What survives (statuses).**
* Theorem D(iii) (a)-(d) for the explicit legs `Phi_k^b`: PROVED. (a) and (d) need nothing; (b) and (c) need D(ii). Reading the
  `Phi_k^b` as the TR legs of the ancestor insertions is the Stage-1 leg dictionary.
* D(ii) in the used form, for all `k`: PROVED, given Lemma V (gauge note) and `Pi_0 = PiC` via `Y^C` (s.3).
* Cor. 4 and Cor. 5: PROVED, given D(ii).
* Dimension claim `d(w) = floor(w/5) + 1 - [5|w]` on both sides: PROVED. VERIFIED for `w <= 4000` (counts) and in every run.
* Falsifier (a): passed. VERIFIED: exact over Q for `w <= 400`; mod 2^61-1 and 2^89-1 for `w <= 1000`; critic re-run byte-identical.
* s.9 consequence: PROVED modulo E9 and the inputs of item 5.

**Suggested status line for the Stage-2 record:**
"D(iii): PROVED for all weights (`stage2/dict/DICT_PROOF.md`; D(ii) in the used form, referee s.3), given the Stage-1 leg
dictionary and `Pi_0 = PiC`; VERIFIED for `w <= 400` (Q) and `w <= 1000` (two primes)."

**Process notes.**
* No existing file was modified. Runs used my scratch folder, at most 4 concurrent processes; none is left running.
* My first mod-prime runs printed spurious FAILs, caused by an unreduced coefficient in my own `x_coords`. I stopped them, fixed the
  bug and re-ran; the mod-prime logs here are from the fixed code.
  * Over Q the fix is a no-op, because Fractions are canonical. The exact `w <= 400` run had started before the fix.
  * The `w <= 160` run, repeated with the final code, gives identical output.
* My first `rf_dii` D3 omitted the `c_b` and failed. Correcting it is what exposed Correction 1.
