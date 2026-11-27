# Falsifier (b) of the Stage-2 route: TR correlators paired with all L-legs and the log leg (2026-11-14)

Folder `allgenus/stage2/logleg/`. AI-written; no human has read it. Labels as in `LOGLEG_PROOF.md`. Task: SYNTHESIS_STAGE01 s.6,
cheapest falsifier (b): compute `omega_(4,2)` and `omega_(3,3)` on the model curve `x = W = p^4/4 + v/p`, `y = p`, pair them with every
leg of the class `L` and with `-log p` in every weight that can be singular, and test the ML_K prediction. Done first, before any proof.

## 0. Verdict

**ML_K SURVIVES. No pairing with an L-leg and negative v-exponent is non-zero; every pure log pairing equals the Gap-Theorem
derivative.** VERIFIED mod 2^89-1 and 2^61-1 (identical supports, every non-zero value rationally reconstructed at 2^89-1 and
confirmed at 2^61-1), for `(h,n) = (4,2), (3,3)` (the target, including the model worker's missing amp(4,2)), `(2,4), (1,5), (0,6)`
(free by-products of the recursion), and every `(h,n)` with `2h - 2 + n <= 6` below them; exactly over `Q(zeta5)` for `(0,3), (0,4),
(0,5), (1,1), (1,2), (1,3), (2,1), (2,2), (3,1)`. Controls behave: off-L single-valued legs give non-zero singular values, functions
of x pair to exactly zero. Since the L-legs `p^m` span `L` mod `Q[v][x]` and a correlator is `const * v^e`, this is an EXHAUSTIVE test
of ML_K at every tested `(h,n)`, given D(iii) in weights `<= 40` (VERIFIED by the critic to weight 240).

Cost: `omega_(4,2)` together with every lower `omega` it needs (`(3,3), (2,4), (1,5), (0,6)`, ...) takes 52-61 s per prime with the
allgenus/tr engine (`tr_rec.TR`, imported read-only); the largest tensors are `omega_(2,4)` (134120 e-basis coefficients) and
`omega_(1,5)` (127080); `omega_(4,2)` has 4320, `omega_(3,3)` 44160. Pairing all singular weights of `(1,5)` takes 50 s, of `(4,2)` 0.3 s.
Wall time of the whole two-prime run for the five new cases: 448 s (one process). The model engine's estimate for amp(4,2) was about
4 CPU hours.

## 1. What is computed

* `omega_(h,n)` at `v = 1` in the basis `e_{b,k} = dp/(p - c_b)^(k+1)` (tr_rec, residues at `c_0` + Z/5 rotation), mod a prime.
* L-leg `p^m`: `Sum_b Res_{c_b} p^m e_{b,k} = binom(m,k) c_b^(m-k)` (= coefficient at infinity; `tr_v0.at_infinity`).
  Log leg `-log p`: `Int_0^oo e_{b,k} = (-1)^k/(k c_b^k)` (the ancestor `tau_0(H)`, LOGLEG_PROOF s.1).
* For every `(h,n)` and every split `n = a + k` (`a` L-legs, `k` log legs): all ordered tuples `(m_1..m_a)`, `m_i >= 1`, with
  `Sum m_i < 5(2h-2+n)` (negative exponent `e = 2 - 2h - n + Sum m/5`), plus the boundary `Sum m = 5(2h-2+n)` (e = 0) for the record.
* Predictions of ML_K: (P1) `a >= 1`, `e < 0` => 0; (P2) `a = 0` => `<(-log p)^k>_h = F_h^(k)(1)`, `F_h = B_2h/(2h(2h-2)) v^(2-2h)`
  (TR normalisation of x = W; frozen = `(-5)^(1-h)` x TR), `F_1 = -(1/12) log v`, `F_0 = (1/2) v^2 log v`; (P3) Z/5: `Sum m = 0 mod 5`.

## 2. Table (per prime; ordered weight tuples; "sing" = number of tuples with e < 0 tested, all found ZERO)

| (h,n) | `a=n,k=0` | `a=n-1,k=1` | `a=n-2,k=2` | `a=n-3,k=3` | `a=n-4,k=4` | `a=n-5,k=5` | `<(-log p)^n>_h` = `F_h^(n)(1)` | frozen value |
|---|---|---|---|---|---|---|---|---|
| **(4,2)** | 133 | 7 | - | - | - | - | -7/240 | 7/30000 (new) |
| **(3,3)** | 986 | 99 | 6 | - | - | - | -5/42 | -1/210 |
| (2,4) | 3445 | 580 | 70 | 5 | - | - | -1/2 | 1/10 |
| (1,5) | 5004 | 1421 | 304 | 46 | 4 | - | -2 | -2 |
| (0,6) | 2128 | 1128 | 452 | 133 | 27 | 3 | -6 | 30 |
| (4,1) | 6 | - | | | | | 1/240 | -1/30000 |
| (3,2) | 70 | 5 | | | | | 5/252 | 1/1260 |
| (2,3) | 304 | 46 | 4 | | | | 1/10 | -1/50 |
| (1,4) | 452 | 133 | 27 | 3 | | | 1/2 | 1/2 |
| (0,5) | 127 | 88 | 42 | 13 | 2 | | 2 | -10 |
| (3,1) | 4 | | | | | | -1/252 | -1/6300 |
| (2,2) | 27 | 3 | | | | | -1/40 | 1/200 |
| (1,3) | 42 | 13 | 2 | | | | -1/6 | -1/6 |
| (0,4) | 4 | 6 | 4 | 1 | | | -1 | 5 |
| (2,1) | 2 | | | | | | 1/120 | -1/600 |
| (1,2) | 4 | 1 | | | | | 1/12 | 1/12 |
| (0,3), (1,1) | 0 | 0 | 0 | | | | 1, -1/12 | -5, -1/12 |

Totals per prime: 17416 singular ordered tuples (4678 of them mixed, `a >= 1` and `k >= 1`), all zero; 0 Z/5 selection violations;
18 pure-log values, all equal to the Gap-Theorem derivative. The sharpest new entries: `<p^m, -log p>_4 = 0` for `m = 5, .., 35`
(the seven mixed genus-4 pairings, exponents -7..-1), `<p^m1, p^m2>_4 = 0` for all 133 ordered pairs with `m1 + m2 <= 35`,
`<p^m1, p^m2, -log p>_3 = 0` (99 pairs) and `<p^m, (-log p)^2>_3 = 0` (m = 5..30).
Boundary (e = 0, regular, recorded): non-zero in every case, e.g. `<p^40, -log p>_4 = -868220928`, `<p^35, (-log p)^2>_3 =
-144609822225/1024`, and `<p^(5(2h-1))>_h` = 2079/128, 3132675/1024, -45541532421/32768 (h = 2, 3, 4) = the v -> 0 limits of
`omega_{h,1}` found by the tr worker (TR_MODEL s.5): an independent cross-check of engine and pairing.

## 3. Controls (`ll_falsifier_*.log`, every (h,n))

* Functions of x (`x`, `x^2` as legs, with the other slots log legs, and `x` with all `p^m`): exactly zero in every case (PASS).
* Single-valued legs outside `L`: `p^-5`, `p^-1`, `p^-2` (n >= 2, with one `p^m` and log legs) and `p^-5`, `p^-10`, `p^-15` (n = 1):
  non-zero values with negative exponent in every `(h,n)`, e.g. `<p^-5, p^25>_4 = -518283675/16`, `<p^-2, p^37>_4 = -2084566248765/1024`,
  `<p^-5>_1 = 15/4`, `<p^-1, p^1, -log p>_0 = -1`. (Some single controls vanish, e.g. `<p^-5>_3 = <p^-5>_4 = 0` and the `p^-1` control
  at (1,2), (4,2); the check requires at least one non-zero singular value per `(h,n)`, which holds. A first run with only `p^-5` at
  n = 1 flagged (3,1), (4,1) as "control failed"; the controls were then extended to `p^-10`, `p^-15`, not the predictions.)
* For the proof steps (LOGLEG_PROOF_2 s.4): the variational formula without, or with the wrong-sign, Lie term fails; the step (V1)
  with the xi-term dropped or sign-flipped fails.

## 4. Comparison with the model worker (`model/MODEL_SINGULAR*.md`, `ms_ktest*.log`)

* (ML1) values: all 17 of the model's `<H^n>_h` (ms_ktest.log, ms_ktest_24.log) are reproduced as `(-5)^(h-1)`-rescaled pure log
  pairings: `-5, 5, -10, 30` (h = 0, n = 3..6), `-1/12, 1/12, -1/6, 1/2, -2` (h = 1), `-1/600, 1/200, -1/50, 1/10` (h = 2),
  `-1/6300, 1/1260, -1/210` (h = 3), `-1/30000` (h = 4); new: `<H^2>_4 = 7/30000 = F_4''` (frozen `F_4 = v^-6/180000`).
* (ML2): the model found 0 non-zero singular K-correlators among 1395 (genus 0-4, incl. (3,3): 246, (2,4): 376, (1,5): 313,
  (0,6): 152). Here the same (h,n) are tested in the leg basis and exhaustively (all weights), and (4,2), which the model could not
  reach, is added. The two bases are related by D(iii) (`K_0` <-> `L` mod functions of x), so the data agree; a value-by-value
  comparison would need the explicit leg of each `kappa(X,k)` and was not made.
* Normalisation: TR pairing = `Lambda^(h-1)` x frozen ancestor correlator, `Lambda = -5` (TR_MODEL s.2.4), confirmed by all 18 values.

## 5. Status and caveats

* Finite checks are evidence, not proof. Labels: (P1)-(P3) VERIFIED in the stated range, two primes (exact Q(zeta5) where stated).
* The translation of the TR test into ML_K uses DNOPS Thm 8.1 (quoted) and D(iii) (VERIFIED, weights <= 240).
* Not done: `omega_(4,3)` (needs `omega_(3,4), (2,5), (1,6), (0,7)`; `omega_(1,6)` alone was computed in the `ll_eo` runs, about 9
  minutes per prime with everything else); not needed for the verdict.
* Scripts: `ll_core.py` (pairings), `ll_falsifier.py`; logs `ll_falsifier_new.log` ((4,2), (3,3), (2,4), (1,5), (0,6), both primes),
  `ll_falsifier_low.log` (13 lower cases, both primes), `ll_falsifier_exact.log` (Q(zeta5)).
