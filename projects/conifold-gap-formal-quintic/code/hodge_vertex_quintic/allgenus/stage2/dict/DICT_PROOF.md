# D(iii) for all weights: `K_0` = the TR leg class `L` modulo functions of x (Stage 2, first lemma; 2026-11-14)

AI-written; no human has checked it. Labels as in `../../MECHANISM.md`: PROVED, VERIFIED, CONJECTURAL, SKETCH, REFUTED.
Conventions of `../../SYNTHESIS_STAGE01.md` s.3 (`H = d/dv`, `K_vc` = the slice `K`, `K_0 = K cap ker Pi_0`, "ML_K", "quintic Theorem R").

## 0. Bottom line

* **Theorem D(iii): PROVED for every weight `w >= 0`** (s.1-s.5). The proof is elementary: integration by parts, a normal form, and
  graded Nakayama over `Q[v]`. It uses D(ii) (PROVED, SYNTHESIS s.6) only to identify `K_0` with the log-free combinations. It uses
  no TR theorem and no quoted literature.
* **Sharper than asked.**
  * `L` restricted to legs with `Q[v]`-coefficients is the ring `Q[p,q]`, `q = v/p`: a leg lies in `L` iff every term `v^m p^j` has `j >= -m`.
  * Every leg, log or not, lies in `M = Q[p,q] (+) log p . Q[v,x]`.
  * The full leg map `S_w -> M_w/Q[v,x]_w` is bijective. So `N = M/Q[v,x]` is a free `Q[v]`-module with basis the legs of the `tau_k(H^b)`.
  * `S_w = K_0 (+) span{v^m tau_k(H)}` (Cor. 4).
* **Falsifier passed.**
  * `critic/cr_dictionary.py 80 320` (scratch copy): ALL CHECKS PASSED, 10400 log-free basis vectors, all in `L` (13 min).
  * `critic/cr_dictionary.py 100 400`: ALL CHECKS PASSED, 16200 log-free basis vectors, all in `L`; D(ii) Test A on 101
    non-zero and 404 zero slots (77 min).
  * No log-free leg outside `L` was found.
* **VERIFIED step by step** by new, independent exact code over `Q` (s.7), with failing controls:
  * weights `<= 400` for the leg-map statements;
  * weights `<= 80` for the monomial-level lemmas.
* The dimension claim `dim K_0 = dim (L/Q(v)[x])_w = floor(w/5) + 1 - [5|w]` was checked first. It is PROVED (Lemma 0, Lemma 1(d)) and VERIFIED for `w <= 400`.

## 1. Set-up and statement

* **Curve, weights.** The model curve is `x = W = p^4/4 + v/p`, `y = p`. Weights: `p` has 1, `v` has 5, `x` has 4, `q := v/p` has 4, `log p` has 0.
  * `dx = (p^3 - v p^-2) dp`.
  * `p` is a root of `T^5 - 4xT + 4v`.
* **Legs.** These are the legs of `critic/cr_dictionary.py`; `dict_xcheck.py` X1 confirms they are equal term by term.
  * `Phi_0^b = p, log p, v/p, p^3, p^2` for `b = 0..4`.
  * `Phi_(k+1)^b = I(Phi_k^b)`, where `I(F) := -Int F dx`, integrated termwise in `p` with the `p`-constant dropped.
  * A different constant convention changes legs only by elements of `Q(v)[x]`: `I` maps `Q(v)[x]` into itself up to constants (Lemma 2(iii)). So every statement modulo `Q(v)[x]` is convention free.
* **Weights of legs, and `S_w`.** `wt(Phi_k^b) = w_b + 4k`, with `(w_0..w_4) = (1, 0, 4, 3, 2)`. Then `S_w` is the `Q`-span of the triples
  `v^m tau_k(H^b)`, `m in Z_(>=0)`, with `w_b + 4k + 5m = w`. `F = (+)_w S_w` is the free graded `Q[v]`-module on the `tau_k(H^b)`.
* **D(ii)** (PROVED, SYNTHESIS s.6; VERIFIED by critic Test A for `k <= 60`, re-run here to `k <= 80` and `100`).
  * The log coefficient of `Phi_k^b` is `c_b Sum_i Pi_i(tau_k(H^b)) (-x)^i/i!`, with `c = (-1, 1, -1, -3, -2)`.
  * This is `cr_dictionary` Test A: `mu_0 = c_b PiC`, together with `Pi_i(tau_(k+1)) = Pi_(i-1)(tau_k)`.
  * Put `Lambda(v^m tau_k(H^b)) := c_b^-1 v^m Phi_k^b` (`Q`-linear).
  * For `s in S_w`: **`Lambda(s)` is log-free iff `Pi_i(s) = 0` for all `i >= 0`, i.e. iff `s in K_0`**, since the `x^i` are independent.
  * Rescaling by the `c_b` changes no statement below.
* **The class `L` and the target space.** Following the critic, `L := Q[v][p] + Q(v)[x]`. xyswap's class is `Q[v][p] + Q[v][x]`, and the
  two agree on legs (Lemma 1). Write `d(w) := floor(w/5) + 1 - [5|w]` and `e(w) := #{(m,n) in Z_(>=0)^2 : 5m + 4n = w} = dim Q[v,x]_w`.

**Theorem D(iii).** For every `w >= 0`:
* (a) `S_w` has `d(w)` triples with `b != 1` and `e(w)` triples with `b = 1`.
* (b) For `s in K_0 cap S_w`, the leg `Lambda(s)` lies in `Q[p,q] = Q[v][p] + Q[v][x]`, inside `L`. This holds exactly, not only modulo functions of x.
* (c) `Lambda : K_0 cap S_w -> (L/Q(v)[x])_w` is a bijection, and both sides have dimension `d(w)`.
* (d) `Lambda : S_w -> N_w := (Q[p,q]/Q[v,x])_w (+) log p . Q[v,x]_w` is a bijection. Hence `N` is a free `Q[v]`-module with basis
  `{Lambda(tau_k(H^b))}`.

**The `Q[v]` coefficients are needed.** Over `Q(v)` the statement fails: `tau_0(H^2)/v` has leg `1/p`, which is not in `L` (control C4).

## 2. Lemma 0 (counting): PROVED

For fixed `m` with `5m <= w`, the four weights `{1, 2, 3, 4}` of `b = 0, 4, 3, 2` meet every class mod 4 exactly once.
* So there is exactly one triple with `b != 1` and this `m`, unless `w - 5m = 0`; then the class-0 candidate `b = 2` would need `w - 5m >= 4`.
* Hence `#{b != 1} = (floor(w/5) + 1) - [5|w] = d(w)`.
* The `b = 1` triples are the pairs `(m,k)` with `5m + 4k = w`, so there are `e(w)` of them. VERIFIED by T1 for `w <= 400`.

## 3. Lemma 1 (the class `L`): PROVED

Here `Q[p,q]`, `q = v/p`, is the subring of `Q[v][p, p^-1]` spanned by the `p^a q^c = v^c p^(a-c)`, i.e. by the `v^m p^j` with `j >= -m`.

* (a) `Q[p,q] = Q[p,x] = Q[v][p] + Q[v][x]`, and `Q[p,x]` is free over `Q[v,x]` with basis `1, p, p^2, p^3, p^4`.
* (b) `Q[v][p^(+-1)] cap Q(v)[x] = Q[v][x]`.
* (c) `L cap Q[v][p^(+-1)] = Q[p,q]`. In words: a Laurent polynomial with `Q[v]` coefficients is in `L` iff every term `v^m p^j` has `j >= -m`.
* (d) `L/Q(v)[x] = Q[p,q]/Q[v,x] = Q[v][p]/Q[v] = (+)_(i=1..4) Q[v,x] p^i`. It is free over `Q[v]`, of dimension `d(w)` in weight `w`.

*Proof.*
* (a)
  * `q = x - p^4/4` gives `Q[p,q] = Q[p,x]`, and `v = pq = px - p^5/4`.
  * The map `Q[v,x][T]/(T^5 - 4xT + 4v) -> Q[p,x]`, `T -> p`, is onto. Its source is `Q[x,T]`, since `v = (4xT - T^5)/4` can be eliminated, so the map is an isomorphism. That gives the basis.
  * `Q[v][p] + Q[v][x]` lies in `Q[p,q]`, since `v = pq` and `x = p^4/4 + q`.
  * Conversely, let `F in Q[p,q]` be homogeneous, and let `v^m p^-n` (`n > 0`, `m >= n`) be its most negative `p`-term.
    * Then `F - c v^(m-n) x^n` lies in `Q[p,q]`.
    * `x^n = Sum_t C(n,t) 4^(t-n) v^t p^(4n-5t)` has lowest term `v^n p^-n`, so the most negative exponent strictly rises.
    * By induction `F` reaches `Q[v][p]`, subtracting only elements of `Q[v][x]`.
* (b) Let `P = Sum_(n<=N) a_n x^n` with `a_N != 0` and `P in Q[v][p^(+-1)]`. Only `x^N` contributes to `p^(4N)`, with coefficient `a_N 4^-N`.
  So `a_N in Q[v]`; induct on `N`.
* (c) The inclusion `Q[p,q] in L` follows from (a). Conversely, if `F = R + P` with `R in Q[v][p]` and `P in Q(v)[x]`, then
  `P = F - R in Q[v][p^(+-1)] cap Q(v)[x] = Q[v][x]` by (b), so `F in Q[p,q]`.
* (d)
  * `Q[p,q] -> L/Q(v)[x]` is onto, since `Q[v][p]` lies in `Q[p,q]`. By (b) its kernel is `Q[p,q] cap Q(v)[x] = Q[v,x]`.
  * By (a), `Q[p,q]/Q[v,x] = (+)_(i>=1) Q[v,x] p^i`.
  * `Q[v][p] cap Q[v][x] = Q[v]`, because the lowest term of `x^n` is `v^n p^-n`. Hence `L/Q(v)[x] = Q[v][p]/Q[v]`, of dimension `floor(w/5) + 1 - [5|w] = d(w)`.
  * Byproduct: `floor(w/4) + 1 = d(w) + e(w)`. VERIFIED by T1.

VERIFIED as follows:
* T2: every monomial `v^m p^j` with `w <= 80` and `0 <= m <= w/5 + w/4 + 3` is tested with an exact linear solve over `Q` (`dict_nf.L_test`), which is independent of the `j >= -m` rule. Both halves of (c) hold; the rejected half is control C1.
* X3: the reduction in (a) is checked for every `p^a q^c` with `a + 4c <= 40`.

## 4. Lemma 2 (closure of `M` under the descendant step): PROVED

Let `M := Q[p,q] (+) log p . Q[v,x]`. The sum is direct, since `1` and `log p` are independent over Laurent polynomials.
* (i) `I(p^a q^c) = -p^(a+4) q^c/(a-c+4) + p^a q^(c+1)/(a-c-1)`. Read a term with zero denominator as `-v^c log p`, resp. `+v^(c+1) log p`.
* (ii) For `g in Q[v,x]`, let `G` be its `x`-primitive with `G(0) = 0`. Then `I(log p . g(x)) = -log p . G(x) + Int G(x) dp/p`. Here
  `Int p^a q^c dp/p = p^a q^c/(a-c)` if `a != c`, and `= v^c log p` if `a = c`.
* (iii) `I` is `Q[v]`-linear, `I(M) in M`, and `I(Q[v,x]) in Q[v,x]`.

So every leg `Phi_k^b` lies in `M`, by induction from the five `Phi_0^b in M`.

*Proof.*
* (i) Compute `-Int v^c p^(a-c)(p^3 - v p^-2) dp`.
* (ii) Integrate by parts: `d(log p . G(x)) = G dp/p + log p . g dx`. Then `G(x(p)) in Q[p,q]` and (i)-type integration apply.
* (iii) By (i) and (ii), `I` maps `Q[p,q]` into `M`; and `-log p . G` lies in `log p . Q[v,x]`.
  * Every exceptional term is `v^c log p in log p . Q[v,x]`.
  * For `g in Q[v,x]`: `g(x) x' = d/dp G(x(p))` has no `p^-1` term, so `I(g) = -G(x) + (its p^0-coefficient) in Q[v,x]`.

**Proof of Theorem (b).** For `s in K_0`, `Lambda(s) in M` by Lemma 2, and it is log-free by D(ii). So `Lambda(s) in Q[p,q]`, which lies in
`L` by Lemma 1(a). □

VERIFIED as follows:
* T3: for every basis monomial of `M_w`, `w <= 80`, the generic route (termwise antiderivative of `p^j (log p)^l`, `dict_core.I_A`) lands in `M` and equals formula (i), and the log part matches (ii).
* T4: all `Phi_k^b`, `k <= 101`, lie in `M`.
* X1: the legs equal the critic's recursion (`k <= 100`).

## 5. Lemma 3 (reduction mod v) and the proof of the Theorem: PROVED

**Lemma 3.** The map `Lambda (x) Q[v]/(v) : F/vF -> N/vN` is an isomorphism in every weight.

*Proof.*
* **The source.** In weight `w`, `F/vF` has basis the `tau_k(H^b)` with `w_b + 4k = w`. That is one element if `4 does not divide w`; `{tau_(n-1)(H^2), tau_n(H)}` if `w = 4n >= 4`; and `{tau_0(H)}` if `w = 0`.
* **The target.** `N/vN = Q[p,q]/(Q[v,x] + pq Q[p,q]) (+) log p . Q[x]`.
  * `Q[p,q]/(pq)` has basis `1, p^a, q^c` (`a, c >= 1`).
  * The image of `Q[v,x]` is spanned by `1` and `x^n = p^(4n)/4^n + q^n` mod `pq`.
  * So the first summand has basis `{p^a : a >= 1}`, with `q^n = -4^-n p^(4n)`.
  * Hence `(N/vN)_w` has basis `{p^w}` (`w >= 1`, `4 does not divide w`), `{p^(4n), log p . x^n}` (`w = 4n >= 4`), or `{log p}` (`w = 0`).
* **`I` descends.** `I` is well defined on `N/vN` by Lemma 2(iii).
* **The images** (with `c_b` dropped; all non-zero):
  * By (i), `I(p^a) = -p^(a+4)/(a+4)` mod `vM` for `a >= 1`. So for `j = 1, 2, 3`, `Phi_k = (-1)^k p^(j+4k)/prod_(i=1..k)(j+4i)`.
  * `q = -p^4/4` mod `Q[v,x]`, so `Phi_k^2 = (-1)^(k+1) p^(4k+4)/(4 prod_(i=1..k)(4+4i))`.
  * By (ii), the `log p`-component of `Phi_k^1` mod `v` is `(-1)^k x^k/k!`.
* **Conclusion.** In weight `4n` the matrix in the bases above is triangular, with diagonal entries non-zero. In the other weights it is
  `1 x 1` and non-zero. □

VERIFIED as follows:
* T5: rank of `Lambda` mod `v` = `dim F/vF = [w>=1] + [4|w]` for every `w <= 400`.
* X2: the three closed images, `k <= 100`.
* C2: replacing `Phi_0^2` by the function `x` makes T5 and T6 fail (first at `w = 4`).

**Proof of Theorem (d).**
* `Lambda` is a graded `Q[v]`-linear map between modules with finite-dimensional pieces, both zero in negative weight. Put `C = coker Lambda`. By right exactness `C/vC = coker(Lambda mod v) = 0` (Lemma 3).
* If `C != 0`, let `w` be minimal with `C_w != 0`. Then `(vC)_w = v C_(w-5) = 0`, so `C_w = (C/vC)_w = 0`, a contradiction (graded Nakayama). So `Lambda` is onto.
* `dim S_w = d(w) + e(w)` (Lemma 0), and `dim N_w = d(w) + e(w)` (Lemma 1(d), together with `dim Q[v,x]_w = e(w)`). So `Lambda_w` is bijective, and `N` is isomorphic to `F`, which is free. □

**Proof of Theorem (c).**
* By (b), `Lambda(K_0 cap S_w)` lies in `Q[p,q]`. By Lemma 1(d), `L/Q(v)[x] = Q[p,q]/Q[v,x]`.
* Injective: if `Lambda(s)` lies in `Q(v)[x]`, then `Lambda(s) in Q[p,q] cap Q(v)[x] = Q[v,x]` (Lemma 1(b)). So `Lambda(s) = 0` in `N`, and `s = 0` by (d).
* Onto: given `y in (Q[p,q]/Q[v,x])_w`, put `s := Lambda^-1(y (+) 0)`. Its leg is log-free, so `s in K_0` by D(ii).
* Both sides have dimension `d(w)`, by Lemma 1(d) and (a). □

VERIFIED as follows:
* T6: rank `Lambda = dim S_w`, for every `w <= 400`.
* T7: `K_0` is computed as the exact kernel of the log map. Then:
  * `dim K_0 = d(w)`;
  * each basis leg is in `L` by the independent linear-solve test, and lies in `Q[p,q]`;
  * the rank modulo `Q(v)[x]` is `d(w)` (normal form of Lemma 1(a): coordinates on `p, .., p^4` over `Q[v,x]`).
* C3: the wrong leg `v^2/p^6` (outside `Q[p,q]`) breaks Lemma 2, and some of its single legs fail `L_test`.

**Corollary 4: PROVED.** `S_w = (K_0 cap S_w) (+) span{v^m tau_k(H) : 5m + 4k = w}`.

*Proof.*
* By Lemma 2(ii), the log component of `Phi_k^1` is `(-1)^k x^k/k!` plus terms of lower `x`-degree. Those terms come from `-log p . G` with `x`-degree below `k`, and from the constants `v^c log p`.
* So the log map restricted to the `b = 1` triples is triangular onto `Q[v,x]_w` (`v^m x^k <-> v^m tau_k(H)`), hence bijective. Its kernel `K_0` therefore has this complement.

VERIFIED: X4 (`k <= 100`). A consequence: every `s in S_w` is `s_0 + Sum a_(m,k) v^m tau_k(H)` with `s_0 in K_0`, and the coefficients `a_(m,k)` are read off from the monodromy (the log component).

## 6. Exactly what is used

* **D(ii)**, in the form "`Lambda(s)` is log-free iff `Pi_i(s) = 0` for all `i`", for the legs normalised by `c_b`. This is used only in Theorem (b) and (c), to pass between "log-free" and `K_0`.
  * Everything else is about the leg map alone. (d) and Cor. 4 do not use D(ii).
* **Gradings.** The weights `p:1, v:5, x:4`. Combinations use `m >= 0` (`Q[v]`, not `Q(v)`; this is sharp, C4). The model's `K_0` is weight-homogeneous, being cut out by the homogeneous `Pi_i`; that is standard in Stage 1 and assumed here.
* **Leg conventions.** The legs of `critic/cr_dictionary.py`, up to `c_b` and up to `Q(v)[x]`.
* **Not used.** The generating-function route of SYNTHESIS s.6 item 2 turned out to be unnecessary. Its content (the principal part at `p = 0` of the monodromy-free part) is Lemma 1(c) together with Lemma 2: monodromy-free legs automatically have pole order at most their `v`-degree.
* No TR, no quoted theorem, nothing about the quintic's own mirror curve.

## 7. Verification record (exact over `Q`; `py -I -B`, stdlib)

| script (this folder) | what | range | result |
|---|---|---|---|
| `dict_core.py` | route A Laurent-log calculus, legs, `M`-test (imports nothing from other folders) | | |
| `dict_nf.py` | route B normal form (`p^5 = 4xp - 4v`), exact rank/kernel, `L_test` by linear solve | | |
| `dict_verify.py 400 80` | T1-T7, C1-C3 | `w <= 400` (T1, T5-T7); `w <= 80` (T2, T3); `k <= 101` (T4) | log `dict_verify_400_80.log` |
| `dict_verify.py 160 60` | the same, smaller (12 s) | `w <= 160` | ALL CHECKS PASSED, 2640 `K_0` vectors |
| `dict_xcheck.py 100` | X1-X4, C4; imports `critic/cr_laurent.py` read-only | `k <= 100` | log `dict_xcheck_100.log` |
| `critic/cr_dictionary.py 80 320` (scratch copy) | the critic's falsifier | `k <= 80`, `w <= 320` | log `cr_dictionary_80_320.log` (summary) |
| `critic/cr_dictionary.py 100 400` (scratch copy) | the same, larger | `k <= 100`, `w <= 400` | log `cr_dictionary_100_400.log` (summary) |

## 8. Run results (2026-11-14; all exact over Q, no primes needed)

* `dict_verify.py 400 80`: ALL CHECKS PASSED (13.5 min).
  * T1: `w <= 400`. T2: 861 monomials accepted (all with `j >= -m`) and 859 rejected (all with `j < -m`).
  * T3: 1042 basis elements. T4: `k <= 101`.
  * T5-T7: 16200 `K_0` basis vectors, failing weights none. C2 and C3 fail as they should.
* `dict_verify.py 160 60`: ALL CHECKS PASSED (12 s), 2640 `K_0` vectors.
* `dict_xcheck.py 100`: ALL CHECKS PASSED (46 s). X1 505/505 legs equal the critic's; X2 505 checks; X3 231; X4 101; C4 fails as it should.
* `cr_dictionary.py 80 320` and `100 400` (scratch copies of the critic script, unchanged): ALL CHECKS PASSED (s.0).
* Logs are in this folder: `dict_verify_400_80.log`, `dict_verify_160_60.log`, `dict_xcheck_100.log`, `cr_dictionary_80_320.log`,
  `cr_dictionary_100_400.log`. The two `cr_*` logs have 90 and 94 lines; their Test A lines are long.
* These finite checks are EVIDENCE for the proof's steps. The theorem's status rests on s.2-s.5.

## 9. Consequences for the rest of Stage 2

* ML on `K_0`: every `K_0` insertion has a leg in `Q[p,q] = Q[v][p] + Q[v][x]`, exactly. So E9 (L-regularity, PROVED modulo
  arXiv:2212.00320 Thm 7.1) applies directly; no `Q(v)` coefficients enter. With D(ii), this closes "ML_K restricted to `K_0`
  insertions" modulo E9.
* **Corollary 5: PROVED** (from Cor. 4 and D(ii)). `K cap S_w = (K_0 cap S_w) (+) Q v^(w/5) tau_0(H)` if `5|w`, and
  `K cap S_w = K_0 cap S_w` otherwise.
  * Proof. By D(ii), `s in K` iff the log component `Sum_i Pi_i(s)(-x)^i/i!` of `Lambda(s)` is `x`-free, i.e. equals `Pi_0(s)`, which lies in `Q v^(w/5)`.
  * Write `s = s_0 + t` with `s_0 in K_0` and `t` in the `b = 1` span (Cor. 4). Then `Lambda(t)` has the same log component as `Lambda(s)`.
  * The leg of `v^(w/5) tau_0(H)` is `v^(w/5) log p`, with log component `v^(w/5)`. The log map is injective on the `b = 1` span, so `t = Pi_0(s) v^(-w/5) . v^(w/5) tau_0(H)`.
  * So on `K` the only legs beyond `Q[p,q]` are multiples of `log p`. This is exactly the setting of (V) ("legs `tau_0(H) <-> -log p`
    keep polynomiality"), which stays open (SKETCH, SYNTHESIS s.6 item 3).
