# Gap cluster G3, part 4 of 4: item G3-3, `lem:tra-E` step 2 (O-8); computations; files

Continues `GAPS_G3.md` - `GAPS_G3_3.md`. Vector form (NOTATION_2 s.10 item 6): `f` is a vector, `eta` the Gram matrix, `eta^-1` its
inverse, `E_R` a bivector written as a matrix `E^(ab)` (first index: first slot).

## 1. The ring (answer to O-8)

**The issue.** The source proof (RD2 s.3.1; DC1 B1 step 2) expands `1/(x + psi') = Sum_(j >= 0) (-psi')^j x^(-j-1)` "in each of the two
terms of the numerator", and TRF s.2 calls this "legitimate in `Q((x^-1))[[psi1]]`". Lemma E is stated for every symplectic
`R = 1 + O(z)` in `gl_5(Q[[u, z]])` and applied to `R = P_0^-1 P`, where the dressing `P` is not assumed polynomial in `z` at any order
of `u` (for `P_can` already `P_30(0, z) = (127/3150) z^3 + (1416472237/725359556250000) z^8 - ..`, GI s.5; the computation of s.4 uses an
`R` with all `z`-powers non-zero at `u = 0`). For such `R` the second term `R^-1(x) eta^-1 R^-1(psi')^T x^(-j-1)` has positive powers of
`x` of unbounded degree and is not an element of `Q((x^-1))[[psi']]`; and `alpha(x) = eta f(-x)` has infinitely many negative powers,
so the residue pairing needs a convergence hypothesis.
**The ring.** Let `AA := Q((x))[[u, psi']]`, the power-series ring in `u, psi'` over the field `Q((x))` of Laurent series with finitely
many negative powers. Then:
* `Q[[u, x, psi']]` is a subring of `AA` (it contains `E_R(x, psi')`, `R^-1(x)`, `R^-1(psi')`, the numerator `N(x, psi')` of `E_R`);
* `Q[x^-1][[u]]` is a subring of `AA`; it contains `alpha(x) = eta f(-x)` for every `f` in `HH^u_- = H (x) z^-1 Q[z^-1][[u]]`;
* `x + psi' = x (1 + x^-1 psi')` is a unit of `AA` (`x` is a unit of `Q((x))`, and `1 + x^-1 psi'` has constant term `1` in the
  power-series ring `AA` over `Q((x))`), with inverse `iota := Sum_(j >= 0) (-psi')^j x^(-j-1)`;
* `Res_x : AA -> Q[[u, psi']]`, the coefficient of `x^-1`, taken coefficientwise in `(u, psi')`, is `Q[[u, psi']]`-linear and commutes
  with sums that converge `(u, psi')`-adically in `AA`.
So the sentence of step 2 becomes: *in `AA`, `E_R = iota N = iota eta^-1 - iota R^-1(x) eta^-1 R^-1(psi')^T`, both terms being elements of
`gl_5(AA)`, and the contraction is `Res_x alpha(x)^T E_R(x, psi')` computed in `AA`.* The hypothesis that makes the residue pairing
converge is `f` in `HH^u_-` (u-adic decay of the negative tail); in the applications `f = z^-k J^C` it holds by `lem:tra-conv` (c).

## 2. Lemma E in this ring, with all steps

**Statement (proposed final form of `lem:tra-E`).** Let `R = 1 + O(z)` lie in `gl_5(Q[[u, z]])` with `R(-z)^T eta R(z) = eta`, and let `f`
lie in `HH^u_-` (vector form). Then `E_R` lies in `H (x) H (x) Q[[u, psi, psi']]`, the contraction of `Omega_s(f, .)` into its first slot
(which pairs `psi^m e` with `eta([z^(-m-1)] f(-z), e)`) converges u-adically, and
        `(Omega_s(f, .) (x) 1) E_R(psi') = - R(psi')^-1 (R f)_(>=0)(psi')`   in `H (x) Q[[u, psi']]`,
where `R f` lies in `HH^u` and `(R f)_(>=0)` in `HH^u_+`. So `E_R` annihilates `Omega_s(f, .)` (in either slot) iff `R f` lies in `HH^u_-`.
**Proof.** *Step 0 (the edge is a power series; `def:tra-edge`).* Put `B := Q[[u]]` and `N(psi, psi') := eta^-1 - R^-1(psi) eta^-1 R^-1(psi')^T`
in `gl_5(B[[psi, psi']])` (`R` is invertible in `gl_5(Q[[u, z]])` because `R(u, 0) = 1`). (Sp) inverted gives `R^-1(psi) eta^-1 R^-1(-psi)^T =
eta^-1`, i.e. `N(psi, -psi) = 0`. The substitution `psi' = w - psi` is a `B`-algebra automorphism of `B[[psi, psi']]` (inverse `w = psi +
psi'`); it turns `N` into `N~(psi, w)` with `N~(psi, 0) = 0`, so `N~ = w G` with `G` in `gl_5(B[[psi, w]])`; transforming back,
`N = (psi + psi') E_R` with `E_R` in `gl_5(B[[psi, psi']])`, unique since `psi + psi'` is not a zero divisor. *Symmetry:* `N(psi, psi')^T =
N(psi', psi)` (`eta^-1` is symmetric), hence `E_R(psi, psi')^T = E_R(psi', psi)`.
*Step 1 (the contraction as a residue).* Put `alpha(x) := eta f(-x)`, an element of `x^-1 Q[x^-1][[u]] (x) H` (only negative powers of `x`).
By definition the contraction is `Sum_m eta([z^(-m-1)] f(-z), .) applied to [psi^m] E_R`, and `eta([z^(-m-1)] f(-z), e) = [x^(-m-1)]
alpha(x)^T e`; so it equals `C(psi')^T := Res_x alpha(x)^T E_R(x, psi')`, the product taken in `AA`. Modulo `u^(M+1)` the vector
`alpha` has finitely many non-zero coefficients, so at each order `u^M psi'^n` only finitely many `m` contribute: this is the
u-adic convergence of the contraction.
*Step 2 (the expansion; O-8).* By s.1, `E_R = iota N` in `gl_5(AA)`, because `(x + psi') E_R = N` holds in `Q[[u, x, psi']]`, a subring of
`AA`, and `x + psi'` is a unit of `AA` with inverse `iota`. Distributivity in the commutative ring `AA` gives
`alpha^T E_R = alpha^T eta^-1 iota - alpha^T R^-1(x) eta^-1 R^-1(psi')^T iota`, an identity of elements of `H (x) AA`.
*Step 3 (the first term).* `Res_x(alpha^T eta^-1 iota) = Sum_j (-psi')^j Res_x(alpha(x)^T x^(-j-1)) eta^-1 = Sum_j (-psi')^j ([x^j] alpha)^T
eta^-1 = 0`, since `alpha` has no non-negative power of `x` (the sum over `j` converges `psi'`-adically, and `Res_x` commutes with it).
*Step 4 (the second term).* Put `beta(x) := R^-1(x)^T alpha(x)`, a product of elements of `Q[[x]][[u]]` and `Q[x^-1][[u]]`, hence in
`Q((x))[[u]] (x) H`, and `beta_(>=0)(y) := Sum_(j >= 0) y^j [x^j] beta`. Since `R^-1(psi')^T` does not involve `x`,
`Res_x(alpha^T R^-1(x) eta^-1 R^-1(psi')^T iota) = Sum_j (-psi')^j ([x^j] beta)^T eta^-1 R^-1(psi')^T = beta_(>=0)(-psi')^T eta^-1 R^-1(psi')^T`.
*Step 5 (symplecticity).* (Sp) gives `R(x)^-1 = eta^-1 R(-x)^T eta`, so `R^-1(x)^T = eta R(-x) eta^-1` (`eta` symmetric) and
`beta(x) = eta R(-x) f(-x) = eta (R f)(-x)`, where `(R f)(-x)` is the element `R f` of `HH^u` with `z -> -x`. Therefore
`[x^j] beta = (-1)^j eta [z^j](R f)` and `beta_(>=0)(-psi') = eta Sum_(j >= 0) psi'^j [z^j](R f) = eta (R f)_(>=0)(psi')`.
Transposing Steps 2-4: `C(psi') = -R^-1(psi') eta^-1 eta (R f)_(>=0)(psi') = -R(psi')^-1 (R f)_(>=0)(psi')`.
*Step 6 ("iff"; second slot).* `R(psi')^-1` is invertible in `gl_5(Q[[u, psi']])`, so `C = 0` iff `(R f)_(>=0) = 0` iff `R f` lies in `HH^u_-`
(the splitting `HH^u = HH^u_+ (+) HH^u_-` of `GAPS_G3.md` s.2.1). The second slot follows from the symmetry of Step 0.
*Step 7 (applications).* For `f = z^-k J^C(v(u))`, `k >= 1`, the hypothesis `f` in `HH^u_-` is `lem:tra-conv` (c); for `R = P_0^-1 P`,
`(R f)_(>=0) = P_0^-1 (P f)_(>=0)` since `P_0` does not involve `z`. QED

## 3. Remarks

1. *Where the source's ring is right.* If `R` is a polynomial in `z` modulo every power of `u` (the hatted rows: `[u^M] Ghat` is a
   polynomial in `zhat`, GH s.1), then modulo `u^(M+1)` all objects lie in `Q((x^-1))[[psi']]` and no decay of `f` is needed; this is the
   hatted contraction lemma (GH s.3, "at each u-order both sides are finite sums"), whose covector does not decay (GH s.2). For the
   note's flat Lemma E neither condition of that setting holds; `AA` with `f` in `HH^u_-` is the correct statement. Both rings contain
   `iota eta^-1`; only `AA` contains `iota R^-1(x)` in the flat case.
2. *The other expansion.* In `Q((psi'))[[u, x]]` one has `1/(x + psi') = Sum_j (-x)^j psi'^(-j-1)`, and the split of `E_R` into two
   terms is different (the first term is then non-zero after `Res_x`); the total is the same, but Step 3 fails term by term. The
   statement must therefore name the ring in which the split is taken: `AA`.
3. *What the hypothesis excludes* (computed in s.4, control C1): for `f` with constant (non-decaying) coefficients and `R(0, z)` a
   full power series, the contraction `Sum_m` does not converge (its `u^0 psi'^0` value changes with every truncation), and `R f` is not
   defined; the decay of `f` is a genuine hypothesis of Lemma E, not a convenience.

## 4. Computations (evidence only; exact over `Q`; scripts and logs in this folder, `py -I -B`, stdlib)

| script (log) | what is computed | result |
|---|---|---|
| `g3_frames.py 80` (`g3_frames_80.log`; also `_40`) | `Y^C` from the closed form of `lem:gau-vanshape` and, independently, as the regular solution of `z dY_b/dv = bmod_b Y_(b+1)`, `Y(0) = -e_1`, by iterated integrals; frames `k <= 80`, `a != 1`, by `Theta^-1 = Sum (-N)^r` (part 2 Step 1-3) and by the source's top-down elimination | the two `Y^C` agree (`n < 96`); 324 frames agree; `Pi_i(frame) = 0`, `0 <= i <= k + 12`; all exponents integral and `>= 1`; the outline's examples reproduced. Control C1 (frames built from `Y^C(+z)`): 304 of 314 non-trivial frames are not annihilated |
| `g3_decomp.py 10 30` (`g3_decomp_10_30.log`; also `_8_24`) | u-adic setting, `v = u tau^(-5/4)`, mod `u^10`: `Theta` invertible; `pr(q)` for a random `q` of `z`-degree 30; frame reconstruction; exact u-orders of all frame coefficients `k <= 90`; the infinite element `Sum_k z^k 1`; a toy bilinear form with the degree bound (`h = 2`, `n = 2`) including singular exponents | `Pi_i(pr q) = 0` (`i <= 33`); `pr(q) = Sum q_(k,a) frame(k,a)` on all coefficients; 3312 u-orders equal `(4(k-j)+nu_a+1)/5`; the `[u^9 z^j H]` coefficients of `pr` stop changing exactly at the predicted truncation; the toy correlator of frame truncations converges with `ord(mu(X^(K+5)) - mu(X^(K))) = 1, 5, 9, 13, 17` (the bound `4K/5 - 3`, sharp). Controls with `v = 1 + u` (`v(0) != 0`): the exact partial sums of `[u^0 z^0 H]` change at every `K = 1 mod 5` up to 60, and the toy differences stay of order 0: no limit, as predicted |
| `g3_lemE.py 8 6` (`g3_lemE_8_6.log`; also `_6_5`) | Lemma E with `R = exp(X)`, `X = Sum_(k<=3) z^k X_k(u)` infinitesimally symplectic, random integers; `R(0, z)` has all 23 computed `z`-powers non-zero; `f = z^-1 J^C(v(u))`, `z^-2 J^C`, and a random `f` with `[z^-n] f = O(u^ceil(n/2))`; both sides mod `u^8`, `psi'^0..psi'^6`, every `z`-truncation exact | (Sp) through `z^23`; LHS = RHS on all coefficients for the three `f` (all 35 components non-zero); divisibility of the numerator. Controls: C3 (sign `psi' -> -psi'`) fails for all three; C1 (non-decaying `f`, `L = 8..28` terms) gives 5 different `u^0 psi'^0` values (no limit), each truncation satisfying the finite identity; C2 (non-symplectic `R`) breaks the divisibility |

Helper modules: `g3_series.py` (truncated series, `v(u) = u tau^(-5/4)`), `g3_mat.py` (matrices, `exp`, random infinitesimally
symplectic `X`). Run from this folder, e.g. `py -I -B g3_lemE.py 8 6`; the logs print no paths. Total run time below one minute.
The computations test the bookkeeping of the proofs (orders, signs, the role of the hypotheses); they are not used by any proof.

## 5. Files of this folder

`GAPS_G3.md` (results, sources searched, conventions, observations), `GAPS_G3_2.md` (G3-1, G3-2g), `GAPS_G3_3.md` (G3-2a-f),
`GAPS_G3_4.md` (G3-3, computations); `g3_frames.py`, `g3_series.py`, `g3_decomp.py`, `g3_mat.py`, `g3_lemE.py`; logs
`g3_frames_40.log`, `g3_frames_80.log`, `g3_decomp_8_24.log`, `g3_decomp_10_30.log`, `g3_lemE_6_5.log`, `g3_lemE_8_6.log`.
No existing file was modified; nothing was written outside this folder and the scratch folder.
