# Step (V) of ML_K, part 2: the proof, the corollary, verification, status (2026-11-14)

Continues `LOGLEG_PROOF.md` (s.0-2; same labels and notation). AI-written; no human has read it.

## 3. Proof of Theorem V

**Recursion (V1) (PROVED given (Q1)).** For `h >= 0`, `a, k >= 0`, `n = a + k`, `2h - 2 + n > 0`, `(h, n) != (0, 2)`, legs `Phi_i in Q[v][p]`:

    <Phi_1..Phi_a, (-log p)^(k+1)>_h = d_v <Phi_1..Phi_a, (-log p)^k>_h - Sum_i <.., d_v Phi_i, .., (-log p)^k>_h
                                      + Sum_{i=1..a} <.., (xi dPhi_i)_+, .., (-log p)^k>_h,

`(xi dPhi)_+` = the quotient of `p Phi'(p)` by `p^5 - v` (s.2.2(a)); for `Phi = p^m`: `m Sum_{q>=0, m-5-5q>=1} v^q p^(m-5-5q)`.
*Proof.* Fix `v_0 != 0`; for `v` near `v_0` choose `R > max_b |c_b(v)|` and a ray `gamma` avoiding all `c_b(v)`. Then
`C(v) := <Phi, (-log p)^k>_h = (2 pi i)^-a oint_{|p_i|=R} Int_{gamma^k} prod_i Phi_i(p_i) omega_{h,n}`, all contours independent of v.
Differentiate under the integrals (the integrand is holomorphic in v on the compact contours; on the rays it is bounded uniformly in v
near v_0 by `const |dq|/(1 + |q|)^2`, by (F1) and continuity of the coefficients): `d_v C = Sum_i <.., d_v Phi_i, ..> + (integral of
prod Phi_i d_v|_p omega_{h,n})`. By s.2.1 and (Q1), `d_v|_p omega_{h,n} = Int_gamma omega_{h,n+1}(., q) + Sum_{l=1..n} Lie_xi_l omega_{h,n}`.
The first term is the pairing with one more log leg (Fubini, (F1)); a Lie term in an L-slot gives `-<.., (xi dPhi_i)_+, ..>`
(s.2.2(a)); a Lie term in a log slot gives 0 (s.2.2(b)). Rearranging gives (V1). QED

**Proof of (i).** Fix `h` and `a >= 1`; induction on `k`. `k = 0`: L-regularity (Q2) (and, for `(h, n) = (0, 3)`, the explicit formula
of s.2.4). Step `k -> k+1`: if `(h, a + k) = (0, 2)` the pairing is a `(0,3)` pairing (`(a, k+1) = (2, 1)` or `(1, 2)`), given explicitly
by s.2.4 and in `Q[v]`. Otherwise every term on the right of (V1) is a pairing with the same `a >= 1` L-legs (legs `d_v Phi_i` and `(xi dPhi_i)_+`
again lie in `Q[v][p]`, since `p^5 - v` is monic) and `k` log legs, or the v-derivative of one: in `Q[v]` by the induction hypothesis.
QED  (The induction never leaves the class `a >= 1`: the L-legs are never consumed, so a log leg never meets a "bare" correlator.)

**Proof of (ii).** For `a = 0`, (V1) reads `<(-log p)^(k+1)>_h = d_v <(-log p)^k>_h`. Bases: `h >= 2`: `<-log p>_h = Int_gamma omega_{h,1}
= d_v F_h` by (Q1) at `n = 0`, and `F_h = B_2h/(2h(2h-2)) v^(2-2h)` by (Q3). `h = 1`: by homogeneity `<-log p>_1 = c/v` with one rational
number `c`, computed exactly over `Q(zeta5)`: `c = -1/12` (`ll_falsifier_exact.log`, (1,1)); this is a proof (one exact number), and
it fixes `F_1 = -(1/12) log v` up to a constant. `h = 0`: `<(-log p)^3>_0 = 1/v` (s.2.4) `= F_0'''` for `F_0 = (1/2) v^2 log v`. QED

**Corollary ML_K (PROVED modulo (Q1)-(Q3) and D(iii)).** Let `X_1, .., X_n in K` be homogeneous (`n >= 1`, `2h - 2 + n > 0`), write
`X_i = Y_i + pi_i H` with `Y_i in K_0`, `pi_i = Pi_0(X_i)` (on the model basis `pi = 1` for `tau_0(H)` and 0 for `kappa(X,k)`). D(iii)
gives legs `Phi_i in L` with `q_Phi_i = Y_i`; by homogeneity a weight-w vector of `K_0` has the leg `Sum_j alpha_j v^j p^(w-5j)` mod
`Q[v][x]` with constant `alpha_j`, so `Phi_i in Q[v][p]` up to functions of x. By multilinearity, s.1 and Theorem V,
`<X_1..X_n>_h = Lambda^(1-h) Sum_{S subset [n]} prod_{i not in S} pi_i <(Phi_i)_{i in S}, (-log p)^(n-|S|)>_h`; the `S = {}` term is
`Lambda^(1-h) F_h^(n) prod pi_i`, which is the frozen `F_h^(n) prod pi_i` (frozen `F_h = c_h v^(2-2h)`, `F_1 = -(1/12) log v`,
`F_0 = -(5/2) v^2 log v`, as in MODEL_SINGULAR s.5), and every other term is in `Q[v]`. This is ML_K exactly as stated in
MODEL_SINGULAR s.5 ((ML1) and (ML2)) and as used in quintic Theorem R. QED

**What the critic's sketch (SYNTHESIS s.6 item 3) needed, and what changed.** The route is the critic's. Made precise here: the
pairings are ancestor correlators (s.1) and `-log p` is exactly the ancestor `tau_0(H)`; the hypotheses of (Q1) (s.2.3); the sign of
the xi-term (`+`, not `-`); the step never removes an L-leg, so the induction is on `k` at fixed `a`; "genus 0 from `omega_{0,2} = B`"
is the case (0,2) of (Q1), proved directly (s.2.4), and the (0,3) pairings are explicit; the bases `h = 1` (`-1/(12v)`) and `h = 0`
(`1/v`) are separate from the Gap Theorem.

## 4. Verification (`ll_eo.py`, `ll_vstep.py`, `ll_falsifier.py`; all exact mod p, rational values reconstructed and confirmed)

| step | identity checked | range | field | control (must fail) | result |
|---|---|---|---|---|---|
| (X) s.1 | `<-log p, d xif_{a,d}> = delta_{(a,d),(1,0)}` | a <= 4, d <= 8 | 2^89-1, 2^61-1 | - | PASS |
| (E0) s.2.4 | `Int omega_{0,3}(p1,p2,.) = Sum_b e_{b,1} e_{b,1}/(5 c_b^3) = d_v B|_x` | exact tensor | both primes | opposite sign | PASS, control fails |
| (E1) (Q1), n = 0 | `Int omega_{h,1} = d_v F_h` | h = 1..4 | both primes | - | PASS |
| (E2) (Q1) | `d_v|_p omega - Sum Lie_xi omega = Int omega_{h,n+1}` as e-basis tensors | 14 cases, (0,3)..(4,1) | both primes | no Lie term; Lie sign flipped | PASS, both controls fail in all 14 |
| (E3) s.2.2(b) | `Int_gamma Lie_xi omega_{h,n+1} = 0` | same 14 | both primes | - | PASS |
| (S) (V1) | the step on EVERY weight tuple `Sum m <= 5(2h-1+n) + 5` | 32 triples (h,a,k), 32559 tuples, h <= 4, n + 1 <= 6 | both primes | xi-term sign flipped; dropped | PASS; each control fails on 2..12210 tuples per triple |
| (D) remark | `P_{k+1} = (2-2h-n) P_k - (1/5) <.., p^5, ..>` | same 32559 | both primes | - | PASS |
| (i), (ii) | conclusions, every possibly singular pairing | `2h-2+n <= 6`, (3,3), (4,1), (4,2) | both primes; exact Q(zeta5) for `2h-2+n <= 3` and (3,1), (2,2), (1,3), (0,5) | off-L legs `p^-1, p^-2, p^-5, p^-10, p^-15` | PASS (FALSIFIER_B.md) |

The 14 cases of (E2) are `(h,n) -> (h,n+1)` for (0,3), (0,4), (0,5), (1,1), (1,2), (1,3), (1,4), (1,5), (2,1), (2,2), (2,3), (3,1),
(3,2), (4,1); the largest tensors are `omega_{1,6}` and `omega_{2,4}`.
Logs: `ll_eo_89.log`, `ll_eo_61.log` (530 s each), `ll_vstep_89.log`, `ll_vstep_61.log` (670 s each; identical counts at both
primes), `ll_falsifier_{new,low,exact}.log` (448 s, 32 s, 49 s).

## 5. Remarks

* **A second, partial route (no (Q1)).** The TR dilaton equation (arXiv:math-ph/0702045 Thm 4.7) with `Phi_ydx = p^5/5 - v log p`
  gives `v <J, -log p>_h = (2-2h-n) <J>_h - (1/5) <J, p^5>_h` for any legs J (VERIFIED, (D)). With (i) at k log legs this shows that
  `<J, -log p>` has at most a SIMPLE pole at `v = 0`; only (Q1) excludes the `v^-1` term. Comparing with (V1): the Euler-dilaton leg
  `-p^5/5` (= `e = tau_1(1) + v tau_0(H)`) acts on L-legs as the polynomial part of `-(4/5) (x/x') d/dp`, i.e. of the vector field
  `-(4/5) x d/dx`, together with the weight `-(Sum m_i)/5` (SKETCH; follows from (V1) and (D)).
* **Scope.** Theorem V concerns the frozen model at the point v (the vertices of the graph sum). Deformed curves never enter: in
  quintic Theorem R each vertex is undeformed, and the edges and leaves enter only through their K-components (Lemma K).
* **What remains for ML_K.** D(iii) for all weights (VERIFIED <= 240; the `dict` worker's task). The quoted inputs (Q1)-(Q3).
  Everything else on the model side of the mechanism is now PROVED modulo quoted theorems.
* **Hard rules.** Python `py -I -B`, stdlib only; the allgenus/tr modules are imported read-only (`tr_field`, `tr_core`, `tr_rec`,
  `tr_v0`, `tr_xi`); no existing file was modified; runs were made with the working directory in the scratch folder.
