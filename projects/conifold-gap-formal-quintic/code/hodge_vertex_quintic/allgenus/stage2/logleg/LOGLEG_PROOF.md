# Step (V) of ML_K: log legs keep polynomiality (Stage 2, logleg worker, 2026-11-14)

Folder `allgenus/stage2/logleg/`. AI-written and AI-checked only; no human has read it. Labels: PROVED (argument written here, or a
quoted theorem with arXiv number and statement number), VERIFIED (exact computation; range, field or primes stated), CONJECTURAL,
SKETCH, REFUTED. Conventions as in `../../SYNTHESIS_STAGE01.md` s.3: `H = d/dv`; "Gap Theorem" = `F_h = c_h v^(2-2h)`; `K_vc` = the
vanishing-covector slice `K` (`K_0 = K cap ker Pi_0`); "ML_K" = the Model Lemma on `K`; "quintic Theorem R". Model curve only:
`x = W = p^4/4 + v/p`, `y = p` (TR on this curve is allowed; nothing here touches the quintic's mirror curve). Companion files:
`FALSIFIER_B.md` (the falsifier, done first), `ll_core.py`, `ll_falsifier.py`, `ll_eo.py`, `ll_vstep.py` and their logs.

## 0. Result

**Theorem V (PROVED modulo three quoted theorems, s.3).** For `h >= 0`, `a, k >= 0`, `n = a + k`, `2h - 2 + n > 0` and legs
`Phi_1, .., Phi_a in Q[v][p]`, the TR pairing `<Phi_1 .. Phi_a, (-log p)^k>_h` (s.1) satisfies
* (i) if `a >= 1`: it lies in `Q[v]`, i.e. it is a polynomial in `v` (no pole at `v = 0`);
* (ii) if `a = 0`: it equals `F_h^(k)(v)`, the k-th derivative of `F_h = B_2h/(2h(2h-2)) v^(2-2h)` (h >= 2, the Gap Theorem in TR
  normalisation), `F_1 = -(1/12) log v` (k >= 1), `F_0 = (1/2) v^2 log v` (k >= 3);
* (iii) the recursion (V1) of s.3 holds, with the sign `+` in front of the xi-term (the critic's sketch had `-`; the `+` sign is
  forced by the derivation and VERIFIED, the opposite sign FAILS on every non-trivial tuple, `ll_vstep_*.log`).

Quoted inputs: (Q1) the variational formula of arXiv:math-ph/0702045, Theorem 5.1 (s.2.3; hypotheses checked; the genus-0 case is
re-proved here, s.2.4); (Q2) L-regularity, E9 of the synthesis, PROVED modulo arXiv:2212.00320 Theorem 7.1; (Q3) the Gap Theorem,
PROVED modulo arXiv:2502.18115 Corollary 1.2 (only for (ii), h >= 2).

**Corollary (ML_K from the TR side).** Translated by DNOPS (arXiv:1509.06954 Thm 8.1, eq. (2-35); `tr/TR_MODEL.md` s.2), Theorem V
says: for `X_1, .., X_n in K` (n >= 1, 2h - 2 + n > 0) the ANCESTOR correlator of the frozen theory at the point v satisfies
`<X_1 .. X_n>_h - F_h^(n)(v) prod_i Pi_0(X_i)` is regular at `v = 0` (frozen `F_h = c_h v^(2-2h)`), provided D(iii) holds
(`K_0 = L` modulo functions of x; VERIFIED weights <= 240 by the critic, not proved here). Status: PROVED modulo (Q1)-(Q3) and D(iii).
This closes item (V) of SYNTHESIS s.6 and leaves D(iii) as the only unproved model-side input of ML_K.

Evidence (VERIFIED, two primes 2^89-1, 2^61-1, cross-checked by rational reconstruction; plus exact Q(zeta5) for 2h-2+n <= 3):
the variational formula at the level of differentials for every `(h,n) -> (h,n+1)` through `(4,1) -> (4,2)`, `(3,2) -> (3,3)`,
`(2,3) -> (2,4)`, `(1,5) -> (1,6)`, `(0,5) -> (0,6)`, 14 cases (`ll_eo_*.log`); the induction step (V1) on every weight tuple in range
(`ll_vstep_*.log`); conclusion (i)-(ii) on every pairing that could be singular for all `2h - 2 + n <= 6` plus (3,3), (4,1), (4,2)
(`FALSIFIER_B.md`). Every check has a failing control.

## 1. The objects

* Curve family `C_v = (P^1, x_v = p^4/4 + v/p, y = p, B = dp1 dp2/(p1-p2)^2)`, `v != 0`. `dx = (p^5 - v) dp/p^2`: five simple branch
  points `c_b = zeta^-b v^(1/5)`, `x''(c_b) = 5 c_b^2 != 0`, `dy = dp != 0` there. TR in the convention of `tr/tr_rec.py`
  (arXiv:1509.06954 (2-26)-(2-29), the Eynard-Orantin kernel); `F_h` (h >= 2) by `(2-2h) F_h = Sum_b Res_{c_b} Phi omega_{h,1}`, `dPhi = y dx`.
* (F1) For `2h - 2 + n > 0`, `omega_{h,n}` is symmetric, has poles only at the branch points, with zero residues (arXiv:math-ph/0702045
  Thm 4.1-4.2, standard). In particular it is holomorphic at `p = 0` and at `p = oo`, where it is `O(dp/p^2)`. No diagonal poles.
* (F2) Homogeneity: `p -> s p`, `v -> s^5 v` gives `omega -> s^(5(2-2h-n)) omega` (TR_MODEL s.5). A pairing with legs `p^m_i` and
  `k` log legs is `const * v^e`, `e = 2 - 2h - n + Sum m_i/5`; Z/5 forces `Sum m_i = 0 mod 5`.
* **Pairings.** `gamma` = a ray from `0` to `oo` avoiding the branch points; `|p| = R` a circle enclosing them.
  - L-leg `Phi in Q[v][p]`: `<.., Phi, ..> := Sum_b Res_{c_b} Phi omega = (1/2 pi i) oint_{|p|=R} Phi omega` (no pole at 0 or inside
    except the c_b; constants pair to 0). Functions of x pair to 0 (linear loop equation), so `L = Q[v][p] + Q[v][x]` acts mod `Q[v][x]`.
  - Log leg: `<.., -log p, ..> := Int_gamma omega(.., p', ..) = Sum_b Res_{c_b} (-log p') omega` (cut along gamma; contour lemma of the
    xyswap referee: zero residues + regularity at 0, oo). Independent of the ray: moving it across `c_b` changes it by `2 pi i Res = 0`.
    In the e-basis `e_{b,k} = dp/(p - c_b)^(k+1)`: `Int_0^oo e_{b,k} = (-1)^k/(k c_b^k)`.
  - Mixed: L-slots on `|p| = R`, log slots on `gamma` (Fubini; no diagonal poles by (F1)).
* **Which correlators.** By DNOPS (2-35) (TR_MODEL s.2, VERIFIED there against MVertex6 on 72120 correlators),
  `omega_{h,n} = Lambda^(h-1) Sum <prod tau_{d_i}(H^{a_i})>_h prod d xif_{a_i,d_i}`, `Lambda = -5` (x = W), the bracket being the
  ANCESTOR correlator (integral over Mbar_{h,n} with its own psi classes) of the frozen CohFT at the point v. Hence
  `<Phi_1..Phi_n>_h = Lambda^(h-1) <q_Phi1 .. q_Phin>_h`, with the ancestor vector of a leg `q_Phi(z) := Sum_{a,d} <Phi, d xif_{a,d}> H^a z^d`.
  - The log leg is EXACTLY the ancestor insertion `tau_0(H)`: `<-log p, d xif_{a,d}> = xif_{a,d}(oo) - xif_{a,d}(0) = delta_{(a,d),(1,0)}`
    (for `d >= 1`, `d/dx = (p^2/(p^5 - v)) d/dp` kills the value at `p = 0`; `xif_{1,0}(0) = -1`, the others vanish at 0). PROVED;
    VERIFIED a <= 4, d <= 8 (`ll_vstep_*.log`, line (X)).
  - Calibration. `K = {q : Pi_i(q) = 0, i >= 1}` uses only the column `S(z) H^2` of the calibration (log free, fixed by homogeneity),
    so `K` and ML_K do not depend on the calibration ambiguity. The model's basis `kappa(X,k) = [S(z) X z^k]_+` (X != H) are the
    ancestor images of the descendants `tau_k(X)`; the bare ancestor insertion `z^k X` is off `K` in general, and its leg carries
    `log p . f(x)` with `f = -Sum_i Pi_i(z^k X)(-x)^i/i!` (D(ii)); the calibration terms `S_j X z^(k-j)` remove exactly this log, which
    is why the xyswap referee's "every descendant acquires log z" (for iterated x-primitives of the natural legs) is no obstruction.
    TR legs themselves need no calibration: an integration constant in an x-primitive changes the leg by a function of x.
  - Descendant form (remark, not used): on the slice `t^H_k = 0` (k >= 1) the descendant correlators of `tau_k(X)`, X != H, and
    `tau_0(H)` equal the ancestor correlators of `kappa(X,k)` and `H` for stable `(h,n) != (0,2)` (ancestor-descendant relation,
    arXiv:math/0108100 s.5); quintic Theorem R uses the ancestor form only.

## 2. Three lemmas

### 2.1 Fixed p versus fixed x (PROVED, elementary)
At fixed `x`, the point `p(x, v)` moves with `d_v p|_x = -(d_v x)/(d_p x) = -(1/p)/((p^5 - v)/p^2) = -p/(p^5 - v)`. For any family of forms,
`d_v|_x = d_v|_p - Lie_xi` in each variable, `xi := p/(p^5 - v) d/dp` (the critic's convention `d_v|_p = d_v|_x + Lie_xi`).
In particular `d_v(y dx)|_x = (d_v p|_x) dx = -dp/p`, and `Int_{q in gamma} B(p, q) = dp [1/(p - q)]_{q=0}^{q=oo} = -dp/p`. So
`d_v(y dx)|_x = Int_gamma B(., q)`, i.e. `Lambda = 1` on `gamma`: a third-kind variation (the residue `-v` of `y dx` at `p = 0`, `+v` at oo).

### 2.2 Integration by parts (PROVED, elementary)
(a) L-slot: `oint Phi Lie_xi w = -oint (xi . dPhi) w` (closed contour; `Lie_xi(Phi w) = d(i_xi(Phi w))`). With `p Phi'(p) = Q (p^5 - v) + r`,
`deg r < 5` (division by the monic `p^5 - v`), `xi . dPhi = Q + r/(p^5 - v)` and `oint (r/(p^5-v)) w = 0` since the integrand is
`O(dp/p^3)` at oo. So `oint Phi Lie_xi w = -<(xi dPhi)_+, w>` with `(xi dPhi)_+ := Q in Q[v][p]`; for `Phi = p^m`:
`(xi dPhi)_+ = m Sum_{q >= 0, m - 5 - 5q >= 1} v^q p^(m-5-5q)` (the constant term pairs to zero).
(b) Log slot: `Int_gamma Lie_xi w = [i_xi w]_0^oo = 0`, since `xi(0) = 0`, `w` holomorphic at 0, and at oo `xi = O(p^-4)`, `w = O(dp/p^2)`.
VERIFIED (E3) at the level of e-basis tensors for all cases of s.0 (`ll_eo_*.log`).

### 2.3 The variational formula (QUOTED, hypotheses checked)
**(Q1)** arXiv:math-ph/0702045 (Eynard-Orantin), s.5.1, Theorem 5.1, in the form used here: for a family of spectral curves with
simple branch points and fixed Bergman kernel such that `d_v(y dx)|_x = Int_gamma B(., q) Lambda(q)`, with `gamma` avoiding the branch
points and `Lambda` holomorphic on it, one has for `2h - 2 + n >= 0`, `(h,n) != (0,0)`:
`d_v|_x omega_{h,n}(J) = Int_gamma omega_{h,n+1}(J, q) Lambda(q)` and `d_v F_h = Int_gamma omega_{h,1} Lambda` (h >= 2).
Sign: for `(h,n) = (0,1)` the formula is the hypothesis itself; our kernel is EO's (TR_MODEL s.1), and the case (0,2) is re-proved
with this sign in s.2.4. Hypotheses: (a) simple branch points `c_b(v)` for all `v != 0`, `dy != 0` there: yes; (b) the variation is at
fixed x and B is the Bergman kernel of the fixed P^1 (genus 0: `dp1 dp2/(p1-p2)^2` in the global coordinate): yes; (c) `d_v(y dx)|_x`
is EXACTLY the third-kind differential `Int_gamma B` (s.2.1; no second-kind part), the "temperature" variation of EO s.5.1 between the
two poles 0 and oo of `y dx`; (d) the integrals converge: `omega_{h,n+1}` is holomorphic at the end points (F1). The residues at the
moving branch points commute with `d_v|_x` (small fixed circles in the x-plane around the moving branch values).
VERIFIED (E2): `d_v|_p omega_{h,n} - Sum_i Lie_xi_i omega_{h,n} = Int_gamma omega_{h,n+1}` exactly as e-basis tensors for all cases of
s.0, both primes; controls (no Lie term; Lie term of the opposite sign) fail in every case. (E1): `Int_gamma omega_{h,1} = (2-2h) F_h(1)
= d_v F_h` for h = 2, 3, 4 and `= -1/12` for h = 1. Independent proof of the n = 0 case (xyswap referee): dilaton + homogeneity +
L-regularity of the polynomial part of `Phi`.

### 2.4 Genus 0 from `omega_{0,2} = B` (PROVED)
`d_v|_p B = 0`, so `d_v|_x B = -(Lie_xi1 + Lie_xi2) B = -d1 d2 [(xi(p1) - xi(p2))/(p1 - p2)]` (xi read as the function `p/(p^5-v)`).
With `xi(p) = Sum_b r_b/(p - c_b)`, `r_b = 1/(5 c_b^3)` (no polynomial part), `(xi(p1)-xi(p2))/(p1-p2) = -Sum_b r_b/((p1-c_b)(p2-c_b))`, so
`d_v|_x B = Sum_b r_b e_{b,1}(p1) e_{b,1}(p2)`. On the other side the recursion gives `omega_{0,3} = -Sum_b e_{b,1} e_{b,1} e_{b,1}/x''(c_b)`
(leading order of the kernel `K = e_{b,1}(p0)/(2 dx) (1 + O(t))` and `B(sigma q, .) = -B(q, .) (1 + O(t))`; this is `-Sum Res BBB/(dx dy)` in
this convention, the sign that gives Witten-Kontsevich on `x = -p^2/2`), and `Int_gamma e_{b,1} = -1/c_b`, `x''(c_b) = 5 c_b^2`. Hence
`Int_gamma omega_{0,3}(p1, p2, .) = Sum_b r_b e_{b,1}(p1) e_{b,1}(p2) = d_v|_x B`: the case (0,2) of (Q1) with the right sign. VERIFIED (E0).
Consequences used below (v general, `Sum_b c_b^N = 5 v^(N/5)` if `5 | N`, else 0):
`<Phi1, Phi2, -log p>_0 = (1/5) Sum_b c_b^-3 Phi1'(c_b) Phi2'(c_b)`, `<Phi1, (-log p)^2>_0 = -(1/5) Sum_b c_b^-4 Phi1'(c_b)`,
`<(-log p)^3>_0 = (1/5) Sum_b c_b^-5 = 1/v`, `<Phi1,Phi2,Phi3>_0 = -(1/5) Sum_b c_b^-2 prod Phi_i'(c_b)`. For `Phi_i in p Q[v][p]` the
powers of `c_b` that occur are `>= -4`, and only multiples of 5 survive the sum, so every pairing with an L-leg is in `Q[v]`;
`<(-log p)^3>_0 = 1/v` is `F_0'''` for `F_0 = (1/2) v^2 log v`. VERIFIED: the e = 0 values `<p, p^4, -log p>_0 = 4`, `<p^5, (-log p)^2>_0 = -5`.
