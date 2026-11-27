# Stage 5 critic: a proof of (CYC), the cyclic form of the formal quintic's system (2026-11-14)

Folder `allgenus/stage5/critic/` (critic of the Stage-5 workflow; new files only). AI-written by the critic; NOT refereed; no human
has read it. Labels as in the brief: PROVED (argument written out here), VERIFIED (exact computation, range stated).
**Why this note.** The pb audit (`../pb/PB_AUDIT.md` s.2.3, `PB_AUDIT_2.md` P1) found that Theorem P(b) uses (CYC) in step 5, that
(CYC) is missing from every hypothesis list (POLY s.6, the poly referee's C6, STAGE4_SYNTHESIS R5, MECHANISM_STAGE4 s.2), and that it
is an implicit input of the whole programme (UR s.0 item 1 labels the cyclic form VERIFIED). pb gave an outline (the kernel of the
factorised operator) and an exact operator check through `u^56`, but no proof and no citable statement number. This note proves it.

## 0. Statement

Notation of `natural_form/fq_conifold.py` (docstring, read-only) and PB_AUDIT s.2.3. `u = 1 - 5^5 q`, `D = q d/dq = -(1-u) d/du`,

    PF4 = D^4 - 5q prod_(k=1..4) (5D + k) = Sum_(k=0..4) A_k D^k,
    A_4 = u, A_3 = -2(1-u), A_2 = -(7/5)(1-u), A_1 = -(2/5)(1-u), A_0 = -(24/625)(1-u).

At the conifold `u = 0` the solutions are `f_0, f_1, f_2` (power series, `f_k = u^k + O(u^3)`) and `y_3 = g log u + h` (`h` a power
series), where `g = f_1 + gamma f_2` is the vanishing period (`gamma = 7/10`). A *conifold frame* is `varpi_0 = y_0` any power-series
solution with `y_0(0) != 0` (fq_conifold: `varpi_0 = f_0 + alpha g + beta f_2`) and `t = g/y_0` (`t_c` of UR), so `Dt` is a unit.

**Theorem CYC (PROVED below, s.1-s.7).** In every conifold frame,

    PF4 = M := (1/y_0) D (1/Dt) D (u y_0^2 (Dt)^2) D (1/Dt) D (1/y_0)        (operators with coefficients in Q((u))).

**Corollary (PROVED).** With `c = (Dt, L^5/(y_0^2 (Dt)^2), Dt, y_0, y_0)`, `L^5 = 1/u`, the cyclic system `z D Y_a = c_a Y_(a+1)`
(`a` in Z/5) of UR s.1 holds for `(Y_0, .., Y_4)` iff `Psi := y_0 Y_0` solves `PF5 = M^5 - 1 - q prod_(k=1..5) (5M + kz)` (`M = zD`) and
`Y_(a+1) = z D Y_a / c_a` for `a = 0..3`. *Proof.* The four definitions and the closing relation give
`Y_0 = z^5 (1/c_4) D (1/c_3) D (1/c_2) D (1/c_1) D (1/c_0) D Y_0`, i.e. `z^5 D M Psi = Psi`. Since `D o q = q o (D + 1)`,
`q (5D + 5) = 5 D q`, so `PF5 = z^5 D PF4 - 1` (PB_AUDIT s.2.2), and `PF5 Psi = 0` reads `z^5 D PF4 Psi = Psi`. Apply the Theorem. QED

So the closing relation `4 -> 0` of the quintic rows (fq_rows, UR s.1) is a theorem, in every frame and for every normalisation of `g`.

## 1. The ring and its constants

`K := Q((u))` with the derivation `D`; `R := K[l]`, `l` an indeterminate with `D l := 1 - 1/u` (`l = log u`). The solutions of PF4 used
below lie in `R`; their Q-span is `S := span(f_0, f_1, f_2, y_3)`, with `y_3 = g l + h`.
*Frobenius structure (PROVED by the recursion).* PF4 maps `Q[[u]]` into itself and `[u^(n-3)] PF4 u^n = n (n-1)^2 (n-2)`, non-zero for
`n >= 3`; so the coefficients of `u^0, u^1, u^2` are free and determine a power-series solution. For `g l + h`, `PF4(g l) = r(g)` has no
`l`-term, and `PF4 h = -r(g)` is solvable in `Q[[u]]` iff the coefficients of `u^-3, u^-2, u^-1` of `r(g)` vanish: finitely many
coefficients of `g`, a linear condition that holds exactly for `g = f_1 + (7/10) f_2` (computed, `cr5_cyc_run`).
**Lemma 0 (constants, PROVED).** `{f in R : D f = 0} = Q`. *Proof.* For `f in K`, `D f = 0` means `f' = 0`, so `f in Q`. For
`f = Sum_(k<=n) a_k l^k` with `a_n != 0`, `n >= 1`: the `l^n` coefficient of `Df` is `D a_n`, so `a_n in Q^x`; the `l^(n-1)` coefficient
is `D a_(n-1) + n a_n (1 - 1/u)`, and `D a = n a_n (1-u)/u` would need `a' = -n a_n/u`, impossible in `Q((u))`. QED

## 2. Self-adjointness

**Lemma C1 (PROVED).** The formal adjoint `PF4* := Sum_k (-D)^k o A_k` equals PF4.
*Proof.* `PF4 = D^4 - 5^5 q o P(D)`, `P(D) = prod_(k=1..4) (D + k/5)`. Since `D o q = q o (D + 1)`, `(q o P(D))* = P(-D) o q =
q o P(-D - 1)`, and `P(-D - 1) = prod_(k=1..4) (-(D + (5-k)/5)) = P(D)` (four sign changes; `k -> 5 - k` permutes the factors). With
`(D^4)* = D^4`: `PF4* = PF4`. QED (VERIFIED exactly on the coefficients, `cr5_cyc_run` (a); the control `PF4 + u D` is not.)

## 3. The form B

**Lemma C2 (PROVED).** For `L = Sum_(k<=4) A_k D^k` put `P(f, g) := Sum_(k=1..4) Sum_(j=0..k-1) (-1)^j D^j(A_k g) D^(k-1-j) f`. Then
(i) `g L f - f L* g = D P(f, g)` for all `f, g` (an identity of differential polynomials);
(ii) for `L = PF4` and `f, g in S`, `B(f, g) := P(f, g)` lies in Q; `B` is bilinear and alternating;
(iii) `P(f, g) = Sum_(a,b) D^a f . P_ab . D^b g` with `P_ab = 0` for `a + b >= 4` and `P_(3-j, j) = (-1)^j u`;
(iv) in the basis `(f_0, f_1, f_2, y_3)`, `Bmat = Wr^T P Wr` with `Wr = (D^a y_i)`, and `det Bmat = 4`; hence `B` and `Wr` are invertible.
*Proof.* (i) For each `k`, with `H = A_k g`: `D[(-1)^j D^j H D^(k-1-j) f] = (-1)^j D^(j+1) H D^(k-1-j) f + (-1)^j D^j H D^(k-j) f`; the sum over
`j = 0..k-1` telescopes to `H D^k f - (-D)^k(H) f`. (ii) By C1 and (i), `D P(f, g) = 0` for solutions, so `P(f, g) in Q` (Lemma 0).
`Q := P(f, g) + P(g, f)` satisfies `D Q = 0` as a differential polynomial in the indeterminates `f, g`. Write `Q = Sum q_ab f_a g_b`
(`f_a = D^a f`). If `s = max{a + b : q_ab != 0}`, the coefficient of `f_(a+1) g_b` in `D Q` (`a + b = s`) is `q_ab + q_(a+1, b-1)`; the
one of `f_0 g_(s+1)` is `q_(0,s)`. All vanish, so inductively every `q_ab` with `a + b = s` is 0, a contradiction; hence `Q = 0` and `B` is
alternating. (iii) The order of `D^j(A_k g) D^(k-1-j) f` is `k - 1 <= 3`; the order-3 part comes from `k = 4`, `A_4 = u`, the term
`(-1)^j u D^j g D^(3-j) f`. (iv) The matrix identity is (iii) evaluated on the basis. `Bmat` is a constant matrix, computed exactly
from the solutions (`cr5_cyc_run` (b)):
`Bmat = [[0, -7/5, 2, -9/10], [7/5, 0, 0, -1], [-2, 0, 0, 0], [9/10, 1, 0, 0]]`, `det = 4`. As `det P = u^4` (anti-diagonal
`u, -u, u, -u`), `det Bmat = (det Wr)^2 u^4`, so `det Wr != 0`. QED

## 4. Griffiths transversality

**Lemma C3 (PROVED).** For every Q-basis `(y_i)` of S: `Sum_(i,j) (Bmat^-1)_ij D^a y_i D^b y_j = (P^-1)_ab`. It vanishes for `a + b <= 2`,
and `(P^-1)_(0,3) = 1/u`. In particular `T_k := Sum (Bmat^-1)_ij y_i D^k y_j` is 0 for `k = 0, 1, 2` and `T_3 = 1/u`.
*Proof.* From C2 (iv), `P^-1 = Wr Bmat^-1 Wr^T`. Let `J` be the 4x4 exchange matrix. `(P J)_(ac) = P_(a, 3-c)` vanishes for `c < a`
(then `a + 3 - c > 3`), so `P J` is upper triangular with diagonal `P_(c, 3-c)`; so is `(P J)^-1 = J P^-1`, with diagonal
`1/P_(c, 3-c)`. Its entries `(J P^-1)_(c a) = (P^-1)_(3-c, a)` vanish for `a < c`, i.e. `(P^-1)_(ba) = 0` for `a + b < 3`; and
`(P^-1)_(0, 3) = 1/P_(3, 0) = 1/u`. QED

## 5. Isotropy of every conifold frame (formal Picard-Lefschetz)

**Lemma C4 (PROVED).** `B(f, g) = 0` for every power-series solution `f`. In particular `B(y_0, g) = 0` in every conifold frame.
*Proof.* Let `tau` be the K-algebra automorphism of R with `tau(l) = l + 1`. It commutes with D (`D l in K`), so `tau(S) = S`,
`tau(P(f, f')) = P(tau f, tau f')`, and hence `B(tau f, tau f') = B(f, f')` (B is a rational number). Now `tau f = f` for power series and
`tau y_3 = y_3 + g`. So `B(f, y_3) = B(f, y_3 + g) = B(f, y_3) + B(f, g)`. QED (Consistent with fq_conifold's recorded check
"B(f_a, g) = 0"; VERIFIED in `cr5_cyc_run` (b): `B(f_i, g) = 0`, `i = 0, 1, 2`; `B(y_3, g) = 1`.)

## 6. Special geometry

**Lemma C5 (PROVED).** Fix a conifold frame `y_0`, `y_1 := g`, `t = y_1/y_0`, and `th := (Dt)^-1 D` (a derivation of R). There are
`y_2, y_3' in S` such that, with `f_i := y_i/y_0` (`f_0 = 1`, `f_1 = t`, `f_3 := y_3'/y_0`) and `Y := th^2 f_2`,
(i) `th^2 f_3 = -t Y`, so `th^2 (y/y_0) in Y . span_Q(1, t)` for every `y in S`; (ii) `Y = 1/(u y_0^2 (Dt)^3)`.
*Proof.* By C4, `span(y_0, y_1)` is B-isotropic, hence Lagrangian in the 4-dimensional symplectic space `(S, B)` (C2). Linear algebra
gives `y_2, y_3'` with `B(y_0, y_3') = B(y_1, y_2) = 1` and `B(y_0, y_2) = B(y_1, y_3') = B(y_2, y_3') = 0`. In this basis `Bmat` is the
signed permutation matrix with `Bmat_03 = Bmat_12 = 1 = -Bmat_30 = -Bmat_21`, orthogonal, so `Bmat^-1 = -Bmat`. Put
`s_r := Sum (Bmat^-1)_ij f_i D^r f_j` and `tau_k := Sum (Bmat^-1)_ij f_i th^k f_j`; `s_0 = tau_0 = 0` (antisymmetry). Leibniz on `y_j = y_0 f_j`:
* `T_1 = y_0^2 s_1`, `T_2 = 2 y_0 (D y_0) s_1 + y_0^2 s_2`, `T_3 = 3 y_0 (D^2 y_0) s_1 + 3 y_0 (D y_0) s_2 + y_0^2 s_3`;
* with `D = Dt . th`: `s_1 = Dt tau_1`, `s_2 = D(Dt) tau_1 + (Dt)^2 tau_2`, `s_3 = D^2(Dt) tau_1 + 3 Dt D(Dt) tau_2 + (Dt)^3 tau_3`.
C3 gives `T_1 = T_2 = 0`, `T_3 = 1/u`; `y_0` and `Dt` are units. Hence `tau_1 = tau_2 = 0` and `tau_3 = 1/(u y_0^2 (Dt)^3)`.
In the adapted basis, `tau_k = -th^k f_3 + f_3 th^k f_0 - t th^k f_2 + f_2 th^k f_1`; for `k >= 2` (`th^k 1 = th^k t = 0`):
`tau_2 = -th^2 f_3 - t Y`, so (i) follows from `tau_2 = 0`; and `tau_3 = -th^3 f_3 - t th Y = (Y + t th Y) - t th Y = Y`. QED

## 7. Proof of Theorem CYC

Since `D (1/Dt) D = Dt . th^2` as operators, `M = (Dt/y_0) th^2 o (u y_0^2 (Dt)^3) o th^2 o (1/y_0)`. For `y = Sum c_i y_i in S`,
C5 gives `u y_0^2 (Dt)^3 th^2 (y/y_0) = c_2 - c_3 t`, which `th^2` kills: `M y = 0` on S. The coefficient of `D^4` in M is
`(1/y_0)(1/Dt)(u y_0^2 (Dt)^2)(1/Dt)(1/y_0) = u = A_4`. So `E := M - PF4 = Sum_(k<=3) e_k D^k`, `e_k in Q((u))`, kills the basis
`(f_0, f_1, f_2, y_3)`: `E y_i = Sum_a e_a Wr_(a i)`, so the row vector `(e_0, .., e_3)` annihilates the invertible 4x4 matrix
`Wr = (D^a y_i)` from the left (C2 (iv)). Hence `E = 0`. QED

*What the proof uses.* Only the explicit PF4, the Frobenius structure at `u = 0` (three power-series solutions and one solution
`g l + h`; computed, `gamma = 7/10` found), and finite algebra. No analysis, no monodromy of actual functions (the automorphism
`tau` replaces it), no Hodge theory, and nothing from the project except the definition of the frame. The Yukawa coupling
`th^3 F` of special geometry is `Y = 1/(u y_0^2 (Dt)^3)` here; its normalisation `1/u` is fixed by C3, i.e. by `A_4 = u`.

## 8. Checks (VERIFIED; `cr5_cyc.py` + `cr5_cyc_run.py`, exact over Q, logs `logs/cr5_cyc_run_N.log`, N = 24, 48, 80)

| check | content | range (N = 80) | result |
|---|---|---|---|
| set-up | PF4 maps power series to power series; indicial `n(n-1)^2(n-2)`; `f_0 = 1 + (2/625)u^3`, `f_1 = u - (19/60)u^3`, `f_2 = u^2 + (37/30)u^3` (as fq_conifold/UR); `gamma = 7/10` FOUND from the `u^-1` obstruction; `y_3` solves PF4 | to `u^73` | PASS |
| (a) C1 | `PF4* = PF4` exactly | exact | PASS; control `PF4 + uD` not self-adjoint |
| (b) C2, C4 | `B(y_i, y_j)` log-free and constant; alternating; `B(f_i, g) = 0` (i = 0, 1, 2); `B(y_3, g) != 0`; `det Bmat = 4` | to `u^72` | PASS; control: top-order part of B alone is not constant |
| (d) C4 in a frame | `B(varpi0, g) = 0` for `varpi0 = f_0 + f_1/3 - 2f_2/7`; the control pairs have `B(varpi0, f_1) = -7/5`, `B(varpi0, f_1 + f_2/2) = -2/5` | exact | PASS |
| (c) C3 | `T_0 = T_1 = T_2 = 0`, `T_3 = 1/u` | to `u^71` | PASS |
| (e) Theorem | `M y = 0` for `f_0, f_1, f_2, y_3`, and `M u^n = PF4 u^n` (`n = 0..4`), frames `(a, b) = (0,0), (1/3,-2/7), (7/11,-3/5)` (`y_0 = f_0 + a f_1 + b f_2`) | to `u^70` | PASS; controls `t = f_1/y_0` and `t = (f_1 + f_2/2)/y_0` (non-isotropic pairs) fail |

Independent of pb's `pb_cyc.py` (which checks the operator identity on monomials through `u^56`): different code, different method
(solutions, the form B, transversality), and each lemma of the proof is tested separately. My own bug: the first run (`N = 16`)
crashed in the control `t = f_2/y_0`, whose `Dt` is not a unit (`f_2 = O(u^2)`); the control was replaced by `t = f_1/y_0` before any
recorded run. A left-over log line still quoted `B(varpi0, f_2)`; it was replaced by check (d) and the three logs were re-run from the
final code (overwriting my own earlier logs of the same runs, same verdicts). No other failure.

## 9. Remarks

* This is classical special geometry for a fourth-order self-adjoint operator, in a frame that is isotropic for the concomitant
  (Strominger; Candelas-de la Ossa; Ceresole-D'Auria-Ferrara-Lerche-Louis, "Picard-Fuchs equations and special geometry",
  arXiv:hep-th/9204035). I have NOT read these sources and cite no statement from them; the proof above is self-contained.
* The frame enters only through `B(y_0, g) = 0`, which holds for EVERY power-series `y_0` (C4). So "(CYC) in every conifold frame"
  is not a frame condition; a non-isotropic pair (e.g. `t = f_1/f_0`) violates it, as the controls show.
* Status: PROVED (critic; AI-written; unrefereed) and VERIFIED twice (pb, critic). It discharges the hypothesis H3 of PB_AUDIT s.1
  and the "VERIFIED" label of UR s.0 item 1.
