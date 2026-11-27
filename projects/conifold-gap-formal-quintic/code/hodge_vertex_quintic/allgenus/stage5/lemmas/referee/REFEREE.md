# Referee report: Stage 5 / lemmas, "the formal stationary-phase lemmas" (2026-11-14)

Adversarial referee of `allgenus/stage5/lemmas/`: `FORMAL_SP_LEMMAS.md` and `_2` to `_5`, the code `sp_core.py`, `sp_iter.py`,
`sp_rand.py`, `chk_util.py`, `chk_*.py`, and the 15 logs. AI-written and AI-checked; no human has read this report. Labels as in the
note: PROVED, VERIFIED, CONJECTURAL, SKETCH, REFUTED. Referee code `referee/rl_*.py`: new, stdlib only, run as `py -I -B` from the
session scratch folder. It imports no project module and no worker module. Logs are in `referee/logs/`. Part 2: `REFEREE_2.md`.
Read in full: the five parts of the note; all worker code and logs; `MECHANISM_STAGE4.md`; `STAGE4_SYNTHESIS.md` s.0-s.3;
`POLY.md` s.5-s.8; the poly referee's `REFEREE.md`; `stage5/chart/CHART_PROOF.md`, `_2`, `_3` (the versions of 19:22-19:24; part 4
was not read); `stage5/pb/PB_AUDIT.md` s.2.0-s.2.7 (read only for how it uses the lemmas, not refereed). DLMF s.2.1, s.2.3, s.5.11
(opened 2026-11-14).

## 0. Verdict: sound, with corrections

* **No mathematical error found.** I checked every PROVED argument line by line (s.1-s.2).
  * The proofs are algebraic over an arbitrary commutative Q-algebra.
  * The one real integral, in (G5), is an exact Gaussian integral, transferred by polynomial identities. That transfer is correct.
  * Proof A of L3 is purely formal and complete. It supplies the formal proof that `CHART_PROOF_2` s.2.8 lists as "not done".
  * Proof B of L3 quotes DLMF (2.3.13)-(2.3.16), 5.11.1 and s.2.1(iii). The quotes match the sources (part 2, s.4).
* **Statements against their users** (s.3). The lemmas cover what the chart proof and Theorem P(b) use, with three exceptions:
  * **(F3), the z-rescaling (homogeneity), is not stated.** Only its special case (H) inside L3 appears. Users:
    * `CHART_PROOF` s.1.7, s.2.8, s.3.9 and s.3.10;
    * Theorem P(b) step 4 (`w/Y_c^4 = tau w`, poly referee s.3; PB_AUDIT s.2.6).
    It is trivial. I give a proof in s.3.2 (PROVED here).
  * **The relabelling `V = s Res_X(...)` of Theorem P(b) step 6** is the note's s.10.4.
    * The worker labels it SKETCH. The label is correct, and the item lies outside FL2-FL3.
    * PB_AUDIT s.2.7 writes the argument out. I did not referee it.
    * Until someone does, this step stays open in the assembled chain R5.
  * **The formal reflection identity of `CHART_PROOF` s.1.7** is not in the note. It identifies POLY's integrand with GAMMA's. It follows
    from L3(iv) and the oddness of `mu` (s.3.3).
* **Corrections C1-C8** (part 2, s.5) concern statement precision and wording. None changes a result.
* **Code.** The code is correct. I re-ran all 15 worker logs in a scratch copy: each is identical to the stored one, timings excepted.
  The checks test what the note says they test (part 2, s.6).
* **Own checks** (VERIFIED: exact over Q, each with failing controls; part 2, s.7):
  * a **toy chart**: a 4-variable analogue of the chart argument with log-type data and the amplitude `1/(x_0 x_1)`. Results:
    * total SP = order (a) = order (b) through `z^6` at `L = 3` and `L = -3/2`;
    * the toy (Gamma): `c_n` has no pole at `L = c_0` and takes the value `h_toy`, through `z^14`;
  * L3 by a third method (Lagrange-Buermann) for `m = -6..6` through `z^40`, and every step of Proof A. This includes the `eps -> z`
    substitution at the moving critical point, computed directly by Wick contraction;
  * FL4 (also iterated) and the shift `sigma -> sigma + z` of Theorem P(a), on Mellin-Barnes-type data with a moving critical point.
* **One observation on controls** (part 2, s.7.4). In the toy, the control "chart amplitude without the Jacobian" agrees with the true
  series at the special point `L = 2` through `z^9`, so a control can be blind at one parameter value. My final runs use generic points.

## 1. The proofs, item by item

| item | note | referee verdict | remarks |
|---|---|---|---|
| Lemma 1.1 (ring `R_n(A)`, exp) | s.1.2 | PROVED (checked) | finiteness per weight gives finite convolution sums; `P^m` has weight `>= m/2` |
| Wick recursion (W), `E_zC` | s.1.3 | PROVED | Leibniz on `(C xi)_i e^(..)`; `E_zC` preserves weight, finitely many inputs per output power |
| Lemma U, (U1), (U2) | s.1.5 | PROVED | the bound "orders `3..6(k - j0)`" is correct but not sharp (`2(k - j0) + 2` suffices, cf. CHART (F1)) |
| (G1)-(G4) | s.2 | PROVED | (G3): the multiplied form `E[t_k g] = z Sum_j C_kj E[d_j g]` is equivalent because `C` is invertible |
| (G5) reweighting | s.2 | PROVED | polynomial identity in the `n(n+1)/2` entries of `C_0`, open set of positive-definite matrices; domination by `e^(abs(t^T k t)/2)` for `norm(k) norm(C_0) < 1`; the formal `det^(-1/2)` is the Taylor series of the positive root |
| (G5') with `z` | s.2 | PROVED, scope remark C5 | `I - C_0 K_1 = H_0^-1 K`, `C_1 = -K^-1` checked. L1b uses it with `K_1 = K(y) - H_xx` depending on the outer variables and with `P` in `R_(n+m)(A)`; same proof, not stated |
| FL1 | s.3.1 | PROVED | via (G3) applied to `g_i = e^(F_3/z) h_i`; the Remark's weight count for Theorem P(a) (`>= (k+1)/2`) is right |
| formal IFT | s.3.2 | PROVED, wording C4 | the bracket is "a matrix with entries in `I`, times `(c - c')`", not "quadratic in `(c, c')`" |
| Lemma T (a), (b) | s.3.2 | PROVED | `dB/ds = Sum_i p_i (d_i B + z^-1 (d_i F)(c+t) B)` re-derived directly from the definition of `b_p` (no `e^(F/z)` heuristics needed); polynomiality in `s` per coefficient is right in both cases |
| Corollary T' | s.3.2 | PROVED | T(a) with amplitude `e^(delta/z) a`, the `F_(0,3)` cancellation, then (G5') with `K_1 = H_c - H_0` |
| L1c (i)-(iii) | s.4 | PROVED | (ii) is (U1) for the Taylor map and the evaluation; (iii) is T' |
| L1c (iv) (T-z) | s.4 | PROVED, wording C4 | the weighted ring with `eps` of weight 1 contains `S`, by (1); substitution `eps -> z` is a continuous ring map commuting with `E_zC`. "every term of `e^(eps g/z)`" should read "of `eps g/z`" |
| FL4 | s.4.5 | PROVED, hypothesis C2 | the final clause `P(zD_q) SP = SP[P(sigma) ..]` needs a `q`-independent amplitude (true in the Mellin-Barnes case) |
| L1a | s.5 | PROVED, constant C3 | linear step by (G2); homotopy `psi_s = y + s chi`; identities (a), (b) re-derived (Jacobi, Piola); FL1 at fixed Hessian. `P_s` must also subtract `G(0)` |
| Piola identity | s.5 | PROVED | closed `(n-1)`-forms `dpsi_1 ^ .. ^ dpsi_n` (one factor omitted); cofactor signs checked |
| FL3 for the SP (CFU) | s.5.4 | PROVED | L1a over `Q[[u]]`, `phi'(0) = Y'(1; u)` a unit, `f''(1) = 5`; `Bhat in Q[X, w][[u]]` is an admissible amplitude over `A = Q[[u]]` |
| L1b (Fubini) | s.6 | PROVED, constant C3 | straightening by L1a with `det D Theta = 1`; `H_xy = 0` afterwards; (G4) then (G5'); the "Rings" paragraph correctly identifies the inner functional with the inner SP over `A_y` |
| L2 (delta lemma) | s.7 | PROVED | `C = [[h, -1], [-1, 0]]`, `C_gg = 0`; the cubic part and the amplitude depend on `g` only; `(rho^2)^(-1/2) = rho^-1` in `1 + I` |
| L3 (i)-(v), Proof A | s.9.1-9.3 | PROVED | (H) power count `r - abs(alpha)/2 = -k`; A1 = FL1 with `h = x^m`; A2 = L1c(iv) with `g = log(x/a)` (all factors lie in the weighted ring, so `eps = z` may be substituted factorwise); A3 injectivity by leading terms; A4 both directions in `m` |
| L3(iv) | s.9.2 | PROVED | each displayed identity re-derived; `D^2` injective on `Y^-1 Q[[1/Y]]`; Q[[D]] acts because `D` raises the order |
| L3, Proof B | s.9.4 | PROVED modulo DLMF | quotes and hypotheses checked (part 2, s.4); left half `b_s = (-1)^s g_s/2` re-derived with orientation |
| Lemma 10.1, FL2, FL3 (residue) | s.10.1-10.3 | PROVED | ring `L` closed under products (bound `min(eps', eps'')abs(j) - M' - M''`); FL2 trivial in `L`; the telescoping identity re-derived; the Taylor shift is a continuous ring map, so the pull-back equals the naive one |
| s.10.4 relabelling | s.10.4 | SKETCH (label correct) | see s.0 and s.3.4 |
| s.8 interface | s.8 | correct | matches the chart proof's steps (s.3.1); warnings (1), (2) are right |

## 2. The delicate points, in detail

* **(G5).** The transfer is legitimate. For a fixed monomial in the entries of `k`, both coefficients are polynomials over Q in the
  entries of `C_0`. For real positive-definite `C_0` they are the Taylor coefficients at `k = 0` of two real-analytic functions that
  coincide for `norm(k) < 1/norm(C_0)`. Integrability and termwise integration both follow from `e^(norm(k) abs(t)^2/2) abs(t^alpha)`
  being integrable against `N_(C_0)`.
* **The weighted rings.**
  * `R_(n+m)(A)` is NOT contained in `R_n(A[[y]])`: e.g. `Sum_k x^(2k) y^(2k) z^(-k)` lies in the first ring but not in the second.
  * So the inner functional of L1b has to be identified with the inner SP coefficientwise, as the "Rings" paragraph of s.6 does. Its
    argument is right: each `y^beta z^k` receives finitely many terms, the same on both sides.
  * The same remark applies to (G5') in L1b (C5).
* **L1c(iv) is the one genuinely new device.**
  * Lemma T' puts `S(eps; z)` into the weighted ring with `eps` of weight 1, since every term of `eps g/z` has weight `>= 1/2`.
  * There `eps -> z` is a continuous ring map that commutes with `E_zC`.
  * The note uses it only in Step A2 of L3. The other two uses listed in s.4 ((S4), the shift of P(a)) are covered by Lemma T(b) as well.
  * I re-did Step A2 numerically without the closed form: the Wick SP at the moving critical point `1 + eps` over `Q[[eps]]`, then
    `eps -> z` (`rl_l3` (c)). The controls fail: without `(Delta_eps/Delta)^(-1/2)`, and with `eps -> 2z`.
  * This answers the objection of `CHART_PROOF_2` s.2.8, that shifting `a` by the formal variable `z` "is not a formal operation".
    The shift is realised as the perturbation `eps log(x/a)` of the phase, and the weight bound makes `eps = z` legitimate.
* **L1a.** Every coefficient of `Phi(s)` is a polynomial in `s`: by Lemma U it involves finitely many Taylor coefficients, each
  polynomial in `s`. So `Phi' = 0` gives `Phi(1) = Phi(0)` over a Q-algebra. Identity (b) needs the Piola identity for `M_s = D psi_s`,
  and the note proves it in exactly that form.
* **Branch rules.** The full-SP rules of L1a, L1b(4) and L2(4) are mutually consistent:
  `det(-J^T H J)^(1/2) = det J det(-H)^(1/2)`, `Delta^(1/2) = Delta_in(y_c)^(1/2) det(-S)^(1/2)`, and
  `(-G_t^2)^(-1/2) = (-g_c^2)^(-1/2) g_c/G_t`. The project normalises by the constant term, so only the z-dependence matters there.
* **L2.** The note proves exactness by a change of variables to `(s, g)`, so it needs L1a. The chart proof (`CHART_PROOF_2` s.2.7)
  counts degrees instead (every `xi` pairs with an `s`) and needs no change of variables. Both proofs are correct.
* **FL2 and FL3 for the residue are trivial in the ring `L`.** The content of Theorem P(b) step 6 is the relabelling of s.10.4, and that
  is not mere bookkeeping, because the two expansions group infinitely many terms differently. At a fixed hatted order `u^N`, the
  hatted residue is an infinite series in `1/w` collecting original-variable terms `u^(N + (1-4r)/5) z^r`, `r = 4, -1, -6, ...`.
  For example, `s Res_X[e^(-f_0/w) dX] = s Sum_l (-1)^(5l+1) v^(5l+1)/((4l+1)! 4^l l!)` contains `u^(4l+1) z^(-5l-1)` for every `l`.
  So `V = s Res_X(...)` needs a summability argument, as in PB_AUDIT s.2.7. The SKETCH label of s.10.4 is right.

## 3. Statements against their users

### 3.1 The chart proof (`stage5/chart/CHART_PROOF*.md`, which has its own s.2 lemmas)

| chart lemma | used in | note | match |
|---|---|---|---|
| (F1) finiteness, rationality | Theorem Gamma (i), s.3.5 | Lemma U | yes (the chart's bound is sharper) |
| (F2) constants in the phase | everywhere | implicit in Def. 1.2 | trivial |
| **(F3) homogeneity** `N[mu Phi, A(., z)](z) = N[Phi, A(., mu z')](z')`, `z' = z/mu` | s.1.7, s.2.8, s.3.9, s.3.10 | **not stated** (only (H) in L3) | missing, trivial; proof in s.3.2 |
| (S1) integration by parts | (L1a), (L1c) | FL1 | yes |
| (PG) perturbed Gaussian | (Fub), (L1c)(iii') | (G5), (G5') | yes |
| (L0) products | steps (a), (d), (e) | (G4) | yes |
| (L1a) | s.1.7, s.3.9, s.3.10 | L1a (any Q-algebra) | yes |
| (Fub), `beta = (det H_vv(w)/det H_vv(w_0))^(-1/2)` | steps (a), (b), (c) | L1b(4), same normalisation | yes |
| (L1c)(i), (ii) | s.3.5 (rationality), s.3.6 (`L = 1`) | L1c(i), (ii) | yes |
| (L1c)(iii') | (iii) | FL4 | yes |
| (L1c)(iii) | step (e) | Corollary T' = L1c(iii) | yes; the note's statement is the cleaner one |
| (L2) | step (b) | L2(3), (4); `beta = X_0(w)/X_0(w_*)` is `g_c/G_t(tau(w), w)` | yes |
| (L3), all `M`, both measures | steps (a), (d), (e) | L3(i)-(v); `Poch_M(w) St_0 = St_M` by (iv) | yes; the note adds a formal proof (Proof A) |
| reflection of Stirling forms | s.1.7 | not stated | follows from L3(iv), s.3.3 |

The hypotheses listed in the note's s.8 are exactly those the chart proof checks:
* step (a): `5L != 0` and `zeta^j - L != 0` (chart s.3.3);
* step (b): `X_0 = 1 - L != 0` (chart s.3.4);
* step (c): `det H_7 = 1/(L^2 (1 - L^5)^2)` (chart s.3.2);
* step (d): the Hessian `diag(-1/5, -1/(zeta^j - 1))` at `L = 1` (chart s.3.6);
* step (e): negative `M` (chart s.3.7).

My toy chart (part 2, s.7.1) runs this interface through steps (a)-(d) on rational data, including the continuation to `L = c_0`.

The note closes two items that `CHART_PROOF_4` s.5 lists as open:
* item 3, the residue part of (FL) and Stirling's functional equation: the note's s.10.1-s.10.3 and L3(iv). The relabelling of s.10.4
  is not covered by this.
* item 4, a formal proof of (L3): Proof A, s.9.3.

### 3.2 The missing z-rescaling lemma (referee; PROVED)

**Lemma H.** Let `mu` be a unit of `A`, `F` a phase with nondegenerate critical point `c`, and `a(t; z)` an amplitude. Then
`SP^norm_c[mu F; a(.; z)](z) = SP^norm_c[F; a(.; mu z')](z')` at `z' = z/mu`, i.e. after the ring isomorphism `A((z)) -> A((z'))`,
`z -> mu z'`. In full form the symbols match:
`[e^(mu F(c)/z)] (2 pi z)^(n/2) [det(-mu H)^(-1/2)] = [e^(F(c)/z')] (2 pi z')^(n/2) [det(-H)^(-1/2)]`, with the branch
`det(-mu H)^(1/2) := mu^(n/2) det(-H)^(1/2)`.
*Proof.* `Hess(mu F) = mu H`, the covariance is `C/mu`, and the cubic part is `mu F_3`. For a monomial,
`E_(z C/mu)[t^alpha z^j] = mu^(-abs(alpha)/2) <t^alpha>_C z^(j + abs(alpha)/2) = <t^alpha>_C mu^j z'^(j + abs(alpha)/2) = E_(z' C)[t^alpha (mu z')^j]`.
Both sides are A-linear and preserve weight, with finitely many inputs per output power. So
`E_(z C/mu)[P(t, z)] = E_(z' C)[P(t, mu z')]` for every `P` in `R_n(A)`. Take `P = e^(mu F_3/z) a(t; z) = e^(F_3/z') a(t; mu z')`. QED
Users:
* the chart proof: `z = L^-4 w` (s.1.7, s.3.10), `x = a t` (s.2.8), `t = zeta t'` (s.3.9);
* Theorem P(b) step 4: `Y = Y_c eta` (L1a, linear) followed by `w -> tau w` (`mu = Y_c^4`);
* P(a): the identification of the flat with the hatted SP.
(H) of s.9.3 is the case `F = a log x - x`.

### 3.3 The formal reflection identity used by CHART_PROOF s.1.7 (referee; PROVED)

From L3(iv) at `Y = -y`: `Gamma_St(1 - y) = (-y) Gamma_St(-y)` in the formal sense of s.9.1. Then
`log Gamma_St(y) + log Gamma_St(1 - y) = log 2 pi + (y - 1/2)(log y - log(-y)) + mu(y) + mu(-y)`, and `mu(-y) = -mu(y)`, because all
exponents `1 - 2k` are odd. This is the chart's identity.
The branch part is not a formal-SP lemma and is not checked here:
* `log y - log(-y) = i pi eps_j`;
* criticality forces `e = 0`.

### 3.4 Theorem P(b) (POLY s.6, poly referee s.3 and C6)

| P(b) step | formal input | note | status |
|---|---|---|---|
| (a) `PF[SP] = SP[PF ..]` | FL4 with a `q`-independent amplitude, also iterated; `Gamma_St(y+1) = y Gamma_St(y)`; the shift `sigma -> sigma + z` | FL4 (C2), L3(iv), Lemma T(b) | PROVED; VERIFIED on Mellin-Barnes-type data (`rl_fl4`) |
| (a) flat to hatted | L1a (`sigma = L/X`) and z-rescaling | L1a, Lemma H (s.3.2) | PROVED (Lemma H here) |
| 2. reduction mod exact forms | FL1 over `A = Q[[u]]` (or `Q[w][[u]]`) | FL1 | PROVED. `Gamma in Y Q[Y, w][[u]]` is ONE amplitude in `A[[t]][[w]]`, so no u-adic series of exact forms has to be summed |
| 3. CFU invariance | L1a over `Q[[u]]`; residue: FL3 | s.5.4, s.10.3 | PROVED |
| 4. model rows | L1a (`Y = Y_c eta`) and `w -> tau w` | L1a, Lemma H | PROVED (Lemma H here) |
| 6. `V = s Res_X(..)`, original versus hatted variables | summability of two regroupings | s.10.4 | **SKETCH** in the note; written out in PB_AUDIT s.2.7 (not refereed by me) |
| 6-7. residue kills exact forms | FL2 | s.10.2 | PROVED |

So the "four standard formal lemmas" of C6 are not the full list. R5 also uses the Stirling functional equation, the shift `T(b)`,
the z-rescaling (Lemma H) and the relabelling of s.10.4.
*Numbering hazard:* PB_AUDIT s.2.0 numbers its own formal lemmas FL1-FL7. Its FL2 (rank-one quotient, invariance) and FL3
(derivations) are NOT the note's FL2 and FL3.
