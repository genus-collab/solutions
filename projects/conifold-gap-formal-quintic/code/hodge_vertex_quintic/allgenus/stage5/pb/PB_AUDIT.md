# Stage 5 / pb: adversarial audit of Theorem P(b) and of the closing logic `D = 1 <=> (Gamma)` (2026-11-14)

Folder `allgenus/stage5/pb/` (task key `pb`; new files only). AI-written; no human has read it. Labels: PROVED (argument written
here, or in the cited project note and re-checked here line by line), VERIFIED (exact computation; range and field or primes stated),
CONJECTURAL, SKETCH, REFUTED. Notation of the brief: colour 0, `u = 1 - 5^5 q`, `s = u^(1/5) = 1/L`, `z = s^4 w`, `c = (0,-1,3,2,1)`,
`eps = (+,-,-,+,+)`, `b_k = B_2k/(2k(2k-1))`, `h = h_(w_0)`. UR = `natural_form/theory/UNIFORM_REDUCTION.md`, GI =
`allgenus/gauge/GAUGE_IDENTITY.md`. Part 2 (own code, results, findings, verdict): `PB_AUDIT_2.md`.
Read in full: POLY.md, the poly referee's REFEREE.md, EVID.md, MECHANISM_STAGE4.md, STAGE4_SYNTHESIS(_2).md, GI s.0-s.6, UR s.0-s.5,
FORMAL_QUINTIC_genus2.md s.1; the gamma referee's s.2.5-2.6. I wrote my own derivation first and compared afterwards. No worker,
referee or critic code was opened. The only project source I read is the docstring of `fq_conifold.py` (definition of the conifold frame).

## 0. Bottom line

* **Theorem P(b) is PROVED given (U), GI Theorem G and the hypotheses H3-H7 of s.1.** These are the cyclic form (CYC), Prop. S, GI
  Lemmas V, Q and 1, and seven formal lemmas. Each is proved in the project, classical, or proved in s.2. **No mathematical gap was found.**
* Three steps were not written out in POLY s.6 or in the poly referee's s.3. They are written out here, and none changes the statement.
  * **(CYC) is used but not listed.** Lemma 2 compares `Phi^SP` and `Phi^can` as solutions of one equation. `Phi^SP` solves PF5 by P(a);
    `Phi^can` solves the cyclic system. That the two equations are the same is (CYC): classical special geometry, s.2.3, and VERIFIED
    here as an exact operator identity.
  * **Two different expansions of one residue.** `V` is expanded in `Lam`, the hatted residue is not (`e^(-1/(Xw))` stays whole). They
    are rearrangements of one summable family, with an injective monomial map (s.2.7).
  * **The key lemma (s.2.7).** The coefficient of `u^E z^b` in `s Res_X` involves the reduction coefficients `pi_m` only at u-orders
    `<= E + (4b-1)/5`, through finitely many terms. So the hatted residue sum is a well-defined element of `Q((z))[[u]]`. This closes the
    gamma referee's G2 for the poly route.
* **The closing logic is sound, and it is an exact identity:** `D_0(z) = Phi^can(z,1)/h(z) = c(z;1)/h(z)` (s.3).
  * It uses UR s.4, formal WKB uniqueness, (QRR) and the polynomiality of `ph^can_n`. **Lemma N is not needed for the all-orders
    statement**; it is used only by the finite route `g <= 69`.
  * `D_i(z) = D_0(zeta^-i z)` covers all colours.
  * `Phi^can` is frame-independent: it follows from P(b), and is VERIFIED in two frames.
* **Own code (s.4).** It is independent of the poly, gamma, evid and critic code and uses different methods. VERIFIED:
  * `Phi^SP = Phi^can` in every coefficient through `w^80`, mod `2^61-1` and `2^89-1`, frames (0,0) and (1/3,-2/7);
  * `V = C0 (Y^C_Q)_0` exactly over Q through `u^40`, two frames;
  * the mechanism itself through `u^20`: CFU plus reduction reproduces the canonical gauge's first row, and the residue side reproduces
    `V`. The gamma referee noted that this test had never been done;
  * (CYC) exactly through `u^56` in three frames, and (S2) through `t^60`.
  * Every control fails where predicted.
* Nothing is REFUTED. One bug in my own test was found and fixed (part 2, s.5).

## 1. Hypotheses used by P(b): complete list

| # | hypothesis | status | where used |
|---|---|---|---|
| H1 | (U) a regular gauge `P_can` exists (UR s.3) | PROVED in UR, AI-audited; hypotheses re-checked for `(A, At)` (s.2.8); VERIFIED here by an independent construction to `z^80` | step 5 |
| H2 | GI Theorem G for every regular gauge | PROVED (GI s.3); the proof uses only regularity, (R1) and GI Lemma 1, re-checked s.2.8 | steps 5-6 |
| H3 | (CYC) the cyclic system of UR s.1 in the conifold frame is PF5 for `Psi = varpi_0 Y_0` | classical; VERIFIED (UR s.1 to `z^6`; here, exactly, s.4) | step 5 (Lemma 2) |
| H4 | UR Prop. S: unique homogeneous model WKB column with `F_0(0) = 1` | PROVED (pivot `4n != 0`) | step 4, (S2) |
| H5 | GI Lemma V (model residues), Lemma Q, Lemma 1 (quintic vanishing solution) | PROVED (GI s.2), re-checked | steps 6-7 |
| H6 | (FL1)-(FL7) formal SP and residue lemmas, Stirling's functional equation (list s.2.0) | standard; short proofs in s.2.0 | (a), steps 2-6 |
| H7 | POLY Lemma 2: power-series solutions of HPF are `Phi^(0) c(z)`, `c in 1 + z^5 Q[[z^5]]` | PROVED (POLY s.2, poly referee s.2) | step 5 |

Used only by the glue (s.3):
* H8, formal WKB uniqueness (standard). In the matrix form it is UR s.4; in the scalar form, the WKB recursion of PF5 fixes `A_n` up to
  `c L`.
* H9, (QRR) in the form `A^QRR_n(q = 0) = h_n`. **Not used anywhere:** Lemma N, Lemma P, the degree bound `deg ph_n <= floor(4n/5)` (only polynomiality
is used), (Gamma), the chart argument.

## 2. Re-derivation of P(b)

### 2.0 Formal lemmas (H6)
Let `A` be a Q-algebra, here `Q[[u]]` or `Q[w][[u]]`. Let `f in A[[t]]` with `f = f_0 + (a/2) t^2 + O(t^3)` and `a in A^x`, and let
`y = t (2(f - f_0)/(a t^2))^(1/2) in t + t^2 A[[t]]` be its Morse coordinate. Define
`SP[g dt] := sqrt(2 pi w/a) e^(-f_0/w) Sum_m (2m-1)!! (w/a)^m [y^2m] (g(t(y)) t'(y))`.
* **(FL1)** `SP[d(gamma e^(-f/w))] = 0` for `gamma in A[[t]][[w]]`. In `y` the form is `(gamma_y - a y gamma/w) dy`. Its moments are
  `Sum_m (2m+1)!! (w/a)^m gamma_(2m+1) - Sum_m (2m-1)!! (w/a)^(m-1) gamma_(2m-1)`, and the two sums cancel after `m -> m+1`. PROVED.
* **(FL2)** The quotient of amplitudes by exact forms is free of rank one over `A[[w]]`, generated by `1`. Indeed `y^(2m+1) == 0` and
  `y^(2m) == (2m-1)!! (w/a)^m`. So `SP` is the unique `A[[w]]`-linear, `w`-adically continuous functional that kills exact forms and is
  normalised on `1`. A change of variables fixing the critical point (derivative a unit of `A`) maps exact forms to exact forms. Hence
  `SP` is invariant, up to the branch of `sqrt` fixed by `phi'(t_c) in 1 + uA`. PROVED. The Wick (Taylor) definition used by every code
  is the same functional, by the same uniqueness.
* **(FL3)** For a parameter-dependent phase, SP commutes with derivations of the parameter, also when the critical point moves. This
  follows from (FL2): differentiate the identity `SP[...] = coefficient of 1` in the rank-one quotient. Standard; used only in P(a).
* **(FL4), (FL5)** The residue: see s.2.7.
* **(FL6)** `Gamma_St(y+1) = y Gamma_St(y)` as formal series in `1/y`. It is the asymptotic expansion of an exact identity, and
  asymptotic expansions are unique; equivalently, it is the generating identity of the Bernoulli numbers. Standard.
* **(FL7)** `SP[g(sigma) - g(sigma+z)] = 0` for `g = e^(S/z) a`. Proof: `g(sigma+z) - g(sigma) = d/dsigma Sum_(k>=1) z^k g^(k-1)/k!`.
  At each `z`-order the sum is infinite but convergent: e.g. its `z^1` part is `a (e^(S') - 1)/S'`. It is therefore an amplitude with
  coefficients analytic at `L`, and (FL1) applies. POLY s.6 (a) states this step without the remark.

### 2.1 The objects (re-derived from the definitions)
* **Phase.** With `sigma = zs` and `x = 1/sigma`, Stirling for `Gamma(1+5s)` and `1/Gamma(1+s-zeta^j/z)` gives the exponent
  `S(sigma)/z`, `S = sigma log q + 5 sigma log 5sigma - Sum_j (sigma - zeta^j) log(sigma - zeta^j)`. The `log z` terms and the linear terms
  cancel because `Sum_j zeta^j = 0`. Then `S(1/x) = -W`, `W = Lam/x + G(x)`, and
  `G = Sum_j (1/x - zeta^j) log(1 - zeta^j x) = Sum_m x^(5m-1)/(m(5m-1))`: only `5 | n` and `5 | n+1` survive the sum over `j`.
* **Amplitude.** `sqrt(5 sigma/(sigma^5 - 1)) dsigma = -sqrt5 (1-x^5)^(-1/2) dx`, and
  `log at = -(1/2) log(1-x^5) + Sum_k b_k (xz)^(2k-1) [5^(1-2k) - 5 Sum_P C(2k-2+5P, 5P) x^(5P)]`, using
  `Sum_j (1-zeta^j x)^(-n) = 5 Sum_(5|i) C(n-1+i, i) x^i`.
* **Hatted.** `x = sX` gives `W/z = What/w` with `What = Sum_m (u^(m-1)/m)[1/X + X^(5m-1)/(5m-1)]` and `B(X) = at(sX, s^4 w)`. At each
  u-order, `B` is a polynomial in `X` and `w`; `B(1)|_(w=0) = (1-u)^(-1/2)`.
  * `What'(X) = -(uX^2)^(-1) log((1-uX^5)/(1-u))`, so the critical point is `X = 1` for every `u`.
  * `What''(1) = 5/(1-u)` and `What(1) = 5 Sum_m u^(m-1)/(5m-1) = (5/4)/tau`.
* **The two functions.**
  * `Phi^SP := SP_(X=1)[e^(-What/w) B dX] / (sqrt(2 pi w/5) e^(-What(1)/w))`. It has `ph_0 = 1`.
  * `Phi^can := C0 x (colour-0 entry of row 0 of the hatted P_can Q^m(v(u)))`, with `P_11(0,z) = 1` and `C0 = varpi_0`.
  * `P^SP := Q^SP (Q^m)^-1`. Here `Q^SP` is built from `Phi^SP`: rows by the cyclic system, colours by covariance.
* **The two functionals.** `V := Res_(x=0) e^(-W/z) at dx = Sum_(k>=1) (-Lam/z)^k/k! [x^(k-1)](e^(-G/z) at)`, which lies in
  `Q[z^-1][[u]]`. The exponential is `e^(-What(1)/w) = e^(u_0/z)` with `D u_0 = L`, the colour-0 exponential (`dW/dx = 0` at `x = s`).

### 2.2 Theorem P(a): PROVED (given FL3, FL6, FL7)
* `PF5 = M^5 - 1 - q Prod_k (5M + kz) = z^5 D PF4 - 1`. Use `q(5D+5) = 5 D q` and `PF4 = D^4 - 5q Prod_(k<=4)(5D+k)`, the u-form of
  `fq_conifold`.
* `PF5[q^s phi] = g(s) - g(s+1)`, `g = q^s((zs)^5 - 1) phi`: the recursion of `phi` re-checked by hand.
* By (FL6) the same holds for the Stirling form, with the shift realised as `sigma -> sigma + z`. By (FL3) and (FL7),
  `PF5[SP] = SP[PF5 integrand] = 0`.
* On `0 < u < 1` the analytic change of variables `sigma -> X` identifies the flat SP with the hatted SP (FL2). The hatted Taylor data
  (`[t^k] What(1+t)`, `1/What''(1) = (1-u)/5`, `B`) lie in `Q[[u]]`, so `[w^n] Phi^SP in Q[[u]]` solves HPF.
* By POLY Lemma 2 (re-checked: `(4n - 5 theta_u)` is diagonal on `u^k`, and `r_n` is polynomial if the lower `ph` are), these
  coefficients are polynomials.
* HPF's pivot structure, re-derived: `L + zD = s^-1[1 + w(1-u)((4/5) theta_w - theta_u)]`, so the `w^(n+1)` coefficient of HPF is
  `(1-u)(4n - 5 theta_u) ph_n + (1-u) r_n`. This is the only part of HPF that P(b) uses.

### 2.3 (CYC): the hypothesis that makes `Phi^can` a solution of the same equation
`zD Y_a = c_a Y_(a+1)` with `c = (Dt, L^5/(varpi_0^2 Dt^2), Dt, varpi_0, varpi_0)` gives `z^5 (1/c_4) D ... (1/c_0) D Y_0 = Y_0`.
With `Psi = varpi_0 Y_0`, the condition `PF5 Psi = 0` reads `z^5 (1/varpi_0) D PF4 Psi = Y_0`. Since `c_4 = varpi_0` and `D o O = 0` forces
`O = 0`, (CYC) is equivalent to
    `PF4 = (1/varpi_0) D (1/Dt) D (u varpi_0^2 Dt^2) D (1/Dt) D (1/varpi_0)`   (both sides have leading term `u D^4`).
With `theta_t = (1/Dt) D` the right side is `(Dt/varpi_0) theta_t^2 (u varpi_0^2 Dt^3) theta_t^2 (1/varpi_0)`. Its kernel is
`varpi_0 (1, t, F_t, 2F - tF_t)` with `F_ttt = 1/(u varpi_0^2 Dt^3)`, i.e. the special-geometry form, with Yukawa coupling `Y_DDD` (up to the
factor 5 of FORMAL_QUINTIC_genus2 s.1) `= 1/u`. This is classical for a fourth-order Calabi-Yau operator in a frame `(varpi_0, varpi_0 t)` with
vanishing symplectic pairing. Every conifold frame has this property: `B(analytic, g) = 0`, as the docstring of `fq_conifold` notes.
Status: classical, used and not listed in POLY s.6; **VERIFIED here exactly** (s.4). The control "`t = f_2/f_0`, a non-isotropic pair"
fails, as it must.

### 2.4 CFU (step 1): PROVED
Put `Y = X y(T)` with `T = X^5`; then `V/Y + Y^4/4 = X^-1 (V/y + T y^4/4)`. The equation `V/y + T y^4/4 = Sum_m (u^(m-1)/m)(1 + T^m/(5m-1))`
at order `u^k` reads `(T-1) Y_k + V_k = R_k(T)`, with `R_k` polynomial in `T` and built from lower orders. So `V_k = R_k(1)` and `Y_k` is a
polynomial; the solution is unique. `X = 1` goes to a critical point of `f = V/Y + Y^4/4`, because `Y'(1) in 1 + uQ[[u]]`. Therefore
`What(1) = (5/4) V^(4/5)` and `V = tau^(-5/4)`. The inverse `X(Y) = Y x(Y^5)` is also polynomial at each u-order. VERIFIED (s.4):
`V = tau^(-5/4)` through `u^20`, and `Y_1 = (1+T)/18`.

### 2.5 Reduction modulo exact forms (step 2): PROVED
* `beta(Y) dY := B(X(Y)) X'(Y) dY` with `beta in Q[Y, w][[u]]`. The identity
  `w d(Y^m e^(-f/w)) = (m w Y^(m-1) + V Y^(m-2) - Y^(m+3)) e^(-f/w) dY` gives `Y^n == (n-3) w Y^(n-4) + V Y^(n-5)` for `n >= 4`.
* **Termination.** Each step lowers the degree by at least 4. Starting from a polynomial, it ends in the exponents `{-1, ..., 3}`; the
  lowest exponent ever produced is `4 - 5 = -1`. The result is `beta e^(-f/w) dY = Sum_(m=-1..3) pi_m Y^m e^(-f/w) dY + d(Gamma e^(-f/w))`,
  with `Gamma in Y Q[Y, w][[u]]` and `pi_m in Q[w][[u]]`. Uniqueness of the `pi_m` is not needed: s.2.6 pins them down.
* **Z/5 grading.** `f(zeta Y) = f(Y)/zeta`, `B(zeta X; w) = B(X; zeta w)` (since `at` depends on `xz` and `x^5`), and CFU commutes with
  `X -> zeta X`. With `deg Y = 1`, `deg w = -1`, the form `beta dY` has degree 1 and the relation preserves degree mod 5. So `pi_m` contains
  only `w^k` with `k = m mod 5`; for `m = c_b` that is `k = -b mod 5`, the grading of the hatted entry `G_b`. VERIFIED (s.4).

### 2.6 Identification (step 4): signs, weights, `kappa = C0`. PROVED
* **(S2)** `SP_(eta=1)[eta^(c_b) e^(-(1/eta+eta^4/4)/t) d eta] / (sqrt(2 pi t/5) e^(-5/(4t))) = eps_b F_b(t)`. **PROVED**, not only VERIFIED
  (the poly referee checked it through `w^50`).
  * The three-line argument of UR s.1 uses only `z d/dv I_m = -I_(m-1)` and `I_(m+3) = m z I_(m-1) + v I_(m-2)`, i.e. (FL3) and (FL1).
    So `(eps_b SP[x^(c_b) e^(lambda/z) dx])_b` solves the model system.
  * It is homogeneous (`x -> kx`, `v -> k^5 v`, `z -> k^4 z`). By Prop. S it is a constant multiple of `(F_b)`, and the constant is 1.
  * VERIFIED through `t^60` at two primes (s.4).
* **Weights.** With `Y = Y_c eta` and `Y_c = V^(1/5) = tau^(-1/4)`: `f(Y_c eta) = Y_c^4 (1/eta + eta^4/4)` and `w/Y_c^4 = tau w`. So
  `SP_(Y_c)[Y^(c_b) e^(-f/w) dY] = N_w eps_b tau^((1-c_b)/4) F_b(tau w)`, with `N_w = sqrt(2 pi w/5) e^(-f(Y_c)/w)`,
  `Y_c^(c_b+1) tau^(1/2) = tau^((1-c_b)/4)` and `f(Y_c) = What(1)`. By (FL2) for CFU and (FL1) for `d(Gamma e^(-f/w))`:
  `Phi^SP = Sum_b eps_b pi_(c_b) tau^((1-c_b)/4) F_b(tau w)`.
* **Gauge side.** `Q^m_b = sigma^(c_b-1) F_b(z/sigma^4)` with `sigma = s tau^(-1/4)`, and `Qhat_a = s^(1-c_a) Q_a`. So
  `Ghat_ab = s^(c_b-c_a) tau^((1-c_b)/4) P_ab`, and row 0 reads `Phi/C0 = Sum_b tau^((1-c_b)/4) G_b F_b(tau w)` with
  `G_b = s^(c_b) P_0b in Q[[u]][[w]]`, graded (only `w^k` with `k = -b`).
* **Triangular uniqueness in `Q[[u]][[w]]`.** At `w^N` the new term is `[w^N] g_(-N)` times a unit times `eps_(-N) != 0`. Both
  `eps_b pi_(c_b)/C0` (by s.2.5) and `G^SP_b` (since `P^SP = s P_can` is graded, s.2.8) are graded, so `G^SP_b = eps_b pi_(c_b)/C0`.
  Only colour 0 is used; the other colours follow by covariance.
* VERIFIED directly (s.4): `eps_b pi_(c_b) = C0 G^can_b` for all `b` through `u^20`, `w^26` (mod `2^61-1`).

### 2.7 The residue (step 6): PROVED, with the bookkeeping written out
* **(FL4)** At each u-order the hatted integrand is `e^(-(1/X + X^4/4)/w)` times a Laurent polynomial in `X` over `Q[w, 1/w]`.
  * Let `R` be the ring of doubly infinite Laurent series in `X` whose coefficients are `1/w`-adically convergent elements of `Q((1/w))`.
    `e^(-V/(Xw)) e^(-X^4/(4w))` lies in `R`: the `w`-powers of the `X^n` coefficient tend to `-oo`.
  * `R` is stable under Laurent polynomials in `X` and under `d/dX`, and `Res := [X^-1]` kills `d/dX(R)`.
  * Equivalently, at fixed `w != 0`, `(1/2 pi i)` times a contour integral on `C^*` (poly referee).
* **(FL5)** Write `Y = X + u h` with `h in Q[X][[u]]` (CFU). Then `F(X+uh)(1+uh') = F + Sum_(j>=1) (u^j/j!) (F^(j-1) h^j)'`, re-checked by
  hand. At each u-order the sum is finite, so `Res_X` of the pull-back equals `Res_Y`. The Taylor shift is a ring homomorphism, so the
  pull-back of `e^(-f/w) beta dY` is the hatted integrand computed directly.
* **Model residues.** `Res_(Y=0)[Y^m e^(-f/w) dY] = s^-(m+1) I^C_m(v(u), z)` termwise: put `x = sY`; then `f/w = (x^4/4 + v/x)/z` with
  `v = uV`, the `v(u)` of UR.
* **Two expansions of one residue.** `V` expands `e^(-Lam/(xz))` in `Lam`; the hatted residue keeps `e^(-1/(Xw))` whole. Both are
  rearrangements of ONE family of monomials.
  * With `x = sX`, `z = s^4 w`, the term `(-Lam/(xz))^k/k! . c_(a,b) x^a z^b`, with `u^(k+r)` taken from `Lam^k`, becomes the monomial
    `u^(r + (a+4b)/5) w^(b-k) X^(a-k)`.
  * For a fixed `(u, w)`-degree only finitely many `(k, r, b)` occur, and `u^N w^beta -> u^(N - 4beta/5) z^beta` is injective on monomials.
  * Hence `V = s Res_X[e^(-What/w) B dX]` as an identity of summable families. PROVED.
* **Key lemma (new).** Consider `s Sum_m pi_m Res_Y[Y^m e^(-f/w) dY] = Sum_b pi_(c_b) s^(-c_b) I^C_(c_b)(v(u), z)`.
  * Combine a term `u^i w^j` of `pi_m` (`j = m mod 5`), the term `(-1/(4z))^(i')(-v/z)^e/(i'! e!)` of `I^C_m` (`e = m + 4i' + 1`) and
    `[u^(e+l)] v^e`. The product lands on `u^E z^beta` with `beta = j - 5i' - m - 1` and `E = i + l - (4beta - 1)/5`.
  * So `(u^E, z^beta)` receives `pi` only from u-orders `i <= E + (4beta-1)/5`, and from finitely many `(j, i')` each (`pi_m` at `u^i` is
    polynomial in `w`).
  * The sum is therefore a well-defined element of `Q((z))[[u]]`, and only `beta = 4 mod 5` occurs. PROVED.
* **Hence** `V = Sum_b eps_b pi_(c_b) s^(-c_b) Y^C_b = C0 Sum_b P^SP_0b Y^C_b = C0 (P^SP Y^C)_0`. This uses `Y^C_b = eps_b I^C_(c_b)`
  (GI Lemma V) and `s^(-c_b) G_b = P_0b`. So `kappa = C0` exactly, in every frame, as the poly referee found.
  * VERIFIED (s.4) from my own `pi_m` through `u^20`, on every admissible key. That includes 45 keys `z^beta`, `beta > 0`, where both sides
    must vanish.
* **Remark.** In step 7 (`E = 1`) the vanishing of `z^(5q-1)` at `u^1` is `s_q = 0`, and it involves the reduction up to u-order `4q`. So
  step 7 is one finite identity per `q`; the structure proves all of them at once.

### 2.8 Regularity, Theorem U and Theorem G (step 5): PROVED given H1, H2, H3
* Both `Phi^SP` and `Phi^can` solve HPF (s.2.2; s.2.3 with H3) and have `ph_0 = 1`. For `Phi^can`: `p_(0,0) = tau^(-1/4)/varpi_0`
  (the diagonal of GI s.5) gives `Ghat_00(u,0) = 1/C0`. By Lemma 2, `Phi^SP = Phi^can s(z)` with `s in 1 + z^5 Q[[z^5]]`
  (`z^5 = u^4 w^5`).
* Rows by the cyclic system and colours by covariance commute with `s(z)`, since `s(zeta z) = s(z)`. So `Q^SP = s Q^can` and
  `P^SP = s P_can`, which is regular and graded.
* **The hypotheses of U for `(A, At)` hold.**
  * `b_0 = b_2 = t_c'` and `b_3 = b_4 = -varpi_0/(1-u)` are units: `t_c'(0) = varpi_0(0) = 1` in every frame.
  * `u b_1` is a unit.
  * `bt = (v', -v'/v, v', -v', -v')`.
  * `Prod b = -1/(u(1-u)^5) = Prod bt` iff `v^(-1/5) dv = u^(-1/5) du/(1-u)` iff `v^(4/5) = u^(4/5)/tau`.
* **Theorem G for `P^SP`.** Its proof needs four things, and all hold for `P^SP`:
  * `Y^C(v(u), z)` solves the pulled-back model (chain rule);
  * `P Y^C in Q((z))[[u]]^5` solves the quintic system;
  * (R1): order `u^-1` with `r = rt = -1`;
  * GI Lemma 1, since `det((n+1)z - B_-1) = ((n+1)z)^5` is a unit of `Q((z))`.

  With `P^SP_11(0,z) = s(z)` this gives `(P^SP Y^C)_0 = s(z) (Y^C_Q)_0`.

### 2.9 Step 7: PROVED
`V = C0 s(z) (Y^C_Q)_0` in `Q((z))[[u]]`.
* At `u^1`, `V_1 = -1/z`: only `k = 1` contributes, and `at(0,z) = 1`, `G(0) = 0`.
* `[u^1] (Y^C_Q)_0 = -1/z`: by GI Lemma Q(a), `ord_u y^(n)_0 = (4n+1)/5`, so only `n = 1` contributes.
* `C0(0) = 1`.

Hence `-1/z = -s(z)/z` in `Q((z))`, so `s = 1`. Alternative ending: the gamma referee's Theorem A(i), `V = C0 (Y^C_Q)_0` without U,
gives `s(z) (Y^C_Q)_0 = (Y^C_Q)_0` directly. QED

## 3. The glue: `D = 1 <=> (Gamma)`, PROVED given P(b), H8, H9
* **UR s.4.** `Q D = P_can Q^m(v(u))`, with `D` constant and diagonal. Its row 0, colour 0, hatted, reads
  `Phi^QRR(w,u) D_0(s^4 w) = Phi^can(w,u)`: the factor `1/C0` is common to both sides.
  * `Phi` is the WKB amplitude of the I-function and is frame-free. `C0` is never evaluated at `q = 0`, where it is singular.
* **The scalar identity is settled near `q = 0`, with no convergence of the gauge needed.**
  * The `ph^can_n` are polynomials (Lemma 2), so `A^can_n = u^(-(4n+1)/5) ph^can_n` is analytic near `q = 0`, like `A^QRR_n`.
  * Both solve the WKB recursion of PF5. That recursion fixes `A_n` up to `c L`, i.e. up to `Phi -> Phi (1 + c z^n)` (H8). So the two
    differ by one constant series `d(z)`, and the colour-0 entry of UR s.4 identifies `d = D_0`.
  * At `u = 1` (`s = 1`, `z = w`), with (QRR) `A^QRR_n(0) = h_n`, this gives `h(z) D_0(z) = Phi^can(z,1)`.
* **With P(b): `D_0(z) = Phi^SP(z,1)/h(z) = c(z;1)/h(z)` EXACTLY.** This is gamma's Corollary A at `L = 1`.
  * By P(a), `c_n(L) = L^(4n) ph^SP_n(L^-5)` is a Laurent polynomial. So "regular at `L = 1`" in (Gamma) is automatic, and the content of
    (Gamma) is the value.
* **All colours.** `D_i(z) = D_0(zeta^-i z)`: `h_(w_i)(z) = h_(w_0)(zeta^-i z)` because `w_i = zeta^i w_0`. So `D_0 = 1` gives
  `D = 1`, and (Gamma) gives (F_g) for every `g`.
* **Same object.** EVID's `Phi^can` (`K_q = 0`) and POLY's (`P_11(0,z) = 1`) coincide, because Lemma N says `K_q = [z^(5q)] P_11(0,z)`.
  VERIFIED: my `Phi^can`, built by the Theorem-U construction, reproduces EVID's `f_5`, `f_10`, `f_15` and `m5` at both primes and in both
  frames (`pb_glue`).
* **Where Lemma N is needed.** Its formula for `K_q` serves only to compute `D_0` order by order (the finite route `g <= 69`). The
  all-orders chain needs only the identity above.
