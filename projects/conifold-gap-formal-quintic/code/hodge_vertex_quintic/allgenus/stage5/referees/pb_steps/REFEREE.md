# Referee report: Stage 5 / pb, PB_AUDIT s.2.6-s.2.9 and s.3 (2026-11-14)

Referee of the unrefereed steps of the second audit of Theorem P(b): `stage5/pb/PB_AUDIT.md` (PB) s.2.6 (identification: signs
`eps_b`, weights, `kappa = C0`, (S2)), s.2.7 (the residue: FL4, FL5, the relabelling `V = s Res_X`, the key lemma), s.2.8-s.2.9
(regularity, Theorem U and Theorem G for `P^SP`, step 7) and s.3 (the glue `D_0 = Phi^can(z,1)/h(z) = c(z;1)/h(z)`). AI-written;
no human has read it. Labels: PROVED, VERIFIED, CONJECTURAL, SKETCH, REFUTED, as in the brief. Part 2: `REFEREE_2.md` (own code and
results, the precise list of corrections, what survives, record).
Read in full: PB_AUDIT.md, PB_AUDIT_2.md; STAGE5_SYNTHESIS.md and `_2` (s.1, s.2, K1-K14, s.6, s.9); POLY.md; the poly referee's
REFEREE.md; FORMAL_SP_LEMMAS_2.md s.4-s.5 and `_4.md` s.10 and `_5.md`; the lemma referee's REFEREE.md and REFEREE_2.md s.5;
UNIFORM_REDUCTION.md (UR) s.0-s.4; GAUGE_IDENTITY.md (GI) s.0-s.5; EVID.md s.2. The pb logs were read as data only. No pb, poly,
evid, gamma, lemma or critic code was opened or imported. Abbreviations: FSL = FORMAL_SP_LEMMAS (+ parts), RL = lemma referee,
Lemma H = RL s.3.2 (z-rescaling), `L` = FSL's ring of s.10.1. Notation of the brief and of PB.

## 0. Verdict: sound, with corrections

* **No mathematical error found in s.2.6-s.2.9 or s.3.** Every displayed identity was re-derived by hand. Every step that admits a
  computation was re-done by new exact code over Q with failing controls (part 2, s.5).
* **s.2.6 is PROVED as claimed.** That includes (S2): the model system (FL4 with a moving critical point, and FL1 used once, at
  `m = 0`), homogeneity (L1a plus Lemma H) and the uniqueness half of UR Prop. S (pivot `4n`). The signs `eps_b`, the weights
  `tau^((1-c_b)/4)` and `kappa = C0` are right.
  * VERIFIED, gauge-free: `Phi^SP = Sum_b eps_b pi_(c_b) tau^((1-c_b)/4) F_b(tau w)`, exact, through `w^25`, `u^25`.
  * VERIFIED: (S2) and the pivots through `t^80`.
* **s.2.7: the relabelling is a complete argument.** It needs two sentences that PB leaves implicit (s.2.2 below): the precise
  family, and why the hatted residue computed in `L[[u]]` is that family's sum in the hatted order.
  * The answers to the brief's three questions are all yes. The fibres are finite in both orders. The regrouping
    `e^(-Lam/(xz)) = e^(-1/(Xw)) e^(-(Lam/u-1)/(Xw))` holds coefficient by coefficient, so in `L[[u]]` and in FSL's ring `B`. The
    monomial map is injective.
  * VERIFIED: the hatted residue computed directly (with `e^(-f0/w)` kept whole) equals the relabelled `V`, exact, hatted orders
    `u^0..u^28`, every w-power that `V` through `u^120` reaches (965 coefficients).
* **The key lemma's bound `i <= E + (4b-1)/5` is right, and it is attained.** The next sentence, "therefore a well-defined element
  of `Q((z))[[u]]`", needs one more fact, which PB does not state. I prove it in s.2.3: a weight bound `m + 4j <= 5i` for the
  monomials `u^i w^j` of `pi_m`. It is not needed by the chain.
* **s.2.8-s.2.9 are correct.** One small redundancy and one implicit fact are noted (s.3).
* **s.3 (the glue) is correct.** It uses formal WKB uniqueness and (QRR) in the form `A^QRR_n(0) = h_n` correctly.
  * One sentence is missing: the ring in which the near-`q = 0` identity and UR s.4 are compared. That ring is `Q[L, 1/L][[z]]`,
    inside `Q((s))[[w]]` (s.4.2).
  * EVID's `Phi^can` (`K_q = 0` for all `q`) and POLY's (`P_11(0,z) = 1`) are literally the same element of `Q[u][[w]]`. The reason
    is Lemma N: `Ghat(Phi) = hatted(P_can c)` with `K_q = [z^(5q)] c`.
* **Corrections:** C1-C10 (part 2, s.6). All concern justification and wording; none changes a statement or a label of the chain.
* **Overall:** I agree with the critic's K5 and K12. K5's point (iv) needs the weight bound of s.2.3.

## 1. PB s.2.6: identification (signs, weights, `kappa = C0`, (S2))

**(S2): PROVED (confirmed).** I re-did the argument with these details.
* (i) *The data.* Put `A = Q[[v-1]]`, `sigma = v^(1/5)` in `A`, `lambda = -(x^4/4 + v/x)`, and `J_m := SP_(x=sigma)[x^m e^(lambda/z) dx]`
  (full SP; Hessian `-lambda''(sigma) = 5 sigma^2`, a unit).
* (ii) *The model system.* FSL FL4 (s.4.5, full-SP form, moving critical point) gives `z dJ_m/dv = -J_(m-1)`. FL1 (FSL s.3.1) for
  `d(x^0 e^(lambda/z))` gives `J_3 = v J_(-2)`; this is the only exact-form relation the system uses. Then
  `Y := (J_0, -J_(-1), -J_3, J_2, J_1)` solves `z d/dv Y_b = bmod_b Y_(b+1)`, `bmod = (1, -1/v, 1, -1, -1)`. Re-derived row by row:
  * `b = 1`: `z d/dv(-J_(-1)) = J_(-2) = J_3/v = (-1/v) Y_2`;
  * the other four rows are the derivative identity alone.
* (iii) *Homogeneity.* L1a over `A` with the linear map `x = sigma x'`, then Lemma H with `mu = sigma^4` (a unit of `A`), give
  `J_(c_b) = sqrt(2 pi z/5) e^(U/z) sigma^(c_b - 1) G_b(z/sigma^4)`. Here `U = -(5/4) sigma^4`, and `G_b(t)` is the normalised SP
  of `eta^(c_b)` at `v = 1`. So `Y_b = eps_b J_(c_b)` has the form of Prop. S with `F_b = eps_b G_b`.
* (iv) *Prop. S.* The chain rule gives `z d/dv[e^(U/z) sigma^(c-1) F(t)] = e^(U/z) sigma^(c-2) T_(c-1) F` (using `U' = -1/sigma`,
  `sigma' = 1/(5 sigma^4)`). So the system is exactly UR Prop. S's T-system:
  `F_1 = T_-1 F_0, F_2 = -T_-2 F_1, F_3 = T_2 F_2, F_4 = -T_1 F_3, F_0 = -T_0 F_4`.
* (v) *Uniqueness.* `F_0` enters the closing relation at `t^(n+1)` with coefficient `-(1/5) Sum_i (k_i - 4n) = 4n`, since
  `Sum k_i = -1 - 2 + 2 + 1 + 0 = 0`. So the solution with `F_0(0) = 1` is unique.
  * This half of Prop. S does not use (H). UR's label "PROVED given (H)" concerns the identification of the symbol, not the
    uniqueness. PB's H4 uses it correctly.
* (vi) *The constant.* `G_0(0) = 1`, so `eps_b G_b = F_b`.
* VERIFIED (part 2, T4): pivots `4n` for `n <= 80`. The SP rows, computed by Lagrange-Buermann, satisfy the T-system and equal `F_b`
  through `t^80`. Controls fail: all signs `+`; `c_b -> c_b + 1`; `T_k -> T_(k+1)`.
* *Minor (C1).* PB s.2.0 says its FL3 is "used only in P(a)". It is also used here.

**Weights: confirmed.**
* The algebra: `f(Y_c eta) = Y_c^4 (1/eta + eta^4/4)` because `V = Y_c^5`; `w/Y_c^4 = tau w`; `f(Y_c) = (5/4)/tau = What(1)`; and
  `Y_c^(c_b+1) tau^(1/2) = tau^((1-c_b)/4)`.
* The inputs: L1a (`Y = Y_c eta`, linear over `Q[[u]]`) and Lemma H with `mu = Y_c^4 = 1/tau`. PB cites neither; RL C1 already
  asked for Lemma H.
* `Phi^SP = Sum_b eps_b pi_(c_b) tau^((1-c_b)/4) F_b(tau w)` also uses three further steps:
  * FSL s.5.4: the SP at `X = 1` equals the SP at `Y_c`, with `Y(1) = V^(1/5)` and `Y'(1)` in `1 + uQ[[u]]`;
  * FL1 over `A = Q[[u]]` for `d(Gamma e^(-f/w))`, `Gamma` in `Y Q[Y,w][[u]]`. This is ONE amplitude in `A[[t]][[w]]` (RL s.3.4);
  * `A[[w]]`-linearity, since `pi_m` lies in `Q[w][[u]]`, inside `Q[[u]][[w]]`.
* VERIFIED independently of any gauge (part 2, T5): exact through `w^25`, `u^25`. Controls: without `eps_b` it first fails at
  `(w^3, u^3)`; without the tau-weights at `(w^0, u^1)`; with `F_b(w)` at `(w^1, u^1)`.

**Gauge side: confirmed.**
* `Qhat_a = Sum_b [s^(c_b - c_a) P_ab] tau^((1-c_b)/4) F_b(tau w)`. So row 0 is `Phi/C0 = Sum_b tau^((1-c_b)/4) G_b F_b(tau w)`,
  with `G_b = s^(c_b) P_0b`.
* Grading: `c_b = -b mod 5`, so `G_b` is in `Q[[u]][[w]]` and contains only `w^k` with `k = -b mod 5`. Regularity of `P` makes
  `[w^k] G_b` divisible by `u^((c_b + 4k)/5)`.

**Triangular uniqueness: confirmed.** Let `Delta_b` be the difference of two graded families with
`Sum_b tau^((1-c_b)/4) Delta_b F_b(tau w) = 0`. At its least w-degree `N`, only `b = -N mod 5` contributes, through
`tau^(..) eps_b [w^N] Delta_b`; so that coefficient vanishes, a contradiction.
* No circularity: s.2.6 uses s.2.8 (`P^SP = s P_can` graded), and s.2.8 uses neither s.2.6 nor s.2.7.
* The direct check `eps_b pi_(c_b) = C0 G^can_b` (pb_mech, mod `2^61-1`) needs the gauge. I did not re-run it. T5 covers the half
  of s.2.6 that needs no gauge.

**`kappa = C0`: confirmed.** In `V = s Sum_b pi_(c_b) Res_Y[...]`:
* the signs cancel: the same `eps_b` occurs in (S2) and in GI Lemma V (`Y^C_b = eps_b I^C_(c_b)`), and `eps_b^2 = 1`;
* the powers of `s` cancel: `s . s^(-(c_b+1)) . s^(c_b) = 1`;
* no power of `tau` survives: the residue side carries no tau-weight, and `G_b` is defined without it.

## 2. PB s.2.7: the residue

### 2.1 FL4, FL5, model residues
* **(FL4) is correct, but the ring is described imprecisely (C2).** "Doubly infinite Laurent series in `X` whose coefficients are
  1/w-adically convergent elements of `Q((1/w))`" should state the actual condition: `ord_(1/w) c_n -> +oo` as `|n| -> oo`.
  * With it, products converge and the set is a ring, by the argument of FSL Lemma 10.1.
  * Simpler: use FSL's ring `L` (linear growth; refereed, RL s.1). All integrands here lie in `L[[u]]`.
* **(FL5) is FSL FL3 (s.10.3, refereed).** The telescoping identity `F(X+uh)(1+uh') = F + Sum_j (u^j/j!)(F^(j-1) h^j)'` was
  re-derived. `h` lies in `Q[X][[u]]` because the CFU coefficients `y_k` are polynomials, of degree `<= k` in `T`.
* **Model residues: an exact termwise identity.** Put `v = s^5 V`, `z = s^4 w`. The `i`-th term of `s^(-(m+1)) I^C_m(v, z)` equals the
  `l = i` term of `Res_Y[Y^m e^(-f/w)]`, which is `(-V/w)^e/e! (-1/(4w))^i/i!` with `e = m + 4i + 1`.
  * Expanding `V^e` in `u` is again a relabelling. It is covered by s.2.2.

### 2.2 The relabelling `V = s Res_X[e^(-What/w) B dX]`: PROVED (PB's argument, completed)
* *The labels.* Monomials are `u^E z^b x^al` (original) and `u^N w^j X^al` (hatted). The substitution `x = sX`, `z = s^4 w`, `u = s^5`
  acts on exponents by `iota(E, b, al) = (E + (4b + al)/5, b, al)`, which is injective.
  * With `dx = s dX`, the claim is coefficientwise: `[u^N w^j] H = [u^(N + (1-4j)/5) z^j] V`, where `H := Res_X[e^(-What/w) B dX]`.
* *The family F.* Write `A := e^(-G/z) at = exp(l)`. The exponent `l` is a sum of atoms:
  * `-x^(5m-1)/(m(5m-1) z)`;
  * `x^(5n)/(2n)`;
  * `b_k c_(k,P) (xz)^(2k-1) x^(5P)`.

  A member of F is a triple (`k`, `r`, a finite multiset of atoms). Its term is
  `[u^r](Lam/u)^k (-1)^k/k! u^(k+r) (xz)^(-k)` times the product of the atoms, divided by the multiplicity factorials.
* *Original order.* A fixed `(u^E, z^b, x^-1)` receives finitely many members. Indeed `k <= E`, and every atom has positive
  x-degree, so the multiset is bounded by the x-degree `k - 1`.
* *Hatted order.* Every atom has hatted weight 0: its hatted u-order is `(a + 4b)/5 >= 0`. The only atom of hatted u-order 0 is
  `-x^4/(4z) -> -X^4/(4w)`, and `u^(k+r)(xz)^(-k) -> u^r (Xw)^(-k)`.
  * Fix `(u^N, w^j, X^-1)`. Then `r <= N`, and at most `N` atoms have positive hatted order, from a finite list.
  * The X- and w-balances `-k + 4l + a' = -1` and `-k - l + b' = j` then fix `k` and the multiplicity `l` of `X^4/(4w)` uniquely
    (determinant 5). So this fibre is finite too.
* *`H` is the hatted-order sum of F.* In `L[[u]]` (FSL s.10.1),
  `e^(-What/w) B = e^(-v/X) e^(-vX^4/4) . exp(O(u) Laurent polynomials)`, and `L`'s product is the convergent sum of products of
  coefficients.
  * `e^(-v/X)` is matched with the `u^0` part of `(Lam/u)^k` by the coefficientwise identity
    `e^(-v(Lam/u)/X) = e^(-v/X) e^(-v(Lam/u - 1)/X)`. At `(u^N, X^-n)` both sides equal `(-v)^n/n! [u^N](Lam/u)^n`, a finite
    expression. So the regrouping holds in `L[[u]]`, and also in FSL's ring `B`.
* *Conclusion.* `V` and `s H` are sums of one family over the fibres of two injective labellings related by `iota`. The fibres are
  finite in both orders, so the identity holds coefficient by coefficient. QED.
  * PB states the hatted finiteness for the coarse index `(k, r, b)`, taking the coefficients of `A` as atoms. That is legitimate:
    `iota` is injective on the monomials of `A`, and the same argument applied to `l` shows that the hatted coefficients of `A` are
    those of `e^(-Ghat/w) B`.
  * RL's objection (infinitely many terms grouped differently) is answered by finite fibres. For example, hatted `u^0` collects
    the original orders `u^1, u^5, u^9, ...`, which is RL's `u^(4l+1) z^(-5l-1)`.
* VERIFIED (T1): `H` computed directly equals the relabelled `V`, exact, on hatted orders `u^0..u^28` and all w-powers reached by
  `V` through `u^120`. That is 965 coefficients, 779 non-zero; `H` has no positive w-power and no w-power outside `4 mod 5`.
  Five controls fail, each at the first compared coefficient: `Lam/u -> 1` (no regrouping factor), no `(1-uX^5)^(-1/2)`,
  `b_k -> -b_k`, an extra `u` per Stirling term, `G` cut to `x^4/4`.

### 2.3 The key lemma: the bound is PROVED; the sentence on `Q((z))[[u]]` needs the weight lemma
* *The bound (re-derived).* Combine a term `u^i w^j` of `pi_m`, the `i'`-th term of `I^C_m`, and `[u^l] V^e` (`v = uV`). The product
  has z-power `b = j - 5i' - m - 1`, which is `4 mod 5`, and u-power `E = i + l + e - (4j + m)/5 = i + l - (4b - 1)/5`. So only
  `i <= E + (4b - 1)/5` occurs, with finitely many `(j, i')` for each `i`. Correct.
* *The bound is attained.* VERIFIED (T3): the top order `i = E + (4b-1)/5` contributes on 744 of the 884 admissible keys. Truncating
  `pi` at `u^24` makes 30 of the 32 keys of the layer `E + (4b-1)/5 = 25` wrong; with `u^25` all 32 are right.
* *What finiteness does not give.* PB writes "the sum is therefore a well-defined element of `Q((z))[[u]]`". Finite fibres give a
  well-defined formal sum only.
  * For fixed `E`, the bound `b >= (1 - 5E)/4` does follow from `i + l >= 0`.
  * But `E >= 0` needs a bound on the w-degree of `pi_m`, and PB gives none.
* **Weight lemma (PROVED here).** Give `u^N Y^n w^j` the weight `omega = 5N - n - 4j`.
  * (a) Every monomial `u^N X^p w^q` of `B` has `5N - p - 4q = 0`, since `x^a z^b -> u^((a+4b)/5) X^a w^b`.
  * The CFU map is `X(1 + terms u^k T^d)` with `d <= k`, by the CFU induction (`deg R_k <= k + 1`).
  * In the original variable it is therefore `x(1 + Sum c u^(k-d) x^(5d))`: a change of variables over `Q[[u]]`, tangent to the
    identity. Its inverse has the same form, so it is `Y(1 + terms u^k T^d)` with `d <= k` as well.
  * So every monomial of `beta` and of `X'(Y)` has `omega >= 0`.
  * (b) The reduction keeps `omega >= 0`. The step `Y^n -> (n-3) w Y^(n-4)` keeps `omega`; the step `Y^n -> V_l u^l Y^(n-5)`
    raises it by `5(l + 1)`.
  * Hence `5i - m - 4j >= 0` for every `u^i w^j` in `pi_m`, and `E = (5i - 4j - m)/5 + e + l >= 0`. QED.
  * Via s.2.6 this is the same statement as the regularity of row 0 of `P^SP`; here it comes from the mechanism alone. The chain
    does not need it: it compares two elements of `Q((z))[[u]]`.
  * VERIFIED (T2, `u^0..u^25`): the a-priori bound `floor((5i+1)/4)` is reached at every order except `i = 4q + 3`.
    * By s.2.6, the top coefficient there is `-C0(0) p_(5q+4,1)(0)`. It vanishes, as (R1) of GI s.3 requires.
    * So the mechanism reproduces (R1) by itself.
* *"Hence `V = Sum_b eps_b pi_(c_b) s^(-c_b) Y^C_b = C0 (P^SP Y^C)_0`": confirmed.*
  * The hatted product `pi_m . Res_Y[...]` has finite fibres, because `pi_m` is polynomial in `w` at each u-order. `iota` is a
    monoid map, so the fibres correspond and the image of the product is the product of the images.
  * Then use s.2.6 (`pi_(c_b) = C0 eps_b s^(c_b) P^SP_0b`) and GI Lemma V.
  * VERIFIED in hatted variables (T2): `H = Sum_m pi_m Res_Y[Y^m e^(-f/w) dY]` at hatted `u^0..u^25`, 4283 coefficients. Controls
    fail: no `X'(Y)`; `V -> 1` in the reduction; no CFU.
  * VERIFIED in original variables (T3): 884 admissible keys, including 200 with `b > 0` and 134 with `E <= 0`, where both sides
    vanish. Controls fail: sign of `pi_-1`; `pi_2` and `pi_3` exchanged; `V -> 1`.
* *Remark (step 7).* `E = 1`, `b = 5q - 1` gives `E + (4b-1)/5 = 4q`. Correct.

## 3. PB s.2.8-s.2.9

* **s.2.8: correct.**
  * `ph^can_0 = 1`: GI s.5 gives `p_(0,0) = tau^(-1/4)/varpi_0`, hence `Ghat_00(u,0) = 1/C0`.
  * Both functions have `[w^n]` in `Q[[u]]` and solve HPF: P(a) for `Phi^SP`; (CYC) for `Phi^can`. So Lemma 2 gives
    `Phi^SP = Phi^can s(z)`, `s` in `1 + z^5 Q[[z^5]]`.
  * `s(zeta z) = s(z)` and `th` commutes with `s(z)`. So `Q^SP = s Q^can` and `P^SP = s P_can`, which is regular and graded.
  * The hypotheses of U, re-derived: `Prod b = Prod bt` iff `v'^5/v = 1/(u(1-u)^5)`; with `v = u + O(u^2)`, this picks the fifth
    root and gives `v^(4/5) = u^(4/5)/tau`.
  * Theorem G for `P^SP` is correct but redundant (C6). `P^SP` is itself a regular gauge (`P^SP_0 = P_can,0`), and GI Theorem G is
    stated for every regular gauge; or multiply GI's identity for `P_can` by `s(z)`.
* **s.2.9: correct.**
  * `V_1 = -1/z`. `[u^1](Y^C_Q)_0 = -1/z`: GI Lemma Q(a), with only `n = 1` and `y^(1)_0 = -t_c`. `C0(0) = 1`.
  * One fact is used without being stated (C7): `[u^0](Y^C_Q)_0 = 0` (`y^(0) = -e_1`). Without it, the `u^1` term of `C0` would
    enter. Then `s = 1`.

## 4. PB s.3: the glue `D_0(z) = Phi^can(z,1)/h(z) = c(z;1)/h(z)`

* **4.1 UR s.4.** `Q D = P_can Q^m(v(u))` holds for the amplitudes, with the SAME exponential stripped on both sides; then `D` is
  constant, diagonal and in `Q[[z]]`.
  * UR writes `e^((u_i - u_i^c)/z)`. The constant `e^(u_i^c/z)` is not a power series, so it must not be absorbed into `D`.
  * PB's amplitudes do this correctly. Row 0, colour 0, hatted: `Phi^QRR D_0(s^4 w) = Phi^can`. Confirmed.
* **4.2 WKB uniqueness and the common ring (PB's argument with the missing sentence, C8).**
  * (a) `A^can_n = L^(4n+1) ph^can_n(L^-5)` lies in `Q[L, 1/L]` (Lemma 2).
  * (b) In `Q[[q]]`, with `L = (1 - 3125q)^(-1/5)`, both `A^can` and `A^QRR` are WKB columns of PF5 with the exponential
    `e^(u_0/z)`. The transport operator has kernel `Q L` there (the characteristic roots `zeta^i L` are distinct). So
    `A^can = A^QRR d(z)`, with `d` in `1 + z Q[[z]]`.
  * (c) Hence `A^QRR_n` lies in `Q[L, 1/L]`: the QRR amplitudes are Laurent polynomials in `L`. This is a by-product, and no degree
    bound is needed.
  * (d) `Q[L, 1/L][[z]]` embeds in `Q((s))[[w]]`, where UR s.4 holds. `Phi^QRR` has constant term 1, so it is not a zero divisor.
    Then `Phi^QRR D_0 = Phi^can = Phi^QRR d` gives `D_0 = d`.
    * The rows of UR s.4 near `u = 0` are, by the definition R0, the QRR WKB columns, the same element as in (b). By (c) there is
      no continuation issue.
  * (e) At `L = 1`, `Phi^QRR(z, 1) = Sum_n A^QRR_n(0) z^n = h(z)` by (QRR). So `h(z) D_0(z) = Phi^can(z, 1)`. PROVED, given P(b) for
    the polynomiality, H8 and H9.
* **4.3 (QRR).** The glue needs exactly H9, `A^QRR_n(q=0) = h_n` for all `n`: the flat colour-0 form of `R(q=0) = diag(h_(w_i))`
  (EVID's u = 1 convention). Passing from `R(q=0)` to H9 belongs to the input (QRR), not to the glue; I did not re-check it.
* **4.4 Same object: yes, literally.** Lemma N (EVID s.2.3) says: for each polynomial-type `Phi`, `Ghat(Phi)` is the hatted form of
  `P_can c`, with `K_q = [z^(5q)] c`.
  * The gauge `P_can c` has `P_11(0, z) = c(z)`, because `P_can,11(0,z) = 1`. So `K_q = 0` for all `q` iff `c = 1` iff
    `Phi = C0 [hatted P_can Q^m]_00`, which is POLY's `Phi^can`: the same element of `Q[u][[w]]`.
  * Lemma N, like Lemma 2 for `Phi^can`, builds rows by the cyclic system and so needs (CYC), now PROVED (u).
* **4.5 The rest: confirmed.**
  * `D_i(z) = D_0(zeta^-i z)` and `h_(w_i)(z) = h_(w_0)(zeta^-i z)`: `Sum_(a in w_i) a^(1-2k) = zeta^(i(1-2k)) Sum_(a in w_0) a^(1-2k)`.
  * `c_n(L) = L^(4n) ph^SP_n(L^-5)` is a Laurent polynomial: two rational functions that agree as expansions at `L = oo` are equal.
  * Lemma N's formula for `K_q` is not used (K3).
