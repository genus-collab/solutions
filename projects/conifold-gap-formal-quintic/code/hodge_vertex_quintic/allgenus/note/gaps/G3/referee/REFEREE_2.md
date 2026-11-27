# Referee report on gap cluster G3 (s04, transport), part 2 of 2: items G3-2d to G3-3, conventions, computations, corrections

Continues `REFEREE.md` (verdict table s.0 there; abbreviations FL, CV, LE as there).

## 6. Item G3-2d: `lem:tra-conv` (d) (CV (d), s.3)

`wt(k,b) = (4k + nu_b - 4)/5 >= 3/5` for `k >= 2` (the minimum is at `b = 1`, `nu_1 = -1`); the degree bound with `n + m` points and
`l_j >= 2` gives `Sum k_i <= 3h - 3 + n - m`, hence `m <= 3h - 3 + n`. Both uses named (`lem:tra-F`; the hatted second proof, where only
the first statement is needed) are correct. **Verdict: accepted.**

## 7. Item G3-2e (added): the vertex resummation for `thm:red-main` steps 8-9 (CV (e), s.4)

* (e1) from (d) with `n = 0`. (e2) For fixed `m`, `<L^m>_g = Sum_(l+j=m) binom(m,l) delta^j <t'^l H^j>_g` (finite, multilinear); frames
  with `k >= 2` and `a != 1` have `wt >= 4/5`, `H` has `wt = -1`; with `ord_u v = 1` (this is where `v'(0) != 0` enters: for negative
  exponents `ord_u v^e = e` needs it) and `ord_u delta^j >= 2j`, every term has u-order `>= 2 - 2g + j + 4l/5`; together with (b2) for the
  frame labels the family is summable, and the two groupings give the two sides. Correct. Note that the individual `<t'^l H^j>_g` do
  not vanish for `l + j > 3g - 3` (`t'` has the `z^0 H`-component `-delta`); only their `m`-sums do, consistently with (e1).
* (e3) For `j >= 1`, `<H^j>_g - (F^m_g)^((j))(v)` lies in `Q[v]` by `ML_K(g)` (i) and is homogeneous of exponent `2 - 2g - j < 0`, hence `0`;
  this sharpens RD2 s.3.2 step (3) ("`+ Q[[u]]`") to an equality. The binomial series in `delta/v` (in `u Q[[u]]`) gives
  `c_g (v + delta)^(2-2g)`; the `l >= 1` part is a summable family of elements of `Q[[u]]`. Correct.
* The shorter route (Remark): the tail `j > 3g - 3` has u-order `>= j + 2 - 2g >= g`. Correct; `r3_resum.py` [E] finds the order of
  `S_(3g-3) - c_g (v + delta)^(2-2g)` equal to `g` exactly for `g = 2, 3, 4` (the bound is sharp), and the controls show that both
  hypotheses are needed: `v = u^2` (E1) and `delta = u/2` (E2) give constant u-order terms and no limit.
* Dependencies (K10): homogeneity of `def:set-modcohft`; `ML_K(g)` (the hypothesis of `thm:red-main` (b)); `lem:tra-T` (d)
  (`delta_P` in `u^2 Q[[u]]`); `v'(0) != 0` (O-2). The prover names all of them. **Verdict: accepted.**

## 8. Item G3-2f (added): dilaton re-expansion and `d(0) = 1` (CV (f), s.5)

* (f1) The dilaton equation applied `m_1` times gives the rising factorial `(N)_(m_1)`, `N = 2h - 2 + n + m_2 >= 1` (every intermediate
  step has `2h - 2 + r > 0`; for `h = 1` the anomaly needs `r = 0`, excluded). `Sum_(m_1) binom(N + m_1 - 1, m_1) c^(m_1) = (1 - c)^(-N) =
  kappa^N` when `ord_u c >= 1`; `m_2 <= 3h - 3 + n` by (d). Correct; `r3_resum.py` [F] confirms the series identity mod `u^24`.
* (f2) **Defect (K3).** The proof shows that for FIXED `m_2` the `m_1`-series has no limit. The "Consequence" is about the Taylor
  series of `Fbar_h` at the base point, whose `m`-th term `T_m = Sum_(m_1 + m_2 = m)` mixes several `m_2`, and a cancellation between
  them is not excluded by the written argument. The claim is true; the missing lines: let `c_0 := c(0) != 0`, `A_(m_2) :=
  <X T^(m_2)>_h/m_2!`, `o := min ord_u A_(m_2)` (finite by hypothesis), `S := {m_2 : ord_u A_(m_2) = o}`, `a_(m_2) := [u^o] A_(m_2)`. Then
  for `m >= max S`, `[u^o] T_m = c_0^m Sum_(m_2 in S) c_0^(-m_2) a_(m_2) binom(N_(m_2) + m - m_2 - 1, N_(m_2) - 1)`, `N_(m_2) = 2h - 2 + n + m_2`
  (`ord_u c >= 0`, `[u^0] c^k = c_0^k`). Each binomial is a polynomial in `m` of degree `N_(m_2) - 1`, the degrees are distinct, so the
  sum is a non-zero polynomial and `[u^o] T_m != 0` for all large `m`: `T_m` does not tend to `0`, the Taylor series has no u-adic limit.
  (`r3_resum.py` [F2]: no zero of `[u^0] T_m` for `m <= 80` in a three-term example.) Add this to (f2), or state (f2) for the
  `m`-grouped series directly.
* Consequence and recommendation (K4). The diagnosis is right and is new: `def:tra-dressing` allows `d(0) != 1` (`D_0 = diag(lambda, mu,
  1/mu, 1/lambda, 1)` is diagonal, graded, satisfies `D_0^T eta D_0 = eta` and keeps (A); then `d -> d/lambda`), and for such `P` the
  vertex tensor "`d^n Fbar_h(q_T)`" of `def:red-flatsum` is not defined as a u-adic Taylor series. Defining it by the homogeneous form
  `d^(2h-2+n) d^n Fbar_h(d q_T)` (whose leaf `d q_T + z 1 = d T_(>=2)` lies in `z^2 HH^u_+`) is the right fix, and `thm:red-main` is
  unaffected (`d^(2g-2) c_g (d T_P)^(2-2g) = c_g T_P^(2-2g)` for any unit `d`). Two precisions for the editor: (i) the same proviso
  sits in `lem:set-dilaton` ("hence `d^n Fbar_h(q/kappa) = kappa^(2h-2+n) d^n Fbar_h(q)`") and in S15 s.1.1 and S15 Lemma 2 ("termwise"):
  as identities of Taylor series at `-z 1` they need `kappa(0) = 1`; S15 uses `kappa = d` with `d(0) = 1` and `kappa = 1`, so S15 is
  fine, but the sentence "it is also what the R-matrix action ... produces ..., cf. S15 Lemmas 1-2" should say "with the vertex at
  the dilaton point of the frame-changed theory (unit `d^-1 1`), which is the homogeneous form", not cite S15 Lemma 2 as a proof for
  `kappa(0) != 1`. (ii) `d(0) = 1` for `P_can` rests on `lem:gau-diag` (`p_(0,0) = tau^(-1/4)/varpi_0`, so `p_(0,0)(0) = 1`, and `d = 1/p_(0,0)`
  by (Sp) at `z^0`); that lemma's step 4 carries the `\draftnote` G8 (`p_(0,a)(0) = 1`), so the dependency must be printed (K10). TRF s.2
  found `d = tau^(1/4) varpi_0` by computation (VERIFIED), consistent.
* Dependencies: `lem:set-dilaton` (its citation is a `\draftnote`, R7 list), `lem:gau-diag`. **Verdict: accepted-with-corrections**
  (K3, K4, K10).

## 9. Item G3-2g (added): frame expansion of bivectors (`GAPS_G3_2.md` s.3)

(ii) => (i) by finite type and linearity; (i) => (ii): columns in `K_0`, FL Step 4 in `psi`, regrouping, the summation lemma for
`1 (x) Pi_i` (finite type coefficientwise in the first slot), uniqueness of frame coordinates, FL Step 4 in `psi'`. The summability of
the double family is right (for a coefficient `psi^k' psi'^l' H^a' (x) H^b'` with `a' = b' = 1` the u-order is at least
`(4(k - k') + 1)/5 + (4(l - l') + 1)/5`). This writes the "duality step" of DC1 B2 step 5. `r3_twoslot.py` (s.12) builds `E =
(pr (x) pr)(B)` for random `B` WITHOUT frames and recovers it from its flat components with `a, b != 1`; the control `(pr (x) 1)(B)`
fails both checks. **Verdict: accepted.**

## 10. Item G3-3: `lem:tra-E` step 2, the ring of O-8 (`GAPS_G3_4.md`)

* *Attribution.* Checked: item O8 of `stage2/STAGE2_SYNTHESIS_2.md` s.5 concerns (E1) as evidence for (Q1) (logleg); the O-8 of the
  outline's `\draftnote` is DC1_6 s.7 O-8. The prover's correction of the pointer is right.
* *The ring.* `AA = Q((x))[[u, psi']]` is a power-series ring over a field; it contains `Q[[u, x, psi']]` and `Q[x^-1][[u]]`; `x + psi' =
  x(1 + x^-1 psi')` is a unit with inverse `iota`; `Res_x` is coefficientwise and commutes with `(u, psi')`-adic sums. The diagnosis of
  TRF's ring is right: `R^-1(x)` has infinitely many positive powers of `x` at `u^0` for a general symplectic `R` (`r3_lemE` [P]: all 28
  computed coefficients non-zero), so `iota R^-1(x)` is not in `Q((x^-1))[[psi']]`. For `P_can` the sources only show non-zero terms
  `z^3, z^8, z^13, z^18` of `P_30(0, z)` (GI, `gi_relations_ab_0_0.log`): that `P_can(0, z)` is not a polynomial is VERIFIED to `z^18`, not
  proved; the correction does not need it, since Lemma E is stated for every symplectic `R`. GH's Lemma L is the finite, polynomial
  setting (GH s.1: `[u^N] Ghat` polynomial in `zhat`), where TRF's ring works and no decay of the covector is needed (GH s.2).
* *The proof in `AA` (LE Steps 0-7).* Step 0: `N(psi, -psi) = 0` from (Sp) inverted; the substitution `psi' = w - psi` is an automorphism
  of `B[[psi, psi']]`, so `N = (psi + psi') E_R` with `E_R` a power series, unique (domain); symmetry. Step 1: the contraction is the
  residue, finite modulo `u^(M+1)` because `alpha` decays. Steps 2-4: `E_R = iota N` in `AA`; the `eta^-1`-term has no residue;
  `beta = R^-1(x)^T alpha` in `Q((x))[[u]] (x) H`. Step 5: `R^-1(x)^T = eta R(-x) eta^-1` gives `beta(x) = eta (R f)(-x)` and
  `C = -R(psi')^-1 (R f)_(>=0)(psi')`. Step 6: iff, second slot. I re-derived every sign and transpose. Correct. `r3_lemE_run.py` checks the
  final identity and the intermediate Step-4 form `-R^-1(psi') eta^-1 beta_(>=0)(-psi')` separately.
* *Defects.* (K5) Remark 2 ("In `Q((psi'))[[u, x]]` ...") names a ring that does not contain `alpha` (negative powers of `x`); the
  comparison should be made in `Q((psi'))((x))[[u]]`, where `x + psi' = psi'(1 + x/psi')` is a unit. Side remark; not used. (K6) Step 7
  applies the vector-form lemma to `f = z^-k J^C` and `P`; (A) is stated for `P Y^C` in lowered components. Write the line `P^vec J^C =
  eta^-1 P eta eta^-1 Y^C = eta^-1 (P Y^C)` (`P^vec = eta^-1 P eta = J P J`), so `(P^vec z^-k J^C)_(>=0) = z^-k eta^-1 (P Y^C)_(>=k) = 0`; this
  is the O-7 convention, but the proof of `cor:tra-edges` step 3 must contain it.
* *Evidence remark (K8 (a), not a defect of the proof).* `g3_lemE.py` computes `[psi^m psi'^n] E_R = Sum_j (-1)^j N_(m+1+j, n-j)` (the
  `psi'/psi` expansion of `N/(psi + psi')`). For `m >= 0` this equals the difference quotient `-Sum_N r_N Sum_(k<N) psi^(N-1-k) (-psi')^k
  eta^-1 R^-1(psi')^T` for EVERY `R`, symplectic or not, since the defect `eta^-1 - R^-1(-psi') eta^-1 R^-1(psi')^T` is free of `psi` and
  contributes only negative powers of `psi` (`r3_lemE_run` [C2]). So that LHS does not test divisibility; the prover's separate [DIV]
  check and control C2 do, and the description in `GAPS_G3_4.md` s.4 lists them separately. Fine as it stands.
* **Verdict: accepted-with-corrections** (K5, K6).

## 11. Conventions (`GAPS_G3.md` s.2) and observations (s.3)

* s.2.1 (two settings, `HH^u`, the splitting, `Omega_s`, the substitution): correct (checked in s.5 of part 1). "Lagrangian" means here
  "isotropic and complementary", which is what is used. s.2.2 (summable families, maps of finite type, summation lemma and its
  proof, regrouping): correct; the summation lemma is the one tool that makes every "u-adically convergent" of the sources precise.
* s.3 item 1: the frame coordinates of FL Step 4 in weight `w` give the basis `v^m frame(k,a)` of `K_0 cap V_w` and `dim = d(w)`
  directly (no `prop:mod-dii`), and `V_w = (K_0 cap V_w) (+) span{v^m z^k H}`. Correct. **(K9)** The parenthesis "DICT Cor. 4, proved
  there through the log map and D(ii)" misattributes: DICT s.6 says "(d) and Cor. 4 do not use D(ii)"; Cor. 4 there is a statement
  about the kernel of the log map, and D(ii) is what identifies that kernel with the `K_0` cut out by the `Pi_i`. Rephrase; the
  conclusion (FL (a) in weight `w` gives the `Pi_i`-version with no D(ii)) stands. The "assumed homogeneity" of DICT s.6 is proved
  by FL Step 0 (0.3) and Step 2, as stated. Items 2-4, 6, 7: correct (item 7 is s.2.4 of part 1). Item 5: no comment.

## 12. Independent computations (own code in this folder; exact over Q; `py -I -B`; working directory the scratch folder)

| script (log) | what it checks, by a route different from the prover's | result |
|---|---|---|
| `r3_frames.py 60` (`r3_frames_60.log`) | `Y^C` from the residue definition (double sum of the two exponentials) at `v = 1, 2`: model system, `Y^C(0) = -e_1`, closed form; frames by generic Gaussian elimination of `Pi_i(frame) = 0` (pairing computed directly; no shift rule, no `Theta`, no top-down rule), `k <= 60`; homogeneity by comparing `v = 2` with `v = 1` | 0 mismatches; 244 frames, all systems regular, `Pi_i(frame) = 0` to `i = k + 12`; trivial frames exactly the ten `k < nu_a + 1` (missing from the eight: `(2,2), (3,2)`); `f_(j0) = eps_a v^(nu_a+1)/(nu_a+1)!`; exponents integral, `>= 1`, `f_j(2) = 2^e f_j(1)`; examples. Controls: `Y^C` without its `z^-n`, `n >= 5`, terms: 214 of 224 frames not annihilated; wrong exponent `(4(k-j)+nu_a)/5` refuted on 1488 coefficients; the rule `k <= nu_a + 1` refuted |
| `r3_uadic.py 10 12` (`r3_uadic_10_12.log`) | `v(u)` from the Pochhammer quotient and the ODE `F g' = (5/4) F' g`; mod `u^10`: (U1) top-down solve = frame formula for a finite `q`; (U2) an INFINITE `q` (`q_(k,b) = (b+1)(1 + k u)`, all `k`): `sigma` from the frame formula, then `Pi_i(q) = (Theta sigma)_i` with `q` infinite; (U3) sharpness of the cut | `v = u + (5/9) u^2 + (220/567) u^3 + ..` (GI s.5); U1 true; U2: 0 failures (`i <= 12`), no frame beyond the predicted cut visible; U3: last changing frame `= j + floor((5m-1)/4)` exactly. Controls: the frame sum cut 2 steps early fails for all 13 `i`; `v = 1 + u`: the `u^0` parts of `sigma_0` and of `Pi_0(q)` change at 48 resp. 61 values of the cut (no u-adic limit) |
| `r3_twoslot.py 7 14` (`r3_twoslot_7_14.log`) | G3-2g: `E = (pr (x) pr)(B)`, `B` random symmetric of bidegree `<= 14`, built without frames, mod `u^7` | annihilated by `Pi_i` (`i <= 17`) in both slots; `E` = frame double sum of its flat components with `a, b != 1`; symmetric. Control `(pr (x) 1)(B)`: 1110 non-zero columns in the second slot, reconstruction fails |
| `r3_lemE_run.py 7 5` with helpers `r3_lemE.py` (`r3_lemE_7_5.log`) | Lemma E with `R` the Cayley transform of a random `X = Sum_(k<=4) z^k X_k(u)`; `E_R` by the difference quotient and by both expansions of `N/(psi+psi')`; `f = z^-1 J^C, z^-3 J^C` (`v = u tau^(-5/4)`) and a random `f` with `[z^-n] f = O(u^floor(n/3))`; mod `u^7`, `psi'^0..psi'^5` | (Sp) to `z^28`, `R^-1(0, x)` non-polynomial (28 of 28); the three `E_R` agree; LHS = RHS and = the Step-4 form for all three `f` (30 of 30 components non-zero). Controls: sign flip fails; an extra `z^0` term in `f` fails; constant `f` gives 4 different `u^0 psi'^0` values for 4 cuts; non-symplectic `R`: the two expansions disagree |
| `r3_resum.py 24` (`r3_resum_24.log`) | G3-2e, G3-2f scalar identities | as in s.7-s.8: tail order `= g` (`g = 2, 3, 4`); E1, E2 fail as predicted; `Sum binom(N+m-1,m) c^m = kappa^N` (`N = 1, 2, 5`); F1 no limit; F2 no zero for `m <= 80` |

Also re-run: the prover's six logs (s.1 of part 1), byte-identical. Run times below 30 s each.

## 13. Consolidated corrections for the editor (none changes a resolution)

* **K1** (G3-1) FL s.1 (c), graded clause, and the first line of FL Step 4: "finitely many non-zero members, each homogeneous" ->
  "finitely supported families in `Q[v]`; on `K_0 cap V_w`, the `(k,a)`-member a rational multiple of `v^((w - 4k - nu_a - 1)/5)`".
* **K2** (G3-2c) "GI s.2 Lemma V (= `lem:gau-vanshape`)" -> "`lem:gau-vanshape` (DB S23; GI s.2, the shape statement following Lemma V;
  GI's Lemma V itself is `lem:gau-vanmod`, NOTATION_2 s.9)", in `GAPS_G3.md` s.0, s.1 (i) and CV (c).
* **K3** (G3-2f) add the polynomial-in-`m` argument of s.8 to (f2), so that the Taylor series grouped by `m` is covered.
* **K4** (G3-2f; `lem:set-dilaton`; S15 s.1.1, Lemma 2) the "hence" clause of `lem:set-dilaton` and S15's termwise identity are
  identities of Taylor series at `-z 1` only when `kappa(0) = 1`; adopt the prover's definition of the vertex tensor by the homogeneous
  form in `def:red-flatsum`, and phrase the S15 reference as in s.8.
* **K5** (G3-3) Remark 2 of `GAPS_G3_4.md` s.3: `Q((psi'))[[u, x]]` -> `Q((psi'))((x))[[u]]`.
* **K6** (G3-3, `cor:tra-edges` step 3) write `P^vec J^C = eta^-1 (P Y^C)`, so (A) in lowered form gives the vector-form hypothesis.
* **K7** (G3-2b) if the note keeps "homogeneity is not needed for the convergence", define `e_0` as the least `v`-exponent of the
  finitely many flat correlators with `Sum k_i <= D`.
* **K8** (a) evidence remark of s.10 (no change needed); (b) state `Theta = 1 + N` and its inverse once, as a lemma of `sec:tra-tower`.
* **K9** (observations s.3 item 1) DICT Cor. 4 does not use D(ii) (DICT s.6); rephrase as in s.11.
* **K10** print the dependencies with the resolutions: G3-2b, G3-2e on `def:set-modcohft` (homogeneity, `\draftnote`); G3-2e on
  `ML_K(g)`, `lem:tra-T` (d), `v'(0) != 0`; G3-2f on `lem:set-dilaton` (`\draftnote`) and, for `P_can`, on `lem:gau-diag` (`\draftnote` G8).
* For `design/SOURCE_ERRATA.md`: RD1 s.1.3 and DC1 D9 list eight trivial frames; there are ten (add `(2,2)`, `(3,2)`); DC1 D10 is right.

## 14. Files (this folder; nothing else was written outside it and the scratch folder; no existing file was modified)

`REFEREE.md`, `REFEREE_2.md`; scripts `r3_frames.py`, `r3_uadic.py`, `r3_twoslot.py`, `r3_lemE.py` (helpers), `r3_lemE_run.py`,
`r3_resum.py`; logs `r3_frames_60.log`, `r3_uadic_10_12.log`, `r3_twoslot_7_14.log`, `r3_lemE_7_5.log`, `r3_resum_24.log`. Run, e.g.,
`py -I -B r3_frames.py 60` with the scratch folder as working directory (give the script's path); the logs print no paths.
