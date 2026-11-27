# Referee report on gap cluster G4, part 2: items 3 and 4, findings F1-F10, own checks, corrections, files

Conventions and abbreviations: REFEREE.md (header) and GAPS_G4.md s.1. Corrections are numbered R-1 .. R-8 across both parts.

## 1. Item 3: `thm:mod-Lreg`, `lem:mod-abdks`, `prop:mod-period`, `lem:mod-vanper` (GAPS_G4_4) -- accepted-with-corrections (minor)

### 1.1 The quoted theorem
* ABDKS v3 (the current version; source numbering): Theorem 7.1 (label `cor:y=z`) is the first statement of s.7: "Let `x(z)` be an
  arbitrary rational function. Then all differentials `omega^(g)_(m,0)` of topological recursion with the spectral curve `(x = x(z),
  y = z)` and the initial form `omega^(g)_(1,0) = -y dx` are given by" (7.9); vertex weight (7.6) `W_i = -exp(-w_i(S(w_i hbar d_(z_i))
  - 1) x(z_i)) dz_i/w_i`; edge weight (7.7); (7.10) the one-point formula; Remark 7.3 on limits to the discriminant. The s.7
  introduction restricts to "rational x such that dx has only simple zeros"; the standing assumptions of s.1.1 (simple critical
  points of `x`, `y` regular with `dy != 0` there) are the TR hypotheses. For `x = W_v`, `v != 0`: `dW_v = (p^5 - v) dp/p^2` has five
  simple zeros and poles at `0`, `oo` only; `y = p`. So Thm 7.1 applies to `x = W_v` directly: **W5 dissolves**, as G4 says.
* Kernel (1.1): `(1/2) Int_z^(sigma z) omega_(2,0)(z_1, .) / (omega_(1,0)(sigma z) - omega_(1,0)(z))` with `omega_(1,0) = -y dx`, i.e.
  `(1/2) Int_z^(sigma z) B/((y(z) - y(sigma z)) dx)`: EO's printed kernel, minus the engine's. So `omega^A_(g,m) = (-1)^m omega_(g,m)`,
  and F4 (XYR l. 37-38 says "minus the kernel of arXiv:math-ph/0702045"; it is equal to it) is a correct source erratum; XYR's
  conclusion `omega^code = (-1)^n omega^A` stands. The sign is also confirmed numerically at `(1,1)`: `-d Phi_1` from (7.10) equals the
  TR `omega_(1,1)` (s.3.3 C3 with s.3.2 B1 below).
* **R-4 (recommended).** GAPS_G4_4 s.1 transcribes (7.7) as `W_(i,j) = exp(hbar^2 w_i w_j S S (z_i - z_j)^-2) - 1`. ABDKS print
  `W_(i,j) = (exp(hbar^2 w_i w_j S(w_i hbar d_(z_i)) S(w_j hbar d_(z_j)) (z_i - z_j)^-2) - 1)/hbar^2 = w_i w_j/((z_i - z_j)^2 - (hbar^2/4)
  (w_i + w_j)^2)` (XY s.1 has the closed form right). Harmless for `lem:mod-abdks`: for fixed `n` the graphs are the finitely many
  simple graphs on `n` numbered vertices; in `W_(i,j) = Sum_k (hbar^2/4)^k (w_i + w_j)^(2k) w_i w_j (z_i - z_j)^(-2k-2)` every `hbar^2`
  comes with `(w_i + w_j)^2` and one more `(z_i - z_j)^-2`, and in `W_i` every `hbar^(2k)` with `w_i^(2k+1)`; so at fixed `hbar^(2h)`
  the `w`-degrees, hence the `r_i`, are bounded. Rewrite step (1) of the proof with the printed weight.

### 1.2 The proofs, step by step
* `lem:mod-abdks` (1): the ring `Q[v][p_i, p_i^-1, (p_i^5 - v)^-1, (p_i - p_j)^-1]` is right (`W_v^((2k))` in `Q[v][p, p^-1]`,
  `D_x(f dp) = d(p^2 f/(p^5 - v))`). (2): the division argument is correct: on the generic point of `{p_i = p_j}` and of `{p_i = 0}`
  the other denominator factors are units, regularity there (EO Thm 4.2) gives divisibility for each fixed `v != 0`, and the
  remainder of the division by a monic divisor over `Q[v][p_k : k != i]` vanishes for all `v != 0`, hence identically; the
  divisions can be done one factor at a time. No negative power of `v` arises. The optional bound `K = 6h - 4 + 2n` is argued from
  (2-35), i.e. from item 1: harmless but it mixes the two routes; the lemma needs no bound on `K` (R-4, same paragraph).
* `thm:mod-Lreg`: correct. In one variable `Sum_b Res_(varsigma_b) = -Res_(p_l = oo)` (the form is holomorphic at `p_l = 0`);
  `Res_(p=oo) p^a dp/(p^5 - v)^K = -C(K+j-1, j) v^j` if `a + 1 = 5(K + j)`, else `0`; iteration is legitimate because after each step
  the result is again polynomial over the remaining `prod (p_m^5 - v)^K`. Legs in `Q(v)[x]` drop by `lem:mod-xfun`, so the pairing of
  `LL = Q[v][p] + Q(v)[x]` lies in `Q[v]` (the decomposition of a leg need not be unique; any one works). Remark (ii) correct.
* `prop:mod-period`. Step 1: `D_x^r(a dp) = d(delta^(r-1)(a/x'))` for `r >= 1`; the `r = 0` term is `-[w^1] e^E dp = 0`; `a_r = -Sum_j
  (1/j!) Sum_(k_1+..+k_j = h) prod eps_(k_i)`, `r = 2h - 1 + j`, `1 <= j <= h`; the engine sign `omega_(h,1) = -d Phi^A_h`: all correct.
  Step 2: `deg eps_k <= 4 - 2k` (`(p^4/4)^((2k))` vanishes for `k >= 3`), every term has degree `<= 5 - 10h` (the bound is attained:
  `rf4_check_c`). Step 3: the pure part is re-derived (`xi = v/p`, `delta_0 = d/dxi`, `-p^2/v = -v/xi^2`, `a^0_r = -(2v)^(-2h) xi^(r+1)
  c_(h,r)`, pure term `v (2v)^(-2h) (r-1)! c_(h,r)`); every correction raises the order (by `>= 5` in fact: `eps_1, eps_2` corrections by
  5, `1/x'` and `delta` by 5); `delta_0` raises the order of every surviving monomial by exactly 1, so no term returns below 0.
  Correct. Step 4: `c_(h,r) = ((-1)^j/j!) [T^h] phit^j` (`E_0 = -t phit(epsl^2 t^2)`), the binomial resummation (the `j = 0` and
  `j > h` terms vanish), `1 + phit(s^2) = artanh(s)/s`, and Lagrange-Buermann for the Laurent function `H(w) = w^(1-2h)` (valid:
  `[s^n] H(w(s)) = (1/n) Res H'(w) tanh(w)^-n dw` holds for Laurent `H`, by the change of variables and one integration by parts),
  `w coth w = Sum 4^n B_2n w^(2n)/(2n)!`: `S_h = -4^h B_2h/(2h)`. Correct; it is TRM s.4.1 (v)'s `A_h` times `4^h`. The Borel
  integral of XY is indeed not needed (F5). Step 5 (contour lemma; branch constants of `log p` drop since `omega_(h,1)` has zero
  residues) and (c) (`<p^5>_h` in `Q[v]` homogeneous of exponent `2 - 2h < 0`): correct.
* Verdict: **PROVED modulo ABDKS Thm 7.1 and EO Thm 4.2** (both in (Q)); accepted with R-4.

## 2. Item 4: base `h = 1` of `thm:mod-V` (GAPS_G4_5) -- accepted

* `lem:mod-w11` (i) re-derived: `Int_(sigma q)^q B(p_0, .) = 2 dp_0 Sum_(m even) B_m zeta^(m+1)/(m+1)`, `y(q) - y(sigma q) = 2 Sum_(k odd)
  y_k zeta^k`, `dx = 2 alpha zeta dzeta`, `B(zeta, -zeta) = -[(4 zeta^2)^-1 + beta_00 + O(zeta^2)] dzeta^2`; the `zeta^-1` coefficient
  gives `C_a`. `beta_00 = (1/6){t; zeta}(0)` is the standard Schwarzian term. (ii): Fubini is legitimate (`K(p_0, q)` is
  `O(p_0^-2)` on the ray uniformly for `q` on the small circles); `Int_(ray) B(p_0, q) = -dq/q`. The Airy anchor (`-(1/16) dp/p^4`,
  then `lem:mod-scale` with `kappa = -1/2`: `(1/8) dp/p^4 = 3!! <tau_1>_1`) re-checked.
* The evaluation at `v = 1`, `a = 1`, re-done by hand: `x(1+t) = 5/4 + (5/2)t^2 + (5/4)t^4 - t^5 + ..` (the `t^3` term vanishes because
  `x'''(varsigma_b) = 6 varsigma_b - 6 v/varsigma_b^4 = 0`), `zeta = t + t^3/4 - t^4/5 + ..`, `t = zeta - zeta^3/4 + zeta^4/5 + ..`,
  `y_1 = 1`, `y_3 = -1/4`, `{t; zeta}(0) = -3/2`, `beta_00 = -1/4`, `-dq/q = -(1 - zeta + zeta^2/4 + ..) dzeta`, `C_1 = -(1/10)(1/4 -
  1/48 - 1/16) = -1/60`. The rotation `rho_b` and the homogeneity step are correct. Total `<-log p>_1 = -1/(12 v)`: **PROVED**.
* Both cross-checks re-done: the closed form with `u = p^5`, `w = u - v` (`G(-v) = -(2/3)/v`, path independence since no `w^-1`
  term) gives `-1/(12 v)`; `prop:mod-period` at `h = 1` gives `-(B_2/2) v^-1`.
* **R-6 (recommended).** "nothing quoted" is right, but the ray-independence used in "All v" needs that `omega_(1,1)` has zero
  residues; say why: by (i) each `B_m(p_0) dp_0` is the `p_0`-differential of a rational function, so `omega_(1,1)` is exact (or cite
  EO (4-26)). Same for "`O(dp/p^2)` at `oo` and holomorphic at 0" in (ii), which also follows from (i).
* Note on the evidence: since `x'''(varsigma_b) = 0` at every branch point, the local involution is `sigma_b(t) = -t + O(t^4)` in
  `t = p - varsigma_b`; a "linearised involution" is therefore NOT a failing control for any `(1,1)` quantity (the referee's first
  control did not fail; replaced by a wrong `sigma = -t + t^2/a`, which gives `-29/96`). G4's controls (`beta_00` dropped: `1/24`;
  `y_3` dropped: `-11/96`) are meaningful.

## 3. The referee's own exact checks (VERIFIED; evidence only; scripts and logs in this folder; stdlib, `py -I -B`)

Code written from scratch (`rf4_field.py`: `Q(zeta5)` in a redundant 5-term basis, power series; `rf4_model.py`: T-system, `M_0`,
(2-31) by residue extraction, local involution by a fixed point); it imports nothing from G4 or the project. Run times 98 s, 2 s, 1 s
(one process each; the order-9 run of check A was repeated and gave an identical log).

### 3.1 `rf4_check_a.py 9` -> `rf4_check_a_9.log` (item 1): see REFEREE.md s.1.4. ALL CHECKS PASSED.

### 3.2 `rf4_check_b.py` -> `rf4_check_b.log` (items 1 (b), 2, 4)
* B1: `omega_(1,1)` of `x = kappa W_1` (`kappa = 1, 4`) by the recursion (2-28) itself = closed form `/kappa` at `p_0 = 3, -2, 1/2, 7/3`.
* B2, the `(1,1)` graph sum used for R-1 (derived here; `X1 := [lam^-1] R^-1` from (2-31), `X1[j][i] = (R^-1_1)^j_i`, `c_i = varsigma_i`):
  `alpha_(1,1)(e_i) = Delta_i^(1/2) + [leaf] (R^-1_1 e_i) psi + [dilaton] T_2 at a genus-1 vertex with two points, pushed forward
  (`pi_* psi_2^2 = kappa_1`, `Int kappa_1 = 1/24`) + [loop] the genus-0 vertex with the degree-0 edge term `-Sum_a (R^-1_1 e_a) (x) e_a`
  (`R^-1_1` is symmetric by symplecticity), `|Aut| = 2`, `Int rho_* 1 = 1`. With the dilaton leaf `s psi (1 - R^-1(psi)) 1`, (2-35) at
  `(1,1)` becomes, after the factors `sqrt(Lambda)` cancel,
  `omega_(1,1)/dp = Sum_i [ A_i (p - c_i)^-2 + (1/24) d/dp( (p - c_i)^-2 / x'(p) ) ]`,
  `A_i = (1/24) Sum_j X1[j][i] c_j/c_i - (s/24) Sum_k X1[i][k] c_i/c_k - (1/2) X1[i][i]`.
  Result: `s = +1` (PPZ, DOSS) matches TR at all 4 points for `kappa = 1` and `kappa = 4`; `s = -1` (DNOPS s.2.3 as printed) fails at
  all 4; control `X1 -> -X1` fails at all 4. E.g. at `p_0 = 3`, `kappa = 1`: TR and `s = +1` give `-25688745/27437936768`, `s = -1`
  gives `-13829535/27437936768`.
* B3: `Res_(varsigma_b) x^k omega_(1,1) = 0` (`b = 0..4`, `k = 0..6`); control `p^k`: 40 of 45 non-zero; `Res_0 x^5 omega_(1,1) = -15/4`
  (evaluated by hand in the script, see REFEREE.md s.2).
* B4: `<-log p>_1 = -1/12` at `v = 1` by pushing the ray integral into the five residues (kernel `(1/2) log(sigma q/q)`, i.e. a route
  that does not use the local formula `C_a`), and from the closed-form integral; control (wrong involution) `-29/96`.

### 3.3 `rf4_check_c.py 7` -> `rf4_check_c_7.log` (item 3)
* C1: `S_h = -4^h B_2h/(2h)` for `h = 1..20` from the finite sum; control (`1/(2k+2)`) agrees for no `h`.
* C2: (7.10) as transcribed in G4_4 Step 1, for `x = A p^4 + V/p`, `(A, V) = (1/4, 1), (1/4, 3), (1, 2)`, `h = 1..7`: no negative power at
  `p = 0`, `Phi_h(0) = V (2V)^(-2h) S_h = -(B_2h/2h) V^(1-2h)`, top exponent at `oo` equal to the bound `5 - 10h`.
* C3: `-d Phi_1` = the closed form of `omega_(1,1)` through `p^80` at `v = 1, 3` (with B1: (7.10) with the engine sign = TR at `(1,1)`).

### 3.4 What G4's own scripts show
`g4_rtr.py 10`, `g4_w11.py`, `g4_period.py 5` were read; their logs are consistent with the claims of G4_2 s.8, G4_3 s.4, G4_4 s.5,
G4_5 s.5, and their controls fail as stated. They and the referee's scripts agree wherever they overlap.

## 4. The findings F1-F10 of GAPS_G4_5 s.6

| F | referee |
|---|---|
| F1 | agree (formal Laplace reading with the finite part; consistent with DOSS's prescription) |
| F2 | agree (residue pairing for `LL`; circles only for legs holomorphic at 0) |
| F3 | agree: `lem:set-M0sympl` closes R7's "`def:set-M0` symplecticity". It does NOT touch the companion R7 item "`def:set-modcohft` homogeneity" |
| F4 | agree (checked in ABDKS (1.1), EO (4-3), (4-9) and XYR l. 37-38) |
| F5 | agree |
| F6 | agree |
| F7 | agree, with R-3 (what remains quoted is DOSS's theorem via DNOPS s.2.6, not "Thm 8.1") |
| F8 | agree |
| F9 | agree (numbers re-counted in the EO source) |
| F10 | agree |
| new | R-1 (dilaton-leaf sign; required), R-2 .. R-8 |

## 5. Corrections for the editor (none changes a statement of the sources' mathematics)

1. **R-1 (required).** Dilaton leaf = `psi (1hat - M_0^-1(psi/(kappa varsigma^4)) 1hat)`, PPZ Def. 2.13's `T`; DNOPS s.2.3's printed
   `-psi (Id - R^-1(psi)) e_1` is a misprint (DOSS, PPZ, TR at `(1,1)`). Edit GAPS_G4_2 s.1 (b), s.4, s.7 bullets 2-3, GAPS_G4.md s.1;
   add the erratum to `design/SOURCE_ERRATA.md`; in `prop:mod-Rtr` print the PPZ form, as the outline does.
2. **R-2.** Replace "(push-forwards preserve degree)" in the transfer by the bookkeeping of REFEREE.md s.1.3.
3. **R-3.** Cite DOSS Thm 4.1 via DNOPS s.2.6 (+ Prop. 4.3 (1), the sentence before Thm 8.1) for "TR = (2-35)", and DOSS's vertex
   weight for the TFT normalisation used in the transfer (i). Key: root `DuninBarkowski:2012bw` (add to BIBPLAN). Since the conditionality must be exact
   (OUTLINE rule 8), either add [DOSS] Thm 4.1 to the (Q) list of OUTLINE s.0.2 (recommended), or let the DNOPS entry of (Q) read
   "s.2.6 (2-30)-(2-35) and the sentence before Thm 8.1, at `a = (0,0,0,0,4v)`" instead of "Thm 8.1".
4. **R-4.** ABDKS (7.7) has `/hbar^2`; rewrite step (1) of `lem:mod-abdks`; drop or reroute the optional bound `K = 6h - 4 + 2n`.
5. **R-5.** `eta(H^a, H^b) = 5 delta_(a + b = 3 mod 5)`; indices of `bmod`, `eps`, `nu` mod 5 (GAPS_G4.md s.1, `lem:set-M0sympl`).
6. **R-6.** Item 4: say why `omega_(1,1)` is residue-free and holomorphic at `0` (from `lem:mod-w11` (i)).
7. **R-7.** `lem:mod-laplace`: the "Meaning" sentence is motivation, not part of the definition (keep it as a remark).
8. **R-8.** Lemma T (ii): choose the discs `U_b` disjoint from `J = {p_2, .., p_n}`.

After R-1, the `\draftnote`s D6/G2, G3, W5 and D5/G1 of OUTLINE_6 can be deleted; `prop:mod-Rtr`, `lem:mod-xfun`, `thm:mod-Lreg` and
`thm:mod-V` (base `h = 1`) leave the R7 list. The conditionality of the chain is unchanged: items 1 and 3 rest on (Q) (DNOPS/DOSS;
ABDKS Thm 7.1; EO Thm 4.2), items 2 and 4 on nothing quoted beyond EO Thms 4.2, 4.6.

## 6. Files in this folder (all new)

| file | content |
|---|---|
| `REFEREE.md`, `REFEREE_2.md` | this report |
| `rf4_field.py`, `rf4_model.py` | referee's own `Q(zeta5)` arithmetic, series, T-system, `M_0`, (2-31) by residues, local involution |
| `rf4_check_a.py` -> `rf4_check_a_9.log` | item 1: `X = M_0(t)^T` for four `kappa`, (2-32), IBP form, symplecticity; controls |
| `rf4_check_b.py` -> `rf4_check_b.log` | TR `omega_(1,1)` vs closed form; the `(1,1)` graph sum with both dilaton signs (R-1); item 2 residues; item 4 by the ray |
| `rf4_check_c.py` -> `rf4_check_c_7.log` | item 3: `S_h` (`h <= 20`), the (7.10) terms at `0` and `oo` (`h <= 7`), `omega_(1,1) = -d Phi_1` |
| `eqcount.py`, `eo_assemble.py`, `eo_appA_count.py` | reading aids (equation numbering from the TeX sources); not evidence |

Run (working directory = a scratch copy of this folder): `py -I -B rf4_check_a.py 9`, `py -I -B rf4_check_b.py`, `py -I -B
rf4_check_c.py 7`; each log ends with "ALL CHECKS PASSED".
