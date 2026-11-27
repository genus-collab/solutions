# Referee report on the gap-closing note of cluster G2 (Theorem U, the canonical gauge, the colours of `D`; s03), part 1

Task key `gaps/G2/referee`. Author-facing; not note text. AI referee; no human has read this. Object: `NOTE/gaps/G2/GAPS_G2.md` ..
`GAPS_G2_6.md`, the scripts `g2_*.py` and their logs.
* This file: verdict (s.0), what was read (s.1), items 1a-4 (s.2.1-s.2.8).
* `REFEREE_2.md`: item 5 (s.2.9-s.2.10), required and recommended corrections (s.3).
* `REFEREE_3.md`: the referee's own exact checks (s.4), the dependency of item 5 on cluster G1 (s.5), refutation attempts (s.6),
  files (s.7).
Labels: PROVED, VERIFIED (exact; range, field or prime stated), CONJECTURAL, SKETCH, OPEN. Verdicts: accepted /
accepted-with-corrections / rejected. Notation of the G2 note (its s.2) unless said otherwise.

## 0. Overall verdict: ACCEPTED WITH CORRECTIONS

* Every proof of s.3-s.9 was checked line by line on paper (s.2 below). **No mathematical error was found.** No step hides
  a finite check, and no hypothesis is used beyond those the note names. The proofs of s.3-s.8 use only (H), (EP) and
  characteristic 0, as claimed. The main new ideas are correct and are real simplifications: `lem:gau-lead`, which gives the
  uniqueness clause in its strong form "every solution in `gl_5(k[[u]])[[z]]` is `P^nat c(z)`", and `lem:gau-degree`, which gives
  graded gauges with no `zeta` and no averaging.
* The resolution "PROVED" is right for items 1a-1h, 2, 3 and 4.
* For item 5 (G6, `cor:gau-colours`) it is right only as a conditional statement. The computation, the covariance of `P_can` and
  `lem:gau-colmod` are PROVED. `lem:gau-colrows` is PROVED from named s02 inputs. G2 calls one of them OPEN (R1). Cluster G1's
  note, written in parallel, now reports these inputs closed: `prop:set-cyclicLR` PROVED, the colour-`i` passage of (QRR) PROVED,
  (M) CITED, `cor:gau-rowsNN` written. Its referee is not yet on file (`G1/referee/` was empty at the time of this report). So
  item 5 is **PROVED given cluster G1 (pending its referee) and `lem:gau-wkb` stated over `Q(zeta)`**, not "PROVED" outright. The
  G2 summary's one-word resolution for G6 and for `lem:gau-colrows` overstates this; the note's own text (s.9.4 "Status", s.11
  item 2) does not. Details in s.5.
* Nothing is FOUND or CITED. The search record (s.1) is confirmed: a grep of the whole package for the key phrases, and a
  reading of UR s.3, GI s.1-s.5, GR, PS s.4, NFS item 2 and `natural_form/verify_theory/`, find no derivation of any of the
  steps. UR s.3 and GI give statements and one-sentence arguments only.
* All six G2 logs reproduce byte for byte (carriage returns aside) from copies run in a scratch folder (`rerun_record.log`).
  This referee also ran three exact checks of its own, on routes independent of the G2 code (s.4). All their controls fail as
  they must:
  * the projection facts for **general** weights, with `pr` computed by honest linear algebra;
  * the uniqueness clause, `lem:gau-lead` and `V(A~, A~) = k[[z]] 1` for **random** one-pole systems mod `2^61-1`, with a new
    failing control. Without the pole the solution module has rank 5, not 1;
  * the colour covariance of the actual row 0 at large radius (Steps 1-3 of `lem:gau-colrows`), by an exact WKB computation of
    PF5 over `Q(zeta)` with constants fixed by `h_(w_i)`.
* Corrections (part 2, s.3) are of four kinds:
  * the status of item 5, as above;
  * notation clashes not registered in G2 s.11: `gamma` (= 7/10), `beta` (frame parameter), `k` (orders vs field), `S`, `V`;
  * two imprecise sentences (s.9.4 Step 2; s.0 "one dependency");
  * two evidence lines that are not tests.

| item (G2 s.0 #) | outline place | G2 resolution | verdict |
|---|---|---|---|
| 1a | `thm:gau-U` step 3 (G1) | PROVED | **accepted** |
| 1b | step 4 (integral commutant) | PROVED | **accepted** |
| 1c | step 5 (G2, projection formula) | PROVED | **accepted** |
| 1d | step 6 (pole bound) | PROVED | **accepted** |
| 1e | step 8 (G3, two projection facts) | PROVED | **accepted** |
| 1f | step 10 (G4, uniqueness clause; `lem:gau-lead`) | PROVED | **accepted** (one optional wording point, s.3 R2) |
| 1g | step 11 (G5, graded clause; `lem:gau-degree`) | PROVED | **accepted** |
| 1h | step 12 (G7, field `Q`; D15) | PROVED | **accepted** |
| 2 | `prop:gau-graded` step 1 [M] | PROVED | **accepted** |
| 3 | `lem:gau-diag` step 4 (G8) | PROVED, needs (N0) | **accepted** (scope statement required, as G2 says) |
| 4 | `thm:gau-G` Cor. (a) (G9) | PROVED | **accepted** |
| 5 | `cor:gau-colours` (G6) incl. `lem:gau-colmod` | PROVED (with an inherited input) | **accepted-with-corrections** (s.5; C1-C4) |
| 5' | `lem:gau-colrows` (input (C2) of G6's computation) | PROVED from four named inputs | **accepted-with-corrections** (s.5; C2, C3) |

## 1. What was read and re-run

* G2 note in full: `GAPS_G2.md` .. `GAPS_G2_6.md`; all seven scripts; all six logs.
* Design: `OUTLINE_4.md` (s03), `OUTLINE_3.md` s02 items 2-10, `OUTLINE_2.md` s.0.3-s.0.4, `NOTATION.md`, `NOTATION_2.md`, `LABELS.md`
  rows 63-102, `design/OUTLINE_REVIEW.md` R1, R7, `design/DECISIONS.md` (D7, D15), `DIGEST_B` D10-D15, S12-S18 and s.7 G1-G11.
* Sources:
  * UR s.1-s.5;
  * GI s.0-s.5 and its referee GR;
  * PS s.3-s.4;
  * FQ2 s.1 ("Genus-0 data", "Rows") and CS s.1 ("Definition", "Covariance");
  * NFS item 2;
  * `natural_form/verify_theory/` (listing; the docstring of `v2_gauge.py`).
* Cluster G1, for the dependency of item 5: `GAPS_G1.md` s.0, `GAPS_G1_2.md` s.2.0-s.2.1, `GAPS_G1_3.md` s.2.10-s.2.11,
  `GAPS_G1_4.md` s.3.1.
* Grep over `formalquintic/hodge_vertex_quintic/` for "regular semisimple", "commutant", "trace form", "uniqueness clause",
  "one constant per order", "every regular gauge". The hits outside the design files, the G2 folder and a G4 remark are
  quotations of the statements (MEC s.0, S01 E1-E2, PB H2, PS2, GI, GR, UR). **Search record confirmed.**
* Re-runs: `g2_linalg.py 6`, `g2_lemma2.py 14`, `g2_gauge.py 15 32`, `g2_G.py 15 32 10`, `g2_unique.py 16 24 2`, `g2_colours.py 10`,
  from copies in the scratch folder. Every log is identical to the stored one (`rerun_record.log`).

## 2. Item by item

### 2.1 Items 1a, 1c, 1d: `lem:gau-comm` (a)-(e), (g) (GAPS_G2_2 s.3). Accepted

* **(a)** The induction `calB_a^(m) b~_(a+m) = calB_a^(m+1)` is correct. So is the formula for `calY^j`.
* **(a), Leibniz argument.** A permutation contributes only if `pi(a)` lies in `{a, a+1}` for every `a`. Once one `pi(a) = a+1`,
  injectivity forces the whole 5-cycle, which is even and contributes `(-1)^5 Bprod`. So the determinant is `x^5 - Bprod`.
* **(b)** Re-derived: `[A~, X] = Sum_a (b~_a xi_(a+1) - xi_a b~_(a+m)) E_(a,a+m+1)`.
* **(c)** The shifts are permuted bijectively by `ad A~`, so commuting is a condition on each shift component separately. On
  `gl^(m)` the condition is the recursion `xi_(a+1) = xi_a b~_(a+m)/b~_a`. It closes up around the cycle, since the product of the
  ratios is 1, so the kernel is the line `F A~^m`.
* **(d)** Trace orthogonality `tr(Z[A~, Y]) = 0` for `Z` in the commutant: correct.
* **(d)** `tr(X A~^(5-m)) = Bprod lambda_m(X)`, using `calB_a^(m) calB_(a+m)^(5-m) = Bprod`: correct.
* **(d)** Dimension count: the image of `ad A~` on `gl^(m-1)` has dimension `5 - 1 = 4`, which equals `dim ker lambda_m`. With
  `lambda_m(A~^m) = 5`, which is not 0 in characteristic 0, this gives the decomposition. Correct.
* **(e)** Both forms of `pr` are right: the `A~^m`-form and the `calY`-form, by `A~^m = Bprod calY^(5-m)`.
* **(e), Gram determinant.** It is `5 (5 Bprod)^4`. The 4x4 block is anti-diagonal, and reversing four indices is an even
  permutation (two transpositions).
* **(g)** The Eisenstein argument after `x = 1/X` is correct. It is a remark the proof does not use.
* **Presentation.** The `m = 0` conventions (`lambda_0 = tr/Bprod`, "read `A~^m` as `A~^5`") are consistent but awkward. Indexing
  shifts by `m = 1..5`, with `gl^(5) = gl^(0)` as UR does, removes both special readings (recommendation R1, s.3).
* **Independent evidence.** `ref_proj_run.log` (s.4.1): the decomposition and the closed formula, against honest linear algebra
  at 5 random points with general weights.

### 2.2 Item 1b: (f), integrality and pole bound. Accepted

* The supports of the `calY^j`, `j = 0..4`, are the five classes `row - column = j` mod 5. They partition the 25 positions.
* The window `{2, .., j+1}` avoids index 1, so `(calB_2^(j))^-1` is a unit. Every `(calB_a^(j))^-1` is integral.
* Both directions of "`c` integral iff every `c_j` integral" follow.
* `k[[u]][calY]` is a ring, because `calY^5 = u beta^-1`.
* The pole bound: at most one `b~_1` occurs among `5 - m <= 4` consecutive weights. Correct.

### 2.3 Item 1e: `lem:gau-lift`, `lem:gau-scal` (GAPS_G2_2 s.4). Accepted

* **`lem:gau-lift`.**
  * `[A~, calZ] = Sum (rho_(a+1) - rho_a) calB_a^(m+1) E_(a,a+m+1)` is correct.
  * Solvable iff `Sum sg_a = lambda_(m+1)(X) = 0`, i.e. iff `pr X = 0`.
  * `sg_a(0) = 0` when the window `{a..a+m}` contains 1.
  * The backward formula `rho_a = rho_1 - Sum_(c=a..0) sg_c` gives `rho_a(0) = 0` for `a` in `{2-m..0}`. Hence `calZ` is integral.
  * The kernel clause comes from (c) and (f).
  * The choice `rho_1 = 0` is harmless: any other lift changes `Z` by an element of the commutant, which `lem:gau-scal` (b)
    absorbs up to the free constant.
* **Projection facts.**
  * First fact: each index occurs in exactly `j` of the five windows of length `j`. Re-derived.
  * Second fact: `xi_(a') = e_(a') (calB_(a'-j)^(j))^-1`, so the closed formula returns `(1/5) Sum_(a') e_(a')`. Re-derived.
* **`kappa^U_j` and the pivots.** `kappa^U_j = j/(5u) + h_j` with `h_j` integral. The pivots are `n + j/5`, never 0.
* **`j = 0`.** `frak d' = -(1/5) tr(dlog Theta_U) frak d`, and `c_0 = frak d (gamma - Int phi_0/frak d)`. Correct; `ker d/du = k`
  in characteristic 0.
* **Independent evidence.** `ref_proj_run.log`: both facts and part (a) for general weights, exact, with the factor `j/4` and a
  non-diagonal matrix as failing controls.

### 2.4 Item 1f: existence, `lem:gau-lead`, uniqueness clause (GAPS_G2_3 s.5.2-s.5.4). Accepted

* **Reduction 5.2.** `A Theta_U = Theta_U A~` holds entrywise. The equation for `Phat = Theta_U^-1 P` is correct at every order.
* **Existence 5.3.** (I_k), (S_k) hold by induction: `X_(k+1)` is integral with `pr X_(k+1) = 0`. Each shift component is lifted
  separately, which is legitimate because `pr` preserves the shift. (S_(k+1)) follows from `lem:gau-scal` (b) with
  `Z = calZ' + dlog calZ`.
* **`lem:gau-lead`.** The order-`z^n` equation puts `Phat_n` in the integral commutant. The order-`z^(n+1)` equation gives
  `pr(Phat_n' + dlog Phat_n) = 0`. `lem:gau-scal` (c) then gives `Phat_n = gamma frak d 1`. No hidden hypothesis: integrality of
  `P` is what makes (f) and the pivots `n + j/5` applicable.
* **Bijectivity.** Injectivity uses `P^nat_0 != 0`. The surjectivity recursion is a correct inverse-limit argument.
* **What the proof uses.** Only `P^nat` in `V(A, A~)` and `P^nat_0 = frak d Theta_U`. So any element with these two properties
  generates `V(A, A~)` freely.
* **Consequence 2 (`V(A~, A~) = k[[z]] 1`).** It follows either way. With `P^nat := 1`, directly. With the constructed
  `P^nat`: `1 = P^nat c`, `c(0) = 1`, so `c` is a unit.
* **What this delivers.** Exactly what `prop:gau-symp` step 3 needs. Invertibility of `T_0` is indeed not needed.
* **Consequence 1** (`P_0 = c(0) frak d Theta_U` for every regular gauge) is correct. So is s.11 item 7 (a regular gauge can be
  defined by `P_0 != 0`).
* **Theorem U (e)** (necessity of (EP)): correct. It even holds for any non-zero solution: the leading term `P_n` satisfies
  `A P_n = P_n A~`. Optional.
* **Independent evidence.** `ref_unique_run_11_20_2.log` (s.4.2). On random (H)+(EP) systems, two seeds, the window dimension
  is the predicted `#{k <= 10 : k = delta mod 5}` in every degree. The coordinates `k < delta` vanish, and `P_delta` spans the
  line of `frak d Theta_U`, which was computed independently from the weights. This is `lem:gau-lead` itself, beyond G2's checks.

### 2.5 Items 1g and 2: `lem:gau-degree`, graded gauges, `prop:gau-graded` step 1, s.6.4 (GAPS_G2_4 s.6). Accepted

* **Degree count.** `A E_(cd) = b_(c-1) E_(c-1,d)` and `E_(cd) A~ = b~_d E_(c,d+1)`. Both raise `k + b - a` by 1, as does
  `z d/du`. So `GEop(P^[delta]) = GEop(P)^[delta+1]` and every degree part of a solution is a solution.
* **(i)** `P^gr = (P^nat)^[0]` has `P^gr_0 = P^nat_0` (diagonal), so it is a graded regular gauge.
* **(ii)** Injectivity of `c -> P c` for `P` regular holds: `P c = P^nat c_P c`, and `k[[z]]` is a domain. With it, the freedom is
  `c(z^5)`.
* **(iii)** Only `k = 0 mod 5` feeds the entry `(1,1)`, and its constant term is `c(0) d^U_1(0) != 0`.
* **6.3** Correct: `P_can` is independent of the choice, and it is the unique regular gauge, graded or not, with
  `P_11(0,z) = 1`.
* **6.4 (i), (ii)** Re-derived:
  * `Rot_i P = Sum_delta zeta^(-i delta) P^[delta]`;
  * `T_i = zeta^i S^-i`;
  * `A T_i = zeta^i T_i A`;
  * `Rot_i(P Y) = (Rot_i P)(Rot_i Y)`;
  * UR's form is `Rot_(-1)`, and averaging gives `P^[0]`.
* **Unsplit vs split.** The referee's unsplit linear system (small box) has a null space of exactly the sum of the five degree
  parts (s.4.2).

### 2.6 Item 1h: the field (GAPS_G2_3 s.5.5). Accepted

* Every construction step is a field operation, a derivative or a `u^-1`-free integration. Only `lambda_m(A~^m) = 5` and the
  pivots `n + j/5` need characteristic 0.
* `frak d` is the binomial series of `(det Theta_U/det Theta_U(0))^(-1/5)`, rational over `Q`.
* `V_k = P^nat k[[z]]` for every `k` containing `Q` follows from (b) over `k`: the rational `P^nat` has the required `z^0` term.
* The rationality of `c` follows from the recursion `c_n = ([z^n] R^(n))_00(0)`, using `d^U_0 = 1` and `frak d(0) = 1`.
* D15 is settled on the gauge side, as claimed. Convergence is indeed never used.

### 2.7 Item 3: `lem:gau-diag` step 4, G8 (GAPS_G2_4 s.7). Accepted

* **General formula.** Re-derived: `prod_a p_(0,a) = c(0)^5 det Theta_U(0) = p_(0,1)(0)^5 prod_a rho^U_a`. Graded is not needed:
  every regular gauge has a diagonal `P_0`.
* **Under (N0).** The values at 0 are `(t_c', u b_1, t_c', -varpi_0, -varpi_0)(0) = (1, -1, 1, -1, -1)`, and
  `(v', -u v'/v, v', -v', -v')(0)` is the same vector. Both lines rest on `t_c'(0) = varpi_0(0) = v'(0) = 1` in the programme
  family (`def:set-frames`, `prop:gau-equalprod`). Hence every `p_(0,a)(0) = 1`.
* **Scope sentence, checked by hand.**
  * `varpi_0 -> w varpi_0`: the ratios at 0 become `(w, 1, w, 1/w, 1/w)`; `rho^U = (1/w, 1, 1, w, 1)` for `a = 0..4`; product 1.
  * `varpi_1 -> lambda varpi_1`: the ratios become `(1/lambda, lambda^2, 1/lambda, 1, 1)`, and the product of the `rho^U_a` is
    `lambda^5`. "Rescaling `varpi_1` changes `r` and the product": confirmed.
* **Second proof.** Correct, including GR's `p_(0,0)^5 = prod_(c<4) (b_c/b~_c)^(4-c) = v/(u varpi_0^5)`. This uses
  `v'^5 = v/(u(1-u)^5)`. The referee recomputed the exponent bookkeeping (`t_c'`: 0; `v'`: -10; `varpi_0`: -5; `(1-u)`: -10).
* **Control.** The deformed-model control `1/16` agrees with the general formula: `rho^U = (1/2, 1, 1/2, 1/2, 1/2)`.
* **Correction.** OUTLINE_4's "`= 1` for `P_can`" must carry the scope (N0), or the general value must be stated, as G2 s.11
  item 5 says.

### 2.8 Item 4: `thm:gau-G` Cor. (a), G9 (GAPS_G2_4 s.8). Accepted

* **Uniqueness.** By 5.1 (b), `P = P_can c(z)`. Then `(c(z) - 1) Y^C_Q = 0`, and the `u^0` coefficient `-e_1` (`lem:gau-vanq` (c))
  forces `c = 1`.
* **Ring.** All products live in `Q((z))[[u]]^5`. `Y^C(v(u), z)` lies in `Q[z^-1][[u]]^5`, so each `u`-order of the products is a
  finite sum.
* **Variant and Cor. (b).** Correct.
* **Control.** `g2_G.log` tests `P_can (1 + z^5/7)`, giving `-1/7` in the slot `[z^5 u^0]`. A genuine failing control.
