# Referee report `cyc_formal`: the Stage-5 critic's proof of (CYC), `critic/CYC_PROOF.md` (2026-11-14)

Folder `allgenus/stage5/referees/cyc_formal/` (new files only). AI-written by the referee; no human has read it. Labels as in the
brief: PROVED, VERIFIED (exact computation, range stated), CONJECTURAL, SKETCH, REFUTED. CYC = `critic/CYC_PROOF.md` (its s.0-s.9);
UR, GI, POLY, PB as in `STAGE5_SYNTHESIS.md`. Part 2 (`REFEREE_2.md`): C4, C5, the Theorem, the Corollary, frames, own code, errata.

**Read in full:** CYC; `critic/cr5_cyc.py`, `critic/cr5_cyc_run.py` and the three logs; `STAGE5_SYNTHESIS.md` and `_2.md`
(s.1, s.2, s.5, s.6, s.9); PB s.0-s.2.3, s.2.8, s.3 and PB_AUDIT_2 (P1, check table); `fq_conifold.py` (docstring and code);
the docstrings of `fq_rows.py` (rows, `kappa`) and `fq_ring.py`; UR s.0-s.2; GI s.0-s.1; POLY s.1-s.3. Looked at, for the frame
convention only: the docstring of `pb/pb_cyc.py` and `pb/pb_core.py` lines 245-248 (`set_frame`). No agent code was imported.

## 0. Verdict: sound-with-corrections

* **Every step of CYC s.0-s.7 is correct.** I found no mathematical error and no gap in the chain
  Lemma 0 -> C1 -> C2 -> C3 -> C4 -> C5 -> Theorem -> Corollary. I re-derived each step by hand (s.1-s.5 here, part 2 s.6-s.8).
* **The two computed inputs are now written out.** They are `gamma = 7/10` and `det Bmat = 4`. Both are finite exact computations of
  quantities the note proves to be constant. I redid both by hand from the first four Taylor coefficients (s.2.3, s.4.4). With this
  report the proof is written out in full.
* **The statement is exactly the one the programme uses.** It is H3 of PB s.1, in both directions:
  * P(b) step 5 and the glue use "cyclic system => PF5";
  * UR s.1 (the quintic rows) and PB's `Q^SP` use "PF5 => closing relation `4 -> 0`".

  The frames are exactly those of UR, GI, POLY and PB (part 2, s.9).
* **Independent evidence (VERIFIED).** My own exact code (part 2, s.10) tests each lemma separately, with failing controls. Highlights:
  * `M = PF4` coefficientwise below `u^77` in five frames;
  * the same identically in the frame parameters `(alpha, beta)` below `u^27`;
  * all 16 entries of the transversality identity;
  * `det Wr = -2/u^2`.
* **Status after this referee.** (CYC) is PROVED. It is AI-written and has been refereed once, not by a human.

**Errata** (none affects validity; part 2, s.11 gives the precise text):
* E1: the uniqueness of `g` is not covered by the cited computation.
* E2: add the hand values that fix `det Bmat`.
* E3: the frame labels in CYC s.8 use another parametrisation.
* E4: wrong file path for `fq_conifold.py`.
* E5: "every frame" should read "every conifold frame".
* E6: notation clashes.
* E7: say where isotropy enters C5.
* E8: name the ring on which the Corollary acts.

## 1. Item (1): PF4, `u` and `D` are the project's

**The project's definitions.** The file is `fq_conifold.py` at the package root.
* The brief and CYC s.0 cite `natural_form/fq_conifold.py`. No such file exists (E4). The two copies under
  `natural_form/genus5_package_review/pkg/` and `.../spot/pkg/` are byte-identical to the root file (checked with `cmp`).
* Its docstring sets `u := 1 - 3125 q` and `D = q d/dq = -(1-u) d/du`, and gives the coefficients
  `p_4 = u`, `p_3 = -2(1-u)`, `p_2 = -7/5 (1-u)`, `p_1 = -2/5 (1-u)`, `p_0 = -24/625 (1-u)`.
* The code agrees: the `pk` dictionary at lines 51-52, and `Dop(f) = -(1-u) f'` at line 94.
* These are CYC's `A_k`.

**Independent derivation (PROVED: an exact polynomial identity).** `cyc_ref_alg` A1 expands `D^4 - 5q prod_(k=1..4)(5D+k)`, with
`q = (1-u)/3125`, by operator composition and gets exactly these `A_k`.
* `D u = -(1-u) = q d/dq (1 - 3125 q)`.
* UR s.1 writes PF5 as `M^5 - 1 - q prod_(k=1..5)(5M + kz)`. Its u-form (also POLY s.1) has the coefficients
  `3, 17/5, 9/5, 274/625, 24/625`. These equal `5^(5-j) e_j(1,..,5)/3125` (A1).

**Cross-check against the project's own code (VERIFIED, exact over Q).** `cyc_ref_xproj` runs scratch copies of `fq_conifold.py`,
`fq_fp.py` and `fq_series.py` (Conifold over QRing, `Nu = 20`). The following coincide with mine:
* `f_0, f_1, f_2`, `gamma`, the log companion `h` and the zero log obstruction, below `u^28`;
* the pairing matrix in its basis `(f_0, g, f_2)`, namely `[[0,0,2],[0,0,0],[-2,0,0]]`, and its flags `Bconst_ok = pairing_vacuous = True`;
* `varpi0` and `tc` in the frames `(0,0)` and `(1/3,-2/7)`.

**Verdict on item (1): correct.**

## 2. Item (2): Lemma 0 and the Frobenius structure

### 2.1 Lemma 0: PROVED as written
`D` extends uniquely to the polynomial ring `R = K[l]` with `D l = 1 - 1/u`.

The proof runs as follows. The `l^n` coefficient of `Df` is `D a_n`, so `a_n` lies in Q. The `l^(n-1)` coefficient then forces
`-(1-u) a_(n-1)' = n a_n (1-u)/u`, i.e. `a_(n-1)' = -n a_n/u`. This is impossible, because a derivative in `Q((u))` has no `u^-1`
term. (Note that `D a` itself can have a `u^-1` term; the argument correctly divides by the unit `1-u` first.)

The essential hypothesis is that `D l` is not in `D(K)`. My control shows that it matters: with `D l' = 1 - u = D(-u)`, the element
`l' + u` is a non-constant D-constant (`cyc_ref_sol` S6).

### 2.2 Recursion and pivot: PROVED
`cyc_ref_alg` A3 computes `PF4 u^n = Sum_(s=-3..1) c_s(n) u^(n+s)` symbolically in `n`, with:
* `c_-3 = n(n-1)^2(n-2)`;
* `c_-2 = (17/5)n - (47/5)n^2 + 10n^3 - 4n^4`;
* `c_-1 = -(7/5)n + (26/5)n^2 - 6n^3 + 6n^4`;
* `c_0 = -24/625 - (2/5)n - (11/5)n^2 - 2n^3 - 4n^4`;
* `c_1 = 24/625 + (2/5)n + (7/5)n^2 + 2n^3 + n^4`.

These are polynomial identities in `n`, not finite checks. By hand, only `A_4 D^4` and `A_3 D^3` reach `u^(n-3)`, giving
`n(n-1)(n-2)(n-3) + 2n(n-1)(n-2) = n(n-1)^2(n-2)`. The values `c_s(n)` vanish whenever `n + s < 0` and `n <= 2`.

So the recursion is triangular, the coefficients of `u^0, u^1, u^2` are free, and there are no conditions at negative powers.
VERIFIED: `f_0, f_1, f_2` solve PF4 below `u^77` (S1). Their first terms are `1 + (2/625)u^3`, `u - (19/60)u^3` and
`u^2 + (37/30)u^3`, as in CYC s.8 and UR.

### 2.3 The log obstruction and `gamma = 7/10`: PROVED (here by hand)
**Set-up.** For a power-series solution `f`, `r(f) := PF4(f l) - l PF4(f) = Sum_k A_k Sum_(m>=1) C(k,m) D^(k-m) f D^m l`. The
powers of `D l` are:
* `D l = -u^-1 + 1`;
* `D^2 l = -u^-2 + u^-1`;
* `D^3 l = -2u^-3 + 3u^-2 - u^-1`;
* `D^4 l = -6u^-4 + 12u^-3 - 7u^-2 + u^-1`.

So `r(f)` lies in `u^-3 Q[[u]]`, and `g l + h` is a solution iff `g` is a power-series solution and `[u^m] r(g) = 0` for
`m = -3, -2, -1`. This is CYC s.1, correct.

**The full 3x3 obstruction matrix.** Rows are `m = -3, -2, -1`; columns are `f_0, f_1, f_2` (VERIFIED, S2):

    [[-2, 0, 0], [17/5, 0, 0], [-7/5, -7/5, 2]]

**Checked by hand.** I used `Df = -1 + u + (19/20)u^2`, `D^2 f = -1 - (9/10)u`, `D^3 f = 9/10` for `f_1`, and `Df = -2u - (17/10)u^2`,
`D^2 f = 2 + (7/5)u`, `D^3 f = -7/5` for `f_2`:
* `[u^-3] r(f) = -6 f(0) + 4 f(0) = -2 f(0)` (from `A_4` and `A_3`);
* `[u^-1] r(f_1) = 0 - 7/5 + 2 - 2 = -7/5` (the contributions of `k = 1, 2, 3, 4`);
* `[u^-1] r(f_2) = 0 + 0 + 4 - 2 = 2`;
* the `u^-2` entries of the `f_1` and `f_2` columns are `0 - 2 + 2 = 0` and `0`.

**Consequence.** The kernel is spanned by `(0, 1, 7/10)`. So the log solution is unique up to scale and up to adding power series,
and `g = f_1 + (7/10) f_2` is canonical up to scale.

**Erratum E1.** CYC's code, like `fq_conifold`, searches only `span(f_1, f_2)`. The `f_0` direction is excluded only by the `u^-3`
row. CYC s.1's "holds exactly for `g = f_1 + (7/10) f_2` (computed, `cr5_cyc_run`)" therefore claims more than the cited computation
shows. The proof itself does not use the uniqueness: it needs only the existence of `y_3`. The definition of the vanishing period
does use it.

**Controls.** `gamma = 1/2` and `g + f_0` are obstructed (S2). VERIFIED: `y_3 = g l + h` solves PF4 below `u^77`, with
`g = u + (7/10)u^2 + (41/75)u^3 + ...` and `h = -(23/360)u^3 - (6397/60000)u^4 + ...`.

## 3. Item (3): C1, self-adjointness: PROVED as written

The note's argument uses two facts:
* `(q o P(D))* = P(-D) o q = q o P(-D-1)`, which uses `D o q = q o (D+1)`;
* `P(-D-1) = P(D)`: the four sign changes cancel, and `k -> 5-k` permutes the factors.

**Checked exactly (`cyc_ref_alg` A4).**
* `D o q = q o (D+1)` and `P(-D-1) = P(D)`, as identities of operators and of polynomials.
* `PF4* = PF4` coefficientwise.
* An equivalent check that I derived. For an operator of order 4, expanding `L* = D^4 A_4 - D^3 A_3 + D^2 A_2 - D A_1 + A_0` shows that
  `L* = L` iff `A_3 = 2 D A_4` and `A_1 = D A_2 - D^3 A_4` (the other two conditions follow). For PF4 these read
  `2 D u = -2(1-u) = A_3` and `-(7/5)(1-u) + (1-u) = -(2/5)(1-u) = A_1`. Both hold.

**Controls: none of the following is self-adjoint.**
* `PF4 + uD`;
* `A_3 -> -3(1-u)`;
* `D^4 - 5 prod(5D+k) o q` (q placed on the right).

## 4. Item (4): C2, the concomitant

### 4.1 (i) The telescoping identity: PROVED as written
Put `T_j = (-1)^j D^j H D^(k-j) f`. Then `D[(-1)^j D^j H D^(k-1-j) f] = T_j - T_(j+1)`, and the sum over `j` is `T_0 - T_k`.

VERIFIED as an exact identity of jets (bilinear differential polynomials in two indeterminates, `cyc_ref_alg` A5):
* `g L f - f L* g - D P(f,g) = 0` for `L = PF4`;
* the same for a random order-4 operator with quadratic coefficients, since (i) is a general identity;
* control: without the signs `(-1)^j` the identity fails.

### 4.2 (ii) Constancy and alternation: PROVED as written
**Constancy.** `D P(f,g) = g PF4 f - f PF4* g = 0` on S. With Lemma 0 this puts `P(f, g)` in Q.

**Alternation.** The note's argument is correct:
* the coefficient of `f_(a+1) g_b` in `DQ` at top order `s + 1` is `q_ab + q_(a+1,b-1)`;
* the coefficient of `f_0 g_(s+1)` is `q_(0,s)`;
* these vanish successively, so every `q_ab` with `a + b = s` is zero.

**VERIFIED, and stronger.**
* `P(f,g) + P(g,f) = 0` holds identically as a jet (A5).
* Control: for `PF4 + uD` the sum `Q` is not zero, and neither is `DQ`.
* On S, `B` is log-free and constant below `u^78`, and alternating (`cyc_ref_sol` S3).

### 4.3 (iii) The jet matrix: PROVED (hand computation)
From the definition, `P_ab = Sum_(k=a+b+1..4) (-1)^(k-1-a) C(k-1-a, b) D^(k-1-a-b) A_k`. With `D(1-u) = 1-u` and `D u = -(1-u)`:

    P = [[0, (2/5)(1-u), 1-u, -u], [-(2/5)(1-u), 0, u, 0], [-(1-u), -u, 0, 0], [u, 0, 0, 0]]

that is, `B(f,g) = (2/5)(1-u) W_01 + (1-u) W_02 + u W_12 - u W_03`, with `W_ab = D^a f D^b g - D^b f D^a g`.

Properties:
* the matrix equals the code's exactly (A6);
* it vanishes for `a + b >= 4`, and `P_(3-j,j) = (-1)^j u`, as CYC (iii) states;
* `det P = u^4`, because the longest permutation of four letters is even (A6).

### 4.4 (iv) `Bmat`, `det Bmat = 4`, and the Wronskian: PROVED (hand values added, E2)
**The matrix.** `Bmat = [[0, -7/5, 2, -9/10], [7/5, 0, 0, -1], [-2, 0, 0, 0], [9/10, 1, 0, 0]]` in the basis `(f_0, f_1, f_2, y_3)`
(VERIFIED, S3; it equals CYC s.3 and the critic's log).

**Only four entries matter.** Since `b_12 = b_23 = 0`, the Pfaffian is `Pf = b_01 b_23 - b_02 b_13 + b_03 b_12 = -b_02 b_13`.
`B` is constant, so each entry is its `u^0` coefficient.

**The power-series entries.** At `u = 0` we have `B(f,g) = (2/5)(f Dg - g Df) + (f D^2 g - g D^2 f)`. The values
`(f, Df, D^2 f)(0)` are `(1,0,0)`, `(0,-1,-1)` and `(0,0,2)` for `f_0`, `f_1`, `f_2`. Hence `b_01 = -7/5`, `b_02 = 2` and `b_12 = 0`.

**The entries with `y_3`.** Let `Y_k` be the `l^0` part of `D^k y_3`:
* `Y_0 = h`;
* `Y_1 = g Dl + Dh = -1 + (3/10)u + ...`;
* `Y_2 = 2 Dg Dl + g D^2 l + D^2 h = u^-1 - 9/10 + ...`;
* `Y_3 = u^-2 - (7/5)u^-1 + 87/100 + ...`.

These come from `g = u + (7/10)u^2 + (41/75)u^3`, `Dg = -1 - (2/5)u`, `D^2 g = 2/5` and `h = -(23/360)u^3`. The four terms of
`B(f_1, y_3)` then contribute `0 + 1 - 1 - 1`, so `b_13 = -1`. All four terms of `B(f_2, y_3)` vanish, so `b_23 = 0`.

So `Pf = 2` and `det Bmat = 4`.

**The Wronskian, an independent route.**
* Abel's identity gives `D det Wr = -(A_3/A_4) det Wr = (2(1-u)/u) det Wr`, so `D(u^2 det Wr) = 0`.
* By Lemma 0, `det Wr = C u^-2` with `C` in Q. VERIFIED: `C = -2` (S4).
* This is consistent with `det Bmat = (det Wr)^2 det P = C^2`.
* So `det Wr` is a unit of K. The note's "`det Wr != 0`" already suffices for the Theorem, because R is a domain.

## 5. Item (5): C3, transversality: PROVED as written

The identity `P^-1 = Wr Bmat^-1 Wr^T` holds in Frac(R). The exchange-matrix argument is correct: `PJ` is upper triangular with
diagonal `P_(c,3-c)`, so `(P^-1)_ab = 0` for `a + b < 3` and `(P^-1)_03 = 1/P_30 = 1/u`. The full inverse, computed by hand, is

    P^-1 = [[0, 0, 0, 1/u], [0, 0, -1/u, -(1-u)/u^2], [0, 1/u, 0, (2/5)(1-u)/u^2], [-1/u, (1-u)/u^2, -(2/5)(1-u)/u^2, 0]]

and `P P^-1 = 1` exactly (A6).

**VERIFIED.** `Sum_ij Binv_ij D^a y_i D^b y_j = (P^-1)_ab` holds for all 16 entries below `u^75` (S5). The critic's check (c) tests
row 0 only. Control: perturbing `Binv` by `E_01 - E_10` breaks the identity.
