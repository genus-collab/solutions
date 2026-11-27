# Referee report: Stage 2, "(V) the log legs + falsifier (b)" (allgenus/stage2/logleg/, 2026-11-14)

Adversarial referee, AI-written; no human has read it. Labels as in the project: PROVED / VERIFIED / CONJECTURAL / SKETCH / REFUTED.
Conventions of `../../SYNTHESIS_STAGE01.md` s.3. Objects read: `LOGLEG_PROOF.md`, `LOGLEG_PROOF_2.md`, `FALSIFIER_B.md`, the worker's
logs (scripts only skimmed; NOTHING of the worker's or of `allgenus/tr`, `allgenus/xyswap` is imported by the referee code).
The printed arXiv:math-ph/0702045 (Eynard-Orantin, "EO" below) was consulted directly, pages 20-36 (s.4.1-5.3).

## 0. Verdict

* **Falsifier (b): SURVIVES, independently reproduced and extended.** With an independent TR engine (different basis, different
  algorithm, s.1) every pairing of `omega_(h,n)` on `x = p^4/4 + v/p`, `y = p` with at least one L-leg `p^m` and negative v-exponent
  is ZERO, mod 2^89-1 and mod 2^61-1, for all `(h,n)` with `g + n <= 6`, `g <= 4` (the worker's range, incl. (4,2), (3,3)) AND for the
  extension `g + n <= 7`, `g <= 6`: (4,3), (3,4), (2,5), (1,6), (0,7), (5,1), (5,2), (6,1). 40739 singular multisets per prime
  (19194 mixed L/log), 0 non-zero, 0 Z/5-selection violations. VERIFIED (two primes).
* **Pure log pairings = Gap-Theorem derivatives:** all 26 `(h,n)` of the extended range, incl. the new frozen `<H^2>_4 = 7/30000`
  and, new here, `<H^3>_4 = -7/3750`, `<H>_5 = -1/82500`, `<H^2>_5 = 3/27500`, `<H>_6 = -691/102375000`. VERIFIED (two primes).
* **Theorem V and the recursion (V1): the mathematics is correct; the `+` sign is right.** (V1) VERIFIED by me on 21733 tuples,
  and the variational formula VERIFIED as a TENSOR identity for 7 cases up to (4,1)->(4,2), (3,2)->(3,3); controls fail.
* **(Q1) as quoted is NOT the printed statement, but it is a correct translation of it** (s.2, finding F1): the worker's kernel is
  the OPPOSITE of EO's and EO's hypothesis carries the opposite sign; the two signs cancel. The note must say this.
* No error that changes a result was found. Six precision / citation findings (s.3), one of which (F3) is a real wording gap in
  the Corollary (homogeneity alone does not put the legs in `Q[v][p]`). Status after correction: Theorem V PROVED modulo (Q1) (as
  translated in s.2), (Q2), (Q3); Corollary ML_K PROVED modulo (Q1)-(Q3), D(iii) in the form of F3, and DNOPS (arXiv:1509.06954
  Thm 8.1) for the ancestor identification.

## 1. Independent re-derivation (referee code `rf_*.py`, `py -I -B`, stdlib)

Engine (`rf_tr.py`, `rf_rec.py`, `rf_build.py`): at each branch point `c_b` (v = 1, c_b = zeta5^b in F_P) the local coordinate
`zeta` with `x - x(c_b) = alpha_b zeta^2` (alpha_b = x''(c_b)/2, no square roots; `t(zeta)` by series reversion, checked
identically); `omega_(g,n)` stored as a SYMMETRIC tensor in the odd basis `dxi_(b,k) = Res B(.,q) zeta^-(2k+1)/(2k+1)` (principal part
`zeta^-(2k+2) dzeta` asserted); kernel `K = Sum_k dxi_(b,k)(p0) zeta^(2k)/(2 alpha (y(zeta)-y(-zeta))) / dzeta`, i.e. orientation
`Int_{sigma q}^{q}` (the worker's, opposite to EO (4-3)); `B(q,p_l)` projected on its odd part in the J-slots. The worker's engine
works in the e-basis `dp/(p-c_b)^(k+1)` with residues at `c_0` and Z/5 rotation; mine has no rotation and different data.
Built-in checks, all passed at both primes: symmetry of EVERY computed tensor (the projection is only legitimate if the result is
symmetric: 0 violations in all 26 tensors), the dimension bound `Sum k <= 3g-3+n` (attained, never exceeded), and two routes to the
log functional (`-(1/(2k+1)) [zeta^2k] t'/(c+t)` versus `Int_0^oo dp/(p-c)^(m+1) = (-1)^m/(m c^m)`): identical.
Sizes: `omega_(4,2)` 790, `omega_(3,3)` 2705, `omega_(2,4)` 3505 coefficients; the whole range to (4,2) in 1.5 s, the run with all
pairings and controls 27 s per prime; extended range to (6,1) 110 s per prime.

| check | range | primes | result | log |
|---|---|---|---|---|
| (P1) L-leg present, e < 0 => 0 | all multisets, all residues mod 5, `g+n <= 6, g <= 4` | 2^89-1, 2^61-1 | 9943 multisets, 0 non-zero | `rf_falsifier_*.log` |
| same, extended | `g+n <= 7, g <= 6` (adds (4,3),(3,4),(2,5),(1,6),(0,7),(5,1),(5,2),(6,1)) | both | 40739, 0 non-zero | `rf_extended_*.log` |
| (P2) `<(-log p)^n>_h = F_h^(n)(1)` | 26 cases, h <= 6 | both | all equal | both |
| F_h from its definition `(2-2h)F_h = Sum Res Phi omega_(h,1)` | h = 2,3,4 | both | `-1/240, 1/1008, -1/1440` = `B_2h/(2h(2h-2))` | both |
| (V1), `+` sign | 5063 / 21733 tuples (worker / extended range) | both | 0 failures | both |
| control: (V1) with `-` sign; without xi-term | same | both | 467 / 1920 failures each (must fail) | both |
| dilaton (D) (EO Thm 4.7, does not use (Q1)) | same tuples | both | 0 failures | both |
| variational formula as a TENSOR identity in the n old slots | (h,1)->(h,2), h=1..4; (h,2)->(h,3), h=1..3; all `m_i <= 10(kmax+1)+10` | both | 0 failures; controls fail on every non-zero tuple | `rf_tensor_*.log` |
| controls: off-L legs `p^-1, p^-2, p^-5, p^-10, p^-15`; a log leg based at `p = 2` (`-log(p-2)`) | every (h,n) | both | non-zero singular values in every (h,n) | both |
| the worker's quoted numbers (13), e.g. `<p^40,-log p>_4 = -868220928`, `<p^-2,p^37>_4 = -2084566248765/1024`, `<p^-5>_3 = <p^-5>_4 = 0`, the e = 0 values 2079/128, 3132675/1024, -45541532421/32768 | - | both | all reproduced | both |

Tuple counts: I enumerate unordered multisets with all residues of `Sum m` mod 5; the worker counts ordered tuples with
`Sum m = 0 mod 5`. I recomputed the worker's ordered count from its definition: 17416 (4678 mixed), as stated. The coverage is the
same (ordered tuples are permutations of multisets; the pairing is symmetric).
Why the tensor test determines the tensor: a differential with poles of order `<= 2 kmax + 2` at the five `c_b` and holomorphic
elsewhere is `N(p) dp/(p^5-1)^(2kmax+2)` with `deg N < 10(kmax+1)`; it is fixed by its first `10(kmax+1)` expansion coefficients at
infinity, i.e. by the pairings with `p^m`, `m <= 10(kmax+1)`; I used `m <= 10(kmax+1) + 10`.
"Functions of x pair to zero": in my coordinate `x = x_b + alpha zeta^2` this is identically true, so it is a consistency check of the
worker's e-basis engine, not a failing control (F4).

## 2. (Q1) against the printed EO paper

Printed statements (arXiv:math-ph/0702045):
* (4-3): `dE_q(p) = (1/2) Int_q^{q-bar} B(xi, p)`, with footnote 3 stating this is the opposite of the notation of other references;
  recursion (4-9) with `dE_q(p)/omega(q)`, `omega(q) = (y(q) - y(q-bar)) dx(q)`. Theorem 4.1 (4-20): `W_3^(0) = Res BBB/(dx dy)`.
* (5-1): `delta_Omega y|_x dx = -Omega`. Rauch (5-4): `delta B = Sum_i Res Omega B B/(dx dy)` (assumes `Omega/dy` has no pole at the
  branch points). (5-13): `Omega(p) = Int_{dOmega} B(p,q) Lambda(q)`, path not meeting small circles around branch points (footnote 8:
  hard edges excluded). **Theorem 5.1**: for `g + k > 1`, `D_Omega W_k^(g)(p_1..p_k)|_{x(p_i)} = Int_{dOmega} W_(k+1)^(g)(.., q) Lambda(q)`
  (5-14); for `g >= 1`, `D_Omega F^(g) = -Int_{dOmega} W_1^(g) Lambda` (5-15). Temperatures: s.5.3, (5-27)-(5-31), `Omega = Int_alpha^alpha' B`,
  `Lambda = 1`, alpha, alpha' two poles of `y dx`. (`D_Omega` = `delta_Omega` for a genus-0 curve, kappa irrelevant.)

Dictionary to the worker's (and my) convention (kernel orientation `q-bar -> q`, ω_{0,3} = `-Sum Res BBB/(dx dy)`, which both engines
have: my `omega_(0,3) = -1/(2 alpha) dxi^3` = `-Sum e_{b,1}^3/x''(c_b)`): each kernel changes sign, so `omega_(g,n) = (-1)^n W_n^(g)`
(2g-2+n kernels) and `F_g = -F^(g)_EO`. The worker's hypothesis `d_v(y dx)|_x = +Int_gamma B` means, in EO's normalisation,
`Omega = -d_v(y dx) = -Int_gamma B`, i.e. `Lambda_EO = -1`. Then (5-14) gives `d_v W_n = -Int W_(n+1)`, hence
`d_v omega_n = (-1)^n (-1) (-1)^(n+1) Int omega_(n+1) = +Int omega_(n+1)`; and (5-15) gives `d_v F^(g)_EO = +Int W_1 = -Int omega_1`, hence
`d_v F_g = +Int omega_(g,1)`. **So the worker's "form used" is correct**, but LOGLEG_PROOF s.2.3 says "our kernel is EO's" - FALSE
(F1) - and omits EO's sign in (5-1); the two discrepancies cancel. Hypotheses: simple branch points with `dy != 0` (yes); `Omega/dy`
regular at the branch points (`Omega = -dp/p`, `dy = dp`: yes); path `[0, oo]` between the two poles of `y dx` avoiding the branch points
(yes; this is EO's temperature variation of s.5.3, not "s.5.1" as the worker writes); `g + k > 1` (the worker's range `(h,n) != (0,0)`,
`2h-2+n >= 0` adds only (0,1), which is the definition, and (h,0), which is (5-15)). I checked the Rauch case independently:
`d_v|_x B = Sum_b e_{b,1} e_{b,1}/(5 c_b^3)` (my derivation agrees with LOGLEG s.2.4), and `Res Omega B B/(dx dy)` with `Omega = +dp/p`
gives the same: EO's (5-4) is consistent with this dictionary. Status of (Q1) as used: PROVED (quoted, EO Thm 5.1, after the
dictionary above), VERIFIED by the referee's tensor test in 7 cases.

## 3. Line-by-line findings on the PROVED claims

**F1 (citation error; result unaffected).** LOGLEG s.2.3 "our kernel is EO's (TR_MODEL s.1)" and the attribution of the "form used"
to EO Thm 5.1 without the dictionary of s.2. Fix: replace by the dictionary (s.2 above) and cite s.5.3 (5-27)-(5-30) for the
temperature variation. Also LOGLEG s.2.4 calls `-Sum Res BBB/(dx dy)` "the sign that gives Witten-Kontsevich on `x = -p^2/2`": fine,
but it is the opposite of EO's Theorem 4.1 sign; say so.

**F2 (evidence overcounted; (ii) for h >= 2 does not need (Q1)).** (E1) "`Int omega_(h,1) = d_v F_h`, h = 2,3,4" is listed as
evidence for (Q1). It is not independent of L-regularity: by the definition `(2-2h) F_h = (1/5)<p^5>_h + v <-log p>_h` and
homogeneity, `<-log p>_h - d_v F_h = -(1/(5v)) <p^5>_h`, so (E1) is equivalent to `<p^5>_h = 0`, a (P1) instance. Consequence (in the
worker's favour): Theorem V(ii) base `h >= 2` is PROVED from (Q2) + (Q3) without (Q1); the worker mentions this route only as a remark.

**F3 (gap in the wording of the Corollary).** LOGLEG_2 s.3: "by homogeneity a weight-w vector of `K_0` has the leg
`Sum_j alpha_j v^j p^(w-5j)` mod `Q[v][x]` with constant `alpha_j`, so `Phi_i in Q[v][p]`". Homogeneity does NOT exclude `j < 0`
(`v^-1 p^(w+5)` is homogeneous of weight w, and `<v^-1 p^(w+5), ..>` can be singular). What is needed is D(iii) in the form
"the leg of a weight-w vector of the model's `Q[v]`-lattice `K_0` lies in `Q[v][p] + Q(v)[x]`" (functions of x may have any
coefficients, they pair to zero), together with the X_i having coefficients in `Q[v]`. SYNTHESIS s.6 writes `K_0 = L/Q(v)[x]`; the
dimension count `floor(w/5) + 1 - [5|w]` there is the dimension of `Q[v][p]_w / constants`, so the critic's computation supports the
needed form, but the Corollary must cite it in that form rather than "by homogeneity".

**F4 (controls).** "Functions of x pair to exactly zero" is a consistency check (trivially true in a local coordinate adapted to
x), not a failing control. The falsifier is vacuous at (0,3) and (1,1) (0 singular tuples with `Sum m = 0 mod 5`), which the table
shows but "every (h,n)" phrasing in the summary obscures. The genuine failing controls (off-L legs, wrong-sign/no-xi steps) are
adequate; I added a further one (log leg based at `p = 2`: non-zero singular values in every (h,n) with n >= 2), which shows that the
endpoint of the log leg at the pole `p = 0` of `y dx` is essential.

**F5 (convention of (Q3)).** Theorem V(ii) uses `F_h = B_2h/(2h(2h-2)) v^(2-2h)` in the worker's sign convention (kernel `q-bar -> q`,
`F` by EO (4-10)), which is `-F_EO`. The note does not state in which convention arXiv:2502.18115 Cor 1.2 is formulated. The values
agree with the TR definition for h = 2, 3, 4 (both engines) and, from the pure-log pairings, for h = 5, 6 (referee). Fix: state the
convention of the quoted corollary explicitly. (I could not check the statement of arXiv:2502.18115 Cor 1.2 beyond its abstract:
it gives `F_g - F_g^vee` for all `g >= 2` via residues at the singularities of x and y, incl. logarithmic ones.)

**F6 (minor wording).** (a) "the critic's sketch had `-`": SYNTHESIS s.6 item 3 writes `C'(v) - (polynomial)` without defining the
polynomial; the `+` in (V1) is right (VERIFIED here), but "the sketch was wrong" is not established. (b) The (Q2) citation:
arXiv:2212.00320 is the x-y swap paper (universal formula; explicit formula for unramified y and rational x, which fits y = p);
I did not check that the statement number is Theorem 7.1. (c) The descendant remark (arXiv:math/0108100 s.5, no statement number) is
correctly marked as not used.

**Checked and correct (no finding).**
* s.1 (F1): poles only at branch points, zero residues, symmetry: EO Thm 4.2, (4-26), Thm 4.6; holomorphic at 0 and oo, `O(dp/p^2)`.
* s.1 (F2) homogeneity `p -> s p, v -> s^5 v`: `y dx -> s^5 y dx`, `omega -> s^(5(2-2h-n)) omega`; exponent formula. Correct.
* s.1 log leg = `Int_0^oo` = `Sum_b Res (-log p) omega`: correct (my two routes agree). Ancestor identity (X):
  `xif_1 = 1/(p^5-1)` has `xif_1(0) = -1`, `xif_1(oo) = 0`; the other `xif_a` vanish at 0 and oo; `d/dx = p^2/(p^5-v) d/dp` kills the
  value at 0 and keeps the decay at oo. PROVED (elementary), given the DNOPS expansion. The calibration paragraph of s.1 I did NOT check.
* s.2.1: `d_v p|_x = -p/(p^5-v)`, `d_v|_x = d_v|_p - Lie_xi`, `d_v(y dx)|_x = -dp/p = Int_0^oo B(., q)`. Re-derived; correct.
* s.2.2 (a), (b): integration by parts on the closed circle; division by the monic `p^5 - v`; `r/(p^5-v)` pairs to zero
  (`O(dp/p^3)` at oo); boundary terms `[i_xi w]_0^oo = 0`. Correct.
* s.2.4: `xi = Sum r_b/(p-c_b)`, `r_b = 1/(5 c_b^3)`; `d_v|_x B = Sum r_b e_{b,1} e_{b,1}`; `omega_(0,3) = -Sum e^3/x''`;
  `Int e_{b,1} = -1/c_b`; the (0,3) pairing formulas and `<(-log p)^3>_0 = 1/v`. Re-derived; correct; values reproduced.
* LOGLEG_2 s.3, (V1) proof: differentiation under fixed contours, the Lie terms (L-slot `-<(xi dPhi)_+>`, log slot 0), the
  rearrangement and the sign `+`. Correct. Proof of (i): induction on k at fixed `a >= 1`, base (Q2) and the explicit (0,3) case,
  the (0,2) exception handled. Correct. Proof of (ii): h = 1 by homogeneity + one exact number (reproduced mod two primes),
  h = 0 explicit, h >= 2 see F2. Correct. The mixed contour configuration (circle and ray intersect) is harmless: no diagonal poles.
* Corollary: the expansion over subsets S, the `S = {}` term `(-5)^(1-h) B_2h/(2h(2h-2)) = |B_2h|/(2h(2h-2) 5^(h-1)) = c_h`, frozen
  `F_0 = -(5/2) v^2 log v`. Correct, subject to F3 and to the DNOPS identification.
* A remark that is right and worth keeping: the dilaton route (D) alone gives `<J, -log p>_h` with at most a simple pole; the only
  case it cannot exclude is exponent `e = -1` (where `v <J,-log p>` has exponent 0); (Q1) is needed exactly there.

## 4. What survives

* Falsifier (b): passed; ML_K not refuted; now VERIFIED (two primes, independent engine) through `g + n <= 7`, genus <= 6.
* Theorem V (i)-(iii): PROVED modulo (Q1) [EO Thm 5.1 via the dictionary of s.2], (Q2) [arXiv:2212.00320, statement number not
  checked by me], (Q3) [arXiv:2502.18115 Cor 1.2, convention to be stated, F5]; (ii) for h >= 2 needs only (Q2), (Q3) (F2).
* Corollary ML_K: PROVED modulo (Q1)-(Q3), DNOPS (arXiv:1509.06954 Thm 8.1, eq. (2-35)) and D(iii) in the form of F3 (D(iii) itself
  VERIFIED to weight 240 only). Not checked here: the calibration paragraph of LOGLEG s.1; the identification of the TR ancestor
  correlators with the project's frozen theory beyond the 26 pure-log values (which all match).
* Required edits (worker's files are not modified by the referee): F1 (dictionary, s.5.3), F3 (D(iii) form), F5 (convention);
  F2, F4, F6 optional.

## 5. Files

`rf_tr.py` (series, local coordinate), `rf_rec.py` (dxi-basis expansions, kernel tensors), `rf_build.py` (recursion, symmetry
check), `rf_pair.py` (leg functionals, contraction, Gap values), `rf_run.py` (falsifier), `rf_checks.py` ((V1), (D), controls,
quoted numbers), `rf_tensor.py` (tensor-level variational formula). Logs: `rf_falsifier_{89,61}.log` (worker's range, 27 s each),
`rf_extended_{89,61}.log` (g + n <= 7, 110 s each), `rf_tensor_{89,61}.log` (55-68 s). Reproduce (working directory = scratch):
`py -I -B <referee>/rf_checks.py 89 12 4,2 3,3`, `... 89 17 6,1`, `py -I -B <referee>/rf_tensor.py 89 12 1,1 2,1 3,1 4,1 1,2 2,2 3,2`
(replace 89 by 61 for the second prime). No process left running; no `__pycache__`.
