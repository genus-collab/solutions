# Referee cyc_adversary: an attempt to break (CYC) and its use (2026-11-14)

Folder `allgenus/stage5/referees/cyc_adversary/` (task key `referees/cyc_adversary`; new files only). AI-written; no human has read
it. Object: `stage5/critic/CYC_PROOF.md` (CYC), and its uses: PB = `stage5/pb/PB_AUDIT.md` s.2.3, s.2.8; POLY =
`stage4/poly/POLY.md` Lemma 2 and s.6; UR = `natural_form/theory/UNIFORM_REDUCTION.md` s.1, s.3, s.4; GI =
`allgenus/gauge/GAUGE_IDENTITY.md` s.1-s.3. Labels as in the brief. Part 2 (`REFEREE_2.md`): frames (c), deformations (d), own code
and logs, the precise list of corrections, compliance.

Read in full: CYC_PROOF.md; STAGE5_SYNTHESIS.md and _2.md; PB_AUDIT.md s.0-s.3; POLY.md; UR s.0-s.4; GI s.0-s.5;
FORMAL_QUINTIC_genus2.md s.1-s.4; the code of `fq_conifold.py` and `fq_rows.py`; in `fq_ring.py` the docstring, `W`, `C1`, the
relation and `subst`. Read for the frame conventions only: the frame constructors of `pb_cyc.py`/`pb_core.py`, `cr5_cyc_run.py`,
`gi_core.py`, `ur_gauge.py`, `pp_gauge.py`, `rf_core.py`, `rf_res.py` (gamma referee), `cr4_core.py`, `ev_core.py`.

## 0. Verdict: sound-with-corrections

I did not manage to break (CYC).
* **The theorem is true and the proof is complete.** I checked every step by hand (s.3). I also re-verified every lemma with my own
  code, exact over Q, by methods different from the critic's and pb's: a closed-form recurrence for the Frobenius basis, and
  operator composition coefficient by coefficient. Every control failed as it had to.
* **It is the statement the programme uses.** I re-derived the cyclic system from the code of fq_rows/fq_ring/fq_conifold (s.2).
  * The project's own row coefficients (`Rows.kinv`), pushed through the project's own conifold map (`Conifold.cf_images`), equal
    CYC's `1/c_a`. VERIFIED, mod 2^89-1, in two frames.
  * Both directions of the Corollary hold on objects built the way the programme builds them. (i) The WKB rows of PF5 satisfy the
    closing relation 4 -> 0. (ii) A solution built from the cyclic system alone gives a solution of PF5.
  * PF5 = z^5 D PF4 - 1 is the project's PF5 and the source of HPF. VERIFIED as an exact operator identity.
* **No hidden analytic input.** The proof uses no convergence, no actual monodromy and no Hodge theory (s.3).
* **Frames.** It covers every conifold frame the programme uses. The hypothesis `varpi_0(0) != 0` is not needed for CYC; it is
  needed downstream, by Theorem U (REFEREE_2 s.4).
* **Deformations.** Each deformation breaks exactly the lemma it should (REFEREE_2 s.5).
* **Corrections K1-K7 (REFEREE_2 s.7).** None changes a statement or a result of the programme. Summary:
  * K1: the critic's frame labels "(1/3,-2/7)", "(7/11,-3/5)" are not the programme's frames of those names.
  * K2: the scope of "every conifold frame" (wider than stated; the downstream normalisations are narrower).
  * K3: the Theorem-U shape of the coefficients does not detect a wrong `gamma`.
  * K4: the Corollary should name the ring the `Y_a` live in.
  * K5: C4 has a further, unlisted use (the `i_3` relation of fq_ring under Cf).
  * K6: the wording on the Yukawa factor 5.
  * K7: the wording of SYN S9.

## 1. Notation

`K = Q((u))`, `D = q d/dq = -(1-u) d/du`, `L^5 = 1/u`, `s = u^(1/5) = 1/L`. PF4 = `Sum A_k D^k` with `A_4 = u`, `A_3 = -2(1-u)`,
`A_2 = -(7/5)(1-u)`, `A_1 = -(2/5)(1-u)`, `A_0 = -(24/625)(1-u)`. Frobenius basis at `u = 0`: `f_0, f_1, f_2 = u^k + O(u^3)`, and
`y_3 = g log u + h` with `h_0 = h_1 = h_2 = 0`. Here `g = f_1 + gamma f_2` and `gamma = 7/10`. For a frame `(y_0, y_1)`, `t = y_1/y_0`
and `M(y_0, y_1) := (1/y_0) D (1/Dt) D (u y_0^2 Dt^2) D (1/Dt) D (1/y_0)`. The programme's frame `(alpha, beta)` is
`varpi_0 = f_0 + alpha g + beta f_2` (fq_conifold), with `y_1 = g` and `t_c = g/varpi_0`.

## 2. Task (a): is the proved statement the one used?

### 2.1 The rows, from the code
* `fq_rows.Rows` (lines 221-243) builds `X_0 = psi_i/c_0` and `X_(a+1) = kinv[a] (zD + zeta^i L) X_a` for `a = 0..3`.
  * `kinv = [c_0^2/W, W^2 L^-5/c_0^2, c_0^2/W, 1/c_0]`, with `W = c_0^2 + c_0 i_1 - c_1 i_0 = c_0^2 C1` (fq_ring lines 9-10, 211-215).
  * Hence `1/kinv[a] = (C1, L^5/(C0^2 C1^2), C1, C0) = C1 kappa_a`, with `kappa = (1, L^5/(C0^2 C1^3), 1, C0/C1, C0/C1)` (docstring).
* The closing coefficient is `C1 kappa_4 = C0`. It appears only in the docstring and in UR s.1; the code never uses it. It is also
  forced: the symbol `l^5 (1 - 5^5 q) = 1` requires `prod c_a = L^5`.
* With `Y_a = e^(u_i/z) X_a` and `D u_i = zeta^i L` this is `zD Y_a = c_a Y_(a+1)`, `a` in Z/5.

### 2.2 The conifold frame, from the code
* `fq_conifold.cf_images` sends `c_j -> D^j varpi_0` and `i_j -> D^j g - j D^(j-1) varpi_0 - Lam D^j varpi_0`.
* By hand: `Cf(W) = varpi_0^2 + varpi_0 (Dg - varpi_0 - Lam D varpi_0) - D varpi_0 (g - Lam varpi_0) = varpi_0 Dg - g D varpi_0
  = varpi_0^2 D t_c`.
  * The `c_0^2` term of `W` cancels the `-j D^(j-1)` part of `Cf(i_1)`, and the `Lam` terms cancel. A `W` without the `c_0^2` term
    would give `C1 -> D t_c - 1`; the code's `W` has it.
* Hence `Cf(C1) = D t_c`, and `Cf(c) = (Dt_c, L^5/(varpi_0^2 Dt_c^2), Dt_c, varpi_0, varpi_0)`. **This is CYC s.0's `c`, in the
  same order.**
* **VERIFIED** (`cy_proj`, (P2)), with the project's objects imported read-only. `Cf(kinv[a]) = 1/c_a` for `a = 0..3`, and
  `Cf(C1) = D t_c`, `Cf(c_0) = varpi_0`. As s-series mod 2^89-1 through s^127, frames (0,0) and (1/3,-2/7), `Lam = 0` and `Lam = 3`.
  The control `1/(5 c_1)` fails.
* UR s.1's `b_a = -c_a/(1-u)` (`th = z d/du`) agrees for all five `a`, exact over Q ((P3)).

### 2.3 Order of the rows
`c_0 = c_2` and `c_3 = c_4`, so M is palindromic. The order still matters: with `c_1` and `c_2` swapped, the closing relation of the
WKB rows fails from z^1 on, and the vanishing solution fails PF5 at `n = 6` (failing controls in `cy_rows`). The programme's order is
CYC's ((P2)).

### 2.4 z versus w, and the hatting
* CYC and its Corollary are stated for `(z, D = q d/dq)` with `z` constant. That is the form of fq_rows (`zD + zeta^i L`) and of UR
  s.1 (`th = z d/du`).
* The hatting `z = s^4 w` enters only afterwards. HPF is PF5 conjugated by `e^(u_0/z)` and rewritten through
  `(zD + L)(s^e F) = s^(e-1) A_e F`.
  * I re-derived this: at fixed `z`, `D Fh(w(u), u) = ((1-u)/u)((4/5) th_w - th_u) Fh` and `D s^e = -(e/5)(1-u) s^e/u`.
  * `z^m N^(5-m) (s^-1 Phi_0) = s^-1 u^(m-1) w^m Phi_(5-m)` with `N = zD + L`.
* The u-form of PF5 that POLY s.1 hats, `u M^5 - 1 - (1-u)(3zM^4 + (17/5)z^2M^3 + (9/5)z^3M^2 + (274/625)z^4M + (24/625)z^5)`, is
  VERIFIED exactly (`cy_run` (G2)).
* No `w` enters CYC. Nothing to correct.

### 2.5 Normalisation of t_c and of g
M is invariant under `y_0 -> mu y_0`, `g -> lambda g`: then `Dt -> (lambda/mu) Dt`, and the factors `(1/Dt)(Dt^2)(1/Dt)` and
`(1/y_0)(y_0^2)(1/y_0)` balance. It is also invariant under `y_1 -> y_1 + a y_0` (`t -> t + a`).
* VERIFIED for `(2 varpi_0, 3 g)`, `(-f_0/2, -5g)` and `(f_0, g + f_0/3)` (`cy_frames`), and on the WKB rows (`cy_rows`).
* So the scale of `g` and of `t_c` does not matter for CYC. It does matter downstream (REFEREE_2 s.4, K2).

### 2.6 The factor 5 (FORMAL_QUINTIC_genus2 s.1)
* There `Y_3 = <<H,H,H>> = 5 L^5/(C0^2 C1^3)`, with `eta(H^a, H^b) = 5 delta_(a+b = 3)`. So `kappa_1 = L^5/(C0^2 C1^3) = Y_3/5`: the
  5 lives in `eta`.
* The row coefficient `c_1 = C1 kappa_1` carries no 5, in the code ((P2)) and in CYC. CYC's "Yukawa coupling"
  `Y = 1/(u y_0^2 Dt^3)` is `kappa_1`, normalised by the Darboux form of B; it equals `Y_3^cf/5`.
* A factor 5 in `c_1` alone would give `prod c_a = 5 L^5`. The closing relation then fails already at z^0, and the vanishing solution
  fails PF5 at `n = 6` (failing controls).
* Moving constants between the `c_a` at fixed product is a constant diagonal gauge, invisible to the closing relation.
* Consistent; only the wording of CYC s.7 needs a clause (K6).

### 2.7 PF5 is the project's PF5
* `PF5 = M^5 - 1 - q prod_(k=1..5)(5M + kz)` is fq_rows' Pcal before conjugation (line 8; coefficients 3125, 9375, 10625, 5625,
  1370, 120 at line 121) and the equation of FORMAL_QUINTIC_genus2 s.1.
* It matches the I-function. `(M^5 - 1) q^(lambda/z + d) = prod_j (lambda - lambda_j + dz) q^(lambda/z + d)` with
  `lambda_j = zeta^j`, and `q prod_k (5M + kz)` on `q^(lambda/z + d - 1)` gives `prod_(k = 5d-4..5d) (5 lambda + kz)`. This is the
  ratio of the d-th and (d-1)-th terms.
* `z^5 D PF4 - 1 = PF5`: VERIFIED as an exact z-graded operator identity (`cy_run` (G1)).

### 2.8 How the programme uses the Corollary, and the test of each use
| use | direction | object it is applied to | own test (exact over Q) |
|---|---|---|---|
| UR s.1 and s.4 (`Q D = P_can Q^m`), hence R0 and the Stage-2 chain | PF5 + four forward relations => closing 4 -> 0 | WKB rows `X_a` of fq_rows after Cf | `cy_rows` (H): colour-0 WKB solution of PF5 by its own recursion; the closing relation holds at z^0..z^12, frames (0,0), (1/3,-2/7), (7/11,-3/5) and a rescaled frame; 4 controls fail |
| PB s.2.8, POLY s.6 step 5, POLY Lemma 2 for `Phi^can` | cyclic system => PF5, hence HPF | column 0 of `P_can Q^m` | `cy_rows` (G): the vanishing solution `Y^C_Q` of GI Lemma Q, built from the five `b_a` alone; `D PF4 Psi_n = Psi_(n-5)` for `n <= 16`, three frames; 4 controls fail at `n = 6` |
| `fq_ring.reduce` under Cf (identities proved in N modulo the relation, e.g. symplecticity) | C4 alone | relation `B(c_0, I_1) = 0`, whose Cf image is `B(varpi_0, g) = 0` | `cy_run` (C3); `cy_proj` (P1): the project's pairing table |

The third use is not in any input table (K5). (G) is the stronger test of the direction P(b) needs: `Y^C_Q` is never built from PF5.
The same tests also re-verify GI Lemma Q(a), (d) in three frames: the grading, `ord_u y^(n)_a = (4n + c_a + 1)/5`, and
`y^(1)_0 = -t_c`.

## 3. Task (b): line-by-line check of the proof; hidden inputs

* **s.1, Frobenius structure.** Correct.
  * My closed form `PF4 u^x = Sum_j P_j(x) u^(x+j)` gives `P_-3 = x(x-1)^2(x-2)`, `P_-2 = -4x^4 + 10x^3 - (47/5)x^2 + (17/5)x`, and
    `P_j(n) = 0` whenever `n + j < 0` and `n <= 2`. So PF4 maps `Q[[u]]` to itself.
  * The obstruction of `g log u` is `O_m = Sum_j P_j'(m-j) b_(m-j)` for `m = -3, -2, -1`.
    * `O_-3 = -2 b_0`, so `g` has no `f_0` part.
    * `O_-2` is automatic, because `P_-3'(1) = 0` (double root at 1).
    * `O_-1 = 0` gives `gamma = -P_-2'(1)/P_-3'(2) = (7/5)/2 = 7/10`.
  * This is a finite exact computation, hence PROVED. The u^-1 obstruction functional equals `B(f_0, .)` on `span(f_1, f_2)`, with
    values `(-7/5, 2)`. So `gamma = 7/10` is fixed twice, by solvability and by isotropy: `B(f_0, f_1 + gamma f_2) = -7/5 + 2 gamma`.
* **Lemma 0.** Correct. It uses characteristic 0, and the fact that `D log u = 1 - 1/u` has a u^-1 term, which no `d/du` of a Laurent
  series has.
* **C1.** Correct. Own criterion, from expanding `(-D)^k o A_k`: `L* = L` iff `A_3 = 2 D A_4` and `A_1 = D A_2 - D^3 A_4`; `A_0` and
  `A_2` are free. PF4 satisfies both (`cy_run` (B1), (B2)).
* **C2.**
  * (i) the telescoping is correct.
  * (ii) the alternation argument is a correct induction on the top total order of `P(f,g) + P(g,f)`.
  * (iii) correct.
  * (iv) `det P = sgn(reversal of 4 letters) u (-u) u (-u) = +u^4`. My Bmat equals CYC's, with Pfaffian 2 (`cy_run` (C1), (C2)).
* **C3.** Correct: `P` is anti-triangular, hence so is `P^-1`, and `(P^-1)_(0,3) = 1/u`. VERIFIED (D1).
* **C4.** Correct.
  * `tau: log u -> log u + 1` is well defined, because `R = K[l]` is a polynomial ring with `D l := 1 - 1/u`.
  * It commutes with `D`: `D o tau` and `tau o D` are tau-derivations that agree on `K` and on `l`. It preserves `S`, and it fixes `B`,
    which is rational.
  * No analytic continuation is used.
* **C5.** Correct: the Leibniz expansions, `Bmat^-1 = -Bmat`, and the `tau_k` formula.
* **s.7.** Correct: `D (1/Dt) D = Dt th^2`, the leading coefficient is `u`, and `E` kills a basis. The Wronskian matrix is
  invertible since `(det Wr)^2 = 4/u^4`.
* **Corollary.** Correct, including the composition order and `q (5D + 5) = 5 D q`.
* **Hidden inputs: none.**
  * The proof uses only `K`, `R = K[l]` and finitely many exact rational numbers: `gamma`, `Bmat`, the first coefficients. The
    constant values of `B` are read off exact truncations, which is legitimate because constancy is proved.
  * No convergence, no actual monodromy (`tau` is formal), no Hodge theory: "Griffiths transversality" C3 is matrix algebra.
  * The hypothesis `y_0(0) != 0` is never used, because `K` is a field (K2).
