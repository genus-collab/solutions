# Referee report on cluster G2, part 3: own checks (s.4), the dependency on cluster G1 (s.5), refutation attempts (s.6), files (s.7)

Continuation of `REFEREE_2.md`.

## 4. Own exact checks (VERIFIED; evidence only; no verdict above rests on them)

All scripts are stdlib Python 3.14, run as `py -I -B <script> <args>` from this folder. They import no project module and
no G2 module. Logs print no paths. Each log ends with `RESULT: ALL CHECKS PASSED`.

### 4.1 `ref_proj.py`, `ref_proj_run.py` -> `ref_proj_run.log` (args 5): `lem:gau-comm` (c)-(e), `lem:gau-scal` (a)

* **Route.**
  * General weights: random quadratic polynomials in `u`, at a random rational point `u0`.
  * Dual numbers, so `d/du` is exact and nothing is truncated.
  * (EP) is imposed by `b_4 := prod b~/(b_0 b_1 b_2 b_3)`.
  * `pr` is computed by solving `X = Sum_j c_j calY^j + [A~, W]` at `u0` (25 equations, 30 unknowns, exact over `Q`). Neither
    `lambda_m` nor `calB` is used. So this checks the closed formula itself, not only its consequences.
* **Checks (5 points, 10 each).**
  * `Theta_U` closes up (value and derivative);
  * rank `[calY^0..calY^4 | ad A~] = 25`;
  * the closed formula of (e) equals the linear-algebra `pr` on a random `X`;
  * `pr((calY^j)') = -(j/5)(Bprod'/Bprod) calY^j` and `pr(dlog Theta_U calY^j) = (tr/5) calY^j` for `j = 0..4`;
  * `pr(c' + dlog Theta_U c) = Sum_j (c_j' + kappa^U_j c_j) calY^j` for a random `c`.
* **Failing controls (all fail as they must).**
  * the factor `j/4`;
  * a non-diagonal matrix in place of `dlog Theta_U`;
  * the Frobenius-orthogonal projection in place of `pr` (the complement matters);
  * unequal products (`Theta_U` does not close).

### 4.2 `ref_unique_core.py`, `ref_unique_run.py` -> `ref_unique_run_11_20_2.log` (args 11 20 2): uniqueness clause, `lem:gau-lead`, `(A~, A~)`

* **Route.**
  * RANDOM one-pole systems mod `p = 2^61-1`, two seeds; not the quintic, not the model.
  * `b_a`, `b~_a` random units for `a != 1`; `b_1`, `b~_1` equal to `u^-1` times random units; `b~_4` fixed by (EP).
  * All equations of `z P' = A P - P A~` that involve only the unknowns `[u^m](P_k)_ab`, `k <= 11`, `m <= 20`, split by the Z/5-degree.
  * The null space is computed as an explicit basis (reduced row echelon form); then the rank of its projection onto the
    window `k <= 10`, `m <= 20 - k - 2`.
  * Own elimination code. G2's `g2_unique.py` uses a ranking formula instead and only the quintic/model pair.
* **Results (both seeds).**
  * Window dimension `#{k <= 10 : k = delta mod 5} = 3, 2, 2, 2, 2` for `delta = 0..4`, exactly as Theorem U (b)-(c) predicts.
  * `lem:gau-lead`: in degree `delta` the coordinates with `k < delta` vanish on the window. Those with `k = delta` span one line,
    which contains the vector of `frak d Theta_U`, computed independently from the weights (binomial series mod `p`).
  * `(A~, A~)`: in every degree the window projection is spanned by the `z^r 1`.
  * Unsplit system (`k <= 4`, `m <= 8`): null-space dimension `125` = the sum of the five degree parts.
* **Failing controls.**
  * **No pole** (`b_1`, `b~_1` units, (EP) kept): dimension `K = 11` in EVERY degree, i.e. a module of rank 5, not 1. This shows
    that the simple pole is what makes the gauge unique up to `c(z)` (recommendation R3).
  * **Unequal products** (`b~_4 -> (1+u) b~_4`): dimension 0 on the box `k <= 1`, `m <= 8`, while equal products give `1, 1, 0, 0, 0`
    there. (Truncation artefacts travel about five `u`-orders per `z`-order here, as G2 also observed.)
  * **Windows touching the truncation** (margins 0, 1): dimensions 5 and 3 instead of 2 in degree 1. So the method does see
    spurious solutions.

### 4.3 `ref_wkb.py`, `ref_wkb_run.py` -> `ref_wkb_run_6_10.log` (args 6 10): `lem:gau-colrows` Steps 1-3, `lem:set-hrot`

* **Route.**
  * The formal WKB amplitudes of PF5 at large radius, colour `i`, as series in `z` with coefficients in `Q(zeta)[[q]]`,
    `z^0..z^6`, `q^0..q^8`.
  * Exponent `D u_i = zeta^i L`, `L = (1 - 5^5 q)^(-1/5)`.
  * Constants of integration fixed at `q = 0` by `h_(w_i)`. These are computed from the weight multisets
    `w_i = {zeta^i - zeta^j}_(j != i) U {-5 zeta^i}`, with exact `Q(zeta)` inversion.
  * The operator is applied literally, as `(zD + zeta^i L)^5 - 1 - q prod_k (5(zD + zeta^i L) + kz)` on truncated series.
* **Checks.**
  * `h_(w_0) = 1 - (3/20) z + (9/800) z^2 - (269/80000) z^3 + ..` is rational. `Sum_(a in w_0) a^-1 = 9/5` and
    `Sum a^-3 = -126/125` reproduce its first coefficients by hand.
  * `lem:set-hrot`, computed from the weights.
  * The transport operator `T(L g) = 5 zeta^(4i) D g` and the vanishing symbol, for every colour.
  * No logarithm at any order.
  * `psi_(0,0) = L`.
  * Direct substitution: `PF5^(u_0)[psi_0] = 0` through `z^6`.
  * The covariance `psi_i(q, z) = psi_0(q, zeta^-i z)` for `i = 1..4`, in every coefficient of the range.
* **Failing controls.**
  * Colour `i` normalised by `h_(w_0)` instead of `h_(w_i)`: covariance fails (`i = 1..4`).
  * `PF5 + q z^2`, which is not invariant under `z -> zeta z`: covariance fails.
* **Scope.** This is evidence for row 0 only. The rows `a >= 1` follow by Step 4, which is algebra given `prop:set-cyclicLR`. No
  computation of the colour-`i` R-matrix rows themselves was made. FQ2 s.1 and CS s.1 record that check (verified, not derived).

### 4.4 Re-runs of the G2 scripts

`rerun_record.log`: all six logs of `GAPS_G2_6` s.10 reproduce identically from copies run in the scratch folder.

## 5. The dependency of item 5 on cluster G1 (status at the time of this report)

| input of `lem:gau-colrows` / `cor:gau-colours` | G2's statement | status now | where |
|---|---|---|---|
| (i) cyclic shape `X_(a+1) = (C_1 kappa_a)^(-1)(zD + zeta^i L) X_a`, `kappa'_a` independent of `z`, `i` | OPEN (R1, G1) | PROVED in G1 (all orders in `q`, from (M) and the genus-0 axioms); FOUND corroboration GJR18 s.5.1 (13), Guo-Ross 1703.06955 Props. 6.1-6.2; unrefereed | `GAPS_G1_2`, `GAPS_G1_3` |
| (ii) `D u_i = zeta^i L` | in the sources | derived in G1 s.2.10 (a); integration constants irrelevant (amplitudes) | `GAPS_G1_3` |
| (iii) row 0 = WKB amplitude of PF5 | (M), G1 | (M) CITED (Givental alg-geom/9603021v2, Cor. 11.7, read in the PDF by G1); WKB form of colour `i` derived (G1 Cor. 2.10 (b)) | `GAPS_G1.md` s.1, `GAPS_G1_3` |
| (iii) `psi_i(q=0) = h_(w_i)` (colour-`i` (QRR)) | `thm:set-qrr` draftnote (a) | PROVED in G1 for every colour, given the quoted (QRR) | `GAPS_G1_4` s.3.1 |
| Step 6: rows in `NN` (`cor:gau-rowsNN`) | R1 item 3 | written in G1 s.2.10 (c) | `GAPS_G1_3` |
| (iv) `lem:gau-wkb` over `Q(zeta)`, every colour | "[T] in another cluster" | to be stated as in C2; proof unchanged (transport operator `5 zeta^(4i) D` on `f/L`, s.2.10) | s03 |
| (iv) `lem:set-hrot` | one line | correct; VERIFIED from the weights (s.4.3) | s02 |

So item 5 is closed **conditionally on the acceptance of cluster G1's items 1, 2 and 3(a) by G1's referee**, and on writing
`lem:gau-wkb` over `Q(zeta)`. This referee has NOT refereed cluster G1. The table records G1's claims, read in its files, not a
verdict on them. Two consistency points between the clusters were checked:
* G1's normalisation `X_a := zeta^(-i) Q^(i)_a`, with `X_0 = L Phi_i/C_0`, is the `QQ^(i)_a = zeta^i X^(i)_a` of G2 s.9.4 Step 5.
* G1's `prop:set-cyclicLR` is stated for every colour, with coefficients `C_1 kappa_a` independent of `z` and `i`. That is
  exactly input (i).

The dependency is not new to item 5. The chain already needs `prop:set-cyclicLR` for `prop:gau-quintic`: the quintic one-pole
system IS the cyclic system of the actual rows. So `cor:gau-colours` adds no input beyond those of `prop:gau-quintic`,
`thm:gau-R0` and the quoted (M), (QRR).

## 6. Attempts to refute; what was looked for and not found

1. **Hidden use of regular semisimplicity, or of an eigenbasis over an extension field.** None. `lem:gau-comm` works
   shift by shift over `F`. (g) is a remark.
2. **Characteristic 0 used silently.** It is used in exactly three places, all named in s.5.6:
   * `lambda_m(A~^m) = 5`;
   * the pivots `n + j/5`;
   * `ker d/du = k` and the binomial series of `frak d`.
3. **Integrality.** The pole bound and the window argument of the lift are the only places where the pole enters. Without
   the pole the uniqueness clause is FALSE (rank 5, s.4.2), so (H) is genuinely needed, and the note uses it where it must.
4. **The uniqueness clause for non-graded gauges and for solutions with `P_0` not invertible.** It is covered by
   `lem:gau-lead`, with no assumption on `P_0`. VERIFIED on random systems in every degree, unsplit and split.
5. **Circularity.**
   * `lem:gau-colrows` Step 6 uses `thm:gau-R0` (iii), colour 0. That in turn uses `lem:gau-wkb`, `lem:gau-poly` and the
     colour-0 row 0 only.
   * `cor:gau-colours` uses R0 for colour 0, or for every colour if C4's simplification is not adopted.
   * `prop:set-cyclicLR` (G1) uses `thm:gau-R0` (iii) only for `cor:gau-rowsNN`.
   * No cycle.
6. **Fields.** `Q(zeta)` enters only in s.9, as claimed. `P_can` is rational. `D_i` is NOT rational for `i != 0` unless
   `D_0` lies in `Q[[z^5]]`, which is not known before `thm:gam-D` (s.2.9). The correction is real, not cosmetic.
7. **Analytic issues** (Stokes phenomena under `z -> zeta^-i z`). Not present: every statement is about formal series. The
   colour covariance is derived from formal WKB uniqueness and the `q = 0` values, not from asymptotics of the I-function, which
   would rotate the direction of `z`.
8. **G8 value 1.** It depends on (N0), i.e. on the scale of `varpi_1` (s.2.7), as G2 says. The outline's unconditional "`= 1` for
   `P_can`" is true in the programme family only.
9. **Evidence that cannot fail.** Two lines of `g2_colours.log` are tautological (R4). Every other central G2 check has a failing
   control, and so has every check of this referee.
10. **Finite checks sold as proofs.** None. G2 labels every computation VERIFIED and says no proof uses them. This was confirmed
    by reading each proof without the logs.

## 7. Files in this folder

* `REFEREE.md`: verdict, item table, items 1a-4.
* `REFEREE_2.md`: item 5, corrections.
* `REFEREE_3.md`: this file.
* `ref_proj.py`, `ref_proj_run.py`, `ref_proj_run.log`: s.4.1.
* `ref_unique_core.py`, `ref_unique_run.py`, `ref_unique_run_11_20_2.log`: s.4.2. About 2 minutes on one process.
* `ref_wkb.py`, `ref_wkb_run.py`, `ref_wkb_run_6_10.log`: s.4.3.
* `rerun_record.log`: s.4.4.

No file outside this folder was created, modified, moved or deleted. The G2 scripts were run only from copies in the
scratch folder.
