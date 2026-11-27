# Referee report on allgenus/xyswap (XYSWAP_MODEL.md, x-y swap and the model lemma) - 2026-11-13

Adversarial referee; AI only. Own code in this folder (`rf_*.py`, pure Python 3.14 stdlib, `py -I -B`, at most two processes,
every run < 1 minute); no worker module imported. Labels: PROVED, VERIFIED (exact; field, points, orders stated), CONJECTURAL,
REFUTED. Verdicts per claim: confirmed / weakened / refuted / unverifiable.

**Overall verdict: SOUND WITH CORRECTIONS.** Every central mathematical claim (applicability of Theorem 7.1 of
arXiv:2212.00320, closed formulas, `F_g = B_2g/(2g(2g-2)) V^(2-2g)` on the whole family, the vanishing period, ML1) survives an
independent re-derivation by a different route (exact Q at generic points with rational branch points, own TR, own closed-formula
code, and a second theorem: Corollary 1.2 of arXiv:2502.18115). Corrections: one arithmetic slip, one sign error in prose, one
control claimed but absent from the files, an overstated piece of evidence, a missed reference, a hard-rule violation (author names).
The mechanism step (M1)/(M2) towards Conjecture N was not done and is correctly labelled CONJECTURAL.

## 1. What was re-done, and how (all logs in this folder)

| file / log | content | result |
|---|---|---|
| `rf_ser.py`, `rf_tr.py` | own direct TR, exact over Q, for `y = z`, `x` a Laurent polynomial whose critical points are all RATIONAL (involution by square root + Lagrange reversion, checked `x(a+s(t)) = x(a+t)` to full order; symmetry of every `omega_{g,n}` checked) | - |
| `rf_fam.py` | family `x = z^4 + c3 z^3 + c2 z^2 + c1 z + V/z` built from five rational branch points with `Sum 1/r_j = 0`; legs | - |
| `rf_closed.py` | own implementation of (7.9) of arXiv:2212.00320, SYMBOLIC in `z` (numeric parameters): `n = 1` all `g`, `n = 2` single edge with the edge expanded in `hbar`, `(0,3)` | - |
| `rf_cor12.py` | `F_g` from Corollary 1.2, eq. (1.3), of arXiv:2502.18115 (residues at the poles of `x`, `y`; no TR) | - |
| `rf_run_a_g3.log` | two generic points `r = (1,-2,3,-6,-3/2)`: `(c3,c2,c1,V) = (22/3,-22,-198,-216)`; `r = (3,6,-2,4,-4)`: `(-28/3,-32,592,2304)`; control `x = z + 9/z` | all PASS, see s.2 |
| `rf_run_b.log` | ML1 OFF the symmetric line along a path to `V = 0` (12 points, exact) | all PASS, see s.3 |
| `rf_run_c.log` | leg dictionary, string-equation sign, Euler-dilaton leg on the family, `1/z` legs, descendant chain | s.4 |
| `rf_run_d.log` | the printed formulas of XYSWAP s.2.2 on the line (`v = 4`, `1/3`) and line one-point values by `Res_inf` | all PASS |
| `rf_run_e.log` | `F_g`, `g = 2..6`, from Corollary 1.2 of arXiv:2502.18115 at five points; double-pole control; the Borel-sum identity `g = 1..10` | all PASS, control behaves |
| re-runs (scratch only) | `xy_run1 qz 2 10`, `xy_run2 qz 1`, `xy_explicit`, `xy_legs qz`, `xy_trgen p89 1 2`, `xy_trgen p61 11 3` | reproduce the worker's logs exactly (timings aside); no `__pycache__` created in the worker folder |

## 2. Theorem 7.1, closed formulas, free energies, vanishing period

* **Theorem 7.1 (arXiv:2212.00320, read: s.1.1-1.5, s.7, pp. 2-12, 46-49).** Statement: for `S = CP^1`, `y = z`, `x` an arbitrary rational
  function (standing assumption: simple critical points of `x`), initial form `-y dx`, all `omega^{(g)}_{m,0}` are given by (7.9). It is
  a corollary of Theorem 1.14 used in the reverse direction; the dual hypotheses (zeros of `dy = dz`: none; `y` regular at the zeros of
  `dx`) hold. `x = z^4 + v/z` is example (7.2) (`r = 5`) after `z -> v^(1/5) z`. **confirmed.** Minor: Remark 7.3 speaks of the
  discriminant of the Hurwitz space; `v = 0` is a boundary point where the pole of `x` disappears (degree drops), not literally the
  discriminant; nothing below depends on it (ML1 is proved for `V != 0` and extends as a polynomial).
* **Kernel convention (PROVED here).** (1.1) of arXiv:2212.00320 has kernel `(1/2) Int_z^{sigma z} B / (omega_{1,0}(sigma z) - omega_{1,0}(z))`
  with `omega_{1,0} = -y dx`, i.e. minus the kernel of arXiv:math-ph/0702045; induction gives `omega^A_{g,n} = (-1)^n omega^code_{g,n}`.
  VERIFIED with own code at two generic points (control: the opposite sign fails).
* **Closed formula = TR (VERIFIED, exact Q, generic points, own codes on both sides).** `n = 1`: `g = 1, 2, 3` at 3 points; `n = 2`: `g = 1, 2`
  (all terms of the single-edge graph, `W_12` expanded to `hbar^4`) at 2 points; `(0,3)` at 2 points; both generic points.
  The printed `omega_{1,1}`, `omega_{2,1}` and their `h`, `k` forms agree with own closed formula on the line (`v = 4`, `1/3`; run D).
* **F_g on the whole family (confirmed; PROVED by two routes mod quoted theorems).**
  - Own exact TR: `F_2`, `F_3` at both generic points and for `x = z + 9/z` equal `B_2g/(2g(2g-2)) V^(2-2g)`; the polynomial part of
    `Phi` contributes exactly 0, the `-V log z` part everything. The GUE control fixes the normalisation: the standard kernel gives
    `F_2 = -1/240 V^-2` with no tuning (the worker's "convention chosen to match `c_2 = -1/240`" is simply the standard one; flipping
    `y -> -y` consistently, i.e. kernel AND `Phi`, leaves `F_g` unchanged).
  - Corollary 1.2 of arXiv:2502.18115 (read pp. 1-5; not in the worker's literature search): for unramified `y = z`, `F_g` (g >= 2) is
    a sum of residues at the poles of `x`, `y`. Own implementation (run E): `g = 2..6` at the line `v = 4`, at `r = (1,-2,3,-6,-3/2)`, at a random rational
    point `x = z^4 + (2/7) z^3 - (5/3) z^2 + (11/2) z + 3 - (9/5)/z` (irrational branch points), and two GUE points: all equal
    `c_g V^(2-2g)`; the residue at infinity is identically 0. Control: with a double pole (`+ W/z^2`) the formula gives `F_2 = 0`, not
    `c_2 V^-2`, so the test is sensitive.
  - **Short proof from (1.3).** At `z = 0`: `1/x' = -(z^2/V)(1 + O(z^2))` (the `z^1` coefficient of `z^2 x'` vanishes), so `D_x^(m-1) z =
    (-1)^(m-1) (m-1)! V^(1-m) z^m (1 + O(z^2))`; `x^(2k) = V (2k)! z^(-2k-1) (1 + O(z^(2k+1)))`. The integrand is `(const/z)(1 + O(z^2)) dz`
    with `const` a function of `V` only and identical to the GUE one. At infinity the integrand has degree `<= -1 - (d+1)(2g-2) < -1`
    (`d = deg P`; `9 - 10g` for `d = 4`). Hence
    `F_g = F_g^GUE(V)` for every polynomial part; the classical GUE value (or the worker's Borel sum, re-checked: `Sum_r (r-1)!
    [eps^2g][t^(r+1)] e^(E_0) = -2^(2g) B_2g/(2g)` exactly for `g = 1..10`) gives the constant. The worker's own proof (weighted degree
    + contour lemma + monomial count) was re-derived step by step and is correct; the contour lemma has a two-line proof
    (`Res_a log z dPhi = -Res_a Phi dz/z`, residue theorem for `Phi dz/z`).
* **Vanishing period (confirmed).** `Int_inf^0 omega_{g,1} = (B_2g/2g) V^(1-2g)` and `Sum_a Res_a log z omega_{g,1} = Phi_g(0) - Phi_g(inf)`:
  VERIFIED exactly, `g = 1, 2, 3`, both generic points and the GUE control; closed-formula primitives `T_r` regular at 0, degree
  `5 - 10g` at infinity, `Phi^A_g(0) = -(B_2g/2g) V^(1-2g)` (`g = 1..3`, generic points).
* **F_1 = -(1/12) log v (weakened).** Rests on the variational formula with its orientation "fixed by consistency". At `g >= 2` the
  orientation is not a choice: dilaton + homogeneity give `<omega_{g,1}, -log z> = -Int_inf^0 omega_{g,1} = d_V F_g` (no free sign);
  applying the same orientation at `g = 1` gives `-(1/12) log v`, which agrees with the GUE and with the project's independent
  genus-one computation (LG_MODEL s.5.4). But `F_1` was not computed from its definition, and off the line it was not treated.

## 3. Model lemma ML1 (confirmed; PROVED mod Theorem 7.1, VERIFIED off the line)

The proof (pairing `= (-1)^n Res_inf .. Res_inf`, `omega_{g,n} = N/prod q(z_i)^K` with `N` in `Q[params][z]` because (7.9) has only
`z_i`, `q(z_i)`, `z_i - z_j` in denominators and `z_i`, `z_i - z_j` are monic) is correct. The worker's tests are on the line only,
where homogeneity reduces regularity to a sign of an exponent, and the "11 nontrivial zeros" are forced by the decay of `omega` at
infinity (weight `< 0`): they test the code, not the lemma (overstated as evidence). **New test (run B):** branch points
`(1, 2, -2/3, s, -s)`, i.e. `c3 = -28/9`, `c2 = -2S`, `c1 = 16/3 + 28S/3`, `V = 16S/3` (`S = s^2`); `V -> 0` with `c3, c1 != 0`, two branch
points colliding with the pole. 12 values of `s`, exact. All 36 `L`-correlators tested (`g = 1, 2` one-point `z^d`, `d <= 15`; eight
`(1,2)` pairs; four `(0,3)` triples) are polynomials in `S` of degree `<= floor(W/2)` (zero when `W < 0`), and their value at `S = 0`
equals direct TR of the limiting curve `x0 = z^4 - (28/9) z^3 + (16/3) z` (three branch points). Controls: `-log z` has a pole of order
1 (`g = 1`), 3 (`g = 2`) in `S`; `z^-5` order 2, 4; mixed `(z^5, z^-5)` order 2; `x^2` (function of `x`) pairs to 0 identically.

## 4. Dictionary, string equation, Euler-dilaton (run C, two generic points)

* `<omega_{1,1}, -Phi> = 5/24` (= rank/24) and `<omega_{2,1}, -Phi> = 2 F_2`: `tau_1(1) <-> -Phi` confirmed on the family.
* **String equation sign (correction).** `<omega_{g,2}, -z, Int z^m dx> = + <omega_{g,1}, z^m>` (`g = 1`: `m = 5, 7, 10`; `g = 2`: `m = 15, 16`):
  the unit leg acts by `+d/dx`, consistent with `tau_k = Int tau_(k-1) dx`. XYSWAP s.0 item 5 says the unit acts by `-d/dx` (with
  that sign its own conclusion would read `<tau_0(1) tau_1(H)> = -<tau_0(H)>`); s.3 says "`-d/dx` of `tau_1(H)` is `-log z`", wrong
  twice (`d/dx tau_1(H) = -log z`), the two errors cancelling. The numerical identity
  `<tau_0(1) tau_1(H)>_g = <tau_0(H)>_g` is confirmed (`g = 1, 2`, both points). The worker's log label `tau1H = x log z - z^4/4 + v/z`
  is MINUS `tau_1(H)` in their dictionary; the note converts the sign correctly.
* **Euler-dilaton (confirmed, extended).** On the whole family `e = -Phi_pol` (log-free part of `-Phi`): `<e>_1 = 1/8`,
  `<omega_{1,2}, e, H> = <omega_{2,2}, e, H> = 0` (two generic points). The edge decomposition `E = A1(x)A1' - (psi A1)(x)1 - 1(x)(psi' A1')`,
  `A1 = psi 1 + t H`, of UNIFORM_REDUCTION s.6.3(ii) was re-derived; `t <-> v` is consistent at genus one (the project's frozen point
  `t = -1` matches the TR point `v = -1`: `<H>_1 = -1/(12 v) = 1/12`, `<tau_1(1)>_1 = 5/24`, `<e>_1 = 1/8`, the project's `-1/8` being
  `-<e>_1` in its sign convention); the full dictionary of LG_MODEL s.5.2 was not re-derived here.
  **Arithmetic slip:** `psi A1 <-> Int(-(4/5) z^5) dx = -(16/45) z^9 + (v/5) z^4`, not `-(2/5) z^9 + (v/5) z^4` (still in `L`;
  conclusion unaffected). The chain `Int^k(-(4/5) z^5) dx` is log-free for `k <= 4` and has `log z` at `k = 5` (confirmed).
* `1/z` legs: `<omega_{0,3}, 1/z, 1/z, 1/z>` has a simple pole in `S` (`~ V^-1`), several others vanish; the note's "1/z gives v^-1"
  is true but no such value is in the worker's files (all their `1/z` pairings were zero). "Mixed correlators are singular" (report
  summary) is too strong: `(z^5, -log z)` vanishes on the whole path, and the worker's own table has regular mixed values.

## 5. Process, files and hard rules

* The control "the polynomial part of `Phi` contributes 0 and dropping `-V log z` destroys the identity" (XYSWAP s.0 item 3, report
  claim 5) has NO code or log in the folder (`xy_trgen.py` only checks the total): unverifiable from the files; its substance is
  VERIFIED here (run A). `xy_trgen` is an independent implementation by the same author with the same algorithm (modest independence).
* **Hard rule violated:** XYSWAP_MODEL.md lines 55, 65-66 give author names of arXiv:2212.00320 and arXiv:2211.08917; the acronym
  of the initials and several eponyms are used throughout. Files are ASCII, under 250 lines, stdlib only; logs reproduce.
* Literature: arXiv:2502.18115 (free energy under the x-y swap; Corollary 1.2) directly covers `F_g` for unramified `y`; it was missed.
* Toda (s.5): true for `log G(1 + v/hbar)`; for the TR/CohFT series the right-hand side is `log(v/hbar)` only up to a constant fixed by
  the genus-0 quadratic terms, which were not identified. Matrix model: sketch only (labelled so).

## 6. What remains open (CONJECTURAL, as the worker says)

(M1)/(M2) were not tested: no column of the regular gauge `P` beyond the leading Euler-dilaton block was translated into legs.
Obstruction to note: every descendant of every primary acquires `log z` after finitely many `x`-integrations (`tau_k(1)`, `k >= 1`;
`tau_2(H^4)`; `tau_3(H^3)`; `tau_k(H)` always), so (M1) needs exact cancellations at all orders, not just at the Euler-dilaton order.
Also ML1 is a statement about the UNDEFORMED family; the sketch in s.4 applies it to "the deformed curve" (after the `R'`-action),
which ML1 does not cover - in the graph-sum form each vertex is undeformed, so the needed input is exactly (M1)/(M2).

## 7. Claim-by-claim (the worker's report)

| # | claim | verdict |
|---|---|---|
| 1 | Theorem 7.1 of arXiv:2212.00320 applies verbatim; example (7.2), `r = 5`; log TR not needed | confirmed (Remark 7.3 scope at `v = 0` overstated, harmless) |
| 2 | own TR = closed formula with `omega^code = (-1)^n omega^A` | confirmed (independent exact check at generic points) |
| 3 | explicit `(0,3)`, `(1,1)`, `(1,2)`, `(2,1)` formulas | confirmed |
| 4 | `F_g = B_2g/(2g(2g-2)) V^(2-2g)`, all `g >= 2`, whole family | confirmed; PROVED by two routes mod quoted theorems |
| 5 | direct-TR `F_2`, `F_3` checks; polynomial part of `Phi` gives 0 (control) | confirmed in substance; the control is not on file |
| 6 | vanishing period `(B_2g/2g) v^(1-2g)`; `F_1 = -(1/12) log v` | period confirmed; `F_1` weakened |
| 7 | ML1 (all legs in `L` => polynomial) | confirmed; line-only tests and "nontrivial zeros" overstated, off-line test added |
| 8 | zero residue / one good leg not sufficient | confirmed (`1/z` example only in this report) |
| 9 | string reconciliation | confirmed numerically; prose sign errors |
| 10 | Euler-dilaton leg `-(4/5) z^5`, `<e,H> = 0`, `<e>_1 = 1/8`, edge legs in `L` | confirmed with a coefficient correction (`-16/45`) |
| 11 | all edges/leaves of `P` in `L` => Conjecture N | not tested by the worker; CONJECTURAL; gap noted in s.6 |
| 12 | Toda / second-difference equation | weakened (constant from the genus-0 normalisation) |
