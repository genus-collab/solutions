# s.1.5: the flat P-action graph sum equals the hatted outer sum, graph by graph (Stage 2, flat worker, 2026-11-14)

Folder `allgenus/stage2/flat/`. AI-written and AI-checked only; no human has read it. Labels: PROVED (argument written here, or a quoted
theorem with arXiv number), VERIFIED (exact computation; range and primes stated), CONJECTURAL, SKETCH, REFUTED. Conventions as in
`../../SYNTHESIS_STAGE01.md` s.3 (`H = d/dv`, `m5 = -43397/807975000`, "(G)", "quintic Theorem R", `K_vc`). Inputs read: MECHANISM.md,
SYNTHESIS_STAGE01.md, `theory/REDUCTION_DRAFT*.md` (s.1.5 is in part 1), `theory/referee/REFEREE.md`, `rf_g2flat.py/.log`,
`natural_form/theory/UNIFORM_REDUCTION.md` (UR) s.2-s.4, `genus3_natural.py`, `fq3_msum.py`, `fq_sser.py` (docstrings: rows, edges, leaves).

## 0. Results

1. **Theorem S15 (PROVED given (F)).** Let `g >= 2`, `Gamma` a stable leg-free graph of genus g, `P = P_can` the canonical Theorem-U
   gauge, `F_Gamma(u)` its contribution to the flat graph sum of `P.Omega^mod_{v(u)}` (REDUCTION_DRAFT s.1.5), and `G_Gamma(u)` its
   contribution to the project's hatted outer sum `R'.(M.omega^)` with `M = M_nat = M_0(tau(u) zhat)` taken exactly. Then
   `G_Gamma(u) = u^(2g-2) F_Gamma(u)`. Input (F): the WKB columns of the quintic satisfy `Q = P Q^m(v(u))` with trivial diagonal
   factor `D = 1` (UR s.4), i.e. (H) purity together with (N) `c^act = 1`. (F) is not new: it is the bracket [(H) and (N)] of the
   mechanism statement (SYNTHESIS s.4), CONJECTURAL to all orders, VERIFIED `zhat <= 24`, `u <= 17` (E4).
   **Consequence.** The separate input "s.1.5" of SYNTHESIS s.4 can be deleted: it is implied by Theorem U and (H), (N), which are
   listed anyway. Theorem R for `P_can` (flat) then IS Conjecture N (hatted), graph by graph, in every genus.
2. **VERIFIED, genus 3** (exact mod `2^89-1` and `2^61-1`; frame (0,0) at both primes, frame (1/3,-2/7) at `2^89-1`): for all 42 graphs
   the flat values equal the hatted values on `u^-4 .. u^3` (8 orders x 42 graphs = 336 coefficients per prime and frame, 112 of them
   non-zero regular coefficients of edge graphs; nothing fitted). The flat values equal the stored record `g3_nat_<p>_nat.json` on its
   valid part `u^-4..u^-1` (the stored pass used `M mod u^4`); beyond it the hatted side was recomputed with `M_nat` exact (s.6).
   All 41 edge graphs are regular in the flat frame (quintic Theorem R's conclusion at genus 3, never computed in the flat frame before).
   Genus 2: own code reproduces the seven values of `rf_g2flat.log` at both primes.
3. **Controls** (s.6): wrong dilaton convention, leaves removed, edge sign dropped, the stored gauge normalisation, and the constant-`M_0`
   split all FAIL as they must; the stored normalisation fails only in regular orders of graphs 1, 2, 5, 6, 7, 8, so the genus-3
   per-graph match is a graph-sum-level test of (N) at `z^5`.
4. **Where it can fail** (s.5): only through (F) (equivalently a non-trivial `D(z)`), through using an `M` other than `M_nat`, or through
   bookkeeping conventions (dilaton factor, leaf `z^1` part, sign of the edge quotient) that the proof fixes and the controls test.

## 1. The two sums

**1.1 Outer sums with a non-unital constant part.** Let `Omega` be a CohFT on `(W, eta, 1)` with ancestor potentials `F_h(q)`, `q` in
`W[[z]]` in dilaton-shifted form (base point `-z 1`). Let `O(z) = Sum_k O_k z^k : W -> W'` be symplectic, `O(-z)^* O(z) = 1` (adjoint for
`eta, eta'`), `O_0` invertible with `O_0^-1 1' = kappa^-1 1` for a scalar unit `kappa` (hypothesis (U)). Put
  `E_O(psi, psi') = (eta^-1 - O^-1(psi) eta'^-1 O^-1(psi')^t)/(psi + psi')`,   `q_O(z) = O^-1(z)(-z 1')`,
and for a stable graph `Gamma`: `Cont_Gamma[E, q; Omega] = |Aut Gamma|^-1 <prod_v d^(n_v) F_(h_v)(q), prod_e E>` (derivatives at `q`, one slot
per half-edge, contracted with the edge bivectors). By symplecticity `E_O` is a power series. If `O_0 = 1` this is the Givental-Teleman
R-action term of `Gamma` (edges, legs-free, translation `T(z) = z(1 - O^-1(z) 1)`, arXiv:1303.1043 s.2). By the dilaton equation
`d^n F_h(q/kappa) = kappa^(2h-2+n) d^n F_h(q)` (valid for `2h - 2 + n > 0`; the genus-1 anomaly is killed by `n >= 1`), the vertex at
`q_O = -(z/kappa) 1 + T` equals `kappa^(2h-2+n) d^n F_h(-z 1 + kappa T)`.

**1.2 The flat sum** (REDUCTION_DRAFT s.1.5, `rf_g2flat.py`, here `s15_core.py`): `F_Gamma := Cont_Gamma[E_P, q_P; Omega^mod_{v(u)}]`, with
`Omega^mod_v` the model CohFT in its flat basis, `P` acting from the model's flat space to the quintic's, `kappa = d = (P_0)_33`. Explicitly:
vertices `d^(2h-2+n) Sum_L prod_{l in L} (dT)_l / mult(L)! <labels, L>_h(v(u))`, `<prod tau_k(phi_a)>_h(v) = C v^e`, `e = 2-2h+Sum w(k,a)`.

**1.3 The hatted sum** (`fq3_msum.GraphSummer6`, `fq_sser` docstring). Rows `Q_a^(i)(z) = eta(H^a, R(z) e~_i)` in the normalised canonical
basis `e~_i` of the quintic TFT `omega^Q`; `R = R' M`; `Q' = Q M^-1`; edges `E'_ij = (delta_ij - Sum Q'_a^(i)(-psi) eta^ab Q'_b^(j)(-psi'))/(psi + psi')`;
leaves `T'^(i)(z) = z - z Q'^(i)_0(-z)` (orders `z^>=2`); vertices `(M.omega^Q)_{g,n}`. Since `R'^-1(z) = R'(-z)^*` (symplectic), the rows
formula is `E_{R'}`, and `R'_0 = Psi_Q` is unital (`G(u, 0) = 1` because `R~(0) = M(0) = 1`), so `T' = q_{R'} + z 1`. Hence
`T_Gamma := Cont_Gamma[E_{R'}, q_{R'}; M.omega^Q]` is the un-hatted outer term, and `G_Gamma = u^(2g-2) T_Gamma` (Lemma 3).

## 2. Four lemmas (all PROVED)

**Lemma 1 (frame change).** For a constant isometry `A : (W~, delta) -> (W, eta)`: `Cont_Gamma[E, q; A_* Omega~] =
Cont_Gamma[(A^-1 (x) A^-1) E, A^-1 q; Omega~]`, where `(A_* Omega~)_{g,n} = Omega~_{g,n} o (A^-1)^(x n)`. *Proof.* Every term is a full
contraction of tensors; insert `A A^-1` on every half-edge and every leaf slot. QED

**Lemma 2 (TFT dilaton invariance).** For a scalar unit `kappa` and `(kappa.Omega)_{g,n} := kappa^(2g-2+n) Omega_{g,n}` (unit `kappa^-1 1`):
`Cont_Gamma[E, q; kappa.Omega] = Cont_Gamma[E, q; Omega]`. *Proof.* `F^(kappa.Omega)_h(q) = kappa^(2h-2) F^Omega_h(kappa q)` termwise, so
`d^n F^(kappa.Omega)_h(q) = kappa^(2h-2+n) (d^n F^Omega_h)(kappa q) = (d^n F^Omega_h)(q)` by the dilaton equation of 1.1. QED

**Lemma 3 (hatting).** The hatted code uses `zhat = z/s^4` in every series, rows `Qhat_a = s^(1-c_a) Q_a`, and `xhat = x/s` in the vertex
TFT. Then `G_Gamma = s^(10(g-1)) T_Gamma = u^(2g-2) T_Gamma`. *Proof.* (i) `S = diag(s^(1-c_a))` satisfies `s^(1-c_a) s^(1-c_b) = 1` for
`a + b = 3 mod 5` (`c = (0,-1,3,2,1)`), so `Sum Q'_a eta^ab Q'_b` and the edges are unchanged; the unit row `a = 0` gets `s`, so every
leaf gets a factor `s`. (ii) Degree grading: a term of `Gamma` has cohomological degree `D_v + Sum_(h at v) m_h + Sum_(leaves) k_l =
3h_v - 3 + n_v + j_v` at each vertex (top degree). Under `z = lambda zhat`, `lambda = s^4`: vertex classes of degree `D` get `lambda^D`
(the R-action is graded by z-degree; for `M(z/s^4)` this is the homogeneity of s.3), edge coefficients `psi^m psi'^n` get `lambda^(m+n+1)`,
leaf coefficients `z^k` get `lambda^(k-1)`: in total `lambda^(Sum_v (3h_v-3+n_v) + E) = lambda^(3g-3)`. (iii) `xhat^(2h-2+n+j)` with `j` leaves
gives `s^-(2h-2+n+j)` per vertex, the leaf factors of (i) cancel the `s^-j`, leaving `s^-(2g-2)`. Total `s^(12(g-1) - 2(g-1))`. QED

**Lemma 4 (R-action and TFT rescaling commute).** `R.(kappa.omega) = kappa.(R.omega)` for `R(0) = 1`. *Proof.* The inner Givental sum of
`R.(kappa.omega)` is `Cont[E_R, R^-1(-z kappa^-1 1); kappa.omega]` with legs; by the computation of Lemma 2 each inner vertex gives
`kappa^(2h_v-2+n_v)` (legs included) and `Sum_v (2h_v - 2 + n_v) = 2g - 2 + n`. QED

## 3. Inputs

* **(I1) Model rows and homogeneity** (PROVED: UR s.2, s.4; LG_MODEL s.3; Teleman arXiv:0712.0160 for the semisimple model). The model
  rows at `v` are `Q^m_b(v; z) = sigma^(c_b-1) R^_b(z/sigma^4)`, `sigma = s tau^(-1/4)`, `R^(zeta) = R^(0) M_0(zeta)`, and
  `z/sigma^4 = tau zhat`. So `Omega^mod_{v(u)} = Psi^m_* (M_nat.omega^m)`, where `Psi^m = Q^m(v; 0)` is the model's constant isometry
  (flat <-> normalised canonical; at `v = 1` the vectors `c[a][i] = NSG[a] zeta^(i(a+1))` of the code) and `omega^m` the model's TFT at
  `v(u)` in that basis. That the quintic symbol (used by the hatted vertex) equals the model's `M_0` (used by the flat vertex) is UR s.2,
  PROVED given (H).
* **(I2) Gauge and WKB columns** (PROVED up to `D`, UR s.4 from Theorem U and uniqueness of formal WKB solutions): `Q D(z) = P Q^m(v(u))`,
  `D` a constant diagonal series, `D_i(z) = D_0(zeta^-i z)` by covariance. **(F)**: `D = 1` for `P = P_can`.
* **(I3) Symplecticity** of `P_can` (Prop. S, E2; PROVED given Theorem U): makes `E_P` a power series and `R' = P Psi^m` symplectic.
* **(I4) Units.** `R'_0 = Psi_Q` maps the unit of `omega^Q` to `1_Q`; `Psi^m` maps the unit of `omega^m` to `1_m`; `P_0^-1 1_Q = d^-1 1_m`.
  Since `R'_0 = P_0 Psi^m` ((F) at `z = 0`; this also identifies the two index sets `i`), the unit of `omega^Q` is `d^-1` times that of
  `omega^m`. A TFT diagonal in an orthonormal basis is fixed by its unit `1 = Sum_i y_i e~_i`: `omega_{g,n}(e~_i, .., e~_i) = y_i^-(2g-2+n)`.
  Hence `omega^Q = d.omega^m` (Lemma 2 notation). [Consistent with the referee's `d = tau^(1/4) varpi_0` and `xhat = varpi_0`.]

## 4. Proof of Theorem S15 (PROVED given (F); every step is one of s.2, s.3)

By (F) the hatted outer rows are `Q' = Q M_nat^-1 = P Q^m(v; z) M_nat^-1 = P Psi^m`, i.e. `R' = P o Psi^m` (z-dependence only through `P`).
Since `Psi^m` is an isometry, `R'^-1 = Psi^m^-1 P^-1` gives `E_{R'} = (Psi^m^-1 (x) Psi^m^-1) E_P` and `q_{R'} = Psi^m^-1 q_P`. Then

    T_Gamma = Cont[E_R', q_R'; M_nat.omega^Q]                 (1.3)
            = Cont[E_R', q_R'; M_nat.(d.omega^m)]             (I4)
            = Cont[E_R', q_R'; d.(M_nat.omega^m)]             (Lemma 4, M_nat(0) = 1)
            = Cont[E_R', q_R'; M_nat.omega^m]                 (Lemma 2)
            = Cont[E_P, q_P; Psi^m_*(M_nat.omega^m)]          (Lemma 1)
            = Cont[E_P, q_P; Omega^mod_{v(u)}] = F_Gamma       (I1),

and `G_Gamma = u^(2g-2) T_Gamma` (Lemma 3). Hypothesis (U) holds with `kappa = 1` for `R'` and `kappa = d` for `P`. QED
*Remarks.* (a) The "dilaton factor `d^(2h-2+n)`" and the leaf `dT` of the flat vertex are not a convention: they are 1.1 applied to `q_P`,
whose `z^1` part is `-(z/d) 1_m`; on the hatted side the `z^1` part of `q_{R'}` is exactly `-z 1` and the TFT carries `d` (Lemma 2 moves it).
(b) Only one property of `P` beyond (F) is used: symplecticity (I3). Flatness in `u` (the u-ODE) is not used, as in quintic Theorem R.
(c) The proof is genus-independent and graph-local: no total potential, no KP/Toda, no TR on any curve.

## 5. Exactly where the identification can fail

1. **(F), i.e. `D(z) != 1`.** `D` diagonal and covariant; (H) makes it a scalar series `c(z)` in `z^5`; (N) says `c = 1` for `P_can`.
   If `c != 1`, the proof gives `T_Gamma = Cont[E, q; ..]` for the gauge `P_can c(z)^-1`, a different flat sum: the graphs change
   individually. Size at genus 3 (control X5, the stored normalisation `P_can c(z)`, `c = 1 + m5 z^5 + ..`, `c` read off by `s15_core`): graphs 1, 2, 5, 6, 7, 8
   change, only in regular orders; e.g. graph 2 by `43397/517104000000` at `u^0`, graphs 1 and 5 by `+-43397/37231488000000`
   (`logs/s15_rational.log`). At genus 4 the stored normalisation makes graphs 1, 2, 3, 4, 7, 8 SINGULAR (E12). So s.1.5 is exactly as
   strong as (N): the per-graph genus-3 match F2 is an independent graph-sum-level confirmation of (N) at `z^5`.
2. **The split `M`.** The theorem holds for `M = M_nat` exactly. For another symplectic `M'` the hatted outer sum is
   `Cont[.., R'' = Q M'^-1; M'.omega^Q]`, a regrouping with the same total and different graphs (control X1: the constant `M_0` split,
   94 of 252 coefficients differ, totals agree). Truncations of `M_nat` (the stored `A3`, `M mod u^4`) agree only on the orders they control.
3. **`E_M` versus `E_out`.** Only the outer data `(E_{R'}, q_{R'})` enter; the `M`-edges live inside the vertex CohFT
   `M_nat.omega^m = Psi^m^*(Omega^mod)` (I1). The regrouping `R~.omega = R'.(M.omega)` (composition of R-actions, arXiv:1303.1043 s.2) is
   needed only to identify the TOTAL with `F_g` of the quintic, not for s.1.5; check F3 confirms it at genus 3.
4. **Dilaton shift and leaves.** Two wrong conventions are tempting and both fail: dropping the vertex factor and the `d` in the leaves
   (X2: 77 coefficients wrong) or dropping the leaves (X3: 88 wrong; this is the referee's leafless value at genus 2). The genus-1 anomaly
   of the dilaton equation never enters (genus-1 vertices have `n >= 1`). Hypothesis (U) holds because `P_0` is diagonal (grading).
5. **Edge formula.** The rows formula of the hatted code equals `E_{R'}` only for symplectic `R'` (I3). For `P_stored` symplecticity fails
   at `z^10` (invisible below genus 4). Dropping the alternating sign of the quotient (X4) breaks 103 coefficients.
6. **Vertex identification** `M_0(symbol) = M_0(LG rows)`: UR s.2, PROVED given (H), VERIFIED `zhat <= 6`; the flat code uses the LG rows
   (`rf_mla_core.lg_M0`), the hatted code the quintic symbol, so F2 also tests this at `zhat <= 6`.
7. **Computation only.** Window and precision (flat series exact through `u^(N-2g+1)`; one extra u-order is needed because `W = v/u` loses
   one, fixed in `s15_core.py`; hatted rows `prec_abs >= 5K + 4`), and `MVertex6`'s psi-deficit `<= 6` (enough for genus 3, not 4).

**Remaining gap.** None beyond (F). A human check of Lemmas 1-4 and of the identification of the code conventions (lowered rows versus
vector form, `eta = 5 x` involution) is outstanding; the conventions are pinned down by the exact genus-2 and genus-3 matches.

## 6. Verification (VERIFIED; exact mod `2^89-1` and `2^61-1`; `logs/s15_compare_*.log`)

Genus-3 index `k` = coefficient of `u^(k-4)`; graphs in Mathematica order (`fq3_graphs.outer_graphs3`, = the stored tables' order).

| check | what | 2^89-1 (0,0) | 2^61-1 (0,0) | 2^89-1 (1/3,-2/7) |
|---|---|---|---|---|
| H1 | own hatted driver, constant `M_0`, = stored `A_M0` (`g3_state_<p>.json`), 42 x 6 | PASS | PASS | PASS |
| H2 | own hatted driver, `M_nat` exact, = stored `A3` (`g3_nat_<p>_nat.json`) on its valid part `k <= 3` | PASS | PASS | PASS |
| H3 | hatted `K = 6` (prec 34) = `K = 8` (prec 44) on `k <= 5` | PASS | PASS | (K8 only) |
| F1 | **flat = stored hatted record**, 42 graphs, `k <= 3` (168 coefficients) | PASS | PASS | PASS |
| F2 | **flat = hatted (`M_nat` exact)**, 42 graphs, `u^-4..u^3` (336; 112 non-zero regular) | PASS | PASS | PASS |
| F3 | sum of flat graphs = stored point-frame total `total_id`, `k <= 5` | PASS | PASS | PASS |
| G2 | own flat code, genus 2 = `rf_g2flat.log` (7 graphs x 4 orders) | PASS | PASS | - |
| X1 | control: flat vs constant-`M_0` split | 94/252 differ (totals agree) | same | same |
| X2 | control: dilaton factor and leaf `d` dropped | 77/336 differ, 35 graphs | same | 77/336 |
| X3 | control: leaves removed | 88/336 differ | same | - |
| X4 | control: edge quotient without `(-1)^t` | 103/336 differ | same | - |
| X5 | control: stored normalisation `P_can c(z)` | 21/336 differ, graphs 1,2,5-8 | same | - |

Rational values (CRT of the two primes; display only, `logs/s15_rational.log`): vertex `u^-4(1/25200 - u/9000 + 407u^2/3780000
- 151u^3/3780000 + 226213751 u^4/62052480000000 + ..)`; first non-zero flat coefficients of edge graphs, e.g. graph 2 (loop on genus 2)
`-202813/38782800000000` at `u^0`, graph 4 `-113/35437500`, graph 9 `127/425250000`, graph 39 `512/7688671875` at `u^2`; graphs 37, 41,
42 vanish through `u^3`. Quintic Theorem R at genus 3 in the flat frame: all 41 edge graphs `O(u^0)`; vertex `= (1/25200) t_c^-4 +
regular` (the stored check `vertex_tc`, here reproduced: same four polar coefficients, `T_P = t_c` to the printed order).

## 7. Files, reuse, cost

| file | content |
|---|---|
| `s15_core.py` | flat engine: `P_can` (via `theory/referee/rf_flat_core`, which imports `natural_form/theory/ur_gauge` read-only), `E_P`, `q_P`, `d`; model correlators from `fq3_mvertex.MVertex6` fed with `rf_mla_core.lg_M0`; vertex tensors in flat labels; Laurent series with tracked precision; variants |
| `s15_flat.py` | flat graph sum for any leg-free graph (`fq3_graphs.outer_graphs`), tensor contraction over psi-decorations and flat colours |
| `s15_hat.py` | hatted outer sum with `M_nat` exact: `fq3_par.setup`, `fq3_msum.GraphSummer6`, `MVertex6` (all read-only); `tau^D`-scaled frozen vertex instead of eps-interpolation |
| `s15_compare.py`, `s15_rational.py`, `s15_run_all.sh` | comparisons H1-X5; display of rationals; the whole production (at most 4 processes) |
| `logs/` | one log per run (`s15_flat_g*_<p>_<variant>[_ab..].log`, `s15_hat_<p>_<mode>_K*.log`), `s15_compare_*.log`, `s15_rational.log`, `s15_run_all_times.log` |

JSON tables are in the session scratch folder `allgenus_s2_flat/` (residues; regenerate with `sh s15_run_all.sh`). Cost: the production `s15_run_all.sh`
(22 single-threaded runs, at most 4 at a time) took 2441 s CPU, 10.5 min wall (`logs/s15_run_all_times.log`): genus-3 flat sum 90-108 s
per run (43 s without leaves), hatted sum 167-190 s per run, genus 2 about 5 s; development runs about another 30 min CPU. Far below 3 h.
Nothing in the project outside this folder was modified; no wolframscript; Python stdlib only, `py -I -B`.
