# DIGEST_C2 (key read_C2): the reduction to the model, part 2 (Stage 2)

Design digest for the all-genus note. Covers: the model as topological recursion (TR) on the model curve; the leg dictionary LD;
D(ii); D(iii) (the lattice `K_0 = L`); L-regularity; the Gap Theorem; Theorem V (log legs); the Model Lemma on K (ML_K) for all
genera; Theorem S15 (flat = hatted, given (F)); the assembled chain to Conjecture N(g) given (F_g); Lemma C and universality.
No LaTeX. Every item carries its source (file + section) and its status in the sources.

Pieces: this file (s.0 sources, s.1 definitions); `_2` (s.2 statements, model side); `_3` (s.2 statements, flat = hatted, chain,
universality); `_4`, `_5` (s.3 proofs); `_6` (s.4 conventions, s.5 evidence); `_7` (s.6 literature, s.7 open issues).

**Source keys** (all paths relative to HVQ/allgenus/ unless stated):
* MS2 = `MECHANISM_STAGE2.md`; SYN2 = `stage2/STAGE2_SYNTHESIS.md`; SYN2b = `stage2/STAGE2_SYNTHESIS_2.md`.
* DP = `stage2/dict/DICT_PROOF.md`; DR = `stage2/dict/referee/REFEREE.md`.
* LL = `stage2/logleg/LOGLEG_PROOF.md`; LL2 = `stage2/logleg/LOGLEG_PROOF_2.md`; FB = `stage2/logleg/FALSIFIER_B.md`;
  LR = `stage2/logleg/referee/REFEREE.md`.
* S15 = `stage2/flat/S15_PROOF.md`; FR = `stage2/flat/referee/REFEREE.md`.
* BV = `stage2/bvar/BVAR.md`; BVT = `stage2/bvar/BV_TABLES.md`; BR = `stage2/bvar/referee/REFEREE.md`; BR4 = `.../REFEREE_S4.md`.
* Stage-1 cross-references (read in part): TRM = `tr/TR_MODEL.md`; XY = `xyswap/XYSWAP_MODEL.md`; XYR = `xyswap/referee/REFEREE.md`;
  MSG = `model/MODEL_SINGULAR.md`; RD = `theory/REDUCTION_DRAFT.md`; GI = `gauge/GAUGE_IDENTITY.md`; S01 = `SYNTHESIS_STAGE01.md`.
* Later stages: MS4 = `MECHANISM_STAGE4.md`; SYN4 = `stage4/STAGE4_SYNTHESIS.md`; MS5 = `MECHANISM_STAGE5.md`;
  SYN5 = `stage5/STAGE5_SYNTHESIS.md`; RS1, RS2 = `stage5/referees/REFEREES_SYNTHESIS.md`, `_2.md`.
* UR = `natural_form/theory/UNIFORM_REDUCTION.md` (cited by the sources; NOT read for this digest).
* [DW] = a check made by the digest writer (not in the project sources), always flagged as such.

Status labels are those of the sources: PROVED (written argument, or a quoted theorem with arXiv number and statement number),
VERIFIED (exact finite computation), CONJECTURAL, SKETCH, REFUTED. "quoted" = PROVED modulo a quoted theorem.

## s.0 Sources read, and what was not read

| source | sections read | remark |
|---|---|---|
| MS2 | s.1-s.6, in full | entry note of Stage 2 |
| SYN2 | s.0-s.2 in full (incl. the LD box and the input table) | the corrected statements of record |
| SYN2b | s.3-s.7 in full (falsifiers, K1-K11, O1-O27, s.6.1-s.6.5, critic record) | adjudications and corrections |
| DP | s.0-s.9 in full | D(iii) |
| DR | s.0-s.6 in full | D(ii) proof of record (s.3); Corrections 1-6 |
| LL, LL2, FB | in full | Theorem V, (V1), ML_K corollary; falsifier (b) |
| LR | s.0-s.5 in full | sign dictionary to printed EO (s.2); findings F1-F6 |
| S15, FR | in full | Theorem S15, Lemmas 1-4; corrections C1-C6 |
| BV, BVT, BR, BR4 | in full | B-variation; pole-order lemma; Lemma C |
| MS5 | s.1-s.5 in full | final assembled statement (supersedes MS4 s.2) |
| TRM | s.0-s.4.3 (lines 1-197) | hypotheses of arXiv:1509.06954 Thm 8.1; dictionary; Gap Theorem proof |
| XY | s.0 items 4-6, s.1, s.2.1, s.2.3, s.3 | ABDKS Thm 7.1 hypotheses; ML1 (= L-regularity); vanishing period |
| XYR | s.1-s.4 | Thm 7.1 statement checked; kernel convention; short proof of the Gap Theorem |
| MSG | s.1, s.5 | definitions of `Pi_i`, `K`, `K_0`, `kappa(X,k)`; statement (ML) |
| RD | s.1.1-s.1.5, s.2 (first part); headings of RD_2, RD_3 by search | Lemma 0 (frames), the flat graph sum s.1.5, ML_Lambda |
| GI | s.2-s.3 (Lemma V, Lemma Q, Lemma 1, (R1), Theorem G, Corollaries) | used by D(ii) and chain step 1-2 |
| S01 | s.3 (C9-C11), s.4, s.6, s.7 | name clashes; the Stage-1 route |
| MS4 s.2, SYN4 s.2, SYN5 s.1-s.2, RS1 s.3 (table), RS2 s.6.2 | as stated | checked that the Stage-2 chain is unchanged; (QRR), (CYC) |
| [DW] arXiv:math-ph/0702045 | printed pp. 21-25, 31-36 of the current arXiv PDF | Thm 4.1-4.7, Cor. 4.1, Thm 5.1, s.5.3 |
| [DW] arXiv:1303.1043 | pp. 19-27 of the current arXiv PDF | Remark 2.1, (6), Def. 2.2-Prop. 2.12 |

**Not read:** UR (all sections: Theorem U, s.2 symbol = `M_0`, s.4 WKB columns, Prop. S); `NATURAL_FORM.md`; `LG_MODEL.md`;
`model/MODEL_SINGULAR_2.md` (Lemma K, model Theorem R); the bodies of `theory/REDUCTION_DRAFT_2.md` (Lemmas E, T, F, quintic
Theorem R) and `_3.md`; GI s.4-s.6 (Prop. S); `gauge/GAUGE_IDENTITY_HAT.md`; all referee reports of Stage 1 except XYR; every
script and log (only their summaries in the notes); the papers arXiv:1509.06954, 2212.00320, 2502.18115 (statements taken from
TRM, XY, XYR). Chain steps 1-5 are Stage-1 results: they are recorded here as statements with their sources, not digested.

## s.1 Definitions (final form; notation as in the sources; clashes flagged in s.7)

**D1. Model curve (TR form).** On P^1 with global coordinate `p`: `x = W_v(p) = p^4/4 + v/p`, `y = p`, `B = dp1 dp2/(p1 - p2)^2`,
`v != 0`. Then `dx = (p^5 - v) dp/p^2`; five simple branch points `c_b = zeta^-b v^(1/5)` (`zeta = e^(2 pi i/5)`), `x''(c_b) = 5 c_b^2`,
`dy = dp != 0` there; `p` is a root of `T^5 - 4xT + 4v`. Weights: `p:1`, `v:5`, `x:4`, `q := v/p : 4`, `log p : 0`.
Source: LL s.1; DP s.1; TRM s.2.1 (the `kappa = 1` form). Other normalisations in the sources: `lambda = p^4 + v/p` with `v = 4u`
(DNOPS form, `kappa = 4`, TRM s.2.1); `x = z^4 + v/z`, `y = z` (ABDKS example (7.2), XY s.1).

**D2. TR correlators (engine convention).** `omega_{0,1} = y dx`, `omega_{0,2} = B`, recursion with kernel
`(1/2) Int_{sigma z}^{z} B/((y(z) - y(sigma z)) dx)` (equal to DNOPS (2-26)-(2-29)); `F_h` for `h >= 2` by
`(2 - 2h) F_h = Sum_b Res_{c_b} Phi omega_{h,1}`, `dPhi = y dx`. Source: LL s.1; TRM s.1; SYN2b K1. Relation to the printed EO
normalisation: s.4 item C1.

**D3. Pairings.** (i) L-leg `Phi in Q[v][p]`: `<.., Phi, ..> := Sum_b Res_{c_b} Phi omega = (1/2 pi i) oint_{|p|=R} Phi omega`
(`R > max |c_b|`; constants pair to 0). (ii) Log leg: `<.., -log p, ..> := Int_gamma omega(.., p', ..)`, `gamma` a ray from 0 to oo
avoiding the `c_b`; equals `Sum_b Res_{c_b}(-log p') omega` (cut along `gamma`); independent of the ray. (iii) Mixed: L-slots on the
circle, log slots on `gamma`. (iv) LD's pairing `<Phi, w> := Sum_c Res_c Phi w` (sum over the five branch points) for any leg analytic
near the branch points (any branch of `log p`). Notation `<Phi_1 .. Phi_a, (-log p)^k>_h`. Source: LL s.1; SYN2 s.2 step 7.1.

**D4. Dual descendant functions.** `xif_a = eps_a v^alpha_a p^(m_a)/(p^5 - v)`, `(m_a, eps_a, alpha_a) = (4,-1,0), (0,1,1), (1,1,0),
(2,-1,0), (3,-1,0)` for `a = 0..4`; `D := d/dx = p^2/(p^5 - v) d/dp`; `xif_{a,d} := D^d xif_a`. At `v = 1` this is the list of TRM s.2.2
(`xif_a` is the xi-function of the dual vector `H^(3-a)/5`); the powers `alpha_a` are forced by homogeneity. Source: SYN2 s.2 step 7.2.

**D5. Ancestor vector of a leg.** `q_Phi(z) := Sum_{a,d} <Phi, d xif_{a,d}> H^a z^d`. Source: LL s.1.

**D6. Flat data and the frozen model.** Flat basis `H^a` (`a in Z/5`), `eta(H^a, H^b) = 5 delta_{a+b = 3}`, `H = H^1`, `H = d/dv`,
`tau_k(H^b) := z^k H^b`, grading `c = (c_0..c_4) = (0,-1,3,2,1)` (RD s.1.1). The frozen model: the CohFT `M_0.omega^` (Givental-Teleman
CohFT of the LG model `x^4 + v/x`, form `dx`, on its Z/5 line; the Hurwitz Frobenius manifold `H_{0;1,4}` at `a = (0,0,0,0,v)`; MSG s.1,
TRM s.1). `<prod tau_{k_i}(H^{a_i})>_h(v)` = ANCESTOR correlators (CohFT integrals over `Mbar_{h,n}` with its own psi classes), homogeneous
`C v^e`, `e = 2 - 2h + Sum_i w(k_i, a_i)`, `w(k, a) = (4k + c_a - 4)/5` (RD s.1.2; other forms of the exponent in s.4 item C7).
Frozen free energies: `F_h = c_h v^(2-2h)` (`h >= 2`), `F_1 = -(1/12) log v`, `F_0 = -(5/2) v^2 log v`, `c_h = |B_2h|/(2h(2h-2) 5^(h-1))`
(SYN2 s.1). `F_h^(n)` = n-th v-derivative.

**D7. Descendant functionals, the slice K.** `Pi_i(q) := (1/5)[z^i] eta(S(-z)H^2, q(z))` (`S(z) = 1 + S_1/z + ..` the calibration,
`z d/dv S = (H *) S`, constants fixed by homogeneity; the column `S(z)H^2` is log-free), `Pi_0 = PiC`, `Pi_i(tau_k(H^b)) =
PiC(tau_(k-i)(H^b))`. Closed form: `PiC(tau_n(H^b)) = (-1)^(n+1) [z^-n] Y^C_b`, `Y^C = (I_0, -I_-1, -I_3, I_2, I_1)`,
`I_m = Res_{p=0} p^m e^(-W/z) dp` (GI Lemma V), `Y^C_b = -(1/5) eta(S(z)H^2, H^b)`. `K := {q : Pi_i(q) = 0 for all i >= 1}`;
`K_0 := K cap ker Pi_0`. Sources: MSG s.1; SYN2 M3; DR s.3. Equal to RD's `T_v (+) Q phi_1` and `T_v` (RD s.1.3; S01 C9 names K `K_vc`).
Bases: `kappa(X, k) := [S(z) X z^k]_+` (`X in {1, H^2, H^3, H^4}`, Q[v] coefficients) span `K_0`, and `K = Q tau_0(H) + K_0` (MSG s.1);
frames `frame(k, a) = z^k H^a + Sum_{j<k} f_j z^j H` (`a != 1`), unique in `K_0`, `f_j = const . v^((4(k-j)+c_a+1)/5)` (RD Lemma 0;
SYN2 M5). Example: `frame(1,0) = e := tau_1(1) + v tau_0(H)`, the Euler-dilaton vector (RD Lemma 0; DR Correction 1).

**D8. Legs.** `Phi_0^b = p, log p, v/p, p^3, p^2` (`b = 0..4`); `Phi_(k+1)^b = I(Phi_k^b)`, `I(F) := -Int F dx`, integrated termwise in `p`
with integration constant 0 (DP s.1; DR s.2: the termwise antiderivative never produces a `p^0` term); so `D Phi_(k+1) = -Phi_k`.
Weights `wt(Phi_k^b) = w_b + 4k`, `(w_0..w_4) = (1,0,4,3,2)` (`w_b = c_b + 1`). Source: DP s.1; SYN2 s.2 step 7.2.

**D9. Weight spaces and the leg map.** `S_w` := Q-span of `v^m tau_k(H^b)` (`m >= 0`) with `w_b + 4k + 5m = w`; `F := (+)_w S_w`, the
free graded `Q[v]`-module on the `tau_k(H^b)`. `Lambda(v^m tau_k(H^b)) := c_b^-1 v^m Phi_k^b`, Q-linear, with the LD constants
`c = (c_0..c_4) = (-1, 1, -1, -3, -2)`. By LD the TR leg of `s in S_w` is `-Lambda(s)`. Source: DP s.1; SYN2 M2.

**D10. Rings and modules.** `Q[p,q]` (`q = v/p`) = span of the `v^m p^j` with `j >= -m`; the leg class `L := Q[v][p] + Q(v)[x]`
(xyswap's `Q[v][p] + Q[v][x]` agrees on legs, DP Lemma 1); `M := Q[p,q] (+) log p . Q[v,x]`; `N := M/Q[v,x]`,
`N_w = (Q[p,q]/Q[v,x])_w (+) log p . Q[v,x]_w`; `d(w) := floor(w/5) + 1 - [5|w]`; `e(w) := #{(m,n) in Z_(>=0)^2 : 5m + 4n = w}
= dim Q[v,x]_w`. "Log-free" = zero coefficient of `log p`. Source: DP s.1, s.3, s.4.

**D11. Variation field.** `xi := p/(p^5 - v) d/dp`; `(xi dPhi)_+` := the quotient `Q` in `p Phi'(p) = Q (p^5 - v) + r`, `deg r < 5`; for
`Phi = p^m`: `(xi dPhi)_+ = m Sum_{q >= 0, m-5-5q >= 1} v^q p^(m-5-5q)`. Source: LL s.2.1-s.2.2.

**D12. Outer sums with a non-unital constant part (S15 s.1.1).** CohFT `Omega` on `(W, eta, 1)` with ancestor potentials `F_h(q)`,
`q in W[[z]]` dilaton-shifted (base point `-z 1`). `O(z) = Sum_k O_k z^k : W -> W'` symplectic (`O(-z)^* O(z) = 1`), `O_0` invertible,
`O_0^-1 1' = kappa^-1 1` for a scalar unit `kappa` (S15 calls this "hypothesis (U)"; rename, s.7). Edge bivector
`E_O(psi, psi') = (eta^-1 - O^-1(psi) eta'^-1 O^-1(psi')^t)/(psi + psi')` (a power series by symplecticity), leaf point
`q_O(z) = O^-1(z)(-z 1')`, `Cont_Gamma[E, q; Omega] := |Aut Gamma|^-1 <prod_v d^(n_v) F_(h_v)(q), prod_e E>`. TFT rescaling
`(kappa.Omega)_{g,n} := kappa^(2g-2+n) Omega_{g,n}` (unit `kappa^-1 1`).

**D13. The flat sum** (RD s.1.5; S15 s.1.2): `F_Gamma := Cont_Gamma[E_P, q_P; Omega^mod_{v(u)}]`, `P = P_can` (canonical Theorem-U gauge,
`P_can := P P_11(0,z)^-1`, GI Cor. (a)), `kappa = d = (P_0)_33`, `d(0) = 1`. Explicitly: vertices `d^(2h-2+n) Sum_L prod_{l in L} (dT)_l /
mult(L)! <labels, L>_h(v(u))`.

**D14. The hatted outer sum** (S15 s.1.3; FR s.1): rows `Q_a^(i)(z) = eta(H^a, R(z) e~_i)` in the normalised canonical basis `e~_i` of
the quintic TFT `omega^Q`; `R = R' M`, `Q' = Q M^-1`; edges `E'_ij = (delta_ij - Sum Q'_a^(i)(-psi) eta^ab Q'_b^(j)(-psi'))/(psi + psi')`;
leaves `T' = z y + q_{R'}` (corrected form, FR s.1; `y_i = Delta_i^(-1/2)`); vertices `(M.omega^Q)_{g,n}`; `M = M_nat := M_0(tau(u) zhat)`
taken exactly, `tau = 1/2F1(1, 4/5; 9/5; u)`, `v = u tau^(-5/4)`. Un-hatted outer term `T_Gamma := Cont_Gamma[E_{R'}, q_{R'}; M.omega^Q]`;
hatted term `G_Gamma`. Hatting: `u = s^5`, `zhat = z/s^4`, rows `Qhat_a = s^(1-c_a) Q_a`, vertex TFT `xhat = x/s` with `x = s varpi_0`.

**D15. WKB columns and (F).** `Q D(z) = P_can Q^m(v(u))` (UR s.4): `Q` the WKB columns of the formal quintic at the conifold
(normalised by the QRR constants at `q = 0`), `Q^m` the model rows, `D` a constant diagonal series, `D_i(z) = D_0(zeta^-i z)`.
(H): `D` is a scalar series `c(z)` in `z^5`; (N): `c^act = 1`; (F) = (H) + (N), i.e. `D = 1`; **(F_g): `D == 1 mod z^(3g-2)`**.
Source: SYN2 s.2 inputs; S15 (I2); SYN2b s.6.1.

**D16. Conjecture N(g)** (as used in SYN2 s.2 step 8): every stable graph with an edge in the hatted outer sum of the formal quintic
is regular at `u = 0`, and the vertex graph is `c_g t_c^(2-2g)` + regular. Engine normalisation: `G_Gamma = u^(2g-2) F_Gamma`; a graph
is regular iff `[u^0 .. u^(2g-3)] G_Gamma = 0` (BV header). `t_c` = conifold period, `[z^-1](Y^C_Q)_0 = -t_c` (GI Lemma Q(d)).

**D17. One-pole systems (universality).** `th Y_a = b_a Y_(a+1)`, `th = z d/du`, `b = (f, B/(f^2 g^2), f, g, g)` (self-dual:
`b_(a-1) = b_(3-a)`), `f = t_c' F(u)`, `g = -(varpi_0/(1-u)) Gm(u)`, `B = -m^20 (Pn/Pd)(u)/(u (1-u)^5)`; model `bmod = (1, -1/v, 1, -1, -1)`;
pull-back condition `v'^5 = -B v`; `1/tau = Sum_n 4 beta5_n/(5n+4) u^n`, `beta5 = (-u B)^(1/5)`; `T'' = c0 Int_0^u f`, `c0 = P_11(0,0)`.
Source: BV s.1; BR s.0, s.3.
