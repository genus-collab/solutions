# DIGEST_C2, part 3: s.2 Statements, continued (log legs, ML_K, flat = hatted, the chain, universality)

Continues `DIGEST_C2_2.md`. Same format.

## s.2.3 (cont.) Log legs and the Model Lemma

**`prop:mlk-V1` (prop) The recursion (V1).** For `h >= 0`, `a, k >= 0`, `n = a + k`, `2h - 2 + n > 0`, legs `Phi_i in Q[v][p]`:
`<Phi_1..Phi_a, (-log p)^(k+1)>_h = d_v <Phi_1..Phi_a, (-log p)^k>_h - Sum_i <.., d_v Phi_i, .., (-log p)^k>_h
  + Sum_{i=1..a} <.., (xi dPhi_i)_+, .., (-log p)^k>_h`.
* Source: LL2 s.3; LL s.0 item (iii). Status: PROVED given (Q1); VERIFIED (32559 tuples worker, 21733 referee; wrong-sign and
  no-xi controls fail). LR: "correct; the `+` sign is right".
* Corrections: K2/O10/LR F6(a): the critic's Stage-1 sketch `C'(v) - (polynomial)` (S01 s.6 item 3) agrees with (V1) once the polynomial
  is `Sum <.., Phi_i Lie_xi, ..>`, since `<Phi, Lie_xi w> = -<(xi dPhi)_+, w>`; the `+` form with `(xi dPhi)_+` is the precise statement.
  LL2 adds "(h,n) != (0,2)", which is implied by `2h - 2 + n > 0`.
* Depends on: `lem:mlk-fixedx`, `lem:mlk-ibp`, `thm:mlk-var`, (F1) properties of `omega` (poles only at the branch points, zero residues,
  symmetric, holomorphic at 0 and oo with `O(dp/p^2)` at oo: EO Thm 4.2, Cor. 4.1 (4-26), Thm 4.6, LR s.3).

**`thm:mlk-logleg` (thm) Theorem V (M9).** For `h >= 0`, `a, k >= 0`, `n = a + k`, `2h - 2 + n > 0` and legs `Phi_1..Phi_a in Q[v][p]`:
(i) if `a >= 1`, `<Phi_1 .. Phi_a, (-log p)^k>_h in Q[v]`;
(ii) if `a = 0`, `<(-log p)^k>_h = (F_h^TR)^(k)(v)` with `F_h^TR = B_2h/(2h(2h-2)) v^(2-2h)` (`h >= 2`), `F_1^TR = -(1/12) log v` (`k >= 1`),
     `F_0^TR = (1/2) v^2 log v` (`k >= 3`);
(iii) (V1) holds, with `+` in front of the xi-term.
* Source: LL s.0, LL2 s.3; SYN2 M9. Status: PROVED mod (Q1) = arXiv:math-ph/0702045 Thm 5.1 through the LR s.2 dictionary,
  (Q2) = `thm:mlk-lreg`, (Q3) = `thm:mlk-gap`. Base `h >= 2` of (ii): only (Q2) + (Q3) (LR F2). Base `h = 0`: the explicit (0,3) of
  `lem:mlk-genus0`. Base `h = 1`: homogeneity plus one exact number `-1/12` computed over `Q(zeta5)` (LL2 s.3) -- a finite computation
  inside a proof (s.7 item G1). VERIFIED (falsifier (b); tensor checks; (V1) tuples).
* Corrections: LR F1 (dictionary, s.5.3), F2 (base `h >= 2` without (Q1); evidence (E1) is equivalent to `<p^5>_h = 0`: O8), F5 (convention
  of (Q3): O13), F6; O11 (the status line must list arXiv:1509.06954 Thm 8.1; D(iii), LD and M5 are now proved as well).

**`thm:mlk-main` (thm) The Model Lemma on K (ML_K), every genus.** (Insertions) Let `h >= 0`, `n >= 1`, `2h - 2 + n > 0`, and let
`X_1..X_n in K` be homogeneous with coefficients in `Q[v]`. Then `<X_1 .. X_n>_h(v) - F_h^(n)(v) prod_i Pi_0(X_i)` lies in `Q[v]`.
(Vertex graph) `F_h = c_h v^(2-2h)` for `h >= 2`. (Frozen `F_h` of D6.)
* Equivalent forms (MSG s.5): (ML1) `<H^n>_h = F_h^(n)(v)`; (ML2) a correlator of K-generators with at least one `kappa(X,k)` has
  v-exponent `>= 0` or vanishes; function form: for `q in -z + K`, `Fbar_h(q) = F_h(Pi_0(q)) + Reg_h(q)`, all Taylor coefficients of
  `Reg_h` at `-z` regular. RD s.2 states the same as ML_Lambda(h) on `T_v (+) Q phi_1` (S01 C9: names `ML1/ML2` clash with XY's).
* Source: SYN2 s.1 (statement and table M1-M10), s.2 step 7; LL2 s.3 Corollary; MSG s.5. Status: PROVED for every `h >= 0`
  modulo arXiv:1509.06954 Thm 8.1, arXiv:2212.00320 Thm 7.1, arXiv:math-ph/0702045 Thm 5.1 with the s.4 properties of TR, and optionally
  arXiv:2502.18115 Cor. 1.2 (the Borel route needs only arXiv:2212.00320 Thm 7.1). "No non-quoted model-side input remains, and nothing
  about the quintic enters" (SYN2 s.1). VERIFIED: genus 0-6 in the leg basis; genus 0-4 in the K basis (1395 singular, all zero).
* Corrections: LR F3/O9 ("by homogeneity ... `Phi_i in Q[v][p]`" is a gap: homogeneity allows `v^-j p^(w+5j)`; closed by D(iii)(b));
  O5, O27, DR Correction 5 (the step also needs M1, LD, M6 and the `Q[v]` coefficients M5); K11.
* Depends on: `thm:trm-ident`, `prop:trm-legdict`, `prop:trm-dii`, `thm:trm-lattice`(b) with `lem:trm-classL`(a), `cor:trm-Ksplit`,
  `prop:trm-coeff`, `lem:trm-xfun`, `thm:mlk-lreg`, `thm:mlk-gap`, `thm:mlk-logleg`.

**`rmk:mlk-dilaton` (rmk) A partial route without (Q1).** The TR dilaton equation (EO Thm 4.7) with `Phi_ydx = p^5/5 - v log p` gives
`v <J, -log p>_h = (2-2h-n) <J>_h - (1/5) <J, p^5>_h` for any legs `J` (VERIFIED (D)); with (i) at `k` log legs, `<J, -log p>_h` has at most
a simple pole at `v = 0`; (Q1) is needed exactly to exclude exponent `-1`. The action of the Euler-dilaton leg as the polynomial part of
`-(4/5) x d/dx`: SKETCH. Source: LL2 s.5; LR s.3 ("right and worth keeping"); SYN2b s.6.4 (a (Q1)-free proof is optional further work:
log-TR swap arXiv:2312.16950, hypermaps/KP).

## s.2.4 Flat versus hatted (Theorem S15); all from S15, refereed by FR ("CONFIRMED")

**`lem:flh-frame` (lem) Frame change (S15 Lemma 1).** For a constant isometry `A : (W~, delta) -> (W, eta)`:
`Cont_Gamma[E, q; A_* Omega~] = Cont_Gamma[(A^-1 (x) A^-1) E, A^-1 q; Omega~]`, `(A_* Omega~)_{g,n} := Omega~_{g,n} o (A^-1)^(x n)`. PROVED.

**`lem:flh-dil` (lem) TFT dilaton invariance (S15 Lemma 2).** For a scalar unit `kappa`: `Cont_Gamma[E, q; kappa.Omega] =
Cont_Gamma[E, q; Omega]`. Read (FR s.1): `q` is absolute, i.e. each CohFT's potential is evaluated at `t = q + z . (its own unit)`; then
`F^(kappa.Omega)_h(q) = kappa^(2h-2) F^Omega_h(kappa q) = F^Omega_h(q)`. Uses the dilaton identity `d^n F_h(q/kappa) = kappa^(2h-2+n) d^n F_h(q)`
(`2h - 2 + n > 0`; at `h = 1` the anomaly is killed by `n >= 1`). PROVED.

**`lem:flh-hat` (lem) Hatting (S15 Lemma 3).** With `zhat = z/s^4` in every series, rows `Qhat_a = s^(1-c_a) Q_a`, vertex TFT `xhat = x/s`:
`G_Gamma = s^(10(g-1)) T_Gamma = u^(2g-2) T_Gamma` for every stable leg-free genus-g graph. PROVED "given the description of the hatted
code" (FR s.1: it matches the docstrings; FR: "a statement about the code, not about the mathematics; the genus-2/3 matches test it").
Corrections: C5/O18 (reference "s.3" should be (I1); `x` is never defined: `x = s varpi_0`). For the note: s.7 item G4.

**`lem:flh-comm` (lem) R-action and TFT rescaling commute (S15 Lemma 4).** `R.(kappa.omega) = kappa.(R.omega)` for `R(0) = 1`. PROVED
(FR: the translation for `kappa.omega` is `kappa^-1 T`, arXiv:1303.1043 Def. 2.5).

**`lem:flh-inputs` (lem) Inputs (I1)-(I4) (S15 s.3; FR s.1 CONFIRMED).**
(I1) Model rows: `Q^m_b(v; z) = sigma^(c_b - 1) R^_b(z/sigma^4)`, `sigma = s tau^(-1/4)`, `R^(zeta) = R^(0) M_0(zeta)`, `z/sigma^4 = tau zhat`;
  hence `Omega^mod_{v(u)} = Psi^m_*(M_nat.omega^m)`, `Psi^m = Q^m(v; 0)` (constant isometry, flat <-> normalised canonical), `omega^m` the
  model TFT at `v(u)`. The quintic symbol (hatted vertex) equals the model's `M_0` (flat vertex): UR s.2, PROVED given (H).
(I2) `Q D(z) = P Q^m(v(u))`, `D` constant diagonal, `D_i(z) = D_0(zeta^-i z)` (UR s.4, from Theorem U and uniqueness of formal WKB
  solutions); (F): `D = 1` for `P = P_can`.
(I3) `P_can` is symplectic (Prop. S, given Theorem U): `E_P` is a power series and `R' = P Psi^m` is symplectic.
(I4) `R'_0 = P_0 Psi^m` ((F) at `z = 0`; this identifies the index sets), `P_0^-1 1_Q = d^-1 1_m`, and a TFT diagonal in an orthonormal
  basis is fixed by its unit (`omega_{g,n}(e~_i, .., e~_i) = y_i^-(2g-2+n)`); hence `omega^Q = d.omega^m`
  (consistent with `d = tau^(1/4) varpi_0`, `xhat = varpi_0`; `d = 1 - u/9 - 23u^2/567 ..`, FR).

**`thm:flh-main` (thm) Theorem S15 (flat = hatted).** Let `g >= 2`, `Gamma` a stable leg-free graph of genus `g`, `P = P_can`,
`F_Gamma(u)` its term in the flat sum of `P.Omega^mod_{v(u)}` (D13), `G_Gamma(u)` its term in the hatted outer sum with
`M = M_nat = M_0(tau(u) zhat)` taken exactly (D14). Assume (F_g). Then `G_Gamma(u) = u^(2g-2) F_Gamma(u)`.
* Window: at genus `g` only `P` through `z^(3g-3)` enters (the single-vertex leaf `z^(3g-2)` needs `P^-1_(3g-3)`), so (F) is used only
  through `z^(3g-3)`, i.e. as (F_g) (FR (I2); SYN2b K9; O20). Genus 3 sees only `c_5`; `c_10 = c_5^2/2` is forced by symplecticity.
* Consequence: the separate input "s.1.5" (S01 s.4) is deleted (K8); Theorem R for `P_can` (flat) is Conjecture N (hatted), graph by
  graph, in every genus.
* Source: S15 s.0-s.4; FR s.0-s.1, s.5; SYN2 s.2 step 6. Status: PROVED given (F) (FR: CONFIRMED). VERIFIED at genus 3 (s.5).
* Remarks (S15 s.4): the dilaton factor `d^(2h-2+n)` and the leaf `dT` of the flat vertex are s.1.1 applied to `q_P` (`z^1` part
  `-(z/d) 1_m`); only symplecticity of `P` beyond (F) is used (not the u-ODE); the proof is genus-independent and graph-local (no total
  potential, no KP/Toda, no TR on any curve). FR: `prod_v d^(2h_v - 2 + n_v) = d^(2g-2)` for every graph (one global unit per graph).
* Corrections: FR C1-C6 = O14-O20 (s.7 item K-list); the leaf formula of S15 s.1.3 is `T' = z y + q_{R'}` (only `d >= 2` used downstream).
* Depends on: `lem:flh-frame`, `lem:flh-dil`, `lem:flh-hat`, `lem:flh-comm`, (I1)-(I4); (U) via UR s.4; Prop. S; arXiv:1303.1043
  Def. 2.2, Prop. 2.4, Def. 2.5, Prop. 2.9 (and [DW] Prop. 2.12, s.7 item L3); the dilaton equation.

**`rmk:flh-fail` (rmk) Where S15 can fail (S15 s.5, FR).** (1) `D != 1`: with (H), `D = c(z)` and the proof gives the flat sum of the gauge
`P_can c(z)^-1`, a different graph decomposition (stored normalisation: genus 3 graphs 1, 2, 5-8 move in regular orders; genus 4 graphs
1, 2, 3, 4, 7, 8 become singular); (2) another split `M'`: same total, different graphs (constant-`M_0` split: 94 of 252 coefficients
differ); truncations of `M_nat` agree only on the orders they control; (3) only outer data `(E_{R'}, q_{R'})` enter; the regrouping
`R~.omega = R'.(M.omega)` (composition of R-actions) is needed only to identify the TOTAL with `F_g` of the quintic; (4) dilaton shift
and leaves; (5) the rows formula equals `E_{R'}` only for symplectic `R'` (true for the actual `R' = R~ M_nat^-1`, with or without (F),
FR s.1); (6) the vertex identification symbol = `M_0`: UR s.2 given (H). Status: remark; the controls are evidence (s.5).

## s.2.5 The assembled chain (SYN2 s.2; unchanged by Stages 4-5: MS4 s.2, SYN4 s.2, SYN5 S8, RS2 s.6.2)

| step | label | statement (final form) | source; status |
|---|---|---|---|
| 1 | `prop:mai-vanish` | `Y^C(v(u))` and `Y^C_Q` lie in `Q[z^-1][[u]]^5`, equal `-e_1` at `u = 0`; `[z^-1](Y^C_Q)_0 = -t_c`; a regular solution is determined by its value at `u = 0` | GI s.2 Lemmas V, Q, 1; PROVED |
| 2 | `thm:mai-G` | uses (U): every regular gauge has `P Y^C(v(u)) = P_11(0,z) Y^C_Q`; `P_can := P P_11(0,z)^-1` satisfies (A) | GI s.3 Theorem G; PROVED given Theorem U |
| 3 | `prop:mai-symp` | uses (U): `P_can(-z)^T eta P_can(z) = eta` | GI s.4 Prop. S (not read); PROVED given Theorem U |
| 4 | `prop:mai-transport` | edges of `P_can` in `K_0 (x) K_0`, as `Q[[u]]`-combinations of frames (M5); leaf point `d q_T in -z + K`, `Pi_0 = d t_c = v + s`, `s in u^2 Q[[u]]`; each graph term a finite sum | RD_2 s.3.1 Lemmas E, T, F, Cor. E' (= model Lemma K); PROVED; uses 2, 3, arXiv:1303.1043 Def. 2.2 |
| 5 | `thm:mai-R` | quintic Theorem R: if ML_K(h), `h <= g-1`, every genus-g stable graph with an edge of the flat sum `P_can.Omega^mod_{v(u)}` lies in `Q[[u]]`; with ML_K(g) the vertex graph is `c_g t_c^(2-2g) + Q[[u]]` | RD_2 s.3.2, MSG_2 s.7; PROVED (refereed twice) |
| 6 | `thm:flh-main` | under (F_g): `G_Gamma = u^(2g-2) F_Gamma` for every stable leg-free genus-g graph | S15; PROVED given (F) |
| 7 | `thm:mlk-main` | ML_K(h) for all `h` | s.2.3; PROVED mod (Q) |
| 8 | `thm:mai-stage2` | Conjecture N(g) (D16), every `g >= 2` | steps 5-7 |

**`thm:mai-stage2` (thm) The Stage-2 theorem.** Assume Theorem U, (F_g) (equivalently (H) and `c^act = 1` through `z^(3g-3)`), and the
quoted theorems (Q): arXiv:1509.06954 Thm 8.1; arXiv:2212.00320 Thm 7.1; arXiv:math-ph/0702045 Thm 5.1 and the linear loop equation of
its s.4 (on the model curve only); arXiv:1303.1043 Def. 2.2, Prop. 2.4, Def. 2.5, Prop. 2.9; optionally arXiv:2502.18115 Cor. 1.2.
Then Conjecture N(g) holds for every `g >= 2`. The s.1.5 identification is not an input (it is step 6).
* Source: SYN2 s.0, s.2 (boxed theorem and input table); MS2 s.1. Status: PROVED modulo (U), (F_g), (Q); AI-written, AI-refereed.
* Later updates (not Stage-2 content, recorded for the main theorem): (QRR) was implicit in "the actual R-matrix" (MS4 s.2, SYN4 s.2);
  (F_g) for every `g` is proved in Stages 4-5 (MS5 s.2: `D_0 = 1` to all orders, modulo (U), (QRR)); (CYC) (PROVED in Stage 5) feeds
  UR s.1, hence Theorem U for the actual quintic and chain steps 1-6 (SYN5 S9; RS1 s.3 table). Final: N(g) for every `g >= 2` modulo
  exactly (U), (Q), (QRR) (MS5 s.2; RS2 s.6.2); (U) is to be PROVED in the note.
* Inconsistency recorded: SYN2 s.0 lists "arXiv:1303.1043 Def. 2.2, Prop. 2.4, Prop. 2.9" without Def. 2.5; SYN2 s.2, MS2 s.1, MS5 s.2
  include Def. 2.5 (used in S15 Lemma 4). Use the full list.

## s.2.6 Universality and Lemma C (plan s.9)

**`lem:uni-pole` (lem) The pole order (BV s.1; BR s.1).** (i) A Theorem-U gauge between a self-dual cyclic system `b` and the
pulled-back model `bt` needs `prod b = prod bt = -v'^5/v` (the `z^0` equation `b_(c-1) q^0_c = bt_(c-1) q^0_(c-1)`, `P_0 = diag(q^0)`
invertible), i.e. `v'^5 = -B v`. (ii) If `v` is a power series with `v(0) = 0` of order `n`, then `ord_0 B = 4n - 5 in {-1, 3, 7, ..}`;
a pole forces `n = 1`, `ord_0 B = -1`, and `v'(0)^4 = -res_0 B`. Status: PROVED; the converse paragraph of BV s.1 (`ord_0 B = -k` for
`k = 2..13`) is SKETCH and not needed (O24). Scope: relative to the fixed model `bmod` (BR s.1).

**`lem:uni-coord` (lem) Lemma C: a change of B is a change of coordinate (BR s.3).** Let `X = (f, g, B)` be self-dual cyclic with
`ord_0 B = -1`, model pull-back `v_X` (`v_X'^5 = -B v_X`, `v_X(0) = 0`, `a := v_X'(0)`, `a^4 = -res_0 B`); quintic `B_q`, `v_q`, `v_q'(0) = 1`;
`phi := v_X^-1 o v_q in w Q[[w]]`, `phi'(0) = 1/a`. Then (a) `X~ := phi^* X` (coefficients `phi' . b^X o phi`) has the self-dual shape and
`prod b~ = B_q`, so it is a univ system with `F~ = phi' f(phi)/t_c'`, `Gm~ = phi' g(phi)/g_q`, `F~(0) = f(0)/a`, `Gm~(0) = -g(0)/a`, and its
model pull-back is `v_q` (its `tau` is the quintic `tau`); (b) `P(phi(w), z)` is a Theorem-U gauge of `X~` with the same `P_11(0,z)`,
`P_00(0,z)`, `P_0(0)`; canonical goes to canonical, stored to stored, with the same `C`, `sigma`, `c0`, `m5`; (c) the pure-`1/z` solutions
pull back; (G) for `X` iff for `X~`; `T''_X~ = T''_X o phi`; (d) `R_X~ = R_X o phi`, `M_X~ = M_X o phi`, `R'_X~ = R'_X o phi`, and every
outer-graph amplitude satisfies `V^X~_Gamma(w) = V^X_Gamma(phi(w))`; (e) a graph is regular for `X` iff for `X~`, and the vertex law
transports. Hatted encoding: `G~(w) = (w/phi(w))^(2g-2) G(phi(w))`. Converse: every univ system and every `phi` (`phi(0) = 0`,
`phi'(0) != 0`) gives a B-perturbed system; the residue of `B` is a coordinate artefact (`phi = c w`).
Status: PROVED (elementary); VERIFIED (BR `rf_reparam`, s.5). Corrections: K7, O21, O22, O25.

**`cor:uni-flat` (cor) Universality, all genera, flat form (SYN2 s.2 Corollary).** For a self-dual cyclic one-pole system with equal
products and its canonical gauge, steps 1-5 and 7 give: every edge graph of the flat sum is regular, and the vertex graph is
`c_g (T/lambda)^(2-2g)` + regular (theory Cor. U; `T = Int_0^u f` (`= Int b_0`), `T'' = c0 T`; `|c0| = m^10/|f(0) g(0)|` for the bvar systems,
the law being even in `c0`). B-perturbed systems are included by Lemma C. Status: PROVED modulo (Q) and Theorem U for the system.
The hatted-engine form follows from step 6 with `D = 1` by construction; matching the univ engine's conventions to Lemma 3: SKETCH.
Corrections: K6/O23 (`c0 = 1/lambda` holds only for `m = 1`, `g(0) = -1`; S5 gives 1024, S6 59049/2); S01 C10 (`lambda = f(0)|g(0)|`).
