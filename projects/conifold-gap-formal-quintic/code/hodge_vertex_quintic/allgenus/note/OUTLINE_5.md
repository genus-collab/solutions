# Outline, part 5: s04 the vanishing tower, the slice K, transport; s05 Theorem R

## s04_transport — target 6 pp

Purpose: fix the two coefficient settings, define the vanishing tower, the slice `K` and its frames, and prove that every edge of an
annihilating symplectic dressing (in particular of `P_can`) lies in `K_0 (x) K_0` and that its leaf point lies in `-z + K`, with every
infinite sum justified u-adically and every graph term finite.

Covers: DC1 s.1 D1-D2, D6-D12, D18-D19, s.2 parts A-B (A1-A6, B1-B7), s.3 A-B, s.4 CV1, CV6, CV13-CV14, s.7 O-2, O-3, O-5..O-8; DC2 D7;
RD1 s.1.1-s.1.4, s.2.1-s.2.2; RD2 s.3.1; RD3 s.4.2; TRF s.2; MSG1 s.1; MSG2 s.6; MRF items 2, 5, 6, 8; ST2 s.2 steps 4, M5; G3 (all pieces)
with G3R, G3R_2 (s.13: K1-K10). Subsections `sec:tra-tower`, `sec:tra-dress`, `sec:tra-lemmas`.

### sec:tra-tower
* `conv:tra-forms` (convention; G3 s.2.1 [G]): the state space `\bbH` with basis `H^a` (REV_2 M13 (e)); weights `wt(v) = 5`, `wt(z^k H^a) =
  4k + nu_a + 1`. Graded setting: `V_w` finite-dimensional, `\cH^gr_+ = \bbH (x) Q[v][z]`, Givental's `\cH = \bbH((z^-1))` with `\cH_+ =
  \bbH[z]`, every sum finite. u-adic setting: `v = v(u)` in `u Q[[u]]`, `\cH^u = \bbH (x) Q((z))[[u]]`, `\cH^u_+ = \bbH (x) Q[[u, z]]`, `\cH^u_- =
  \bbH (x) z^-1 Q[z^-1][[u]]`, the splitting `f = f_(>=0) + f_(<0)`, `\Omega_s(f, g) = Res_z eta(f(-z), g(z)) dz` (`Q[[u]]`-bilinear; `\cH^u_+`,
  `\cH^u_-` Lagrangian); the substitution `v -> v(u)`. The earlier "`\cH = H((z^-1))`, `\cH_+ = H[[z]]`" is inconsistent and is not used
  (G3 s.2.1, "Correction to the outline"). Insertions `psi^k H^b` = vectors `z^k H^b`; lowered components `Y_a = eta(H^a, J)`, `<Y, n> := Sum_a
  Y_a n^a = eta(J, n)`; vector form `P^vec = J P J` (`J`: `a <-> 3 - a`, `eta = 5 J`); `d := (P_0)_33 = (P_0^vec)_00`; (A) and `Y'` lowered,
  Lemmas E and T in vector form (O-7). Summable families and maps of finite type: G3 s.2.2 [G].
* `lem:tra-sum` (lem; new; G3 s.2.2 [G]): the summation lemma (a map of finite type commutes with summable sums), with its examples.
* `lem:tra-theta` (lem; new; G3_2 s.2 Step 1 [G]; G3R_2 s.13 K8 (b)): the operator `Theta = 1 + N` through which the shift rule acts on the
  sequences `(Pi_i(q))_(i >= 0)`, and its inverse `Sum_r (-N)^r`, in both settings; stated once in this subsection.
* `def:tra-tower` (defn; DEC D1): `Pi_i(q) := -[z^i] <Y^C(-z), q(z)>` (`i >= 0`), `Pi_0` the vanishing covector; shift rule `Pi_i(z^k H^b) =
  Pi_0(z^(k-i) H^b)`, zero for `k < i`; `Pi_0(z^k H^b) = (-1)^(k+1) [z^-k] Y^C_b`, zero unless `b + k = 1 mod 5`; first values (`Pi_0 = 1` on
  `H`, `-v` on `z 1`, `-v^2/2` on `z^2 H^4`, ..) from `def:gau-van` (DC1 D7; ST2 M3; S01 C4); in the graded setting `Pi_i` maps `V_w` into
  `Q v^((w - 4i)/5)` (G3_2 s.2 Steps 0, 2). No `\draftnote` (the definition through `S(z)` is not used).
* `def:tra-slice` (defn; DC1 D8): `K := {q : Pi_i(q) = 0, i >= 1}`, `K_0 := K cap ker Pi_0 = {n : <Y^C(-z), n(z)> in z^-1 Q[z^-1][[u]]}` in
  the u-adic setting (G3 s.3 item 6), `K = K_0 (+) Q H`; the affine slice `-z + K = -z 1 + (K_0 (+) Q H)`, `Pi_0(-z 1 + n + mu H) = v + mu`
  (`mu` a declared local scalar; REV R8 (a)).
* `lem:tra-frames` (lem; DC1 A1, theory's "Lemma 0"): statement G3_2 s.1 (a)-(d) [G] with K1: (a) `\cH^A_+ = K_0^A (+) A{z} H`, graded:
  `V_w = (K_0 cap V_w) (+) (A{z} H cap V_w)`; (b) for `(k, a)`, `a != 1`, the unique `frame(k,a) = z^k H^a + Sum_(j<k) f_j z^j H` in `K_0`,
  `f_j = c_(k,a,j) v^((4(k-j) + nu_a + 1)/5)`, `c_(k,a,j)` in `Q`, zero unless the exponent is an integer, which is then `>= 1`; no frame with
  `a = 1`; (c) frame coordinates are bijective (graded: finitely supported families in `Q[v]`, on `K_0 cap V_w` the `(k,a)`-member a rational
  multiple of `v^((w - 4k - nu_a - 1)/5)`), the inverse by summable families; (d) the u-adic frames are the graded ones at `v = v(u)`.
  Trivial frames: `frame(k,a) = z^k H^a` iff `k < nu_a + 1`, ten of them (G3 s.3 item 7). Examples (RD1 s.1.3): `frame(1,0) = e = z 1 + v H`,
  `frame(2,4) = z^2 H^4 + (v^2/2) H`, `frame(3,3) = z^3 H^3 + (v^3/6) H`, `frame(4,2) = z^4 H^2 - (v^4/24) H`, `frame(6,0) = z^6 1 + v z^5 H -
  (v^5/120) H`. Proof G3_2 s.2 Steps 0-5 [G] (Step 2 replaces the source's "top power" argument, which holds only for polynomial
  combinations; Step 4 is the former [S] step). Local letters for G3's `s(z)`, `sigma(q)`, `pr`. Delete the [S] `\draftnote`.
* `cor:tra-slice` (cor; A2): the description of `-z + K`, from `Pi_0(-z 1) = v`, `Pi_i(-z 1) = 0`, `Pi_0(H) = 1`, `Pi_i(H) = 0` (`i >= 1`). G3_2
  s.2 Corollary (ii) [G]; RD1 s.1.3 [W].
* `lem:tra-coeff` (lem; A3 = ST2 M5): over `Q[[u]]`, `K_0` is spanned u-adically by the frames, every coefficient has `v`-exponent `>= 1`;
  every insertion of Theorem R is a `Q[[u]]`-combination of homogeneous elements of `K_0` and of `H`. Proof G3_2 s.2 Step 5 [G] (the count
  `wt(k,a) - wt(j,1) = (4(k-j) + nu_a + 1)/5 >= 1`; TRF s.2; ST2 M5). The basis statement for `K_0 cap V_w` is `cor:mod-basis` (s06).
* `lem:tra-conv` (lem; DC1 O-3): statement G3_3 s.1 (a)-(f) [G]: (a) frame families are summable, with the u-order bound, and give the
  decomposition of `lem:tra-frames` (c) (G3_3 s.2); (b) correlators extend to symmetric `Q[[u]]`-multilinear maps of finite type on
  `(\cH^u_+)^n`, frame expansions inside a correlator are summable and agree with flat expansions, and no flat component `z^j H` (`j >= 1`)
  of a frame occurs as an insertion of its own (G3_3 s.2 (b1)-(b4); K7: `e_0` is the least `v`-exponent of the finitely many flat
  correlators with `Sum k_i <= D`); (c) `[z^-n] Y^C(v(u))` by the closed form of `lem:gau-vanshape` (K2), of u-order `>= 4n/5`, so every
  pairing and product of s04 is a ring operation in `Q((z))[[u]]`, resp. `Q((x))[[u, psi']]`; (d) a flat insertion `z^k H^b` with `k >= 2` has
  `wt >= 3/5`, and slot degrees at a vertex of genus `g_v` with `m` such insertions are `<= 3g_v - 3 + n - m` (G3_3 s.3); (e) the vertex
  resummation `Sum_j delta_P^j (F^m_g)^((j))(v)/j! = c_g (v + delta_P)^(2-2g)` and the leaf sum of `thm:red-main` steps 8-9 (G3_3 s.4; its
  Remark gives the finite route); (f) the Taylor expansion of `Fbar_(g_v)` at `-z 1` in the direction `q_T + z 1` converges u-adically iff
  `d(0) = 1`, and then equals the homogeneous form of `def:red-flatsum` (G3_3 s.5 with the polynomial-in-`m` argument of G3R_2 s.8, K3).
  Print the dependencies with the lemma (K10): `lem:set-modhom` ((b), (e)), `lem:set-dilaton` ((f)), `lem:gau-diag` ((f), for `P_can`),
  ML_K(g) ((e), as the hypothesis of `thm:red-main` (b)), `lem:tra-T` (d) and `v'(0) != 0` ((e)). Delete the "argued" `\draftnote`.
* `lem:tra-twoslot` (lem; new; G3_2 s.3 [G]): a bivector annihilated by every `Pi_i` in each slot is the summable double frame series of its
  flat components, and conversely; used in `cor:tra-edges` step 5.
* `rmk:tra-consistency` (rmk; A5 (a)-(c) only): the string operator preserves `-z + K`; off `K` no function of `Pi_0` can be the singular
  part (`<psi H, 1>_h = <H>_h`): the necessity of `K`; the dilaton and Euler-dilaton consistency, values computed (appendix C). [W] RD1
  s.2.1, RD3 s.4.2, MSG2 s.9(d), TRF s.2. Dropped (D8): `rmk:tra-isotropic`, `rmk:tra-descendant`, items (d)-(e).

### sec:tra-dress — dressings, edges, the leaf point
* `def:tra-dressing` (defn; DC1 D11): a dressing of the model at `v(u)` (`v(0) = 0`, `v'(0) != 0`): `P = Sum z^k P_k(u)`, `P_k` in `gl_5(Q[[u]])`,
  graded, `P_0` invertible, acting on lowered components; (Sp) `P(-z)^T eta P(z) = eta`; (A) `Y' := P Y^C(v(u))` has no positive power of `z`;
  `T_P := -[z^-1] Y'_0`; `d := (P_0)_33`; `R := P_0^-1 P`. Examples: the canonical gauge of any self-dual one-pole system with the model's
  product; `D P_A` with `D` flat-regular, graded, symplectic, `(D Y^C_A)_+ = 0`.
* `def:tra-edge` (defn; DC1 D12): `E_R(psi, psi') := (eta^-1 - R^-1(psi) eta^-1 R^-1(psi')^T)/(psi + psi')`, `E_P := E_R`; over `Q[[u]]` the
  numerator is divisible by `psi + psi'` (one line, G3_4 s.2 Step 0 [G]; G3 s.3 item 4), so `E_R` lies in `\bbH (x) \bbH (x) Q[[u, psi, psi']]`;
  leaf point `q_T := P^-1(-z 1)` (vector form), `qt := d q_T`; `Omega^P := P.Omega^mod_(v(u))`: [PPZ15] Def. 2.2 for `R`, preceded by the constant
  frame change `P_0`, with the translation `\ttr_R(z) = z(1 - R^-1(z) 1)` of [PPZ15] Def. 2.13 (= `T_b` of Prop. 2.12).

### sec:tra-lemmas — Lemmas E, T, F and the canonical gauge
* `lem:tra-E` (lem; B1, "Lemma E"): statement and proof G3_4 s.2 [G]: `R = 1 + O(z)` in `gl_5(Q[[u, z]])` symplectic, `f` in `\cH^u_-` (vector
  form); the contraction of `\Omega_s(f, .)` into the first slot of `E_R` converges u-adically, and `(\Omega_s(f, .) (x) 1) E_R(psi') =
  -R(psi')^-1 (R f)_(>=0)(psi')` in `\bbH (x) Q[[u, psi']]`; so `E_R` annihilates `\Omega_s(f, .)` (in either slot) iff `R f` lies in `\cH^u_-`. The
  expansion ring is `Q((x))[[u, psi']]` (G3_4 s.1: the source's `Q((x^-1))[[psi']]` is right only in the hatted setting); the side remark
  with `Q((psi'))((x))[[u]]` (K5); a local letter for `(x + psi')^-1` (not `iota`). Delete the O-8 `\draftnote`.
* `cor:tra-edges` (cor; B2, "Corollary E'"): if `P = P_0 R` satisfies (Sp), (A), the edges annihilate every `Pi_i` in each slot; `E_P` lies in
  `K_0 (x) K_0` (completed), as `Q[[u]]`-combinations of frames; and (A) iff the edges annihilate `Pi_0` iff they annihilate the whole
  tower. Proof B2 steps 1-5: 1-4 [W] (RD2 s.3.1; TRF s.2), step 3 with `P^vec J^C = eta^-1 (P Y^C)` (G3R_2 s.13 K6); step 5 by
  `lem:tra-twoslot` (G3_2 s.3), with the duality step of MSG1 s.1 (triangularity `Pi_i(z^k H) = Pi_0(z^(k-i) H)`, `Pi_0(H) = 1`).
* `lem:tra-T` (lem; B3, "Lemma T"): for `P` with (Sp), (A): (a) `\Omega_s(f, q_T) = \Omega_s(P f, -z 1)`; (b) `Pi_0(q_T) = T_P`, `Pi_i(q_T) = 0`
  (`i >= 1`), so `q_T` lies in `K`; (c) `q_T = -(z/d) 1 + T_(>=2)`, `T_(>=2)` in `z^2 \cH^u_+`, so `qt = d q_T` lies in `-z + K` with `Pi_0(qt) =
  d T_P =: v + delta_P` (REV R8 (a)); (d) `delta_P` lies in `u^2 Q[[u]]` (`(P_0)_00 (P_0)_33 = 1` by (Sp) at `z^0`). Proof B3 steps 1-6 [W]/[W-c]
  (RD2 s.3.1; TRF s.2; MSG2 s.6); step 6: write the exclusion of the other `u^1` terms from the grading. Step 7 (`delta_P = u^2/30 + ..`,
  frame independence) [V]: appendix C only. Never write "the leaf point has `Pi_0 = t_c`" without (b)-(c).
* `lem:tra-F` (lem; B4, "Lemma F"): every `Cont_Gamma` is a finite sum of products of model correlators at `v(u)` (monomials `C v(u)^e`) and
  entries of `P_k`, `(P^-1)_k`, `d^(+-1)`. Proof B4 steps 1-3 [W]/[W-c] (RD2 s.3.1; TRF s.2): step 1 with the homogeneous vertex tensor of
  `def:red-flatsum` (DEC Q5); where `lem:set-dilaton` is applied, `d(0) = 1` is stated, and then the vertex tensor is the Taylor expansion at
  the base point (`lem:tra-conv` (f)); step 1 at `g_v = 1` [W-c] (the log term is killed by `n >= 1` derivatives). Correction O-6:
  "monomials in `v(u)`", not "Laurent monomials in `u`".
* `rmk:tra-K`: dropped (D8). Omitted (B6, theory's Fock-space form): its quoted formula has no statement number (rule 4); no label.
* `prop:tra-can` (prop; B7): `P_can` is a dressing with (Sp) (`prop:gau-symp`) and (A) (`thm:gau-G`), `Y' = Y^C_Q`, `T_P = t_c` (`lem:gau-vanq`
  (d)), `d(0) = 1` (`lem:gau-diag`); hence its edges lie in `K_0 (x) K_0` (frames over `Q[[u]]`), `q_T` lies in `K` with `Pi_0(q_T) = t_c`, `d q_T`
  lies in `-z + K` with `Pi_0 = v + delta_P`, `delta_P` in `u^2 Q[[u]]`, and every graph term is a finite sum. Proof [W] ST2 s.2 steps 1-4;
  RD1 s.1.4 example (i).

Length 6 pp.

## s05_reduction — target 3 pp

Purpose: reduce regularity of every flat graph term, and the vertex law, to the Model Lemma on `K` (Theorem R), using only the
transport of s04; no gauge equation, no quintic-specific datum beyond `T_P`.

Covers: DC1 s.1 D4, D5, D12, D14, s.2 part C (C1-C6), s.3 C, s.7 O-2, O-3 (e), (g), O-4, O-9; RD2 s.3.2-s.3.4; TRF s.2; MSG2 s.7;
MRF item 8; ST2 s.2 step 5; G3 s.3 item 2; G3_3 s.4-s.5. Subsections `sec:red-sum`, `sec:red-thm`.

* `def:red-flatsum` (defn; DC1 D12, RD1 s.1.5; DEC Q5): for a dressing `P` and a leg-free stable graph `Gamma` of genus `g`, `Cont_Gamma`
  (including `|Aut Gamma|^-1`, `def:set-graphsum`) contracts `E_P` along the edges with the vertex tensors defined by the HOMOGENEOUS form
  `d^(2g_v - 2 + n_v) \dd^(n_v) Fbar_(g_v)(d q_T)` (the form in which the flat sum is computed, S15 s.1.1; G3 s.3 item 2); it equals the Taylor
  expansion `\dd^(n_v) Fbar_(g_v)(q_T)` at the base point iff `d(0) = 1` (`lem:tra-conv` (f)), which holds for `P_can` (`lem:gau-diag`; GS
  s.2.3); `F^P_g = Sum_Gamma Cont_Gamma`; vertex graph, edge graph; regular = `Cont_Gamma` in `Q[[u]]`; for `P = P_can` the flat term
  `F_Gamma` (s07). `Fbar_h` and its homogeneity: `lem:set-dilaton`. The `n`-th derivative is written `\dd^n` (roman), the dilaton factor `d`.
* `def:red-ML` (defn; DC1 D14, final form): ML_K(h), `h >= 0`: (i) for `n >= 1`, `2h - 2 + n > 0`, `X_1..X_n` in `K` homogeneous over `Q[v]`:
  `<X_1 .. X_n>_h(v) - (F^m_h)^((n))(v) prod_i Pi_0(X_i)` lies in `Q[v]`; (ii) for `h >= 2`, `F^m_h = c_h v^(2-2h)`. Equivalent function form
  `Fbar_h | (-z + K) = F^m_h o Pi_0 + Reg_h` (RD1 s.2; MSG1 s.5) as one sentence.
* `lem:red-genus` (lem; O-3 (e)): in a stable graph with an edge, every vertex of a non-zero term has `n >= 1` half-edges and genus
  `g_v <= g - 1` (a genus-`g` vertex forces a tree of genus-0 vertices; a genus-0 vertex with `<= 2` half-edges vanishes by dimension since
  the leaves are `O(z^2)`). Proof [W] RD2 s.3.2 step (2), first sentence; TRF s.2; MRF item 8 (a), (d).
* `thm:red-main` (thm; C1, "Theorem R"): fix `g >= 2` and a dressing `P` with (Sp), (A). (a) If ML_K(h) holds for `0 <= h <= g-1`, then
  `Cont_Gamma` lies in `Q[[u]]` for every genus-`g` graph with an edge. (b) If moreover ML_K(g) holds, `Cont_vertex = c_g T_P^(2-2g) + Q[[u]]`.
  The gauge equation is not used. Proof DC1 s.3 C1 steps 1-10: 1 [W-c] (frames for every half-edge, `lem:tra-coeff`); 2 `lem:red-genus`;
  3-5 [W] (RD2 s.3.2 step (2); MSG2 s.7; MRF item 8 (b)); 6 `lem:tra-conv` (a)-(d) [G]; 7 [W] (RD2 s.3.2 steps (2)-(3)); 8-9 `lem:tra-conv`
  (e) [G] (G3_3 s.4; or its finite route); 10 [W]. Hypothesis used in (b): `v(0) = 0`, `v'(0) != 0` (O-2). Corrections: name; ML for
  `h <= g` suffices (not MSG2's "every `h >= 0`"); no claim that each hypothesis is needed.
* `cor:red-can` (cor): for `P = P_can`: given ML_K(h <= g-1) every edge graph of `P_can.Omega^mod_(v(u))` is regular; given ML_K(g) the
  vertex graph is `c_g t_c^(2-2g) + Q[[u]]` (`prop:tra-can`, `thm:red-main`). [W] ST2 s.2 step 5.
* `rmk:red-hyp` (rmk; C2 (i)-(v)): what of ML is used. [W] RD2 s.3.2 "Hypotheses needed"; TRF s.2.
* `rmk:red-leaf` (rmk; C3, "Corollary L"; kept under D8: mechanism): with the leafless normalisation `c_g vhat^(2-2g)`, `vhat = v/d`, the
  leaves turn the vertex into `c_g t_c^(2-2g) + ..`, first correction `-(g-1) c_g (2/vhat^(2g-1))(t_c - vhat)`. The input `t_c - vhat = u^2/30
  + O(u^3)` is computed only: `\draftnote{no written derivation of the u^2-coefficient of delta_P; the frame independence 1/15 = 2(7/10) -
  4/3 uses the leaf coefficient -p_{1,4}(0), whose derivation is not written either}` (outside the chain); numbers in appendix C.
* Omitted (D8): Theorem R(c), Prop. I, the deformation route (C4-C5); Cor. N, N_9 are superseded (C6).

Length 3 pp.
