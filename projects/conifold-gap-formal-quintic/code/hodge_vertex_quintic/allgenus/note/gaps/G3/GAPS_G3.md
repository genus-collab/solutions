# Gap cluster G3 (s04, transport): lem:tra-frames step 4, lem:tra-conv, lem:tra-E step 2 (part 1 of 4)

Author-facing working file of the gap-closing phase; not part of the note. Pieces: this file (s.0 results, s.1 what was searched,
s.2 the standing conventions every proof below uses, s.3 observations outside the cluster); `GAPS_G3_2.md` (item 1: the frame
lemma, all steps, both coefficient settings; the two-slot frame expansion of edges); `GAPS_G3_3.md` (item 2: lem:tra-conv, parts
(a)-(d), and two added parts (e), (f)); `GAPS_G3_4.md` (item 3: lem:tra-E step 2, the ring of O-8, Lemma E in that ring; the
computations with failing controls; files). Labels of the statuses as in the binding rules: PROVED (argument written here, in the
note's notation, dependencies named), VERIFIED (exact computation; range and field stated), FOUND, CITED, OPEN. Every
computation below is evidence only; no proof uses it.

## 0. Results (one resolution per item)

| item | statement (outline: OUTLINE_5 s04) | resolution | where |
|---|---|---|---|
| G3-1 | `lem:tra-frames` step 4: `HH_+ = K_0 (+) Q[[z]] H` for general (infinite) elements, frame coordinates = flat components with `a != 1`; the ring stated | PROVED, over the graded ring `Q[v]` and over `Q[[u]]` | part 2 s.2 (Steps 1-5) |
| G3-1x | (found while proving G3-1) `lem:tra-frames` step 3 ("top power") is valid only for polynomial combinations; for power series in `z` it is replaced by the u-adic invertibility of one operator `Theta = 1 + N` | PROVED (replacement) | part 2 s.2, Step 1-2 |
| G3-2a | `lem:tra-conv` (a): the decomposition converges u-adically | PROVED | part 2 s.2 Step 4; part 3 s.1 |
| G3-2b | `lem:tra-conv` (b): frame expansions inside a model correlator converge u-adically, agree with the flat expansion, never isolate `psi^k H` | PROVED (no homogeneity needed for the convergence) | part 3 s.2 |
| G3-2c | `lem:tra-conv` (c): `[z^-n] Y^C = O(v^(n - floor(n/5)))`, so the pairings of s04 converge | PROVED (the decay itself is FOUND, exactly, in GI s.2 Lemma V = `lem:gau-vanshape`; the convergence of the pairings is the ring statement of s.2 below) | part 3 s.1 |
| G3-2d | `lem:tra-conv` (d): each leaf of psi-degree `>= 2` raises the exponent by `>= 3/5`; slot degrees `<= 3h - 3 + n - m` | PROVED (elementary) | part 3 s.3 |
| G3-2e | (added; "argued" in RD2 s.3.2 step (3), TRF s.2) the vertex resummation `Sum_j delta_P^j F^m_g^(j)(v)/j! = c_g (v + delta_P)^(2-2g)` and the extension/re-ordering of the leaf sum used by `thm:red-main` steps 8-9 | PROVED | part 3 s.4 |
| G3-2f | (added; flat referee "unstated: d(0) = 1") the Taylor expansion of `Fbar_h` at `-z 1` in the direction `q_T + z 1` converges u-adically iff `d(0) = 1`, and then equals the homogeneous form of `lem:tra-F` step 1 | PROVED; recommendation for `def:red-flatsum` | part 3 s.5 |
| G3-2g | (added; DC1 B2 step 5 "the duality step is implicit: write it") an edge bivector annihilated by every `Pi_i` in each slot is the summable frame double series of its flat components | PROVED | part 2 s.3 |
| G3-3 | `lem:tra-E` step 2 (O-8): the ring of the expansion of `1/(x + psi')` and of the residue pairing | PROVED: the ring is `Q((x))[[u, psi']]`, with `f` in `HH^u_- = H (x) z^-1 Q[z^-1][[u]]`; the source's `Q((x^-1))[[psi']]` is wrong in the flat frame | part 4 s.1-s.2 |

No item of the cluster is left OPEN. Inputs used by the proofs, by label: `lem:gau-vanmod`, `lem:gau-vanshape` (the closed form of
`Y^C`), `def:tra-tower`, `def:tra-slice`, `def:tra-dressing`, `def:tra-edge`, `def:set-modcohft` (correlators of the model at `v(u)`;
the CohFT degree bound), `lem:set-dilaton` (only in G3-2f), `def:red-ML` (only in G3-2e, as the hypothesis ML_K(g) of
`thm:red-main` (b)), `lem:tra-T` (d) (only in G3-2e: `delta_P` in `u^2 Q[[u]]`), and (Sp) for `lem:tra-E`. No quoted theorem
is used. Two of these inputs carry their own `\draftnote`s in the outline (homogeneity in `def:set-modcohft`; the citation of the
dilaton equation in `lem:set-dilaton`). From `def:set-modcohft` the proofs use: in G3-2b only that each flat model correlator is a
Laurent polynomial in `v` (so that it can be evaluated at `v(u)`), and the CohFT degree bound; in G3-2e the exponent formula
`e = 2 - 2h + Sum wt(k_i, a_i)` (homogeneity). The dilaton equation enters only G3-2f.

**Attribution of O-8.** The task text locates "O-8" in `stage2/STAGE2_SYNTHESIS_2.md` s.5. Item O8 there concerns the logleg
note ((E1) listed as evidence for (Q1)); it has nothing to do with Lemma E. The O-8 meant by the outline's `\draftnote` at
`lem:tra-E` step 2 is DIGEST_C1_6 s.7 O-8 ("The ring of Lemma E's expansion", citing TRF s.2's phrase "legitimate in
`Q((x^-1))[[psi1]]`"). G3-3 resolves that item.

## 1. What was searched (FOUND / not found)

Read in full: `theory/REDUCTION_DRAFT.md`, `_2.md`, `_3.md` (RD1-RD3); `theory/referee/REFEREE.md` (TRF); `gauge/GAUGE_IDENTITY.md`
(GI), `gauge/GAUGE_IDENTITY_HAT.md` (GH); `stage2/STAGE2_SYNTHESIS.md`, `_2.md`. Read in part: `model/MODEL_SINGULAR_2.md` s.6-s.7,
`model/MODEL_SINGULAR.md` s.1 (MS1, MS2); `model/referee/REFEREE.md` items 5-10 (MRF); `stage2/flat/S15_PROOF.md` s.0-s.5 and its
referee's table s.1 and C5; `stage2/dict/DICT_PROOF.md` (headings, Cor. 4-5, s.6); the docstrings of `theory/th_core.py`; text search of
`stage2/logleg/*.md` (nothing on s04 there); the design files OUTLINE_5, NOTATION, NOTATION_2, LABELS (s04),
DIGEST_C1, _2, _4, _6, DIGEST_B_3 S22-S24, DECISIONS, OUTLINE_REVIEW R7-R8, OUTLINE_REVIEW_2 M13. Text searches (all of `HVQ`,
`*.md`): "u-adic", "converg", "Q((z))[[u]]", "L[[u]]", "summable", "completed tensor", "top power", "direct sum".

* **FOUND.** (i) The exact form of the decay in G3-2c: GI s.2 Lemma V ("`Y^C_b` is a sum of monomials `v^j z^-n` with `n = 5i + c_b + 1`,
  `5j = 4n + c_b + 1`: at order `u^N` ... only `n <= (5N - c_b - 1)/4`, so `Y^C(v(u), z)` is in `Q[z^-1][[u]]^5`"), digested as S23 =
  `lem:gau-vanshape`. Since `0 <= c_b + 1 <= 4`, `i = floor(n/5)` and `j = n - i`: the "`O(v^(n - floor(n/5)))`" of RD2 s.3.1 is an
  exact monomial. (ii) The ring `L[[u]]`, `L = Q((z))`, is already the ring of the gauge identity (GI s.3, `thm:gau-G`: "in
  `Q((z))[[u]]^5`"); s.2 below makes it the ring of all of s04, which keeps s03 and s04 consistent. (iii) The hatted analogue of
  Lemma E (GH s.3, Lemma L) states its finiteness ("at each u-order both sides are finite sums", via the polynomiality of the
  hatted rows); this is the setting where TRF's `Q((x^-1))[[psi1]]` is right (part 4 s.3), not the flat one.
* **Not found anywhere** (hence proved here): the decomposition for infinite elements (G3-1); a convergence statement with a
  topology or a ring (all sources: "converges u-adically", "argued, not formalised" (MS2 s.10)); the precise ring of O-8; the
  role of `d(0) = 1` (only the flat referee's one sentence).

## 2. Standing conventions (proposed content of `conv:tra-forms`; every proof of parts 2-4 uses them)

**2.1 Two coefficient settings.** Weights: `wt(v) = 5`, `wt(z^k H^a) = 4k + nu_a + 1` (the weight of `def:mod-weights`: `V_w` is the
span of `v^m z^k H^b` with `w_b + 4k + 5m = w`, `w_b = nu_b + 1`); for scalar series `wt(v^m z^k) = 5m + 4k`, `k` in `Z`.
* *Graded setting.* `V_w` is finite-dimensional, and its elements are polynomials in `z` (`4k <= w`). Put `HH^gr_+ := (+)_w V_w =
  H (x) Q[v][z]`. The vanishing solution is homogeneous: `Y^C_b` has scalar weight `nu_b + 1` (`lem:gau-vanshape`: `[z^-n] Y^C_b`
  is a rational multiple of `v^(n - floor(n/5))`, and `5(n - floor(n/5)) - 4n = n - 5 floor(n/5) = nu_b + 1`). Hence for `q` in `V_w`
  the pairing `<Y^C(-z), q(z)>` is a homogeneous element of `Q[v]((z^-1))` of scalar weight `w` (finitely many positive powers),
  and `Pi_i(q)` is a rational multiple of `v^((w - 4i)/5)`, zero unless that exponent is a non-negative integer. In this setting
  Givental's convention `HH = H((z^-1))` with `HH_+ = H[z]` is consistent, and every sum is finite.
* *u-adic setting.* `v = v(u)` in `u Q[[u]]`, `v != 0` (s04 uses only `ord_u v >= 1`; `v'(0) != 0` enters only in G3-2e and in
  `lem:tra-T` (d)). The loop space is
  `HH^u := H (x) Q((z))[[u]]`, `HH^u_+ := H (x) Q[[z]][[u]] = H (x) Q[[u, z]]`, `HH^u_- := H (x) z^-1 Q[z^-1][[u]]`.
  `Q((z))[[u]]` is the power-series ring in `u` over the field `Q((z))`. Equivalently its elements are the two-sided series
  `Sum_(k in Z) c_k z^k`, `c_k` in `Q[[u]]`, with `c_k -> 0` u-adically as `k -> -oo`; the product is the Cauchy product, a
  u-adically convergent sum at each `z^k` (both tails are small), and it is finite at each order of `u`. Every element splits uniquely,
  coefficientwise in `u`, as `f = f_(>=0) + f_(<0)` with `f_(>=0)` in `HH^u_+`, `f_(<0)` in `HH^u_-`. `Omega_s(f, g) = Res_z eta(f(-z),
  g(z)) dz` is a `Q[[u]]`-bilinear antisymmetric form on `HH^u` (`z -> -z` is an automorphism of `Q((z))`, `Res_z` is taken
  coefficientwise in `u`), and `HH^u_+`, `HH^u_-` are Lagrangian.
* *Substitution.* `v -> v(u)` sends each `V_w` into `HH^u_+`, and each homogeneous element of the negative part, `Sum_n c_n v^((w +
  4n)/5) z^-n`, into `HH^u_-` (the `z^-n` coefficient has u-order `(w + 4n)/5 -> oo`). In particular `Y^C(v(u), z)` lies in
  `Q[z^-1][[u]]^5` and `J^C = eta^-1 Y^C`, `z^-k J^C` (`k >= 1`) lie in `HH^u`, resp. `HH^u_-`.
* *Correction to the outline.* `conv:tra-forms` as drafted ("`HH = H((z^-1))`, `HH_+ = H[[z]]`") is inconsistent: `H[[z]]` is not a
  subspace of `H((z^-1))` (their intersection is `H[z]`), and `z^-1 H[[z^-1]]` cannot be paired with `H[[z]]` without a convergence
  condition. The edges and leaf points of s04 are power series in `z` with infinitely many terms, and `Y^C` has infinitely many
  negative powers, so neither Givental's `H((z^-1))` nor `H((z))` contains all objects of s04. `HH^u` above does (the gauge `P`
  in `gl_5(Q[[u, z]])`, `Y^C(v(u))`, `P Y^C`, `q_T = P^-1(-z 1)`, `E_P`). Write `conv:tra-forms` with the two settings above.

**2.2 Summable families; maps of finite type.** A family `(x_l)` in `HH^u_+` is *summable* if for every flat index `(k, b)` and every
`N` only finitely many `x_l` have `[z^k] x_l^b` outside `u^N Q[[u]]`; the sum is taken coefficientwise (finite modulo each `u^N`).
The same definition is used in `H (x) H (x) Q[[psi, psi']][[u]]` (edges; indices `(k, l, a, b)`). A family `(y_l)` in `Q((u))` is
summable if all `y_l` lie in one `u^(N_0) Q[[u]]` and for every `N` only finitely many lie outside `u^N Q[[u]]`. Any regrouping of
a summable family (a partition into subfamilies) consists of summable families, and the sum of their sums is the total sum
(coefficientwise this is a statement about finite sums). A summable family of elements of `Q[[u]]` has its sum in `Q[[u]]`.
A `Q[[u]]`-linear map `phi : HH^u_+ -> Q((u))` is *of finite type* if `phi(HH^u_+)` lies in some `u^(N_0) Q[[u]]` and for every `N`
there are a finite set `F_N` of flat indices and an integer `M_N` such that `phi(x)` lies in `u^N Q[[u]]` whenever `[z^k] x^b` lies
in `u^(M_N) Q[[u]]` for all `(k, b)` in `F_N` (the other coefficients of `x` being arbitrary).
**Summation lemma.** If `phi` is of finite type and `(x_l)` is summable, then `(phi(x_l))` is summable and `phi(Sum x_l) = Sum phi(x_l)`.
*Proof.* Fix `N`, take `F_N`, `M_N`. Only a finite set `L_0` of indices `l` has a coefficient `[z^k] x_l^b`, `(k, b)` in `F_N`, outside
`u^(M_N) Q[[u]]`; for `l` not in `L_0`, `phi(x_l)` lies in `u^N Q[[u]]`. The element `Sum_(l not in L_0) x_l` has all its
`F_N`-coefficients in `u^(M_N) Q[[u]]` (each is a coefficientwise sum of such), so `phi(Sum_l x_l) - Sum_(l in L_0) phi(x_l)` lies in
`u^N Q[[u]]`. As `N` is arbitrary, the claim follows. QED
*Examples* (used below): `x -> [z^k] x^b`; `Pi_i` (part 3 s.1); a model correlator in one slot, the others fixed (part 3 s.2);
`(Pi_i (x) 1)` and the contraction of `lem:tra-E` in the first slot of a bivector (coefficientwise in the second slot).
*Remark.* On `HH^u_+ = H (x) Q[[u, z]]` one may use instead the `(u, z)`-adic topology; the frame families below are summable for
both (part 2 s.2 Step 4). The coefficientwise notion is weaker, and is all that the proofs need.

## 3. Observations outside the cluster (for the editor; nothing here is needed by G3)

1. **`cor:mod-basis` needs no input from s06.** By G3-1 in the graded setting (part 2, Step 5), `n -> ([z^k] n^a)_(a != 1)` is an
   isomorphism from `K_0 cap V_w` onto the span of the triples `v^m z^k H^a`, `a != 1`, of weight `w`. So the `v^m frame(k,a)` of weight
   `w` are a basis of `K_0 cap V_w` and `dim K_0 cap V_w = d(w)` (the number of triples with `b != 1`, `thm:mod-lattice` (a)), with no
   use of `thm:mod-lattice` (c) or `prop:mod-dii` (the source routes, DR O2, needed D(ii)). Likewise `V_w = (K_0 cap V_w) (+)
   span{v^m z^k H}` (DICT Cor. 4, proved there through the log map and D(ii)) is G3-1 (a) in weight `w`, and the sentence of
   `stage2/dict/DICT_PROOF.md` s.6 ("Gradings": "the model's `K_0` is weight-homogeneous, being cut out by the homogeneous `Pi_i`;
   that is standard in Stage 1 and assumed here") is proved by part 2, Step 0 (0.3) and Step 2 (`Pi_i` maps `V_w` to `Q v^((w-4i)/5)`).
2. **`lem:tra-F` step 1 and `def:red-flatsum`.** For a dressing with `d(0) != 1` the vertex tensor `d^n Fbar_h(q_T)` is not defined by
   the Taylor expansion at the base point (G3-2f). Define the vertex tensor by the homogeneous form `d^(2h-2+n) d^n Fbar_h(d q_T)`
   (the form in which the flat sum is computed, S15 s.1.1), and state that it equals the Taylor expansion when `d(0) = 1`, which
   holds for `P_can` (`lem:gau-diag`: `p_(0,0)(0) = 1`, `d = 1/p_(0,0)`).
3. **`thm:red-main` step 9** can be written without extending the leaf sum: part 3 s.4, Remark.
4. **`def:tra-edge`** ("a power series by (Sp)"): over the coefficient ring `Q[[u]]` the divisibility of the numerator by `psi + psi'` is
   one line (part 4 s.2 step 0); it is worth printing, since `E_R` is used with `Q[[u]]` coefficients throughout.
5. **M13 (e) of OUTLINE_REVIEW_2** (the state space `H = Q^5` against the class `H`): the proofs below keep `H` for both, as NOTATION
   does; they read unambiguously (`H^a` is a basis vector, `z^k H` a flat insertion).
6. **The sources' `K_0 = {n : <Y^C(-z), n(z)> in z^-1 Q[[z^-1]]}`** must read `z^-1 Q[z^-1][[u]]` in the u-adic setting (s.2.1).
7. **Trivial frames (erratum to DC1 D9 / RD1 s.1.3).** `frame(k,a) = z^k H^a` iff `k < nu_a + 1`: ten frames, the eight listed in the
   sources plus `(2,2)` and `(3,2)` (`frame(2,2) = z^2 H^2`, `frame(3,2) = z^3 H^2`). Proof in part 2, Step 3; computed for `k < 40`.
   This agrees with MS1 s.1's membership list for `K` (`psi^k H^b` in `K` iff `k <= nu_b + 1`, `b != 1`), the case `k = nu_b + 1` being
   `K \ K_0`. For `design/SOURCE_ERRATA.md` (the editor's file; not written by me).
