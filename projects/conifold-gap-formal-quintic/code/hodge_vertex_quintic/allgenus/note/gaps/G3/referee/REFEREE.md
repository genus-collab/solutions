# Referee report on gap cluster G3 (s04, transport), part 1 of 2: verdict, items G3-1, G3-2a, G3-2b, G3-2c

Author-facing; not part of the note (key `gaps/G3/referee`). Object, read in full: `NOTE/gaps/G3/GAPS_G3.md`, `GAPS_G3_2.md`,
`GAPS_G3_3.md`, `GAPS_G3_4.md`, the scripts `g3_frames.py`, `g3_series.py`, `g3_decomp.py`, `g3_mat.py`, `g3_lemE.py` and their six
logs. Part 2 (`REFEREE_2.md`): items G3-2d to G3-3, the conventions and observations of `GAPS_G3.md` s.2-s.3, my computations, the
consolidated corrections K1-K10, files. Statuses as in the binding rules. No verdict below rests on a computation; the computations
of s.12 (part 2) are evidence with failing controls. Citations "FL Step n" = `GAPS_G3_2.md` s.2 Step n; "CV (x)" = `GAPS_G3_3.md`,
part (x) of the proposed `lem:tra-conv`; "LE Step n" = `GAPS_G3_4.md` s.2 Step n.

## 0. Verdict

**Overall: accepted-with-corrections.** Every item of the cluster is closed; I found no gap, no wrong statement in a proof of
record and no hidden hypothesis that the prover does not name. The corrections K1-K10 (part 2, s.13) are a misstated clause of
the graded bijection (K1), a label (K2), one proof step that covers only half of a claim (K3, I supply the missing argument), one
wrongly named ring in a side remark (K5), one application step that must be written (K6), an attribution (K9), and wording. None
changes a resolution; no item becomes OPEN.

| item | prover | verdict | corrections |
|---|---|---|---|
| G3-1 `lem:tra-frames` step 4 [S] (with G3-1x, the replacement of step 3, and the trivial-frame erratum) | PROVED | accepted-with-corrections | K1, K8 (b) |
| G3-2a `lem:tra-conv` (a) | PROVED | accepted | - |
| G3-2b `lem:tra-conv` (b) | PROVED | accepted-with-corrections | K7, K10 |
| G3-2c `lem:tra-conv` (c) | FOUND (decay) + PROVED (ring) | accepted-with-corrections | K2 |
| G3-2d `lem:tra-conv` (d) | PROVED | accepted | - |
| G3-2e (added) vertex resummation | PROVED | accepted | K10 |
| G3-2f (added) dilaton re-expansion, `d(0) = 1` | PROVED | accepted-with-corrections | K3, K4, K10 |
| G3-2g (added) two-slot frame expansion | PROVED | accepted | - |
| G3-3 `lem:tra-E` step 2 (O-8) | PROVED | accepted-with-corrections | K5, K6 |

No item is CITED: the cluster quotes no literature statement (`GAPS_G3.md` s.0, "No quoted theorem is used"); I confirmed that
none of the four files invokes a theorem outside the note's own labels. So there was no full text to open.

## 1. What was read, what was re-run

* Design: `OUTLINE_5.md` (s04, s05), `NOTATION.md`, `NOTATION_2.md` (tables 5-6, s.9 dictionary, s.10), `LABELS.md` (s02-s05 rows),
  `LABELS_2.md` (`cor:mod-basis`), `OUTLINE_3.md` (`def:set-modcohft`, `lem:set-dilaton`), `OUTLINE_4.md` (`lem:gau-diag`,
  `lem:gau-vanmod`, `lem:gau-vanshape`), `OUTLINE_6.md` (`def:mod-weights`, `thm:mod-lattice`, `cor:mod-basis`), `design/OUTLINE_REVIEW.md`
  (R1, R7, R8), `design/DECISIONS.md`; digests DC1 (D8-D11), DC1_4 (A1-A4, B1-B5), DC1_6 (s.6, O-3, O-5, O-8), DB_3 (S22-S24).
* Sources: `theory/REDUCTION_DRAFT.md` s.1 (RD1), `REDUCTION_DRAFT_2.md` s.3.1-s.3.2 (RD2), `theory/referee/REFEREE.md` s.1-s.2 (TRF),
  `gauge/GAUGE_IDENTITY.md` s.0-s.5 (GI) and its log `gi_relations_ab_0_0.log` (the `z^3, z^8, z^13, z^18` terms of `P_30(0,z)`),
  `gauge/GAUGE_IDENTITY_HAT.md` s.0-s.4 (GH), `model/MODEL_SINGULAR_2.md` s.6-s.7, s.10 (MS2), `stage2/dict/DICT_PROOF.md` Cor. 4 and
  s.6, `stage2/flat/S15_PROOF.md` s.1-s.4, `stage2/STAGE2_SYNTHESIS_2.md` s.5 (its item O8).
* Re-run (copies in my scratch folder, `py -I -B`): `g3_frames.py 40`, `80`; `g3_decomp.py 8 24`, `10 30`; `g3_lemE.py 6 5`, `8 6`.
  All six outputs are byte-identical to the logs in `NOTE/gaps/G3/`. The logs say what the table of `GAPS_G3_4.md` s.4 says, with one
  nuance: "stop changing exactly at the predicted truncation" holds in the `10_30` run; in the `8_24` run the change stops earlier
  than the bound (`K = 6` resp. `9` against `8` resp. `11`), i.e. the prediction is an upper bound, which is all that is claimed.

## 2. Item G3-1: `lem:tra-frames` step 4 [S] (`GAPS_G3_2.md` s.1-s.2)

**2.1 Statement.** The proposed statement (a)-(d) is the right one for the note: two coefficient settings (`A = Q[v]`, polynomial
in `z`; `A = Q[[u]]`, `v = v(u)` in `u Q[[u]]`, power series in `z`), the direct sum `HH^A_+ = K_0^A (+) A{z} H`, the frames with their
exponents, the frame coordinates as a bijection onto all families (u-adic), and the substitution statement (d).

**2.2 Line-by-line check.**
* *Step 0.* From `I^C_m = Sum_i (-1/(4z))^i/i! (-v/z)^(m+4i+1)/(m+4i+1)!` (residue of `x^m e^(-W_v/z)`): the `z^-n` term has `n = 5i + m + 1`,
  `v`-power `n - i`, sign `(-1)^n`; with `0 <= nu_b + 1 <= 4`, `i = floor(n/5)`. With `Pi_0(z^m H^b) = (-1)^(m+1) [z^-m] Y^C_b` (from the
  definition `Pi_i(q) = -[z^i] <Y^C(-z), q(z)>`; I re-derived the sign) one gets `Pi_0(z^m H^b) = -eps_b v^(m-q)/(4^q q! (m-q)!)` for
  `b + m = 1 mod 5` (`nu_b + 1 = 1 - b mod 5`). (0.1)-(0.3) and the five examples are correct (`Pi_0(z^5 H) = v^4/96`,
  `Pi_0(z^6 1) = -v^5/480` checked by hand).
* *Step 1.* The shift rule and `(Theta s)_i = s_i + Sum_(q>=1) pi_q v^(4q) s_(i+5q)` are correct. Graded: `N` lowers the support by 5, so
  `N` is nilpotent on finitely supported sequences. u-adic: `N` raises the u-order of every entry by at least 4, uniformly in `i`. In
  one line: the module `prod_(j>=0) Q[[u]]` is complete for the filtration `u^m prod Q[[u]]` and `N` maps the `m`-th step into the
  `(m+4)`-th, so `1 + N` is invertible with inverse `Sum_r (-N)^r`; the prover's argument (summability of `((-N)^r s)_r`, `N` commuting
  with the sum, telescoping on both sides) is this, written out. The weight bookkeeping (`v^(4q) s_(i+5q)` of weight `w - 4i`) is right.
  `Pi_i` of finite type: from (0.3), modulo `u^N` only `k < i + 5N/4` matter. Correct.
* *Step 2.* `pr(q) = q - sigma(q)(z) H`, `sigma = Theta^-1 (Pi_i(q))_i`: `Pi_i(pr q) = 0` uses `Pi_i(s(z) H) = (Theta s)_i` for an infinite
  `s`, which holds because `Pi_i` is the u-adically convergent sum of Step 1 (ring product of `GAPS_G3.md` s.2.1). Directness from the
  injectivity of `Theta`. Graded homogeneity of `sigma(q)` from Step 1. Correct.
* *Step 3.* Support of `sigma(z^k H^a)` in `[0, k)` (vanishing for `i > k` by the shift rule, at `i = k` by (0.2), and `Theta^-1`
  preserves supports `[0, L)`); uniqueness; no frame with `a = 1`; exponent `(4(k-j) + nu_a + 1)/5 >= 1` for `j <= k - 1`, `nu_a >= 0`.
  (d) by uniqueness after the ring map `Q[v] -> Q[[u]]`. The examples, including `frame(6,0)` (`f_0 = v^5/480 - v^5/96 = -v^5/120`),
  are correct. Correct.
* *Step 4 (the [S] step).* Summability from the u-order `(4(k-j) + nu_a + 1)/5` of the `z^j H`-coefficient; `n(c)` in `K_0` by the
  summation lemma; injectivity/surjectivity by (a); the projection formula `sigma(q)_j = q_(j,1) - Sum_(k>j, a!=1) q_(k,a) f_j^(k,a)`. The
  `(u, z)`-adic remark (`frame(k,a)` in `(u,z)^floor(4k/5)`) is right. Correct; this is the missing step of DC1 A1 step 4, with its ring.
* *Step 5, Corollaries (i)-(ii).* Correct (`Pi_i(H) = 0`, `i >= 1`; `Pi_0(H) = 1`; `Pi_0(-z 1) = v`, `Pi_1(-z 1) = Pi_0(H^0) = 0`).
* *Defect (K1).* Statement (c), graded clause: "onto the families with finitely many non-zero members, each homogeneous" is not
  right for `K_0^A = K_0 cap H (x) Q[v][z]`, whose elements need not be homogeneous; the image is the finitely supported families in
  `Q[v]`, and on `K_0 cap V_w` the families whose `(k,a)`-member is a rational multiple of `v^((w - 4k - nu_a - 1)/5)`. (Since each `Pi_i`
  maps `V_w` to `Q v^((w-4i)/5)`, `K_0 cap HH^gr_+ = (+)_w K_0 cap V_w`, so the two forms are equivalent.) The proof is unaffected.

**2.3 G3-1x (replacement of the source's step 3).** Accepted. RD1 s.1.3 ("look at its top power") is valid for polynomials only;
modulo `u^N` an element of `Q[[u, z]] H` is still a power series in `z`, so even an order-by-order reading has no top power. The
injectivity of `Theta` is the correct replacement (equivalently: if `m` is the least u-order among the `s_j`, then `(Theta s)_i = s_i`
modulo `u^(m+1)`, so `Theta s = 0` forces `s = 0`).

**2.4 The trivial-frame erratum.** Confirmed, three ways. (i) The proof (Step 3): at `j_0 = k - nu_a - 1` the recursion has no
correction because `j_0 + 5q >= k + 4 - nu_a > k` (`nu_a <= 3`), and `f_(j_0) = eps_a v^(nu_a+1)/(nu_a+1)! != 0`; for `k < nu_a + 1` no index
`j < k` is admissible. (ii) Directly: `Pi_i(z^2 H^2)`, `Pi_i(z^3 H^2)` vanish for all `i` (`2 + m` is never `1 mod 5` for `m <= 3`).
(iii) `r3_frames.py` (independent route, s.12): for `k <= 60`, at `v = 1` and `v = 2`, the trivial frames are exactly the ten
`{(k,a) : k < nu_a + 1}`. The prover's own control count is consistent (`g3_frames_80.log`: 314 non-trivial among the 320 frames with
`k >= 1`, i.e. six trivial ones with `k >= 1`). One more source point for the editor: DC1 D10 (MS1 s.1: "`tau_k(H^b)` lies in `K` iff
... `b = 2, k <= 4` ...") already implies the ten; the digest is internally inconsistent (D9 against D10), so `design/SOURCE_ERRATA.md`
should record RD1 s.1.3 and DC1 D9 as incomplete lists.

**2.5 Verdict G3-1: accepted-with-corrections** (K1; K8 (b) recommends stating `Theta` once as a lemma, since Step 2, Step 4,
G3-2g and Corollary (i) all use it).

## 3. Item G3-2a: `lem:tra-conv` (a) (CV (a), FL Step 4)

The u-order bound of the `z^j H`-coefficient of `frame(k,a)` at `v(u)` is FL Step 3 with `ord_u v >= 1`; summability for arbitrary
coefficients in `Q[[u]]` and the decomposition of infinite elements are FL Step 4; the two-slot version is G3-2g. The prover's T4
(3312 exact u-orders, `k <= 90`) and my U2-U3 (s.12) agree with the bound, which is sharp (`r3_uadic` U3: the last frame index that
changes `[u^m] sigma_j` equals `j + floor((5m-1)/4)` in all eight cases tested). **Verdict: accepted.**

## 4. Item G3-2b: `lem:tra-conv` (b) (CV (b), s.2 (b1)-(b4))

* (b1) The multilinear extension through the finitely many flat slots with `Sum k_i <= D = 3h - 3 + n` is the correct reading of a
  vertex contracted with power-series insertions (the class has non-negative degree; `dim Mbar_(h,n) = D`). Finite type in each
  slot with `F_N = {(k,a) : k <= D}`, `M_N = N - e_0`. Correct.
* (b2) A frame with `kappa(beta) > D` meets a slot `k_i <= D` only through an `H`-coefficient of u-order `>= (4(kappa(beta) - D) + 1)/5`;
  the multi-family over `B^n` is summable; the iterated application of the summation lemma equals the unordered sum. Every term is
  a correlator of homogeneous elements of `K` (frames, `H`), so `ML_K(h)` applies termwise and no `z^j H`, `j >= 1`, is isolated.
  Correct.
* (b3)-(b4) Leaves (`L = d q_T + z 1` in `K^u`, `L = t' + delta_P H`, `t'_(k,a) = 0` for `k <= 1`) and graph terms (edge by edge, each
  edge expanded by G3-2g, the bound of (b2) at both ends) are correct; this is the precise form of `thm:red-main` steps 1 and 6.
* *Defects.* (K7) The parenthesis "homogeneity is not needed for the convergence" is true, but the written bound
  `e_0 = (2 - 2h - n) ord_u v` uses homogeneity (`wt >= -1`); without it, `e_0` must be the least `v`-exponent of the finitely many flat
  correlators with `Sum k_i <= D`. Say which form is printed. (K10) The proof uses `def:set-modcohft` (Laurent-polynomiality, the
  degree bound, and in (b2) for the termwise use of `ML_K` the homogeneity of frames): the printed resolution should read "PROVED,
  given `def:set-modcohft`", whose homogeneity carries its own `\draftnote` (OUTLINE_3; R7 list).
* **Verdict: accepted-with-corrections** (K7, K10; wording only).

## 5. Item G3-2c: `lem:tra-conv` (c) (CV (c); `GAPS_G3.md` s.2.1)

* *The FOUND part (decay).* Checked in the source: GI s.2, paragraph "Lemma V (model; PROVED)", gives `I^C_m` in closed form as the
  formal residue and states that `Y^C_b` is a sum of monomials `v^j z^-n`, `n = 5i + c_b + 1`, `5j = 4n + c_b + 1`, and that
  `Y^C(v(u), z)` lies in `Q[z^-1][[u]]^5`; DB S23 (= `lem:gau-vanshape`) digests it with the closed form `y_b(t)`. The argument is
  complete (a direct reading of the double series). The exponent `j = n - floor(n/5)` is exact, so "`O(v^(n - floor(n/5)))`" of RD2
  s.3.1 is an exact monomial, as the prover says. `r3_frames.py` [Y] re-derives the closed form from the residue definition and checks
  the model system `z dY_b/dv = bmod_b Y_(b+1)` and `Y^C(0) = -e_1` (`n < 80`, `v = 1, 2`).
* *Label (K2).* `GAPS_G3.md` s.0, s.1 (i) and CV (c) write "GI s.2 Lemma V (= `lem:gau-vanshape`)". By NOTATION_2 s.9 and LABELS,
  GI's Lemma V is `lem:gau-vanmod` (DB S22); the shape and closed form used here are `lem:gau-vanshape` (DB S23, sources GI s.2,
  GH s.1-s.2). Cite `lem:gau-vanshape` for the closed form and `lem:gau-vanmod` for `Y^C` being the residue solution.
* *The PROVED part (the ring).* The two-sided description of `Q((z))[[u]]` is right (an element `Sum_N u^N a_N(z)`, `a_N` in `Q((z))`,
  is `Sum_k c_k z^k` with `c_k -> 0` u-adically as `k -> -oo`, and conversely); the Cauchy product converges at each `z^k` and is
  finite modulo each `u^N`; the splitting `( )_(>=0) + ( )_(<0)` is coefficientwise; `z -> -z` and `Res_z` are coefficientwise, and
  `Omega_s` is antisymmetric (`Res_z h(-z) = -Res_z h(z)`, `eta` symmetric). `Pi_i` on `HH^u_+`, `(P z^-k J^C)_(>=0) = z^-k (P J^C)_(>=k)`,
  and `Omega_s(f, P^-1 g) = Omega_s(P f, g)` (from (Sp): `P(-z)^T eta = eta P(z)^-1`) are ring identities there. I checked each.
* *The correction of `conv:tra-forms`.* Right: `H[[z]]` is not contained in `H((z^-1))` (the intersection is `H[z]`); the objects of
  s04 (`Y^C(v(u))` with infinitely many negative powers, `q_T`, `E_P` with infinitely many positive powers, `P_can(0, z)` non-polynomial,
  GI log: non-zero `z^3, z^8, z^13, z^18` terms of `P_30(0, z)`) need `H (x) Q((z))[[u]]`, which is GI's ring `L[[u]]` (GI s.1).
  The graded setting with Givental's `H((z^-1))`, `HH_+ = H[z]` is consistent.
* **Verdict: accepted-with-corrections** (K2 only).
