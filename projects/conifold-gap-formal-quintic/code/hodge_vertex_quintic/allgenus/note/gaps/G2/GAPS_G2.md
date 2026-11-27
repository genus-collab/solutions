# Gap closing, cluster G2: Theorem U, the graded and canonical gauge, the colours of `D` (section s03)

Author-facing working file of the gap-closing phase (key gaps/G2); not note text. Pieces: this file (s.0 results, s.1 search
record, s.2 setting and local notation); `GAPS_G2_2.md` (s.3 `lem:gau-comm`, s.4 `lem:gau-lift`, `lem:gau-scal`);
`GAPS_G2_3.md` (s.5 Theorem U: existence, the leading-term lemma, the uniqueness clause, the field); `GAPS_G2_4.md` (s.6 graded
gauges and `prop:gau-graded`, s.7 `lem:gau-diag` step 4, s.8 `thm:gau-G` Cor. (a)); `GAPS_G2_5.md` (s.9 `cor:gau-colours`);
`GAPS_G2_6.md` (s.10 evidence, s.11 corrections and notes for the editor, s.12 dependencies, files).
Labels used for results here: PROVED (complete argument written in this folder), VERIFIED (exact computation; range, field or
prime stated), OPEN. Finite checks are evidence only; no proof below uses them. Notation: `NOTATION.md`, `NOTATION_2.md`; labels:
`LABELS.md`, `LABELS_2.md`; new local symbols and new labels are proposed in s.2 and s.11 and must be registered by the editor.
Sources (abbreviations of `design/DIGEST_B.md`): UR s.3-s.5, GI s.1-s.5, GR s.2, NFS item 2, PS s.3-s.4, CS s.1, FQ2 s.1, TH s.1.

## 0. Results

| # | item (outline place) | gap | resolution | where |
|---|---|---|---|---|
| 1a | `thm:gau-U` step 3 | G1: why `A~` is regular semisimple; commutant `= KK[calY]`; `gl_5(KK) = KK[calY] (+) im ad A~` | **PROVED** (elementary, graded; regular semisimplicity is a remark, not needed) | s.3 (a)-(d), (g) |
| 1b | `thm:gau-U` step 4 | `KK[calY] cap gl_5(Q[[u]]) = Q[[u]][calY]` (disjoint supports) | **PROVED** | s.3 (f) |
| 1c | `thm:gau-U` step 5 | G2: the projection formula (trace form) | **PROVED** | s.3 (c)-(e) |
| 1d | `thm:gau-U` step 6 | pole bound of `pr(X)` | **PROVED** (transcribed, checked) | s.3 (f) |
| 1e | `thm:gau-U` step 8 | G3: `pr((calY^j)') = -(j/5)(log Bprod)' calY^j`, `pr(dlog Theta_U calY^j) = (tr dlog Theta_U/5) calY^j` | **PROVED** | s.4.2 |
| 1f | `thm:gau-U` step 10 | G4: uniqueness clause | **PROVED**, in the stronger form "every solution of the gauge equation in `gl_5(Q[[u]])[[z]]` is `P c(z)`", via a new leading-term lemma `lem:gau-lead` | s.5.3-s.5.4 |
| 1g | `thm:gau-U` step 11 | G5: graded clause (existence of a graded gauge; freedom `c(z^5)`) | **PROVED** via the Z/5-degree of the gauge equation (`lem:gau-degree`); no `zeta`, no averaging needed | s.6.1-s.6.2 |
| 1h | `thm:gau-U` step 12 | G7: the field `Q` | **PROVED** (every step is rational; statement over any field `k` of characteristic 0) | s.5.5 |
| 2 | `prop:gau-graded` step 1 [M] | existence of graded gauges, `S A S^-1 = zeta A` | **PROVED** (s.6.2); the `S`-conjugation form is shown equivalent to gradedness (s.6.4) | s.6 |
| 3 | `lem:gau-diag` step 4 | G8: `prod_a p_(0,a) = 1` for `P_can` | **PROVED**, with its exact hypothesis (programme normalisation `b_a(0) = b~_a(0)`, `r = r~`); also a new two-line proof of steps 1-5 from `lem:gau-lead` | s.7 |
| 4 | `thm:gau-G` Cor. (a) | G9: `P_can` is the unique regular gauge with `P Y^C = Y^C_Q` | **PROVED** (uses the uniqueness clause 1f and (G) for `P_can`) | s.8 |
| 5 | `cor:gau-colours` | G6: `D_i(z) = D_0(zeta^-i z)`; colour covariance of `P_can` | **PROVED**: the computation, the covariance of `P_can` (= gradedness) and of the model columns in full; the covariance of the actual columns is PROVED from four named inputs, one of which (the cyclic shape of the actual rows, part of `prop:set-cyclicLR`) is the open item R1 of cluster G1 | s.9 |

Nothing in this cluster is left OPEN on its own account. Item 5 inherits one dependency from cluster G1 (s.9.4: the cyclic shape of the actual rows and, through it, their membership in
the differential ring; s.11 item 2).
Side findings for the editor (s.11): (i) `D_i` lies in `Q(zeta)[[z]]`, not `Q[[z]]`, for `i != 0` (correction to `thm:gau-R0`
(iv), DB D14, PS s.4.1); (ii) the colour argument is on amplitudes and needs no normalisation of the integration constants of the
`u_i` (s.11 item 4); (iii) every regular gauge has `P_0 = gamma frak d Theta_U` (diagonal); (iv) new labels and local symbols.

**Evidence (VERIFIED; s.10).** Independent of the proofs, all exact, every central check with a failing control: the linear
algebra of s.3 at random rational points; the two projection facts of s.4.2 in truncated series (frame (0,0), `u^13`); the graded
gauge built along the proof of s.5 for frames `(0,0)`, `(1/3,-2/7)` (`z^0..z^15`, `u^19` or better), reproducing the sources'
canonical values (`pi(z) = (127/3150) z^3 + (1416472237/725359556250000) z^8`, `p_(1,b)(0) = (8/45, 0, 8/45, -1/15, -1/15)`);
G8; (G) for `P_can` and the G9 control; the uniqueness clause INCLUDING non-graded gauges by a route independent of the
construction (truncated linear system mod `2^61-1`, all five Z/5-degrees, both frames, `k <= 15`); the model colour covariance
over `Q(zeta)` to `t^10`.

## 1. Search record

The task asked first whether any step is written in a source the design missed. Searched (read unless stated): UR s.0-s.5 (s.3 in
full); NFS items 1-14 and s.2 item 2 (the only written audit of Theorem U: "audited line by line by a second AI and re-implemented
along the proof route", no audit text exists); RN (written before Theorem U); `natural_form/verify_theory/` (logs and code only;
`v2_gauge.py` re-implements the proof route in GRADED form and asserts `Sum sg_a = 0` and holomorphy at every order, which is
evidence, not a proof; `v2_run.log` read); `natural_form/critic/` (logs only); GI s.0-s.9; GR s.0-s.5 (paper checks of (G), Prop.
S, the diagonal; it quotes the uniqueness clause, does not prove it); TH s.1 (property (A) quoted); PS s.3-s.4 (s.4.5 says
`D_i(z) = D_0(zeta^-i z)` is "confirmed", with only the `h_w` covariance written); CS s.1 and FQ2 s.1 (the covariance of the
actual rows is stated, and checked numerically, not derived); DIGEST_A D6 (the same); EVID s.2.1-s.2.3 (colours "by covariance",
i.e. by definition, in that note); UNIVERSALITY s.2 and its referee (a different gauge code; free constants only at `k = 0 mod 5`,
evidence); MS2, MS4, MS5 (quote Theorem U as input (U)); NATURAL_FORM s.0-s.1 ("ARGUED"). Grep over the whole package for "regular
semisimple", "commutant", "trace form", "uniqueness clause", "one constant per order", "every regular gauge": no derivation of any
of the steps 1a-1h, 2-5 exists. Hence no item is resolved as FOUND; all are resolved by proofs written here.

## 2. Setting, standing hypotheses, local notation

**Rings.** `k` a field of characteristic 0 (in the note `k = Q`; `Q(zeta)` in s.9); `k[[u]]`, `k((u))`; for `k = Q`,
`k((u)) = KK` of `def:gau-conc`. `' := d/du`. Matrix units `E_ab`, indices `a, b` in Z/5. Gauges are formal series
`P = Sum_(k>=0) z^k P_k(u)`, `P_k` in `gl_5(k[[u]])`, i.e. elements of `gl_5(k[[u]])[[z]] = gl_5(k[[u, z]])`.

**One-pole systems (`def:gau-onepole`).** `A = Sum_a b_a E_(a,a+1)`, `A~ = Sum_a b~_a E_(a,a+1)` with
(H) `b_a, b~_a` in `k[[u]]^x` for `a != 1`, and `u b_1`, `u b~_1` in `k[[u]]^x` (residues `r := (u b_1)(0)`, `r~ := (u b~_1)(0)`);
(EP) equal products `Bprod := prod_a b_a = prod_a b~_a` (in `k((u))`); by (H), `Bprod = beta/u` with `beta` in `k[[u]]^x`.
The gauge equation (written `GE(A, A~)`): `z P' = A P - P A~`. It says that `Y = P Y~` maps solutions of `z d_u Y~ = A~ Y~` to
solutions of `z d_u Y = A Y`. Write `V(A, A~)` for the set of ALL solutions `P` in `gl_5(k[[u]])[[z]]` of `GE(A, A~)` (no condition
on `P_0`); a **regular gauge** (`def:gau-onepole`) is an element of `V(A, A~)` with `P_0` invertible. `V(A, A~)` is a module over
`k[[z]]`: for `c(z)` in `k[[z]]` (constant in `u`, scalar) `z (P c)' = (A P - P A~) c = A (P c) - (P c) A~`.

**Local notation (proposed; to be declared at the top of `s03` as local to `sec:gau-U`).**
* `calB_a^(m) := b~_a b~_(a+1) .. b~_(a+m-1)` (`m >= 0`, indices mod 5; `calB_a^(0) = 1`, `calB_a^(5) = Bprod`): partial products of
  the weights of `A~`; calligraphic B as `Bprod` (NOTATION table 4), which is the case `m = 5`.
* "shift": `gl^(m) := (+)_a KK E_(a,a+m)`, `m` in Z/5; an element `Sum_a xi_a E_(a,a+m)` has shift `m` (UR: "degree m").
* `calY := A~^-1` (UR: `y`; renamed because `y_i`, `y^(n)` are solutions in s03).
* `Theta_U = diag(d^U_a)` (registered in NOTATION table 4), `d^U_0 := 1`, `d^U_(a+1) := d^U_a b~_a/b_a`; `dlog Theta_U :=
  Theta_U^-1 Theta_U'` (UR: `Dl`); `frak d(u) := (det Theta_U(u)/det Theta_U(0))^(-1/5)` (binomial series; `frak d` in
  `1 + u k[[u]]`), so that `frak d' = -(1/5) tr(dlog Theta_U) frak d`.
* `calZ` the lift of `lem:gau-lift` (registered: "UR's lift Pi -> local calligraphic Z").
* `pr` the projection of `lem:gau-comm` (d).
* Z/5-degree of series (s.6): `deg(z^k E_ab) := k + b - a` mod 5; `P^[delta]` the degree-`delta` part.
* `gr := diag(0, 1, 2, 3, 4)` and `T_i := zeta^(i(gr + 1)) = diag(zeta^(i(a+1)))`; the colour operator `Rot_i` (s.6.4, s.9):
  `(Rot_i Y)(u, z) := T_i Y(u, zeta^-i z)` on vectors, `(Rot_i P)(u, z) := T_i P(u, zeta^-i z) T_i^-1` on matrices.

**Programme normalisation (used only in s.7 and s.9).** For the quintic in a frame of the programme family and the model pulled back
by `v(u)` (`prop:gau-quintic`, `prop:gau-equalprod`): `b = (t_c', -1/(u varpi_0^2 (1-u)^3 t_c'^2), t_c', -varpi_0/(1-u),
-varpi_0/(1-u))`, `b~ = (v', -v'/v, v', -v', -v')`, and (N0) `b_a(0) = b~_a(0)` for `a != 1` (values `1, 1, -1, -1` at `a = 0, 2,
3, 4`) and `r = r~ = -1`, because `t_c'(0) = varpi_0(0) = 1` and `v = u + O(u^2)`. (H) and (EP) hold there (`prop:gau-quintic`,
`prop:gau-equalprod`, outside this cluster), and all coefficients lie in `Q[[u]]` (`varpi_0`, `t_c` in `Q[[u]]` with
`varpi_0(0) = 1`; `tau^(-5/4)` in `1 + u Q[[u]]` by the binomial series).
