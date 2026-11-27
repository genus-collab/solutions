# Review rf4: s04_transport (parts 1-4) and s05_reduction (part 1 of 2: verdict, issues, fixes, open items)

Author-facing; not part of the note (key `review/rf4`). Object: `NOTE/sections/s04_transport.tex`, `_2`, `_3`, `_4`,
`NOTE/sections/s05_reduction.tex`. Part 2 (`REVIEW_2.md`): the line-by-line record of what was checked and how, the changes
that affect other sections and the design files, the compile record. Only the five section files above were edited.

## 0. Verdict

**s04_transport: sound-after-fixes. s05_reduction: sound-after-fixes.**

No step of the chain (C7, C8 of OUTLINE s.0.2) is wrong or open. Every all-orders step of the frames, of the convergence lemma,
of Lemmas E, T, F, of the edges in `K_0 (x) K_0`, of the leaf point, of the homogeneous vertex tensor and of Theorem R is
written in sentences and follows from statements proved earlier in the note (s02, s03, earlier parts of s04) or from the pinned
definitions of [PPZ15]. No proof cites a finite computation. No `\draftnote` sits on a step of the chain (before and after the
review). I found one wrong sub-statement (F2), one scope error in the dependency structure (F1), one imprecise identification
with the literature (F3), one misattributed hypothesis (F4), one remark whose "known only by computation" claim is false because
the note proves it (F5), two forward references inside proofs (F7, F8), and smaller points (F6, F9-F11). All are fixed in the
section files; none changes a statement used by another section (see `REVIEW_2.md` s.3).

## 1. Issues found and fixed (location = file after the fix; "old" quotes the reviewed text)

**F1 (scope error; the leaf identity used outside the hypotheses of the lemma that proved it).**
`lem:tra-F` is stated for a dressing with (Sp) only, and so are `def:red-flatsum` and `eq:red-vertex`; but `lem:tra-F` Steps 1
and 3, the paragraph after `def:red-flatsum`, the paragraph after `def:tra-edge` and the remark after `lem:tra-conv` took the
identity `qt_T + z1 = d T_(>=2) = ttr_(R^vec)(z)`, and the very definition of `T_(>=2)`, from `lem:tra-T` (c), whose
hypotheses are (Sp) **and (A)**. The identity is pure algebra (grading of `P_0`, `(P_0^vec)_00 = d`, `R^vec = 1 + O(z)`).
*Fix:* a new paragraph "The leaf point" after the divisibility paragraph of `def:tra-edge` (s04_transport_4) derives, for
every dressing with (Sp) and without using (Sp) or (A), the auxiliary display `eq:tra-leaf`:
`q_T = -(z/d)(R^vec)^-1(z) 1 = -(z/d) 1 + T_(>=2)`, `T_(>=2) := d^-1 (z1 - z (R^vec)^-1(z) 1)` in `z^2 HH^u_+`, hence
`qt_T = -z1 + d T_(>=2)`, `d T_(>=2) = ttr_(R^vec)(z) = z1 - z d (P^vec)^-1(z) 1`. `lem:tra-T` (c) (statement and proof),
`lem:tra-F` (Steps 1, 3), `eq:red-vertex` and the paragraph after it, `thm:red-main` Step 3 and the remark after
`lem:tra-conv` now cite `eq:tra-leaf`. Content unchanged; the dependency is now correct.

**F2 (wrong sub-statement; `lem:tra-T` (c), s04_transport_4).** Old: "`d T_(>=2) = y° + delta_P H` with `y°` in `K_0` and
`y°_(k,a) = 0` for `k <= 1`". With the note's notation `q_(k,b) = [z^k] q^b` this includes `(k,a) = (0,1)`, and
`y°_(0,1) = -delta_P`, which is non-zero in general (`delta_P = u^2/30 + ..` for `P_can`; G3R_2 s.7 records exactly this:
"`t'` has the `z^0 H`-component `-delta`"). *Fix:* "`y°` in `K_0`, whose frame coordinates `y°_(k,a)` (`a != 1`) vanish for
`k <= 1`"; the proof now says why (`y° = d T_(>=2) - delta_P H` has the same components with `a != 1` as `d T_(>=2)`). The uses
(`thm:red-main` Step 4, `lem:tra-conv` (b3), (e)) only ever used the frame coordinates, so nothing downstream changes.

**F3 (imprecise identification with [PPZ15]; `def:tra-edge`, last two sentences).** Old: "The theory `Omega^P := P.Omega^mod`
is the R-matrix action of [PPZ15, Def. 2.2] for `R^vec`, ..., with the translation `ttr_(R^vec)` of [PPZ15, Def. 2.13] ...,
preceded by the constant frame change `P_0`. Its leg-free graph sum, in the form used in this note, is the flat sum". Read
literally, the leg-free graph sum of `(P_0)_*(R^vec.Omega^mod)` is `Sum_Gamma Cont_Gamma[E_(R^vec), q_(R^vec); Omega^mod]`,
which differs from the flat sum graph by graph by `Prod_v d^(2g_v-2+n_v) = d^(2g-2)`: the source description (RD1 s.1.5,
DC1 D12) omits the rescaling of the classes by `d^(2h-2+n)` ("absorbed in the vertex TFT", RD1 s.1.5). The symbol `Omega^P`
is used nowhere else in the note. *Fix:* `def:tra-edge` now names "the flat sum of `P.Omega^mod_(v(u))`" (the only graph sum
attached to `P` that the note uses; it is `def:red-flatsum`) and refers to the exact relation, stated after `eq:tra-leaf`:
`E_P = E_(R^vec)` and `qt_T = q_(R^vec)` are the edge and leaf point of `R^vec` in the sense of `def:set-graphsum`, i.e. of the
unit-preserving action `R^vec.Omega^mod` [PPZ15, Defs. 2.2, 2.13]; the vertex tensors of the flat sum are those of this
action multiplied by `d^(2g_v-2+n_v)`; for `P_0 = 1` the flat sum is its leg-free graph sum. This was already the content of
the old second paragraph; only the inaccurate first sentence went. No proof used it.

**F4 (misattributed hypothesis; `thm:red-main` Step 10, s05).** Old: "The hypotheses `v(0) = 0` and `v'(0) != 0` enter
through `lem:tra-conv` (e) and `lem:tra-T` (d)". The proof of `lem:tra-T` (d) uses only `v` in `u Q[[u]]` (every term
`k >= 1` of `T_P - p_00 v` has `v`-exponent `>= 2`); `v'(0) != 0` enters only through `lem:tra-conv` (e) (`ord_u v = 1` for
negative exponents and `delta/v` in `u Q[[u]]`) and makes `T_P` equal to `u` times a unit. *Fix:* the sentence now says so.

**F5 (remark claiming "known only by computation" for a proved formula; `rmk:red-leaf`, s05).** Old: "The value
`t_c - vhat = u^2/30 + O(u^3)` is known only by computation" with the `\draftnote` "No derivation of the `u^2`-coefficient of
`delta_(P_can)` is written; its frame independence ... uses the leaf coefficient `-p_(1,4)(0)` ...". The note proves
`d = 1/p_(0,0)` (`lem:tra-T` (d)), `p_(0,0) = tau^(-1/4)/varpi_0` (`lem:gau-diag`, the G2 s.7 result, already used in
s10_main_3 l.98) and `T_(P_can) = t_c` (`prop:tra-can`); with `varpi_0 t_c = varpi_1` (`def:set-frames`) this gives
`delta_(P_can) = d t_c - v = tau^(1/4) varpi_1 - u tau^(-5/4)`, which is the formula of the sources (TRF s.2 "new remark";
DC1 B3; DC1_5 CV14) and is manifestly frame independent; with `tau = 1 - (4/9) u + O(u^2)` (`def:set-tau`) and
`varpi_1 = u + (7/10) u^2 + O(u^3)` (`lem:set-frobbasis`) its `u^2`-coefficient is `7/10 - 1/9 - 5/9 = 1/30`. *Fix:* the
remark now writes this derivation, keeps "not used", and points to appendix C only for the next coefficient; the
`\draftnote` is removed (it was outside the chain). Sanity check, evidence only: my exact series computation of
`tau^(1/4) varpi_1 - u tau^(-5/4)` gives `u^2/30 + (127/3150) u^3 + ..`, the value recorded in c_evidence_3 l.32.

**F6 (self-contradictory clause; `lem:tra-frames` (c), s04_transport_2).** Old: "The map ... is bijective. In the graded setting
its image is the set of finitely supported families" (G3R K1 asked for the graded clause; the merge kept "bijective" for both
settings). *Fix:* "is injective. In the u-adic setting it is bijective; in the graded setting its image is the set of finitely
supported families in `Q[v]`, and ...". The proof already covered both cases.

**F7 (forward reference inside a proof; `lem:tra-coeff`, s04_transport_2).** Old proof ended "The basis statement for
`K_0 cap V_w` is `cor:mod-basis`" (s06). Not circular (`cor:mod-basis` is deduced from `lem:tra-frames`), but rule 1/the task
forbid forward pointers inside proofs. *Fix:* moved after the proof as a sentence saying where the basis is stated and that it
is deduced from `lem:tra-frames`.

**F8 (forward reference inside a proof; `lem:tra-conv` (c), s04_transport_3).** The statement of (c) included "the contraction
of `lem:tra-E`" and its proof ended "for the contraction of `lem:tra-E` see its proof" (a later lemma). *Fix:* the statement
lists the one-variable pairings only and adds, in parentheses, that the contraction of `lem:tra-E` is treated in the proof of
that lemma, in `Q((psi))[[u, psi']]`; the forward pointer is gone from the proof. (`lem:tra-E` Step 1 proves its own
convergence and uses only `f` in `HH^u_-`, which for `f = z^-k J^C` is the first sentence of (c): no circularity before either.)

**F9 (citation list; s04 introduction).** Old: "the definitions of the R-matrix action in [PPZ15, Def. 2.2 and Def. 2.13]"
although the section also cites [PPZ15, Prop. 2.12] (the `T_b`). *Fix:* "[PPZ15, Def. 2.2, Prop. 2.12, Def. 2.13], recalled
after `def:tra-edge`", and one sentence on where the dilaton equation enters (`lem:tra-conv` (f), and through it the last
statement of `lem:tra-F`). Statement numbers confirmed on the public arXiv HTML rendering of 1303.1043 (`REVIEW_2.md` s.1.4).

**F10 (notation; all five files).** The unit was `\mathbf 1` in s04-s05 only; s02 (`def:set-graphsum`, `lem:set-dilaton`), s07,
appendix B and the comment of the macro `\ttr` write `1` (`z1`, `zR^-1(z)1`). Rule 6 (one notation): harmonised to `1`
(62 occurrences); where a bare `1` could be misread I wrote `z^6\,1`, `z^2\,1`, `\eta\,1` and `(P_0^vec)^-1 1 = d^-1 1`.

**F11 (wording, precision; no change of content).** "A shorter route" -> "A shorter argument" (rule 5 lists "route");
"the homogeneous form, the form in which the flat sum is computed" -> "..., which makes sense for every dressing with (Sp): its
leaf `d T_(>=2)` has psi-degree at least 2" (no hint at computations); `rmk:red-hyp` (v) "only finitely many correlators are
tested" -> "the Model Lemma is used for finitely many correlators only"; `rmk:red-hyp` (ii) now mentions the derivatives of the
gap used in (iv); `prop:tra-can` now cites `thm:gau-G` (a), `lem:gau-vanq` (c) (the shape that gives (A)) and the places where
(N0) holds (`def:set-frames`, `prop:gau-equalprod`), and `lem:tra-T` (d) for `d = 1/p_(0,0)`; `lem:tra-conv` (f) cites
`lem:set-dilaton` (ii) (the flat unit of the model, needed to apply (i)); `def:red-flatsum` cites `def:set-modcohft`,
`def:set-graphsum` for `Fbar_h` (where it is defined); the s05 introduction cites `sec:set-graph`, `sec:tra-dress` for the
[PPZ15] definitions; `conv:tra-forms` no longer says that `Q((z))[[u]]` "contains every object of this section" (the edges are
two-variable objects, treated in `lem:tra-twoslot`, `lem:tra-E`).

## 2. Open items (none on the chain; none is a gap of Theorem A)

* O1. `rmk:tra-consistency` keeps its `\draftnote` (the string equation of the model is not stated, proved or pinned). The remark
  says it is not used; this is O2 of OUTLINE s.0.2 and stays so. The sources (RD1 s.2.1) only use the string equation, they do
  not prove it; a one-paragraph proof from `Omega_(h,n+1)(.., 1) = p^* Omega_(h,n)` and the comparison of psi classes is
  standard but not in the sources, so I did not write it (rule 1).
* O2. Local letters that also occur elsewhere with another meaning: `pi_r = 1/(4^r r! (4r)!)` (s04) against `pi_m` (s08,
  NOTATION_2 table, reduction coefficients); `A` (the coefficient ring) against `A_(m_2)` in `lem:tra-conv` (f2); `T_m` (the
  m-th term in (f2)) against `T_Gamma`, `T_P`. All are declared locally; harmless; the editor may rename `pi_r`.
* O3. Forward references in statements (not in proofs), by design: `conv:tra-forms` -> `def:mod-weights` (s06; a consistency
  remark about the same weight spaces), `lem:tra-conv` (e3) -> `def:red-ML` (s05; used as a hypothesis, OUTLINE_5 and G3R K10),
  `def:tra-dressing` (ii) -> `rmk:mai-annih` (s10; an example). Not circular.
* O4. `def:tra-dressing` example (i) for a general self-dual one-pole system leaves the pull-back `v_A(u)` of the model
  implicit (it is made precise in s10, `cor:mai-univ`). Example only; not used.
* O5. Labels: the new display label `eq:tra-leaf` carries the section prefix, which `LABELS_2.md` ("Equation labels fixed now
  (others: section prefix, free)") allows without registration; the editor may add it to that list. No other label was added,
  renamed or removed. I did not edit design files.
