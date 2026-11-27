# Draft notes (part 2 of 2)

Continuation of `DRAFT_NOTES.md`: s.0 is the chain of Theorem A, and s.1-s.3 list the 12 notes by category.

## 4. Notes removed or added by the trim (16 -> 12)

| Place in the full version | Label | Note | Action |
|---|---|---|---|
| s02_setup_4.tex:77 | rmk:set-qrrderiv | Mumford's relation and the genus-one values used without a pinned statement | removed with the remark (cut) |
| s04_transport_3.tex:70 | rmk:tra-consistency | the string equation for the model is not stated, proved or pinned | removed with the remark (cut) |
| s09_gamma_4.tex:88 | rmk:gam-sevenvar | sketch only; the formal delta lemma is not in the stationary-phase appendix | removed with the remark (cut) |
| s09_gamma_4.tex:107 | rmk:gam-Mcomplex | not proved (formal stationary phase with a non-integer exponent) | removed with the remark (cut) |
| s10_main_3.tex:68 | rmk:mai-tau | not proved; uniqueness of $\tau$ | removed with the remark (cut) |
| s10_main_4.tex:79 | rmk:mai-analogy | the genus-3 computation is not recorded | removed with the remark (cut) |
| (new) s10_main_2.tex:46 | rmk:mai-graphlocal | genus-4 computation no longer in the evidence tables | added (trim T1) |
| (new) s10_main_2.tex:51 | rmk:mai-necessity | the computations behind the first sentence no longer in the evidence tables | added (assembly) |

The six removed notes are exactly those inside cut remarks (the trim plan counted seven; `rmk:set-qrrderiv`,
`rmk:tra-consistency`, `rmk:gam-sevenvar`, `rmk:gam-Mcomplex`, `rmk:mai-tau`, `rmk:mai-analogy` carried one each). Every
other note of the full version is kept, with unchanged text; their line numbers changed (s01_intro.tex 38 -> 34, 40 -> 36;
s01_intro_2.tex 74 -> 72; s06_model_7.tex 99 -> 89; s09_gamma_4.tex 117 -> 83; s10_main_2.tex 55 -> 41).

If the authors keep decision D8 for the cut items (the 7-variable chart, LMI for complex $M_0$, uniqueness of $\tau$), the
cheapest fallback is one sentence each, marked unproved where the source remark was, with the draftnotes of
`rmk:gam-Mcomplex` and `rmk:mai-tau` restored (TRIM_PLAN_2 s.7 item 1).

## 5. Points recorded without a draftnote (for the authors)

None of these is a gap in the chain of Theorem A.

* **Pins to ratify** (unchanged by the trim). [PPZ15, (8)] (pull-back of $\psi$ classes, used in `lem:set-pushins` (i));
  the PPZ lists §0.2 (1), (5), Def. 0.3, the conventions after Defs. 2.2 and 2.5, "Prop. 2.14, with Props. 2.4, 2.8 and 2.9
  through its proof"; [Eynard:2007kz, Thm 4.7, (4-30)] and (5-15) under "cited but not used"; the reading pins in the sentence
  "The inputs are read with ..." of `sec:int-inputs`.
* **The role of (T)**: (T) identifies the $R$-matrix in Theorems A and C alike. Human check: the proof of Theorem C uses this
  matrix only through `def:set-rows` (b), $R(0)=1$ and (QRR).
* **DNOPS §2.3** prints its action without $1/|\mathrm{Aut}\,\Gamma|$; the note reads it as PPZ's action, which has it.
* **Evidence pointers.** The evidence appendix now keeps only the tables (`sec:evi-graphs`, `sec:evi-normal`,
  `sec:evi-model`). Besides the two draftnoted remarks, the text cites it in `s01_intro_3.tex` l.21 and l.51 (generic pointers)
  and `s06_model.tex` l.40, `s06_model_10.tex` l.48 (`sec:evi-model`, both marked "not used"). The full appendix is in
  `archive_full_162pp/`; an ancillary file is an author decision.
* **Bibliography.** The entry MR1728879 is no longer cited (it was cited only in `rmk:set-qrrderiv`); BibTeX drops it from
  the printed list. BibTeX warns "missing pages in ABDKS25" (published online first).
* **Further cuts, second pass (done).** `prop:con-deg` (with `eq:con-PL`), `lem:con-u0`, `lem:con-weight`,
  `cor:con-frame`, `rmk:con-710`, `rmk:gam-necessity`, `rmk:fsp-classical`, `cor:gau-Gstar`, `cor:mod-basis`,
  `cor:mod-primary`, `lem:mai-coord` (the coordinate-change lemma of the Complements) and the subsections `sec:gau-symp`,
  `sec:sgn-qrr` were removed; none is used by a proof of the chain (`review/trim/CHECK2.md`). All are in
  `archive_full_162pp/`. Whether to keep `lem:mai-coord` is an author decision.
* **Earlier referee points, unchanged.** rf2: (a) translation (d) of the pin of (M) uses the string equation for stable maps
  without a pin; (b) the weights of the virtual normal bundle stay invertible at $\lambda_i=\zeta^i$ (part of the quoted
  definition); (c) in `def:set-rows` the properties listed after "convention:" do not determine $R$, (T) and (QRR) fix it. rf3:
  $\Phi^\circ$ of `lem:gau-poly` (iii) is used nowhere; `def:gau-van` (ii) cites `lem:gau-vanq` (b) ahead of it (not
  circular); some local letters reuse global ones. rf4: three forward references in statements by design
  (`conv:tra-forms`, `lem:tra-conv` (e3), now in the appendix, and `def:tra-dressing` (ii)). rf6: whether
  $\Dn_0(z)\Dn_0(-z)=1$ can be proved without Theorem Gamma is open and not needed. rf7: the stationary-phase appendix
  writes Jacobians with $D$; only `lem:fsp-generic` uses $\Dq$.
* **Letters with several meanings**, each declared where used: $q$, $J$, $x$, $\lambda$, $\mu$, $\xi$.
* **Typography** (119-page build): no overfull hbox is left; the widths were fixed by rewording or by displays.

## 6. Corrections for the design files

These files were not edited: the design, digests and gap notes are outside the assembly.

* **NOTATION.** Hatting matrices $\Lambda_s$, $\Lambda_\varsigma$, bare $\Lambda$ the leg map; HPF operators
  $\mathscr O_e$; the Mellin variable is $\sigma$; $x=s\varpi_0$; moduli spaces $\overline{\mathcal M}$. The constants
  $\mathsf a_j$ of `rmk:gau-Uremarks` are gone with that remark.
* **LABELS.** Add `lem:set-pushins`, `eq:tra-leaf`, and the trim's new labels `sec:set-lrproof`, `eq:set-McycLR`,
  `sec:gau-cycproof`, `sec:tra-convproof`, `sec:cmp`; remove the cut labels (TRIM_PLAN s.4, `lem:mod-vanper`,
  `sec:mai-rmk`, `sec:con-finite`, `lem:con-N`, `thm:con-conn`, the six removed `sec:evi-*` subsections).
* **BIBPLAN.** PPZ15 pins as in s.5; the reading pins (MR2276766, EO, DOSS); DOSS §3.1, §3.5; LP20 §2.4.3 (17), §2.4.4
  (18)-(20); Hock25 also in `tab:sgn-tr`; EO Thm 4.7 (4-30), (5-15) recalled, not inputs; MR1728879 no longer cited.
* **OUTLINE s.0.2.** Item (T) should name Theorem C; the PPZ pin list lacks (5), Prop. 2.8, §0.2 (1) and (8); step C9 uses
  `lem:set-pushins`; row C10 uses PPZ Def. 0.3, Rem. 2.1, (6), Defs. 2.2, 2.5, 2.13, Prop. 2.12 and the dilaton lemma. The
  outline's layout must also follow the trim: section 11 Complements, appendices A-F (A-C new, D-F the old A-C).
* **design/TITLE_ABSTRACT.md** still holds the old abstract.
* **design/SOURCE_ERRATA.md**: some sources carry the symplecticity of the rows over to the conifold expansion, unproved
  (the note proves $R'(-z)^T\eta^{-1}R'(z)=M(\Dn(z)\Dn(-z))^{-1}M^{-1}$); some take the colour-0 amplitude rational before
  (QRR) is used (the note proves it in `thm:con-glue`).
* **DECISIONS D8** (remarks about two pages): after the trim about 2.3 pp of remarks remain; the cut items D8 named are
  listed in s.4.
