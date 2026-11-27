# Review rf5: section s06 (sections/s06_model.tex, s06_model_2.tex .. s06_model_10.tex)

Author-facing review file; not part of the note. Continued in `REVIEW_2.md` (statement-by-statement record, open items,
changes affecting other sections, compilation).

## 0. Verdict: **sound-after-fixes**

Every definition, statement and proof of s06 was checked line by line against OUTLINE.md s.0 and s.0.2 (row C9), OUTLINE_6.md,
NOTATION(_2), LABELS_2, BIBPLAN, DECISIONS (Q4), GAPS_SYNTHESIS_2 s.3 and GAPS_SYNTHESIS_3 (the OUTLINE_6 items), the gap note G4
(`GAPS_G4.md`, `_2` .. `_5`) with its referee reports (R-1 .. R-8), and against the definitions of s02-s05 that s06 uses
(`def:set-model`, `prop:set-modrows`, `def:set-M0`, `lem:set-M0sympl`, `def:set-modcohft`, `lem:set-modhom`, `lem:set-dilaton`,
`def:gau-van`, `lem:gau-vanshape`, `conv:tra-forms`, `def:tra-tower`, `def:tra-slice`, `lem:tra-frames`, `lem:tra-coeff`,
`lem:tra-conv`, `def:red-ML`). The quoted statements were re-read in the arXiv sources of the pinned versions (downloaded to a
scratch folder): EO math-ph/0702045v4 (s.4 Defs. 4.1-4.3, Thms 4.1-4.7, Cor. 4.1; s.5 (5-1)-(5-32), Thm 5.1), DNOPS 1509.06954v2
(s.2.3, s.2.5, s.2.6, s.4 Thm 4.1, Prop. 4.3, s.8 (8-1) and the sentence before Thm 8.1), DOSS 1211.4021v1 (s.2.3, s.3.1, s.3.5,
Thm 4.1), ABDKS 2212.00320v3 (s.1.1 (1.1), s.7 Thm 7.1, (7.6)-(7.12)), Hock 2502.18115v2 (Thm 1.1, Cor. 1.2, (1.1), (2.1)).
Equation numbers of EO s.4 and s.5 were recounted from the source (equation / eqnarray counter): (4-3), (4-9), (4-10), (4-20),
(4-23)-(4-30), (5-1), (5-4), (5-13), (5-14), (5-15), (5-27)-(5-30) are as cited.

The mathematics of s06 is correct and complete modulo the quoted inputs: the identification with topological recursion (DOSS
Thm 4.1 at `f = 4 W_v`, `a = (0,0,0,0,4v)`, transferred to `x = W_v` by `lem:mod-scale`), `prop:mod-Rtr` to all orders,
`lem:mod-xfun`, the leg dictionary, the residue identity, the lattice theorem, L-regularity, the Gap Theorem by the closed formula of
ABDKS Thm 7.1 and the exact identity `Sigma_h = -4^h B_2h/(2h)` (binomial series and Lagrange-Buermann), Theorem V with the bases
`(h,a) = (0,1), (0,2)`, `h = 0`, `h = 1` (local formula, `-1/60` per branch point) and `h >= 2` (without the variational formula), and
the Model Lemma on `K` in every genus. Topological recursion is used only on the model curve (and on `x = kappa W_v`, the same curve
up to scaling); the Airy curve and `x = p + 1/p` appear only as sign anchors used in no proof. No proof relies on a finite
computation; the only computed facts (`def:mod-pair`, `rmk:mod-offK`) are labelled as computed and unused. **No `\draftnote` was
added; no step of the chain of Theorem A is left open in s06.** The only `\draftnote` of s06 is the one in `rmk:mod-hock`
(outside the chain, prescribed by OUTLINE_6).

The defects found were: one wrong citation (a section of DOSS), one misquoted theorem (EO Thm 4.7), an inexact list of quoted
inputs (incomplete for DNOPS; two items listed that no proof uses), one wrong item reference, three references in proofs to
later statements, an imprecise justification of the DNOPS compatibility condition, a notation clash (`Lambda_kappa` versus the leg
map `Lambda`), several notation and wording slips, and a misparaphrase in a remark. All are fixed in the section files.

## 1. Issues found and fixed (location; what; fix)

1. **Wrong citation** (`s06_model_2.tex`, paragraph after `lem:mod-laplace`). "this is the prescription of [DOSS, s.2.3] (the
   singular part `dw dw'/(w - w')^2` gives the Kronecker delta, the regular part is transformed)". DOSS s.2.3 is "Expressions in
   terms of graphs" (leaf, dilaton-leaf, edge and vertex weights of the Givental side). The decomposition `B^(i,j) = delta_ij dz dz'/
   (z - z')^2 + regular part` is DOSS s.3.1, and its Laplace transform `Bcheck^(i,j)(u,v) = delta_ij uv/(u + v) + ..`, with the
   functions `f^i_d(u,j) = delta_ij (-1)^d u^d - ..` that carry the Kronecker delta, is DOSS s.3.5 (read in the v1 source). Fixed to
   `\cite[\S3.1, \S3.5]{DuninBarkowski:2012bw}`. This is a new pin for BIBPLAN (REVIEW_2 s.3).
2. **Misquoted theorem** (`s06_model_4.tex`, proof of `cor:mod-trF`). "This is the case `n = 0` of the dilaton equation [EO,
   Thm 4.7 (4-30)]". EO state Thm 4.7 "for `k >= 1`" and remark after it that for `k = 0`, `g >= 2` it holds by definition if
   `W_0^(g) := -F^(g)`. Rewritten: only the definition `(2-2h) F^TR_h := Sum_b Res Phi_ydx omega_(h,1)` is used on the side of the
   curve; a parenthetical sentence records that it is EO's Def. 4.3 (4-10) up to the sign `F^TR_h = -F^(h)` and, as remarked after
   EO Thm 4.7, extends to `n = 0` the dilaton equation (4-30), stated there for `n >= 1`. The translated form `Sum_b Res Phi_ydx
   omega_(h,n+1) = (2 - 2h - n) omega_(h,n)` is right: EO print the factor `2g + k - 2` for `W`, and `W_n = (-1)^n omega_(h,n)`;
   independently, on the anchor `x = -p^2/2` of `conv:mod-kernel` it gives `-omega_(0,3)` from `omega_(0,4)` and `-omega_(1,1)` from
   `omega_(1,2)`, as it must.
3. **Inexact list of quoted inputs** (`s06_model.tex`, second paragraph; and the sentence after `thm:mod-tr`).
   (a) The DNOPS entry lacked `(2-30)-(2-35)`, the action of DNOPS s.2.3 with the dilaton sign of [PPZ15] Def. 2.13, the
   compatibility condition Prop. 4.3 (1), the sentence preceding Thm 8.1 and (8-1) at `a = (0,0,0,0,4v)`, all of which Step 1 of
   the proof of `thm:mod-tr` uses; the ABDKS entry lacked (7.6), (7.7), (7.9), (7.10) and the kernel (1.1). Now the same as the (Q)
   entry of `sec:int-inputs` (s01). The sentence after `thm:mod-tr` ("Its only quoted input is DOSS Thm 4.1 ...") now also names
   Prop. 4.3 (1) and the sentence before Thm 8.1.
   (b) The list said that EO Thm 4.7 (4-30) and EO (5-15) "enter here". No proof of the note uses them: `cor:mod-trF` uses the
   definition of `F^TR_h` (item 2), and `thm:mod-var` (ii), the only place where (5-15) enters, is used in no proof (`thm:mod-V` (ii)
   for `h >= 2` goes through `thm:mod-Lreg` and `thm:mod-gap`). The paragraph now lists the inputs that enter and says that these
   two items of (Q) are recalled in `cor:mod-trF` and `thm:mod-var` (ii) but are used in no proof (cross-section consequence:
   REVIEW_2 s.3, item 1).
4. **Inexact list of inputs of Theorem V** (`s06_model_9.tex`, sentence before the proof of `thm:mod-V`). It named EO Thm 5.1 and
   ABDKS Thm 7.1 only; the proof also rests on EO Thms 4.2 and 4.6 (poles only at the branch points, symmetry: the paragraph after
   `conv:mod-kernel`, `def:mod-pair`, `lem:mod-xfun`), and from EO Thm 5.1 only (5-14) enters. Fixed.
5. **Wrong item reference** (`s06_model.tex`, `def:mod-curve`): "`lem:mod-scale`~(b)"; the lemma has items (i), (ii) (first-level
   `enumerate` is labelled (i), (ii) by the preamble). Fixed to (ii).
6. **Reference to a later statement in a proof** (`s06_model_4.tex`, `prop:mod-LD`, Step 1): it cited `lem:mod-closure` "below" for
   the branch independence (no circularity, but OUTLINE s.0 requires that a proof cite earlier results). Replaced by a self-contained
   argument from `def:mod-legs`: every leg is a polynomial in `log p` over `Q[v][p, p^-1]`; the identity `d phi^b_(k+1) = -phi^b_k dx`
   is an identity of such expressions and survives `log p -> log p + 2 pi i N`; so the change `delta^b_k` satisfies `d delta^b_(k+1) =
   -delta^b_k dx` with `delta^b_0` constant, hence is a polynomial in `x`; such functions pair to zero by the integration by parts
   already in Step 1.
7. **Reference to a later statement in a proof** (`s06_model_8.tex`, proof of `thm:mod-var`): "the case (0,2) is proved directly in
   `lem:mod-genus0`". The case is covered by EO (5-14) (`g + k > 1`), so the sentence is a remark; moved after the proof ("also
   proved directly, in `lem:mod-genus0` (iii) below"), with "Part (ii) is used in no proof".
8. **Indirect reference to a later proof; missing justification** (`s06_model_3.tex`, proof of `prop:mod-Rtr` (ii)): "the graph
   formula of the theorem quoted in the proof of `thm:mod-tr`" (a proof that comes later). Now cites DOSS Thm 4.1 as recalled in
   DNOPS s.2.6 directly, and adds the justification of the dilaton sign transcribed from G4R R-1: in DOSS s.2.3 a dilaton leaf
   carrying `psi^k`, `k >= 2`, has the weight `[z^(k-1)](-R(-z))` on the unit, which sums to `psi (1 - R(-psi))` on the unit, and
   DOSS's `R(-z)` (which also decorates their ordinary leaves) is DNOPS's `R^-1(z)` (checked in both sources).
9. **Imprecise justification** (`s06_model_3.tex`, proof of `thm:mod-tr`, Step 1): "the compatibility of [DNOPS Prop. 4.3 (1)]
   holds, as `dy/df` ... has poles only at the branch points and no residue at `p = oo`". In the source, Prop. 4.3 gives `eta = 0`
   if (1) the curve is rational with `x(oo) = oo` in the global coordinate, or (2) `dy/dx` has poles only at the zeros of `dx`; the
   condition quoted was that of (2), and "no residue at oo" of a function has no meaning. Rewritten: (1) holds since the curve is
   rational and `f` has a pole at `p = oo`; moreover `dy/df = p^2/(4(p^5 - v))` has poles only at the zeros of `df`. (Both hold;
   (2) is what the proof of Prop. 4.3 uses.)
10. **Notation clash** (`s06_model_3.tex`, proof of `prop:mod-Rtr`, 15 occurrences): `Lambda_kappa := -5 kappa` clashes with the
    leg map `Lambda` of `def:mod-weights`; NOTATION_2 (renames: "TR's `Lambda = -5 kappa` -> the number `-5 kappa`") and GS_3
    ask for the number. Replaced by `sqrt(-5 kappa)` and `-kappa/(-5 kappa) = 1/5`.
11. **Wording** (`s06_model_3.tex`, Step 3 of the proof of `prop:mod-Rtr`): "By (Mod1), `F_b = eps_b J_(nu_b)(1; varsigma^4
    vartheta)`" attributes to (Mod1) what is the homogeneity `J_m(v; vartheta) = varsigma^(m-1) J_m(1; varsigma^4 vartheta)` of
    `lem:mod-laplace` (iii); fixed.
12. **Notation** (several files): `\Rhat^{\mathrm m}{}^{(j)}_b` -> `\Rhat^{\mathrm m,(j)}_b` as in `def:set-M0` (s02 and s03, s07);
    `\overline{M}_{h,n}` -> `\overline{\mathcal M}_{h,n}` as in s02, s04, s05; `\Omega^{\mathrm{mod}}_{v,h,0}` ->
    `\Omega^{\mathrm{mod}}_{v;h,0}` (`def:set-modcohft`); the local constants `kappa_b` of `prop:mod-LD` Step 5 (clash with the
    scaling constant `kappa` of s06 and the global `kappa_a`; OUTLINE_6 asks to rename them) -> `gamma_b`; the local letter `y` in
    the proof of `thm:mod-lattice` (c) (clash with `y = p`) -> `phi`; in the statement of `prop:mod-LD`, `q_(phi^b_k) = -varrho_b
    psi^k H^b` -> `-varrho_b z^k H^b` (`q_phi` is a vector of `H[z]`, as in its proof; the insertion `psi^k H^b` is identified with
    `z^k H^b` by `conv:tra-forms`, now said).
13. **Typo** (`s06_model_6.tex`, proof of `lem:mod-classL` (iv)): `(+)_(i >= 1) Q[v,x] p^i` -> `(+)_(i=1..4)` (basis `1, p, .., p^4`).
14. **Misparaphrase in a remark** (`s06_model_7.tex`, `rmk:mod-hock`): the residues of [Hock25] Cor. 1.2 were said to be taken
    "at the points where `y dx` or `x dy` has a non-zero residue". By Hock's Thm 1.1 they are taken at the logarithmic singularities
    of `x` and at the other singularities of `x` and `y`; for `x = W_v`, `y = p` these are `p = 0` and `p = oo`, so the conclusion is
    unchanged. Fixed; "This proof is not used" -> "This argument is only sketched here and is not used" (OUTLINE rule 10).
15. **References** (`s06_model_4.tex`, statement of `cor:mod-trF`): `F^m_h` was redefined; it is now recalled as the free energy of
    `def:set-modcohft`, and the potential `Fbar_h` is attributed to `def:set-modcohft`, `def:set-graphsum` (where it is defined; the
    same correction was made in s05 by its reviewer), not to `lem:set-dilaton`.
16. **References and typography** (`s06_model.tex`, `def:mod-pair`; `s06_model_4.tex`, `lem:mod-xfun`; `s06_model_5.tex`,
    `lem:mod-resid`; `s06_model_6.tex`, `cor:mod-split`). `def:mod-pair` used `lem:mod-xfun` (stated later) inside a justification:
    the holomorphy at `0`, `oo` is now attributed to the paragraph after `conv:mod-kernel` (EO Thm 4.2), and the zero residues to
    `lem:mod-xfun` "proved below; its proof does not use pairings" (true: only its "Consequences" paragraph mentions pairings, by
    their definition). In `lem:mod-xfun` the symmetry is attributed (EO Thm 4.6, paragraph after `conv:mod-kernel`). In
    `lem:mod-resid` the unnecessary citation of `lem:gau-vanmod` is dropped and two formulas are displayed (overfull lines of 38 pt
    and 54 pt removed); `cor:mod-split` is reworded (overfull line of 39 pt removed).

No label was added, removed or renamed; no statement changed its mathematical content.
