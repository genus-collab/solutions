# Review rf2: `sections/s02_setup*.tex` (part 1 of 2)

Author-facing; not part of the note. Object: the eight parts `s02_setup.tex` .. `s02_setup_8.tex`, checked line by line against
OUTLINE s.0 and s.0.2, OUTLINE_3 items 1-9, OUTLINE_10 items 10-17, NOTATION (+`_2`), LABELS, BIBPLAN, DECISIONS (Q1, Q3, Q7),
GAPS_SYNTHESIS s.2.1-s.2.4, s.2.10, GAPS_SYNTHESIS_2 s.3, the gap notes G1 (all pieces) with G1R and G1R_2, and the s03/s07-s10 statements
that s02 points to. Part 2 (`REVIEW_2.md`): what remains open, the changes that affect other sections, the checks on record.

## 0. Verdict: **sound-after-fixes**

* The mathematics of s02 is correct as transcribed from G1 (with the G1R corrections), DB, DD1/DD2 and G5, with one real error, now
  fixed: `def:set-outer` defined the outer edge by the formula of [PPZ15, Def. 2.2] and justified its being a power series by the
  symplecticity of the conifold-expanded rows, which the note had not proved (s.1, item S1).
* Four forward references inside proofs (to `def:set-frob` and `def:set-rows`) were removed by moving one computation and two proof
  paragraphs (S2).
* One new `\draftnote` (S11): the identification of the contraction `Cont_Gamma[E_R, q_R; Omega]` with the term of `Gamma` in the action
  of [PPZ15, Defs. 2.2, 2.5, 2.13] needs the comparison of the pulled-back classes `psi` (push-forward form of the translation) with the
  classes `psi` of the bigger moduli space (insertion form). The argument is two lines and is written inside the draftnote; neither the
  sources nor PPZ15 (Def. 2.5 and Remark 2.6 read in the v3 text) state it. It does not enter the proof of Theorem A as written in s05
  and s07 (there `T_Gamma` and `F_Gamma` are both contractions); it underlies `thm:set-GT` in integrated form, hence Corollary B, and
  possibly the identification in `thm:mod-tr` (s06). No input is added to the conditionality.
* Compiles with the current versions of all other sections: no error, no undefined reference, no overfull box in s02 (159 pages).

## 1. Issues, location, fix

**S1 (error; fixed). `def:set-outer` (`s02_setup_8.tex`), the symplecticity of `R'`.** The old text said "The matrix `R'` is symplectic
(`lem:flh-inputs` (iii)); so the numerator of `E_(R')` vanishes at `psi' = -psi`". `lem:flh-inputs` (iii) (old version) derived this from
`def:set-rows` (a), an identity of power series in `q`. The conifold expansion is defined on the presented ring `Nring (x) Q(zeta)`; to
carry the identity there one needs it in `Nring (x) Q(zeta)`, which is not proved (the injectivity of `ev` is not claimed, OUTLINE O4).
The fallback sentence ("written as `(N(psi,psi') - N(-psi',psi'))/(psi+psi')` the edge is a power series without this input") gave a
non-symmetric bivector, while `Cont_Gamma` needs a symmetric one. Decision, settled consistently with s07 and s10 (see `REVIEW_2.md` s.3):
* The edge is now *defined* without division problems and symmetrised:
  `E_(R') := (N(psi,psi') - N(-psi',psi')/2 - N(psi,-psi)/2)/(psi + psi')`, `N := 1 - (R'^-1(psi) (x) R'^-1(psi')) eta^-1`, `1` the inverse
  bivector of the colour space. It is a symmetric power series for every `R'` with `R'(0)` invertible; its `psi^a psi'^b` coefficient
  involves `R'^-1` only through `z^(a+b+1)`; it is the edge of [PPZ15, Def. 2.2] when `R'` is symplectic. This is the form that s10
  (`def:mai-Fcf`) now uses for `E_X`.
* A new paragraph proves the exact symplecticity defect, from `cor:gau-colours`, `lem:flh-inputs` (i), `prop:gau-symp`,
  `def:set-modcohft` (isometry `Psi^m`) and `lem:set-M0sympl`:
  `R'(-z)^T eta^-1 R'(z) = M(z) diag((D_i(z) D_i(-z))^-1) M(-z)^T` (equal to s07's `M(z)(D(z)D(-z))^-1 M(z)^-1`, as `M(-z)^T = M(z)^-1`).
  Hence `R'` is symplectic modulo `z^(3g-2)` under (F_g), and symplectic because `D = 1` (`thm:gam-D`).
* Non-circularity: the proof of `thm:gam-D` uses `thm:con-glue`, `thm:con-Pb`, `prop:gam-hatflat` (`thm:con-Pa` (iv), `thm:gam-main`
  (ii), `thm:gam-lmi`, appendix A) and `cor:gau-colours`; a grep of s03, s08, s09 and appendix A finds no reference to `def:set-outer`,
  `def:set-hat`, `lem:flh-*`, `R'` or `T_Gamma`.
* Where symplecticity of `R'` is actually needed: nowhere in the proof of Theorem A. `thm:flh-main` uses only `R' = R'' mod z^(3g-2)`
  (`lem:flh-inputs` (ii)) and the symplecticity of `R'' = P_can Psi^m`; with the new edge its Step 1 holds verbatim. It is needed only to
  read the outer graph sum as the action `R'.(M.TFT)` of [PPZ15, Def. 2.13], i.e. in the regrouping of Corollary B, where `cor:mai-gap`
  derives it from `D = 1`, as above.
* The constant term: "`R'(0)` is the isometry with lowered entries `eta(H^a, e~_i)`, which maps `1_TFT` to `1`" was a `q`-series
  statement. Now justified by the explicit form of `Psi_Q` in `Nring (x) Q(zeta)` (`lem:mai-Psi` (i)): `Psi_Q^T eta^-1 Psi_Q = 1` and
  `Sum_i Delta_i^(-1/2) (Psi_Q)_(ai) = (L/C_0)(Psi_Q)_(a0) Sum_i zeta^(i(a+2)) = 5 delta_(a3)`, identities in the ring, hence for the expansion.

**S2 (forward references in proofs; fixed).** (a) `lem:set-grading` cited `def:set-frob` (later) for `eta(H^(a+b), H^c) = 5 delta`;
`lem:set-unitfrob` and the proof of `cor:set-cyclicPF5` (a) used the same values. The localisation computation of `eta` and of
`eta^(cc')` is now in `def:set-fq`, where `eta` is defined; `def:set-frob` keeps the display with a pointer back. (b) The proof of
`prop:set-cyclicLR` (iii) used `def:set-rows` and `cor:set-cyclicPF5` (a), both later, and the last sentence of the proof of
`cor:set-cyclicPF5` (b) used `def:set-rows`. Both are now proved in a proof environment right after `def:set-rows`; the proof of (iii)
works with `w = R(z) e~_i` and the relation `z d_T w + mu_i w = H * w` of `def:set-rows` (b), so it needs no formal exponential; the WKB
statement keeps the formal exponential, with its derivation rule stated. The statements are unchanged (the proposition now says "defined
in `sec:set-rows`"). (c) The clause "large-radius analogue of `M_cyc` of `def:gau-conc`" moved from the proof to the paragraph after it.

**S3 (gap in an argument; fixed). `cor:set-cyclicPF5` (a).** "`Sum_i e_i = 1` gives `x_0 = 1/5`" tacitly assumed the same `x_0` for all
`i`. Now: the coefficient of `H^a` gives `Sum_i x_0^(i) zeta^(-ia) = delta_(a0)`, whose unique solution is `x_0^(i) = 1/5` (Fourier matrix).

**S4-S5 (minor; fixed).** `lem:set-periods` "Constants": the case `N = 0` (constants of `D` on `Q[[q]]`) is stated. Proof of
`prop:set-cyclicLR` (ii): commutativity of `*` and the constancy of `eta(H^a, .)` are named.

**S6 (notation; fixed). `lem:set-modhom` (a)** used `h_v` for the vertex genus, the source name that NOTATION renames to `g_v`.

**S7 (notation; fixed). `lem:set-dilaton` (i)** used `p` and `pi` for the same forgetful map; now `pi` throughout, and the position of
the unit insertion is fixed by symmetry.

**S8 (scope; fixed). `lem:set-dilaton` (iii)** was stated for complete local `Q`-algebras only, but `lem:flh-dil` applies it over
`A = Q(zeta)((s))` (a field) with `varkappa = d` in `1 + uQ[[u]]`. The hypothesis now reads: complete local with `varkappa` in `1 + m`, or
`Q(zeta)((s))` with the `s`-adic topology and `varkappa` in `1 + uQ[[u]]`; the proof is unchanged (`eps^b -> 0`). The remark after the
lemma names both rings.

**S9 (missing justification; fixed). `def:set-fq`** claims `F_g` in `Q[[Q]]`. Added one sentence: the localisation value is symmetric
in the weights (the permutation argument of (G0)(a), G1_2 s.2.1), so its value at `lambda_i = zeta^i` is Galois-invariant, hence rational.

**S10 (rule 4; fixed).** The Teleman `\draftnote` in `rmk:set-qrrderiv` Step 1 was a second copy; the first citation of [MR2917177],
with the draftnote, is in s01. The Givental v2 draftnote stays (first citation of [MR1408320] is in s02).

**S11 (new `\draftnote`; open). `def:set-graphsum`**, the sentence "for `g >= 2`, `Cont_Gamma[E_R, q_R; Omega]` is the term of `Gamma`
in `Int (R.Omega)_(g,0)`". [PPZ15, Def. 2.5] translates by push-forward along the forgetful maps (read in the v3 text, p.23: "(T Omega)
_(g,n)(v_1..v_n) = Sum_m (1/m!) (p_m)_* Omega_(g,n+m)(v_1..v_n, T(psi_(n+1)), ..)"; Remark 2.6, p.24), so the half-edges carry pulled-back
classes `psi`; `Cont_Gamma` inserts the points and uses the classes of `Mbar_(g_v, n_v+m)`. They agree: `psi_j - p_m^* psi_j` is the sum
of the divisors `D_(j,S)` (`S` a non-empty set of translation points on a rational component carrying only `j`, `S` and the node), and on
`D_(j,S)`, of rational factor of dimension `|S| - 1`, the classes `psi_l^(m_l)`, `m_l >= 2`, of the points of `S` vanish. Not in the
sources or in PPZ15; rule 1 forbids filling it, so it is a draftnote (with the argument). Used by: `thm:set-GT` ("`F_g` = graph sum"),
`def:mai-Fcf` (third property, `ev(F_g^Nring) = F_g`), `cor:mai-gap` Steps 1-2; to be checked in `thm:mod-tr` (s06) and `def:tra-edge`.

**S12 (scope; fixed). `def:set-graphsum`:** one sentence that the definitions and propositions of [PPZ15] used are identities between
finite sums and hold over every commutative `Q`-algebra (the note applies them over `Q(zeta)[[q]]`, `Q(zeta)[varsigma^(+-1)]`,
`Q(zeta)((s))`); it matches the parenthetical of `cor:mai-gap`.

**S13 (wording; fixed). `lem:set-canon`:** the canonical coordinates in the conifold chart are characterised only by `D u_i = zeta^i L`
(they are not identified with the large-radius ones by continuation); `cor:set-cyclicPF5` (a) now says why they are determined up to
constants at large radius (constants of `D` on `Q(zeta)[[q]][logq]`), pointing back to `lem:set-periods`.

**S14 (typography; fixed).** `rmk:set-cg`: a 40 pt overfull line; the formula for `c_h` is displayed.

**S15 (intro; fixed).** The opening paragraph of s02 now lists every later statement that the section points to (`lem:gau-confexp`,
`cor:gau-rowsNN`, `lem:gau-colrows`, the symplecticity of `R'` by `thm:gam-D`), and says that `def:set-confexp` uses the two operator
identities `lem:gau-adj`, `lem:gau-conc` (i) of s03, whose proofs use only the coefficients of PF4 (checked: no s02 input).

Continued in `REVIEW_2.md`.
