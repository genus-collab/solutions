# APPLIED_1: log of the application of AMEND_1 (2026-11-18)

Applied to `main.tex` and `sections/*.tex` with exact-string edits. Nothing compiled. `macros.tex`, `preamble.tex`,
`refs.bib` unchanged (no new macro, no new bibliography key). Snapshot before the edits: `archive_pre_amend_2026-11-18/`.
Decisions followed: A1 recommended variant with the lemma; A2 new part (iii); A3.4 version (a); B2 power-series
repair; B4 polynomial-identity argument; C3 replaced by A4-02 (AMEND_4_cross_lit.md gives a complete re-proof);
D9 not added.

Counts: 65 entries. APPLIED 62, APPLIED WITH ADAPTATION 2 (A1.16, C3.1), SKIPPED 1 (D9, by decision).

## A1 (flatness as part of (QRR); all-genus axioms in (T))

| item | status | file |
|---|---|---|
| A1.1 | APPLIED (whole item (QRR)) | s01_intro_3.tex |
| A1.2 | APPLIED (whole item (T)); "symplectic and" kept, see note 1 | s01_intro_3.tex |
| A1.3 | APPLIED (both replacements) | s01_intro_3.tex |
| A1.4 | APPLIED (both places) | s01_intro_3.tex, s10_main.tex |
| A1.5 | APPLIED | s01_intro_3.tex |
| A1.6 | APPLIED (both replacements, plus the sentence referring to `lem:set-Runique`) | s02_setup_4.tex |
| A1.7 | APPLIED (title, statement, sentence after the theorem) | s02_setup_4.tex |
| A1.8 | APPLIED | s02_setup_2.tex |
| A1.9 | APPLIED | s02_setup_3.tex |
| A1.10 | APPLIED | s02_setup_5.tex |
| A1.11 | APPLIED | s08_connection_8.tex |
| A1.12 | APPLIED | s09_gamma_3.tex |
| A1.13 | APPLIED | s10_main.tex |
| A1.14 | APPLIED | main.tex |
| A1.15 | APPLIED | b_signs.tex |
| A1.16 | APPLIED WITH ADAPTATION: the lemma and its proof are placed after the proof of `lem:set-hrot` (end of Section 2.4, before "The graph sum"), not after the proof at l.23-25, so that `\eqref{eq:set-hw}` and $h_{\bw_i}$ are backward references. Text unchanged. The sentence of A1.6 is then a forward reference inside `def:set-rows`. | s02_setup_4.tex |

Note 1. AMEND_1 says that with the lemma the words "symplectic and" "may be dropped" from A1.2 and A1.10. Not done:
A1.10 does not contain these words literally, and A1.6 (applied) says "(a) is part of the quoted input (T)". The text
now says both that (a) is part of (T) and that (a) follows from (QRR) by the lemma; this is redundant, not contradictory.
No statement "(T) is needed for Corollary B only" was added: AMEND_1 gives no text for it.

## A2 (edge condition)

| item | status | file |
|---|---|---|
| A2.1 | APPLIED (new item (iii) of `prop:mod-Rtr`, `eq:mod-edge`) | s06_model_2.tex |
| A2.2 | APPLIED (new paragraph "(iii)" at the end of the proof) | s06_model_3.tex |
| A2.3 | APPLIED | s06_model_3.tex |
| A2.4 | APPLIED | s06_model_3.tex |
| A2.5 | APPLIED (both replacements in (Q)) | s01_intro_3.tex |
| A2.6 | APPLIED | s06_model_2.tex |
| A2.7 | APPLIED | s06_model_2.tex |
| A2.8 | APPLIED | b_signs.tex |
| A2.9 | APPLIED (l.52 two replacements, l.58 one) | b_signs.tex |

The optional simplification after A2.4 (cut Step 2) was not made, as recommended there.

## A3 (formal change of frame)

| item | status | file |
|---|---|---|
| A3.1 | APPLIED (three replacements in the abstract) | main.tex |
| A3.2 | APPLIED | s02_setup_3.tex |
| A3.3 | APPLIED (new item 3 of "What is not claimed"; the later items move down by one; no cross-reference uses their numbers) | s01_intro_3.tex |
| A3.4 | APPLIED, version (a): new paragraph "Membership" before "Differences"; the rest of the old "Differences" line (the sentence on admissibility) is kept after "by Corollary B." | b_signs_2.tex |

## B

| item | status | file |
|---|---|---|
| B1.1 | APPLIED | a_formalsp_5.tex |
| B1.2 | APPLIED | a_formalsp_5.tex |
| B1.3 | APPLIED (both replacements in the proof) | a_formalsp_5.tex |
| B2.1 | APPLIED (both replacements) | s08_connection_4.tex |
| B2.2 | APPLIED (Step 3 replaced; `eq:con-bz` kept, `eq:con-factor` new). Also the optional rewriting in Step 1: `\chi A_{\mathrm{St}}/A_{\mathrm{St}}(L,0)` -> `\chi a_0` (one occurrence; the other one is in the line replaced by B2.1) | s08_connection_4.tex |
| B3.1 | APPLIED | a_formalsp.tex |
| B3.2 | APPLIED (statement and proof) | a_formalsp_2.tex |
| B3.3 | APPLIED | a_formalsp_3.tex |
| B3.4 | APPLIED | a_formalsp_3.tex |
| B3.5 | APPLIED | a_formalsp_3.tex |
| B3.6 | APPLIED (three places) | a_formalsp_4.tex (two), s08_connection_5.tex |
| B3.7 | APPLIED (statement (iv), end of the proof of (iv), proof of `lem:fsp-fubini`) | a_formalsp.tex, a_formalsp_4.tex |
| B4.1 | APPLIED | s02_setup_4.tex |
| B4.2 | APPLIED | s02_setup_5.tex |
| B4.3 | APPLIED | s10_main_2.tex |
| B4.4 | APPLIED (two deletions and the last sentence; done in one edit with D8) | s01_intro_3.tex |
| B4.5 | APPLIED (one edit with D8) | s10_main_2.tex |
| B5.1 | APPLIED | s10_main.tex |
| B5.2 | APPLIED (optional; "apply everything") | s03_gauge_8.tex |

## C

| item | status | file |
|---|---|---|
| C1 | APPLIED | s02_setup_4.tex |
| C2.1 | APPLIED | s03_gauge.tex |
| C2.2 | APPLIED (statement l.20 and proof l.31) | s03_gauge_2.tex |
| C3.1 | APPLIED WITH ADAPTATION: A4-02 applied INSTEAD of C3.1 (decision 5: A4-02 gives a complete re-proof of (ii) from `prop:mod-period`). Differences from C3.1: the proof is split into "(i) ..." and a separate paragraph "(ii) This part is not used. ..."; the sentence "on $\bbP^1$ there is no period matrix and no parameter $\kappa$, so the covariant variation $D_\Omega$ ... is $\delta_\Omega$" of A4-02 is included; the last sentence of the old proof ("The hypotheses hold: ...") is unchanged, as A4-02 asks. | s06_model_8.tex |
| C3.2 | APPLIED (identical in C3.2 and A4-02) | b_signs_2.tex |
| C4 | APPLIED; draft note deleted | s10_main_2.tex |
| C5 | APPLIED; draft note deleted | s10_main_2.tex |
| C6 | APPLIED | s01_intro_2.tex |

## D

| item | status | file |
|---|---|---|
| D1 | APPLIED (two places) | s02_setup_8.tex, s07_comparison_2.tex |
| D2 | APPLIED | b_signs.tex |
| D3 | APPLIED | s02_setup_8.tex |
| D4 | APPLIED | s07_comparison_2.tex |
| D5 | APPLIED (three places) | s10_main.tex, s01_intro_2.tex, s02_setup_8.tex |
| D6 | APPLIED | s08_connection_5.tex |
| D7 | APPLIED (Step 3 and Step 4) | s08_connection_5.tex |
| D8 | APPLIED (both lists) | s10_main_2.tex, s01_intro_3.tex |
| D9 | SKIPPED (decision 6: not to be added) | - |

"Inputs after the amendments" (AMEND_1_chain_4, last section): contains no edit of its own; it restates A1, A2, B4.
The texts of (Q), (QRR), (T) in `s01_intro_3.tex` now agree with it.

## Consequential edits (not listed in AMEND_1)

1. `s06_model_3.tex`, proof of `thm:mod-tr`, Step 1: `as recalled in \cite[\S2.6, (2-30)--(2-35)]{...}` ->
   `(2-30)--(2-32), (2-34), (2-35)`, to agree with the new (Q) (A2.5), from which (2-33) is removed. Minimal change.

Looked at and left unchanged:
* `s09_gamma_3.tex` l.82, last sentence ("only through \cref{def:set-rows}~(b), $R(0)=1$ and (QRR)"): (b) is now
  part of (QRR), so the sentence is redundant but true. A1.5 rewrote the parallel sentence of `s01_intro_3.tex`; no
  text is given for this one.
* `s01_intro.tex` l.16, `s01_intro_2.tex` l.36, l.46 ((QRR) named for the value only): AMEND_1 says they stay.
* `s02_setup.tex` l.14 and the first clause of `s01_intro_3.tex` l.17 ("(M) and (G0) enter only in ..."): left for
  A4-05 (README, "Duplicates").
* `s10_main.tex` l.31 "\cref{prop:mod-Rtr}~(ii)" and `s06_model_3.tex` "unit ... (\cref{prop:mod-Rtr}~(ii))": still correct.
* `s01_intro_3.tex` l.23 (list of citations that enter no proof, with (5-15)): now true, as A4-02 says.

## New labels

`lem:set-Runique` (s02_setup_4.tex), `eq:mod-edge` (s06_model_2.tex), `eq:con-factor` (s08_connection_4.tex).
Checked by script: no duplicate label in `main.tex` and `sections/`; every `\cref`, `\Cref`, `\eqref` target exists;
every citation key is in `refs.bib`; counts of braces, `$`, `\begin`/`\end`, `\[`/`\]`, `\bigl`/`\bigr` and parentheses
are balanced in each modified file exactly as before. Commands used in the new text that are not in `macros.tex`:
`\Bbbk` (amssymb, already used in appendix D), `\mathscr` (eucal), `\varrho`, `\mathsf`; all available.

## Files modified (28)

main.tex; sections: a_formalsp, a_formalsp_2, a_formalsp_3, a_formalsp_4, a_formalsp_5, b_signs, b_signs_2,
s01_intro_2, s01_intro_3, s02_setup_2, s02_setup_3, s02_setup_4, s02_setup_5, s02_setup_8, s03_gauge, s03_gauge_2,
s03_gauge_8, s06_model_2, s06_model_3, s06_model_8, s07_comparison_2, s08_connection_4, s08_connection_5,
s08_connection_8, s09_gamma_3, s10_main, s10_main_2 (.tex).

## Mathematics: what I rechecked while typing, and doubts

Rechecked by hand, no error found:
* A1.16: the recursion $(\Dq Y_k+\Theta Y_k)_{ji}=(\zeta^j-\zeta^i)L(Y_{k+1})_{ji}$ and $z\Dq N=[\Dq U,N]$.
* A2.2: the formula for $\mathsf D\omega_{0,2}$, the limit $f'''/6$, the cancellation of the simple poles, the principal
  part, the identity for $\mathsf D(\dd w\,\dd w'/(w-w')^2)$, and the signs in Step C.
* A3.4 (a): $\mathsf W=c_0^2+c_0\pi_1$, $\pi_2=(\Dq\mathsf W-2c_0c_1)/c_0$, the action of the shift on $x$, $\pi_j$, and
  the weights $\zeta^{1-n}$, $\zeta^{a+1-n}$, $\Delta_i^{1/2}\mapsto\Delta_{i+1}^{1/2}$.
* B1.3: the equation for $\mathsf T_{\mathsf t}(E)$ and the second solution.
* B2.2: \eqref{eq:con-factor}, the decomposition of $\mathsf S_0$, $\mathsf S_1$, $\mathsf S_{2k}$ into the six factors, the
  extra term $(\mathsf T_z-1)\mathsf l_{\rho_{\mathrm N}}$, and the product of the multipliers.
* B3.3 (the common exponent), C1 (the equation and $u^{4/5}$), A4-02 (ii) (against `prop:mod-period` (ii), (iii)).

Doubts and points for the author (nothing was changed on their account):
1. B4.2 against its use. The new sentence in `lem:set-pushins` (iii) says that over a ring other than $\bbQ$
   [PPZ15, Prop. 2.12] "is used only for the topological field theories of `def:set-graphsum`". In the proof of
   Corollary B (`s10_main_2.tex` l.27) (iii) is applied to $M.\Cexp(\TFT)$ and the series $X$, and $M.\Cexp(\TFT)$ is
   not a topological field theory. The density argument of B4.1 covers this ($X.(M.\TFT)=(XM).\TFT$ by Prop. 2.14,
   and all conclusions are polynomial identities in $y$, $r$), but the wording "only for the topological field
   theories" should perhaps read "only for theories $R.\TFT_y$ obtained from the topological field theories of ...".
2. A1.2, A1.6, A1.10 with A1.16: symplecticity is at the same time quoted as part of (T) and proved from (QRR) in
   `lem:set-Runique` (note 1 above). The author may want to drop it from (T).
3. A1.16: the proof uses silently that the $\tilde e_i$ form a basis over $\bbQ(\zeta)\fps q$ (so that $\Theta$ and $Y$
   have entries there) and that $L$ is a unit of $\bbQ(\zeta)\fps q$; both hold by `cor:set-cyclicPF5` (a), but are
   not said. Existence ("exactly one") is by solving the recursion, which the proof shows only implicitly
   ("the solutions correspond to the families of constant terms").
4. A1.10: the all-genus divisor equation is declared a quoted input, yet the parenthesis "(as in genus 0, ...:
   $E_{g,n+1,d}$ is the pull-back of $E_{g,n,d}$ and $\int_dH=d$)" reads as a proof sketch. The notation $E_{g,n,d}$
   occurs in the note only in genus 0 (`thm:set-G0`), so it is used here without definition for $g>0$. In the
   intro item (T) the coefficient ring $\bbQ(\zeta)\fps{\Qn}$ of A1.10 is not mentioned.
5. A2.2, Step A: "holomorphic at $p_1=0$ ... ($f$ has a double zero at $0$)". The form is holomorphic at $p_1=0$
   simply because $f$ is regular there; the double zero is not needed. Harmless.
6. A2.6: "(2-32) is $R^{-1}$ applied to the unit [DNOPS, Section 3]" and A2.1 "[DOSS, (4-1), with (3-2) and
   (3-57)]", A2.7 "[DOSS, (2-23)]", A2.8 "[DOSS, (2-26)]": pins typed as given, not checked by me.
7. A3.1: the abstract now uses the phrase "the holomorphic limit in the conifold frame", which is not defined in
   the note. "No analytic continuation is involved" is the cautious wording asked for in AMEND_1.
8. A3.4 (a): the step "$\varrho$ permutes the colours and fixes $F_g^{\Nring}$" uses that $\varrho$ fixes $\zeta$ and that
   the whole graph sum (edges, leaves, vertices) is built from the rows, $\Delta_i^{1/2}$ and $\Psi_{\mathrm Q}$ only; it
   is the same argument as for $\sigma_k$ in `def:mai-Fcf`, and I did not find a gap, but it is new and uses
   Theorem C. The definitions of [CGYZ26] are those recalled at l.100 of `b_signs_2.tex` (from an HTML rendering).
9. D6: the letter $a_0$ is used in the proof of (ii) but is "local to the proof of (i)" (B2.1), and the list of
   local letters of the proof of (ii) ($\phi$, $J$, $a$) was not extended.
10. D7, Step 4: "the sums over $j$, which lie in $\bbQ[s,X,(1-s^5X^5)^{-1}]$": not checked.
11. D1: "the matrix $\eta^{-1}=\eta/25$": consistent with $\eta$ being 5 times a permutation matrix; not rechecked
    against `def:set-frob`.
12. C2.2: "for the rows it is obtained in `thm:gau-R0` (iv) and `cor:gau-colours`": typed as given; `thm:gau-R0` (iv)
    does mention the closing relation; `cor:gau-colours` not re-read.
13. C4, C5: the computations remain unrecorded; the wording is only softened, as AMEND_1 says.
