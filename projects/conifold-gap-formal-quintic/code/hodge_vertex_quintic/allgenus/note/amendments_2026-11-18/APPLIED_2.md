# APPLIED_2: log of the application of AMEND_2 (2026-11-18)

Applied to `sections/*.tex` with exact-string edits, after AMEND_1 (`APPLIED_1.md`). Nothing compiled. `main.tex`,
`macros.tex`, `preamble.tex`, `refs.bib` unchanged (no new macro, label or bibliography key). All AI work.
Every old fragment of AMEND_2 was found verbatim in the sources as they stood after AMEND_1 (line numbers shifted in
`s02_setup_4.tex` only), except A2-38, which another pass had already changed.

Counts (items A2-01..A2-39): APPLIED 32, ADAPTED 2 (A2-06 optional part, A2-39), ALREADY DONE 4 (A2-01, A2-02, A2-03,
A2-38), SKIPPED 1 (A2-11, reserved) plus one sub-item of A2-39 (comment in `macros.tex`). Rejected R-1..R-7: not applied.

## P1

| item | status | file |
|---|---|---|
| A2-01 | ALREADY DONE (AMEND_1 C1; not touched) | s02_setup_4.tex |
| A2-02 | ALREADY DONE (AMEND_1 C2.1; not touched) | s03_gauge.tex |
| A2-03 | ALREADY DONE (AMEND_1 C2.2; not touched) | s03_gauge_2.tex |
| A2-04 | APPLIED (three places: proof after `def:set-rows`; Step 1 of `lem:gau-colrows`; proof of `thm:gau-R0` (iii)) | s02_setup_4.tex, s03_gauge_7.tex, s03_gauge_8.tex |
| A2-05 | APPLIED | a2_cyclicPF4.tex |
| A2-06 | APPLIED (statement (ii)). Optional part ADAPTED: the proof steps (1)-(5) are renamed `\emph{Step 1.}` ... `\emph{Step 5.}`, and "using (2)" became "using Step~2" (the item gives the names of the steps, not the sentence). No other file refers to these step numbers (grep). | s03_gauge_6.tex |
| A2-07 | APPLIED (both places) | s04_transport.tex, s04_transport_2.tex |
| A2-08 | APPLIED | s05_reduction.tex |
| A2-09 | APPLIED (both places) | a3_uconv.tex |
| A2-10 | APPLIED | s05_reduction.tex |
| A2-11 | SKIPPED (README: reserved for A4-05, a later pass) | s02_setup.tex |

## P2

| item | status | file |
|---|---|---|
| A2-12 | APPLIED | s04_transport.tex |
| A2-13 | APPLIED (both places) | s04_transport_2.tex |
| A2-14 | APPLIED (the first of the two occurrences of `(\cref{lem:tra-frames}~(c))`, the one after "$(4(k-k')+1)/5$"; the second, in "(i)$\Rightarrow$(ii)", is a correct reference and is unchanged) | s04_transport_2.tex |
| A2-15 | APPLIED | a3_uconv.tex |
| A2-16 | APPLIED (both places) | a3_uconv.tex, s05_reduction.tex |
| A2-17 | APPLIED | s03_gauge_6.tex |
| A2-18 | APPLIED | s02_setup_7.tex |
| A2-19 | APPLIED | s02_setup_8.tex |
| A2-20 | APPLIED | s02_setup_4.tex |

## P3

| item | status | file |
|---|---|---|
| A2-21 | APPLIED | s02_setup_6.tex |
| A2-22 | APPLIED (three places) | s02_setup_2.tex, s02_setup.tex, a1_cyclicLR.tex |
| A2-23 | APPLIED (inside the parenthesis) | s02_setup.tex |
| A2-24 | APPLIED (both places) | s02_setup_3.tex, s02_setup_4.tex |
| A2-25 | APPLIED (pairing in two files; range of the sum; symmetry of $E$) | s02_setup_4.tex, s02_setup_5.tex |
| A2-26 | APPLIED (both places) | s03_gauge_2.tex |
| A2-27 | APPLIED (both places) | s03_gauge_4.tex |
| A2-28 | APPLIED | s03_gauge_5.tex |
| A2-29 | APPLIED | s03_gauge_5.tex |
| A2-30 | APPLIED | s03_gauge_7.tex |
| A2-31 | APPLIED (three places) | s03_gauge_8.tex |
| A2-32 | APPLIED | s04_transport.tex |
| A2-33 | APPLIED | s04_transport_2.tex |
| A2-34 | APPLIED | s04_transport_4.tex |
| A2-35 | APPLIED (both places) | s04_transport_4.tex, b_signs.tex |
| A2-36 | APPLIED (both replacements) | a3_uconv.tex |
| A2-37 | APPLIED (`thm:red-main` (b) only; `cor:red-can` unchanged, as the item says) | s05_reduction.tex |
| A2-38 | ALREADY DONE: when I reached the line it read "(\cref{lem:set-dilaton}~(iii); for the flat model correlators this is \cref{lem:tra-conv}~(f))", changed by another pass during this session; same content as A2-38. Not touched. | s07_comparison.tex |
| A2-39 | ADAPTED, see below | a3_uconv.tex, s02_setup_7.tex |

A2-39, details.
* `a3_uconv.tex`: `D := 3h-3+n` -> `\mathsf D_{h,n} := 3h-3+n`, and every later `D` with this meaning in the statement
  and in the proofs of (b1)-(b3) -> `\mathsf D_{h,n}` (11 occurrences in all); in the proof of (b4)
  `D_w := 3g_w-3+n_w`, `D_{w'}` -> `\mathsf D_w`, `\mathsf D_{w'}`. The explicit expressions "$3h-3+n$" in (d), (f) are unchanged.
* `s02_setup_7.tex`, proof of `lem:set-modhom` (b): the item changes only "The degree is $D=$". The same sentence
  continues "$2h-2+\sum_j\nu_{b_j}-4D=5e$" with the same letter; I changed it to `-4\mathsf D` as well, so that the
  sentence does not use two letters for one number. The statement (a) of the lemma and the proof of (a) keep the
  letter $D$ ("for every cohomological degree $D$"), as the item does not mention them; (b) now applies (a) with
  $D=\mathsf D$.
* `macros.tex` l.61, comment of `\cR`: SKIPPED, the item says "Comment only" and gives no new text.

## New labels

None. Checked by script after the edits: no duplicate label in `main.tex` and `sections/`; every `\cref`, `\Cref`,
`\eqref` target of the 22 modified files exists; in each modified file the counts of braces, `$`, parentheses,
brackets, `\bigl`/`\bigr`, `\begin`/`\end` are balanced exactly as in `archive_pre_amend_2026-11-18/`, each line has
balanced braces and `$`, and the number of lines is unchanged (one paragraph per line kept; `s02_setup_4.tex` differs
from the archive only by the lines added by AMEND_1). Commands in the new text that are not in `macros.tex`:
`\mathsf`, `\mathscr`, `\Bbbk`, `\mathrm{GL}`, `\operatorname`; all already used in the note.

## Files modified (22, all in `sections/`)

s02_setup, s02_setup_2, s02_setup_3, s02_setup_4, s02_setup_5, s02_setup_6, s02_setup_7, s02_setup_8, s03_gauge_2,
s03_gauge_4, s03_gauge_5, s03_gauge_6, s03_gauge_7, s03_gauge_8, s04_transport, s04_transport_2, s04_transport_4,
s05_reduction, a1_cyclicLR, a2_cyclicPF4, a3_uconv, b_signs (.tex).

## Mathematics: what I rechecked while typing, and doubts

Rechecked by hand, no error found:
* A2-08: $\sum_w(2g_w-2+n_w)=2g-2$ for a connected leg-free graph of genus $g$.
* A2-10: $2h-2+n>0$ gives $n\geq3$, $n\geq1$, all $n$ for $h=0$, $h=1$, $h\geq2$.
* A2-30: with $L=s^{-1}$ the stated $\Dq$ gives $\Dq L=(L^6-L)/5$.
* A2-04: $A^{\mathrm{act}}=C_0\Qr^{(0)}_0=L\Phi_0$ by the definition of $\Phi_0$ in `def:set-rows`.
* A2-06: the coefficient of $\zh^0$ is preserved by $\mathscr O_e$ (formula of Step 1 at $m=0$).
* A2-28: $d^U_0:=1$ in `s03_gauge_2.tex` l.60, consistent with $\cn(0)=p_{0,0}(0)$.

Doubts and points for the author (nothing was changed on their account):
1. A2-04 (`s02_setup_4.tex`). The new text calls $\mathsf B$ "a differential $\bbQ\fps q$-algebra in which $z$ is a
   constant"; the coefficients are in $\bbQ(\zeta)\fps q$ and $\Dq$ is not $\bbQ\fps q$-linear, so "differential
   algebra over the differential ring $\bbQ\fps q$" is meant. The step "$C_0v_0=\zeta^iL\Phi_i\,\mathsf E$ solves
   $\PFfive$ by `cor:set-cyclicPF5` (b)" applies that corollary over $\mathsf B$; I did not check that its statement
   allows this coefficient algebra. In `s03_gauge_7.tex` and `s03_gauge_8.tex` the new words "the proof after
   `def:set-rows`" point to an unnumbered proof.
2. A2-05. "The right side has degree at most $1$ in the indeterminate $\ell$ over $\bbK$" uses that $\ell$ is
   transcendental over $\bbK$ (or that the appendix treats it as an indeterminate); not checked. The conclusion is
   stated "for $f\in\mathcal S\cap\bbK$", while the sentence before it has no such restriction on $f$.
3. A2-06. The new clause "$1-u$ being non-zero in the field $\bbQ\lser s$" fixes the ring to $\bbQ\lser s$; A2-20 says that
   (i), (ii) "hold in the same way for coefficients in $\bbQ(\zeta)\fps q$", where $1-u=5^5q$ is not a unit but is a
   non-zero-divisor, which is all that is needed. Consistent, but said in two places in two ways.
4. A2-09. In (e2) the member $(l,j)=(0,0)$ is now excluded from "expanded in frames by (b2)"; the members with
   $(l,j)\neq(0,0)$ are correlators with $n=l+j\geq1$ slots, as (b) requires. In (b4) the new reference
   `lem:red-genus` (Section 5) is cited from the statement of a lemma used in Section 4 and 5; the label exists, and
   I did not check that this lemma states $n_w\geq1$ for graphs with an edge.
5. A2-16. The new list of dependencies of `lem:tra-conv` and the claim "which precede `sec:tra-dress` and do not use
   this lemma" are typed as given (AMEND_2 says it was checked by grep); not rechecked.
6. A2-19. "$X\coloneqq R'(0)^{-1}R'$ (symplectic, $X(0)=1$)" and the reference `lem:set-pushins` (iii): typed as given.
   The sentence is inside the clause "If $R'$ is symplectic, ...", which is where $X$ is symplectic.
7. A2-21. The sign pattern $\sum_a(-1)^{a+1}f_ag_{3-a}$ and "the coefficients of $\PFfour$ below the leading one are
   multiples of $1-u$" were not recomputed; `lem:gau-conc` (iii) was not re-read.
8. A2-22. "(G0)(b)" as the name of the divisor equation, the ring $\bbQ[T]\fps{\re^T}$, and "see (d) below" (the
   item (d) of the paragraph "Translation" of `s02_setup.tex`; the same file also has an `\item[(d)]` in (G0), so the
   pointer is not unique) are typed as given.
9. A2-23. The list of localisation weights is typed as given; AMEND_2 itself says it was not checked against a
   source. The classes "$(\lambda_i-\lambda_j)/d-\psi$" are invertible for a different reason ($\psi$ nilpotent,
   $\lambda_i\neq\lambda_j$) than the convex-combination argument quoted for the weights.
10. A2-25. The pairing is now "perfect symmetric" in `def:set-graphsum` and `lem:set-dilaton`; the sentence added
    by AMEND_1 B4.1 in the same definition says "with a perfect pairing". Consistent. The restriction of the sum to
    $2h-2+n>0$ in $\overline F^\Omega_h$ was not propagated to the parallel formula for $\dd^n\overline F^\Omega_h$ in the
    same sentence (the item does not ask for it).
11. A2-34. The four citations added to example (i) of `def:tra-dressing` are typed as given; not checked.
12. A2-39. See the adaptation above; the letter `\mathsf D` is also used by AMEND_1 A2.2 for an operator in
    `s06_model_3.tex` (a different file, local there).
