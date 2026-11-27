# Review rf7, part 2: appendix A lemma by lemma; edits; changes touching other sections; compilation

Continuation of `REVIEW.md` (s.0 verdict, s.1 LMI step 2, s.2 s09 statement by statement).

## 3. Appendix A against FSL (with RL-C1..C8), CP2 and the other sources

Every proof was re-done by hand. "RL-Cn" are the corrections of RL2 s.5; all eight are applied where they apply.

| label | source of record | result |
|---|---|---|
| `def:fsp-ring`, `lem:fsp-ring` (i)-(v) | FSL s.1.2-s.1.3, Lemma 1.1, (W), (G2); RL-C5 (outer variables) | correct; (iv) the finiteness count for `E^t` (terms with `beta <= delta`, `j + |alpha|/2 = k`, weight `<= k + |delta|/2`) checked; (v) holds for any `J` (homogeneous substitution) |
| `def:fsp-wick` | FSL Def. 1.2-1.3; RL-C7; REV_2 M11 (b) (linearity, input of LMI) | correct; constant term `a_0(c)` from (eq:fsp-univ) at `k = 0` |
| `rmk:fsp-analytic` | FSL s.1.4 "Meaning"; no Hoermander (rule 4) | correct, not used |
| `lem:fsp-univ` | FSL Lemma U with CP (F1)'s sharp bound | correct: `r <= 2(k-m)`, phase orders `<= 2(k-m)+2`, amplitude orders `<= 2(k-m)` re-derived |
| `lem:fsp-const`, `lem:fsp-homog` | CP (F2); RL s.3.2 Lemma H (RL-C1) | correct (moment homogeneity of degree `|alpha|/2` in `C`; symbols) |
| `lem:fsp-ibp` | FSL2 s.3.1, (G3) | correct (the equivalence with `E[t_k g] = z Sum_j C_kj E[d_j g]` via `C H = -1`) |
| `lem:fsp-pg` (i)-(iii), second proof | CP2 s.2.2 (algebraic, relied on); FSL (G5), (G5') with RL-C4 (wording of the ideal), RL-C5 (`P` in `R_(n+m)(A)`) | correct: `dC_s/ds = C_s K' C_s`, the double Wick step, `d log d_s/ds = tr(C_s K')/2`, the `I^k` bound and the evaluation at `s = 1`; the integrability condition `||k|| ||C_0|| < 1` of the second proof is right |
| `lem:fsp-prod` | FSL (G4); CP2 s.2.3 | correct |
| `lem:fsp-ift` | FSL2 s.3.2 with RL-C4 ("a matrix with entries in `I`, times `(c - c')`") | correct (difference quotients, contraction, closedness of `I`, invertibility mod `I`) |
| `lem:fsp-transl` (a), (b) | FSL2 s.3.2 Lemma T | correct (polynomiality in `s`, `dB/ds` formula re-derived) |
| `cor:fsp-movebase` | FSL2 Corollary T'; CP2 s.2.6 (iii) | correct (cancellation of `F_(0,>=3)`, (iii) of `lem:fsp-pg` over `A_0[[q]]`, branch rule) |
| `lem:fsp-param` (i)-(iv) | FSL2 s.4 (i)-(iv), RL-C4 in (iv); DEC Q2 | correct; (i) asserts holomorphy only ALONG a given holomorphic family (no holomorphic IFT), as DEC Q2 requires; (iv) the weighted ring with `eps` of weight 1 and the substitution `eps -> z` checked |
| `lem:fsp-dq` | FSL2 s.4.5 with RL-C2 (`D a = 0` for the iterated clause) | correct |
| `lem:fsp-generic` (i)-(iii) | G5 s.1.3 Steps 1-2 [G] | correct: `AL = Q[L, delta^-1]`, `r = 5 sigma/(sigma^5 - 1)` a unit with `delta r = L^2` (sign irrelevant, only ratios enter), the ODE uniqueness, the two-homomorphism arguments, the derivative `(L/z)N + (DL/L)N + DN = N^sigma/z`, and the necessity example `z (D f) L SPn[..]` |
| `lem:fsp-piola` | FSL2 s.5 (Piola) | correct (signs `(-1)^(k-1)`, `(-1)^(i+1)`) |
| `lem:fsp-cov` | FSL2 s.5 L1a with RL-C3 (`P_s` with `- F(0)`); CP2 s.2.4 alternative | correct: (a), (b), Jacobi's formula, Piola, `lem:fsp-ibp` over `A[s]` |
| `lem:fsp-cfu` | FSL2 s.5.4 | correct (`Y'(X_0) = psi(1;u) + 5 d_T psi(1;u)`) |
| `lem:fsp-fubini` (i)-(iii) | FSL3 s.6 with "Rings" (RL s.2) and RL-C3, RL-C5; CP2 s.2.5 | correct; the counterexample `Sum_k v_1^(2k) w_1^(2k) z^-k` checked |
| `def:fsp-stirling` | DD2 D2; FSL4 s.9.1 | correct (`beta_(2k-1)(0) = bk_k`, `beta_(2k)(0) = 0`, `Poch_(M+1) = (1 + Mw) Poch_M` for all `M` in Z) |
| `lem:fsp-stirfe` | FSL4 s.9.2 and the bookkeeping of s.9.1 | correct; every identity of the `D_Y^2` computation re-derived |
| `lem:fsp-stirling` (i)-(iv) | FSL4 s.9.3, Proof A (formal, the proof the chain uses, ST5b K2); FSL's (v) is the note's (iv) | correct: the homogeneity count `a^(r - Sum l_i + k + r) = a^(-k)`, Step 1 (`lem:fsp-ibp` with `h = x^m`), Step 2 (`lem:fsp-param` (iv), `F_eps`, Hessian ratio `a/(a + eps)`, the weighted-ring factors), Step 3 (injectivity of `f(Y) - f(Y+1)`), Step 4 |
| `rmk:fsp-classical` | CP2 s.2.8; FSL4 s.9.4; RF-C6 | correct, not used; DLMF 2.3.8 is Watson's lemma (CP2 s.2.8 step 4 quotes it), 5.11.1 Stirling, 5.11.8 the shifted form |
| `lem:fsp-dlmf5118` | LF s.4.2 = LA s.2.3 | correct; both sides in `w^-1 Q[theta][[1/w]]`, the shift identity via `rho(w + theta) = 0` and `1 + (theta+1)/w = (1 + theta/w)(1 + 1/(w + theta))`, `B_k(theta+1) - B_k(theta) = k theta^(k-1)` from 24.2.3 |
| `lem:fsp-refl` | RL s.3.3; branch part not formal (REV_2 M11 (a)) | correct |
| `lem:fsp-res` (i)-(iii) | FSL4 s.10.1-s.10.3 (PS-C2) | correct; the relabelling is `lem:con-relab` (s08), as stated |

The first paragraph names both non-formal places (the exact Gaussian integral of the second proof of `lem:fsp-pg`, the branch part of
`lem:fsp-refl`), as REV_2 M11 (a) asks. `lem:fsp-delta` and `lem:fsp-L3C` are not written, as LABELS_2 prescribes.

## 4. Edits made (all in my files; no label added, removed or renamed)

* **E1** `s09_gamma.tex`, header comment: local letters `chi, mu, Phi_2, Psi_Q, p(L)` declared; discs are `\bbD`.
* **E2** `def:gam-lmidata`: "the Stirling form of `Gamma(M_* + a/z)` without `(2 pi z)^(1/2)` and its phase" replaced by the exact object:
  the expansion of `z^(M_* + a/z) Gamma_St(M_* + a/z)` (`lem:fsp-stirling` (iii), `lem:fsp-stirfe`) without `(2 pi z)^(1/2)` and without the
  phase factor `e^((a log a - a)/z)`. (Without the factor `z^(M_* + a/z)` the expansion carries `sigma`-dependent `log z` terms; LA-C1.)
* **E3** The disc `D` of `def:gam-lmidata`, of the proof of `thm:gam-lmi` (steps 0, 2, holomorphy), of `thm:gam-main` (ii) and of
  `thm:gam-ginf` renamed `\bbD` (as `\bbD_L` in `def:gam-c`): `D` is the operator `q d/dq` (`\Dq`, NOTATION s.11; rule 6).
* **E4** Step 2: `Psi_Q` (and `Psi_0`) declared local to the proof; `Psi_0` is also the matrix of `def:set-M0`.
* **E5** Step 1: the factor `FF` is carried along by the `C[[sigma - L]][[z]]`-linearity of the inner `SPn` (the text said
  `C[[z]]`-linearity, but `FF(sigma, z)` depends on `sigma`); the outer amplitude is now displayed (it overflowed by 24pt).
* **E6** Step 2: the claim "the Taylor expansion of `p` at `s_*` is one" justified (Taylor expansion commutes with composition), and the
  sentence added that this uniqueness is not needed: the continuation uses only the explicit family `p`, so no implicit function
  theorem (formal or holomorphic) and no convergence of a family in `Q` enters (REVIEW s.1).
* **E7** Step 3: the constant `phi(s_*)` of the phase removed by `lem:fsp-const` before `lem:gam-deltasrc` is applied (the lemma's phase
  is exactly `P(s) + a_* s m`).
* **E8** The input list after the proof: `lem:fsp-stirfe` added (form of `A`); where `lem:fsp-ift` enters (inside `lem:fsp-fubini`; the
  uniqueness sentence); the linearity ring; "no analytic theorem" made precise (elementary holomorphy only).
* **E9** `thm:gam-main` (iv): "in the notation there" after `m = 5, k = 1` (the `k` of [GJR18], not the Bernoulli index).
* **E10** Proof of `thm:gam-main`: "their proofs use (i), (iii) and `thm:gam-lmi` only" was inaccurate (`thm:gam-ginf` uses
  `lem:fsp-dlmf5118`; `cor:gam-poly` uses `thm:con-Pa` (iv) for its further statements); corrected. No circularity either way.
* **E11** (merged into E3: the disc of `thm:gam-main` (ii).)
* **E12** `thm:gam-D`: "In every conifold frame (`def:set-frames`)" added (rule 7: `D_0` is defined in a frame, `thm:gau-R0`); the proof
  fixes a frame and notes that `Phi^SP(z, 1)/h(z)` does not depend on it (`def:con-PhiSP`, `thm:con-Pb`).
* **E13** Proof of `thm:gam-ginf`: "the formal version of [DLMF, 5.11.8]" removed; BIBPLAN allows 5.11.8 only in `rmk:fsp-classical`, to
  which the text now points ("no asymptotic statement is used"). 24.4.1 stays (derived in the proof of `lem:fsp-dlmf5118`).
* **E14** The display `c(z; 0) = prod_j Gamma_St(zeta^j/z + 1/5)/[(2 pi)^(5/2) z^(3/2) e^(phi(0)/z) e^(-phi'(0)/5) FF(0,0)]` was stated without
  proof (LF-C5 gives the normalisation only): the branch condition and a one-line derivation from `log V_j` added
  (`log Gamma_St(w + theta) = log V_j + (w + theta - 1/2)(log zeta^j - log z) - w + (1/2) log 2 pi`; I checked every term).
* **E15** After `cor:gam-poly`: "`ph^SP_n` has no negative power of `u`, that is `deg c_n <= 4n`, is `thm:con-Pa` (iv)" replaced by the
  precise statement (`ph^SP_n` lies in `Q[[u]]` by `def:con-PhiSP` and is a polynomial by `thm:con-Pa` (iv); with `c_n = L^(4n) ph^SP_n(L^-5)`
  this is `deg c_n <= 4n`). Same content as DD2 GAM-8 and RF-C3.
* **E16** `cor:gam-free`: the free coefficients of `Phi^act` are defined only once `Phi^act` has coefficients in `Q[[u]]`, which holds
  because `Phi^act = Phi^can` (`thm:gam-D`); said in the statement.
* **E17** `rmk:gam-blind`: "the series changes away from `L = 1`" is known only by computation; now said so, "not used" (rule 3).
* **E18** `rmk:gam-necessity` (a): the weights `lambda_j`, `l, d, e, Q, S` declared local (`lambda` is the logarithm of `-q` in s08-s09).
* **E19** `rmk:gam-general`: the condition that the `sigma`-dependent part of the powers of `z` be absorbed into `lambda` added (LA-C1;
  GEN's `Q = (-1)^(N + l_conc) z^(Sum c_f) q` in DD2 D6); "this is not used" added (rule 10).
* **E20** `\cite[s.~10]`, `\cite[s.~7, (11)]` written `\cite[\S10]`, `\cite[\S7, (11)]`, as in every other section.
* **E21** `a_formalsp.tex`, first paragraph: "the proof of record" (project wording) replaced by "the one the note relies on".

## 5. Open points

* None in the chain. P1 of DESIGN_CHANGES_2 s.7 can be closed: the formal route of LMI step 2 is written in full (REVIEW s.1).
* Outside the chain, unchanged and planned: the `\draftnote`s of `rmk:gam-sevenvar` (sketch; full proof would need `lem:fsp-delta`),
  `rmk:gam-Mcomplex` (non-integer exponent not in appendix A), `rmk:gam-prior` (neutral description of [MR1901075, \S10]).
* Not changed, for the editor: appendix A writes Jacobian matrices `D phi`, `D psi`, `D v_c` with the letter of the operator `D`
  (different lemmas from `lem:fsp-generic`, the only user of `\Dq` there); standard and unambiguous, but strictly a second use of `D`.

## 6. Changes that touch other sections (no edit made there)

* `thm:gam-D` now states its frame scope. Its users: `thm:int-D` (s01; its "equivalently, in every conifold frame" already matches),
  `def:gau-Fg` (s03), s07 (intro, Step 1), `thm:mai-main` proof (s10: "Fix `g >= 2` and the frame. By `thm:gam-D` ..."), Corollary B and
  `cor:mai-natfact`, `cor:mai-actual` (s10). All consistent; nothing to change.
* No other statement, definition or label of s09 or appendix A changed meaning. `[DLMF, 5.11.8]` now appears only in
  `rmk:fsp-classical`, as BIBPLAN says; s01 (`sec:int-inputs`) names only 24.4.1, consistent.
* Optional for the writer of s01: `s01_intro_3` says "the implicit function theorem only in its formal version (`lem:fsp-ift`)"; still
  true (Fubini; `lem:fsp-generic` in Theorem P(a)); one may add that LMI step 2 needs none.
* s08: `cor:con-A2` points to `cor:gam-equiv` inside its statement (a pointer, not a proof step); no circularity.

## 7. Compilation and checks

* Scratch copy of `main.tex`, `preamble.tex`, `macros.tex`, `refs.bib`, `sections/`; pdflatex, bibtex, pdflatex twice: no error, no
  undefined reference, no multiply defined label, no rerun request; 153 pages. Overfull boxes 14 (15 before; the two in s09 are gone,
  appendix A has none). BibTeX: one warning, "missing pages in ABDKS25", pre-existing (BIBPLAN, "cosmetic").
* Own exact arithmetic (Python, Fractions; evidence for this review only, not for the note): `B_6(1/5)/6 = 4537/3937500`; the printed
  `B_6(theta)`; `beta_(2k-1)(0) = bk_k`, `beta_(2k)(0) = 0` for `n <= 8`; `5 beta_(5l)(1/5) = (-1)^(l+1) B_(5l+1)(1/5)/(l(5l+1))` for `l <= 3`;
  `h_1 = bk_1 T_1 = -3/20`.
* Every file edited stays below 250 lines (s09: 90, 87, 99, 118; appendix A: 54, 107, 91, 84, 106).
