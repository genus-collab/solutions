# Draft notes

Every `\draftnote` in the note, by category, with file and line (line numbers as of the assembled version). The notes are visible in the PDF in red; set `\showdraftnotesfalse` in `main.tex` to hide them. Remove each note when the point is settled.

## 1. Citations to confirm in print

### [BS] pointers (cited by content; the label in the [BS] source to be checked)

| File:line | Nearest label | Pointer |
|---|---|---|
| s01_intro.tex:97 | cor:int-conj | conj: rigidity |
| s01_intro.tex:111 | sec:int-inputs | prop:ResConifold, eq:resconGWref |
| s01_intro.tex:150 | sec:int-BS | defn: refined GW invariant; conj: rigidity; prop: rigidity unref limit; prop:ResConifold |
| s02_setup.tex:18 | eq:set-wts | compare the data of [BS, Appendix A] for the tangent weights |
| s02_setup.tex:64 | def:set-inv | defn: refined GW invariant, defn: refined GW series |
| s02_setup.tex:209 | prop:set-frob | compare BS App. A for the Frobenius data (the note derives them from the genus-0 potential, lem:ide-g0) |
| s05_identification.tex:76 | lem:ide-fixed | defn: refined GW invariant (definition by localisation on the proper fixed locus) |
| s05_identification.tex:375 | prop:ide-rows | compare BS, Appendix A, for the Frobenius data at generic weights |
| s05_identification.tex:499 | rmk:ide-scope | lem: GW contravariance |
| s06_regularity.tex:247 | rmk:reg-betavalues | eq:betapol |
| s07_rigidity.tex:62 | lem:rig-poly | defn: rigid; the sentence that the resolved conifold is equirigid |
| s07_rigidity.tex:70 | lem:rig-poly | lem: shape of refGW series; lem: GW contravariance |
| s07_rigidity.tex:210 | lem:rig-g1 | prop:ResConifold, eq:resconGWref |
| s07_rigidity.tex:272 | (statement of the conjecture, before cor:rig-conj) | conj: rigidity; defn: rigid; defn: refined GW series |
| s07_rigidity.tex:293 | cor:rig-conj | eq:resconGWref |
| s07_rigidity.tex:293 | cor:rig-conj | lem: GW contravariance; the identification R_G = R_T^W used in the definition of the invariants |
| s07_rigidity.tex:300 | rmk:rig-BS | prop:ResConifold |
| s07_rigidity.tex:305 | rmk:rig-BS | prop: rigidity unref limit; eq: res conifold unref |
| a_evidence.tex:55 | sec:evi-loc | the displayed five-torus invariants of degree one in the resolved-conifold section, eq: res con GW 2 |
| a_evidence.tex:67 | sec:evi-rig | eq:betapol |

### Other citations

| File:line | Nearest label | Note |
|---|---|---|
| s01_intro.tex:109 | sec:int-inputs | Teleman is cited with the numbering of the arXiv version; the journal numbering has not been compared. |
| s01_intro.tex:110, s04_leaves.tex:47, s06_regularity.tex:8 | sec:int-inputs, rmk:lea-cy, sec:reg | LP2 numbering provisional (Theorem 8.1, Corollaries 8.2-8.4, Theorem 9.1, Theorem 9.3). |
| s05_identification.tex:78 | lem:ide-fixed | The fundamental-class and splitting axioms of the virtual class are used in standard form; a reference (Behrend-Manin) is to be added to the bibliography. |
| s06_regularity.tex:376 | rmk:reg-genus | A proof of M_irr(2) independent of the Marino-Vafa formula exists in the sources but is not among the sources of this note; the authors may wish to cite it, so that genus 2 rests on thm:reg-VI alone. |
| s07_rigidity.tex:367 | rmk:rig-F0 | The value of F_g(0) on T_F is quoted from the sources (derived there from [GP99] and [FP00, Theorem 4]); the derivation is not reproduced. |

## 2. Author decisions

| File:line | Nearest label | Decision |
|---|---|---|
| s01_intro.tex:5 | sec:int | Working title; to be settled. |
| s07_rigidity.tex:356 | prop:rig-flop | The description of the flop and its tangent weights is stated as in the source; add a reference or a line of justification for the weights of the two small resolutions. |
| s07_rigidity.tex:362 | prop:rig-flop (iii) | The identification of the graph sum at 1/Q = 0 with the Q' = 0 value of the frame-B sum of the flop Z' is one sentence in the source; the note spells out why the factors depend on the eigenvalues and the pairing only through the rows, the norms Delta_i and the orthogonality of the eigen-covectors. To be checked. |

## 3. Wording and steps to check (written out beyond the source, or where sections had to be matched)

| File:line | Nearest label | Note |
|---|---|---|
| s03_rows.tex:566 | lem:row-betaphase | The Stirling step is written out (the source gives it in condensed form): the series for log Gamma(a/z) + log Gamma(b/z) - log Gamma((a+b)/z), its tail, and the identification with the exponent of h_{(-a,-b,a+b)} through the oddness of 1 - 2k. |
| s03_rows.tex:601 | rmk:row-betaroutes | The relation between the two normalised expansions (amplitude 1 in s versus amplitude 1 in y = log U) is a consequence of the two statements as proved, not a step of either proof; the sources do not state it, which is why prop:row-beta is proved separately. |
| s05_identification.tex:74 | lem:ide-fixed | The sources say that R pi_* f^* N is moving on the fixed loci without listing the weights; the list given is the normalisation-sequence computation read for O(-1)_{alpha_1} + O(-1)_{alpha_2}. To be checked. |
| s05_identification.tex:410 | prop:ide-rows (ii) | The sources conclude Row 2 "at every Q with distinct roots not in {0, 1, 1/Q}"; the passage to the function field is the principle (P) with the parameter t in place of tau. Check that prop:row-P2 is stated with a general parameter (it is: prop:row-P2 is stated for an arbitrary analytic parameter). |
| s05_identification.tex:453 | prop:ide-frameB | The sources state "Hg_w (x) omega = h_w . omega (h_w is scalar)" in one phrase; the projection-formula sentence expands it. To be checked. |
| s07_rigidity.tex:70 | lem:rig-poly, last paragraph of the proof | The sources state the restriction of the invariants to a subtorus for T_CY (and, in the three-torus argument, for T_3, where the normal weights of the zero section are the four fibre weights of N_cf). The note writes out the general statement for a subtorus T' of the five-torus, with the hypothesis that no weight of R pi_* f^* N_cf on the fixed loci vanishes on T' (the weights are (fibre weight at p_pm) + j alpha_0/k, 0 <= j <= k); the four fibre weights alone do not suffice when alpha_0 != 0 on T' (alpha_1 + j alpha_0/k can vanish although alpha_1, B do not). This hypothesis is the one used in cor:rig-conj and in s01 (cor:int-conj) for the subtori of T_CY. To be checked. |
| s07_rigidity.tex:240 | eq:rig-loc | The degree-one localisation formula is the standard one with the listed ingredients; its normalisation agrees with the genus-zero value. Confirm the sign conventions of the node and edge factors. |
| s07_rigidity.tex:380 | lem:rig-3Tid | The source specialises the five-weight identity in two steps (alpha_0 -> 0, then epsilon_2 -> alpha_1 + alpha_2 - epsilon_1); the note's Theorem B is over K_CY and the ring map alpha_0 -> 0 applied to it is the same specialisation. Check that the wording of prop:ide-main and rmk:ide-alt supports this. |
| s07_rigidity.tex:383 | lem:rig-3Tid | The restriction of the invariants to T_3 for d >= 1 is stated in the source as a consequence of polynomiality without comment; the sentence on the T_3-fixed locus follows the d = 0 argument and should be checked. |

## 4. Source slips corrected in the note (to be back-ported to the source files)

The first three slips were adjudicated by the fidelity check of 2026-11-03 (the note is right); their draft notes were removed from the text and the rows below are the record. The fourth (s02_setup.tex:140) still carries its draft note. The fifth was found by the fidelity check and corrected in the text without a draft note.

| File:line | Nearest label | Slip in the source | Corrected to |
|---|---|---|---|
| s04b_closedforms.tex:108 (draft note removed) | thm:lea-loops, eq:lea-loops | LC_CY_PROOF.md, Thm 11.4_CY: the second (tau-variable) form of the loop integral is printed without the factor (1 - v) | e_{il}(t) = -(il eps1 eps2 B/(E^2 t^2)) Int_0^tau v (1 - v) s^*(v) L_2(v) (dA_i/dv)(dA_l/dv) dv: one (1 - v) from each A' = (1 - v) dA/dv and one (1 - v)^{-1} from du = dv/(1 - v) |
| s04b_closedforms.tex:155 (draft note removed) | rmk:lea-low | LC_CY_PROOF.md, s.11.4: g_1 printed as -(u/2) lam_1 lam_2 | g_1 = (u/(2E))(2 lam_2 - lam_1) = +(u/2) lam_1 lam_2 by (LAM5) |
| s04b_closedforms.tex:178 (draft note removed) | cor:lea-a0 | LC_CY_PROOF.md, Cor 12.1: the leg rescaling printed with the exponent alpha_2^{-2-b} | f_{j,b}(psi) = kappa alpha_2^{-1-b} f^F_{j,b}(psi/alpha_2), from E_{1j}(a,b) = kappa alpha_2^{-1-a-b} E^F_{1j}(a,b) (step (6) of the proof) |
| s02_setup.tex:140 | lem:set-crit (i) | LC_CY_PROOF.md (C2'): the multiplied-out critical equation is printed with an overall minus sign | (1 - QU)(1 - U) U W_U = +P(U); the sign is immaterial for the critical equation |
| s04_leaves.tex:208 (no draft note) | prop:lea-row2, the sentence after it | LC_CY_PROOF.md s.4.6, the F-form of Row 2 after (4.7): printed as z (1 - tau) (PsiF1)_tau = PsiF1 - a1 a2 F, i.e. without the coefficient a2 of PsiF1 | z (1 - tau) (PsiF1)_tau = a2 PsiF1 - a1 a2 F, because z (1 - tau) (PsiF1)_tau - a2 PsiF1 + a1 a2 F = (1 - tau) z^2 Gauss[F]; the coefficient is -[omega - z (1 - tau)/tau - (E/tau - a0)] = A - a0 = a2 = L_2 + dphi'. Re-derived and checked numerically on 2026-11-03; the statement is not used later in the note |

## 5. Removed by the assembler

- main.tex abstract: the placeholder note "Abstract to be written by the assembler" was replaced by the abstract.
- s01_intro.tex: the stub note was replaced by the section.
- s02_setup.tex, former subsection "The main statements": the note "the assembler may keep one copy only" and the four duplicated theorems (thm:set-A ... thm:set-D) were removed; the statements now live in s01 (thm:int-leaf, thm:int-ident, thm:int-reg, thm:int-rig, cor:int-conj) and the subsection is a pointer paragraph.

## 6. Closed by the fidelity check (2026-11-03)

- s04b_closedforms.tex: the three draft notes on the source slips (loop integral factor (1 - v), sign of g_1, leg exponent -1-b) were removed; the mathematics stands as written (see section 4).
- s07_rigidity.tex:272: the request to check the normalisation epsilon_pm = (epsilon_1 pm epsilon_2)/2 was removed from the draft note; [BS] defines epsilon_pm = (epsilon_1 pm epsilon_2)/2, and the note's wording matches.
- Not applied: the proposal to strengthen the draft note at eq:rig-loc (s07_rigidity.tex:240) with "the g = 1 output of the formula agrees with the genus-1 part of the degree-one five-torus invariant displayed in [BS]". In the [BS] source of record the displayed five-torus invariants of degree one are commented out, and the genus-1 display there reads (1/(12 e1 e2))(-eps_-^2 + eps_+(eps_+ + a0 + a1 + a2)) = e1 e2/12 + (e1 + e2)(a0 + a1 + a2)/24, whereas the formula of the note gives e1 e2/12 - (e1 + e2)(a0 + a1 + a2)/24 (the value that restricts correctly to [BS]'s formula on T_F). The two agree only up to the sign convention of epsilon_1, epsilon_2 (or of the alpha_i), so the comparison does not confirm the normalisation; the draft note stands as it was.
