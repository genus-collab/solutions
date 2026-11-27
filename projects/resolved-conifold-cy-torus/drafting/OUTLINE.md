# Outline of the note

Working file for the section writers. LCY = cy_torus/LC_CY_PROOF.md, RIG = cy_torus/RIGIDITY_PROOF.md, LSP = leaf_proof/leaf_structure_proof.md, LAP = lemmaA/LEMMA_A_PROOF.md, LP2 = lp2_explore/LP2_REGULARITY_PROOF_draft.md, CP = companion LaTeX paper lp2_explore/paper/ (main.tex, sections/*.tex, NOTATION.md, LABELS.md), FEAS = cy_torus_explore/feasibility_result.json, BS = BS_v2_2026-11-02/main.tex (read only). Notation: NOTATION.md. Labels: LABELS.md.

## 0. Binding rules

1. Transcribe mathematics, invent none. Every statement and every proof step comes from LCY, RIG, LSP, LAP (statement of Theorem L only), LP2 s.7.1 and s.7.6 (statements and the proofs of (H0), (S), (C)), or the companion paper CP (LaTeX of the rank-free arguments). If a step is missing from the sources, write the surrounding text and put `\draftnote{...}` saying exactly what is missing. Do not fill the gap.
2. Every all-orders proof step is written out in sentences, not sketched. "Verbatim" in LCY means: transcribe the LSP step with the substitutions of s.0.3 below; it never means "cite LSP". The note is self-contained apart from the cited literature and the two cited theorems of [LP2].
3. Finite-order checks appear only in the appendix, as evidence. No proof step cites the appendix. No proof step says "checked numerically".
4. Cite with pinned theorem numbers (LABELS.md, bibliography table). [LP2] numbers are provisional: mark the first citation in each section with `\draftnote{LP2 numbering provisional}`.
5. No project bookkeeping: no rounds, referees, verdicts, tags [P] [C] [E], status tables, file names, script names, code names (rc_lib.m, hvc_lib.m, wickZ, LofU, cconj, rmatrixAtClosed, graphValueC, ffcClosed, FEAS, ...), "proof of record", "the brief", "the critic", AI, dates of writing. Where LCY says "the code's row", write "the point-frame row" or "the normalised expansion".
6. BS: cite by content, never by line number; `\draftnote{BS label: ...}` after each pointer (NOTATION.md s.10). No errata of BS. Where the note's reading differs from BS, neutral wording (NOTATION.md s.9). The note proves the conjecture of BS for the resolved conifold; say what is covered and what is not (RIG s.5.6) without evaluating BS.
7. Plain prose; one notation (NOTATION.md); macros from main.tex; new macros only at the top of the section file. `\draftnote` for anything the authors must check (choice of wording, a citation to confirm, a convention).
8. The note proves the leaf structure at generic weights (Theorem A). The fibrewise theorem of LSP is not stated separately: it is the case alpha_0 = 0, alpha_i = epsilon_i, mentioned in `rmk:lea-fibrewise`.
9. Hypotheses are stated with their exact scope: "at the generic point of T_CY, or at any weight point with (ND)"; "per colouring"; "as formal series in z, w, psi to all orders in t". The regularity statements (s06) are statements about the decoration sum as defined in s02; their identification with F_g is s05 and is used only in s07.
10. Conditionality. The note's results are conditional only on: Teleman's classification [Tel12]; [PPZ15] for the R-action; [GP99]; Mumford's relation and formula [Mum83], [FP00]; [FP00, Theorem 4] and p. 3; [CG07, App. 1]; Stirling's series with remainder [DLMF 5.11.1]; the two theorems of the companion paper [LP2] (Theorem L, Theorem VI), which rest on the Marino-Vafa formula [LLZ03]; and, for genus 1 and d >= 2 only, the formula of [BS] on the subtori alpha_i = epsilon_j. Say so in s01 and at each place where an input enters. Nothing uses topological recursion, a mirror theorem, or a finite-order computation.

### 0.1 Global dependency graph

    s02 definitions and genus-0 data (uses nothing)
    s03 Euler integral and closed form: uses s02 (lem:set-N, lem:set-eig, def:set-rows, def:set-rho)
    s04 leaf structure: uses s02, s03 (thm:row-main, cor:row-F2, cor:row-reg, prop:row-transport, lem:row-hess, lem:row-betaphase, lem:row-canon)
    s05 identification: uses s02, s03 (Laplace calculus prop:row-P2, prop:row-P3, def:row-Lnorm, lem:row-betaphase), s04 only for the remark that the frame-B rows are unitary (thm:lea-qde); not Theorem A
    s06 regularity: uses s02 (def:set-fpot), s04 (thm:lea-main, cor:lea-alllegs, rmk:lea-use), [LP2]; not s05
    s07 rigidity: uses s05 (prop:ide-main, lem:ide-terms), s06 (cor:reg-main), s02 (lem:set-disc, lem:set-limits), s04 (cor:lea-a0 for the three-torus remark)
    appendix: points to the above, proves nothing

Keep visible: Theorem A (s03-s04) uses one analytic input, Stirling's series through `prop:row-beta`; the regularity theorem (s06) uses Theorem A and the two theorems of [LP2]; Teleman, PPZ, GP enter only through s05, i.e. only in the identification with F_g and hence in Theorem D.

### 0.2 Order of the sections and of the statements

Results are used only after they are stated and proved. Items that change place relative to the sources are marked **moved**.

### 0.3 The substitution dictionary (from LCY header and s.13.1; apply to every LSP step transcribed)

LSP is written on T_F with parameters x, c = 1 + x. LCY gives for each step whether it is used verbatim, by substitution, or with a new argument. The note writes every step at generic weights. Dictionary:

| LSP | Note |
|---|---|
| (x, 1, 1) in the Euler phase log s + x log(1-s) - log(1 - tau s) | (alpha_1, alpha_2, A): phi_tau(s) = A log s + alpha_1 log(1 - s) - alpha_2 log(1 - tau s) |
| c | E |
| x/tau, c/tau in the quantum product and A | alpha_1 alpha_2/tau, E/tau - alpha_0 |
| -x c^3 | -epsilon_1 epsilon_2 E^3 |
| x^3 (Delta_2(0) c) | epsilon_1 epsilon_2 alpha_1 alpha_2 A B/E |
| class (c, -x, -1) | (E, -epsilon_1, -epsilon_2) |
| rows h_w^{-1} r^pt | h_w^{-1} h_nw r^pt = h_{(alpha_1, A, -E)} r^pt (`lem:set-hmult`) |
| 1 - s = x s (1 - tau s) (LSP Lemma 3.1(b)(ii)) | the two forms of `lem:set-N`(ii): (1 - tau) alpha_2 s = (A - B s)(1 - tau s) and (E - B tau s)(1 - s) = alpha_1 (1 - tau s) |
| double zeros (1 + jz)^2 of the numerator of F_n | simple zeros at -alpha_2/j and -A/j; only degree and leading coefficient are used (`lem:lea-taucoef`) |
| x^2 - c x L_2 = -x^3 tau s*^2 | alpha_1 alpha_2 - E L_2 = -B tau s^* L_2 (`lem:set-Delta2`) |
| h_nw = 1 | the multiset identity (3.9a) (`lem:set-hmult`) |
| delta^2 = c^2 - 4 x tau | delta^2 = (E + alpha_0 tau)^2 - 4AB tau |
| s* = 2/(c + delta), s~ = (c + delta)/(2 x tau) | s^* = (E + alpha_0 tau - delta)/(2B tau) = 2A/(E + alpha_0 tau + delta), s~ = (E + alpha_0 tau + delta)/(2B tau) |
| U_2(0) = -1/x | U_2(0) = -alpha_2/B |
| mu = c/tau - 1 | mu = E/tau - A |
| Theta = -x^3 tau s*^2 | Theta = -epsilon_1 epsilon_2 B tau s^* L_2 |
| P_k = prod (k + 1 - i c)^2 | p_k = prod_{i<=k} ((k+1) alpha_2 - iE)((k+1) A - iE) |
| K_x = Q(x) | K_alpha (or K_CY) |
| v_2(z; 0) = h_w^{-2}(1, x/c) | r_2(z; 0) = h_reg(z)(1, alpha_1 alpha_2/E) |

Every place where LSP used a_2 = A or alpha_i = epsilon_i is listed in LCY s.13.1 ("No step uses ..."): the four replacements are in the table. A writer who finds another such place stops and puts a `\draftnote`.

## s01_intro (written last, by the assembler)

Covers: RIG s.0.3, s.0.4 (neutral), s.5.6; LCY s.0.2, Remark 0.1, s.0.5; the inputs list of rule 10; BS's conjecture by content.

Contents, in order.
1. Setting in words: Z = X x C^2, the five-torus, the Calabi-Yau four-torus T_CY, the series F_g, the conifold point Q = 1, frame B and the Hodge-vertex graph sum.
2. Statements: `thm:int-leaf` (Theorem A, with (F1)-(F3), (L1)-(L3) in words and the remark that the Calabi-Yau condition is not used), `thm:int-ident` (Theorem B), `thm:int-reg` (Theorem C: every coloured graph with an edge regular; polar part of the vertex graph exactly beta_{2g-2} t^{2-2g}; hence the polar part of F_g), `thm:int-rig` (Theorem D: F_g - F_g(0) = P_g Li_{3-2g}(Q)/(epsilon_1 epsilon_2), g >= 2, equivalently epsilon_1 epsilon_2 GW_{g,d} = P_g d^{2g-3}, at every point of Lie(T_CY)), `cor:int-conj` (with genus 0 and 1: the rigidity conjecture of [BS] for the resolved conifold, every d >= 1, for T_CY and every connected group acting through it). Each with exactly the hypotheses of its body version.
3. Inputs, listed once (rule 10), and which results use which: Theorem A uses Stirling; Theorem C uses Theorem A and [LP2]; Theorems B and D add Teleman, PPZ, GP, Mumford, FP, CG; genus 1 adds [BS] at d >= 2.
4. What is not claimed: nothing off T_CY (the regularity fails there; evidence only, appendix); the constant term F_g(0) (depends on alpha; not determined); the conjecture for group actions not conjugate into the five-torus; the full genus sum at d = 0 (not part of the conjecture); no human-independent check of the cited theorems. The proofs of the cited theorems are not reproduced. Not used: topological recursion, mirror theorems, any finite-order computation, any degeneration argument.
5. Strategy: Euler integral -> Gauss's equation -> KEY FORMULA -> leaf structure (s03-s04); Teleman at generic weights with Q = 0 semisimple on T_CY -> frame change (s05); string reduction and the vanishing theorem -> regularity (s06); poles of a rational function of Q: {1, infinity, Disc}, the alpha_0-adic argument on Disc, regularity at infinity, comparison with Li (s07). The whole polar part sits in the single-vertex graph. The flop remark and the three-torus.
6. Plan of the note.

`\draftnote` at Teleman: numbers from arXiv v3, journal numbering not compared (as in CP). `\draftnote` at [LP2]: numbering provisional. `\draftnote` on the title (working title).

Length: 15-20k characters.

## s02_setup

Covers: RIG s.0.1, s.1.1-s.1.5; LCY s.1.1-s.1.5 ((C1')-(C14'), s.1.3, s.1.4, s.1.5 as far as it is definitions), Lemma 3.1_CY(a), (b), Lemma 4.1_CY, Lemma 4.2_CY(a)-(d), Prop 4.3_CY(a)-(c) and the algebraic parts of (d), (e); LSP s.1.4 (rings, principles, leaf identities), Lemma 1.1, s.0.1 (graph sum), s.13.0 (definitions of the Hodge side, not (H0)); LAP s.1 for Xi and dg; CP s02_setup.tex def:set-fpot and s11_conifold.tex "Set-up" block for the LaTeX.

Order.
1. Target, five-torus, tangent weights `eq:set-wts`; T_CY, (ND) `def:set-ND`; the symmetries `lem:set-sym` (RIG s.1.1, with the reason: automorphisms of Z; the invariants are invariant). Fields `eq:set-fields`.
2. Invariants by localisation and F_g, d >= 0 `def:set-inv`, `rmk:set-d0` (RIG s.0.1; the T_CY-fixed loci are the alpha_0-fixed loci of the proper moduli space of maps to P^1). The constant-map term is named here and discussed in s07.
3. Superpotential, critical equation, roots, (1.3), L at the critical points, W_2 `lem:set-crit` (LCY (C1'), (C2'), Lemma 4.2_CY(a), (b); RIG s.1.3). Write the multiplication-out of the critical equation as in LCY (C2').
4. Frobenius data `prop:set-frob` (LCY Lemma 4.1_CY, with the full computation from the genus-0 potential PrepCl + Li_3/(epsilon_1 epsilon_2); say that the potential itself is derived in s05 `lem:ide-g0` and that only its third derivatives are used here; the classical part by localisation at the two fixed points as in LSP proof of Lemma 4.1(a)). State the convention: covectors paired by eta^{-1}, d/dt<a,b> = <a',b> + <a,b'>.
5. Idempotents and Delta_i `prop:set-idem` (LCY Prop 4.3_CY(a)-(c): the four forms, with the proof of (c) written out; the evaluation in the quadratic algebra by the two characters). The quadratic algebra and its involution `def:set-QA` (LCY (C3'), without code: the involution U -> U_1 + U_2 - U of K(Q)[U]/(P) exchanges the two characters; every factor of idempotent 2 is the image of that of idempotent 1 under the involution, because the scalar constants are the same at both idempotents, `lem:set-hmult` and `lem:ide-norm`).
6. Algebraic functions of tau `def:set-alg` (LCY s.1.3, including r_D and the case alpha_0 = 0), N(s) and its identities `lem:set-N` (LCY Lemma 3.1_CY(a), (b) with proofs), eigenvalues in terms of s^* `lem:set-eig` (LCY Lemma 4.2_CY(c), (d), with the expansion and l(0)), Delta_2 in terms of s^* and Delta_2'(0) = 0 `lem:set-Delta2` (LCY Prop 4.3_CY(d), (e), the parts that use only `lem:set-N`; the part of (d) with alg^2 is `prop:lea-key`).
7. The discriminant `lem:set-disc` (RIG s.1.4 (i)-(v) with proofs) and the behaviour at Q = 1, 0, infinity `lem:set-limits` (RIG s.1.5; the exact unit at Q = 1 is forward-referenced to (F1)).
8. Rows `def:set-rows`: the point-frame row through the Wick polynomials (LSP (C7) formula and the characterisation of wick_m, including the count "a term with j factors and v-degree 2n has z-order n - j"; LSP Lemma 3.9's homogeneity statement is in s03), the H-row (H), the covector r_i; frame B: r_{0,i} := h_{(alpha_1, A, -E)} r^pt_{0,i}, with `lem:set-hmult` showing that this equals h_w^{-1} h_nw r^pt_{0,i} (LCY (3.9a), including the one-line uniqueness of the class) and the forward reference to s05 for the meaning of h_nw (the Q = 0 constant). Dil_1(1) = -1.
9. Factors `def:set-factors` (LCY (C10')-(C12')), `lem:set-edgegen` (LSP Lemma 1.1 with proof), rho and the leaves `def:set-rho`, scope `rmk:set-scope`.
10. Coloured stable graphs and the graph sum `def:set-graphsum`, `eq:set-graphsum` (RIG (3.1); decorations as in LSP s.0.1 and (C10): exponents on half-edges, ordered tuples of dilaton degrees d >= 2 with weight 1/m!, the term d = 1 as the dilaton shift). Colour 1 singular, colour 2 regular.
11. Rings and principles `lem:set-rings` (LSP s.1.4 transcribed: O, M, C((z)), R_z, iota, refl, TRANSFER, (P) with LCY's list of denominators at tau = 0), leaf spaces and identities `lem:set-leafid` ((LAM1)-(LAM6) with their one-line proofs, (PF) with its proof, (T)).
12. Hodge side `def:set-fpot`, `def:set-mirr` (LSP s.13.0 definitions: Hg_w, base-point derivatives, Dom, F^pot_h on Dom, the dilaton-shifted form and the overlap, q_t in Dom at idempotent 1 (forward reference to `thm:lea-L1`(iv) and (F1)), not at idempotent 2; Xi_h, dg, Z[h, n], M_irr(h)). Reuse CP def:set-fpot LaTeX with E in place of E and t in place of t_c.

Dropped: every sentence about code, libraries, scripts, references to BS App. A beyond one `\draftnote{compare BS App. A}` at the Frobenius data.

Depends on: nothing. Used by: all.

Length: 30-40k characters.

## s03_rows

Covers: LCY s.3 entire (s.3.1-s.3.11, Remark 3.20 optional) and Lemma 4.2_CY(e); LSP s.3.4 entire (transcribed), proofs of LSP Prop 3.10, 3.12, Lemma 3.13, Thm 3.14, Prop 3.15, Cor 3.16, Thm 3.17 (transcribed with the dictionary); LSP Lemma 13.3's analytic paragraph (as `lem:row-betaphase`).

Order: `lem:row-euler` -> `lem:row-hess` -> `lem:row-wr` -> the Laplace calculus (`def:row-Acal`, `lem:row-gauss`, `def:row-Lnorm` with (D1)-(D5), `prop:row-P3`, `prop:row-P2`, `prop:row-P1`, `lem:row-laplace`, `lem:row-homog`; LSP s.3.4 in LSP's order, each proof in full) -> `prop:row-wickrow` -> `lem:row-amp` -> `prop:row-gaussode` -> `lem:row-unique` -> `thm:row-closed` -> `lem:row-betaphase` -> `prop:row-beta` -> `cor:row-reg` -> `thm:row-main` -> `cor:row-F2` -> `lem:row-eik` -> `prop:row-transport` -> `rmk:row-first` -> `lem:row-canon` -> `rmk:row-eulerint` (optional).

Instructions.
- `prop:row-wickrow`: LCY says "word for word the proof of LSP Prop 3.10"; transcribe LSP s.3.5 steps (1)-(2) (the row in y = log U as Lnorm[W(e^y), 1; y_1], the coordinate change hc(s) = log(1 - tau s), (P1)(iii), (D3); (ii) by Lemma 3.9 with lambda = tau, then `lem:row-wr`, then lambda = -1 and (D1) at tau = 0; (iii) by (D2)). Replace "the code's row" by "the point-frame row".
- `prop:row-gaussode`: transcribe LSP Prop 3.12's proof ((P2) twice, (D3), the amplitude identity, (P3), conjugation by Exp[dphi/z](f_2/f_2(0))^{-1/2}).
- `lem:row-unique`: the recursion (k + 1)(k + c) Y_{k+1} = (k + a)(k + b) Y_k with k + c = (E + (k + 2)z)/z a unit because E != 0; the sentence on the second exponent 1 - c = -1 - E/z (not used).
- `lem:row-betaphase` (new placement): state for a, b in a field of characteristic 0 with a b (a + b) != 0: the normalised expansion of the phase a log U + b log(1 - U) in y = log U at the critical point U = a/(a + b) equals h_{(-a, -b, a+b)}(z). Proof from LSP Lemma 13.3: both sides rational in (a, b) by (D2); for a, b > 0, f(y) = a y + b log(1 - e^y) strictly concave on (-infinity, 0), the Beta integral B(a/z, b/z + 1) [DLMF 5.12.1], `lem:row-laplace` with amplitude 1, Stirling's series with remainder [DLMF 5.11.1, s.5.11(ii)] for log B(a/z, b/z), B(a/z, b/z + 1) = B(a/z, b/z) b/(a + b), uniqueness of asymptotic expansions. Then rationality: an identity between rational functions of (a, b) valid on an open set of the real quadrant holds identically (LSP Prop 3.15's "then for all x by rationality"). `\draftnote{the Stirling step: write the series for log Gamma(a/z) + log Gamma(b/z) - log Gamma((a+b)/z) and identify its tail with -Sum_k bk_k z^{2k-1}(a^{1-2k} + b^{1-2k} - (a+b)^{1-2k}); the exponent of h_{(-a,-b,a+b)}}`. This lemma is the only analytic input of s03-s04, and it is also used in s04 (`prop:lea-regrows`(a)) and s05 (`lem:ide-norm`), s07 (`prop:rig-flop`).
- `prop:row-beta`: either transcribe LSP Prop 3.15 directly at the parameters (A, alpha_1): the Beta integral B(A/z + 1, alpha_1/z + 1) = Int_0^1 s^{A/z}(1 - s)^{alpha_1/z} ds, Laplace at s = A/E, Stirling, the elementary factor E/(E + z) from B(a + 1, b + 1) = B(a, b) ab/((a + b)(a + b + 1)); or derive it from `lem:row-betaphase` and (D4) as LCY does (phi_0 = A ph with ph the phase at x = alpha_1/A; note that LSP Prop 3.15's phase is in the variable s, not log U, so the writer must check which lemma applies and say so; `\draftnote` if the relation between the two phases is not written in the sources). Include LCY's z^1 hand check as a remark.
- `thm:row-main`: the three equalities r_{0,1} = h_{(alpha_1,A,-E)} r^pt_{0,1} = h_{(alpha_1,A,-E)} N = h_{(alpha_1,A,-E)} A_beta F X alg = (E/(E + z)) F X alg, and the coefficientwise meaning (`cor:row-reg`, `prop:row-wickrow`(iii), (T)).
- `lem:row-eik`, `prop:row-transport`: LCY s.3.10 in full, including the derivation of (3.19), the eikonal discriminant computation, Eik(dphi_tau) = 0 via N(s^*) = 0, P_tr = delta, R_tr and R_tr(0), and the hand proof of R_tr = delta (log alg)_tau ending at `lem:set-N`(ii).
- No pole at z = 0 has two proofs (`cor:row-reg` and `prop:row-transport`); keep both, as LCY does.
- Analytic/formal: mark each statement [formal] or [analytic] as LCY does, in words ("as an identity of formal series"; "as analytic functions on D").

Depends on: s02. Used by: s04 (all), s05 (Laplace calculus, `lem:row-betaphase`), s07 (`lem:row-betaphase`).

Length: 55-70k characters.

## s04_leaves

Covers: LCY s.0.2-0.3 (statements), Remark 0.1, s.4.4-s.4.8, s.5-s.11, s.12 (Cor 12.1), s.10 (assembly); the LSP steps that LCY marks verbatim, transcribed with the dictionary.

Order: statements `thm:lea-main`, `rmk:lea-cy`, `thm:lea-qde`, `thm:lea-closed` (the closed forms may be stated here and proved in the last subsection) -> `prop:lea-key` -> `lem:lea-op` -> `prop:lea-row2` -> `prop:lea-A` -> `prop:lea-K1`, `rmk:lea-wronsk` -> `lem:lea-taucoef` -> `lem:lea-tcoef` -> `rmk:lea-nopole` -> `prop:lea-pf` -> `thm:lea-L1` -> `prop:lea-sigma` -> `prop:lea-proj` -> `eq:lea-Ebar` -> `lem:lea-Yqde` -> `thm:lea-flow` -> `eq:lea-int` -> `thm:lea-L3` -> `cor:lea-rhoe` -> `lem:lea-conj` -> `lem:lea-reg` -> `lem:lea-regdata` -> `thm:lea-row2reg` -> `cor:lea-orth` -> `prop:lea-regrows` -> `cor:lea-K2` -> `rmk:lea-E22` -> `lem:lea-leg` -> `prop:lea-L2` -> `prop:lea-E12` -> `prop:lea-transp` -> `cor:lea-alllegs` -> `rmk:lea-use` -> proof of Theorem A and `thm:lea-qde`, `eq:lea-table`, `rmk:lea-inputs` -> closed forms `thm:lea-Aj`, `cor:lea-sumrule`, `cor:lea-top`, `thm:lea-loops`, `cor:lea-anti`, `rmk:lea-legs`, `rmk:lea-low` -> `cor:lea-a0`, `rmk:lea-fibrewise`.

Instructions.
- Each "verbatim [PR]" of LCY is a transcription task: the LSP proof with the dictionary. The list, with LSP locations: Prop 4.4(a)-(c) (LSP s.4.4); Prop 4.6 (LSP s.4.6, the conjugation (4.6) and steps (i)-(iii)); Prop 4.8 (LSP s.4.8, (a)-(d), the division by psi legitimate since (A - L_1) chi_1 = 0); Remarks 4.9, 4.10 (LSP s.4.8); Lemma 5.1(b), (c), Lemma 5.2, Prop 5.3, Thm 5.4 with remarks (LSP s.5); Prop 6.1 (LSP s.6.2, with LCY's k = 0 computation); s.7.1, Lemma 7.1, Thm 7.2 (i), (ii), s.7.4, Thm 7.3, Cor 7.4 (LSP s.7); Lemma 8.1 (LSP s.8.1, rephrased: no code), Thm 8.4 Steps A-D (LSP s.8.4), Cor 8.5, Cor 8.7, Remarks 8.8-8.10 (LSP s.8); Lemma 9.1, Prop 9.2, 9.3, 9.4, Cor 9.5, s.9.6 (LSP s.9); Thm 11.4, Cor 11.2(f) (LSP s.11). Everything else is written out in LCY.
- `prop:lea-key`: (F1) for Delta_1 is here (t^2 rho^2 Delta_1 = -epsilon_1 epsilon_2 E^3 exactly), together with the alg^2 form of Delta_1 (LCY Prop 4.3_CY(d) first form, Prop 4.4_CY(a)).
- `thm:lea-row2reg`: LCY gives the proof via Gauss's operator at s~ (Steps A-D) and notes that `prop:ide-rows`(ii) proves the same via the Picard-Fuchs operator at indeterminate weights. Write the Gauss route in full here; one sentence pointing forward to s05 as an alternative.
- `prop:lea-regrows`(a) is an application of `lem:row-betaphase` with (a, b) = (-alpha_2, E): say so and do not re-prove.
- Sections 5-7 are over K_alpha ((P) of `lem:set-rings`); say so at the start of the subsection, as LCY does.
- `thm:lea-L1` remark (iv) gives (F3) for Dil_1 and the value q_0 at t = 0: keep.
- `rmk:lea-use` (LSP s.9.6): transcribe; it is what s06 step (ii) refers to.
- `eq:lea-table`: the table of (F3) with the values at t = 0 (LCY s.10) as a tabular; `rmk:lea-inputs`: where each input and each factor of (ND) is used (LCY s.10 last paragraph).
- Closed forms: LCY s.11 is written out except Thm 11.4 (transcribe LSP Thm 11.4) and Cor 11.2(f); the items LCY lists as "not re-derived" (branch-free forms of rho and alg, closed form of e^{dphi}, higher tables) are omitted, with one sentence saying which closed forms are not given.
- `cor:lea-a0`: LCY Cor 12.1 with proof (1)-(7) and remarks (i)-(iv); then `rmk:lea-fibrewise`: at alpha_0 = 0, alpha_i = epsilon_i the theorem is the fibrewise statement, with the dictionary of s.0.3 read backwards (c = 1 + x, etc.). Remark (iv) of Cor 12.1 (the direct proof covers alpha_0 = 0; only the citations of the Q = 0 constant need alpha_0 != 0 and extend by rationality) is important for s07's three-torus corollary: state it clearly.
- Sharpness statements: keep, they are proved (Cor 11.3, 11.5).
- Every "[Evidence: ...]" bracket of LCY is dropped here; its content goes to the appendix.

Depends on: s02, s03. Used by: s06 (`thm:lea-main`, `cor:lea-alllegs`, `rmk:lea-use`), s05 (`thm:lea-qde` in a remark), s07 (`cor:lea-a0`).

Length: 100-125k characters (the longest section; if the file passes 110k, split the closed forms s.11 into a second file `s04b_closedforms.tex` input from s04_leaves.tex, same label prefix).

## s05_identification

Covers: RIG s.3 entire (statement, inputs, proof read at T_CY, description of the coloured terms, alternative route); LSP Prop 13.2 statement, Lemma 13.3, proof (a) (including "Field of definition", "Descent", "Flatness" (F1)-(F9), the remark on Prop 7.7), (b)(i)-(iv), (c), (e) transcribed at indeterminate weights; LSP Remarks (ii), (v) after Prop 13.2; CP s11_conifold.tex prop:con-frame and lem:con-norm (LaTeX to reuse), CP s04_generic.tex prop:gen-flat, prop:gen-teleman, lem:gen-free, lem:gen-ambig, lem:gen-mumford, lem:gen-ch (LaTeX of the rank-free steps, to be specialised to rank 2 as LSP (F1)-(F9) prescribe).

Order: statement `prop:ide-main`, `eq:ide-main` (with the two sentences of RIG s.3 on what Kw is and why part (d) of LSP is not needed) -> `lem:ide-fixed` -> `lem:ide-g0` -> `lem:ide-ss` -> `prop:ide-teleman` with `lem:ide-free` -> `lem:ide-ambig` -> `prop:ide-flat` (Fl1)-(Fl9) -> `rmk:ide-tel77` (optional) -> `lem:ide-norm` -> `prop:ide-rows` (i)-(iv), `eq:ide-PF` -> `lem:ide-mumford` -> `lem:ide-ch` -> `prop:ide-const` -> `prop:ide-frameB` -> `lem:ide-terms` -> proof of `prop:ide-main` assembled -> `rmk:ide-alt` -> `rmk:ide-scope` -> `rmk:ide-frame`.

Instructions.
- The proof is LSP Prop 13.2 (a), (b), (c), (e) at the five indeterminate weights (field K_alpha), followed by the observation (RIG s.3) that on T_CY (field K_CY) every hypothesis is a non-vanishing which holds: list them exactly as RIG does ((a): fixed loci unchanged, alpha_0 != 0; semisimplicity L_1(0) - L_2(0) = alpha_0; the field has characteristic 0 and cardinality 2^aleph_0; (Fl8) the ratio Delta_+(0)/Delta_-(0) non-constant; (b): the limiting pairs (-alpha_2, A), (alpha_0, alpha_1) with a, b, a + b non-zero; (c): the weights non-zero; (e): w has non-zero entries). Write the proof once, over K_CY directly, and remark that the same text holds over K_alpha; or write it over K_alpha and specialise (RIG's alternative route `rmk:ide-alt`). Choose the first (RIG's main route) and give the alternative as a remark.
- `prop:ide-flat`: LSP's (F1)-(F9) give, for each step, the reason it does not see the rank and the rank-2 data; CP prop:gen-flat gives the full text of each step in rank 3. Write each step in full in rank 2 (the edgeless form, the pushforward identities on the universal curve, the bubble term, the divisor equation for P^1 at every d including d = 0, the defect matrix and its g-independence, stability [Tel12, 2.18, 2.19, (2.15)] with the degree count g > 4n + 1, skewness of M_0 from Psi^T eta Psi = 1 and symplecticity, non-degeneracy (Fl8), the equation and its row form z D_i r_i = (A - L_i) r_i). Rename LSP's (F1)-(F9) to (Fl1)-(Fl9) to avoid the clash with Theorem A's (F1)-(F3).
- `prop:ide-teleman`: the embedding of K_CY((Q)) into C, the base change, the descent by Aut(C/K_CY((Q))) with the freeness lemma (R_m diagonal, odd, entries killed by Mumford and ch_{2k-1} != 0; the lemma is stated in general and its proof uses `lem:ide-mumford`, `lem:ide-ch`, which therefore are stated before it: **moved** forward relative to LSP, as CP does).
- `lem:ide-norm`: LSP Lemma 13.3 with the limits of W_k(U_+-) at Q = 0 written out, then `lem:row-betaphase` for the value; the multiset identity h_{wt_+} h_{(-alpha_2, A, -alpha_0)} = h_{wt_-} h_{(alpha_0, alpha_1, -B)} = h_nw.
- `prop:ide-rows`(ii): Euler's identity PF e^{W/z} = theta_U(k_m e^{W/z}) with k_m = -z alpha_1 Q(1 - U)/(1 - QU), "multiply out": write the multiplication out (`\draftnote` if it does not close in one page: the sources say it is a direct identity in U and Q). Then (P2) twice, (D3), (P3) give Row 2 "at every Q with distinct roots not in {0, 1, 1/Q}".
- `prop:ide-const`: LSP (c) in full, including the passage from M_k to M_{k,1} by the flat unit and pi_*(psi_1 pi^* y) = (2k - 2) y.
- `prop:ide-frameB`: LSP (e) in full: R'.(Hg_w (x) omega) = (R' h_w).omega by [PPZ15, Prop 2.14]; the exponent -1 is derived.
- `lem:ide-terms`: RIG's description, with the denominators listed (W_2(U_i), 1 - U_i, 1 - QU_i, y, BQ, constants), the deck involution, the colour sum in K_CY(Q), regularity at both points over Q = 0 (unramified by `lem:set-disc`(v)), expansion Sum GW_{g,d} Q^d. This lemma is what s07 `lem:rig-poles` cites.
- Remarks: RIG's bracketed references to BS App. A become one `\draftnote{compare BS App. A}`; LSP's remarks (i), (iii), (iv) on the code and on T_F are dropped except the mathematical content of (iv) (on T_F and on alpha_0 = 0, Q = 0 is not semisimple: the identification there is `lem:rig-3Tid`).
- Theorem A is not used here; `thm:lea-qde` is cited in `rmk:ide-frame` only.

Depends on: s02, s03 (Laplace calculus), `lem:row-betaphase`. Used by: s07.

Length: 45-60k characters.

## s06_regularity

Covers: LSP s.13.0 (H0) with proof, Lemma 13.0; LP2 s.7.1 (S_P2), (C_P2) with proofs, (L0_P2) as a remark, s.7.2 (remark); LAP Thm L, Cor L1, Cor L3, s.7.5 (statement in the w-convention); LP2 s.7.6 Theorem VI (statement); RIG s.4 entire (Steps 0-6, specialisation, hand check of (0.1)); RIG (0.1); LSP s.13.1 Theorem R (i)-(v) and s.13.2 "By genus"; CP s10_regularity.tex lem:reg-H0, lem:reg-string, lem:reg-C, lem:reg-L0, lem:reg-range and s11_conifold.tex cor:con-reg proof (LaTeX to reuse, with c -> E, x c^3 -> epsilon_1 epsilon_2 E^3).

Order: `lem:reg-H0` (statement and the three-step proof) -> `lem:reg-range` -> `lem:reg-string`, `lem:reg-C`, `rmk:reg-L0` -> `rmk:reg-nol0` -> `thm:reg-L` (cited; statement in the note's convention: for every h >= 2 and every non-empty multiset M of integers >= 2 with dg(M) <= 2h - 3, Xi_h(M) = 0, where Xi_h is built from Hg_w, w = (E, -epsilon_1, -epsilon_2), E = epsilon_1 + epsilon_2, at the generic point of T_CY and at every point with epsilon_1 epsilon_2 E != 0, the diagonal epsilon_1 = epsilon_2 included; the scale law Xi_h(U) = E^{-(3h-3+2r)} Xi^tau_h(U; tau_w = -epsilon_1/E) and sigma_w = -epsilon_1 epsilon_2/E^2 != 0; cite [LP2, Theorem 8.1, Corollaries 8.2, 8.3, 8.4]; one sentence: the leaf scale E is the first weight of w, which is what the theorem requires, and this is the one place where the Calabi-Yau condition is used) -> `thm:reg-VI` (cited: [LP2, Theorem 9.1]: (-epsilon_1 epsilon_2 E^3)^{g-1} Z[g, 0] = beta_{2g-2}(epsilon_1, epsilon_2) for g >= 2, identically in (epsilon_1, epsilon_2) and at every specialisation with epsilon_1 epsilon_2 E != 0, because Z_g is a polynomial in 1/sigma_w (RIG s.4 Step 5 explanation)) -> `lem:reg-beta`, `eq:reg-beta` (RIG (0.1) with the proof "both sides are (2g-3)!/(epsilon_1 epsilon_2) times [t^{2g}] S(epsilon_1 t) S(epsilon_2 t)"; P_g symmetric polynomial in epsilon_1^2, epsilon_2^2; the hand values at g = 2, 3 as an example) -> `thm:reg-main` (statement with hypotheses M_irr(h) for the stated range, as LSP Theorem R, so that the use of each instance is visible; proof Steps (i)-(v) from RIG s.4 Steps 0-5 and LSP s.13.1, every sentence; Step 0 says that q_t and rho q_t lie in Dom by `thm:lea-L1` and (F1)) -> `cor:reg-main` (a)-(c) (the hypotheses discharged by `thm:reg-L`; Step 6 for (c): each coloured term lies in K_CY(C) by `lem:ide-terms`... **careful**: (c) for F_g uses `prop:ide-main`; state (c) in two halves: for the graph sum (no s05) and for F_g (with s05), as RIG does) -> `rmk:reg-special` -> `rmk:reg-genus`.

Instructions.
- The statements of (S) and (C) are at scale E with class Hg_w; LP2 s.7.1 is written exactly in this convention (E, epsilon). Transcribe including the string-equation computation S lam_j = (j/E) lam_j, S T = -psi lam_1.
- `lem:reg-H0`: LSP s.13.0 proof (Steps 1-3) or LP2 s.7.1 (identical); CP lem:reg-H0 has the LaTeX. Use the formal variable name `\mathsf{r}` for the auxiliary variable of Step 3 (CP uses it), not vartheta.
- `thm:reg-main` Step (ii) refers to `rmk:lea-use` for how the factors attach; Step (iii) uses `lem:reg-string`, `lem:reg-C`, and M_irr(h); Step (iv) is the counting with the two sentences on genus-0 and genus-1 singular vertices and regular vertices ((F1), (F3)); Step (v) the single vertex.
- Say once (RIG s.4 Status): as a statement about the decoration sum, Theorem C uses Theorem A and the two cited theorems of [LP2]; its identification with F_g is `prop:ide-main`.
- "Regularity at t = 0 after summing over decorations" means: the Laurent series in t of Cont(Gamma, c) has no negative powers; define this once.
- Dropped: RIG's status lines, the parenthetical references to the feasibility reports, LP2's audit s.7.7, LSP Remark 13.1 (the Reduction Lemma), everything on TR.

Depends on: s02 (`def:set-fpot`), s04, [LP2]; s05 only for the F_g half of (c). Used by: s07.

Length: 25-35k characters.

## s07_rigidity

Covers: RIG s.5 entire, s.6.1-6.3, s.7 (adapted), s.0.4 (neutral), s.5.6; LSP Prop 13.2(a) Aspinwall-Morrison for genus 0; CP prop:gen-poly for the polynomiality LaTeX.

Order: statement `thm:rig-main`, `eq:rig-main` -> `sec:rig-plan` (RIG s.5.1: the four steps and which use what; G_g := F_g - F_g(0)) -> `lem:rig-poles` (RIG Lemma 5.1, from `lem:ide-terms`; the paragraph "What is actually proved about the rows") -> `lem:rig-poly` (RIG Lemma 5.2 in full: H^0(C, f^* O(-1)) = 0, the rank of R^1, the class e_T(...)^{-1} as (1/(epsilon_1 epsilon_2)) times an honest class, [GP99] on the proper space with a generic one-parameter subgroup, degree 2g, symmetries, restriction to T_CY with the reason no denominator vanishes; BS's shape lemma and contravariance cited by content with `\draftnote{BS label: lem: shape of refGW series; lem: GW contravariance}`) -> `rmk:rig-d0` (the constant-map term as a rational function with denominators alpha_0 alpha_1 alpha_2 epsilon_1 epsilon_2 and alpha_0 A B epsilon_1 epsilon_2; constant in Q) -> `lem:rig-root` (RIG Lemma 5.3, with the simple-root lifting and the first coefficients) -> `lem:rig-disc` (RIG Lemma 5.4 steps 1-3 verbatim, including the sentence on why the evaluation is legitimate) -> `lem:rig-inf` (RIG Lemma 5.5 in full, both points over Q = infinity, the values (1.5)) -> `lem:rig-li` (RIG Lemma 5.6 in full: A_n recursion, degree, Li(1/Q), the expansion at t = 0 and why g >= 2 is needed) -> proof of `thm:rig-main` (RIG s.5.5, the four bullets and the extension to the excluded loci) -> `prop:rig-g0` (genus 0: GW_{0,d} = 1/(d^3 epsilon_1 epsilon_2) from `lem:ide-g0`) -> `lem:rig-g1` with `eq:rig-loc` (RIG Lemma 5.7 and "The case d = 1 by direct localisation": the fixed loci, the five factors, V_+-(1), the sum; the [BS] input for d >= 2 stated plainly; the sentence that a direct computation of L_d would remove it, as a remark) -> `cor:rig-conj` (RIG s.5.6: the conjecture stated by content, with `\draftnote{BS label: conj: rigidity}`; what is implied: X the resolved conifold, every d >= 1, G = T_CY and every connected G whose maximal torus acts through a subtorus of T_CY, via contravariance and R_G = R_T^W; the explicit value Sum_g P_g d^{2g-3}/(epsilon_1 epsilon_2) = (1/d)(2 sinh(d epsilon_1/2))^{-1}(2 sinh(d epsilon_2/2))^{-1}; not implied: actions not conjugate into the five-torus; the constant-map term) -> `rmk:rig-BS` (RIG s.6.3, neutral: the subtori alpha_i = epsilon_j by degeneration; the case E = 0 reached by polynomial continuation and the unrefined formula there; one sentence that BS verified the formula on T_CY in low degree and genus by computer; no evaluation) -> `prop:rig-flop` (RIG s.6.1: F_g(infinity) = F_g(0) from `thm:rig-main` and `lem:rig-li`; the limiting phases at s = 1/Q = 0 at both colours, their normalised expansions by `lem:row-betaphase` ((a, b) = (alpha_1 - alpha_2, A) and (-alpha_2, alpha_1)), multiplication by h_nw, the rows diag(h_{wt'_-}, h_{wt'_+}), Delta_i(infinity) = e(wt'_-+), the identification with the constant-map term of the flop Z' by `prop:ide-const` read for Z'; the geometric description of the flop; the resulting identity among constant-map lambda-class integrals, "for which no direct proof is given here". `\draftnote{the analyticity of the W_k(U_i) at s = 0 and the limits of the phases are to be written from Lemma 5.5's computation; check that nothing beyond lem:row-betaphase and lem:rig-inf is needed}`) -> `rmk:rig-F0` -> `lem:rig-3Tid` (RIG s.7 proof (4): the ring map alpha_0 -> 0 then epsilon_2 -> alpha_1 + alpha_2 - epsilon_1 applied to the identity at indeterminate weights; GW_{g,d} for d >= 1 in the ring of `lem:rig-poly`; GW_{g,0} regular at alpha_0 = 0 by the residue cancellation written out; the T_3-fixed locus of the degree-0 space and the localisation on P^1; the graph side regular along alpha_0 = 0 at generic Q with the list of reasons; "h_nw -> h_{(alpha_1, alpha_2, -epsilon_1, -epsilon_2)}") -> `cor:rig-3T` (RIG s.7 statement and proof (1) (= `cor:lea-a0`), (2) (Theorem A holds at alpha_0 = 0 directly, by `rmk:lea-fibrewise` and LCY Remark (iv) of Cor 12.1), (3) (`thm:reg-main` applies verbatim; the vertex class carries (epsilon_1, epsilon_2), the rows (alpha_1, alpha_2), and M_irr needs no relation between them); the Status paragraph's last sentence as a remark: rigidity on T_3 is the restriction of Theorem D by polynomiality; the alpha_0-adic step is not available on T_3 itself).

Instructions.
- RIG's "[P mod ...]" tags become sentences naming the input at the point of use.
- RIG s.0.4 and s.6.3 contain comparisons with BS; keep only the mathematical content, neutrally, without line numbers.
- The letter N for the numerator in `lem:rig-disc` clashes with the critical polynomial N(s) of s02; the critical polynomial does not occur in s07, but say "N here denotes the numerator" or use the subscripted macro of NOTATION.md s.8.
- The local variable s = 1/Q of `lem:rig-inf` and `prop:rig-flop` is `\mathsf{s}`; the Euler variable does not occur in s07.
- Evidence sentences of RIG ("[E: ...]", "Checks s.9", "Numerical confirmation") are dropped; their content goes to the appendix.

Depends on: s02, s04 (`cor:lea-a0`, `rmk:lea-fibrewise`), s05, s06. Used by: s01.

Length: 40-50k characters.

## a_evidence

Covers: LCY s.13.3 (the two scripts' checks, described without names), LCY s.13.1 second referee item (the independent check to t^14, z^23 and the four failing controls, described as an independent computation), RIG s.9 A-F and the "Not re-run here" list, FEAS reports[0].tested, reports[1].tested, crit.tested (the lists are printed in the JSON; read them there), FEAS crit "NEW" items.

Opening sentence (`sec:evi-scope`): nothing in this appendix is a proof step; the proofs of the note are complete without it; the checks are exact rational or F_p arithmetic at the stated orders and weight points; a check that fails is listed when it is a control that must fail.

Subsections and contents.
- `sec:evi-rows`: identities of s02-s03 at random rational points (LCY check1: the amplitude identity, the operator identity, the identities of `lem:set-N`, the three forms of Theta, self-adjointness, R_tr = delta (log alg)_tau, the multiset identity as power sums, first-order values); the rows (LCY check2: KEY FORMULA through tau^8 z^11 at two CY points, r^pt_{0,1}(z; 0), L1 degree test k <= 8, the A_j closed form k <= 7, F2, Row 2, K_1 = 1, E_11 values, Theta, the flow, L3 bidegree, alpha_0 = 0 rescaling); the independent computation through t^14, z^23 at nine weight points (KEY FORMULA, F1-F3, L1 k <= 13, L2 k <= 9, L3 k <= 9, flow, QDE, symplecticity, regular colour, alpha_0 = 0, non-CY points) with the four controls that fail (class (E + 1, -epsilon_1, -epsilon_2); no constant h_nw; A -> alpha_2 in the rows; lam_j = 1/(E + j psi)) and the two sharpness observations (u = 0, odd k; E = 9 alpha_2/7, k = 8).
- `sec:evi-graphs`: FEAS: coloured graphs with an edge regular at g = 2 (20/20), 3 (350/350), 4 (8596/8596), 5 (284806/284806, one sample); vertex graph with only t^{2-2g}; tallies (pole order, bound) at g = 2, 3, 4; special samples (u = 0, epsilon_1 = epsilon_2, alpha_1 = alpha_2, 3 alpha_2 = E, ...); controls: class (E, -alpha_1, -A) and the point frame give singular edge graphs with unchanged totals; frames S1, S3; wrong constants.
- `sec:evi-loc`: direct localisation epsilon_1 epsilon_2 GW_{g,d} = P_g d^{2g-3}: g <= 7 at d = 1; g <= 5 at d <= 3; g <= 4 at d <= 4; E = 0 and alpha_1 = alpha_2 samples; genus 0 and 1; non-CY control fails for g >= 1. Also the referee's d = 1, g <= 2 check (RIG s.5.5 bracket), described as an independent computation.
- `sec:evi-rational`: rational reconstruction in Q at g = 2, 3: totals (1 - Q)^{-2}, (1 - Q)^{-4} with numerator degrees 2, 4; edge graphs with Disc^{-4}, Disc^{-8} poles cancelling in the total; regular at infinity; the point frame at g = 2.
- `sec:evi-rig`: RIG s.9 A-F (S(y) expansion, P_2, P_3, beta values; Disc identities; the quadratic algebra identities; Q = infinity limits in floating point; the small root to O(alpha_0^8); Li_{-n} identities; the power-sum identities of the multisets).
- `sec:evi-flop`: F_g(0) = F_g^flop(0) for g = 2..6 at two CY samples, fails off CY; F_g(0) depends on alpha (g = 2, 3, 4, three samples).
- `sec:evi-offCY`: off T_CY with class (E_alpha, -epsilon_1, -epsilon_2): the row still freezes, Delta_1 t^2 -> -E_alpha^3 epsilon_1 epsilon_2, edge graphs regular at g = 2 but the vertex graph gains t^{-1}; at g = 3 four edge graphs singular; at g = 4 fifty-four; the total fails.
- `sec:evi-limits`: finite genus and orders; F_p zero tests probabilistic; what was not tested (the stability argument (Fl6), the descent, the analytic lemma).

No file names, no script names, no referee names, no "report 1/2", no "critic". Write "an independent computation" or "a second implementation".

Length: 12-18k characters.

## Dropped from all sections

Rounds, referees and verdicts, tags, status tables (LCY s.13.1-13.2, RIG s.8), provenance, the dictionary to earlier reports, "proof of record", file and script names, code names and line numbers, "the brief", "the critic", "the feasibility reports", mention of AI, dates of writing, errata of BS, BS line numbers, every "[Evidence: ...]" bracket in the body (content to the appendix), every "[P]/[C]/[E]" tag, "Not used" paragraphs except where they prevent a misreading (one sentence), LSP s.12, s.14, Appendices A-B, LAP beyond the statements, LP2 beyond s.7.1 and s.7.6 statements.

## Flags for all writers

- Macro clashes to avoid (see NOTATION.md): A, B (weights) are `\wA`, `\wB`, never bare A, B (B_{2k} are Bernoulli numbers; `\Aop` is the QDE matrix; `\Rinv` is R^{-1}(psi)). Gauss's parameters are `\ga`, `\gb`, `\gc`, never alpha, beta, gamma. `\tau` is 1 - 1/Q only; the scale-1 parameter of [LP2] is `\tauw`. `\Lsp` is the leaf span, `\Lnorm` the Laplace expansion, `\Nlap` the expansion N(z, tau), `\Ns` the critical polynomial, `\Nb` the normal bundle. `\Rm` is the R-matrix; the t-coefficients of zG are `\Rk_k`. `u` (sans, `\uq`) is the quadratic (E - 2alpha_1)(E - 2alpha_2); canonical coordinates are italic u_i. `E` is the weight; the Hodge bundle is `\Hdg`; edge factors are `\Edge`; `\Ephi` is Exp[-dphi/z]. `\Pg` is P_g, `\Psh_k` the sharpness product, `\Pq` the critical polynomial in U. `\Dop_i` and `\Dd` both print as sans D: the first always carries an index i, the second never.
- Renames relative to the sources: LCY's A, B -> `\wA`, `\wB`; u -> `\uq`; N(z, tau) -> `\Nlap`; Abeta -> `\Abeta`; alg -> `\alg`; dphi -> `\dphi`; Ephi -> `\Ephi`; X_n -> `\Ephi_n`; Psi_n -> `\Kcoef_n`; R_k -> `\Rk_k`; Lnorm -> `\Lnorm`; l -> `\lrem`; the QDE matrix A -> `\Aop`; pi -> `\pim`; P_k -> `\Psh_k`; C_j -> `\Cj_j`; e_i (idempotents) -> `\eid_i`; v_i -> `\vv_i`; r^pt -> `\rpt`; r' (frame B) -> `\rr`; Omega_w -> `\Hg_{\bw}`; R (multiset) -> `\msM`; D(R) -> `\dg`; Fpot -> `\Fpot`; Dom -> `\Dom`; T (dilaton) -> `\Tdil`; RIG's Kw -> `\KCY`, K -> `\Kred`, P(U) -> `\Pq`, G -> `\sG_g`, C -> `\Ccov`, y -> `\ydc`, A_n -> `\Aeul_n`, S(y) -> `\Ssh`.
- Where LSP steps are transcribed with substitutions: s03 (Laplace calculus verbatim, no weights; Props 3.10, 3.12, Lemma 3.13, Thm 3.14, Prop 3.15, Cor 3.16, Thm 3.17 with the dictionary), s04 (the list in the s04 instructions), s05 (Prop 13.2 (a)-(c), (e) and Lemma 13.3 at indeterminate weights; no substitution, but the hypotheses listed for T_CY), s06 (s.13.0, s.13.1 with c -> E, -x c^3 -> -epsilon_1 epsilon_2 E^3, Q(x) -> K_CY).
- Length expectations: s01 15-20k, s02 30-40k, s03 55-70k, s04 100-125k, s05 45-60k, s06 25-35k, s07 40-50k, appendix 12-18k characters. Equations in `equation`/`align` environments with labels only where referred to; displays without labels otherwise.
- Each section file begins with the comment line `% sNN_name. Section file; new macros at the top; label prefix: xxx.` and its `\section` and `\label{sec:xxx}` exactly as in the stub.
- Compile check: copy note/ to the scratch folder and run pdflatex there (never inside the Drive folder); fix every undefined reference to a label of LABELS.md before handing the section back; a reference to a label of another section that is not yet written is allowed and reported.
