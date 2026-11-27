# Label registry

Every named result of the note, its unique `\label`, the section file that proves it and its source location. LCY = cy_torus/LC_CY_PROOF.md, RIG = cy_torus/RIGIDITY_PROOF.md, LSP = leaf_proof/leaf_structure_proof.md, LAP = lemmaA/LEMMA_A_PROOF.md, LP2 = lp2_explore/LP2_REGULARITY_PROOF_draft.md, CP = companion LaTeX paper lp2_explore/paper/sections/*.tex (reusable LaTeX), FEAS = cy_torus_explore/feasibility_result.json. Use these labels exactly. A writer adds further labels (equations, auxiliary lemmas) only with the prefix of its own section file, and refers to another section's results only through labels listed here. The "note name" is the name to use in the text if the statement needs one; the source tag in the last column never appears in the note.

Prefixes: s01 `int`, s02 `set`, s03 `row`, s04 `lea`, s05 `ide`, s06 `reg`, s07 `rig`, appendix `evi`. Form: `thm:row-closed`, `eq:lea-key`, `rmk:set-d0`. Section labels (already in the stubs): `sec:int`, `sec:set`, `sec:row`, `sec:lea`, `sec:ide`, `sec:reg`, `sec:rig`, `sec:evi`.

Named parts inside a statement keep the source's tags: (F1)-(F3), (L1)-(L3) in the leaf theorem; (a)-(e) in the identification; (i)-(v) in the regularity theorem; (G0)-(G2), (D1)-(D5), (P1)-(P3) in the Laplace calculus; (LAM1)-(LAM6), (PF), (T), (P) in the conventions; (H0), (S), (C), (L0) on the Hodge side; (F1)-(F9) in the flatness proof (these clash in name with the leaf theorem's (F1)-(F3); inside s05 write them (Fl1)-(Fl9) and say so once).

## s01_intro (statements only; each proved later)

| Note name | Label | Restates |
|---|---|---|
| Theorem A (leaf structure at generic weights) | `thm:int-leaf` | `thm:lea-main` |
| Theorem B (the graph sum is the potential) | `thm:int-ident` | `prop:ide-main` |
| Theorem C (regularity and the polar part) | `thm:int-reg` | `cor:reg-main` |
| Theorem D (rigidity) | `thm:int-rig` | `thm:rig-main` |
| Corollary (the conjecture of [BS] for the resolved conifold) | `cor:int-conj` | `cor:rig-conj` |
| Subsections | `sec:int-setting`, `sec:int-statements`, `sec:int-inputs`, `sec:int-notclaimed`, `sec:int-strategy`, `sec:int-plan` | |

## s02_setup

| Note name | Label | Source |
|---|---|---|
| The target, the five-torus, the tangent weights | `eq:set-wts` | RIG s.0.1, s.1.1; LSP Prop 13.2 statement; BS App. A (comparison only) |
| The Calabi-Yau torus; non-degeneracy (ND) | `def:set-ND` | LCY s.0.2, s.1.1; RIG s.1.1 |
| Symmetries of the weights (exchange of p_+ and p_-, alpha_1 <-> alpha_2, epsilon_1 <-> epsilon_2) | `lem:set-sym` | RIG s.1.1 |
| Invariants by localisation; the series F_g with d >= 0 | `def:set-inv`, `rmk:set-d0` | RIG s.0.1; BS defn (by content) |
| Fields K_alpha, K_CY, K | `eq:set-fields` | LCY s.1.1; RIG s.0.1, s.1.1 |
| The superpotential and its derivatives; the critical equation; the roots and their limits | `lem:set-crit` | LCY (C1'), (C2'), Lemma 4.2_CY(a), (b); RIG s.1.3 incl. (1.3) |
| Frobenius data in the flat basis: eta, the product, L_1 + L_2, L_1 L_2, commutativity | `prop:set-frob` | LCY Lemma 4.1_CY (proof from the genus-0 potential); RIG s.1.2 |
| Idempotents and Delta_i: the four forms of (4.1), the product identity, 1/Delta_1 + 1/Delta_2 | `prop:set-idem` | LCY Prop 4.3_CY(a)-(c) and the algebraic identities of (d); RIG (1.1), (1.4) |
| The quadratic algebra and its involution | `def:set-QA` | LCY (C3'), Lemma 8.1_CY (first sentence, without code) |
| Algebraic functions of tau: delta, s^*, s~, the discs, U_1 = 1 - tau s^*, U_2 = 1 - tau s~ | `def:set-alg` | LCY s.1.3 |
| The critical polynomial N(s): roots, the identities (i)-(iv), ds^*/dtau, the envelope formula | `lem:set-N` | LCY Lemma 3.1_CY(a), (b) (the parts (c), (d) are `lem:row-hess`) |
| Eigenvalues in terms of s^*: L_1 = A/(tau s^*) - alpha_0, L_2 = B s^* - alpha_0, sum, product, difference delta/tau, the expansion L_1 = E/t + l(t), l(0) | `lem:set-eig` | LCY Lemma 4.2_CY(c), (d) (part (e) is `lem:row-canon`) |
| Delta_2 = epsilon_1 epsilon_2 B delta s^* L_2; Delta_1 Delta_2; Delta_2(0), Delta_2'(0) = 0 | `lem:set-Delta2` | LCY Prop 4.3_CY(d) second and third formulas, (e) for Delta_2 |
| The discriminant Disc(Q): (i)-(v) | `lem:set-disc` | RIG s.1.4 |
| Behaviour at Q = 1, Q = 0, Q = infinity | `lem:set-limits` | RIG s.1.5 |
| Point-frame row (Wick polynomials), frame-B rows, H-row, the covector r_i | `def:set-rows` | LCY (C7')-(C9'); LSP (C7), (C9) for the Wick formula and the H-row |
| The multiset identity h_nw h_{(E,-alpha_1,-A)} = h_w and uniqueness of the class | `lem:set-hmult` | LCY s.3.9 "The constant line", (3.9a) |
| Dilaton factors, q_t, the pairing P_ij, edge factors, legs | `def:set-factors` | LCY (C10')-(C12') |
| Edge factors as generating functions | `lem:set-edgegen` | LSP Lemma 1.1 (weight-free) |
| rho, the leaves, Lam_m | `def:set-rho` | LCY (C13') |
| Scope: the formal series | `rmk:set-scope` | LCY (C14') |
| Coloured stable graphs, Cont(Gamma, c), the graph sum G_g | `def:set-graphsum`, `eq:set-graphsum` | RIG s.3 (3.1); LSP s.0.1; LSP (C10) for the dilaton-shifted form |
| Rings O, M, C((z)), R_z, iota, refl, TRANSFER, the passage (P) | `lem:set-rings` | LSP s.1.4 with LCY s.1.4 (denominators at tau = 0) |
| Leaf spaces; (LAM1)-(LAM6); (PF); (T) | `lem:set-leafid` | LSP s.1.4 with c -> E |
| Hodge side: Hg_w, R, base-point derivatives, Dom, F^pot_h on Dom, the dilaton-shifted form and the overlap, Xi_h, dg, Z[h,n], M_irr(h) | `def:set-fpot`, `def:set-mirr` | LSP s.13.0 (definitions only; (H0) is `lem:reg-H0`); LAP s.1; CP s02_setup.tex def:set-fpot |

## s03_rows

| Note name | Label | Source |
|---|---|---|
| The Euler coordinate (3.1), (3.2) | `lem:row-euler` | LCY s.3.1 |
| Hessians (3.3), alg^2 (3.4), (3.4'), the tensor rule (3.5) | `lem:row-hess` | LCY Lemma 3.1_CY(c), (d); LSP Lemma 3.1(d) |
| The rescaled derivatives w_k = tau^k W_k(U_1): (3.6), (3.7), (3.8) | `lem:row-wr` | LCY s.3.3, Lemma 3.2_CY |
| The ring A and the Gaussian functional; (G0)-(G2) | `def:row-Acal`, `lem:row-gauss` | LSP s.3.4, Lemma 3.3 |
| Normalised formal Laplace expansions; (D1)-(D3), (D4), (D5) | `def:row-Lnorm` | LSP Def 3.4; LCY s.3.4 (D4), (D5) |
| Total derivatives (P3) | `prop:row-P3` | LSP Prop 3.5 |
| Parameter derivative with moving expansion point (P2) | `prop:row-P2` | LSP Prop 3.6 |
| Coordinate invariance (P1) | `prop:row-P1` | LSP Prop 3.7 |
| Laplace's method | `lem:row-laplace` | LSP Lemma 3.8 (proved in full there) |
| Homogeneity of the Wick coefficients | `lem:row-homog` | LSP Lemma 3.9 |
| The point-frame row is the normalised expansion of the Euler integrand | `prop:row-wickrow` | LCY Prop 3.10_CY with the proof of LSP Prop 3.10 (steps (1)-(2)); LSP Remark 3.11 |
| The amplitude identity (3.16) | `lem:row-amp` | LCY s.3.6 (proof written out there) |
| Gauss's equation for the formal expansion | `prop:row-gaussode` | LCY Prop 3.12_CY with the proof of LSP Prop 3.12 |
| Uniqueness | `lem:row-unique` | LCY Lemma 3.13_CY with the recursion of LSP Lemma 3.13 |
| The closed form (3.17) | `thm:row-closed` | LCY Thm 3.14_CY, LSP Thm 3.14 |
| The Beta-phase lemma (normalised expansion of a log U + b log(1 - U) at its critical point equals h_{(-a,-b,a+b)}) | `lem:row-betaphase` | LSP Lemma 13.3, the paragraph "This limit is h_{(-a,-b,a+b)}" (Laplace + Stirling); CP s11 lem:con-norm proof |
| The constant A_beta(z) = h_{(E,-alpha_1,-A)} E/(E + z) | `prop:row-beta` | LCY Prop 3.15_CY; LSP Prop 3.15 (proof transcribed at (alpha_1, A); alternatively the rescaling (D4)); LCY's z^1 hand check |
| Regularity of the tau-coefficients | `cor:row-reg` | LCY Cor 3.16_CY, LSP Cor 3.16 |
| The frame-B singular row (3.18) | `thm:row-main` | LCY Thm 3.17_CY, LSP Thm 3.17 (coefficientwise meaning) |
| (F2) at t = 0 | `cor:row-F2` | LCY Cor 3.18_CY |
| Eikonal and transport decomposition (3.19); Eik(dphi_tau) = 0; P_tr = delta; R_tr | `lem:row-eik` | LCY s.3.10 |
| The transport equations; Phi_0 = 1/alg; Ephi F = Phi | `prop:row-transport` | LCY Prop 3.19_CY (the identity R_tr = delta (log alg)_tau proved by hand there) |
| First orders: F_1, X_1, H_1, [tau^1] alg | `rmk:row-first` | LCY s.3.11 |
| Euler's integral (motivation, used nowhere) | `rmk:row-eulerint` | LCY Remark 3.20 (optional) |
| Canonical coordinates u_i, L_1 = mu + dphi', L_2 = mu + alpha_2 - L_1 | `lem:row-canon` | LCY Lemma 4.2_CY(e) |

## s04_leaves

| Note name | Label | Source |
|---|---|---|
| Theorem A: leaf structure at generic weights, (F1)-(F3), (L1)-(L3), sharpness, large-z form | `thm:lea-main` | LCY Thm LC_CY (s.0.2) |
| Where the Calabi-Yau condition is not used | `rmk:lea-cy` | LCY Remark 0.1 |
| Quantum differential equation and symplecticity at both idempotents | `thm:lea-qde` | LCY Thm Q_CY (s.0.3) |
| Closed forms | `thm:lea-closed` | LCY Thm CF_CY (s.0.3, s.11) |
| The Hessian in the Euler coordinate, rho, the KEY FORMULA | `prop:lea-key`, `eq:lea-key` | LCY Prop 4.4_CY; LSP Prop 4.4 (a)-(c) |
| The operator identity (4.4) | `lem:lea-op` | LCY Lemma 4.5_CY |
| Row 2 of the quantum differential equation at the singular idempotent | `prop:lea-row2` | LCY Prop 4.6_CY; LSP Prop 4.6 |
| The operator A, its eigen-covectors, self-adjointness, (4.8) | `prop:lea-A` | LCY Prop 4.7_CY |
| The symplectic constant K_1 = 1; sigma; (4.9) | `prop:lea-K1` | LCY Prop 4.8_CY; LSP Prop 4.8 |
| K_1 = 1 as a Wronskian identity | `rmk:lea-wronsk` | LCY Remark 4.9 |
| The tau-coefficients F_n, X_m, H_n: poles and values at infinity | `lem:lea-taucoef` | LCY Lemma 5.1_CY; LSP Lemma 5.1(b), (c) |
| The t-coefficients of zG: (5.1), (5.2) | `lem:lea-tcoef` | LCY Lemma 5.2_CY; LSP Lemma 5.2 |
| No pole at z = 0 | `rmk:lea-nopole` | LCY s.5.3 |
| Partial fractions (5.3) | `prop:lea-pf` | LCY Prop 5.3_CY; LSP Prop 5.3 |
| (L1) for all k; g_k; the sum rule | `thm:lea-L1` | LCY Thm 5.4_CY with remarks (i)-(iv); LSP Thm 5.4 |
| The H-row in leaf form; sigma^{(k)}; (6.3), (6.4); k = 0 by hand | `prop:lea-sigma` | LCY Prop 6.1_CY; LSP Prop 6.1 |
| The projection on the regular eigen-covector | `prop:lea-proj` | LCY Prop 6.2_CY |
| E_11 in terms of sigma: (7.1), symmetry | `eq:lea-Ebar` | LCY s.7.1; LSP s.7.1 |
| The equation for Y | `lem:lea-Yqde` | LCY Lemma 7.1_CY |
| The flow (7.2), (7.3); expansion of Theta | `thm:lea-flow` | LCY Thm 7.2_CY; LSP Thm 7.2 (i), (ii) |
| Integration from t = 0; (F3) for E_11 | `eq:lea-int` | LCY s.7.4; LSP s.7.4 |
| (L3) for all k; (7.5), (7.6); rank one at t = 0 | `thm:lea-L3` | LCY Thm 7.3_CY; LSP Thm 7.3 |
| Powers of rho; (L2) for j = 1; (7.7) | `cor:lea-rhoe` | LCY Cor 7.4_CY; LSP Cor 7.4 |
| The involution exchanges the idempotents; the same scalar at both | `lem:lea-conj` | LCY Lemma 8.1_CY (rephrased without code) |
| Analyticity at the regular idempotent; (F1) for Delta_2; (F3) for Dil_2, E_22 | `lem:lea-reg` | LCY Lemma 8.2_CY |
| Data at the regular critical point (8.1) | `lem:lea-regdata` | LCY Lemma 8.3_CY |
| Row 2 at the regular idempotent (8.2) | `thm:lea-row2reg` | LCY Thm 8.4_CY with LSP Thm 8.4 Steps A-D |
| Orthogonality | `cor:lea-orth` | LCY Cor 8.5_CY; LSP Cor 8.5 |
| The regular rows at t = 0: (a) analytic input via `lem:row-betaphase`, (b), (c) | `prop:lea-regrows` | LCY Prop 8.6_CY |
| K_2 = 1 | `cor:lea-K2` | LCY Cor 8.7_CY |
| E_22 at t = 0 is not of rank one; the regular column is determined by the singular data | `rmk:lea-E22` | LCY after Prop 8.6_CY; LSP Remarks 8.8-8.10 |
| Extraction of a leg | `lem:lea-leg` | LSP Lemma 9.1 (weight-free) |
| (L2) for j = 2 and (F3) for E_12 | `prop:lea-L2` | LCY Prop 9.2_CY; LSP Prop 9.2 |
| Leaf form of the E_12 edge (9.6), (9.7) | `prop:lea-E12` | LCY Prop 9.3_CY; LSP Prop 9.3 |
| E_21 = E_12^T; regularity of E_21 | `prop:lea-transp` | LCY Prop 9.4_CY; LSP Prop 9.4 |
| All legs; unrescaled ends | `cor:lea-alllegs` | LCY Cor 9.5_CY; LSP Cor 9.5 |
| How the factors attach to a singular vertex | `rmk:lea-use` | LCY s.9.6; LSP s.9.6 |
| Proof of Theorem A and of `thm:lea-qde`; the table of (F3) with values at t = 0 | `eq:lea-table` | LCY s.10 |
| Where each input is used; where (ND) is used | `rmk:lea-inputs` | LCY s.10 last paragraph |
| Dilaton leaf coefficients A_j(t) | `thm:lea-Aj` | LCY Thm 11.1_CY |
| Sum rule, Theta_k, S_j, dilaton consistency | `cor:lea-sumrule` | LCY Cor 11.2_CY(a)-(e) |
| Top and subleading coefficients; sharpness | `cor:lea-top` | LCY Cor 11.3_CY |
| The loops in closed form | `thm:lea-loops` | LCY Thm 11.4_CY; LSP Thm 11.4 |
| Anti-diagonal and corner | `cor:lea-anti` | LCY Cor 11.5_CY |
| The legs: top coefficients, E_12 at t = 0 | `rmk:lea-legs` | LCY s.11.3 |
| Low orders by hand | `rmk:lea-low` | LCY s.11.4 |
| alpha_0 = 0: rescaling to the fibrewise rows | `cor:lea-a0` | LCY Cor 12.1 with Remarks (i)-(iv) |
| The fibrewise theorem as the special case | `rmk:lea-fibrewise` | LCY s.0 header dictionary; LSP Thm LC (statement only, as the case alpha_0 = 0, alpha_i = epsilon_i) |

## s05_identification

| Note name | Label | Source |
|---|---|---|
| Theorem B: the frame-B formula on T_CY | `prop:ide-main`, `eq:ide-main` | RIG s.3 statement; LSP Prop 13.2 statement |
| Fixed loci, the twisted theory, CohFT axioms, the pairing | `lem:ide-fixed` | LSP Prop 13.2(a) first paragraph; CP s04_generic.tex lem:gen-fixed, lem:gen-axioms (rank-free parts) |
| The genus-0 potential: classical part, Aspinwall-Morrison term, H o H | `lem:ide-g0` | LSP Prop 13.2(a), (b)(i); LCY Lemma 4.1_CY proof |
| Semisimplicity at Q = 0 on T_CY | `lem:ide-ss` | RIG s.3 (a); RIG s.1.5 |
| Existence of R; field of definition; descent; freeness | `prop:ide-teleman`, `lem:ide-free` | LSP Prop 13.2(a) "Field of definition", "Descent"; CP prop:gen-teleman, lem:gen-free |
| Ambiguity of R | `lem:ide-ambig` | LSP Prop 13.2(a) (the X argument); CP lem:gen-ambig |
| Flatness: (Fl1)-(Fl9) | `prop:ide-flat` | LSP Prop 13.2(a) "Flatness" (F1)-(F9); CP prop:gen-flat for the full text of each step, specialised to rank 2 |
| The route through Teleman's Proposition 7.7 (not used) | `rmk:ide-tel77` | LSP Prop 13.2(a) final remark (optional, neutral) |
| Normalisation constants at Q = 0: r^pt_{0,+-}(z;0) and h_nw | `lem:ide-norm` | LSP Lemma 13.3; CP s11 lem:con-norm; uses `lem:row-betaphase` |
| The rows: (i) Frobenius data (13.2), (ii) Row 2 via Euler's identity (13.3), (iii) Q = 0, (iv) symplecticity | `prop:ide-rows`, `eq:ide-PF` | LSP Prop 13.2(b); CP prop:con-frame step (b) |
| Hodge classes as an R-action: Hg_v = h_v.1 | `lem:ide-mumford` | LSP Prop 13.2(c) first sentences; CP lem:gen-mumford |
| Non-vanishing of ch_{2k-1} on M_{k,1} | `lem:ide-ch` | LSP Prop 13.2(c); CP lem:gen-ch |
| The constants: R = R_Z | `prop:ide-const` | LSP Prop 13.2(c) |
| Frame B: R' = R_Z h_w^{-1} | `prop:ide-frameB` | LSP Prop 13.2(e); CP prop:con-frame step (e) |
| The coloured terms as elements of K_CY(C); regularity off {0, 1, infinity, Disc = 0}; the colour sum in K_CY(Q); regularity at Q = 0 | `lem:ide-terms` | RIG s.3 "The description of the coloured terms" |
| The alternative route by specialisation epsilon_2 -> E - epsilon_1 | `rmk:ide-alt` | RIG s.3 parenthesis |
| What is and is not an analytic continuation; the T_CY-invariants as restrictions | `rmk:ide-scope` | RIG s.5.2 paragraph after Lemma 5.1; RIG s.3 last sentence of the list |
| Unitarity and t-independence of the frame-B rows; the dilaton shift | `rmk:ide-frame` | LSP Remarks (ii), (v) after Prop 13.2; CP prop:con-frame remarks (1), (4) |

## s06_regularity

| Note name | Label | Source |
|---|---|---|
| Dilaton homogeneity (H0) on Dom | `lem:reg-H0`, `eq:reg-H0` | LSP s.13.0 (H0) with proof; LP2 s.7.1 (H0); CP s10 lem:reg-H0 |
| Range: every g_v <= g - 1 in a graph with an edge | `lem:reg-range` | LSP Lemma 13.0; CP lem:reg-range |
| String identity (S) | `lem:reg-string` | LP2 s.7.1 (S_P2); CP lem:reg-string |
| Vanishing on 1^p (C) | `lem:reg-C` | LP2 s.7.1 (C_P2); CP lem:reg-C |
| The lam_0 identity (L0) (remark; not used) | `rmk:reg-L0` | LP2 s.7.1 (L0_P2), s.7.2; CP lem:reg-L0 |
| Theorem L in the convention of this note (cited) | `thm:reg-L` | LAP Thm L, Cor L1, Cor L3, s.7.5; [LP2, Theorem 8.1, Corollaries 8.2-8.4] |
| Theorem VI (cited) | `thm:reg-VI` | LP2 s.7.6; [LP2, Theorem 9.1] |
| beta_{2g-2} = (2g-3)! P_g/(epsilon_1 epsilon_2); hand values at g = 2, 3 | `lem:reg-beta`, `eq:reg-beta` | RIG (0.1), s.4 hand check |
| Theorem R on T_CY: (a) with M_irr(h), h <= g - 1; (b) with M_irr(g) and Theorem VI | `thm:reg-main` | RIG s.4 Steps 0-5; LSP s.13.1 Theorem R (i)-(v); CP s11 cor:con-reg proof |
| Theorem C: Corollary R_CY (a)-(c) | `cor:reg-main` | RIG s.0.3, s.4 Step 6 |
| Specialisation to weight points; the value gap_g | `rmk:reg-special` | RIG s.4 "Specialisation" |
| Which instances of M_irr are used, by genus | `rmk:reg-genus` | LSP s.13.2 "By genus"; RIG s.4 Step 4 |
| Why a constant leaf would be fatal and does not occur | `rmk:reg-nol0` | LP2 s.7.2 (generic form: `prop:lea-sigma`, `thm:lea-L3`, `prop:lea-L2` have no constants) |

## s07_rigidity

| Note name | Label | Source |
|---|---|---|
| Theorem D: rigidity for g >= 2 | `thm:rig-main`, `eq:rig-main` | RIG s.0.3, s.5 |
| Reduction and plan | `sec:rig-plan` | RIG s.5.1 |
| Possible poles of F_g | `lem:rig-poles` | RIG Lemma 5.1 |
| Polynomiality of epsilon_1 epsilon_2 GW_{g,d}, d >= 1; restriction to T_CY | `lem:rig-poly` | RIG Lemma 5.2; CP prop:gen-poly for the LaTeX of the argument; BS by content (shape lemma, contravariance) |
| The constant-map term GW_{g,0} | `rmk:rig-d0` | RIG remark after Lemma 5.2 |
| The small root Q_1 | `lem:rig-root` | RIG Lemma 5.3 |
| No pole on Disc = 0 | `lem:rig-disc` | RIG Lemma 5.4 (steps 1-3) |
| Regularity at Q = infinity | `lem:rig-inf` | RIG Lemma 5.5 |
| The polylogarithms Li_{-n} | `lem:rig-li` | RIG Lemma 5.6 |
| Proof of Theorem D; extension to the excluded loci | end of `sec:rig-proof` | RIG s.5.5 |
| Genus 0 | `prop:rig-g0` | RIG s.5.5 (g = 0); LSP Prop 13.2(a) (Aspinwall-Morrison) |
| Genus 1: epsilon_1 epsilon_2 GW_{1,d} = P_1/d | `lem:rig-g1` | RIG Lemma 5.7 (d >= 2 with the formula of [BS] on alpha_i = epsilon_j; d = 1 by localisation) |
| Degree-1 localisation formula (L) | `eq:rig-loc` | RIG s.5.5 "The case d = 1" |
| What is implied for the conjecture of [BS] | `cor:rig-conj` | RIG s.5.6 |
| Relation to the results of [BS] (neutral) | `rmk:rig-BS` | RIG s.6.3 |
| The flop: F_g(infinity) = F_g(0) and its identification with the constant-map term of the flop | `prop:rig-flop` | RIG s.6.1 (limiting rows via `lem:row-betaphase`, written out) |
| The constant-map term depends on alpha (not determined here) | `rmk:rig-F0` | RIG s.6.2 (statement only; evidence in appendix) |
| The three-torus: identification on T_3 by the ring map alpha_0 -> 0 | `lem:rig-3Tid` | RIG s.7 proof (4) |
| Corollary 3T: regularity and the polar part on T_3; rigidity there by restriction | `cor:rig-3T` | RIG s.7 statement, proof (1)-(3) (with (2) now `thm:lea-main` at alpha_0 = 0), Status paragraph |

## a_evidence

Subsections: `sec:evi-scope` (opening: nothing here is a proof step; arithmetic exact unless stated), `sec:evi-rows` (LCY s.13.3: identities, rows, KEY FORMULA, L1, L3, flow, K_1, alpha_0 = 0 rescaling; the independent check through t^14, z^23 and the four failing controls, LCY s.13.1 second referee item, described as checks), `sec:evi-graphs` (FEAS reports[0], reports[1] "tested": coloured graphs with an edge regular at g = 2..5; vertex graph; tallies; frame and class controls that fail), `sec:evi-loc` (FEAS crit: direct localisation through g = 7 at d = 1), `sec:evi-rational` (FEAS crit: rational reconstruction in Q at g = 2, 3; Disc poles cancel), `sec:evi-rig` (RIG s.9 A-F), `sec:evi-flop` (flop identity g = 2..6; alpha-dependence of F_g(0)), `sec:evi-offCY` (controls off T_CY: polar terms at g = 2 vertex, g = 3 edges), `sec:evi-limits` (what the checks do not cover).

## Pinned numbers for the companion paper [LP2] (draft of 2026-11-02; mark each first use with \draftnote{LP2 numbering provisional})

| Companion label | Number | Content |
|---|---|---|
| thm:van-main | Theorem 8.1 | Theorem L at scale 1 |
| cor:van-spec | Corollary 8.2 | specialisation to any tau_0 with sigma != 0 |
| cor:van-parts | Corollary 8.3 | parts >= 1 |
| cor:van-scale | Corollary 8.4 | general scale, eq. (8.?) `eq:van-scale` |
| thm:val-main | Theorem 9.1 | Theorem VI: (-eps1 eps2 E^3)^{g-1} Z[g,0] = beta_{2g-2} |
| thm:val-tau | Theorem 9.3 | Theorem VI at scale 1 |
| lem:reg-H0 | Lemma 10.1 | (H0) |
| lem:reg-string | Lemma 10.2 | (S) |
| lem:reg-C | Lemma 10.3 | (C) |
| lem:reg-L0 | Lemma 10.4 | (L0) |
| lem:reg-range | Lemma 10.6 | range |
| thm:reg-main, cor:reg-sv | Theorem 10.8, Corollary 10.9 | regularity for local P2 (comparison only) |
| prop:gen-poly | Proposition 4.5 | polynomiality by Graber-Pandharipande (comparison only) |
| lem:gen-mumford, lem:gen-ch | Lemmas 4.13, 4.14 | Mumford; ch_{2k-1} != 0 (comparison only; the note proves its own) |
| lem:gen-free, prop:gen-teleman, prop:gen-flat, lem:gen-ambig | Lemma 4.15, Proposition 4.16, Proposition 4.23, Lemma 4.29 | rank-3 versions (comparison only; the note proves its own in rank 2) |
| thm:con-leaf, cor:con-reg, prop:con-frame, lem:con-norm | Theorem 11.1, Corollary 11.4, Proposition 11.7, Lemma 11.8 | the fibrewise conifold section (comparison only) |

Only Theorems 8.1, 9.1 and Corollaries 8.2-8.4 are inputs; everything else of [LP2] is "compare" or is re-proved here.

## Bibliography keys (refs.bib)

| Key | Work | Pinned numbers (form: `\cite[Theorem~(2)]{Tel12}`) |
|---|---|---|
| BS | the paper whose conjecture is proved here | by statement, not number; `\draftnote{BS label: ...}` |
| LP2 | companion paper (placeholder entry) | Theorem 8.1, Corollaries 8.2-8.4, Theorem 9.1 |
| Tel12 | Teleman | Theorem (2) of the Introduction; (1.15); Prop 8.7(ii); 2.17, 2.18, 2.19; (2.15). In remarks only: Def 7.1, Prop 7.7, (7.7.a), Ex 7.17, Prop 8.2(ii), (8.6), Prop 8.7(i), Prop 8.4 |
| PPZ15 | Pandharipande-Pixton-Zvonkine | Defs 2.2, 2.13; Prop 2.14 |
| GP99 | Graber-Pandharipande | the localisation formula |
| Mum83 | Mumford | c(E) c(E^vee) = 1; Mumford's formula |
| FP00 | Faber-Pandharipande | s.1.1 (Mumford's formula), Theorem 4 (lambda_g lambda_{g-1} lambda_{g-2}), p. 3 (Aspinwall-Morrison) |
| LLZ03 | Liu-Liu-Zhou | eq. (1)-(3): only through [LP2] |
| OP04 | Okounkov-Pandharipande | Theorem 1 (remark only, if at all) |
| CG07 | Coates-Givental | App. 1 (iv), (v), (vi) |
| Giv01 | Givental | (1), (6), (7); s.1.3 (3)-(5); s.2.2 Theorem 1 (remark); s.2.3 Lemma |
| Loo96, MW07, Har85, Iva90 | stability | as reported by Teleman 2.18, 2.19 |
| DLMF | | 5.11.1 and s.5.11(ii) (Stirling with remainder); 5.12.1 (Beta integral) |
| Liu13 | Liu | background on localisation, introduction only |
| HK10, CKK14 | background in the introduction | none pinned |

## Labels added while writing the sections

(Writers append here: label, file, one line.)

Appended by the assembler. Every `\label` of the assembled note not listed in the tables above, by file, with the environment or section it belongs to. Line numbers refer to the assembled version. New labels keep the section prefixes; no label is defined twice.

### s01_intro.tex

| Label | Line | Belongs to |
|---|---|---|
| `sec:int-BS` | 148 | section: Relation to the results of \cite{BS} |

### s02_setup.tex

| Label | Line | Belongs to |
|---|---|---|
| `sec:set-target` | 10 | section: The target, the tori and the weights |
| `sec:set-inv` | 60 | section: The invariants and the series |
| `sec:set-crit` | 95 | section: The superpotential and the critical equation |
| `eq:set-P` | 112 | equation: Pq(U)\coloneqq\wB\cQ U^2-\big(\alpha_0+\cQ(\alpha_1-\alpha_2)\big)U- |
| `eq:set-Wk` | 127 | equation: thU^k\log(1-aU)=-\Pol_{k-1}(p_a)\quad(k\geq2),\qquad \cW_k=-\alpha_1 |
| `sec:set-frob` | 155 | section: Frobenius data |
| `eq:set-pairing` | 181 | equation: langle a,b\rangle\coloneqq\sum\eta^{a'b'}a_{a'}b_{b'}=\epsilon_1\eps |
| `sec:set-idem` | 213 | section: Idempotents, the vertex weights, the quadratic algebra |
| `eq:set-Delta` | 230 | equation: begin{split} \Delta_i&=\frac{\alpha_1\alpha_2\wA\wB\epsilon_1\epsilo |
| `sec:set-alg` | 260 | section: Algebraic functions of $\tau$; the critical polynomial |
| `eq:set-envelope` | 295 | equation: dphi_\tau=\frac{\alpha_2\sast}{1-\tau\sast}=\frac{\wA-\wB\sast}{1-\t |
| `sec:set-disc` | 372 | section: The discriminant and the double cover; behaviour at $\cQ=1,0,\infty$ |
| `eq:set-Dinf` | 413 | equation: begin{gathered} \Delta_1(\infty)=\epsilon_1\epsilon_2(\alpha_1-\alph |
| `sec:set-rows` | 435 | section: The rows |
| `eq:set-h` | 439 | equation: h_{\mathsf{v}}(z)\coloneqq\exp\Big[-\sum_{k\geq1}\bk_k\,z^{2k-1}\sum_ |
| `eq:set-rpt` | 453 | equation in defn (Rows) |
| `eq:set-rowB` | 468 | equation: rr_{0,i}(z;t)\coloneqq h_{(\alpha_1,\wA,-E)}(z)\,\rpt_{0,i}(z;t)=\hw |
| `eq:set-Hrow` | 474 | equation: rr_{1,i}(z)\coloneqq(\Lv_i+z\Dop_i)\,\rr_{0,i}(z),\qquad \Dop_i\colo |
| `eq:set-hmult` | 484 | equation in lem (The multiset identity) |
| `sec:set-factors` | 498 | section: The factors |
| `eq:set-Pij` | 511 | equation: Pij_{ij}(z,w)\coloneqq-\langle\vv_i(-z),\vv_j(-w)\rangle=\sum_{p,q}\ |
| `eq:set-edge` | 517 | equation: Edge_{ij}(z,w)\coloneqq\frac{\Pij_{ij}(z,w)-\Pij_{ij}(-w,w)}{z+w}=\s |
| `eq:set-EdN` | 529 | equation in lem (Edge factors as generating functions) |
| `rmk:set-edgesym` | 539 | rmk |
| `sec:set-graphsum` | 562 | section: Coloured stable graphs and the graph sum |
| `eq:set-Hg` | 566 | equation: Hg_{\bw}=\prod_{a=1}^{3}\frac{\Lamv(\bw_a)}{\bw_a^{h}},\qquad \bw=(E |
| `eq:set-Gsum` | 581 | equation: Gsum_g\coloneqq\sum_{\Gamma}\frac{1}{\|\Aut\Gamma\|}\sum_{\col:\Ver(\G |
| `sec:set-rings` | 592 | section: Rings, general principles, leaf identities |
| `sec:set-hodge` | 636 | section: The Hodge side |
| `sec:set-statements` | 676 | section: The main statements |

### s03_rows.tex

| Label | Line | Belongs to |
|---|---|---|
| `sec:row-euler` | 20 | section: The Euler coordinate |
| `eq:row-phase` | 24 | equation: phi_\tau(s)\coloneqq\wA\log s+\alpha_1\log(1-s)-\alpha_2\log(1-\tau |
| `eq:row-phaseder` | 29 | equation: phi_\tau'(s)=\frac{\wA}{s}-\frac{\alpha_1}{1-s}+\frac{\alpha_2\tau}{ |
| `eq:row-dphi` | 34 | equation: dphi(\tau)=\phi_\tau(\sast(\tau))-\phi_0(\wA/E),\qquad \phi_0(\wA/E) |
| `eq:row-dphider` | 39 | equation: dphi_\tau=\frac{\alpha_2\sast}{1-\tau\sast}=\frac{\wA-\wB\sast}{1-\t |
| `eq:row-euler1` | 52 | equation: cW(1-\tau s)&=\alpha_1\log\big(-\cQ\tau(1-s)\big)-\alpha_2\log(1-\ta |
| `eq:row-euler2` | 55 | equation: re^{\cW/z}\,\frac{\dd U}{U}&=\kph(z)\,\cQ^{\alpha_1/z}\,\tau^{E/z+1} |
| `sec:row-hess` | 69 | section: Hessians and the tensor rule |
| `eq:row-algdef` | 73 | equation: alg(\tau)\coloneqq(1-\tau\sast)\Big(\frac{\phi_\tau''(\sast)}{\phi_0 |
| `eq:row-hess` | 84 | equation: phi_\tau''(\sast)=-\frac{\delta}{\sast(1-\sast)U_1},\qquad \phi_\tau |
| `eq:row-alg` | 91 | equation: alg^2=\frac{U_1^2\,\phi_\tau''(\sast)}{\phi_0''(\wA/E)}=\frac{\alpha |
| `eq:row-algp` | 96 | equation: alg^2=\frac{\wA\,\delta\,(E-\wB\tau\sast)}{E^3\,\sast},\qquad E-\wB\ |
| `eq:row-tensor` | 107 | equation: cW_2(U_i)=\frac{U_i^2\,\phi_\tau''(s_i)}{\tau^2}\quad((s_1,s_2)\colo |
| `rmk:row-hessfib` | 126 | rmk |
| `sec:row-wr` | 132 | section: The rescaled derivatives $\tau^k\cW_k(U_1)$ |
| `eq:row-phik` | 136 | equation: begin{gathered} \phi_\tau^{(k)}(s)=(k-1)!\Big[(-1)^{k-1}\wA s^{-k}-\ |
| `eq:row-Wk` | 148 | equation: cW_k=-\alpha_1\Pol_{k-1}(p_\cQ)-\wA\Pol_{k-1}(p_1)\qquad(k\geq2). \e |
| `eq:row-wr` | 157 | equation in lem (Rescaled derivatives; formal) |
| `sec:row-laplace` | 183 | section: The formal Laplace calculus |
| `eq:row-gaussfun` | 191 | equation: gauss{K}{a}\coloneqq\sum_{p\ \mathrm{even},\,q}K_{pq}\,(p-1)!!\,\Big |
| `eq:row-P2` | 249 | equation: de_\tau\mathsf{n}+\Big(\frac{f_{0,\tau}}{z}-\frac{f_{2,\tau}}{2f_2}\ |
| `eq:row-P1` | 282 | equation: Lnorm[f\circ h,(\amp\circ h)\,h';v_0]=h_1\,\Lnorm[f,\amp;u_0]. \end{ |
| `eq:row-laplace` | 296 | equation: int_I\re^{f(u)/z}\amp(u)\,\dd u=\re^{f(u_0)/z}\Big(\frac{2\pi z}{-f' |
| `eq:row-homog` | 323 | equation: wick_m(\lambda^2\cW_2,\lambda^3\cW_3,\dots,\lambda^{2m+2}\cW_{2m+2}) |
| `sec:row-wickrow` | 337 | section: The point-frame row is the normalised expansion of the Euler integrand |
| `eq:row-Ndef` | 345 | equation: begin{gathered} \Nlap(z,\tau)\coloneqq\Lhat[\phi_\tau,(1-\tau s)^{-1 |
| `eq:row-wickrow` | 357 | equation in prop (formal) |
| `eq:row-wick13` | 362 | equation in prop (formal) |
| `rmk:row-unif` | 385 | rmk: The degenerate coordinate; a uniform version without continuation |
| `eq:row-unif` | 393 | equation: cW(\re^{-\tau\tilde y})=\phi_\tau(\varsigma_\tau(\tilde y))+E\log\ta |
| `eq:row-unif2` | 398 | equation: (-\tau)^k\cW_k(U_1)=(-1)^k\wrk_k(\tau)=\frac{\dd^k}{\dd\tilde y^k}\bi |
| `sec:row-gauss` | 406 | section: Gauss's equation for the formal expansion |
| `eq:row-gaussdef` | 410 | equation: begin{gathered} \ga\coloneqq1+\frac{\alpha_2}{z},\qquad \gb\coloneqq |
| `eq:row-gaussop` | 419 | equation: begin{gathered} \Gauss[Y]\coloneqq\tau(1-\tau)Y_{\tau\tau}+\big(\gc- |
| `eq:row-amp` | 435 | equation: tau(1-\tau)\Dam^2\amp_0+\big(\gc-(\ga+\gb+1)\tau\big)\Dam\amp_0-\ga\ |
| `sec:row-closed` | 478 | section: Uniqueness and the closed form |
| `eq:row-closed` | 499 | equation in thm (The closed form; formal) |
| `sec:row-beta` | 510 | section: The Beta constant |
| `eq:row-betaphase` | 518 | equation: Lnorm[f,1;y_0]=h_{(-a,-b,a+b)}(z)\qquad\text{in }\bbQ(a,b)\fps{z}. \ |
| `eq:row-betalap` | 532 | equation: Bfn\Big(\frac az,\frac bz+1\Big)=\re^{f^\ast/z}\Big(\frac{2\pi z\,b} |
| `eq:row-stirling` | 538 | equation: log\Gamma(u)=\Big(u-\frac12\Big)\log u-u+\frac12\log(2\pi)+\sum_{k=1 |
| `eq:row-beta` | 571 | equation in prop (The Beta constant; analytic input) |
| `rmk:row-betaroutes` | 600 | rmk |
| `ex:row-beta1` | 607 | example: The coefficient of $z$, by hand |
| `sec:row-main` | 617 | section: Regularity, the frame-B row, and (F2) |
| `eq:row-main` | 634 | equation in thm (The frame-B singular row; formal given \eqref{eq:row-beta}) |
| `sec:row-eik` | 655 | section: A purely formal proof of the regularity: eikonal and transport equations |
| `eq:row-eik` | 663 | equation in lem (Eikonal and transport decomposition; formal) |
| `eq:row-transport` | 722 | equation: Ephi\,\Fhyp=\Phia,\qquad \Phia\in\Oan(\Ddisc)\fps{z},\qquad \coef{z^ |
| `sec:row-first` | 755 | section: First orders |
| `sec:row-canon` | 776 | section: Canonical coordinates |

### s04_leaves.tex

| Label | Line | Belongs to |
|---|---|---|
| `sec:lea-statements` | 25 | section: Statements |
| `sec:lea-qde` | 75 | section: The metric in the Euler coordinate, the key formula, the quantum differential equation and the symplectic constant at the singular idempotent |
| `eq:lea-rho` | 86 | equation: t^2\rho^2\Delta_1 = -\epsilon_1\epsilon_2E^3\quad\text{exactly}; \end |
| `eq:lea-op` | 122 | equation: (\Dd+\omq)^2 - \Bigl(\frac{E}{\tau}-\alpha_0\Bigr)(\Dd+\omq) + \frac{ |
| `eq:lea-qde` | 170 | equation: z\,\leaVn' = (\Aop-\Lv_1)\leaVn . \end{equation} \end{prop}  \begin{p |
| `eq:lea-conj` | 178 | equation: Lv_1 + z\Dop_1 = (t\rho)^{-1}\bigl(z\,\tfrac{d}{dt} + \Lv_1\bigr)(t\ |
| `eq:lea-Hrow` | 205 | equation: rr_{1,1} = \rr_{0,1}\Bigl[\frac{E+z}{\tau} - (\wA+z) + z(1-\tau)\fra |
| `eq:lea-proj` | 227 | equation: langle(\Aop-\Lv_1)a,b\rangle = \frac{\Lv_2-\Lv_1}{\Delta_2}\pim(a)\p |
| `eq:lea-E11sigma` | 263 | equation: rho^2\Edge_{11}(\psi,\psi') = \frac{-\epsilon_1\epsilon_2E^3/t^2 - \ |
| `rmk:lea-flat` | 300 | rmk: flat sections; not used |
| `sec:lea-L1` | 305 | section: The dilaton leaves: (L1) |
| `eq:lea-Rk` | 333 | equation in lem (the $t$-coefficients of $z\Gk$) |
| `eq:lea-Rinf` | 338 | equation: Rk_k(\infty) = E\,\delta_{k0}. \end{equation} \end{lem}  \begin{proo |
| `eq:lea-pf` | 362 | equation in prop ((L1) in the variable $z$: partial fractions) |
| `eq:lea-gk` | 381 | equation: gl_k = \sum_{j=1}^{k+1}g_j^{(k)}\lf_j = -\frac{1}{E}\sum_{j=1}^{k+1} |
| `sec:lea-sigma` | 397 | section: The H-row in leaf form |
| `eq:lea-sigma1` | 401 | equation: begin{gathered} \rho\,\rr_{1,1}(z) = \Bigl(\Lv_1+\frac{z}{t}\Bigr)\G |
| `eq:lea-L1exp` | 410 | equation: Lv_1 = \frac{E}{t} + \lrem(t),\qquad \lrem = \frac{\wA}{\tau\sast} - |
| `eq:lea-sigk` | 421 | equation in prop (the leaf structure of $\brow$) |
| `eq:lea-Sj` | 427 | equation: brow^{(k)} = \sum_{j=1}^{k+1}S_j^{(k)}\lf_j,\qquad S_j^{(k)} := \Big |
| `sec:lea-flow` | 465 | section: The $\Edge_{11}$ flow: (F3) for $\Edge_{11}$, (L3), and (L2) for $j = 1$ |
| `eq:lea-flow` | 491 | equation in thm (the $\Edge_{11}$ flow) |
| `eq:lea-Theta` | 493 | equation in thm (the $\Edge_{11}$ flow) |
| `eq:lea-L3` | 522 | equation in thm ((L3), all $k$) |
| `eq:lea-anti` | 531 | equation: e^{(k)}_{k_1+1,k_2+1} = (-1)^{k+1}\,\epsilon_1\epsilon_2\alpha_1\alph |
| `eq:lea-L2one` | 555 | equation: coef{t^k}(\rho f_{1,b}) = \sum_{i=1}^{k+1}\Bigl(\sum_{m=0}^{k}\coef{ |
| `sec:lea-regular` | 572 | section: The regular colour |
| `eq:lea-f0` | 612 | equation: f_0' = (1-\tau)\frac{\alpha_2\stil}{1-\tau\stil} = \wA - \wB\stil = \ |
| `eq:lea-row2reg` | 634 | equation in thm (Row 2 of the quantum differential equation at the regular idempotent) |
| `eq:lea-stepA` | 643 | equation: rpt_{0,2}(z;t) = U_2\,\leaLt(\tau),\qquad \leaLt(\tau) := \Lnorm[\ph |
| `eq:lea-stepB` | 650 | equation: (\Lv_2 + z\Dop_2)\circ U_2 = U_2\circ(z\leaDt + \omq). \end{equation} |
| `eq:lea-stepC` | 657 | equation: Gauss_{\stil}[\leaLt] := \tau(1-\tau)\leaDtau^2\leaLt + (\gc - (\ga+ |
| `eq:lea-stepD` | 668 | equation: (z\leaDt+\omq)^2 - \Bigl(\frac{E}{\tau}-\alpha_0\Bigr)(z\leaDt+\omq) |
| `sec:lea-legs` | 769 | section: The legs: (L2) for $j = 2$, (F3) for $\Edge_{12}$ and $\Edge_{21}$ |
| `eq:lea-legext` | 775 | equation in lem (extraction of a leg) |
| `eq:lea-N12` | 801 | equation: (\psi+w)\,\rho\Edge_{12}(\psi,w) = \leaNtw(\psi,w) - \leaNtw(-w,w). \ |
| `eq:lea-N12coef` | 806 | equation: coef{t^k}\leaNtw(\psi,w) = -\sum_{k_1=-1}^{k}\langle\brow^{(k_1)}(\p |
| `eq:lea-legrat` | 813 | equation: begin{gathered} \coef{t^k}(\rho f_{2,b})(\psi) = \sum_{n=0}^{b}(-1)^ |
| `eq:lea-sigleaf` | 834 | equation in prop (leaf form of the $\Edge_{12}$ edge) |
| `eq:lea-E12leaf` | 836 | equation: rho\Edge_{12}(\psi,w) &= -\sum_{j\geq 1}j\,\langle S_j(t),\vv_2(-w;t |
| `eq:lea-E12coef` | 841 | equation: coef{t^k}(\rho f_{2,b})(\psi) = -\sum_{j=1}^{k+1}j\,\lf_j(\psi)\sum_ |
| `sec:lea-assembly` | 899 | section: Proof of Theorems A and \ref{thm:lea-qde}; (F1)--(F3) in one table |

### s04b_closedforms.tex

| Label | Line | Belongs to |
|---|---|---|
| `sec:lea-closed` | 5 | section: Closed forms |
| `eq:lea-Aj` | 16 | equation in thm (the dilaton leaf coefficients) |
| `eq:lea-loops` | 97 | equation in thm (the loops in closed form) |
| `eq:lea-antidiag` | 115 | equation in cor (anti-diagonal and corner) |
| `sec:lea-a0` | 165 | section: The case $\alpha_0 = 0$: rescaling to the fibrewise rows |

### s05_identification.tex

| Label | Line | Belongs to |
|---|---|---|
| `eq:ide-Delta` | 104 | equation: Delta_i=\epsilon_1\epsilon_2\big(\alpha_0\alpha_1\alpha_2+2\alpha_1\ |
| `eq:ide-mumford` | 126 | equation in lem (Hodge classes as an R-action) |
| `eq:ide-flat` | 198 | equation in lem (Ambiguity of flat symplectic solutions) |
| `eq:ide-roweq` | 217 | equation: z\Dop_i\vv_i=(\Aop-\Lv_i)\vv_i,\qquad \Dop_i=\frac{\dd}{\dd t}-\frac{ |
| `eq:ide-Fl1` | 232 | equation: begin{gathered} \Omega_{g,n}(x_1,\dots,x_n)=\sum_{i=\pm}\ideP_i\,\id |
| `eq:ide-Fl3` | 253 | equation: p_*\Omega_{g,n+1}(x,H)=\sum_i\ideP_i\ideN_i\sum_{l=1}^n\fid_i^*\Big(\ |
| `eq:ide-Fl5` | 270 | equation: sum_i\ideP_i\ideN_i\Big[(\Mdef(\psi_1)\mathrm{a}_1)^i\,\mathrm{a}_2^ |
| `eq:ide-euler` | 363 | equation: idePF\,\re^{\cW/z}=\thU\big(\idekm\,\re^{\cW/z}\big),\qquad\idekm\co |

### s06_regularity.tex

| Label | Line | Belongs to |
|---|---|---|
| `sec:reg-hodge` | 15 | section: The Hodge side: domain, dilaton homogeneity, range of genera |
| `sec:reg-string` | 79 | section: String and dilaton identities |
| `eq:reg-S` | 87 | equation: frac{\lf_j-\lf_j(0)}{\psi}=\frac{1}{\psi}\Big(\frac{1}{E-j\psi}-\fra |
| `sec:reg-cited` | 137 | section: The two cited theorems |
| `lem:reg-scale` | 152 | lem: Scale law |
| `eq:reg-scale` | 155 | equation: Xiz^{\mathsf{v}}_h(\msU)=E^{-(3h-3+2r)}\,\Xit_h(\msU)\big\|_{\tau=\ma |
| `eq:reg-laurent` | 173 | equation: Xit_h(\msU;\tau)=\sum_{k=0}^h\Big(-\frac{1}{\st}\Big)^k\,\Xi^{(k)}_h |
| `eq:reg-VI` | 198 | equation: (-\epsilon_1\epsilon_2E^3)^{g-1}\,\ZZ[g,0]=\beta_{2g-2}(\epsilon_1,\e |
| `rmk:reg-betavalues` | 239 | rmk: The first values |
| `sec:reg-thmR` | 251 | section: Theorem R: regularity and the polar part of the decoration sums |
| `eq:reg-F1` | 255 | equation: Delta_1\,\rho^2=-\frac{\epsilon_1\epsilon_2E^3}{t^2}. \end{equation} |
| `eq:reg-pole` | 266 | equation: (-\epsilon_1\epsilon_2E^3)^{g-1}\,\Fpot_g(\qo)\,t^{2-2g}=(-\epsilon_1 |
| `eq:reg-vertex` | 276 | equation: Delta_i^{\,h-1}\,\dd^n\Fpot_h\big(\qt^{(i)}\big)[f_1,\dots,f_n], \en |
| `eq:reg-resc` | 285 | equation: Delta_1^{\,h-1}\,\dd^n\Fpot_h(\qt)[f_1,\dots,f_n] =\Big(\frac{-\epsi |
| `sec:reg-thmC` | 335 | section: Theorem C |
| `eq:reg-gap` | 344 | equation: gapg\,t^{2-2g},\qquad \gapg\coloneqq(-\epsilon_1\epsilon_2E^3)^{g-1} |

### s07_rigidity.tex

| Label | Line | Belongs to |
|---|---|---|
| `rmk:rig-rows` | 44 | rmk: what is proved about the rows |
| `eq:rig-NG` | 107 | equation: sG_g=\frac{N}{(1-\cQ)^n\,\Disc^m},\qquad N\in\KCY[\cQ],\quad n,m\geq |
| `eq:rig-Dinf` | 145 | equation: end{align} both non-zero in $\KCY$. The same values follow from $\De |
| `rmk:rig-g1input` | 262 | rmk |
| `rmk:rig-3Trig` | 415 | rmk |

### Changes to the registered labels

- `thm:set-A`, `thm:set-B`, `thm:set-C`, `thm:set-D` (the duplicated statements at the end of s02) were removed; `sec:set-statements` is now a pointer paragraph. The statements are `thm:int-leaf`, `thm:int-ident`, `thm:int-reg`, `thm:int-rig`, `cor:int-conj` in s01, as prescribed above.
- s01 has the additional subsection label `sec:int-BS` (relation to the results of [BS]); the three-torus corollary is not a separate numbered statement in s01 (it is stated in words with a pointer to `cor:rig-3T`).
- `lem:set-N` has outer parts (a), (b) and inner identities (i)-(iv) under (b), so that `\cref{lem:set-N}(b)(iii)` and "identity (ii) of `lem:set-N`" read as the writers intended.
- `def:set-inv` now also defines $P_g$, $\beta_n$ and $\mathrm{Li}_{-n}$ (previously defined only inside the duplicated Theorem C of s02; s06 and s07 recall them from `def:set-inv`).

## Macros defined in the section files

All section-local macros are `\providecommand`s (s05: `\newcommand`s) at the top of their file; none clashes with main.tex. Three macros that two sections defined independently were moved to main.tex: `\Dam` (s03 and s04, same meaning), `\amp` (s03, now also used in s04 for the amplitude of the regular-row Laplace expansion), and the local coordinate $1/\mathcal{Q}$ at infinity, which s02 called `\sinv` and s07 `\sQ`: it is now `\sQ` everywhere (s02, s07, appendix).

| Macro | File | Definition | Meaning |
|---|---|---|---|
| `\amp` | main.tex | `\mathsf{m}` | amplitude of a Laplace expansion (s03, s04) |
| `\Dam` | main.tex | `\mathsf{D}_{\mathrm{am}}` | d_tau q + q (d_tau f)/z (s03, s04) |
| `\sQ` | main.tex | `\mathsf{s}` | 1/Q, the local coordinate at Q = infinity (s02, s07, appendix); not the Euler variable |
| `\pc` | s02_setup.tex | `\mathsf{p}` | pc_{k,i}: coefficients of the polynomials Pol_k |
| `\kE` | s03_rows.tex | `k_E` | s(1-s)(1-tau s)^{-2} |
| `\kph` | s03_rows.tex | `\kappa_{\mathrm{ph}}` | constant of the phase, used nowhere |
| `\Lt` | s03_rows.tex | `\Lnorm_{\tau}` | Lnorm[phi_tau, (1 - tau s)^{-1}; s^*] |
| `\Ngauss` | s03_rows.tex | `\Nlap_{\mathrm{G}}` | the conjugated expansion that solves Gauss's equation |
| `\Spht` | s03_rows.tex | `\mathsf{S}_{\phi,\tau}` | d/dtau of the eikonal phase (avoids the double subscript \Sph_\tau) |
| `\Sphtt` | s03_rows.tex | `\mathsf{S}_{\phi,\tau\tau}` | second tau-derivative of the eikonal phase |
| `\Abetam{m}` | s03_rows.tex | `\mathsf{A}_{\beta,m}` | z^m coefficient of \Abeta (avoids the double subscript \Abeta_m) |
| `\Bfn` | s03_rows.tex | `\mathrm{B}` | Euler's Beta function |
| `\leaVn` | s04_leaves.tex | `\mathbf{n}` | normalised vector t rho v_1 |
| `\leaPr` | s04_leaves.tex | `\mathrm{Pr}` | projectors Pr_1, Pr_2 onto the eigenlines |
| `\leaPhiF` | s04_leaves.tex | `\boldsymbol{\Phi}_{\Fhyp}` | (F, (D_z + omega) F) in the Wronskian remark |
| `\leaLt` | s04_leaves.tex | `\Lnorm_{2}` | Lnorm[phi_tau, (1 - tau s)^{-1}; s~] |
| `\leaDt` | s04_leaves.tex | `\mathsf{D}_{\stil}` | d/dt + f_0'/z - f_2'/(2 f_2) |
| `\leaDtau` | s04_leaves.tex | `\mathsf{D}_{\stil,\tau}` | (1 - tau)^{-1} D_{s~} |
| `\leaNtw` | s04_leaves.tex | `\mathsf{N}_{12}` | rho P_{12} |
| `\leaDisc` | s04_leaves.tex | `\Ddisc_{1}` | a disc around tau_0 with closure in D'' |
| `\TFsup` | s04_leaves.tex | `\mathrm{F}` | superscript for the fibrewise point (0, x, 1, x, 1) |
| `\ideKt` | s05_identification.tex | `\widetilde{\bbK}_{\mathrm{CY}}` | K_CY with e(wt_+)^{1/2}, e(wt_-)^{1/2} adjoined |
| `\ideFz` | s05_identification.tex | `\bbF_{\circ}` | ideKt((Q)), the field of the flatness proof |
| `\ideone` | s05_identification.tex | `\mathbf{1}` | the unit of the state space; also the trivial rank-one CohFT |
| `\ideV` | s05_identification.tex | `\mathcal{V}` | state space K_CY 1 + K_CY H |
| `\idePF` | s05_identification.tex | `\mathsf{PF}` | the second-order operator of Euler's identity |
| `\idekm` | s05_identification.tex | `\mathsf{k}_{\mathrm{m}}` | its amplitude -z alpha_1 Q(1-U)/(1-QU) |
| `\ideCh` | s05_identification.tex | `\mathcal{X}` | ideCh_h = sum_k b_k (2k)!/B_{2k} ch_{2k-1}(E) |
| `\ideCurve` | s05_identification.tex | `\mathfrak{C}` | universal curve over the open moduli space |
| `\idesec` | s05_identification.tex | `\varsigma` | its sections |
| `\ideLam` | s05_identification.tex | `\boldsymbol{\Lambda}` | diag(L_+, L_-) |
| `\idecan` | s05_identification.tex | `\boldsymbol{u}` | diag(u_+, u_-), canonical coordinates |
| `\ideS` | s05_identification.tex | `\mathscr{S}` | fundamental solution R Psi e^{u/z} |
| `\idelg{i}` | s05_identification.tex | `\mathrm{a}^{i}` | leg functions a^i(psi; x) = f_i^*(B(psi) x) |
| `\ideP` | s05_identification.tex | `\mathfrak{P}` | P_i = Delta_i^{(2g-2+n)/2} |
| `\ideN` | s05_identification.tex | `\mathfrak{N}` | kappa series of the translation legs |
| `\ideG` | s05_identification.tex | `\mathfrak{G}` | kappa series of the pushforward |
| `\ideW` | s05_identification.tex | `\mathsf{W}_{\Psi}` | Psi^{-1} theta_Q Psi |
| `\ideRf` | s05_identification.tex | `\mathsf{R}^{\mathrm{f}}` | Psi^{-1} R Psi |
| `\ideY` | s05_identification.tex | `\mathsf{Y}` | the pairing <v_+(z), v_-(-z)> in the symplecticity step |
| `\ideZ` | s05_identification.tex | `\mathsf{Z}` | Psi^{-1} R_1^{-1} R_2 Psi in the ambiguity lemma |
| `\idesp` | s05_identification.tex | `\mathfrak{s}` | the scalar series s_p = exp sum_k a_k^{(p)} z^{2k-1} |
| `\ideQU` | s05_identification.tex | `\mathrm{V}` | V = Q U_-, the regular coordinate at Q = 0 |
| `\ideT` | s05_identification.tex | `\mathsf{T}_{\!\Rm}` | translation z(1 - R^{-1}(z) 1) |
| `\regbra{f}{U}` | s06_regularity.tex | `\langle f\mid U\rangle` | one extra insertion f(psi_0) in Xi_h(U) |
| `\regD` | s06_regularity.tex | `\dg` | total leaf degree at a vertex, regD_v |
