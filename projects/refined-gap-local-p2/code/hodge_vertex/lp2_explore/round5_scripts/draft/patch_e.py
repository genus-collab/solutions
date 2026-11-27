# Round-5 patch, part E: s.8, s.9, Appendix.
PE = []
PE.append((737, 769, '    Statement ', r'''    Statement                                        Status      Rests on / source
    ------------------------------------------------ ----------- ------------------------------------------------
    Lemma 1.2(a),(b) (t_c, ρ, ρ^2 Delta_0)           P           algebra
    Lemma 1.2(c) (t_c = t_cf)                        P mod C     BS 5145, eq:pfinf (read by round-4 critics);
                                                                 used only for Corollary G
    Lemma 1.4 (code bookkeeping)                     P           reading hv_lib3.m, hv_lib4.m, lp2_lib.m, lp2_percolour.wls
    Theorem A, Cor 2.2 (PF, uniqueness)              P           algebra; R3 Part I; Theorem A re-derived by round-5 referee
    SP1-SP5 (formal stationary phase)                S (written) R3 Part C; interface only
    Lemma 2.1 (critical points)                      P           R2:Lemma 1; interface only
    Lemma 2.3 (analyticity)                          P           R2:Lemma 2; interface only
    Lemma 2.4 (PF under the integral)                P           SP2, SP3; interface only
    Lemma 2.5 (WKB uniqueness)                       P           R2:Lemma 4; interface only
    Prop 2.6 (normalisation; eq:normNLR derived)     P mod C     (C1), (C2); R2:Prop 5; interface only
    Thm 2.7 (continuity, closed form, Rec = Wick)    P mod C     (C1)-(C3); R2:Thm 6, Prop 8; R3:Cor B; interface only;
                                                                 reading of BS 5257-5267 re-checked by referee
    Cor 2.8 = Lemma CL (corrected)                   P           algebra; rationality at colour 0, frame B only
    Cor 2.9: F2, L1, Form, R, Z-infinity             P           closed form + Theorem A
    Prop 2.10 (symplecticity, K = 1, orthogonality)  P           Theorem A only; NEW in round 5, one referee pass
    Prop 2.11 (Givental-Teleman on T_F)              P mod C     (C1)-(C3), (C5); not used for Theorems LC_P2, R_P2
    Prop 2.12 (lem:modrel, both relations)           P           no citation; NEW in round 5, one referee pass;
                                                                 interface only
    Lemma D0 (L1, F3 for Dil_0)                      P           Cor 2.9; written here
    Prop 3.1 (closure, general c1)                   P           Theorem A; not used (its c1 = 1 form is S1 of Prop 2.10)
    Lemmas 4.1-4.3, Prop H2, Lemma G                 P           s.3.1, Prop 2.10
    Lemma R, Theorem E00                             P           s.3.1, H2 Steps 1-2, Prop 2.10
    Theorem M                                        P           Cor 2.9(R), Prop 2.10, H2
    Theorem LC_P2                                    P           for the decoration sum as defined; no BS input
    Fpot_h on 𝔇, (H0)                                P           dilaton equation; NEW in this edit, unrefereed
    (S_P2), (C_P2), (L0_P2)                          P           string/dilaton equations
    Lemma 7.1 (range)                                P           stability
    Lemma A_P2(a),(c)                                P           degree count
    Lemma A_P2(b) (vanishing)                        C           lemmaA/LEMMA_A_PROOF.md Theorem L, P there modulo (MV) and
                                                                 textbook facts; becomes P mod (MV) after its open point O1;
                                                                 see s.7.4 and 8.2
    Theorem R_P2                                     P mod C     Theorem LC_P2, Lemma A_P2(b) for h <= g-1
    Theorem VI (vertex identity (VI))                P mod C     LLZ eqs. (1), (5), boson-fermion, string/dilaton, (7.2);
                                                                 NEW in round 5, one referee pass
    Corollary SV_P2                                  P mod C     above + Lemma A_P2(b) at h = g, Theorem VI
    Proposition 13.2_P2 (G_g = GW^cf_g)              O           s.8.3
    Corollary G (BS rmk:cfgap)                       Conditional on Proposition 13.2_P2. NOT CLAIMED.

"Interface only" means: used for Proposition 2.11, Proposition 13.2_P2 or Corollary G, and not in the proofs of Theorems LC_P2, R_P2 or Corollary SV_P2. Round-4 critics asked for H2 Step 3, E00 and (M-d) to be relabelled P mod (C1)-(C3); since Proposition 2.10 is now P, that relabelling is unnecessary. All of these theorems are about the decoration sum (0.1) as defined.'''))

PE.append((771, 779, '### 8.2 Remaining inputs', r'''### 8.2 Remaining inputs for Theorem R_P2 (and Corollary SV_P2)

The following is the outcome of a dependency reading, not of a full referee pass over the document.
1. **Lemma A: LEMMA_A_PROOF.md (s.7.4).** Literature: (MV), the Marino-Vafa formula, LLZ eq. (1) (arXiv v3 read, JDG numbering not confirmed), with OP Theorem 1 as a second published proof; and the textbook facts listed in s.7.4. The 2D Toda/mKP equations are not a citation in the proof of record, since the P_k are derived in s.4 of that file. The matching of the MV series with Xi_h is no longer an open item (that file's Proposition 2.1, s.3 and Proposition 8.1). Open there: O1, the assembled text has not had its own referee; O2, citation details, none mathematical; O3, the proof of (MV) in LLZ was not re-checked.
2. **Theorem VI, for Corollary SV_P2 only.** LLZ eq. (1) and Theorem 1 eq. (5); the boson-fermion correspondence (vertex-operator form of ψ, ψ^*); the principal specialisation. The string and dilaton equations and (7.2) are standard or proved here. The ABDKS/Eynard route is no longer used.
3. **BS.** No BS input enters Theorems LC_P2 and R_P2, or Corollary SV_P2 on the recorded route. This rests on the dependency reading of s.3-s.7 by the round-5 modrel author and referee; the η^{ab} and the operator eq:pf enter as definitions. Identification with GW still needs (C1)-(C5) and Proposition 13.2_P2, s.8.3.

Unconfirmed citations are collected in s.7.4 and s.7.8: JDG numbering of LLZ; theorem numbers for the mKP bilinear identity in Date-Jimbo-Kashiwara-Miwa and Miwa-Jimbo-Date; Macdonald I.3 Ex. 2. F2 has **not** been taken from thm:cfasym: since round 2 it follows from the closed form, and eq:Rcf is then a consistency check.'''))

PE.append((794, 794, '2. *N4.*', r'''2. *N4: closed.* lem:modrel, including the second relation E(y) = (D - X)/2, is proved in Proposition 2.12, and it supplies the dictionary Δ_c A^2 = 1/D, d2 = e_1(α*) E(y) identifying Definition 1.3 at c1 = A with eq:RHb.'''))

PE.append((799, 799, '**Remark on N4', r'''**Remark on N4 and Proposition 2.10.** In round 4, Proposition 2.10 step 3 compared with BS's rows at c1 = A, which made lem:modrel a load-bearing input of Theorem LC_P2, and the draft wrongly attributed c2 = D to eq:ADEFXdef (it is the first relation of lem:modrel). Both issues are gone: Proposition 2.10 is now proved from Theorem A, and lem:modrel is Proposition 2.12. The audit of BS's own use of lem:modrel is Remark 2.13.

Round-1 items N1-N8 are recorded in journal wf_e928f771-9e3 (only N4 appears in ROUND1_STATUS_2026-11-02.md).'''))

PE.append((824, 824, '11. **lem:modrel', r'''11. **lem:modrel (2695)** E = (D - X)/2: BS's proof is a sketch ("We omit the details"), and the first relation is cited. Fix: Proposition 2.12, whose proof can replace 2697-2749 verbatim. Also: eq:SHb (about 3299) cites eq:allc, whose printed row 2 already has modrel applied; it should cite eq:ssH1 (Remark 2.13(a)). BS's derivation of Delta_i^{-1} passes through modrel-simplified formulas; a modrel-free derivation is in Remark 2.13(c).'''))

PE.append((827, 827, '14. **Line references.**', r'''14. **Line references.** eq:pf label at 2628 (displayed about 2627-2632); eq:Wadm at 3225; eq:Rcf at 5334; thm:cfasym at 5301; lem:cfpole at 5407; rmk:cfgap at 5422; eq:Ifun at 2621; eq:sympl at 3084; eq:SHb at 3305; eq:RHb at 3573; prop:ResConifold at 1708 ("Prop. 3.1" is the numbering of arXiv:2410.00118, not a v2 label); prop:cfgraphsum at 6366. Earlier round notes quote slightly different numbers.'''))

PE.append((830, 830, '## Appendix.', r'''## Appendix. Computer checks (evidence only; none is a proof step)

The round-4 critic's objection, that the items S1, S3 and proof (ii) of S2 made finite-order evidence load-bearing, no longer applies: those steps now have written all-order proofs (s.7.4, s.7.8), and gap_identity_check.pl, dexrel_check.pl and the data comparisons are demoted to evidence.'''))

PE.append((846, 846, '- dexrel_check.pl:', r'''- dexrel_check.pl: E(y) = (D - X)/2 to y^40. Evidence only since Proposition 2.12.'''))

PE.append((863, 863, '**Lemma A data', r'''**Lemma A data (hodge_vertex_conifold/).** hvc_multiprobe_out.txt, hvc_multiprobe_g4_out.txt, hvc_multiprobe_g4b.wls (M_irr(3), M_irr(4)); lemmaA/E_formulas/mirr_g7.txt (all 194 genus-7 points zero, 77 controls non-zero); lemmaA/round2/critic/gap_identity_check.pl (44 comparisons through h = 12).

**Round 5 (SP/lp2_round5/). Exact rational arithmetic, or F_p where stated. Checks by the MVmatch and ThmL owners and referees are recorded in LEMMA_A_PROOF.md and journal wf_9d7fd94e-8d7.**
- modrel/modrel_series.py: L2b A = 0, (K), L2b P_I = 0, A^2 D = X, 2E(y) = D - X and eq:limitconv through y^60; eq:Ian coefficients reproduced from eq:Ifun for d <= 6 at a generic rational α.
- modrel/pairing_check.py: Π_ij = δ_ij on the closed form alone, all i <= j, through z^6 at (1,1), (5,2), and through z^5 at (7,-3), (1/2,3/5). Controls (r_2 perturbed, η^00 doubled) fail.
- modrel/closure_frob_check.py: d2 = e_1 E(y) and Δ_c c1^2 = 1/D at c1 = A through t^24 (Laurent truncation); Definition 1.3 = eq:RHb through z^4 up to t^4 at colour 0 and (5,2), (1,1); floating-point test of the modrel-free Phi_2i.
- modrel_ref1/thm3_fp.py (F_p, p = 2^61 - 1, random section c1, random t0): Π_ij = δ_ij for all colour pairs through z^14 at 9 weight pairs; Proposition 3.1 for general c1 through z^13. Controls (d2 doubled, r_3 perturbed, η^00 doubled) fail.
- modrel_ref1/thm1_direct.py: I-function expanded directly from eq:Ifun at 4 α, including T_F; A^2 D = X, 2E(y) = D - X, eq:DEXrel through y^22; (K) through y^200. Control fails.
- modrel_ref1/generic_alpha.py (F_p): at generic α, η_ab from Remark 2.13(c), and an independent WKB solution of eq:pf with rows from eq:RHb: B_ij = 0 for i != j, θ B_ii = 0 through z^5. Controls fail.
- VI/vi_check.py: (*) on the Schur side, p_1^n, n <= 7 at τ = 2/7, -5/3; Z_h extracted from MV data, h = 2..6, equal to β_{2h-2}; (E1) with dual numbers in τ.
- VI/vi_grade.py: (*) in monomials x^a w^b, x-order <= 12: grade <= 1 vanishes with Z_h = β_{2h-2}/σ^{h-1}. Perturbing Z_3 gives non-zero grade 1.
- VI_ref1/faber.py, h5.py: MV-free values (Mumford/Faber recursion and DVV) of σ^{h-1} Z_h = β_{2h-2}(-τ, 1+τ) for h = 2, 3, 4 at τ = 3/11, -7/4, 5/2, and for h = 5 at τ = 3/11. Control β·σ fails.
- VI_ref1/lastchecks.py: slice normalisation (B1) against Faber integrals for (h, a) up to (3, 2); (E1) with no mismatch for N^1 down to N^{-33} at τ = 3/11, -7/4, 13/5. This closes the top-order flag of VI/vi_check.py as a truncation artefact. Kernel control: adding (5/3) f_* ((1+τ)t)^4, which breaks τ-symmetry, still solves (E1), so the symmetry step is genuinely needed.
- VI_ref1/schur_modp.py (mod 2^61 - 1): (*) on p_1^n, n <= 12; Z_h for h = 2..8 at τ = 3/11, -7/4, 13/5 by the w-polynomial fit, equal to β.
- VI_ref1/grade.py: with B_s and higher grades as arbitrary random data, grade 0 of (*) vanishes and grade 1 equals 𝓛[B_s] - RHS(E1) exactly. Controls: perturbing the log coefficient or the genus-1 constant fails, and perturbing the -1/4 of G_0 changes nothing, as 𝓛[N^2/σ] = 0 predicts.
- draft/check_transcription.py: this text's formulas for (E1), 𝔏 and the (H0) base-point identity, re-evaluated exactly with the notation of s.7.8.'''))
