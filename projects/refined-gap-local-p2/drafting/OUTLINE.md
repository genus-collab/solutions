# Outline of the paper

Working file for the section writers. FB = FRAMEB_PROOF.md, DR = LP2_REGULARITY_PROOF_draft.md, LA = LEMMA_A_PROOF.md, ST = ROUND6_STATUS_2026-11-02.md, BS = BS_v2_2026-11-02/main.tex, CP = leaf_structure_proof.md. Notation: NOTATION.md. Labels: LABELS.md.

The order of sections is fixed. Results are used only after they are stated and proved, and to achieve this a few source items change place (marked **moved**). Finite-order computation appears only in the appendix. Every all-order proof step of a source is written out.

## Global dependency graph

    s02 definitions
    s03 closed form: Thm PF → Cor unique; SP calculus → Lemmas crit, analytic, pfint, wkb; Lemma closed12 → Cor carry; Cor closed; Prop sympl; Prop modrel
    s04 uses s03 only in prop:gen-osc
    s05 uses s03, s04
    s06 uses s03 (Prop sympl, Prop modrel), s04 (Lemma mumford), s05
    s07 uses s02, s03 only (no s04-s06)
    s08 uses nothing from s02-s07 beyond the definition of Ξ
    s09 uses s08: prop:van-real, lem:van-dil, eq:van-bil at k = 2, cor:van-scale; nothing else
    s10 uses s07, s08, s09 (statements about the decoration sum), then s06 for the corollaries
    s11 points to the above

The regularity statements (s07, s08, s09, first half of s10) do not depend on the mirror theorem or on Teleman. Keep that visible: state it in s01 and at the start of s07.

## s01_intro

Covers: FB s.0.2 (statements only), s.0.3 (scope, not proved, dependency, end-to-end inputs), DR s.0.2, s.0.4 (not claimed), ST s.1, s.3.2. BS's def:dirint and rmk:cfgap as background.

Contents, in order.
1. Set-up in words: Z = K_P2 x C², fibrewise torus TF, generating series ℱ_g^{TF}, conifold substitution Cf, t_cf.
2. Statements: `thm:int-frame`, `thm:int-regularity` with the polar part, `cor:int-gap`, `cor:int-gapkh`. Each with exactly the hypotheses and scope of its body version: fibrewise torus, per colouring where the statement is per colouring. The meaning of O(1) in full: the difference is holomorphic in a neighbourhood of t_cf = 0, Taylor coefficients in 𝕂, value not determined; either convention for d = 0. Genus 2 display from BS's remark.
3. Inputs, listed once: FB s.0.3 end-to-end list, as "the proofs use the following published results". [CCIT15, Theorem 31]; the genus-0 differential equation is derived from [Pan98, Lemma 1] and the boundary structure (not an input); [Tel12, Theorem (2), Prop 8.7(ii), 2.18, 2.19] with [Loo96], [MW07], [Har85], [Iva90]; [PPZ15, Defs 2.2, 2.5, 2.13, Prop 2.14]; [GP99]; virtual-class axioms with [CG07, App. 1 (iv)-(vi)]; [FP00, s.1.1, Theorem 4] for Mumford; [DLMF, 15.6.7], Barnes's second lemma, Gauss's values; [LLZ03, eq. (1)] with textbook facts (string, dilaton, genus-0 multinomial formula, Frobenius, Murnaghan-Nakayama, principal specialisation, Lagrange inversion). Say which results use which: regularity and the polar value use only [LLZ03] and textbook facts; the mirror theorem enters through the frame formula only (FB head, dependency partition).
4. What is not claimed: FB s.0.3 list, DR s.0.4: nothing off TF; the value of the regular part; ℱ_{k,h} in Q[X^{±1},F] (only ℱ_g^{TF} in 𝕂[y^{±1},ξ^{-1}][F]); Orbifold regularity; holomorphic anomaly equations; uniqueness in direct integration. The proofs of the cited theorems are not reproduced. Not used: topological recursion on the mirror curve of local P2, remodelling, BS's thm:cfasym, lem:cfpole, prop:cfgraphsum, and any finite-order computation.
5. Strategy: three steps, generic weights → ring map → frame change; closed form → leaf structure; Mariño-Vafa → vanishing; counting. The whole polar part sits in the single-vertex graph.
6. Plan of the paper.

\draftnote at the Barnes citation: the form used should be compared with a printed source before submission (FB s.1.2). \draftnote at Teleman: numbers from arXiv v3; journal numbering not compared. Same for CCIT15 print.

## s02_setup

Covers: DR s.0.1, s.1.1-s.1.6, excluding Lemma 1.2(c) and Lemma 1.4; FB s.0.1; BS definitions (I-function, A, D, E(y), F, X, I_n^{[a]}, mirror map, η, J, S, 𝒲, dΩ, Asym, Θ, Υ) as cited definitions.

Order.
1. Target, tori, tangent weights, pairing (BS); invariants by localisation; ℱ_g with d ≥ 0: `rmk:set-d0`.
2. Mirror data: `eq:set-I`, `eq:set-pf`, mirror map (with g(y) independent of α), `eq:set-gens`, L̂, W, Y_i. Write-out of FB s.0.1 Base.
3. Colours on TF: t_i, Y_i, q_i, v_i, V_i, W(Y_i) = E q_i², Δ_i on TF and fixed Δ_i^{1/2} := sqrt(3ε1ε2E)/v_i. Remove the sentence about edges_sing.pl and code names.
4. Leaves and (LAM1)-(LAM4). **λ** normalisation warning dropped (historical).
5. Conifold coordinate: 𝓛_2, w_cf, t_c, c_cf, c_1, F_c, d_2. `lem:set-coord` with proof, (a), (b) only. Say that (c) is `lem:fra-tcf`. DR's remark after Lemma 1.2 on rescaling t_c: keep its mathematical content (normalisation enters through ρ(0) = 1).
6. Special functions: Θ, two-variable Θ, h_w, N(z), (Rec), S, S^{[m]}. Polynomiality of [z^k]S in v.
7. Rows: `def:set-rows`, normalised rows, the two remarks. Cf is defined in s06: only say "at c_1 = c_cf". Pairing η^{ab} on TF.
8. Factors: Δ_i, Dil̂_i, q_t, edge `eq:set-edge`, edge symmetry (forward reference to `prop:row-sympl` for the δ form; write the definition without it, and state the δ form after that proposition, in s03). `rmk:set-signs`.
9. Stable graphs and decorations: `def:set-graphsum`, `eq:set-graphsum`. Definition of Cont(Γ,𝔠), 𝒢_g. Finite decoration sum. Colour 0 is singular.
10. Hodge side: Hg, 𝐑, 𝔡^s, 𝔇, F^pot, overlap of the two definitions, Ξ_h, dg, Z[h,n], `def:set-mirr`, `eq:set-vertex`. `rmk:set-frame` (DR s.1.6 without pointers to other notes).
11. Givental action recalled for the reader: legs R^{-1}(ψ), T_R(z) = z(1 - R^{-1}(z)1), edge bivector, as used in FB A2.7, A2.8, C Step 1, citing [PPZ15, Defs 2.2, 2.5, 2.13, Prop 2.14]. Only what the sources state; \draftnote if a formula of PPZ is needed beyond them.

Dropped: DR Lemma 1.4 and all code references; the definition (0.1) stands on its own. DR s.2.5 (status of code) is dropped entirely.

## s03_rows

Covers: DR s.2.1; s.2.2 (calculus); s.2.3 (Lemmas 2.1, 2.3, 2.4); s.2.4 (Lemma 2.5 with its algebraic form); Lemma 4.1 **moved** from DR s.4; s.2.6 (Cor 2.8, Cor 2.9, Prop 2.10); s.2.7 (Prop 2.12, Dictionary).

Not included: Prop 2.6, Thm 2.7 (their content on the route is FB A3.5, B4), Prop 2.11 (now `thm:spe-main`(d)), Remark 2.13 (audit of BS; its (a) is in FB A2.1), DR "Scope" paragraph beyond one neutral sentence. Lemma 2.3 defines Φ and Jet_m here; FB Theorem FB_B(a) then only refers to it.

Order: Thm PF, Cor unique → calculus (G1, SP2, FP, SP3, SP1, SP4, SP5 in the source's logical order; FP uses G1, SP3 uses FP and SP2, SP1 uses FP and SP2) → crit → analytic → pfint → wkb → closed12 → carry → closed → sympl (statement, S0-S6, then the δ form of the edge and edge symmetry) → modrel → dictionary.

Depends on: s02. Used by: s04 (osc), s05 (B2, B3, B4), s06 (sympl, modrel), s07 (all of closed, carry, sympl, closed12, thm PF). Note in the section opening that thm PF, closed12, carry, closed, sympl, modrel are statements about the closed form only, and that the calculus and lemmas crit, analytic, pfint, wkb are used only for the frame formula.

## s04_generic

Covers: FB s.2 entire; FB s.1.2 bullets on eq:J=I and on the QDE; B5 part. Statement `thm:gen-main` first. FB A2.2 as a convention sentence.

Order (**moved** items make this differ from FB's numbering).
1. Statement of the theorem.
2. A1.1, A1.2.
3. **moved** from B5: `prop:gen-poly`: GP for T; d ≥ 1; d = 0, g ≥ 2; passage to TF, including e_3 not in the prime. A1.3 is this proposition. The rings 𝕆_α, 𝔮 are defined here. The constant term stays in s05.
4. A1.4, A1.5, A1.6.
5. **moved** from FB s.1.2: `prop:gen-qde`: the psi relation with its T-equivariant form; boundary restriction; TRR (called the genus-0 topological recursion relation, "the psi-class relation on Mbar_{0,n}", not Eynard-Orantin); V(a,b); string; divisor; restriction to the small locus; commutativity. The double bracket convention.
6. **moved** from FB s.1.2: `thm:gen-mirror`, with proof: scope (fan, semi-projective), torus, dictionary, substitution, mirror map, divisor equation, comparison with A1.2. Sentence: [Coa09] gives it on TF; Stage cannot be moved to TF because the product at 𝒬 = 0 is not semisimple there. No "Item O1". Not-adopted sources may be named in one remark with `Iri`, `CFK`.
7. A1.8 needs A1.7, which needs A3.2 and the first paragraph of A3.3: so **moved** forward: `lem:gen-mumford`, `lem:gen-ch`, then `lem:gen-free`, `prop:gen-teleman`.
8. A2.1-A2.7. In A2.1 the row-2 paragraph is the proof of C^c_{11}. Do not cite DR Remark 2.13.
9. `prop:gen-flat`: A2.8 in full: set-up, pushforward, defect, stability, equation. Keep the remarks: kappa terms; comparison with Teleman's Prop 7.7 (labelled not used); first-order sign check; PPZ (12) check; attribution. Remove "Rewritten in round 6c" and every mention of re-derivations.
10. A2.9, A3.1, rest of A3.3 (colour isolation, `thm:gen-const`), A3.4, A3.5 (uses s03), A3.6. Rank-3 remarks as one remark. Drop "Discharged for the draft".

Remove FB's parenthetical on the report's defects in A1.3. In A1.7 and A3.3 colour isolation appears twice; keep one proof and cite it.

Scope reminder for the writer: the T_max statement is an intermediate step; nothing is claimed about GW theory at generic weights beyond it.

## s05_specialisation

Covers: FB s.3 entire; Theorem FB_B statement from FB s.0.2. Part (f) becomes one neutral remark: the argument needs no rational tori, no limit along a family, no algebraic independence; no history.

Order: statement; B1 (`lem:spe-ring`), colours, reduction; B2; B3 with (Pair), 𝖥 free; B4 with chart, block, s_3-integral, recursion, conclusion; `rmk:spe-analytic`; `prop:spe-const` (constant term, uses `prop:gen-poly`); B6 (`lem:spe-numden`), with "why no limit is taken"; proof of (e): coefficient ring, `prop:spe-transc`. Rank-3 specifics as a remark. Part (b) statement refers to `prop:gen-poly` and `prop:spe-const`. Drop the paragraph "Against draft Thm 2.7" and the paragraph on draft Prop 2.11's φ_{k,m}. BFGW paragraph: one sentence, only for statements involving A.

Depends on: `thm:gen-main`, `prop:gen-osc`, `lem:gen-unit`, `lem:gen-rows`; `thm:row-pf`, `cor:row-unique`, `lem:row-pfint`, SP lemmas, `prop:row-sympl`, `prop:row-modrel`, `rmk:row-dict`.

## s06_frame

Covers: FB s.4 Notation, C Steps 0-6 and the statements of FB s.0.2: Theorem FB_C (C1), (C2), Definition (Cf), Lemma (t_cf), Prop 13.2_P2. DR Lemma 1.2(c) is the same lemma.

Order: notation; `lem:fra-setup`; `lem:fra-unfold`; `thm:fra-change`; `thm:fra-dict`; `lem:fra-LR`; `lem:fra-tcf` (**moved** before the Proposition, since the Proposition uses its logarithmic form; it is independent of s04, s05; constant part uses Barnes); `def:fra-cf`; `lem:fra-cfwell` (i)-(iv); `prop:fra-frame`, including the d > 0 form and the constant.

Remarks to keep, neutrally: the other branch of log y (statement in the lemma; it's part of it); A_cf = i√3(τ_ℍ - 1)A only if the writer wants; rule 7: do not mention "contradicts". "The round-5 critic" paragraph and C Step 2's sentence on the critic and "draft s.8.3" are dropped. C Steps 7-9 go to s10.

Dependencies: (C1) needs only `prop:row-sympl`, `lem:gen-mumford`, PPZ; this is worth saying. Identification with 𝒯_g needs `thm:spe-main`, `lem:spe-red`, `lem:spe-numden`. ev_A on 𝖠 needs [BFGW21, Lemma 4.2]; not used in the Proposition.

## s07_leaf

Covers: DR s.0.2 statement of Theorem LC_P2; s.3 (Remark E-scaling, s.3.1, s.3.2; s.3.3 reduced to the on-shell remark without [E]); s.4 except Lemma 4.1 (**moved** to s03); s.5; s.6.

Order: statement `thm:lea-main`; `rmk:lea-scale` (drop r_k evidence; if kept, in the appendix); s.3.1 list; `lem:lea-D0`; `lem:lea-conj`; `lem:lea-Dc`; `prop:lea-H2`; `rmk:lea-gauge`; `lem:lea-res`; `thm:lea-E00`; `thm:lea-M`; proof of the theorem (DR Cor 6.1 table in sentences); table of values at t_c = 0. DR's Prop 3.1 (closure for general c_1) is not used: omit and give only the first-order system at c_1 = 1, which is `eq:row-S1`. Remarks that refer to round 1, 3 and to scripts are dropped; their numeric content goes to the appendix.

Depends on: `lem:set-coord`, `cor:row-carry`, `cor:row-closed`, `prop:row-sympl`, `thm:row-pf`, `lem:row-closed12`. No BS input and no s04-s06. State this.

## s08_vanishing

Covers: LA s.0-s.7, plus the unstable-term conventions of s.1. From s.8 nothing is needed by the proofs. Optional complement, if a writer wants it: s.8.5 last bullets (the non-vanishing of Ξ_h(∅) is stated in s10 as a remark and is not used). LA s.7.4 and s.7.5 are the transfers: put the local-P2 one in s10 (`lem:reg-A`), the conifold one in s11. Local to this section: all integrals at scale 1 in 𝖰(τ), `\Xit`; specialise only in `cor:van-spec`.

Order: statement of `thm:van-main` and its corollaries; notation; LLZ statement `eq:van-LLZ` by [LLZ03, eq. (1)-(3)]; `prop:van-real`; second derivation through [OP04, Theorem 1] as a remark, with the caveat on the printed 1/|Aut μ|, wording neutral; scaled variables; s.3 lemmas; s.4 (F1)-(F6) `prop:van-bil`; `thm:van-Mk`; proof of the theorem; base case by hand as a remark; `cor:van-spec`, `cor:van-parts`, `lem:van-part01`, `cor:van-scale`.

Drop: LA s.9-s.12 as text (s.9 informs citations; s.12 goes to the appendix), provenance, dictionary to old reports.

## s09_value

Covers: DR s.7.6 statement of Theorem VI; s.7.8 Parts A, B. Part A uses `cor:van-scale` in place of "(7.2)" (they are the same formula at r = 0; `lem:reg-A`(a) in s10 cites the same corollary). Inputs: `prop:van-real` at μ = 1^a only, `lem:van-dil`, `eq:van-bil` at k = 2. It uses nothing from `thm:van-main`: say so. Conversion to R2 Cor 7 is dropped as a pointer to another route. Sanity limits: one remark. Status paragraph dropped; evidence to the appendix.

The notation of this section differs most from the source (𝖥 ⇒ ℱ_MV; D ⇒ ϑ; w ⇒ ϖ; ε ⇒ x̌; y ⇒ 𝗒). See NOTATION s.7.

## s10_regularity

Covers: DR s.7.1-s.7.7 and FB C Steps 7-9. Statements of Theorem R_P2 and Cor SV_P2 exactly as in DR s.0.2, with their hypotheses M_irr(h) for the stated range of h, followed by the sentence that `lem:reg-A`(b) gives M_irr(h) for every h ≥ 2, so that the hypothesis holds, the hypothesis being kept to show which h are used. Keep the scope: a statement about the decoration sum, per colouring. `thm:int-regularity` in s01 is the combination with `lem:reg-A`(b), and it names [LLZ03] as the input.

Order: `lem:reg-H0`; `lem:reg-string`, `lem:reg-C`, `lem:reg-L0`; `rmk:reg-nol0`; `lem:reg-range`; `lem:reg-A` (a) identification (b) vanishing (c) range; `thm:reg-main` (steps (i)-(iv)); `cor:reg-sv`, using `thm:val-main`; `cor:reg-gap`, `cor:reg-gapkh`, proofs from FB C Step 7, including the Taylor-coefficient paragraph and the Π_{γ_cf} paragraph from FB's statement; identification of ℱ_{0,g} remark (C Step 7 end, not used); `rmk:reg-compare` (C Step 8; DR "Comparison at g = 2" and "Consistency with BS", neutral); `rmk:reg-beta` (C Step 9); `rmk:reg-rank` (DR s.7.7, items (a)-(h), dropping (i) and pointers); `rmk:reg-dict`, optional, or move to s11. Non-vanishing of Z[h,0] does not specialise: part of `rmk:reg-rank`(f).

## s11_conifold

Covers: CP s.13 statements: s.13.0 (port), s.13.1, s.13.2 Corollary R with its input list in neutral form, Prop 13.2 statement of s.13.3, s.13.3 (5) specialisation of the vertex identity; LA s.7.4; DR s.7.0 dictionary if not in s10. Short. State results with their inputs, sketch where the methods of this paper apply. CP's Theorem LC (rank-2 leaf structure) has its own proof by a different route (Stirling, DLMF 5.11.1): state it as such with a pointer, not reproduced. \draftnote: full text to be decided. Use `\mathsf{x}` and `\mathsf{c} = 1 + \mathsf{x}` for x, c of CP, since x is the MV variable here. Do not import the repair history; state only "Teleman's inputs as in s04, rank 2". Item: non-degeneracy (F8).

## a_evidence

Covers: FB Appendix A.1-A.4 (A.5 dropped); DR Appendix; LA s.12; ST s.3.3 (gap at g ≥ 5). Opening sentence: nothing here is a proof step; arithmetic exact unless stated. No file paths or script names, no referee names. Describe each check: what was compared, at which weights and orders, controls that fail. Keep the limits of FB s.A.4: genus 2 sees only F_cf(0) at O(1); rank-one sign controls; R vs R^*; scale law; coverage statement. Note the sampling warning at α = (2,-3,7).

## Dropped from all sections

Round numbers, referee names and verdicts, tags [P] [C] [E] [S] [O], status tables (DR s.8, FB s.1.3-1.4, A.5, LA s.10-11), provenance, errata lists (FB s.5, DR s.9), user rulings, "proof of record", file paths, Lemma 1.4, code reading, mention of AI. Readings that differ from BS: NOTATION s.8.
