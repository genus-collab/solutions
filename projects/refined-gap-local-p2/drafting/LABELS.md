# Label registry

Every named result of the paper, its unique `\label`, the section file that proves it and its source location. FB = FRAMEB_PROOF.md, DR = LP2_REGULARITY_PROOF_draft.md, LA = LEMMA_A_PROOF.md, CP = leaf_structure_proof.md. Use these labels exactly. A writer adds further labels (equations, auxiliary lemmas) only with the prefix of its own section file, and refers to another section's results only through labels listed here. The paper name is the name to use in the text if the statement needs a name; the source tag in the last column must not appear in the paper.

Prefixes: s01 `int`, s02 `set`, s03 `row`, s04 `gen`, s05 `spe`, s06 `fra`, s07 `lea`, s08 `van`, s09 `val`, s10 `reg`, s11 `con`, appendix `evi`. Form: `thm:row-pf`, `eq:row-rec`, and so on. Section labels: `sec:int`, `sec:set`, `sec:row`, `sec:gen`, `sec:spe`, `sec:fra`, `sec:lea`, `sec:van`, `sec:val`, `sec:reg`, `sec:con`, `sec:evi` (already in the stubs).

## s01_intro (statements only; each proved later)

| Paper name | Label | Restates |
|---|---|---|
| Theorem (frame formula) | `thm:int-frame` | `prop:fra-frame` |
| Theorem (regularity of graphs with an edge) | `thm:int-regularity` | `thm:reg-main`, with `cor:reg-sv` |
| Corollary (refined conifold gap) | `cor:int-gap` | `cor:reg-gap` |
| Corollary (coefficients) | `cor:int-gapkh` | `cor:reg-gapkh` |

## s02_setup

| Paper name | Label | Source |
|---|---|---|
| Convention on d = 0 | `rmk:set-d0` | FB s.0.2 (after graph sum), s.5.1 item 1 |
| Equivariant Picard-Fuchs equation, I-function | `eq:set-I`, `eq:set-pf` | BS eq:Ifun, eq:pf |
| Generators | `eq:set-gens` | BS eq:ADEFXdef |
| Pairing | `eq:set-eta` | BS eq:etaflat; DR s.1.3 on TF |
| Fundamental solution; small J | `eq:set-fundsol`, `eq:set-J` | BS eq:fundsol, eq:Jred |
| Superpotential | `eq:set-spot` | BS eq:spotP2 |
| Leaf identities (LAM1)-(LAM4) | `eq:set-lam1` ... `eq:set-lam4` | DR s.1.1 |
| Conifold coordinate and the unit ρ | `lem:set-coord` | DR Lemma 1.2(a), (b) |
| Recursion (Rec) | `eq:set-rec` | DR s.1.3 |
| Frame-B rows | `def:set-rows` | DR Def 1.3 |
| Edge factor | `eq:set-edge` | DR (1.1) |
| Decoration sum | `eq:set-graphsum`, `def:set-graphsum` | DR (0.1), s.0.1 |
| Remark on signs of square roots | `rmk:set-signs` | DR s.1.4, Remark (signs) |
| Hodge class, F^pot on 𝔇, Ξ_h | `def:set-fpot` | DR s.1.5 |
| The vanishing statement M_irr(h) | `def:set-mirr` | DR s.1.5 |
| Vertex factor | `eq:set-vertex` | DR (1.2) |
| Frame and mode remark | `rmk:set-frame` | DR s.1.6 |

## s03_rows

| Paper name | Label | Source |
|---|---|---|
| Picard-Fuchs for the closed form | `thm:row-pf` | DR Theorem A |
| Uniqueness of solutions | `cor:row-unique` | DR Cor 2.2 |
| Gaussian integration by parts | `lem:row-G1` | DR G1 |
| Total derivatives vanish | `lem:row-SP2` | DR SP2 |
| Fixed critical point | `lem:row-FP` | DR FP |
| Moving critical point | `lem:row-SP3` | DR SP3 |
| Coordinate invariance | `lem:row-SP1` | DR SP1 |
| Fubini | `lem:row-SP4` | DR SP4 |
| Formal Stirling | `lem:row-SP5` | DR SP5 |
| Critical points | `lem:row-crit` | DR Lemma 2.1 |
| Analyticity and the algebraic row Φ | `lem:row-analytic`, `eq:row-Phi` | DR Lemma 2.3; formula of Φ, Jet_m from FB Thm FB_B head |
| Picard-Fuchs under the integral | `lem:row-pfint` | DR Lemma 2.4 |
| WKB uniqueness | `lem:row-wkb` | DR Lemma 2.5, with its algebraic form |
| Closed forms of rows 1, 2 | `lem:row-closed12` | DR Lemma 4.1 (moved here: Cor 2.8 uses it) |
| Laurent form; carrying identities to the conifold | `cor:row-carry` | DR Cor 2.8 |
| Consequences of the closed form (F2), (L1), (Form), (R), (Z-infinity) | `cor:row-closed`, `eq:row-zinf` | DR Cor 2.9 |
| Symplecticity | `prop:row-sympl`, `eq:row-S1` | DR Prop 2.10 (S0-S6) |
| The relations A²D = X, 2E(y) = D - X | `prop:row-modrel` | DR Prop 2.12 |
| Dictionary at c_1 = A | `rmk:row-dict` | DR s.2.7, Dictionary |

## s04_generic

| Paper name | Label | Source |
|---|---|---|
| Main theorem at generic weights | `thm:gen-main`, parts (A1), (A2), (A3); `eq:gen-graphsum` | FB Theorem FB_A |
| Fixed loci | `lem:gen-fixed` | FB A1.1 |
| Weights move | `lem:gen-weights` | FB A1.2 |
| Polynomiality and regularity on the diagonal | `prop:gen-poly` | FB A1.3 with B5 (d ≥ 1, d = 0, GP for T, passage to TF; moved here) |
| CohFT axioms | `lem:gen-axioms` | FB A1.4 |
| Divisor equation | `lem:gen-divisor` | FB A1.5 |
| Semisimplicity | `lem:gen-ss` | FB A1.6 |
| Genus-0 differential equation | `prop:gen-qde` | FB s.1.2, second bullet |
| Mirror theorem at generic weights | `thm:gen-mirror` | FB s.1.2, first bullet, from [CCIT15, Theorem 31] |
| Hodge classes as an R-action | `lem:gen-mumford`, `eq:gen-mumford` | FB A3.2 (moved before A1.7) |
| Non-vanishing of ch_{2k-1} | `lem:gen-ch` | FB A3.3, first paragraph (moved before A1.7) |
| Freeness of the action | `lem:gen-free` | FB A1.7 |
| Existence and field of definition of R | `prop:gen-teleman` | FB A1.8 |
| Rows are formal | `lem:gen-rows` | FB A2.1 |
| Formal WKB | `lem:gen-wkb` | FB A2.3 |
| Comparison with the flat rows | `lem:gen-compare` | FB A2.4 |
| Y_j and Δ_j | `lem:gen-Y` | FB A2.5 |
| Unitarity | `lem:gen-unit` | FB A2.6 |
| Dictionary | `lem:gen-dict` | FB A2.7 |
| Flatness of R | `prop:gen-flat` | FB A2.8 |
| Ambiguity | `lem:gen-ambig` | FB A2.9 |
| Reduction at 𝒬 = 0 | `lem:gen-red` | FB A3.1 |
| The constants | `thm:gen-const` | FB A3.3, rest |
| Tensor structure | `lem:gen-tensor` | FB A3.4 |
| Oscillating integral | `prop:gen-osc` | FB A3.5 |
| Graph sum | proof ending at `eq:gen-graphsum` | FB A3.6 |

## s05_specialisation

| Paper name | Label | Source |
|---|---|---|
| Specialisation theorem | `thm:spe-main`, parts (a)-(e) | FB Theorem FB_B |
| The ring 𝒞 | `lem:spe-ring` | FB B1 |
| Row 0 at generic weights | `lem:spe-row0` | FB B2 |
| Rows 1, 2 and the pairing | `lem:spe-pair`, `eq:spe-pair` | FB B3, (Pair) |
| Reduction of the rows | `lem:spe-red` | FB B4 |
| Analytic corollary | `rmk:spe-analytic` | FB B4 |
| Constant term | `prop:spe-const` | FB B5, constant term |
| Numerator-denominator step | `lem:spe-numden` | FB B6 |
| Transcendence of F | `prop:spe-transc` | FB s.3, proof of (e) |

## s06_frame

| Paper name | Label | Source |
|---|---|---|
| Abstract set-up on TF | `lem:fra-setup` | FB C Step 0 |
| Unfolding | `lem:fra-unfold` | FB C Step 1 |
| Frame change | `thm:fra-change` | FB Theorem FB_C (C1), C Step 2 |
| Dictionary | `thm:fra-dict` | FB Theorem FB_C (C2) |
| Large-radius form | `lem:fra-LR` | FB C Step 3 |
| The conifold substitution | `def:fra-cf` | FB Definition (Cf) |
| Well-definedness | `lem:fra-cfwell` | FB C Step 4 |
| t_cf = t_c | `lem:fra-tcf`, with `eq:fra-log` for the logarithmic form | FB Lemma (t_cf), C Step 6; DR Lemma 1.2(c) |
| Frame formula | `prop:fra-frame` | FB Prop 13.2_P2, C Step 5 |

## s07_leaf

| Paper name | Label | Source |
|---|---|---|
| Leaf structure | `thm:lea-main`, parts (F1)-(F3), (L1)-(L3) | DR Theorem LC_P2 |
| Scaling | `rmk:lea-scale` | DR s.3 Remark (E-scaling) |
| Dilaton factors at colour 0 | `lem:lea-D0` | DR Lemma D0 |
| Conjugation by ρ | `lem:lea-conj`, (I1), (I2) | DR Lemma 4.2 |
| Recursion for coefficients | `lem:lea-Dc` | DR Lemma 4.3 |
| Structure of the singular row | `prop:lea-H2`, parts (a)-(d) | DR Prop H2 |
| Gauge remark | `rmk:lea-gauge` | DR Lemma G |
| Residue lemma | `lem:lea-res` | DR Lemma R |
| The edge at colour 0 | `thm:lea-E00`, parts (i)-(v) | DR Theorem E00 |
| Mixed and regular edges | `thm:lea-M`, parts (M-a)-(M-e) | DR Theorem M |
| Proof of the leaf structure | end of section | DR Cor 6.1 |

## s08_vanishing

| Paper name | Label | Source |
|---|---|---|
| Vanishing theorem | `thm:van-main` | LA Theorem L |
| Specialisation | `cor:van-spec` | LA Cor L1 |
| Parts ≥ 1 | `cor:van-parts` | LA Cor L2 |
| General scale | `cor:van-scale` | LA Cor L3 (DR (7.2)) |
| Mariño-Vafa formula | `eq:van-LLZ` | LA s.2.1 |
| Real form | `prop:van-real`, `eq:van-real` | LA Prop 2.1, (2.1) |
| Scaled variables | `eq:van-scaled`, `eq:van-G` | LA (2.2)-(2.4) |
| Operators | `lem:van-ops` | LA Lemma 3.1 |
| Dilaton leaves as unit legs | `lem:van-dil`, `eq:van-dil` | LA Lemma 3.2, (3.1) |
| Genus 0 | `lem:van-g0` | LA Lemma 3.3 |
| Leading coefficient | `lem:van-lead`, `eq:van-lead` | LA Lemma 3.4, (3.2) |
| String | `cor:van-string`, `eq:van-string` | LA Cor 3.5, (S) |
| Bilinear equations | `prop:van-bil`, `eq:van-bil` | LA s.4, (F1)-(F6), (4.2), (4.3) |
| Relations at the dilaton point | `thm:van-Mk`, `eq:van-Mk` | LA Theorem M_k, (5.1) |
| Parts 1 and 0 | `lem:van-part01` | LA s.7.2 |

## s09_value

| Paper name | Label | Source |
|---|---|---|
| Vertex identity | `thm:val-main` | DR Theorem VI |
| Vertex identity at scale 1 | `thm:val-tau` | DR Theorem VI_τ |
| Closed equation on the p_1 slice | `eq:val-star` | DR s.7.8 (*) |
| Equation for the Z_h | `eq:val-E1` | DR s.7.8 (E1) |

## s10_regularity

| Paper name | Label | Source |
|---|---|---|
| Dilaton homogeneity | `lem:reg-H0`, `eq:reg-H0` | DR s.7.1 (H0), (7.1) |
| String identity | `lem:reg-string` | DR s.7.1 (S_P2) |
| Vanishing on 1^p | `lem:reg-C` | DR s.7.1 (C_P2) |
| The **λ**_0 identity | `lem:reg-L0` | DR s.7.1 (L0_P2) |
| Absence of **λ**_0 | `rmk:reg-nol0` | DR s.7.2 |
| Range | `lem:reg-range` | DR Lemma 7.1 |
| Transfer of the vanishing | `lem:reg-A` | DR Lemma A_P2 |
| Regularity theorem | `thm:reg-main` | DR Theorem R_P2, s.7.5 |
| The single-vertex graph | `cor:reg-sv` | DR Cor SV_P2 |
| Refined conifold gap | `cor:reg-gap` | FB Cor G, C Step 7 |
| Coefficients | `cor:reg-gapkh` | FB Cor G_{k,h}, C Step 7 |
| Comparison at low genus and limits | `rmk:reg-compare` | FB C Step 8, DR s.7.6 |
| Graph sum for β | `rmk:reg-beta` | FB C Step 9 |
| Rank and parameter audit | `rmk:reg-rank` | DR s.7.7 |
| Comparison with the conifold (dictionary) | `rmk:reg-dict` | DR s.7.0 |

## s11_conifold

| Paper name | Label | Source |
|---|---|---|
| Frame formula for the conifold | `prop:con-frame` | CP Prop 13.2 |
| Regularity for the conifold | `cor:con-reg` | CP Cor R |
| Range (conifold) | `lem:con-range` | CP Lemma 13.0 |

## a_evidence

Subsections `sec:evi-gen` (FB s.A.1), `sec:evi-spe` (FB s.A.2), `sec:evi-fra` (FB s.A.3), `sec:evi-rows` (DR Appendix), `sec:evi-van` (LA s.12), `sec:evi-limits` (FB s.A.4; ST s.3.3 coverage).

## Bibliography keys (refs.bib)

| Key | Work | Pinned numbers (form: `\cite[Theorem 31]{CCIT15}`) |
|---|---|---|
| BS | Brini-Schuler, revised paper | by statement, not number |
| CCIT15 | Coates-Corti-Iritani-Tseng | Theorem 31, Definition 28, eq. (19), eq. (4), Remark 3; s.2.2, s.3.1; proof structure Theorem 41, Props 45-47, Theorem 2 |
| Coa09 | Coates | Prop 4.1, Cor 4.2 (covers TF) |
| Tel12 | Teleman | Theorem (2) of the Introduction; Prop 8.7(ii); 2.17, 2.18, 2.19; (2.15). In remarks only: Def 7.1, (7.5), (7.6), Prop 7.7, (7.7.a), Ex 7.17, Prop 8.2(ii), Prop 8.7(i) |
| PPZ15 | Pandharipande-Pixton-Zvonkine | Defs 2.2, 2.5, 2.13; Prop 2.14; (5); s.3.2. In remarks only: Props 2.8, 2.9, eq. (12), footnote 10, Remark 1.6 |
| GP99 | Graber-Pandharipande | formula (1), Apps A, B, C |
| Mum83 | Mumford | (5.2), Cor (5.3), (5.4), (5.5) |
| FP00 | Faber-Pandharipande | Theorem (Mumford) s.1.1; Theorem 4; s.2 before (22) |
| LLZ03 | Liu-Liu-Zhou | eq. (1)-(3), called there the Mariño-Vafa formula |
| OP04 | Okounkov-Pandharipande | Theorem 1 (second published proof) |
| Pan98 | Pandharipande | Lemma 1; s.1.2 (6), s.1.3 Prop 2 |
| CG07 | Coates-Givental | App. 1 (iv), (v), (vi) |
| Giv01 | Givental | (3)-(6), Proposition (a)-(d), s.2.3 Lemma (comparison) |
| Giv03 | Givental | Theorem 1 (consistency only) |
| Loo96, MW07, Har85, Iva90 | stability | Looijenga Props (2.1), (2.2); Madsen-Weiss Theorem 1.1.1. Harer and Ivanov as reported by Looijenga s.1 |
| Wah, RW | ranges | Wahl s.1, Theorems 1.1, 1.2; Randal-Williams s.1.3 (comparison of range only) |
| DLMF | | Eq. 15.6.7; 5.11.1, s.5.11(ii) (conifold) |
| Bar10, Bai35, AAR99, Jan | Barnes's second lemma | Bailey s.6.2; AAR Theorem 2.4.3; Jantzen eq. (2) |
| BFGW21 | Bousseau-Fan-Guo-Wu | Lemma 4.2, Conj 1.9, s.4.1 |
| LP18 | Lho-Pandharipande | (19), (24), (27), (29), (30), (36), Lemma 8 (comparison only) |
| Liu13 | Liu | eq. (51) |
| Mac95, Sta99 | textbook | Macdonald I.3 Ex. 2 and pp. 10-11, Stanley Cor 7.21.3: numbers from memory, mark with \draftnote |
| CFK, Iri | not used, only in a remark | CFK Theorem 1.2.2, Cors 1.4.1, 1.4.2; Iritani Prop 2.4 |
| JM83, HO15 | corroboration of the bilinear equations, optional | Jimbo-Miwa Appendix 1 p. 996; Harnad-Orlov (2.2) |
| HK10, CKK14, HKR08, CI18, ZZ08 | background in the introduction | none pinned |

## Labels added while writing the sections

These are defined in the section files in addition to the plan above.

- `a_evidence`: `sec:evi-con`
- `s01_intro`: `sec:int-inputs`, `sec:int-notclaimed`, `sec:int-plan`, `sec:int-problem`, `sec:int-statements`, `sec:int-strategy`
- `s02_setup`: `eq:set-DeltaTF`, `eq:set-Gsum`, `eq:set-Hg`, `eq:set-Ian`, `eq:set-LW`, `eq:set-S`, `eq:set-WTF`, `eq:set-asym`, `eq:set-colours`, `eq:set-dil`, `eq:set-etaTF`, `eq:set-hw`, `eq:set-lam2`, `eq:set-lam3`, `eq:set-qprod`, `eq:set-series`, `rmk:set-norm`, `sec:set-colours`, `sec:set-coord`, `sec:set-factors`, `sec:set-givental`, `sec:set-graphsum`, `sec:set-hodge`, `sec:set-leaves`, `sec:set-mirror`, `sec:set-rows`, `sec:set-special`, `sec:set-target`
- `s03_rows`: `eq:row-DEX`, `eq:row-conv`, `eq:row-edgedelta`, `eq:row-pfcoef`, `eq:row-thetav`, `eq:row-vS`, `sec:row-closed12`, `sec:row-conseq`, `sec:row-crit`, `sec:row-modrel`, `sec:row-pf`, `sec:row-sp`, `sec:row-wkb`
- `s04_generic`: `eq:gen-rowrec`, `rmk:gen-attr`, `rmk:gen-deps`, `rmk:gen-kappa`, `rmk:gen-mirror`, `rmk:gen-ppz12`, `rmk:gen-rank3`, `rmk:gen-sign`, `rmk:gen-teleman77`
- `s05_specialisation`: `eq:spe-identity`, `rmk:spe-cfgraph`, `rmk:spe-noinput`, `rmk:spe-nolimit`, `rmk:spe-rank3`, `sec:spe-statement`
- `s06_frame`: `eq:fra-pair`, `rmk:fra-branch`, `rmk:fra-cfread`, `rmk:fra-leaves`
- `s07_leaf`: `eq:lea-Et`, `eq:lea-K1`, `eq:lea-Zinf`, `eq:lea-gaugecomb`, `eq:lea-v`, `eq:lea-vexp`, `rmk:lea-onshell`, `sec:lea-proof`, `sec:lea-row0`
- `s08_vanishing`: `eq:van-Gd`, `eq:van-LLZV`, `eq:van-LLZdef`, `eq:van-Xi`, `eq:van-bilc`, `eq:van-bilid`, `eq:van-block`, `eq:van-exp`, `eq:van-scale`, `rmk:van-scale`, `rmk:van-table`, `sec:van-MV`, `sec:van-Mk`, `sec:van-bil`, `sec:van-cors`, `sec:van-dil`, `sec:van-notation`, `sec:van-proof`, `sec:van-statements`
- `s09_value`: `eq:val-beta`, `eq:val-slice`, `sec:val-proof`
- `s10_regularity`: `eq:reg-scale`, `rmk:reg-unref`, `sec:reg-gap`
- `s11_conifold`: `eq:con-Delta`, `eq:con-frame`, `eq:con-value`, `lem:con-norm`, `thm:con-leaf`
