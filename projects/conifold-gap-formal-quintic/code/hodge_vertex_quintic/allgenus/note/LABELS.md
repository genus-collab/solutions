# Label registry (part 1 of 2: scheme, s01-s05)

Every named result of the note, its unique `\label`, and its source (digest item; primary source; gap note where it closes a step).
Revised after the gap phase (record: `design/DESIGN_CHANGES.md`). Part 2 (`LABELS_2.md`): s06-s10, appendices, equation labels. Use these
labels exactly. A writer adds auxiliary labels (equations, sublemmas) only with the prefix of its own section file and refers to another
section only through labels listed here. The "note name" is the name used in the text when the statement needs one; the source column
never appears in the note (rule 5). Rows marked **dropped (D8)** are not written (DEC D8); their labels must not be used.

**Prefixes.** s01 `int`, s02 `set`, s03 `gau`, s04 `tra`, s05 `red`, s06 `mod`, s07 `flh`, s08 `con`, s09 `gam`, s10 `mai`; appendix A `fsp`,
B `sgn`, C `evi`. Form `type:prefix-name`, types `thm`, `prop`, `lem`, `cor`, `defn` -> `def`, `rmk`, `conv`, `eq`, `sec`, `tab`. Exception
(DEC Q6): `def:gam-c` and `prop:gam-rational` keep their names but are placed in s08 (`sec:con-form`). Section labels (in the stubs):
`sec:int`, `sec:set`, `sec:gau`, `sec:tra`, `sec:red`, `sec:mod`, `sec:flh`, `sec:con`, `sec:gam`, `sec:mai`, `sec:fsp`, `sec:sgn`, `sec:evi`.
Named parts inside a statement keep the sources' item letters ((i), (ii), (a)-(d)); the input and hypothesis names (Q), (QRR), (M), (G0),
(T), (Sp), (A), (F_g), (G), (unit), (N0), (EP) are fixed by `OUTLINE.md` s.0.2 and `NOTATION_2.md` (never (S1)-(S3), REV R2).

## s01_intro (statements only; each proved later)

| note name | label | restates |
|---|---|---|
| Theorem A (graph-local conifold regularity, Property N(g) for all `g >= 2`) | `thm:int-main` | `thm:mai-main` |
| Corollary B (the conifold gap of the formal quintic) | `cor:int-gap` | `cor:mai-gap` |
| Theorem C (the normalisation `D = 1`; the actual gauge is canonical) | `thm:int-D` | `thm:gam-D`, `cor:mai-actual` |
| Theorem Gamma (the value at `q = 0` and at `q = oo`) | `thm:int-gamma` | `thm:gam-main` (ii), (iv) |
| Remark (the constant `c_g` and the quintic's conjectured constant) | `rmk:int-constant` | `rmk:mai-constant` |
| subsections | `sec:int-setting`, `sec:int-gap`, `sec:int-results`, `sec:int-inputs`, `sec:int-notclaimed`, `sec:int-strategy`, `sec:int-plan` | |

## s02_setup (subsection in brackets; REV_2 M14)

| note name | label | source |
|---|---|---|
| The formal quintic; invariants, `F_g` [fq] | `def:set-fq` | DA D1; LP18 s.0.2 (5); LP20 (5), (10); FQ2 s.1 |
| Genus-0 data, PF5, PF4 [fq] | `def:set-genus0` | DA D2; DB D2; FQ2 s.1; CYC s.0 |
| `PF5 = z^5 D PF4 - 1` [fq] | `lem:set-pf5` | DB S1; CYC Cor.; CF2 s.8.2 |
| The genus-0 mirror theorem (quoted, (M)) [fq] | `thm:set-mirror` | G1 s.1.1-s.1.4; GS s.2.1; MR1408320 Thm 11.6 |
| The genus-0 axioms of the twisted theory (quoted, (G0)) [fq] | `thm:set-G0` | G1_2 s.2.1; GJR18 s.4.1, s.5.2 (14); PPZ15 Def. 0.3 |
| The large-radius cyclic structure (flat-basis product) [fq] | `prop:set-cyclicLR` | G1_2 s.2.0, G1_3 s.2.8-s.2.9; REV R1 |
| Lemma: the Z/5 grading (G1 "Lemma A") [fq] | `lem:set-grading` | G1_2 s.2.2 |
| Lemma: unit and Frobenius property (G1 "Lemma B") [fq] | `lem:set-unitfrob` | G1_2 s.2.3 |
| Lemma: the Birkhoff recursion (G1 "Lemma C") [fq] | `lem:set-birkhoff` | G1_2 s.2.4 |
| Lemma: the periods of PF4 (G1 "Lemma D") [fq] | `lem:set-periods` | G1_2 s.2.5 |
| Lemma: the Wronskian (G1 "Lemma E") [fq] | `lem:set-wronskian` | G1_2 s.2.6 |
| Lemma: operators of order at most 3 (G1 "Lemma F") [fq] | `lem:set-order3` | G1_2 s.2.7 |
| Eigenvalues, `Delta_i`, `D u_i`; every QDE solution closes on PF5 [fq] | `cor:set-cyclicPF5` | G1_3 s.2.10 (a)-(b) |
| Points and variables `u, s, L, zhat` [fq] | `def:set-points` | DA D3; DD1 D1.1 |
| Frobenius data at the zeta point; the TFT `\TFT` [fq] | `def:set-frob` | DA D4; C2-C3; REV R6 |
| The R-matrix and its rows; `Phi_i`, `A_n` [rows] | `def:set-rows` | DA D6; S15 s.1.3; OU2 s.1; G1_4 s.3.1 |
| (QRR) (quoted) [rows] | `thm:set-qrr` | DD2 QRR-1; GJR18 (17); GJR17 (19)-(20); G1_4 s.3.1 |
| The twist sign [rows] | `rmk:set-qrrsign` | DD2 QRR-2; DA A5 |
| How (QRR) follows from the CohFT [rows] | `rmk:set-qrrderiv` | G1_4 s.3.2 with G1R_2 s.4; DEC Q1 |
| The vertex-level derivation of (QRR) | `rmk:set-qrrvertex` | **dropped** (subsumed by `rmk:set-qrrderiv`) |
| `h_(w_i)(z) = h_(w_0)(zeta^-i z)` [rows] | `lem:set-hrot` | DD1 cor:glu-colours; CP s.0.2 |
| The Givental-Teleman graph sum [graph] | `def:set-graphsum` | DA D9; PPZ15 Defs. 2.2, 2.5, 2.13; G1R_2 s.7 |
| The potential is the graph sum (quoted, (T)) [graph] | `thm:set-GT` | DA A6; MR2917177 Theorem (2); GJR18 Thm 4.1 |
| Conifold frames [frames] | `def:set-frames` | DA D11; DB D4 |
| The Frobenius basis at the conifold point; `gamma = 7/10` [frames] | `lem:set-frobbasis` | DB S2; CYC s.1; CF s.2.2-s.2.3; CA s.3 |
| The differential ring and the conifold expansion [frames] | `def:set-confexp` | DA D11; G5_2 s.2.2; G5R F2.1 |
| `tau(u)`, `v(u)`, `varsigma` [frames] | `def:set-tau` | DA D14; DB D7 |
| The canonical coordinates at the conifold point [frames] | `lem:set-canon` | DB D6; UR s.1; RN s.3(a) |
| The model `W_v`, its integrals and rows [model] | `def:set-model` | DB D8; DA D16, C9 |
| The model system and the T-system (uniqueness half) [model] | `prop:set-modrows` | DB S18 (a)-(b); UR s.1-s.2; PS s.1 |
| `M_0`, the model R-matrix in Teleman form [model] | `def:set-M0` | DB D15, D18; DD1 D1.5 |
| Symplecticity of `M_0` [model] | `lem:set-M0sympl` | G4_2 s.5 = G1_5 s.5.1; DEC Q3 |
| The model CohFT `Omega^mod_v`; the model TFT [model] | `def:set-modcohft` | DC1 D3; DC2 D6; S15 (I1); GS s.2.10 |
| Homogeneity of the model CohFT [model] | `lem:set-modhom` | G1_5 s.5.2; G1R_2 s.7 |
| Dilaton equation and homogeneity of the model potential [model] | `lem:set-dilaton` | G1_4 s.4; G1R_2 s.5; RD1 s.1.2; S15 Lemma 2 |
| The hatted theory [hat] | `def:set-hat` | DA D12; OU2 s.1; DC2 D14 |
| The outer graph sum; `T_Gamma`; `G_Gamma`; vertex graph [hat] | `def:set-outer` | DA D13, D18; DC2 D14, D16; REV R6 |
| Polar parts by reversion [hat] | `rmk:set-reversion` | DA A8; OU2 s.1; NF s.6 |
| Property N(g) [N] | `def:set-N` | DA A14; NF s.1 item 7 |
| The values of `c_g` [N] | `rmk:set-cg` | DA A16; S01 C2 |
| subsections | `sec:set-fq`, `sec:set-rows`, `sec:set-graph`, `sec:set-frames`, `sec:set-model`, `sec:set-hat`, `sec:set-N` | |

## s03_gauge

| note name | label | source |
|---|---|---|
| Objects of the cyclic-form theorem | `def:gau-conc` | DB D3, D17; CYC s.0 |
| Constants of `KK[ell]` | `lem:gau-const` | DB S3; CYC s.1; CF s.2.1 |
| PF4 is self-adjoint | `lem:gau-adj` | DB S4; CYC s.2; CF s.3 |
| The concomitant | `lem:gau-conc` | DB S5; CYC s.3; CF s.4 |
| Transversality | `lem:gau-transv` | DB S6; CYC s.4; CF s.5 |
| Isotropy (formal Picard-Lefschetz) | `lem:gau-iso` | DB S7; CYC s.5; CF2 s.6 |
| The conifold expansion is a ring map | `lem:gau-confexp` | FQ2 s.1; RS2 item 44; CA-K5; G5_2 s.2.2 |
| Special geometry | `lem:gau-sg` | DB S8; CYC s.6; CF2 s.7 |
| Theorem (the cyclic form of PF4, "(CYC)") | `thm:gau-cyc` | DB S9; CYC s.7; RS1 A1 |
| The cyclic system is PF5 | `cor:gau-cycsys` | DB S10; CYC Cor.; CF2 s.8.2 |
| Remarks on the cyclic form | `rmk:gau-cyc` | **dropped (D8)** |
| Cyclic systems with one simple pole; regular gauges | `def:gau-onepole` | DB D10-D11; UR s.3; GI s.1; G2 s.2 |
| The quintic cyclic system | `prop:gau-quintic` | DB D5, S12; UR s.1; PB s.2.8; GS s.2.9 |
| Equal products and `v(u)` | `prop:gau-equalprod` | DB S12; PS s.3; G2_3 s.5.1 (e) |
| Theorem U | `thm:gau-U` | DB S13; UR s.3; G2_3 s.5 |
| The commutant | `lem:gau-comm` | DB S14; UR s.3 (ii); G2_2 s.3 |
| Lifting | `lem:gau-lift` | DB S15; UR s.3 (iii); G2_2 s.4.1 |
| The scalar equation | `lem:gau-scal` | DB S16; UR s.3 (iv); G2_2 s.4.2 |
| The leading term of a solution of the gauge equation | `lem:gau-lead` | G2_3 s.5.4 (new) |
| The Z/5-degree of the gauge equation | `lem:gau-degree` | G2_4 s.6.1 (new) |
| Graded gauges | `prop:gau-graded` | DB S17; UR s.3-s.4; G2_4 s.6.2-s.6.4 |
| Remarks on Theorem U (what is used; the role of the pole) | `rmk:gau-Uremarks` | UR s.3 remarks; G2_3 s.5.6; G2R_2 s.3 R3 |
| The residue relation (R1) | `lem:gau-R1` | DB S21; UR s.5; GI s.3 |
| The diagonal of the gauge | `lem:gau-diag` | DB S27; GI s.5; GR s.2; G2_4 s.7 |
| Vanishing solutions | `def:gau-van` | DB D13; GI s.2 |
| The model's vanishing solution | `lem:gau-vanmod` | DB S22 (GI Lemma V) |
| Shape of `Y^C` | `lem:gau-vanshape` | DB S23 |
| The quintic's vanishing solution | `lem:gau-vanq` | DB S24 (GI Lemma Q) |
| Uniqueness of regular solutions | `lem:gau-uniq` | DB S25 (GI Lemma 1) |
| The vanishing cycle | `rmk:gau-monodromy` | **dropped (D8)** |
| The gauge identity (G) | `thm:gau-G` | DB S28; GI s.3; G2_4 s.8 (Cor. (a)) |
| The canonical gauge | `def:gau-can` | DB D12; G2_4 s.6.3 |
| Symplecticity | `prop:gau-symp` | DB S29; GI s.4; G2_3 s.5.4 |
| The dual identity | `cor:gau-Gstar` | DB S30 |
| Residue relations to all orders | `rmk:gau-res` | **dropped (D8)** |
| HPF | `lem:gau-HPF` | DD1 lem:pol-HPF; POLY s.1 |
| Power-series solutions of HPF are polynomial | `lem:gau-poly` | DD1 lem:pol-poly; POLY s.2 |
| Formal WKB uniqueness (every colour) | `lem:gau-wkb` | DD1 lem:glu-wkb; PS-C10; G2R_2 s.2.10, s.3 C2 |
| Colour covariance of the model columns | `lem:gau-colmod` | G2_5 s.9.3 (new) |
| Colour covariance of the actual rows | `lem:gau-colrows` | G2_5 s.9.4 (new); G2R_2 s.2.10 |
| The normalisation series `D` (WKB comparison, colour 0) | `thm:gau-R0` | DB S20; DD1 thm:glu-main (a)-(d); UR s.4; PS s.4 |
| The rows lie in the differential ring | `cor:gau-rowsNN` | G1_3 s.2.10 (c); G2_5 s.9.4 Step 6; REV R1 (new) |
| Colours of `D`; R0 in every colour | `cor:gau-colours` | DB S17 step 3; DD1 cor:glu-colours; G2_5 s.9.1-s.9.2 |
| (F_g) | `def:gau-Fg` | DD1 D1.8; ST2 s.2 |
| subsections | `sec:gau-cyc`, `sec:gau-systems`, `sec:gau-U`, `sec:gau-van`, `sec:gau-symp`, `sec:gau-R0` | |

## s04_transport

| note name | label | source |
|---|---|---|
| Loop spaces, the two settings, forms, insertions | `conv:tra-forms` | DC1 D1-D2, CV6; G3 s.2.1-s.2.2 |
| The summation lemma | `lem:tra-sum` | G3 s.2.2 (new) |
| The operator `Theta = 1 + N` | `lem:tra-theta` | G3_2 s.2 Step 1; G3R_2 K8 (new) |
| The vanishing tower | `def:tra-tower` | DC1 D7; DEC D1 |
| The slice `K` and `K_0` | `def:tra-slice` | DC1 D8; G3 s.3 item 6 |
| Frames (theory's "Lemma 0") | `lem:tra-frames` | DC1 A1; RD1 s.1.3; G3_2 s.1-s.2 |
| The affine slice | `cor:tra-slice` | DC1 A2; G3_2 s.2 Cor. (ii) |
| Coefficients of frames | `lem:tra-coeff` | DC1 A3; ST2 M5; G3_2 s.2 Step 5 |
| u-adic convergence | `lem:tra-conv` | DC1 O-3; G3_3 s.1-s.5 |
| Frame expansion of bivectors | `lem:tra-twoslot` | G3_2 s.3 (new) |
| Consistency checks | `rmk:tra-consistency` | DC1 A5 (a)-(c) |
| The isotropic sector | `rmk:tra-isotropic` | **dropped (D8)** |
| Descendant reading of the tower | `rmk:tra-descendant` | **dropped (D8)** |
| Dressings, (Sp), (A) | `def:tra-dressing` | DC1 D11 |
| Edges and the leaf point | `def:tra-edge` | DC1 D12; G3_4 s.2 Step 0 |
| Lemma E | `lem:tra-E` | DC1 B1; RD2 s.3.1; G3_4 s.1-s.2 |
| Corollary E' (edges in `K_0 (x) K_0`) | `cor:tra-edges` | DC1 B2; G3_2 s.3 |
| Lemma T (the leaf point) | `lem:tra-T` | DC1 B3 |
| Lemma F (finiteness) | `lem:tra-F` | DC1 B4 |
| Second proof by generating functions | `rmk:tra-K` | **dropped (D8)** |
| The canonical gauge as a dressing | `prop:tra-can` | DC1 B7; ST2 s.2 step 4 |
| subsections | `sec:tra-tower`, `sec:tra-dress`, `sec:tra-lemmas` | |

## s05_reduction

| note name | label | source |
|---|---|---|
| The flat sum `F^P_g`, flat term `F_Gamma` (homogeneous vertex tensor) | `def:red-flatsum` | DC1 D12; RD1 s.1.5; DEC Q5; G3 s.3 item 2 |
| The Model Lemma ML_K(h) | `def:red-ML` | DC1 D14; ST2 s.1 |
| Genus bound for vertices of edge graphs | `lem:red-genus` | DC1 O-3 (e); RD2 s.3.2; TRF s.2 |
| Theorem R | `thm:red-main` | DC1 C1; RD2 s.3.2; MSG2 s.7; G3_3 s.4 |
| Theorem R for the canonical gauge | `cor:red-can` | DC1 B7, C1; ST2 s.2 step 5 |
| What of the Model Lemma is used | `rmk:red-hyp` | DC1 C2 |
| The leaf bookkeeping | `rmk:red-leaf` | DC1 C3; RD2 s.3.3 |
| subsections | `sec:red-sum`, `sec:red-thm` | |
