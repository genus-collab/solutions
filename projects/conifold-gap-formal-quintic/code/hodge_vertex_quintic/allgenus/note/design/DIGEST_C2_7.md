# DIGEST_C2, part 7: s.6 Literature used; s.7 Open issues for the writers

Continues `DIGEST_C2_6.md`.

## s.6 Literature (statement numbers; what is used; whether the sources read the statement in the full text)

| arXiv | statements | used for | read in full text? | root refs.bib key |
|---|---|---|---|---|
| 1509.06954 (DNOPS) | Thm 8.1 (Hurwitz case, superpotential (8-1)); (2-26)-(2-29) TR conventions; (2-30)-(2-35) dictionary ((2-31) Laplace transform of `B`, (2-32) of `dy`, (2-34) `xi`, (2-35) the expansion of `omega_{g,k}`); (2-7) metric; Thm 5.1 (general Dubrovin curve) | `thm:trm-ident` (M1) | yes: TRM s.1 "read s.1-8 from the PDF" (v2); hypotheses checked there (TRM s.1) | `Dunin-Barkowski:2015caa` (root entry: J. Inst. Math. Jussieu, year 2017, pp. 1-49; TRM: 18 (2019); pin it) |
| 2212.00320 (x-y swap) | Thm 7.1, closed formula (7.9), example (7.2) (`r = 5`), Remark 7.3; Thm 1.8 / 1.14; kernel (1.1); (1.13) | `thm:mlk-lreg` (M7); Gap route 2; `lem:mlk-vanper` | yes: XY s.1 (s.1, s.6-s.7); XYR s.2 (s.1.1-1.5, s.7, pp. 2-12, 46-49), "statement and hypotheses checked". LR F6(b), DR s.2: number not re-checked by them | none: new key |
| math-ph/0702045 (EO) | Thm 5.1, (5-14) (`g + k > 1`), (5-15) (`g >= 1`); s.5.3 temperatures (5-27)-(5-31); (5-1); Rauch (5-4); (5-13) with footnote 8; Def. 4.2 (4-6)-(4-9) and the kernel (4-3), footnote 3; Def. 4.3 (4-10); Thm 4.1 (4-20); Thm 4.2; Cor. 4.1 (4-26)-(4-29); Thm 4.4 (4-23)-(4-24) (linear loop equation); Thm 4.5 (4-25); Thm 4.6; Thm 4.7 (4-30) (dilaton) | (Q1) `thm:mlk-var`; (F1) properties; M6 `lem:trm-xfun`; `rmk:mlk-dilaton` | yes: LR consulted printed pp. 20-36 (s.4.1-5.3). The linear loop equation's number was never pinned (SYN2 M6). [DW] read pp. 21-25 and 31-36: Thm 4.4 is the linear loop equation; Thm 5.1 and s.5.3 as LR quotes them | `Eynard:2007kz` |
| 2502.18115 (Hock) | Cor. 1.2, formula (1.3); Thm 1.1 / 3.2; Prop. 3.4-3.5; the example of s.3.2 has a sign slip | Gap route 1 (optional) | yes: TRM s.1 (v2: s.1-3.2, Prop. 3.4-3.5); XYR pp. 1-5; LR abstract only (F5: convention not checked by LR) | none: new key |
| 1303.1043 (R-matrix action) | (5) symplectic condition; Remark 2.1 and (6) (the edge series is a power series); Def. 2.2 (R-matrix action: vertices, legs `R^-1(psi)`, edges); Prop. 2.3; Prop. 2.4 (left group action; red/blue proof); Def. 2.5 (translation); Props. 2.7, 2.8; Prop. 2.9 (`T_a = R T_b` implies `T_a R Omega = R T_b Omega`); Remark 2.10; Prop. 2.12 (unit-preserving combination, `T_a = z[R(1) - 1]`, `T_b = z[1 - R^-1(1)]`) | chain step 4 (Def. 2.2); S15 (Prop. 2.4, Prop. 2.9; Def. 2.5 in Lemma 4; Remark 2.1, (6) for `E_O`) | yes: FR located Def. 2.2, Prop. 2.4, Def. 2.5, Prop. 2.9 in the PDF; [DW] read pp. 19-27 of the current arXiv PDF (version not recorded) | none: new key |
| 0712.0160 (Teleman) | (no number in Stage 2) | background to (I1) | FR: title and abstract only | `MR2917177` |
| 1211.4021 (DOSS) | source of DNOPS (2-30)-(2-35) | background (TRM s.1) | not stated | `DuninBarkowski:2012bw` |
| math/0108100 (Givental) s.5 | ancestor-descendant relation (no number) | remark only, not used (LL s.1; LR F6(c)) | no | none |
| 2312.16950 (log TR) | abstract | candidate for a (Q1)-free proof (SYN2b s.6.4); not used | abstract (XY s.1) | none |
| 2211.12259 | metric (1.6), s.1-2.1 | sign comparison of metrics (TRM s.1) | TRM: s.1-2.1 | none |
| 2211.08917, 2304.03032, 2408.02608 | abstracts | not needed (XY s.1) | abstracts | none |

Topological recursion is applied ONLY to the model curve (D1 and the families of XY); nothing in this digest uses the quintic's own
mirror curve (SYN2 s.1 "nothing about the quintic enters"; SYN5 s.1; RS2 s.6.2).

## s.7 Open issues for the writers

### Gaps against the binding rules (each needs a `\draftnote` or an author decision)
* **G1 (rule 3). Theorem V(ii), base `h = 1`.** LL2 s.3 proves `<-log p>_1 = -1/(12v)` by homogeneity plus ONE exact computation over
  `Q(zeta5)`. A written alternative exists in the sources: `lem:mlk-vanper` at `g = 1` (XY s.2.3, PROVED mod arXiv:2212.00320 Thm 7.1)
  with the orientation of XYR s.2. No source combines the two; `\draftnote` in s.6 of the note. (Stage 1 had `h = 0, 1` VERIFIED only,
  S01 s.4; Stage 2 made `h = 0` PROVED, LL s.2.4.)
* **G2 (rule 3). M1 to all orders.** TRM s.2.3 gets `R_TR = M_0` from the `(U,V)`-uniqueness plus the agreement of `R_1`, the latter
  from the `z^8` check (VERIFIED); the sign of `kappa` is "ARGUED". A hand derivation of `R_1` (or `V`) on both sides is missing.
  Stage-1 material (TRM); coordinate with the digest that covers it; `\draftnote` in s.5.
* **G3 (rule 2). M6 (functions of x pair to zero).** Only "standard; linear loop equation" in the sources. [DW] suggestion for the authors,
  NOT in the sources: by EO Thm 4.4 (4-23) the sum over the preimages of `x` of `omega/dx` vanishes (stable `omega`), and the only
  points of P^1 where `f(x) omega` can have residues besides the `c_b` are `p = 0, oo`, the preimages of `x = oo`; so
  `Sum_b Res_{c_b} f(x) omega = -Res_{x=oo} f(x) Sum_i omega(p^i)/dx(p^i) dx = 0`. Needs an author's check; `\draftnote`.
* **G4 (rules 3, 5). S15 Lemma 3 is a statement about the code** (FR). The note must define the hatted outer sum mathematically (rows
  `Qhat_a = s^(1-c_a) Q_a` at `zhat = z/s^4`, `xhat = x/s`, `x = s varpi_0`, leaves from row 0) from NATURAL_FORM s.5-s.6 / UR (not read
  here) and prove Lemma 3 for that definition; "pinned down by the genus-2/3 matches" cannot be cited.
* **G5. Which definition of K.** MSG s.1 defines `Pi_i` through the calibration column `S(z)H^2` and identifies
  `Y^C_b = -(1/5) eta(S(z)H^2, H^b)` by code; LR did NOT check LL s.1's calibration paragraph. RD s.1.3 defines `T_v` directly from `Y^C`.
  [DW]: `V^(k) = (-1)^(k-1) Pi_(k-1)` on the basis, so `K`, `K_0` can be defined from `Y^C` alone (then D(ii) needs only GI Lemma V);
  the descendant reading of `Pi_i` becomes a remark. Author decision; `\draftnote` in s.3 or s.5.
* **G6. LD has no separate referee** (proved by the critic in SYN2 s.2 step 7.2; VERIFIED). First item of the human-check list.
* **G7. SKETCH items** to keep out of proofs: the hatted-engine form of the universality corollary (SYN2 s.2); BV s.1's converse
  paragraph (O24); `rmk:mlk-dilaton`'s Euler-dilaton action; TRM's sign-of-`kappa` argument.
* **G8.** Use LR F2 for the base `h >= 2` of Theorem V(ii) (no (Q1)); (Q1) is still needed for (i) and for the steps `k -> k+1`.

### Writing details
* **W1.** State the kernel convention and the dictionary once (appendix A); cite EO Thm 5.1 with s.5.3 (5-27)-(5-30), not "s.5.1" (O7, F1);
  say that the engine `omega_{0,3}` sign is opposite to EO Thm 4.1 (O13).
* **W2.** State the convention of (Q3) (F5, O13).
* **W3.** Use (F_g) (`D == 1 mod z^(3g-2)`) wherever S15 says (F) (K9, O20).
* **W4.** S15 slips: leaf formula `T' = z y + q_{R'}` (O18); define `x`; reference (I1) for homogeneity; `d(0) = 1` (u-adic convergence) (C5).
* **W5.** Transfer from `x = z^4 + V/z` (XY), `x = P(z) + v/(z-a)` (TRM s.4.1), `lambda = p^4 + v/p` (DNOPS) to `x = W` (D1): scaling of
  TR under `x -> kappa x` (TRM s.2.4 gives it for the CohFT side); not written for the closed formula and the leg classes. One paragraph.
* **W6.** Evidence counts: 108 + 4 (not 112; K4, O14); X4 "fails", no count (K5, O16); F1's 168 = 164 zeros + 4 (C4, O17); "16200 vectors of
  the raw-leg kernel" (O1); E1 of LL is not evidence for (Q1) (O8); "functions of x pair to zero" is a consistency check (O12, F4).
* **W7.** LL2: "(h,n) != (0,2)" is implied; write the bases of Theorem V(i) for `(h, a) = (0,1), (0,2)` as the explicit (0,3) pairings.
* **W8.** The (Q) list: SYN2 s.0 omits arXiv:1303.1043 Def. 2.5; SYN2 s.2, MS2 s.1, MS5 s.2 include it. Use the full list; add PPZ
  Prop. 2.12 if the unit-preserving form is cited (L3).
* **W9.** LR s.2 writes `Omega = -dp/p` in its hypothesis list but `+dp/p` in its Rauch check; `+dp/p` is right ([DW], EO (5-1)); harmless.
* **W10.** Erratum to this digest: in `DIGEST_C2_4.md`, `thm:trm-ident` step 1, "`dig2` entry" means the hypotheses listed under
  `thm:trm-ident` in `DIGEST_C2_2.md`.
* **W11.** Notation clashes: s.4 items C6, C7, C13.

### Literature to pin
* **L1.** EO's linear loop equation: Thm 4.4, (4-23)-(4-24) ([DW]); cite it for M6.
* **L2.** arXiv:2212.00320 Thm 7.1: number confirmed by XYR (Stage 1) only.
* **L3.** The translation `T(z) = z(1 - R^-1(z) 1)` of RD s.1.5 and S15 s.1.1 is `T_b` of PPZ Prop. 2.12 ([DW]); the sources cite "s.2".
* **L4.** Hock Cor. 1.2: convention through TRM s.1 only.
* **L5.** New BibTeX keys for arXiv:2212.00320, 2502.18115, 1303.1043 (and math/0108100 if the remark is kept); reuse
  `Dunin-Barkowski:2015caa`, `Eynard:2007kz`, `MR2917177`, `DuninBarkowski:2012bw`.

### Adjudicated conflicts (all settled in SYN2b s.4; apply silently, list in DRAFT_NOTES)
K1 kernel orientation (C1); K2 sign of the xi-term (`+`); K3 Cor. 4 uses D(ii); K4 108 + 4; K5 X4 count orientation-dependent;
K6 `|c0| = m^10/|f(0) g(0)|`; K7 `B` is a coordinate (Lemma C); K8 s.1.5 is not an input; K9 window `z^(3g-3)`; K10 LD of record;
K11 D(iii) PROVED. Overclaims O1-O27 (SYN2b s.5): applied in s.2 of this digest where they touch a statement (O1-O7, O9-O11, O13,
O18, O20-O27); O8, O12, O14-O17, O19 concern evidence (W6).

### Human-check list of the sources (for DRAFT_NOTES only; SYN2b s.6.2, still open per MS5 s.5 item 1.7)
1. LD, D(ii) (DR s.3), D(iii) (DP s.2-s.5). 2. Theorem V with the sign dictionary (EO Thm 5.1; differentiation under the contours).
3. Quintic Theorem R and theory Lemmas 0, E, T, F. 4. Theorem S15, Lemma 3 against the definition of the hatted outer sum.
5. Theorem U, Theorem G, Prop. S. 6. Hypotheses of the quoted theorems (DNOPS Thm 8.1 case (8-1), `n = 4`, `a = (0,0,0,0,v)`;
ABDKS Thm 7.1 with `y = z`, rational `x`, simple critical points; EO Thm 5.1; Hock Cor. 1.2; PPZ s.2).

### Suggested placement (write-up plan SYN2b s.6.3, updated by MS5 s.5 item 2)
* s.5 (model as TR): D1-D5, D7-D10; `thm:trm-ident`, `lem:trm-xfun`, `prop:trm-legdict`, `lem:trm-lemG`, `lem:trm-resid`, `prop:trm-dii`,
  DP Lemmas 0-3, `thm:trm-lattice`, `cor:trm-split`, `cor:trm-Ksplit`, `prop:trm-coeff`, `rmk:trm-h2leg`.
* s.6: `thm:mlk-lreg`, `thm:mlk-gap`, `lem:mlk-vanper`, `lem:mlk-fixedx`, `lem:mlk-ibp`, `thm:mlk-var`, `lem:mlk-genus0`, `prop:mlk-V1`,
  `thm:mlk-logleg`, `thm:mlk-main`, `rmk:mlk-dilaton`.
* s.7: D12-D15; S15 Lemmas 1-4, (I1)-(I4), `thm:flh-main`, `rmk:flh-fail`.
* s.8: the chain (s.2.5), merged with the Stage 4-5 discharge of (F_g); the final conditionality (Q), (QRR), with (U) proved.
* s.9: `lem:uni-pole`, `lem:uni-coord`, `cor:uni-flat`; the failing controls as necessity examples (appendix B).
* Appendix A (signs): s.4 C1-C13. Appendix B (evidence): s.5 E1-E20. Appendix C (errata to the working notes) is project
  bookkeeping and conflicts with binding rule 5: keep it in DRAFT_NOTES, not in the note.
