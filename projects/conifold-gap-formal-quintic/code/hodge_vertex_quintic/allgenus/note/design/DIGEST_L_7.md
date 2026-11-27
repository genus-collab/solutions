# DIGEST_L_7: s.3.5 entries B; s.5 evidence; s.6 literature used; s.7 open issues

Continuation of `DIGEST_L_6.md`.

## 3.5 Entries, B (topological recursion, Gamma factors, DLMF; optional appendix references)

```bibtex
@article{ABDKS25,
  author = {Alexandrov, Alexander and Bychkov, Boris and Dunin-Barkowski, Petr and Kazarian, Maxim and Shadrin, Sergey},
  title = {A universal formula for the {$x$}--{$y$} swap in topological recursion},
  journal = {J. Eur. Math. Soc. (JEMS)}, year = {2025}, doi = {10.4171/JEMS/1615},
  eprint = {2212.00320}, archivePrefix = {arXiv}, primaryClass = {math-ph},
  note = {arXiv:2212.00320v3; numbers cited from this version. Volume and pages to be added (the arXiv page gives journal,
          DOI and year only)}}

@article{Hock25,
  author = {Hock, Alexander},
  title = {Symplectic (non-)invariance of the free energy in topological recursion},
  journal = {Comm. Math. Phys.}, volume = {406}, year = {2025}, pages = {192}, doi = {10.1007/s00220-025-05373-8},
  eprint = {2502.18115}, archivePrefix = {arXiv}, primaryClass = {math-ph},
  note = {arXiv:2502.18115v2; numbers cited from this version; 192 is the article number}}

@article{CvdLPS23,
  author = {Carlet, G. and van de Leur, J. and Posthuma, H. and Shadrin, S.},
  title = {Enumeration of hypermaps and {H}irota equations for extended rationally constrained {KP}},
  journal = {Commun. Number Theory Phys.}, volume = {17}, number = {3}, year = {2023}, pages = {643--708},
  doi = {10.4310/CNTP.2023.v17.n3.a3},
  eprint = {2211.12259}, archivePrefix = {arXiv}, primaryClass = {math-ph}, note = {arXiv:2211.12259v1. Optional}}

@article{GGI16,
  author = {Galkin, Sergey and Golyshev, Vasily and Iritani, Hiroshi},
  title = {Gamma classes and quantum cohomology of {F}ano manifolds: {G}amma conjectures},
  journal = {Duke Math. J.}, volume = {165}, number = {11}, year = {2016}, pages = {2005--2077},
  doi = {10.1215/00127094-3476593},
  eprint = {1404.6407}, archivePrefix = {arXiv}, primaryClass = {math.AG},
  note = {arXiv:1404.6407v4; numbers cited from this version}}

@misc{DLMF,
  key = {DLMF},
  title = {{NIST} Digital Library of Mathematical Functions},
  howpublished = {\url{https://dlmf.nist.gov/}},
  note = {Version 1.2.8, release date 2026-09-15. Cited: Eqs. 2.3.8, 2.3.13--2.3.18, \S2.3(iv), 5.11.1, 5.11.8, 24.4.1}}

@misc{Lei24,
  author = {Lei, Patrick},
  title = {Higher-genus {G}romov--{W}itten theory of one-parameter {C}alabi--{Y}au threefolds {I}: {P}olynomiality},
  year = {2024}, eprint = {2409.11659}, archivePrefix = {arXiv}, primaryClass = {math.AG}, note = {arXiv:2409.11659v2. Optional}}

@misc{Eti24,
  author = {Etingof, Pavel},
  title = {Mathematical ideas and notions of quantum field theory},
  year = {2024}, eprint = {2409.03117}, archivePrefix = {arXiv}, primaryClass = {math-ph},
  note = {arXiv:2409.03117v2 (2026). Optional (appendix context)}}

@article{Kir10,
  author = {Kirwin, William D.},
  title = {Higher asymptotics of {L}aplace's approximation},
  journal = {Asymptot. Anal.}, volume = {70}, number = {3-4}, year = {2010}, pages = {231--248},
  eprint = {0810.1700}, archivePrefix = {arXiv}, primaryClass = {math.CA}, note = {arXiv:0810.1700v2. Optional}}

@incollection{Pol05,
  author = {Polyak, Michael},
  title = {Feynman diagrams for pedestrians and mathematicians},
  booktitle = {Graphs and patterns in mathematics and theoretical physics}, series = {Proc. Sympos. Pure Math.}, volume = {73},
  publisher = {Amer. Math. Soc., Providence, RI}, year = {2005}, pages = {15--42},
  eprint = {math/0406251}, archivePrefix = {arXiv}, primaryClass = {math.GT},
  note = {arXiv:math/0406251v1. Volume and pages: arXiv journal reference; book title as in the root entry MR2131017 of the
          same volume. Optional}}

@misc{BGRT,
  author = {Bar-Natan, Dror and Garoufalidis, Stavros and Rozansky, Lev and Thurston, Dylan},
  title = {The {A}arhus integral of rational homology 3-spheres {II}: {I}nvariance and universality},
  year = {1998}, eprint = {math/9801049}, archivePrefix = {arXiv}, primaryClass = {math.QA},
  note = {arXiv:math/9801049v4; journal version to be added. Optional}}
```
Not on arXiv, only if GJR18's proof of Lemma 7.7 is retold: Hermite, J. Reine Angew. Math. 115 (1895) 201-208, and Nemes, Integral
Transforms Spec. Funct. 24(8) (2013) 597-606 (both as listed in GJR18's references [22], [30]; not read). DLMF 5.11.8 suffices.

## 5. Evidence for the appendix (literature side only)
* I ran no mathematics code. The one literature-side exact check is the project's own: `stage5/lit/lit_vertex_check.py` (LIT s.1.5):
  from LPfq (26) at `lambda_i = zeta^i`, `N_(2,0) = -19/9600` and `N_(3,0) = 1/224000` (the project's localisation values), and the
  control with twist weight `+5` gives `-97/28800`, `1/81000`; exact over `Q(zeta_5)`; re-run identical (MECHANISM_STAGE5 s.4.1).
  It supports the sign of the twist weight in (QRR) (LIT s.1.4); it is evidence, not a proof step (rule 3).
* Computations in the literature that the introduction may mention as context (not evidence for the note): HKQ09 s.3.6, s.5.1 (quintic to
  genus 51, assuming the gap); LR22 (genus 53, assuming the gap); CI21 (local `P^2` gap checked to genus 7, per Br25 s.1.2); the GJR
  verification for the quintic to genus 5 (known only through CGYZ26 s.1 and Remark 1.9).

## 6. Literature used: statement numbers, use, how read
"Me": my reading ([A], [H], [P], [abs] as in s.0). "Sources": the project's reading (LIT, referees).

| paper | statements | use in the note | me | sources |
|---|---|---|---|---|
| GJR18 | (17) in Lemma 5.2; Rem. 7.2; Lemma 7.7, (37); Thm 4.1; Rem. 1.6 | (QRR); prior art Gamma, Ginf; intro | [A] exact | LIT [P] |
| GJR17 | (19)-(20) in Lemma 4.6 | (QRR), second source | [A] exact | LIT [P] |
| Giv01 | s.2.2 Thm 1 and Remark | (QRR) basis | [A] | LIT [P] |
| MR1901075 | 4.1, 4.2, 9.1, 9.2(a), s.10, 10.4 | (QRR) vertex derivation; prior art of the chart | [A], [P] | LIT [P] |
| MR2276766 | Thm 1, Rem. (2), (9), Cor. 1, (11), Thm 2 | (QRR) twist sign; Euler integrals | [A], [P] | LIT [P]/ar5iv |
| MR1728879, MR717614 | s.1.1 Mumford; Prop. 2; Thm 4; (22) | (QRR) vertex derivation; evidence | [A] / not read | LIT ar5iv |
| MR2917177 | Thm (2) in (1.15); 8.4; (8.6); 8.7(i)-(ii) (v3) | uniqueness of R; (QRR) | [A] + [P] v2, v3 | summariser only |
| PPZ15 | Defs 2.2, 2.5, 2.13; Props 2.4, 2.9, 2.14; Thm 5 | R-action, translations (Q) | [A], [P] | flat referee [P] |
| Eynard:2007kz | Thm 4.4 (4-23)-(4-24); Thm 5.1 (5-37)-(5-38); Defs 4.1-4.3; Thms 4.1, 4.7 | (Q): TR on the model | [A], [P] | logleg referee (printed) |
| Dunin-Barkowski:2015caa | (8-1), Thm 8.1 | (Q): the model CohFT is TR | [A], [P] | TR_MODEL s.1 |
| ABDKS25 | Thm 7.1, (7.9) | (Q): L-regularity; gap route 2 | [A], [P] | xyswap referee |
| Hock25 | Cor. 1.2, (1.3) | (Q) optional: model gap | [H], [P] | TR_MODEL s.4 |
| CvdLPS23 | Thm 2; Thms 14, 17 | context only (KP; hypermaps) | [A], [P] | "read in part" |
| ZZ08 | (1), Lemma 1, (5), Thms 3-4, (39) | degree bound; normalisation | [A], [P] | LIT [P] |
| MR2553377, GGI16 | Lemma 4.20; Thms 1.1-1.2; Prop. 5.2.4 | context (Gamma structures) | [A] | LIT [P] |
| DLMF | 5.11.1, 5.11.8, 24.4.1, 2.3.x | Stirling; formal 5.11.8 | read | referees |
| CGYZ26, Br25, FCCZ26, LR22, CGL21, CGL26, HK07, HKQ09, Ghoshal:1995wm, Bershadsky:1993cx | as in s.1.1-1.3 | introduction only | [H]/[P]/[abs] | not read by any source |
| LP18, LPfq | (5); Thms 1-2; s.0.4; (26) | definition of the formal quintic; context | [A], [P] | LIT [P]/ar5iv |

## 7. Open issues for the writers (conflicts, gaps, `\draftnote` items, placement)
1. **Teleman numbering (conflict; settled here).** The project cites "arXiv:0712.0160v3 Theorem (2) (s.1.9) and Prop. 8.5(ii) (s.8.2)"
   (MECHANISM_STAGE5 s.2; REFEREES_SYNTHESIS_2 s.5; STAGE5_SYNTHESIS_2 s.7 item 5). These are the numbers of the ar5iv rendering. The v3
   PDF has Theorem (2) in paragraph (1.15), 8.4 Proposition ([D, G1]: `R` unique up to `exp(H_odd)`), paragraph (8.6), and 8.7
   Proposition (ii) (classification by `E exp h`). The conifold note's LABELS already uses v3 numbering. Use v3; `\draftnote{numbers
   from arXiv v3; journal numbering not compared}`. The (QRR) step "the CohFT fixes `h`" should cite 8.4 and 8.7 together.
2. **New prior art (not in any project source).** CGYZ26 (quintic gap proved; Remark 1.9: formal quintic gap open), Br25, FCCZ26,
   LR22, CI21. The introduction must cite CGYZ26 and Br25 (s.2.2). CGYZ26 is a preprint (v2 of 2026-10-31) whose Remark 1.9 may change:
   `\draftnote{re-check arXiv:2609.13921 (latest version) and arXiv:2509.19298 (journal) before submission}`.
3. **Placement.** s01 intro: background (L1-L4), results on the gap (L5-L13), the formal quintic (L14-L19), main theorem and
   corollary, method (L20-L42), what is not claimed (s.2.3). The (QRR) subsection: L19, L20-L24, L27 with exact pins (GJR18 (17) first).
   The topological-recursion section: L28-L31 (state there that TR is applied to `x^4 + v/x` only). The Theorem Gamma section: one
   remark with L35-L39 (Givental s.10, GJR18 Rem. 7.2, MR2553377 Lemma 4.20, MR2276766 (11)). The `Ginf` subsection: L40 (cite GJR18
   Lemma 7.7 as the first proof of the formula for the LG series). The formal stationary-phase appendix: L41-L42.
4. **The corollary "gap of the total" needs its derivation written** (rule 1): the sources state N(g) and verify the total at genus
   2-5, and MECHANISM s.2 names the gap of the total as the weaker property; the step from N(g) to `F_g^c = c_g t_c^(2-2g) + O(1)` needs
   (a) the reconstruction of the formal quintic CohFT from its (QRR)-normalised R-matrix (MR2917177 Thm (2); GJR18 Thm 4.1), (b) the
   conifold-frame map Cf (FORMAL_QUINTIC_genus2 s.1) and (c) Theorem S15 (flat = hatted). If a source does not contain the sentence,
   the writers put a `\draftnote` rather than add an argument. Also compare the note's definition of the conifold expansion with CGYZ26
   Defs 1.1-1.2 and the admissible family of its Remark 1.4.
5. **Corrections to the project's literature notes** (for DRAFT_NOTES.md; the notes themselves stay unchanged): GGI16 Prop. 5.2.4 is
   for `Psi(t)` with exponent `t^(-(N-1)/2)` (LIT_2 s.3.4 has `phi`, `t^((N-1)/2)`); MR2553377 Lemma 4.20 is at `z = -1` (LIT_2 s.3.3:
   `z = 1`); the three illegible items of LIT_3 s.5 in GJR18 s.7 are resolved (L39-L40); LIT_3 s.4 item 1, "Teleman not opened", is
   superseded (L24); the sign `-5 h_alpha` of GJR18 (17), reconstructed in LIT s.1.2, is confirmed in the authors' LaTeX (L27).
6. **Bibliography data.** Root entry Dunin-Barkowski:2015caa outdated; MR717614 wrong entry type; Coates:2012vs lacks its journal
   (s.3.3). To be confirmed: LPfq journal version; CGL26 (Duke 2026, from one reference list); Br25 journal; ABDKS25 volume and pages;
   BGRT journal; Giv01 DOI (from the conifold note).
7. **Keys.** This note uses root keys (rule 9): MR2917177, MR2276766, MR1728879, MR1666787, MR717614 instead of the conifold note's
   Tel12, CG07, FP00, GP99, Mum83. Giv01 and PPZ15 are kept from the conifold note (not in the root file).
8. **Not re-checked by me** (keep the project's wording or add a `\draftnote` if repeated): Hock25's printed example s.3.2 "sign slip"
   (TR_MODEL s.4.1); the `R^-1` convention of arXiv:1712.02528 (LIT s.1.3); the statement numbers of Eti24, Kir10, Pol05, BGRT (LIT_2
   s.2); the link of CvdLPS23 Thm 2 to the `H^2` tower (MECHANISM s.1); BCOV's role (only through HKQ09's attribution).
9. **The constant.** `c_g` = the constant of CGYZ26 Conj. 1.5 in the same normalisation (s.2.1). FORMAL_QUINTIC_genus2 s.4 made no
   statement about HKQ's constant; with CGYZ26 a printed comparison now exists. A remark is optional; `\draftnote` on normalisations.
10. **Do not cite** arXiv:2610.38473 (withdrawn); Lei's lecture notes and the HKUST theses (not read) only by the authors' decision.
11. **Author decisions** (`\draftnote` or DRAFT_NOTES.md): whether to cite the authors' companion notes (resolved conifold, local `P^2`),
    and whether to quote Br25's sentence on R-matrix methods (s.2.3 item 2).
