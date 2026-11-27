# DIGEST_L_5: s.3 the bibliography (part 1): keys, and entries under root keys

Continuation of `DIGEST_L_4.md`. Part 2 (`DIGEST_L_6.md`): complete BibTeX entries for the papers not in the root bibliography.

## 3.1 Rules applied
* Binding rule 9: `NOTE/refs.bib` is the note's own file; every paper already in `refGW/refs.bib` keeps its ROOT key (searched by eprint,
  author and title, s.0.1). The conifold note used short keys for some of these (Tel12, CG07, FP00, GP99, Mum83); this note must not.
* For papers not in the root file, new keys in the conifold note's style (authors' initials and year). Two of them are the conifold
  note's own keys for the same papers (Giv01 = arXiv:math/0008067, PPZ15 = arXiv:1303.1043), reused for consistency.
* Every entry carries `eprint` (the version from which statement numbers are cited goes into `note`, as in the conifold note). Journal
  data: from the arXiv abstract page when listed there ([abs]); otherwise from the reference list of a paper read in full ([bib], named in
  the `note` field) and marked "to be confirmed" when only one secondary source gives them.

## 3.2 Key table

| key | work | arXiv (version cited) | root? | action | pinned numbers used |
|---|---|---|---|---|---|
| Bershadsky:1993cx | BCOV, Kodaira-Spencer theory, CMP 1994 | hep-th/9309140v1 | yes | copy | none (attribution as in HKQ09 s.1) |
| Bershadsky:1993ta | BCOV, holomorphic anomalies, NPB 1993 | hep-th/9302103v1 | yes | copy (optional) | none |
| Ghoshal:1995wm | Ghoshal-Vafa, c = 1 string and the conifold | hep-th/9506122v1 | yes | copy | none |
| HK07 | Huang-Klemm, holomorphic anomaly, gauge theories | hep-th/0605195v2 | no | new | s.2.3 |
| HKQ09 | Huang-Klemm-Quackenbush, compact CY | hep-th/0612125v1 | no | new | (1.2), (2.29), (3.79), s.3.5-3.6, s.5.1 |
| YY04 | Yamaguchi-Yau, polynomials | hep-th/0406078v2 | no | new (optional) | none |
| GJR18 | Guo-Janda-Ruan, structure, quintic | 1812.11908v1 | no | new | abstract; Rem. 1.6; s.3.1; Thm 4.1; Lemma 5.2, (17); Lemma 5.5; Thm 7.1; Rem. 7.2; Cor. 7.3-7.4; Prop. 7.5; Lemma 7.7, (37) |
| GJR17 | Guo-Janda-Ruan, genus two | 1709.07392v1 | no | new | proof of Lemma 4.6, (19)-(20) |
| CGL21 | Chang-Guo-Li, polynomial structure (Ann. Math.) | 1809.11058v1 | no | new | none (finite generation, as summarised in CGYZ26 s.1.1) |
| CGL26 | Chang-Guo-Li, BCOV Feynman rule | 1810.00394v2 | no | new | none |
| CGLL21 | Chang-Guo-Li-Li, NMSP fields | 1809.08806v2 | no | new (optional) | none |
| CGLZ20 | Chang-Guo-Li-Zhou, genus one via MSP | 1711.10118v1 | no | new (optional) | none |
| LR22 | Liu-Ruan, Castelnuovo bound | 2210.13411v2 | no | new | abstract only |
| CGYZ26 | Chang-Guo-You-Zhang, conifold gap, quintic | 2609.13921v2 | no | new | Defs 1.1-1.2, 1.8; Ex. 1.3; Rem. 1.4; Conj. 1.5; Thms 1-2; Rem. 1.9; Prop. 3.1 |
| Br25 | Brini, conifold gap, local P^2 | 2509.19298v3 | no | new | (1.5); Theorem (= Cor. 4.4); s.1.1, footnote 1 |
| CI21 | Coates-Iritani, local P^2 and modular forms | 1804.03292v2 | no | new (optional) | none |
| FCCZ26 | Fang-Chen-Chen-Zong, gap theorem for TR | 2608.11960v1 | no | new | abstract only |
| LP18 | Lho-Pandharipande, stable quotients and HAE | 1702.06096v4 | no | new | s.0.2, (5); s.0.3, (6); Thms 3-4 |
| LPfq | Lho-Pandharipande, HAE for the formal quintic | 1803.01409v2 | no | new | (5); Thms 1-2; s.0.4; (18); (26) |
| ZZ08 | Zagier-Zinger, hypergeometric series | 0710.0889v1 | no | new | (1); Lemma 1; (5); Thms 3-4; (39) |
| Giv01 | Givental, semisimple Frobenius structures | math/0008067v4 | no | new (conifold note's key) | s.2.2, Thm 1 and the Remark after it; s.2.3 |
| MR1901075 | Givental, quantization of quadratic hamiltonians | math/0108100v2 | yes | copy + eprint | 4.1 Thm; 4.2 Rem.; 9.1 Thm; 9.2 Rem. (a); s.10 (10.1, 10.4) |
| MR1408320 | Givental, equivariant GW invariants | alg-geom/9603021v2 | yes | copy (optional) | none |
| MR2276766 | Coates-Givental, quantum Riemann-Roch | math/0110142v2 | yes | copy + eprint | Thm 1, (7); Rem. (2); (9); Cor. 1; s.7, (11); Thm 2; App. 1 |
| MR1728879 | Faber-Pandharipande, Hodge integrals | math/9810173v1 | yes | copy + eprint | s.1.1 Theorem (Mumford); Props 1-2; Thm 4; s.2.1, (22) |
| MR717614 | Mumford 1983 | (none) | yes | copy, fix entry type | Mumford's relation and formula |
| MR1666787 | Graber-Pandharipande, localization | alg-geom/9708001v2 | yes | copy + eprint | the localisation formula |
| MR2917177 | Teleman, 2D semisimple field theories | 0712.0160v3 | yes | copy + note | Thm (2) in (1.15); 8.4 Prop.; (8.6); 8.7 Prop. (i)-(ii) (v3 numbering) |
| PPZ15 | Pandharipande-Pixton-Zvonkine, 3-spin | 1303.1043v3 | no | new (conifold note's key) | Thm 5; Defs 2.2, 2.5, 2.13; Props 2.3, 2.4, 2.9, 2.14 |
| Eynard:2007kz | Eynard-Orantin, invariants of algebraic curves | math-ph/0702045v4 | yes | copy + note | Defs 4.1-4.3, (4-3), (4-9), (4-10); Thms 4.1, 4.4, 4.7, 5.1; (4-23)-(4-24); (5-37)-(5-38) |
| Dunin-Barkowski:2015caa | DNOPS, superpotential as global spectral curve | 1509.06954v2 | yes | key kept, data corrected | (8-1); Thm 8.1 |
| ABDKS25 | Alexandrov-Bychkov-Dunin-Barkowski-Kazarian-Shadrin, x-y swap | 2212.00320v3 | no | new | Thm 7.1, (7.9); Thm 1.14 |
| Hock25 | Hock, symplectic (non-)invariance of F_g | 2502.18115v2 | no | new | Thm 1.1; Cor. 1.2, (1.3) |
| CvdLPS23 | Carlet-van de Leur-Posthuma-Shadrin, hypermaps, KP | 2211.12259v1 | no | new (optional) | Thm 2; Thms 14, 17 |
| DuninBarkowski:2012bw | DOSS, Givental formula = TR | 1211.4021v1 | yes | copy (optional) | none |
| dunin2016primary | DNOPS, primary invariants | 1605.07644v2 | yes | copy (optional) | none |
| MR2553377 | Iritani, integral structure | 0903.1463v3 | yes | copy + eprint | Thms 1.1-1.2 (= 4.11, 4.14); Lemma 4.20 |
| GGI16 | Galkin-Golyshev-Iritani, Gamma conjectures | 1404.6407v4 | no | new | Thm 1.3.1 (= 5.0.1); (5.2.3); Prop. 5.2.4 |
| Coates:2012vs | Coates-Iritani, convergence of Givental's formula | 1203.4193v2 | yes | key kept, journal added (optional) | none |
| DLMF | NIST DLMF, Version 1.2.8 | (none) | no | new | 5.11.1, 5.11.8, 24.4.1, 2.3.8, 2.3.13-2.3.18, s.2.3(iv) |
| Eti24, Kir10, Pol05, BGRT02 | formal Gaussian integrals (appendix context) | see `DIGEST_L_6.md` | no | new (optional) | as LIT_2 s.2 (not re-checked) |

Not to be cited: arXiv:2610.38473 (withdrawn). Lecture notes and theses cited by CGYZ26 (Lei 2025, You 2025, Zhang 2025): only if the
authors decide (`\draftnote`); not read.

## 3.3 Entries under root keys: copy from `refGW/refs.bib`, with these changes (checked against the arXiv pages)
* **Dunin-Barkowski:2015caa** (root: "J. Inst. Math. Jussieu, pages 1--49, 2017", the online-first data). Use:
  ```
  @article{Dunin-Barkowski:2015caa,
    author = {Dunin-Barkowski, P. and Norbury, P. and Orantin, N. and Popolitov, A. and Shadrin, S.},
    title = {Dubrovin's superpotential as a global spectral curve},
    journal = {J. Inst. Math. Jussieu}, volume = {18}, number = {3}, year = {2019}, pages = {449--497},
    doi = {10.1017/S147474801700007X}, eprint = {1509.06954}, archivePrefix = {arXiv}, primaryClass = {math-ph},
    note = {arXiv:1509.06954v2; numbers cited from this version}}
  ```
* **MR2917177** (Teleman): keep the root fields; add `note = {arXiv:0712.0160v3; numbers cited from the v3 PDF (the ar5iv rendering
  renumbers sections 1 and 8); journal numbering not compared}`; the DOI on the arXiv page is 10.1007/s00222-011-0352-5.
* **MR717614** (Mumford): the root entry is typed `@article` with `journal = {Progr. Math.}`; use `@incollection` with `booktitle =
  {Arithmetic and geometry, Vol. II}`, `series = {Progr. Math.}`, `volume = {36}`, `pages = {271--328}`, `publisher = {Birkh\"{a}user
  Boston, Boston, MA}`, `year = {1983}` (the conifold note's Mum83 has exactly these fields).
* **Coates:2012vs**: add `journal = {Michigan Math. J.}`, `volume = {64}`, `number = {3}`, `year = {2015}`, `pages = {587--631}`,
  `doi = {10.1307/mmj/1441116660}` ([abs]).
* **Add `eprint`, `archivePrefix` and a version note** to MR1901075 (math/0108100v2), MR2276766 (math/0110142v2), MR1728879
  (math/9810173v1), MR1666787 (alg-geom/9708001v2), MR2553377 (0903.1463v3), MR1408320 (alg-geom/9603021v2), and `note =
  {arXiv:math-ph/0702045v4; numbers cited from this version}` to Eynard:2007kz. Their journal data in the root file agree with the
  arXiv pages where those list them (MR2553377: Adv. Math. 222 (2009), no. 3, 1016-1079; the others list none).
* Bershadsky:1993cx, Bershadsky:1993ta, Ghoshal:1995wm, DuninBarkowski:2012bw, dunin2016primary: copy unchanged (arXiv pages agree:
  CMP 165 (1994) 311-428; NPB 405 (1993) 279-304 (the arXiv page lists S. Katz as a fifth author; the root entry has four); NPB 453 (1995) 121-128; CMP 328 (2014)
  669-700; dunin2016primary may gain `booktitle = {Topological recursion and its influence in analysis, geometry, and topology}`,
  `series = {Proc. Sympos. Pure Math.}`, `volume = {100}`, `pages = {297--331}`, `publisher = {Amer. Math. Soc.}`, `year = {2018}`,
  `doi = {10.1090/pspum/100/01768}` [abs], if it is cited at all).
