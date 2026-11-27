# Bibliography plan (part 1 of 2: rules, key table, statement numbers, changes to root entries)

`NOTE/refs.bib` is built from this file and `BIBPLAN_2.md` (complete new entries). Source of the data: DL s.3 (DIGEST_L_5..L_7, every
new entry checked by DL against the arXiv abstract page), the statement numbers read in the full text in the gap phase (G1 s.1.2, s.7;
G1R_2 s.10; G4R s.1; G5 s.8), corrected below where the designer or the editor re-checked them; revised after the gap phase (record:
`design/DESIGN_CHANGES.md`).

## Rules
1. Rule 9: every paper present in `refGW/refs.bib` keeps its ROOT key, copied verbatim with the changes of the last table (never the
   conifold note's short keys Tel12, CG07, FP00, GP99, Mum83). Papers not in the root file get new keys (authors' initials + year, the
   conifold note's style); Giv01 and PPZ15 are the conifold note's own keys for papers absent from the root file and are reused.
2. Every arXiv item carries `eprint`, `archivePrefix`, and in `note` the version from which statement numbers are cited; provenance
   sentences ("journal data from ..", "to be confirmed") go into `annote`, which amsalpha does not print.
3. Statement numbers are cited exactly as below (arXiv version named); the journal numbering is not compared: `\draftnote` at the first
   citation of [MR2917177] (v3 PDF numbering differs from the ar5iv rendering) and of [MR1408320] (arXiv v2 numbering; IMRN not compared).
4. Not cited: arXiv:2610.38473 (withdrawn); Hoermander, ALPDO I (cited from memory in the sources); Lei's lecture notes and the theses
   cited by [CGYZ26] (not read; only by the authors' decision); arXiv:hep-th/9204035 and classical special geometry (not read); the
   Chester-Friedman-Ursell paper (named in the sources without a reference); GJR18's [22] Hermite and [30] Nemes (DLMF 5.11.8 suffices);
   Madsen-Weiss (wrongly attributed in the derivation of (QRR), G1R_2 s.4.5); Witten 1991 (the genus-1 logarithmic term is not printed);
   [MR1653024] (not needed with the pin on [MR1408320]); the appendix-A context items Eti24, Kir10, Pol05, BGRT (dropped under D8).

## Key table (cited works; "root" = key of `refGW/refs.bib`)

| key | work (arXiv version) | root? | statements cited | where |
|---|---|---|---|---|
| Eynard:2007kz | Eynard-Orantin, Invariants of algebraic curves and topological expansion (math-ph/0702045v4) | yes | Def. 4.1 with the kernel (4-3), Def. 4.2 (4-9), Def. 4.3 (4-10); Thm 4.1 (4-20) (sign comparison only); Thm 4.2; Thm 4.6; Thm 4.7 (4-30); Thm 5.1 (5-14)-(5-15) in s.5.1; s.5.3 (5-27)-(5-30); (5-1), (5-4), (5-13) with footnote 8; optional pins Thm 4.4 (4-23)-(4-24), Thm 4.5 (4-25), Cor. 4.1 (4-26)-(4-29), App. A (A-1), (A-5), (A-25) | s06, App. B |
| DuninBarkowski:2012bw | DOSS, Identification of the Givental formula with the spectral curve topological recursion procedure (1211.4021v1, the only version) | yes | s.2.3 (leaf, dilaton-leaf, edge and vertex weights; the graph formula); Thm 4.1 | s06 ((Q)), App. B |
| Dunin-Barkowski:2015caa | DNOPS, Dubrovin's superpotential as a global spectral curve (1509.06954v2) | yes (data corrected) | Remark 2.1; s.2.3 (the action; its dilaton-leaf sign is misprinted); s.2.5 (2-26)-(2-29); s.2.6 (2-30)-(2-35); Prop. 4.3 (1); s.8 (8-1) at `a = (0,0,0,0,4v)`, the sentence before Thm 8.1; Thm 8.1 only as context | s06, App. B |
| ABDKS25 | Alexandrov-Bychkov-Dunin-Barkowski-Kazarian-Shadrin, x-y swap (2212.00320v3) | no | Thm 7.1 with (7.6), (7.7) (factor `1/hbar^2`), (7.9), (7.10); kernel (1.1); example (7.2); Remark 7.3; Thm 1.14 (context) | s06, App. B |
| Hock25 | Hock, Symplectic (non-)invariance of the free energy (2502.18115v2) | no | Cor. 1.2, formula (1.3); kernel (2.1); `F_g` (1.1); only in `rmk:mod-hock` and `conv:mod-kernel` (not an input, DEC Q4) | s06 |
| PPZ15 | Pandharipande-Pixton-Zvonkine, Relations via 3-spin structures (1303.1043v3) | no | Def. 0.3 (iii); Rem. 2.1 and (6); Def. 2.2 with the convention sentence after it (p.20); Prop. 2.3; Prop. 2.4; Def. 2.5; Prop. 2.9; Prop. 2.12; Def. 2.13; Prop. 2.14; Thm 5 (context) | s02, s04, s06, s07, s10, App. B |
| GJR18 | Guo-Janda-Ruan, Structure of higher genus GW invariants of quintic 3-folds (1812.11908v1) | no | abstract; Rem. 1.6; s.3.1 eq. (9); s.4.1 (the twisted CohFT), Thm 4.1; s.5.2 eq. (14); Lemma 5.2 with eq. (17); Rem. 7.2; Lemma 7.7 eq. (37); NOT eq. (13) (misprinted) | s01, s02, s09 |
| GJR17 | Guo-Janda-Ruan, A mirror theorem for genus two GW invariants of quintic threefolds (1709.07392v1) | no | proof of Lemma 4.6, eqs. (19)-(20) | s02 |
| Giv01 | Givental, Semisimple Frobenius structures at higher genus (math/0008067v4) | no | s.2.2 Thm 1 and the Remark after it | s02, App. B |
| MR1408320 | Givental, Equivariant Gromov-Witten invariants (alg-geom/9603021v2) | yes | s.9 set-up (p.23-27); Prop. 11.4; Thm 11.6 (the pin of (M)); Cor. 11.7 only for its misprint | s02, App. B |
| MR2276766 | Coates-Givental, Quantum Riemann-Roch, Lefschetz and Serre (math/0110142v2) | yes | Thm 1 (7); Rem. (2); s.4 (9), Cor. 1; s.6 (10) (normalisation of `J`); s.7 (11); Thm 2 only as the modern form of (M) (non-equivariant: not the pin); no Cor. 5.1 exists | s02, s09, App. B |
| MR2917177 | Teleman, The structure of 2D semi-simple field theories (0712.0160v3) | yes | Theorem (2) in paragraph (1.15); 8.2 Proposition; 8.4 Proposition; paragraph (8.6); 8.7 Proposition (i)-(ii); the Riemann-Roch identity of p.40; 8.10 Proposition (context) (v3 PDF numbering) | s02, s10 |
| MR1901075 | Givental, Gromov-Witten invariants and quantization of quadratic hamiltonians (math/0108100v2) | yes | 4.1 Theorem; 4.2 Remark; 9.1 Theorem; 9.2 Remarks (a); s.10: 10.1 Theorem, 10.4 Lemma | s01, s09, App. B |
| MR1728879 | Faber-Pandharipande, Hodge integrals and GW theory (math/9810173v1) | yes | s.2.1 Mumford's relation, (22): only `rmk:set-qrrderiv` (with its `\draftnote`) | s02 |
| MR717614 | Mumford, Towards an enumerative geometry of the moduli space of curves (1983) | yes (type fixed) | Mumford's relation, only if `rmk:set-qrrderiv` cites it directly | s02 (remark) |
| MR1666787 | Graber-Pandharipande, Localization of virtual classes (alg-geom/9708001v2) | yes | the virtual localisation formula (definition of the invariants) | s02, App. C |
| AC96 | Arbarello-Cornalba, Combinatorial and algebro-geometric cohomology classes on the moduli spaces of curves (alg-geom/9406008v1) | no | formula (1.7); `kappa_0 = 2g - 2 + n` (p.5) | s02 (`lem:set-dilaton`) |
| LP18 | Lho-Pandharipande, Stable quotients and the holomorphic anomaly equation (1702.06096v4) | no | s.0.2 eq. (5); s.0.3 eq. (6) | s01, s02 |
| LP20 | Lho-Pandharipande, Holomorphic anomaly equations for the formal quintic (1803.01409v2; key `LPfq` in the digests and gap notes, "LP18" in G1) | no | s.0.2 eq. (5); (10); Thms 1-2; s.0.4; s.2.4.3 (17), (19)-(20) (corroboration only); s.3.2 eq. (26) | s01, s02, App. C |
| ZZ08 | Zagier-Zinger, Some properties of hypergeometric series associated with mirror symmetry (0710.0889v1) | no | eq. (1); Lemma 1; expansion (5); Thm 4 (i)-(ii) | s01, s02, App. B |
| GR19 | Guo-Ross, The genus-one global mirror theorem for the quintic threefold (1703.06955v1) | no | Props. 6.1-6.2 (corroboration of `prop:set-cyclicLR`; optional sentence) | s02 (optional) |
| LLY97 | Lian-Liu-Yau, Mirror principle I (alg-geom/9712011v1) | no | only if the writer of s02 reads it and pins a statement (DEC Q7) | s02 (conditional) |
| MR2553377 | Iritani, An integral structure in quantum cohomology (0903.1463v3) | yes (add eprint) | Lemma 4.20 (at `z = -1`) | s09 (`rmk:gam-prior`) |
| CGYZ26 | Chang-Guo-You-Zhang, Conifold gap for the quintic threefold (2609.13921v2) | no | Defs. 1.1-1.2; Ex. 1.3; Rem. 1.4; Conj. 1.5; Def. 1.8 (sign warning only); Thm 1; Thm 2 (the quintic only, REV R5); Rem. 1.9; Prop. 3.1; s.4.2 Prop.-Def. 4.5 (1) | s01, s10, App. B |
| Br25 | Brini, Conifold gap and all-genus mirror symmetry for local `P^2` (2509.19298v3) | no | (1.5); the main Theorem (= Cor. 4.4); s.1.1 footnote 1 | s01, App. B |
| FCCZ26 | Fang-Chen-Chen-Zong, Conifold gap theorem for topological recursion (2608.11960v1) | no | abstract only (no statement used) | s01 |
| HKQ09 | Huang-Klemm-Quackenbush, Topological string theory on compact CY (hep-th/0612125v1) | no | (1.2); s.3.5 (3.79); s.5.1 | s01, App. B |
| HK07 | Huang-Klemm, Holomorphic anomaly in gauge theories and matrix models (hep-th/0605195v2) | no | s.2.3 | s01 |
| Bershadsky:1993cx | BCOV, Kodaira-Spencer theory of gravity (hep-th/9309140) | yes | none (attribution through HKQ09 s.1) | s01 |
| Ghoshal:1995wm | Ghoshal-Vafa, c = 1 string as the topological theory of the conifold (hep-th/9506122) | yes | none (abstract) | s01 |
| CGL21 | Chang-Guo-Li, Polynomial structure of GW potential of quintic 3-folds (1809.11058v1) | no | none (as summarised in CGYZ26 s.1.1) | s01 |
| CGL26 | Chang-Guo-Li, BCOV Feynman rule of the quintic threefold (1810.00394v2) | no | none | s01 |
| LR22 | Liu-Ruan, Castelnuovo bound (2210.13411v2) | no | abstract only | s01 |
| CI21 | Coates-Iritani, GW invariants of local `P^2` and modular forms (1804.03292v2) | no | none (as reported in Br25 s.1.2) | s01 (optional) |
| CvdLPS23 | Carlet-van de Leur-Posthuma-Shadrin, hypermaps and Hirota equations (2211.12259v1) | no | the name `H_(0;1,4)` only (`rmk:mod-H2` is dropped) | App. B |
| GGI16 | Galkin-Golyshev-Iritani, Gamma conjectures (1404.6407v4) | no | Thm 1.3.1 (= 5.0.1); Prop. 5.2.4 | s01 (optional context) |
| DLMF | NIST Digital Library of Mathematical Functions, Version 1.2.8 | no | 24.4.1 (from 24.2.3) in `lem:fsp-dlmf5118`; 5.11.1, 5.11.8, 2.3.8 only in `rmk:fsp-classical` | App. A, s09 |
| Optional, s01 | Bershadsky:1993ta, YY04, CGLL21, CGLZ20, dunin2016primary, Lei24 | mixed | none | s01 if the authors want them |

## Changes to root entries when copying into `NOTE/refs.bib` (DL s.3.3; gap phase)
* **Dunin-Barkowski:2015caa**: replace the online-first data by J. Inst. Math. Jussieu 18 (2019), no. 3, 449-497, DOI
  10.1017/S147474801700007X; add `note = {arXiv:1509.06954v2; numbers cited from this version}` (full entry in `BIBPLAN_2.md`).
* **DuninBarkowski:2012bw**: add `note = {arXiv:1211.4021v1 (the only version); numbers cited from this version}`.
* **MR1408320**: add `eprint = {alg-geom/9603021}`, `archivePrefix = {arXiv}`, `note = {arXiv:alg-geom/9603021v2; numbers cited from this
  version; IMRN numbering not compared}`.
* **MR2917177**: add `note = {arXiv:0712.0160v3; numbers cited from the v3 PDF (the ar5iv rendering renumbers sections 1 and 8);
  journal numbering not compared}`; DOI 10.1007/s00222-011-0352-5 (arXiv page).
* **MR717614**: type `@incollection`, `booktitle = {Arithmetic and geometry, Vol. II}`, `series = {Progr. Math.}`, `volume = {36}`,
  `pages = {271--328}`, `publisher = {Birkh\"auser Boston, Boston, MA}`, `year = {1983}`.
* **Eynard:2007kz**: add `note = {arXiv:math-ph/0702045v4; numbers cited from this version}` (numbers re-checked in the v4 PDF, OUTLINE_2
  s.0.4 item 1).
* **MR1901075, MR2276766, MR1728879, MR1666787, MR2553377**: add `eprint` (math/0108100, math/0110142, math/9810173, alg-geom/9708001,
  0903.1463), `archivePrefix = {arXiv}` and the version note (v2, v2, v1, v2, v3); journal data of the root file agree with the arXiv pages.
* **Bershadsky:1993cx, Ghoshal:1995wm** (and Bershadsky:1993ta, dunin2016primary if cited): copy unchanged (Bershadsky:1993ta: the arXiv
  page lists S. Katz as a fifth author, the root entry four: `\draftnote` if cited). **Coates:2012vs**: not cited.

## Points to confirm before submission (`DRAFT_NOTES.md`)
LP20 journal version (unconfirmed lead: Peking Math. J., DOI 10.1007/s42543-018-0008-0); AC96 journal data (J. Algebraic Geom. 5 (1996)
705-749, from a reader's memory; the arXiv page lists none); LLY97 journal data, only if cited; CGL26 journal data (one reference list
only); Br25 journal (to appear); ABDKS25 volume and pages; Giv01 DOI (from the conifold note); the IMRN numbering of [MR1408320] and the
journal numbering of [MR2917177]; the latest versions of CGYZ26 (v2) and Br25 (v3), with CGYZ26 Defs. 1.1-1.2, Ex. 1.3, Rem. 1.4, Conj. 1.5
and Prop.-Def. 4.5 (1) re-read in the PDF; the published numbering of every journal item (statement numbers are from the arXiv versions).
Cosmetic, from the root data copied unchanged: MR1901075's `Mrnumber` field already contains "MR", so amsalpha prints "MR MR1901075";
`dunin2016primary` (optional) has no journal (BibTeX warns only if it is cited); `ABDKS25` has no pages yet (BibTeX warning).
