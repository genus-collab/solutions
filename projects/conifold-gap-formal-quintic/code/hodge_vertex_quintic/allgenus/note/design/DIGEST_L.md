# DIGEST_L: literature positioning for the introduction, and the bibliography (key read_L)

Design digest for the note in `allgenus/note/` (no LaTeX here). Pieces: this file (s.0 sources), `DIGEST_L_2.md` (s.1, part 1:
the conifold gap and the formal quintic), `DIGEST_L_3.md` (s.1, part 2: Givental-Teleman and (QRR), topological recursion and the
model curve, Mellin-Barnes and Gamma; what is new), `DIGEST_L_4.md` (s.2 the novelty statement; s.4 conventions of the literature),
`DIGEST_L_5.md`, `DIGEST_L_6.md` and `DIGEST_L_7.md` s.3.5 (s.3 the bibliography: keys, corrections of root entries, complete new
entries), `DIGEST_L_7.md` (also s.5 evidence, s.6 table of literature used, s.7 open issues). The sections s.0-s.3 follow the task's layout for this key; s.4-s.7 are those of the generic digest format. The generic s.1-s.3
(definitions, statements, proofs) do not apply to a literature digest: here s.1 is "what is known, claim by claim" (claims L1-L48),
s.2 the novelty statement, s.3 the bibliography.

Status labels. For project statements: as in the sources (PROVED, quoted, VERIFIED, CONJECTURAL, SKETCH). For a literature item:
"proved there", "stated without proof", "conjecture/prediction (physics)", "assumed", and how I read it:
* [A] full text in the ar5iv rendering, converted by my own stdlib script that replaces every formula by the authors' LaTeX source
  (the `alttext` of the MathML), so signs and Greek letters are exact. ar5iv renders the latest arXiv version. CAUTION: for papers
  whose theorem counter is shared with numbered paragraphs (Teleman), ar5iv's statement numbers differ from the PDF (s.7 item 1);
* [H] full text in arXiv's own HTML (same conversion), for recent papers;
* [P] full text in my own text extraction of the arXiv PDF (minus signs and Greek letters can be lost; used for numbering and for
  papers without a usable HTML rendering);
* [abs] the arXiv abstract page only (title, authors, versions, journal reference, DOI, comments);
* [bib] bibliographic data taken from the reference list of another paper read in full (named).

## 0. Sources read

### 0.1 Project files (read-only)
* `allgenus/MECHANISM_STAGE5.md`: all (s.1-s.5). The assembled statement s.2: (U), (Q), (QRR); s.3.6 literature; s.5 item 2 plan.
* `allgenus/MECHANISM_STAGE4.md`: all. `allgenus/MECHANISM_STAGE2.md`: all. `allgenus/MECHANISM.md`: all (s.1 cites arXiv:2211.12259
  Thm 2; s.2 "no argument about total potentials ... can prove Conjecture N").
* `allgenus/SYNTHESIS_STAGE01.md`: s.4 table end and s.5 only (lines 155-179).
* `allgenus/stage5/lit/LIT.md`, `LIT_2.md`, `LIT_3.md`: all (the project's own literature notes; every claim the introduction could
  repeat is re-checked in s.1 below, marked "LIT confirmed" or "LIT corrected"). `stage5/lit/pdftext.py` read; a copy was used.
* `NATURAL_FORM.md` (package root): s.0-s.3 and s.5 (lines 1-165): literature remarks are s.0 (Hurwitz Frobenius manifold `M_{0;3,0}`,
  `d = 1/2`), s.2 (DNOPS-type superpotential with a simple pole; arXiv:2211.12259 "read only in part, not used").
* `FORMAL_QUINTIC_genus2.md` s.1-s.5 (set-up citing arXiv:1702.06096 s.0.2; conifold frame "HKQ-style"; the genus-2 total gap and
  the remark "no statement about the HKQ constant is made", s.4).
* `allgenus/stage2/STAGE2_SYNTHESIS.md` s.0-s.2 (lines 1-115: the quoted theorems (Q) and where each is used, table M1-M10);
  `stage2/STAGE2_SYNTHESIS_2.md` s.5 (O7, O13) and s.6 (write-up plan s.6.3, human checks s.6.2).
* `allgenus/stage5/STAGE5_SYNTHESIS_2.md` s.5 (K6-K14, in particular K7, K9), s.6 (O11, O15, O17), s.7-s.8 (the literature tasks:
  Teleman s.8.2 verbatim; GJR (17) and s.7 typeset; cite GJR Remark 7.2, Lemma 7.7 and Givental s.10).
* `allgenus/stage5/referees/REFEREES_SYNTHESIS*.md`: searched for literature items only (DLMF 5.11.8, 24.4.1; the (Q), (QRR) lists).
* The conifold note (template, read-only): `OUTLINE.md` s.0 (binding rules) and s01_intro plan; `LABELS.md` "Bibliography keys";
  `refs.bib` (all); `main.tex` lines 1-120 (preamble, `\draftnote`). Not read: its `NOTATION.md`, `DRAFT_NOTES.md`, section files.
* Root bibliography `refGW/refs.bib` (read-only): searched by eprint number, author and title for every paper below (s.3).

### 0.2 Literature read on the web (arXiv and DLMF only; generic user agent; no personal data in any request)
Versions are those served on 2026-11-14/15.
* Conifold gap, quintic and formal quintic:
  * hep-th/0612125v1 HKQ [P] s.1, s.2.2, s.3.5-s.3.6, s.4, references; [abs]. hep-th/0605195v2 Huang-Klemm [P] s.2.3; [abs].
  * hep-th/9506122 Ghoshal-Vafa [abs]. hep-th/9309140, hep-th/9302103 BCOV [abs] only.
  * 1702.06096v4 Lho-Pandharipande [A] s.0.2-s.0.6 headings; [P] searched; [abs].
  * 1803.01409v2 Lho-Pandharipande [A] s.0.2-s.0.4, s.3.2 (26); [P] searched ("conifold", "gap": no occurrence); [abs].
  * 1812.11908v1 GJR [A] abstract, s.1 (Remarks 1.6-1.7), s.3.1, s.4.1, Lemma 5.2 with (17), s.7 (Thm 7.1 to Lemma 7.9), references;
    [P] cross-checked numbering; [abs]. 1709.07392v1 GJR [A] Lemma 4.6 with (17)-(22); [abs].
  * 1809.11058v1, 1810.00394v2, 1809.08806v2, 1711.10118v1 (Chang-Guo-Li(-Li), Chang-Guo-Li-Zhou) [abs] only.
  * 2210.13411v2 Liu-Ruan [abs]. 2409.11659v2 Lei [abs].
  * **2609.13921v2 Chang-Guo-You-Zhang [H]: s.1 in full (Defs 1.1-1.2, 1.7-1.8, Remarks 1.4, 1.6, 1.9, Conj. 1.5, Thms 1-3), s.3.1
    (Prop. 3.1, Lemma 3.2), s.4 opening, Thm 4, Prop.-Def. 4.5, Thm 5, references; [P] and [abs].** Not read: s.2 beyond Thm 3, s.5,
    appendix, the proofs in s.3.2 and s.4.3.
  * 2509.19298v3 Brini [H] s.1.1-s.1.2 (main theorem, Conjecture (1.5)), references; [abs]. 2608.11960v1 Fang-Chen-Chen-Zong [abs].
  * 1804.03292v2 Coates-Iritani [abs]. 2610.38473 Hecht [abs] (v2 withdrawn on 2026-10-30; not to be cited).
* Givental-Teleman and (QRR):
  * math/0008067v4 Givental [A] s.2.2 (Thm 1, Remark), headings; [P] numbering. math/0108100v2 Givental [A] s.4, s.9, s.10; [P].
  * math/0110142v2 Coates-Givental [A] s.3 (Thm 1, Remarks), s.4 (Cor. 1), s.7 (11), Thm 2; [P] numbering.
  * math/9810173v1 Faber-Pandharipande [A] s.0 (Thms 1, 4), s.1.1 (Mumford), Prop. 2, s.2.1 (22); [P] numbering.
  * **0712.0160 Teleman [A] (ar5iv) s.1.9 and s.8.1-s.8.2, and [P] of v2 and v3 for the numbering of s.1 and s.8 (they differ from
    ar5iv: s.7 item 1).** 1303.1043v3 PPZ [A] s.1 (Thm 5), s.2 (Defs 2.2, 2.5, 2.13; Props 2.3, 2.4, 2.7-2.9, 2.14); [P] numbering.
  * alg-geom/9708001v2 Graber-Pandharipande, alg-geom/9603021v2 Givental [abs] only.
* Topological recursion and the model:
  * math-ph/0702045v4 Eynard-Orantin [A] s.2.2 (list of theorems), s.4.1-s.4.4 (Defs 4.1-4.3, Thms 4.1-4.7, Cor. 4.1), s.5.1 (Thm 5.1);
    [P] numbering. 1509.06954v2 DNOPS [A] s.6 (Thm 6.1), s.8 (Thm 8.1, Remark 8.2, proof start); [P] numbering.
  * 2212.00320v3 ABDKS [A] s.7 (Thm 7.1, Remarks 7.2-7.4); [P] numbering. 2502.18115v2 Hock [H] s.1 (Thm 1.1, Cor. 1.2); [P].
  * 2211.12259v1 Carlet-van de Leur-Posthuma-Shadrin [A] introduction, Thm 2 with proof start; [P] numbering of Thms 14, 17.
  * 1211.4021 DOSS, 1605.07644 DNOPS, 1110.2949 Eynard [abs] only.
* Mellin-Barnes, stationary phase, Gamma:
  * 0903.1463v3 Iritani [A] s.1 (Thms 1.1-1.2), Lemma 4.20 with proof; headings of s.4. 1404.6407v4 GGI [A] Thm 1.3.1, s.5.2.
  * 0710.0889v1 Zagier-Zinger [A] s.1 (eq. (1), Lemma 1, (5), Thms 3-4), s.3 Conjecture (39); [P] numbering.
  * DLMF (Version 1.2.8, release 2026-09-15): s.5.11 (5.11.1, 5.11.8 with its errata and generalisation), 24.4.1, s.2.3 (2.3.8,
    2.3.13-2.3.18) [read in full on dlmf.nist.gov].
  * 2409.03117v2 Etingof, 0810.1700v2 Kirwin, math/0406251 Polyak, math/9801049v4 Bar-Natan-Garoufalidis-Rozansky-Thurston [abs] only
    (their statement numbers are LIT_2 s.2's, read there through a summariser; not re-checked).

### 0.3 Not read, and limits
* Not on an allowed site, so not read: Lei's lecture notes (GitHub), the HKUST theses of L. You and H. Zhang, GJR's "log GLSM
  localization formula under development" (unpublished), Hermite (1895), Nemes (2013), Jahnke-Emde, Mumford (1983), Olver; the journal
  versions of all papers (no journal page was opened; journal data come from arXiv pages or from [bib], each marked in s.3).
* arXiv:1812.11908v1 is the only version of GJR; arXiv:2609.13921 may be revised again (v2 is of 2026-10-31): s.7 item 2.
* I ran no mathematics code. Scripts (stdlib, `py -I -B`, in the session scratch folder): a copy of `pdftext.py`, `ar5iv2txt.py`
  (HTML to text with LaTeX `alttext`), `absinfo.py` (abstract-page fields). Downloads were kept in separate scratch folders.
* Nothing in the project was modified; the only files written are `note/design/DIGEST_L*.md`.
