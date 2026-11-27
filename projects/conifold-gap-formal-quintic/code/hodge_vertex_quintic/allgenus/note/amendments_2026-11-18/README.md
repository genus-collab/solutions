# Suggested amendments after the audit of 2026-11-17/18

Review of `../audit_2026-11-18/` (README, 44 reports, 8 judgements). Every amendment below was checked against the current
TeX (the quoted old text was matched in the sources) and is written as an exact change, old text -> new text. Nothing in
the note has been edited and nothing was compiled. All of this is AI work; no source paper was re-read for this review,
so corrected statement numbers rest on the audit referees' readings.

## Files

| file | content | items |
|---|---|---|
| `AMEND_1_chain.md`, `_2`, `_3`, `_4` | the audit's three must-fix points (A1-A3), five routine missing steps (B1-B5), statements false as written (C1-C6), wording from E8 (D1-D9), and the inputs as they read after the amendments | 23 items, about 60 edits |
| `AMEND_2_units_01_08.md`, `_2` | minor findings: setup, cyclic systems, Theorem U, canonical gauge, colours, transport, Theorem R | P1 11, P2 9, P3 19, rejected 7 |
| `AMEND_3_units_09_18.md`, `_2`, `_3` | minor findings: the model, comparison theorem, Theorem P, the glue, Theorem Gamma, formal stationary phase, main theorem, Complements | P1 12, P2 13, P3 38, rejected 13 |
| `AMEND_4_cross_lit.md`, `_2`, `_3` | cross-cutting reports (circularity, interfaces, signs, rings, ledger of inputs) and literature reports; table of 18 pins to correct | P1 10, P2 14, P3 8, rejected 8 |

P1 = a statement or proof step is wrong or incomplete as written; P2 = citation, pin or cross-reference; P3 = wording
or notation, optional.

## Recommended order

1. **AMEND_1, A1-A3** (must fix), then **B1-B5**. These change what the note declares as inputs and close the five
   routine steps. Apply them first: several later items touch the same lines.
2. **All P1 items** of AMEND_2, AMEND_3, AMEND_4 (33 before removing duplicates, about 28 after). Each is a one-phrase
   or one-paragraph repair; none changes a main result.
3. **P2 items and the pin table** (AMEND_4, rows 1-18), after a human has opened the sources for the single-referee pins
   listed under "Not verified" in each file.
4. **P3 items** at the authors' discretion.

## Duplicates between the files (apply once)

| topic | keep | drop |
|---|---|---|
| `rmk:set-qrrsign` (ph_1 not forced without (QRR)) | AMEND_1 C1 | A2-01, A4-03 |
| `cor:gau-cycsys` (ii), closing relation in `prop:gau-quintic` | AMEND_1 C2 | A2-02, A2-03 |
| `thm:mod-var` (ii) and EO (5-15) | A4-02 (re-proof from `prop:mod-period`) if (ii) is kept; else AMEND_1 C3 | A3-13 |
| remarks resting on unrecorded computations | AMEND_1 C4-C6 | A3-11, A4-09 |
| ledger sentences (`s02_setup` l.14, `s01_intro_3` l.17) | A4-05, after AMEND_1 A1 | A2-11 |
| PPZ lists after Corollary B | AMEND_1 B4 | A3-24, A3-59 apply with it |
| `conv:mod-kernel` and DNOPS (2-26)-(2-29) | A4-16 | A3-15 |

## Points where the amendments depart from the audit's own repairs

1. **Flatness of R (E4).** The judge redefines R by "flat, with the (QRR) value at q = 0". AMEND_1 A1 keeps the
   definition, declares flatness as part of (QRR), pinned to GJR18 Lemma 5.2 (16), and adds the uniqueness lemma
   `lem:set-Runique` as optional. Recommendation: include the lemma. With it Theorems A and C are statements about an
   explicitly characterised R, and (T) is needed for Corollary B only.
2. **Edge condition (E5).** DNOPS Prop. 4.3 is redundant, as referees U09 and L3 say; the edge condition (2-33) is a
   separate hypothesis and the printed formula is false. A1's reader re-derived the judge's proof and wrote it as a new
   part (iii) of `prop:mod-Rtr` (AMEND_1 A2). Recommendation: add it; the citations of DNOPS Prop. 4.3, the sentence
   before Thm 8.1 and (8-1) then leave (Q).
3. **`eq:con-bz` (E1).** A shorter repair inside power series replaces the judge's ring with logarithm symbols
   (AMEND_1 B2). Needs a careful human read, since it is a new argument.
4. **PPZ over rings (E6).** With the polynomial-identity argument only the statements of Props. 2.12 and 2.14 are
   used; Props. 2.4, 2.8, 2.9 leave (Q) (AMEND_1 B4). This reverses the earlier decision to list Prop. 2.8.
5. **Formal frame change (E7).** New wording says "no analytic continuation is involved"; the note does not prove
   that the expansion differs from the continuation. For the CGYZ comparison, version (a) of AMEND_1 A3.4 adds the
   short membership argument; version (b) restricts the claim. Injectivity of ev (D9) is not added: it needs
   Beukers-Heckman and is used nowhere.

## The declared inputs after the amendments (AMEND_1_chain_4, last section)

Same five: (Q), (QRR), (M), (G0), (T), with [AC96, (1.7)]. (QRR) gains the flatness of R. (T) gains the all-genus CohFT
axioms with flat unit and the all-genus divisor equation. (Q) loses DNOPS (2-33), Prop. 4.3, the sentence before
Thm 8.1, (8-1), and PPZ Props. 2.4, 2.8, 2.9. No new source is cited.

## Not verified

* Source papers: only spot-checks in the audit's stored text extracts. Single-referee pins are listed in AMEND_4.
  Not read by the audit at all: the journal numbering of Teleman and Givental, Hock25, LP20 (17)-(20) and (26), GR19,
  CGYZ26.
* The computations behind `rmk:mai-graphlocal` and `rmk:mai-necessity`; the wording is only softened.
* None of the new LaTeX has been compiled.
