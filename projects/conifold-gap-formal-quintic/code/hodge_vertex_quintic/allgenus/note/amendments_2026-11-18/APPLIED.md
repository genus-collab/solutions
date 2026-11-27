# Amendments applied to the TeX sources, 2026-11-18

All of `AMEND_1` to `AMEND_4` was applied to `main.tex`, `refs.bib` and `sections/*.tex` (instruction: "apply everything
and rebuild"). The sources as they were before are in `../archive_pre_amend_2026-11-18/` (with the 119-page `main.pdf`).
`main.pdf` now has 126 pages; the build has no undefined reference or citation, no multiply defined label, no overfull box.
All of this is AI work. No source paper was re-read; the corrected statement numbers rest on the audit referees' readings.

## Logs, item by item

| log | list | result |
|---|---|---|
| `APPLIED_1.md` | AMEND_1 (A1-A3, B1-B5, C1-C6, D1-D9) | 65 entries: 62 as written, 2 adapted, D9 not added (by decision) |
| `APPLIED_2.md` | AMEND_2 (sections 1-5) | 39 items: 32 applied, 2 adapted, 4 already done, A2-11 left to A4-05 |
| `APPLIED_3.md` | AMEND_3 (sections 6-10, appendices) | 63 items: 56 applied, 3 adapted, 3 duplicates, A3-54 has no text |
| `APPLIED_4.md` | AMEND_4 (cross-cutting, literature, 18 pins) | 32 items: 18 applied, 9 adapted, 2 already done, 2 duplicates, A4-24 has no text; pins: 12 applied, 2 adapted, 4 already done |

Rejected (R) items of the lists were not applied.

## Choices made where a list offered alternatives

* A1: definition of R kept; flatness declared in (QRR); `lem:set-Runique` added (after `lem:set-hrot`).
* A2: edge condition proved as `prop:mod-Rtr` (iii), `eq:mod-edge`.
* A3.4: version (a), the membership argument.
* `thm:mod-var` (ii): re-proved from `prop:mod-period` (A4-02) instead of weakened (C3.1).
* `def:set-fq`, non-vanishing of the localisation weights: the A2-23 text stands; A4-06 was not substituted.

## Edits made after the four passes (main session)

* `s02_setup_5.tex`, `lem:set-pushins`: "used only when Omega is one of the topological field theories of ...".
* `s10_main_4.tex`, `cor:mai-univ`: pointer "(compare `lem:mai-pole`)" restored (the lemma was otherwise cited nowhere).
* `s10_main_2.tex`: "Def. 0.3, in particular (iii)", as in (Q).
* `s01_intro_3.tex`: the recursion DNOPS (2-26)-(2-29), "to which the recollection of DOSS Thm 4.1 in DNOPS s.2.6 refers".
* `main.tex`: the words "(the holomorphic limit in the conifold frame)" removed from the abstract (the phrase is not
  defined in the note); `\emergencystretch=2em` added.
* `s08_connection_4.tex`: one display split over two lines.

## Open after the amendments (reported by the passes; nothing changed on their account)

New mathematics that no referee has read: `lem:set-Runique`; `prop:mod-Rtr` (iii); the membership argument of A3.4;
the repair of `eq:con-bz` (B2); `thm:mod-var` (ii); the new statement of `cor:mai-univ` (A3-62).

Points to look at:
1. `lem:set-Runique` uses without saying so that the `e~_i` are a basis over Q(zeta)[[q]] and that L is a unit;
   symplecticity is now both quoted in (T) and proved in the lemma.
2. `thm:set-GT`: the all-genus divisor equation is a quoted input, but its parenthesis reads as a proof sketch and uses
   `E_{g,n,d}`, which the note defines in genus 0 only.
3. A2-04: the algebra B is called a differential Q[[q]]-algebra although its coefficients are in Q(zeta)[[q]].
4. A2-05: "degree at most 1 in l over K" needs l transcendental over K.
5. A2-25: the restriction 2h-2+n > 0 was added to one of two parallel formulas only.
6. A3-03: uniformity of the bound in `prop:mod-V1` in the circle variables and in v is not stated.
7. A3-25: `prop:gau-graded` assumes (EP), which `cor:mai-univ` does not restate.
8. A4-01: `thm:set-mirror` now refers forward to `thm:set-G0` (b). A4-17: `s06_model_2.tex` uses `R_TR` before it is
   introduced. A4-23: "whose closed form" has no clear antecedent.
9. Pins typed as given from a single referee: DNOPS (2-28); DOSS (3-57), (3-64), (4-1), (2-23), (2-26); ABDKS Thm 1.14
   and the journal numbers; EO (4-7), (5-27); GJR18 Cor. 7.6 and (37); DLMF 5.11.8.
10. The computations behind `rmk:mai-graphlocal` and `rmk:mai-necessity` are still unrecorded.
11. Draft notes left: 10 (title, abstract, D11, version re-checks, journal numbering of Givental and Teleman, remarks).

## Referee pass on the new mathematics (2026-11-18, later the same day)

Five AI referees read the six new pieces; reports in `referee/REF_1` to `REF_5`, fixes applied (log `referee/APPLIED_REF.md`,
30 findings, 41 replacements). `main.pdf` rebuilt: 127 pages, clean.

| piece | verdict | checked by computation |
|---|---|---|
| `lem:set-Runique` | sound; 5 wording fixes | recursion and symplecticity, generic data, through z^5, q^6 |
| `prop:mod-Rtr` (iii), `eq:mod-edge` | sound; wording (the DNOPS misprint described exactly; coordinates stated) | all 25 pairs, degree <= 5, five parameter sets, floating point |
| membership argument (`sec:sgn-lit`) | sound; wording (admissibility of the pair is quoted from CGYZ26, not proved) | weights and identities, exact |
| `eq:con-bz`, `eq:con-factor` | sound; 6 wording fixes | exact series to total order 13 at rational L |
| `thm:mod-var` (ii) | sound; pins of EO corrected | h = 1..5 exact |
| `cor:mai-univ` | one minor error fixed: lambda := -b_0(0) b_3(0)/v_A'(0)^2 (was without the sign, giving -1 for the quintic) | hand checks |

Still open: the CGYZ26 pins (Def. 1.8 (2), Ex. 1.3) rest on an automatic reading of the HTML of version 2; points 2-10 of the
list above, except those covered by the referee fixes (point 1).
