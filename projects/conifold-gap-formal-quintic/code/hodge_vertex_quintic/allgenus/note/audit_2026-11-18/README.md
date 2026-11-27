# Audit of the all-genus proof (2026-11-17/18)

Independent referee audit of the all-genus argument of the note in `../` (`main.tex`, `sections/`, the 119-page build).
Scope: the all-genus proof only; no genus 2-5 computation was rerun. All referees were AI agents (Claude Fable 5.2 in run 1,
Claude Opus 5.5 in runs 2-3); none of this is a human check.

## Verdict

No fatal error and no gap that could not be closed was found, conditional on the declared inputs. Three points need a
text change to make the declared inputs and the claims literally accurate; five places omit a routine step (repairs are
written out); the rest is wording, references and notation.

**Must fix**
1. *Undeclared input* (`verify/E4.md`): property (b) of `def:set-rows` (the columns R e_i e^{u_i/z} solve the QDE) is
   called a convention but is used on the chain of Theorems A, B, C (`prop:set-cyclicLR` (iii), `cor:set-cyclicPF5` (b),
   `lem:gau-colrows`, `thm:gau-R0` (iii), `thm:con-glue`). It is not in Teleman Thm (2) of (1.15) nor GJR18 Thm 4.1;
   pin it to GJR18 Lemma 5.2 (16) and Teleman v3 Props 7.7, 8.2(ii), 8.4, 8.7(ii). Add to (T) the all-genus divisor
   equation and the all-genus CohFT axioms with flat unit (GJR18 Sec. 4.1).
2. *Misprinted source formula* (`verify/E5.md`): DNOPS (2-33) is false as printed; Step 1 of `thm:mod-tr`
   (`sections/s06_model_3.tex` l.28) leans on it. The correct identity (DOSS (4-1)) holds for the model curve; the
   escalation judge gives a short proof. Also: DNOPS (8-1) tuple is a_1..a_5; Prop. 4.3 (2), not (1); DNOPS Sec. 2.3
   omits 1/|Aut|. (Referees U09 and L3 think the note's direct proof of (2-32) already makes the condition redundant.)
3. *Wording of what is proved* (`verify/E7.md`): the conifold expansion is a formal frame change, not analytic
   continuation. Reword the abstract, `def:set-points` l.28 and `sec:sgn-lit` l.106-108 (the "not needed" at l.108 is
   wrong for the CGYZ comparison, but the needed membership holds by a short argument using Theorem C). E7 and X4 also
   argue that ev is injective (differential Galois group Sp_4), which would make Corollary B intrinsic; external theorems,
   not needed by the note.

**Routine missing steps (repairs in the files named)**
- `lem:fsp-res` (iii), exponential under X -> X + chi, used by `prop:con-res` (`verify/E3.md`).
- `eq:con-bz` proved "factor by factor" in a ring that does not contain the factors (`verify/E1.md`).
- Appendix D scope: undefined q-dependent Gaussian functional (`cor:fsp-movebase`); `lem:fsp-param` (i) false over
  subrings not containing Q; `lem:fsp-cov`/`cfu` proved only for power-series amplitudes (`verify/E2.md`).
- PPZ applied over rings other than Q: the "holds verbatim" sentence is not a proof; replacement given (`verify/E6.md`).
- Galois descent in `def:mai-Fcf` cites `lem:gau-colrows` for rationality; cite `thm:con-glue` (`reports/X4.md`, `U06_*`).

**False as written, used by no proof**: `rmk:set-qrrsign` (ph_1 forced independently of (QRR)); `cor:gau-cycsys` (ii)
literal reading; `thm:mod-var` (ii) uses EO (5-15), declared unused; `rmk:mai-graphlocal`, `rmk:mai-necessity` and
`s01_intro_2.tex` l.25 rest on computations not in the evidence appendix.

## Files

- `reports/`: referee reports.
  - Units (two lenses A = line by line, B = adversarial; C = one combined referee):
    U01 setup I, U02 setup II, U03 App. A/B + s03_gauge, U04 Theorem U, U05 gauge identity/canonical gauge,
    U06 WKB/colours/thm:gau-R0, U07 transport, U08 App. C + Theorem R, U09-U11 the model (TR), U12 comparison,
    U13-U14 Theorem P and the glue, U15 Theorem Gamma/C, U16-U17 App. D (U1617B = second referee for all of App. D),
    U18 main theorem, Corollary B, Complements. `U03B`, `U05B`, `U07B` are the second referees of run 3.
  - Cross-cutting: X1 circularity, X2 interfaces (incl. the constant c_g), X3 conventions and signs, X4 rings and
    expansions, X5 ledger of inputs and unproved assertions.
  - Literature (checked against the arXiv versions of `../refs.bib`): L1 PPZ, L2 Eynard-Orantin, L3 DOSS/DNOPS,
    L4 ABDKS x-y swap, L5 QRR/Teleman/G0, L6 mirror theorem and remaining formulas.
- `verify/`: E1-E8, independent judges of the minor findings that touch the proof chain (verdict, repair, text changes).
- `data/`: structured results of the three runs (`done_wave1.json`, `wave2_result.json`, `wave3_result.json`),
  `wave2_digest.txt` (all summaries and minor findings of run 2), and the label dependency index
  (`deps.txt`, `items.json`, made by `depgraph.py`).
- `work/`: the referees' own small check scripts and outputs (wolframscript / Python); downloaded PDFs of the sources
  were not kept.
