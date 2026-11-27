# Conifold regularity of the formal quintic, graph by graph

**Strand:** solutions
**Status:** Active
**Proposal issue:** #5
**Credited to:** the GENUS collaboration

## Scope

The formal quintic of Lho and Pandharipande (the O(5)-twisted equivariant Gromov–Witten theory of P⁴
specialised at the fifth roots of unity), expanded at the conifold point u = 1 − 5⁵q = 0 in a conifold frame.
The expanded R-matrix is factorised as R = R′ M₀(τ(u) ẑ), with M₀ the R-matrix of the Landau–Ginzburg model
p⁴/4 + v/p (primitive form dp) and τ(u)ẑ the natural variable. The project shows (with the status given below)
that for every g ≥ 2, in every conifold frame of a two-parameter family, the term of every leg-free stable
graph of genus g with at least one edge is regular at the conifold point, and the single-vertex graph
contributes c_g t_c^{2−2g} plus a regular series, c_g = (−1)^{g−1} B_{2g}/(2g(2g−2) 5^{g−1}) (Theorem A).
Summing over graphs gives the conifold gap for this expansion of the formal quintic (Corollary B).

This project covers graph regularity and the conifold gap only. The later work on the regular terms of F_g at
the conifold (the exploration of 23–27 November 2026 and the genus 4–6 cluster runs for them) is not part of it.

Outputs:

- [`output/main.pdf`](output/main.pdf) (116 pages; source `main.tex`, `preamble.tex`, `macros.tex`, `sections/`,
  `refs.bib`) and the abridged version [`output/main_summary.pdf`](output/main_summary.pdf) (30 pages; source
  `main_summary.tex`, `summary/`).
- `code/`: the genus 2–5 graph sums at the conifold (genus 5 on a cluster), the natural factorisation and
  Theorem U, the all-genus mechanism (stages 0–5), and the working material of the note (design, gaps, reviews,
  the audit and amendments of 18 November 2026). See [`code/README.md`](code/README.md).
- `data/`: where the results are, and the large files in [genus-collab/data](https://github.com/genus-collab/data).
  See [`data/README.md`](data/README.md).

Related projects: [solutions/refined-gap-local-p2](../refined-gap-local-p2/) and
[solutions/resolved-conifold-cy-torus](../resolved-conifold-cy-torus/), where the analogous factorisation is the
Hodge one; the natural form here is the analogue of that Hodge trick.

## Contributors

Listed for the record; workshop PDFs carry no individual author names
(git history is the durable record of participation).

- Andrea Brini (University of Sheffield)
- Shuai Guo (Peking University)
- Lei You (Peking University)
- Yannik Schuler (ETH Zurich)

## Methodology

1. **Tools used.**
   - Claude Code (Anthropic), November 2026, with the Claude Opus 5 / 5.5 and Claude Fable 5.1 / 5.2
     models, including multi-agent workflows (separate worker, referee, critic and auditor agents).
   - Python 3 (standard library; exact arithmetic over Q and modular arithmetic modulo 2^61 − 1 and 2^89 − 1,
     one referee also 2^521 − 1) for the graph sums and most checks; Mathematica 13 (`wolframscript`) for
     independent checks.
   - The Stanage HPC cluster (University of Sheffield) for the genus-5 graph sums.
2. **What they were used for.** The AI tools wrote the code, ran the computations, drafted the arguments and
   the LaTeX note, and checked them in several rounds of independent agent passes: referees and critics for each
   stage of the all-genus mechanism, a review of the note, and an audit of every statement and citation
   (18 November 2026) with amendments. The research direction and the rulings on which inputs are allowed were
   set by A. Brini.
3. **Human verification performed.** No claim has yet been checked by a contributor independent of its
   originator, so no claim is labelled *proved (human-checked)*. The proofs have had AI referee and audit passes
   only. The exact and modular computations in genus 2–5 (Appendix F of `main.pdf`) agree with the statements;
   they are labelled *numerically supported*. Contributors who check a statement: please record it in the table
   below and change its label in the LaTeX.
4. **Archived AI transcripts:** see [`transcripts/README.md`](transcripts/README.md).

There is no stigma attached to extensive AI use, nor to none.

## Claim status

Every theorem, proposition, lemma and corollary of `main.pdf` (163) and `main_summary.pdf` (53) is labelled
*AI-proposed, unchecked*; Appendix F of `main.pdf` (and Section 10 of `main_summary.pdf`) is labelled
*numerically supported*. The main claims:

| Claim (label in `main.tex`) | Status | Verifier / method |
|---|---|---|
| Theorem A, regularity of every graph with an edge and the vertex term c_g t_c^{2−2g} (Thm 1.1 = Thm 10.1, `thm:int-main`, `thm:mai-main`) | AI-proposed, unchecked | AI referee, critic and audit passes; checked graph by graph in genus 2–5 (`code/`, Appendix F) |
| Corollary B, the conifold gap of the formal quintic (Cor 1.2 = Cor 10.4, `cor:int-gap`, `cor:mai-gap`) | AI-proposed, unchecked | follows from Theorem A; gap checked in genus 2–5 |
| Theorem C, the normalisation is canonical (Thm 1.3 = Thm 9.7, `thm:int-D`, `thm:gam-D`) | AI-proposed, unchecked | AI referee passes; the finite identity checked for 2 ≤ g ≤ 9 |
| Theorem Γ (Thm 1.4 = Thm 9.4, `thm:int-gamma`, `thm:gam-main`) and Theorem LMI (Thm 9.3) | AI-proposed, unchecked | AI referee passes (stage 5 of the mechanism) |
| Theorem R, reduction to the model (Thm 5.4, `thm:red-main`) | AI-proposed, unchecked | two refereed write-ups (stage 0–1) |
| Finite computations (Appendix F) | numerically supported | exact over Q and Q(ζ), and modulo the primes named (a zero test modulo a prime is probabilistic) |

Allowed statuses: conjectural · numerically supported · AI-proposed, unchecked ·
proved (human-checked) · formally verified.

## Reproducing

See [`code/README.md`](code/README.md). The PDFs are built with `latexmk -pdf main.tex` and
`latexmk -pdf main_summary.tex`, with `genus.sty` from
[genus-collab/meta](https://github.com/genus-collab/meta/tree/main/latex) in this folder.
