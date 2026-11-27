# Refined conifold gap and regularity of stable graphs for local P² × C²

**Strand:** solutions
**Status:** Active
**Proposal issue:** #3
**Credited to:** the GENUS collaboration

## Scope

The refined Gromov–Witten theory of K_{P²} × C² on the fibrewise torus, which acts trivially on P², with
ε₁ε₂(ε₁+ε₂) ≠ 0. In every genus g ≥ 2 the project shows (with the status given below) that the generating
series, after the conifold substitution in its polynomial form, equals an explicit sum over stable graphs
with a triple Hodge class at the vertices (Theorem 1.1); that each graph with at least one edge is analytic at
the conifold point, for each colouring, and that the polar part of the one-vertex graph is the single term
(−3)^{g−1} β_{2g−2}(ε₁, ε₂) t_c^{2−2g} (Theorem 1.2). Together these give the refined conifold gap on this
torus (Corollary 1.3). The regularity statement uses the Mariño–Vafa formula and not the mirror theorem; the
identification with Gromov–Witten theory uses the mirror theorem for toric stacks and Teleman's
classification at generic weights. No topological recursion on the mirror curve is used.

Outputs:

- [`output/main.pdf`](output/main.pdf) (94 pages; source `main.tex`, `sections/`, `refs.bib`) and the abridged
  version [`output/main_summary.pdf`](output/main_summary.pdf) (30 pages; source `main_summary.tex`, `summary/`).
- `code/`: R-matrices and the conifold gap in genus 2–4 (notebooks and the twisted R-matrix), the Hodge-vertex
  graph sum in genus 2–4 (genus 4 on a cluster), the per-graph comparison with the resolved conifold, and the
  scripts, reports and referee notes of the six proof rounds. See [`code/README.md`](code/README.md).
- `data/`: where the result tables are, and the large files in
  [genus-collab/data](https://github.com/genus-collab/data). See [`data/README.md`](data/README.md).

Related project: [solutions/resolved-conifold-cy-torus](../resolved-conifold-cy-torus/) (the resolved conifold
× C² on the Calabi–Yau torus). It cites Theorems 8.1 and 9.1 and Corollaries 8.2–8.4 of `main.pdf` in its
Theorem C. The working notes of this project (`code/hodge_vertex/lp2_explore/`) use its Lemma A
(`code/hodge_vertex_conifold/lemmaA/LEMMA_A_PROOF.md` there).

## Contributors

Listed for the record; workshop PDFs carry no individual author names
(git history is the durable record of participation).

- Andrea Brini (University of Sheffield)
- Shuai Guo (Peking University)
- Lei You (Peking University)
- Yannik Schuler (ETH Zurich)

## Methodology

1. **Tools used.**
   - Claude Code (Anthropic), September–November 2026, with the Claude Opus 5 / 5.5 and Claude Fable
     5.1 / 5.2 models, including multi-agent workflows (separate author, critic and referee agents).
   - Mathematica 13.0 (`wolframscript`), for the R-matrices, the graph sums and most checks.
   - Perl with `Math::BigRat`, and Python 3 with exact rational arithmetic (`fractions`), for independent checks.
   - The Stanage HPC cluster (University of Sheffield) for the genus-4 computations.
2. **What they were used for.** The AI tools wrote most of the code, ran the computations, drafted the
   proofs and the LaTeX note, and refereed the proofs in six rounds of independent agent passes (the status
   notes, reports and referee files under `code/hodge_vertex/lp2_explore/`). The research direction and the
   rulings on which inputs are allowed (recorded in `code/hodge_vertex/lp2_explore/ALLGENUS_PATHWAY_2026-11-01.md`
   and the round status notes) were set by A. Brini.
3. **Human verification performed.** No claim has yet been checked by a contributor independent of its
   originator, so no claim is labelled *proved (human-checked)*. The proofs have had AI referee passes only.
   The exact computations in genus 2–4 (Appendix A of `main.pdf`) agree with the statements; they are labelled
   *numerically supported*. Contributors who check a statement: please record it in the table below and
   change its label in the LaTeX.
4. **Archived AI transcripts:** see [`transcripts/README.md`](transcripts/README.md).

There is no stigma attached to extensive AI use, nor to none.

## Claim status

Every theorem, proposition, lemma and corollary of `main.pdf` (102) and `main_summary.pdf` (48) is labelled
*AI-proposed, unchecked*; Appendix A of `main.pdf` (and the evidence part of Section 11 of `main_summary.pdf`)
is labelled *numerically supported*. The main claims:

| Claim (label in `main.tex`) | Status | Verifier / method |
|---|---|---|
| Theorem 1.1, frame formula (`thm:int-frame`) | AI-proposed, unchecked | AI referee passes; graph sum = F_g exactly in genus 2, 3 and 4 (genus 4: 379 stable graphs, symbolically; `code/hodge_vertex/`) |
| Theorem 1.2, regularity of graphs with an edge, value of the one-vertex graph (`thm:int-regularity`) | AI-proposed, unchecked | AI referee passes; checked graph by graph in genus 2, 3, 4 |
| Corollary 1.3, refined conifold gap (`cor:int-gap`); Corollary 1.4, coefficients (`cor:int-gapkh`) | AI-proposed, unchecked | follow from Theorems 1.1 and 1.2; gap checked in genus 2–4 (also via the twisted R-matrix, `code/cluster_genus4/`) |
| Theorem 8.1 and Corollaries 8.2–8.4, vanishing of Hodge integrals (`thm:van-main`, `cor:van-*`) | AI-proposed, unchecked | AI referee passes; used by solutions/resolved-conifold-cy-torus |
| Theorem 9.1, value of the vertex (`thm:val-main`) | AI-proposed, unchecked | AI referee passes (two rounds) |
| Section 11, the resolved conifold | AI-proposed, unchecked | see solutions/resolved-conifold-cy-torus |
| Finite-order evidence (Appendix A) | numerically supported | exact arithmetic over Q, Q(ω) and finite fields; several independent implementations |

Allowed statuses: conjectural · numerically supported · AI-proposed, unchecked ·
proved (human-checked) · formally verified.

## Reproducing

See [`code/README.md`](code/README.md). The PDFs are built with `latexmk -pdf main.tex` and
`latexmk -pdf main_summary.tex`, with `genus.sty` from
[genus-collab/meta](https://github.com/genus-collab/meta/tree/main/latex) in this folder.
