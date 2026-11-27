# Refined Gromov–Witten theory of the resolved conifold on the Calabi–Yau torus: conifold regularity and rigidity

**Strand:** solutions
**Status:** Active
**Proposal issue:** #1
**Credited to:** the GENUS collaboration

## Scope

The equivariant Gromov–Witten theory of the resolved conifold times C² on its Calabi–Yau four-torus, written
as a Givental–Teleman graph sum with the normalised Hodge class of weights (E, −ε₁, −ε₂) at every vertex. At
the conifold point Q = 1 one of the two idempotents is singular. The project proves (with the status given
below) that every factor of the graph sum at the singular idempotent has a leaf structure (Theorem A); that the
graph sum equals the Gromov–Witten potential on the Calabi–Yau torus (Theorem B); that every coloured stable
graph with an edge is regular at the conifold point, the whole polar part sitting in the single-vertex graph
(Theorem C); and that for g ≥ 2, ε₁ε₂ GW_{g,d} = P_g d^{2g−3} for every d ≥ 1 (Theorem D). With genus 0 and 1
this gives the rigidity conjecture of Brini–Schuler (arXiv:2410.00118) for the resolved conifold. The
constant-map term is not determined.

Outputs:

- [`output/main.pdf`](output/main.pdf) (100 pages; source `main.tex`, `sections/`, `refs.bib`) and the abridged
  version [`output/main_summary.pdf`](output/main_summary.pdf) (31 pages; source `main_summary.tex`, `summary/`).
- `code/`: the Hodge-vertex graph sum in genus 2–4 (Mathematica), the genus-4 cluster package, the checks of
  Lemma A, of the leaf structure and of the rigidity statements (Mathematica, Perl, Python), and the
  working notes of the proofs. See [`code/README.md`](code/README.md).
- `data/`: where the result tables are, and the large files in
  [genus-collab/data](https://github.com/genus-collab/data). See [`data/README.md`](data/README.md).

Related project: the companion project on local P² × C² (cited as [GEN26] in the PDF) proves the two
theorems on Hodge integrals used in Theorem C (Theorems 6.8 and 6.9 of `main.pdf`), and itself uses Lemma A of
this project (`code/hodge_vertex_conifold/lemmaA/LEMMA_A_PROOF.md`). It is to be added to the GENUS
repositories; the link will be added here.

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
   - Mathematica 13.0 (`wolframscript`), for the graph sums and most checks.
   - Perl with `Math::BigRat`, and Python 3 with exact rational arithmetic (`fractions`), for independent checks.
   - The Stanage HPC cluster (University of Sheffield) for the genus-4 graph sum (120 workers, about 14 CPU-hours).
2. **What they were used for.** The AI tools wrote most of the code, ran the computations, drafted the
   proofs and the LaTeX note, and refereed the proofs in several rounds of independent agent passes (the
   `report.md` / `referee*.md` files under `code/hodge_vertex_conifold/`). The research direction, the choice of
   frame and the rulings on which inputs are allowed (recorded in `code/hodge_vertex_conifold/lemmaA/RESUME.md`
   and the status notes) were set by A. Brini.
3. **Human verification performed.** No claim has yet been checked by a contributor independent of its
   originator, so no claim is labelled *proved (human-checked)*. The proofs have had AI referee passes only.
   The exact computations in genus 2–4 (Appendix A of `main.pdf`) agree with the statements; they are labelled
   *numerically supported*. Contributors who check a statement: please record it in the table below and
   change its label in the LaTeX.
4. **Archived AI transcripts:** see [`transcripts/README.md`](transcripts/README.md).

There is no stigma attached to extensive AI use, nor to none.

## Claim status

Every theorem, proposition, lemma and corollary of `main.pdf` (114) and `main_summary.pdf` (67) is labelled
*AI-proposed, unchecked*; Appendix A of `main.pdf` (and Section 8 of `main_summary.pdf`) is labelled
*numerically supported*. The main claims:

| Claim (label in `main.tex`) | Status | Verifier / method |
|---|---|---|
| Theorem A, leaf structure (Thm 1.1 = Thm 4.1, `thm:int-leaf`, `thm:lea-main`) | AI-proposed, unchecked | AI referee passes; exact checks of the leaf coefficients (Appendix A.2–A.3) |
| Theorem B, graph sum = potential (Thm 1.2 = Prop 5.1, `thm:int-ident`, `prop:ide-main`) | AI-proposed, unchecked | AI referee passes; graph sum = F_g exactly in genus 2–4 (`code/hodge_vertex_conifold/`) |
| Theorem C, regularity and polar part (Thm 1.3 = Cor 6.13, `thm:int-reg`, `cor:reg-main`) | AI-proposed, unchecked | AI referee passes; regularity of every graph with an edge checked exactly in genus 2, 3, 4 (379 stable graphs in genus 4) |
| Theorems 6.8, 6.9 (`thm:reg-L`, `thm:reg-VI`), cited from the local P² × C² project | AI-proposed, unchecked | status as in that project |
| Lemma A (`code/hodge_vertex_conifold/lemmaA/LEMMA_A_PROOF.md`) | AI-proposed, unchecked | two AI referee passes (round 6) |
| Theorem D, rigidity (Thm 1.4, `thm:int-rig`) | AI-proposed, unchecked | AI referee passes; closed form checked exactly in genus 2–4 |
| Corollary 1.5, rigidity conjecture of [BS] for the resolved conifold (`cor:int-conj`) | AI-proposed, unchecked | follows from Theorem D and genus 0, 1 |
| Finite-order checks (Appendix A) | numerically supported | exact rational arithmetic and computations over F_p, several independent implementations |

Allowed statuses: conjectural · numerically supported · AI-proposed, unchecked ·
proved (human-checked) · formally verified.

## Reproducing

See [`code/README.md`](code/README.md). The PDFs are built with `latexmk -pdf main.tex` and
`latexmk -pdf main_summary.tex`, with `genus.sty` from
[genus-collab/meta](https://github.com/genus-collab/meta/tree/main/latex) in this folder.
