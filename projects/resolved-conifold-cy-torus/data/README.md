# Data

Plain-text formats only (CSV, JSON or similar). Each table needs a header or a
companion `*.md` file stating conventions, normalisations and how it was computed.
Tables of general interest should also be contributed to the `atlas` repo, with
provenance linking back to this project.

## Where the data of this project are

The result tables stay next to the code that produced and reads them, in `../code/` (Mathematica
`InputForm`, plain text). Conventions: `../code/hodge_vertex_conifold/README.txt` (weights, frames A/B,
modes, file formats).

| Data | File(s) |
|---|---|
| F₁₁, F₂₀, F₃₀ on the general five-torus (exact rational functions of Q and the weights) | `../code/resolvedconifold_results.m` (also `../code/hodge_vertex_conifold/ref/`) |
| Hodge-vertex graph sum, genus 2 and 3, modes A, B and point: per stable graph and per decorated graph | `../code/hodge_vertex_conifold/hvc_results_g{2,3}_{Am,Bm,point}.m`, `hvc_terms_g*_c*.m`, tables `hvc_table_*.txt`, `hvc_summary_*.txt`, `hvc_analysis_*.txt` |
| Same, genus 4, mode B (379 stable graphs, 27172 decorated graphs) | `../code/hodge_vertex_conifold/stanage_genus4/hvc_results_g4_Bm.m` (per stable graph; per decorated graph without the term), `hvc_table_g4_Bm.txt`, `hvc_summary_g4_Bm.txt`, `hvc_analysis_g4_Bm.txt` |
| Per-graph conifold expressions and CY-torus polar parts, genus 2 and 3 (point vertices) | `../code/conifold_graphterms_genus2_3.m`, `../code/conifold_CY_graphterms_polar.m` |
| Lemma A data (Hodge integrals, moments) | `../code/hodge_vertex_conifold/lemmaA/data/`, `lemmaA/E_formulas/` |

## Large files: genus-collab/data

Too large for this repository, in release
[`solutions/resolved-conifold-cy-torus/v1`](https://github.com/genus-collab/data/releases/tag/solutions/resolved-conifold-cy-torus/v1)
of [genus-collab/data](https://github.com/genus-collab/data), with checksums and unpacking instructions in
[`projects/solutions/resolved-conifold-cy-torus/`](https://github.com/genus-collab/data/tree/main/projects/solutions/resolved-conifold-cy-torus)
there:

- `hvc_terms_g4_Bm.tar.gz`: the 120 genus-4 chunk files `stanage_genus4/hvc_terms_g4_Bm_c<k>of120.m`
  (the term of every decorated graph);
- `engine_h5_r7.txt.gz`: `lemmaA/E_formulas/engine_h5_r7.txt`.

Unpack them into `../code/hodge_vertex_conifold/`.
