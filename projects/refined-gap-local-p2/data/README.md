# Data

Plain-text formats only (CSV, JSON or similar). Each table needs a header or a
companion `*.md` file stating conventions, normalisations and how it was computed.
Tables of general interest should also be contributed to the `atlas` repo, with
provenance linking back to this project.

## Where the data of this project are

The result tables stay next to the code that produced and reads them, in `../code/` (Mathematica
`InputForm`, plain text). Conventions: `../code/hodge_vertex/README.txt` and the cluster READMEs.

| Data | File(s) |
|---|---|
| F₂₀, F₃₀ on the general torus and the R-matrix to z⁶ | `../code/F20_F30.m`, `F20_result.m`, `F30_result.m`, `F20redv15.m`, `F30red.txt`, `cPP_z6.m`, `cPPu6.m` |
| Twisted R-matrix and the conifold gap, genus 2–4 | `../code/genRmat_*.m`, `genRmat_*.txt`, `genRmat_gap_results.m`; genus 4: `../code/cluster_genus4/genRmat_F40_cluster_result.m` |
| Hodge-vertex graph sum, genus 2 and 3 | `../code/hodge_vertex/hodge_vertex_genus2_results.m`, `hodge_vertex_genus3_results.m` |
| Hodge-vertex graph sum, genus 4 (379 stable graphs, 27172 decorated graphs) | `../code/hodge_vertex/hodge_F40.m` (total), `hodge_stablegraphs_genus4.m` (per stable graph), `hodge_collect_genus4.txt`, `hodge_vertex_genus4_summary.txt`, `hodge_vertex_genus4_analysis.m`; per decorated graph: data release below |
| Per-graph local P² terms and polar parts, genus 2 and 3 | `../code/localP2_graphterms_genus3.m`, `polar_check_genus2.m` |

## Large files: genus-collab/data

Too large for this repository, in release
[`solutions/refined-gap-local-p2/v1`](https://github.com/genus-collab/data/releases/tag/solutions/refined-gap-local-p2/v1)
of [genus-collab/data](https://github.com/genus-collab/data), with checksums and unpacking instructions in
[`projects/solutions/refined-gap-local-p2/`](https://github.com/genus-collab/data/tree/main/projects/solutions/refined-gap-local-p2)
there:

- `hv_terms_g4.tar.gz`: the 120 genus-4 chunk files `hodge_vertex/stanage_run_2026-10-30/hv_terms_g4_c<k>of120.m`;
- `hodge_graphterms_genus4.m.gz`: `hodge_vertex/hodge_graphterms_genus4.m` (every genus-4 decorated graph);
- `StableGraph-7.bin.gz`: `StableGraph-7.bin`.

Unpack them into `../code/`.
