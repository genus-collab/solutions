# Frozen copies of resolved-conifold results

Copied on 2026-11-27 from the resolved-conifold project (resolved conifold x C^2, Hodge-vertex graph sum, mode B), so that the
comparison in `..` runs inside this project. They are reference data for the comparison only, not results of this project;
the originals and the code that produced them stay in the resolved-conifold project.

| file | copied from (resolved-conifold project) | bytes | sha256 (first 16) |
|---|---|---|---|
| `hvc_results_g2_Bm.m` | `hodge_vertex_conifold/hvc_results_g2_Bm.m` | 33292 | dd035ab0f04744be |
| `hvc_results_g3_Bm.m` | `hodge_vertex_conifold/hvc_results_g3_Bm.m` | 2727770 | c63ea83afe7e4953 |
| `hvc_results_g4_Bm.m` | `hodge_vertex_conifold/stanage_genus4/hvc_results_g4_Bm.m` (genus-4 cluster run) | 17860474 | 40427a87223bf4bf |

The genus-4 comparison previously read the binary `hvc_results_g4_Bm.mx` (219 MB); on 2026-11-27 the comparison was rerun
from the plain-text `.m` copy here and reproduced `../hvc_compare_localP2_g4.txt` byte for byte.
