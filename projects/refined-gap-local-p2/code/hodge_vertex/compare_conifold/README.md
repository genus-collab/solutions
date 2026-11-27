# Local P2 versus resolved conifold, decorated graph by decorated graph (mode B)

Compares the polar parts of the Hodge-vertex graph sum of K_P2 x C^2 (conifold frame, t_c at 1 + 27 y = 0; results in `..`)
with those of the resolved conifold x C^2 (t = Log Q at Q = 1), decorated graph by decorated graph, at genus 2, 3, 4.

Until 2026-11-27 this comparison lived in the resolved-conifold project (`hodge_vertex_conifold/`) and read the local P2
results from this folder by a relative path across the two projects. It was moved here so that both projects are standalone:
this project now holds the comparison and frozen copies of the conifold results it needs (`ref_conifold/`, with provenance).

| file | what it is |
|---|---|
| `hvc_compare_localP2.wls` | the comparison; reads `../hodge_vertex_genus2_results.m`, `../hodge_vertex_genus3_results.m`, `../hodge_graphterms_genus4.m` and `ref_conifold/hvc_results_g<g>_Bm.m` |
| `lp2_tc_section.m` | section 6 of `../lp2_lib.m` (conifold-frame t_c expansion: `tcSetup`, `polarOfTerm`, `gapCoefficient`), used at genus 2 |
| `hvc_compare_localP2_g2.txt`, `_g3.txt`, `_g4.txt` | stored outputs |
| `nb_hvc_compare.m` | the comparison cells formerly in the resolved-conifold notebook `resolvedconifold_hodge.nb` (record; they need the conifold results loaded in the notebook session and do not run on their own) |
| `ref_conifold/` | frozen copies of the conifold mode-B results, see its README |

Usage (from any directory): `wolframscript -file hvc_compare_localP2.wls <genus>`, genus 2, 3 or 4 (about 2 s, 3 s, 40 s).
Rerun on 2026-11-27 from this folder: all three outputs identical to the stored ones.
