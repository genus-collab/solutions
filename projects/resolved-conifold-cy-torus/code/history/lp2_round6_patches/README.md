# Round-6 patch scripts from the local P2 work (record)

One-off Python scripts written during round 6 of the local P2 regularity work (2026-11-02) that edited documents of this
project in place (`hodge_vertex_conifold/leaf_proof/leaf_structure_proof.md`, `proof_status_mixed_scale.md`,
`conifold_regularity_note.md`, `lemmaA/LEMMA_A_PROOF.md`, `lemmaA/RESUME.md`). They were moved here from
`hodge_vertex/lp2_explore/round6_scripts/` of the local P2 project on 2026-11-27, keeping their subfolders, so that neither
project reaches into the other; only their path line was changed (it now climbs to this project's folder).

All of them have already been applied. They are a record of how the documents were edited and must not be rerun: each
replaces fixed strings and would fail or double-apply on the current text.

Three further round-6 scripts edited documents of both projects and stay in the local P2 project
(`round6c/apply_minor/edit_conifold.py`, `round6c/integrate/patch_conifold.py`, `round6d/repair/apply_edits.py`); they
now take this project's folder from the environment variable `GENUS_RESOLVEDCONIFOLD` and stop if it is not set.
