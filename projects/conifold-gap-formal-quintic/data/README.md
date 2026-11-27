# Data

Plain-text formats only (CSV, JSON or similar). Each table needs a header or a
companion `*.md` file stating conventions, normalisations and how it was computed.
Tables of general interest should also be contributed to the `atlas` repo, with
provenance linking back to this project.

## Where the data of this project are

The results stay next to the code that produced and reads them, in `../code/` (JSON and plain text).
Conventions: `../code/hodge_vertex_quintic/OU_EXTENSION_genus3.md` s.1 and `NATURAL_FORM.md`.

| Data | File(s) |
|---|---|
| Genus 2 and 3 graph sums at the conifold, by frame and prime | `../code/hodge_vertex_quintic/*.json` (e.g. `conifold_symbol_result.json`, `g3_nat_*_nat.json`, `g3_state_*.json`), logs in `logs/` |
| Genus 4 (edge splitting, frames and primes) | `../code/hodge_vertex_quintic/g4_split_*.json`, `logs/genus4_split_*.log`, partial sums in `parts/` |
| Genus 5 (Stanage) | `../code/cluster_genus5_quintic/results/`, `parts/`, `logs/`; the audit in `../code/hodge_vertex_quintic/natural_form/genus5_cluster/` (census, controls, independent recomputation) |

## Large files: genus-collab/data

Too large for this repository, in release
[`solutions/conifold-gap-formal-quintic/v1`](https://github.com/genus-collab/data/releases/tag/solutions/conifold-gap-formal-quintic/v1)
of [genus-collab/data](https://github.com/genus-collab/data), with checksums and unpacking instructions in
[`projects/solutions/conifold-gap-formal-quintic/`](https://github.com/genus-collab/data/tree/main/projects/solutions/conifold-gap-formal-quintic)
there:

- `g5_evidence.tar.gz`: `logs/` and `parts/` of the five genus-5 cluster runs; unpack into `../code/g5ev/`;
- `parts_g5_89_nat_local.json.gz`: `hodge_vertex_quintic/natural_form/genus5_cluster/local/out/parts_g5_89_nat_local.json`.
