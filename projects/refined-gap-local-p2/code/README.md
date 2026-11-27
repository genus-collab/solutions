# Code

- **Language and version:** Mathematica 13.0 (`wolframscript` and notebooks); Perl 5 with `Math::BigRat`;
  Python 3 (standard library, exact `fractions`); bash and SLURM for the genus-4 cluster runs.
- **Dependencies / packages:** none beyond the above. Stable graphs are read from `StableGraph-<g>.bin`
  (genus 2–6 here; `StableGraph-7.bin`, 72 MB, is in the data release, see `../data/README.md`).
- **Folder layout.** This folder keeps the layout of the original working folder, so relative paths inside it
  work as they are:

  | item | what it is |
  |---|---|
  | `localp2final_CY5_v15.nb` | genus-2 base of everything below (R-matrix coefficients cPP, derivations) |
  | `localp2final_CY5_v16.nb`, `localp2CY5_genus3.nb`, `localp2CY5_optimised.nb` | genus 3: R-matrix to z^6, the genus-3 pipeline and checks; data `cPP_z6.m`, `cPPu6.m`, `PFu6.m`, `v15data.m`, `moncycproj.m`, `F20_*`, `F30_*`, `F20_F30.m`, `F20redv15.m`, `F30red.txt`; logs `genus3_checks.log`, `partA_*.log`, `partB_run_g2.log` |
  | `localp2CY5_genRmat.nb`, `localp2CY5_genRmat_v3.nb`, `localp2CY5_genRmat_genus3.nb`, `run_genRmat_genus3.wls`, `genRmat_*` | twisted (generalised) R-matrix and the conifold gap in genus 2–4; `v3_F20_twisted_reference.m` |
  | `cluster_genus4/` | genus-4 twisted R-matrix on the cluster: the conifold gap fixes q1, q2, q3 (`README.txt`) |
  | `hodge_vertex/` | Hodge-vertex (Guo–Lam) graph sum, genus 2–4 (`README.txt` there: conventions, results, files); notebook `localp2CY5_hodge_genus3.nb`; genus-4 results `hodge_F40.m`, `hodge_stablegraphs_genus4.m`; the cluster run `stanage_run_2026-10-30/` (logs, SLURM outputs, per-chunk partial sums; genus 2 and 3 term files) |
  | `hodge_vertex/lp2_explore/` | the all-genus programme (rounds 1–6): `FRAMEB_PROOF.md`, `LP2_REGULARITY_PROOF_draft*.md`, `ROUND*_STATUS_*.md`, `ALLGENUS_PATHWAY_2026-11-01.md`, and the scripts, reports and referee notes of each round (`round*_scripts/`) |
  | `hodge_vertex/compare_conifold/` | decorated-graph comparison with the resolved conifold × C² (genus 2–4), with frozen copies of the conifold results in `ref_conifold/` (provenance README there) |
  | `cluster_genus4_hodge/` | genus-4 Hodge-vertex cluster package (`README.txt`) and the copied-back worker logs |
  | `cluster_genus4_polar/`, `cluster_genus4_polar_validation_2026-10-29/` | genus-4 per-graph polar parts of local P² and of the resolved conifold, their comparison, and the validation logs of the cluster run |
  | `localP2_graphterms_genus3.m`, `polar_check_*` | per-graph local P² terms and the 2026-09 comparison with the resolved conifold (genus 2, 3; CY torus) |

  The resolved-conifold data used in the comparisons are frozen copies from the GENUS project
  [solutions/resolved-conifold-cy-torus](../../resolved-conifold-cy-torus/). The LaTeX note that lived in
  `hodge_vertex/lp2_explore/paper/` is now `main.tex`, `sections/`, `main_summary.tex` and `summary/` at the top
  of the project; its working notes are in `../drafting/`. The `*.md` files are working notes, AI reports and
  referee reports; `*.out`, `*.log` and `*_out.txt` files are stored outputs.

  Three round-6 scripts in `hodge_vertex/lp2_explore/round6_scripts/` (`round6c/apply_minor/edit_conifold.py`,
  `round6c/integrate/patch_conifold.py`, `round6d/repair/apply_edits.py`) also edited documents of the
  resolved-conifold project; they read its location from the environment variable `GENUS_RESOLVEDCONIFOLD`.
  All round-6 patch scripts are already applied and must not be rerun.

- **Not included here:**
  - Mathematica `.mx` files (binary `DumpSave` files, tied to the Mathematica version). The 13 at the top of
    this folder that the notebooks read are included as plain-text `.m` files with the same name, written with
    `Save` and checked to give identical definitions when read back: `Get["name.m"]` instead of
    `Get["name.mx"]`. The others (about 7,600: the sample caches `hv_samples_*`, `hv3_samples_B.mx`, the orbit
    tables `hv_orbits_*` and the cluster chunk files `hv_chunk_*`) are intermediate results, rebuilt by the
    `prep` / `sample` stages of the scripts; the notebooks recompute them when they are missing.
  - The 120 genus-4 term files `hodge_vertex/stanage_run_2026-10-30/hv_terms_g4_c<k>of120.m`,
    `hodge_vertex/hodge_graphterms_genus4.m` and `StableGraph-7.bin`: see `../data/README.md`.
  - `localp2CY5_genRmat_2.nb`, which is damaged (blocks of binary garbage, from before the project moved), and
    a duplicate of `hodge_vertex/localp2CY5_hodge_genus3.nb` that sat at the top of the folder.
  - Text extracted from published papers that the agents used for reference (arXiv:math/0306434,
    Liu–Liu–Zhou; Okounkov–Pandharipande, Hodge integrals and invariants of the unknot, Geom. Topol. 8 (2004); arXiv:1407.7800, Harnad–Orlov): files
    `round5_scripts/MVmatch/lit/{llz,op}.txt`, `round5_scripts/ho.txt`, `round5_scripts/VI_ref1/llz*.txt`.
    Scripts that compare against these texts need a local copy of the papers.
  - Snapshot tarballs of the working folder.

- **To reproduce the main outputs** (from `hodge_vertex/`; see `hodge_vertex/README.txt` for timings):

  ```
  wolframscript -file hv_genus2.wls checks
  wolframscript -file hv_genus2.wls calib
  wolframscript -file hv_genus2.wls sample point
  wolframscript -file hv_genus2.wls sample A
  wolframscript -file hv_genus2.wls sample B
  wolframscript -file hv_genus2.wls recon
  wolframscript -file hv_genus3.wls checks
  wolframscript -file hv_genus3.wls prep
  wolframscript -file hv_genus3.wls sample
  wolframscript -file hv_genus3.wls recon
  ```

  (`sample` of genus 3 can be split over workers as `sample <w> <nw>`.)

  Genus 4 needs a cluster: see `cluster_genus4_hodge/README.txt` (Hodge-vertex graph sum) and
  `cluster_genus4/README.txt` (twisted R-matrix). The notebook `hodge_vertex/localp2CY5_hodge_genus3.nb` runs
  headless with `wolframscript -file run_nb_headless.wls localp2CY5_hodge_genus3.nb` (from `hodge_vertex/`);
  its genus-4 section needs `hodge_graphterms_genus4.m` from the data release. Checked on 2026-11-27 on this
  folder without any `.mx` file and with the data release unpacked: all 30 input cells run without messages in
  about 7 minutes (the genus-3 samples are recomputed) and give the stored outputs of `nb_headless_g4.out`, except
  that the second half of cell 13 (comparison of the genus-4 evaluator with the genus-3 one on the decorated
  graphs of genus 3) is skipped: it runs only when `hv_orbits_g3.mx` exists or the genus-3 cells have already
  been evaluated.

Mathematica notebooks (`.nb`) are accompanied by a plain-text export (`.m`, input cells as code, text cells as
comments, outputs omitted), so that changes are diff-readable.
