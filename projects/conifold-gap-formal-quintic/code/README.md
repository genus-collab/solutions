# Code

- **Language and version:** Python ≥ 3.8 (standard library only; exact arithmetic over Q and modular arithmetic),
  Mathematica 13 (`wolframscript`, `.wl` / `.wls`) for independent checks, bash and SLURM for the genus-5 cluster runs.
- **Dependencies / packages:** none beyond the above. Stable graphs are read from
  `hodge_vertex_quintic/graph_data/StableGraph-{2,3,4}.bin`; the genus-5 package ships its own graph lists.
- **Folder layout.** This folder keeps the layout of the original working folder, so relative paths inside it
  work as they are:

  | item | what it is |
  |---|---|
  | `hodge_vertex_quintic/` | the genus 2–4 programme. Notes of record: `NATURAL_FORM.md` (the natural form M = M₀(τ(u) ẑ); Conjecture N, edge graphs regular and vertex term c_g t_c^{2−2g}, checked for g = 2–5), `FORMAL_QUINTIC_genus2.md`, `CONIFOLD_SYMBOL_genus2.md`, `OU_EXTENSION_genus{2,3}.md`, the earlier scans (`DIAG_`, `SCALAR_`, `SYMPLECTIC_SCAN_genus2.md`) and plans. Code `fq*.py`, `gsum_q.py`, drivers `genus*_*.py`; results `*.json`; logs `logs/`; partial sums `parts/`; independent referees `referee/`, `referee_u3/`, `referee_u3b/` |
  | `hodge_vertex_quintic/natural_form/` | the natural form: Landau–Ginzburg model (`lg/`), Theorem U (`theory/UNIFORM_REDUCTION.md`), referee and critic (`referee/`, `critic/`, `SYNTHESIS.md`), verifiers (`verify_*`), the genus-4 referee, the review of the genus-5 package and the audit of the genus-5 evidence (`genus5_cluster/`) |
  | `hodge_vertex_quintic/allgenus/` | the all-genus mechanism for Conjecture N, stages 0–5: entry points `MECHANISM.md`, `MECHANISM_STAGE2.md`, `MECHANISM_STAGE4.md`, `MECHANISM_STAGE5.md`; worker notes, referee reports, critic code |
  | `hodge_vertex_quintic/allgenus/note/` | the working material of the LaTeX note (now `main.tex`, `sections/`, `summary/` at the top of the project): outlines, labels, notation, bibliography plan, `design/`, `gaps/`, `review/`, the audit of 18 November 2026 (`audit_2026-11-18/`: reports per statement, work folders) and its `amendments_2026-11-18/`, and two archived earlier versions (`archive_*`) |
  | `cluster_genus5_quintic/` | the genus-5 Stanage package (`README.txt`: what is computed, files, steps), with its validation data, logs and results |

  The `*.md` files are working notes, AI reports, referee reports and syntheses; `*.log`, `*.out` and
  `*.console` files are stored outputs.

- **Not included here:**
  - The exploration of 23–27 November 2026 on the regular conifold coefficients of the edge graphs
    (`allgenus/note/explore_2026-11-23/`) and the genus 4–6 cluster packages for the regular terms
    (`cluster_ladder_*`): outside the scope of this project.
  - Copies of published papers that the agents of the audit downloaded to check citations (arXiv sources,
    HTML pages and text extractions in `audit_2026-11-18/work/{L1,L2,L3,L4,L5,L6,U09B,U09_B,U10,U10_B,U11}`):
    among them Teleman; Pandharipande–Pixton–Zvonkine; Givental; Graber–Pandharipande; Lho–Pandharipande;
    Zagier–Zinger; Arbarello–Cornalba; Guo–Janda–Ruan; Eynard–Orantin; Dunin-Barkowski et al.; Alexandrov et al.;
    pages of the DLMF.
    Scripts that compare against them need a local copy of the papers.
  - Five small Mathematica `.mx` files of the audit (`U09`, `U09_B`, `U14_A`), replaced by plain-text `.m` files
    with the same name, written with `Save` and checked to give identical definitions when read back.
  - `natural_form/genus5_cluster/local/out/parts_g5_89_nat_local.json` (23 MB) and the archive of the logs and
    partial results of the five genus-5 cluster runs (`g5_evidence.tar.gz`): see `../data/README.md`.
  - Snapshot tarballs and Python bytecode.

- **To reproduce the main outputs.** Each note of record names its drivers and logs; the genus-2 note
  `hodge_vertex_quintic/FORMAL_QUINTIC_genus2.md`, for example, lists `genus0_check.py`, `rows_check.py`,
  `loc_check.py`, `calib_check.py` and `genus2_frames.py`, run from `hodge_vertex_quintic/` with
  `python <script>`. Genus 3 and 4: `OU_EXTENSION_genus3.md` and `NATURAL_FORM.md` (s.5, the edge-splitting
  engine `fq4_split.py`). Genus 5: `cluster_genus5_quintic/README.txt`
  (`python genus5_run.py --validate`, then the SLURM arrays and `--collect`).
