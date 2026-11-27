# RESUME: package the second-prime run in the frame (1/3, -2/7) for Stanage

Written 2026-11-13, before a machine restart. Nothing was in flight: no job, no agent and no background process was running when this was
written, and no file of the package was modified for this task yet. The genus-5 note (`../../NATURAL_FORM.md` s.6) and the memory are
already up to date with the audit of the copied-back evidence; what follows is only the packaging of the one remaining cluster job.

## What the user asked for

The second-prime natural run in the second conifold frame: tag **`g5_61_nat_ab1o3_m2o7`**, full `--maxout`, so that the second frame's 628
rows rest on two primes instead of one (s.6, "What is not established at genus 5", first bullet). About 52 CPU hours; 100 one-core array
tasks; the heaviest single inner graph sets the wall time at about 0.85 h on Stanage (`f = 0.94`, measured).

## What is already settled (do not re-derive)

* **No code change is needed.** `genus5_run.py` takes `--prime 61 --variant nat --ab 1/3,-2/7`; `g5_checks.tag_of` then gives the tag above.
* **Additivity will run this time** if the point-frame parts of the same prime *and* frame are present. `genus5_run.py` builds
  `idtag = tag_of(genus, prime, 'id', ab, None)` = `g5_61_id_ab1o3_m2o7` and runs the check when `len(idp) == 4555` and `--maxout` is unset;
  otherwise it logs `(no complete ... run: the additivity check is skipped)` and `failed` stays `[]`. That is why it was skipped in the
  p89 frame-2 run.
* **Those parts exist and are audited**: `frame2/out/parts/g5_61_id_ab1o3_m2o7/local.json` (963 kB, `tag` field correct, 4555 inner graphs)
  and the p89 one beside it (1.27 MB). They were made and checked on 2026-11-13 (`FRAME2.md`): 4555 graphs, GAP pass, `c_5 = 1/660000`, and
  all 4555 additivity identities against the cluster's stored pieces.
* **The SLURM plumbing is adequate as it stands**: `S_array.sh` has `--time=48:00:00`, `--mem=2G`, one core; `C_collect.sh` 12 h and 8 GB;
  `env.sh` loads any Python >= 3.8 (the cluster used 3.9.25). Observed per-task wall time in the p89 frame-2 run: 1756 / 1900 / 3025 s.
* The package lives in `~/formalquintic/cluster_genus5_quintic` on Stanage, and the cluster still holds its own `parts/` there.

## The plan that was about to run (4 agents, ultracode)

1. **verify** - copy the package to the scratchpad (never smoke-test inside the real package: a stray `task<I>of<N>.json` would be merged by
   `read_parts`, whose regex accepts any splitting) and there (a) run `V_validate.sh`/`--validate`, (b) run `--genus 4 --prime 61 --variant
   nat --ab 1/3,-2/7 --run` (4.5 min) and compare with the genus-4 frame-2 p61 table of record in `../../g4_split_*.json`, (c) smoke-test the
   genus-5 code path itself with one cheap graph (`--task` of a large `--ntasks`, not task 0: LPT puts the heaviest graph there), (d) confirm
   the tag and that `read_parts` finds the shipped id parts.
2. **plan-commands** - the exact copy-paste sequence (scp, extract, submit, `tail -n`), the wall/CPU estimate from the measured p89 frame-2
   run, and a list of every path the supplementary tarball contains, to show it cannot clobber the cluster's own results.
3. **package** - new files only, nothing existing modified: `submit_frame2_p61.sh` (wrapping `sbatch --array=0-99 S_array.sh 100 --prime 61
   --variant nat --ab 1/3,-2/7` plus the `afterany` collect job), `README_FRAME2_P61.txt`, the two id part directories copied into
   `parts/`, and a new standalone `g5_prov.py` implementing s.6's own "for the next run" recommendation (record every check as
   pass/fail/**skipped**, assert `min(valO) == 1` *together with* a non-zero entry of `E_out`, and stamp host, interpreter, CPU seconds,
   part-file inventory and a digest of the `.py` files into `results/<tag>.prov.json`) - `g5_checks.py` and `genus5_run.py` must not be
   touched.
4. **critic** - adversarial pass: a frame- or prime-dependent input that is not shipped, the `assert part['tag'] == tag` on the shipped
   parts, tag collisions with the five completed runs, the resumability claim, absolute paths, and that the user's email appears only in the
   established `#SBATCH --mail-user` lines.

Optional extra, flagged but **not** part of what the user approved: a `nat_t6 --maxout 2` control (about 15 CPU h), which would be the first
test of `tau_6` at any genus. Mention it; do not bundle it into the default submit target.

## Deliverables

`cluster_genus5_quintic.tar.gz` rebuilt in the main folder, plus a small overlay tarball with only the new files (safer to extract over the
cluster copy), the copy-paste commands as one fenced `bash` block each, and a short entry in s.6a of the note once the job has run.
