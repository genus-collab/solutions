Genus-4 cluster package: Hodge-vertex graph sum of the resolved conifold x C^2 (fibrewise torus, mode B)
=====================================================================================================

What it computes. The 379 stable graphs / 27172 decorated graphs of genus 4 of the Givental-Teleman graph sum for F_{4,0} of
X x C^2, X = O(-1)+O(-1) -> P^1, on the fibrewise torus alpha0 = 0, alpha1 = eps1 = x, alpha2 = eps2 = 1 (opposite weights on each
O + O(-1) pair), with the local-P2 triple Hodge class as vertex (mode B: weights {eps1+eps2, -eps1, -eps2}, GW theory of C^3) and
R' = R_Z h_w^(-1) (hodge_vertex_conifold/README.txt in the project folder explains the frame; genus 2 and 3 were done locally).
Everything is exact: Q and x symbolic; for every decorated graph the term (rational function of Q and x) and its polar part at
Q = 1 in t = Log Q (coefficients of t^-6, ..., t^-1). Expected (from genus 2, 3 and from local P2 at genus 4): only the genus-4
vertex graph (97 decorated graphs) is singular at Q = 1 and carries exactly the gap polar part of P_4 Li_{-5}(Q)/(eps1 eps2),
(381 x^8 + 620 x^6 + 686 x^4 + 620 x^2 + 381)/(3870720 x) at t^-6 [= local P2 gap / (-3)^3]; singular decorated graphs only in the 32
stable graphs with a vertex of genus >= 2; pole order <= Sum_v max(2 g_v - 2, 0).

Files
  hvc_lib.m, hvc_run.wls, hvc_analysis.wls, hvc_summary.wls, hvc_rank4_check.wls      code (identical to hodge_vertex_conifold/)
  rc_lib.m, graphs.m, orbits.m, StableGraph-2/3/4.bin                                  shared conifold / graph code
  hvc_validate.wls      compares a fresh genus-2/3 run with the stored local runs ref/hvc_results_g{2,3}_Bm.m
  ref/                  resolvedconifold_results.m (F11, F20, F30 on the general torus), conifold_graphterms_genus2_3.m (point-vertex
                        terms 2026-09), hvc_results_g2_Bm.m, hvc_results_g3_Bm.m (local mode-B runs of 2026-10-30)
  env.sh                sourced by every job script and by run_workers.sh: module load Mathematica/13.2.1
  run_workers.sh        robust launcher (cluster_genus4_hodge/): one single-kernel wolframscript per core of an array task, staggered
                        starts, kernel-start timeout with process-tree kill, retries until the final "finished:" line, exit 1 on failure
  V_validate.sh         validation job (checks; genus 2 and 3 with 6 workers; comparison with ref/)
  C1_prep.sh            genus-4 factors (R' to z^9 with Q, x symbolic) + decorated graphs, one kernel, 64 GB
  C2_orbits.sh          array 1-20 x 6 workers = 120 chunks of the 27172 decorated graphs, 72 GB per task, restartable
  C3_collect.sh         collect, totals and gap, per stable graph, analysis, brief summary (one kernel, 64 GB)
  submit_genus4.sh      submits C1 -> C2 -> C3 with afterok dependencies

How to run (on Stanage, in ~/localP2genus4 or wherever the tarball was unpacked)
    sbatch V_validate.sh              # about 20-40 min on 6 cores
    cat validate_summary.log          # must end with "ALL VALIDATIONS OK -- ready for submit_genus4.sh" and contain no VALIDATION FAIL
    bash submit_genus4.sh             # C1 (hours), C2 (array), C3
Progress of C2: count the records written so far
    cat hvc_terms_g4_Bm_c*of120.m | grep -c "^{{"           (27172 when complete)
    grep -c "orbit {" hvc.log                               (every finished decorated graph; hvc.log is shared by all stages/genera)
If an array task dies or hits the time limit, simply resubmit C2 (sbatch C2_orbits.sh): the workers skip the decorated graphs already
in their chunk file (a torn last record makes the file unreadable; the worker then moves it aside and recomputes that chunk).
Then sbatch C3_collect.sh by hand (or with --dependency=afterok:<C2 job id>).

Results to copy back (INTO A DEDICATED EMPTY SUBFOLDER of the project, never into hodge_vertex_conifold/ itself):
    hvc_results_g4_Bm.m        per stable graph (term, polar part, counts) and per decorated graph (term, polar part); large
    hvc_table_g4_Bm.txt        per stable graph table
    hvc_analysis_g4_Bm.txt     structural analysis (regular/singular graphs, vertex-genus statement, pole-order bound, denominators)
    hvc_summary_g4_Bm.txt      polar parts of every decorated graph by stable graph (brief form, no terms)
    hvc.log, hvc_validate.log, validate_summary.log, C*.out, V_validate.*.out
The decorated-graph comparison with local P2 at genus 4 is part of the local P2 project (hodge_vertex/compare_conifold/ there,
since 2026-11-27), which holds a frozen copy of hvc_results_g4_Bm.m.

Reference at genus 4. ref/resolvedconifold_results.m has F_{g,0} only up to genus 3, so C3 compares the total with the closed form
F_4(0) + P_4 Li_{-5}(Q)/(eps1 eps2), F_g(0) = -(eps1^(2g-2) + eps2^(2g-2)) |B_2g B_2g-2|/(4g (2g-2) (2g-2)!) (the constant-map term;
reproduces the independently computed genus-2 and genus-3 values). This is a prediction, flagged as such in hvc.log; C3 also logs
separately whether total - P_4 Li_{-5}(Q)/(eps1 eps2) is Q-independent (true if and only if the Q-dependence is the refined conifold
formula) and the value of that constant.

SIZING (measured locally on the desktop, Mathematica 13.0, one kernel, 2026-10-30; Stanage was 1.5x slower per kernel for local P2):
  C1: R' factors to z^9 with Q and x symbolic 529 s (LeafCount of the largest edge factor 6049; file 1.7 MB); decorated graphs
      (379 graphs, 71198 labelled decorations, 27172 orbits) 132 s (file 11 MB). Expect about 20 min on the cluster.
  C2: 17 sample decorated graphs with one labelled decoration each (hodge_vertex_conifold/hvc_time_g4.txt): 0.0-0.3 s for graphs with
      <= 4 edges, 0.7-2.5 s with 5-7 edges, 4.5-4.8 s for the 9-edge six-vertex graphs (64 colourings), 0.6 s for the vertex decoration
      with nine dilaton leaves; the cost is proportional to the number of labelled decorations of the orbit (2.6 on average). Estimate
      for all 27172 decorated graphs: 20-100 CPU-hours, i.e. 10-50 min wall-clock on 120 workers (plus a minute per worker to load
      the factors and the orbit table). The 48 h limit and the restartable chunk files are precautions.
  C3: genus 3 collect + analysis took about 2 min locally for 693 decorated graphs; expect well under an hour at genus 4.
  Validation (V_validate.sh): genus 2 seconds, genus 3 about 5 min with 6 workers (2 min for 347 decorated graphs per chunk locally).
Memory: the genus-4 factor file hvc_factors_g4_Bm.mx (1.7 MB) and the orbit table (11 MB) are loaded by every worker; 72 GB per
6-worker task leaves 12 GB per kernel, far more than the exact route needs (the genus-3 workers stayed below 1 GB).
