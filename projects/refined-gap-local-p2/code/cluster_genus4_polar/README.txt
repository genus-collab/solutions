Genus-4 cluster pipeline: refined local P2 x C2 and resolved conifold x C2, per decorated stable graph
==================================================================================================

What this computes
------------------
1. Local P2 x C2 on the fibrewise torus T_F (untwisted R-matrix from cPP_z9.m): the contribution of every decorated stable
   graph (orbit of psi/dilaton decorations) of genus 4 to F_{4,0}, exactly in y, XX = y c1'/c1, eps1, eps2, by sampling and
   reconstruction (same conventions as localp2CY5_genus3.nb and cluster_genus4/), and the polar part of each in the conifold
   coordinate t_c. The sum is checked against the q -> 0 part of the twisted cluster result genRmat_F40_cluster_result.m
   (an independent R-matrix normalisation path) and the gap; GV invariants, refined BPS numbers, the anomaly equation and the
   refined split follow in lp2_checks.m.
2. Resolved conifold x C2 on the general five-dimensional torus: F11, F20, F30, F40 as exact rational functions of Q and the
   weights (sampling at integer weights, exact Q-reconstruction, symmetric-polynomial interpolation; degree-zero term symbolic),
   with the CY-restriction (rigidity) and degree-1 checks.
3. Resolved conifold at alpha0 = 0, alpha_i = eps_i: every decorated graph of genus 4 as a rational function of Q and
   x = eps1/eps2 and its polar part at Q = 1 (t = Log Q).
4. Resolved conifold on the CY torus alpha0 = eps1+eps2-alpha1-alpha2 with alpha1, alpha2 free: the polar part of every
   decorated graph of genus 4 as exact rational functions of alpha1, alpha2, eps1, eps2 (Laurent-series sampling at integer
   weights with the closed-form QRR normalisation, multivariate fit verified on spare points).
5. The decorated-graph-by-decorated-graph comparison of the polar parts (local P2 versus conifold).

Files
-----
graphs.m, orbits.m                  stable graphs, decorations, Witten-Kontsevich numbers; decorated graphs = Aut-orbits (shared)
lp2_lib.m                           local P2 library (core routines of genRmat_g4_lib.m, untwisted; per-orbit sampling; fits; t_c)
lp2_prep.m / lp2_sample.m / lp2_reshape.m / lp2_recon.m / lp2_collect.m / lp2_checks.m     local P2 stages
rc_lib.m                            conifold library (resolvedconifold_genus3.nb cells + closed-form QRR constants, alpha0 = 0
                                    R-matrix, Laurent-series layer, CY-torus fit)
rc_prep.m / rc_prep_x.m / rc_sample.m / rc_const.m / rc_collect.m                            general torus and alpha0 = 0 factors
rc_polar_x.m / rc_polar_collect.m                                                            per-graph polar parts at alpha0 = 0
rc_polarCY_sample.m / rc_polarCY_reshape.m / rc_polarCY_fit.m / rc_polarCY_collect.m        CY-torus per-graph polar parts
rc_check_ffc.m                      closed-form QRR constants versus the Q -> 0 limit (part of the validation)
make_ref_y.wls                      one-off local utility that produced ref/F20_F30_y.m from the project's F20/F30_result.mx (already run)
compare_g4.m                        the comparison
cPP_z6.m, cPP_z9.m                  local P2 R-matrix identity column (from cluster_genus4/)
StableGraph-2/3/4.bin               stable graphs
ref/                                stored genus-2/3 results used by the validation stage, and genRmat_F40_cluster_result.m
env.sh                              sourced by every job script and by run_workers.sh: initialises the module system if needed and
                                    runs "module load Mathematica/13.2.1" (nothing has to be loaded by hand before sbatch)
run_workers.sh                      launches one single-kernel wolframscript per core of an array task (6 workers x 20 tasks = 120)
V_validate.sh                       validation job (genus 2 and 3, all stages, 6 workers)
S1_prep.sh ... S10_compare.sh       genus-4 stages;  submit_genus4.sh submits them with dependencies
All text files are ASCII; intermediate data on the cluster are DumpSave .mx files (one Mathematica version), final results .m.
Tested locally with Mathematica 13.0 (wolframscript) at genus 2 (local P2: full stack; conifold: see validation); the cluster
validation stage repeats genus 2 and 3 for every stage with 13.2.1.

How to run
----------
Copy the whole directory to the cluster, then in it:
    sbatch V_validate.sh                     # about 1-2 h on 6 cores; wait for it
    cat validate_summary.log                 # must show only "VALIDATION OK" lines ("ALL VALIDATIONS OK -- ready for submit_genus4.sh")
    bash submit_genus4.sh                    # submits S1..S10 with afterok dependencies (see below)
Results to copy back to the project folder: localP2_F40.m, localP2_graphterms_genus4.m, localP2_polar_genus4.m,
lp2_checks_g4.log, resolvedconifold_F40.m, rc_collect.log, conifold_graphterms_genus4.m (+ _full.m), conifold_CY_graphterms_polar_genus4.m,
polar_check_genus4.log, polar_check_genus4_table.txt, and the *.log files.

Stages (genus 4)
----------------
S1  prep (1 core):           lp2_prep.m 4 (orbits), rc_prep.m 4, rc_prep_x.m 4 9 (alpha0 = 0 factors with Q, x symbolic to z^9),
                             rc_const.m 4 9 (degree-zero term, symbolic in the five weights)
S2  lp2 sampling (20 x 6):   lp2_sample.m 4 9 w 120 -- 714 points (34 t x 21 eps-ratios), per point the XX^0..XX^9 coefficient
                             vectors of all decorated graphs (PREC 350, Rationalize). Expected about 7 min per point.
S3a lp2 reshape (1 core, 64 GB): regroup by decorated graph into 120 chunks;   S3b lp2 reconstruction (20 x 6): exact terms and
                             t_c polar parts (NC = 12)
S4  lp2 collect + checks (1 core): sum == q -> 0 twisted F40, gap, per-graph files; lp2_checks.m 4 (GV, BPS, anomaly, split)
S5  rc sampling (20 x 6):    rc_sample.m 4 9 w 120 closed -- 187 integer weight samples, per sample the R-matrix to z^9, 18
                             Q-values, graph sums for genus 1 (one leg), 2, 3, 4, Q-reconstruction. Cost unmeasured at genus 4
                             (genus 3: 75 s per sample including all genera); the log line per sample reports the time per genus.
S6  rc collect (1 core):     assemble genus 1..4 (genus 1-3 must reproduce ref/resolvedconifold_results.m), CY and degree-1 checks
S7a rc per-graph at alpha0 = 0 (20 x 6):  rc_polar_x.m 4 c 120 (exact rational functions of Q, x + polar parts); S7b collector
S8  rc CY-torus sampling (20 x 6):  rc_polarCY_sample.m 4 9 w 120 -- 460 integer points, per point the polar parts of all
                             decorated graphs by Laurent-series arithmetic to absolute order t^(6g-2) = t^22 (a product with E
                             edge factors keeps its coefficients up to t^(M-2E); the vertex factors Delta^(g_v-1) are expanded from
                             exact rational functions, since series inversion of Delta loses precision). Measured on the cluster:
                             7-9 s per point at genus 3 with order t^12; expect roughly 10-20 min per point at genus 4.
S9a reshape (1 core, 64 GB); S9b fit (20 x 6): rc_polarCY_fit.m 4 c 120 (ansatz per order: denominator eps1^6 eps2^6 (eps1+eps2)^(2(n-1)),
                             numerator degree 18 + 2(n-1), alpha-degree <= 2(n-1), escalated if a coefficient fails to verify);
                             S9c collector (structure checks, totals, output file)
S10 comparison (1 core):     compare_g4.m 4 -> polar_check_genus4.log, polar_check_genus4_table.txt
Dependencies: S1 -> S2 -> S3a -> S3b -> S4;  S1 -> S5 -> S6;  S1 -> S7a -> S7b (also after S6);  S1 -> S8 -> S9a -> S9b -> S9c;
S4, S7b, S9c -> S10.  Expected total: about 1000-1500 core-hours, 12-20 h wall-clock (S5 and S8 dominate).

Restartability and sizing
-------------------------
Every worker writes each finished item to its own output file (the samplers one small .mx file per point, the other stages one
record per item appended to a text file) and skips items already present, so a killed array can simply be resubmitted; different worker counts can be mixed as long as every item is eventually covered (items are assigned round-robin
by index). 120 concurrent kernels need 120 licence tokens (180 were used in the 2026-09 run). Memory: 48 GB per 6-worker task,
64 GB for the single-core reshape/collect jobs. If S5 or S8 turn out slower than estimated, the arrays can be resubmitted after
raising --array (and NTASKS) -- the item assignment changes, but finished items are still skipped by index.

Validation criteria (V_validate.sh, genus 2 and 3)
--------------------------------------------------
- local P2: sum of the decorated graphs == F20, F30 of ref/F20_F30_y.m (y, XX, eps form); conifold gap; genus-3 per-graph terms and polar parts ==
  ref/localP2_graphterms_genus3.m; lp2_checks_g3.log: anomaly equation, GV table rows g <= 3, gap.
- conifold: degree-zero terms, F11/F20/F30 == ref/resolvedconifold_results.m; CY restriction; per-graph polar parts and rational
  functions at alpha0 = 0 == ref/conifold_graphterms_genus2_3.m; CY-torus polar parts == ref/conifold_CY_graphterms_polar.m;
  closed-form QRR constants == Q -> 0 limit to z^9 at three general weight points.
