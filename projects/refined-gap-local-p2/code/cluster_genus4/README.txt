Cluster computation of the twisted (generalised R-matrix) genus-4 free energy of local P2 x C2
==============================================================================================

What this computes
------------------
F_{g,0} of K_P2 x C^2 on the fibrewise torus T_F (weights alpha_i = (eps1+eps2)/3) with the identity row of the
R-matrix right-multiplied by ff[z] = Exp[q0 z + q1 z^3 + q2 z^5 + q3 z^7 + q4 z^9], q_k = Qrrc[0,k] of
localp2CY5_genRmat_v3.nb, kept symbolic. The Givental-Teleman graph sum (StableGraph-<g>.bin, automated decorations,
Witten-Kontsevich vertex numbers, conventions of localp2CY5_genus3.nb) is evaluated exactly at rational sample points
(t0 = y^(1/3), x0 = eps1/eps2) with PREC-digit arithmetic keeping C1i = 1/c1, C1p = c1' and q0..q4 symbolic, and
F_g is then reconstructed exactly as a polynomial in XX = y c1'/c1 and in the q's, with coefficients rational in y
and homogeneous in (eps1, eps2). Same method and same code as localp2CY5_genRmat_genus3.nb (whose results were
verified against localp2CY5_genRmat_v3.nb); only the sum over the 3^V colourings is done by variable elimination
instead of enumeration (identical results, checked at genus 2 and 3), and the work is split into independent jobs.

Files
-----
genRmat_g4_lib.m            library (sections 1-3 verbatim from the notebook, section 4 the cluster extensions)
genRmat_g4_sample.m         one sampling job:  wolframscript -file genRmat_g4_sample.m genus KMAX jobIndex nJobs [nKernels] [PREC]
genRmat_g4_collect.m        collector:         wolframscript -file genRmat_g4_collect.m genus
genRmat_g4_partA_pf.m       R-matrix, part A1: wolframscript -file genRmat_g4_partA_pf.m KMAX     (Picard-Fuchs ansatz, levels 0..KMAX)
genRmat_g4_partA_limit.m    R-matrix, part A2: wolframscript -file genRmat_g4_partA_limit.m KMAX  (QRR normalisation + delta -> 0 limit)
cPP_z6.m                    R-matrix identity column to z^6 (= cPPu6.mx of localp2final_CY5_v16.nb, portable text form)
PFu6.m                      Picard-Fuchs ansatz coefficients to level 6 (= PFu6.mx, portable): lets part A1 start at level 7
StableGraph-2/3/4.bin       stable graphs (binary, portable)
All files are plain text except the .bin graph files; nothing depends on the platform or on .mx files. Tested with
Mathematica 13.0 (wolframscript). Copy the whole directory to the cluster.

Step 1: R-matrix to z^9 (needed for genus 4; genus <= 3 runs with cPP_z6.m)
----------------------------------------------------------------------------
    wolframscript -file genRmat_g4_partA_pf.m 9        # writes PFu7.m, PFu8.m, PFu9.m (checkpoints per level; restartable)
    wolframscript -file genRmat_g4_partA_limit.m 9     # writes cPP_z9.m (and cPPu9_all.m with the intermediate data)
Single job, one core, a few GB of memory. Measured level times on the reference machine: level 2: 11 s, level 3: 65 s,
level 4: 154 s (levels 5, 6 previously took up to ~2 h in total); expect levels 7-9 to take hours up to about a day.
The logs partA_pf_num.log / partA_limit_u.log must say "cancellation of polar terms (all zero?): True".
Check: with KMAX = 4 the scripts reproduce cPP_z6.m[[1;;5]] exactly (done).

Step 2: sampling jobs (embarrassingly parallel)
-----------------------------------------------
Genus 4: 27 t-values x 20 eps-ratios = 540 sample points (pointList[4]); KMAX = 9; 379 graphs, 38,633 labelled
decorations, 33 q-monomials. Job j of N computes the points j, j+N, j+2N, ... and appends each finished point to
samples_g4_job<j>of<N>.m, so a job can be killed and restarted (it skips the points already in its file), and
different N can be mixed as long as every point is eventually covered.

    # SLURM array example, 60 jobs of ~9 points each, one core per job:
    #SBATCH --array=1-60
    wolframscript -file genRmat_g4_sample.m 4 9 $SLURM_ARRAY_TASK_ID 60

    # with local subkernels on a multi-core node (needs a licence for that many kernels):
    wolframscript -file genRmat_g4_sample.m 4 9 $SLURM_ARRAY_TASK_ID 60 8

Per-point cost (single core, PREC 350): genus 3 about 11-16 s; genus 4: see the line "genus 4 per point" at the end of
this file. Memory per job: a few GB. Each job logs one line per point in genRmat_g4_sample.log (shared, append-only);
"WARNING unexpected monomials" would indicate a precision problem (raise PREC).

Step 3: collect (locally or on the cluster)
-------------------------------------------
    wolframscript -file genRmat_g4_collect.m 4
Reads all samples_g4_job*.m, uses every t-value whose 20 eps-ratio points are present (lists what is missing), fits the
eps-structure (degree 2g+w palindromic polynomial, 2 spare points per fit), reconstructs the y-dependence (rational
reconstruction verified on >= 5 spare t-values) and writes
    genRmat_F40_cluster_result.m  = {Fg, Fgred, mons, FgQrrc, FgredQrrc}
    genRmat_F40_cluster_result.txt, genRmat_g4_collect.log
A partial set of jobs can be collected as soon as >= ~18 complete t-values exist (all 27 give the safety margin).
Then upload genRmat_F40_cluster_result.m (and the samples files, log) to the project folder and run, from the project
folder,  wolframscript -file genRmat_g4_analysis.wls 4  for the conifold-frame expansion and the gap conditions
(t_c^-5 .. t_c^-1, solved for Qrrc[0,3], Qrrc[0,4] with the genus-2/3 relations substituted, consistency check, and
the substitution observation of the genus-3 run).

Validation done before shipping
-------------------------------
- genus 2 through the whole pipeline (2 jobs, one with 2 subkernels, collector): exactly the twisted F20 of the notebook.
- genus 3 through the whole pipeline (1 job with 6 subkernels, collector): exactly the twisted F30 of the notebook.
- variable-elimination colouring sum == enumeration at a genus-3 sample point to 10^-291.
- Part A scripts at KMAX = 4 reproduce cPP_z6.m to order z^4.

genus 4 per point (measured on the cluster run of 2026-09-07, Mathematica 13.2.1, all 540 points, PREC 350)
------------------------------------------------------------------------------------------------------------
One kernel on one core: 31-34 min per point (mean 1919 s over 100 points), i.e. about 290 core-hours for the 540 points.
Eight subkernels on a 2-core allocation (cpus-per-task=2, kernel argument 8): mean 1.9 h per point (32 min to 2.6 h over
540 points). Oversubscribing the cores with subkernels is several times slower than one kernel per job: prefer more
single-kernel jobs (e.g. --array=1-180 with no kernel argument).
Wall-clock of the whole run: R-matrix levels 7, 8, 9 took 15, 25, 37 min and the delta -> 0 limit 3 min (single core);
sampling 14:51-19:43 (a 180-job single-kernel array, replaced after 100 points by a 50-job array with 8 kernels);
collection 1 min. Collector: 540 of 540 points, 27 complete t-values, 0 warnings, F40 LeafCount 35767.
