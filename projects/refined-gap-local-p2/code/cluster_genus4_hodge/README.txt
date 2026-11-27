Genus-4 Hodge-vertex pipeline: K_P2 x C^2 on the fibrewise torus, F_{4,0} per decorated stable graph with Hodge vertices
=====================================================================================================================

What this computes
------------------
The Givental-Teleman graph sum for F_{g,0}(K_P2 x C^2, T_F) written with three copies of the Hodge CohFT (equivariant GW theory
of C^3, vertex classes prod_a Lambda_g^vee(w_a)/w_a) carrying the e_T weights w = {eps1+eps2, -eps1, -eps2} of the normal bundle
O(-3) + O + O of the fixed P^2, and the R-matrix R' = R_Z h_w^(-1), h_w(z) = Exp[-Sum_k B_2k/(2k(2k-1)) z^(2k-1) Sum_a w_a^(1-2k)]
(R_Z = the R-matrix of localp2CY5_genus3.nb / cluster_genus4_polar, identity column cPP_z9.m). The total equals the point-vertex
F_{g,0}; the pipeline delivers the contribution of every stable graph and every decorated graph (Aut-orbit of psi/dilaton
decorations) exactly in y, XX = y c1'/c1, eps1, eps2, with the polar parts at the conifold (t_c^-(2g-2), ..., t_c^-1).
Genus 2 and 3 (computed locally, hodge_vertex/ in the project folder): all stable graphs except the single genus-g vertex graph
are regular at the conifold, and that graph carries the whole polar part (= the gap). Genus 4 tests this with 379 stable graphs.

Hodge integrals (hv_lib4.m, valid in genus <= 4 without closed forms): every lambda class is a polynomial in lambda_1 and ch_3(E)
(Newton with ch_2 = ch_4 = 0: lambda_2 = lambda_1^2/2, lambda_3 = lambda_1^3/6 + 2 ch_3, lambda_4 = lambda_1^4/24 + 2 lambda_1 ch_3);
Int lambda_1^k1 ch_3^k3 psi^a by Mumford's formula for ch_1 and ch_3 (kappa_k = pi_* psi^(k+1), boundary restrictions additive,
subset sums organised by sub-multisets, memoised). Checks (hv_checks.m): lambda_1 on M_{1,1}; lambda_1 psi^3, lambda_2 psi^2,
lambda_1 lambda_2 psi, lambda_1 lambda_2, ch_3 = -1/34560 on M_2; lambda_1^6 = 1/90720 on M_3; lambda_3 psi^4, lambda_3 lambda_2 psi^2,
lambda_1 lambda_2 lambda_3 (genus 3); lambda_1^9 = 1/113400 on M_4, lambda_4 psi^6 = 127/154828800, lambda_4 lambda_3 psi^3 = 1/3225600,
lambda_4 lambda_3 psi^2 psi^2 = 1/276480 (Getzler-Pandharipande), lambda_4 lambda_3 lambda_2 = 1/87091200; vanishing above the top
degree of the lambda ring computed without the cut-off; lambda_g and lambda_g lambda_{g-1} closed forms versus the recursion;
string and dilaton equations; PPZ scalar R-action of h on the point theory = normalised Hodge class on M_g (genus 2 symbolic,
genus 3 and 4 at exact rational weights, genus 4 over the 379 graphs).

Files
-----
graphs.m, orbits.m, lp2_lib.m       shared with cluster_genus4_polar/ (stable graphs, WK numbers, orbits; R-matrix rows, colour sums, fits)
hv_lib4.m                           Hodge integrals, Hodge-vertex decorations, R' rows, sampling, reconstruction with
                                    (eps1 eps2)^a (eps1+eps2)^b denominators (a, b <= 1 at genus 2, 3; AMAX = BMAX = 5)
hv_prep.m / hv_sample.m / hv_reshape.m / hv_recon.m / hv_collect.m / hv_checks.m      the stages (usage in each file's header)
cPP_z6.m, cPP_z9.m                  R-matrix identity column to z^6 / z^9;  StableGraph-2/3/4.bin  stable graphs
ref/F20_F30_y.m                     F20, F30 (y, XX, eps);  ref/genRmat_F40_cluster_result.m  F40 (q -> 0 part = untwisted F40)
ref/hodge_vertex_genus2_results.m, ref/hodge_vertex_genus3_results.m    stored local runs (per decorated graph) for the validation
env.sh, run_workers.sh              Mathematica/13.2.1 module; 6 single-kernel workers per array task (20 tasks = 120 workers)
V_validate.sh                       genus 2 and 3, all stages, versus the stored runs; then hv_checks.m 4
H1_prep.sh ... H4_collect.sh        genus-4 stages;  submit_genus4.sh submits them with afterok dependencies

How to run (on the cluster, in this directory)
----------------------------------------------
    sbatch V_validate.sh                     # genus 2 + 3 end to end (about 1 h on 6 cores) and the Hodge checks
    cat validate_summary.log                 # must show only "VALIDATION OK" lines ("ALL VALIDATIONS OK -- ready for submit_genus4.sh")
    bash submit_genus4.sh                    # H1 prep -> H2 sampling (20 x 6 workers) -> H3a reshape -> H3b reconstruction -> H4 collect
Results to copy back: hodge_F40.m (total), hodge_graphterms_genus4.m (every decorated graph with its polar part),
hodge_stablegraphs_genus4.m (per stable graph: term, polar part, number of singular decorated graphs), hodge_collect_genus4.txt
(table), hv_collect.log, hv_sample.log, hv_prep.log, hv_reshape.log, hv_recon.log.

Stages (genus 4)
----------------
H1  prep (1 core):           hv_prep.m 4 -- decorated graphs with Hodge vertices: 379 graphs, 71,198 labelled decorations, 27,172
                             decorated graphs (5,647 dilaton-free), 2 min (the point version has 38,633 labelled decorations);
                             hv_orbits_g4.mx is already included in the package, so H1 only re-creates it
H2  sampling (20 x 6):       hv_sample.m 4 9 w 120 -- 680 points (34 t x 20 eps-ratios), per point the XX^0..XX^9 coefficient
                             vectors of all decorated graphs (PREC 350, exact rationals), each point checked against F40 (q -> 0 of
                             the twisted cluster result). Measured locally (hv_timing4.wls, Mathematica 13.0, one kernel): 182 s per
                             point (R-matrix rows to z^9: 2 s; evalPoint 0.1 s; the rest is the 27,172 graph sums), and the orbit sum
                             equals F40 at the point. 680 points on 120 workers: about 6 points per worker, well under an hour of wall
                             time even at half the local speed. Restartable, one file per point.
H3a reshape (1 core, 64 GB): regroup by decorated graph into 120 chunks (VALIDATION line: orbit sums = F40 at every point)
H3b reconstruction (20 x 6): hv_recon.m 4 c 120 -- exact terms and t_c polar parts (tcSetup[12])
H4  collect (1 core):        total = F40, gap, per-graph files, table of regular / singular stable graphs and decorated graphs
Each stage is restartable (a resubmitted array skips the items already in the output files).

Local test: wolframscript -file hv_prep.m 2 ; hv_sample.m 2 3 1 1 ; hv_reshape.m 2 4 ; hv_recon.m 2 c 4 (c = 1..4) ; hv_collect.m 2
(about 5 min), and the same with genus 3 (KMAX 6, about 25 min on two workers: hv_sample.m 3 6 w 2).

Run of 2026-10-29 on Stanage and the robust launcher (2026-10-30)
-----------------------------------------------------------------
V_validate.sh passed. H2 (20 x 6 workers) sampled 626 points in 33 minutes (18:04 - 18:37; 273 s per point on average, at most
330 s, i.e. 1.5 x the desktop), but 9 of the 120 kernels never started: their worker_hv_sample_w*.out hold only the harmless
"Failed to open configuaration file" lines and, for four of them, "LinkConnect::linkc" -- the WSTP link to the kernel was never
made when 120 kernels were launched within 5 s -- and the wolframscript processes hung for more than 15 h, keeping 8 array
tasks alive (run_workers.sh waits for its 6 workers) and holding the afterok chain. The 54 points of those 9 workers
(w = 8, 11, 29, 39, 54, 56, 73, 99, 104: the points w, w + 120, ...) are the only missing ones.
run_workers.sh is now robust: launches STAGGER = 5 s apart; a worker with no kernel output within START_TIMEOUT = 600 s is
killed (process tree) and relaunched; a worker that exits without its final line ("job finished" / "finished: ...") is
relaunched (the stages skip what is on disk); MAX_TRIES = 4, then the task exits 1 and the chain stops. Existing worker
*.out files are moved to old_worker_out/. Tested locally with a stub wolframscript (hang, early exit, persistent failure).
    bash submit_resume.sh        # scancel the stalled H2/H3a/H3b/H4 by name; H2r_resume (9 workers, missing points only,
                                 # about 30 min) -> H3a -> H3b (robust launcher) -> H4
Counting progress: hv_sample.log is shared by all genera (the validation wrote 1664 point lines before genus 4), so count the
genus-4 point files:  ls hv_samples_g4_p*.mx | wc -l   (680 when complete).
Fix 2026-10-30: hv_collect.m's regQ now canonicalises (Together) before testing a summed polar part for zero; the genus-4 run's
log line "regular at the conifold: 348; singular: 31" was an artefact of this (after Together: 378 regular, only the genus-4
vertex graph singular). Results and analysis: hodge_vertex/ (README there), hv_analysis_g4.wls, hodge_vertex_genus4_summary.txt.
