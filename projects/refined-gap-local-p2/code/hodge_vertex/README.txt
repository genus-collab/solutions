Hodge-vertex form of the genus-2 graph sum for K_P2 x C^2 on the fibrewise torus
===============================================================================

Idea (suggested by Guo and Lam). The Givental-Teleman reconstruction of the T_F-equivariant theory of Z = K_P2 x C^2,
    Omega_Z = R_Z . (point CohFT)^{+3}                       (localp2CY5_genus3.nb, cluster_genus4_polar/lp2_*),
is rewritten with the three idempotents carrying the GW theory of C^3 (Hodge CohFT  prod_a Lambda_g^vee(w_a)/w_a) with the
three fibre-direction weights  w = {-eps1-eps2, eps1, eps2}  (the normal weights of the fixed P^2, Calabi-Yau):
    Omega_Z = R' . (Hodge_w)^{+3},        R' = R_Z h_w^(-1),
    h_w(z) = Exp[ -Sum_{k>=1} B_{2k}/(2k(2k-1)) z^(2k-1) Sum_a w_a^(1-2k) ]      (scalar; h_w . point = Hodge_w by Mumford/QRR).
Both expansions give the same F_{2,0} (group-action property) but distribute it differently over the stable graphs.
Vertices: Hodge integrals  Int_{M_{g,n}} prod_a Lambda_g^vee(w_a)/w_a^g psi^a  (genus <= 2: Mumford's formula for lambda_1,
the lambda_g and lambda_g lambda_{g-1} closed forms, lambda_1^2 = 2 lambda_2), each psi degree <= dimension contributes.
Legs/edges/dilaton leaves: as before with R' in place of R_Z (rows of R' = rows of R_Z times the coefficient list of h^(-1)).

Modes computed
  point : Witten-Kontsevich vertices, R_Z          (the original expansion; 17 decorated graphs)
  A     : Hodge vertices, weights {-eps1-eps2, eps1, eps2}, R' = R_Z h^(-1)     (27 decorated graphs)
  B     : Hodge vertices, weights {eps1+eps2, -eps1, -eps2} (= paper/notebook convention e_T(K fibre) = +3 alpha,
          C^2 weights -eps1, -eps2), R' = R_Z h^(+1)  [h_{-w} = 1/h_w]; Lambda^vee(-w) = (-1)^g Lambda(w)
The sign s in R' = R_Z h^s is not a convention choice: with the geometric vertex for the stated weights exactly one sign
reproduces F_{2,0} (hv_calib.m; s = -1 for A), the other does not.
INTENDED CONVENTION (confirmed by the user 2026-10-29): mode B, the e_T weights of the normal bundle O(-3)+O+O of the fixed P^2
(twisted metric: e_T(K fibre) = 3 alpha - 3H = eps1 + eps2 - 3H, CY => C^2 weights -eps1, -eps2). Mode A is kept for comparison.

Main result at genus 2 (mode B): every one of the 20 decorated graphs with at least one edge is regular at the conifold
(hodge_vertex_B_polar_orbits.txt: polar parts {0, 0} orbit by orbit; the simple (1+27y)^-1 poles of the XX^0 and XX^1
coefficients cancel because XX = -2/9 at the conifold), and the 7 decorations of the single genus-2 vertex carry the whole polar
part, summing to the gap coefficient -(7 + 10 x^2 + 7 x^4)/(1920 x) at t_c^-2 and 0 at t_c^-1. The edge graphs are also much
simpler as functions (e.g. the (0,1) graph is (eps1^2 + eps2^2)(1/216 + 5 XX/144 + XX^2/16), y-independent).

Files
  hv_lib.m             Hodge integrals and checks, Hodge-vertex decorations, modified R-matrix rows, sampling, fits with
                       (eps1 eps2)^a (eps1+eps2)^b denominators   (loads lp2_lib.m -> graphs.m, orbits.m)
  hv_genus2.wls        driver:  checks | calib | sample <mode> [w nw] | recon
  lp2_lib.m, graphs.m, orbits.m, cPP_z6.m, StableGraph-2.bin, ref/F20_F30_y.m    copies from cluster_genus4_polar/
  hv_calib.m           the sign s and the weights
  hv_samples_<mode>_p<idx>.mx   one file per sample point (26 t x 34 eps-ratios), with the exact total-minus-F20 check
  hodge_vertex_genus2_results.m   Association: "point"/"A"/"B" -> <|graphTerms (7 stable graphs, order of StableGraph-2.bin),
                       orbitTerms (per decorated graph, keys {psi powers on half-edges, dilaton leaves}), polar {t_c^-2, t_c^-1},
                       total, totalOK, descriptions, weights|>, "graphs", "calibration", "hCoeffsA", "readme"
  hodge_vertex_genus2_summary.txt  human-readable per-graph terms and polar parts
  hodge_vertex_genus2.log          all checks

Checks performed (hodge_vertex_genus2.log)
  Hodge integrals: known values (1/24, 1/480, 7/5760, 1/2880, 1/5760); genus-1 Mumford = lambda_g formula; lambda_1^2 = 2 lambda_2
  on M_{2,n}, n <= 3, from Mumford's formula; string and dilaton equations of the Hodge correlators; the PPZ R-action of the
  scalar h on the point theory reproduces the normalised Hodge class on M_2 for arbitrary weights.
  Graph sums: sum over decorated graphs = F_{2,0} exactly at every sample point and symbolically after reconstruction, for all
  three modes; total conifold polar part = gap coefficient.

Run:  wolframscript -file hv_genus2.wls checks ; ... calib ; ... sample point ; ... sample A ; ... sample B ; ... recon

Genus 3 (mode B only) -- hv_lib3.m, hv_genus3.wls, and the runnable notebook localp2CY5_hodge_genus3.nb
----------------------------------------------------------------------------------------------------------
hv_lib3.m is the genus-general (g <= 3) version of hv_lib.m: Int lambda_1^k psi^a by k-fold Mumford (LpowRaw/Lpow), the
lambda_g formula (LgTop), lambda_g lambda_{g-1} (Lg1; closed form for a_i >= 1, zeros removed by the string equation),
lambda_1 lambda_2 lambda_3 on M_3 (Lg2), the normalised Hodge class expanded in lambda monomials by the Mumford relations
(hodgeClass, lamInt, hodgeNormU in the inverse weights u[a] = 1/w_a), the PPZ scalar-R-action check in any genus (ppzScalar), and
the reconstruction with a first-sample (a,b) search reused across t (xfitH with start).
Stages:  wolframscript -file hv_genus3.wls checks | prep | sample [w nw] | recon
  checks  known values; Mumford recursion vs closed forms (lambda_1^5 = 16 lambda_2 lambda_3, lambda_1^6 = 16 lambda_1 lambda_2 lambda_3
          with psi insertions, lambda_1^6 on M_3 = 1/90720); string/dilaton; PPZ in genus 2 (symbolic) and genus 3 (42 graphs, exact)
  prep    42 graphs, 1383 labelled decorations, 693 decorated graphs (orbits) -> hv_orbits_g3.mx            (1 s)
  sample  780 points (26 t x 30 eps-ratios), R' rows to z^6; about 3 s per point per worker (two workers: ~20 min);
          one file hv3_samples_B_p<idx>.mx per point with the exact total-minus-F30 check; bundle_samples.wls -> hv3_samples_B.mx
  recon   per-orbit reconstruction, 42 stable-graph terms, total = F30 (numeric + symbolic), gap, polar parts {t_c^-4..t_c^-1} per
          graph and per decorated graph -> hodge_vertex_genus3_results.m, hodge_vertex_genus3_summary.txt
RESULTS at genus 3 (2026-10-29; hodge_vertex_genus3_results.m / _summary.txt, hodge_vertex_genus3.log):
  * 693 decorated graphs sum to F_{3,0} exactly at all 780 points and symbolically; total polar part
    {-(1+x^2)(31+18x^2+31x^4)/(17920 x), 0, 0, 0} = the gap (leading coefficient (-3)^2 3! [t^6] t^2/(4 sinh(e1 t/2) sinh(e2 t/2))).
  * 41 of the 42 stable graphs are regular at the conifold; the single genus-3 vertex graph carries the whole polar part, i.e.
    the gap is a property of the genus-3 Hodge vertex dressed with the dilaton leaves of R'.
  * Decorated graphs: those whose vertices all have genus <= 1 (38 stable graphs) are regular one by one; the stable graphs
    with a genus >= 2 vertex (1 = the vertex, 2 = genus-2 vertex with a loop, 5 = (1)-(2), 6 = (0 with loop)-(2)) have singular
    decorated graphs (188 of 693 in all) which cancel only in the full sum over the decorations of the graph (grouping by psi
    assignment or by dilaton pattern alone is not enough). hv_group3.wls does this analysis.
Notebook: localp2CY5_hodge_genus3.nb (project folder; copy here) = the code files inlined as Input cells + driver cells
  (nb_cell_*.m; regenerate with make_hodge_notebook.wls): checks, genus 2 (recomputed in a few minutes and compared with the
  stored run), genus 3 (decorated graphs from hv_orbits_g3.mx, sample records from hv3_samples_B.mx -- delete it to resample --
  reconstruction, F30 and gap checks, per-graph polar table, singular decorated graphs). Data directory = the notebook's folder if
  it holds cPP_z6.m, else its subfolder hodge_vertex/. run_nb_headless.wls evaluates the notebook's Input cells with wolframscript
  (used as the test that the notebook runs).

Genus 4 (mode B) -- cluster_genus4_hodge/ on Stanage, results and analysis in this folder (2026-10-30)
---------------------------------------------------------------------------------------------------
Pipeline: cluster_genus4_hodge/ (hv_lib4.m: every lambda class as a polynomial in lambda_1 and ch_3, Int lambda_1^k1 ch_3^k3 psi^a
by Mumford's formula; stages prep/sample/reshape/recon/collect; see its README). Run: V_validate.sh (genus 2 + 3 reproduce the
stored runs here), then H1-H4: 680 sample points (34 t x 20 eps-ratios, 273 s per point on the cluster), 379 stable graphs,
71,198 labelled decorations = 27,172 decorated graphs, reconstruction 1 s per decorated graph, 0 failures. First attempt: 9 of
the 120 sampling kernels never started (LinkConnect::linkc) and hung their array tasks for 15 h; run_workers.sh is now robust
(staggered launches, start-up timeout, retries) and H2r_resume.sh completed the 54 missing points.
Files here: hodge_F40.m (total), hodge_stablegraphs_genus4.m (per stable graph: genera, edges, Aut, decoratedGraphs,
singularDecoratedGraphs, term, polar), hodge_graphterms_genus4.m (36 MB: every decorated graph with key, description, term and
polar part), hodge_collect_genus4.txt (table; rewritten by hv_table_g4.wls with canonical regularity), hv_analysis_g4.wls ->
hodge_vertex_genus4_summary.txt, hodge_vertex_genus4_analysis.m (Association), hodge_vertex_genus4_analysis.log;
stanage_run_2026-10-30/ = the complete copy of the cluster directory (logs, SLURM outputs, all .mx sample and chunk files, 4 GB).
RESULTS at genus 4:
  * The 379 stable-graph terms sum to F_{4,0} exactly (at all 680 points and symbolically); the total polar part is
    {-(381 + 620 x^2 + 686 x^4 + 620 x^6 + 381 x^8)/(143360 x), 0, 0, 0, 0, 0} = the gap ((-3)^3 5! [t^8] t^2/(4 sinh sinh)).
  * Exactly one stable graph is singular at the conifold: graph 1, the single genus-4 vertex (97 decorated graphs, 96 singular),
    and its polar part is exactly the gap (t_c^-6 only). All 378 other stable graphs are regular. (The collector's log line
    "regular: 348, singular: 31" tested unsimplified sums; after Together the 30 extra graphs are regular -- fixed in hv_collect.m.)
  * The 32 stable graphs with a vertex of genus >= 2 (28 with a genus-2 vertex, 3 with a genus-3 vertex, the vertex graph) are
    exactly the stable graphs containing singular decorated graphs; the 19,912 decorated graphs of the 347 graphs whose vertices
    all have genus <= 1 are regular one by one. Inside the 32 graphs 7,225 of 7,260 decorated graphs are singular; of the 35
    regular ones 26 are identically zero (vanishing Hodge integrals) and 9 (all in graph 11 = (1 with loop)-(2)) are nonzero.
  * Pole orders: every decorated graph has conifold pole order <= Sum_v max(2 g_v - 2, 0) over its vertices, i.e. each Hodge
    vertex of genus h >= 2 contributes at most the order 2h - 2 of the genus-h gap and the orders add (graph 7 = (2)-(2) reaches
    t_c^-4). The bound is attained in 31 of the 32 graphs; graph 118 (genera {0,0,1,2}, edges {{2,2},{1,2},{1,3},{1,4}}) has only
    simple poles. Odd orders do occur in individual decorated graphs (609 decorated graphs of genus-2-vertex graphs have a simple
    pole) and cancel in the sums.
  * Terms: all 379 stable-graph terms are polynomials in XX (degree <= 9) with denominators y^al (al <= 3) (1+27y)^be (be <= 6, so
    the polar lists t_c^-6..t_c^-1 are complete); 223 of the 378 regular graphs have no (1+27y) in the denominator at all, 155 have
    off-shell poles of their XX-coefficients that cancel on the curve XX = XX(y). Only 11 stable-graph terms have eps-denominators
    (eps1 eps2 at most, never eps1 + eps2).
Notebook: localp2CY5_hodge_genus3.nb section "Genus 4 (cluster results)" reads these files from the data directory
(ref/genRmat_F40_cluster_result.m is needed for the F40 comparison) and reproduces the gap, the regular/singular table, the
vertex-genus statement and the pole-order bound.
Copy-back pitfall (2026-10-30): the cluster directory was copied into this folder wholesale, overwriting README.txt and
hv_orbits_g3.mx (the cluster's .mx defines preparedH/OTH, not the notebook's preparedH3/OTH3, which silently broke the genus-3
section) -- both restored from hodge_vertex.tar.gz; copy cluster results into a dedicated subfolder.

Companion (2026-10-30): the same check for the resolved conifold x C^2 on the fibrewise torus with opposite weights on each
O + O(-1) pair, with the SAME triple Hodge class (mode B) as vertex: folder hodge_vertex_conifold/ of the resolved-conifold project
(README there), notebook resolvedconifold_hodge.nb. The decorated-graph-by-decorated-graph comparison of the two is compare_conifold/
in this folder (moved here 2026-11-27; frozen copies of the conifold results in compare_conifold/ref_conifold/). Result: the local P2 picture is reproduced exactly at genus 2 and 3 (only the vertex graph singular,
decorated graphs match one by one, universal ratios (-3)^(g-1) and (-3)^(g-1)(-(eps1+eps2)^2/(9(eps1-eps2)^2)) at the two leading orders).
