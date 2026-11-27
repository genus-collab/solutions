Hodge-vertex form of the Givental-Teleman graph sum for the resolved conifold x C^2 on the fibrewise torus
==========================================================================================================

Geometry and torus. X = O(-1) + O(-1) -> P^1, 5-fold Z = X x C^2. Conventions of resolvedconifold_genus3.nb / rc_lib.m: tangent
weights at the two fixed points {alpha0, alpha1, alpha2, -eps1, -eps2} and {-alpha0, alpha0+alpha1, alpha0+alpha2, -eps1, -eps2},
Calabi-Yau condition alpha0 + alpha1 + alpha2 = eps1 + eps2. Here: the FIBREWISE torus (alpha0 = 0, trivial action on the base P^1)
with OPPOSITE WEIGHTS ON EACH O + O(-1) PAIR, alpha1 = eps1, alpha2 = eps2, i.e. normal weights (eps1, eps2, -eps1, -eps2) along the
fixed P^1 (the specialisation used in the point-vertex comparison of 2026-09, conifold_graphterms_genus2_3.m). eps1 = x, eps2 = 1.

The frame (user's choice, 2026-10-30): the SAME triple Hodge class as for local P2, not the rank-4 class of the actual normal bundle.
    Omega_Z = R' . (Hodge_w)^{+2},     vertex = Int_{M_{g,n}} prod_{a=1}^3 Lambda_g^vee(w_a)/w_a^g prod psi^a   (GW theory of C^3),
    mode B:  w = {eps1+eps2, -eps1, -eps2}      mode A:  w = {-eps1-eps2, eps1, eps2}     (overall sign; Lambda^vee <-> Lambda),
    R' = R_Z h_w^s,  h_w(z) = Exp[-Sum_k B_2k/(2k(2k-1)) z^(2k-1) Sum_a w_a^(1-2k)],  s = -1 (the sign that reproduces F_{g,0};
    s = +1 does not, for either sign of the weights -- exactly as for local P2).
Mode B is the class that gave the local P2 results of hodge_vertex/ (in the code's eps conventions, where the C^2 weights are
-eps1, -eps2 for both theories; in the user's notation it is Lambda^vee(-eps1-eps2) Lambda^vee(eps1) Lambda^vee(eps2)). Mode A is
kept for comparison. Mode "point" = Witten-Kontsevich vertices with R_Z (the original expansion; regression against the 2026-09 files).
R_Z = conifold R-matrix of resolvedconifold_genus3.nb from the Wick expansion of the superpotential; at alpha0 = 0, alpha_i = eps_i
it needs no QRR normalisation constants: they are h_{(eps1, eps2, -eps1, -eps2)} = 1, because the rank-4 Hodge class of the normal
bundle O(-1)+O(-1)+O+O is trivial (Mumford's c(E) c(E^vee) = 1; checked: it integrates to the WK numbers for all 811 psi-tuples of
genus 1-3, hvc_rank4_check.wls). So with the rank-4 geometric class the Hodge frame would coincide with the point frame, where every
stable graph is singular at Q = 1 (genus 2: 7 of 7, genus 3: 42 of 42); the local-P2 triple class is a genuinely different frame.
Everything is exact: Q symbolic, x symbolic; terms are rational functions of Q and x; polar parts at Q = 1 in t = Log Q (coefficients
of t^-(2g-2), ..., t^-1), where F_{g,0} has its gap singularity P_g Li_{3-2g}(Q)/(eps1 eps2). Note: on this torus the two idempotents
collide at the discriminant Q_c = -4x/(1-x)^2 of U^2 - S U + P, not at Q = 1; at Q = 1 the roots are U = 1 (the logarithmic
singularity of the superpotential) and U = -1/x. EVERY stable-graph term (all modes) has a pole of order 4g-4 at Q_c (factor
(Q(1-x)^2 + 4x)^(4g-4) in the denominator; decorated graphs have orders between 0 and 4g-4), and these poles cancel only in the sum
over all stable graphs: the total F_{g,0} is regular at Q_c. [Corrected 2026-10-30 after the review: the first version of
hvc_analysis.wls tested divisibility with PolynomialRemainder, whose remainder with x-dependent coefficients does not simplify to 0,
and wrongly reported no discriminant poles; the detector now uses FactorList and checks that no unlisted factor remains.]

RESULTS
-------
Genus 2 (7 stable graphs; 27 decorated graphs with Hodge vertices, 17 with point vertices; all three modes sum to F_{2,0} exactly
and their total polar part is the gap (7 + 10 x^2 + 7 x^4)/(5760 x) at t^-2, 0 at t^-1):
  mode B : ONLY the genus-2 vertex graph is singular at Q = 1 and its polar part is EXACTLY the gap; all six graphs with edges are
           regular, and all 20 decorated graphs with an edge are regular individually (pole-order tally {0: 20, 2: 7}); the 7
           decorations of the vertex carry the polar part (leading coefficients not proportional to the gap polynomial).
           Decorated graph by decorated graph this is the local P2 mode-B picture: same 20 regular decorated graphs, and for the 7
           vertex decorations localP2/conifold = -3 at the leading order (t_c^2 = -3 t^2) and (eps1+eps2)^2/(3 (eps1-eps2)^2) at the
           next order -- the two universal factors of the 2026-09 point-vertex comparison (hvc_compare_localP2_g2.txt, now in the local P2 project, hodge_vertex/compare_conifold/).
  mode A : every stable graph singular (as mode A for local P2); 27 of 27 decorated graphs singular, all of order 2.
  point  : every stable graph singular; all 17 decorated graphs and polar parts agree with conifold_graphterms_genus2_3.m.
  Structure: mode-B stable-graph terms have denominators (1-Q)^2 (vertex graph only) Q (Q(1-x)^2+4x)^4 x^{<=1} (1+x)^2 (numerator and
  denominator of degree 7 in Q); every stable-graph term has a simple pole at Q = 0 (the canonical coordinates degenerate there,
  P = -1/(Q x)) and a pole of order 4 at the discriminant point Q_c, both cancelling only in the sum over all graphs (total at Q = 0 =
  F_{2,0}(0), total regular at Q_c); decorated graphs: discriminant orders {0: 1, 1: 2, 2: 7, 4: 17}.
Genus 3 (42 stable graphs; 693 decorated graphs with Hodge vertices, 378 with point vertices; all three modes sum to F_{3,0}
exactly (symbolically in Q and x) and their total polar part is the gap -(1+x^2)(31+18x^2+31x^4)/(161280 x) at t^-4, 0 below):
  mode B : ONLY the genus-3 vertex graph (30 decorated graphs) is singular at Q = 1 and its polar part is EXACTLY the gap; the 41
           other stable graphs are regular. 188 of the 693 decorated graphs are singular individually, and they sit exactly in the
           4 stable graphs with a vertex of genus >= 2 (graph 1 = (3): 29 of 30, graph 2 = (2) with a loop: 59 of 59, graph 5 =
           (1)-(2): 75 of 75, graph 6 = (0 with loop)-(2): 25 of 25); all 504 decorated graphs of the 38 graphs whose vertices have
           genus <= 1 are regular one by one. Pole order of every decorated graph <= Sum_v max(2 g_v - 2, 0), attained in all four
           graphs (tally {order, bound}: {0,0}: 504, {0,4}: 1, {2,2}: 159, {4,4}: 29); no decorated graph has an odd leading order.
           The single decorated graph with an identically vanishing term is the vertex with one dilaton leaf psi^2 (key {{}, {{2}}}),
           as for local P2. These are EXACTLY the local P2 mode-B numbers of hodge_vertex/ (41 of 42 regular; graphs 1, 2, 5, 6 with
           188 singular decorated graphs). Decorated graph by decorated graph (hvc_compare_localP2_g3.txt, now in the local P2 project, hodge_vertex/compare_conifold/): identical sets of 505
           regular decorated graphs; for the 29 (28) singular ones localP2/conifold = 9 = (-3)^2 at t^-4 and -(eps1+eps2)^2/(eps1-eps2)^2
           at t^-3, universal; at t^-2 (188 graphs) 33 distinct ratios, at t^-1 137 distinct ratios and 3 decorated graphs of graph 5
           with a nonzero coefficient on one side only -- the same "two universal leading orders" structure as in the point frame.
  mode A : every stable graph singular; 692 of 693 decorated graphs singular, all of order 4 (the vanishing term is again {{}, {{2}}}).
  point  : every stable graph singular; all 378 decorated graphs and polar parts agree with conifold_graphterms_genus2_3.m.
  Structure (mode B): stable-graph terms are rational in Q with denominators (1-Q)^4 (vertex graph only) Q^2 (Q(1-x)^2+4x)^8 x^{<=1}
  (1+x)^6, numerator and denominator degree <= 14 in Q; every stable graph has a double pole at Q = 0 and a pole of order 8 at the
  discriminant point Q_c = -4x/(1-x)^2 (same orders in modes A and point), all cancelling only in the sum over the 42 graphs (total
  regular at Q = 0 and at Q_c, equal to F_{3,0}(0) = -(eps1^4 + eps2^4)/1451520 at Q = 0 with eps2 = 1); decorated graphs: discriminant
  orders {0: 2, 2: 4, 3: 21, 5: 68, 6: 220, 8: 378} (the point frame: all 378 of order 8). The Laurent coefficient at t^-3 is nonzero
  for 28 vertex decorations although the stable-graph sum has none.
  Timings (Mathematica 13.0, desktop): R-matrix factors to z^6 with Q, x symbolic 65 s; 693 decorated graphs 4 min on two kernels;
  collect 45 s.
Notebook resolvedconifold_hodge.nb (project folder + copy here; 51 cells, 34 Input): core code inlined, checks, genus 2 recomputed
live for the three modes (seconds), genus 3 read from the stored hvc_results_g3_*.m (recompute3 = True recomputes, about 5 min per
mode on one kernel) and display cells. Evaluated headlessly with run_nb_headless.wls: every cell
produced the expected output, no messages (nb_headless.out, 2026-10-30). On 2026-11-27 the comparison cells with local P2 were removed
(they read files of the local P2 project; they are kept there, in the local P2 project, hodge_vertex/compare_conifold//nb_hvc_compare.m) and the notebook was rebuilt with
make_hvc_notebook.wls; headless rerun nb_headless_2026-11-27.out: all 34 Input cells evaluate, output identical to the 2026-10-30 run
for every cell kept, apart from the printed data directory.

FILES
-----
  hvc_lib.m            library: Get["rc_lib.m"] (conifold Frobenius data, quadratic algebra, Wick R-matrix, graph factors, graph
                       sum, polar parts) + sections 1-4 of hodge_vertex/hv_lib4.m verbatim (Hodge integrals by Mumford's formula
                       for lambda_1 and ch_3, normalised Hodge class, checks, decorated graphs with Hodge vertices) + conifold
                       specific: modeWeights, hodgeClassN (any number of weights), twistRows (rows of R_Z times the coefficient list
                       of h_w^s), buildFactors (R' at alpha0 = 0, alpha_i = eps_i, Q and x symbolic), loadOrbitsHVC, orbitTerm,
                       FrefHVC (F_{g,0} of ref/resolvedconifold_results.m specialised), gapHVC (polar part of P_g Li_{3-2g}/(e1 e2))
  hvc_run.wls          driver:  checks | g2sign | prep <g> <mode> <m|p> | orbits <g> <mode> <m|p> <chunk> <nchunks> | collect <g> <mode> <m|p>
                       (the sign is passed as the letter m or p: a numeric -1 would be swallowed by wolframscript as its own option)
  hvc_analysis.wls     structural analysis of one (genus, mode) run -> hvc_analysis_g<g>_<tag>.txt
  hvc_summary.wls      per stable graph / per decorated graph listing of terms and polar parts -> hvc_summary_g<g>_<tag>.txt
  (hvc_compare_localP2.wls, the decorated-graph comparison with local P2, and lp2_tc_section.m moved on 2026-11-27 to the local P2 project, hodge_vertex/compare_conifold/)
  hvc_rank4_check.wls  triviality of the rank-4 class at the level of integrals;  explore_pointframe.wls  the 2026-09 point-frame tables
  rc_lib.m, graphs.m, orbits.m, StableGraph-2/3.bin   copies from cluster_genus4_polar/ (identical to the hodge_vertex/ copies)
  ref/resolvedconifold_results.m   {F11, F20, F30} on the general torus;  ref/conifold_graphterms_genus2_3.m   point-vertex terms 2026-09;
  (ref/hv_lib4_localP2.m, a byte-identical copy of the local P2 library hv_lib4.m that the Hodge sections were taken from, removed 2026-11-27)
  hvc_factors_g<g>_<tag>.mx        graph factors (DumpSave, Mathematica 13.0);  hvc_orbits_g<g>_{hodge,point}.mx   decorated graphs
  hvc_terms_g<g>_<tag>_c<c>of<n>.m per decorated graph {id {graph, orbit}, key {psi on half-edges, dilaton leaves}, description,
                                   term (rational function of Q, x), polar {t^-(2g-2), ..., t^-1}, symmetric}
  hvc_results_g<g>_<tag>.m         Association: stable (per stable graph: genera, edges, Aut, decoratedGraphs, singularDecoratedGraphs,
                                   term, polar), entries, total, totalOK, polar, gapOK, readme
  hvc_table_g<g>_<tag>.txt         per stable graph table;  hvc.log  run log with VALIDATION lines
  make_hvc_notebook.wls + nb_hvc_*.m   generator of the notebook resolvedconifold_hodge.nb (project folder + copy here);
  run_nb_headless.wls  evaluates a notebook's Input cells with wolframscript (test)
Tags: Bm = mode B with s = -1, Am = mode A with s = -1, point.

Genus 4 (mode B, Stanage run of 2026-10-30: V_validate all OK, C1 prep 10 min, C2 array 18 min wall / 14 CPU-h, 120 chunk files copied back
into stanage_genus4/ and assembled locally; results stanage_genus4/hvc_results_g4_Bm.mx (full, DumpSave 219 MB) and .m (18 MB, without the
decorated-graph terms), hvc_table_g4_Bm.txt, hvc_analysis_g4_Bm.txt, hvc_summary_g4_Bm.txt (polar parts of all decorated graphs), and
hvc_compare_localP2_g4.txt, now in the local P2 project, hodge_vertex/compare_conifold/):
  the 27172 decorated graphs (379 stable graphs) sum EXACTLY (symbolic in Q, x) to the closed form F_4(0) + P_4 Li_{-5}(Q)/(eps1 eps2):
  total - P_4 Li_{-5}/(eps1 eps2) is Q-independent and equals -(eps1^6 + eps2^6)/87091200 = -|B_8 B_6|/(16 6 6!) (eps1^6+eps2^6), the
  predicted constant-map term; total polar part = gap (381 x^8 + 620 x^6 + 686 x^4 + 620 x^2 + 381)/(3870720 x) at t^-6, 0 below.
  ONLY the genus-4 vertex graph (97 decorated graphs, 96 singular) is singular at Q = 1 and its polar part is EXACTLY the gap; the 378
  other stable graphs are regular. 7234 decorated graphs are singular individually, exactly in the 32 stable graphs with a vertex of genus
  >= 2 (the same 32 graphs as for local P2); all 19912 decorated graphs of the other 347 graphs are regular; inside the 32 graphs the 26
  regular decorations are identically zero terms (graphs 1, 2, 6, 8). Pole order <= Sum_v max(2 g_v - 2, 0) for every decorated graph,
  attained in all 32 graphs (local P2: 31 of 32); tally {0: 19938, 2: 6069, 4: 1069, 6: 96}; NO decorated graph has an odd leading
  order (local P2 had 609 with t_c^-1). Denominators of the stable-graph terms: (1-Q)^6 (vertex graph only) Q^3 (Q(1-x)^2+4x)^12 x^{<=1}
  (1+x)^8, degree <= 21 in Q; the discriminant power 4g-4 = 12 in every stable graph and the Q = 0 poles (order 3) cancel only in the
  total, which is regular at Q = 0 (= F_4(0)) and at Q_c.
  Comparison with local P2 (hodge_graphterms_genus4.m, same 27172 keys): localP2/conifold = -27 = (-3)^3 at t^-6 for all 96 singular
  vertex decorations and 3(1+x)^2/(1-x)^2 = (-3)^3 (-(eps1+eps2)^2/(9(eps1-eps2)^2)) at t^-5 for all 94, universal; no universal ratio
  below (100, 482, 1268, 4886 distinct ratios at t^-4..t^-1), and the sets of regular decorated graphs differ slightly (conifold 19938,
  local P2 19947: the 609 local-P2 simple poles and the 9 regular nonzero decorations of graph 11 have no counterpart here), while the
  stable-graph-level statements coincide. Local assembly: collect 6.5 min, analysis + summary + comparison 10 min.
  PITFALL: writing the full results file with Put takes hours (1 MB/min for 27172 rational functions): collect now DumpSaves the full
  Association to .mx under the global symbol hvcResults and writes the slim .m; the chunk files keep the terms.
Mechanism and proof status (2026-10-30, user's question): see conifold_regularity_note.md in this folder (building blocks at Q = 1: the
singular idempotent is frozen in frame B, R~' -> diag(1, h_w^-2); only Delta_1 ~ t^-2 is singular; the vertex functional Z[h, n] with leaves
1/(c - psi) and dilaton leaves -psi^2/(c - psi) vanishes for n >= 1 (PROVED for all h >= 2 by the string equation, Theorem C of the note; three independent agent routes, six verifiers) and gives the gap for n = 0; scripts hvc_asymptotics*.wls,
hvc_normalised_R.wls, hvc_vertex_identity2.wls, hvc_percolouring.wls).
Genus 4 package: cluster package cluster_genus4_hodge_conifold/ in the project folder (tarball there; README inside): the same code
with restartable chunk files (orbits stage appends every finished decorated graph; hvc_validate.wls compares a fresh genus-2/3 run with
ref/), SLURM scripts V_validate.sh, C1_prep.sh, C2_orbits.sh (120 workers), C3_collect.sh, submit_genus4.sh. Local sizing
2026-10-30: factors to z^9 with Q, x symbolic 529 s, 27172 decorated graphs (71198 labelled decorations), 0-5 s per decorated graph
(hvc_time_g4.wls -> hvc_time_g4.txt), i.e. tens of CPU-hours in total. Copy the results back into a dedicated empty subfolder
(e.g. stanage_genus4/). The comparison with local P2 at genus 4 is run in the local P2 project, hodge_vertex/compare_conifold/.

Pipeline (local, Mathematica 13.0; genus 3 takes about 4 minutes per mode on two kernels):
    wolframscript -file hvc_run.wls checks ; wolframscript -file hvc_run.wls g2sign
    wolframscript -file hvc_run.wls prep 3 B m ; wolframscript -file hvc_run.wls orbits 3 B m 1 2 & ... orbits 3 B m 2 2 ; collect 3 B m
    wolframscript -file hvc_analysis.wls 3 Bm ; hvc_summary.wls 3 Bm

Mixed-scale vanishing and the all-genus programme (2026-10-30/31). Summary: conifold_regularity_note.md sections 7 and 7a.
  Findings and proofs: mixed_scale_findings.md (single leaf of another scale), proof_status_mixed_scale.md (synthesis of the proof workflow,
  with an addendum), proof_M_hurwitz.md, proof_M_residue.md, proof_M_framing.md (stub; that route ran out of output), referee_hurwitz_1.md,
  referee_residue_1.md, leaf_structure.md. Perl/BigRat checks of the workflow: review_checks/ (status_*.pl, referee1_*.pl).
  Scripts (wolframscript -file <script>; all exact in x):
    hvc_mixed_scale.wls       one leaf 1/(cp - psi) in the Theorem-C configuration (zeros cp = c/k, k = 1..2h-2)
    hvc_multiprobe.wls        Conjecture M on lattice points, genus 2-3, r <= 3, plus the lambda_h (x -> 0) limit   -> hvc_multiprobe_out.txt
    hvc_multiprobe_g4.wls     genus 4, r <= 2                                                                     -> hvc_multiprobe_g4_out.txt
    hvc_multiprobe_g4b.wls    genus 4, the irreducible points with r = 3, 4, 5 and D = 6 controls (185 s)         -> hvc_multiprobe_g4b_out.txt
    hvc_which_class.wls       which Hodge classes have the zeros j = 1..2h-2                                      -> hvc_which_class_out.txt
    hvc_check61.wls           closed one-leaf formula (6.1) against Xi_h(j), h = 2, 3, 4
    hvc_tr_tau.wls <tau> <gmax> <nmax> [hodge]   exact EO recursion on C_tau: x = v^tau/(1+v)^(tau+1), y = -v; Lemma D test and,
                              with "hodge", the dictionary P2 against the Hodge side (genus <= 4)  -> hvc_tr_tau<tau>_g<gmax>n<nmax>.txt
                              (gmax = 3, nmax = 3: about 10 s; genus 4 needs omega_{0,11} and is slow in this implementation)
    hvc_moment_test.wls       tau-dependence of the first non-vanishing moment (refutes the naive shear heuristic)
    hvc_leaf_structure.wls    leaf decomposition of the t-corrections (args: none, "book", "kt5")                -> hvc_leaf_structure_out*.txt
  hvc_mp_header.m is the shared header (Lval, msets, Xi) extracted from hvc_multiprobe.wls.

Leaf structure to all orders (Lemma C), proved without topological recursion (2026-10-31). Summary: conifold_regularity_note.md section 7b.
  The singular row is exactly hypergeometric: rho r_{0,1}(z; t) = (c/(c+z)) (tau/t) 2F1(1/z+1, 1/z+1; c/z+2; tau) Exp[-dphi(tau)/z], tau = 1 - e^-t
  (Euler integral, formal Laplace calculus, Gauss's equation, Stirling for the tau = 0 constant); the loops follow from the exact flow
  d/dt [t^2 rho^2 E_11] = Theta (d/dt t g)(d/dt t g'), Theta = -x^3 tau s*^2; the E_12 legs by rationality. Consequence: with Theorem R
  (proof_status_mixed_scale.md s.3) the all-genus regularity theorem is conditional only on Lemma A, M_irr(h) for h >= 5.
  leaf_proof/                 the proof: leaf_structure_proof.md (assembled reference text, 1890 lines), outline.md (plan and verified facts),
                              S1_thimble.md, S2_metric_qde.md, S3_leaves_loops.md, S4_legs_closed.md (section files, refereed "sound" and
                              revised), S*_checks*.pl with outputs (Perl/BigRat checks of the sections), referee_scripts/ (the referees'
                              independent Perl checks, the assembly spot checks and the critic's crit_rank.pl, with outputs)
  hvc_leaf_proof_check.wls [KT] [KX]   V1-V8 against the stored z^9 factors; closed form to t^KX   -> hvc_leaf_proof_check_out.txt (6 12: 112 s)
  hvc_pf_gauss.wls                     Picard-Fuchs = Gauss, and Delta_1 = x^2 W_2 mod the critical equation (prints)
  hvc_leaf_closed_forms.wls [KX]       closed forms A_j(t), e_il(t) against the closed form and the stored fits -> hvc_leaf_closed_forms_out.txt (12: 10 s)
  hvc_leaf_proof_check2.wls [KT] [KT2] requested checks against the stored factors: Q1, Q2, M8, M10b, M13-M15 -> hvc_leaf_proof_check2_out.txt (6 5: 63 s)
  hvc_leaf_proof_check3.wls            requested symbolic checks, no stored factors: M2-M7, M9-M11, Q3            -> hvc_leaf_proof_check3_out.txt (62 s)
  hvc_leaf_proof_check4.wls x0 [ZM] [KT]  M1: factors rebuilt from scratch to z^ZM at rational x0; KEY FORMULA, QDE row 2, symplecticity
                                       -> hvc_leaf_proof_check4_out_<x0>.txt (2 12 8 and 3/5 12 8: 40 s each)
  Logs leaf_proof_check*_run*.log are run logs only. All checks pass.
  Errata applied the same day: conifold_regularity_note.md F3 (E_22 at t = 0 is (x^3/c)[1 - h_w(z)^2 h_w(w)^2]/(z + w), not of rank one; E_22
  has no double pole in any frame) and s.5a (the Reduction Lemma's factorisation fails for E_11 edges between two singular vertices of
  genus >= 2; Theorem R is the valid reduction); hvc_vertex_identity.wls line 14 and hvc_asymptotics*.wls line 24 now read EdN[2,2] at U1.
  PITFALL (conjugate encoding): rc_lib.m stores every idempotent-2 quantity as cconj of the idempotent-1 one, so EVERY stored factor must be
  read at U1; read at U2, EdN[2,2] is the E_11 table. The graph sums are cconj-invariant (cSymQ) and unaffected.
