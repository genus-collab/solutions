(* genus 4 is computed on the cluster with cluster_genus4_hodge/ (379 stable graphs; the same code: hv_lib4.m generalises the Hodge
   integrals to genus 4 with Int lambda_1^k1 ch_3^k3 psi^a by Mumford's formula). This section reads the results copied back
   into the data directory: hodge_stablegraphs_genus4.m (per stable graph), hodge_graphterms_genus4.m (per decorated graph). *)
g4files = {"hodge_stablegraphs_genus4.m", "hodge_graphterms_genus4.m"};
If[! AllTrue[g4files, FileExistsQ], Print["genus-4 result files not present in ", Directory[], " -- run cluster_genus4_hodge/ on the cluster and copy back ", g4files], stable4 = Get["hodge_stablegraphs_genus4.m"]; entries4 = Rest[Get["hodge_graphterms_genus4.m"]]; Print[Length[stable4], " stable graphs, ", Length[entries4], " decorated graphs"]]
(* total polar part {t_c^-6, ..., t_c^-1} versus the gap coefficient (-3)^3 5! [t^8] t^2/(4 sinh(e1 t/2) sinh(e2 t/2)) *)
If[ValueQ[stable4], {Factor /@ Total[stable4[[All, "polar"]]], Factor[gapCoefficient[4] /. {e1 -> x, e2 -> 1}]}]
(* which stable graphs are regular at the conifold, and which contain singular decorated graphs *)
If[ValueQ[stable4], regular4 = Select[Range[Length[stable4]], Union[Together /@ stable4[[#, "polar"]]] === {0} &];   (* canonicalise the summed polar parts before testing for zero *)
  {Length[regular4], "regular of", Length[stable4], "; singular stable graphs:", Table[{k, stable4[[k, "genera"]], stable4[[k, "edges"]], Factor /@ stable4[[k, "polar"]]}, {k, Complement[Range[Length[stable4]], regular4]}]}]
If[ValueQ[stable4], Grid[Prepend[Table[{k, stable4[[k, "genera"]], stable4[[k, "edges"]], stable4[[k, "decoratedGraphs"]], stable4[[k, "singularDecoratedGraphs"]]}, {k, Select[Range[Length[stable4]], stable4[[#, "singularDecoratedGraphs"]] > 0 &]}], {"graph", "genera", "edges", "decorated graphs", "singular decorated graphs"}], Frame -> All]]
(* the stable graphs with a vertex of genus >= 2 are exactly those containing singular decorated graphs, and the conifold pole
   order of every decorated graph is bounded by Sum_v max(2 g_v - 2, 0) over its vertices: {sets equal, their number, bound holds,
   tally of {pole order, bound}} *)
If[ValueQ[stable4], Module[{ge2 = Select[Range[Length[stable4]], Max[stable4[[#, "genera"]]] >= 2 &], sg = Select[Range[Length[stable4]], stable4[[#, "singularDecoratedGraphs"]] > 0 &], ord},
   ord = Table[{Module[{p = Position[e[[5]], _?(# =!= 0 &), {1}, 1, Heads -> False]}, If[p === {}, 0, 7 - p[[1, 1]]]], Total[Max[2 # - 2, 0] & /@ stable4[[e[[1, 1]], "genera"]]]}, {e, entries4}];
   {ge2 === sg, Length[ge2], AllTrue[ord, #[[1]] <= #[[2]] &], Sort[Tally[ord]]}]]
(* the total F_{4,0} agrees with the untwisted genus-4 cluster result (q -> 0 part of genRmat_F40_cluster_result.m) *)
If[ValueQ[stable4] && FileExistsQ["ref/genRmat_F40_cluster_result.m"], Module[{ref4 = Get["ref/genRmat_F40_cluster_result.m"][[1]] /. {q0 -> 0, q1 -> 0, q2 -> 0, q3 -> 0, q4 -> 0}, tot4 = Together[Total[stable4[[All, "term"]]]]},
   SeedRandom[4]; Table[Together[(tot4 - ref4) /. {y -> RandomInteger[{1, 50}]/RandomInteger[{51, 97}], XX -> RandomInteger[{-30, 30}]/RandomInteger[{1, 17}], e1 -> RandomInteger[{1, 40}], e2 -> RandomInteger[{1, 40}]}], {5}]]]
