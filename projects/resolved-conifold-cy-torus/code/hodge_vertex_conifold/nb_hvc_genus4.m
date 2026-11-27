(* genus 4 (mode B) is computed on Stanage with cluster_genus4_hodge_conifold/ (379 stable graphs, 27172 decorated graphs, exact in Q and x)
   and copied back into the subfolder stanage_genus4/ of the data directory: hvc_results_g4_Bm.m (per stable graph and per decorated
   graph), hvc_table_g4_Bm.txt, hvc_analysis_g4_Bm.txt. This section reads them if present. *)
g4file = FileNameJoin[{"stanage_genus4", "hvc_results_g4_Bm.m"}];
If[! FileExistsQ[g4file], Print["genus-4 results not present (", g4file, ") -- run cluster_genus4_hodge_conifold/ on Stanage and copy hvc_results_g4_Bm.m back into stanage_genus4/"], res4 = Get[g4file]; Print[Length[res4["stable"]], " stable graphs, ", Length[res4["entries"]], " decorated graphs; total = F40 (closed form): ", res4["totalOK"], "; total polar part = gap: ", res4["gapOK"]]]
(* total polar part {t^-6, ..., t^-1} versus the gap polar part of P_4 Li_{-5}(Q)/(eps1 eps2) *)
If[ValueQ[res4], {Factor /@ res4["polar"], gapHVC[4, x]}]
(* which stable graphs are singular at Q = 1 (expected: only the genus-4 vertex graph, with exactly the gap) *)
If[ValueQ[res4], Module[{sing = Select[Range[Length[res4["stable"]]], ! regQ[res4["stable"][[#]]["polar"]] &]}, {Length[res4["stable"]] - Length[sing], "regular of", Length[res4["stable"]], "; singular:", Table[{k, res4["stable"][[k]]["genera"], res4["stable"][[k]]["edges"], Factor /@ res4["stable"][[k]]["polar"]}, {k, sing}]}]]
(* the structural statements: singular decorated graphs only in the stable graphs with a vertex of genus >= 2, pole order <= Sum_v max(2 g_v - 2, 0) *)
If[ValueQ[res4], structure[res4]]
If[ValueQ[res4], Grid[Prepend[Table[With[{st = res4["stable"][[k]]}, {k, st["genera"], st["edges"], st["decoratedGraphs"], st["singularDecoratedGraphs"]}], {k, Select[Range[Length[res4["stable"]]], res4["stable"][[#]]["singularDecoratedGraphs"] > 0 &]}], {"graph", "genera", "edges", "decorated graphs", "singular decorated graphs"}], Frame -> All]]
