(* genus 3: 42 stable graphs, 693 decorated graphs with Hodge vertices (378 with point vertices). The stored results of hvc_run.wls
   (hvc_results_g3_Bm.m, _Am.m, _point.m) are read; set recompute3 = True to recompute everything here (about 5 minutes per mode
   on one kernel, exact in Q and x) *)
recompute3 = False;
res3[mode_] := res3[mode] = Module[{f = "hvc_results_g3_" <> modeTag[mode, -1] <> ".m"}, If[recompute3 || ! FileExistsQ[f], runMode[3, mode, -1], Get[f]]];
Table[{m, Length[res3[m]["entries"]], "total = F30" -> res3[m]["totalOK"], "total polar part = gap" -> res3[m]["gapOK"], Factor /@ res3[m]["polar"]}, {m, {"B", "A", "point"}}] // TableForm
(* mode B per stable graph: only the genus-3 vertex graph is singular at Q = 1 and it carries exactly the gap *)
tableOf[res3["B"]]
{Factor /@ res3["B"]["stable"][[1]]["polar"], gapHVC[3, x]}
(* mode B: the stable graphs with a vertex of genus >= 2 are exactly those containing singular decorated graphs, and the conifold
   pole order of every decorated graph is bounded by Sum_v max(2 g_v - 2, 0) *)
structure[res3["B"]]
(* the genus-3 vertex graph: decomposition of its polar part over its decorated graphs (dilaton-leaf patterns) *)
Grid[Prepend[Table[{e[[1]], e[[2, 2]], Factor /@ e[[5]]}, {e, Select[res3["B"]["entries"], #[[1, 1]] == 1 &]}], {"id", "dilaton leaves", "polar part {t^-4, ..., t^-1}"}], Frame -> All]
(* the other stable graphs with a genus-2 vertex: singular decorated graphs that cancel only in the sum over the decorations of the graph *)
Table[With[{st = res3["B"]["stable"][[k]]}, {k, st["genera"], st["edges"], st["decoratedGraphs"], st["singularDecoratedGraphs"], If[regQ[st["polar"]], "regular", Factor /@ st["polar"]]}], {k, Select[Range[42], Max[res3["B"]["stable"][[#]]["genera"]] >= 2 &]}] // TableForm
(* modes A and point: every stable graph singular *)
tableOf[res3["A"]]
tableOf[res3["point"]]
