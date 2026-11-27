(* hv_table_g4.wls -- rewrite hodge_collect_genus4.txt from hodge_stablegraphs_genus4.m with the polar parts canonicalised
   (the collector's regQ tested the unsimplified sum of per-orbit polar parts, so 30 regular graphs were printed as
   "{0, 0, 0, 0, 0, 0}" instead of "regular"). Also prints the corrected counts. Run in hodge_vertex/. *)
SetDirectory[DirectoryName[$InputFileName]];
stable = Get["hodge_stablegraphs_genus4.m"];
regQ[v_] := Union[Together /@ v] === {0};
ngr = Length[stable];
regular = Select[Range[ngr], regQ[stable[[#]]["polar"]] &];
Print["stable graphs: ", ngr, "; regular at the conifold: ", Length[regular], "; singular: ", Complement[Range[ngr], regular]];
Print["stable graphs containing singular decorated graphs: ", Count[stable, s_ /; s["singularDecoratedGraphs"] > 0], "; with a vertex of genus >= 2: ", Count[stable, s_ /; Max[s["genera"]] >= 2]];
tab = StringRiffle[Prepend[Table[StringRiffle[{k, stable[[k]]["genera"], stable[[k]]["edges"], stable[[k]]["Aut"], stable[[k]]["decoratedGraphs"], stable[[k]]["singularDecoratedGraphs"],
      If[regQ[stable[[k]]["polar"]], "regular", ToString[InputForm[Factor /@ stable[[k]]["polar"]]]]}, "  "], {k, ngr}],
   "graph  genera  edges  |Aut|  decorated graphs  singular decorated graphs  polar part {t_c^-(2g-2), ..., t_c^-1} (eps1 = x, eps2 = 1)"], "\n"];
Export["hodge_collect_genus4.txt", tab, "Text"];
Print["rewritten hodge_collect_genus4.txt"];
