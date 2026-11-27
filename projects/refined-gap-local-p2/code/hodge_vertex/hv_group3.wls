SetDirectory[DirectoryName[$InputFileName]];
r = Get["hodge_vertex_genus3_results_nb.m"]; Get["hv_orbits_g3.mx"];
po = r["polarOrbits"]; desc = r["descriptions"];
regQ[v_] := Union[v] === {0};
Print["stable graphs with singular decorated graphs: ", Select[Range[42], Function[k, AnyTrue[Select[Keys[po], #[[1]] == k &], ! regQ[po[#]] &]]]];
(* grouping within each stable graph: by psi assignment on half-edges (summing over dilaton-leaf patterns), and by dilaton pattern *)
byPsi = Table[Module[{ids = Select[Keys[po], #[[1]] == k &], groups}, groups = GatherBy[ids, desc[#]["psiOnHalfEdges"] &]; {k, Length[groups], Count[groups, g_ /; ! regQ[Total[po /@ g]]]}], {k, 42}];
byDil = Table[Module[{ids = Select[Keys[po], #[[1]] == k &], groups}, groups = GatherBy[ids, desc[#]["dilatonLeaves"] &]; {k, Length[groups], Count[groups, g_ /; ! regQ[Total[po /@ g]]]}], {k, 42}];
Print["{graph, #psi-classes, #singular psi-classes} for graphs with any singular class: ", Select[byPsi, #[[3]] > 0 &]];
Print["{graph, #dilaton-classes, #singular dilaton-classes} for graphs with any singular class: ", Select[byDil, #[[3]] > 0 &]];
(* dilaton-free decorated graphs *)
df = Select[Keys[po], desc[#]["dilatonFree"] &];
Print["dilaton-free decorated graphs: ", Length[df], ", singular among them: ", Count[df, id_ /; ! regQ[po[id]]]];
(* orders of the poles of the singular decorated graphs of the edge graphs *)
Print["edge-graph singular orbits: max pole order t_c^-n present: ", Max[Table[Max[Position[Reverse[po[id]] /. (0 -> Nothing) /. {} -> {0}, _?(# =!= 0 &)][[All, 1]]], {id, Select[Keys[po], #[[1]] > 1 && ! regQ[po[#]] &]}]]];
Print["genus-3 vertex graph: ", Length[Select[Keys[po], #[[1]] == 1 &]], " decorated graphs, singular ", Count[Select[Keys[po], #[[1]] == 1 &], id_ /; ! regQ[po[id]]]];
Print["graph 2: ", preparedH3[[2]]["g"], preparedH3[[2]]["edges"], "; graph 5: ", preparedH3[[5]]["g"], preparedH3[[5]]["edges"], "; graph 6: ", preparedH3[[6]]["g"], preparedH3[[6]]["edges"]];
Print["list of graphs with singular orbits with genera/edges: ", Table[{k, preparedH3[[k]]["g"], preparedH3[[k]]["edges"]}, {k, Select[Range[42], Function[k, AnyTrue[Select[Keys[po], #[[1]] == k &], ! regQ[po[#]] &]]]}]];
