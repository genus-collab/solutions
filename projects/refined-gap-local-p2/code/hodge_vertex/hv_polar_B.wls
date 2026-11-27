(* hv_polar_B.wls -- conifold polar parts of every decorated graph (orbit) in mode B, and the residue mechanism *)
SetDirectory[DirectoryName[$InputFileName]];
Get["lp2_lib.m"]; tcSetup[8];
res = Get["hodge_vertex_genus2_results.m"];
e1 = \[Epsilon][1]; e2 = \[Epsilon][2];
Print["XX at the conifold (y = -1/27): ", Normal[XXcf] /. tcv -> 0, ";  XX series: ", Normal[Series[XXcf, {tcv, 0, 2}]]];
Do[Module[{B = res["B"], ids},
  ids = Keys[B["orbitTerms"]];
  Print["mode B, per-orbit polar parts {t_c^-2, t_c^-1} (eps1 = x, eps2 = 1); graph k = id[[1]]:"];
  Do[Module[{d = B["descriptions"][id], p = polarOfTerm[B["orbitTerms"][id], 2]},
    p = p /. SeriesCoefficient[a_, _] /; FreeQ[a, tcv] :> 0;
    Print["  ", id, "  genera ", d["genera"], " edges ", d["edges"], " psi ", d["psiOnHalfEdges"], " dil ", d["dilatonLeaves"], " -> ", InputForm[Factor /@ p]]], {id, ids}];
  Print["sum over the 7 orbits of the vertex graph: ", InputForm[Factor /@ Total[Table[polarOfTerm[B["orbitTerms"][id], 2], {id, Select[ids, #[[1]] == 1 &]}]]]];
  ], {1}];
