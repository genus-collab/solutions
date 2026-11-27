(* explore_pointframe.wls -- structure of the stored conifold references and the per-graph regularity in the point frame
   at alpha0 = 0, alpha_i = eps_i (where the rank-4 Hodge vertex is trivial). *)
SetDirectory[DirectoryName[$InputFileName]];
tt = AbsoluteTime[];
res = Get["ref/resolvedconifold_results.m"];
Print["resolvedconifold_results.m: Head ", Head[res], ", Length ", Length[res], ", variables ", Union[Cases[res, _Symbol, {0, Infinity}] /. {Plus -> Nothing, Times -> Nothing, Power -> Nothing, Rational -> Nothing, Integer -> Nothing, List -> Nothing}] , ", symbols with heads ", Union[Cases[res, (\[Alpha] | \[Epsilon])[_], {0, Infinity}]]];
Print["LeafCounts: ", LeafCount /@ res];
e1 = \[Epsilon][1]; e2 = \[Epsilon][2]; a0 = \[Alpha][0]; a1 = \[Alpha][1]; a2 = \[Alpha][2];
(* CY torus, fibrewise, alpha1 = al: check F20 - P2 Li_{-1} structure *)
Pg[g_] := Module[{tv}, SeriesCoefficient[Series[e1 e2 tv^2/(4 Sinh[e1 tv/2] Sinh[e2 tv/2]), {tv, 0, 2 g}], 2 g]];
Do[Module[{F = res[[g]] /. {a0 -> 0, a1 -> al, a2 -> e1 + e2 - al}, inst},
  inst = Together[F - (F /. Q -> 0)];
  Print["genus ", g, ": F_g(Q) - F_g(0) on the fibrewise CY torus minus P_g Li_{3-2g}(Q)/(e1 e2): ", Together[inst - Pg[g] PolyLog[3 - 2 g, Q]/(e1 e2)] === 0,
    ";  F_g(0) = ", Factor[F /. Q -> 0]]], {g, 2, 3}];
(* the stored per-graph terms at alpha0 = 0, alpha_i = eps_i *)
gt = Get["ref/conifold_graphterms_genus2_3.m"];
Print["conifold_graphterms_genus2_3.m: ", Head[gt], " of length ", Length[gt], "; part 2 keys: ", Keys[Association[gt[[2]]]]];
dat = Association[gt[[2]]];
Do[Module[{ents = dat[g], ngr, stable, pol, reg, singDec, lead, ords},
  ngr = Max[ents[[All, 1, 1]]];
  Print["genus ", g, ": ", Length[ents], " decorated graphs, ", ngr, " stable graphs; entry structure: ", Length /@ ents[[1]], "; description keys: ", Keys[Association[ents[[1, 2]]]]];
  pol = Table[Together /@ Total[Select[ents, #[[1, 1]] == k &][[All, 4]]], {k, ngr}];
  reg = Select[Range[ngr], Union[pol[[#]]] === {0} &];
  singDec = Table[Count[Select[ents, #[[1, 1]] == k &], e_ /; Union[Together /@ e[[4]]] =!= {0}], {k, ngr}];
  ords = Table[Module[{p = Position[Together /@ e[[4]], _?(# =!= 0 &), {1}, 1, Heads -> False]}, If[p === {}, 0, 2 g - 2 - (p[[1, 1]] - 1)]], {e, ents}];
  Print["  stable graphs regular at Q = 1: ", Length[reg], " of ", ngr, " -> ", reg];
  Print["  singular decorated graphs per stable graph: ", singDec];
  Print["  pole-order tally over decorated graphs: ", Sort[Tally[ords]]];
  Print["  total polar part: ", Factor /@ Together /@ Total[pol], "; expected ", Factor /@ Table[SeriesCoefficient[Series[Pg[g] PolyLog[3 - 2 g, Q]/(e1 e2) /. Q -> Exp[t], {t, 0, -1}], -n], {n, 2 g - 2, 1, -1}]];
  Print["  graph 1 (genus-", g, " vertex): genera ", "genera" /. ents[[1, 2]], ", polar ", Factor /@ pol[[1]]];
  Print["  per stable graph {k, genera, edges, leading polar coefficient}: "];
  Do[Print["    ", {k, "genera" /. First[Select[ents, #[[1, 1]] == k &]][[2]], "edges" /. First[Select[ents, #[[1, 1]] == k &]][[2]], Factor[pol[[k, 1]]]}], {k, ngr}]], {g, 2, 3}];
Print["done in ", Round[AbsoluteTime[] - tt], " s"];
