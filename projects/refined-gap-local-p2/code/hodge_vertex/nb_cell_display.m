(* the genus-3 vertex graph: decomposition of its polar part over its decorated graphs (dilaton-leaf patterns) *)
Grid[Prepend[Table[{id, OTH3["descriptions"][id]["dilatonLeaves"], Factor /@ polarOrbits3[id]}, {id, Select[Keys[terms3], #[[1]] == 1 &]}], {"id", "dilaton leaves", "polar part"}], Frame -> All]
(* which stable graphs contain singular decorated graphs (expected: exactly those with a vertex of genus >= 2), and whether the
   cancellation inside such a graph needs the full sum over its decorations (grouping by psi assignment or by dilaton pattern alone) *)
regQ[v_] := Union[v] === {0};
singGraphs = Select[Range[Length[preparedH3]], Function[k, AnyTrue[Select[Keys[polarOrbits3], #[[1]] == k &], ! regQ[polarOrbits3[#]] &]]];
{singGraphs, Table[{k, preparedH3[[k]]["g"], preparedH3[[k]]["edges"]}, {k, singGraphs}], Select[Range[Length[preparedH3]], Max[preparedH3[[#]]["g"]] >= 2 &]}
Table[Module[{ids = Select[Keys[polarOrbits3], #[[1]] == k &]}, {k, Count[GatherBy[ids, OTH3["descriptions"][#]["psiOnHalfEdges"] &], g_ /; ! regQ[Total[polarOrbits3 /@ g]]], Count[GatherBy[ids, OTH3["descriptions"][#]["dilatonLeaves"] &], g_ /; ! regQ[Total[polarOrbits3 /@ g]]]}], {k, singGraphs}]
(* the terms of a few edge graphs as polynomials in XX (E1 = eps1, E2 = eps2) *)
showTerm[term_] := Module[{tm = Together[term]}, Column[Table[Row[{"XX^", j, " : ", Factor[Coefficient[Numerator[tm], XX, j]/Denominator[tm]] /. {e1 -> E1, e2 -> E2}}], {j, 0, Exponent[Numerator[tm], XX]}]]];
Table[{k, preparedH3[[k]]["g"], preparedH3[[k]]["edges"], showTerm[graphTerms3[[k]]]}, {k, Select[Range[Length[preparedH3]], Length[preparedH3[[#]]["g"]] >= 3 &]}] // TableForm
