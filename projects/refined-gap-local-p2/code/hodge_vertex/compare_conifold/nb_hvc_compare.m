(* comparison with local P2 x C2 in the same Hodge-vertex frame (mode B, hodge_vertex/hodge_vertex_genus{2,3}_results.m): the stable
   graphs and the orbit keys {psi powers on half-edges, dilaton leaves} are identical on both sides; local P2 polar parts are in the
   conifold coordinate t_c (t_c^2 = -3 t^2 at the leading order), the conifold ones in t = Log Q *)
lp2dir = FileNameJoin[{ParentDirectory[Directory[]], "hodge_vertex"}];
conMap[gg_] /; gg <= 3 := conMap[gg] = Association[Table[{e[[1, 1]], e[[2, 1]], e[[2, 2]]} -> (Together /@ e[[5]]), {e, If[gg == 2, res2["B"], res3["B"]]["entries"]}]];
compareWith[gg_, lpMap_] := Module[{cm = conMap[gg], common, rows},
   common = Intersection[Keys[cm], Keys[lpMap]];
   rows = Table[Module[{both, one, rat}, both = Select[common, cm[#][[j]] =!= 0 && lpMap[#][[j]] =!= 0 &]; one = Select[common, (cm[#][[j]] === 0) =!= (lpMap[#][[j]] === 0) &];
      rat = Union[Table[Factor[lpMap[k][[j]]/cm[k][[j]]], {k, both}]]; {"t^-" <> ToString[2 gg - 1 - j], Length[both], Length[one], rat}], {j, 2 gg - 2}];
   <|"decorated graphs" -> {Length[cm], Length[lpMap], "same keys" -> (Sort[Keys[cm]] === Sort[Keys[lpMap]])}, "regular on both sides" -> {Count[common, k_ /; Union[cm[k]] === {0}], Count[common, k_ /; Union[lpMap[k]] === {0}], Sort[Select[common, Union[cm[#]] === {0} &]] === Sort[Select[common, Union[lpMap[#]] === {0} &]]},
     "per order {order, both nonzero, exactly one nonzero, distinct ratios localP2/conifold}" -> rows, "expected leading ratio (-3)^(g-1)" -> (-3)^(gg - 1)|>];
(* genus 2: the local P2 polar parts of the decorated graphs are computed here from the stored orbit terms (t_c expansion of lp2_lib.m) *)
If[FileExistsQ[FileNameJoin[{lp2dir, "hodge_vertex_genus2_results.m"}]], Module[{r2 = Get[FileNameJoin[{lp2dir, "hodge_vertex_genus2_results.m"}]]["B"]}, tcSetup[8];
   lpMap2 = Association[Table[Module[{d = r2["descriptions"][id]}, {id[[1]], d["psiOnHalfEdges"], d["dilatonLeaves"]} -> (polarOfTerm[r2["orbitTerms"][id], 2] /. SeriesCoefficient[a_, _] /; FreeQ[a, tcv] :> 0)], {id, Keys[r2["orbitTerms"]]}]];
   compareWith[2, lpMap2]], "hodge_vertex/hodge_vertex_genus2_results.m not found"]
(* genus 3 *)
If[FileExistsQ[FileNameJoin[{lp2dir, "hodge_vertex_genus3_results.m"}]], Module[{r3 = Get[FileNameJoin[{lp2dir, "hodge_vertex_genus3_results.m"}]]},
   lpMap3 = Association[Table[Module[{d = r3["descriptions"][id]}, {id[[1]], d["psiOnHalfEdges"], d["dilatonLeaves"]} -> (Together /@ r3["polarOrbits"][id])], {id, Keys[r3["polarOrbits"]]}]];
   compareWith[3, lpMap3]], "hodge_vertex/hodge_vertex_genus3_results.m not found"]
(* the two universal ratios: -3 = t_c^2/t^2 at the leading order and (eps1+eps2)^2/(3 (eps1-eps2)^2) at the next order, as in the
   point-vertex comparison of 2026-09 (polar_check_genus2.log) *)
{(-3)^(2 - 1), Factor[(-3) (-(e1 + e2)^2/(9 (e1 - e2)^2)) /. {e1 -> x, e2 -> 1}]}
(* comparison with local P2 at genus 4 (hodge_vertex/hodge_graphterms_genus4.m, polar parts in t_c): same regular sets, ratios (-3)^3 = -27 at t^-6
   and (-3)^3 (-(eps1+eps2)^2/(9 (eps1-eps2)^2)) at t^-5 expected *)
If[ValueQ[res4] && FileExistsQ[FileNameJoin[{lp2dir, "hodge_graphterms_genus4.m"}]], Module[{r4 = Rest[Get[FileNameJoin[{lp2dir, "hodge_graphterms_genus4.m"}]]]},
   conMap[4] = Association[Table[{e[[1, 1]], e[[2, 1]], e[[2, 2]]} -> (Together /@ e[[5]]), {e, res4["entries"]}]];
   lpMap4 = Association[Table[{e[[1, 1]], e[[2, 1]], e[[2, 2]]} -> (Together /@ e[[5]]), {e, r4}]];
   compareWith[4, lpMap4]]]
