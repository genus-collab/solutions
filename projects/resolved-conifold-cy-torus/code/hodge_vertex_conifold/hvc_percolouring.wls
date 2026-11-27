(* hvc_percolouring.wls -- mode B, genus 3: is the sum over all decorations of a stable graph with a genus-2 vertex regular at Q = 1
   for each colouring separately (evaluated at the singular root, i.e. colour 1 = singular idempotent)? *)
SetDirectory[DirectoryName[$InputFileName]];
Get["hvc_lib.m"]; Get["hvc_factors_g3_Bm.mx"]; loadOrbitsHVC[3, "B"];
S0 = Together[SsumQ]; P0 = Together[PprodQ]; sq = ((1 + x)/x) Sqrt[1 + 4 x (1/Q - 1)/(1 + x)^2]; U1 = (S0 + sq)/2;
usub = usubMode["B", x];
graphValueCol[gr_, col_] := Module[{V = gr["V"], ed = gr["edges"], H = gr["H"], hsAt = gr["hsAt"], total = C[0, 0], hp, facs},
   Do[hp = ConstantArray[0, Length[H]];
    Do[Do[hp[[hsAt[[v, r]]]] = dec[[v, 1, r]], {r, Length[hsAt[[v]]]}], {v, V}];
    facs = Join[Table[deltaPow[col[[v]], gr["g"][[v]] - 1], {v, V}], Flatten[Table[Table[LfN[col[[v]]][b], {b, dec[[v, 2]]}], {v, V}]],
      Flatten[Table[Table[DilN[col[[v]]][d], {d, dec[[v, 3]]}], {v, V}]], Table[EdN[col[[ed[[e, 1]]]], col[[ed[[e, 2]]]]][hp[[2 e - 1]], hp[[2 e]]], {e, Length[ed]}]];
    total = cadd[total, cscale[Product[dec[[v, 4]], {v, V}], cmul @@ facs]], {dec, gr["decs"]}];
   cscale[1/gr["aut"], total]];
polarAt1[C[a_, b_]] := Module[{s = Series[(a + b U1) /. Q -> Exp[t], {t, 0, -1}]}, If[Head[s] === SeriesData, Table[Factor[SeriesCoefficient[s, -n]], {n, 4, 1, -1}], {0, 0, 0, 0}]];
Do[Module[{gr = prepared[[k]], V, cols},
   V = Length[gr["g"]]; cols = Tuples[{1, 2}, V];
   Print["graph ", k, " genera ", gr["g"], " edges ", gr["edges"], ":"];
   Do[Module[{val = C[0, 0]},
     Do[val = cadd[val, graphValueCol[numGraph[OT["graphs"][id], usub], col]], {id, Select[OT["ids"], #[[1]] == k &]}];
     Print["   colouring ", col, ": polar part at the singular root {t^-4..t^-1} = ", polarAt1[val]]], {col, cols}]], {k, {2, 5, 6}}];
