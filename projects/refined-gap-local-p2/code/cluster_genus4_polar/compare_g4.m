(* compare_g4.m -- decorated-graph-by-decorated-graph comparison of the conifold polar parts: local P2 (t_c) versus the resolved
   conifold at alpha0 = 0, alpha_i = eps_i (t = Log Q) and on the CY torus with alpha1, alpha2 free.
   Usage:  wolframscript -file compare_g4.m genus
   Input:  localP2_polar_genus<genus>.m, conifold_graphterms_genus<genus>.m, conifold_CY_graphterms_polar_genus<genus>.m (optional)
   Output: polar_check_genus<genus>.log (summary) and polar_check_genus<genus>_table.txt (per-graph table). *)
SetDirectory[DirectoryName[$InputFileName]];
gg = ToExpression[$ScriptCommandLine[[2]]];
e1 = \[Epsilon][1]; e2 = \[Epsilon][2]; A1 = \[Alpha][1]; A2 = \[Alpha][2];
out = {}; log[m_] := (AppendTo[out, m]; Print[m]);
hom[f_] := Factor[e2^(2 gg - 2) (f /. x -> e1/e2)];
P2 = Association[Table[e[[1]] -> e, {e, Rest[Get["localP2_polar_genus" <> ToString[gg] <> ".m"]]}]];
CF = Association[Table[e[[1]] -> e, {e, Rest[Get["conifold_graphterms_genus" <> ToString[gg] <> ".m"]]}]];
keys = Sort[Keys[P2]];
log["genus " <> ToString[gg] <> ": " <> ToString[Length[keys]] <> " decorated graphs (local P2), " <> ToString[Length[CF]] <> " (conifold); keys identical: " <> ToString[Keys[CF] === Keys[P2] || Sort[Keys[CF]] === keys] <> "; keys agree: " <> ToString[And @@ (P2[#][[2]] === CF[#][[2]] & /@ keys)]];
orders = Table[2 gg - 1 - n, {n, 1, 2 gg - 2}];   (* t^-(2g-2) ... t^-1 *)
ratio[p_, c_] := Which[c === 0 && p === 0, "0/0", c === 0, "conifold 0", True, Factor[p/c]];
rat = Table[ratio[hom[P2[k][[4, n]]], hom[CF[k][[4, n]]]], {k, keys}, {n, 2 gg - 2}];
Do[Module[{rs = Union[rat[[All, n]]]},
   log["order t^-" <> ToString[orders[[n]]] <> ": distinct localP2/conifold ratios over all decorated graphs: " <> ToString[Length[rs]] <> If[Length[rs] <= 3, "  " <> ToString[InputForm /@ rs], "  e.g. " <> ToString[InputForm /@ Take[rs, 2]]]]], {n, 2 gg - 2}];
log["expected universal factors: leading (-3)^(g-1) = " <> ToString[(-3)^(gg - 1)] <> ", next (-3)^(g-1) (-(e1+e2)^2/(9 (e1-e2)^2)) = " <> ToString[InputForm[(-3)^(gg - 1) (-(e1 + e2)^2/(9 (e1 - e2)^2))]]];
dilfree = Select[keys, ("dilatonFree" /. P2[#][[3]]) === True &];
log["dilaton-free decorated graphs: " <> ToString[Length[dilfree]]];
Do[Module[{rs = Union[Table[rat[[Position[keys, k][[1, 1]], n]], {k, dilfree}]]}, log["   order t^-" <> ToString[orders[[n]]] <> ", dilaton-free only: " <> ToString[Length[rs]] <> " distinct ratios"]], {n, 2 gg - 2}];
Do[log["dilaton-free: [t^-" <> ToString[orders[[n]]] <> "]/[t^-" <> ToString[2 gg - 2] <> "] local P2: " <> ToString[InputForm[Union[Table[Factor[hom[P2[k][[4, n]]]/hom[P2[k][[4, 1]]]], {k, dilfree}]]]] <> "   conifold: " <> ToString[InputForm[Union[Table[Factor[hom[CF[k][[4, n]]]/hom[CF[k][[4, 1]]]], {k, dilfree}]]]]], {n, 2, Min[3, 2 gg - 2]}];
(* regrouping by underlying stable graph *)
graphsOf = Union[keys[[All, 1]]];
Do[Module[{rs = Union[Table[ratio[Together[Total[Table[hom[P2[k][[4, n]]], {k, Select[keys, #[[1]] == gk &]}]]], Together[Total[Table[hom[CF[k][[4, n]]], {k, Select[keys, #[[1]] == gk &]}]]]], {gk, graphsOf}]]},
   log["sum over all decorations of each underlying graph (" <> ToString[Length[graphsOf]] <> " graphs), order t^-" <> ToString[orders[[n]]] <> ": " <> ToString[Length[rs]] <> " distinct ratios"]], {n, 2 gg - 2}];
totP2 = Factor /@ Total[Table[hom /@ P2[k][[4]], {k, keys}]]; totCF = Factor /@ Total[Table[hom /@ CF[k][[4]], {k, keys}]];
log["totals local P2: " <> ToString[InputForm[totP2]]];
log["totals conifold: " <> ToString[InputForm[totCF]]];
(* CY torus with free alphas *)
If[FileExistsQ["conifold_CY_graphterms_polar_genus" <> ToString[gg] <> ".m"],
  Module[{CY = Association[Table[e[[1]] -> e, {e, Rest[Get["conifold_CY_graphterms_polar_genus" <> ToString[gg] <> ".m"]]}]], ratCY, uu = (e1 + e2 - 2 A1) (e1 + e2 - 2 A2)},
   ratCY = Table[ratio[hom[P2[k][[4, n]]], CY[k][[4, n]]], {k, keys}, {n, 2 gg - 2}];
   Do[Module[{rs = Union[ratCY[[All, n]]]}, log["CY torus (alphas free), order t^-" <> ToString[orders[[n]]] <> ": distinct localP2/conifold ratios: " <> ToString[Length[rs]] <> If[Length[rs] <= 3, "  " <> ToString[InputForm /@ rs], "  e.g. " <> ToString[InputForm /@ Take[rs, 2]]]]], {n, 2 gg - 2}];
   log["expected next-order factor on the CY torus: (-3)^(g-1) (e1+e2)^2/(9 u) = " <> ToString[InputForm[(-3)^(gg - 1) (e1 + e2)^2/(9 uu)]]];
   Do[log["CY torus, order t^-" <> ToString[orders[[n]]] <> ": alpha-independent for " <> ToString[Count[Table[CY[k][[4, n]], {k, keys}], c_ /; FreeQ[c, A1 | A2]]] <> " graphs; divisible by u for " <> ToString[Count[Table[CY[k][[4, n]], {k, keys}], c_ /; Together[c /. A1 -> (e1 + e2)/2] === 0]]], {n, 2 gg - 2}]]];
Export["polar_check_genus" <> ToString[gg] <> ".log", StringRiffle[out, "\n"], "Text"];
(* table *)
tab = OpenWrite["polar_check_genus" <> ToString[gg] <> "_table.txt"];
WriteString[tab, "Genus ", gg, ": per decorated graph, local P2 polar coefficients {t_c^-", 2 gg - 2, ", ..., t_c^-1}, conifold (alpha0 = 0, alpha_i = eps_i) {t^-", 2 gg - 2, ", ..., t^-1}, ratios (eps1 = x, eps2 = 1)\n\n"];
Do[WriteString[tab, ToString[k], "  ", ToString[P2[k][[2]]], "  g=", ToString["genera" /. P2[k][[3]]], " edges=", ToString["edges" /. P2[k][[3]]], "\n   localP2:  ", ToString[InputForm[P2[k][[4]]]], "\n   conifold: ", ToString[InputForm[CF[k][[4]]]], "\n   ratios:   ", ToString[InputForm[rat[[Position[keys, k][[1, 1]]]]]], "\n"], {k, keys}];
Close[tab];
Print["written polar_check_genus", gg, ".log and polar_check_genus", gg, "_table.txt"];
