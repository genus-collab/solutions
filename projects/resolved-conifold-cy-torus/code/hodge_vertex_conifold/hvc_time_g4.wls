(* hvc_time_g4.wls -- timing of a sample of genus-4 decorated graphs (mode B) with the exact route (Q and x symbolic), to size the
   cluster run. Needs hvc_factors_g4_Bm.mx and hvc_orbits_g4_hodge.mx (from hvc_run.wls prep 4 B m).
   Usage: wolframscript -file hvc_time_g4.wls [limitSeconds]     Output: hvc_time_g4.txt *)
SetDirectory[DirectoryName[$InputFileName]];
limit = If[Length[$ScriptCommandLine] >= 2, ToExpression[$ScriptCommandLine[[2]]], 900];
Get["hvc_lib.m"];
tt = AbsoluteTime[]; Get["hvc_factors_g4_Bm.mx"]; loadOrbitsHVC[4, "B"];
Print["loaded factors and ", Length[OT["ids"]], " decorated graphs in ", Round[AbsoluteTime[] - tt, 0.1], " s; LeafCount EdN[1,2][8,0] = ", LeafCount[EdN[1, 2][8, 0]], ", EdN[1,1][4,4] = ", LeafCount[EdN[1, 1][4, 4]]];
usub = usubMode["B", x];
(* sample: the first orbit of a spread of stable graphs, the last orbit of graph 1 (largest dilaton pattern), and the graphs with most vertices/edges *)
nV = Table[Length[prepared[[k]]["g"]], {k, Length[prepared]}]; nE = Table[Length[prepared[[k]]["edges"]], {k, Length[prepared]}];
pick = Join[Table[First[Select[OT["ids"], #[[1]] == k &]], {k, {1, 2, 3, 7, 11, 20, 50, 100, 150, 200, 250, 300, 350, 379}}], {Last[Select[OT["ids"], #[[1]] == 1 &]]},
   Table[First[Select[OT["ids"], #[[1]] == k &]], {k, Take[Ordering[nE, -3], 3]}]];
pick = DeleteDuplicates[pick];
out = {}; emit[s_] := (AppendTo[out, s]; Print[s]);
emit["genus 4, mode B, exact route (Q, x symbolic): " <> ToString[Length[pick]] <> " sample decorated graphs, limit " <> ToString[limit] <> " s each"];
Do[Module[{gr = OT["graphs"][id], t0 = AbsoluteTime[], r},
   r = TimeConstrained[orbitTerm[gr, usub, 4], limit, $Aborted];
   emit[ToString[id] <> " genera " <> ToString[gr["g"]] <> " edges " <> ToString[Length[gr["edges"]]] <> " key " <> ToString[OT["orbits"][[id[[1]], id[[2]]]]] <> " decorations " <> ToString[Length[gr["decs"]]] <> ": " <>
     If[r === $Aborted, "ABORTED after " <> ToString[limit] <> " s", ToString[Round[AbsoluteTime[] - t0, 0.1]] <> " s, LeafCount " <> ToString[LeafCount[r[[1]]]] <> ", symmetric " <> ToString[r[[3]]] <> ", polar " <> ToString[InputForm[r[[2]]]]]]], {id, pick}];
Export["hvc_time_g4.txt", StringRiffle[out, "\n"], "Text"];
