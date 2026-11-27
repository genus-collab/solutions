(* hvc_summary.wls -- human-readable summary of one (genus, mode) run: per stable graph the polar part and the term, per decorated
   graph the key, the polar part and the term (coefficients as rational functions of Q and x, eps1 = x, eps2 = 1).
   Usage:  wolframscript -file hvc_summary.wls <genus> <tag> [brief]      Output: hvc_summary_g<genus>_<tag>.txt
   With "brief" the terms of the decorated graphs are omitted (polar parts only; for genus 4 the full terms would be gigabytes). *)
SetDirectory[DirectoryName[$InputFileName]];
{gg, tag} = {ToExpression[$ScriptCommandLine[[2]]], $ScriptCommandLine[[3]]};
brief = Length[$ScriptCommandLine] >= 4 && $ScriptCommandLine[[4]] == "brief";
Get["hvc_lib.m"];
resBase = "hvc_results_g" <> ToString[gg] <> "_" <> tag;
resDir = If[FileExistsQ[resBase <> ".m"] || FileExistsQ[resBase <> ".mx"], ".", "stanage_genus4"];   (* cluster results copied back *)
If[FileExistsQ[FileNameJoin[{resDir, resBase <> ".mx"}]], Get[FileNameJoin[{resDir, resBase <> ".mx"}]]; res = hvcResults, res = Get[FileNameJoin[{resDir, resBase <> ".m"}]]];
If[Head[res] =!= Association, Print["could not load ", resBase]; Exit[1]];
stable = res["stable"]; entries = res["entries"]; ngr = Length[stable];
out = {}; emit[s_] := AppendTo[out, s];
regQ[v_] := Union[Together /@ v] === {0};
emit[res["readme"]];
emit["Polar parts are the coefficients of t^-" <> ToString[2 gg - 2] <> ", ..., t^-1 at Q = 1 in t = Log Q, with eps1 = x, eps2 = 1. Terms are exact rational functions of Q and x."];
emit["total = F" <> ToString[gg] <> "0: " <> ToString[res["totalOK"]] <> "; total polar part " <> ToString[InputForm[Factor /@ res["polar"]]] <> "; equal to the gap P_g Li_{3-2g} polar part: " <> ToString[res["gapOK"]]];
emit[""];
Do[Module[{st = stable[[k]], es = Select[entries, #[[1, 1]] == k &]},
   emit["=== stable graph " <> ToString[k] <> ": genera " <> ToString[st["genera"]] <> ", edges " <> ToString[st["edges"]] <> ", |Aut| = " <> ToString[st["Aut"]] <> ", " <> ToString[st["decoratedGraphs"]] <> " decorated graphs (" <> ToString[st["singularDecoratedGraphs"]] <> " singular at Q = 1)"];
   emit["  polar part = " <> If[regQ[st["polar"]], "regular", ToString[InputForm[Factor /@ st["polar"]]]]];
   If[! brief, emit["  term = " <> ToString[InputForm[Factor[st["term"]]]]]];
   Do[emit["  -- " <> ToString[e[[1]]] <> " psi " <> ToString[e[[2, 1]]] <> " dil " <> ToString[e[[2, 2]]] <> "  polar " <> If[regQ[e[[5]]], "regular", ToString[InputForm[Factor /@ e[[5]]]]]];
      If[! brief, emit["       term " <> ToString[InputForm[Factor[e[[4]]]]]]], {e, es}];
   emit[""]], {k, ngr}];
Export[FileNameJoin[{resDir, "hvc_summary_g" <> ToString[gg] <> "_" <> tag <> ".txt"}], StringRiffle[out, "\n"], "Text"];
Print["written ", FileNameJoin[{resDir, "hvc_summary_g" <> ToString[gg] <> "_" <> tag <> ".txt"}], " (", Length[out], " lines)"];
