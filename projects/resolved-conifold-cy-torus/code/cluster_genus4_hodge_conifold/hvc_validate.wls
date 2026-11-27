(* hvc_validate.wls -- compare a fresh run (hvc_results_g<genus>_Bm.m in this directory) with the stored local run ref/hvc_results_g<genus>_Bm.m:
   per decorated graph (matched by stable graph and orbit key {psi on half-edges, dilaton leaves}) the exact term and the polar part,
   and per stable graph the term and the polar part.  Usage:  wolframscript -file hvc_validate.wls <genus>
   Appends VALIDATION OK / VALIDATION FAIL lines to hvc_validate.log. *)
SetDirectory[DirectoryName[$InputFileName]];
gg = ToExpression[$ScriptCommandLine[[2]]];
log[m_] := Module[{st = OpenAppend["hvc_validate.log"]}, WriteString[st, DateString[{"Year", "-", "Month", "-", "Day", " ", "Hour", ":", "Minute", ":", "Second"}], "  [g=", gg, "] ", m, "\n"]; Close[st]; Print[m]];
new = Get["hvc_results_g" <> ToString[gg] <> "_Bm.m"]; ref = Get["ref/hvc_results_g" <> ToString[gg] <> "_Bm.m"];
keyOf[e_] := {e[[1, 1]], e[[2, 1]], e[[2, 2]]};
refMap = Association[Table[keyOf[e] -> e, {e, ref["entries"]}]];
nt = 0; np = 0; nmiss = 0;
Do[Module[{k = keyOf[e], r},
   If[! KeyExistsQ[refMap, k], nmiss++, r = refMap[k];
    If[Together[e[[4]] - r[[4]]] === 0, nt++];
    If[Union[Together /@ (e[[5]] - r[[5]])] === {0}, np++]]], {e, new["entries"]}];
ok1 = nt == Length[new["entries"]] && np == Length[new["entries"]] && nmiss == 0 && Length[ref["entries"]] == Length[new["entries"]];
log[If[ok1, "VALIDATION OK: all " <> ToString[Length[new["entries"]]] <> " decorated graphs (terms and polar parts) agree with the stored local run",
   "VALIDATION FAIL: decorated graphs: terms " <> ToString[nt] <> ", polar parts " <> ToString[np] <> " of " <> ToString[Length[new["entries"]]] <> " agree (stored " <> ToString[Length[ref["entries"]]] <> ", unmatched keys " <> ToString[nmiss] <> ")"]];
ok2 = Length[new["stable"]] == Length[ref["stable"]] && And @@ Table[Together[new["stable"][[k]]["term"] - ref["stable"][[k]]["term"]] === 0 && Union[Together /@ (new["stable"][[k]]["polar"] - ref["stable"][[k]]["polar"])] === {0}, {k, Length[new["stable"]]}];
log[If[ok2, "VALIDATION OK: all " <> ToString[Length[new["stable"]]] <> " stable-graph terms and polar parts agree with the stored local run", "VALIDATION FAIL: stable-graph terms or polar parts differ from the stored local run"]];
log[If[new["totalOK"] && new["gapOK"], "VALIDATION OK: fresh run reports total = F" <> ToString[gg] <> "0 and the gap", "VALIDATION FAIL: fresh run reports totalOK = " <> ToString[new["totalOK"]] <> ", gapOK = " <> ToString[new["gapOK"]]]];
