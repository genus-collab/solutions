(* lp2_prep.m -- enumerate the decorated stable graphs (orbits) of genus <genus> for the local P2 pipeline and store them.
   Usage:  wolframscript -file lp2_prep.m genus
   Output: lp2_orbits_g<genus>.mx  (prepared, OT) read by the sampling and reconstruction jobs; counts in lp2_prep.log. *)
SetDirectory[DirectoryName[$InputFileName]];
gg = ToExpression[$ScriptCommandLine[[2]]];
Get["lp2_lib.m"];
log[m_] := Module[{s = OpenAppend["lp2_prep.log"]}, WriteString[s, DateString[{"Year", "-", "Month", "-", "Day", " ", "Hour", ":", "Minute", ":", "Second"}], "  [g=", gg, "] ", m, "\n"]; Close[s]; Print[m]];
t0 = AbsoluteTime[];
setupOrbits[gg];
DumpSave["lp2_orbits_g" <> ToString[gg] <> ".mx", {prepared, OT}];
log[ToString[Length[prepared]] <> " graphs, " <> ToString[Total[Length[#["decs"]] & /@ prepared]] <> " labelled decorations, " <> ToString[Length[OT["ids"]]] <> " decorated graphs, " <> ToString[Count[Values[OT["descriptions"]], d_ /; d["dilatonFree"]]] <> " dilaton-free; " <> ToString[Length[pointList[gg]]] <> " sample points (" <> ToString[Length[tsamplesFor[gg]]] <> " t x " <> ToString[Length[xsamplesFor[gg]]] <> " eps-ratios); " <> ToString[Round[AbsoluteTime[] - t0, 0.1]] <> " s"];
