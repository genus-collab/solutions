(* rc_prep.m -- enumerate the decorated stable graphs (orbits) of genus <genus> for the conifold pipeline and store them.
   Usage:  wolframscript -file rc_prep.m genus
   Output: rc_orbits_g<genus>.mx (prepared, OT); counts in rc_prep.log. (Same enumeration as lp2_prep.m, so the keys match.) *)
SetDirectory[DirectoryName[$InputFileName]];
gg = ToExpression[$ScriptCommandLine[[2]]];
Get["rc_lib.m"];
log[m_] := Module[{s = OpenAppend["rc_prep.log"]}, WriteString[s, DateString[{"Year", "-", "Month", "-", "Day", " ", "Hour", ":", "Minute", ":", "Second"}], "  [g=", gg, "] ", m, "\n"]; Close[s]; Print[m]];
t0 = AbsoluteTime[];
setupOrbits[gg];
DumpSave["rc_orbits_g" <> ToString[gg] <> ".mx", {prepared, OT}];
log[ToString[Length[prepared]] <> " graphs, " <> ToString[Total[Length[#["decs"]] & /@ prepared]] <> " labelled decorations, " <> ToString[Length[OT["ids"]]] <> " decorated graphs, " <> ToString[Count[Values[OT["descriptions"]], d_ /; d["dilatonFree"]]] <> " dilaton-free; " <> ToString[Round[AbsoluteTime[] - t0, 0.1]] <> " s"];
