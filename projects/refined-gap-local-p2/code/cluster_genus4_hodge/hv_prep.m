(* hv_prep.m -- enumerate the decorated stable graphs with Hodge vertices (orbits) of genus <genus> and store them.
   Usage:  wolframscript -file hv_prep.m genus
   Output: hv_orbits_g<genus>.mx (preparedH, OTH, hodgeClassG) read by the sampling, reconstruction and collection stages;
           counts in hv_prep.log. *)
SetDirectory[DirectoryName[$InputFileName]];
gg = ToExpression[$ScriptCommandLine[[2]]];
Get["hv_lib4.m"];
log[m_] := Module[{s = OpenAppend["hv_prep.log"]}, WriteString[s, DateString[{"Year", "-", "Month", "-", "Day", " ", "Hour", ":", "Minute", ":", "Second"}], "  [g=", gg, "] ", m, "\n"]; Close[s]; Print[m]];
t0 = AbsoluteTime[];
preparedH = prepareGraphH4 /@ graphsWithLegs[gg, 0]; OTH = orbitTable[preparedH]; hodgeClassG = hodgeClass4[gg];
DumpSave["hv_orbits_g" <> ToString[gg] <> ".mx", {preparedH, OTH, hodgeClassG}];
log[ToString[Length[preparedH]] <> " graphs, " <> ToString[Total[Length[#["decs"]] & /@ preparedH]] <> " labelled decorations, " <> ToString[Length[OTH["ids"]]] <> " decorated graphs (orbits), " <> ToString[Count[Values[OTH["descriptions"]], d_ /; d["dilatonFree"]]] <> " dilaton-free; " <> ToString[Length[pointListH[gg]]] <> " sample points (" <> ToString[Length[tsamplesFor[gg]]] <> " t x " <> ToString[nxFor[gg]] <> " eps-ratios); " <> ToString[Round[AbsoluteTime[] - t0, 0.1]] <> " s; written hv_orbits_g" <> ToString[gg] <> ".mx"];
