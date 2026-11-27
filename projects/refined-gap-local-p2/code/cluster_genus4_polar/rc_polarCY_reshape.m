(* rc_polarCY_reshape.m -- regroup the CY-torus sampling output by decorated graph: one chunk file per fitting worker.
   Usage:  wolframscript -file rc_polarCY_reshape.m genus nChunks
   Input:  rc_polarCY_g<genus>_w*.mx (any number of workers); points whose orbit sum failed the P_g Li check are dropped.
   Output: rc_polarCY_chunk_g<genus>_c<c>of<n>.mx with chunkRecs[point] = matrix rows of the chunk's orbits, chunkIdx, ptsAll; log in rc_polarCY_reshape.log *)
SetDirectory[DirectoryName[$InputFileName]];
{gg, nchunks} = ToExpression /@ $ScriptCommandLine[[2 ;; 3]];
Get["rc_lib.m"];
log[m_] := Module[{s = OpenAppend["rc_polarCY_reshape.log"]}, WriteString[s, DateString[{"Year", "-", "Month", "-", "Day", " ", "Hour", ":", "Minute", ":", "Second"}], "  [g=", gg, "] ", m, "\n"]; Close[s]; Print[m]];
loadOrbits[gg]; norb = Length[OT["ids"]];
files = FileNames["rc_polarCY_g" <> ToString[gg] <> "_w*_p*.mx"];
If[files === {}, log["no sample files found"]; Exit[1]];
allRecs = <||>; nbad = 0;
Do[rec = Null; Get[f];
   If[rec === Null, log["ERROR: " <> f <> " holds no record"]; Continue[]];
   If[Dimensions[rec[[2]]] =!= {norb, 2 gg - 2}, log["ERROR: record " <> ToString[rec[[1]]] <> " in " <> f <> " has dimensions " <> ToString[Dimensions[rec[[2]]]]]; Continue[]];
   If[! rec[[3]], nbad++; Continue[]];
   allRecs[rec[[1]]] = rec[[2]], {f, files}];
ptsAll = Keys[allRecs];
log[ToString[Length[files]] <> " files, " <> ToString[Length[ptsAll]] <> " usable points, " <> ToString[nbad] <> " dropped (orbit sum check failed)"];
Do[chunkIdx = Select[Range[norb], Mod[# - 1, nchunks] + 1 == c &];
   chunkRecs = Association[Table[p -> allRecs[p][[chunkIdx]], {p, ptsAll}]];
   DumpSave["rc_polarCY_chunk_g" <> ToString[gg] <> "_c" <> ToString[c] <> "of" <> ToString[nchunks] <> ".mx", {chunkRecs, chunkIdx, ptsAll}], {c, nchunks}];
log["wrote " <> ToString[nchunks] <> " chunk files"];
