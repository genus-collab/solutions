(* lp2_reshape.m -- regroup the sampling output by decorated graph: one chunk file per reconstruction worker.
   Usage:  wolframscript -file lp2_reshape.m genus nChunks
   Input:  lp2_samples_g<genus>_w*.mx (any number of workers); a t-value is kept only if all its eps-ratio points are present.
   Output: lp2_chunk_g<genus>_c<c>of<n>.mx with chunkRecs[{t0, x0}] = vectors of the chunk's orbits, chunkIdx = their positions in
           OT["ids"], tsOK = complete t-values, xsAll = eps-ratios; log in lp2_reshape.log. Orbit i goes to chunk Mod[i-1, n]+1. *)
SetDirectory[DirectoryName[$InputFileName]];
args = Rest[$ScriptCommandLine];
{gg, nchunks} = ToExpression /@ args[[1 ;; 2]];
Get["lp2_lib.m"];
log[m_] := Module[{s = OpenAppend["lp2_reshape.log"]}, WriteString[s, DateString[{"Year", "-", "Month", "-", "Day", " ", "Hour", ":", "Minute", ":", "Second"}], "  [g=", gg, "] ", m, "\n"]; Close[s]; Print[m]];
loadOrbits[gg];
norb = Length[OT["ids"]];
files = FileNames["lp2_samples_g" <> ToString[gg] <> "_w*_p*.mx"];
If[files === {}, log["no sample files found"]; Exit[1]];
allRecs = <||>; nbad = 0; imax = 0;
Do[rec = Null; Get[f];
   If[rec === Null, log["ERROR: " <> f <> " holds no record"]; Continue[]];
   If[Length[rec[[2]]] != norb, log["ERROR: record " <> ToString[rec[[1]]] <> " in " <> f <> " has " <> ToString[Length[rec[[2]]]] <> " orbits instead of " <> ToString[norb]]; Continue[]];
   allRecs[rec[[1]]] = rec[[2]]; imax = Max[imax, "imag" /. rec[[3]]]; If[("bad" /. rec[[3]]) =!= {}, nbad++], {f, files}];
xsAll = xsamplesFor[gg];
tsOK = Select[tsamplesFor[gg], And @@ Table[KeyExistsQ[allRecs, {#, x0}], {x0, xsAll}] &];
log[ToString[Length[files]] <> " files, " <> ToString[Length[allRecs]] <> " of " <> ToString[Length[pointList[gg]]] <> " points; max |Im| " <> ToString[imax] <> "; points with unexpected monomials: " <> ToString[nbad] <> "; complete t-values: " <> ToString[Length[tsOK]] <> " of " <> ToString[Length[tsamplesFor[gg]]] <> "; incomplete: " <> ToString[Complement[tsamplesFor[gg], tsOK]]];
If[Length[tsOK] < 8, log["too few complete t-values"]; Exit[1]];
Do[chunkIdx = Select[Range[norb], Mod[# - 1, nchunks] + 1 == c &];
   chunkRecs = Association[Table[{t0, x0} -> allRecs[{t0, x0}][[chunkIdx]], {t0, tsOK}, {x0, xsAll}]];
   DumpSave["lp2_chunk_g" <> ToString[gg] <> "_c" <> ToString[c] <> "of" <> ToString[nchunks] <> ".mx", {chunkRecs, chunkIdx, tsOK, xsAll}], {c, nchunks}];
log["wrote " <> ToString[nchunks] <> " chunk files (" <> ToString[Ceiling[norb/nchunks]] <> " orbits each at most)"];
