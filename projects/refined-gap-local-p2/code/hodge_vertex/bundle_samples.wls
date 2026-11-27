(* bundle_samples.wls -- collect the per-point genus-3 sample files hv3_samples_B_p<idx>.mx into one sidecar hv3_samples_B.mx
   holding the list recs3 (in the order of the point list), as read by the notebook. Usage: wolframscript -file bundle_samples.wls *)
SetDirectory[DirectoryName[$InputFileName]];
files = FileNames["hv3_samples_B_p*.mx"];
idx[f_] := ToExpression[StringReplace[FileNameTake[f], {"hv3_samples_B_p" -> "", ".mx" -> ""}]];
files = SortBy[files, idx];
recs3 = Table[Get[f]; rec, {f, files}];
Print[Length[recs3], " records; indices ", Min[idx /@ files], "..", Max[idx /@ files], "; all total-ref zero: ", Union[Flatten[("total-ref" /. #[[4]]) & /@ recs3]] === {0}];
DumpSave["hv3_samples_B.mx", recs3];
Print["written hv3_samples_B.mx (", Round[FileByteCount["hv3_samples_B.mx"]/1024], " kB)"];
