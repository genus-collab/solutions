(* lp2_recon.m -- reconstruction worker: exact contribution of every decorated graph of one chunk and its polar part
   in the conifold coordinate t_c (eps1 = x, eps2 = 1).
   Usage:  wolframscript -file lp2_recon.m genus chunk nChunks
   Input:  lp2_chunk_g<genus>_c<chunk>of<n>.mx (from lp2_reshape.m), lp2_orbits_g<genus>.mx
   Output: lp2_terms_g<genus>_c<chunk>of<n>.m   = list of {id, key, description, term(y, XX, eps1, eps2), polar {t_c^-(2g-2), ..., t_c^-1}}
           lp2_partial_g<genus>_c<chunk>of<n>.m = {sum of the terms, sum of the polar parts, list of failed ids}
           log in lp2_recon.log *)
SetDirectory[DirectoryName[$InputFileName]];
args = Rest[$ScriptCommandLine];
{gg, cc, nchunks} = ToExpression /@ args[[1 ;; 3]];
Get["lp2_lib.m"];
tag = "g" <> ToString[gg] <> "_c" <> ToString[cc] <> "of" <> ToString[nchunks];
log[m_] := Module[{s = OpenAppend["lp2_recon.log"]}, WriteString[s, DateString[{"Year", "-", "Month", "-", "Day", " ", "Hour", ":", "Minute", ":", "Second"}], "  [", tag, "] ", m, "\n"]; Close[s]; Print[m]];
loadOrbits[gg];
Get["lp2_chunk_" <> tag <> ".mx"];   (* chunkRecs, chunkIdx, tsOK, xsAll *)
log["start: " <> ToString[Length[chunkIdx]] <> " decorated graphs, " <> ToString[Length[tsOK]] <> " t-values x " <> ToString[Length[xsAll]] <> " eps-ratios"];
tcSetup[2 gg + 4];
entries = {}; failed = {}; termSum = 0; polarSum = ConstantArray[0, 2 gg - 2];
Do[Module[{i = chunkIdx[[pos]], id, samp, term, pol, tt = AbsoluteTime[]},
   id = OT["ids"][[i]];
   samp = Association[Table[t0 -> Association[Table[x0 -> chunkRecs[{t0, x0}][[pos]], {x0, xsAll}]], {t0, tsOK}]];
   term = reconstructOrbit[samp, tsOK, xsAll, gg];
   If[Head[term] === List && term[[1]] === $Failed, AppendTo[failed, {id, term[[2]]}]; log["RECONSTRUCTION FAILED for " <> ToString[id] <> ": " <> ToString[term[[2]]]]; Continue[]];
   pol = polarOfTerm[term, gg];
   termSum += term; polarSum += pol;
   AppendTo[entries, {id, OT["orbits"][[id[[1]], id[[2]]]], Normal[OT["descriptions"][id]], term, pol}];
   log["orbit " <> ToString[id] <> " " <> ToString[OT["orbits"][[id[[1]], id[[2]]]]] <> " done (" <> ToString[Round[AbsoluteTime[] - tt, 0.1]] <> " s), LeafCount " <> ToString[LeafCount[term]] <> ", polar " <> ToString[InputForm[pol]]]], {pos, Length[chunkIdx]}];
Put[entries, "lp2_terms_" <> tag <> ".m"];
Put[{Together[termSum], polarSum, failed}, "lp2_partial_" <> tag <> ".m"];
log["finished: " <> ToString[Length[entries]] <> " reconstructed, " <> ToString[Length[failed]] <> " failed"];
