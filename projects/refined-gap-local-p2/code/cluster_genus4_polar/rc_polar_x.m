(* rc_polar_x.m -- per-decorated-graph worker for the resolved conifold at alpha0 = 0, alpha1 = eps1 = x, alpha2 = eps2 = 1:
   exact contribution as a rational function of Q and x and its polar part at Q = 1 in t = Log Q.
   (A Laurent-series variant with x symbolic was tried and is far slower than the exact route; it was removed.)
   Usage:  wolframscript -file rc_polar_x.m genus chunk nChunks
   Input:  rc_factors_g<genus>_x.mx (from rc_prep_x.m), rc_orbits_g<genus>.mx (from rc_prep.m). Orbit i belongs to chunk Mod[i-1, n]+1.
   Output: rc_polar_x_g<genus>_c<chunk>of<n>.m = list of {id, key, description, polar {t^-(2g-2), ..., t^-1}, rational function or Null}
           rc_polar_x_partial_g<genus>_c<chunk>of<n>.m = {sum of polar parts, sum of rational functions (or Null), failed ids}; log in rc_polar_x.log *)
SetDirectory[DirectoryName[$InputFileName]];
args = Rest[$ScriptCommandLine];
{gg, cc, nchunks} = ToExpression /@ args[[1 ;; 3]];
Get["rc_lib.m"];
tag = "g" <> ToString[gg] <> "_c" <> ToString[cc] <> "of" <> ToString[nchunks];
log[m_] := Module[{s = OpenAppend["rc_polar_x.log"]}, WriteString[s, DateString[{"Year", "-", "Month", "-", "Day", " ", "Hour", ":", "Minute", ":", "Second"}], "  [", tag, "] ", m, "\n"]; Close[s]; Print[m]];
Get["rc_factors_g" <> ToString[gg] <> "_x.mx"];
loadOrbits[gg];
chunkIdx = Select[Range[Length[OT["ids"]]], Mod[# - 1, nchunks] + 1 == cc &];
log["start: " <> ToString[Length[chunkIdx]] <> " decorated graphs, ZMAX " <> ToString[ZMAX]];
entries = {}; failed = {}; polarSum = ConstantArray[0, 2 gg - 2]; rqSum = 0;
Do[Module[{id = OT["ids"][[i]], v, rq = Null, pol, tt = AbsoluteTime[]},
   v = graphValueC[OT["graphs"][id]];
   If[! cSymQ[v], AppendTo[failed, id]; log["NON-SYMMETRIC value for orbit " <> ToString[id]]; Continue[]];
   rq = Together[v[[1]]]; pol = polarOfRational[rq, gg]; rqSum += rq;
   pol = Factor /@ pol; polarSum += pol;
   AppendTo[entries, {id, OT["orbits"][[id[[1]], id[[2]]]], Normal[OT["descriptions"][id]], pol, rq}];
   log["orbit " <> ToString[id] <> " " <> ToString[OT["orbits"][[id[[1]], id[[2]]]]] <> " (" <> ToString[Round[AbsoluteTime[] - tt, 0.1]] <> " s): polar " <> ToString[InputForm[pol]]]], {i, chunkIdx}];
Put[entries, "rc_polar_x_" <> tag <> ".m"];
Put[{polarSum, Together[rqSum], failed}, "rc_polar_x_partial_" <> tag <> ".m"];
log["finished: " <> ToString[Length[entries]] <> " orbits, " <> ToString[Length[failed]] <> " failed"];
