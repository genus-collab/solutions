(* rc_sample.m -- sampling worker for the resolved conifold on the general five-dimensional torus: at every integer weight
   sample assigned to this worker, the exact rational functions of Q of F11 (one leg), F20, ..., F<genus>0.
   Usage:  wolframscript -file rc_sample.m genus ZMAX worker nWorkers [method]
     genus 2, 3 or 4 (ZMAX >= 3 genus - 3); worker w of nWorkers takes the samples w, w + nWorkers, ... of the list
     wsamplesAll[[;; nSamplesFor[genus]]] (same generator and seed as resolvedconifold_genus3.nb);
     method "closed" (default: closed-form QRR constants) or "limit" (constants from the Q -> 0 limit, as in the notebook).
   Output: rc_samples_g<genus>_w<w>of<n>.m, one record {i, wsub, {{N11, a11}, {N20, a20}, ...}} per sample (appended; a
   restarted job skips the samples already present); log lines in rc_sample.log. *)
SetDirectory[DirectoryName[$InputFileName]];
args = Rest[$ScriptCommandLine];
If[Length[args] < 4, Print["usage: wolframscript -file rc_sample.m genus ZMAX worker nWorkers [method]"]; Exit[1]];
{gg, zmax, w, nw} = ToExpression /@ args[[1 ;; 4]];
method = If[Length[args] >= 5, args[[5]], "closed"];
Get["rc_lib.m"];
tag = "g" <> ToString[gg] <> "_w" <> ToString[w] <> "of" <> ToString[nw];
log[m_] := Module[{s = OpenAppend["rc_sample.log"]}, WriteString[s, DateString[{"Year", "-", "Month", "-", "Day", " ", "Hour", ":", "Minute", ":", "Second"}], "  [", tag, "] ", m, "\n"]; Close[s]; Print[m]];
log["start: genus " <> ToString[gg] <> ", ZMAX " <> ToString[zmax] <> ", method " <> method <> ", " <> $Version];
setupWick[zmax];
preps = Table[prepareGraph /@ graphsWithLegs[gi, If[gi == 1, 1, 0]], {gi, 1, gg}];
log["graphs: " <> ToString[Table[{gi, Length[preps[[gi]]], Total[Length[#["decs"]] & /@ preps[[gi]]]}, {gi, 1, gg}]]];
ws = wsamplesAll[[;; nSamplesFor[gg]]]; qs = qsFor[gg];
outfile = "rc_samples_" <> tag <> ".m";
mine = Range[w, Length[ws], nw];
done = If[FileExistsQ[outfile], Quiet[Check[ReadList[outfile][[All, 1]], {}]], {}];
todo = Select[mine, ! MemberQ[done, #] &];
log[ToString[Length[mine]] <> " samples assigned (of " <> ToString[Length[ws]] <> "), " <> ToString[Length[done]] <> " already done, " <> ToString[Length[todo]] <> " to do; " <> ToString[Length[qs]] <> " Q-values"];
Do[Module[{wsub = wsubOf[ws[[i]]], t0 = AbsoluteTime[], tR, vals, recs, times = {}, sym = True},
   If[method == "limit", rmatrixAt[wsub], rmatrixAtClosed[wsub]];
   tR = Round[AbsoluteTime[] - t0, 0.1];
   vals = Table[evalQ[q];
     Prepend[Table[Module[{t1 = AbsoluteTime[], v = graphSumC[preps[[gi]]]}, AppendTo[times, {gi, Round[AbsoluteTime[] - t1, 0.1]}]; If[! cSymQ[v], sym = False]; v[[1]]], {gi, 1, gg}], q], {q, qs}];
   If[! sym, log["sample " <> ToString[i] <> " " <> ToString[ws[[i]]] <> ": NON-SYMMETRIC graph sum, skipped"]; Continue[]];
   recs = Table[reconQ[Table[{v[[1]], v[[gi + 1]]}, {v, vals}]], {gi, 1, gg}];
   If[MemberQ[recs, $Failed], log["sample " <> ToString[i] <> " " <> ToString[ws[[i]]] <> ": Q-reconstruction FAILED for genus " <> ToString[Flatten[Position[recs, $Failed]]] <> ", skipped"]; Continue[]];
   PutAppend[{i, wsub, recs}, outfile];
   log["sample " <> ToString[i] <> " " <> ToString[ws[[i]]] <> ": pole orders " <> ToString[recs[[All, 2]]] <> ", numerator degrees " <> ToString[Exponent[#[[1]], Q] & /@ recs] <> "  (R-matrix " <> ToString[tR] <> " s, graph sums per genus at the last Q " <> ToString[Table[Last[Select[times, #[[1]] == gi &]][[2]], {gi, 1, gg}]] <> " s, total " <> ToString[Round[AbsoluteTime[] - t0, 0.1]] <> " s)"]], {i, todo}];
log["job finished"];
