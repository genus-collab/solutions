(* Sampling job for the twisted Givental-Teleman graph sum of K_P2 x C^2 on the fibrewise torus.
   Usage:  wolframscript -file genRmat_g4_sample.m genus KMAX jobIndex nJobs [nKernels] [PREC]
     genus     2, 3 or 4  (KMAX must be >= 3 genus - 3: 3, 6, 9)
     jobIndex  1 .. nJobs; job j computes the sample points j, j + nJobs, j + 2 nJobs, ... of pointList[genus]
     nKernels  optional: > 1 uses that many local subkernels (ParallelMap) for the points of this job
     PREC      optional working precision in digits (default 350)
   Output: samples_g<genus>_job<jobIndex>of<nJobs>.m in the script directory, one expression
     {{t0, x0}, matrix, diagnostics} per finished sample point (appended after every point; the job can be
     restarted and skips the points already in the file), and a log line per point in genRmat_g4_sample.log.
   Needs in the same directory: genRmat_g4_lib.m, cPP_z<KMAX>.m, StableGraph-<genus>.bin. *)
SetDirectory[DirectoryName[$InputFileName]];
args = Rest[$ScriptCommandLine];
If[Length[args] < 4, Print["usage: wolframscript -file genRmat_g4_sample.m genus KMAX jobIndex nJobs [nKernels] [PREC]"]; Exit[1]];
{gg, KMAX, job, njobs} = ToExpression /@ args[[1 ;; 4]];
nk = If[Length[args] >= 5, ToExpression[args[[5]]], 1];
Get["genRmat_g4_lib.m"];
If[Length[args] >= 6, PREC = ToExpression[args[[6]]]];
tag = "g" <> ToString[gg] <> "_job" <> ToString[job] <> "of" <> ToString[njobs];
log[m_] := Module[{s = OpenAppend["genRmat_g4_sample.log"]}, WriteString[s, DateString[{"Year", "-", "Month", "-", "Day", " ", "Hour", ":", "Minute", ":", "Second"}], "  [", tag, "] ", m, "\n"]; Close[s]; Print[m]];
log["start: genus " <> ToString[gg] <> ", KMAX " <> ToString[KMAX] <> ", PREC " <> ToString[PREC] <> ", kernels " <> ToString[nk] <> ", " <> $Version];
setupGenus[gg, KMAX];
outfile = "samples_" <> tag <> ".m";
mine = pointList[gg][[job ;; ;; njobs]];
done = If[FileExistsQ[outfile], Quiet[Check[ReadList[outfile][[All, 1]], {}]], {}];
todo = Select[mine, ! MemberQ[done, #] &];
log[ToString[Length[mine]] <> " points assigned, " <> ToString[Length[done]] <> " already done, " <> ToString[Length[todo]] <> " to do"];
save[res_] := (PutAppend[res, outfile]; log["point " <> ToString[res[[1]]] <> " done in " <> ToString["seconds" /. res[[3]]] <> " s; max |Im| " <> ToString["imag" /. res[[3]]] <> If[("bad" /. res[[3]]) =!= {}, "; WARNING unexpected monomials " <> ToString["bad" /. res[[3]]], ""]]);
If[nk > 1,
  LaunchKernels[nk]; DistributeDefinitions["Global`"];
  log["launched " <> ToString[Length[Kernels[]]] <> " subkernels"];
  Do[save /@ ParallelMap[computeSample[#[[1]], #[[2]]] &, batch], {batch, Partition[todo, UpTo[Length[Kernels[]]]]}];
  CloseKernels[],
  Do[save[computeSample[p[[1]], p[[2]]]], {p, todo}]];
log["job finished"];
