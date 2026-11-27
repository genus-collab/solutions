(* lp2_sample.m -- sampling worker: exact XX-coefficient vectors of every decorated stable graph of local P2 x C2 (T_F)
   at the sample points (t0 = y^(1/3), x0 = eps1/eps2) assigned to this worker.
   Usage:  wolframscript -file lp2_sample.m genus KMAX worker nWorkers [PREC]
     genus 2, 3 or 4 (KMAX >= 3 genus - 3); worker w of nWorkers computes the points w, w + nWorkers, ... of pointList[genus].
   Output: lp2_samples_g<genus>_w<w>of<n>.mx holding `recs` = list of {{t0, x0}, vectors, diagnostics} (re-saved after every
   point; a restarted job skips the points already present); one log line per point in lp2_sample.log.
   Needs: lp2_lib.m, graphs.m, orbits.m, cPP_z<K>.m, StableGraph-<genus>.bin, lp2_orbits_g<genus>.mx (from lp2_prep.m). *)
SetDirectory[DirectoryName[$InputFileName]];
args = Rest[$ScriptCommandLine];
If[Length[args] < 4, Print["usage: wolframscript -file lp2_sample.m genus KMAX worker nWorkers [PREC]"]; Exit[1]];
{gg, KMAX, w, nw} = ToExpression /@ args[[1 ;; 4]];
Get["lp2_lib.m"];
If[Length[args] >= 5, PREC = ToExpression[args[[5]]]];
tag = "g" <> ToString[gg] <> "_w" <> ToString[w] <> "of" <> ToString[nw];
log[m_] := Module[{s = OpenAppend["lp2_sample.log"]}, WriteString[s, DateString[{"Year", "-", "Month", "-", "Day", " ", "Hour", ":", "Minute", ":", "Second"}], "  [", tag, "] ", m, "\n"]; Close[s]; Print[m]];
log["start: genus " <> ToString[gg] <> ", KMAX " <> ToString[KMAX] <> ", PREC " <> ToString[PREC] <> ", " <> $Version];
loadOrbits[gg]; loadRmatrix[KMAX];
(* one DumpSave file per finished point (no growing re-saves): lp2_samples_<tag>_p<index in pointList>.mx holding `rec` *)
allpts = pointList[gg];
mine = Range[w, Length[allpts], nw];
pfile[idx_] := "lp2_samples_" <> tag <> "_p" <> ToString[idx] <> ".mx";
todo = Select[mine, ! FileExistsQ[pfile[#]] &];
log[ToString[Length[mine]] <> " points assigned, " <> ToString[Length[mine] - Length[todo]] <> " already done, " <> ToString[Length[todo]] <> " to do"];
Do[Module[{p = allpts[[idx]], r},
   r = computeSampleOrbits[p[[1]], p[[2]], gg];
   rec = r; DumpSave[pfile[idx], {rec}];
   log["point " <> ToString[p] <> " done in " <> ToString["seconds" /. r[[3]]] <> " s; max |Im| " <> ToString["imag" /. r[[3]]] <> If[("bad" /. r[[3]]) =!= {}, "; WARNING unexpected monomials in orbits " <> ToString["bad" /. r[[3]]], ""]]], {idx, todo}];
log["job finished"];
