(* hv_sample.m -- sampling worker: exact XX-coefficient vectors of every decorated graph (Hodge vertices, weights
   w = {eps1+eps2, -eps1, -eps2}, R' = R_Z h^(-1)) at the sample points (t0 = y^(1/3), x0 = eps1/eps2) assigned to this worker.
   Usage:  wolframscript -file hv_sample.m genus KMAX worker nWorkers [PREC]
     genus 2, 3 or 4 (KMAX >= 3 genus - 3); worker w of nWorkers computes the points w, w + nWorkers, ... of pointListH[genus].
   Output: hv_samples_g<genus>_p<idx>.mx holding `rec` = {{t0, x0}, vectors (order OTH["ids"]), diagnostics, "total-ref" -> vector}
           (one file per point; a restarted job skips the points already present); one log line per point in hv_sample.log.
   Needs: hv_lib4.m, lp2_lib.m, graphs.m, orbits.m, cPP_z<K>.m, StableGraph-<genus>.bin, hv_orbits_g<genus>.mx (from hv_prep.m),
          ref/F20_F30_y.m (genus 2, 3) or ref/genRmat_F40_cluster_result.m (genus 4: q -> 0 part) for the per-point total check. *)
SetDirectory[DirectoryName[$InputFileName]];
args = Rest[$ScriptCommandLine];
If[Length[args] < 4, Print["usage: wolframscript -file hv_sample.m genus KMAX worker nWorkers [PREC]"]; Exit[1]];
{gg, KMAX, w, nw} = ToExpression /@ args[[1 ;; 4]];
Get["hv_lib4.m"];
If[Length[args] >= 5, PREC = ToExpression[args[[5]]]];
tag = "g" <> ToString[gg] <> "_w" <> ToString[w] <> "of" <> ToString[nw];
log[m_] := Module[{s = OpenAppend["hv_sample.log"]}, WriteString[s, DateString[{"Year", "-", "Month", "-", "Day", " ", "Hour", ":", "Minute", ":", "Second"}], "  [", tag, "] ", m, "\n"]; Close[s]; Print[m]];
log["start: genus " <> ToString[gg] <> ", KMAX " <> ToString[KMAX] <> ", PREC " <> ToString[PREC] <> ", " <> $Version];
Get["hv_orbits_g" <> ToString[gg] <> ".mx"]; loadRmatrix[KMAX];
wB = {e1s + e2s, -e1s, -e2s}; usubB = usubOf[wB];
rows = rowsTimes[hCoeffs[wB, -1, KMAX], KMAX];
log["rows of R' = R_Z h^(-1) built to z^" <> ToString[KMAX] <> "; LeafCounts " <> ToString[LeafCount /@ rows]];
Fref = If[gg <= 3, Get["ref/F20_F30_y.m"][[gg - 1]], Get["ref/genRmat_F40_cluster_result.m"][[1]] /. {q0 -> 0, q1 -> 0, q2 -> 0, q3 -> 0, q4 -> 0}];
allpts = pointListH[gg];
mine = Range[w, Length[allpts], nw];
pfile[idx_] := "hv_samples_g" <> ToString[gg] <> "_p" <> ToString[idx] <> ".mx";
todo = Select[mine, ! FileExistsQ[pfile[#]] &];
log[ToString[Length[mine]] <> " points assigned, " <> ToString[Length[mine] - Length[todo]] <> " already done, " <> ToString[Length[todo]] <> " to do"];
Do[Module[{p = allpts[[idx]], r, ref},
   r = computeSampleH[p[[1]], p[[2]], gg, OTH, rows, usubB /. {e1s -> p[[2]], e2s -> 1}];
   ref = refVector[Fref, p[[1]], p[[2]], gg];
   rec = Append[r, "total-ref" -> (Total[r[[2]]] - ref)];
   DumpSave[pfile[idx], {rec}];
   log["point " <> ToString[idx] <> " " <> ToString[p] <> " done in " <> ToString["seconds" /. r[[3]]] <> " s; total-ref " <> ToString[Total[r[[2]]] - ref] <> "; max |Im| " <> ToString["imag" /. r[[3]]] <> If[("bad" /. r[[3]]) =!= {}, "; WARNING unexpected monomials in orbits " <> ToString["bad" /. r[[3]]], ""]]], {idx, todo}];
log["job finished"];
