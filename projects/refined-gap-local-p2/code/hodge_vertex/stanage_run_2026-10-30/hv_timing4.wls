(* hv_timing4.wls -- local timing of one genus-4 sample point (R-matrix rows to z^9, evalPoint, graph sums over the 27172 decorated
   graphs) and the per-point total check against F40. Usage: wolframscript -file hv_timing4.wls *)
SetDirectory[DirectoryName[$InputFileName]];
Get["hv_lib4.m"];
log[m_] := Module[{s = OpenAppend["hv_timing4.log"]}, WriteString[s, DateString[{"Year", "-", "Month", "-", "Day", " ", "Hour", ":", "Minute", ":", "Second"}], "  ", m, "\n"]; Close[s]; Print[m]];
gg = 4; KMAX = 9;
tt = AbsoluteTime[]; Get["hv_orbits_g4.mx"]; log["orbits loaded: " <> ToString[Length[OTH["ids"]]] <> " (" <> ToString[Round[AbsoluteTime[] - tt, 0.1]] <> " s)"];
tt = AbsoluteTime[]; loadRmatrix[KMAX]; log["buildRmatrix[9]: " <> ToString[Round[AbsoluteTime[] - tt, 0.1]] <> " s"];
wB = {e1s + e2s, -e1s, -e2s}; usubB = usubOf[wB];
tt = AbsoluteTime[]; rows = rowsTimes[hCoeffs[wB, -1, KMAX], KMAX]; log["rows of R': " <> ToString[Round[AbsoluteTime[] - tt, 0.1]] <> " s; LeafCounts " <> ToString[LeafCount /@ rows]];
tt = AbsoluteTime[]; Fref = Get["ref/genRmat_F40_cluster_result.m"][[1]] /. {q0 -> 0, q1 -> 0, q2 -> 0, q3 -> 0, q4 -> 0}; log["reference F40 loaded (" <> ToString[Round[AbsoluteTime[] - tt, 0.1]] <> " s)"];
{t0, x0} = {1/2, 3};
tt = AbsoluteTime[]; data = evalPointRows[t0, x0, KMAX, rows]; log["evalPoint at (1/2, 3): " <> ToString[Round[AbsoluteTime[] - tt, 0.1]] <> " s"];
tt = AbsoluteTime[]; r = computeSampleH[t0, x0, gg, OTH, rows, usubB /. {e1s -> x0, e2s -> 1}];
log["full point (evalPoint + 27172 graph sums): " <> ToString[Round[AbsoluteTime[] - tt, 0.1]] <> " s; diagnostics " <> ToString[r[[3]]]];
ref = refVector[Fref, t0, x0, gg];
log["sum of decorated graphs minus F40 at the point: " <> ToString[Total[r[[2]]] - ref]];
log[If[Union[Total[r[[2]]] - ref] === {0}, "VALIDATION OK: genus-4 orbit sum = F40 at one point", "VALIDATION FAIL: genus-4 orbit sum differs from F40"]];
