(* hv_genus2.wls -- F_{2,0} of K_P2 x C^2 (fibrewise torus) from three copies of the Hodge CohFT of C^3 with the fibre
   weights w = {-eps1-eps2, eps1, eps2} (mode A) or the opposite sign (mode B), compared with the point-vertex sum (mode point).
   Usage:
     wolframscript -file hv_genus2.wls checks                  Hodge-integral consistency checks and the PPZ scalar R-action check
     wolframscript -file hv_genus2.wls calib                   determine the sign s in R' = R_Z h^s from the total at two sample points
     wolframscript -file hv_genus2.wls sample <mode> [w nw]    sampling worker (mode point | A | B), points w, w+nw, ... of the list
     wolframscript -file hv_genus2.wls recon                   reconstruction, per-graph terms, checks, output files
   Output: hv_samples_<mode>_p<idx>.mx, hv_calib.m, hodge_vertex_genus2_results.m, hodge_vertex_genus2.log *)
SetDirectory[DirectoryName[$InputFileName]];
args = Rest[$ScriptCommandLine];
stage = If[Length[args] >= 1, args[[1]], "checks"];
Get["hv_lib.m"];
gg = 2; KMAX = kmaxFor[gg]; NX = 34;
log[m_] := Module[{s = OpenAppend["hodge_vertex_genus2.log"]}, WriteString[s, DateString[{"Year", "-", "Month", "-", "Day", " ", "Hour", ":", "Minute", ":", "Second"}], "  [", stage, "] ", m, "\n"]; Close[s]; Print[m]];
wA = {-(e1s + e2s), e1s, e2s}; wB = -wA;
weightsOf[mode_] := Switch[mode, "A", wA, "B", wB, _, None];
Fref = Get["ref/F20_F30_y.m"][[1]];
ts = tsamplesFor[gg]; xs = xsamplesH[NX];
allpts = Flatten[Table[{t0, x0}, {t0, ts}, {x0, xs}], 1];
pfile[mode_, idx_] := "hv_samples_" <> mode <> "_p" <> ToString[idx] <> ".mx";

setupTables[] := (
   prepared = prepareGraph /@ graphsWithLegs[gg, 0]; OT = orbitTable[prepared];
   preparedH = prepareGraphH /@ graphsWithLegs[gg, 0]; OTH = orbitTable[preparedH];
   log["point vertices: " <> ToString[Length[prepared]] <> " graphs, " <> ToString[Total[Length[#["decs"]] & /@ prepared]] <> " labelled decorations, " <> ToString[Length[OT["ids"]]] <> " orbits; Hodge vertices: " <>
     ToString[Total[Length[#["decs"]] & /@ preparedH]] <> " labelled decorations, " <> ToString[Length[OTH["ids"]]] <> " orbits"]);
modeData[mode_] := Module[{w = weightsOf[mode], s},   (* -> {OTx, rows, hesubSymbolic} *)
   If[mode == "point", Return[{OT, None, {}}]];
   s = If[FileExistsQ["hv_calib.m"], Get["hv_calib.m"]["s"], Print["no hv_calib.m -- run calib first"]; Exit[1]];
   {OTH, rowsTimes[hCoeffs[w, s, KMAX], KMAX], heSub[w]}];
numsub[hesub_, x0_] := hesub /. {e1s -> x0, e2s -> 1};

Which[
 stage == "checks",
 Module[{ok, rep, w3 = {w1, w2, w3}, ppz, expct},
  {ok, rep} = hodgeChecks[]; Do[log[r], {r, rep}];
  ppz = ppzScalarG2[w3]; expct = hodgeNorm[2, {}] /. heSub[w3];
  log["PPZ scalar R-action of h on the point theory, genus 2, general weights: " <> ToString[InputForm[Factor[ppz]]] <> "; normalised Hodge class on M_2: " <> ToString[InputForm[Factor[expct]]] <> "; equal: " <> ToString[Together[ppz - expct] === 0]];
  log["normalised Hodge class on M_2 for the CY weights wA: " <> ToString[InputForm[Factor[hodgeNorm[2, {}] /. heSub[wA]]]] <> " (geometric vertex prod Lambda^vee(w)/w integrates to 1/5760 for CY weights)"];
  log["h(z) for wA to z^3: " <> ToString[InputForm[Factor /@ hCoeffs[wA, 1, 3]]]];
  log[If[ok && Together[ppz - expct] === 0, "CHECKS OK", "CHECKS FAILED"]]],

 stage == "calib",
 Module[{pts = {{1/2, 3}, {2/3, 5}}, rowsP, rowsM, res},
  setupTables[]; loadRmatrix[KMAX];
  rowsP = rowsTimes[hCoeffs[wA, 1, KMAX], KMAX]; rowsM = rowsTimes[hCoeffs[wA, -1, KMAX], KMAX];
  res = Table[Module[{t0 = p[[1]], x0 = p[[2]], ref, tp, tA1, tA2},
     ref = refVector[Fref, t0, x0, gg];
     tp = Total[computeSampleH[t0, x0, gg, OT, None, {}][[2]]];
     tA1 = Total[computeSampleH[t0, x0, gg, OTH, rowsP, numsub[heSub[wA], x0]][[2]]];
     tA2 = Total[computeSampleH[t0, x0, gg, OTH, rowsM, numsub[heSub[wA], x0]][[2]]];
     log["point " <> ToString[p] <> ": point vertices - ref = " <> ToString[tp - ref] <> "; Hodge(wA) with R_Z h^(+1) - ref = " <> ToString[N[tA1 - ref, 6]] <> "; with R_Z h^(-1) - ref = " <> ToString[N[tA2 - ref, 6]]];
     {tp - ref, tA1 - ref, tA2 - ref}], {p, pts}];
  Which[
   Union[Flatten[res[[All, 3]]]] === {0} && Union[Flatten[res[[All, 1]]]] === {0}, s = -1,
   Union[Flatten[res[[All, 2]]]] === {0} && Union[Flatten[res[[All, 1]]]] === {0}, s = 1,
   True, log["CALIBRATION FAILED: neither sign reproduces F20"]; Exit[1]];
  Put[<|"s" -> s, "weightsA" -> wA, "weightsB" -> wB, "note" -> "R' = R_Z h^s with h the Bernoulli series of the vertex weights (mode A: wA, mode B: wB uses h_wB = h_wA^(-1), same s)"|>, "hv_calib.m"];
  log["CALIBRATION OK: s = " <> ToString[s] <> " (rows of R' = rows of R_Z times h^s); written hv_calib.m"]],

 stage == "sample",
 Module[{mode = args[[2]], w = 1, nw = 1, OTx, rows, hesub, mine, todo},
  If[Length[args] >= 4, {w, nw} = ToExpression /@ args[[3 ;; 4]]];
  setupTables[]; loadRmatrix[KMAX];
  {OTx, rows, hesub} = modeData[mode];
  mine = Range[w, Length[allpts], nw]; todo = Select[mine, ! FileExistsQ[pfile[mode, #]] &];
  log["mode " <> mode <> ": " <> ToString[Length[mine]] <> " points assigned, " <> ToString[Length[todo]] <> " to do"];
  Do[Module[{p = allpts[[idx]], r, ref},
    r = computeSampleH[p[[1]], p[[2]], gg, OTx, rows, numsub[hesub, p[[2]]]];
    ref = refVector[Fref, p[[1]], p[[2]], gg];
    rec = Append[r, "total-ref" -> (Total[r[[2]]] - ref)];
    DumpSave[pfile[mode, idx], {rec}];
    log["mode " <> mode <> " point " <> ToString[idx] <> " " <> ToString[p] <> " in " <> ToString["seconds" /. r[[3]]] <> " s; total-ref " <> ToString[Total[r[[2]]] - ref] <> If[("bad" /. r[[3]]) =!= {}, "; WARNING unexpected monomials " <> ToString["bad" /. r[[3]]], ""]]], {idx, todo}];
  log["mode " <> mode <> " finished"]],

 stage == "recon",
 Module[{modes = {"point", "A", "B"}, results = <||>, samp, recs, terms, graphTerms, total, polar, ok, cal},
  setupTables[]; cal = Get["hv_calib.m"]; tcSetup[8];
  Do[Module[{OTx = If[mode == "point", OT, OTH], files, missing},
    files = Table[pfile[mode, idx], {idx, Length[allpts]}]; missing = Select[files, ! FileExistsQ[#] &];
    If[missing =!= {}, log["mode " <> mode <> ": " <> ToString[Length[missing]] <> " sample files missing -- skipped"]; Continue[]];
    recs = Table[Get[f]; rec, {f, files}];
    ok = Union[Flatten[("total-ref" /. #[[4]]) & /@ recs]] === {0};
    log["mode " <> mode <> ": " <> ToString[Length[recs]] <> " points; sum of orbits = F20 at every point: " <> ToString[ok]];
    Clear[samp]; Do[samp[r[[1, 1]]][r[[1, 2]]] = r[[2]], {r, recs}];
    terms = Association[Table[Module[{sv, o = oo}, Clear[sv]; Do[sv[t0][x0] = samp[t0][x0][[o]], {t0, ts}, {x0, xs}];
        OTx["ids"][[o]] -> reconstructOrbitH[sv, ts, xs, gg]], {oo, Length[OTx["ids"]]}]];
    If[Select[Values[terms], ListQ] =!= {}, log["mode " <> mode <> ": FAILED reconstructions " <> ToString[Select[Keys[terms], ListQ[terms[#]] &]]]];
    graphTerms = Table[Together[Total[Table[terms[id], {id, Select[OTx["ids"], #[[1]] == k &]}]]], {k, Length[prepared]}];
    total = Together[Total[graphTerms]];
    log["mode " <> mode <> ": total - F20 = " <> ToString[Together[total - Fref] === 0]];
    polar = Table[polarOfTerm[gt, gg], {gt, graphTerms}];
    log["mode " <> mode <> ": per-graph polar parts {t_c^-2, t_c^-1} (eps1 = x, eps2 = 1): " <> ToString[InputForm[polar]]];
    log["mode " <> mode <> ": total polar part " <> ToString[InputForm[Factor /@ Total[polar]]] <> "; gap coefficient " <> ToString[InputForm[Factor[gapCoefficient[gg] /. {e1s -> x, e2s -> 1}]]]];
    results[mode] = <|"graphTerms" -> graphTerms, "orbitTerms" -> terms, "polar" -> polar, "total" -> total, "totalOK" -> (Together[total - Fref] === 0),
      "descriptions" -> OTx["descriptions"], "weights" -> weightsOf[mode]|>], {mode, modes}];
  results["graphs"] = Table[<|"genera" -> prepared[[k]]["g"], "adjacency" -> prepared[[k]]["A"], "edges" -> prepared[[k]]["edges"], "Aut" -> prepared[[k]]["aut"]|>, {k, Length[prepared]}];
  results["calibration"] = cal;
  results["hCoeffsA"] = hCoeffs[wA, 1, KMAX];
  results["readme"] = "F_{2,0} of K_P2 x C^2 (fibrewise torus) by the Givental-Teleman graph sum with three kinds of vertices: 'point' (Witten-Kontsevich, R-matrix R_Z as in localp2CY5_genus3.nb), 'A' (Hodge CohFT of C^3 with weights {-eps1-eps2, eps1, eps2}, R' = R_Z h^s), 'B' (weights {eps1+eps2, -eps1, -eps2}). graphTerms: the 7 stable graphs of genus 2 in the order of 'graphs' (StableGraph-2.bin), as functions of y, XX = y c1'/c1, eps1, eps2; orbitTerms: per decorated graph (orbit key as in orbits.m); polar: {t_c^-2, t_c^-1} with eps1 = x, eps2 = 1. All three totals equal F20.";
  Put[results, "hodge_vertex_genus2_results.m"];
  log["written hodge_vertex_genus2_results.m"]],

 True, Print["unknown stage ", stage]];
