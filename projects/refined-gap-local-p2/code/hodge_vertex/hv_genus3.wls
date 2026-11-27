(* hv_genus3.wls -- F_{3,0} of K_P2 x C^2 (fibrewise torus) from three copies of the Hodge CohFT of C^3 with the e_T weights of
   the normal bundle, w = {eps1+eps2, -eps1, -eps2} (mode B), R' = R_Z h_w^(-1); contribution of every stable / decorated graph.
   Usage:
     wolframscript -file hv_genus3.wls checks            Hodge-integral checks, PPZ scalar checks in genus 2 (symbolic) and 3 (exact numeric)
     wolframscript -file hv_genus3.wls prep              decorated graphs with Hodge vertices -> hv_orbits_g3.mx
     wolframscript -file hv_genus3.wls sample [w nw]     sampling worker (points w, w+nw, ...) -> hv3_samples_B_p<idx>.mx
     wolframscript -file hv_genus3.wls recon             reconstruction, per-graph terms, F30 and gap checks, polar parts, output files
   Output: hodge_vertex_genus3_results.m, hodge_vertex_genus3_summary.txt, hodge_vertex_genus3.log *)
SetDirectory[DirectoryName[$InputFileName]];
args = Rest[$ScriptCommandLine];
stage = If[Length[args] >= 1, args[[1]], "checks"];
Get["hv_lib3.m"];
gg = 3; KMAX = kmaxFor[gg]; NX = 30;
log[m_] := Module[{s = OpenAppend["hodge_vertex_genus3.log"]}, WriteString[s, DateString[{"Year", "-", "Month", "-", "Day", " ", "Hour", ":", "Minute", ":", "Second"}], "  [", stage, "] ", m, "\n"]; Close[s]; Print[m]];
wB = {e1s + e2s, -e1s, -e2s};
usubB = usubOf[wB];
Fref = Get["ref/F20_F30_y.m"][[2]];
ts = tsamplesFor[gg]; xs = xsamplesH[NX];
allpts = Flatten[Table[{t0, x0}, {t0, ts}, {x0, xs}], 1];
pfile[idx_] := "hv3_samples_B_p" <> ToString[idx] <> ".mx";
fmt[expr_] := ToString[InputForm[expr /. {e1s -> E1, e2s -> E2}]];

Which[
 stage == "checks",
 Module[{ok, rep, w3 = {w1, w2, w3}, ppz2, exp2, ppz3, exp3, wnum, tt},
  tt = AbsoluteTime[];
  {ok, rep} = hodgeChecks3[]; Do[log[r], {r, rep}];
  ppz2 = ppzScalar[2, w3]; exp2 = hodgeNormU[2, {}] /. usubOf[w3];
  log["PPZ scalar R-action, genus 2, general weights: " <> ToString[InputForm[Factor[ppz2]]] <> " = normalised Hodge class on M_2: " <> ToString[Together[ppz2 - exp2] === 0]];
  SeedRandom[7];
  Do[wnum = Table[RandomInteger[{-9, 9}]/RandomInteger[{1, 5}], {3}]; If[MemberQ[wnum, 0], Continue[]];
   ppz3 = ppzScalar[3, wnum]; exp3 = hodgeNormU[3, {}] /. usubOf[wnum];
   log["PPZ scalar R-action, genus 3 (42 graphs), weights " <> ToString[wnum] <> ": " <> ToString[ppz3] <> " = normalised Hodge class on M_3 " <> ToString[exp3] <> ": " <> ToString[ppz3 == exp3]]; ok = ok && ppz3 == exp3, {2}];
  log["h_wB(z) to z^6: " <> fmt[Factor /@ hCoeffs[wB, 1, 6]]];
  log["normalised genus-3 Hodge class for wB (coefficients of lambda monomials {i,j,k}): " <> fmt[Table[{c[[1]], Factor[c[[2]] /. usubB]}, {c, hodgeClass[3]}]]];
  log[If[ok, "CHECKS OK (" <> ToString[Round[AbsoluteTime[] - tt]] <> " s)", "CHECKS FAILED"]]],

 stage == "prep",
 Module[{tt = AbsoluteTime[]},
  preparedH3 = prepareGraphH /@ graphsWithLegs[gg, 0]; OTH3 = orbitTable[preparedH3];
  hodgeClass3 = hodgeClass[3];
  DumpSave["hv_orbits_g3.mx", {preparedH3, OTH3, hodgeClass3}];
  log["genus 3 with Hodge vertices: " <> ToString[Length[preparedH3]] <> " graphs, " <> ToString[Total[Length[#["decs"]] & /@ preparedH3]] <> " labelled decorations, " <> ToString[Length[OTH3["ids"]]] <> " decorated graphs (orbits); " <> ToString[Length[allpts]] <> " sample points (" <> ToString[Length[ts]] <> " t x " <> ToString[NX] <> " eps-ratios); " <> ToString[Round[AbsoluteTime[] - tt]] <> " s; written hv_orbits_g3.mx"]],

 stage == "sample",
 Module[{w = 1, nw = 1, rows, mine, todo, hc},
  If[Length[args] >= 3, {w, nw} = ToExpression /@ args[[2 ;; 3]]];
  Get["hv_orbits_g3.mx"]; loadRmatrix[KMAX];
  hc = hCoeffs[wB, -1, KMAX]; rows = rowsTimes[hc, KMAX];
  log["worker " <> ToString[w] <> "/" <> ToString[nw] <> ": rows of R' = R_Z h^(-1) built to z^" <> ToString[KMAX] <> "; " <> $Version];
  mine = Range[w, Length[allpts], nw]; todo = Select[mine, ! FileExistsQ[pfile[#]] &];
  log["worker " <> ToString[w] <> ": " <> ToString[Length[mine]] <> " points assigned, " <> ToString[Length[todo]] <> " to do"];
  Do[Module[{p = allpts[[idx]], r, ref},
    r = computeSampleH[p[[1]], p[[2]], gg, OTH3, rows, usubB /. {e1s -> p[[2]], e2s -> 1}];
    ref = refVector[Fref, p[[1]], p[[2]], gg];
    rec = Append[r, "total-ref" -> (Total[r[[2]]] - ref)];
    DumpSave[pfile[idx], {rec}];
    log["worker " <> ToString[w] <> " point " <> ToString[idx] <> " " <> ToString[p] <> " in " <> ToString["seconds" /. r[[3]]] <> " s; total-ref " <> ToString[Total[r[[2]]] - ref] <> If[("bad" /. r[[3]]) =!= {}, "; WARNING unexpected monomials " <> ToString["bad" /. r[[3]]], ""]]], {idx, todo}];
  log["worker " <> ToString[w] <> " finished"]],

 stage == "recon",
 Module[{files, missing, recs, ok, samp, terms, failed, graphTerms, total, polar, polarOrbits, numchk, symchk, gapok, expected, results, out = {}, emit, tt = AbsoluteTime[]},
  Get["hv_orbits_g3.mx"]; tcSetup[10];
  files = Table[pfile[idx], {idx, Length[allpts]}]; missing = Select[files, ! FileExistsQ[#] &];
  If[missing =!= {}, log[ToString[Length[missing]] <> " sample files missing -- abort"]; Exit[1]];
  recs = Table[Get[f]; rec, {f, files}];
  ok = Union[Flatten[("total-ref" /. #[[4]]) & /@ recs]] === {0};
  log[ToString[Length[recs]] <> " points; sum of decorated graphs = F30 at every point: " <> ToString[ok]];
  Clear[samp]; Do[samp[r[[1, 1]]][r[[1, 2]]] = r[[2]], {r, recs}];
  terms = Association[Table[Module[{sv, o = oo}, Clear[sv]; Do[sv[t0][x0] = samp[t0][x0][[o]], {t0, ts}, {x0, xs}];
      OTH3["ids"][[o]] -> reconstructOrbitH[sv, ts, xs, gg]], {oo, Length[OTH3["ids"]]}]];
  failed = Select[Keys[terms], ListQ[terms[#]] &];
  log["reconstruction: " <> ToString[Length[terms] - Length[failed]] <> " of " <> ToString[Length[terms]] <> " decorated graphs; failed: " <> ToString[failed] <> "; " <> ToString[Round[AbsoluteTime[] - tt]] <> " s"];
  If[failed =!= {}, Exit[1]];
  graphTerms = Table[Together[Total[Table[terms[id], {id, Select[OTH3["ids"], #[[1]] == k &]}]]], {k, Length[preparedH3]}];
  total = Together[Total[graphTerms]];
  SeedRandom[4];
  numchk = Table[Together[(total - Fref) /. {y -> RandomInteger[{1, 50}]/RandomInteger[{51, 97}], XX -> RandomInteger[{-30, 30}]/RandomInteger[{1, 17}], e1s -> RandomInteger[{1, 40}], e2s -> RandomInteger[{1, 40}]}], {5}];
  symchk = TimeConstrained[Together[total - Fref] === 0, 1800, "timeout"];
  log["total minus F30 at 5 random rational points: " <> ToString[numchk] <> "; symbolic: " <> ToString[symchk]];
  log[If[Union[numchk] === {0}, "VALIDATION OK: sum of the decorated graphs = F30", "VALIDATION FAIL: total differs from F30"]];
  polar = Table[polarOfTerm[gt, gg] /. SeriesCoefficient[a_, _] /; FreeQ[a, tcv] :> 0, {gt, graphTerms}];
  expected = gapCoefficient[gg] /. {e1s -> x, e2s -> 1};
  gapok = Together[Total[polar][[1]] - expected] === 0 && Union[Together /@ Rest[Total[polar]]] === {0};
  log["total polar part {t_c^-4, ..., t_c^-1} (eps1 = x, eps2 = 1): " <> ToString[InputForm[Factor /@ Total[polar]]] <> "; expected leading " <> ToString[InputForm[Factor[expected]]] <> "; gap: " <> ToString[gapok]];
  polarOrbits = Association[Table[id -> (polarOfTerm[terms[id], gg] /. SeriesCoefficient[a_, _] /; FreeQ[a, tcv] :> 0), {id, Keys[terms]}]];
  log["per stable graph: polar parts: " <> ToString[InputForm[Table[{k, preparedH3[[k]]["g"], preparedH3[[k]]["edges"], Factor /@ polar[[k]]}, {k, Length[preparedH3]}]]]];
  log["decorated graphs with a nonzero polar part: " <> ToString[Select[Keys[polarOrbits], Union[polarOrbits[#]] =!= {0} &]]];
  results = <|"graphTerms" -> graphTerms, "orbitTerms" -> terms, "polar" -> polar, "polarOrbits" -> polarOrbits, "total" -> total, "totalOK" -> (Union[numchk] === {0}), "symbolicOK" -> symchk, "gapOK" -> gapok,
    "descriptions" -> OTH3["descriptions"], "weights" -> wB, "hCoeffs" -> hCoeffs[wB, 1, KMAX], "hodgeClass3" -> hodgeClass3,
    "graphs" -> Table[<|"genera" -> preparedH3[[k]]["g"], "adjacency" -> preparedH3[[k]]["A"], "edges" -> preparedH3[[k]]["edges"], "Aut" -> preparedH3[[k]]["aut"]|>, {k, Length[preparedH3]}],
    "readme" -> "F_{3,0} of K_P2 x C^2 (fibrewise torus) by the Givental-Teleman graph sum with Hodge vertices (GW theory of C^3, e_T weights {eps1+eps2, -eps1, -eps2} of the normal bundle of the fixed P^2) and R' = R_Z h^(-1). graphTerms: the 42 stable graphs of genus 3 in the order of StableGraph-3.bin (see 'graphs'), as functions of y, XX = y c1'/c1, eps1, eps2; orbitTerms: per decorated graph (keys as in orbits.m); polar / polarOrbits: coefficients of t_c^-4, ..., t_c^-1 with eps1 = x, eps2 = 1. The terms sum to F30 (ref/F20_F30_y.m)."|>;
  Put[results, "hodge_vertex_genus3_results.m"];
  (* summary *)
  emit[s_] := AppendTo[out, s];
  emit["F_{3,0} of K_P2 x C^2 (fibrewise torus): Hodge-vertex expansion (weights {eps1+eps2, -eps1, -eps2}, R' = R_Z h^(-1)); E1 = eps1, E2 = eps2, XX = y c1'/c1"];
  emit["total = F30: " <> ToString[Union[numchk] === {0}] <> " (symbolic: " <> ToString[symchk] <> "); total polar part " <> fmt[Factor /@ Total[polar]] <> "; gap: " <> ToString[gapok]];
  emit[""];
  Do[Module[{gr = preparedH3[[k]], tm = Together[graphTerms[[k]]], ids = Select[OTH3["ids"], #[[1]] == k &]},
    emit["=== stable graph " <> ToString[k] <> ": genera " <> ToString[gr["g"]] <> ", edges " <> ToString[gr["edges"]] <> ", |Aut| = " <> ToString[gr["aut"]] <> ", " <> ToString[Length[ids]] <> " decorated graphs"];
    emit["  polar {t_c^-4, t_c^-3, t_c^-2, t_c^-1} = " <> fmt[Factor /@ polar[[k]]]];
    Do[Module[{c = Factor[Coefficient[Numerator[tm], XX, j]/Denominator[tm]]}, If[c =!= 0, emit["    XX^" <> ToString[j] <> " : " <> fmt[c]]]], {j, 0, Exponent[Numerator[tm], XX]}];
    Do[Module[{d = OTH3["descriptions"][id], to = Together[terms[id]]},
      emit["  -- " <> ToString[id] <> " psi " <> ToString[d["psiOnHalfEdges"]] <> " dil " <> ToString[d["dilatonLeaves"]] <> "  polar " <> fmt[Factor /@ polarOrbits[id]]];
      Do[Module[{c = Factor[Coefficient[Numerator[to], XX, j]/Denominator[to]]}, If[c =!= 0, emit["       XX^" <> ToString[j] <> " : " <> fmt[c]]]], {j, 0, Exponent[Numerator[to], XX]}]], {id, ids}];
    emit[""]], {k, Length[preparedH3]}];
  Export["hodge_vertex_genus3_summary.txt", StringRiffle[out, "\n"], "Text"];
  log["written hodge_vertex_genus3_results.m and hodge_vertex_genus3_summary.txt; " <> ToString[Round[AbsoluteTime[] - tt]] <> " s"]],

 True, Print["unknown stage ", stage]];
