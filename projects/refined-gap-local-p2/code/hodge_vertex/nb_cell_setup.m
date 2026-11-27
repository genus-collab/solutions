(* data directory: this notebook's folder if it holds the sidecars (cPP_z6.m, StableGraph-2/3.bin, ref/F20_F30_y.m, hv_orbits_g3.mx,
   hv3_samples_B.mx, hodge_vertex_genus2_results.m), otherwise its subfolder hodge_vertex/ *)
datadir = If[FileExistsQ[FileNameJoin[{NotebookDirectory[], "cPP_z6.m"}]], NotebookDirectory[], FileNameJoin[{NotebookDirectory[], "hodge_vertex"}]];
SetDirectory[datadir]; Print["data directory: ", Directory[]];
e1 = \[Epsilon][1]; e2 = \[Epsilon][2];
wB = {e1 + e2, -e1, -e2};          (* e_T weights of the normal bundle O(-3) + O + O of the fixed P^2 *)
usubB = usubOf[wB];
{Fref20, Fref30} = Get["ref/F20_F30_y.m"];   (* F_{2,0}, F_{3,0} of localp2CY5_genus3.nb in y, XX, eps *)
(* reconstruction of every decorated graph from a list of sample records {{t0, x0}, vectors, diagnostics, ...} *)
reconstructAll[recs_, OTx_, ts_, xs_, gg_] := Module[{samp, sv},
   Do[samp[r[[1, 1]]][r[[1, 2]]] = r[[2]], {r, recs}];
   Association[Table[Module[{o = oo}, Clear[sv]; Do[sv[t0][x0] = samp[t0][x0][[o]], {t0, ts}, {x0, xs}];
      OTx["ids"][[o]] -> reconstructOrbitH[sv, ts, xs, gg]], {oo, Length[OTx["ids"]]}]]];
graphSums[terms_, OTx_, ngraphs_] := Table[Together[Total[Table[terms[id], {id, Select[OTx["ids"], #[[1]] == k &]}]]], {k, ngraphs}];
polarOf[term_, gg_] := polarOfTerm[term, gg] /. SeriesCoefficient[a_, _] /; FreeQ[a, tcv] :> 0;
sampleAll[pts_, gg_, OTx_, rows_] := (If[$KernelCount == 0, Quiet[LaunchKernels[]]];
   If[$KernelCount > 0, DistributeDefinitions["Global`"];
    ParallelTable[computeSampleH[p[[1]], p[[2]], gg, OTx, rows, usubB /. {e1 -> p[[2]], e2 -> 1}], {p, pts}],
    Table[computeSampleH[p[[1]], p[[2]], gg, OTx, rows, usubB /. {e1 -> p[[2]], e2 -> 1}], {p, pts}]]);
