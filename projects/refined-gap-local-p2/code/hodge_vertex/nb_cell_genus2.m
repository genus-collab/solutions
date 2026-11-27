(* genus 2: 7 stable graphs, 27 decorated graphs with Hodge vertices; 26 t x 34 eps-ratio sample points (a few minutes) *)
gg = 2; KMAX = kmaxFor[gg];
preparedH2 = prepareGraphH /@ graphsWithLegs[2, 0]; OTH2 = orbitTable[preparedH2];
{Length[preparedH2], Total[Length[#["decs"]] & /@ preparedH2], Length[OTH2["ids"]]}
loadRmatrix[KMAX]; rows2 = rowsTimes[hCoeffs[wB, -1, KMAX], KMAX];   (* rows of R' = R_Z h^(-1) *)
ts2 = tsamplesFor[2]; xs2 = xsamplesH[34]; pts2 = Flatten[Table[{t0, x0}, {t0, ts2}, {x0, xs2}], 1];
AbsoluteTiming[recs2 = sampleAll[pts2, 2, OTH2, rows2];]
(* the sum over the decorated graphs equals F_{2,0} exactly at every sample point *)
Union[Flatten[Table[Total[r[[2]]] - refVector[Fref20, r[[1, 1]], r[[1, 2]], 2], {r, recs2}]]]
terms2 = reconstructAll[recs2, OTH2, ts2, xs2, 2]; Select[Keys[terms2], ListQ[terms2[#]] &]   (* failed reconstructions: none expected *)
graphTerms2 = graphSums[terms2, OTH2, Length[preparedH2]]; Together[Total[graphTerms2] - Fref20] === 0
(* agreement with the stored run (hodge_vertex_genus2_results.m, mode B) *)
stored2 = Get["hodge_vertex_genus2_results.m"]["B"]["graphTerms"]; Table[Together[graphTerms2[[k]] - stored2[[k]]] === 0, {k, Length[preparedH2]}]
(* conifold polar parts {t_c^-2, t_c^-1} (eps1 = x, eps2 = 1): everything sits in the genus-2 vertex graph *)
tcSetup[8]; Grid[Prepend[Table[{k, preparedH2[[k]]["g"], preparedH2[[k]]["edges"], Factor /@ polarOf[graphTerms2[[k]], 2]}, {k, Length[preparedH2]}], {"graph", "genera", "edges", "polar part"}], Frame -> All]
{Factor /@ Total[Table[polarOf[graphTerms2[[k]], 2], {k, Length[preparedH2]}]], Factor[gapCoefficient[2] /. {e1 -> x, e2 -> 1}]}
