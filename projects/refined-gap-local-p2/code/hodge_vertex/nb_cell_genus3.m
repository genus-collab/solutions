(* genus 3: 42 stable graphs; decorated graphs with Hodge vertices from the sidecar hv_orbits_g3.mx (or recomputed) *)
gg = 3; KMAX = kmaxFor[gg];
If[FileExistsQ["hv_orbits_g3.mx"], Get["hv_orbits_g3.mx"], preparedH3 = prepareGraphH /@ graphsWithLegs[3, 0]; OTH3 = orbitTable[preparedH3]];
{Length[preparedH3], Total[Length[#["decs"]] & /@ preparedH3], Length[OTH3["ids"]]}
loadRmatrix[KMAX]; rows3 = rowsTimes[hCoeffs[wB, -1, KMAX], KMAX];
ts3 = tsamplesFor[3]; xs3 = xsamplesH[30]; pts3 = Flatten[Table[{t0, x0}, {t0, ts3}, {x0, xs3}], 1];
(* sampling: 780 points; the sidecar hv3_samples_B.mx holds the records of the stored run (delete it to recompute; see the
   README for the timing per point) *)
If[FileExistsQ["hv3_samples_B.mx"], Get["hv3_samples_B.mx"], AbsoluteTiming[recs3 = sampleAll[pts3, 3, OTH3, rows3];]; DumpSave["hv3_samples_B.mx", recs3]];
Length[recs3]
(* the sum over the decorated graphs equals F_{3,0} exactly at every sample point *)
Union[Flatten[Table[Total[r[[2]]] - refVector[Fref30, r[[1, 1]], r[[1, 2]], 3], {r, recs3}]]]
AbsoluteTiming[terms3 = reconstructAll[recs3, OTH3, ts3, xs3, 3];]
Select[Keys[terms3], ListQ[terms3[#]] &]   (* failed reconstructions: none expected *)
graphTerms3 = graphSums[terms3, OTH3, Length[preparedH3]]; F30hodge = Together[Total[graphTerms3]];
SeedRandom[4]; Table[Together[(F30hodge - Fref30) /. {y -> RandomInteger[{1, 50}]/RandomInteger[{51, 97}], XX -> RandomInteger[{-30, 30}]/RandomInteger[{1, 17}], e1 -> RandomInteger[{1, 40}], e2 -> RandomInteger[{1, 40}]}], {5}]
TimeConstrained[Together[F30hodge - Fref30] === 0, 1800, "timeout"]
(* conifold frame: polar parts {t_c^-4, t_c^-3, t_c^-2, t_c^-1} (eps1 = x, eps2 = 1) *)
tcSetup[10]; polar3 = Table[polarOf[graphTerms3[[k]], 3], {k, Length[preparedH3]}];
{Factor /@ Total[polar3], Factor[gapCoefficient[3] /. {e1 -> x, e2 -> 1}]}   (* gap: leading = (-3)^2 3! [t^6] t^2/(4 sinh sinh), the rest 0 *)
Grid[Prepend[Table[{k, preparedH3[[k]]["g"], preparedH3[[k]]["edges"], preparedH3[[k]]["aut"], Length[Select[OTH3["ids"], #[[1]] == k &]], If[Union[polar3[[k]]] === {0}, "regular", Factor /@ polar3[[k]]]}, {k, Length[preparedH3]}], {"graph", "genera", "edges", "|Aut|", "decorated graphs", "polar part"}], Frame -> All]
polarOrbits3 = Association[Table[id -> polarOf[terms3[id], 3], {id, Keys[terms3]}]];
singular3 = Select[Keys[polarOrbits3], Union[polarOrbits3[#]] =!= {0} &]; {Length[singular3], singular3}
Put[<|"graphTerms" -> graphTerms3, "orbitTerms" -> terms3, "polar" -> polar3, "polarOrbits" -> polarOrbits3, "total" -> F30hodge, "descriptions" -> OTH3["descriptions"], "weights" -> wB|>, "hodge_vertex_genus3_results_nb.m"]
