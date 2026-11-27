(* orbits.m -- decorated stable graphs = orbits of labelled decorations under Aut(Gamma). Shared by both pipelines
   (identical prepareGraph structure on both sides, so the keys match across theories).
   A decorated graph is the orbit of {psi powers on the two half-edges of each edge, sorted dilaton lists per vertex}
   under vertex permutations preserving genera/legs/adjacency, permutations of parallel edges and flips of loops.
   Canonical representative = lexicographically smallest image. *)

autosOf[gr_] := Module[{gv = gr["g"], mv = gr["m"], A = gr["A"], ed = gr["edges"], V = gr["V"], perms, out = {}},
   perms = Select[Permutations[Range[V]], gv[[#]] == gv && mv[[#]] == mv && A[[#, #]] == A &];
   Do[Module[{img = Map[Sort[sigma[[#]]] &, ed], classes, choices},
     classes = GatherBy[Range[Length[ed]], img[[#]] &];
     choices = Tuples[Table[Module[{cl = c, targets = Flatten[Position[ed, img[[c[[1]]]], {1}]]}, Map[Thread[cl -> #] &, Permutations[targets]]], {c, classes}]];
     Do[Do[AppendTo[out, {sigma, Flatten[ch], flips}], {flips, Tuples[Table[If[ed[[e, 1]] == ed[[e, 2]], {0, 1}, {0}], {e, Length[ed]}]]}], {ch, choices}]], {sigma, perms}]; out];
decData[gr_, dec_] := Module[{ed = gr["edges"], hsAt = gr["hsAt"], psi, dils},
   psi = Table[{0, 0}, Length[ed]];
   Do[Do[psi[[gr["H"][[hsAt[[v, r]], 1]], gr["H"][[hsAt[[v, r]], 2]]]] = dec[[v, 1, r]], {r, Length[hsAt[[v]]]}], {v, gr["V"]}];
   dils = Table[Sort[dec[[v, 3]]], {v, gr["V"]}]; {psi, dils}];
applyAut[{psi_, dils_}, {sigma_, emap_, flips_}, ed_] := Module[{psi2 = Table[{0, 0}, Length[ed]], dils2 = dils},
   Do[Module[{e2 = e /. emap, a = psi[[e]], tgt, ord}, tgt = {sigma[[ed[[e, 1]]]], sigma[[ed[[e, 2]]]]};
     ord = If[tgt == ed[[e2]], {1, 2}, {2, 1}]; If[ed[[e2, 1]] == ed[[e2, 2]] && flips[[e]] == 1, ord = Reverse[ord]];
     psi2[[e2, ord[[1]]]] = a[[1]]; psi2[[e2, ord[[2]]]] = a[[2]]], {e, Length[ed]}];
   Do[dils2[[sigma[[v]]]] = dils[[v]], {v, Length[dils]}]; {psi2, dils2}];
canonOf[gr_, dec_, auts_] := First[Sort[Table[applyAut[decData[gr, dec], a, gr["edges"]], {a, auts}]]];
(* {orbitKeys, orbits}: orbitKeys[[k]] = canonical key of every labelled decoration of graph k; orbits[[k]] = sorted distinct keys *)
computeOrbits[prep_] := Module[{orbitKeys, orbits},
   orbitKeys = Table[Module[{auts = autosOf[gr]}, Table[canonOf[gr, dec, auts], {dec, gr["decs"]}]], {gr, prep}];
   orbits = Table[Union[orbitKeys[[k]]], {k, Length[prep]}];
   {orbitKeys, orbits}];
orbitDecs[prep_, orbitKeys_, k_, okey_] := prep[[k]]["decs"][[Flatten[Position[orbitKeys[[k]], okey, {1}]]]];
orbitDescription[gr_, okey_] := <|"genera" -> gr["g"], "legs" -> gr["m"], "adjacency" -> gr["A"], "edges" -> gr["edges"], "psiOnHalfEdges" -> okey[[1]], "dilatonLeaves" -> okey[[2]], "Aut" -> gr["aut"], "dilatonFree" -> (Union[Flatten[okey[[2]]]] === {})|>;
(* full orbit table for a list of prepared graphs: orbitIds = {{k, j}, ...}, orbitGraphs[{k, j}] = graph restricted to the orbit *)
orbitTable[prep_] := Module[{orbitKeys, orbits, ids, graphs},
   {orbitKeys, orbits} = computeOrbits[prep];
   ids = Flatten[Table[{k, j}, {k, Length[prep]}, {j, Length[orbits[[k]]]}], 1];
   graphs = Association[Table[id -> <|prep[[id[[1]]]], "decs" -> orbitDecs[prep, orbitKeys, id[[1]], orbits[[id[[1]], id[[2]]]]]|>, {id, ids}]];
   <|"orbitKeys" -> orbitKeys, "orbits" -> orbits, "ids" -> ids, "graphs" -> graphs,
     "descriptions" -> Association[Table[id -> orbitDescription[prep[[id[[1]]]], orbits[[id[[1]], id[[2]]]]], {id, ids}]]|>];
