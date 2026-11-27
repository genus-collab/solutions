(* data directory: this notebook's folder if it holds hvc_lib.m (StableGraph-2/3.bin, ref/resolvedconifold_results.m, the stored
   hvc_results_g*.m), otherwise its subfolder hodge_vertex_conifold/.  The code above is inlined from graphs.m, orbits.m, rc_lib.m
   and hvc_lib.m; their Get[...] lines were removed. *)
datadir = If[FileExistsQ[FileNameJoin[{NotebookDirectory[], "hvc_lib.m"}]], NotebookDirectory[], FileNameJoin[{NotebookDirectory[], "hodge_vertex_conifold"}]];
SetDirectory[datadir]; Print["data directory: ", Directory[]];
wB = modeWeights["B"]; wA = modeWeights["A"];   (* {eps1+eps2, -eps1, -eps2} and its negative; e1, e2 are the symbols Epsilon[1], Epsilon[2] *)
(* one complete run: R' = R_Z h_w^s with Q and x = eps1/eps2 symbolic (alpha0 = 0, alpha1 = eps1 = x, alpha2 = eps2 = 1), decorated
   graphs with Hodge (or point) vertices, exact term and polar part at Q = 1 of every decorated graph; returns an Association
   with the same keys as the stored hvc_results_g<g>_<tag>.m files *)
runMode[gg_, mode_, s_] := Module[{nfo, usub, entries, stable, total, Fref, pol, gap},
   buildFactors[3 gg - 3, mode, s, x]; nfo = loadOrbitsHVC[gg, mode]; usub = usubMode[mode, x];
   entries = Table[Module[{r = orbitTerm[OT["graphs"][id], usub, gg]}, {id, OT["orbits"][[id[[1]], id[[2]]]], Normal[OT["descriptions"][id]], r[[1]], r[[2]], r[[3]]}], {id, OT["ids"]}];
   stable = Table[Module[{es = Select[entries, #[[1, 1]] == k &]}, <|"graph" -> k, "genera" -> prepared[[k]]["g"], "edges" -> prepared[[k]]["edges"], "Aut" -> prepared[[k]]["aut"], "decoratedGraphs" -> Length[es],
        "singularDecoratedGraphs" -> Count[es, e_ /; Union[Together /@ e[[5]]] =!= {0}], "term" -> Together[Total[es[[All, 4]]]], "polar" -> (Together /@ Total[es[[All, 5]]])|>], {k, Length[prepared]}];
   total = Together[Total[stable[[All, "term"]]]]; Fref = FrefHVC[gg, x]; pol = Together /@ Total[stable[[All, "polar"]]]; gap = gapHVC[gg, x];
   <|"genus" -> gg, "mode" -> mode, "s" -> s, "counts" -> nfo, "stable" -> stable, "entries" -> entries, "total" -> total, "totalOK" -> (Together[total - Fref] === 0), "polar" -> pol, "gapOK" -> (Union[Together /@ (pol - gap)] === {0}), "symmetric" -> And @@ entries[[All, 6]]|>];
regQ[v_] := Union[Together /@ v] === {0};
tableOf[res_] := Grid[Prepend[Table[With[{st = res["stable"][[k]]}, {k, st["genera"], st["edges"], st["Aut"], st["decoratedGraphs"], st["singularDecoratedGraphs"], If[regQ[st["polar"]], "regular", Factor /@ st["polar"]]}], {k, Length[res["stable"]]}],
    {"graph", "genera", "edges", "|Aut|", "decorated graphs", "singular decorated graphs", "polar part {t^-(2g-2), ..., t^-1} at Q = 1 (eps1 = x, eps2 = 1)"}], Frame -> All];
(* the structural statements of a run: singular decorated graphs versus vertex genera, pole-order bound Sum_v max(2 g_v - 2, 0) *)
structure[res_] := Module[{st = res["stable"], en = res["entries"], gg = res["genus"], ge2, sg, ord},
   ge2 = Select[Range[Length[st]], Max[st[[#]]["genera"]] >= 2 &]; sg = Select[Range[Length[st]], st[[#]]["singularDecoratedGraphs"] > 0 &];
   ord = Table[{Module[{p = Position[Together /@ e[[5]], _?(# =!= 0 &), {1}, 1, Heads -> False]}, If[p === {}, 0, 2 gg - 1 - p[[1, 1]]]], Total[Max[2 # - 2, 0] & /@ st[[e[[1, 1]]]]["genera"]]}, {e, en}];
   <|"graphs with a vertex of genus >= 2" -> ge2, "graphs containing singular decorated graphs" -> sg, "equal" -> (ge2 === sg), "pole order <= bound for every decorated graph" -> AllTrue[ord, #[[1]] <= #[[2]] &], "tally of {pole order, bound}" -> Sort[Tally[ord]]|>];
