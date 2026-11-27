(* graphs.m -- stable graphs, decorations and Witten-Kontsevich numbers (shared by lp2_lib.m and rc_lib.m).
   Verbatim from genRmat_g4_lib.m / localp2CY5_genus3.nb / resolvedconifold_genus3.nb.
   graph = {genera, legs, adjacency matrix (diagonal = loops)}; StableGraph-<g>.bin must be in the current directory. *)

(* ---------- Witten-Kontsevich numbers ---------- *)
df[n_] := If[n <= 0, 1, n df[n - 2]];
splitSum[g_, k_, rest_] := Module[{n = Length[rest], idx, tot = 0},
   Do[idx = Subsets[Range[n]][[m]];
    tot += Sum[df[2 r + 1] df[2 (k - 1 - r) + 1] Sum[wk[g1, Join[{r}, rest[[idx]]]] wk[g - g1, Join[{k - 1 - r}, rest[[Complement[Range[n], idx]]]]], {g1, 0, g}], {r, 0, k - 1}], {m, 2^n}]; tot];
wk[g_, ds_List] := wk[g, ds] = Module[{d = Sort[ds], n, k, rest, tot}, n = Length[d];
   If[g < 0 || n == 0 || Min[d] < 0 || Total[d] != 3 g - 3 + n || 2 g - 2 + n <= 0, Return[0]];
   If[g == 0 && d == {0, 0, 0}, Return[1]]; If[g == 1 && d == {1}, Return[1/24]];
   If[MemberQ[d, 0], rest = Delete[d, First[Position[d, 0]]]; Return[Sum[wk[g, ReplacePart[rest, j -> rest[[j]] - 1]], {j, Length[rest]}]]];
   k = d[[1]] - 1; rest = Rest[d];
   tot = Sum[df[2 k + 2 rest[[j]] + 1]/df[2 rest[[j]] - 1] wk[g, ReplacePart[rest, j -> rest[[j]] + k]], {j, Length[rest]}];
   tot += 1/2 Sum[df[2 r + 1] df[2 (k - 1 - r) + 1] wk[g - 1, Join[{r, k - 1 - r}, rest]], {r, 0, k - 1}];
   tot += 1/2 splitSum[g, k, rest]; tot/df[2 k + 3]];

(* ---------- stable graphs ---------- *)
tabtoadj[diag_, up_] := Module[{K = Length[diag], A}, A = ConstantArray[0, {K, K}];
   Do[A[[j, j]] = diag[[j]]; Do[A[[i, j]] = A[[j, i]] = up[[j, i - j]], {i, j + 1, K}], {j, 1, K}]; A];
graphlist[raw_] := Module[{count = 1, K, g, m, l, atab, a, out = {}, adj},
   While[count <= Length[raw],
    K = raw[[count]]; g = Take[raw, {count + 1, count + K}]; m = Take[raw, {count + K + 1, count + 2 K}];
    l = Take[raw, {count + 2 K + 1, count + 3 K}]; a = {};
    If[K > 1, atab = Take[raw, {count + 3 K + 1, count + 3 K + K (K - 1)/2}];
     Do[AppendTo[a, Take[atab, K - j]]; atab = Drop[atab, K - j], {j, 1, K - 1}]];
    adj = tabtoadj[l, a]; AppendTo[out, {g, m, adj}];
    count += 3 K + Length[Flatten[a]] + 1]; out];
loadGraphs[gg_] := graphlist[DeleteCases[Delete[Delete[BinaryReadList["StableGraph-" <> ToString[gg] <> ".bin"], 1], -1], 13]];
(* hand-coded graphs with legs for the genus <= 1 cases *)
graphsWithLegs[0, 3] = {{{0}, {3}, {{0}}}};
graphsWithLegs[0, 4] = {{{0}, {4}, {{0}}}, {{0, 0}, {2, 2}, {{0, 1}, {1, 0}}}};
graphsWithLegs[0, 5] = {{{0}, {5}, {{0}}}, {{0, 0}, {2, 3}, {{0, 1}, {1, 0}}}, {{0, 0, 0}, {2, 1, 2}, {{0, 1, 0}, {1, 0, 1}, {0, 1, 0}}}};
graphsWithLegs[1, 1] = {{{1}, {1}, {{0}}}, {{0}, {1}, {{1}}}};
graphsWithLegs[g_, 0] := loadGraphs[g];

edgesOf[A_] := Join[Flatten[Table[ConstantArray[{i, i}, A[[i, i]]], {i, Length[A]}], 1],
   Flatten[Table[ConstantArray[{i, j}, A[[i, j]]], {i, Length[A]}, {j, i + 1, Length[A]}], 2]];
autCount[gv_, mv_, A_] := Module[{V = Length[gv], perms},
   perms = Select[Permutations[Range[V]], gv[[#]] == gv && mv[[#]] == mv && A[[#, #]] == A &];
   Length[perms] Product[A[[i, j]]!, {i, 1, V}, {j, i + 1, V}] Product[2^A[[i, i]] A[[i, i]]!, {i, 1, V}]];

(* vertex decoration options: for a vertex with genus g, hs half-edges, m legs:
   list of {hpowers, legpowers, dilpowers, weight}, weight = wk * 1/(m! nd!) *)
vertexOptions[g_, hs_, m_] := Module[{res = {}, dim, nd, dils, hp, lp, rest, budget},
   Do[dim = 3 g - 3 + hs + m + nd;
    If[dim < 2 nd, Continue[]];
    Do[budget = dim - Total[dils];
     Do[Do[
       AppendTo[res, {hp, lp, dils, wk[g, Join[hp, lp, dils]]/(m! nd!)}],
       {lp, Select[Tuples[Range[0, budget - Total[hp]], m], Total[#] == budget - Total[hp] &]}],
      {hp, Select[Tuples[Range[0, budget], hs], Total[#] <= budget &]}],
     {dils, If[nd == 0, {{}}, Flatten[Table[Permutations /@ Select[IntegerPartitions[s, {nd}], Min[#] >= 2 &], {s, 2 nd, dim}], 2]]}],
    {nd, 0, Max[0, 3 g - 3 + hs + m]}];
   Select[res, #[[4]] =!= 0 &]];

(* prepare a graph: half-edge bookkeeping and decoration list *)
prepareGraph[{gv_, mv_, A_}] := Module[{V = Length[gv], ed, H, hv, hsAt, opts, decs, aut},
   ed = edgesOf[A]; H = Flatten[Table[{{e, 1}, {e, 2}}, {e, Length[ed]}], 1];
   hv = Table[ed[[h[[1]], h[[2]]]], {h, H}];
   hsAt = Table[Flatten[Position[hv, v]], {v, V}];   (* indices into H of half-edges at v *)
   opts = Table[vertexOptions[gv[[v]], Length[hsAt[[v]]], mv[[v]]], {v, V}];
   decs = Tuples[opts];
   aut = autCount[gv, mv, A];
   <|"g" -> gv, "m" -> mv, "A" -> A, "edges" -> ed, "H" -> H, "hsAt" -> hsAt, "decs" -> decs, "aut" -> aut, "V" -> V|>];
