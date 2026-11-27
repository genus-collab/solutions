(* hv_analysis_g4.wls -- structural analysis of the genus-4 Hodge-vertex decomposition of F_{4,0}(K_P2 x C^2, T_F) from the
   cluster results hodge_stablegraphs_genus4.m (379 stable graphs) and hodge_graphterms_genus4.m (27172 decorated graphs).
   Usage: wolframscript -file hv_analysis_g4.wls      (run in hodge_vertex/, where lp2_lib.m, cPP_z6.m and the result files are)
   Output: hodge_vertex_genus4_summary.txt (human-readable), hodge_vertex_genus4_analysis.m (Association), hodge_vertex_genus4_analysis.log *)
SetDirectory[DirectoryName[$InputFileName]];
Get["lp2_lib.m"];
gg = 4; e1 = \[Epsilon][1]; e2 = \[Epsilon][2];
strm = OpenWrite["hodge_vertex_genus4_analysis.log"];
log[m_] := (Print[m]; WriteString[strm, DateString[{"Hour", ":", "Minute", ":", "Second"}], "  ", m, "\n"]);
tt0 = AbsoluteTime[];
stable = Get["hodge_stablegraphs_genus4.m"];
entries = Rest[Get["hodge_graphterms_genus4.m"]];
ngr = Length[stable];
log["loaded " <> ToString[ngr] <> " stable graphs and " <> ToString[Length[entries]] <> " decorated graphs (" <> ToString[Round[AbsoluteTime[] - tt0]] <> " s)"];
log["entry format: id " <> ToString[InputForm[entries[[1, 1]]]] <> ", key " <> ToString[InputForm[entries[[1, 2]]]] <> ", description " <> ToString[InputForm[entries[[1, 3]]]]];
tcSetup[2 gg + 4];
gap = gapCoefficient[gg] /. {e1 -> x, e2 -> 1};
res = <|"readme" -> "Genus-4 Hodge-vertex decomposition of F_{4,0}(K_P2 x C^2, T_F), mode B weights {eps1+eps2, -eps1, -eps2}; analysis of the cluster results of 2026-10-30."|>;

(* ---------- 1. polar parts of the stable graphs, canonicalised ---------- *)
pol = Table[Together /@ st["polar"], {st, stable}];
nonzero = Select[Range[ngr], Union[pol[[#]]] =!= {0} &];
g1ok = Together[pol[[1, 1]] - gap] === 0 && Union[Rest[pol[[1]]]] === {0};
log["stable graphs with a nonzero polar part (orders t_c^-6..t_c^-1, after Together): " <> ToString[nonzero] <> " of " <> ToString[ngr]];
log["graph 1 (genus-4 vertex) polar part = {gap coefficient, 0, 0, 0, 0, 0}: " <> ToString[g1ok] <> "; gap = " <> ToString[InputForm[Factor[gap]]]];
totpol = Together /@ Total[pol];
log["sum of all stable-graph polar parts = {gap, 0, 0, 0, 0, 0}: " <> ToString[Together[totpol[[1]] - gap] === 0 && Union[Rest[totpol]] === {0}]];
res["nonzeroPolarGraphs"] = nonzero; res["graph1PolarIsGap"] = g1ok; res["gap"] = Factor[gap];

(* ---------- 2. total = F40 ---------- *)
F40ref = Get["ref/genRmat_F40_cluster_result.m"][[1]] /. {q0 -> 0, q1 -> 0, q2 -> 0, q3 -> 0, q4 -> 0};
tt = AbsoluteTime[];
Ftot = Together[Total[stable[[All, "term"]]]];
totok = Together[Ftot - F40ref] === 0;
log["sum of the 379 stable-graph terms = F40 (q -> 0 part of the twisted cluster result), symbolically: " <> ToString[totok] <> " (" <> ToString[Round[AbsoluteTime[] - tt]] <> " s)"];
res["totalIsF40"] = totok;

(* ---------- 3. vertex genera versus singular decorated graphs ---------- *)
maxg = Table[Max[st["genera"]], {st, stable}];
sing = Table[st["singularDecoratedGraphs"], {st, stable}];
ndec = Table[st["decoratedGraphs"], {st, stable}];
setG = Select[Range[ngr], maxg[[#]] >= 2 &]; setS = Select[Range[ngr], sing[[#]] > 0 &];
log["stable graphs with a vertex of genus >= 2: " <> ToString[Length[setG]] <> "; with singular decorated graphs: " <> ToString[Length[setS]] <> "; the same set: " <> ToString[setG === setS]];
byg = Table[{h, Count[maxg, h], Total[Pick[ndec, maxg, h]], Total[Pick[sing, maxg, h]]}, {h, 0, 4}];
log["by maximal vertex genus {h, stable graphs, decorated graphs, singular decorated graphs}: " <> ToString[byg]];
res["byMaxVertexGenus"] = byg; res["singularSetEqualsGenusGe2Set"] = setG === setS;
(* regular decorated graphs inside graphs with a vertex of genus >= 2 *)
regIn = Select[entries, MemberQ[setG, #[[1, 1]]] && Union[#[[5]]] === {0} &];
log[ToString[Length[regIn]] <> " regular decorated graphs inside the " <> ToString[Length[setG]] <> " stable graphs with a vertex of genus >= 2; by graph: " <> ToString[Tally[regIn[[All, 1, 1]]]]];
log["their keys {graph, psi on half-edges, dilaton leaves} (first 40): " <> ToString[InputForm[Take[Table[{e[[1, 1]], e[[2]]}, {e, regIn}], UpTo[40]]]]];
res["regularDecorationsInGenusGe2Graphs"] = Table[{e[[1]], e[[2]]}, {e, regIn}];

(* ---------- 4. pole orders of the decorated graphs: leading order versus the maximal vertex genus ---------- *)
leadOrder[p_] := Module[{pos = Position[p, _?(# =!= 0 &), {1}, 1, Heads -> False]}, If[pos === {}, 0, 2 gg - 2 - (pos[[1, 1]] - 1)]];   (* 0 = regular; k = pole t_c^-k *)
lead = Table[leadOrder[e[[5]]], {e, entries}];
gOf = Table[maxg[[e[[1, 1]]]], {e, entries}];
prof = Table[{h, Tally[Pick[lead, gOf, h]] // Sort}, {h, 0, 4}];
log["leading pole order of the decorated graphs, by maximal vertex genus h: {h, {{order, count}, ...}}: " <> ToString[prof]];
(* the pole order of a decorated graph is bounded by the sum over its vertices of (2 g_v - 2)_+ : each Hodge vertex of genus
   h >= 2 contributes a conifold singularity of order 2h - 2 (the order of the genus-h gap), genus 0 and 1 vertices none *)
sumBound = Table[Total[Max[2 # - 2, 0] & /@ stable[[e[[1, 1]]]]["genera"]], {e, entries}];
bound = And @@ Table[lead[[i]] <= sumBound[[i]], {i, Length[entries]}];
log["every decorated graph has pole order <= Sum_v max(2 g_v - 2, 0) over its vertices: " <> ToString[bound]];
attained = Table[{k, Max[Pick[lead, entries[[All, 1, 1]], k]], Total[Max[2 # - 2, 0] & /@ stable[[k]]["genera"]]}, {k, setG}];
log["{graph, max pole order among its decorations, Sum_v (2 g_v - 2)_+} for the graphs with a vertex of genus >= 2; bound attained in every case: " <> ToString[And @@ Table[a[[2]] == a[[3]], {a, attained}]] <> "; " <> ToString[attained]];
oddzero = And @@ Table[Union[e[[5, {2, 4, 6}]]] === {0}, {e, entries}];
log["odd orders t_c^-5, t_c^-3, t_c^-1 vanish for every decorated graph: " <> ToString[oddzero]];
nz = Count[regIn, e_ /; Together[e[[4]]] === 0];
log["of the " <> ToString[Length[regIn]] <> " regular decorated graphs inside graphs with a vertex of genus >= 2, identically zero terms: " <> ToString[nz]];
res["poleOrderProfile"] = prof; res["poleOrderBound"] = bound; res["poleOrderAttained"] = attained; res["oddOrdersVanish"] = oddzero; res["regularDecorationsZero"] = nz;

(* ---------- 5. denominators: powers of (1 + 27 y), y, eps1 eps2, eps1 + eps2 ---------- *)
powOf[d_, f_, v_] := Module[{k = 0, q = d}, While[Together[q /. First[Solve[f == 0, v]]] === 0 && k < 40, q = PolynomialQuotient[q, f, v]; k++]; k];
denData[t_] := Module[{d = Denominator[Together[t]]},
   {powOf[d, 1 + 27 y, y], Exponent[d, y, Min] /. -Infinity -> 0, Exponent[d, e1, Min], Exponent[d, e2, Min], powOf[d, e1 + e2, e1]}];
tt = AbsoluteTime[];
dst = Table[denData[st["term"]], {st, stable}];
log["stable-graph denominators {power of (1+27y), min power of y, of eps1, of eps2, of (eps1+eps2)}: max over graphs " <> ToString[Max /@ Transpose[dst]] <> "; graphs with eps-denominators: " <> ToString[Count[dst, {_, _, a_, b_, c_} /; a + b + c > 0]] <> " (" <> ToString[Round[AbsoluteTime[] - tt]] <> " s)"];
log["pole order at the conifold is bounded by the power of (1+27y) in the denominator; max over stable graphs: " <> ToString[Max[dst[[All, 1]]]] <> " (<= 6 means the polar lists t_c^-6..t_c^-1 are complete)"];
res["stableDenominators"] = dst;
regTally = Sort[Tally[dst[[Complement[Range[ngr], nonzero], 1]]]];
log["regular stable graphs: tally {power of (1+27y) in the denominator of the term, graphs}: " <> ToString[regTally] <> " -- off-shell poles of the XX-polynomial coefficients that cancel on the curve XX = XX(y)"];
res["regularGraphsConifoldPowerTally"] = regTally;
tt = AbsoluteTime[];
dpow = Table[powOf[Denominator[Together[e[[4]]]], 1 + 27 y, y], {e, entries}];
log["decorated graphs: max power of (1+27y) in the denominators: " <> ToString[Max[dpow]] <> "; tally: " <> ToString[Sort[Tally[dpow]]] <> " (" <> ToString[Round[AbsoluteTime[] - tt]] <> " s)"];
res["decoratedConifoldPowerTally"] = Sort[Tally[dpow]];
polyXX = Table[PolynomialQ[st["term"], XX], {st, stable}];
degXX = Table[Exponent[Together[st["term"]], XX], {st, stable}];
log["stable-graph terms polynomial in XX: " <> ToString[Count[polyXX, True]] <> " of " <> ToString[ngr] <> "; XX-degree tally: " <> ToString[Sort[Tally[degXX]]]];
yfree = Select[Range[ngr], FreeQ[stable[[#]]["term"], y] &];
log["stable-graph terms independent of y: " <> ToString[Length[yfree]] <> " graphs " <> ToString[yfree]];
res["yIndependentGraphs"] = yfree;

(* ---------- 6. the genus-4 vertex graph: its 97 decorated graphs ---------- *)
g1 = Select[entries, #[[1, 1]] == 1 &];
g1tab = Table[{e[[2]], Factor[e[[5, 1]]], Factor[e[[5, 3]]], Factor[e[[5, 5]]]}, {e, g1}];
log["graph 1: " <> ToString[Length[g1]] <> " decorated graphs; regular ones: " <> ToString[InputForm[Select[g1, Union[#[[5]]] === {0} &][[All, 2]]]]];
ratios = DeleteDuplicates[Table[Together[e[[5, 1]]/gap], {e, g1}]];
log["graph 1: t_c^-6 coefficients proportional to the gap polynomial? distinct ratios: " <> ToString[Length[ratios]] <> If[Length[ratios] <= 6, " " <> ToString[InputForm[ratios]], ""]];
log["graph 1: sum of t_c^-6 coefficients = gap: " <> ToString[Together[Total[g1[[All, 5, 1]]] - gap] === 0] <> "; t_c^-4 and t_c^-2 sums vanish: " <> ToString[Together[Total[g1[[All, 5, 3]]]] === 0 && Together[Total[g1[[All, 5, 5]]]] === 0]];
res["graph1Decorations"] = g1tab;

(* ---------- 7. the cancelling graphs: pole orders present among their decorated graphs ---------- *)
canc = Table[{k, stable[[k]]["genera"], stable[[k]]["edges"], ndec[[k]], sing[[k]], Max[Pick[lead, entries[[All, 1, 1]], k]]}, {k, Complement[setS, {1}]}];
log["stable graphs with singular decorated graphs but regular sum {graph, genera, edges, decorated, singular, max pole order}: " <> ToString[canc]];
res["cancellingGraphs"] = canc;

(* ---------- 8. simplest terms ---------- *)
lc = Table[LeafCount[st["term"]], {st, stable}];
small = Ordering[lc, 6];
log["six simplest stable-graph terms (graph, genera, edges, term): " <> ToString[InputForm[Table[{k, stable[[k]]["genera"], stable[[k]]["edges"], Factor[stable[[k]]["term"]]}, {k, small}]]]];
res["simplestTerms"] = Table[{k, stable[[k]]["genera"], stable[[k]]["edges"], Factor[stable[[k]]["term"]]}, {k, small}];

Put[res, "hodge_vertex_genus4_analysis.m"];
summary = StringRiffle[{
  "Genus-4 Hodge-vertex decomposition of F_{4,0}(K_P2 x C^2, T_F) -- mode B weights {eps1+eps2, -eps1, -eps2}, R' = R_Z h^(-1)",
  "Cluster run Stanage 2026-10-29/30: 680 sample points, 379 stable graphs, 27172 decorated graphs, 0 failed reconstructions.",
  "",
  "total of the 379 stable-graph terms = F40 (symbolic): " <> ToString[totok],
  "total polar part = {gap, 0, 0, 0, 0, 0}, gap = " <> ToString[InputForm[Factor[gap]]] <> ": " <> ToString[Together[totpol[[1]] - gap] === 0 && Union[Rest[totpol]] === {0}],
  "stable graphs with a nonzero polar part: " <> ToString[nonzero] <> " (graph 1 = the single genus-4 vertex); its polar part is exactly the gap: " <> ToString[g1ok],
  "stable graphs with a vertex of genus >= 2: " <> ToString[Length[setG]] <> " = stable graphs containing singular decorated graphs: " <> ToString[setG === setS],
  "by maximal vertex genus {h, stable graphs, decorated graphs, singular decorated graphs}: " <> ToString[byg],
  "leading pole order of decorated graphs by maximal vertex genus: " <> ToString[prof],
  "pole order <= Sum_v (2 g_v - 2)_+ for every decorated graph: " <> ToString[bound] <> "; bound attained in every graph with a vertex of genus >= 2: " <> ToString[And @@ Table[a[[2]] == a[[3]], {a, attained}]] <> "; odd orders vanish: " <> ToString[oddzero],
  "regular decorated graphs inside graphs with a vertex of genus >= 2: " <> ToString[Length[regIn]] <> ", of which identically zero: " <> ToString[nz],
  "max power of (1+27y) in the denominators: stable graphs " <> ToString[Max[dst[[All, 1]]]] <> ", decorated graphs " <> ToString[Max[dpow]] <> "; regular stable graphs by that power: " <> ToString[regTally],
  "stable-graph terms with eps-denominators: " <> ToString[Count[dst, {_, _, a_, b_, c_} /; a + b + c > 0]] <> "; max powers {eps1, eps2, eps1+eps2}: " <> ToString[Max /@ Transpose[dst[[All, 3 ;; 5]]]],
  "y-independent stable-graph terms: " <> ToString[yfree],
  "graph 1: regular decorated graphs " <> ToString[InputForm[Select[g1, Union[#[[5]]] === {0} &][[All, 2]]]] <> "; distinct ratios of t_c^-6 coefficients to the gap: " <> ToString[Length[ratios]],
  "cancelling graphs {graph, genera, edges, decorated, singular, max pole order}: " <> ToString[canc]}, "\n"];
Export["hodge_vertex_genus4_summary.txt", summary, "Text"];
log["written hodge_vertex_genus4_analysis.m and hodge_vertex_genus4_summary.txt; total " <> ToString[Round[AbsoluteTime[] - tt0]] <> " s"];
Close[strm];
