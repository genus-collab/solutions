(* hv_collect.m -- collector for the Hodge-vertex per-decorated-graph pipeline.
   Usage:  wolframscript -file hv_collect.m genus
   Input:  hv_partial_g<genus>_c*.m and hv_terms_g<genus>_c*.m (all chunks), hv_orbits_g<genus>.mx,
           reference totals ref/F20_F30_y.m (genus 2, 3) or ref/genRmat_F40_cluster_result.m (genus 4, q -> 0),
           stored results ref/hodge_vertex_genus<2|3>_results.m (validation of the per-graph terms at genus 2 and 3).
   Output: hodge_F<genus>0.m (the total), hodge_graphterms_genus<genus>.m (all decorated graphs with polar parts),
           hodge_stablegraphs_genus<genus>.m (per stable graph: term, polar part), hodge_collect_genus<genus>.txt (table),
           hv_collect.log with VALIDATION lines. *)
SetDirectory[DirectoryName[$InputFileName]];
gg = ToExpression[$ScriptCommandLine[[2]]];
Get["hv_lib4.m"];
log[m_] := Module[{s = OpenAppend["hv_collect.log"]}, WriteString[s, DateString[{"Year", "-", "Month", "-", "Day", " ", "Hour", ":", "Minute", ":", "Second"}], "  [g=", gg, "] ", m, "\n"]; Close[s]; Print[m]];
e1 = \[Epsilon][1]; e2 = \[Epsilon][2];
Get["hv_orbits_g" <> ToString[gg] <> ".mx"];
pfiles = FileNames["hv_partial_g" <> ToString[gg] <> "_c*.m"];
If[pfiles === {}, log["no partial files found"]; Exit[1]];
parts = Get /@ pfiles;
F = Together[Total[parts[[All, 1]]]]; polar = Total[parts[[All, 2]]]; failed = Join @@ parts[[All, 3]];
log[ToString[Length[pfiles]] <> " chunk files; failed reconstructions: " <> ToString[Length[failed]] <> If[failed =!= {}, " " <> ToString[failed], ""]];
Put[F, "hodge_F" <> ToString[gg] <> "0.m"];
(* ---- total against the reference ---- *)
ref = If[gg <= 3, Get["ref/F20_F30_y.m"][[gg - 1]], Get["ref/genRmat_F40_cluster_result.m"][[1]] /. {q0 -> 0, q1 -> 0, q2 -> 0, q3 -> 0, q4 -> 0}];
SeedRandom[4];
pts = Table[{y -> RandomInteger[{1, 50}]/RandomInteger[{51, 97}], XX -> RandomInteger[{-30, 30}]/RandomInteger[{1, 17}], e1 -> RandomInteger[{1, 40}], e2 -> RandomInteger[{1, 40}]}, {5}];
numchk = Table[Together[(F - ref) /. p], {p, pts}];
symchk = TimeConstrained[Together[F - ref] === 0, 1200, "timeout"];
log["total minus reference at 5 random rational points: " <> ToString[numchk] <> "; symbolic: " <> ToString[symchk]];
log[If[Union[numchk] === {0}, "VALIDATION OK: sum of decorated graphs = reference F" <> ToString[gg] <> "0", "VALIDATION FAIL: total differs from the reference"]];
(* ---- gap ---- *)
expected = gapCoefficient[gg] /. {e1 -> x, e2 -> 1};
gapok = Together[polar[[1]] - expected] === 0 && Union[Together /@ Rest[polar]] === {0};
log["total polar part {t_c^-" <> ToString[2 gg - 2] <> ", ..., t_c^-1} (eps1 = x, eps2 = 1): " <> ToString[InputForm[Factor /@ polar]] <> "; expected leading " <> ToString[InputForm[Factor[expected]]]];
log[If[gapok, "VALIDATION OK: conifold gap", "VALIDATION FAIL: conifold gap"]];
(* ---- per-graph files ---- *)
tfiles = FileNames["hv_terms_g" <> ToString[gg] <> "_c*.m"];
entries = SortBy[Join @@ (Get /@ tfiles), First];
readme = "Local P2 x C2, fibrewise torus, genus " <> ToString[gg] <> ": the " <> ToString[Length[entries]] <> " decorated stable graphs of the Givental-Teleman sum for F_{" <> ToString[gg] <> ",0} with Hodge vertices (GW theory of C^3 with the e_T weights {eps1+eps2, -eps1, -eps2} of the normal bundle of the fixed P^2) and R' = R_Z h^(-1). " <>
  "Entry: {id {graph, orbit}, key {psi powers on the two half-edges of each edge, dilaton-leaf psi powers per vertex}, description (genera, legs, adjacency with loops on the diagonal, edges, Aut, dilatonFree), " <>
  "contribution as a function of y, XX = y c1'(y)/c1(y), eps1, eps2 (weight 1/|Aut| times the sum over labelled decorations and colourings), " <>
  "polar part in the conifold coordinate t_c = Delta + 11/18 Delta^2 + ... (Delta = 1 + 27 y) with eps1 = x, eps2 = 1: coefficients of t_c^-(2g-2), ..., t_c^-1. The contributions sum to F_{g,0}.";
Put[Prepend[entries, readme], "hodge_graphterms_genus" <> ToString[gg] <> ".m"];
regQ[v_] := Union[v] === {0};
ngr = Length[preparedH];
stable = Table[Module[{es = Select[entries, #[[1, 1]] == k &]}, <|"graph" -> k, "genera" -> preparedH[[k]]["g"], "edges" -> preparedH[[k]]["edges"], "Aut" -> preparedH[[k]]["aut"], "decoratedGraphs" -> Length[es],
     "singularDecoratedGraphs" -> Count[es, e_ /; ! regQ[e[[5]]]], "term" -> Together[Total[es[[All, 4]]]], "polar" -> Total[es[[All, 5]]]|>], {k, ngr}];
Put[stable, "hodge_stablegraphs_genus" <> ToString[gg] <> ".m"];
regular = Select[Range[ngr], regQ[stable[[#]]["polar"]] &];
log[ToString[ngr] <> " stable graphs; regular at the conifold: " <> ToString[Length[regular]] <> "; singular: " <> ToString[Table[{k, stable[[k]]["genera"], stable[[k]]["edges"], InputForm[Factor /@ stable[[k]]["polar"]]}, {k, Complement[Range[ngr], regular]}]]];
log["stable graphs containing singular decorated graphs: " <> ToString[Table[{k, stable[[k]]["genera"], stable[[k]]["singularDecoratedGraphs"], stable[[k]]["decoratedGraphs"]}, {k, Select[Range[ngr], stable[[#]]["singularDecoratedGraphs"] > 0 &]}]] <> "; graphs with a vertex of genus >= 2: " <> ToString[Select[Range[ngr], Max[preparedH[[#]]["g"]] >= 2 &]]];
tab = StringRiffle[Prepend[Table[StringRiffle[{k, stable[[k]]["genera"], stable[[k]]["edges"], stable[[k]]["Aut"], stable[[k]]["decoratedGraphs"], stable[[k]]["singularDecoratedGraphs"], If[regQ[stable[[k]]["polar"]], "regular", ToString[InputForm[Factor /@ stable[[k]]["polar"]]]]}, "  "], {k, ngr}],
   "graph  genera  edges  |Aut|  decorated graphs  singular decorated graphs  polar part {t_c^-(2g-2), ..., t_c^-1} (eps1 = x, eps2 = 1)"], "\n"];
Export["hodge_collect_genus" <> ToString[gg] <> ".txt", tab, "Text"];
log["wrote hodge_graphterms_genus" <> ToString[gg] <> ".m (" <> ToString[Length[entries]] <> " entries), hodge_stablegraphs_genus" <> ToString[gg] <> ".m, hodge_collect_genus" <> ToString[gg] <> ".txt"];
(* ---- genus-2/3 validation against the stored local runs (per decorated graph, matched by graph index + key) ---- *)
tcSetup[2 gg + 4];
If[gg <= 3 && FileExistsQ["ref/hodge_vertex_genus" <> ToString[gg] <> "_results.m"],
  Module[{res = Get["ref/hodge_vertex_genus" <> ToString[gg] <> "_results.m"], refTerms, refDesc, refKey, nt = 0, np = 0, nmiss = 0, refPolar},
   If[gg == 2, res = res["B"]];
   refTerms = res["orbitTerms"]; refDesc = res["descriptions"];
   refKey = Association[Table[{id[[1]], refDesc[id]["psiOnHalfEdges"], refDesc[id]["dilatonLeaves"]} -> id, {id, Keys[refTerms]}]];
   refPolar = If[KeyExistsQ[res, "polarOrbits"], res["polarOrbits"], Association[Table[id -> polarOfTerm[refTerms[id], gg] /. SeriesCoefficient[a_, _] /; FreeQ[a, tcv] :> 0, {id, Keys[refTerms]}]]];
   Do[Module[{key = {e[[1, 1]], e[[2, 1]], e[[2, 2]]}},
     If[! KeyExistsQ[refKey, key], nmiss++; Continue[]];
     If[Together[e[[4]] - refTerms[refKey[key]]] === 0, nt++];
     If[And @@ Table[Together[e[[5, n]] - refPolar[refKey[key]][[n]]] === 0, {n, 2 gg - 2}], np++]], {e, entries}];
   log[If[nt == Length[entries] && np == Length[entries] && nmiss == 0 && Length[entries] == Length[refTerms], "VALIDATION OK: all " <> ToString[Length[entries]] <> " genus-" <> ToString[gg] <> " decorated graphs and polar parts agree with the stored local run",
     "VALIDATION FAIL: genus-" <> ToString[gg] <> " per-graph agreement: terms " <> ToString[nt] <> ", polar parts " <> ToString[np] <> " of " <> ToString[Length[entries]] <> " (stored " <> ToString[Length[refTerms]] <> ", unmatched keys " <> ToString[nmiss] <> ")"]]]];
log["done"];
