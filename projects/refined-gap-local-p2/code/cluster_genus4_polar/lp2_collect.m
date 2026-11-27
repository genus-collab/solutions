(* lp2_collect.m -- collector for the local P2 per-decorated-graph pipeline.
   Usage:  wolframscript -file lp2_collect.m genus
   Input:  lp2_partial_g<genus>_c*.m and lp2_terms_g<genus>_c*.m (all chunks), lp2_orbits_g<genus>.mx,
           reference totals ref/F20_F30.m (genus 2, 3) or ref/genRmat_F40_cluster_result.m (genus 4, q -> 0),
           reference per-graph data ref/localP2_graphterms_genus3.m (genus 3 validation).
   Output: localP2_F<genus>0.m (the total), localP2_graphterms_genus<genus>.m (all decorated graphs with polar parts),
           localP2_polar_genus<genus>.m (polar parts only), lp2_collect.log with VALIDATION lines. *)
SetDirectory[DirectoryName[$InputFileName]];
gg = ToExpression[$ScriptCommandLine[[2]]];
Get["lp2_lib.m"];
log[m_] := Module[{s = OpenAppend["lp2_collect.log"]}, WriteString[s, DateString[{"Year", "-", "Month", "-", "Day", " ", "Hour", ":", "Minute", ":", "Second"}], "  [g=", gg, "] ", m, "\n"]; Close[s]; Print[m]];
e1 = \[Epsilon][1]; e2 = \[Epsilon][2];
loadOrbits[gg];
pfiles = FileNames["lp2_partial_g" <> ToString[gg] <> "_c*.m"];
If[pfiles === {}, log["no partial files found"]; Exit[1]];
parts = Get /@ pfiles;
F = Together[Total[parts[[All, 1]]]]; polar = Total[parts[[All, 2]]]; failed = Join @@ parts[[All, 3]];
log[ToString[Length[pfiles]] <> " chunk files; failed reconstructions: " <> ToString[Length[failed]] <> If[failed =!= {}, " " <> ToString[failed], ""]];
Put[F, "localP2_F" <> ToString[gg] <> "0.m"];
(* ---- total against the reference ---- *)
ref = If[gg <= 3, Get["ref/F20_F30_y.m"][[gg - 1]], Get["ref/genRmat_F40_cluster_result.m"][[1]] /. {q0 -> 0, q1 -> 0, q2 -> 0, q3 -> 0, q4 -> 0}];   (* F20, F30 in y, XX, eps *)
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
tfiles = FileNames["lp2_terms_g" <> ToString[gg] <> "_c*.m"];
entries = SortBy[Join @@ (Get /@ tfiles), First];
readme = "Local P2 x C2, fibrewise torus, genus " <> ToString[gg] <> ": the " <> ToString[Length[entries]] <> " decorated stable graphs of the Givental-Teleman sum for F_{" <> ToString[gg] <> ",0}. " <>
  "Entry: {id {graph, orbit}, key {psi powers on the two half-edges of each edge, dilaton-leaf psi powers per vertex}, description (genera, legs, adjacency with loops on the diagonal, edges, Aut, dilatonFree), " <>
  "contribution as a function of y, XX = y c1'(y)/c1(y), eps1, eps2 (conventions of localp2CY5_genus3.nb; weight 1/|Aut| times the sum over labelled decorations and colourings), " <>
  "polar part in the conifold coordinate t_c = Delta + 11/18 Delta^2 + ... (Delta = 1 + 27 y) with eps1 = x, eps2 = 1: coefficients of t_c^-(2g-2), ..., t_c^-1 (homogeneous of degree 2g-2 in eps after restoring eps2). The contributions sum to F_{g,0}.";
Put[Prepend[entries, readme], "localP2_graphterms_genus" <> ToString[gg] <> ".m"];
Put[Prepend[Table[{e[[1]], e[[2]], e[[3]], e[[5]]}, {e, entries}], readme], "localP2_polar_genus" <> ToString[gg] <> ".m"];
log["wrote localP2_graphterms_genus" <> ToString[gg] <> ".m (" <> ToString[Length[entries]] <> " entries) and localP2_polar_genus" <> ToString[gg] <> ".m"];
(* ---- genus-3 validation against the stored per-graph data ---- *)
If[gg == 3 && FileExistsQ["ref/localP2_graphterms_genus3.m"],
  Module[{refE = Rest[Get["ref/localP2_graphterms_genus3.m"]], refA, nt = 0, np = 0, hom},
   hom[f_] := Factor[e2^(2 gg - 2) (f /. x -> e1/e2)];
   refA = Association[Table[r[[1]] -> r, {r, refE}]];
   Do[If[KeyExistsQ[refA, e[[1]]],
      If[Together[e[[4]] - refA[e[[1]]][[3]]] === 0, nt++];
      If[And @@ Table[Together[hom[e[[5, n]]] - refA[e[[1]]][[4, n]]] === 0, {n, 2 gg - 2}], np++]], {e, entries}];
   log[If[nt == Length[refE] && np == Length[refE] && Length[entries] == Length[refE], "VALIDATION OK: all " <> ToString[Length[refE]] <> " genus-3 decorated graphs and polar parts agree with the stored results",
     "VALIDATION FAIL: genus-3 per-graph agreement: terms " <> ToString[nt] <> ", polar parts " <> ToString[np] <> " of " <> ToString[Length[refE]] <> " (computed " <> ToString[Length[entries]] <> ")"]]]];
log["done"];
