(* rc_polar_collect.m -- collector for the per-decorated-graph conifold data at alpha0 = 0, alpha_i = eps_i.
   Usage:  wolframscript -file rc_polar_collect.m genus
   Input:  rc_polar_x_g<genus>_c*.m and rc_polar_x_partial_g<genus>_c*.m, resolvedconifold_F<genus>0.m (from rc_collect.m, optional),
           ref/conifold_graphterms_genus2_3.m (genus 2, 3 validation).
   Output: conifold_graphterms_genus<genus>.m (polar parts, eps1 = x, eps2 = 1), conifold_graphterms_genus<genus>_full.m (rational
           functions of Q and x, if computed), rc_polar_collect.log with VALIDATION lines. *)
SetDirectory[DirectoryName[$InputFileName]];
gg = ToExpression[$ScriptCommandLine[[2]]];
Get["rc_lib.m"];
log[m_] := Module[{s = OpenAppend["rc_polar_collect.log"]}, WriteString[s, DateString[{"Year", "-", "Month", "-", "Day", " ", "Hour", ":", "Minute", ":", "Second"}], "  [g=", gg, "] ", m, "\n"]; Close[s]; Print[m]];
hom[f_] := Factor[e2^(2 gg - 2) (f /. x -> e1/e2)];
pfiles = FileNames["rc_polar_x_partial_g" <> ToString[gg] <> "_c*.m"];
If[pfiles === {}, log["no partial files found"]; Exit[1]];
parts = Get /@ pfiles;
polar = Total[parts[[All, 1]]]; failed = Join @@ parts[[All, 3]];
exact = ! MemberQ[parts[[All, 2]], Null];
rqtot = If[exact, Together[Total[parts[[All, 2]]]], Null];
log[ToString[Length[pfiles]] <> " chunk files; failed orbits: " <> ToString[Length[failed]] <> If[failed =!= {}, " " <> ToString[failed], ""] <> "; full rational functions: " <> ToString[exact]];
(* total polar part = polar part of P_g Li_{3-2g}(Q)/(eps1 eps2) *)
expected = cyPolar[gg] /. {e1 -> x, e2 -> 1};
gapok = Union[Together /@ (polar - expected)] === {0};
log["total polar part: " <> ToString[InputForm[Factor /@ polar]] <> "; expected " <> ToString[InputForm[Factor /@ expected]]];
log[If[gapok, "VALIDATION OK: sum of the polar parts = polar part of P_g Li_{3-2g}", "VALIDATION FAIL: total polar part"]];
(* total rational function = F<g>0 of the general torus at alpha0 = 0, alpha_i = eps_i (eps1 = x, eps2 = 1) *)
If[exact && FileExistsQ["resolvedconifold_F" <> ToString[gg] <> "0.m"],
  Module[{Fref = Last[Get["resolvedconifold_F" <> ToString[gg] <> "0.m"]] /. {a0 -> 0, a1 -> x, a2 -> 1, e1 -> x, e2 -> 1}, ok},
   ok = Together[rqtot - Fref] === 0;
   log[If[ok, "VALIDATION OK: sum of the decorated graphs = F" <> ToString[gg] <> "0 at alpha0 = 0, alpha_i = eps_i", "VALIDATION FAIL: sum of the decorated graphs differs from F" <> ToString[gg] <> "0 specialised"]]],
  log["(no comparison of the full sum with F" <> ToString[gg] <> "0: " <> If[exact, "resolvedconifold_F" <> ToString[gg] <> "0.m not found", "series method, no rational functions"] <> ")"]];
(* per-graph files *)
tfiles = FileNames["rc_polar_x_g" <> ToString[gg] <> "_c*.m"];
entries = SortBy[Join @@ (Get /@ tfiles), First];
readme = "Resolved conifold x C2 at alpha0 = 0, alpha1 = eps1, alpha2 = eps2, genus " <> ToString[gg] <> ": the " <> ToString[Length[entries]] <> " decorated stable graphs of the Givental-Teleman sum (same enumeration, keys and descriptions as the local P2 files). " <>
  "Entry: {id {graph, orbit}, key, description, polar part at Q = 1 in t = Log Q with eps1 = x, eps2 = 1: coefficients of t^-(2g-2), ..., t^-1 (restore eps2 by homogeneity of degree 2g-2)}. Weight 1/|Aut| times the sum over labelled decorations and colourings.";
Put[Prepend[Table[e[[1 ;; 4]], {e, entries}], readme], "conifold_graphterms_genus" <> ToString[gg] <> ".m"];
If[exact, Put[Prepend[Table[{e[[1]], e[[2]], e[[5]]}, {e, entries}], readme <> " This file: {id, key, contribution as a rational function of Q and x = eps1/eps2 with eps2 = 1}."], "conifold_graphterms_genus" <> ToString[gg] <> "_full.m"]];
log["wrote conifold_graphterms_genus" <> ToString[gg] <> ".m (" <> ToString[Length[entries]] <> " entries)" <> If[exact, " and the _full.m file", ""]];
(* validation against the stored genus-2/3 per-graph polar parts *)
If[gg <= 3 && FileExistsQ["ref/conifold_graphterms_genus2_3.m"],
  Module[{ref = Get["ref/conifold_graphterms_genus2_3.m"], refE, refA, np = 0, nf = 0},
   refE = gg /. ref[[2]];
   refA = Association[Table[r[[1]] -> r, {r, refE}]];
   Do[If[KeyExistsQ[refA, e[[1]]],
      If[And @@ Table[Together[hom[e[[4, n]]] - refA[e[[1]]][[4, n]]] === 0, {n, 2 gg - 2}], np++];
      If[exact && Together[hom[e[[5]]] - refA[e[[1]]][[3]]] === 0, nf++]], {e, entries}];
   log[If[np == Length[refE] && Length[entries] == Length[refE] && (! exact || nf == Length[refE]), "VALIDATION OK: all " <> ToString[Length[refE]] <> " genus-" <> ToString[gg] <> " decorated graphs agree with the stored results (polar parts" <> If[exact, " and rational functions)", ")"],
     "VALIDATION FAIL: genus-" <> ToString[gg] <> " per-graph agreement: polar " <> ToString[np] <> ", full " <> ToString[nf] <> " of " <> ToString[Length[refE]] <> " (computed " <> ToString[Length[entries]] <> ")"]]]];
log["done"];
