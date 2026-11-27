(* hvc_compare_localP2.wls -- decorated-graph-by-decorated-graph comparison of the conifold polar parts (mode B, t = Log Q at Q = 1)
   with the local P2 mode-B polar parts (t_c at 1 + 27 y = 0) of hodge_vertex/ (same stable-graph enumeration, same orbit keys).
   Moved here from the resolved-conifold project on 2026-11-27 (see README.md); the conifold results are the frozen copies in ref_conifold/.
   Usage:  wolframscript -file hvc_compare_localP2.wls <genus>        Output: hvc_compare_localP2_g<genus>.txt *)
SetDirectory[DirectoryName[$InputFileName]];
gg = ToExpression[$ScriptCommandLine[[2]]];
Get["lp2_tc_section.m"];                     (* tcSetup, polarOfTerm, gapCoefficient of lp2_lib.m (section 6), for the genus-2 orbit terms *)
out = {}; emit[s_] := (AppendTo[out, s]; Print[s]);
(* the conifold results: frozen copies in ref_conifold/ (see ref_conifold/README.md) *)
resFile = First[Select[{"ref_conifold/hvc_results_g" <> ToString[gg] <> "_Bm.m"}, FileExistsQ], None];
If[resFile === None, Print["no hvc_results_g", gg, "_Bm found"]; Exit[1]];
If[StringEndsQ[resFile, ".mx"], Get[resFile]; res = hvcResults, res = Get[resFile]];
If[Head[res] =!= Association, Print["could not load ", resFile]; Exit[1]];
entries = res["entries"];
conMap = Association[Table[{e[[1, 1]], e[[2, 1]], e[[2, 2]]} -> (Together /@ e[[5]]), {e, entries}]];
(* local P2 mode B: genus 2 from the stored orbit terms (t_c expansion), genus 3 from polarOrbits, genus 4 from the cluster file *)
Which[gg == 2,
  r2 = Get["../hodge_vertex_genus2_results.m"]["B"]; tcSetup[8];
  lpMap = Association[Table[Module[{d = r2["descriptions"][id]}, {id[[1]], d["psiOnHalfEdges"], d["dilatonLeaves"]} -> (polarOfTerm[r2["orbitTerms"][id], 2] /. SeriesCoefficient[a_, _] /; FreeQ[a, tcv] :> 0)], {id, Keys[r2["orbitTerms"]]}]],
  gg == 3,
  r3 = Get["../hodge_vertex_genus3_results.m"];
  lpMap = Association[Table[Module[{d = r3["descriptions"][id]}, {id[[1]], d["psiOnHalfEdges"], d["dilatonLeaves"]} -> (Together /@ r3["polarOrbits"][id])], {id, Keys[r3["polarOrbits"]]}]],
  gg == 4,
  r4 = Rest[Get["../hodge_graphterms_genus4.m"]];   (* entries {id, key {psi, dils}, description, term, polar {t_c^-6, ..., t_c^-1}} *)
  lpMap = Association[Table[{e[[1, 1]], e[[2, 1]], e[[2, 2]]} -> (Together /@ e[[5]]), {e, r4}]]];
keys = Keys[conMap];
emit["genus " <> ToString[gg] <> ", mode B (weights eps1+eps2, -eps1, -eps2; R' = R_Z h^-1): " <> ToString[Length[keys]] <> " conifold decorated graphs, " <> ToString[Length[lpMap]] <> " local P2 decorated graphs; key sets identical: " <> ToString[Sort[keys] === Sort[Keys[lpMap]]]];
common = Intersection[keys, Keys[lpMap]];
Do[Module[{n = 2 gg - 2 - (j - 1), both, one, rat},
   both = Select[common, conMap[#][[j]] =!= 0 && lpMap[#][[j]] =!= 0 &];
   one = Select[common, (conMap[#][[j]] === 0) =!= (lpMap[#][[j]] === 0) &];
   rat = Union[Table[Factor[lpMap[k][[j]]/conMap[k][[j]]], {k, both}]];
   emit["order t^-" <> ToString[n] <> " (local P2: t_c^-" <> ToString[n] <> "): both nonzero for " <> ToString[Length[both]] <> " decorated graphs, exactly one nonzero for " <> ToString[Length[one]] <> " " <> ToString[one[[;; Min[6, Length[one]]]]] <> "; distinct ratios localP2/conifold: " <> ToString[Length[rat]] <> " " <> ToString[InputForm[rat[[;; Min[6, Length[rat]]]]]]]], {j, 2 gg - 2}];
(* expected: leading (-3)^(g-1) since t_c^2 = -3 t^2 *)
emit["expected leading ratio (-3)^(g-1) = " <> ToString[(-3)^(gg - 1)]];
(* regularity pattern per decorated graph agrees? *)
regC = Select[common, Union[conMap[#]] === {0} &]; regL = Select[common, Union[lpMap[#]] === {0} &];
emit["regular decorated graphs: conifold " <> ToString[Length[regC]] <> ", local P2 " <> ToString[Length[regL]] <> ", same set: " <> ToString[Sort[regC] === Sort[regL]]];
(* per stable graph polar parts *)
stC = Table[Together /@ Total[Table[conMap[k], {k, Select[common, #[[1]] == g1 &]}]], {g1, Max[common[[All, 1]]]}];
stL = Table[Together /@ Total[Table[lpMap[k], {k, Select[common, #[[1]] == g1 &]}]], {g1, Max[common[[All, 1]]]}];
emit["stable graphs with nonzero polar part: conifold " <> ToString[Select[Range[Length[stC]], Union[stC[[#]]] =!= {0} &]] <> ", local P2 " <> ToString[Select[Range[Length[stL]], Union[stL[[#]]] =!= {0} &]]];
emit["vertex graph (graph 1), decorated graphs {key, conifold polar, local P2 polar, ratios}: "];
Do[emit["  " <> ToString[InputForm[{k[[2 ;;]], conMap[k], lpMap[k], Table[If[conMap[k][[j]] === 0, If[lpMap[k][[j]] === 0, "0/0", "0 vs nonzero"], Factor[lpMap[k][[j]]/conMap[k][[j]]]], {j, 2 gg - 2}]}]]], {k, Sort[Select[common, #[[1]] == 1 &]]}];
Export["hvc_compare_localP2_g" <> ToString[gg] <> ".txt", StringRiffle[out, "\n"], "Text"];
Print["written hvc_compare_localP2_g" <> ToString[gg] <> ".txt"];
