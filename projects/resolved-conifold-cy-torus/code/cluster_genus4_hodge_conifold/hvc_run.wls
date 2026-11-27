(* hvc_run.wls -- driver for the Hodge-vertex graph sum of the resolved conifold x C^2 (fibrewise torus, alpha_i = eps_i).
   Usage:
     wolframscript -file hvc_run.wls checks                          Hodge-integral checks, PPZ scalar check, rank-4 class triviality
     wolframscript -file hvc_run.wls g2sign                          genus 2, all (mode, s) combinations: which reproduce F20
     wolframscript -file hvc_run.wls prep <g> <mode> <s>             factors (Q, x symbolic) -> hvc_factors_g<g>_<tag>.mx; orbits -> hvc_orbits_g<g>_*.mx
     wolframscript -file hvc_run.wls orbits <g> <mode> <s> <c> <n>   decorated graphs c, c+n, ... -> hvc_terms_g<g>_<tag>_c<c>of<n>.m
     wolframscript -file hvc_run.wls collect <g> <mode> <s>          totals, gap, per stable graph -> hvc_results_g<g>_<tag>.m, hvc_table_g<g>_<tag>.txt
   mode = B (weights eps1+eps2, -eps1, -eps2), A (weights -eps1-eps2, eps1, eps2) or point; s = m (h_w^-1) or p (h_w^+1), the exponent
   of h_w (a numeric "-1" would be swallowed by wolframscript as one of its own options); tag = B/A + m/p.
   Log: hvc.log *)
SetDirectory[DirectoryName[$InputFileName]];
args = Rest[$ScriptCommandLine];
stage = If[Length[args] >= 1, args[[1]], "checks"];
Get["hvc_lib.m"];
sOf[mode_] := If[mode == "point", 0, Switch[args[[4]], "m" | "-1", -1, "p" | "1" | "+1", 1, _, Print["bad sign argument ", args[[4]], " (use m or p)"]; Exit[1]]];
log[m_] := Module[{st = OpenAppend["hvc.log"]}, WriteString[st, DateString[{"Year", "-", "Month", "-", "Day", " ", "Hour", ":", "Minute", ":", "Second"}], "  [", stage, If[Length[args] >= 2, " " <> StringRiffle[args[[2 ;;]], " "], ""], "] ", m, "\n"]; Close[st]; Print[m]];
fmt[expr_] := ToString[InputForm[expr /. {e1s -> E1, e2s -> E2}]];
tt = AbsoluteTime[];

Which[
 stage == "checks",
 Module[{ok, rep, w3 = {w1, w2, w3}, ppz2, exp2, ppz3, exp3, wnum, cl4, wB, wA},
  {ok, rep} = hodgeChecks4[]; Do[log[r], {r, rep}];
  ppz2 = ppzScalar[2, w3]; exp2 = hodgeNormU4[2, {}] /. usubOf[w3];
  log["PPZ scalar R-action (rows h_w^-1 on point vertices), genus 2, general weights: " <> ToString[InputForm[Factor[ppz2]]] <> " = normalised Hodge class on M_2: " <> ToString[Together[ppz2 - exp2] === 0]];
  SeedRandom[7];
  Do[wnum = Table[RandomInteger[{-9, 9}]/RandomInteger[{1, 5}], {3}]; If[MemberQ[wnum, 0], Continue[]];
   ppz3 = ppzScalar[3, wnum]; exp3 = hodgeNormU4[3, {}] /. usubOf[wnum];
   log["PPZ scalar R-action, genus 3 (42 graphs), weights " <> ToString[wnum] <> ": " <> ToString[ppz3] <> " = normalised Hodge class on M_3 " <> ToString[exp3] <> ": " <> ToString[ppz3 == exp3]]; ok = ok && ppz3 == exp3, {2}];
  (* the rank-4 class of the actual normal bundle O(-1)+O(-1)+O+O with weights (eps1, eps2, -eps1, -eps2) is trivial *)
  cl4 = Table[{g, Table[{c[[1]], Factor[c[[2]] /. {u[1] -> 1/e1s, u[2] -> 1/e2s, u[3] -> -1/e1s, u[4] -> -1/e2s}]}, {c, hodgeClassN[g, 4]}]}, {g, 1, 4}];
  log["rank-4 Hodge class Lambda^vee(eps1) Lambda^vee(eps2) Lambda^vee(-eps1) Lambda^vee(-eps2)/(eps1 eps2)^(2g), lambda monomials {k1, k3} with coefficients, genus 1-4: " <> fmt[cl4] <> "; trivial: " <> ToString[And @@ Table[DeleteCases[c[[2]], {_, 0}] === {{{0, 0}, 1}}, {c, cl4}]]];
  wB = modeWeights["B"]; wA = modeWeights["A"];
  log["h_wB(z) to z^6: " <> fmt[Factor /@ hCoeffs[wB, 1, 6]]];
  log["h_wA(z) = h_wB(z)^-1: " <> ToString[Union[Together /@ (hCoeffs[wA, 1, 6] - hCoeffs[wB, -1, 6])] === {0}]];
  log["normalised triple Hodge class, mode B, genus 2 and 3 (lambda monomials {k1, k3}): " <> fmt[Table[{g, Table[{c[[1]], Factor[c[[2]] /. usubOf[wB]]}, {c, hodgeClass4[g]}]}, {g, 2, 3}]]];
  log[If[ok, "CHECKS OK (" <> ToString[Round[AbsoluteTime[] - tt]] <> " s)", "CHECKS FAILED"]]],

 stage == "g2sign",
 Module[{gg = 2, Fref, res = {}, nfo, terms, tot, diff, pol, gap},
  Fref = FrefHVC[gg, x];
  Do[Module[{mode = ms[[1]], s = ms[[2]], t0 = AbsoluteTime[]},
    buildFactors[3, mode, s, x];
    nfo = loadOrbitsHVC[gg, mode];
    terms = Table[orbitTerm[OT["graphs"][id], usubMode[mode, x], gg], {id, OT["ids"]}];
    tot = Together[Total[terms[[All, 1]]]]; diff = Together[tot - Fref];
    pol = Together /@ Total[terms[[All, 2]]]; gap = gapHVC[gg, x];
    log["mode " <> mode <> ", s = " <> ToString[s] <> ": " <> ToString[nfo] <> " (graphs, labelled decorations, decorated graphs); all values symmetric: " <> ToString[And @@ terms[[All, 3]]] <>
      "; sum of decorated graphs - F20 = " <> ToString[InputForm[diff]] <> "; total polar part " <> ToString[InputForm[Factor /@ pol]] <> " (gap " <> ToString[InputForm[gap]] <> "); " <> ToString[Round[AbsoluteTime[] - t0, 0.1]] <> " s"];
    AppendTo[res, {mode, s, diff === 0}]], {ms, {{"point", 0}, {"B", -1}, {"B", 1}, {"A", -1}, {"A", 1}}}];
  log["genus-2 sign test: {mode, s, reproduces F20}: " <> ToString[res]]],

 stage == "prep",
 Module[{gg = ToExpression[args[[2]]], mode = args[[3]], s, tag, nfo, t0 = AbsoluteTime[]},
  s = sOf[mode]; tag = modeTag[mode, s];
  buildFactors[3 gg - 3, mode, s, x];
  log["factors for genus " <> ToString[gg] <> ", mode " <> mode <> ", s = " <> ToString[s] <> " built to z^" <> ToString[ZMAX] <> " (" <> ToString[Round[AbsoluteTime[] - t0, 0.1]] <> " s); LeafCount of EdN[1,2][" <> ToString[ZMAX - 1] <> ",0]: " <> ToString[LeafCount[EdN[1, 2][ZMAX - 1, 0]]]];
  DumpSave["hvc_factors_g" <> ToString[gg] <> "_" <> tag <> ".mx", {rN, DeltaN, DilN, LfN, EdN, ratCN, DeltaCN, Ssum, Pprod, SsumQ, PprodQ, etainvCur, genSym, ZMAX}];
  t0 = AbsoluteTime[]; nfo = loadOrbitsHVC[gg, mode];
  log["decorated graphs (" <> If[mode == "point", "point", "Hodge"] <> " vertices): " <> ToString[nfo] <> " (graphs, labelled decorations, decorated graphs), " <> ToString[Round[AbsoluteTime[] - t0, 0.1]] <> " s; files hvc_factors_g" <> ToString[gg] <> "_" <> tag <> ".mx, " <> orbitFile[gg, mode]]],

 stage == "orbits",
 (* restartable: every finished decorated graph is appended at once to the chunk file; on a restart the ids already present are skipped *)
 Module[{gg = ToExpression[args[[2]]], mode = args[[3]], s, tag, cc, nch, idx, usub, nonsym = 0, outfile, done, todo, ndone = 0},
  s = sOf[mode]; tag = modeTag[mode, s];
  {cc, nch} = ToExpression /@ args[[5 ;; 6]];
  Get["hvc_factors_g" <> ToString[gg] <> "_" <> tag <> ".mx"]; loadOrbitsHVC[gg, mode];
  usub = usubMode[mode, x];
  idx = Select[Range[Length[OT["ids"]]], Mod[# - 1, nch] + 1 == cc &];
  outfile = "hvc_terms_g" <> ToString[gg] <> "_" <> tag <> "_c" <> ToString[cc] <> "of" <> ToString[nch] <> ".m";
  done = If[FileExistsQ[outfile], Module[{r = readChunk[outfile]}, If[r === $Failed, Module[{bad = outfile <> ".corrupt." <> DateString[{"Hour", "Minute", "Second"}]}, RenameFile[outfile, bad]; log["chunk file " <> outfile <> " unreadable (torn record?), moved to " <> bad <> "; recomputing the chunk"]]; {}, r[[All, 1]]]], {}];
  todo = Select[idx, ! MemberQ[done, OT["ids"][[#]]] &];
  log["chunk " <> ToString[cc] <> " of " <> ToString[nch] <> ": " <> ToString[Length[idx]] <> " decorated graphs, " <> ToString[Length[idx] - Length[todo]] <> " already in " <> outfile <> ", " <> ToString[Length[todo]] <> " to do; ZMAX " <> ToString[ZMAX]];
  Do[Module[{id = OT["ids"][[i]], r, t0 = AbsoluteTime[]},
    r = orbitTerm[OT["graphs"][id], usub, gg];
    If[! r[[3]], nonsym++];
    PutAppend[{id, OT["orbits"][[id[[1]], id[[2]]]], Normal[OT["descriptions"][id]], r[[1]], r[[2]], r[[3]]}, outfile]; ndone++;
    log["orbit " <> ToString[id] <> " " <> ToString[OT["orbits"][[id[[1]], id[[2]]]]] <> " (" <> ToString[Round[AbsoluteTime[] - t0, 0.1]] <> " s): polar " <> ToString[InputForm[r[[2]]]] <> If[r[[3]], "", "  NON-SYMMETRIC"]]], {i, todo}];
  log["finished: " <> ToString[ndone] <> " decorated graphs computed now (" <> ToString[Length[idx]] <> " in the chunk), " <> ToString[nonsym] <> " non-symmetric, " <> ToString[Round[AbsoluteTime[] - tt]] <> " s"]],

 stage == "collect",
 (* optional 5th argument: a subdirectory holding the chunk files (e.g. stanage_genus4, copied back from the cluster); the results,
    the table and the log lines are then written there as well *)
 Module[{gg = ToExpression[args[[2]]], mode = args[[3]], s, tag, files, entries, ngr, stable, total, Fref, diff, numchk, pol, gap, gapok, regQ, regular, tab, ref, refKey, nt, np, nmiss, ords, results, dir},
  s = sOf[mode]; tag = modeTag[mode, s];
  dir = If[Length[args] >= 5, args[[5]], "."];
  loadOrbitsHVC[gg, mode];
  files = FileNames["hvc_terms_g" <> ToString[gg] <> "_" <> tag <> "_c*of*.m", dir];
  If[dir =!= ".", log["chunk files read from " <> dir <> "; outputs written there"]];
  (* chunk files hold either one list of records (old Put format) or one record per expression (PutAppend format) *)
  entries = readChunk /@ files;
  If[MemberQ[entries, $Failed], log["unreadable chunk files (rerun their orbits stage, which moves them aside and recomputes): " <> ToString[Pick[files, entries, $Failed]]]; entries = DeleteCases[entries, $Failed]];
  entries = SortBy[Join @@ entries, First] /. SeriesCoefficient[a_, _] /; FreeQ[a, t] :> 0;
  entries = DeleteDuplicatesBy[entries, First];   (* a chunk restarted after a crash during PutAppend could hold an id twice *)
  log[ToString[Length[files]] <> " chunk files, " <> ToString[Length[entries]] <> " decorated graphs of " <> ToString[Length[OT["ids"]]] <> "; non-symmetric values: " <> ToString[Count[entries, e_ /; ! e[[6]]]]];
  If[Length[entries] != Length[OT["ids"]] || Union[entries[[All, 1]]] =!= Union[OT["ids"]], log["INCOMPLETE -- missing ids: " <> ToString[Complement[OT["ids"], entries[[All, 1]]][[;; Min[20, Length[Complement[OT["ids"], entries[[All, 1]]]]]]]] <> " -- abort"]; Exit[1]];
  ngr = Length[prepared];
  regQ[v_] := Union[Together /@ v] === {0};
  stable = Table[Module[{es = Select[entries, #[[1, 1]] == k &]}, <|"graph" -> k, "genera" -> prepared[[k]]["g"], "edges" -> prepared[[k]]["edges"], "Aut" -> prepared[[k]]["aut"], "decoratedGraphs" -> Length[es],
       "singularDecoratedGraphs" -> Count[es, e_ /; ! regQ[e[[5]]]], "term" -> Together[Total[es[[All, 4]]]], "polar" -> (Together /@ Total[es[[All, 5]]])|>], {k, ngr}];
  total = Together[Total[stable[[All, "term"]]]]; Fref = FrefHVC[gg, x];
  diff = TimeConstrained[Together[total - Fref], 3600, "timeout"];
  SeedRandom[4]; numchk = Table[Together[(total - Fref) /. {Q -> RandomInteger[{1, 50}]/RandomInteger[{51, 97}], x -> RandomInteger[{2, 40}]/RandomInteger[{1, 7}]}], {5}];
  log["sum of decorated graphs minus F" <> ToString[gg] <> "0 [reference: " <> FrefKind[gg] <> "] (symbolic in Q, x): " <> ToString[InputForm[diff]] <> "; at 5 random rational points: " <> ToString[numchk]];
  log[If[diff === 0, "VALIDATION OK: sum of decorated graphs = F" <> ToString[gg] <> "0", "VALIDATION FAIL: total differs from the reference"]];
  (* independent of the constant-map term: total - P_g Li_{3-2g}(Q)/(eps1 eps2) must be Q-independent, and its value is F_g(0) *)
  Module[{c0 = Together[total - (Pg[gg] PolyLog[3 - 2 gg, Q]/(e1 e2) /. {e1 -> x, e2 -> 1})]},
   log["total minus P_g Li_{3-2g}(Q)/(eps1 eps2) is Q-independent: " <> ToString[Together[D[c0, Q]] === 0] <> "; its value (constant-map term F_g(0)): " <> ToString[InputForm[Factor[c0 /. Q -> 2]]] <> "; closed form -(eps1^(2g-2)+eps2^(2g-2))|B_2g B_2g-2|/(4g(2g-2)(2g-2)!): " <> ToString[InputForm[Factor[F0closed[gg, x]]]]]];
  pol = Together /@ Total[stable[[All, "polar"]]]; gap = gapHVC[gg, x];
  gapok = Union[Together /@ (pol - gap)] === {0};
  log["total polar part {t^-" <> ToString[2 gg - 2] <> ", ..., t^-1} (eps1 = x, eps2 = 1): " <> ToString[InputForm[Factor /@ pol]] <> "; expected " <> ToString[InputForm[gap]] <> "; " <> If[gapok, "VALIDATION OK: gap", "VALIDATION FAIL: gap"]];
  regular = Select[Range[ngr], regQ[stable[[#]]["polar"]] &];
  log[ToString[ngr] <> " stable graphs; regular at Q = 1: " <> ToString[Length[regular]] <> " -> " <> ToString[regular] <> "; singular: " <> ToString[Complement[Range[ngr], regular]]];
  log["singular stable graphs {k, genera, edges, polar part}: " <> ToString[InputForm[Table[{k, stable[[k]]["genera"], stable[[k]]["edges"], Factor /@ stable[[k]]["polar"]}, {k, Complement[Range[ngr], regular]}]]]];
  log["stable graphs containing singular decorated graphs {k, genera, singular, total}: " <> ToString[Table[{k, stable[[k]]["genera"], stable[[k]]["singularDecoratedGraphs"], stable[[k]]["decoratedGraphs"]}, {k, Select[Range[ngr], stable[[#]]["singularDecoratedGraphs"] > 0 &]}]] <> "; graphs with a vertex of genus >= 2: " <> ToString[Select[Range[ngr], Max[prepared[[#]]["g"]] >= 2 &]]];
  ords = Table[Module[{p = Position[Together /@ e[[5]], _?(# =!= 0 &), {1}, 1, Heads -> False]}, If[p === {}, 0, 2 gg - 2 - (p[[1, 1]] - 1)]], {e, entries}];
  log["pole-order tally over decorated graphs {order, count}: " <> ToString[Sort[Tally[ords]]] <> "; zero terms: " <> ToString[Count[entries, e_ /; Together[e[[4]]] === 0]]];
  tab = StringRiffle[Prepend[Table[StringRiffle[{k, stable[[k]]["genera"], stable[[k]]["edges"], stable[[k]]["Aut"], stable[[k]]["decoratedGraphs"], stable[[k]]["singularDecoratedGraphs"], If[regQ[stable[[k]]["polar"]], "regular", ToString[InputForm[Factor /@ stable[[k]]["polar"]]]]}, "  "], {k, ngr}],
     "graph  genera  edges  |Aut|  decorated graphs  singular decorated graphs  polar part {t^-(2g-2), ..., t^-1} at Q = 1, t = Log Q (eps1 = x, eps2 = 1)"], "\n"];
  Export[FileNameJoin[{dir, "hvc_table_g" <> ToString[gg] <> "_" <> tag <> ".txt"}], tab, "Text"];
  (* point mode: regression against the stored per-graph terms at alpha0 = 0, alpha_i = eps_i *)
  If[mode == "point" && FileExistsQ["ref/conifold_graphterms_genus2_3.m"],
   ref = Association[Get["ref/conifold_graphterms_genus2_3.m"][[2]]][gg];
   refKey = Association[Table[{e[[1, 1]], "psiOnHalfEdges" /. e[[2]], "dilatonLeaves" /. e[[2]]} -> e, {e, ref}]];
   nt = 0; np = 0; nmiss = 0;
   Do[Module[{key = {e[[1, 1]], e[[2, 1]], e[[2, 2]]}, r},
     If[! KeyExistsQ[refKey, key], nmiss++, r = refKey[key];
      If[Together[e[[4]] - (r[[3]] /. {e1 -> x, e2 -> 1})] === 0, nt++];
      If[Union[Together /@ (e[[5]] - (r[[4]] /. {e1 -> x, e2 -> 1}))] === {0}, np++]]], {e, entries}];
   log[If[nt == Length[entries] && np == Length[entries] && nmiss == 0 && Length[ref] == Length[entries], "VALIDATION OK: all " <> ToString[Length[entries]] <> " point-vertex decorated graphs and polar parts agree with conifold_graphterms_genus2_3.m",
     "VALIDATION FAIL: point-vertex regression: terms " <> ToString[nt] <> ", polar " <> ToString[np] <> " of " <> ToString[Length[entries]] <> " (stored " <> ToString[Length[ref]] <> ", unmatched " <> ToString[nmiss] <> ")"]]];
  results = <|"genus" -> gg, "mode" -> mode, "s" -> s, "weights" -> If[mode == "point", None, modeWeights[mode]], "stable" -> stable, "entries" -> entries, "total" -> total, "totalOK" -> (diff === 0), "polar" -> pol, "gapOK" -> gapok,
    "readme" -> "Resolved conifold x C2, fibrewise torus alpha0 = 0, alpha1 = eps1 = x, alpha2 = eps2 = 1, genus " <> ToString[gg] <> ", mode " <> mode <> If[mode == "point", " (Witten-Kontsevich vertices, R_Z)", " (triple Hodge vertex with weights " <> ToString[InputForm[modeWeights[mode]]] <> ", R' = R_Z h_w^(" <> ToString[s] <> "))"] <>
      ". stable: per stable graph (order of StableGraph-" <> ToString[gg] <> ".bin) genera, edges, |Aut|, numbers of decorated graphs, term (rational function of Q and x), polar part {t^-(2g-2), ..., t^-1} at Q = 1 in t = Log Q. entries: per decorated graph {id {graph, orbit}, key {psi on half-edges, dilaton leaves}, description, term, polar, symmetric}. Weight convention 1/|Aut| times the sum over labelled decorations and the 2^V colourings."|>;
  (* genus >= 4: the full results (with the 27172 decorated-graph terms) go to a fast binary sidecar .mx (DumpSave; same Mathematica
     version and platform), and the portable .m file omits the per-decorated-graph terms (they stay in the chunk files); writing the
     full .m in InputForm would take hours *)
  If[gg >= 4,
   hvcResults = results;   (* a global symbol: DumpSave of the Module-local `results` would store it under its localised name *)
   DumpSave[FileNameJoin[{dir, "hvc_results_g" <> ToString[gg] <> "_" <> tag <> ".mx"}], hvcResults];
   Put[Append[KeyDrop[results, "entries"], "entries" -> (ReplacePart[#, 4 -> "term in the chunk files / the .mx sidecar"] & /@ entries)], FileNameJoin[{dir, "hvc_results_g" <> ToString[gg] <> "_" <> tag <> ".m"}]];
   log["written " <> FileNameJoin[{dir, "hvc_results_g" <> ToString[gg] <> "_" <> tag <> ".mx"}] <> " (full) and the .m file without the decorated-graph terms, and hvc_table_g" <> ToString[gg] <> "_" <> tag <> ".txt (" <> ToString[Round[AbsoluteTime[] - tt]] <> " s)"],
   Put[results, FileNameJoin[{dir, "hvc_results_g" <> ToString[gg] <> "_" <> tag <> ".m"}]];
   log["written " <> FileNameJoin[{dir, "hvc_results_g" <> ToString[gg] <> "_" <> tag <> ".m"}] <> " and hvc_table_g" <> ToString[gg] <> "_" <> tag <> ".txt (" <> ToString[Round[AbsoluteTime[] - tt]] <> " s)"]]],

 True, Print["unknown stage ", stage]];
