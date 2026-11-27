(* hvc_analysis.wls -- structural analysis of the per-decorated-graph results of one (genus, mode) run of the resolved conifold
   Hodge-vertex graph sum (same questions as hv_analysis_g4.wls for local P2):
     which stable graphs are regular at Q = 1; does the genus-g vertex graph carry exactly the gap; which decorated graphs are
     singular and how this correlates with the vertex genera; pole-order bound Sum_v max(2 g_v - 2, 0); odd orders; the
     structure of the terms as rational functions of Q and x (poles at Q = 1, Q = 0 and at the discriminant point
     Q_c = -4x/(1-x)^2 where the two idempotents collide; x-denominators).
   Usage:  wolframscript -file hvc_analysis.wls <genus> <tag>      (tag = Bm, Am, point)
   Output: hvc_analysis_g<genus>_<tag>.txt (and the same lines on stdout). *)
SetDirectory[DirectoryName[$InputFileName]];
{gg, tag} = {ToExpression[$ScriptCommandLine[[2]]], $ScriptCommandLine[[3]]};
Get["hvc_lib.m"];
(* results in this folder, or (copied back from the cluster) in stanage_genus4/; the full binary sidecar .mx (genus 4) is preferred over
   the .m file, which omits the decorated-graph terms; the output goes next to the results file *)
resBase = "hvc_results_g" <> ToString[gg] <> "_" <> tag;
resDir = If[FileExistsQ[resBase <> ".m"] || FileExistsQ[resBase <> ".mx"], ".", "stanage_genus4"];
If[FileExistsQ[FileNameJoin[{resDir, resBase <> ".mx"}]], Get[FileNameJoin[{resDir, resBase <> ".mx"}]]; res = hvcResults, res = Get[FileNameJoin[{resDir, resBase <> ".m"}]]];
If[Head[res] =!= Association, Print["could not load ", resBase]; Exit[1]];
stable = res["stable"]; entries = res["entries"]; ngr = Length[stable];
out = {}; emit[s_] := (AppendTo[out, s]; Print[s]);
regQ[v_] := Union[Together /@ v] === {0};
emit[res["readme"]];
emit[""];
(* ---- totals ---- *)
gap = gapHVC[gg, x];
emit["total = F" <> ToString[gg] <> "0 (symbolic in Q, x): " <> ToString[res["totalOK"]] <> "; total polar part {t^-" <> ToString[2 gg - 2] <> ", ..., t^-1} = gap: " <> ToString[res["gapOK"]] <> "; gap = " <> ToString[InputForm[gap]]];
(* ---- stable graphs ---- *)
regular = Select[Range[ngr], regQ[stable[[#]]["polar"]] &]; singular = Complement[Range[ngr], regular];
emit[ToString[ngr] <> " stable graphs; regular at Q = 1: " <> ToString[Length[regular]] <> "; singular: " <> ToString[Table[{k, stable[[k]]["genera"], stable[[k]]["edges"]}, {k, singular}]]];
emit["genus-" <> ToString[gg] <> " vertex graph (graph 1): polar part = " <> ToString[InputForm[Factor /@ stable[[1]]["polar"]]] <> "; equal to the gap: " <> ToString[Union[Together /@ (stable[[1]]["polar"] - gap)] === {0}] <> "; polar part of the sum of all other graphs: " <> ToString[InputForm[Factor /@ Together /@ Total[stable[[2 ;;, "polar"]]]]]];
(* ---- decorated graphs ---- *)
sing = Table[stable[[k]]["singularDecoratedGraphs"], {k, ngr}];
maxg = Table[Max[stable[[k]]["genera"]], {k, ngr}];
setG = Select[Range[ngr], maxg[[#]] >= 2 &]; setS = Select[Range[ngr], sing[[#]] > 0 &];
emit[ToString[Length[entries]] <> " decorated graphs, " <> ToString[Total[sing]] <> " singular individually; stable graphs containing singular decorated graphs: " <> ToString[setS] <> "; stable graphs with a vertex of genus >= 2: " <> ToString[setG] <> "; equal: " <> ToString[setS === setG]];
emit["per stable graph {k, genera, decorated graphs, singular}: " <> ToString[Table[{k, stable[[k]]["genera"], stable[[k]]["decoratedGraphs"], sing[[k]]}, {k, Select[Range[ngr], sing[[#]] > 0 || maxg[[#]] >= 2 &]}]]];
leadOrder[p_] := Module[{pos = Position[Together /@ p, _?(# =!= 0 &), {1}, 1, Heads -> False]}, If[pos === {}, 0, 2 gg - 2 - (pos[[1, 1]] - 1)]];
lead = leadOrder /@ entries[[All, 5]];
sumBound = Table[Total[Max[2 # - 2, 0] & /@ stable[[e[[1, 1]]]]["genera"]], {e, entries}];
emit["pole-order tally over decorated graphs {order, count}: " <> ToString[Sort[Tally[lead]]] <> "; bound order <= Sum_v max(2 g_v - 2, 0) holds for all: " <> ToString[And @@ Thread[lead <= sumBound]] <> "; tally of {order, bound}: " <> ToString[Sort[Tally[Transpose[{lead, sumBound}]]]]];
If[setG =!= {}, emit["bound attained {graph, max order among its decorated graphs, bound}: " <> ToString[Table[{k, Max[Pick[lead, entries[[All, 1, 1]], k]], Total[Max[2 # - 2, 0] & /@ stable[[k]]["genera"]]}, {k, setG}]]]];
emit["odd pole orders occur in single decorated graphs: " <> ToString[Count[lead, _?OddQ]] <> " (t^-(2g-3), t^-(2g-5), ...: " <> ToString[Table[Count[lead, 2 gg - 2 - j], {j, 1, 2 gg - 3, 2}]] <> ")"];
nonzeroCoeffs = Table[Count[entries, e_ /; Together[e[[5, n]]] =!= 0], {n, 2 gg - 2}];
emit["number of decorated graphs with a nonzero coefficient at t^-" <> ToString[Range[2 gg - 2, 1, -1]] <> ": " <> ToString[nonzeroCoeffs]];
zeroTerms = Select[entries, Together[#[[4]]] === 0 &];
emit["decorated graphs with identically vanishing term: " <> ToString[Length[zeroTerms]] <> " " <> ToString[zeroTerms[[All, 1]]]];
regIn = Select[entries, MemberQ[setG, #[[1, 1]]] && regQ[#[[5]]] &];
emit["regular decorated graphs inside the stable graphs with a vertex of genus >= 2: " <> ToString[Length[regIn]] <> ", of which identically zero: " <> ToString[Count[regIn, e_ /; Together[e[[4]]] === 0]] <> "; nonzero ones: " <> ToString[Select[regIn, Together[#[[4]]] =!= 0 &][[All, 1]]]];
(* decorated graphs of the vertex graph: leading coefficients versus the gap polynomial *)
v1 = Select[entries, #[[1, 1]] == 1 &];
ratios = Union[Table[Factor[e[[5, 1]]/gap[[1]]], {e, v1}]];
emit["vertex graph: " <> ToString[Length[v1]] <> " decorated graphs; leading coefficients proportional to the gap polynomial: " <> ToString[AllTrue[ratios, FreeQ[#, x] &]] <> " (" <> ToString[Length[ratios]] <> " distinct ratios); sums of their coefficients at t^-" <> ToString[Range[2 gg - 2, 1, -1]] <> ": " <> ToString[InputForm[Factor /@ Together /@ Total[v1[[All, 5]]]]]];
emit["vertex graph, per decorated graph {key {psi on half-edges, dilaton leaves}, polar part}: " <> ToString[InputForm[Table[{e[[2]], e[[5]]}, {e, v1}]]]];
(* ---- structure of the terms: denominators ---- *)
disc = Q (1 - x)^2 + 4 x;   (* discriminant of U^2 - S U + P at alpha0 = 0, alpha_i = eps_i (times x^2 Q); zero at Q_c = -4x/(1-x)^2 *)
(* the multiplicity of an irreducible factor f in a polynomial d; the residual is what remains after dividing out all the listed factors
   (it must be a nonzero constant for the claimed denominator structure to hold -- checked explicitly, so that no other factor, e.g. the
   discriminant written in another normalisation or an unexpected polynomial in Q and x, can hide in the denominators) *)
(* multiplicities from FactorList (a polynomial remainder with x-dependent coefficients does not reliably simplify to 0, which hid the
   discriminant factor in the first version of this script); a factor is identified with a listed one when their ratio is a number *)
knownFactors = {1 - Q, Q, disc, x, 1 + x, 1 - x};
denData[t_] := Module[{fl = FactorList[Denominator[Together[t]]], pw = ConstantArray[0, Length[knownFactors]], res = {}, hit},
   Do[If[NumberQ[fl[[i, 1]]], Continue[]];
    hit = Position[knownFactors, f_ /; NumberQ[Together[fl[[i, 1]]/f]], {1}, 1, Heads -> False];
    If[hit === {}, AppendTo[res, fl[[i]]], pw[[hit[[1, 1]]]] += fl[[i, 2]]], {i, Length[fl]}];
   Append[pw, If[res === {}, 1, Times @@ (Power @@@ res)]]];
dst = denData /@ stable[[All, "term"]];
emit["denominators of the stable-graph terms {power of (1-Q), power of Q, power of the discriminant Q(1-x)^2+4x, power of x, power of (1+x), power of (1-x), residual factor (1 = none)}:"];
emit["  tally over the " <> ToString[ngr] <> " stable graphs: " <> ToString[InputForm[Sort[Tally[dst]]]]];
emit["  max powers: " <> ToString[Max /@ Transpose[dst[[All, ;; 6]]]] <> "; graphs with a discriminant pole: " <> ToString[Select[Range[ngr], dst[[#, 3]] > 0 &]] <> "; graphs with a residual (unlisted) denominator factor: " <> ToString[Select[Range[ngr], dst[[#, 7]] =!= 1 &]]];
den = denData /@ entries[[All, 4]];
emit["denominators of the decorated-graph terms: max powers " <> ToString[Max /@ Transpose[den[[All, ;; 6]]]] <> "; with a discriminant pole: " <> ToString[Count[den, d_ /; d[[3]] > 0]] <> "; with a residual (unlisted) denominator factor: " <> ToString[Count[den, d_ /; d[[7]] =!= 1]] <> If[Count[den, d_ /; d[[7]] =!= 1] > 0, " " <> ToString[InputForm[Union[Select[den, #[[7]] =!= 1 &][[All, 7]]]]], ""] <> "; (1-Q)-power tally: " <> ToString[Sort[Tally[den[[All, 1]]]]] <> "; x-power tally: " <> ToString[Sort[Tally[den[[All, 4]]]]] <> "; (1+x)-power tally: " <> ToString[Sort[Tally[den[[All, 5]]]]]];
emit["VALIDATION " <> If[dst[[All, 7]] === ConstantArray[1, ngr] && Count[den, d_ /; d[[7]] =!= 1] == 0, "OK", "FAIL"] <> ": every denominator is a product of powers of (1-Q), Q, x, (1+x), (1-x) and the discriminant only"];
(* the discriminant poles (collision of the two idempotents at Q_c = -4x/(1-x)^2) of the individual terms must cancel in the total *)
emit["discriminant-power tally over the stable graphs {power, count}: " <> ToString[Sort[Tally[dst[[All, 3]]]]] <> "; over the decorated graphs: " <> ToString[Sort[Tally[den[[All, 3]]]]] <> "; total has a discriminant pole: " <> ToString[denData[tot0][[3]] > 0] <> " (its denominator: " <> ToString[InputForm[Factor[Denominator[tot0]]]] <> ")"];
emit["discriminant power versus the number of edges of the stable graph {edges, max power}: " <> ToString[Table[{ne, Max[Pick[dst[[All, 3]], Table[Length[stable[[k]]["edges"]], {k, ngr}], ne]]}, {ne, Union[Table[Length[stable[[k]]["edges"]], {k, ngr}]]}]]];
emit["degree in Q of numerator/denominator of the stable-graph terms (max): " <> ToString[Max[Exponent[Numerator[Together[#]], Q] & /@ stable[[All, "term"]]]] <> "/" <> ToString[Max[Exponent[Denominator[Together[#]], Q] & /@ stable[[All, "term"]]]]];
(* ---- regular stable graphs: pole at Q = 1 of the term itself (regular means the Laurent coefficients t^-k vanish; the term may still have (1-Q) in its denominator) ---- *)
emit["regular stable graphs whose term still has (1-Q) in the denominator (cancellation of the pole only in the Laurent expansion is impossible for a rational function, so these must be 0): " <> ToString[Select[regular, dst[[#, 1]] > 0 &]]];
emit["(1-Q)-power of the stable-graph terms versus the leading pole order of the polar part {k, (1-Q)-power, pole order}: " <> ToString[Table[{k, dst[[k, 1]], leadOrder[stable[[k]]["polar"]]}, {k, singular}]]];
(* ---- Q = 0 (constant maps): individual terms have a pole at Q = 0 (the canonical coordinates degenerate at Q = 0, P = -1/(Qx));
   the singular coefficients cancel in the sum and the total at Q = 0 is the constant-map term F_g(0) ---- *)
tot0 = Together[Total[stable[[All, "term"]]]];
resQ0 = Table[Factor[Residue[stable[[k]]["term"], {Q, 0}]], {k, ngr}];
emit["stable graphs with a pole at Q = 0: " <> ToString[Count[dst, d_ /; d[[2]] > 0]] <> " of " <> ToString[ngr] <> " (maximal order " <> ToString[Max[dst[[All, 2]]]] <> "); sum of the residues (coefficient of Q^-1): " <> ToString[Factor[Together[Total[resQ0]]]] <> "; residue of the vertex graph: " <> ToString[InputForm[resQ0[[1]]]]];
emit["the total is regular at Q = 0: " <> ToString[FreeQ[Denominator[tot0] /. Q -> 0, 0] && (Denominator[tot0] /. Q -> 0) =!= 0]];
emit["total at Q = 0 (constant-map contribution F_g(0)): " <> ToString[InputForm[Factor[tot0 /. Q -> 0]]] <> "; reference F_g(0): " <> ToString[InputForm[Factor[FrefHVC[gg, x] /. Q -> 0]]] <> "; agree: " <> ToString[Together[(tot0 /. Q -> 0) - (FrefHVC[gg, x] /. Q -> 0)] === 0]];
Export[FileNameJoin[{resDir, "hvc_analysis_g" <> ToString[gg] <> "_" <> tag <> ".txt"}], StringRiffle[out, "\n"], "Text"];
Print["written ", FileNameJoin[{resDir, "hvc_analysis_g" <> ToString[gg] <> "_" <> tag <> ".txt"}]];
