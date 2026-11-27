(* rc_collect.m -- collector for the general-torus conifold sampling: assembles F11, F20, ..., F<genus>0 exactly and runs the checks.
   Usage:  wolframscript -file rc_collect.m genus ZMAX
   Input:  rc_samples_g<genus>_w*.m (all workers), rc_const_g<genus>.m (degree-zero term; computed here if absent),
           ref/resolvedconifold_results.m (genus 1-3 reference).
   Output: resolvedconifold_F<genus>0.m = {F11, F20, ..., F<genus>0} (exact rational functions of Q and the five weights),
           rc_collect.log with the checks and VALIDATION lines. *)
SetDirectory[DirectoryName[$InputFileName]];
{gg, zmax} = ToExpression /@ $ScriptCommandLine[[2 ;; 3]];
Get["rc_lib.m"];
log[m_] := Module[{s = OpenAppend["rc_collect.log"]}, WriteString[s, DateString[{"Year", "-", "Month", "-", "Day", " ", "Hour", ":", "Minute", ":", "Second"}], "  [g=", gg, "] ", m, "\n"]; Close[s]; Print[m]];
files = FileNames["rc_samples_g" <> ToString[gg] <> "_w*.m"];
If[files === {}, log["no sample files found"]; Exit[1]];
allSamples = <||>;
Do[Do[allSamples[r[[1]]] = {r[[2]], r[[3]]}, {r, ReadList[f]}], {f, files}];
log[ToString[Length[files]] <> " files, " <> ToString[Length[allSamples]] <> " weight samples (needed: " <> ToString[Length[orbitBasis[2 gg]] + 3] <> " for the degree-" <> ToString[2 gg] <> " basis of " <> ToString[Length[orbitBasis[2 gg]]] <> " symmetric monomials)"];
res = <||>; okall = True;
Do[Module[{r = assemble[gi, 2 gi]}, If[r === $Failed, okall = False; log["assembly FAILED for genus " <> ToString[gi]], res[gi] = r; log["genus " <> ToString[gi] <> " instanton part assembled: pole order " <> ToString[r[[2]]] <> ", LeafCount " <> ToString[LeafCount[r[[1]]]]]]], {gi, 1, gg}];
If[! okall, Exit[1]];
(* degree-zero terms *)
setupWick[zmax];
F0[1] = F11const;
Do[F0[gi] = If[gi == gg && FileExistsQ["rc_const_g" <> ToString[gg] <> ".m"], Get["rc_const_g" <> ToString[gg] <> ".m"],
     If[gi <= 3 && FileExistsQ["ref/resolvedconifold_results.m"], Get["ref/resolvedconifold_results.m"][[gi]] /. Q -> 0, const0[gi]]], {gi, 2, gg}];
Do[Fall[gi] = res[gi][[1]] + F0[gi], {gi, 1, gg}];
Put[Table[Fall[gi], {gi, 1, gg}], "resolvedconifold_F" <> ToString[gg] <> "0.m"];
(* validation against the stored genus-1..3 results *)
If[FileExistsQ["ref/resolvedconifold_results.m"],
  Module[{ref = Get["ref/resolvedconifold_results.m"], ok},
   Do[ok = Together[Fall[gi] - ref[[gi]]] === 0;
      log[If[ok, "VALIDATION OK: genus " <> ToString[gi] <> " agrees with the stored result", "VALIDATION FAIL: genus " <> ToString[gi] <> " differs from the stored result"]], {gi, 1, Min[gg, 3]}]]];
(* checks at the top genus *)
CYr = {a0 -> e1 + e2 - a1 - a2};
Finst = Fall[gg] - F0[gg];
cyok = Together[(Finst /. CYr) - Pg[gg] PolyLog[3 - 2 gg, Q]/(e1 e2)] === 0;
log["CY torus: F" <> ToString[gg] <> "0 - F" <> ToString[gg] <> "0(0) = P_" <> ToString[gg] <> " Li_{" <> ToString[3 - 2 gg] <> "}(Q)/(eps1 eps2), independent of alpha1, alpha2: " <> ToString[cyok]];
log[If[cyok, "VALIDATION OK: CY restriction (rigidity) at genus " <> ToString[gg], "VALIDATION FAIL: CY restriction at genus " <> ToString[gg]]];
ep = (e1 + e2)/2; em = (e1 - e2)/2; c1 = e1 + e2 - a0 - a1 - a2;
deg1 = Factor[SeriesCoefficient[Series[Fall[gg] e1 e2, {Q, 0, 1}], 1]];
log["degree-1 term eps1 eps2 [Q^1] F_" <> ToString[gg] <> " = " <> ToString[InputForm[deg1]]];
log["   its CY part: " <> ToString[InputForm[Factor[deg1 /. CYr]]]];
nc = Together[(deg1 - (deg1 /. CYr))/(c1 ep)];
log["   (non-CY part)/(c1 e+) = " <> ToString[InputForm[Factor[nc]]] <> "; polynomial in the weights: " <> ToString[PolynomialQ[Numerator[nc], {a0, a1, a2, e1, e2}] && NumberQ[Denominator[nc]]]];
quo = Together[(Finst - (Finst /. CYr))/(c1 ep)];
log["(F - F|_CY)/(c1 e+) has no weight denominators besides eps1 eps2: " <> ToString[FreeQ[Factor[Denominator[quo]], a0 | a1 | a2]]];
log["unrefined limit (eps2 -> -eps1) of the instanton part: " <> ToString[InputForm[Factor[Finst /. e2 -> -e1]]]];
log["NS limit (eps2 -> 0) of eps1 eps2 x instanton part: " <> ToString[InputForm[Factor[Limit[e1 e2 Finst, e2 -> 0]]]]];
log["F" <> ToString[gg] <> "0 instanton part = " <> ToString[InputForm[Collect[Together[Finst], Q, Factor]]]];
log["F" <> ToString[gg] <> "0 constant term  = " <> ToString[InputForm[F0[gg]]]];
log["done"];
