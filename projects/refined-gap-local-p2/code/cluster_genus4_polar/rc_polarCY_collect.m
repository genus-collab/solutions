(* rc_polarCY_collect.m -- collector for the CY-torus per-decorated-graph polar parts (alpha1, alpha2 free).
   Usage:  wolframscript -file rc_polarCY_collect.m genus
   Input:  rc_polarCY_fit_g<genus>_c*.m (all chunks), ref/conifold_CY_graphterms_polar.m (genus 2, 3 validation).
   Output: conifold_CY_graphterms_polar_genus<genus>.m, rc_polarCY_collect.log with structure checks and VALIDATION lines. *)
SetDirectory[DirectoryName[$InputFileName]];
gg = ToExpression[$ScriptCommandLine[[2]]];
Get["rc_lib.m"];
log[m_] := Module[{s = OpenAppend["rc_polarCY_collect.log"]}, WriteString[s, DateString[{"Year", "-", "Month", "-", "Day", " ", "Hour", ":", "Minute", ":", "Second"}], "  [g=", gg, "] ", m, "\n"]; Close[s]; Print[m]];
files = FileNames["rc_polarCY_fit_g" <> ToString[gg] <> "_c*.m"];
If[files === {}, log["no fit files found"]; Exit[1]];
data = Get /@ files;
entries = SortBy[Join @@ data[[All, 1]], First]; failed = Join @@ data[[All, 2]];
log[ToString[Length[files]] <> " chunk files, " <> ToString[Length[entries]] <> " decorated graphs; failed fits: " <> ToString[failed]];
ok = Select[entries, FreeQ[#[[4]], $Failed] &];
(* (u, alpha0) form: symmetric in alpha1 <-> alpha2, so a function of s1 = a1 + a2 = e1 + e2 - a0 and p1 = a1 a2 = (u + (e1+e2)^2 - 2 (e1+e2) a0)/4 *)
uform[f_] := Module[{g = Together[f], num, den, red},
   num = Numerator[g]; den = Denominator[g];
   If[! FreeQ[den, a1 | a2], Return[f]];
   red = SymmetricReduction[num, {a1, a2}, {s1, p1}];
   If[red[[2]] =!= 0, Return[f]];
   Factor[(red[[1]] /. {s1 -> e1 + e2 - \[Alpha][0], p1 -> (u + (e1 + e2)^2 - 2 (e1 + e2) \[Alpha][0])/4})/den]];
(* structure checks *)
Do[Module[{cs = Table[e[[4, n]], {e, ok}], nind, ndivu},
   nind = Count[cs, c_ /; FreeQ[c, a1 | a2]];
   ndivu = Count[cs, c_ /; Together[c /. a1 -> (e1 + e2)/2] === 0 && Together[c /. a2 -> (e1 + e2)/2] === 0];
   log["order t^-" <> ToString[2 gg - 1 - n] <> ": independent of alpha1, alpha2 for " <> ToString[nind] <> " of " <> ToString[Length[cs]] <> " graphs; divisible by u = (eps1+eps2-2 alpha1)(eps1+eps2-2 alpha2) for " <> ToString[ndivu] <> "; identically zero for " <> ToString[Count[cs, 0]]]], {n, 1, 2 gg - 2}];
tot = Together /@ Total[Table[e[[4]], {e, ok}]];
expected = cyPolar[gg];
totok = Length[ok] == Length[entries] && (Together /@ (tot - expected)) === ConstantArray[0, 2 gg - 2];
log["sum over all decorated graphs minus the polar part of P_g Li_{3-2g}(Q)/(eps1 eps2): " <> ToString[InputForm[Factor /@ (tot - expected)]]];
log[If[totok, "VALIDATION OK: CY-torus totals", "VALIDATION FAIL: CY-torus totals (or failed fits)"]];
(* output *)
readme = "Resolved conifold x C2 on the CY torus alpha0 = eps1 + eps2 - alpha1 - alpha2 (alpha1, alpha2 free), genus " <> ToString[gg] <> ": polar parts at Q = 1 in t = Log Q of the " <> ToString[Length[entries]] <> " decorated stable graphs " <>
  "(same enumeration, keys and descriptions as the other per-graph files). Entry: {id, key, description, coefficients of t^-(2g-2), ..., t^-1 as rational functions of alpha1, alpha2, eps1, eps2 (homogeneous of degree 2g-2), " <>
  "the same written in u = (eps1+eps2-2 alpha1)(eps1+eps2-2 alpha2) and alpha0 = eps1+eps2-alpha1-alpha2 (symbols u and Alpha[0]), fit ansatz {numerator degree, eps-exponent, (eps1+eps2)-exponent, alpha-degree} per order}. " <>
  "Weight 1/|Aut| times the sum over labelled decorations and colourings; reconstructed from exact integer weight samples and verified on spare points.";
Put[Prepend[Table[{e[[1]], e[[2]], e[[3]], e[[4]], uform /@ e[[4]], e[[5]]}, {e, entries}], readme], "conifold_CY_graphterms_polar_genus" <> ToString[gg] <> ".m"];
log["wrote conifold_CY_graphterms_polar_genus" <> ToString[gg] <> ".m"];
(* validation against the stored genus-2/3 results *)
If[gg <= 3 && FileExistsQ["ref/conifold_CY_graphterms_polar.m"],
  Module[{ref = Get["ref/conifold_CY_graphterms_polar.m"], refE, refA, np = 0},
   refE = gg /. ref[[2]];
   refA = Association[Table[r[[1]] -> r, {r, refE}]];
   Do[If[KeyExistsQ[refA, e[[1]]] && And @@ Table[Together[e[[4, n]] - refA[e[[1]]][[3, n]]] === 0, {n, 2 gg - 2}], np++], {e, ok}];
   log[If[np == Length[refE] && Length[entries] == Length[refE], "VALIDATION OK: all " <> ToString[Length[refE]] <> " genus-" <> ToString[gg] <> " CY-torus polar parts agree with the stored results",
     "VALIDATION FAIL: genus-" <> ToString[gg] <> " CY-torus polar parts: " <> ToString[np] <> " of " <> ToString[Length[refE]] <> " agree (computed " <> ToString[Length[entries]] <> ")"]]]];
log["done"];
