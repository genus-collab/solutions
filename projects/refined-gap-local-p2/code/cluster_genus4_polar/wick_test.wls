(* wick_test.wls -- regression and timing of the direct Wick expansion: equality with the x-series method at z^3 and z^6,
   timing of wickDirect[9], of the W_k derivatives at numeric weights to k = 20, and of the z^9 R-matrix at numeric weights. *)
SetDirectory[DirectoryName[$InputFileName]];
Get["rc_lib.m"];
oldWick[zmax_] := Module[{ser, tot},
   ser = Series[Exp[Sum[Wk[k] x^k/(k! z), {k, 3, 2 zmax + 2}]], {x, 0, 6 zmax}];
   tot = Sum[SeriesCoefficient[ser, 2 n] (2 n)! (-z)^n/(n! 2^n Wk[2]^n), {n, 0, 3 zmax}];
   Table[Expand[Coefficient[Normal[Series[tot, {z, 0, zmax}]], z, m]], {m, 0, zmax}]];
Do[Module[{t0 = AbsoluteTime[], a, b},
   a = wickDirect[zm]; b = oldWick[zm];
   Print["ZMAX = ", zm, ": direct == series method: ", And @@ Table[Expand[a[[m + 1]] - b[[m + 1]]] === 0, {m, 0, zm}], "; terms per order ", Length /@ (List @@@ a), " (", Round[AbsoluteTime[] - t0, 0.1], " s)"]], {zm, {3, 6}}];
t0 = AbsoluteTime[]; w9 = wickDirect[9];
Print["wickDirect[9]: ", Round[AbsoluteTime[] - t0, 0.1], " s; terms per order ", Length /@ (List @@@ w9)];
t0 = AbsoluteTime[]; setupWick[9]; Print["setupWick[9]: ", Round[AbsoluteTime[] - t0, 0.1], " s"];
ws = {a0 -> 5, a1 -> 3, a2 -> -7, e1 -> 2, e2 -> 11};
t0 = AbsoluteTime[]; computeWU[ws]; Print["W_k to k = 20 at numeric weights: ", Round[AbsoluteTime[] - t0, 0.1], " s; LeafCount W_20 = ", LeafCount[WU[20]]];
t0 = AbsoluteTime[]; rmatrixAtClosed[ws]; Print["rmatrixAtClosed at z^9: ", Round[AbsoluteTime[] - t0, 0.1], " s; LeafCount ratCN[1,2][9] = ", LeafCount[ratCN[1, 2][9]]];
t0 = AbsoluteTime[]; evalQ[3/7]; Print["evalQ at a numeric Q: ", Round[AbsoluteTime[] - t0, 0.1], " s"];
(* genus-3 graph sums at this Q as a sanity check of the z^9 factors (must equal the genus-3 result) and timing *)
prep3 = prepareGraph /@ graphsWithLegs[3, 0];
t0 = AbsoluteTime[]; v3 = graphSumC[prep3]; Print["genus-3 graph sum at z^9 factors: ", Round[AbsoluteTime[] - t0, 0.1], " s; symmetric ", cSymQ[v3]];
F30 = Get["ref/resolvedconifold_results.m"][[3]];
Print["genus-3 graph sum minus stored F30 at this point: ", Together[v3[[1]] - (F30 /. ws /. Q -> 3/7)]];
(* genus-4 graph sum at one Q: the per-Q cost of S5 *)
prep4 = prepareGraph /@ graphsWithLegs[4, 0];
Print["genus 4: ", Length[prep4], " graphs, ", Total[Length[#["decs"]] & /@ prep4], " labelled decorations"];
t0 = AbsoluteTime[]; v4 = graphSumC[prep4]; Print["genus-4 graph sum at one Q: ", Round[AbsoluteTime[] - t0, 0.1], " s; symmetric ", cSymQ[v4]];
