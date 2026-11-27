(* review_verify_proof_C.wls -- consistency of the genus-4 machinery used above: (i) gap formula (-e1 e2 c^3)^(h-1) Z[h,0] = P_h (2h-3)!/(e1 e2), h = 2,3,4,
   P_h = [t^(2h)] e1 e2 t^2/(4 sinh(e1 t/2) sinh(e2 t/2)), at e2 = 1, e1 = x, c = 1+x; (ii) Lemma 1 (string equation) at genus 4 with symbolic weights;
   (iii) Z[4,3] with symbolic weights and cv (time-guarded). Read-only on project files. *)
SetDirectory[DirectoryName[$InputFileName]];
tt = AbsoluteTime[]; el[] := Round[AbsoluteTime[] - tt];
Get["hvc_lib.m"];
msets[k_, lo_, smax_] := Module[{rec}, rec[kk_, l_, s_] := If[kk == 0, {{}}, Join @@ Table[Prepend[#, a] & /@ rec[kk - 1, a, s - a], {a, l, s}]]; rec[k, lo, smax]];
nOrd[lst_] := Length[lst]!/Times @@ (Last /@ Tally[lst]!);
ZW[h_, F_List, cv_, clf_] := Module[{n = Length[F], tot = 0, dim, legs, w},
   Do[dim = 3 h - 3 + n + m;
    If[dim < 2 m || 2 h - 2 + n + m <= 0, Continue[]];
    legs = Select[Tuples[Range[0, dim], n], Total[#] <= dim &];
    Do[w = Times @@ MapThread[#1[#2] &, {F, A}];
     If[w === 0, Continue[]];
     Do[tot += w (nOrd[D]/m!) (Times @@ (-cv^(1 - #) & /@ D)) clf[h, Join[A, D]], {D, msets[m, 2, dim - Total[A]]}],
     {A, legs}],
    {m, 0, 3 h - 3 + n}];
   Together[tot]];
phi[cv_] := Function[k, cv^(-k - 1)]; one = Function[k, If[k == 0, 1, 0]];
Z[h_, n_, cv_, clf_] := Z[h, n, cv, clf] = ZW[h, ConstantArray[phi[cv], n], cv, clf];
clGen = Function[{h, a}, hodgeNormU4[h, a]];
clCY = Function[{h, a}, hodgeNormU4[h, a] /. usubOf[{1 + x, -x, -1}]];
budget = 150; guard[label_, expr_] := If[el[] > budget, Print[label, ": skipped (time budget)"], Print[label, " = ", TimeConstrained[expr, budget - el[], "timed out"], "  t = ", el[], " s"]];
Ph[h_] := SeriesCoefficient[Series[x t^2/(4 Sinh[x t/2] Sinh[t/2]), {t, 0, 2 h}], 2 h];
Print["=== (i) gap formula at CY weights, e1 = x, e2 = 1, c = 1+x: (-x c^3)^(h-1) Z[h,0] - P_h (2h-3)!/x, h = 2, 3, 4 ==="];
Do[guard["h = " <> ToString[h], Together[(-x (1 + x)^3)^(h - 1) Z[h, 0, 1 + x, clCY] - Ph[h] (2 h - 3)!/x]], {h, 2, 4}];
Print["  (stated t^-6 coefficient at genus 4: (381+620x^2+686x^4+620x^6+381x^8)/(3870720 x); P_4 3!/x = ", Factor[Ph[4] 3!/x], ")"];
Print["=== (ii) Lemma 1 at genus 4, symbolic weights, random polynomial insertions ==="];
corr[h_, F_List, clf_] := Module[{n = Length[F], dim = 3 h - 3 + Length[F], tot = 0, w},
   Do[w = Times @@ MapThread[#1[#2] &, {F, A}]; If[w =!= 0, tot += w clf[h, A]], {A, Select[Tuples[Range[0, dim], n], Total[#] <= dim &]}]; Together[tot]];
poly[cl_List] := Function[k, If[k < Length[cl], cl[[k + 1]], 0]];
Dop[cl_List] := Function[k, If[k + 1 < Length[cl], cl[[k + 2]], 0]];
f1 = {3, 2, 5, 1, 7, 2, 1, 4, 6, 2, 1, 1, 3}; f2 = {1, -1, 0, 4, 2, 3, 5, 1, 1, 2, 2, 1, 1};
guard["h=4, N=1: <1,f1> - <Df1>", Together[corr[4, {one, poly[f1]}, clGen] - corr[4, {Dop[f1]}, clGen]]];
guard["h=4, N=2: <1,f1,f2> - <Df1,f2> - <f1,Df2>", Together[corr[4, {one, poly[f1], poly[f2]}, clGen] - corr[4, {Dop[f1], poly[f2]}, clGen] - corr[4, {poly[f1], Dop[f2]}, clGen]]];
Print["=== (iii) Z[4,3], symbolic weights and cv ==="];
guard["Z[4,3] generic", Z[4, 3, cv, clGen]];
Print["total time ", el[], " s"];
