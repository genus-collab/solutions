(* review_verify_proof_B.wls -- genus-4 (and further genus-3) checks of the string-equation proof with SYMBOLIC weights u[1..3] and symbolic cv,
   the A/B string decomposition at genus 4, and a four-weight Hodge class at genus 2 (flat-unit generality). Time-guarded. Read-only on project files. *)
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
phi[cv_] := Function[k, cv^(-k - 1)]; one = Function[k, If[k == 0, 1, 0]]; zphi[cv_] := Function[k, If[k >= 1, cv^(-k), 0]];
Z[h_, n_, cv_, clf_] := Z[h, n, cv, clf] = ZW[h, ConstantArray[phi[cv], n], cv, clf];
clGen = Function[{h, a}, hodgeNormU4[h, a]];
clWK = Function[{h, a}, HI[h, 0, 0, a]];
hodgeNormN[g_, nw_, a_List] := hodgeNormN[g, nw, a] = Expand[Total[Table[cl[[2]] HI[g, cl[[1, 1]], cl[[1, 2]], a], {cl, hodgeClassN[g, nw]}]]];
cl4 = Function[{h, a}, hodgeNormN[h, 4, a]];                 (* four symbolic weights u[1..4] *)
budget = 150; guard[label_, expr_] := If[el[] > budget, Print[label, ": skipped (time budget)"], Print[label, " = ", TimeConstrained[expr, budget - el[], "timed out"], "  t = ", el[], " s"]];
Print["=== four symbolic weights (rank-one Hodge CohFT of C^4-type weights), genus 2 ==="];
guard["4 weights: Z[2,1], Z[2,2]", {Z[2, 1, cv, cl4], Z[2, 2, cv, cl4]}];
guard["4 weights: cv Z[1,1] - Int_{M_{1,1}} Omega", Together[cv Z[1, 1, cv, cl4] - hodgeNormN[1, 4, {0}]]];
Print["=== genus 3, symbolic weights: Z[3,3] and the A/B decomposition at (3,1) ==="];
guard["Z[3,3] generic", Z[3, 3, cv, clGen]];
guard["(3,1): A + B (= cv Z[3,1]) and A, B separately nonzero?", Module[{Ap = ZW[3, {one}, cv, clGen], Bp = ZW[3, {zphi[cv]}, cv, clGen]}, {Together[Ap + Bp], Ap =!= 0, Bp =!= 0}]];
Print["=== genus 4 (the prover's unrun check): symbolic weights, symbolic cv ==="];
guard["WK: Z^WK[4,1]", Z[4, 1, cv, clWK]];
guard["generic: Z[4,1]", Z[4, 1, cv, clGen]];
guard["generic: (4,1) A + B and A nonzero?", Module[{Ap = ZW[4, {one}, cv, clGen], Bp = ZW[4, {zphi[cv]}, cv, clGen]}, {Together[Ap + Bp], Ap =!= 0}]];
guard["generic: Z[4,2]", Z[4, 2, cv, clGen]];
guard["generic: Z[4,0] (gap, nonzero expected)", Factor[Z[4, 0, cv, clGen]]];
Print["total time ", el[], " s"];
