(* review_verify_proof_A.wls -- independent verification of the string-equation proof of Z[h,n] = 0 (h >= 2, n >= 1).
   Z[h,n] is re-implemented from the task's definition (phi_k = c^{-k-1}, T_d = -c^{1-d} for d >= 2, 1/m!, legs distinguishable).
   Weights u[1], u[2], u[3] are kept SYMBOLIC (no CY condition) and the leg parameter cv is symbolic. Read-only on project files. *)
SetDirectory[DirectoryName[$InputFileName]];
tt = AbsoluteTime[]; el[] := Round[AbsoluteTime[] - tt];
Get["hvc_lib.m"];
Print["library loaded, t = ", el[], " s"];
msets[k_, lo_, smax_] := Module[{rec}, rec[kk_, l_, s_] := If[kk == 0, {{}}, Join @@ Table[Prepend[#, a] & /@ rec[kk - 1, a, s - a], {a, l, s}]]; rec[k, lo, smax]];
nOrd[lst_] := Length[lst]!/Times @@ (Last /@ Tally[lst]!);
(* ZW[h, F, cv, clf] = Sum_m 1/m! Int_{M_{h,n+m}} Omega prod_i F_i(psi_i) prod_j T(psi_{n+j});  F = list of n coefficient functions (distinguishable legs),
   clf[h, a] = Int_{M_{h,|a|}} Omega psi^a; dilaton leaves T = -psi^2/(cv - psi) as multisets with ordering factors *)
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
phi[cv_] := Function[k, cv^(-k - 1)];                 (* 1/(cv - psi) *)
one = Function[k, If[k == 0, 1, 0]];                  (* tau_0 *)
zphi[cv_] := Function[k, If[k >= 1, cv^(-k), 0]];     (* psi/(cv - psi) *)
Z[h_, n_, cv_, clf_] := Z[h, n, cv, clf] = ZW[h, ConstantArray[phi[cv], n], cv, clf];
clGen = Function[{h, a}, hodgeNormU4[h, a]];                                  (* Hodge class, symbolic u[1..3] *)
clCY = Function[{h, a}, hodgeNormU4[h, a] /. usubOf[{1 + x, -x, -1}]];       (* task's CY weights, e2 = 1, e1 = x *)
clWK = Function[{h, a}, HI[h, 0, 0, a]];                                     (* trivial CohFT *)
kappaGen = hodgeNormU4[1, {0}];
Print["Int_{M_{1,1}} Hodge_1 = ", kappaGen, "  (expected -(u1+u2+u3)/24: ", Together[kappaGen + (u[1] + u[2] + u[3])/24] === 0, ")"];

Print["=== A. conventions: CY weights {1+x,-x,-1}, c = 1+x, against the values stated in the task ==="];
cx = 1 + x;
Print["Z[2,0] + (7+10x^2+7x^4)/(5760 x^2 (1+x)^3) = ", Together[Z[2, 0, cx, clCY] + (7 + 10 x^2 + 7 x^4)/(5760 x^2 (1 + x)^3)]];
Print["Z[1,1..3] minus stated = ", {Together[Z[1, 1, cx, clCY] - (1 + x + x^2)/(24 x (1 + x)^2)], Together[Z[1, 2, cx, clCY] - (1 + x + x^2)/(24 x (1 + x)^4)], Together[Z[1, 3, cx, clCY] - (1 + x + x^2)/(12 x (1 + x)^6)]}];
Print["Z[2,1], Z[2,2], Z[3,1] (CY) = ", {Z[2, 1, cx, clCY], Z[2, 2, cx, clCY], Z[3, 1, cx, clCY]}, "  t = ", el[], " s"];

Print["=== B. Lemma 1 (string equation, operator form) for random polynomial insertions, symbolic weights ==="];
corr[h_, F_List, clf_] := Module[{n = Length[F], dim = 3 h - 3 + Length[F], tot = 0, w},
   Do[w = Times @@ MapThread[#1[#2] &, {F, A}]; If[w =!= 0, tot += w clf[h, A]], {A, Select[Tuples[Range[0, dim], n], Total[#] <= dim &]}]; Together[tot]];
poly[cl_List] := Function[k, If[k < Length[cl], cl[[k + 1]], 0]];
Dop[cl_List] := Function[k, If[k + 1 < Length[cl], cl[[k + 2]], 0]];
f1 = {3, 2, 5, 1, 7, 2, 1, 4, 6}; f2 = {1, -1, 0, 4, 2, 3, 5, 1, 1};
Print["h=2, N=2: <1,f1,f2> - <Df1,f2> - <f1,Df2> = ", Together[corr[2, {one, poly[f1], poly[f2]}, clGen] - corr[2, {Dop[f1], poly[f2]}, clGen] - corr[2, {poly[f1], Dop[f2]}, clGen]]];
Print["h=3, N=1: <1,f1> - <Df1> = ", Together[corr[3, {one, poly[f1]}, clGen] - corr[3, {Dop[f1]}, clGen]]];
Print["h=3, N=2: <1,f1,f2> - <Df1,f2> - <f1,Df2> = ", Together[corr[3, {one, poly[f1], poly[f2]}, clGen] - corr[3, {Dop[f1], poly[f2]}, clGen] - corr[3, {poly[f1], Dop[f2]}, clGen]]];
Print["h=2,3, N=0: <1>_{h,1} = Int Hodge_h on M_{h,1} = ", {hodgeNormU4[2, {0}], hodgeNormU4[3, {0}]}, ";  h=1: <1,f1>_{1,2} - <Df1>_{1,1} = ", Together[corr[1, {one, poly[f1]}, clGen] - corr[1, {Dop[f1]}, clGen]], "  t = ", el[], " s"];

Print["=== C. generic (non-CY) symbolic weights u[1..3], symbolic cv: the theorem's statements ==="];
Print["Z[2,1] = ", Z[2, 1, cv, clGen], ";  Z[2,2] = ", Z[2, 2, cv, clGen], ";  Z[2,3] = ", TimeConstrained[Z[2, 3, cv, clGen], 60, "timed out"], "  t = ", el[], " s"];
Print["Z[2,0] = ", Factor[Z[2, 0, cv, clGen]]];
Print["genus 1 recursion: cv Z[1,1] - kappa = ", Together[cv Z[1, 1, cv, clGen] - kappaGen], ";  cv Z[1,n] - (n-1)/cv Z[1,n-1], n = 2,3,4: ", Table[Together[cv Z[1, n, cv, clGen] - (n - 1)/cv Z[1, n - 1, cv, clGen]], {n, 2, 4}]];
Print["Z[1,n] - (n-1)! kappa cv^(1-2n), n = 1..4: ", Table[Together[Z[1, n, cv, clGen] - (n - 1)! kappaGen cv^(1 - 2 n)], {n, 1, 4}], "  t = ", el[], " s"];

Print["=== D. Lemma 3 pieces: A = Sum 1/m! <1, phi^(n-1), T^m>, B = Sum 1/m! <psi phi, phi^(n-1), T^m>; check A + B = cv Z[h,n] and A = (n-1)/cv Z[h,n-1] - B + delta kappa ==="];
Do[Module[{Ap, Bp}, Ap = ZW[h, Prepend[ConstantArray[phi[cv], n - 1], one], cv, clGen]; Bp = ZW[h, Prepend[ConstantArray[phi[cv], n - 1], zphi[cv]], cv, clGen];
   Print["(h,n) = ", {h, n}, ": A = ", Factor[Ap], ";  B = ", Factor[Bp], ";  A + B - cv Z = ", Together[Ap + Bp - cv Z[h, n, cv, clGen]], ";  A - (n-1)/cv Z[h,n-1] + B - delta kappa = ", Together[Ap - (n - 1)/cv Z[h, n - 1, cv, clGen] + Bp - If[h == 1 && n == 1, kappaGen, 0]]]],
  {h, 1, 2}, {n, 1, 2}];
Print["t = ", el[], " s"];

Print["=== E. trivial CohFT (Witten-Kontsevich), symbolic cv ==="];
Print["Z^WK[2,0] + 1/(1440 cv^3) = ", Together[Z[2, 0, cv, clWK] + 1/(1440 cv^3)], ";  Z^WK[2,1], Z^WK[2,2], Z^WK[3,1], Z^WK[3,2] = ", {Z[2, 1, cv, clWK], Z[2, 2, cv, clWK], Z[3, 1, cv, clWK], Z[3, 2, cv, clWK]}];
Print["Z^WK[1,n], n=1..3 (anomaly 0 expected) = ", Table[Z[1, n, cv, clWK], {n, 1, 3}]];
(* the m-decomposition of Z^WK[2,1] quoted in Remark 1: (5 - 78 + 273 - 340 + 140)/(5760 cv^5) *)
termM[m_] := Module[{dim = 4 + m, tot = 0}, Do[Do[tot += cv^(-k - 1) (nOrd[D]/m!) (Times @@ (-cv^(1 - #) & /@ D)) HI[2, 0, 0, Join[{k}, D]], {D, msets[m, 2, dim - k]}], {k, 0, dim}]; Together[tot]];
Print["m-terms of Z^WK[2,1] times 5760 cv^5: ", Table[Together[5760 cv^5 termM[m]], {m, 0, 4}], "  (expected {5, -78, 273, -340, 140})"];
Print["=== F. genus 3, generic symbolic weights and cv (time-guarded) ==="];
budget = 150; guard[label_, expr_] := If[el[] > budget, Print[label, ": skipped (time budget)"], Print[label, " = ", TimeConstrained[expr, budget - el[], "timed out"], "  t = ", el[], " s"]];
guard["Z[3,1] generic", Z[3, 1, cv, clGen]];
guard["Z[3,2] generic", Z[3, 2, cv, clGen]];
Print["total time ", el[], " s"];
