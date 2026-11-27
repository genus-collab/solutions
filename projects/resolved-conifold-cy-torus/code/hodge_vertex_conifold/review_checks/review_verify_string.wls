(* review_verify_string.wls -- independent verifier checks of the string-equation proof of Z[h,n] = 0 (h >= 2, n >= 1).
   Own enumeration (not the prover's ZG). Read-only on the project files. *)
SetDirectory[DirectoryName[$InputFileName]];
tt = AbsoluteTime[];
Get["hvc_lib.m"];
Print["library loaded in ", Round[AbsoluteTime[] - tt], " s"];

(* correlator <prod psi^a_i>_h with an arbitrary class function cl[h, a] *)
(* --- my own enumeration: all ORDERED tuples (small cases) --- *)
phiC[k_, cc_] := 1/cc^(k + 1);
TC[d_, cc_] := If[d >= 2, -1/cc^(d - 1), 0];
(* ordered tuples of n legs (>=0) and m dilaton leaves (>=2) with total psi-degree dim *)
Ztuples[h_, n_, cc_, cl_] := Module[{tot = 0, dim, legs, dils},
   Do[dim = 3 h - 3 + n + m;
    If[dim < 2 m || 2 h - 2 + n + m <= 0, Continue[]];
    legs = Select[Tuples[Range[0, dim], n], Total[#] <= dim &];
    Do[dils = Select[Tuples[Range[2, dim - Total[A]], m], Total[#] <= dim - Total[A] &];
     Do[tot += (Times @@ (phiC[#, cc] & /@ A)) (Times @@ (TC[#, cc] & /@ D))/m! cl[h, Join[A, D]], {D, dils}], {A, legs}], {m, 0, 3 h - 3 + n}];
   Together[tot]];
(* --- multiset enumeration for larger cases --- *)
sortedLists[k_, lo_, smax_] := Module[{rec}, rec[kk_, l_, s_] := If[kk == 0, {{}}, Join @@ Table[Prepend[#, a] & /@ rec[kk - 1, a, s - a], {a, l, s}]]; rec[k, lo, smax]];
nOrd[lst_] := Multinomial @@ (Last /@ Tally[lst]);
Zms[h_, n_, cc_, cl_] := Module[{tot = 0, dim},
   Do[dim = 3 h - 3 + n + m;
    If[dim < 2 m || 2 h - 2 + n + m <= 0, Continue[]];
    Do[Do[If[Total[A] + Total[D] <= dim,
        tot += nOrd[A] nOrd[D]/m! (Times @@ (phiC[#, cc] & /@ A)) (Times @@ (TC[#, cc] & /@ D)) cl[h, Join[A, D]]],
       {D, sortedLists[m, 2, dim - Total[A]]}], {A, sortedLists[n, 0, dim]}], {m, 0, 3 h - 3 + n}];
   Together[tot]];

clGen = Function[{h, a}, hodgeNormU4[h, a]];                       (* generic weights, u[a] = 1/w_a symbolic *)
clCY = Function[{h, a}, hodgeNormU4[h, a] /. usubOf[{1 + x, -x, -1}]];
clWK = Function[{h, a}, hodgeNormU4[h, a] /. u[_] -> 0];

Print["=== A. conventions: CY weights {1+x,-x,-1}, c = 1+x, against the task's values (ordered-tuple enumeration) ==="];
Print["Z[1,1] - (1+x+x^2)/(24x(1+x)^2) = ", Together[Ztuples[1, 1, 1 + x, clCY] - (1 + x + x^2)/(24 x (1 + x)^2)]];
Print["Z[1,2] - (1+x+x^2)/(24x(1+x)^4) = ", Together[Ztuples[1, 2, 1 + x, clCY] - (1 + x + x^2)/(24 x (1 + x)^4)]];
Print["Z[1,3] - (1+x+x^2)/(12x(1+x)^6) = ", Together[Ztuples[1, 3, 1 + x, clCY] - (1 + x + x^2)/(12 x (1 + x)^6)]];
Print["Z[2,0] + (7+10x^2+7x^4)/(5760x^2(1+x)^3) = ", Together[Ztuples[2, 0, 1 + x, clCY] + (7 + 10 x^2 + 7 x^4)/(5760 x^2 (1 + x)^3)]];
Print["Z[2,1], Z[2,2] = ", {Ztuples[2, 1, 1 + x, clCY], Ztuples[2, 2, 1 + x, clCY]}];
Print["multiset vs tuple enumeration agree (2,2),(1,3): ", {Together[Zms[2, 2, 1 + x, clCY] - Ztuples[2, 2, 1 + x, clCY]], Together[Zms[1, 3, 1 + x, clCY] - Ztuples[1, 3, 1 + x, clCY]]}];
Print["time ", Round[AbsoluteTime[] - tt], " s"];

Print["=== B. Lemma 1 ingredients from the library: <tau_0>_2 = ", clGen[2, {0}], ";  <tau_0>_1 = ", Factor[clGen[1, {0}]], " (expect -(u1+u2+u3)/24)"];
Print["string equation on Hodge correlators (generic u), genus 2, a = {1,2},{3},{2,2,1}: ",
  Table[Together[clGen[2, Append[a, 0]] - Sum[If[a[[i]] >= 1, clGen[2, ReplacePart[a, i -> a[[i]] - 1]], 0], {i, Length[a]}]], {a, {{1, 2}, {3}, {2, 2, 1}, {4, 1}}}]];

Print["=== C. the key step (S') with power-series insertions t_s = (c s - z^2)/(c - z), genus 2 and 1 (generic u, symbolic c, s) ==="];
(* <A_1(psi_1) ... A_N(psi_N)>_h for power series given as coefficient functions, truncated at the dimension *)
corrSeries[h_, fs_List, cl_] := Module[{N = Length[fs], dim = 3 h - 3 + Length[fs], tot = 0},
   If[2 h - 2 + N <= 0, Return[0]];
   Do[If[Total[a] <= dim, tot += (Times @@ MapThread[#1[#2] &, {fs, a}]) cl[h, a]], {a, Tuples[Range[0, dim], N]}];
   Together[tot]];
one[k_] := If[k == 0, 1, 0];
tsC[k_] := If[k == 0, s, If[k == 1, s/c, s/c^k - 1/c^(k - 1)]];           (* coefficients of (c s - z^2)/(c - z) = s + s z/c + Sum_{k>=2} (s/c^k - c^(1-k)) z^k *)
Print["check t_s coefficients: ", Simplify[Sum[tsC[k] z^k, {k, 0, 12}] - Normal[Series[(c s - z^2)/(c - z), {z, 0, 12}]]]];
lowC[k_] := If[k == 0, 1 + (s - c)/c, (s - c)/c^(k + 1)];                (* coefficients of (t_s(z) - t_s(0))/z = 1 + (s - c) phi(z) *)
Print["check Lemma 2 coefficients: ", Simplify[Sum[lowC[k] z^k, {k, 0, 12}] - Normal[Series[((c s - z^2)/(c - z) - s)/z, {z, 0, 12}]]]];
Do[Print["genus ", h, ", M = ", M, ": <1 t_s^M>_{h,M+1} - M <[(t_s-s)/psi] t_s^(M-1)>_{h,M} = ",
   Together[corrSeries[h, Prepend[ConstantArray[tsC, M], one], clGen] - M corrSeries[h, Prepend[ConstantArray[tsC, M - 1], lowC], clGen]]], {h, 1, 2}, {M, 1, 3}];
Print["genus 1, M = 0: <1>_{1,1} = ", clGen[1, {0}], " (the anomaly)"];
Print["time ", Round[AbsoluteTime[] - tt], " s"];

Print["=== D. predicted universality: generic non-CY weights u1,u2,u3 and unrelated symbolic c ==="];
Print["Z[2,1] = ", Zms[2, 1, c, clGen], ";  Z[2,2] = ", Zms[2, 2, c, clGen], ";  Z[2,3] = ", Zms[2, 3, c, clGen]];
Do[Print["Z[1,", n, "] - (n-1)! c^(1-2n) (-(u1+u2+u3)/24) = ", Together[Zms[1, n, c, clGen] - (n - 1)! c^(1 - 2 n) (-(u[1] + u[2] + u[3])/24)]], {n, 1, 4}];
Print["Z[2,0] generic = ", Factor[Zms[2, 0, c, clGen]]];
Print["time ", Round[AbsoluteTime[] - tt], " s"];

Print["=== E. Witten-Kontsevich vertex (u = 0), symbolic c ==="];
Print["Z^WK[2,1], Z^WK[2,2], Z^WK[3,1], Z^WK[3,2] = ", {Zms[2, 1, c, clWK], Zms[2, 2, c, clWK], Zms[3, 1, c, clWK], Zms[3, 2, c, clWK]}];
Print["Z^WK[1,1], Z^WK[1,2] = ", {Zms[1, 1, c, clWK], Zms[1, 2, c, clWK]}, ";  Z^WK[2,0], Z^WK[3,0] = ", {Zms[2, 0, c, clWK], Zms[3, 0, c, clWK]}];
Print["time ", Round[AbsoluteTime[] - tt], " s"];

Print["=== F. genus 3, Hodge class with random non-CY rational weights and unrelated rational c ==="];
SeedRandom[7];
wR = {3/2, -5/7, 11/3}; cR = 13/4;
clR = Function[{h, a}, hodgeNormU4[h, a] /. usubOf[wR]];
Print["weights ", wR, " (sum ", Total[wR], "), c = ", cR];
Print["Z[3,1] = ", Zms[3, 1, cR, clR], "   (", Round[AbsoluteTime[] - tt], " s)"];
Print["Z[3,2] = ", Zms[3, 2, cR, clR], "   (", Round[AbsoluteTime[] - tt], " s)"];
Print["Z[3,0] = ", Zms[3, 0, cR, clR]];
Print["total time ", Round[AbsoluteTime[] - tt], " s"];
