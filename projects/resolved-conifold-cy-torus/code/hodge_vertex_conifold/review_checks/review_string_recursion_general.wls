(* review_string_recursion_general.wls -- the vanishing Z[h,n] = 0 (h >= 2, n >= 1) and the genus-1 anomaly c Z[1,1] = Int_{M_{1,1}} Omega_{1,1}(1)
   use neither the CY condition, nor the value of c, nor the Hodge class: they follow from the string equation alone (flat unit).
   Tests: Witten-Kontsevich vertex (point class), three non-CY weights, the CY triple with an unrelated (numeric or symbolic) c,
   a rank-one Hodge class with four weights.  Read-only on the project files. *)
SetDirectory[DirectoryName[$InputFileName]];
tt = AbsoluteTime[];
Get["hvc_lib.m"];
msets[k_, lo_, smax_] := Module[{rec}, rec[kk_, l_, s_] := If[kk == 0, {{}}, Join @@ Table[Prepend[#, a] & /@ rec[kk - 1, a, s - a], {a, l, s}]]; rec[k, lo, smax]];
ordFactor[lst_] := Length[lst]!/Times @@ (Last /@ Tally[lst]!);
(* ZG[h, n, cv, clf] = Sum_m (1/m!) Int_{M_{h,n+m}} (class) prod_{i<=n} phi(psi_i) prod_{j<=m} T(psi_{n+j}), phi = 1/(cv - psi), T = -psi^2/(cv - psi),
   clf[h, a] = Int_{M_{h,Length[a]}} (class) psi^a; legs and dilaton leaves as multisets weighted by the number of orderings *)
ZG[h_, n_, cv_, clf_] := ZG[h, n, cv, clf] = Module[{tot = 0, dim},
   Do[dim = 3 h - 3 + n + m; If[dim < 2 m || 2 h - 2 + n + m <= 0, Continue[]];
    Do[If[Total[aa] > dim, Continue[]];
     Do[tot += (Times @@ (cv^(-# - 1) & /@ aa)) ordFactor[aa] (ordFactor[dd]/m!) (Times @@ (-cv^(1 - #) & /@ dd)) clf[h, Join[aa, dd]], {dd, msets[m, 2, dim - Total[aa]]}],
     {aa, msets[n, 0, dim]}], {m, 0, 3 h - 3 + n}];
   Together[tot]];
hodgeNormN[g_, nw_, a_List] := hodgeNormN[g, nw, a] = Expand[Total[Table[cl[[2]] HI[g, cl[[1, 1]], cl[[1, 2]], a], {cl, hodgeClassN[g, nw]}]]];
clWK = Function[{h, a}, HI[h, 0, 0, a]];                                                   (* point class: Witten-Kontsevich numbers *)
usCY = usubOf[{1 + x, -x, -1}]; clCY = Function[{h, a}, hodgeNormU4[h, a] /. usCY];       (* the CY triple of the task *)
usN = usubOf[{2, 3, -7}]; clN = Function[{h, a}, hodgeNormU4[h, a] /. usN];                (* three non-CY weights *)
us4 = {u[1] -> 1/2, u[2] -> -1/3, u[3] -> 1/5, u[4] -> 1/7}; cl4 = Function[{h, a}, hodgeNormN[h, 4, a] /. us4];   (* four weights (2, -3, 5, 7) *)
Print["regression (CY triple, c = 1+x): Z[2,1], Z[3,1] = ", {ZG[2, 1, 1 + x, clCY], ZG[3, 1, 1 + x, clCY]}, ";  Z[2,0] = ", Factor[ZG[2, 0, 1 + x, clCY]]];
Print["Witten-Kontsevich vertex, c = 1+x: Z[1,1], Z[1,2] = ", {Factor[ZG[1, 1, 1 + x, clWK]], Factor[ZG[1, 2, 1 + x, clWK]]}, " (anomaly Int_{M_{1,1}} 1 = 0);  Z[2,0..2] = ", Table[Factor[ZG[2, n, 1 + x, clWK]], {n, 0, 2}], ";  Z[3,0..2] = ", Table[Factor[ZG[3, n, 1 + x, clWK]], {n, 0, 2}]];
Print["non-CY weights (2, 3, -7), c = 1+x: Z[2,1], Z[2,2], Z[3,1], Z[3,2] = ", {ZG[2, 1, 1 + x, clN], ZG[2, 2, 1 + x, clN], ZG[3, 1, 1 + x, clN], ZG[3, 2, 1 + x, clN]}, ";  Z[2,0] = ", Factor[ZG[2, 0, 1 + x, clN]]];
Print["  genus 1: c Z[1,1] = ", Factor[(1 + x) ZG[1, 1, 1 + x, clN]], " = Int_{M_{1,1}} Hodge_1 = -(1/24)(1/2 + 1/3 - 1/7): ", Together[(1 + x) ZG[1, 1, 1 + x, clN] + (1/24) (1/2 + 1/3 - 1/7)] === 0, ";  c^2 Z[1,2] = Z[1,1]: ", Together[(1 + x)^2 ZG[1, 2, 1 + x, clN] - ZG[1, 1, 1 + x, clN]] === 0];
Print["CY triple but c = 5: Z[2,1], Z[2,2], Z[3,1] = ", {ZG[2, 1, 5, clCY], ZG[2, 2, 5, clCY], ZG[3, 1, 5, clCY]}, ";  c = y symbolic: Z[2,1], Z[2,2], Z[3,1] = ", {ZG[2, 1, y, clCY], ZG[2, 2, y, clCY], ZG[3, 1, y, clCY]}];
Print["  c = y symbolic, genus 1: y Z[1,1] - Int_{M_{1,1}} Hodge_1 = ", Together[y ZG[1, 1, y, clCY] - clCY[1, {0}]], ";  y^2 Z[1,2] - Z[1,1] = ", Together[y^2 ZG[1, 2, y, clCY] - ZG[1, 1, y, clCY]], ";  Z[1,1] = ", Factor[ZG[1, 1, y, clCY]]];
Print["four weights (2, -3, 5, 7), c = y: Z[2,1], Z[2,2], Z[3,1] = ", {ZG[2, 1, y, cl4], ZG[2, 2, y, cl4], ZG[3, 1, y, cl4]}, ";  Z[2,0] = ", Factor[ZG[2, 0, y, cl4]], ";  y Z[1,1] = Int_{M_{1,1}} Hodge_1 = -(1/24) Sum 1/w_a: ", Together[y ZG[1, 1, y, cl4] + (1/24) (1/2 - 1/3 + 1/5 + 1/7)] === 0];
Print["done in ", Round[AbsoluteTime[] - tt], " s"];
