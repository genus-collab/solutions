(* review_string_recursion.wls -- numerical check of the string-equation proof of Z[h, n] = 0 (h >= 2, n >= 1).
   Z[h, n] = Sum_m (1/m!) Int_{M_{h,n+m}} Hodge_h prod_{i<=n} phi(psi_i) prod_{j<=m} T(psi_{n+j}), phi = 1/(c - psi), T = -psi^2/(c - psi),
   Hodge_h with weights (c, -eps1, -eps2), eps2 = 1, eps1 = x, c = 1 + x.  Proof mechanism: c phi(psi) = 1 + psi phi(psi), so
   c Z[h, n+1] = A + B with A = Sum_m (1/m!) Int Hodge 1 prod phi prod T and B = Sum_m (1/m!) Int Hodge (psi phi)(psi_{n+1}) prod phi prod T;
   the string equation (forgetting the point with leaf 1) gives A = (n/c) Z[h, n] - B + delta_{h,1} delta_{n,0} Int_{M_{1,1}} Hodge_1,
   hence c^2 Z[h, n+1] = n Z[h, n] for h >= 2 (and for h = 1, n >= 1) and c Z[h, 1] = 0 for h >= 2.  Read-only on the project files. *)
SetDirectory[DirectoryName[$InputFileName]];
tt = AbsoluteTime[];
Get["hvc_lib.m"];
Print["library loaded in ", Round[AbsoluteTime[] - tt], " s"];
cc = 1 + x;                                   (* eps2 = 1, eps1 = x, c = eps1 + eps2 ; the symbol c is an iterator inside hodgeNormU4, so cc *)
usub = usubOf[{cc, -x, -1}];                  (* weights (c, -eps1, -eps2) of mode B *)
TC[d_] := -1/cc^(d - 1);                      (* psi^d coefficient (d >= 2) of T(psi) = -psi^2/(c - psi) *)
leafPhi[k_] := 1/cc^(k + 1);                  (* phi(psi) = 1/(c - psi) *)
leafOne[k_] := If[k == 0, 1, 0];              (* the constant leaf 1 *)
leafPsiPhi[k_] := If[k >= 1, 1/cc^k, 0];      (* psi phi(psi) = Sum_{k>=1} psi^k / c^k = (T(psi) - T(0))/(-psi) *)
msets[k_, lo_, smax_] := Module[{rec}, rec[kk_, l_, s_] := If[kk == 0, {{}}, Join @@ Table[Prepend[#, a] & /@ rec[kk - 1, a, s - a], {a, l, s}]]; rec[k, lo, smax]];
ordFactor[lst_] := Length[lst]!/Times @@ (Last /@ Tally[lst]!);
(* IntT[h, leaves, m] = Int_{M_{h,n+m}} Hodge_h prod_i leaf_i(psi_i) prod_{j<=m} T(psi_{n+j}) with m labelled dilaton points
   (legs as ordered degree tuples, dilaton degrees as multisets weighted by the number of their orderings) *)
IntT[h_, leaves_List, m_] := IntT[h, leaves, m] = Module[{n = Length[leaves], tot = 0, dim, wA},
   dim = 3 h - 3 + n + m;
   If[dim < 2 m || 2 h - 2 + n + m <= 0, Return[0]];
   Do[If[Total[aa] > dim, Continue[]];
    wA = Times @@ MapThread[#1[#2] &, {leaves, aa}];
    If[wA === 0, Continue[]];
    Do[tot += wA ordFactor[dd] (Times @@ (TC /@ dd)) (hodgeNormU4[h, Join[aa, dd]] /. usub), {dd, msets[m, 2, dim - Total[aa]]}],
    {aa, Tuples[Range[0, dim], n]}];
   Together[tot]];
V[h_, leaves_List] := V[h, leaves] = Together[Sum[IntT[h, leaves, m]/m!, {m, 0, Max[0, 3 h - 3 + Length[leaves]]}]];
Z[h_, n_] := Z[h, n] = V[h, ConstantArray[leafPhi, n]];
Aterm[h_, n_] := Aterm[h, n] = V[h, Append[ConstantArray[leafPhi, n], leafOne]];
Bterm[h_, n_] := Bterm[h, n] = V[h, Append[ConstantArray[leafPhi, n], leafPsiPhi]];

Print["=== genus 1 ==="];
Zexp1 = {(1 + x + x^2)/(24 x (1 + x)^2), (1 + x + x^2)/(24 x (1 + x)^4), (1 + x + x^2)/(12 x (1 + x)^6)};
Print["Z[1,n], n = 1..3: ", Table[Factor[Z[1, n]], {n, 1, 3}], ";  equal to the values quoted in the task: ", Table[Together[Z[1, n] - Zexp1[[n]]] === 0, {n, 1, 3}]];
anom = hodgeNormU4[1, {0}] /. usub;
Print["Int_{M_{1,1}} Hodge_1 = ", Factor[anom], ";  = -(1/24)(1/c - 1/eps1 - 1/eps2): ", Together[anom + (1/24) (1/cc - 1/x - 1)] === 0, ";  c Z[1,1] = Int_{M_{1,1}} Hodge_1: ", Together[cc Z[1, 1] - anom] === 0];
Print["recursion c^2 Z[1,n+1] = n Z[1,n] for n = 1, 2: ", Table[Together[cc^2 Z[1, n + 1] - n Z[1, n]] === 0, {n, 1, 2}]];
Print["closed form Z[1,n] = (n-1)! c^(-2(n-1)) Z[1,1], n = 2, 3: ", Table[Together[Z[1, n] - (n - 1)! cc^(-2 (n - 1)) Z[1, 1]] === 0, {n, 2, 3}]];

Print["=== genus 2 ==="];
Print["Z[2,0] = ", Factor[Z[2, 0]], ";  (-x c^3) Z[2,0] = gap_2 leading coefficient: ", Together[(-x cc^3) Z[2, 0] - gapHVC[2, x][[1]]] === 0];
Pcoef[h_] := SeriesCoefficient[x t^2/(4 Sinh[x t/2] Sinh[t/2]), {t, 0, 2 h}];
Print["(-eps1 eps2 c^3)^(h-1) Z[h,0] = (2h-3)! P_h/(eps1 eps2) at h = 2: ", Together[(-x cc^3) Z[2, 0] - 1! Pcoef[2]/x] === 0];
Print["Z[2,n], n = 1..3: ", Table[Factor[Z[2, n]], {n, 1, 3}]];
Print["recursion c^2 Z[2,n+1] - n Z[2,n], n = 0..2: ", Table[Together[cc^2 Z[2, n + 1] - n Z[2, n]], {n, 0, 2}]];

Print["=== the string-equation step, c Z[h,n+1] = A + B, A = (n/c) Z[h,n] - B (+ anomaly at h = n+1 = 1) ==="];
Do[Print["h = ", h, ", n = ", n, ":  c Z[h,n+1] = A + B: ", Together[cc Z[h, n + 1] - Aterm[h, n] - Bterm[h, n]] === 0,
   ";  A - (n/c) Z[h,n] + B - anomaly = ", Together[Aterm[h, n] - (n/cc) Z[h, n] + Bterm[h, n] - If[h == 1 && n == 0, anom, 0]],
   ";  A = ", Factor[Aterm[h, n]], ",  B = ", Factor[Bterm[h, n]]], {h, 1, 2}, {n, 0, 2}];

Print["=== term by term in the number m of dilaton leaves (h = 2, n = 0, i.e. Z[2,1]): string equation IntT[1; m] = (0/c) IntT[phi^0; m] - m IntT[psi phi; m-1] ==="];
Do[Print["  m = ", m, ":  Int_{M_{2,1+m}} Hodge 1 prod T = ", Factor[IntT[2, {leafOne}, m]], ";  -m Int_{M_{2,m}} Hodge (psi phi) prod T = ", Factor[-m IntT[2, {leafPsiPhi}, m - 1]],
   ";  Int_{M_{2,1+m}} Hodge (psi phi) prod T / m! = ", Factor[IntT[2, {leafPsiPhi}, m]/m!]], {m, 0, 5}];
Print["  A_m = IntT[1; m]/m! summed: ", Factor[Aterm[2, 0]], ";  B = ", Factor[Bterm[2, 0]], ";  A + B = c Z[2,1] = ", Together[Aterm[2, 0] + Bterm[2, 0]]];
Print["  the same for h = 2, n = 1 (Z[2,2]): IntT[phi,1; m] - (1/c) IntT[phi; m] + m IntT[phi, psi phi; m-1] for m = 0..5: ", Table[Together[IntT[2, {leafPhi, leafOne}, m] - (1/cc) IntT[2, {leafPhi}, m] + m IntT[2, {leafPhi, leafPsiPhi}, m - 1]], {m, 0, 5}]];

Print["=== genus 3 ==="];
Print["Z[3,0] = ", Factor[Z[3, 0]], ";  (-x c^3)^2 Z[3,0] = gap_3: ", Together[(-x cc^3)^2 Z[3, 0] - gapHVC[3, x][[1]]] === 0, ";  = 3! P_3/(eps1 eps2): ", Together[(-x cc^3)^2 Z[3, 0] - 3! Pcoef[3]/x] === 0];
Print["Z[3,1] = ", Factor[Z[3, 1]], ";  Z[3,2] = ", Factor[Z[3, 2]], " (", Round[AbsoluteTime[] - tt], " s)"];
Do[Print["h = 3, n = ", n, ":  c Z[3,n+1] = A + B: ", Together[cc Z[3, n + 1] - Aterm[3, n] - Bterm[3, n]] === 0, ";  A - (n/c) Z[3,n] + B = ", Together[Aterm[3, n] - (n/cc) Z[3, n] + Bterm[3, n]], ";  B = ", Factor[Bterm[3, n]]], {n, 0, 1}];
TimeConstrained[Print["Z[3,3] = ", Factor[Z[3, 3]], ";  recursion c^2 Z[3,3] - 2 Z[3,2] = ", Together[cc^2 Z[3, 3] - 2 Z[3, 2]], " (", Round[AbsoluteTime[] - tt], " s)"], 50, Print["Z[3,3] skipped (time)"]];
TimeConstrained[Print["Z[4,1] = ", Factor[Z[4, 1]], " (", Round[AbsoluteTime[] - tt], " s)"], 50, Print["Z[4,1] skipped (time)"]];
Print["done in ", Round[AbsoluteTime[] - tt], " s"];
