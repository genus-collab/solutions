(* review_string.wls -- checks for the string-equation proof of Z[h, n] = 0 (h >= 2, n >= 1).
   Z[h, n] = Sum_m 1/m! < prod_{i<=n} phi(psi_i) prod_{j<=m} T(psi_j) >_h with phi = 1/(c - psi), T = -psi^2/(c - psi),
   for an arbitrary rank-one vertex theory corr[h, a_List] (Hodge class with arbitrary weights, or Witten-Kontsevich).
   Read-only on the project files. *)
SetDirectory[DirectoryName[$InputFileName]];
tt = AbsoluteTime[];
Get["hvc_lib.m"];
Print["library loaded in ", Round[AbsoluteTime[] - tt], " s"];

msets[k_, lo_, hi_, smax_] := Module[{rec}, rec[kk_, l_, s_] := If[kk == 0, {{}}, Join @@ Table[Prepend[#, a] & /@ rec[kk - 1, a, s - a], {a, l, Min[hi, s]}]]; rec[k, lo, smax]];
ordFactor[lst_] := Length[lst]!/Times @@ (Last /@ Tally[lst]!);
(* the vertex functional for any correlator function corr[h, a] and any leaf parameter cc *)
ZG[h_, n_, corr_, cc_] := Module[{tot = 0, dim, phiC, TC},
   phiC[k_] := 1/cc^(k + 1); TC[d_] := -1/cc^(d - 1);
   Do[dim = 3 h - 3 + n + m; If[dim < 2 m, Continue[]];
    Do[Module[{sa = Total[A]},
      Do[tot += ordFactor[A] (ordFactor[D]/m!) (Times @@ (phiC /@ A)) (Times @@ (TC /@ D)) corr[h, Join[A, D]], {D, msets[m, 2, dim - sa, dim - sa]}]], {A, msets[n, 0, dim, dim]}], {m, 0, 3 h - 3 + n}];
   Together[tot]];
hodgeCorr[w_List][h_, a_List] := hodgeNormU4[h, a] /. usubOf[w];
wkCorr[h_, a_List] := If[2 h - 2 + Length[a] > 0 && Total[a] == 3 h - 3 + Length[a], wk[h, a], 0];

Print["=== 0. the two algebraic lemmas ==="];
ts[z_] := (c s - z^2)/(c - z);
Print["Lemma 2: (t_s(z) - t_s(0))/z - (1 + (s - c)/(c - z)) = ", Simplify[(ts[z] - ts[0])/z - (1 + (s - c)/(c - z))]];
Print["t_s = c s phi + T: ", Simplify[ts[z] - (c s/(c - z) - z^2/(c - z))]];
qT = -c z/(c - z); ord = 8;
flow = Normal[Series[Exp[eps/z] qT, {z, 0, ord}]];
flowPlus = Sum[Coefficient[flow, z, k] z^k, {k, 0, ord}];     (* nonnegative part in z *)
qs = Normal[Series[c (c (1 - Exp[eps/c]) - z)/(c - z), {z, 0, ord}]];
Print["string orbit: [e^(eps/z) q_T]_+ = q_s with s = c(1 - e^(eps/c)), to z^", ord, ": ", Simplify[flowPlus - qs] === 0];

Print["=== 1. Calabi-Yau weights {1+x, -x, -1}, c = 1 + x (the case of the claim) ==="];
wB = {1 + x, -x, -1}; cB = 1 + x;
tau0 = -(1/24) Total[1/wB];
Print["<tau_0>_1 = Int_{M_{1,1}} Hodge_1 = ", Factor[tau0], "; from the library: ", Factor[hodgeCorr[wB][1, {0}]], "; equal: ", Together[tau0 - hodgeCorr[wB][1, {0}]] === 0];
Do[Print["Z[1, ", n, "] = ", Factor[ZG[1, n, hodgeCorr[wB], cB]], ";  (n-1)! c^(1-2n) <tau_0>_1 = ", Factor[(n - 1)! cB^(1 - 2 n) tau0], ";  equal: ", Together[ZG[1, n, hodgeCorr[wB], cB] - (n - 1)! cB^(1 - 2 n) tau0] === 0], {n, 1, 4}];
Do[Print["Z[2, ", n, "] = ", ZG[2, n, hodgeCorr[wB], cB]], {n, 1, 3}];
Pg[h_, ee1_, ee2_] := SeriesCoefficient[Series[ee1 ee2 t^2/(4 Sinh[ee1 t/2] Sinh[ee2 t/2]), {t, 0, 2 h}], 2 h];
gapFormula[h_, ee1_, ee2_] := (-1)^(h - 1) (2 h - 3)! Pg[h, ee1, ee2]/((ee1 ee2)^h (ee1 + ee2)^(3 h - 3));
Print["Z[2, 0] = ", Factor[ZG[2, 0, hodgeCorr[wB], cB]], ";  (-1)^(h-1)(2h-3)! P_h/((e1 e2)^h c^(3h-3)) at h = 2 = ", Factor[gapFormula[2, x, 1]], ";  equal: ", Together[ZG[2, 0, hodgeCorr[wB], cB] - gapFormula[2, x, 1]] === 0];
Print["time ", Round[AbsoluteTime[] - tt], " s"];

Print["=== 2. generic weights {w1, w2, w3} (not Calabi-Yau) and an independent leaf parameter c ==="];
wG = {w1, w2, w3}; tau0G = -(1/24) Total[1/wG];
Do[Print["Z[1, ", n, "] - (n-1)! c^(1-2n) <tau_0>_1 = ", Together[ZG[1, n, hodgeCorr[wG], c] - (n - 1)! c^(1 - 2 n) tau0G]], {n, 1, 3}];
Do[Print["Z[2, ", n, "] = ", ZG[2, n, hodgeCorr[wG], c]], {n, 1, 2}];
Print["Z[2, 0] (generic weights) = ", Factor[ZG[2, 0, hodgeCorr[wG], c]]];
Print["time ", Round[AbsoluteTime[] - tt], " s"];

Print["=== 3. Witten-Kontsevich vertex (trivial CohFT), symbolic c ==="];
Do[Print["Z^WK[1, ", n, "] = ", ZG[1, n, wkCorr, c]], {n, 1, 3}];
Do[Print["Z^WK[2, ", n, "] = ", ZG[2, n, wkCorr, c]], {n, 1, 2}];
Do[Print["Z^WK[3, ", n, "] = ", ZG[3, n, wkCorr, c]], {n, 1, 2}];
Print["Z^WK[2, 0] = ", ZG[2, 0, wkCorr, c], ",  Z^WK[3, 0] = ", ZG[3, 0, wkCorr, c], "  (the orbit invariants of the point theory)"];
Print["total time ", Round[AbsoluteTime[] - tt], " s"];
