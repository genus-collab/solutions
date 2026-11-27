(* review_string3.wls -- genus-3 checks of the string-equation argument (Hodge vertex with random non-Calabi-Yau weights and random c),
   plus the n = 0 constants: Int_{M_h} Hodge_h closed form and the gap formula at genus 3 (Calabi-Yau weights). Read-only on the project files. *)
SetDirectory[DirectoryName[$InputFileName]];
tt = AbsoluteTime[];
Get["hvc_lib.m"];
Print["library loaded in ", Round[AbsoluteTime[] - tt], " s"];
msets[k_, lo_, hi_, smax_] := Module[{rec}, rec[kk_, l_, s_] := If[kk == 0, {{}}, Join @@ Table[Prepend[#, a] & /@ rec[kk - 1, a, s - a], {a, l, Min[hi, s]}]]; rec[k, lo, smax]];
ordFactor[lst_] := Length[lst]!/Times @@ (Last /@ Tally[lst]!);
ZG[h_, n_, corr_, cc_] := Module[{tot = 0, dim, phiC, TC},
   phiC[k_] := 1/cc^(k + 1); TC[d_] := -1/cc^(d - 1);
   Do[dim = 3 h - 3 + n + m; If[dim < 2 m, Continue[]];
    Do[Module[{sa = Total[A]},
      Do[tot += ordFactor[A] (ordFactor[D]/m!) (Times @@ (phiC /@ A)) (Times @@ (TC /@ D)) corr[h, Join[A, D]], {D, msets[m, 2, dim - sa, dim - sa]}]], {A, msets[n, 0, dim, dim]}], {m, 0, 3 h - 3 + n}];
   Together[tot]];
hodgeCorr[w_List][h_, a_List] := hodgeNormU4[h, a] /. usubOf[w];

Print["=== 0. string orbit: [e^(eps/z) q_T]_+ = q_s with s = c (1 - e^(eps/c)), q_T = -c z/(c - z), q_s = c (s - z)/(c - z) (double series) ==="];
K = 10; ord = 6;
qTc[j_] := If[j == 0, 0, -c^(1 - j)];                                         (* coefficient of z^j in -c z/(c - z) = -Sum_{j>=1} c^(1-j) z^j *)
flowCoeff[k_] := Sum[eps^i/i! qTc[k + i], {i, 0, K}];                           (* coefficient of z^k in [e^(eps/z) q_T]_+, truncated at eps^K *)
qsCoeff[k_] := Normal[Series[SeriesCoefficient[Series[c (c (1 - Exp[eps/c]) - z)/(c - z), {z, 0, ord}], k], {eps, 0, K}]];
Print["agreement of the z^k coefficients as series in eps to eps^", K, ", k = 0..", ord, ": ", Table[Simplify[Normal[Series[flowCoeff[k] - qsCoeff[k], {eps, 0, K}]]] === 0, {k, 0, ord}]];

Print["=== 1. Int_{M_h} Hodge_h versus (-1)^h (w1 w2 w3)^(1-h) |B_2h B_2h-2| / (4h (2h-2) (2h-2)!) for Calabi-Yau weights, h = 2, 3 ==="];
wB = {1 + x, -x, -1}; cB = 1 + x;
Ch[h_, w_] := (-1)^h (Times @@ w)^(1 - h) Abs[BernoulliB[2 h] BernoulliB[2 h - 2]]/(4 h (2 h - 2) (2 h - 2)!);
Do[Print["h = ", h, ": ", Together[hodgeCorr[wB][h, {}] - Ch[h, wB]] === 0], {h, 2, 3}];

Print["=== 2. genus 3, Hodge vertex, random non-Calabi-Yau rational weights and random c ==="];
SeedRandom[20261030];
wR = RandomInteger[{1, 40}, 3]/RandomInteger[{1, 7}, 3] {1, -1, 1}; cR = RandomInteger[{1, 30}]/RandomInteger[{1, 5}];
Print["weights ", wR, " (sum ", Total[wR], "), c = ", cR];
Print["Z[3, 1] = ", ZG[3, 1, hodgeCorr[wR], cR], "   (", Round[AbsoluteTime[] - tt], " s)"];
Print["Z[3, 2] = ", ZG[3, 2, hodgeCorr[wR], cR], "   (", Round[AbsoluteTime[] - tt], " s)"];
Print["Z[1, 4] - 3! c^-7 <tau_0>_1 = ", ZG[1, 4, hodgeCorr[wR], cR] - 3! cR^(-7) (-(1/24) Total[1/wR])];

Print["=== 3. genus 3, Calabi-Yau weights: Z[3, 0] versus (-1)^(h-1) (2h-3)! P_h/((e1 e2)^h c^(3h-3)) ==="];
Pg[h_, ee1_, ee2_] := SeriesCoefficient[Series[ee1 ee2 t^2/(4 Sinh[ee1 t/2] Sinh[ee2 t/2]), {t, 0, 2 h}], 2 h];
gapFormula[h_, ee1_, ee2_] := (-1)^(h - 1) (2 h - 3)! Pg[h, ee1, ee2]/((ee1 ee2)^h (ee1 + ee2)^(3 h - 3));
z30 = ZG[3, 0, hodgeCorr[wB], cB];
Print["Z[3, 0] = ", Factor[z30], ";  formula: ", Factor[gapFormula[3, x, 1]], ";  equal: ", Together[z30 - gapFormula[3, x, 1]] === 0, "   (", Round[AbsoluteTime[] - tt], " s)"];
Print["genus-1 orbit invariant Z[1, 0] (CY weights) = ", Factor[ZG[1, 0, hodgeCorr[wB], cB]]];
Print["total time ", Round[AbsoluteTime[] - tt], " s"];
