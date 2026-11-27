(* review_string_cor.wls -- corollaries and necessity checks for the string-equation proof. Read-only on the project files. *)
SetDirectory[DirectoryName[$InputFileName]];
tt = AbsoluteTime[];
Get["hvc_lib.m"];
msets[k_, lo_, hi_, smax_] := Module[{rec}, rec[kk_, l_, s_] := If[kk == 0, {{}}, Join @@ Table[Prepend[#, a] & /@ rec[kk - 1, a, s - a], {a, l, Min[hi, s]}]]; rec[k, lo, smax]];
ordFactor[lst_] := Length[lst]!/Times @@ (Last /@ Tally[lst]!);
(* generalised functional: n legs 1/(u - psi), M string insertions 1, dilaton leaves -psi^2/(c - psi) *)
ZGen[h_, n_, M_, corr_, u_, cc_] := Module[{tot = 0, dim, phiC, TC},
   phiC[k_] := 1/u^(k + 1); TC[d_] := -1/cc^(d - 1);
   Do[dim = 3 h - 3 + n + M + m; If[dim < 2 m, Continue[]];
    Do[Module[{sa = Total[A]},
      Do[tot += ordFactor[A] (ordFactor[D]/m!) (Times @@ (phiC /@ A)) (Times @@ (TC /@ D)) corr[h, Join[A, ConstantArray[0, M], D]], {D, msets[m, 2, dim - sa, dim - sa]}]], {A, msets[n, 0, dim, dim]}], {m, 0, 3 h - 3 + n + M}];
   Together[tot]];
hodgeCorr[w_List][h_, a_List] := hodgeNormU4[h, a] /. usubOf[w];
wkCorr[h_, a_List] := If[2 h - 2 + Length[a] > 0 && Total[a] == 3 h - 3 + Length[a], wk[h, a], 0];
wB = {1 + x, -x, -1}; cB = 1 + x;

Print["=== 1. corollary of string + dilaton: Sum_m 1/m! <tau_0^M prod T>_h = (-1)^M (2h-3+M)!/((2h-3)! c^M) Z[h,0], i.e. F_h(a + T) = (1 + a/c)^(2-2h) Z[h,0] ==="];
Do[Print["h = 2, M = ", M, ": ", Together[ZGen[2, 0, M, hodgeCorr[wB], cB, cB] - (-1)^M (1 + M)!/(1! cB^M) ZGen[2, 0, 0, hodgeCorr[wB], cB, cB]] === 0], {M, 1, 3}];
Print["h = 3, M = 1: ", Together[ZGen[3, 0, 1, hodgeCorr[wB], cB, cB] - (-1) (4)!/(3! cB) ZGen[3, 0, 0, hodgeCorr[wB], cB, cB]] === 0, "   (", Round[AbsoluteTime[] - tt], " s)"];
Print["mixed: Sum_m 1/m! <phi tau_0 prod T>_2 = ", ZGen[2, 1, 1, hodgeCorr[wB], cB, cB], ",  <phi^2 tau_0 prod T>_2 = ", ZGen[2, 2, 1, hodgeCorr[wB], cB, cB], ",  <phi tau_0^2 prod T>_2 = ", ZGen[2, 1, 2, hodgeCorr[wB], cB, cB]];

Print["=== 2. necessity of matching parameters: legs 1/(u - psi) with u != c, dilaton leaves -psi^2/(c - psi) ==="];
Print["genus 2, n = 1, Hodge CY weights: ", Factor[ZGen[2, 1, 0, hodgeCorr[wB], u, cB]]];
Print["genus 2, n = 1, WK: ", Factor[ZGen[2, 1, 0, wkCorr, u, c]], ";  genus 2, n = 2, WK: ", Factor[ZGen[2, 2, 0, wkCorr, u, c]]];

Print["=== 3. winding-1 one-leg vertex: <phi>_h = c^(1-h) (e1 e2)^(-h) (2^(2h-1)-1) |B_2h| / (2^(2h-1) (2h)!)  (Marino-Vafa p_1 coefficient = 1/(2 sin(lambda/2)), i.e. lambda_h psi^(2h-2)) ==="];
J1[h_, w_, cc_] := Together[hodgeCorr[w][h, #] /. {} -> 0 & /@ {}; Module[{tot = 0}, Do[tot += hodgeCorr[w][h, {k}]/cc^(k + 1), {k, 0, 3 h - 2}]; Together[tot]]];
Do[Print["h = ", h, ": ", Together[J1[h, wB, cB] - cB^(1 - h) (x)^(-h) (2^(2 h - 1) - 1) Abs[BernoulliB[2 h]]/(2^(2 h - 1) (2 h)!)] === 0], {h, 1, 3}];
Print["generic weights, h = 2: <phi>_2 depends on the weights only through (w1 w2 w3)^-2 and the leg weight w1 = c: ", Factor[J1[2, {c, w2, w3}, c]]];

Print["=== 4. the orbit invariant of the point theory, Z^WK[h, 0], h = 2..4 ==="];
Do[Print["h = ", h, ": ", ZGen[h, 0, 0, wkCorr, c, c]], {h, 2, 4}];
Print["total time ", Round[AbsoluteTime[] - tt], " s"];
