(* hvc_tcorrections.wls -- z-structure of the t^k corrections (k = 1, 2, 3) of the identity row r_{0,1}(z; t) of R' at the singular idempotent
   (mode B, genus-4 factors, z^9); expectation from the Euler-integral representation of the thimble integral,
   tau^(a+b+1) B(a+1,b+1) 2F1(b+1, b+1; a+b+2; tau), a = eps1/z, b = eps2/z, tau = 1 - e^-t: poles only at z = -c/j, j = 1..k+1. *)
SetDirectory[DirectoryName[$InputFileName]];
Get["hvc_lib.m"]; Get["hvc_factors_g4_Bm.mx"];
S0 = Together[SsumQ]; P0 = Together[PprodQ]; sq = ((1 + x)/x) Sqrt[1 + 4 x (1/Q - 1)/(1 + x)^2]; U1 = (S0 + sq)/2;
c = 1 + x;
r0 = Table[Series[(ratCN[0, 1][m][[1]] + ratCN[0, 1][m][[2]] U1) /. Q -> Exp[t], {t, 0, 3}], {m, 0, 9}];
Do[Module[{co, gen, pade, den},
   co = Table[If[m == 0, 0, Factor[SeriesCoefficient[r0[[m + 1]], k]]], {m, 0, 9}];
   gen = Sum[co[[m + 1]] z^m, {m, 0, 9}];
   pade = Factor[PadeApproximant[gen, {z, 0, {4, 4}}]];
   den = FactorList[Denominator[pade]];
   Print["t^", k, ": Pade [4/4] denominator factors ", Select[den, ! FreeQ[#[[1]], z] &], "; z^9 reproduced: ", Simplify[SeriesCoefficient[Series[pade, {z, 0, 9}], 9] - co[[10]]] === 0, "; numerator ", Numerator[pade]]], {k, 1, 3}];
(* the H-row regular part and the dilaton series: same poles? *)
r1 = Table[Series[(ratCN[1, 1][m][[1]] + ratCN[1, 1][m][[2]] U1) /. Q -> Exp[t], {t, 0, 2}], {m, 1, 9}];
Do[Module[{co = Table[Factor[SeriesCoefficient[r1[[m]], k]], {m, 1, 9}], gen, pade},
   gen = Sum[co[[m]] z^m, {m, 1, 9}]; pade = Factor[PadeApproximant[gen, {z, 0, {4, 4}}]];
   Print["H-row, t^", k, ": denominator factors ", Select[FactorList[Denominator[pade]], ! FreeQ[#[[1]], z] &], "; z^9 ok: ", Simplify[SeriesCoefficient[Series[pade, {z, 0, 9}], 9] - co[[9]]] === 0]], {k, 0, 2}];
