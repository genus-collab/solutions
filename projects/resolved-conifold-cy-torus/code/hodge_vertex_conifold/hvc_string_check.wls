(* hvc_string_check.wls -- (1) the string-equation argument for Z[h,n] = 0: Phi_h(s) = Sum_n (c s)^n Z[h,n]/n! obeys (c - s) Phi_h'(s) = c delta_{h1} A_1,
   A_1 = Int_{M_{1,1}} Hodge_1 = -(1/24) Sum_a 1/w_a; check the genus-1 prediction Z[1,n] = (n-1)! A_1 / c^(2n-1) against the computed values,
   and the identities behind the argument.  (2) the O(t), O(t^2) corrections of the identity row r_{0,1}(z; t) at the singular idempotent
   (genus-4 factors, z^9): Pade reconstruction of their z-dependence. *)
SetDirectory[DirectoryName[$InputFileName]];
Get["hvc_lib.m"];
c = 1 + x;
phi[z_] := 1/(c - z); T[z_] := -z^2/(c - z); ts[s_, z_] := c s phi[z] + T[z];
Print["t_s(z) = (c s - z^2)/(c - z): ", Simplify[ts[s, z] - (c s - z^2)/(c - z)] === 0, "; d/ds t_s = c phi = 1 + z phi: ", Simplify[D[ts[s, z], s] - 1 - z phi[z]] === 0,
  "; D t_s + z phi = s phi with D f = (f(z) - f(0))/z: ", Simplify[(ts[s, z] - ts[s, 0])/z + z phi[z] - s phi[z]] === 0];
A1 = -(1/24) (1/c - 1/x - 1);
Print["A_1 = Int_{M_{1,1}} Hodge_1 = ", Factor[A1], " ; from hodgeNormU4: ", Factor[hodgeNormU4[1, {0}] /. usubOf[{c, -x, -1}]], " equal: ", Simplify[A1 - (hodgeNormU4[1, {0}] /. usubOf[{c, -x, -1}])] === 0];
Print["genus-1 prediction Z[1,n] = (n-1)! A_1/c^(2n-1) versus the computed Z[1,1], Z[1,2], Z[1,3]: ",
  Table[Factor[(n - 1)! A1/c^(2 n - 1)], {n, 1, 3}], " vs ", {(1 + x + x^2)/(24 x (1 + x)^2), (1 + x + x^2)/(24 x (1 + x)^4), (1 + x + x^2)/(12 x (1 + x)^6)}];
Print["agree: ", And @@ Table[Simplify[(n - 1)! A1/c^(2 n - 1) - {(1 + x + x^2)/(24 x (1 + x)^2), (1 + x + x^2)/(24 x (1 + x)^4), (1 + x + x^2)/(12 x (1 + x)^6)}[[n]]] === 0, {n, 1, 3}]];
(* (2) corrections of the identity row *)
Get["hvc_factors_g4_Bm.mx"];
S0 = Together[SsumQ]; P0 = Together[PprodQ]; sq = ((1 + x)/x) Sqrt[1 + 4 x (1/Q - 1)/(1 + x)^2]; U1 = (S0 + sq)/2;
r0 = Table[Series[(ratCN[0, 1][m][[1]] + ratCN[0, 1][m][[2]] U1) /. Q -> Exp[t], {t, 0, 2}], {m, 0, 9}];
Do[Module[{co = Table[Factor[SeriesCoefficient[r0[[m + 1]], k]], {m, 0, 9}], pade},
   Print["t^", k, " coefficients of r_{0,1}[z^m], m = 0..9: ", co];
   pade = Factor[PadeApproximant[Sum[co[[m + 1]] z^m, {m, 0, 9}], {z, 0, {4, 4}}]];
   Print["   Pade [4/4] in z: ", pade, "; reproduces z^9 coefficient: ", Simplify[SeriesCoefficient[Series[pade, {z, 0, 9}], 9] - co[[10]]] === 0]], {k, 1, 2}];
