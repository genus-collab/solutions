(* hvc_normalised_R.wls -- the R-matrix in the normalised idempotent basis at the conifold point, from the flat-component rows:
   R(z) e_i = sum_j R_ji(z) e_j with e_j = (1, L_j) in the flat basis {1, H}; Rtilde_ji = R_ji Sqrt[Delta_j/Delta_i].
   Checks the symplectic condition and the t -> 0 limits in the point frame and in mode B.  Usage: ... <Bm|point> *)
SetDirectory[DirectoryName[$InputFileName]];
tag = $ScriptCommandLine[[2]];
Get["hvc_lib.m"]; Get["hvc_factors_g3_" <> tag <> ".mx"];
S0 = Together[SsumQ]; P0 = Together[PprodQ]; sq = ((1 + x)/x) Sqrt[1 + 4 x (1/Q - 1)/(1 + x)^2]; U1 = (S0 + sq)/2; U2 = (S0 - sq)/2;
ev[C[a_, b_], root_] := a + b root;
NT = 4; ser[e_] := Series[e /. Q -> Exp[t], {t, 0, NT}];
(* rows at the two roots: r0[i][m], r1[i][m] with i = 1 (singular), 2 (regular); L_i = r1[i][0]; Delta_i *)
Do[r0[i] = Table[ser[ev[ratCN[0, 1][m], {U1, U2}[[i]]]], {m, 0, 6}]; r1[i] = Table[ser[ev[ratCN[1, 1][m], {U1, U2}[[i]]]], {m, 0, 6}]; L[i] = r1[i][[1]]; Dl[i] = ser[ev[DeltaCN[1], {U1, U2}[[i]]]], {i, 1, 2}];
Print["L_1, L_2 leading: ", Normal[L[1]] // Factor, " ; ", Normal[L[2]] // Factor];
(* idempotent components of R(z) e_i, as series in z (coefficients series in t) *)
Rc[1, i_] := Table[Together /@ ((r1[i][[m + 1]] - L[2] r0[i][[m + 1]])/(L[1] - L[2])), {m, 0, 6}];
Rc[2, i_] := Table[Together /@ ((L[1] r0[i][[m + 1]] - r1[i][[m + 1]])/(L[1] - L[2])), {m, 0, 6}];
Do[Rm[j, i] = Rc[j, i], {j, 1, 2}, {i, 1, 2}];
ordS[s_] := Module[{nn = Normal[s]}, If[nn === 0 || Together[nn] === 0, "zero", Exponent[Together[nn] t^40, t, Min] - 40]];
lead0[s_] := Factor[Normal[Series[s, {t, 0, 0}]]];
Print["=== ", tag, " ===  orders in t of the idempotent-basis entries R_ji[z^m], m = 0..6:"];
Do[Print["  R_", j, i, ": ", ordS /@ Rm[j, i]], {j, 1, 2}, {i, 1, 2}];
(* normalised entries: Rtilde_ji = R_ji Sqrt[Delta_j/Delta_i]; Delta_1 ~ A/t^2, Delta_2 ~ B *)
rat12 = Series[Sqrt[Dl[1]/Dl[2] t^2], {t, 0, NT}]/t;   (* Sqrt[Delta_1/Delta_2] = (..)/t *)
rat21 = Series[Sqrt[Dl[2]/Dl[1] /t^2], {t, 0, NT}] t;
Print["Sqrt[Delta_1/Delta_2] leading: ", Normal[rat12] // Factor];
Rt[1, 1] = Rm[1, 1]; Rt[2, 2] = Rm[2, 2]; Rt[1, 2] = Table[Series[Rm[1, 2][[m]] rat12, {t, 0, NT - 1}], {m, 1, 7}]; Rt[2, 1] = Table[Series[Rm[2, 1][[m]] rat21, {t, 0, NT - 1}], {m, 1, 7}];
Print["normalised entries at t = 0, as series in z (coefficients z^0..z^6):"];
Do[Print["  Rtilde_", j, i, "(z; t=0) = ", Table[lead0[Rt[j, i][[m]]], {m, 1, 7}]], {j, 1, 2}, {i, 1, 2}];
hB = hCoeffs[{1 + x, -x, -1}, 1, 6];
Print["  h_wB(z) coefficients: ", Factor /@ hB, ";  h_wB^-2: ", Factor /@ CoefficientList[Normal[Series[(Sum[hB[[k + 1]] z^k, {k, 0, 6}])^-2, {z, 0, 6}]], z]];
(* symplectic check R(z) R(-z)^T = 1 for the full normalised matrix, order by order in z, to O(t^2) *)
RtS[j_, i_] := Sum[Rt[j, i][[m + 1]] z^m, {m, 0, 6}];
M = Table[RtS[j, i], {j, 1, 2}, {i, 1, 2}];
prod = Table[Normal[Series[Sum[M[[j, k]] (M[[i, k]] /. z -> -z), {k, 1, 2}] - KroneckerDelta[i, j], {z, 0, 6}]], {j, 1, 2}, {i, 1, 2}];
Print["symplectic defect R(z)R(-z)^T - 1 (coefficient series in t truncated), max |order|: ", Map[If[# === 0, "0", Normal[Series[#, {t, 0, 2}]] // Simplify] &, prod, {2}]];
