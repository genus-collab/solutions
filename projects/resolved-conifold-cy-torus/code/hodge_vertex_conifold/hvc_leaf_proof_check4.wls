(* hvc_leaf_proof_check4.wls -- check M1 (and M16 at rational x) of leaf_proof/leaf_structure_proof.md s.12.3: the frame-B factors are
   rebuilt from scratch by the code (setupWick, rmatrixAt0, twistRows of rc_lib.m / hvc_lib.m; not the stored .mx files) at ZMAX = ZM
   and a rational x = x0, and then
   (a) rho r_{0,1}[m] = [z^m] of the KEY FORMULA (c/(c + z)) (tau/t) 2F1(1/z + 1, 1/z + 1; c/z + 2; tau) Exp[-dphi/z], m <= ZM, t^0..t^KT,
       with rho = (Delta_1 t^2/(-x c^3))^(-1/2) from the rebuilt Delta_1;
   (b) V5, the second QDE row (D_i r_{1,i})[m-1] = L_j (r_{1,i}[m] - L_i r_{0,i}[m]) at both idempotents, m = 1..ZM, t^0..t^KT;
   (c) V6, symplecticity Sum_{p+q=n} (-1)^q <v_i[p], v_j[q]> = delta_{n0} delta_ij Delta_i, exactly in Q, n <= ZM, all i, j.
   Rows of idempotent 2 are evaluated at U1 (conjugate encoding). Usage: wolframscript -file hvc_leaf_proof_check4.wls x0 [ZM] [KT]
   Output: hvc_leaf_proof_check4_out_<x0>.txt (rewritten after every line) *)
SetDirectory[DirectoryName[$InputFileName]];
args = Rest[$ScriptCommandLine]; x0 = ToExpression[args[[1]]];
ZM = If[Length[args] >= 2, ToExpression[args[[2]]], 12]; KT = If[Length[args] >= 3, ToExpression[args[[3]]], 8];
Get["hvc_lib.m"]; t0 = AbsoluteTime[];
outFile = "hvc_leaf_proof_check4_out_" <> StringReplace[ToString[InputForm[x0]], {"/" -> "_", "-" -> "m"}] <> ".txt";
out = {}; emit[s_] := (AppendTo[out, s]; Print[s]; Export[outFile, StringRiffle[out, "\n"], "Text"]);
tm[] := " [" <> ToString[Round[AbsoluteTime[] - t0]] <> " s]";
c0 = 1 + x0;
emit["x = " <> ToString[InputForm[x0]] <> ", ZMAX = " <> ToString[ZM] <> ", t-orders 0.." <> ToString[KT]];
setupWick[ZM]; emit["wickDirect ready" <> tm[]];
rmatrixAt0[{a1 -> x0, a2 -> 1, e1 -> x0, e2 -> 1}]; emit["R-matrix at alpha0 = 0 ready" <> tm[]];
twistRows[hCoeffs[{c0, -x0, -1}, -1, ZM]]; emit["rows twisted by h_w^-1 (frame B)" <> tm[]];
tidy[s_SeriesData] := MapAt[Together, s, 3]; tidy[s_] := Together[s];
sc[s_SeriesData, k_] := SeriesCoefficient[s, k]; sc[s_, k_] := If[k == 0, s, 0];
evC[cc_C] := cc[[1]] + cc[[2]] U;
S0 = Together[SsumQ]; sq = (c0/x0) Sqrt[1 + 4 x0 (1/Q - 1)/c0^2]; U1 = (S0 + sq)/2;
ser[cc_C, n_] := tidy[Series[evC[cc] /. U -> U1 /. Q -> Exp[t], {t, 0, n}]];

(* (a) the KEY FORMULA against the rebuilt singular row *)
tauS = Normal[Series[1 - Exp[-t], {t, 0, KT + 3}]];
sst = Normal[tidy[Series[(1 - (U1 /. Q -> Exp[t]))/(1 - Exp[-t]), {t, 0, KT + 2}]]];
dphi = Integrate[Normal[tidy[Series[1 - x0 sst, {t, 0, KT + 1}]]] /. t -> uu, {uu, 0, t}];
F21 = Sum[(tauS/z)^k Product[(1 + j z)^2, {j, 1, k}]/(k! Product[c0 + j z, {j, 2, k + 1}]), {k, 0, KT + 1}];
Gz = MapAt[Together, Series[(c0/(c0 + z)) (tauS/t) F21 Exp[-dphi/z], {t, 0, KT}], 3];
Gk = Table[Together[SeriesCoefficient[Gz, k]], {k, 0, KT}];
del1 = ser[DeltaCN[1], KT + 1]; del2 = ser[DeltaCN[2], KT + 1];
rho = tidy[Series[(Normal[ser[DeltaCN[1], KT + 2]] t^2/(-x0 c0^3))^(-1/2), {t, 0, KT}]];
r01 = Table[ser[ratCN[0, 1][m], KT], {m, 0, ZM}];
okA = Table[Module[{l = tidy[Series[Normal[rho] Normal[r01[[m + 1]]], {t, 0, KT}]]},
     And @@ Table[Together[sc[l, k] - Coefficient[Normal[Series[Gk[[k + 1]], {z, 0, ZM}]], z, m]] === 0, {k, 0, KT}]], {m, 0, ZM}];
emit["(a) rho r_{0,1}[m] = [z^m] KEY FORMULA, all t^k (k <= " <> ToString[KT] <> "), m = 0.." <> ToString[ZM] <> ": " <> ToString[okA] <>
  "; sample [t^1] rho r_{0,1}[" <> ToString[ZM] <> "] = " <> ToString[InputForm[Together[sc[tidy[Series[Normal[rho] Normal[r01[[ZM + 1]]], {t, 0, KT}]], 1]]]] <> tm[]];

(* (b) V5 at both idempotents *)
L1s = ser[LCN[1], KT + 1]; L2s = ser[LCN[2], KT + 1];
rows[i_] := rows[i] = Table[{ser[ratCN[0, i][m], KT], ser[ratCN[1, i][m], KT]}, {m, 0, ZM}];
Dop[f_, dl_] := tidy[D[f, t] - f D[dl, t]/(2 dl)];
qde[i_] := Module[{dl = If[i == 1, del1, del2], Li = If[i == 1, L1s, L2s], Lj = If[i == 1, L2s, L1s], rr = rows[i]},
   Table[Normal[tidy[Dop[rr[[m, 2]], dl] - Lj (rr[[m + 1, 2]] - Li rr[[m + 1, 1]])]] === 0, {m, 1, ZM}]];
emit["(b) second QDE row, idempotent 1 (m = 1.." <> ToString[ZM] <> "): " <> ToString[qde[1]] <> "; idempotent 2: " <> ToString[qde[2]] <> tm[]];

(* (c) V6 exactly in Q *)
pairC[i_, j_, p_, q_] := cadd @@ Flatten[Table[cscale[etainvCur[[a + 1, b + 1]], cmul[ratCN[a, i][p], ratCN[b, j][q]]], {a, 0, 1}, {b, 0, 1}]];
symp = Table[Module[{s = cadd @@ Table[cscale[(-1)^q, pairC[i, j, n - q, q]], {q, 0, n}], target},
     target = If[n == 0 && i == j, DeltaCN[i], C[0, 0]];
     Together[s[[1]] - target[[1]]] === 0 && Together[s[[2]] - target[[2]]] === 0], {i, 1, 2}, {j, 1, 2}, {n, 0, ZM}];
emit["(c) symplecticity <v_i(w), v_j(-w)> = delta_ij Delta_i, exact in Q, w^n for n <= " <> ToString[ZM] <> ", all i, j: " <> ToString[Union[Flatten[symp]]] <> tm[]];
emit["total time " <> ToString[Round[AbsoluteTime[] - t0]] <> " s"];
