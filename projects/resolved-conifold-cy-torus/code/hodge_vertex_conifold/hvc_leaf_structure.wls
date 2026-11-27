(* hvc_leaf_structure.wls -- leaf structure of the t-expansion of the singular column (mode B, genus-4 factors, z^9).
   rho(t) = (Delta_1 t^2/(-x c^3))^(-1/2); q_t(psi) = Sum_{d>=1} Dil_1(d; t) psi^d; legs f_{j,b}(psi) = Sum_a E_{1j}(a, b) psi^a; loops E11(psi, psi').
   For k = 0..KT the order-t^k coefficients of rho q_t - q_0, rho f_{j,b}, rho^2 E11 are fitted as combinations of lam_j(psi) = 1/(c - j psi),
   j = 0..k+1 (products lam_i(psi) lam_l(psi') for E11), from the psi-orders 0..k+1; all remaining available psi-orders are then checked exactly.
   Usage: wolframscript -file hvc_leaf_structure.wls [book | kt5]; "book" adds stage 2 (leaf+Xi bookkeeping versus the direct vertex functional,
   genus 2 and 3, KT = 3 only); "kt5" extends the decompositions to t^5 (loops then only with the restricted ansatz).
   Output: hvc_leaf_structure_out.txt / hvc_leaf_structure_out_kt5.txt (printed), hvc_leaf_structure_data(_kt5).m (decompositions). *)
SetDirectory[DirectoryName[$InputFileName]];
Get["hvc_lib.m"]; Get["hvc_factors_g4_Bm.mx"]; tt0 = AbsoluteTime[];
S0 = Together[SsumQ]; P0 = Together[PprodQ]; sq = ((1 + x)/x) Sqrt[1 + 4 x (1/Q - 1)/(1 + x)^2]; U1 = (S0 + sq)/2;
c = 1 + x; KT = If[MemberQ[$ScriptCommandLine, "kt5"], 5, 3]; ZM = ZMAX;   (* "kt5": orders t^4, t^5 too (genus-4 range) *)
ser[cc_C, n_] := Series[(cc[[1]] + cc[[2]] U1) /. Q -> Exp[t], {t, 0, n}];
tco[cc_C, n_] := Module[{s = ser[cc, n]}, If[Head[s] === SeriesData, Table[Together[SeriesCoefficient[s, k]], {k, 0, n}], Prepend[ConstantArray[0, n], Together[s]]]];
(* rho *)
dl = ser[DeltaCN[1], KT - 2];
rhoS = Series[(Normal[dl] t^2/(-x c^3))^(-1/2), {t, 0, KT}]; rhoC = Table[Together[SeriesCoefficient[rhoS, k]], {k, 0, KT}];
rho2C = Table[Together[SeriesCoefficient[rhoS^2, k]], {k, 0, KT}];
Print["rho = ", Factor /@ rhoC];
(* leaf fits *)
lamco[j_, n_] := If[j == 0, If[n == 0, 1/c, 0], j^n/c^(n + 1)];
Mlam[K_] := Table[lamco[j, n], {n, 0, K}, {j, 0, K}];
fit1[v_List, K_] := Module[{N = Length[v] - 1, sol, res},
   If[N < K, Return[{Missing["too few orders"], 0, False}]];
   sol = Together /@ LinearSolve[Mlam[K], v[[1 ;; K + 1]]];
   res = Table[Together[Sum[sol[[j + 1]] lamco[j, n], {j, 0, K}] - v[[n + 1]]], {n, K + 1, N}];
   {Factor /@ sol, Length[res], And @@ (# === 0 & /@ res)}];
fit2[vf_, K_, avail_] := Module[{Mi = Inverse[Mlam[K]], V, sol, res},
   V = Table[vf[a, b], {a, 0, K}, {b, 0, K}];
   sol = Map[Together, Mi . V . Transpose[Mi], {2}];
   res = Table[Together[Sum[sol[[i + 1, l + 1]] lamco[i, p[[1]]] lamco[l, p[[2]]], {i, 0, K}, {l, 0, K}] - vf @@ p], {p, Select[avail, Max[#] > K &]}];
   {Map[Factor, sol, {2}], Length[res], And @@ (# === 0 & /@ res)}];
conv[a_List, b_List, k_] := Together[Sum[a[[i + 1]] b[[k - i + 1]], {i, 0, k}]];
(* (i) dilaton *)
qC = Table[tco[DilN[1][d], KT], {d, 1, ZM + 1}];
q0co[d_] := If[d == 0, 0, -1/c^(d - 1)];
dilRes = Table[Module[{v = Table[If[d == 0, 0, conv[rhoC, qC[[d]], k]] - If[k == 0, q0co[d], 0], {d, 0, ZM + 1}]}, fit1[v, k + 1]], {k, 0, KT}];
Do[Print["(i) rho q_t - q_0, t^", k, ": coefficients of lam_0..lam_", k + 1, " = ", dilRes[[k + 1, 1]], "; checks ", dilRes[[k + 1, 2]], " ok: ", dilRes[[k + 1, 3]]], {k, 0, KT}];
Print["time ", Round[AbsoluteTime[] - tt0], " s"];
(* (ii) legs *)
eC[j_, a_, b_] := eC[j, a, b] = tco[EdN[1, j][a, b], KT];
legRes = Table[Module[{fk},
     fk = Table[Table[conv[rhoC, eC[j, a, b], k], {a, 0, ZM - 1 - b}], {k, 0, KT}];
     Table[fit1[fk[[k + 1]], k + 1], {k, 0, KT}]], {j, 1, 2}, {b, 0, ZM - 1}];
Do[Print["(ii) leg j = ", j, ", b = ", b, ", t^", k, ": lam_0..lam_", k + 1, " = ", legRes[[j, b + 1, k + 1, 1]], "; checks ", legRes[[j, b + 1, k + 1, 2]], " ok: ", legRes[[j, b + 1, k + 1, 3]]],
  {j, 1, 2}, {b, 0, ZM - 1}, {k, 0, KT}];
Print["time ", Round[AbsoluteTime[] - tt0], " s"];
(* (iii) loops *)
availE = Select[Tuples[Range[0, ZM - 1], 2], Total[#] <= ZM - 1 &];
loopRes = Table[Module[{vf}, vf[a_, b_] := conv[rho2C, eC[1, a, b], k]; fit2[vf, k + 1, availE]], {k, 0, Min[KT, 3]}];
Do[Print["(iii) rho^2 E11, t^", k, ": matrix (lam_i(psi) lam_l(psi'), i, l = 0..", k + 1, ") = ", loopRes[[k + 1, 1]], "; checks ", loopRes[[k + 1, 2]], " ok: ", loopRes[[k + 1, 3]]], {k, 0, Min[KT, 3]}];
(* restricted symmetric ansatz for the loops: Sum_{1 <= i <= l, (i-1)+(l-1) <= k} e_il (lam_i(psi) lam_l(psi') + (i<->l)), solved from ALL available (a, a'),
   a + a' <= 8 (overdetermined); works also for k = 4, 5 where the full (k+2)^2 ansatz cannot be fitted *)
fitRestricted[vf_, k_, avail_] := Module[{pairs, M, rhs, rows = {}, sol, res},
   pairs = Select[Tuples[Range[1, k + 1], 2], #[[1]] <= #[[2]] && Total[#] - 2 <= k &];
   M = Table[Table[If[p[[1]] == p[[2]], p[[1]]^q[[1]] p[[2]]^q[[2]], p[[1]]^q[[1]] p[[2]]^q[[2]] + p[[2]]^q[[1]] p[[1]]^q[[2]]], {p, pairs}], {q, avail}];
   rhs = Table[Together[c^(Total[q] + 2) vf @@ q], {q, avail}];
   Do[If[MatrixRank[M[[Append[rows, i]]]] > Length[rows], AppendTo[rows, i]], {i, Length[avail]}];
   If[Length[rows] < Length[pairs], Return[{"underdetermined", Length[rows], Length[pairs]}]];
   sol = Together /@ LinearSolve[M[[rows]], rhs[[rows]]];
   res = Together /@ (M . sol - rhs);
   {Thread[pairs -> Factor /@ sol], Length[avail] - Length[pairs], And @@ (# === 0 & /@ res)}];
loopResR = Table[Module[{vf}, vf[a_, b_] := conv[rho2C, eC[1, a, b], k]; fitRestricted[vf, k, availE]], {k, 0, KT}];
Do[Print["(iii') rho^2 E11, t^", k, ", restricted ansatz (no lam_0, D <= k): ", loopResR[[k + 1, 1]], "; extra equations ", loopResR[[k + 1, 2]], " all satisfied: ", loopResR[[k + 1, 3]]], {k, 0, KT}];
Put[{rhoC, dilRes, legRes, loopRes, loopResR}, If[KT == 3, "hvc_leaf_structure_data.m", "hvc_leaf_structure_data_kt5.m"]];
Print["total ", Round[AbsoluteTime[] - tt0], " s"];

(* ===================== stage 2: bookkeeping (leaf expansion + Xi versus the direct vertex functional) ===================== *)
usub = usubMode["B", x];
msets[k_, lo_, hi_, smax_] := Module[{rec}, rec[kk_, l_, s_] := If[kk == 0, {{}}, Join @@ Table[Prepend[#, a] & /@ rec[kk - 1, a, s - a], {a, l, Min[hi, s]}]]; rec[k, lo, smax]];
ordFactor[lst_] := Length[lst]!/Times @@ (Last /@ Tally[lst]!);
hodgeN[h_, a_List] := hodgeN[h, Sort[a]] /; a =!= Sort[a];
hodgeN[h_, a_List] := hodgeN[h, a] = Together[hodgeNormU4[h, a] /. usub];
TC0[d_] := -1/c^(d - 1);
Lval[h_, a_List] := Lval[h, a] = Module[{tot = 0, r = Length[a], sa = Total[a], dim},
   Do[dim = 3 h - 3 + r + m;
    Do[tot += (ordFactor[D]/m!) (Times @@ (TC0 /@ D)) hodgeN[h, Join[a, D]], {D, msets[m, 2, dim, dim - sa]}], {m, 0, 3 h - 3 + r - sa}];
   Together[tot]];
(* Xi_h(js) by direct summation (r <= 3) *)
XiDirect[h_, js_List] := XiDirect[h, js] = Module[{r = Length[js], amax, tot = 0},
   amax = 3 h - 3 + r;
   Do[If[Total[a] <= amax, tot += (Times @@ Table[lamco[js[[i]], a[[i]]], {i, r}]) Lval[h, Sort[a]]], {a, Tuples[Range[0, amax], r]}];
   Together[tot]];
(* Xi_h with the lam_1 leaves removed by the string equation (proved): Xi(1^p, rest) = prod_{i=0}^{p-1} (J_rest + i)/c^2 Xi(rest) *)
XiRed[h_, js0_List] := Module[{js = Sort[js0], p, rest},
   p = Count[js, 1]; rest = DeleteCases[js, 1];
   If[MemberQ[rest, 0], Return[Missing["lam_0 leaf"]]];
   Together[Product[(Total[rest] + i)/c^2, {i, 0, p - 1}] If[rest === {}, Lval[h, {}], XiDirect[h, rest]]]];
Print["string-reduction sanity check (genus 2, 3; (1,2), (1,1,2), (1,3), (1,2,2)): ",
  Table[Together[XiRed[h, js] - XiDirect[h, js]] === 0, {h, 2, 3}, {js, {{1, 2}, {1, 1, 2}, {1, 3}, {1, 2, 2}}}]];
(* leaf-based prediction of the coefficients of t^(2-2h+k), k = 0..KT *)
deltaS = Sum[tt^k Sum[dilRes[[k + 1, 1, j + 1]] LL[j], {j, 0, k + 1}], {k, 1, KT}];
legS[j_, b_] := Sum[tt^k Sum[legRes[[j, b + 1, k + 1, 1, jp + 1]] LL[jp], {jp, 0, k + 1}], {k, 0, KT}];
loopS = Sum[tt^k Sum[loopRes[[k + 1, 1, i + 1, l + 1]] LL[i] LL[l], {i, 0, k + 1}, {l, 0, k + 1}], {k, 0, KT}];
truncT[e_] := Normal[Series[e, {tt, 0, KT}]];
monoToXi[h_, poly_] := Module[{cr = CoefficientRules[Expand[poly], Array[LL, 5, 0]], maxD = 0, tot = 0},
   Do[Module[{js = Flatten[Table[ConstantArray[i - 1, rr[[1, i]]], {i, 5}]], xi},
      xi = XiRed[h, js]; maxD = Max[maxD, Total[js - 1]];
      tot += rr[[2]] xi], {rr, cr}];
   {Together[tot], maxD}];
vertexLeaf[h_, legs_List, nloops_] := Module[{gen, out},
   gen = truncT[(Times @@ (legS @@@ legs)) loopS^nloops Sum[deltaS^m/m!, {m, 0, KT}]];
   out = Table[monoToXi[h, Coefficient[gen, tt, k]], {k, 0, KT}];
   {Factor /@ ((-x c^3)^(h - 1) out[[All, 1]]), out[[All, 2]]}];
(* direct vertex functional from the factor series (no leaf decomposition) *)
dilT[d_] := dilT[d] = tco[DilN[1][d], KT];
smul[a_List, b_List] := Table[Sum[a[[i + 1]] b[[k - i + 1]], {i, 0, k}], {k, 0, KT}];
one = PadRight[{1}, KT + 1];
dPow[h_] := Module[{s = Series[(Normal[dl])^(h - 1), {t, 0, 2 - 2 h + KT}]}, Table[Together[SeriesCoefficient[s, 2 - 2 h + k]], {k, 0, KT}]];
vertexDirect[h_, legs_List, nloops_] := Module[{tot = ConstantArray[0, KT + 1], nl = Length[legs], nhalf, dim},
   nhalf = nl + 2 nloops;
   Do[dim = 3 h - 3 + nhalf + m; If[dim < 2 m, Continue[]];
    Do[Module[{pd = Fold[smul, one, dilT /@ D], loc = ConstantArray[0, KT + 1], sD = Total[D]},
      Do[Do[Module[{lp = Join[A, Flatten[L], D], w},
          If[Total[lp] > dim, Continue[]];
          w = hodgeN[h, lp]; If[w === 0, Continue[]];
          loc += w Fold[smul, one, Join[Table[eC[legs[[i, 1]], A[[i]], legs[[i, 2]]], {i, nl}], Table[eC[1, L[[q, 1]], L[[q, 2]]], {q, nloops}]]]],
         {L, Tuples[Table[Tuples[Range[0, dim - sD], 2], {nloops}]]}],
        {A, Tuples[Range[0, dim - sD], nl]}];
      tot += (ordFactor[D]/m!) smul[pd, Together /@ loc]], {D, msets[m, 2, dim, dim]}],
    {m, 0, 3 h - 3 + nhalf}];
   Factor /@ smul[dPow[h], Together /@ tot]];
compare[h_, legs_, nloops_] := Module[{t1 = AbsoluteTime[], vl, vd},
   vl = vertexLeaf[h, legs, nloops]; vd = vertexDirect[h, legs, nloops];
   Print["genus ", h, ", legs ", legs, ", loops ", nloops, ": direct {t^", 2 - 2 h, "..t^", 2 - 2 h + KT, "} = ", vd,
    "; leaf+Xi prediction agrees: ", Table[Together[vl[[1, k]] - vd[[k]]] === 0, {k, KT + 1}], "; max leaf degree D per order ", vl[[2]], "   (", Round[AbsoluteTime[] - t1, 0.1], " s)"]];
If[MemberQ[$ScriptCommandLine, "book"] && KT == 3,
  Print["gap_2 = ", gapHVC[2, x][[1]], ", gap_3 = ", gapHVC[3, x][[1]]];
  Print["Conjecture M instances needed after string reduction (all leaves j >= 2, D <= 2h-3): genus 2 {2}: ", XiDirect[2, {2}], "; genus 3 {2},{3},{4},{2,2},{2,3},{2,2,2}: ", Table[XiDirect[3, js], {js, {{2}, {3}, {4}, {2, 2}, {2, 3}, {2, 2, 2}}}], "; first non-zero ones (D = 2h-2): genus 2 {3}, {2,2} nonzero: ", Together[XiDirect[2, #]] =!= 0 & /@ {{3}, {2, 2}}, "; genus 3 {5}, {2,4}, {3,3}, {2,2,3}, {2,2,2,2} nonzero: ", Together[XiRed[3, #]] =!= 0 & /@ {{5}, {2, 4}, {3, 3}, {2, 2, 3}}];
  compare[2, {}, 0]; compare[2, {{1, 0}}, 0]; compare[2, {{2, 0}}, 0]; compare[2, {{2, 1}}, 0]; compare[2, {}, 1]; compare[2, {{1, 0}, {2, 0}}, 0];
  Print["time ", Round[AbsoluteTime[] - tt0], " s"];
  compare[3, {}, 0]; compare[3, {{1, 0}}, 0]; compare[3, {{2, 0}}, 0]; compare[3, {{1, 1}}, 0]; compare[3, {{2, 1}}, 0];
  Print["time ", Round[AbsoluteTime[] - tt0], " s"];
  compare[3, {}, 1]; compare[3, {{2, 0}, {2, 0}}, 0]; compare[3, {{1, 0}, {2, 0}}, 0];
  Print["total ", Round[AbsoluteTime[] - tt0], " s"]];
