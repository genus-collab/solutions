(* rc_lib.m -- library for the resolved conifold x C^2 pipeline (general five-dimensional torus), genus <= 4.
   Load with Get after SetDirectory to the data directory (graphs.m, orbits.m, StableGraph-<g>.bin must be there).
   Sections 1-4 are verbatim copies of the Input cells of resolvedconifold_genus3.nb (setup, quadratic algebra, numeric
   R-matrix and graph factors, graph sum, degree-zero term, Q-reconstruction and weight interpolation); section 5 adds
   the closed-form QRR normalisation, the alpha0 = 0 R-matrix, the Laurent-series layer at Q = 1, the sample lists and
   the multivariate fit on the CY torus.  ASCII only; the weights are the symbols Alpha[0..2], Epsilon[1..2]
   (typed as backslash-bracket-Alpha / -Epsilon). *)

Get["graphs.m"]; Get["orbits.m"];

(* ============================== 1. setup: Frobenius structure, canonical coordinates, Delta (verbatim) ============================== *)
a0 = \[Alpha][0]; a1 = \[Alpha][1]; a2 = \[Alpha][2]; e1 = \[Epsilon][1]; e2 = \[Epsilon][2];
CYearly = False;
CYrepl = {\[Alpha][0] -> \[Epsilon][1] + \[Epsilon][2] - \[Alpha][1] - \[Alpha][2]};
cy[x_] := If[CYearly, x /. CYrepl, x];
Replp2 = -Expand[p (p + a0)] + p^2;
exp3 = Simplify[Expand[Expand[Expand[(t[0] + t[1] p)^3] /. p^4 -> p^2 Replp2] /. p^3 -> p Replp2] /. p^2 -> Replp2];
PrepCl = cy[Simplify[Normal[Series[((exp3/a0/(-p + a2)/(-p + a1)/e1/e2 /. p -> 0) + (exp3/(-a0)/(-p + a2)/(-p + a1)/e1/e2 /. p -> -a0))/6, {t[0], 0, 3}, {t[1], 0, 3}]]]];
eta = Table[Simplify[D[D[D[PrepCl, t[0]], t[i]], t[j]]], {i, 0, 1}, {j, 0, 1}];
etainv = Simplify[Inverse[eta]];
Prep = PrepCl + PolyLog[3, Exp[t[1]]]/(e1 e2);
Do[cch[i, j, k] = Simplify[Sum[etainv[[i + 1, m + 1]] D[D[D[Prep, t[m]], t[j]], t[k]], {m, 0, 1}] /. t[1] -> Log[Q]], {i, 0, 1}, {j, 0, 1}, {k, 0, 1}];
Leq = Together[L^2 - cch[0, 1, 1] - cch[1, 1, 1] L];
Ssum = Together[cch[1, 1, 1]];
Pprod = Together[-cch[0, 1, 1]];
RedLeq[arg_, var_] := Block[{out = Together[arg]},
   Together[PolynomialRemainder[Numerator[out], Leq /. L -> var, var]/PolynomialRemainder[Denominator[out], Leq /. L -> var, var]]];
LDer[Lv_] := Together[-D[Leq, Q]/D[Leq, L] /. L -> Lv];
Module[{Psi = {{1, L1}, {1, L2}}, etacan},
  etacan = Simplify[Transpose[Inverse[Psi]] . eta . Inverse[Psi]];
  DeltaL = RedLeq[1/etacan[[1, 1]] /. L2 -> Ssum - L1 /. L1 -> L, L]];

(* ============================== 2. arithmetic in the quadratic algebra (verbatim, series-aware cTog) ============================== *)
cTog[x_SeriesData] := x;
cTog[x_] := Together[x];
cadd[C[a_, b_], C[c_, d_]] := C[cTog[a + c], cTog[b + d]];
cadd[x__] := Fold[cadd, C[0, 0], {x}];
cscale[k_, C[a_, b_]] := C[cTog[k a], cTog[k b]];
cmul[C[a_, b_], C[c_, d_]] := C[cTog[a c - b d Pprod], cTog[a d + b c + b d Ssum]];
cmul[x__] := Fold[cmul, C[1, 0], {x}];
cconj[C[a_, b_]] := C[cTog[a + b Ssum], cTog[-b]];
cnorm[C[a_, b_]] := cTog[a^2 + a b Ssum + b^2 Pprod];
cinv[C[a_, b_]] := Module[{n = cnorm[C[a, b]]}, C[cTog[(a + b Ssum)/n], cTog[-b/n]]];
cpow[c_, 0] := C[1, 0];
cpow[c_, n_Integer?Positive] := cmul @@ ConstantArray[c, n];
cpow[c_, n_Integer?Negative] := cpow[cinv[c], -n];
genSym = L;
toC[e_] := Module[{r, num, den, a, b, c, d, v = genSym},
   r = Together[e]; num = PolynomialRemainder[Numerator[r], v^2 - Ssum v + Pprod, v];
   den = PolynomialRemainder[Denominator[r], v^2 - Ssum v + Pprod, v];
   a = Coefficient[num, v, 0]; b = Coefficient[num, v, 1]; c = Coefficient[den, v, 0]; d = Coefficient[den, v, 1];
   cmul[C[a, b], cinv[C[c, d]]]];
symPow[sym_, n_, symvals_] := symPow[sym, n, symvals] = cpow[symvals[sym], n];
termEval[t_, symvals_] := Module[{facs = If[Head[t] === Times, List @@ t, {t}], coef = 1, fac = C[1, 0]},
   Do[Which[KeyExistsQ[symvals, f], fac = cmul[fac, symvals[f]],
      MatchQ[f, Power[_, _Integer]] && KeyExistsQ[symvals, f[[1]]], fac = cmul[fac, symPow[f[[1]], f[[2]], symvals]],
      True, coef = coef f], {f, facs}];
   cscale[coef, fac]];
cEval[poly_, symvals_Association] := Module[{ex = Expand[poly], terms},
   terms = If[Head[ex] === Plus, List @@ ex, {ex}];
   cadd @@ (termEval[#, symvals] & /@ terms)];
cExpr[C[a_, b_]] := a + b genSym;
cSymQ[C[a_, b_]] := Module[{bb = If[Head[b] === SeriesData, Normal[b], b]}, bb === 0 || Together[bb] === 0];

(* ============================== 3. R-matrix at numeric weights, graph factors (verbatim + closed-form variants) ============================== *)
(* symbolic pieces computed once for a given ZMAX.
   wickZ[[m+1]] = coefficient of z^m of sum_n (2n)! (-z)^n/(n! 2^n W2^n) [x^2n] Exp[sum_{k>=3} W_k x^k/(k! z)] as a polynomial in the
   symbols Wk[k]: a term with j factors W_{k_1} ... W_{k_j}, sum k_i = 2n, has z-order n - j, so the z^m coefficient is the sum over
   n = m .. 3m and over the partitions of 2n into exactly n - m parts in [3, 2 ZMAX + 2], with the multiset factor
   prod_k 1/((k!)^{m_k} m_k!).  Enumerating the partitions directly avoids the x-series of the exponential (hopeless at z^9).
   The derivatives W_k of the superpotential are computed by computeWU only after the weights are substituted (the fully
   symbolic W_k explode for k > 14). *)
wickDirect[zmax_] := Table[Module[{tot = 0, j, parts},
     Do[j = n - m;
      parts = If[j == 0, If[n == 0, {{}}, {}], Select[IntegerPartitions[2 n, {j}], Min[#] >= 3 && Max[#] <= 2 zmax + 2 &]];   (* parts in [3, 2 zmax + 2]; the explicit part-list form of IntegerPartitions is very slow here *)
      Do[tot += (2 n)! (-1)^n/(n! 2^n Wk[2]^n) Times @@ (Wk[#[[1]]]^#[[2]]/((#[[1]]!)^#[[2]] #[[2]]!) & /@ Tally[p]), {p, parts}], {n, m, 3 m}];
     Expand[tot]], {m, 0, zmax}];
setupWick[zmax_] := (ZMAX = zmax;
   wickZ = wickDirect[zmax];
   spotU = a1 Log[1 - Q U] - a2 Log[U] + (a0 + a2) Log[1 - U];
   LofU = -(Q (a0 + a1) U + Q a2)/(1 - Q);
   etaSym = {{(a0 + a1 + a2)/(a1 (a0 + a1) a2 (a0 + a2) e1 e2), 1/((a0 + a1) (a0 + a2) e1 e2)}, {1/((a0 + a1) (a0 + a2) e1 e2), -(a0/((a0 + a1) (a0 + a2) e1 e2))}};
   etainvSym = Simplify[Inverse[etaSym]];
   wtSym[1] = {a0, a1, a2, -e1, -e2}; wtSym[2] = {-a0, a0 + a1, a0 + a2, -e1, -e2};
   "wick ready");
QrrcSym[i_, k_] := (BernoulliB[2 k + 2]/(2 k + 2)/(2 k + 1)) Sum[1/w^(2 k + 1), {w, wtSym[i]}];
(* derivatives W_k = (U d/dU)^k spot at the given weights (rational functions of U and Q; cheap even for k = 20) *)
computeWU[wsub_] := (WU[1] = Together[U D[spotU /. wsub, U]];
   Do[WU[k] = Together[U D[WU[k - 1], U]], {k, 2, 2 ZMAX + 2}]; "W ready");
(* QRR normalisation constants from the Q -> 0 limit (verbatim); returns {ffc list root 1, ffc list root 2} (ffc[0] = 1) *)
ffcLimit[wsub_, SU_, LeqU_] := Module[{U1ser, U2ser, lim, eqs, sol, out},
   computeWU[wsub];
   Module[{c, sol1}, sol1 = First[Solve[Thread[Rest[CoefficientList[Normal[Series[LeqU /. U -> Sum[c[j] Q^j, {j, 0, 4}], {Q, 0, 4}]], Q]] == 0] /. c[0] -> (-a2/a0 /. wsub), Table[c[j], {j, 1, 4}]]];
     U1ser = Sum[c[j] Q^j, {j, 0, 4}] /. c[0] -> (-a2/a0 /. wsub) /. sol1];
   U2ser = Together[SU - U1ser];
   out = Table[Module[{Wser, ff},
      Wser = Association[Table[Wk[k] -> (WU[k] /. U -> If[i == 1, U1ser, U2ser]), {k, 2, 2 ZMAX + 2}]];
      lim = Table[Normal[Series[wickZ[[m + 1]] /. Wser, {Q, 0, 0}]], {m, 0, ZMAX}];
      If[! FreeQ[lim, Q], Print["WARNING: Q->0 limit not finite at root ", i, ": ", lim]];
      ff = Normal[Series[Exp[-Sum[QrrcSym[i, k] z^(2 k + 1) /. wsub, {k, 0, ZMAX}]], {z, 0, ZMAX}]];
      eqs = CoefficientList[Normal[Series[Sum[lim[[m + 1]] z^m, {m, 0, ZMAX}] (1 + Sum[ffc[k] z^k, {k, 1, ZMAX}]), {z, 0, ZMAX}]] - ff, z];
      sol = First[Solve[Thread[eqs == 0], Table[ffc[k], {k, 1, ZMAX}]]];
      Prepend[Table[ffc[k] /. sol, {k, 1, ZMAX}], 1]], {i, 1, 2}];
   out];
(* closed (Stirling) form of the same constants, identical for both roots; verified against ffcLimit to z^6 at general weights *)
ffcClosed[wsub_] := Module[{S, lst},
   S[k_] := (a1^-(2 k - 1) + (a0 + a2)^-(2 k - 1) - e1^-(2 k - 1) - e2^-(2 k - 1)) /. wsub;
   lst = CoefficientList[Normal[Series[Exp[-Sum[BernoulliB[2 k]/(2 k (2 k - 1)) z^(2 k - 1) S[k], {k, 1, ZMAX}]], {z, 0, ZMAX}]], z];
   PadRight[lst, ZMAX + 1]];
(* core: R-matrix at the weights wsub (Q symbolic) with given normalisation lists; sets ratCN[a,i][m], DeltaCN[i], SsumQ, PprodQ, etainvCur *)
rmatrixCore[wsub_, ffc1_, ffc2_] := Module[{SU, PU, LeqU, dUC, etainvN, DeltaLN, m, i, k},
   genSym = U;
   SU = (a0 + Q (a1 - a2))/(Q (a0 + a1)) /. wsub; PU = -a2/(Q (a0 + a1)) /. wsub;
   SsumQ = Together[SU]; PprodQ = Together[PU]; Ssum = SsumQ; Pprod = PprodQ;
   LeqU = Q (a0 + a1) U^2 - (a0 + Q (a1 - a2)) U - a2 /. wsub;
   computeWU[wsub];
   Do[WcN[k] = toC[WU[k]], {k, 2, 2 ZMAX + 2}];
   wsymN = Association[Table[Wk[k] -> WcN[k], {k, 2, 2 ZMAX + 2}]];
   Clear[symPow];
   symPow[sym_, n_, symvals_] := symPow[sym, n, symvals] = cpow[symvals[sym], n];
   Do[wickCN[m] = cEval[wickZ[[m + 1]], wsymN], {m, 0, ZMAX}];
   Do[ratCN[0, 1][m] = cadd @@ Table[cscale[ffc1[[k + 1]], wickCN[m - k]], {k, 0, m}], {m, 0, ZMAX}];
   Do[ratCN[0, 2][m] = cadd @@ Table[cscale[ffc2[[k + 1]], cconj[wickCN[m - k]]], {k, 0, m}], {m, 0, ZMAX}];
   dUC = toC[-D[LeqU, Q]/D[LeqU, U]];
   cDerN[C[a_, b_]] := cadd[C[Together[D[a, Q]], Together[D[b, Q]]], cscale[b, dUC]];
   etainvN = Together[etainvSym /. wsub];
   LCN[1] = toC[LofU /. wsub]; LCN[2] = cconj[LCN[1]];
   DeltaLN = DeltaL /. wsub /. L -> (LofU /. wsub);
   DeltaCN[1] = toC[DeltaLN]; DeltaCN[2] = cconj[DeltaCN[1]];
   Do[dlogDeltaN[i] = cmul[cDerN[DeltaCN[i]], cinv[cscale[2, DeltaCN[i]]]], {i, 1, 2}];
   Do[ratCN[1, i][0] = cmul[LCN[i], ratCN[0, i][0]];
      Do[ratCN[1, i][m] = cadd[cmul[LCN[i], ratCN[0, i][m]], cscale[Q, cadd[cDerN[ratCN[0, i][m - 1]], cscale[-1, cmul[ratCN[0, i][m - 1], dlogDeltaN[i]]]]]], {m, 1, ZMAX}], {i, 1, 2}];
   etainvCur = etainvN;
   "ok"];
(* the three entry points *)
rmatrixAt[wsub_] := Module[{SU = (a0 + Q (a1 - a2))/(Q (a0 + a1)) /. wsub, LeqU = Q (a0 + a1) U^2 - (a0 + Q (a1 - a2)) U - a2 /. wsub, ffs},
   ffs = ffcLimit[wsub, SU, LeqU]; ffcLast = ffs; rmatrixCore[wsub, ffs[[1]], ffs[[2]]]];              (* QRR constants from the Q -> 0 limit *)
rmatrixAtClosed[wsub_] := Module[{ff = ffcClosed[wsub]}, ffcLast = {ff, ff}; rmatrixCore[wsub, ff, ff]];   (* closed-form constants *)
rmatrixAt0[wsub_] := rmatrixCore[Join[{a0 -> 0}, wsub], PadRight[{1}, ZMAX + 1], PadRight[{1}, ZMAX + 1]];  (* alpha0 = 0, alpha_i = eps_i: no constants *)
(* everything at a numeric (or symbolic) Q *)
evalQ[Q0_] := Module[{sub = {Q -> Q0}, i, j, m, n, k, pre},
   Ssum = Together[SsumQ /. sub]; Pprod = Together[PprodQ /. sub];
   Do[Do[rN[a, i][m] = ratCN[a, i][m] /. sub, {m, 0, ZMAX}], {a, 0, 1}, {i, 1, 2}];
   Do[DeltaN[i] = DeltaCN[i] /. sub, {i, 1, 2}];
   Do[DilN[i][n] = cscale[(-1)^n, rN[0, i][n - 1]], {i, 1, 2}, {n, 1, ZMAX + 1}];
   Do[LfN[i][n] = cscale[(-1)^n, rN[1, i][n]], {i, 1, 2}, {n, 0, ZMAX}];
   Do[pre = Table[cadd @@ Flatten[Table[cscale[-(-1)^(p + q) etainvCur[[a + 1, b + 1]], cmul[rN[a, i][p], rN[b, j][q]]], {a, 0, 1}, {b, 0, 1}]], {p, 0, ZMAX}, {q, 0, ZMAX}];
      Do[EdN[i, j][m, n] = cadd @@ Table[cscale[(-1)^k, pre[[m + k + 2, n - k + 1]]], {k, 0, n}], {m, 0, ZMAX - 1}, {n, 0, ZMAX - 1 - m}], {i, 1, 2}, {j, 1, 2}];
   "ok"];
(* the same with every factor expanded as a Laurent series in t = Log Q at Q = 1 to absolute order M.
   The vertex factors Delta_i^(g_v - 1), g_v - 1 = -1 .. gmax - 1, are computed exactly as rational functions first and expanded
   afterwards (series inversion of Delta loses precision and mislabels the valuation of Delta^-1), via deltaPow[i, k]. *)
toS[expr_] := Series[expr /. Q -> Exp[t], {t, 0, M}];
deltaPow[i_, k_] := cpow[DeltaN[i], k];                    (* rational mode: plain power; series mode: redefined by evalQseries *)
evalQseries[Mord_, gmax_: 4] := Module[{i, j, m, n, k, pre, ex},
   M = Mord;
   Ssum = Together[SsumQ]; Pprod = Together[PprodQ];       (* exact powers of Delta with the rational algebra data *)
   Do[ex = cpow[DeltaCN[i], k]; deltaPow[i, k] = C[toS[ex[[1]]], toS[ex[[2]]]], {i, 1, 2}, {k, -1, gmax - 1}];
   Ssum = toS[SsumQ]; Pprod = toS[PprodQ];
   Do[Do[rN[a, i][m] = C[toS[ratCN[a, i][m][[1]]], toS[ratCN[a, i][m][[2]]]], {m, 0, ZMAX}], {a, 0, 1}, {i, 1, 2}];
   Do[DeltaN[i] = C[toS[DeltaCN[i][[1]]], toS[DeltaCN[i][[2]]]], {i, 1, 2}];
   Do[DilN[i][n] = cscale[(-1)^n, rN[0, i][n - 1]], {i, 1, 2}, {n, 1, ZMAX + 1}];
   Do[LfN[i][n] = cscale[(-1)^n, rN[1, i][n]], {i, 1, 2}, {n, 0, ZMAX}];
   Do[pre = Table[cadd @@ Flatten[Table[cscale[-(-1)^(p + q) etainvCur[[a + 1, b + 1]], cmul[rN[a, i][p], rN[b, j][q]]], {a, 0, 1}, {b, 0, 1}]], {p, 0, ZMAX}, {q, 0, ZMAX}];
      Do[EdN[i, j][m, n] = cadd @@ Table[cscale[(-1)^k, pre[[m + k + 2, n - k + 1]]], {k, 0, n}], {m, 0, ZMAX - 1}, {n, 0, ZMAX - 1 - m}], {i, 1, 2}, {j, 1, 2}];
   "ok"];

(* ============================== 4. graph sum, degree-zero term, Q-reconstruction, weight interpolation (verbatim) ============================== *)
graphValueC[gr_] := Module[{V = gr["V"], ed = gr["edges"], H = gr["H"], hsAt = gr["hsAt"], total = C[0, 0], hp, facs},
   Do[Do[hp = ConstantArray[0, Length[H]];
     Do[Do[hp[[hsAt[[v, r]]]] = dec[[v, 1, r]], {r, Length[hsAt[[v]]]}], {v, V}];
     facs = Join[
       Table[deltaPow[c[[v]], gr["g"][[v]] - 1], {v, V}],
       Flatten[Table[Table[LfN[c[[v]]][b], {b, dec[[v, 2]]}], {v, V}]],
       Flatten[Table[Table[DilN[c[[v]]][d], {d, dec[[v, 3]]}], {v, V}]],
       Table[EdN[c[[ed[[e, 1]]]], c[[ed[[e, 2]]]]][hp[[2 e - 1]], hp[[2 e]]], {e, Length[ed]}]];
     total = cadd[total, cscale[Product[dec[[v, 4]], {v, V}], cmul @@ facs]], {dec, gr["decs"]}], {c, Tuples[{1, 2}, V]}];
   cscale[1/gr["aut"], total]];
graphSumC[prepared_] := cadd @@ (graphValueC /@ prepared);
(* degree-zero (constant map) part, symbolic in the five weights *)
const0[gg_] := Module[{ffs, Del0, res = 0, prep},
   Do[ffs[i] = CoefficientList[Normal[Series[Exp[-Sum[QrrcSym[i, k] z^(2 k + 1), {k, 0, ZMAX}]], {z, 0, ZMAX}]], z];
      Del0[i] = Together[DeltaL /. Q -> 0 /. L -> {0, -a0}[[i]]], {i, 1, 2}];
   prep = prepareGraph /@ graphsWithLegs[gg, 0];
   Do[Module[{i = col, tot = 0},
      Do[Dil0[n] = (-1)^n ffs[i][[n]], {n, 1, ZMAX + 1}];
      Do[Ed0[m, n] = Del0[i] Sum[(-1)^k (-(-1)^(m + k + 1 + n - k) ffs[i][[m + k + 2]] ffs[i][[n - k + 1]]), {k, 0, n}], {m, 0, ZMAX - 1}, {n, 0, ZMAX - 1 - m}];
      Do[Module[{V = gr["V"], ed = gr["edges"], H = gr["H"], hsAt = gr["hsAt"], hp},
        Do[hp = ConstantArray[0, Length[H]];
         Do[Do[hp[[hsAt[[v, r]]]] = dec[[v, 1, r]], {r, Length[hsAt[[v]]]}], {v, V}];
         tot += (1/gr["aut"]) Product[dec[[v, 4]], {v, V}] Product[Dil0[d], {v, V}, {d, dec[[v, 3]]}] Product[Ed0[hp[[2 e - 1]], hp[[2 e]]], {e, Length[ed]}] Del0[i]^(Total[gr["g"]] - V),
         {dec, gr["decs"]}]], {gr, prep}];
      res += Together[tot]], {col, 1, 2}];
   Together[res]];
F11const = (-(\[Alpha][0]^3*\[Epsilon][1]) - \[Alpha][0]^2*\[Alpha][1]*\[Epsilon][1] - \[Alpha][0]^2*\[Alpha][2]*\[Epsilon][1] - \[Alpha][0]*\[Alpha][1]*\[Alpha][2]*\[Epsilon][1] - \[Alpha][0]^3*\[Epsilon][2] - \[Alpha][0]^2*\[Alpha][1]*\[Epsilon][2] - \[Alpha][0]^2*\[Alpha][2]*\[Epsilon][2] - \[Alpha][0]*\[Alpha][1]*\[Alpha][2]*\[Epsilon][2] + \[Alpha][0]^2*\[Epsilon][1]*\[Epsilon][2] - \[Alpha][1]*\[Alpha][2]*\[Epsilon][1]*\[Epsilon][2])/(24*(\[Alpha][0] + \[Alpha][1])*(\[Alpha][0] + \[Alpha][2])*\[Epsilon][1]*\[Epsilon][2]);
(* reconstruction of F(Q) from exact samples {q, v}: F = P(Q)/(1-Q)^a, verified on 5 spare points *)
reconQ[pts_, amax_: 16] := Module[{n = Length[pts], qs = pts[[All, 1]], vs = pts[[All, 2]], w, P, ok, res = $Failed},
   If[Union[vs] === {0}, Return[{0, 0}]];
   Do[w = Table[vs[[i]] (1 - qs[[i]])^a, {i, n}];
    P = InterpolatingPolynomial[Transpose[{qs[[1 ;; n - 5]], w[[1 ;; n - 5]]}], Q];
    ok = And @@ Table[(P /. Q -> qs[[i]]) == w[[i]], {i, n - 4, n}];
    If[ok, res = {Expand[P], a}; Return[]], {a, 0, amax}];
   res];
(* interpolation of the instanton part over the weight samples: eps1 eps2 [Q^j] F_g homogeneous of degree 2g, symmetric under a1<->a2, e1<->e2 *)
wvars = {a0, a1, a2, e1, e2};
orbitBasis[deg_] := Module[{mons, orb},
   mons = Select[Tuples[Range[0, deg], 5], Total[#] == deg &];
   orb[m_] := Union[{m, m[[{1, 3, 2, 4, 5}]], m[[{1, 2, 3, 5, 4}]], m[[{1, 3, 2, 5, 4}]]}];
   Union[Table[Total[Times @@@ (wvars^# & /@ orb[m])], {m, mons}]]];
fitPoly[values_Association, deg_] := Module[{basis = orbitBasis[deg], keys = Keys[values], A, rhs, sol, n, res},
   n = Length[basis];
   If[Length[keys] < n + 3, Print["not enough samples: ", Length[keys], " < ", n + 3]; Return[$Failed]];
   A = Table[basis /. allSamples[k][[1]], {k, keys[[1 ;; n]]}]; rhs = Table[values[k], {k, keys[[1 ;; n]]}];
   sol = LinearSolve[A, rhs];
   res = Table[(basis /. allSamples[k][[1]]) . sol - values[k], {k, keys[[n + 1 ;;]]}];
   If[Union[res] =!= {0}, Print["fitPoly: verification FAILED on ", Count[res, Except[0]], " of ", Length[res], " spare samples"]; Return[$Failed]];
   Expand[basis . sol]];
(* assemble the instanton part of genus gi (record index gi in allSamples[k][[2]]) as N(Q)/(eps1 eps2 (1-Q)^a) *)
assemble[gi_, gdeg_] := Module[{a = Max[Table[allSamples[k][[2, gi, 2]], {k, Keys[allSamples]}]], cj, vals, Npoly = 0},
   Do[vals = Association[Table[k -> With[{NN = Expand[allSamples[k][[2, gi, 1]] (1 - Q)^(a - allSamples[k][[2, gi, 2]])], w = allSamples[k][[1]]},
        Coefficient[NN - (NN /. Q -> 0) (1 - Q)^a, Q, j] (e1 e2 /. w)], {k, Keys[allSamples]}]];
      cj = fitPoly[vals, gdeg];
      If[cj === $Failed, Return[$Failed]];
      Npoly += cj Q^j, {j, 1, a}];
   {Npoly/(e1 e2 (1 - Q)^a), a}];

(* ============================== 5. sample lists, orbits, polar parts, CY-torus fit ============================== *)
zmaxFor[gg_] := 3 gg - 3;
okQ[w_] := And @@ (# != 0 & /@ {w[[1]], w[[2]], w[[3]], w[[1]] + w[[2]], w[[1]] + w[[3]], w[[4]], w[[5]], w[[1]] + w[[2]] + w[[3]], w[[4]] + w[[5]], w[[4]] + w[[5]] - w[[1]] - w[[2]] - w[[3]]});
(* integer weight samples for the general torus (same generator and seed as resolvedconifold_genus3.nb: the first 100 are its samples) *)
wsamplesAll := wsamplesAll = Module[{ws = {}}, SeedRandom[20260906]; While[Length[ws] < 400, Module[{w = RandomInteger[{-9, 9}, 5]}, If[okQ[w] && ! MemberQ[ws, w], AppendTo[ws, w]]]]; ws];
nSamplesFor[gg_] := Length[orbitBasis[2 gg]] + 12;                     (* unknowns of the degree-2g symmetric basis plus spare samples *)
qsFor[gg_] := Table[1/(3 + 2 s) + 1/(50 + 3 s), {s, 0, If[gg <= 3, 15, 17]}];
wsubOf[w_] := {a0 -> w[[1]], a1 -> w[[2]], a2 -> w[[3]], e1 -> w[[4]], e2 -> w[[5]]};
(* orbit table (prepared graphs + orbits), stored by rc_prep.m in rc_orbits_g<gg>.mx *)
setupOrbits[gg_] := (prepared = prepareGraph /@ graphsWithLegs[gg, 0]; OT = orbitTable[prepared];
   Print["genus ", gg, ": ", Length[prepared], " graphs, ", Total[Length[#["decs"]] & /@ prepared], " labelled decorations, ", Length[OT["ids"]], " decorated graphs (orbits)"]);
loadOrbits[gg_] := Module[{f = "rc_orbits_g" <> ToString[gg] <> ".mx"},
   If[FileExistsQ[f], Get[f], Print["orbit file ", f, " not found: computing"]; setupOrbits[gg]; DumpSave[f, {prepared, OT}]]];
(* polar part at Q = 1 in t = Log Q of a rational function of Q: coefficients of t^-(2g-2) .. t^-1 *)
polarOfRational[rq_, gg_] := Module[{ser = Series[rq /. Q -> Exp[t], {t, 0, -1}]}, Table[SeriesCoefficient[ser, -n], {n, 2 gg - 2, 1, -1}]];
polarOfSeries[ser_, gg_] := Table[SeriesCoefficient[ser, -n], {n, 2 gg - 2, 1, -1}];
(* the known total on the CY torus: P_g Li_{3-2g}(Q)/(eps1 eps2) (plus a constant) *)
Pg[g_] := Module[{tv}, SeriesCoefficient[Series[e1 e2 tv^2/(4 Sinh[e1 tv/2] Sinh[e2 tv/2]), {tv, 0, 2 g}], 2 g]];
cyPolar[gg_] := polarOfRational[Pg[gg] PolyLog[3 - 2 gg, Q]/(e1 e2), gg];
(* multivariate fit on the CY torus (alpha0 = eps1+eps2-alpha1-alpha2, eps2 = 1): coefficient = N/D with D = eps1^da eps2^da (eps1+eps2)^db,
   N homogeneous of degree d in (alpha1, alpha2, eps1, eps2), symmetric under alpha1<->alpha2 and eps1<->eps2, alpha-degree <= Da *)
symA[j_] := Table[If[a == j - a, a1^a a2^a, a1^a a2^(j - a) + a1^(j - a) a2^a], {a, 0, Floor[j/2]}];
symE[j_] := Table[If[a == j - a, e1^a e2^a, e1^a e2^(j - a) + e1^(j - a) e2^a], {a, 0, Floor[j/2]}];
cyBasis[d_, Da_] := Flatten[Table[Outer[Times, symA[j], symE[d - j]], {j, 0, Min[Da, d]}]];
cySub[pt_] := {a1 -> pt[[1]], a2 -> pt[[2]], e1 -> pt[[3]], e2 -> 1};
(* line probe: along the line pt = base + s dir (integer s) a coefficient is a univariate rational function of s; the generic
   reconstruction gives its numerator/denominator degrees, which fixes d, da, db when D is of the assumed form *)
lineRecon[vals_] := Module[{pts = Transpose[{Range[Length[vals]], vals}]}, reconstructU[pts]];
reconstructU[pts_] := Module[{n = Length[pts], us = pts[[All, 1]], vs = pts[[All, 2]], res = $Failed, T, dd, dn, m, A, rhs, sol, N0, D0, ok},
   If[Union[vs] === {0}, Return[0]];
   Do[Do[dn = T - dd; m = dn + dd + 1; If[m + 4 > n, Continue[]];
     A = Table[Join[Table[us[[s]]^j, {j, 0, dn}], Table[-vs[[s]] us[[s]]^j, {j, 0, dd - 1}]], {s, 1, m}];
     rhs = Table[vs[[s]] us[[s]]^dd, {s, 1, m}];
     sol = Quiet[Check[LinearSolve[A, rhs], $Failed]]; If[sol === $Failed, Continue[]];
     N0 = Sum[sol[[j + 1]] u^j, {j, 0, dn}]; D0 = u^dd + Sum[sol[[dn + 2 + j]] u^j, {j, 0, dd - 1}];
     ok = And @@ Table[((N0 /. u -> us[[s]]) - vs[[s]] (D0 /. u -> us[[s]])) == 0, {s, m + 1, n}];
     If[ok, res = Together[N0/D0]; Return[]], {dd, 0, T}]; If[res =!= $Failed, Return[]], {T, 0, n - 5}];
   res];
(* group fit: vals[key] = list of values at the points pts (random points first, probe points last); returns association key -> rational function *)
fitGroup[vals_Association, pts_, d_, da_, db_, DaMax_, nspare_: 8] := Module[{keys = Keys[vals], Dp, W, res = <||>, still, Da = 0, bs, nU, Ball, lsf, done = False},
   Dp = e1^da e2^da (e1 + e2)^db;
   W = Association[Table[k -> Table[vals[k][[s]] (Dp /. cySub[pts[[s]]]), {s, Length[pts]}], {k, keys}]];
   still = keys;
   While[! done && Da <= DaMax,
    bs = cyBasis[d, Da]; nU = Length[bs];
    If[nU + nspare > Length[pts], Print["fitGroup: not enough points for alpha-degree ", Da, " (", nU, " unknowns, ", Length[pts], " points)"]; Break[]];
    Ball = Table[bs /. cySub[pt], {pt, pts}];
    lsf = Quiet[Check[LinearSolve[Ball[[;; nU]]], $Failed]];
    If[lsf === $Failed, Da++; Continue[]];
    still = Select[still, Module[{sol = lsf[W[#][[;; nU]]]}, If[Ball[[nU + 1 ;;]] . sol == W[#][[nU + 1 ;;]], res[#] = Factor[(sol . bs)/Dp]; False, True]] &];
    If[still === {}, done = True, Da++]];
   {res, still, Da}];
