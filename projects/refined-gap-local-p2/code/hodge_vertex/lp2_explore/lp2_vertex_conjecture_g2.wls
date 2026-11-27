(* lp2_vertex_conjecture_g2.wls -- CHECK 4: Conjecture V (steps (i)-(iii) of Theorem R) for K_P2 x C^2 at genus 2.
   Port of hodge_vertex_conifold/hvc_vertex_conjecture_g2.wls to local P2: fibrewise torus T_F, frame B (R' = R_Z h_w^(-1),
   w = {eps1+eps2, -eps1, -eps2}, mode-B triple Hodge vertex), conifold polarisation XX -> XXcf(t_c), singular colour = colour 0
   on the REAL branch t = y^(1/3) -> -1/3 (y -> -1/27).  Exact arithmetic: eps2 = 1, eps1 = x0 rational with sq = Sqrt[x0/(1+x0)^3]
   rational (x0 = 1/3: sq = 3/8; x0 = 4/5: sq = 10/27); every factor is an exact Laurent series in t_c with coefficients in Q(omega),
   omega = om (om^2 + om + 1 = 0) for the regular colours 1, 2 (t = om^i y^(1/3)).

   Normalised factors (per half-edge / dilaton one sqrt(Delta), as in rc_lib.m):
     NR[a,i][k] = sqrtDelta_i * (R' row a at colour i)[z^k]  (a = 0, 1, 2 = rows P0, P1, P2 of lp2_lib.m, twisted by h_w^(-1));
     Del[i] = Delta_i = sqrtDelta_i^2;  DilN[i][d] = (-1)^d NR[0,i][d-1];
     EN[i,j][m,n] = sqrtDelta_i sqrtDelta_j Ed[i,j][m,n] = Sum_k (-1)^(k+m+n) Sum_ab etainv_ab NR[a,i][m+k+1] NR[b,j][n-k].
   Vertex kernel  Vker[i,h,A] = Delta_i^(h-1) Sum_m 1/m! Sum_{d_1..d_m >= 2} Int_{M_{h,n+m}} Hodge_w psi^A prod psi^d DilN[i][d].
   Vertex functional (Conjecture V)  V_h(legs, loops) = Sum_A Vker[0,h,A] prod_legs EN[0,j_l][A_l, b_l] prod_loops EN[0,0][A, A'],
   a leg = {j, b}: the edge column to colour j (0 = another singular vertex, 1/2 = regular colours) with external psi-power b.

   Steps.
   (A) calibration: C1i-independence of every EN (C1p -> XX/(y C1i), C1i = 1 vs 2); the frozen frame-B identity row
       (leading t_c coefficients of NR[0,0][k] / NR[0,0][0] = (-1/E)^k, E = 1+x0, i.e. E/(E+z));
       the colour-summed genus-2 stable graphs (all 7) rebuilt from these factors versus the stored exact terms
       cluster_genus4_hodge/hodge_stablegraphs_genus2.m and F20 (hodge_F20.m), through t_c^0.
   (B) valuations of the normalised factors at the singular colour; rank of the t_c^0 matrices of EN[0,0], EN[0,1] (F3 rank-one test).
   (C) Conjecture V at genus 2 (same cases as the conifold script): bare vertex; one leg (j = 0, 1, 2; b = 0..4); two legs;
       one loop; loop + leg; two loops (a genus-4 graph's vertex); two loops + leg.

   EXPECTED OUTPUT (if the conifold mechanism transfers):
     (A) "C1i-independent: True", "frozen row E/(E+z): True", every "graph g: diff t^-2..t^0 = {0,0,0}" and "F20 diff {0,0,0}";
     (B) valuations >= 0 for DilN[0], EN[0,0], EN[0,1]; Del[0] valuation -2; rank 1 for the t_c^0 matrices (if F3 transfers);
     (C) bare vertex {t^-2, t^-1} = {gap_2, 0} with gap_2 = gapCoefficient[2] = stored polar part of graph 1;
         every case with a leg or a loop: {t^-2, t^-1} = {0, 0} at both x0 (printed as "OK").
     Final O-order of every printed series must be >= 1 (checked; a WARNING line otherwise).
   Cost: a few minutes on one kernel (rows to z^9: seconds; series: ~1 min per weight; two-loop cases the largest).
   Inputs (read only): cluster_genus4_hodge/ (hv_lib4.m, lp2_lib.m, graphs.m, orbits.m, cPP_z9.m, hodge_stablegraphs_genus2.m, hodge_F20.m).
   Output: this directory, lp2_vertex_conjecture_g2_out.txt (also printed). Run: wolframscript -file lp2_vertex_conjecture_g2.wls *)

$HistoryLength = 0;
outDir = DirectoryName[ExpandFileName[$InputFileName]];
dataDir = FileNameJoin[{ParentDirectory[ParentDirectory[outDir]], "cluster_genus4_hodge"}];
outFile = FileNameJoin[{outDir, "lp2_vertex_conjecture_g2_out.txt"}];
If[FileExistsQ[outFile], DeleteFile[outFile]];
emit[s_] := Module[{st = OpenAppend[outFile]}, WriteString[st, s, "\n"]; Close[st]; Print[s]];
T0 = AbsoluteTime[]; el[] := ToString[Round[AbsoluteTime[] - T0]] <> " s";

SetDirectory[dataDir];
Get["hv_lib4.m"];
KZ = 9; loadRmatrix[KZ];
NT = 16; tcSetup[NT];
KS = 13;   (* precision of the t_c-series of t and XX *)
KF = 8;    (* truncation of the normalised factors: O(t_c^KF) *)
tS = -(1/3) (1 - invcf)^(1/3) + O[tcv]^KS;      (* real cube root of y = (Del - 1)/27 *)
yS = (invcf - 1)/27 + O[tcv]^KS;
XXS = XXcf + O[tcv]^KS;
emit["lp2_vertex_conjecture_g2: " <> $Version <> "; t(t_c) = " <> ToString[InputForm[Normal[tS + O[tcv]^3]]] <> "; XX(t_c) = " <> ToString[InputForm[Normal[XXS + O[tcv]^3]]]];
rowsSym = rowsTimes[hCoeffs[{e1s + e2s, -e1s, -e2s}, -1, KZ], KZ];
emit["frame-B rows built (" <> el[] <> ")"];

(* ---------- series utilities ---------- *)
redOm[p_] := PolynomialRemainder[Expand[p], om^2 + om + 1, om];
normS[s_SeriesData] := Normal[s] + O[tcv]^(s[[5]]/s[[6]]);
normS[s_] := s;
tidy[s_SeriesData] := normS[SeriesData[s[[1]], s[[2]], redOm /@ s[[3]], s[[4]], s[[5]], s[[6]]]];
tidy[s_] := redOm[s];
trunc[s_SeriesData] := If[s[[5]]/s[[6]] > KF, s + O[tcv]^KF, s];
trunc[s_] := s;
ordS[s_SeriesData] := s[[5]]/s[[6]]; ordS[s_] := Infinity;
valS[s_SeriesData] := If[s[[3]] === {}, Infinity, s[[4]]/s[[6]] + (First[FirstPosition[s[[3]], c_ /; c =!= 0, {Infinity}, {1}]] - 1)/s[[6]]];
valS[s_] := If[s === 0, Infinity, 0];
coefS[s_SeriesData, k_] := If[k >= ordS[s], Missing["precision"], redOm[SeriesCoefficient[s, k]]];
coefS[s_, k_] := If[k == 0, redOm[s], 0];
Clear[powT, powY, powX];
powT[k_] := powT[k] = If[k == 0, 1, tS^k]; powY[k_] := powY[k] = If[k == 0, 1, yS^k]; powX[k_] := powX[k] = If[k == 0, 1, XXS^k];
serPoly[p_, v1_, pw1_] := Module[{rules = CoefficientRules[p, {v1, XX}]},
   If[rules === {}, 0, normS[Total[(#[[2]] pw1[#[[1, 1]]] powX[#[[1, 2]]]) & /@ rules]]]];
serRat[F_, v1_, pw1_] := Module[{f = Together[F], n, d}, n = serPoly[Expand[Numerator[f]], v1, pw1]; d = serPoly[Expand[Denominator[f]], v1, pw1];
   If[Head[d] === SeriesData && Normal[d] === 0, emit["ERROR: denominator vanishes to the working precision"]; Abort[]];
   normS[n/d]];
(* a rational function of t (and XX) at colour i: t = om^i tt *)
colourSeries[F_, 0] := serRat[F /. t -> tt, tt, powT];
colourSeries[F_, i_] := Module[{f = Together[F /. t -> om^i tt], n, d, A, B, n2, d2},
   n = redOm[Numerator[f]]; d = redOm[Denominator[f]];
   A = Coefficient[d, om, 0]; B = Coefficient[d, om, 1];
   n2 = redOm[n ((A - B) - om B)]; d2 = Expand[A^2 - A B + B^2];
   If[! FreeQ[d2, om], emit["ERROR: om in the norm"]; Abort[]];
   tidy[serRat[Coefficient[n2, om, 0]/d2, tt, powT] + om serRat[Coefficient[n2, om, 1]/d2, tt, powT]]];
irrQ[ex_] := ! FreeQ[ex, Power[_Integer | _Rational, _Rational] | \[Epsilon][_] | sq | C1i | C1p | c1pp | LLL[_][_]];

(* ---------- normalised factors at the weight x0 (C1i = c1v) ---------- *)
msets[0, lo_, smax_] := {{}};
msets[k_, lo_, smax_] := Join @@ Table[Prepend[#, a1] & /@ msets[k - 1, a1, smax - a1], {a1, lo, Floor[smax/k]}];
ordFactor[lst_] := Length[lst]!/Times @@ ((Last /@ Tally[lst])!);
setupFactors[x0_, c1v_] := Module[{sqv, sub, nr, delSym},
   Clear[NRS, Del, DelInv, ENf, Vker, hodgeN]; defineMemos[];
   sqv = Sqrt[x0/(1 + x0)^3];
   If[! MatchQ[sqv, _Rational | _Integer], emit["ERROR: sq not rational at x0 = " <> ToString[x0]]; Abort[]];
   sub = {e1s -> x0, e2s -> 1, sq -> sqv};
   etai = etainvF /. sub; Ecur = 1 + x0; xcur = x0;
   usubCur = usubOf[{x0 + 1, -x0, -1}];
   Do[nr = Together[((sqDsym rowsSym[[a + 1, k + 1]]) /. C1p -> XX/(t^3 C1i)) /. C1i -> c1v /. sub];
    If[irrQ[nr], emit["ERROR: normalised row entry not rational in t, XX: a = " <> ToString[a] <> ", k = " <> ToString[k]]; Abort[]];
    Do[NRS[a, i][k] = trunc[colourSeries[nr, i]], {i, 0, 2}], {a, 0, 2}, {k, 0, KZ}];
   delSym = Together[sqDsym^2 /. sub];
   If[irrQ[delSym], emit["ERROR: Delta not rational"]; Abort[]];
   Do[Del[i] = trunc[colourSeries[delSym, i]]; DelInv[i] = trunc[colourSeries[1/delSym, i]], {i, 0, 2}];
   "ok"];
defineMemos[] := (
ENf[i_, j_][m_, n_] := ENf[i, j][m, n] = If[m + n > KZ - 1, (emit["ERROR: edge index out of range " <> ToString[{m, n}]]; Abort[]),
    tidy[Sum[(-1)^(k + m + n) Sum[If[etai[[a + 1, b + 1]] === 0, 0, etai[[a + 1, b + 1]] NRS[a, i][m + k + 1] NRS[b, j][n - k]], {a, 0, 2}, {b, 0, 2}], {k, 0, n}]]];
hodgeN[h_, lst_] := hodgeN[h, lst] = (hodgeNormU4[h, Sort[lst]] /. usubCur);
Vker[i_, h_, A_List] := Vker[i, h, A] = Module[{n = Length[A], sA = Total[A], tot = 0, dim, w},
   If[3 h - 3 + n < 0, emit["ERROR: unstable vertex"]; Abort[]];
   Do[dim = 3 h - 3 + n + m;
    Do[w = hodgeN[h, Join[A, Dl]]; If[w =!= 0, tot += w ordFactor[Dl]/m! (Times @@ (DilN[i] /@ Dl))], {Dl, msets[m, 2, dim - sA]}],
    {m, 0, 3 h - 3 + n - sA}];
   If[tot === 0, 0, tidy[Which[h == 0, DelInv[i], h == 1, 1, True, Del[i]^(h - 1)] tot]]]);
DilN[i_][d_] := If[1 <= d <= KZ + 1, (-1)^d NRS[0, i][d - 1], (emit["ERROR: dilaton index " <> ToString[d]]; Abort[])];
ConjV[i_, h_, legs_List, nloops_] := Module[{nl = Length[legs], n, smax, tot = 0, vk},
   n = nl + 2 nloops; smax = 3 h - 3 + n;
   Do[If[Total[A] <= smax, vk = Vker[i, h, Sort[A]];
     If[vk =!= 0, tot += vk Product[ENf[i, legs[[l, 1]]][A[[l]], legs[[l, 2]]], {l, nl}] Product[ENf[i, i][A[[nl + 2 q - 1]], A[[nl + 2 q]]], {q, nloops}]]],
    {A, Tuples[Range[0, smax], n]}];
   tidy[tot]];
(* colour-summed stable graph from the factors: genera, edge list {v1, v2} (loops {v, v}), |Aut| *)
graphColSum[gen_, edges_, aut_] := Module[{V = Length[gen], ne = Length[edges], hv, hsAt, perV, tot = 0, p},
   hv = Flatten[edges]; hsAt = Table[Flatten[Position[hv, v]], {v, V}];
   perV = Table[Select[Tuples[Range[0, 3 gen[[v]] - 3 + Length[hsAt[[v]]]], Length[hsAt[[v]]]], Total[#] <= 3 gen[[v]] - 3 + Length[hsAt[[v]]] &], {v, V}];
   Do[Do[p = ConstantArray[0, 2 ne];
      Do[If[hsAt[[v]] =!= {}, p[[hsAt[[v]]]] = pv[[v]]], {v, V}];
      tot += Product[Vker[c[[v]], gen[[v]], Sort[pv[[v]]]], {v, V}] Product[ENf[c[[edges[[e, 1]]]], c[[edges[[e, 2]]]]][p[[2 e - 1]], p[[2 e]]], {e, ne}],
      {pv, Tuples[perV]}], {c, Tuples[{0, 1, 2}, V]}];
   tidy[tot/aut]];
storedSeries[term_, x0_] := serRat[term /. {e1s -> x0, e2s -> 1}, y, powY];
laurentList[s_, kmin_, kmax_] := Table[coefS[s, k], {k, kmin, kmax}];
precOK[s_, k_] := If[ordS[s] <= k, emit["WARNING: precision O(t_c^" <> ToString[ordS[s]] <> ") insufficient"]; False, True];

(* ============================== (A) calibration ============================== *)
xlist = {1/3, 4/5};
x0 = xlist[[1]];
setupFactors[x0, 2];
chkIdx = Select[Tuples[Range[0, 4], 2], Total[#] <= 4 &];
T2 = Table[ENf[0, j][mn[[1]], mn[[2]]], {j, 0, 1}, {mn, chkIdx}];
setupFactors[x0, 1];
T1 = Table[ENf[0, j][mn[[1]], mn[[2]]], {j, 0, 1}, {mn, chkIdx}];
emit["(A) x0 = " <> ToString[x0] <> ": EN[0,0], EN[0,1] (m+n <= 4) C1i-independent (C1i = 1 vs 2): " <> ToString[And @@ Flatten[Map[Normal[tidy[#]] === 0 &, T1 - T2, {2}]]] <> " (" <> el[] <> ")"];
Module[{v0 = valS[NRS[0, 0][0]], lead},
  lead = Table[coefS[NRS[0, 0][k], v0]/coefS[NRS[0, 0][0], v0], {k, 0, KZ}];
  emit["(A) frozen frame-B identity row at the singular colour (valuation " <> ToString[v0] <> "): leading NR[0,0][k]/NR[0,0][0] = (-1/E)^k, k <= 9 (E/(E+z)): " <>
    ToString[Simplify[lead == Table[(-1/Ecur)^k, {k, 0, KZ}]]]]];
stored2 = Get[FileNameJoin[{dataDir, "hodge_stablegraphs_genus2.m"}]];
F20 = Get[FileNameJoin[{dataDir, "hodge_F20.m"}]];
Module[{tot = 0},
  Do[Module[{mine = graphColSum[gr["genera"], gr["edges"], gr["Aut"]], ref = storedSeries[gr["term"], x0], d},
     tot += mine; d = tidy[mine - ref];
     emit["(A) genus-2 graph " <> ToString[gr["graph"]] <> " genera " <> ToString[gr["genera"]] <> " edges " <> ToString[gr["edges"]] <> ": factors - stored, t^-2..t^1 = " <>
       ToString[InputForm[laurentList[d, -2, 1]]] <> "; own polar {t^-2, t^-1} = " <> ToString[InputForm[laurentList[mine, -2, -1]]] <> " (" <> el[] <> ")"]], {gr, stored2}];
  emit["(A) sum of the 7 graphs - F20 (hodge_F20.m), t^-2..t^1: " <> ToString[InputForm[laurentList[tidy[tot - storedSeries[F20, x0]], -2, 1]]]]];

(* ============================== (B) factor valuations at the singular colour ============================== *)
emit["(B) valuations: Del[0] " <> ToString[valS[Del[0]]] <> "; DilN[0][d], d = 1..10: " <> ToString[Table[valS[DilN[0][d]], {d, 1, KZ + 1}]]];
emit["(B) valuations EN[0,0][m,n], m,n = 0..4: " <> ToString[Table[valS[ENf[0, 0][m, n]], {m, 0, 4}, {n, 0, 4}]]];
emit["(B) valuations EN[0,1][m,n], m,n = 0..4: " <> ToString[Table[valS[ENf[0, 1][m, n]], {m, 0, 4}, {n, 0, 4}]]];
Module[{M00 = Table[coefS[ENf[0, 0][m, n], 0], {m, 0, 3}, {n, 0, 3}], M01 = Table[coefS[ENf[0, 1][m, n], 0], {m, 0, 3}, {n, 0, 3}]},
  emit["(B) t_c^0 matrix of EN[0,0] (m,n <= 3): rank " <> ToString[MatrixRank[M00 /. om -> (-1 + I Sqrt[3])/2]] <> "; first row " <> ToString[InputForm[M00[[1]]]] <> "; first column " <> ToString[InputForm[M00[[All, 1]]]]];
  emit["(B) t_c^0 matrix of EN[0,1] (m,n <= 3): rank " <> ToString[MatrixRank[M01 /. om -> (-1 + I Sqrt[3])/2]] <> "; first row " <> ToString[InputForm[M01[[1]]]]]];

(* ============================== (C) Conjecture V at genus 2 ============================== *)
caseList = Join[
   {{"bare vertex", {}, 0}},
   Flatten[Table[{"one leg j=" <> ToString[j] <> ", b=" <> ToString[b], {{j, b}}, 0}, {j, 0, 2}, {b, 0, 4}], 1],
   {{"two legs (0,0)(0,0)", {{0, 0}, {0, 0}}, 0}, {"two legs (0,0)(0,1)", {{0, 0}, {0, 1}}, 0}, {"two legs (0,1)(0,2)", {{0, 1}, {0, 2}}, 0}, {"two legs (0,2)(0,3)", {{0, 2}, {0, 3}}, 0},
    {"two legs (1,0)(1,0)", {{1, 0}, {1, 0}}, 0}, {"two legs (1,1)(1,2)", {{1, 1}, {1, 2}}, 0}, {"two legs (1,0)(2,0)", {{1, 0}, {2, 0}}, 0}, {"two legs (0,1)(1,0)", {{0, 1}, {1, 0}}, 0},
    {"one loop", {}, 1}, {"loop + leg j=1, b=0", {{1, 0}}, 1}, {"loop + leg j=0, b=1", {{0, 1}}, 1}},
   {{"two loops", {}, 2}, {"two loops + leg j=1, b=0", {{1, 0}}, 2}}];
heavy = {"two loops", "two loops + leg j=1, b=0"};
runCases[x0_, which_] := Module[{gap = gapCoefficient[2] /. {e1s -> x0, e2s -> 1}, ref1 = (Select[stored2, #["graph"] == 1 &][[1]]["polar"]) /. x -> x0, t0b = {}},
   emit["(C) x0 = " <> ToString[x0] <> ": gap_2 = gapCoefficient[2] = " <> ToString[InputForm[gap]] <> "; stored polar of graph 1 = " <> ToString[InputForm[Together[ref1]]]];
   Do[Module[{tt0 = AbsoluteTime[], v, l, ok},
      v = ConjV[0, 2, cs[[2]], cs[[3]]]; precOK[v, 0];
      l = laurentList[v, -2, 0];
      ok = If[cs[[1]] === "bare vertex", l[[1]] === Together[gap] && l[[2]] === 0, l[[1]] === 0 && l[[2]] === 0];
      If[StringStartsQ[cs[[1]], "one leg"], AppendTo[t0b, {cs[[2, 1]], l[[3]]}]];
      emit["(C) x0 = " <> ToString[x0] <> ", " <> cs[[1]] <> ": {t^-2, t^-1} = " <> ToString[InputForm[l[[1 ;; 2]]]] <> If[ok, "  OK", "  FAIL"] <>
        "; t^0 = " <> ToString[InputForm[l[[3]]]] <> "; O(t_c^" <> ToString[ordS[v]] <> ") (" <> ToString[Round[AbsoluteTime[] - tt0, 0.1]] <> " s)"]],
    {cs, Select[caseList, which]}];
   Do[Module[{vals = Select[t0b, #[[1]] == j &][[All, 2]]},
      If[Length[vals] == 5 && vals[[1]] =!= 0, emit["(C) x0 = " <> ToString[x0] <> ": one leg j=" <> ToString[j] <> ", t^0(b)/t^0(0) for b = 0..4: " <> ToString[InputForm[FullSimplify[(#/vals[[1]]) /. om -> (-1 + I Sqrt[3])/2] & /@ vals]] <>
          "  (compare (-1/E)^b = " <> ToString[InputForm[Table[(-1/(1 + x0))^b, {b, 0, 4}]]] <> ")"]]], {j, 0, 2}]];
runCases[xlist[[1]], True &];
emit["(C) done at x0 = " <> ToString[xlist[[1]]] <> " (" <> el[] <> ")"];
Do[setupFactors[xx, 1]; runCases[xx, ! MemberQ[heavy, #[[1]]] &]; emit["(C) done at x0 = " <> ToString[xx] <> " (" <> el[] <> ")"], {xx, Rest[xlist]}];
emit["total " <> el[]];
