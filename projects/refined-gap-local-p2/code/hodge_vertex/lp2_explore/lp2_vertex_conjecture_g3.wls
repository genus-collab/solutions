(* lp2_vertex_conjecture_g3.wls -- CHECK 4 at genus 3: Conjecture V (steps (i)-(iii) of Theorem R) for K_P2 x C^2.
   Port of lp2_vertex_conjecture_g2.wls (same frame and conventions) to genus 3, with the valuation readout fixed.

   SETTING (unchanged from the genus-2 script)
   Fibrewise torus T_F, frame B (R' = R_Z h_w^(-1), w = {eps1+eps2, -eps1, -eps2}, mode-B triple Hodge vertex), conifold
   polarisation XX -> XXcf(t_c), singular colour = colour 0 on the REAL branch t = y^(1/3) -> -1/3 (y -> -1/27).
   Exact arithmetic: eps2 = 1, eps1 = x0 rational with sq = Sqrt[x0/(1+x0)^3] rational (x0 = 1/3: sq = 3/8; x0 = 4/5:
   sq = 10/27). Every factor is an exact Laurent series in t_c. Its coefficients lie in Q(omega), omega = om, om^2 + om + 1 = 0,
   for the regular colours 1 and 2 (t = om^i y^(1/3)).
   Normalised factors (one sqrt(Delta) per half-edge or dilaton, as in rc_lib.m):
     NR[a,i][k] = sqrtDelta_i (R' row a at colour i)[z^k];  DelS[i] = Delta_i = sqrtDelta_i^2;  DilN[i][d] = (-1)^d NR[0,i][d-1];
     EN[i,j][m,n] = sqrtDelta_i sqrtDelta_j Ed[i,j][m,n] = Sum_k (-1)^(k+m+n) Sum_ab etainv_ab NR[a,i][m+k+1] NR[b,j][n-k].
   Vertex kernel  Vker[i,h,A] = Delta_i^(h-1) Sum_m 1/m! Sum_{d_1..d_m >= 2} Int_{M_{h,n+m}} Hodge_w psi^A prod psi^d DilN[i][d].
   Conjecture V at genus 3:  V_3(legs, loops) = Sum_A Vker[0,3,A] prod_legs EN[0,j_l][A_l, b_l] prod_loops EN[0,0][A, A'].
     It should be O(t_c^0) for every configuration with n >= 1 legs or loops. The bare vertex should have polar part
     gap_3 t_c^-4, i.e. {t^-4, t^-3, t^-2, t^-1} = {gap_3, 0, 0, 0}, with gap_3 = gapCoefficient[3]
     = 9 * 3! [t^6] t^2/(4 sinh(eps1 t/2) sinh(eps2 t/2)).
   A leg {j, b} is the edge column to colour j (0 = another singular vertex; 1 and 2 = the regular colours) with external
   psi-power b.

   CHANGES FROM THE GENUS-2 SCRIPT
   1. Valuation readout fixed. The genus-2 helper
        valS[s] = s[[4]] + (First[FirstPosition[s[[3]], c_ /; c =!= 0, {Infinity}, {1}]] - 1)
      used FirstPosition with the default Heads -> True. At level {1} the head List (position {0}) matches c =!= 0, so it
      reported nmin - 1 for every genuine SeriesData. That explains Delta_0 at -3 (true -2) and edges and DilN[d >= 2] at -1
      (true 0). DilN[0][1] = -NR[0,0][0] is an exact constant, not a SeriesData, so it read 0 correctly.
      New valS: read the coefficients off Normal[s] and test each one in canonical form (Together, then reduction mod
      om^2+om+1). (B) prints the old and the new readout side by side, and compares the new one with lp2_percolour_g2.log
      (eps1 = 1/3): sqrtDelta_0 -1, unnormalised P0 rows 1, unnormalised Dil 1, E_00 2, E_01 = E_02 = E_10 1, E_11 0.
      Normalised, that means DelS[0] -2, NR[0,0][k] 0, DilN[0][d] 0, EN[0,0] 0, EN[0,1] 0.
   2. Every factor and product goes through canon[] (canonical coefficients, leading zeros trimmed), so that the series
      precision is propagated from the true valuations.
   3. Vertex Hodge integrals: the string and dilaton equations of the normalised Hodge class (it is pulled back under the
      forgetful maps) reduce every integral to psi-lists with entries >= 2. Only those reach hv_lib4.m's Mumford recursion
      (hodgeNormU4, genus <= 4). In genus 3 such a list has at most 6 points. Calibrated in (A0) against the direct
      hodgeNormU4 on sample lists.
   4. Precision: KS = 24 (t and XX series), NT = 28 (tcSetup), KF = 10 (factor truncation). The bare genus-3 vertex is
      DelS^2 (valuation -4) times a regular sum, so it comes out O(t_c^6). A WARNING line is printed if a case comes out
      at O(t_c^k) with k <= 0.
   5. R-matrix range (IMPORTANT). EN[i,j][m,n] needs R' rows to z^(m+n+1), and DilN[i][d] needs them to z^(d-1).
      The script uses the largest cPP_z<K>.m in the data directory (now K = 9). At a genus-h vertex with n = #legs + 2 #loops
      the psi powers reach smax = 3h-3+n, so a case needs
        KZ >= max(smax, smax + 1 + max_l b_l [if legs], smax + 1 [if loops]).
      Genus 3 with K = 9 therefore allows only: the bare vertex (6); one leg with b <= 1 (8+b); two legs with b = (0,0) (9);
      one loop (9). One leg with b >= 2, two legs with some b > 0, three legs (10 + max b), loop + leg (10 + b), two loops
      (11) and two loops + leg (12 + b) are printed as "SKIPPED (needs R' to z^K)". Truncating there would change the
      answer, because the Hodge weights with large psi powers are nonzero. To run every case, compute cPP_z12.m with
      cluster_genus4/genRmat_g4_partA_limit.m 12 (Picard-Fuchs route) and put it in the data directory. The script then
      uses it automatically.

   STEPS
   (A0) hodgeR (string/dilaton-reduced) = direct hodgeNormU4 on sample lists, genus 1-3.
   (A)  C1i-independence of EN (C1p -> XX/(y C1i), C1i = 1 vs 2). Frozen frame-B identity row at colour 0
        (leading NR[0,0][k]/NR[0,0][0] = (-1/E)^k, E = 1+x0). All 42 genus-3 stable graphs, colour-summed and rebuilt
        from the factors, against the stored exact terms (cluster_genus4_hodge/hodge_stablegraphs_genus3.m), t^-4..t^1.
        Their sum against F30 (hodge_F30.m).
   (B)  Factor valuations at colour 0 (fixed readout), with the old readout alongside, the comparison with lp2_percolour,
        and the ranks of the t_c^0 matrices of EN[0,0] and EN[0,1].
   (C)  Conjecture V at genus 3. Cases: bare vertex; one leg (j = 0,1,2; b = 0..4); two legs; three legs; one loop;
        loop + leg; two loops; two loops + leg (heaviest last). Each prints its timing. Weight x0 = 1/3 for everything,
        x0 = 4/5 for the cases with n <= 3.

   USAGE   wolframscript -file lp2_vertex_conjecture_g3.wls [KZ]
           (optional KZ = R-matrix order to use, at most the largest cPP_z<K>.m available; default the largest)
   INPUTS  (read only) cluster_genus4_hodge/: hv_lib4.m, lp2_lib.m, graphs.m, orbits.m, cPP_z9.m,
           hodge_stablegraphs_genus3.m, hodge_F30.m.
   OUTPUT  lp2_vertex_conjecture_g3_out.txt in this directory (lp2_explore), also printed to stdout.

   EXPECTED OUTPUT (if the conifold mechanism transfers)
   (A0) "hodgeR = direct hodgeNormU4 on N sample lists: True".
   (A)  "C1i-independent: True"; "frozen row E/(E+z): True"; for every one of the 42 graphs
        "factors - stored, t^-4..t^1 = {0, 0, 0, 0, 0, 0}"; own polar parts {0,0,0,0} except graph 1 (the gap);
        "sum of the 42 graphs - F30 ... {0, 0, 0, 0, 0, 0}".
   (B)  new valuations: DelS[0] -2, DelS[1] = DelS[2] 0; DilN[0][d] 0; EN[0,0], EN[0,1] all 0 (>= 0);
        "consistent with lp2_percolour_g2.log: True"; old readout = new - 1 for every SeriesData entry;
        rank 1 for the t_c^0 matrices (as in genus 2).
   (C)  bare vertex: {t^-4..t^-1} = {gap_3, 0, 0, 0} OK (x0 = 1/3: gap_3 = -169/27216; x0 = 4/5: gap_3 = -1414951/224000000), and equal to
        the stored polar part of graph 1. Every case with a leg or a loop: {t^-4..t^-1} = {0, 0, 0, 0} OK at both x0.
        With cPP_z9 the SKIPPED lines described in point 5 above appear. The final summary line counts
        OK / FAIL / SKIPPED.
   COST (estimate; not run): setup of the factors a few seconds to ~1 min per weight (three setups); the 42-graph
        calibration well under a minute; each feasible (C) case under a few seconds. Total about 2-6 minutes on one kernel
        with cPP_z9. With cPP_z12 the two-loop cases (n = 4, 5; psi tuples up to 12^5) would add a few minutes. *)

$HistoryLength = 0;
outDir = DirectoryName[ExpandFileName[$InputFileName]];
dataDir = FileNameJoin[{ParentDirectory[ParentDirectory[outDir]], "cluster_genus4_hodge"}];
outFile = FileNameJoin[{outDir, "lp2_vertex_conjecture_g3_out.txt"}];
If[FileExistsQ[outFile], DeleteFile[outFile]];
emit[s_] := Module[{st = OpenAppend[outFile]}, WriteString[st, s, "\n"]; Close[st]; Print[s]];
fmt[e_] := ToString[e, InputForm];
T0 = AbsoluteTime[]; el[] := ToString[Round[AbsoluteTime[] - T0]] <> " s";

SetDirectory[dataDir];
Get["hv_lib4.m"];
availK = Select[ToExpression[StringReplace[FileNameTake /@ FileNames["cPP_z*.m"], {"cPP_z" -> "", ".m" -> ""}]], IntegerQ];
args = Rest[$ScriptCommandLine];
KZ = If[Length[args] >= 1 && IntegerQ[ToExpression[args[[1]]]], Min[ToExpression[args[[1]]], Max[availK]], Max[availK]];
loadRmatrix[KZ];
NT = 28; tcSetup[NT];
KS = 24;   (* precision of the t_c-series of t, y and XX *)
KF = 10;   (* truncation of the normalised factors: O(t_c^KF) *)
HG = 3;    (* genus of the vertex *)
tS = -(1/3) (1 - invcf)^(1/3) + O[tcv]^KS;      (* real cube root of y = (Del - 1)/27 *)
yS = (invcf - 1)/27 + O[tcv]^KS;
XXS = XXcf + O[tcv]^KS;
emit["lp2_vertex_conjecture_g3: " <> $Version <> "; " <> DateString[] <> "; R' rows to z^" <> ToString[KZ] <> " (available cPP_z: " <> fmt[availK] <>
   "); KS = " <> ToString[KS] <> ", KF = " <> ToString[KF] <> "; t(t_c) = " <> fmt[Normal[tS + O[tcv]^3]] <> "; XX(t_c) = " <> fmt[Normal[XXS + O[tcv]^3]]];
rowsSym = rowsTimes[hCoeffs[{e1s + e2s, -e1s, -e2s}, -1, KZ], KZ];
emit["frame-B rows built (" <> el[] <> ")"];

(* ---------- series utilities (canonical coefficients; fixed valuation readout) ---------- *)
redOm[p_] := PolynomialRemainder[Expand[p], om^2 + om + 1, om];
cz[c_] := redOm[Together[c]];
ordS[s_SeriesData] := s[[5]]/s[[6]]; ordS[s_] := Infinity;
chkDen[s_SeriesData] := If[s[[6]] =!= 1, emit["ERROR: fractional powers of t_c in a series"]; Abort[]];
normS[s_SeriesData] := Normal[s] + O[tcv]^(s[[5]]/s[[6]]);
normS[s_] := s;
canon[s_SeriesData] := Module[{p = ordS[s], nrm, lo},
   chkDen[s]; nrm = Expand[Normal[s]];
   If[nrm === 0, Return[SeriesData[tcv, 0, {}, p, p, 1]]];
   lo = Min[Exponent[nrm, tcv, Min], p];
   Sum[cz[Coefficient[nrm, tcv, k]] tcv^k, {k, lo, p - 1}] + O[tcv]^p];
canon[s_] := cz[s];
tidy[s_] := canon[s];
trunc[s_SeriesData] := If[ordS[s] > KF, s + O[tcv]^KF, s];
trunc[s_] := s;
(* valuation = exponent of the first coefficient that is nonzero in canonical form; Infinity = zero to the computed order *)
valS[s_SeriesData] := Module[{nrm, ks},
   chkDen[s]; nrm = Expand[Normal[s]];
   If[nrm === 0, Return[Infinity]];
   ks = Select[Range[Exponent[nrm, tcv, Min], ordS[s] - 1], cz[Coefficient[nrm, tcv, #]] =!= 0 &, 1];
   If[ks === {}, Infinity, First[ks]]];
valS[s_] := Which[cz[s] === 0, Infinity, FreeQ[s, tcv], 0, True, emit["ERROR: valS of a non-series expression in tcv"]; Abort[]];
(* the genus-2 readout, kept verbatim for the diagnostic in (B) only *)
oldValS[s_SeriesData] := If[s[[3]] === {}, Infinity, s[[4]]/s[[6]] + (First[FirstPosition[s[[3]], c_ /; c =!= 0, {Infinity}, {1}]] - 1)/s[[6]]];
oldValS[s_] := If[s === 0, Infinity, 0];
coefS[s_SeriesData, k_] := If[k >= ordS[s], Missing["precision"], cz[Coefficient[Expand[Normal[s]], tcv, k]]];
coefS[s_, k_] := If[k == 0, cz[s], 0];
Clear[powT, powY, powX];
powT[k_] := powT[k] = If[k == 0, 1, tS^k]; powY[k_] := powY[k] = If[k == 0, 1, yS^k]; powX[k_] := powX[k] = If[k == 0, 1, XXS^k];
serPoly[p_, v1_, pw1_] := Module[{rules = CoefficientRules[p, {v1, XX}]},
   If[rules === {}, 0, normS[Total[(#[[2]] pw1[#[[1, 1]]] powX[#[[1, 2]]]) & /@ rules]]]];
serRat[F_, v1_, pw1_] := Module[{f = Together[F], n, d}, n = serPoly[Expand[Numerator[f]], v1, pw1]; d = serPoly[Expand[Denominator[f]], v1, pw1];
   If[Head[d] === SeriesData && Normal[d] === 0, emit["ERROR: denominator vanishes to the working precision"]; Abort[]];
   canon[normS[n/d]]];
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
   Clear[NRS, DelS, DelInv, ENf, Vker, hodgeN, hodgeR, dilProd]; defineMemos[];
   sqv = Sqrt[x0/(1 + x0)^3];
   If[! MatchQ[sqv, _Rational | _Integer], emit["ERROR: sq not rational at x0 = " <> fmt[x0]]; Abort[]];
   sub = {e1s -> x0, e2s -> 1, sq -> sqv};
   etai = etainvF /. sub; Ecur = 1 + x0; xcur = x0;
   usubCur = usubOf[{x0 + 1, -x0, -1}];
   Do[nr = Together[((sqDsym rowsSym[[a + 1, k + 1]]) /. C1p -> XX/(t^3 C1i)) /. C1i -> c1v /. sub];
    If[irrQ[nr], emit["ERROR: normalised row entry not rational in t, XX: a = " <> ToString[a] <> ", k = " <> ToString[k]]; Abort[]];
    Do[NRS[a, i][k] = trunc[colourSeries[nr, i]], {i, 0, 2}], {a, 0, 2}, {k, 0, KZ}];
   delSym = Together[sqDsym^2 /. sub];
   If[irrQ[delSym], emit["ERROR: Delta not rational"]; Abort[]];
   Do[DelS[i] = trunc[colourSeries[delSym, i]]; DelInv[i] = trunc[colourSeries[1/delSym, i]], {i, 0, 2}];
   "ok"];
defineMemos[] := (
ENf[i_, j_][m_, n_] := ENf[i, j][m, n] = If[m + n > KZ - 1, (emit["ERROR: edge index out of range " <> fmt[{m, n}]]; Abort[]),
    trunc[tidy[Sum[(-1)^(k + m + n) Sum[If[etai[[a + 1, b + 1]] === 0, 0, etai[[a + 1, b + 1]] NRS[a, i][m + k + 1] NRS[b, j][n - k]], {a, 0, 2}, {b, 0, 2}], {k, 0, n}]]]];
hodgeN[h_, lst_] := hodgeN[h, lst] = (hodgeNormU4[h, Sort[lst]] /. usubCur);
(* string / dilaton reduction of the normalised Hodge correlator; lst sorted ascending *)
hodgeR[h_, lst_List] := hodgeR[h, lst] = Module[{n = Length[lst], rest},
   Which[
    Total[lst] > 3 h - 3 + n, 0,
    n >= 1 && First[lst] == 0 && 2 h - 2 + n - 1 > 0,
     rest = Rest[lst]; Sum[If[rest[[j]] >= 1, hodgeR[h, Sort[ReplacePart[rest, j -> rest[[j]] - 1]]], 0], {j, Length[rest]}],
    n >= 1 && First[lst] == 1 && 2 h - 2 + n - 1 > 0, (2 h - 2 + n - 1) hodgeR[h, Rest[lst]],
    True, hodgeN[h, lst]]];
dilProd[i_, Dl_List] := dilProd[i, Dl] = If[Dl === {}, 1, trunc[tidy[Times @@ (DilN[i] /@ Dl)]]];
Vker[i_, h_, A_List] := Vker[i, h, A] = Module[{n = Length[A], sA = Total[A], tot = 0, dim, w},
   If[3 h - 3 + n < 0, emit["ERROR: unstable vertex"]; Abort[]];
   Do[dim = 3 h - 3 + n + m;
    Do[w = hodgeR[h, Sort[Join[A, Dl]]]; If[w =!= 0, tot += w ordFactor[Dl]/m! dilProd[i, Dl]], {Dl, msets[m, 2, dim - sA]}],
    {m, 0, 3 h - 3 + n - sA}];
   If[tot === 0, 0, trunc[tidy[Which[h == 0, DelInv[i], h == 1, 1, True, DelS[i]^(h - 1)] tot]]]]);
DilN[i_][d_] := If[1 <= d <= KZ + 1, (-1)^d NRS[0, i][d - 1], (emit["ERROR: dilaton index " <> ToString[d]]; Abort[])];
(* R-matrix order needed by a case at a genus-h vertex (see the header, point 5) *)
reqKZ[h_, legs_List, nloops_] := Module[{n = Length[legs] + 2 nloops, smax},
   smax = 3 h - 3 + n;
   Max[smax, If[legs === {}, 0, smax + 1 + Max[legs[[All, 2]]]], If[nloops > 0, smax + 1, 0]]];
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

(* ============================== (A0) Hodge reduction calibration ============================== *)
xlist = {1/3, 4/5};
x0 = xlist[[1]];
setupFactors[x0, 2];
Module[{samples, res},
  samples = {{1, {1}}, {1, {0, 1}}, {1, {1, 1}}, {1, {0, 0, 3}}, {2, {}}, {2, {0}}, {2, {1}}, {2, {4}}, {2, {0, 5}}, {2, {1, 1, 2}}, {2, {0, 2, 3}}, {2, {0, 0, 6}},
    {3, {}}, {3, {0}}, {3, {1}}, {3, {2}}, {3, {7}}, {3, {0, 7}}, {3, {1, 6}}, {3, {2, 5}}, {3, {1, 1, 2}}, {3, {0, 2, 4}}, {3, {2, 2, 3}}, {3, {0, 0, 9}}, {3, {0, 1, 8}}, {3, {0, 0, 3}}};
  res = Table[Together[hodgeR[s[[1]], Sort[s[[2]]]] - hodgeN[s[[1]], s[[2]]]] === 0, {s, samples}];
  emit["(A0) hodgeR (string/dilaton-reduced) = direct hodgeNormU4 on " <> ToString[Length[samples]] <> " sample lists (genus 1-3): " <> ToString[And @@ res] <>
    If[And @@ res, "", "; failures " <> fmt[Pick[samples, res, False]]] <> " (" <> el[] <> ")"]];

(* ============================== (A) calibration ============================== *)
chkIdx = Select[Tuples[Range[0, 4], 2], Total[#] <= 4 &];
T2 = Table[ENf[0, j][mn[[1]], mn[[2]]], {j, 0, 1}, {mn, chkIdx}];
setupFactors[x0, 1];
T1 = Table[ENf[0, j][mn[[1]], mn[[2]]], {j, 0, 1}, {mn, chkIdx}];
emit["(A) x0 = " <> fmt[x0] <> ": EN[0,0], EN[0,1] (m+n <= 4) C1i-independent (C1i = 1 vs 2): " <> ToString[And @@ Flatten[Map[Normal[tidy[#]] === 0 &, T1 - T2, {2}]]] <> " (" <> el[] <> ")"];
emit["(A) factor precision: min O(t_c^k) over NR[a,i][k], k <= " <> ToString[KZ] <> ": " <> ToString[Min[Table[ordS[NRS[a, i][k]], {a, 0, 2}, {i, 0, 2}, {k, 0, KZ}]]] <>
   "; over DelS, DelInv: " <> ToString[Min[Table[{ordS[DelS[i]], ordS[DelInv[i]]}, {i, 0, 2}]]] <> " (KF = " <> ToString[KF] <> "; lower values reduce the final orders)"];
Module[{v0 = valS[NRS[0, 0][0]], lead},
  lead = Table[coefS[NRS[0, 0][k], v0]/coefS[NRS[0, 0][0], v0], {k, 0, KZ}];
  emit["(A) frozen frame-B identity row at the singular colour (valuation " <> ToString[v0] <> "): leading NR[0,0][k]/NR[0,0][0] = (-1/E)^k, k <= " <> ToString[KZ] <> " (E/(E+z)): " <>
    ToString[Simplify[lead == Table[(-1/Ecur)^k, {k, 0, KZ}]]]]];
stored3 = Get[FileNameJoin[{dataDir, "hodge_stablegraphs_genus3.m"}]];
F30 = Get[FileNameJoin[{dataDir, "hodge_F30.m"}]];
Module[{tot = 0},
  Do[Module[{tg = AbsoluteTime[], mine = graphColSum[gr["genera"], gr["edges"], gr["Aut"]], ref = storedSeries[gr["term"], x0], d},
     tot += mine; d = tidy[mine - ref];
     emit["(A) genus-3 graph " <> ToString[gr["graph"]] <> " genera " <> fmt[gr["genera"]] <> " edges " <> fmt[gr["edges"]] <> ": factors - stored, t^-4..t^1 = " <>
       fmt[laurentList[d, -4, 1]] <> "; own polar {t^-4..t^-1} = " <> fmt[laurentList[mine, -4, -1]] <> "; O(t_c^" <> ToString[ordS[mine]] <> ") (" <>
       ToString[Round[AbsoluteTime[] - tg, 0.1]] <> " s)"]], {gr, stored3}];
  emit["(A) sum of the " <> ToString[Length[stored3]] <> " graphs - F30 (hodge_F30.m), t^-4..t^1: " <> fmt[laurentList[tidy[tot - storedSeries[F30, x0]], -4, 1]] <> " (" <> el[] <> ")"]];

(* ============================== (B) factor valuations at the singular colour ============================== *)
Module[{vDel = Table[valS[DelS[i]], {i, 0, 2}], vDil, vNR, vE00, vE01, vE02, vE11, uP0, uE00, uE01, uE02, uE10, uE11, ok, nmn},
  nmn = Min[4, KZ - 1];
  vDil = Table[valS[DilN[0][d]], {d, 1, KZ + 1}];
  vNR = Table[valS[NRS[a, 0][k]], {a, 0, 2}, {k, 0, KZ}];
  vE00 = Table[valS[ENf[0, 0][m, n]], {m, 0, nmn}, {n, 0, nmn}];
  vE01 = Table[valS[ENf[0, 1][m, n]], {m, 0, nmn}, {n, 0, nmn}];
  vE02 = Table[valS[ENf[0, 2][m, n]], {m, 0, nmn}, {n, 0, nmn}];
  vE11 = Table[valS[ENf[1, 1][m, n]], {m, 0, nmn}, {n, 0, nmn}];
  emit["(B) FIXED readout. valuations: DelS[0..2] " <> fmt[vDel] <> " (sqrtDelta_i: " <> fmt[vDel/2] <> "); DilN[0][d], d = 1.." <> ToString[KZ + 1] <> ": " <> fmt[vDil]];
  emit["(B) NR[a,0][k], k = 0.." <> ToString[KZ] <> " (rows a = 0, 1, 2; rows 1, 2 depend on the C1i gauge, only row 0 is gauge-free): " <> fmt[vNR]];
  emit["(B) EN[0,0][m,n], m,n = 0.." <> ToString[nmn] <> ": " <> fmt[vE00]];
  emit["(B) EN[0,1][m,n]: " <> fmt[vE01] <> "; EN[0,2][m,n]: " <> fmt[vE02] <> "; EN[1,1][m,n]: " <> fmt[vE11]];
  (* unnormalised valuations, to compare with lp2_percolour_g2.log at eps1 = 1/3 (indices <= 3, m+n <= 2) *)
  uP0 = Table[valS[tidy[NRS[0, 0][k]^2 DelInv[0]]]/2, {k, 0, 3}];
  uE00 = Table[valS[tidy[ENf[0, 0][m, n] DelInv[0]]], {m, 0, 2}, {n, 0, 2 - m}];
  uE01 = Table[valS[ENf[0, 1][m, n]] - vDel[[1]]/2 - vDel[[2]]/2, {m, 0, 2}, {n, 0, 2 - m}];
  uE02 = Table[valS[ENf[0, 2][m, n]] - vDel[[1]]/2 - vDel[[3]]/2, {m, 0, 2}, {n, 0, 2 - m}];
  uE10 = Table[valS[ENf[1, 0][m, n]] - vDel[[2]]/2 - vDel[[1]]/2, {m, 0, 2}, {n, 0, 2 - m}];
  uE11 = Table[valS[ENf[1, 1][m, n]] - vDel[[2]], {m, 0, 2}, {n, 0, 2 - m}];
  ok = {vDel[[1]] == -2, vDel[[2]] == 0, vDel[[3]] == 0, uP0 == {1, 1, 1, 1}, Union[Flatten[uE00]] == {2}, Union[Flatten[uE01]] == {1}, Union[Flatten[uE02]] == {1},
    Union[Flatten[uE10]] == {1}, Union[Flatten[uE11]] == {0}};
  emit["(B) unnormalised (as in lp2_percolour): sqrtDelta_0 " <> fmt[vDel[[1]]/2] <> "; P0[k], k <= 3: " <> fmt[uP0] <> " (Dil[n] = (-1)^n P0[n-1]); E_00 " <> fmt[uE00] <>
    "; E_01 " <> fmt[uE01] <> "; E_02 " <> fmt[uE02] <> "; E_10 " <> fmt[uE10] <> "; E_11 " <> fmt[uE11]];
  emit["(B) consistent with lp2_percolour_g2.log (eps1 = 1/3: sqrtDelta_0 -1, sqrtDelta_1,2 0, P0 {1,1,1,1}, Dil {1,1,1,1}, E_00 2, E_01 = E_02 = E_10 1, E_11 0): " <>
    ToString[And @@ ok] <> If[And @@ ok, "", "; individual checks " <> fmt[ok]]];
  emit["(B) old genus-2 readout (FirstPosition, Heads -> True) vs fixed: DelS[0] " <> fmt[{oldValS[DelS[0]], vDel[[1]]}] <> "; DilN[0][1] " <> fmt[{oldValS[DilN[0][1]], valS[DilN[0][1]]}] <>
    "; DilN[0][2] " <> fmt[{oldValS[DilN[0][2]], valS[DilN[0][2]]}] <> "; EN[0,0][0,0] " <> fmt[{oldValS[ENf[0, 0][0, 0]], valS[ENf[0, 0][0, 0]]}] <>
    "; EN[0,1][1,1] " <> fmt[{oldValS[ENf[0, 1][1, 1]], valS[ENf[0, 1][1, 1]]}] <> "; Head of NR[0,0][0]: " <> ToString[Head[NRS[0, 0][0]]]]];
Module[{M00 = Table[coefS[ENf[0, 0][m, n], 0], {m, 0, 3}, {n, 0, 3}], M01 = Table[coefS[ENf[0, 1][m, n], 0], {m, 0, 3}, {n, 0, 3}]},
  emit["(B) t_c^0 matrix of EN[0,0] (m,n <= 3): rank " <> ToString[MatrixRank[M00 /. om -> (-1 + I Sqrt[3])/2]] <> "; first row " <> fmt[M00[[1]]] <> "; first column " <> fmt[M00[[All, 1]]]];
  emit["(B) t_c^0 matrix of EN[0,1] (m,n <= 3): rank " <> ToString[MatrixRank[M01 /. om -> (-1 + I Sqrt[3])/2]] <> "; first row " <> fmt[M01[[1]]]]];
emit["(B) done (" <> el[] <> ")"];

(* ============================== (C) Conjecture V at genus 3 ============================== *)
lg[j_, b_] := {j, b};
legName[legs_] := StringJoin[("(" <> ToString[#[[1]]] <> "," <> ToString[#[[2]]] <> ")") & /@ legs];
caseList = Join[
   {{"bare vertex", {}, 0}},
   Flatten[Table[{"one leg j=" <> ToString[j] <> ", b=" <> ToString[b], {lg[j, b]}, 0}, {j, 0, 2}, {b, 0, 4}], 1],
   Table[{"two legs " <> legName[L], L, 0}, {L, {{lg[0, 0], lg[0, 0]}, {lg[0, 0], lg[1, 0]}, {lg[0, 0], lg[2, 0]}, {lg[1, 0], lg[1, 0]}, {lg[1, 0], lg[2, 0]}, {lg[2, 0], lg[2, 0]},
       {lg[0, 0], lg[0, 1]}, {lg[0, 1], lg[0, 2]}, {lg[1, 1], lg[1, 2]}, {lg[0, 1], lg[1, 0]}}}],
   {{"one loop", {}, 1}},
   Table[{"loop + leg " <> legName[L], L, 1}, {L, {{lg[1, 0]}, {lg[0, 0]}, {lg[0, 1]}}}],
   Table[{"three legs " <> legName[L], L, 0}, {L, {{lg[0, 0], lg[0, 0], lg[0, 0]}, {lg[0, 0], lg[1, 0], lg[2, 0]}, {lg[1, 0], lg[1, 0], lg[2, 0]}, {lg[0, 1], lg[1, 0], lg[2, 0]}}}],
   {{"two loops", {}, 2}, {"two loops + leg (1,0)", {lg[1, 0]}, 2}}];
cheapQ[cs_] := Length[cs[[2]]] + 2 cs[[3]] <= 3;
tally = <|"OK" -> 0, "FAIL" -> 0, "SKIPPED" -> 0, "WARN" -> 0|>;
runCases[x0_, which_] := Module[{gap = Together[gapCoefficient[HG] /. {e1s -> x0, e2s -> 1}], ref1, t0b = {}},
   ref1 = Together[(Select[stored3, #["graph"] == 1 &][[1]]["polar"]) /. x -> x0];
   emit["(C) x0 = " <> fmt[x0] <> ": gap_3 = gapCoefficient[3] = " <> fmt[gap] <> "; stored polar {t^-4..t^-1} of graph 1 = " <> fmt[ref1] <>
     "; stored = {gap_3, 0, 0, 0}: " <> ToString[ref1 === {gap, 0, 0, 0}]];
   Do[Module[{tt0 = AbsoluteTime[], v, l, ok, need = reqKZ[HG, cs[[2]], cs[[3]]], pOK},
      If[need > KZ,
       tally["SKIPPED"] = tally["SKIPPED"] + 1;
       emit["(C) x0 = " <> fmt[x0] <> ", " <> cs[[1]] <> ": SKIPPED (needs R' to z^" <> ToString[need] <> ", have z^" <> ToString[KZ] <> ")"],
       v = ConjV[0, HG, cs[[2]], cs[[3]]]; pOK = precOK[v, 0]; If[! pOK, tally["WARN"] = tally["WARN"] + 1];
       l = laurentList[v, -4, 0];
       ok = valS[v] >= -4 && If[cs[[1]] === "bare vertex", l[[1 ;; 4]] === {gap, 0, 0, 0}, l[[1 ;; 4]] === {0, 0, 0, 0}];
       If[ok, tally["OK"] = tally["OK"] + 1, tally["FAIL"] = tally["FAIL"] + 1];
       If[StringStartsQ[cs[[1]], "one leg"], AppendTo[t0b, {cs[[2, 1, 1]], cs[[2, 1, 2]], l[[5]]}]];
       emit["(C) x0 = " <> fmt[x0] <> ", " <> cs[[1]] <> ": {t^-4, t^-3, t^-2, t^-1} = " <> fmt[l[[1 ;; 4]]] <> If[ok, "  OK", "  FAIL"] <>
         "; t^0 = " <> fmt[l[[5]]] <> "; valuation " <> ToString[valS[v]] <> "; O(t_c^" <> ToString[ordS[v]] <> ")" <>
         If[cs[[1]] === "bare vertex", "; = stored graph-1 polar: " <> ToString[l[[1 ;; 4]] === ref1], ""] <>
         " (" <> ToString[Round[AbsoluteTime[] - tt0, 0.1]] <> " s)"]]],
    {cs, Select[caseList, which]}];
   Do[Module[{vals = SortBy[Select[t0b, #[[1]] == j &], #[[2]] &]},
      If[Length[vals] >= 2 && vals[[1, 2]] == 0 && vals[[1, 3]] =!= 0,
       emit["(C) x0 = " <> fmt[x0] <> ": one leg j=" <> ToString[j] <> ", t^0(b)/t^0(0) for b = " <> fmt[vals[[All, 2]]] <> ": " <>
         fmt[FullSimplify[(#/vals[[1, 3]]) /. om -> (-1 + I Sqrt[3])/2] & /@ vals[[All, 3]]] <> "  (compare (-1/E)^b = " <> fmt[(-1/(1 + x0))^vals[[All, 2]]] <> ")"]]], {j, 0, 2}]];
runCases[xlist[[1]], True &];
emit["(C) done at x0 = " <> fmt[xlist[[1]]] <> " (" <> el[] <> ")"];
Do[setupFactors[xx, 1]; runCases[xx, cheapQ]; emit["(C) done at x0 = " <> fmt[xx] <> " (" <> el[] <> ")"], {xx, Rest[xlist]}];
emit["SUMMARY (C): OK " <> ToString[tally["OK"]] <> ", FAIL " <> ToString[tally["FAIL"]] <> ", SKIPPED " <> ToString[tally["SKIPPED"]] <>
   " (R' to z^" <> ToString[KZ] <> "), precision warnings " <> ToString[tally["WARN"]]];
emit["total " <> el[]];
