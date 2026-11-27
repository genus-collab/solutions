(* hvc_vertex_conjecture_g2.wls -- all-order vertex conjecture V at genus 2 (mode B, exact in Q, x symbolic; genus-4 factors, z^9):
   V = Delta_1(t) * Sum_m 1/m! < prod_legs (E_{1j}[a_i, b_i]) prod_loops E_11[a, a'] prod_dilaton Dil_1[d_j] >_{h=2} evaluated at the singular root
   and expanded at Q = 1, t = Log Q. Conjecture: the t^-2 and t^-1 coefficients vanish whenever there is at least one leg or loop; for the bare
   vertex t^-2 = gap_2, t^-1 = 0.  External series on a leg = monomial w^b (b = external psi-power), other end at idempotent j (1 singular, 2 regular). *)
SetDirectory[DirectoryName[$InputFileName]];
Get["hvc_lib.m"]; Get["hvc_factors_g4_Bm.mx"]; tt0 = AbsoluteTime[];
S0 = Together[SsumQ]; P0 = Together[PprodQ]; sq = ((1 + x)/x) Sqrt[1 + 4 x (1/Q - 1)/(1 + x)^2]; U1 = (S0 + sq)/2;
h = 2; c = 1 + x; usub = usubMode["B", x];
msets[k_, lo_, hi_, smax_] := Module[{rec}, rec[kk_, l_, s_] := If[kk == 0, {{}}, Join @@ Table[Prepend[#, a] & /@ rec[kk - 1, a, s - a], {a, l, Min[hi, s]}]]; rec[k, lo, smax]];
ordFactor[lst_] := Length[lst]!/Times @@ (Last /@ Tally[lst]!);
hodgeN[a_List] := hodgeN[a] = hodgeNormU4[h, a] /. usub;
(* legs: list of {j, b} (other end idempotent j, external power b); loops: number of E_11 loops *)
vertexV[legs_List, nloops_] := Module[{tot = C[0, 0], nl = Length[legs], nhalf, dim, fac},
   nhalf = nl + 2 nloops;
   Do[dim = 3 h - 3 + nhalf + m; If[dim < 2 m, Continue[]];
    Do[Module[{apows = Tuples[Range[0, dim - Total[D]], nl], sa},
      Do[sa = Total[A] + Total[D];
       Do[Module[{lp = Join[A, Flatten[L], D], w},
         If[Total[lp] > dim, Continue[]];
         w = hodgeN[lp]; If[w === 0, Continue[]];
         fac = Module[{fs = Join[Table[EdN[1, legs[[i, 1]]][A[[i]], legs[[i, 2]]], {i, nl}], Table[EdN[1, 1][L[[q, 1]], L[[q, 2]]], {q, nloops}], Table[DilN[1][d], {d, D}]]}, If[fs === {}, C[1, 0], cmul @@ fs]];
         tot = cadd[tot, cscale[w ordFactor[D]/m!, fac]]], {L, Tuples[Table[Tuples[Range[0, dim - sa], 2], {nloops}]]}], {A, apows}]], {D, msets[m, 2, dim, dim]}], {m, 0, 3 h - 3 + nhalf}];
   cmul[DeltaCN[1], tot]];
laurent[C[a_, b_]] := Module[{s = Series[(a + b U1) /. Q -> Exp[t], {t, 0, 0}]}, If[Head[s] === SeriesData, Table[Factor[SeriesCoefficient[s, k]], {k, -2, 0}], {0, 0, Factor[Normal[s]]}]];
report[label_, legs_, nloops_] := Module[{t0 = AbsoluteTime[], v, l}, v = vertexV[legs, nloops]; l = laurent[v];
   Print[label, ": {t^-2, t^-1, t^0} = ", l, "   (", Round[AbsoluteTime[] - t0, 0.1], " s)"]];
report["bare vertex (gap expected at t^-2: " <> ToString[InputForm[gapHVC[2, x][[1]]]] <> ")", {}, 0];
(* leg and loop cases run above on 2026-10-30; rerun here only the bare vertex and the two-loop case *)
Print["total ", Round[AbsoluteTime[] - tt0], " s"];
report["two loops (genus-2 vertex with two loops = a genus-4 graph, colouring 1)", {}, 2];
report["two loops + one leg via E12, b = 0", {{2, 0}}, 2];
