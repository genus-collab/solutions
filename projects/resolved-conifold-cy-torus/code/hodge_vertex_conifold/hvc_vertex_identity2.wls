(* hvc_vertex_identity2.wls -- (1) exact-in-Q structure tests: is the identity row of R' at the singular idempotent geometric in z for all Q,
   are the edge terms of rank one for all Q;  (2) the leading-order vertex functional Z[h, n] for genus 3 (and 4 if feasible) by multiset enumeration. *)
SetDirectory[DirectoryName[$InputFileName]];
Get["hvc_lib.m"]; Get["hvc_factors_g3_Bm.mx"];
S0 = Together[SsumQ]; P0 = Together[PprodQ]; sq = ((1 + x)/x) Sqrt[1 + 4 x (1/Q - 1)/(1 + x)^2]; U1 = (S0 + sq)/2; U2 = (S0 - sq)/2;
c = 1 + x;
(* exact-in-Q tests are done in the quadratic algebra: an algebra element C[a,b] is zero iff a = b = 0 *)
zeroC[C[a_, b_]] := Together[a] === 0 && Together[b] === 0;
r[m_] := ratCN[0, 1][m];
Print["identity row geometric for all Q (r[m] r[1] = r[m+1], m = 1..5): ", Table[zeroC[cadd[cmul[r[m], r[1]], cscale[-1, r[m + 1]]]], {m, 1, 5}]];
Print["  r[2] - r[1]^2 series in t at the singular root: ", Series[(cadd[r[2], cscale[-1, cmul[r[1], r[1]]]] /. C[a_, b_] :> a + b U1) /. Q -> Exp[t], {t, 0, 2}] // Factor];
rk1[i_, j_] := Table[zeroC[cadd[cmul[EdN[i, j][m, n], EdN[i, j][0, 0]], cscale[-1, cmul[EdN[i, j][m, 0], EdN[i, j][0, n]]]]], {m, 1, 3}, {n, 1, 3}];
Print["edge terms rank one for all Q: E11 ", rk1[1, 1], ", E12 ", rk1[1, 2], ", E22 ", rk1[2, 2]];
Print["  E11[1,1]E11[0,0] - E11[1,0]E11[0,1] series in t at the singular root: ", Series[(cadd[cmul[EdN[1, 1][1, 1], EdN[1, 1][0, 0]], cscale[-1, cmul[EdN[1, 1][1, 0], EdN[1, 1][0, 1]]]] /. C[a_, b_] :> a + b U1) /. Q -> Exp[t], {t, 0, 2}] // Factor];
(* vertex functional with multiset enumeration *)
usub = usubOf[{c, -x, -1}];
phiC[k_] := 1/c^(k + 1); TC[d_] := -1/c^(d - 1);
multisets[k_, lo_, hi_, smax_] := If[k == 0, {{}}, Join @@ Table[Join @@ Table[Prepend[#, a] & /@ multisets[k - 1, a, hi, smax - a], {a, lo, Min[hi, smax]}], {1}]];
(* multisets as sorted lists a1 <= a2 <= ... with sum <= smax *)
msets[k_, lo_, hi_, smax_] := Module[{rec}, rec[kk_, l_, s_] := If[kk == 0, {{}}, Join @@ Table[Prepend[#, a] & /@ rec[kk - 1, a, s - a], {a, l, Min[hi, s]}]]; rec[k, lo, smax]];
ordFactor[lst_] := Length[lst]!/Times @@ (Last /@ Tally[lst]!);
Zm[h_, n_] := Zm[h, n] = Module[{tot = 0, dim},
   Do[dim = 3 h - 3 + n + m; If[dim < 2 m, Continue[]];
    Do[Module[{sa = Total[A]},
      Do[tot += ordFactor[A] (ordFactor[D]/m!) (Times @@ (phiC /@ A)) (Times @@ (TC /@ D)) (hodgeNormU4[h, Join[A, D]] /. usub), {D, msets[m, 2, dim - sa, dim - sa]}]], {A, msets[n, 0, dim, dim]}], {m, 0, 3 h - 3 + n}];
   Together[tot]];
Print["consistency with the brute-force values: Z[2,0] gap: ", Together[(-x c^3) Zm[2, 0] - gapHVC[2, x][[1]]] === 0, ", Z[2,1] = ", Zm[2, 1], ", Z[2,2] = ", Zm[2, 2], ", Z[1,1] = ", Factor[Zm[1, 1]]];
Print["Z[3, 0] reproduces the gap: ", Together[(-x c^3)^2 Zm[3, 0] - gapHVC[3, x][[1]]] === 0];
Do[Print["Z[3, ", n, "] = ", Factor[Zm[3, n]]], {n, 1, 3}];
Print["genus 4: gap check (-x c^3)^3 Z[4,0] - gap_4 = ", Together[(-x c^3)^3 Zm[4, 0] - gapHVC[4, x][[1]]]];
Print["Z[4, 1] = ", Factor[Zm[4, 1]]];
Print["Z[4, 2] = ", Factor[Zm[4, 2]]];
