(* Appendix A checks in q, with formal log l *)
N0 = 7;
tr[f_] := Expand[Normal[Series[f, {q, 0, N0}]]];
DD[f_] := tr[q D[f, q] + D[f, l]];
C0 = Sum[(5 d)!/(d!)^5 q^d, {d, 0, N0}];
It1 = 5 Sum[(5 d)!/(d!)^5 Sum[1/m, {m, d + 1, 5 d}] q^d, {d, 1, N0}];
I1 = C0 l + It1;
PF4[f_] := Module[{d1, d2, d3, d4}, d1 = DD[f]; d2 = DD[d1]; d3 = DD[d2]; d4 = DD[d3];
   tr[(1 - 3125 q) d4 - 6250 q d3 - 4375 q d2 - 1250 q d1 - 120 q f]];
Print["PF4 C0: ", PF4[C0]];
Print["PF4 I1: ", PF4[I1]];
(* periods from P(w,q) *)
P = Sum[q^d Product[5 w + k, {k, 1, 5 d}]/Product[(w + k)^5 - w^5, {k, 1, d}], {d, 0, N0}];
Pser = Normal[Series[P, {w, 0, 5}]];
pt = Table[Expand[tr[Coefficient[Pser, w, m]]], {m, 0, 5}];
pk = Table[Expand[Sum[l^j/j! pt[[k - j + 1]], {j, 0, k}]], {k, 0, 5}];
Print["p0-C0: ", Expand[pk[[1]] - C0], "  p1-I1: ", Expand[pk[[2]] - I1]];
Print["PF4 p_k k=0..4: ", Table[PF4[pk[[k + 1]]], {k, 0, 4}]];
(* Wronskian *)
Wm = Table[Nest[DD, pk[[k + 1]], j], {j, 0, 3}, {k, 0, 3}];
Wr = tr[Det[Wm]];
L5 = tr[1/(1 - 3125 q)];
Print["Wr - L^10: ", Expand[Wr - tr[L5^2]]];
(* Birkhoff *)
C1 = tr[DD[tr[I1/C0]]];
inv[f_] := tr[1/f];
dT[f_] := tr[inv[C1] DD[f]];
p0 = Table[tr[pk[[k + 1]] inv[C0]], {k, 0, 5}];
k1 = dT[dT[p0[[3]]]];
Print["kappa1 - L^5/(C0^2 C1^3): ", Expand[tr[k1 - L5 inv[C0^2 C1^3]]]];
k2 = dT[tr[dT[dT[p0[[4]]]] inv[k1]]];
Print["kappa2 - 1: ", Expand[k2 - 1]];
k3 = dT[dT[tr[dT[dT[p0[[5]]]] inv[k1]]]];
Print["kappa3 - C0/C1: ", Expand[tr[k3 - C0 inv[C1]]]];
