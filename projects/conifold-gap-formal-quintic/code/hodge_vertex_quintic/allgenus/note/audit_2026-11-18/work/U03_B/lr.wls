(* Large-radius sanity checks, appendix A (a1_cyclicLR) *)
M = 7; (* q-order *)
ser[e_] := Normal[Series[e, {q, 0, M}]];
Dq[e_] := Expand[q D[e, q] + D[e, lq]];  (* lq = formal log q, Dq lq = 1 *)
Pw = Sum[q^d Product[5 w + k, {k, 1, 5 d}]/Product[(w + k)^5 - w^5, {k, 1, d}], {d, 0, M}];
pt[m_] := pt[m] = Expand[SeriesCoefficient[Pw, {w, 0, m}]];
p[k_] := p[k] = Expand[Sum[lq^j/j! pt[k - j], {j, 0, k}]];
C0 = pt[0]; 
I1 = p[1];
Print["C0 check: ", Expand[C0 - Sum[(5 d)!/(d!)^5 q^d, {d, 0, M}]]];
It1 = 5 Sum[(5 d)!/(d!)^5 Sum[1/m, {m, d + 1, 5 d}] q^d, {d, 1, M}];
Print["I1 check: ", Expand[pt[1] - It1]];
PF4[f_] := Expand[Dq[Dq[Dq[Dq[f]]]] - 5 q Expand[(5 Dq[#] + 1 #) &[(5 Dq[#] + 2 #) &[(5 Dq[#] + 3 #) &[(5 Dq[f] + 4 f)]]]]];
trunc[e_] := Expand[Normal[Series[e, {q, 0, M - 1}]]];
Do[Print["PF4 p_", k, " = ", trunc[PF4[p[k]]]], {k, 0, 4}];
(* explicit form *)
PF4b[f_] := Expand[(1 - 5^5 q) Dq[Dq[Dq[Dq[f]]]] - 6250 q Dq[Dq[Dq[f]]] - 4375 q Dq[Dq[f]] - 1250 q Dq[f] - 120 q f];
Print["PF4 forms agree on p4: ", trunc[PF4[p[4]] - PF4b[p[4]]]];
T = I1/C0; C1 = ser[Dq[ser[T]] ];
C1 = ser[(Dq[I1] C0 - I1 Dq[C0])/C0^2];
Print["C1 has no lq: ", FreeQ[Expand[C1], lq], "  C1 = ", Series[C1, {q, 0, 2}]];
dT[f_] := ser[Dq[f]/C1];
L = (1 - 5^5 q)^(-1/5);
p0k[k_] := ser[p[k]/C0];
kap1 = dT[dT[p0k[2]]];
Print["kappa1 - L^5/(C0^2 C1^3) = ", trunc[kap1 - L^5/(C0^2 C1^3)]];
kap2 = dT[ser[dT[dT[p0k[3]]]/kap1]];
Print["kappa2 - 1 = ", trunc[kap2 - 1]];
kap3 = dT[dT[ser[dT[dT[p0k[4]]]/kap1]]];
Print["kappa3 - C0/C1 = ", trunc[kap3 - C0/C1]];
Print["p0_1 - T = ", trunc[p0k[1] - T]];
(* Wronskian *)
W = Det[Table[Nest[Dq, p[k], j], {j, 0, 3}, {k, 0, 3}]];
Print["Wr - L^10 = ", trunc[W - L^10]];
(* DL *)
Print["Dq L - (L^6-L)/5 = ", trunc[q D[L, q] - (L^6 - L)/5]];
