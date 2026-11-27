(* sanity checks for appendix A: periods, Birkhoff recursion, kappa^LR *)
NQ = 6;
P = Sum[q^d Product[5 w + k, {k, 1, 5 d}]/Product[(w + k)^5 - w^5, {k, 1, d}], {d, 0, NQ}];
pt[m_] := pt[m] = Normal[Series[SeriesCoefficient[P, {w, 0, m}], {q, 0, NQ}]];
C0 = Sum[(5 d)!/(d!)^5 q^d, {d, 0, NQ}];
I1t = 5 Sum[(5 d)!/(d!)^5 Sum[1/m, {m, d + 1, 5 d}] q^d, {d, 1, NQ}];
Print["pt0-C0: ", Expand[pt[0] - C0], "  pt1-I1t: ", Expand[pt[1] - I1t]];
(* p_k with formal log l ; Dq = q d/dq + d/dl *)
p[k_] := Sum[l^j/j! pt[k - j], {j, 0, k}];
Dq[f_] := q D[f, q] + D[f, l];
Dn[f_, 0] := f; Dn[f_, n_] := Dq[Dn[f, n - 1]];
PF4[f_] := Dn[f, 4] - 5 q (5 Dq[#] + 1 # &)[(5 Dq[#] + 2 # &)[(5 Dq[#] + 3 # &)[(5 Dq[#] + 4 # &)[f]]]];
tr[f_] := Normal[Series[f, {q, 0, NQ}]];
Do[Print["PF4 p_", k, " = ", Expand[tr[PF4[p[k]]]]], {k, 0, 5}];
(* explicit form of PF4 *)
PF4b[f_] := (1 - 5^5 q) Dn[f, 4] - 6250 q Dn[f, 3] - 4375 q Dn[f, 2] - 1250 q Dn[f, 1] - 120 q f;
Print["PF4 forms agree: ", Expand[PF4[G[q, l]] - PF4b[G[q, l]]] ];
(* Wronskian *)
L = (1 - 5^5 q)^(-1/5);
Wr = Det[Table[Dn[p[k], j], {j, 0, 3}, {k, 0, 3}]];
Print["Wr - L^10 = ", Simplify[tr[Wr - L^10]]];
(* mirror map *)
T = l + I1t/C0; C1 = Dq[T];
dT[f_] := Dq[f]/C1;
p0[k_] := p[k]/C0;
k1 = dT[dT[p0[2]]];
k2 = dT[dT[dT[p0[3]]]/k1];
k3 = dT[dT[dT[dT[p0[4]]]/k1]];
Print["p0_1 - T: ", Simplify[tr[p0[1] - T]]];
Print["kappa1 - L^5/(C0^2 C1^3): ", Simplify[tr[k1 - L^5/(C0^2 C1^3)]]];
Print["kappa2 - 1: ", Simplify[tr[k2 - 1]]];
Print["kappa3 - C0/C1: ", Simplify[tr[k3 - C0/C1]]];
Print["kappa1 = ", Simplify[tr[k1]]];
(* kappa4 from J_4: p^{(4)}_5 derivative should be kappa4 = C0/C1 as well *)
p1[k_] := dT[p0[k]]; p2[k_] := dT[p1[k]]/k1; p3[k_] := dT[p2[k]]; p4[k_] := dT[p3[k]]/k3;
Print["kappa4 - C0/C1 (from dT p4_5): ", Simplify[tr[dT[p4[5]] - C0/C1]]];
