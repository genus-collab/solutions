(* Check of the edge identity (2-33) off the diagonal, kappa = 1, v = 1:
   (psi1+psi2) L_{c_i} x L_{c_j}[B] = - (Rinv(psi1) Rinv(psi2)^T)_{ij},  Rinv(psi) = M0(-psi)^T (row = Laplace point).
   Also tests the transposed reading. *)
Get["m0.mx"];
Nd = 3; prec = 50; kappa = 1;
sqk = Sqrt[-5 kappa];
x[p_] := kappa (p^4/4 + 1/p);
cs = Table[N[z5^(-j), prec], {j, 0, 4}];
pser[c_, w_] := Module[{a, pw, eq, s},
  pw = c + w/(sqk c) + Sum[a[k] w^k, {k, 2, 2 Nd + 3}];
  s = {};
  Do[eq = Chop[SeriesCoefficient[Series[x[pw //. s] - x[c] + w^2/2, {w, 0, k + 1}], k + 1], 10^-40];
     s = Join[s, {a[k] -> (a[k] /. Solve[eq == 0, a[k]][[1]])}], {k, 2, 2 Nd + 3}];
  Chop[(pw //. s) /. a[_] -> 0, 10^-40]];
G[m_, th_] := Which[OddQ[m], 0, m >= 0, (m - 1)!! th^(-m/2), True, (-1)^(-m/2) th^(-m/2)/((-m - 1)!!)];
Rinv[psi_] := Transpose[M0 /. t -> -psi/kappa];
res = {};
Do[If[i != j,
   Module[{p1 = pser[cs[[i]], w1], p2 = pser[cs[[j]], w2], Bf, ser, LL, rhs, rhsT, d1, d2},
    Bf = D[p1, w1] D[p2, w2]/(p1 - p2)^2;
    ser = Normal[Series[Bf, {w1, 0, 2 Nd}, {w2, 0, 2 Nd}]];
    LL = Sum[Coefficient[Coefficient[ser, w1, m], w2, n] G[m, 1/s1] G[n, 1/s2], {m, 0, 2 Nd}, {n, 0, 2 Nd}];
    rhs = -(Rinv[s1].Transpose[Rinv[s2]])[[i, j]];
    rhsT = -(Transpose[Rinv[s1]].Rinv[s2])[[i, j]];
    d1 = Expand[(s1 + s2) LL - rhs];
    d2 = Expand[(s1 + s2) LL - rhsT];
    (* keep total degree <= Nd *)
    trn[e_] := Sum[Coefficient[Coefficient[e, s1, a], s2, b] s1^a s2^b, {a, 0, Nd}, {b, 0, Nd - a}];
    AppendTo[res, {i, j, Max[Abs[Flatten[CoefficientList[trn[d1], {s1, s2}]]]], Max[Abs[Flatten[CoefficientList[trn[d2], {s1, s2}]]]]}]]],
  {i, 1, 5}, {j, 1, 5}];
Print["(i,j, maxdiff PPZ reading, maxdiff transposed reading) through total degree ", Nd, ":"];
Print[N[res, 5]];
