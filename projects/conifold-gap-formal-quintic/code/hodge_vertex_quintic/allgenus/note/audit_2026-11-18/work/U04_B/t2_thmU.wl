(* U04_B sanity test 2: Theorem U by brute-force linear algebra on a generic pair of one-pole systems with (EP) *)
N0 = 9; K = 3;
bt = {1, -1/u, 1, -1, -1};
b0 = 1 + u; b2 = 1 - 2 u + 3 u^2; b3 = -1 + 3 u - u^3; b4 = -2 + u + u^2;
b1 = Normal[Series[(-1/u)/(b0 b2 b3 b4), {u, 0, N0 + 1}]];
b = {b0, b1, b2, b3, b4};
mk[c_] := Table[If[Mod[i, 5] + 1 == j, c[[i]], 0], {i, 5}, {j, 5}];
A = mk[b]; At = mk[bt];
P[k_] := Table[Sum[p[k, a, bb, n] u^n, {n, 0, N0}], {a, 5}, {bb, 5}];
Pm1 = ConstantArray[0, {5, 5}];
eqs = {};
Do[
  prev = If[k == 0, Pm1, D[P[k - 1], u]];
  ex = Expand[u (prev - A . P[k] + P[k] . At)];
  Do[AppendTo[eqs, Coefficient[ex[[i, j]], u, m]], {i, 5}, {j, 5}, {m, 0, N0}],
  {k, 0, K}];
vars = Flatten[Table[p[k, a, bb, n], {k, 0, K}, {a, 5}, {bb, 5}, {n, 0, N0}]];
M = Normal[CoefficientArrays[eqs, vars][[2]]];
ns = NullSpace[M];
Print["unknowns ", Length[vars], " equations ", Length[eqs], " nullity ", Length[ns]];
pos[k_, a_, bb_, n_] := Position[vars, p[k, a, bb, n]][[1, 1]];
(* rank of the projection to (P_k, u^j), k<=kmax, j<=J *)
Do[
  idx = Flatten[Table[pos[k, a, bb, n], {k, 0, kmax}, {a, 5}, {bb, 5}, {n, 0, J}]];
  Print["kmax=", kmax, " J=", J, " rank of projection = ", MatrixRank[ns[[All, idx]]]],
  {kmax, 0, K - 1}, {J, 0, 4}];
(* P0 of a solution versus dfrak*Theta_U *)
d = FoldList[Normal[Series[#1*#2, {u, 0, 6}]] &, 1, Table[Normal[Series[bt[[a]]/b[[a]], {u, 0, 6}]], {a, 1, 4}]];
dfr = Normal[Series[(Times @@ d/(Times @@ d /. u -> 0))^(-1/5), {u, 0, 4}]];
pred = Table[Normal[Series[dfr d[[a]], {u, 0, 4}]], {a, 5}];
idx0 = Flatten[Table[pos[0, a, bb, n], {a, 5}, {bb, 5}, {n, 0, 4}]];
proj = RowReduce[ns[[All, idx0]]];
row = proj[[1]];
P0 = Table[Sum[row[[(a - 1)*25 + (bb - 1)*5 + n + 1]] u^n, {n, 0, 4}], {a, 5}, {bb, 5}];
Print["P0 offdiag zero: ", Union[Flatten[Table[If[a != bb, P0[[a, bb]], 0], {a, 5}, {bb, 5}]]]];
Print["P0 diag / pred: ", Table[Series[P0[[a, a]]/pred[[a]], {u, 0, 4}], {a, 5}]];
