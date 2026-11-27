(* E5: leg identity L_{c_j}[B(.,p)] = sum_k R^{-1}(psi)^j_k sum_d psi^d d((d/dx)^d xi^k)(p), and its transpose *)
Get["edge233_core.wl"];
pt = 7/10 - 3 I/10;
xi[k_, p_] := -1/(s5 c[k] (p - c[k]));
Dx[e_] := Simplify[D[e, p]/D[x[p], p]];
legForm[k_, d_] := D[Nest[Dx, xi[k, p], d], p];   (* coefficient of dp *)
LB[j_] := Module[{P = Ps[[j + 1]], ser, cs},
  ser = Normal[Series[D[P, w] /(P - pt)^2, {w, 0, 2 nmax + 1}]];
  cs = CoefficientList[ser, w];
  Sum[cs[[m + 1]] G[m, psi], {m, 0, Length[cs] - 1}]];
Do[
  lhs = trunc[LB[j], nmax];
  rhs = trunc[Sum[Rinv[[j + 1, k + 1]] Sum[psi^d (legForm[k, d] /. p -> pt), {d, 0, nmax}], {k, 0, 4}], nmax];
  rhsT = trunc[Sum[Rinv[[k + 1, j + 1]] Sum[psi^d (legForm[k, d] /. p -> pt), {d, 0, nmax}], {k, 0, 4}], nmax];
  Print["j=", j, " max|leg residual| = ", N[Max[Abs[CoefficientList[Expand[lhs - rhs], psi]]], 5],
    "   transposed reading: ", N[Max[Abs[CoefficientList[Expand[lhs - rhsT], psi]]], 5]], {j, 0, 4}];
