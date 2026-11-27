(* E5: sanity check of the corrected DNOPS (2-33) = DOSS (4-1) on the model curve x = W_v (kappa=1), y = p, v = 1.
   Also checks the bidifferential identity D(B) := d1(B/dx1)+d2(B/dx2) = sum_k dxi^k (x) dxi^k exactly. *)
prec = 60;
v = 1; kap = 1;
zeta = Exp[2 Pi I/5];
c[j_] := zeta^(-j);              (* varsigma = 1 *)
s5 = I Sqrt[5 kap];              (* one root sqrt(-5 kappa) for all j *)
x[p_] := kap (p^4/4 + v/p);
(* --- exact identity D(B) = sum_k dxi^k dxi^k (symbolic in p1,p2, exact zeta) --- *)
f[p_] := 1/D[x[q], q] /. q -> p;
Bt[p1_, p2_] := 1/(p1 - p2)^2;
DB = D[f[p1] Bt[p1, p2], p1] + D[f[p2] Bt[p1, p2], p2];
dxi[k_, p_] := 1/(s5 c[k] (c[k] - p)^2);   (* coefficient of dp in dxi^k *)
S = Sum[dxi[k, p1] dxi[k, p2], {k, 0, 4}];
Print["D(B) - sum dxi dxi at random points: ",
  N[(DB - S) /. {p1 -> 3/10 + 7 I/10, p2 -> -11/10 + 2 I/5}, 30],
  "  ", N[(DB - S) /. {p1 -> 2 + I/3, p2 -> 1/7 - 5 I/4}, 30]];
(* --- local inversion p_j(w), x = x(c_j) - w^2/2 --- *)
ord = 14;
pser[j_] := Module[{cj = N[c[j], prec], ser, wt, inv},
  ser = Series[-2 (x[cj + t] - x[cj])/(t^2 (-5 kap cj^2)), {t, 0, ord}];
  wt = N[s5, prec] cj t Sqrt[ser];
  inv = InverseSeries[wt, w];
  Expand[cj + Normal[inv]]];
Ps = Table[pser[j], {j, 0, 4}];
G[m_, psi_] := Which[OddQ[m], 0, m >= 0, (m - 1)!! psi^(m/2), m == -2, -1/psi, True, (-1)^(-m/2) psi^(-m/2)/((-m - 1)!!)];
nmax = 3;
(* single Laplace of g(p) dp at c_j, truncated at psi^nmax *)
Lap[g_, j_] := Module[{ser = Series[(g /. p -> Ps[[j + 1]]) D[Ps[[j + 1]], w], {w, 0, 2 nmax + 1}], cs},
  cs = CoefficientList[Normal[ser] w^4, w];  (* shift by w^4 to keep negative powers *)
  Sum[cs[[m + 5]] G[m, psi], {m, -4, Length[cs] - 5}]];
Rinv = Table[Expand[-psi Lap[dxi[i, p], j]], {j, 0, 4}, {i, 0, 4}];   (* Rinv[[j+1,i+1]] = R^{-1}(psi)^j_i *)
trunc[e_, n_] := Normal[Series[e, {psi, 0, n}]];
Rinv = Map[trunc[#, nmax] &, Rinv, {2}];
Print["R^{-1}(0) = I ? ", Chop[(Rinv /. psi -> 0) - IdentityMatrix[5], 10^-40] == 0 IdentityMatrix[5]];
(* double Laplace of a bidifferential given by a bivariate series in w1,w2 *)
LL[ser_] := Module[{cl = CoefficientList[ser, {w1, w2}]},
  Sum[cl[[a + 1, b + 1]] G[a, p1s] G[b, p2s], {a, 0, Length[cl] - 1}, {b, 0, Length[cl[[1]]] - 1}]];
Breg[i_, j_] := Module[{P1 = Ps[[i + 1]] /. w -> w1, P2 = Ps[[j + 1]] /. w -> w2, Q, Nn, res},
  If[i != j,
    res = Normal[Series[D[P1, w1] D[P2, w2]/(P1 - P2)^2, {w1, 0, 2 nmax}, {w2, 0, 2 nmax}]],
    Q = Cancel[(P1 - P2)/(w1 - w2)];
    Nn = Expand[D[P1, w1] D[P2, w2] - Q^2];
    Nn = PolynomialQuotient[Nn, (w1 - w2)^2, w1];
    res = Normal[Series[Nn/Q^2, {w1, 0, 2 nmax}, {w2, 0, 2 nmax}]]];
  res];
check[i_, j_] := Module[{lhs, rhs, rhsPrinted, diff, diffP},
  lhs = LL[Breg[i, j]];
  rhs = (KroneckerDelta[i, j] - Sum[(Rinv[[i + 1, k + 1]] /. psi -> p1s) (Rinv[[j + 1, k + 1]] /. psi -> p2s), {k, 0, 4}]);
  rhsPrinted = Sum[(Rinv[[i + 1, k + 1]] /. psi -> p1s) (Rinv[[j + 1, k + 1]] /. psi -> p2s), {k, 0, 4}];
  (* compare (p1s+p2s)*lhs with rhs, total degree <= nmax *)
  diff = Expand[(p1s + p2s) lhs - rhs];
  diffP = Expand[(p1s + p2s) lhs - rhsPrinted];
  {Max[Abs[Flatten[Table[Coefficient[Coefficient[diff, p1s, a], p2s, b] /. {p1s -> 0, p2s -> 0}, {a, 0, nmax}, {b, 0, nmax - a}]]]],
   Max[Abs[Flatten[Table[Coefficient[Coefficient[diffP, p1s, a], p2s, b] /. {p1s -> 0, p2s -> 0}, {a, 0, nmax}, {b, 0, nmax - a}]]]]}];
Do[Print["(i,j)=", {i, j}, "  max|corrected residual|, max|printed residual| (total degree <= ", nmax, "): ", N[check[i, j], 8]], {i, 0, 4}, {j, i, 4}];
