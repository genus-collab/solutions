(* Test of lem:con-cfu: Y = X(1+Sum u^k Y_k(X^5)), What = V/Y + Y^4/4 *)
K = 5;
yser = 1 + Sum[u^k Ycoef[k][T], {k, 1, K}];
Vser = 1 + Sum[u^k Vc[k], {k, 1, K}];
rhs = Sum[u^(m - 1)/m (1 + T^m/(5 m - 1)), {m, 1, K + 1}];
lhs = Vser/yser + T yser^4/4;
sol = {};
Do[
  eq = Expand[SeriesCoefficient[(lhs - rhs) /. sol, {u, 0, k}]];
  (* Ycoef[k][T] unknown polynomial of degree <= k+1 *)
  poly = Sum[c[k, j] T^j, {j, 0, k + 2}];
  eq2 = Expand[eq /. Ycoef[k][T] -> poly];
  cs = CoefficientList[eq2, T];
  s = Solve[cs == 0, Join[Table[c[k, j], {j, 0, k + 2}], {Vc[k]}]];
  If[Length[s] != 1, Print["no unique solution at k=", k, " ", s]];
  newY = poly /. s[[1]];
  AppendTo[sol, Ycoef[k][T] -> newY];
  AppendTo[sol, Ycoef[k] -> Function[{TT}, Evaluate[newY /. T -> TT]]];
  AppendTo[sol, Vc[k] -> (Vc[k] /. s[[1]])];
  Print["k=", k, " Y_k=", Factor[newY], "  deg=", Exponent[newY, T], "  V_k=", Vc[k] /. s[[1]]],
  {k, 1, K}];
tau = 1/Sum[4 u^n/(5 n + 4), {n, 0, K + 2}];
Print["tau^(-5/4) series: ", Series[tau^(-5/4), {u, 0, K}]];
Print["V series:          ", Series[Vser /. sol, {u, 0, K}]];
(* Y'(1) and Y(1) *)
Yfun = X (yser /. sol) /. T -> X^5;
Yc = Series[Yfun /. X -> 1, {u, 0, K}];
Print["Y(1) - tau^(-1/4): ", Series[(Yfun /. X -> 1) - tau^(-1/4), {u, 0, K}]];
Yp1 = D[Yfun, X] /. X -> 1;
Print["Y'(1) Y(1) - (1-u)^(-1/2): ", Series[Yp1 (Yfun /. X -> 1) - (1 - u)^(-1/2), {u, 0, K}]];
(* check identity in X *)
What = Sum[u^(m - 1)/m (1/X + X^(5 m - 1)/(5 m - 1)), {m, 1, K + 1}];
chk = Series[(Vser /. sol)/Yfun + Yfun^4/4 - What, {u, 0, K}];
Print["W_V(Y(X)) - What: ", Simplify[Normal[chk]]];
