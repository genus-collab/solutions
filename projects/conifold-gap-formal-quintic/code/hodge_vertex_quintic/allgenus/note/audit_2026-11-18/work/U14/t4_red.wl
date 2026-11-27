(* Test lem:con-red + lem:con-id(b): PhiSP = Sum_b eps_b pi_{nu_b} tau^{(1-nu_b)/4} Fs_b(tau zh), to u^Ku, zh^Kz *)
Ku = 2; Kz = 6;
bk[k_] := BernoulliB[2 k]/(2 k (2 k - 1));
tr[e_] := Normal[Series[e, {u, 0, Ku}]];
(* uniformising map to order Ku *)
yser = 1 + Sum[u^k Ycoef[k][T], {k, 1, Ku}];
Vser = 1 + Sum[u^k Vc[k], {k, 1, Ku}];
rhs = Sum[u^(m - 1)/m (1 + T^m/(5 m - 1)), {m, 1, Ku + 1}];
lhs = Vser/yser + T yser^4/4;
sol = {};
Do[eq = Expand[SeriesCoefficient[(lhs - rhs) /. sol, {u, 0, k}]];
  poly = Sum[c[k, j] T^j, {j, 0, k + 2}];
  s = Solve[CoefficientList[Expand[eq /. Ycoef[k][T] -> poly], T] == 0, Join[Table[c[k, j], {j, 0, k + 2}], {Vc[k]}]];
  newY = poly /. s[[1]];
  AppendTo[sol, Ycoef[k][T] -> newY]; AppendTo[sol, Vc[k] -> (Vc[k] /. s[[1]])], {k, 1, Ku}];
Vv = tr[Vser /. sol];
Yof[X_] := tr[X (yser /. sol) /. T -> X^5];
(* inverse map X(Y) = Y(1 + Sum u^k Xk(Y^5)) *)
Xinv = Y (1 + Sum[u^k Sum[d[k, j] Y^(5 j), {j, 0, k}], {k, 1, Ku}]);
eqinv = Expand[tr[Yof[Xinv] - Y]];
sinv = Solve[Flatten[CoefficientList[eqinv, {u, Y}]] == 0, Flatten[Table[d[k, j], {k, 1, Ku}, {j, 0, k}]]];
Xinv = Xinv /. sinv[[1]];
Print["X(Y) = ", Xinv];
(* hatted amplitude *)
logah = -1/2 Log[1 - u X^5] + Sum[bk[k] (u X zh)^(2 k - 1) (5^(1 - 2 k) - 5 Sum[Binomial[2 k - 2 + 5 P, 5 P] u^P X^(5 P), {P, 0, Ku}]), {k, 1, Ku}];
ah = tr[Exp[logah]];
fY = Expand[tr[(ah /. X -> Xinv) D[Xinv, Y]]];
Print["f_Y polynomial in Y? exponents: ", Union[Exponent[fY, Y, List]] // Short];
(* reduction: Y^n -> (n-3) zh Y^(n-4) + V Y^(n-5) for n >= 4 *)
reduce[expr_] := Module[{e = Expand[expr], n},
  While[(n = Exponent[e, Y]) >= 4,
    cc = Coefficient[e, Y, n];
    e = Expand[tr[e - cc Y^n + cc ((n - 3) zh Y^(n - 4) + Vv Y^(n - 5))]]];
  e];
red = reduce[fY];
pim[m_] := Coefficient[red, Y, m];
Print["pi_m grading check (k - m mod 5 for each zh^k): ",
  Table[{m, Union[Mod[Exponent[Expand[pim[m]], zh, List] - m, 5]]}, {m, -1, 3}]];
(* T-system Fs_b *)
Tsys[k_, F_] := -F + t/5 (k F - 4 t D[F, t]);
F0 = 1 + Sum[a[n] t^n, {n, 1, Kz + 2}];
F1 = Tsys[-1, F0]; F2 = -Tsys[-2, F1]; F3t = Tsys[2, F2]; F4 = -Tsys[1, F3t];
clos = Expand[F0 + Tsys[0, F4]];
sa = Solve[Table[Coefficient[clos, t, n] == 0, {n, 1, Kz + 2}], Table[a[n], {n, 1, Kz + 2}]][[1]];
Fs = Map[Normal[Series[# /. sa, {t, 0, Kz}]] &, {F0, F1, F2, F3t, F4}];
nu = {0, -1, 3, 2, 1}; eps = {1, -1, -1, 1, 1};
tau = tr[1/Sum[4 u^n/(5 n + 4), {n, 0, Ku + 1}]];
formula = Sum[eps[[b + 1]] pim[nu[[b + 1]]] tr[tau^((1 - nu[[b + 1]])/4)] (Fs[[b + 1]] /. t -> tau zh), {b, 0, 4}];
formula = Normal[Series[tr[Expand[formula]], {zh, 0, Kz}]];
(* direct PhiSP by Wick at X = 1 *)
What = Sum[u^(m - 1)/m (1/X + X^(5 m - 1)/(5 m - 1)), {m, 1, Ku + 1}];
Fser = tr[Normal[Series[-(What /. X -> 1 + x), {x, 0, 2 Kz + 2}]]];
c2 = Coefficient[Fser, x, 2]; Cov = tr[-1/(2 c2)];
F3 = Expand[Fser - (Fser /. x -> 0) - Coefficient[Fser, x, 1] x - c2 x^2];
Print["F'(1) = ", Coefficient[Fser, x, 1], "   -1/F''(1) = ", Cov];
mom[j_] := If[OddQ[j], 0, (j - 1)!! Cov^(j/2)];
amp = Normal[Series[ah /. X -> 1 + x, {x, 0, 2 Kz + 2}]];
e = Exp[(F3 /. x -> w tt)/w^2] (amp /. x -> w tt /. zh -> w^2);
s = Normal[Series[e, {w, 0, 2 Kz}]];
s = Expand[tr[Expand[s]]];
s = s /. tt^j_. :> mom[j]; s = s /. tt -> 0;
direct = Expand[tr[Expand[s /. w -> Sqrt[zh]] tr[(1 - u)^(1/2)]]];
direct = Normal[Series[direct, {zh, 0, Kz}]];
Print["PhiSP direct (zh^0..3): ", Table[Coefficient[direct, zh, n], {n, 0, 3}]];
Print["difference formula - direct: ", Expand[tr[formula - direct]]];
