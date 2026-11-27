(* independent sanity check of def:gam-c: compute c_n(L) by Wick contraction *)
NN = ToExpression[$ScriptCommandLine[[2]]];
Tmax = 3 NN + 2;
bk[k_] := BernoulliB[2 k]/(2 k (2 k - 1));
s0pp[x_] := 5/(x (1 - x^5));
s1p[x_] := 1/(2 x) + 5 x^4/(2 (1 - x^5));
f[y_] := 5 y^4/(1 - y^5);
pn[x_, n_] := D[f[x], {x, n - 1}]/(n - 1)!;
g = Normal[Series[s0pp[L + t], {t, 0, Tmax}]];
h2 = s0pp[L];
S0c = Sum[Coefficient[g, t, k] t^(k + 2)/((k + 1) (k + 2)), {k, 1, Tmax}]; (* order >=3 part *)
g1 = Normal[Series[s1p[L + t], {t, 0, Tmax}]];
S1t = Sum[Coefficient[g1, t, k] t^(k + 1)/(k + 1), {k, 0, Tmax}];
S2k[k_] := Normal[Series[bk[k] ((5 (L + t))^(1 - 2 k) + Evaluate[pn[x, 2 k - 1] /. x -> L + t]), {t, 0, Tmax}]];
logamp = S1t + Sum[z^(2 k - 1) S2k[k], {k, 1, Ceiling[NN/2] + 1}];
(* substitute t = r tau, z = r^2 *)
expo = (S0c /. t -> r tau)/r^2 + (logamp /. {t -> r tau, z -> r^2});
ser = Normal[Series[Exp[expo], {r, 0, 2 NN}]];
CC = -1/h2;
wick[expr_] := Module[{e = Expand[expr]}, Sum[Coefficient[e, tau, 2 j] (2 j - 1)!! CC^j, {j, 0, 3 NN + 2}] + (e /. tau -> 0) - (e /. tau -> 0)];
res = Table[Together[Factor[wick[Coefficient[ser, r, 2 n] ] /. {0^0 -> 1}]], {n, 0, NN}];
(* careful: j=0 term *)
Do[Print["c_", n, " = ", Factor[res[[n + 1]]]], {n, 0, NN}];
Print["values at L=1: ", Table[Limit[res[[n + 1]], L -> 1], {n, 0, NN}]];
hser = Series[Exp[Sum[bk[k] (5^(1 - 2 k) + Sum[(Exp[2 Pi I j/5] - 1)^(1 - 2 k), {j, 1, 4}]) z^(2 k - 1), {k, 1, NN}]], {z, 0, NN}];
Print["h_n: ", Table[FullSimplify[SeriesCoefficient[hser, n]], {n, 0, NN}]];
Print["values at L=0: ", Table[Limit[res[[n + 1]], L -> 0], {n, 0, NN}]];
