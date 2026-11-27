bk[k_] := BernoulliB[2 k]/(2 k (2 k - 1));
NZ = 3; KT = 2 NZ + 3; KK = NZ + 1;
p[k_] := p[k] = Together[SeriesCoefficient[5 (L + t)^4/(1 - (L + t)^5), {t, 0, k - 1}]];
F3 = Sum[(5 (-1)^(k + 1) L^(-k) + p[k])/(k (k + 1)) t^(k + 1), {k, 2, 2 NZ + 2}];
logA = Sum[((-1)^(k + 1) L^(-k) + p[k])/(2 k) t^k, {k, 1, 2 NZ + 2}] +
   Sum[z^(2 k - 1) Sum[bk[k] Binomial[1 - 2 k, r] (5^(1 - 2 k) L^(1 - 2 k - r) + (-1)^r p[2 k - 1 + r]) t^r, {r, 0, 2 NZ + 2}], {k, 1, KK}];
CC = -L (1 - L^5)/5;
integrand = Exp[(F3 /. t -> h tau)/h^2] Exp[logA /. {t -> h tau, z -> h^2}];
ser = Normal[Series[integrand, {h, 0, 2 NZ}]];
ser = Expand[ser];
gauss[expr_] := Module[{e = Expand[expr]}, e /. tau^n_. :> If[EvenQ[n], (n - 1)!! CC^(n/2), 0]];
res = gauss[ser /. tau -> tau];
(* careful: constant term tau^0 *)
res = Expand[ser] /. {tau^n_ :> If[EvenQ[n], (n - 1)!! CC^(n/2), 0]};
res = res /. tau -> 0;
cs = Table[Together[Coefficient[res, h, 2 n]], {n, 0, NZ}];
Print["c_n: ", cs];
(* check odd powers vanish *)
Print["odd: ", Table[Together[Coefficient[res, h, 2 n + 1]], {n, 0, NZ - 1}]];
(* ph_n(u) = L^{-4n} c_n with L = u^{-1/5} *)
Do[Print["ph_", n, " = ", Factor[Together[(L^(-4 n) cs[[n + 1]]) /. L -> u^(-1/5)]]], {n, 0, NZ}];
(* PF5^L check *)
cser = Sum[cs[[n + 1]] z^n, {n, 0, NZ}];
DD[f_] := D[f, L] (L^6 - L)/5;
op[f_] := z DD[f] + L f;
f0 = L cser;
pow[f_, n_] := Nest[op, f, n];
q = (1 - L^(-5))/5^5;
prodterm = Fold[Function[{acc, k}, 5 op[acc] + k z acc], f0, Range[5]];
pf = pow[f0, 5] - f0 - q prodterm;
Print["PF5^L(Lc) coefficients: ", Table[Together[Coefficient[Expand[pf], z, n]], {n, 0, NZ}]];
