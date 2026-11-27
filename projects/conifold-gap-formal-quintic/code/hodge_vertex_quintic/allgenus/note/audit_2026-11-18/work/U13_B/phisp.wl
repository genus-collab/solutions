bk[k_] := BernoulliB[2 k]/(2 k (2 k - 1));
NZ = 3; UU = 5;
W[X_] := Sum[u^(m - 1)/m (1/X + X^(5 m - 1)/(5 m - 1)), {m, 1, UU + 1}];
(* checks of lem:con-form (c) *)
Print["W'(1): ", Series[D[W[X], X] /. X -> 1, {u, 0, UU - 1}]];
Print["W''(1) - 5/(1-u): ", Series[(D[W[X], {X, 2}] /. X -> 1) - 5/(1 - u), {u, 0, UU - 1}]];
tauinv = Sum[4 u^n/(5 n + 4), {n, 0, UU + 1}];
Print["W(1) - 5/4 tau^-1: ", Series[W[1] - 5/4 tauinv, {u, 0, UU - 1}]];
Print["W' closed form: ", Series[D[W[X], X] + 1/(u X^2) Log[(1 - u X^5)/(1 - u)], {u, 0, UU - 1}] // Simplify];
loga = -1/2 Log[1 - u X^5] + Sum[bk[k] (u X zh)^(2 k - 1) (5^(1 - 2 k) - 5 Sum[Binomial[2 k - 2 + 5 P, 5 P] u^P X^(5 P), {P, 0, UU}]), {k, 1, NZ + 1}];
(* stationary phase at X=1 with phase F=-W *)
Ft = Normal[Series[-W[1 + t] + W[1], {t, 0, 2 NZ + 2}]];
H = Coefficient[Ft, t, 2] 2;
F3 = Ft - Coefficient[Ft, t, 2] t^2;
CC = -1/H;
amp = Exp[loga /. X -> 1 + t] / ((1 - u)^(-1/2));
integrand = Exp[(F3 /. t -> h tau)/h^2] (amp /. {t -> h tau, zh -> h^2});
ser = Normal[Series[integrand, {h, 0, 2 NZ}]];
res = Expand[ser] /. {tau^n_ :> If[EvenQ[n], (n - 1)!! CC^(n/2), 0]};
res = res /. tau -> 0;
Do[Print["ph_", n, ": ", Series[Coefficient[res, h, 2 n], {u, 0, UU - 2}]], {n, 0, NZ}];
