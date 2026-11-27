(* downstream use of stirfe in Theorem P(a)(i) Step 3: shift identity of the Stirling-form summand, branch-free form *)
K = 5; zeta = N[Exp[2 Pi I/5], 60];
bk[k_] := BernoulliB[2 k]/(2 k (2 k - 1));
L1p[x_] := Sum[(-1)^(n + 1) x^n/n, {n, 1, K + 3}]; (* log(1+x) truncated *)
sv = 7/3;
cj = Table[zeta^j - sv, {j, 0, 4}];
(* (S0(s+z)-S0(s))/z minus the constant 5log5s - sum log(zeta^j - s) *)
ex1 = (5 (sv + z) L1p[z/sv] - 5 z + Sum[(cj[[j + 1]] - z) L1p[-z/cj[[j + 1]]] + z, {j, 0, 4}])/z;
ex2 = L1p[z/sv]/2 - Sum[L1p[-z/cj[[j + 1]]], {j, 0, 4}]/2;
S2k[k_, s_] := bk[k] ((5 s)^(1 - 2 k) + Sum[(zeta^j - s)^(1 - 2 k), {j, 0, 4}]);
ex3 = Sum[z^(2 k - 1) (S2k[k, sv + z] - S2k[k, sv]), {k, 1, K + 1}];
const = (5 sv)^5/(1 - sv^5);
lhs = const Exp[ex1 + ex2 + ex3];
rhs = -Product[5 sv + k z, {k, 1, 5}]/((sv + z)^5 - 1);
d = Series[lhs - rhs, {z, 0, K}];
Print[Chop[Normal[d], 10^-40]]; Print[Table[Abs[SeriesCoefficient[d, n]], {n, 0, K}]];
Print["rhs check ", N[Normal[Series[rhs, {z, 0, 2}]], 10]];
