(* Unit U05 sanity check: canonical gauge between two self-dual cyclic systems of the shape of
   prop:gau-quintic / prop:gau-equalprod, with GENERIC rational series in place of the periods
   (the product of the b_a does not depend on them, so (EP) holds). Checks to low order:
   lem:gau-diag (p_{0,0} = tau^(-1/4)/varpi_0), lem:gau-R1, thm:gau-G (G), prop:gau-symp. *)
NN = 11; KK = 9;
SeedRandom[7];
rr[] := RandomInteger[{-5, 5}]/RandomInteger[{1, 7}];
w0 = 1 + Sum[rr[] u^n, {n, 1, NN + 3}];
tc = u + Sum[rr[] u^n, {n, 2, NN + 4}];
tau = Normal@Series[1/Sum[4 u^n/(5 n + 4), {n, 0, NN + 4}], {u, 0, NN + 3}];
vv = Normal@Series[u tau^(-5/4), {u, 0, NN + 3}];
b = {D[tc, u], -1/(u w0^2 (1 - u)^3 D[tc, u]^2), D[tc, u], -w0/(1 - u), -w0/(1 - u)};
bt = {D[vv, u], -D[vv, u]/vv, D[vv, u], -D[vv, u], -D[vv, u]};
bser = Series[#, {u, 0, NN + 1}] & /@ b; btser = Series[#, {u, 0, NN + 1}] & /@ bt;
Print["EP check (prod b - prod bt): ", Series[Times @@ b - Times @@ bt, {u, 0, NN - 2}]];
bc[a_, i_] := bc[a, i] = SeriesCoefficient[bser[[Mod[a, 5] + 1]], i];
btc[a_, i_] := btc[a, i] = SeriesCoefficient[btser[[Mod[a, 5] + 1]], i];
X[k_, c_, n_] := If[k < 0 || n < 0 || n > NN, 0, x[k, Mod[c, 5], n]];
eqs = Flatten@Table[
    Sum[bc[c + k - 1, i] X[k, c, n - i], {i, -1, n}] - 
      Sum[btc[c - 1, i] X[k, c - 1, n - i], {i, -1, n}] - (n + 1) X[k - 1, c, n + 1] == 0,
    {k, 0, KK}, {c, 0, 4}, {n, -1, NN - 1}];
norm = Join[{x[0, 1, 0] == 1}, Table[x[5 m, 1, 0] == 0, {m, 1, Floor[KK/5]}]];
vars = Flatten@Table[x[k, c, n], {k, 0, KK}, {c, 0, 4}, {n, 0, NN}];
sol = Quiet@Solve[Join[DeleteCases[eqs, True], norm], vars];
Print["number of solutions: ", Length[sol]];
sol = First[sol];
det[k_, c_, n_] := NumericQ[x[k, c, n] /. sol];
(* trusted range: largest n0(k) with all x[k,c,n], n<=n0 determined *)
n0[k_] := Module[{n = -1}, While[n < NN && And @@ Table[det[k, c, n + 1], {c, 0, 4}], n++]; n];
Print["determined u-orders per k: ", Table[n0[k], {k, 0, KK}]];
p[k_, c_] := Sum[(x[k, Mod[c, 5], n] /. sol) u^n, {n, 0, n0[k]}];
(* lem:gau-diag *)
Print["p_{0,a}(0): ", Table[p[0, c] /. u -> 0, {c, 0, 4}], "   prod p_{0,a} - 1 = ", 
  Series[Product[p[0, c], {c, 0, 4}] - 1, {u, 0, n0[0]}]];
Print["p_{0,0} - tau^(-1/4)/varpi0 = ", Series[p[0, 0] - tau^(-1/4)/w0, {u, 0, n0[0]}]];
(* matrices P_k: entry (c+k, c) = p_{k,c}; 1-based indices *)
P[k_] := Normal@SparseArray[Table[{Mod[c + k, 5] + 1, c + 1} -> p[k, c], {c, 0, 4}], {5, 5}];
eta = Table[If[Mod[a + bb, 5] == 3, 5, 0], {a, 0, 4}, {bb, 0, 4}];
(* lem:gau-R1 at u=0 *)
P0z = Sum[z^k (P[k] /. u -> 0), {k, 0, KK - 1}];
Print["P(0,z) (mod z^", KK, "): column 1 = ", P0z[[All, 2]], "  row 2 = ", P0z[[3, All]]];
(* symplecticity *)
nmin = Min[Table[n0[k], {k, 0, KK - 1}]];
sym = Table[Series[Sum[(-1)^k Transpose[P[k]].eta.P[m - k], {k, 0, m}] - If[m == 0, eta, 0 eta], {u, 0, nmin}], {m, 0, KK - 1}];
Print["symplecticity P(-z)^T eta P(z) - eta, orders z^0..z^", KK - 1, ", mod u^", nmin + 1, ": ", 
  Union[Flatten[Normal /@ sym]]];
(* vanishing solutions *)
nu = {0, -1, 3, 2, 1}; ep = {1, -1, -1, 1, 1};
MM = 13;
YCm[n_] := Table[If[Mod[n - nu[[a + 1]] - 1, 5] == 0, 
    With[{m = Floor[n/5]}, ep[[a + 1]] (-1)^n vv^(n - m)/(4^m m! (n - m)!)], 0], {a, 0, 4}];
yQ[0] = {0, -1, 0, 0, 0};
yQ[n_] := yQ[n] = Table[
    Module[{ig = Normal@Series[b[[a + 1]] yQ[n - 1][[Mod[a + 1, 5] + 1]], {u, 0, NN}]},
      If[(ig /. u -> 0) === ComplexInfinity || ! PolynomialQ[ig, u], Print["POLE in integrand ", n, a]];
      Integrate[ig, {u, 0, u}]], {a, 0, 4}];
Print["[z^-1] Y^C_Q + tc e_0 = ", Series[yQ[1] + {tc, 0, 0, 0, 0}, {u, 0, NN}]];
(* (G): coefficient of z^j of Pc.Y^C(v(u)) minus that of Y^C_Q, for u-orders that are exact *)
res = Table[
   Module[{nex, lhs, rhs},
    nex = Min[nmin, Ceiling[4 (KK - j)/5] - 1];
    lhs = Sum[If[k - j >= 0 && k - j <= MM, P[k].YCm[k - j], 0 {1, 1, 1, 1, 1}], {k, 0, KK - 1}];
    rhs = If[j <= 0, yQ[-j], {0, 0, 0, 0, 0}];
    {j, nex, Union@Flatten@Normal@Series[lhs - rhs, {u, 0, nex}]}], {j, -8, 6}];
Print["(G) residuals {j, exact u-order checked, residual set}: ", res];
