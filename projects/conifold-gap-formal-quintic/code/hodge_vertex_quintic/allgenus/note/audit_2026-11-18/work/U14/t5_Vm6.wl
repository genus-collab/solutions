(* Test of prop:con-res + thm:con-Pb at z^{-6}: V_{z^-6} = varpi0 * y0^(6) (frame (0,0)); any free factor s(z)=1+s5 z^5 would add s5*varpi0*y0^(11) = O(u^9) *)
K = 12;
th[f_] := -(1 - u) D[f, u];
prodop[f_] := Module[{y = f}, Do[y = th[y] + k/5 y, {k, 4, 1, -1}]; y];
PF4[f_] := Nest[th, f, 4] - (1 - u) prodop[f];
ff = a[0] + Sum[a[n] u^n, {n, 1, K + 3}];
res = Expand[PF4[ff]];
eqs = Table[Coefficient[res, u, n] == 0, {n, 0, K}];
frob[k_] := Module[{fix, s},
  fix = Table[a[j] -> If[j == k, 1, 0], {j, 0, 2}];
  s = Solve[eqs /. fix, Table[a[n], {n, 3, K + 3}]];
  If[Length[s] != 1, Print["frob problem ", k, s]];
  Normal[Series[ff /. fix /. s[[1]], {u, 0, K}]]];
f0 = frob[0]; f1 = frob[1]; f2 = frob[2];
Print["f0 = ", Series[f0, {u, 0, 5}]];
vp0 = f0; vp1 = f1 + 7/10 f2;
tc = Normal[Series[vp1/vp0, {u, 0, K}]];
Print["tc = ", Series[tc, {u, 0, 4}]];
tcp = D[tc, u];
bq = {tcp, -1/(u vp0^2 (1 - u)^3 tcp^2), tcp, -vp0/(1 - u), -vp0/(1 - u)};
bq = Map[Normal[Series[#, {u, 0, K}]] &, bq];
(* note bq[[2]] has a 1/u pole *)
y = {0, -1, 0, 0, 0}; (* y^(0) = -e_1 ; components a = 0..4 at positions 1..5 *)
ys = {y};
Do[
  ynew = Table[Module[{integrand},
     integrand = Normal[Series[bq[[a + 1]] ys[[-1]][[Mod[a + 1, 5] + 1]], {u, 0, K}]];
     Integrate[integrand, {u, 0, u}]], {a, 0, 4}];
  AppendTo[ys, Expand[ynew]], {n, 1, 11}];
Print["y^(1)_0 = ", Series[ys[[2]][[1]], {u, 0, 4}]];
y60 = ys[[7]][[1]]; y110 = ys[[12]][[1]];
Print["ord y^(6)_0: ", Exponent[y60, u, Min], "   ord y^(11)_0: ", Exponent[y110, u, Min]];
rhs6 = Normal[Series[vp0 y60, {u, 0, K}]];
(* flat residue V at z^-6 *)
bk[k_] := BernoulliB[2 k]/(2 k (2 k - 1));
Xi = Sum[u^n/n, {n, 1, K + 1}];
G = Sum[x^(5 m - 1)/(m (5 m - 1)), {m, 1, 3}];
logf = Sum[x^(5 n)/(2 n), {n, 1, 3}] + Sum[bk[k] (x z)^(2 k - 1) (5^(1 - 2 k) KroneckerDelta[P, 0] - 5 Binomial[2 k - 2 + 5 P, 5 P]) x^(5 P), {k, 1, 7}, {P, 0, 2}];
ef = Normal[Series[Exp[-G/z + logf], {x, 0, K}]];
V6 = Sum[(-Xi)^k/k! Coefficient[Expand[Coefficient[ef, x, k - 1] /. (x -> 0)], z, k - 6], {k, 1, K}];
(* careful with k-1=0 coefficient *)
efc[j_] := If[j == 0, ef /. x -> 0, Coefficient[ef, x, j]];
V6 = Sum[(-Xi)^k/k! Coefficient[Expand[efc[k - 1] z^(-k)], z, -6], {k, 1, K}];
V6 = Normal[Series[V6, {u, 0, K}]];
Print["V_{z^-6}       = ", V6];
Print["varpi0 y0^(6)  = ", rhs6];
Print["difference     = ", Expand[V6 - rhs6]];
Print["varpi0 y0^(11) = ", Normal[Series[vp0 y110, {u, 0, K}]]];
