(* End-to-end low-order check of lem:gau-R1, lem:gau-diag, thm:gau-G, prop:gau-symp in a conifold frame *)
NU = 9; KK = 7; MMs = NU + 5;
al = 1/3; be = -2/7;  (* a generic rational frame *)
Dq[f_] := -(1 - u) D[f, u];
PF4[f_] := u Nest[Dq, f, 4] - (1 - u) (2 Nest[Dq, f, 3] + 7/5 Nest[Dq, f, 2] + 2/5 Dq[f] + 24/625 f);
ans = Sum[c[j] u^j, {j, 0, MMs}];
cl = CoefficientList[Expand[PF4[ans]], u];
solc = Solve[Table[cl[[m + 1]] == 0, {m, 0, MMs - 3}], Table[c[j], {j, 3, MMs}]][[1]];
fgen = ans /. solc;
f0 = fgen /. {c[0] -> 1, c[1] -> 0, c[2] -> 0};
f1 = fgen /. {c[0] -> 0, c[1] -> 1, c[2] -> 0};
f2 = fgen /. {c[0] -> 0, c[1] -> 0, c[2] -> 1};
w1 = f1 + 7/10 f2; w0 = f0 + al w1 + be f2;
ser[x_] := Series[x, {u, 0, NU + 2}];
tc = Normal[ser[w1/w0]];
Print["tc (frame 0,0 would be u+7/10u^2+41/75u^3): ", Normal[Series[(f1 + 7/10 f2)/f0, {u, 0, 3}]]];
tcp = D[tc, u];
tau = Normal[ser[1/Sum[4 u^n/(5 n + 4), {n, 0, NU + 4}]]];
vv = Normal[ser[u tau^(-5/4)]]; vp = D[vv, u];
(* u*b_a and u*tilde b_a as polynomials mod u^(NU+2) *)
ub = Normal /@ ser /@ {u tcp, -1/(w0^2 (1 - u)^3 tcp^2), u tcp, -u w0/(1 - u), -u w0/(1 - u)};
utb = Normal /@ ser /@ {u vp, -u vp/vv, u vp, -u vp, -u vp};
Print["EP check: ", Series[Times @@ ub - Times @@ utb, {u, 0, NU}]];
md[a_] := Mod[a, 5];
p[k_, b_] := If[k < 0, 0, Sum[x[k, b, j] u^j, {j, 0, NU}]];
eqs = Flatten[Table[
    CoefficientList[Expand[u D[p[k - 1, b], u] - ub[[md[b + k - 1] + 1]] p[k, b] + utb[[md[b - 1] + 1]] p[k, md[b - 1]]] + u^(NU + 5), u][[1 ;; NU + 1]],
    {k, 0, KK}, {b, 0, 4}]];
vars = Flatten[Table[x[k, b, j], {k, 0, KK}, {b, 0, 4}, {j, 0, NU}]];
eqs = Join[eqs, {x[0, 1, 0] - 1, x[5, 1, 0]}];
sol = Quiet[Solve[Thread[eqs == 0], vars]][[1]];
P[k_, b_] := p[k, b] /. sol;
free = Complement[vars, First /@ sol];
Print["number of free parameters: ", Length[free]];
freeQ[e_] := And @@ (FreeQ[e, #] & /@ free);
(* diag *)
Print["p00 - tau^(-1/4)/w0 : ", Series[P[0, 0] - tau^(-1/4)/w0, {u, 0, 6}]];
Print["p0a(0): ", Table[P[0, a] /. u -> 0, {a, 0, 4}]];
Print["prod p0a: ", Series[Product[P[0, a], {a, 0, 4}], {u, 0, 5}]];
(* R1: column 1 at u=0: p_{k,1}(0)=0 for k>=1; row 2: p_{k,2-k}(0) = delta; P22(0,z) = P11(0,z) *)
Print["p_{k,1}(0), k=0..KK: ", Table[P[k, 1] /. u -> 0, {k, 0, KK}]];
Print["p_{k,2-k}(0), k=0..KK: ", Table[P[k, md[2 - k]] /. u -> 0, {k, 0, KK}]];
(* vanishing solutions *)
eps = {1, -1, -1, 1, 1}; nu = {0, -1, 3, 2, 1};
YCc[b_, n_] := If[Mod[n - nu[[b + 1]] - 1, 5] == 0, With[{m = Floor[n/5]}, eps[[b + 1]] (-1)^n vv^(n - m)/(4^m m! (n - m)!)], 0];
bq = {ub[[1]]/u, ub[[2]]/u, ub[[3]]/u, ub[[4]]/u, ub[[5]]/u} // Expand;
nmax = 14;
y[0] = {0, -1, 0, 0, 0};
Do[y[n + 1] = Table[Integrate[Normal[Series[bq[[a + 1]] y[n][[md[a + 1] + 1]], {u, 0, NU}]], u], {a, 0, 4}], {n, 0, nmax}];
Print["z^-1 coefficient of YCQ: ", Series[y[1] + {tc, 0, 0, 0, 0}, {u, 0, 6}]];
(* (G): coefficient of z^{-n} (n may be negative) of (P YC)_a = sum_k p_{k,a-k} [z^{-n-k}] YC_{a-k} *)
PY[a_, n_] := Sum[If[n + k >= 0, P[k, md[a - k]] YCc[md[a - k], n + k], 0], {k, 0, KK}];
YQ[a_, n_] := If[n >= 0, y[n][[a + 1]], 0];
ordlim[n_] := Ceiling[4 (n + KK + 1)/5] - 1;  (* trustworthy u-order: terms k>KK start at u^ceil(4(n+KK+1)/5) *)
res = Table[
   With[{o = Min[ordlim[n], NU - 1]},
    If[o >= 0, Normal[Series[PY[a, n] - YQ[a, n], {u, 0, o}]], 0]],
   {n, -KK, 6}, {a, 0, 4}];
Print["(G) residuals zero: ", res // Flatten // Union];
(* model is Y^C_{tilde A}: check recursion gives YC(v(u)) *)
tb = Expand[utb/u];
yt[0] = {0, -1, 0, 0, 0};
Do[yt[n + 1] = Table[Integrate[Normal[Series[tb[[a + 1]] yt[n][[md[a + 1] + 1]], {u, 0, NU}]], u], {a, 0, 4}], {n, 0, 8}];
Print["Y^C_{tilde A} = YC(v(u)) : ", Union[Flatten[Table[Normal[Series[yt[n][[b + 1]] - YCc[b, n], {u, 0, NU}]], {n, 0, 8}, {b, 0, 4}]]]];
(* symplecticity *)
eta = Table[If[Mod[a + b, 5] == 3, 5, 0], {a, 0, 4}, {b, 0, 4}];
Pmat[zz_] := Sum[zz^k Sum[P[k, b] SparseArray[{{md[b + k] + 1, b + 1} -> 1}, {5, 5}], {b, 0, 4}], {k, 0, KK}];
S = Transpose[Pmat[-z]] . eta . Pmat[z] - eta // Normal;
Sc = Table[Normal[Series[Coefficient[S, z, m], {u, 0, 3}]], {m, 0, 4}];
Print["symplecticity residual (z^0..z^4, u<=3): ", Union[Flatten[Sc]]];
Print["which entries involve free params at z^5,z^6: ", Table[freeQ[Normal[Series[Coefficient[S, z, m], {u, 0, 1}]]], {m, 5, 6}], "  values: ", Table[Union[Flatten[Normal[Series[Coefficient[S, z, m], {u, 0, 1}]]]], {m, 5, 6}]];
