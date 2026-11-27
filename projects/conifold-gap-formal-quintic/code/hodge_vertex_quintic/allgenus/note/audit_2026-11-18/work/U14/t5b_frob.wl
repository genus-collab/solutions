K = 8;
th[f_] := -(1 - u) D[f, u];
prodop[f_] := Module[{y = f}, Do[y = th[y] + k/5 y, {k, 4, 1, -1}]; y];
PF4[f_] := Nest[th, f, 4] - (1 - u) prodop[f];
ff = Sum[a[n], {n, 0, 0}] + Sum[a[n] u^n, {n, 1, K + 4}];
res = Expand[PF4[ff]];
Do[Print["u^", n, ": ", Collect[Coefficient[res, u, n], Table[a[j], {j, 0, n + 4}]]], {n, 0, 2}];
