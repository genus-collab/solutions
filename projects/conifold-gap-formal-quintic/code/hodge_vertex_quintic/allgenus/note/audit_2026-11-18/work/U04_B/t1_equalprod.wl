(* check prop:gau-equalprod (ii): v = u tau^(-5/4) satisfies v'^5/v = 1/(u(1-u)^5) *)
Nn = 12;
tau = 1/Sum[4 u^n/(5 n + 4), {n, 0, Nn + 2}];
v = Normal[Series[u*tau^(-5/4), {u, 0, Nn}]];
lhs = Series[D[v, u]^5/v, {u, 0, Nn - 3}];
rhs = Series[1/(u (1 - u)^5), {u, 0, Nn - 3}];
Print["EP difference: ", Normal[lhs - rhs] // Simplify];
Print["tau = ", Series[tau, {u, 0, 3}]];
Print["v = ", Series[v, {u, 0, 4}]];
(* 2F1 check *)
Print["2F1 check: ", Series[1/Hypergeometric2F1[1, 4/5, 9/5, u] - tau, {u, 0, 8}]];
