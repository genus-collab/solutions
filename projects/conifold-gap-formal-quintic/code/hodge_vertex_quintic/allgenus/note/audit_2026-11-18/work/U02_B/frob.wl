(* Check lem:set-frobbasis *)
Dq[f_] := -(1 - u) D[f, u];
PF4[f_] := u Dq[Dq[Dq[Dq[f]]]] - (1 - u) (2 Dq[Dq[Dq[f]]] + 7/5 Dq[Dq[f]] + 2/5 Dq[f] + 24/625 f);
(* P_j(n): coefficient of u^(n+j) in PF4 u^n *)
res = Expand[PF4[u^n]/u^n] // Simplify;
coefs = Table[Factor[Coefficient[Expand[res], u, j]], {j, -3, 1}];
Print["P_{-3..1}(n) = ", coefs];
Print["P_{-3} check: ", Simplify[coefs[[1]] - n (n - 1)^2 (n - 2)]];
Print["P_{-2} check: ", Simplify[coefs[[2]] + n (n - 1) (4 n^2 - 6 n + 17/5)]];
Print["P_{-1}(0), P_{-1}'(0): ", {coefs[[3]] /. n -> 0, D[coefs[[3]], n] /. n -> 0}];
(* solutions *)
Nmax = 8;
sol[b0_, b1_, b2_] := Module[{g, eqs, vars, s},
  g = b0 + b1 u + b2 u^2 + Sum[c[k] u^k, {k, 3, Nmax + 3}];
  eqs = Table[Coefficient[Expand[PF4[g]], u, k] == 0, {k, 0, Nmax}];
  vars = Table[c[k], {k, 3, Nmax + 3}];
  s = Solve[eqs, vars][[1]];
  Series[g /. s, {u, 0, Nmax}] // Normal];
f0 = sol[1, 0, 0]; f1 = sol[0, 1, 0]; f2 = sol[0, 0, 1];
Print["f0 = ", Series[f0, {u, 0, 5}]];
Print["f1 = ", Series[f1, {u, 0, 5}]];
Print["f2 = ", Series[f2, {u, 0, 5}]];
w1 = f1 + 7/10 f2;
Print["varpi1 = ", Series[w1, {u, 0, 5}]];
Print["tc(0,0) = ", Series[w1/f0, {u, 0, 5}]];
(* negative-power check: PF4 of low terms *)
Print["PF4[f0] low: ", Series[PF4[f0], {u, 0, 4}]];
(* tau *)
tau = 1/Sum[4 u^k/(5 k + 4), {k, 0, 10}];
Print["tau = ", Series[tau, {u, 0, 4}]];
v = u Series[tau, {u, 0, 8}]^(-5/4);
Print["v = ", v];
(* EP: v'^5/v = 1/(u(1-u)^5) *)
vn = Normal[v];
Print["EP check: ", Series[D[vn, u]^5/vn - 1/(u (1 - u)^5), {u, 0, 5}]];
(* concomitant of f0, w1 etc using Jet *)
Jet = {{0, 2/5 (1 - u), 1 - u, -u}, {-2/5 (1 - u), 0, u, 0}, {-(1 - u), -u, 0, 0}, {u, 0, 0, 0}};
jet[f_] := NestList[Dq, f, 3];
conc[f_, g_] := Series[jet[f].Jet.jet[g], {u, 0, 4}];
Print["Conc(f0,w1) = ", conc[f0, w1], "  Conc(f2,w1) = ", conc[f2, w1], "  Conc(f0,f1)=", conc[f0, f1], " Conc(f0,f2)=", conc[f0, f2]];
