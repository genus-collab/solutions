(* Conifold sanity checks, appendix B (a2_cyclicPF4) *)
NN = 14;
Dq[e_] := -(1 - u) D[e, u];     (* ell = Log[u]: Dq ell = 1 - 1/u automatically *)
A = {-24/625 (1 - u), -2/5 (1 - u), -7/5 (1 - u), -2 (1 - u), u}; (* A0..A4 *)
PF4[f_] := Sum[A[[k + 1]] Nest[Dq, f, k], {k, 0, 4}];
(* compare with q-form *)
th[e_] := Dq[e];
PF4q[f_] := Nest[th, f, 4] - 5 ((1 - u)/5^5) ((5 th[#] + #) &[(5 th[#] + 2 #) &[(5 th[#] + 3 #) &[(5 th[f] + 4 f)]]]);
tf = Exp[3 u] Sin[u] + u^2 Log[u];
Print["u-form vs q-form: ", Simplify[PF4[tf] - PF4q[tf]]];
(* power series solutions *)
sol[k_] := Module[{c, f, eqs, s},
  f = u^k + Sum[c[n] u^n, {n, 3, NN}];
  eqs = CoefficientList[Normal[Series[PF4[f], {u, 0, NN - 3}]], u];
  s = Solve[Thread[eqs == 0], Table[c[n], {n, 3, NN}]];
  f /. First[s]];
f0 = sol[0]; f1 = sol[1]; f2 = sol[2];
Print["f0 = ", Series[f0, {u, 0, 5}]]; Print["f1 = ", Series[f1, {u, 0, 5}]]; Print["f2 = ", Series[f2, {u, 0, 5}]];
(* log solution: g = f1 + gam f2 *)
Clear[gam, r];
g = f1 + gam f2; h = Sum[r[n] u^n, {n, 3, NN}];
res = Normal[Series[PF4[g Log[u] + h] /. Log[u] -> ell, {u, 0, NN - 4}]];
resnolog = res /. ell -> 0;
obst = Table[Coefficient[resnolog, u, m], {m, -3, -1}];
Print["obstruction coefficients u^-3,u^-2,u^-1: ", Simplify[obst]];
solg = Solve[obst[[3]] == 0, gam]; Print["gamma: ", solg];
gamv = gam /. First[solg];
eqs2 = Table[Coefficient[resnolog /. gam -> gamv, u, m], {m, 0, NN - 4}];
sr = Solve[Thread[eqs2 == 0], Table[r[n], {n, 3, NN - 1}]];
rho3 = (h /. r[NN] -> 0) /. First[sr];
w1 = f1 + gamv f2; y3 = w1 Log[u] + rho3;
Print["w1 = ", Series[w1, {u, 0, 4}]];
Print["PF4 y3 residual order: ", Series[PF4[y3] /. Log[u] -> ell, {u, 0, NN - 6}]];
(* Jet matrix from the general formula *)
Jet = Table[Sum[(-1)^(k - 1 - a) Binomial[k - 1 - a, b] Nest[Dq, A[[k + 1]], k - 1 - a - b], {k, a + b + 1, 4}], {a, 0, 3}, {b, 0, 3}];
Print["Jet = ", Simplify[Jet] // MatrixForm];
Print["det Jet = ", Simplify[Det[Jet]]];
(* concomitant from the definition *)
ConcDef[f_, gg_] := Sum[Sum[(-1)^j Nest[Dq, A[[k + 1]] gg, j] Nest[Dq, f, k - 1 - j], {j, 0, k - 1}], {k, 1, 4}];
ConcJ[f_, gg_] := Sum[Nest[Dq, f, a] Jet[[a + 1, b + 1]] Nest[Dq, gg, b], {a, 0, 3}, {b, 0, 3}];
Print["Conc def vs Jet: ", Simplify[ConcDef[tf, Cos[u] + Log[u]] - ConcJ[tf, Cos[u] + Log[u]]]];
(* Lagrange identity *)
ff = Exp[2 u] + Log[u] u; gg0 = Sin[u] + u^3 Log[u]^2;
Print["Lagrange: ", Simplify[gg0 PF4[ff] - ff PF4[gg0] - Dq[ConcDef[ff, gg0]]]];
basis = {f0, f1, f2, y3};
Gram = Table[Normal[Series[ConcJ[basis[[i]], basis[[j]]] /. Log[u] -> ell, {u, 0, 2}]], {i, 4}, {j, 4}];
Print["Gram (should be constants) = ", Simplify[Gram] // MatrixForm];
Wr = Table[Nest[Dq, basis[[i]], a], {a, 0, 3}, {i, 4}];
Print["u^2 det Wr = ", Series[u^2 Det[Wr] /. Log[u] -> ell, {u, 0, 3}]];
(* cyclic operator *)
Mcyc[y0_, y1_][f_] := Module[{t = y1/y0, td}, td = Dq[t];
   (1/y0) Dq[(1/td) Dq[(u y0^2 td^2) Dq[(1/td) Dq[f/y0]]]]];
al = 3/7; be = -5/11;
w0 = f0 + al w1 + be f2;
tests = {1, u, u^2, u^3, u^4, u^5, Log[u], u Log[u]};
Do[Print["Mcyc - PF4 on ", tt, ": ", Series[(Mcyc[w0, w1][tt] - PF4[tt]) /. Log[u] -> ell, {u, 0, 5}]], {tt, tests}];
(* non-isotropic plane: should fail *)
Do[Print["Mcyc(f0, f1) - PF4 on ", tt, ": ", Series[(Mcyc[f0, f1][tt] - PF4[tt]) /. Log[u] -> ell, {u, 0, 2}]], {tt, {u^2, u^3}}];
(* general isotropic plane with y0 = w1, y1 = f0 *)
Do[Print["Mcyc(w1, f0) - PF4 on ", tt, ": ", Series[(Mcyc[w1, f0][tt] - PF4[tt]) /. Log[u] -> ell, {u, 0, 3}]], {tt, {1, u^2, u^3, Log[u]}}];
