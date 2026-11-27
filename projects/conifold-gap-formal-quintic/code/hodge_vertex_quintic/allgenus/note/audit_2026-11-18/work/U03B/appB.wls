(* Appendix B checks in u; D = -(1-u) d/du; ell formal log u with D ell = 1 - 1/u *)
N0 = 14;
tr[f_] := Expand[Normal[Series[f, {u, 0, N0}]]];
DD[f_] := Expand[-(1 - u) D[f, u] + (1 - 1/u) D[f, ell]];
A = {-24/625 (1 - u), -2/5 (1 - u), -7/5 (1 - u), -2 (1 - u), u};
PF4[f_] := Module[{ds}, ds = NestList[DD, f, 4]; Expand[Sum[A[[k + 1]] ds[[k + 1]], {k, 0, 4}]]];
(* power series solutions *)
sol[b0_, b1_, b2_] := Module[{cs, g, eqs, ss},
  cs = Table[c[n], {n, 3, N0 + 3}];
  g = b0 + b1 u + b2 u^2 + Sum[c[n] u^n, {n, 3, N0 + 3}];
  eqs = Table[Coefficient[PF4[g], u, m] == 0, {m, 0, N0}];
  ss = Solve[eqs, cs][[1]];
  Expand[Normal[Series[g /. ss, {u, 0, N0}]]]];
f0 = sol[1, 0, 0]; f1 = sol[0, 1, 0]; f2 = sol[0, 0, 1];
Print["PF4 f0 low: ", Table[Coefficient[PF4[f0], u, m], {m, -3, N0 - 4}]];
w1 = Expand[f1 + 7/10 f2];
(* log solution *)
rs = Sum[r[n] u^n, {n, 3, N0 + 3}];
y3g = w1 ell + rs;
eq3 = PF4[y3g];
eqs3 = Table[Coefficient[eq3 /. ell -> 0, u, m] == 0, {m, -3, N0 - 4}];
ss3 = Solve[eqs3, Table[r[n], {n, 3, N0 + 3}]];
Print["log sol solvable: ", Length[ss3]];
rho3 = Expand[Normal[Series[rs /. ss3[[1]], {u, 0, N0}]]];
y3 = w1 ell + rho3;
(* concomitant, by definition *)
Conc[f_, g_] := Expand[Sum[Sum[(-1)^j Nest[DD, A[[k + 1]] g, j] Nest[DD, f, k - 1 - j], {j, 0, k - 1}], {k, 1, 4}]];
(* jet matrix: use symbolic jets *)
Clear[F, G];
(* jet matrix from formula *)
Jf[a_, b_] := Sum[(-1)^(k - 1 - a) Binomial[k - 1 - a, b] Nest[DD, A[[k + 1]], k - 1 - a - b], {k, a + b + 1, 4}];
J = Table[Expand[Jf[a, b]], {a, 0, 3}, {b, 0, 3}];
Print["J = ", J];
Print["det J = ", Factor[Det[J]]];
Print["J^-1 = ", Simplify[Inverse[J]]];
(* check C(f,g) = sum D^a f J_ab D^b g on test functions *)
tf = 1 + 3 u + 5 u^3 + u^5; tg = u^2 - 7 u^4 + ell u^3;
lhs = Conc[tf, tg]; rhs = Expand[Sum[Nest[DD, tf, a] J[[a + 1, b + 1]] Nest[DD, tg, b], {a, 0, 3}, {b, 0, 3}]];
Print["jet check: ", Expand[lhs - rhs]];
(* Gram entries: constants; take low-order coefficient and check nonconstant parts vanish up to truncation *)
cc[f_, g_] := Module[{c = Conc[f, g]}, {Coefficient[c /. ell -> 0, u, 0], Table[Coefficient[c /. ell -> 0, u, m], {m, -3, 6}], Coefficient[c, ell]}];
Print["C(f0,f1): ", cc[f0, f1]];
Print["C(f0,f2): ", cc[f0, f2]];
Print["C(f1,f2): ", cc[f1, f2]];
Print["C(f0,w1): ", cc[f0, w1]];
Print["C(f0,y3): ", cc[f0, y3]];
Print["C(f1,y3): ", cc[f1, y3]];
Print["C(f2,y3): ", cc[f2, y3]];
(* Wronskian standard basis *)
Wf = Table[Nest[DD, {f0, f1, f2, y3}[[i]], a], {a, 0, 3}, {i, 1, 4}];
dw = Expand[u^2 Det[Wf]];
Print["u^2 det Wr_f low terms: ", Table[Coefficient[dw /. ell -> 0, u, m], {m, -4, 3}], " ell-part: ", Table[Coefficient[Coefficient[dw, ell], u, m], {m, -4, 3}]];
(* frame (0,0) tc *)
tc = Normal[Series[w1/f0, {u, 0, 5}]];
Print["tc frame(0,0): ", tc];
(* cyclic operator check: compose operators as action on generic jets; compare coefficient of operators *)
Mcyc[y0_, y1_, h_] := Module[{t, td, s1, s2, s3, s4},
  t = Normal[Series[y1/y0, {u, 0, N0 - 1}]];
  td = Expand[DD[t]];
  s1 = DD[h/y0];
  s2 = DD[s1/td];
  s3 = DD[u y0^2 td^2 s2];
  s4 = DD[s3/td];
  s4/y0];
(* apply to generic polynomial test function and compare with PF4 *)
Do[
  y0 = f0 + al w1 + be f2;
  test = u^3 + 2 u^5 - u^6;
  diff = Normal[Series[Mcyc[y0, w1, test] - PF4[test], {u, 0, 4}]];
  Print["Mcyc-PF4 on test, (al,be)=", {al, be}, ": ", Simplify[diff]],
  {al, {0, 3/7}}, {be, {0, -2}}];
(* Note Mcyc applied to h with ell: *)
y0 = f0; testl = ell u^2 + u;
Print["Mcyc-PF4 on ell test: ", Simplify[Normal[Series[Mcyc[y0, w1, testl] - PF4[testl], {u, 0, 3}]]]];
