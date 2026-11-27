(* U04_B sanity test 1: the quintic one-pole system, the pulled-back model, (EP), (N0), Mcyc = PF4 *)
NN = 14;
Dq[f_] := -(1 - u) D[f, u];
PF4[f_] := u Nest[Dq, f, 4] - (1 - u) (2 Nest[Dq, f, 3] + 7/5 Nest[Dq, f, 2] + 2/5 Dq[f] + 24/625 f);
(* power-series solutions *)
bs = Array[bb, NN + 4, 0];
gen = Sum[bb[n] u^n, {n, 0, NN + 3}];
eqs = CoefficientList[Normal[Series[PF4[gen], {u, 0, NN}]], u];
sol = Solve[Thread[eqs == 0], Table[bb[n], {n, 3, NN + 3}]][[1]];
fgen = gen /. sol;
f0 = fgen /. {bb[0] -> 1, bb[1] -> 0, bb[2] -> 0};
f1 = fgen /. {bb[0] -> 0, bb[1] -> 1, bb[2] -> 0};
f2 = fgen /. {bb[0] -> 0, bb[1] -> 0, bb[2] -> 1};
Print["f0 = ", Series[f0, {u, 0, 5}]];
Print["f1 = ", Series[f1, {u, 0, 5}]];
Print["f2 = ", Series[f2, {u, 0, 5}]];
al = 1/3; be = -2/7;
w1 = f1 + 7/10 f2; w0 = f0 + al w1 + be f2;
tc = Normal[Series[w1/w0, {u, 0, NN}]];
Print["tc = ", Series[tc, {u, 0, 4}]];
tcp = D[tc, u];
b = {tcp, -1/(u w0^2 (1 - u)^3 tcp^2), tcp, -w0/(1 - u), -w0/(1 - u)};
Print["prod b * u(1-u)^5 = ", Series[Times @@ b*u (1 - u)^5, {u, 0, NN - 3}]];
tau = 1/Sum[4 u^n/(5 n + 4), {n, 0, NN + 2}];
v = Normal[Series[u tau^(-5/4), {u, 0, NN}]];
Print["v = ", Series[v, {u, 0, 4}]];
vp = D[v, u];
bt = {vp, -vp/v, vp, -vp, -vp};
Print["prod bt * u(1-u)^5 = ", Series[Times @@ bt*u (1 - u)^5, {u, 0, NN - 3}]];
Print["b(0) a=0,2,3,4: ", (b[[{1, 3, 4, 5}]] /. u -> 0), "  bt(0): ", (bt[[{1, 3, 4, 5}]] /. u -> 0)];
Print["res b1 = ", Limit[u b[[2]], u -> 0], "  res bt1 = ", Limit[u bt[[2]], u -> 0]];
(* Mcyc(w0,w1) applied to a test series vs PF4 *)
DtcQ = Dq[tc];
Mcyc[f_] := (1/w0) Dq[(1/DtcQ) Dq[u w0^2 DtcQ^2 Dq[(1/DtcQ) Dq[f/w0]]]];
test = 1 + 3 u - 2 u^2 + 5 u^3 + u^4/7 - u^5 + 2 u^6;
Print["Mcyc - PF4 on test: ", Series[Mcyc[test] - PF4[test], {u, 0, NN - 6}]];
(* Theta_U, and p00 = tau^(-1/4)/w0 *)
d = FoldList[#1*#2 &, 1, Table[bt[[a]]/b[[a]], {a, 1, 4}]];
dser = Normal[Series[#, {u, 0, NN - 3}]] & /@ d;
Print["d_a(0) = ", dser /. u -> 0];
Print["closing d4*bt4/b4 = ", Series[d[[5]] bt[[5]]/b[[5]], {u, 0, NN - 4}]];
detT = Times @@ d;
dd = Series[(detT)^(-1/5), {u, 0, NN - 4}];
Print["dfrak - tau^(-1/4)/w0 = ", dd - Series[tau^(-1/4)/w0, {u, 0, NN - 4}]];
