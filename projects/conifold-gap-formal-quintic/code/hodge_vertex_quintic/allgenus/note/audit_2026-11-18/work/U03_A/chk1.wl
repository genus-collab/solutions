(* sanity checks for appendix B: Jet matrix, Frobenius basis, concomitant, cyclic form *)
Dq[f_] := -(1 - u) D[f, u];
Dn[f_, 0] := f; Dn[f_, n_] := Dq[Dn[f, n - 1]];
A = {-24/625 (1 - u), -2/5 (1 - u), -7/5 (1 - u), -2 (1 - u), u}; (* A0..A4 *)
PF4[f_] := Sum[A[[k + 1]] Dn[f, k], {k, 0, 4}];
(* Jet by formula *)
Jet = Table[Sum[(-1)^(k - 1 - a) Binomial[k - 1 - a, b] Dn[A[[k + 1]], k - 1 - a - b], {k, a + b + 1, 4}], {a, 0, 3}, {b, 0, 3}] // Simplify;
Print["Jet = ", Jet];
Print["Jet antisym: ", Simplify[Jet + Transpose[Jet]] ];
Print["det Jet = ", Simplify[Det[Jet]]];
(* Jet from definition of concomitant with generic functions *)
conc[f_, g_] := Sum[(-1)^j Dn[A[[k + 1]] g, j] Dn[f, k - 1 - j], {k, 1, 4}, {j, 0, k - 1}];
(* Lagrange identity: g PF4 f - f PF4 g = Dq conc *)
Print["Lagrange: ", Simplify[G[u] PF4[F[u]] - F[u] PF4[G[u]] - Dq[conc[F[u], G[u]]]]];
(* compare conc with Jet form *)
Print["conc vs Jet: ", Simplify[conc[F[u], G[u]] - Sum[Dn[F[u], a] Jet[[a + 1, b + 1]] Dn[G[u], b], {a, 0, 3}, {b, 0, 3}]]];
(* power series solutions *)
NN = 14;
sol[k_] := Module[{c, f, eqs, s},
  f = u^k + Sum[c[n] u^n, {n, 3, NN + 3}];
  eqs = CoefficientList[Normal[Series[PF4[f], {u, 0, NN - 1}]], u];
  s = Solve[Thread[eqs == 0], Table[c[n], {n, 3, NN + 2}]][[1]];
  (f /. s) /. c[NN + 3] -> 0];
f0 = sol[0]; f1 = sol[1]; f2 = sol[2];
Print["f0 = ", Series[f0, {u, 0, 5}]];
Print["f1 = ", Series[f1, {u, 0, 5}]];
Print["f2 = ", Series[f2, {u, 0, 5}]];
w1 = f1 + 7/10 f2;
Print["w1 = ", Series[w1, {u, 0, 4}]];
(* conc on power series *)
Do[Print["Conc(f", i, ",f", j, ") = ", Series[conc[{f0, f1, f2}[[i + 1]], {f0, f1, f2}[[j + 1]]], {u, 0, 6}]], {i, 0, 2}, {j, i + 1, 2}];
(* tc in frame (0,0) *)
Print["tc(0,0) = ", Series[w1/f0, {u, 0, 4}]];
(* cyclic operator *)
y0 = f0 + al w1 + be f2; t = w1/y0; td = Dq[t];
Mcyc[F_] := (1/y0) Dq[(1/td) Dq[(u y0^2 td^2) Dq[(1/td) Dq[F/y0]]]];
diff = Mcyc[F[u]] - PF4[F[u]];
cf = Table[Series[Coefficient[diff, Derivative[j][F][u]], {u, 0, 5}], {j, 0, 4}] // Simplify;
Print["Mcyc - PF4 coefficients (should be O(u^6)): ", cf];
(* a non-isotropic plane: y0=f0, y1=f1 should fail *)
y0b = f0; tb = f1/f0; tdb = Dq[tb];
Mb[F_] := (1/y0b) Dq[(1/tdb) Dq[(u y0b^2 tdb^2) Dq[(1/tdb) Dq[F/y0b]]]];
cfb = Table[Series[Coefficient[Mb[F[u]] - PF4[F[u]], Derivative[j][F][u]], {u, 0, 2}], {j, 0, 4}];
Print["non-isotropic plane (f0,f1): ", cfb];
(* log solution: y3 = w1 ell + rho3, with Dq ell = 1 - 1/u ; treat ell as log u *)
rho = Sum[r[n] u^n, {n, 3, NN + 2}];
y3 = w1 Log[u] + rho;
res = Simplify[PF4[y3]];
resl = Coefficient[res, Log[u]]; res0 = Simplify[res - resl Log[u]];
Print["log coefficient order: ", Series[resl, {u, 0, 8}]];
eqs = CoefficientList[Normal[Series[u^3 res0, {u, 0, NN - 2 + 3}]], u];
Print["obstruction coefficients u^-3,u^-2,u^-1: ", eqs[[1 ;; 3]]];
sr = Solve[Thread[eqs[[4 ;; NN + 2]] == 0], Table[r[n], {n, 3, NN + 1}]][[1]];
y3s = y3 /. sr /. r[NN + 2] -> 0;
Print["rho3 = ", Series[(rho /. sr /. r[NN + 2] -> 0), {u, 0, 5}]];
Do[Print["Conc(f", i, ",y3) = ", Series[Simplify[conc[{f0, f1, f2}[[i + 1]], y3s] /. Log[u] -> LL], {u, 0, 4}]], {i, 0, 2}];
Print["Conc(y3,w1) = ", Series[Simplify[conc[y3s, w1] /. Log[u] -> LL], {u, 0, 4}]];
W = Table[Dn[{f0, f1, f2, y3s}[[i]], a], {a, 0, 3}, {i, 1, 4}];
Print["u^2 det Wr_f = ", Series[Simplify[u^2 Det[W] /. Log[u] -> LL], {u, 0, 4}]];
