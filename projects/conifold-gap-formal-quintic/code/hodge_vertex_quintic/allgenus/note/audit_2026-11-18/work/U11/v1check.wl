nn = 4;
T[k_, F_] := Expand[-F + (t/5) (k F - 4 t D[F, t])];
F0 = 1 + Sum[f[n] t^n, {n, 1, nn}];
F1 = T[-1, F0]; F2 = -T[-2, F1]; F3 = T[2, F2]; F4 = -T[1, F3];
cl = Expand[F0 + T[0, F4]];
sol = First[Solve[Table[Coefficient[cl, t, n] == 0, {n, 1, nn}], Table[f[n], {n, 1, nn}]]];
Fs = Normal[Series[#, {t, 0, nn - 1}]] & /@ ({F0, F1, F2, F3, F4} /. sol);
eps = {1, -1, -1, 1, 1}; nu = {0, -1, 3, 2, 1};
z5 = Exp[2 Pi I/5];
M0 = Table[Sum[(1/5) eps[[b + 1]] z5^(-i (b + 1)) z5^(j (b + 1)) (Fs[[b + 1]] /. t -> z5^(-j) t), {b, 0, 4}], {i, 0, 4}, {j, 0, 4}];
M01 = Coefficient[Expand[M0], t, 1];
s = 7/5; v = s^5;
d = Table[z5^(4 i) s, {i, 0, 4}];
R1 = M01/s^4;
one = Table[1/d[[i]], {i, 1, 5}];
vec[a_] := Table[s^(nu[[a + 1]] - 1) eps[[a + 1]] z5^(i (a + 1)), {i, 0, 4}];
(* genus-0 four-point, psi^0 everywhere, degree-1 part *)
g04[ph_] := Module[{l, pr},
  Sum[Sum[-(R1.ph[[l]])[[i]] Product[If[m == l, 1, ph[[m]][[i]]], {m, 1, 4}], {l, 1, 4}] d[[i]]^2 + Product[ph[[m]][[i]], {m, 1, 4}] (R1.one)[[i]] d[[i]]^3, {i, 1, 5}] +
  Sum[With[{pa = pr[[1]], pb = pr[[2]]}, Sum[ph[[pa[[1]]]][[i]] ph[[pa[[2]]]][[i]] d[[i]] R1[[i, j]] ph[[pb[[1]]]][[j]] ph[[pb[[2]]]][[j]] d[[j]], {i, 1, 5}, {j, 1, 5}]], {pr, {{{1, 2}, {3, 4}}, {{1, 3}, {2, 4}}, {{1, 4}, {2, 3}}}}]];
(* TFT four-point *)
t04[ph_] := Sum[Product[ph[[m]][[i]], {m, 1, 4}] d[[i]]^2, {i, 1, 5}];
(* q_phi via residues at infinity *)
Dx[f_] := Together[p^2/(p^5 - v) D[f, p]];
xi[a_] := Module[{m, ep, al}, {m, ep, al} = {{4, -1, 0}, {0, 1, 1}, {1, 1, 0}, {2, -1, 0}, {3, -1, 0}}[[a + 1]]; ep v^al p^m/(p^5 - v)];
xiad[a_, dd_] := Nest[Dx, xi[a], dd];
pair[phi_, a_, dd_] := Module[{form = phi D[xiad[a, dd], p]}, -(-SeriesCoefficient[form /. p -> 1/w, {w, 0, -1}] (-1))];
(* residue at infinity of f(p)dp = -coefficient of 1/w in f(1/w)/w^2 ... compute carefully below *)
resinf[f_] := -SeriesCoefficient[(f /. p -> 1/w)/w^2, {w, 0, -1}];
pairR[phi_, a_, dd_] := -resinf[phi D[xiad[a, dd], p]];
phi = p^10;
qphi = Table[{a, dd, pairR[phi, a, dd]}, {a, 0, 4}, {dd, 0, 3}];
Print["q_{p^10} components (a,d,coef): ", Select[Flatten[qphi, 1], #[[3]] =!= 0 &]];
H = vec[1];
total = Sum[With[{a = c[[1]], dd = c[[2]], co = c[[3]]}, co If[dd == 0, g04[{vec[a], H, H, H}], If[dd == 1, t04[{vec[a], H, H, H}], 0]]], {c, Flatten[qphi, 1]}];
Print["<q_{p^10} H H H>_0 = ", FullSimplify[ComplexExpand[total]], "   V1 prediction (-5)*(-60) = 300"];
(* also check q_{p^5} = -5(z 1 + v H) *)
Print["q_{p^5}: ", Select[Flatten[Table[{a, dd, pairR[p^5, a, dd]}, {a, 0, 4}, {dd, 0, 3}], 1], #[[3]] =!= 0 &], "  v = ", v];
