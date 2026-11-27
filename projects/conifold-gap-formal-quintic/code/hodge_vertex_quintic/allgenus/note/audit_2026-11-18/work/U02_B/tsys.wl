Nt = 9;
z5 = Exp[2 Pi I/5];
T[i_, k_][F_] := -zz^i F + (t/5) (k F - 4 t D[F, t]);
(* colour-0 system, zz -> 1 *)
F0 = 1 + Sum[a[n] t^n, {n, 1, Nt}];
chain[i_, F_] := Module[{F1, F2, F3, F4, F0b},
  F1 = T[i, -1][F]; F2 = -T[i, -2][F1]; F3 = T[i, 2][F2]; F4 = -T[i, 1][F3]; F0b = -T[i, 0][F4];
  {F, F1, F2, F3, F4, F0b}];
ch = chain[0, F0] /. zz -> 1;
eqs = Table[Coefficient[Expand[ch[[6]] - ch[[1]]], t, n] == 0, {n, 1, Nt + 1}];
sol = Solve[eqs, Table[a[n], {n, 1, Nt}]];
Print["num sols: ", Length[sol]];
Fs = Table[Series[ch[[b]] /. sol[[1]], {t, 0, Nt}] // Normal, {b, 1, 5}];
Print["F_0 = ", Series[Fs[[1]], {t, 0, 5}]];
Print["F_b(0) = ", Fs /. t -> 0];
(* QRR h for comparison *)
bk[k_] := BernoulliB[2 k]/(2 k (2 k - 1));
w0 = Join[Table[1 - z5^j, {j, 1, 4}], {-5}];
h = Series[Exp[-Sum[bk[k] t^(2 k - 1) Total[w0^(1 - 2 k)], {k, 1, 5}]], {t, 0, 6}] // FullSimplify;
Print["h = ", h];
(* M0 *)
eps = {1, -1, -1, 1, 1};
Rhat = Table[z5^(i (b + 1)) (Fs[[b + 1]] /. t -> z5^(-i) t), {b, 0, 4}, {i, 0, 4}];
Psi0 = Rhat /. t -> 0;
Print["Psi0 check: ", FullSimplify[Psi0 - Table[eps[[b + 1]] z5^(i (b + 1)), {b, 0, 4}, {i, 0, 4}]] // Flatten // Union];
etainv = Table[If[Mod[a + b, 5] == 3, 1/5, 0], {a, 0, 4}, {b, 0, 4}];
Print["Psi0^T etainv Psi0 = 1? ", FullSimplify[Transpose[Psi0].etainv.Psi0 - IdentityMatrix[5]] // Flatten // Union];
nn = Transpose[Rhat /. t -> -t].etainv.Rhat;
nnS = Map[Chop[N[Normal[Series[#, {t, 0, Nt - 1}]], 30], 10^-20] &, nn, {2}];
Print["Rhat(-t)^T etainv Rhat(t) (numerical, through t^", Nt - 1, ") = ", nnS];
