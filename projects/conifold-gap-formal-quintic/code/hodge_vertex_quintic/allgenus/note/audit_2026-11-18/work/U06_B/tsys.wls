(* colour-0 T-system of prop:set-modrows(b); compare F_0 with the u=0 values of the HPF solution *)
NN = 7;
T[k_][F_] := -F + (t/5) (k F - 4 t D[F, t]);
F0 = 1 + Sum[a[n] t^n, {n, 1, NN}];
F1 = T[-1][F0]; F2 = -T[-2][F1]; F3 = T[2][F2]; F4 = -T[1][F3];
cl = Expand[-T[0][F4] - F0];
sol = Solve[Table[Coefficient[cl, t, n] == 0, {n, 1, NN + 1}] /. t -> 0, Table[a[n], {n, 1, NN}]];
Print[Table[a[n], {n, 1, NN}] /. sol];
(* colour covariance check of lem:gau-colmod (i) for i=1 symbolically with w=zeta *)
Fs = (F0 /. sol[[1]]);
Fsb = {Fs, T[-1][Fs], -T[-2][T[-1][Fs]]};
w = Exp[2 Pi I/5];
Ti[k_][F_] := -w F + (t/5) (k F - 4 t D[F, t]);
G0 = Fs /. t -> t/w; G1 = w (Fsb[[2]] /. t -> t/w);
Print["G1 - T^(1)_{-1} G0 (low order): ", Simplify[Normal[Series[G1 - Ti[-1][G0], {t, 0, 4}]]]];
