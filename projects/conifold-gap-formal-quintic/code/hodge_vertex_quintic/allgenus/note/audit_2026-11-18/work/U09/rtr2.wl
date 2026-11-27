(* Exact check of prop:mod-Rtr (i) at branch point j=0 (c=1), kappa=1, v=1, and of (2-32) at i=0 *)
z5 = Exp[2 Pi I/5];
kap = 1;
xk[p_] := kap (p^4/4 + 1/p);
sq = Sqrt[-5 kap];
ord = 9;
c = 1;
pp = Module[{a, q, ser, sol},
   q = c + w/(sq c) + Sum[a[k] w^k, {k, 2, ord}];
   ser = Series[xk[c] - xk[q] - w^2/2, {w, 0, ord + 1}];
   sol = Solve[Table[SeriesCoefficient[ser, k] == 0, {k, 3, ord + 1}], Table[a[k], {k, 2, ord}]][[1]];
   Simplify[q /. sol]];
xiF[i_][p_] := -1/(sq z5^(-i) (p - z5^(-i)));
G[m_] := Which[OddQ[m], 0, m >= 0, (m - 1)!! s^(m/2), True, (-1)^(-m/2) s^(m/2)/((-m - 1)!!)];  (* s = 1/theta *)
Lap[form_] := Module[{ser}, ser = Series[form, {w, 0, ord - 3}];
   Sum[SeriesCoefficient[ser, m] G[m], {m, -2, ord - 3}]];
RinvRow = Table[Expand[-s Lap[(D[xiF[i][y], y] /. y -> pp) D[pp, w]]], {i, 0, 4}];
(* model *)
N0 = 6;
Tk[k_, F_] := -F + t/5 (k F - 4 t D[F, t]);
F0 = 1 + Sum[cc[n] t^n, {n, 1, N0}];
F1 = Tk[-1, F0]; F2 = -Tk[-2, F1]; F3 = Tk[2, F2]; F4 = -Tk[1, F3];
eq = Series[F0 + Tk[0, F4], {t, 0, N0 + 1}];
sol = Solve[Table[SeriesCoefficient[eq, n] == 0, {n, 2, N0 + 1}], Table[cc[n], {n, 1, N0}]][[1]];
Fs = Normal[Series[#, {t, 0, N0}]] & /@ ({F0, F1, F2, F3, F4} /. sol);
Rhat[tt_] := Table[z5^(i (b + 1)) (Fs[[b + 1]] /. t -> z5^(-i) tt), {b, 0, 4}, {i, 0, 4}];
Psi0 = Rhat[0];
Mt = Inverse[Psi0] . Rhat[t];
(* compare R^{-1}(s)^0_i (j=0 row) with M0(t)_{i0} at t=-s (note) and with M0_{0i} (transposed) *)
note = Table[N[Coefficient[RinvRow[[i + 1]], s, k] - (-1)^k Coefficient[Mt[[i + 1, 1]], t, k], 25], {i, 0, 4}, {k, 0, 3}];
trans = Table[N[Coefficient[RinvRow[[i + 1]], s, k] - (-1)^k Coefficient[Mt[[1, i + 1]], t, k], 25], {i, 0, 4}, {k, 0, 3}];
Print["note reading residuals (should be 0): ", Chop[note, 10^-20]];
Print["transposed reading residuals: ", Chop[trans, 10^-20]];
(* (2-32) at i=0: sum_k R^{-1 0}_k Delta_k^{-1/2} = L_0[dp] ; Delta_k^{-1/2} = 1/(sq varsigma_k) *)
lhs = Expand[Sum[RinvRow[[k + 1]]/(sq z5^(-k)), {k, 0, 4}]];
rhs = Expand[Lap[D[pp, w]]];
Print["(2-32) residuals: ", Chop[Table[N[Coefficient[lhs - rhs, s, k], 25], {k, 0, 3}], 10^-20]];
