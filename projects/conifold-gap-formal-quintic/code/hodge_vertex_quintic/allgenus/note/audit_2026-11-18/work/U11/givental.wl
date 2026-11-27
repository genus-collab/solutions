(* colour-0 T-system *)
nn = 4;
T[k_, F_] := Expand[-F + (t/5) (k F - 4 t D[F, t])];
F0 = 1 + Sum[f[n] t^n, {n, 1, nn}];
F1 = T[-1, F0]; F2 = -T[-2, F1]; F3 = T[2, F2]; F4 = -T[1, F3];
cl = Expand[F0 + T[0, F4]];
eqs = Table[Coefficient[cl, t, n] == 0, {n, 1, nn + 1}];
sol = First[Solve[Take[eqs, nn], Table[f[n], {n, 1, nn}]]];
Fs = {F0, F1, F2, F3, F4} /. sol;
Fs = Normal[Series[#, {t, 0, nn - 1}]] & /@ Fs;
Print["Fs_0 = ", Fs[[1]]];
eps = {1, -1, -1, 1, 1};
z5 = Exp[2 Pi I/5];
(* M0(t)_{ij} = sum_b (1/5) eps_b zeta^{-i(b+1)} zeta^{j(b+1)} F_b(zeta^{-j} t) *)
M0 = Table[Sum[(1/5) eps[[b + 1]] z5^(-i (b + 1)) z5^(j (b + 1)) (Fs[[b + 1]] /. t -> z5^(-j) t), {b, 0, 4}], {i, 0, 4}, {j, 0, 4}];
M01 = Simplify[Coefficient[Expand[M0], t, 1]];
M01 = FullSimplify[ComplexExpand[M01]];
Print["M0_1 symmetric? ", Simplify[M01 - Transpose[M01]] == 0 * M01];
sympl = Simplify[Expand[(M0 /. t -> -t)\[Transpose].M0]];
Print["symplecticity check coefficient t^1, t^2: ", Simplify[Coefficient[sympl, t, 1]], Simplify[Coefficient[sympl, t, 2]]];
(* model data at given varsigma = s (real positive), v = s^5 *)
s = 7/5; v = s^5;
d = Table[z5^(4 i) s, {i, 0, 4}];
R1 = M01/s^4;
one = Table[1/d[[i + 1]], {i, 0, 4}];
phiH = Table[-s^(-2) z5^(2 i), {i, 0, 4}];
genus1[phi_] := (1/24) Sum[-(R1.phi)[[i]] d[[i]] + phi[[i]] d[[i]]^2 (R1.one)[[i]], {i, 1, 5}] + (1/2) Sum[phi[[i]] R1[[i, i]] d[[i]], {i, 1, 5}];
g1 = FullSimplify[ComplexExpand[genus1[phiH]]];
Print["<H>_1 = ", g1, "   -1/(12 v) = ", -1/(12 v)];
g03 = FullSimplify[ComplexExpand[Sum[phiH[[i]]^3 d[[i]], {i, 1, 5}]]];
Print["<HHH>_0 = ", g03, "   -5/v = ", -5/v];
g04 = Sum[4 (-(R1.phiH)[[i]] phiH[[i]]^3 d[[i]]^2) + phiH[[i]]^4 (R1.one)[[i]] d[[i]]^3, {i, 1, 5}] + 3 Sum[phiH[[i]]^2 d[[i]] R1[[i, j]] phiH[[j]]^2 d[[j]], {i, 1, 5}, {j, 1, 5}];
g04 = FullSimplify[ComplexExpand[g04]];
Print["<HHHH>_0 = ", g04, "   5/v^2 = ", 5/v^2];
