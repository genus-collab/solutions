(* Independent check of thm:mod-tr at (h,n)=(0,3),(1,1): model side, v=1, varsigma=1 *)
z5 = Exp[2 Pi I/5];
nu = {0, -1, 3, 2, 1}; eps = {1, -1, -1, 1, 1};
N0 = 6;
Tk[k_, F_] := -F + t/5 (k F - 4 t D[F, t]);
F0 = 1 + Sum[c[n] t^n, {n, 1, N0}];
F1 = Tk[-1, F0]; F2 = -Tk[-2, F1]; F3 = Tk[2, F2]; F4 = -Tk[1, F3];
eq = Series[F0 + Tk[0, F4], {t, 0, N0 + 1}];
sol = Solve[Table[SeriesCoefficient[eq, n] == 0, {n, 2, N0 + 1}], Table[c[n], {n, 1, N0}]][[1]];
Fs = Normal[Series[#, {t, 0, N0}]] & /@ ({F0, F1, F2, F3, F4} /. sol);
Print["Fs0 = ", Fs[[1]]];
Rhat[tt_] := Table[z5^(i (b + 1)) (Fs[[b + 1]] /. t -> z5^(-i) tt), {b, 0, 4}, {i, 0, 4}];
Psi0 = Rhat[0];
M0[tt_] := Simplify[Inverse[Psi0] . Rhat[tt]];
Mt = Expand[N[M0[t], 40]];
chk = Chop[N[Normal[Series[Transpose[Mt /. t -> -t] . Mt, {t, 0, 4}]], 30], 10^-20];
Print["symplectic check: ", chk == IdentityMatrix[5]];
M1 = Chop[N[Coefficient[Mt, t, 1], 30]];
Print["M1 symmetric: ", Max[Abs[M1 - Transpose[M1]]]];
(* r1 = coefficient of z in R^{-1}(z) = M0(-z)^T ; varsigma=1 *)
r1 = -Transpose[M1];
delta = Table[z5^(4 i), {i, 0, 4}];
one = Table[1/delta[[i]], {i, 1, 5}];
(* <x>_1 for x in colour basis *)
om11[x_] := Sum[delta[[i]] x[[i]], {i, 5}];
om12[x_, y_] := Sum[delta[[i]]^2 x[[i]] y[[i]], {i, 5}];
om03[x_, y_, w_] := Sum[delta[[i]] x[[i]] y[[i]] w[[i]], {i, 5}];
E00 = -r1; (* bivector: E00[[l,k]] = -(r1)_{lk} *)
g1[x_] := 1/24 om11[r1 . x] + 1/24 om12[x, -r1 . one] + 1/2 Sum[E00[[l, k]] om03[x, UnitVector[5, l], UnitVector[5, k]], {l, 5}, {k, 5}];
PsiInvH[a_] := Table[Psi0[[a + 1, i + 1]], {i, 0, 4}];
Print["<H^a>_1 (a=0..4): ", Chop[N[Table[g1[PsiInvH[a]], {a, 0, 4}], 20]]];
Print["<psi H^a>_1 (a=0..4): ", Chop[N[Table[1/24 om11[PsiInvH[a]], {a, 0, 4}], 20]]];
(* genus-0 3-point in flat basis *)
g03 = Table[Chop[N[om03[PsiInvH[a], PsiInvH[b], PsiInvH[c]], 20]], {a, 0, 4}, {b, 0, 4}, {c, 0, 4}];
Print["nonzero <H^aH^bH^c>_0: ", Select[Flatten[Table[{a, b, c, g03[[a + 1, b + 1, c + 1]]}, {a, 0, 4}, {b, 0, 4}, {c, 0, 4}], 2], #[[4]] != 0 &]];
Save["model_data.mx", {Fs, M1}];
g03r = Round[Re[g03]];
Put[Hold[g03 = Evaluate[g03r]] // ReleaseHold; Unevaluated[g03 = g03r], "g03tmp.wl"];
Export["g03.wl", "g03 = " <> ToString[g03r, InputForm] <> ";", "Text"];
Print["g03 exported"];
