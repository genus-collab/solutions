(* U12_B sanity checks of displayed identities used in s07 *)
z5 = Exp[2 Pi I/5];
nu = {0, -1, 3, 2, 1}; eps = {1, -1, -1, 1, 1};
eta = Table[If[Mod[a + b, 5] == 3, 5, 0], {a, 0, 4}, {b, 0, 4}];
etai = Inverse[eta];
Print["eta^-1 == eta/25: ", Simplify[etai == eta/25]];
Psi0 = Table[eps[[b + 1]] z5^(i (b + 1)), {b, 0, 4}, {i, 0, 4}];
Print["Psi0^T eta^-1 Psi0 == 1: ", FullSimplify[Transpose[Psi0].etai.Psi0] == IdentityMatrix[5]];
(* T-system colour 0 *)
nmax = 7;
Tk[k_, F_] := -F + (t/5) (k F - 4 t D[F, t]);
F0 = 1 + Sum[a[n] t^n, {n, 1, nmax}];
F1 = Tk[-1, F0]; F2 = -Tk[-2, F1]; F3 = Tk[2, F2]; F4 = -Tk[1, F3]; cl = -Tk[0, F4] - F0;
eqs = Table[Coefficient[Expand[cl], t, n] == 0, {n, 1, nmax + 1}];
sol = Solve[eqs, Table[a[n], {n, 1, nmax}]][[1]];
Fs = {F0, F1, F2, F3, F4} /. sol;
Fs = Map[Normal[Series[#, {t, 0, nmax - 1}]] &, Fs];
Print["Fs0 = ", Fs[[1]]];
Rhat[tt_] := Table[z5^(i (b + 1)) (Fs[[b + 1]] /. t -> z5^(-i) tt), {b, 0, 4}, {i, 0, 4}];
M0[tt_] := Inverse[Psi0].Rhat[tt];
symp = Expand[Transpose[M0[-x]].M0[x]];
symp = Map[Normal[Series[#, {x, 0, nmax - 2}]] &, symp, {2}];
Print["M0(-t)^T M0(t) - 1 (mod t^", nmax - 1, ") vanishes: ", Simplify[FullSimplify[symp - IdentityMatrix[5]] == 0 IdentityMatrix[5]]];
(* Psi^m with formal varsigma *)
Lv = DiagonalMatrix[Table[vs^(1 - nu[[b + 1]]), {b, 0, 4}]];
Psim = Inverse[Lv].Psi0;
Print["Psi^m isometry: ", FullSimplify[Transpose[Psim].etai.Psim] == IdentityMatrix[5]];
unitm = Table[vs^-1 z5^i, {i, 0, 4}];
Print["Psi^m unit -> lowered comps of H^0: ", FullSimplify[Psim.unitm]];
(* Psi_Q with m_a; P0 = Psi_Q Psi^m^-1 diagonal? *)
mm = {L/C0, L^2/(C0 C1), C0 C1/L^2, C0/L, 1};
PsiQ = Table[z5^(j (a + 1)) mm[[a + 1]], {a, 0, 4}, {j, 0, 4}];
P0 = FullSimplify[PsiQ.Inverse[Psim]];
Print["P0 = Psi_Q (Psi^m)^-1 = ", P0];
Print["Psi_Q isometry: ", FullSimplify[Transpose[PsiQ].etai.PsiQ] == IdentityMatrix[5]];
unitQ = Table[z5^j L/C0, {j, 0, 4}];
Print["Psi_Q 1_TFT: ", FullSimplify[PsiQ.unitQ]];
