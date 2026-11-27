(* Check 1: T-system colour 0, M0, and direct Laplace series (prop:mod-Rtr) at v=1, kappa=1 *)
ord = 4;
Tk[k_, F_] := -F + (t/5) (k F - 4 t D[F, t]);
F0 = 1 + Sum[a[n] t^n, {n, 1, ord}];
F1 = Tk[-1, F0]; F2 = -Tk[-2, F1]; F3 = Tk[2, F2]; F4 = -Tk[1, F3];
cl = CoefficientList[Expand[F0 + Tk[0, F4]], t];
sol = Solve[Table[cl[[n + 1]] == 0, {n, 1, ord + 1}], Table[a[n], {n, 1, ord}]][[1]];
Fs = Normal[Series[{F0, F1, F2, F3, F4} /. sol, {t, 0, ord}]];
Print["Fs = ", Fs];
eps = {1, -1, -1, 1, 1};
z5 = Exp[2 Pi I/5];
Rhat[tt_] := Table[z5^(i (b + 1)) (Fs[[b + 1]] /. t -> z5^(-i) tt), {b, 0, 4}, {i, 0, 4}];
Psi0 = Table[eps[[b + 1]] z5^(i (b + 1)), {b, 0, 4}, {i, 0, 4}];
M0[tt_] := Inverse[Psi0] . Rhat[tt];
M0ser = Table[N[Coefficient[M0[tt], tt, n], 40], {n, 0, ord}];
Print["M0(0) = 1? ", Chop[M0ser[[1]] - IdentityMatrix[5], 10^-30] // Flatten // Union];
(* symplectic check M0(-t)^T M0(t) = 1 to order ord *)
Mm[tt_] := Sum[M0ser[[n + 1]] tt^n, {n, 0, ord}];
sy = CoefficientList[Transpose[Mm[-tt]] . Mm[tt] // Expand, tt];
Print["sympl defect (orders 1..ord): ", Table[Max[Abs[Flatten[Coefficient[Expand[Transpose[Mm[-tt]] . Mm[tt]], tt, n]]]], {n, 1, ord}]];
(* direct Laplace series *)
s5 = Sqrt[-5]; W[p_] := p^4/4 + 1/p;
sig[j_] := z5^(-j);
wo = 2 ord + 6;
eser[j_] := Module[{ws},
   ws = s5 sig[j] e Sqrt[Series[-2 (W[sig[j] + e] - W[sig[j]])/(-5 sig[j]^2 e^2), {e, 0, wo}]];
   InverseSeries[N[ws, 60], w]];
G[m_] := Which[OddQ[m], 0, m >= 0, (m - 1)!! th^(-m/2), True, (-1)^(-m/2) th^(-m/2)/((-m - 1)!!)];
Lap[ser_] := Module[{cl, m0}, 
   Sum[SeriesCoefficient[ser, m] G[m], {m, -4, 2 ord}]];
RinvTR = Table[
   Module[{es = eser[j], form},
    form = (1/(s5 sig[i])) D[Normal[es], w]/(sig[i] - sig[j] - es)^2;
    form = Series[form, {w, 0, 2 ord + 1}];
    Expand[-th^(-1) Lap[form]]], {j, 0, 4}, {i, 0, 4}];
(* compare RinvTR[[j,i]] with M0(t)_{ij}, t=-1/th *)
cmp = Table[RinvTR[[j + 1, i + 1]] - Sum[M0ser[[n + 1, i + 1, j + 1]] (-1/th)^n, {n, 0, ord}], {j, 0, 4}, {i, 0, 4}];
Print["max defect R^-1_TR(j,i) - M0(t)_{ij}: ", Max[Abs[Flatten[Table[Coefficient[cmp[[j, i]], th, -n], {j, 5}, {i, 5}, {n, 0, ord}]]]]];
cmpT = Table[RinvTR[[j + 1, i + 1]] - Sum[M0ser[[n + 1, j + 1, i + 1]] (-1/th)^n, {n, 0, ord}], {j, 0, 4}, {i, 0, 4}];
Print["max defect vs transposed M0(t)_{ji}: ", Max[Abs[Flatten[Table[Coefficient[cmpT[[j, i]], th, -n], {j, 5}, {i, 5}, {n, 0, ord}]]]]];
(* antisymmetric part of M0 coefficients *)
Print["antisym part norms of M0_n: ", Table[Max[Abs[Flatten[M0ser[[n + 1]] - Transpose[M0ser[[n + 1]]]]]], {n, 0, ord}]];
(* (2-32): sum_k Rinv^i_k Delta_k^{-1/2} = L_i[dy] *)
lhs232 = Table[Sum[RinvTR[[i + 1, k + 1]]/(s5 sig[k]), {k, 0, 4}], {i, 0, 4}];
rhs232 = Table[Module[{es = eser[i]}, Expand[Lap[Series[D[Normal[es], w], {w, 0, 2 ord + 1}]]]], {i, 0, 4}];
Print["(2-32) defect: ", Max[Abs[Flatten[Table[Coefficient[lhs232[[i]] - rhs232[[i]], th, -n], {i, 5}, {n, 0, ord}]]]]];
DumpSave["m0.mx", {M0ser, Fs}];
