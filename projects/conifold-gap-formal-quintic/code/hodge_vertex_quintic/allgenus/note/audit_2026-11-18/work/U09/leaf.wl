(* Leaf localisation test: with legs R^{-1}(psi) (PPZ) and R^{-1}(psi)^j_i = M0(-psi)_{ij} (manuscript reading),
   the leaf of vertex colour j, sum_i sum_k [psi^{D-k}] Rinv_{ji} d(D_x^k xi~^i), should have poles only at vs_j.
   Compare with the transposed reading. v=1, vs=1, kappa=1. *)
Tk[k_, F_] := -F + t/5 (k F - 4 t D[F, t]);
F0 = 1 + Sum[al[j] t^j, {j, 1, 6}];
F1 = Tk[-1, F0]; F2 = -Tk[-2, F1]; F3 = Tk[2, F2]; F4 = -Tk[1, F3]; F0c = -Tk[0, F4];
eqs = Table[SeriesCoefficient[F0c - F0, {t, 0, j}] == 0, {j, 1, 6}];
alsol = Solve[eqs, Table[al[j], {j, 1, 5}]][[1]];
Fs = Normal[Series[{F0, F1, F2, F3, F4} /. alsol, {t, 0, 4}]];
zeta = Exp[2 Pi I/5]; eps = {1, -1, -1, 1, 1};
Psi0inv = Table[(1/5) eps[[b + 1]] zeta^(-i (b + 1)), {i, 0, 4}, {b, 0, 4}];
Rhat = Table[zeta^(i (b + 1)) (Fs[[b + 1]] /. t -> zeta^(-i) t), {b, 0, 4}, {i, 0, 4}];
M0 = Psi0inv . Rhat;
Mk[n_] := Mk[n] = N[Table[SeriesCoefficient[M0[[i, j]], {t, 0, n}], {i, 5}, {j, 5}], 30];
(* M0(-psi) coefficient of psi^n = (-1)^n Mk[n] *)
sv = Table[zeta^(-i), {i, 0, 4}];
xit[i_] := -1/(sv[[i + 1]] (p - sv[[i + 1]]));
Dx[f_] := p^2/(p^5 - 1) D[f, p];
form[i_, k_] := D[Nest[Dx, xit[i], k], p];
leafA[j_, DD_] := Sum[(-1)^(DD - k) Mk[DD - k][[i + 1, j + 1]] form[i, k], {i, 0, 4}, {k, 0, DD}];  (* manuscript *)
leafB[j_, DD_] := Sum[(-1)^(DD - k) Mk[DD - k][[j + 1, i + 1]] form[i, k], {i, 0, 4}, {k, 0, DD}];  (* transposed *)
pp[f_, m_] := Module[{s}, s = Series[f /. p -> N[sv[[m + 1]], 30] + ee, {ee, 0, -1}]; Max[Abs[Table[SeriesCoefficient[s, -r], {r, 1, 12}]]]];
Do[
 Print["D=", DD, " j=", j, "  max principal part off vs_j: manuscript=", Max[Table[If[m == j, 0, pp[leafA[j, DD], m]], {m, 0, 4}]],
   "  transposed=", Max[Table[If[m == j, 0, pp[leafB[j, DD], m]], {m, 0, 4}]]],
 {DD, 1, 3}, {j, 0, 1}];
