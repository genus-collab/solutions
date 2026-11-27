(* colour-i T-system; check uniqueness/scaling F^{(i)}_b(t) = zeta^{i b} Fs_b(zeta^{-i} t) to order t^6 *)
N0 = 7;
T[i_, k_][F_] := Expand[-z^i F + t/5 (k F - 4 t D[F, t])];
solve[i_] := Module[{a, F0, F1, F2, F3, F4, eqs, sol},
  F0 = 1 + Sum[a[n] t^n, {n, 1, N0}];
  F1 = T[i, -1][F0]; F2 = -T[i, -2][F1]; F3 = T[i, 2][F2]; F4 = -T[i, 1][F3];
  eqs = Table[Coefficient[Expand[F0 + T[i, 0][F4]] /. z^5 -> 1, t, n] == 0, {n, 1, N0}];
  sol = Solve[eqs /. {z^5 -> 1}, Table[a[n], {n, 1, N0}]];
  {F0, F1, F2, F3, F4} /. First[sol]];
FsZ = solve[0];
Print["Fs_0 to t^3: ", Series[FsZ[[1]], {t, 0, 3}]];
z = Exp[2 Pi I/5];
F1sys = solve[1];
chk = Table[Simplify[Series[F1sys[[b + 1]] - z^b (FsZ[[b + 1]] /. t -> t/z), {t, 0, N0 - 1}] // Normal], {b, 0, 4}];
Print["scaling check (should be zeros): ", Simplify[chk]];
