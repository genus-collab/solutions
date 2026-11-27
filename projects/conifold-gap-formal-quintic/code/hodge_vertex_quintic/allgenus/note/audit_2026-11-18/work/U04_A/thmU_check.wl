(* Sanity check of the construction in the proof of Theorem U(a) on a random pair of one-pole systems with equal products. *)
NN = 14; KK = 3;
SeedRandom[20261117];
ser[f_] := Series[f, {u, 0, NN}];
ru := Sum[RandomInteger[{-3, 3}] u^i, {i, 1, 6}];
unit := RandomChoice[{1, 2, -1, 3, -2}] + ru;
bt = Table[If[a == 1, unit/u, unit], {a, 0, 4}];
g = Table[unit, {a, 0, 3}];
b = Table[If[a <= 3, bt[[a + 1]] g[[a + 1]], bt[[5]]/(Times @@ g)], {a, 0, 4}];
idx[a_] := Mod[a, 5] + 1;
E5[a_, c_] := SparseArray[{idx[a], idx[c]} -> 1, {5, 5}];
A = Normal[Sum[b[[idx[a]]] E5[a, a + 1], {a, 0, 4}]];
At = Normal[Sum[bt[[idx[a]]] E5[a, a + 1], {a, 0, 4}]];
Print["EP: ", ser[Times @@ b - Times @@ bt]];
(* Theta_U *)
dU = Table[1, {5}];
Do[dU[[a + 2]] = Together[dU[[a + 1]] bt[[a + 1]]/b[[a + 1]]], {a, 0, 3}];
Print["closing d: ", ser[dU[[5]] bt[[5]]/b[[5]] - 1]];
Theta = DiagonalMatrix[dU];
dlog = DiagonalMatrix[Table[Normal[ser[D[dU[[a]], u]/dU[[a]]]], {a, 1, 5}]];
detT = Times @@ dU;
frakd = Normal[ser[(detT/(detT /. u -> 0))^(-1/5)]];
Print["A Theta = Theta At: ", Normal[ser[A . Theta - Theta . At]] // Flatten // Union];
Bm[m_, a_] := Product[bt[[idx[a + i]]], {i, 0, m - 1}];
Bfull = Bm[5, 0];
Ycal = Inverse[At];
tr[M_] := Map[Normal[ser[#]] &, M, {2}];
(* projection: coefficients phi_j of Y^j *)
prcoef[X_] := Module[{phi = Table[0, {5}]},
   Do[
    With[{j = Mod[5 - m, 5]},
     phi[[j + 1]] = Normal[ser[(1/5) Sum[X[[idx[a], idx[a + m]]] Bm[5 - m, a + m], {a, 0, 4}]]]],
    {m, 1, 5}];
   phi];
lift[X_] := Module[{Z = ConstantArray[0, {5, 5}], xi, sig, rho},
   Do[
    sig = Table[X[[idx[a], idx[a + m + 1]]]/Bm[m + 1, a], {a, 0, 4}];
    rho = Table[0, {5}];
    rho[[idx[1]]] = 0;
    Do[rho[[idx[a + 1]]] = rho[[idx[a]]] + sig[[idx[a]]], {a, 1, 4}];
    Do[Z[[idx[a], idx[a + m]]] = Normal[ser[rho[[idx[a]]] Bm[m, a]]], {a, 0, 4}],
    {m, 0, 4}];
   Z];
trd = Normal[ser[Tr[dlog]]];
chi[j_] := -(j/5) D[Bfull, u]/Bfull + trd/5;
solvew[j_, phi_] := Module[{c, w, eq, sol},
   w = Sum[c[n] u^n, {n, 0, NN - 2}];
   If[j == 0,
    Normal[ser[frakd (0 - Integrate[Normal[ser[phi/frakd]], u])]],
    eq = Normal[Series[u (D[w, u] + chi[j] w + phi), {u, 0, NN - 2}]];
    sol = Solve[Thread[CoefficientList[eq, u] == 0], Table[c[n], {n, 0, NN - 2}]];
    w /. First[sol]]];
Pch = {frakd IdentityMatrix[5]};
Do[
  X = tr[D[Pch[[k + 1]], u] + dlog . Pch[[k + 1]]];
  Print["k=", k, " pr(X) orders: ", Map[Series[#, {u, 0, NN - 2 k - 3}] &, prcoef[X]]];
  Z = lift[X];
  Print["  lift integral? min exponents: ", Min[Flatten[Map[If[# === 0, Infinity, Exponent[#, u, Min]] &, Z, {2}]]]];
  Print["  [At,Z]-X: ", Union[Flatten[Map[Series[#, {u, 0, NN - 2 k - 4}] &, At . Z - Z . At - X, {2}]]]];
  ZZ = tr[D[Z, u] + dlog . Z];
  phi = prcoef[ZZ];
  Print["  phi pole orders: ", Map[If[# === 0, Infinity, Exponent[#, u, Min]] &, phi]];
  w = Table[solvew[j, phi[[j + 1]]], {j, 0, 4}];
  Wm = Sum[w[[j + 1]] MatrixPower[Ycal, j], {j, 0, 4}];
  AppendTo[Pch, tr[Z + Wm]],
  {k, 0, KK - 1}];
P = Table[tr[Theta . Pch[[k + 1]]], {k, 0, KK}];
Print["P0 diag: ", Series[Diagonal[P[[1]]], {u, 0, 2}]];
Do[
  res = If[k == 0, A . P[[1]] - P[[1]] . At, A . P[[k + 1]] - P[[k + 1]] . At - D[P[[k]], u]];
  Print["gauge eq at z^", k, ": ", Union[Flatten[Map[Series[#, {u, 0, NN - 2 k - 5}] &, res, {2}]]]];
  Print["   min exponent of P_k entries: ", Min[Flatten[Map[If[# === 0, Infinity, Exponent[#, u, Min]] &, P[[k + 1]], {2}]]]],
  {k, 0, KK}];
