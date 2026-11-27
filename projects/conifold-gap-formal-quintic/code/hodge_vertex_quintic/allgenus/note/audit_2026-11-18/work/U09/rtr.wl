(* Direct check of prop:mod-Rtr (i): R^{-1}(theta^{-1})^j_i := -theta^{-1} L_{c_j}[d xi^i] equals M0(t)_{ij}, t=-1/(kappa varsigma^4 theta), kappa=1, v=1 *)
$MaxExtraPrecision = 200;
prec = 50;
z5 = Exp[2 Pi I/5];
kap = 1;
xk[p_] := kap (p^4/4 + 1/p);
cs = Table[N[z5^(-b), prec], {b, 0, 4}];
sq = N[Sqrt[-5 kap], prec];
ord = 9;
(* p as series in w at c: x(c) - x(p) = w^2/2, dp/dw(0) = 1/(sq c) *)
pOfW[c_] := Module[{a, pp, ser, sol},
  pp = c + w/(sq c) + Sum[a[k] w^k, {k, 2, ord}];
  ser = Series[xk[c] - xk[pp] - w^2/2, {w, 0, ord + 1}];
  sol = Solve[Table[SeriesCoefficient[ser, k] == 0, {k, 3, ord + 1}], Table[a[k], {k, 2, ord}]][[1]];
  pp /. sol];
pws = Table[pOfW[cs[[j]]], {j, 1, 5}];
xiF[i_][p_] := -1/(sq cs[[i]] (p - cs[[i]]));
G[m_] := Which[OddQ[m], 0, m >= 0, (m - 1)!! th^(-m/2), True, (-1)^(-m/2) th^(-m/2)/((-m - 1)!!)];
Lap[j_, form_] := Module[{ser, cl}, (* form = f(w) dw given as f *)
  ser = Series[form, {w, 0, ord - 2}];
  cl = Table[SeriesCoefficient[ser, m], {m, -4, ord - 2}];
  Sum[cl[[m + 5]] G[m], {m, -4, ord - 2}]];
Rinv = Table[Module[{pp = pws[[j]], f},
    f = (D[xiF[i][y], y] /. y -> pp) D[pp, w];
    Expand[-th^(-1) Lap[j, f]]], {j, 1, 5}, {i, 1, 5}];
(* model M0(t) *)
N0 = 6;
Tk[k_, F_] := -F + t/5 (k F - 4 t D[F, t]);
F0 = 1 + Sum[c[n] t^n, {n, 1, N0}];
F1 = Tk[-1, F0]; F2 = -Tk[-2, F1]; F3 = Tk[2, F2]; F4 = -Tk[1, F3];
eq = Series[F0 + Tk[0, F4], {t, 0, N0 + 1}];
sol = Solve[Table[SeriesCoefficient[eq, n] == 0, {n, 2, N0 + 1}], Table[c[n], {n, 1, N0}]][[1]];
Fs = Normal[Series[#, {t, 0, N0}]] & /@ ({F0, F1, F2, F3, F4} /. sol);
Rhat[tt_] := Table[z5^(i (b + 1)) (Fs[[b + 1]] /. t -> z5^(-i) tt), {b, 0, 4}, {i, 0, 4}];
Psi0 = Rhat[0];
Mt = Expand[N[Inverse[Psi0] . Rhat[t], prec]];
(* compare: Rinv[[j,i]] vs Mt[[i,j]] with t = -1/th, i.e. th^{-1} = -t *)
diff = Table[Chop[Expand[(Rinv[[j, i]] /. th -> -1/t) - Mt[[i, j]]], 10^-25], {j, 1, 5}, {i, 1, 5}];
diffT = Table[Chop[Normal[Series[diff[[j, i]], {t, 0, 3}]], 10^-25], {j, 1, 5}, {i, 1, 5}];
Print["max discrepancy through t^3: ", Max[Abs[Flatten[Table[Coefficient[diffT[[j, i]], t, k], {j, 5}, {i, 5}, {k, 0, 3}]]]]];
(* alternative (transposed) reading *)
diffAlt = Table[Chop[Normal[Series[Expand[(Rinv[[j, i]] /. th -> -1/t) - Mt[[j, i]]], {t, 0, 3}]], 10^-25], {j, 1, 5}, {i, 1, 5}];
Print["transposed reading discrepancy through t^3: ", Max[Abs[Flatten[Table[Coefficient[diffAlt[[j, i]], t, k], {j, 5}, {i, 5}, {k, 0, 3}]]]]];
(* (2-32): sum_k Rinv^i_k Delta_k^{-1/2} = L_i[dy] *)
lhs = Table[Sum[Rinv[[i, k]]/(sq cs[[k]]), {k, 5}], {i, 5}];
rhs = Table[Expand[Lap[i, D[pws[[i]], w]]], {i, 5}];
Print["(2-32) discrepancy: ", Max[Abs[Flatten[Table[Coefficient[Expand[lhs[[i]] - rhs[[i]]], th, -k], {i, 5}, {k, 0, 3}]]]]];
