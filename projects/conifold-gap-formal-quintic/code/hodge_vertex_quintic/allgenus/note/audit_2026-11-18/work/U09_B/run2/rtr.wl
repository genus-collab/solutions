(* Independent check of prop:mod-Rtr (i),(ii) at v = 1 (varsigma = 1), for kappa = 1 and kappa = 4:
   R^{-1}_TR(psi)^j_i := -psi L_{c_j}[d xi^i]  versus  M0(t)_{ij}, t = -psi/kappa. *)
Get["m0.mx"];
Nt = 6; prec = 50;
sq5 = Sqrt[-5];
G[m_] := Which[OddQ[m], 0, m >= 0, (m - 1)!! th^(-m/2), True, (-1)^(-m/2) th^(-m/2)/((-m - 1)!!)];
check[kappa_] := Module[{x, cs, res, sqk, Rinv, dif, Lap, dyL, lhs232},
  sqk = Sqrt[-5 kappa];
  x[p_] := kappa (p^4/4 + 1/p);
  cs = Table[N[z5^(-j), prec], {j, 0, 4}];
  (* local inverse p(w) at c: x(p(w)) - x(c) = -w^2/2, p'(0) = 1/(sqk c) *)
  pser[c_] := Module[{a, pw, eq, s},
    pw = c + w/(sqk c) + Sum[a[k] w^k, {k, 2, 2 Nt + 4}];
    s = {};
    Do[eq = Chop[SeriesCoefficient[Series[x[pw //. s] - x[c] + w^2/2, {w, 0, k + 1}], k + 1], 10^-40];
       s = Join[s, {a[k] -> (a[k] /. Solve[eq == 0, a[k]][[1]])}], {k, 2, 2 Nt + 4}];
    Chop[(pw //. s) /. a[_] -> 0, 10^-40]];
  Lap[form_] := Module[{ser, cl, lo},
    ser = Series[form, {w, 0, 2 Nt + 1}];
    lo = ser[[4]]; cl = ser[[3]];
    Sum[cl[[k]] G[lo + k - 1], {k, 1, Length[cl]}]];
  res = Table[0, {5}, {5}];
  Do[Module[{pw = pser[cs[[j + 1]]], dp},
      dp = D[pw, w];
      Do[res[[j + 1, i + 1]] = Expand[-(1/th) Lap[(1/(sqk cs[[i + 1]])) dp/(pw - cs[[i + 1]])^2]], {i, 0, 4}];
      (* (2-32): RHS L_{c_j}[dy] *)
      dyL[j] = Expand[Lap[dp]];
    ], {j, 0, 4}];
  (* res[[j,i]] is a polynomial in 1/th = psi *)
  Rinv = res /. th -> 1/psi;
  dif = Table[Expand[Normal[Series[Rinv[[j, i]] - (M0[[i, j]] /. t -> -psi/kappa), {psi, 0, Nt - 1}]]], {j, 1, 5}, {i, 1, 5}];
  Print["kappa=", kappa, ": max coefficient of R^{-1}_TR(psi)^j_i - M0(-psi/kappa)_{ij} through psi^", Nt - 1, ": ",
    Max[Abs[Flatten[Map[CoefficientList[#, psi] &, dif, {2}]]]], "   transposed reading: ", Max[Abs[Flatten[Table[CoefficientList[Expand[Normal[Series[Rinv[[j, i]] - (M0[[j, i]] /. t -> -psi/kappa), {psi, 0, Nt - 1}]]], psi], {j, 1, 5}, {i, 1, 5}]]]]];
  Print["  sample R^{-1}_TR(psi)^1_2 = ", Chop[N[Normal[Series[Rinv[[1, 2]], {psi, 0, 2}]], 8]]];
  (* (2-32): sum_k R^{-1 i}_k Delta_k^{-1/2} = L_{c_i}[dy], Delta_k^{-1/2} = 1/(sqk c_k) *)
  lhs232 = Table[Expand[Sum[Rinv[[i, k]]/(sqk cs[[k]]), {k, 1, 5}]], {i, 1, 5}];
  Print["  (2-32) max diff through psi^", Nt - 1, ": ",
    Max[Abs[Flatten[Table[CoefficientList[Normal[Series[lhs232[[i]] - (dyL[i - 1] /. th -> 1/psi), {psi, 0, Nt - 1}]], psi], {i, 1, 5}]]]]];
 ];
check[1];
check[4];
