(* Independent check of prop:mod-Rtr (i) at v=1 (varsigma=1), kappa=1, to order t^nn:
   R_TR^{-1}(1/th)^j_i := -(1/th) L_{c_j}[d xi^i]  (DNOPS (2-31) read through lem:mod-laplace)
   should equal M0(t)_{ij}, t = -1/(kappa th).  Also check M0(-t)^T M0(t) = 1 and Fs_0. *)
nn = 3; prec = 40; z5 = N[Exp[2 Pi I/5], prec];
Tk[k_, F_] := -F + (t/5) (k F - 4 t D[F, t]);
truncK[e_] := Normal[Series[e, {t, 0, nn + 1}]];
trunc[e_] := Normal[Series[e, {t, 0, nn}]];
F0 = 1 + Sum[f[n] t^n, {n, 1, nn + 1}];
F1 = truncK[Tk[-1, F0]]; F2 = truncK[-Tk[-2, F1]]; F3 = truncK[Tk[2, F2]]; F4 = truncK[-Tk[1, F3]];
eqs = CoefficientList[truncK[-Tk[0, F4] - F0], t];
sol = Solve[Thread[Take[eqs, {3, nn + 2}] == 0], Table[f[n], {n, 1, nn}]][[1]];
Fs = trunc /@ Expand[{F0, F1, F2, F3, F4} /. sol];
Print["Fs0 = ", Fs[[1]]];
Rhat = Table[z5^(i (b + 1)) (Fs[[b + 1]] /. t -> z5^(-i) t), {b, 0, 4}, {i, 0, 4}];
Psi0 = Rhat /. t -> 0;
M0 = Map[trunc, Inverse[Psi0].Rhat, {2}];
symp = Map[trunc, Transpose[M0 /. t -> -t].M0, {2}];
Print["max |M0(-t)^T M0(t) - 1|: ", Max[Abs[Flatten[CoefficientList[#, t] & /@ Flatten[symp - IdentityMatrix[5]]]]]];
s5 = N[Sqrt[-5], prec]; xW[p_] := p^4/4 + 1/p;
Lap[ser_] := Module[{e = Expand[Normal[ser]], lo = -4, hi = 2 nn},
   Sum[With[{a = Coefficient[e, w, m]},
      Which[EvenQ[m] && m >= 0, a (m - 1)!! th^(-m/2), EvenQ[m] && m < 0, a (-1)^(-m/2) th^(-m/2)/(-m - 1)!!, True, 0]],
     {m, lo, hi}]];
pw[j_] := Module[{c = z5^(-j), wp},
   wp = Series[s5 c (p - c) Sqrt[Series[-2 (xW[p] - xW[c])/(-5 c^2 (p - c)^2), {p, c, 2 nn + 3}]], {p, c, 2 nn + 3}];
   InverseSeries[wp, w]];
Rinv = Table[Module[{c = z5^(-j), ci = z5^(-i), pser, form},
     pser = pw[j];
     form = Series[(1/(s5 ci))/(ci - pser)^2 D[pser, w], {w, 0, 2 nn + 1}];
     Expand[-(1/th) Lap[form]]], {j, 0, 4}, {i, 0, 4}];
RinvT = Map[Expand[# /. th -> -1/t] &, Rinv, {2}];
diff = Table[trunc[RinvT[[j + 1, i + 1]]] - M0[[i + 1, j + 1]], {j, 0, 4}, {i, 0, 4}];
Print["max |R_TR^{-1}(1/th)^j_i - M0(t)_{ij}| (t^0..t^", nn, "): ",
  Max[Abs[Flatten[CoefficientList[#, t] & /@ Flatten[diff]]]]];
Print["transpose test max |R_TR^{-1}(1/th)^j_i - M0(t)_{ji}|: ",
  Max[Abs[Flatten[CoefficientList[#, t] & /@ Flatten[Table[trunc[RinvT[[j + 1, i + 1]]] - M0[[j + 1, i + 1]], {j, 0, 4}, {i, 0, 4}]]]]]];
Print["M0_{00} = ", Chop[N[M0[[1, 1]], 10]], " ; M0_{12} = ", Chop[N[M0[[2, 3]], 10]]];
