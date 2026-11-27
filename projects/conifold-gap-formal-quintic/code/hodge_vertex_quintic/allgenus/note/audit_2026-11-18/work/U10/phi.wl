x[p_] := p^4/4 + v/p;
eps[k_, p_] := -D[x[p], {p, 2 k}]/(4^k (2 k + 1)!);
Dx[f_] := Simplify[D[f, p]/D[x[p], p]];
(* a_r = -[hbar^{2h}][w^{r+1}] e^E, E = sum_k w^{2k+1} hbar^{2k} eps_k *)
PhiA[h_] := Module[{E0, ser, tot = 0, ar, r, w, hb},
  E0 = Sum[w^(2 k + 1) hb^(2 k) eps[k, p], {k, 1, h}];
  ser = Normal[Series[Exp[E0], {hb, 0, 2 h}]];
  ser = Coefficient[Expand[ser], hb, 2 h];
  Do[ar = -Coefficient[ser, w, r + 1];
   If[ar =!= 0, tot += Nest[Dx, ar/D[x[p], p], r - 1]], {r, 1, 4 h}];
  Together[tot]];
Do[ph = PhiA[h];
 Print[h, ": at0 ", Series[ph, {p, 0, 1}], "  expected ", -BernoulliB[2 h]/(2 h) v^(1 - 2 h), "  atInf ", Limit[ph, p -> Infinity]], {h, 1, 4}];
