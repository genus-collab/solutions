(* genus-1 two-point pairings, test of (V1) at (1,1)->(1,2) *)
ord = 9; prec = 40;
setup[v_] := Module[{},
  sig = Table[N[Exp[-2 Pi I b/5] v^(1/5), prec], {b, 0, 4}];
  x[p_] := p^4/4 + v/p;
  invol = Table[Module[{c = sig[[bb]], sAns, eqs, sol},
     sAns = -t + Sum[a[k] t^k, {k, 2, ord}];
     eqs = CoefficientList[Normal[Series[x[c + sAns] - x[c + t], {t, 0, ord + 1}]], t];
     sol = Quiet[Solve[Thread[Drop[eqs, 3] == 0], Table[a[k], {k, 2, ord}]][[1]]];
     c + (sAns /. sol)], {bb, 1, 5}];
  (* omega_{1,1}(y) dy coefficient: sum_c Res_t K_c(y,q) w02(q,iq) *)
  w11[y_] = Sum[Module[{c = sig[[bb]], iq = invol[[bb]], q, diq, xp},
      q = c + t; diq = D[iq, t]; xp = D[x[p], p] /. p -> q;
      SeriesCoefficient[Series[(1/2) (1/(y - q) - 1/(y - iq)) diq/((q - iq)^3 xp), {t, 0, 2}], -1]],
     {bb, 1, 5}];
  ];
pair12[f0_, f1_] := Module[{tot = 0, d, A1},
  d[f_, p0_] := (D[f[p], p] /. p -> p0);
  A1 = Table[d[f1, sig[[b]]]/(5 sig[[b]]^2), {b, 1, 5}];
  Do[Module[{c = sig[[bb]], iq = invol[[bb]], q, diq, xp, br},
    q = c + t; diq = D[iq, t]; xp = D[x[p], p] /. p -> q;
    br = -Sum[1/(q - sig[[b]])^2 diq/(iq - sig[[b]])^2 A1[[b]], {b, 1, 5}] +
      d[f1, q] w11[iq] diq + w11[q] d[f1, iq] diq;
    tot += SeriesCoefficient[Series[(1/2) (f0[q] - f0[iq])/((q - iq) xp) br, {t, 0, 2}], -1]],
   {bb, 1, 5}];
  Chop[tot, 10^-20]];
lg = (-Log[#]) &;
pair1[f0_] := Chop[Sum[Module[{c = sig[[bb]], iq = invol[[bb]], q, diq, xp},
     q = c + t; diq = D[iq, t]; xp = D[x[p], p] /. p -> q;
     SeriesCoefficient[Series[(1/2) (f0[q] - f0[iq]) diq/((q - iq)^3 xp), {t, 0, 2}], -1]], {bb, 1, 5}], 10^-20];
Do[setup[vv];
 Print["v=", vv, ": <-log p>_1=", pair1[lg], " (expect ", N[-1/(12 vv)], ")"];
 Print["   <(-log p)^2>_1=", pair12[lg, lg], " (expect ", N[1/(12 vv^2)], ")"];
 Print["   <p^5>_1=", pair1[#^5 &], "  <p^10>_1=", pair1[#^10 &]];
 Print["   <p^5,-log p>_1=", pair12[#^5 &, lg], " (expect 0)"];
 Print["   <p^10,-log p>_1=", pair12[lg, #^10 &], " (V1 predicts ", pair1[#^10 &]/vv + 10 pair1[#^5 &], ")"],
 {vv, {13/10, 7/10, 2}}];
