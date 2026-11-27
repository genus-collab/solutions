(* genus-0 four-point pairings by the recursion, to test thm:mod-var (0,3)->(0,4) via thm:mod-V *)
ord = 8;
pair4[v_, legs_] := Module[{sig, tot = 0, c, x, sAns, eqs, sol, iq, q, diq, xp, T, d, res, f0, others, splits, br},
  sig = Table[N[Exp[-2 Pi I b/5] v^(1/5), 40], {b, 0, 4}];
  x[p_] := p^4/4 + v/p;
  d[f_, p0_] := (D[f[p], p] /. p -> p0);
  (* T_{jk}(q') coefficient of dq' *)
  T[fj_, fk_, qq_] := -Sum[d[fj, sig[[b]]] d[fk, sig[[b]]]/((qq - sig[[b]])^2 5 sig[[b]]^2), {b, 1, 5}];
  f0 = legs[[1]]; others = Rest[legs];
  Do[
   c = sig[[bb]];
   sAns = -t + Sum[a[k] t^k, {k, 2, ord}];
   eqs = CoefficientList[Normal[Series[x[c + sAns] - x[c + t], {t, 0, ord + 1}]], t];
   sol = Solve[Thread[Drop[eqs, 3] == 0], Table[a[k], {k, 2, ord}]][[1]];
   iq = c + (sAns /. sol); q = c + t; diq = D[iq, t];
   xp = D[x[p], p] /. p -> q;
   br = Sum[
     With[{i = ii, jk = Delete[Range[3], ii]},
      d[others[[i]], q] T[others[[jk[[1]]]], others[[jk[[2]]]], iq] diq +
       T[others[[jk[[1]]]], others[[jk[[2]]]], q] d[others[[i]], iq] diq], {ii, 1, 3}];
   res = SeriesCoefficient[Series[(1/2) (f0[q] - f0[iq])/((q - iq) xp) br, {t, 0, 2}], -1];
   tot += res, {bb, 1, 5}];
  Chop[tot, 10^-20]];
lg = (-Log[#]) &;
Print["<p^6 p^4 (-log p)^2>_0 at v=1.7 (expect 48): ", pair4[17/10, {#^6 &, #^4 &, lg, lg}]];
Print["<p^6 p^4 (-log p)^2>_0 slots permuted: ", pair4[17/10, {lg, #^4 &, lg, #^6 &}]];
Print["<p^7 p^3 (-log p)^2>_0 at v=0.9: ", pair4[9/10, {#^7 &, #^3 &, lg, lg}], " expect ", N[21 + 7*3 ]];
