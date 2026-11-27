(* ============================== 6. conifold frame ============================== *)
(* t_c(Delta) from the analytic solution of theta^2 + 3y(3theta+1)(3theta+2) at Delta = 1+27y = 0; XX and y as series in tcv *)
tcSetup[NC_] := Module[{L2b, wc, ceq, csol, wcDel},
   L2b[f_] := Module[{th = y D[#, y] &}, (th[th[f]] + 3 y (9 th[th[f]] + 9 th[f] + 2 f))];
   wc = Sum[a[n] Del^n, {n, 0, NC}] /. Del -> 1 + 27 y;
   ceq = CoefficientList[Normal[Series[L2b[wc] /. y -> (Del - 1)/27, {Del, 0, NC - 1}]], Del];
   csol = Solve[Thread[ceq == 0] /. a[0] -> 1, Table[a[n], {n, 1, NC}]];
   wcDel = Sum[a[n] Del^n, {n, 0, NC}] /. First[csol] /. a[0] -> 1;
   tc = Module[{ser = Normal[Series[wcDel/(Del - 1), {Del, 0, NC}]]}, Sum[Coefficient[ser, Del, n] Del^(n + 1)/(n + 1), {n, 0, NC}]];
   tc = Normal[Series[tc/SeriesCoefficient[tc, {Del, 0, 1}], {Del, 0, NC}]];
   invcf = InverseSeries[Series[tc, {Del, 0, NC}], tcv];
   C1cf = Normal[y 27 D[tc, Del]] /. Del -> 1 + 27 y;
   XXcf = Series[y D[C1cf, y]/C1cf /. y -> (Del - 1)/27 /. Del -> invcf, {tcv, 0, NC - 2}];
   ycf = (invcf - 1)/27;
   NCcur = NC];
(* polar part t_c^-(2g-2) .. t_c^-1 of a term in y, XX, eps1, eps2 with eps1 = x, eps2 = 1 *)
polarOfTerm[term_, gg_] := Module[{ser},
   ser = Series[term /. \[Epsilon][1] -> x /. \[Epsilon][2] -> 1 /. XX -> XXcf /. y -> ycf, {tcv, 0, -1}];
   Table[Factor[SeriesCoefficient[ser, -n]], {n, 2 gg - 2, 1, -1}]];
(* the known gap coefficient: (-3)^(g-1) (2g-3)! [t^(2g)] t^2/(4 sinh(e1 t/2) sinh(e2 t/2)) *)
gapCoefficient[gg_] := Module[{tv}, (-3)^(gg - 1) (2 gg - 3)! SeriesCoefficient[Series[tv^2/(4 Sinh[\[Epsilon][1] tv/2] Sinh[\[Epsilon][2] tv/2]), {tv, 0, 2 gg}], 2 gg]];
