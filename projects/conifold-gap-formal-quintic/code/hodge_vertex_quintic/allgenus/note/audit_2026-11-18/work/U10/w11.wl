v = 1; x[p_] := p^4/4 + v/p;
ord = 7;
inv[a_, t_] := Module[{f, g, gi, u},
  f = Normal[Series[x[a + t] - x[a], {t, 0, ord + 2}]];
  g = Series[Sqrt[Expand[f/t^2]] t, {t, 0, ord + 1}];
  gi = InverseSeries[g, u];
  Normal[ComposeSeries[gi, -g]]];
w11[p0_] := Module[{tot = 0, a, s, q, expr, t},
  Do[a = Exp[-2 Pi I b/5];
   s = a + inv[a, t];
   q = a + t;
   expr = (1/2) (1/(p0 - q) - 1/(p0 - s)) D[s, t]/((q - s)^3 (D[x[y], y] /. y -> q));
   tot += SeriesCoefficient[Series[expr, {t, 0, 1}], -1], {b, 0, 4}];
  N[tot, 20]];
PhiA1[p_] := Module[{y}, (1/D[x[y], y]) D[D[x[y], {y, 2}]/(24 D[x[y], y]), y] /. y -> p];
dPhi[p0_] := Module[{y}, D[PhiA1[y], y] /. y -> p0];
Print[Table[{p0, w11[p0], N[-dPhi[p0], 20]}, {p0, {37/100 + 21 I/100, 23/10}}]];
