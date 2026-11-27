(* check DNOPS (2-33) index placement for the model curve, v=1, kappa=1 *)
xf[p_] := p^4/4 + 1/p;
z = Exp[2 Pi I/5]; s5 = I Sqrt[5];
ord = 6;
pofw[c_] := Module[{ww, ser, inv},
  ser = Series[s5 c e Sqrt[Normal[Series[-2 (xf[c + e] - xf[c])/(-5 c^2 e^2), {e, 0, ord}]]], {e, 0, ord}];
  c + Normal[InverseSeries[ser, w]]];
c0 = 1; c1 = z^(-1);
p1 = pofw[c0] /. w -> w1; p2 = pofw[c1] /. w -> w2;
(* sanity: x(p(w)) - x(c) = -w^2/2 *)
Print["check local coord: ", N[Normal[Series[xf[p1] - xf[c0] + w1^2/2, {w1, 0, 5}]], 10] // Chop];
Bf = D[p1, w1] D[p2, w2]/(p1 - p2)^2;
ser = Normal[Series[Bf, {w1, 0, 3}, {w2, 0, 3}]];
b[m_, n_] := N[Coefficient[Coefficient[ser, w1, m], w2, n] /. {w1 -> 0, w2 -> 0}, 20];
Print["b00 = ", Chop[b[0, 0]]];
Print["L L order psi1: b20 = ", Chop[b[2, 0]], "   order psi2: b02 = ", Chop[b[0, 2]]];
