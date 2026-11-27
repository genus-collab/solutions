(* Independent check: omega_{1,1} by direct TR with the kernel K_b of conv:mod-kernel,
   for x = p^4/4 + v/p, y = p, versus the ABDKS (7.10) prediction with omega = -omega^A. *)
X[p_] := p^4/4 + v/p;
nmax = 6;
sig = -s + Sum[c[k] s^k, {k, 2, nmax}];
eq = Series[(X[a + sig] - X[a + s]) /. v -> a^5, {s, 0, nmax + 1}];
coefs = Table[SeriesCoefficient[eq, k], {k, 3, nmax + 1}];
sol = {};
Do[ sol = Join[sol, Solve[(coefs[[k - 2]] /. sol) == 0, c[k - 1]][[1]]], {k, 3, nmax + 1}];
sigS = sig /. sol;
integrand = (1/2) (1/(p0 - a - s) - 1/(p0 - a - sigS)) D[sigS, s]/((s - sigS)^3 (D[X[q], q] /. q -> a + s /. v -> a^5));
res = SeriesCoefficient[Series[integrand, {s, 0, 1}], -1] // Together;
w11 = RootSum[(#^5 - v) &, Function[aa, Evaluate[res /. a -> aa]]] // Together;
target = -(1/24) D[(-3 p0^10 - 24 v p0^5 + 2 v^2)/(p0^5 - v)^3, p0] // Together;
Print["direct TR omega_11/dp: ", Factor[w11]];
Print["ABDKS (7.10) with omega=-omega^A: ", Factor[target]];
Print["difference: ", Simplify[w11 - target]];
