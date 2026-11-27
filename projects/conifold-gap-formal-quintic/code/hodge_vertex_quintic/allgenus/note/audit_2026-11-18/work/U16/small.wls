(* univ: Stirling example z^1 coefficient; subring counterexample; pg(i) 1D; homog 1D *)
mom[n_, C_] := If[EvenQ[n], (n - 1)!! C^(n/2), 0];
SPgen[F3_, amp_, C_, K_] := Module[{ser, w, tau},
  ser = Normal[Series[Exp[F3[w tau]/w^2] amp[w tau, w^2], {w, 0, 2 K}]];
  Expand[Expand[ser] /. tau^n_. :> mom[n, C] /. w -> Sqrt[z]]];
(* Stirling: F = log(1+t)-t, C = 1 *)
f3st[t_] := Normal[Series[Log[1 + t] - t + t^2/2, {t, 0, 12}]];
Print["Stirling SP: ", SPgen[f3st, (1 &), 1, 3], "   expected ", Series[Exp[z/12 - z^3/360], {z, 0, 3}]];
(* subring counterexample: F = -t^2/2 + t^3, amplitude 1, C = 1 *)
Print["Z-data SP: ", SPgen[(#^3 &), (1 &), 1, 2]];
(* pg(i): <exp(k t^2/2) t^2>_1 vs (1-k)^(-3/2) *)
Print["pg(i): ", Series[Sum[k^j/(2^j j!) mom[2 j + 2, 1], {j, 0, 6}] - (1 - k)^(-3/2), {k, 0, 6}]];
(* homog: SP[mu F; a(z)](z) vs SP[F; a(mu z')](z') at z' = z/mu, F3 = t^3 + t^4, a = 1 + z t, C = 1 for F, C/mu for mu F *)
lhs = SPgen[(mu (#^3 + #^4) &), (1 + #2 #1 &), 1/mu, 3];
rhs = SPgen[((#^3 + #^4) &), (1 + mu #2 #1 &), 1, 3] /. z -> z/mu;
Print["homog residual: ", Simplify[Series[lhs - rhs, {z, 0, 3}]]];
