(* sanity check of cor:gau-cycsys (i),(ii): chain identity with zD+w, non-constant w *)
Dq[f_] := -(1 - u) D[f, u];
Dn[f_, 0] := f; Dn[f_, n_] := Dq[Dn[f, n - 1]];
A = {-24/625 (1 - u), -2/5 (1 - u), -7/5 (1 - u), -2 (1 - u), u};
PF4[f_] := Sum[A[[k + 1]] Dn[f, k], {k, 0, 4}];
NN = 16;
sol[k_] := Module[{c, f, eqs, s},
  f = u^k + Sum[c[n] u^n, {n, 3, NN + 3}];
  eqs = CoefficientList[Normal[Series[PF4[f], {u, 0, NN - 1}]], u];
  s = Solve[Thread[eqs == 0], Table[c[n], {n, 3, NN + 2}]][[1]];
  (f /. s) /. c[NN + 3] -> 0];
f0 = sol[0]; f1 = sol[1]; f2 = sol[2]; w1 = f1 + 7/10 f2;
y0 = f0 + 1/3 w1 - 2/7 f2; t = w1/y0; td = Dq[t];
kap = {td, 1/(u y0^2 td^2), td, y0, y0};
q = (1 - u)/5^5;
w = 2/3 + 5/7 u - 3/11 u^2 + u^3/13;   (* non-constant w *)
F = 1 + 3 u - 1/2 u^2 + 7/5 u^3 + 2 u^4 - u^5/9 + u^6/4 + 5 u^7;
op[f_] := z Dq[f] + w f;
(* left: kappa-chain applied to F/y0 *)
lhs = Fold[Normal[Series[op[#1]/#2, {u, 0, 9}]] &, Normal[Series[F/y0, {u, 0, 14}]], kap];
(* right: y0^{-1} (sigma(PF5)+1) F, PF5 = (zD)^5 - 1 - q prod (5 zD + k z) *)
p5 = Nest[op, F, 5];
pr = Fold[(5 op[#1] + #2 z #1) &, F, {1, 2, 3, 4, 5}];
rhs = (p5 - q pr)/y0;
Print["chain identity (ii), should vanish to the computed order: ", Simplify[Series[lhs - rhs, {u, 0, 4}]]];
(* also check (i): w=0 *)
w = 0;
lhs0 = Fold[Normal[Series[op[#1]/#2, {u, 0, 9}]] &, Normal[Series[F/y0, {u, 0, 14}]], kap];
p50 = Nest[op, F, 5]; pr0 = Fold[(5 op[#1] + #2 z #1) &, F, {1, 2, 3, 4, 5}];
Print["chain identity (i): ", Simplify[Series[lhs0 - (p50 - q pr0)/y0, {u, 0, 4}]]];
(* and PF5 = z^5 D PF4 - 1 *)
Print["PF5 = z^5 D PF4 - 1: ", Simplify[(p50 - F - q pr0) - (z^5 Dq[PF4[F]] - F)]];
