(* 1-D normalised stationary phase via t = h tau, z = h^2 *)
sp[F3_, amp_, C_, K_] := Module[{ser, poly, res},
  (* F3[t] cubic part, amp[t,z] amplitude, C covariance, K order in z *)
  ser = Series[Exp[F3[h tau]/h^2] amp[h tau, h^2], {h, 0, 2 K}];
  poly = Normal[ser];
  poly = Expand[poly];
  res = poly /. tau^n_. :> If[EvenQ[n], (n - 1)!! C^(n/2), 0];
  res = res /. tau -> 0;
  Expand[Normal[Series[res, {h, 0, 2 K}]] /. h -> Sqrt[z]]
];
(* Test 1: Stirling *)
K = 5;
bk[k_] := BernoulliB[2 k]/(2 k (2 k - 1));
poch[M_, w_] := If[M >= 0, Product[1 + i w, {i, 0, M - 1}], Product[1/(1 - i w), {i, 1, -M}]];
Do[
  F3s[t_] := a Log[1 + t/a] - t + t^2/(2 a);
  amps[t_, zz_] := (a + t)^(M - 1);
  lhs = sp[F3s, amps, a, K];
  rhs = Normal[Series[a^(M - 1) poch[M, z/a] Exp[Sum[bk[k] (z/a)^(2 k - 1), {k, 1, K}]], {z, 0, K}]];
  Print["Stirling M=", M, " diff: ", Simplify[lhs - rhs]],
  {M, {-2, -1, 0, 1, 2, 3}}];
