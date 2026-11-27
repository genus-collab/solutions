sp[F3_, amp_, C_, K_] := Module[{ser, poly, res},
  ser = Series[Exp[F3[h tau]/h^2] amp[h tau, h^2], {h, 0, 2 K}];
  poly = Expand[Normal[ser]];
  res = poly /. tau^n_. :> If[EvenQ[n], (n - 1)!! C^(n/2), 0];
  res = res /. tau -> 0;
  Expand[Normal[Series[res, {h, 0, 2 K}]] /. h -> Sqrt[z]]
];
K = 3;
delta[s_] := (s^6 - s)/5;
r[s_] := 5 s/(s^5 - 1);
Spp[s_] := -1/delta[s];
(* cubic part of S at L: integrate twice the Taylor series of S'' minus its constant term *)
ser2 = Normal[Series[Spp[L + t], {t, 0, 6 K + 2}]];
F3expr = Integrate[Integrate[ser2 - (ser2 /. t -> 0), {t, 0, t1}] /. t1 -> t, {t, 0, t2}] /. t2 -> t;
F3expr = Expand[F3expr];
F3g[tt_] := F3expr /. t -> tt;
CC = delta[L];
Do[
  ampN[tt_, zz_] := Sqrt[r[L + tt]/r[L]] b[L + tt, zz];
  ampNs[tt_, zz_] := (L + tt) Sqrt[r[L + tt]/r[L]] b[L + tt, zz];
  N0 = sp[F3g, ampN, CC, K];
  Ns = sp[F3g, ampNs, CC, K];
  lhs = z delta[L] D[L N0, L] + L (L N0);
  diff = Normal[Series[lhs - L Ns, {z, 0, K}]];
  Print["b=", b[s, zz], "  diff: ", Simplify[diff]],
  {b, {Function[{s, zz}, 1], Function[{s, zz}, s^2 + zz/s], Function[{s, zz}, 1/(s^5 - 1) + zz^2 s]}}];
(* sanity: wrong sign variant *)
b = Function[{s, zz}, 1];
ampN[tt_, zz_] := Sqrt[r[L + tt]/r[L]];
ampNs[tt_, zz_] := (L + tt) Sqrt[r[L + tt]/r[L]];
N0 = sp[F3g, ampN, CC, K]; Ns = sp[F3g, ampNs, CC, K];
Print["wrong sign variant (-zD+L): ", Simplify[Normal[Series[-z delta[L] D[L N0, L] + L^2 N0 - L Ns, {z, 0, 2}]]]];
