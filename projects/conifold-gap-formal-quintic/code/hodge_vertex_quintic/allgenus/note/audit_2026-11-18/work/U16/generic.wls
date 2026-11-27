(* Independent check of eq:fsp-generic and of lem:fsp-univ on the Stirling example *)
delta[s_] := (s^6 - s)/5;
KK = 4; ord = 2 KK + 4;
S2[s_] := -1/delta[s];
F3[t_] := Sum[(D[S2[s], {s, k - 2}] /. s -> L)/k! t^k, {k, 3, ord}];
r[s_] := 5 s/(s^5 - 1);
aaPoly = Normal[Series[Sqrt[r[L + tt]/r[L]], {tt, 0, ord}]]; aa[t_] := aaPoly /. tt -> t;
mom[n_, C_] := If[EvenQ[n], (n - 1)!! C^(n/2), 0];
SPn[amp_] := Module[{ser, w, tau, poly},
  ser = Normal[Series[Exp[F3[w tau]/w^2] amp[w tau, w^2], {w, 0, 2 KK}]];
  poly = Expand[ser];
  poly = poly /. tau^n_. :> mom[n, delta[L]];
  Simplify[Expand[poly /. w -> Sqrt[z]]]];
(* test 1: b = 1 *)
ampA[t_, z_] := aa[t];
ampS[t_, z_] := (L + t) aa[t];
N0 = SPn[ampA]; Ns = SPn[ampS];
DD[f_] := delta[L] D[f, L];
chk = Series[z DD[L N0] + L (L N0) - L Ns, {z, 0, KK - 1}];
Print["generic b=1 residual: ", Simplify[chk]];
(* test 2: b = 1 + z/sigma + z^2 sigma^2 *)
bb[t_, z_] := 1 + z/(L + t) + z^2 (L + t)^2;
ampA2[t_, z_] := aa[t] bb[t, z];
ampS2[t_, z_] := (L + t) aa[t] bb[t, z];
N2 = SPn[ampA2]; Ns2 = SPn[ampS2];
chk2 = Series[z DD[L N2] + L (L N2) - L Ns2, {z, 0, KK - 1}];
Print["generic b=1+z/s+z^2 s^2 residual: ", Simplify[chk2]];
Print["N0 to z^2: ", Simplify[Series[N0, {z, 0, 2}]]];
