delta[s_] := (s^6 - s)/5;
KK = 3; ord = 2 KK + 6;
S2[s_] := -1/delta[s];
Scoef[k_] := (D[S2[s], {s, k - 2}] /. s -> L)/k!;
SL[t_] := Sum[Scoef[k] t^k, {k, 2, ord + 4}];   (* S(L+t)-S(L) *)
F3[t_] := Sum[Scoef[k] t^k, {k, 3, ord}];
mom[n_, C_] := If[EvenQ[n], (n - 1)!! C^(n/2), 0];
SPn[ampfun_] := Module[{ser, w, tau, poly},
  ser = Normal[Series[Exp[F3[w tau]/w^2] ampfun[w tau, w^2], {w, 0, 2 KK}]];
  poly = Expand[ser] /. tau^n_. :> mom[n, delta[L]];
  Simplify[Expand[poly /. w -> Sqrt[z]]]];
(* negative control: amplitude 1 instead of sqrt(r/r(L)) in generic identity *)
DD[f_] := delta[L] D[f, L];
one[t_, z_] := 1; sig[t_, z_] := L + t;
N0 = SPn[one]; Ns = SPn[sig];
Print["negative control residual (should be nonzero): ", Simplify[Series[z DD[L N0] + L (L N0) - L Ns, {z, 0, 2}]]];
(* transl (b) with p = z, amplitude a(sigma)=sigma^2: b_z(t) = exp((S(L+t+z)-S(L+t))/z) a(L+t+z) *)
bz[t_, z_] := Exp[(SL[t + z] - SL[t])/z] (L + t + z)^2;
aplain[t_, z_] := (L + t)^2;
lhs = SPn[bz]; rhs = SPn[aplain];
Print["transl(b) residual: ", Simplify[Series[lhs - rhs, {z, 0, KK - 1}]]];
