(* independent check: c_n(L) from def:gam-c via Wick, n<=N0 *)
N0 = 3;
S0pp[s_] := 5/(s (1 - s^5));
P1[s_] := 5 s^4/(1 - s^5);
psum[r_, s_] := (D[P1[x], {x, r - 1}]/(r - 1)!) /. x -> s;
bk[k_] := BernoulliB[2 k]/(2 k (2 k - 1));
S2k[k_, s_] := bk[k] ((5 s)^(1 - 2 k) + psum[2 k - 1, s]);
S1p[s_] := 1/(2 s) + P1[s]/2;
CC = Together[-1/S0pp[L]];
ser0 = Series[S0pp[L + t], {t, 0, 2 N0}];
F3 = Sum[Together[SeriesCoefficient[ser0, j - 2]] t^j/(j (j - 1)), {j, 3, 2 N0 + 2}];
ser1 = Series[S1p[L + t], {t, 0, 2 N0}];
logA = Sum[Together[SeriesCoefficient[ser1, k - 1]] t^k/k, {k, 1, 2 N0}] +
   Sum[z^(2 k - 1) Normal[Series[S2k[k, L + t], {t, 0, 2 N0}]], {k, 1, N0}];
expo = (F3/z + logA) /. {t -> e y, z -> e^2};
ser = Normal[Series[Exp[expo], {e, 0, 2 N0}]];
wick[expr_] := Module[{ex = Expand[expr]}, ex /. y^p_. :> If[EvenQ[p], (p - 1)!! CC^(p/2), 0]];
res = Table[Factor[Together[wick[Coefficient[ser, e, 2 n]] /. y -> 0]], {n, 1, N0}];
Print["c_n(L): ", res];
Print["at L=1: ", Table[Limit[res[[n]], L -> 1], {n, 1, N0}]];
h = Exp[-Sum[bk[k] z^(2 k - 1) (Sum[(1 - Exp[2 Pi I j/5])^(1 - 2 k), {j, 1, 4}] + (-5)^(1 - 2 k)), {k, 1, N0}]];
Print["h_n: ", Table[FullSimplify[SeriesCoefficient[h, {z, 0, n}]], {n, 1, N0}]];
