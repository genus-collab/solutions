bk[k_] := BernoulliB[2 k]/(2 k (2 k - 1));
f[x_] := 5 x^4/(1 - x^5);
P[r_, x_] := (D[f[y], {y, r - 1}]/(r - 1)!) /. y -> x;
S2k[k_, x_] := bk[k] ((5 x)^(1 - 2 k) + P[2 k - 1, x]);
Sp[x_] := lq + 5 Log[5 x] - Log[x^5 - 1];
NN = 7; KK = 4;
lhs = Sum[(D[Sp[x], {x, m - 1}] /. x -> sig) z^(m - 1)/m!, {m, 1, NN + 2}] +
   (1/2 Log[(sig + z)/sig] - 1/2 Log[(1 - (sig + z)^5)/(1 - sig^5)]) +
   Sum[z^(2 k - 1) (S2k[k, sig + z] - S2k[k, sig]), {k, 1, KK}];
rhs = lq + Sum[Log[5 sig + k z], {k, 1, 5}] - Log[(sig + z)^5 - 1];
d = Series[lhs - rhs, {z, 0, NN}];
Print[Simplify[Normal[d], sig > 2]];
(* also check Taylor data claims *)
p[k_, L_] := SeriesCoefficient[5 (L + t)^4/(1 - (L + t)^5), {t, 0, k - 1}];
S0pp[s_] := 5/s + Sum[1/(Exp[2 Pi I j/5] - s), {j, 0, 4}];
Print[Simplify[S0pp[s] - 5/(s (1 - s^5))]];
(* S0(L+t)-S0(L) coefficient check: S0'' = 5/s + sum 1/(zeta^j - s); m-th Taylor coeff *)
Do[Print[m, " ", Simplify[(D[5/s + 5 s^4/(1 - s^5), {s, m - 2}]/m! /. s -> L) - (5 (-1)^(m) L^(-(m - 1)) + p[m - 1, L])/(m (m - 1))]], {m, 2, 7}];
