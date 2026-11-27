(* sanity: Phi^SP from def:con-PhiSP vs c_n(L)=L^{4n} ph_n(L^{-5}) ; and Stirling lemma *)
N0 = 3; Mu = 5;
bk[k_] := BernoulliB[2 k]/(2 k (2 k - 1));
W[X_] := Sum[u^(m - 1)/m (1/X + X^(5 m - 1)/(5 m - 1)), {m, 1, Mu + 1}];
Wser = Normal[Series[W[1 + t], {t, 0, 2 N0 + 2}]];
F3 = -Expand[Wser - (Wser /. t -> 0) - Coefficient[Wser, t, 2] t^2];
CC = (1 - u)/5;
logA = Normal[Series[-(1/2) Log[1 - u (1 + t)^5] + (1/2) Log[1 - u], {t, 0, 2 N0}]] +
   Sum[bk[k] (u (1 + t) zh)^(2 k - 1) (5^(1 - 2 k) - 5 Sum[Binomial[2 k - 2 + 5 P, 5 P] u^P (1 + t)^(5 P), {P, 0, Mu}]), {k, 1, N0}];
expo = Expand[(F3/zh + logA) /. {t -> e y, zh -> e^2}];
ser = Normal[Series[Exp[expo], {e, 0, 2 N0}]];
wick[expr_] := Module[{ex = Expand[expr]}, (ex /. y^p_. :> If[EvenQ[p], (p - 1)!! CC^(p/2), 0]) /. y -> 0];
ph = Table[Normal[Series[wick[Coefficient[ser, e, 2 n]], {u, 0, 3}]], {n, 1, N0}];
Print["ph^SP_n (to u^3): ", Expand[ph]];
Print["expected ph3: ", Expand[(11571 - 16128 u + 4288 u^2)/80000]];
(* Stirling lemma check: SPn_a[a log x - x; x^(M-1)] *)
st[M_, a_] := Module[{F3s, sers, CCs = a},
  F3s = Sum[(-1)^(l + 1) t^l/(l a^(l - 1)), {l, 3, 10}];
  sers = Normal[Series[Exp[Expand[(F3s/z) /. {t -> e y, z -> e^2}]] ((a + e y)^(M - 1)), {e, 0, 8}]];
  Table[(Expand[Coefficient[sers, e, 2 n]] /. y^p_. :> If[EvenQ[p], (p - 1)!! CCs^(p/2), 0]) /. y -> 0, {n, 0, 4}]];
pred[M_, a_] := Table[SeriesCoefficient[a^(M - 1) If[M >= 0, Product[1 + i z/a, {i, 0, M - 1}], Product[1/(1 - i z/a), {i, 1, -M}]] Exp[Sum[bk[k] (z/a)^(2 k - 1), {k, 1, 4}]], {z, 0, n}], {n, 0, 4}];
Do[Print["Stirling M=", M, ": ", st[M, 3/2] == pred[M, 3/2]], {M, {-2, 0, 1, 3}}];
