(* independent test of Theorem LMI (eq:gam-lmi) on generic data *)
N0 = 3; ss = 0; as = 3; M = 2;
phi[s_] := s/5 + s^2/7 + 2 s^3/3 - s^4/11;
FF[s_, z_] := (1 + 2 s) Exp[z (1 + s^2) + z^2 s];
bk[k_] := BernoulliB[2 k]/(2 k (2 k - 1));
a[s_] := as (s - ss);
Poch[m_, w_] := If[m >= 0, Product[1 + i w, {i, 0, m - 1}], Product[1/(1 - i w), {i, 1, -m}]];
(* phase S'' and higher derivs at L: phi^(k) + (a log a - a)^(k) *)
Spp[s_] := phi''[s] + as^2/a[s];
CC = Together[-1/Spp[L]];
F3 = Sum[Together[(D[phi[x], {x, j}] /. x -> L)/j! + (-1)^j as/(j (j - 1) (L - ss)^(j - 1))] t^j, {j, 3, 2 N0 + 2}];
ampl = (1 + t/(L - ss))^(M - 1/2) Poch[M, z/a[L + t]] Exp[Sum[bk[k] (z/a[L + t])^(2 k - 1), {k, 1, N0}]] FF[L + t, z]/FF[L, 0];
expo = Expand[F3/z] /. {t -> e y, z -> e^2};
amp2 = ampl /. {t -> e y, z -> e^2};
ser = Normal[Series[Exp[expo] amp2, {e, 0, 2 N0}]];
wick[expr_] := Module[{ex = Expand[expr]}, (ex /. y^p_. :> If[EvenQ[p], (p - 1)!! CC^(p/2), 0]) /. y -> 0];
res = Table[Together[wick[Coefficient[ser, e, 2 n]]], {n, 1, N0}];
lim = Table[Simplify[Limit[res[[n]], L -> ss]], {n, 1, N0}];
Print["limits L->s*: ", lim];
del = -M z/as;
phit[x_] := phi[ss + x] - phi[ss] - phi'[ss] x;
pred = Series[Exp[phit[del]/z] FF[ss + del, z]/FF[ss, 0], {z, 0, N0}];
Print["prediction: ", Table[SeriesCoefficient[pred, n], {n, 1, N0}]];
Print["c_1(L) sample: ", Simplify[res[[1]]]];
