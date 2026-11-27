(* toy test of eq:gam-lmi: s*=0, a*=1, M*=2, phi=sig^3, F=1+sig z + sig^2 z^0 ... *)
NN = 3; Tmax = 3 NN + 3;
bk[k_] := BernoulliB[2 k]/(2 k (2 k - 1));
M = 2; as = 1;
phi[x_] := x^3 + x^4/2;
FF[x_, z_] := Exp[x] (1 + x z + x^2 z^2);
Sfun[x_] := phi[x] + as x Log[as x] - as x;  (* lambda sigma dropped: only order >=2 matters *)
Afun[x_, z_] := (as x)^(M - 1/2) (1 + z/(as x)) Exp[Sum[bk[k] (z/(as x))^(2 k - 1), {k, 1, NN + 1}]] FF[x, z];
S2 = D[Sfun[x], {x, 2}] /. x -> L;
S3 = Normal[Series[Sfun[L + t], {t, 0, Tmax}]];
cub = Expand[S3 - (S3 /. t -> 0) - Coefficient[S3, t, 1] t - Coefficient[S3, t, 2] t^2];
amp = Normal[Series[Afun[L + t, z]/Afun[L, 0] /. z -> zz, {t, 0, Tmax}]];
expo = (cub /. t -> r tau)/r^2;
full = Normal[Series[Exp[expo] (amp /. {t -> r tau, zz -> r^2}), {r, 0, 2 NN}]];
CC = -1/S2;
wick[e0_] := Module[{e = Expand[e0]}, (e /. tau -> 0) + Sum[Coefficient[e, tau, 2 j] (2 j - 1)!! CC^j, {j, 1, 3 NN + 3}]];
res = Table[Simplify[wick[Coefficient[full, r, 2 n]]], {n, 0, NN}];
Print["limits L->0: ", Table[Limit[res[[n + 1]], L -> 0], {n, 0, NN}]];
delta = -M z/as;
pred = Series[Exp[(phi[delta] - phi[0] - phi'[0] delta)/z] FF[delta, z]/FF[0, 0], {z, 0, NN}];
Print["predicted: ", Table[SeriesCoefficient[pred, n], {n, 0, NN}]];
