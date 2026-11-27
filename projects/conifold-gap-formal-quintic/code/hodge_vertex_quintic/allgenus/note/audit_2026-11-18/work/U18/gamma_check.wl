(* sanity check of Theorem Gamma (a): c(z;1)=h(z) to low order; and c_n(zeta L)=zeta^-n c_n(L) *)
NN = 2; K = 6*NN + 1;
pgen = Series[5 (L + t)^4/(1 - (L + t)^5), {t, 0, K}];
p[k_] := p[k] = Together[SeriesCoefficient[pgen, k - 1]];
a = Together[5/L + p[1]];
f = Sum[(5 (-1)^(k + 1) L^(-k) + p[k])/(k (k + 1)) t^(k + 1), {k, 2, K}]; (* cubic and higher *)
s1 = Sum[((-1)^(k + 1) L^(-k) + p[k])/(2 k) t^k, {k, 1, K}];
bk[k_] := BernoulliB[2 k]/(2 k (2 k - 1));
p1sig[sig_] := 5 sig^4/(1 - sig^5);
S2k[k_, sig_] := bk[k] ((5 sig)^(1 - 2 k) + Sum[(Exp[2 Pi I j/5] - sig)^(1 - 2 k), {j, 0, 4}]);
(* use generating function for sum over j of (zeta^j - sig)^(-m): coefficient extraction *)
psum[m_, sig_] := Module[{tt}, SeriesCoefficient[5 (sig + tt)^4/(1 - (sig + tt)^5), {tt, 0, m - 1}]];
S2kT[k_] := bk[k] (Normal[Series[(5 (L + t))^(1 - 2 k), {t, 0, K}]] + Normal[Series[psum[2 k - 1, L + t], {t, 0, K}]]);
logA = s1 + Sum[z^(2 k - 1) S2kT[k], {k, 1, NN}];
expo = f/z + logA;
(* expand exp in t up to degree K, z handled as Laurent *)
ser = Normal[Series[Exp[expo /. t -> eps t], {eps, 0, K}]] /. eps -> 1;
ser = Expand[ser];
rule = t^n_. :> If[OddQ[n], 0, (n - 1)!! (-z/a)^(n/2)];
moment[expr_] := Module[{e = Expand[expr]}, e /. {t^n_ :> If[OddQ[n], 0, (n - 1)!! (-z/a)^(n/2)], t -> 0}];
cs = moment[ser];
cser = Series[cs, {z, 0, NN}];
c0 = Together[SeriesCoefficient[cser, 0]];
Print["c0 = ", Simplify[c0]];
Do[cn = Together[SeriesCoefficient[cser, n]/c0]; cn = Factor[cn]; Print["c", n, "(L) = ", cn]; Print["  limit L->1: ", Limit[cn, L -> 1]], {n, 1, NN}];
hser = Series[Exp[Sum[bk[k] (5^(1 - 2 k) + Sum[(Exp[2 Pi I j/5] - 1)^(1 - 2 k), {j, 1, 4}]) z^(2 k - 1), {k, 1, NN}]], {z, 0, NN}];
Print["h = ", FullSimplify[Normal[hser]]];
