(* independent check of Theorem P(a)(i),(iv) to low order: compute c_n(L) by Wick contraction *)
nmax = 6;
f[x_] := 5 x^4/(1 - x^5);
pw[m_, x_] := D[f[y], {y, m - 1}]/(m - 1)! /. y -> x;
bk[k_] := BernoulliB[2 k]/(2 k (2 k - 1));
S2k[k_, x_] := bk[k] ((5 x)^(1 - 2 k) + pw[2 k - 1, x]);
S0pp[x_] := 5/(x (1 - x^5));
tdeg = 2 nmax + 2;
(* Taylor of S0 at L: S0(L+t)-S0(L) = int_0^t (t-u) S0''(L+u) du *)
s0ser = Normal[Series[S0pp[L + u], {u, 0, tdeg}]];
S0t = Expand[Integrate[(t - u) s0ser, {u, 0, t}]];
H = Coefficient[S0t, t, 2]*2;  (* S0''(L) *)
Cc = -1/H;
F3 = S0t - H t^2/2;
amp = Exp[(1/2) Log[((L + t)/L) (1 - L^5)/(1 - (L + t)^5)] + Sum[z^(2 k - 1) S2k[k, L + t], {k, 1, nmax}]];
(* substitute t = w tau, z = w^2 *)
expr = Exp[(F3 /. t -> w tau)/w^2] (amp /. {t -> w tau, z -> w^2});
ser = Normal[Series[expr, {w, 0, 2 nmax}]];
ser = Expand[ser];
mom[n_] := If[OddQ[n], 0, (n - 1)!! Cc^(n/2)];
res = Sum[Coefficient[ser, tau, n] mom[n], {n, 0, Exponent[ser, tau]}] /. tau -> 0;
res = Expand[res];
cs = Table[Together[Simplify[Coefficient[res, w, 2 n]]], {n, 0, nmax}];
cs = Together[cs/cs[[1]]];
Print["c_n at L=1: ", cs /. L -> 1];
Print["c_n as Laurent in L: ", Table[Factor[cs[[n + 1]]], {n, 0, nmax}]];
(* PF5^L (L c) mod z^(nmax+1) *)
cz = Sum[cs[[n + 1]] z^n, {n, 0, nmax}];
DD[g_] := D[g, L] (L^6 - L)/5;
op[g_] := z DD[g] + L g;
q = (1 - L^-5)/5^5;
lhs = Nest[op, L cz, 5] - L cz - q Fold[5 op[#1] + #2 z #1 &, L cz, Range[5]];
lhs = Series[Together[lhs], {z, 0, nmax}];
Print["PF5^L(Lc) mod z^", nmax + 1, ": ", Simplify[Normal[lhs]]];
