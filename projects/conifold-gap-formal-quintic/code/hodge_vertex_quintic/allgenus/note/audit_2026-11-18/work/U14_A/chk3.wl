(* U14_A sanity checks: lem:con-cfu, lem:con-red, lem:con-S2, lem:con-id (b), thm:con-Pa (ii)/(iv), residue identity V = Sum pi_m s^-m IC_m *)
U = 12;   (* u-truncation *)
NZ = 4;  (* zh-truncation for PhiSP *)
trunc[expr_, n_] := Normal[Series[expr /. u -> eps*u, {eps, 0, n}]] /. eps -> 1;
bk[k_] := BernoulliB[2 k]/(2 k (2 k - 1));

(* ---------- lem:con-cfu ---------- *)
ys = 1 + Sum[u^k Sum[a[k, d] T^d, {d, 0, k}], {k, 1, U}];
Vs = 1 + Sum[vv[k] u^k, {k, 1, U}];
rhs = Sum[u^(m - 1)/m (1 + T^m/(5 m - 1)), {m, 1, U + 1}];
eq = trunc[Vs/ys + T ys^4/4 - rhs, U];
sol = {};
Do[
  ck = Expand[Coefficient[eq /. sol, u, k]];
  vars = Join[Table[a[k, d], {d, 0, k}], {vv[k]}];
  s1 = Solve[Thread[CoefficientList[ck, T] == 0], vars];
  If[Length[s1] != 1, Print["cfu: no/multiple solution at k=", k, " ", s1]];
  sol = Join[sol, s1[[1]]], {k, 1, U}];
Ysol = Expand[ys /. sol]; Vsol = Expand[Vs /. sol];
Print["Y_1(T) = ", Coefficient[Ysol, u, 1], "   V_1 = ", Coefficient[Vsol, u, 1]];
tau = 1/Sum[4 u^n/(5 n + 4), {n, 0, U + 1}];
Print["V - tau^(-5/4) mod u^(U+1): ", Simplify[trunc[Vsol - tau^(-5/4), U]]];
Print["deg Y_k <= k: ", And @@ Table[Exponent[Coefficient[Ysol, u, k], T] <= k, {k, 1, U}]];
YofX = X (Ysol /. T -> X^5);
What = Sum[u^(m - 1)/m (1/X + X^(5 m - 1)/(5 m - 1)), {m, 1, U + 1}];
Print["What - W_V(Y(X)) mod u^(U+1): ", Simplify[trunc[What - (Vsol/YofX + YofX^4/4), U]]];
(* inverse map X(Y) = Y (1 + sum u^k Xk(Y^5)) *)
xs = 1 + Sum[u^k Sum[b[k, d] T^d, {d, 0, k + 1}], {k, 1, U}];
XofY = Y (xs /. T -> Y^5);
eqi = trunc[(XofY /. Y -> YofX) - X, U];
soli = {};
Do[ck = Expand[Coefficient[eqi /. soli, u, k]];
  vars = Table[b[k, d], {d, 0, k + 1}];
  s1 = Solve[Thread[CoefficientList[ck, X] == 0], vars];
  If[Length[s1] != 1, Print["inverse: problem at k=", k]];
  soli = Join[soli, s1[[1]]], {k, 1, U}];
XofY = Expand[XofY /. soli];
Print["deg X_k <= k: ", And @@ Table[Exponent[Coefficient[xs /. soli, u, k], T] <= k, {k, 1, U}]];
Print["Y(1) - tau^(-1/4): ", Simplify[trunc[(YofX /. X -> 1) - tau^(-1/4), U]]];

(* ---------- hatted amplitude ---------- *)
loga = 1/2 Sum[(u X^5)^n/n, {n, 1, U + 1}] +
   Sum[bk[k] (u X zh)^(2 k - 1) (5^(1 - 2 k) - 5 Sum[Binomial[2 k - 2 + 5 P, 5 P] u^P X^(5 P), {P, 0, U}]), {k, 1, U}];
ahat = Expand[trunc[Exp[loga], U]];

(* ---------- lem:con-red ---------- *)
fY = Expand[trunc[(ahat /. X -> XofY) D[XofY, Y], U]];
reduce1[poly_] := Module[{p = Expand[poly], n, c, done = False},
  While[! done,
   n = Exponent[p, Y];
   If[n < 4, done = True,
    c = Coefficient[p, Y, n];
    p = Expand[trunc[p - c Y^n + c ((n - 3) zh Y^(n - 4) + Vsol Y^(n - 5)), U]]]];
  p];
red = reduce1[fY];
pis = Table[m -> Expand[Coefficient[red, Y, m]], {m, -1, 3}];
Print["exponents left after reduction: ", Union[Exponent[#, Y] & /@ (List @@ Expand[red])], " min ", Exponent[red, Y, Min]];
(* grading: pi_m contains zh^k only with k = m mod 5 *)
grad = And @@ Flatten[Table[
     Module[{p = m /. pis, cl}, cl = CoefficientList[p, zh];
      Table[If[cl[[kk + 1]] =!= 0, Mod[kk - m, 5] == 0, True], {kk, 0, Length[cl] - 1}]], {m, -1, 3}]];
Print["grading of pi_m: ", grad];

bk[k_] := BernoulliB[2 k]/(2 k (2 k - 1));
(* --- residue identity --- *)
Emax = 11;
trunc[expr_, n_] := Normal[Series[expr /. u -> eps*u, {eps, 0, n}]] /. eps -> 1;
Xi = Sum[u^n/n, {n, 1, Emax}];
Gx = Sum[x^(5 m - 1)/(m (5 m - 1)), {m, 1, Emax}];
logafl = 1/2 Sum[x^(5 n)/n, {n, 1, Emax}] + Sum[bk[k] (x z)^(2 k - 1) (5^(1 - 2 k) - 5 Sum[Binomial[2 k - 2 + 5 P, 5 P] x^(5 P), {P, 0, Emax}]), {k, 1, Emax}];
eGa = Normal[Series[Exp[-Gx/z + logafl], {x, 0, Emax}]];
Vflat = Expand[trunc[Sum[(-Xi/z)^k/k! Coefficient[eGa, x, k - 1], {k, 1, Emax}], Emax]];
(* Sigma *)
Vpow[e_] := Vpow[e] = Expand[trunc[Vsol^e, Emax + 2]];
Sig = 0;
Do[
  pm = Expand[(m /. pis)];
  terms = If[pm === 0, {}, CoefficientRules[pm, {u, zh}]];
  Do[{i, j} = tr[[1]]; c = tr[[2]];
   ip = 0;
   While[True,
    e = m + 4 ip + 1;
    beta = j - 5 ip - m - 1;
    E0 = i - (4 beta - 1)/5;
    If[E0 > Emax, Break[]];
    If[e >= 0,
     Do[Sig += c (-1/4)^ip/ip! (-1)^e/e! Coefficient[Vpow[e], u, l] uu^(E0 + l) z^beta, {l, 0, Emax - E0}]];
    ip++], {tr, terms}], {m, -1, 3}];
Sig = Expand[Sig];
(* valid range: target (E,beta) needs pi to order i <= E + (4 beta - 1)/5 <= U *)
K = 80;
rules[expr_] := Module[{p = Expand[expr z^K]}, If[p === 0, {}, ({#[[1, 1]], #[[1, 2]] - K} -> #[[2]]) & /@ CoefficientRules[p, {uu, z}]]];
rS = rules[Sig]; rV = rules[Vflat /. u -> uu];
valid[{E_, b_}] := (E + (4 b - 1)/5 <= U) && E <= Emax;
keys = Union[First /@ Join[rS, rV]];
keysv = Select[keys, valid];
getc[r_, k_] := Module[{s = Select[r, #[[1]] == k &]}, If[s === {}, 0, s[[1, 2]]]];
bad = Select[keysv, getc[rS, #] =!= getc[rV, #] &];
Print["monomials compared: ", Length[keysv], "; of which nonzero in Vflat: ", Length[Select[keysv, getc[rV, #] =!= 0 &]], "; with beta >= 0: ", Length[Select[keysv, #[[2]] >= 0 &]]];
Print["mismatches: ", Length[bad], " ", Take[bad, UpTo[5]]];
Print["[z^-1] Vflat: ", Sum[getc[rV, {E, -1}] uu^E, {E, 0, Emax}]];
