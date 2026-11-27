(* U14_A: Pa (ii)/(iv) at low order (c_n(L) from def:gam-c Taylor data vs L^{4n} ph_n(L^-5)),
   and the residue identity V = Sum_m pi_m s^{-m} IC_m(v(u),z) (con-relab + fsp-res + con-red + model residues) *)
Get["state.mx"];
bk[k_] := BernoulliB[2 k]/(2 k (2 k - 1));
Print["ph^SP: ", phs];
(* --- c_n(L) from def:gam-c --- *)
nmax = 3; ord = 2 nmax + 4;
s0pp = Normal[Series[5/((L + tt) (1 - (L + tt)^5)), {tt, 0, ord}]];
S0t = Integrate[Integrate[s0pp, {tt, 0, tt1}] /. tt1 -> tt, {tt, 0, tt2}] /. tt2 -> tt;
S0t = Expand[S0t];
(* power sums sum_j (sigma - zeta^j)^(-r) *)
psum[r_, sig_] := (-1)^(r - 1)/(r - 1)! D[5 sg^4/(sg^5 - 1), {sg, r - 1}] /. sg -> sig;
S2k[k_, sig_] := bk[k] ((5 sig)^(1 - 2 k) + (-1)^(1 - 2 k) psum[2 k - 1, sig]);
rr[sig_] := 5 sig/(1 - sig^5);
ampSt = Sqrt[rr[L + tt]/rr[L]] Exp[Sum[zz^(2 k - 1) S2k[k, L + tt], {k, 1, nmax + 1}]];
ampt = Normal[Series[ampSt, {tt, 0, ord}, {zz, 0, nmax}]];
cov = -1/(5/(L (1 - L^5)));
F3t = S0t - Coefficient[S0t, tt, 2] tt^2;
integrand = Exp[(F3t /. tt -> Sqrt[zz] TT)/zz] (ampt /. tt -> Sqrt[zz] TT);
ser = Normal[Series[integrand, {zz, 0, nmax}]];
res = Expand[ser] /. TT^n_ :> If[OddQ[n], 0, (n - 1)!! cov^(n/2)];
res = res /. TT -> 0;
Do[cn = Together[Coefficient[res, zz, n]];
  target = L^(4 n) (phs[[n + 1]] /. u -> L^(-5));
  Print["c_", n, "(L) = ", Factor[cn], "   minus L^{4n} ph_n(L^-5): ", Simplify[cn - target]], {n, 1, nmax}];

(* --- residue identity --- *)
Emax = 5;
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
