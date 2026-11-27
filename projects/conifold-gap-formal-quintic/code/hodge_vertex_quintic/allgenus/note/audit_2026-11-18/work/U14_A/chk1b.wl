(* U14_A sanity checks: lem:con-cfu, lem:con-red, lem:con-S2, lem:con-id (b), thm:con-Pa (ii)/(iv), residue identity V = Sum pi_m s^-m IC_m *)
U = 7;   (* u-truncation *)
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

(* ---------- model rows Fs_b from T-system ---------- *)
NT = NZ + 2;
Tk[k_, F_] := -F + t/5 (k F - 4 t D[F, t]);
F0 = 1 + Sum[f[n] t^n, {n, 1, NT}];
F1 = Tk[-1, F0]; F2 = -Tk[-2, F1]; F3 = Tk[2, F2]; F4 = -Tk[1, F3];
cl = Expand[Normal[Series[F0 + Tk[0, F4], {t, 0, NT + 1}]]];
solF = Solve[Table[Coefficient[cl, t, n] == 0, {n, 2, NT + 1}], Table[f[n], {n, 1, NT}]];
Fs = Table[Normal[Series[(F /. solF[[1]]), {t, 0, NT}]], {F, {F0, F1, F2, F3, F4}}];
Print["Fs_0 = ", Fs[[1]]];
nu = {0, -1, 3, 2, 1}; epsl = {1, -1, -1, 1, 1};

(* ---------- stationary phase helper (one variable) ---------- *)
(* SPn_c[F;a] with F given by Taylor series in t (F(c+t)-F(c)), amplitude a(t,z); returns series in z to order nz *)
spn[Ft_, at_, z_, nz_, cov_] := Module[{F3t, integrand, ser, res},
  F3t = Ft - Coefficient[Ft, tt, 2] tt^2;
  integrand = Exp[(F3t /. tt -> Sqrt[z] TT)/z] (at /. tt -> Sqrt[z] TT);
  ser = Expand[Normal[Series[integrand, {z, 0, nz}]]];
  res = ser /. TT^n_ :> If[OddQ[n], 0, (n - 1)!! cov^(n/2)];
  res = res /. TT -> 0;
  Expand[res]];

(* ---------- lem:con-S2: SPn_1[-W1; p^nu_b](t) = eps_b Fs_b(t) ---------- *)
W1t = Normal[Series[(1 + tt)^4/4 + 1/(1 + tt), {tt, 0, 2 NZ + 4}]] - 5/4;
Do[
  amp = Normal[Series[(1 + tt)^nu[[b]], {tt, 0, 2 NZ + 4}]];
  g = spn[-W1t, amp, t, NZ, 1/5];
  Print["con-S2 b=", b - 1, ": ", Expand[g - epsl[[b]] Normal[Series[Fs[[b]], {t, 0, NZ}]]]], {b, 1, 5}];

(* ---------- PhiSP from hatted data ---------- *)
Wt = Expand[trunc[Normal[Series[What /. X -> 1 + tt, {tt, 0, 2 NZ + 4}]] - (What /. X -> 1), U]];
at = Expand[trunc[Normal[Series[(ahat /. X -> 1 + tt), {tt, 0, 2 NZ + 4}]] Sqrt[1 - u], U]];
PhiSPser = spn[-Wt, at, zh, NZ, (1 - u)/5];
PhiSPser = Expand[trunc[PhiSPser, U]];
phs = Table[Coefficient[PhiSPser, zh, n], {n, 0, NZ}];
Print["ph^SP_n (trunc u^", U, "): ", phs];

(* ---------- lem:con-id (b) ---------- *)
rhsb = Sum[epsl[[b]] (nu[[b]] /. pis) tau^((1 - nu[[b]])/4) (Fs[[b]] /. t -> tau zh), {b, 1, 5}];
rhsb = Expand[trunc[Normal[Series[rhsb, {zh, 0, NZ}]], U]];
Print["con-id (b) difference: ", Simplify[Expand[PhiSPser - rhsb]]];

(* save for later *)

DumpSave["state.mx", {pis, Vsol, phs, U, NZ, Fs, tau}];
