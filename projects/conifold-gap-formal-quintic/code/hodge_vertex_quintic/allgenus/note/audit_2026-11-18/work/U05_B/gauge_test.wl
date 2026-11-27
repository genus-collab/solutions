(* Independent low-order test of: lem:gau-R1, lem:gau-diag, lem:gau-vanmod/vanshape/vanq, thm:gau-G, prop:gau-symp,
   and (for the true quintic periods) that Phi^can solves HPF.
   mode = "quintic" (periods of PF4, frame (al,be)) or "generic" (random self-dual one-pole system with (N0)). *)
M = 22;          (* truncation order in u for the unknowns *)
K = 6;           (* gauge computed for z^0..z^K (level K+1 used for solvability) *)
OO = M + 4;      (* working order of the data series *)
ser[f_] := Series[f, {u, 0, OO}];
mode = If[Length[$ScriptCommandLine] >= 2, $ScriptCommandLine[[2]], "quintic"];
al = If[Length[$ScriptCommandLine] >= 3, ToExpression[$ScriptCommandLine[[3]]], 0];
be = If[Length[$ScriptCommandLine] >= 4, ToExpression[$ScriptCommandLine[[4]]], 0];
Print["mode = ", mode, "  (alpha,beta) = ", {al, be}];

(* ---------- data ---------- *)
If[mode === "quintic",
  Dop[f_] := -(1 - u) D[f, u];
  PF4[f_] := u Nest[Dop, f, 4] - (1 - u) (2 Nest[Dop, f, 3] + 7/5 Nest[Dop, f, 2] + 2/5 Dop[f] + 24/625 f);
  nf = OO + 3;
  fs = Sum[a[j] u^j, {j, 0, nf}];
  cl = CoefficientList[PF4[fs], u];
  solf = First@Solve[Thread[cl[[1 ;; nf - 2]] == 0], Table[a[j], {j, 3, nf}]];
  fgen = fs /. solf;
  f0 = fgen /. {a[0] -> 1, a[1] -> 0, a[2] -> 0};
  f1 = fgen /. {a[0] -> 0, a[1] -> 1, a[2] -> 0};
  f2 = fgen /. {a[0] -> 0, a[1] -> 0, a[2] -> 1};
  w1 = ser[f1 + 7/10 f2]; w0 = ser[f0] + al w1 + be ser[f2];
  tc = w1/w0;
  Print["tc = ", Series[Normal[tc], {u, 0, 3}]];
  ,
  SeedRandom[12345];
  w0 = ser[1 + Sum[RandomInteger[{-5, 5}]/RandomInteger[{1, 7}] u^j, {j, 1, OO}]];
  tc = ser[u + Sum[RandomInteger[{-5, 5}]/RandomInteger[{1, 7}] u^j, {j, 2, OO + 1}]];
  ];
tcp = D[tc, u];
bq = {tcp, -1/(u w0^2 (1 - u)^3 tcp^2), tcp, -w0/(1 - u), -w0/(1 - u)};   (* b_0..b_4 *)
tau = 1/ser[Sum[4 u^n/(5 n + 4), {n, 0, OO + 2}]];
v = u tau^(-5/4); vp = D[v, u];
bm = {vp, -vp/v, vp, -vp, -vp};                                           (* \tilde b_0..\tilde b_4 *)
Print["(EP) residual: ", Series[Times @@ bq - Times @@ bm, {u, 0, 5}]];
B[a_] := bq[[Mod[a, 5] + 1]]; Bt[a_] := bm[[Mod[a, 5] + 1]];

(* ---------- graded gauge by undetermined coefficients ---------- *)
p[k_, c_] := Sum[x[k, Mod[c, 5], j] u^j, {j, 0, M - 1}];
p[-1, c_] := 0;
eqs = Flatten@Table[
    CoefficientList[
      Normal@Series[u (D[p[k - 1, c], u] - B[c + k - 1] p[k, c] + Bt[c - 1] p[k, c - 1]), {u, 0, M - 1}] + u^M, u][[1 ;; M]],
    {k, 0, K + 1}, {c, 0, 4}];
vars = Flatten@Table[x[k, c, j], {k, 0, K + 1}, {c, 0, 4}, {j, 0, M - 1}];
norm = {x[0, 1, 0] == 1, x[5, 1, 0] == 0};
Off[Solve::svars];
solP = First@Solve[Join[Thread[eqs == 0], norm], vars];
pp[k_, c_] := p[k, c] /. solP;
(* trusted order: largest J such that all coefficients x[k,c,j], k<=K, j<=J are parameter-free *)
free[k_, j_] := Union@Flatten@Table[Variables[x[k, c, j] /. solP], {c, 0, 4}];
Jt = Table[Module[{j = 0}, While[j <= M - 1 && free[k, j] === {}, j++]; j - 1], {k, 0, K}];
Print["trusted u-order of p_k, k=0..K: ", Jt];
J = Min[Jt]; Print["J = ", J];
tr[f_] := Series[f, {u, 0, J}];
P[k_, c_] := tr[pp[k, c]];

(* ---------- T1: (R1) ---------- *)
Print["(R1) column 1: p_{k,1}(0), k=1..K: ", Table[pp[k, 1] /. u -> 0, {k, 1, K}]];
Print["(R1) row 2: p_{k,2-k}(0), k=1..K: ", Table[pp[k, 2 - k] /. u -> 0, {k, 1, K}]];
Print["P_11(0,z) coefficients k=0,5: ", {pp[0, 1] /. u -> 0, pp[5, 1] /. u -> 0},
  "   P_22(0,z): ", {pp[0, 2] /. u -> 0, pp[5, 2] /. u -> 0}];

(* ---------- T2: lem:gau-diag ---------- *)
Print["p_{0,a}(0): ", Table[pp[0, a] /. u -> 0, {a, 0, 4}]];
Print["p_{0,0} - tau^(-1/4)/w0: ", tr[pp[0, 0]] - tr[tau^(-1/4)/w0]];
Print["prod_a p_{0,a} - 1: ", tr[Product[pp[0, a], {a, 0, 4}]] - 1];
Print["p_{0,0} p_{0,3} - 1, p_{0,1} p_{0,2} - 1, p_{0,4}^2 - 1: ",
  {tr[pp[0, 0] pp[0, 3]] - 1, tr[pp[0, 1] pp[0, 2]] - 1, tr[pp[0, 4]^2] - 1}];

(* ---------- vanishing solutions ---------- *)
nmaxY = 3 J + 12;
IC[m_, n_] := Module[{i = (n - m - 1)/5},       (* coefficient of z^-n of I^C_m as a series in v; from the definition *)
   If[IntegerQ[i] && i >= 0 && m + 4 i + 1 >= 0,
    (-1/4)^i/i! (-1)^(m + 4 i + 1) vv^(m + 4 i + 1)/(m + 4 i + 1)!, 0]];
nu = {0, -1, 3, 2, 1}; eps = {1, -1, -1, 1, 1};
YCmod[n_, b_] := eps[[b + 1]] IC[nu[[b + 1]], n] /. vv -> v;      (* [z^-n] Y^C_b(v(u),z) *)
int[f_] := Integrate[Normal[f], u];
yrec[bb_] := Module[{y, tab}, y[0] = {0, -1, 0, 0, 0};
   Do[y[n + 1] = Table[int[Series[bb[[a + 1]] y[n][[Mod[a + 1, 5] + 1]], {u, 0, OO - 2}]], {a, 0, 4}], {n, 0, nmaxY}];
   Table[y[n], {n, 0, nmaxY + 1}]];
yQ = yrec[bq]; yM = yrec[bm];
Print["Y^C_Atilde - Y^C(v(u)) (n<=12): ",
  Union@Flatten@Table[Normal[Series[yM[[n + 1, b + 1]] - YCmod[n, b], {u, 0, J}]], {n, 0, 12}, {b, 0, 4}]];
Print["orders of y^(n)_a for the quintic vs (4n+nu_a+1)/5 (n<=10): ",
  Union@Flatten@Table[
     If[Mod[n - 1 + a, 5] == 0,
      Exponent[Normal[Series[yQ[[n + 1, a + 1]], {u, 0, OO - 3}]], u, Min] - (4 n + nu[[a + 1]] + 1)/5,
      If[Normal[yQ[[n + 1, a + 1]]] === 0, 0, bad]], {n, 0, 10}, {a, 0, 4}]];
Print["[z^-1] Y^C_Q + tc e_0: ", Normal[Series[yQ[[2]] + {tc, 0, 0, 0, 0}, {u, 0, J}]]];

(* ---------- T3: the gauge identity (G) ---------- *)
(* (P Y)_a at z^(k-n) : p_{k,a-k} Ytilde^(n)_{a-k} *)
Gcoef[jz_] := Table[
   Sum[If[n >= 0 && n <= nmaxY, Normal@Series[pp[k, a - k] yM[[n + 1, Mod[a - k, 5] + 1]], {u, 0, J}], 0] /. n -> k - jz,
     {k, 0, K}] /. n -> 0, {a, 0, 4}];
Gc[jz_] := Table[Sum[With[{n = k - jz}, If[n >= 0 && n <= nmaxY,
       Normal@Series[pp[k, a - k] yM[[n + 1, Mod[a - k, 5] + 1]], {u, 0, J}], 0]], {k, 0, K}], {a, 0, 4}];
(* the truncation k<=K neglects P_k ytilde^(k-jz), k>K, of u-order >= (4(k-jz))/5 ; report only below that order *)
Do[
  With[{omit = Ceiling[4 (K + 1 - jz)/5]},
   Module[{lhs = Gc[jz], rhs, d},
    rhs = If[jz <= 0, Normal@Series[yQ[[-jz + 1]], {u, 0, J}], {0, 0, 0, 0, 0}];
    d = Normal@Series[lhs - rhs, {u, 0, Min[J, omit - 1]}];
    Print["(G) at z^", jz, ": residual mod u^", Min[J, omit - 1] + 1, " = ", d]]],
  {jz, -6, 3}];

(* ---------- T4: symplecticity ---------- *)
eta = Table[If[Mod[a + b, 5] == 3, 5, 0], {a, 0, 4}, {b, 0, 4}];
Pmat[k_] := Table[If[Mod[a - b - k, 5] == 0, pp[k, b], 0], {a, 0, 4}, {b, 0, 4}];
Pz[zz_] := Sum[zz^k Pmat[k], {k, 0, K}];
S = Transpose[Pz[-z]].eta.Pz[z];
Sc = Table[Normal@Series[Coefficient[S, z, k], {u, 0, J}], {k, 0, K}];
Print["symplecticity: [z^k](P(-z)^T eta P(z)) - delta_{k0} eta, k=0..K, all zero? ",
  Table[Union@Flatten[Sc[[k + 1]] - If[k == 0, eta, 0 eta]] === {0}, {k, 0, K}]];

(* ---------- T5: Phi^can solves HPF (true quintic only) ---------- *)
If[mode === "quintic",
  nF = K + 2;
  Top[k_][f_] := -f + t/5 (k f - 4 t D[f, t]);
  F0u = 1 + Sum[ff[n] t^n, {n, 1, nF}];
  F1u = Top[-1][F0u]; F2u = -Top[-2][F1u]; F3u = Top[2][F2u]; F4u = -Top[1][F3u]; cl0 = F0u + Top[0][F4u];
  solF = First@Solve[Thread[CoefficientList[cl0, t][[3 ;; nF + 2]] == 0], Table[ff[n], {n, 1, nF}]];
  Fs = {F0u, F1u, F2u, F3u, F4u} /. solF;
  Print["F_0 = ", Series[Fs[[1]], {t, 0, 3}], "   F_b(0) = ", Fs /. t -> 0];
  PhiCan = Sum[
     zh^k Normal@Series[w0 tau^((1 - nu[[Mod[-k, 5] + 1]])/4) u^((4 k + nu[[Mod[-k, 5] + 1]])/5) pp[k, -k], {u, 0, J}]*
      (Fs[[Mod[-k, 5] + 1]] /. t -> Normal[Series[tau, {u, 0, J}]] zh), {k, 0, K}];
  phc[n_] := Normal@Series[Coefficient[PhiCan, zh, n], {u, 0, J}];
  Print["ph^can_n, n=0..3: ", Table[phc[n], {n, 0, 3}]];
  PhiT = Sum[zh^n phc[n], {n, 0, K}];
  thu[f_] := u D[f, u]; thz[f_] := zh D[f, zh];
  Oop[e_][f_] := f + zh (1 - u) (4/5 thz[f] - e/5 f - thu[f]);
  PhiB[0] = PhiT; Do[PhiB[j] = Expand[Oop[-j][PhiB[j - 1]]], {j, 1, 5}];
  hpf = Expand[PhiB[5] - PhiT - (1 - u) (3 zh PhiB[4] + 17/5 u zh^2 PhiB[3] + 9/5 u^2 zh^3 PhiB[2] +
        274/625 u^3 zh^4 PhiB[1] + 24/625 u^4 zh^5 PhiT)];
  Print["HPF residual, coefficients of zh^1..zh^(K+1), mod u^", J + 1, ": ",
   Table[Normal@Series[Coefficient[hpf, zh, n], {u, 0, J}], {n, 1, K + 1}]];
  Print["degrees of ph^can_n in u (polynomiality), n=0..K: ", Table[Exponent[phc[n], u], {n, 0, K}]];
  ];
