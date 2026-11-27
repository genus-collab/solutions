(* E1: independent checks of the factor-by-factor step behind eq:con-bz.  Own code; no project script. *)
K = 8; (* number of Stirling terms *)
bk[k_] := BernoulliB[2 k]/(2 k (2 k - 1));
mu[w_] := Sum[bk[k] w^(2 k - 1), {k, 1, K}];
ord = 14;

(* (A) functional equation rho = 0 and its 5-fold iterate in Q[[W]], W = 1/Y *)
rho = (1/W + 1/2) Log[1 + W] - 1 + mu[W/(1 + W)] - mu[W];
rhoN = (1/W + 11/2) Log[1 + 5 W] - 5 + mu[W/(1 + 5 W)] - mu[W] - Sum[Log[1 + k W], {k, 1, 5}];
Print["(A) rho  series: ", Normal[Series[rho, {W, 0, ord}]] // Simplify];
Print["(A) rhoN series: ", Normal[Series[rhoN, {W, 0, ord}]] // Simplify];

(* (B) change of base point: T_z (t -> t+z) applied to the images of W and of log Y under eps_N, eps_D *)
sig = L + t;
epsNW = z/(5 sig);  tauNW = W/(1 + 5 W);
Print["(B) T_z eps_N(W) - eps_N(tau_5 W): ",
  Simplify[(epsNW /. t -> t + z) - (tauNW /. W -> epsNW)]];
epsDW = z/(a + t); tauDW = W/(1 + W);
Print["(B) T_z eps_D(W) - eps_D(tau_1 W): ",
  Simplify[(epsDW /. t -> t + z) - (tauDW /. W -> epsDW)]];
n = 6;
Print["(B) log-symbol N: ",
  Normal[Series[(Log[1 + (t + z)/L]) - (Log[1 + t/L] + Log[1 + 5 W] /. W -> epsNW), {t, 0, n}, {z, 0, n}]] // Simplify];
Print["(B) log-symbol D: ",
  Normal[Series[(Log[1 + (t + z)/a]) - (Log[1 + t/a] + Log[1 + W] /. W -> epsDW), {t, 0, n}, {z, 0, n}]] // Simplify];

(* (D) factor by factor in the symbol ring: lamN, lam, ell are free symbols (log(5L), log(L - zeta^j), log z) *)
FN = (5 sig/z + 1/2) (lamN + Log[1 + t/L] - ell) - 5 sig/z + Sum[bk[k] (z/(5 sig))^(2 k - 1), {k, 1, K}];
claimN = 5 (lamN - ell) + Sum[Log[1 + (5 t + k z)/(5 L)], {k, 1, 5}];
dN = (FN /. t -> t + z) - FN - claimN;
Print["(D) N-factor: ", Normal[Series[dN, {t, 0, 5}, {z, 0, 5}]] // Simplify];
FD = -((((a + t)/z) + 1/2) (lam + Log[1 + t/a] - ell) - (a + t)/z + Sum[bk[k] (z/(a + t))^(2 k - 1), {k, 1, K}]);
claimD = -(lam - ell) - Log[1 + (t + z)/a];
dD = (FD /. t -> t + z) - FD - claimD;
Print["(D) D-factor: ", Normal[Series[dD, {t, 0, 5}, {z, 0, 5}]] // Simplify];
F0 = sig lamq/z;
Print["(D) q-factor: ", Simplify[(F0 /. t -> t + z) - F0 - lamq]];

(* (C) end to end: the logarithmic form of eq:con-bz over Q(L), all rational/log data normalised at L *)
q = (1 - L^-5)/5^5;
S0p[s_] := 5 Log[s/L] - Log[(1 - s^5)/(1 - L^5)];            (* S_0'(sigma), S_0'(L) = 0 *)
S1[s_] := (1/2) Log[s/L] - (1/2) Log[(1 - s^5)/(1 - L^5)];   (* S_1(sigma) - S_1(L) *)
pr[r_, s_] := (1/(r - 1)!) D[5 x^4/(1 - x^5), {x, r - 1}] /. x -> s;  (* sum_j (zeta^j - s)^(-r) *)
S2k[k_, s_] := bk[k] ((5 s)^(1 - 2 k) + pr[2 k - 1, s]);
N0 = 6;
Eexp = Sum[(D[S0p[x], {x, m - 1}] /. x -> sig) z^(m - 1)/m!, {m, 1, N0 + 1}];
amp[s_] := S1[s] + Sum[z^(2 k - 1) S2k[k, s], {k, 1, 4}];
lhs = Eexp + (amp[sig + z] - amp[sig]);
rhs = Log[q Product[5 sig + k z, {k, 1, 5}]/((sig + z)^5 - 1)];
ser = Series[lhs - rhs, {t, 0, 4}, {z, 0, N0}];
Print["(C) log eq:con-bz, symbolic L, t^<=4 z^<=6: ", Simplify[Normal[ser]]];
