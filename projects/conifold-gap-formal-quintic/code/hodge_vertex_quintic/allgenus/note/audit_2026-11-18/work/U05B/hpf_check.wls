(* Independent check of lem:gau-HPF (i): PF5^(0)(L Phi) = s^-1 * HPF(Phi), Phi = sum_n zh^n ph_n(s) *)
NN = 3;
Dq[f_] := -(1 - s^5)/(5 s^4) D[f, s];          (* D = -(1-u) d/du, u = s^5, at fixed z *)
M[f_] := z Dq[f] + f/s;                          (* zD + L *)
Mpow[f_, k_] := Nest[M, f, k];
sig = {3, 17/5, 9/5, 274/625, 24/625};
Phi = Sum[(z/s^4)^n ph[n][s], {n, 0, NN}];
LPhi = Phi/s;
lhs = s^5 Mpow[LPhi, 5] - LPhi - (1 - s^5) Sum[sig[[k]] z^k Mpow[LPhi, 5 - k], {k, 1, 5}];
(* HPF side: work with coefficient lists in zh; theta_u at fixed zh = (s/5) d/ds on coefficients *)
thu[f_] := (s/5) D[f, s];
Oe[e_, c_List] := Module[{n = Length[c]},  (* c = {F_0,...,F_{n-1}} coefficient list; returns list of length n+1 *)
  Table[If[m <= n - 1, c[[m + 1]], 0] + If[m >= 1, (1 - s^5) ((4 (m - 1) - e)/5 c[[m]] - thu[c[[m]]]), 0], {m, 0, n}]];
c0 = Table[ph[n][s], {n, 0, NN}];
P[0] = c0;
Do[P[j] = Oe[-j, P[j - 1]], {j, 1, 5}];
toZ[c_List] := Sum[(z/s^4)^(m - 1) c[[m]], {m, 1, Length[c]}];
hpf = toZ[P[5]] - toZ[P[0]] - (1 - s^5) (3 (z/s^4) toZ[P[4]] + 17/5 s^5 (z/s^4)^2 toZ[P[3]] + 9/5 s^10 (z/s^4)^3 toZ[P[2]] + 274/625 s^15 (z/s^4)^4 toZ[P[1]] + 24/625 s^20 (z/s^4)^5 toZ[P[0]]);
diff = Simplify[Expand[lhs - hpf/s]];
Print["difference PF5^(0)(L Phi) - s^-1 HPF(Phi): ", diff];
(* r_1 with ph0 = 1 and the pivot (4n - 5 theta_u) *)
sub = {ph[0] -> Function[x, 1]};
coef1 = Simplify[Coefficient[Expand[hpf /. sub], z, 2] s^8];  (* zh^2 coefficient *)
Print["zh^2 coefficient (ph0=1): ", Collect[coef1, {ph[1][s], ph[1]'[s]}]];
(* zh^1 coefficient with ph0=1 *)
coefA = Simplify[Coefficient[Expand[hpf /. sub], z, 1] s^4];
Print["zh^1 coefficient (ph0=1): ", Collect[Expand[coefA], {ph[1][s]}]];
