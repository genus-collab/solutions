(* Check lem:gau-HPF: PF5^{(0)}(L Phi) = s^{-1} * LHS(HPF), symbolic ph_n(u), up to z^NN *)
NN = 4;
(* z-s coordinates *)
Dq[f_] := -(1 - s^5)/(5 s^4) D[f, s];
M[f_] := z Dq[f] + f/s;
Phizs = Sum[z^n s^(-4 n) ph[n][s^5], {n, 0, NN}];
LPhi = Phizs/s;
q = (1 - s^5)/5^5;
Mpow[f_, k_] := Nest[M, f, k];
term1 = Mpow[LPhi, 5] - LPhi;
(* prod_{k=1}^5 (5M + k z) applied *)
ap[f_, k_] := 5 M[f] + k z f;
prodop = Fold[ap, LPhi, {1, 2, 3, 4, 5}];
pf5 = term1 - q prodop;
pf5ser = Table[Coefficient[Expand[pf5], z, m], {m, 0, NN}];
(* hat coordinates *)
Oe[e_][F_] := F + zh (1 - u) (4/5 zh D[F, zh] - e/5 F - u D[F, u]);
Phih = Sum[zh^n ph[n][u], {n, 0, NN}];
PhiJ[0] = Phih;
Do[PhiJ[j] = Expand[Oe[-j][PhiJ[j - 1]]], {j, 1, 5}];
lhs = PhiJ[5] - Phih - (1 - u) (3 zh PhiJ[4] + 17/5 u zh^2 PhiJ[3] + 9/5 u^2 zh^3 PhiJ[2] + 274/625 u^3 zh^4 PhiJ[1] + 24/625 u^4 zh^5 Phih);
lhszs = (lhs /. {zh -> z/s^4}) /. u -> s^5;
lhsser = Table[Coefficient[Expand[lhszs/s], z, m], {m, 0, NN}];
Print["HPF identity check (should be all 0): ", Simplify[pf5ser - lhsser]];
(* coefficient formula (ii): coefficient of zh^{n+1} = (1-u)[(4n - 5 theta_u) ph_n + r_n] *)
cf[n_] := Coefficient[Expand[lhs], zh, n + 1];
Print["zh^0 coefficient: ", Simplify[Coefficient[Expand[lhs], zh, 0]]];
Do[
  rn = Simplify[cf[n]/(1 - u) - (4 n ph[n][u] - 5 u ph[n]'[u])];
  Print["r_", n, " free of ph_n: ", FreeQ[rn, ph[n]], "  ; r_n with ph0=1: ", Simplify[rn /. {ph[0] -> (1 &)}]],
  {n, 0, 2}];
(* solve HPF at u=0 with polynomial ph_n: ph0=1, ph1=-3/20; get ph2 *)
sol1 = -3/20;
r2 = Simplify[(cf[2]/(1 - u) - (8 ph[2][u] - 5 u ph[2]'[u])) /. {ph[0] -> (1 &), ph[1] -> (-3/20 &)}];
Print["r_2 with ph0=1, ph1=-3/20: ", r2];
(* (8 - 5 theta) ph2 = -r2 *)
c2 = CoefficientList[r2, u];
ph2 = Sum[-c2[[k + 1]]/(8 - 5 k) u^k, {k, 0, Length[c2] - 1}];
Print["ph_2 = ", ph2];
(* T-system closing relation for model F0 *)
T[k_][F_] := -F + t/5 (k F - 4 t D[F, t]);
F0 = 1 + a1 t + a2 t^2 + a3 t^3;
cl = Series[F0 + T[0][T[1][T[2][T[-2][T[-1][F0]]]]], {t, 0, 4}];
Print["T-system: ", Solve[Table[Coefficient[Normal[cl], t, k] == 0, {k, 1, 3}], {a1, a2}]];
