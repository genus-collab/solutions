(* Check of lem:gau-HPF (i): PF5^{(0)}(L Phi) = s^{-1} * (left side of eq:gau-HPF), as an exact identity
   for Phi = sum_{n<=nmax} zh^n ph_n(u) with symbolic ph_n. Variables: s (u = s^5), z; zh = z/s^4. *)
nmax = 3;
Dq[f_] := -(1 - s^5)/(5 s^4) D[f, s];
Mp[f_] := z Dq[f] + f/s;                (* z Dq + L *)
q = (1 - s^5)/5^5;
PF50[f_] := Module[{g = f, h = f, k},
   Do[g = Mp[g], {5}];
   Do[h = 5 Mp[h] + k z h, {k, 1, 5}];
   g - f - q h];
PhiS = Sum[(z/s^4)^n ph[n][s^5], {n, 0, nmax}];
lhs = PF50[PhiS/s];

(* hatted side: variables zh, u *)
thu[f_] := u D[f, u]; thz[f_] := zh D[f, zh];
Oop[e_][f_] := f + zh (1 - u) (4/5 thz[f] - e/5 f - thu[f]);
PhiH = Sum[zh^n ph[n][u], {n, 0, nmax}];
PhiB[0] = PhiH; Do[PhiB[j] = Oop[-j][PhiB[j - 1]], {j, 1, 5}];
hpf = PhiB[5] - PhiH - (1 - u) (3 zh PhiB[4] + 17/5 u zh^2 PhiB[3] + 9/5 u^2 zh^3 PhiB[2] +
      274/625 u^3 zh^4 PhiB[1] + 24/625 u^4 zh^5 PhiH);
rhs = (hpf /. {zh -> z/s^4, u -> s^5})/s;
Print["HPF (i) identity residual: ", Simplify[lhs - rhs]];

(* (ii): coefficient of zh^(n+1) is (1-u)[(4n-5 theta_u) ph_n + r_n]; check the ph_n-part and r_1 = 3/5 for ph_0 = 1 *)
Do[
  c = Coefficient[hpf, zh, n + 1];
  part = c /. {ph[m_] /; m != n :> (0 &)};
  Print["n=", n, "  ph_n part minus (1-u)(4n-5 theta_u)ph_n: ",
   Simplify[part - (1 - u) (4 n ph[n][u] - 5 u ph[n]'[u])]], {n, 0, nmax}];
c1 = Coefficient[hpf, zh, 2] /. {ph[0] -> (1 &), ph[1] -> (0 &)};
Print["r_1 for ph_0 = 1: ", Simplify[c1/(1 - u)]];
