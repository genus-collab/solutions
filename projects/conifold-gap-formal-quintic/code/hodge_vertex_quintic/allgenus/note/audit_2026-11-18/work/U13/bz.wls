(* check eq:con-bz: exp((S(s+z)-S(s))/z) ((s+z)^5-1) A(s+z,z)/A(s,z) == q prod_k (5 s + k z), with exp(S'(s)) = -q (5s)^5/(1-s^5) *)
N0 = 9;
f[x_] := 5 x^4/(1 - x^5);              (* sum_j 1/(zeta^j - x) *)
pw[m_, x_] := D[f[y], {y, m - 1}]/(m - 1)! /. y -> x;  (* sum_j (zeta^j-x)^(-m) *)
bk[k_] := BernoulliB[2 k]/(2 k (2 k - 1));
S2[k_, x_] := bk[k] ((5 x)^(1 - 2 k) + pw[2 k - 1, x]);
S2pp[x_] := 5/(x (1 - x^5));
(* (S(s+z)-S(s))/z minus S'(s) = sum_{n>=2} S^(n)(s) z^(n-1)/n! *)
Srest[s_, z_] := Sum[ (D[S2pp[y], {y, n - 2}] /. y -> s) z^(n - 1)/n!, {n, 2, N0 + 2}];
dS1[s_, z_] := (1/2) Log[((s + z)/s) (1 - s^5)/(1 - (s + z)^5)];
Lhs[s_, z_] := (5 s)^5/(s^5 - 1) ((s + z)^5 - 1) Exp[Srest[s, z] + dS1[s, z] + Sum[z^(2 k - 1) (S2[k, s + z] - S2[k, s]), {k, 1, 5}]];
Rhs[s_, z_] := Product[5 s + k z, {k, 1, 5}];
(* q cancels: LHS has factor -q * (5s)^5/(1-s^5) = q (5s)^5/(s^5-1) *)
Do[ser = Series[Lhs[s0, z] - Rhs[s0, z], {z, 0, 8}]; Print[s0, " -> ", Simplify[Normal[ser]]], {s0, {3/7, -2/5, 11/3}}];
(* symbolic sigma, low order *)
ser = Series[Lhs[s, z]/Rhs[s, z], {z, 0, 5}];
Print[Simplify[ser]];
