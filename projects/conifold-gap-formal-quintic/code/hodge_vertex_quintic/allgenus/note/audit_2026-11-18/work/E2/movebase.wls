(* Independent check of cor:fsp-movebase in one variable (q-free amplitude, and a q-dependent one
   with non-uniform z-order, read with the extended definition SPn_c[F;sum q^d a_d] := sum q^d SPn_c[F;a_d]). *)
K = 3; Qo = 3;
mom[n_, C_] := If[OddQ[n], 0, (n - 1)!! C^(n/2)];
Eop[expr_, C_] := Module[{e = Expand[expr], cl},
  cl = CoefficientList[e, \[Tau]];
  Sum[cl[[i]] mom[i - 1, C], {i, Length[cl]}]];
F0[t_] := -t^2/2 + t^3/3 + t^4/5;
F1[t_] := q (t + t^2/2) + q^2 t^3;
F[t_] := F0[t] + F1[t];
(* LHS: functional at c0 = 0 with C0 = 1, amplitude e^{F1/z} a, coefficientwise in q *)
lhs[a_] := Module[{int, s},
  int = Normal[Series[Exp[(F0[r \[Tau]] + (r \[Tau])^2/2)/r^2] Exp[F1[r \[Tau]]/r^2] a[r \[Tau], r^2], {q, 0, Qo}]];
  s = Normal[Series[int, {r, 0, 2 K}]];
  Expand[Eop[s, 1]]];
(* RHS: critical point c(q), Hessian, cubic part, prefactors *)
csol = Module[{cc = 0}, Do[cc = Normal[Series[cc - (D[F[x], x] /. x -> cc)/(D[F[x], {x, 2}] /. x -> 0 /. q -> 0) , {q, 0, Qo}]], {12}]; cc];
Print["check crit: ", Simplify[Normal[Series[D[F[x], x] /. x -> csol, {q, 0, Qo}]]]];
Hc = Normal[Series[D[F[x], {x, 2}] /. x -> csol, {q, 0, Qo}]];
Cc = Normal[Series[-1/Hc, {q, 0, Qo}]];
rhs[a_] := Module[{Fc3, int, s, spn, pref},
  Fc3[t_] := Normal[Series[F[csol + t] - F[csol] - Hc t^2/2, {q, 0, Qo}]];
  int = Exp[Fc3[r \[Tau]]/r^2] a[csol + r \[Tau], r^2];
  s = Normal[Series[int, {r, 0, 2 K + 2 Qo}]];
  spn = Normal[Series[Eop[s, Cc], {q, 0, Qo}]];
  pref = Normal[Series[Exp[(F[csol] - F0[0])/r^2] (Hc/(-1))^(-1/2), {q, 0, Qo}]];
  Expand[Normal[Series[pref spn, {q, 0, Qo}]]]];
cmp[a_] := Module[{L = lhs[a], R = rhs[a], d},
  d = Expand[L - R];
  Table[{j, k, Coefficient[Coefficient[d, q, j], r, 2 k]}, {j, 0, Qo}, {k, -Qo, K}]];
a1[t_, z_] := 1 + t + z t^2;
Print["q-free amplitude: ", Union[Flatten[cmp[a1][[All, All, 3]]]]];
(* q-dependent amplitude with non-uniform z-order: a = (1+t) sum_d q^d z^{-d} *)
a2[t_, z_] := (1 + t) Sum[q^d z^(-d), {d, 0, Qo}];
Print["q-dependent amplitude: ", Union[Flatten[cmp[a2][[All, All, 3]]]]];
(* negative control: drop the Hessian prefactor *)
rhsBad[a_] := Module[{Fc3, int, s, spn, pref},
  Fc3[t_] := Normal[Series[F[csol + t] - F[csol] - Hc t^2/2, {q, 0, Qo}]];
  int = Exp[Fc3[r \[Tau]]/r^2] a[csol + r \[Tau], r^2];
  s = Normal[Series[int, {r, 0, 2 K + 2 Qo}]];
  spn = Normal[Series[Eop[s, Cc], {q, 0, Qo}]];
  pref = Normal[Series[Exp[(F[csol] - F0[0])/r^2], {q, 0, Qo}]];
  Expand[Normal[Series[pref spn, {q, 0, Qo}]]]];
dbad = Expand[lhs[a1] - rhsBad[a1]];
Print["negative control nonzero: ", Coefficient[Coefficient[dbad, q, 1], r, 0]];
