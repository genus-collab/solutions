(* Test: z^{-1}-coefficient of the flat residue V equals -varpi_1 (consequence of prop:con-res + thm:con-Pb with P = Pcan) *)
K = 12;
bk[k_] := BernoulliB[2 k]/(2 k (2 k - 1));
g = Normal[Series[Exp[Sum[bk[j] (5^(1 - 2 j) - 5) w^(2 j - 1), {j, 1, K}]], {w, 0, K}]];
Xi = Sum[u^n/n, {n, 1, K + 1}];
Vm1 = Normal[Series[Sum[(-Xi)^k/k! Coefficient[g, w, k - 1] /. (Coefficient[g, w, 0] -> 1), {k, 1, K}], {u, 0, K}]];
(* careful: Coefficient[g,w,0] *)
g0 = g /. w -> 0;
Vm1 = Normal[Series[Sum[(-Xi)^k/k! If[k == 1, g0, Coefficient[g, w, k - 1]], {k, 1, K}], {u, 0, K}]];
Print["V_{z^-1} = ", Vm1];
th[f_] := -(1 - u) D[f, u];
PF4[f_] := th[th[th[th[f]]]] - (1 - u) (th[f] + f/5) // Function[{y}, y];
PF4b[f_] := Module[{a = f, r},
  r = Nest[th, f, 4] - (1 - u) Fold[#1 + 0 &, 0, {}];
  r];
(* apply theta^4 - (1-u) prod_{k=1..4} (theta + k/5) *)
prodop[f_] := Module[{y = f}, Do[y = th[y] + k/5 y, {k, 4, 1, -1}]; y];
pf = Nest[th, -Vm1, 4] - (1 - u) prodop[-Vm1];
Print["PF4(-V_{z^-1}) = ", Normal[Series[Expand[pf], {u, 0, K - 5}]]];
