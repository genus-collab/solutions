(* Check lem:mod-closure, lem:mod-genfn, lem:mod-resid, prop:mod-dii on small cases *)
X[p_] := p^4/4 + v/p;
(* termwise antiderivative of a polynomial in p, 1/p, Log[p] *)
antider[F_] := Module[{e = Expand[F], terms},
  terms = If[Head[e] === Plus, List @@ e, {e}];
  Total[Map[Integrate[#, p] &, terms]]];
Iop[F_] := Expand[antider[Expand[-F (p^3 - v p^-2)]]];
(* remove a constant term p^0 Log^0 if Integrate produced one -- it should not *)
phi0 = {p, Log[p], v/p, p^3, p^2};
rho = {-1, 1, -1, -3, -2};
nu = {0, -1, 3, 2, 1}; eps = {1, -1, -1, 1, 1};
K = 7;
legs = Table[NestList[Iop, phi0[[b + 1]], K], {b, 0, 4}];
(* check: max power of Log *)
Print["max log power: ", Max[Table[Exponent[legs[[b, k]] /. Log[p] -> LL, LL], {b, 1, 5}, {k, 1, K + 1}]]];
(* check: constant term p^0 Log^0 absent *)
Print["constant terms: ", Union[Flatten[Table[Coefficient[Coefficient[Expand[legs[[b, k]] /. Log[p] -> LL], LL, 0], p, 0], {b, 1, 5}, {k, 2, K + 1}]]]];
(* check D_x phi_{k+1} = - phi_k *)
Print["Dx check: ", Union[Flatten[Table[Simplify[p^2/(p^5 - v) D[legs[[b, k + 1]], p] + legs[[b, k]]], {b, 1, 5}, {k, 1, K}]]]];
(* Pi_0(z^m H^b) from (eq:tra-values) *)
Pi0[m_, b_] := If[Mod[b + m, 5] == 1, With[{r = Floor[m/5]}, -eps[[b + 1]] v^(m - r)/(4^r r! (m - r)!)], 0];
(* also directly from residue definition *)
Pi0res[m_, b_] := -eps[[b + 1]]/m! Coefficient[Expand[p^(nu[[b + 1]]) X[p]^m p], p, 0];
Print["Pi0 consistency: ", Union[Flatten[Table[Simplify[Pi0[m, b] - Pi0res[m, b]], {m, 0, 12}, {b, 0, 4}]]]];
(* R^b_n = Res W^n dphi0 / n! *)
Rb[b_, n_] := Coefficient[Expand[X[p]^n D[phi0[[b + 1]], p] p], p, 0]/n!;
Print["resid identity: ", Union[Flatten[Table[Simplify[Rb[b, n] - rho[[b + 1]] Pi0[n, b]], {b, 0, 4}, {n, 0, 12}]]]];
(* log p coefficient vs formula *)
logcoef[F_] := Coefficient[Expand[F /. Log[p] -> LL], LL, 1];
dii = Table[Simplify[Expand[logcoef[legs[[b + 1, k + 1]]] - rho[[b + 1]] Sum[Pi0[k - i, b] (-X[p])^i/i!, {i, 0, k}]]], {b, 0, 4}, {k, 0, K}];
Print["dii check: ", Union[Flatten[dii]]];
(* the non-log part lies in Q[p, pbar]: every term v^m p^j has j >= -m *)
nonlog[F_] := Coefficient[Expand[F /. Log[p] -> LL], LL, 0];
okQ[F_] := Module[{e = Expand[F], terms}, If[e === 0, True,
   terms = If[Head[e] === Plus, List @@ e, {e}];
   And @@ Map[(Exponent[#, p] + Exponent[#, v] >= 0 || (Exponent[#, p] >= 0)) &, terms]]];
(* careful: Exponent of p for negative powers *)
okQ2[F_] := Module[{e = Expand[F], cl}, If[e === 0, True,
   cl = CoefficientRules[e /. p -> pp, {pp, v}] /. {}; True]];
minexp[F_] := Module[{e = Expand[F]}, If[e === 0, Infinity,
   Min[Map[(Exponent[#, p, Min] + Exponent[#, v]) &, If[Head[e] === Plus, List @@ e, {e}]]]]];
Print["min (j+m) over non-log terms: ", Min[Flatten[Table[minexp[nonlog[legs[[b, k]]]], {b, 1, 5}, {k, 1, K + 1}]]]];
Print["phi^2_1 = ", legs[[3, 2]]];
Print["phi^1_1 = ", legs[[2, 2]]];
Print["phi^1_2 = ", legs[[2, 3]]];
