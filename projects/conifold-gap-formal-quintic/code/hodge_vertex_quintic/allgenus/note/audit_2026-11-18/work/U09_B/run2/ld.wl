(* Check of the leg dictionary prop:mod-LD: <phi^b_k, d xi_{a,d}> = -rho_b delta_ab delta_dk, at v = 2. *)
vv = 2;
(* legs: expressions in p and L = log p *)
phi0 = {p, L, vv/p, p^3, p^2};
Iop[F_] := Module[{G},
  (* termwise antiderivative of -F (p^3 - v p^-2) dp, with d L = dp/p, no constant term *)
  G = Integrate[-(F /. L -> Log[p]) (p^3 - vv p^-2), p];
  G = Expand[G /. Log[p] -> L];
  G - (G /. {p -> 0, L -> 0} /. {} ) // Expand];
(* safer: remove pure constant term *)
Iop2[F_] := Module[{G = Expand[Integrate[-(F /. L -> Log[p]) (p^3 - vv p^-2), p] /. Log[p] -> L]},
   Expand[G - Coefficient[Coefficient[G /. p -> pp, L, 0] /. pp -> p, p, 0] ]];
phi[b_, k_] := Nest[Iop2, phi0[[b + 1]], k];
mm = {4, 0, 1, 2, 3}; ep = {-1, 1, 1, -1, -1}; al = {0, 1, 0, 0, 0};
xi[a_] := ep[[a + 1]] vv^al[[a + 1]] p^mm[[a + 1]]/(p^5 - vv);
Dx[f_] := p^2/(p^5 - vv) D[f, p];
xid[a_, d_] := Nest[Dx, xi[a], d];
cs = Table[N[Exp[-2 Pi I j/5] vv^(1/5), 40], {j, 0, 4}];
pair[F_, a_, d_] := Module[{form = (F /. L -> Log[p]) D[xid[a, d], p]},
  Sum[SeriesCoefficient[form /. p -> c + s, {s, 0, -1}], {c, cs}]];
rho = {-1, 1, -1, -3, -2};
res = Table[Chop[pair[phi[b, k], a, d], 10^-25] - If[a == b && d == k, -rho[[b + 1]], 0], {b, 0, 4}, {k, 0, 2}, {a, 0, 4}, {d, 0, 2}];
Print["max |<phi^b_k, d xi_{a,d}> + rho_b delta delta| over b,a in Z/5, k,d <= 2: ", Max[Abs[Flatten[res]]]];
Print["phi^1_1 = ", phi[1, 1]];
Print["phi^0_1 = ", phi[0, 1]];
