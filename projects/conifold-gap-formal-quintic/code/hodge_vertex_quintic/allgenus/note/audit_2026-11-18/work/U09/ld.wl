(* Check eq:mod-LD numerically: <phi^b_k, dxi_{a,d}> = -rho_b delta_ab delta_dk, at v=2 (generic), with log p principal branch near each branch point *)
vv = 2; zeta = Exp[2 Pi I/5]; vs = vv^(1/5);
Lg = Symbol["Lg"]; (* stands for log p *)
(* antiderivative operator I on expressions polynomial in Lg with Laurent coefficients in p *)
Iop[F_] := Module[{G, terms, res = 0},
  G = Expand[-F (p^3 - vv p^-2)];
  terms = If[Head[G] === Plus, List @@ G, {G}];
  Do[Module[{c, j, l, tt = term},
     l = Exponent[tt, Lg]; j = Exponent[tt /. Lg -> 1, p]; c = Simplify[tt/(p^j Lg^l)];
     res += c If[j == -1, Lg^(l + 1)/(l + 1),
        (* int p^j Lg^l dp *)
        Sum[(-1)^r l!/(l - r)! Lg^(l - r) p^(j + 1)/(j + 1)^(r + 1), {r, 0, l}]]], {term, terms}];
  Expand[res]];
phi0 = {p, Lg, vv/p, p^3, p^2};
phi[b_, 0] := phi0[[b + 1]];
phi[b_, k_] := phi[b, k] = Iop[phi[b, k - 1]];
(* check D_x phi_{k+1} = -phi_k *)
DxE[f_] := Module[{g = f /. Lg -> Log[p]}, Simplify[p^2/(p^5 - vv) D[g, p]]];
Print["DxCheck: ", Table[Simplify[DxE[phi[b, k + 1]] + (phi[b, k] /. Lg -> Log[p])], {b, 0, 4}, {k, 0, 2}]];
eps1 = {-1, 1, 1, -1, -1}; mm = {4, 0, 1, 2, 3}; alp = {0, 1, 0, 0, 0};
xi[a_] := eps1[[a + 1]] vv^alp[[a + 1]] p^mm[[a + 1]]/(p^5 - vv);
xid[a_, d_] := Nest[Simplify[p^2/(p^5 - vv) D[#, p]] &, xi[a], d];
pairing[b_, k_, a_, d_] := Module[{f = phi[b, k] /. Lg -> Log[p], w = D[xid[a, d], p]},
   Sum[Module[{c = zeta^(-j) vs}, SeriesCoefficient[Series[(f w) /. p -> c + e, {e, 0, 0}], -1]], {j, 0, 4}]];
rho = {-1, 1, -1, -3, -2};
errs = Table[{b, k, a, d, Chop[N[pairing[b, k, a, d] + rho[[b + 1]] If[a == b && d == k, 1, 0], 20], 10^-12]}, {b, 0, 4}, {k, 0, 2}, {a, 0, 4}, {d, 0, 2}];
bad = Select[Flatten[errs, 3], #[[5]] =!= 0 &];
Print["number of failures: ", Length[bad], "  ", Take[bad, Min[5, Length[bad]]]];
