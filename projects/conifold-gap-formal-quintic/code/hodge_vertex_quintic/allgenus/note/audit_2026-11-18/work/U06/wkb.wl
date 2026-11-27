(* check T^w_1 and z^0 coefficient of PF5^{(i)} ; c = zeta^i with c^5 = 1 *)
L = (1 - 5^5 q)^(-1/5);
w = c L;
Mop[f_] := z q D[f, q] + w f;
pow[f_, 0] := f; pow[f_, n_] := Mop[pow[f, n - 1]];
prodop[f_] := Fold[Function[{g, k}, 5 Mop[g] + k z g], f, {5, 4, 3, 2, 1}];
PF5i[f_] := pow[f, 5] - f - q prodop[f];
expr = PF5i[L g[q] + z h[q]];
c0 = Simplify[Coefficient[Expand[expr], z, 0] /. c^5 -> 1];
c1 = Simplify[(Coefficient[Expand[expr], z, 1] /. {c^5 -> 1, c^6 -> c, c^9 -> c^4, c^10 -> 1}) ];
Print["z^0 coefficient: ", Simplify[c0 /. c^5 -> 1]];
(* z^1 coefficient = T1(h) + T2? no: z^1 coeff of PF5i(L g + z h) = T1(L g) + T0(h) where T0 = 0 *)
Print["z^1 coefficient minus 5 c^4 q g': ", Simplify[PowerExpand[c1 - 5 c^4 q g'[q]] /. {c^5 -> 1, c^6 -> c, c^9 -> c^4}]];
