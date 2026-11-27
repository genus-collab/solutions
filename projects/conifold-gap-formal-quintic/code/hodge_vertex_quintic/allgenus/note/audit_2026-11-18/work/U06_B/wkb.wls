(* check lem:gau-wkb (i): z^0 and z^1 coefficients of PF5^(i) *)
L = (1 - 5^5 q)^(-1/5);
c = w L;  (* w = zeta^i, w^5 = 1 *)
M[g_] := z q D[g, q] + c g;
prodOp[g_] := Module[{h = g}, Do[h = 5 M[h] + k z h, {k, 5, 1, -1}]; h];
PF5i[g_] := Nest[M, g, 5] - g - q prodOp[g];
expr = Expand[PF5i[f[q]]];
c0 = Coefficient[expr, z, 0];
c1 = Coefficient[expr, z, 1];
red[x_] := Simplify[PowerExpand[x /. w^5 -> 1 /. w^n_ :> w^Mod[n, 5]]];
Print["z^0 coefficient: ", Simplify[c0 /. w^5 -> 1]];
Tw1 = w^4 (5 L^(-1) q f'[q] + (L^(-1) - L^4) f[q]);
Print["z^1 coefficient minus claimed T^w_1 f: ", Simplify[(c1 - Tw1) /. w^5 -> 1 /. w^9 -> w^4 /. w^6 -> w]];
(* T^w_1 (L g) = 5 w^4 D g *)
TwOn[h_] := w^4 (5 L^(-1) q D[h, q] + (L^(-1) - L^4) h);
Print["T^w_1(L g) - 5 w^4 D g: ", Simplify[TwOn[L g[q]] - 5 w^4 q g'[q]]];
(* also check commutation of the five factors (5M + k z) *)
A1[g_] := 5 M[g] + 1 z g; A2[g_] := 5 M[g] + 2 z g;
Print["[A1,A2] f: ", Simplify[A1[A2[f[q]]] - A2[A1[f[q]]]]];
