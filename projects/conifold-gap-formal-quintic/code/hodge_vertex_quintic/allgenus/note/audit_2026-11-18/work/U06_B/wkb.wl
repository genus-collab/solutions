(* lem:gau-wkb (i): the z^0 and z^1 coefficients of PF5^{(i)}; w = zeta^i, w^5 = 1 *)
dq[f_] := D[f, L]*(L^6 - L)/5;
q = (1 - L^(-5))/5^5;
M[f_] := z*dq[f] + w*L*f;
Mpow[f_, n_] := Nest[M, f, n];
prod[f_] := Module[{g = f}, Do[g = 5*M[g] + k*z*g, {k, 1, 5}]; g];
pf5i[f_] := Mpow[f, 5] - f - q*prod[f];
res = Expand[pf5i[g[L]] /. w^n_ :> w^Mod[n, 5]];
c0 = Simplify[Coefficient[res, z, 0] /. w^5 -> 1];
c1 = Coefficient[res, z, 1];
claimT1 = w^4*(5*L^(-1)*dq[g[L]] + (L^(-1) - L^4)*g[L]);
Print["z^0 coefficient: ", c0];
Print["z^1 coefficient minus claim: ", Simplify[c1 - claimT1]];
Print["max z-degree: ", Exponent[res, z]];
(* T1(L g) = 5 w^4 Dq g *)
resLg = Expand[pf5i[L*g[L]] /. w^n_ :> w^Mod[n, 5]];
Print["T1(L g) - 5 w^4 Dq g: ", Simplify[Coefficient[resLg, z, 1] - 5*w^4*dq[g[L]]]];
(* sigma_i intertwining: check on a test series with generic coefficients a0..a3 in L *)
A[zz_] := Sum[a[n][L]*zz^n, {n, 0, 3}];
dqs[f_] := (D[f, L]*(L^6 - L)/5);
M0[f_] := z*dqs[f] + L*f; Mi[f_] := z*dqs[f] + w*L*f;
p0[f_] := Module[{g = Nest[M0, f, 5], h = f}, Do[h = 5*M0[h] + k*z*h, {k, 1, 5}]; g - f - q*h];
pI[f_] := Module[{g = Nest[Mi, f, 5], h = f}, Do[h = 5*Mi[h] + k*z*h, {k, 1, 5}]; g - f - q*h];
lhs = (p0[A[z]]) /. z -> z/w;     (* sigma_i(PF5^(0) A): z -> zeta^{-i} z *)
rhs = pI[A[z/w]];
diff = Expand[(lhs - rhs)*w^10] /. w^n_ :> w^Mod[n, 5];
Print["sigma_i PF5^(0) - PF5^(i) sigma_i : ", Simplify[diff]];
