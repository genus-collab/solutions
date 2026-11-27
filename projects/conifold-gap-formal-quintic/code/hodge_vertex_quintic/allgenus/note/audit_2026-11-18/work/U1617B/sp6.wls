(* independent 2-D check of cor:fsp-movebase (with lem:fsp-pg (iii), lem:fsp-transl (a)) *)
Kz = 2; Kq = 2;
mom2[i_, j_, C_] := If[OddQ[i + j], 0, Sum[If[EvenQ[i - k] && EvenQ[j - k],
    Binomial[i, k] Binomial[j, k] k! Max[1, (i - k - 1)!!] Max[1, (j - k - 1)!!] C[[1, 1]]^((i - k)/2) C[[2, 2]]^((j - k)/2) C[[1, 2]]^k, 0], {k, 0, Min[i, j]}]];
trunc[e_] := Normal[Series[e, {q, 0, Kq}]];
spn[P3_, amp_, C_, K_] := Module[{poly, res},
  poly = Expand[Normal[Series[Exp[P3[h t1, h t2]/h^2] amp[h t1, h t2, h^2], {h, 0, 2 K}]]];
  poly = Expand[poly /. q^n_ /; n > Kq -> 0];
  res = poly /. {t1^i_. t2^j_. :> mom2[i, j, C], t1^i_. :> mom2[i, 0, C], t2^j_. :> mom2[0, j, C]};
  Expand[trunc[Normal[Series[res, {h, 0, 2 K}]] /. h -> Sqrt[z]]]];
F0[x1_, x2_] := -x1^2/2 - x2^2 + x1 x2/3 + x1^3/5 + x1 x2^2/2 - x2^4/7 + x1^4/11;
G[x1_, x2_] := x1 - 2 x2 + x1 x2 + x2^3/3;
am[x1_, x2_, zz_] := 1 + x1 + zz x2^2 - x1 x2/2;
H0 = {{-1, 1/3}, {1/3, -2}}; C0 = -Inverse[H0];
P30[t1_, t2_] := F0[t1, t2] - ({t1, t2}.H0.{t1, t2})/2;
lhs = Expand[Sum[q^d/d! z^(-d) spn[P30, Function[{a1, a2, zz}, G[a1, a2]^d am[a1, a2, zz]], C0, Kz + d], {d, 0, Kq}]];
Print["lhs done"];
Ff[x1_, x2_] := F0[x1, x2] + q G[x1, x2];
sol = {};
Do[eqs = Thread[Coefficient[Normal[Series[{D[Ff[x1, x2], x1], D[Ff[x1, x2], x2]} /. {x1 -> Sum[a1[k] q^k, {k, 1, Kq + 1}], x2 -> Sum[a2[k] q^k, {k, 1, Kq + 1}]} /. sol, {q, 0, kk}]], q, kk] == 0];
   sol = Join[sol, First[Solve[eqs, {a1[kk], a2[kk]}]]], {kk, 1, Kq + 1}];
cq = {Sum[a1[k] q^k, {k, 1, Kq + 1}], Sum[a2[k] q^k, {k, 1, Kq + 1}]} /. sol;
Print["c(q) = ", cq];
FcPoly = Expand[Ff[cq[[1]] + t1, cq[[2]] + t2]];
FcPoly = Expand[trunc[FcPoly]];
grad = {Coefficient[Coefficient[FcPoly, t1, 1], t2, 0], Coefficient[Coefficient[FcPoly, t2, 1], t1, 0]};
Print["gradient at c (should be O(q^3)): ", trunc[grad]];
Hc = {{2 Coefficient[Coefficient[FcPoly, t1, 2], t2, 0], Coefficient[Coefficient[FcPoly, t1, 1], t2, 1]},
      {Coefficient[Coefficient[FcPoly, t1, 1], t2, 1], 2 Coefficient[Coefficient[FcPoly, t2, 2], t1, 0]}};
Cc = Map[trunc, -Inverse[Hc], {2}];
crit = trunc[(FcPoly /. {t1 -> 0, t2 -> 0}) - F0[0, 0]];
P3cPoly = Expand[FcPoly - (FcPoly /. {t1 -> 0, t2 -> 0}) - grad[[1]] t1 - grad[[2]] t2 - ({t1, t2}.Hc.{t1, t2})/2];
P3c[s1_, s2_] := P3cPoly /. {t1 -> s1, t2 -> s2};
ampc[s1_, s2_, zz_] := Expand[am[cq[[1]] + s1, cq[[2]] + s2, zz]];
spc = spn[P3c, ampc, Cc, Kz + Kq]; (* exp(crit/z) has z^-k, so spc is needed to z^(Kz+Kq) *)
detr = trunc[(Det[Hc]/Det[H0])^(-1/2)];
rhs = Expand[trunc[Expand[Series[Exp[crit/z], {q, 0, Kq}] // Normal] detr spc]];
diff = Expand[lhs - rhs];
Print["coefficients of diff, q^0..q^2 (rows), z^-2..z^2 (cols): ", Table[Simplify[Coefficient[Coefficient[diff, q, d], z, k]], {d, 0, Kq}, {k, -Kq, Kz}]];
Print["sample lhs coefficients q^2: ", Table[Coefficient[Coefficient[lhs, q, 2], z, k], {k, -2, 2}]];
