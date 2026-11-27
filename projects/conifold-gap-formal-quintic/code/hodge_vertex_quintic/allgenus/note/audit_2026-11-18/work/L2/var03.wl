(* (0,3): d/dv|_x omega03(p0,p1,p2) = int_{p3 in ray} omega04(p0,p1,p2,p3) *)
prec = 80; ord = 10;
xF[p_, v_] := p^4/4 + v/p;
xp[p_, v_] := p^3 - v/p^2;
invol[c_, v_] := Module[{a, s, eq, sol = {}, F, u},
   F[w_] := Normal[Series[xF[c + u, v] - xF[c, v], {u, 0, ord + 3}]] /. u -> w;
   s = -t + Sum[a[k] t^k, {k, 2, ord}];
   eq = CoefficientList[Normal[Series[F[s] - F[t], {t, 0, ord + 1}]], t];
   Do[sol = Join[sol, Solve[(eq[[k + 2]] /. sol) == 0, a[k]][[1]]], {k, 2, ord}];
   Expand[s /. sol]];
branch[v_] := Table[N[v^(1/5), prec] Exp[-2 Pi I b/5], {b, 0, 4}];
kern[c_, v_, s_, p0_] := 1/2 (1/(p0 - (c + t)) - 1/(p0 - (c + s)))/(((c + t) - (c + s)) xp[c + t, v]);
res[expr_] := SeriesCoefficient[expr, {t, 0, -1}];
w02[a_, b_] := 1/(a - b)^2;
w03[v_][a_, b_, c_] := Module[{cs = branch[v]}, -Sum[1/((a - cb)^2 (b - cb)^2 (c - cb)^2)/(5 cb^2), {cb, cs}]];
v0 = 13/10; p0 = 7/10 + 2/5 I; p1 = -3/10 + 9/10 I; p2 = 1/5 - 6/5 I;
cs = branch[v0];
(* LHS: d/dv|_p w03 - sum_l d/dp_l (g w03) in each slot *)
h = 10^-20;
dv = (w03[v0 + h][p0, p1, p2] - w03[v0 - h][p0, p1, p2])/(2 h);
g[P_] := P/(P^5 - v0);
lie = (D[g[P] w03[v0][P, p1, p2], P] /. P -> p0) + (D[g[P] w03[v0][p0, P, p2], P] /. P -> p1) + (D[g[P] w03[v0][p0, p1, P], P] /. P -> p2);
lhs = dv - lie;
(* RHS: int over p3 on the ray of omega04(p0,p1,p2,p3), integral pushed inside:
   int_r w02(a,p3) = -1/a ; int_r w03(a,b,p3) = -sum_cb e(a)e(b)(-1/cb)/(5cb^2) *)
I02[a_] := -1/a;
I03[a_, b_] := -Sum[1/((a - cb)^2 (b - cb)^2) (-1/cb)/(5 cb^2), {cb, cs}];
W03 = w03[v0];
rhs = 0;
Do[s = invol[c, v0]; ds = D[s, t]; q = c + t; iq = c + s;
  (* bracket of omega04(p0; p1,p2,p3): sum over splittings I1|I2 of {p1,p2,p3}, both nonempty;
     slot p3 integrated over the ray *)
  br = ds (
     w02[q, p1] I03[iq, p2] + w02[q, p2] I03[iq, p1] + I02[q] W03[iq, p1, p2] +
     w02[iq, p1] I03[q, p2] + w02[iq, p2] I03[q, p1] + I02[iq] W03[q, p1, p2]);
  rhs += res[Series[kern[c, v0, s, p0] br, {t, 0, 6}]], {c, cs}];
Print["(0,3): lhs = ", N[lhs, 20]];
Print["(0,3): rhs = ", N[rhs, 20]];
Print["(0,3): diff = ", N[Abs[lhs - rhs], 5]];
