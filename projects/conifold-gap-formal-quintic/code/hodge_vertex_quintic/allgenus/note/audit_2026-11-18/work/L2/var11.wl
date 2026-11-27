(* Sanity check of thm:mod-var (i) for (h,n)=(1,1) and (0,3), in the convention of conv:mod-kernel:
   d/dv|_x omega_{h,n}(J) = int_{q in ray} omega_{h,n+1}(J,q).
   Curve x = p^4/4 + v/p, y = p, omega02 = dp1 dp2/(p1-p2)^2,
   K_b(p0,q) = 1/2 (1/(p0-q) - 1/(p0-iq)) / ((y(q)-y(iq)) x'(q))  [coefficient of dp0/dq]. *)
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
(* kernel series at branch point c, as function of p0 *)
kern[c_, v_, s_, p0_] := 1/2 (1/(p0 - (c + t)) - 1/(p0 - (c + s)))/(((c + t) - (c + s)) xp[c + t, v]);
res[expr_] := SeriesCoefficient[expr, {t, 0, -1}];
(* omega11 coefficient function at parameter v, point p0 *)
f11[v_, p0_] := Module[{cs = branch[v], tot = 0, s, ds},
   Do[s = invol[c, v]; ds = D[s, t];
     tot += res[Series[kern[c, v, s, p0] ds/(t - s)^2, {t, 0, 2}]], {c, cs}];
   tot];
(* omega11 as exact partial-fraction-free symbolic function of P (needed inside the recursion) *)
f11sym[v_] := Module[{P}, Function @@ {P, f11[v, P]}];
v0 = 13/10; p0 = 7/10 + 2/5 I;
cs = branch[v0];
F11 = f11sym[v0];
(* LHS: d/dv|_x omega11 = d/dv|_p f - d/dp (g f), g = p/(p^5-v) *)
h = 10^-20;
dvf = (f11[v0 + h, p0] - f11[v0 - h, p0])/(2 h);
gf[P_] := P/(P^5 - v0) F11[P];
lhs = dvf - (D[gf[P], P] /. P -> p0);
(* RHS: int over p1 on the ray of omega12(p0,p1), pushed inside the residues:
   int_r dp1/(p1-c)^2 = -1/c, int_r dp1/(a-p1)^2 = -1/a *)
rhs = 0;
Do[s = invol[c, v0]; ds = D[s, t]; q = c + t; iq = c + s;
  w03int = -Sum[1/((q - cb)^2 (iq - cb)^2) (-1/cb)/(5 cb^2), {cb, cs}];
  br = ds (w03int + (-1/q) F11[iq] + F11[q] (-1/iq));
  rhs += res[Series[kern[c, v0, s, p0] br, {t, 0, 6}]], {c, cs}];
Print["(1,1): lhs = ", N[lhs, 25]];
Print["(1,1): rhs = ", N[rhs, 25]];
Print["(1,1): diff = ", N[Abs[lhs - rhs], 5]];
(* residue of omega11 at a branch point, should vanish *)
Print["Res omega11 at c0: ", N[Chop[SeriesCoefficient[Series[F11[cs[[1]] + t], {t, 0, 2}], {t, 0, -1}], 10^-40], 10]];
(* also ray integral of omega11 vs -1/(12 v) (thm:mod-V proof) *)
ray = Exp[I 3/10];
Print["int_ray omega11 = ", N[NIntegrate[F11[r ray] ray, {r, 0, Infinity}, WorkingPrecision -> 40], 20], "  vs -1/(12v) = ", N[-1/(12 v0), 20]];
