N0 = 12;
DD[f_] := Expand[-(1 - u) D[f, u] + (1 - 1/u) D[f, ell]];
A = {-24/625 (1 - u), -2/5 (1 - u), -7/5 (1 - u), -2 (1 - u), u};
PF4[f_] := Module[{ds}, ds = NestList[DD, f, 4]; Expand[Sum[A[[k + 1]] ds[[k + 1]], {k, 0, 4}]]];
sol[b0_, b1_, b2_] := Module[{cs, g, eqs, ss}, cs = Table[c[n], {n, 3, N0 + 3}];
  g = b0 + b1 u + b2 u^2 + Sum[c[n] u^n, {n, 3, N0 + 3}];
  eqs = Table[Coefficient[PF4[g], u, m] == 0, {m, 0, N0}]; ss = Solve[eqs, cs][[1]];
  Expand[Normal[Series[g /. ss, {u, 0, N0}]]]];
f0 = sol[1, 0, 0]; f1 = sol[0, 1, 0]; f2 = sol[0, 0, 1]; w1 = Expand[f1 + 7/10 f2];
rs = Sum[r[n] u^n, {n, 3, N0 + 3}]; eq3 = PF4[w1 ell + rs];
ss3 = Solve[Table[Coefficient[eq3 /. ell -> 0, u, m] == 0, {m, -3, N0 - 4}], Table[r[n], {n, 3, N0 + 3}]];
y3 = w1 ell + Expand[Normal[Series[rs /. ss3[[1]], {u, 0, N0}]]];
y0 = f0 + 2/3 w1 - 5 f2;
t = Normal[Series[w1/y0, {u, 0, N0}]]; td = DD[t];
dt[h_] := Normal[Series[DD[h]/td, {u, 0, N0 - 3}]];
Yuk = Normal[Series[1/(u y0^2 td^3), {u, 0, N0 - 6}]];
Do[ r2 = Normal[Series[dt[dt[Normal[Series[y/y0, {u, 0, N0}]]]]/Yuk, {u, 0, 4}]];
  (* fit r2 = a + b t *)
  sl = Solve[Table[Coefficient[Normal[Series[r2 - (aa + bb t), {u, 0, 4}]], u, m] == 0, {m, 0, 1}], {aa, bb}];
  Print[y /. {ell -> "ell"} // Short, "  -> a,b: ", sl, "  residual: ", Simplify[Normal[Series[(r2 - (aa + bb t)) /. sl[[1]], {u, 0, 4}]]]],
  {y, {f1, f2, y3, f0}}];
