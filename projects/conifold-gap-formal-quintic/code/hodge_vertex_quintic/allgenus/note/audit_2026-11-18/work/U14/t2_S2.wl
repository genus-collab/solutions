(* Test lem:con-S2: SPn_{p=1}[-W_1; p^nu_b](t) = eps_b Fs_b(t) *)
K = 5; (* order in t *)
W1[p_] := p^4/4 + 1/p;
Fser = Normal[Series[-W1[1 + x] + 5/4, {x, 0, 2 K + 2}]];
c2 = Coefficient[Fser, x, 2]; (* = -5/2 *)
Cov = -1/(2 c2);
F3 = Fser - c2 x^2;
mom[j_] := If[OddQ[j], 0, (j - 1)!! Cov^(j/2)];
SPn[amp_] := Module[{e, s, cl},
  e = Exp[(F3 /. x -> w tt)/w^2] (amp /. x -> w tt);
  s = Normal[Series[e, {w, 0, 2 K}]];
  s = Expand[s];
  s = s /. tt^j_. :> mom[j];
  s = s /. tt -> 0;
  Expand[s /. w -> Sqrt[t]]];
nu = {0, -1, 3, 2, 1}; eps = {1, -1, -1, 1, 1};
(* T-system for colour 0 *)
Tsys[k_, F_] := -F + t/5 (k F - 4 t D[F, t]);
F0 = 1 + Sum[a[n] t^n, {n, 1, K + 2}];
F1 = Tsys[-1, F0]; F2 = -Tsys[-2, F1]; F3t = Tsys[2, F2]; F4 = -Tsys[1, F3t];
clos = Expand[F0 + Tsys[0, F4]];
eqs = Table[Coefficient[clos, t, n] == 0, {n, 1, K + 2}];
sa = Solve[eqs, Table[a[n], {n, 1, K + 2}]][[1]];
Fs = Map[Normal[Series[# /. sa, {t, 0, K}]] &, {F0, F1, F2, F3t, F4}];
Print["Fs_0 = ", Fs[[1]]];
Do[
  g = SPn[(1 + x)^nu[[b + 1]]];
  diff = Expand[Normal[Series[g - eps[[b + 1]] Fs[[b + 1]], {t, 0, K}]]];
  Print["b=", b, " SPn = ", Normal[Series[g, {t, 0, 3}]], "   diff to t^", K, ": ", diff],
  {b, 0, 4}];
