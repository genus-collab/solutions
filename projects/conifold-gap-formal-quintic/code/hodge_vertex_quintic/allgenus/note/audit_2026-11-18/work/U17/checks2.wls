(* U17: change of variables (lem:fsp-cov) and iterated stationary phase (lem:fsp-fubini) on explicit examples *)
ord = 3;
(* 1D normalised stationary phase at 0: Gexpr in tt with G'(0)=0; ampexpr in tt, zz *)
spn1[Gexpr_, ampexpr_, ord_] := Module[{s, tau, H, C, ser, res, G0, G1},
  H = D[Gexpr, {tt, 2}] /. tt -> 0; C = -1/H; G0 = Gexpr /. tt -> 0; G1 = D[Gexpr, tt] /. tt -> 0; If[G1 =!= 0, Print["dropping linear term of w-order ", Exponent[G1 /. w -> 1/ww, ww, Min]]];
  ser = Normal[Series[Exp[((Gexpr /. tt -> s tau) - G0 - G1 s tau - H (s tau)^2/2)/s^2] (ampexpr /. {tt -> s tau, zz -> s^2}), {s, 0, 2 ord}]];
  ser = Expand[ser];
  res = ser /. tau^n_. :> If[OddQ[n], 0, (n-1)!! C^(n/2)];
  res = res /. tau -> 0;
  Expand[res /. s -> Sqrt[z]]];

(* --- change of variables, 1D --- *)
F1[x_] := -x^2/2 + x^3/3 + x^4/7 - x^5/2;
a1[x_, zz_] := 1 + x + x^2/3 + zz x;
phi[y_] := 2 y + y^2 - y^3/5 + y^4/3;
J = phi'[0];
lhs = spn1[F1[tt], a1[tt, zz], 4];
rhs = spn1[Normal[Series[F1[phi[tt]], {tt, 0, 14}]], Normal[Series[a1[phi[tt], zz] phi'[tt]/J, {tt, 0, 14}]], 4];
Print["cov 1D: lhs - rhs = ", Simplify[lhs - rhs]];
Print["cov 1D lhs = ", lhs];

(* --- 2D total SPn --- *)
F2[v_, w_] := -v^2/2 - w^2 + v w/3 + v^3/5 + v^2 w/2 - v w^2/7 + w^3/3 + v^4/11 + v^2 w^2/13 + w^4/9 - v^3 w/4;
a2[v_, w_, zz_] := 1 + v + w^2 + v w/2 + zz (1 + w) + zz^2 v^2;
Hm = {{D[F2[v, w], v, v], D[F2[v, w], v, w]}, {D[F2[v, w], w, v], D[F2[v, w], w, w]}} /. {v -> 0, w -> 0};
Cm = -Inverse[Hm];
mom[i_, j_] := mom[i, j] = (D[Exp[{x1, x2}.Cm.{x1, x2}/2], {x1, i}, {x2, j}] /. {x1 -> 0, x2 -> 0});
spn2 = Module[{s, t1, t2, ser, quad},
  quad = {s t1, s t2}.Hm.{s t1, s t2}/2;
  ser = Normal[Series[Exp[(F2[s t1, s t2] - F2[0, 0] - quad)/s^2] a2[s t1, s t2, s^2], {s, 0, 2 ord}]];
  ser = Expand[ser];
  ser = Expand[ser /. {t1^i_. t2^j_. :> mom[i, j]}];
  ser = ser /. {t1^i_. :> mom[i, 0], t2^j_. :> mom[0, j]};
  Expand[ser /. s -> Sqrt[z]]];
Print["2D total SPn = ", spn2];

(* --- iterated: inner v, outer w --- *)
K = 10;
vc = Sum[c[k] w^k, {k, 1, K}];
eqs = CoefficientList[Normal[Series[D[F2[v, w], v] /. v -> vc, {w, 0, K}]], w];
sol = Solve[Thread[eqs == 0], Table[c[k], {k, 1, K}]][[1]];
vcw = vc /. sol;
Fflat = Normal[Series[F2[vcw, w], {w, 0, K}]];
Hin = Normal[Series[D[F2[v, w], v, v] /. v -> vcw, {w, 0, K}]];
Hin0 = Hin /. w -> 0;
bfac = Normal[Series[(Hin/Hin0)^(-1/2), {w, 0, K}]];
Bw = spn1[Expand[F2[vcw + tt, w]], Expand[a2[vcw + tt, w, zz]], ord];
Bw = Expand[Normal[Series[Bw, {w, 0, 2 ord + 2}]]];
outerAmp = Expand[Normal[Series[bfac Bw, {w, 0, 2 ord + 2}]]];
it = spn1[Fflat /. w -> tt, outerAmp /. {w -> tt, z -> zz}, ord];
it = Expand[Normal[Series[it, {z, 0, ord}]]];
tot = Expand[Normal[Series[spn2, {z, 0, ord}]]];
Print["iterated = ", it];
Print["fubini 2D: total - iterated = ", Simplify[tot - it]];
Print["det H = ", Det[Hm], "  det Hvv * S = ", Hm[[1, 1]] (Hm[[2, 2]] - Hm[[2, 1]] Hm[[1, 2]]/Hm[[1, 1]]), "  S = ", (Hm[[2, 2]] - Hm[[2, 1]] Hm[[1, 2]]/Hm[[1, 1]]), "  Fflat''(0) = ", D[Fflat, w, w] /. w -> 0];
