(* TR side at v=1 for x=p^4/4+1/p, y=p, kernel of conv:mod-kernel *)
$MaxExtraPrecision = 200;
prec = 60;
z5 = Exp[2 Pi I/5];
xf[p_] := p^4/4 + 1/p;
cs = Table[N[z5^(-b), prec], {b, 0, 4}];
Nord = 12;
(* local involution at c: iota(eps) = c + s(eps), s = -eps + ... *)
iotaSeries[c_] := Module[{s, a, eqs, sol, ser},
  s = -e + Sum[a[k] e^k, {k, 2, Nord}];
  ser = Series[xf[c + s] - xf[c + e], {e, 0, Nord + 1}];
  sol = Solve[Table[SeriesCoefficient[ser, k] == 0, {k, 3, Nord + 1}], Table[a[k], {k, 2, Nord}]][[1]];
  s /. sol];
iotas = Table[iotaSeries[cs[[b]]], {b, 1, 5}];
(* omega_{1,1}(p0)/dp0 *)
om11TR[p0_] := Sum[Module[{c = cs[[b]], s = iotas[[b]], q, iq, integrand},
    q = c + e; iq = c + s;
    integrand = 1/2 (1/(p0 - q) - 1/(p0 - iq)) D[s, e]/((q - iq)^3 (D[xf[w], w] /. w -> q));
    SeriesCoefficient[Series[integrand, {e, 0, 1}], -1]], {b, 1, 5}];
(* model prediction with v=1 *)
DxF[f_] := p^2/(p^5 - 1) D[f, p];
xi0 = -p^4/(p^5 - 1); xi1 = 1/(p^5 - 1); xi2 = p/(p^5 - 1); xi3 = -p^2/(p^5 - 1); xi4 = -p^3/(p^5 - 1);
xis = {xi0, xi1, xi2, xi3, xi4};
pred11[p0_] := (5/24 D[DxF[xi0], p] - 1/12 D[xi1, p]) /. p -> p0;
pts = {N[3/10 + 7/10 I, prec], N[-1/2 + 1/5 I, prec], N[2, prec]};
Print["omega_{1,1}: TR vs model:"];
Do[Print[N[om11TR[p0], 15], "   ", N[pred11[p0], 15]], {p0, pts}];
(* omega_{0,3} *)
om03TR[p0_, p1_, p2_] := Sum[Module[{c = cs[[b]], s = iotas[[b]], q, iq, integrand},
    q = c + e; iq = c + s;
    integrand = 1/2 (1/(p0 - q) - 1/(p0 - iq))/((q - iq) (D[xf[w], w] /. w -> q)) *
      (D[s, e]/((q - p1)^2 (iq - p2)^2) + D[s, e]/((q - p2)^2 (iq - p1)^2));
    SeriesCoefficient[Series[integrand, {e, 0, 1}], -1]], {b, 1, 5}];
(* genus-0 3-pt of model in flat basis at v=1: from model.wl output, <H^aH^bH^c>_0 = 5 * sign *)
Get["g03.wl"];
pred03[p0_, p1_, p2_] := (-1/5) Sum[g03[[a + 1, b + 1, cc + 1]] (D[xis[[a + 1]], p] /. p -> p0) (D[xis[[b + 1]], p] /. p -> p1) (D[xis[[cc + 1]], p] /. p -> p2), {a, 0, 4}, {b, 0, 4}, {cc, 0, 4}];
Print["omega_{0,3}: TR vs model:"];
Print[N[om03TR[pts[[1]], pts[[2]], pts[[3]]], 15], "   ", N[pred03[pts[[1]], pts[[2]], pts[[3]]], 15]];
