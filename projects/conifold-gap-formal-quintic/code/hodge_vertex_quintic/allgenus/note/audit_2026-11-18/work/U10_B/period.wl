(* Check prop:mod-period Steps 1-4 for the formula as transcribed *)
X[p_] := p^4/4 + v/p;
eps[k_] := -D[X[p], {p, 2 k}]/(4^k (2 k + 1)!);
xp = D[X[p], p];
Dx[f_] := Together[D[f, p]/xp];
(* a_r = -[hbar^{2h}][w^{r+1}] e^E, E = Sum_k w^{2k+1} hb^{2k} eps_k *)
PhiA[h_] := Module[{E0, ser, tot = 0},
  E0 = Sum[w^(2 k + 1) hb^(2 k) eps[k], {k, 1, h}];
  ser = Normal[Series[Exp[E0], {hb, 0, 2 h}]];
  ser = Expand[Coefficient[ser, hb, 2 h]];
  Do[With[{ar = -Coefficient[ser, w, r + 1]},
     If[ar =!= 0, tot += Nest[Dx, Together[ar/xp], r - 1]]], {r, 1, 3 h}];
  Together[tot]];
Do[Module[{f = PhiA[h], den, at0, atinf},
   den = Factor[Denominator[f]];
   at0 = Series[f, {p, 0, 1}];
   atinf = Limit[f, p -> Infinity];
   Print["h=", h, " denominator=", den, "  Phi(0)=", Simplify[SeriesCoefficient[at0, 0]],
     " expected=", Simplify[-BernoulliB[2 h]/(2 h) v^(1 - 2 h)], " ord0>=0: ",
     SeriesCoefficient[Series[f, {p, 0, -1}], -1] === 0, " Phi(inf)=", atinf]], {h, 1, 4}];
(* Step 4 identity *)
aa[x_] := Sum[x^k/(2 k + 1), {k, 1, 20}];
Do[Module[{S1, S2},
   S1 = Sum[(-1)^j (2 h - 2 + j)!/j! SeriesCoefficient[aa[x]^j, {x, 0, h}], {j, 1, h}];
   S2 = (2 h - 2)! SeriesCoefficient[(1 + aa[x])^(1 - 2 h), {x, 0, h}];
   Print["h=", h, " Sigma=", S1, " binom=", S2, " target=", -4^h BernoulliB[2 h]/(2 h)]], {h, 1, 8}];
(* Airy check of the transcribed formula: x=z^2/2 *)
