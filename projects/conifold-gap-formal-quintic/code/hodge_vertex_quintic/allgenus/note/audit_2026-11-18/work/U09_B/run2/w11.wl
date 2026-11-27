(* omega_{1,1} of x = p^4/4 + 1/p, y = p (v = 1) from the recursion with the kernel of conv:mod-kernel,
   versus the genus-1 one-point graph sum of M0(z).omega^m (unit-preserving action, PPZ sign of the dilaton leaf). *)
Get["m0.mx"];
prec = 50; Nw = 8;
x[p_] := p^4/4 + 1/p;
cs = Table[N[z5^(-j), prec], {j, 0, 4}];
sq = Sqrt[-5];
pser[c_] := Module[{a, pw, eq, s},
  pw = c + w/(sq c) + Sum[a[k] w^k, {k, 2, Nw}];
  s = {};
  Do[eq = Chop[SeriesCoefficient[Series[x[pw //. s] - x[c] + w^2/2, {w, 0, k + 1}], k + 1], 10^-40];
     s = Join[s, {a[k] -> (a[k] /. Solve[eq == 0, a[k]][[1]])}], {k, 2, Nw}];
  Chop[(pw //. s) /. a[_] -> 0, 10^-40]];
(* TR: omega11(p0)/dp0 = sum_b Res_{w=0} (1/2)[1/(p0-p(w)) - 1/(p0-p(-w))] p'(w) p'(-w) / (w (p(w)-p(-w))^3) *)
p0 = N[3/10 + 7 I/10, prec];
tr = Sum[Module[{pw = pser[c], pm, f},
    pm = pw /. w -> -w;
    f = (1/2) (1/(p0 - pw) - 1/(p0 - pm)) D[pw, w] (D[pw, w] /. w -> -w) /(w (pw - pm)^3);
    SeriesCoefficient[Series[f, {w, 0, 1}], -1]], {c, cs}];
(* the factor D[pw,w] /. w->-w equals p'(-w) since pm' = -p'(-w): B(q, iota q) = p'(w) d(p(-w)) /(...)^2, d(p(-w)) = -p'(-w) dw;
   K has dx = -w dw; the two signs cancel, as in the formula above *)
Print["TR omega_{1,1}(p0)/dp0 = ", Chop[N[tr, 15]]];
(* Model side *)
R1 = Map[Coefficient[#, t, 1] &, M0, {2}];
sig = cs; (* varsigma_i = zeta^{-i} at v = 1, i = 0..4 *)
corr0 = Table[-(1/24) Sum[R1[[j, i]] sig[[j]], {j, 1, 5}] + (1/2) R1[[i, i]] sig[[i]] + (1/24) Sum[R1[[i, k]] sig[[i]]^2/sig[[k]], {k, 1, 5}], {i, 1, 5}];
corr0wrong = Table[-(1/24) Sum[R1[[j, i]] sig[[j]], {j, 1, 5}] + (1/2) R1[[i, i]] sig[[i]] - (1/24) Sum[R1[[i, k]] sig[[i]]^2/sig[[k]], {k, 1, 5}], {i, 1, 5}];
corr1 = Table[sig[[i]]/24, {i, 1, 5}];
xit[i_, p_] := -1/(sig[[i]] (p - sig[[i]]));
Dx[f_] := p^2/(p^5 - 1) D[f, p];
model[c0_] := Sum[c0[[i]] D[xit[i, p], p] + corr1[[i]] D[Dx[xit[i, p]], p], {i, 1, 5}] /. p -> p0;
Print["model omega_{1,1}(p0)/dp0 (PPZ dilaton sign) = ", Chop[N[model[corr0], 15]]];
Print["model with the opposite dilaton sign         = ", Chop[N[model[corr0wrong], 15]]];
