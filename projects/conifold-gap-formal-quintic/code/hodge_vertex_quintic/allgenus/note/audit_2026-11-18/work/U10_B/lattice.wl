(* Check thm:mod-lattice (d) and (c) dimensions in small weights *)
X[p_] := p^4/4 + v/p;
antider[F_] := Module[{e = Expand[F], terms},
  If[e === 0, 0, terms = If[Head[e] === Plus, List @@ e, {e}];
  Total[Map[Integrate[#, p] &, terms]]]];
Iop[F_] := Expand[antider[Expand[-F (p^3 - v p^-2)]]];
phi0 = {p, Log[p], v/p, p^3, p^2};
rho = {-1, 1, -1, -3, -2};
wb = {1, 0, 4, 3, 2};
eps = {1, -1, -1, 1, 1};
K = 9;
legs = Table[NestList[Iop, phi0[[b + 1]], K], {b, 0, 4}];
Pi0[m_, b_] := If[Mod[b + m, 5] == 1, With[{r = Floor[m/5]}, -eps[[b + 1]] v^(m - r)/(4^r r! (m - r)!)], 0];
Lam[{m_, k_, b_}] := v^m legs[[b + 1, k + 1]]/rho[[b + 1]];
basisV[w_] := Select[Flatten[Table[{m, k, b}, {m, 0, Floor[w/5]}, {k, 0, Floor[w/4]}, {b, 0, 4}], 2], 5 #[[1]] + 4 #[[2]] + wb[[#[[3]] + 1]] == w &];
(* coordinates: non-log part as Laurent poly in p (v determined by weight), log part also *)
coords[F_, w_] := Module[{e = Expand[F /. Log[p] -> LL], A, B},
  A = Coefficient[e, LL, 0]; B = Coefficient[e, LL, 1];
  Join[Table[Coefficient[A /. v -> 1, p, w - 5 m], {m, 0, Floor[w/4]}],
       Table[Coefficient[B /. v -> 1, p, w - 5 m], {m, 0, Floor[w/4]}]]];
(* Q[v,x]_w elements (non-log) *)
qvx[w_] := Select[Flatten[Table[{m, n}, {m, 0, Floor[w/5]}, {n, 0, Floor[w/4]}], 1], 5 #[[1]] + 4 #[[2]] == w &];
res = Table[Module[{bV = basisV[w], rows, extra, dimP, e, d},
   rows = Map[coords[Lam[#], w] &, bV];
   extra = Map[coords[Expand[v^#[[1]] X[p]^#[[2]]], w] &, qvx[w]];
   e = Length[qvx[w]]; d = Length[Select[bV, #[[3]] != 1 &]];
   {w, Length[bV], d, e, MatrixRank[Join[rows, extra]], Floor[w/4] + 1 + e,
    (* K0 dimension: kernel of all Pi_i *)
    Length[bV] - MatrixRank[Table[(Expand[If[#[[2]] >= i, v^#[[1]] Pi0[#[[2]] - i, #[[3]]], 0]] /. v -> 1) & /@ bV, {i, 0, Floor[w/4]}]]}],
  {w, 0, 36}];
Print[TableForm[res, TableHeadings -> {None, {"w", "dimV", "d", "e", "rank[Lam;Qvx]", "dimQ[p,pb]_w+e", "dimK0"}}]];
Print["all full rank: ", And @@ Map[#[[5]] == #[[6]] && #[[5]] == #[[2]] + #[[4]] &, res]];
Print["K0 dim = d: ", And @@ Map[#[[7]] == #[[3]] &, res]];
