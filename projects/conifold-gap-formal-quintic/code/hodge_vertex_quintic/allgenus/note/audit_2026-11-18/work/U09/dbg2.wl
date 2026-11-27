Get["rtr.wl"];
cR[j_, i_, k_] := Coefficient[Rinv[[j, i]] /. th -> 1/s, s, k];
cM[i_, j_, k_] := (-1)^k Coefficient[Mt[[i, j]], t, k];
Print["max |R^{-1 j}_i - M0_{ij}| coeffs k=0..4: ", Table[Max[Table[Abs[SetPrecision[cR[j, i, k] - cM[i, j, k], MachinePrecision]], {j, 5}, {i, 5}]], {k, 0, 3}]];
Print["max |R^{-1 j}_i - M0_{ji}| coeffs k=0..4 (transposed reading): ", Table[Max[Table[Abs[SetPrecision[cR[j, i, k] - cM[j, i, k], MachinePrecision]], {j, 5}, {i, 5}]], {k, 0, 3}]];
lhs2 = Table[Expand[Sum[Rinv[[i, k]]/(sq cs[[k]]), {k, 5}] /. th -> 1/s], {i, 5}];
rhs2 = Table[Expand[Lap[i, D[pws[[i]], w]] /. th -> 1/s], {i, 5}];
Print["(2-32) max discrepancy k=0..4: ", Table[Max[Table[Abs[SetPrecision[Coefficient[lhs2[[i]] - rhs2[[i]], s, k], MachinePrecision]], {i, 5}]], {k, 0, 3}]];
