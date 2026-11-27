Get["rtr.wl"];
cR[j_, i_, k_] := Coefficient[Rinv[[j, i]] /. th -> 1/s, s, k];
cM[i_, j_, k_] := (-1)^k Coefficient[Mt[[i, j]], t, k];
vals = Table[{j, i, k, Chop[N[cR[j, i, k], 12], 10^-10], Chop[N[cM[i, j, k], 12], 10^-10], Chop[N[cM[j, i, k], 12], 10^-10]}, {k, 2, 3}, {j, 1, 5}, {i, 1, 5}];
dd = Flatten[vals, 2];
Print["max |R - M_ij| k=2,3: ", Max[Abs[Chop[#[[4]] - #[[5]], 10^-9] & /@ dd]]];
Print["max |R - M_ji| k=2,3: ", Max[Abs[Chop[#[[4]] - #[[6]], 10^-9] & /@ dd]]];
Print[Take[dd, 4]];
