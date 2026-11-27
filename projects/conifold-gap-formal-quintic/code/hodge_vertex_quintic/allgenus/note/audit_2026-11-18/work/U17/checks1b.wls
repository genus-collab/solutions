bk[k_] := BernoulliB[2k]/(2k(2k-1));
beta[n_, th_] := (-1)^(n+1) Expand[BernoulliB[n+1, th]]/(n (n+1));
lhs4 = Series[(1/v + th - 1/2) Log[1 + th v] - th + Sum[bk[k] (v/(1 + th v))^(2k-1), {k, 1, 8}], {v, 0, 12}];
rhs4 = Sum[beta[n, th] v^n, {n, 1, 12}];
Print["DLMF check: ", Expand[Normal[lhs4] - rhs4]];
Print["beta_1(th) = ", beta[1,th], "  beta_{2k-1}(0)=bk: ", Table[beta[2k-1,0]-bk[k],{k,1,5}], " beta_{2k}(0): ", Table[beta[2k,0],{k,1,5}]];
