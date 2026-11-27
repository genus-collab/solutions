zeta = Exp[2 Pi I/5];
w0 = Join[Table[1 - zeta^j, {j, 1, 4}], {-5}];
bk[k_] := BernoulliB[2 k]/(2 k (2 k - 1));
pk[k_] := RootReduce[Total[w0^(1 - 2 k)]];
Print[Table[FullSimplify[pk[k]], {k, 1, 3}]];
h = Exp[-Sum[bk[k] z^(2 k - 1) pk[k], {k, 1, 3}]];
Print[Series[h, {z, 0, 5}] // FullSimplify];
