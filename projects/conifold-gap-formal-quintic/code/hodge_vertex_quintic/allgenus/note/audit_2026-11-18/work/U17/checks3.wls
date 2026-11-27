bk[k_] := BernoulliB[2k]/(2k(2k-1));
s0[w_] := Exp[Sum[bk[k] w^(2k-1), {k, 1, 8}]];
Print["Step2 functional eq: ", Series[s0[w] - Exp[(1/w + 1/2) Log[1 + w] - 1] s0[w/(1 + w)], {w, 0, 14}]];
ord = 6; (* s-order in SPn *)
spnA[ampT_] := Module[{s, tau, ser},
  ser = Normal[Series[Exp[(a Log[1 + s tau/a] - s tau + (s tau)^2/(2 a))/s^2] (ampT /. t -> s tau), {s, 0, 2 ord}]];
  ser = Expand[ser] /. tau^n_. :> If[OddQ[n], 0, (n-1)!! a^(n/2)];
  Expand[(ser /. tau -> 0) /. s -> Sqrt[z]]];
direct = Sum[e^d/d! z^(-d) spnA[Log[1 + t/a]^d (1 + t/a)^(-1)], {d, 0, 4}];
closed = Expand[Normal[Series[Exp[((a + e) Log[1 + e/a] - e)/z] (1 + e/a)^(-1/2) s0[z/(a + e)], {e, 0, 4}]]];
res = Table[Simplify[Coefficient[Expand[Normal[Series[Coefficient[closed, e, d], {z, 0, 2}]] - Coefficient[direct, e, d]], z, k]], {d, 0, 4}, {k, -2, 2}];
Print["Step2 S_eps check (zeros expected): ", res];
