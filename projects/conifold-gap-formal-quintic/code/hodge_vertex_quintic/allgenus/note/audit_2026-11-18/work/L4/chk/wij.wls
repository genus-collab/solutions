(* Check (7.7): (exp(h^2 wi wj S(wi h d_i) S(wj h d_j) (zi-zj)^-2) - 1)/h^2 = sum_k (h^2/4)^k (wi+wj)^(2k) wi wj (zi-zj)^(-2k-2) *)
S[u_] := Sinh[u/2]/(u/2);
N0 = 4;
(* S(w h d) acting on f: sum_k (w h)^(2k) d^(2k) f/(4^k (2k+1)!) *)
Sop[w_, var_, f_] := Sum[(w hb)^(2 k) D[f, {var, 2 k}]/(4^k (2 k + 1)!), {k, 0, N0}];
A = hb^2 wi wj Sop[wj, zj, Sop[wi, zi, (zi - zj)^-2]];
lhs = Normal[Series[(Exp[A] - 1)/hb^2, {hb, 0, 2 N0 - 2}]];
rhs = Sum[(hb^2/4)^k (wi + wj)^(2 k) wi wj (zi - zj)^(-2 k - 2), {k, 0, N0 - 1}];
Print[Simplify[Normal[Series[lhs - rhs, {hb, 0, 2 N0 - 2}]]]];
