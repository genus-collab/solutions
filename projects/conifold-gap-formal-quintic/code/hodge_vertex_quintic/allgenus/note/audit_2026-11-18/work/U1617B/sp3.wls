K = 9;
bk[k_] := BernoulliB[2 k]/(2 k (2 k - 1));
(* dlmf5118 in variable e = 1/w *)
lhs = (1/e + th - 1/2) Log[1 + th e] - th + Sum[bk[k] e^(2 k - 1) (1 + th e)^(1 - 2 k), {k, 1, K}];
rhs = Sum[(-1)^(n + 1) BernoulliB[n + 1, th]/(n (n + 1)) e^n, {n, 1, K}];
Print["dlmf5118: ", Simplify[Normal[Series[lhs - rhs, {e, 0, K}]]]];
(* stirfe rho(Y) = 0, variable e = 1/Y *)
mu[e_] := Sum[bk[k] e^(2 k - 1), {k, 1, K}];
rho = (1/e + 1/2) Log[1 + e] - 1 + mu[e/(1 + e)] - mu[e];
Print["stirfe rho: ", Simplify[Normal[Series[rho, {e, 0, K}]]]];
(* St_m recursion *)
St[m_, w_] := Exp[(1/w + m - 1/2) Log[1 + m w] - m + Sum[bk[k] w^(2 k - 1) (1 + m w)^(1 - 2 k), {k, 1, K}]];
Do[Print["St rec m=", m, ": ", Simplify[Normal[Series[St[m + 1, w] - (1 + m w) St[m, w], {w, 0, 7}]]]], {m, {-3, -1, 0, 2}}];
(* Ginf first coefficient *)
Print["B6(1/5)/6 = ", BernoulliB[6, 1/5]/6];
(* Step 2 equation of the Stirling lemma: s0(w) = exp((1/w+1/2)log(1+w)-1) s0(w/(1+w)) with s0 = St_0 *)
s0[w_] := Exp[Sum[bk[k] w^(2 k - 1), {k, 1, K}]];
Print["step2: ", Simplify[Normal[Series[s0[w] - Exp[(1/w + 1/2) Log[1 + w] - 1] s0[w/(1 + w)], {w, 0, 8}]]]];
