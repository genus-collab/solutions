(* U17 sanity checks: Stirling functional equation, St_m, stationary phase of a log x - x, DLMF 5.11.8 form *)
bk[k_] := BernoulliB[2k]/(2k(2k-1));
NN = 14;
(* 1. rho(Y) = 0 in Q[[1/Y]]: use y = 1/Y *)
mu[yy_, K_] := Sum[bk[k] yy^(2k-1), {k, 1, K}];   (* mu as function of 1/Y *)
rho = Series[(1/y + 1/2) Log[1 + y] - 1 + mu[y/(1+y), NN] - mu[y, NN], {y, 0, 2NN-2}];
Print["rho series (should be O): ", rho];

(* 2. St_m(w) vs Poch_m(w) St_0(w) *)
St[m_, w_, K_] := Exp[(1/w + m - 1/2) Log[1 + m w] - m + Sum[bk[k] w^(2k-1) (1 + m w)^(1-2k), {k, 1, K}]];
Poch[M_, w_] := If[M >= 0, Product[1 + i w, {i, 0, M-1}], Product[1/(1 - i w), {i, 1, -M}]];
Do[
  d = Series[St[m, w, 8] - Poch[m, w] St[0, w, 8], {w, 0, 12}];
  Print["m=", m, " St_m - Poch_m St_0: ", d], {m, -3, 3}];

(* 3. SPn_a[a log x - x; (x/a)^(m-1)] by Wick, compare with St_m(z/a) *)
spnStir[m_, ord_] := Module[{s, tau, ser, res, F3},
  (* F(a+t)-F(a) = a log(1+t/a) - t; H = -1/a; C = a *)
  ser = Normal[Series[Exp[(a Log[1 + s tau/a] - s tau + (s tau)^2/(2 a))/s^2] (1 + s tau/a)^(m-1), {s, 0, 2 ord}]];
  ser = Expand[ser];
  res = ser /. tau^n_. :> If[OddQ[n], 0, (n-1)!! a^(n/2)];
  res = res /. tau -> 0;
  Expand[res /. s -> Sqrt[z]]];
Do[
  lhs = spnStir[m, 5];
  rhs = Normal[Series[St[m, z/a, 6], {z, 0, 5}]];
  Print["m=", m, " SPn - St_m(z/a) = ", Simplify[lhs - rhs]], {m, -2, 3}];

(* 4. DLMF 5.11.8 formal: in Q[theta][[1/w]], with v=1/w *)
beta[n_, th_] := (-1)^(n+1) BernoulliPolynomial[n+1, th]/(n (n+1));
lhs4 = Series[(1/v + th - 1/2) Log[1 + th v] - th + Sum[bk[k] (v/(1 + th v))^(2k-1), {k, 1, 8}], {v, 0, 12}];
rhs4 = Sum[beta[n, th] v^n, {n, 1, 12}];
Print["DLMF check: ", Simplify[Normal[lhs4] - rhs4]];

(* 5. reflection: mu(-y) = -mu(y) trivial; check Gst(1-y) formula via rho at Y=-y *)
(* log Gst(1-y) := (1/2-y)(log(-y)+log(1-1/y)) - (1-y) + 1/2 log 2pi + mu(1-y);
   claim = (1/2-y) log(-y) + y + 1/2 log 2pi + mu(-y). Difference in terms of v=1/y: *)
diff5 = Series[(1/2 - 1/v) Log[1 - v] - 1 + mu[-v/(1 - v), NN] - mu[-v, NN], {v, 0, 2NN-2}];
Print["refl auxiliary (should be O): ", diff5];
