(* chain check: two-variable SP (Euler integral) = one-variable SP with Stirling amplitude (eq gam-step1) *)
K = 3;
mom2[i_, j_, C_] := Module[{x1, x2}, D[Exp[({x1, x2}.C.{x1, x2})/2], {x1, i}, {x2, j}] /. {x1 -> 0, x2 -> 0}];
sp2[P3_, amp_, C_, K_] := Module[{poly, res},
  poly = Expand[Normal[Series[Exp[P3[h t1, h t2]/h^2] amp[h t1, h t2, h^2], {h, 0, 2 K}]]];
  res = poly /. {t1^i_. t2^j_. :> mom2[i, j, C], t1^i_. :> mom2[i, 0, C], t2^j_. :> mom2[0, j, C]};
  Expand[Normal[Series[res, {h, 0, 2 K}]] /. h -> Sqrt[z]]];
sp1[F3_, amp_, C_, K_] := Module[{poly, res},
  poly = Expand[Normal[Series[Exp[F3[h tau]/h^2] amp[h tau, h^2], {h, 0, 2 K}]]];
  res = poly /. tau^n_. :> If[EvenQ[n], (n - 1)!! C^(n/2), 0]; res = res /. tau -> 0;
  Expand[Normal[Series[res, {h, 0, 2 K}]] /. h -> Sqrt[z]]];
bk[k_] := BernoulliB[2 k]/(2 k (2 k - 1));
poch[M_, w_] := If[M >= 0, Product[1 + i w, {i, 0, M - 1}], Product[1/(1 - i w), {i, 1, -M}]];
as = 2; aL = 3; c2 = 1/5; c3 = -2/7; c4 = 1/3;
phit[s_] := c2 s^2 + c3 s^3 + c4 s^4;   (* phi(L+s)-phi(L)-phi'(L)s *)
Fz[s_, zz_] := 1 + zz s + zz^2 s^2/2 - s/3;
Do[
 (* two-variable: variables (s, m) *)
 Hm = {{2 c2, as}, {as, -aL}};
 C2 = -Inverse[Hm];
 P3two[s_, m_] := phit[s] - c2 s^2 - aL (Exp[m] - 1 - m - m^2/2);
 amp2[s_, m_, zz_] := Exp[Ms m] Fz[s, zz];
 lhs = sp2[P3two, amp2, C2, K];
 lhs = Expand[Normal[Series[lhs/(lhs /. z -> 0), {z, 0, K}]]];
 (* one-variable *)
 Sfull[s_] := phit[s] + (aL + as s) Log[1 + as s/aL] - as s;
 H1 = 2 c2 + as^2/aL; C1 = -1/H1;
 F3one[s_] := Normal[Series[Sfull[x], {x, 0, 12}]] /. x -> s // (# - H1 s^2/2 &);
 amp1[s_, zz_] := (1 + as s/aL)^(Ms - 1/2) poch[Ms, zz/(aL + as s)] Exp[Sum[bk[k] (zz/(aL + as s))^(2 k - 1), {k, 1, K}]] Fz[s, zz];
 rhs = sp1[F3one, amp1, C1, K];
 rhs = Expand[Normal[Series[rhs/(rhs /. z -> 0), {z, 0, K}]]];
 Print["M*=", Ms, " diff: ", Simplify[lhs - rhs], "   lhs z^1: ", Coefficient[lhs, z, 1]],
 {Ms, {0, 1, 2}}];
