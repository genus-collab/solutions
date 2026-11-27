(* review_adv_lemma2.wls -- adversarial verification of the string-equation proof of Z[h,n] = 0 (h >= 2, n >= 1).
   Independent enumeration code (ordered tuples, pruned by the support of each leaf series). Read-only on the project files.
   Checks: (1) the string equation (F2) for the Hodge correlators hodgeNormU4 (random tuples, genus 1-3);
           (2) Lemma 2 term by term, I_m(1) = (n/c) J_m - m I_{m-1}(z phi), for genus 1-2, n = 0..2, all m, with the CY weights and c = 1+x,
               and with a symbolic c unrelated to the weights;
           (3) the by-hand table for (h,n) = (2,1); (4) the recursion c^2 Z[h,n+1] = n Z[h,n] and the genus-1 anomaly;
           (5) Z[2,n] and Z[1,n] for the CY triple and for non-CY numeric weights / arbitrary c; (6) genus 3, n = 0,1 if time permits. *)
SetDirectory[DirectoryName[$InputFileName]];
t0 = AbsoluteTime[];
Get["hvc_lib.m"];
Print["library loaded: ", Round[AbsoluteTime[] - t0], " s"];
el[] := Round[AbsoluteTime[] - t0];

(* correlator Int_{M_{h,N}} Omega prod_i f_i(psi_i); legs = list of coefficient functions k -> coefficient of psi^k (0 outside the support).
   Omega is of mixed degree, so every tuple with Total[a] <= dim contributes (hodgeNormU4 picks the right degree component of Omega). *)
corrLE[h_, legs_List, usub_] := Module[{N = Length[legs], dim, tot = 0, rec},
   dim = 3 h - 3 + N;
   If[2 h - 2 + N <= 0, Return[0]];
   rec[i_, prefix_, budget_, coef_] := If[i > N, tot += coef (hodgeNormU4[h, prefix] /. usub),
      Do[Module[{cf = legs[[i]][k]}, If[cf =!= 0, rec[i + 1, Append[prefix, k], budget - k, coef cf]]], {k, 0, budget}]];
   rec[1, {}, dim, 1];
   Together[tot]];

(* leaf series with parameter cc *)
phiF[cc_] := Function[k, 1/cc^(k + 1)];
TF[cc_] := Function[k, If[k >= 2, -1/cc^(k - 1), 0]];
oneF = Function[k, If[k == 0, 1, 0]];
zphiF[cc_] := Function[k, If[k >= 1, 1/cc^k, 0]];
(* J_m, I_m(f) and Z[h,n] *)
Jm[h_, n_, m_, cc_, usub_] := corrLE[h, Join[ConstantArray[phiF[cc], n], ConstantArray[TF[cc], m]], usub];
Ileaf[h_, n_, m_, f_, cc_, usub_] := corrLE[h, Join[ConstantArray[phiF[cc], n], {f}, ConstantArray[TF[cc], m]], usub];
Zf[h_, n_, cc_, usub_] := Together[Sum[Jm[h, n, m, cc, usub]/m!, {m, 0, 3 h - 3 + n}]];

usCY = usubOf[{1 + x, -x, -1}];       (* the task's weights, e2 = 1, e1 = x *)
usN = usubOf[{2, 5, -7}];             (* three non-CY numeric weights *)

(* ---------- (1) string equation for the Hodge correlators, independent random check ---------- *)
SeedRandom[7];
strOK = And @@ Table[Module[{g = RandomInteger[{1, 3}], n, a},
    n = RandomInteger[{1, If[g == 3, 3, 4]}]; a = RandomInteger[{0, 3 g - 3 + n}, n];
    Together[hodgeNormU4[g, Append[a, 0]] - Sum[If[a[[i]] >= 1, hodgeNormU4[g, ReplacePart[a, i -> a[[i]] - 1]], 0], {i, n}]] === 0], {30}];
Print["(1) string equation Int Omega psi^a psi_{n+1}^0 = Sum_i Int Omega psi^{a - e_i}, 30 random (g, a), genus 1-3, symbolic u: ", strOK, "  [", el[], " s]"];

(* ---------- (2) Lemma 2 term by term ---------- *)
lemma2[h_, n_, cc_, usub_] := Table[Module[{lhs, rhs},
    lhs = Ileaf[h, n, m, oneF, cc, usub];
    rhs = (n/cc) Jm[h, n, m, cc, usub] - If[m >= 1, m Ileaf[h, n, m - 1, zphiF[cc], cc, usub], 0];
    {m, Together[lhs - rhs] === 0, lhs =!= 0 || rhs =!= 0}], {m, 0, 3 h - 2 + n}];
Do[Print["(2) Lemma 2, h = ", h, ", n = ", n, ", CY weights, c = 1+x: {m, identity holds, nontrivial} = ", lemma2[h, n, 1 + x, usCY], "  [", el[], " s]"], {h, 1, 2}, {n, If[h == 1, 1, 0], 2}];
Print["(2') Lemma 2 with symbolic c unrelated to the weights, h = 2, n = 1: ", lemma2[2, 1, cc, usCY], "  [", el[], " s]"];
Print["(2'') Lemma 2 fails only at (h,n,m) = (1,0,0) as the proof says: ", lemma2[1, 0, 1 + x, usCY], "  [", el[], " s]"];

(* ---------- (3) by-hand table (h,n) = (2,1): I_m(1), I_m(z phi) ---------- *)
Print["(3) h=2, n=1 (so Lemma 2 has n = 0: the proof's table): I_m(1) = ", Table[Factor[Ileaf[2, 0, m, oneF, 1 + x, usCY]], {m, 0, 5}]];
Print["    I_m(z phi) = ", Table[Factor[Ileaf[2, 0, m, zphiF[1 + x], 1 + x, usCY]], {m, 0, 4}]];
Print["    expected (proof): I_0(z phi) = 7/(5760 x^2); I_1(1) = -7/(5760 x^2); I_2(z phi)/2 = 5(1+3x+x^2)/(288 x (1+x)^4); I_3(z phi)/6 = -7/(288 (1+x)^4); I_4(1) = 7/(12 (1+x)^4)"];
AA = Together[Sum[Ileaf[2, 0, m, oneF, 1 + x, usCY]/m!, {m, 0, 5}]]; BB = Together[Sum[Ileaf[2, 0, m, zphiF[1 + x], 1 + x, usCY]/m!, {m, 0, 4}]];
Print["    A = ", Factor[AA], ", B = ", Factor[BB], ", A + B = ", Together[AA + BB], " (expected A = (7+10x^2+7x^4)/(2880 x^2 (1+x)^4) = -B)  [", el[], " s]"];

(* ---------- (4), (5) recursion, anomaly, vanishing ---------- *)
Z2 = Table[Zf[2, n, 1 + x, usCY], {n, 0, 3}];
Print["(4) CY triple, c = 1+x: Z[2,0..3] = ", Factor /@ Z2, ";  (-x c^3) Z[2,0] - gap_2 = ", Together[-x (1 + x)^3 Z2[[1]] - gapHVC[2, x][[1]]], "  [", el[], " s]"];
Z1 = Table[Zf[1, n, 1 + x, usCY], {n, 1, 4}];
Print["    genus 1: Z[1,1..4] = ", Factor /@ Z1, ";  c^2 Z[1,n+1] - n Z[1,n] = ", Table[Together[(1 + x)^2 Z1[[n + 1]] - n Z1[[n]]], {n, 1, 3}],
   ";  c Z[1,1] - Int_{M_{1,1}} Omega = ", Together[(1 + x) Z1[[1]] - (hodgeNormU4[1, {0}] /. usCY)], ";  Int_{M_{1,1}} Omega = -(1/24) Sum 1/w_a: ", Together[(hodgeNormU4[1, {0}] /. usCY) + (1/24) (1/(1 + x) - 1/x - 1)] === 0];
Print["    quoted values: Z[1,1] = (1+x+x^2)/(24 x (1+x)^2): ", Together[Z1[[1]] - (1 + x + x^2)/(24 x (1 + x)^2)] === 0, ", Z[1,2]: ", Together[Z1[[2]] - (1 + x + x^2)/(24 x (1 + x)^4)] === 0, ", Z[1,3]: ", Together[Z1[[3]] - (1 + x + x^2)/(12 x (1 + x)^6)] === 0, "  [", el[], " s]"];
Print["(5) generality: Z[2,1], Z[2,2] with symbolic c unrelated to the CY weights: ", {Zf[2, 1, cc, usCY], Zf[2, 2, cc, usCY]},
   ";  non-CY weights (2,5,-7), c = 3: Z[2,1], Z[2,2] = ", {Zf[2, 1, 3, usN], Zf[2, 2, 3, usN]}, ", Z[2,0] = ", Zf[2, 0, 3, usN],
   ";  genus 1, non-CY, c = 3: c^2 Z[1,2] - Z[1,1] = ", Together[9 Zf[1, 2, 3, usN] - Zf[1, 1, 3, usN]], ", c Z[1,1] - Int Omega = ", Together[3 Zf[1, 1, 3, usN] - (hodgeNormU4[1, {0}] /. usN)], "  [", el[], " s]"];
Print["    Witten-Kontsevich vertex (u -> 0), c = 1+x: Z[1,1], Z[1,2], Z[2,1], Z[2,2] = ", {Zf[1, 1, 1 + x, {u[_] -> 0}], Zf[1, 2, 1 + x, {u[_] -> 0}], Zf[2, 1, 1 + x, {u[_] -> 0}], Zf[2, 2, 1 + x, {u[_] -> 0}]}, "  [", el[], " s]"];

(* ---------- (6) genus 3 if time permits ---------- *)
If[el[] < 90,
  Print["(6) genus 3, Lemma 2 with n = 0 (all m): ", lemma2[3, 0, 1 + x, usCY], "  [", el[], " s]"];
  Z31 = Zf[3, 1, 1 + x, usCY]; Print["    Z[3,1] = ", Z31, "  [", el[], " s]"];
  If[el[] < 130, Print["    Lemma 2 with n = 1, genus 3 (m <= 4 only): ", Table[Module[{lhs, rhs}, lhs = Ileaf[3, 1, m, oneF, 1 + x, usCY]; rhs = (1/(1 + x)) Jm[3, 1, m, 1 + x, usCY] - If[m >= 1, m Ileaf[3, 1, m - 1, zphiF[1 + x], 1 + x, usCY], 0]; {m, Together[lhs - rhs] === 0}], {m, 0, 4}], "  [", el[], " s]"]],
  Print["(6) skipped (time)"]];
Print["total ", el[], " s"];
