(* review_verifier_A.wls -- independent verifier checks of the string-equation proof of Z[h,n] = 0 (h >= 2, n >= 1).
   Own enumeration (generic typed correlator), tests Lemma 2 TERM BY TERM in m, the recursion cc^2 Z[h,n+1] = n Z[h,n],
   the genus-1 anomaly, the trivial CohFT, the quoted values and the quoted h = 2 table. Read-only on the project files. *)
SetDirectory[DirectoryName[$InputFileName]];
tt = AbsoluteTime[];
Get["hvc_lib.m"];
Print["library loaded in ", Round[AbsoluteTime[] - tt], " s"];
lap[] := (Print["   [t = ", Round[AbsoluteTime[] - tt], " s]"]);

(* multisets (sorted lists) of length k, entries >= lo, total <= smax *)
ms[k_, lo_, smax_] := Module[{rec}, rec[kk_, l_, s_] := If[kk == 0, {{}}, Join @@ Table[Prepend[#, a] & /@ rec[kk - 1, a, s - a], {a, l, s}]]; rec[k, lo, smax]];
nord[lst_] := Length[lst]!/Times @@ (Last /@ Tally[lst]!);      (* number of distinct orderings of a multiset *)

(* Corr[h, V, types]: Int_{M_{h,N}} (class) prod_points f_i(psi_i); types = {{f, number of points, minimal exponent}, ...},
   f[k] = psi^k coefficient of the leaf series; V[h, a] = Int_{M_{h,Length[a]}} (class) psi^a (zero unless degrees match) *)
Corr[h_, V_, types_List] := Module[{N = Total[types[[All, 2]]], dim, tot = 0, rec},
   If[2 h - 2 + N <= 0, Return[0]];
   dim = 3 h - 3 + N;
   rec[{}, acc_, budget_] := (tot += (Times @@ (Last /@ acc)) V[h, Join @@ (First /@ acc)]);
   rec[{{f_, k_, lo_}, rest___}, acc_, budget_] := Do[rec[{rest}, Append[acc, {A, nord[A] (Times @@ (f /@ A))}], budget - Total[A]], {A, ms[k, lo, budget]}];
   rec[types, {}, dim];
   Together[tot]];

phiF[cc_] := Function[k, cc^(-k - 1)];                       (* phi = 1/(cc - psi) *)
TF[cc_] := Function[d, -cc^(1 - d)];                         (* T = -psi^2/(cc - psi), d >= 2 *)
oneF = Function[k, If[k == 0, 1, 0]];                        (* constant leaf 1 *)
zphiF[cc_] := Function[k, If[k >= 1, cc^(-k), 0]];           (* psi phi(psi), k >= 1 *)

J[h_, n_, m_, V_, cc_] := Corr[h, V, {{phiF[cc], n, 0}, {TF[cc], m, 2}}];
I1[h_, n_, m_, V_, cc_] := Corr[h, V, {{phiF[cc], n, 0}, {oneF, 1, 0}, {TF[cc], m, 2}}];
Izp[h_, n_, m_, V_, cc_] := If[m < 0, 0, Corr[h, V, {{phiF[cc], n, 0}, {zphiF[cc], 1, 1}, {TF[cc], m, 2}}]];
Z[h_, n_, V_, cc_] := Z[h, n, V, cc] = Together[Sum[J[h, n, m, V, cc]/m!, {m, 0, 3 h - 3 + n}]];

(* vertex classes *)
Vgen = Function[{h, a}, hodgeNormU4[h, a]];                                  (* generic weights u[1], u[2], u[3] (no CY relation) *)
usCY = usubOf[{1 + x, -x, -1}]; Vcy = Function[{h, a}, hodgeNormU4[h, a] /. usCY];
usN = usubOf[{3, -5, 11}]; Vnum = Function[{h, a}, hodgeNormU4[h, a] /. usN];  (* numeric non-CY weights *)
Vwk = Function[{h, a}, wk[h, a]];                                           (* trivial CohFT *)

Print["=== 1. Lemma 2 term by term: I_m(1) - (n/cc) J_m + m I_{m-1}(z phi) == 0, generic weights, symbolic cc ==="];
Do[Print["(h,n) = ", {h, n}, ": ", Table[Together[I1[h, n, m, Vgen, cc] - (n/cc) J[h, n, m, Vgen, cc] + m Izp[h, n, m - 1, Vgen, cc]] === 0, {m, 0, 3 h - 2 + n}]],
  {h, 1, 2}, {n, If[h == 1, 1, 0], 2}]; lap[];
Print["    genus 3 (numeric non-CY weights 3, -5, 11; symbolic cc): "];
Do[Print["(h,n) = ", {3, n}, ": ", Table[Together[I1[3, n, m, Vnum, cc] - (n/cc) J[3, n, m, Vnum, cc] + m Izp[3, n, m - 1, Vnum, cc]] === 0, {m, 0, 7 + n}]], {n, 0, 1}]; lap[];
Print["    the genus-1, n = 0 exception: I_0(1) = Int_{M_{1,1}} Omega = -(u1+u2+u3)/24: ", Together[I1[1, 0, 0, Vgen, cc] + (u[1] + u[2] + u[3])/24] === 0,
  ";  Lemma 2 for m >= 1 at (h,n) = (1,0): ", Table[Together[I1[1, 0, m, Vgen, cc] + m Izp[1, 0, m - 1, Vgen, cc]] === 0, {m, 1, 3}]];

Print["=== 2. Lemma 1: cc Z[h,n+1] == Sum_m (I_m(1) + I_m(z phi))/m!, generic weights, symbolic cc ==="];
Do[Print["(h,n) = ", {h, n}, ": ", Together[cc Z[h, n + 1, Vgen, cc] - Sum[(I1[h, n, m, Vgen, cc] + Izp[h, n, m, Vgen, cc])/m!, {m, 0, 3 h - 2 + n}]] === 0], {h, 1, 2}, {n, 0, 1}]; lap[];

Print["=== 3. Recursion cc^2 Z[h,n+1] == n Z[h,n] and the genus-1 anomaly ==="];
Print["generic weights, symbolic cc, h = 1, n = 1..3: ", Table[Together[cc^2 Z[1, n + 1, Vgen, cc] - n Z[1, n, Vgen, cc]] === 0, {n, 1, 3}]];
Print["generic weights, symbolic cc, h = 2, n = 0..2: ", Table[Together[cc^2 Z[2, n + 1, Vgen, cc] - n Z[2, n, Vgen, cc]] === 0, {n, 0, 2}]];
Print["anomaly: cc Z[1,1] + (u1+u2+u3)/24 = ", Together[cc Z[1, 1, Vgen, cc] + (u[1] + u[2] + u[3])/24], ";  Z[1,1] = ", Factor[Z[1, 1, Vgen, cc]]];
Print["genus 2 generic: Z[2,1] = ", Z[2, 1, Vgen, cc], ", Z[2,2] = ", Z[2, 2, Vgen, cc], ", Z[2,3] = ", Z[2, 3, Vgen, cc], ";  Z[2,0] = ", Factor[Z[2, 0, Vgen, cc]]]; lap[];
Print["genus 3, CY weights, cc = 1+x: Z[3,1] = ", Z[3, 1, Vcy, 1 + x], ", Z[3,2] = ", Z[3, 2, Vcy, 1 + x]];
Print["genus 3, numeric non-CY weights, symbolic cc: Z[3,1] = ", Z[3, 1, Vnum, cc], ", Z[3,2] = ", Z[3, 2, Vnum, cc], ";  Z[3,0] = ", Factor[Z[3, 0, Vnum, cc]]]; lap[];

Print["=== 4. Trivial CohFT (Witten-Kontsevich vertex), symbolic cc ==="];
Print["Z[1,1], Z[1,2] = ", {Z[1, 1, Vwk, cc], Z[1, 2, Vwk, cc]}, ";  Z[2,1], Z[2,2], Z[3,1], Z[3,2] = ", {Z[2, 1, Vwk, cc], Z[2, 2, Vwk, cc], Z[3, 1, Vwk, cc], Z[3, 2, Vwk, cc]}, ";  Z[2,0], Z[3,0] = ", {Z[2, 0, Vwk, cc], Z[3, 0, Vwk, cc]}]; lap[];

Print["=== 5. Quoted values (CY weights {1+x, -x, -1}, cc = 1+x) ==="];
Print["Z[1,1] - (1+x+x^2)/(24x(1+x)^2) = ", Together[Z[1, 1, Vcy, 1 + x] - (1 + x + x^2)/(24 x (1 + x)^2)], ";  Z[1,2] - (1+x+x^2)/(24x(1+x)^4) = ", Together[Z[1, 2, Vcy, 1 + x] - (1 + x + x^2)/(24 x (1 + x)^4)],
  ";  Z[1,3] - (1+x+x^2)/(12x(1+x)^6) = ", Together[Z[1, 3, Vcy, 1 + x] - (1 + x + x^2)/(12 x (1 + x)^6)]];
Print["Z[2,0] + (7+10x^2+7x^4)/(5760x^2(1+x)^3) = ", Together[Z[2, 0, Vcy, 1 + x] + (7 + 10 x^2 + 7 x^4)/(5760 x^2 (1 + x)^3)], ";  Z[2,1], Z[2,2] = ", {Z[2, 1, Vcy, 1 + x], Z[2, 2, Vcy, 1 + x]}];
Print["gap check (-x c^3)^(h-1) Z[h,0] vs gapHVC, h = 2, 3: ", Table[Together[(-x (1 + x)^3)^(h - 1) Z[h, 0, Vcy, 1 + x] - gapHVC[h, x][[1]]] === 0, {h, 2, 3}]];

Print["=== 6. The quoted h = 2 table (n = 0 legs; one special leg; cc = 1+x) ==="];
Do[Print["m = ", m, ": I_m(1) = ", Factor[I1[2, 0, m, Vcy, 1 + x]], ";  I_m(z phi)/m! = ", Factor[Izp[2, 0, m, Vcy, 1 + x]/m!]], {m, 0, 5}];
Print["A = ", Factor[Sum[I1[2, 0, m, Vcy, 1 + x]/m!, {m, 0, 5}]], ";  B = ", Factor[Sum[Izp[2, 0, m, Vcy, 1 + x]/m!, {m, 0, 4}]]];
Print["done in ", Round[AbsoluteTime[] - tt], " s"];
