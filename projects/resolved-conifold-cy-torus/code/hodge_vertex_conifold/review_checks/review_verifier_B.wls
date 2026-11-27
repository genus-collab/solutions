(* review_verifier_B.wls -- Lemma 2 term by term at genus 3 (n = 2) and genus 4 (n = 0, 1), CY weights {1+x, -x, -1}, cc = 1+x. Read-only. *)
SetDirectory[DirectoryName[$InputFileName]];
tt = AbsoluteTime[];
Get["hvc_lib.m"];
lap[] := (Print["   [t = ", Round[AbsoluteTime[] - tt], " s]"]);
ms[k_, lo_, smax_] := Module[{rec}, rec[kk_, l_, s_] := If[kk == 0, {{}}, Join @@ Table[Prepend[#, a] & /@ rec[kk - 1, a, s - a], {a, l, s}]]; rec[k, lo, smax]];
nord[lst_] := Length[lst]!/Times @@ (Last /@ Tally[lst]!);
Corr[h_, V_, types_List] := Module[{N = Total[types[[All, 2]]], dim, tot = 0, rec},
   If[2 h - 2 + N <= 0, Return[0]];
   dim = 3 h - 3 + N;
   rec[{}, acc_, budget_] := (tot += (Times @@ (Last /@ acc)) V[h, Join @@ (First /@ acc)]);
   rec[{{f_, k_, lo_}, rest___}, acc_, budget_] := Do[rec[{rest}, Append[acc, {A, nord[A] (Times @@ (f /@ A))}], budget - Total[A]], {A, ms[k, lo, budget]}];
   rec[types, {}, dim];
   Together[tot]];
phiF[cc_] := Function[k, cc^(-k - 1)]; TF[cc_] := Function[d, -cc^(1 - d)]; oneF = Function[k, If[k == 0, 1, 0]]; zphiF[cc_] := Function[k, If[k >= 1, cc^(-k), 0]];
J[h_, n_, m_, V_, cc_] := Corr[h, V, {{phiF[cc], n, 0}, {TF[cc], m, 2}}];
I1[h_, n_, m_, V_, cc_] := Corr[h, V, {{phiF[cc], n, 0}, {oneF, 1, 0}, {TF[cc], m, 2}}];
Izp[h_, n_, m_, V_, cc_] := If[m < 0, 0, Corr[h, V, {{phiF[cc], n, 0}, {zphiF[cc], 1, 1}, {TF[cc], m, 2}}]];
Z[h_, n_, V_, cc_] := Z[h, n, V, cc] = Together[Sum[J[h, n, m, V, cc]/m!, {m, 0, 3 h - 3 + n}]];
usCY = usubOf[{1 + x, -x, -1}]; Vcy = Function[{h, a}, hodgeNormU4[h, a] /. usCY];
c = 1 + x;
Print["genus 3, n = 2, Lemma 2 term by term (m = 0..9): ", Table[Together[I1[3, 2, m, Vcy, c] - (2/c) J[3, 2, m, Vcy, c] + m Izp[3, 2, m - 1, Vcy, c]] === 0, {m, 0, 9}]]; lap[];
Print["genus 3: Z[3,3] = ", Z[3, 3, Vcy, c], ";  c^2 Z[3,3] - 2 Z[3,2] = ", Together[c^2 Z[3, 3, Vcy, c] - 2 Z[3, 2, Vcy, c]]]; lap[];
Print["genus 4, n = 0, Lemma 2 term by term (m = 0..10): ", Table[Together[I1[4, 0, m, Vcy, c] + m Izp[4, 0, m - 1, Vcy, c]] === 0, {m, 0, 10}]]; lap[];
Print["genus 4: Z[4,1] = ", Z[4, 1, Vcy, c], ";  gap check (-x c^3)^3 Z[4,0] - gap_4: ", Together[(-x c^3)^3 Z[4, 0, Vcy, c] - gapHVC[4, x][[1]]]]; lap[];
Print["genus 4, n = 1, Lemma 2 term by term (m = 0..11): ", Table[Together[I1[4, 1, m, Vcy, c] - (1/c) J[4, 1, m, Vcy, c] + m Izp[4, 1, m - 1, Vcy, c]] === 0, {m, 0, 11}]]; lap[];
Print["genus 4: Z[4,2] = ", Z[4, 2, Vcy, c]];
Print["done in ", Round[AbsoluteTime[] - tt], " s"];
