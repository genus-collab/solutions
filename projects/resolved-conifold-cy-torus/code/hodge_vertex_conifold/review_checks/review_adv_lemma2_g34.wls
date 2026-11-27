(* review_adv_lemma2_g34.wls -- adversarial verification, part 2: Lemma 2 (string step) term by term in genus 3 (n = 0, 1, 2, all m)
   and genus 4 (n = 0, 1), plus Z[3,2], Z[3,3], Z[4,1], Z[4,2]; time-guarded. Independent enumeration code; read-only on the project files. *)
SetDirectory[DirectoryName[$InputFileName]];
t0 = AbsoluteTime[];
Get["hvc_lib.m"];
el[] := Round[AbsoluteTime[] - t0];
corrLE[h_, legs_List, usub_] := Module[{N = Length[legs], dim, tot = 0, rec},
   dim = 3 h - 3 + N;
   If[2 h - 2 + N <= 0, Return[0]];
   rec[i_, prefix_, budget_, coef_] := If[i > N, tot += coef (hodgeNormU4[h, prefix] /. usub),
      Do[Module[{cf = legs[[i]][k]}, If[cf =!= 0, rec[i + 1, Append[prefix, k], budget - k, coef cf]]], {k, 0, budget}]];
   rec[1, {}, dim, 1];
   Together[tot]];
phiF[cc_] := Function[k, 1/cc^(k + 1)];
TF[cc_] := Function[k, If[k >= 2, -1/cc^(k - 1), 0]];
oneF = Function[k, If[k == 0, 1, 0]];
zphiF[cc_] := Function[k, If[k >= 1, 1/cc^k, 0]];
Jm[h_, n_, m_, cc_, usub_] := corrLE[h, Join[ConstantArray[phiF[cc], n], ConstantArray[TF[cc], m]], usub];
Ileaf[h_, n_, m_, f_, cc_, usub_] := corrLE[h, Join[ConstantArray[phiF[cc], n], {f}, ConstantArray[TF[cc], m]], usub];
Zf[h_, n_, cc_, usub_] := Together[Sum[Jm[h, n, m, cc, usub]/m!, {m, 0, 3 h - 3 + n}]];
usCY = usubOf[{1 + x, -x, -1}];
usN = usubOf[{2, 5, -7}];
lemma2m[h_, n_, m_, cc_, usub_] := Module[{lhs, rhs},
   lhs = Ileaf[h, n, m, oneF, cc, usub];
   rhs = (n/cc) Jm[h, n, m, cc, usub] - If[m >= 1, m Ileaf[h, n, m - 1, zphiF[cc], cc, usub], 0];
   {m, Together[lhs - rhs] === 0, lhs =!= 0 || rhs =!= 0}];
Print["genus 3, Lemma 2, CY weights, c = 1+x: {m, holds, nontrivial}"];
Do[If[el[] < 100, Print["  n = ", n, ": ", Table[lemma2m[3, n, m, 1 + x, usCY], {m, 0, 3 3 - 2 + n}], "  [", el[], " s]"]], {n, 0, 2}];
If[el[] < 100, Print["  n = 1, symbolic c unrelated to the weights: ", Table[lemma2m[3, 1, m, cc, usCY], {m, 0, 7}], "  [", el[], " s]"]];
If[el[] < 110, Print["Z[3,1], Z[3,2], Z[3,3] (CY, c = 1+x) = ", {Zf[3, 1, 1 + x, usCY], Zf[3, 2, 1 + x, usCY], Zf[3, 3, 1 + x, usCY]}, ";  (-x c^3)^2 Z[3,0] - gap_3 = ", Together[(x (1 + x)^3)^2 Zf[3, 0, 1 + x, usCY] - gapHVC[3, x][[1]]], "  [", el[], " s]"]];
If[el[] < 110, Print["Z[3,1] with symbolic c: ", Zf[3, 1, cc, usCY], ";  Z[3,1] non-CY weights (2,5,-7), c = 3: ", Zf[3, 1, 3, usN], "  [", el[], " s]"]];
Print["genus 4, Lemma 2 with n = 0 (the step that gives Z[4,1] = 0), m = 0..10, CY weights, c = 1+x:"];
Do[If[el[] < 120, Print["  ", lemma2m[4, 0, m, 1 + x, usCY], "  [", el[], " s]"]], {m, 0, 10}];
If[el[] < 130, Print["Z[4,1] = ", Zf[4, 1, 1 + x, usCY], "  [", el[], " s]"]];
Do[If[el[] < 150, Print["  genus 4, n = 1: ", lemma2m[4, 1, m, 1 + x, usCY], "  [", el[], " s]"]], {m, 0, 11}];
If[el[] < 160, Print["Z[4,2] = ", Zf[4, 2, 1 + x, usCY], "  [", el[], " s]"]];
Print["total ", el[], " s"];
