(* review_rigour_check3.wls -- adversarial verifier (rigour lens), genus 4 with generic weights u[1..3] and generic cc (the case the prover did not run),
   plus a direct spot check of Lemma 1 (string equation) for the normalised Hodge correlators at genus 3-4 with 4-5 marked points.
   Same independent enumeration as review_rigour_check1.wls. Read-only on the project files. *)
SetDirectory[DirectoryName[$InputFileName]];
t0 = AbsoluteTime[]; el[] := Round[AbsoluteTime[] - t0, 0.1];
Get["hvc_lib.m"]; Print["lib loaded ", el[], " s"];
phiK[cc_][k_] := cc^(-k - 1); TK[cc_][k_] := If[k >= 2, -cc^(1 - k), 0]; oneK[k_] := If[k == 0, 1, 0]; zphiK[cc_][k_] := If[k >= 1, cc^(-k), 0];
sl[k_, lo_, smax_] := Module[{rec}, rec[kk_, l_, s_] := If[kk == 0, {{}}, Join @@ Table[Prepend[#, a] & /@ rec[kk - 1, a, s - a], {a, l, s}]]; rec[k, lo, smax]];
nord[l_] := Multinomial @@ (Last /@ Tally[l]);
Wm[h_, f_, n_, cc_, cls_] := Module[{N0 = If[f === None, 0, 1], res = {}, dim, tot, ks, k0},
  Do[dim = 3 h - 3 + N0 + n + m; tot = 0;
   If[2 h - 2 + N0 + n + m > 0 && dim >= 2 m,
    ks = If[f === None, {Null}, Range[0, dim]];
    Do[k0 = If[k === Null, 0, k];
     Do[Do[tot += If[k === Null, 1, f[k]] nord[A] nord[D] (Times @@ (phiK[cc] /@ A)) (Times @@ (TK[cc] /@ D))/m! cls[h, Join[If[k === Null, {}, {k}], A, D]],
       {D, sl[m, 2, dim - k0 - Total[A]]}], {A, sl[n, 0, dim - k0]}], {k, ks}]];
   AppendTo[res, Together[tot]], {m, 0, 3 h - 3 + N0 + n}];
  res];
W[h_, f_, n_, cc_, cls_] := Together[Total[Wm[h, f, n, cc, cls]]];
Zc[h_, n_, cc_, cls_] := W[h, None, n, cc, cls];
hodgeGen[h_, a_] := hodgeNormU4[h, a];
budget = 160;
guard[label_, expr_] := If[AbsoluteTime[] - t0 > budget, Print[label, ": SKIPPED (time budget)"],
   Print[label, " = ", TimeConstrained[expr, budget - (AbsoluteTime[] - t0), "TIMED OUT"], "   [", el[], " s]"]];

Print["=== 1. Lemma 1 (string equation) for the normalised Hodge correlators, generic weights, random tuples with 4-5 points, genus 3 and 4 ==="];
SeedRandom[20261030];
strOK[h_, a_List] := Together[hodgeNormU4[h, Append[a, 0]] - Sum[If[a[[i]] >= 1, hodgeNormU4[h, ReplacePart[a, i -> a[[i]] - 1]], 0], {i, Length[a]}]] === 0;
randTuples[h_, n_, k_] := Table[RandomInteger[{0, 3 h - 3 + n}, n], {k}];
guard["genus 3, 4 points, 6 random tuples", And @@ (strOK[3, #] & /@ randTuples[3, 4, 6])];
guard["genus 4, 4 points, 4 random tuples", And @@ (strOK[4, #] & /@ randTuples[4, 4, 4])];
guard["genus 4, 5 points, 3 random tuples", And @@ (strOK[4, #] & /@ randTuples[4, 5, 3])];
guard["genus 4, <1>_{4,1} = 0 and <1, psi^k>_{4,2} = <psi^(k-1)>_{4,1}, k = 1..10", And @@ Prepend[Table[strOK[4, {k}], {k, 1, 10}], hodgeNormU4[4, {0}] === 0]];

Print["=== 2. genus 4, generic weights u[1..3] and generic cc (not run by the prover) ==="];
guard["Z[4,1]", Zc[4, 1, cc, hodgeGen]];
Module[{Ap, Bp},
  If[AbsoluteTime[] - t0 > budget, Print["pieces (4,1): SKIPPED (time budget)"],
   Ap = W[4, oneK, 0, cc, hodgeGen]; Bp = W[4, zphiK[cc], 0, cc, hodgeGen];
   Print["(h,n) = (4,1): A nonzero? ", Ap =!= 0, ";  B nonzero? ", Bp =!= 0, ";  A + B = ", Together[Ap + Bp], ";  A at u = (1/2,1/3,1/5), cc = 7: ", Ap /. {u[1] -> 1/2, u[2] -> 1/3, u[3] -> 1/5, cc -> 7}, "   [", el[], " s]"]]];
guard["Z[4,2]", Zc[4, 2, cc, hodgeGen]];
guard["Z[4,0] at u = (1/2, 1/3, 1/5), cc = 7 (nonzero expected)", Zc[4, 0, 7, Function[{h, a}, hodgeNormU4[h, a] /. {u[1] -> 1/2, u[2] -> 1/3, u[3] -> 1/5}]]];
Print["total ", el[], " s"];
