(* review_rigour_check2.wls -- adversarial verifier (rigour lens), genus 3 with generic weights u[1..3] and generic cc; trivial CohFT at genus 3 and 4.
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
trivWK[h_, a_] := HI[h, 0, 0, a];                            (* trivial CohFT directly from the WK numbers (DVV recursion of graphs.m) *)
budget = 165;
guard[label_, expr_] := If[AbsoluteTime[] - t0 > budget, Print[label, ": SKIPPED (time budget)"],
   Print[label, " = ", TimeConstrained[expr, budget - (AbsoluteTime[] - t0), "TIMED OUT"], "   [", el[], " s]"]];

Print["=== 1. trivial CohFT (WK numbers), genus 3 and 4 ==="];
guard["Z^WK[3,1], Z^WK[3,2], Z^WK[3,3]", Table[Zc[3, n, cc, trivWK], {n, 1, 3}]];
guard["Z^WK[3,0] (nonzero expected)", Zc[3, 0, cc, trivWK]];
guard["Z^WK[4,1], Z^WK[4,2]", Table[Zc[4, n, cc, trivWK], {n, 1, 2}]];
guard["Z^WK[4,0] (nonzero expected)", Zc[4, 0, cc, trivWK]];

Print["=== 2. genus 3, generic weights u[1..3] and generic cc ==="];
guard["Z[3,1]", Zc[3, 1, cc, hodgeGen]];
guard["Z[3,2]", Zc[3, 2, cc, hodgeGen]];
guard["Z[3,0] at u = (1/2, 1/3, 1/5), cc = 7 (nonzero expected)", Zc[3, 0, 7, Function[{h, a}, hodgeNormU4[h, a] /. {u[1] -> 1/2, u[2] -> 1/3, u[3] -> 1/5}]]];
Print["=== 3. Lemma 3 pieces at genus 3, generic weights and cc ==="];
Do[Module[{Ap, Bp, Zn, Zp},
   If[AbsoluteTime[] - t0 > budget, Print["(h,n) = ", {3, n}, ": SKIPPED (time budget)"],
    Ap = W[3, oneK, n - 1, cc, hodgeGen]; Bp = W[3, zphiK[cc], n - 1, cc, hodgeGen]; Zn = Zc[3, n, cc, hodgeGen]; Zp = If[n == 1, 0, Zc[3, n - 1, cc, hodgeGen]];
    Print["(h,n) = ", {3, n}, ": A nonzero? ", Ap =!= 0, ";  B nonzero? ", Bp =!= 0, ";  cc Z - (A + B) = ", Together[cc Zn - Ap - Bp],
     ";  A - (n-1)/cc Z[h,n-1] + B = ", Together[Ap - (n - 1)/cc Zp + Bp], "   [", el[], " s]"]]], {n, 1, 2}];
guard["Z[3,3] generic", Zc[3, 3, cc, hodgeGen]];
Print["total ", el[], " s"];
