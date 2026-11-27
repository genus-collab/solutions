(* review_string_proof_g3.wls -- genus-3 checks of the string-equation proof with generic weights and generic cc (time-guarded). *)
SetDirectory[DirectoryName[$InputFileName]];
tt = AbsoluteTime[];
Get["hvc_lib.m"];
msets[k_, lo_, smax_] := Module[{rec}, rec[kk_, l_, s_] := If[kk == 0, {{}}, Join @@ Table[Prepend[#, a] & /@ rec[kk - 1, a, s - a], {a, l, s}]]; rec[k, lo, smax]];
ordFactor[lst_] := Length[lst]!/Times @@ (Last /@ Tally[lst]!);
phiC[k_, cc_] := 1/cc^(k + 1); TC[d_, cc_] := -1/cc^(d - 1);
Zp[h_, n_, f1_, usub_, cc_] := Module[{tot = 0, dim},
   Do[dim = 3 h - 3 + n + m; If[dim < 2 m, Continue[]];
    Do[If[f1[k] === 0, Continue[]];
     Do[Module[{sa = k + Total[A]},
       Do[tot += f1[k] ordFactor[A] (Times @@ (phiC[#, cc] & /@ A)) (ordFactor[D]/m!) (Times @@ (TC[#, cc] & /@ D)) (hodgeNormU4[h, Join[{k}, A, D]] /. usub),
        {D, msets[m, 2, dim - sa]}]], {A, msets[n - 1, 0, dim - k]}], {k, 0, dim}], {m, 0, 3 h - 3 + n}];
   Together[tot]];
one[k_] := If[k == 0, 1, 0]; psiphi[cc_][k_] := If[k >= 1, 1/cc^k, 0];
triv = {u[_] -> 0}; gen = {};
budget = 150;
guard[label_, expr_] := If[AbsoluteTime[] - tt > budget, Print[label, ": skipped (time budget)"],
   Print[label, " = ", TimeConstrained[expr, budget - (AbsoluteTime[] - tt), "timed out"], "   [t = ", Round[AbsoluteTime[] - tt], " s]"]];
guard["trivial CohFT, Z[3,1]", Zp[3, 1, Function[k, phiC[k, cc]], triv, cc]];
guard["generic weights, Z[3,1]", Zp[3, 1, Function[k, phiC[k, cc]], gen, cc]];
guard["generic weights, Apiece(3,1) + Bpiece(3,1)", Together[Zp[3, 1, one, gen, cc] + Zp[3, 1, psiphi[cc], gen, cc]]];
guard["generic weights, Apiece(3,1) alone (nonzero expected)", Factor[Zp[3, 1, one, gen, cc]]];
guard["trivial CohFT, Z[3,2]", Zp[3, 2, Function[k, phiC[k, cc]], triv, cc]];
guard["generic weights, Z[3,2]", Zp[3, 2, Function[k, phiC[k, cc]], gen, cc]];
Print["total time ", Round[AbsoluteTime[] - tt], " s"];
