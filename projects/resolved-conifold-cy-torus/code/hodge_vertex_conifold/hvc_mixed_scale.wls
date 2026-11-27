(* hvc_mixed_scale.wls -- mixed-scale vertex functionals: legs phi_c = 1/(c - psi) (scale c) plus ONE leg phi_{cp} = 1/(cp - psi) of another
   scale cp (symbolic), dilaton leaves T_c = -psi^2/(c - psi), Hodge class with weights (c, -x, -1), c = 1 + x.
   Y[h, n] := Sum_m 1/m! < phi_cp phi_c^n T_c^m >_h.  Theorem C: Y with cp = c is Z[h, n+1] = 0.  Question: Y[h,n] = 0 for all cp? *)
SetDirectory[DirectoryName[$InputFileName]];
Get["hvc_lib.m"];
c = 1 + x; usub = usubOf[{c, -x, -1}];
phiC[k_] := 1/c^(k + 1); phiP[k_] := 1/cp^(k + 1); TC[d_] := -1/c^(d - 1);
msets[k_, lo_, hi_, smax_] := Module[{rec}, rec[kk_, l_, s_] := If[kk == 0, {{}}, Join @@ Table[Prepend[#, a] & /@ rec[kk - 1, a, s - a], {a, l, Min[hi, s]}]]; rec[k, lo, smax]];
ordFactor[lst_] := Length[lst]!/Times @@ (Last /@ Tally[lst]!);
Y[h_, n_] := Y[h, n] = Module[{tot = 0, dim},
   Do[dim = 3 h - 3 + n + 1 + m; If[dim < 2 m, Continue[]];
    Do[Do[Do[Module[{lp = Join[{ap}, A, D]}, If[Total[lp] <= dim, tot += phiP[ap] ordFactor[A] (ordFactor[D]/m!) (Times @@ (phiC /@ A)) (Times @@ (TC /@ D)) (hodgeNormU4[h, lp] /. usub)]],
       {D, msets[m, 2, dim, dim - ap - Total[A]]}], {A, msets[n, 0, dim, dim - ap]}], {ap, 0, dim}], {m, 0, 3 h - 3 + n + 1}];
   Factor[Together[tot]]];
Print["genus 2: Y[2,0] (one leg of scale cp, no scale-c legs) = ", Y[2, 0]];
Print["genus 2: Y[2,1] = ", Y[2, 1], ";  Y[2,2] = ", Y[2, 2]];
Print["check cp -> c: ", Factor[Y[2, 0] /. cp -> c], " ", Factor[Y[2, 1] /. cp -> c]];
Print["genus 3: Y[3,0] = ", Y[3, 0]];
Print["genus 3: Y[3,1] = ", Y[3, 1]];
Print["genus 1 for comparison: Y[1,0] = ", Y[1, 0], ";  Y[1,1] = ", Y[1, 1]];
