(* hvc_which_class.wls -- which Hodge classes give the mixed-scale zeros?  Single probe lam_j = 1/(cl - j psi), dilaton leaves -psi^2/(cl - psi),
   class prod_a Sum_i (-1)^i lambda_i u_a^i with u_a = 1/w_a (u_a = 0 drops the factor).  Prints the factorisation of Xi_h(j) in j. *)
SetDirectory[DirectoryName[$InputFileName]];
Get["hvc_lib.m"];
msets[k_, lo_, hi_, smax_] := Module[{rec}, rec[kk_, l_, s_] := If[kk == 0, {{}}, Join @@ Table[Prepend[#, a] & /@ rec[kk - 1, a, s - a], {a, l, Min[hi, s]}]]; rec[k, lo, smax]];
ordFactor[lst_] := Length[lst]!/Times @@ (Last /@ Tally[lst]!);
XiOne[h_, cl_, us_] := Module[{tot = 0, TC, dim},
   TC[d_] := -1/cl^(d - 1);
   Do[Do[dim = 3 h - 2 + m;
      Do[tot += jv^a/cl^(a + 1) (ordFactor[D]/m!) (Times @@ (TC /@ D)) (hodgeNormU4[h, Join[{a}, D]] /. us), {D, msets[m, 2, dim, dim - a]}],
      {m, 0, 3 h - 2 - a}], {a, 0, 3 h - 2}];
   Factor[Together[tot]]];
cases = {
   {"trivial class (KW)", 1, {u[1] -> 0, u[2] -> 0, u[3] -> 0}},
   {"single Lambda(cl), cl = 1", 1, {u[1] -> 1, u[2] -> 0, u[3] -> 0}},
   {"Lambda(cl) Lambda(w2), cl = 1, w2 generic", 1, {u[1] -> 1, u[2] -> 1/w2, u[3] -> 0}},
   {"Lambda(cl) Lambda(-cl) (Mumford, trivial up to sign)", 1, {u[1] -> 1, u[2] -> -1, u[3] -> 0}},
   {"triple, cl = w1 = 1, w2, w3 generic (non-CY)", 1, {u[1] -> 1, u[2] -> 1/w2, u[3] -> 1/w3}},
   {"CY triple (1, w2, -1-w2), leaf scale cl = 1 = w1 (control)", 1, {u[1] -> 1, u[2] -> 1/w2, u[3] -> 1/(-1 - w2)}},
   {"CY triple (1, w2, -1-w2), leaf scale cl = 2 (not a weight)", 2, {u[1] -> 1, u[2] -> 1/w2, u[3] -> 1/(-1 - w2)}},
   {"triple (w1, w2, w3) generic, leaf scale cl generic", cl, {u[1] -> 1/w1, u[2] -> 1/w2, u[3] -> 1/w3}}};
Do[Print["genus ", h, ", ", cs[[1]], ":  Xi(j) = ", XiOne[h, cs[[2]], cs[[3]]]], {h, 2, 3}, {cs, If[h == 2, cases, cases[[;; 7]]]}];
