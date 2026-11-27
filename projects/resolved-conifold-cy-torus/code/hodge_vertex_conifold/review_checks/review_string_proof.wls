(* review_string_proof.wls -- checks of the string-equation proof of Z[h,n] = 0 (h >= 2, n >= 1), genus <= 2 and genus 1.
   Generic weights u[1], u[2], u[3] (NOT assumed Calabi-Yau) and generic leg parameter cc; also the trivial CohFT (u = 0)
   and the project's CY weights {1+x, -x, -1} with cc = 1+x as a convention check. Read-only on the project files. *)
SetDirectory[DirectoryName[$InputFileName]];
tt = AbsoluteTime[];
Get["hvc_lib.m"];
Print["library loaded in ", Round[AbsoluteTime[] - tt], " s"];
(* multisets as sorted lists of length k, entries >= lo, sum <= smax; ordFactor = number of distinct orderings *)
msets[k_, lo_, smax_] := Module[{rec}, rec[kk_, l_, s_] := If[kk == 0, {{}}, Join @@ Table[Prepend[#, a] & /@ rec[kk - 1, a, s - a], {a, l, s}]]; rec[k, lo, smax]];
ordFactor[lst_] := Length[lst]!/Times @@ (Last /@ Tally[lst]!);
phiC[k_, cc_] := 1/cc^(k + 1);                    (* psi^k coefficient of phi = 1/(cc - psi) *)
TC[d_, cc_] := -1/cc^(d - 1);                     (* psi^d coefficient (d >= 2) of T = -psi^2/(cc - psi) *)
(* Zp[h, n, f1, usub, cc]: Sum_m 1/m! Int Omega f1(psi_1) prod_{i=2}^n phi(psi_i) prod_{j<=m} T(psi_j); f1[k] = psi^k coefficient of the first leg *)
Zp[h_, n_, f1_, usub_, cc_] := Module[{tot = 0, dim},
   Do[dim = 3 h - 3 + n + m; If[dim < 2 m, Continue[]];
    Do[If[f1[k] === 0, Continue[]];
     Do[Module[{sa = k + Total[A]},
       Do[tot += f1[k] ordFactor[A] (Times @@ (phiC[#, cc] & /@ A)) (ordFactor[D]/m!) (Times @@ (TC[#, cc] & /@ D)) (hodgeNormU4[h, Join[{k}, A, D]] /. usub),
        {D, msets[m, 2, dim - sa]}]], {A, msets[n - 1, 0, dim - k]}], {k, 0, dim}], {m, 0, 3 h - 3 + n}];
   Together[tot]];
Zfull[h_, n_, usub_, cc_] := If[n == 0, Zp[h, 1, Function[k, If[k == 0, 1, 0]], usub, cc] /. {} , Zp[h, n, Function[k, phiC[k, cc]], usub, cc]];
(* n = 0 needs its own enumeration (no legs) *)
Z0[h_, usub_, cc_] := Module[{tot = 0, dim},
   Do[dim = 3 h - 3 + m; If[dim < 2 m, Continue[]];
    Do[tot += (ordFactor[D]/m!) (Times @@ (TC[#, cc] & /@ D)) (hodgeNormU4[h, D] /. usub), {D, msets[m, 2, dim]}], {m, 0, 3 h - 3}];
   Together[tot]];
Zany[h_, n_, usub_, cc_] := If[n == 0, Z0[h, usub, cc], Zfull[h, n, usub, cc]];
one[k_] := If[k == 0, 1, 0];                       (* tau_0 at the first point *)
psiphi[cc_][k_] := If[k >= 1, 1/cc^k, 0];          (* psi phi(psi) at the first point *)
kappa = -(u[1] + u[2] + u[3])/24;                  (* Int_{M_{1,1}} Omega_{1,1} = -(1/24) Sum_a 1/w_a *)
cy = usubOf[{1 + x, -x, -1}]; triv = {u[_] -> 0}; gen = {};

Print["=== A. convention check against the task's stated values (CY weights {1+x,-x,-1}, cc = 1+x) ==="];
Print["Z[2,0] + (7+10x^2+7x^4)/(5760 x^2 (1+x)^3) = ", Together[Zany[2, 0, cy, 1 + x] + (7 + 10 x^2 + 7 x^4)/(5760 x^2 (1 + x)^3)]];
Print["Z[1,1] - (1+x+x^2)/(24x(1+x)^2) = ", Together[Zany[1, 1, cy, 1 + x] - (1 + x + x^2)/(24 x (1 + x)^2)], ",  Z[1,2] - (1+x+x^2)/(24x(1+x)^4) = ", Together[Zany[1, 2, cy, 1 + x] - (1 + x + x^2)/(24 x (1 + x)^4)]];
Print["Z[2,1], Z[2,2] (CY) = ", {Zany[2, 1, cy, 1 + x], Zany[2, 2, cy, 1 + x]}];

Print["=== B. new prediction 1: the Hodge class is irrelevant -- trivial CohFT (u = 0), symbolic cc ==="];
Print["trivial CohFT: Z[2,1], Z[2,2], Z[2,3] = ", {Zany[2, 1, triv, cc], Zany[2, 2, triv, cc], Zany[2, 3, triv, cc]}, ";  Z[2,0] = ", Zany[2, 0, triv, cc], " (nonzero expected)"];
Print["trivial CohFT, genus 1: Z[1,n], n = 1..3 = ", Table[Zany[1, n, triv, cc], {n, 1, 3}], " (zero expected: kappa = 0 for the point)"];

Print["=== C. new prediction 2: arbitrary (non-CY) weights u[1],u[2],u[3] and arbitrary cc ==="];
Print["generic: Z[2,1] = ", Zany[2, 1, gen, cc]];
Print["generic: Z[2,2] = ", Zany[2, 2, gen, cc]];
Print["generic: Z[2,0] = ", Factor[Zany[2, 0, gen, cc]]];

Print["=== D. new prediction 3: genus 1, Z[1,n] = (n-1)! kappa / cc^(2n-1) with kappa = -(u1+u2+u3)/24, generic u and cc ==="];
Do[Print["n = ", n, ": Z[1,n] - (n-1)! kappa/cc^(2n-1) = ", Together[Zany[1, n, gen, cc] - (n - 1)! kappa/cc^(2 n - 1)], "   (Z[1,n] = ", Factor[Zany[1, n, gen, cc]], ")"], {n, 1, 4}];

Print["=== E. the two string-equation pieces separately: Apiece = Sum 1/m! <1, phi^(n-1), T^m>, Bpiece = Sum 1/m! <psi phi, phi^(n-1), T^m> ==="];
Print["   identity used: cc Z[h,n] = Apiece + Bpiece (cc phi = 1 + psi phi) and string equation: Apiece = (n-1)/cc Z[h,n-1] - Bpiece + delta_{h1} delta_{n1} kappa"];
Do[Module[{Ap = Zp[h, n, one, gen, cc], Bp = Zp[h, n, psiphi[cc], gen, cc], Zprev = If[n == 1, Zany[h, 0, gen, cc], Zany[h, n - 1, gen, cc]]},
   Print["(h,n) = ", {h, n}, ": Apiece = ", Factor[Ap], ";  Bpiece = ", Factor[Bp], ";  Apiece + Bpiece - (n-1)/cc Z[h,n-1] - delta kappa = ",
    Together[Ap + Bp - (n - 1)/cc Zprev - If[h == 1 && n == 1, kappa, 0]]]], {h, 1, 2}, {n, 1, 2}];
Print["total time ", Round[AbsoluteTime[] - tt], " s"];
