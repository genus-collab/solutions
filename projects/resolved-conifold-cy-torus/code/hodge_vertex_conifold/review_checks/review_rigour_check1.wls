(* review_rigour_check1.wls -- adversarial verifier (rigour lens), independent check of the string-equation proof of Z[h,n] = 0 (h >= 2, n >= 1).
   Own enumeration (sorted exponent lists with multinomial weights; every tuple with total psi-degree <= dim contributes, the class
   correlator picks the matching degree component). Generic weights u[1], u[2], u[3] (no CY relation) and generic leg parameter cc.
   Read-only on the project files. *)
SetDirectory[DirectoryName[$InputFileName]];
t0 = AbsoluteTime[]; el[] := Round[AbsoluteTime[] - t0, 0.1];
Get["hvc_lib.m"]; Print["lib loaded ", el[], " s"];

(* leaf series as coefficient functions k -> coefficient of psi^k *)
phiK[cc_][k_] := cc^(-k - 1);                      (* phi = 1/(cc - psi) *)
TK[cc_][k_] := If[k >= 2, -cc^(1 - k), 0];         (* T = -psi^2/(cc - psi) *)
oneK[k_] := If[k == 0, 1, 0];                      (* 1 *)
zphiK[cc_][k_] := If[k >= 1, cc^(-k), 0];          (* psi phi *)
(* sorted exponent lists of length k, entries >= lo, sum <= smax *)
sl[k_, lo_, smax_] := Module[{rec}, rec[kk_, l_, s_] := If[kk == 0, {{}}, Join @@ Table[Prepend[#, a] & /@ rec[kk - 1, a, s - a], {a, l, s}]]; rec[k, lo, smax]];
nord[l_] := Multinomial @@ (Last /@ Tally[l]);
(* Wm[h, f, n, cc, cls] = list over m of (1/m!) < f, phi^n, T^m >_{h, 1+n+m}; f = None: no distinguished leg, (1/m!) < phi^n, T^m >_{h, n+m}.
   cls[h, a] = Int_{M_{h,|a|}} Omega_h psi^a (mixed-degree class: all tuples with Total[a] <= dim contribute). *)
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
hodgeGen[h_, a_] := hodgeNormU4[h, a];                       (* generic weights u[1], u[2], u[3] *)
hodgeCY[h_, a_] := hodgeNormU4[h, a] /. usubOf[{1 + x, -x, -1}];   (* the project's CY weights, eps2 = 1, eps1 = x *)
triv[h_, a_] := hodgeNormU4[h, a] /. u[_] -> 0;              (* trivial CohFT (Witten-Kontsevich) *)
kappaGen = hodgeNormU4[1, {0}];
Print["kappa = <1>_{1,1} for generic weights: ", kappaGen, "   (expected -(u1+u2+u3)/24)"];

budget = 165;
guard[label_, expr_] := If[AbsoluteTime[] - t0 > budget, Print[label, ": SKIPPED (time budget)"],
   Print[label, " = ", TimeConstrained[expr, budget - (AbsoluteTime[] - t0), "TIMED OUT"], "   [", el[], " s]"]];

Print["=== 1. trivial CohFT, genus 2 (independent of the Hodge machinery: pure WK numbers) ==="];
guard["Z^WK[2,1] by dilaton order m (hand values 5,-78,273,-340,140 over 5760 cc^5)", 5760 cc^5 Wm[2, None, 1, cc, triv]];
guard["Z^WK[2,1], Z^WK[2,2]", {Zc[2, 1, cc, triv], Zc[2, 2, cc, triv]}];
guard["Z^WK[2,0] (expected -1/(1440 cc^3))", Zc[2, 0, cc, triv]];

Print["=== 2. genus 1, generic weights and cc: Z[1,n] vs (n-1)! kappa cc^(1-2n) ==="];
Do[guard["n = " <> ToString[n] <> ": Z[1,n] - (n-1)! kappa cc^(1-2n)", Together[Zc[1, n, cc, hodgeGen] - (n - 1)! kappaGen cc^(1 - 2 n)]], {n, 1, 3}];

Print["=== 3. genus 2, generic weights and cc ==="];
guard["Z[2,1]", Zc[2, 1, cc, hodgeGen]];
guard["Z[2,2]", Zc[2, 2, cc, hodgeGen]];
guard["Z[2,0] (nonzero expected)", Factor[Zc[2, 0, cc, hodgeGen]]];

Print["=== 4. Lemma 3 piece by piece, generic weights and cc: A = Sum 1/m! <1, phi^(n-1), T^m>, B = Sum 1/m! <psi phi, phi^(n-1), T^m> ==="];
Do[Module[{Ap, Bp, Zn, Zp},
   If[AbsoluteTime[] - t0 > budget, Print["(h,n) = ", {h, n}, ": SKIPPED (time budget)"],
    Ap = W[h, oneK, n - 1, cc, hodgeGen]; Bp = W[h, zphiK[cc], n - 1, cc, hodgeGen]; Zn = Zc[h, n, cc, hodgeGen]; Zp = If[n == 1, 0, Zc[h, n - 1, cc, hodgeGen]];
    Print["(h,n) = ", {h, n}, ": A = ", Factor[Ap], ";  B = ", Factor[Bp],
     ";  cc Z[h,n] - (A + B) = ", Together[cc Zn - Ap - Bp],
     ";  A - (n-1)/cc Z[h,n-1] + B - delta kappa = ", Together[Ap - (n - 1)/cc Zp + Bp - If[h == 1 && n == 1, kappaGen, 0]], "   [", el[], " s]"]]],
  {h, 1, 2}, {n, 1, 2}];

Print["=== 5. convention check with the project's CY weights (1+x, -x, -1), cc = 1+x, against the task's stated values ==="];
guard["Z[1,1] - (1+x+x^2)/(24x(1+x)^2), Z[1,2] - (1+x+x^2)/(24x(1+x)^4), Z[1,3] - (1+x+x^2)/(12x(1+x)^6)",
  Together /@ {Zc[1, 1, 1 + x, hodgeCY] - (1 + x + x^2)/(24 x (1 + x)^2), Zc[1, 2, 1 + x, hodgeCY] - (1 + x + x^2)/(24 x (1 + x)^4), Zc[1, 3, 1 + x, hodgeCY] - (1 + x + x^2)/(12 x (1 + x)^6)}];
guard["Z[2,0] + (7+10x^2+7x^4)/(5760 x^2 (1+x)^3)", Together[Zc[2, 0, 1 + x, hodgeCY] + (7 + 10 x^2 + 7 x^4)/(5760 x^2 (1 + x)^3)]];
Print["total ", el[], " s"];
