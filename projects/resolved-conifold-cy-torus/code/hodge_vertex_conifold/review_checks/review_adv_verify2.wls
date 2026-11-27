(* review_adv_verify2.wls -- genus-3 adversarial checks (symbolic weights u[a], symbolic c, s). Read-only on the project files. *)
SetDirectory[DirectoryName[$InputFileName]];
tt = AbsoluteTime[];
Get["hvc_lib.m"];
corrU[h_, a_List] := hodgeNormU4[h, a];
tuples[g_, nmax_] := Join @@ Table[Tuples[Range[0, 3 g - 3 + n], n], {n, 1, nmax}];
strOK[g_, a_] := Together[corrU[g, Append[a, 0]] - Sum[If[a[[i]] >= 1, corrU[g, ReplacePart[a, i -> a[[i]] - 1]], 0], {i, Length[a]}]] === 0;
Print["(S) genus 3, N <= 2, symbolic weights: ", And @@ (strOK[3, #] & /@ tuples[3, 2]), "   (", Round[AbsoluteTime[] - tt], " s)"];
otuples[k_, lo_, hi_, smax_] := If[k == 0, {{}}, Join @@ Table[Prepend[#, a] & /@ otuples[k - 1, lo, hi, smax - a], {a, lo, Min[hi, smax]}]];
phiCoef[k_, cc_] := 1/cc^(k + 1);
TCoef[d_, cc_] := If[d >= 2, -1/cc^(d - 1), 0];
Zord[h_, n_, corr_, cc_] := Module[{tot = 0, dim},
   Do[dim = 3 h - 3 + n + m;
    If[2 h - 2 + n + m <= 0 || dim < 2 m, Continue[]];
    Do[Do[tot += (Times @@ (phiCoef[#, cc] & /@ A)) (Times @@ (TCoef[#, cc] & /@ D)) corr[h, Join[A, D]]/m!,
       {D, otuples[m, 2, dim - Total[A], dim - Total[A]]}], {A, otuples[n, 0, dim, dim]}], {m, 0, 3 h - 3 + n}];
   Together[tot]];
tsCoef[k_] := s c^(-k) - If[k >= 2, c^(1 - k), 0];
gCoef[k_] := If[k == 0, 1, 0] + (s - c) c^(-k - 1);
one[k_] := If[k == 0, 1, 0];
corrIns[h_, Acoef_, Bcoef_, M_] := Module[{dim = 3 h - 3 + M + 1, tot = 0},
   If[2 h - 2 + M + 1 <= 0, Return[0]];
   Do[Do[tot += Acoef[a0] (Times @@ (Bcoef /@ B)) corrU[h, Prepend[B, a0]], {B, otuples[M, 0, dim - a0, dim - a0]}], {a0, 0, dim}];
   Together[tot]];
Print["genus 3: <tau_0 t_s^M>_{3,M+1} - M <(1+(s-c)phi) t_s^(M-1)>_{3,M}, M = 1..3: ", Table[Together[corrIns[3, one, tsCoef, M] - M corrIns[3, gCoef, tsCoef, M - 1]], {M, 1, 3}], "   (", Round[AbsoluteTime[] - tt], " s)"];
Print["Z[3,1] (symbolic u, c) = ", Zord[3, 1, corrU, c], "   (", Round[AbsoluteTime[] - tt], " s)"];
Print["Z[3,0] at CY weights: ", Factor[Zord[3, 0, corrU, c] /. usubOf[{1 + x, -x, -1}] /. c -> 1 + x], "   (", Round[AbsoluteTime[] - tt], " s)"];
Print["Z[3,2] (symbolic u, c) = ", Zord[3, 2, corrU, c], "   (", Round[AbsoluteTime[] - tt], " s)"];
Print["total time ", Round[AbsoluteTime[] - tt], " s"];
