(* review_adv_verify1.wls -- independent adversarial checks of the string-equation proof of Z[h,n] = 0 (h >= 2, n >= 1).
   Independent implementation: ordered-tuple enumeration, symbolic weights u[a] = 1/w_a (NOT specialised to Calabi-Yau), symbolic c and s.
   Read-only on the project files. *)
SetDirectory[DirectoryName[$InputFileName]];
tt = AbsoluteTime[];
Get["hvc_lib.m"];
Print["library loaded in ", Round[AbsoluteTime[] - tt], " s"];
corrU[h_, a_List] := hodgeNormU4[h, a];           (* <prod tau_{a_i}>_h as a polynomial in u[1], u[2], u[3]; 0 outside the stable range *)
sumU = u[1] + u[2] + u[3];

Print["=== A. Lemma 1: string equation (S) for the Hodge correlators with SYMBOLIC weights ==="];
tuples[g_, nmax_] := Join @@ Table[Tuples[Range[0, 3 g - 3 + n], n], {n, 1, nmax}];
strOK[g_, a_] := Together[corrU[g, Append[a, 0]] - Sum[If[a[[i]] >= 1, corrU[g, ReplacePart[a, i -> a[[i]] - 1]], 0], {i, Length[a]}]] === 0;
Print["(S) genus 1, N <= 4 points before adding tau_0: ", And @@ (strOK[1, #] & /@ tuples[1, 4])];
Print["(S) genus 2, N <= 3: ", And @@ (strOK[2, #] & /@ tuples[2, 3]), "   (", Round[AbsoluteTime[] - tt], " s)"];
Print["<tau_0>_2 = ", corrU[2, {0}], ";  <tau_0>_3 = ", corrU[3, {0}], ";  <tau_0>_4 = ", corrU[4, {0}]];
Print["anomaly: <tau_0>_1 + (1/24)(u1+u2+u3) = ", Together[corrU[1, {0}] + sumU/24]];
Print["CY check of A_1: ", Together[(corrU[1, {0}] /. usubOf[{1 + x, -x, -1}]) - (1 + x + x^2)/(24 (1 + x) x)] === 0];

Print["=== B. independent evaluation of Z[h,n] (ordered tuples), symbolic u and c ==="];
otuples[k_, lo_, hi_, smax_] := If[k == 0, {{}}, Join @@ Table[Prepend[#, a] & /@ otuples[k - 1, lo, hi, smax - a], {a, lo, Min[hi, smax]}]];
phiCoef[k_, cc_] := 1/cc^(k + 1);                       (* 1/(c - psi) = Sum psi^k / c^(k+1) *)
TCoef[d_, cc_] := If[d >= 2, -1/cc^(d - 1), 0];        (* -psi^2/(c - psi) = -Sum_{d>=2} psi^d / c^(d-1) *)
Zord[h_, n_, corr_, cc_] := Module[{tot = 0, dim},
   Do[dim = 3 h - 3 + n + m;
    If[2 h - 2 + n + m <= 0 || dim < 2 m, Continue[]];
    Do[Do[tot += (Times @@ (phiCoef[#, cc] & /@ A)) (Times @@ (TCoef[#, cc] & /@ D)) corr[h, Join[A, D]]/m!,
       {D, otuples[m, 2, dim - Total[A], dim - Total[A]]}], {A, otuples[n, 0, dim, dim]}], {m, 0, 3 h - 3 + n}];
   Together[tot]];
A1 = -sumU/24;
Do[Print["Z[1,", n, "] - (n-1)! c^(1-2n) A_1 = ", Together[Zord[1, n, corrU, c] - (n - 1)! c^(1 - 2 n) A1]], {n, 1, 4}];
Do[Print["Z[2,", n, "] (symbolic u, c) = ", Zord[2, n, corrU, c]], {n, 1, 3}];
Print["Z[2,0] (symbolic) = ", Factor[Zord[2, 0, corrU, c]]];
Print["Z[2,0] at CY weights vs -(7+10x^2+7x^4)/(5760 x^2 (1+x)^3): ", Together[(Zord[2, 0, corrU, c] /. usubOf[{1 + x, -x, -1}] /. c -> 1 + x) + (7 + 10 x^2 + 7 x^4)/(5760 x^2 (1 + x)^3)] === 0];
Print["genus-1 CY values: ", Table[Factor[Zord[1, n, corrU, c] /. usubOf[{1 + x, -x, -1}] /. c -> 1 + x], {n, 1, 3}]];
Print["   (", Round[AbsoluteTime[] - tt], " s)"];

Print["=== C. the key step of the Theorem: (S') with the series insertion t_s, symbolic s, c, u ==="];
(* coefficient of psi^k in t_s = c s phi + T, in 1 + (s - c) phi, in 1 (tau_0) *)
tsCoef[k_] := s c^(-k) - If[k >= 2, c^(1 - k), 0];
gCoef[k_] := If[k == 0, 1, 0] + (s - c) c^(-k - 1);
(* <A(psi_0) prod_{i=1}^M B(psi_i)>_h by ordered enumeration, Acoef/Bcoef the coefficient functions *)
corrIns[h_, Acoef_, Bcoef_, M_] := Module[{dim = 3 h - 3 + M + 1, tot = 0},
   If[2 h - 2 + M + 1 <= 0, Return[0]];
   Do[Do[tot += Acoef[a0] (Times @@ (Bcoef /@ B)) corrU[h, Prepend[B, a0]], {B, otuples[M, 0, dim - a0, dim - a0]}], {a0, 0, dim}];
   Together[tot]];
one[k_] := If[k == 0, 1, 0];
Do[Print["genus ", h, ": <tau_0 t_s^M>_{h,M+1} - M <(1+(s-c)phi) t_s^(M-1)>_{h,M}, M = 1..", If[h == 1, 5, 4], ": ",
   Table[Together[corrIns[h, one, tsCoef, M] - M corrIns[h, gCoef, tsCoef, M - 1]], {M, 1, If[h == 1, 5, 4]}]], {h, 1, 2}];
(* the telescoped identity itself: Psi_1(s) = <tau_0>_h + Psi_1(s) + (s-c) Psi_phi(s), i.e. (s - c) Psi_phi(s) + <tau_0>_h = 0, truncated consistently:
   Psi_phi(s) = Sum_M 1/M! <phi t_s^M>_{h,M+1}; only M <= 3h-2+2n contribute to s^n, so truncating at Mmax is exact up to s^((Mmax-3h+2)/2) *)
PsiPhi[h_, Mmax_] := Sum[corrIns[h, phiCoef[#, c] &, tsCoef, M]/M!, {M, 0, Mmax}];
Do[Module[{Mmax = If[h == 1, 7, 6], nex}, nex = Floor[(Mmax - 3 h + 2)/2];
   Print["genus ", h, ": (s-c) Psi_phi(s) + <tau_0>_h, exact to s^", nex, ": ", Normal[Series[Together[(s - c) PsiPhi[h, Mmax] + corrU[h, {0}]], {s, 0, nex}]]]], {h, 1, 2}];
Print["   (", Round[AbsoluteTime[] - tt], " s)"];

Print["=== D. Witten-Kontsevich vertex (u -> 0 is NOT the point vertex; use wk directly), symbolic c ==="];
wkCorr[h_, a_List] := If[2 h - 2 + Length[a] > 0 && Total[a] == 3 h - 3 + Length[a], wk[h, a], 0];
Print["Z^WK[1,n], n=1..3: ", Table[Zord[1, n, wkCorr, c], {n, 1, 3}], ";  Z^WK[2,n], n=1..2: ", Table[Zord[2, n, wkCorr, c], {n, 1, 2}], ";  Z^WK[2,0] = ", Zord[2, 0, wkCorr, c]];
Print["total time ", Round[AbsoluteTime[] - tt], " s"];
