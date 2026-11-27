(* hvc_vertex_identity.wls -- the leading-order vertex functional at the conifold point (mode B, singular idempotent):
   leaves on half-edges phi(psi) = 1/(c - psi), dilaton leaves T(psi) = -psi^2/(c - psi), c = eps1 + eps2 (= 1 + x, eps2 = 1), Hodge class
   with weights {c, -x, -1}.  Z[h, n] = Sum_m 1/m! Int_{M_{h,n+m}} Hodge prod_{i<=n} phi(psi_i) prod_{j<=m} T(psi_j).
   Checks: (i) edge terms at t = 0 factorise as E11 = -(x^3/2) phi(z) phi(w) (and the form of E12, E22); (ii) Z[h,0] times the leading
   coefficient (-x c^3)^(h-1) of Delta_1^(h-1) reproduces the gap at genus h; (iii) Z[h, n] for n >= 1. *)
SetDirectory[DirectoryName[$InputFileName]];
Get["hvc_lib.m"]; Get["hvc_factors_g3_Bm.mx"];
S0 = Together[SsumQ]; P0 = Together[PprodQ]; sq = ((1 + x)/x) Sqrt[1 + 4 x (1/Q - 1)/(1 + x)^2]; U1 = (S0 + sq)/2; U2 = (S0 - sq)/2;
val0[C[a_, b_], root_] := Factor[SeriesCoefficient[Series[(a + b root) /. Q -> Exp[t], {t, 0, 1}], 0]];
c = 1 + x;
Print["E11[m,n](t=0) = -(x^3/2)/c^(m+n+2) for m,n <= 4 (m+n <= 5): ", And @@ Flatten[Table[If[m + n <= 5, Together[val0[EdN[1, 1][m, n], U1] + (x^3/2)/c^(m + n + 2)] === 0, True], {m, 0, 4}, {n, 0, 4}]]];
Print["E12[m,n](t=0) (1 = singular, 2 = regular), m,n = 0..2: ", Table[val0[EdN[1, 2][m, n], U1], {m, 0, 2}, {n, 0, 2}]];
Print["E12[m,n](t=0) / E12[0,n](t=0) (z-dependence of the singular leg): ", Table[Factor[val0[EdN[1, 2][m, n], U1]/val0[EdN[1, 2][0, n], U1]], {m, 0, 3}, {n, 0, 2}]];
(* every stored factor is read at U1: rc_lib.m stores the idempotent-2 data as cconj of the idempotent-1 data, so EdN[2,2] read at U2 is the E11 table (erratum 2026-10-31, leaf_proof/leaf_structure_proof.md Remark 8.8) *)
Print["E22[m,n](t=0), m,n = 0..2 (not of rank one; = (x^3/c)[z^m w^n](1 - h_w(z)^2 h_w(w)^2)/(z + w)): ", Table[val0[EdN[2, 2][m, n], U1], {m, 0, 2}, {n, 0, 2}]];
Print["Delta_1 leading coefficient (t^-2): ", Factor[SeriesCoefficient[Series[(DeltaCN[1][[1]] + DeltaCN[1][[2]] U1) /. Q -> Exp[t], {t, 0, 0}], -2]]];
(* the vertex functional *)
usub = usubOf[{c, -x, -1}];
phiC[k_] := 1/c^(k + 1);            (* coefficient of psi^k in 1/(c - psi) *)
TC[d_] := -1/c^(d - 1);             (* coefficient of psi^d (d >= 2) in -psi^2/(c - psi) *)
Z[h_, n_] := Z[h, n] = Module[{tot = 0, dim, tuples},
   Do[dim = 3 h - 3 + n + m;
    If[dim < 2 m, Continue[]];
    tuples = Select[Tuples[Join[Table[Range[0, dim], {n}], Table[Range[2, dim], {m}]]], Total[#] <= dim &];
    Do[tot += (Times @@ (phiC /@ tp[[;; n]])) (Times @@ (TC /@ tp[[n + 1 ;;]])) (hodgeNormU4[h, tp] /. usub), {tp, tuples}];
    , {m, 0, 3 h - 3 + n}];
   Together[tot/1]];
(* note: 1/m! for the dilaton leaves is supplied by summing over ordered tuples divided by m! -- do it explicitly *)
Zm[h_, n_] := Zm[h, n] = Module[{tot = 0, dim, tuples},
   Do[dim = 3 h - 3 + n + m;
    If[dim < 2 m, Continue[]];
    tuples = Select[Tuples[Join[Table[Range[0, dim], {n}], Table[Range[2, dim], {m}]]], Total[#] <= dim &];
    tot += Sum[(Times @@ (phiC /@ tp[[;; n]])) (Times @@ (TC /@ tp[[n + 1 ;;]])) (hodgeNormU4[h, tp] /. usub), {tp, tuples}]/m!, {m, 0, 3 h - 3 + n}];
   Together[tot]];
Do[Print["Z[", h, ", 0] = ", Factor[Zm[h, 0]], ";  (-x c^3)^(h-1) Z[h,0] = ", Factor[(-x c^3)^(h - 1) Zm[h, 0]], ";  gap_h leading coefficient = ", gapHVC[h, x][[1]], ";  equal: ", Together[(-x c^3)^(h - 1) Zm[h, 0] - gapHVC[h, x][[1]]] === 0], {h, 2, 3}];
Do[Print["Z[", h, ", ", n, "] = ", Factor[Zm[h, n]]], {h, 1, 3}, {n, 1, If[h == 3, 2, 3]}];
Print["genus 0: Z[0, n] for n = 3, 4: ", {Factor[Zm[0, 3]], Factor[Zm[0, 4]]}];
