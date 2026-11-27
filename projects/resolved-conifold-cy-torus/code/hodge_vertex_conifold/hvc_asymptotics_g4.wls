(* hvc_asymptotics.wls -- Laurent behaviour at the conifold point Q = 1 (t = Log Q) of the building blocks of the graph sum (genus-3 factors,
   z^6): rows of R, Delta, dilaton leaves, edge terms, at the singular root U_1 -> 1 and the regular root U_2 -> -1/x, with the square
   root of the discriminant written so that Series expands it cleanly.  Usage: wolframscript -file hvc_asymptotics.wls <Bm|Am|point> *)
SetDirectory[DirectoryName[$InputFileName]];
tag = $ScriptCommandLine[[2]];
Get["hvc_lib.m"]; Get["hvc_factors_g4_" <> tag <> ".mx"];
S0 = Together[SsumQ]; P0 = Together[PprodQ];
(* S^2 - 4P = ((1+x)^2 + 4x (1/Q - 1))/x^2 ; take the branch with sqrt -> (1+x)/x at Q = 1 *)
sq = ((1 + x)/x) Sqrt[1 + 4 x (1/Q - 1)/(1 + x)^2];
Print["check: sq^2 - (S^2 - 4P) = ", Simplify[sq^2 - (S0^2 - 4 P0)]];
U1 = (S0 + sq)/2; U2 = (S0 - sq)/2;
ser[C[a_, b_], root_, n_: 6] := Series[(a + b root) /. Q -> Exp[t], {t, 0, n}];
ordOf[s_] := Module[{nn = Normal[s]}, If[nn === 0 || Together[nn] === 0, "zero", Exponent[Together[nn] t^30, t, Min] - 30]];
lead[s_] := Module[{o = ordOf[s]}, If[o === "zero", 0, Factor[SeriesCoefficient[s, o]]]];
val0[s_] := Factor[SeriesCoefficient[s, 0]];
Print["=== frame ", tag, " ===  roots at Q = 1: ", {U1, U2} /. Q -> 1 // Simplify];
Do[Module[{s = ser[DeltaCN[1], r, 4]}, Print["Delta at ", If[r === U1, "singular", "regular"], " root: order ", ordOf[s], ", leading coefficient ", lead[s]]], {r, {U1, U2}}];
Do[Module[{ss = Table[ser[ratCN[a, 1][m], r, 4], {m, 0, 9}]},
   Print["row ", a, " at the ", If[r === U1, "singular", "regular"], " root: orders ", ordOf /@ ss, "; leading coefficients ", lead /@ ss]], {a, 0, 1}, {r, {U1, U2}}];
Print["row 0 at the singular root at t = 0, versus (-1/(1+x))^m: ", Table[Factor[val0[ser[ratCN[0, 1][m], U1, 2]] - (-1/(1 + x))^m], {m, 0, 9}]];
Print["row 0 at the regular root at t = 0: ", Table[val0[ser[ratCN[0, 1][m], U2, 2]], {m, 0, 9}]];
Print["dilaton leaves Dil[1][n] (n = 2..6) orders at the singular root: ", Table[ordOf[ser[DilN[1][n], U1, 4]], {n, 2, 9}]];
(* every stored factor is read at U1: rc_lib.m stores the idempotent-2 data as cconj of the idempotent-1 data, so EdN[2,2] read at U2 would print the orders of E11 (erratum 2026-10-31, leaf_proof/leaf_structure_proof.md Remark 8.8) *)
Do[Module[{rr = {U1, U1, U1}[[Position[{{1, 1}, {1, 2}, {2, 2}}, ij][[1, 1]]]]},
   Print["edge terms Ed", ij, " (root labelling 1 = singular): orders m, n = 0..2: ", Table[ordOf[ser[EdN[ij[[1]], ij[[2]]][m, n], U1, 4]], {m, 0, 4}, {n, 0, 4}],
     If[ij === {1, 1}, "; leading coefficients (m,n = 0..1): " <> ToString[InputForm[Table[lead[ser[EdN[1, 1][m, n], U1, 4]], {m, 0, 1}, {n, 0, 1}]]], ""]]], {ij, {{1, 1}, {1, 2}, {2, 2}}}];
Print["L (canonical-coordinate component, row 1 z^0) at the two roots: ", Table[ordOf[ser[ratCN[1, 1][0], r, 3]], {r, {U1, U2}}], "; leading: ", Table[lead[ser[ratCN[1, 1][0], r, 3]], {r, {U1, U2}}]];
