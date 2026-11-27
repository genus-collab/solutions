(* lp2_checks.m -- checks on the reconstructed local P2 free energies up to genus GMAX (port of the Checks section of
   localp2CY5_genus3.nb): Q-expansions, unrefined GV invariants against the literature table, refined BPS numbers at low
   degree, the modular anomaly equation (normalisation fixed at genus 2), the conifold expansion (gap), the refined split
   (5.54) and the unrefined / NS limits.
   Usage:  wolframscript -file lp2_checks.m GMAX
   Input:  ref/F20_F30.m (F20, F30 in y, XX, eps), localP2_F40.m (if GMAX = 4; from lp2_collect.m). Output: lp2_checks_g<GMAX>.log *)
SetDirectory[DirectoryName[$InputFileName]];
GMAX = ToExpression[$ScriptCommandLine[[2]]];
log[m_] := Module[{s = OpenAppend["lp2_checks_g" <> ToString[GMAX] <> ".log"]}, WriteString[s, DateString[{"Hour", ":", "Minute", ":", "Second"}], "  ", m, "\n"]; Close[s]; Print[m]];
e1 = \[Epsilon][1]; e2 = \[Epsilon][2];
reflist = Get["ref/F20_F30_y.m"];   (* {F20, F30} in y, XX, eps (from F20_result.mx, F30_result.mx) *)
F[2] = reflist[[1]]; F[3] = reflist[[2]];
If[GMAX >= 4, F[4] = Get["localP2_F40.m"]];
log["loaded F_{g,0} for g = 2.." <> ToString[GMAX] <> "; LeafCounts " <> ToString[Table[LeafCount[F[g]], {g, 2, GMAX}]]];
F11 = -(4 + (e1^2 + e2^2)/(e1 e2 (1 + 27 y)) + 12 XX)/(24 cc);
C1[y_] = Hypergeometric2F1[1/3, 2/3, 1, -27 y];
NQ = 10;
II10 = Log[y] + Sum[3 (3 d - 1)!/(d!)^3 (-y)^d, {d, 1, 2 NQ}];
invmm = InverseSeries[Series[Exp[II10], {y, 0, 2 NQ}] /. y -> Q];
F03 = -(1/(3 (1 + 27 y) cc^3 e1 e2));
XXser = Series[y C1'[y]/C1[y], {y, 0, 2 NQ}];
Qser[f_] := Series[Normal[Series[f /. XX -> XXser /. cc -> C1[y], {y, 0, 2 NQ}]] /. y -> invmm, {Q, 0, NQ}];
F03Q = Qser[F03];
F0Q = Sum[SeriesCoefficient[F03Q, d]/d^3 Q^d, {d, 1, NQ}];
F1Q = Sum[SeriesCoefficient[Qser[F11], d]/d Q^d, {d, 1, NQ}];
FgQ[0] = F0Q; FgQ[1] = F1Q;
Do[tt = AbsoluteTime[]; FgQ[g] = Normal[Qser[F[g]]]; log["Q-expansion of F" <> ToString[g] <> " to Q^" <> ToString[NQ] <> " done (" <> ToString[Round[AbsoluteTime[] - tt, 0.1]] <> " s)"], {g, 2, GMAX}];
Do[log["F" <> ToString[g] <> " Q-series (eps1 eps2 times, first 4): " <> ToString[InputForm[Simplify[Table[SeriesCoefficient[FgQ[g] e1 e2, {Q, 0, d}], {d, 1, 4}]]]]], {g, 2, GMAX}];
(* ---------- unrefined GV invariants ---------- *)
Funref[g_] := If[g == 0, FgQ[0] e1 e2 /. e2 -> -e1 /. e1 -> 1 // Simplify, If[g == 1, Simplify[(FgQ[1] /. e2 -> -e1)], Simplify[(FgQ[g] /. e2 -> -e1)/(-e1^2)^(g - 1)]]];
Do[Fu[g] = Funref[g], {g, 0, GMAX}];
gvSeries[lam_] := Sum[n[g, d]/k (2 Sin[k lam/2])^(2 g - 2) Q^(d k), {g, 0, GMAX}, {d, 1, NQ}, {k, 1, Floor[NQ/d]}];
eqs = Flatten[Table[SeriesCoefficient[Series[gvSeries[lam], {lam, 0, 2 GMAX - 2}], {lam, 0, 2 g - 2}] - Fu[g] // Expand // CoefficientList[#, Q] &, {g, 0, GMAX}]];
nsol = Solve[Thread[DeleteCases[eqs, 0] == 0], Flatten[Table[n[g, d], {g, 0, GMAX}, {d, 1, NQ}]]];
gvtab = Table[n[g, d], {g, 0, GMAX}, {d, 1, NQ}] /. First[nsol];
log["unrefined GV invariants n^g_d (rows g = 0.." <> ToString[GMAX] <> ", d = 1.." <> ToString[NQ] <> "): " <> ToString[gvtab]];
(* literature table for local P2 (Katz-Klemm-Vafa / Aganagic-Mari\~no-Vafa), d = 1..7 *)
lit = {{3, -6, 27, -192, 1695, -17064, 188454}, {0, 0, -10, 231, -4452, 80948, -1438086}, {0, 0, 0, -102, 5430, -194022, 5784837}, {0, 0, 0, 15, -3672, 290853, -15363990}, {0, 0, 0, 0, 1386, -290400, 29056614}};
Do[log["genus " <> ToString[g] <> " GV d = 1..7 versus the literature table: " <> If[gvtab[[g + 1, 1 ;; 7]] === lit[[g + 1]], "agree", "DIFFER: computed " <> ToString[gvtab[[g + 1, 1 ;; 7]]] <> " table " <> ToString[lit[[g + 1]]]]], {g, 0, Min[GMAX, 4]}];
(* ---------- refined BPS numbers ---------- *)
chi[j_, q_] := Sum[q^(2 m), {m, -j, j}];
GWtot = Sum[FgQ[g], {g, 0, GMAX}]; GWtot = Normal[GWtot] - (Normal[GWtot] /. Q -> 0);
spins[dd_] := Select[Flatten[Table[{jL, jR}, {jL, 0, (dd - 1) (dd - 2)/4, 1/2}, {jR, 0, dd (dd + 3)/4, 1/2}], 1], Mod[2 #[[1]] + 2 #[[2]] - dd - 1, 2] == 0 &];
Do[Module[{bpsAnsatz, diffser, beqs, bunknowns, bsol},
  bpsAnsatz = Sum[1/k Sum[(-1)^(2 sp[[1]] + 2 sp[[2]]) NN[dd, sp[[1]], sp[[2]]] chi[sp[[1]], Exp[k (e1 - e2)/2]] chi[sp[[2]], Exp[k (e1 + e2)/2]], {sp, spins[dd]}]/((Exp[k e1/2] - Exp[-k e1/2]) (Exp[k e2/2] - Exp[-k e2/2])) Q^(dd k), {dd, 1, DREF}, {k, 1, Floor[DREF/dd]}];
  diffser = Series[(bpsAnsatz - Normal[GWtot] /. Q^p_ /; p > DREF -> 0) e1 e2 /. e2 -> s x2 /. e1 -> s x1, {s, 0, 2 GMAX}];
  beqs = DeleteCases[Flatten[Table[CoefficientList[Expand[SeriesCoefficient[diffser, {s, 0, m}]], {Q, x1, x2}], {m, 0, 2 GMAX}]], 0];
  bunknowns = Flatten[Table[NN[dd, sp[[1]], sp[[2]]], {dd, 1, DREF}, {sp, spins[dd]}]];
  bsol = Solve[Thread[beqs == 0], bunknowns];
  If[bsol === {}, log["refined BPS: INCONSISTENT at d <= " <> ToString[DREF]],
   log["refined BPS numbers determined (nonzero), d <= " <> ToString[DREF] <> ": " <> ToString[InputForm[Select[Thread[bunknowns -> (bunknowns /. First[bsol])], #[[2]] =!= 0 && NumberQ[#[[2]]] &]]]];
   log["   undetermined at this genus: " <> ToString[Length[Select[bunknowns /. First[bsol], ! NumberQ[#] &]]]]]], {DREF, 1, 4}];
(* ---------- modular anomaly equation ---------- *)
QdQ[f_] := Together[(y/cc) (D[f, y] + D[f, XX] (-(XX^2 + (27 XX y + 6 y)/(1 + 27 y))/y) + D[f, cc] (XX cc/y))];
rhs[2] = cc^2 (QdQ[F11] + F11^2);
ratio2 = Together[D[F[2], XX]/rhs[2]];
log["genus 2 anomaly: d F2/dXX / (cc^2 [QdQ F11 + F11^2]) = " <> ToString[InputForm[ratio2]]];
If[GMAX >= 3, rhs[3] = cc^2 (QdQ[QdQ[F[2]]] + 2 F11 QdQ[F[2]]);
  log["genus 3 anomaly holds with the genus-2 normalisation: " <> ToString[Together[D[F[3], XX] - ratio2 rhs[3]] === 0]]];
If[GMAX >= 4, rhs[4] = cc^2 (QdQ[QdQ[F[3]]] + 2 F11 QdQ[F[3]] + QdQ[F[2]]^2);
  log["genus 4 anomaly holds with the genus-2 normalisation: " <> ToString[Together[D[F[4], XX] - ratio2 rhs[4]] === 0]]];
(* ---------- conifold expansion ---------- *)
Get["lp2_lib.m"]; tcSetup[2 GMAX + 4];
cfexp[f_] := Series[f /. XX -> XXcf /. cc -> C1cf /. y -> ycf, {tcv, 0, 1}];
Do[Module[{s = cfexp[F[g]]}, log["genus " <> ToString[g] <> " conifold expansion t_c^-" <> ToString[2 g - 2] <> " .. t_c^1: " <> ToString[InputForm[Table[Factor[SeriesCoefficient[s, n]], {n, -(2 g - 2), 1}]]]];
   log["   gap and leading coefficient (-3)^(g-1)(2g-3)! [t^2g] t^2/(4 sinh sinh): " <> ToString[Union[Table[Together[SeriesCoefficient[s, n]], {n, -(2 g - 3), -1}]] === {0} && Together[SeriesCoefficient[s, -(2 g - 2)] - gapCoefficient[g]] === 0]]], {g, 2, GMAX}];
(* ---------- refined split (5.54) ---------- *)
Do[Module[{poly, red, split},
   poly = Expand[F[g] e1 e2]; red = SymmetricReduction[poly, {e1, e2}, {ss, pp}][[1]];
   split = Table[Coefficient[Coefficient[red, ss, 2 g1], pp, g - g1], {g1, 0, g}];
   Do[log["F_{" <> ToString[g1] <> "," <> ToString[g - g1] <> "} = " <> ToString[InputForm[Collect[Factor[split[[g1 + 1]] /. y -> (1 - LL^3)/(27 LL^3)], XX, Factor]]]], {g1, 0, g}]], {g, 2, GMAX}];
log["unrefined limit (eps2 -> -eps1) of eps1 eps2 F_g: " <> ToString[InputForm[Table[Factor[F[g] e1 e2 /. e2 -> -e1], {g, 2, GMAX}]]]];
log["NS limit (eps2 -> 0) of eps1 eps2 F_g: " <> ToString[InputForm[Table[Factor[Limit[F[g] e1 e2, e2 -> 0]], {g, 2, GMAX}]]]];
log["done"];
