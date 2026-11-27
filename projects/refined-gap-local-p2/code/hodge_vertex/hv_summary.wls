(* hv_summary.wls -- human-readable summary of hodge_vertex_genus2_results.m *)
SetDirectory[DirectoryName[$InputFileName]];
res = Get["hodge_vertex_genus2_results.m"];
e1 = \[Epsilon][1]; e2 = \[Epsilon][2];
modes = Select[{"point", "A", "B"}, KeyExistsQ[res, #] &];
out = {};
emit[s_] := AppendTo[out, s];
fmt[expr_] := ToString[InputForm[expr /. {e1 -> E1, e2 -> E2}]];   (* E1, E2 = eps1, eps2 for readability *)
byXX[expr_] := Module[{t = Together[expr], k, c}, StringJoin[Table[c = Factor[Coefficient[Numerator[t], XX, k]/Denominator[t]]; If[c === 0, "", "    XX^" <> ToString[k] <> " : " <> fmt[c] <> "\n"], {k, 0, Exponent[Numerator[t], XX]}]]];
emit["F_{2,0} of K_P2 x C^2 (fibrewise torus T_F): contribution of each stable graph of genus 2 for three kinds of vertices"];
emit["  point : Witten-Kontsevich vertices, R-matrix R_Z (localp2CY5_genus3.nb conventions)"];
emit["  A     : Hodge CohFT of C^3, weights {-eps1-eps2, eps1, eps2}, R' = R_Z h^(" <> ToString[res["calibration"]["s"]] <> ")"];
emit["  B     : Hodge CohFT of C^3, weights {eps1+eps2, -eps1, -eps2}, R' = R_Z h^(" <> ToString[-res["calibration"]["s"]] <> ")"];
emit["Variables: y (1 + 27 y = 0 conifold), XX = y c1'(y)/c1(y), E1 = eps1, E2 = eps2. Polar parts: coefficients of t_c^-2, t_c^-1 with eps1 = x, eps2 = 1."];
emit["h_wA(z) = " <> fmt[Sum[res["hCoeffsA"][[k + 1]] z^k, {k, 0, Length[res["hCoeffsA"]] - 1}]]];
emit[""];
Do[emit["mode " <> m <> ": total = F20: " <> ToString[res[m]["totalOK"]] <> "; total polar part: " <> fmt[Factor /@ Total[res[m]["polar"]]]], {m, modes}];
emit[""];
Do[Module[{gr = res["graphs"][[k]]},
   emit["=== stable graph " <> ToString[k] <> ": genera " <> ToString[gr["genera"]] <> ", edges " <> ToString[gr["edges"]] <> " (adjacency " <> ToString[gr["adjacency"]] <> "), |Aut| = " <> ToString[gr["Aut"]]];
   Do[emit["  [" <> m <> "]  polar {t_c^-2, t_c^-1} = " <> fmt[Factor /@ res[m]["polar"][[k]]]];
    emit[byXX[res[m]["graphTerms"][[k]]]], {m, modes}];
   If[MemberQ[modes, "A"] && MemberQ[modes, "point"],
    emit["  [A - point]"]; emit[byXX[res["A"]["graphTerms"][[k]] - res["point"]["graphTerms"][[k]]]]];
   emit[""]], {k, Length[res["graphs"]]}];
(* decorated graphs *)
Do[emit["=== mode " <> m <> ": " <> ToString[Length[res[m]["orbitTerms"]]] <> " decorated graphs (orbits)"];
  Do[Module[{d = res[m]["descriptions"][id]},
    emit["  " <> ToString[id] <> "  genera " <> ToString[d["genera"]] <> ", edges " <> ToString[d["edges"]] <> ", psi on half-edges " <> ToString[d["psiOnHalfEdges"]] <> ", dilaton leaves " <> ToString[d["dilatonLeaves"]]];
    emit[byXX[res[m]["orbitTerms"][id]]]], {id, Keys[res[m]["orbitTerms"]]}];
  emit[""], {m, modes}];
Export["hodge_vertex_genus2_summary.txt", StringRiffle[out, "\n"], "Text"];
Print["written hodge_vertex_genus2_summary.txt (", Length[out], " lines)"];
