SetDirectory[DirectoryName[$InputFileName]];
r = Get["hodge_vertex_genus3_results.m"]; e1 = \[Epsilon][1]; e2 = \[Epsilon][2];
ab[t_] := Module[{d = Denominator[Together[t]]}, {Exponent[d, e1], Exponent[d /. e2 -> 1 /. e1 -> (s - 1), s] - Exponent[d /. e2 -> 1, e1] + 0}];
(* exponent of (e1+e2): substitute e2 -> 1, e1 -> s - 1 so that e1 + e2 = s; count the power of s in the denominator *)
bexp[t_] := Module[{d = Factor[Denominator[Together[t]]]}, Exponent[d /. {e2 -> 1, e1 -> s - 1} // Factor, s] - 0];
res = Table[Module[{d = Factor[Denominator[Together[r["orbitTerms"][id]]]]}, {id, Exponent[d, e1], Exponent[Factor[d /. {e2 -> 1, e1 -> s - 1}], s]}], {id, Keys[r["orbitTerms"]]}];
Print["genus 3 mode B: max power of eps1 in denominators: ", Max[res[[All, 2]]], "; max power of (eps1+eps2): ", Max[res[[All, 3]]]];
Print["distribution {a, b} -> count: ", Tally[res[[All, 2 ;; 3]]]];
