(* genus 2 (a few seconds per mode, everything exact in Q and x): point vertices with R_Z; Hodge vertices of mode B (weights
   eps1+eps2, -eps1, -eps2) and of mode A (the opposite sign), each with R' = R_Z h_w^-1 *)
res2 = Association[Table[m -> runMode[2, m, -1], {m, {"point", "B", "A"}}]];
Table[{m, res2[m]["counts"], "all values symmetric" -> res2[m]["symmetric"], "total = F20" -> res2[m]["totalOK"], "total polar part = gap" -> res2[m]["gapOK"]}, {m, Keys[res2]}] // TableForm
(* the opposite exponent s = +1 does not reproduce F20 (shown for mode B) *)
Module[{r = runMode[2, "B", 1]}, {r["totalOK"], Factor[r["total"] - FrefHVC[2, x]]}]
(* mode B per stable graph: only the genus-2 vertex graph is singular at Q = 1 and it carries exactly the gap *)
tableOf[res2["B"]]
{Factor /@ res2["B"]["stable"][[1]]["polar"], gapHVC[2, x]}
(* mode A and the point expansion: every stable graph is singular *)
tableOf[res2["A"]]
tableOf[res2["point"]]
(* mode B at the level of decorated graphs: the 7 decorations of the vertex carry the polar part, all 20 decorated graphs with an edge
   are regular individually; the structural statements *)
Grid[Prepend[Table[{e[[1]], e[[2, 2]], Factor /@ e[[5]]}, {e, Select[res2["B"]["entries"], #[[1, 1]] == 1 &]}], {"id", "dilaton leaves", "polar part {t^-2, t^-1}"}], Frame -> All]
{Count[Select[res2["B"]["entries"], #[[1, 1]] > 1 &], e_ /; regQ[e[[5]]]], structure[res2["B"]]}
(* the terms of the mode-B edge graphs are simple rational functions of Q and x *)
Table[{k, res2["B"]["stable"][[k]]["genera"], res2["B"]["stable"][[k]]["edges"], Factor[res2["B"]["stable"][[k]]["term"]]}, {k, 2, 7}] // TableForm
(* individual terms have a simple pole at Q = 0 (the canonical coordinates degenerate there) whose residues cancel; the sum at Q = 0
   is the constant-map term F_{2,0}(0) *)
{Factor /@ Table[Residue[res2["B"]["stable"][[k]]["term"], {Q, 0}], {k, 7}], Factor[res2["B"]["total"] /. Q -> 0], Factor[FrefHVC[2, x] /. Q -> 0]}
(* regression against the stored runs of hvc_run.wls *)
Table[If[FileExistsQ["hvc_results_g2_" <> t <> ".m"], Module[{st = Get["hvc_results_g2_" <> t <> ".m"], m = StringTake[t, 1]}, If[m == "p", m = "point"]; {t, And @@ Table[Together[res2[m]["stable"][[k]]["term"] - st["stable"][[k]]["term"]] === 0, {k, 7}]}], {t, "no stored file"}], {t, {"Bm", "Am", "point"}}]
