(* mode B, genus 3: the terms of the stable graphs with three or four vertices as rational functions of Q and x *)
Table[With[{st = res3["B"]["stable"][[k]]}, {k, st["genera"], st["edges"], Factor[st["term"]]}], {k, Select[Range[42], Length[res3["B"]["stable"][[#]]["genera"]] >= 3 &]}] // TableForm
(* denominators of the mode-B stable-graph terms at genus 3 from FactorList: powers of (1-Q), Q, the discriminant Q(1-x)^2+4x (the two
   idempotents collide at Q_c = -4x/(1-x)^2: every stable graph has such a pole, which cancels only in the sum over all graphs), x, (1+x),
   and whatever else remains (1 = nothing) *)
denPowers[t_] := Module[{fl = FactorList[Denominator[Together[t]]], known = {1 - Q, Q, Q (1 - x)^2 + 4 x, x, 1 + x}, pw, res = {}, hit}, pw = ConstantArray[0, 5];
   Do[If[NumberQ[fl[[i, 1]]], Continue[]]; hit = Position[known, f_ /; NumberQ[Together[fl[[i, 1]]/f]], {1}, 1, Heads -> False]; If[hit === {}, AppendTo[res, fl[[i]]], pw[[hit[[1, 1]]]] += fl[[i, 2]]], {i, Length[fl]}];
   Append[pw, If[res === {}, 1, Times @@ (Power @@@ res)]]];
{Sort[Tally[denPowers /@ res3["B"]["stable"][[All, "term"]]]], "total:", denPowers[res3["B"]["total"]]}
(* residues at Q = 0 of the stable-graph terms cancel; the total at Q = 0 is F_{3,0}(0) *)
{Factor[Total[Table[Residue[res3["B"]["stable"][[k]]["term"], {Q, 0}], {k, 42}]]], Factor[res3["B"]["total"] /. Q -> 0], Factor[FrefHVC[3, x] /. Q -> 0]}
