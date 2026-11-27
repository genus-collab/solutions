(* rc_cy_debug.wls -- at one CY-torus point and genus 3: polar parts of every decorated graph from the Laurent-series layer
   versus the exact rational-function route; lists the orbits that differ and tests whether a higher series order cures them. *)
SetDirectory[DirectoryName[$InputFileName]];
Get["rc_lib.m"];
gg = 3; setupWick[6]; loadOrbits[gg];
p = {5, 6, -4}; wsub = {a0 -> p[[3]] + 1 - p[[1]] - p[[2]], a1 -> p[[1]], a2 -> p[[2]], e1 -> p[[3]], e2 -> 1};
Print["point ", p, " weights ", wsub];
t0 = AbsoluteTime[]; rmatrixAtClosed[wsub]; Print["R-matrix: ", Round[AbsoluteTime[] - t0, 0.1], " s"];
t0 = AbsoluteTime[]; evalQ[Q];
exact = Association[Table[id -> polarOfRational[Together[graphValueC[OT["graphs"][id]][[1]]], gg], {id, OT["ids"]}]];
Print["exact route: ", Round[AbsoluteTime[] - t0, 0.1], " s; total minus P_3 Li_{-3} polar part: ", Together /@ (Total[Values[exact]] - (cyPolar[gg] /. wsub))];
Do[t0 = AbsoluteTime[]; evalQseries[M0];
   ser = Association[Table[id -> polarOfSeries[graphValueC[OT["graphs"][id]][[1]], gg], {id, OT["ids"]}]];
   diffs = Select[OT["ids"], exact[#] =!= ser[#] &];
   Print["series order ", M0, ": ", Round[AbsoluteTime[] - t0, 0.1], " s; total minus expected: ", Together /@ (Total[Values[ser]] - (cyPolar[gg] /. wsub)), "; orbits differing from exact: ", Length[diffs]];
   Do[Print["   ", id, " ", OT["orbits"][[id[[1]], id[[2]]]], " genera ", OT["descriptions"][id]["genera"], " edges ", OT["descriptions"][id]["edges"], "\n      exact  ", exact[id], "\n      series ", ser[id]], {id, Take[diffs, Min[4, Length[diffs]]]}], {M0, {12, 20}}];
