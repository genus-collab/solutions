(* Headless driver for localp2CY5_genRmat_genus3.nb: evaluates its Input cells in order.
   Usage: wolframscript -file run_genRmat_genus3.wls [N]   (N = evaluate only the first N Input cells; default all) *)
SetDirectory[DirectoryName[$InputFileName]];
nb = Get["localp2CY5_genRmat_genus3.nb"];
cells = Cases[nb, Cell[s_String, "Input", ___] :> s, Infinity];
nmax = If[Length[$ScriptCommandLine] >= 2, ToExpression[$ScriptCommandLine[[2]]], Length[cells]];
Print["evaluating ", nmax, " of ", Length[cells], " Input cells of localp2CY5_genRmat_genus3.nb"];
Do[ToExpression[cells[[k]]], {k, nmax}];
Print["done"];
