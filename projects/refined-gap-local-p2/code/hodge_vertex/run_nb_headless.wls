(* run_nb_headless.wls -- evaluate the Input cells of a notebook in order with wolframscript (no front end), as a test that the
   notebook runs. Usage: wolframscript -file run_nb_headless.wls <notebook.nb> [firstCell lastCell]
   NotebookDirectory[] is stubbed to the notebook's folder. *)
args = Rest[$ScriptCommandLine];
nbfile = args[[1]];
nbdir = DirectoryName[AbsoluteFileName[nbfile]];
Unprotect[NotebookDirectory]; NotebookDirectory[] := nbdir; Protect[NotebookDirectory];
nb = Get[nbfile];
inputs = Cases[nb, Cell[s_String, "Input"] :> s, Infinity];
{i1, i2} = If[Length[args] >= 3, ToExpression /@ args[[2 ;; 3]], {1, Length[inputs]}];
Print[Length[inputs], " Input cells; evaluating ", i1, "..", i2];
Do[Module[{tt = AbsoluteTime[], res},
   Print["\n=== Input cell ", k, " ===\n", StringTake[inputs[[k]], Min[160, StringLength[inputs[[k]]]]], If[StringLength[inputs[[k]]] > 160, " ...", ""]];
   res = ToExpression[inputs[[k]]];
   If[res =!= Null, Print["--> ", ToString[res, InputForm] // StringTake[#, Min[3000, StringLength[#]]] &]];
   Print["    (", Round[AbsoluteTime[] - tt, 0.1], " s)"]], {k, i1, i2}];
Print["\nheadless run finished"];
