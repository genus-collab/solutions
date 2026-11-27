(* E_crosscheck_g4r2.wls -- requested check of Route E: the symbolic two- and three-leaf functionals Xi[4,2] and Xi[3,3]
   (hvc_mp_header.m conventions: c-normalisation, symbol x, leaves 1/(c - jj[i] psi), dilaton leaves -psi^2/(c - psi); NO TR),
   exported in InputForm so that compare_g4r2.pl can compare them with the Perl engine (engine_h4_r2.txt, engine_h3_r3.txt),
   which is the data behind the shortened-Norlund and Newton-regular structures of the Route E report.
   Usage (from this folder): wolframscript -file E_crosscheck_g4r2.wls  ->  E_crosscheck_g4r2_out.txt  ; then perl compare_g4r2.pl *)
SetDirectory[ParentDirectory[ParentDirectory[DirectoryName[$InputFileName]]]];   (* HVC, where hvc_mp_header.m and hvc_lib.m live *)
Get["hvc_mp_header.m"];
t0 = AbsoluteTime[];
lines = {};
Do[Module[{e = Xi @@ hr}, AppendTo[lines, "genus " <> ToString[hr[[1]]] <> " r=" <> ToString[hr[[2]]] <> " symbolic: " <> ToString[InputForm[Together[e]]]];
   Print["done ", hr, " [", Round[AbsoluteTime[] - t0], " s]"]], {hr, {{4, 2}, {3, 3}}}];
Export[FileNameJoin[{DirectoryName[$InputFileName], "E_crosscheck_g4r2_out.txt"}], StringRiffle[lines, "\n"], "Text"];
