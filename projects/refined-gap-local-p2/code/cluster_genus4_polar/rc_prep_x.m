(* rc_prep_x.m -- conifold R-matrix and graph factors at alpha0 = 0, alpha1 = eps1 = x, alpha2 = eps2 = 1, with Q and x symbolic.
   Usage:  wolframscript -file rc_prep_x.m genus ZMAX        (ZMAX >= 3 genus - 3)
   Output: rc_factors_g<genus>_x.mx with the factors rN, DeltaN, DilN, LfN, EdN (C-elements with coefficients rational in Q, x),
           Ssum, Pprod (= SsumQ, PprodQ), etainvCur, genSym, ZMAX; timings in rc_prep_x.log.
   No QRR normalisation constants are needed at this specialisation (they vanish identically there). *)
SetDirectory[DirectoryName[$InputFileName]];
{gg, zmax} = ToExpression /@ $ScriptCommandLine[[2 ;; 3]];
Get["rc_lib.m"];
log[m_] := Module[{s = OpenAppend["rc_prep_x.log"]}, WriteString[s, DateString[{"Year", "-", "Month", "-", "Day", " ", "Hour", ":", "Minute", ":", "Second"}], "  [g=", gg, "] ", m, "\n"]; Close[s]; Print[m]];
t0 = AbsoluteTime[];
setupWick[zmax];
log["Wick expansion to z^" <> ToString[zmax] <> " ready (" <> ToString[Round[AbsoluteTime[] - t0, 0.1]] <> " s)"];
t1 = AbsoluteTime[];
rmatrixAt0[{a1 -> x, a2 -> 1, e1 -> x, e2 -> 1}];
log["R-matrix at alpha0 = 0 with x symbolic built (" <> ToString[Round[AbsoluteTime[] - t1, 0.1]] <> " s); S = " <> ToString[InputForm[SsumQ]] <> ", P = " <> ToString[InputForm[PprodQ]]];
t2 = AbsoluteTime[];
evalQ[Q];
log["factors at symbolic Q ready (" <> ToString[Round[AbsoluteTime[] - t2, 0.1]] <> " s); LeafCount of EdN[1,2][" <> ToString[zmax - 1] <> ",0]: " <> ToString[LeafCount[EdN[1, 2][zmax - 1, 0]]]];
DumpSave["rc_factors_g" <> ToString[gg] <> "_x.mx", {rN, DeltaN, DilN, LfN, EdN, ratCN, DeltaCN, Ssum, Pprod, SsumQ, PprodQ, etainvCur, genSym, ZMAX}];
log["saved rc_factors_g" <> ToString[gg] <> "_x.mx"];
