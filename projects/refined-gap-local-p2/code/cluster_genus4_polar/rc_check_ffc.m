(* rc_check_ffc.m -- check of the closed (Stirling) form of the QRR normalisation constants against the Q -> 0 limit at general weights.
   Usage:  wolframscript -file rc_check_ffc.m ZMAX     (ZMAX = 9 for genus 4; the limit method takes a few minutes per point at z^9)
   Output: rc_check_ffc.log with a VALIDATION line. *)
SetDirectory[DirectoryName[$InputFileName]];
zmax = ToExpression[$ScriptCommandLine[[2]]];
Get["rc_lib.m"];
log[m_] := Module[{s = OpenAppend["rc_check_ffc.log"]}, WriteString[s, DateString[{"Year", "-", "Month", "-", "Day", " ", "Hour", ":", "Minute", ":", "Second"}], "  ", m, "\n"]; Close[s]; Print[m]];
setupWick[zmax];
okall = True;
Do[Module[{t0 = AbsoluteTime[], SU, LeqU, lim, cl, diff},
   SU = (a0 + Q (a1 - a2))/(Q (a0 + a1)) /. ws; LeqU = Q (a0 + a1) U^2 - (a0 + Q (a1 - a2)) U - a2 /. ws;
   lim = ffcLimit[ws, SU, LeqU]; cl = ffcClosed[ws];
   diff = Table[Together[lim[[i]] - cl], {i, 1, 2}];
   If[Union[Flatten[diff]] =!= {0}, okall = False];
   log["weights " <> ToString[ws] <> ": limit minus closed form (roots 1, 2) = " <> ToString[diff] <> " (" <> ToString[Round[AbsoluteTime[] - t0, 0.1]] <> " s)"]],
  {ws, {{a0 -> 5, a1 -> 3, a2 -> -7, e1 -> 2, e2 -> 11}, {a0 -> 2 + 9 - 3 - 4, a1 -> 3, a2 -> 4, e1 -> 2, e2 -> 9}, {a0 -> -5 + 7 + 6 - 1, a1 -> 6, a2 -> -1, e1 -> -5, e2 -> 7}}}];
log[If[okall, "VALIDATION OK: closed-form QRR normalisation constants to z^" <> ToString[zmax], "VALIDATION FAIL: closed-form QRR normalisation constants"]];
