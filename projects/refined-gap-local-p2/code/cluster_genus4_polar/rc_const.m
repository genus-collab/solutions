(* rc_const.m -- degree-zero (constant map) term F_{g,0}(Q = 0) of the resolved conifold, symbolic in the five weights.
   Usage:  wolframscript -file rc_const.m genus ZMAX
   Output: rc_const_g<genus>.m (the rational function of alpha0, alpha1, alpha2, eps1, eps2); log in rc_const.log.
   For genus 2 and 3 the result is compared with the Q = 0 part of ref/resolvedconifold_results.m. *)
SetDirectory[DirectoryName[$InputFileName]];
{gg, zmax} = ToExpression /@ $ScriptCommandLine[[2 ;; 3]];
Get["rc_lib.m"];
log[m_] := Module[{s = OpenAppend["rc_const.log"]}, WriteString[s, DateString[{"Year", "-", "Month", "-", "Day", " ", "Hour", ":", "Minute", ":", "Second"}], "  [g=", gg, "] ", m, "\n"]; Close[s]; Print[m]];
t0 = AbsoluteTime[];
setupWick[zmax];
F0 = const0[gg];
log["F_{" <> ToString[gg] <> ",0}(Q = 0) computed (" <> ToString[Round[AbsoluteTime[] - t0, 0.1]] <> " s), LeafCount " <> ToString[LeafCount[F0]]];
Put[F0, "rc_const_g" <> ToString[gg] <> ".m"];
If[gg <= 3 && FileExistsQ["ref/resolvedconifold_results.m"],
  Module[{ref = Get["ref/resolvedconifold_results.m"][[gg]] /. Q -> 0},
   log[If[Together[F0 - ref] === 0, "VALIDATION OK: degree-zero term agrees with the stored genus-" <> ToString[gg] <> " result", "VALIDATION FAIL: degree-zero term differs from the stored result"]]]];
