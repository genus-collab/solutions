(* hv_checks.m -- checks of the Hodge-integral evaluator and of the scalar R-action (PPZ) up to genus GMAX.
   Usage:  wolframscript -file hv_checks.m GMAX      (GMAX = 3: fast; GMAX = 4 adds the 379-graph genus-4 PPZ check, minutes)
   Output: hv_checks.log with a VALIDATION line. *)
SetDirectory[DirectoryName[$InputFileName]];
GMAX = ToExpression[$ScriptCommandLine[[2]]];
Get["hv_lib4.m"];
log[m_] := Module[{s = OpenAppend["hv_checks.log"]}, WriteString[s, DateString[{"Year", "-", "Month", "-", "Day", " ", "Hour", ":", "Minute", ":", "Second"}], "  ", m, "\n"]; Close[s]; Print[m]];
tt = AbsoluteTime[];
{ok, rep} = hodgeChecks4[]; Do[log[r], {r, rep}];
log["hodgeChecks4 done in " <> ToString[Round[AbsoluteTime[] - tt]] <> " s"];
ppz2 = ppzScalar[2, {w1, w2, w3}]; c = Together[ppz2 - (hodgeNormU4[2, {}] /. usubOf[{w1, w2, w3}])] === 0; log["PPZ genus 2 symbolic = normalised Hodge class on M_2: " <> ToString[c]]; ok = ok && c;
SeedRandom[11];
Do[Module[{wn, gr, c2},
   gr = prepareGraph /@ graphsWithLegs[g, 0];
   Do[wn = Table[RandomInteger[{-9, 9}]/RandomInteger[{1, 5}], {3}]; If[MemberQ[wn, 0], Continue[]];
    tt = AbsoluteTime[]; c2 = ppzScalar[g, wn, gr] == (hodgeNormU4[g, {}] /. usubOf[wn]);
    log["PPZ genus " <> ToString[g] <> " (" <> ToString[Length[gr]] <> " graphs) at weights " <> ToString[wn] <> ": " <> ToString[c2] <> " (" <> ToString[Round[AbsoluteTime[] - tt]] <> " s)"]; ok = ok && c2, {If[g == 4, 4, 2]}]], {g, 3, GMAX}];
log[If[ok, "VALIDATION OK: Hodge integrals and scalar R-action up to genus " <> ToString[GMAX], "VALIDATION FAIL: Hodge integrals / scalar R-action"]];
