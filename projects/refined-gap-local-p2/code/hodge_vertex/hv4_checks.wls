(* hv4_checks.wls -- checks of the genus-general Hodge evaluator (hv_lib4.m): known values and relations, agreement with the
   closed-form implementation of hv_lib3.m on all genus-2 and genus-3 decoration weights, PPZ scalar checks (genus 2 symbolic,
   genus 3 and 4 exact numeric). Usage: wolframscript -file hv4_checks.wls [ppz4]   (ppz4: also the 379-graph genus-4 check) *)
SetDirectory[DirectoryName[$InputFileName]];
Get["hv_lib3.m"]; Get["hv_lib4.m"];
log[m_] := Module[{s = OpenAppend["hv4_checks.log"]}, WriteString[s, DateString[{"Year", "-", "Month", "-", "Day", " ", "Hour", ":", "Minute", ":", "Second"}], "  ", m, "\n"]; Close[s]; Print[m]];
tt = AbsoluteTime[];
{ok, rep} = hodgeChecks4[]; Do[log[r], {r, rep}];
log["hodgeChecks4 done in " <> ToString[Round[AbsoluteTime[] - tt]] <> " s"];
(* agreement with hv_lib3 (closed forms) on every decoration weight of genus 2 and 3 *)
tt = AbsoluteTime[];
Get["hv_orbits_g3.mx"];
decs3 = Union[Flatten[Table[Table[{gr["g"][[v]], Join[dec[[v, 1]], dec[[v, 2]], dec[[v, 3]]]}, {v, gr["V"]}], {gr, preparedH3}, {dec, gr["decs"]}], 2]];
agree = And @@ Table[Together[hodgeNormU4[d[[1]], d[[2]]] - hodgeNormU[d[[1]], d[[2]]]] === 0, {d, decs3}];
log["hodgeNormU4 = hodgeNormU (hv_lib3) on the " <> ToString[Length[decs3]] <> " distinct vertex correlators of the genus-3 decorations: " <> ToString[agree] <> " (" <> ToString[Round[AbsoluteTime[] - tt]] <> " s)"];
ok = ok && agree;
(* genus-4 vertex correlators: timing and a sample *)
tt = AbsoluteTime[];
samp4 = {hodgeNormU4[4, {}], hodgeNormU4[4, {2, 2, 2, 3}], hodgeNormU4[4, {2, 2, 2, 2, 2, 2, 2, 2, 2}], hodgeNormU4[4, {0, 0, 4, 5}]};
log["genus-4 correlators: hodgeNormU4[4, {}] = " <> ToString[InputForm[Factor[samp4[[1]]]]] <> "; [4,{2,2,2,3}] = " <> ToString[InputForm[Factor[samp4[[2]]]]] <> "; [4,{2^9}] = " <> ToString[InputForm[samp4[[3]]]] <> "; [4,{0,0,4,5}] = " <> ToString[InputForm[Factor[samp4[[4]]]]] <> " (" <> ToString[Round[AbsoluteTime[] - tt, 0.1]] <> " s)"];
wB = {e1s + e2s, -e1s, -e2s};
log["normalised genus-4 Hodge class for wB: " <> ToString[InputForm[Table[{c[[1]], Factor[c[[2]] /. usubOf[wB]]}, {c, hodgeClass4[4]}]]]];
(* PPZ scalar checks *)
ppz2 = ppzScalar[2, {w1, w2, w3}]; log["PPZ genus 2 symbolic = normalised class: " <> ToString[Together[ppz2 - (hodgeNormU4[2, {}] /. usubOf[{w1, w2, w3}])] === 0]];
SeedRandom[11];
Do[wn = Table[RandomInteger[{-9, 9}]/RandomInteger[{1, 5}], {3}]; If[MemberQ[wn, 0], Continue[]];
  c = ppzScalar[3, wn] == (hodgeNormU4[3, {}] /. usubOf[wn]); log["PPZ genus 3 at weights " <> ToString[wn] <> ": " <> ToString[c]]; ok = ok && c, {2}];
If[MemberQ[Rest[$ScriptCommandLine], "ppz4"],
  tt = AbsoluteTime[];
  graphs4 = prepareGraph /@ graphsWithLegs[4, 0];
  log["genus-4 point graphs: " <> ToString[Length[graphs4]] <> " graphs, " <> ToString[Total[Length[#["decs"]] & /@ graphs4]] <> " labelled decorations (" <> ToString[Round[AbsoluteTime[] - tt]] <> " s)"];
  Do[wn = Table[RandomInteger[{-9, 9}]/RandomInteger[{1, 5}], {3}]; If[MemberQ[wn, 0], Continue[]];
   tt = AbsoluteTime[]; c = ppzScalar[4, wn, graphs4] == (hodgeNormU4[4, {}] /. usubOf[wn]);
   log["PPZ genus 4 (379 graphs) at weights " <> ToString[wn] <> ": " <> ToString[c] <> " (" <> ToString[Round[AbsoluteTime[] - tt]] <> " s)"]; ok = ok && c, {4}]];
log[If[ok, "HV4 CHECKS OK", "HV4 CHECKS FAILED"]];
