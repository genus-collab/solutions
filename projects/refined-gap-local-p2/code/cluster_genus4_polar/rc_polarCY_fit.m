(* rc_polarCY_fit.m -- fitting worker: exact polar coefficients on the CY torus (alpha1, alpha2, eps1, eps2 symbolic) for the
   decorated graphs of one chunk, from the sampled values.
   Ansatz per order index n = 1..2g-2 (n = 1 is the leading t^-(2g-2)): coefficient = N/(eps1^da eps2^da (eps1+eps2)^db) with
   da = 2g-2, db = 2(n-1), N homogeneous of degree d = 3(2g-2) + 2(n-1), symmetric under alpha1<->alpha2, eps1<->eps2, alpha-degree
   <= 2(n-1) (the structure found at genus 2 and 3); on failure the exponents are escalated (db + 2, then da + 1) and the
   alpha-degree bound raised, everything verified on >= 8 spare points.
   Usage:  wolframscript -file rc_polarCY_fit.m genus chunk nChunks
   Input:  rc_polarCY_chunk_g<genus>_c<chunk>of<n>.mx, rc_orbits_g<genus>.mx
   Output: rc_polarCY_fit_g<genus>_c<chunk>of<n>.m = list of {id, key, description, polar list (rational functions), ansatz info}
           and the failed ids; log in rc_polarCY_fit.log *)
SetDirectory[DirectoryName[$InputFileName]];
{gg, cc, nchunks} = ToExpression /@ $ScriptCommandLine[[2 ;; 4]];
Get["rc_lib.m"];
tag = "g" <> ToString[gg] <> "_c" <> ToString[cc] <> "of" <> ToString[nchunks];
log[m_] := Module[{s = OpenAppend["rc_polarCY_fit.log"]}, WriteString[s, DateString[{"Year", "-", "Month", "-", "Day", " ", "Hour", ":", "Minute", ":", "Second"}], "  [", tag, "] ", m, "\n"]; Close[s]; Print[m]];
loadOrbits[gg];
Get["rc_polarCY_chunk_" <> tag <> ".mx"];   (* chunkRecs, chunkIdx, ptsAll *)
log["start: " <> ToString[Length[chunkIdx]] <> " decorated graphs, " <> ToString[Length[ptsAll]] <> " points"];
result = <||>; info = <||>; failedAll = {};
Do[Module[{vals, zero, pending, d, da, db, DaMax, tries = 0, fit, tt = AbsoluteTime[]},
   vals = Association[Table[chunkIdx[[pos]] -> Table[chunkRecs[p][[pos, n]], {p, ptsAll}], {pos, Length[chunkIdx]}]];
   zero = Select[Keys[vals], Union[vals[#]] === {0} &];
   Do[result[{i, n}] = 0; info[{i, n}] = "zero", {i, zero}];
   pending = KeyDrop[vals, zero];
   da = 2 gg - 2; db = 2 (n - 1); d = 3 (2 gg - 2) + 2 (n - 1); DaMax = 2 (n - 1) + 2;
   While[Length[pending] > 0 && tries < 6,
    fit = fitGroup[pending, ptsAll, d, da, db, DaMax];
    Do[result[{i, n}] = fit[[1]][i]; info[{i, n}] = {d, da, db, fit[[3]]}, {i, Keys[fit[[1]]]}];
    pending = KeyTake[pending, fit[[2]]];
    If[Length[pending] > 0, tries++; If[OddQ[tries], db += 2; d += 2, da += 1; db -= 2; d -= 0]; DaMax += 2;
      log["order " <> ToString[n] <> ": " <> ToString[Length[pending]] <> " coefficients not verified, escalating to d = " <> ToString[d] <> ", da = " <> ToString[da] <> ", db = " <> ToString[db] <> ", DaMax = " <> ToString[DaMax]]]];
   If[Length[pending] > 0, AppendTo[failedAll, {n, Keys[pending]}]];
   log["order index " <> ToString[n] <> " (t^-" <> ToString[2 gg - 1 - n] <> "): " <> ToString[Length[chunkIdx] - Length[zero] - Length[pending]] <> " fitted, " <> ToString[Length[zero]] <> " zero, " <> ToString[Length[pending]] <> " FAILED; alpha-degrees " <> ToString[Tally[Table[If[ListQ[info[{i, n}]], info[{i, n}][[4]], -1], {i, chunkIdx}]]] <> " (" <> ToString[Round[AbsoluteTime[] - tt, 0.1]] <> " s)"]], {n, 1, 2 gg - 2}];
entries = Table[Module[{id = OT["ids"][[i]]}, {id, OT["orbits"][[id[[1]], id[[2]]]], Normal[OT["descriptions"][id]], Table[Lookup[result, Key[{i, n}], $Failed], {n, 1, 2 gg - 2}], Table[Lookup[info, Key[{i, n}], "failed"], {n, 1, 2 gg - 2}]}], {i, chunkIdx}];
Put[{entries, failedAll}, "rc_polarCY_fit_" <> tag <> ".m"];
log["finished; failed: " <> ToString[failedAll]];
