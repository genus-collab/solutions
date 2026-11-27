(* rc_polarCY_sample.m -- sampling worker on the CY torus alpha0 = eps1 + eps2 - alpha1 - alpha2 with alpha1, alpha2 free:
   at every integer weight point assigned to this worker, the polar part at Q = 1 (t = Log Q) of every decorated stable graph,
   from Laurent-series arithmetic (closed-form QRR normalisation, eps2 = 1).
   Usage:  wolframscript -file rc_polarCY_sample.m genus ZMAX worker nWorkers
   Output: rc_polarCY_g<genus>_w<w>of<n>.mx holding `recs` = list of {point {a1, a2, e1}, matrix (orbits x orders t^-(2g-2)..t^-1),
           totalOK, seconds} (re-saved after every point; a restarted job skips the points already present); log in rc_polarCY_sample.log.
   Needs: rc_lib.m, graphs.m, orbits.m, StableGraph-<genus>.bin, rc_orbits_g<genus>.mx. *)
SetDirectory[DirectoryName[$InputFileName]];
args = Rest[$ScriptCommandLine];
If[Length[args] < 4, Print["usage: wolframscript -file rc_polarCY_sample.m genus ZMAX worker nWorkers"]; Exit[1]];
{gg, zmax, w, nw} = ToExpression /@ args[[1 ;; 4]];
Get["rc_lib.m"];
tag = "g" <> ToString[gg] <> "_w" <> ToString[w] <> "of" <> ToString[nw];
log[m_] := Module[{s = OpenAppend["rc_polarCY_sample.log"]}, WriteString[s, DateString[{"Year", "-", "Month", "-", "Day", " ", "Hour", ":", "Minute", ":", "Second"}], "  [", tag, "] ", m, "\n"]; Close[s]; Print[m]];
(* integer points (a1, a2, e1), eps2 = 1, alpha0 = e1 + 1 - a1 - a2, all weights and weight sums generic *)
nCYFor[g_] := {40, 160, 460}[[g - 1]];
cyOK[{A1_, A2_, E1_}] := Module[{A0 = E1 + 1 - A1 - A2}, And @@ (# != 0 & /@ {A1, A2, E1, E1 + 1, A0, A0 + A1, A0 + A2, E1 + 1 - 2 A1, E1 + 1 - 2 A2, A1 - E1, A2 - 1, A0 + A2 - 1, A0 + A1 - E1})];
cyPoints[g_] := Module[{ps = {}}, SeedRandom[20260907];
   While[Length[ps] < nCYFor[g], Module[{p = RandomInteger[{-12, 12}, 3]}, If[cyOK[p] && ! MemberQ[ps, p], AppendTo[ps, p]]]]; ps];
cyWsub[p_] := {a0 -> p[[3]] + 1 - p[[1]] - p[[2]], a1 -> p[[1]], a2 -> p[[2]], e1 -> p[[3]], e2 -> 1};
log["start: genus " <> ToString[gg] <> ", ZMAX " <> ToString[zmax] <> ", " <> $Version];
setupWick[zmax]; loadOrbits[gg];
pts = cyPoints[gg];
(* one DumpSave file per finished point: rc_polarCY_<tag>_p<index>.mx holding `rec` = {point, matrix, totalOK, seconds} *)
mine = Range[w, Length[pts], nw];
pfile[idx_] := "rc_polarCY_" <> tag <> "_p" <> ToString[idx] <> ".mx";
todo = Select[mine, ! FileExistsQ[pfile[#]] &];
log[ToString[Length[mine]] <> " points assigned (of " <> ToString[Length[pts]] <> "), " <> ToString[Length[mine] - Length[todo]] <> " already done, " <> ToString[Length[todo]] <> " to do; " <> ToString[Length[OT["ids"]]] <> " decorated graphs"];
(* absolute series order: a product with E edge factors keeps its coefficients up to t^(M - 2E) (measured at genus 3, with the vertex
   factors expanded from exact rational functions); E <= 3g - 3, so M = 6g - 2 leaves a margin of four orders above t^-1 *)
Mord = 6 gg - 2;
Do[Module[{p = pts[[idx]], wsub, t0 = AbsoluteTime[], mat, tot, expected, ok, nonsym = 0},
   wsub = cyWsub[p];
   rmatrixAtClosed[wsub]; evalQseries[Mord, gg];
   mat = Table[Module[{v = graphValueC[OT["graphs"][id]]}, If[! cSymQ[v], nonsym++]; polarOfSeries[v[[1]], gg]], {id, OT["ids"]}];
   tot = Total[mat]; expected = cyPolar[gg] /. wsub;
   ok = (Together /@ (tot - expected)) === ConstantArray[0, 2 gg - 2] && nonsym == 0;
   rec = {p, mat, ok, Round[AbsoluteTime[] - t0, 0.1]}; DumpSave[pfile[idx], {rec}];
   log["point " <> ToString[p] <> " done in " <> ToString[Round[AbsoluteTime[] - t0, 0.1]] <> " s; sum of orbits = P_g Li polar part: " <> ToString[ok] <> If[nonsym > 0, "; WARNING " <> ToString[nonsym] <> " non-symmetric orbit values", ""]]], {idx, todo}];
log["job finished"];
