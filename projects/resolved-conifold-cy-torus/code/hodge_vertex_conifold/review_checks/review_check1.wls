(* review_check1.wls -- independent review checks (read-only on the project data; writes only to scratch/ in this folder) *)
SetDirectory[DirectoryName[$InputFileName]];
scr = FileNameJoin[{Directory[], "scratch"}]; If[!DirectoryQ[scr], CreateDirectory[scr]];
tt = AbsoluteTime[];
Get["hvc_lib.m"];
Print["=== 1. notebook (14:25) inlined code cells versus the current files ==="];
nb = Get["resolvedconifold_hodge.nb"];
inputs = Cases[nb, Cell[s_String, "Input"] :> s, Infinity];
strip[f_] := StringTrim[StringReplace[Import[f, "Text"], StartOfLine ~~ "Get[\"" ~~ Except["\n"] ... ~~ "\n" -> ""]];
Export[FileNameJoin[{scr, "review_nb_cell4.txt"}], inputs[[4]], "Text"];
Export[FileNameJoin[{scr, "review_cur_hvclib.txt"}], strip["hvc_lib.m"], "Text"];
Print["cell 4 (hvc_lib.m) identical to current: ", inputs[[4]] === strip["hvc_lib.m"], "; lengths ", StringLength[inputs[[4]]], " vs ", StringLength[strip["hvc_lib.m"]]];
Print["cell 1 = graphs.m: ", inputs[[1]] === StringTrim[Import["graphs.m", "Text"]], "; cell 2 = orbits.m: ", inputs[[2]] === StringTrim[Import["orbits.m", "Text"]], "; cell 3 = rc_lib.m: ", inputs[[3]] === strip["rc_lib.m"]];

Print["=== 2. QRR constants (closed form) at alpha0 = 0, alpha_i = eps_i ==="];
setupWick[3];
Print[ffcClosed[{a0 -> 0, a1 -> x, a2 -> 1, e1 -> x, e2 -> 1}]];

Print["=== 3. polarOfRational ordering: 1/(1-Q)^2 at g = 2 -> ", polarOfRational[1/(1 - Q)^2, 2], " (expect {1, -1}); 1/(1-Q)^4 at g = 3 -> ", polarOfRational[1/(1 - Q)^4, 3]];

(* check 4 (gap normalisation local P2 / conifold) removed on 2026-11-27: it needed lp2_tc_section.m of the local P2 project; the
   local P2 comparison now lives in that project (hodge_vertex/compare_conifold/). Check 1 no longer tests a cell 5 for the same reason. *)

Print["=== 5. stored genus-2 factors: Bm rows = point rows (Cauchy) hCoeffs[wB, -1, 3] ? ==="];
Get["hvc_factors_g2_point.mx"];
Do[pt[a, i, m] = ratCN[a, i][m], {a, 0, 1}, {i, 1, 2}, {m, 0, ZMAX}];
Clear[ratCN]; Get["hvc_factors_g2_Bm.mx"];
hcm = hCoeffs[modeWeights["B"] /. {e1s -> x, e2s -> 1}, -1, 3];
hcp = hCoeffs[modeWeights["B"] /. {e1s -> x, e2s -> 1}, 1, 3];
twist[hc_] := Table[Module[{tw = cadd @@ Table[cscale[hc[[m - j + 1]], pt[a, i, j]], {j, 0, m}], st = ratCN[a, i][m]}, {Together[tw[[1]] - st[[1]]], Together[tw[[2]] - st[[2]]]}], {a, 0, 1}, {i, 1, 2}, {m, 0, 3}];
Print["with h^-1 (s = -1): all rows agree: ", Union[Flatten[twist[hcm]]] === {0}, ";  with h^+1 instead: ", Union[Flatten[twist[hcp]]] === {0}];
Print["point row 0, root 1, z^1 coefficient: ", ratCN[0, 1][1] /. C[a_, b_] :> {Factor[a], Factor[b]} // InputForm];

Print["=== 6. scalar R-action in the rc_lib graph conventions (Dil = (-1)^n R_{n-1}, Ed = Sum_k (-1)^k pre, pre = -(-1)^(p+q) R_p R_q): which sign s of h_w^s on point vertices gives the Hodge class? ==="];
coniScalar[gg_, w_, s_] := Module[{KMAX = 3 gg - 3, A, pre, Ed, Dil, grs, tot = 0, hp, vfac, efac},
   A = hCoeffs[w, s, KMAX];
   pre[p_, q_] := -(-1)^(p + q) A[[p + 1]] A[[q + 1]];
   Ed[m_, n_] := Ed[m, n] = Together[Sum[(-1)^k pre[m + k + 1, n - k], {k, 0, n}]];
   Dil[n_] := (-1)^n A[[n]];
   grs = prepareGraph /@ graphsWithLegs[gg, 0];
   Do[Do[hp = ConstantArray[0, Length[gr["H"]]];
     Do[Do[hp[[gr["hsAt"][[v, r]]]] = dec[[v, 1, r]], {r, Length[gr["hsAt"][[v]]]}], {v, gr["V"]}];
     vfac = Product[dec[[v, 4]] Product[Dil[d], {d, dec[[v, 3]]}], {v, gr["V"]}];
     efac = Product[Ed[hp[[2 e - 1]], hp[[2 e]]], {e, Length[gr["edges"]]}];
     tot += vfac efac/gr["aut"], {dec, gr["decs"]}], {gr, grs}];
   Together[tot]];
hod2 = hodgeNormU4[2, {}] /. usubOf[{w1, w2, w3}];
Print["genus 2, general weights: s = +1 gives the Hodge class: ", Together[coniScalar[2, {w1, w2, w3}, 1] - hod2] === 0, ";  s = -1 gives it: ", Together[coniScalar[2, {w1, w2, w3}, -1] - hod2] === 0, ";  value with s = +1: ", Factor[coniScalar[2, {w1, w2, w3}, 1]] // InputForm];
wn = {8/5, 3/2, -7/2}; hod3 = hodgeNormU4[3, {}] /. usubOf[wn];
Print["genus 3 at ", wn, ": s = +1: ", coniScalar[3, wn, 1] == hod3, ";  s = -1: ", coniScalar[3, wn, -1] == hod3];
Print["(PPZ-convention check of the pipeline for comparison, ppzScalar with h^-1, genus 2): ", Together[ppzScalar[2, {w1, w2, w3}] - hod2] === 0];

Print["=== 7. stored genus-2 Bm results: internal consistency ==="];
res2 = Get["hvc_results_g2_Bm.m"]; en = res2["entries"]; st = res2["stable"];
Print["entries: ", Length[en], "; stored polar = polarOfRational[stored term] for all: ", And @@ Table[Module[{rq = Together[e[[4]]], p}, p = If[rq === 0, ConstantArray[0, 2], polarOfRational[rq, 2]]; Union[Together /@ (p - e[[5]])] === {0}], {e, en}]];
Print["stable terms = sums of entry terms: ", And @@ Table[Together[st[[k]]["term"] - Total[Select[en, #[[1, 1]] == k &][[All, 4]]]] === 0, {k, 7}]];
Print["total = F20: ", Together[Total[st[[All, "term"]]] - FrefHVC[2, x]] === 0, "; graph-1 polar = gap: ", Union[Together /@ (st[[1]]["polar"] - gapHVC[2, x])] === {0}];
Print["regular (polar list zero) <=> no (1-Q) in the denominator, per entry: ", And @@ Table[(Union[Together /@ e[[5]]] === {0}) === (PolynomialRemainder[Denominator[Together[e[[4]]]], 1 - Q, Q] =!= 0), {e, en}]];

Print["=== 8. stored genus-3 Bm results: polar = polarOfRational[term] on graphs 1, 2, 5, 6 (time-limited) ==="];
res3 = Get["hvc_results_g3_Bm.m"]; en3 = Select[res3["entries"], MemberQ[{1, 2, 5, 6}, #[[1, 1]]] &];
t0 = AbsoluteTime[]; nchk = 0; nbad = 0;
TimeConstrained[Do[Module[{rq = Together[e[[4]]], p}, p = If[rq === 0, ConstantArray[0, 4], polarOfRational[rq, 3]]; nchk++; If[Union[Together /@ (p - e[[5]])] =!= {0}, nbad++; Print["MISMATCH ", e[[1]]]]], {e, en3}], 80];
Print[nchk, " of ", Length[en3], " checked in ", Round[AbsoluteTime[] - t0], " s; mismatches: ", nbad];
Print["regular <=> no (1-Q) in the denominator, all 693 entries: ", And @@ Table[(Union[Together /@ e[[5]]] === {0}) === (PolynomialRemainder[Denominator[Together[e[[4]]]], 1 - Q, Q] =!= 0), {e, res3["entries"]}]];
Print["done in ", Round[AbsoluteTime[] - tt], " s"];
