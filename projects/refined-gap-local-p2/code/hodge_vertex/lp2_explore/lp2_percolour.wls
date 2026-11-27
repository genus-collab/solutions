(* lp2_percolour.wls -- CHECK 3 for K_P2 x C^2 (fibrewise torus T_F, mode B = Hodge vertices with the weights {E, -eps1, -eps2},
   E = eps1 + eps2, frame B: R' = R_Z h_wB^(-1), conifold polarisation c1 = C1cf of lp2_lib.m section 6).

   WHAT IT CHECKS
   (A) Factor table.  For each colour i = 0, 1, 2 (colour i = the idempotent evaluated at t -> w^i t, w = exp(2 pi i/3), on the
       REAL branch t = -(1/3)(1 - Delta)^(1/3), Delta = 1 + 27 y, so colour 0 is the singular idempotent at the conifold) the
       conifold valuation (order in t_c) of every per-colour factor of the graph sum: sqrtDelta_i, the rows P0, P1, P2 of R',
       the dilaton leaves Dil_i[n], the leg factors Leaf_i[n] and the edge factors E_ij[m, n] for all ordered colour pairs.
       Question: are the factors of the two regular idempotents (colours 1, 2) analytic at the conifold?
   (B) Per-colouring polar parts.  For every stable graph of genus gg and every colouring c of its vertices by {0, 1, 2}:
       the sum over ALL decorations of the graph, colouring fixed, expanded at t_c = 0; its polar coefficients t_c^(-P)..t_c^(-1).
       Also the sums over colourings with the same singular set S = c^(-1)(0) (i.e. summed over the regular colours only).
       Question: are graphs with an edge regular colouring by colouring (as for the resolved conifold, Theorem B/R), or only
       after the colour sum?
   (C) Self-check: the sum over all 3^V colourings of graph k equals, as a series in t_c through t_c^1, the stored colour-summed
       graph term (hodge_vertex_genus2_results.m ["B"]["graphTerms"], resp. hodge_vertex_genus3_results.m ["graphTerms"]) with
       XX -> XXcf, y -> ycf.  Also XX = y c1'/c1 for c1 = C1cf reproduces XXcf (XX(0) = -2/9).

   Arithmetic: exact (Q(sqrt(-3)) coefficients; eps2 = 1, eps1 = x with x in a list of rationals for which
   sq = sqrt(x/(1+x)^3) is rational: x = 1/(n^2-1), default {1/3, 1/8}; pass "x" for symbolic eps1).
   Usage:  wolframscript -file lp2_percolour.wls [gg = 2] [xlist = {1/3, 1/8} | x] [graphs = All | {k1, k2, ...}]
   Cost: genus 2 seconds to a minute per x; genus 3 (all 42 graphs, 3^V colourings, 1383 labelled decorations) estimated a few
   minutes to tens of minutes per x (untested).
   Output: lp2_percolour_g<gg>_out.txt (and .m with the data) in this folder (lp2_explore).

   EXPECTED OUTPUT
   - "self-check XX: ... True" and, for every graph, "colour-sum = stored term through t_c^1: True".  If False, the per-colour
     machinery (branch, c1 normalisation, colour labelling) is wrong and nothing below is meaningful.
   - Vertex graph (graph 1): colouring {0} carries the whole gap (genus 2: -(7+10x^2+7x^4)/(1920x) at t_c^-2, 0 at t_c^-1),
     colourings {1}, {2} are 0.
   - Factor table: if the regular idempotents are analytic, every colour-1/2 factor (and E_ij with i, j in {1, 2}) has valuation >= 0.
     Colour 0: sqrtDelta_0 is expected to have a simple pole (Delta_0 ~ t_c^-2); P1, P2 rows and edges touching colour 0 may have poles.
   - (B): the answer to CHECK 3 (unknown before running): "REGULAR colouring by colouring" or a list of singular colourings. *)

args = Rest[$ScriptCommandLine];
gg = If[Length[args] >= 1, ToExpression[args[[1]]], 2];
xlist = If[Length[args] >= 2, If[args[[2]] === "x", {"x"}, ToExpression[args[[2]]]], {1/3, 1/8}];
graphSel = If[Length[args] >= 3, ToExpression[args[[3]]], All];
here = DirectoryName[$InputFileName]; base = ParentDirectory[here];
SetDirectory[base];
Get["hv_lib3.m"];          (* loads lp2_lib.m, graphs.m, orbits.m; Hodge vertices for g <= 3 *)
outF = FileNameJoin[{here, "lp2_percolour_g" <> ToString[gg] <> "_out.txt"}];
If[FileExistsQ[outF], DeleteFile[outF]];
emit[m__] := Module[{s = OpenAppend[outF], str = StringJoin[ToString /@ {m}]}, WriteString[s, str, "\n"]; Close[s]; Print[str]];
fmt[e_] := ToString[InputForm[e]];

wB = {e1s + e2s, -e1s, -e2s};
KMAX = 3 gg - 3;
loadRmatrix[KMAX];
rows = rowsTimes[hCoeffs[wB, -1, KMAX], KMAX];      (* frame B rows (as hv_genus3.wls, stage sample) *)
ORD = 4 gg + 4;                                       (* absolute t_c order of every factor series *)
NCT = ORD + 4 gg + 8;
tcSetup[NCT];
Del0 = Normal[invcf];                                 (* Delta = 1 + 27 y as a polynomial in tcv *)
C1D = (Del - 1) D[tc, Del]; C1pD = 27 D[C1D, Del];    (* c1 = C1cf = y dt_c/dy, c1' = dc1/dy, as functions of Delta *)
treal = -(1/3) (1 - Del)^(1/3);
wwE = (-1 + I Sqrt[3])/2;
xxChk = Normal[Series[((Del - 1)/27) C1pD/C1D /. Del -> Del0, {tcv, 0, 4}]] - Normal[Series[XXcf, {tcv, 0, 4}]];
emit["genus ", gg, ", KMAX = ", KMAX, ", factor order ORD = ", ORD, ", t_c series order ", NCT, "; ", DateString[]];
emit["self-check XX: y c1'/c1 with c1 = C1cf equals XXcf through t_c^4: ", Expand[xxChk] === 0, "; XX(0) = ", Normal[XXcf] /. tcv -> 0,
  "; c1(0) = ", C1D /. Del -> 0];

csN[e_] := Expand[ComplexExpand[e]];
csS[e_] := Module[{c = ComplexExpand[e]}, Together[ComplexExpand[Re[c]]] + I Together[ComplexExpand[Im[c]]]];
cleanS[s_SeriesData, cs_] := Module[{l = cs /@ s[[3]]}, SeriesData[s[[1]], s[[2]], l, s[[4]], s[[5]], s[[6]]]];
cleanS[s_, cs_] := cs[s];
valS[s_SeriesData] := Module[{p = Position[s[[3]], c_ /; c =!= 0, {1}, 1, Heads -> False]}, If[p === {}, Infinity, (s[[4]] + p[[1, 1]] - 1)/s[[6]]]];
valS[0] = Infinity; valS[s_] := 0;
precS[s_SeriesData] := s[[5]]/s[[6]]; precS[_] := Infinity;
coefS[s_SeriesData, n_] := SeriesCoefficient[s, n]; coefS[s_, n_] := If[n == 0, s, 0];

results = <||>;
Do[Module[{xv = xr, sqv, cs, etai, fser, rowS, sqDS, DilS, LfS, EdS, preS, sqDp, prepared, usubN, graphs, stored,
    factorTable, perCol, tt0 = AbsoluteTime[]},
  If[xv === "x", xv = x; sqv = sq; cs = csS, sqv = Sqrt[xv/(1 + xv)^3]; cs = csN];
  emit[""]; emit["==================== eps1 = ", fmt[xv], ", eps2 = 1, sq = ", fmt[sqv], " ===================="];
  etai = etainvF /. {e1s -> xv, e2s -> 1};
  fser[expr_, i_] := cleanS[Series[expr /. {e1s -> xv, e2s -> 1, sq -> sqv} /. {C1i -> 1/C1D, C1p -> C1pD} /. t -> wwE^i treal /. Del -> Del0, {tcv, 0, ORD}], cs];
  Do[rowS[i, a, k] = fser[rows[[a + 1, k + 1]], i], {i, 0, 2}, {a, 0, 2}, {k, 0, KMAX}];
  Do[sqDS[i] = fser[sqDsym, i], {i, 0, 2}];
  Do[DilS[i, n] = (-1)^n rowS[i, 0, n - 1], {i, 0, 2}, {n, 1, KMAX + 1}];
  Do[LfS[i, n] = (-1)^n rowS[i, 1, n], {i, 0, 2}, {n, 0, KMAX}];
  Do[preS[i, j, m, n] = cleanS[-(-1)^(m + n) Sum[rowS[i, a, m] etai[[a + 1, b + 1]] rowS[j, b, n], {a, 0, 2}, {b, 0, 2}], cs], {i, 0, 2}, {j, 0, 2}, {m, 0, KMAX}, {n, 0, KMAX}];
  Do[EdS[i, j, m, n] = cleanS[Sum[(-1)^k preS[i, j, m + k + 1, n - k], {k, 0, n}], cs], {i, 0, 2}, {j, 0, 2}, {m, 0, KMAX - 1}, {n, 0, KMAX - 1 - m}];
  sqDp[i_, p_] := sqDp[i, p] = cleanS[sqDS[i]^p, cs];
  emit["factor series built in ", Round[AbsoluteTime[] - tt0], " s"];

  (* ---------- (A) factor table ---------- *)
  emit["(A) conifold valuations (order in t_c; Infinity = identically 0 to the order computed) of the per-colour factors"];
  factorTable = Table[Module[{vr, vd, vl},
     vr = Table[valS[rowS[i, a, k]], {a, 0, 2}, {k, 0, KMAX}];
     vd = Table[valS[DilS[i, n]], {n, 1, KMAX + 1}];
     vl = Table[valS[LfS[i, n]], {n, 0, KMAX}];
     emit["  colour ", i, ": sqrtDelta ", valS[sqDS[i]], "; rows P0 ", vr[[1]], ", P1 ", vr[[2]], ", P2 ", vr[[3]], "; Dil[1..", KMAX + 1, "] ", vd, "; Leaf[0..", KMAX, "] ", vl];
     <|"colour" -> i, "sqrtDelta" -> valS[sqDS[i]], "rows" -> vr, "Dil" -> vd, "Leaf" -> vl|>], {i, 0, 2}];
  Do[Module[{tab = Table[valS[EdS[i, j, m, n]], {m, 0, KMAX - 1}, {n, 0, KMAX - 1 - m}]},
     emit["  edge E_", i, j, "[m, n] (rows m = 0..", KMAX - 1, ", n = 0..KMAX-1-m): ", tab, "   min ", Min[Flatten[tab]]]], {i, 0, 2}, {j, 0, 2}];
  Module[{regmin = Min[Flatten[{Table[{valS[sqDS[i]], Table[valS[rowS[i, a, k]], {a, 0, 2}, {k, 0, KMAX}]}, {i, 1, 2}],
        Table[valS[EdS[i, j, m, n]], {i, 1, 2}, {j, 1, 2}, {m, 0, KMAX - 1}, {n, 0, KMAX - 1 - m}]}]]},
   emit["  regular idempotents (colours 1, 2): minimal valuation of sqrtDelta, P0, P1, P2 (hence Dil, Leaf), E_ij (i, j in {1,2}) = ", regmin,
     If[regmin >= 0, "  => ANALYTIC at the conifold (to z^" <> ToString[KMAX] <> ")", "  => NOT analytic"],
     "; sqrtDelta_1, sqrtDelta_2 at t_c = 0: ", fmt[{coefS[sqDS[1], 0], coefS[sqDS[2], 0]}]]];

  (* ---------- (B) per-colouring polar parts ---------- *)
  prepared = prepareGraphH /@ graphsWithLegs[gg, 0];
  usubN = usubOf[wB] /. {e1s -> xv, e2s -> 1};
  graphs = If[graphSel === All, Range[Length[prepared]], graphSel];
  stored = If[gg == 2, Get[FileNameJoin[{base, "hodge_vertex_genus2_results.m"}]]["B"]["graphTerms"],
    Get[FileNameJoin[{base, "hodge_vertex_genus3_results.m"}]]["graphTerms"]];
  emit["(B) per-colouring polar parts {t_c^-P .. t_c^-1}, sum over all decorations of the stable graph with the colouring fixed"];
  perCol = Table[Module[{gr = numGraph[prepared[[k]], usubN], V, cols, vals, tot, sto, chk, polar, Pmin, bySet, out},
     V = gr["V"]; cols = Tuples[{0, 1, 2}, V];
     vals = Association[Table[col -> Module[{ed = gr["edges"], H = gr["H"], hsAt = gr["hsAt"], tot0 = 0, hp},
          Do[hp = ConstantArray[0, Length[H]];
           Do[Do[hp[[hsAt[[v, r]]]] = dec[[v, 1, r]], {r, Length[hsAt[[v]]]}], {v, V}];
           tot0 += Product[dec[[v, 4]], {v, V}] *
             Product[sqDp[col[[v]], 2 gr["g"][[v]] - 2 + Length[hsAt[[v]]] + gr["m"][[v]] + Length[dec[[v, 3]]]] Product[DilS[col[[v]], d], {d, dec[[v, 3]]}], {v, V}] *
             Product[EdS[col[[ed[[e, 1]]]], col[[ed[[e, 2]]]], hp[[2 e - 1]], hp[[2 e]]], {e, Length[ed]}], {dec, gr["decs"]}];
          cleanS[tot0/gr["aut"], cs]], {col, cols}]];
     tot = cleanS[Total[Values[vals]], cs];
     sto = Series[stored[[k]] /. {e1s -> xv, e2s -> 1} /. XX -> Normal[XXcf] /. y -> ycf, {tcv, 0, 1}];
     chk = (And @@ Table[cs[coefS[tot, n] - SeriesCoefficient[sto, n]] === 0, {n, -(2 gg + 2), 1}]) && precS[tot] >= 2;
     Pmin = Min[0, Min[valS /@ Values[vals]]];
     polar = Association[Table[col -> Table[If[precS[vals[col]] > -n, Factor[cs[coefS[vals[col], -n]]], "undetermined"], {n, -Pmin, 1, -1}], {col, cols}]];
     bySet = Association[Table[S -> Module[{sel = Select[cols, Flatten[Position[#, 0]] == S &]},
          Table[Factor[cs[Total[Table[If[precS[vals[c]] > -n, coefS[vals[c], -n], Indeterminate], {c, sel}]]]], {n, -Pmin, 1, -1}]], {S, Subsets[Range[V]]}]];
     emit["  graph ", k, ": genera ", gr["g"], ", edges ", gr["edges"], ", Aut ", gr["aut"], ", decorations ", Length[gr["decs"]],
       "; colour-sum = stored term through t_c^1: ", chk, " (precision t_c^", precS[tot], ")"];
     Do[If[Union[polar[col]] =!= {0} && polar[col] =!= {}, emit["      colouring ", col, " (prec t_c^", precS[vals[col]], "): ", fmt[polar[col]]]], {col, cols}];
     out = Select[cols, polar[#] =!= {} && Union[polar[#]] =!= {0} &];
     emit["      singular colourings: ", Length[out], " of ", Length[cols], If[out === {}, "  => REGULAR colouring by colouring", ""],
       If[Length[gr["edges"]] > 0 && out =!= {}, "  [graph with an edge NOT regular per colouring]", ""]];
     Do[If[bySet[S] =!= {} && Union[bySet[S]] =!= {0}, emit["      singular set S = ", S, " summed over the regular colours: ", fmt[bySet[S]]]], {S, Subsets[Range[V]]}];
     emit["      colour sum: ", fmt[Table[Factor[cs[coefS[tot, -n]]], {n, -Pmin, 1, -1}]]];
     <|"graph" -> k, "genera" -> gr["g"], "edges" -> gr["edges"], "check" -> chk, "polar" -> polar, "bySingularSet" -> bySet,
       "singularColourings" -> out|>], {k, graphs}];
  emit["SUMMARY eps1 = ", fmt[xv], ": self-checks all True: ", And @@ (#["check"] & /@ perCol),
    "; graphs with an edge that are singular for some colouring: ", #["graph"] & /@ Select[perCol, Length[#["edges"]] > 0 && #["singularColourings"] =!= {} &],
    "; graphs with an edge singular for some singular set S (after summing the regular colours): ",
    #["graph"] & /@ Select[perCol, Function[pc, Length[pc["edges"]] > 0 && Select[Values[pc["bySingularSet"]], # =!= {} && Union[#] =!= {0} &] =!= {}]],
    "; time ", Round[AbsoluteTime[] - tt0], " s"];
  results[xv] = <|"factorTable" -> factorTable, "perColouring" -> perCol|>], {xr, xlist}];
Put[results, FileNameJoin[{here, "lp2_percolour_g" <> ToString[gg] <> "_out.m"}]];
emit["written ", outF];
