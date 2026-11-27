(* lp2_lib.m -- library for the per-decorated-graph Givental-Teleman graph sum of K_P2 x C^2 on the fibrewise torus T_F
   (untwisted R-matrix), genus <= 4.  Load with Get after SetDirectory to the data directory (graphs.m, orbits.m,
   cPP_z<K>.m, StableGraph-<g>.bin must be there).
   Sections 1 (core routines), 3 (rational reconstruction) and the colouring sum by variable elimination are verbatim
   from genRmat_g4_lib.m (cluster_genus4/, itself verbatim from localp2CY5_genus3.nb); the twist is removed.
   Conventions: y = t^3; eps1, eps2 are the symbols Epsilon[1], Epsilon[2] (typed as backslash-bracket-Epsilon);
   C1i = 1/c1, C1p = dc1/dy; XX = y C1p/c1; PREC digits. ASCII only. *)

PREC = 350;
ww = Exp[2 Pi I/3];
Get["graphs.m"]; Get["orbits.m"];

(* ============================== 1. core routines ============================== *)
alpha = (\[Epsilon][1] + \[Epsilon][2])/3;
eeF[1] = 3 alpha; eeF[2] = 3 alpha^2; eeF[3] = alpha^3;
LLt0 = -alpha/(1 + 3 t);                         (* LLL[0][y] on T_F with y = t^3 *)
sqrule = Power[(\[Epsilon][1] \[Epsilon][2])/(\[Epsilon][1] + \[Epsilon][2])^3, r_] :> sq^(2 r);
etaF = -{{-(1/(3 eeF[3] \[Epsilon][1] \[Epsilon][2])), 0, 0}, {0, 0, 1/(3 \[Epsilon][1] \[Epsilon][2])}, {0, 1/(3 \[Epsilon][1] \[Epsilon][2]), -(eeF[1]/(3 \[Epsilon][1] \[Epsilon][2]))}};
etainvF = Simplify[Inverse[etaF]];
Dy[expr_] := D[expr, t]/(3 t^2) - C1p C1i^2 D[expr, C1i] + c1pp D[expr, C1p];

buildRmatrix[KMAX_] := Module[{P0, P1, P2, c2, d2, k},
   cPPtk[k_] := cPPtk[k] = Together[cPP[0, 0][k] /. LLL[0][y] -> LLt0 /. y -> t^3 /. sqrule];
   P0 = Sum[cPPtk[k] z^k, {k, 0, KMAX}];
   P1 = Normal[Series[(LLt0 P0 + z t^3 Dy[P0]) C1i, {z, 0, KMAX}]];
   c2 = C1i^2/(1 + 27 t^3); d2 = (\[Epsilon][1] + \[Epsilon][2]) (C1i^2 - 1)/(2 (1 + 27 t^3));
   P2 = Normal[Series[(LLt0 P1 + z t^3 Dy[P1] - d2 P1)/c2, {z, 0, KMAX}]];
   If[! FreeQ[{P1, P2}, c1pp], Print["WARNING: c1'' appeared"]];
   PPsym[0] = Table[Together[Coefficient[P0, z, k]], {k, 0, KMAX}];
   PPsym[1] = Table[Together[Coefficient[P1, z, k]], {k, 0, KMAX}];
   PPsym[2] = Table[Together[Coefficient[P2, z, k]], {k, 0, KMAX}];
   sqDsym = Sqrt[27] sq (3 eeF[3] + eeF[2] LLt0);   (* sqrtDelta[0] = Sqrt[eps1 eps2/ee3] (3 ee3 + ee2 L) *)
   KMAXcur = KMAX;
   Print["R-matrix rows built to z^", KMAX, "; LeafCounts row0/1/2: ", LeafCount /@ {PPsym[0], PPsym[1], PPsym[2]}];
   ];

evalPoint[t0_, x0_, KMAX_] := Module[{sub, Pnum, i, a, k, pre, Dil, Lf, Ed, sqD, etainv, sqv, tt},
   sqv = N[Sqrt[x0/(1 + x0)^3], PREC];
   etainv = N[etainvF /. {\[Epsilon][1] -> x0, \[Epsilon][2] -> 1}, PREC];
   Do[tt = N[ww^i t0, PREC];
    sub = {t -> tt, \[Epsilon][1] -> N[x0, PREC], \[Epsilon][2] -> N[1, PREC], sq -> sqv};
    Do[Pnum[a, i] = Expand[PPsym[a] /. sub], {a, 0, 2}];
    sqD[i] = Expand[sqDsym /. sub];
    Do[Dil[i][n] = (-1)^n Pnum[0, i][[n]], {n, 1, KMAX + 1}];
    Do[Lf[i][n] = (-1)^n Pnum[1, i][[n + 1]], {n, 0, KMAX}];
    , {i, 0, 2}];
   Do[
    pre = Table[KroneckerDelta[i, j] KroneckerDelta[m, 0] KroneckerDelta[n, 0] -
       Sum[(-1)^(m + n) Pnum[a, i][[m + 1]] etainv[[a + 1, b + 1]] Pnum[b, j][[n + 1]], {a, 0, 2}, {b, 0, 2}],
      {m, 0, KMAX}, {n, 0, KMAX}];
    Do[Ed[i, j][m, n] = Expand[Sum[(-1)^k pre[[m + k + 2, n - k + 1]], {k, 0, n}]],
     {m, 0, KMAX - 1}, {n, 0, KMAX - 1 - m}];
    , {i, 0, 2}, {j, 0, 2}];
   etainvN = etainv; <|"Dil" -> Dil, "Leaf" -> Lf, "Edge" -> Ed, "sqD" -> sqD, "Pnum" -> Pnum|>];

(* graph value by enumeration of the 3^V colourings (reference implementation) *)
graphValue[gr_, data_] := Module[{V = gr["V"], ed = gr["edges"], H = gr["H"], hsAt = gr["hsAt"], Dil = data["Dil"], Lf = data["Leaf"], Ed = data["Edge"], sqD = data["sqD"], total = 0, hp, vfac, efac},
   Do[Do[hp = ConstantArray[0, Length[H]];
     Do[Do[hp[[hsAt[[v, r]]]] = dec[[v, 1, r]], {r, Length[hsAt[[v]]]}], {v, V}];
     vfac = Product[dec[[v, 4]] sqD[c[[v]]]^(2 gr["g"][[v]] - 2 + Length[hsAt[[v]]] + gr["m"][[v]] + Length[dec[[v, 3]]])
        Product[Lf[c[[v]]][b], {b, dec[[v, 2]]}] Product[Dil[c[[v]]][d], {d, dec[[v, 3]]}], {v, V}];
     efac = Product[Ed[c[[ed[[e, 1]]]], c[[ed[[e, 2]]]]][hp[[2 e - 1]], hp[[2 e]]], {e, Length[ed]}];
     total += vfac efac, {dec, gr["decs"]}], {c, Tuples[{0, 1, 2}, V]}];
   Expand[total/gr["aut"]]];

(* colouring sum by variable elimination (verbatim from genRmat_g4_lib.m section 4) *)
contractColourings[V_, ed_, vf_, ef_] := Module[{factors, remaining = Range[V], v, cpos, cand, U, keep, posU, idx, keepidx, newassoc, key, val, sizes},
   factors = Join[
     Table[{{w}, Association[Table[{c} -> vf[[w, c + 1]], {c, 0, 2}]]}, {w, V}],
     Table[If[ed[[e, 1]] == ed[[e, 2]],
        {{ed[[e, 1]]}, Association[Table[{c} -> ef[[e, c + 1, c + 1]], {c, 0, 2}]]},
        {ed[[e]], Association[Flatten[Table[{c1, c2} -> ef[[e, c1 + 1, c2 + 1]], {c1, 0, 2}, {c2, 0, 2}]]]}], {e, Length[ed]}]];
   While[remaining =!= {},
    sizes = Table[Length[Union @@ Select[factors[[All, 1]], MemberQ[#, w] &]], {w, remaining}];
    v = remaining[[First[Ordering[sizes, 1]]]];
    remaining = DeleteCases[remaining, v];
    cpos = Flatten[Position[factors[[All, 1]], l_ /; MemberQ[l, v], {1}, Heads -> False]];
    cand = factors[[cpos]];
    U = Union @@ cand[[All, 1]]; keep = DeleteCases[U, v];
    posU = AssociationThread[U -> Range[Length[U]]];
    idx = Table[Lookup[posU, f[[1]]], {f, cand}];
    keepidx = Lookup[posU, keep];
    newassoc = Association[Table[k -> 0, {k, Tuples[{0, 1, 2}, Length[keep]]}]];
    Do[val = Expand[Times @@ Table[cand[[j, 2]][asg[[idx[[j]]]]], {j, Length[cand]}]];
     key = asg[[keepidx]];
     newassoc[key] += val, {asg, Tuples[{0, 1, 2}, Length[U]]}];
    factors = Append[Delete[factors, List /@ cpos], {keep, newassoc}]];
   Expand[Times @@ Table[f[[2]][{}], {f, factors}]]];
graphValueVE[gr_, data_] := Module[{V = gr["V"], ed = gr["edges"], H = gr["H"], hsAt = gr["hsAt"], Dil = data["Dil"], Lf = data["Leaf"], Ed = data["Edge"], sqD = data["sqD"], total = 0, hp, vf, ef},
   Do[
    hp = ConstantArray[0, Length[H]];
    Do[Do[hp[[hsAt[[v, r]]]] = dec[[v, 1, r]], {r, Length[hsAt[[v]]]}], {v, V}];
    vf = Table[Expand[dec[[v, 4]] sqD[c]^(2 gr["g"][[v]] - 2 + Length[hsAt[[v]]] + gr["m"][[v]] + Length[dec[[v, 3]]])
        Product[Lf[c][b], {b, dec[[v, 2]]}] Product[Dil[c][d], {d, dec[[v, 3]]}]], {v, V}, {c, 0, 2}];
    ef = Table[Ed[c1, c2][hp[[2 e - 1]], hp[[2 e]]], {e, Length[ed]}, {c1, 0, 2}, {c2, 0, 2}];
    total += contractColourings[V, ed, vf, ef], {dec, gr["decs"]}];
   Expand[total/gr["aut"]]];

(* ============================== 3. generic rational reconstruction (verbatim) ============================== *)
reconstruct[pts_] := Module[{n = Length[pts], us = pts[[All, 1]], vs = pts[[All, 2]], res = $Failed, T, dd, dn, m, A, rhs, sol, N0, D0, ok},
   If[Union[vs] === {0}, Return[0]];
   Do[
    Do[dn = T - dd; m = dn + dd + 1;
     If[m + 5 > n, Continue[]];
     A = Table[Join[Table[us[[s]]^j, {j, 0, dn}], Table[-vs[[s]] us[[s]]^j, {j, 0, dd - 1}]], {s, 1, m}];
     rhs = Table[vs[[s]] us[[s]]^dd, {s, 1, m}];
     sol = Quiet[Check[LinearSolve[A, rhs], $Failed]];
     If[sol === $Failed, Continue[]];
     N0 = Sum[sol[[j + 1]] u^j, {j, 0, dn}]; D0 = u^dd + Sum[sol[[dn + 2 + j]] u^j, {j, 0, dd - 1}];
     ok = And @@ Table[((N0 /. u -> us[[s]]) - vs[[s]] (D0 /. u -> us[[s]])) == 0, {s, m + 1, n}];
     If[ok, res = Together[N0/D0]; Return[]],
     {dd, 0, T}];
    If[res =!= $Failed, Return[]],
    {T, 0, n - 6}];
   res];

(* ============================== 4. R-matrix file, sample points, per-orbit sampling ============================== *)
loadRmatrix[KMAX_] := Module[{avail, file, lst},
   avail = Select[ToExpression[StringReplace[FileNameTake /@ FileNames["cPP_z*.m"], {"cPP_z" -> "", ".m" -> ""}]], IntegerQ[#] && # >= KMAX &];
   If[avail === {}, Print["ERROR: no cPP_z<K>.m with K >= ", KMAX, " in ", Directory[]]; Abort[]];
   file = "cPP_z" <> ToString[Min[avail]] <> ".m"; lst = Get[file];
   If[! ListQ[lst] || Length[lst] < KMAX + 1, Print["ERROR: ", file, " does not contain ", KMAX + 1, " coefficients"]; Abort[]];
   Print["R-matrix coefficients from ", file];
   Do[cPP[0, 0][k] = lst[[k + 1]], {k, 0, KMAX}];
   buildRmatrix[KMAX]];

kmaxFor[gg_] := 3 gg - 3;
tsamplesFor[gg_] := If[gg <= 3,
   {1/2, 1/3, 2/3, 1/4, 3/4, 1/5, 2/5, 3/5, 4/5, 1/6, 5/6, 1/7, 2/7, 3/7, 4/7, 5/7, 6/7, 1/8, 3/8, 5/8, 7/8, 1/9, 2/9, 1/10, 3/10, 7/10},
   {1/2, 1/3, 2/3, 1/4, 3/4, 1/5, 2/5, 3/5, 4/5, 1/6, 5/6, 1/7, 2/7, 3/7, 4/7, 5/7, 6/7, 1/8, 3/8, 5/8, 7/8, 1/9, 2/9, 1/10, 3/10, 7/10, 9/10, 1/11, 2/11, 3/11, 4/11, 5/11, 6/11, 7/11}];
(* the XX^k coefficient of a decorated graph is N(eps)/(eps1 eps2)^a with a <= 2g-2 and N homogeneous of degree 2g-2+2a:
   3(2g-2)+1 eps-ratios suffice for the largest a, plus 2 spare *)
xsamplesFor[gg_] := Table[Prime[r + 1], {r, 1, 3 (2 gg - 2) + 3}];
pointList[gg_] := Flatten[Table[{t0, x0}, {t0, tsamplesFor[gg]}, {x0, xsamplesFor[gg]}], 1];

(* orbit table of genus gg (prepared graphs + orbits); computed by lp2_prep.m and stored in lp2_orbits_g<gg>.mx *)
setupOrbits[gg_] := (prepared = prepareGraph /@ graphsWithLegs[gg, 0]; OT = orbitTable[prepared];
   Print["genus ", gg, ": ", Length[prepared], " graphs, ", Total[Length[#["decs"]] & /@ prepared], " labelled decorations, ", Length[OT["ids"]], " decorated graphs (orbits)"]);
loadOrbits[gg_] := Module[{f = "lp2_orbits_g" <> ToString[gg] <> ".mx"},
   If[FileExistsQ[f], Get[f], Print["orbit file ", f, " not found: computing"]; setupOrbits[gg]; DumpSave[f, {prepared, OT}]]];

ratl[v_] := Rationalize[Re[v], 10^-(PREC - 60)];
(* one sample point: for every orbit the coefficients of XX^0..XX^(3g-3); returns {{t0, x0}, list of vectors in the order OT["ids"], diagnostics} *)
computeSampleOrbits[t0_, x0_, gg_] := Module[{tt = AbsoluteTime[], data, wmax = 3 gg - 3, vecs, imax = 0, bad = {}},
   data = evalPoint[t0, x0, KMAXcur];
   vecs = Table[Module[{val = graphValueVE[OT["graphs"][id], data], rules, asc},
      rules = CoefficientRules[val, {C1i, C1p}]; asc = Association[rules];
      If[rules =!= {}, imax = Max[imax, Max[Abs[Im[rules[[All, 2]]]]]]];
      If[Select[rules, Abs[#[[2]]] > 10^-(PREC - 80) && (#[[1, 1]] != #[[1, 2]] || #[[1, 1]] > wmax) &] =!= {}, AppendTo[bad, id]];
      Table[ratl[Lookup[asc, Key[{k, k}], 0]/t0^(3 k)], {k, 0, wmax}]], {id, OT["ids"]}];
   {{t0, x0}, vecs, {"imag" -> N[imax, 5], "bad" -> bad, "seconds" -> Round[AbsoluteTime[] - tt, 0.1]}}];

(* ============================== 5. structured fits for one decorated graph ============================== *)
(* eps-direction: x^a v(x) is a polynomial of degree 2g-2+2a (eps2 = 1, eps1 = x); minimal a, verified on the spare points.
   Returns {a, coefficient list of length 2g-2+2a+1} or $Failed. *)
xfit[vals_, xs_, gg_] := Module[{res = $Failed, NX = Length[xs], deg, pts, pol, ok},
   If[Union[vals] === {0}, Return[{0, ConstantArray[0, 2 gg - 1]}]];
   Do[deg = 2 gg - 2 + 2 a; If[deg + 2 > NX, Break[]];
    pts = Table[{xs[[s]], xs[[s]]^a vals[[s]]}, {s, deg + 1}];
    pol = InterpolatingPolynomial[pts, xv];
    ok = And @@ Table[(pol /. xv -> xs[[s]]) == xs[[s]]^a vals[[s]], {s, deg + 2, NX}];
    If[ok, res = {a, CoefficientList[Expand[pol], xv, deg + 1]}; Return[]], {a, 0, 2 gg - 2 + 1}];
   res];
(* y-direction, several sequences (over the t-samples ts) sharing a denominator y^al (1+27y)^be: vectorised interpolation,
   5 spare t-values. Returns {al, be, list of rational functions of y} or $Failed. *)
yfitC[seqs_, ts_, gg_] := Module[{res = $Failed, NT = Length[ts], ys = ts^3, cands, nfit, A, W, sol, ok},
   nfit = NT - 5;
   cands = SortBy[Flatten[Table[{al, be}, {al, 0, 2 gg}, {be, 0, 3 gg - 3}], 1], {Total[#], #[[1]]} &];
   A = Table[ys[[s]]^d, {s, nfit}, {d, 0, nfit - 1}];
   Do[Module[{al = c[[1]], be = c[[2]]},
     W = Table[seqs[[q, s]] ys[[s]]^al (1 + 27 ys[[s]])^be, {q, Length[seqs]}, {s, NT}];
     sol = Transpose[LinearSolve[A, Transpose[W[[All, ;; nfit]]]]];
     ok = And @@ Flatten[Table[Sum[sol[[q, d + 1]] ys[[s]]^d, {d, 0, nfit - 1}] == W[[q, s]], {q, Length[seqs]}, {s, nfit + 1, NT}]];
     If[ok, res = {al, be, Table[Together[Sum[sol[[q, d + 1]] y^d, {d, 0, nfit - 1}]/(y^al (1 + 27 y)^be)], {q, Length[seqs]}]}; Return[]]], {c, cands}];
   res];
(* reconstruct one decorated graph from its sample matrix samp[t0][x0] = coefficient vector (XX^0..XX^(3g-3)); returns the
   exact term in y, XX, eps1, eps2 or {$Failed, info} *)
reconstructOrbit[samp_, ts_, xs_, gg_] := Module[{term = 0, e1 = \[Epsilon][1], e2 = \[Epsilon][2], wmax = 3 gg - 3, xf, avals, amax, deg, seqs, yf, ok = True, info = {}},
   Do[xf = Table[xfit[Table[samp[t0][x0][[k + 1]], {x0, xs}], xs, gg], {t0, ts}];
    If[MemberQ[xf, $Failed], ok = False; AppendTo[info, {k, "x"}]; Continue[]];
    avals = Union[xf[[All, 1]]]; amax = Max[avals]; deg = 2 gg - 2 + 2 amax;
    seqs = Table[Table[Module[{a = xf[[s, 1]], cl = xf[[s, 2]]}, Coefficient[Expand[xv^(amax - a) Sum[cl[[i + 1]] xv^i, {i, 0, Length[cl] - 1}]], xv, jj]], {s, Length[ts]}], {jj, 0, deg}];
    If[Union[Flatten[seqs]] === {0}, Continue[]];
    yf = yfitC[seqs, ts, gg];
    If[yf === $Failed, ok = False; AppendTo[info, {k, "y"}]; Continue[]];
    term += XX^k Sum[yf[[3, jj + 1]] e1^jj e2^(deg - jj), {jj, 0, deg}]/(e1 e2)^amax, {k, 0, wmax}];
   If[ok, term, {$Failed, info}]];

(* ============================== 6. conifold frame ============================== *)
(* t_c(Delta) from the analytic solution of theta^2 + 3y(3theta+1)(3theta+2) at Delta = 1+27y = 0; XX and y as series in tcv *)
tcSetup[NC_] := Module[{L2b, wc, ceq, csol, wcDel},
   L2b[f_] := Module[{th = y D[#, y] &}, (th[th[f]] + 3 y (9 th[th[f]] + 9 th[f] + 2 f))];
   wc = Sum[a[n] Del^n, {n, 0, NC}] /. Del -> 1 + 27 y;
   ceq = CoefficientList[Normal[Series[L2b[wc] /. y -> (Del - 1)/27, {Del, 0, NC - 1}]], Del];
   csol = Solve[Thread[ceq == 0] /. a[0] -> 1, Table[a[n], {n, 1, NC}]];
   wcDel = Sum[a[n] Del^n, {n, 0, NC}] /. First[csol] /. a[0] -> 1;
   tc = Module[{ser = Normal[Series[wcDel/(Del - 1), {Del, 0, NC}]]}, Sum[Coefficient[ser, Del, n] Del^(n + 1)/(n + 1), {n, 0, NC}]];
   tc = Normal[Series[tc/SeriesCoefficient[tc, {Del, 0, 1}], {Del, 0, NC}]];
   invcf = InverseSeries[Series[tc, {Del, 0, NC}], tcv];
   C1cf = Normal[y 27 D[tc, Del]] /. Del -> 1 + 27 y;
   XXcf = Series[y D[C1cf, y]/C1cf /. y -> (Del - 1)/27 /. Del -> invcf, {tcv, 0, NC - 2}];
   ycf = (invcf - 1)/27;
   NCcur = NC];
(* polar part t_c^-(2g-2) .. t_c^-1 of a term in y, XX, eps1, eps2 with eps1 = x, eps2 = 1 *)
polarOfTerm[term_, gg_] := Module[{ser},
   ser = Series[term /. \[Epsilon][1] -> x /. \[Epsilon][2] -> 1 /. XX -> XXcf /. y -> ycf, {tcv, 0, -1}];
   Table[Factor[SeriesCoefficient[ser, -n]], {n, 2 gg - 2, 1, -1}]];
(* the known gap coefficient: (-3)^(g-1) (2g-3)! [t^(2g)] t^2/(4 sinh(e1 t/2) sinh(e2 t/2)) *)
gapCoefficient[gg_] := Module[{tv}, (-3)^(gg - 1) (2 gg - 3)! SeriesCoefficient[Series[tv^2/(4 Sinh[\[Epsilon][1] tv/2] Sinh[\[Epsilon][2] tv/2]), {tv, 0, 2 gg}], 2 gg]];
