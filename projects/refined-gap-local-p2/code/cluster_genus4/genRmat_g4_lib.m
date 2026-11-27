(* genRmat_g4_lib.m -- library for the cluster computation of the twisted (generalised R-matrix) Givental-Teleman
   graph sum of K_P2 x C^2 on the fibrewise torus T_F, genus <= 4.  Load with Get after SetDirectory to the data directory.
   Sections 1-3 are verbatim copies of the verified Input cells of localp2CY5_genRmat_genus3.nb (core routines of
   localp2CY5_genus3.nb, the twisted buildRmatrix, the rational reconstruction); section 4 adds the genus-4 extensions.
   Conventions: y = t^3; eps1, eps2 are the symbols Epsilon[1], Epsilon[2] (typed as backslash-bracket-Epsilon); C1i = 1/c1, C1p = dc1/dy; XX = y C1p/c1; PREC digits. *)



(* ============================== 1. core routines (verbatim) ============================== *)
(* Part B core: fibrewise (T_F) Givental graph sum for K_P2 x C^2 by exact-structure sampling.
   Load after defining cPP[0, i_][k] for k = 0..KMAX (rational functions of y, LLL[i][y], eps).
   Everything is evaluated at numeric sample points (t0 = y^(1/3), x0 = eps1/eps2) in high precision,
   keeping C1i = 1/c1, C1p = c1' symbolic so that the dependence on XX = y c1'/c1 is separated. *)

PREC = 350;   (* default working precision for the cluster runs; 300 was used at genus <= 3 *)
ww = Exp[2 Pi I/3];

(* ---------- symbolic setup in t (i = 0 only; other idempotents by t -> w^i t) ---------- *)
alpha = (\[Epsilon][1] + \[Epsilon][2])/3;
eeF[1] = 3 alpha; eeF[2] = 3 alpha^2; eeF[3] = alpha^3;
LLt0 = -alpha/(1 + 3 t);                         (* LLL[0][y] on T_F with y = t^3 *)
sqrule = Power[(\[Epsilon][1] \[Epsilon][2])/(\[Epsilon][1] + \[Epsilon][2])^3, r_] :> sq^(2 r);
etaF = -{{-(1/(3 eeF[3] \[Epsilon][1] \[Epsilon][2])), 0, 0}, {0, 0, 1/(3 \[Epsilon][1] \[Epsilon][2])}, {0, 1/(3 \[Epsilon][1] \[Epsilon][2]), -(eeF[1]/(3 \[Epsilon][1] \[Epsilon][2]))}};
etainvF = Simplify[Inverse[etaF]];

(* total y-derivative acting on expressions in t, C1i = 1/c1[y], C1p = c1'[y] *)
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

(* ---------- numeric evaluation at a sample point ---------- *)
(* returns association with Dil[i][n], Leaf[i][n] (a=1), Edge[i,j][m,n], sqD[i] as polynomials in C1i, C1p *)
evalPoint[t0_, x0_, KMAX_] := Module[{sub, Pnum, i, a, k, pre, Dil, Lf, Ed, sqD, etainv, sqv, tt, zz, wv, cf},
   sqv = N[Sqrt[x0/(1 + x0)^3], PREC];
   etainv = N[etainvF /. {\[Epsilon][1] -> x0, \[Epsilon][2] -> 1}, PREC];
   Do[tt = N[ww^i t0, PREC];
    sub = {t -> tt, \[Epsilon][1] -> N[x0, PREC], \[Epsilon][2] -> N[1, PREC], sq -> sqv};
    Do[Pnum[a, i] = Expand[PPsym[a] /. sub], {a, 0, 2}];   (* list of z-coefficients *)
    sqD[i] = Expand[sqDsym /. sub];
    Do[Dil[i][n] = (-1)^n Pnum[0, i][[n]], {n, 1, KMAX + 1}];   (* Dilaton[i,n] = coeff z^(n-1) of -P0(-z) *)
    Do[Lf[i][n] = (-1)^n Pnum[1, i][[n + 1]], {n, 0, KMAX}];   (* Leaf[i,n,1] = coeff z^n of P1(-z) *)
    , {i, 0, 2}];
   Do[
    pre = Table[KroneckerDelta[i, j] KroneckerDelta[m, 0] KroneckerDelta[n, 0] -
       Sum[(-1)^(m + n) Pnum[a, i][[m + 1]] etainv[[a + 1, b + 1]] Pnum[b, j][[n + 1]], {a, 0, 2}, {b, 0, 2}],
      {m, 0, KMAX}, {n, 0, KMAX}];   (* pre[[m+1,n+1]] = coeff of w^m z^n in PreEdgeTerm[i,j] *)
    Do[Ed[i, j][m, n] = Expand[Sum[(-1)^k pre[[m + k + 2, n - k + 1]], {k, 0, n}]],
     {m, 0, KMAX - 1}, {n, 0, KMAX - 1 - m}];
    , {i, 0, 2}, {j, 0, 2}];
   etainvN = etainv; <|"Dil" -> Dil, "Leaf" -> Lf, "Edge" -> Ed, "sqD" -> sqD, "Pnum" -> Pnum|>];

(* ---------- Witten-Kontsevich numbers ---------- *)
df[n_] := If[n <= 0, 1, n df[n - 2]];
splitSum[g_, k_, rest_] := Module[{n = Length[rest], idx, tot = 0},
   Do[idx = Subsets[Range[n]][[m]];
    tot += Sum[df[2 r + 1] df[2 (k - 1 - r) + 1] Sum[wk[g1, Join[{r}, rest[[idx]]]] wk[g - g1, Join[{k - 1 - r}, rest[[Complement[Range[n], idx]]]]], {g1, 0, g}], {r, 0, k - 1}], {m, 2^n}]; tot];
wk[g_, ds_List] := wk[g, ds] = Module[{d = Sort[ds], n, k, rest, tot}, n = Length[d];
   If[g < 0 || n == 0 || Min[d] < 0 || Total[d] != 3 g - 3 + n || 2 g - 2 + n <= 0, Return[0]];
   If[g == 0 && d == {0, 0, 0}, Return[1]]; If[g == 1 && d == {1}, Return[1/24]];
   If[MemberQ[d, 0], rest = Delete[d, First[Position[d, 0]]]; Return[Sum[wk[g, ReplacePart[rest, j -> rest[[j]] - 1]], {j, Length[rest]}]]];
   k = d[[1]] - 1; rest = Rest[d];
   tot = Sum[df[2 k + 2 rest[[j]] + 1]/df[2 rest[[j]] - 1] wk[g, ReplacePart[rest, j -> rest[[j]] + k]], {j, Length[rest]}];
   tot += 1/2 Sum[df[2 r + 1] df[2 (k - 1 - r) + 1] wk[g - 1, Join[{r, k - 1 - r}, rest]], {r, 0, k - 1}];
   tot += 1/2 splitSum[g, k, rest]; tot/df[2 k + 3]];

(* ---------- stable graphs ---------- *)
(* graph = {genera, legs, adjacency matrix (diagonal = loops)} *)
tabtoadj[diag_, up_] := Module[{K = Length[diag], A}, A = ConstantArray[0, {K, K}];
   Do[A[[j, j]] = diag[[j]]; Do[A[[i, j]] = A[[j, i]] = up[[j, i - j]], {i, j + 1, K}], {j, 1, K}]; A];
graphlist[raw_] := Module[{count = 1, K, g, m, l, atab, a, out = {}, adj},
   While[count <= Length[raw],
    K = raw[[count]]; g = Take[raw, {count + 1, count + K}]; m = Take[raw, {count + K + 1, count + 2 K}];
    l = Take[raw, {count + 2 K + 1, count + 3 K}]; a = {};
    If[K > 1, atab = Take[raw, {count + 3 K + 1, count + 3 K + K (K - 1)/2}];
     Do[AppendTo[a, Take[atab, K - j]]; atab = Drop[atab, K - j], {j, 1, K - 1}]];
    adj = tabtoadj[l, a]; AppendTo[out, {g, m, adj}];
    count += 3 K + Length[Flatten[a]] + 1]; out];
loadGraphs[gg_] := graphlist[DeleteCases[Delete[Delete[BinaryReadList["StableGraph-" <> ToString[gg] <> ".bin"], 1], -1], 13]];
(* hand-coded graphs with legs for the genus <= 1 checks *)
graphsWithLegs[0, 3] = {{{0}, {3}, {{0}}}};
graphsWithLegs[0, 4] = {{{0}, {4}, {{0}}}, {{0, 0}, {2, 2}, {{0, 1}, {1, 0}}}};
graphsWithLegs[0, 5] = {{{0}, {5}, {{0}}}, {{0, 0}, {2, 3}, {{0, 1}, {1, 0}}}, {{0, 0, 0}, {2, 1, 2}, {{0, 1, 0}, {1, 0, 1}, {0, 1, 0}}}};
graphsWithLegs[1, 1] = {{{1}, {1}, {{0}}}, {{0}, {1}, {{1}}}};
graphsWithLegs[g_, 0] := loadGraphs[g];

edgesOf[A_] := Join[Flatten[Table[ConstantArray[{i, i}, A[[i, i]]], {i, Length[A]}], 1],
   Flatten[Table[ConstantArray[{i, j}, A[[i, j]]], {i, Length[A]}, {j, i + 1, Length[A]}], 2]];
autCount[gv_, mv_, A_] := Module[{V = Length[gv], perms},
   perms = Select[Permutations[Range[V]], gv[[#]] == gv && mv[[#]] == mv && A[[#, #]] == A &];
   Length[perms] Product[A[[i, j]]!, {i, 1, V}, {j, i + 1, V}] Product[2^A[[i, i]] A[[i, i]]!, {i, 1, V}]];

(* vertex decoration options: for vertex with genus g, hs = number of half-edges, m = legs:
   list of {hpowers, legpowers, dilpowers, weight}, weight = wk * 1/(m! nd!) *)
vertexOptions[g_, hs_, m_] := Module[{res = {}, dim, nd, dils, hp, lp, rest, budget},
   Do[dim = 3 g - 3 + hs + m + nd;
    If[dim < 2 nd, Continue[]];
    Do[budget = dim - Total[dils];
     Do[Do[
       AppendTo[res, {hp, lp, dils, wk[g, Join[hp, lp, dils]]/(m! nd!)}],
       {lp, Select[Tuples[Range[0, budget - Total[hp]], m], Total[#] == budget - Total[hp] &]}],
      {hp, Select[Tuples[Range[0, budget], hs], Total[#] <= budget &]}],
     {dils, If[nd == 0, {{}}, Flatten[Table[Permutations /@ Select[IntegerPartitions[s, {nd}], Min[#] >= 2 &], {s, 2 nd, dim}], 2]]}],
    {nd, 0, Max[0, 3 g - 3 + hs + m]}];
   Select[res, #[[4]] =!= 0 &]];

(* prepare a graph: half-edge bookkeeping and decoration list *)
prepareGraph[{gv_, mv_, A_}] := Module[{V = Length[gv], ed, H, hv, hsAt, opts, decs, aut},
   ed = edgesOf[A]; H = Flatten[Table[{{e, 1}, {e, 2}}, {e, Length[ed]}], 1];
   hv = Table[ed[[h[[1]], h[[2]]]], {h, H}];
   hsAt = Table[Flatten[Position[hv, v]], {v, V}];   (* indices into H of half-edges at v *)
   opts = Table[vertexOptions[gv[[v]], Length[hsAt[[v]]], mv[[v]]], {v, V}];
   decs = Tuples[opts];
   aut = autCount[gv, mv, A];
   <|"g" -> gv, "m" -> mv, "A" -> A, "edges" -> ed, "H" -> H, "hsAt" -> hsAt, "decs" -> decs, "aut" -> aut, "V" -> V|>];

(* evaluate one prepared graph at a sample point (data from evalPoint) *)
graphValue[gr_, data_] := Module[{V = gr["V"], ed = gr["edges"], H = gr["H"], hsAt = gr["hsAt"], Dil = data["Dil"], Lf = data["Leaf"], Ed = data["Edge"], sqD = data["sqD"], total = 0, hp, c, term, vfac, efac},
   Do[  (* colourings *)
    Do[ (* decorations: dec = list over vertices of {hp, lp, dils, w} *)
     hp = ConstantArray[0, Length[H]];
     Do[Do[hp[[hsAt[[v, r]]]] = dec[[v, 1, r]], {r, Length[hsAt[[v]]]}], {v, V}];
     vfac = Product[dec[[v, 4]] sqD[c[[v]]]^(2 gr["g"][[v]] - 2 + Length[hsAt[[v]]] + gr["m"][[v]] + Length[dec[[v, 3]]])
        Product[Lf[c[[v]]][b], {b, dec[[v, 2]]}] Product[Dil[c[[v]]][d], {d, dec[[v, 3]]}], {v, V}];
     efac = Product[Ed[c[[ed[[e, 1]]]], c[[ed[[e, 2]]]]][hp[[2 e - 1]], hp[[2 e]]], {e, Length[ed]}];
     total += vfac efac;
     , {dec, gr["decs"]}];
    , {c, Tuples[{0, 1, 2}, V]}];
   Expand[total/gr["aut"]]];

(* full graph sum F_{g,n} at a sample point *)
graphSum[prepared_, data_] := Expand[Total[graphValue[#, data] & /@ prepared]];

(* ============================== 2. twisted R-matrix (verbatim; qsyms etc. are redefined in section 4) ============================== *)
(* The twist. ff[z] = Exp[q0 z + q1 z^3 + q2 z^5 + ...] (so ff[z] ff[-z] = 1), with q_k = Qrrc[0,k] of
   localp2CY5_genRmat_v3.nb kept as symbols throughout the numeric evaluation. Only the identity row of the
   R-matrix is multiplied by ff[z] (truncated at z^KMAX); the rows 1, 2 follow from the flatness recursion exactly
   as in the untwisted buildRmatrix above, which this definition replaces. *)
qsyms = {q0, q1, q2};
qrule = Table[qsyms[[k + 1]] -> Qrrc[0, k], {k, 0, 2}];
vars = Join[{C1i, C1p}, qsyms];
(* q-monomials q0^a q1^b q2^c of weight a + 3b + 5c <= wmax, as {{a,b,c}, weight} *)
qmonos[wmax_] := SortBy[Select[Flatten[Table[{{a, b, c}, a + 3 b + 5 c}, {a, 0, wmax}, {b, 0, Floor[wmax/3]}, {c, 0, Floor[wmax/5]}], 2], #[[2]] <= wmax &], {#[[2]], #[[1]]} &];
qmono[e_List] := Times @@ (qsyms^e);
coeffs[expr_] := Module[{e = Expand[expr]}, If[e === 0, {}, If[Head[e] === Plus, List @@ e, {e}] /. Thread[vars -> 1]]];   (* numeric coefficients of the monomials in vars (negative powers of C1i allowed) *)
maxabs[expr_] := If[Expand[expr] === 0, 0, Max[Abs[coeffs[expr]]]];

buildRmatrix[KMAX_] := Module[{P0, P1, P2, c2, d2, k, ffz},
   cPPtk[k_] := cPPtk[k] = Together[cPP[0, 0][k] /. LLL[0][y] -> LLt0 /. y -> t^3 /. sqrule];
   ffz = Normal[Series[Exp[Sum[qsyms[[k + 1]] z^(2 k + 1), {k, 0, Floor[(KMAX - 1)/2]}]], {z, 0, KMAX}]];
   P0 = Normal[Series[Sum[cPPtk[k] z^k, {k, 0, KMAX}] ffz, {z, 0, KMAX}]];
   P1 = Normal[Series[(LLt0 P0 + z t^3 Dy[P0]) C1i, {z, 0, KMAX}]];
   c2 = C1i^2/(1 + 27 t^3); d2 = (\[Epsilon][1] + \[Epsilon][2]) (C1i^2 - 1)/(2 (1 + 27 t^3));
   P2 = Normal[Series[(LLt0 P1 + z t^3 Dy[P1] - d2 P1)/c2, {z, 0, KMAX}]];
   If[! FreeQ[{P1, P2}, c1pp], Print["WARNING: second derivative of c1 appeared"]];
   PPsym[0] = Table[Together[Coefficient[P0, z, k]], {k, 0, KMAX}];
   PPsym[1] = Table[Together[Coefficient[P1, z, k]], {k, 0, KMAX}];
   PPsym[2] = Table[Together[Coefficient[P2, z, k]], {k, 0, KMAX}];
   sqDsym = Sqrt[27] sq (3 eeF[3] + eeF[2] LLt0);
   KMAXcur = KMAX;
   Print["twisted R-matrix rows built to z^", KMAX, " (ff[z] to q", Floor[(KMAX - 1)/2], "); LeafCounts row0/1/2: ", LeafCount /@ {PPsym[0], PPsym[1], PPsym[2]}];
   ];

(* ============================== 3. rational reconstruction (verbatim) ============================== *)
(* generic rational reconstruction from exact samples {u_s, v_s}: N(u)/D(u), D monic, verified on >= 5 spare points *)
reconstruct[pts_] := Module[{n = Length[pts], us = pts[[All, 1]], vs = pts[[All, 2]], res = $Failed, T, dd, dn, m, A, rhs, sol, N0, D0, ok},
   If[Union[vs] === {0}, Return[0]];
   Do[
    Do[dn = T - dd; m = dn + dd + 1;
     If[m + 5 > n, Continue[]];
     (* unknowns: n0..n_dn, d0..d_{dd-1}; equations N(u_s) - v_s (u_s^dd + sum d_j u_s^j) = 0 *)
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

(* ============================== 4. cluster extensions ============================== *)
(* ======================================================================================================
   Extensions for the cluster version (genus 4 and below):
   - five twist coefficients q0..q4 (q_k multiplies z^(2k+1) in ff[z]; to order z^9 no further one enters),
   - the R-matrix identity column is read from the portable text file cPP_z<KMAX>.m (list of cPPL[k], k = 0..KMAX,
     functions of y, LLL[0][y], eps, sq; cPP_z6.m is the content of cPPu6.mx, cPP_z9.m must be produced with
     genRmat_g4_partA_pf.m / genRmat_g4_partA_limit.m),
   - the sum over the 3^V colourings of a decorated graph is done by variable elimination (contractColourings)
     instead of enumeration: identical result, far fewer polynomial multiplications for V = 5, 6,
   - sample-point bookkeeping for independent cluster jobs (computeSample, pointList).
   ====================================================================================================== *)
qsyms = {q0, q1, q2, q3, q4};
qweights = {1, 3, 5, 7, 9};
qrule = Table[qsyms[[k + 1]] -> Qrrc[0, k], {k, 0, 4}];
vars = Join[{C1i, C1p}, qsyms];
qmonos[wmax_] := SortBy[Select[Flatten[Table[{{a, b, c, d, e}, a + 3 b + 5 c + 7 d + 9 e}, {a, 0, wmax}, {b, 0, Floor[wmax/3]}, {c, 0, Floor[wmax/5]}, {d, 0, Floor[wmax/7]}, {e, 0, Floor[wmax/9]}], 4], #[[2]] <= wmax &], {#[[2]], #[[1]]} &];
qweight[e_List] := e . qweights;
qmono[e_List] := Times @@ (qsyms^e);

loadRmatrix[KMAX_] := Module[{avail, file, lst},
   (* use cPP_z<K>.m with the smallest K >= KMAX available (cPP_z6.m serves genus <= 3, cPP_z9.m is needed for genus 4) *)
   avail = Select[ToExpression[StringReplace[FileNameTake /@ FileNames["cPP_z*.m"], {"cPP_z" -> "", ".m" -> ""}]], IntegerQ[#] && # >= KMAX &];
   If[avail === {}, Print["ERROR: no cPP_z<K>.m with K >= ", KMAX, " in ", Directory[]]; Abort[]];
   file = "cPP_z" <> ToString[Min[avail]] <> ".m"; lst = Get[file];
   If[! ListQ[lst] || Length[lst] < KMAX + 1, Print["ERROR: ", file, " does not contain ", KMAX + 1, " coefficients"]; Abort[]];
   Print["R-matrix coefficients from ", file];
   Do[cPP[0, 0][k] = lst[[k + 1]], {k, 0, KMAX}];
   buildRmatrix[KMAX]];

(* ---------- colouring sum by variable elimination ----------
   vf[[v, c+1]]: vertex factor of vertex v with colour c;  ef[[e, c1+1, c2+1]]: factor of edge e = {v1, v2} (v1 <= v2; loops v1 == v2 use the diagonal).
   Returns sum over all colourings c in {0,1,2}^V of Product[vf[[v, c_v]]] Product[ef[[e, c_v1, c_v2]]]. *)
contractColourings[V_, ed_, vf_, ef_] := Module[{factors, remaining = Range[V], v, cpos, cand, U, keep, posU, idx, keepidx, newassoc, key, val, sizes},
   factors = Join[
     Table[{{w}, Association[Table[{c} -> vf[[w, c + 1]], {c, 0, 2}]]}, {w, V}],
     Table[If[ed[[e, 1]] == ed[[e, 2]],
        {{ed[[e, 1]]}, Association[Table[{c} -> ef[[e, c + 1, c + 1]], {c, 0, 2}]]},
        {ed[[e]], Association[Flatten[Table[{c1, c2} -> ef[[e, c1 + 1, c2 + 1]], {c1, 0, 2}, {c2, 0, 2}]]]}], {e, Length[ed]}]];
   While[remaining =!= {},
    (* greedy: eliminate the vertex whose merged factor has the fewest variables *)
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
    total += contractColourings[V, ed, vf, ef];
    , {dec, gr["decs"]}];
   Expand[total/gr["aut"]]];
graphSumVE[prepared_, data_] := Expand[Total[graphValueVE[#, data] & /@ prepared]];

(* ---------- sample points and jobs ---------- *)
tsamplesFor[gg_] := If[gg <= 3,
   {1/2, 1/3, 2/3, 1/4, 3/4, 1/5, 2/5, 3/5, 4/5, 1/6, 5/6, 1/7, 2/7, 3/7, 4/7, 5/7, 6/7, 1/8, 3/8, 5/8, 7/8, 1/9, 2/9},
   {1/2, 1/3, 2/3, 1/4, 3/4, 1/5, 2/5, 3/5, 4/5, 1/6, 5/6, 1/7, 2/7, 3/7, 4/7, 5/7, 6/7, 1/8, 3/8, 5/8, 7/8, 1/9, 2/9, 1/10, 3/10, 7/10, 9/10}];
xsamplesFor[gg_] := Table[Prime[r + 1], {r, 1, 2 gg + (3 gg - 3) + 3}];   (* 2g + w + 1 needed for q-weight w <= 3g-3, plus 2 checks *)
pointList[gg_] := Flatten[Table[{t0, x0}, {t0, tsamplesFor[gg]}, {x0, xsamplesFor[gg]}], 1];

setupGenus[gg_, KMAX_] := (GG = gg; loadRmatrix[KMAX]; wmax = 3 gg - 3; mons = qmonos[wmax];
   prepared = prepareGraph /@ graphsWithLegs[gg, 0];
   Print["genus ", gg, ": ", Length[prepared], " graphs, ", Total[Length[#["decs"]] & /@ prepared], " labelled decorations, ", Length[mons], " q-monomials of weight <= ", wmax, ", ", Length[pointList[gg]], " sample points"]);

ratl[v_] := Rationalize[Re[v], 10^-(PREC - 60)];
(* one sample point: returns {{t0, x0}, matrix (rows k = 0..wmax = power of XX, columns = q-monomials in the order of mons), diagnostics} *)
computeSample[t0_, x0_] := Module[{tt = AbsoluteTime[], data, val, rules, asc, imax, bad, mat},
   data = evalPoint[t0, x0, KMAXcur];
   val = graphSumVE[prepared, data];
   rules = CoefficientRules[val, vars]; asc = Association[rules];
   imax = If[rules === {}, 0, Max[Abs[Im[rules[[All, 2]]]]]];
   bad = Select[rules, Abs[#[[2]]] > 10^-(PREC - 80) && (#[[1, 1]] != #[[1, 2]] || #[[1, 1]] > wmax || qweight[#[[1, 3 ;;]]] > wmax) &][[All, 1]];
   mat = Table[ratl[Lookup[asc, Key[Join[{k, k}, mons[[m, 1]]]], 0]/t0^(3 k)], {k, 0, wmax}, {m, Length[mons]}];
   {{t0, x0}, mat, {"imag" -> N[imax, 5], "bad" -> bad, "seconds" -> Round[AbsoluteTime[] - tt, 0.1]}}];