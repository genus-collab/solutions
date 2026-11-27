(* Plain-text export of localp2CY5_hodge_genus3.nb (2026-11-27): Input cells as code, text cells as comments; outputs omitted. *)

(* ::Title:: Local P2 x C2: Hodge-vertex graph sum on the fibrewise torus, genus 2, 3 (and 4 from the cluster) *)

(* ::Text:: Successor of localp2CY5_genus3.nb. The T_F-equivariant Gromov-Witten theory of Z = K_P2 x C^2 is reconstructed from three copies of the Hodge CohFT (the equivariant GW theory of C^3, vertex classes prod_a Lambda_g^vee(w_a)/w_a) with the e_T weights of the normal bundle O(-3) + O + O of the fixed P^2, w = {eps1+eps2, -eps1, -eps2}, instead of three copies of the GW theory of the point: Omega_Z = R'.(Hodge_w)^{+3} with R' = R_Z h_w^(-1), h_w(z) = Exp[-Sum_k B_2k/(2k(2k-1)) z^(2k-1) Sum_a w_a^(1-2k)] the scalar Bernoulli series (h_w . point = Hodge_w by Mumford's formula / Givental's quantum Riemann-Roch). R_Z is the R-matrix of localp2CY5_genus3.nb (identity column cPP_z6.m, rows by the T_R ansatz), so the totals F_{2,0}, F_{3,0} are unchanged (group-action property) while their distribution over the stable graphs changes. Evaluate top to bottom. Data files (in the notebook directory or its subfolder hodge_vertex/): cPP_z6.m, StableGraph-2.bin, StableGraph-3.bin, ref/F20_F30_y.m, hv_orbits_g3.mx (decorated graphs of genus 3, recomputed if absent), hv3_samples_B.mx (the 780 genus-3 sample records of the stored run; delete it to recompute, several hours on one kernel), hodge_vertex_genus2_results.m (stored genus-2 run for the regression). *)

(* ::Text:: Conventions. Vertices carry idempotent labels i in {0,1,2}; on T_F the canonical coordinates are LLL[i][y] = -alpha/(1 + 3 omega^i t), t = y^(1/3), alpha = (eps1+eps2)/3. Graph factors as in localp2CY5_genus3.nb (Dilaton[i,n], Leaf, EdgeTerm, sqrtDelta[i] = Sqrt[eps1 eps2/e3](3 e3 + e2 LLL[i])) but built from the rows of R' = R_Z h^(-1) (coefficient lists of R_Z times the coefficient list of h^(-1)), and the vertex factor sqrtDelta^(2g_v-2+|H(v)|) <prod tau> is replaced by sqrtDelta^(2g_v-2+|H(v)|) Int_{M_{g_v,n_v}} prod_a Lambda^vee_{g_v}(w_a)/w_a^{g_v} prod psi^a, so every psi assignment of degree at most the dimension contributes (the lambda classes make up the difference). A decorated graph is an orbit of labelled psi/dilaton decorations under Aut(Gamma); each stable graph contributes 1/|Aut(Gamma)| times the sum over labelled decorations and 3^V colourings, with 1/nd_v! for the dilaton leaves. *)

(* ::Text:: Hodge integrals (genus <= 3). The normalised class prod_a Lambda_g^vee(w_a)/w_a^g is expanded in lambda monomials with the Mumford relations lambda_1^2 = 2 lambda_2, lambda_2^2 = 2 lambda_1 lambda_3, lambda_3^2 = 0; the integrals against psi classes are Int lambda_1^k psi^a by k-fold Mumford's formula lambda_1 = (kappa_1 - Sum psi_i + delta)/12 (kappa_1 = pi_* psi^2, <kappa_1 X> = <tau_2 X>, lambda_1 restricts to the components of the boundary), the lambda_g formula, the lambda_g lambda_{g-1} closed form (2g-3+n)! |B_2g|/(2^(2g-1) (2g)! prod (2a_i-1)!!) for a_i >= 1 with zeros removed by the string equation, and lambda_1 lambda_2 lambda_3 on M_3 = |B_4 B_6|/(2 4 6 4!). Checks: known values, the recursion against the closed forms (lambda_1^5 = 16 lambda_2 lambda_3, lambda_1^6 = 16 lambda_1 lambda_2 lambda_3 with psi insertions, lambda_1^6 on M_3 = 1/90720), string and dilaton equations, and the R-action of the scalar h on the point theory in the PPZ conventions reproducing the normalised Hodge class on M_2 (symbolically) and M_3 (42 graphs, exact rational weights). *)

(* ::Text:: Numerics. As in localp2CY5_genus3.nb: rational sample points (t0, x0 = eps1/eps2) with 350-digit arithmetic, C1i = 1/c1 and C1p = c1' symbolic so that the XX-dependence is separated; exact rational results by Rationalize, reconstructed as polynomials in XX with coefficients rational in y (denominators y^a (1+27y)^b) and homogeneous rational functions of (eps1, eps2) with denominators (eps1 eps2)^a (eps1+eps2)^b, every fit verified on spare points. The sum over the decorated graphs is compared with F_{2,0}, F_{3,0} at every point and symbolically; the total conifold polar part is compared with the gap coefficient. *)

(* ::Section:: Core routines *)

(* ::Subsection:: Stable graphs, decorations, Witten-Kontsevich numbers (graphs.m) *)

(* graphs.m -- stable graphs, decorations and Witten-Kontsevich numbers (shared by lp2_lib.m and rc_lib.m).
   Verbatim from genRmat_g4_lib.m / localp2CY5_genus3.nb / resolvedconifold_genus3.nb.
   graph = {genera, legs, adjacency matrix (diagonal = loops)}; StableGraph-<g>.bin must be in the current directory. *)

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
(* hand-coded graphs with legs for the genus <= 1 cases *)
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

(* vertex decoration options: for a vertex with genus g, hs half-edges, m legs:
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

(* ::Subsection:: Decorated graphs as Aut-orbits (orbits.m) *)

(* orbits.m -- decorated stable graphs = orbits of labelled decorations under Aut(Gamma). Shared by both pipelines
   (identical prepareGraph structure on both sides, so the keys match across theories).
   A decorated graph is the orbit of {psi powers on the two half-edges of each edge, sorted dilaton lists per vertex}
   under vertex permutations preserving genera/legs/adjacency, permutations of parallel edges and flips of loops.
   Canonical representative = lexicographically smallest image. *)

autosOf[gr_] := Module[{gv = gr["g"], mv = gr["m"], A = gr["A"], ed = gr["edges"], V = gr["V"], perms, out = {}},
   perms = Select[Permutations[Range[V]], gv[[#]] == gv && mv[[#]] == mv && A[[#, #]] == A &];
   Do[Module[{img = Map[Sort[sigma[[#]]] &, ed], classes, choices},
     classes = GatherBy[Range[Length[ed]], img[[#]] &];
     choices = Tuples[Table[Module[{cl = c, targets = Flatten[Position[ed, img[[c[[1]]]], {1}]]}, Map[Thread[cl -> #] &, Permutations[targets]]], {c, classes}]];
     Do[Do[AppendTo[out, {sigma, Flatten[ch], flips}], {flips, Tuples[Table[If[ed[[e, 1]] == ed[[e, 2]], {0, 1}, {0}], {e, Length[ed]}]]}], {ch, choices}]], {sigma, perms}]; out];
decData[gr_, dec_] := Module[{ed = gr["edges"], hsAt = gr["hsAt"], psi, dils},
   psi = Table[{0, 0}, Length[ed]];
   Do[Do[psi[[gr["H"][[hsAt[[v, r]], 1]], gr["H"][[hsAt[[v, r]], 2]]]] = dec[[v, 1, r]], {r, Length[hsAt[[v]]]}], {v, gr["V"]}];
   dils = Table[Sort[dec[[v, 3]]], {v, gr["V"]}]; {psi, dils}];
applyAut[{psi_, dils_}, {sigma_, emap_, flips_}, ed_] := Module[{psi2 = Table[{0, 0}, Length[ed]], dils2 = dils},
   Do[Module[{e2 = e /. emap, a = psi[[e]], tgt, ord}, tgt = {sigma[[ed[[e, 1]]]], sigma[[ed[[e, 2]]]]};
     ord = If[tgt == ed[[e2]], {1, 2}, {2, 1}]; If[ed[[e2, 1]] == ed[[e2, 2]] && flips[[e]] == 1, ord = Reverse[ord]];
     psi2[[e2, ord[[1]]]] = a[[1]]; psi2[[e2, ord[[2]]]] = a[[2]]], {e, Length[ed]}];
   Do[dils2[[sigma[[v]]]] = dils[[v]], {v, Length[dils]}]; {psi2, dils2}];
canonOf[gr_, dec_, auts_] := First[Sort[Table[applyAut[decData[gr, dec], a, gr["edges"]], {a, auts}]]];
(* {orbitKeys, orbits}: orbitKeys[[k]] = canonical key of every labelled decoration of graph k; orbits[[k]] = sorted distinct keys *)
computeOrbits[prep_] := Module[{orbitKeys, orbits},
   orbitKeys = Table[Module[{auts = autosOf[gr]}, Table[canonOf[gr, dec, auts], {dec, gr["decs"]}]], {gr, prep}];
   orbits = Table[Union[orbitKeys[[k]]], {k, Length[prep]}];
   {orbitKeys, orbits}];
orbitDecs[prep_, orbitKeys_, k_, okey_] := prep[[k]]["decs"][[Flatten[Position[orbitKeys[[k]], okey, {1}]]]];
orbitDescription[gr_, okey_] := <|"genera" -> gr["g"], "legs" -> gr["m"], "adjacency" -> gr["A"], "edges" -> gr["edges"], "psiOnHalfEdges" -> okey[[1]], "dilatonLeaves" -> okey[[2]], "Aut" -> gr["aut"], "dilatonFree" -> (Union[Flatten[okey[[2]]]] === {})|>;
(* full orbit table for a list of prepared graphs: orbitIds = {{k, j}, ...}, orbitGraphs[{k, j}] = graph restricted to the orbit *)
orbitTable[prep_] := Module[{orbitKeys, orbits, ids, graphs},
   {orbitKeys, orbits} = computeOrbits[prep];
   ids = Flatten[Table[{k, j}, {k, Length[prep]}, {j, Length[orbits[[k]]]}], 1];
   graphs = Association[Table[id -> <|prep[[id[[1]]]], "decs" -> orbitDecs[prep, orbitKeys, id[[1]], orbits[[id[[1]], id[[2]]]]]|>, {id, ids}]];
   <|"orbitKeys" -> orbitKeys, "orbits" -> orbits, "ids" -> ids, "graphs" -> graphs,
     "descriptions" -> Association[Table[id -> orbitDescription[prep[[id[[1]]]], orbits[[id[[1]], id[[2]]]]], {id, ids}]]|>];

(* ::Subsection:: R-matrix rows, colour sums, fits, conifold frame (lp2_lib.m) *)

(* lp2_lib.m -- library for the per-decorated-graph Givental-Teleman graph sum of K_P2 x C^2 on the fibrewise torus T_F
   (untwisted R-matrix), genus <= 4.  Load with Get after SetDirectory to the data directory (graphs.m, orbits.m,
   cPP_z<K>.m, StableGraph-<g>.bin must be there).
   Sections 1 (core routines), 3 (rational reconstruction) and the colouring sum by variable elimination are verbatim
   from genRmat_g4_lib.m (cluster_genus4/, itself verbatim from localp2CY5_genus3.nb); the twist is removed.
   Conventions: y = t^3; eps1, eps2 are the symbols Epsilon[1], Epsilon[2] (typed as backslash-bracket-Epsilon);
   C1i = 1/c1, C1p = dc1/dy; XX = y C1p/c1; PREC digits. ASCII only. *)

PREC = 350;
ww = Exp[2 Pi I/3];

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

(* ::Subsection:: Hodge integrals, Hodge vertices, R' = R_Z h^(-1) (hv_lib3.m) *)

(* hv_lib3.m -- genus-general (g <= 3) Hodge-vertex machinery for the Givental-Teleman graph sum of K_P2 x C^2 on T_F.
   Omega_Z = R' . (Hodge_w)^{+3}, R' = R_Z h_w^(-1), h_w(z) = Exp[-Sum_k B_2k/(2k(2k-1)) z^(2k-1) Sum_a w_a^(1-2k)].
   Vertex = normalised Hodge class prod_a Lambda_g^vee(w_a)/w_a^g expanded in lambda monomials (Mumford relations) with
   coefficients polynomial in the inverse weights u[a] = 1/w_a; the lambda integrals with psi insertions come from
   (i) Int lambda_1^k psi^a by k-fold Mumford's formula lambda_1 = (kappa_1 - Sum psi_i + delta)/12,
   (ii) the lambda_g formula, (iii) lambda_g lambda_{g-1} (closed form for a_i >= 1, zeros by the string equation),
   (iv) lambda_g lambda_{g-1} lambda_{g-2} in genus 3 (dilaton from Int_{M_3} = |B_4 B_6|/(2 4 6 4!)).
   Built on lp2_lib.m (R-matrix rows, colour sums, fits) and graphs.m / orbits.m.  ASCII only. *)

e1s = \[Epsilon][1]; e2s = \[Epsilon][2];
multin[a_List] := Total[a]!/Times @@ (a!);
bg[g_] := Abs[BernoulliB[2 g]];
lamTop[g_] := g (g + 1)/2;

(* ============================== 1. lambda integrals ============================== *)
(* Int_{M_{g,n}} lambda_1^k psi^a by Mumford: <kappa_1 X> = <tau_2 X> (kappa_1 = pi_* psi_{n+1}^2 and lambda_1 pulls back);
   boundary: 1/2 irreducible (genus g-1, two new points, lambda_1 restricts to lambda_1) + every separating divisor once
   (= 1/2 sum over ordered (h,S); lambda_1 restricts to lambda_1' + lambda_1''); unstable components give 0 *)
LpowRaw[g_, k_, a_List] /; a =!= Sort[a] := LpowRaw[g, k, Sort[a]];
LpowRaw[g_, k_, a_List] := LpowRaw[g, k, a] = Module[{n = Length[a], tot, Sc},
   Which[g < 0 || k < 0, 0, 2 g - 2 + n <= 0, 0, Total[a] + k != 3 g - 3 + n, 0, k == 0, wk[g, a],
    True,
    tot = LpowRaw[g, k - 1, Append[a, 2]] - Sum[LpowRaw[g, k - 1, ReplacePart[a, i -> a[[i]] + 1]], {i, n}] + LpowRaw[g - 1, k - 1, Join[a, {0, 0}]]/2;
    tot += 1/2 Sum[Sc = Complement[Range[n], S];
       Sum[Binomial[k - 1, j] LpowRaw[h, j, Append[a[[S]], 0]] LpowRaw[g - h, k - 1 - j, Append[a[[Sc]], 0]], {j, 0, k - 1}], {h, 0, g}, {S, Subsets[Range[n]]}];
    tot/12]];
(* with the vanishing lambda_1^k = 0 for k > g(g+1)/2 built in (used by the vertex weights) *)
Lpow[g_, k_, a_List] := If[k > lamTop[g], 0, LpowRaw[g, k, a]];
(* lambda_g psi^a (Faber-Pandharipande): multinomial (2g-3+n; a) (2^(2g-1)-1)/2^(2g-1) |B_2g|/(2g)! *)
LgTop[g_, a_List] := If[g >= 1 && 2 g - 2 + Length[a] > 0 && Total[a] == 2 g - 3 + Length[a], multin[a] (2^(2 g - 1) - 1)/2^(2 g - 1) bg[g]/(2 g)!, 0];
(* lambda_g lambda_{g-1} psi^a (psi degree g-2+n): (2g-3+n)! |B_2g| / (2^(2g-1) (2g)! prod (2a_i-1)!!) for all a_i >= 1;
   a zero is removed by the string equation whenever the target space is stable *)
Lg1[g_, a_List] := Module[{n = Length[a], z, rest},
   If[g < 1 || 2 g - 2 + n <= 0 || Total[a] != g - 2 + n, Return[0]];
   z = Position[a, 0, {1}, Heads -> False];
   If[z =!= {} && 2 g - 2 + n - 1 > 0, rest = Delete[a, z[[1]]];
    Return[Sum[If[rest[[j]] >= 1, Lg1[g, ReplacePart[rest, j -> rest[[j]] - 1]], 0], {j, Length[rest]}]]];
   (2 g - 3 + n)! bg[g]/(2^(2 g - 1) (2 g)! Times @@ (df[2 # - 1] & /@ a))];
(* lambda_3 lambda_2 lambda_1 psi^a on M_{3,n} (psi degree n): all a_i = 1 by the string equation, then dilaton *)
Lg2[3, a_List] := Module[{n = Length[a], z, rest},
   If[Total[a] != n, Return[0]];
   z = Position[a, 0, {1}, Heads -> False];
   If[z =!= {}, rest = Delete[a, z[[1]]]; Return[Sum[If[rest[[j]] >= 1, Lg2[3, ReplacePart[rest, j -> rest[[j]] - 1]], 0], {j, Length[rest]}]]];
   (n + 3)!/3!/1451520];
Lg2[g_, a_List] := 0 /; g != 3;

(* ============================== 2. the normalised Hodge class ============================== *)
(* prod_a Lambda_g^vee(w_a)/w_a^g = prod_a Sum_j (-1)^j lambda_j u_a^j reduced by c(E)c(E^vee) = 1:
   lambda_1^2 = 2 lambda_2, lambda_2^2 = 2 lambda_1 lambda_3, lambda_3^2 = 0 (genus 3), lambda_1^2 = 2 lambda_2, lambda_2^2 = 0 (genus 2),
   lambda_1^2 = 0 (genus 1).  Result: list of {exponents {i,j,k} of lambda_1^i lambda_2^j lambda_3^k, coefficient in u[1], u[2], u[3]} *)
mumfordRules[g_] := Switch[g,
   0, {lam[1] -> 0, lam[2] -> 0, lam[3] -> 0},
   1, {lam[1]^k_ /; k >= 2 -> 0, lam[2] -> 0, lam[3] -> 0},
   2, {lam[1]^k_ /; k >= 2 :> 2 lam[2] lam[1]^(k - 2), lam[2]^k_ /; k >= 2 -> 0, lam[3] -> 0},
   3, {lam[1]^k_ /; k >= 2 :> 2 lam[2] lam[1]^(k - 2), lam[2]^k_ /; k >= 2 :> 2 lam[1] lam[3] lam[2]^(k - 2), lam[3]^k_ /; k >= 2 -> 0}];
hodgeClass[g_] := hodgeClass[g] = Module[{cl, rules = mumfordRules[g]},
   cl = Expand[Product[Sum[(-1)^j lam[j] u[a]^j, {j, 0, g}], {a, 3}] /. lam[0] -> 1];
   cl = FixedPoint[Expand[# /. rules] &, cl];
   {#[[1]], Factor[#[[2]]]} & /@ CoefficientRules[cl, {lam[1], lam[2], lam[3]}]];
(* integral of a lambda monomial against psi^a *)
lamInt[g_, {0, 0, 0}, a_] := wk[g, a];
lamInt[g_, {1, 0, 0}, a_] := Lpow[g, 1, a];
lamInt[g_, {0, 1, 0}, a_] := Lpow[g, 2, a]/2;        (* lambda_2 = lambda_1^2/2 *)
lamInt[g_, {0, 0, 1}, a_] := LgTop[g, a] /; g == 3;
lamInt[g_, {1, 1, 0}, a_] := Lpow[g, 3, a]/2;        (* lambda_1 lambda_2 = lambda_1^3/2 *)
lamInt[g_, {1, 0, 1}, a_] := Lpow[g, 4, a]/8;        (* lambda_1 lambda_3 = lambda_2^2/2 = lambda_1^4/8 *)
lamInt[g_, {0, 1, 1}, a_] := Lg1[g, a] /; g == 3;
lamInt[g_, {1, 1, 1}, a_] := Lg2[g, a] /; g == 3;
lamInt[g_, m_, a_] := (Print["lamInt: unexpected monomial ", m, " in genus ", g]; Abort[]);
(* normalised Hodge correlator, polynomial in u[1], u[2], u[3] *)
hodgeNormU[g_, a_List] := hodgeNormU[g, a] = Expand[Total[Table[c[[2]] lamInt[g, c[[1]], a], {c, hodgeClass[g]}]]];
usubOf[w_List] := Table[u[a] -> Together[1/w[[a]]], {a, 3}];
(* the scalar Bernoulli series h(z)^s of the weights, coefficient list to z^KMAX *)
hCoeffs[w_List, s_, KMAX_] := Module[{p, ex},
   p[m_] := Total[w^(-m)];
   ex = -Sum[BernoulliB[2 k]/(2 k (2 k - 1)) z^(2 k - 1) p[2 k - 1], {k, 1, Ceiling[(KMAX + 1)/2]}];
   Table[Together[SeriesCoefficient[Series[Exp[s ex], {z, 0, KMAX}], j]], {j, 0, KMAX}]];

(* ============================== 3. checks ============================== *)
hodgeChecks3[] := Module[{ok = True, c, rep = {}, tuples},
   tuples[g_, nmax_] := Join @@ Table[Tuples[Range[0, 3 g - 3 + n], n], {n, 1, nmax}];
   c = {Lpow[1, 1, {0}] == 1/24, Lpow[2, 1, {3}] == 1/480, Lpow[2, 2, {2}]/2 == 7/5760, Lpow[2, 3, {1}]/2 == 1/2880, Lpow[2, 3, {}]/2 == 1/5760,
     LpowRaw[3, 6, {}] == 1/90720, LgTop[3, {4}] == 31/967680, Lg1[3, {2}] == 1/120960, Lg2[3, {}] == 1/1451520};
   AppendTo[rep, "known values (lambda_1 on M_{1,1}; lambda_1 psi^3, lambda_2 psi^2, lambda_1 lambda_2 psi on M_{2,1}; lambda_1 lambda_2 on M_2; lambda_1^6 on M_3 by six-fold Mumford = 1/90720; lambda_3 psi^4 on M_{3,1}; lambda_3 lambda_2 psi^2 on M_{3,1}; lambda_1 lambda_2 lambda_3 on M_3): " <> ToString[c]]; ok = ok && And @@ c;
   (* Mumford recursion versus the closed forms and the vanishing relations *)
   c = And @@ Flatten[{Table[LpowRaw[1, 2, a] == 0, {a, Select[tuples[1, 4], Total[#] == Length[#] - 2 &]}],
       Table[Lpow[1, 1, a] == LgTop[1, a], {a, Select[tuples[1, 4], Total[#] == Length[#] - 1 &]}],
       Table[Lpow[2, 2, a]/2 == LgTop[2, a], {a, Select[tuples[2, 4], Total[#] == Length[#] + 1 &]}],
       Table[Lpow[2, 3, a]/2 == Lg1[2, a], {a, Select[tuples[2, 4], Total[#] == Length[#] &]}],
       Table[LpowRaw[2, 4, a] == 0, {a, Select[tuples[2, 3], Total[#] == Length[#] - 1 &]}],
       Table[Lpow[3, 5, a] == 16 Lg1[3, a], {a, Select[tuples[3, 3], Total[#] == Length[#] + 1 &]}],
       Table[Lpow[3, 6, a] == 16 Lg2[3, a], {a, Select[tuples[3, 3], Total[#] == Length[#] &]}]}];
   AppendTo[rep, "Mumford recursion: lambda_1^2 = 0 (g=1), lambda_1 = lambda_g formula (g=1), lambda_1^2 = 2 lambda_2 and lambda_1^3 = 2 lambda_1 lambda_2 (g=2), lambda_1^4 = 0 (g=2), lambda_1^5 = 16 lambda_2 lambda_3 and lambda_1^6 = 16 lambda_1 lambda_2 lambda_3 with psi insertions (g=3): " <> ToString[c]]; ok = ok && c;
   (* string and dilaton equations of the Hodge correlators *)
   c = And @@ Flatten[Table[{Together[hodgeNormU[g, Append[a, 0]] - Sum[If[a[[i]] >= 1, hodgeNormU[g, ReplacePart[a, i -> a[[i]] - 1]], 0], {i, Length[a]}]] === 0,
        Together[hodgeNormU[g, Append[a, 1]] - (2 g - 2 + Length[a]) hodgeNormU[g, a]] === 0},
      {g, 1, 3}, {a, tuples[g, If[g == 3, 2, 3]]}]];
   AppendTo[rep, "string and dilaton equations for the Hodge correlators (genus 1-2: n <= 3, genus 3: n <= 2): " <> ToString[c]]; ok = ok && c;
   AppendTo[rep, "normalised Hodge class, genus 3: " <> ToString[InputForm[hodgeClass[3]]]];
   {ok, rep}];

(* R-action of the scalar h on the point theory in the PPZ conventions, genus gg, summed over the stable graphs with point (WK)
   vertices; equals hodgeNormU[gg, {}] (normalised Hodge class on M_g). Exact for the given (symbolic or rational) weights. *)
ppzScalar[gg_, w_List] := Module[{KMAX = 3 gg - 3, A, Bc, E1, Dil1, graphs, tot = 0, hp, vfac, efac},
   A = hCoeffs[w, -1, KMAX];
   Bc[i_, j_] := If[i == 0 && j == 0, 0, -A[[i + 1]] A[[j + 1]]];
   E1[m_, n_] := E1[m, n] = Together[Sum[(-1)^k Bc[m + 1 + k, n - k], {k, 0, n}]];
   Dil1[d_] := -A[[d]];
   graphs = prepareGraph /@ graphsWithLegs[gg, 0];
   Do[Do[hp = ConstantArray[0, Length[gr["H"]]];
     Do[Do[hp[[gr["hsAt"][[v, r]]]] = dec[[v, 1, r]], {r, Length[gr["hsAt"][[v]]]}], {v, gr["V"]}];
     vfac = Product[dec[[v, 4]] Product[Dil1[d], {d, dec[[v, 3]]}], {v, gr["V"]}];
     efac = Product[E1[hp[[2 e - 1]], hp[[2 e]]], {e, Length[gr["edges"]]}];
     tot += vfac efac/gr["aut"], {dec, gr["decs"]}], {gr, graphs}];
   Together[tot]];

(* ============================== 4. decorated graphs with Hodge vertices ============================== *)
vertexOptionsH[g_, hs_, m_] := Module[{res = {}, dim, nd, dils, hp, lp, budget, wgt},
   Do[dim = 3 g - 3 + hs + m + nd;
    If[dim < 2 nd, Continue[]];
    Do[budget = dim - Total[dils];
     Do[Do[
       wgt = hodgeNormU[g, Join[hp, lp, dils]]/(m! nd!);
       If[wgt =!= 0, AppendTo[res, {hp, lp, dils, wgt}]],
       {lp, Select[Tuples[Range[0, budget - Total[hp]], m], Total[#] <= budget - Total[hp] &]}],
      {hp, Select[Tuples[Range[0, budget], hs], Total[#] <= budget &]}],
     {dils, If[nd == 0, {{}}, Flatten[Table[Permutations /@ Select[IntegerPartitions[s, {nd}], Min[#] >= 2 &], {s, 2 nd, dim}], 2]]}],
    {nd, 0, Max[0, 3 g - 3 + hs + m]}];
   res];
prepareGraphH[{gv_, mv_, A_}] := Module[{V = Length[gv], ed, H, hv, hsAt, opts, decs, aut},
   ed = edgesOf[A]; H = Flatten[Table[{{e, 1}, {e, 2}}, {e, Length[ed]}], 1];
   hv = Table[ed[[h[[1]], h[[2]]]], {h, H}];
   hsAt = Table[Flatten[Position[hv, v]], {v, V}];
   opts = Table[vertexOptionsH[gv[[v]], Length[hsAt[[v]]], mv[[v]]], {v, V}];
   decs = Tuples[opts];
   aut = autCount[gv, mv, A];
   <|"g" -> gv, "m" -> mv, "A" -> A, "edges" -> ed, "H" -> H, "hsAt" -> hsAt, "decs" -> decs, "aut" -> aut, "V" -> V|>];
numGraph[gr_, sub_] := <|gr, "decs" -> (gr["decs"] /. sub)|>;   (* u[a] -> numbers *)

(* ============================== 5. modified R-matrix rows and sampling ============================== *)
rowsTimes[hc_List, KMAX_] := Table[Table[Together[Sum[PPsym[a][[j + 1]] hc[[k - j + 1]], {j, 0, k}]], {k, 0, KMAX}], {a, 0, 2}];
evalPointRows[t0_, x0_, KMAX_, rows_] := Block[{PPsym}, PPsym[0] = rows[[1]]; PPsym[1] = rows[[2]]; PPsym[2] = rows[[3]]; evalPoint[t0, x0, KMAX]];
computeSampleH[t0_, x0_, gg_, OTx_, rows_, usub_] := Module[{tt = AbsoluteTime[], data, wmax = 3 gg - 3, vecs, imax = 0, bad = {}},
   data = If[rows === None, evalPoint[t0, x0, KMAXcur], evalPointRows[t0, x0, KMAXcur, rows]];
   vecs = Table[Module[{val = graphValueVE[numGraph[OTx["graphs"][id], usub], data], rules, asc},
      rules = CoefficientRules[val, {C1i, C1p}]; asc = Association[rules];
      If[rules =!= {}, imax = Max[imax, Max[Abs[Im[rules[[All, 2]]]]]]];
      If[Select[rules, Abs[#[[2]]] > 10^-(PREC - 80) && (#[[1, 1]] != #[[1, 2]] || #[[1, 1]] > wmax) &] =!= {}, AppendTo[bad, id]];
      Table[ratl[Lookup[asc, Key[{k, k}], 0]/t0^(3 k)], {k, 0, wmax}]], {id, OTx["ids"]}];
   {{t0, x0}, vecs, {"imag" -> N[imax, 5], "bad" -> bad, "seconds" -> Round[AbsoluteTime[] - tt, 0.1]}}];
refVector[Fref_, t0_, x0_, gg_] := Table[Together[Coefficient[Fref, XX, k] /. {y -> t0^3, e1s -> x0, e2s -> 1}], {k, 0, 3 gg - 3}];

(* ============================== 6. reconstruction with (eps1 eps2)^a (eps1+eps2)^b denominators ============================== *)
AMAX = 10; BMAX = 16;
xsamplesH[nx_] := Table[Prime[r + 1], {r, 1, nx}];
xtry[vals_, xs_, gg_, a_, b_] := Module[{NX = Length[xs], deg = 2 gg - 2 + 2 a + b, pts, pol, ok},
   If[deg + 2 > NX, Return[$Failed]];
   pts = Table[{xs[[s]], xs[[s]]^a (1 + xs[[s]])^b vals[[s]]}, {s, deg + 1}];
   pol = InterpolatingPolynomial[pts, xv];
   ok = And @@ Table[(pol /. xv -> xs[[s]]) == xs[[s]]^a (1 + xs[[s]])^b vals[[s]], {s, deg + 2, NX}];
   If[ok, {a, b, CoefficientList[Expand[pol], xv, deg + 1]}, $Failed]];
(* search (a, b) by increasing a + b; `start` = {a0, b0} is tried first *)
xfitH[vals_, xs_, gg_, start_: None] := Module[{res = $Failed, a, r},
   If[Union[vals] === {0}, Return[{0, 0, {0}}]];
   If[start =!= None, r = xtry[vals, xs, gg, start[[1]], start[[2]]]; If[r =!= $Failed, Return[r]]];
   Do[Do[a = tot - b; If[a < 0 || a > AMAX, Continue[]];
     r = xtry[vals, xs, gg, a, b]; If[r =!= $Failed, res = r; Return[]], {b, 0, Min[tot, BMAX]}];
    If[res =!= $Failed, Return[]], {tot, 0, AMAX + BMAX}];
   res];
reconstructOrbitH[samp_, ts_, xs_, gg_] := Module[{term = 0, wmax = 3 gg - 3, xf, x1, amax, bmax, deg, seqs, yf, ok = True, info = {}},
   Do[x1 = xfitH[Table[samp[ts[[1]]][x0][[k + 1]], {x0, xs}], xs, gg];
    If[x1 === $Failed, ok = False; AppendTo[info, {k, "x"}]; Continue[]];
    xf = Prepend[Table[xfitH[Table[samp[t0][x0][[k + 1]], {x0, xs}], xs, gg, x1[[1 ;; 2]]], {t0, Rest[ts]}], x1];
    If[MemberQ[xf, $Failed], ok = False; AppendTo[info, {k, "x"}]; Continue[]];
    amax = Max[xf[[All, 1]]]; bmax = Max[xf[[All, 2]]]; deg = 2 gg - 2 + 2 amax + bmax;
    seqs = Table[Table[Module[{a = xf[[s, 1]], b = xf[[s, 2]], cl = xf[[s, 3]]},
         Coefficient[Expand[xv^(amax - a) (1 + xv)^(bmax - b) Sum[cl[[i + 1]] xv^i, {i, 0, Length[cl] - 1}]], xv, jj]], {s, Length[ts]}], {jj, 0, deg}];
    If[Union[Flatten[seqs]] === {0}, Continue[]];
    yf = yfitC[seqs, ts, gg];
    If[yf === $Failed, ok = False; AppendTo[info, {k, "y"}]; Continue[]];
    term += XX^k Sum[yf[[3, jj + 1]] e1s^jj e2s^(deg - jj), {jj, 0, deg}]/((e1s e2s)^amax (e1s + e2s)^bmax), {k, 0, wmax}];
   If[ok, term, {$Failed, info}]];

(* ::Subsection:: Genus-general Hodge evaluator used by the cluster package for genus 4 (hv_lib4.m): every lambda class as a polynomial in lambda_1 and ch_3, integrals by Mumford's formula *)

(* hv_lib4.m -- genus-general Hodge-vertex machinery for the Givental-Teleman graph sum of K_P2 x C^2 on T_F, valid in genus <= 4
   (no closed forms needed): every lambda class is a polynomial in lambda_1 and ch_3(E) (Newton's identities with ch_2 = ch_4 = 0:
   lambda_2 = lambda_1^2/2, lambda_3 = lambda_1^3/6 + 2 ch_3, lambda_4 = lambda_1^4/24 + 2 lambda_1 ch_3), and
   HI[g, k1, k3, a] = Int_{M_{g,n}} lambda_1^k1 ch_3^k3 psi^a is computed by Mumford's formula for ch_1 = lambda_1 and ch_3:
     ch_k = B_{k+1}/(k+1)! ( kappa_k - Sum_i psi_i^k + 1/2 Sum_boundary iota_* Sum_{j=0}^{k-1} (-1)^j psi'^j psi''^(k-1-j) ),
   kappa_k = pi_*(psi_{n+1}^(k+1)) so that <kappa_k X> = <tau_{k+1} X> for pulled-back X; ch_k restricts additively to the
   components of the boundary. Vertex = normalised Hodge class prod_a Lambda_g^vee(w_a)/w_a^g, u[a] = 1/w_a.
   The rest (decorations, R' = R_Z h^(-1) rows, sampling, fits) as in hv_lib3.m. Built on lp2_lib.m, graphs.m, orbits.m. ASCII only. *)

e1s = \[Epsilon][1]; e2s = \[Epsilon][2];

(* ============================== 1. Hodge integrals by Mumford's formula ============================== *)
(* sum over the subsets S of the marked points of f[a_S, a_S^c], organised by sub-multisets of the sorted list a *)
subsetSum[a_List, f_] := Module[{tal = Tally[a], tot = 0},
   Do[tot += (Times @@ MapThread[Binomial, {tal[[All, 2]], sel}]) f[Flatten[MapThread[ConstantArray, {tal[[All, 1]], sel}]], Flatten[MapThread[ConstantArray, {tal[[All, 1]], tal[[All, 2]] - sel}]]],
    {sel, Tuples[Range[0, #] & /@ tal[[All, 2]]]}];
   tot];
(* lambda monomials of degree > dmax[g] vanish on every M_{g,n}: they are pulled back from M_{g,1} (g = 1) or M_g (g >= 2) *)
dmax[g_] := If[g <= 0, 0, If[g == 1, 1, Min[g (g + 1)/2, 3 g - 3]]];
(* the Mumford recursion; F is the function called recursively (HI with the vanishing cut-off, HIraw without) *)
hiBody[g_, k1_, k3_, a_List, F_] := Module[{n = Length[a], tot},
   Which[
    g < 0 || k1 < 0 || k3 < 0, 0,
    2 g - 2 + n <= 0, 0,
    k1 + 3 k3 + Total[a] != 3 g - 3 + n, 0,
    g == 0 && (k1 > 0 || k3 > 0), 0,
    k1 == 0 && k3 == 0, wk[g, a],
    k1 > 0,
    tot = F[g, k1 - 1, k3, Append[a, 2]] - Sum[F[g, k1 - 1, k3, ReplacePart[a, i -> a[[i]] + 1]], {i, n}] + F[g - 1, k1 - 1, k3, Join[a, {0, 0}]]/2;
    tot += 1/2 subsetSum[a, Function[{aS, aSc}, Sum[Binomial[k1 - 1, j1] Binomial[k3, j3] F[h, j1, j3, Append[aS, 0]] F[g - h, k1 - 1 - j1, k3 - j3, Append[aSc, 0]], {h, 0, g}, {j1, 0, k1 - 1}, {j3, 0, k3}]]];
    tot/12,
    True,
    tot = F[g, 0, k3 - 1, Append[a, 4]] - Sum[F[g, 0, k3 - 1, ReplacePart[a, i -> a[[i]] + 3]], {i, n}] + 1/2 Sum[(-1)^j F[g - 1, 0, k3 - 1, Join[a, {j, 2 - j}]], {j, 0, 2}];
    tot += 1/2 subsetSum[a, Function[{aS, aSc}, Sum[(-1)^j Binomial[k3 - 1, j3] F[h, 0, j3, Append[aS, j]] F[g - h, 0, k3 - 1 - j3, Append[aSc, 2 - j]], {h, 0, g}, {j, 0, 2}, {j3, 0, k3 - 1}]]];
    BernoulliB[4]/4! tot]];
HI[g_, k1_, k3_, a_List] /; a =!= Sort[a] := HI[g, k1, k3, Sort[a]];
HI[g_, k1_, k3_, a_List] := HI[g, k1, k3, a] = If[k1 + 3 k3 > dmax[g], 0, hiBody[g, k1, k3, a, HI]];
HIraw[g_, k1_, k3_, a_List] /; a =!= Sort[a] := HIraw[g, k1, k3, Sort[a]];
HIraw[g_, k1_, k3_, a_List] := HIraw[g, k1, k3, a] = hiBody[g, k1, k3, a, HIraw];

(* ============================== 2. the normalised Hodge class ============================== *)
lamSub = {lam[1] -> L1, lam[2] -> L1^2/2, lam[3] -> L1^3/6 + 2 C3, lam[4] -> L1^4/24 + 2 L1 C3};
(* {{k1, k3}, coefficient in u[1], u[2], u[3]} of lambda_1^k1 ch_3^k3 in prod_a Lambda_g^vee(w_a)/w_a^g *)
hodgeClass4[g_] /; g > 4 := (Print["hodgeClass4: genus ", g, " needs lambda_5 and the ch_5 Mumford branch (not implemented)"]; Abort[]);
hodgeClass4[g_] := hodgeClass4[g] = Module[{cl},
   cl = Expand[Product[Sum[(-1)^j lam[j] u[a]^j, {j, 0, g}], {a, 3}] /. lam[0] -> 1 /. lamSub];
   Select[{#[[1]], Factor[#[[2]]]} & /@ CoefficientRules[cl, {L1, C3}], #[[1, 1]] + 3 #[[1, 2]] <= dmax[g] &]];
hodgeNormU4[g_, a_List] := hodgeNormU4[g, a] = Expand[Total[Table[c[[2]] HI[g, c[[1, 1]], c[[1, 2]], a], {c, hodgeClass4[g]}]]];
usubOf[w_List] := Table[u[a] -> Together[1/w[[a]]], {a, 3}];
hCoeffs[w_List, s_, KMAX_] := Module[{p, ex},
   p[m_] := Total[w^(-m)];
   ex = -Sum[BernoulliB[2 k]/(2 k (2 k - 1)) z^(2 k - 1) p[2 k - 1], {k, 1, Ceiling[(KMAX + 1)/2]}];
   Table[Together[SeriesCoefficient[Series[Exp[s ex], {z, 0, KMAX}], j]], {j, 0, KMAX}]];

(* ============================== 3. checks ============================== *)
multin[a_List] := Total[a]!/Times @@ (a!);
bg[g_] := Abs[BernoulliB[2 g]];
LgTop[g_, a_List] := If[g >= 1 && 2 g - 2 + Length[a] > 0 && Total[a] == 2 g - 3 + Length[a], multin[a] (2^(2 g - 1) - 1)/2^(2 g - 1) bg[g]/(2 g)!, 0];
Lg1[g_, a_List] := Module[{n = Length[a], z, rest},
   If[g < 1 || 2 g - 2 + n <= 0 || Total[a] != g - 2 + n, Return[0]];
   z = Position[a, 0, {1}, Heads -> False];
   If[z =!= {} && 2 g - 2 + n - 1 > 0, rest = Delete[a, z[[1]]];
    Return[Sum[If[rest[[j]] >= 1, Lg1[g, ReplacePart[rest, j -> rest[[j]] - 1]], 0], {j, Length[rest]}]]];
   (2 g - 3 + n)! bg[g]/(2^(2 g - 1) (2 g)! Times @@ (df[2 # - 1] & /@ a))];
(* lambda_g, lambda_g lambda_{g-1} and lambda_g lambda_{g-1} lambda_{g-2} as polynomials in lambda_1, ch_3, integrated by HI *)
lamPoly[g_, mono_List] := Expand[Product[lam[j], {j, mono}] /. lamSub];
lamIntHI[g_, mono_List, a_List] := Module[{p = lamPoly[g, mono]}, Total[Table[c[[2]] HI[g, c[[1, 1]], c[[1, 2]], a], {c, CoefficientRules[p, {L1, C3}]}]]];
hodgeChecks4[] := Module[{ok = True, c, rep = {}, tuples},
   tuples[g_, nmax_] := Join @@ Table[Tuples[Range[0, 3 g - 3 + n], n], {n, 1, nmax}];
   (* known values *)
   c = {HI[1, 1, 0, {0}] == 1/24, HI[2, 1, 0, {3}] == 1/480, HI[2, 2, 0, {2}]/2 == 7/5760, HI[2, 3, 0, {1}]/2 == 1/2880, HI[2, 3, 0, {}]/2 == 1/5760,
     HI[2, 0, 1, {}] == -1/34560, HIraw[3, 6, 0, {}] == 1/90720, lamIntHI[3, {3}, {4}] == 31/967680, lamIntHI[3, {3, 2}, {2}] == 1/120960, lamIntHI[3, {1, 2, 3}, {}] == 1/1451520,
     HIraw[4, 9, 0, {}] == 1/113400, lamIntHI[4, {4}, {6}] == 127/154828800, lamIntHI[4, {4, 3}, {3}] == 1/3225600, lamIntHI[4, {4, 3}, {2, 2}] == 1/276480, lamIntHI[4, {4, 3, 2}, {}] == 1/87091200};
   AppendTo[rep, "known values: genus 1-2 (lambda_1; lambda_1 psi^3, lambda_2 psi^2, lambda_1 lambda_2 psi, lambda_1 lambda_2; ch_3 on M_2 = -1/34560), genus 3 (lambda_1^6 = 1/90720 raw; lambda_3 psi^4; lambda_3 lambda_2 psi^2; lambda_1 lambda_2 lambda_3), genus 4 (lambda_1^9 = 1/113400 raw; lambda_4 psi^6 = 127/154828800; lambda_4 lambda_3 psi^3 = 1/3225600; lambda_4 lambda_3 psi^2 psi^2 = 1/276480 [GP]; lambda_4 lambda_3 lambda_2 = 1/87091200): " <> ToString[c]]; ok = ok && And @@ c;
   (* vanishing relations without the cut-off: ch_3 = 0 in genus 1, lambda_1^4 = 0 in genus 2, lambda_1^7 = 0 in genus 3, degree-10 monomials in genus 4 *)
   c = And @@ Flatten[{Table[HIraw[1, 0, 1, a] == 0, {a, Select[tuples[1, 4], Total[#] == Length[#] - 3 &]}],
       Table[HIraw[2, 4, 0, a] == 0, {a, Select[tuples[2, 3], Total[#] == Length[#] - 1 &]}],
       Table[HIraw[3, 7, 0, a] == 0, {a, Select[tuples[3, 2], Total[#] == Length[#] - 1 &]}],
       HIraw[3, 4, 1, {0}] == 0, HIraw[4, 10, 0, {0}] == 0, HIraw[4, 7, 1, {0}] == 0, HIraw[4, 4, 2, {0}] == 0, HIraw[4, 1, 3, {0}] == 0}];
   AppendTo[rep, "vanishing above the top degree of the lambda ring, computed without the cut-off: " <> ToString[c]]; ok = ok && c;
   (* Mumford's relation c(E) c(E^vee) = 1 is NOT built into lamSub (which only uses ch_2 = ch_4 = 0): lambda_g^2 = 0 and
      lambda_{g-1}^2 = 2 lambda_g lambda_{g-2} are independent tests of the ch_3 sector (in genus 4 they pin Int_{M_4} ch_3^3) *)
   c = {lamIntHI[4, {3, 3, 3}, {}] == 1/43545600, lamIntHI[3, {2, 2, 2}, {}] == 1/725760,
     And @@ Flatten[Table[If[2 g + Total[a] == 3 g - 3 + Length[a], lamIntHI[g, {g, g}, a] == 0, True], {g, 2, 4}, {a, tuples[g, 2]}]],
     And @@ Flatten[Table[If[6 + Total[a] == 9 + Length[a], lamIntHI[4, {3, 3}, a] == 2 lamIntHI[4, {4, 2}, a], True], {a, tuples[4, 3]}]],
     And @@ Flatten[Table[If[4 + Total[a] == 6 + Length[a], lamIntHI[3, {2, 2}, a] == 2 lamIntHI[3, {3, 1}, a], True], {a, tuples[3, 3]}]]};
   AppendTo[rep, "Mumford relations not encoded in lamSub (lambda_3^3 = 2 lambda_4 lambda_3 lambda_2 = 1/43545600 on M_4, lambda_2^3 = 1/725760 on M_3, lambda_g^2 = 0 (g = 2-4), lambda_3^2 = 2 lambda_4 lambda_2 on M_{4,n}, lambda_2^2 = 2 lambda_3 lambda_1 on M_{3,n}): " <> ToString[c]]; ok = ok && And @@ c;
   (* lambda_g and lambda_g lambda_{g-1} closed forms versus the Mumford recursion with psi insertions *)
   c = And @@ Flatten[Table[{If[Total[a] == 2 g - 3 + Length[a], lamIntHI[g, {g}, a] == LgTop[g, a], True], If[Total[a] == g - 2 + Length[a] && g >= 2, lamIntHI[g, {g, g - 1}, a] == Lg1[g, a], True]},
       {g, 2, 4}, {a, tuples[g, If[g == 4, 2, 3]]}]];
   AppendTo[rep, "lambda_g formula and lambda_g lambda_{g-1} closed form (zeros by string) versus the recursion, genus 2-4 (genus 4: n <= 2): " <> ToString[c]]; ok = ok && c;
   (* string and dilaton equations *)
   c = And @@ Flatten[Table[{Together[hodgeNormU4[g, Append[a, 0]] - Sum[If[a[[i]] >= 1, hodgeNormU4[g, ReplacePart[a, i -> a[[i]] - 1]], 0], {i, Length[a]}]] === 0,
        Together[hodgeNormU4[g, Append[a, 1]] - (2 g - 2 + Length[a]) hodgeNormU4[g, a]] === 0},
      {g, 1, 4}, {a, tuples[g, If[g >= 3, 2, 3]]}]];
   AppendTo[rep, "string and dilaton equations of the Hodge correlators, genus 1-4: " <> ToString[c]]; ok = ok && c;
   {ok, rep}];
(* R-action of the scalar h on the point theory in the PPZ conventions, genus gg, from prepared point graphs (prepareGraph) *)
ppzScalar[gg_, w_List, graphs_: Automatic] := Module[{KMAX = 3 gg - 3, A, Bc, E1, Dil1, grs, tot = 0, hp, vfac, efac},
   A = hCoeffs[w, -1, KMAX];
   Bc[i_, j_] := If[i == 0 && j == 0, 0, -A[[i + 1]] A[[j + 1]]];
   E1[m_, n_] := E1[m, n] = Together[Sum[(-1)^k Bc[m + 1 + k, n - k], {k, 0, n}]];
   Dil1[d_] := -A[[d]];
   grs = If[graphs === Automatic, prepareGraph /@ graphsWithLegs[gg, 0], graphs];
   Do[Do[hp = ConstantArray[0, Length[gr["H"]]];
     Do[Do[hp[[gr["hsAt"][[v, r]]]] = dec[[v, 1, r]], {r, Length[gr["hsAt"][[v]]]}], {v, gr["V"]}];
     vfac = Product[dec[[v, 4]] Product[Dil1[d], {d, dec[[v, 3]]}], {v, gr["V"]}];
     efac = Product[E1[hp[[2 e - 1]], hp[[2 e]]], {e, Length[gr["edges"]]}];
     tot += vfac efac/gr["aut"], {dec, gr["decs"]}], {gr, grs}];
   Together[tot]];

(* ============================== 4. decorated graphs with Hodge vertices ============================== *)
vertexOptionsH4[g_, hs_, m_] := Module[{res = {}, dim, nd, dils, hp, lp, budget, wgt},
   Do[dim = 3 g - 3 + hs + m + nd;
    If[dim < 2 nd, Continue[]];
    Do[budget = dim - Total[dils];
     Do[Do[
       wgt = hodgeNormU4[g, Join[hp, lp, dils]]/(m! nd!);
       If[wgt =!= 0, AppendTo[res, {hp, lp, dils, wgt}]],
       {lp, Select[Tuples[Range[0, budget - Total[hp]], m], Total[#] <= budget - Total[hp] &]}],
      {hp, Select[Tuples[Range[0, budget], hs], Total[#] <= budget &]}],
     {dils, If[nd == 0, {{}}, Flatten[Table[Permutations /@ Select[IntegerPartitions[s, {nd}], Min[#] >= 2 &], {s, 2 nd, dim}], 2]]}],
    {nd, 0, Max[0, 3 g - 3 + hs + m]}];
   res];
prepareGraphH4[{gv_, mv_, A_}] := Module[{V = Length[gv], ed, H, hv, hsAt, opts, decs, aut},
   ed = edgesOf[A]; H = Flatten[Table[{{e, 1}, {e, 2}}, {e, Length[ed]}], 1];
   hv = Table[ed[[h[[1]], h[[2]]]], {h, H}];
   hsAt = Table[Flatten[Position[hv, v]], {v, V}];
   opts = Table[vertexOptionsH4[gv[[v]], Length[hsAt[[v]]], mv[[v]]], {v, V}];
   decs = Tuples[opts];
   aut = autCount[gv, mv, A];
   <|"g" -> gv, "m" -> mv, "A" -> A, "edges" -> ed, "H" -> H, "hsAt" -> hsAt, "decs" -> decs, "aut" -> aut, "V" -> V|>];
numGraph[gr_, sub_] := <|gr, "decs" -> (gr["decs"] /. sub)|>;

(* ============================== 5. R' rows, sampling, reconstruction ============================== *)
rowsTimes[hc_List, KMAX_] := Table[Table[Together[Sum[PPsym[a][[j + 1]] hc[[k - j + 1]], {j, 0, k}]], {k, 0, KMAX}], {a, 0, 2}];
evalPointRows[t0_, x0_, KMAX_, rows_] := Block[{PPsym}, PPsym[0] = rows[[1]]; PPsym[1] = rows[[2]]; PPsym[2] = rows[[3]]; evalPoint[t0, x0, KMAX]];
computeSampleH[t0_, x0_, gg_, OTx_, rows_, usub_] := Module[{tt = AbsoluteTime[], data, wmax = 3 gg - 3, vecs, imax = 0, bad = {}},
   data = If[rows === None, evalPoint[t0, x0, KMAXcur], evalPointRows[t0, x0, KMAXcur, rows]];
   vecs = Table[Module[{val = graphValueVE[numGraph[OTx["graphs"][id], usub], data], rules, asc},
      rules = CoefficientRules[val, {C1i, C1p}]; asc = Association[rules];
      If[rules =!= {}, imax = Max[imax, Max[Abs[Im[rules[[All, 2]]]]]]];
      If[Select[rules, Abs[#[[2]]] > 10^-(PREC - 80) && (#[[1, 1]] != #[[1, 2]] || #[[1, 1]] > wmax) &] =!= {}, AppendTo[bad, id]];
      Table[ratl[Lookup[asc, Key[{k, k}], 0]/t0^(3 k)], {k, 0, wmax}]], {id, OTx["ids"]}];
   {{t0, x0}, vecs, {"imag" -> N[imax, 5], "bad" -> bad, "seconds" -> Round[AbsoluteTime[] - tt, 0.1]}}];
refVector[Fref_, t0_, x0_, gg_] := Table[Together[Coefficient[Fref, XX, k] /. {y -> t0^3, e1s -> x0, e2s -> 1}], {k, 0, 3 gg - 3}];
(* eps-ratios: mode-B terms have denominators (eps1 eps2)^a (eps1+eps2)^b with a, b <= 1 at genus 2 and 3; allow up to AMAX, BMAX *)
AMAX = 5; BMAX = 5;
nxFor[gg_] := Switch[gg, 2, 34, 3, 30, 4, 20];
xsamplesH[nx_] := Table[Prime[r + 1], {r, 1, nx}];
pointListH[gg_] := Flatten[Table[{t0, x0}, {t0, tsamplesFor[gg]}, {x0, xsamplesH[nxFor[gg]]}], 1];
xtry[vals_, xs_, gg_, a_, b_] := Module[{NX = Length[xs], deg = 2 gg - 2 + 2 a + b, pts, pol, ok},
   If[deg + 2 > NX, Return[$Failed]];
   pts = Table[{xs[[s]], xs[[s]]^a (1 + xs[[s]])^b vals[[s]]}, {s, deg + 1}];
   pol = InterpolatingPolynomial[pts, xv];
   ok = And @@ Table[(pol /. xv -> xs[[s]]) == xs[[s]]^a (1 + xs[[s]])^b vals[[s]], {s, deg + 2, NX}];
   If[ok, {a, b, CoefficientList[Expand[pol], xv, deg + 1]}, $Failed]];
xfitH[vals_, xs_, gg_, start_: None] := Module[{res = $Failed, a, r},
   If[Union[vals] === {0}, Return[{0, 0, {0}}]];
   If[start =!= None, r = xtry[vals, xs, gg, start[[1]], start[[2]]]; If[r =!= $Failed, Return[r]]];
   Do[Do[a = tot - b; If[a < 0 || a > AMAX, Continue[]];
     r = xtry[vals, xs, gg, a, b]; If[r =!= $Failed, res = r; Return[]], {b, 0, Min[tot, BMAX]}];
    If[res =!= $Failed, Return[]], {tot, 0, AMAX + BMAX}];
   res];
reconstructOrbitH[samp_, ts_, xs_, gg_] := Module[{term = 0, wmax = 3 gg - 3, xf, x1, amax, bmax, deg, seqs, yf, ok = True, info = {}},
   Do[x1 = xfitH[Table[samp[ts[[1]]][x0][[k + 1]], {x0, xs}], xs, gg];
    If[x1 === $Failed, ok = False; AppendTo[info, {k, "x"}]; Continue[]];
    xf = Prepend[Table[xfitH[Table[samp[t0][x0][[k + 1]], {x0, xs}], xs, gg, x1[[1 ;; 2]]], {t0, Rest[ts]}], x1];
    If[MemberQ[xf, $Failed], ok = False; AppendTo[info, {k, "x"}]; Continue[]];
    amax = Max[xf[[All, 1]]]; bmax = Max[xf[[All, 2]]]; deg = 2 gg - 2 + 2 amax + bmax;
    seqs = Table[Table[Module[{a = xf[[s, 1]], b = xf[[s, 2]], cl = xf[[s, 3]]},
         Coefficient[Expand[xv^(amax - a) (1 + xv)^(bmax - b) Sum[cl[[i + 1]] xv^i, {i, 0, Length[cl] - 1}]], xv, jj]], {s, Length[ts]}], {jj, 0, deg}];
    If[Union[Flatten[seqs]] === {0}, Continue[]];
    yf = yfitC[seqs, ts, gg];
    If[yf === $Failed, ok = False; AppendTo[info, {k, "y"}]; Continue[]];
    term += XX^k Sum[yf[[3, jj + 1]] e1s^jj e2s^(deg - jj), {jj, 0, deg}]/((e1s e2s)^amax (e1s + e2s)^bmax), {k, 0, wmax}];
   If[ok, term, {$Failed, info}]];

(* ::Section:: Setup *)

(* data directory: this notebook's folder if it holds the sidecars (cPP_z6.m, StableGraph-2/3.bin, ref/F20_F30_y.m, hv_orbits_g3.mx,
   hv3_samples_B.mx, hodge_vertex_genus2_results.m), otherwise its subfolder hodge_vertex/ *)
datadir = If[FileExistsQ[FileNameJoin[{NotebookDirectory[], "cPP_z6.m"}]], NotebookDirectory[], FileNameJoin[{NotebookDirectory[], "hodge_vertex"}]];
SetDirectory[datadir]; Print["data directory: ", Directory[]];
e1 = \[Epsilon][1]; e2 = \[Epsilon][2];
wB = {e1 + e2, -e1, -e2};          (* e_T weights of the normal bundle O(-3) + O + O of the fixed P^2 *)
usubB = usubOf[wB];
{Fref20, Fref30} = Get["ref/F20_F30_y.m"];   (* F_{2,0}, F_{3,0} of localp2CY5_genus3.nb in y, XX, eps *)

(* reconstruction of every decorated graph from a list of sample records {{t0, x0}, vectors, diagnostics, ...} *)
reconstructAll[recs_, OTx_, ts_, xs_, gg_] := Module[{samp, sv},
   Do[samp[r[[1, 1]]][r[[1, 2]]] = r[[2]], {r, recs}];
   Association[Table[Module[{o = oo}, Clear[sv]; Do[sv[t0][x0] = samp[t0][x0][[o]], {t0, ts}, {x0, xs}];
      OTx["ids"][[o]] -> reconstructOrbitH[sv, ts, xs, gg]], {oo, Length[OTx["ids"]]}]]];
graphSums[terms_, OTx_, ngraphs_] := Table[Together[Total[Table[terms[id], {id, Select[OTx["ids"], #[[1]] == k &]}]]], {k, ngraphs}];
polarOf[term_, gg_] := polarOfTerm[term, gg] /. SeriesCoefficient[a_, _] /; FreeQ[a, tcv] :> 0;
sampleAll[pts_, gg_, OTx_, rows_] := (If[$KernelCount == 0, Quiet[LaunchKernels[]]];
   If[$KernelCount > 0, DistributeDefinitions["Global`"];
    ParallelTable[computeSampleH[p[[1]], p[[2]], gg, OTx, rows, usubB /. {e1 -> p[[2]], e2 -> 1}], {p, pts}],
    Table[computeSampleH[p[[1]], p[[2]], gg, OTx, rows, usubB /. {e1 -> p[[2]], e2 -> 1}], {p, pts}]]);

(* ::Section:: Checks of the Hodge integrals and of the scalar R-action *)

(* Hodge integrals: known values, Mumford recursion against the closed forms, string and dilaton equations *)
{hodgeOK, hodgeReport} = hodgeChecks3[]; Column[Append[hodgeReport, "all checks passed: " <> ToString[hodgeOK]]]

(* the scalar R-action of h on the point theory (PPZ conventions) reproduces the normalised Hodge class on M_2 for arbitrary weights ... *)
ppz2 = ppzScalar[2, {w1, w2, w3}]; {Factor[ppz2], Together[ppz2 - (hodgeNormU[2, {}] /. usubOf[{w1, w2, w3}])] === 0}

(* ... and on M_3 (42 stable graphs) at exact rational weights *)
SeedRandom[7]; Table[Module[{wn = Table[RandomInteger[{-9, 9}]/RandomInteger[{1, 5}], {3}]}, If[MemberQ[wn, 0], Nothing, {wn, ppzScalar[3, wn] == (hodgeNormU[3, {}] /. usubOf[wn])}]], {2}]

(* the normalised Hodge class of genus 3 for the weights wB, as coefficients of lambda_1^i lambda_2^j lambda_3^k *)
Table[{c[[1]], Factor[c[[2]] /. usubB]}, {c, hodgeClass[3]}] // TableForm

(* the scalar h_wB(z) to z^6 *)
Factor /@ hCoeffs[wB, 1, 6]

(* the genus-general evaluator (hv_lib4.m): known values up to genus 4 (lambda_1^9 on M_4 = 1/113400 by nine-fold Mumford, the
   Getzler-Pandharipande value lambda_4 lambda_3 psi^2 psi^2 = 1/276480, lambda_4 lambda_3 lambda_2 = 1/87091200, ...), vanishing above
   the top degree, closed forms versus the recursion, string/dilaton; and agreement with hv_lib3's correlators on the genus-3 decorations *)
{hodge4OK, hodge4Report} = hodgeChecks4[]; Column[Append[hodge4Report, "all checks passed: " <> ToString[hodge4OK]]]
If[ValueQ[preparedH3] || FileExistsQ["hv_orbits_g3.mx"], If[! ValueQ[preparedH3], Get["hv_orbits_g3.mx"]];
  Module[{d3 = Union[Flatten[Table[Table[{gr["g"][[v]], Join[dec[[v, 1]], dec[[v, 2]], dec[[v, 3]]]}, {v, gr["V"]}], {gr, preparedH3}, {dec, gr["decs"]}], 2]]},
   {Length[d3], And @@ Table[Together[hodgeNormU4[d[[1]], d[[2]]] - hodgeNormU[d[[1]], d[[2]]]] === 0, {d, d3}]}]]

(* the normalised genus-4 Hodge class for wB as coefficients of lambda_1^k1 ch_3^k3 *)
Table[{c[[1]], Factor[c[[2]] /. usubOf[wB]]}, {c, hodgeClass4[4]}] // TableForm

(* ::Section:: Genus 2 (regression against the stored run) *)

(* genus 2: 7 stable graphs, 27 decorated graphs with Hodge vertices; 26 t x 34 eps-ratio sample points (a few minutes) *)
gg = 2; KMAX = kmaxFor[gg];
preparedH2 = prepareGraphH /@ graphsWithLegs[2, 0]; OTH2 = orbitTable[preparedH2];
{Length[preparedH2], Total[Length[#["decs"]] & /@ preparedH2], Length[OTH2["ids"]]}
loadRmatrix[KMAX]; rows2 = rowsTimes[hCoeffs[wB, -1, KMAX], KMAX];   (* rows of R' = R_Z h^(-1) *)
ts2 = tsamplesFor[2]; xs2 = xsamplesH[34]; pts2 = Flatten[Table[{t0, x0}, {t0, ts2}, {x0, xs2}], 1];
AbsoluteTiming[recs2 = sampleAll[pts2, 2, OTH2, rows2];]

(* the sum over the decorated graphs equals F_{2,0} exactly at every sample point *)
Union[Flatten[Table[Total[r[[2]]] - refVector[Fref20, r[[1, 1]], r[[1, 2]], 2], {r, recs2}]]]
terms2 = reconstructAll[recs2, OTH2, ts2, xs2, 2]; Select[Keys[terms2], ListQ[terms2[#]] &]   (* failed reconstructions: none expected *)
graphTerms2 = graphSums[terms2, OTH2, Length[preparedH2]]; Together[Total[graphTerms2] - Fref20] === 0

(* agreement with the stored run (hodge_vertex_genus2_results.m, mode B) *)
stored2 = Get["hodge_vertex_genus2_results.m"]["B"]["graphTerms"]; Table[Together[graphTerms2[[k]] - stored2[[k]]] === 0, {k, Length[preparedH2]}]

(* conifold polar parts {t_c^-2, t_c^-1} (eps1 = x, eps2 = 1): everything sits in the genus-2 vertex graph *)
tcSetup[8]; Grid[Prepend[Table[{k, preparedH2[[k]]["g"], preparedH2[[k]]["edges"], Factor /@ polarOf[graphTerms2[[k]], 2]}, {k, Length[preparedH2]}], {"graph", "genera", "edges", "polar part"}], Frame -> All]
{Factor /@ Total[Table[polarOf[graphTerms2[[k]], 2], {k, Length[preparedH2]}]], Factor[gapCoefficient[2] /. {e1 -> x, e2 -> 1}]}

(* ::Section:: Genus 3 *)

(* genus 3: 42 stable graphs; decorated graphs with Hodge vertices from the sidecar hv_orbits_g3.mx (or recomputed) *)
gg = 3; KMAX = kmaxFor[gg];
If[FileExistsQ["hv_orbits_g3.mx"], Get["hv_orbits_g3.mx"], preparedH3 = prepareGraphH /@ graphsWithLegs[3, 0]; OTH3 = orbitTable[preparedH3]];
{Length[preparedH3], Total[Length[#["decs"]] & /@ preparedH3], Length[OTH3["ids"]]}
loadRmatrix[KMAX]; rows3 = rowsTimes[hCoeffs[wB, -1, KMAX], KMAX];
ts3 = tsamplesFor[3]; xs3 = xsamplesH[30]; pts3 = Flatten[Table[{t0, x0}, {t0, ts3}, {x0, xs3}], 1];

(* sampling: 780 points; the sidecar hv3_samples_B.mx holds the records of the stored run (delete it to recompute; see the
   README for the timing per point) *)
If[FileExistsQ["hv3_samples_B.mx"], Get["hv3_samples_B.mx"], AbsoluteTiming[recs3 = sampleAll[pts3, 3, OTH3, rows3];]; DumpSave["hv3_samples_B.mx", recs3]];
Length[recs3]

(* the sum over the decorated graphs equals F_{3,0} exactly at every sample point *)
Union[Flatten[Table[Total[r[[2]]] - refVector[Fref30, r[[1, 1]], r[[1, 2]], 3], {r, recs3}]]]
AbsoluteTiming[terms3 = reconstructAll[recs3, OTH3, ts3, xs3, 3];]
Select[Keys[terms3], ListQ[terms3[#]] &]   (* failed reconstructions: none expected *)
graphTerms3 = graphSums[terms3, OTH3, Length[preparedH3]]; F30hodge = Together[Total[graphTerms3]];
SeedRandom[4]; Table[Together[(F30hodge - Fref30) /. {y -> RandomInteger[{1, 50}]/RandomInteger[{51, 97}], XX -> RandomInteger[{-30, 30}]/RandomInteger[{1, 17}], e1 -> RandomInteger[{1, 40}], e2 -> RandomInteger[{1, 40}]}], {5}]
TimeConstrained[Together[F30hodge - Fref30] === 0, 1800, "timeout"]

(* conifold frame: polar parts {t_c^-4, t_c^-3, t_c^-2, t_c^-1} (eps1 = x, eps2 = 1) *)
tcSetup[10]; polar3 = Table[polarOf[graphTerms3[[k]], 3], {k, Length[preparedH3]}];
{Factor /@ Total[polar3], Factor[gapCoefficient[3] /. {e1 -> x, e2 -> 1}]}   (* gap: leading = (-3)^2 3! [t^6] t^2/(4 sinh sinh), the rest 0 *)
Grid[Prepend[Table[{k, preparedH3[[k]]["g"], preparedH3[[k]]["edges"], preparedH3[[k]]["aut"], Length[Select[OTH3["ids"], #[[1]] == k &]], If[Union[polar3[[k]]] === {0}, "regular", Factor /@ polar3[[k]]]}, {k, Length[preparedH3]}], {"graph", "genera", "edges", "|Aut|", "decorated graphs", "polar part"}], Frame -> All]
polarOrbits3 = Association[Table[id -> polarOf[terms3[id], 3], {id, Keys[terms3]}]];
singular3 = Select[Keys[polarOrbits3], Union[polarOrbits3[#]] =!= {0} &]; {Length[singular3], singular3}
Put[<|"graphTerms" -> graphTerms3, "orbitTerms" -> terms3, "polar" -> polar3, "polarOrbits" -> polarOrbits3, "total" -> F30hodge, "descriptions" -> OTH3["descriptions"], "weights" -> wB|>, "hodge_vertex_genus3_results_nb.m"]

(* ::Section:: Display *)

(* the genus-3 vertex graph: decomposition of its polar part over its decorated graphs (dilaton-leaf patterns) *)
Grid[Prepend[Table[{id, OTH3["descriptions"][id]["dilatonLeaves"], Factor /@ polarOrbits3[id]}, {id, Select[Keys[terms3], #[[1]] == 1 &]}], {"id", "dilaton leaves", "polar part"}], Frame -> All]

(* which stable graphs contain singular decorated graphs (expected: exactly those with a vertex of genus >= 2), and whether the
   cancellation inside such a graph needs the full sum over its decorations (grouping by psi assignment or by dilaton pattern alone) *)
regQ[v_] := Union[v] === {0};
singGraphs = Select[Range[Length[preparedH3]], Function[k, AnyTrue[Select[Keys[polarOrbits3], #[[1]] == k &], ! regQ[polarOrbits3[#]] &]]];
{singGraphs, Table[{k, preparedH3[[k]]["g"], preparedH3[[k]]["edges"]}, {k, singGraphs}], Select[Range[Length[preparedH3]], Max[preparedH3[[#]]["g"]] >= 2 &]}
Table[Module[{ids = Select[Keys[polarOrbits3], #[[1]] == k &]}, {k, Count[GatherBy[ids, OTH3["descriptions"][#]["psiOnHalfEdges"] &], g_ /; ! regQ[Total[polarOrbits3 /@ g]]], Count[GatherBy[ids, OTH3["descriptions"][#]["dilatonLeaves"] &], g_ /; ! regQ[Total[polarOrbits3 /@ g]]]}], {k, singGraphs}]

(* the terms of a few edge graphs as polynomials in XX (E1 = eps1, E2 = eps2) *)
showTerm[term_] := Module[{tm = Together[term]}, Column[Table[Row[{"XX^", j, " : ", Factor[Coefficient[Numerator[tm], XX, j]/Denominator[tm]] /. {e1 -> E1, e2 -> E2}}], {j, 0, Exponent[Numerator[tm], XX]}]]];
Table[{k, preparedH3[[k]]["g"], preparedH3[[k]]["edges"], showTerm[graphTerms3[[k]]]}, {k, Select[Range[Length[preparedH3]], Length[preparedH3[[#]]["g"]] >= 3 &]}] // TableForm

(* ::Section:: Genus 4 (cluster results) *)

(* genus 4 is computed on the cluster with cluster_genus4_hodge/ (379 stable graphs; the same code: hv_lib4.m generalises the Hodge
   integrals to genus 4 with Int lambda_1^k1 ch_3^k3 psi^a by Mumford's formula). This section reads the results copied back
   into the data directory: hodge_stablegraphs_genus4.m (per stable graph), hodge_graphterms_genus4.m (per decorated graph). *)
g4files = {"hodge_stablegraphs_genus4.m", "hodge_graphterms_genus4.m"};
If[! AllTrue[g4files, FileExistsQ], Print["genus-4 result files not present in ", Directory[], " -- run cluster_genus4_hodge/ on the cluster and copy back ", g4files], stable4 = Get["hodge_stablegraphs_genus4.m"]; entries4 = Rest[Get["hodge_graphterms_genus4.m"]]; Print[Length[stable4], " stable graphs, ", Length[entries4], " decorated graphs"]]

(* total polar part {t_c^-6, ..., t_c^-1} versus the gap coefficient (-3)^3 5! [t^8] t^2/(4 sinh(e1 t/2) sinh(e2 t/2)) *)
If[ValueQ[stable4], {Factor /@ Total[stable4[[All, "polar"]]], Factor[gapCoefficient[4] /. {e1 -> x, e2 -> 1}]}]

(* which stable graphs are regular at the conifold, and which contain singular decorated graphs *)
If[ValueQ[stable4], regular4 = Select[Range[Length[stable4]], Union[Together /@ stable4[[#, "polar"]]] === {0} &];   (* canonicalise the summed polar parts before testing for zero *)
  {Length[regular4], "regular of", Length[stable4], "; singular stable graphs:", Table[{k, stable4[[k, "genera"]], stable4[[k, "edges"]], Factor /@ stable4[[k, "polar"]]}, {k, Complement[Range[Length[stable4]], regular4]}]}]
If[ValueQ[stable4], Grid[Prepend[Table[{k, stable4[[k, "genera"]], stable4[[k, "edges"]], stable4[[k, "decoratedGraphs"]], stable4[[k, "singularDecoratedGraphs"]]}, {k, Select[Range[Length[stable4]], stable4[[#, "singularDecoratedGraphs"]] > 0 &]}], {"graph", "genera", "edges", "decorated graphs", "singular decorated graphs"}], Frame -> All]]

(* the stable graphs with a vertex of genus >= 2 are exactly those containing singular decorated graphs, and the conifold pole
   order of every decorated graph is bounded by Sum_v max(2 g_v - 2, 0) over its vertices: {sets equal, their number, bound holds,
   tally of {pole order, bound}} *)
If[ValueQ[stable4], Module[{ge2 = Select[Range[Length[stable4]], Max[stable4[[#, "genera"]]] >= 2 &], sg = Select[Range[Length[stable4]], stable4[[#, "singularDecoratedGraphs"]] > 0 &], ord},
   ord = Table[{Module[{p = Position[e[[5]], _?(# =!= 0 &), {1}, 1, Heads -> False]}, If[p === {}, 0, 7 - p[[1, 1]]]], Total[Max[2 # - 2, 0] & /@ stable4[[e[[1, 1]], "genera"]]]}, {e, entries4}];
   {ge2 === sg, Length[ge2], AllTrue[ord, #[[1]] <= #[[2]] &], Sort[Tally[ord]]}]]

(* the total F_{4,0} agrees with the untwisted genus-4 cluster result (q -> 0 part of genRmat_F40_cluster_result.m) *)
If[ValueQ[stable4] && FileExistsQ["ref/genRmat_F40_cluster_result.m"], Module[{ref4 = Get["ref/genRmat_F40_cluster_result.m"][[1]] /. {q0 -> 0, q1 -> 0, q2 -> 0, q3 -> 0, q4 -> 0}, tot4 = Together[Total[stable4[[All, "term"]]]]},
   SeedRandom[4]; Table[Together[(tot4 - ref4) /. {y -> RandomInteger[{1, 50}]/RandomInteger[{51, 97}], XX -> RandomInteger[{-30, 30}]/RandomInteger[{1, 17}], e1 -> RandomInteger[{1, 40}], e2 -> RandomInteger[{1, 40}]}], {5}]]]
