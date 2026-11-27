(* Plain-text export of resolvedconifold_hodge.nb (2026-11-27): Input cells as code, other cells as comments. *)

(* ::Title:: Resolved conifold x C2: Hodge-vertex graph sum on the fibrewise torus, genus 2 and 3 *)

(* ::Text:: Companion of localp2CY5_hodge_genus3.nb. The equivariant Gromov-Witten theory of X x C^2, X = O(-1) + O(-1) -> P^1, on the fibrewise torus with opposite weights on each O + O(-1) pair (alpha0 = 0, alpha1 = eps1, alpha2 = eps2 in the conventions of resolvedconifold_genus3.nb: tangent weights {0, eps1, eps2, -eps1, -eps2} along the fixed P^1) is reconstructed from two copies of the SAME triple Hodge class as for local P2, the equivariant GW theory of C^3 with weights w = {eps1+eps2, -eps1, -eps2} (mode B; mode A is the opposite sign), instead of two copies of the GW theory of the point: Omega_Z = R'.(Hodge_w)^{+2} with R' = R_Z h_w^(-1), h_w(z) = Exp[-Sum_k B_2k/(2k(2k-1)) z^(2k-1) Sum_a w_a^(1-2k)] (h_w . point = Hodge_w by Mumford/QRR; a scalar, so the totals F_{g,0} are unchanged by the group-action property). R_Z is the conifold R-matrix of resolvedconifold_genus3.nb built from the Wick expansion of the superpotential; at alpha0 = 0, alpha_i = eps_i it needs no QRR normalisation constants (they are h_{(eps1, eps2, -eps1, -eps2)} = 1, the Bernoulli series of the actual normal bundle O(-1)+O(-1)+O+O of P^1, whose rank-4 Hodge class is trivial by c(E) c(E^vee) = 1). Everything is exact: Q symbolic, eps1 = x symbolic, eps2 = 1; polar parts are taken at Q = 1 in t = Log Q (the gap singularity P_g Li_{3-2g}(Q)/(eps1 eps2), where one canonical coordinate runs into the logarithmic singularity of the superpotential). Evaluate top to bottom. Data files (in the notebook directory or its subfolder hodge_vertex_conifold/): StableGraph-2.bin, StableGraph-3.bin, ref/resolvedconifold_results.m (F_{1,1}, F_{2,0}, F_{3,0} on the general torus), ref/conifold_graphterms_genus2_3.m (point-vertex terms of 2026-09), the stored hvc_results_g3_*.m (genus 3 is read from them unless recompute3 = True), and for the comparison ../hodge_vertex/hodge_vertex_genus{2,3}_results.m. *)

(* ::Text:: Conventions. Vertices carry idempotent labels i in {1, 2} (the two roots U_1, U_2 of U^2 - S U + P = 0, S = (eps1 - eps2)/eps1, P = -eps2/(Q eps1), elements of the quadratic algebra C[a, b] = a + b U); graph factors as in resolvedconifold_genus3.nb (Delta_i^(g_v - 1) per vertex, dilaton leaves (-1)^n R'_{0,i}[z^(n-1)], leaves (-1)^n R'_{1,i}[z^n], edge terms from the rows of R' and eta^-1), with the Witten-Kontsevich vertex <prod tau> replaced by Int_{M_{g_v,n_v}} prod_a Lambda^vee_{g_v}(w_a)/w_a^{g_v} prod psi^a, so every psi assignment of degree at most the dimension contributes. A decorated graph is an orbit of labelled psi/dilaton decorations under Aut(Gamma); each stable graph contributes 1/|Aut(Gamma)| times the sum over labelled decorations and 2^V colourings, with 1/nd_v! for the dilaton leaves. Every orbit value must have zero U-component (Galois symmetry), which is checked. *)

(* ::Text:: Hodge integrals (genus <= 4, hv_lib4.m): every lambda class is a polynomial in lambda_1 and ch_3(E) (Newton's identities with ch_2 = ch_4 = 0), and Int lambda_1^k1 ch_3^k3 psi^a is computed by Mumford's formula for ch_1 and ch_3. Checks: known values up to genus 4, vanishing above the top degree, Mumford relations, closed forms, string and dilaton equations, and the R-action of the scalar h_w^-1 on the point theory in the PPZ conventions reproducing the normalised Hodge class on M_2 (symbolically) and M_3 (42 graphs, exact rational weights). The exponent s in R' = R_Z h_w^s is fixed by the requirement that the sum of the decorated graphs equals F_{g,0}: s = -1 for both signs of the weights, as for local P2. *)

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

(* ::Subsection:: Conifold Frobenius structure, quadratic algebra, Wick R-matrix, graph factors and graph sum (rc_lib.m) *)

(* rc_lib.m -- library for the resolved conifold x C^2 pipeline (general five-dimensional torus), genus <= 4.
   Load with Get after SetDirectory to the data directory (graphs.m, orbits.m, StableGraph-<g>.bin must be there).
   Sections 1-4 are verbatim copies of the Input cells of resolvedconifold_genus3.nb (setup, quadratic algebra, numeric
   R-matrix and graph factors, graph sum, degree-zero term, Q-reconstruction and weight interpolation); section 5 adds
   the closed-form QRR normalisation, the alpha0 = 0 R-matrix, the Laurent-series layer at Q = 1, the sample lists and
   the multivariate fit on the CY torus.  ASCII only; the weights are the symbols Alpha[0..2], Epsilon[1..2]
   (typed as backslash-bracket-Alpha / -Epsilon). *)


(* ============================== 1. setup: Frobenius structure, canonical coordinates, Delta (verbatim) ============================== *)
a0 = \[Alpha][0]; a1 = \[Alpha][1]; a2 = \[Alpha][2]; e1 = \[Epsilon][1]; e2 = \[Epsilon][2];
CYearly = False;
CYrepl = {\[Alpha][0] -> \[Epsilon][1] + \[Epsilon][2] - \[Alpha][1] - \[Alpha][2]};
cy[x_] := If[CYearly, x /. CYrepl, x];
Replp2 = -Expand[p (p + a0)] + p^2;
exp3 = Simplify[Expand[Expand[Expand[(t[0] + t[1] p)^3] /. p^4 -> p^2 Replp2] /. p^3 -> p Replp2] /. p^2 -> Replp2];
PrepCl = cy[Simplify[Normal[Series[((exp3/a0/(-p + a2)/(-p + a1)/e1/e2 /. p -> 0) + (exp3/(-a0)/(-p + a2)/(-p + a1)/e1/e2 /. p -> -a0))/6, {t[0], 0, 3}, {t[1], 0, 3}]]]];
eta = Table[Simplify[D[D[D[PrepCl, t[0]], t[i]], t[j]]], {i, 0, 1}, {j, 0, 1}];
etainv = Simplify[Inverse[eta]];
Prep = PrepCl + PolyLog[3, Exp[t[1]]]/(e1 e2);
Do[cch[i, j, k] = Simplify[Sum[etainv[[i + 1, m + 1]] D[D[D[Prep, t[m]], t[j]], t[k]], {m, 0, 1}] /. t[1] -> Log[Q]], {i, 0, 1}, {j, 0, 1}, {k, 0, 1}];
Leq = Together[L^2 - cch[0, 1, 1] - cch[1, 1, 1] L];
Ssum = Together[cch[1, 1, 1]];
Pprod = Together[-cch[0, 1, 1]];
RedLeq[arg_, var_] := Block[{out = Together[arg]},
   Together[PolynomialRemainder[Numerator[out], Leq /. L -> var, var]/PolynomialRemainder[Denominator[out], Leq /. L -> var, var]]];
LDer[Lv_] := Together[-D[Leq, Q]/D[Leq, L] /. L -> Lv];
Module[{Psi = {{1, L1}, {1, L2}}, etacan},
  etacan = Simplify[Transpose[Inverse[Psi]] . eta . Inverse[Psi]];
  DeltaL = RedLeq[1/etacan[[1, 1]] /. L2 -> Ssum - L1 /. L1 -> L, L]];

(* ============================== 2. arithmetic in the quadratic algebra (verbatim, series-aware cTog) ============================== *)
cTog[x_SeriesData] := x;
cTog[x_] := Together[x];
cadd[C[a_, b_], C[c_, d_]] := C[cTog[a + c], cTog[b + d]];
cadd[x__] := Fold[cadd, C[0, 0], {x}];
cscale[k_, C[a_, b_]] := C[cTog[k a], cTog[k b]];
cmul[C[a_, b_], C[c_, d_]] := C[cTog[a c - b d Pprod], cTog[a d + b c + b d Ssum]];
cmul[x__] := Fold[cmul, C[1, 0], {x}];
cconj[C[a_, b_]] := C[cTog[a + b Ssum], cTog[-b]];
cnorm[C[a_, b_]] := cTog[a^2 + a b Ssum + b^2 Pprod];
cinv[C[a_, b_]] := Module[{n = cnorm[C[a, b]]}, C[cTog[(a + b Ssum)/n], cTog[-b/n]]];
cpow[c_, 0] := C[1, 0];
cpow[c_, n_Integer?Positive] := cmul @@ ConstantArray[c, n];
cpow[c_, n_Integer?Negative] := cpow[cinv[c], -n];
genSym = L;
toC[e_] := Module[{r, num, den, a, b, c, d, v = genSym},
   r = Together[e]; num = PolynomialRemainder[Numerator[r], v^2 - Ssum v + Pprod, v];
   den = PolynomialRemainder[Denominator[r], v^2 - Ssum v + Pprod, v];
   a = Coefficient[num, v, 0]; b = Coefficient[num, v, 1]; c = Coefficient[den, v, 0]; d = Coefficient[den, v, 1];
   cmul[C[a, b], cinv[C[c, d]]]];
symPow[sym_, n_, symvals_] := symPow[sym, n, symvals] = cpow[symvals[sym], n];
termEval[t_, symvals_] := Module[{facs = If[Head[t] === Times, List @@ t, {t}], coef = 1, fac = C[1, 0]},
   Do[Which[KeyExistsQ[symvals, f], fac = cmul[fac, symvals[f]],
      MatchQ[f, Power[_, _Integer]] && KeyExistsQ[symvals, f[[1]]], fac = cmul[fac, symPow[f[[1]], f[[2]], symvals]],
      True, coef = coef f], {f, facs}];
   cscale[coef, fac]];
cEval[poly_, symvals_Association] := Module[{ex = Expand[poly], terms},
   terms = If[Head[ex] === Plus, List @@ ex, {ex}];
   cadd @@ (termEval[#, symvals] & /@ terms)];
cExpr[C[a_, b_]] := a + b genSym;
cSymQ[C[a_, b_]] := Module[{bb = If[Head[b] === SeriesData, Normal[b], b]}, bb === 0 || Together[bb] === 0];

(* ============================== 3. R-matrix at numeric weights, graph factors (verbatim + closed-form variants) ============================== *)
(* symbolic pieces computed once for a given ZMAX.
   wickZ[[m+1]] = coefficient of z^m of sum_n (2n)! (-z)^n/(n! 2^n W2^n) [x^2n] Exp[sum_{k>=3} W_k x^k/(k! z)] as a polynomial in the
   symbols Wk[k]: a term with j factors W_{k_1} ... W_{k_j}, sum k_i = 2n, has z-order n - j, so the z^m coefficient is the sum over
   n = m .. 3m and over the partitions of 2n into exactly n - m parts in [3, 2 ZMAX + 2], with the multiset factor
   prod_k 1/((k!)^{m_k} m_k!).  Enumerating the partitions directly avoids the x-series of the exponential (hopeless at z^9).
   The derivatives W_k of the superpotential are computed by computeWU only after the weights are substituted (the fully
   symbolic W_k explode for k > 14). *)
wickDirect[zmax_] := Table[Module[{tot = 0, j, parts},
     Do[j = n - m;
      parts = If[j == 0, If[n == 0, {{}}, {}], Select[IntegerPartitions[2 n, {j}], Min[#] >= 3 && Max[#] <= 2 zmax + 2 &]];   (* parts in [3, 2 zmax + 2]; the explicit part-list form of IntegerPartitions is very slow here *)
      Do[tot += (2 n)! (-1)^n/(n! 2^n Wk[2]^n) Times @@ (Wk[#[[1]]]^#[[2]]/((#[[1]]!)^#[[2]] #[[2]]!) & /@ Tally[p]), {p, parts}], {n, m, 3 m}];
     Expand[tot]], {m, 0, zmax}];
setupWick[zmax_] := (ZMAX = zmax;
   wickZ = wickDirect[zmax];
   spotU = a1 Log[1 - Q U] - a2 Log[U] + (a0 + a2) Log[1 - U];
   LofU = -(Q (a0 + a1) U + Q a2)/(1 - Q);
   etaSym = {{(a0 + a1 + a2)/(a1 (a0 + a1) a2 (a0 + a2) e1 e2), 1/((a0 + a1) (a0 + a2) e1 e2)}, {1/((a0 + a1) (a0 + a2) e1 e2), -(a0/((a0 + a1) (a0 + a2) e1 e2))}};
   etainvSym = Simplify[Inverse[etaSym]];
   wtSym[1] = {a0, a1, a2, -e1, -e2}; wtSym[2] = {-a0, a0 + a1, a0 + a2, -e1, -e2};
   "wick ready");
QrrcSym[i_, k_] := (BernoulliB[2 k + 2]/(2 k + 2)/(2 k + 1)) Sum[1/w^(2 k + 1), {w, wtSym[i]}];
(* derivatives W_k = (U d/dU)^k spot at the given weights (rational functions of U and Q; cheap even for k = 20) *)
computeWU[wsub_] := (WU[1] = Together[U D[spotU /. wsub, U]];
   Do[WU[k] = Together[U D[WU[k - 1], U]], {k, 2, 2 ZMAX + 2}]; "W ready");
(* QRR normalisation constants from the Q -> 0 limit (verbatim); returns {ffc list root 1, ffc list root 2} (ffc[0] = 1) *)
ffcLimit[wsub_, SU_, LeqU_] := Module[{U1ser, U2ser, lim, eqs, sol, out},
   computeWU[wsub];
   Module[{c, sol1}, sol1 = First[Solve[Thread[Rest[CoefficientList[Normal[Series[LeqU /. U -> Sum[c[j] Q^j, {j, 0, 4}], {Q, 0, 4}]], Q]] == 0] /. c[0] -> (-a2/a0 /. wsub), Table[c[j], {j, 1, 4}]]];
     U1ser = Sum[c[j] Q^j, {j, 0, 4}] /. c[0] -> (-a2/a0 /. wsub) /. sol1];
   U2ser = Together[SU - U1ser];
   out = Table[Module[{Wser, ff},
      Wser = Association[Table[Wk[k] -> (WU[k] /. U -> If[i == 1, U1ser, U2ser]), {k, 2, 2 ZMAX + 2}]];
      lim = Table[Normal[Series[wickZ[[m + 1]] /. Wser, {Q, 0, 0}]], {m, 0, ZMAX}];
      If[! FreeQ[lim, Q], Print["WARNING: Q->0 limit not finite at root ", i, ": ", lim]];
      ff = Normal[Series[Exp[-Sum[QrrcSym[i, k] z^(2 k + 1) /. wsub, {k, 0, ZMAX}]], {z, 0, ZMAX}]];
      eqs = CoefficientList[Normal[Series[Sum[lim[[m + 1]] z^m, {m, 0, ZMAX}] (1 + Sum[ffc[k] z^k, {k, 1, ZMAX}]), {z, 0, ZMAX}]] - ff, z];
      sol = First[Solve[Thread[eqs == 0], Table[ffc[k], {k, 1, ZMAX}]]];
      Prepend[Table[ffc[k] /. sol, {k, 1, ZMAX}], 1]], {i, 1, 2}];
   out];
(* closed (Stirling) form of the same constants, identical for both roots; verified against ffcLimit to z^6 at general weights *)
ffcClosed[wsub_] := Module[{S, lst},
   S[k_] := (a1^-(2 k - 1) + (a0 + a2)^-(2 k - 1) - e1^-(2 k - 1) - e2^-(2 k - 1)) /. wsub;
   lst = CoefficientList[Normal[Series[Exp[-Sum[BernoulliB[2 k]/(2 k (2 k - 1)) z^(2 k - 1) S[k], {k, 1, ZMAX}]], {z, 0, ZMAX}]], z];
   PadRight[lst, ZMAX + 1]];
(* core: R-matrix at the weights wsub (Q symbolic) with given normalisation lists; sets ratCN[a,i][m], DeltaCN[i], SsumQ, PprodQ, etainvCur *)
rmatrixCore[wsub_, ffc1_, ffc2_] := Module[{SU, PU, LeqU, dUC, etainvN, DeltaLN, m, i, k},
   genSym = U;
   SU = (a0 + Q (a1 - a2))/(Q (a0 + a1)) /. wsub; PU = -a2/(Q (a0 + a1)) /. wsub;
   SsumQ = Together[SU]; PprodQ = Together[PU]; Ssum = SsumQ; Pprod = PprodQ;
   LeqU = Q (a0 + a1) U^2 - (a0 + Q (a1 - a2)) U - a2 /. wsub;
   computeWU[wsub];
   Do[WcN[k] = toC[WU[k]], {k, 2, 2 ZMAX + 2}];
   wsymN = Association[Table[Wk[k] -> WcN[k], {k, 2, 2 ZMAX + 2}]];
   Clear[symPow];
   symPow[sym_, n_, symvals_] := symPow[sym, n, symvals] = cpow[symvals[sym], n];
   Do[wickCN[m] = cEval[wickZ[[m + 1]], wsymN], {m, 0, ZMAX}];
   Do[ratCN[0, 1][m] = cadd @@ Table[cscale[ffc1[[k + 1]], wickCN[m - k]], {k, 0, m}], {m, 0, ZMAX}];
   Do[ratCN[0, 2][m] = cadd @@ Table[cscale[ffc2[[k + 1]], cconj[wickCN[m - k]]], {k, 0, m}], {m, 0, ZMAX}];
   dUC = toC[-D[LeqU, Q]/D[LeqU, U]];
   cDerN[C[a_, b_]] := cadd[C[Together[D[a, Q]], Together[D[b, Q]]], cscale[b, dUC]];
   etainvN = Together[etainvSym /. wsub];
   LCN[1] = toC[LofU /. wsub]; LCN[2] = cconj[LCN[1]];
   DeltaLN = DeltaL /. wsub /. L -> (LofU /. wsub);
   DeltaCN[1] = toC[DeltaLN]; DeltaCN[2] = cconj[DeltaCN[1]];
   Do[dlogDeltaN[i] = cmul[cDerN[DeltaCN[i]], cinv[cscale[2, DeltaCN[i]]]], {i, 1, 2}];
   Do[ratCN[1, i][0] = cmul[LCN[i], ratCN[0, i][0]];
      Do[ratCN[1, i][m] = cadd[cmul[LCN[i], ratCN[0, i][m]], cscale[Q, cadd[cDerN[ratCN[0, i][m - 1]], cscale[-1, cmul[ratCN[0, i][m - 1], dlogDeltaN[i]]]]]], {m, 1, ZMAX}], {i, 1, 2}];
   etainvCur = etainvN;
   "ok"];
(* the three entry points *)
rmatrixAt[wsub_] := Module[{SU = (a0 + Q (a1 - a2))/(Q (a0 + a1)) /. wsub, LeqU = Q (a0 + a1) U^2 - (a0 + Q (a1 - a2)) U - a2 /. wsub, ffs},
   ffs = ffcLimit[wsub, SU, LeqU]; ffcLast = ffs; rmatrixCore[wsub, ffs[[1]], ffs[[2]]]];              (* QRR constants from the Q -> 0 limit *)
rmatrixAtClosed[wsub_] := Module[{ff = ffcClosed[wsub]}, ffcLast = {ff, ff}; rmatrixCore[wsub, ff, ff]];   (* closed-form constants *)
rmatrixAt0[wsub_] := rmatrixCore[Join[{a0 -> 0}, wsub], PadRight[{1}, ZMAX + 1], PadRight[{1}, ZMAX + 1]];  (* alpha0 = 0, alpha_i = eps_i: no constants *)
(* everything at a numeric (or symbolic) Q *)
evalQ[Q0_] := Module[{sub = {Q -> Q0}, i, j, m, n, k, pre},
   Ssum = Together[SsumQ /. sub]; Pprod = Together[PprodQ /. sub];
   Do[Do[rN[a, i][m] = ratCN[a, i][m] /. sub, {m, 0, ZMAX}], {a, 0, 1}, {i, 1, 2}];
   Do[DeltaN[i] = DeltaCN[i] /. sub, {i, 1, 2}];
   Do[DilN[i][n] = cscale[(-1)^n, rN[0, i][n - 1]], {i, 1, 2}, {n, 1, ZMAX + 1}];
   Do[LfN[i][n] = cscale[(-1)^n, rN[1, i][n]], {i, 1, 2}, {n, 0, ZMAX}];
   Do[pre = Table[cadd @@ Flatten[Table[cscale[-(-1)^(p + q) etainvCur[[a + 1, b + 1]], cmul[rN[a, i][p], rN[b, j][q]]], {a, 0, 1}, {b, 0, 1}]], {p, 0, ZMAX}, {q, 0, ZMAX}];
      Do[EdN[i, j][m, n] = cadd @@ Table[cscale[(-1)^k, pre[[m + k + 2, n - k + 1]]], {k, 0, n}], {m, 0, ZMAX - 1}, {n, 0, ZMAX - 1 - m}], {i, 1, 2}, {j, 1, 2}];
   "ok"];
(* the same with every factor expanded as a Laurent series in t = Log Q at Q = 1 to absolute order M.
   The vertex factors Delta_i^(g_v - 1), g_v - 1 = -1 .. gmax - 1, are computed exactly as rational functions first and expanded
   afterwards (series inversion of Delta loses precision and mislabels the valuation of Delta^-1), via deltaPow[i, k]. *)
toS[expr_] := Series[expr /. Q -> Exp[t], {t, 0, M}];
deltaPow[i_, k_] := cpow[DeltaN[i], k];                    (* rational mode: plain power; series mode: redefined by evalQseries *)
evalQseries[Mord_, gmax_: 4] := Module[{i, j, m, n, k, pre, ex},
   M = Mord;
   Ssum = Together[SsumQ]; Pprod = Together[PprodQ];       (* exact powers of Delta with the rational algebra data *)
   Do[ex = cpow[DeltaCN[i], k]; deltaPow[i, k] = C[toS[ex[[1]]], toS[ex[[2]]]], {i, 1, 2}, {k, -1, gmax - 1}];
   Ssum = toS[SsumQ]; Pprod = toS[PprodQ];
   Do[Do[rN[a, i][m] = C[toS[ratCN[a, i][m][[1]]], toS[ratCN[a, i][m][[2]]]], {m, 0, ZMAX}], {a, 0, 1}, {i, 1, 2}];
   Do[DeltaN[i] = C[toS[DeltaCN[i][[1]]], toS[DeltaCN[i][[2]]]], {i, 1, 2}];
   Do[DilN[i][n] = cscale[(-1)^n, rN[0, i][n - 1]], {i, 1, 2}, {n, 1, ZMAX + 1}];
   Do[LfN[i][n] = cscale[(-1)^n, rN[1, i][n]], {i, 1, 2}, {n, 0, ZMAX}];
   Do[pre = Table[cadd @@ Flatten[Table[cscale[-(-1)^(p + q) etainvCur[[a + 1, b + 1]], cmul[rN[a, i][p], rN[b, j][q]]], {a, 0, 1}, {b, 0, 1}]], {p, 0, ZMAX}, {q, 0, ZMAX}];
      Do[EdN[i, j][m, n] = cadd @@ Table[cscale[(-1)^k, pre[[m + k + 2, n - k + 1]]], {k, 0, n}], {m, 0, ZMAX - 1}, {n, 0, ZMAX - 1 - m}], {i, 1, 2}, {j, 1, 2}];
   "ok"];

(* ============================== 4. graph sum, degree-zero term, Q-reconstruction, weight interpolation (verbatim) ============================== *)
graphValueC[gr_] := Module[{V = gr["V"], ed = gr["edges"], H = gr["H"], hsAt = gr["hsAt"], total = C[0, 0], hp, facs},
   Do[Do[hp = ConstantArray[0, Length[H]];
     Do[Do[hp[[hsAt[[v, r]]]] = dec[[v, 1, r]], {r, Length[hsAt[[v]]]}], {v, V}];
     facs = Join[
       Table[deltaPow[c[[v]], gr["g"][[v]] - 1], {v, V}],
       Flatten[Table[Table[LfN[c[[v]]][b], {b, dec[[v, 2]]}], {v, V}]],
       Flatten[Table[Table[DilN[c[[v]]][d], {d, dec[[v, 3]]}], {v, V}]],
       Table[EdN[c[[ed[[e, 1]]]], c[[ed[[e, 2]]]]][hp[[2 e - 1]], hp[[2 e]]], {e, Length[ed]}]];
     total = cadd[total, cscale[Product[dec[[v, 4]], {v, V}], cmul @@ facs]], {dec, gr["decs"]}], {c, Tuples[{1, 2}, V]}];
   cscale[1/gr["aut"], total]];
graphSumC[prepared_] := cadd @@ (graphValueC /@ prepared);
(* degree-zero (constant map) part, symbolic in the five weights *)
const0[gg_] := Module[{ffs, Del0, res = 0, prep},
   Do[ffs[i] = CoefficientList[Normal[Series[Exp[-Sum[QrrcSym[i, k] z^(2 k + 1), {k, 0, ZMAX}]], {z, 0, ZMAX}]], z];
      Del0[i] = Together[DeltaL /. Q -> 0 /. L -> {0, -a0}[[i]]], {i, 1, 2}];
   prep = prepareGraph /@ graphsWithLegs[gg, 0];
   Do[Module[{i = col, tot = 0},
      Do[Dil0[n] = (-1)^n ffs[i][[n]], {n, 1, ZMAX + 1}];
      Do[Ed0[m, n] = Del0[i] Sum[(-1)^k (-(-1)^(m + k + 1 + n - k) ffs[i][[m + k + 2]] ffs[i][[n - k + 1]]), {k, 0, n}], {m, 0, ZMAX - 1}, {n, 0, ZMAX - 1 - m}];
      Do[Module[{V = gr["V"], ed = gr["edges"], H = gr["H"], hsAt = gr["hsAt"], hp},
        Do[hp = ConstantArray[0, Length[H]];
         Do[Do[hp[[hsAt[[v, r]]]] = dec[[v, 1, r]], {r, Length[hsAt[[v]]]}], {v, V}];
         tot += (1/gr["aut"]) Product[dec[[v, 4]], {v, V}] Product[Dil0[d], {v, V}, {d, dec[[v, 3]]}] Product[Ed0[hp[[2 e - 1]], hp[[2 e]]], {e, Length[ed]}] Del0[i]^(Total[gr["g"]] - V),
         {dec, gr["decs"]}]], {gr, prep}];
      res += Together[tot]], {col, 1, 2}];
   Together[res]];
F11const = (-(\[Alpha][0]^3*\[Epsilon][1]) - \[Alpha][0]^2*\[Alpha][1]*\[Epsilon][1] - \[Alpha][0]^2*\[Alpha][2]*\[Epsilon][1] - \[Alpha][0]*\[Alpha][1]*\[Alpha][2]*\[Epsilon][1] - \[Alpha][0]^3*\[Epsilon][2] - \[Alpha][0]^2*\[Alpha][1]*\[Epsilon][2] - \[Alpha][0]^2*\[Alpha][2]*\[Epsilon][2] - \[Alpha][0]*\[Alpha][1]*\[Alpha][2]*\[Epsilon][2] + \[Alpha][0]^2*\[Epsilon][1]*\[Epsilon][2] - \[Alpha][1]*\[Alpha][2]*\[Epsilon][1]*\[Epsilon][2])/(24*(\[Alpha][0] + \[Alpha][1])*(\[Alpha][0] + \[Alpha][2])*\[Epsilon][1]*\[Epsilon][2]);
(* reconstruction of F(Q) from exact samples {q, v}: F = P(Q)/(1-Q)^a, verified on 5 spare points *)
reconQ[pts_, amax_: 16] := Module[{n = Length[pts], qs = pts[[All, 1]], vs = pts[[All, 2]], w, P, ok, res = $Failed},
   If[Union[vs] === {0}, Return[{0, 0}]];
   Do[w = Table[vs[[i]] (1 - qs[[i]])^a, {i, n}];
    P = InterpolatingPolynomial[Transpose[{qs[[1 ;; n - 5]], w[[1 ;; n - 5]]}], Q];
    ok = And @@ Table[(P /. Q -> qs[[i]]) == w[[i]], {i, n - 4, n}];
    If[ok, res = {Expand[P], a}; Return[]], {a, 0, amax}];
   res];
(* interpolation of the instanton part over the weight samples: eps1 eps2 [Q^j] F_g homogeneous of degree 2g, symmetric under a1<->a2, e1<->e2 *)
wvars = {a0, a1, a2, e1, e2};
orbitBasis[deg_] := Module[{mons, orb},
   mons = Select[Tuples[Range[0, deg], 5], Total[#] == deg &];
   orb[m_] := Union[{m, m[[{1, 3, 2, 4, 5}]], m[[{1, 2, 3, 5, 4}]], m[[{1, 3, 2, 5, 4}]]}];
   Union[Table[Total[Times @@@ (wvars^# & /@ orb[m])], {m, mons}]]];
fitPoly[values_Association, deg_] := Module[{basis = orbitBasis[deg], keys = Keys[values], A, rhs, sol, n, res},
   n = Length[basis];
   If[Length[keys] < n + 3, Print["not enough samples: ", Length[keys], " < ", n + 3]; Return[$Failed]];
   A = Table[basis /. allSamples[k][[1]], {k, keys[[1 ;; n]]}]; rhs = Table[values[k], {k, keys[[1 ;; n]]}];
   sol = LinearSolve[A, rhs];
   res = Table[(basis /. allSamples[k][[1]]) . sol - values[k], {k, keys[[n + 1 ;;]]}];
   If[Union[res] =!= {0}, Print["fitPoly: verification FAILED on ", Count[res, Except[0]], " of ", Length[res], " spare samples"]; Return[$Failed]];
   Expand[basis . sol]];
(* assemble the instanton part of genus gi (record index gi in allSamples[k][[2]]) as N(Q)/(eps1 eps2 (1-Q)^a) *)
assemble[gi_, gdeg_] := Module[{a = Max[Table[allSamples[k][[2, gi, 2]], {k, Keys[allSamples]}]], cj, vals, Npoly = 0},
   Do[vals = Association[Table[k -> With[{NN = Expand[allSamples[k][[2, gi, 1]] (1 - Q)^(a - allSamples[k][[2, gi, 2]])], w = allSamples[k][[1]]},
        Coefficient[NN - (NN /. Q -> 0) (1 - Q)^a, Q, j] (e1 e2 /. w)], {k, Keys[allSamples]}]];
      cj = fitPoly[vals, gdeg];
      If[cj === $Failed, Return[$Failed]];
      Npoly += cj Q^j, {j, 1, a}];
   {Npoly/(e1 e2 (1 - Q)^a), a}];

(* ============================== 5. sample lists, orbits, polar parts, CY-torus fit ============================== *)
zmaxFor[gg_] := 3 gg - 3;
okQ[w_] := And @@ (# != 0 & /@ {w[[1]], w[[2]], w[[3]], w[[1]] + w[[2]], w[[1]] + w[[3]], w[[4]], w[[5]], w[[1]] + w[[2]] + w[[3]], w[[4]] + w[[5]], w[[4]] + w[[5]] - w[[1]] - w[[2]] - w[[3]]});
(* integer weight samples for the general torus (same generator and seed as resolvedconifold_genus3.nb: the first 100 are its samples) *)
wsamplesAll := wsamplesAll = Module[{ws = {}}, SeedRandom[20260906]; While[Length[ws] < 400, Module[{w = RandomInteger[{-9, 9}, 5]}, If[okQ[w] && ! MemberQ[ws, w], AppendTo[ws, w]]]]; ws];
nSamplesFor[gg_] := Length[orbitBasis[2 gg]] + 12;                     (* unknowns of the degree-2g symmetric basis plus spare samples *)
qsFor[gg_] := Table[1/(3 + 2 s) + 1/(50 + 3 s), {s, 0, If[gg <= 3, 15, 17]}];
wsubOf[w_] := {a0 -> w[[1]], a1 -> w[[2]], a2 -> w[[3]], e1 -> w[[4]], e2 -> w[[5]]};
(* orbit table (prepared graphs + orbits), stored by rc_prep.m in rc_orbits_g<gg>.mx *)
setupOrbits[gg_] := (prepared = prepareGraph /@ graphsWithLegs[gg, 0]; OT = orbitTable[prepared];
   Print["genus ", gg, ": ", Length[prepared], " graphs, ", Total[Length[#["decs"]] & /@ prepared], " labelled decorations, ", Length[OT["ids"]], " decorated graphs (orbits)"]);
loadOrbits[gg_] := Module[{f = "rc_orbits_g" <> ToString[gg] <> ".mx"},
   If[FileExistsQ[f], Get[f], Print["orbit file ", f, " not found: computing"]; setupOrbits[gg]; DumpSave[f, {prepared, OT}]]];
(* polar part at Q = 1 in t = Log Q of a rational function of Q: coefficients of t^-(2g-2) .. t^-1 *)
polarOfRational[rq_, gg_] := Module[{ser = Series[rq /. Q -> Exp[t], {t, 0, -1}]}, Table[SeriesCoefficient[ser, -n], {n, 2 gg - 2, 1, -1}]];
polarOfSeries[ser_, gg_] := Table[SeriesCoefficient[ser, -n], {n, 2 gg - 2, 1, -1}];
(* the known total on the CY torus: P_g Li_{3-2g}(Q)/(eps1 eps2) (plus a constant) *)
Pg[g_] := Module[{tv}, SeriesCoefficient[Series[e1 e2 tv^2/(4 Sinh[e1 tv/2] Sinh[e2 tv/2]), {tv, 0, 2 g}], 2 g]];
cyPolar[gg_] := polarOfRational[Pg[gg] PolyLog[3 - 2 gg, Q]/(e1 e2), gg];
(* multivariate fit on the CY torus (alpha0 = eps1+eps2-alpha1-alpha2, eps2 = 1): coefficient = N/D with D = eps1^da eps2^da (eps1+eps2)^db,
   N homogeneous of degree d in (alpha1, alpha2, eps1, eps2), symmetric under alpha1<->alpha2 and eps1<->eps2, alpha-degree <= Da *)
symA[j_] := Table[If[a == j - a, a1^a a2^a, a1^a a2^(j - a) + a1^(j - a) a2^a], {a, 0, Floor[j/2]}];
symE[j_] := Table[If[a == j - a, e1^a e2^a, e1^a e2^(j - a) + e1^(j - a) e2^a], {a, 0, Floor[j/2]}];
cyBasis[d_, Da_] := Flatten[Table[Outer[Times, symA[j], symE[d - j]], {j, 0, Min[Da, d]}]];
cySub[pt_] := {a1 -> pt[[1]], a2 -> pt[[2]], e1 -> pt[[3]], e2 -> 1};
(* line probe: along the line pt = base + s dir (integer s) a coefficient is a univariate rational function of s; the generic
   reconstruction gives its numerator/denominator degrees, which fixes d, da, db when D is of the assumed form *)
lineRecon[vals_] := Module[{pts = Transpose[{Range[Length[vals]], vals}]}, reconstructU[pts]];
reconstructU[pts_] := Module[{n = Length[pts], us = pts[[All, 1]], vs = pts[[All, 2]], res = $Failed, T, dd, dn, m, A, rhs, sol, N0, D0, ok},
   If[Union[vs] === {0}, Return[0]];
   Do[Do[dn = T - dd; m = dn + dd + 1; If[m + 4 > n, Continue[]];
     A = Table[Join[Table[us[[s]]^j, {j, 0, dn}], Table[-vs[[s]] us[[s]]^j, {j, 0, dd - 1}]], {s, 1, m}];
     rhs = Table[vs[[s]] us[[s]]^dd, {s, 1, m}];
     sol = Quiet[Check[LinearSolve[A, rhs], $Failed]]; If[sol === $Failed, Continue[]];
     N0 = Sum[sol[[j + 1]] u^j, {j, 0, dn}]; D0 = u^dd + Sum[sol[[dn + 2 + j]] u^j, {j, 0, dd - 1}];
     ok = And @@ Table[((N0 /. u -> us[[s]]) - vs[[s]] (D0 /. u -> us[[s]])) == 0, {s, m + 1, n}];
     If[ok, res = Together[N0/D0]; Return[]], {dd, 0, T}]; If[res =!= $Failed, Return[]], {T, 0, n - 5}];
   res];
(* group fit: vals[key] = list of values at the points pts (random points first, probe points last); returns association key -> rational function *)
fitGroup[vals_Association, pts_, d_, da_, db_, DaMax_, nspare_: 8] := Module[{keys = Keys[vals], Dp, W, res = <||>, still, Da = 0, bs, nU, Ball, lsf, done = False},
   Dp = e1^da e2^da (e1 + e2)^db;
   W = Association[Table[k -> Table[vals[k][[s]] (Dp /. cySub[pts[[s]]]), {s, Length[pts]}], {k, keys}]];
   still = keys;
   While[! done && Da <= DaMax,
    bs = cyBasis[d, Da]; nU = Length[bs];
    If[nU + nspare > Length[pts], Print["fitGroup: not enough points for alpha-degree ", Da, " (", nU, " unknowns, ", Length[pts], " points)"]; Break[]];
    Ball = Table[bs /. cySub[pt], {pt, pts}];
    lsf = Quiet[Check[LinearSolve[Ball[[;; nU]]], $Failed]];
    If[lsf === $Failed, Da++; Continue[]];
    still = Select[still, Module[{sol = lsf[W[#][[;; nU]]]}, If[Ball[[nU + 1 ;;]] . sol == W[#][[nU + 1 ;;]], res[#] = Factor[(sol . bs)/Dp]; False, True]] &];
    If[still === {}, done = True, Da++]];
   {res, still, Da}];

(* ::Subsection:: Hodge integrals, Hodge vertices, twisted rows R' = R_Z h^s, references (hvc_lib.m) *)

(* hvc_lib.m -- Hodge-vertex graph sum for the resolved conifold x C^2 on the fibrewise torus alpha0 = 0, alpha1 = eps1,
   alpha2 = eps2 (opposite weights on each O + O(-1) pair), genus <= 3 locally (the Hodge machinery is valid to genus 4).
   Omega_Z = R' . (Hodge_w)^{+2} with the SAME triple Hodge class as for local P2: vertex = Int prod_{a=1}^3 Lambda_g^vee(w_a)/w_a^g psi^a
   (GW theory of C^3), mode B: w = {eps1+eps2, -eps1, -eps2}, mode A: w = {-eps1-eps2, eps1, eps2} (A and B differ by the overall sign
   of the weights, Lambda^vee <-> Lambda); R' = R_Z h_w^s with h_w(z) = Exp[-Sum_k B_2k/(2k(2k-1)) z^(2k-1) Sum_a w_a^(1-2k)] and
   s = -1 or +1 (fixed by the requirement that the sum of the decorated graphs is F_{g,0}; for local P2, s = -1 in these conventions).
   R_Z is the conifold R-matrix of resolvedconifold_genus3.nb / rc_lib.m at alpha0 = 0, alpha_i = eps_i, where no QRR normalisation
   constants arise (rmatrixAt0). Everything is exact: Q symbolic, eps1 = x symbolic, eps2 = 1; polar parts at Q = 1 in t = Log Q.
   Sections 1-4 are verbatim from hv_lib4.m (hodge_vertex/, local P2); section 5 is conifold specific. ASCII only. *)

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


(* ============================== 5. conifold specific: modes, twisted rows, factors, orbits, references ============================== *)
modeWeights["B"] = {e1s + e2s, -e1s, -e2s};
modeWeights["A"] = {-e1s - e2s, e1s, e2s};
(* the normalised Hodge class for an arbitrary list of nw weights (u[a] = 1/w_a), same truncation as hodgeClass4 *)
hodgeClassN[g_, nw_] := Module[{cl},
   cl = Expand[Product[Sum[(-1)^j lam[j] u[a]^j, {j, 0, g}], {a, nw}] /. lam[0] -> 1 /. lamSub];
   Select[{#[[1]], Factor[#[[2]]]} & /@ CoefficientRules[cl, {L1, C3}], #[[1, 1]] + 3 #[[1, 2]] <= dmax[g] &]];
(* multiply every row of R_Z (both roots) by the scalar series with coefficient list hc: ratCN[a,i][m] -> Sum_j ratCN[a,i][j] hc[[m-j+1]] *)
twistRows[hc_List] := Module[{old},
   Do[old[a, i, m] = ratCN[a, i][m], {a, 0, 1}, {i, 1, 2}, {m, 0, ZMAX}];
   Do[ratCN[a, i][m] = cadd @@ Table[cscale[hc[[m - j + 1]], old[a, i, j]], {j, 0, m}], {a, 0, 1}, {i, 1, 2}, {m, 0, ZMAX}];
   "rows twisted"];
(* R' = R_Z h_w^s at alpha0 = 0, alpha1 = eps1 = x0, alpha2 = eps2 = 1 (x0 symbolic or rational), then all graph factors at symbolic Q *)
buildFactors[zmax_, mode_, s_, x0_] := (setupWick[zmax];
   rmatrixAt0[{a1 -> x0, a2 -> 1, e1 -> x0, e2 -> 1}];
   If[mode =!= "point", twistRows[hCoeffs[modeWeights[mode] /. {e1s -> x0, e2s -> 1}, s, zmax]]];
   evalQ[Q]; "factors ready");
usubMode[mode_, x0_] := If[mode === "point", {}, usubOf[modeWeights[mode] /. {e1s -> x0, e2s -> 1}]];
modeTag[mode_, s_] := mode <> If[mode === "point", "", If[s < 0, "m", "p"]];
(* decorated graphs: Hodge vertices (prepareGraphH4) or point vertices (prepareGraph) *)
orbitFile[gg_, mode_] := "hvc_orbits_g" <> ToString[gg] <> If[mode === "point", "_point", "_hodge"] <> ".mx";
loadOrbitsHVC[gg_, mode_] := Module[{f = orbitFile[gg, mode]},
   If[FileExistsQ[f], Get[f],
    prepared = If[mode === "point", prepareGraph, prepareGraphH4] /@ graphsWithLegs[gg, 0]; OT = orbitTable[prepared];
    DumpSave[f, {prepared, OT}]];
   {Length[prepared], Total[Length[#["decs"]] & /@ prepared], Length[OT["ids"]]}];
(* records of a chunk file: either one list of records (Put format of the first local runs) or one record per expression (PutAppend) *)
readChunk[f_] := Module[{r = Quiet[Check[ReadList[f], $Failed]]},
   If[r === $Failed, Return[$Failed]];                      (* unreadable (e.g. a record torn by a kill during PutAppend) *)
   r = DeleteCases[r, $Failed];
   If[Length[r] == 1 && ListQ[r[[1]]] && r[[1]] =!= {} && ListQ[r[[1, 1]]] && Length[r[[1, 1]]] == 6 && MatchQ[r[[1, 1, 1]], {_Integer, _Integer}], r[[1]], r]];
(* one decorated graph: {rational function of Q (and x), polar part {t^-(2g-2), ..., t^-1}, symmetric?} *)
orbitTerm[gr_, usub_, gg_] := Module[{v = graphValueC[If[usub === {}, gr, numGraph[gr, usub]]], rq, pol},
   rq = Together[v[[1]]];
   pol = If[rq === 0, ConstantArray[0, 2 gg - 2], Factor /@ polarOfRational[rq, gg]];   (* Series[0] is not a SeriesData: an identically vanishing term has no polar part *)
   {rq, pol, cSymQ[v]}];
(* references: F_{g,0} of resolvedconifold_results.m (genus <= 3, independent computation on the general torus) at alpha0 = 0,
   alpha_i = eps_i, eps1 = x0, eps2 = 1. For genus >= 4 the closed form F_g = F_g(0) + P_g Li_{3-2g}(Q)/(eps1 eps2) is used, with the
   constant-map term F_g(0) = -(eps1^(2g-2) + eps2^(2g-2)) |B_2g B_{2g-2}| / (4g (2g-2) (2g-2)!), which reproduces the genus-2 and
   genus-3 values -(eps1^2+eps2^2)/5760 and -(eps1^4+eps2^4)/1451520 of the stored results (a prediction at genus 4, flagged in the log). *)
F0closed[gg_, x0_] := -(x0^(2 gg - 2) + 1) Abs[BernoulliB[2 gg] BernoulliB[2 gg - 2]]/(4 gg (2 gg - 2) (2 gg - 2)!);
FrefHVC[gg_, x0_] := If[gg <= 3, Get["ref/resolvedconifold_results.m"][[gg]] /. {a0 -> 0, a1 -> e1, a2 -> e2} /. {e1 -> x0, e2 -> 1},
   F0closed[gg, x0] + (Pg[gg] PolyLog[3 - 2 gg, Q]/(e1 e2) /. {e1 -> x0, e2 -> 1})];
FrefKind[gg_] := If[gg <= 3, "resolvedconifold_results.m (independent computation)", "closed form F_g(0) + P_g Li_{3-2g}(Q)/(eps1 eps2) (prediction)"];
gapHVC[gg_, x0_] := Factor /@ (cyPolar[gg] /. {e1 -> x0, e2 -> 1});

(* ::Section:: Setup *)

(* data directory: this notebook's folder if it holds hvc_lib.m (StableGraph-2/3.bin, ref/resolvedconifold_results.m, the stored
   hvc_results_g*.m), otherwise its subfolder hodge_vertex_conifold/.  The code above is inlined from graphs.m, orbits.m, rc_lib.m
   and hvc_lib.m; their Get[...] lines were removed. *)
datadir = If[FileExistsQ[FileNameJoin[{NotebookDirectory[], "hvc_lib.m"}]], NotebookDirectory[], FileNameJoin[{NotebookDirectory[], "hodge_vertex_conifold"}]];
SetDirectory[datadir]; Print["data directory: ", Directory[]];
wB = modeWeights["B"]; wA = modeWeights["A"];   (* {eps1+eps2, -eps1, -eps2} and its negative; e1, e2 are the symbols Epsilon[1], Epsilon[2] *)

(* one complete run: R' = R_Z h_w^s with Q and x = eps1/eps2 symbolic (alpha0 = 0, alpha1 = eps1 = x, alpha2 = eps2 = 1), decorated
   graphs with Hodge (or point) vertices, exact term and polar part at Q = 1 of every decorated graph; returns an Association
   with the same keys as the stored hvc_results_g<g>_<tag>.m files *)
runMode[gg_, mode_, s_] := Module[{nfo, usub, entries, stable, total, Fref, pol, gap},
   buildFactors[3 gg - 3, mode, s, x]; nfo = loadOrbitsHVC[gg, mode]; usub = usubMode[mode, x];
   entries = Table[Module[{r = orbitTerm[OT["graphs"][id], usub, gg]}, {id, OT["orbits"][[id[[1]], id[[2]]]], Normal[OT["descriptions"][id]], r[[1]], r[[2]], r[[3]]}], {id, OT["ids"]}];
   stable = Table[Module[{es = Select[entries, #[[1, 1]] == k &]}, <|"graph" -> k, "genera" -> prepared[[k]]["g"], "edges" -> prepared[[k]]["edges"], "Aut" -> prepared[[k]]["aut"], "decoratedGraphs" -> Length[es],
        "singularDecoratedGraphs" -> Count[es, e_ /; Union[Together /@ e[[5]]] =!= {0}], "term" -> Together[Total[es[[All, 4]]]], "polar" -> (Together /@ Total[es[[All, 5]]])|>], {k, Length[prepared]}];
   total = Together[Total[stable[[All, "term"]]]]; Fref = FrefHVC[gg, x]; pol = Together /@ Total[stable[[All, "polar"]]]; gap = gapHVC[gg, x];
   <|"genus" -> gg, "mode" -> mode, "s" -> s, "counts" -> nfo, "stable" -> stable, "entries" -> entries, "total" -> total, "totalOK" -> (Together[total - Fref] === 0), "polar" -> pol, "gapOK" -> (Union[Together /@ (pol - gap)] === {0}), "symmetric" -> And @@ entries[[All, 6]]|>];
regQ[v_] := Union[Together /@ v] === {0};
tableOf[res_] := Grid[Prepend[Table[With[{st = res["stable"][[k]]}, {k, st["genera"], st["edges"], st["Aut"], st["decoratedGraphs"], st["singularDecoratedGraphs"], If[regQ[st["polar"]], "regular", Factor /@ st["polar"]]}], {k, Length[res["stable"]]}],
    {"graph", "genera", "edges", "|Aut|", "decorated graphs", "singular decorated graphs", "polar part {t^-(2g-2), ..., t^-1} at Q = 1 (eps1 = x, eps2 = 1)"}], Frame -> All];

(* the structural statements of a run: singular decorated graphs versus vertex genera, pole-order bound Sum_v max(2 g_v - 2, 0) *)
structure[res_] := Module[{st = res["stable"], en = res["entries"], gg = res["genus"], ge2, sg, ord},
   ge2 = Select[Range[Length[st]], Max[st[[#]]["genera"]] >= 2 &]; sg = Select[Range[Length[st]], st[[#]]["singularDecoratedGraphs"] > 0 &];
   ord = Table[{Module[{p = Position[Together /@ e[[5]], _?(# =!= 0 &), {1}, 1, Heads -> False]}, If[p === {}, 0, 2 gg - 1 - p[[1, 1]]]], Total[Max[2 # - 2, 0] & /@ st[[e[[1, 1]]]]["genera"]]}, {e, en}];
   <|"graphs with a vertex of genus >= 2" -> ge2, "graphs containing singular decorated graphs" -> sg, "equal" -> (ge2 === sg), "pole order <= bound for every decorated graph" -> AllTrue[ord, #[[1]] <= #[[2]] &], "tally of {pole order, bound}" -> Sort[Tally[ord]]|>];

(* ::Section:: Checks *)

(* Hodge integrals of the genus-general evaluator (hv_lib4.m): known values up to genus 4, vanishing above the top degree, Mumford
   relations, closed forms versus the recursion, string and dilaton equations *)
{hodgeOK, hodgeReport} = hodgeChecks4[]; Column[Append[hodgeReport, "all checks passed: " <> ToString[hodgeOK]]]

(* the scalar R-action of h_w^-1 on the point theory (PPZ conventions) reproduces the normalised Hodge class on M_2 for arbitrary weights ... *)
ppz2 = ppzScalar[2, {w1, w2, w3}]; {Factor[ppz2], Together[ppz2 - (hodgeNormU4[2, {}] /. usubOf[{w1, w2, w3}])] === 0}

(* ... and on M_3 (42 stable graphs) at exact rational weights *)
SeedRandom[7]; Table[Module[{wn = Table[RandomInteger[{-9, 9}]/RandomInteger[{1, 5}], {3}]}, If[MemberQ[wn, 0], Nothing, {wn, ppzScalar[3, wn] == (hodgeNormU4[3, {}] /. usubOf[wn])}]], {2}]

(* the rank-4 Hodge class of the actual normal bundle O(-1)+O(-1)+O+O of P^1 in the 5-fold, weights (eps1, eps2, -eps1, -eps2), is
   trivial by Mumford's c(E) c(E^vee) = 1: it integrates to the Witten-Kontsevich numbers for every psi assignment (as a polynomial in
   lambda_1, ch_3 it is not literally 1, since lambda_4 = 0 in genus 3 and lambda_g^2 = 0 are not in the substitution table) *)
usub4 = {u[1] -> 1/e1, u[2] -> 1/e2, u[3] -> -1/e1, u[4] -> -1/e2};
hodgeNormN[g_, nw_, a_List] := Expand[Total[Table[c[[2]] HI[g, c[[1, 1]], c[[1, 2]], a], {c, hodgeClassN[g, nw]}]]];
Table[Module[{as = Select[Join @@ Table[Tuples[Range[0, 3 g - 3 + n], n], {n, 1, If[g == 3, 3, 4]}], Total[#] <= 3 g - 3 + Length[#] &]}, {g, Length[as], "disagreements with Int psi^a:", Count[as, a_ /; Together[(hodgeNormN[g, 4, a] /. usub4) - wk[g, a]] =!= 0]}], {g, 1, 3}]

(* the triple Hodge class of local P2 (mode B) in genus 2 and 3, as coefficients of lambda_1^k1 ch_3^k3, and the scalar h_wB(z) to z^6 *)
Table[{g, Table[{c[[1]], Factor[c[[2]] /. usubOf[wB]]}, {c, hodgeClass4[g]}]}, {g, 2, 3}] // TableForm
Factor /@ hCoeffs[wB, 1, 6]

(* the conifold Frobenius data at alpha0 = 0, alpha_i = eps_i: the quadratic equation of the canonical coordinates U^2 - S U + P = 0,
   whose roots are U = 1 and U = -1/x at Q = 1 (one canonical coordinate runs into the logarithmic singularity of the superpotential),
   and the discriminant Q (1-x)^2 + 4x (times x^2 Q), which vanishes at Q_c = -4x/(1-x)^2, not at Q = 1 *)
{Together[(a0 + Q (a1 - a2))/(Q (a0 + a1)) /. {a0 -> 0, a1 -> x, a2 -> 1}], Together[-a2/(Q (a0 + a1)) /. {a0 -> 0, a1 -> x, a2 -> 1}], Factor[Discriminant[Q x U^2 - Q (x - 1) U - 1, U]], Solve[(U^2 - (x - 1)/x U - 1/x) == 0, U]}

(* ::Section:: Genus 2 *)

(* genus 2 (a few seconds per mode, everything exact in Q and x): point vertices with R_Z; Hodge vertices of mode B (weights
   eps1+eps2, -eps1, -eps2) and of mode A (the opposite sign), each with R' = R_Z h_w^-1 *)
res2 = Association[Table[m -> runMode[2, m, -1], {m, {"point", "B", "A"}}]];
Table[{m, res2[m]["counts"], "all values symmetric" -> res2[m]["symmetric"], "total = F20" -> res2[m]["totalOK"], "total polar part = gap" -> res2[m]["gapOK"]}, {m, Keys[res2]}] // TableForm

(* the opposite exponent s = +1 does not reproduce F20 (shown for mode B) *)
Module[{r = runMode[2, "B", 1]}, {r["totalOK"], Factor[r["total"] - FrefHVC[2, x]]}]

(* mode B per stable graph: only the genus-2 vertex graph is singular at Q = 1 and it carries exactly the gap *)
tableOf[res2["B"]]
{Factor /@ res2["B"]["stable"][[1]]["polar"], gapHVC[2, x]}

(* mode A and the point expansion: every stable graph is singular *)
tableOf[res2["A"]]
tableOf[res2["point"]]

(* mode B at the level of decorated graphs: the 7 decorations of the vertex carry the polar part, all 20 decorated graphs with an edge
   are regular individually; the structural statements *)
Grid[Prepend[Table[{e[[1]], e[[2, 2]], Factor /@ e[[5]]}, {e, Select[res2["B"]["entries"], #[[1, 1]] == 1 &]}], {"id", "dilaton leaves", "polar part {t^-2, t^-1}"}], Frame -> All]
{Count[Select[res2["B"]["entries"], #[[1, 1]] > 1 &], e_ /; regQ[e[[5]]]], structure[res2["B"]]}

(* the terms of the mode-B edge graphs are simple rational functions of Q and x *)
Table[{k, res2["B"]["stable"][[k]]["genera"], res2["B"]["stable"][[k]]["edges"], Factor[res2["B"]["stable"][[k]]["term"]]}, {k, 2, 7}] // TableForm

(* individual terms have a simple pole at Q = 0 (the canonical coordinates degenerate there) whose residues cancel; the sum at Q = 0
   is the constant-map term F_{2,0}(0) *)
{Factor /@ Table[Residue[res2["B"]["stable"][[k]]["term"], {Q, 0}], {k, 7}], Factor[res2["B"]["total"] /. Q -> 0], Factor[FrefHVC[2, x] /. Q -> 0]}

(* regression against the stored runs of hvc_run.wls *)
Table[If[FileExistsQ["hvc_results_g2_" <> t <> ".m"], Module[{st = Get["hvc_results_g2_" <> t <> ".m"], m = StringTake[t, 1]}, If[m == "p", m = "point"]; {t, And @@ Table[Together[res2[m]["stable"][[k]]["term"] - st["stable"][[k]]["term"]] === 0, {k, 7}]}], {t, "no stored file"}], {t, {"Bm", "Am", "point"}}]

(* ::Section:: Genus 3 *)

(* genus 3: 42 stable graphs, 693 decorated graphs with Hodge vertices (378 with point vertices). The stored results of hvc_run.wls
   (hvc_results_g3_Bm.m, _Am.m, _point.m) are read; set recompute3 = True to recompute everything here (about 5 minutes per mode
   on one kernel, exact in Q and x) *)
recompute3 = False;
res3[mode_] := res3[mode] = Module[{f = "hvc_results_g3_" <> modeTag[mode, -1] <> ".m"}, If[recompute3 || ! FileExistsQ[f], runMode[3, mode, -1], Get[f]]];
Table[{m, Length[res3[m]["entries"]], "total = F30" -> res3[m]["totalOK"], "total polar part = gap" -> res3[m]["gapOK"], Factor /@ res3[m]["polar"]}, {m, {"B", "A", "point"}}] // TableForm

(* mode B per stable graph: only the genus-3 vertex graph is singular at Q = 1 and it carries exactly the gap *)
tableOf[res3["B"]]
{Factor /@ res3["B"]["stable"][[1]]["polar"], gapHVC[3, x]}

(* mode B: the stable graphs with a vertex of genus >= 2 are exactly those containing singular decorated graphs, and the conifold
   pole order of every decorated graph is bounded by Sum_v max(2 g_v - 2, 0) *)
structure[res3["B"]]

(* the genus-3 vertex graph: decomposition of its polar part over its decorated graphs (dilaton-leaf patterns) *)
Grid[Prepend[Table[{e[[1]], e[[2, 2]], Factor /@ e[[5]]}, {e, Select[res3["B"]["entries"], #[[1, 1]] == 1 &]}], {"id", "dilaton leaves", "polar part {t^-4, ..., t^-1}"}], Frame -> All]

(* the other stable graphs with a genus-2 vertex: singular decorated graphs that cancel only in the sum over the decorations of the graph *)
Table[With[{st = res3["B"]["stable"][[k]]}, {k, st["genera"], st["edges"], st["decoratedGraphs"], st["singularDecoratedGraphs"], If[regQ[st["polar"]], "regular", Factor /@ st["polar"]]}], {k, Select[Range[42], Max[res3["B"]["stable"][[#]]["genera"]] >= 2 &]}] // TableForm

(* modes A and point: every stable graph singular *)
tableOf[res3["A"]]
tableOf[res3["point"]]

(* ::Section:: Comparison with local P2 *)

(* ::Text:: The decorated-graph-by-decorated-graph comparison with local P2 in the same frame (genus 2, 3, 4) belongs to the local P2 project since 2026-11-27 (its folder hodge_vertex/compare_conifold/, with frozen copies of the mode-B results hvc_results_g<g>_Bm.m of this project). This notebook no longer reads any file of that project. *)

(* ::Section:: Display *)

(* mode B, genus 3: the terms of the stable graphs with three or four vertices as rational functions of Q and x *)
Table[With[{st = res3["B"]["stable"][[k]]}, {k, st["genera"], st["edges"], Factor[st["term"]]}], {k, Select[Range[42], Length[res3["B"]["stable"][[#]]["genera"]] >= 3 &]}] // TableForm

(* denominators of the mode-B stable-graph terms at genus 3 from FactorList: powers of (1-Q), Q, the discriminant Q(1-x)^2+4x (the two
   idempotents collide at Q_c = -4x/(1-x)^2: every stable graph has such a pole, which cancels only in the sum over all graphs), x, (1+x),
   and whatever else remains (1 = nothing) *)
denPowers[t_] := Module[{fl = FactorList[Denominator[Together[t]]], known = {1 - Q, Q, Q (1 - x)^2 + 4 x, x, 1 + x}, pw, res = {}, hit}, pw = ConstantArray[0, 5];
   Do[If[NumberQ[fl[[i, 1]]], Continue[]]; hit = Position[known, f_ /; NumberQ[Together[fl[[i, 1]]/f]], {1}, 1, Heads -> False]; If[hit === {}, AppendTo[res, fl[[i]]], pw[[hit[[1, 1]]]] += fl[[i, 2]]], {i, Length[fl]}];
   Append[pw, If[res === {}, 1, Times @@ (Power @@@ res)]]];
{Sort[Tally[denPowers /@ res3["B"]["stable"][[All, "term"]]]], "total:", denPowers[res3["B"]["total"]]}

(* residues at Q = 0 of the stable-graph terms cancel; the total at Q = 0 is F_{3,0}(0) *)
{Factor[Total[Table[Residue[res3["B"]["stable"][[k]]["term"], {Q, 0}], {k, 42}]]], Factor[res3["B"]["total"] /. Q -> 0], Factor[FrefHVC[3, x] /. Q -> 0]}

(* ::Section:: Genus 4 (cluster results) *)

(* genus 4 (mode B) is computed on Stanage with cluster_genus4_hodge_conifold/ (379 stable graphs, 27172 decorated graphs, exact in Q and x)
   and copied back into the subfolder stanage_genus4/ of the data directory: hvc_results_g4_Bm.m (per stable graph and per decorated
   graph), hvc_table_g4_Bm.txt, hvc_analysis_g4_Bm.txt. This section reads them if present. *)
g4file = FileNameJoin[{"stanage_genus4", "hvc_results_g4_Bm.m"}];
If[! FileExistsQ[g4file], Print["genus-4 results not present (", g4file, ") -- run cluster_genus4_hodge_conifold/ on Stanage and copy hvc_results_g4_Bm.m back into stanage_genus4/"], res4 = Get[g4file]; Print[Length[res4["stable"]], " stable graphs, ", Length[res4["entries"]], " decorated graphs; total = F40 (closed form): ", res4["totalOK"], "; total polar part = gap: ", res4["gapOK"]]]

(* total polar part {t^-6, ..., t^-1} versus the gap polar part of P_4 Li_{-5}(Q)/(eps1 eps2) *)
If[ValueQ[res4], {Factor /@ res4["polar"], gapHVC[4, x]}]

(* which stable graphs are singular at Q = 1 (expected: only the genus-4 vertex graph, with exactly the gap) *)
If[ValueQ[res4], Module[{sing = Select[Range[Length[res4["stable"]]], ! regQ[res4["stable"][[#]]["polar"]] &]}, {Length[res4["stable"]] - Length[sing], "regular of", Length[res4["stable"]], "; singular:", Table[{k, res4["stable"][[k]]["genera"], res4["stable"][[k]]["edges"], Factor /@ res4["stable"][[k]]["polar"]}, {k, sing}]}]]

(* the structural statements: singular decorated graphs only in the stable graphs with a vertex of genus >= 2, pole order <= Sum_v max(2 g_v - 2, 0) *)
If[ValueQ[res4], structure[res4]]
If[ValueQ[res4], Grid[Prepend[Table[With[{st = res4["stable"][[k]]}, {k, st["genera"], st["edges"], st["decoratedGraphs"], st["singularDecoratedGraphs"]}], {k, Select[Range[Length[res4["stable"]]], res4["stable"][[#]]["singularDecoratedGraphs"] > 0 &]}], {"graph", "genera", "edges", "decorated graphs", "singular decorated graphs"}], Frame -> All]]
