(* hvc_lib.m -- Hodge-vertex graph sum for the resolved conifold x C^2 on the fibrewise torus alpha0 = 0, alpha1 = eps1,
   alpha2 = eps2 (opposite weights on each O + O(-1) pair), genus <= 3 locally (the Hodge machinery is valid to genus 4).
   Omega_Z = R' . (Hodge_w)^{+2} with the SAME triple Hodge class as for local P2: vertex = Int prod_{a=1}^3 Lambda_g^vee(w_a)/w_a^g psi^a
   (GW theory of C^3), mode B: w = {eps1+eps2, -eps1, -eps2}, mode A: w = {-eps1-eps2, eps1, eps2} (A and B differ by the overall sign
   of the weights, Lambda^vee <-> Lambda); R' = R_Z h_w^s with h_w(z) = Exp[-Sum_k B_2k/(2k(2k-1)) z^(2k-1) Sum_a w_a^(1-2k)] and
   s = -1 or +1 (fixed by the requirement that the sum of the decorated graphs is F_{g,0}; for local P2, s = -1 in these conventions).
   R_Z is the conifold R-matrix of resolvedconifold_genus3.nb / rc_lib.m at alpha0 = 0, alpha_i = eps_i, where no QRR normalisation
   constants arise (rmatrixAt0). Everything is exact: Q symbolic, eps1 = x symbolic, eps2 = 1; polar parts at Q = 1 in t = Log Q.
   Sections 1-4 are verbatim from hv_lib4.m (hodge_vertex/, local P2); section 5 is conifold specific. ASCII only. *)

Get["rc_lib.m"];                      (* conifold Frobenius data, quadratic algebra C[a, b], Wick R-matrix, graph factors, graphValueC, polar parts *)
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
