(* Plain-text export of localp2CY5_genRmat_genus3.nb (2026-11-27): Input cells as code, text cells as comments; outputs omitted. *)

(* ::Title:: Local P2 x C2: genus 3 free energy with a generalised (twisted) R-matrix *)

(* ::Text:: Extension to genus 3 of localp2CY5_genRmat_v3.nb. There the identity row of the fibrewise R-matrix of K_P2 x C^2 (torus T_F, weights alpha_i = (eps1+eps2)/3) is right-multiplied by a scalar ff[z] = Exp[Sum_k Qrrc[0,k] z^(2k+1)], so that ff[z] ff[-z] = 1 with otherwise arbitrary y-independent coefficients, and the twisted genus-2 free energy F20 is computed symbolically (17 hand-written graphs) and expanded in the conifold frame; the vanishing of its t_c^-1 term fixes Qrrc[0,1] in terms of Qrrc[0,0]. That symbolic route does not scale to genus 3 (42 stable graphs, 840 labelled decorations), so this notebook applies the twist inside the numeric sampling pipeline of localp2CY5_genus3.nb: the graph sum is evaluated at rational sample points (t0 = y^(1/3), x0 = eps1/eps2) with 300-digit arithmetic, keeping C1i = 1/c1, C1p = c1' and the twist coefficients q0, q1, q2 = Qrrc[0,0], Qrrc[0,1], Qrrc[0,2] (for z, z^3, z^5; to order z^6 no further coefficient enters) symbolic, and the result is reconstructed exactly. Evaluate top to bottom (genus 2: seconds; genus 3: about 1.5 hours), or headless with run_genRmat_genus3.wls, which evaluates the Input cells of this notebook in order. Requires cPPu6.mx (R-matrix to z^6, from localp2final_CY5_v16.nb), StableGraph-2.bin, StableGraph-3.bin, v3_F20_twisted_reference.m (the F20xy series, the conifold expansion and the Qrrc[0,1] solution extracted from the output cells of localp2CY5_genRmat_v3.nb), and F20_result.mx, F30_result.mx (untwisted results of localp2CY5_genus3.nb) in the notebook directory. Log: genRmat_run.log. *)

(* ::Text:: Structure used for the reconstruction (every assumption is verified on spare sample points and logged). F_g is a polynomial in XX = y c1'/c1 of degree <= 3g-3 and in q0, q1, q2 of weighted degree <= 3g-3 (weights 1, 3, 5): 5 monomials at genus 2, 14 at genus 3. Since the z^k coefficient of the R-matrix is homogeneous of degree -k in eps and a q-monomial of weight w replaces z^w, the coefficient of a q-monomial of weight w in F_g satisfies: eps1 eps2 * coefficient = symmetric homogeneous polynomial of degree 2g+w in (eps1, eps2), so with eps2 = 1, x = eps1 it is a palindromic polynomial of degree 2g+w divided by x, needing 2g+w+1 values of x (plus 2 checks): 10 eps-ratio samples at genus 2, 15 at genus 3, times 23 values of t. The y-dependence is recovered by generic rational reconstruction verified on >= 5 spare points. *)

(* ::Text:: Checks. (i) The twisted R-matrix is still symplectic (ff[z] ff[-z] = 1) and the edge terms symmetric; (ii) the twisted genus-2 result reproduces v3's F20 (as the double series in eps1+eps2, eps1 eps2 stored by v3, and as the conifold expansion) and v3's Qrrc[0,1] relation; (iii) at q = 0 both genera reproduce the untwisted F20, F30 of localp2CY5_genus3.nb; (iv) the genus-3 sum at q = 0 was also checked to agree with the untwisted graph sum to 10^-295 at a sample point when this pipeline was set up. *)

(* ::Section:: Core routines *)

(* ::Text:: Verbatim copy of the core of localp2CY5_genus3.nb (conventions: Dilaton[i,n] = [z^(n-1)] (-P_0i(-z)), Leaf[i,n,a] = [z^n] P_ai(-z), EdgeTerm[i,j,m,n] = [w^m z^n] (delta_ij - sum_ab P_ai(-w) eta^{-1}_ab P_bj(-z))/(w+z), sqrtDelta[i] = Sqrt[eps1 eps2/e3] (3 e3 + e2 LLL[i]); a graph contributes 1/|Aut| times the sum over all labelled psi-assignments, dilaton leaves and 3^V colourings with the Witten-Kontsevich vertex numbers). Its buildRmatrix is replaced by the twisted one in the next cell. *)

(* Part B core: fibrewise (T_F) Givental graph sum for K_P2 x C^2 by exact-structure sampling.
   Load after defining cPP[0, i_][k] for k = 0..KMAX (rational functions of y, LLL[i][y], eps).
   Everything is evaluated at numeric sample points (t0 = y^(1/3), x0 = eps1/eps2) in high precision,
   keeping C1i = 1/c1, C1p = c1' symbolic so that the dependence on XX = y c1'/c1 is separated. *)

PREC = 300;
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

(* ::Section:: Twisted driver *)

(* Twisted driver: F_{g,0} on T_F as a polynomial in the twist coefficients q0, q1, q2, by sampling + reconstruction.
   Structure used for the reconstruction (verified on spare sample points): the coefficient of a q-monomial of
   weight w in F_g is a polynomial in XX of degree <= 3g-3 whose coefficients are rational in y, and
   eps1 eps2 * (that coefficient) is a symmetric homogeneous polynomial of degree 2g + w in (eps1, eps2)
   (each z^k of the R-matrix carries eps^-k, and a q-monomial of weight w replaces z^w).
   Output: genRmat_F<g>0_result.mx with Fgtw (in y, XX, eps, q), Fgredtw (in LL = (1+27y)^(-1/3), XX, eps, q), monstw;
   genRmat_F<g>0_result.txt with Fgred in the Qrrc[0,k] notation. Returns Fg. *)
runGenusTw[gg_, KMAX_, cppfile_] := Module[{wmax, mons, prepared, tsamples, xsamples, NT, NX, sample, aCoef, Fg, Fgred, failed, tstart = AbsoluteTime[]},
log[m_] := Module[{s = OpenAppend["genRmat_run.log"]}, WriteString[s, DateString[{"Hour", ":", "Minute", ":", "Second"}], "  [g=", gg, "] ", m, "\n"]; Close[s]; Print[m]];
Get[cppfile]; Do[cPP[0, 0][k] = cPPL[[k + 1]], {k, 0, KMAX}];
buildRmatrix[KMAX];
wmax = 3 gg - 3; mons = qmonos[wmax];
log["genus " <> ToString[gg] <> ", KMAX = " <> ToString[KMAX] <> ": " <> ToString[Length[mons]] <> " q-monomials of weight <= " <> ToString[wmax] <> ": " <> ToString[qmono /@ mons[[All, 1]], InputForm]];

(* R-matrix checks at a sample point: symplecticity sum_i P_bi(-z) etainv_ac P_ci(z) = delta_ab to O(z^KMAX) (holds for the twisted matrix because ff[z] ff[-z] = 1), edge symmetry *)
Module[{data = evalPoint[1/3, 5, KMAX], Pn, symp, edgesym},
  Pn[a_, i_] := data["Pnum"][a, i];
  symp = Max[Table[maxabs[Sum[Sum[(-1)^m Pn[b, i][[m + 1]] etainvN[[a + 1, c + 1]] Pn[c, i][[n - m + 1]], {m, 0, n}], {i, 0, 2}, {c, 0, 2}] - KroneckerDelta[a, b] KroneckerDelta[n, 0]], {a, 0, 2}, {b, 0, 2}, {n, 0, KMAX}]];
  edgesym = Max[Table[maxabs[data["Edge"][i, j][m, n] - data["Edge"][j, i][n, m]], {i, 0, 2}, {j, 0, 2}, {m, 0, KMAX - 1}, {n, 0, KMAX - 1 - m}]];
  log["twisted R-matrix symplecticity residual (should be ~0): " <> ToString[N[symp, 5], InputForm] <> ";  edge-term symmetry residual: " <> ToString[N[edgesym, 5], InputForm]]];
prepared = prepareGraph /@ graphsWithLegs[gg, 0];
log[ToString[Length[prepared]] <> " graphs, labelled decorations " <> ToString[Total[Length[#["decs"]] & /@ prepared]]];

tsamples = {1/2, 1/3, 2/3, 1/4, 3/4, 1/5, 2/5, 3/5, 4/5, 1/6, 5/6, 1/7, 2/7, 3/7, 4/7, 5/7, 6/7, 1/8, 3/8, 5/8, 7/8, 1/9, 2/9};
xsamples = Table[Prime[r + 1], {r, 1, 2 gg + wmax + 3}];   (* eps1/eps2 values 3,5,7,...: 2g+w+1 needed for weight w, +2 checks *)
NT = Length[tsamples]; NX = Length[xsamples];
log[ToString[NT] <> " t-samples x " <> ToString[NX] <> " eps-ratio samples"];
ratl[v_] := Rationalize[Re[v], 10^-(PREC - 60)];
sample = <||>;
Do[Module[{tt = AbsoluteTime[], data, val, rules, asc, bad},
   data = evalPoint[t0, x0, KMAX];
   val = graphSum[prepared, data];
   rules = CoefficientRules[val, vars]; asc = Association[rules];
   If[Max[Abs[Im[rules[[All, 2]]]]] > 10^-(PREC - 80), log["WARNING: imaginary part at " <> ToString[{t0, x0}]]];
   (* structure: only (C1i C1p)^k, k <= 3g-3, and q-monomials of weight <= 3g-3 *)
   bad = Select[rules, Abs[#[[2]]] > 10^-(PREC - 80) && (#[[1, 1]] != #[[1, 2]] || #[[1, 1]] > wmax || #[[1, 3]] + 3 #[[1, 4]] + 5 #[[1, 5]] > wmax) &][[All, 1]];
   If[bad =!= {}, log["WARNING: unexpected monomials (exponents of C1i, C1p, q0, q1, q2) " <> ToString[bad] <> " at " <> ToString[{t0, x0}]]];
   (* sample[{t0,x0}][[k+1, m]] = coefficient of XX^k * (q-monomial m) *)
   sample[{t0, x0}] = Table[ratl[Lookup[asc, Key[Join[{k, k}, mons[[m, 1]]]], 0]/t0^(3 k)], {k, 0, wmax}, {m, Length[mons]}];
   If[x0 == xsamples[[-1]], log["t0 = " <> ToString[t0] <> " done (" <> ToString[NX] <> " eps-ratios, " <> ToString[Round[AbsoluteTime[] - tt, 0.1]] <> " s for the last one; elapsed " <> ToString[Round[(AbsoluteTime[] - tstart)/60, 0.1]] <> " min)"]];
   ], {t0, tsamples}, {x0, xsamples}];
samplestw = sample; tsamplestw = tsamples; xsamplestw = xsamples; monstw = mons;   (* global names, so that Get restores them *)
DumpSave["genRmat_F" <> ToString[gg] <> "0_samples.mx", {samplestw, tsamplestw, xsamplestw, monstw}];
log["all samples done"];

(* eps-structure: for the coefficient of XX^k and q-monomial m of weight w:  x * coeff(x) = sum_{j=0}^{2g+w} a_j x^j (eps2 = 1), palindromic *)
aCoef = <||>;
Do[Module[{w = mons[[m, 2]], DD, A, rhs, sol, chk},
   DD = 2 gg + w;
   A = Table[x0^j, {x0, xsamples[[1 ;; DD + 1]]}, {j, 0, DD}];
   rhs = Table[x0 sample[{t0, x0}][[k + 1, m]], {x0, xsamples[[1 ;; DD + 1]]}];
   sol = LinearSolve[A, rhs];
   chk = Table[Sum[sol[[j + 1]] x0^j, {j, 0, DD}] - x0 sample[{t0, x0}][[k + 1, m]], {x0, xsamples[[DD + 2 ;;]]}];
   If[DeleteCases[chk, 0] =!= {}, log["WARNING: eps-polynomial of degree " <> ToString[DD] <> " fails at t=" <> ToString[t0] <> " k=" <> ToString[k] <> " q-monomial " <> ToString[qmono[mons[[m, 1]]], InputForm] <> ": " <> ToString[N[Max[Abs[chk]], 3]]]];
   If[sol =!= Reverse[sol], log["WARNING: eps-polynomial not palindromic at t=" <> ToString[t0] <> " k=" <> ToString[k] <> " q-monomial " <> ToString[qmono[mons[[m, 1]]], InputForm]]];
   Do[aCoef[{t0, k, m, j}] = sol[[j + 1]], {j, 0, DD}]], {t0, tsamples}, {k, 0, wmax}, {m, Length[mons]}];
log["eps structure done"];

(* y-dependence: reconstruct a_{k,m,j}(y) as rational functions of y = t^3 *)
Fg = 0; failed = {};
Do[Module[{w = mons[[m, 2]], pts, r},
   pts = Table[{t0^3, aCoef[{t0, k, m, j}]}, {t0, tsamples}];
   r = reconstruct[pts] /. u -> y;
   If[r === $Failed, AppendTo[failed, {k, m, j}], Fg += r XX^k \[Epsilon][1]^j \[Epsilon][2]^(2 gg + w - j) qmono[mons[[m, 1]]]]], {k, 0, wmax}, {m, Length[mons]}, {j, 0, 2 gg + mons[[m, 2]]}];
If[failed =!= {}, log["RECONSTRUCTION FAILED for (k, m, j) = " <> ToString[failed]]];
Fg = Fg/(\[Epsilon][1] \[Epsilon][2]);
Fgred = Collect[Expand[Together[Fg/(\[Epsilon][1] \[Epsilon][2]) /. y -> (1 - LL^3)/(27 LL^3)]], Join[{XX}, qsyms], Factor];
Fgtw = Fg; Fgredtw = Fgred; monstw = mons; DumpSave["genRmat_F" <> ToString[gg] <> "0_result.mx", {Fgtw, Fgredtw, monstw}];   (* global names: Fgtw (y, XX, eps, q), Fgredtw (LL, XX, eps, q) *)
Export["genRmat_F" <> ToString[gg] <> "0_result.txt", ToString[InputForm[Fgred /. qrule]], "Text"];
log["F" <> ToString[gg] <> "0 reconstructed in " <> ToString[Round[(AbsoluteTime[] - tstart)/60, 0.1]] <> " min; denominators: " <> ToString[Union[Cases[Fg, Power[b_, e_?Negative] :> b, Infinity]]] <> "; LeafCount " <> ToString[LeafCount[Fg]]];
Fg]

(* ::Section:: Genus 2: regression against the symbolic twisted computation of localp2CY5_genRmat_v3.nb *)

If[$Notebooks, SetDirectory[NotebookDirectory[]]];
F20tw = runGenusTw[2, 3, "cPPu6.mx"];
{refdesc, F20xyRef, F20cfRef, Q01Ref, cfmmRef} = Get["v3_F20_twisted_reference.m"];
e1 = \[Epsilon][1]; e2 = \[Epsilon][2];
(* 1. against v3's F20xy: F20 with eps1, eps2 -> (xx -+ Sqrt[xx^2 - 4 yy])/2 expanded in xx = eps1 + eps2, yy = eps1 eps2 *)
replxy = First[Solve[{xx == e1 + e2, yy == e1 e2}, {e1, e2}]];
F20c1 = F20tw /. XX -> y Derivative[1][c1][y]/c1[y] /. qrule;
diffxy = Simplify[Normal[Series[F20c1 /. replxy, {xx, 0, 8}, {yy, 0, 5}] - F20xyRef]];
log["genus 2 twisted: difference from the v3 F20xy double series (should be 0): " <> ToString[diffxy, InputForm]];
(* 2. at q = 0 against the untwisted F20 of localp2CY5_genus3.nb (itself checked against v15) *)
Module[{}, Get["F20_result.mx"]; F20untw = Fg];
log["genus 2 twisted at q = 0 minus untwisted F20_result.mx (should be 0): " <> ToString[Simplify[(F20tw /. Thread[qsyms -> 0]) - F20untw], InputForm]];
Module[{}, Get["genRmat_F20_result.mx"]; F20redtw = Fgredtw /. qrule];
Print["F20red (twisted; F20 = eps1 eps2 F20red, LL = (1+27y)^(-1/3)) = ", InputForm[F20redtw]];

(* ::Section:: Genus 3 *)

F30tw = runGenusTw[3, 6, "cPPu6.mx"];
Module[{}, Get["F30_result.mx"]; F30untw = Fg];
log["genus 3 twisted at q = 0 minus untwisted F30_result.mx (should be 0): " <> ToString[Simplify[(F30tw /. Thread[qsyms -> 0]) - F30untw], InputForm]];

(* ::Section:: Conifold frame and gap conditions *)

(* Conifold frame (as in localp2CY5_genus3.nb): t_c from the analytic solution of theta^2 + 3y(3theta+1)(3theta+2) at Delta = 1 + 27 y = 0,
   C1cf = theta t_c, and the expansion of F_g /. c1 -> C1cf in t_c to O(t_c^1). Then the gap conditions: v3 fixes Qrrc[0,1] from the
   vanishing t_c^-1 term at genus 2; here the t_c^-3, t_c^-2, t_c^-1 terms at genus 3 (each linear in Qrrc[0,2]) are examined. *)
GMAX = 3; NC = 2 GMAX + 4;
L2b[f_] := Module[{th = y D[#, y] &}, (th[th[f]] + 3 y (9 th[th[f]] + 9 th[f] + 2 f))];
wc = Sum[a[n] Del^n, {n, 0, NC}] /. Del -> 1 + 27 y;
ceq = CoefficientList[Normal[Series[L2b[wc] /. y -> (Del - 1)/27, {Del, 0, NC - 1}]], Del];
csol = Solve[Thread[ceq == 0] /. a[0] -> 1, Table[a[n], {n, 1, NC}]];
wcDel = Sum[a[n] Del^n, {n, 0, NC}] /. First[csol] /. a[0] -> 1;
tc = Module[{ser = Normal[Series[wcDel/(Del - 1), {Del, 0, NC}]]}, Sum[Coefficient[ser, Del, n] Del^(n + 1)/(n + 1), {n, 0, NC}]];
tc = Normal[Series[tc/SeriesCoefficient[tc, {Del, 0, 1}], {Del, 0, NC}]];
log["conifold flat coordinate t_c(Delta) = " <> ToString[tc, InputForm] <> "; agrees with v3's cfmm to Delta^4: " <> ToString[Simplify[Normal[Series[tc, {Del, 0, 4}]] - Normal[cfmmRef /. \[CapitalDelta] -> Del]] === 0]];
invcf = InverseSeries[Series[tc, {Del, 0, NC}], tcv];
C1cf = Normal[y 27 D[tc, Del]] /. Del -> 1 + 27 y;
cfexp[f_] := Series[f /. XX -> y D[C1cf, y]/C1cf /. cc -> C1cf /. y -> (Del - 1)/27 /. Del -> invcf, {tcv, 0, 1}];

(* genus 2: against v3's twisted conifold expansion, and its Qrrc[0,1] solution *)
cf[2] = cfexp[F20tw /. qrule];
log["genus 2 twisted conifold expansion minus v3 (should be 0): " <> ToString[Simplify[Normal[cf[2] - F20cfRef]], InputForm]];
q1sol = Solve[SeriesCoefficient[cf[2], -1] == 0, Qrrc[0, 1]];
log["genus 2 gap (t_c^-1 = 0) solved for Qrrc[0,1]: " <> ToString[FullSimplify[q1sol], InputForm] <> ";  minus v3's solution: " <> ToString[Simplify[(Qrrc[0, 1] /. First[q1sol]) - (Qrrc[0, 1] /. First[Q01Ref])], InputForm]];

(* genus 3 *)
cf[3] = cfexp[F30tw /. qrule];
Do[log["genus 3 twisted conifold expansion, coefficient of t_c^" <> ToString[n] <> ": " <> ToString[Factor[SeriesCoefficient[cf[3], n]], InputForm]], {n, -4, 1}];
Do[cs[n] = Factor[SeriesCoefficient[cf[3], n] /. First[q1sol]], {n, -4, 1}];
Do[log["  with Qrrc[0,1] from the genus-2 gap, coefficient of t_c^" <> ToString[n] <> ": " <> ToString[cs[n], InputForm]], {n, -4, 1}];
Do[q2sol[n] = Solve[cs[n] == 0, Qrrc[0, 2]]; log["genus 3 gap condition t_c^" <> ToString[n] <> " = 0 solved for Qrrc[0,2]: " <> ToString[FullSimplify[q2sol[n]], InputForm]], {n, -3, -1}];
log["consistency of the three genus-3 gap conditions (differences of the Qrrc[0,2] solutions from t_c^-3, t_c^-2 against t_c^-1; 0 if consistent): " <> ToString[Simplify[Table[(Qrrc[0, 2] /. First[q2sol[n]]) - (Qrrc[0, 2] /. First[q2sol[-1]]), {n, -3, -2}]], InputForm]];
log["genus 3 leading coefficient t_c^-4 with both gap solutions: " <> ToString[Factor[cs[-4] /. First[q2sol[-1]]], InputForm]];
log["genus 3 regular part (t_c^0, t_c^1) with both gap solutions: " <> ToString[Factor /@ {cs[0] /. First[q2sol[-1]], cs[1] /. First[q2sol[-1]]}, InputForm]];
DumpSave["genRmat_conifold.mx", {tc, invcf, C1cf, cf, q1sol, q2sol}];

(* Observations on the gap-constrained expansions (Qrrc[0,1], Qrrc[0,2] as fixed by the gap conditions):
   (a) every coefficient of the conifold expansions computed (t_c^-2 .. t_c^1 at genus 2, t_c^-4 .. t_c^1 at genus 3)
       equals the untwisted coefficient with (eps1 + eps2)^2 -> (eps1 + eps2)^2 + 12 Qrrc[0,0] eps1 eps2 (eps1 + eps2) at fixed eps1 eps2
       (written in sig = eps1 + eps2, pp = eps1 eps2 the untwisted coefficients are even in sig, so the substitution is sig -> Sqrt[sig^2 + 12 Qrrc[0,0] pp sig]);
   (b) the pure-Qrrc[0,0] parts of the gap solutions are the quantum Riemann-Roch coefficients of one extra weight w = -1/(12 Qrrc[0,0]):
       Qrrc[0,k] ~ -B_{2k+2}/((2k+2)(2k+1)) w^-(2k+1), i.e. -24/5 Qrrc[0,0]^3 and 6912/35 Qrrc[0,0]^5. *)
toSP[f_] := Module[{g = Together[f]}, SymmetricReduction[Numerator[g], {\[Epsilon][1], \[Epsilon][2]}, {sig, pp}][[1]]/SymmetricReduction[Denominator[g], {\[Epsilon][1], \[Epsilon][2]}, {sig, pp}][[1]]];
gapsub = Join[First[q1sol], First[q2sol[-1]]];
subcheck[c_] := Module[{tw = Together[c /. gapsub], un}, un = tw /. Qrrc[0, 0] -> 0; Simplify[(toSP[un] /. sig -> Sqrt[sig^2 + 12 Qrrc[0, 0] pp sig]) - toSP[tw]]];
log["(a) twisted-with-gap minus [untwisted with (eps1+eps2)^2 -> (eps1+eps2)^2 + 12 Qrrc[0,0] eps1 eps2 (eps1+eps2)], genus 2 t_c^-2..t_c^1: " <> ToString[Table[subcheck[SeriesCoefficient[cf[2], n]], {n, -2, 1}], InputForm] <> "; genus 3 t_c^-4..t_c^1: " <> ToString[Table[subcheck[SeriesCoefficient[cf[3], n]], {n, -4, 1}], InputForm]];
qrr[k_] := -BernoulliB[2 k + 2]/((2 k + 2) (2 k + 1)) (-12 Qrrc[0, 0])^(2 k + 1);
log["(b) pure-Qrrc[0,0] part of the gap solutions minus the QRR coefficients of the weight -1/(12 Qrrc[0,0]): " <> ToString[Simplify[{Coefficient[Qrrc[0, 1] /. gapsub, Qrrc[0, 0], 3] Qrrc[0, 0]^3 - qrr[1], Coefficient[Qrrc[0, 2] /. gapsub, Qrrc[0, 0], 5] Qrrc[0, 0]^5 - qrr[2]}], InputForm]];
log["gap solutions collected in Qrrc[0,0]: Qrrc[0,1] = " <> ToString[Collect[Together[Qrrc[0, 1] /. gapsub], Qrrc[0, 0], Factor], InputForm] <> ";  Qrrc[0,2] = " <> ToString[Collect[Together[Qrrc[0, 2] /. gapsub], Qrrc[0, 0], Factor], InputForm]];

(* ::Section:: Results *)

(* ::Text:: Results of the headless run of 7 September 2026 (genRmat_run.log; genus 2 in 0.3 min, genus 3 in 78.6 min for 23 x 15 sample points at about 13 s each, no warnings; symplecticity residual of the twisted R-matrix 10^-296). All regression checks returned exactly 0: the twisted genus-2 F20 against v3 (F20xy double series and conifold expansion), the Qrrc[0,1] relation against v3, and both twisted free energies at q = 0 against the untwisted F20_result.mx, F30_result.mx.

Files: genRmat_F20_result.mx, genRmat_F30_result.mx (Fgtw in y, XX, eps, q0..q2; Fgredtw in LL = (1+27y)^(-1/3); monstw), genRmat_F30_result.txt (Fgred with F30 = eps1 eps2 Fgred in the Qrrc notation), genRmat_F30_samples.mx (exact samples), genRmat_conifold.mx (t_c, its inverse, C1cf, the expansions cf[2], cf[3], q1sol, q2sol), genRmat_gap_results.m.

Genus 3 twisted conifold expansion: the coefficients of t_c^-4 .. t_c^1 are polynomials in Qrrc[0,0], Qrrc[0,1], Qrrc[0,2] (weighted degree <= 6) with coefficients homogeneous in eps (see the log). With Qrrc[0,1] fixed by the genus-2 gap condition, each of the three genus-3 gap coefficients (t_c^-3, t_c^-2, t_c^-1) factorises as (eps1+eps2)^3 times a q0-dependent polynomial in eps times one common factor linear in Qrrc[0,2], and the three conditions are mutually consistent: they give the same solution

Qrrc[0,1] = -24/5 Qrrc[0,0]^3 - Qrrc[0,0] (eps1^2 + eps2^2)/(10 eps1^2 eps2^2) - 6/5 Qrrc[0,0]^2 (eps1^2 + eps1 eps2 + eps2^2)/(eps1 eps2 (eps1 + eps2)),

Qrrc[0,2] = 6912/35 Qrrc[0,0]^5 + 576/7 Qrrc[0,0]^4 (eps1^2 + eps1 eps2 + eps2^2)/(eps1 eps2 (eps1+eps2)) + 48/7 Qrrc[0,0]^3 (2 eps1^4 + 4 eps1^3 eps2 + 5 eps1^2 eps2^2 + 4 eps1 eps2^3 + 2 eps2^4)/(eps1^2 eps2^2 (eps1+eps2)^2) + 4/7 Qrrc[0,0]^2 (2 eps1^4 + 2 eps1^3 eps2 + 3 eps1^2 eps2^2 + 2 eps1 eps2^3 + 2 eps2^4)/(eps1^3 eps2^3 (eps1+eps2)) + Qrrc[0,0] (eps1^2 - eps1 eps2 + eps2^2)(eps1^2 + eps1 eps2 + eps2^2)/(21 eps1^4 eps2^4).

With both relations imposed the expansion has the gap form c_-4 t_c^-4 + c_0 + c_1 t_c + O(t_c^2) with
c_-4 = -(eps1^2 + eps2^2 + 12 Qrrc[0,0] eps1 eps2 (eps1+eps2)) (31 (eps1^2 + eps2^2 + 12 Qrrc[0,0] eps1 eps2 (eps1+eps2))^2 - 44 eps1^2 eps2^2)/(17920 eps1 eps2),
and at genus 2 the leading term is -(7 (eps1^2 + eps2^2 + 12 Qrrc[0,0] eps1 eps2 (eps1+eps2))^2 - 4 eps1^2 eps2^2)/(1920 eps1 eps2).

Two observations (verified symbolically in the last Input cell): (a) with the gap relations imposed, every coefficient computed (t_c^-2 .. t_c^1 at genus 2, t_c^-4 .. t_c^1 at genus 3) equals the untwisted coefficient with (eps1+eps2)^2 replaced by (eps1+eps2)^2 + 12 Qrrc[0,0] eps1 eps2 (eps1+eps2) at fixed eps1 eps2; (b) the pure-Qrrc[0,0] parts of the two relations, -24/5 Qrrc[0,0]^3 and 6912/35 Qrrc[0,0]^5, are the quantum Riemann-Roch coefficients -B_{2k+2}/((2k+2)(2k+1)) w^-(2k+1) of a single extra weight w = -1/(12 Qrrc[0,0]); the eps-dependent terms are cross terms with the weights of the geometry. *)
