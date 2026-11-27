(* hv_lib4.m -- genus-general Hodge-vertex machinery for the Givental-Teleman graph sum of K_P2 x C^2 on T_F, valid in genus <= 4
   (no closed forms needed): every lambda class is a polynomial in lambda_1 and ch_3(E) (Newton's identities with ch_2 = ch_4 = 0:
   lambda_2 = lambda_1^2/2, lambda_3 = lambda_1^3/6 + 2 ch_3, lambda_4 = lambda_1^4/24 + 2 lambda_1 ch_3), and
   HI[g, k1, k3, a] = Int_{M_{g,n}} lambda_1^k1 ch_3^k3 psi^a is computed by Mumford's formula for ch_1 = lambda_1 and ch_3:
     ch_k = B_{k+1}/(k+1)! ( kappa_k - Sum_i psi_i^k + 1/2 Sum_boundary iota_* Sum_{j=0}^{k-1} (-1)^j psi'^j psi''^(k-1-j) ),
   kappa_k = pi_*(psi_{n+1}^(k+1)) so that <kappa_k X> = <tau_{k+1} X> for pulled-back X; ch_k restricts additively to the
   components of the boundary. Vertex = normalised Hodge class prod_a Lambda_g^vee(w_a)/w_a^g, u[a] = 1/w_a.
   The rest (decorations, R' = R_Z h^(-1) rows, sampling, fits) as in hv_lib3.m. Built on lp2_lib.m, graphs.m, orbits.m. ASCII only. *)

Get["lp2_lib.m"];
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
