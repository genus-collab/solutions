(* hv_lib3.m -- genus-general (g <= 3) Hodge-vertex machinery for the Givental-Teleman graph sum of K_P2 x C^2 on T_F.
   Omega_Z = R' . (Hodge_w)^{+3}, R' = R_Z h_w^(-1), h_w(z) = Exp[-Sum_k B_2k/(2k(2k-1)) z^(2k-1) Sum_a w_a^(1-2k)].
   Vertex = normalised Hodge class prod_a Lambda_g^vee(w_a)/w_a^g expanded in lambda monomials (Mumford relations) with
   coefficients polynomial in the inverse weights u[a] = 1/w_a; the lambda integrals with psi insertions come from
   (i) Int lambda_1^k psi^a by k-fold Mumford's formula lambda_1 = (kappa_1 - Sum psi_i + delta)/12,
   (ii) the lambda_g formula, (iii) lambda_g lambda_{g-1} (closed form for a_i >= 1, zeros by the string equation),
   (iv) lambda_g lambda_{g-1} lambda_{g-2} in genus 3 (dilaton from Int_{M_3} = |B_4 B_6|/(2 4 6 4!)).
   Built on lp2_lib.m (R-matrix rows, colour sums, fits) and graphs.m / orbits.m.  ASCII only. *)

Get["lp2_lib.m"];
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
