(* hv_lib.m -- Hodge-vertex form of the Givental-Teleman graph sum for K_P2 x C^2 on the fibrewise torus T_F.
   The semisimple CohFT is written as  Omega_Z = R' . (Hodge^{+3})  instead of  R_Z . (point^{+3}):
   every vertex carries the equivariant GW theory of C^3 (Hodge CohFT  prod_a Lambda_g^vee(w_a)/w_a) with the three
   fibre-direction weights w = {w1, w2, w3} (the same at the three idempotents), and the R-matrix is R' = R_Z h^(-1) with
   h(z) = Exp[-Sum_k B_2k/(2k(2k-1)) z^(2k-1) Sum_a w_a^(1-2k)] the scalar Bernoulli series of the weights, so that
   h . point = Hodge (Mumford / Givental QRR) and R' . Hodge = (R' h) . point = R_Z . point by the group-action property.
   Built on lp2_lib.m (R-matrix rows, colour sums, fits) and graphs.m / orbits.m.  ASCII only. *)

Get["lp2_lib.m"];
e1s = \[Epsilon][1]; e2s = \[Epsilon][2];

(* ============================== 1. Hodge integrals in genus <= 2 ============================== *)
multin[a_List] := Total[a]!/Times @@ (a!);
(* lambda_1 on M_{g,n} by Mumford's formula lambda_1 = (kappa_1 - Sum psi_i + delta)/12, kappa_1 = pi_*(psi_{n+1}^2):
   <kappa_1 prod tau_a> = <tau_2 prod tau_a>; boundary = 1/2 irreducible (genus g-1, two extra points) + every separating
   divisor once (= 1/2 the sum over ordered (h,S)); unstable components give wk = 0 *)
L1[g_, a_List] := Module[{n = Length[a], tot, Sc},
   If[g < 1 || Total[a] != 3 g - 3 + n - 1 || 2 g - 2 + n <= 0, Return[0]];
   tot = wk[g, Append[a, 2]] - Sum[wk[g, ReplacePart[a, i -> a[[i]] + 1]], {i, n}] + wk[g - 1, Join[a, {0, 0}]]/2;
   tot += 1/2 Sum[Sc = Complement[Range[n], S]; wk[h, Append[a[[S]], 0]] wk[g - h, Append[a[[Sc]], 0]], {h, 0, g}, {S, Subsets[Range[n]]}];
   tot/12];
(* closed forms: lambda_g formula (Faber-Pandharipande) and lambda_g lambda_{g-1} (Getzler-Pandharipande) *)
L1closed1[a_List] := If[Length[a] >= 1 && Total[a] == Length[a] - 1, multin[a]/24, 0];                 (* lambda_1 on M_{1,n} *)
L2g2[a_List] := If[Total[a] == Length[a] + 1, 7 multin[a]/5760, 0];                                      (* lambda_2 on M_{2,n} *)
(* lambda_1 lambda_2 on M_{2,n} (psi degree n): the closed form (2g-3+n)! |B_2g| / (2^(2g-1) (2g)! prod (2a_i-1)!!) holds for
   a_i >= 1 (here: all a_i = 1, value (n+1)!/5760); zeros are removed by the string equation *)
L12g2[a_List] := Module[{n = Length[a], z, rest},
   If[Total[a] != n, Return[0]];
   z = Position[a, 0, {1}, Heads -> False];
   If[z =!= {}, rest = Delete[a, z[[1]]]; Return[Sum[If[rest[[j]] >= 1, L12g2[ReplacePart[rest, j -> rest[[j]] - 1]], 0], {j, Length[rest]}]]];
   (n + 1)!/5760];
(* lambda_1^2 on M_{2,n} by Mumford's formula applied to lambda_1 . lambda_1 (lambda_1 restricts to the components); must equal 2 lambda_2 *)
L1sq[a_List] := Module[{n = Length[a], tot, Sc},
   If[Total[a] != n + 1, Return[0]];
   tot = L1[2, Append[a, 2]] - Sum[L1[2, ReplacePart[a, i -> a[[i]] + 1]], {i, n}] + L1[1, Join[a, {0, 0}]]/2;
   tot += 1/2 Sum[Sc = Complement[Range[n], S];
      Which[h == 1, L1[1, Append[a[[S]], 0]] wk[1, Append[a[[Sc]], 0]] + wk[1, Append[a[[S]], 0]] L1[1, Append[a[[Sc]], 0]],
       h == 0, wk[0, Append[a[[S]], 0]] L1[2, Append[a[[Sc]], 0]],
       h == 2, L1[2, Append[a[[S]], 0]] wk[0, Append[a[[Sc]], 0]]], {h, 0, 2}, {S, Subsets[Range[n]]}];
   tot/12];

(* normalised Hodge class  prod_a Lambda_g^vee(w_a)/w_a^g  (degree-0 part 1) integrated against psi^a, in terms of the
   elementary symmetric functions he1, he2, he3 of the inverse weights 1/w_a (Mumford relations lambda_1^2 = 2 lambda_2,
   lambda_2^2 = 0 used in genus 2):  1 - he1 lambda_1 + he1^2 lambda_2 - (he1 he2 - he3) lambda_1 lambda_2 *)
hodgeNorm[0, a_List] := wk[0, a];
hodgeNorm[1, a_List] := wk[1, a] - he1 L1[1, a];
hodgeNorm[2, a_List] := wk[2, a] - he1 L1[2, a] + he1^2 L2g2[a] - (he1 he2 - he3) L12g2[a];
heSub[w_List] := {he1 -> Together[Total[1/w]], he2 -> Together[Total[Times @@@ Subsets[1/w, {2}]]], he3 -> Together[Times @@ (1/w)]};

(* the scalar Bernoulli series h(z)^s of the weights, coefficient list to z^KMAX *)
hCoeffs[w_List, s_, KMAX_] := Module[{p, ex},
   p[m_] := Total[w^(-m)];
   ex = -Sum[BernoulliB[2 k]/(2 k (2 k - 1)) z^(2 k - 1) p[2 k - 1], {k, 1, Ceiling[(KMAX + 1)/2]}];
   Table[Together[SeriesCoefficient[Series[Exp[s ex], {z, 0, KMAX}], j]], {j, 0, KMAX}]];

(* ============================== 2. consistency checks of the Hodge integrals ============================== *)
hodgeChecks[] := Module[{ok = True, c, rep = {}},
   (* known values *)
   c = {L1[1, {0}] == 1/24, L1[2, {3}] == 1/480, L2g2[{2}] == 7/5760, L12g2[{1}] == 1/2880, L12g2[{}] == 1/5760};
   AppendTo[rep, "known values (lambda_1 on M_{1,1}, lambda_1 psi^3, lambda_2 psi^2, lambda_1 lambda_2 psi on M_{2,1}, lambda_1 lambda_2 on M_2): " <> ToString[c]]; ok = ok && And @@ c;
   (* Mumford in genus 1 against the closed form *)
   c = And @@ Table[L1[1, a] == L1closed1[a], {a, Join @@ Table[Select[Tuples[Range[0, n], n], Total[#] == n - 1 &], {n, 1, 4}]}];
   AppendTo[rep, "genus-1 Mumford = (1/24) multinomial: " <> ToString[c]]; ok = ok && c;
   (* lambda_1^2 = 2 lambda_2 in genus 2 *)
   c = And @@ Table[L1sq[a] == 2 L2g2[a], {a, Join @@ Table[Select[Tuples[Range[0, n + 1], n], Total[#] == n + 1 &], {n, 1, 3}]}];
   AppendTo[rep, "genus-2 Mumford: lambda_1^2 = 2 lambda_2 on M_{2,n}, n <= 3: " <> ToString[c]]; ok = ok && c;
   (* string and dilaton equations for the Hodge correlators (lambda classes are pulled back from M_{g,n}) *)
   c = And @@ Flatten[Table[{Together[hodgeNorm[g, Append[a, 0]] - Sum[If[a[[i]] >= 1, hodgeNorm[g, ReplacePart[a, i -> a[[i]] - 1]], 0], {i, Length[a]}]] === 0,
        Together[hodgeNorm[g, Append[a, 1]] - (2 g - 2 + Length[a]) hodgeNorm[g, a]] === 0},
      {g, 1, 2}, {a, Join @@ Table[Tuples[Range[0, 3 g - 3 + n], n], {n, 1, 3}]}]];
   AppendTo[rep, "string and dilaton equations for Hodge correlators, genus 1-2, n <= 3: " <> ToString[c]]; ok = ok && c;
   {ok, rep}];

(* R-action of the scalar h on the point theory in the PPZ conventions (legs R^-1(psi), edges (1 - R^-1(psi')R^-1(psi''))/(psi'+psi''),
   translation z(1 - R^-1(z))), summed over the stable graphs of genus 2 with point (WK) vertices: must give the normalised
   Hodge class integrated over M_2, i.e. hodgeNorm[2, {}] = -(he1 he2 - he3)/5760, for arbitrary weights *)
ppzScalarG2[w_List] := Module[{KMAX = 3, hc, A, Bc, E1, Dil1, graphs, tot = 0, hp, vfac, efac},
   hc = hCoeffs[w, 1, KMAX]; A = hCoeffs[w, -1, KMAX];
   Bc[i_, j_] := If[i == 0 && j == 0, 0, -A[[i + 1]] A[[j + 1]]];
   E1[m_, n_] := Together[Sum[(-1)^k Bc[m + 1 + k, n - k], {k, 0, n}]];
   Dil1[d_] := -A[[d]];
   graphs = prepareGraph /@ graphsWithLegs[2, 0];
   Do[Do[hp = ConstantArray[0, Length[gr["H"]]];
     Do[Do[hp[[gr["hsAt"][[v, r]]]] = dec[[v, 1, r]], {r, Length[gr["hsAt"][[v]]]}], {v, gr["V"]}];
     vfac = Product[dec[[v, 4]] Product[Dil1[d], {d, dec[[v, 3]]}], {v, gr["V"]}];
     efac = Product[E1[hp[[2 e - 1]], hp[[2 e]]], {e, Length[gr["edges"]]}];
     tot += vfac efac/gr["aut"], {dec, gr["decs"]}], {gr, graphs}];
   Together[tot]];

(* ============================== 3. decorated graphs with Hodge vertices ============================== *)
(* as vertexOptions, but the psi degree may be smaller than the dimension (the lambda classes make up the difference) *)
vertexOptionsH[g_, hs_, m_] := Module[{res = {}, dim, nd, dils, hp, lp, budget, wgt},
   Do[dim = 3 g - 3 + hs + m + nd;
    If[dim < 2 nd, Continue[]];
    Do[budget = dim - Total[dils];
     Do[Do[
       wgt = hodgeNorm[g, Join[hp, lp, dils]]/(m! nd!);
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
numGraph[gr_, sub_] := <|gr, "decs" -> (gr["decs"] /. sub)|>;   (* decoration weights he1, he2, he3 -> numbers *)

(* ============================== 4. modified R-matrix rows and sampling ============================== *)
(* rows of R' = R_Z h^s (coefficient lists to z^KMAX); hc = coefficient list of h^s *)
rowsTimes[hc_List, KMAX_] := Table[Table[Together[Sum[PPsym[a][[j + 1]] hc[[k - j + 1]], {j, 0, k}]], {k, 0, KMAX}], {a, 0, 2}];
evalPointRows[t0_, x0_, KMAX_, rows_] := Block[{PPsym}, PPsym[0] = rows[[1]]; PPsym[1] = rows[[2]]; PPsym[2] = rows[[3]]; evalPoint[t0, x0, KMAX]];
(* one sample point for an orbit table OTx with vertex weights depending on he1, he2, he3 (hesub numeric) and R-matrix rows *)
computeSampleH[t0_, x0_, gg_, OTx_, rows_, hesub_] := Module[{tt = AbsoluteTime[], data, wmax = 3 gg - 3, vecs, imax = 0, bad = {}},
   data = If[rows === None, evalPoint[t0, x0, KMAXcur], evalPointRows[t0, x0, KMAXcur, rows]];
   vecs = Table[Module[{val = graphValueVE[numGraph[OTx["graphs"][id], hesub], data], rules, asc},
      rules = CoefficientRules[val, {C1i, C1p}]; asc = Association[rules];
      If[rules =!= {}, imax = Max[imax, Max[Abs[Im[rules[[All, 2]]]]]]];
      If[Select[rules, Abs[#[[2]]] > 10^-(PREC - 80) && (#[[1, 1]] != #[[1, 2]] || #[[1, 1]] > wmax) &] =!= {}, AppendTo[bad, id]];
      Table[ratl[Lookup[asc, Key[{k, k}], 0]/t0^(3 k)], {k, 0, wmax}]], {id, OTx["ids"]}];
   {{t0, x0}, vecs, {"imag" -> N[imax, 5], "bad" -> bad, "seconds" -> Round[AbsoluteTime[] - tt, 0.1]}}];
(* XX-coefficient vector of the reference total at a point *)
refVector[Fref_, t0_, x0_, gg_] := Table[Together[Coefficient[Fref, XX, k] /. {y -> t0^3, e1s -> x0, e2s -> 1}], {k, 0, 3 gg - 3}];

(* ============================== 5. reconstruction with (eps1 eps2)^a (eps1+eps2)^b denominators ============================== *)
AMAX = 8; BMAX = 10;
xsamplesH[nx_] := Table[Prime[r + 1], {r, 1, nx}];
xfitH[vals_, xs_, gg_] := Module[{res = $Failed, NX = Length[xs], deg, pts, pol, ok, a},
   If[Union[vals] === {0}, Return[{0, 0, {0}}]];
   Do[Do[a = tot - b; If[a < 0 || a > AMAX, Continue[]];
     deg = 2 gg - 2 + 2 a + b; If[deg + 2 > NX, Continue[]];
     pts = Table[{xs[[s]], xs[[s]]^a (1 + xs[[s]])^b vals[[s]]}, {s, deg + 1}];
     pol = InterpolatingPolynomial[pts, xv];
     ok = And @@ Table[(pol /. xv -> xs[[s]]) == xs[[s]]^a (1 + xs[[s]])^b vals[[s]], {s, deg + 2, NX}];
     If[ok, res = {a, b, CoefficientList[Expand[pol], xv, deg + 1]}; Return[]], {b, 0, Min[tot, BMAX]}];
    If[res =!= $Failed, Return[]], {tot, 0, AMAX + BMAX}];
   res];
reconstructOrbitH[samp_, ts_, xs_, gg_] := Module[{term = 0, wmax = 3 gg - 3, xf, amax, bmax, deg, seqs, yf, ok = True, info = {}},
   Do[xf = Table[xfitH[Table[samp[t0][x0][[k + 1]], {x0, xs}], xs, gg], {t0, ts}];
    If[MemberQ[xf, $Failed], ok = False; AppendTo[info, {k, "x"}]; Continue[]];
    amax = Max[xf[[All, 1]]]; bmax = Max[xf[[All, 2]]]; deg = 2 gg - 2 + 2 amax + bmax;
    seqs = Table[Table[Module[{a = xf[[s, 1]], b = xf[[s, 2]], cl = xf[[s, 3]]},
         Coefficient[Expand[xv^(amax - a) (1 + xv)^(bmax - b) Sum[cl[[i + 1]] xv^i, {i, 0, Length[cl] - 1}]], xv, jj]], {s, Length[ts]}], {jj, 0, deg}];
    If[Union[Flatten[seqs]] === {0}, Continue[]];
    yf = yfitC[seqs, ts, gg];
    If[yf === $Failed, ok = False; AppendTo[info, {k, "y"}]; Continue[]];
    term += XX^k Sum[yf[[3, jj + 1]] e1s^jj e2s^(deg - jj), {jj, 0, deg}]/((e1s e2s)^amax (e1s + e2s)^bmax), {k, 0, wmax}];
   If[ok, term, {$Failed, info}]];
