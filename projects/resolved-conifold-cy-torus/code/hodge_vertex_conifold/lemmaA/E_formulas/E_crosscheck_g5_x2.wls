(* E_crosscheck_g5.wls -- independent Mathematica cross-check (requested check of Route E) of the Perl engine
   lemmaA/E_formulas/Hodge.pm + HodgePt.pm at GENUS 5 (where hvc_lib.m's hodgeNormU4 is not available: it needs the ch_5 branch).
   Self-contained: Witten-Kontsevich numbers (verbatim from graphs.m), Mumford's GRR recursion for Int ch_1^k1 ch_3^k3 ch_5^k5 psi^a
   (the same recursion as hvc_lib.m's HI, extended to ch_5), the CY Hodge class prod_a Lambda^vee_5(w_a)/w_a^5, w = (c,-x,-1), c = 1+x,
   leaves 1/(c - j psi), dilaton leaves T_c = -psi^2/(c - psi) (as in hvc_multiprobe_g4b.wls).  NO topological recursion.
   Compares Xi_5(R) (c-normalisation) with the engine predictions in expected_g5.m (written by expected_g5.pl).
   Lemma A points: {2},{3},{8},{2,2},{2,7},{3,6},{2,2,2},{2,2,6} (predicted 0); controls {9},{2,8},{2,2,7} (predicted non-zero).
   Cheaper variant: set xval = 2 (or any rational) below to work at a numeric point.
   Usage: wolframscript -file E_crosscheck_g5.wls  ->  E_crosscheck_g5_out.txt *)
SetDirectory[DirectoryName[$InputFileName]];
xval = 2;   (* main session: numeric run at x = 2 *)
out = {}; emit[s_] := (AppendTo[out, s]; Print[s]; Export["E_crosscheck_g5_out.txt", StringRiffle[out, "\n"], "Text"]);
(* ---------- Witten-Kontsevich numbers (graphs.m) ---------- *)
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
(* ---------- Mumford recursion with ch_1, ch_3, ch_5 ---------- *)
chl = {1, 3, 5};
subsetSum[a_List, f_] := Module[{tal = Tally[a], tot = 0},
   Do[tot += (Times @@ MapThread[Binomial, {tal[[All, 2]], sel}]) f[Flatten[MapThread[ConstantArray, {tal[[All, 1]], sel}]], Flatten[MapThread[ConstantArray, {tal[[All, 1]], tal[[All, 2]] - sel}]]],
    {sel, Tuples[Range[0, #] & /@ tal[[All, 2]]]}];
   tot];
HI5[g_, ks_List, a_List] /; a =!= Sort[a] := HI5[g, ks, Sort[a]];
HI5[g_, ks_List, a_List] := HI5[g, ks, a] = Module[{n = Length[a], dch = ks . chl, li, l, kp, tot},
   Which[
    g < 0, 0,
    n > 0 && Min[a] < 0, 0,
    2 g - 2 + n <= 0, 0,
    dch + Total[a] != 3 g - 3 + n, 0,
    dch == 0, wk[g, a],
    g == 0, 0,
    g == 1 && dch > 1, 0,
    True,
    li = Max[Flatten[Position[ks, _?Positive]]]; l = chl[[li]]; kp = ReplacePart[ks, li -> ks[[li]] - 1];
    tot = HI5[g, kp, Append[a, l + 1]] - Sum[HI5[g, kp, ReplacePart[a, i -> a[[i]] + l]], {i, n}]
      + 1/2 Sum[(-1)^p HI5[g - 1, kp, Join[a, {p, l - 1 - p}]], {p, 0, l - 1}];
    tot += 1/2 subsetSum[a, Function[{aS, aSc},
        Sum[(Times @@ MapThread[Binomial, {kp, js}]) (-1)^p HI5[h, js, Append[aS, p]] HI5[g - h, kp - js, Append[aSc, l - 1 - p]],
         {h, 0, g}, {js, Tuples[Range[0, #] & /@ kp]}, {p, 0, l - 1}]]];
    BernoulliB[l + 1]/(l + 1)! tot]];
(* ---------- the CY Hodge class at genus g in the ch basis ---------- *)
c = 1 + x; If[xval =!= None, c = 1 + xval];
xx = If[xval === None, x, xval];
hodgeRules[g_] := hodgeRules[g] = Module[{lam, om},
   lam = CoefficientList[Normal[Series[Exp[eps ch[1] + 2 eps^3 ch[3] + 24 eps^5 ch[5]], {eps, 0, g}]], eps];
   om = Expand[Product[Sum[(-1)^k lam[[k + 1]] w^(-k), {k, 0, g}], {w, {1 + xx, -xx, -1}}]];
   CoefficientRules[om, {ch[1], ch[3], ch[5]}]];
hodgeNorm[g_, a_List] := hodgeNorm[g, Sort[a]] = Together[Total[(#[[2]] HI5[g, #[[1]], Sort[a]]) & /@ Select[hodgeRules[g], #[[1]] . chl + Total[a] == 3 g - 3 + Length[a] &]]];
(* ---------- dilaton leaves and lattice-point evaluation (as hvc_multiprobe_g4b.wls) ---------- *)
TC[d_] := -1/c^(d - 1);
msets[k_, lo_, hi_, smax_] := Module[{rec}, rec[kk_, l_, s_] := If[kk == 0, {{}}, Join @@ Table[Prepend[#, aa] & /@ rec[kk - 1, aa, s - aa], {aa, l, Min[hi, s]}]]; rec[k, lo, smax]];
ordFactor[lst_] := Length[lst]!/Times @@ (Last /@ Tally[lst]!);
Lval[h_, a_List] := Lval[h, a] = Module[{tot = 0, r = Length[a], sa = Total[a], dim},
   Do[dim = 3 h - 3 + r + m;
    Do[tot += (ordFactor[D]/m!) (Times @@ (TC /@ D)) hodgeNorm[h, Join[a, D]], {D, msets[m, 2, dim, dim - sa]}],
    {m, 0, 3 h - 3 + r - sa}];
   Together[tot]];
XiAt[h_, jl_List] := Module[{r = Length[jl], amax, tot = 0},
   amax = 3 h - 3 + r;
   Do[tot += (Times @@ (1/c^(a + 1))) (Total[Table[Times @@ (jl^p), {p, Permutations[a]}]]) Lval[h, a], {a, msets[r, 0, amax, amax]}];
   Together[tot]];
(* ---------- sanity checks of the genus-5 Hodge part: Faber, Getzler-Pandharipande, lambda_g ---------- *)
lamInt[g_, mono_List, a_List] := Module[{lam, p},
   lam = CoefficientList[Normal[Series[Exp[eps ch[1] + 2 eps^3 ch[3] + 24 eps^5 ch[5]], {eps, 0, g}]], eps];
   p = Expand[Product[lam[[k + 1]], {k, mono}]];
   Total[(#[[2]] HI5[g, #[[1]], a]) & /@ CoefficientRules[p, {ch[1], ch[3], ch[5]}]]];
emit["genus 5 lambda_5 lambda_4 lambda_3 = " <> ToString[lamInt[5, {5, 4, 3}, {}]] <> " (Faber: " <> ToString[Abs[BernoulliB[8]] Abs[BernoulliB[10]]/(2 8! 8 10)] <> ")"];
emit["genus 5 psi^4 lambda_5 lambda_4 = " <> ToString[lamInt[5, {5, 4}, {4}]] <> " (Getzler-Pandharipande: " <> ToString[Abs[BernoulliB[10]]/(2^9 df[9] 10)] <> ")"];
(* ---------- the cross-check ---------- *)
Get["expected_g5.m"];
t0 = AbsoluteTime[];
Do[Module[{p = pr[[1]], e = pr[[2]], v, d},
   v = XiAt[5, p]; e = If[xval === None, e, e /. x -> xval]; d = Together[v - e];
   emit["genus 5, R = " <> ToString[p] <> ": engine prediction " <> If[d === 0, "CONFIRMED", "DISAGREES (difference " <> ToString[InputForm[Factor[d]]] <> ")"] <>
     "; value " <> If[v === 0, "0", "non-zero"] <> "   [" <> ToString[Round[AbsoluteTime[] - t0]] <> " s]"]], {pr, expected5}];
emit["done"];
