(* hvc_tr_tau.wls -- exact Eynard-Orantin recursion on the base curve C_tau of the dictionary P2:
     X = log x = tau log v - (tau+1) log(1+v),  y = -v,  B = dv dv'/(v-v')^2,  single branch point v = tau,
   local coordinate z at v = tau: X = X(tau) - z^2/(2 tau (1+tau)); basis e_k(p) = Res_{q->tau} B(p,q) z(q)^(-2k-1)
   (= (2k+1) dz/z^(2k+2) + regular).  Recursion in ABCD form with precomputed residue tables.
   Output: Xi^TR_g(mu) := sigma^(g-1) [x^mu] F_{g,n} / prod binom(mu_i(tau+1)-1, mu_i-1)  (c = 1 normalisation, P2),
   sigma = tau(1+tau); with "hodge", for g <= 4 compared with the exact Hodge functional (hvc_mp_header.m) at x = -tau/(1+tau).
   Usage: wolframscript -file hvc_tr_tau.wls <tau> <gmax> <nmax> [hodge]      Output: hvc_tr_tau<tau>_g<gmax>n<nmax>.txt *)
SetDirectory[DirectoryName[$InputFileName]];
{tauv, gmax, nmax} = ToExpression /@ $ScriptCommandLine[[2 ;; 4]];
doHodge = Length[$ScriptCommandLine] >= 5 && $ScriptCommandLine[[5]] == "hodge";
tag = "hvc_tr_tau" <> StringReplace[ToString[tauv, InputForm], {"/" -> "_", "-" -> "m"}] <> "_g" <> ToString[gmax] <> "n" <> ToString[nmax];
out = {}; emit[s_] := (AppendTo[out, s]; Print[s]; Export[tag <> ".txt", StringRiffle[out, "\n"], "Text"]);
t0 = AbsoluteTime[];
kmax = 3 gmax - 3 + nmax; ordB = 2 kmax + 4; ordV = 2 kmax + ordB + 8;
sig = tauv (1 + tauv);
trunc[p_, d_] := Module[{e = Expand[p]}, If[e === 0, 0, Total[Select[If[Head[e] === Plus, List @@ e, {e}], Exponent[#, z] <= d &]]]];
(* local coordinate *)
Xloc = tauv Log[1 + u/tauv] - (tauv + 1) Log[1 + u/(1 + tauv)];
qq = Normal[Series[-2 sig Xloc/u^2, {u, 0, ordV + 2}]];
zOfU = Series[u Sqrt[qq], {u, 0, ordV + 2}];
vz = Expand[tauv + Normal[InverseSeries[zOfU, z]]];
cv = CoefficientList[vz, z]; uz = Expand[vz - tauv]; dvz = D[vz, z];
Yhat = Expand[(vz - (vz /. z -> -z))/(z sig)];                        (* Y = z^2 Yhat, Y dz = omega01(q) - omega01(sigma q) *)
invY = Expand[Normal[Series[1/Yhat, {z, 0, ordV - 2}]]];
(* basis differentials in the local coordinate: e_k(p) = Sum_{j=0}^{2k} (j+1) A_{k,j} dp/(p-tau)^(j+2), A_{k,j} = [z^(2k)](v'(z) u(z)^j),
   hence E_k(z) := e_k/dz at p = v(z) = v'(z) Sum_j (j+1) A_{k,j} u(z)^(-j-2) *)
Akj[k_, j_] := Coefficient[trunc[dvz uz^j, 2 k], z, 2 k];
zOverU = Normal[Series[z/uz, {z, 0, ordB + 2 kmax + 6}]];
zuPow[p_] := zuPow[p] = trunc[zOverU^p, ordB + 2 kmax + 6];
EE[k_] := EE[k] = trunc[Expand[dvz Sum[(j + 1) Akj[k, j] z^(-j - 2) zuPow[j + 2], {j, 0, 2 k}]], ordB];
res[expr_] := Coefficient[Expand[expr], z, -1];
(* residue tables *)
rhoAB[k0_, a_, b_] := rhoAB[k0, a, b] = -res[z^(2 k0 - 1) EE[a] (EE[b] /. z -> -z) invY]/(2 k0 + 1);
rhoC[k0_, l_, b_] := rhoC[k0, l, b] = -res[z^(2 k0 - 1) z^(2 l) (EE[b] + (EE[b] /. z -> -z)) invY]/(2 k0 + 1);
emit["tau = " <> ToString[tauv, InputForm] <> ", sigma = " <> ToString[sig, InputForm] <> "; v(z) = " <> ToString[InputForm[Take[cv, 5]]] <> "...; E_1(z) = " <> ToString[InputForm[trunc[EE[1], 1]]] <> "; setup " <> ToString[Round[AbsoluteTime[] - t0, 0.1]] <> " s"];
(* initial data *)
W = <||>;
W[{0, 3}] = Association[{0, 0, 0} -> -2 Coefficient[invY, z, 0]];
B11 = trunc[Normal[Series[-dvz (dvz /. z -> -z)/(vz - (vz /. z -> -z))^2, {z, 0, 6}]], 6];
W[{1, 1}] = Association[Select[Table[{k0} -> res[z^(2 k0 - 1) B11 invY]/(2 k0 + 1), {k0, 0, 1}], Last[#] =!= 0 &]];
emit["omega_{0,3} = " <> ToString[InputForm[W[{0, 3}]]] <> ", omega_{1,1} = " <> ToString[InputForm[W[{1, 1}]]]];
stableQ[g_, n_] := 2 g - 2 + n > 0;
SetAttributes[addTo, HoldFirst];
addTo[assoc_, key_, val_] := If[val =!= 0, assoc[key] = Lookup[assoc, Key[key], 0] + val];
compute[g_, np1_] := Module[{n = np1 - 1, acc = <||>, kk0max = 3 g - 3 + np1},
   (* A: omega_{g-1,n+2}(q, sigma q, J) *)
   If[g >= 1 && stableQ[g - 1, n + 2],
    KeyValueMap[Function[{key, cf}, Do[addTo[acc, Prepend[key[[3 ;;]], k0], cf rhoAB[k0, key[[1]], key[[2]]]], {k0, 0, kk0max}]], W[{g - 1, n + 2}]]];
   (* B: stable products *)
   Do[Module[{I = Select[Range[n], BitGet[mask, # - 1] == 1 &], Ic},
      Ic = Complement[Range[n], I];
      If[stableQ[h, 1 + Length[I]] && stableQ[g - h, 1 + Length[Ic]],
       KeyValueMap[Function[{k1, c1}, KeyValueMap[Function[{k2, c2}, Module[{kJ = ConstantArray[0, n]},
             kJ[[I]] = k1[[2 ;;]]; kJ[[Ic]] = k2[[2 ;;]];
             Do[addTo[acc, Prepend[kJ, k0], c1 c2 rhoAB[k0, k1[[1]], k2[[1]]]], {k0, 0, kk0max}]]], W[{g - h, 1 + Length[Ic]}]]], W[{h, 1 + Length[I]}]]]],
    {h, 0, g}, {mask, 0, 2^n - 1}];
   (* C: omega_{0,2}(q, p_j) omega_{g,n}(sigma q, J\j) + (q <-> sigma q) *)
   If[n >= 1 && stableQ[g, n],
    KeyValueMap[Function[{key, cf}, Do[Do[addTo[acc, Prepend[Insert[key[[2 ;;]], l, j], k0], cf rhoC[k0, l, key[[1]]]], {l, 0, key[[1]] + 1}, {k0, 0, key[[1]] + 1 - l}], {j, 1, n}]], W[{g, n}]]];
   W[{g, np1}] = Select[acc, # =!= 0 &]];
(* all (g,n) needed: n <= nmax + 2 (gmax - g) *)
order = SortBy[Select[Flatten[Table[{g, n}, {g, 0, gmax}, {n, 1, nmax + 2 (gmax - g)}], 1], stableQ @@ # && ! MemberQ[{{0, 1}, {0, 2}, {0, 3}, {1, 1}}, #] &], {2 #[[1]] - 2 + #[[2]], #[[1]]} &];
Do[compute @@ gn, {gn, order}];
emit["recursion done: " <> ToString[Length[order]] <> " correlators, " <> ToString[Round[AbsoluteTime[] - t0, 0.1]] <> " s"];
symBad[as_] := Count[KeyValueMap[Function[{k, v}, Lookup[as, Key[Sort[k]], 0] - v], as], x_ /; x =!= 0];
(* x-expansion at v = oo:  [x^mu] xi_k = -((2k+1)/mu) Sum_{m=0}^{mu-1} [z^(2k+1)](v^(m+1)) binom((tau+1) mu, mu-m-1) *)
vpow[p_] := vpow[p] = trunc[vz^p, 2 kmax + 3];
xiCoef[k_, mu_] := xiCoef[k, mu] = -((2 k + 1)/mu) Sum[Coefficient[vpow[m + 1], z, 2 k + 1] Binomial[(tauv + 1) mu, mu - m - 1], {m, 0, mu - 1}];
xCoef[gn_, mu_List] := Total[KeyValueMap[Function[{k, v}, v Times @@ MapThread[xiCoef, {k, mu}]], W[gn]]];
XiTR[g_, mu_List] := sig^(g - 1) xCoef[{g, Length[mu]}, mu]/Times @@ (Binomial[# (tauv + 1) - 1, # - 1] & /@ mu);
If[doHodge, Get["hvc_mp_header.m"];
  XiAt[h_, jl_List] := Module[{r = Length[jl], amax, tot = 0}, amax = 3 h - 3 + r;
     Do[tot += (Times @@ (1/c^(a + 1))) (Total[Table[Times @@ (jl^p), {p, Permutations[a]}]]) Lval[h, a], {a, msets[r, 0, amax, amax]}]; Together[tot]];
  xH = -tauv/(1 + tauv);
  XiH[g_, mu_List] := (c^(3 g - 3 + 2 Length[mu]) XiAt[g, mu]) /. x -> xH];
Do[If[gn[[1]] >= 2 && gn[[2]] <= nmax,
   Module[{g = gn[[1]], n = gn[[2]], lim, tups, bad = {}, cmpBad = {}, firstNZ = {}, ncmp = 0},
    lim = 2 g - 2 + n + If[n == 1, 3, 1];
    tups = Select[Tuples[Range[1, lim], n], OrderedQ[#] && Total[#] <= lim &];
    Do[Module[{val = XiTR[g, mu], pred = Total[mu] <= 2 g - 3 + n},
       If[(val === 0) =!= pred, AppendTo[bad, mu -> val]];
       If[! pred && Total[mu] == 2 g - 2 + n && val =!= 0 && Length[firstNZ] < 3, AppendTo[firstNZ, mu -> val]];
       If[doHodge && g <= 4 && Total[mu] <= 2 g - 2 + n + 1, ncmp++; Module[{hv = XiH[g, mu]}, If[Together[hv - val] =!= 0, AppendTo[cmpBad, {mu, val, hv}]]]]], {mu, tups}];
    emit["(" <> ToString[g] <> "," <> ToString[n] <> "): " <> ToString[Length[W[gn]]] <> " terms, symmetry violations " <> ToString[symBad[W[gn]]] <>
      "; lattice points |mu| <= " <> ToString[lim] <> ": " <> ToString[Length[tups]] <> ", disagreements with Lemma D / Conjecture M: " <> ToString[InputForm[bad]] <>
      "; first layer " <> ToString[InputForm[firstNZ]] <> If[doHodge && g <= 4, "; TR vs Hodge (P2) at " <> ToString[ncmp] <> " points, mismatches: " <> ToString[InputForm[cmpBad]], ""] <>
      "  [" <> ToString[Round[AbsoluteTime[] - t0, 0.1]] <> " s]"]]], {gn, order}];
