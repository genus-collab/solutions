(* lemmaA_data.wls -- exact Hodge-side tables for the Lemma A workflow (no topological recursion anywhere).
   Hodge class Omega = prod_{a=1}^3 Lambda^vee_g(w_a)/w_a^g with w = (c, -x, -1), c = 1 + x (hvc_lib.m hodgeNormU4, Mumford's formula),
   leaves lam_j(psi) = 1/(c - j psi). Output files in lemmaA/data/ (one record per line, fields separated by " | ", values in InputForm):
   psi_g<g>.txt      g | n | {a_1..a_n} (sorted) | <Omega psi_1^a_1 ... psi_n^a_n>_g        for every sorted tuple with sum <= 3g-3+n
   mv_g<g>_l<l>.txt  g | l | P(m1..ml) = <Omega prod_i lam_{m_i}(psi_i)>_g as a polynomial in symbols m1..ml (no dilaton leaves)
   q_g<g>_r<r>.txt   g | r | Q(s; j1..jr) = Sum_k s^k/k! <Omega prod lam_{j_i} Tt^k>_g,  Tt(psi) = psi^2/(c - psi),
                     a polynomial in s and j1..jr; s = -1 gives Xi_g(j1..jr) of hvc_multiprobe.wls (dilaton leaves T_c = -Tt).
   Usage: wolframscript -file lemmaA_data.wls *)
SetDirectory[DirectoryName[$InputFileName]];
Get[FileNameJoin[{ParentDirectory[], "hvc_lib.m"}]];
t0 = AbsoluteTime[];
c = 1 + x; usub = usubOf[{c, -x, -1}];
fmt[e_] := ToString[InputForm[e]];
sortedTuples[n_, smax_] := Select[Tuples[Range[0, smax], n], OrderedQ[#] && Total[#] <= smax &];
om[g_, a_List] := om[g, a] = Together[hodgeNormU4[g, Sort[a]] /. usub];
(* 3. s-interpolated dilaton leaves: Q = Sum_k s^k/k! <prod lam_j Tt^k> *)
TCs[d_] := s/c^(d - 1);
msets[k_, lo_, hi_, smax_] := Module[{rec}, rec[kk_, l_, sm_] := If[kk == 0, {{}}, Join @@ Table[Prepend[#, aa] & /@ rec[kk - 1, aa, sm - aa], {aa, l, Min[hi, sm]}]]; rec[k, lo, smax]];
ordFactor[lst_] := Length[lst]!/Times @@ (Last /@ Tally[lst]!);
Lvals[h_, a_List] := Lvals[h, a] = Module[{tot = 0, r = Length[a], sa = Total[a], dim},
   Do[dim = 3 h - 3 + r + m;
    Do[tot += (ordFactor[DD]/m!) (Times @@ (TCs /@ DD)) om[h, Join[a, DD]], {DD, msets[m, 2, dim, dim - sa]}],
    {m, 0, 3 h - 3 + r - sa}];
   Together[tot]];
jsym = {j1, j2, j3};
Do[Module[{h = hr[[1]], r = hr[[2]], pol = 0, amax},
   amax = 3 h - 3 + r;
   Do[If[Total[a] <= amax, pol += (Times @@ Table[jsym[[i]]^a[[i]]/c^(a[[i]] + 1), {i, r}]) Lvals[h, Sort[a]]], {a, Tuples[Range[0, amax], r]}];
   pol = Together[pol];
   Export["data/q_g" <> ToString[h] <> "_r" <> ToString[r] <> ".txt", StringRiffle[{ToString[h], ToString[r], fmt[pol]}, " | "], "Text"];
   Print["Q table genus ", h, ", r = ", r, " [", Round[AbsoluteTime[] - t0], " s]"]],
  {hr, {{2, 1}, {2, 2}, {3, 1}}}];
Print["total ", Round[AbsoluteTime[] - t0], " s"];
