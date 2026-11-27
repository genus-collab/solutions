(* hvc_multiprobe.wls -- multi-probe mixed-scale vertex functionals.
   Leaves lam_j(psi) = 1/(c - j psi) (j = 0: constant 1/c), dilaton leaves T_c = -psi^2/(c - psi), Hodge weights (c, -x, -1), c = 1 + x.
   Xi[h, r] = Sum_m 1/m! < lam_{j1} ... lam_{jr} T_c^m >_h as a polynomial in symbolic jj[1..r].
   Conjecture M: Xi = 0 whenever the j_i are integers >= 0 with 1 <= J = Sum j_i <= 2h - 3 + r.
   lambda_h limit (x -> 0): Xi ~ prod_{l=1}^{2h-3+r} (J - l) (Norlund B_N^(N+1)(J) = (J-1)...(J-N)). *)
SetDirectory[DirectoryName[$InputFileName]];
Get["hvc_lib.m"];
c = 1 + x; usub = usubOf[{c, -x, -1}];
TC[d_] := -1/c^(d - 1);
msets[k_, lo_, hi_, smax_] := Module[{rec}, rec[kk_, l_, s_] := If[kk == 0, {{}}, Join @@ Table[Prepend[#, a] & /@ rec[kk - 1, a, s - a], {a, l, Min[hi, s]}]]; rec[k, lo, smax]];
ordFactor[lst_] := Length[lst]!/Times @@ (Last /@ Tally[lst]!);
Lval[h_, a_List] := Lval[h, a] = Module[{tot = 0, r = Length[a], sa = Total[a], dim},
   Do[dim = 3 h - 3 + r + m;
    Do[tot += (ordFactor[D]/m!) (Times @@ (TC /@ D)) (hodgeNormU4[h, Join[a, D]] /. usub), {D, msets[m, 2, dim, dim - sa]}],
    {m, 0, 3 h - 3 + r - sa}];
   Together[tot]];
Xi[h_, r_] := Xi[h, r] = Module[{tot = 0, amax = 3 h - 3 + r},
   Do[If[Total[a] <= amax, tot += (Times @@ Table[jj[i]^a[[i]]/c^(a[[i]] + 1), {i, r}]) Lval[h, Sort[a]]], {a, Tuples[Range[0, amax], r]}];
   Together[tot]];
latticeReport[h_, r_, jmax_] := Module[{pts, bad = {}, xi = Xi[h, r]},
   pts = Select[Tuples[Range[0, jmax], r], OrderedQ[#] && Total[#] <= 2 h - 1 + r &];
   Do[Module[{v = Together[xi /. Thread[Array[jj, r] -> p]], pred = 1 <= Total[p] <= 2 h - 3 + r},
      If[(v === 0) =!= pred, AppendTo[bad, {p, v === 0, pred}]]], {p, pts}];
   Print["genus ", h, ", r = ", r, ": lattice points tested ", Length[pts], " (all ordered tuples with J <= 2h-1+r); disagreements with Conjecture M: ", bad]];
lamLimit[h_, r_] := Module[{xi = Xi[h, r], lead, J = Total[Array[jj, r]], norl},
   lead = Factor[Limit[x^h xi, x -> 0]];
   norl = Product[J - l, {l, 1, 2 h - 3 + r}];
   Print["genus ", h, ", r = ", r, ": x->0 leading coefficient / prod_{l=1}^{2h-3+r}(J-l) = ", Factor[lead/norl]]];
t0 = AbsoluteTime[];
latticeReport[4, 1, 9]; lamLimit[4, 1];
Print["genus 4 r=1 symbolic: ", Factor[Xi[4, 1]]];
Print["time so far ", Round[AbsoluteTime[] - t0], " s"];
latticeReport[4, 2, 9]; lamLimit[4, 2];
Print["total time ", Round[AbsoluteTime[] - t0], " s"];
