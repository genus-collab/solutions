(* hvc_rank4_check.wls -- the rank-4 Hodge class of the actual normal bundle O(-1)+O(-1)+O+O of P^1 in the 5-fold, weights
   (eps1, eps2, -eps1, -eps2), integrates to the Witten-Kontsevich numbers: Int prod_a Lambda_g^vee(w_a)/w_a^g psi^a = Int psi^a
   for all psi-degrees (the class is 1 by Mumford's c(E) c(E^vee) = 1; as a polynomial in lambda_1, ch_3 it is not literally 1
   because lambda_4 = 0 (genus 3) and lambda_g^2 = 0 are not in the substitution table, but these monomials integrate to zero). *)
SetDirectory[DirectoryName[$InputFileName]];
Get["hvc_lib.m"];
usub4 = {u[1] -> 1/e1s, u[2] -> 1/e2s, u[3] -> -1/e1s, u[4] -> -1/e2s};
hodgeNormN[g_, nw_, a_List] := Expand[Total[Table[c[[2]] HI[g, c[[1, 1]], c[[1, 2]], a], {c, hodgeClassN[g, nw]}]]];
tuples[g_, nmax_] := Join @@ Table[Tuples[Range[0, 3 g - 3 + n], n], {n, 1, nmax}];
res = Table[{g, Module[{as = Select[tuples[g, If[g == 3, 3, 4]], Total[#] <= 3 g - 3 + Length[#] &], bad},
     bad = Select[as, Together[(hodgeNormN[g, 4, #] /. usub4) - wk[g, #]] =!= 0 &];
     {Length[as], "psi-tuples checked; disagreements with Int psi^a:", Length[bad], If[bad =!= {}, bad[[;; Min[5, Length[bad]]]], ""]}]}, {g, 1, 3}];
Print[res];
Print["rank-4 class integrates to the WK numbers for all tuples: ", And @@ Table[r[[2, 3]] == 0, {r, res}]];
