(* Hodge integrals of the genus-general evaluator (hv_lib4.m): known values up to genus 4, vanishing above the top degree, Mumford
   relations, closed forms versus the recursion, string and dilaton equations *)
{hodgeOK, hodgeReport} = hodgeChecks4[]; Column[Append[hodgeReport, "all checks passed: " <> ToString[hodgeOK]]]
(* the scalar R-action of h_w^-1 on the point theory (PPZ conventions) reproduces the normalised Hodge class on M_2 for arbitrary weights ... *)
ppz2 = ppzScalar[2, {w1, w2, w3}]; {Factor[ppz2], Together[ppz2 - (hodgeNormU4[2, {}] /. usubOf[{w1, w2, w3}])] === 0}
(* ... and on M_3 (42 stable graphs) at exact rational weights *)
SeedRandom[7]; Table[Module[{wn = Table[RandomInteger[{-9, 9}]/RandomInteger[{1, 5}], {3}]}, If[MemberQ[wn, 0], Nothing, {wn, ppzScalar[3, wn] == (hodgeNormU4[3, {}] /. usubOf[wn])}]], {2}]
(* the rank-4 Hodge class of the actual normal bundle O(-1)+O(-1)+O+O of P^1 in the 5-fold, weights (eps1, eps2, -eps1, -eps2), is
   trivial by Mumford's c(E) c(E^vee) = 1: it integrates to the Witten-Kontsevich numbers for every psi assignment (as a polynomial in
   lambda_1, ch_3 it is not literally 1, since lambda_4 = 0 in genus 3 and lambda_g^2 = 0 are not in the substitution table) *)
usub4 = {u[1] -> 1/e1, u[2] -> 1/e2, u[3] -> -1/e1, u[4] -> -1/e2};
hodgeNormN[g_, nw_, a_List] := Expand[Total[Table[c[[2]] HI[g, c[[1, 1]], c[[1, 2]], a], {c, hodgeClassN[g, nw]}]]];
Table[Module[{as = Select[Join @@ Table[Tuples[Range[0, 3 g - 3 + n], n], {n, 1, If[g == 3, 3, 4]}], Total[#] <= 3 g - 3 + Length[#] &]}, {g, Length[as], "disagreements with Int psi^a:", Count[as, a_ /; Together[(hodgeNormN[g, 4, a] /. usub4) - wk[g, a]] =!= 0]}], {g, 1, 3}]
(* the triple Hodge class of local P2 (mode B) in genus 2 and 3, as coefficients of lambda_1^k1 ch_3^k3, and the scalar h_wB(z) to z^6 *)
Table[{g, Table[{c[[1]], Factor[c[[2]] /. usubOf[wB]]}, {c, hodgeClass4[g]}]}, {g, 2, 3}] // TableForm
Factor /@ hCoeffs[wB, 1, 6]
(* the conifold Frobenius data at alpha0 = 0, alpha_i = eps_i: the quadratic equation of the canonical coordinates U^2 - S U + P = 0,
   whose roots are U = 1 and U = -1/x at Q = 1 (one canonical coordinate runs into the logarithmic singularity of the superpotential),
   and the discriminant Q (1-x)^2 + 4x (times x^2 Q), which vanishes at Q_c = -4x/(1-x)^2, not at Q = 1 *)
{Together[(a0 + Q (a1 - a2))/(Q (a0 + a1)) /. {a0 -> 0, a1 -> x, a2 -> 1}], Together[-a2/(Q (a0 + a1)) /. {a0 -> 0, a1 -> x, a2 -> 1}], Factor[Discriminant[Q x U^2 - Q (x - 1) U - 1, U]], Solve[(U^2 - (x - 1)/x U - 1/x) == 0, U]}
