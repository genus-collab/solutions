(* Hodge integrals: known values, Mumford recursion against the closed forms, string and dilaton equations *)
{hodgeOK, hodgeReport} = hodgeChecks3[]; Column[Append[hodgeReport, "all checks passed: " <> ToString[hodgeOK]]]
(* the scalar R-action of h on the point theory (PPZ conventions) reproduces the normalised Hodge class on M_2 for arbitrary weights ... *)
ppz2 = ppzScalar[2, {w1, w2, w3}]; {Factor[ppz2], Together[ppz2 - (hodgeNormU[2, {}] /. usubOf[{w1, w2, w3}])] === 0}
(* ... and on M_3 (42 stable graphs) at exact rational weights *)
SeedRandom[7]; Table[Module[{wn = Table[RandomInteger[{-9, 9}]/RandomInteger[{1, 5}], {3}]}, If[MemberQ[wn, 0], Nothing, {wn, ppzScalar[3, wn] == (hodgeNormU[3, {}] /. usubOf[wn])}]], {2}]
(* the normalised Hodge class of genus 3 for the weights wB, as coefficients of lambda_1^i lambda_2^j lambda_3^k *)
Table[{c[[1]], Factor[c[[2]] /. usubB]}, {c, hodgeClass[3]}] // TableForm
(* the scalar h_wB(z) to z^6 *)
Factor /@ hCoeffs[wB, 1, 6]
(* the genus-general evaluator (hv_lib4.m): known values up to genus 4 (lambda_1^9 on M_4 = 1/113400 by nine-fold Mumford, the
   Getzler-Pandharipande value lambda_4 lambda_3 psi^2 psi^2 = 1/276480, lambda_4 lambda_3 lambda_2 = 1/87091200, ...), vanishing above
   the top degree, closed forms versus the recursion, string/dilaton; and agreement with hv_lib3's correlators on the genus-3 decorations *)
{hodge4OK, hodge4Report} = hodgeChecks4[]; Column[Append[hodge4Report, "all checks passed: " <> ToString[hodge4OK]]]
If[ValueQ[preparedH3] || FileExistsQ["hv_orbits_g3.mx"], If[! ValueQ[preparedH3], Get["hv_orbits_g3.mx"]];
  Module[{d3 = Union[Flatten[Table[Table[{gr["g"][[v]], Join[dec[[v, 1]], dec[[v, 2]], dec[[v, 3]]]}, {v, gr["V"]}], {gr, preparedH3}, {dec, gr["decs"]}], 2]]},
   {Length[d3], And @@ Table[Together[hodgeNormU4[d[[1]], d[[2]]] - hodgeNormU[d[[1]], d[[2]]]] === 0, {d, d3}]}]]
(* the normalised genus-4 Hodge class for wB as coefficients of lambda_1^k1 ch_3^k3 *)
Table[{c[[1]], Factor[c[[2]] /. usubOf[wB]]}, {c, hodgeClass4[4]}] // TableForm
