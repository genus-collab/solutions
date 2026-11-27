(* hvc_leaf_closed_forms.wls -- closed forms of all leaf coefficients (leaf_proof/, 2026-10-31).
   Claim: rho q_t(psi) = c + Sum_{j>=1} A_j(t) lam_j(psi), lam_j = 1/(c - j psi), exactly, with
     A_j(t) = (-1)^j c^2 (tau/t) Exp[j dphi(tau)/c] tau^(j-1) Pochhammer[1 - j/c, j-1]^2/(j ((j-1)!)^2) Hypergeometric2F1[j x/c, j x/c, j, tau],
   tau = 1 - e^-t, dphi = phi_tau(sst) - phi_0(1/c) (see leaf_proof/outline.md), and the E11 loop coefficients
     t^2 rho^2 E11 = Sum_{i,l} e_il(t) lam_i(psi) lam_l(psi'),  e_il(t) = (i l/c^2) Int_0^t Theta d/dt(t A_i) d/dt(t A_l),  Theta = -x^3 tau sst^2.
   Checks: (a) against the closed form G of the singular row (hvc_leaf_proof_check.wls, matched to the stored factors) for k <= KX;
   (b) Sum_j A_j = -c^2 (the row vanishes at psi = 0); (c) against the stored fits of leaf_structure.md (hvc_leaf_structure_data_kt5.m),
   dilaton k <= 5 and loops k <= 5 (restricted ansatz).
   Usage: wolframscript -file hvc_leaf_closed_forms.wls [KX]     Output: hvc_leaf_closed_forms_out.txt *)
SetDirectory[DirectoryName[$InputFileName]];
args = Rest[$ScriptCommandLine]; KX = If[Length[args] >= 1, ToExpression[args[[1]]], 10];
t0 = AbsoluteTime[]; out = {}; emit[s_] := (AppendTo[out, s]; Print[s]);
c = 1 + x; tidy[s_SeriesData] := MapAt[Together, s, 3]; tidy[s_] := Together[s];
tauS = Normal[Series[1 - Exp[-t], {t, 0, KX + 2}]];
U1 = ((x - 1)/x + ((1 + x)/x) Sqrt[1 + 4 x (Exp[-t] - 1)/(1 + x)^2])/2;
U1x = Normal[tidy[Series[U1, {t, 0, KX + 4}]]];
sst = Normal[tidy[Series[(1 - U1x)/tauS, {t, 0, KX + 2}]]];
dphi = Normal[MapAt[Simplify[PowerExpand[#], x > 0] &, Series[(Log[s] + x Log[1 - s] - Log[1 - tauS s] /. s -> sst) - (Log[1/c] + x Log[x/c]), {t, 0, KX + 2}], 3]];
(* closed form of the row, as in hvc_leaf_proof_check.wls *)
F21 = Sum[(tauS/z)^k Product[(1 + j z)^2, {j, 1, k}]/(k! Product[c + j z, {j, 2, k + 1}]), {k, 0, KX + 1}];
zG = tidy[Series[(c z/(c + z)) (tauS/t) F21 Exp[-dphi/z], {t, 0, KX}]];
(* A_j(t) as t-series *)
hyp[j_] := Sum[Pochhammer[j x/c, m]^2/(Pochhammer[j, m] m!) tauS^m, {m, 0, KX + 1}];
Aj[j_] := Aj[j] = tidy[Series[(-1)^j c^2 (tauS/t) Exp[j dphi/c] tauS^(j - 1) Pochhammer[1 - j/c, j - 1]^2/(j ((j - 1)!)^2) hyp[j], {t, 0, KX}]];
pred = tidy[Series[c + Sum[Normal[Aj[j]]/(c + j z), {j, 1, KX + 1}], {t, 0, KX}]];
emit["(a) z G = c + Sum_j A_j(t)/(c + j z) through t^" <> ToString[KX] <> ": " <> ToString[Table[Together[SeriesCoefficient[zG, k] - SeriesCoefficient[pred, k]] === 0, {k, 0, KX}]]];
emit["(b) Sum_j A_j(t) = -c^2 through t^" <> ToString[KX] <> ": " <> ToString[Normal[tidy[Sum[Aj[j], {j, 1, KX + 1}] + c^2]] === 0]];
emit["    order of vanishing of A_j at t = 0 (expected j - 1): " <> ToString[Table[Module[{s = Aj[j]}, If[Normal[s] === 0, "all", s[[4]]]], {j, 1, Min[KX + 1, 8]}]]];
Print["time ", Round[AbsoluteTime[] - t0], " s"];
(* (c) the stored fits: dilRes[[k+1, 1]] = coefficients of lam_0..lam_{k+1} of rho q_t - q_0 at t^k *)
{rhoC, dilRes5, legRes5, loopRes5, loopResR5} = Get["hvc_leaf_structure_data_kt5.m"];
q0lam = Table[If[j == 1, -c^2, 0], {j, 0, 7}];
dilOK = Table[Module[{fit = dilRes5[[k + 1, 1]]}, And @@ Table[Together[fit[[j + 1]] - (If[j == 0, 0, SeriesCoefficient[Aj[j], k]] - If[k == 0, q0lam[[j + 1]], 0])] === 0, {j, 0, k + 1}]], {k, 0, 5}];
emit["(c1) A_j against the stored dilaton fits (rho q_t - q_0, lam_0..lam_{k+1}), k = 0..5: " <> ToString[dilOK]];
(* loops: e_il(t) *)
Theta = Normal[tidy[Series[-x^3 tauS sst^2, {t, 0, KX + 1}]]];
eil[i_, l_] := eil[i, l] = Module[{integrand = Normal[tidy[Series[Theta D[t Normal[Aj[i]], t] D[t Normal[Aj[l]], t], {t, 0, 7}]]]},
    (i l/c^2) Integrate[integrand, {t, 0, t}]];
loopOK = Table[Module[{sol = loopResR5[[k + 1, 1]]},
    And @@ Table[Module[{p = rule[[1]], val = rule[[2]], pred},
       pred = Together[Coefficient[Expand[eil[p[[1]], p[[2]]]], t, k + 2]];
       Together[pred - val] === 0], {rule, sol}]], {k, 0, 5}];
emit["(c2) loop coefficients e_il from the flow against the stored restricted fits (rho^2 E11 at t^k, pairs i <= l, (i-1)+(l-1) <= k), k = 0..5: " <> ToString[loopOK]];
emit["total time " <> ToString[Round[AbsoluteTime[] - t0]] <> " s"];
Export["hvc_leaf_closed_forms_out.txt", StringRiffle[out, "\n"], "Text"];
