(* Part A2 (u-family version): fibrewise R-matrix identity column cPP[0,0][k], k <= KMAX, by the v15 route
   (T_R ansatz + QRR normalisation over all five weights + delta -> 0 limit along alpha_l = alpha/(1 - w^l delta)).
   Input: PFu<KMAX>.mx from partA_pf_u.wl (ansatz coefficients on the curve, a = 1, rational in u = delta^3).
   Both the PF solution (Laurent in u) and the y-independent normalisation series C(delta,z) (Laurent in delta)
   are expanded, and the finite part of the product is assembled; all polar terms are checked to cancel.
   Usage: wolframscript -file genRmat_g4_partA_limit.m KMAX   (this file: the Part A2 cell of localp2final_CY5_v16.nb, unchanged except for the argument handling and the portable input/output PFu<K>.m, cPP_z<K>.m)
   Output: cPPu<KMAX>.mx with cPPt[k] (functions of t = y^(1/3), eps, sq) and cPPL[k] (in y, LLL[0][y]). *)

KMAX = ToExpression[$ScriptCommandLine[[2]]];
SetDirectory[DirectoryName[$InputFileName]];
PFu = Get["PFu" <> ToString[KMAX] <> ".m"];
log[m_] := Module[{s = OpenAppend["partA_limit_u.log"]}, WriteString[s, DateString[{"Hour", ":", "Minute", ":", "Second"}], "  ", m, "\n"]; Close[s]; Print[m]];

alpha = (\[Epsilon][1] + \[Epsilon][2])/3;
ee2u = 3 alpha^2/(1 - u); ee3u = alpha^3/(1 - u);
ppu[k_][l_][m_] := ppu[k][l][m] = Together[alpha^(8 k - l) (pp[k][l][m] /. PFu)];
S0u[k_][yy_, X_] := Sum[ppu[k][l][m] X^l yy^m, {l, 0, 2}, {m, 0, 2 k}]/(3 ee3u + ee2u X)^(k + 2)/
   ((1 + 27 yy) (-ee2u^3 + 27 (1 + 27 yy) ee3u^2))^k;

(* ---------- constants C(delta, z) ---------- *)
wred[e_] := Expand[Expand[e] /. w^n_ :> w^Mod[n, 3] /. w^2 -> -1 - w];
alphaD[l_] := alpha/(1 - w^l d);
invwt[1] = (1 - w d) (1 - d)/(alpha d) (-(w + 2)/3);
invwt[2] = (1 - w^2 d) (1 - d)/(alpha d) ((w - 1)/3);
Qrrc[k_] := Qrrc[k] = Together[-BernoulliB[2 k + 2]/((2 k + 2) (2 k + 1)) (
     (3 alphaD[0])^(-(2 k + 1)) + (-\[Epsilon][1])^(-(2 k + 1)) + (-\[Epsilon][2])^(-(2 k + 1)) +
      wred[invwt[1]^(2 k + 1) + invwt[2]^(2 k + 1)])];
If[! FreeQ[Table[Qrrc[k], {k, 0, KMAX}], w], log["ERROR: Qrrc not real"]; Abort[]];
ffz = Normal[Series[Exp[Sum[Qrrc[k] z^(2 k + 1), {k, 0, KMAX}]], {z, 0, KMAX}]];
Sk0[k_] := Sk0[k] = Together[S0u[k][yy, -alphaD[0]] /. yy -> 0 /. u -> d^3];
PPun0 = Sum[Sk0[k] z^k, {k, 0, KMAX}];
sqrtDcl = Sqrt[27] sq Sqrt[1 - d^3] (3 ee3u - ee2u alphaD[0]) /. u -> d^3;
normLRz = Normal[Series[1/(sqrtDcl PPun0), {z, 0, KMAX}]];
Cz = Normal[Series[normLRz ffz, {z, 0, KMAX}]];
cb = Table[Together[Coefficient[Cz, z, b]], {b, 0, KMAX}];
pord = Table[Module[{s = Series[cb[[b + 1]], {d, 0, 0}]}, If[s === 0, 0, -s[[4]]]], {b, 0, KMAX}];
log["pole orders of C_b(delta): " <> ToString[pord]];
J = Max[0, Ceiling[Max[pord]/3]];

(* ---------- Laurent expansion in u of PPunnorm along the curve ---------- *)
q = (1 + 27 t^3)/(27 t^3);
Jneg = 0;
Do[Module[{s = Series[S0u[a][t^3, alpha/(-3 t (1 - u q)^(1/3) - 1)], {u, 0, 0}]}, If[s =!= 0, Jneg = Max[Jneg, -s[[4]]]]], {a, 0, KMAX}];
log["u-pole order of PPunnorm: " <> ToString[Jneg] <> "; delta-pole order of C: " <> ToString[Max[pord]] <> "; J = " <> ToString[J]];
cbrt = -3 t Normal[Series[(1 - u q)^(1/3), {u, 0, J + Jneg}]];
X0u = alpha/(cbrt - 1);
Qja = Table[Module[{s = Series[S0u[a][t^3, X0u], {u, 0, J}]},
    Association[Table[j -> Together[SeriesCoefficient[s, j]], {j, -Jneg, J}]]], {a, 0, KMAX}];
log["u-Laurent expansion of PPunnorm done"];
cLaurent = Table[Module[{s = Series[cb[[b + 1]], {d, 0, 3 Jneg}]},
    Association[Table[m -> Together[SeriesCoefficient[s, m]], {m, -3 J, 3 Jneg}]]], {b, 0, KMAX}];
log["delta-Laurent expansion of C done"];
cL[b_][m_] := If[m < -3 J || m > 3 Jneg, 0, cLaurent[[b + 1]][m]];
Qj[a_][j_] := If[j < -Jneg || j > J, 0, Qja[[a + 1]][j]];

(* ---------- assemble: coefficient of delta^p is sum_{a+b=k} sum_j Qj[a][j] cL[b][p-3j] ---------- *)
coeffP[k_, p_] := Together[Sum[Sum[Qj[a][j] cL[k - a][p - 3 j], {j, -Jneg, J}], {a, 0, k}]];
cPPt = Table[coeffP[k, 0], {k, 0, KMAX}];
checks = Table[coeffP[k, p], {k, 0, KMAX}, {p, -3 J - 3 Jneg, -1}];
log["cancellation of polar terms (all zero?): " <> ToString[Union[Flatten[checks]] === {0}]];
If[Union[Flatten[checks]] =!= {0}, log["nonzero polar terms at {k+1, index}: " <> ToString[Position[checks, Except[0], {2}, Heads -> False]]]];
cPPL = Table[Together[cPPt[[k + 1]] /. t -> -(LLL[0][y] + alpha)/(3 LLL[0][y])], {k, 0, KMAX}];
Put[cPPL, "cPP_z" <> ToString[KMAX] <> ".m"];   (* portable: the list cPPL[[k+1]] = cPP[0,0][k], k = 0..KMAX, read by genRmat_g4_lib.m *)
Put[{cPPt, cPPL, J, Jneg, pord}, "cPPu" <> ToString[KMAX] <> "_all.m"];
log["saved cPP_z" <> ToString[KMAX] <> ".m; LeafCounts " <> ToString[LeafCount /@ cPPt]];

(* numeric comparison with v15 closed forms where available *)
If[FileExistsQ["v15data.mx"],
  Module[{sqrule, pt, ok},
   Get["v15data.mx"];
   sqrule = Power[(\[Epsilon][1] \[Epsilon][2])/(\[Epsilon][1] + \[Epsilon][2])^3, r_] :> sq^(2 r);
   pt = {t -> 1/2, \[Epsilon][1] -> 2, \[Epsilon][2] -> 1, sq -> Sqrt[2/27]};
   ok = Table[Module[{ref = cPP[0, 0][k] /. LLL[0][y] -> -alpha/(1 + 3 t) /. y -> t^3 /. sqrule},
       N[(cPPt[[k + 1]] - ref) /. pt, 30]], {k, 0, Min[3, KMAX]}];
   log["numeric difference from v15 cPP[0,0][0..3] at a sample point: " <> ToString[ok]]]];