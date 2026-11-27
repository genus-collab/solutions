(* Plain-text export of localp2final_CY5_v16.nb (2026-11-27): Input cells as code, text cells as comments; outputs omitted. *)

(* ::Title:: Local P2 x C2: the fibrewise R-matrix to order z^6 *)

(* ::Text:: Derivation of the identity column cPP[0,i][k], k = 0..6, of the canonical R-matrix of K_P2 x C^2 restricted to the fibrewise torus T_F (alpha_i = (eps1+eps2)/3), following the route of localp2final_CY5_v15.nb: (i) the ansatz solution of the equivariant Picard-Fuchs equation on the rational sub-torus T_R, (ii) normalisation at y = 0 and the quantum Riemann-Roch materialisation constants Qrrc over all five tangent weights, (iii) the delta -> 0 limit along alpha_l = (eps1+eps2)/(3 (1 - omega^l delta)). This notebook supersedes the subsections 'Picard-Fuchs calculations' and 'The R-matrix' of v15 (whose k <= 3 results it reproduces exactly); the rest of v15 is unchanged. *)

(* ::Text:: The organisation differs from v15 in two respects that make the computation cheap:
1. Along the delta-curve, e1 = 3 alpha/(1-u), e2 = 3 alpha^2/(1-u), e3 = alpha^3/(1-u) with u = delta^3, so the PF solution depends on delta only through u, and the cubic for the canonical coordinates is explicit: (X + alpha)^3 = X^3 (u - 27 y (1 - u)). The ansatz is therefore solved on this one-parameter family (with alpha = 1, restored by homogeneity pp[k][l][m] ~ alpha^(8k-l)), at rational sample values of u, imposing the identity in (y, X) modulo the cubic by sampling rational X and y = y(X) on the cubic; each coefficient is then reconstructed as an exact rational function of u.
2. The y = 0 normalisation and the QRR factor form a y-independent series C(delta, z). The limit is assembled as R_k^{TF}(y) = sum_{a+b=k} sum_j [u^j] PPunnorm_a(y, X0(y;u); u) [delta^(-3j)] C_b(delta), with all other polar terms checked to cancel. This replaces the expansion in delta of one very large rational function. *)

(* ::Text:: Evaluate top to bottom. Timings on the original machine: levels 0-3 of the ansatz about 1 minute, level 4 about 7 minutes, levels 5-6 substantially longer (see the log printed by the cells). The results are also stored in the sidecar files PFu6.mx (ansatz coefficients) and cPPu6.mx (R-matrix coefficients), which the genus-3 notebook loads. *)

(* ::Section:: The Picard-Fuchs ansatz on the delta-family *)

(* Part A1 (numeric-sampling variant): v15's PF ansatz on the delta-curve family
     ee1 = 3/(1-u), ee2 = 3/(1-u), ee3 = 1/(1-u)  (alpha = 1),
   solved at rational sample values u_s (all lower levels numeric at the same u_s), then each pp[k][l][m](u)
   is reconstructed as P(u)/(u^a (1-u)^b) by polynomial interpolation with verification on extra samples.
   Usage: wolframscript -file partA_pf_num.wl KMAX
   Output: PFu<k>.mx (PFu rules, exact rational functions of u) after each level k. *)

KMAX = 6;
SetDirectory[NotebookDirectory[]];
log[m_] := Module[{s = OpenAppend["partA_pf_num.log"]}, WriteString[s, DateString[{"Hour", ":", "Minute", ":", "Second"}], "  ", m, "\n"]; Close[s]; Print[m]];

ee1 = 3/(1 - u); ee2 = 3/(1 - u); ee3 = 1/(1 - u);
polPFL[yy_, X_] := ee1 X^2 + ee2 X + ee3 + (1 + 27 yy) X^3;
ySol[X_] := -(ee1 X^2 + ee2 X + ee3 + X^3)/(27 X^3);
Module[{Pf, r1, r2, r3, Xf},
  Pf[a_, b_] := polPFL[a, b];
  r1[x_, X_] = Together[-D[Pf[x, X], x]/D[Pf[x, X], X]];
  r2[x_, X_] = Together[D[r1[x, Xf[x]], x] /. Xf'[x] -> r1[x, Xf[x]] /. Xf[x] -> X];
  r3[x_, X_] = Together[D[r2[x, Xf[x]], x] /. Xf'[x] -> r1[x, Xf[x]] /. Xf[x] -> X];
  Xd1[x_, X_] = r1[x, X]; Xd2[x_, X_] = r2[x, X]; Xd3[x_, X_] = r3[x, X];];

ZMAX = KMAX;
IAs[y_] = Exp[mm[y]/z] Sum[S0[k][y] z^k, {k, 0, ZMAX + 1}];
pf = Expand[Exp[-mm[y]/z] Block[{coeff = {ee3, ee2 z, ee1 z^2, z^3},
      coeff2 = CoefficientList[Series[(3 z DY) (3 z DY + z) (3 z DY + 2 z), {DY, 0, 3}], DY]},
     Sum[coeff[[n + 1]] D[IAs[Exp[Y]], {Y, n}], {n, 0, 3}] + y Sum[coeff2[[n + 1]] D[IAs[Exp[Y]], {Y, n}], {n, 0, 3}]] /. Y -> Log[y]
   /. {mm'[y] -> XX[y]/y, mm''[y] -> D[XX[y]/y, y], mm'''[y] -> D[XX[y]/y, {y, 2}]}];
pfk[k_] := pfk[k] = Coefficient[pf, z, k];
log["pf built"];

S0ans[k_][yy_, X_] := Sum[pp[k][l][m] X^l yy^m, {l, 0, 2}, {m, 0, 2 k}]/(3 ee3 + ee2 X)^(k + 2)/((1 + 27 yy) (-ee2^3 + 27 (1 + 27 yy) ee3^2))^k;
Sder[k_, n_] := Sder[k, n] = Module[{f = S0ans[k][yy, Xf[yy]]}, D[f, {yy, n}] /. {Derivative[3][Xf][yy] -> X3, Derivative[2][Xf][yy] -> X2, Derivative[1][Xf][yy] -> X1, Xf[yy] -> X0}];
subsAnsatz[k_] := subsAnsatz[k] = Module[{e = pfk[k]},
   e = e /. {Derivative[n_][S0[j_]][y] :> Sder[j, n], S0[j_][y] :> Sder[j, 0]};
   e /. {Derivative[n_][XX][y] :> {X1, X2, X3}[[n]], XX[y] -> X0, y -> yy}];

samplesX = Table[-2/3 - 1/(7 + 3 s), {s, 0, 80}];
usamples = Table[1/(4 + 3 s) + 1/(100 + 7 s), {s, 0, 140}];    (* rational u in (0,1) *)
NUS[k_] := Max[24, 2 (7 k - 1) + 8];  (* initial number of u samples per level; increased if reconstruction fails *)
known = <||>;
(* reuse the highest available PFu<k>.mx (each contains all levels <= k) *)
Module[{avail = Select[Range[0, KMAX], FileExistsQ["PFu" <> ToString[#] <> ".mx"] &]},
  If[avail =!= {}, Get["PFu" <> ToString[Max[avail]] <> ".mx"];
   Do[known[k] = Select[PFu, #[[1, 0, 0, 1]] == k &], {k, 0, Max[avail]}]]];
Print["levels loaded from disk: ", Keys[known]];

(* numeric values: val[k][u_s] = association pp[k][l][m] -> rational *)
val[k_][us_] := val[k][us] = If[KeyExistsQ[known, k], Thread[known[k][[All, 1]] -> (known[k][[All, 2]] /. u -> us)], Module[{unknowns, nU, e, eqs = {}, sol, lower},
   unknowns = Flatten[Table[pp[k][l][m], {l, 0, 2}, {m, 0, 2 k}]]; nU = Length[unknowns];
   lower = Flatten[Table[val[j][us], {j, 0, k - 1}]];
   e = subsAnsatz[k + 1] /. u -> us /. lower;
   Do[Module[{Xs = samplesX[[s]], ys, vals}, ys = ySol[Xs] /. u -> us;
     vals = {yy -> ys, X0 -> Xs, X1 -> Xd1[ys, Xs] /. u -> us, X2 -> Xd2[ys, Xs] /. u -> us, X3 -> Xd3[ys, Xs] /. u -> us};
     AppendTo[eqs, e /. vals]], {s, 1, nU + 3}];
   sol = Solve[Join[Thread[eqs == 0], {pp[k][0][0] == KroneckerDelta[k, 0]}], unknowns];
   If[Length[sol] != 1 || Length[Flatten[sol]] != nU, log["level " <> ToString[k] <> " u=" <> ToString[us] <> ": bad solve " <> ToString[Length[sol]]]; Abort[]];
   Flatten[sol]]];

(* rational reconstruction from exact samples {u_s, v_s}: N(u)/D(u) with deg N + deg D + 1 = n - 5 (balanced),
   one NullSpace computation; verified on the 5 spare points *)
reconstruct[pts_] := Module[{n = Length[pts], us = pts[[All, 1]], vs = pts[[All, 2]], m, dd, dn, A, ns, v, N0, D0, ok},
   If[Union[vs] === {0}, Return[0]];
   m = n - 5; dd = Floor[(m - 1)/2]; dn = m - 1 - dd;
   A = Table[Join[Table[us[[s]]^j, {j, 0, dn}], Table[-vs[[s]] us[[s]]^j, {j, 0, dd}]], {s, 1, m}];
   ns = NullSpace[A]; If[ns === {}, Return[$Failed]];
   v = First[ns]; N0 = Sum[v[[j + 1]] u^j, {j, 0, dn}]; D0 = Sum[v[[dn + 2 + j]] u^j, {j, 0, dd}];
   If[D0 === 0, Return[$Failed]];
   ok = And @@ Table[((N0 /. u -> us[[s]]) - vs[[s]] (D0 /. u -> us[[s]])) == 0, {s, m + 1, n}];
   If[ok, Together[N0/D0], $Failed]];

PFu = {};
Do[
  Module[{t0 = AbsoluteTime[], keys, nus = NUS[k], recon, done = False},
   keys = Flatten[Table[pp[k][l][m], {l, 0, 2}, {m, 0, 2 k}]];
   If[KeyExistsQ[known, k], PFu = Join[PFu, known[k]]; log["level " <> ToString[k] <> " loaded from disk"]; Continue[]];
   While[! done,
    Do[val[k][usamples[[s]]], {s, 1, nus}];
    recon = Table[key -> reconstruct[Table[{usamples[[s]], key /. val[k][usamples[[s]]]}, {s, 1, nus}]], {key, keys}];
    If[FreeQ[recon, $Failed], done = True,
     log["level " <> ToString[k] <> ": reconstruction failed with " <> ToString[nus] <> " u-samples; increasing"]; nus += 8;
     If[nus > Length[usamples], log["too many samples needed -- ABORT"]; Abort[]]]];
   PFu = Join[PFu, recon];
   DumpSave["PFu" <> ToString[k] <> ".mx", PFu];
   log["level " <> ToString[k] <> " done: " <> ToString[Length[keys]] <> " coefficients, " <> ToString[nus] <> " u-samples, " <>
     ToString[Round[AbsoluteTime[] - t0, 0.1]] <> " s; max LeafCount " <> ToString[Max[LeafCount /@ recon[[All, 2]]]]];
   ], {k, 0, KMAX}];
log["saved PFu" <> ToString[KMAX] <> ".mx"];

If[FileExistsQ["v15data.mx"],
  Module[{ok, save = PFu},
   Get["v15data.mx"]; PFu = save;
   ok = Table[Module[{v15v = r[[2]] /. {ee[2] -> 3/(1 - u), ee[3] -> 1/(1 - u)}}, Simplify[(r[[1]] /. PFu) - v15v] === 0],
     {r, Select[PFAsRepl, MemberQ[PFu[[All, 1]], #[[1]]] &]}];
   log["agreement with v15 PFAsRepl on the curve: " <> ToString[Union[ok]] <> " (" <> ToString[Length[ok]] <> " coefficients)"]]];

(* ::Section:: Normalisation, QRR constants and the delta -> 0 limit *)

(* ::Text:: Requires PFu6.mx from the previous section (or from the sidecar). Produces cPPt (functions of t = y^(1/3)) and cPPL (functions of y and LLL[0][y]) and compares with the stored v15 closed forms for k <= 3. *)

(* Part A2 (u-family version): fibrewise R-matrix identity column cPP[0,0][k], k <= KMAX, by the v15 route
   (T_R ansatz + QRR normalisation over all five weights + delta -> 0 limit along alpha_l = alpha/(1 - w^l delta)).
   Input: PFu<KMAX>.mx from partA_pf_u.wl (ansatz coefficients on the curve, a = 1, rational in u = delta^3).
   Both the PF solution (Laurent in u) and the y-independent normalisation series C(delta,z) (Laurent in delta)
   are expanded, and the finite part of the product is assembled; all polar terms are checked to cancel.
   Usage: wolframscript -file partA_limit_u.wl KMAX
   Output: cPPu<KMAX>.mx with cPPt[k] (functions of t = y^(1/3), eps, sq) and cPPL[k] (in y, LLL[0][y]). *)

KMAX = 6;
SetDirectory[NotebookDirectory[]];
Get["PFu" <> ToString[KMAX] <> ".mx"];
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
DumpSave["cPPu" <> ToString[KMAX] <> ".mx", {cPPt, cPPL, J, Jneg, pord}];
log["saved cPPu" <> ToString[KMAX] <> ".mx; LeafCounts " <> ToString[LeafCount /@ cPPt]];

(* numeric comparison with v15 closed forms where available *)
If[FileExistsQ["v15data.mx"],
  Module[{sqrule, pt, ok},
   Get["v15data.mx"];
   sqrule = Power[(\[Epsilon][1] \[Epsilon][2])/(\[Epsilon][1] + \[Epsilon][2])^3, r_] :> sq^(2 r);
   pt = {t -> 1/2, \[Epsilon][1] -> 2, \[Epsilon][2] -> 1, sq -> Sqrt[2/27]};
   ok = Table[Module[{ref = cPP[0, 0][k] /. LLL[0][y] -> -alpha/(1 + 3 t) /. y -> t^3 /. sqrule},
       N[(cPPt[[k + 1]] - ref) /. pt, 30]], {k, 0, Min[3, KMAX]}];
   log["numeric difference from v15 cPP[0,0][0..3] at a sample point: " <> ToString[ok]]]];

(* ::Section:: Results *)

(* ::Text:: The coefficients cPP[0,i][k] = cPPL[[k+1]] /. LLL[0][y] -> LLL[i][y] in the mixed basis of v15 (flat basis for the argument, normalised canonical basis for the image), with sq = Sqrt[eps1 eps2/(eps1+eps2)^3]. *)

Get["cPPu6.mx"]; Do[cPP[0, i_][k] = cPPL[[k + 1]] /. LLL[0][y] -> LLL[i][y], {k, 0, 6}]; Table[LeafCount[cPP[0, 0][k]], {k, 0, 6}]

cPP[0, 0][4]

cPP[0, 0][5]

cPP[0, 0][6]
